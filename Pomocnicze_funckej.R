library(forecast)
library(tseries)
library(fracdiff)
library(fpc)
library(klaR)

find.freq <- function(x){
  n <- length(x)
  spec <- spec.ar(c(na.contiguous(x)),plot=FALSE)
  if(max(spec$spec)>10) # Arbitrary threshold chosen by trial and error.
  {
    period <- round(1/spec$freq[which.max(spec$spec)])
    if(period==Inf) # Find next local maximum
    {
      j <- which(diff(spec$spec)>0)
      if(length(j)>0)
      {
        nextmax <- j[1] + which.max(spec$spec[j[1]:500])
        if(nextmax <= length(spec$freq))
          period <- round(1/spec$freq[nextmax])
        else
          period <- 1
      }
      else
        period <- 1
    }
  }
  else
    period <- 1
  
  return(period)
}

##########

decomp <- function(x,transform=TRUE){
  require(forecast)
  # Transform series
  if(transform & min(x,na.rm=TRUE) >= 0)
  {
    lambda <- BoxCox.lambda(na.contiguous(x))
    x <- BoxCox(x,lambda)
  }
  else
  {
    lambda <- NULL
    transform <- FALSE
  }
  # Seasonal data
  if(frequency(x)>1)
  {
    x.stl <- stl(x,s.window="periodic",na.action=na.contiguous)
    trend <- x.stl$time.series[,2]
    season <- x.stl$time.series[,1]
    remainder <- x - trend - season
  }
  else #Nonseasonal data
  {
    require(mgcv)
    tt <- 1:length(x)
    trend <- rep(NA,length(x))
    trend[!is.na(x)] <- fitted(gam(x ~ s(tt)))
    season <- NULL
    remainder <- x - trend
  }
  return(list(x=x,trend=trend,season=season,remainder=remainder,
              transform=transform,lambda=lambda))
}

#######

# f1 maps [0,infinity) to [0,1]
f1 <- function(x,a,b){
  eax <- exp(a*x)
  if (eax == Inf)
    f1eax <- 1
  else
    f1eax <- (eax-1)/(eax+b)
  return(f1eax)
}

# f2 maps [0,1] onto [0,1]
f2 <- function(x,a,b){
  eax <- exp(a*x)
  ea <- exp(a)
  return((eax-1)/(eax+b)*(ea+b)/(ea-1))
}

#######

measures <- function(x){
  require(forecast)
  
  N <- length(x)
  freq <- find.freq(x)
  #
  if (2*freq >= N){freq <- 1}
  #
  fx <- c(frequency=(exp((freq-1)/50)-1)/(1+exp((freq-1)/50)))
  x <- ts(x,f=freq)
  
  # Decomposition
  decomp.x <- decomp(x)
  
  # Adjust data
  if(freq > 1)
    fits <- decomp.x$trend + decomp.x$season
  else # Nonseasonal data
    fits <- decomp.x$trend
  adj.x <- decomp.x$x - fits + mean(decomp.x$trend, na.rm=TRUE)
  
  # Backtransformation of adjusted data
  if(decomp.x$transform)
    tadj.x <- InvBoxCox(adj.x,decomp.x$lambda)
  else
    tadj.x <- adj.x
  
  # Trend and seasonal measures
  v.adj <- var(adj.x, na.rm=TRUE)
  if(freq > 1)
  {
    detrend <- decomp.x$x - decomp.x$trend
    deseason <- decomp.x$x - decomp.x$season
    trend <- ifelse(var(deseason,na.rm=TRUE) < 1e-10, 0, 
                    max(0,min(1,1-v.adj/var(deseason,na.rm=TRUE))))
    season <- ifelse(var(detrend,na.rm=TRUE) < 1e-10, 0,
                     max(0,min(1,1-v.adj/var(detrend,na.rm=TRUE))))
  }
  else #Nonseasonal data
  {
    trend <- ifelse(var(decomp.x$x,na.rm=TRUE) < 1e-10, 0,
                    max(0,min(1,1-v.adj/var(decomp.x$x,na.rm=TRUE))))
    season <- 0
  }
  
  m <- c(fx,trend,season)
  
  # Measures on original data
  xbar <- mean(x,na.rm=TRUE)
  s <- sd(x,na.rm=TRUE)
  
  # Serial correlation
  Q <- Box.test(x,lag=10)$statistic/(N*10)
  fQ <- f2(Q,7.53,0.103)
  
  # Nonlinearity
  p <- terasvirta.test(na.contiguous(x))$statistic
  fp <- f1(p,0.069,2.304)
  
  # Skewness
  skewness <- abs(mean((x-xbar)^3,na.rm=TRUE)/s^3)
  fs <- f1(skewness,1.510,5.993)
  
  # Kurtosis
  k <- mean((x-xbar)^4,na.rm=TRUE)/(s^4)
  fk <- f1(k,2.273,11567)
  
  # Hurst=d+0.5 where d is fractional difference.
  H <- fracdiff(na.contiguous(x),0,0)$d + 0.5
  
  # Lyapunov Exponent
  if(freq > N-10)
    stop("Insufficient data")
  Ly <- numeric(N-freq)
  for(i in 1:(N-freq))
  {
    idx <- order(abs(x[i] - x))
    idx <- idx[idx < (N-freq)]
    j <- idx[2]
    Ly[i] <- log(abs((x[i+freq] - x[j+freq])/(x[i]-x[j])))/freq
    if(is.na(Ly[i]) | Ly[i]==Inf | Ly[i]==-Inf)
      Ly[i] <- NA
  }
  Lyap <- mean(Ly,na.rm=TRUE)
  fLyap <- exp(Lyap)/(1+exp(Lyap))
  
  m <- c(m,fQ,fp,fs,fk,H,fLyap)
  
  # Measures on adjusted data
  xbar <- mean(tadj.x, na.rm=TRUE)
  s <- sd(tadj.x, na.rm=TRUE)
  
  # Serial
  Q <- Box.test(adj.x,lag=10)$statistic/(N*10)
  fQ <- f2(Q,7.53,0.103)
  
  # Nonlinearity
  p <- terasvirta.test(na.contiguous(adj.x))$statistic
  fp <- f1(p,0.069,2.304)
  
  # Skewness
  skewness <- abs(mean((tadj.x-xbar)^3,na.rm=TRUE)/s^3)
  fs <- f1(skewness,1.510,5.993)
  
  # Kurtosis
  k <- mean((tadj.x-xbar)^4,na.rm=TRUE)/s^4
  fk <- f1(k,2.273,11567)
  
  m <- c(m,fQ,fp,fs,fk)
  names(m) <- c("frequency", "trend","seasonal",
                "autocorrelation","non_linear","skewness","kurtosis",
                "Hurst","Lyapunov",
                "dc_autocorrelation","dc_non_linear","dc_skewness","dc_kurtosis")
  
  return(m)
}

####
#funkcja do ekstrakcji cech z macierzy szeregów czasowych
####
get_measure_df <- function(data_frame){
  measure_matrix <- matrix(NA, nrow = nrow(data_frame), ncol = 13)
  for (i in 1:nrow(data_frame)){
    measure_matrix[i,] <- measures(as.numeric(data_frame[i,]))
    print(i)
  }
  measure_matrix <- as.data.frame(measure_matrix)
  names(measure_matrix) <- c("frequency", "trend","seasonal",
                                 "autocorrelation","non_linear","skewness","kurtosis",
                                 "Hurst","Lyapunov",
                                 "dc_autocorrelation","dc_non_linear","dc_skewness","dc_kurtosis")
  
  return(measure_matrix)
}

###
#funkcja do obliczania macierzy odległości acf
###
acf_dist <- function(df){
  n <- nrow(df)
  dist_matrix <- matrix(nrow = n, ncol = n)
  for (i in 1:n){
    for (j in 1:n){
      dist_matrix[i,j] <- diss.ACF(as.numeric(df[i,]), as.numeric(df[j,]))
    }
  }
  return(dist_matrix)
}

###
#wybór cech przy pomocy ładunków PCA
###
pca_select_features <- function(data, min_var_explained=0.9){
  pca_model <- prcomp(data, scale. = TRUE)
  var_explained <- cumsum(pca_model$sdev^2) / sum(pca_model$sdev^2)
  n_components <- sum(var_explained >= min_var_explained)
  n_components <- ncol(data) + 1 - n_components
  rotations <- pca_model$rotation[,1:n_components]
  indeksy_zmiennych <- numeric(n_components)
  
  for (i in 1:n_components) {
    obciazenia_i <- abs(rotations[, i])
    indeks <- which.max(obciazenia_i)
    
    while (indeks %in% indeksy_zmiennych) {
      obciazenia_i[indeks] <- 0
      indeks <- which.max(obciazenia_i)
    }
    
    indeksy_zmiennych[i] <- indeks
  }
  
  return(indeksy_zmiennych)
}

###
#funkcja wyznaczająca kolumny ramki danych, które mają więcej unikalnych wartości niż liczba klastrów
###
unique_values <- function(df, k, digits=7){
  columns <- c()
  for (i in 1:ncol(df)){
    uniq_len <- length(unique(signif(df[,i], digits = digits)))
    if (uniq_len > k) {
      columns <- c(columns,i)
    }
  }
  return(columns)
}

###
#forward search (metoda nadzorowana)
###
forward <- function(df, krange, itter=1000, discrete = TRUE, track = TRUE){
  if (discrete == TRUE){
    n <- ncol(df)
    index <- c(1:n)
    char <- colnames(df)
  }
  else {
    x <- c("frequency","Hurst","dc_kurtosis")
    index <- c(1,8,13)
    del <- c()
    for (i in 1:length(index)){
      if ((sum(x[i] == colnames(df))) == 1){del <- c(del,index[i])}
    }
    df <- df[,-del]
    n <- ncol(df)
    index <- c(1:n)
    char <- colnames(df)
  }
  picked <- c()
  jaccard_index <- c(0)
  for (j in 1:n){
    jaccard <- c()
    for (i in 1:length(index)){
      stab <- mean(clusterboot(df[,c(picked,index[i])], B=itter, 
                              bootmethod = "subset", subtuning = floor(0.7*nrow(df)),
                             count = TRUE,clustermethod = kmeansCBI, krange = krange)$subsetmean)
      jaccard <- c(jaccard, stab)
    }
    names(jaccard) <- char
    best <- which.max(jaccard)
    print(jaccard[best])
    if (jaccard[best] <= jaccard_index[length(jaccard_index)]){
      break
    }
    else {
      picked <- c(picked,best)
      char <- char[-best]
      index <- index[-best]
      jaccard_index <- c(jaccard_index, jaccard[best])
    }
    if (track == TRUE){
      print(picked)
    }
  }
  if (discrete == FALSE){
    ind <- c(2,3,4,5,6,7,9,10,11,12)
    index <- ind[index]
  }
  
  return(list(picked=picked,jaccard_index=jaccard_index[-1], to_del_index=index))
}
