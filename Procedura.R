library(TSA)
library(npreg)
library(forecast)
library(stats)
library(fracdiff)
library(tseriesChaos)

######
#periodicity
######
#mod_ss <- ss(AirPassengers,nknots = 3)
########

remove_trend <- function(ts_data, nknots = 3) {
  # Przekształcenie szeregu czasowego na macierz dwukolumnową
  ts_matrix <- cbind(time(ts_data), ts_data)
  
  # Dopasowanie trendu przy użyciu smooth spline
  trend_fit <- ss(ts_data, nknots = nknots)
  
  # Usunięcie dopasowanego trendu z szeregu czasowego
  detrended_ts <- ts_data - trend_fit$y
  
  # Zwrócenie przetworzonego szeregu czasowego
  return(detrended_ts)
}

#plot(remove_trend(AirPassengers))
#acf(remove_trend(AirPassengers))

find_peaks_and_troughs <- function(ts_data) {
  acf_data <- acf(ts_data, lag.max = floor(length(ts_data)/3), plot = FALSE)
  peaks <- c()
  troughs <- c()
  if (acf_data$acf[1] > acf_data$acf[2]){peaks <- c(peaks, 1)}
  else if (acf_data$acf[1] < acf_data$acf[2]){troughs <- c(troughs, 1)}
  if (acf_data$acf[length(acf_data$acf)] > acf_data$acf[length(acf_data$acf)-1]){
    peaks <- c(peaks, length(acf_data$acf))
  }
  else if (acf_data$acf[length(acf_data$acf)] < acf_data$acf[length(acf_data$acf)-1]){
    troughs <- c(troughs, length(acf_data$acf))
  }
  
  for(i in 2:(length(acf_data$acf) - 1)) {
    if(acf_data$acf[i] > acf_data$acf[i-1] && acf_data$acf[i] > acf_data$acf[i+1]) {
      peaks <- c(peaks, i)
    } else if(acf_data$acf[i] < acf_data$acf[i-1] && acf_data$acf[i] < acf_data$acf[i+1]) {
      troughs <- c(troughs, i)
    }
  }
  peaks <- sort(peaks)
  troughs <- sort(troughs)
  return(list(peaks = peaks, troughs = troughs))
}

find_seasonality <- function(ts_data) {
  #detrend time series
  detrended_ts <- remove_trend(ts_data)
  
  # calculate autocorrelation function for lags up to 1/3 of series length
  acf_data <- acf(detrended_ts, lag.max = floor(length(detrended_ts)/3), 
                  plot = FALSE)
  
  # find peaks and troughs in autocorrelation function
  peaks <- find_peaks_and_troughs(acf_data$acf)$peaks
  troughs <- find_peaks_and_troughs(acf_data$acf)$troughs
  if (peaks[1] < troughs[1]){k <- 2}
  else {k <- 1}
  
  # check if there is a peak that satisfies conditions (a), (b), and (c)
  for (i in k:length(peaks)) {
    if (length(troughs[troughs<peaks[i]]) >= 1 & 
        abs(acf_data$acf[peaks[i]] - acf_data$acf[max(troughs[troughs<peaks[i]])]) >= 0.1 &
        acf_data$acf[peaks[i]] > 0) {
      return(peaks[i]) # return frequency of seasonality
    }
  }
  
  # if no peak satisfies the conditions, set frequency to 1 (non-seasonal)
  return(1)
}

#find_seasonality(AirPassengers)
#find_peaks_and_troughs(remove_trend(AirPassengers))

####
#dekompozycja szeregu
####

stl_decomp <- function(ts_data, freq){
  stl_result <- stl(ts_data, s.window = "periodic", t.window = freq, robust = TRUE)
  
  return(list(seasonal = stl_result$time.series[,1], trend = stl_result$time.series[,2],
              remainder = stl_result$time.series[,3]))
}

# ts_after_BoxCox <- BoxCox(AirPassengers, lambda = "auto")
# ts_after_stl <- stl_decomp(ts_after_BoxCox,find_seasonality(AirPassengers))
# length(ts_after_stl$remainder)
# 
# X_t <- ts_after_BoxCox - ts_after_stl$trend
# Z_t <- ts_after_BoxCox - ts_after_stl$seasonal
# R_t <- ts_after_BoxCox - ts_after_stl$trend - ts_after_stl$seasonal
# 
# trend_measure <- 1-var(R_t)/var(Z_t)
# seasonality_measure <- 1 - var(R_t)/var(X_t)

Trend_seasonality_measure <- function(ts_data, lambda = "auto"){
  ts_minimum <- min(ts_data)
  if (ts_minimum > 0){
    ts_data <- BoxCox(ts_data, lambda = lambda)
  }
  else if (ts_minimum == 0){
    ts_data <- ts_data + 0.001*max(ts_data)
    ts_data <- BoxCox(ts_data, lambda = lambda)
  }
  
  freq <- find_seasonality(ts_data = ts_data)
  if (freq > 1){
    stl_ts <- stl_decomp(ts_data = ts_data, freq=freq)
    X_t <- ts_data - stl_ts$trend
    Z_t <- ts_data - stl_ts$seasonal
    R_t <- ts_data - stl_ts$trend - stl_ts$seasonal
  }
  else if (freq == 1){
    ss_ts_trend <- ss(ts_data, nknots = 3)$y
    X_t <- ts_data - ss_ts_trend
    Z_t <- ts_data
    R_t <- ts_data - ss_ts_trend
  }
  trend_measure <- 1 - var(R_t)/var(Z_t)
  seasonality_measure <- 1 - var(R_t)/var(X_t)
  
  return(list(trend_measure=trend_measure, seasonality_measure=seasonality_measure))
}

#Trend_seasonality_measure(AirPassengers, lambda = 0.5)

####
#serial corelation
####
Box_Pierce_stats <- function(ts_data){
  #value of stattistic for RAW data
  Box_Pierce_RAW <- Box.test(ts_data, lag = 20)$statistic
  #TSA data using stl(when freq > 1) or ss(when freq = 1)
  freq <- find_seasonality(ts_data = ts_data)
  if (freq > 1){
    stl_ts <- stl_decomp(ts_data = ts_data, freq=freq)
    R_t <- ts_data - stl_ts$trend - stl_ts$seasonal
  }
  else if (freq == 1){
    ss_ts_trend <- ss(ts_data, nknots = 3)$y
    R_t <- ts_data - ss_ts_trend
  }
  Box_Pierce_TSA <- Box.test(R_t, lag = 20)$statistic
  
  return(list(Box_Pierce_RAW = Box_Pierce_RAW, Box_Pierce_TSA = Box_Pierce_TSA))
}

#Box_Pierce_stats(AirPassengers)

####
#Non-linear AR structure
####
#Tsay.test()

####
#Skewness
####
#gotowa funkcja z pakietu TSA
TSA::skewness(ts_data)
#ze wzoru z artykułu characteristic-based...
skewness_manual <- function(ts_data){
  skewness <- 1/(length(ts_data)*sd(ts_data)^3)*sum((ts_data-mean(ts_data))^3)
  return(skewness=skewness)
}

####
#Kurtosis (heavy-tails)
####
#gotowa funkcja z pakietu TSA
TSA::kurtosis(ts_data)
kurtosis_manual <- function(ts_data){
  kurtosis <- 1/(length(ts_data)*sd(ts_data)^4)*sum((ts_data-mean(ts_data))^4)-3
  return(kurtosis = kurtosis)
}

####
#self-dependency
####
self_similarity <- function(ts_data){
  d <- fracdiff(ts_data)$d
  Hurst <- d + 0.5
  return(Hurst=Hurst)
}
self_similarity(AirPassengers)

####
#chaos (dynamic systems)
####
output <-lyap_k(AirPassengers, m=1, d=2, s=10, t=10, ref=100, k=2, eps=4)
plot(output)
lyap(output, 0.2, 0.4)
