library(class)
library(clv)
library(FactoMineR)
###
#Procedura
###
Procedura_grupowania <- function(train, test, train_class, test_class, k, method = "all",
                                 hclust_method = "complete", B=100, min_var=0.7,
                                 grid_dim = c(6,6), topo = "hexagonal", subset=0.7,
                                 discrete = TRUE){
  #macierz odległości (euklidesowa oraz acf) do grupowania metodami referencyjnymi
  eucl_dist_matrix <- dist(test, diag = TRUE, upper = TRUE)
  acf_dist_matrix <- acf_dist(test)
  
  #grupowanie metodami referencyjnymi
  #euklidesowa
  kmean_eucl_ref <- kmeans(test, centers = k)
  pam_eucl_ref <- pam(eucl_dist_matrix, k = k, metric = "euclidean")
  hclust_eucl_ref <- hclust(eucl_dist_matrix, method = hclust_method)
  hclust_eucl_ref_clust <- cutree(hclust_eucl_ref, k=k)
  grid <- kohonen::somgrid(grid_dim[1], grid_dim[2], topo=topo)
  train_scaled <- scale(train)
  test_scaled <- scale(test, center = attr(train_scaled, "scaled:center"),
                        scale = attr(train_scaled, "scaled:scale"))
  training_data <- list(measurements = train_scaled,
                        class = as.factor(train_class))
  test_data <- list(measurements = test_scaled, class = as.factor(test_class))
  som_ref <- kohonen::supersom(training_data, grid = grid)
  som_pred_ref <- predict(som_ref, newdata = test_data[1])
  som_tab_ref <- table(test_class, som_pred_ref$predictions[["class"]])
  som_tab_ref_permut <- matchClasses(som_tab_ref, method = "exact", verbose = FALSE)
  which.na <- c(attr(na.omit(as.integer(som_pred_ref$predictions[["class"]])),"na.action"))
  l <- length(which.na)
  print(which.na)
  print(l)
  test_row <- nrow(test)
  #acf
  pam_acf_ref <- pam(acf_dist_matrix, k = k, metric = "euclidean")
  hclust_acf_ref <- hclust(as.dist(acf_dist_matrix), method = hclust_method)
  hclust_acf_ref_clust <- cutree(hclust_acf_ref, k=k)  
  print("grupowanie metodami ref ok")
  
  #ocena metod referencyjnych
  #accuracy
  kmean_eucl_ref_acc <- compareMatchedClasses(kmean_eucl_ref$cluster, test_class, method = "exact", verbose = FALSE)$diag
  pam_eucl_ref_acc <- compareMatchedClasses(pam_eucl_ref$clustering, test_class, method = "exact", verbose = FALSE)$diag
  hclust_eucl_ref_acc <- compareMatchedClasses(hclust_eucl_ref_clust, test_class, method = "exact",verbose = FALSE)$diag
  som_ref_acc <- sum(diag(som_tab_ref[,som_tab_ref_permut]))/test_row
  
  pam_acf_ref_acc <- compareMatchedClasses(pam_acf_ref$clustering, test_class, method = "exact", verbose = FALSE)$diag
  hclust_acf_ref_acc <- compareMatchedClasses(hclust_acf_ref_clust, test_class, method = "exact",verbose = FALSE)$diag
  
  #spójność (connectivity)
  kmean_eucl_conn <- connectivity.diss.mx(as.matrix(eucl_dist_matrix), kmean_eucl_ref$cluster, neighbour.num = 5)
  pam_eucl_conn <- connectivity.diss.mx(as.matrix(eucl_dist_matrix), pam_eucl_ref$clustering, neighbour.num = 5)
  hclust_eucl_conn <- connectivity.diss.mx(as.matrix(eucl_dist_matrix), hclust_eucl_ref_clust, neighbour.num = 5)
  if(l != 0){
    som_ref_conn <- connectivity(as.matrix(test)[-which.na,], c(na.omit(as.integer(som_pred_ref$predictions[["class"]]))), 
                                 neighbour.num = 5)
  }
  else{
    som_ref_conn <- connectivity(as.matrix(test), c(na.omit(as.integer(som_pred_ref$predictions[["class"]]))), 
                                 neighbour.num = 5)
  }
  
  pam_acf_conn <- connectivity.diss.mx(as.matrix(acf_dist_matrix), pam_acf_ref$clustering, neighbour.num = 5)
  hclust_acf_conn <- connectivity.diss.mx(as.matrix(acf_dist_matrix), hclust_acf_ref_clust, neighbour.num = 5)
  
  #indeks silhouette
  kmean_eucl_silh <- summary(silhouette(kmean_eucl_ref$cluster, dist = eucl_dist_matrix))$avg.width
  pam_eucl_silh <- summary(silhouette(pam_eucl_ref))$avg.width
  hclust_eucl_silh <- summary(silhouette(hclust_eucl_ref_clust, dist = eucl_dist_matrix))$avg.width
  if(l != 0){
    som_ref_silh <- summary(silhouette(na.omit(as.integer(som_pred_ref$predictions[["class"]])), 
                                       dist = dist(test_data$measurements[-which.na,], diag = TRUE,upper = TRUE)))$avg.width
  }
  else{
    som_ref_silh <- summary(silhouette(na.omit(as.integer(som_pred_ref$predictions[["class"]])), 
                                       dist = dist(test_data$measurements, diag = TRUE,upper = TRUE)))$avg.width
  }
  
  pam_acf_silh <- summary(silhouette(pam_acf_ref))$avg.width
  hclust_acf_silh <- summary(silhouette(hclust_acf_ref_clust, dist = acf_dist_matrix))$avg.width
  
  #stability
  kmean_eucl_stab <- mean(clusterboot(test, B=B, bootmethod = "subset", subtuning = floor(subset*nrow(test)),
                                      noisetuning = c(0,0), count = FALSE, clustermethod = kmeansCBI, krange=k, seed = 249730)$subsetmean)
  pam_eucl_stab <- mean(clusterboot(test, B=B, bootmethod = "subset", subtuning = floor(subset*nrow(test)),
                                    noisetuning = c(0,0), count = FALSE, clustermethod = claraCBI, k=k, seed = 249730)$subsetmean)
  hclust_eucl_stab <- mean(clusterboot(test, B=B, bootmethod = "subset", subtuning = floor(subset*nrow(test)),
                                      noisetuning = c(0,0), count = FALSE, clustermethod = hclustCBI, 
                                      method=hclust_method, k=k, seed = 249730)$subsetmean)
  
  pam_acf_stab <- mean(clusterboot(as.dist(acf_dist_matrix), B=B, bootmethod = "subset", subtuning = floor(subset*nrow(test)),
                                    noisetuning = c(0,0), count = FALSE, clustermethod = claraCBI, k=k, seed = 249730)$subsetmean)
  hclust_acf_stab <- mean(clusterboot(as.dist(acf_dist_matrix), B=B, bootmethod = "subset", subtuning = floor(subset*nrow(test)),
                                       noisetuning = c(0,0), count = FALSE, clustermethod = hclustCBI, 
                                       method=hclust_method, k=k, seed = 249730)$subsetmean)
  
  #tabelka
  wskazniki_eucl_ref <- c(kmean_eucl_ref_acc, kmean_eucl_conn, kmean_eucl_silh, kmean_eucl_stab,
                          pam_eucl_ref_acc, pam_eucl_conn, pam_eucl_silh, pam_eucl_stab,
                          hclust_eucl_ref_acc, hclust_eucl_conn, hclust_eucl_silh, hclust_eucl_stab,
                          som_ref_acc, som_ref_conn, som_ref_silh, NA)

  wskazniki_acf_ref <- c(pam_acf_ref_acc, pam_acf_conn, pam_acf_silh, pam_acf_stab,
                         hclust_acf_ref_acc, hclust_acf_conn, hclust_acf_silh, hclust_acf_stab)
  
  table_eucl_ref_method <- as.data.frame(matrix(wskazniki_eucl_ref, nrow = 4, ncol = 4, byrow = TRUE))
  colnames(table_eucl_ref_method) <- c("Accuracy", "Connectivity", "Silhouette", "Stability")
  rownames(table_eucl_ref_method) <- c("Kmeans", "PAM", "Hclust", "SOM")
  table_acf_ref_method <- as.data.frame(matrix(wskazniki_acf_ref, nrow = 2, ncol = 4, byrow = TRUE))
  colnames(table_acf_ref_method) <- c("Accuracy", "Connectivity", "Silhouette", "Stability")
  rownames(table_acf_ref_method) <- c("PAM", "Hclust")
  
  #ekstrakcja cech
  train_meas <- get_measure_df(train)
  test_meas <- get_measure_df(test)
  
  train_meas_raw <- train_meas[,1:9]
  train_meas_tsa <- train_meas[,10:13]
  test_meas_raw <- test_meas[,1:9]
  test_meas_tsa <- test_meas[,10:13]
  #usuwamy kolumny z mniejsza liczbą unikalnych wartości niz liczba klastrow k
  train_meas <- train_meas[, unique_values(train_meas, k=k)]
  test_meas <- test_meas[, unique_values(test_meas, k=k)]
  
  train_meas_raw <- train_meas_raw[, unique_values(train_meas_raw, k=k)]
  train_meas_tsa <- train_meas_tsa[, unique_values(train_meas_tsa, k=k)]
  test_meas_raw <- test_meas_raw[, unique_values(test_meas_raw, k=k)]
  test_meas_tsa <- test_meas_tsa[, unique_values(test_meas_tsa, k=k)]
  
  #grupowanie oparte na ekstrakcji cech
  #macierz odległości do grupowania metodami 
  #opartymi na ekstrakcji cech
  eucl_dist_meas_matrix <- dist(test_meas, diag = TRUE, upper = TRUE)
  #wszystkie cechy
  kmean_eucl_char <- kmeans(test_meas, centers = k)
  pam_eucl_char <- pam(eucl_dist_meas_matrix, k = k, metric = "euclidean")
  hclust_eucl_char <- hclust(eucl_dist_meas_matrix, method = hclust_method)
  hclust_eucl_char_clust <- cutree(hclust_eucl_char, k=k)
  train_scaled <- scale(train_meas)
  test_scaled <- scale(test_meas, center = attr(train_scaled, "scaled:center"),
                       scale = attr(train_scaled, "scaled:scale"))
  training_data <- list(measurements = train_scaled,
                        class = as.factor(train_class))
  test_data <- list(measurements = test_scaled, class = as.factor(test_class))
  som_char <- kohonen::supersom(training_data, grid = grid)
  som_pred_char <- predict(som_char, newdata = test_data[1])
  som_tab_char <- table(test_class, som_pred_char$predictions[["class"]])
  som_tab_char_permut <- matchClasses(som_tab_char, method = "exact", verbose = FALSE)
  which.na <- c(attr(na.omit(as.integer(som_pred_char$predictions[["class"]])),"na.action"))
  l <- length(which.na)
  print(which.na)
  print(l)
  
  #ocena metod opartych na ekstrakcji cech (wszystkie cechy)
  #accuracy
  kmean_eucl_char_acc <- compareMatchedClasses(kmean_eucl_char$cluster, test_class, method = "exact", verbose = FALSE)$diag
  pam_eucl_char_acc <- compareMatchedClasses(pam_eucl_char$clustering, test_class, method = "exact", verbose = FALSE)$diag
  hclust_eucl_char_acc <- compareMatchedClasses(hclust_eucl_char_clust, test_class, method = "exact",verbose = FALSE)$diag
  som_char_acc <- sum(diag(som_tab_char[,som_tab_char_permut]))/test_row

  #spójność (connectivity)
  kmean_eucl_conn_char <- connectivity.diss.mx(as.matrix(eucl_dist_meas_matrix), kmean_eucl_char$cluster, neighbour.num = 5)
  pam_eucl_conn_char <- connectivity.diss.mx(as.matrix(eucl_dist_meas_matrix), pam_eucl_char$clustering, neighbour.num = 5)
  hclust_eucl_conn_char <- connectivity.diss.mx(as.matrix(eucl_dist_meas_matrix), hclust_eucl_char_clust, neighbour.num = 5)
  if(l != 0){
    som_char_conn <- connectivity(as.matrix(test_meas)[-which.na,], c(na.omit(as.integer(som_pred_char$predictions[["class"]]))), 
                                 neighbour.num = 5)
  }
  else{
    som_char_conn <- connectivity(as.matrix(test_meas), c(na.omit(as.integer(som_pred_char$predictions[["class"]]))), 
                                 neighbour.num = 5)
  }

  #indeks silhouette
  kmean_eucl_silh_char <- summary(silhouette(kmean_eucl_char$cluster, dist = eucl_dist_meas_matrix))$avg.width
  pam_eucl_silh_char <- summary(silhouette(pam_eucl_char))$avg.width
  hclust_eucl_silh_char <- summary(silhouette(hclust_eucl_char_clust, dist = eucl_dist_meas_matrix))$avg.width
  if(l != 0){
    som_char_silh <- summary(silhouette(na.omit(as.integer(som_pred_char$predictions[["class"]])), 
                                       dist = dist(test_data$measurements[-which.na,], diag = TRUE,upper = TRUE)))$avg.width
  }
  else{
    som_char_silh <- summary(silhouette(na.omit(as.integer(som_pred_char$predictions[["class"]])), 
                                       dist = dist(test_data$measurements, diag = TRUE,upper = TRUE)))$avg.width
  }
  
  #stability
  kmean_char_stab <- mean(clusterboot(test_meas, B=B, bootmethod = "subset", subtuning = floor(subset*nrow(test_meas)),
                                      noisetuning = c(0,0), count = FALSE, clustermethod = kmeansCBI, krange=k, seed = 249730)$subsetmean)
  pam_char_stab <- mean(clusterboot(test_meas, B=B, bootmethod = "subset", subtuning = floor(subset*nrow(test_meas)),
                                    noisetuning = c(0,0), count = FALSE, clustermethod = claraCBI, k=k, seed = 249730)$subsetmean)
  hclust_char_stab <- mean(clusterboot(test_meas, B=B, bootmethod = "subset", subtuning = floor(subset*nrow(test_meas)),
                                       noisetuning = c(0,0), count = FALSE, clustermethod = hclustCBI, 
                                       method=hclust_method, k=k, seed = 249730)$subsetmean)
  print("oceny po ekstrakcji")
  
  #tabelka
  wskazniki_eucl_char <- c(kmean_eucl_char_acc, kmean_eucl_conn_char, kmean_eucl_silh_char, kmean_char_stab,
                          pam_eucl_char_acc, pam_eucl_conn_char, pam_eucl_silh_char, pam_char_stab,
                          hclust_eucl_char_acc, hclust_eucl_conn_char, hclust_eucl_silh_char, hclust_char_stab,
                          som_char_acc, som_char_conn, som_char_silh, NA)

  table_eucl_char_method <- as.data.frame(matrix(wskazniki_eucl_char, nrow = 4, ncol = 4, byrow = TRUE))
  colnames(table_eucl_char_method) <- c("Accuracy", "Connectivity", "Silhouette", "Stability")
  rownames(table_eucl_char_method) <- c("Kmeans", "PAM", "Hclust", "SOM")
  
  #grupowanie oparte na ekstrakcji cech (wybor zmiennych metoda ładunków PCA)
  #macierz odległości do grupowania metodami 
  #opartymi na ekstrakcji cech (PCA)
  #pca
  pca_ind <- pca_select_features(test_meas, min_var_explained = min_var)
  pca_names <- names(test_meas[,pca_ind])
  #
  eucl_dist_pca_matrix <- dist(test_meas[,pca_ind], diag = TRUE, upper = TRUE)
  kmean_eucl_pca <- kmeans(test_meas[,pca_ind], centers = k)
  pam_eucl_pca <- pam(eucl_dist_pca_matrix, k = k, metric = "euclidean")
  hclust_eucl_pca <- hclust(eucl_dist_pca_matrix, method = hclust_method)
  hclust_eucl_pca_clust <- cutree(hclust_eucl_pca, k=k)
  train_scaled <- scale(train_meas[,pca_ind])
  test_scaled <- scale(test_meas[,pca_ind], center = attr(train_scaled, "scaled:center"),
                       scale = attr(train_scaled, "scaled:scale"))
  training_data <- list(measurements = train_scaled,
                        class = as.factor(train_class))
  test_data <- list(measurements = test_scaled, class = as.factor(test_class))
  som_pca <- kohonen::supersom(training_data, grid = grid)
  som_pred_pca <- predict(som_pca, newdata = test_data[1])
  som_tab_pca <- table(test_class, som_pred_pca$predictions[["class"]])
  som_tab_pca_permut <- matchClasses(som_tab_pca, method = "exact", verbose = FALSE)
  which.na <- c(attr(na.omit(as.integer(som_pred_pca$predictions[["class"]])),"na.action"))
  l <- length(which.na)
  print(which.na)
  print(l)
  #ocena metod opartych na ekstrakcji cech (cechy wybrane metodą ładunków PCA)
  #accuracy
  kmean_eucl_pca_acc <- compareMatchedClasses(kmean_eucl_pca$cluster, test_class, method = "exact", verbose = FALSE)$diag
  pam_eucl_pca_acc <- compareMatchedClasses(pam_eucl_pca$clustering, test_class, method = "exact", verbose = FALSE)$diag
  hclust_eucl_pca_acc <- compareMatchedClasses(hclust_eucl_pca_clust, test_class, method = "exact",verbose = FALSE)$diag
  som_pca_acc <- sum(diag(som_tab_pca[,som_tab_pca_permut]))/test_row

  #spójność (connectivity)
  kmean_eucl_conn_pca <- connectivity.diss.mx(as.matrix(eucl_dist_pca_matrix), kmean_eucl_pca$cluster, neighbour.num = 5)
  pam_eucl_conn_pca <- connectivity.diss.mx(as.matrix(eucl_dist_pca_matrix), pam_eucl_pca$clustering, neighbour.num = 5)
  hclust_eucl_conn_pca <- connectivity.diss.mx(as.matrix(eucl_dist_pca_matrix), hclust_eucl_pca_clust, neighbour.num = 5)
  if(l != 0){
    som_pca_conn <- connectivity(as.matrix(test_meas[,pca_ind])[-which.na,], c(na.omit(as.integer(som_pred_pca$predictions[["class"]]))), 
                                neighbour.num = 5)
  }
  else{
    som_pca_conn <- connectivity(as.matrix(test_meas[,pca_ind]), c(na.omit(as.integer(som_pred_pca$predictions[["class"]]))), 
                                 neighbour.num = 5)
  }

  #indeks silhouette
  kmean_eucl_silh_pca <- summary(silhouette(kmean_eucl_pca$cluster, dist = eucl_dist_pca_matrix))$avg.width
  pam_eucl_silh_pca <- summary(silhouette(pam_eucl_pca))$avg.width
  hclust_eucl_silh_pca <- summary(silhouette(hclust_eucl_pca_clust, dist = eucl_dist_pca_matrix))$avg.width
  if(l != 0){
    som_pca_silh <- summary(silhouette(na.omit(as.integer(som_pred_pca$predictions[["class"]])), 
                                      dist = dist(test_data$measurements[-which.na,], diag = TRUE,upper = TRUE)))$avg.width
  }
  else{
    som_pca_silh <- summary(silhouette(na.omit(as.integer(som_pred_pca$predictions[["class"]])), 
                                       dist = dist(test_data$measurements, diag = TRUE,upper = TRUE)))$avg.width
  }
  #stability
  kmean_pca_stab <- mean(clusterboot(test_meas[,pca_ind], B=B, bootmethod = "subset", subtuning = floor(subset*nrow(test_meas)),
                                      noisetuning = c(0,0), count = FALSE, clustermethod = kmeansCBI, krange=k, seed = 249730)$subsetmean)
  pam_pca_stab <- mean(clusterboot(test_meas[,pca_ind], B=B, bootmethod = "subset", subtuning = floor(subset*nrow(test_meas)),
                                    noisetuning = c(0,0), count = FALSE, clustermethod = claraCBI, k=k, seed = 249730)$subsetmean)
  hclust_pca_stab <- mean(clusterboot(test_meas[,pca_ind], B=B, bootmethod = "subset", subtuning = floor(subset*nrow(test_meas)),
                                       noisetuning = c(0,0), count = FALSE, clustermethod = hclustCBI, 
                                       method=hclust_method, k=k, seed = 249730)$subsetmean)
  
  #tabelka
  wskazniki_eucl_pca <- c(kmean_eucl_pca_acc, kmean_eucl_conn_pca, kmean_eucl_silh_pca, kmean_pca_stab,
                           pam_eucl_pca_acc, pam_eucl_conn_pca, pam_eucl_silh_pca, pam_pca_stab,
                           hclust_eucl_pca_acc, hclust_eucl_conn_pca, hclust_eucl_silh_pca, hclust_pca_stab,
                          som_pca_acc, som_pca_conn, som_pca_silh, NA)

  table_eucl_pca_method <- as.data.frame(matrix(wskazniki_eucl_pca, nrow = 4, ncol = 4, byrow = TRUE))
  colnames(table_eucl_pca_method) <- c("Accuracy", "Connectivity", "Silhouette", "Stability")
  rownames(table_eucl_pca_method) <- c("Kmeans", "PAM", "Hclust", "SOM")
  
  #grupowanie po ekstrakcji cech wylacznie na cechach z danych RAW
  eucl_dist_raw_matrix <- dist(test_meas_raw, diag = TRUE, upper = TRUE)

  kmean_eucl_raw <- kmeans(test_meas_raw, centers = k)
  pam_eucl_raw <- pam(eucl_dist_raw_matrix, k = k, metric = "euclidean")
  hclust_eucl_raw <- hclust(eucl_dist_raw_matrix, method = hclust_method)
  hclust_eucl_raw_clust <- cutree(hclust_eucl_raw, k=k)
  train_scaled <- scale(train_meas_raw)
  test_scaled <- scale(test_meas_raw, center = attr(train_scaled, "scaled:center"),
                       scale = attr(train_scaled, "scaled:scale"))
  training_data <- list(measurements = train_scaled,
                        class = as.factor(train_class))
  test_data <- list(measurements = test_scaled, class = as.factor(test_class))
  som_raw <- kohonen::supersom(training_data, grid = grid)
  som_pred_raw <- predict(som_raw, newdata = test_data[1])
  som_tab_raw <- table(test_class, som_pred_raw$predictions[["class"]])
  som_tab_raw_permut <- matchClasses(som_tab_raw, method = "exact", verbose = FALSE)
  which.na <- c(attr(na.omit(as.integer(som_pred_raw$predictions[["class"]])),"na.action"))
  l <- length(which.na)
  print(which.na)
  print(l)
  #ocena metod opartych na ekstrakcji cech (cechy RAW)
  #accuracy
  kmean_eucl_raw_acc <- compareMatchedClasses(kmean_eucl_raw$cluster, test_class, method = "exact", verbose = FALSE)$diag
  pam_eucl_raw_acc <- compareMatchedClasses(pam_eucl_raw$clustering, test_class, method = "exact", verbose = FALSE)$diag
  hclust_eucl_raw_acc <- compareMatchedClasses(hclust_eucl_raw_clust, test_class, method = "exact",verbose = FALSE)$diag
  som_raw_acc <- sum(diag(som_tab_raw[,som_tab_raw_permut]))/test_row

  #spójność (connectivity)
  kmean_eucl_conn_raw <- connectivity.diss.mx(as.matrix(eucl_dist_raw_matrix), kmean_eucl_raw$cluster, neighbour.num = 5)
  pam_eucl_conn_raw <- connectivity.diss.mx(as.matrix(eucl_dist_raw_matrix), pam_eucl_raw$clustering, neighbour.num = 5)
  hclust_eucl_conn_raw <- connectivity.diss.mx(as.matrix(eucl_dist_raw_matrix), hclust_eucl_raw_clust, neighbour.num = 5)
  if(l != 0){
    som_raw_conn <- connectivity(as.matrix(test_meas_raw)[-which.na,], c(na.omit(as.integer(som_pred_raw$predictions[["class"]]))), 
                                 neighbour.num = 5)
  }
  else{
    som_raw_conn <- connectivity(as.matrix(test_meas_raw), c(na.omit(as.integer(som_pred_raw$predictions[["class"]]))), 
                                 neighbour.num = 5)
  }

  #indeks silhouette
  kmean_eucl_silh_raw <- summary(silhouette(kmean_eucl_raw$cluster, dist = eucl_dist_raw_matrix))$avg.width
  pam_eucl_silh_raw <- summary(silhouette(pam_eucl_raw))$avg.width
  hclust_eucl_silh_raw <- summary(silhouette(hclust_eucl_raw_clust, dist = eucl_dist_raw_matrix))$avg.width
  if(l != 0){
    som_raw_silh <- summary(silhouette(na.omit(as.integer(som_pred_raw$predictions[["class"]])), 
                                       dist = dist(test_data$measurements[-which.na,], diag = TRUE,upper = TRUE)))$avg.width
  }
  else{
    som_raw_silh <- summary(silhouette(na.omit(as.integer(som_pred_raw$predictions[["class"]])), 
                                       dist = dist(test_data$measurements, diag = TRUE,upper = TRUE)))$avg.width
  }
  
  #stability
  kmean_raw_stab <- mean(clusterboot(test_meas_raw, B=B, bootmethod = "subset", subtuning = floor(subset*nrow(test_meas)),
                                     noisetuning = c(0,0), count = FALSE, clustermethod = kmeansCBI, krange=k, seed = 249730)$subsetmean)
  pam_raw_stab <- mean(clusterboot(test_meas_raw, B=B, bootmethod = "subset", subtuning = floor(subset*nrow(test_meas)),
                                   noisetuning = c(0,0), count = FALSE, clustermethod = claraCBI, k=k, seed = 249730)$subsetmean)
  hclust_raw_stab <- mean(clusterboot(test_meas_raw, B=B, bootmethod = "subset", subtuning = floor(subset*nrow(test_meas)),
                                      noisetuning = c(0,0), count = FALSE, clustermethod = hclustCBI, 
                                      method=hclust_method, k=k, seed = 249730)$subsetmean)
  
  #tabelka
  wskazniki_eucl_raw <- c(kmean_eucl_raw_acc, kmean_eucl_conn_raw, kmean_eucl_silh_raw, kmean_raw_stab,
                          pam_eucl_raw_acc, pam_eucl_conn_raw, pam_eucl_silh_raw, pam_raw_stab,
                          hclust_eucl_raw_acc, hclust_eucl_conn_raw, hclust_eucl_silh_raw, hclust_raw_stab,
                          som_raw_acc, som_raw_conn, som_raw_silh, NA)

  table_eucl_raw_method <- as.data.frame(matrix(wskazniki_eucl_raw, nrow = 4, ncol = 4, byrow = TRUE))
  colnames(table_eucl_raw_method) <- c("Accuracy", "Connectivity", "Silhouette", "Stability")
  rownames(table_eucl_raw_method) <- c("Kmeans", "PAM", "Hclust", "SOM")
  
  #grupowanie po ekstrakcji cech wylacznie na cechach z danych TSA
  eucl_dist_tsa_matrix <- dist(test_meas_tsa, diag = TRUE, upper = TRUE)
  
  kmean_eucl_tsa <- kmeans(test_meas_tsa, centers = k)
  pam_eucl_tsa <- pam(eucl_dist_tsa_matrix, k = k, metric = "euclidean")
  hclust_eucl_tsa <- hclust(eucl_dist_tsa_matrix, method = hclust_method)
  hclust_eucl_tsa_clust <- cutree(hclust_eucl_tsa, k=k)
  train_scaled <- scale(train_meas_tsa)
  test_scaled <- scale(test_meas_tsa, center = attr(train_scaled, "scaled:center"),
                       scale = attr(train_scaled, "scaled:scale"))
  training_data <- list(measurements = train_scaled,
                        class = as.factor(train_class))
  test_data <- list(measurements = test_scaled, class = as.factor(test_class))
  som_tsa <- kohonen::supersom(training_data, grid = grid)
  som_pred_tsa <- predict(som_tsa, newdata = test_data[1])
  som_tab_tsa <- table(test_class, som_pred_tsa$predictions[["class"]])
  som_tab_tsa_permut <- matchClasses(som_tab_tsa, method = "exact", verbose = FALSE)
  which.na <- c(attr(na.omit(as.integer(som_pred_tsa$predictions[["class"]])),"na.action"))
  l <- length(which.na)
  print(which.na)
  print(l)
  #ocena metod opartych na ekstrakcji cech (cechy TSA)
  #accuracy
  kmean_eucl_tsa_acc <- compareMatchedClasses(kmean_eucl_tsa$cluster, test_class, method = "exact", verbose = FALSE)$diag
  pam_eucl_tsa_acc <- compareMatchedClasses(pam_eucl_tsa$clustering, test_class, method = "exact", verbose = FALSE)$diag
  hclust_eucl_tsa_acc <- compareMatchedClasses(hclust_eucl_tsa_clust, test_class, method = "exact",verbose = FALSE)$diag
  som_tsa_acc <- sum(diag(som_tab_tsa[,som_tab_tsa_permut]))/test_row

  #spójność (connectivity)
  kmean_eucl_conn_tsa <- connectivity.diss.mx(as.matrix(eucl_dist_tsa_matrix), kmean_eucl_tsa$cluster, neighbour.num = 5)
  pam_eucl_conn_tsa <- connectivity.diss.mx(as.matrix(eucl_dist_tsa_matrix), pam_eucl_tsa$clustering, neighbour.num = 5)
  hclust_eucl_conn_tsa <- connectivity.diss.mx(as.matrix(eucl_dist_tsa_matrix), hclust_eucl_tsa_clust, neighbour.num = 5)
  if(l != 0){
    som_tsa_conn <- connectivity(as.matrix(test_meas_tsa)[-which.na,], c(na.omit(as.integer(som_pred_tsa$predictions[["class"]]))), 
                                 neighbour.num = 5)
  }
  else{
    som_tsa_conn <- connectivity(as.matrix(test_meas_tsa), c(na.omit(as.integer(som_pred_tsa$predictions[["class"]]))), 
                                 neighbour.num = 5)
  }

  #indeks silhouette
  kmean_eucl_silh_tsa <- summary(silhouette(kmean_eucl_tsa$cluster, dist = eucl_dist_tsa_matrix))$avg.width
  pam_eucl_silh_tsa <- summary(silhouette(pam_eucl_tsa))$avg.width
  hclust_eucl_silh_tsa <- summary(silhouette(hclust_eucl_tsa_clust, dist = eucl_dist_tsa_matrix))$avg.width
  if(l != 0){
    som_tsa_silh <- summary(silhouette(na.omit(as.integer(som_pred_tsa$predictions[["class"]])), 
                                       dist = dist(test_data$measurements[-which.na,], diag = TRUE,upper = TRUE)))$avg.width
  }
  else{
    som_tsa_silh <- summary(silhouette(na.omit(as.integer(som_pred_tsa$predictions[["class"]])), 
                                       dist = dist(test_data$measurements, diag = TRUE,upper = TRUE)))$avg.width
  }
  
  #stability
  kmean_tsa_stab <- mean(clusterboot(test_meas_tsa, B=B, bootmethod = "subset", subtuning = floor(subset*nrow(test_meas)),
                                     noisetuning = c(0,0), count = FALSE, clustermethod = kmeansCBI, krange=k, seed = 249730)$subsetmean)
  pam_tsa_stab <- mean(clusterboot(test_meas_tsa, B=B, bootmethod = "subset", subtuning = floor(subset*nrow(test_meas)),
                                   noisetuning = c(0,0), count = FALSE, clustermethod = claraCBI, k=k, seed = 249730)$subsetmean)
  hclust_tsa_stab <- mean(clusterboot(test_meas_tsa, B=B, bootmethod = "subset", subtuning = floor(subset*nrow(test_meas)),
                                      noisetuning = c(0,0), count = FALSE, clustermethod = hclustCBI, 
                                      method=hclust_method, k=k, seed = 249730)$subsetmean)
  
  #tabelka
  wskazniki_eucl_tsa <- c(kmean_eucl_tsa_acc, kmean_eucl_conn_tsa, kmean_eucl_silh_tsa, kmean_tsa_stab,
                          pam_eucl_tsa_acc, pam_eucl_conn_tsa, pam_eucl_silh_tsa, pam_tsa_stab,
                          hclust_eucl_tsa_acc, hclust_eucl_conn_tsa, hclust_eucl_silh_tsa, hclust_tsa_stab,
                          som_tsa_acc, som_tsa_conn, som_tsa_silh, NA)

  table_eucl_tsa_method <- as.data.frame(matrix(wskazniki_eucl_tsa, nrow = 4, ncol = 4, byrow = TRUE))
  colnames(table_eucl_tsa_method) <- c("Accuracy", "Connectivity", "Silhouette", "Stability")
  rownames(table_eucl_tsa_method) <- c("Kmeans", "PAM", "Hclust", "SOM")
  
  #grupowanie na podstawie wyboru cech metodą forward search
  #(nadzorowana) funkcja stepclass
  #wybór cech
  sc_object <- stepclass(train_meas, train_class, "sknn", direction = "forward",
                         criterion = "CR", improvement = 0.005, kn = 1)
  names_sc <- sc_object$model[,"name"]
  if (length(sc_object$model[,"nr"]) == 1){
    print("Stosując nadzorowaną metodę wyboru zmiennych na wszystkich ekstrachowanych cechach\n
          wybrano tylko jedną zmienną. Stosujemy metodę biorąc pod uwagę wyłącznie\n
          zmienne typu ciągłego.")
    sc_object <- stepclass(train_meas[,-1], train_class, "sknn", direction = "forward",
                           criterion = "CR", improvement = 0.005, kn = 1)
    names_sc <- sc_object$model[,"name"]
  }
  train_meas_sc <- train_meas[,c(sc_object$model[,"nr"])]
  test_meas_sc <- test_meas[,c(sc_object$model[,"nr"])]
  #grupowanie po ekstrakcji cech wylacznie na cechach po stepclass
  eucl_dist_sc_matrix <- dist(test_meas_sc, diag = TRUE, upper = TRUE)
  #
  kmean_eucl_sc <- kmeans(test_meas_sc, centers = k)
  pam_eucl_sc <- pam(eucl_dist_sc_matrix, k = k, metric = "euclidean")
  hclust_eucl_sc <- hclust(eucl_dist_sc_matrix, method = hclust_method)
  hclust_eucl_sc_clust <- cutree(hclust_eucl_sc, k=k)
  train_scaled <- scale(train_meas_sc)
  test_scaled <- scale(test_meas_sc, center = attr(train_scaled, "scaled:center"),
                       scale = attr(train_scaled, "scaled:scale"))
  training_data <- list(measurements = train_scaled,
                        class = as.factor(train_class))
  test_data <- list(measurements = test_scaled, class = as.factor(test_class))
  som_sc <- kohonen::supersom(training_data, grid = grid)
  som_pred_sc <- predict(som_sc, newdata = test_data[1])
  som_tab_sc <- table(test_class, som_pred_sc$predictions[["class"]])
  som_tab_sc_permut <- matchClasses(som_tab_sc, method = "exact", verbose = FALSE)
  which.na <- c(attr(na.omit(as.integer(som_pred_sc$predictions[["class"]])),"na.action"))
  l <- length(which.na)
  print(which.na)
  print(l)
  #ocena metod opartych na ekstrakcji cech (stepclass)
  #accuracy
  kmean_eucl_sc_acc <- compareMatchedClasses(kmean_eucl_sc$cluster, test_class, method = "exact", verbose = FALSE)$diag
  pam_eucl_sc_acc <- compareMatchedClasses(pam_eucl_sc$clustering, test_class, method = "exact", verbose = FALSE)$diag
  hclust_eucl_sc_acc <- compareMatchedClasses(hclust_eucl_sc_clust, test_class, method = "exact",verbose = FALSE)$diag
  som_sc_acc <- sum(diag(som_tab_sc[,som_tab_sc_permut]))/test_row
  
  #spójność (connectivity)
  kmean_eucl_conn_sc <- connectivity.diss.mx(as.matrix(eucl_dist_sc_matrix), kmean_eucl_sc$cluster, neighbour.num = 5)
  pam_eucl_conn_sc <- connectivity.diss.mx(as.matrix(eucl_dist_sc_matrix), pam_eucl_sc$clustering, neighbour.num = 5)
  hclust_eucl_conn_sc <- connectivity.diss.mx(as.matrix(eucl_dist_sc_matrix), hclust_eucl_sc_clust, neighbour.num = 5)
  if(l != 0){
    som_sc_conn <- connectivity(as.matrix(test_meas_sc)[-which.na,], c(na.omit(as.integer(som_pred_sc$predictions[["class"]]))), 
                                 neighbour.num = 5)
  }
  else{
    som_sc_conn <- connectivity(as.matrix(test_meas_sc), c(na.omit(as.integer(som_pred_sc$predictions[["class"]]))), 
                                 neighbour.num = 5)
  }
  
  #indeks silhouette
  kmean_eucl_silh_sc <- summary(silhouette(kmean_eucl_sc$cluster, dist = eucl_dist_sc_matrix))$avg.width
  pam_eucl_silh_sc <- summary(silhouette(pam_eucl_sc))$avg.width
  hclust_eucl_silh_sc <- summary(silhouette(hclust_eucl_sc_clust, dist = eucl_dist_sc_matrix))$avg.width
  if(l != 0){
    som_sc_silh <- summary(silhouette(na.omit(as.integer(som_pred_sc$predictions[["class"]])), 
                                       dist = dist(test_data$measurements[-which.na,], diag = TRUE,upper = TRUE)))$avg.width
  }
  else{
    som_sc_silh <- summary(silhouette(na.omit(as.integer(som_pred_sc$predictions[["class"]])), 
                                       dist = dist(test_data$measurements, diag = TRUE,upper = TRUE)))$avg.width
  }
  
  #stability
  kmean_sc_stab <- mean(clusterboot(test_meas_sc, B=B, bootmethod = "subset", subtuning = floor(subset*nrow(test_meas)),
                                     noisetuning = c(0,0), count = FALSE, clustermethod = kmeansCBI, krange=k, seed = 249730)$subsetmean)
  pam_sc_stab <- mean(clusterboot(test_meas_sc, B=B, bootmethod = "subset", subtuning = floor(subset*nrow(test_meas)),
                                   noisetuning = c(0,0), count = FALSE, clustermethod = claraCBI, k=k, seed = 249730)$subsetmean)
  hclust_sc_stab <- mean(clusterboot(test_meas_sc, B=B, bootmethod = "subset", subtuning = floor(subset*nrow(test_meas)),
                                      noisetuning = c(0,0), count = FALSE, clustermethod = hclustCBI, 
                                      method=hclust_method, k=k, seed = 249730)$subsetmean)
  
  #tabelka
  wskazniki_eucl_sc <- c(kmean_eucl_sc_acc, kmean_eucl_conn_sc, kmean_eucl_silh_sc, kmean_sc_stab,
                          pam_eucl_sc_acc, pam_eucl_conn_sc, pam_eucl_silh_sc, pam_sc_stab,
                          hclust_eucl_sc_acc, hclust_eucl_conn_sc, hclust_eucl_silh_sc, hclust_sc_stab,
                          som_sc_acc, som_sc_conn, som_sc_silh, NA)
  table_eucl_sc_method <- as.data.frame(matrix(wskazniki_eucl_sc, nrow = 4, ncol = 4, byrow = TRUE))
  colnames(table_eucl_sc_method) <- c("Accuracy", "Connectivity", "Silhouette", "Stability")
  rownames(table_eucl_sc_method) <- c("Kmeans", "PAM", "Hclust", "SOM")
  
  #grupowanie po metodzie forward search (unsupervised)
  #wybór cech
  fs_object <- forward(test_meas, krange = k, 
                       itter=100, discrete = discrete, track = FALSE)
  names_fs <- fs_object$picked
  if (length(fs_object$picked) == 1){
    print("Stosując nienadzorowaną metodę wyboru zmiennych na wszystkich ekstrachowanych cechach\n
          wybrano tylko jedną zmienną. Stosujemy metodę biorąc pod uwagę wyłącznie\n
          zmienne typu ciągłego.")
    fs_object <- forward(test_meas, krange = k, 
                         itter=100, discrete = FALSE, track = FALSE)
    names_fs <- fs_object$picked
  }
  train_meas_fs <- train_meas[,-fs_object$to_del_index]
  test_meas_fs <- test_meas[,-fs_object$to_del_index]
  #grupowanie po ekstrakcji cech wylacznie na cechach po forward search
  eucl_dist_fs_matrix <- dist(test_meas_fs, diag = TRUE, upper = TRUE)
  #
  kmean_eucl_fs <- kmeans(test_meas_fs, centers = k)
  pam_eucl_fs <- pam(eucl_dist_fs_matrix, k = k, metric = "euclidean")
  hclust_eucl_fs <- hclust(eucl_dist_fs_matrix, method = hclust_method)
  hclust_eucl_fs_clust <- cutree(hclust_eucl_fs, k=k)
  train_scaled <- scale(train_meas_fs)
  test_scaled <- scale(test_meas_fs, center = attr(train_scaled, "scaled:center"),
                       scale = attr(train_scaled, "scaled:scale"))
  training_data <- list(measurements = train_scaled,
                        class = as.factor(train_class))
  test_data <- list(measurements = test_scaled, class = as.factor(test_class))
  som_fs <- kohonen::supersom(training_data, grid = grid)
  som_pred_fs <- predict(som_fs, newdata = test_data[1])
  som_tab_fs <- table(test_class, som_pred_fs$predictions[["class"]])
  som_tab_fs_permut <- matchClasses(som_tab_fs, method = "exact", verbose = FALSE)
  which.na <- c(attr(na.omit(as.integer(som_pred_fs$predictions[["class"]])),"na.action"))
  l <- length(which.na)
  print(which.na)
  print(l)
  #ocena metod opartych na ekstrakcji cech (unsupervised forward search)
  #accuracy
  kmean_eucl_fs_acc <- compareMatchedClasses(kmean_eucl_fs$cluster, test_class, method = "exact", verbose = FALSE)$diag
  pam_eucl_fs_acc <- compareMatchedClasses(pam_eucl_fs$clustering, test_class, method = "exact", verbose = FALSE)$diag
  hclust_eucl_fs_acc <- compareMatchedClasses(hclust_eucl_fs_clust, test_class, method = "exact",verbose = FALSE)$diag
  som_fs_acc <- sum(diag(som_tab_fs[,som_tab_fs_permut]))/test_row
  
  #spójność (connectivity)
  kmean_eucl_conn_fs <- connectivity.diss.mx(as.matrix(eucl_dist_fs_matrix), kmean_eucl_fs$cluster, neighbour.num = 5)
  pam_eucl_conn_fs <- connectivity.diss.mx(as.matrix(eucl_dist_fs_matrix), pam_eucl_fs$clustering, neighbour.num = 5)
  hclust_eucl_conn_fs <- connectivity.diss.mx(as.matrix(eucl_dist_fs_matrix), hclust_eucl_fs_clust, neighbour.num = 5)
  if(l != 0){
    som_fs_conn <- connectivity(as.matrix(test_meas_fs)[-which.na,], c(na.omit(as.integer(som_pred_fs$predictions[["class"]]))), 
                                neighbour.num = 5)
  }
  else{
    som_fs_conn <- connectivity(as.matrix(test_meas_fs), c(na.omit(as.integer(som_pred_fs$predictions[["class"]]))), 
                                neighbour.num = 5)
  }
  
  #indeks silhouette
  kmean_eucl_silh_fs <- summary(silhouette(kmean_eucl_fs$cluster, dist = eucl_dist_fs_matrix))$avg.width
  pam_eucl_silh_fs <- summary(silhouette(pam_eucl_fs))$avg.width
  hclust_eucl_silh_fs <- summary(silhouette(hclust_eucl_fs_clust, dist = eucl_dist_fs_matrix))$avg.width
  if(l != 0){
    som_fs_silh <- summary(silhouette(na.omit(as.integer(som_pred_fs$predictions[["class"]])), 
                                      dist = dist(test_data$measurements[-which.na,], diag = TRUE,upper = TRUE)))$avg.width
  }
  else{
    som_fs_silh <- summary(silhouette(na.omit(as.integer(som_pred_fs$predictions[["class"]])), 
                                      dist = dist(test_data$measurements, diag = TRUE,upper = TRUE)))$avg.width
  }
  
  #stability
  kmean_fs_stab <- mean(clusterboot(test_meas_fs, B=B, bootmethod = "subset", subtuning = floor(subset*nrow(test_meas)),
                                    noisetuning = c(0,0), count = FALSE, clustermethod = kmeansCBI, krange=k, seed = 249730)$subsetmean)
  pam_fs_stab <- mean(clusterboot(test_meas_fs, B=B, bootmethod = "subset", subtuning = floor(subset*nrow(test_meas)),
                                  noisetuning = c(0,0), count = FALSE, clustermethod = claraCBI, k=k, seed = 249730)$subsetmean)
  hclust_fs_stab <- mean(clusterboot(test_meas_fs, B=B, bootmethod = "subset", subtuning = floor(subset*nrow(test_meas)),
                                     noisetuning = c(0,0), count = FALSE, clustermethod = hclustCBI, 
                                     method=hclust_method, k=k, seed = 249730)$subsetmean)
  
  #tabelka
  wskazniki_eucl_fs <- c(kmean_eucl_fs_acc, kmean_eucl_conn_fs, kmean_eucl_silh_fs, kmean_fs_stab,
                         pam_eucl_fs_acc, pam_eucl_conn_fs, pam_eucl_silh_fs, pam_fs_stab,
                         hclust_eucl_fs_acc, hclust_eucl_conn_fs, hclust_eucl_silh_fs, hclust_fs_stab,
                         som_fs_acc, som_fs_conn, som_fs_silh, NA)
  table_eucl_fs_method <- as.data.frame(matrix(wskazniki_eucl_fs, nrow = 4, ncol = 4, byrow = TRUE))
  colnames(table_eucl_fs_method) <- c("Accuracy", "Connectivity", "Silhouette", "Stability")
  rownames(table_eucl_fs_method) <- c("Kmeans", "PAM", "Hclust", "SOM")
  
  #grupowanie stosując składowe głowne z metody PCA
  train_comp_model <- prcomp(train_meas, scale. = TRUE)
  test_comp_model <- prcomp(test_meas, scale. = TRUE)
  var_explained <- cumsum(test_comp_model$sdev^2) / sum(test_comp_model$sdev^2)
  n_components <- sum(var_explained >= min_var)
  n_components <- ncol(test_meas) + 1 - n_components
  train_meas_comp <- train_comp_model$x[,1:n_components]
  test_meas_comp <- test_comp_model$x[,1:n_components]
  #grupowanie po ekstrakcji cech na składowych głownych (pca)
  eucl_dist_comp_matrix <- dist(test_meas_comp, diag = TRUE, upper = TRUE)
  #
  kmean_eucl_comp <- kmeans(test_meas_comp, centers = k)
  pam_eucl_comp <- pam(eucl_dist_comp_matrix, k = k, metric = "euclidean")
  hclust_eucl_comp <- hclust(eucl_dist_comp_matrix, method = hclust_method)
  hclust_eucl_comp_clust <- cutree(hclust_eucl_comp, k=k)
  train_scaled <- scale(train_meas_comp)
  test_scaled <- scale(test_meas_comp, center = attr(train_scaled, "scaled:center"),
                       scale = attr(train_scaled, "scaled:scale"))
  training_data <- list(measurements = train_scaled,
                        class = as.factor(train_class))
  test_data <- list(measurements = test_scaled, class = as.factor(test_class))
  som_comp <- kohonen::supersom(training_data, grid = grid)
  som_pred_comp <- predict(som_comp, newdata = test_data[1])
  som_tab_comp <- table(test_class, som_pred_comp$predictions[["class"]])
  som_tab_comp_permut <- matchClasses(som_tab_comp, method = "exact", verbose = FALSE)
  which.na <- c(attr(na.omit(as.integer(som_pred_comp$predictions[["class"]])),"na.action"))
  l <- length(which.na)
  print(which.na)
  print(l)
  #ocena metod opartych na składowych głównych PCA
  #accuracy
  kmean_eucl_comp_acc <- compareMatchedClasses(kmean_eucl_comp$cluster, test_class, method = "exact", verbose = FALSE)$diag
  pam_eucl_comp_acc <- compareMatchedClasses(pam_eucl_comp$clustering, test_class, method = "exact", verbose = FALSE)$diag
  hclust_eucl_comp_acc <- compareMatchedClasses(hclust_eucl_comp_clust, test_class, method = "exact",verbose = FALSE)$diag
  som_comp_acc <- sum(diag(som_tab_comp[,som_tab_comp_permut]))/test_row
  
  #spójność (connectivity)
  kmean_eucl_conn_comp <- connectivity.diss.mx(as.matrix(eucl_dist_comp_matrix), kmean_eucl_comp$cluster, neighbour.num = 5)
  pam_eucl_conn_comp <- connectivity.diss.mx(as.matrix(eucl_dist_comp_matrix), pam_eucl_comp$clustering, neighbour.num = 5)
  hclust_eucl_conn_comp <- connectivity.diss.mx(as.matrix(eucl_dist_comp_matrix), hclust_eucl_comp_clust, neighbour.num = 5)
  if(l != 0){
    som_comp_conn <- connectivity(as.matrix(test_meas_comp)[-which.na,], c(na.omit(as.integer(som_pred_comp$predictions[["class"]]))), 
                                neighbour.num = 5)
  }
  else{
    som_comp_conn <- connectivity(as.matrix(test_meas_comp), c(na.omit(as.integer(som_pred_comp$predictions[["class"]]))), 
                                neighbour.num = 5)
  }
  
  #indeks silhouette
  kmean_eucl_silh_comp <- summary(silhouette(kmean_eucl_comp$cluster, dist = eucl_dist_comp_matrix))$avg.width
  pam_eucl_silh_comp <- summary(silhouette(pam_eucl_comp))$avg.width
  hclust_eucl_silh_comp <- summary(silhouette(hclust_eucl_comp_clust, dist = eucl_dist_comp_matrix))$avg.width
  if(l != 0){
    som_comp_silh <- summary(silhouette(na.omit(as.integer(som_pred_comp$predictions[["class"]])), 
                                      dist = dist(test_data$measurements[-which.na,], diag = TRUE,upper = TRUE)))$avg.width
  }
  else{
    som_comp_silh <- summary(silhouette(na.omit(as.integer(som_pred_comp$predictions[["class"]])), 
                                      dist = dist(test_data$measurements, diag = TRUE,upper = TRUE)))$avg.width
  }
  
  #stability
  kmean_comp_stab <- mean(clusterboot(test_meas_comp, B=B, bootmethod = "subset", subtuning = floor(subset*nrow(test_meas)),
                                    noisetuning = c(0,0), count = FALSE, clustermethod = kmeansCBI, krange=k, seed = 249730)$subsetmean)
  pam_comp_stab <- mean(clusterboot(test_meas_comp, B=B, bootmethod = "subset", subtuning = floor(subset*nrow(test_meas)),
                                  noisetuning = c(0,0), count = FALSE, clustermethod = claraCBI, k=k, seed = 249730)$subsetmean)
  hclust_comp_stab <- mean(clusterboot(test_meas_comp, B=B, bootmethod = "subset", subtuning = floor(subset*nrow(test_meas)),
                                     noisetuning = c(0,0), count = FALSE, clustermethod = hclustCBI, 
                                     method=hclust_method, k=k, seed = 249730)$subsetmean)
  
  #tabelka
  wskazniki_eucl_comp <- c(kmean_eucl_comp_acc, kmean_eucl_conn_comp, kmean_eucl_silh_comp, kmean_comp_stab,
                         pam_eucl_comp_acc, pam_eucl_conn_comp, pam_eucl_silh_comp, pam_comp_stab,
                         hclust_eucl_comp_acc, hclust_eucl_conn_comp, hclust_eucl_silh_comp, hclust_comp_stab,
                         som_comp_acc, som_comp_conn, som_comp_silh, NA)
  table_eucl_comp_method <- as.data.frame(matrix(wskazniki_eucl_comp, nrow = 4, ncol = 4, byrow = TRUE))
  colnames(table_eucl_comp_method) <- c("Accuracy", "Connectivity", "Silhouette", "Stability")
  rownames(table_eucl_comp_method) <- c("Kmeans", "PAM", "Hclust", "SOM")
  
  
  return(list(ref=list(eucl_ref_method = table_eucl_ref_method, acf_ref_method = table_acf_ref_method),
              char = table_eucl_char_method,
              pca = list(table=table_eucl_pca_method, pca_names=pca_names),
              raw = table_eucl_raw_method,
              tsa = table_eucl_tsa_method,
              sc = list(table=table_eucl_sc_method, names_sc = names_sc),
              fs = list(table=table_eucl_fs_method, names_fs = names_fs),
              comp = list(table=table_eucl_comp_method, n_components = n_components)))
}

