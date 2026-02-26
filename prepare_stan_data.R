
prepare_stan_data <- function(f = as.formula(ch_stunt ~ age_of_child + sex_of_child + breastfeeding +  maternal_age + 
                                               maternal_education + wealth_index),
                              y = y, sf_joined = sf_joined, 
                              coord_m = coord_m, 
                              adminGroup_vars = c("admi_NAME"), 
                              cluster_vars_idz = c("DHSCLUST", "Precipitation_2020"),
                              urban_rural_index12 = as.numeric(unique(y[, c("DHSCLUST","URBAN_RURA")])$URBAN_RURA)) {
  library(dplyr)
  library(sf)
  library(geodist)
  
  ################################################################################
  ## Create layer coordinates m & cluster coordinate points J
  ################################################################################
  
  coord_m <- data.frame(coord_m)
  coord_J <- data.frame(coord_J)
  names(coord_m) <- names(coord_J)[2:3]
  
  d <- geodist::geodist(coord_J[,2:3], coord_m, measure = "geodesic") / 1000 # in KM
  
  ################################################################################
  ## Function to create W matrix
  ################################################################################
  fnc_w <- function(cluster_J, layer_m) {
    W_Jxm <- table(cluster_J, layer_m)
    J <- nrow(W_Jxm)
    m <- ncol(W_Jxm)
    W <- matrix(0, nrow = J, ncol = J)
    for (i in 1:m) {
      ck <- W_Jxm
      dimnames(ck)[[1]] <- 1:dim(ck)[1]
      number <- as.numeric(names(ck[, i][ck[, i] == 1]))
      W[number, number] <- 1
    }
    return(W)
  }
  W <- fnc_w(cluster_J = coord_J[,1], layer_m = sf_joined[[adminGroup_vars]])
  
  ################################################################################
  ## Prepare cluster-level data
  ################################################################################
  z <- unique(y[, cluster_vars_idz])
  
  ################################################################################
  ## Create IMAT matrix - initial
  ################################################################################
  Ind <- data.frame(clst=y[, cluster_vars_idz[1]])
  for(i in z[,1]){
    Ind[[paste0("clst", i)]] <- as.numeric(Ind$clst == i)
  }
  IMAT <- as.matrix(Ind[, -1]) # N x J
  
  ################################################################################
  ## Urban-rural index
  ################################################################################
  ur <- urban_rural_index12
  
  ################################################################################
  ## MW matrix - inital
  ################################################################################
  MW <- (1 - z[,2] %*% solve(t(z[,2]) %*% z[,2]) %*% t(z[,2])) %*% 
    W %*% 
    (1 - z[,2] %*% solve(t(z[,2]) %*% z[,2]) %*% t(z[,2]))
  
  ################################################################################
  ## Prepare QQ and eigen decomposition - initial
  ################################################################################
  x <- model.matrix(f, y)
  phi <- max(d)
  QQ <- t(MW %*% (1 - (d/phi)^2)^2) %*% (diag(W)*1 - W) %*% (MW %*% (1 - (d/phi)^2)^2)
  eigen_QQ <- eigen(QQ)
  
  ################################################################################
  ## Initial zstar
  ################################################################################
  J <- length(unique(sf_joined[[cluster_vars_idz[1]]]))
  M <- length(unique(sf_joined[[adminGroup_vars[1]]]))
  M_JxJ <- eigen(MW)$vectors
  D_JxM <- d
  Dstar <- d
  Sigma_diag <- abs(diag(eigen_QQ$vectors))
  
  sigma_r <- 1 / rgamma(100, 10, scale = 1/25)
  sigma_u <- 1 / rgamma(100, 20, scale = 1/20)
  
  for(m in 1:M){
    for(j in 1:J){
      if(ur[j] == 1){ # rural
        sigma_r_val <- sd(rnorm(length(sigma_r), mean = D_JxM[j, m], sd = sigma_r))
        Dstar[j, m] <- rnorm(1, D_JxM[j, m], sigma_r_val)
      } else { # urban
        sigma_u_val <- sd(rnorm(length(sigma_u), mean = D_JxM[j, m], sd = sigma_u))
        Dstar[j, m] <- rnorm(1, D_JxM[j, m], sigma_u_val)
      }
    }
  }
  
  Psi <- M_JxJ %*% (1 - (Dstar / phi)^2)^2
  Psi_eta <- Psi %*% Sigma_diag
  zstar <- z[,2] + Psi_eta
  z$zstar <- zstar
  
  ################################################################################
  ## Scale z for Stan - initial
  ################################################################################
  z_scale <- scale(cbind(z[,2], z[, "zstar"]))
  
  ################################################################################
  ## Create Stan data list
  ################################################################################
  data_list <- list(
    N = nrow(y),
    y = model.frame(f, y)[,1],
    K = ncol(x),
    Q = 1,
    J = J,
    M = M,
    x = x,
    z = as.matrix(z_scale[,1]),
    zstar = as.matrix(z_scale[,2]),
    ur = ur,
    D_JxM = D_JxM,
    Dstar = Dstar,
    IMAT = IMAT,
    MJJ = M_JxJ,
    Sigma_diag = as.matrix(Sigma_diag),
    phi = phi
  )
  
  return(data_list)
}
