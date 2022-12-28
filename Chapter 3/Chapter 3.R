### CHAPTER 3 ###

###  Includes (amongst others) the following objects:                   ###
###  The first order Maximum Mean Discrepancy test statistic for the    ###
###  signature kernel, the second order Maximum Mean Discrepancy,       ###
###  the experiment in Section 3.4, the Maximum Mean Discrepancy for    ###
###  different kernels                                                  ###

### install packages ###

install.packages("pracma")
install.packages("moments")
install.packages("einsum")
install.packages("abind")
install.packages("xtable")
library(pracma)
library(moments)
library(einsum)
library(abind)
library(xtable)

### The functions sig and BM from Chapter 2 must be initialized before ###
### executing this code.                                               ###

### MMD with the signature kernel for the normalized signature ###
### (followed the approach in Github Marketsimulator)          ###

SigMMD <- function(set1,set2,truncorder){
  
  m <- dim(set1)[3]
  n <- dim(set2)[3]
  
  #X
  d <- dim(set1)[1]
  X <- c()  #rep(0,m*(((d*((d^truncorder)-1))/(d-1))+1))
  for (i in (1:m)) {
    #path <- t(set1[,,i])
    
    path <- set1[,,i]
    ifelse(is.null(dim(path)),dimpath<-1,dimpath<-dim(path)[1])
    ifelse(is.null(dim(path)),lenpath<-length(path),lenpath<-dim(path)[2])
    
    #function Phi
    #signa <- c(1,sig(path,truncorder))
    
    signa <- c(1,sig(array(path,c(dimpath,lenpath)),truncorder))
    
    #function get_keys
    tuples <- as.vector(0)
    for (j in 1:truncorder) {
      tuples<-c(tuples,rep(j,d^(j)))
    }
    #keys <- tuples
    
    #function phi
    a<- signa^2
    
    #function psi
    psi <- norm(signa,"2")^2
    if (psi > 4){psi <- 4+((4^(1+1))*((4^(-1))-(psi^(-1))))/1}
    
    a[1] <- a[1] - psi
    f<-function(x) {     
      z<-x^(2*(0:(length(a)-1)))
      return(as.numeric(a%*%z))
    }
    
    phi_x <- brent(f,0,5)$root 
    
    lambda <- phi_x^(tuples)
    Phi<-lambda*signa
    #X[(1+(i-1)*(((d*((d^truncorder)-1))/(d-1))+1)):(i*(((d*((d^truncorder)-1))/(d-1))+1))] <- Phi
    X<-c(X,Phi)
  }
  
  #Y
  d <- dim(set2)[1]
  Y <- c()   #rep(0,n*(((d*((d^truncorder)-1))/(d-1))+1))
  for (i in (1:n)) {
    #path <- t(set2[,,i])
    
    path <- set2[,,i]
    ifelse(is.null(dim(path)),dimpath<-1,dimpath<-dim(path)[1])
    ifelse(is.null(dim(path)),lenpath<-length(path),lenpath<-dim(path)[2])
    
    #function Phi
    #signa <- c(1,sig(path,truncorder))
    
    signa <- c(1,sig(array(path,c(dimpath,lenpath)),truncorder))
    
    #function get_keys
    tuples <- as.vector(0)
    for (j in 1:truncorder) {
      tuples<-c(tuples,rep(j,d^(j)))
    }
    #keys <- tuples
    
    #function phi
    a<- signa^2
    
    #function psi
    psi <- norm(signa,"2")^2
    if (psi > 4){psi <- 4+((4^(1+1))*((4^(-1))-(psi^(-1))))/1}
    
    a[1] <- a[1] - psi
    f<-function(x) {     
      z<-x^(2*(0:(length(a)-1)))
      return(as.numeric(a%*%z))
    }
    
    phi_x <- brent(f,0,5)$root 
    
    lambda <- phi_x^(tuples)
    Phi<-lambda*signa
    #Y[(1+(i-1)*(((d*((d^truncorder)-1))/(d-1))+1)):(i*(((d*((d^truncorder)-1))/(d-1))+1))] <- Phi
    Y<-c(Y,Phi)
  }
  
  
  TU <- (1/(m^2))*sum(X%*%t(X)) - (2/(m*n))*sum(X%*%t(Y)) + (1/(n^2))*sum(Y%*%t(Y))
  return(TU)
}

### Associated test statistic (paragraph 3.2.5. in "A data-driven market  ###
### simulator for small data environments")                               ###
### with adjusted threshold from "A Kernel Two-Sample Test" (Corollary 9) ###

SigMMD_test <- function(set1,set2,truncorder,conflevel){
  m <- dim(set1)[3]
  #c_alpha<-4*sqrt(-log(conflevel)/m)
  K <- 4 + 4/1
  c_alpha <- sqrt((2*K)/m)*(1+sqrt(2*log(1/conflevel)))
  return(paste("The samples come from the same distribution:", SigMMD(set1,set2,truncorder)<=c_alpha))
}

### Second order MMD                                  ###
### followed the approach from Github higherOrderKME  ### 
### commented out is the original Python code         ###

algo1 <- function(G_static){
  #cdef int A = G_static.shape[0]
  #cdef int B = G_static.shape[1]
  #cdef int M = G_static.shape[2]
  #cdef int N = G_static.shape[3]
  #cdef int i, j, l, m
  A <- dim(G_static)[3]
  B <- dim(G_static)[4]
  M <- dim(G_static)[1]
  N <- dim(G_static)[2]
  
  #cdef double[:,:,:,:] K = np.zeros((A,B,M+1,N+1), dtype=np.float64)
  K <- array(0,c(M+1,N+1,A,B))
  
  #for l in range(A):
  #  for m in range(B):
  #     for i in range(M+1):
  #       K[l,m,i,0] = 1.
  #     for j in range(N+1):
  #       K[l,m,0,j] = 1.
  #     for i in range(M):
  #       for j in range(N):
  #         K[l,m,i+1,j+1] = (K[l,m,i+1,j] + K[l,m,i,j+1])*(1. + 0.5*G_static[l,m,i,j]+(1./12)*G_static[l,m,i,j]**2) - K[l,m,i,j]*(1. - (1./12)*G_static[l,m,i,j]**2)
  
  for (l in 1:A) {
    for (m in 1:B) {
      for (i in (1:(M+1))) {
        K[i,1,l,m] = 1
      }
      for (j in (1:(N+1))) {
        K[1,j,l,m] = 1
      }
      for (i in 1:M) {
        for (j in 1:N) {
          K[i+1,j+1,l,m] = (K[i+1,j,l,m] + K[i,j+1,l,m])*(1 + 0.5*G_static[i,j,l,m] + (1/12)*G_static[i,j,l,m]^2) - K[i,j,l,m]*(1 - (1/12)*G_static[i,j,l,m]^2)
        }
      }
    }
  }
  
  return(K)
}

algo2 <- function(X,Y,sigma,return_sol_grid=FALSE){
  
  #A, B, M, N = X.shape[0], Y.shape[0], X.shape[1], Y.shape[1]
  #X.shape[0] = number of paths, X.shape[1] = length of path, X.shape[2] = dimension of path
  A <- dim(X)[3]
  B <- dim(Y)[3]
  M <- dim(X)[2]
  N <- dim(Y)[2]
  
  
  #Xs = torch.sum(X**2, dim=2)
  #Ys = torch.sum(Y**2, dim=2)
  Xs <- colSums(X^2,dims = 1)
  Ys <- colSums(Y^2,dims = 1)
  #dist = -2.*torch.einsum('ipk,jqk->ijpq', X, Y)
  #dist += torch.reshape(Xs, (A, 1, M, 1)) + torch.reshape(Ys, (1, B, 1, N))
  
  dist <- array(0,c(M,N,A,B)) 
  for (i in 1:M) {
    for (j in 1:N) {
      for (k in 1:A) {
        for (l in 1:B) {
          dist[i,j,k,l] <- Xs[i,k] + Ys[j,l]
        }
      }
    }
  }
  dist <- dist - 2*einsum("kpi,kqj -> pqij",X,Y)
  
  #G_static = torch.exp(-dist/self.sigma)
  G_static <- exp(-dist/sigma)
  
  
  #G_static_ = G_static[:, :, 1:, 1:] + G_static[:, :, :-1, :-1] - G_static[:, :, 1:, :-1] - G_static[:, :, :-1, 1:] 
  #G_static_ = tile(tile(G_static_, 2, 2**dyadic_order)/float(2**dyadic_order), 3, 2**dyadic_order)/float(2**dyadic_order)
  G_static_ <- G_static[2:M,2:N,,] + G_static[(1:(M-1)),(1:(N-1)),,] - G_static[(1:(M-1)),2:N,,] - G_static[2:M,(1:(N-1)),,]
  
  init_dim <- dim(G_static_)[1] 
  repeat_idx <- rep(1,length(dim(G_static_)))
  repeat_idx[3] <- 2^2
  a <- abind(G_static_,G_static_,G_static_,G_static_,along = 1)
  order_index <- c()
  for (i in (0:(init_dim-1))) { 
    order_index <- c(order_index,init_dim*(0:3) + i) 
  }
  tile1 <- array(0,c(1,dim(a)[2],dim(a)[3],dim(a)[4]))
  for (i in order_index) {
    tile1 <- abind(tile1,a[i+1,,,],along = 1)
  }
  tile1 <- tile1[2:(dim(tile1)[1]),,,]
  
  tile1 <- tile1/(2^2)
  init_dim <- dim(tile1)[2]
  repeat_idx <- rep(1,length(dim(tile1)))
  repeat_idx[4] <- 2^2
  a <- abind(tile1,tile1,tile1,tile1,along = 2)
  order_index <- c()
  for (i in (0:(init_dim-1))) { 
    order_index <- c(order_index,init_dim*(0:3) + i) 
  }
  tile <- array(0,c(dim(a)[1],1,dim(a)[3],dim(a)[4]))
  for (i in order_index) {
    tile <- abind(tile,a[,i+1,,],along = 2) 
  }
  tile <- tile[,2:(dim(tile)[2]),,]
  G_static_ <- tile/(2^2)
  
  
  #G = torch.tensor(sig_kernel_Gram_varpar(G_static_.detach().numpy(), sym, _naive_solver), dtype=G_static.dtype, device=G_static.device)
  #sig_kernel_Gram_varpar = algo1
  
  G <- algo1(G_static = G_static_)
  
  #if not return_sol_grid:
  #  return G[:, :, -1, -1]
  #else:
  #  return G[:, :, ::2**dyadic_order, ::2**dyadic_order]
  
  if (return_sol_grid == TRUE){return(G[seq(1,dim(G)[1],2^2),seq(1,dim(G)[2],2^2),,])}
  else {return(G[dim(G)[1],dim(G)[2],,])}
  
}

algo4 <- function(K_XX,K_XY,K_YY,sigma,lambda){
  
  #G_base = innerprodCKME(K_XX, K_XY, K_YY, lambda_, static_kernel, sym=sym)    # <---- this is the only change compared to order 1
  #G_base_ = G_base[:, :, 1:, 1:] + G_base[:, :, :-1, :-1] - G_base[:, :, 1:, :-1] - G_base[:, :, :-1, 1:] 
  #G_base_ = tile(tile(G_base_, 2, 2**dyadic_order)/float(2**dyadic_order), 3, 2**dyadic_order)/float(2**dyadic_order)
  #innerprodCKME = algo6
  
  G_base <- algo6(K_XX,K_XY,K_YY,sigma = sigma,lambda = lambda)
  G_base_ <- G_base[(2:(dim(G_base)[1])),(2:(dim(G_base)[2])),,] + G_base[(1:(dim(G_base)[1]-1)),(1:(dim(G_base)[2]-1)),,] - G_base[(1:(dim(G_base)[1]-1)),(2:(dim(G_base)[2])),,] - G_base[(2:(dim(G_base)[1])),(1:(dim(G_base)[2]-1)),,]
  
  init_dim <- dim(G_base_)[1] 
  repeat_idx <- rep(1,length(dim(G_base_)))
  repeat_idx[3] <- 2^2
  a <- abind(G_base_,G_base_,G_base_,G_base_,along = 1)
  order_index <- c()
  for (i in (0:(init_dim-1))) { 
    order_index <- c(order_index,init_dim*(0:3) + i) 
  }
  tile1 <- array(0,c(1,dim(a)[2],dim(a)[3],dim(a)[4]))
  for (i in order_index) {
    tile1 <- abind(tile1,a[i+1,,,],along = 1)
  }
  tile1 <- tile1[2:(dim(tile1)[1]),,,]
  
  tile1 <- tile1/(2^2)
  init_dim <- dim(tile1)[2]
  repeat_idx <- rep(1,length(dim(tile1)))
  repeat_idx[4] <- 2^2
  a <- abind(tile1,tile1,tile1,tile1,along = 2)
  order_index <- c()
  for (i in (0:(init_dim-1))) { 
    order_index <- c(order_index,init_dim*(0:3) + i) 
  }
  tile <- array(0,c(dim(a)[1],1,dim(a)[3],dim(a)[4]))
  for (i in order_index) {
    tile <- abind(tile,a[,i+1,,],along = 2) 
  }
  tile <- tile[,2:(dim(tile)[2]),,]
  G_base_ <- tile/(2^2)
  
  #G = torch.tensor(sig_kernel_Gram_varpar(G_base_.detach().numpy(), sym, _naive_solver), dtype=G_base.dtype, device=G_base.device)
  #sig_kernel_Gram_varpar = algo1
  
  G <- algo1(G_static = G_base_)
  
  #return G[:, :, -1, -1]
  return(G[dim(G)[1],dim(G)[2],,])
}

algo5 <- function(X,Y,sigma,lambda){
  
  #K_XX_1 = self.compute_Gram(X, X, sym=True, return_sol_grid=True)   # shape (batch_X, batch_X, length_X, length_X)
  #K_YY_1 = self.compute_Gram(Y, Y, sym=True, return_sol_grid=True)   # shape (batch_Y, batch_Y, length_Y, length_Y)
  #K_XY_1 = self.compute_Gram(X, Y, sym=False, return_sol_grid=True)  # shape (batch_X, batch_Y, length_X, length_Y)
  #compute_Gram = algo2(X,Y,sigma,True)
  
  K_XX_1 <- algo2(X,X,sigma = sigma,return_sol_grid = TRUE)
  K_YY_1 <- algo2(Y,Y,sigma = sigma,return_sol_grid = TRUE)
  K_XY_1 <- algo2(X,Y,sigma = sigma,return_sol_grid = TRUE)
  
  #K_XX_2 = self.compute_HigherOrder_Gram(K_XX_1, K_XX_1, K_XX_1, lambda_, sym=True)       # shape (batch_X, batch_X)
  #K_YY_2 = self.compute_HigherOrder_Gram(K_YY_1, K_YY_1, K_YY_1, lambda_, sym=True)       # shape (batch_Y, batch_Y)
  #K_XY_2 = self.compute_HigherOrder_Gram(K_XX_1, K_XY_1, K_YY_1, lambda_, sym=False)      # shape (batch_X, batch_Y)
  #compute_HigherOrder_Gram = algo4
  
  K_XX_2 <- algo4(K_XX_1,K_XX_1,K_XX_1,sigma = sigma,lambda = lambda)
  K_YY_2 <- algo4(K_YY_1,K_YY_1,K_YY_1,sigma = sigma,lambda = lambda)
  K_XY_2 <- algo4(K_XX_1,K_XY_1,K_YY_1,sigma = sigma,lambda = lambda)
  
  #return torch.mean(K_XX_2) + torch.mean(K_YY_2) - 2.*torch.mean(K_XY_2)
  return(mean(K_XX_2) + mean(K_YY_2) - 2*mean(K_XY_2))
}

algo6 <- function(K_XX,K_XY,K_YY,sigma,lambda){
  
  #m  = K_XX.shape[0]
  #n  = K_YY.shape[0]
  #LX = K_XX.shape[2]
  #LY = K_YY.shape[2]
  
  m <- dim(K_XX)[3]
  n <- dim(K_YY)[3]
  LX <- dim(K_XX)[1]
  LY <- dim(K_YY)[1]
  
  #H_X = torch.eye(m, device=K_XX.device, dtype=K_XX.dtype) - (1./m)*torch.ones((m, m), device=K_XX.device, dtype=K_XX.dtype)       # centering matrix
  #H_Y = torch.eye(n, device=K_YY.device, dtype=K_YY.dtype) - (1./n)*torch.ones((n, n), device=K_YY.device, dtype=K_YY.dtype)       # centering matrix
  
  H_X <- diag(1,m,m) - matrix(1/m,m,m)
  H_Y <- diag(1,n,n) - matrix(1/n,n,n)
  
  #G = torch.zeros((m, n, LX, LY), dtype=K_XX.dtype, device=K_XX.device)
  
  G <- array(0,c(LX,LY,m,n))
  
  #to_inv_XX = torch.zeros((m, m, LX), dtype=K_XX.dtype, device=K_XX.device)
  #for p in range(LX):
  #  to_inv_XX[:, :, p] = K_XX[:, :, p, p]
  #to_inv_XX += lambda_*torch.eye(m, dtype=K_XX.dtype, device=K_XX.device)[:, :, None]
  #to_inv_XX = to_inv_XX.T
  #inv_X = torch.linalg.inv(to_inv_XX)
  #inv_X = inv_X.T
  
  to_inv_XX <- array(0,c(m,LX,m))
  dia <- array(0,c(m,1,m))
  for (i in 1:LX) {
    to_inv_XX[,i,] <- K_XX[i,i,,]
    dia <- abind(dia,array(diag(lambda,m,m),c(m,1,m)),along = 2)
  }
  dia <- dia[,2:dim(dia)[2],]
  to_inv_XX <- to_inv_XX + dia
  to_inv_XX <- aperm(to_inv_XX,perm = c(1,3,2))
  inv_X <- array(0,c(dim(to_inv_XX)[1],dim(to_inv_XX)[2],dim(to_inv_XX)[3]))
  for (i in (1:dim(to_inv_XX)[3])) {
    inv_X[,,i] <- inv(to_inv_XX[,,i]) 
  }
  inv_X <- aperm(inv_X,perm = c(1,3,2))
  
  #to_inv_YY = torch.zeros((n, n, LY), dtype=K_YY.dtype, device=K_YY.device)
  #for q in range(LY):
  #  to_inv_YY[:, :, q] = K_YY[:, :, q, q]
  #to_inv_YY += lambda_*torch.eye(n, dtype=K_YY.dtype, device=K_YY.device)[:, :, None]
  #to_inv_YY = to_inv_YY.T
  #inv_Y = torch.linalg.inv(to_inv_YY)
  #inv_Y = inv_Y.T
  
  to_inv_YY <- array(0,c(n,LY,n))
  dia <- array(0,c(n,1,n))
  for (i in 1:LY) {
    to_inv_YY[,i,] <- K_YY[i,i,,]
    dia <- abind(dia,array(diag(lambda,n,n),c(n,1,n)),along = 2)
  }
  dia <- dia[,2:dim(dia)[2],]
  to_inv_YY <- to_inv_YY + dia
  to_inv_YY <- aperm(to_inv_YY,perm = c(1,3,2))
  inv_Y <- array(0,c(dim(to_inv_YY)[1],dim(to_inv_YY)[2],dim(to_inv_YY)[3]))
  for (i in (1:dim(to_inv_YY)[3])) {
    inv_Y[,,i] <- inv(to_inv_YY[,,i]) 
  }
  inv_Y <- aperm(inv_Y,perm = c(1,3,2))
  
  #for p in range(LX): # TODO: to optimize (e.g. when X=Y)
  #   WX = inv_X[:, :, p]
  #   WX_ = torch.matmul(K_XX[:, :, p, p].t(), WX)
  #   for q in range(LY):
  #     WY = inv_Y[:, :, q]  
  #     WY_ = torch.matmul(WY, K_YY[:, :, q, q])
  #     G_cross  = -2*torch.matmul(WX_, torch.matmul(K_XY[:, :, -1, -1], WY_))
  #     G_row = torch.diag(torch.matmul(WX_, torch.matmul(K_XX[:, :, -1, -1], WX_.t())))[:, None]
  #     G_col = torch.diag(torch.matmul(WY_.t(), torch.matmul(K_YY[:, :, -1, -1],WY_)))[None, :]
  #     dist  = G_cross + G_row + G_col
  #     G[:,:,p,q] = torch.exp(-dist/static_kernel.sigma)
  
  for (p in 1:LX) {
    WX <- inv_X[,p,]
    WX_ <- t(K_XX[p,p,,])%*%WX
    for (q in 1:LY) {
      WY <- inv_Y[,q,]
      WY_ <- WY%*%K_YY[q,q,,]
      G_cross <- -2*(WX_%*%(K_XY[dim(K_XY)[1],dim(K_XY)[2],,]%*%WY_))
      G_row <- diag(WX_%*%(K_XX[dim(K_XX)[1],dim(K_XX)[2],,]%*%t(WX_)))
      G_col <- diag(t(WY_)%*%(K_YY[dim(K_YY)[1],dim(K_YY)[2],,]%*%WY_))
      dist <- matrix(0,dim(G_cross)[1],dim(G_cross)[2])
      for (i in 1:dim(dist)[1]) {
        for (j in 1:dim(dist)[2]) {
          dist[i,j] <- G_cross[i,j] + G_row[i] + G_col[j]
        }
      }
      G[p,q,,] <- exp(-dist/sigma)
    }
  }
  
  return(G)
}

secondorderMMD <- function(X,Y,sigma = 0.5,lambda = 0.00001){
  return(algo5(X,Y,sigma,lambda))
}

firstordersigMMD <- function(X,Y,sigma=0.5){
  
  #K_XX = self.compute_Gram(X, X, sym=True)
  #K_YY = self.compute_Gram(Y, Y, sym=True)
  #K_XY = self.compute_Gram(X, Y, sym=False)
  #compute_Gram = algo2(X,Y,sigma)
  
  K_XX <- algo2(X,X,sigma = sigma,return_sol_grid = FALSE)
  K_YY <- algo2(Y,Y,sigma = sigma,return_sol_grid = FALSE)
  K_XY <- algo2(X,Y,sigma = sigma,return_sol_grid = FALSE)
  
  
  #return torch.mean(K_XX) + torch.mean(K_YY) - 2.*torch.mean(K_XY)
  return(mean(K_XX) + mean(K_YY) - 2*mean(K_XY))
}


### The experiment from Section 3.4 ###

GBM2vola <- function(n,T,S0,drift,diffusion1,diffusion2){
  time_grid <- seq(0,T,T/(n-1))
  #time_increment <- time_grid[2]
  Brownian <- BM(T,1,n-1)
  #Brownianincrement <- diff(Brownian)
  gbm <- matrix(S0,2,n)
  #for (t in (1:(n-1))) {
  #gbm[t+1] = gbm[t] + gbm[t]*(drift*time_increment + diffusion*Brownianincrement[t])
  #}
  for (t in (2:n)) {
    gbm[1,t] = S0*exp((drift-0.5*(diffusion1)^2)*time_grid[t] + diffusion1*Brownian[t])
    gbm[2,t] = S0*exp((drift-0.5*(diffusion2)^2)*time_grid[t] + diffusion2*Brownian[t])
  }
  return(gbm)
}

### Table in the thesis ### 

r <- 0
S0 <- 34
sigma1 <- 0.2
length_of_paths <- 15
paths_to_be_simulated <- 50
conflevel <- 0.95
c_alpha <- sqrt((2*8)/paths_to_be_simulated)*(1+sqrt(2*log(1/conflevel)))

TableSec34 <- matrix(0,12,5)
TableSec34[1,3:5] <- "?" #"Do the samples come from the same distributions?"
TableSec34[2,1:5] <- c("n","d_BS(sigma_1,sigma_2)","K-S","MMD_1","MMD_2")
dbs <- c(1,10,50,100,200,400,600,1000,2000,4000) 
T <- c(5/252,63/252,252/252)
for (i in (1:10)) {
  n <- dbs[i]
  sigma2 <- sigma1 + 1/n
  TableSec34[i+2,1] <- n
  TableSec34[i+2,2] <- 1/n
  for (j in (3:3)) {
    t <- T[j]
    stock_prices1 <- array(0,c(2,length_of_paths,paths_to_be_simulated))
    stock_prices2 <- array(0,c(2,length_of_paths,paths_to_be_simulated))
    for (l in 1:paths_to_be_simulated) {
      stocks <- GBM2vola(length_of_paths,t,S0,r,sigma1,sigma2)
      stock_prices1[1,,l] <- seq(0,t,t/(length_of_paths-1))
      stock_prices1[2,,l] <- stocks[1,]
      stock_prices2[1,,l] <- seq(0,t,t/(length_of_paths-1))
      stock_prices2[2,,l] <- stocks[2,]
    }
    
    TableSec34[2+i,(3*1)] <- ifelse((ks.test(stock_prices1[2,,],stock_prices2[2,,])$p.value >= conflevel),"Yes","No") 
    TableSec34[2+i,(3*1)+1] <- ifelse((SigMMD(stock_prices1,stock_prices2,5) >= c_alpha),"No","Yes")
    TableSec34[2+i,(3*1)+2] <- ifelse((secondorderMMD(stock_prices1,stock_prices2,0.5) >= c_alpha),"No","Yes")
  }
}
TableSec34

print(xtable(TableSec34,type="latex"),file = "sec34.tex")


### MMD for the linear, rbf and laplace kernel ###

#First we have to define the different kernels

rbf <- function(x,y,sigma=1,norm="2"){
  return(exp((-norm(x-y,type=norm)^2)/(2*sigma^2)))
}

laplace <- function(x,y,sigma=1,norm="2"){
  return(exp((-norm(x-y,type=norm))/sigma))
}

matern32 <- function(x,y,sigma=1,norm="2"){
  return((1 + (sqrt(3)*norm(x-y,type=norm))/(sigma^2))*exp(-(sqrt(3)*norm(x-y,type=norm))/(sigma^2)))
}

#for now don't change the type of norm... because it only works for matrices and 
#vectors for "2"
#Norm() allows other norms for vectors but not for matrices...

linear <- function(x,y){
  return(t(x)%*%y)
}


MMD <- function(X,Y,k,sigma=1,signa=FALSE,truncorder=4){
  A <- X
  B <- Y
  numpathA <- dim(A)[3]
  numpathB <- dim(B)[3]
  
  if (signa == TRUE){
    dimpath<-dim(X)[1]
    lenpath<-dim(X)[2]
    A <- array(0,c(1,length(c(1,sig(array(X[,,1],c(dimpath,lenpath)),truncorder))),numpathA))
    B <- array(0,c(1,length(c(1,sig(array(X[,,1],c(dimpath,lenpath)),truncorder))),numpathA))
    for (i in (1:numpathA)) {
      A[,,i] <- c(1,sig(array(X[,,i],c(dimpath,lenpath)),truncorder))
      B[,,i] <- c(1,sig(array(Y[,,i],c(dimpath,lenpath)),truncorder))
    }
  }
  
  if (k=="rbf"){
    TU <- 0
    for (i in (1:numpathA)) {
      for (j in (1:numpathA)) {
        TU <- TU + (rbf(A[,,i],A[,,j],sigma) - 2*rbf(A[,,i],B[,,j],sigma) + rbf(B[,,i],B[,,j],sigma))
      }
    }
  }
  else if (k=="laplace"){
    TU <- 0
    for (i in (1:numpathA)) {
      for (j in (1:numpathA)) {
        TU <- TU + (laplace(A[,,i],A[,,j],sigma) - 2*laplace(A[,,i],B[,,j],sigma) + laplace(B[,,i],B[,,j],sigma))
      }
    }
  }
  else if (k=="linear"){
    TU <- 0
    for (i in (1:numpathA)) {
      for (j in (1:numpathA)) {
        TU <- TU + (linear(A[,,i],A[,,j]) - 2*linear(A[,,i],B[,,j]) + linear(B[,,i],B[,,j]))
      }
    }
  }
  else if (k=="matern"){
    TU <- 0
    for (i in (1:numpathA)) {
      for (j in (1:numpathA)) {
        TU <- TU + (matern32(A[,,i],A[,,j],sigma) - 2*matern32(A[,,i],B[,,j],sigma) + matern32(B[,,i],B[,,j],sigma))
      }
    }
  }
  else {stop("Kernel not defined.")}
  #other norm than 2 is atm not implemented
  return((1/(numpathA^2))*TU)
}

### Other metrics that we did not continue to elaborate as the thesis ###
### was orientated differently                                        ###

### Sig-W1 presented in "Conditional Sig-Wasserstein GANs ### 
### for Time Series Generation"                           ###
### Function exp_sig from Chapter 2 needs to be           ### 
### initialized first                                     ###

SigW1 <- function(X,Y,truncorder){
  dimpath <- dim(X)[1]
  lenpath <- dim(X)[2]
  numpath <- dim(X)[3]
  sigX <- array(0,c(1,length(c(1,sig(array(X[,,1],c(dimpath,lenpath)),truncorder))),numpath))
  sigY <- array(0,c(1,length(c(1,sig(array(X[,,1],c(dimpath,lenpath)),truncorder))),numpath))
  for (i in (1:numpath)) {
    sigX[,,i] <- c(1,sig(array(X[,,i],c(dimpath,lenpath)),truncorder))
    sigY[,,i] <- c(1,sig(array(Y[,,i],c(dimpath,lenpath)),truncorder))
  }
  #maybe use theoretical value for length of sig -> gives problem when d=1
  return(norm(exp_sig(sigX) - exp_sig(sigY),type = "2"))
}

### Kullback-Leibler Divergence ###

KL <- function(Finv,f,g,n){
  X <- Finv(runif(n,0,1))
  return((1/n)*sum(log(f(X)/g(X))))
}
# at the moment this function still has a problem if used with the uniform distribution
# because f and g are then constant and f(X) and g(X) give a constant then and not a 
# vector of constants (atm one would need to multiply the result by n...)

# Example for the exponential law with lambda=1 for X~F and lambda=2 for X~G

Finv <- function(x){
  return(-log(1-x))
}
f <- function(x){
  return(exp(-x))
}
g <- function(x){
  return(2*exp(-2*x))
}

Ginv <- function(x){
  return((-1/2)*log(1-x))
}

KL(Finv,f,g,1000000) #theoretical value D(F||G)=1-ln(2)
KL(Ginv,g,f,1000000) #theoretical value D(G||F)=ln(2)-0.5
1-log(2)
log(2)-0.5

### Metric on marginal distribution ###

marginaldist <- function(X,Y){
  X <- X[,,1] 
  Y <- Y[,,1]
  #for now only with one path. 
  
  ifelse(is.null(dim(X)),d<-1,d<-dim(X)[1])
  sum <- 0
  for (i in (1:d)) {
    firstpoint <- min(min(X[i,]),min(Y[i,]))
    lastpoint <- max(max(X[i,]),max(Y[i,]))
    densX <- hist(X[i,],breaks = seq(firstpoint,lastpoint,(lastpoint - firstpoint)/100),plot = FALSE)$density
    densY <- hist(Y[i,],breaks = seq(firstpoint,lastpoint,(lastpoint - firstpoint)/100),plot = FALSE)$density
    sum <- sum + sum(abs(densX-densY))
  }
  
  return((1/d)*sum)
}

### Moment matching ###

moment_test <- function(X,Y,threshold){
  if((dim(X)[1] > 1) || (dim(Y)[1] > 1)){stop("d > 1.")} 
  X <- X[,,1] 
  Y <- Y[,,1] 
  #for now only for the first path
  ord <- 1
  
  while((abs(moment(X,order = ord) - moment(Y,order = ord)) <= threshold) && (ord <=20)){
    ord <- ord + 1
  }
  return(paste("The moments match up to order",ord-1))
}


### Absolute difference of lag-1 auto-correlation ###

lag1autocorr <- function(X,Y){
  X <- X[,,1] 
  Y <- Y[,,1]
  #for now only with one path. 
  
  ifelse(is.null(dim(X)),d<-1,d<-dim(X)[1])
  ifelse(is.null(dim(X)),lenpathX<-length(X),lenpathX<-dim(X)[2])
  ifelse(is.null(dim(Y)),lenpathY<-length(Y),lenpathY<-dim(Y)[2])
  
  sum <- 0
  for (i in 1:d) {
    Xi <- X[i,]
    rhoreal0 <- (1/lenpathX)*sum((Xi-mean(Xi))^2)
    rhoreal1 <- (1/(lenpathX-1))*sum((Xi[1:(lenpathX-1)]-mean(Xi))*(Xi[2:lenpathX]-mean(Xi)))
    
    Yi <- Y[i,]
    rhogen0 <- (1/lenpathY)*sum(Yi[2:(lenpathY+1)]*Yi[1:lenpathY]) - (1/lenpathY)*sum(Yi[2:(lenpathY+1)])*(1/lenpathY)*sum(Yi[1:lenpathY])
    rhogen1 <- (1/lenpathY)*sum(Yi[2:(lenpathY+1)]^2) - ((1/lenpathY)^2)*(sum(Yi[2:(lenpathY+1)])^2)
    
    sum <- sum + abs((rhoreal1/rhoreal0)-(rhogen1/rhogen0))
  }
  
  return((1/d)*sum)
}#obvious problem with lenpathY+1... need to fully understand the theory first


### Metric on the correlation ###
### did not continue this     ###

corrmetric <- function(X,Y){
  X <- X[,,1] 
  Y <- Y[,,1]
  #for now only with one path.  
  
  ifelse(is.null(dim(X)),d<-1,d<-dim(X)[1])
  ifelse(is.null(dim(X)),lenpathX<-length(X),lenpathX<-dim(X)[2])
  ifelse(is.null(dim(Y)),lenpathY<-length(Y),lenpathY<-dim(Y)[2])
  
  covreal <- matrix(0,d,d)
  for (i in 1:d) {
    for (j in 1:d) {
      covreal[i,j] <- (1/lenpathX)*sum(X[i,]*X[j,]) - (1/lenpathX)*sum(X[i,])*(1/lenpathX)*sum(X[j,])
    }
  }
  
  
  
  corrreal <- matrix(0,d,d)
  corrgen <- matrix(0,d,d)
  for (i in 1:d) {
    for (j in 1:d) {
      corrreal <- covreal[i,j]/sqrt(covreal[i,i]*covreal[j,j])
      corrgen <- covgen[i,j]/sqrt(covgen[i,i]*covgen[j,j])
    }
  }
  
  return(norm(corrreal - corrgen,type = "1"))
}
