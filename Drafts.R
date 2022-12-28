### NOT FOR PRESENTATION ###
### NOT RELEVANT         ###


### Signature of a path ###


#sig(array(c(1,1,3,10,20,30),c(3,2)),4)

## Generate BM
#BM = function(T,n){
#  delta = T/n;
#  W = c(0, cumsum(sqrt(delta)*rnorm(n,0,1)));  
#  return(W)
#}
#BM = function(T,n){
#  W <- rep(0,n)
#  delta = T/n
#  for (t in 2:n) {
#    W[t] <- W[t-1] + rnorm(1,0,delta)
#  }
#  return(W)
#}

#p <- 1 #number of paths
#T <- 1
#n <- 1000
#W= matrix(0,n+1,p)
#for (i in 1:p){
#  W[,i] = BM(T,n) #each column of the matrix contains one path of a Brownian motion
#}

# signature

m <- 200
t <- 4
#for (t in 2:N) {
expsig <- rep(0,length(sig(array(BM(T,d,n),c(d,n+1)),t)))
for (j in 1:length(sig(array(BM(T,d,n),c(d,n+1)),t))) {
  for (i in 1:m) {
    expsig[j] <- expsig[j] + sig(array(BM(T,d,n),c(d,n+1)),t)[j]
  }
  expsig[j] <- expsig[j]/m
}
expsig <- c(1,expsig)
expsig
matplot(seq(0,t,1),abs(expsig),"l")
#normexpsig <- 1 + sum(abs(expsig))

#}
#factorial decay  
res <- rep(0,t+1)
for (i in 2:(t+1)) {
  res[i] <- (1/(i-1))*sum(abs(expsig[1:i])) 
}
res[1] <- 1

matplot(0:t,res,"l")


#stopping <- function(l,truncorder){
#  return(paste(exp_shuffle(paste("-",paste(shuffle(l,l),1,sep=""),sep=""),truncorder),2,sep=""))
#}

simplify <- function(x){
  y <- x
  count <- 0
  for (i in (1:length(x))) {
    for (j in ((i+1):length(x))) {
      if (!(is.null(summation(x[i],x[j])))){
        y[i] <- summation(x[i],x[j]) 
        #x <- x[-j]
        count <- count+1
      }
    }
  }
  return(y)#[1:(length(y)-count)])
}

simplify <- function(x){#work to do here
  y <- x
  count <- 0
  for (i in (1:(length(x)-1))) {
    for (j in ((i+1):length(x))) {
      if (!(is.null(summation(y[i],x[j])))){
        y[i] <- summation(y[i],x[j]) 
        x <- x[-j]
        y <- y[-j]
        count <- count+1
      }
    }
  }
  return(y)#[1:(length(y)-count)])
}

summation <- function(a,b){
  ifelse(((strsplit(a," ")[[1]][2]==strsplit(b," ")[[1]][2]) || ((is.na(strsplit(a," ")[[1]][2])) && (is.na(strsplit(b," ")[[1]][2])))),return(paste(as.numeric(strsplit(a," ")[[1]][1]) + as.numeric(strsplit(b," ")[[1]][1]),strsplit(a," ")[[1]][2])),stop("Not the same word."))
}#check condition, problem when a or b is "l " because then strsplit(a," ")[[1]][2]=NA
#and NA can not be compared(gives NA)
#maybe try with normal if conditions...

summation <- function(a,b){
  if ((is.na(strsplit(a," ")[[1]][2]))){
    if ((is.na(strsplit(b," ")[[1]][2]))){
      return(paste(as.numeric(strsplit(a," ")[[1]][1]) + as.numeric(strsplit(b," ")[[1]][1]),gsub("NA","",paste(strsplit(a," ")[[1]][2],"",sep=""))))
    }
    else {return("0")}#stop("Not the same word.")}
  }
  else {
    if ((is.na(strsplit(b," ")[[1]][2]))){return("0")
      #stop("Not the same word.")
    }
    else if (strsplit(a," ")[[1]][2]==strsplit(b," ")[[1]][2]){
      return(paste(as.numeric(strsplit(a," ")[[1]][1]) + as.numeric(strsplit(b," ")[[1]][1]),gsub("NA","",paste(strsplit(a," ")[[1]][2],"",sep=""))))
    }
    else {return("0")}#stop("Not the same word.")}
  }
}


exp_shuffle <-function(l,truncorder){
  sexp <- ""
  term <- ""
  fact <- c(1)
  if (l >= 0){
    for (r in (1:truncorder)) {
      term <- shuffle(term,l)
      sexp <- c(sexp,term)
      fact <- c(fact,rep(1/factorial(r),length(term)))
    }
  }
  else if (l < 0){
    for (r in (1:truncorder)) {
      term <- shuffle(term,substr(l,2,nchar(l)))
      sexp <- c(sexp,term)
      fact <- c(fact,rep(((-1)^r)/factorial(r),length(term)))
    }
  }
  return(paste(fact,sexp))
} #(any count duplicates function to include in the loop?)
#table function looks very interesting
#does this function work? maybe depends on if it is N or Z or R in Proposition 2 in Outline 



stopping <- function(l,truncorder){
  y <- c()
  arg_for_exp <- paste("-",paste(shuffle(l,l),1,sep=""),sep="")
  sexp <- matrix("",length(arg_for_exp),length(exp_shuffle(arg_for_exp[1],truncorder)))
  for (i in (1:length(arg_for_exp))) {
    sexp[i,] <- exp_shuffle(arg_for_exp[i],truncorder)
  }
  for (j in (1:(dim(sexp)[2]))) {
    for (k in (1:(dim(sexp)[2]))) {
      y <- c(y,concatenation(sexp[1,j],sexp[2,k]))
    }
  }
  return(paste(y,2,sep=""))
}#the function is ok but is only made for exp(aaa + aaa) i.e. for exp with 2 terms
#needs to be generalised to i terms



### Metrics ###


#KL1 <- function(Finv,f,g,n){
#  divergence <- 0
#  for (i in 1:n) {
#    X <- Finv(runif(1))
#    divergence <- divergence + log(f(X)/g(X)) 
#  }
#  return((1/n)*divergence)
#}


#test
set1<-array(0,c(1,1000+1,50))
set2<-array(0,c(1,1000+1,50))
for (i in 1:50) {
  set1[,,i]<-array(BM(1,1,1000),c(1,1000+1))
  set2[,,i]<-array(BM(1,1,1000),c(1,1000+1))
}

SigMMD(set1,set2,4)
#value of SigMMD varies a lot... + problem with higher dimensions


set1<-array(0,c(1,1000,50))
set2<-array(0,c(1,1000,50))
for (i in 1:50) {
  set1[,,i]<-array(rexp(1000),c(1,1000))
  set2[,,i]<-array(rnorm(1000),c(1,1000))
}
#This is a problem... Something does not work as it should...


#test
T<-1
d<-2
n<-1000
m<-100
c <- c()
b <- c()
for (i in (1:m)) {
  c <- c(c,BM(T,d,n))
  b <- c(b,BM(T,d,n))
}
X <- array(c,c(d,n+1,m))
Y <- array(b,c(d,n+1,m))
SigW1(X,X,4)


### Signature kernel to implement SigMMD without normalization ###

#sigkernel <- function(x,y,truncorder=4){
#  ifelse(is.null(dim(x)),dimpath<-1,dimpath<-dim(x)[1])
#  ifelse(is.null(dim(x)),lenpath<-length(x),lenpath<-dim(x)[2])
#  #dimpath <- dim(x)[1]
#  #lenpath <- dim(x)[2]
#  sigx <- c(1,sig(array(x,c(dimpath,lenpath)),truncorder))
#  sigy <- c(1,sig(array(y,c(dimpath,lenpath)),truncorder))
#  return(t(sigx)%*%sigy) #is this scalar product correct? as we are in T((R^d))... maybe it is here where the normalization constant comes into play?
#} #again it is assumed that x and y have the same dimensions


#MMD <- function(X,Y,k,sigma=1,truncorder=4){
#  m <- dim(X)[3]
#  n <- dim(Y)[3]
#  if (k=="rbf"){
#    TU <- 0
#    for (i in (1:m)) {
#      for (j in (1:m)) {
#        TU <- TU + (rbf(X[,,i],X[,,j],sigma) - 2*rbf(X[,,i],Y[,,j],sigma) + rbf(Y[,,i],Y[,,j],sigma))
#      }
#    }
#  }
#  else if (k=="laplace"){
#    TU <- 0
#    for (i in (1:m)) {
#      for (j in (1:m)) {
#        TU <- TU + (laplace(X[,,i],X[,,j],sigma) - 2*laplace(X[,,i],Y[,,j],sigma) + laplace(Y[,,i],Y[,,j],sigma))
#      }
#    }
#  }
#  else if (k=="linear"){
#    TU <- 0
#    for (i in (1:m)) {
#      for (j in (1:m)) {
#        TU <- TU + (linear(X[,,i],X[,,j]) - 2*linear(X[,,i],Y[,,j]) + linear(Y[,,i],Y[,,j]))
#      }
#    }
#  }
#  else if (k=="sigkernel"){
#    TU <- 0
#    for (i in (1:m)) {
#      for (j in (1:m)) {
#        TU <- TU + (sigkernel(X[,,i],X[,,j],truncorder) - 2*sigkernel(X[,,i],Y[,,j],truncorder) + sigkernel(Y[,,i],Y[,,j],truncorder))
#      }
#    }
#  }
#  else {stop("Kernel not defined.")}
#  #does it make sense to consider m and n different?
#  #other norm than 2 is atm not implemented
#  return((1/(m^2))*TU)
#}



SigMMD_test(set1,set2,4,0.95)


#einsum("kpi,kqj -> pqij",X,Y)
#array(Xs,c(M,1,A)) + array(Ys,c(N,1,B))
#array(ys,c(4,1,3,3))
#array(abind(array(ys,c(4,1,3)),array(ys,c(4,1,3)),array(ys,c(4,1,3)),array(ys,c(4,1,3)),along = 2),c(4,4,3,3)) + einsum("kpi,kqj -> pqij",x,y)

#XS <- c()
#for (i in 1:M) {
#  XS <- abind(XS,array(Xs,c(M,1,A)),along = 2)
#}
#XS <- array(XS,c(M,M,A,A))
#YS <- c()
#for (i in 1:N) {
#  YS <- abind(YS,array(Ys,c(1,N,B)),along = 2)
#}
#YS <- array(YS,c(N,N,B,B))
#dist <- -2*einsum("kpi,kqj -> pqij",X,Y) + XS + YS


G_base_ <- G_base[(2:(dim(G_base)[1])),(2:(dim(G_base)[2])),,] + G_base[(1:(dim(G_base)[1]-1)),(1:(dim(G_base)[2]-1)),,] - G_base[(1:(dim(G_base)[1]-1)),(2:(dim(G_base)[2])),,] - G_base[(2:(dim(G_base)[1])),(1:(dim(G_base)[2]-1)),,]

init_dim <- dim(G_base_)[2] #is it 2 in R too?
repeat_idx <- rep(1,length(dim(G_base_)))
repeat_idx[2] <- 2^2
a <-
  order_index <- c()
for (i in (0:(init_dim-1))) { #0:init_dim-1 or 1:init_dim
  order_index <- c(order_index,init_dim*(0:3) + i) #0:3 or 1:4? #i or i-1?
}
tile1 <- 
  
  tile1 <- tile1/(2^2)
init_dim <- dim(tile1)[3] #is it 3 in R too?
repeat_idx <- rep(1,length(dim(tile1)))
repeat_idx[3] <- 2^2
a <-
  order_index <- c()
for (i in (0:(init_dim-1))) { #0:init_dim-1 or 1:init_dim
  order_index <- c(order_index,init_dim*(0:3) + i) #0:3 or 1:4? #i or i-1?
}
tile <-
  G_base_ <- tile/(2^2)


ifelse(is.null(dim(X)),A <- 1,ifelse(is.na(dim(X)[3]),ifelse(dim(X)[1]>dim(X)[2],A <- dim(X)[2],A <- dim(X)[1]),A <- dim(X)[3])) #only works as long as the length of the path is greater than the number of dimensions or paths (temporary solution)
ifelse(is.null(dim(Y)),B <- 1,ifelse(is.na(dim(Y)[3]),ifelse(dim(Y)[1]>dim(Y)[2],B <- dim(Y)[2],B <- dim(Y)[1]),B <- dim(Y)[3]))
ifelse(is.null(dim(X)),M <- length(X),ifelse(is.na(dim(X)[3]),ifelse(dim(X)[1]>dim(X)[2],M <- dim(X)[1],M <- dim(X)[2]),M <- dim(X)[3]))
ifelse(is.null(dim(Y)),N <- length(Y),ifelse(is.na(dim(Y)[3]),ifelse(dim(Y)[1]>dim(Y)[2],N <- dim(Y)[1],N <- dim(Y)[2]),N <- dim(Y)[3]))



### LSMC approach ###



#CallBS <- function(St,K,T,t,r,sigma){
#    d1 <- (log(St/K) + (r+sigma^2/2)*(T-t))/sigma*sqrt(T-t)
#    d2 <- (log(St/K) + (r-sigma^2/2)*(T-t))/sigma*sqrt(T-t)
#    return(St*pnorm(d1,0,1) - K*exp(-r*(T-t))*pnorm(d2,0,1))
#}

#call_option_value_BS <- CallBS(St = St,
#                               K = strike,
#                               T = maturity,
#                               t = t,
#                               r = r,
#                               sigma = vola)


#diff_degree <- rep(0,10)
#for (d in 1:10) {
#  diff_degree[d] <- CallBS(St = St,K = strike,T = maturity,t = t,r = r,sigma = vola) - 
#                  LSM_american_option(state_variables = stock_prices, payoff = stock_prices, K = strike, dt = delta,rf = r,call = TRUE,orthogonal = "Laguerre",degree = d,verbose = FALSE)
#}
#matplot(1:10,diff_degree,"l")


#diff_delta <- rep(0,20)
#for (i in seq(10,200,10)) {
#  dt <- (maturity - t)/(250*i)
#  diff_delta[i/10] <- CallBS(St = St,K = strike,T = maturity,t = t,r = r,sigma = vola) - 
#    LSM_american_option(state_variables = stock_prices, payoff = stock_prices, K = strike, dt = dt,rf = r,call = TRUE,orthogonal = "Laguerre",degree = 3,verbose = FALSE)
#}
#lines(seq(10,200,10),diff_delta,col=palette.colors(7)[3])


#PutBS <- function(St,K,T,t,r,sigma){
#  d1 <- (log(St/K) + (r+sigma^2/2)*(T-t))/sigma*sqrt(T-t)
#  d2 <- (log(St/K) + (r-sigma^2/2)*(T-t))/sigma*sqrt(T-t)
#  return(K*exp(-r*(T-t))*pnorm(-d2,0,1) - St*pnorm(-d1,0,1))
#}

#put_option_value_BS <- PutBS(St = St,
#                             K = strike,
#                             T = maturity,
#                             t = t,
#                             r = r,
#                             sigma = vola)



### Rough Bergomi model ###

### parameters ###

m <- 10
rho <- 0.9
H <- 0.07
n<-100
xi <- 0.09
eta <- 1.9
S0 <- 100

### Computing the covariance matrix for the Volterra process Wtilde and ###
### a Brownian motion Z                                                 ###

CovMatrix <- matrix(0,m,m)
for (i in 1:m) {
  for (j in 1:m) {
    if (i > j) {
      CovMatrix[i,j] = rho*(sqrt(2*H)/H+0.5)*(i^(H+0.5)-(i-j)^(H+0.5))
    }
    else {
      CovMatrix[i,j] = rho*(sqrt(2*H)/H+0.5)*(i^(H+0.5))
    }
  }
}

### Get the Cholesky factorisation of CovMatrix ###

CovMatrix <- chol(CovMatrix)
CovMatrix <- t(CovMatrix)

### For each m, generate iid normal random vectors of length 2n ###

B <- matrix(rnorm(2*n*m,0,1),2*n,m)

### Generate the paths of Wtilde and Z with the correct joint marginals ###

pathsmatrix <- B %*% CovMatrix


V <- function(t){
  Vt <- rep(0,n)
  for (i in 1:n) {
    Vt[i] <- xi*exp(eta*pathsmatrix[i,t] - 0.5*eta^(2)*t^(2*H))
  }
  return(Vt)
}

S <- function(t){
  St <- rep(0,n)
  for (i in 1:n) {
    St[i] <- S0*exp(-0.5*V(t)[i]*t + sqrt(V(t)[i])*pathsmatrix[n+i,t])
  }
  return(St)
}

rbpath <- rep(0,m)
for (t in 1:m) {
  rbpath <- c(rbpath,S(t)[1])
}


### Generate fBM ###

#fbm <- function(hurst, n){
#  delta <- 1/n
#  r <- numeric(n+1)
#  r[1] <- 1
#  for(k in 1:n)
#    r[k+1] <- 0.5 * ((k+1)^(2*hurst) - 2*k^(2*hurst) + (k-1)^(2*hurst))
#  r <- c(r, r[seq(length(r)-1, 2)])
#  lambda <- Re((fft(r)) / (2*n))
#  W <- fft(sqrt(lambda) * (rnorm(2*n) + rnorm(2*n)*1i))
#  W <- n^(-hurst) * cumsum(Re(W[1:(n+1)]))
#  X <- ts(W, start=0, deltat=delta)
#  return(X)
#}

#matplot(seq(0,1,1/10000),fbm(0.3,10000),"l")


FBM(100,0.5,1)
N <- 1000
T <- 1
matplot(seq(0,T,T/N),c(0,FBM(N,0.5,T)),"l",col=palette.colors(3)[2])
lines(seq(0,T,T/N),BM(T,1,N),"l",col=palette.colors(3)[3])
#matplot(seq(0,1,1/N),BM(1,1,N),"l",col=palette.colors(3)[3])
#matplot(seq(0,2+1/N,1/N),BM(1,2,N),"l",col=palette.colors(3)[3])


numpaths <- 100
lenpaths <- (maturity - t)/((maturity - t)/(250*10))
delta<- (maturity - t)/(250*10)
stock_prices_GBM <- GBM_simulate(n = numpaths,
                                 t = maturity-t,
                                 mu = drift,
                                 sigma = vola,
                                 S0 = St,
                                 dt = delta)

stock_prices_FBM <- matrix(0,lenpaths+1,numpaths)
for (i in 1:numpaths) {
  stock_prices_FBM[,i] <- FBM(N=lenpaths+1,H=0.3,t=maturity-t)
}

### BM gives GBM as price process, it remains to clearify what fBM gives as 
### price process?


#stock_prices_GBM <- stock_prices_GBM[1:2500]
#stock_prices_FBM <- FBM(N=n,H=0.3,t=maturity-t)


H<-c(0.07,0.05)
T<-1
rho<-c(-0.9,1,-1)
xi0<-c(0.04,0.09)
nu<-c(1.7,1.23,1.52,1)


stock_prices_rBerg1 <- t(rough_bergomi(grid_points = 50,
                                       M=4000,
                                       H=setting1[1],
                                       T=maturity-t,
                                       rho=setting1[2],
                                       xi0=setting1[3],
                                       nu=setting1[4],
                                       S0=St))

LSM_american_option(state_variables = stock_prices_rBerg1, 
                    payoff = stock_prices_rBerg1, 
                    K = 90,
                    dt = (maturity - t)/50,
                    rf = r,
                    call = FALSE,
                    orthogonal = "Laguerre",
                    degree = 3,
                    verbose = FALSE)




### Fundamental Theorem of Derivatives Trading ###

#Vi_increments <- diff(optionprocess(S,K,T,r,sigma_i,call))
#Vh_increments <- diff(optionprocess(S,K,T,r,sigma_h,call))

#Pi <- rep(optionprocess(S,K,T,r,sigma_h,call)[1],length(S))

#for (t in (1:(length(S)-1))) {
#  Pi[t+1] = Pi[t] + Vi_increments[t] - Vh_increments[t] + (-r[t]*(optionprocess(S,K,T,r,sigma_i,call)[t] - optionprocess(S,K,T,r,sigma_h,call)[t]) + 0.5*((S[t])^2)*(((sigma_r[t])^2) - ((sigma_h[t])^2))*optiongamma(S,K,T,r,sigma_h)[t])*time_increment
#}

#pnl <- 0


#ru <- sum(r[1:t]*time_increment) 
#ru <- 0 
#for (u in 1:t) { 
#  ru <- ru + r[u]*time_increment
#}
#pnl <- pnl + (exp(-ru))*Pi_increments[t]
#pnl <- pnl + (exp(-sum(r[1:t]*time_increment)))*Pi_increments[t]
