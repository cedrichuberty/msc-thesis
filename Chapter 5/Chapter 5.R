### CHAPTER 5 ###

###  Includes (amongst others) the following objects:                   ###
###  The Least-Squares Monte Carlo algorithm with the code from         ###
###  subsection 5.2.2, Example 11, Fractional Brownian motion, rough    ###
###  Bergomi model, exponential shuffle                                 ###


### install packages ###

install.packages("xtable")
install.packages("hypergeo")
library(xtable)
library(hypergeo)

### The functions BM and BS from Chapter 2 and Chapter 4 must be ###
### initialized before executing this code.                      ###

### Own LSMC algorithm ###

### Simulate GBM ###

GBM <- function(n,T,S0,drift,diffusion){
  time_grid <- seq(0,T,T/(n-1))
  #time_increment <- time_grid[2]
  Brownian <- BM(T,1,n-1)
  #Brownianincrement <- diff(Brownian)
  gbm <- rep(S0,n)
  #for (t in (1:(n-1))) {
  #gbm[t+1] = gbm[t] + gbm[t]*(drift*time_increment + diffusion*Brownianincrement[t])
  #}
  for (t in (2:n)) {
    gbm[t] = S0*exp((drift-0.5*(diffusion)^2)*time_grid[t] + diffusion*Brownian[t])
  }
  return(gbm)
}

payoff <- function(x,strike,call=TRUE){
  return(ifelse(call,max(x-strike,0),max(strike-x,0)))
}

### The LSMC algorithm ###

LSMC <- function(paths,strike,r,basis="Monomial",call=TRUE){
  T <- dim(paths)[2] #this should be ex_opp...
  numpaths <- dim(paths)[1]
  #A <- round(seq(T-1,2,-((T-1)/(ex_opp-1))))
  scenario <- matrix(0,numpaths,T)
  #scenario <- matrix(0,numpaths,length(A)+2)
  for (i in (1:numpaths)) {
    scenario[i,T] <- payoff(paths[i,T],strike,call)
    #scenario[i,length(A)+2] <- payoff(paths[i,T],strike,call)
  }
  pb <- txtProgressBar(min = 0, max = T-3, style = 3, width = 50, char = "=")
  #pb <- txtProgressBar(min = 0, max = A[1]-A[length(A)], style = 3, width = 50, char = "=")
  for (t in ((T-1):2)) { #seq(T-1,2,-((T-1-2)/ex_opp))
    #for (t in round(seq(T-1,2,-((T-1)/(ex_opp-1))))){
    #for (t in A){
    X <- rep(0,numpaths)
    Y <- rep(0,numpaths)
    for (i in 1:numpaths) {
      if (payoff(paths[i,t],strike,call)<=0){scenario[i,t] <- 0}
      else{
        X[i] <- paths[i,t]
        Y[i] <- max(exp(-r*(1:(T-t)))*scenario[i,(t+1):T]) #still needs to be adapted for every t
        #Y[i] <- max(exp(-r*(1:((length(A)+2)-t)))*scenario[i,(t+1):(length(A)+2)])
        #Y[i] <- max(exp(-r*(-diff(c(rev(A[(match(t,rev(A))+1):length(A)]),T))))*scenario[i,c(A[(match(t,A)+1):length(A)],length(A)+2)]) #problem here, need to go through this with example
      }
    }
    Yprime <- Y[!X==0]
    Xprime <- X[!X==0]
    if (basis == "Monomial"){
      if (length(Xprime)>=9){regression <- as.numeric(lsfit(matrix(c(rep(1,length(Xprime)),Xprime,Xprime^2,Xprime^3,Xprime^4,Xprime^5,Xprime^6,Xprime^7),length(Xprime),8),Yprime)$coef)[1:8]}
      if (length(Xprime)==8){regression <- as.numeric(lsfit(matrix(c(rep(1,length(Xprime)),Xprime,Xprime^2,Xprime^3,Xprime^4,Xprime^5,Xprime^6),length(Xprime),7),Yprime)$coef)[1:7]}
      if (length(Xprime)==7){regression <- as.numeric(lsfit(matrix(c(rep(1,length(Xprime)),Xprime,Xprime^2,Xprime^3,Xprime^4,Xprime^5),length(Xprime),6),Yprime)$coef)[1:6]}
      if (length(Xprime)==6){regression <- as.numeric(lsfit(matrix(c(rep(1,length(Xprime)),Xprime,Xprime^2,Xprime^3,Xprime^4),length(Xprime),5),Yprime)$coef)[1:5]}
      if (length(Xprime)==5){regression <- as.numeric(lsfit(matrix(c(rep(1,length(Xprime)),Xprime,Xprime^2,Xprime^3),length(Xprime),4),Yprime)$coef)[1:4]}
      if (length(Xprime)==4){regression <- as.numeric(lsfit(matrix(c(rep(1,length(Xprime)),Xprime,Xprime^2),length(Xprime),3),Yprime)$coef)[1:3]}
      if (length(Xprime)==3){regression <- as.numeric(lsfit(matrix(c(rep(1,length(Xprime)),Xprime),length(Xprime),2),Yprime)$coef)[1:2]}
      if (length(Xprime)==2){regression <- as.numeric(lsfit(matrix(rep(1,length(Xprime)),length(Xprime),1),Yprime)$coef)[1]}
      #else if (length(Xprime)==1){regression <- Yprime/Xprime}
      keep <- rep(0,numpaths)
      for (i in (1:numpaths)){
        if (payoff(paths[i,t],strike,call)>0){
          if (length(Xprime)>=9){keep[i] <- max(as.numeric(regression%*%c(1,X[i],(X[i])^2,(X[i])^3,(X[i])^4,(X[i])^5,(X[i])^6,(X[i])^7)),0)} #round yes or no? makes a huge difference
          if (length(Xprime)==8){keep[i] <- max(as.numeric(regression%*%c(1,X[i],(X[i])^2,(X[i])^3,(X[i])^4,(X[i])^5,(X[i])^6)),0)}
          if (length(Xprime)==7){keep[i] <- max(as.numeric(regression%*%c(1,X[i],(X[i])^2,(X[i])^3,(X[i])^4,(X[i])^5)),0)} 
          if (length(Xprime)==6){keep[i] <- max(as.numeric(regression%*%c(1,X[i],(X[i])^2,(X[i])^3,(X[i])^4)),0)}
          if (length(Xprime)==5){keep[i] <- max(as.numeric(regression%*%c(1,X[i],(X[i])^2,(X[i])^3)),0)}
          if (length(Xprime)==4){keep[i] <- max(as.numeric(regression%*%c(1,X[i],(X[i])^2)),0)} 
          if (length(Xprime)==3){keep[i] <- max(as.numeric(regression%*%c(1,X[i])),0)}
          if (length(Xprime)==2){keep[i] <- max(regression,0)}
          #else if (length(Xprime)==3){keep[i] <- max(as.numeric(round(regression,2)%*%c(1,X[i])),0)}
          #else if (length(Xprime)==2){keep[i] <- max(round(regression,2),0)}
          if (length(Xprime) >= 2){
            if (keep[i]>payoff(paths[i,t],strike,call)){
              scenario[i,t] <- 0 #keep[i]  
            }
            else{
              scenario[i,t] <- max(keep[i],payoff(paths[i,t],strike,call)) #max here necessary? could be payoff() directly
              scenario[i,(t+1):T] <- 0
              #scenario[i,(t+1):(length(A)+2)] <- 0
              #scenario[i,c(A[(match(t,A)+1):length(A)],length(A)+2)] <- 0
            }
          }
        }
      }
    }
    else if (basis == "Laguerre"){
      #if (length(Xprime)>=7){regression <- as.numeric(lsfit(matrix(c(exp(-Xprime/2),exp(-Xprime/2)*(1-Xprime),exp(-Xprime/2)*(1-2*Xprime+((Xprime)^2)/2),exp(-Xprime/2)*(1-3*Xprime+(3/2)*(Xprime)^2-(1/6)*(Xprime)^3),exp(-Xprime/2)*(1-4*Xprime+3*(Xprime)^2-(2/3)*(Xprime)^3+(1/24)*(Xprime)^4),exp(-Xprime/2)*(1-5*Xprime+5*(Xprime)^2-(5/3)*(Xprime)^3+(5/24)*(Xprime)^4-(1/120)*(Xprime)^5)),length(Xprime),6),Yprime)$coef)[1:6]}
      #else if (length(Xprime)==6){regression <- as.numeric(lsfit(matrix(c(exp(-Xprime/2),exp(-Xprime/2)*(1-Xprime),exp(-Xprime/2)*(1-2*Xprime+((Xprime)^2)/2),exp(-Xprime/2)*(1-3*Xprime+(3/2)*(Xprime)^2-(1/6)*(Xprime)^3),exp(-Xprime/2)*(1-4*Xprime+3*(Xprime)^2-(2/3)*(Xprime)^3+(1/24)*(Xprime)^4)),length(Xprime),5),Yprime)$coef)[1:5]}
      #else if (length(Xprime)==5){regression <- as.numeric(lsfit(matrix(c(exp(-Xprime/2),exp(-Xprime/2)*(1-Xprime),exp(-Xprime/2)*(1-2*Xprime+((Xprime)^2)/2),exp(-Xprime/2)*(1-3*Xprime+(3/2)*(Xprime)^2-(1/6)*(Xprime)^3)),length(Xprime),4),Yprime)$coef)[1:4]}
      #else if (length(Xprime)==4){regression <- as.numeric(lsfit(matrix(c(exp(-Xprime/2),exp(-Xprime/2)*(1-Xprime),exp(-Xprime/2)*(1-2*Xprime+((Xprime)^2)/2)),length(Xprime),3),Yprime)$coef)[1:3]}
      #else if (length(Xprime)==3){regression <- as.numeric(lsfit(matrix(c(exp(-Xprime/2),exp(-Xprime/2)*(1-Xprime)),length(Xprime),2),Yprime)$coef)[1:2]}
      #else if (length(Xprime)==2){regression <- as.numeric(lsfit(matrix(exp(-Xprime/2),length(Xprime),1),Yprime)$coef)[1]}
      ##else if (length(Xprime)==1){regression <- Yprime/Xprime}
      if (length(Xprime)>=7){regression <- as.numeric(lsfit(matrix(c(rep(1,length(Xprime)),1-Xprime,1-2*Xprime+((Xprime)^2)/2,1-3*Xprime+(3/2)*(Xprime)^2-(1/6)*(Xprime)^3,1-4*Xprime+3*(Xprime)^2-(2/3)*(Xprime)^3+(1/24)*(Xprime)^4,1-5*Xprime+5*(Xprime)^2-(5/3)*(Xprime)^3+(5/24)*(Xprime)^4-(1/120)*(Xprime)^5),length(Xprime),6),Yprime)$coef)[1:6]}
      else if (length(Xprime)==6){regression <- as.numeric(lsfit(matrix(c(rep(1,length(Xprime)),1-Xprime,1-2*Xprime+((Xprime)^2)/2,1-3*Xprime+(3/2)*(Xprime)^2-(1/6)*(Xprime)^3,1-4*Xprime+3*(Xprime)^2-(2/3)*(Xprime)^3+(1/24)*(Xprime)^4),length(Xprime),5),Yprime)$coef)[1:5]}
      else if (length(Xprime)==5){regression <- as.numeric(lsfit(matrix(c(rep(1,length(Xprime)),1-Xprime,1-2*Xprime+((Xprime)^2)/2,1-3*Xprime+(3/2)*(Xprime)^2-(1/6)*(Xprime)^3),length(Xprime),4),Yprime)$coef)[1:4]}
      else if (length(Xprime)==4){regression <- as.numeric(lsfit(matrix(c(rep(1,length(Xprime)),1-Xprime,1-2*Xprime+((Xprime)^2)/2),length(Xprime),3),Yprime)$coef)[1:3]}
      else if (length(Xprime)==3){regression <- as.numeric(lsfit(matrix(c(rep(1,length(Xprime)),1-Xprime),length(Xprime),2),Yprime)$coef)[1:2]}
      else if (length(Xprime)==2){regression <- as.numeric(lsfit(matrix(rep(1,length(Xprime)),length(Xprime),1),Yprime)$coef)[1]}
      keep <- rep(0,numpaths)
      for (i in (1:numpaths)){
        if (payoff(paths[i,t],strike,call)>0){
          #if (length(Xprime)>=7){keep[i] <- max(as.numeric(regression%*%c(exp(-(X[i])/2),exp(-(X[i])/2)*(1-X[i]),exp(-(X[i])/2)*(1-2*X[i]+((X[i])^2)/2),exp(-(X[i])/2)*(1-3*X[i]+(3/2)*(X[i])^2-(1/6)*(X[i])^3),exp(-(X[i])/2)*(1-4*X[i]+3*(X[i])^2-(2/3)*(X[i])^3+(1/24)*(X[i])^4),exp(-(X[i])/2)*(1-5*X[i]+5*(X[i])^2-(5/3)*(X[i])^3+(5/24)*(X[i])^4-(1/120)*(X[i])^5))),0)} #round yes or no? makes a huge difference
          #else if (length(Xprime)==6){keep[i] <- max(as.numeric(regression%*%c(exp(-(X[i])/2),exp(-(X[i])/2)*(1-X[i]),exp(-(X[i])/2)*(1-2*X[i]+((X[i])^2)/2),exp(-(X[i])/2)*(1-3*X[i]+(3/2)*(X[i])^2-(1/6)*(X[i])^3),exp(-(X[i])/2)*(1-4*X[i]+3*(X[i])^2-(2/3)*(X[i])^3+(1/24)*(X[i])^4))),0)}
          #else if (length(Xprime)==5){keep[i] <- max(as.numeric(regression%*%c(exp(-(X[i])/2),exp(-(X[i])/2)*(1-X[i]),exp(-(X[i])/2)*(1-2*X[i]+((X[i])^2)/2),exp(-(X[i])/2)*(1-3*X[i]+(3/2)*(X[i])^2-(1/6)*(X[i])^3))),0)}
          #else if (length(Xprime)==4){keep[i] <- max(as.numeric(regression%*%c(exp(-(X[i])/2),exp(-(X[i])/2)*(1-X[i]),exp(-(X[i])/2)*(1-2*X[i]+((X[i])^2)/2))),0)}
          #else if (length(Xprime)==3){keep[i] <- max(as.numeric(regression%*%c(exp(-(X[i])/2),exp(-(X[i])/2)*(1-X[i]))),0)}
          #else if (length(Xprime)==2){keep[i] <- max(as.numeric(regression%*%c(exp(-(X[i])/2))),0)}
          if (length(Xprime)>=7){keep[i] <- max(as.numeric(regression%*%c(1,1-X[i],1-2*X[i]+((X[i])^2)/2,1-3*X[i]+(3/2)*(X[i])^2-(1/6)*(X[i])^3,1-4*X[i]+3*(X[i])^2-(2/3)*(X[i])^3+(1/24)*(X[i])^4,1-5*X[i]+5*(X[i])^2-(5/3)*(X[i])^3+(5/24)*(X[i])^4-(1/120)*(X[i])^5)),0)} #round yes or no? makes a huge difference
          else if (length(Xprime)==6){keep[i] <- max(as.numeric(regression%*%c(1,1-X[i],1-2*X[i]+((X[i])^2)/2,1-3*X[i]+(3/2)*(X[i])^2-(1/6)*(X[i])^3,1-4*X[i]+3*(X[i])^2-(2/3)*(X[i])^3+(1/24)*(X[i])^4)),0)}
          else if (length(Xprime)==5){keep[i] <- max(as.numeric(regression%*%c(1,1-X[i],1-2*X[i]+((X[i])^2)/2,1-3*X[i]+(3/2)*(X[i])^2-(1/6)*(X[i])^3)),0)}
          else if (length(Xprime)==4){keep[i] <- max(as.numeric(regression%*%c(1,1-X[i],1-2*X[i]+((X[i])^2)/2)),0)}
          else if (length(Xprime)==3){keep[i] <- max(as.numeric(regression%*%c(1,1-X[i])),0)}
          else if (length(Xprime)==2){keep[i] <- max(as.numeric(regression%*%c(1)),0)}
          if (length(Xprime) >= 2){
            if (keep[i]>payoff(paths[i,t],strike,call)){
              scenario[i,t] <- 0 #keep[i]  
            }
            else{
              scenario[i,t] <- max(keep[i],payoff(paths[i,t],strike,call))
              scenario[i,(t+1):T] <- 0
              #scenario[i,(t+1):(length(A)+2)] <- 0
              #scenario[i,c(A[(match(t,A)+1):length(A)],length(A)+2)] <- 0
            }
          }
        }
      }
    }
    setTxtProgressBar(pb, T-1-t)
    #setTxtProgressBar(pb, A[1]-t)
  }
  close(pb)
  for (i in 1:numpaths) {
    scenario[i,1] <- max(exp(-r*(1:(T-1)))*scenario[i,2:T])
    #scenario[i,1] <- max(exp(-r*(1:((length(A)+2)-1)))*scenario[i,2:(length(A)+2)])
  }
  return((1/numpaths)*sum(scenario[,1]))
  #return(scenario)
}

### Part presented in subsection 5.2.2. of the thesis ### 

S0 <- 100
K <- 95
T <- 1
r <- 0.06
sigma <- 0.3
paths_to_be_simulated <- 10
n <- 4 
strike <- K
call <- TRUE

BM(T,1,n)

GBM(n,T,S0,r,sigma)

stock_prices <- matrix(0,paths_to_be_simulated,n)
for (i in 1:paths_to_be_simulated) {
  stock_prices[i,] <- GBM(n,T,S0,r,sigma)
}

print(xtable(round(stock_prices,2),type="latex"),file = "stock_prices.tex")

#Hardcoding the paths I got so that one can redo the steps:

#stock_prices_1 <- matrix(c(100,100,100,100,100,100,100,100,100,100,
#                           132.47,134.17,111.60,77.99,83.48,108.81,85.92,105.11,94.81,84.68,
#                           131.15,141.91,90.05,104.71,84.24,90.22,103.37,117.81,75.88,60.91,
#                           107.52,165.03,93.26,107.09,99.38,90.26,113.34,97.37,63.06,46.00),10,4)

paths <- round(stock_prices,2)

T <- dim(paths)[2] 
numpaths <- dim(paths)[1]
scenario <- matrix(0,numpaths,T)


for (i in (1:numpaths)) {
  scenario[i,T] <- payoff(paths[i,T],strike,call)
}

print(xtable(scenario,type="latex"),file = "scenariot3.tex")

t<-T-1

X <- rep(0,numpaths)
Y <- rep(0,numpaths)
for (i in 1:numpaths) {
  if (payoff(paths[i,t],strike,call)<=0){scenario[i,t] <- 0}
  else{
    X[i] <- paths[i,t]
    Y[i] <- max(exp(-r*(1:(T-t)))*scenario[i,(t+1):T]) 
  }
}

print(xtable(matrix(c(X,Y),10,2),type="latex"),file = "XY3.tex")

Yprime <- Y[!X==0]
Xprime <- X[!X==0]

if (length(Xprime)>=4){regression <- as.numeric(lsfit(matrix(c(rep(1,length(Xprime)),Xprime,Xprime^2),length(Xprime),3),Yprime)$coef)[1:3]}
if (length(Xprime)==3){regression <- as.numeric(lsfit(matrix(c(rep(1,length(Xprime)),Xprime),length(Xprime),2),Yprime)$coef)[1:2]}
if (length(Xprime)==2){regression <- as.numeric(lsfit(matrix(rep(1,length(Xprime)),length(Xprime),1),Yprime)$coef)[1]}

regression

keep <- rep(0,numpaths)
for (i in (1:numpaths)){
  if (payoff(paths[i,t],strike,call)>0){
    if (length(Xprime)>=4){keep[i] <- max(as.numeric(regression%*%c(1,X[i],(X[i])^2)),0)} 
    if (length(Xprime)==3){keep[i] <- max(as.numeric(regression%*%c(1,X[i])),0)}
    if (length(Xprime)==2){keep[i] <- max(regression,0)}
    if (length(Xprime) >= 2){
      if (keep[i]>payoff(paths[i,t],strike,call)){
        scenario[i,t] <- 0  
      }
      else{
        scenario[i,t] <- payoff(paths[i,t],strike,call) 
        scenario[i,(t+1):T] <- 0
      }
    }
  }
}

print(xtable(matrix(keep,10,1),type="latex"),file = "keep3.tex")
print(xtable(scenario,type="latex"),file = "scenariot2.tex")

t <- t-1

X <- rep(0,numpaths)
Y <- rep(0,numpaths)
for (i in 1:numpaths) {
  if (payoff(paths[i,t],strike,call)<=0){scenario[i,t] <- 0}
  else{
    X[i] <- paths[i,t]
    Y[i] <- max(exp(-r*(1:(T-t)))*scenario[i,(t+1):T]) 
  }
}

print(xtable(matrix(c(X,Y),10,2),type="latex"),file = "XY2.tex")

Yprime <- Y[!X==0]
Xprime <- X[!X==0]

if (length(Xprime)>=4){regression <- as.numeric(lsfit(matrix(c(rep(1,length(Xprime)),Xprime,Xprime^2),length(Xprime),3),Yprime)$coef)[1:3]}
if (length(Xprime)==3){regression <- as.numeric(lsfit(matrix(c(rep(1,length(Xprime)),Xprime),length(Xprime),2),Yprime)$coef)[1:2]}
if (length(Xprime)==2){regression <- as.numeric(lsfit(matrix(rep(1,length(Xprime)),length(Xprime),1),Yprime)$coef)[1]}

regression

keep <- rep(0,numpaths)
for (i in (1:numpaths)){
  if (payoff(paths[i,t],strike,call)>0){
    if (length(Xprime)>=4){keep[i] <- max(as.numeric(regression%*%c(1,X[i],(X[i])^2)),0)} 
    if (length(Xprime)==3){keep[i] <- max(as.numeric(regression%*%c(1,X[i])),0)}
    if (length(Xprime)==2){keep[i] <- max(regression,0)}
    if (length(Xprime) >= 2){
      if (keep[i]>payoff(paths[i,t],strike,call)){
        scenario[i,t] <- 0  
      }
      else{
        scenario[i,t] <- payoff(paths[i,t],strike,call) 
        scenario[i,(t+1):T] <- 0
      }
    }
  }
}

print(xtable(matrix(keep,10,1),type="latex"),file = "keep2.tex")
print(xtable(scenario,type="latex"),file = "scenariot1.tex")

for (i in 1:numpaths) {
  scenario[i,1] <- max(exp(-r*(1:(T-1)))*scenario[i,2:T])
}

print(xtable(matrix(scenario[,1],10,1),type="latex"),file = "scenariot0.tex")

(1/numpaths)*sum(scenario[,1])

BS(S0,K,T,0,r,sigma,TRUE)

### Price given by the LSMC for the Netflix stock (Example 11 in the thesis) ###
### one first needs to calculate the historic volatility in Chapter 4     ###

S0 <- 240.13
K <- 300
r <- 0.0467214 #0.02263
sigma <- historic_vola
T <- 64/251
n <- 4 #8 
paths_to_be_simulated <- 100000 

synthetic_prices <- matrix(0,paths_to_be_simulated,n)
for (i in 1:paths_to_be_simulated) {
  synthetic_prices[i,] <- GBM(n,T,S0,r,sigma)
}

LSMC(synthetic_prices,K,r,"Laguerre",TRUE)

LSMC(synthetic_prices,K,r,"Laguerre",FALSE)


### Fractional Brownian Motion (taken from the Github fractionalBM) ###

FBM <- function(N, H, t=1, method='davies-harte') {
  
  # Hurst parameter check
  if((H < 0) | (H > 1)) {
    stop("Hurst parameter must be between 0 and 1.")
  }
  
  # autocovariance function
  gamma <- function(k, h) {
    return(0.5*( (abs(k-1)^(2*h)) - (2*(abs(k)^(2*h))) + (abs(k+1)^(2*H))))
  }
  
  # Davies Harte method:
  if(method == 'davies-harte') {
    
    # Get eigenvalues
    g = c()
    for(k in 0:(N-1)){g <- c(g, gamma(k,H))}
    r = c(g, 0, rev(g)[1:(N-1)])
    j = seq(0, ((2*N)-1))
    K = (2*N)-1
    i = complex(real=0, imaginary=1)
    lk = rev(fft(r*exp(2*pi*i*K*j*(1/(2*N)))))
    
    # Generate random variables
    Vj <- cbind(seq(0,0,length.out=2*N), seq(0,0,length.out=2*N))
    Vj[1,1] <- rnorm(1)
    Vj[N+1, 1] <- rnorm(1)
    Vj1 <- rnorm(N-1)
    Vj2 <- rnorm(N-1)
    Vj[2:N,1] <- Vj1
    Vj[2:N,2] <- Vj2
    Vj[(N+2):(2*N),1] <- rev(Vj1)
    Vj[(N+2):(2*N),2] <- rev(Vj2)
    
    # Compute Z (fractional Gaussian Noise)
    wk = seq(0,0,length.out=2*N)
    wk[1] <- sqrt(lk[1]/(2*N))*Vj[1,1]
    wk[2:N] <- sqrt(lk[2:N]/(4*N))*(Vj[2:N,1] + i*Vj[2:N,2])
    wk[N+1] <- sqrt(lk[N+1]/(2*N))*Vj[N+1,1]
    wk[(N+2):(2*N)] <- sqrt(lk[(N+2):(2*N)]/(4*N))*(Vj[(N+2):(2*N),1] - i*Vj[(N+2):(2*N),2])
    Z = fft(wk)
    fGn = Z[1:N]
    fBm = cumsum(fGn)*(N^(-H))
    return(Re((t^H)*fBm))
    
    
  } else if(method == 'hosking') {
    
    # Starting values for the recursion
    X = c(rnorm(1)) # vector containing fractional Gaussian Noise
    mu = c(gamma(1,H)*X[1])
    sigsq = c(1-(gamma(1,H)^2))
    tau = c(gamma(1,H)^2)
    d = c(gamma(1,H))
    
    for(n in 2:N){
      Fmat = apply(diag(n),2,rev)
      cn = c()
      for(k in 0:(n-1)){cn <- c(cn, gamma(k+1,H))}
      
      # Generate sig(n+1)^2
      s = sigsq[n-1] - ((gamma(n,H) - tau[n-1])^2)/sigsq[n-1]
      sigsq = c(sigsq, s)
      
      # Generate d(n+1)
      phi = (gamma(n,H) - tau[n-1])/sigsq[n-1]
      d = d - phi*rev(d)
      d = c(d, phi)
      
      # mu(n+1) and tau(n+1)
      Xnplus1 = mu[n-1] + sigsq[n-1]*rnorm(1)
      X = c(X, Xnplus1)
      mu = c(mu, d %*% rev(X))
      tau = c(tau, cn %*% Fmat %*% d)
    }
    
    fBm = cumsum(X)*(N^(-H))
    return((t^H)*fBm)
    
    
  } else if(method == 'cholesky') {
    
    # Starting values for recursion
    L = matrix(0, N,N)
    X = seq(0,0,length.out=N)
    V = rnorm(N)
    
    L[1,1] = 1.0
    X[1] = V[1]
    
    L[2,1] = gamma(1,H)
    L[2,2] = sqrt(1 - (L[2,1]^2))
    X[2] = sum(L[2,1:2])
    
    # Generate Cholesky matrix by row
    for(i in 3:N) {
      
      # First row value
      L[i,1] = gamma(i-1,H)
      
      # Middle values
      for(j in 2:(i-1)){
        L[i,j] = (1/L[j,j])*(gamma(i-j,H) - (L[i,1:(j-1)]%*%(L[j,1:(j-1)])))
      }
      
      # Last value
      L[i,i] = sqrt(1 - sum(L[i,1:(i-1)]^2))
      
      # Use row to compute X(n+1)
      X[i] = L[i,1:i] %*% V[1:i]
    }
    
    fBm = cumsum(X)*(N^(-H))
    return((t^H)*fBm)
    
  } else {stop("Incorrect method")}
}

### Plots in the thesis ###

#plot(seq(0,1,1/(250-1)),c(0,FBM(250-1,1/3,1)),type="l",col="blue",ylim = c(-2,2),main = "fBM - BM",col.main="black",xlab = "time t",ylab = "")
#lines(seq(0,1,1/(250-1)),c(0,FBM(250-1,2/3,1)),col="purple")
#lines(seq(0,1,1/(250-1)),BM(1,1,250-1),col="red")
#for (i in 2:3) {
#  lines(seq(0,1,1/(250-1)),BM(1,1,250-1),col="red")
#  lines(seq(0,1,1/(250-1)),c(0,FBM(250-1,1/3,1)),col="blue") 
#  lines(seq(0,1,1/(250-1)),c(0,FBM(250-1,2/3,1)),col="purple") 
#}

plot(seq(0,1,1/(250-1)),c(0,FBM(250-1,1/3,1)),type="l",col="blue",ylim=c(-2,2),main = "Fractional Brownian motion with H=1/3",col.main="black",xlab = "time t",ylab = "")
for (i in 2:3) {
  lines(seq(0,1,1/(250-1)),c(0,FBM(250-1,1/3,1)),col="blue") 
}

plot(seq(0,1,1/(250-1)),c(0,FBM(250-1,2/3,1)),type="l",col="purple",ylim=c(-2,2),main = "Fractional Brownian motion with H=2/3",col.main="black",xlab = "time t",ylab = "")
for (i in 2:3) {
  lines(seq(0,1,1/(250-1)),c(0,FBM(250-1,2/3,1)),col="purple") 
}

plot(seq(0,1,1/(250-1)),BM(1,1,250-1),type="l",col="red",ylim=c(-2,2),main = "Brownian motion",col.main="black",xlab = "time t",ylab = "")
for (i in 2:3) {
  lines(seq(0,1,1/(250-1)),BM(1,1,250-1),col="red") 
}

### rough Bergomi model (followed the approach from the Github Marketsimulator) ###

volterra_BM_path_chol <- function(grid_points,M,H,T=1,rho){
  
  X <- seq(0,T,T/(grid_points-1))
  
  X <- X[2:grid_points]
  
  size <- 2*(grid_points-1)
  Sigma <- matrix(0,size,size)
  
  for (j in (1:(grid_points-1))) {
    for (i in (1:(grid_points-1))) {
      if (i==j){
        Sigma[i,j]=(X[i]^(2*H))/(2*H)
      } else{
        s<-min(X[i],X[j])
        t<-max(X[i],X[j])
        Sigma[i,j] <- ((t-s)^(H-0.5))*(1/(H+0.5))*(s^(0.5+H))*hypergeo(0.5-H,0.5+H,1.5+H,(-s)/(t-s)) 
      }
    }
  }
  
  for (j in (1:(grid_points-1))) {
    for (i in (1:(grid_points-1))) {
      Sigma[i,j+(grid_points-1)] = rho*(1/(H+0.5))*((X[i]^(H+0.5))-((X[i]-min(X[i],X[j]))^(H+0.5)))
      Sigma[i+(grid_points-1),j] = rho*(1/(H+0.5))*((X[j]^(H+0.5))-((X[j]-min(X[i],X[j]))^(H+0.5)))
    }
  }
  
  for (j in (1:(grid_points-1))) {
    for (i in (1:(grid_points-1))) {
      Sigma[i+(grid_points-1),j+(grid_points-1)]=min(X[i],X[j])
    }
  }
  
  P<-t(chol(Re(Sigma)))
  
  #Z<-matrix(0,M,2*(grid_points-1)) 
  #for (i in 1:M) {
  #  for (j in 1:(2*(grid_points-1))) {
  #    Z[i,j] <- rnorm(1,0,1)
  #  }
  #}
  Z<-matrix(rnorm(M*2*(grid_points-1),0,1),M,2*(grid_points-1))
  
  V <- matrix(0,M,grid_points)
  W <- matrix(0,M,grid_points)
  
  for (i in (1:M)) {
    aux <- P%*%Z[i,] 
    V[i,2:grid_points] <- aux[1:(grid_points-1)]
    W[i,2:grid_points] <- aux[grid_points:(2*(grid_points-1))]
  }
  return(array(cbind(V,W),c(M,grid_points,2))) 
}

rough_bergomi <- function(grid_points,M,H,T=1,rho,xi0,nu,S0){
  
  path <- volterra_BM_path_chol(grid_points,M,H,T,rho)
  V_path <- path[,,1]
  W_path <- path[,,2]
  
  time_grid <- seq(0,T,T/(grid_points-1))
  C_H <- sqrt((2*H*gamma(3/2-H))/(gamma(H+0.5)*gamma(2-2*H)))
  Variance <- xi0*exp(2*nu*C_H*V_path - (((nu*C_H)^2)*((time_grid)^(2*H))*(1/H)))
  log_stock <- matrix(1,M,grid_points)*log(S0)
  time_increment <- time_grid[2]
  
  if (M>1){
    brownian_increments <- matrix(0,M,grid_points-1)
    for (i in (1:M)) {
      brownian_increments[i,] <- diff(W_path[i,])
    }
    
    for (i in (1:(grid_points-1))) {
      log_stock[,i+1] <- log_stock[,i] - 0.5*Variance[,i]*time_increment + sqrt(Variance[,i])*brownian_increments[,i]
    }
  }else{
    brownian_increments<-diff(W_path)
    
    for (i in (1:(grid_points-1))) {
      log_stock[i+1] <- log_stock[i] - 0.5*Variance[i]*time_increment + sqrt(Variance[i])*brownian_increments[i]
    }
  }
  
  return(exp(log_stock))
}

### Plots in the thesis ###

setting1<-c(0.07,-0.9,0.09,1.23) #setting 1 in "Pricing under rough volatility"
setting2<-c(0.05,-0.9,0.09,1.52) #setting 2 in "Pricing under rough volatility"

rBergomi <- t(rough_bergomi(grid_points = 250,
                            M=5,
                            H=setting1[1],
                            T=1,
                            rho=setting1[2],
                            xi0=setting1[3],
                            nu=setting1[4],
                            S0=100))

#plot(rBergomi[,1],type="l",col="blue",ylim=c(50,150))
#lines(GBM(250,1,100,0,0.3),col="red")
#for (i in 2:3) {
#  lines(GBM(250,1,100,0,0.3),col="red")
#  lines(rBergomi[,i],col="blue") 
#}

plot(seq(0,1,1/(250-1)),rBergomi[,1],type="l",col="blue",ylim=c(50,150),main = "Simulated price process in the rough Bergomi model",col.main="black",xlab = "time t",ylab = "Stock price")
for (i in 2:5) {
  lines(seq(0,1,1/(250-1)),rBergomi[,i],col="blue") 
}

plot(seq(0,1,1/(250-1)),GBM(250,1,100,0,0.3),type="l",col="red",ylim=c(50,150),main = "Price process simulated by geometric Brownian motion",col.main="black",xlab = "time t",ylab = "Stock price")
for (i in 2:5) {
  lines(seq(0,1,1/(250-1)),GBM(250,1,100,0,0.3),col="red") 
}


### For the next part, the functions from Chapter 2 need to be run first ###

### Exponential shuffle ###

exponential_shuffle <- function(l,truncorder,simple = FALSE){ 
  res <- paste(1,"")
  pb <- txtProgressBar(min = 0, max = length(l), style = 3, width = 50, char = "=")
  for (i in 1:length(l)) {
    argument <- l[i]
    #truncorder <- nchar(strsplit(argument," ")[[1]][2]) #maybe this is wrong, but don't need to change anything, just comment it out and put truncorder in the arguments of the function
    word <- ""
    int <- c()
    for (j in 1:truncorder) {
      word <- shuffle(word,strsplit(argument," ")[[1]][2]) 
      int <- c(int,paste(((as.numeric(strsplit(argument," ")[[1]][1]))^j)/factorial(j),word)) #if (nchar(word) <= truncorder) {int <- ...} maybe this is unnecessary if I have a condition before concatenation already
    }
    int <- c(paste(1,""),int)
    int_res <- res
    t <- c()
    for (k in 1:length(int_res)) {
      for (h in 1:length(int)) {
        lenA <- ifelse(is.na(strsplit(int_res[k]," ")[[1]][2]),0,nchar(strsplit(int_res[k]," ")[[1]][2])) 
        lenB <- ifelse(is.na(strsplit(int[h]," ")[[1]][2]),0,nchar(strsplit(int[h]," ")[[1]][2]))
        if ((lenA + lenB) <= truncorder){
          t <- c(t,concatenation(int_res[k],int[h])) 
        }
      }
    }
    res <- t
    setTxtProgressBar(pb, i)
  }
  close(pb)
  
  ifelse(simple,return(simplify(res)),return(res))
}

### LHS in the pairing in Corollary 6.6 in ###
### "Optimal stopping with signatures"     ###

stopping <- function(l,truncorder){
  stopp <- shuffle2(l,l)
  for (i in 1:length(stopp)) {
    stopp[i] <- paste(stopp[i],1,sep="")
    stopp[i] <- paste((-1)*as.numeric(strsplit(stopp[i]," ")[[1]][1]),strsplit(stopp[i]," ")[[1]][2])
  }
  stopp <- exponential_shuffle(simplify(stopp),truncorder)
  for (i in 1:length(stopp)) {
    stopp[i] <- paste(stopp[i],2,sep="")
  }
  return(stopp)
}

### Example in the thesis ###

a <- paste(1,"")
length(stopping(a,1))
length(stopping(a,2))
length(stopping(a,3))
length(stopping(a,4))
length(stopping(a,5))
length(stopping(a,6))
length(stopping(a,7))

### l1-norm for words ###

l1_words <- function(a){
  l1 <- 0
  for (w in 1:length(a)) {
    l1 <- l1 + abs(as.numeric(strsplit(a[w]," ")[[1]][1]))
  }
  return(l1)
}

### The beginning of how to compute Corollary 6.6 in "Optimal stopping  ###
### with signatures" but we did not elaborate this as it quickly became ###
### clear that it won't be feasible                                     ###

optimalstopping <- function(pricepath,truncorder){
  dimpath <- dim(pricepath)[1]
  lenpath <- dim(pricepath)[2]
  numpath <- dim(pricepath)[3]
  ifelse(dimpath==1,signa <- array(0,c(1,1+truncorder,numpath)),signa <- array(0,c(1,((dimpath*((dimpath^truncorder)-1))/(dimpath-1))+1,numpath)))
  for (i in (1:numpath)) {
    signa[,,i] <- c(1,sig(array(pricepath[,,i],c(dimpath,lenpath)),truncorder))
  }
  expsigna <- exp_sig(signa)
  
  #optimizing (not elaborated)
  for (K in 1:2) {
    for (N in 1:4) {
      word <- c()
      for (l in keys(3,N)) {
        a <- keys(3,N)
        for (i in 1:length(keys(3,N))) {
          for (k in seq(-K,K,1/10)) {
            l <- paste(k,a[i])
          }
          word <- c(word,l)
          #rest of bazar
        }
        #word <- c(word,paste(runif(1,-K,K),l))
      }
      if ((l1_words(word)+degree_of_word(word)) <= K){
        pairing(stopping(word,truncorder),expsigna,keys(3,N))
      } #or something like use the distributivity of scalar product and take out coefficents so that you end up with only the coefficients as unknowns and then use Linear Programming for each l(=combination of keys without coefficients)
    }
  }
  pairing(stopping(l,truncorder),expsigna,key)
  return()
}


### Comparison of the rough Bergomi model with the tables 4 and 5 (paragraph 8.2) ###
### in "Machine Learning for Pricing American Options in High-Dimensional         ###
### Markovian and non-Markovian models"                                           ###

N <- c(50,100)
K <- c(70,80,90,100,110,120,130,140)
P <- c(500,1000,2000,4000,8000)
Ndescr <- c("50 time steps","100 time steps")
Kdescr <- c("Strike=70","Strike=80","Strike=90","Strike=100","Strike=110","Strike=120","Strike=130","Strike=140")
Pdescr <- c("500 paths","1000 paths","2000 paths","4000 paths","8000 paths")
res_comp <- array(0,c(8,5,2),dimnames = list(Kdescr,Pdescr,Ndescr))
for (k in (1:2)) {
  for (i in (1:8)) {
    for (j in (1:5)) {
      stock_prices_rBerg1 <- t(rough_bergomi(grid_points = N[k],
                                             M=P[j],
                                             H=setting1[1],
                                             T=1,
                                             rho=setting1[2],
                                             xi0=setting1[3],
                                             nu=setting1[4],
                                             S0=100))
      
      res_comp[i,j,k] <- LSM_american_option(state_variables = stock_prices_rBerg1, 
                                             payoff = stock_prices_rBerg1, 
                                             K = K[i],
                                             dt = 1/N[k],
                                             rf = 0.05,
                                             call = FALSE,
                                             orthogonal = "Laguerre",
                                             degree = 3,
                                             verbose = FALSE)
    }
  }
}
round(res_comp,3)

GPR_Tree <- array(c(1.87,3.18,5.24,8.36,13.04,20.19,30.00,40.00,
                    1.88,3.19,5.24,8.37,13.06,20.19,30.00,40.00,
                    1.88,3.20,5.25,8.37,13.08,20.20,30.00,40.00,
                    1.86,3.20,5.26,8.39,13.12,20.22,30.00,40.00,
                    1.87,3.17,5.25,8.42,13.15,20.21,30.00,40.00,
                    1.86,3.18,5.28,8.45,13.18,20.24,30.00,40.00,
                    1.86,3.19,5.26,8.42,13.16,20.22,30.00,40.00,
                    1.86,3.19,5.28,8.46,13.20,20.23,30.00,40.00)
                  ,c(8,4,2),dimnames = list(Kdescr,Pdescr[1:4],Ndescr))
GPR_EI <- array(c(1.82,3.14,5.19,8.30,13.05,20.19,30.00,40.00,
                  1.84,3.16,5.22,8.33,13.07,20.20,30.00,40.00,
                  1.85,3.18,5.24,8.36,13.10,20.21,30.00,40.00,
                  1.85,3.17,5.24,8.38,13.10,20.21,30.00,40.00,
                  1.86,3.22,5.29,8.44,13.20,20.24,30.00,40.00,
                  1.88,3.24,5.30,8.46,13.18,20.21,30.00,40.00,
                  1.87,3.21,5.28,8.45,13.17,20.21,30.00,40.00,
                  1.88,3.22,5.29,8.45,13.17,20.21,30.00,40.00)
                ,c(8,4,2),dimnames = list(Kdescr,Pdescr[2:5],Ndescr))
Bayer <- array(c(1.88,3.22,5.31,8.50,13.23,20.00,30.00,40.00),c(8,1,1),dimnames = list(Kdescr))

round(res_comp[,1:4,]-GPR_Tree,3)
round(abs(res_comp[,1:4,]-GPR_Tree)*100/GPR_Tree,1)
round(res_comp[,2:5,]-GPR_EI,3)
round(abs(res_comp[,2:5,]-GPR_EI)*100/GPR_EI,1)
round(res_comp-array(Bayer,c(8,5,2)),3)
round(abs(res_comp-array(Bayer,c(8,5,2)))*100/array(Bayer,c(8,5,2)),1)
