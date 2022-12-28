### CHAPTER 4 ###

###  Includes (amongst others) the following objects:                     ###
###  The Fundamental Theorem of Derivatives Trading (the version from     ###
###  Section 4.2) in the one-dimensional case, the Black-Scholes formula, ###
###  Delta, Gamma, Implied volatility, the code for Figure 4.1 and        ###
###  Example 8 with real market data                                      ###

### install packages ###

install.packages("pracma")
library(pracma)


### The functions BM and GBM from Chapter 2 and Chapter 5 must be ###
### initialized before executing this code.                       ###


### FUNDAMENTAL THEOREM OF DERIVATIVES TRADING ###

### Case for a single risky asset ###

riskfreerate <- function(S,T = 1){
  #t <- seq(0,T,T/(length(S)-1))
  r <- rep(0.02,length(S)-1)
  return(r)
}

sigma_r <- function(Stilde,T = 1){
  #t <- seq(0,T,T/(length(Stilde)-1))
  sigma_r <- rep(0.4,length(Stilde)-1)
  return(sigma_r)
}

sigma_i <- function(S,T = 1){
  #t <- seq(0,T,T/(length(S)-1))
  sigma_i <- rep(0.5,length(S)-1)
  return(sigma_i)
}

sigma_h <- function(S,T = 1){
  #t <- seq(0,T,T/(length(S)-1))
  sigma_h <- rep(0.5,length(S)-1)
  return(sigma_h)
}

optionprocess <- function(S,K,T,r,sigma,call = TRUE){
  ifelse(call,payoff <- max(S[length(S)]-K,0),payoff <- max(K-S[length(S)],0))
  t <- seq(0,T,T/(length(S)-1))[1:(length(S)-1)]
  S <- S[1:(length(S)-1)]
  d1 <- (log(S/K) + (r + (sigma^2)/2)*(T-t))/(sigma*sqrt(T-t))
  d2 <- (log(S/K) + (r - (sigma^2)/2)*(T-t))/(sigma*sqrt(T-t))
  if (call==TRUE){
    return(c(S*pnorm(d1,0,1) - K*exp(-r*(T-t))*pnorm(d2,0,1),payoff))
  }
  else if (call==FALSE){
    return(c(K*exp(-r*(T-t))*pnorm(-d2,0,1) - S*pnorm(-d1,0,1),payoff))
  }
} 

optiongamma <- function(S,K,T,r,sigma){
  t <- seq(0,T,T/(length(S)-1))[1:(length(S)-1)]
  S <- S[1:(length(S)-1)]
  d1 <- (log(S/K) + (r + (sigma^2)/2)*(T-t))/(sigma*sqrt(T-t))
  return(c(dnorm(d1,0,1)/(S*sigma*sqrt(T-t)),0))
}

pnlprocess <- function(S,K,T,r,sigma_r,sigma_i,sigma_h,call = TRUE,whole = TRUE){
  time_grid <- seq(0,T,T/(length(S)-1))
  time_increment <- time_grid[2]
  Vi <- optionprocess(S,K,T,r,sigma_i,call)
  Vh <- optionprocess(S,K,T,r,sigma_h,call)
  Vi_increments <- diff(Vi)
  Vh_increments <- diff(Vh)
  GammaV <- optiongamma(S,K,T,r,sigma_h)
  Pi <- rep(0,length(S))
  for (t in (1:(length(S)-1))) {
    Pi[t+1] = Pi[t] + Vi_increments[t] - Vh_increments[t] + (-r[t]*(Vi[t] - Vh[t]) + 0.5*((S[t])^2)*(((sigma_r[t])^2) - ((sigma_h[t])^2))*GammaV[t])*time_increment
  }
  
  if (whole == TRUE){return(Pi)}
  else{
    Pi_increments <- diff(Pi)
    pnl <- rep(0,length(S)-1)
    for (t in (1:(length(S)-1))) {
      pnl[t] <- (exp(-sum(r[1:t]*time_increment)))*Pi_increments[t]
    }
    return(sum(pnl))
  }
}

### Example from the paper "The Fundamental Theorem of Derivative Trading - ###
### Expositions, Extensions, and Experiments" (Section 2.2)                 ###

S0 <- 100
K <- 120
T <- 1
drift <- 0.3

S <- matrix(0,10,5000)
for (i in 1:10) {
  S[i,] <- GBM(5000,T,S0,drift,0.4)
}
plot(pnlprocess(S[1,],K,T,riskfreerate(S[1,]),sigma_r(S[1,]),sigma_i(S[1,]),sigma_h(S[1,]),whole=TRUE),type="l",col=palette.colors(9)[1],ylim = c(-8,0),main = "Hedging performance when choosing the implied volatility",xlab = "time t in years",ylab = "P&L",xaxt="n")
axis(1, at=seq(0,5000,5000/5), labels=seq(0,T,T/5))
for (i in 2:10) {
  #lines(pnlprocess(S[i,],K,T,riskfreerate(S[i,]),sigma_r(S[i,]),sigma_i(S[i,]),sigma_h(S[i,]),whole=TRUE),col=palette.colors(9)[i-1])
  lines(pnlprocess(S[i,],K,T,riskfreerate(S[i,]),sigma_r(S[i,]),sigma_i(S[i,]),sigma_h(S[i,]),whole=TRUE))
  print(pnlprocess(S[i,],K,T,riskfreerate(S[i,]),sigma_r(S[i,]),sigma_i(S[i,]),sigma_h(S[i,]),whole=FALSE))
}


### The Black-Scholes formula for European options ###

BS <- function(St,K,T,t,r,sigma,call=TRUE){
  d1 <- (log(St/K) + (r + (sigma^2)/2)*(T-t))/(sigma*sqrt(T-t))
  d2 <- (log(St/K) + (r - (sigma^2)/2)*(T-t))/(sigma*sqrt(T-t))
  if (call==TRUE){
    return(St*pnorm(d1,0,1) - K*exp(-r*(T-t))*pnorm(d2,0,1))
  }
  else if (call==FALSE){
    return(K*exp(-r*(T-t))*pnorm(-d2,0,1) - St*pnorm(-d1,0,1))
  }
}

### Implied volatility ###

#at t=0

get_impliedvol <- function(V,S,K,T,r,call = TRUE){
  S <- S[1]
  V <- V[1]
  #d1 <- (log(S/K) + (r + (sigma^2)/2)*T)/(sigma*sqrt(T))
  #d2 <- (log(S/K) + (r - (sigma^2)/2)*T)/(sigma*sqrt(T))
  #if (call==TRUE){
  #  return(S*pnorm(d1,0,1) - K*exp(-r*(T-t))*pnorm(d2,0,1))
  #}
  #else if (call==FALSE){
  #  return(K*exp(-r*(T-t))*pnorm(-d2,0,1) - S*pnorm(-d1,0,1))
  #}
  f <- function(x){
    return(V - BS(S,K,T,0,r,x,call))
  }
  return(brent(f,0,1)$root)
}

#and the whole process

impliedvol <- function(V,S,K,T,r,call = TRUE){
  if (length(V) != length(S)) {stop("The length of the option's price path differs from the length of the underlying's price path.")}
  else{
    f <- function(x){return(V - optionprocess(S,K,T,r,rep(x,length(S)-1),call))}
    iv_path <- c()
    for (i in 1:(length(V)-1)) {
      g <- function(x){return(f(x)[i])}
      iv_path <- c(iv_path,brent(g,-30,30)$root) #c(iv_path,brent(g,0,1)$root)
    }
    return(iv_path)
  }
}

### Delta ###

optiondelta <- function(S,K,T,r,sigma){
  t <- seq(0,T,T/(length(S)-1))[1:(length(S)-1)]
  S_T <- S[length(S)] 
  mat <- S_T > K
  S <- S[1:(length(S)-1)]
  d1 <- (log(S/K) + (r + (sigma^2)/2)*(T-t))/(sigma*sqrt(T-t))
  if (mat == TRUE){return(c(pnorm(d1,0,1),1))}
  else {return(c(pnorm(d1,0,1),0))}
}

### EXPERIMENTS WITH REAL MARKET DATA (Example 8 in the thesis)               ###
### THE STOCK OF NETFLIX, INC. (NFLX) WAS CHOSEN AND THE DATA WAS TAKEN FROM  ###
### YAHOO FINANCE                                                             ###


### 18.09.2017 - 16.09.2022 ###

### Reading in data (that was pre-processed in the Excel sheet ### 
### "NFLX Data 5y until16.9.22")                               ###
### Further information can be found in the sheet              ###

dates <- readClipboard()
dates <- dates[2:length(dates)]

prices <- readClipboard()
prices <- as.numeric(prices[2:length(prices)])

### Plotting the stock prices over the last 5 years ###

plot(prices,type="l",ylim=c(0,750),main = "Netflix, Inc. (NFLX)",col.main="black",xlab = "Date",ylab = "Close price in USD",xaxt="n")
axis(1,at=c(1,252,502,755,1007,1259),labels = dates[c(1,252,502,755,1007,1259)],las=1)

### Implied vola ###

get_impliedvol(9.15,240.13,300,64/251,0.0467214,TRUE) #0.02263

### Calculation historic volatility ###

logreturns <- log(prices[2:length(prices)]/prices[1:(length(prices)-1)])
expected_logreturn <- mean(logreturns)
vola_over_period <- sqrt((1/(length(logreturns)-1))*sum((logreturns-expected_logreturn)^2))
period_considered <- 5/length(logreturns)
historic_vola <- (1/sqrt(period_considered))*vola_over_period
round(historic_vola*100,1)

### Calculating V_0^h ###

Vh <- BS(prices[length(prices)],300,64/251,0,0.0467214,historic_vola,call=TRUE) #0.01063 #0.02263
round(Vh,2)


### 16.09.2022 - 16.12.2022 ###

### Reading in data (that was pre-processed in the Excel sheet ### 
### "NFLX Data 16.9.22 - 16.12.22")                            ###
### Further information can be found in the sheet              ###

dates3M <- readClipboard()
dates3M <- dates3M[2:length(dates3M)]

prices3M <- readClipboard()
prices3M <- as.numeric(prices3M[2:length(prices3M)])

### Reading in data (that was pre-processed in the Excel sheet ### 
### "Netflix Option Prices")                                   ###
### Further information can be found in the sheet              ###

optionprices <- readClipboard()
optionprices <- as.numeric(optionprices[2:length(optionprices)])
optionprices <- rev(optionprices)

### Plotting the stock prices ###

plot(prices3M,type="l",ylim=c(0,350),main = "The spot price of Netflix, Inc. (NFLX) \nover the lifetime of the option",col.main="black",xlab = "Date",ylab = "Close price in USD",xaxt="n")
lines(rep(300,length(prices3M)),col="red")
axis(1,at=c(1,22,44,65),labels = dates3M[c(1,22,44,65)],las=1)

### Calculation realized volatility ###

logreturns3M <- log(prices3M[2:length(prices3M)]/prices3M[1:(length(prices3M)-1)])
expected_logreturn3M <- mean(logreturns3M)
vola_over_period3M <- sqrt((1/(length(logreturns3M)-1))*sum((logreturns3M-expected_logreturn3M)^2))
period_considered3M <- 0.25/length(logreturns3M)
realized_vola <- (1/sqrt(period_considered3M))*vola_over_period3M
round(realized_vola*100,1)

plot(c(prices,prices3M[2:length(prices3M)]),type="l",ylim=c(0,750),main = "Netflix, Inc. (NFLX)",col.main="black",xlab = "Date",ylab = "Close price in USD",xaxt="n")
#lines(rep(300,length(c(prices,prices3M[2:length(prices3M)]))),col="red")
axis(1,at=c(1,252,502,755,1007,1259,1323),labels = c(substr(dates,3,10),substr(dates3M[2:length(dates3M)],3,10))[c(1,252,502,755,1007,1259,1323)],las=2)
#axis(1,at=c(1,252,502,755,1007,1259,1323),labels = FALSE)
#text(x=c(1,252,502,755,1007,1259,1323),y=rep(-120,7),xpd=NA,labels = c(dates,dates3M[2:length(dates3M)])[c(1,252,502,755,1007,1259,1323)],srt=55)

### Calculating hedging performance ###

### Reading in data (that was pre-processed in the Excel sheet ### 
### "Netflix Option Prices")                                   ###
### Further information can be found in the sheet              ###

riskfreerate <- readClipboard()
riskfreerate <- as.numeric(riskfreerate[2:(length(riskfreerate))])
riskfreerate <- (1/100)*rev(riskfreerate)
riskfreerate <- riskfreerate[1:(length(riskfreerate)-1)]

sigma_r <- rep(realized_vola,length(prices3M)-1)

sigma_i <- impliedvol(optionprices,prices3M,300,64/251,riskfreerate,TRUE)

sigma_h <- rep(historic_vola,length(prices3M)-1)

pnlprocess(prices3M,300,64/251,riskfreerate,sigma_r,sigma_i,sigma_h,call = TRUE,whole = FALSE)

plot(pnlprocess(prices3M,300,64/251,riskfreerate,sigma_r,sigma_i,sigma_h,call = TRUE,whole = TRUE),type = "l",ylim=c(-10,11),main = expression(paste("Value process of the hedge portfolio ",Pi)),col.main="black",xlab = "Date",ylab = "Value in USD",xaxt="n")
axis(1,at=c(1,22,44,65),labels = dates3M[c(1,22,44,65)],las=1)
plot(optionprocess(prices3M,300,64/251,riskfreerate,sigma_r,call = TRUE),type = "l",ylim=c(-1,30),main = "Value process of the option according to B-S \nversus the quoted option prices at NASDAQ",col.main="black",col="blue",xlab = "Date",ylab = "Value in USD",xaxt="n")
lines(optionprices,col="red")
legend(1,30,legend = c("B-S","NASDAQ"),fill = c("blue","red"))
axis(1,at=c(1,22,44,65),labels = dates3M[c(1,22,44,65)],las=1)
#plot(optiongamma(prices3M,300,63/241,riskfreerate,sigma_r),type = "l",ylim=c(0,0.030),main = "Gamma of the option",col.main="black",xlab = "Date",ylab = expression(Gamma),xaxt="n")
#axis(1,at=c(1,22,44,65),labels = dates3M[c(1,22,44,65)],las=1)
par(mar=c(5, 4, 4, 6) + 0.1)
plot(optiondelta(prices3M,300,64/251,riskfreerate,sigma_r),type = "l",ylim=c(-0.1,1),col="green",main = "Delta and Gamma of the option",col.main="black",xlab = "",ylab = "",xaxt = "n",yaxt = "n")
axis(1,at=c(1,22,44,65),labels = dates3M[c(1,22,44,65)],las=1)
mtext("Date",side = 1, line = 3)
axis(2,ylim=c(-0.1,1),las=1,col.axis="green",col="green")
mtext(expression(Delta),side = 2, line = 3,col = "green")
par(new=TRUE)
plot(optiongamma(prices3M,300,64/251,riskfreerate,sigma_r),type = "l",ylim=c(0,0.030),col="purple",xlab = "",ylab = "",axes=FALSE)
axis(4,ylim=c(0,0.030),las=1,col.axis="purple",col="purple")
mtext(expression(Gamma),side=4,line=4,col="purple")
legend(2,0.030,legend = c("Delta","Gamma"),fill = c("green","purple"))
#plot(optionprices,type="l",ylim=c(0,25),main = "Option prices at NASDAQ",col.main="black",xlab = "Date",ylab = "Close price in USD",xaxt="n")
#axis(1,at=c(1,22,44,65),labels = dates3M[c(1,22,44,65)],las=1)
plot(riskfreerate*100,type="l",ylim=c(4,6),main = "USD Libor 12 months",col.main="black",xlab = "Date",ylab = "Rate in %",xaxt="n")
axis(1,at=c(1,22,44,65),labels = dates3M[c(1,22,44,65)],las=1)
plot(sigma_i*100,type="l",ylim=c(-1,65),main = "Implied volatility of the option",col.main="black",xlab = "Date",ylab = "Implied volatility in %",xaxt="n")
axis(1,at=c(1,22,44,65),labels = dates3M[c(1,22,44,65)],las=1)
