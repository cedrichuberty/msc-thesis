### SIGNATURES, SHUFFLE PRODUCT, APPROACHES IN "OPTIMAL EXECUTION WITH ROUGH PATH ###
### SIGNATURES" AND "OPTIMAL STOPPING WITH SIGNATURES"                            ###

### install packages ###

install.packages("BBmisc")
install.packages("lpSolve")
library(BBmisc)
library(lpSolve)

### Signature ### (taken from Github iisignature)

#Takes a path in R^d given by an array of n points with dims c(d,n)
#and returns its signature up to level m excluding level 0
#e.g. path=array(c(1,1,3,10,20,30),c(3,2))
#Note that the order of path's dimensions is the other 
#way around from iisignature.sig.
sig=function(path,m, flat=TRUE){
  n=dim(path)[2]
  d=dim(path)[1]
  diffs=path[,-1,drop=FALSE]-path[,-n,drop=FALSE]
  if (n<2){
    o=sapply(1:m,function(x)rep(0,d**x))
    if (flat) return (do.call("c",o))
    return (o)
  }
  r=lapply(1:(n-1),function(x)Reduce(kronecker,rep(list(diffs[,x]),m),accumulate=TRUE))
  facts=lapply(1:m,factorial)
  r=lapply(r,function(x)mapply("/",x,facts))
  chen=function(x,y) c(list(x[[1]]+y[[1]]),
                       lapply(2:m,
                              function(z)x[[z]]+y[[z]]
                              +Reduce("+",mapply(kronecker,x[1:z-1],rev(y[1:z-1]),
                                                 SIMPLIFY=FALSE))))
  o=Reduce(chen,r)
  if (flat) return (do.call("c",o))
  o
}

### Brownian motion ###

BM = function(T,d,n){
  W <- rep(0,n*d)
  delta <- T/n
  for (t in seq(d+1,n*d+1,d)) {
    for (i in 0:(d-1)) {
      W[t+i] <- W[t+i-d] + sqrt(delta)*rnorm(1,0,1)  
    }
  }
  return(W)
}

d <- 3 #dimension
N <- 5 #order of truncation
T <- 1 #time horizon
n <- 1000

sigofBM <- c(1,sig(array(BM(T,d,n),c(d,n+1)),N))
  
### Expected signature ###  
  
exp_sig <- function(sig){#also works for paths - not only for signatures
  expsig <- rep(0,dim(sig)[2])
  for (j in (1:dim(sig)[2])) {
    for (i in (1:dim(sig)[3])) {
      expsig[j] <- expsig[j] + sig[,j,i]
    }
    expsig[j] <- expsig[j]/dim(sig)[3]
  }
  return(expsig)
}  


### Shuffle product ###

shuffle <- function(a,b){
  if(nchar(a)==0){return(b)}
  else if (nchar(b)==0){return(a)} 
  else {return(rev(c(paste(shuffle(substr(a,1,nchar(a)-1),b),substr(a,nchar(a),nchar(a)),sep=""),paste(shuffle(a,substr(b,1,nchar(b)-1)),substr(b,nchar(b),nchar(b)),sep=""))))}
}

simplify(paste(1,shuffle("1","2")))

shuffle2 <- function(a,b,simple = FALSE){
  t <- c()
  for (i in 1:length(a)) {
    for (j in 1:length(b)) {
      coeffa <- as.numeric(strsplit(a[i]," ")[[1]][1])
      coeffb <- as.numeric(strsplit(b[j]," ")[[1]][1])
      worda <- ifelse(is.na(strsplit(a[i]," ")[[1]][2]),"",strsplit(a[i]," ")[[1]][2])
      wordb <- ifelse(is.na(strsplit(b[j]," ")[[1]][2]),"",strsplit(b[j]," ")[[1]][2])
      t <- c(t,paste(coeffa*coeffb,shuffle(worda,wordb)))
    }
  }
  ifelse(simple,return(simplify(t)),return(t))
}

shuffle3 <- function(a,b,simple = FALSE){
  t <- c()
  for (i in 1:length(a)) {
    for (j in 1:length(b)) {
      coeffa <- strsplit(a[i]," ")[[1]][1]
      coeffb <- strsplit(b[j]," ")[[1]][1]
      worda <- ifelse(is.na(strsplit(a[i]," ")[[1]][2]),"",strsplit(a[i]," ")[[1]][2])
      wordb <- ifelse(is.na(strsplit(b[j]," ")[[1]][2]),"",strsplit(b[j]," ")[[1]][2])
      t <- c(t,paste(paste("(",coeffa,")","*","(",coeffb,")",sep=""),shuffle(worda,wordb)))
    }
  }
  ifelse(simple,return(simplify(t)),return(t))
}

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


### The scalar product for words and signatures ###

scalarproduct <- function(a,b,keys){
  aux <- rep("",length(b))
  for (i in (1:length(b))) {
    aux[i] <- paste(b[i],keys[i])
  }
  #find terms a and aux with same keys and add them
  sum <- 0
  for (j in 1:length(a)) {
    for (k in 1:length(aux)) {
      sum = sum + as.numeric(strsplit(multiplication(a[j],aux[k])," ")[[1]][1]) 
    }
  }
  return(sum)
}

### LHS in the scalar product in Proposition 4.2 in ###
### "Optimal execution with rough path signatures"  ###

#execution <- function(l,lambda=0.001,q0=1,alpha,phi){
#  return(c(paste(1,paste(shuffle("2",l),1,sep="")),
#           paste(1,paste(shuffle("",l),1,sep="")),
#           paste(-lambda,paste(shuffle(l,l),1,sep="")),
#           paste(-(q0^2)*phi,(shuffle("","1"))),
#           paste(2*q0*phi,(paste(l,shuffle("","11"),sep=""))),
#           paste(-phi,paste(shuffle(paste(l,1,sep=""),paste(l,1,sep="")),1,sep="")),
#           paste(-(q0^2)*alpha,shuffle("","")),
#           paste(2*q0*alpha,shuffle("",paste(l,1,sep=""))),
#           paste(-alpha,shuffle(paste(l,1,sep=""),paste(l,1,sep=""))),
#           paste(q0,shuffle("2","")),
#           paste(q0,shuffle("","")),
#           paste(-q0*lambda,shuffle("",l)),
#           paste(-1,shuffle(paste(l,1,sep=""),"2")),
#           paste(-1,paste(l,1,sep="")),
#           paste(lambda,shuffle(paste(l,1,sep=""),l)))
#  )
#}

execution <- function(a,lambda=0.001,q0=1,alpha,phi){
  return(c(paste(q0-(q0^2)*alpha,""),
           paste(shuffle2(paste(1,"2"),a),"1",sep=""),
           concatenation(shuffle2(a,a),paste(-lambda,"1")), 
           paste(-(q0^2)*phi,"1"),
           concatenation(paste(a,"11",sep=""),paste(2*q0*phi,"")),
           concatenation(shuffle2(paste(a,"1",sep=""),paste(a,"1",sep="")),paste(-phi,"1")),
           concatenation(paste(a,"1",sep=""),paste(2*q0*alpha,"")),
           concatenation(shuffle2(paste(a,"1",sep=""),paste(a,"1",sep="")),paste(-alpha,"")),
           paste(q0,"2"),
           concatenation(a,paste(-q0*lambda,"")),
           shuffle2(paste(a,"1",sep=""),paste(-1,"2")),
           concatenation(shuffle2(paste(a,1,sep=""),a),paste(lambda,"")))
  )
}

execution2 <- function(x){
  words <- keys(2,truncorder)
  a <- c()
  for (w in (1:length(words))) {
    a <- c(a,paste(x[w],words[w]))
  }
  execution(a,lambda,q0,alpha,phi)
}

as.numeric(y)*3
quote(3 * 2)

words <- keys(2,truncorder)
a <- c()
x <- c()
for (w in (1:length(words))) {
  x <- c(x,paste("x_",w-1,sep=""))
  a <- c(a,paste(paste("x_",w-1,sep=""),words[w]))
}

execution3 <- function(a,lambda=0.001,q0=1,alpha,phi){ 
  return(c(paste(q0-(q0^2)*alpha,""),
           paste(shuffle3(paste(1,"2"),a),"1",sep=""),
           concatenation3(shuffle3(a,a),paste(-lambda,"1")), 
           paste(-(q0^2)*phi,"1"),
           concatenation3(paste(a,"11",sep=""),paste(2*q0*phi,"")),
           concatenation3(shuffle3(paste(a,"1",sep=""),paste(a,"1",sep="")),paste(-phi,"1")),
           concatenation3(paste(a,"1",sep=""),paste(2*q0*alpha,"")),
           concatenation3(shuffle3(paste(a,"1",sep=""),paste(a,"1",sep="")),paste(-alpha,"")),
           paste(q0,"2"),
           concatenation3(a,paste(-q0*lambda,"")),
           shuffle3(paste(a,"1",sep=""),paste(-1,"2")),
           concatenation3(shuffle3(paste(a,1,sep=""),a),paste(lambda,"")))
  )
} #if strsplit...[2] > truncorder, delete, don't keep

### Looking for the supremum from Proposition 4.2 in ###
### "Optimal execution with rough path signatures"   ###
###                                                  ###
### very naive approach for now                      ###

#optimalexecution <- function(pricepath,truncorder,lambda = 0.001,q0 = 1,alpha,phi){
#  dimpath <- dim(pricepath)[1]
#  lenpath <- dim(pricepath)[2]
#  numpath <- dim(pricepath)[3]
#  signa <- array(0,c(1,((dimpath*((dimpath^truncorder)-1))/(dimpath-1))+1,numpath))
#  for (i in (1:numpath)) {
#    signa[,,i] <- c(1,sig(array(pricepath[,,i],c(dimpath,lenpath)),truncorder))
#  }
#  words <- keys(2,truncorder)
#  keys <- keys(dimpath,truncorder)
#  expsigna <- exp_sig(signa)
#  supremum <- scalarproduct(execution(words[1],lambda,q0,alpha,phi),expsigna,keys)
#  l <- words[1]
  
#  pb <- txtProgressBar(min = 0, max = length(words), style = 3, width = 50, char = "=")
  
#  for (i in (2:length(words))) {
#    candidate <- scalarproduct(execution(words[i],lambda,q0,alpha,phi),expsigna,keys)
#    if (candidate >= supremum){
#      supremum <- candidate
#      l <- words[i]
#    }
#    setTxtProgressBar(pb, i)
#  }
#  close(pb)
  
#  return(c(supremum,l))
#}

optimalexecution <- function(pricepath,truncorder,lambda = 0.001,q0 = 1,alpha,phi){
  dimpath <- dim(pricepath)[1]
  lenpath <- dim(pricepath)[2]
  numpath <- dim(pricepath)[3]
  signa <- array(0,c(1,((dimpath*((dimpath^truncorder)-1))/(dimpath-1))+1,numpath))
  pb <- txtProgressBar(min = 1, max = numpath, style = 3, width = 50, char = "=")
  for (i in (1:numpath)) {
    signa[,,i] <- c(1,sig(array(pricepath[,,i],c(dimpath,lenpath)),truncorder))
    setTxtProgressBar(pb, i)
  }
  close(pb)
  expsigna <- exp_sig(signa)
  key <- keys(dimpath,truncorder)
  
  words <- keys(2,truncorder)
  #a <- c()
  #x <- c()
  #for (w in (1:length(words))) {
  #  x <- c(x,paste("x_",w-1,sep=""))
  #  a <- c(a,paste(paste("x_",w-1,sep=""),words[w]))
  #}
  
  costfunction <- function(x){
    a <- c()
    for (w in (1:length(words))) {
      a <- c(a,paste(x[w],words[w]))
    }
    return(scalarproduct(execution(a,lambda,q0,alpha,phi),expsigna,key))
  }
  
  #scalarproduct(execution(a,lambda,q0,alpha,phi),expsigna,key)
  
  optimum <- optim(rep(0,length(words)),costfunction,control = list(fnscale = -1))
  
  return(c(paste("Cost:",optimum$value,"Coefficients:"),optimum$par))
}

truncorder <- 1
#execute everything above costfunction 
a <- c()
for (w in (1:length(words))) {
  a <- c(a,paste(1,words[w]))
}
scalarproduct(execution(a,lambda,q0,alpha,phi),expsigna,key)

### Example in the thesis ###

T <-1/250
n <-8*60
paths_to_be_simulated <- 10000
augmentedstock <- array(0,c(2,n,paths_to_be_simulated))
pb <- txtProgressBar(min = 1, max = paths_to_be_simulated, style = 3, width = 50, char = "=")
for (i in 1:paths_to_be_simulated) {
  augmentedstock[1,,i] <- seq(0,T,T/(n-1))
  augmentedstock[2,,i] <- 0.2*BM(T,1,n-1)
  setTxtProgressBar(pb, i)
}
close(pb)

truncorder <- 2
lambda <- 0.001
q0 <- 0
alpha <- 100000
phi <- 0.0001

optimalexecution(augmentedstock,truncorder,lambda,q0,alpha,phi)

sol <- c(paste(-4.63154540495257*10^(-8),""),paste(-0.0515163150843151,"1"),paste(-0.0154001681435968,"2"),paste(0.0126901683658508,"11"),paste(0.034538207873199,"12"),paste(0.0515311274900253,"21"),paste(0.0742716522713664,"22"))


### Figure in the thesis ###

averageprice <- rep(0,dim(augmentedstock)[2])
for (j in (1:dim(augmentedstock)[2])) {
  for (i in (1:dim(augmentedstock)[3])) {
    averageprice[j] <- averageprice[j] + augmentedstock[2,j,i]
  }
  averageprice[j] <- averageprice[j]/dim(augmentedstock)[3]
}
averageaugmentedprice <- array(0,c(2,n))
averageaugmentedprice[1,] <- augmentedstock[1,,1]
averageaugmentedprice[2,] <- averageprice

tradstrataverage <- c()
costaverage <- c()
for (i in 1:dim(averageaugmentedprice)[2]) {
  s <- c(1,sig(array(averageaugmentedprice[,1:i],c(2,i)),truncorder))
  tradstrataverage <- c(tradstrataverage,scalarproduct(sol,s,keys(2,truncorder)))
  costaverage <- c(costaverage,scalarproduct(execution(sol,lambda,q0,alpha,phi),s,keys(2,truncorder)))
}
wealthaverage <- (averageaugmentedprice[2,] - lambda*tradstrataverage)*tradstrataverage
wealthaverage <- cumsum(wealthaverage)
inventoryaverage <- -cumsum(tradstrataverage)

#plot(augmentedstock[2,,1],type = "l",xaxt="n")
#axis(1, at = seq(0,n,n/8), labels = c("9:00","10:00","11:00","12:00","13:00","14:00","15:00","16:00","17:00"))


chosenpaths <- round(runif(5,1,10000)) 
j <- chosenpaths[1]
tradstrat <- c()
cost <- c()
for (i in 1:dim(augmentedstock)[2]) {
  s <- c(1,sig(array(augmentedstock[,1:i,j],c(2,i)),truncorder))
  tradstrat <- c(tradstrat,scalarproduct(sol,s,keys(2,truncorder)))
  cost <- c(cost,scalarproduct(execution(sol,lambda,q0,alpha,phi),s,keys(2,truncorder)))
}
wealth <- (augmentedstock[2,,j] - lambda*tradstrat)*tradstrat
wealth <- cumsum(wealth)
inventory <- -cumsum(tradstrat)
#plot(tradstrat,type = "l",col="blue",xaxt="n")
#plot(wealth,type = "l",col="blue",xaxt="n")
plot(cost,type="l",col="blue",xaxt="n")
#plot(inventory,type = "l",col="blue",xaxt="n",ylim = c(-0.04,0.1))
for (j in chosenpaths[2:length(chosenpaths)]) {
  #tradstrat <- c()
  cost <- c()
  for (i in 1:dim(augmentedstock)[2]) {
    s <- c(1,sig(array(augmentedstock[,1:i,j],c(2,i)),truncorder))
    #tradstrat <- c(tradstrat,scalarproduct(sol,s,keys(2,truncorder)))
    cost <- c(cost,scalarproduct(execution(sol,lambda,q0,alpha,phi),s,keys(2,truncorder)))
  }
  #wealth <- (augmentedstock[2,,j] - lambda*tradstrat)*tradstrat
  #wealth <- cumsum(wealth)
  #inventory <- -cumsum(tradstrat)
  #lines(tradstrat,col="blue")
  #lines(wealth,col="blue")
  lines(cost,col="blue")
  #lines(inventory,col="blue")
}
#lines(tradstrataverage,col="red")
#lines(wealthaverage,col="red")
lines(costaverage,col="red")
#lines(inventoryaverage,col="red")
axis(1, at = seq(0,n,n/8), labels = c("9:00","10:00","11:00","12:00","13:00","14:00","15:00","16:00","17:00"))

plot(augmentedstock[2,,j],type="l")
par(new=TRUE)
plot(tradstrat,col="red",type="l")
par(new=TRUE)
plot(wealth,col="blue",type="l")
par(new=TRUE)
plot(cost,col="purple",type="l")
par(new=TRUE)
plot(inventory,col="green",type = "l")

### Test ###

T<-1
d<-2
n<-100
m<-1000
c <- c()
b <- c()
for (i in (1:m)) {
  c <- c(c,BM(T,d,n))
  b <- c(b,BM(T,d,n))
}
X <- array(c,c(d,n+1,m))

optimalexecution(X,4,0.001,1,0.1,0.01)

### LHS in the scalar product in Corollary 6.6 in ###
### "Optimal stopping with signatures"            ###

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




optimalstopping <- function(pricepath,truncorder){
  dimpath <- dim(pricepath)[1]
  lenpath <- dim(pricepath)[2]
  numpath <- dim(pricepath)[3]
  ifelse(dimpath==1,signa <- array(0,c(1,1+truncorder,numpath)),signa <- array(0,c(1,((dimpath*((dimpath^truncorder)-1))/(dimpath-1))+1,numpath)))
  for (i in (1:numpath)) {
    signa[,,i] <- c(1,sig(array(pricepath[,,i],c(dimpath,lenpath)),truncorder))
  }
  expsigna <- exp_sig(signa)
  
  #optimizing
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
        scalarproduct(stopping(word,truncorder),expsigna,keys(3,N))
      } #or something like use the distributivity of scalar product and take out coefficents so that you end up with only the coefficients as unknowns and then use Linear Programming for each l(=combination of keys without coefficients)
    }
  }
  scalarproduct(stopping(l,truncorder),expsigna,key)
  return()
}

### Get the words corresponding to the terms in the signature ###

keys <- function(d,truncorder){ #applied brute force for now
  keys <- c("")
  for (a in 1:d) {
    keys <- c(keys,a)
  }
  if (truncorder >= 2){
    for (a in 1:d) {
      for (b in 1:d) {
        comb <- paste(a,b,sep="")
        keys <- c(keys,comb)
      }
    }
  }
  if (truncorder >= 3){
    for (a in 1:d) {
      for (b in 1:d) {
        for (c in 1:d) {
          comb <- paste(paste(a,b,sep=""),c,sep="")
          keys <- c(keys,comb)   
        }
      }
    }
  }
  if (truncorder >= 4){
    for (a in 1:d) {
      for (b in 1:d) {
        for (c in 1:d) {
          for (e in 1:d) {
            comb <- paste(paste(paste(a,b,sep=""),c,sep=""),e,sep="")
            keys <- c(keys,comb)   
          }
        }
      }
    }
  }
  if (truncorder >= 5){
    for (a in 1:d) {
      for (b in 1:d) {
        for (c in 1:d) {
          for (e in 1:d) {
            for (f in 1:d) {
              comb <- paste(paste(paste(paste(a,b,sep=""),c,sep=""),e,sep=""),f,sep="")
              keys <- c(keys,comb)
            }
          }
        }
      }
    }
  }
  if (truncorder >= 6){
    for (a in 1:d) {
      for (b in 1:d) {
        for (c in 1:d) {
          for (e in 1:d) {
            for (f in 1:d) {
              for (g in 1:d) {
                comb <- paste(paste(paste(paste(paste(a,b,sep=""),c,sep=""),e,sep=""),f,sep=""),g,sep="")
                keys <- c(keys,comb)
              }
            }
          }
        }
      }
    }
  }
  if (truncorder >= 7){
    for (a in 1:d) {
      for (b in 1:d) {
        for (c in 1:d) {
          for (e in 1:d) {
            for (f in 1:d) {
              for (g in 1:d) {
                for (h in 1:d) {
                  comb <- paste(paste(paste(paste(paste(paste(a,b,sep=""),c,sep=""),e,sep=""),f,sep=""),g,sep=""),h,sep="")
                  keys <- c(keys,comb)
                }
              }
            }
          }
        }
      }
    }
  }
  if (truncorder >= 8){
    for (a in 1:d) {
      for (b in 1:d) {
        for (c in 1:d) {
          for (e in 1:d) {
            for (f in 1:d) {
              for (g in 1:d) {
                for (h in 1:d) {
                  for (i in 1:d) {
                    comb <- paste(paste(paste(paste(paste(paste(paste(a,b,sep=""),c,sep=""),e,sep=""),f,sep=""),g,sep=""),h,sep=""),i,sep="")
                    keys <- c(keys,comb)
                  }
                }
              }
            }
          }
        }
      }
    }
  }
  if (truncorder >= 9){
    for (a in 1:d) {
      for (b in 1:d) {
        for (c in 1:d) {
          for (e in 1:d) {
            for (f in 1:d) {
              for (g in 1:d) {
                for (h in 1:d) {
                  for (i in 1:d) {
                    for (j in 1:d) {
                      comb <- paste(paste(paste(paste(paste(paste(paste(paste(a,b,sep=""),c,sep=""),e,sep=""),f,sep=""),g,sep=""),h,sep=""),i,sep=""),j,sep="")
                      keys <- c(keys,comb)
                    }
                  }
                }
              }
            }
          }
        }
      }
    }
  }
  if (truncorder >= 10){
    for (a in 1:d) {
      for (b in 1:d) {
        for (c in 1:d) {
          for (e in 1:d) {
            for (f in 1:d) {
              for (g in 1:d) {
                for (h in 1:d) {
                  for (i in 1:d) {
                    for (j in 1:d) {
                      for (k in 1:d) {
                        comb <- paste(paste(paste(paste(paste(paste(paste(paste(paste(a,b,sep=""),c,sep=""),e,sep=""),f,sep=""),g,sep=""),h,sep=""),i,sep=""),j,sep=""),k,sep="")
                        keys <- c(keys,comb)
                      }
                    }
                  }
                }
              }
            }
          }
        }
      }
    }
  }
  if (truncorder >= 11){stop("Order of truncation to high.")}
  return(keys)
}

### Auxiliary functions ###

#concatenation <- function(a,b){
#  return(paste(as.numeric(strsplit(a," ")[[1]][1])*as.numeric(strsplit(b," ")[[1]][1]),gsub("NA","",paste(strsplit(a," ")[[1]][2],strsplit(b," ")[[1]][2],sep=""))))
#}

concatenation <- function(a,b){
  concat <- c()
  for (w in (1:length(a))) {
    for (v in (1:length(b))) {
      concat <- c(concat,paste(as.numeric(strsplit(a[w]," ")[[1]][1])*as.numeric(strsplit(b[v]," ")[[1]][1]),gsub("NA","",paste(strsplit(a[w]," ")[[1]][2],strsplit(b[v]," ")[[1]][2],sep=""))))
    }
  }
  return(concat)
}

concatenation3 <- function(a,b){
  concat <- c()
  for (w in (1:length(a))) {
    for (v in (1:length(b))) {
      concat <- c(concat,paste(paste("(",strsplit(a[w]," ")[[1]][1],")","*","(",strsplit(b[v]," ")[[1]][1],")",sep=""),gsub("NA","",paste(strsplit(a[w]," ")[[1]][2],strsplit(b[v]," ")[[1]][2],sep=""))))
    }
  }
  return(concat)
}


multiplication <- function(a,b){
  if ((is.na(strsplit(a," ")[[1]][2]))){
    if ((is.na(strsplit(b," ")[[1]][2]))){
      return(paste(as.numeric(strsplit(a," ")[[1]][1]) * as.numeric(strsplit(b," ")[[1]][1]),gsub("NA","",paste(strsplit(a," ")[[1]][2],"",sep=""))))
    }
    else {return("0")}#stop("Not the same word.")}
  }
  else {
    if ((is.na(strsplit(b," ")[[1]][2]))){return("0")
      #stop("Not the same word.")
    }
    else if (strsplit(a," ")[[1]][2]==strsplit(b," ")[[1]][2]){
      return(paste(as.numeric(strsplit(a," ")[[1]][1]) * as.numeric(strsplit(b," ")[[1]][1]),gsub("NA","",paste(strsplit(a," ")[[1]][2],"",sep=""))))
    }
    else {return("0")}#stop("Not the same word.")}
  }
}


summation <- function(a,b){
  if ((is.na(strsplit(a," ")[[1]][2]))){
    if ((is.na(strsplit(b," ")[[1]][2]))){
      return(paste(as.numeric(strsplit(a," ")[[1]][1]) + as.numeric(strsplit(b," ")[[1]][1]),gsub("NA","",paste(strsplit(a," ")[[1]][2],"",sep=""))))
    }
    #else {return("0")}#stop("Not the same word.")}
  }
  else {
    #if ((is.na(strsplit(b," ")[[1]][2]))){return("0")
      #stop("Not the same word.")
    #}
    #else 
    if ((is.na(strsplit(b," ")[[1]][2]))==FALSE){
      if (strsplit(a," ")[[1]][2]==strsplit(b," ")[[1]][2]){
        return(paste(as.numeric(strsplit(a," ")[[1]][1]) + as.numeric(strsplit(b," ")[[1]][1]),gsub("NA","",paste(strsplit(a," ")[[1]][2],"",sep=""))))
      }
    }
    #else {return("0")}#stop("Not the same word.")}
  }
}

summation3 <- function(a,b){
  if ((is.na(strsplit(a," ")[[1]][2]))){
    if ((is.na(strsplit(b," ")[[1]][2]))){
      return(paste(paste("((",strsplit(a," ")[[1]][1],")","+","(",strsplit(b," ")[[1]][1],"))",sep=""),gsub("NA","",paste(strsplit(a," ")[[1]][2],"",sep=""))))
    }
    #else {return("0")}#stop("Not the same word.")}
  }
  else {
    #if ((is.na(strsplit(b," ")[[1]][2]))){return("0")
    #stop("Not the same word.")
    #}
    #else 
    if ((is.na(strsplit(b," ")[[1]][2]))==FALSE){
      if (strsplit(a," ")[[1]][2]==strsplit(b," ")[[1]][2]){
        return(paste(paste("((",strsplit(a," ")[[1]][1],")","+","(",strsplit(b," ")[[1]][1],"))",sep=""),gsub("NA","",paste(strsplit(a," ")[[1]][2],"",sep=""))))
      }
    }
    #else {return("0")}#stop("Not the same word.")}
  }
}

degree_of_word <- function(a){
  deg <- 0 
  for (w in 1:length(a)) {
    if (strsplit(a[w]," ")[[1]][1] != 0){
      if (!(is.na(strsplit(a[w]," ")[[1]][2]))){
        if (deg < nchar(strsplit(a[w]," ")[[1]][2])){deg <- nchar(strsplit(a[w]," ")[[1]][2])}
      }
    }
  }
  return(as.numeric(deg))
}

l1_words <- function(a){
  l1 <- 0
  for (w in 1:length(a)) {
    l1 <- l1 + abs(as.numeric(strsplit(a[w]," ")[[1]][1]))
  }
  return(l1)
}

simplify <- function(a){
  maxt <- 1 
  maxd <- 1
  for (w in 1:length(a)) {
    if (!(is.na(strsplit(a[w]," ")[[1]][2]))){
      if (maxt < nchar(strsplit(a[w]," ")[[1]][2])){maxt <- nchar(strsplit(a[w]," ")[[1]][2])}
      for (v in (1:(nchar(strsplit(a[w]," ")[[1]][2])))) {
        if (maxd < substr(strsplit(a[w]," ")[[1]][2],v,v)){maxd <- substr(strsplit(a[w]," ")[[1]][2],v,v)}
      }
    }
  }
  key <- keys(as.numeric(maxd),as.numeric(maxt))
  aux <- rep("",length(key))
  for (i in (1:length(key))) {
    aux[i] <- paste(0,key[i])
  }
  res <- c()
  for (j in 1:length(aux)) {
    sum <- aux[j]
    for (k in 1:length(a)) { 
      if (is.null(summation(sum,a[k]))==FALSE){
        sum <- summation(sum,a[k])
      }
    }
    if (strsplit(sum," ")[[1]][1] != "0"){
      res <- c(res,sum)
    }
  }
  return(res)
}

simplify3 <- function(a){
  maxt <- 1 
  maxd <- 1
  for (w in 1:length(a)) {
    if (!(is.na(strsplit(a[w]," ")[[1]][2]))){
      if (maxt < nchar(strsplit(a[w]," ")[[1]][2])){maxt <- nchar(strsplit(a[w]," ")[[1]][2])}
      for (v in (1:(nchar(strsplit(a[w]," ")[[1]][2])))) {
        if (maxd < substr(strsplit(a[w]," ")[[1]][2],v,v)){maxd <- substr(strsplit(a[w]," ")[[1]][2],v,v)}
      }
    }
  }
  key <- keys(as.numeric(maxd),as.numeric(maxt))
  aux <- rep("",length(key))
  for (i in (1:length(key))) {
    aux[i] <- paste(0,key[i])
  }
  res <- c()
  for (j in 1:length(aux)) {
    sum <- aux[j]
    for (k in 1:length(a)) { 
      if (is.null(summation3(sum,a[k]))==FALSE){
        sum <- summation3(sum,a[k])
      }
    }
    if (strsplit(sum," ")[[1]][1] != "0"){
      res <- c(res,sum)
    }
  }
  return(res)
}
