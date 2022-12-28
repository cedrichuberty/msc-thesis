### CHAPTER 2 ###

###  Includes (amongst others) the following objects:                   ###
###  Pairing, degree, shuffle product, concatenation product, signature ###
###  Brownian motion and the experiment in Section 2.6                  ###


### First running the auxiliary functions at the end is best ###
### They are necessary to work on the space of words         ###


### The pairing for words and signatures ###

pairing <- function(a,b,keys){
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

### The degree of a word ###

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

### Shuffle product ###

shuffle <- function(a,b){
  if(nchar(a)==0){return(b)}
  else if (nchar(b)==0){return(a)} 
  else {return(rev(c(paste(shuffle(substr(a,1,nchar(a)-1),b),substr(a,nchar(a),nchar(a)),sep=""),paste(shuffle(a,substr(b,1,nchar(b)-1)),substr(b,nchar(b),nchar(b)),sep=""))))}
}

shuffle("12","3")
shuffle("12","23")

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

shuffle2(c(paste(1,"1"),paste(-1,"2")),paste(3,"2"))

### Concatenation product ###

concatenation <- function(a,b){
  concat <- c()
  for (w in (1:length(a))) {
    for (v in (1:length(b))) {
      concat <- c(concat,paste(as.numeric(strsplit(a[w]," ")[[1]][1])*as.numeric(strsplit(b[v]," ")[[1]][1]),gsub("NA","",paste(strsplit(a[w]," ")[[1]][2],strsplit(b[v]," ")[[1]][2],sep=""))))
    }
  }
  return(concat)
}

### Signature ### (taken from Github iisignature)

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

# Example #

d <- 3 #dimension
N <- 5 #order of truncation
T <- 1 #time horizon
n <- 1000

sigofBM <- c(1,sig(array(BM(T,d,n),c(d,n+1)),N))


### The experiment for the signature method applied to the optimal execution ###
### problem in Section 2.6                                                   ###

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

### LHS in the pairing in Proposition 4.2 in        ###
### "Optimal execution with rough path signatures"  ###

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

### Looking for the supremum from Proposition 4.2 in ###
### "Optimal execution with rough path signatures"   ###

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
    return(pairing(execution(a,lambda,q0,alpha,phi),expsigna,key))
  }
  
  #pairing(execution(a,lambda,q0,alpha,phi),expsigna,key)
  
  optimum <- optim(rep(0,length(words)),costfunction,control = list(fnscale = -1))
  
  return(c(paste("Cost:",optimum$value,"Coefficients:"),optimum$par))
}

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

#plugged this in manually
sol <- c(paste(-4.63154540495257*10^(-8),""),paste(-0.0515163150843151,"1"),paste(-0.0154001681435968,"2"),paste(0.0126901683658508,"11"),paste(0.034538207873199,"12"),paste(0.0515311274900253,"21"),paste(0.0742716522713664,"22"))

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
  tradstrataverage <- c(tradstrataverage,pairing(sol,s,keys(2,truncorder)))
  costaverage <- c(costaverage,pairing(execution(sol,lambda,q0,alpha,phi),s,keys(2,truncorder)))
}
wealthaverage <- (averageaugmentedprice[2,] - lambda*tradstrataverage)*tradstrataverage
wealthaverage <- cumsum(wealthaverage)
inventoryaverage <- -cumsum(tradstrataverage)
wealthaverage

### Possible figures ###

chosenpaths <- round(runif(5,1,10000)) 
j <- chosenpaths[1]
tradstrat <- c()
cost <- c()
for (i in 1:dim(augmentedstock)[2]) {
  s <- c(1,sig(array(augmentedstock[,1:i,j],c(2,i)),truncorder))
  tradstrat <- c(tradstrat,pairing(sol,s,keys(2,truncorder)))
  cost <- c(cost,pairing(execution(sol,lambda,q0,alpha,phi),s,keys(2,truncorder)))
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
    #tradstrat <- c(tradstrat,pairing(sol,s,keys(2,truncorder)))
    cost <- c(cost,pairing(execution(sol,lambda,q0,alpha,phi),s,keys(2,truncorder)))
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


### Auxiliary functions ###

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


### Get the words corresponding to the terms in the signature ###

keys <- function(d,truncorder){ #applied brute force for now (any help on how to do this in a more elegant way is welcome)
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
