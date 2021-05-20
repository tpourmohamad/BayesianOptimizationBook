####################################################
### Appendix A.1 -- Getting Started

library(GA)
library(tgp)
library(laGP)
library(mvtnorm)
library(ggplot2)
library(gridExtra)
library(CompModels)
library(scatterplot3d)

####################################################
### Chapter 4 Code


### Section 4.2.1 ###
### Constrained Expected Improvement

### Figure 4.1

# The objective function
obj = function(x1,x2){
  value <- 4*x1^2 - x1 - x2 - 2.5
  return(value)
}

# The constraint functions
con1 = function(x1,x2){
  value <- -x2^2 + 1.5*x1^2 - 2*x1 + 1
  return(value)
}

con2 = function(x1,x2){
  value <- x2^2 + 3*x1^4 - 2*x1 - 4.25
    return(value)
}

# Objective and constraint function for use 
# with laGP package optimization algorithms
bbox <- function(X){
  obj <- 4*X[1]^2 - X[1] - X[2] - 2.5
  con1 <- -X[2]^2 + 1.5*X[1]^2 - 2*X[1] + 1
  con2 <- X[2]^2 + 3*X[1]^4 - 2*X[1] - 4.25
  con <- c(con1, con2)
  return(list(obj = obj, c = con))
}

n <- 200
x1 <- seq(-1.5, 2.5, len = n)
x2 <- seq(-3, 3, len = n)

x <- expand.grid(x1, x2)
obj <- rep(NA, nrow(x))
con <- matrix(NA, nrow = nrow(x), ncol = 2)

for(i in 1:nrow(x)){
  temp <- unlist(bbox(x[i,]))
  obj[i] <- temp[1]
  con[i,] <- c(temp[2], temp[3])
}

y <- obj
y[con[,1] > 0 | con[,2] > 0] <- NA

z <- obj
z[!(con[,1] > 0 | con[,2] > 0)] <- NA

par(ps=15)
plot(0, 0, type = "n", xlim = c(-1.5, 2.5), ylim = c(-3, 3), 
     xlab = expression(x[1]), ylab = expression(x[2]))
c1 <- matrix(con[,1], ncol = n)
contour(x1, x2, c1, nlevels = 1, levels = 0, drawlabels = FALSE, add = TRUE, lwd = 2)
c2 <- matrix(con[,2], ncol = n)
contour(x1, x2, c2, nlevels = 1, levels = 0, drawlabels = FALSE, add = TRUE, lwd = 2, lty = 2)
contour(x1, x2, matrix(y, ncol = n), nlevels = 10, add = TRUE, col = "forestgreen")
contour(x1, x2, matrix(z, ncol = n), nlevels = 25, add = TRUE, col = 2, lty = 2)
points(0.1707747, 2.141688, pch = 21, bg = "deepskyblue")


### Figures 4.2 -- 4.3

B <- matrix(c(-1.5, -3, 2.5, 3), ncol=2)
x <- lhs(n = 10, rect = B)
out.CEI <- optim.efi(bbox, B, fhat = TRUE, end=60, Xstart = x)

par(ps=15)
plot(0, 0, type = "n", xlim = c(-1.5, 2.5), ylim = c(-3, 3), 
     xlab = expression(x[1]), ylab = expression(x[2]),
     main = "CEI")
c1 <- matrix(con[,1], ncol = n)
contour(x1, x2, c1, nlevels = 1, levels = 0, drawlabels = FALSE, add = TRUE, lwd = 2)
c2 <- matrix(con[,2], ncol = n)
contour(x1, x2, c2, nlevels = 1, levels = 0, drawlabels = FALSE, add = TRUE, lwd = 2, lty = 2)
contour(x1, x2, matrix(y, ncol = n), nlevels = 10, add = TRUE, col = "forestgreen")
contour(x1, x2, matrix(z, ncol = n), nlevels = 25, add = TRUE, col = 2, lty = 2)

best <- which.min(out.CEI$prog)

points(out.CEI$X[-c(1:10),1],out.CEI$X[-c(1:10),2],pch=19)
points(out.CEI$X[1:10,1],out.CEI$X[1:10,2],pch="x", col ="purple")
points(out.CEI$X[best,1],out.CEI$X[best,2],pch=21,bg="deepskyblue",cex=1.5)
legend("bottomleft", c("Initial LHS", "Data", "Minimum Found"),
       pch = c(4, 19, 19), col = c("purple", "black", "deepskyblue"),
       pt.cex = c(2,1,1), bg = "white")


S <- 50 # The number of additional inputs to select
n.results = 30 # Number of Monte Carlo runs

out.CEI <- matrix(NA,nrow=n.results,ncol=(10+S))
B <- matrix(c(-1.5, -3, 2.5, 3), ncol=2)

for(result.iter in 1:n.results){
  
  set.seed(result.iter)
  x <- lhs(n = 10, rect = B)
  out.CEI[result.iter,] <- optim.efi(bbox, B, fhat = TRUE, end=60, Xstart = x)$prog
  
}

E.CEI <- apply(out.CEI,2,mean)

par(ps=15)
matplot(t(out.CEI),  col="grey", lty = 1,  type="l",
        xlab="black-box evalulations (n)",
        ylab="best valid objective (f)",
        main = "CEI",
        ylim=c(-5,1))
lines(E.CEI, lwd = 2, col ="red")
legend("topright",c("Average Solution","Global Solution"),
       lty=c(1,2),lwd=c(2,2),bty="n", col=c("red","black"))
abline(h=-4.6958,lty=2,lwd=2,col="black")



####################################################
### Section 4.2.2 ###
### Asymmetric Entropy

### Figures 4.4 -- 4.5

n <- 200
x1 <- seq(-1.5, 2.5, len = n)
x2 <- seq(-3, 3, len = n)

x <- expand.grid(x1, x2)
obj <- rep(NA, nrow(x))
con <- matrix(NA, nrow = nrow(x), ncol = 2)

for(i in 1:nrow(x)){
  temp <- unlist(bbox(x[i,]))
  obj[i] <- temp[1]
  con[i,] <- c(temp[2], temp[3])
}

y <- obj
y[con[,1] > 0 | con[,2] > 0] <- NA

z <- obj
z[!(con[,1] > 0 | con[,2] > 0)] <- NA

par(ps=15)
plot(0, 0, type = "n", xlim = c(-1.5, 2.5), ylim = c(-3, 3), 
     xlab = expression(x[1]), ylab = expression(x[2]),
     main = "AE")
c1 <- matrix(con[,1], ncol = n)
contour(x1, x2, c1, nlevels = 1, levels = 0, drawlabels = FALSE, add = TRUE, lwd = 2)
c2 <- matrix(con[,2], ncol = n)
contour(x1, x2, c2, nlevels = 1, levels = 0, drawlabels = FALSE, add = TRUE, lwd = 2, lty = 2)
contour(x1, x2, matrix(y, ncol = n), nlevels = 10, add = TRUE, col = "forestgreen")
contour(x1, x2, matrix(z, ncol = n), nlevels = 25, add = TRUE, col = 2, lty = 2)


# The objective function
obj = function(x1,x2){
  value <- 4*x1^2 - x1 - x2 - 2.5
  return(value)
}

# The constraint functions
con1 = function(x1,x2){
  value <- -x2^2 + 1.5*x1^2 - 2*x1 + 1
  return(value)
}

con2 = function(x1,x2){
  value <- x2^2 + 3*x1^4 - 2*x1 - 4.25
  return(value)
}


B <- matrix(c(-1.5, -3, 2.5, 3), ncol=2)
x <- lhs(n = 10, rect = B)

z = obj(x[,1],x[,2])
w = con1(x[,1],x[,2])
v = con2(x[,1],x[,2])
xx = lhs(n=1000,rect=B)


for(i in 1:50){
  
  model.f = laGP(xx, 6, 8+i, x, z)   
  model.c1 = laGP(xx, 6, 8+i, x, w)  
  model.c2 = laGP(xx, 6, 8+i, x, v)  
  
  mu.f =  model.f$mean
  mu.c1 = model.c1$mean
  mu.c2 = model.c2$mean
  
  sigma.f = model.f$s2
  sigma.c1 = model.c1$s2
  sigma.c2 = model.c2$s2
  
  p <- pnorm(0,mean=mu.c1,sd = sqrt(sigma.c1))*
    pnorm(0,mean=mu.c2,sd = sqrt(sigma.c2))
  
  fmin <- min(z[w<=0 & v <= 0])
  EI <- (fmin - mu.f) * pnorm((fmin - mu.f) / sqrt(sigma.f)) +
         sqrt(sigma.f) * dnorm((fmin - mu.f) / sqrt(sigma.f))
  ww <- 2 / 3
  Sx <- (2 * p * (1 - p)) / (p - 2 * ww * p + ww^2)
  acq <- (EI)*(Sx^5)        
  
  idx = which.max(acq)
  
  x = rbind(x,xx[idx,])
  z = c(z,obj(xx[idx,1],xx[idx,2]))
  w = c(w,con1(xx[idx,1],xx[idx,2]))
  v = c(v,con2(xx[idx,1],xx[idx,2]))
  
  xx = lhs(n=1000,rect=B)
  print(paste(i))
  
}

fmin <- min(z[w<=0 & v <= 0])
best = which(z == fmin)

points(x[-c(1:10),1],x[-c(1:10),2],pch=19)
points(x[1:10,1],x[1:10,2],pch="x", col = "purple")
points(x[best,1],x[best,2],pch=21,bg="deepskyblue",cex=1.5)
legend("bottomleft", c("Initial LHS", "Data", "Minimum Found"),
       pch = c(4, 19, 19), col = c("purple", "black", "deepskyblue"),
       pt.cex = c(2,1,1), bg = "white")

RESULTS.f = matrix(NA,nrow=n.results,ncol=(10+S))
RESULTS.c1 = matrix(NA,nrow=n.results,ncol=(10+S))
RESULTS.c2 = matrix(NA,nrow=n.results,ncol=(10+S))
RESULTS.x = matrix(NA,nrow=n.results,ncol=(10+S))

for(result.iter in 1:n.results){
  
  set.seed(result.iter) 
  
  x = lhs(n=10,rect=B) 
  z = obj(x[,1],x[,2])
  w = con1(x[,1],x[,2])
  v = con2(x[,1],x[,2])
  xx = lhs(n=1000,rect=B)
  
  for(i in 1:50){
    
    model.f = laGP(xx, 6, 8+i, x, z)   
    model.c1 = laGP(xx, 6, 8+i, x, w)  
    model.c2 = laGP(xx, 6, 8+i, x, v)  
    
    mu.f =  model.f$mean
    mu.c1 = model.c1$mean
    mu.c2 = model.c2$mean
    
    sigma.f = model.f$s2
    sigma.c1 = model.c1$s2
    sigma.c2 = model.c2$s2
    
    p <- pnorm(0,mean=mu.c1,sd = sqrt(sigma.c1))*
         pnorm(0,mean=mu.c2,sd = sqrt(sigma.c2))
    
    fmin <- min(z[w<=0 & v <= 0])
    EI <- (fmin - mu.f) * pnorm((fmin - mu.f) / sqrt(sigma.f)) +
          sqrt(sigma.f) * dnorm((fmin - mu.f) / sqrt(sigma.f))
    ww <- 2 / 3
    Sx <- (2 * p * (1 - p)) / (p - 2 * ww * p + ww^2)
    acq <- (EI)*(Sx^5)        
    
    idx = which.max(acq)
    
    x = rbind(x,xx[idx,])
    z = c(z,obj(xx[idx,1],xx[idx,2]))
    w = c(w,con1(xx[idx,1],xx[idx,2]))
    v = c(v,con2(xx[idx,1],xx[idx,2]))
    
    xx = lhs(n=1000,rect=B)
    print(paste(i,result.iter))
    
  }
  
  RESULTS.f[result.iter,] = z
  RESULTS.c1[result.iter,] = w
  RESULTS.c2[result.iter,] = v
  
}

# Calculate the best feasible objective value over the function evalutions 
fx = RESULTS.f
for(i in 1:nrow(fx)){
  
  idx = which(RESULTS.c1[i,]<=0 & RESULTS.c2[i,]<=0)
  if(idx[1] != 1){
    fx[i,1:(idx[1]-1)] = Inf
  }
  hx = rep(0,ncol(fx))
  hx[idx] = RESULTS.f[i,idx]
  
  for(j in 1:(ncol(fx)-1)){
    
    if(RESULTS.f[i,j+1]<fx[i,j] & hx[j+1]!=0 ){
      fx[i,j+1] = RESULTS.f[i,j+1]
    }else{
      fx[i,j+1] = fx[i,j]
    }
    
  }
  
}

E.AE = apply(fx,2,mean)

par(ps=15)
matplot(t(fx),  col="grey", lty = 1,  type="l",
        xlab="black-box evalulations (n)",
        ylab="best valid objective (f)",
        main = "AE", ylim=c(-5,1))
lines(E.AE, lwd = 2, col = "red")
axis(1)
axis(2)
legend("topright",c("Average Solution","Global Solution"),
       lty=c(1,2),lwd=c(2,2),bty="n", col=c("red","black"))
abline(h=-4.6958,lty=2,lwd=2,col="black")



####################################################
### Section 4.2.3 ###
### Augmented Lagrangian

### Figures 4.6 -- 4.7

n <- 200
x1 <- seq(-1.5, 2.5, len = n)
x2 <- seq(-3, 3, len = n)

x <- expand.grid(x1, x2)
obj <- rep(NA, nrow(x))
con <- matrix(NA, nrow = nrow(x), ncol = 2)

for(i in 1:nrow(x)){
  temp <- unlist(bbox(x[i,]))
  obj[i] <- temp[1]
  con[i,] <- c(temp[2], temp[3])
}

y <- obj
y[con[,1] > 0 | con[,2] > 0] <- NA

z <- obj
z[!(con[,1] > 0 | con[,2] > 0)] <- NA

B <- matrix(c(-1.5, -3, 2.5, 3), ncol=2)
x <- lhs(n = 10, rect = B)
out.AL <- optim.auglag(bbox, B, fhat = TRUE, end=60, Xstart = x)

par(ps=15)
plot(0, 0, type = "n", xlim = c(-1.5, 2.5), ylim = c(-3, 3), 
     xlab = expression(x[1]), ylab = expression(x[2]),
     main = "AL")
c1 <- matrix(con[,1], ncol = n)
contour(x1, x2, c1, nlevels = 1, levels = 0, drawlabels = FALSE, add = TRUE, lwd = 2)
c2 <- matrix(con[,2], ncol = n)
contour(x1, x2, c2, nlevels = 1, levels = 0, drawlabels = FALSE, add = TRUE, lwd = 2, lty = 2)
contour(x1, x2, matrix(y, ncol = n), nlevels = 10, add = TRUE, col = "forestgreen")
contour(x1, x2, matrix(z, ncol = n), nlevels = 25, add = TRUE, col = 2, lty = 2)

best <- which.min(out.AL$prog)

points(out.AL$X[-c(1:10),1],out.AL$X[-c(1:10),2],pch=19)
points(out.AL$X[1:10,1],out.AL$X[1:10,2],pch="x", col ="purple")
points(out.AL$X[best,1],out.AL$X[best,2],pch=21,bg="deepskyblue",cex=1.5)
legend("bottomleft", c("Initial LHS", "Data", "Minimum Found"),
       pch = c(4, 19, 19), col = c("purple", "black", "deepskyblue"),
       pt.cex = c(2,1,1), bg = "white")

S <- 50 # The number of additional inputs to select
n.results = 30 # Number of Monte Carlo runs

out.AL <- matrix(NA,nrow=n.results,ncol=(10+S))
B <- matrix(c(-1.5, -3, 2.5, 3), ncol=2)

for(result.iter in 1:n.results){
  
  set.seed(result.iter)
  x <- lhs(n = 10, rect = B)
  out.AL[result.iter,] <- optim.auglag(bbox, B, fhat = TRUE, end=60, Xstart = x)$prog
  
}

E.AL <- apply(out.AL,2,mean)

par(ps=15)
matplot(t(out.AL),  col="grey", lty = 1,  type="l",
        xlab="black-box evalulations (n)",
        ylab="best valid objective (f)",
        main = "AL",
        ylim=c(-5,1))
lines(E.AL, lwd = 2, col = "red")
legend("topright",c("Average Solution","Global Solution"),
       lty=c(1,2),lwd=c(2,2),bty="n", col=c("red","black"))
abline(h=-4.6958,lty=2,lwd=2,col="black")



####################################################
### Section 4.2.4 ###
### Barrier Method

### Figure 4.8

### Figures 4.9 -- 4.10

n <- 200
x1 <- seq(-1.5, 2.5, len = n)
x2 <- seq(-3, 3, len = n)

x <- expand.grid(x1, x2)
obj <- rep(NA, nrow(x))
con <- matrix(NA, nrow = nrow(x), ncol = 2)

for(i in 1:nrow(x)){
  temp <- unlist(bbox(x[i,]))
  obj[i] <- temp[1]
  con[i,] <- c(temp[2], temp[3])
}

y <- obj
y[con[,1] > 0 | con[,2] > 0] <- NA

z <- obj
z[!(con[,1] > 0 | con[,2] > 0)] <- NA

par(ps=15)
plot(0, 0, type = "n", xlim = c(-1.5, 2.5), ylim = c(-3, 3), 
     xlab = expression(x[1]), ylab = expression(x[2]),
     main = "BM")
c1 <- matrix(con[,1], ncol = n)
contour(x1, x2, c1, nlevels = 1, levels = 0, drawlabels = FALSE, add = TRUE, lwd = 2)
c2 <- matrix(con[,2], ncol = n)
contour(x1, x2, c2, nlevels = 1, levels = 0, drawlabels = FALSE, add = TRUE, lwd = 2, lty = 2)
contour(x1, x2, matrix(y, ncol = n), nlevels = 10, add = TRUE, col = "forestgreen")
contour(x1, x2, matrix(z, ncol = n), nlevels = 25, add = TRUE, col = 2, lty = 2)


# The objective function
obj = function(x1,x2){
  value <- 4*x1^2 - x1 - x2 - 2.5
  return(value)
}

# The constraint functions
con1 = function(x1,x2){
  value <- -x2^2 + 1.5*x1^2 - 2*x1 + 1
  return(value)
}

con2 = function(x1,x2){
  value <- x2^2 + 3*x1^4 - 2*x1 - 4.25
  return(value)
}

B <- matrix(c(-1.5, -3, 2.5, 3), ncol=2)
x <- lhs(n = 10, rect = B)

z = obj(x[,1],x[,2])
w = con1(x[,1],x[,2])
v = con2(x[,1],x[,2])
xx = lhs(n=1000,rect=B)


for(i in 1:50){
  
  model.f = laGP(xx, 6, 8+i, x, z)   
  model.c1 = laGP(xx, 6, 8+i, x, w)  
  model.c2 = laGP(xx, 6, 8+i, x, v)  
  
  mu.f =  model.f$mean
  mu.c1 = model.c1$mean
  mu.c2 = model.c2$mean
  
  sigma.f = model.f$s2
  sigma.c1 = model.c1$s2
  sigma.c2 = model.c2$s2
  
  gamma = 1/sigma.f
  E = mu.f - (1/gamma)*((log(-mu.c1) + sigma.c1/(2*mu.c1^2)) + 
                          (log(-mu.c2) + sigma.c2/(2*mu.c2^2)))
  
  idx = which.min(E)
  
  x = rbind(x,xx[idx,])
  z = c(z,obj(xx[idx,1],xx[idx,2]))
  w = c(w,con1(xx[idx,1],xx[idx,2]))
  v = c(v,con2(xx[idx,1],xx[idx,2]))
  
  xx = lhs(n=1000,rect=B)
  print(paste(i))
  
}

fmin <- min(z[w<=0 & v <= 0])
best = which(z == fmin)

points(x[-c(1:10),1],x[-c(1:10),2],pch=19)
points(x[1:10,1],x[1:10,2],pch=4, col = "purple")
points(x[best,1],x[best,2],pch=21,bg="deepskyblue",cex=1.5)
legend("bottomleft", c("Initial LHS", "Data", "Minimum Found"),
       pch = c(4, 19, 19), col = c("purple", "black", "deepskyblue"),
       pt.cex = c(2,1,1), bg = "white")

RESULTS.f = matrix(NA,nrow=n.results,ncol=(10+S))
RESULTS.c1 = matrix(NA,nrow=n.results,ncol=(10+S))
RESULTS.c2 = matrix(NA,nrow=n.results,ncol=(10+S))
RESULTS.x = matrix(NA,nrow=n.results,ncol=(10+S))

for(result.iter in 1:n.results){
  
  set.seed(result.iter) 
  
  x = lhs(n=10,rect=B) 
  z = obj(x[,1],x[,2])
  w = con1(x[,1],x[,2])
  v = con2(x[,1],x[,2])
  xx = lhs(n=1000,rect=B)
  
  for(i in 1:50){
    
    model.f = laGP(xx, 6, 8+i, x, z)   
    model.c1 = laGP(xx, 6, 8+i, x, w)  
    model.c2 = laGP(xx, 6, 8+i, x, v)  
    
    mu.f =  model.f$mean
    mu.c1 = model.c1$mean
    mu.c2 = model.c2$mean
    
    sigma.f = model.f$s2
    sigma.c1 = model.c1$s2
    sigma.c2 = model.c2$s2
    
    gamma = 1/sigma.f
    E = mu.f - (1/gamma)*((log(-mu.c1) + sigma.c1/(2*mu.c1^2)) + 
                            (log(-mu.c2) + sigma.c2/(2*mu.c2^2)))
    
    idx = which.min(E)
    
    x = rbind(x,xx[idx,])
    z = c(z,obj(xx[idx,1],xx[idx,2]))
    w = c(w,con1(xx[idx,1],xx[idx,2]))
    v = c(v,con2(xx[idx,1],xx[idx,2]))
    
    xx = lhs(n=1000,rect=B)
    print(paste(i,result.iter))
    
  }
  
  RESULTS.f[result.iter,] = z
  RESULTS.c1[result.iter,] = w
  RESULTS.c2[result.iter,] = v
  
}

# Calculate the best feasible objective value over the function evalutions 
fx = RESULTS.f
for(i in 1:nrow(fx)){
  
  idx = which(RESULTS.c1[i,]<=0 & RESULTS.c2[i,]<=0)
  if(idx[1] != 1){
    fx[i,1:(idx[1]-1)] = Inf
  }
  hx = rep(0,ncol(fx))
  hx[idx] = RESULTS.f[i,idx]
  
  for(j in 1:(ncol(fx)-1)){
    
    if(RESULTS.f[i,j+1]<fx[i,j] & hx[j+1]!=0 ){
      fx[i,j+1] = RESULTS.f[i,j+1]
    }else{
      fx[i,j+1] = fx[i,j]
    }
    
  }
  
}

E.BM = apply(fx,2,mean)

par(ps=15)
matplot(t(fx),  col="grey", lty = 1,  type="l",
        xlab="black-box evalulations (n)",
        ylab="best valid objective (f)",
        main = "BM", ylim=c(-5,1))
lines(E.BM, lwd = 2, col = "red")
axis(1)
axis(2)
legend("topright",c("Average Solution","Global Solution"),
       lty=c(1,2),lwd=c(2,2),bty="n", col=c("red","black"))
abline(h=-4.6958,lty=2,lwd=2,col="black")



### Figure 4.11
par( ps = 15)
plot(E.CEI, col=1, lwd=2, type="l",
     xlab="black-box evaluations (n)", ylab="average best valid objective value",
     ylim = c(-5, -3))
lines(E.AE, col="red", lwd=2)
lines(E.AL, col="green", lwd=2)
lines(E.BM, col="deepskyblue", lwd=2)
abline(v=10, lty=3)
abline(h=-4.6958,lty=2,lwd=2,col="black")
legend("topright", c(expression(a[CEI]*"(x)"), expression(a[AE]*"(x)"), 
                     expression(a[AL]*"(x)"), expression(a[BM]*"(x)"),
                     "Initial n", "Global Solution"),
       col=c("black", "red", "green", "deepskyblue","black","black"), lwd=c(2,2,2,2,1,2), 
       lty = c(1,1,1,1,3,2), bty="n")


####################################################
### Section 4.3 ###

### Figure 4.12 -- 4.13


