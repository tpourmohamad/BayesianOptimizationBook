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
### Chapter 3 Code


### Section 3.1 ###
### Figure 3.1

f <- function(x) -(cos(5 * x) + 2 * sin(x))
x <- seq(0, 10, .01)  

par( ps = 15 )
plot(x, f(x), type = "l", ylab = "Amount of Contaminants, f(x)")
points(4.42124, f(4.42124), pch = 21, bg = "red")



####################################################
### Section 3.3.1 ###
### Figures 3.2 -- 3.3

f <- function(x) exp(-1.4 * x) * cos(7 * pi * x / 2)
n <- 6
x <- lhs(n, c(0, 1))  
fx <- f(x)

da <- darg(list(mle = TRUE), x)
gp <- newGP(x, fx, da$start, g = .001*var(fx), dK = TRUE)
xx <- matrix(seq(0, 1, .001), ncol = 1)
preds <- predGP(gp, xx)
deleteGP(gp)

fmin <- min(fx)
mu <- preds$mean
sig <- sqrt(diag(preds$Sigma))
probimpr <- pnorm((fmin - mu) / sig)

par(mfrow = c(1, 2), ps = 15)
plot(x, fx, pch = 21, bg = "red", ylim = c(-1, 1), 
     xlab = "x", ylab = "f(x)", main = "Predictive Surface")
lines(xx, mu, col = "red")
lines(xx, mu + 1.96 * sig, lty = 2, col = "red")
lines(xx, mu - 1.96 * sig, lty = 2, col = "red")
legend("topright", c("Mean", "95% PI", "Data"), col = "red",
       lty = c(1, 2, NA), pch = c(NA, NA, 19), bty = "n")
plot(xx, probimpr, type = "l", xlab = "x", ylab = expression(a[PI]*"(x)"), 
     main = "Probability of Improvement")

# Next best point to choose
next.point <- xx[which.max(probimpr)]
x <- matrix(c(x, next.point), ncol = 1)
fx <- c(fx, f(next.point))

da <- darg(list(mle = TRUE), x)
gp <- newGP(x, fx, da$start, g = .001*var(fx), dK = TRUE)
xx <- matrix(seq(0, 1, .001), ncol = 1)
preds <- predGP(gp, xx)
deleteGP(gp)

fmin <- min(fx)
mu <- preds$mean
sig <- sqrt(diag(preds$Sigma))
probimpr <- pnorm((fmin - mu) / sig)

par(mfrow = c(1, 2), ps = 15)
plot(x, fx, pch = 21, bg = "red", ylim = c(-1, 1), 
     xlab = "x", ylab = "f(x)", main = "Predictive Surface")
lines(xx, mu, col = "red")
lines(xx, mu + 1.96 * sig, lty = 2, col = "red")
lines(xx, mu - 1.96 * sig, lty = 2, col = "red")
legend("topright", c("Mean", "95% PI", "Data"), col = "red",
       lty = c(1, 2, NA), pch = c(NA, NA, 19), bty = "n")
plot(xx, probimpr, type = "l", xlab = "x", ylab = expression(a[PI]*"(x)"), 
     main = "Probability of Improvement")


### Figures 3.4 -- 3.6
f <- function(x) (cos(5 * x) + 2 * sin(x))
n <- 10

xx <- matrix(seq(0, 10, .01), ncol = 1)
x <- lhs(n, c(0, 10))  
fx <- f(x)

da <- darg(list(mle = TRUE), x)
gp <- newGP(x, fx, da$start, g = .0001*var(fx), dK = TRUE)
preds <- predGP(gp, xx)
deleteGP(gp)

fmin <- min(fx)
mu <- preds$mean
sig <- sqrt(diag(preds$Sigma))
PI <- pnorm((fmin - mu) / sig)

par( ps = 15, mfrow = c(1,2) )
plot(xx, f(xx), type = "l", xlab = "x", ylab = "Amount of Contaminants, f(x)",
     main = "Predictive Surface for n = 10")
lines(xx, mu, col = "red")
lines(xx, mu + 1.96 * sig, lty = 2, col = "red")
lines(xx, mu - 1.96 * sig, lty = 2, col = "red")
points(x, fx, pch = 21, bg = "red", cex = 1.2)

next.point <- xx[which.max(PI),1]
x <- matrix(c(x, next.point), ncol = 1)
fx <- c(fx, f(next.point)) 

points(next.point, f(next.point), pch = 21, bg = "deepskyblue", cex = 1.2)

plot(xx, PI, type = "l", xlab = "x", ylab = expression(a[PI]*"(x)"),
     main = "Probability of Improvement")

for(i in 1:3){
  ga <- garg(list(mle = TRUE), fx)
  gp <- newGP(x, fx, da$start, g = .0001*var(fx), dK = TRUE)
  preds <- predGP(gp, xx)
  deleteGP(gp)
  
  fmin <- min(fx)
  mu <- preds$mean
  sig <- sqrt(diag(preds$Sigma))
  PI <- pnorm((fmin - mu) / sig)
  
  next.point <- xx[which.max(PI),1]
  x <- matrix(c(x, next.point), ncol = 1)
  fx <- c(fx, f(next.point)) 
  
}

par( ps = 15, mfrow = c(1,2) )
plot(xx, f(xx), type = "l", xlab = "x", ylab = "Amount of Contaminants, f(x)",
     main = "Predictive Surface for n = 13")
lines(xx, mu, col = "red")
lines(xx, mu + 1.96 * sig, lty = 2, col = "red")
lines(xx, mu - 1.96 * sig, lty = 2, col = "red")
points(x, fx, pch = 21, bg = "red", cex = 1.2)

next.point <- xx[which.max(PI),1]
x <- matrix(c(x, next.point), ncol = 1)
fx <- c(fx, f(next.point)) 

points(next.point, f(next.point), pch = 21, bg = "deepskyblue", cex = 1.2)

plot(xx, PI, type = "l", xlab = "x", ylab = expression(a[PI]*"(x)"),
     main = "Probability of Improvement")

for(i in 1:2){
  da <- darg(list(mle = TRUE), x)
  gp <- newGP(x, fx, da$start, g = .0001*var(fx), dK = TRUE)
  preds <- predGP(gp, xx)
  deleteGP(gp)
  
  fmin <- min(fx)
  mu <- preds$mean
  sig <- sqrt(diag(preds$Sigma))
  PI <- pnorm((fmin - mu) / sig)
  
  next.point <- xx[which.max(PI),1]
  x <- matrix(c(x, next.point), ncol = 1)
  fx <- c(fx, f(next.point)) 
  
}

par( ps = 15, mfrow = c(1,2) )
plot(xx, f(xx), type = "l", xlab = "x", ylab = "Amount of Contaminants, f(x)",
     main = "Predictive Surface for n = 15")
lines(xx, mu, col = "red")
lines(xx, mu + 1.96 * sig, lty = 2, col = "red")
lines(xx, mu - 1.96 * sig, lty = 2, col = "red")
points(x, fx, pch = 21, bg = "red", cex = 1.2)

next.point <- xx[which.max(PI),1]
x <- matrix(c(x, next.point), ncol = 1)
fx <- c(fx, f(next.point)) 

points(next.point, f(next.point), pch = 21, bg = "deepskyblue")

plot(xx, PI, type = "l", xlab = "x", ylab = expression(a[PI]*"(x)"),
     main = "Probability of Improvement", cex = 1.2)


### Figure 3.7
bov <- function(y, end = length(y)){ # Calculates best overall value
  prog <- rep(min(y), end)
  prog[1:min(end, length(y))] <- y[1:min(end, length(y))]
  for(i in 2:end)
    if(is.na(prog[i]) || prog[i] > prog[i-1]) prog[i] <- prog[i-1]
  return(prog)
}

f <- function(x) (cos(5 * x) + 2 * sin(x))
n <- 10
xx <- matrix(seq(0, 10, .01), ncol = 1)

runs <- 30
budget <- 50
solutions <- matrix(NA, nrow = runs, ncol = budget)

for(j in 1:runs){
  x <- lhs(n, c(0, 10))  
  fx <- f(x)
  for(i in 1:(budget-length(x))){
    ga <- garg(list(mle = TRUE), fx)
    gp <- newGP(x, fx, da$start, g = .0001*var(fx), dK = TRUE)
    preds <- predGP(gp, xx)
    deleteGP(gp)

    fmin <- min(fx)
    mu <- preds$mean
    sig <- sqrt(diag(preds$Sigma))
    PI <- pnorm((fmin - mu) / sig)

    next.point <- xx[which.max(PI),1]
    x <- matrix(c(x, next.point), ncol = 1)
    fx <- c(fx, f(next.point)) 
    print(c(j,i))
  }
  solutions[j, ] <- fx
}

best <- apply(solutions, 1, bov)
avg <- apply(best, 1, mean)
  
par( ps = 15 )
matplot(1:50,best, col = "gray", lwd = 0.5, lty = 1, type = "l",
        xlab = "black-box evaluations (n)", ylab = "best objective value",
        main = "Intializing with a Sample Size of n = 10")
lines(1:50, avg, lwd = 2)
abline(h = f(4.42124), lty = 2, lwd = 2)
legend("topright", c("Average Solution","Global Solution"), lty = 1:2, bty = "n")


### Figure 3.8
f <- function(x) (cos(5 * x) + 2 * sin(x))
n <- 15
xx <- matrix(seq(0, 10, .01), ncol = 1)

runs <- 30
budget <- 50
solutions <- matrix(NA, nrow = runs, ncol = budget)

for(j in 1:runs){
  x <- lhs(n, c(0, 10))  
  fx <- f(x)
  for(i in 1:(budget-length(x))){
    da <- darg(list(mle = TRUE), x)
    gp <- newGP(x, fx, da$start, g = .0001*var(fx), dK = TRUE)
    preds <- predGP(gp, xx)
    deleteGP(gp)
    
    fmin <- min(fx)
    mu <- preds$mean
    sig <- sqrt(diag(preds$Sigma))
    PI <- pnorm((fmin - mu) / sig)
    
    next.point <- xx[which.max(PI),1]
    x <- matrix(c(x, next.point), ncol = 1)
    fx <- c(fx, f(next.point)) 
    print(c(j,i))
  }
  solutions[j, ] <- fx
}

best <- apply(solutions, 1, bov)
avg <- apply(best, 1, mean)

par( ps = 15 )
matplot(1:50,best, col = "gray", lwd = 0.5, lty = 1, type = "l",
        xlab = "black-box evaluations (n)", ylab = "best objective value",
        main = "Intializing with a Sample Size of n = 15")
lines(1:50, avg, lwd = 2)
abline(h = f(4.42124), lty = 2, lwd = 2)
legend("topright", c("Average Solution","Global Solution"), lty = 1:2, bty = "n")



