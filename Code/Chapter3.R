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



####################################################
### Section 3.3.2 ###
### Figures 3.10 -- 3.13

f <- function(x) -(1.4 - 3 * x) * sin(18 * x)
n <- 10
x <- lhs(n, c(0, 1.2))  
fx <- f(x)

da <- darg(list(mle = TRUE), x)
gp <- newGP(x, fx, da$start, g = .0001*var(fx), dK = TRUE)
xx <- matrix(seq(0, 1.2, .001), ncol = 1)
preds <- predGP(gp, xx)
deleteGP(gp)

fmin <- min(fx)
mu <- preds$mean
sig <- sqrt(diag(preds$Sigma))
ei <- (fmin - mu) * pnorm((fmin - mu) / sig) + sig * dnorm ((fmin - mu) / sig)

par(mfrow = c(1, 2), ps = 15)
plot(x, fx, pch = 21, bg = "red", ylim = c(-1.5, 2.5), xlim = c(0, 1.2), 
     xlab = "x", ylab = "f(x)", main = "Predictive Surface")
lines(xx, mu, col = "red")
lines(xx, mu + 1.96 * sig, lty = 2, col = "red")
lines(xx, mu - 1.96 * sig, lty = 2, col = "red")
legend("topleft", c("Mean", "95% PI", "Data"), col = "red",
       lty = c(1, 2, NA), pch = c(NA, NA, 19), bty = "n")
plot(xx, ei, type = "l", xlab = "x", ylab = expression(a[EI]*"(x)"), 
     main = "Expected Improvement")

par( ps = 15 )
xgrid <- seq(0, 1.2, .001)
plot(xgrid, f(xgrid), xlab = "x", ylab = "f(x)", type = "l",
     main = "True Objective Function")
points(x, f(x), bg = "red", pch = 21)
legend("toplef", "Data", pch = 19, col = "red", bty = "n")

# Next best point to choose
next.point <- xx[which.max(ei)]
x <- matrix(c(x, next.point), ncol = 1)
fx <- c(fx, f(next.point))

da <- darg(list(mle = TRUE), x)
gp <- newGP(x, fx, da$start, g = .0001*var(fx), dK = TRUE)
xx <- matrix(seq(0, 1.2, .001), ncol = 1)
preds <- predGP(gp, xx)
deleteGP(gp)

fmin <- min(fx)
mu <- preds$mean
sig <- sqrt(diag(preds$Sigma))
ei <- (fmin - mu) * pnorm((fmin - mu) / sig) + sig * dnorm ((fmin - mu) / sig)

par(mfrow = c(1, 2), ps = 15)
plot(x, fx, pch = 21, bg = "red", ylim = c(-1.5, 2.5), 
     xlab = "x", ylab = "f(x)", main = "Predictive Surface")
lines(xx, mu, col = "red")
lines(xx, mu + 1.96 * sig, lty = 2, col = "red")
lines(xx, mu - 1.96 * sig, lty = 2, col = "red")
legend("topleft", c("Mean", "95% PI", "Data"), col = "red",
       lty = c(1, 2, NA), pch = c(NA, NA, 19), bty = "n")
plot(xx, ei, type = "l", xlab = "x", ylab = expression(a[EI]*"(x)"), 
     main = "Expected Improvement")

for(i in 1:9){
  next.point <- xx[which.max(ei)]
  x <- matrix(c(x, next.point), ncol = 1)
  fx <- c(fx, f(next.point))
  
  da <- darg(list(mle = TRUE), x)
  gp <- newGP(x, fx, da$start, g = .0001*var(fx), dK = TRUE)
  xx <- matrix(seq(0, 1.2, .001), ncol = 1)
  preds <- predGP(gp, xx)
  deleteGP(gp)
  
  fmin <- min(fx)
  mu <- preds$mean
  sig <- sqrt(diag(preds$Sigma))
  ei <- (fmin - mu) * pnorm((fmin - mu) / sig) + sig * dnorm ((fmin - mu) / sig)
}

par(mfrow = c(1, 2), ps = 15)
plot(x, fx, pch = 21, bg = "red", ylim = c(-1.5, 2.5), 
     xlab = "x", ylab = "f(x)", main = "Predictive Surface")
lines(xx, mu, col = "red")
lines(xx, mu + 1.96 * sig, lty = 2, col = "red")
lines(xx, mu - 1.96 * sig, lty = 2, col = "red")
legend("topleft", c("Mean", "95% PI", "Data"), col = "red",
       lty = c(1, 2, NA), pch = c(NA, NA, 19), bty = "n")
plot(xx, ei, type = "l", xlab = "x", ylab = expression(a[EI]*"(x)"), 
     main = "Expected Improvement")



### Figures 3.14 -- 3.16
f <- function(x1, x2) -(cos((x1 - .1) * x2))^2 - x1 * sin(3 * x1 + x2)
n <- 10
x <- lhs(n, rbind(c(-2.25, 2.5), c(-2.5, 1.75)))
fx <- f(x[,1], x[,2])

da <- darg(list(mle = TRUE), x)
ga <- garg(list(mle = TRUE), fx)
gp <- newGP(x, fx, da$start, ga$start, dK = TRUE)

x1 <- seq(-2.25, 2.5, .05)
x2 <- seq(-2.5, 1.75, .05)
xx <- expand.grid(x1, x2)
preds <- predGP(gp, xx, nonug=TRUE)
deleteGP(gp)

fmin <- min(fx)
mu <- preds$mean
sig <- sqrt(diag(preds$Sigma))
probimpr <- pnorm((fmin - mu) / sig)
ei <- (fmin - mu) * pnorm((fmin - mu) / sig) + sig * dnorm ((fmin - mu) / sig)

obj <- function(x) f(x[1], x[2])
mini <- optim(c(2.49, 0.09), obj, method = "L-BFGS-B", 
              lower = c(-2.25, -2.5), upper = c(2.5, 1.75))

x <- data.frame(x1 = x[,1], x2 = x[,2])

x1 <- seq(-2.25, 2.5, .01)
x2 <- seq(-2.5, 1.75, .01)
out <- outer(x1, x2, f)

X <- expand.grid(x1, x2)

temp <- data.frame(x = X[,1], y = X[,2], z = c(out))

ggplot(temp, aes(x, y, z = z)) + 
  geom_tile(aes(fill=z)) + 
  geom_contour(color = "black") + 
  scale_x_continuous(expand = c(0, 0)) + 
  scale_y_continuous(expand = c(0, 0)) + 
  scale_fill_gradientn(colors = c("yellow","red"), name = "") + 
  labs(x = expression(x[1]), y = expression(x[2]),
       title = "Objective Function") + 
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 10),
        plot.title = element_text(hjust = 0.5))  

temp <- data.frame(x = xx[,1], y = xx[,2], z = mu)

ggplot(temp, aes(x, y, z = z)) + 
  geom_tile(aes(fill=z)) + 
  geom_contour(color = "black") + 
  geom_point(data = x, aes(x = x1, y = x2, z = NULL), size = 2) + 
  scale_x_continuous(expand = c(0, 0)) + 
  scale_y_continuous(expand = c(0, 0)) + 
  scale_fill_gradientn(colors = c("yellow","red"), name = "") + 
  labs(x = expression(x[1]), y = expression(x[2]),
       title = "Mean") + 
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 10),
        plot.title = element_text(hjust = 0.5)) 

temp <- data.frame(x = xx[,1], y = xx[,2], z = sig)

ggplot(temp, aes(x, y, z = z)) + 
  geom_tile(aes(fill=z)) + 
  geom_contour(color = "black") + 
  geom_point(data = x, aes(x = x1, y = x2, z = NULL), size = 2) + 
  scale_x_continuous(expand = c(0, 0)) + 
  scale_y_continuous(expand = c(0, 0)) + 
  scale_fill_gradientn(colors = c("yellow","red"), name = "") + 
  labs(x = expression(x[1]), y = expression(x[2]),
       title = "Standard Deviation") + 
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 10),
        plot.title = element_text(hjust = 0.5)) 

temp <- data.frame(x = X[,1], y = X[,2], z = c(out))
next.point.pi <- data.frame(xx[which.max(probimpr),])
next.point.ei <- data.frame(xx[which.max(ei),])
optimum <- data.frame(x1 = mini$par[1], x2 = mini$par[2])

ggplot(temp, aes(x, y, z = z)) + 
  geom_tile(aes(fill=z)) + 
  geom_contour(color = "black") + 
  geom_point(data = next.point.pi, aes(x = Var1, y = Var2, z = NULL),
             size = 3, bg = "green", pch = 21) + 
  geom_point(data = next.point.ei, aes(x = Var1, y = Var2, z = NULL),
             size = 3, bg = "deepskyblue", pch = 21) + 
  geom_point(data = optimum, aes(x = x1, y = x2, z = NULL),
             size = 3, bg = "blue", pch = 21) +
  scale_x_continuous(expand = c(0, 0)) + 
  scale_y_continuous(expand = c(0, 0)) + 
  scale_fill_gradientn(colors = c("yellow","red"), name = "") + 
  labs(x = expression(x[1]), y = expression(x[2]),
       title = "Objective Function") + 
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 10),
        plot.title = element_text(hjust = 0.5))  


temp <- data.frame(x = xx[,1], y = xx[,2], z = probimpr)

ggplot(temp, aes(x, y, z = z)) + 
  geom_tile(aes(fill=z)) + 
  geom_contour(color = "black") + 
  geom_point(data = next.point.pi, aes(x = Var1, y = Var2, z = NULL),
             size = 3, bg = "green", pch = 21) + 
  scale_x_continuous(expand = c(0, 0)) + 
  scale_y_continuous(expand = c(0, 0)) + 
  scale_fill_gradientn(colors = c("yellow","red"), name = "") + 
  labs(x = expression(x[1]), y = expression(x[2]),
       title = expression(a[PI]*"(x)")) + 
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 10),
        plot.title = element_text(hjust = 0.5))  

temp <- data.frame(x = xx[,1], y = xx[,2], z = ei)

ggplot(temp, aes(x, y, z = z)) + 
  geom_tile(aes(fill=z)) + 
  geom_contour(color = "black") + 
  geom_point(data = next.point.ei, aes(x = Var1, y = Var2, z = NULL),
             size = 3, bg = "deepskyblue", pch = 21) +
scale_x_continuous(expand = c(0, 0)) + 
  scale_y_continuous(expand = c(0, 0)) + 
  scale_fill_gradientn(colors = c("yellow","red"), name = "") + 
  labs(x = expression(x[1]), y = expression(x[2]),
       title = expression(a[EI]*"(x)")) + 
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 10),
        plot.title = element_text(hjust = 0.5)) 


# Picture of the 2D problem

obj = function(x1,x2){
  value = -(cos((x1-.1)*x2))^2 - x1*sin(3*x1+x2)
  return(value)
}

x1 = seq(-2.25,2.5,.015)
x2 = seq(-2.5,1.75,.015)

out = outer(x1,x2,obj)
infeasible = outer(x1,x2,obj)

my.heat.colors <- function(...) heat.colors(..., rev=TRUE)

par(ps = 15, mfrow = c(1,1))
p <- persp3D(x1,x2,out,col.palette = my.heat.colors, phi = 5, theta = 20,
            box = TRUE, border = NA, shade = 0.4,zlab="f(x1,x2)",zlim=c(-3.3,2.5))
cl = contourLines(x1,x2,out)
levs = factor(sapply(cl, `[[`, "level"))
Map(
  function(cl,col) lines(trans3d(cl$x, cl$y, -3.3, pmat=p$persp), col=col),
  cl,
  heat.colors(length(levs), rev = TRUE)
)
mypoints <- trans3d(mini$par[1], mini$par[2], min(out,na.rm=TRUE), pmat=p$persp)
points(mypoints, pch=21, bg="blue")



### Figures 3.17 -- 3.19
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
EI <- (fmin - mu) * pnorm((fmin - mu) / sig) + sig * dnorm ((fmin - mu) / sig)

par( ps = 15, mfrow = c(1,2) )
plot(xx, f(xx), type = "l", xlab = "x", ylab = "Amount of Contaminants, f(x)",
     main = "Predictive Surface for n = 10")
lines(xx, mu, col = "red")
lines(xx, mu + 1.96 * sig, lty = 2, col = "red")
lines(xx, mu - 1.96 * sig, lty = 2, col = "red")
points(x, fx, pch = 21, bg = "red", cex = 1.2)

next.point <- xx[which.max(EI),1]
x <- matrix(c(x, next.point), ncol = 1)
fx <- c(fx, f(next.point)) 

points(next.point, f(next.point), pch = 21, bg = "deepskyblue", cex = 1.2)

plot(xx, EI, type = "l", xlab = "x", ylab = expression(a[EI]*"(x)"),
     main = "Expected Improvement")

for(i in 1:3){
  da <- darg(list(mle = TRUE), x)
  gp <- newGP(x, fx, da$start, g = .0001*var(fx), dK = TRUE)
  preds <- predGP(gp, xx)
  deleteGP(gp)
  
  fmin <- min(fx)
  mu <- preds$mean
  sig <- sqrt(diag(preds$Sigma))
  EI <- (fmin - mu) * pnorm((fmin - mu) / sig) + sig * dnorm ((fmin - mu) / sig)
  
  next.point <- xx[which.max(EI),1]
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

next.point <- xx[which.max(EI),1]
x <- matrix(c(x, next.point), ncol = 1)
fx <- c(fx, f(next.point)) 

points(next.point, f(next.point), pch = 21, bg = "deepskyblue", cex = 1.2)

plot(xx, EI, type = "l", xlab = "x", ylab = expression(a[EI]*"(x)"),
     main = "Expected Improvement")

for(i in 1:2){
  da <- darg(list(mle = TRUE), x)
  gp <- newGP(x, fx, da$start, g = .0001*var(fx), dK = TRUE)
  preds <- predGP(gp, xx)
  deleteGP(gp)
  
  fmin <- min(fx)
  mu <- preds$mean
  sig <- sqrt(diag(preds$Sigma))
  EI <- (fmin - mu) * pnorm((fmin - mu) / sig) + sig * dnorm ((fmin - mu) / sig)
  
  next.point <- xx[which.max(EI),1]
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

next.point <- xx[which.max(EI),1]
x <- matrix(c(x, next.point), ncol = 1)
fx <- c(fx, f(next.point)) 

points(next.point, f(next.point), pch = 21, bg = "deepskyblue")

plot(xx, EI, type = "l", xlab = "x", ylab = expression(a[EI]*"(x)"),
     main = "Expected Improvement", cex = 1.2)


### Figure 3.20
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
    da <- darg(list(mle = TRUE), x)
    gp <- newGP(x, fx, da$start, g = .0001*var(fx), dK = TRUE)
    preds <- predGP(gp, xx)
    deleteGP(gp)
    
    fmin <- min(fx)
    mu <- preds$mean
    sig <- sqrt(diag(preds$Sigma))
    EI <- (fmin - mu) * pnorm((fmin - mu) / sig) + sig * dnorm ((fmin - mu) / sig)

    
    next.point <- xx[which.max(EI),1]
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
    EI <- pnorm((fmin - mu) / sig)
    
    next.point <- xx[which.max(EI),1]
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



####################################################
### Section 3.3.3 ###
### Figures 3.21 -- 3.23

f <- function(x) -(1.4 - 3 * x) * sin(18 * x)
n <- 10
x <- lhs(n, c(0, 1.2))  
fx <- f(x)

da <- darg(list(mle = TRUE), x)
gp <- newGP(x, fx, da$start, g = .001*var(fx), dK = TRUE)
xx <- matrix(seq(0, 1.2, .001), ncol = 1)
preds <- predGP(gp, xx)
deleteGP(gp)

mu <- preds$mean
sig <- sqrt(diag(preds$Sigma))
beta <- 3
lcb <- mu - beta * sig

# Next best point to choose
next.point <- xx[which.min(lcb)]
x <- matrix(c(x, next.point), ncol = 1)
fx <- c(fx, f(next.point))

par(ps = 15)
plot(x, fx, pch = 21, bg = "red", ylim = c(-3, 3), 
     xlab = "x", ylab = "f(x)", main = "Predictive Surface")
lines(xx, mu, col = "red")
lines(xx, mu + 1.96 * sig, lty = 2, col = "red")
lines(xx, mu - 1.96 * sig, lty = 2, col = "red")
lines(xx, lcb, col = "blue")
points(next.point, f(next.point), pch = 21, bg = "blue")
legend("topleft", c("Mean", "95% PI", expression(LCB*"(x;"*beta*")"),"Data"), 
       col = c("red", "red", "blue", "red"), lty = c(1, 2, 1, NA), 
       pch = c(NA, NA, NA, 19), bty = "n")

da <- darg(list(mle = TRUE), x)
gp <- newGP(x, fx, da$start, g = .001*var(fx), dK = TRUE)
xx <- matrix(seq(0, 1.2, .001), ncol = 1)
preds <- predGP(gp, xx)
deleteGP(gp)

mu <- preds$mean
sig <- sqrt(diag(preds$Sigma))
lcb <- mu - beta * sig

next.point <- xx[which.min(lcb)]

par(ps = 15)
plot(x, fx, pch = 21, bg = "red", ylim = c(-3, 3), 
     xlab = "x", ylab = "f(x)", main = "Predictive Surface")
lines(xx, mu, col = "red")
lines(xx, mu + 1.96 * sig, lty = 2, col = "red")
lines(xx, mu - 1.96 * sig, lty = 2, col = "red")
lines(xx, lcb, col = "blue")
points(next.point, f(next.point), pch = 21, bg = "blue")
legend("topleft", c("Mean", "95% PI", expression(LCB*"(x;"*beta*")"),"Data"), 
       col = c("red", "red", "blue", "red"), lty = c(1, 2, 1, NA), 
       pch = c(NA, NA, NA, 19), bty = "n")

for(i in 1:9){
  next.point <- xx[which.min(lcb)]
  x <- matrix(c(x, next.point), ncol = 1)
  fx <- c(fx, f(next.point))
  
  da <- darg(list(mle = TRUE), x)
  gp <- newGP(x, fx, da$start, g = .001*var(fx), dK = TRUE)
  xx <- matrix(seq(0, 1.2, .001), ncol = 1)
  preds <- predGP(gp, xx)
  deleteGP(gp)
  
  mu <- preds$mean
  sig <- sqrt(diag(preds$Sigma))
  lcb <- mu - beta * sig
}

next.point <- xx[which.min(lcb)]

par(ps = 15)
plot(x, fx, pch = 21, bg = "red", ylim = c(-3, 3), 
     xlab = "x", ylab = "f(x)", main = "Predictive Surface")
lines(xx, mu, col = "red")
lines(xx, mu + 1.96 * sig, lty = 2, col = "red")
lines(xx, mu - 1.96 * sig, lty = 2, col = "red")
lines(xx, lcb, col = "blue")
points(next.point, f(next.point), pch = 21, bg = "blue")
legend("topleft", c("Mean", "95% PI", expression(LCB*"(x;"*beta*")"),"Data"), 
       col = c("red", "red", "blue", "red"), lty = c(1, 2, 1, NA), 
       pch = c(NA, NA, NA, 19), bty = "n")



### Figure 3.24
bov <- function(y, end = length(y)){ # Calculates best overall value
  prog <- rep(min(y), end)
  prog[1:min(end, length(y))] <- y[1:min(end, length(y))]
  for(i in 2:end)
    if(is.na(prog[i]) || prog[i] > prog[i-1]) prog[i] <- prog[i-1]
  return(prog)
}

f <- function(x) -(1.4 - 3 * x) * sin(18 * x)
beta <- c(0, .5, 1, 2, 4)
S <- 30
solution <- matrix(NA, nrow = 30*length(beta), ncol = 30)
iter <- 0

for(k in 1:length(beta)){
  set.seed(43)
  for(j in 1:S){
    n <- 10
    x <- lhs(n, c(0, 1.2))  
    fx <- f(x)
    for(i in 1:(ncol(solution) - n)){
      
      da <- darg(list(mle = TRUE), x)
      gp <- newGP(x, fx, da$start, g = .001*var(fx), dK = TRUE)
      xx <- matrix(seq(0, 1.2, .001), ncol = 1)
      preds <- predGP(gp, xx)
      deleteGP(gp)
      
      mu <- preds$mean
      sig <- sqrt(diag(preds$Sigma))
      lcb <- mu - beta[k] * sig
      
      next.point <- xx[which.min(lcb)]
      x <- matrix(c(x, next.point), ncol = 1)
      fx <- c(fx, f(next.point)) 
    }
    iter <- iter + 1
    solution[iter,] <- bov(fx)
    print(iter)
  } 
}

par(ps = 15)
plot(colMeans(solution[1:30,]), col = "black", lwd = 2, type = "l",
     ylim = c(-1.5, 0), xlab="black-box evaluations (n)", 
     ylab="average best objective value", axes = FALSE)
lines(colMeans(solution[30 + (1:30),]), col="green", lwd=2)
lines(colMeans(solution[2*30 + (1:30),]), col="orange", lwd=2)
lines(colMeans(solution[3*30 + (1:30),]), col="blue", lwd=2)
lines(colMeans(solution[4*30 + (1:30),]), col="red", lwd=2)
abline(h = f(.96609), lty = 2)
legend("topright", c(expression(beta*" = 0"), expression(beta*" = 0.5"),
                     expression(beta*" = 1"), expression(beta*" = 2"), 
                     expression(beta*" = 4"), "Global Solution"),
       col = c("black", " green", "orange", "blue", "red", "black"), 
       lty = c(1, 1, 1, 1, 1, 2), lwd = c(2, 2, 2, 2, 2, 1), bty = "n")
axis(1)
axis(2)


### Figures 3.25 -- 3.27
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
beta <- 3
lcb <- mu - beta * sig

par( ps = 15, mfrow = c(1,2) )
plot(xx, f(xx), type = "l", xlab = "x", ylab = "Amount of Contaminants, f(x)",
     main = "Predictive Surface for n = 10")
lines(xx, mu, col = "red")
lines(xx, mu + 1.96 * sig, lty = 2, col = "red")
lines(xx, mu - 1.96 * sig, lty = 2, col = "red")
points(x, fx, pch = 21, bg = "red", cex = 1.2)

next.point <- xx[which.min(lcb),1]
x <- matrix(c(x, next.point), ncol = 1)
fx <- c(fx, f(next.point)) 

points(next.point, f(next.point), pch = 21, bg = "deepskyblue", cex = 1.2)

plot(xx, lcb, type = "l", xlab = "x", ylab = expression(LCB*"(x;"*beta*")"),
     main = "Lower Confidence Bound")

for(i in 1:3){
  da <- darg(list(mle = TRUE), x)
  gp <- newGP(x, fx, da$start,g = .0001*var(fx), dK = TRUE)
  preds <- predGP(gp, xx)
  deleteGP(gp)
  
  fmin <- min(fx)
  mu <- preds$mean
  sig <- sqrt(diag(preds$Sigma))
  beta <- 3
  lcb <- mu - beta * sig
  
  next.point <- xx[which.min(lcb),1]
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

next.point <- xx[which.min(lcb),1]
x <- matrix(c(x, next.point), ncol = 1)
fx <- c(fx, f(next.point)) 

points(next.point, f(next.point), pch = 21, bg = "deepskyblue", cex = 1.2)

plot(xx, lcb, type = "l", xlab = "x", ylab = expression(LCB*"(x;"*beta*")"),
     main = "Lower Confidence Bound")

for(i in 1:2){
  da <- darg(list(mle = TRUE), x)
  gp <- newGP(x, fx, da$start, g = .0001*var(fx), dK = TRUE)
  preds <- predGP(gp, xx)
  deleteGP(gp)
  
  fmin <- min(fx)
  mu <- preds$mean
  sig <- sqrt(diag(preds$Sigma))
  beta <- 3
  lcb <- mu - beta * sig
  
  next.point <- xx[which.min(lcb),1]
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

next.point <- xx[which.min(lcb),1]
x <- matrix(c(x, next.point), ncol = 1)
fx <- c(fx, f(next.point)) 

points(next.point, f(next.point), pch = 21, bg = "deepskyblue")

plot(xx, lcb, type = "l", xlab = "x", ylab = expression(LCB*"(x;"*beta*")"),
     main = "Lower Confidence Bound", cex = 1.2)


### Figure 3.28
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
    da <- darg(list(mle = TRUE), x)
    gp <- newGP(x, fx, da$start, g = .0001*var(fx), dK = TRUE)
    preds <- predGP(gp, xx)
    deleteGP(gp)
    
    fmin <- min(fx)
    mu <- preds$mean
    sig <- sqrt(diag(preds$Sigma))
    beta <- 3
    lcb <- mu - beta * sig    
    
    next.point <- xx[which.min(lcb),1]
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


