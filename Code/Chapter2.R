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
### Chapter 2 Code


### Section 2.1 ###
### Figure 2.1
x <- seq(0, 5, .01)
tau <- as.matrix(dist(x, diag = T, upper = T))
color <- c("black","red","blue","green")

par( ps = 15)
plot(x, x, type = 'n', ylim = c(-3, 3), xlim = c(0, 5), 
     ylab = "Y(x)", main = "")
R <- exp((-abs(tau)^2)) # Squared Exponential

for(i in 1:4){
  y <- rmvnorm(1, rep(0, length(x)), R, method = "svd")
  lines(x, y, col = color[i])
}

plot(x, x, type = 'n', ylim = c(-3, 3), xlim = c(0, 5),
     ylab = "Y(x)", main = "")
R <- exp((-abs(tau))) # Powered Exponential

for(i in 1:4){
  y <- rmvnorm(1, rep(0, length(x)), R, method = "svd")
  lines(x, y, col = color[i])
}

plot(x, x, type = 'n', ylim = c(-3, 3), xlim = c(0, 5),
     ylab = "Y(x)", main = "")
R <- exp(tau * besselJ(tau, 1)) # Matern
diag(R) <- 1

for(i in 1:4){
  y <- rmvnorm(1, rep(0, length(x)), R, method = "svd")
  lines(x, y, col = color[i])
}


### Figure 2.2
f <- function(x) exp(-1.4 * x) * cos(7 * pi * x / 2)
n <- 7
x <- lhs(n, c(0, 1))  
fx <- f(x)

gp <- newGP(x, fx, 0.015, 1e-6, dK = TRUE)
xx <- matrix(seq(0, 1, .001), ncol = 1)
preds <- predGP(gp, xx)
deleteGP(gp)

N <- 100
y <- rmvt(N, preds$Sigma, preds$df)
y <- y + t(matrix(rep(preds$mean, N), ncol = N))

par( ps = 15 )
matplot(xx, t(y), col = "gray", lwd = 0.5, lty = 1, type = "l",
        xlab = "x", ylab = "f(x)")
lines(xx, f(xx), lwd = 2)
lines(xx, preds$mean, col = "red", lwd = 2, lty = 2)
points(x, fx, pch = 21, bg = "red")
legend("topright", c("f(x)","Mean", "GP Realizations", "Data"),
       lty = c(1, 2, 1, NA), pch = c(NA, NA, NA, 19), 
       col = c("black", "red", "gray", "red"), bty = "n")


### Figure 2.3
f <- function(x1, x2) -(cos((x1 - .1) * x2))^2 - x1 * sin(3 * x1 + x2)
n <- 20
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

out <- outer(x1, x2, f)
par(ps=15)
image(x1, x2, out, xlab = expression(x[1]), ylab = expression(x[2]),
      main = "True Function")
points(x[,1], x[,2])
par(ps=15)
image(x1, x2, matrix(preds$mean, nrow = length(x1)), xlab = expression(x[1]), 
      ylab = expression(x[2]), main = "Mean")
par(ps=15)
image(x1, x2, matrix(sqrt(diag(preds$Sigma)), nrow = length(x1)),
      xlab = expression(x[1]), ylab = expression(x[2]),
      main = "Standard Deviation")
points(x[,1], x[,2])


### Figure 2.4
f <- function(x) exp(-1.4 * x) * cos(7 * pi * x / 2)
n <- 7
x <- lhs(n, c(0, 1))  
fx <- f(x)

gp <- newGP(x, fx, d = 0.05, g = 1e-6, dK = TRUE)
xx <- matrix(seq(0, 1, .001), ncol = 1)
preds <- predGP(gp, xx)
deleteGP(gp)

N <- 100
y <- rmvnorm(N, preds$mean, preds$Sigma)

par( ps = 15 )
matplot(xx, t(y), col = "gray", lwd = 0.5, lty = 1, type = "l",
        xlab = "x", ylab = "f(x)")
lines(xx, f(xx), lwd = 2)
lines(xx, preds$mean, col = "red", lwd = 2, lty = 2)
points(x, fx, pch = 21, bg = "red")
legend("topright", c("f(x)","Mean", "GP Realizations", "Data"),
       lty = c(1, 2, 1, NA), pch = c(NA, NA, NA, 19), 
       col = c("black", "red", "gray", "red"), bty = "n")


### Section 2.2 ###



### Section 2.3 ###
### Figure 2.8
gauss <- function(r, c) exp(-(c * r)^2)
multiquad <- function(r, c) sqrt(r^2 + c^2)
iquad <- function(r, c) 1 / (r^2 + c^2)
imultiquad <- function(r, c) 1 / sqrt(r^2 + c^2)
polyharm <- function(r, c) ifelse(rep(c, length(r)) %% 2 != 0, r^c, r^c * log(abs(r)))

r <- seq(-6, 6, .01)
plot(r, gauss(r, c = 0.35), type = "l", ylab = expression(varphi*"(r)"), 
     xlab = "r", col = "blue")
lines(r, gauss(r, c = 0.5), col = "green4")
lines(r, gauss(r, c = 1), col = "red")
lines(r, gauss(r, c = 3.5), col = "deepskyblue")
lines(r, gauss(r, c = 10), col = "purple")
legend("topright", c("c = 0.35", "c = 0.5", "c = 1", "c = 3.5", "c = 10"), 
       lty = 1, col = c("blue", "green4", "red", "deepskyblue", "purple"),
       bty = "n")

plot(r, multiquad(r, c = 0.25), type = "l", ylab = expression(varphi*"(r)"), 
     xlab = "r", col = "blue", ylim = c(0, 8))
lines(r, multiquad(r, c = 1), col = "green4")
lines(r, multiquad(r, c = 2.5), col = "red")
lines(r, multiquad(r, c = 4), col = "deepskyblue")
lines(r, multiquad(r, c = 5.5), col = "purple")
legend("bottomright", c("c = 0.25", "c = 1", "c = 2.5", "c = 4", "c = 5.5"), 
       lty = 1, col = c("blue", "green4", "red", "deepskyblue", "purple"),
       bty = "n")

plot(r, iquad(r, c = 0.6), type = "l", ylab = expression(varphi*"(r)"), 
     xlab = "r", col = "blue", ylim = c(0, 3))
lines(r, iquad(r, c = 0.8), col = "green4")
lines(r, iquad(r, c = 1), col = "red")
lines(r, iquad(r, c = 1.5), col = "deepskyblue")
lines(r, iquad(r, c = 2), col = "purple")
legend("topright", c("c = 0.6", "c = 0.8", "c = 1", "c = 1.5", "c = 2"), 
       lty = 1, col = c("blue", "green4", "red", "deepskyblue", "purple"),
       bty = "n")

r <- seq(-2, 2, .01)
plot(r, polyharm(r, c = 1), type = "l", ylab = expression(varphi*"(r)"), 
     xlab = "r", col = "blue", ylim = c(-4, 4))
lines(r, polyharm(r, c = 2), col = "green4")
lines(r, polyharm(r, c = 3), col = "red")
lines(r, polyharm(r, c = 4), col = "deepskyblue")
lines(r, polyharm(r, c = 5), col = "purple")
legend("bottomright", c("c = 1", "c = 2", "c = 3", "c = 4", "c = 5"), 
       lty = 1, col = c("blue", "green4", "red", "deepskyblue", "purple"),
       bty = "n")


### Figure 2.9
set.seed(4)
f <- function(x) exp(-1.4 * x) * cos(7 * pi * x / 2)
n <- 5
x <- lhs(n, c(0, 1))  
fx <- f(x)

# Gaussian RBF
c <- .1
r <- as.matrix(dist(x, diag = T, upper = T))
R <- exp(-(c*r)^2) 
lambda <- solve(R) %*% fx

xx <- seq(0, 1, .01)
r <- as.matrix(dist(c(x, xx), diag = T, upper = T))[1:n, (n + 1):(n + length(xx))]
s <- c(t(lambda) %*% exp(-(c*r)^2))

c <- 1
r <- as.matrix(dist(x, diag = T, upper = T))
R <- exp(-(c*r)^2) 
lambda <- solve(R) %*% fx

r <- as.matrix(dist(c(x, xx), diag = T, upper = T))[1:n, (n + 1):(n + length(xx))]
s2 <- c(t(lambda) %*% exp(-(c*r)^2))

c <- 10
r <- as.matrix(dist(x, diag = T, upper = T))
R <- exp(-(c*r)^2) 
lambda <- solve(R) %*% fx

r <- as.matrix(dist(c(x, xx), diag = T, upper = T))[1:n, (n + 1):(n + length(xx))]
s3 <- c(t(lambda) %*% exp(-(c*r)^2))


par( ps = 15)
plot(xx, f(xx), type = "l", xlab = "x", ylab = "f(x)")
points(x, fx, pch = 21, bg = "red")
lines(xx, s, col = "red", lty = 2)
lines(xx, s2, col ="blue", lty = 2)
lines(xx, s3, col ="green4", lty = 2)
legend("topright", c("f(x)", "c = 0.1", "c = 1",
                     "c = 10", "Data"),
       col = c("black", "red", "blue", "green4", "red"), 
       lty = c(1, 2, 2, 2, NA), pch = c(NA, NA, NA, NA, 19), bty = "n")


# Multiquadratic RBF
c <- .1
r <- as.matrix(dist(x, diag = T, upper = T))
R <- sqrt(r^2 + c^2)
lambda <- solve(R) %*% fx

r <- as.matrix(dist(c(x, xx), diag = T, upper = T))[1:n, (n + 1):(n + length(xx))]
s <- c(t(lambda) %*% sqrt(r^2 + c^2))

c <- 1
r <- as.matrix(dist(x, diag = T, upper = T))
R <- sqrt(r^2 + c^2)
lambda <- solve(R) %*% fx

r <- as.matrix(dist(c(x, xx), diag = T, upper = T))[1:n, (n + 1):(n + length(xx))]
s2 <- c(t(lambda) %*% sqrt(r^2 + c^2))

c <- 10
r <- as.matrix(dist(x, diag = T, upper = T))
R <- sqrt(r^2 + c^2)
lambda <- solve(R) %*% fx

r <- as.matrix(dist(c(x, xx), diag = T, upper = T))[1:n, (n + 1):(n + length(xx))]
s3 <- c(t(lambda) %*% sqrt(r^2 + c^2))

par( ps = 15)
plot(xx, f(xx), type = "l", xlab = "x", ylab = "f(x)")
points(x, fx, pch = 21, bg = "red")
lines(xx, s, col = "red", lty = 2)
lines(xx, s2, col ="blue", lty = 2)
lines(xx, s3, col ="green4", lty = 2)
legend("topright", c("f(x)", "c = 0.1", "c = 1",
                     "c = 10", "Data"),
       col = c("black", "red", "blue", "green4", "red"), 
       lty = c(1, 2, 2, 2, NA), pch = c(NA, NA, NA, NA, 19), bty = "n")


# Polyharmonic Spline RBF
c <- 1
r <- as.matrix(dist(x, diag = T, upper = T))
R <- r^c 
lambda <- solve(R) %*% fx

r <- as.matrix(dist(c(x, xx), diag = T, upper = T))[1:n, (n + 1):(n + length(xx))]
s <- c(t(lambda) %*% r^c )

c <- 3
r <- as.matrix(dist(x, diag = T, upper = T))
R <- r^c 
lambda <- solve(R) %*% fx

r <- as.matrix(dist(c(x, xx), diag = T, upper = T))[1:n, (n + 1):(n + length(xx))]
s2 <- c(t(lambda) %*% r^c )

c <- 5
r <- as.matrix(dist(x, diag = T, upper = T))
R <- r^c 
lambda <- solve(R) %*% fx

r <- as.matrix(dist(c(x, xx), diag = T, upper = T))[1:n, (n + 1):(n + length(xx))]
s3 <- c(t(lambda) %*% r^c )


par( ps = 15)
plot(xx, f(xx), type = "l", xlab = "x", ylab = "f(x)")
points(x, fx, pch = 21, bg = "red")
lines(xx, s, col = "red", lty = 2)
lines(xx, s2, col ="blue", lty = 2)
lines(xx, s3, col ="green4", lty = 2)
legend("topright", c("f(x)", "c = 1", "c = 3",
                     "c = 5", "Data"),
       col = c("black", "red", "blue", "green4", "red"), 
       lty = c(1, 2, 2, 2, NA), pch = c(NA, NA, NA, NA, 19), bg = "white", box.col = "white")
