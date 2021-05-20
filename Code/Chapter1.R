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
### Chapter 1 Code


### Section 1.2 ###
### Figure 1.7
n <- 10000
X <- lhs(n, rbind(c(0, 90), c(0, 90), c(2e-06, 4e-06), c(0.1, 0.2),
                  c(0.01, 0.02), c(0.01, 0.02), c(1, 2), c(5, 10)))

consumption <- rep(NA, n)
speed <- rep(NA, n)
range <- rep(NA, n)
for(i in 1:n){
  out <- sprinkler(X[i,1], X[i,2], X[i,3], X[i,4],
                   X[i,5], X[i,6], X[i,7], X[i,8])
  
  consumption[i] <- out$obj[1]
  speed[i] <- out$obj[2]
  range[i] <- out$obj[3]
  
}
scatterplot3d(speed, range, consumption, zlim = c(0, 12), 
              pch = ".", cex.symbols = 2,
              xlab = "Speed", ylab = "Range", zlab = "Consumption")


### Figure 1.8
`%notin%` <- Negate(`%in%`)
Z <- NULL

lower <- c(0, 0, 2e-6, 0.1, 0.01, 0.01, 1, 5)
upper <- c(90, 90, 4e-6, 0.2, 0.02, 0.02, 2, 10)
midpoint <- (upper + lower) / 2

for(i in 1:8){
  for(j in 1:8){
    
    Y <- expand.grid(seq(lower[i], upper[i], len = 100),
                     seq(lower[j], upper[j], len = 100))
    X <- matrix(NA, nrow = nrow(Y), ncol = 8)
    X[,(1:8)[(1:8) %notin% c(i,j)]] <- rep(midpoint[-c(i,j)], each = nrow(Y))
    X[, i] = Y[,1]
    X[, j] = Y[,2]
    
    n <- nrow(X)
    range <- rep(NA, n)
    for(k in 1:n){
      out <- sprinkler(X[k,1], X[k,2], X[k,3], X[k,4],
                       X[k,5], X[k,6], X[k,7], X[k,8])
      range[k] <- out$obj[3]
    } # End of k
    
    if(i == j){
      temp <- data.frame(x = X[,i], y = X[,j], z = NA, type = paste(i,j))
      Z <- rbind(Z, temp)
    }else{
      temp <- data.frame(x = X[,i], y = X[,j], z = range, type = paste(i,j))
      Z <- rbind(Z, temp)
    } # end of else
    
  } # End of j
} # End of i

k <- 1
plots <- list()
for(i in 1:8){
  for(j in 1:8){
    temp <- Z[Z$type == paste(j, i, sep = " "),]
    p <- ggplot(temp, aes(x, y, z = z)) + 
      geom_tile(aes(fill=z)) + 
      geom_contour(color = "black") + 
      scale_x_continuous(expand = c(0, 0)) + 
      scale_y_continuous(expand = c(0, 0)) + 
      scale_fill_gradientn(colors = c("yellow", "red"), name = "") + 
      labs(x = "", y = "") + 
      theme(axis.title.x=element_blank(),
            axis.text.x=element_blank(),
            axis.ticks.x=element_blank(),
            axis.title.y=element_blank(),
            axis.text.y=element_blank(),
            axis.ticks.y=element_blank()) + 
      guides(fill=FALSE)
    plots[[k]] <- p
    k <- k + 1
  }
}

out <- do.call(grid.arrange,plots)


####################################################
### Section 1.3 ###
### Figure 1.9
n <- 4
x1 <- (c(1, 4, 3, 2) - .5) / n
x2 <- (c(3, 1, 2, 4) - .5) / n
X  <- cbind(x1, x2)

par(ps = 15)
plot(X, xlab = expression(x[1]), ylab = expression(x[2]), xlim = c(0, 1),
     main = "Latin Hypercube Design", pch = "X", ylim = c(0, 1), cex = 3)
abline(v = seq(0, 1, length = (n + 1)), lty = 2)
abline(h = seq(0, 1, length = (n + 1)), lty = 2)

x1 <- (c(4, 3, 2, 1) - .5)/n
x2 <- (c(1, 2, 3, 4) - .5)/n
X  <- cbind(x1, x2)

par(ps = 15)
plot(X, xlab = expression(x[1]), ylab = expression(x[2]), xlim = c(0, 1),
     main = "Latin Hypercube Design", pch = "X", ylim = c(0, 1), cex = 3)
abline(v = seq(0, 1, length = (n + 1)), lty = 2)
abline(h = seq(0, 1, length = (n + 1)), lty = 2)


### Figure 1.10
n <- 10
X <- lhs(n, rbind(c(0, 1), c(0, 1)))
par(ps = 15)
plot(X, xlab = expression(x[1]), ylab = expression(x[2]),
     main = "Latin Hypercube Sample", pch = 21, bg = "red")
abline(v = seq(0, 1, length = (n + 1)), lty = 2)
abline(h = seq(0, 1, length = (n + 1)), lty = 2)


### Figure 1.11
n <- 10
x1 <- runif(n)
x2 <- runif(n)
par(ps=15)
plot(x1, x2, xlab = expression(x[1]), ylab = expression(x[2]),
     main = "Simple Random Sample", pch = 21, bg = "red")
