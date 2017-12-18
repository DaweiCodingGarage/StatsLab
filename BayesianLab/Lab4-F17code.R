library('Flury')
data(midge)
y <- midge[midge[, "Species"] == "Af", 
           "Wing.Length"]

hist(y, xlab = "y", freq = FALSE,
     main = "Histogram of Wing Lengths")

sigma.sq.0 <- 0.01
nu.0 <- 1
mu.0 <- 1.9
kappa.0 <- 1


S <- 500

y.bar <- mean(y)
s.sq <- var(y)
n <- length(y)

nu.n <- nu.0 + n
kappa.n <- kappa.0 + n
mu.n <- (kappa.0*mu.0 + n*y.bar)/kappa.n
sigma.sq.n <- (nu.0*sigma.sq.0 + (n - 1)*s.sq + kappa.0*n*(y.bar - mu.0)^2/kappa.n)/nu.n

sigma.sq.samps <- 1/rgamma(S, nu.n/2, nu.n*sigma.sq.n/2)
theta.samps <- rnorm(S, mu.n, sqrt(sigma.sq.samps/kappa.n))

theta.samps.mcmc <- sigma.sq.samps.mcmc <- numeric(S)

theta.start <- y.bar

for (i in 1:S) {
  
  if (i == 1) {
    theta.samps.mcmc[1] <- theta.start
  } else {
    theta.samps.mcmc[i] <- rnorm(1, mu.n,
                                 sqrt(sigma.sq.samps.mcmc[i - 1]/kappa.n))
  }
  
  a.n <- (nu.0 + n + 1)/2
  b.n <- (nu.0*sigma.sq.0 + 
            sum((y - theta.samps.mcmc[i])^2) + 
            kappa.0*(theta.samps.mcmc[i] - mu.0)^2)/2
  
  sigma.sq.samps.mcmc[i] <- 1/rgamma(1, a.n, b.n)
}


par(mfrow=c(1,2))
hist(theta.samps, xlab = expression(theta), 30,
     freq = FALSE, 
     main = expression(paste("Hist. of Marginal Post. Dist. of ", 
                             theta, sep = "")))
hist(sigma.sq.samps, xlab=expression(sigma^2), 30,
     freq=FALSE,
     main = expression(paste("Hist. of Marginal Post. Dist. of ", 
                             sigma^2, sep = "")))


par(mfrow=c(1,2))
hist(theta.samps.mcmc, xlab = expression(theta),30,
     freq = FALSE, 
     main = expression(paste("Hist. of Gibbs samples of ", 
                             theta, sep = "")))
hist(sigma.sq.samps.mcmc, xlab=expression(sigma^2),30,
     freq=FALSE,
     main = expression(paste("Hist. of Gibbs samples of ", 
                             sigma^2, sep = "")))
