library(Flury)
data(midge)
y <- midge[midge$Species == "Af", "Wing.Length"]

# Summarize data
n <- length(y)
y.bar <- mean(y)

# Set seed for reproducible approximations
set.seed(1)

```{r}
# Parameters for semi-congugate priors
mu.0 <- 1.9; tau.sq.0 <- 0.95^2;
nu.0 <- 1; sigma.sq.0 <- 0.01;

# Initialize MCMC structures
S <- 1000
theta.samps <- sigma.sq.samps <- numeric(S)

# Perform Gibbs sampling
for (i in 1:S) {
  # Update theta
  if (i == 1) {
    theta.samps[1] <- y.bar
  } else {
    d.n <- (1/tau.sq.0 + n/sigma.sq.samps[i-1])
    c.n <- mu.0/tau.sq.0 + n*y.bar/sigma.sq.samps[i-1]
    theta.samps[i] <- rnorm(1, c.n/d.n, sqrt(1/d.n))
  }

  # Update sigma.sq
  a.n <- 0.5*(nu.0 + n)
  b.n <- 0.5*(nu.0*sigma.sq.0 + sum((y - theta.samps[i])^2))
  sigma.sq.samps[i] <- 1/rgamma(1, a.n, b.n)
}

par(mfrow=c(1,2))
plot(theta.samps, type="l")
plot(sigma.sq.samps, type="l")

par(mfrow=c(1,2))
acf(theta.samps)
acf(sigma.sq.samps)


library(coda)
c(effectiveSize(mcmc(theta.samps)),
  effectiveSize(mcmc(sigma.sq.samps)))

hist(theta.samps, freq = FALSE, ylim = c(0, 9),
     xlab = expression(theta),
     main = expression(paste("Approx. of p(",theta,"|y)", sep = "")))
lines(density(theta.samps), col = "blue")

########################## Latent Class membership
data(midge)
y <- midge[, "Wing.Length"]
n <- length(y)

# Prior parameters
a <- 3; b <-3;
S <- 1000
X <- matrix(NA,n,S)
p <- numeric(S)
# Random initialization to groups
X[,1] <- sample(1:0, n, prob=c(0.5,0.5), replace=T)

for (i in 2:S){
  # Update for p
  n.1 <- sum(X[,i-1]==1)
  n.0 <- sum(X[,i-1]==0)
  p[i] <- rbeta(1, a + n.1, b + n.0)

  # Update for X's (could be vectorized)
  for (j in 1:n){
    w1 <- p[i]*dnorm(y[j], 1.8, 0.13)
    w0 <- (1-p[i])*dnorm(y[j], 1.9, 0.08)
    w  <- w1/(w1+w0)
    X[j,i] <- sample(1:0, 1, prob=c(w, 1-w))
  }
}



