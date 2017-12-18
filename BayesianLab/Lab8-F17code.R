# Load libraries and read in data
set.seed(123)
library(mvtnorm)
library(coda)

data <- read.table("http://www.stat.washington.edu/~pdhoff/Book/Data/hwdata/swim.dat",
                   header = FALSE)

# Set up response variable and design matrix
y <- as.matrix(data)
X <- cbind(1, seq(0, 10, by = 2))


# Set priors for Beta -- Note should probably have less variability then we are giving
beta.0 <- c(23, 0)
Sigma.0 <- matrix(c(0.25, 0, 0, 0.01), 
                  nrow = 2, ncol = 2)
as.numeric(pmvnorm(lower = c(22, -0.2),
                   upper = c(24, 0.2), 
                   mean = beta.0, sigma = Sigma.0))

# Set prior on sigma.sq 
sigma.sq.0 <- 0.002
nu.0 <- 4
plot(density(sqrt(1/rgamma(10000, shape=nu.0/2, rate=nu.0*sigma.sq.0/2))),main="")


# Function to sample full conditional for beta
samp.beta <- function(y, X, sigma.sq, beta.0, Sigma.0) {
  beta.v <- solve(solve(Sigma.0) + t(X)%*%X/sigma.sq)
  beta.e <- beta.v%*%(solve(Sigma.0)%*%beta.0 + 
                        t(X)%*%y/sigma.sq)
  return(t(rmvnorm(1, mean = beta.e, sigma = beta.v)))
}

# Function to sample full conditional for sigma.sq
samp.sigma.sq <- function(y, X, beta, 
                          nu.0, sigma.sq.0) {
  sigma.sq.inv.a <- (nu.0 + length(y))/2
  sigma.sq.inv.b <- (nu.0*sigma.sq.0 + sum(y - X%*%beta)^2)/2
  return(1/rgamma(1, sigma.sq.inv.a, sigma.sq.inv.b))
}

# Function to sample from data Y given X, beta, sigma.sq
samp.y <- function(X, beta, sigma.sq) {
  return(rmvnorm(1, 
                 mean = X%*%beta, 
                 sigma = sigma.sq*diag(nrow(X))))
}

# Gibb's sampler
Xpred <- matrix(c(1,12),1,2)  # predict at week 12

S <- 1000
betas <- array(NA, c(4, S, 2))
sigma.sqs <- matrix(NA, 4, S)
y.fits <- array(NA, c(4, S, 6))
y.preds <- matrix(NA, 4, S)

# Loop for each swimmer
for(i in 1:4){
  y.swimmer <- y[i,]
  # Draws for each posterior
  for (k in 1:S) {
    if(k==1){
      beta <- c(23,0)
      sigma.sq <- 0.1
    }
    # Draw beta
    betas[i, k, ] <- beta <-
      samp.beta(y.swimmer, X, sigma.sq, beta.0, Sigma.0)
    # Draw sigma.sq
    sigma.sqs[i, k] <- sigma.sq <-
      samp.sigma.sq(y.swimmer, X, beta, nu.0, sigma.sq.0)
    # Find mean and find extrapolated mean
    y.fits[i, k, ] <- X%*%beta
    y.preds[i, k] <- samp.y(Xpred, beta, sigma.sq)
  }
}

# Diagnostics
apply(betas, 1, effectiveSize)
apply(sigma.sqs, 1, effectiveSize)

par(mfrow = c(1, 3))
plot(betas[1, , 1], xlab = expression(beta[1]), type="l")
plot(betas[1, , 2], xlab = expression(beta[2]), type="l")
plot(sigma.sqs[1, ], xlab = expression(sigma^2), type="l")

# Posterior Means for each swimmer
par(mfrow = c(1, 1))
matplot(X[,2], t(y), xlab = "Weeks", ylab = "Time", ylim=c(22.4, 23.8))
matpoints(X[,2], t(apply(y.fits,c(1,3),mean)), type="l", lty=1, lwd=2)
matpoints(X[,2], t(apply(y.fits, c(1,3), quantile, probs=0.025)), type="l", lty = 2)
matpoints(X[,2], t(apply(y.fits, c(1,3), quantile, probs=0.975)), type="l", lty = 2)

# Which swimmer is fastest in week 12 (extrapolation)
boxplot(t(y.preds), xlab = "Swimmer", ylab = "Predicted Time")
table(apply(y.preds,2,which.min))/S
