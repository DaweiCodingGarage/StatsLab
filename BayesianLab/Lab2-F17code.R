y <- 9; n <- 5 
a <- 1.25; b <- 0.5
a.y <- a + y; b.y <- b + n

#Gamma Prior
Theta <- seq(0.1, 5, by = 0.1)
ptheta.y <- dgamma(Theta, a.y, b.y)
plot(Theta, ptheta.y, type = "l")

citheta.y <- qgamma(c(0.025, 0.975), 
                      a.y, b.y)
citheta.y

#95% CI
plot(Theta, ptheta.y, type = "l")
abline(v = citheta.y, lty = 2)

#========Find HPD Interval by discretizing
# Discretize theta
Theta <- seq(0.01, 5, by = 0.01) 
ptheta.y <- dgamma(Theta, a.y, b.y)

# Start with a horizontal line on top of posterior
hpd.cutoff <- max(ptheta.y)

# Find intersection of first line w/posterior (will be one point)
hptheta.y <- range(Theta[ptheta.y >= hpd.cutoff])
# Find area between line and posterior (will be 0)
hpd.p <- pgamma(hptheta.y[2], a.y, b.y) - pgamma(hptheta.y[1], a.y, b.y)
while(hpd.p <= 0.95 & hpd.cutoff > 0) {
  hpd.cutoff <- hpd.cutoff - 0.005 # Move hline down
  # Find intersection of line w/posterior
  hptheta.y <- range(Theta[ptheta.y >= hpd.cutoff])
  # Compute area between line and posterior
  hpd.p <- pgamma(hptheta.y[2], a.y, b.y) - pgamma(hptheta.y[1], a.y, b.y)
}
hptheta.y

#Plot 95% CI and 95% HPD Interval
plot(Theta, ptheta.y, type = "l",
     xlab = "theta", ylab = "p(theta | y)")
abline(h = hpd.cutoff, col = "red")
abline(v = citheta.y, lty = 2)
abline(v = hptheta.y, lty = 2, col = "red")
legend("topright", lty = c(2, 2), col = c("red", "black"),
       legend = c("95% HPDI", "95% CI"))

#======================SENSITIVITY ANALYSIS

#Range of parameters
a <- seq(0.1, 2, by = 0.1)
b <- seq(0.1, 1, by = 0.1)

#Posterior mean
etheta.y <- matrix(nrow = length(a), ncol = length(b))
for (i in 1:length(a)) {
  for (j in 1:length(b)) {
    etheta.y[i, j] <- (y + a[i])/(n + b[j])
  }
}

contour(a, b, etheta.y, 
        levels = seq(1, 2, by = 0.05), 
        xlab = "a", ylab = "b")

#Consider a different sample size
y.1 <- 36
n.1 <- 20
etheta.y.1 <- matrix(nrow = length(a), ncol = length(b))
for (i in 1:length(a)) {
  for (j in 1:length(b)) {
    etheta.y.1[i, j] <- (y.1 + a[i])/(n.1 + b[j])
    etheta.y[i, j] <- (y + a[i])/(n + b[j])
  }
}

#Compare two scenarios' contour plots
par(mfrow = c(1, 2))
contour(a, b, etheta.y, 
        xlab = "a", ylab = "b",
        main = "n = 5")
contour(a, b, etheta.y.1,
        xlab = "a", ylab = "b",
        main = "n = 20")



#========MC simulations of posterior distribution
set.seed(1)

y.1 <- 36
n.1 <- 20

# Discrete prior 
Theta <- seq(0.1, 5, by = 0.1)
py.theta.1 <- dpois(y.1, n.1*Theta)
ptheta <- rep(1/length(Theta), length(Theta)) #prior
pytheta.1 <- py.theta.1*ptheta    
ptheta.y.1 <- pytheta.1/sum(pytheta.1)  #posterior

#Sample 100 draws
stheta.y.1.mc100 <- sample(x = Theta, 
                           size = 100, 
                           replace = TRUE, 
                           prob = ptheta.y.1)

#Plot histogram of 100 draws
par(mfrow=c(1,1))
hist(stheta.y.1.mc100, freq = FALSE, main = "Histogram of 100 Draws")


#Sample 50k draws
stheta.y.1.mc50000 <- sample(x = Theta, 
                             size = 50000, 
                             replace = TRUE, 
                             prob = ptheta.y.1)


hist(stheta.y.1.mc50000, freq = TRUE, main = "Histogram of 50,000 Draws")

mean(stheta.y.1.mc100)
mean(stheta.y.1.mc50000)

var(stheta.y.1.mc100)
var(stheta.y.1.mc50000)

quantile(stheta.y.1.mc100, c(0.025, 0.975))
quantile(stheta.y.1.mc50000, c(0.025, 0.975))

n.2 <- 40
y.2 <- 72

py.theta.2 <- dpois(y.2, n.2*Theta)
ptheta <- rep(1/length(Theta), length(Theta))
pytheta.2 <- py.theta.2*ptheta
ptheta.y.2 <- pytheta.2/sum(pytheta.2)

stheta.y.2.mc50000 <- sample(x = Theta, 
                             size = 50000, 
                             replace = TRUE, 
                             prob = ptheta.y.2)

mean(stheta.y.2.mc50000 <= stheta.y.1.mc50000)

hist(stheta.y.2.mc50000/stheta.y.1.mc50000, main = "Hist. of theta2 / theta1")
abline(v = mean(stheta.y.2.mc50000/stheta.y.1.mc50000), col = "red")

