library('mvtnorm')

rmvnorm(n = 3, mean = rep(0, 2), sigma = diag(2))

nu0 <- 2
Sigma0 <- diag(2)
rWishart(1, df = nu0, Sigma = Sigma0)[, , 1]

data<-read.table("http://stat.duke.edu/~nln6/sta601labs/lab8data.dat",
                 header=TRUE, na.strings="NA")
data[1:3,]
data[47:50,]

#Single imputation
imp1 <- data
miss.A <- which(is.na(data$yA))
miss.B <- which(is.na(data$yB))

theta.A.hat <- mean(data$yA, na.rm = TRUE)
theta.B.hat <- mean(data$yB, na.rm = TRUE)
sigma.sq.A.hat <- var(data$yA, na.rm = TRUE)
sigma.sq.B.hat <- var(data$yB, na.rm = TRUE)
rho.hat <- cor(data$yA, data$yB, use = "complete.obs")

#impute
for (i in miss.A) {
  imp1[i, "yA"] <- theta.A.hat +
    rho.hat*sqrt(sigma.sq.A.hat/sigma.sq.B.hat)*(imp1[i, "yB"] - theta.B.hat)
}
for (i in miss.B) {
  imp1[i, "yB"] <- theta.B.hat +
    rho.hat*sqrt(sigma.sq.B.hat/sigma.sq.A.hat)*(imp1[i, "yA"] - theta.A.hat)
}

#paired t-test
t.test(imp1[, "yA"], imp1[, "yB"], paired = TRUE)

#####################
#####Gibbs sampler - multiple imputation
# Set Prior Parameters
y.bar <- apply(data, 2, mean, na.rm = TRUE) #prior mean for theta

complete <- which(complete.cases(data))
#Prior for Sigma
S <- (t(data[complete, ]) - y.bar)%*%t(t(data[complete, ]) - y.bar)/length(complete)
nu.0 <- nrow(S) + 2 #nu0
n <- nrow(data)

samps <- 5000

#Store results
y.A.samps <- y.B.samps <- matrix(nrow = samps, ncol = n)
theta.samps <- matrix(nrow = samps, ncol = 2)

# Set Starting values
Sigma <- S
theta <- y.bar
Y <- as.matrix(imp1) #Starting values are from the simple imputation method
#Can also start from sample means

# Y<-data
# for(j in 1:2)
# {
#   Y[is.na(Y[,j]),j]<-mean(Y[,j],na.rm=TRUE)
# }


for (i in 1:samps) {

  # Update theta
  y.bar.samp <- apply(Y, 2, mean)
  theta <- rmvnorm(1, y.bar.samp, Sigma/(n + 1))
  theta.samps[i, ] <- theta

  # Update Sigma
  Sn<- S + ( t(Y)-c(theta) )%*%t( t(Y)-c(theta) )
  Sigma<-solve(rWishart(1, nu.0+n, solve(Sn))[, , 1])
  rho <- cov2cor(Sigma)[1, 2]
  sigma.sq.A <- Sigma[1, 1]
  sigma.sq.B <- Sigma[2, 2]

  # Update Missing Data
  for (j in miss.A) {
    Y[j, 1] <- rnorm(1,
                     theta[1] + (rho*sqrt(sigma.sq.A/sigma.sq.B))*(Y[j, "yB"] - theta[2]),
                     sqrt(sigma.sq.A*(1 - rho)))
  }
  for (j in miss.B) {
    Y[j, 2] <- rnorm(1,
                     theta[2] + (rho*sqrt(sigma.sq.B/sigma.sq.A))*(Y[j, "yA"] - theta[1]),
                     sqrt(sigma.sq.B*(1 - rho)))
  }
  y.A.samps[i, ] <- Y[, 1]
  y.B.samps[i, ] <- Y[, 2]
}

#Mean and quantiles
mean(theta.samps[, 1] - theta.samps[, 2])
quantile(theta.samps[,1 ] - theta.samps[, 2], c(0.025, 0.975))

#Compare with t-test
t.test(imp1[, "yA"], imp1[, "yB"], paired = TRUE)

imp2 <- cbind(apply(y.A.samps, 2, mean),
              apply(y.B.samps, 2, mean))

plot(imp1[31:50, 1], imp1[31:50, 2],ylim=c(6,18), pch=16,xlab = expression(y[A]),
     ylab = expression(y[B]),
     main = "Comparison of Imputed Values")
points(imp2[31:50, 1], imp2[31:50, 2], pch=17, col = "red")
abline(a = 0, b = 1)
legend("bottomright", pch = c(16, 17), col = c("black", "red"), legend = c("Simple", "Multiple"))
