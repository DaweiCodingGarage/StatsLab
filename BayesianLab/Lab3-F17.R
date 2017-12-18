y <- 36; n <- 20
a <- 1.25; b <- 0.5
a.y <- a + y; b.y <- b + n
Theta <- seq(0.1, 5, by = 0.1)
ptheta.y <- dgamma(Theta, a.y, b.y)
plot(Theta, ptheta.y, type = "l", main="Posterior")

S <- 10000
stheta.y.mc10000 <- sy.y.mc10000 <- numeric(S)
for (s in 1:S) {
  stheta.y.mc10000[s] <- rgamma(1, a.y, b.y)
  sy.y.mc10000[s] <- rpois(1, 20*stheta.y.mc10000[s])
}

hist(sy.y.mc10000, main = "Hist. of Post. Pred. Samples")
abline(v = mean(sy.y.mc10000), col = "blue")
abline(v = y, col = "red")

nesteggs <- read.table(
  "http://www2.stat.duke.edu/courses/Fall16/sta601.001//files/nesteggs.dat",
  header=TRUE)
y.i <- nesteggs$Eggs

t <- mean(y.i >= 5)

S <- 10000
stheta.y.mc10000 <- st.y.mc10000 <- numeric(S)
for (s in 1:S) {
  stheta.y.mc10000[s] <- rgamma(1, a.y, b.y)
  sy.y <- rpois(20, stheta.y.mc10000[s])
  st.y.mc10000[s] <- mean(sy.y >= 5)
}

hist(st.y.mc10000, main = "Hist. of Post. Pred. Samples")
abline(v = mean(st.y.mc10000), col = "blue")
abline(v = t, col = "red")

hist(y.i, main = "Hist. of Observed Data")
