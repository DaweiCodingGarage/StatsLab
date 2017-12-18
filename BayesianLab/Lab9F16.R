library('coda')
numit<-11000      #burn in 1000
#vector to store results
theta.start<-0.5;  #initial value for theta

dist<-function(x){
  y<-(exp(-0.5*(x^2))+0.5*exp(-0.5*((x-3)^2)))
  return(y);
}
msamp<-function(numit,theta.start,sigcand){
  THETA<-rep(0,numit);
  accept<-0;        #acceptance rate
  theta.now<-theta.start;
  for(i in 1: numit){
    theta.new<-rnorm(1,theta.now, sigcand)
    priornew<-log(dist(theta.new))
    priornow<-log(dist(theta.now))
    lograte<- priornew - priornow
    logu<-log(runif(1))
    
    if(lograte > logu){
      theta.now<-theta.new
      accept<-accept+1
    }
    THETA[i]<-theta.now
  }
  accrate = accept/numit; #Rate of acceptance
  return(list(THETA,accrate))
}


msample<-msamp(numit,theta.start,sigcand=4)
THETA<-msample[[1]]
accrate<-msample[[2]]
accrate


#Plots
xtest<-seq(-3,8,length=1000);
thetaanat<-(2/(3*sqrt(2*pi)))*(exp(-0.5*(xtest^2))+0.5*exp(-0.5*((xtest-3)^2)));
thetasam<-density(THETA[1001:11000])

par(mfrow=c(1,2))
plot(xtest,thetaanat, col="red",type="l",main="True vs Sample",ylim=c(-0.01,0.35) )
lines(thetasam,col="blue")
legend('topright',c("Analytical","Sample"),col=c('red','blue'),lty=c(1,1))

traceplot(mcmc(THETA[1001:2000]))


#Change sigma_c
msample1<-msamp(numit,theta.start,sigcand=0.05)
THETA1<-msample1[[1]]
accrate1<-msample1[[2]]
accrate1

#Change sigma_c
msample2<-msamp(numit,theta.start,sigcand=100)
THETA2<-msample2[[1]]
accrate2<-msample2[[2]]
accrate2

#ACF
par(mfrow=c(1,3))
acf(ts(THETA),main=expression(paste('THETA,',sigma[c],"=4")))
acf(ts(THETA1),main=expression(paste('THETA1,',sigma[c],"=0.05")))
acf(ts(THETA2),main=expression(paste('THETA2,',sigma[c],"=100")))
