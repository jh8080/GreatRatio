
rm(list = ls())
library(urca)     # Install urca package
n=100             # Sample size
nit=10000         # Number of Monte Carlo Trials

# Start simulation
set.seed(1234)
stat=matrix(NA,nrow=nit,ncol=3)
for (i in 1:nit){
  y=cumsum(rnorm(n))       # Generate a random walk
  # Test with no intercept, tau
  M=ur.df(y, type = "none", lags = 0)
  stat[i,1]=M@teststat[1]
  # Test with  intercept, tau_mu
  M=ur.df(y, type = "drift", lags = 0)
  stat[i,2]=M@teststat[1]
  # Test with intercept and trend,  tau_tau
  M=ur.df(y, type = "trend", lags = 0)
  stat[i,3]=M@teststat[1]
}
# End of Simulation 


#Plots
par(mfrow=c(1,3))
plot(density(stat[,1]),col="red",xlim=c(-5,3),lwd=2,
     main="tau (red) vs. N(0,1)")
curve(dnorm(x, mean=0, sd=1), 
      col="black", lwd=2, add=TRUE, yaxt="n")
plot(density(stat[,2]),col="red",xlim=c(-5,3),lwd=2,
     main="tau_mu (red) vs. N(0,1)")
curve(dnorm(x, mean=0, sd=1), 
      col="black", lwd=2, add=TRUE, yaxt="n")
plot(density(stat[,3]),col="red",xlim=c(-5,3),lwd=2,
     main="tau_tau (red) vs. N(0,1)")
curve(dnorm(x, mean=0, sd=1), 
      col="black", lwd=2, add=TRUE, yaxt="n")


# Calculate the critical values
prob=c(0.01,0.05, 0.10, 0.15,0.20,0.25,0.30,0.35,0.40,0.45,0.50)
CRV1=round(quantile(stat[,1],probs = prob),3)
CRV2=round(quantile(stat[,2],probs = prob),3)
CRV3=round(quantile(stat[,3],probs = prob),3)

# Print out the critical values
cat("Critical values for model with no intercept: Tau")
print(CRV1)
cat("Critical values for model with intercept: Tau_mu")
print(CRV2)
cat("Critical values for model with intercept: Tau_tau")
print(CRV3)


# Figure 2
par(mfrow=c(1,1))
plot(density(stat[,3]),col="red",xlim=c(-5,3),lwd=2,
     main="",xlab="")
curve(dnorm(x, mean=-0.1/0.05, sd=1), 
      col="black", lwd=2, add=TRUE, yaxt="n")
curve(dnorm(x, mean=0, sd=1), 
      col="grey", lwd=2, add=TRUE, yaxt="n")
#curve(dnorm(x, mean=-0.1/0.033, sd=1), 
#      col="blue", lwd=2, add=TRUE, yaxt="n")
abline(v=-3.45,col="red")
abline(v=-0.1/0.05,col="black")


# Figure 3 and unit root testing results for real GDP
data("npext")
y=ts(na.omit(npext$realgnp),freq=1,start=1909)
plot(y,ylab="")
summary(ur.df(y,type = "trend"))

# Figure 4
setwd("~/")    # Specify your working directory and place he data set in there
gratio=ts(read.csv(file="US.csv",header = TRUE),freq=4,start=c(1995,1))
par(mfrow=c(2,1))
x1=exp(gratio[,1])
plot(x1,ylim=c(0.55,0.75),main="consumption-income ratio",ylab="")
x2=exp(gratio[,2])
plot(x2,ylim=c(0.1,0.3),main="investment-income ratio",ylab="")


# ADF and ERS test statistics for great ratio
ADF=ur.df(x1,type="drift",lags=8,selectlags="BIC")
summary(ADF)

ADF=ur.ers(x1,type="DF-GLS",model="constant",lag.max=1)
summary(ADF)

ADF=ur.df(x2,type="drift",lags=8,selectlags="BIC")
summary(ADF)

ADF=ur.ers(x2,type="DF-GLS",model="constant",lag.max=1)
summary(ADF)

# The t-statistics and optimal level of significance
t=(-0.062+0.05)/0.030
library(OptSig)
Opt.sig.norm.test(ncp=(-0.02+0.05)/0.030,n=103,alternative="greater")

t=(-0.037661+0.05)/0.016
Opt.sig.norm.test(ncp=(-0.02+0.05)/0.016,n=103,alternative="greater")
