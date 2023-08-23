# This R code reproduces Figures 5, 6, 7 of the paper (DF-GLS test)
rm(list=ls(all=TRUE))
# Install urca package and fUnitroots package
library("fUnitRoots");library("urca")
options(warn=-1)
set.seed(5)

# number of Monte Carlo Trials
nit=10000 

# parameter settings for data
a0=0; a1=0; b=0.95; x0=1.5; sig = 1
g = a0*(1-b) + a1*b; d=a1*(1-b)

n1vec=c(80, 100, 120);   # Sample size

del=0.01;
alphavec=seq(0.01,0.99,del)  # alpha values
statmat=numeric()
for(jj in 1:length(n1vec)){
n1=n1vec[jj]
cr=qunitroot(alphavec,N=n1,trend="nc")
mat=numeric()
reject=matrix(NA,nrow=nit,ncol=length(cr))


# Monte Carlo simulation
for( i in 1:nit){
u= arima.sim(list(order = c(1,0,0), ar = 0.5), n = n1)
y = a0+sig*x0

# generate data
for(i1 in 2:n1){
y1= g+d*i+b*y[i1-1]+u[i1]
y=c(y,y1)}

# statistic
ADF=ur.ers(y,type="DF-GLS",model="constant",lag.max=1)
adft=ADF@teststat

# Power and beta calculation
for (i2 in 1:length(cr)){
if(adft < cr[i2]) reject[i,i2] = 1
if(adft >= cr[i2]) reject[i,i2] = 0
}
power1=colMeans(reject)
beta1=1-power1
}

statmat = cbind(statmat,beta1)

}

# plots
P=0.5
dat=statmat
plot(dat[,1],alphavec,xlim=c(0,1),ylim=c(0,1),type="l",col=1,lwd=2,ylab="alpha",xlab="beta",main="DF-GLS"); #text(x=0.9,y=0.5,"n=40",col=2)
points(dat[,2],alphavec,type="l",col=1,lwd=2); #text(x=0.9,y=0.4,"n=60",col=4)
points(dat[,3],alphavec,type="l",col=1,lwd=2); #text(x=0.9,y=0.4,"n=80",col=4)
abline(h=0.05,col=2,lwd=2)
axis(side=2,at=seq(0.1,1,0.1))
abline(v=seq(0.1,1,0.1), col="lightgray", lty="dotted")
abline(h=seq(0.1,1,0.1), col="lightgray", lty="dotted")

# Points of minimizing the expected loss
alpha1=alphavec; 
for (i in 1:ncol(dat)){
beta1=dat[,i]
dd=cbind(alpha1,beta1,abs(P*alpha1+(1-P)*beta1))
dd1= dd[dd[,1] == 0.05,]
dd2= dd[dd[,3] == min(dd[,3]),]
#print(dd1); 
print(dd2)
points(dd2[2],dd2[1],col=4,pch=15,cex=1.5)
}
