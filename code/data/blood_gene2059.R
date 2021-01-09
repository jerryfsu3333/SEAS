truncate<-function(x,f,delta.n){
  y<-f(x)
  y[y<delta.n]<-delta.n
  y[y>1-delta.n]<-1-delta.n
  y
}


getnorm3<-function(x,y){
  f0<-ecdf(x[y==0])
  f1<-ecdf(x[y==1])
  n0<-sum(y==0)
  n1<-sum(y==1)
  n<-n0+n1
  delta.n0<-1/n0^2
  delta.n1<-1/n1^2
  
  v0<-rep(0,n)
  v0[y==0]<-qnorm(truncate(x[y==0],f0,delta.n0))
  v0[y==1]<-qnorm(truncate(x[y==1],f0,delta.n0))
  mu0.hat<-mean(v0[y==1])
  
  v1<-rep(0,n)
  v1[y==0]<-qnorm(truncate(x[y==0],f1,delta.n1))
  v1[y==1]<-qnorm(truncate(x[y==1],f1,delta.n1))
  mu1.hat<-mean(v1[y==0])
  mu.hat<-n0/n*mu0.hat-n1/n*mu1.hat  
  
  transform<-function(t){
    n0*qnorm(truncate(t,f0,delta.n0))/n+n1*(qnorm(truncate(t,f1,delta.n1))+mu.hat)/n
  }
  x.norm<-transform(x)
  
  list(x.norm=x.norm,f0=f0,f1=f1,mu.hat=mu.hat,transform=transform)
}


dat<-scan("~/Documents/GitHub/ssdr/Real_dataset/blood.txt")
x<-matrix(dat,ncol=71)
w<-x[2059,]
w<-scale(w)
y<-c(rep(0,22),rep(1,71-22))
w.den0<-density(w[y==0])
w.den1<-density(w[y==1])


f<-ecdf(w[y==1])
v<-f(w)
v[v<1/49^2]<-1/49^2
v[v>1-1/49^2]<-1-1/49^2
v<-qnorm(v)
v.den0<-density(v[y==0],bw=0.7)
v.den1<-density(v[y==1],bw=0.7)

w<-c(w)
x.norm<-getnorm3(w,y)$x.norm
x.norm.den0<-density(x.norm[y==0],bw=0.55)
x.norm.den1<-density(x.norm[y==1],bw=0.4)

par(mfrow=c(1,3))
plot(w.den0$x,w.den0$y,xlim=c(-2,3),ylim=c(0,4),type="l",lty=2,xlab="x",ylab="Density",main="Gene IRF1: Standardized",lwd=2)
lines(w.den1$x,w.den1$y,lwd=2)
legend("topright",legend=c("Positive","Negative"),lty=c(1,2))
plot(v.den0$x,v.den0$y,xlim=c(-6,3),ylim=c(0,0.5),type="l",lty=2,xlab="x",ylab="Density",main="Gene IRF1: After Naive transformation",lwd=2)
lines(v.den1$x,v.den1$y,lwd=2)
legend("topright",legend=c("Positive","Negative"),lty=c(1,2))
plot(x.norm.den0$x,x.norm.den0$y,xlim=c(-2.5,6),ylim=c(0,0.55),type="l",lty=2,xlab="x",ylab="Density",main="Gene IRF1: After Pooled transformation",lwd=2)
lines(x.norm.den1$x,x.norm.den1$y,lwd=2)
legend("topright",legend=c("Positive","Negative"),lty=c(1,2))


postscript("blood_gene2059_transformed.ps")
plot(v.den0$x,v.den0$y,xlim=c(-3,3),ylim=c(0,0.5),type="l",lty=2,xlab="x",ylab="Density",main="Density for Transformed Data")
lines(v.den1$x,v.den1$y)
legend("topright",legend=c("Positive","Negative"),lty=c(1,2))
dev.off()

w.den0<-density(w[y==0],bw=0.1)
w.den1<-density(w[y==1],bw=0.5)

density(w[y==1])$bw
density(v[y==0])$bw


postscript("blood_gene2059_raw.ps")
plot(w.den0$x,w.den0$y,xlim=c(-3,3),ylim=c(0,4),type="l",lty=2,xlab="x",ylab="Density",main="Density for Raw Data")
lines(w.den1$x,w.den1$y)
legend("topright",legend=c("Positive","Negative"),lty=c(1,2))
dev.off()

postscript("blood_gene2059.ps")
par(mfrow=c(1,2))
plot(w.den0$x,w.den0$y,xlim=c(-3,3),ylim=c(0,3.5),type="l",lty=2,xlab="x",ylab="Density",main="Density for Raw Data")
lines(w.den1$x,w.den1$y)
legend("topright",legend=c("Positive","Negative"),lty=c(1,2))
plot(v.den0$x,v.den0$y,xlim=c(-3,3),ylim=c(0,0.45),type="l",lty=2,xlab="x",ylab="Density",main="Density for Transformed Data")
lines(v.den1$x,v.den1$y)
legend("topright",legend=c("Positive","Negative"),lty=c(1,2))
dev.off()

