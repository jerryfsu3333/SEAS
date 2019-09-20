library(GEOquery)
library(randomForest)
library(energy)

dat<-getGEO("GSE5680")

x<-exprs(dat[[1]])

"1389163_at" %in% featureNames(dat)

y<-2^x[featureNames(dat)=="1389163_at",]

x<-t(x)
#x<-x[,1:31042]

x.exp<-round(2^x,digits=3)

x.exp.max<-apply(x.exp,2,max)
x.exp.min<-apply(x.exp,2,min)

id2<-which(x.exp.max/x.exp.min>2.0005)

id1<-which(x.exp.max>quantile(c(x.exp))[2])

id3<-1:31042

id<-intersect(id1,id2)

id<-id[id<31042]

#id3<-which(x.exp.max/x.exp.min>4)

#id<-intersect(id,id3)

x.r<-x.exp[,featureNames(dat)!="1389163_at"]
x.r<-x.exp[,id]
x.r<-x.r[y>270&y<400,]
y<-y[y>270&y<400]

p<-ncol(x.r)
n<-nrow(x.r)

# nruns<-10
# dc.MSE<-rep(0,nruns)
# ks.MSE<-rep(0,nruns)
# n.train<-90
# 
# dc.res<-matrix(0,n-n.train,nruns)
# ks.res<-matrix(0,n-n.train,nruns)
# 
# 
# for(iruns in 1:nruns){
#   id.train<-sample(1:n,n.train,replace=F)
#   x.train<-log(x.r[id.train,])
#   y.train<-log(y[id.train])
#   x.test<-log(x.r[-id.train,])
#   y.test<-log(y[-id.train])
#   
#   dc.stat<-rep(0,p)
#   
#   for(i in 1:p){
#   	dc.stat[i]<-dcor(exp(x.train[,i]),exp(y.train))
#   }
#   
#   
#   id.scr<-which(dc.stat>=sort(dc.stat,decreasing=T)[n.train/log(n.train)])
#   r<-randomForest(x.train[,id.scr],y.train)
#   pred<-predict(r,x.test[,id.scr])
#   dc.MSE[iruns]<-median(abs(pred-y.test))
#   dc.res[,iruns]<-abs(pred-y.test)
#   
#   ks.stat.max<-rep(0,p)
#   
#   K.max<-ceiling(log(n.train))-2
#   for(K in 3:(2+K.max)){
#     y.dm<-cut(y.train,quantile(y.train,seq(0,1,1/K)),labels=c(1:K),include.lowest=T)
#     ks.stat<-matrix(0,p,K*(K-1)/2)
#     
#     nclass<-0
#     for(j in 1:(K-1)){
#       for(l in (j+1):K){
#         nclass<-nclass+1
#         for(i in 1:p){
#           ks.stat[i,nclass]<-ks.test(x.train[y.dm==j,i],x.train[y.dm==l,i])$statistic}
#       }
#     }  
#     ks.stat.max0<-apply(ks.stat,1,max)
#     ks.stat.max<-ks.stat.max+ks.stat.max0
#   }
#         
#   id.scr<-which(ks.stat.max>=sort(ks.stat.max,decreasing=T)[n.train/log(n.train)])
#   r<-randomForest(x.train[,id.scr],y.train)
#   pred<-predict(r,x.test[,id.scr])
#   ks.MSE[iruns]<-median(abs(pred-y.test))
#   ks.res[,iruns]<-abs(pred-y.test)
# }
