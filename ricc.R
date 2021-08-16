library(wle)
library(leaps)
library(MASS)

n.y<-50
n.s<-1
n.p<-30
devi<- 1

b<-5

for (rou.1 in c(-0.5,0,0.5))
{
output<-NULL
output<-matrix(0,10,7)
lossoutput<-matrix(0,10,7)
for (n.c in c(1:10)){

                      
                      beta.1<-logical(10)
                      beta.1[1:n.c]<-TRUE
                      beta.2<-rep(beta.1,(n.p/10))
                      true.model<-beta.2
                      beta<-beta.2*b
                      
                      Cov.Matrix<-matrix(0,n.p,n.p) 
                      X<-matrix(0, n.y, n.p)
                      mu1<-seq(0,0,length=n.p)
                      
                      loss.aic<-seq(0,length=n.s) 
                      loss.bic<-seq(0,length=n.s) 
                      loss.ric<-seq(0,length=n.s) 
                      loss.mric<-seq(0,length=n.s)
                      loss.cic<-seq(0,length=n.s)
                      loss.ricc<-seq(0,length=n.s)
                      
                      true.aic<-logical(n.s)
                      true.bic<-logical(n.s)
                      true.ric<-logical(n.s)
                      true.mric<-logical(n.s)
                      true.cic<-logical(n.s)
                      true.ricc<-logical(n.s)

                      
                      for (i in 1:n.p) 
                                      {
                                      Cov.Matrix[i,]<-seq(0,0,length=n.p) 
                                      for (j in 1:n.p) 
                                                      {
                                                      if (i==j)  { Cov.Matrix[i,j]<-1} 
                                                      else 
                                                      {Cov.Matrix[i,j]<-(rou.1)^abs(i-j)} 
                                                      }
                                      }
                      
                      
                      
                      
                      for (i in 1:n.s)  
                                        {
                                        #i<-1
                                        X <- mvrnorm(n.y, mu1, Cov.Matrix)
                                        X.1 <- cbind(1, X) 
                                        
                                        
                                        mu<-X%*%beta
                                        Y <- mu + rnorm(n.y,0,(devi)^2)
                                        
                                        
                                        
                                        mm <- leaps(X,Y,nbest=1,method="r2")
                                        model<-mm$which
                                        row.names(model)<-NULL
                                        md<-mm$size
                                        rss<-(n.y-1)*var(Y)*(1-mm$r2)
                                        aic<-rss+2*(devi)^2*md
                                        bic<-rss+log(n.y)*md*(devi)^2 
                                        ric<-rss+2*log(n.p)*md*(devi)^2 
                                        mric<-rss+2*devi^2*(log(n.p)*md-log(factorial(md))) 
                                        cic<-rss+4*devi^2*(log(n.p)*md-log(factorial(md))) 
                                        ricc<-rss+2*(log(n.p)+log(log(n.p)))*md*(devi)^2

                                        
                                        
                                        if (sum(mm$which[order(aic)[1],]-true.model==0)==n.p) true.aic[i]<-TRUE
                                        if (sum(mm$which[order(bic)[1],]-true.model==0)==n.p) true.bic[i]<-TRUE
                                        if (sum(mm$which[order(ric)[1],]-true.model==0)==n.p)  true.ric[i]<-TRUE
                                        if (sum(mm$which[order(mric)[1],]-true.model==0)==n.p) true.mric[i]<-TRUE
                                        if (sum(mm$which[order(cic)[1],]-true.model==0)==n.p)  true.cic[i]<-TRUE
                                        if (sum(mm$which[order(ricc)[1],]-true.model==0)==n.p)  true.ricc[i]<-TRUE
                                        
                                        loss.aic[i]<-sum((lm(Y~X[,mm$which[order(aic)[1],]])$fit-mu)^2)
                                        loss.bic[i]<-sum((lm(Y~X[,mm$which[order(bic)[1],]])$fit-mu)^2)
                                        loss.ric[i]<-sum((lm(Y~X[,mm$which[order(ric)[1],]])$fit-mu)^2)
                                        loss.mric[i]<-sum((lm(Y~X[,mm$which[order(mric)[1],]])$fit-mu)^2)
                                        loss.cic[i]<-sum((lm(Y~X[,mm$which[order(cic)[1],]])$fit-mu)^2)
                                        loss.ricc[i]<-sum((lm(Y~X[,mm$which[order(ricc)[1],]])$fit-mu)^2)
                                        }
                      
                      
                      proportion.aic<-sum(true.aic)/n.s
                      proportion.bic<-sum(true.bic)/n.s
                      proportion.ric<-sum(true.ric)/n.s
                      proportion.mric<-sum(true.mric)/n.s
                      proportion.cic<-sum(true.cic)/n.s
                      proportion.ricc<-sum(true.ricc)/n.s
                      
                      
                      meanloss.aic<-mean(loss.aic)
                      meanloss.bic<-mean(loss.bic)
                      meanloss.ric<-mean(loss.ric)
                      meanloss.ricc<-mean(loss.ricc)
                      output[n.c,]<-c(n.c*3,proportion.aic,proportion.bic,proportion.ric,proportion.mric,proportion.cic,proportion.ricc)
                      
                      }
                      
                      
                      
output<-data.frame(output)
names(output)<-c("md","AIC","BIC","RIC","MRIC","CIC","RICc")
file.name<-paste("C:/Documents and Settings/yongli/My Documents/TechReport/RICc/RICc_",rou.1,"_",n.y,".csv",sep="")

write.csv(output, file=file.name,row.names=F,na="NA")

}

output