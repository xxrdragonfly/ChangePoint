#generate data
H0Data<-function(p,omega,d,n,L){
  q<-c(rep(omega/d,d),rep((1-omega)/(p-d),p-d))
  x<-rmultinom(L,n,q)
  return(x)
}
#frequency under H_0
qH0<-function(x){
  fre<-apply(x,1,sum)/N
  return(fre)
}
#estimate d
dhat<-function(x){
  qh0s<-sort(x,decreasing = TRUE)
  d<-numeric(0)
  for (i in 1:(p-1)) {
    d[i]<-(-1-qh0s[i]*qh0s[i+1])/(sqrt(1+(qh0s[i])^2)
                                  *sqrt(1+(qh0s[i+1])^2))
  }
  destimate<-which.max(d)
  return(destimate)
}

likelihood<-function(data){
  L<-nrow(data)
  if(ncol(data)==2){
    data<-cbind(data[,1],data[,2]-data[,1],data[,2])
    s<-2
  }else{
    s<-ncol(data)-1
  }
  m<-data[,1:s];n<-data[,s+1];
  
  p_hat<-lapply(as.list(1:s), function(i) sum(data[,i])/sum(n))
  p_hat<-unlist(p_hat)
  p_hat<-p_hat[p_hat!=0]
  h<-(-p_hat)%*%log(p_hat)
  h1<-numeric(L-1)
  h2<-numeric(L-1)
  for (i in 1:(L-1)) {
    p1c_hat<-lapply(as.list(1:s), function(j)
      sum(data[1:i,j])/sum(n[1:i]));
    p1c_hat<-unlist(p1c_hat)
    p2c_hat<-lapply(as.list(1:s), function(j)
      sum(data[(i+1):L,j])/sum(n[(i+1):L]));
    p2c_hat<-unlist(p2c_hat)
    p1c_hat<-p1c_hat[p1c_hat!=0]
    h1[i]<-(-p1c_hat)%*%log(p1c_hat)
    p2c_hat<-p2c_hat[p2c_hat!=0]
    h2[i]<-(-p2c_hat)%*%log(p2c_hat)
  }
  M<-lapply(as.list(1:(L-1)),
            function(i) h-(sum(n[1:i])/sum(n))*h1[i]
            -((sum(n[(i+1):L]))/sum(n))*h2[i] )
  M<-unlist(M)
  LL<-2*sum(n)*M
  jqLL<-numeric(0)
  for (i in 1:(L-1)) {
    jqLL[i]<-LL[i]*sum(n[1:i])*sum(n[(i+1):L])/N^2
  }  
  maxL<-max(LL)
  sumL<-sum(jqLL)/L
  taoestimate<-min(which.max(LL)) 
  return(list(maxL=maxL,sumL=sumL,taoestimate=taoestimate))
}

sA<-function(data,destimate,p){
  if(destimate>1){
    datasA<-apply(data[1:destimate,],2,sum)
  }else{
    datasA<-data[1,]
  }
  return(datasA)
}
sB<-function(data,destimate,p){
  if(destimate<p-1){
    datasB<-apply(data[(destimate+1):p,],2,sum)
  }else{
    datasB<-data[p,]
  }
  return(datasB)
}

p<-300;omega<-0.5;d<-6;#parameters
n<-100 #number of experiments
L<-300 #sample size
N<-L*n;
ep<-p;
destimate<-numeric(0);LA<-numeric(0);LB<-numeric(0);Q<-numeric(0);
taoe<-numeric(0); taoez<-numeric(0);taop<-numeric(0);

for (i in 1:1000) {
  data<-H0Data(p,omega,d,n,L)
  qyuan<-qH0(data)
  destimate<-dhat(qyuan)
  datasA<-sA(data,destimate,p)
  datasB<-sB(data,destimate,p)
  dataB<-rbind(data[(destimate+1):p,],datasA)
  dataA<-rbind(data[1:destimate,],datasB)
  dataB<-rbind(dataB,apply(dataB, 2, sum))
  dataA<-rbind(dataA,apply(dataA, 2, sum))
  
  resultA<-likelihood(t(dataA))
  LA[i]<-resultA$maxL
  resultB<-likelihood(t(dataB))
  LB[i]<-resultB$sumL
  rp<-(2*log(log(N))+(destimate/2)*log(log(log(N))))^2
  Q[i]<-(LB[i]+ep*I(LA[i]>rp)-(p-destimate+1)/6)
  /sqrt((p-destimate+1)/45)
}
sum(abs(Q)>1.96)/1000