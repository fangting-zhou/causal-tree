generate.polynomial=function(x) {
  power.series=sapply(x,function(z) z^(0:5))
  if(length(x)==1) poly.basis=power.series else {
    poly.basis=apply(expand.grid(power.series[,1],power.series[,2]),1,prod)
    if(length(x)>=3) for(j in 3:length(x)) poly.basis=apply(expand.grid(poly.basis,power.series[,j]),1,prod)
  }
  
  return(poly.basis)
}

generate.data=function(n,p,g0) {
  x=matrix(0,n,p)
  for(j in 1:p) {
    if(sum(g0[j,])==0) x[,j]=rnorm(n,0,1) else {
      if(sum(g0[j,])==1) degree=(0:5) else {
        degree=apply(expand.grid(0:5,0:5),1,sum)
        if(sum(g0[j,])>=3) degree=apply(expand.grid(degree,0:5),1,sum)
      }
      
      if(sum(g0[j,])==1) sel=rep(0,6) else {
        sel=apply(expand.grid(c(1,rep(0,5)),c(1,rep(0,5))),1,sum)
        if(sum(g0[j,])>=3) sel=apply(expand.grid(sel,c(1,rep(0,5))),1,sum)
      }
      
      scalex=as.matrix(x[,which(g0[j,]!=0)],ncol=sum(g0[j,]))
      for(k in 1:ncol(scalex)) scalex[,k]=2*((scalex[,k]-min(scalex[,k]))/(max(scalex[,k])-min(scalex[,k]))-0.5)
                   
      poly.basis=t(apply(scalex,1,generate.polynomial))
      poly.coeff=runif(ncol(poly.basis),0.5,1)*(2*(runif(ncol(poly.basis))<0.5)-1)
      if(sum(g0[j,])>=2) poly.coeff[sel==(sum(g0[j,])-1)]=(runif(1)<0.5)*poly.coeff[sel==(sum(g0[j,])-1)]
      
      x[,j]=poly.basis%*%poly.coeff+sd(poly.basis%*%poly.coeff)*rnorm(n,0,1)
    }
  }
  
  return(list(g0=g0,x=x))
}
