source('Requirements.R')

set.seed(123)

n=100; p=6; d=2
g0=matrix(0,p,p); g0[lower.tri(g0)]=as.numeric(runif(p*(p-1)/2)<d/(p-1))
data=generate.data(n,p,g0); x=as.data.frame(data$x)

############# tree-GSBIC Experiment #############
order.list=NULL; score.list=NULL; num.run=20
for(j in 1:num.run) {
  initial.order=sample(1:ncol(x))
  candidate.order=treeGSBIC(x,4,initial.order)

  order.list=c(order.list,list(candidate.order$order))
  score.list=c(score.list,candidate.order$score)
}

gTreeGSBIC=prune.from.order(x,order.list[[which.min(score.list)]])

shd(dag2cpdag(as(t(gTreeGSBIC),"graphNEL")),dag2cpdag(as(t(g0),"graphNEL")))
structIntervDist(g0,gTreeGSBIC)$sid
mcc(g0,gTreeGSBIC)

############# CAM Experiment #############
gCAM=CAM(x,variableSel=TRUE,pruning=TRUE)$Adj

shd(dag2cpdag(as(gCAM,"graphNEL")),dag2cpdag(as(t(g0),"graphNEL")))
structIntervDist(g0,t(gCAM))$sid
mcc(g0,t(gCAM))

############# GP-RESIT experiment #############
gRESIT=ICML(x,alpha=0.05,model=train_gp)

shd(dag2cpdag(as(gRESIT,"graphNEL")),dag2cpdag(as(t(data$g0),"graphNEL")))
structIntervDist(data$g0,t(gRESIT))$sid
mcc(data$g0,t(gRESIT))