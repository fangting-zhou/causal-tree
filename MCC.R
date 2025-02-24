mcc=function(trueDAG,estDAG) {
  tp=sum(estDAG[trueDAG==1])
  fp=sum(estDAG[trueDAG==0])
  tn=sum(1-estDAG[trueDAG==0])
  fn=sum(1-estDAG[trueDAG==1])
  
  mcc=(tp*tn-fp*fn)/sqrt((tp+fp)*(tp+fn)*(tn+fp)*(tn+fn))
  
  return(mcc)
}