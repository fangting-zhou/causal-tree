find.order=function(graph) {
  nodes=1:ncol(graph)
  order=which(rowSums(graph)==0)
  
  nodes=nodes[!(nodes%in%order)]
  while(length(nodes)>1) {
    candidate.node=which(rowSums(graph[,-order])==0)
    order=c(order,candidate.node[!(candidate.node%in%order)])
    nodes=nodes[!(nodes%in%order)]
  }
  
  if(length(nodes)==1) order=c(order,nodes)
  
  return(order)
}

regression.cart=function(x) {
  xtree<<-x;num.node=ncol(xtree);num.sample=nrow(xtree)
  varname=colnames(xtree)
  
  regression.formula<<-as.formula(paste(paste(varname[num.node],"~"),paste(varname[1:(num.node-1)],collapse="+")))
  new.regression.tree=tree(regression.formula,data=xtree,split="deviance")
  
  if(length(unique(new.regression.tree$where))==1) {
    tree.depth=1;tree.value=new.regression.tree$frame$yval[new.regression.tree$where]
  } else {
    cv.regression.tree=try(cv.tree(new.regression.tree,,prune.tree,10),silent=TRUE)
    if(!any(class(cv.regression.tree)=="try-error")) {
      tree.depth=min(cv.regression.tree$size[which(cv.regression.tree$dev==min(cv.regression.tree$dev))]) 
    } else tree.depth=1
    
    prune.regression.tree=prune.tree(new.regression.tree,best=tree.depth,method="deviance")
    tree.value=prune.regression.tree$frame$yval[prune.regression.tree$where]
  }
  
  residual=xtree[,num.node]-tree.value
  
  return(residual)
}

treeRESIT=function(x) {
  current.node=1:ncol(x)
  current.order=NULL
  
  while(length(current.node)>1) {
    pvalues=rep(0,length(current.node))
    for(j in 1:length(current.node)) {
      tempz=x[,current.node[-j]]
      tempx=x[,c(current.node[-j],current.node[j])]
      
      temp.residual=regression.cart(tempx)
      pvalues[j]=dhsic.test(tempz,temp.residual)$p.value
    }
    
    current.order=c(current.node[which.max(pvalues)],current.order)
    current.node=current.node[-which.max(pvalues)]
  }
  
  current.order=c(current.node,current.order)
  
  return(current.order)
}

prune.from.order=function(x,order) {
  xtree<<-x[,order];varname=colnames(xtree);num.node=ncol(xtree);num.sample=nrow(xtree)
  estDAG=matrix(0,num.node,num.node)
  
  for(j in 2:num.node) {
    regression.formula<<-as.formula(paste(paste(varname[j],"~"),paste(varname[1:(j-1)],collapse="+")))
    new.regression.tree=tree(regression.formula,data=xtree,split="deviance")
    
    if(length(unique(new.regression.tree$where))==1) tree.depth=1 else {
      cv.regression.tree=try(cv.tree(new.regression.tree,,prune.tree,10),silent=TRUE)
      if(!any(class(cv.regression.tree)=="try-error")) {
        tree.depth=min(cv.regression.tree$size[which(cv.regression.tree$dev==min(cv.regression.tree$dev))]) 
      } else tree.depth=1
      
      prune.regression.tree=prune.tree(new.regression.tree,best=tree.depth,method="deviance")
    }
    
    if(tree.depth!=1) estDAG[order[j],order[which(colnames(xtree)%in%prune.regression.tree$frame[,1])]]=1
  }
  
  return(estDAG)
}

score.cart=function(x,order) {
  xtree<<-x[,order];varname=colnames(xtree);num.sample=nrow(xtree)
  
  score=num.sample*log( mean((xtree[,1]-mean(xtree[,1]))^2))+log(num.sample)
  for(j in 2:ncol(xtree)) {
    regression.formula<<-as.formula(paste(paste(varname[j],"~"),paste(varname[1:(j-1)],collapse="+")))
    new.regression.tree=tree(regression.formula,data=xtree,split="deviance")
    
    if(length(unique(new.regression.tree$where))==1) {
      tree.depth=1;tree.value=new.regression.tree$frame$yval[new.regression.tree$where]
    } else {
      cv.regression.tree=try(cv.tree(new.regression.tree,,prune.tree,10),silent=TRUE)
      if(!any(class(cv.regression.tree)=="try-error")) {
        tree.depth=min(cv.regression.tree$size[which(cv.regression.tree$dev==min(cv.regression.tree$dev))]) 
      } else tree.depth=1
      
      prune.regression.tree=prune.tree(new.regression.tree,best=tree.depth,method="deviance")
      tree.value=prune.regression.tree$frame$yval[prune.regression.tree$where]
    }
    
    score=score+log(num.sample)*tree.depth+num.sample*log(mean((xtree[,j]-tree.value)^2))
  }
  
  return(score)
}

treeBruteForce=function(x) {
  list.order=permn(ncol(x))
  list.score=rep(0,length(list.order))
  
  for(k in 1:length(list.order)) list.score[k]=score.cart(x,list.order[[k]])
    
  return(list(order=list.order,score=list.score))
}

orderInfo=function(x,order) {
  xtree<<-x[,order];varname=colnames(xtree);num.sample=nrow(xtree);num.node=ncol(xtree)
  graph=matrix(0,num.node,num.node)
  
  score=num.sample*log(mean((xtree[,1]-mean(xtree[,1]))^2))+log(num.sample)
  for(j in 2:ncol(xtree)) {
    regression.formula<<-as.formula(paste(paste(varname[j],"~"),paste(varname[1:(j-1)],collapse="+")))
    new.regression.tree=tree(regression.formula,data=xtree,split="deviance")
    
    if(length(unique(new.regression.tree$where))==1) {
      tree.depth=1;tree.value=new.regression.tree$frame$yval[new.regression.tree$where]
    } else {
      cv.regression.tree=try(cv.tree(new.regression.tree,,prune.tree,10),silent=TRUE)
      if(!any(class(cv.regression.tree)=="try-error")) {
        tree.depth=min(cv.regression.tree$size[which(cv.regression.tree$dev==min(cv.regression.tree$dev))]) 
      } else tree.depth=1
      
      prune.regression.tree=prune.tree(new.regression.tree,best=tree.depth,method="deviance")
      tree.value=prune.regression.tree$frame$yval[prune.regression.tree$where]
    }
    
    if(tree.depth!=1) graph[order[j],order[which(colnames(xtree)%in%prune.regression.tree$frame[,1])]]=1
    score=score+log(num.sample)*tree.depth+num.sample*log(mean((xtree[,j]-tree.value)^2))
  }
  
  return(list(score=score,graph=graph))
}

search.candidate=function(x,graph,order.list,search.depth) {
  candiate.order=NULL;candidate.score=NULL;candidate.graph=NULL;wait.list=NULL;iter=1;num.node=ncol(x)
  for(j in 1:num.node) {
    for(l in 1:num.node) {
      if(graph[j,l]==1&all(graph[j,-l]==graph[l,-l])) {
        new.graph=graph;new.graph[j,l]=0;new.graph[l,j]=1
        new.order=find.order(new.graph)
        
        if(!any(order.list%in%list(new.order))) wait.list=c(wait.list,list(new.order))}}}
  
  while(length(wait.list)!=0&iter<search.depth) {
    new.list=NULL
    for(k in 1:length(wait.list)) {
      current.result=orderInfo(x,wait.list[[k]])
      candidate.score=c(candidate.score,current.result$score)
      candidate.graph=c(candidate.graph,list(current.result$graph))
      
      current.graph=current.result$graph
      for(j in 1:num.node) {
        for(l in 1:num.node) {
          if(current.graph[j,l]==1&all(current.graph[j,-l]==current.graph[l,-l])) {
            new.graph=current.graph;new.graph[j,l]=0;new.graph[l,j]=1
            new.order=find.order(new.graph)
            
            if(!any(c(order.list,candiate.order,wait.list,new.list)%in%list(new.order))) 
              new.list=c(new.list,list(new.order))}}}
    }
    
    candiate.order=c(candiate.order,wait.list)
    wait.list=new.list;iter=iter+1
  }
  
  if(length(wait.list)!=0) {
     candiate.order=c(candiate.order,wait.list)
     for(k in 1:length(wait.list)) {
       current.result=orderInfo(x,wait.list[[k]])
       current.order=c(candidate.score,current.result$score)
       current.graph=c(candidate.score,list(current.result$graph))
     }
  }
  
  return(list(order=candiate.order,score=candidate.score,graph=candidate.graph))
}

treeGSBIC=function(x,search.depth,initial.order) {
  current.result=orderInfo(x,initial.order)
  
  current.order=initial.order
  current.score=current.result$score
  current.graph=current.result$graph
  
  terminate=0;order.list=c(list(current.order))
  while(terminate==0) {
    candidate.result=search.candidate(x,current.graph,order.list,search.depth)
    candidate.order=candidate.result$order
    candidate.score=candidate.result$score
    candidate.graph=candidate.result$graph
    
    order.list=c(order.list,candidate.order)
    if(sum(candidate.score<current.score)!=0) {
      current.order=candidate.order[[which.min(candidate.score)]]
      current.graph=candidate.graph[[which.min(candidate.score)]]
      current.score=candidate.score[which.min(candidate.score)]
    } else terminate=1
  }
  
  return(list(order=current.order,score=current.score))
}