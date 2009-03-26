"modifyPartition"<-
function(X, memSubjects, memGenes, geneNames=NULL, alpha=0.05, eps=1e-06)
{
  nGenes<-nrow(X)
  size<-table(memGenes)
  pos<-which(size<2)
  if(length(pos)>0)
  {
    res<-iniMemGenesTestFunc(X, memSubjects, myTtest, geneNames, alpha, eps=eps)
    return(res)
  }

  tmpMemGenes<-memGenes
  tmp<-getPara(X, memSubjects, tmpMemGenes, eps=eps)
  para<-tmp$para
  llkh<-tmp$llkh

  tmpMemGenes2<-memGenes
  tmpMemGenes2[memGenes==1]<-1
  tmpMemGenes2[memGenes==2]<-3
  tmpMemGenes2[memGenes==3]<-2
  tmp2<-getPara(X, memSubjects, tmpMemGenes2, eps=eps)
  para2<-tmp2$para
  llkh2<-tmp2$llkh

  if(llkh2>llkh && para2[2]>=para2[1] && para2[2]>=para2[3])
  { llkh<-llkh2
    para<-para2
    tmpMemGenes<-tmpMemGenes2
  }

  tmpMemGenes3<-memGenes
  tmpMemGenes3[memGenes==1]<-2
  tmpMemGenes3[memGenes==2]<-3
  tmpMemGenes3[memGenes==3]<-1
  tmp3<-getPara(X, memSubjects, tmpMemGenes3, eps=eps)
  para3<-tmp3$para
  llkh3<-tmp3$llkh

  if(llkh3>llkh && para3[2]>=para3[1] && para3[2]>=para3[3])
  { llkh<-llkh3
    para<-para3
    tmpMemGenes<-tmpMemGenes3
  }

  tmpMemGenes4<-memGenes
  tmpMemGenes4[memGenes==1]<-2
  tmpMemGenes4[memGenes==2]<-1
  tmpMemGenes4[memGenes==3]<-3
  tmp4<-getPara(X, memSubjects, tmpMemGenes4, eps=eps)
  para4<-tmp4$para
  llkh4<-tmp4$llkh

  if(llkh4>llkh && para4[2]>=para4[1] && para4[2]>=para4[3])
  { llkh<-llkh4
    para<-para4
    tmpMemGenes<-tmpMemGenes4
  }

  tmpMemGenes5<-memGenes
  tmpMemGenes5[memGenes==1]<-3
  tmpMemGenes5[memGenes==2]<-1
  tmpMemGenes5[memGenes==3]<-2
  tmp5<-getPara(X, memSubjects, tmpMemGenes5, eps=eps)
  para5<-tmp5$para
  llkh5<-tmp5$llkh

  if(llkh5>llkh && para5[2]>=para5[1] && para5[2]>=para5[3])
  { llkh<-llkh5
    para<-para5
    tmpMemGenes<-tmpMemGenes5
  }

  tmpMemGenes6<-memGenes
  tmpMemGenes6[memGenes==1]<-3
  tmpMemGenes6[memGenes==2]<-2
  tmpMemGenes6[memGenes==3]<-1
  tmp6<-getPara(X, memSubjects, tmpMemGenes6, eps=eps)
  para6<-tmp6$para
  llkh6<-tmp6$llkh

  if(llkh6>llkh && para6[2]>=para6[1] && para6[2]>=para6[3])
  { llkh<-llkh6
    para<-para6
    tmpMemGenes<-tmpMemGenes6
  }

  if(sum(is.null(geneNames)))
  { geneNames<-paste("gene", 1:nGenes, sep="") }

  names(tmpMemGenes)<-geneNames

  memGenes<-tmpMemGenes
  memGenes2<-rep(1, nGenes)
  memGenes2[memGenes==2]<-0

  if(llkh== -Inf || para[2]<para[1] || para[2]<para[3])
  {
    res<-iniMemGenesTestFunc(X, memSubjects, myTtest, geneNames, alpha, eps=eps)
    return(res)
  } else {
    return(list(para=para, llkh=llkh, memGenes=memGenes, memGenes2=memGenes2))
  }
}

