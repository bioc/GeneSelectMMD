"myWilcox" <-
function(x, memSubjects, alpha =0.05, ...)
{
  xc<-x[memSubjects==1]
  xn<-x[memSubjects==0]

  m<-sum(memSubjects==1)
  res<-wilcox.test(x=xc, y=xn, conf.level = 1-alpha, ...)
  res2<-c(res$p.value, res$statistic-m*(m+1)/2)
  names(res2)<-c("p.value", "statistic")

  return(res2)
}

