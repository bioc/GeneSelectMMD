"myTtest" <-
function(x, memSubjects, alpha = 0.05, ...)
{
  xc<-x[memSubjects==1]
  xn<-x[memSubjects==0]

  res<-t.test(x=xc, y=xn, conf.level=1-alpha, ...)

  res2<-c(res$p.value, res$statistic)
  names(res2)<-c("p.value", "statistic")
  return(res2)
}

