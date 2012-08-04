
# this was in wiFun.R - the wiFun function has been converted to Fortran
# so checkPara is in its own file

checkPara<-function(para, eps=1.0e-6)
{
  if(length(para) !=TNumParaRP)
  { stop("Number of parameters is not correct!\n") }

  # mixture proportions
  pi.1<-para[1]; pi.2<-para[2]; 
  pi.3<- 1-pi.1-pi.2

  if(pi.1>=1 || pi.1<=0)
  { cat("pi.1=", pi.1, " pi.2=", pi.2, " pi.3=",pi.3,"\n") 
    stop("pi.1 should be in (0, 1)!\n") }
  if(pi.2>=1 || pi.2<=0)
  { cat("pi.1=", pi.1, " pi.2=", pi.2, " pi.3=",pi.3,"\n") 
    stop("pi.2 should be in (0, 1)!\n") }
  if(pi.3>=1 || pi.3<=0)
  { cat("pi.1=", pi.1, " pi.2=", pi.2, " pi.3=",pi.3,"\n") 
    stop("pi.3 should be in (0, 1)!\n") }
  if(abs((pi.1+pi.2+pi.3)-1)>eps)
  { cat("pi.1=", pi.1, " pi.2=", pi.2, " pi.3=",pi.3,"\n") 
    stop("pi.1+pi.2+pi.3 should be equal to 1!\n") }

}


