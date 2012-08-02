!     converted to Fortran 77, June 2012
!
!# modified on Sept. 28, 2011
!#  (1) added 'na.rm=TRUE' to function 'sum'
!#
!# sum_{i=1}^{nGenes} log10(pi.1*f1(xi)+pi.2*f2(xi)+pi.3*f3(xi))

!# negative of the M-step's objective function - the part relating to pi's
!# Q(Psi)=sum_{i=1}^{nGenes} wi*V(pi)+\sum_{i=1}^{nGenes) wi*U(theta)


!     NOTE: removed dat argument
!      function negqfunc(theta, dat, nc, nn, n, &

      SUBROUTINE negqfunc(negq, theta, nc, nn, n,
     +      sumw1, sumw2, sumw3,
     +      sumw1xcTxc, sumw1xcT1, sumw1xcT1_sq,
     +      sumw1xnTxn, sumw1xnT1, sumw1xnT1_sq,
     +      sumw2xTx, sumw2xT1, sumw2xT1_sq,
     +      sumw3xcTxc, sumw3xcT1, sumw3xcT1_sq,
     +      sumw3xnTxn, sumw3xnT1, sumw3xnT1_sq  )

      implicit none

      double precision theta(15)
      double precision negq

      double precision nc, nn, n

      double precision sumw1, sumw2, sumw3
      double precision sumw1xcTxc, sumw1xcT1, sumw1xcT1_sq
      double precision sumw1xnTxn, sumw1xnT1, sumw1xnT1_sq
      double precision sumw2xTx, sumw2xT1, sumw2xT1_sq
      double precision sumw3xcTxc, sumw3xcT1, sumw3xcT1_sq
      double precision sumw3xnTxn, sumw3xnT1, sumw3xnT1_sq

      double precision mu_c1, tau_c1, r_c1, rho_c1, delta_n1, mu_n1
      double precision tau_n1, r_n1, rho_n1, mu_2, tau_2, r_2, rho_2
      double precision mu_c3, tau_c3, r_c3, rho_c3, delta_n3, mu_n3
      double precision tau_n3, r_n3, rho_n3

      double precision sumw1acTac,sumw1acT1_sq,sumw1anTan,sumw1anT1_sq
      double precision sumw2aTa,sumw2aT1_sq,sumw3acTac,sumw3acT1_sq
      double precision sumw3anTan,sumw3anT1_sq

      double precision tt, part1, part2, part3
      double precision part1_1, part1_2, part1_3
      double precision part1_4, part1_5, part1_6
      double precision part2_1, part2_2, part2_3
      double precision part3_1, part3_2, part3_3, part3_4
      double precision part3_5, part3_6

      logical risnan

      double precision pi

      pi = 2d0*asin(1d0)


!     mean expression level for cluster 1 for diseased subjects
      mu_c1 = theta(1)
!     log(variance) of expression levels for cluster 1 for diseased subjects
      tau_c1 = theta(2) 
!     modified logit of correlation among expression levels for cluster 1 for diseased subjects
      r_c1 = theta(3)
 
      rho_c1 = (exp(r_c1)-1d0/(nc-1d0))/(1d0+exp(r_c1))

!     mean expression level for cluster 1 for normal subjects
      delta_n1 = theta(4)
      mu_n1 =  mu_c1-exp(delta_n1)
!     log(variance) of expression levels for cluster 1 for normal subjects
      tau_n1 = theta(5)
!     modified logit of correlation among expression levels for cluster 1 for normal subjects
      r_n1 =  theta(6)
      rho_n1 = (exp(r_n1)-1d0/(nn-1d0))/(1d0+exp(r_n1))
!     mean expression level for cluster 2
      mu_2 = theta(7)
!     log(variance) of expression levels for cluster 2
      tau_2 = theta(8)
!     modified logit of correlation among expression levels for cluster 2
      r_2 = theta(9)
      rho_2 = (exp(r_2)-1d0/(n-1d0))/(1d0+exp(r_2))

!     mean expression level for cluster 3 for diseased subjects
      mu_c3 =  theta(10)
!     log(variance) of expression levels for cluster 3 for diseased subjects
      tau_c3 = theta(11)
!     modified logit of correlation among expression levels for cluster 3 for diseased subjects
      r_c3 = theta(12)

      rho_c3 = (exp(r_c3)-1d0/(nc-1d0))/(1d0+exp(r_c3))

!     mean expression level for cluster 3 for normal subjects
      delta_n3 = theta(13)
      mu_n3 = mu_c3+exp(delta_n3)
!     log(variance) of expression levels for cluster 3 for normal subjects
      tau_n3 = theta(14)
!     modified logit of correlation among expression levels for cluster 3 for normal subjects
      r_n3 = theta(15)
      rho_n3 = (exp(r_n3)-1d0/(nn-1d0))/(1d0+exp(r_n3))

!     ###########
      sumw1acTac = sumw1xcTxc-2d0*mu_c1*sumw1xcT1+nc*mu_c1**2*sumw1
      sumw1acT1_sq = sumw1xcT1_sq-2d0*nc*mu_c1*sumw1xcT1+
     +        nc**2*mu_c1**2*sumw1

      sumw1anTan = sumw1xnTxn-2d0*mu_n1*sumw1xnT1+nn*mu_n1**2*sumw1
      sumw1anT1_sq = sumw1xnT1_sq-2d0*nn*mu_n1*sumw1xnT1+
     +        nn**2*mu_n1**2*sumw1
 
      sumw2aTa = sumw2xTx-2d0*mu_2*sumw2xT1+n*mu_2**2*sumw2
      sumw2aT1_sq = sumw2xT1_sq-2d0*n*mu_2*sumw2xT1+
     +        n**2*mu_2**2*sumw2

      sumw3acTac = sumw3xcTxc-2d0*mu_c3*sumw3xcT1+nc*mu_c3**2*sumw3
      sumw3acT1_sq = sumw3xcT1_sq-2d0*nc*mu_c3*sumw3xcT1+
     +        nc**2*mu_c3**2*sumw3
 
      sumw3anTan = sumw3xnTxn-2d0*mu_n3*sumw3xnT1+nn*mu_n3**2*sumw3
      sumw3anT1_sq = sumw3xnT1_sq-2d0*nn*mu_n3*sumw3xnT1+
     +        nn**2*mu_n3**2*sumw3

!     ##########

      part1_1 = -nc*log(2d0*pi)/2d0-nc*tau_c1/2d0-nc*log(nc)/2d0+
     +        (nc-1d0)*log(nc-1d0)/2d0

      tt = exp(r_c1)
      if (tt .lt. 1.0d308) then
            part1_1 = part1_1+nc*log(1d0+exp(r_c1))/2d0-r_c1/2d0
      else
            part1_1 = part1_1+nc*r_c1/2d0-r_c1/2d0
      endif
      part1_1 = part1_1 * sumw1

      part1_2 = -(nc-1d0)*(exp(-tau_c1)+
     +        exp(r_c1-tau_c1))/(2d0*nc)*sumw1acTac

      part1_3 = ((nc-2d0)*exp(-tau_c1)+exp(r_c1-tau_c1)*(nc-1d0)-
     +        exp(-r_c1-tau_c1))/(2d0*nc**2)*sumw1acT1_sq

      part1_4 = -nn*log(2d0*pi)/2d0-nn*tau_n1/2d0-nn*log(nn)/2d0+
     +        (nn-1d0)*log(nn-1d0)/2d0

      tt = exp(r_n1)
      if (tt .lt. 1.0d308) then
            part1_4 = part1_4+nn*log(1d0+exp(r_n1))/2d0-r_n1/2d0
      else
            part1_4 = part1_4+nn*r_n1/2d0-r_n1/2d0
      endif
      part1_4 = part1_4 * sumw1

      tt = exp(r_n1-tau_n1)
      if (tt .lt. 1.0d308) then
            part1_5 = -(nn-1d0)*(exp(-tau_n1)+
     +        exp(r_n1-tau_n1))/(2d0*nn)*sumw1anTan

            part1_6 = ((nn-2d0)*exp(-tau_n1)+exp(r_n1-tau_n1)*(nn-1d0)-
     +        exp(-r_n1-tau_n1))/(2d0*nn**2)*sumw1anT1_sq

            part1 = part1_1+part1_2+part1_3+part1_4+part1_5+part1_6
      else
            part1 = part1_1+part1_2+part1_3+part1_4
      endif

!     #
!     #  part1<-part1.1+part1.2+part1.3+part1.4+part1.5+part1.6
!     #
!     ##########
      part2_1 = -n*log(2d0*pi)/2d0-n*tau_2/2d0-n*log(n)/2d0+
     +            (n-1d0)*log(n-1d0)/2d0
      tt = exp(r_2)
      if (tt .lt. 1.0d308) then
            part2_1 = part2_1+n*log(1d0+exp(r_2))/2d0-r_2/2d0
      else
            part2_1 = part2_1+n*r_2/2d0-r_2/2d0
      endif
      part2_1 = part2_1 * sumw2

      part2_2 = -(n-1d0)*(exp(-tau_2)+exp(r_2-tau_2))/(2d0*n)*sumw2aTa
      part2_3 = ((n-2d0)*exp(-tau_2)+exp(r_2-tau_2)*(n-1d0)-
     +        exp(-r_2-tau_2))/(2d0*n**2)*sumw2aT1_sq

      if (.not. risnan(part2_2)) then
            part2 = part2_1+part2_2+part2_3
      else
            part2 = part2_1+part2_3
      endif


!     ##########

      part3_1 = -nc*log(2d0*pi)/2d0-nc*tau_c3/2d0-
     +        nc*log(nc)/2d0+(nc-1d0)*log(nc-1d0)/2d0

      tt = exp(r_c3)
      if (tt .lt. 1.0d308) then
            part3_1 = part3_1+nc*log(1d0+exp(r_c3))/2d0-r_c3/2d0
      else
            part3_1 = part3_1+nc*r_c3/2d0-r_c3/2d0
      endif
      part3_1 = part3_1 * sumw3

      part3_2 = -(nc-1d0)*(exp(-tau_c3)+
     +        exp(r_c3-tau_c3))/(2d0*nc)*sumw3acTac
      part3_3 = ((nc-2d0)*exp(-tau_c3)+exp(r_c3-tau_c3)*(nc-1d0)-
     +        exp(-r_c3-tau_c3))/(2d0*nc**2)*sumw3acT1_sq

      part3_4 = -nn*log(2d0*pi)/2d0-nn*tau_n3/2d0-
     +        nn*log(nn)/2d0+(nn-1d0)*log(nn-1d0)/2d0

      tt = exp(r_n3)
      if (tt .lt. 1.0d308) then
            part3_4 = part3_4+nn*log(1d0+exp(r_n3))/2d0-r_n3/2d0
      else
            part3_4 = part3_4+nn*r_n3/2d0-r_n3/2d0
      endif
      part3_4 = part3_4 * sumw3

      part3_5 = -(nn-1d0)*(exp(-tau_n3)+
     +        exp(r_n3-tau_n3))/(2d0*nn)*sumw3anTan
      part3_6 = ((nn-2d0)*exp(-tau_n3)+exp(r_n3-tau_n3)*(nn-1d0)-
     +        exp(-r_n3-tau_n3))/(2d0*nn**2)*sumw3anT1_sq

      part3 = part3_1+part3_2+part3_3+part3_4+part3_5+part3_6


!     #Qfunc<-part0+part1+part2+part3
      negq = part1+part2+part3

      if (.not. risnan(negq)) then
            negq = -negq
      else
            negq = (1.0d308)
      endif


      return
      end




!#    first derivatives of the negative of the objective function
!#    Q(Psi)=sum_{i=1}^{nGenes} wi*V(pi)+\sum_{i=1}^{nGenes) wi*U(theta)


!     NOTE: removed dat argument
!      function dnegqfunc(theta, dat, nc, nn, n, &

      SUBROUTINE dnegqfunc(dnegq, theta, nc, nn, n, 
     +      sumw1, sumw2, sumw3, 
     +      sumw1xcTxc, sumw1xcT1, sumw1xcT1_sq, 
     +      sumw1xnTxn, sumw1xnT1, sumw1xnT1_sq, 
     +      sumw2xTx, sumw2xT1, sumw2xT1_sq, 
     +      sumw3xcTxc, sumw3xcT1, sumw3xcT1_sq, 
     +      sumw3xnTxn, sumw3xnT1, sumw3xnT1_sq  )

      implicit none

      double precision theta(15)
      double precision dnegq(15)

      double precision nc, nn, n

      double precision sumw1, sumw2, sumw3
      double precision sumw1xcTxc, sumw1xcT1, sumw1xcT1_sq
      double precision sumw1xnTxn, sumw1xnT1, sumw1xnT1_sq
      double precision sumw2xTx, sumw2xT1, sumw2xT1_sq
      double precision sumw3xcTxc, sumw3xcT1, sumw3xcT1_sq
      double precision sumw3xnTxn, sumw3xnT1, sumw3xnT1_sq

      double precision mu_c1, tau_c1, r_c1, rho_c1, delta_n1, mu_n1
      double precision tau_n1, r_n1, rho_n1, mu_2, tau_2, r_2, rho_2
      double precision mu_c3, tau_c3, r_c3, rho_c3, delta_n3, mu_n3
      double precision tau_n3, r_n3, rho_n3

      double precision sumw1acTac,sumw1acT1_sq,sumw1anTan,sumw1anT1_sq
      double precision sumw2aTa,sumw2aT1_sq,sumw3acTac,sumw3acT1_sq
      double precision sumw3anTan,sumw3anT1_sq

      double precision dqdmu_c1,dqddelta_n1,dqdmu_2,dqdmu_c3,dqddelta_n3
      double precision dqdtau_c1,dqdtau_n1,dqdtau_2,dqdtau_c3,dqdtau_n3
      double precision dqdr_c1,dqdr_n1,dqdr_2,dqdr_c3,dqdr_n3
      double precision tt

      logical risnan

      integer i


!     mean expression level for cluster 1 for diseased subjects
      mu_c1 = theta(1)
!     log(variance) of expression levels for cluster 1 for diseased subjects
      tau_c1 = theta(2)
!     modified logit of correlation among expression levels for cluster 1 for diseased subjects
      r_c1 = theta(3)
      rho_c1 = (exp(r_c1)-1d0/(nc-1d0))/(1d0+exp(r_c1))

!     mean expression level for cluster 1 for normal subjects
      delta_n1 = theta(4)
      mu_n1 = mu_c1-exp(delta_n1)
!     log(variance) of expression levels for cluster 1 for normal subjects
      tau_n1 = theta(5)
!     modified logit of correlation among expression levels for cluster 1 for normal subjects
      r_n1 = theta(6)
      rho_n1 = (exp(r_n1)-1d0/(nn-1d0))/(1d0+exp(r_n1))

!     mean expression level for cluster 2
      mu_2 = theta(7)
!     log(variance) of expression levels for cluster 2
      tau_2 = theta(8)
!     modified logit of correlation among expression levels for cluster 2
      r_2 = theta(9)
      rho_2 = (exp(r_2)-1d0/(n-1d0))/(1d0+exp(r_2))

!     mean expression level for cluster 3 for diseased subjects
      mu_c3 = theta(10)
!     log(variance) of expression levels for cluster 3 for diseased subjects
      tau_c3 = theta(11)
!     modified logit of correlation among expression levels for cluster 3 for diseased subjects
      r_c3 = theta(12)
      rho_c3 = (exp(r_c3)-1d0/(nc-1d0))/(1d0+exp(r_c3))

!     mean expression level for cluster 3 for normal subjects
      delta_n3 = theta(13)
      mu_n3 = mu_c3+exp(delta_n3)
!     log(variance) of expression levels for cluster 3 for normal subjects
      tau_n3 = theta(14)
!     modified logit of correlation among expression levels for cluster 3 for normal subjects
      r_n3 = theta(15)
      rho_n3 = (exp(r_n3)-1d0/(nn-1d0))/(1d0+exp(r_n3))

!     ###########
      sumw1acTac = sumw1xcTxc-2d0*mu_c1*sumw1xcT1+nc*mu_c1**2*sumw1
      sumw1acT1_sq = sumw1xcT1_sq-2d0*nc*mu_c1*sumw1xcT1+
     +        nc**2*mu_c1**2*sumw1
 
      sumw1anTan = sumw1xnTxn-2d0*mu_n1*sumw1xnT1+nn*mu_n1**2*sumw1
      sumw1anT1_sq = sumw1xnT1_sq-2d0*nn*mu_n1*sumw1xnT1+
     +        nn**2*mu_n1**2*sumw1
 
      sumw2aTa = sumw2xTx-2d0*mu_2*sumw2xT1+n*mu_2**2*sumw2
      sumw2aT1_sq = sumw2xT1_sq-2d0*n*mu_2*sumw2xT1+
     +        n**2*mu_2**2*sumw2
 
      sumw3acTac = sumw3xcTxc-2d0*mu_c3*sumw3xcT1+nc*mu_c3**2*sumw3
      sumw3acT1_sq = sumw3xcT1_sq-2d0*nc*mu_c3*sumw3xcT1+
     +        nc**2*mu_c3**2*sumw3
 
      sumw3anTan = sumw3xnTxn-2d0*mu_n3*sumw3xnT1+nn*mu_n3**2*sumw3
      sumw3anT1_sq = sumw3xnT1_sq-2d0*nn*mu_n3*sumw3xnT1+
     +        nn**2*mu_n3**2*sumw3


!     ##########
      dqdmu_c1 = (exp(-tau_c1)+exp(-r_c1-tau_c1))*
     +        (sumw1xcT1-nc*mu_c1*sumw1)/nc
      dqdmu_c1 = dqdmu_c1+(exp(-tau_n1)+exp(-r_n1-tau_n1))*
     +        (sumw1xnT1-nn*mu_n1*sumw1)/nn

      dqddelta_n1 = -exp(delta_n1)*(exp(-tau_n1)+exp(-r_n1-tau_n1))*
     +        (sumw1xnT1-nn*mu_n1*sumw1)/nn

      dqdmu_2 = (exp(-tau_2)+exp(-r_2-tau_2))*(sumw2xT1-n*mu_2*sumw2)/n

      dqdmu_c3 = (exp(-tau_c3)+exp(-r_c3-tau_c3))*
     +        (sumw3xcT1-nc*mu_c3*sumw3)/nc

      dqdmu_c3 = dqdmu_c3+(exp(-tau_n3)+exp(-r_n3-tau_n3))*
     +        (sumw3xnT1-nn*mu_n3*sumw3)/nn

      dqddelta_n3 = exp(delta_n3)*(exp(-tau_n3)+exp(-r_n3-tau_n3))*
     +        (sumw3xnT1-nn*mu_n3*sumw3)/nn

!     #####
      dqdtau_c1 = -nc*sumw1/2d0+(nc-1d0)*
     +            (exp(-tau_c1)+exp(r_c1-tau_c1))*
     +            sumw1acTac/(2d0*nc)

      dqdtau_c1 = dqdtau_c1-((nc-2d0)*exp(-tau_c1)+(nc-1d0)*
     +         exp(r_c1-tau_c1)-exp(-tau_c1-r_c1))*
     +         sumw1acT1_sq/(2d0*nc**2)

      dqdtau_n1 = -nn*sumw1/2d0+(nn-1d0)*
     +        (exp(-tau_n1)+exp(r_n1-tau_n1))*sumw1anTan/(2d0*nn)

      dqdtau_n1 = dqdtau_n1-((nn-2d0)*exp(-tau_n1)+(nn-1d0)*
     +            exp(r_n1-tau_n1)-exp(-tau_n1-r_n1))*
     +            sumw1anT1_sq/(2d0*nn**2)

      dqdtau_2 = -n*sumw2/2d0+(n-1d0)*(exp(-tau_2)+exp(r_2-tau_2))*
     +        sumw2aTa/(2d0*n)

      dqdtau_2 = dqdtau_2 - ((n-2d0)*exp(-tau_2)+(n-1d0)*exp(r_2-tau_2)-
     +        exp(-tau_2-r_2))*sumw2aT1_sq/(2d0*n**2)


      dqdtau_c3 = -nc*sumw3/2d0+(nc-1d0)*(exp(-tau_c3)+
     +        exp(r_c3-tau_c3))*sumw3acTac/(2d0*nc)

      dqdtau_c3 = dqdtau_c3 - ((nc-2d0)*exp(-tau_c3)+(nc-1d0)*
     +            exp(r_c3-tau_c3)-exp(-tau_c3-r_c3))*
     +            sumw3acT1_sq/(2d0*nc**2)

      dqdtau_n3 = -nn*sumw3/2d0+(nn-1d0)*
     +            (exp(-tau_n3)+exp(r_n3-tau_n3))*
     +            sumw3anTan/(2d0*nn)

      dqdtau_n3 = dqdtau_n3 - ((nn-2d0)*exp(-tau_n3)+(nn-1d0)*
     +            exp(r_n3-tau_n3)-exp(-tau_n3-r_n3))*
     +            sumw3anT1_sq/(2d0*nn**2)


!     #####
      tt = exp(r_c1)
      if (tt .lt. 1.0d308) then
            dqdr_c1 = (nc*exp(r_c1)/(2d0*(1d0+exp(r_c1)))-0.5d0)*sumw1
      else
            dqdr_c1 = (nc/(2d0*(1d0+exp(-r_c1)))-0.5d0)*sumw1
      endif
      dqdr_c1 = dqdr_c1-(nc-1d0)*exp(r_c1-tau_c1)*sumw1acTac/(2d0*nc)
      dqdr_c1 = dqdr_c1+(exp(r_c1-tau_c1)*(nc-1d0)+exp(-r_c1-tau_c1))*
     +        sumw1acT1_sq/(2d0*nc**2)

!     #####
      tt = exp(r_n1)
      if (tt .lt. 1.0d308) then
            dqdr_n1 = (nn*exp(r_n1)/(2d0*(1+exp(r_n1)))-0.5d0)*sumw1
      else
            dqdr_n1 = (nn/(2d0*(exp(-r_n1)+1d0))-0.5d0)*sumw1
      endif
      dqdr_n1 = dqdr_n1-(nn-1d0)*exp(r_n1-tau_n1)*sumw1anTan/(2d0*nn)

      dqdr_n1 = dqdr_n1+(exp(r_n1-tau_n1)*(nn-1d0)+exp(-r_n1-tau_n1))*
     +        sumw1anT1_sq/(2d0*nn**2)


!     #####
      tt = exp(r_2)
      if (tt .lt. 1.0d308) then
            dqdr_2 = (n*exp(r_2)/(2d0*(1d0+exp(r_2)))-0.5d0)*sumw2
      else
            dqdr_2 = (n/(2d0*(exp(-r_2)+1d0))-0.5d0)*sumw2
      endif
      dqdr_2 = dqdr_2-(n-1d0)*exp(r_2-tau_2)*sumw2aTa/(2d0*n)
      dqdr_2 = dqdr_2+(exp(r_2-tau_2)*(n-1d0)+exp(-r_2-tau_2))*
     +        sumw2aT1_sq/(2d0*n**2)

!     #####
      tt = exp(r_c3)
      if (tt .lt. 1.0d308) then
            dqdr_c3 = (nc*exp(r_c3)/(2d0*(1d0+exp(r_c3)))-0.5d0)*sumw3
      else
            dqdr_c3 = (nc/(2d0*(1d0+exp(-r_c3)))-0.5d0)*sumw3
      endif
      dqdr_c3 = dqdr_c3-(nc-1d0)*exp(r_c3-tau_c3)*sumw3acTac/(2d0*nc)
      dqdr_c3 = dqdr_c3+(exp(r_c3-tau_c3)*(nc-1d0)+exp(-r_c3-tau_c3))*
     +        sumw3acT1_sq/(2d0*nc**2)

!     #####
      tt = exp(r_n3)
      if (tt .lt. 1.0d308) then
            dqdr_n3 = (nn*exp(r_n3)/(2d0*(1d0+exp(r_n3)))-0.5d0)*sumw3
      else
            dqdr_n3 = (nn/(2d0*(1d0+exp(-r_n3)))-0.5d0)*sumw3
      endif
      dqdr_n3 = dqdr_n3-(nn-1d0)*exp(r_n3-tau_n3)*sumw3anTan/(2*nn)
      dqdr_n3 = dqdr_n3+(exp(r_n3-tau_n3)*(nn-1d0)+exp(-r_n3-tau_n3))*
     +        sumw3anT1_sq/(2d0*nn**2)

!     #####

      dnegq(1) = dqdmu_c1
      dnegq(2) = dqdtau_c1
      dnegq(3) = dqdr_c1
      dnegq(4) = dqddelta_n1
      dnegq(5) = dqdtau_n1
      dnegq(6) = dqdr_n1
      dnegq(7) = dqdmu_2
      dnegq(8) = dqdtau_2
      dnegq(9) = dqdr_2
      dnegq(10)= dqdmu_c3
      dnegq(11)= dqdtau_c3
      dnegq(12)= dqdr_c3
      dnegq(13)= dqddelta_n3
      dnegq(14)= dqdtau_n3
      dnegq(15)= dqdr_n3

!     no NaN or infinity values
      do 100 i=1,15
            if ((abs(dnegq(i)) .gt. 1.0d308) .or. risnan(dnegq(i))) then
                 dnegq(i) = 0d0
            else
                 dnegq(i) = -dnegq(i)
            endif
  100 continue


      return
      end




