! converted to Fortran July, 2012
!
! modified on Sept. 28, 2011
!  (1) added 'na.rm=TRUE' to function 'sum'
!
! get wij=pi.j*f.j(xi)/[pi.1*f.1(xi)+pi.2*f.2(xi)+pi.3*f.3(xi)]
! j=1,2,3
! Psi.m=(pi.1, pi.2, 
!        mu.c1, tau.c1, r.c1, 
!        delta.n1, tau.n1, r.n1, 
!        mu.2, tau.2, r.2, 
!        mu.c3, tau.c3, r.c3, 
!        delta.n3, tau.n3, r.n3) 

      subroutine wifun(wi1arr, wi2arr, wi3arr,
     +                  x, psi_m, memsubjects, eps,
     +                  ngenes, n, nc, nn)

      implicit none

      double precision wi1arr(ngenes), wi2arr(ngenes), wi3arr(ngenes)
      double precision x(ngenes,n), psi_m(17)
      double precision eps

      integer memsubjects(n), ngenes, n, nc, nn
      integer i, j, ci, ni
      integer pos

      character*85 errmess

      double precision pi_1, pi_2, pi_3
      double precision wi1, wi2, wi3

      double precision mu_c1, tau_c1, r_c1, rho_c1, delta_n1, mu_n1
      double precision tau_n1, r_n1, rho_n1, mu_2, tau_2, r_2, rho_2
      double precision mu_c3, tau_c3, r_c3, rho_c3, delta_n3, mu_n3
      double precision tau_n3, r_n3, rho_n3
      double precision sigma2_c1, sigma2_n1, sigma2_2
      double precision sigma2_c3, sigma2_n3

      double precision xci(nc), xni(nn), xi(n)
      double precision xci_sq(nc), xni_sq(nn), xi_sq(n)

      double precision xiTxi_c, xiT1_c, aiTai_c1, aiT12_c1
      double precision xiTxi_n, xiT1_n, aiTai_n1, aiT12_n1
      double precision part_c1, part_n1, delta1, log_detSigma1
      double precision xiTxi, xiT1, aiTai_2, aiT12_2
      double precision delta2, log_detSigma2

      double precision aiTai_c3, aiT12_c3, aiTai_n3, aiT12_n3
      double precision part_c3, part_n3, delta3, log_detSigma3

      double precision logf2df1, f2df1, logf3df1, f3df1, denom_1
      double precision logf1df2, f1df2, logf3df2, f3df2, denom_2
      double precision logf1df3, f1df3, logf2df3, f2df3, denom_3


      pos = 1
!      eps = 1.0d-6

!     check if parameters are in appropriate ranges
!     not checking dimension of psi_m though
!     this was previously handled by the R routine checkPara

!     mixture proportions
      pi_1 = psi_m(1)
      pi_2 = psi_m(2)
      pi_3 = 1d0 - pi_1 - pi_2

!     if parameters are out of the allowed range then exit with error message
      if((pi_1 .ge. 1d0) .or. (pi_1 .le. 0d0)) then
            write(errmess,'(A,D10.4,A,D10.4,A,D10.4,A)')
     +                        'pi.1=',pi_1,' pi.2=',pi_2,' pi.3=',pi_3,
     +                        '    pi.1 should be in (0, 1)          '
            call rexit(errmess)
      else if((pi_2 .ge. 1d0) .or. (pi_2 .le. 0d0)) then
            write(errmess,'(A,D10.4,A,D10.4,A,D10.4,A)')
     +                        'pi.1=',pi_1,' pi.2=',pi_2,' pi.3=',pi_3,
     +                        '    pi.2 should be in (0, 1)          '
            call rexit(errmess)
      else if((pi_3 .ge. 1d0) .or. (pi_3 .le. 0d0)) then
            write(errmess,'(A,D10.4,A,D10.4,A,D10.4,A)')
     +                        'pi.1=',pi_1,' pi.2=',pi_2,' pi.3=',pi_3,
     +                        '    pi.3 should be in (0, 1)          '
            call rexit(errmess)
      else if(abs((pi_1+pi_2+pi_3)-1d0) .gt. eps) then
            write(errmess,'(A,D10.4,A,D10.4,A,D10.4,A)')
     +                        'pi.1=',pi_1,' pi.2=',pi_2,' pi.3=',pi_3,
     +                        '   pi.1+pi.2+pi.3 should be equal to 1'
            call rexit(errmess)
      end if


!     mixture proportions calculated above

!     mean expression level for cluster 1 for diseased subjects
      mu_c1 = psi_m(3)
!     variance of expression levels for cluster 1 for diseased subjects
      tau_c1 = psi_m(4)
      sigma2_c1 = exp(tau_c1)
!     modified logit of correlation among expression levels for cluster 1 for diseased subjects
      r_c1 = psi_m(5)
      rho_c1 = (exp(r_c1)-1d0/(nc-1d0))/(1d0+exp(r_c1))

!     mu.n1=mu.c1-exp(delta.n1)
!      =mean expression level for cluster 1 for normal subjects
      delta_n1 = psi_m(6)
      mu_n1 = mu_c1-exp(delta_n1)
!     variance of expression levels for cluster 1 for normal subjects
      tau_n1 = psi_m(7)
      sigma2_n1 = exp(tau_n1)
!     modified logit of correlation among expression levels for cluster 1 for normal subjects
      r_n1 = psi_m(8)
      rho_n1 = (exp(r_n1)-1d0/(nn-1d0))/(1d0+exp(r_n1))

!     mean expression level for cluster 2
      mu_2 = psi_m(9)
!     variance of expression levels for cluster 2
      tau_2 = psi_m(10)
      sigma2_2 = exp(tau_2)
!     modified logit of correlation among expression levels for cluster 2
      r_2 = psi_m(11)
      rho_2 = (exp(r_2)-1d0/(n-1d0))/(1d0+exp(r_2))

!     mean expression level for cluster 3 for diseased subjects
      mu_c3 = psi_m(12)
!     variance of expression levels for cluster 3 for diseased subjects
      tau_c3 = psi_m(13)
      sigma2_c3 = exp(tau_c3)
!     modified logit of correlation among expression levels for cluster 3 for diseased subjects
      r_c3 = psi_m(14)
      rho_c3 = (exp(r_c3)-1d0/(nc-1d0))/(1d0+exp(r_c3))

!     mu.n3=mu.c3+exp(delta.n3)
!     =mean expression level for cluster 3 for normal subjects
      delta_n3 = psi_m(15)
      mu_n3 = mu_c3 + exp(delta_n3)
!     variance of expression levels for cluster 3 for normal subjects
      tau_n3 = psi_m(16)
      sigma2_n3 = exp(tau_n3)
!     modified logit of correlation among expression levels for cluster 3 for normal subjects
      r_n3 = psi_m(17)
      rho_n3 = (exp(r_n3)-1d0/(nn-1d0))/(1d0+exp(r_n3))


!     if values are out of the allowed range then exit with error message
      if(nc .ge. n) then
            errmess = 'Number of cases >= Total number of patients!'
            call rexit(errmess)
      else if(nn .ge. n) then
            errmess = 'Number of controls >= Total number of patients!'
            call rexit(errmess)
      end if



!     loop over ngenes to find each wi value
      do 100 j=1,ngenes

            ci = 1
            ni = 1
            do 20 i=1,n
                  xi(i) = x(j,i)
                  xi_sq(i) = x(j,i)**2
!                 expression levels of the gene for diseased subjects
                  if (memsubjects(i) .eq. 1) then
                        xci(ci) = x(j,i)
                        xci_sq(ci) = x(j,i)**2
                        ci = ci + 1
!                 expression levels of the gene for non-diseased subjects
                  else if (memsubjects(i) .eq. 0) then
                        xni(ni) = x(j,i)
                        xni_sq(ni) = x(j,i)**2
                        ni = ni + 1
                  end if
 20         continue



!           ###
!           # density for genes in cluster 1 (over-expressed)
!           ###

            xiTxi_c = sum(xci_sq)
            xiT1_c = sum(xci)
            aiTai_c1 = xiTxi_c-2d0*mu_c1*xiT1_c+nc*mu_c1**2
            aiT12_c1 = xiT1_c**2-2d0*nc*mu_c1*xiT1_c+nc**2*mu_c1**2

            xiTxi_n = sum(xni_sq)
            xiT1_n = sum(xni)
            aiTai_n1 = xiTxi_n-2d0*mu_n1*xiT1_n+nn*mu_n1**2
            aiT12_n1 = xiT1_n**2-2d0*nn*mu_n1*xiT1_n+nn**2*mu_n1**2

            part_c1 = (aiTai_c1-rho_c1*aiT12_c1/(1d0+(nc-1d0)*rho_c1))
     +                  / (sigma2_c1 * (1d0-rho_c1)) 
            part_n1 = (aiTai_n1-rho_n1*aiT12_n1/(1d0+(nn-1d0)*rho_n1))
     +                  / (sigma2_n1 * (1d0-rho_n1)) 
            delta1 = (part_c1+part_n1)

            log_detSigma1 = nc*log(sigma2_c1)+
     +            ((nc-1d0)*log(1d0-rho_c1)+log(1d0+(nc-1d0)*rho_c1))+
     +                  nn*log(sigma2_n1)+
     +            ((nn-1d0)*log(1d0-rho_n1)+log(1d0+(nn-1d0)*rho_n1))


!           ###
!           # density for genes in cluster 2 (non-expressed)
!           ###
            xiTxi = sum(xi_sq)
            xiT1 = sum(xi)
            aiTai_2 = xiTxi-2d0*mu_2*xiT1+n*mu_2**2
            aiT12_2 = xiT1**2-2d0*n*mu_2*xiT1+n**2*mu_2**2

            delta2 = ( aiTai_2 - rho_2* aiT12_2 / (1d0+(n-1d0)*rho_2) )
     +                  / (sigma2_2 * (1d0-rho_2)) 

            log_detSigma2 = n*log(sigma2_2)+
     +                  ((n-1d0)*log(1d0-rho_2)+log(1d0+(n-1d0)*rho_2))


!           ###
!           # density for genes in cluster 3 (under-expressed)
!           ###

            aiTai_c3 = xiTxi_c-2d0*mu_c3*xiT1_c+nc*mu_c3**2
            aiT12_c3 = xiT1_c**2-2d0*nc*mu_c3*xiT1_c+nc**2*mu_c3**2

            aiTai_n3 = xiTxi_n-2d0*mu_n3*xiT1_n+nn*mu_n3**2
            aiT12_n3 = xiT1_n**2-2d0*nn*mu_n3*xiT1_n+nn**2*mu_n3**2

            part_c3 = (aiTai_c3-rho_c3*aiT12_c3/(1d0+(nc-1d0)*rho_c3))
     +                  / (sigma2_c3 * (1d0-rho_c3)) 
            part_n3 = (aiTai_n3-rho_n3*aiT12_n3/(1d0+(nn-1d0)*rho_n3))
     +                  / (sigma2_n3 * (1d0-rho_n3)) 

            delta3 = (part_c3+part_n3)

            log_detSigma3 = nc*log(sigma2_c3)+
     +            (nc-1d0)*log(1d0-rho_c3)+log(1d0+(nc-1d0)*rho_c3)+
     +                  nn*log(sigma2_n3)+
     +            (nn-1d0)*log(1d0-rho_n3)+log(1d0+(nn-1d0)*rho_n3)


            if ((delta1 .lt. delta2) .and. (delta1 .lt. delta3)) then
                  pos = 1
            else if ((delta2 .lt. delta1).and.(delta2 .lt. delta3)) then
                  pos = 2
            else if ((delta3 .lt. delta2).and.(delta3 .lt. delta1)) then
                  pos = 3
            end if

!           ################


            if(pos .eq. 1) then
!                 # wi1=pi.1*f1/(pi.1*f1+pi.2*f2+pi.3*f3)
!                 #    = pi.1/(pi.1+pi.2*f2/f1 + pi.3 * f3/f1)
!                 # wi2=pi.2*f2/f1/[pi.1+pi.2*f2/f1 + pi.3 * f3/f1]
!                 # wi3=pi.3*f3/f1/[pi.1+pi.2*f2/f1 + pi.3 * f3/f1]
!                 # f2/f1
!                 # log(f2/f1)=[log(|Sigma1|)-log(|Sigma2|)+delta1-delta2]/2
                  logf2df1 = (log_detSigma1-
     +                        log_detSigma2+delta1-delta2)/2d0
                  f2df1 = exp(logf2df1)
    
!                 # f3/f1
!                 # log(f3/f1)=[log(|Sigma1|)-log(|Sigma3|)+delta1-delta3]/2
                  logf3df1 = (log_detSigma1-
     +                        log_detSigma3+delta1-delta3)/2d0
                  f3df1 = exp(logf3df1)

                  denom_1 = (pi_1+pi_2*f2df1+pi_3*f3df1)
  
                  wi1 = pi_1/denom_1
                  wi2 = pi_2*f2df1/denom_1
                  wi3 = pi_3*f3df1/denom_1

            else if(pos .eq. 2) then
!                 #wi2=pi.2*f2/(pi.1*f1+pi.2*f2+pi.3*f3)
!                 #   = pi.2/(pi.1*f1/f2+pi.2 + pi.3 * f3/f2)
!                 #wi1=pi.1 * f1/f2 / [pi.1*f1/f2+pi.2 + pi.3 * f3/f2]
!                 #wi3=pi.3 * f3/f2 / [pi.1*f1/f2+pi.2 + pi.3 * f3/f2]
!                 # f1/f2
!                 # log(f1/f2)=(log|Sigma2|-log|Sigma1|+detla2-delta1)/2
                  logf1df2 = (log_detSigma2-
     +                        log_detSigma1+delta2-delta1)/2d0
                  f1df2 = exp(logf1df2)
    
!                 # f3/f2
!                 # log(f3/f2)=(log|Sigma2|-log|Sigma3|+delta2-detla3)/2
                  logf3df2 = (log_detSigma2-
     +                        log_detSigma3+delta2-delta3)/2d0
                  f3df2 = exp(logf3df2)

                  denom_2 = (pi_1*f1df2+pi_2+pi_3*f3df2)
  
                  wi1 = pi_1*f1df2/denom_2
                  wi2 = pi_2/denom_2
                  wi3 = pi_3*f3df2/denom_2

            else
!                 #wi3=pi.3*f3/(pi.1*f1+pi.2*f2+pi.3*f3)
!                 #   = pi.3/(pi.1*f1/f3+pi.2*f2/f3 + pi.3)
!                 #wi1=pi.1*f1/f3 / [pi.1*f1/f3+pi.2*f2/f3 + pi.3]
!                 #wi2=pi.2*f2/f3 / [pi.1*f1/f3+pi.2*f2/f3 + pi.3]
!                 # f1/f3
!                 # log(f1/f3)=(log|Sigma3|-log|Sigma1|+delta3-delta1)/2
                  logf1df3 = (log_detSigma3-
     +                        log_detSigma1+delta3-delta1)/2d0
                  f1df3 = exp(logf1df3)
 
!                 # f2/f3
!                 # log(f2/f3)=(log|Sigma3|-log|Sigma2|+delta3-delta2)/2
                  logf2df3 = (log_detSigma3-
     +                        log_detSigma2+delta3-delta2)/2d0
                  f2df3 = exp(logf2df3)
  
                  denom_3 = (pi_1*f1df3+pi_2*f2df3+pi_3)

                  wi1 = pi_1*f1df3/denom_3
                  wi2 = pi_2*f2df3/denom_3
                  wi3 = pi_3/denom_3
            end if


            wi1arr(j) = wi1
            wi2arr(j) = wi2
            wi3arr(j) = wi3

!     end loop over the genes (rows of x)
100   continue

      return
      end


