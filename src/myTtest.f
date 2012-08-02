
!     provides the t-test data for two matrices of data
!     results are provided for each row (gene) in the matrices
      subroutine myttest(xc, xn, nc, nn, ngenes, pvals, tstats)

      implicit none

      integer nc, nn, ngenes
      integer i, j
      double precision xc(ngenes,nc), xn(ngenes,nn)
      double precision pvals(ngenes), tstats(ngenes)

      double precision tmpxc(nc), tmpxn(nn)
      double precision sumxc, meanxc, sumxn, meanxn
      double precision varxc, varxn
      double precision sumd
      double precision tmpt, df, rpt

      integer nnc, nnn
      logical risnan


!     loop through the genes
!     exclude the NA values in the xc and xn matrices
      do 100 i=1,ngenes
         nnc = 0
         do 20 j=1,nc
            if (.not. risnan(xc(i,j))) then
              tmpxc(j) = xc(i,j)
              nnc = nnc + 1
            end if
 20      continue

         nnn = 0
         do 25 j=1,nn
            if (.not. risnan(xn(i,j))) then
              tmpxn(j) = xn(i,j)
              nnn = nnn + 1
            end if
 25      continue

!        use the nc and nn values for non-NA values (nnc, nnn)

!        find the mean and sample variance
         sumxc = sumd(tmpxc,nnc)
         meanxc = sumxc/nnc
         varxc = 0d0
         do 30 j=1,nnc
            varxc = varxc + (tmpxc(j) - meanxc)**2
 30      continue
         varxc = varxc/(nnc-1d0)

!        find the mean and sample variance
         sumxn = sumd(tmpxn,nnn)
         meanxn = sumxn/nnn
         varxn = 0d0
         do 40 j=1,nnn
            varxn = varxn + (tmpxn(j) - meanxn)**2
 40      continue
         varxn = varxn/(nnn-1d0)

!        calculate test statistic
         tmpt = (meanxc-meanxn)/sqrt((varxc/nnc)+(varxn/nnn))

!        calculate the degrees of freedom
         df = ((varxc/nnc)+(varxn/nnn))**2
     +           /((varxc/nnc)**2/(nnc-1d0)+(varxn/nnn)**2/(nnn-1d0))

!        calculate the p-value using call to c wrapper for R function pt()
         pvals(i) = 2d0*rpt(-abs(tmpt),df,1,0)
         tstats(i) = tmpt

 100   continue

      return
      end

