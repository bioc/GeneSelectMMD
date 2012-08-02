
!     the code in the loop from paraEst.R converted to Fortran

      subroutine paraestloop(w1, w2, w3, x, psi_m, memsubjects,
     +                  xctxc, xct1,
     +                  xntxn, xnt1,
     +                  xtx, xt1,
     +                  ngenes, n, nc, nn,
     +                  itmax, eps, quiet, loop)

      implicit none

      integer itmax, loop
      integer ngenes, n, nc, nn
      integer memsubjects(n)

      double precision x(ngenes,n), psi_m(17)
      double precision xctxc(ngenes), xct1(ngenes)
      double precision xntxn(ngenes), xnt1(ngenes)
      double precision xtx(ngenes), xt1(ngenes)
      double precision eps
      logical quiet

      double precision w1(ngenes), w2(ngenes), w3(ngenes)
      double precision sumw1, sumw2, sumw3

      double precision sumw1xctxc, sumw1xct1, sumw1xct1_sq
      double precision sumw1xntxn, sumw1xnt1, sumw1xnt1_sq
      double precision sumw2xtx, sumw2xt1, sumw2xt1_sq
      double precision sumw3xctxc, sumw3xct1, sumw3xct1_sq
      double precision sumw3xntxn, sumw3xnt1, sumw3xnt1_sq

      double precision const_b1, const_b2, const_b3
      double precision sumb
      double precision pivec_new(2)
      double precision psi_m_new(17)

      double precision xct1_sq(ngenes), xnt1_sq(ngenes), xt1_sq(ngenes)
      double precision sumtwo, sumd

      double precision thetaini(15)

      character*60 task
      integer isave(44)
      double precision dsave(29)

      integer i
      logical done

      done = .TRUE.

      const_b1 = 2d0
      const_b2 = 2d0
      const_b3 = 2d0
      sumb = const_b1 + const_b2 + const_b3


      do 10 i=1,ngenes
        xct1_sq(i) = xct1(i)**2
        xnt1_sq(i) = xnt1(i)**2
        xt1_sq(i) = xt1(i)**2
 10   continue


!     loop until convergence and stop if max number of loops is reached
      do 100 loop=0,itmax

!       check for user interrupt
        call rchkusr()

        if(.not. quiet) then
          call intpr('********** loop=', -1, loop, 1)
!          write(6,*)'************'
!          write(6,*)'loop=',loop
        end if


!       E-step
!       obtain weights w_{ij}(\Psi^{(m)}), i.e. E(Z_{ij} | x_i, Psi^{(m)})

!       checkpara was called within wiFun.R, however this will all be handled
!       in Fortran, so checking the dimension of Psi.m will not be included
!       in the future

        call wifun(w1, w2, w3, x, psi_m, memsubjects, 
     +            1.0d-6, ngenes, n, nc, nn       )


        sumw1 = sumd(w1,ngenes)
        sumw2 = sumd(w2,ngenes)
        sumw3 = sumd(w3,ngenes)

!        if (.not. quiet) then
!          write(6,*)'sumw1=',sumw1,' sumw2=',sumw2,' sumw3=',sumw3
!        end if


        sumw1xctxc = sumtwo(w1, xctxc, ngenes)
        sumw1xct1 = sumtwo(w1, xct1, ngenes)
        sumw1xct1_sq = sumtwo(w1, xct1_sq, ngenes)

        sumw1xntxn = sumtwo(w1, xntxn, ngenes)
        sumw1xnt1 = sumtwo(w1, xnt1, ngenes)
        sumw1xnt1_sq = sumtwo(w1, xnt1_sq, ngenes)

        sumw2xtx = sumtwo(w2, xtx, ngenes)
        sumw2xt1 = sumtwo(w2, xt1, ngenes)
        sumw2xt1_sq = sumtwo(w2, xt1_sq, ngenes)

        sumw3xctxc = sumtwo(w3, xctxc, ngenes)
        sumw3xct1 = sumtwo(w3, xct1, ngenes)
        sumw3xct1_sq = sumtwo(w3, xct1_sq, ngenes)
 
        sumw3xntxn = sumtwo(w3, xntxn, ngenes)
        sumw3xnt1 = sumtwo(w3, xnt1, ngenes)
        sumw3xnt1_sq = sumtwo(w3, xnt1_sq, ngenes)


!       M-step
!       mixture proportions
        pivec_new(1) = (sumw1+const_b1-1d0)/(ngenes+sumb-3d0)
        pivec_new(2) = (sumw2+const_b2-1d0)/(ngenes+sumb-3d0)


        do 20 i=1,15
          thetaini(i) = psi_m(i+2)
 20     continue

!        if(.not. quiet) then
!          write(6,*)
!          write(6,*)'thetaIni>>'
!          write(6,30)'     muc1= ',thetaini(1),
!     +               '   tau.c1= ',thetaini(2),
!     +               '     r.c1= ',thetaini(3),
!     +               ' delta.n1= ',thetaini(4),
!     +               '   tau.n1= ',thetaini(5)
!          write(6,30)'     r.n1= ',thetaini(6),
!     +               '     mu.2= ',thetaini(7),
!     +               '    tau.2= ',thetaini(8),
!     +               '      r.2= ',thetaini(9),
!     +               '    mu.c3= ',thetaini(10)
!          write(6,30)'   tau.c3= ',thetaini(11),
!     +               '     r.c3= ',thetaini(12),
!     +               ' delta.n3= ',thetaini(13),
!     +               '   tau.n3= ',thetaini(14),
!     +               '     r.n3= ',thetaini(15)
!        end if
! 30     format(A,D10.4,A,D10.4,A,D10.4,A,D10.4,A,D10.4)

        
        call lbfgsbdriver(task, isave, dsave, thetaini,
     +                  nc, nn, n,
     +                  sumw1, sumw2, sumw3,
     +                  sumw1xctxc, sumw1xct1, sumw1xct1_sq,
     +                  sumw1xntxn, sumw1xnt1, sumw1xnt1_sq,
     +                  sumw2xtx, sumw2xt1, sumw2xt1_sq,
     +                  sumw3xctxc, sumw3xct1, sumw3xct1_sq,
     +                  sumw3xntxn, sumw3xnt1, sumw3xnt1_sq       )

!        if(.not. quiet) then
!          write(6,*)
!          if(task(1:4) .eq. 'CONV') then
!            write(6,*)'message= ',task
!            write(6,*)'number of calls to fn and gr>> ',isave(34)
!          end if
!          write(6,*)
!        end if
!!        negQ<-res$dsave[2]


!       # update parameters
        psi_m_new(1) = pivec_new(1)
        psi_m_new(2) = pivec_new(2)
        do 40 i=1,15
          psi_m_new(i+2) = thetaini(i)
 40     continue

!       names(Psi.m.new)<-paraNamesRP

!        if(.not. quiet) then
!          write(6,*)'Psi.m.new>>'
!          write(6,50)'      pi1= ',psi_m_new(1),
!     +               '      pi2= ',psi_m_new(2)
!          write(6,55)'     muc1= ',psi_m_new(3),
!     +               '   tau.c1= ',psi_m_new(4),
!     +               '     r.c1= ',psi_m_new(5),
!     +               ' delta.n1= ',psi_m_new(6),
!     +               '   tau.n1= ',psi_m_new(7)
!          write(6,55)'     r.n1= ',psi_m_new(8),
!     +               '     mu.2= ',psi_m_new(9),
!     +               '    tau.2= ',psi_m_new(10),
!     +               '      r.2= ',psi_m_new(11),
!     +               '    mu.c3= ',psi_m_new(12)
!          write(6,55)'   tau.c3= ',psi_m_new(13),
!     +               '     r.c3= ',psi_m_new(14),
!     +               ' delta.n3= ',psi_m_new(15),
!     +               '   tau.n3= ',psi_m_new(16),
!     +               '     r.n3= ',psi_m_new(17)
!          write(6,*)
!        end if
! 50     format(A,D10.4,A,D10.4)
! 55     format(A,D10.4,A,D10.4,A,D10.4,A,D10.4,A,D10.4)


        done = .TRUE.

        do 60 i=1,17
          if(abs(psi_m_new(i)-psi_m(i)) .gt. eps) done = .FALSE.
 60     continue

        if(done) then
!         # the following statement is to 
!         # make sure eta.new is an update from current wMat
!         # so that d Q / d eta = 0
          do 70 i=1,17
            psi_m(i) = psi_m_new(i)
 70       continue
          return
        end if


        do 80 i=1,17
          psi_m(i) = psi_m_new(i)
 80     continue


 100  continue

      if(loop .ge. itmax) then
        call rwarn('***** Warning! ITMAX exceeded *****')
        call rexit('EM algorithm did not converge!')
      end if

      return
      end




! ************************************************************************************
!     function to sum an array of size n

      function sumd(arr1, n)

      implicit none

      integer n
      double precision sumd
      double precision arr1(n)
      integer i

      sumd = 0d0

      do 100 i=1,n
        sumd = sumd + arr1(i)
 100  continue

      return
      end



! ************************************************************************************
!     function to sum an array that is the element by element
!     product of two arrays (arr1, arr2) 
      function sumtwo(arr1, arr2, n)

      implicit none

      integer n
      double precision sumtwo, sumd
      double precision arr1(n), arr2(n), arr3(n)
      integer i

      do 100 i=1,n
        arr3(i) = arr1(i)*arr2(i)
 100  continue

      sumtwo=sumd(arr3,n)

      return
      end



! ************************************************************************************
! find the position of the maximum element of memmat

      subroutine maxposfun(memgenes, memmat, 
     +                  ngenes, maxflag, thrspostprob)

      implicit none

      integer ngenes, memgenes(ngenes)
      double precision memmat(ngenes,3)
      logical maxflag
      double precision thrspostprob

      integer i

      if (maxflag) then
        do 100 i=1,ngenes
          if ((memmat(i,1) .gt. memmat(i,2)) .and. 
     +                   (memmat(i,1) .gt. memmat(i,3))) then
            memgenes(i) = 1
          else if ((memmat(i,2) .gt. memmat(i,1)) .and. 
     +                   (memmat(i,2) .gt. memmat(i,3))) then
            memgenes(i) = 2
          else if ((memmat(i,3) .gt. memmat(i,1)) .and. 
     +                   (memmat(i,3) .gt. memmat(i,2))) then
            memgenes(i) = 3
          end if
 100    continue
        return
      end if

      do 200 i=1,ngenes
          if ((memmat(i,1) .le. thrspostprob) .and.
     +      (memmat(i,2) .le. thrspostprob) .and.
     +      (memmat(i,3) .le. thrspostprob)   ) then
!           non-differentially expressed
            memgenes(i) = 2
          end if
 200  continue
      return

      end


