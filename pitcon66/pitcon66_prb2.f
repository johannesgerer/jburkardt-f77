      program main

c*********************************************************************72
c
cc MAIN is the main program for the aircraft stability problem.
c
c  The variables X(I) are control parameters for an aircraft under
c  the guidance of a pilot:
c
c
c  X(1) = Roll
c  X(2) = Pitch
c  X(3) = Yaw
c  X(4) = Angle of attack
c  X(5) = Sideslip
c  X(6) = Elevator
c  X(7) = Aileron
c  X(8) = Rudder
c
c
c  The function F is of the following form:
c
c
c  For I=1 to 5,
c
c  F(I)=Sum(J=1,8) Sum(K=1,8) B(I,J)*X(J)+PHI(I,J,K)*X(J)*X(K)
c
c  F(6)=X(IFIX1)-VAL1
c  F(7)=X(IFIX2)-VAL2
c
c
c  Options:
c
c
c  ICHOOZ=1  IFIX1=6, VAL1=-.050, Rheinboldt BARRAY.
c  ICHOOZ=2  IFIX1=6, VAL1=-.008, Rheinboldt BARRAY.
c  ICHOOZ=3  IFIX1=6, VAL1=0.000, Rheinboldt BARRAY.
c  ICHOOZ=4  IFIX1=6, VAL1=0.050, Rheinboldt BARRAY.
c  ICHOOZ=5  IFIX1=6, VAL1=0.100, Rheinboldt BARRAY.
c  ICHOOZ=6  IFIX1=6, VAL1=-0.01250, Rheinboldt BARRAY.  (Bifurcation)
c
c  ICHOOZ=11  IFIX1=6, VAL1=-.050, Melhem BARRAY.
c  ICHOOZ=12  IFIX1=6, VAL1=-.008, Melhem BARRAY.
c  ICHOOZ=13  IFIX1=6, VAL1=0.000, Melhem BARRAY.
c  ICHOOZ=14  IFIX1=6, VAL1=0.050, Melhem BARRAY.
c  ICHOOZ=15  IFIX1=6, VAL1=0.100, Melhem BARRAY.
c  ICHOOZ=16  IFIX1=6, VAL1=-0.01250, Melhem BARRAY.  (Bifurcation)
c
c  For all choices of ICHOOZ, IFIX2=8, VAL2=0.0.
c  If ICHOOZ=6 or 16, LIM=0, otherwise LIM=7.
c
c
c  Limit points:
c
c
c  Melhem lists the following limit points in X7.  Note that Melhem
c  has BARRAY(4,1)=1.0, BARRAY(4,2)=0.0:
c
c
c     X1       X2        X3        X4        X5      X6     X7      X8
c
c  -2.9691  0.83074  -0.072748  0.41029  -0.26880  -.05  0.50919    0.0
c  -2.8158 -0.17481  -0.089469  0.026319  0.070951 -.008 0.20442    0.0
c  -3.7571 -0.64911  -0.39350   0.091810  0.19685  -.008 -.0038238  0.0
c  -4.1637  0.092284 -0.092610  0.022402 -0.017106 -.008 0.37823    0.0
c  -2.5839 -0.22128  -0.054079  0.013524  0.090871 0.0   0.18608    0.0
c  -3.9007 -1.1421   -0.57863   0.13284   0.32685  0.0  -0.50703    0.0
c  -2.3610 -0.72360   0.032739 -0.039108  0.29347  0.05  0.29272    0.0
c  -2.2982  1.4033    0.063244 -0.079383  0.58336  0.10  0.58336    0.0
c
c
c  Rheinboldt lists the following limit points, using
c  BARRAY(4,1)=0.0, BARRAY(4,2)=1.0:
c
c
c     X1       X2        X3        X4        X5      X6     X7      X8
c
c   2.9648  0.82556   0.07366   0.041309  0.26734 -0.050 -.050481   0.0
c   2.8173 -0.17628   0.08992   0.026429 -0.07147 -0.008 -.204973   0.0
c   3.7579 -0.65541   0.38658   0.092520 -0.19867 -0.008 0.006200   0.0
c   4.1638  0.08913   0.09480   0.022888  0.16232 -0.008 -.377660   0.0
c   2.5873 -0.22354   0.05468   0.013676 -0.09168  0.000 -.186908   0.0
c   3.9005 -1.14815   0.58156   0.133516 -0.32858   0.000  .510158  0.0
c   2.3639 -0.72974  -0.31604  -0.038785 -0.29583   0.050 -.295772  0.0
c   2.2992 -1.41023  -0.06184  -0.079009 -0.58629   0.100 -.689717  0.0
c
c
c  Bifurcation points:
c
c
c  Rheinboldt, using BARRAY(4,1)=0.0, BARRAY(4,2)=1.0, lists:
c
c
c   4.482   0.1632    0.02373   0.006205  0.03527 -0.0006177 -0.3986 0.0
c   3.319  -0.1869    0.1605    0.04379  -0.06888 -0.01250   -0.2374 0.0
c   4.466   0.1467    0.04045   0.009777  0.03089 -0.006129  -0.3995 0.0
c  -3.325   0.1880   -0.1614    0.04395   0.06911 -0.01247    0.2367 0.0
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    22 August 2011
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Raman Mehra, William Kessel, James Carroll,
c    Global stability and contral analysis of aircraft at high angles of attack,
c    Technical Report CR-215-248-1, -2, -3,
c    Office of Naval Research, June 1977.
c
c    Albert Schy, Margery Hannah,
c    Prediction of Jump Phenomena in Roll-coupled Maneuvers
c    of Airplanes,
c    Journal of Aircraft,
c    Volume 14, 1977,  pages 375-382.
c
c    John Young, Albert Schy, Katherine Johnson,
c    Prediction of Jump Phenomena in Aircraft Maneuvers, Including 
c    Nonlinear Aerodynamic Effects,
c    Journal of Guidance and Control,
c    Volume 1, 1978, pages 26-31.
c
      implicit none

      integer liw
      integer lrw
      integer nvar

      parameter (nvar=8)
      parameter (liw=29+nvar)
      parameter (lrw=29+nvar*(nvar+6))

      external denslv
      external fpair
      external fxair
      external pitcon

      double precision barray(5,8)
      double precision hmax
      integer i
      integer ichooz
      integer ierror
      integer ifix1
      integer ifix2
      integer ipar(2)
      integer iwork(liw)
      integer j
      integer k
      integer lim
      character name*12
      double precision rpar(42)
      double precision rwork(lrw)
      double precision val1
      double precision val2
      double precision xr(nvar)
c
c  Set options.
c
      ichooz=1
      ifix1=6
      ifix2=8
c
c  Initialize work arrays.
c
      do i=1,liw
        iwork(i)=0
      enddo
 
      do i=1,lrw
        rwork(i)=0.0
      enddo
c
c  Set VAL1, VAL2 based on ICHOOZ
c
      if(ichooz.eq.1.or.ichooz.eq.11)val1=-.050
      if(ichooz.eq.2.or.ichooz.eq.12)val1=-.008
      if(ichooz.eq.3.or.ichooz.eq.13)val1=0.000
      if(ichooz.eq.4.or.ichooz.eq.14)val1=0.050
      if(ichooz.eq.5.or.ichooz.eq.15)val1=0.100
      if(ichooz.eq.6.or.ichooz.eq.16)val1=-0.01250
      val2=0.0
c
c  Set some parameters
c
      ipar(1)=ifix1
      ipar(2)=ifix2
      if(ichooz.ne.6.and.ichooz.ne.16)then
        hmax=0.4
      else
        hmax=0.2
        endif
      if(ichooz.ne.6.and.ichooz.ne.16)then
        lim=7
      else
        lim=0
      endif
c
c  Set starting point
c
      do i=1,nvar
        xr(i)=0.0
      enddo
 
      xr(ifix1)=val1
      xr(ifix2)=val2
c
c  Set BARRAY
c
      barray(1,1)=-3.933
      barray(2,1)=0.0
      barray(3,1)=0.002
      if(ichooz.lt.10)then
        barray(4,1)=0.0
      else
        barray(4,1)=1.0
      endif
      barray(5,1)=0.0
      barray(1,2)=0.107
      barray(2,2)=-0.987
      barray(3,2)=0.0
      if(ichooz.lt.10)then
        barray(4,2)=1.0
      else
        barray(4,2)=0.0
      endif
      barray(5,2)=0.0
      barray(1,3)=0.126
      barray(2,3)=0.0
      barray(3,3)=-0.235
      barray(4,3)=0.0
      barray(5,3)=-1.0
      barray(1,4)=0.0
      barray(2,4)=-22.95
      barray(3,4)=0.0
      barray(4,4)=-1.0
      barray(5,4)=0.0
      barray(1,5)=-9.99
      barray(2,5)=0.0
      barray(3,5)=5.67
      barray(4,5)=0.0
      barray(5,5)=-0.196
      barray(1,6)=0.0
      barray(2,6)=-28.37
      barray(3,6)=0.0
      barray(4,6)=-0.168
      barray(5,6)=0.0
      barray(1,7)=-45.83
      barray(2,7)=0.0
      barray(3,7)=-0.921
      barray(4,7)=0.0
      barray(5,7)=-0.0071
      barray(1,8)=-7.64
      barray(2,8)=0.0
      barray(3,8)=-6.51
      barray(4,8)=0.0
      barray(5,8)=0.0
c
c  Copy BARRAY, VAL1, VAL2 into RPAR.
c
      k=0
      do i=1,5
        do j=1,nvar
          k=k+1
          rpar(k)=barray(i,j)
        enddo
      enddo
      k=k+1
      rpar(k)=val1
      k=k+1
      rpar(k)=val2
c
c  Set work array input
c
c  IWORK(1)=0 ; This is a startup
c  IWORK(2)=7 ; Use X(7) for initial parameter
c  IWORK(3)=0 ; Program may choose parameter index
c  IWORK(4)=0 ; Update jacobian every newton step
c  IWORK(5)=0 ; Seek no target values.
c  IWORK(6)=LIM ; Seek limit points in X(LIM)
c  IWORK(7)=1 ; Control amount of output.
c  IWORK(9)=0 ; Jacobian choice.
c
      iwork(1)=0
      iwork(2)=7
      iwork(3)=0
      iwork(4)=0
      iwork(5)=0
      iwork(6)=lim
      iwork(7)=1
      iwork(9)=0
c
c  RWORK(1)=0.0001 ; Absolute error tolerance
c  RWORK(2)=0.0001 ; Relative error tolerance
c  RWORK(3)=0.0001 ; Minimum stepsize
c  RWORK(4)=HMAX   ; Maximum stepsize
c  RWORK(5)=0.1    ; Starting stepsize
c  RWORK(6)=-1.0   ; Starting direction
c  RWORK(7)=0.0    ; Target value
c
      rwork(1)=0.0001
      rwork(2)=0.0001
      rwork(3)=0.0001
      rwork(4)=hmax
      rwork(5)=0.1
      rwork(6)=-1.0
      rwork(7)=0.0

      write ( *, '(a)' ) ' '
      call timestamp ( )
      write(*,*)' '
      write(*,*)'PITCON_PRB2:'
      write ( *, '(a)' ) '  FORTRAN77 version'
      write(*,*)'  pitcon sample program.'
      write(*,*)' '
      write(*,*)'  the aircraft stability problem'
      write(*,*)'  using option ',ichooz
      write(*,*)'  fix variable ',ifix1,' at ',val1
      write(*,*)'  fix variable ',ifix2,' at ',val2
      if(ichooz.lt.10)then
        write(*,*)'  rheinboldt version of barray'
      else
        write(*,*)'  melhem version of barray'
      endif
      write(*,*)' '
      write(*,*)'step  type of point     '
      write(*,*)' '
      i=0
      name='start point  '
      write(*,'(1x,i3,2x,a12,2x,8f7.3)')i,name,(xr(j),j=1,nvar)
 
      do i=1,50
        call pitcon(fpair,rpar,fxair,ierror,ipar,iwork,liw,
     &  nvar,rwork,lrw,xr,denslv)
 
        if(iwork(1).eq.1)then
          name='corrected    '
        elseif(iwork(1).eq.2)then
          name='continuation '
        elseif(iwork(1).eq.3)then
          name='target point '
        elseif(iwork(1).eq.4)then
          name='limit point  '
        elseif(iwork(1).lt.0)then
          name='jacobian   '
        endif
 
        write(*,'(1x,i3,2x,a12,2x,8f7.3)')i,name,(xr(j),j=1,nvar)
 
        if(ierror.ne.0)then
          write(*,*)'iteration halted, ierror=',ierror
          go to 70
        endif
 
      enddo
 
70    continue
 
      write(*,*)' '
      write(*,*)'jacobians:     ',iwork(19)
      write(*,*)'factorizations:',iwork(20)
      write(*,*)'solves:        ',iwork(21)
      write(*,*)'functions:     ',iwork(22)
      write ( *, '(a)' ) ' '
      call timestamp ( )
      stop
      end
      subroutine fxair(nvar,rpar,ipar,x,fx,ierror)

c*********************************************************************72
c
cc FXAIR evaluates the function.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
      integer nvar

      double precision fx(nvar)
      integer i
      integer ierror
      integer ifix1
      integer ifix2
      integer ipar(2)
      integer j
      integer k
      double precision phi
      double precision rpar(42)
      double precision val1
      double precision val2
      double precision x(nvar)
c
c  Compute linear terms
c
      k=0
      do i=1,5
        fx(i)=0.0
        do j=1,8
          k=k+1
          fx(i)=fx(i)+rpar(k)*x(j)
        enddo
      enddo
c
c  Compute nonlinear terms
c
      phi=-.727*x(2)*x(3)+8.39*x(3)*x(4)-684.4*x(4)*x(5)+63.5*x(4)*x(7)
      fx(1)=fx(1)+phi
 
      phi=.949*x(1)*x(3)+.173*x(1)*x(5)
      fx(2)=fx(2)+phi
 
      phi=-.716*x(1)*x(2)-1.578*x(1)*x(4)+1.132*x(4)*x(7)
      fx(3)=fx(3)+phi
 
      phi=-x(1)*x(5)
      fx(4)=fx(4)+phi
 
      phi=x(1)*x(4)
      fx(5)=fx(5)+phi
c
c  Two function values restrict two variables:
c
      ifix1=ipar(1)
      val1=rpar(41)
      fx(6)=x(ifix1)-val1
 
      ifix2=ipar(2)
      val2=rpar(42)
      fx(7)=x(ifix2)-val2
 
      return
      end
      subroutine fpair(nvar,rpar,ipar,x,fprime,ierror)

c*********************************************************************72
c
cc FPAIR evaluates the jacobian matrix.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
      implicit none

      integer nvar

      double precision fprime(nvar,nvar)
      integer i
      integer ierror
      integer ifix1
      integer ifix2
      integer ipar(2)
      integer j
      integer k
      double precision rpar(42)
      double precision x(nvar)
c
c  Linear part of F
c
      k=0
      do i=1,5
        do j=1,8
          k=k+1
          fprime(i,j)=rpar(k)
        enddo
      enddo
 
      do i=6,7
        do j=1,8
          fprime(i,j)=0.0
        enddo
      enddo
 
      ifix1=ipar(1)
      ifix2=ipar(2)
      fprime(6,ifix1)=1.0
      fprime(7,ifix2)=1.0
c
c  Nonlinear part of F
c
      fprime(1,2)=fprime(1,2)-.727*x(3)
      fprime(1,3)=fprime(1,3)-.727*x(2)+8.39*x(4)
      fprime(1,4)=fprime(1,4)+8.39*x(3)-684.4*x(5)+63.5*x(7)
      fprime(1,5)=fprime(1,5)-684.4*x(4)
      fprime(1,7)=fprime(1,7)+63.5*x(4)
 
      fprime(2,1)=fprime(2,1)+.949*x(3)+.173*x(5)
      fprime(2,3)=fprime(2,3)+.949*x(1)
      fprime(2,5)=fprime(2,5)+.173*x(1)
 
      fprime(3,1)=fprime(3,1)-.716*x(2)-1.578*x(4)
      fprime(3,2)=fprime(3,2)-.716*x(1)
      fprime(3,4)=fprime(3,4)-1.578*x(1)+1.132*x(7)
      fprime(3,7)=fprime(3,7)+1.132*x(4)
 
      fprime(4,1)=fprime(4,1)-x(5)
      fprime(4,5)=fprime(4,5)-x(1)
 
      fprime(5,1)=fprime(5,1)+x(4)
      fprime(5,4)=fprime(5,4)+x(1)
 
      return
      end
