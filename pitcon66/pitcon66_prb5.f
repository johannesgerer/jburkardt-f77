      program main

c*********************************************************************72
c
cc MAIN is the main program for PITCON66_PRB5.
c
c  Discussion:
c
c    Discretized two point boundary value problem with parameter LAMBDA:
c
c      Y''+LAMBDA*EXP(Y)=0
c      Y(0)=0 Y'(1)=0
c
c    This test illustrates the use of 
c
c    the full Newton method, which updates the Jacobian on every 
c    Newton iteration,
c
c    the modified Newton method, which updates the Jacobian only at
c    specified steps in the Newton iteration,
c
c    the "cheap" Newton method, which updates the Jacobian only 
c    when convergence fails.
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
c    Herbert Keller,
c    Numerical Methods for Two-point Boundary Value Problems,
c    Dover, 1992,
c    ISBN: 0486669254,
c    LC: QA372.K42.
c
      implicit none

      integer nvar
      parameter (nvar=22)

      integer lrw
      parameter (lrw=29+(nvar+6)*nvar)

      integer liw
      parameter (liw=nvar+29)

      external denslv
      external fpfull
      external fxbend
      external pitcon

      double precision fpar(1)
      integer i
      integer ierror
      integer ipar(2)
      integer itry
      integer iwork(liw)
      integer j
      integer modcon
      character name*12
      integer nit
      double precision rwork(lrw)
      double precision xr(nvar)

      itry=0

      write ( *, '(a)' ) ' '
      call timestamp ( )
      write(*,*)' '
      write(*,*)'PCPRB5'
      write ( *, '(a)' ) '  FORTRAN77 version'
      write(*,*)'  PITCON sample program.'
      write(*,*)'  Two point BVP, limit point in Lambda'
      write(*,*)' '
      write(*,*)'This program will be run three times.'

10    continue

      itry=itry+1
      write(*,*)' '
      write(*,*)'This is run number ',itry
      nit=0

      write(*,*)' '

      do i=1,liw
        iwork(i)=0
      enddo

      do i=1,lrw
        rwork(i)=0.0
      enddo

      ierror=0
c
c  Set input quantities
c
      if(itry.eq.1)modcon=0
      if(itry.eq.2)modcon=1
      if(itry.eq.3)modcon=2
c
c  IWORK(1)=0      ; This is a startup
c  IWORK(2)=NVAR   ; Use X(NVAR) for initial parameter
c  IWORK(3)=0      ; Program may choose parameter index
c  IWORK(4)=MODCON ; Control frequency of jacobian update in newton iteration
c  IWORK(5)=NVAR   ; Seek target values for X(NVAR)
c  IWORK(6)=NVAR   ; Seek limit points in X(NVAR)
c  IWORK(7)=1      ; Control amount of output.
c  IWORK(9)=0      ; Jacobian choice.
c
      iwork(1)=0
      iwork(2)=nvar
      iwork(3)=0
      iwork(4)=modcon
      iwork(5)=nvar
      iwork(6)=nvar
      iwork(7)=1
      iwork(9)=0
c
c  Jacobian upper/lower bandwidths:
c
      ipar(1)=1
      ipar(2)=1
c
c  RWORK(1)=0.0001; Absolute error tolerance
c  RWORK(2)=0.0001; Relative error tolerance
c  RWORK(3)=0.01   ; Minimum stepsize
c  RWORK(4)=0.25   ; Maximum stepsize
c  RWORK(5)=0.05   ; Starting stepsize
c  RWORK(6)=1.0    ; Starting direction
c  RWORK(7)=0.80   ; Target value (Seek solution with X(NVAR)=0.80
c
      rwork(1)=0.0001
      rwork(2)=0.0001
      rwork(3)=0.01
      rwork(4)=0.25
      rwork(5)=0.05
      rwork(6)=1.0
      rwork(7)=0.80
c
c  Set starting point
c
      do i=1,nvar
        xr(i)=0.0
      enddo

      if(itry.eq.1)then
        write(*,*)' '
        write(*,*)'Using the full Newton method.'
        write(*,*)'The jacobian is updated on every iteration.'
      elseif(itry.eq.2)then
        write(*,*)' '
        write(*,*)'Using the modified Newton method.'
        write(*,*)'The jacobian is held fixed while correcting a '
        write(*,*)'point.'
      else
        write(*,*)' '
        write(*,*)'Using the "cheap" Newton method.'
        write(*,*)'The jacobian is held fixed as long as possible,'
        write(*,*)'perhaps over multiple points, until'
        write(*,*)'convergence fails.'
      endif

      write(*,*)' '
      write(*,*)'Step  Type of point     Lambda'
      write(*,*)' '

      i=0
      name='Start point  '
      write(*,'(1x,i3,2x,a12,2x,g14.6)')i,name,xr(nvar)

      do i=1,50

        call pitcon(fpfull,fpar,fxbend,ierror,ipar,iwork,liw,
     &    nvar,rwork,lrw,xr,denslv)

        if(iwork(1).eq.1)then
          name='Corrected    '
        elseif(iwork(1).eq.2)then
          name='Continuation '
        elseif(iwork(1).eq.3)then
          name='Target point '
        elseif(iwork(1).eq.4)then
          name='Limit point  '
        elseif(iwork(1).lt.0)then
          name='Jacobian   '
        endif

        write(*,'(1x,i3,2x,a12,2x,3g14.6)')i,name,xr(nvar)

        if(iwork(1).eq.3)then
          nit=nit+1
          write(*,*)' '
          write(*,*)'Complete value for target point ',nit
          write(*,*)' '
          write(*,'(1x,5g14.6)')(xr(j),j=1,nvar)
          if(nit.ge.2)go to 60
        endif

        if(ierror.ne.0)then
          write(*,*)'Iteration halted, IERROR=',ierror
          go to 60
        endif

      enddo

60    continue

      write(*,*)' '
      write(*,*)'Jacobians      ',iwork(19)
      write(*,*)'Factorizations ',iwork(20)
      write(*,*)'Solves         ',iwork(21)
      write(*,*)'Functions      ',iwork(22)

      if(itry.lt.3)go to 10

      write ( *, '(a)' ) ' '
      call timestamp ( )
      stop
      end
      subroutine fxbend(nvar,fpar,ipar,x,fx,ierror)

c*********************************************************************72
c
cc FXBEND evaluates the function.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
      implicit none

      intrinsic exp

      integer nvar

      double precision fpar(1)
      double precision fx(nvar)
      double precision hsq
      integer i
      integer ierror
      integer ipar(2)
      double precision x(nvar)

      hsq=1.0/(nvar-2)**2

      fx(1)=x(1)

      do i=2,nvar-2
        fx(i)=x(i-1)-2.0*x(i)+x(i+1)+x(nvar)*exp(x(i))*hsq
      enddo

      fx(nvar-1)=x(nvar-1)-x(nvar-2)

      return
      end
      subroutine fpfull(nvar,fpar,ipar,x,fprime,ierror)

c*********************************************************************72
c
cc FPFULL evaluates the jacobian.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
      implicit none

      intrinsic exp

      integer nvar

      double precision fpar(1)
      double precision fprime(nvar,nvar)
      double precision hsq
      integer i
      integer ierror
      integer ipar(2)
      double precision x(nvar)

      hsq=1.0/(nvar-2)**2

      fprime(1,1)=1.0

      do i=2,nvar-2
        fprime(i,i-1)=1.0
        fprime(i,i)=-2.0+x(nvar)*exp(x(i))*hsq
        fprime(i,i+1)=1.0
        fprime(i,nvar)=exp(x(i))*hsq
      enddo

      fprime(nvar-1,nvar-2)=-1.0
      fprime(nvar-1,nvar-1)=1.0

      return
      end
