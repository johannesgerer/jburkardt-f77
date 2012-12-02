      program main

c*********************************************************************72
c
cc MAIN is the main program for PITCON66_PRB3.
c
c  A second order boundary value problem with a parameter:
c
c    Y''+LAMBDA*EXP(Y)=0
c    Y(0)=0  Y'(1)=0
c
c  We use a finite difference approximation, with 21 equally spaced 
c  nodes.  The value of Y at node (I) is X(I), for I=1 to 21.
c  X(22) is the value of LAMBDA.
c
c  We expect a limit point in LAMBDA at roughly LAMBDA=0.878.
c
c  NVAR-1 is the number of grid points between 0 and 1.  This problem may
c  be set up with arbitrarily many grid points, simply by increasing the
c  value of NVAR.  No other change is required.
c
c
c  This problem is solved six times:
c
c    With full storage user jacobian,
c    With full storage forward difference jacobian,
c    With full storage central difference jacobian.
c    With band storage user jacobian,
c    With band storage forward difference jacobian,
c    With band storage central difference jacobian.
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

      integer liw
      integer lrw
      integer nvar

      parameter (nvar=22)
      parameter (lrw=29+(nvar+6)*nvar)
      parameter (liw=nvar+29)

      external banslv
      external denslv
      external fpband
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
      integer jac
      character name*12
      integer nit
      double precision rwork(lrw)
      double precision xr(nvar)

      write ( *, '(a)' ) ' '
      call timestamp ( )
      write(*,*)' '
      write(*,*)'PCPRB3'
      write ( *, '(a)' ) '  FORTRAN77 version'
      write(*,*)'  PITCON sample program.'
      write(*,*)' '
      write(*,*)'  Solve a two point boundary value problem with a'
      write(*,*)'  parameter lambda, and seek a limit point in lambda.'
      write(*,*)' '
      write(*,*)'  This problem will be run six times.'
      write(*,*)' '

      itry=0
10    continue

      itry=itry+1
      write(*,*)' '
      write(*,*)'This is run number ',itry
      write(*,*)' '

      do i=1,liw
        iwork(i)=0
      enddo
 
      do i=1,lrw
        rwork(i)=0.0
      enddo
 
      ierror=0
 
      nit=0
c
c  For runs 1 and 4, supply Jacobian.
c  For runs 2 and 5, use forward difference approximation,
c  For runs 3 and 6, use central difference approximation.
c
      if(mod(itry,3).eq.1)jac=0
      if(mod(itry,3).eq.2)jac=1
      if(mod(itry,3).eq.0)jac=2
c
c  IWORK(1)=0    ; This is a startup
c  IWORK(2)=NVAR ; Use X(NVAR) for initial parameter
c  IWORK(3)=0    ; Program may choose parameter index
c  IWORK(4)=1    ; Cut down evaluations of jacobian on newton steps.
c  IWORK(5)=NVAR ; Seek target values for X(NVAR)
c  IWORK(6)=NVAR ; Seek limit points in X(NVAR)
c  IWORK(7)=1    ; Control amount of output.
c  IWORK(9)=*    ; Jacobian choice.  0=user, 1=forward, 2=central.
c
      iwork(1)=0
      iwork(2)=nvar
      iwork(3)=0
      iwork(4)=1
      iwork(5)=nvar
      iwork(6)=nvar
      iwork(7)=1
      iwork(9)=jac
c
c  Pass the lower and upper bandwidths of the Jacobian
c  in the integer parameter array IPAR.
c
      ipar(1)=1
      ipar(2)=1
c
c  RWORK(1)=0.0001 ; Absolute error tolerance
c  RWORK(2)=0.0001 ; Relative error tolerance
c  RWORK(3)=0.01   ; Minimum stepsize
c  RWORK(4)=0.25   ; Maximum stepsize
c  RWORK(5)=0.05   ; Starting stepsize
c  RWORK(6)=1.0    ; Starting direction
c  RWORK(7)=0.80   ; Target value (Seek solution with X(NVAR)=0.80)
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
 
      if(itry.le.3)then
        write(*,*)'This run uses the full linear solver DENSLV.'
      else
        write(*,*)'This run uses the banded linear solver BANSLV.'
      endif

      if(jac.eq.0)then
        write(*,*)'This run assumes that the user supplies the'
        write(*,*)'jacobian matrix via a subroutine.'
      elseif(jac.eq.1)then
        write(*,*)'This run assumes that PITCON will approximate'
        write(*,*)'the jacobian using forward differences.'
      elseif(jac.eq.2)then
        write(*,*)'This run assumes that PITCON will approximate'
        write(*,*)'the jacobian using central differences.'
      endif
 
      write(*,*)' '
      write(*,*)'Step  Type of point     Lambda'
      write(*,*)' '
      i=0
      name='Start point  '
      write(*,'(1x,i3,2x,a12,2x,g14.6)')i,name,xr(nvar)
c
c  Take steps along the curve
c
      do i=1,50
 
        if(itry.le.3)then
          call pitcon(fpfull,fpar,fxbend,ierror,ipar,iwork,liw,
     &      nvar,rwork,lrw,xr,denslv)
        else
          call pitcon(fpband,fpar,fxbend,ierror,ipar,iwork,liw,
     &      nvar,rwork,lrw,xr,banslv)
        endif
 
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
 
        write(*,'(1x,i3,2x,a12,2x,g14.6)')i,name,xr(nvar)
 
        if(iwork(1).eq.3)then
          write(*,*)' '
          nit=nit+1
          write(*,*)'Value of target point ',nit
          write(*,'(1x,5g14.6)')(xr(j),j=1,nvar)
          write(*,*)' '
          if(nit.ge.2)go to 60
        endif
 
        if(ierror.ne.0)then
          write(*,*)'An error occurred, ierror=',ierror
          write(*,*)'We terminate this portion of the computation.'
          go to 60
        endif
 
      enddo
 
60    continue
      write(*,*)' '
      write(*,*)'Jacobians:      ',iwork(19)
      write(*,*)'Factorizations: ',iwork(20)
      write(*,*)'Solves:         ',iwork(21)
      write(*,*)'Functions:      ',iwork(22)
      if(itry.lt.6)go to 10

      write ( *, '(a)' ) ' '
      call timestamp ( )
      stop
      end
      subroutine fxbend(nvar,fpar,ipar,x,fx,ierror)

c*********************************************************************72
c
cc FXBEND evaluates the function.
c
c  Discussion:
c
c    Nonlinear function evaluation for two point BVP.
c
c    Equation 1:  (boundary condition)
c
c    Y(0)=0  becomes FX(1)=X(1)
c
c    Equations 2 through NVAR-2 are the differential equation:
c
c    Y'' + LAMBDA*EXP(Y) = 0  becomes
c
c    FX(I)= X(I-1) -2*X(I) + X(I+1)
c          +X(NVAR)*EXP(X(I))*HSQ = 0
c
c    where HSQ is the square of the distance between nodes.
c
c    Equation NVAR-1:  (boundary condition)
c
c    Y'(1)=0  becomes  FX(NVAR-1)=X(NVAR-1)-X(NVAR-2)=0
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
      implicit none

      intrinsic exp

      integer nvar

      double precision fpar(*)
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
cc FPFULL evaluates the jacobian matrix using full storage.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
      implicit none

      intrinsic exp

      integer nvar

      double precision fpar(*)
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
      subroutine fpband(nvar,fpar,ipar,x,fprime,ierror)

c*********************************************************************72
c
cc FPBAND evaluates the jacobian matrix using band storage.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
      implicit none

      intrinsic exp

      integer nvar

      double precision fpar(*)
      double precision fprime(*)
      double precision hsq
      integer i
      integer ierror
      integer indx
      integer ipar(2)
      integer ml
      integer mu
      integer nband
      double precision x(nvar)

      hsq=1.0/(nvar-2)**2

      ml=ipar(1)
      mu=ipar(2)
      nband=2*ml+mu+1
      indx=1+1*(nband-1)-ml
      fprime(indx)=1.0
 
      do i=2,nvar-2
        indx=i+(i-1)*(nband-1)-ml
        fprime(indx)=1.0
        indx=i+i*(nband-1)-ml
        fprime(indx)=-2.0+x(nvar)*exp(x(i))*hsq
        indx=i+(i+1)*(nband-1)-ml
        fprime(indx)=1.0
        indx=(nvar-1)*nband+i
        fprime(indx)=exp(x(i))*hsq
      enddo
 
      indx=nvar-1+(nvar-2)*(nband-1)-ml
      fprime(indx)=-1.0
      indx=nvar-1+(nvar-1)*(nband-1)-ml
      fprime(indx)=1.0

      return
      end
