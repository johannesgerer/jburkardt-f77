      program main

c*********************************************************************72
c
cc MAIN is the main program for PITCON66_PRB1.
c
c  Discussion:
c
c    PITCON66_PRB1 treats a system based on the Freudenstein-Roth function.
c
c    The function F(X) is of the form
c
c      FX(1) = X1 - X2**3 + 5*X2**2 -  2*X2 - 13 + 34*(X3-1)
c      FX(2) = X1 + X2**3 +   X2**2 - 14*X2 - 29 + 10*(X3-1)
c
c    Starting from the point (15,-2,0), the program is required to produce
c    solution points along the curve until it reaches a solution point
c    (*,*,1).  It also may be requested to look for limit points in the
c    first or third components.
c
c    The correct value of the solution at X3=1 is (5,4,1).
c
c    Limit points in the first variable occur at:
c
c      (14.28309, -1.741377,  0.2585779)
c      (61.66936,  1.983801, -0.6638797)
c
c    Limit points for the third variable occur at:
c
c      (20.48586, -0.8968053, 0.5875873)
c      (61.02031,  2.230139, -0.6863528)
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
c    F Freudenstein, B Roth,
c    Numerical Solutions of Nonlinear Equations,
c    Journal of the Association for Computing Machinery,
c    Volume 10, 1963, Pages 550-556.
c
      implicit none

      integer liw
      integer lrw
      integer nvar

      parameter (nvar=3)
      parameter (liw=nvar+29)
      parameter (lrw=29+(6+nvar)*nvar)

      external fxroth
      external dfroth
      external denslv
      external pitcon

      double precision fpar(1)
      integer i
      integer ierror
      integer ipar(1)
      integer iwork(liw)
      integer j
      character name*12
      double precision rwork(lrw)
      double precision xr(nvar)
c
c  Set work arrays to zero:
c
      do i=1,liw
        iwork(i)=0
      enddo

      do i=1,lrw
        rwork(i)=0.0
      enddo
c
c  Set some entries of work arrays.
c
c  IWORK(1)=0 ; This is a startup
c  IWORK(2)=2 ; Use X(2) for initial parameter
c  IWORK(3)=0 ; Program may choose parameter index
c  IWORK(4)=0 ; Update jacobian every newton step
c  IWORK(5)=3 ; Seek target values for X(3)
c  IWORK(6)=1 ; Seek limit points in X(1)
c  IWORK(7)=1 ; Control amount of output.
c  IWORK(9)=0 ; Jacobian choice.
c
      iwork(1)=0
      iwork(2)=2
      iwork(3)=0
      iwork(4)=0
      iwork(5)=3
      iwork(6)=1
      iwork(7)=1
      iwork(9)=0
c
c  RWORK(1)=0.00001; Absolute error tolerance
c  RWORK(2)=0.00001; Relative error tolerance
c  RWORK(3)=0.01   ; Minimum stepsize
c  RWORK(4)=10.0   ; Maximum stepsize
c  RWORK(5)=0.3    ; Starting stepsize
c  RWORK(6)=1.0    ; Starting direction
c  RWORK(7)=1.0    ; Target value (Seek solution with X(3)=1)
c
      rwork(1)=0.00001
      rwork(2)=0.00001
      rwork(3)=0.01
      rwork(4)=10.0
      rwork(5)=0.3
      rwork(6)=1.0
      rwork(7)=1.0
c
c  Set starting point.
c
      xr(1)=15.0
      xr(2)=-2.0
      xr(3)=0.0

      write ( *, '(a)' ) ' '
      call timestamp ( )
      write(*,*)' '
      write(*,*)'PCPRB1:'
      write ( *, '(a)' ) '  FORTRAN77 version'
      write(*,*)'  PITCON test problem'
      write(*,*)'  Freudenstein-Roth function'
      write(*,*)'  2 equations, 3 variables.'
      write(*,*)' '
      write(*,*)'Step  Type of point     '//
     & 'X(1)         X(2)         X(3)'
      write(*,*)' '
      i=0
      name='Start point '
      write(*,'(1x,i3,2x,a12,2x,3g14.6)')i,name,(xr(j),j=1,nvar)

      do i=1,30

        call pitcon(dfroth,fpar,fxroth,ierror,ipar,iwork,liw,
     &    nvar,rwork,lrw,xr,denslv)

        if(iwork(1).eq.1)then
          name='Corrected   '
        elseif(iwork(1).eq.2)then
          name='Continuation'
        elseif(iwork(1).eq.3)then
          name='Target point'
        elseif(iwork(1).eq.4)then
          name='Limit point '
        elseif(iwork(1).lt.0)then
          name='Jacobian    '
        endif

        write(*,'(1x,i3,2x,a12,2x,3g14.6)')i,name,(xr(j),j=1,nvar)

        if(iwork(1).eq.3)then
          write(*,*)' '
          write(*,*)'We have reached the point we wanted.'
          write(*,*)'The code may stop now.'
          write ( *, '(a)' ) ' '
          call timestamp ( )
          stop
        endif

        if(ierror.ne.0)then
          write(*,*)' '
          write(*,*)'An error occurred.'
          write(*,*)'We will terminate the computation now.'
          write ( *, '(a)' ) ' '
          call timestamp ( )
          stop
        endif

      enddo

      write(*,*)' '
      write(*,*)'PITCON did not reach the point of interest.'
      write ( *, '(a)' ) ' '
      call timestamp ( )
      stop
      end
      subroutine fxroth(nvar,fpar,ipar,x,f,ierror)

c*********************************************************************72
c
cc FXROTH evaluates the function.
c
c  ( X1 - ((X2-5.0)*X2+2.0)*X2 - 13.0 + 34.0*(X3-1.0)  )
c  ( X1 + ((X2+1.0)*X2-14.0)*X2 - 29.0 + 10.0*(X3-1.0) )
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
      implicit none

      integer nvar

      double precision f(*)
      double precision fpar(*)
      integer ierror
      integer ipar(*)
      double precision x(nvar)

      f(1)=x(1)
     &     -((x(2)-5.0)*x(2)+2.0)*x(2)
     &     -13.0
     &     +34.0*(x(3)-1.0)

      f(2)=x(1)
     &    +((x(2)+1.0)*x(2)-14.0)*x(2)
     &    -29.0
     &    +10.0*(x(3)-1.0)

      return
      end
      subroutine dfroth(nvar,fpar,ipar,x,fjac,ierror)

c*********************************************************************72
c
cc DFROTH evaluates the Jacobian matrix.
c
c  ( 1.0   (-3.0*X(2)+10.0)*X(2)-2.0      34.0   )
c  ( 1.0   (3.0*X(2)+2.0)*X(2)-14.0       10.0   )
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
      implicit none

      integer nvar

      double precision fjac(nvar,nvar)
      double precision fpar(*)
      integer ierror
      integer ipar(*)
      double precision x(nvar)

      fjac(1,1)=1.0
      fjac(1,2)=(-3.0*x(2)+10.0)*x(2)-2.0
      fjac(1,3)=34.0

      fjac(2,1)=1.0
      fjac(2,2)=(3.0*x(2)+2.0)*x(2)-14.0
      fjac(2,3)=10.0

      return
      end
