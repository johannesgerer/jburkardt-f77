      program main

c*********************************************************************72
c
cc MAIN is the main program for PITCON66_PRB6.
c
c  Discussion:
c
c    Demonstrate the options available to approximate the jacobian,
c    and to check its accuracy.
c
c    RWORK(18) is used to determine the size of the finite difference
c    increments used to approximate the jacobian.  The continuation
c    program will set a default value for this, but the user may 
c    override it.
c
c    IWORK(1) can be used to check the user supplied jacobian routine
c    against a forward or central difference approximations, and to
c    specify whether the maximum value of the difference, or the 
c    entire matrix of differences is to be printed.
c
c    The problem solved is the Freudenstein-Roth function, although 
c    the solution procedure is only carried out for five steps.  In 
c    this example, we're more interested in the output of the 
c    jacobian checker than in the solution of the problem.
c
c    On the second try of the problem, the Jacobian is "polluted".
c    The jacobian checker points this out.  Note, however, that the
c    program is able to compute points on the curve, even with a bad
c    Jacobian.  This is typical.  The good news is that small errors 
c    in the jacobian aren't fatal.  The bad news is that you lose 
c    quadratic convergence in Newton's method, and may not realize 
c    that something's wrong, since you still get answers!
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
      implicit none

      integer liw
      integer lrw
      integer nvar

      parameter (nvar=3)
      parameter (liw=nvar+29)
      parameter (lrw=29+(6+nvar)*nvar)

      external denslv
      external dfroth
      external fxroth
      external pitcon

      double precision fpar(1)
      integer i
      integer ierror
      integer ipar(1)
      integer isave
      integer itry
      integer iwork(liw)
      integer j
      character name*12
      double precision rwork(lrw)
      double precision xr(nvar)

      write ( *, '(a)' ) ' '
      call timestamp ( )
      write(*,*)' '
      write(*,*)'PCPRB6'
      write ( *, '(a)' ) '  FORTRAN77 version'
      write(*,*)'  PITCON sample program.'
      write(*,*)'  Freudenstein-Roth function.'
      write(*,*)' '
      write(*,*)'This test demonstrates the use of IWORK(1)'
      write(*,*)'and RWORK(18) to approximate the jacobian,'
      write(*,*)'compare the user jacobian to an approximation,'
      write(*,*)'choose forward or centered differences,'
      write(*,*)'choose the size of the difference increment,'
      write(*,*)'print user, approximate jacobian or difference,'
      write(*,*)'print full matrix, or maximum entry.'
      write(*,*)' '

      itry=0
10    continue
      itry=itry+1

      write(*,*)' '

      if(itry.eq.1)then
        write(*,*)'Run 1: standard run for comparison.'
        write(*,*)'Don''t call jacobian checker.'
      elseif(itry.eq.2)then
        write(*,*)'Run 2: run with bad jacobian.'
        write(*,*)'Call jacobian checker at third step.'
        write(*,*)'Use check value IWORK(1)=-1.'
      elseif(itry.eq.3)then
        write(*,*)'Run 3: run with bad jacobian.'
        write(*,*)'Call jacobian checker at third step.'
        write(*,*)'Use check value iwork(1)=-2.'
      elseif(itry.eq.4)then
        write(*,*)'Run 4: run with bad jacobian.'
        write(*,*)'Call jacobian checker at third step.'
        write(*,*)'Use check value iwork(1)=-3.'
      elseif(itry.eq.5)then
        write(*,*)'Run 5: run with bad jacobian.'
        write(*,*)'Call jacobian checker at third step.'
        write(*,*)'Use check value iwork(1)=-4.'
      elseif(itry.eq.6)then
        write(*,*)'Run 6: run with good jacobian.'
        write(*,*)'Call jacobian checker at third step.'
        write(*,*)'Use check value iwork(1)=-4.'
        write(*,*)'Use default value of rwork(18)'
      elseif(itry.eq.7)then
        write(*,*)'Run 7: run with good jacobian.'
        write(*,*)'Call jacobian checker at third step.'
        write(*,*)'Use check value iwork(1)=-4.'
        write(*,*)'Use finite difference increment rwork(18)=0.1'
      elseif(itry.eq.8)then
        write(*,*)'Run 8: run with good jacobian.'
        write(*,*)'Call jacobian checker at third step.'
        write(*,*)'Use check value iwork(1)=-4.'
        write(*,*)'Use finite difference increment rwork(18)=0.01'
      elseif(itry.eq.9)then
        write(*,*)'Run 9: run with good jacobian.'
        write(*,*)'Call jacobian checker at third step.'
        write(*,*)'Use check value iwork(1)=-4.'
        write(*,*)'Finite difference increment rwork(18)=0.0001'
      elseif(itry.eq.10)then
        write(*,*)'Run 10: run with good jacobian.'
        write(*,*)'Call jacobian checker at third step.'
        write(*,*)'Use check value iwork(1)=-5.'
      elseif(itry.eq.11)then
        write(*,*)'Run 11: run with good jacobian.'
        write(*,*)'Call jacobian checker at third step.'
        write(*,*)'Use check value iwork(1)=-6.'
      elseif(itry.eq.12)then
        write(*,*)'Run 12: run with good jacobian.'
        write(*,*)'Call jacobian checker at third step.'
        write(*,*)'Use check value iwork(1)=-7.'
      elseif(itry.eq.13)then
        write(*,*)'Run 13: run with good jacobian.'
        write(*,*)'Call jacobian checker at third step.'
        write(*,*)'Use check value iwork(1)=-8.'
      elseif(itry.eq.14)then
        write(*,*)'Run 14: run with good jacobian.'
        write(*,*)'Call jacobian checker at third step.'
        write(*,*)'Use check value iwork(1)=-9.'
      elseif(itry.eq.15)then
        write(*,*)'Run 15: run with good jacobian.'
        write(*,*)'Call jacobian checker at third step.'
        write(*,*)'Use check value iwork(1)=-10.'
      endif

      write(*,*)' '
c
c  Use IPAR(1) as a signal to Jacobian routine.  If IPAR(1)=1,
c  then use a 'bad' jacobian.
c
      if(itry.eq.1.or.itry.ge.6)then
        ipar(1)=0
      else
        ipar(1)=1
      endif

      do i=1,liw
        iwork(i)=0
      enddo

      do i=1,lrw
        rwork(i)=0.0
      enddo
c
c  IWORK(1)=0 ; This is a startup
c  IWORK(2)=2 ; Use index 2 for first parameter
c  IWORK(3)=0 ; Program may choose index
c  IWORK(4)=0 ; Update jacobian every newton step
c  IWORK(5)=3 ; Seek target values for index 3
c  IWORK(6)=1 ; Seek limit points in index 1
c  IWORK(7)=0 ; small amount of output
c  IWORK(9)=0 ; Use user's jacobian routine
c
      iwork(1)=0
      iwork(2)=2
      iwork(3)=0
      iwork(4)=0
      iwork(5)=3
      iwork(6)=1
      iwork(7)=0
      iwork(9)=0
c
c  RWORK(1)=0.00001; Absolute error tolerance
c  RWORK(2)=0.0001 ; Relative error tolerance
c  RWORK(3)=0.01   ; Minimum stepsize
c  RWORK(4)=20.0   ; Maximum stepsize
c  RWORK(5)=0.3    ; Starting stepsize
c  RWORK(6)=1.0    ; Starting direction
c  RWORK(7)=1.0    ; Target value (Seek solution with X(3)=1)
c
      rwork(1)=0.00001
      rwork(2)=0.0001
      rwork(3)=0.01
      rwork(4)=20.0
      rwork(5)=0.3
      rwork(6)=1.0
      rwork(7)=1.0
c
c  Set the parameter that determines the size of the increment used for
c  finite difference approximations.
c
      if(itry.eq.7)then
        rwork(18)=0.1
      elseif(itry.eq.8)then
        rwork(18)=0.01
      elseif(itry.eq.9)then
        rwork(18)=0.0001
      else
        rwork(18)=0.0
      endif
c
c  Set starting point.
c
      xr(1)=15.0
      xr(2)=-2.0
      xr(3)=0.0

      write(*,*)' '
      write(*,*)'Step  Type of point     '//
     &  'X(1)         X(2)         X(3)'
      write(*,*)' '
      i=0
      name='Start point  '
      write(*,'(1x,i3,2x,a12,2x,3g14.6)')i,name,(xr(j),j=1,nvar)
c
c  Take another step.  This loop normally runs for 40 steps, but we'll cut it
c  to just 5.
c
      do 40 i=1,5
c
c  Check jacobian on third step
c
        if(i.eq.3.and.itry.gt.1)then

          isave=iwork(1)

          if(itry.eq.2)then
            iwork(1)=-1
          elseif(itry.eq.3)then
            iwork(1)=-2
          elseif(itry.eq.4)then
            iwork(1)=-3
          elseif(itry.eq.5)then
            iwork(1)=-4
          elseif(itry.eq.10)then
            iwork(1)=-5
          elseif(itry.eq.11)then
            iwork(1)=-6
          elseif(itry.eq.12)then
            iwork(1)=-7
          elseif(itry.eq.13)then
            iwork(1)=-8
          elseif(itry.eq.14)then
            iwork(1)=-9
          elseif(itry.eq.15)then
            iwork(1)=-10
          else
            iwork(1)=-4
          endif

          write(*,*)' '
          write(*,*)'Turning on jacobian check option!'
          write(*,*)' '
        endif

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
          name='Jacobian'
        endif
c
c  After Jacobian check, restore value of IWORK(1).
c
        if(i.eq.3.and.itry.gt.1)then
          iwork(1)=isave
        endif

        write(*,'(1x,i3,2x,a12,2x,3g14.6)')i,name,(xr(j),j=1,nvar)

        if(iwork(1).eq.3)then
          write(*,*)' '
          write(*,*)'We have reached the point we wanted.'
          go to 50
        endif

        if(ierror.ne.0)then
          write(*,*)' '
          write(*,*)'An error occurred.'
          write(*,*)'We will terminate the computation now.'
          go to 50
        endif

        if(ierror.ne.0)then
          write(*,*)'Iteration halted, ierror=',ierror
          go to 50
        endif

40      continue

50    continue
      if(itry.lt.15)go to 10
      write(*,*)' '
      write(*,*)'Normal conclusion of tests.'
      write ( *, '(a)' ) ' '
      call timestamp ( )
      stop
      end
      subroutine fxroth(nvar,fpar,ipar,x,f,ierror)

c*********************************************************************72
c
cc FXROTH evaluates the function.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
      implicit none

      integer nvar

      double precision f(nvar)
      double precision fpar(1)
      integer ierror
      integer ipar(1)
      double precision x(nvar)

      f(1)=x(1)-((x(2)-5.0)*x(2)+2.0)*x(2)-13.0+34.0*(x(3)-1.0)

      f(2)=x(1)+((x(2)+1.0)*x(2)-14.0)*x(2)-29.0+10.0*(x(3)-1.0)

      return
      end
      subroutine dfroth(nvar,fpar,ipar,x,fjac,ierror)

c*********************************************************************72
c
cc DFROTH evaluates the jacobian matrix.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
      implicit none

      integer nvar

      double precision fjac(nvar,nvar)
      double precision fpar(1)
      integer ierror
      integer ipar(1)
      double precision x(nvar)

      fjac(1,1)=1.0
      fjac(1,2)=(-3.0*x(2)+10.0)*x(2)-2.0
      fjac(1,3)=34.0

      fjac(2,1)=1.0
      fjac(2,2)=(3.0*x(2)+2.0)*x(2)-14.0
      fjac(2,3)=10.0
c
c  If IPAR(1)=1, put some incorrect values in jacobian
c
      if(ipar(1).eq.1)then
        fjac(1,1)=fjac(1,1)+0.125
        fjac(1,2)=fjac(1,2)+0.250
        fjac(2,2)=fjac(2,2)+0.500
      endif

      return
      end
