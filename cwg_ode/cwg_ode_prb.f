      program main

c*********************************************************************72
c
cc MAIN is the main program for CWG_ODE_PRB.
c
c  Discussion:
c
c    This is a sample calling program for the CWG_ODE code.
c
c    the ode solved here has the following form
c
c      x'' = -x, x(0)=0, x'(0)=1    (solution x(t)=sin(t)).
c
c    this must first be transformed to a pair of first order equations,
c    which results in the following system-
c
c      y1' = y2, y2' = -y1, y1(0)=0, y2(0)=1, (solution y1(t)=sin(t),
c      y2(t)=cos(t)).
c
c    srkode is called first, then sprode with both options of mf,
c    and then smvode with its three mf options.  more accurate results
c    could be obtained if hmin or eps were decreased.
c
c    note the three subroutines which are declared external in this
c    program.
c
c    the subroutine sypone evaluates the right hand side of the odes.
c    in this example, therefore, it assigns the following values
c    to the dy vector-
c    dy(1)=y(2), dy(2)=-y(1), corresponding to the ode above.
c
c    the subroutine spyone, which is only needed when calling smvode
c    with the mf=1 option, takes the partial derivative of right hand
c    side number i with respect to component j of the y vector, and
c    stores it in the entry dermat(i,j), thus doing the following-
c    dermat(1,1)=0,  dermat(1,2)=1,
c    dermat(2,1)=-1, dermat(2,2)=0.
c
c    finally, the subroutine pdummy is a dummy routine which could have
c    been used in place of spyone when calling smvode with mf=0 or
c    mf=2.  since these options don't need a matrix of partial
c    derivatives, passing this dummy name is all that would be
c    required.  in this example, though, we did not use pdummy in the
c    call statement since all three calls (mf=0, 1, and 2) are made
c    with the same statement.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    12 January 2007
c
c  Author:
c
c    John Burkardt
c
      implicit none

      write ( *, '(a)' ) ' '
      call timestamp ( )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'CWG_ODE_PRB'
      write ( *, '(a)' ) '  FORTRAN77 version.'
      write ( *, '(a)' ) '  Test the CWG_ODE library.'

      call test01 ( )
      call test02 ( )
      call test03 ( )
c
c  Terminate.
c
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'CWG_ODE_PRB'
      write ( *, '(a)' ) '  Normal end of execution.'

      write ( *, '(a)' ) ' '
      call timestamp ( )

      stop
      end
      subroutine test01 ( )

c*********************************************************************72
c
cc TEST01 tests SRKODE.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    12 January 2007
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n_max
      parameter ( n_max = 2 )

      real eps
      real error(n_max)
      real hmin
      integer i
      integer j
      integer jstart
      integer kflag
      integer maxstp
      integer n
      external pdummy
      real pi4
      external spyone
      external sypone
      real t
      real tout
      real tstep
      real y(n_max)

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST01'
      write ( *, '(a)' ) '  SRKODE is a Runge-Kutta ODE solver'
c
c  set number of steps to take
c
      maxstp = 12
c
c  set pi/4
c
      pi4 = atan ( 1.0 )
c
c  call srkode
c
      n = 2
      t = 0.0
      tstep = pi4 / 3.0
      y(1) = 0.0
      y(2) = 1.0
      hmin = 0.0001
      eps = 0.0001
      jstart = 0

      write ( *, '(a)' ) ' '
      write ( *, '(a)' )
     &  '              T          Y1            Y2'
      write ( *, '(a)' ) ' '
      i = 0
      write ( *, '(2x,i4,2x,f10.4,2x,5g14.6)' ) i,t,(y(j),j=1,n)

      do i = 1, maxstp

        tout = real ( i ) * tstep

        call srkode ( n, t, tout, y, hmin, eps, error, kflag, jstart, 
     &    sypone )

        write ( *, '(2x,i4,2x,f10.4,2x,5g14.6)' )
     &  i, tout, ( y(j), j = 1, n )

        if ( kflag .ne. 1 ) then
          write ( *, 1050 ) kflag
          go to 20
        end if

      end do

   20 continue

      return
 1050 format(' return flag =',i6)
      end
      subroutine test02 ( )

c*********************************************************************72
c
cc TEST02 tests SPRODE.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    12 January 2007
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n_max
      parameter ( n_max = 2 )

      real eps
      real error(n_max)
      real hmin
      integer i
      integer j
      integer jstart
      integer kflag
      integer maxstp
      integer mf
      integer n
      external pdummy
      real pi4
      external spyone
      external sypone
      real t
      real tout
      real tstep
      real y(n_max)

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST02'
      write ( *, '(a)' ) '  SPRODE is a polyomial/rational function'
      write ( *, '(a)' ) '  approximation ODE solver'
c
c  Set number of steps to take
c
      maxstp = 12
c
c  Set pi/4
c
      pi4 = atan ( 1.0 )
c
c  Call sprode
c
      mf = 0

   30 continue

      n = 2
      t = 0.0
      tstep = pi4 / 3.0
      y(1) = 0.0
      y(2) = 1.0
      hmin = 0.0001
      eps = 1.0e-04
      jstart = 0
      write ( *, '(a)' ) ' '
      write ( *, '(a,i8)' ) '  Using MF = ', mf
      write ( *, '(a)' ) ' '
      write ( *, '(a)' )
     &  '              T          Y1            Y2'
      write ( *, '(a)' ) ' '
      i = 0
      write ( *, '(2x,i4,2x,f10.4,2x,5g14.6)' ) i,t,(y(j),j=1,n)
      do i = 1, maxstp
        tout = t + tstep
        call sprode ( n, t, tout, y, hmin, eps, mf, error, kflag,
     &    jstart, sypone )
        write ( *, '(2x,i4,2x,f10.4,2x,5g14.6)' ) i,t,(y(j),j=1,n)
        if ( kflag .ne. 1 ) then
          write(6,1050)kflag
          go to 50
        end if
      end do

   50 continue
c
c  if mf=0 was just run with sprode, run again with mf=1
c
      if(mf.ne.0) go to 60
      mf=1
      go to 30
   60 continue

      return
 1050 format(' return flag =',i6)
      end
      subroutine test03 ( )

c*********************************************************************72
c
cc TEST03 tests SMVODE.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    12 January 2007
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n_max
      parameter ( n_max = 2 )

      real eps
      real error(n_max)
      real hmin
      integer i
      integer j
      integer jstart
      integer kflag
      integer maxstp
      integer mf
      integer n
      external pdummy
      real pi4
      external spyone
      external sypone
      real t
      real tout
      real tstep
      real y(n_max)

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST03'
      write ( *, '(a)' ) '  SMVODE is a multi-value ODE solver'
c
c  set number of steps to take
c
      maxstp=12
c
c  set pi/4
c
      pi4=atan(1.0)
c
c  call smvode
c
      do mf = 0, 2

        n = 2
        t = 0.0
        tstep = pi4 / 3.0
        y(1) = 0.0
        y(2) = 1.0
        hmin = 0.0001
        eps = 0.00001
        jstart = 0
        write ( *, '(a)' ) ' '
        write ( *, '(a,i8)' ) '  Using MF = ', mf
        write ( *, '(a)' ) ' '
        write ( *, '(a)' )
     &  '              T          Y1            Y2'
        write ( *, '(a)' ) ' '
        i = 0
        write ( *, '(2x,i4,2x,f10.4,2x,5g14.6)' ) i,t,(y(j),j=1,n)

        do i = 1, maxstp

          tout = t + tstep

          call smvode(n,t,tout,y,hmin,eps,mf,error,kflag,
     &      jstart,sypone,spyone)

          write ( *, '(2x,i4,2x,f10.4,2x,5g14.6)' ) i,t,(y(j),j=1,n)

          if(kflag.ne.1)then
            write(6,1050)kflag
           go to 90
          end if

        end do

   90   continue

      end do

      return
 1050 format(' return flag =',i6)
      end
      subroutine sypone(t,y,dy,n)

c*********************************************************************72
c
cc SYPONE
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    12 January 2007
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n

      real dy(n)
      real t
      real y(n)

      dy(1)=y(2)
      dy(2)=-y(1)

      return
      end
      subroutine spyone(t,y,dermat,n)

c*********************************************************************72
c
cc SPYONE
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    12 January 2007
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n

      real dermat(n,n)
      real t
      real y(n)

      dermat(1,1) =  0.0
      dermat(1,2) =  1.0
      dermat(2,1) = -1.0
      dermat(2,2) =  0.0

      return
      end
