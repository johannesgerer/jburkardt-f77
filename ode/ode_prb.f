      program main

c*********************************************************************72
c
cc MAIN is the main program for ODE_PRB.
c
c  Discussion:
c
c    ODE_PRB calls a set of problems for ODE.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    31 May 2008
c
c  Author:
c
c    John Burkardt
c
      implicit none

      call timestamp ( )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'ODE_PRB'
      write ( *, '(a)' ) '  FORTRAN77 version'
      write ( *, '(a)' ) '  Test the ODE library.'

      call test01 ( )
      call test02 ( )
c
c  Terminate.
c
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'ODE_PRB'
      write ( *, '(a)' ) '  Normal end of execution.'

      write ( *, '(a)' ) ' '
      call timestamp ( )

      stop
      end
      subroutine test01 ( )

c*********************************************************************72
c
cc TEST01 tests ODE.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    02 February 2012
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer neqn
      parameter ( neqn = 2 )

      double precision abserr
      external f01
      integer i
      integer iflag
      integer iwork(5)
      double precision pi
      parameter ( pi = 3.141592653589793D+00 )
      double precision relerr
      integer step_num
      parameter ( step_num = 12 )
      double precision t
      double precision tout
      double precision work(100+21*neqn)
      double precision y(neqn)

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST01'
      write ( *, '(a)' ) '  ODE solves a system of ordinary'
      write ( *, '(a)' ) '  differential equations.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '      T           Y(1)         Y(2)'
      write ( *, '(a)' ) ' '

      abserr = 0.00001D+00
      relerr = 0.00001D+00

      iflag = 1

      t = 0.0D+00
      y(1) = 1.0D+00
      y(2) = 0.0D+00

      write ( *, '(2x,f8.4,2x,2g14.6)' ) t, y(1:neqn)

      do i = 1, step_num

        tout = dble ( i ) * 2.0D+00 * pi / dble ( step_num )

        call ode ( f01, neqn, y, t, tout, relerr, abserr, iflag,
     &    work, iwork )

        if ( iflag .ne. 2 ) then
          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) 'TEST01 - Fatal errorc'
          write ( *, '(a,i8)' ) '  ODE returned IFLAG = ', iflag
          exit
        end if

        write ( *, '(2x,f8.4,2x,2g14.6)' ) t, y(1:neqn)

      end do

      return
      end
      subroutine test02 ( )

c*********************************************************************72
c
cc TEST02 tests ODE by integrating in the NEGATIVE time direction.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    02 November 2006
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer neqn
      parameter ( neqn = 2 )

      double precision abserr
      external f01
      integer i
      integer iflag
      integer iwork(5)
      double precision pi
      parameter ( pi = 3.141592653589793D+00 )
      double precision relerr
      integer step_num
      parameter ( step_num = 12 )
      double precision t
      double precision tout
      double precision work(100+21*neqn)
      double precision y(neqn)

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST02'
      write ( *, '(a)' ) '  ODE solves a system of ordinary'
      write ( *, '(a)' ) '  differential equations.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  In this example, we integrate in the'
      write ( *, '(a)' ) '  negative time direction.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '      T           Y(1)         Y(2)'
      write ( *, '(a)' ) ' '

      abserr = 0.00001D+00
      relerr = 0.00001D+00

      iflag = 1

      t = 0.0D+00
      y(1) = 1.0D+00
      y(2) = 0.0D+00

      write ( *, '(2x,f8.4,2x,2g14.6)' ) t, y(1:neqn)

      do i = 1, step_num

        tout = - dble ( i ) * 2.0D+00 * pi / dble ( step_num )

        call ode ( f01, neqn, y, t, tout, relerr, abserr, iflag,
     &    work, iwork )

        if ( iflag /= 2 ) then
          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) 'TEST02 - Fatal errorc'
          write ( *, '(a,i8)' ) '  ODE returned IFLAG = ', iflag
          exit
        end if

        write ( *, '(2x,f8.4,2x,2g14.6)' ) t, y(1:neqn)

      end do

      return
      end
      subroutine f01 ( t, y, yp )

c*********************************************************************72
c
cc F01 supplies the right hand side of the ODE for problem 1.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    02 February 2012
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, double precision T, the time.
c
c    Input, double precision Y(), the dependent variable.
c
c    Output, double precision YP(), the value of the derivative.
c
      implicit none

      double precision t
      double precision y(2)
      double precision yp(2)

      yp(1) =   y(2)
      yp(2) = - y(1)

      return
      end
