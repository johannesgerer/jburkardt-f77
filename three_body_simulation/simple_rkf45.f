      program main

c*********************************************************************72
c
cc MAIN is the main program for SIMPLE_RKF45.
c
c  Discussion:
c
c    SIMPLE_RKF45 uses RKF45 as an integrator for the simple version
c    of the three-body problem.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    07 October 2012
c
c  Author:
c
c    John Burkardt
c
      implicit none

      call timestamp ( )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'SIMPLE_RKF45'
      write ( *, '(a)' ) '  FORTRAN77 version'
      write ( *, '(a)' ) '  Simulate the behavior of three bodies which'
      write ( *, '(a)' ) '  are constrained to lie in a plane, moving'
      write ( *, '(a)' ) '  under the influence of gravity.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Use RKF45 for the ODE integrator.'

      call simple_rkf45_run ( )
c
c  Terminate.
c
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'SIMPLE_RKF45'
      write ( *, '(a)' ) '  Normal end of execution.'

      write ( *, '(a)' ) ' '
      call timestamp ( )

      stop
      end
      subroutine get_unit ( unit )

c*********************************************************************72
c
cc GET_UNIT returns a free FORTRAN unit number.
c
c  Discussion:
c
c    A "free" FORTRAN unit number is a value between 1 and 99 which
c    is not currently associated with an I/O device.  A free FORTRAN unit
c    number is needed in order to open a file with the OPEN command.
c
c    If UNIT = 0, then no free FORTRAN unit could be found, although
c    all 99 units were checked (except for units 5, 6 and 9, which
c    are commonly reserved for console I/O).
c
c    Otherwise, UNIT is a value between 1 and 99, representing a
c    free FORTRAN unit.  Note that GET_UNIT assumes that units 5 and 6
c    are special, and will never return those values.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    26 October 2008
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, integer UNIT, the free unit number.
c
      implicit none

      integer i
      integer unit

      unit = 0

      do i = 1, 99

        if ( i .ne. 5 .and. i .ne. 6 .and. i .ne. 9 ) then

          open ( unit = i, err = 10, status = 'scratch' )
          close ( unit = i )

          unit = i

          return
        end if

10      continue

      end do

      return
      end
      subroutine r8mat_write ( output_filename, m, n, table )

c*********************************************************************72
c
cc R8MAT_WRITE writes a R8MAT file.
c
c  Discussion:
c
c    An R8MAT is an array of R8's.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    22 October 2009
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, character * ( * ) OUTPUT_FILENAME, the output file name.
c
c    Input, integer M, the spatial dimension.
c
c    Input, integer N, the number of points.
c
c    Input, double precision TABLE(M,N), the data.
c
      implicit none

      integer m
      integer n

      integer j
      character * ( * ) output_filename
      integer output_unit
      character * ( 30 ) string
      double precision table(m,n)
c
c  Open the file.
c
      call get_unit ( output_unit )

      open ( unit = output_unit, file = output_filename,
     &  status = 'replace' )
c
c  Create the format string.
c
      if ( 0 .lt. m .and. 0 .lt. n ) then

        write ( string, '(a1,i8,a1,i8,a1,i8,a1)' )
     &    '(', m, 'g', 24, '.', 16, ')'
c
c  Write the data.
c
        do j = 1, n
          write ( output_unit, string ) table(1:m,j)
        end do

      end if
c
c  Close the file.
c
      close ( unit = output_unit )

      return
      end
      subroutine simple_rkf45_run ( )

c*********************************************************************72
c
cc SIMPLE_RKF45_RUN runs the simple three body ODE system.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    07 October 2012
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer neqn
      parameter ( neqn = 12 )
      integer step_num
      parameter ( step_num = 630 )

      double precision abserr
      integer flag
      integer i
      double precision m0
      double precision m1
      double precision m2
      double precision relerr
      external simple_f
      integer step
      double precision t
      character * ( 80 ) t_filename
      double precision t_out
      double precision t_start
      double precision t_stop
      double precision ts(0:step_num)
      double precision y(neqn)
      character * ( 80 ) y_filename
      double precision yp(neqn)
      double precision ys(neqn,0:step_num)

      t_filename = 'simple_rkf45_t.txt'
      y_filename = 'simple_rkf45_y.txt'

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'SIMPLE_RKF45_RUN'
      write ( *, '(a)' ) 
     &  '  Simulate the planar three-body problem as an ODE system'
      write ( *, '(a)' ) '  using RKF45 for the ODE integration.'

      m0 = 5.0D+00
      m1 = 3.0D+00
      m2 = 4.0D+00

      abserr = 1.0D-10
      relerr = 1.0D-10

      flag = 1

      t_start = 0.0D+00
      t_stop = 63.0D+00

      t = 0.0D+00
      t_out = 0.0D+00

      y(1)  =  1.0D+00
      y(2)  = -1.0D+00
      y(3)  =  0.0D+00
      y(4)  =  0.0D+00
      y(5)  =  1.0D+00
      y(6)  =  3.0D+00
      y(7)  =  0.0D+00
      y(8)  =  0.0D+00
      y(9)  = -2.0D+00
      y(10) = -1.0D+00
      y(11) =  0.0D+00
      y(12) =  0.0D+00

      call simple_f ( t, y, yp )

      do i = 1, neqn
        ys(i,0) = y(i)
      end do
      ts(0) = t

      do step = 1, step_num

        t = ( dble ( step_num - step + 1 ) * t_start 
     &      + dble (            step - 1 ) * t_stop ) 
     &      / dble ( step_num            )

        t_out = ( dble ( step_num - step ) * t_start 
     &          + dble (            step ) * t_stop ) 
     &          / dble ( step_num        )

        call r8_rkf45 ( simple_f, neqn, y, yp, t, t_out, relerr, 
     &    abserr, flag )

        if ( abs ( flag ) .ne. 2 ) then
          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) 'SIMPLE_RKF45_RUN - Warningc'
          write ( *, '(a,i4)' ) '  Output value of FLAG = ', flag
        end if

        do i = 1, neqn
          ys(i,step) = y(i)
        end do
        ts(step) = t_out

      end do

      call r8mat_write ( t_filename, 1, step_num + 1, ts )
      call r8mat_write ( y_filename, neqn, step_num + 1, ys )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'SIMPLE_RKF45_RUN:'
      write ( *, '(a)' ) 
     &  '  Time data written to "' // trim ( t_filename ) // '".'
      write ( *, '(a)' ) 
     &  '  Solution data written to "' // trim ( y_filename ) // '".'

      return
      end
      subroutine simple_f ( t, y, yp )

c*********************************************************************72
c
cc SIMPLE_F returns the right hand side of the three body ode system.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    03 April 2011
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, double precision T, the value of the independent variable.
c
c    Input, double precision Y(NEQN), the value of the dependent variable.
c
c    Output, double precision YP(NEQN), the value of the derivative
c    dY(1:NEQN)/dT.
c
      implicit none

      integer neqn
      parameter ( neqn = 12 )

      double precision m0
      double precision m1
      double precision m2
      double precision n0
      double precision n1
      double precision n2
      double precision t
      double precision x0
      double precision x1
      double precision x2
      double precision y(neqn)
      double precision y0
      double precision y1
      double precision y2
      double precision yp(neqn)

      m0 = 5.0D+00
      m1 = 3.0D+00
      m2 = 4.0D+00

      x0 = y(1)
      y0 = y(2)

      x1 = y(5)
      y1 = y(6)

      x2 = y(9)
      y2 = y(10)

      n0 = sqrt ( ( ( x2 - x1 )**2 + ( y2 - y1 )**2 )**3 ) 
      n1 = sqrt ( ( ( x0 - x2 )**2 + ( y0 - y2 )**2 )**3 ) 
      n2 = sqrt ( ( ( x1 - x0 )**2 + ( y1 - y0 )**2 )**3 ) 

      yp(1)  =  y(3)
      yp(2)  =  y(4)
      yp(3)  = - m1 * ( x0 - x1 ) / n2 - m2 * ( x0 - x2 ) / n1
      yp(4)  = - m1 * ( y0 - y1 ) / n2 - m2 * ( y0 - y2 ) / n1
      yp(5)  =  y(7)
      yp(6)  =  y(8)
      yp(7)  = - m2 * ( x1 - x0 ) / n0 - m0 * ( x1 - x2 ) / n2
      yp(8)  = - m2 * ( y1 - y0 ) / n0 - m0 * ( y1 - y2 ) / n2
      yp(9)  = y(11)
      yp(10) = y(12)
      yp(11) = - m0 * ( x2 - x0 ) / n1 - m1 * ( x2 - x1 ) / n0
      yp(12) = - m0 * ( y2 - y0 ) / n1 - m1 * ( y2 - y1 ) / n0

      return
      end
