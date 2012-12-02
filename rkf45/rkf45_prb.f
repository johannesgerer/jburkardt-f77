      program main

c*********************************************************************72
c
cc MAIN is the main program for RKF45_PRB.
c
c  Discussion:
c
c    RKF45_PRB tests the RKF45 ODE integrator.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    17 June 2006
c
c  Author:
c
c    John Burkardt
c
      implicit none

      call timestamp ( )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'RKF45_PRB'
      write ( *, '(a)' ) '  FORTRAN77 version'
      write ( *, '(a)' ) '  Test the RKF45 library.'

      call test01 ( )
      call test02 ( )
      call test03 ( )
      call test04 ( )
      call test05 ( )
      call test06 ( )
c
c  Terminate.
c
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'RKF45_PRB'
      write ( *, '(a)' ) '  Normal end of execution.'

      write ( *, '(a)' ) ' '
      call timestamp ( )

      stop
      end
      subroutine test01 ( )

c*********************************************************************72
c
cc TEST01 solves a scalar ODE.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    17 June 2006
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer neqn
      parameter ( neqn = 1 )

      real abserr
      integer flag
      integer i_step
      integer iwork(5)
      integer n_step
      real r4_epsilon
      external r4_f1
      real r4_y1x
      real relerr
      real t
      real t_out
      real t_start
      real t_stop
      real work(3+6*neqn)
      real y(neqn)

      write ( *, '(a)') ' '
      write ( *, '(a)' ) 'TEST01'
      write ( *, '(a)' ) '  Solve a scalar equation using RKF_S:'
      write ( *, '(a)') ' '
      write ( *, '(a)' ) '  Y'' = 0.25 * Y * ( 1 - Y / 20 )'
      write ( *, '(a)') ' '

      abserr = sqrt ( r4_epsilon ( ) )
      relerr = sqrt ( r4_epsilon ( ) )

      flag = 1

      t_start = 0.0E+00
      t_stop = 20.0E+00

      n_step = 5

      t_out = 0.0E+00
      t = t_out
      y(1) = 1.0E+00

      write ( *, '(a)' ) ' '
      write ( *, '(a,a)' ) '  FLAG     T             Y            ',
     &  'Y_Exact         Error'
      write ( *, '(a)' ) ' '
      write ( *, '(i4,2x,5g14.6)' ) flag, t, y(1), r4_y1x ( t ), 
     &  y(1) - r4_y1x ( t )

      do i_step = 1, n_step

        t = ( real ( n_step - i_step + 1 ) * t_start  
     &      + real (          i_step - 1 ) * t_stop ) 
     &      / real ( n_step              )

        t_out = ( real ( n_step - i_step ) * t_start  
     &          + real (          i_step ) * t_stop ) 
     &          / real ( n_step          )

        call r4_rkf45 ( r4_f1, neqn, y, t, t_out, relerr, abserr, 
     &    flag, work, iwork )

        write ( *, '(i4,2x,5g14.6)' ) flag, t, y(1), 
     &    r4_y1x ( t ), y(1) - r4_y1x ( t )

      end do

      return
      end
      subroutine test02 ( )

c*********************************************************************72
c
cc TEST02 solves a vector ODE.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    17 June 2006
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer neqn
      parameter ( neqn = 2 )

      real abserr
      integer flag
      integer i_step
      integer iwork(5)
      integer n_step
      real r4_epsilon
      external r4_f2
      real relerr
      real t
      real t_out
      real t_start
      real t_stop
      real work(3+6*neqn)
      real y(neqn)

      write ( *, '(a)') ' '
      write ( *, '(a)' ) 'TEST02'
      write ( *, '(a)' ) '  Solve a vector equation using RKF_S:'
      write ( *, '(a)') ' '
      write ( *, '(a)' ) '  Y''(1) =  Y(2)'
      write ( *, '(a)' ) '  Y''(2) = -Y(1)'
      write ( *, '(a)') ' '
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  This system is equivalent to the following'
      write ( *, '(a)' ) '  second order system:'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Z" = - Z.'

      abserr = sqrt ( r4_epsilon ( ) )
      relerr = sqrt ( r4_epsilon ( ) )

      flag = 1

      t_start = 0.0E+00
      t_stop = 2.0E+00 * 3.14159265E+00

      n_step = 12

      t = 0.0E+00
      t_out = 0.0E+00

      y(1) = 1.0E+00
      y(2) = 0.0E+00

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  FLAG       T          Y(1)          Y(2)'
      write ( *, '(a)' ) ' '
      write ( *, '(i4,2x,4g14.6)' ) flag, t, y(1), y(2)

      do i_step = 1, n_step

        t = ( real ( n_step - i_step + 1 ) * t_start 
     &     + real (           i_step - 1 ) * t_stop ) 
     &      / real ( n_step              )

        t_out = ( real ( n_step - i_step ) * t_start 
     &          + real (          i_step ) * t_stop ) 
     &          / real ( n_step          )

        call r4_rkf45 ( r4_f2, neqn, y, t, t_out, relerr, abserr, 
     &    flag, work, iwork )

        write ( *, '(i4,2x,4g14.6)' ) flag, t, y(1), y(2)

      end do

      return
      end
      subroutine test03 ( )

c*********************************************************************72
c
cc TEST03 solves a scalar ODE and uses one-step integration.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    17 June 2006
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer neqn
      parameter ( neqn = 1 )

      real abserr
      integer flag
      integer i_step
      integer iwork(5)
      integer n_step
      real r4_epsilon
      external r4_f1
      real r4_y1x
      real relerr
      real t
      real t_out
      real t_start
      real t_stop
      real work(3+6*neqn)
      real y(neqn)

      write ( *, '(a)') ' '
      write ( *, '(a)' ) 'TEST03'
      write ( *, '(a)' ) '  Solve a scalar equation using RKF_S:'
      write ( *, '(a)') ' '
      write ( *, '(a)' ) '  Y'' = 0.25 * Y * ( 1 - Y / 20 )'
      write ( *, '(a)') ' '
      write ( *, '(a)' ) '  Use the special SINGLE_STEP mode'
      write ( *, '(a)' ) '  which returns after every step.'

      abserr = sqrt ( r4_epsilon ( ) )
      relerr = sqrt ( r4_epsilon ( ) )

      flag = -1

      t_start = 0.0E+00
      t_stop = 20.0E+00

      n_step = 5

      t = 0.0E+00
      t_out = 0.0E+00
      y(1) = 1.0E+00

      write ( *, '(a)' ) ' '
      write ( *, '(a,a)' ) '  FLAG     T             Y           ',
     &  'Y_Exact        Error'
      write ( *, '(a)' ) ' '
      write ( *, '(i4,2x,5g14.6)' ) flag, t, y(1), r4_y1x ( t ), 
     &  y(1) - r4_y1x ( t )

      do i_step = 1, n_step

        t = ( real ( n_step - i_step + 1 ) * t_start 
     &      + real (          i_step - 1 ) * t_stop ) 
     &      / real ( n_step              )

        t_out = ( real ( n_step - i_step ) * t_start  
     &          + real (          i_step ) * t_stop ) 
     &          / real ( n_step          )
c
c  As long as FLAG is negative, we are heading towards T_OUT, but
c  have not reached itc
c
10      continue

        if ( flag < 0 ) then

          call r4_rkf45 ( r4_f1, neqn, y, t, t_out, relerr, abserr, 
     &    flag, work, iwork )

          write ( *, '(i4,2x,5g14.6)' ) flag, t, y(1), 
     &      r4_y1x ( t ), y(1) - r4_y1x ( t )

          go to 10

        end if
c
c  FLAG is returned as +2 when we reach T_OUT.  Reset it to -2
c  to continue to the next T_OUT in one step mode.
c
        flag = -2

      end do

      return
      end
      subroutine test04 ( )

c*********************************************************************72
c
cc TEST04 solves a scalar ODE.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    17 June 2006
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer neqn
      parameter ( neqn = 1 )

      double precision abserr
      integer flag
      integer i_step
      integer iwork(5)
      integer n_step
      double precision r8_epsilon
      external r8_f1
      double precision r8_y1x
      double precision relerr
      double precision t
      double precision t_out
      double precision t_start
      double precision t_stop
      double precision work(3+6*neqn)
      double precision y(neqn)

      write ( *, '(a)') ' '
      write ( *, '(a)' ) 'TEST04'
      write ( *, '(a)' ) '  Solve a scalar equation using RKF_D:'
      write ( *, '(a)') ' '
      write ( *, '(a)' ) '  Y'' = 0.25 * Y * ( 1 - Y / 20 )'
      write ( *, '(a)') ' '

      abserr = sqrt ( r8_epsilon ( ) )
      relerr = sqrt ( r8_epsilon ( ) )

      flag = 1

      t_start = 0.0D+00
      t_stop = 20.0D+00

      n_step = 5

      t_out = 0.0D+00
      t = t_out
      y(1) = 1.0D+00

      write ( *, '(a)' ) ' '
      write ( *, '(a,a)' ) '  FLAG     T             Y            ',
     &  'Y_Exact         Error'
      write ( *, '(a)' ) ' '
      write ( *, '(i4,2x,5g14.6)' ) flag, t, y(1), r8_y1x ( t ), 
     &  y(1) - r8_y1x ( t )

      do i_step = 1, n_step

        t = ( dble ( n_step - i_step + 1 ) * t_start 
     &      + dble (          i_step - 1 ) * t_stop ) 
     &      / dble ( n_step              )

        t_out = ( real ( n_step - i_step ) * t_start 
     &          + real (          i_step ) * t_stop ) 
     &          / real ( n_step          )

        call r8_rkf45 ( r8_f1, neqn, y, t, t_out, relerr, abserr, 
     &    flag, work, iwork )

        write ( *, '(i4,2x,5g14.6)' ) flag, t, y(1),  
     &    r8_y1x ( t ), y(1) - r8_y1x ( t )

      end do

      return
      end
      subroutine test05 ( )

c*********************************************************************72
c
cc TEST05 solves a vector ODE.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    17 June 2006
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer neqn
      parameter ( neqn = 2 )

      double precision abserr
      integer flag
      integer i_step
      integer iwork(5)
      integer n_step
      double precision r8_epsilon
      external r8_f2
      double precision relerr
      double precision t
      double precision t_out
      double precision t_start
      double precision t_stop
      double precision work(3+6*neqn)
      double precision y(neqn)

      write ( *, '(a)') ' '
      write ( *, '(a)' ) 'TEST05'
      write ( *, '(a)' ) '  Solve a vector equation using RKF_D:'
      write ( *, '(a)') ' '
      write ( *, '(a)' ) '  Y''(1) =  Y(2)'
      write ( *, '(a)' ) '  Y''(2) = -Y(1)'
      write ( *, '(a)') ' '

      abserr = sqrt ( r8_epsilon ( ) )
      relerr = sqrt ( r8_epsilon ( ) )

      flag = 1

      t_start = 0.0D+00
      t_stop = 2.0D+00 * 3.14159265D+00

      n_step = 12

      t = 0.0D+00
      t_out = 0.0D+00

      y(1) = 1.0D+00
      y(2) = 0.0D+00

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  FLAG       T          Y(1)          Y(2)'
      write ( *, '(a)' ) ' '
      write ( *, '(i4,2x,4g14.6)' ) flag, t, y(1), y(2)

      do i_step = 1, n_step

        t = ( dble ( n_step - i_step + 1 ) * t_start 
     &      + dble (          i_step - 1 ) * t_stop ) 
     &      / dble ( n_step              )

        t_out = ( real ( n_step - i_step ) * t_start 
     &          + real (          i_step ) * t_stop ) 
     &          / real ( n_step          )

        call r8_rkf45 ( r8_f2, neqn, y, t, t_out, relerr, abserr, 
     &    flag, work, iwork )

        write ( *, '(i4,2x,4g14.6)' ) flag, t, y(1), y(2)

      end do

      return
      end
      subroutine test06 ( )

c*********************************************************************72
c
cc TEST06 solves a scalar ODE and uses one-step integration.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    17 June 2006
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer neqn
      parameter ( neqn = 1 )

      double precision abserr
      integer flag
      integer i_step
      integer iwork(5)
      integer n_step
      double precision r8_epsilon
      external r8_f1
      double precision r8_y1x
      double precision relerr
      double precision t
      double precision t_out
      double precision t_start
      double precision t_stop
      double precision work(3+6*neqn)
      double precision y(neqn)

      write ( *, '(a)') ' '
      write ( *, '(a)' ) 'TEST06'
      write ( *, '(a)' ) '  Solve a scalar equation using RKF_D:'
      write ( *, '(a)') ' '
      write ( *, '(a)' ) '  Y'' = 0.25 * Y * ( 1 - Y / 20 )'
      write ( *, '(a)') ' '
      write ( *, '(a)' ) '  Use the special SINGLE_STEP mode'
      write ( *, '(a)' ) '  which returns after every step.'

      abserr = sqrt ( r8_epsilon ( ) )
      relerr = sqrt ( r8_epsilon ( ) )

      flag = -1

      t_start = 0.0D+00
      t_stop = 20.0D+00

      n_step = 5

      t = 0.0D+00
      t_out = 0.0D+00
      y(1) = 1.0D+00

      write ( *, '(a)' ) ' '
      write ( *, '(a,a)' ) '  FLAG     T             Y           ',
     &  'Y_Exact        Error'
      write ( *, '(a)' ) ' '
      write ( *, '(i4,2x,5g14.6)' ) flag, t, y(1), r8_y1x ( t ), 
     &  y(1) - r8_y1x ( t )

      do i_step = 1, n_step

        t = ( dble ( n_step - i_step + 1 ) * t_start 
     &      + dble (          i_step - 1 ) * t_stop ) 
     &      / dble ( n_step              )

        t_out = ( real ( n_step - i_step ) * t_start 
     &          + real (          i_step ) * t_stop ) 
     &          / real ( n_step          )
c
c  As long as FLAG is negative, we are heading towards T_OUT, but
c  have not reached itc
c
10      continue

        if ( flag < 0 ) then

          call r8_rkf45 ( r8_f1, neqn, y, t, t_out, relerr, abserr, 
     &    flag, work, iwork )

          write ( *, '(i4,2x,4g14.6)' ) flag, t, y(1), 
     &      r8_y1x ( t ), y(1) - r8_y1x ( t )

          go to 10

        end if
c
c  FLAG is returned as +2 when we reach T_OUT.  Reset it to -2
c  to continue to the next T_OUT in one step mode.
c
        flag = -2

      end do

      return
      end
      subroutine r4_f1 ( t, y, yp )

c*********************************************************************72
c
cc R4_F1 evaluates the derivative for the ODE.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    26 March 2004
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, real T, the value of the independent variable.
c
c    Input, real Y(NEQN), the value of the dependent variable.
c
c    Output, real YP(NEQN), the value of the derivative
c    dY(1:NEQN)/dT.
c
      implicit none

      real t
      real y(1)
      real yp(1)

      yp(1) = 0.25E+00 * y(1) * ( 1.0E+00 - y(1) / 20.0E+00 )

      return
      end
      function r4_y1x ( t )

c*********************************************************************72
c
cc R4_Y1X evaluates the exact solution of the ODE.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    26 March 2004
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, real T, the value of the independent variable.
c
c    Output, real Y1X_S, the exact solution.
c
      implicit none

      real t
      real r4_y1x

      r4_y1x = 20.0E+00 / ( 1.0E+00 + 19.0E+00 
     &  * exp ( - 0.25E+00 * t ) )

      return
      end
      subroutine r4_f2 ( t, y, yp )

c*********************************************************************72
c
cc R4_F2 evaluates the derivative for the ODE.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    26 March 2004
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, real T, the value of the independent variable.
c
c    Input, real Y(NEQN), the value of the dependent variable.
c
c    Output, real YP(NEQN), the value of the derivative
c    dY(1:NEQN)/dT.
c
      implicit none

      real t
      real y(2)
      real yp(2)

      yp(1) =  y(2)
      yp(2) = -y(1)

      return
      end
      subroutine r8_f1 ( t, y, yp )

c*********************************************************************72
c
cc R8_F1 evaluates the derivative for the ODE.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    26 March 2004
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

      double precision t
      double precision y(1)
      double precision yp(1)

      yp(1) = 0.25D+00 * y(1) * ( 1.0D+00 - y(1) / 20.0D+00 )

      return
      end
      function r8_y1x ( t )

c*********************************************************************72
c
cc R8_Y1X evaluates the exact solution of the ODE.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    26 March 2004
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, double precision T, the value of the independent variable.
c
c    Output, double precision Y1X_D, the exact solution.
c
      implicit none

      double precision t
      double precision r8_y1x

      r8_y1x = 20.0D+00 / ( 1.0D+00 + 19.0D+00 
     &  * exp ( - 0.25D+00 * t ) )

      return
      end
      subroutine r8_f2 ( t, y, yp )

c*********************************************************************72
c
cc R8_F2 evaluates the derivative for the ODE.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    26 March 2004
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

      double precision t
      double precision y(2)
      double precision yp(2)

      yp(1) =  y(2)
      yp(2) = -y(1)

      return
      end
