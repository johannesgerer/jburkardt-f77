      subroutine fd1d_wave_alpha ( x_num, x1, x2, t_num, t1, t2, c, 
     &  alpha )

c*********************************************************************72
c
cc FD1D_WAVE_ALPHA computes ALPHA for the 1D wave equation.
c
c  Discussion:
c
c    The explicit timestepping procedure uses the quantity ALPHA, which
c    is determined by this function.
c
c    If the spatial region bounds are X1 <= X <= X2, containing X_NUM equally
c    spaced nodes, including the endpoints, and the time domain similarly
c    extends from T1 <= T <= T2 containing T_NUM equally spaced time values,
c    then
c
c      ALPHA = C * DT / DX
c            = C * ( ( T2 - T1 ) / ( T_NUM - 1 ) )
c                / ( ( X2 - X1 ) / ( X_NUM - 1 ) ).
c
c    For a stable computation, it must be the case that ALPHA < 1.
c
c    If ALPHA is greater than 1, then the middle coefficient 1-C^2 DT^2 / DX^2 
c    is negative, and the sum of the magnitudes of the three coefficients 
c    becomes unbounded.  In such a case, the user must reduce the time step 
c    size appropriately.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    24 January 2012
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    George Lindfield, John Penny,
c    Numerical Methods Using MATLAB,
c    Second Edition,
c    Prentice Hall, 1999,
c    ISBN: 0-13-012641-1,
c    LC: QA297.P45.
c
c  Parameters:
c
c    Input, integer X_NUM, the number of nodes in the X direction.
c
c    Input, double precision X1, X2, the first and last X coordinates.
c
c    Input, integer T_NUM, the number of time steps, including the 
c    initial condition.
c
c    Input, double precision T1, T2, the first and last T coordinates.
c
c    Input, double precision C, a parameter which gives the speed of waves.
c
c    Output, double precision ALPHA, the stability coefficient.
c
      implicit none

      double precision alpha
      double precision c
      double precision t_delta
      integer t_num
      double precision t1
      double precision t2
      double precision x_delta
      integer x_num
      double precision x1
      double precision x2

      t_delta = ( t2 - t1 ) / dble ( t_num - 1 )
      x_delta = ( x2 - x1 ) / dble ( x_num - 1 )
      alpha = c * t_delta / x_delta

      write ( *, '(a)' ) ' '
      write ( *, '(a,g14.6)' ) 
     &  '  Stability condition ALPHA = C * DT / DX = ', alpha

      if ( 1.0D+00 .lt. abs ( alpha ) ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'FD1D_WAVE_ALPHA - Warning!'
        write ( *, '(a)' ) 
     &    '  The stability condition |ALPHA| <= 1 fails.'
        write ( *, '(a)' ) 
     &    '  Computed results are liable to be inaccurate.'
      end if

      return
      end
      subroutine fd1d_wave_start ( x_num, x_vec, t, t_delta, alpha, 
     &  u_x1, u_x2, ut_t1, u1, u2 )

c*********************************************************************72
c
cc FD1D_WAVE_START takes the first step for the wave equation.
c
c  Discussion:
c
c    This program solves the 1D wave equation of the form:
c
c      Utt = c^2 Uxx
c
c    over the spatial interval [X1,X2] and time interval [T1,T2],
c    with initial conditions:
c
c      U(T1,X)  = U_T1(X),
c      Ut(T1,X) = UT_T1(X),
c
c    and boundary conditions of Dirichlet type:
c
c      U(T,X1) = U_X1(T),
c      U(T,X2) = U_X2(T).
c
c    The value C represents the propagation speed of waves.
c
c    The program uses the finite difference method, and marches
c    forward in time, solving for all the values of U at the next
c    time step by using the values known at the previous two time steps.
c
c    Central differences may be used to approximate both the time
c    and space derivatives in the original differential equation.
c
c    Thus, assuming we have available the approximated values of U
c    at the current and previous times, we may write a discretized
c    version of the wave equation as follows:
c
c      Uxx(T,X) = ( U(T,   X+dX) - 2 U(T,X) + U(T,   X-dX) ) / dX^2
c      Utt(T,X) = ( U(T+dt,X   ) - 2 U(T,X) + U(T-dt,X   ) ) / dT^2
c
c    If we multiply the first term by C^2 and solve for the single
c    unknown value U(T+dt,X), we have:
c
c      U(T+dT,X) =        (     C^2 * dT^2 / dX^2 ) * U(T,   X+dX)
c                  +  2 * ( 1 - C^2 * dT^2 / dX^2 ) * U(T,   X   )
c                  +      (     C^2 * dT^2 / dX^2 ) * U(T,   X-dX)
c                  -                                  U(T-dT,X   )
c
c    (Equation to advance from time T to time T+dT, except for FIRST step!)
c
c    However, on the very first step, we only have the values of U
c    for the initial time, but not for the previous time step.
c    In that case, we use the initial condition information for dUdT
c    which can be approximated by a central difference that involves
c    U(T+dT,X) and U(T-dT,X):
c
c      dU/dT(T,X) = ( U(T+dT,X) - U(T-dT,X) ) / ( 2 * dT )
c
c    and so we can estimate U(T-dT,X) as
c
c      U(T-dT,X) = U(T+dT,X) - 2 * dT * dU/dT(T,X)
c
c    If we replace the "missing" value of U(T-dT,X) by the known values
c    on the right hand side, we now have U(T+dT,X) on both sides of the
c    equation, so we have to rearrange to get the formula we use
c    for just the first time step:
c
c      U(T+dT,X) =   1/2 * (     C^2 * dT^2 / dX^2 ) * U(T,   X+dX)
c                  +       ( 1 - C^2 * dT^2 / dX^2 ) * U(T,   X   )
c                  + 1/2 * (     C^2 * dT^2 / dX^2 ) * U(T,   X-dX)
c                  +  dT *                         dU/dT(T,   X   )
c
c    (Equation to advance from time T to time T+dT for FIRST step.)
c
c    It should be clear now that the quantity ALPHA = C * DT / DX will affect
c    the stability of the calculation.  If it is greater than 1, then
c    the middle coefficient 1-C^2 DT^2 / DX^2 is negative, and the
c    sum of the magnitudes of the three coefficients becomes unbounded.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    24 January 2012
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    George Lindfield, John Penny,
c    Numerical Methods Using MATLAB,
c    Second Edition,
c    Prentice Hall, 1999,
c    ISBN: 0-13-012641-1,
c    LC: QA297.P45.
c
c  Parameters:
c
c    Input, integer X_NUM, the number of nodes in the X direction.
c
c    Input, double precision X_VEC(X_NUM), the coordinates of the nodes.
c
c    Input, double precision T, the time after the first step has been taken.
c    In other words, T = T1 + T_DELTA.
c
c    Input, double precision T_DELTA, the time step.
c
c    Input, double precision ALPHA, the stability coefficient, computed 
c    by FD1D_WAVE_ALPHA.
c
c    Input, double precision U_X1(T), U_X2(T), functions for the left and 
c    right boundary conditions.
c
c    Input, double precision UT_T1(X), the function that evaluates dUdT at the 
c    initial time.
c
c    Input, real U1(X_NUM), the initial condition.
c
c    Output, real U2(X_NUM), the solution at the first time step.
c
      implicit none

      integer x_num

      double precision alpha
      integer i
      double precision t
      double precision t_delta
      external u_x1
      external u_x2
      external ut_t1
      double precision u1(x_num)
      double precision u2(x_num)
      double precision ut(x_num)
      double precision x_vec(x_num)

      call ut_t1 ( x_num, x_vec, ut )

      call u_x1 ( t, u2(1) )

      do i = 2, x_num - 1
        u2(i) =               alpha**2   * u1(i+1) / 2.0D+00 
     &          + ( 1.0D+00 - alpha**2 ) * u1(i) 
     &          +             alpha**2   * u1(i-1) / 2.0D+00 
     &          +             t_delta    * ut(i)
      end do

      call u_x2 ( t, u2(x_num) )

      return
      end
      subroutine fd1d_wave_step ( x_num, t, alpha, u_x1, u_x2, u1, 
     &  u2, u3 )

c*********************************************************************72
c
cc FD1D_WAVE_STEP computes a step of the 1D wave equation.
c
c  Discussion:
c
c    This program solves the 1D wave equation of the form:
c
c      Utt = c^2 Uxx
c
c    over the spatial interval [X1,X2] and time interval [T1,T2],
c    with initial conditions:
c
c      U(T1,X)  = U_T1(X),
c      Ut(T1,X) = UT_T1(X),
c
c    and boundary conditions of Dirichlet type:
c
c      U(T,X1) = U_X1(T),
c      U(T,X2) = U_X2(T).
c
c    The value C represents the propagation speed of waves.
c
c    The program uses the finite difference method, and marches
c    forward in time, solving for all the values of U at the next
c    time step by using the values known at the previous two time steps.
c
c    Central differences may be used to approximate both the time
c    and space derivatives in the original differential equation.
c
c    Thus, assuming we have available the approximated values of U
c    at the current and previous times, we may write a discretized
c    version of the wave equation as follows:
c
c      Uxx(T,X) = ( U(T,   X+dX) - 2 U(T,X) + U(T,   X-dX) ) / dX^2
c      Utt(T,X) = ( U(T+dt,X   ) - 2 U(T,X) + U(T-dt,X   ) ) / dT^2
c
c    If we multiply the first term by C^2 and solve for the single
c    unknown value U(T+dt,X), we have:
c
c      U(T+dT,X) =        (     C^2 * dT^2 / dX^2 ) * U(T,   X+dX)
c                  +  2 * ( 1 - C^2 * dT^2 / dX^2 ) * U(T,   X   )
c                  +      (     C^2 * dT^2 / dX^2 ) * U(T,   X-dX)
c                  -                                  U(T-dT,X   )
c
c    (Equation to advance from time T to time T+dT, except for FIRST step!)
c
c    However, on the very first step, we only have the values of U
c    for the initial time, but not for the previous time step.
c    In that case, we use the initial condition information for dUdT
c    which can be approximated by a central difference that involves
c    U(T+dT,X) and U(T-dT,X):
c
c      dU/dT(T,X) = ( U(T+dT,X) - U(T-dT,X) ) / ( 2 * dT )
c
c    and so we can estimate U(T-dT,X) as
c
c      U(T-dT,X) = U(T+dT,X) - 2 * dT * dU/dT(T,X)
c
c    If we replace the "missing" value of U(T-dT,X) by the known values
c    on the right hand side, we now have U(T+dT,X) on both sides of the
c    equation, so we have to rearrange to get the formula we use
c    for just the first time step:
c
c      U(T+dT,X) =   1/2 * (     C^2 * dT^2 / dX^2 ) * U(T,   X+dX)
c                  +       ( 1 - C^2 * dT^2 / dX^2 ) * U(T,   X   )
c                  + 1/2 * (     C^2 * dT^2 / dX^2 ) * U(T,   X-dX)
c                  +  dT *                         dU/dT(T,   X   )
c
c    (Equation to advance from time T to time T+dT for FIRST step.)
c
c    It should be clear now that the quantity ALPHA = C * DT / DX will affect
c    the stability of the calculation.  If it is greater than 1, then
c    the middle coefficient 1-C^2 DT^2 / DX^2 is negative, and the
c    sum of the magnitudes of the three coefficients becomes unbounded.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    24 January 2012
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    George Lindfield, John Penny,
c    Numerical Methods Using MATLAB,
c    Second Edition,
c    Prentice Hall, 1999,
c    ISBN: 0-13-012641-1,
c    LC: QA297.P45.
c
c  Parameters:
c
c    Input, integer X_NUM, the number of nodes in the X direction.
c
c    Input, double precision T, the new time, that is, the current 
c    time + T_DELTA.
c
c    Input, double precision ALPHA, the stability coefficient, computed 
c    by FD1D_WAVE_ALPHA.
c
c    Input, double precision U_X1(T), U_X2(T), functions for the left and 
c    right boundary conditions.
c
c    Input, double precision U1(X_NUM), the solution at the old time.
c
c    Input, double precision U2(X_NUM), the solution at the current time.
c
c    Output, double precision U3(X_NUM), the solution at the new time.
c
      implicit none

      integer x_num

      double precision alpha
      integer i
      double precision t
      external u_x1
      external u_x2
      double precision u1(x_num)
      double precision u2(x_num)
      double precision u3(x_num)

      call u_x1 ( t, u3(1) )

      do i = 2, x_num - 1
        u3(i) =                         alpha**2   * u2(i+1) 
     &          + 2.0D+00 * ( 1.0D+00 - alpha**2 ) * u2(i)
     &          +                       alpha**2   * u2(i-1) 
     &          -                                    u1(i)
      end do

      call u_x2 ( t, u3(x_num) )


      return
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
      subroutine piecewise_linear ( nd, xd, yd, nv, xv, yv )

c*********************************************************************72
c
cc PIECEWISE_LINEAR evaluates a piecewise linear spline.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    24 January 2012
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer ND, the number of data points.
c
c    Input, double precision XD(ND), YD(ND), the data values.
c
c    Input, integer NV, the number of evaluation points.
c
c    Input, double precision XV(NV), the evaluation arguments.
c
c    Output, double precision YV(NV), the values.
c
      implicit none

      integer nd
      integer nv

      integer id
      integer iv
      double precision xd(nd)
      double precision xv(nv)
      double precision yd(nd)
      double precision yv(nv)

      do iv = 1, nv

        if ( xv(iv) .lt. xd(1) ) then
          yv(iv) = yd(1)
        else if ( xd(nd) .lt. xv(iv) ) then
          yv(iv) = yd(nd)
        else 

          do id = 2, nd
            if ( xv(iv) .lt. xd(id) ) then
              yv(iv) = ( ( xd(id) - xv(iv)            ) * yd(id-1) 
     &                 + (          xv(iv) - xd(id-1) ) * yd(id) ) 
     &                 / ( xd(id)          - xd(id-1) )
              go to 10
            end if
          end do

10        continue

        end if

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
      subroutine r8vec_linspace ( n, a_first, a_last, a )

c*********************************************************************72
c
cc R8VEC_LINSPACE creates a vector of linearly spaced values.
c
c  Discussion:
c
c    An R8VEC is a vector of R8's.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    14 March 2011
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the number of entries in the vector.
c
c    Input, double precision A_FIRST, A_LAST, the first and last entries.
c
c    Output, double precision A(N), a vector of linearly spaced data.
c
      implicit none

      integer n

      double precision a(n)
      double precision a_first
      double precision a_last
      integer i

      if ( n .eq. 1 ) then

        a(1) = ( a_first + a_last ) / 2.0D+00

      else

        do i = 1, n
          a(i) = ( dble ( n - i     ) * a_first 
     &           + dble (     i - 1 ) * a_last )
     &           / dble ( n     - 1 )
        end do

      end if

      return
      end
      subroutine timestamp ( )

c*********************************************************************72
c
cc TIMESTAMP prints out the current YMDHMS date as a timestamp.
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
c  Parameters:
c
c    None
c
      implicit none

      character * ( 8 ) ampm
      integer d
      character * ( 8 ) date
      integer h
      integer m
      integer mm
      character * ( 9 ) month(12)
      integer n
      integer s
      character * ( 10 ) time
      integer y

      save month

      data month /
     &  'January  ', 'February ', 'March    ', 'April    ',
     &  'May      ', 'June     ', 'July     ', 'August   ',
     &  'September', 'October  ', 'November ', 'December ' /

      call date_and_time ( date, time )

      read ( date, '(i4,i2,i2)' ) y, m, d
      read ( time, '(i2,i2,i2,1x,i3)' ) h, n, s, mm

      if ( h .lt. 12 ) then
        ampm = 'AM'
      else if ( h .eq. 12 ) then
        if ( n .eq. 0 .and. s .eq. 0 ) then
          ampm = 'Noon'
        else
          ampm = 'PM'
        end if
      else
        h = h - 12
        if ( h .lt. 12 ) then
          ampm = 'PM'
        else if ( h .eq. 12 ) then
          if ( n .eq. 0 .and. s .eq. 0 ) then
            ampm = 'Midnight'
          else
            ampm = 'AM'
          end if
        end if
      end if

      write ( *,
     &  '(i2,1x,a,1x,i4,2x,i2,a1,i2.2,a1,i2.2,a1,i3.3,1x,a)' )
     &  d, month(m), y, h, ':', n, ':', s, '.', mm, ampm

      return
      end
