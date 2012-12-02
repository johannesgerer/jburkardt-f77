      subroutine fd1d_heat_explicit ( x_num, x, t, dt, cfl, rhs, bc, 
     &  h, h_new )

c*********************************************************************72
c
cc FD1D_HEAT_EXPLICIT: Finite difference solution of 1D heat equation.
c
c  Discussion:
c
c    This program takes one time step to solve the 1D heat equation 
c    with an explicit method.
c
c    This program solves
c
c      dUdT - k * d2UdX2 = F(X,T)
c
c    over the interval [A,B] with boundary conditions
c
c      U(A,T) = UA(T),
c      U(B,T) = UB(T),
c
c    over the time interval [T0,T1] with initial conditions
c
c      U(X,T0) = U0(X)
c
c    The code uses the finite difference method to approximate the
c    second derivative in space, and an explicit forward Euler approximation
c    to the first derivative in time.
c
c    The finite difference form can be written as
c
c      U(X,T+dt) - U(X,T)                  ( U(X-dx,T) - 2 U(X,T) + U(X+dx,T) )
c      ------------------  = F(X,T) + k *  ------------------------------------
c               dt                                   dx * dx
c
c    or, assuming we have solved for all values of U at time T, we have
c
c      U(X,T+dt) = U(X,T) 
c        + cfl * ( U(X-dx,T) - 2 U(X,T) + U(X+dx,T) ) + dt * F(X,T) 
c
c    Here "cfl" is the Courant-Friedrichs-Loewy coefficient:
c
c      cfl = k * dt / dx / dx
c
c    In order for accurate results to be computed by this explicit method,
c    the CFL coefficient must be less than 0.5.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    25 January 2012
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer X_NUM, the number of points to use in the 
c    spatial dimension.
c
c    Input, double precision X(X_NUM), the coordinates of the nodes.
c
c    Input, double precision T, the current time.
c
c    Input, double precision DT, the size of the time step.
c
c    Input, double precision CFL, the Courant-Friedrichs-Loewy coefficient,
c    computed by FD1D_HEAT_EXPLICIT_CFL.
c
c    Input, double precision H(X_NUM), the solution at the current time.
c
c    Input, external RHS, the function which evaluates the right hand side.
c
c    Input, external BC, the function which evaluates the boundary conditions.
c
c    Output, double precision H_NEW(X_NUM), the solution at time T+DT.
c
      implicit none

      integer x_num

      double precision cfl
      double precision f(x_num)
      double precision dt
      double precision h(x_num)
      double precision h_new(x_num)
      integer j
      double precision t
      double precision x(x_num)

      call rhs ( x_num, x, t, f )

      h_new(1) = 0.0D+00

      do j = 2, x_num - 1
        h_new(j) = h(j) + dt * f(j) 
     &    + cfl * (             h(j-1) 
     &              - 2.0D+00 * h(j) 
     &              +           h(j+1) )
      end do

      h_new(x_num) = 0.0D+00

      call bc ( x_num, x, t + dt, h_new )

      return
      end
      subroutine fd1d_heat_explicit_cfl ( k, t_num, t_min, t_max, 
     &  x_num, x_min, x_max, cfl )

c*********************************************************************72
c
cc FD1D_HEAT_EXPLICIT_CFL: compute the Courant-Friedrichs-Loewy coefficient.
c
c  Discussion:
c
c    The equation to be solved has the form:
c
c      dUdT - k * d2UdX2 = F(X,T)
c
c    over the interval [X_MIN,X_MAX] with boundary conditions
c
c      U(X_MIN,T) = U_X_MIN(T),
c      U(X_MIN,T) = U_X_MAX(T),
c
c    over the time interval [T_MIN,T_MAX] with initial conditions
c
c      U(X,T_MIN) = U_T_MIN(X)
c
c    The code uses the finite difference method to approximate the
c    second derivative in space, and an explicit forward Euler approximation
c    to the first derivative in time.
c
c    The finite difference form can be written as
c
c      U(X,T+dt) - U(X,T)                  ( U(X-dx,T) - 2 U(X,T) + U(X+dx,T) )
c      ------------------  = F(X,T) + k *  ------------------------------------
c               dt                                   dx * dx
c
c    or, assuming we have solved for all values of U at time T, we have
c
c      U(X,T+dt) = U(X,T) 
c        + cfl * ( U(X-dx,T) - 2 U(X,T) + U(X+dx,T) ) + dt * F(X,T) 
c
c    Here "cfl" is the Courant-Friedrichs-Loewy coefficient:
c
c      cfl = k * dt / dx / dx
c
c    In order for accurate results to be computed by this explicit method,
c    the CFL coefficient must be less than 0.5c
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    25 January 2012
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
c    Input, double precision K, the heat conductivity coefficient.
c
c    Input, integer T_NUM, the number of time values, including 
c    the initial value.
c
c    Input, double precision T_MIN, T_MAX, the minimum and maximum times.
c
c    Input, integer X_NUM, the number of nodes.
c
c    Input, double precision X_MIN, X_MAX, the minimum and maximum spatial 
c    coordinates.
c
c    Output, double precision CFL, the Courant-Friedrichs-Loewy coefficient.
c
      implicit none

      double precision cfl
      double precision dx
      double precision dt
      double precision k
      double precision t_max
      double precision t_min
      integer t_num
      double precision x_max
      double precision x_min
      integer x_num

      dx = ( x_max - x_min ) / dble ( x_num - 1 )
      dt = ( t_max - t_min ) / dble ( t_num - 1 )
c
c  Check the CFL condition, print out its value, and quit if it is too large.
c
      cfl = k * dt / dx / dx

      write ( *, '(a)' ) ' '
      write ( *, '(a,g14.6)' ) '  CFL stability criterion value = ', cfl

      if ( 0.5D+00 .le. cfl ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'FD1D_HEAT_EXPLICIT_CFL - Fatal error!'
        write ( *, '(a)' ) '  CFL condition failed.'
        write ( *, '(a)' ) '  0.5 <= K * dT / dX / dX = CFL.'
        stop
      end if

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
      subroutine r8vec_write ( output_filename, n, x )

c*********************************************************************72
c
cc R8VEC_WRITE writes an R8VEC file.
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
c    10 July 2011
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, character * ( * ) OUTPUT_FILENAME, the output file name.
c
c    Input, integer N, the number of points.
c
c    Input, double precision X(N), the data.
c
      implicit none

      integer m
      integer n

      integer j
      character * ( * ) output_filename
      integer output_unit
      double precision x(n)
c
c  Open the file.
c
      call get_unit ( output_unit )

      open ( unit = output_unit, file = output_filename,
     &  status = 'replace' )
c
c  Create the format string.
c
      if ( 0 .lt. n ) then
c
c  Write the data.
c
        do j = 1, n
          write ( output_unit, '(2x,g24.16)' ) x(j)
        end do

      end if
c
c  Close the file.
c
      close ( unit = output_unit )

      return
      end
      subroutine timestamp ( )

c*********************************************************************72
c
cc TIMESTAMP prints out the current YMDHMS date as a timestamp.
c
c  Discussion:
c
c    This FORTRAN77 version is made available for cases where the
c    FORTRAN90 version cannot be used.
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
