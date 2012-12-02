      program main

c*********************************************************************72
c
cc MAIN is the main program for FD1D_HEAT_IMPLICIT.
c
c  Discussion:
c
c    FD1D_HEAT_IMPLICIT solves the 1D heat equation with an implicit method.
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
c    second derivative in space, and an implicit backward Euler approximation
c    to the first derivative in time.
c
c    The finite difference form can be written as
c
c      U(X,T+dt) - U(X,T)                  ( U(X-dx,T+dt) - 2 U(X,T+dt) + U(X+dx,T+dt) )
c      ------------------ = F(X,T+dt) + k *  --------------------------------------
c               dt                                   dx * dx
c
c    so that we have the following linear system for the values of U at time T+dt:
c
c            -     k * dt / dx / dx   * U(X-dt,T+dt)
c      + ( 1 + 2 * k * dt / dx / dx ) * U(X,   T+dt)
c            -     k * dt / dx / dx   * U(X+dt,T+dt)
c      =               dt             * F(X,   T+dt)
c      +                                U(X,   T)
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    30 May 2009
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer t_num
      integer x_num

      parameter ( t_num = 51 )
      parameter ( x_num = 11 )

      double precision a(3,x_num)
      double precision b(x_num)
      double precision fvec(x_num)
      integer i
      integer info
      integer j
      integer job
      double precision k
      double precision t(t_num)
      double precision t_delt
      character * ( 80 ) t_file
      double precision t_max
      double precision t_min
      integer t_unit
      double precision u(x_num,t_num)
      character * ( 80 ) u_file
      integer u_unit
      double precision w
      double precision x(x_num)
      double precision x_delt
      character * ( 80 ) x_file
      double precision x_max
      double precision x_min
      integer x_unit

      call timestamp ( );
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'FD1D_HEAT_IMPLICIT'
      write ( *, '(a)' ) '  FORTRAN77 version'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Finite difference solution of'
      write ( *, '(a)' ) '  the time dependent 1D heat equation'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '    Ut - k * Uxx = F(x,t)'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 
     &  '  for space interval A <= X <= B with boundary conditions'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '    U(A,t) = UA(t)'
      write ( *, '(a)' ) '    U(B,t) = UB(t)'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' )
     &  '  and time interval T0 <= T <= T1 with initial condition'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '    U(X,T0) = U0(X).'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 
     &  '  A second order difference approximation is used for Uxx.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 
     &  '  A first order backward Euler difference approximation'
      write ( *, '(a)' ) '  is used for Ut.'

      k = 5.0D-07
c
c  Set X values.
c
      x_min = 0.0D+00
      x_max = 0.3D+00
      x_delt = ( x_max - x_min ) / dble ( x_num - 1 )

      do i = 1, x_num
        x(i) = ( dble ( x_num - i     ) * x_min   
     &         + dble (         i - 1 ) * x_max ) 
     &         / dble ( x_num     - 1 )
      end do
c 
c  Set T values.
c
      t_min = 0.0D+00
      t_max = 22000.0D+00
      t_delt = ( t_max - t_min ) / dble ( t_num - 1 )

      do j = 1, t_num

        t(j) = ( dble ( t_num - j     ) * t_min   
     &         + dble (         j - 1 ) * t_max ) 
     &         / dble ( t_num     - 1 )
      end do
c
c  Set the initial data, for time T_MIN.
c
      call u0 ( x_min, x_max, t_min, x_num, x, u(1:x_num,1) )
c
c  The matrix A does not change with time.  We can set it once,
c  factor it once, and solve repeatedly.
c
      w = k * t_delt / x_delt / x_delt

      a(1,1) = 0.0D+00

      a(2,1) = 1.0D+00
      a(1,2) = 0.0D+00

      do i = 2, x_num - 1
        a(3,i-1) =                   - w
        a(2,i  ) = 1.0D+00 + 2.0D+00 * w
        a(1,i+1) =                   - w
      end do

      a(3,x_num-1) = 0.0D+00
      a(2,x_num) = 1.0D+00

      a(3,x_num) = 0.0D+00
!
!  Factor the matrix.
!
      call r83_np_fa ( x_num, a, info )

      do j = 2, t_num
!
!  Set the right hand side B.
!
        call ua ( x_min, x_max, t_min, t(j), b(1) )

        call f ( x_min, x_max, t_min, t(j), x_num, x, fvec )

        b(2:x_num-1) = u(2:x_num-1,j-1) + t_delt * fvec(2:x_num-1)

        call ub ( x_min, x_max, t_min, t(j), b(x_num) )

        job = 0
        call r83_np_sl ( x_num, a, b, job )

        u(1:x_num,j) = b(1:x_num)

      end do

      x_file = 'x.txt'
      call get_unit ( x_unit )
      open ( unit = x_unit, file = x_file, status = 'replace' )

      do i = 1, x_num
        write ( x_unit, '(g14.6)' ) x(i)
      end do

      close ( unit = x_unit )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 
     &  '  X data written to "' // trim ( x_file ) // '".'

      t_file = 't.txt'
      call get_unit ( t_unit )
      open ( unit = t_unit, file = t_file, status = 'replace' )

      do j = 1, t_num
        write ( t_unit, '(g14.6)' ) t(j)
      end do

      close ( unit = t_unit )

      write ( *, '(a)' ) 
     &  '  T data written to "' // trim ( t_file ) // '".'

      u_file = 'u.txt'
      call get_unit ( u_unit )
      open ( unit = u_unit, file = u_file, status = 'replace' )

      do j = 1, t_num
        write ( u_unit, '(11g14.6)' ) u(1:x_num,j)
      end do

      close ( unit = u_unit )

      write ( *, '(a)' ) 
     &  '  U data written to "' // trim ( u_file ) // '".'

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'FD1D_HEAT_IMPLICIT'
      write ( *, '(a)' ) '  Normal end of execution.'
      write ( *, '(a)' ) 
      call timestamp ( );

      return
      end
      subroutine f ( a, b, t0, t, n, x, value )

c*********************************************************************72
c
cc F returns the right hand side of the heat equation.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    30 May 2009
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, double precision A, B, the left and right endpoints.
c
c    Input, double precision T0, the initial time.
c
c    Input, double precision T, the current time.
c
c    Input, integer N, the number of points.
c
c    Input, double precision X(N), the current spatial positions.
c
c    Output, double precision VALUE(N), the prescribed value of U(X(:),T0).
c
      implicit none

      integer n

      double precision a
      double precision b
      integer i
      double precision t
      double precision t0
      double precision value(n)
      double precision x(n)

      do i = 1, n
        value(i) = 0.0D+00
      end do

      return
      end
      subroutine get_unit ( unit )

c*********************************************************************72
c
cc GET_UNIT returns a free FORTRAN unit number.
c
c  Discussion:
c
c    A "free" FORTRAN unit number is an integer between 1 and 99 which
c    is not currently associated with an I/O device.  A free FORTRAN unit
c    number is needed in order to open a file with the OPEN command.
c
c    If UNIT = 0, then no free FORTRAN unit could be found, although
c    all 99 units were checked (except for units 5, 6 and 9, which
c    are commonly reserved for console I/O).
c
c    Otherwise, UNIT is an integer between 1 and 99, representing a
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
      subroutine r83_np_fa ( n, a, info )

c*********************************************************************72
c
cc R83_NP_FA factors an R83 matrix without pivoting.
c
c  Discussion:
c
c    The R83 storage format is used for a tridiagonal matrix.
c    The superdiagonal is stored in entries (1,2:N), the diagonal in
c    entries (2,1:N), and the subdiagonal in (3,1:N-1).  Thus, the
c    original matrix is "collapsed" vertically into the array.
c
c    Because this routine does not use pivoting, it can fail even when
c    the matrix is not singular, and it is liable to make larger
c    errors.
c
c    R83_NP_FA and R83_NP_SL may be preferable to the corresponding
c    LINPACK routine SGTSL for tridiagonal systems, which factors and solves
c    in one step, and does not save the factorization.
c
c  Example:
c
c    Here is how an R83 matrix of order 5 would be stored:
c
c       *  A12 A23 A34 A45
c      A11 A22 A33 A44 A55
c      A21 A32 A43 A54  *
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    24 May 2009
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the order of the matrix.
c    N must be at least 2.
c
c    Input/output, double precision A(3,N).
c    On input, the tridiagonal matrix.  On output, factorization information.
c
c    Output, integer INFO, singularity flag.
c    0, no singularity detected.
c    nonzero, the factorization failed on the INFO-th step.
c
      implicit none

      integer n

      double precision a(3,n)
      integer i
      integer info

      info = 0

      do i = 1, n-1

        if ( a(2,i) .eq. 0.0D+00 ) then
          info = i
          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) 'R83_NP_FA - Fatal error!'
          write ( *, '(a,i8)' ) '  Zero pivot on step ', info
          return
        end if
c
c  Store the multiplier in L.
c
        a(3,i) = a(3,i) / a(2,i)
c
c  Modify the diagonal entry in the next column.
c
        a(2,i+1) = a(2,i+1) - a(3,i) * a(1,i+1)

      end do

      if ( a(2,n) .eq. 0.0D+00 ) then
        info = n
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'R83_NP_FA - Fatal error!'
        write ( *, '(a,i8)' ) '  Zero pivot on step ', info
        return
      end if

      return
      end
      subroutine r83_np_sl ( n, a_lu, b, job )

c*********************************************************************72
c
cc R83_NP_SL solves an R83 system factored by R83_NP_FA.
c
c  Discussion:
c
c    The R83 storage format is used for a tridiagonal matrix.
c    The superdiagonal is stored in entries (1,2:N), the diagonal in
c    entries (2,1:N), and the subdiagonal in (3,1:N-1).  Thus, the
c    original matrix is "collapsed" vertically into the array.
c
c  Example:
c
c    Here is how an R83 matrix of order 5 would be stored:
c
c       *  A12 A23 A34 A45
c      A11 A22 A33 A44 A55
c      A21 A32 A43 A54  *
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    24 May 2009
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the order of the matrix.
c    N must be at least 2.
c
c    Input, double precision A_LU(3,N), the LU factors from R83_NP_FA.
c
c    Input/output, double precision B(N).
c    On input, B contains the right hand side of the linear system.
c    On output, B contains the solution of the linear system.
c
c    Input, integer JOB, specifies the system to solve.
c    0, solve A * x = b.
c    nonzero, solve A' * x = b.
c
      implicit none

      integer n

      double precision a_lu(3,n)
      double precision b(n)
      integer i
      integer job

      if ( job .eq. 0 ) then
c
c  Solve L * Y = B.
c
        do i = 2, n
          b(i) = b(i) - a_lu(3,i-1) * b(i-1)
        end do
c
c  Solve U * X = Y.
c
        do i = n, 1, -1
          b(i) = b(i) / a_lu(2,i)
          if ( 1 .lt. i ) then
            b(i-1) = b(i-1) - a_lu(1,i) * b(i)
          end if
        end do

      else
c
c  Solve U' * Y = B
c
        do i = 1, n
          b(i) = b(i) / a_lu(2,i)
          if ( i .lt. n ) then
            b(i+1) = b(i+1) - a_lu(1,i+1) * b(i)
          end if
        end do
c
c  Solve L' * X = Y.
c
        do i = n-1, 1, -1
          b(i) = b(i) - a_lu(3,i) * b(i+1)
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
      subroutine u0 ( a, b, t0, n, x, value )

c*********************************************************************72
c
cc U0 returns the initial condition at the starting time.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    30 May 2009
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, double precision A, B, the left and right endpoints
c
c    Input, double precision T0, the initial time.
c
c    Input, double precision T, the current time.
c
c    Input, integer N, the number of points where initial data is needed.
c
c    Input, double precision X(N), the positions where initial data is needed.
c
c    Output, double precision VALUE(N), the prescribed value of U(X,T0).
c
      implicit none

      integer n

      double precision a
      double precision b
      integer i
      double precision t
      double precision t0
      double precision value(n)
      double precision x(n)

      do i = 1, n
        value(i) = 100.0D+00
      end do

      return
      end
      subroutine ua ( a, b, t0, t, value )

c*********************************************************************72
c
cc UA returns the Dirichlet boundary condition at the left endpoint.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    30 May 2009
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, double precision A, B, the left and right endpoints
c
c    Input, double precision T0, the initial time.
c
c    Input, double precision T, the current time.
c
c    Output, double precision VALUE, the prescribed value of U(A,T).
c
      implicit none

      double precision a
      double precision b
      double precision t
      double precision t0
      double precision value

      value = 20.0D+00

      return
      end
      subroutine ub ( a, b, t0, t, value )

c*********************************************************************72
c
cc UB returns the Dirichlet boundary condition at the right endpoint.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    30 May 2009
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, double precision A, B, the left and right endpoints
c
c    Input, double precision T0, the initial time.
c
c    Input, double precision T, the current time.
c
c    Output, double precision VALUE, the prescribed value of U(B,T).
c
      implicit none

      double precision a
      double precision b
      double precision t
      double precision t0
      double precision value

      value = 20.0D+00

      return
      end
