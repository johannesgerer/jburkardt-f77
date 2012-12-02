      program main

c*********************************************************************72
c
cc FD1D_BURGERS_LEAP solves the nonviscous Burgers equation using leapfrogging.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    18 August 2010
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

      integer n
      parameter ( n = 21 )

      double precision a
      double precision b
      double precision dt
      double precision dx
      integer i
      integer ihi
      integer ilo
      integer step
      integer step_num
      double precision t
      double precision t_init
      double precision t_last
      double precision uc(n)
      double precision un(n)
      double precision uo(n)
      double precision x(n)

      call timestamp ( )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'FD1D_BURGERS_LEAP:'
      write ( *, '(a)' ) '  FORTRAN77 version'
      write ( *, '(a)' ) 
     &  '  Solve the non-viscous time-dependent Burgers equation,'
      write ( *, '(a)' ) '  using the leap-frog method.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Equation to be solved:'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '    du/dt + u * du/dx = 0'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 
     &  '  for x in [ a, b ], for t in [t_init, t_last]'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  with initial conditions:'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '    u(x,o) = u_init'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  and boundary conditions:'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '    u(a,t) = u_a(t), u(b,t) = u_b(t)'
c
c  Set and report the problem parameters.
c

      a = -1.0D+00
      b = +1.0D+00
      dx = ( b - a ) / dble ( n - 1 )
      step_num = 30
      t_init = 0.0D+00
      t_last = 3.0D+00
      dt = ( t_last - t_init ) / dble ( step_num )

      write ( *, '(a)' ) ' '
      write ( *, '(2x,g14.6,a,g14.6)' ) a, ' <= X <= ', b
      write ( *, '(a,i8)' ) '  Number of nodes = ', n
      write ( *, '(a,g14.6)' ) '  DX = ', dx
      write ( *, '(a)' ) ' '
      write ( *, '(2x,g14.6,a,g14.6)' ) t_init, ' <= T <= ', t_last
      write ( *, '(a,i8)' ) '  Number of time steps = ', step_num
      write ( *, '(a,g14.6)' ) '  DT = ', dt

      call r8vec_even ( n, a, b, x )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  X:'
      write ( *, '(a)' ) ' '
      do ilo = 1, n, 5
        ihi = min ( ilo + 4, n )
        write ( *, '(2x,5g14.6)' ) x(ilo:ihi)
      end do
c
c  Set the initial condition,
c  and apply boundary conditions to first and last entries.
c
      step = 0
      t = t_init
      call u_init ( n, x, t, un )
      call u_a ( x(1), t, un(1) )
      call u_b ( x(n), t, un(n) )

      call report ( step, step_num, n, x, t, un )
c
c  Use Euler's method to get the first step.
c
      step = 1
      t = ( dble ( step_num - step ) * t_init   
     &    + dble (            step ) * t_last ) 
     &    / dble ( step_num        )

      do i = 1, n
        uc(i) = un(i)
      end do

      do i = 2, n - 1
        un(i) = uc(i) - dt * uc(i) * ( uc(i+1) - uc(i-1) ) 
     &    / 2.0D+00 / dx
      end do
      call u_a ( x(1), t, un(1) )
      call u_b ( x(n), t, un(n) )

      call report ( step, step_num, n, x, t, un )
c
c  Subsequent steps use the leapfrog method.
c
      do step = 2, step_num
     
        t = ( dble ( step_num - step ) * t_init   
     &      + dble (            step ) * t_last ) 
     &      / dble ( step_num        )

        do i = 1, n
          uo(i) = uc(i)
          uc(i) = un(i)
        end do

        do i = 2, n - 1
          un(i) = uo(i) - dt * uc(i) * ( uc(i+1) - uc(i-1) ) / dx
        end do

        call u_a ( x(1), t, un(1) )
        call u_b ( x(n), t, un(n) )

        call report ( step, step_num, n, x, t, un )

      end do
c
c  Terminate.
c
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'FD1D_BURGERS_LEAP:'
      write ( *, '(a)' ) '  Normal end of execution.'

      write ( *, '(a)' ) ' '
      call timestamp ( )

      stop
      end
      subroutine r8vec_even ( n, alo, ahi, a )

c*********************************************************************72
c
cc R8VEC_EVEN returns an R8VEC of evenly spaced values.
c
c  Discussion:
c
c    An R8VEC is a vector of R8's.
c
c    If N is 1, then the midpoint is returned.
c
c    Otherwise, the two endpoints are returned, and N-2 evenly
c    spaced points between them.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    09 December 2004
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the number of values.
c
c    Input, double precision ALO, AHI, the low and high values.
c
c    Output, double precision A(N), N evenly spaced values.
c    Normally, A(1) = ALO and A(N) = AHI.
c    However, if N = 1, then A(1) = 0.5*(ALO+AHI).
c
      implicit none

      integer n

      double precision a(n)
      double precision ahi
      double precision alo
      integer i

      if ( n .eq. 1 ) then

        a(1) = 0.5D+00 * ( alo + ahi )

      else

        do i = 1, n
          a(i) = ( dble ( n - i     ) * alo   
     &           + dble (     i - 1 ) * ahi ) 
     &           / dble ( n     - 1 )
        end do

      end if

      return
      end
      subroutine report ( step, step_num, n, x, t, u )

c*********************************************************************72
c
cc REPORT prints or plots or saves the data at the current time step.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    15 August 2010
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer STEP, the index of the current step,
c    between 0 and STEP_NUM.
c
c    Input, integer STEP_NUM, the number of steps to take.
c
c    Input, integer N, the number of nodes.
c
c    Input, double precision X(N), the coordinates of the nodes.
c
c    Input, double precision T, the current time.
c
c    Input, double precision U(N), the initial values U(X,T).
c
      implicit none

      integer n

      integer ihi
      integer ilo
      integer step
      integer step_num
      double precision t
      double precision u(n)
      double precision x(n)

      write ( *, '(a)' ) ' '
      write ( *, '(a,i8)' ) '  STEP = ', step
      write ( *, '(a,g14.6)' ) '  TIME = ', t
      write ( *, '(a)' ) ' '
      do ilo = 1, n, 5
        ihi = min ( ilo + 4, n )
        write ( *, '(2x,5g14.6)' ) u(ilo:ihi)
      end do

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
      subroutine u_a ( x, t, ua )

c*********************************************************************72
c
cc U_A sets the boundary condition for U at A.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    15 August 2010
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, double precision X, T, the position and time.
c
c    Output, double precision UA, the prescribed value of U(X,T).
c
      implicit none

      double precision t
      double precision ua
      double precision x

      ua = + 0.5D+00

      return
      end
      subroutine u_b ( x, t, ub )

c*********************************************************************72
c
cc U_B sets the boundary condition for U at B.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    15 August 2010
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, double precision X, T, the position and time.
c
c    Output, double precision UB, the prescribed value of U(X,T).
c
      implicit none

      double precision t
      double precision ub
      double precision x

      ub = - 0.5D+00

      return
      end
      subroutine u_init ( n, x, t, u )

c*********************************************************************72
c
cc U_INIT sets the initial condition for U.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    18 August 2010
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the number of nodes.
c
c    Input, double precision X(N), the coordinates of the nodes.
c
c    Input, double precision T, the current time.
c
c    Output, double precision U(N), the initial values U(X,T).
c
      implicit none

      integer n

      integer i
      double precision pi
      parameter ( pi = 3.141592653589793D+00 )
      double precision q
      double precision r
      double precision s
      double precision t
      double precision u(n)
      double precision ua
      double precision ub
      double precision x(n)

      call u_a ( x(1), t, ua )
      call u_b ( x(n), t, ub )

      q = 2.0D+00 * ( ua - ub ) / pi
      r = ( ua + ub ) / 2.0D+00
c
c  S can be varied.  It is the slope of the initial condition at the midpoint.
c
      s = 1.0D+00

      do i = 1, n
        u(i) = - q * atan ( s * ( 2.0D+00 * x(i) - x(1) - x(n) ) 
     &    / ( x(n) - x(1) ) ) + r
      end do

      return
      end

