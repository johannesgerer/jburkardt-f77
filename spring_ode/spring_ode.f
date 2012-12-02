      program main

c*********************************************************************72
c
cc MAIN is the main program for SPRING_ODE.
c
c  Discussion:
c
c    This is a simple example of how to plot when you don't have a plotter.
c    This is a particular kind of "ASCII graphics", or "typewriter graphics"
c    or "lineprinter graphics", and shows you how valuable an illustration 
c    can be, even when it's as crude as this example.
c
c    Hooke's law for a spring observes that the restoring force is
c    proportional to the displacement: F = - k x
c
c    Newton's law relates the force to acceleration: F = m a
c
c    Putting these together, we have
c
c      m * d^2 x/dt^2 = - k * x
c
c    We can add a damping force with coefficient c:
c
c      m * d^2 x/dt^2 = - k * x - c * dx/dt
c
c    If we write this as a pair of first order equations for (x,v), we have
c
c          dx/dt = v
c      m * dv/dt = - k * x - c * v
c
c    and now we can approximate these values for small time steps.
c
c    Note that the plotting assumes that the value of X will always be
c    between -1 and +1.  If the initial condition uses V = 0, and X starts
c    between -1 and +1, then this will be OK.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    15 May 2012
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

      double precision c
      double precision dt
      integer i
      integer j
      double precision k
      double precision m
      integer n
      integer p
      double precision t
      double precision t_final
      double precision v
      double precision v_old
      double precision x
      double precision x_old
      character z(21)

      call timestamp ( )
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'SPRING_ODE'
      write ( *, '(a)' ) '  FORTRAN77 version'
      write ( *, '(a)' ) 
     &  '  Approximate the solution of a spring equation.'
      write ( *, '(a)' ) 
     &  '  Display the solution with line printer graphics.'
      write ( *, '(a)' ) ' '
c
c  Data
c
      m = 1.0
      k = 1.0
      c = 0.3
      t_final = 20.0
      n = 100
      dt = t_final / dble ( n )
c
c  Initial conditions.
c
      x = 1.0
      v = 0.0
c
c  Compute the approximate solution at equally spaced times.
c
      do i = 0, n

        x_old = x
        v_old = v

        t = dble ( i ) * t_final / dble ( n )
        x = x_old + dt * v_old
        v = v_old + ( dt / m ) * ( - k * x_old - c * v_old )
c
c  Approximate the position of X in [-1,+1] to within 1/10.
c
        p = int ( ( (  1 * ( 1.0 - x ) ) + 21 * ( 1.0 + x ) ) / 2.0 ) 
        p = max ( p, 1 )
        p = min ( p, 21 )
c
c  Fill in the next line of the plot, placing 'x' in the p position.
c
        if ( mod ( i, 10 ) .eq. 0 ) then
          z(1:21) = '-'
        else
          z(1:21) = ' '
        end if
        z(1) = '|'
        z(6) = '.'
        z(11) = '+'
        z(16) = '.'
        z(21) = '|'

        z(p) = 'x'

        write ( *, '(21a)' ) z(1:21)

      end do
c
c  Terminate.
c
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'SPRING_ODE:'
      write ( *, '(a)' ) '  Normal end of execution.'
      write ( *, '(a)' ) ' '
      call timestamp ( )

      stop
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
