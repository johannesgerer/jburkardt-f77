      program main

c*********************************************************************72
c
cc MAIN is the main program for SPRING_ODE2.
c
c  Discussion:
c
c    This is a revision of the SPRING_ODE code.
c
c    In this revision of the program, we want to use vectors (C arrays) to 
c    store the data, and we want to write the data out to a file in a form 
c    that Gnuplot (or other plotting programs) can use.
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
c    26 June 2012
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
      parameter ( n = 100 )

      double precision c
      double precision dt
      integer i
      double precision k
      double precision m
      double precision t(0:n)
      double precision t_final
      double precision v(0:n)
      double precision x(0:n)

      call timestamp ( )
      write ( *, '(a)' ) '#'
      write ( *, '(a)' ) '#SPRING_ODE2'
      write ( *, '(a)' ) '#  FORTRAN77 version'
      write ( *, '(a)' ) 
     &  '#  Approximate the solution of a spring equation.'
      write ( *, '(a)' ) '#  Write data to a file for use by gnuplot.'
      write ( *, '(a)' ) '#'
c
c  Data
c
      m = 1.0D+00
      k = 1.0D+00
      c = 0.3D+00
      t_final = 20.0D+00
      dt = t_final / dble ( n )
c
c  Initial conditions.
c
      t(0) = 0.0D+00
      x(0) = 1.0D+00
      v(0) = 0.0D+00
c
c  Compute the approximate solution at equally spaced times.
c
      do i = 1, n

        t(i) = dble ( i ) * t_final / dble ( n )
        x(i) = x(i-1) + dt * v(i-1)
        v(i) = v(i-1) + ( dt / m ) * ( - k * x(i-1) - c * v(i-1) )

      end do
c
c  Write the data to a file for plotting, possibly by gnuplot.
c  gnuplot expects T, X, and V to be columns of output.
c
      do i = 0, n
        write ( *, '(g14.6,2x,g14.6,2x,g14.6)' ) t(i), x(i), v(i)
      end do
c
c  Terminate.
c
      write ( *, '(a)' ) '#'
      write ( *, '(a)' ) '#SPRING_ODE2:'
      write ( *, '(a)' ) '#  Normal end of execution.'
      write ( *, '(a)' ) '#'
      call timestamp ( )

      stop
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
     &  '(a1,1x,i2,1x,a,1x,i4,2x,i2,a1,i2.2,a1,i2.2,a1,i3.3,1x,a)' ) 
     &  '#', d, month(m), y, h, ':', n, ':', s, '.', mm, ampm

      return
      end
