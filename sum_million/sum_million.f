      program main

c*********************************************************************72
c
cc MAIN is the main program for SUM_MILLION.
c
c  Discussion:
c
c    This code estimates the power of a computer by summing the integers
c    from 1 to 1,000,000.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    18 July 2008
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n
      parameter ( n = 1000000 )

      double precision ctime
      double precision error
      double precision exact
      parameter ( exact = 500000500000.0D+00 )
      integer i
      double precision mflops
      double precision total
      double precision x(n)

      write ( *, '(a)' ) ' '
      call timestamp ( )
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'SUM_MILLION'
      write ( *, '(a)' ) '  FORTRAN77 version'
      write ( *, '(a)' ) '  Sum the integers from 1 to 1,000,000.'
      write ( *, '(a)' ) '  Correct answer is 500000500000.'

      call set_up ( n, x )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 
     &  '         N      CPU time        MFLOPS     ERROR'
      write ( *, '(a)' ) '                (seconds)'
      write ( *, '(a)' ) ' '

      do i = 1, 10

        call sum_up ( n, x, total, ctime )

        mflops = dble ( n ) / 1000000.0D+00 / ctime

        error = total - exact

        write ( *, '(2x,i8,2x,g14.6,2x,g14.6,2x,g14.6)' ) 
     &    n, ctime, mflops, error

      end do

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'SUM_MILLION:'
      write ( *, '(a)' ) '  Normal end of execution.'
      write ( *, '(a)' ) ' '
      call timestamp ( )

      stop
      end
      subroutine set_up ( n, x )

c*********************************************************************72
c
cc SET_UP sets up the data for the SUM_MILLION program.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    24 July 2008
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the number of values to define.
c
c    Output, double precision X(N), a vector which contains the values 1 through N.
c
      implicit none

      integer n

      integer i
      double precision x(n)

      x(1) = 1.0
      do i = 2, n
        x(i) = x(i-1) + 1.0
      end do

      return
      end
      subroutine sum_up ( n, x, total, ctime )

c*********************************************************************72
c
cc SUM_UP carries out the sum for the SUM_MILLION program.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    24 July 2008
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the number of values to define.
c
c    Input, double precision X(N), the data to be summed.
c
c    Output, double precision TOTAL, the sum of the data.
c
c    Output, double precision CTIME, the cpu time required to sum the data.
c
      implicit none

      integer n

      double precision ctime
      double precision ctime1
      double precision ctime2
      integer i
      double precision total
      double precision x(n)

      call cpu_time ( ctime1 )

      total = 0.0
      do i = 1, n
        total = total + x(i)
      end do

      call cpu_time ( ctime2 )

      ctime = ctime2 - ctime1

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
