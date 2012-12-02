      subroutine prime_number ( n, total )

c*********************************************************************72
c
cc PRIME_NUMBER returns the number of primes between 1 and N.
c
c  Discussion:
c
c    A naive algorithm is used.
c
c    Mathematica can return the number of primes less than or equal to N
c    by the command PrimePi[N].
c
c                N  PRIME_NUMBER
c
c                1           0
c               10           4
c              100          25
c            1,000         168
c           10,000       1,229
c          100,000       9,592
c        1,000,000      78,498
c       10,000,000     664,579
c      100,000,000   5,761,455
c    1,000,000,000  50,847,534
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    20 May 2009
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the maximum number to check.
c
c    Output, integer TOTAL, the number of prime numbers up to N.
c
      implicit none

      integer i
      integer j
      integer n
      integer prime
      integer total

      total = 0

      do i = 2, n

        prime = 1

        do j = 2, i - 1
          if ( mod ( i, j ) == 0 ) then
            prime = 0
            go to 10
          end if
        end do

10      continue

        total = total + prime

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
