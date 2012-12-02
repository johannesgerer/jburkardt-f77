      program main

c*********************************************************************72
c
cc MAIN is the main program for DIVISION.
c
c  Discussion:
c
c    DIVISION demonstrates that division should be done to full precision.
c
c    Especially when carrying out division, but also in other arithmetic
c    operations, it is important that constants be specified as double
c    precision values, by using the "D" identifier in the exponent field.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    04 December 20080
c
c  Author:
c
c    John Burkardt
c
      implicit none

      double precision a(4,4)
      double precision b(4,4)
      integer i
      integer j

      call timestamp ( )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'DIVISION:'
      write ( *, '(a)' ) '  FORTRAN77 version.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Demonstrate 16 ways to compute 1/3,'
      write ( *, '(a)' ) '  expecting a double precision real result.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  The results, stored as a table, are:'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 
     &  '  1       / 3   1       / 3.0   1       / 3.0E+00' //
     &  '  1       / 3.0D+00'
      write ( *, '(a)' ) 
     &  '  1.0     / 3   1.0     / 3.0   1.0     / 3.0E+00' //
     &  '  1.0     / 3.0D+00'
      write ( *, '(a)' ) 
     &  '  1.0E+00 / 3   1.0E+00 / 3.0   1.0E+00 / 3.0E+00' //
     &  '  1.0E+00 / 3.0D+00'
      write ( *, '(a)' ) 
     &  '  1.0D+00 / 3   1.0D+00 / 3.0   1.0E+00 / 3.0E+00' //
     &  '  1.0D+00 / 3.0D+00'

      a(1,1) = 1 / 3
      a(1,2) = 1 / 3.0
      a(1,3) = 1 / 3.0E+00
      a(1,4) = 1 / 3.0D+00

      a(2,1) = 1.0 / 3
      a(2,2) = 1.0 / 3.0
      a(2,3) = 1.0 / 3.0E+00
      a(2,4) = 1.0 / 3.0D+00

      a(3,1) = 1.0E+00 / 3
      a(3,2) = 1.0E+00 / 3.0
      a(3,3) = 1.0E+00 / 3.0E+00
      a(3,4) = 1.0E+00 / 3.0D+00

      a(4,1) = 1.0D+00 / 3
      a(4,2) = 1.0D+00 / 3.0
      a(4,3) = 1.0D+00 / 3.0E+00
      a(4,4) = 1.0D+00 / 3.0D+00

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  The matrix of results:'
      write ( *, '(a)' ) ' '

      do i = 1, 4
        write ( *, '(4g20.14)' ) ( a(i,j), j = 1, 4 )
      end do

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Defects might be more obvious if we'
      write ( *, '(a)' ) '  multiply by 3:'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  B = 3 * A:'
      write ( *, '(a)' ) ' '

      do j = 1, 4
        do i = 1, 4
          b(i,j) = 3 * a(i,j)
        end do
      end do

      do i = 1, 4
        write ( *, '(4g20.14)' ) ( b(i,j), j = 1, 4 )
      end do

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  B = 3.0 * A:'
      write ( *, '(a)' ) ' '

      do j = 1, 4
        do i = 1, 4
          b(i,j) = 3.0 * a(i,j)
        end do
      end do

      do i = 1, 4
        write ( *, '(4g20.14)' ) ( b(i,j), j = 1, 4 )
      end do

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  B = 3.0E+00 * A:'
      write ( *, '(a)' ) ' '

      do j = 1, 4
        do i = 1, 4
          b(i,j) = 3.0E+00 * a(i,j)
        end do
      end do

      do i = 1, 4
        write ( *, '(4g20.14)' ) ( b(i,j), j = 1, 4 )
      end do

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  B = 3.0D+00 * A:'
      write ( *, '(a)' ) ' '

      do j = 1, 4
        do i = 1, 4
          b(i,j) = 3.0D+00 * a(i,j)
        end do
      end do

      do i = 1, 4
        write ( *, '(4g20.14)' ) ( b(i,j), j = 1, 4 )
      end do

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'DIVISION:'
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
