      program main

c*********************************************************************72
c
cc MAIN is the main program for ARGS.
c
c  Discussion:
c
c    ARGS demonstrates the use of the (semi-standard) GETARGS utility.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    24 January 2009
c
c  Author:
c
c    John Burkardt
c
c  Modified:
c
c    26 April 2007
c
c  Author:
c
c    John Burkardt
c
c  Usage:
c
c    args arg1 arg2 arg3 ...
c
      implicit none

      character * ( 80 ) arg
      integer i
      integer iargc
      integer numarg

      call timestamp ( )
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'ARGS'
      write ( *, '(a)' ) '  FORTRAN77 version'
      write ( *, '(a)' ) '  Demonstrate the use of the command line'
      write ( *, '(a)' ) '  argument routines GETARG and IARGC.'

      numarg = iargc ( )

      write ( *, '(a)' ) ' '
      write ( *, '(a,i8,a)' ) 
     &  '  ARGS was called with IARGC() = ', numarg, ' arguments.'

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  CALL GETARG(I,ARG) returns the arguments:'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  I     ARG '
      write ( *, '(a)' ) ' '

      do i = 0, numarg
        call getarg ( i, arg )
        write ( *, '(2x,i3,2x,a)' ) i, arg
      end do

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'ARGS:'
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
c    24 January 2009
c
c  Author:
c
c    John Burkardt
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
