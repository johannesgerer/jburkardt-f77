      program main

c*********************************************************************72
c
cc DOOMSDAY_PRB tests DOOMSDAY.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    28 May 2012
c
c  Author:
c
c    John Burkardt
c
      implicit none

      call timestamp ( )
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'DOOMSDAY_PRB:'
      write ( *, '(a)' ) '  FORTRAN77 version'
      write ( *, '(a)' ) '  Test the DOOMSDAY library.'

      call test01 ( )
      call test02 ( )
c
c  Terminate.
c
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'DOOMSDAY_PRB:'
      write ( *, '(a)' ) '  Test the DOOMSDAY library.'
      write ( *, '(a)' ) ' '
      call timestamp ( )

      stop
      end
      subroutine test01 ( )

c*********************************************************************72
c
cc TEST01 tests DOOMSDAY against a couple of test dates.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    27 May 2012
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer d
      integer m
      integer n_data
      integer w
      character * ( 10 ) s1
      character * ( 10 ) s2
      integer y

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST01'
      write ( *, '(a)' ) '  Try a couple selected dates.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  YYYY  MM  DD  Weekday    Weekday'
      write ( *, '(a)' ) '                Tabulated  Computed'
      write ( *, '(a)' ) ' '

      y = 1989
      m = 7
      d = 13
      call doomsday_gregorian ( y, m, d, w )
      call weekday_to_name_common ( w, s1 )
      s2 = 'Thursday'
      write ( *, '(2x,i4,2x,i2,2x,i2,2x,a10,2x,a10)' ) y, m, d, s1, s2

      y = 2012
      m = 5
      d = 26
      call doomsday_gregorian ( y, m, d, w )
      call weekday_to_name_common ( w, s1 )
      s2 = 'Saturday'
      write ( *, '(2x,i4,2x,i2,2x,i2,2x,a10,2x,a10)' ) y, m, d, s1, s2

      return
      end
      subroutine test02 ( )

c*********************************************************************72
c
cc TEST02 tests DOOMSDAY against a number of known values.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    28 May 2012
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer d
      integer m
      integer n_data
      integer w1
      integer w2
      character * ( 10 ) s1
      character * ( 10 ) s2
      integer y

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST02'
      write ( *, '(a)' ) 
     &  '  WEEKDAY_VALUES supplies a list of dates and weekdays.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  YYYY  MM  DD  Weekday    Weekday'
      write ( *, '(a)' ) '                Tabulated  Computed'
      write ( *, '(a)' ) ' '

      n_data = 0

      do

        call weekday_values ( n_data, y, m, d, w1 )

        if ( n_data <= 0 ) then
          exit
        end if
c
c  The transition from Julian to Gregorian calendars occurred in 1582
c  (for some people).  The data in "WEEKDAY_VALUES" before the transition
c  is stored in Julian format, which DOOMSDAY_GREGORIAN can't handle.
c  So let's just refuse to handle 1582 or earlierc
c
        if ( y <= 1582 ) then
          cycle
        end if

        call doomsday_gregorian ( y, m, d, w2 )

        call weekday_to_name_common ( w1, s1 )
        call weekday_to_name_common ( w2, s2 )

        write ( *, '(2x,i4,2x,i2,2x,i2,2x,a10,2x,a10)' ) y, m, d, s1, s2

      end do

      return
      end
