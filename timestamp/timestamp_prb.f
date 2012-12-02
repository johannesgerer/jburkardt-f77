       program main

c*********************************************************************72
c
cc MAIN is the main program for TIMESTAMP_PRB.
c
c  Discussion:
c
c    TIMESTAMP_PRB demonstrates the use of TIMESTAMP.
c
c  Modified:
c
c    12 January 2007
c
c  Author:
c
c    John Burkardt
c
       implicit none

       call timestamp ( )

       write ( *, '(a)' ) ' '
       write ( *, '(a)' ) 'TIMESTAMP_PRB'
       write ( *, '(a)' ) '  FORTRAN77 version'
       write ( *, '(a)' ) '  Test the TIMESTAMP library.'

       call test01 ( )
       call test02 ( )
c
c  Terminate.
c
       write ( *, '(a)' ) ' '
       write ( *, '(a)' ) 'TIMESTAMP_PRB'
       write ( *, '(a)' ) '  Normal end of execution.'

       write ( *, '(a)' ) ' '
       call timestamp ( )

       stop
       end
       subroutine test01 ( )

c*********************************************************************72
c
cc TEST01 demonstrates the use of TIMESTAMP.
c
c  Modified:
c
c    12 January 2007
c
c  Author:
c
c    John Burkardt
c
       implicit none

       write ( *, '(a)' ) ' '
       write ( *, '(a)' ) 'TEST01'
       write ( *, '(a)' ) '  TIMESTAMP prints out the current wallclock'
       write ( *, '(a)' ) '  time, including the year, month, day, '
       write ( *, '(a)' ) '  hours, minutes, seconds, thousandths of a '
       write ( *, '(a)' ) '  second, and AM/PM.'
       write ( *, '(a)' ) ' '
       write ( *, '(a)' ) '  This can be useful in keeping track of the'
       write ( *, '(a)' ) '  date of execution of a particular program'
       write ( *, '(a)' ) '  or to give a rough idea of the length of'
       write ( *, '(a)' ) '  time required to run a program.'

       write ( *, '(a)' ) ' '
       call timestamp ( )

       return
       end
       subroutine test02 ( )

c*********************************************************************72
c
cc TEST02 demonstrates the use of TIMESTRING.
c
c  Modified:
c
c    12 January 2007
c
c  Author:
c
c    John Burkardt
c
      implicit none

      character * (40) string

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST02'
      write ( *, '(a)' )
     &  '  TIMESTRING returns the current wallclock time,'
      write ( *, '(a)' )
     &  '  including the year, month, day, hours, minutes,'
      write ( *, '(a)' )
     &  '  seconds, thousandths of a second, and AM/PM'
      write ( *, '(a)' )
     &  '  in a string, which the user may print or manipulate.'

      call timestring ( string )

      write ( *, '(a)' ) ' '
      write ( *, '(a,a,a)' ) '  TIMESTRING returned the value "'
     &   , string, '".'

      return
      end
