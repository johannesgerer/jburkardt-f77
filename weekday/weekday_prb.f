      program main

c*********************************************************************72
c
cc MAIN is the main program for WEEKDAY_PRB.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    26 May 2012
c
c  Author:
c
c    John Burkardt
c
      implicit none

      call timestamp ( )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'WEEKDAY_PRB:'
      write ( *, '(a)' ) '  FORTRAN77 version'
      write ( *, '(a)' ) '  Test the WEEKDAY library.'

      call test01 ( )
c
c  Terminate.
c
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'WEEKDAY_PRB:'
      write ( *, '(a)' ) '  Noraml end of execution.'

      write ( *, '(a)' ) ' '
      call timestamp ( )

      stop
      end
      subroutine test01 ( )

c*********************************************************************72
c
cc TEST01 tests YMD_TO_WEEKDAY_COMMON.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    26 May 2012
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer d
      integer m
      integer n_data
      character * ( 9 ) s1
      character * ( 9 ) s2
      character * ( 20 ) s3
      integer w1
      integer w2
      integer y

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST01'
      write ( *, '(a)' ) '  For dates in the Common calendar:'
      write ( *, '(a)' ) 
     &  '  YMD_TO_WEEKDAY_COMMON returns the day of the week.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  YMD                   Weekday    Weekday'
      write ( *, '(a)' ) '                        Tabulated  Computed'
      write ( *, '(a)' ) ' '

      do

        call weekday_values ( n_data, y, m, d, w1 )

        if ( n_data .eq. 0 ) then
          exit
        end if
     
        call ymd_to_s_common ( y, m, d, s3 ) 
        call ymd_to_weekday_common ( y, m, d, w2 )
        call weekday_to_name_common ( w1, s1 )
        call weekday_to_name_common ( w2, s2 )

        write ( *, '(2x,a20,2x,a9,2x,a9)' ) s3, s1, s2

      end do

      return
      end
