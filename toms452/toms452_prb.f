      program main

c*********************************************************************72
c
cc TOMS452_PRB tests NXCBN.
c
c  Modified:
c
c    12 January 2006
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer i
      integer ic(10)
      integer j
      integer m
      integer n
      integer number

      call timestamp ( )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TOMS452_PRB'
      write ( *, '(a)' ) '  FORTRAN77 version'
      write ( *, '(a)' ) '  Test the TOMS452 library.'

      n = 10
      m = 3
      number = 5

      do i = 1, n
        if ( i <= m ) then
          ic(i) = 1
        else
          ic(i) = 0
        end if
      end do

      write ( *, '(a)' ) ' '
      write ( *, '(a,i2)' ) '  Let N = ', n
      write ( *, '(a,i2)' ) '  and M = ', m
      write ( *, '(a,i2,a)' ) '  Generate ', number, ' combinations.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Initial combination:'
      write ( *, '(a)' ) ' '
      write ( *, '(2x,10i1)' ) ( ic(i), i = 1, n )
      write ( *, '(a)' ) ' '

      do j = 1, number
        call nxcbn ( n, m, ic )
        write ( *, '(2x,10i1)' ) ( ic(i), i = 1, n )
      end do

      n = 5
      m = 2
      number = 15

      do i = 1, n
        if ( i <= m ) then
          ic(i) = 1
        else
          ic(i) = 0
        end if
      end do

      write ( *, '(a)' ) ' '
      write ( *, '(a,i2)' ) '  Let N = ', n
      write ( *, '(a,i2)' ) '  and M = ', m
      write ( *, '(a,i2,a)' ) '  Generate ', number, ' combinations.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Initial combination:'
      write ( *, '(a)' ) ' '
      write ( *, '(2x,10i1)' ) ( ic(i), i = 1, n )
      write ( *, '(a)' ) ' '

      do j = 1, number
        call nxcbn ( n, m, ic )
        write ( *, '(2x,10i1)' ) ( ic(i), i = 1, n )
      end do
c
c  Terminate.
c
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TOMS452_PRB'
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
c  Modified:
c
c    16 September 2005
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

      character ( len = 8 ) date
      character ( len = 10 ) time

      call date_and_time ( date, time )

      write ( *, '(a8,2x,a10)' ) date, time

      return
      end
