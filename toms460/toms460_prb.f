      program main

c*********************************************************************72
c
cc TOMS460_PRB tests ADIP.
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

      integer n
      parameter ( n = 10 )

      double precision a
      double precision b
      double precision c
      double precision d
      integer i
      integer ier
      integer iopt
      integer itns
      double precision dmu
      double precision omeh(n)
      double precision omev(n)

      call timestamp ( )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TOMS460_PRB'
      write ( *, '(a)' ) '  FORTRAN77 version'
      write ( *, '(a)' ) '  Test the TOMS460 library.'
      write ( *, '(a)' ) ' '

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Given ITNS, compute DMU:'
      write ( *, '(a)' ) ' '

      a = 0.5D+00
      b = 2.0D+00
      c = 0.25D+00
      d = 0.75D+00
      iopt = 1

      do itns = 1, 5

        call adip ( a, b, c, d, iopt, n, itns, dmu, omeh,
     &    omev, ier )

        write ( *, '(2x,i6,2x,g14.6,2x,i6)' ) itns, dmu, ier

      end do

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Given DMU, compute ITNS:'
      write ( *, '(a)' ) ' '

      a = 0.5D+00
      b = 2.0D+00
      c = 0.25D+00
      d = 0.75D+00
      iopt = 2
      dmu = 1.0D+00

      do i = 1, 8

        dmu = dmu / 10.0D+00

        call adip ( a, b, c, d, iopt, n, itns, dmu, omeh,
     &    omev, ier )

         write ( *, '(2x,g14.6,2x,i6,2x,i8)' ) dmu, itns, ier

      end do
c
c  Terminate.
c
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TOMS460_PRB'
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
