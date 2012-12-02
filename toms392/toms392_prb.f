      program main

c*********************************************************************72
c
cc TOMS392_PRB tests TOMS392.
c
c  Modified:
c
c    11 January 2006
c
c  Author:
c
c    John Burkardt
c
      implicit none

      call timestamp ( )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TOMS392_PRB'
      write ( *, '(a)' ) '  FORTRAN77 version'
      write ( *, '(a)' ) '  Test the TOMS392 library.'

      call testch ( )
c
c  Terminate.
c
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TOMS392_PRB'
      write ( *, '(a)' ) '  Normal end of execution.'

      write ( *, '(a)' ) ' '
      call timestamp ( )

      stop
      end
      subroutine testch ( )

c*********************************************************************72
c
cc TESTCH is the test program for CHARAC.
c
      implicit none

      real data(4,81)
      real fm
      integer i
      integer ifail
      integer m
      integer n

      n = 20
c
c  generate initial data.
c
      m = 4 * n + 1
      fm = 4.0e+00 * float ( n )
      do i = 1, m
        data(1,i) = float ( i - 1 ) / fm
        data(2,i) = 0.0e+00
        data(3,i) = 0.0e+00
        data(4,i) = 2.0e+00 * exp ( data(1,i) )
      end do

      ifail = 0
      write ( *, 900 )
900   format ( 1h1 )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' )
     &  '       x               y               u               v'
      write ( *, '(a)' ) ' '

200   continue

      do i = 1, m
        write ( *, 910 ) data(1,i), data(2,i), data(3,i), data(4,i)
      end do

910   format ( 2x,g14.6,2x,g14.6,2x,g14.6,2x,g14.6 )

      if ( m .le. 1 ) go to 300
      if ( ifail .ne. 0 ) go to 300

      call charac ( data, m, ifail )
      write ( *, 900 )

      go to 200

300   continue
      write ( *, 920 ) m, ifail
920   format ( 1x, 3hm =,i2,8h ifail =,i2 )
      return
      end
      subroutine chcoef ( coeff, xyuv )

c*********************************************************************72
c
cc CHCOEF computes coefficients a1, a2, a3, a4, h1, b1, b2, b3, b4, h2,
c  and stores them sequentially in coeff.
c
      implicit none

      real coeff(10)
      real xyuv(4)

      coeff(1) = 1.0e+00 - xyuv(3)**2
      coeff(2) = -xyuv(3) * xyuv(4)
      coeff(3) = -xyuv(3) * xyuv(4)
      coeff(4) = 1.0e+00 - xyuv(4)**2
      coeff(5) = -4.0e+00 * xyuv(3) * exp ( xyuv(1) )**2
      coeff(6) = 0.0e+00
      coeff(7) = 1.0e+00
      coeff(8) = -1.0e+00
      coeff(9) = 0.0e+00
      coeff(10) = 0.0e+00

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
