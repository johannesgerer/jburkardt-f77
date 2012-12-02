      program main

c*********************************************************************72
c
cc DISLIN_EX01 demonstrates the use of CURVE to plot (X,Y) data.
c
c  Modified:
c
c    09 April 2011
c
c  Reference:
c
c    Helmut Michels,
c    The Data Plotting Software DISLIN - version 10.4,
c    Shaker Media GmbH, January 2010,
c    ISBN13: 978-3-86858-517-9.
c
      implicit none

      integer n
      parameter ( n = 301 )

      integer i
      real fpi
      real pi
      parameter ( pi = 3.1415926 )
      real step
      real x
      real xray(n)
      real y1ray(n)
      real y2ray(n)

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'DISLIN_EX01:'
      write ( *, '(a)' ) '  FORTRAN77 version'
      write ( *, '(a)' ) '  Demonstrate the use of CURVE to plot '
      write ( *, '(a)' ) '  (X,Y) data.'

      fpi = pi / 180.0E+00
      step = 360.0E+00 / real ( n - 1 )

      do i = 1, n
        xray(i) = real ( i - 1 ) * step
        x = xray(i) * fpi
        y1ray(i) = sin ( x )
        y2ray(i) = cos ( x )
      end do
c
c  Specify the format of the output file.
c
      call metafl ( 'png' )
c
c  Indicate that new data overwrites old data.
c
      call filmod ( 'delete' )
c
c  Specify the name of the output graphics file.
c
      call setfil ( 'dislin_ex01.png' )
c
c  Choose the page size and orientation.
c
      call setpag ( 'usal' )
c
c  For PNG output, reverse the default black background to white.
c
      call scrmod ( 'reverse' )
c
c  Open DISLIN.
c
      call disini ( )
c
c  Plot a border around the page.
c
      call pagera ( )
c
c  Use the COMPLEX font.
c
      call complx ( )

      call axspos ( 450, 1800 )
      call axslen ( 2200, 1200 )

      call name ( 'X-axis', 'X' )
      call name ( 'Y-axis', 'Y' )

      call labdig ( -1, 'X' )
      call ticks ( 10, 'XY' )

      call titlin ( 'Demonstration of CURVE', 1 )
      call titlin ( 'SIN(X), COS(X)', 3 )

      call graf ( 0.0, 360.0, 0.0, 90.0, -1.0, 1.0, -1.0, 0.5 )
      call title ( )

      call color ( 'RED' )
      call curve ( xray, y1ray, n )
      call color ( 'GREEN' )
      call curve ( xray, y2ray, n )

      call color ( 'FORE' )
      call dash ( )
      call xaxgit ( )
c
c  Close DISLIN.
c
      call disfin ( )
c
c  Terminate.
c
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'DISLIN_EX01:'
      write ( *, '(a)' ) '  Normal end of execution.'

      stop
      end
