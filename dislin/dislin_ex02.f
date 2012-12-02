      program main

c*********************************************************************72
c
cc DISLIN_EX02 demonstrates the use of POLAR to plot (R,Theta) data.
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

      integer m
      integer n
      parameter ( m = 10 )
      parameter ( n = 300 )

      real a
      integer i
      real pi
      real step
      real x2(m)
      real xray(n)
      real y2(m)
      real yray(n)

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'DISLIN_EX02:'
      write ( *, '(a)' ) '  FORTRAN77 version'
      write ( *, '(a)' ) '  Demonstrate the use of POLAR to plot '
      write ( *, '(a)' ) '  (R,Theta) data.'

      pi = 3.1415927
      step = 360.0 / real ( n - 1 )

      do i = 1, n
        a = real ( i - 1 ) * step
        a = a * pi / 180
        yray(i) = a
        xray(i) = sin ( 5 * a )
      end do

      do i = 1, m
        x2(i) = real ( i )
        y2(i) = real ( i )
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
      call setfil ( 'dislin_ex02.png' )
c
c  Choose the page size and orientation.
c
      call setpag ( 'usap' )
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
c  Use the HARDWARE font.
c
      call hwfont ( )

      call titlin ( 'Polar Plots', 2 )
      call ticks ( 3, 'Y' )
      call axends ( 'NOENDS', 'X' )
      call labdig ( -1, 'Y' )
      call axslen ( 1000, 1000 )
      call axsorg ( 1050, 900 )

      call polar ( 1.0, 0.0, 0.2, 0.0, 30.0 )
      call curve ( xray, yray, n )
      call htitle ( 50 )
      call title ( )
      call endgrf ( )

      call labdig ( -1, 'X' )
      call axsorg ( 1050, 2250 )
      call labtyp ( 'VERT', 'Y' )
      call barwth ( 5.0 )
      
      call polar ( 10.0, 0.0, 2.0, 0.0, 30.0 )
      call barwth ( -5.0 )
      call polcrv ( 'FBARS' )
      call curve ( x2, y2, m )
c
c  Close DISLIN.
c
      call disfin ( )
c
c  Terminate.
c
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'DISLIN_EX02:'
      write ( *, '(a)' ) '  Normal end of execution.'

      stop
      end
