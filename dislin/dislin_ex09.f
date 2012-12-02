      program main

c*********************************************************************72
c
cc DISLIN_EX09 demonstrates 3D color plots.
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
      parameter ( n = 100 )

      real fpi
      integer i
      integer j
      real step
      real x
      real y
      real zmat(n,n)

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'DISLIN_EX09:'
      write ( *, '(a)' ) '  FORTRAN77 version'
      write ( *, '(a)' ) '  Demonstrate the creation of a '
      write ( *, '(a)' ) '  3D color plot.'

      fpi = 3.1415927 / 180.0
      step = 360.0 / real ( n - 1 )

      do i = 1, n
        x = real ( i - 1 ) * step
        do j = 1, n
          y = real ( j - 1 ) * step
          zmat(i,j) = 2.0 * sin ( x * fpi ) * sin ( y * fpi )
        end do
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
      call setfil ( 'dislin_ex09.png' )
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
c  Use the HARDWARE font.
c
      call hwfont ( )

      call titlin ( '3-D Color Plot of the Function', 1 )
      call titlin ( 'F(X,Y) = 2 * SIN(X) * SIN(Y)', 3 )

      call name ( 'X-axis', 'X' )
      call name ( 'Y-axis', 'Y' )
      call name ( 'Z-axis', 'Z' )

      call intax ( )
      call autres ( n, n )
      call axspos ( 300, 1850 )
      call ax3len ( 2200, 1400, 1400 )

      call graf3 ( 0.0, 360.0, 0.0, 90.0, 0.0, 360.0, 0.0, 90.0,
     &  -2.0, 2.0, -2.0, 1.0 )

      call crvmat ( zmat, n, n, 1, 1 )

      call height ( 50 )
      call title ( )
      call mpaepl ( 3 )
c
c  Close DISLIN.
c
      call disfin ( )
c
c  Terminate.
c
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'DISLIN_EX09:'
      write ( *, '(a)' ) '  Normal end of execution.'

      stop
      end
