      program main

c*********************************************************************72
c
cc QUICKPLOT_COLOR demonstrates the DISLIN quickplot command QPLCLR.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    13 May 2012
c
c  Author:
c
c    John Burkardt
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
      real pi
      parameter ( pi = 3.1415926 )
      real step
      real x
      real y
      real z(n,n)

      fpi = pi / 180.0

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'QUICKPLOT_COLOR:'
      write ( *, '(a)' ) '  FORTRAN77 version'
      write ( *, '(a)' ) '  Demonstrate the DISLIN "quickplot" command'
      write ( *, '(a)' ) 
     &  '  QPLCLR to make a color plot of a matrix of data.'
c
c  Set up the X and Y data for the plot.
c
      step = 360.0 / real ( n - 1 )
      do i = 1, n
        x = real ( i - 1 ) * step
        do j = 1, n
          y = real ( j - 1 ) * step
          z(i,j) = 2.0 * sin ( x * fpi ) * sin ( y * fpi )
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
      call setfil ( 'quickplot_color.png' )
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
c  Label the axes and the plot.
c
      call name ( "X-axis", "x" );
      call name ( "Y-axis", "y" );
      call titlin ( 'Quick plot by QPLCLR', 2 )
c
c  Draw the curve.
c
      call qplclr ( z, n, n )
c
c  Close DISLIN.
c
      call disfin ( )
c
c  Terminate.
c
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'DISLIN_COLOR:'
      write ( *, '(a)' ) '  Normal end of execution.'

      stop
      end
