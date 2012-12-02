      program main

c*********************************************************************72
c
cc QUICKPLOT_CURVE demonstrates the DISLIN quickplot command QPLOT.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    13 April 2011
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

      integer i
      real pi
      parameter ( pi = 3.1415926 )
      real xray(n)
      real yray(n)

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'QUICKPLOT_CURVE:'
      write ( *, '(a)' ) '  FORTRAN77 version'
      write ( *, '(a)' ) '  Demonstrate the DISLIN "quickplot" command'
      write ( *, '(a)' ) '  QPLOT to plot a curve.'
c
c  Set up the X and Y data for the plot.
c
      do i = 1, n
        xray(i) = real (  i - 1 ) * 360.0 / real ( n - 1 )
      end do

      do i = 1, n
        yray(i) = sin ( pi * xray(i) / 180.0 )
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
      call setfil ( 'quickplot_curve.png' )
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
      call name ( '<-- Angle in Degrees -->', 'X' )
      call name ( '<-- Sine (angle) -->', 'Y' )
      call titlin ( 'Quick plot by QPLOT', 2 )
c
c  Draw the curve.
c
      call qplot ( xray, yray, n )
c
c  Close DISLIN.
c
      call disfin ( )
c
c  Terminate.
c
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'DISLIN_CURVE:'
      write ( *, '(a)' ) '  Normal end of execution.'

      stop
      end
