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
c    30 April 2011
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

      integern
      parameter ( n = 14 )

      integer i
      real xray(n)

      save xray

      data xray /
     &   1.0, 15.0, 38.0, 22.0, 16.0, 
     &  16.0, 26.0, 55.0, 50.0, 40.0, 
     &  16.0,  3.0,  0.0,  1.0 /
        
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'QUICKPLOT_BAR:'
      write ( *, '(a)' ) '  FORTRAN77 version'
      write ( *, '(a)' ) '  Demonstrate the DISLIN "quickplot" command'
      write ( *, '(a)' ) '  QPLBAR to plot a bar chart.'
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
      call setfil ( 'quickplot_bar.png' )
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
      call name ( '<-- Minutes -->', 'X' )
      call name ( '<-- Frequency -->', 'Y' )
      call titlin ( 'Quick plot by QPLBAR', 2 )
c
c  Draw the curve.
c
      call qplbar ( xray, n )
c
c  Close DISLIN.
c
      call disfin ( )
c
c  Terminate.
c
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'QUICKPLOT_BAR:'
      write ( *, '(a)' ) '  Normal end of execution.'

      stop
      end
