      program main

c*********************************************************************72
c
cc QUICKPLOT_PIE demonstrates the DISLIN quickplot command QPLOT.
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
      parameter ( n = 5 )

      integer i
      real xray(n)

      save xray

      data xray / 10.0, 20.0, 15.0, 5.0, 50.0 /

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'QUICKPLOT_PIE:'
      write ( *, '(a)' ) '  FORTRAN77 version'
      write ( *, '(a)' ) '  Demonstrate the DISLIN "quickplot" command'
      write ( *, '(a)' ) '  QPLPIE to plot a pie chart.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 
     &  '  Here, we plot 10 percent luck, 20 percent skill,'
      write ( *, '(a)' ) 
     &  '  15 percent concentrated power of will, '
      write ( *, '(a)' ) '  5 percent pleasure,'
      write ( *, '(a)' ) '  50 percent pain.'
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
      call setfil ( 'quickplot_pie.png' )
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
      call titlin ( 'Quick plot by QPLPIE', 2 )
c
c  Draw the curve.
c
      call qplpie ( xray, n )
c
c  Close DISLIN.
c
      call disfin ( )
c
c  Terminate.
c
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'DISLIN_PIE:'
      write ( *, '(a)' ) '  Normal end of execution.'

      stop
      end
