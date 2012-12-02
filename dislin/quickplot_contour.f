      program main

c*********************************************************************72
c
cc QUICKPLOT_CONTOUR demonstrates the DISLIN quickplot command QPLCON.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    14 May 2012
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

      integer m
      parameter ( m = 100 )
      integer n
      parameter ( n = 100 )

      integer i
      integer j
      integer levels
      real x
      real y
      real zmat(m,n)

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'QUICKPLOT_CONTOUR:'
      write ( *, '(a)' ) '  FORTRAN77 version'
      write ( *, '(a)' ) '  Demonstrate the DISLIN "quickplot" command'
      write ( *, '(a)' ) '  QPLCON to make a contour plot of data'
      write ( *, '(a)' ) '  stored as a matrix.'
c
c  Set up the data.
c
      do i = 1, m
        x = 1.6 * real ( i - 1 ) / real ( m - 1 )
        do j = 1, n
          y = 1.6 * real ( j - 1 ) / real ( n - 1 )
          zmat(i,j) = ( x * x - 1.0 )**2 + ( y * y - 1.0 )**2
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
      call setfil ( 'quickplot_contour.png' )
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
      call name ( '<-- X -->', 'X' )
      call name ( '<-- Y -->', 'Y' )
      call titlin ( 'Quick plot by QPLCON', 2 )
c
c  Draw the curve.
c
      levels = 20
      call qplcon ( zmat, m, n, levels )
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
