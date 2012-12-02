      program main

c*********************************************************************72
c
cc MAIN is the main program for TOMS626_PRB1.
c
c  Discussion:
c
c    This test demonstrates the TRICP routines with MODE = 0 and 2.
c
c  Modified:
c
c    11 January 2007
c
c  Author:
c
c    Albrecht Preusser
c
c  Reference:
c
c    Albrecht Preusser,
c    Algorithm 626:
c    TRICP: A Contour Plot Program for Triangular Meshes,
c    ACM Transactions on Mathematical Software,
c    Volume 10, Number 4, December 1984, pages 473-475.
c
c  Local Parameters:
c
c    Local, integer ND, the number of data points
c
c    Local, integer NT, the number of triangles.
c
      implicit none

      integer nc1
      parameter ( nc1 = 8 )
      integer nc2
      parameter ( nc2 = 20 )
      integer nd
      parameter ( nd = 30 )

      real angle
      real c1(nc1)
      real c2(nc2)
      real cmsc
      real height
      integer i
      integer iwk(31*nd)
      integer j
      integer mode
      integer n
      integer nh
      integer nt
      real scale
      real value
      real wk(5*nd)
      real x(3)
      real xd(nd)
      real xmax
      real xmin
      real xpage
      real y(3)
      real yd(nd)
      real ymax
      real ymin
      real ypage
      real zd(nd)
c
c  Example 1 is taken from algorithm 526 by Akima,
c  figure 1E of the corresponding paper.
c  The contour line with value c=15 is plotted in addition.
c
      data ( xd(i), i = 1, nd ) /
     &  11.16E+00, 24.20E+00, 19.85E+00, 10.35E+00, 19.72E+00,
     &   0.00E+00, 20.87E+00, 19.99E+00, 10.28E+00,  4.51E+00,
     &   0.00E+00, 16.70E+00,  6.08E+00, 25.00E+00, 14.90E+00, 
     &   0.00E+00,  9.66E+00,  5.22E+00, 11.77E+00, 15.10E+00,
     &  25.00E+00, 25.00E+00, 14.59E+00, 15.20E+00,  5.23E+00, 
     &   2.14E+00,  0.51E+00, 25.00E+00, 21.67E+00,  3.31E+00 /

      data ( yd(i), i = 1, nd ) /
     &   1.24E+00, 16.23E+00, 10.72E+00,  4.11E+00,  1.39E+00,
     &  20.00E+00, 20.00E+00,  4.62E+00, 15.16E+00, 20.00E+00,
     &   4.48E+00, 19.65E+00,  4.58E+00, 11.87E+00,  3.12E+00, 
     &   0.00E+00, 20.00E+00, 14.66E+00, 10.47E+00, 17.19E+00, 
     &   3.87E+00,  0.00E+00,  8.71E+00,  0.00E+00, 10.72E+00,
     &  15.03E+00,  8.37E+00, 20.00E+00, 14.36E+00,  0.13E+00 /

      data ( zd(i), i = 1, nd ) /
     &  22.15E+00,  2.83E+00,  7.97E+00, 22.33E+00, 16.83E+00, 
     &  34.60E+00,  5.74E+00, 14.72E+00, 21.59E+00, 15.61E+00,
     &  61.77E+00,  6.31E+00, 35.74E+00,  4.40E+00, 21.70E+00, 
     &  58.20E+00,  4.73E+00, 40.36E+00, 13.62E+00, 12.57E+00,
     &   8.74E+00, 12.00E+00, 14.81E+00, 21.60E+00, 26.50E+00,
     &  53.10E+00, 49.43E+00,  0.60E+00,  5.52E+00, 44.08E+00 /

      data ( c1(i), i = 1, nc1 ) / 
     &   5.0E+00, 10.0E+00, 15.0E+00, 20.0E+00, 30.0E+00,
     &  40.0E+00, 50.0E+00, 60.0E+00 /

      write ( *, '(a)' ) ' '
      call timestamp ( )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TOMS626_PRB1:'
      write ( *, '(a)' ) '  FORTRAN77 version'
      write ( *, '(a)' ) '  An example calling program for TOMS626.'
c
c  Set CMSC to 2.54 when using an inch-calibrated plotter
c  for adjusting the height of symbols.
c
      cmsc = 2.54E+00
c
c  Determine the limits of the physical data to be plotted.
c
      xmin = xd(1)
      xmax = xd(1)
      ymin = yd(1)
      ymax = yd(1)

      do i = 1, nd
        xmin = min ( xmin, xd(i) )
        xmax = max ( xmax, xd(i) )
        ymin = min ( ymin, yd(i) )
        ymax = max ( ymax, yd(i) )
      end do

      write ( *, '(a)' ) ' '
      write ( *, '(g14.6,a,g14.6)' ) xmin, ' <= X <= ', xmax
      write ( *, '(g14.6,a,g14.6)' ) ymin, ' <= Y <= ', ymax
c
c  Rescale the data so that 
c    0.5 <= X <= 8
c    0.5 <= Y <= 10.5
c
      scale = min (  7.5E+00 / ( xmax - xmin ), 
     &              10.0E+00 / ( ymax - ymin ) )

      do i = 1, nd
        xd(i) = 0.5E+00 + scale * ( xd(i) - xmin )
        yd(i) = 0.5E+00 + scale * ( yd(i) - ymin )
      end do
c
c  Initialize plotting.
c
      call plots ( 0, 0, 1 )
c
c  Example 1.
c  Plot contour lines (MODE=0)
c
      mode = 0
      call tricp ( xd, yd, zd, nd, c1, nc1, wk, iwk, mode )
c
c  Mark the data points.
c
      do i = 1, nd

        xpage = xd(i)
        ypage = yd(i)
        height = 0.35E+00 / cmsc
        angle = 0.0E+00

        call symbol ( xpage, ypage, height, 1, angle, -1 )

      end do
c
c  Close the plot.
c
      xpage = 0.0E+00
      ypage = 0.0E+00

      call plot ( xpage, ypage, 999 )
c
c  Example 2.
c
c  Same triangles, but different ZD values (MODE=2)
c
      do i = 1, nd
        zd(i) = real ( 5 + mod ( i, 2 ) * 10 )
      end do
c
c  Define new C values.
c
      do i = 1, nc2
        c2(i) = real ( i )
      end do

      call plots ( 0, 0, 1 )
c
c  Plot contour lines (MODE=2)
c
      mode = 2
      call tricp ( xd, yd, zd, nd, c2, nc2, wk, iwk, mode )
c
c  Plot the triangles.
c
      nt = iwk(1)

      do i = 1, nt
c
c  Load the vertices.
c
        do j = 1, 3
          n = 3 * i + 2 - j
          nh = iwk(n)
          x(j) = xd(nh)
          y(j) = yd(nh)
        end do
c
c  Draw the triangle.
c
        call plot ( x(1), y(1), 3 )
        call plot ( x(2), y(2), 2 )
        call plot ( x(3), y(3), 2 )
        call plot ( x(1), y(1), 2 )
c
c  Label the triangle with its index.
c
        xpage = ( x(1) + x(2) + x(3) ) / 3.0E+00 - 0.3E+00 / cmsc
        ypage = ( y(1) + y(2) + y(3) ) / 3.0E+00 - 0.15E+00 / cmsc
        height = 0.35E+00 / cmsc
        value = real ( i )
        angle = 0.0E+00

        call number ( xpage, ypage, height, value, angle, -1 )

      end do
c
c  Annotate data points.
c
      do i = 1, nd

        xpage = xd(i) - 0.3E+00 / cmsc
        ypage = yd(i) + 0.1E+00 / cmsc
        height = 0.35E+00 / cmsc
        value = real ( i )
        angle = 0.0E+00

        call number ( xpage, ypage, height, value, angle, -1 )

        xpage = xd(i) + 0.1E+00 / cmsc
        ypage = yd(i) - 0.1E+00 / cmsc
        height = 0.35E+00 / cmsc
        value = zd(i)
        angle = 0.0E+00

        call number ( xpage, ypage, height, value, angle, 1 )

      end do
c
c  Close the plot.
c
      xpage = 0.0E+00
      ypage = 0.0E+00

      call plot ( xpage, ypage, 999 )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TOMS626_PRB1:'
      write ( *, '(a)' ) '  Normal end of execution.'

      write ( *, '(a)' ) ' '
      call timestamp ( )

      stop
      end
