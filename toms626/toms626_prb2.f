      program main

c*********************************************************************72
c
cc MAIN is the main program for TOMS626_PRB2.F
c
c  Discussion:
c
c    This test demonstrates the TRICP routines with MODE = 0 and 2.
c
c  Modified:
c
c    10 January 2007
c
c  Author:
c
c    Albrecht Preusser
C
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
c    Local, integer ND, the number of data points.
c
c    Local, integer NT, the number of triangles.
c
      implicit none

      integer nc
      parameter ( nc = 8 )
      integer nd
      parameter ( nd = 47 )
      integer nt_copy
      parameter ( nt_copy = 53)

      real ai
      real angle
      real c(nc)
      real height
      integer i
      integer iwk(3*nt_copy+4*nd+1)
      integer j
      integer jp(11)
      integer jpp
      integer mode
      integer n
      integer nh
      integer nt
      real pd(5,47)
      real scale
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

      data ( xd(i), i = 1, nd ) /
     &  21.95E+00, 21.95E+00, 19.70E+00, 19.70E+00, 17.50E+00,
     &  17.50E+00, 15.30E+00, 15.30E+00, 13.00E+00, 13.00E+00,
     &  10.80E+00, 10.80E+00,  8.60E+00,  8.60E+00,  7.04E+00,
     &   6.77E+00,  6.48E+00,  6.01E+00,  5.44E+00,  5.59E+00,
     &   5.37E+00,  5.51E+00,  4.84E+00,  4.73E+00,  4.77E+00,
     &   4.32E+00,  3.55E+00,  3.56E+00,  3.80E+00,  3.49E+00,
     &   3.09E+00,  3.06E+00,  2.57E+00,  2.58E+00,  2.42E+00,
     &   2.44E+00,  2.64E+00,  2.03E+00,  1.24E+00,  1.43E+00,
     &   1.96E+00,  1.36E+00,  0.53E+00,  1.28E+00,  0.68E+00, 
     &   0.90E+00,  0.29E+00 /

      data ( yd(i), i = 1, nd ) /
     &   5.10E+00, 7.25E+00, 5.10E+00, 7.25E+00, 5.10E+00,
     &   7.25E+00, 5.10E+00, 7.25E+00, 5.10E+00, 7.25E+00,
     &   5.10E+00, 7.25E+00, 5.10E+00, 7.25E+00, 5.10E+00,
     &   6.17E+00, 7.29E+00, 5.10E+00, 5.73E+00, 6.55E+00,
     &   7.42E+00, 4.82E+00, 6.66E+00, 7.87E+00, 4.10E+00,
     &   5.45E+00, 7.05E+00, 8.80E+00, 3.45E+00, 4.46E+00,
     &   5.03E+00, 6.15E+00, 7.06E+00, 7.84E+00, 8.83E+00,
     &   3.52E+00, 4.76E+00, 7.13E+00, 8.22E+00, 3.83E+00,
     &   4.58E+00, 6.96E+00, 7.37E+00, 4.40E+00, 6.78E+00,
     &   4.30E+00, 6.68E+00 /

      data ( zd(i), i = 1, nd ) /
     &   0.1E+00,  0.1E+00,  0.1E+00,  0.1E+00,  1.0E+00,
     &   1.0E+00,  2.1E+00,  2.1E+00,  3.6E+00,  3.6E+00,
     &   5.0E+00,  5.0E+00,  6.8E+00,  7.0E+00,  8.5E+00,
     &   0.2E+00,  8.2E+00,  8.5E+00,  1.0E+00,  2.7E+00,
     &   6.0E+00,  8.0E+00,  3.0E+00,  1.8E+00,  4.5E+00,
     &   0.9E+00,  2.2E+00,  0.6E+00,  7.0E+00,  4.5E+00,
     &   3.0E+00,  2.2E+00,  9.8E+00,  1.8E+00,  1.9E+00,
     &   7.0E+00,  4.8E+00,  9.5E+00,  1.0E+00,  2.2E+00,
     &   8.2E+00,  1.8E+00,  0.5E+00, 15.0E+00,  1.0E+00,
     &   4.0E+00,  0.0E+00 /
c
c     iwk(2...3*nt+1) point numbers for triangles
c
      data ( iwk(i), i = 1, 160 ) /
     &  53,
     &   1,  2,  3,
     &   3,  2,  4,
     &   5,  3,  4,
     &   5,  4,  6,
     &   7,  5,  6,
     &   7,  6,  8,
     &   9,  7,  8,
     &   9,  8, 10,
     &  11,  9, 10,
     &  11, 10, 12,
     &  13, 11, 12,
     &  13, 12, 14,
     &  13, 14, 16,
     &  15, 13, 16,
     &  16, 14, 17,
     &  18, 15, 16,
     &  18, 16, 19,
     &  19, 16, 20,
     &  20, 16, 17,
     &  20, 17, 21,
     &  22, 18, 19,
     &  26, 22, 19,
     &  26, 19, 23,
     &  23, 19, 20,
     &  23, 20, 21,
     &  23, 21, 24,
     &  26, 23, 27,
     &  32, 26, 27,
     &  27, 23, 24,
     &  33, 32, 27,
     &  33, 27, 34,
     &  34, 27, 28,
     &  27, 24, 28,
     &  38, 33, 34,
     &  34, 28, 35,
     &  38, 34, 39,
     &  39, 34, 35,
     &  42, 38, 39,
     &  43, 42, 39,
     &  43, 45, 42,
     &  47, 45, 43,
     &  46, 40, 44,
     &  40, 41, 44,
     &  40, 36, 41,
     &  36, 37, 41,
     &  36, 30, 37,
     &  36, 29, 30,
     &  30, 31, 37,
     &  29, 25, 30,
     &  25, 22, 26,
     &  25, 26, 30,
     &  30, 26, 31,
     &  31, 26, 32 /
c
c  Point numbers for which derivatives will be set to zero.
c
      data jp /
     &   7,  8,  9, 10, 11,
     &  12, 13, 14, 16, 19, 
     &  26 /
c
c  Contour levels
c
      data ( c(i), i = 1, nc ) / 
     &   0.0E+00,  2.0E+00,  4.0E+00,  6.0E+00,  8.0E+00, 
     &  10.0E+00, 12.0E+00, 14.0E+00 /

      write ( *, '(a)' ) ' '
      call timestamp ( )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TOMS626_PRB2:'
      write ( *, '(a)' ) '  FORTRAN77 version'
      write ( *, '(a)' ) '  An example calling program for TOMS626.'
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
c  Example 3.
c
c  Plot contour lines (MODE = 1)
c
      mode = 1
      call tricp ( xd, yd, zd, nd, c, nc, pd, iwk, mode )
c
c  Plot the triangle sides.
c
      nt = iwk(1)

      do i = 1, nt

        do j = 1, 3
          n = 3 * i + 2 - j
          nh = iwk(n)
          x(j) = xd(nh)
          y(j) = yd(nh)
        end do

        call plot ( x(1), y(1), 3 )
        call plot ( x(2), y(2), 2 )
        call plot ( x(3), y(3), 2 )
        call plot ( x(1), y(1), 2 )

      end do
c
c  Plot point numbers
c
      do i = 1, nd

        xpage = xd(i) - 0.3E+00
        ypage = yd(i) + 0.1E+00
        height = 0.2E+00
        ai = real ( i )
        angle = 0.0E+00

        call number ( xpage, ypage, height, ai, angle, -1 )

      end do
c
c  Close the plot.
c
      xpage = 0.0E+00
      ypage = 0.0E+00

      call plot ( xpage, ypage, 999 )
c
c  Change partial derivatives (slopes and curvatures)
c  at some data points and redraw the plot with mode= 3
c
c  Set the slopes and curvatures at the jp-points to zero
c
      do j = 1, 11
        jpp = jp(j)
        do i = 1, 5
          pd(i,jpp) = 0.0E+00
        end do
      end do
c
c  Define X and Y slopes for points 8, 10, 12, 14.
c
      do j = 8, 14, 2
        pd(1,j) = -1.0E+00
        pd(2,j) = 8.5+00
      end do
c
c  Define X and Y slopes for points 7, 9, 11, 13.
c
      do j = 7, 13, 2
        pd(1,j) = -1.0E+00
        pd(2,j) = -8.5E+00
      end do
c
c  Define the Y slope at point 26.
c
      pd(2,26) = 2.0E+00
c
c  Example 4.
c
c  Initialize plotting.
c
      call plots ( 0, 0, 1 )
c
c  Plot the contour lines (MODE = 3)
c
      mode = 3
      call tricp ( xd, yd, zd, nd, c, nc, pd, iwk, mode )
c
c  Plot the triangle.
c
      do i = 1, nt

        do j = 1, 3
          n = 3 * i + 2 - j
          nh = iwk(n)
          x(j) = xd(nh)
          y(j) = yd(nh)
        end do

        call plot ( x(1), y(1), 3 )
        call plot ( x(2), y(2), 2 )
        call plot ( x(3), y(3), 2 )
        call plot ( x(1), y(1), 2 )

      end do
c
c  Terminate the plot.
c
      xpage = 0.0E+00
      ypage = 0.0E+00

      call plot ( xpage, ypage, 999 )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TOMS626_PRB2:'
      write ( *, '(a)' ) '  Normal end of execution.'

      write ( *, '(a)' ) ' '
      call timestamp ( )

      stop
      end
