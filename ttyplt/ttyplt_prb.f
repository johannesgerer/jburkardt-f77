      program main

c*********************************************************************72
c
cc MAIN is the main program for TTYPLT_PRB.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    21 May 2008
c
c  Author:
c
c    John Burkardt
c
      implicit none
c
c  You need to declare the common blocks here, or the data they hold
c  may not be preserved.
c
      character*1 ibak
      integer igrid
      character*1 ihoriz
      character*1 islant
      character*1 istar
      character*1 iverti
      character*1 lgray(30)
      integer margl
      integer ncolr
      integer ngray
      character*1 node
      integer nrowr
      real zmaxw
      real zminw

      common /chdata/ ibak,ihoriz,islant,iverti,node
      common /indata/ igrid
      common /gracom/ margl,ncolr,ngray,nrowr,zmaxw,zminw
      common /chrgra/ lgray

      call timestamp ( )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TTYPLT_PRB'
      write ( *, '(a)' ) '  FORTRAN77 version.'
      write ( *, '(a)' ) '  Test the TTYPLT library.'

      call test01 ( )
      call test02 ( )
      call test03 ( )
      call test04 ( )
      call test05 ( )
      call test06 ( )
      call test07 ( )
      call test08 ( )
c
c  Terminate
c
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TTYPLT_PRB'
      write ( *, '(a)' ) '  Normal end of execution.'
      write ( *, '(a)' ) ' '
      call timestamp ( )

      stop
      end
      subroutine test01 ( )

c*********************************************************************72
c
cc TEST01 makes a dot plot of a spiral.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    21 May 2008
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer nval
      parameter ( nval = 52 )

      integer i
      character*1 ibak
      character*1 ichars(nval)
      integer igrid
      integer ixmax
      integer iymax
      integer letter
      integer margel
      real t
      real xvalue(nval)
      real yvalue(nval)

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST01'
      write ( *, '(a)' ) '  DOTPLT is used to create a dot plot'
      write ( *, '(a)' ) '  of a spiral.'
c
c  Set defaults.
c
      call setdef
c
c  Set left margin to 5.
c
      margel = 5
      call setmar ( margel )
c
c  Set the grid option to no labeling at all.
c
      igrid = 0
      call setigr ( igrid )
c
c  Set the background character to a period.
c
      ibak = '.'
      call setbak ( ibak )
c
c  Set up the data for the spiral.
c
c  Here, the example is a little complicated because
c  we are going to specify the particular letter to print
c  at each point.  We could let the program choose one
c  fixed character to use, which would be simpler.
c
      ixmax = 41
      iymax = 20

      do i = 1, nval
        t = real ( i ) / 4.0E+00
        yvalue(i) = t * sin ( t )
        xvalue(i) = t * cos ( t )
        letter = i - 1 + ichar ( 'A' )
        if ( i .ge. 27) letter = i - 27 + ichar ( 'A' )
        ichars(i) = char ( letter )
      end do

      call dotplt ( ichars, ixmax, iymax, nval, nval, xvalue, yvalue )

      return
      end
      subroutine test02 ( )

c*********************************************************************72
c
cc TEST02 is a dot plot of a sine curve.
c
c  Discussion:
c
C    Make background be 'S' and dots themselves blanks.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    21 May 2008
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer nval
      parameter ( nval = 101 )

      integer i
      character*1 ibak
      character*1 ichars(1)
      integer igrid
      integer ixmax
      integer iymax
      integer nchar
      character*1 node
      real pi
      parameter ( pi = 3.14159265 )
      real t
      real xvalue(nval)
      real yvalue(nval)

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST02'
      write ( *, '(a)' ) '  Use DOTPLT to draw a sine curve.'
      write ( *, '(a)' ) '  We can easily draw it in "reverse video".'
c
c  Restore defaults.
c
      call setdef
c
c  Set grid option.
c
      igrid = 1
      call setigr ( igrid )
c
c  Use letter 'S' for background.
c
      ibak = 'S'
      call setbak ( ibak )
c
c  Use blank for dots.
c
      node = ' '
      call setnod ( node )
c
c  Set up sine curve.
c
      ixmax = 71
      iymax = 20
      do i = 1, nval
        t = 0.0 + real ( i - 1 ) * 2.0 * pi / real ( nval - 1 )
        xvalue(i) = t
        yvalue(i) = sin(t)
      end do
      nchar = 0
 
      call dotplt ( ichars, ixmax, iymax, nchar, nval, xvalue, yvalue )

      return
      end
      subroutine test03 ( )

c*********************************************************************72
c
cc TEST03 plots a doghouse.
c
c  Discussion:
c
c    Plot of 'doghouse' using no grid lines, 'D' for nodes,
c    '-' for horizontal lines,
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    21 May 2008
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer nval
      parameter ( nval = 8 )

      integer igrid
      integer ixmax
      integer iymax
      character*1 node
      real x1(nval)
      real x2(nval)
      real y1(nval)
      real y2(nval)

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST03'
      write ( *, '(a)' ) '  Use DOTPLT to draw a doghouse.'
c
c  Restore defaults.
c
      call setdef
c
c  Cancel grid lines.
c
      igrid = 0
      call setigr ( igrid )
c
c  Set nodes to 'D'.
c
      node = 'D'
      call setnod ( node )
c
c  Set up data.
c
      x1(1) = 0.75
      y1(1) = 3.0
      x2(1) = 1.0
      y2(1) = 3.0
      x1(2) = 1.0
      y1(2) = 3.0
      x2(2) = 1.0
      y2(2) = 2.0
      x1(3) = 1.0
      y1(3) = 2.0
      x2(3) = 0.0
      y2(3) = 1.0
      x1(4) = 0.0
      y1(4) = 1.0
      x2(4) = 2.0
      y2(4) = 1.0
      x1(5) = 1.0
      y1(5) = 2.0
      x2(5) = 2.0
      y2(5) = 1.0
      x1(6) = 0.0
      y1(6) = 1.0
      x2(6) = 0.0
      y2(6) = -1.0
      x1(7) = 2.0
      y1(7) = 1.0
      x2(7) = 2.0
      y2(7) = -1.0
      x1(8) = 0.0
      y1(8) = -1.0
      x2(8) = 2.0
      y2(8) = -1.0

      iymax = 20
      ixmax = 61

      call grdplt ( ixmax, iymax, nval, x1, x2, y1, y2 )

      return
      end
      subroutine test04 ( )

c*********************************************************************72
c
cc TEST04 draws a cube with grid points, and different letters for lines.l
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    21 May 2008
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer nval
      parameter ( nval = 12 )

      character*1 arrow
      character*1 buck
      integer igrid
      character*1 itic
      integer ixmax
      integer iymax
      character*1 wow
      real x1(nval)
      real x2(nval)
      real y1(nval)
      real y2(nval)

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST04'
      write ( *, '(a)' ) '  Use GRDPLT to draw a cube with grid points.'
c
c  Restore defaults.
c
      call setdef
c
c  Set grid option.
c
      igrid = 1
      call setigr ( igrid )
c
c  Reset letters used to draw lines.
c
      arrow = '>'
      buck = '$'
      wow = '!'

      call setlet ( arrow, buck, wow )
c
c  Reset letter used to draw nodes.
c
      itic = '#'
      call setnod ( itic )
c
c  Store data describing beginning and end points of lines to be drawn.
c
      x1(1) = 0.0
      y1(1) = 0.0
      x2(1) = 2.0
      y2(1) = 0.0
      x1(2) = 0.0
      y1(2) = 2.0
      x2(2) = 2.0
      y2(2) = 2.0
      x1(3) = 1.0
      y1(3) = 1.0
      x2(3) = 3.0
      y2(3) = 1.0
      x1(4) = 1.0
      y1(4) = 3.0
      x2(4) = 3.0
      y2(4) = 3.0
      x1(5) = 0.0
      y1(5) = 0.0
      x2(5) = 0.0
      y2(5) = 2.0
      x1(6) = 2.0
      y1(6) = 0.0
      x2(6) = 2.0
      y2(6) = 2.0
      x1(7) = 1.0
      y1(7) = 1.0
      x2(7) = 1.0
      y2(7) = 3.0
      x1(8) = 3.0
      y1(8) = 1.0
      x2(8) = 3.0
      y2(8) = 3.0
      x1(9) = 0.0
      y1(9) = 0.0
      x2(9) = 1.0
      y2(9) = 1.0
      x1(10) = 2.0
      y1(10) = 0.0
      x2(10) = 3.0
      y2(10) = 1.0
      x1(11) = 0.0
      y1(11) = 2.0
      x2(11) = 1.0
      y2(11) = 3.0
      x1(12) = 2.0
      y1(12) = 2.0
      x2(12) = 3.0
      y2(12) = 3.0

      call grdplt ( ixmax, iymax, nval, x1, x2, y1, y2 )

      return
      end
      subroutine test05 ( )

c*********************************************************************72
c
cc TEST05 draws a star.
c
c  Discussion:
c
c    Note that since the star can be drawn
c    with a single series of connected line segments,
c    we can use just X1, Y1 to store both start and end
c    points of lines, by passing NVAL+1 values in X1, Y1,
c    and position X1(2), Y1(2) as the vectors X2 and Y2
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    21 May 2008
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer nval
      parameter ( nval = 5 )

      real delth
      integer i
      character*1 iblank
      integer igrid
      integer ihi
      character*1 istar
      integer ixmax
      integer iymax
      real pi
      real theta
      real x1(nval+1)
      real y1(nval+1)

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST05'
      write ( *, '(a)' ) '  Use GRDPLT to draw a star.'
c
c  Restore defaults.
c
      call setdef
c
c  No grid lines, please.
c
      igrid = 0
      call setigr ( igrid )
c
c  Set line letters to star.
c
      istar = '*'
      call setlet ( istar, istar, istar )
c
c  Blank out node dots.
c
      iblank = ' '
      call setnod(iblank)
c
c  Set up data for star.
c
      pi = 3.14159265
      theta = pi / 10.0
      delth = 8.0 * pi / 10.0
      ihi = nval + 1
      do i = 1, ihi
        x1(i) = cos ( theta )
        y1(i) = sin ( theta )
        theta = theta + delth
      end do
c
c  x1(1) through x1(5) are beginning points,
c  x1(2) through x1(6) are end points.
c  similarly for y1.
c
      ixmax = 41
      iymax = 20

      call grdplt ( ixmax, iymax, nval, x1(1), x1(2), y1(1), y1(2) )

      return
      end
      subroutine test06 ( )

c*********************************************************************72
c
cc TEST06 draws a gray plot for Z(X,Y) defined by a subroutine.
c
c  Discussion:
C
C    Draw gray plot using a subroutine zfunc,
c    setting maximum and minimum values, declaring
c    height and width of region, and calling grafun once.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    21 May 2008
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer igray
      character*1 lgray(1)
      integer margel
      integer ncolr
      integer ngray
      integer nrowr
      real winmax
      real winmin
      real xmax
      real xmin
      real ymax
      real ymin
      external zfunc

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST06'
      write ( *, '(a)' ) '  Use GRAFUN to draw a gray scale plot of'
      write ( *, '(a)' ) '  the contours of a cubic Z(X,Y), whose'
      write ( *, '(a)' ) '  values are defined by a subroutine.'
c
c  Restore defaults.
c
      call setdef
c
c  Define size of region.
c
      nrowr = 20
      ncolr = 41
      call setreg ( ncolr, nrowr )
c
c  Set left margin.
c
      margel = 5
      call setmar ( margel )
c
c  Request first gray scale.
c
      igray = 1
      call setgra ( igray, lgray, ngray ) 
c
c  Set maximum and minimum window values.
c
      winmax = 0.7
      winmin = -0.7
      call setwin ( winmax, winmin )
c
c  Set x and y boundaries.
c
      xmin = -2.0
      xmax = 2.0
      ymin = -2.0
      ymax = 2.0
c
c  Call gray plotter.
c
      call grafun ( xmax, xmin, ymax, ymin, zfunc )

      return
      end
      subroutine test07 ( )

c*********************************************************************72
c
cc TEST07 draws a gray plot for Z(X,Y) by storing values in an array.
c
c  Discussion:
C
C    Create gray plot by storing entries in an array zarray,
c    declare number of rows and columns in array zarray,
c    choose maximum and minimum window values,
c    declaring height and width of region, and calling graray once.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    21 May 2008
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer ncolz
      parameter ( ncolz = 41 )
      integer ngray
      parameter ( ngray = 12 )
      integer nrowz
      parameter ( nrowz = 20 )

      integer i
      integer igray
      integer j
      character*1 lgray(ngray)
      integer margel
      integer ncolr
      integer nrowr
      real xval
      real yval
      real zarray(nrowz,ncolz)
      real zmax
      real zmin

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST07'
      write ( *, '(a)' ) '  Use GRARAY to draw a gray scale plot'
      write ( *, '(a)' ) '  that looks like water waves.'
c
c  Restore defaults.
c
      call setdef
c
c  Define size of region.
c
      nrowr = 20
      ncolr = 41
      call setreg ( ncolr, nrowr )
c
c  Set left margin.
c
      margel = 0
      call setmar ( margel )
c
c  Set up a special gray scale.
c
      igray = 0

      lgray(1) = ' '
      lgray(2) = '1'
      lgray(3) = ' '
      lgray(4) = '3'
      lgray(5) = ' '
      lgray(6) = '5'
      lgray(7) = ' '
      lgray(8) = '7'
      lgray(9) = ' '
      lgray(10) = '9'
      lgray(11) = ' '
      lgray(12) = '*'
 
      call setgra ( igray, lgray, ngray ) 
c
c  Set values in zarray.
c
      do i = 1, nrowr
        yval = 10.0 + 2.0 * real ( i - 1 ) * ( -10.0 ) / ( nrowr - 1 )
        do j = 1, ncolr
          xval = - 10.0 
     &      + 2.0 * real ( j - 1 ) * ( 10.0 ) / ( ncolr - 1 )
          zarray(i,j) = xval * xval + yval * yval
        end do
      end do
c
c  Set windows on Z values.
c
      zmax = 100.0
      zmin = 0.0
      call setwin ( zmax, zmin )
c
c  Call gray plotter.
c
      call graray ( ncolz, nrowz, zarray )

      return
      end
      subroutine test08 ( )

c*********************************************************************72
c
cc TEST08 draws a gray plot of Z(X,Y) by supplying one vector of data at a time.
c
c  Discussion:
c
c    Create gray plot by storing entries, one row per call,
c    in a vector ZVEC.
c
c    declare maximum and minimum, since they can not be found,
c    declare height and width of region,
c    and call gravec.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    21 May 2008
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer nvec
      parameter ( nvec = 80 )

      integer i
      integer igray
      integer j
      character*1 lgray(1)
      integer margel
      integer ncolr
      integer ngray
      integer nrowr
      real winmax
      real winmin
      real x
      real y
      real zvec(nvec)

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST08'
      write ( *, '(a)' ) '  Use GRAVEC to draw a gray scale plot'
      write ( *, '(a)' ) '  that looks like alphabet soup.'
c
c  Restore defaults.
c
      call setdef
c
c  Set maximum and minimum of range.
c
      winmax = 400.0
      winmin = -300.0
      call setwin ( winmax, winmin )
c
c  Define size of region.
c
      nrowr = 20
      ncolr = 2 * nrowr - 1
      call setreg ( ncolr, nrowr )
c
c  Set left margin.
c
      margel = 10
      call setmar ( margel )
c
c  Request third gray scale.
c
      igray = 3
      call setgra ( igray, lgray, ngray )
c
c  Send Z data, one row at a time, to GRAVEC.
c
      do i = 1, nrowr
        y = real ( 2 * i - 26 )
        do j = 1, ncolr
          x = real ( j - 25 )
          zvec(j) =  x * x + y * y - 400.0
        end do
        call gravec ( nvec, zvec )
      end do
c
c  Print out the gray scale.
c
      call getscl

      return
      end
      subroutine zfunc ( xval, yval, zval )

c*********************************************************************72
c
cc ZFUNC defines a function whose gray plot is desired.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    21 May 2008
c
c  Author:
c
c    John Burkardt
c
      implicit none

      real xval
      real yval
      real zval

      zval = xval * ( xval - 1.0 ) * ( xval + 1.0 ) - yval

      return
      end
