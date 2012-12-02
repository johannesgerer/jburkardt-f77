      subroutine do_back ( ixmin, ixmax, icolum )

c*********************************************************************72
c
cc DO_BACK fills in the background.
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
c  Parameters:
c
c    Input, integer IXMIN, IXMAX, the maximum and minimum indices in 
c    COLUMN used for the plot.
c
c    Output, character*1 ICOLUM(IXMAX), the initialized graphics line.
c
      implicit none

      integer ixmax

      character*1 ibak
      character*1 icolum(ixmax)
      character*1 ihoriz
      character*1 islant
      character*1 iverti
      integer ixmin
      integer j
      character*1 node

      common /chdata/ ibak,ihoriz,islant,iverti,node

      do j = ixmin, ixmax
        icolum(j) = ibak
      end do

      return
      end
      subroutine do_grid ( igrid, icolum, im1, ixlin, ixmax, ixmin )

c*********************************************************************72
c
cc DO_GRID places grid lines and axis ticks on the graphics line.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    19 May 2008
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer IGRID, controls the drawing of the grid.
c    0 or less, no grid is drawn.
c    1 or more, grid nodes are drawn, and numeric tick marks.
c    2 or more, grid lines are drawn.
c
c    Input/output, character*1 ICOLUM(IXMAX), the graphics line.
c
c    Input, integer IXMAX, IXMIN, the maximum and minimum indices in 
c    COLUMN used for the plot.
c
      implicit none

      integer ixmax

      character*1 icolum(ixmax)
      integer igrid
      integer im1
      character*1 iminus
      character*1 iplus
      integer ixlin
      integer ixmin
      integer j

      if ( igrid .le. 0 ) then
        return
      end if
c
c  Search for grid lines
c
      do j = ixmin, ixmax

        if ( 1 .lt. igrid ) then

          if ( mod ( im1, 5 ) .eq. 0 ) then
            icolum(j) = '-'
            ixlin = max ( ixlin, j )
          end if

          if ( mod ( j - 1, 10 ) .eq. 0 ) then
            icolum(j) = ':'
            ixlin = max ( ixlin, j )
          end if

        end if

        if(mod(im1,5).eq.0.and.mod(j-1,10).eq.0) then
          icolum(j) = '+'
          ixlin = max0(ixlin,j)
        end if

      end do

      return
      end
      subroutine dotplt ( ichars, ixmax, iymax, nchar, nval, xvalue,
     &  yvalue )

c*********************************************************************72
c
cc DOTPLT makes an (X,Y) "dot plot" based on lists of X and Y values.
c
c  THE SUBROUTINES FOR USE WITH DOTPLT AND GRDPLT ARE-
c
c  SUBROUTINE SETLET(IHORIZ,ISLANT,IVERTI)
c
c  FOR SUBROUTINE GRDPLT, HORIZONTAL, SLANTED, AND VERTICAL
c  LINES ARE EACH MARKED WITH A SEPARATE character.  IT IS
c  POSSIBLE TO ALTER EACH OF THESE characterS BY CALLING
c  SETLET.  THE DEFAULT VALUES ARE
c  IHORIZ = '-', ISLANT = '*', IVERTI = '!'.
c  IHORIZ, ISLANT, AND IVERTI ARE OF TYPE character*1.
c
c
c
c  SUBROUTINES FOR USE WITH GRAFUN, GRARAY AND GRAVEC-
c
c
c  SUBROUTINE SETGRA(IGRAY,LGRAY,NGRAY)
c
c  SETGRA IS USED TO CHOOSE THE NUMBER AND TYPE OF SYMBOLS
c  TO BE USED IN THE PLOT.  NOTE THAT THE FIRST character IN LGRAY
c  IS USED FOR VALUES BELOW THE MINIMUM, AND THE LAST FOR
c  VALUES ABOVE THE MAXIMUM.  LGRAY IS OF TYPE character*1.
c
c  IF IGRAY = 0, THE USER HAS DIMENSIONED THE ARRAY LGRAY
c  TO AT LEAST THE VALUE NGRAY, AND HAS STORED IN LGRAY
c  THE NGRAY characterS TO BE USED.  NGRAY MUST BE AT LEAST 3.
c
c  IF IGRAY = 1, THE USER REQUESTS THAT THE 14 SYMBOL GRAY
c  SCALE WHICH GOES FROM LIGHT TO DARK IS TO BE USED.
c  (LGRAY(1) = LGRAY(16) = ' '.)
c  THE INPUT LGRAY AND NGRAY ARE IGNORED.
c
c  IF IGRAY = 2, THE USER REQUESTS THAT THE 10 SYMBOLS
c  0 THROUGH 9 BE USED (WITH LGRAY(1) = LGRAY(12) = ' ').
c  THE INPUT LGRAY AND NGRAY ARE IGNORED.
c
c  IF IGRAY = 3, THE 26 characterS A THROUGH Z ARE USED
c  (WITH LGRAY(1) = LGRAY(28) = ' ').
c  THE INPUT LGRAY AND NGRAY ARE IGNORED.
c
c  SUBROUTINE SETREG(NCOLR,NROWR)
c
c  THE USER MUST CALL THIS ROUTINE BEFORE CALLING THE MAIN
c  GRAY SCALE PLOTTERS.  THE USER MUST INPUT THE NUMBER OF
c  COLUMNS AND THE NUMBER OF ROWS IN THE REGION TO BE DRAWN.
c  HOWEVER, THE NUMBER OF ROWS IS NOT NEEDED FOR CALLING GRAVEC,
c  EVEN IF A REGION OF SEVERAL ROWS IS TO BE DRAWN,
c  SINCE THE USER TAKES CARE OF THIS BY CALLING GRAVEC ONCE
c  FOR EACH ROW TO BE DRAWN.
c
c  SUBROUTINE SETWIN(WINMAX,WINMIN)
c
c  THIS SUBROUTINE MUST BE CALLED BEFORE CALLING ANY OF THE
c  MAIN GRAY SCALE PLOTTERS.  IT INFORMS THE PLOTTER OF
c  THE MAXIMUM AND MINIMUM VALUES TO BE PLOTTED.  ANY VALUES
c  OUTSIDE THIS RANGE WILL BE REPRESENTED BY A SPECIAL character,
c  THOSE BELOW WINMIN BY LGRAY(1), THOSE ABOVE WINMAX BY LGRAY(NGRAY).
c  WINMAX MUST BE GREATER THAN WINMIN.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    19 May 2008
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c  ichars - a character*1 array of dimension nchar, which if nchar.gt.0,
c           contains characters to be used to mark the location
c           of the dots.
c           if nchar.le.0, ichars is ignored, and the user need
c           not even set it up.
c           if nchar.le.0, the default value supplied by getnod
c             is used.  this value can be altered by calling
c             setnod.
c           if nchar.eq.1, the single character in ichars(1) is
c             used.  this is an alternative to calling setnod.
c           if nchar.ge.nval, the characters contained in ichars
c             are assumed to correspond to the nodes
c             so that ichars(i) is the character to be printed
c             at the node (xvalue(i),yvalue(i)).
c           if nchar.gt.1 and nchar.lt.nval, only the first
c             character in ichars is used.
c
c  ixmax  - number of columns wide the graph is to be.
c           ixmax must be between 1 and 100.  note that if igrid.ne.0,
c           12 more columns are needed for labeling, so that for
c           a terminal with 72 columns, ixmax = 60 is appropriate,
c           and for a line printer, ixmax = 120 is appropriate.
c
c  iymax  - number of lines high the graph is to be.
c           this is limited only by the amount of paper you have.
c           iymax must be at least 1.
c
c  nchar  - the number of characters in ichars to be used.
c           permissible values are-
c           nchar = 0, no characters are to be used.  use the
c             default symbol supplied by getnod.
c           nchar = 1, use ichars(1) as the character to be printed
c             at each node.
c           nchar = nval, each node has its own character in ichars.
c             print ichars(i) at node (xvalue(i),yvalue(i)).
c           if nchar.gt.1 and nchar.lt.nval, the code will behave
c             as though nchar = 1.
c
c    Input, integer NVAL, the number of points to be plotted.
c    Normally, 0 < NVAL, but this is not required.
c
c  xvalue - a vector of length nval containing the x coordinates
c           of the points to be plotted.
c
c  yvalue - a vector of length nval containing the y coordinates
c           of the points to be plotted.
c
      implicit none

      integer nval

      integer i
      character*1 ibak
      character*1 iblank
      character*1 ichara
      character*1 ichars(1)
      character*1 icolum(100)
      integer igrid
      character*1 ihoriz
      integer im1
      character*1 islant
      character*10 istore
      character*1 iverti
      integer ix
      integer ixlin
      integer ixmax
      integer ixmin
      integer iy
      integer iy1
      integer iymax
      integer iymin
      integer j
      character*1 jstore(9)
      integer k
      character*1 lgray(30)
      integer margl
      integer nchar
      integer ncolr
      integer ngray
      integer nlabel
      integer nlen
      character*1 node
      character*1 node_use
      integer nrowr
      integer nstore
      character*100 output
      real r4mat_max
      real r4mat_min
      real x
      real xlabel(15)
      real xmax
      real xmin
      real xvalue(nval)
      real y
      real y1
      real ymax
      real ymin
      real yvalue(nval)
      real zmaxw
      real zminw

      common /chdata/ ibak,ihoriz,islant,iverti,node
      common /indata/ igrid
      common /gracom/ margl,ncolr,ngray,nrowr,zmaxw,zminw
      common /chrgra/ lgray

      iblank = ' '
      ichara = 'a'
      ixmin = 1
      iymin = 1
c
c  Check IXMAX, IYMAX
c
      if ( ixmax .lt. 1 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'DOTPLT - Warning!'
        write ( *, '(a)' ) '  IXMAX < 1.'
        write ( *, '(a)' ) '  This call is ignored.'
        return
      end if

      if ( 100 .lt. ixmax ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'DOTPLT - Warning!'
        write ( *, '(a)' ) '  100 < IXMAX.'
        write ( *, '(a)' ) '  This call is ignored.'
        return
      end if

      if ( iymax .lt. 1 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'DOTPLT - Warning!'
        write ( *, '(a)' ) '  IYMAX < 1.'
        write ( *, '(a)' ) '  This call is ignored.'
        return
      end if
c
c  Find ranges of X and Y.
c
      xmin = r4mat_min ( nval, 1, xvalue )
      xmax = r4mat_max ( nval, 1, xvalue )
      ymin = r4mat_min ( nval, 1, yvalue )
      ymax = r4mat_max ( nval, 1, yvalue )
c
c  IF(NCHAR.LE.0) GET LETTER TO USE FOR DOTS
c
      if ( nchar .gt. 1 .and. nchar .lt. nval ) nchar = 1
      if ( nchar .le. 0 ) node_use = node
      if ( nchar .eq. 1 ) node_use = ichars(1)
c
c  WRITING GRAPH HEADINGS
c
c  IF(IGRID.NE.0) COMPUTE AXIS LABELS
c
      if ( 0 .lt. igrid ) then
        nlabel = ((ixmax-1)/10)+1
        do i = 1,nlabel
          ix = 10*(i-1)+ixmin
          call ntor ( ix, ixmax, ixmin, x, xmax, xmin )
          xlabel(i) = x
        end do
      end if
c
c  Drawing graph
c
      do i = 1,iymax

        iy = iymax+1-i
        im1 = i-1
        call ntor ( iy, iymax, iymin, y, ymax, ymin )

        if ( ibak .eq. ' ' ) then
          ixlin = 0
        else
          ixlin = ixmax
        end if
 
        call do_back ( ixmin, ixmax, icolum )

        call do_grid ( igrid, icolum, im1, ixlin, ixmax, ixmin )

        do k = 1,nval
          y1 = yvalue(k)
          call rton ( iy1, iymax, iymin, y1, ymax, ymin, 0.5 )
          if ( iy1 .eq. iy ) then
            x = xvalue(k)
            call rton ( ix, ixmax, ixmin, x, xmax, xmin, 0.5 )
            if(nchar.gt.1)node_use = ichars(k)
            icolum(ix) = node_use
            if(ix.gt.ixlin)ixlin = ix
          end if
        end do

        if ( igrid .le. 0 ) then
          nstore = 0
        else
          nstore = 9
          write(istore, '(g9.3)' )y
          read(istore,'(9a1)')jstore
        end if

        output(1:margl) = ' '

        nlen = margl + 1

        do j = 1,nstore
          nlen = nlen+1
          output(nlen:nlen) = jstore(j)
        end do

        do j = 1,ixlin
          nlen = nlen+1
          output(nlen:nlen) = icolum(j)
        end do

        write ( *, '(a)' ) output(1:nlen)

      end do
c
c  Writing bottom axis
c
      if(igrid.gt.0)then
        write ( *, '(10x,15(a1,9x))' )  (ichara,i = 1,nlabel)
        write ( *, '(5x,15(1x,g9.3))' )  (xlabel(i),i = 1,nlabel)
      end if

      return
      end
      subroutine getscl ( )

c*********************************************************************72
c
cc GETSCL prints out the gray scaled used in contour plots.
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
c  Parameters:
c
c    None
c
      implicit none

      integer i
      character*1 ibak
      integer icopy
      character*1 ihoriz
      integer ip1
      character*1 islant
      character*1 iverti
      character*1 lgray(30)
      integer margl
      integer ncolr
      integer ngray
      integer nrowr
      character*1 node
      character output*100
      real zhi
      real zlo
      real zmaxw
      real zminw

      common /chdata/ ibak,ihoriz,islant,iverti,node
      common /gracom/ margl,ncolr,ngray,nrowr,zmaxw,zminw
      common /chrgra/ lgray

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Gray scale equivalents:'
      write ( *, '(a)' ) '  Symbol  Low value     High value'
      write ( *, '(a)' ) ' '
      write ( *, '(3x,a1,2x,14x,g14.6)' ) lgray(1), zminw

      do i = 2,ngray-1
        icopy = i
        call ntor ( icopy, ngray, 2, zlo, zmaxw, zminw )
        ip1 = i + 1
        call ntor ( ip1, ngray, 2, zhi, zmaxw, zminw )
        write ( *, '(3x,a1,2x,2g14.6)' ) lgray(i), zlo, zhi
      end do

      write ( *, '(3x,a1,2x,g14.6)') lgray(ngray), zmaxw

      return
      end
      subroutine grafun ( xmax, xmin, ymax, ymin, zname )

c*********************************************************************72
c
cc GRAFUN makes a grayscale plot of Z(X,Y), data from subroutine.
c
c  Discussion:
c
c    before calling grafun-
c
c    write a subroutine of the form subroutine zfunc(x,y,z)
c    which accepts values x and y, and evaluates z(x,y).
c
c    decide on the maximum and minimum z value you want to
c    display, and call setwin(winmax,winmin) to tell the
c    program about this.
c
c    decide how large the region is that you want to display
c    in columns and rows, and call setreg(ncolr,nrowr).
c
c    choose a gray scale, and if igray = 1, 2, or 3, just call
c    setgra(igray,lgray,ngray).  if you want to define your own
c    gray scale, then store up to 30 characters, one per entry,
c    in the vector lgray, set ngray to the number of characters
c    in lgray, and then call setgra(igray,lgray,ngray).
c    lgray is of type character*1.
c
c    then decide the x and y boundaries of the your region, declare
c    your subroutine external and call grafun(xmax,xmin,ymax,ymin,zfunc)
c
c    afterwards you can call getscl to see what the gray scale
c    numeric conversion table is.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    19 May 2008
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, real XMAX, XMIN, the maximum and minimum X values to plot.
c
c    Input, real YMAX, YMIN, the maximum and minimum Y values to plot.
c
c  ZNAME      - THE NAME OF THE USER WRITTEN SUBROUTINE WHICH
c               HAS THE FORM SUBROUTINE ZFUNC(X,Y,Z),
c               WHICH EVALUATES Z(X,Y).  THE NAME MUST BE
c               DECLARED EXTERNAL IN THE CALLING PROGRAM.
c
      implicit none

      integer i
      integer iback
      character*1 ibak
      character*1 ihoriz
      character*1 islant
      character*1 iverti
      integer j
      integer jcopy
      character*1 lgray(30)
      integer margl
      integer ncolr
      integer ngray
      character*1 node
      integer nrowr
      real xmax
      real xmin
      real xval
      real ymax
      real ymin
      real yval
      real zmaxw
      real zminw
      external zname
      real zval
      real zvec(100)

      common /chdata/ ibak,ihoriz,islant,iverti,node
      common /gracom/ margl,ncolr,ngray,nrowr,zmaxw,zminw
      common /chrgra/ lgray

      do i = 1, nrowr
        iback = nrowr+1-i
        call ntor ( iback, nrowr, 1, yval, ymax, ymin )
        do j = 1, ncolr
          jcopy = j
          call ntor ( jcopy, ncolr, 1, xval, xmax, xmin )
          call zname ( xval, yval, zval )
          zvec(j) = zval
        end do
        call gravec ( ncolr, zvec )
      end do

      return
      end
      subroutine graray ( ncolz, nrowz, zarray )

c*********************************************************************72
c
cc GRARAY makes a grayscale plot of Z(X,Y), data from array.
c
c  Discussion:
c
c    before calling graray-
c
c    dimension an array zarray of some sufficient size
c    nrowz rows by ncolz columns, where nrowz should
c    be no less than nrowr, and ncolz no less than ncolr,
c    the sizes of the region to be plotted.
c    store the values of z(x,y) in this array.
c
c    decide on the maximum and minimum z value you want to
c    display, and call setwin(winmax,winmin) to tell the
c    program about this.
c    if you want to use the maximum and minimum values
c    in the array, you can call getmax(ncolz,nrowz,zarray,zmax,zmin)
c    to get these values.
c
c    decide how large the region is that you want to display
c    in columns and rows, and call setreg(ncolr,nrowr).
c
c    choose a gray scale, and if igray = 1, 2, or 3, just call
c    setgra(igray,lgray,ngray).  if you want to define your own
c    gray scale, then store up to 30 characters, one per entry,
c    in the vector lgray, set ngray to the number of characters
c    in lgray, and then call setgra(igray,lgray,ngray).
c    lgray is of type character*1.
c
c    then call graray(ncolz,nrowz,zarray).
c  
c    afterwards you can call getscl to see what the gray scale
c    numeric conversion table is.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    19 May 2008
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer NOUNIT, the dimension of IOUNIT.  This will be either 2 or 3.
c
      implicit none

      integer ncolz
      integer nrowz

      integer i
      integer j
      integer margl
      integer ncolr
      integer ngray
      integer nrowr
      real zarray(nrowz,ncolz)
      real zmaxw
      real zminw
      real zvec(100)

      common /gracom/ margl,ncolr,ngray,nrowr,zmaxw,zminw

      do i = 1, nrowr
        do j = 1, ncolr
          zvec(j) = zarray(i,j)
        end do
        call gravec ( ncolr, zvec )
      end do

      return
      end
      subroutine gravec ( nvec, zvec )

c*********************************************************************72
c
cc GRAVEC makes a grayscale plot of Z(X,Y), data from vector.
c
c  Discussion:
c
c    BEFORE CALLING GRAVEC-
c
c    Dimension a vector zvec of some sufficient size nvec.
c    store the values for a single row of your plot in this
c    vector.  several rows can be drawn by storing one
c    row at a time and calling gravec.
c
c    Decide on the maximum and minimum z value you want to
c    display, and call setwin(winmax,winmin) to tell the
c    program about this.
c
c    Decide how large the region is that you want to display
c    in columns and rows, and call setreg(ncolr,nrowr).
c
c    Choose a gray scale, and if igray = 1, 2, or 3, just call
c    setgra(igray,lgray,ngray).  if you want to define your own
c    gray scale, then store up to 30 characters, one per entry,
c    in the vector lgray, set ngray to the number of characters
c    in lgray, and then call setgra(igray,lgray,ngray).
c    lgray is of type character*1.
c
c    Then call gravec(nvec,zvec).
c
c    Afterwards you can call getscl to see what the gray scale
c    numeric conversion table is.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    19 May 2008
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer NVEC, the size of the data vector.
c
      implicit none

      integer nvec

      character*1 ibak
      integer icol
      character*1 ihoriz
      integer ilen
      integer irow
      character*1 islant
      character*1 iverti
      integeriz
      integer j
      integer jhi
      character lgray*1
      character line(100)*1
      integer margl
      integer ncolr
      integer ngray
      character*1 node
      integer nrowr
      character output*100
      real zvec(nvec)
      real zmaxw
      real zminw
      real zval

      common /chdata/ ibak,ihoriz,islant,iverti,node
      common /gracom/ margl,ncolr,ngray,nrowr,zmaxw,zminw
      common /chrgra/ lgray(30)

      ilen = 0
      irow = 1
      jhi = margl + ncolr

      call do_back ( 1, jhi, line )

      do j = 1, ncolr

        icol = j
        zval = zvec(icol)

        if ( zval .ge. zminw .or. lgray(1) .ne. ' ' ) then

          if ( zval .le. zmaxw .or. lgray(ngray) .ne.' ' ) then
            ilen = j+margl
            call rton ( iz, ngray, 2, zval, zmaxw, zminw, 0.0 )
            line(ilen) = lgray(iz)
          end if

        end if

      end do

      write ( *, '(100a1)' ) ( line(j), j = 1, ilen )

      return
      end
      subroutine grdplt ( ixmax, iymax, nval, x1, x2, y1, y2 )

c*********************************************************************72
c
cc GRDPLT plots a collection of line segments.
c
c  Discussion:
c
c    The routine plots, on a nongraphic terminal, a
c    graph of the network of points and lines supplied by
c    the user.  these points might represent a grid used
c    in a finite element problem, or an electric circuit,
c    or an arbitrary mathematic graph.
c
c    the program requests only the x,y coordinates of the
c    two endpoints of each line, the number of such lines, and
c    the size (in rows and columns of output) of the total
c    graph to be produced. 
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    19 May 2008
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c  ixmax  - the requested width, in columns, of the output graph.
c           ixmax should be between 1 and 100.
c
c  iymax  - the requested height, in rows, of the output graph.
c           iymax should be greater than or equal to 1.
c
c           be sure that your requested value of ixmax
c           is compatible with your output device.  many terminals
c           can handle at most 72 or 80 columns.  the lineprinter
c           can handle up to 100.
c
c           also, note that a square 4 rows by 4 columns actually
c           appears about twice as tall as wide.  hence, you might
c           get less distorted results using ixmax = 2*iymax.
c
c  nval   - the number of line segments to be drawn.  there
c           should be this many valid entries in all of x1,
c           x2, y1 and y2.  nval may not be less than 2.
c
c  x1     - a vector of length at least nval, containing the
c           x coordinates of the starting points of the line.
c
c  x2     - a vector of length at least nval, containing the
c           x coordinates of the end points of the line.
c
c  y1     - a vector of length at least nval, containing the
c           y coordinates of the starting points of the line.
c
c  y2     - a vector of length at least nval, containing the
c           y coordinates of the end points of the line.
c
      implicit none

      integer nval

      integer i
      character*1 ibak
      character*1 iblank
      character*1 icolum(100)
      integer igrid
      character*1 ihoriz
      integer im1
      character*1 islant
      character*1 iverti
      integer ix
      integer ix1
      integer ix2
      integer ixhi
      integer ixlin
      integer ixlo
      integer ixmax
      integer ixmin
      integer iy
      integer iy1
      integer iy2
      integer iymax
      integer iymin
      integer j
      character*1 lgray(30)
      integer margl
      integer ncolr
      integer ngray
      character*1 node
      integer nrowr
      character*100 output
      real r4mat_max
      real r4mat_min
      real x1(nval)
      real x1j
      real x1max
      real x1min
      real x2(nval)
      real x2max
      real x2min
      real x2j
      real xi
      real xmax
      real xmin
      real y1(nval)
      real y1j
      real y1max
      real y1min
      real y2(nval)
      real y2j
      real y2max
      real y2min
      real yi
      real ymax
      real ymin
      real zmaxw
      real zminw

      common /chdata/ ibak,ihoriz,islant,iverti,node
      common /gracom/ margl,ncolr,ngray,nrowr,zmaxw,zminw
      common /chrgra/ lgray

      iblank = ' '
      ixmin = 1
      iymin = 1
c
c  Check NVAL, IXMAX, IYMAX
c
      if ( nval .lt. 2 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'GRDPLT - Warning!'
        write ( *, '(a)' ) '  NVAL < 2.'
        write ( *, '(a)' ) '  This call is ignored.'
        return
      end if

      if ( ixmax .lt. 1 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'GRDPLT - Warning!'
        write ( *, '(a)' ) '  IXMAX < 1.'
        write ( *, '(a)' ) '  This call is ignored.'
        return
      end if

      if ( 100 .lt. ixmax ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'GRDPLT - Warning!'
        write ( *, '(a)' ) '  100 < IXMAX.'
        write ( *, '(a)' ) '  This call is ignored.'
        return
      end if

      if ( iymax .lt. 1 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'GRDPLT - Warning!'
        write ( *, '(a)' ) '  IYMAX < 1.'
        write ( *, '(a)' ) '  This call is ignored.'
        return
      end if
c
c  Calculating upper and lower limits on Y
c  and on X.
c
      x1min = r4mat_min ( nval, 1, x1 )
      x1max = r4mat_max ( nval, 1, x1 )
      y1min = r4mat_min ( nval, 1, y1 )
      y1max = r4mat_max ( nval, 1, y1 )

      x2min = r4mat_min ( nval, 1, x2 )
      x2min = r4mat_max ( nval, 1, x2 )
      y2min = r4mat_min ( nval, 1, y2 )
      y2min = r4mat_max ( nval, 1, y2 )

      xmax = max ( x1max, x2max )
      xmin = min ( x1min, x2min )
      ymax = max ( y1max, y2max )
      ymin = min ( y1min, y2min )
c
c  Begin at maximum Y
c
      do i = 1, iymax

        iy = iymax+1-i
        im1 = i-1
        call ntor ( iy, iymax, iymin, yi, ymax, ymin )

        if ( ibak .eq. ' ' ) then 
          ixlin = 0
        else
          ixlin = ixmax
        end if

        call do_back ( ixmin, ixmax, icolum )

        call do_grid ( igrid, icolum, im1, ixlin, ixmax, ixmin )
c
c  Check if any lines cross this line
c
        do j = 1, nval

          x1j = x1(j)
          y1j = y1(j)
          x2j = x2(j)
          y2j = y2(j)

          call rton ( ix1, ixmax, ixmin, x1j, xmax, xmin, 0.5 )
          call rton ( iy1, iymax, iymin, y1j, ymax, ymin, 0.5 )
          call rton ( ix2, ixmax, ixmin, x2j, xmax, xmin, 0.5 )
          call rton ( iy2, iymax, iymin, y2j, ymax, ymin, 0.5 )

          if(iy.lt.iy1.and.iy.lt.iy2)go to 80
          if(iy.gt.iy1.and.iy.gt.iy2)go to 80
          if(iy.eq.iy1.and.iy.eq.iy2)go to 60

          call icross ( xi, xmax, xmin, x1j, x2j, yi, y1j, y2j )
          call rton ( ix, ixmax, ixmin, xi, xmax, xmin, 0.5 )

          if(ix1.eq.ix2)icolum(ix) = iverti
          if(ix1.ne.ix2)icolum(ix) = islant
          if(ix.gt.ixlin)ixlin = ix
          go to 80
c
c  Horizontal line
c
   60     continue
          ixlo = min ( ix1, ix2 )
          ixhi = max ( ix1, ix2 )
          if(ixhi.gt.ixlin)ixlin = ixhi
          do ix = ixlo,ixhi
            icolum(ix) = ihoriz
          end do
   80     continue

        end do
c
c  See if any nodes are on this line.
c
        do j = 1,nval

          y1j = y1(j)
          call rton ( iy1, iymax, iymin, y1j, ymax, ymin, 0.5 )

          if(iy1.eq.iy)then
            x1j = x1(j)
            call rton ( ix1, ixmax, ixmin, x1j, xmax, xmin, 0.5 )
            icolum(ix1) = node
            if(ix1.gt.ixlin)ixlin = ix1
          end if

          y2j = y2(j)
          call rton ( iy2, iymax, iymin, y2j, ymax, ymin, 0.5 )

          if(iy2.eq.iy) then
            x2j = x2(j)
            call rton ( ix2, ixmax, ixmin, x2j, xmax, xmin, 0.5 )
            icolum(ix2) = node
            if(ix2.gt.ixlin)ixlin = ix2
          end if

        end do

        write ( *, '(100a1)' ) 
     &    (iblank,j = 1,margl), (icolum(j),j = 1,ixlin)

      end do

      return
      end
      subroutine icross ( xc, xmax, xmin, x1, x2, yc, y1, y2 )

c*********************************************************************72
c
cc ICROSS determines the X coordinate where a line has a given Y coordinate.
c
c  Discussion:
c
c    The line is described as:
c
c      ( XC - X2 ) * ( Y1 - Y2 ) = ( X1 - X2 ) * ( YC - Y2 )
c
c    where the points (X1,Y1) and (X2,Y2) are known, the value
c    YC is given, and the value of XC is desired.
c
c    If necessary, XC is increased to XMIN or decreased to XMAX.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    19 May 2008
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, real XC, the computed X coordinate.
c
c    Input, real XMAX, XMIN, the maximum and minimum X values allowed.
c
c    Input, real X1, X2, the X coordinates of two points on the line.
c
c    Input, real YC, the Y coordinate of the unknown point.
c
c    Input, real Y1, Y2, the Y coordinates of two points on the line.
c
      implicit none

      real r4_huge
      real x1
      real x2
      real xc
      real xmax
      real xmin
      real y1
      real y2
      real yc

      if ( y1 .eq. y2 ) then

        if ( yc .ne. y1 ) then
          xc = r4_huge ( )
          return
        end if

        xc = x1

      else

        xc = x2 + ( x1 - x2 ) * ( ( yc - y2 ) / ( y1 - y2 ) )

      end if
c
c  Force XC to be within [XMIN,XMAX].
c
      xc = max ( xc, xmin )
      xc = min ( xc, xmax )

      return
      end
      subroutine ntor ( ix, ixmax, ixmin, x, xmax, xmin )

c*********************************************************************72
c
cc NTOR maps an integer to a real.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    19 May 2008
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
      implicit none

      integer ix
      integer ixmax
      integer ixmin
      real rmax
      real rmin
      real x
      real xmax
      real xmin

      rmax = real(ix-ixmin)*xmax
      rmin = real(ixmax-ix)*xmin
      x = (rmin+rmax)/real(ixmax-ixmin)
      if(x.lt.xmin)x = xmin
      if(x.gt.xmax)x = xmax

      return
      end
      function r4_huge ( )

c*********************************************************************72
c
cc R4_HUGE returns a "huge" R4.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    13 April 2004
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, real R4_HUGE, a huge number.
c
      implicit none

      real r4_huge

      r4_huge = 1.0E+30

      return
      end
      function r4mat_max ( m, n, a )

c*********************************************************************72
c
cc R4MAT_MAX returns the maximum entry of an R4MAT.
c
c  Discussion:
c
c    An R4MAT is an array of real values.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    19 May 2008
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer M, the number of rows in A.
c
c    Input, integer N, the number of columns in A.
c
c    Input, real A(M,N), the matrix.
c
c    Output, real R4MAT_MAX, the maximum entry of A.
c
      implicit none

      integer m
      integer n

      real a(m,n)
      integer i
      integer j
      real r4mat_max
      real value

      value = a(1,1)
      do j = 1, n
        do i = 1, m
          value = max ( value, a(i,j) )
        end do
      end do

      r4mat_max = value

      return
      end
      function r4mat_min ( m, n, a )

c*********************************************************************72
c
cc R4MAT_MIN returns the minimum entry of an R4MAT.
c
c  Discussion:
c
c    An R4MAT is an array of real values.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    19 May 2008
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer M, the number of rows in A.
c
c    Input, integer N, the number of columns in A.
c
c    Input, real A(M,N), the matrix.
c
c    Output, real R4MAT_MIN, the minimum entry of A.
c
      implicit none

      integer m
      integer n

      real a(m,n)
      integer i
      integer j
      real r4mat_min
      real value

      value = a(1,1)
      do j = 1, n
        do i = 1, m
          value = min ( value, a(i,j) )
        end do
      end do

      r4mat_min = value

      return
      end
      subroutine rton ( ix, ixmax, ixmin, x, xmax, xmin, xround )

c*********************************************************************72
c
cc RTON maps a real to an integer.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    19 May 2008
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, integer IX, the corresponding integer value.
c
c    Input, integer IXMAX, IXMIN, the maximum and mininum values
c    in the integer interval.
c
c    Input, real X, the real value to be mapped.
c
c    Input, real XMAX, XMIN, the maximum and minimum values
c    in the real interval.
c
c    Input, real XROUND, a rounding value for real numbers.
c
      implicit none

      integer ix
      integer ixmax
      integer ixmin
      real ratio
      real rmax
      real rmin
      real x
      real xave
      real xmax
      real xmin
      real xround

      if(xmax.gt.xmin)go to 10
      xave = 0.5*(xmin+xmax)
      xmin = xave-1.0
      xmax = xave+1.0
   10 continue
      rmin = (xmax-x)*real(ixmin)
      rmax = (x-xmin)*real(ixmax)
      ratio = (rmin+rmax)/(xmax-xmin)
      ix = ifix(ratio+xround)
      if(ix.ge.ixmin)go to 22
      if(ix.le.1.and.ixmin.eq.2.and.xround.eq.0.0)go to 21
      ix = ixmin
      go to 22
21    continue
      ix = 1
22    continue
      if(ix.gt.ixmax)ix = ixmax

      return
      end
      subroutine setbak ( ibak_new )

c*********************************************************************72
c
cc SETBAK sets the value of the background character.
c
c  Discussion:
c
c    All plots have a background.  This background is simply a character
c    which is assigned as the initial value to each plot "pixel".
c    The user may replace the default value by calling this routine.
c    The default value is the blank character.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    19 May 2008
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, character*1 IBAK_NEW, the new value of the background character.
c
      implicit none

      character*1 ibak
      character*1 ibak_new
      character*1 ihoriz
      character*1 islant
      character*1 iverti
      character*1 node

      common /chdata/ ibak,ihoriz,islant,iverti,node

      ibak = ibak_new

      return
      end
      subroutine setdef

c*********************************************************************72
c
cc SETDEF sets all internal configuration data to default values.
c
c  Discussion:
c
c    Before calling any TTYPLT graphics routines, the user should call
c    SETDEF to set default values for the internal data.
c
c  IBAK = ' '
c  IGRID = 0
c  IHORIZ = '-'
c  ISLANT = '*'
c  IVERTI = '!'
c  NODE = 'N'
c
c  LGRAY(2) THROUGH LGRAY(15) CONTAIN THE FOLLOWING characterS
c  .,:;- = +/YX&@$#
c  WITH LGRAY(1) AND LGRAY(16) BLANK.
c  MARGL = 0
c  NCOLR = 0
c  NGRAY = 16
c  NROWR = 0
c  ZMAXW = 0.0
c  ZMINW = 0.0
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    19 May 2008
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    None
c
      implicit none

      character*1 ibak
      integer igray
      integer igrid
      character*1 ihoriz
      character*1 islant
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

      ibak = ' '
      igrid = 0
      ihoriz = '-'
      islant = '*'
      iverti = '!'
      node = 'n'
      ngray = 16
      igray = 1
      call setgra ( igray, lgray, ngray )
      margl = 0
      ncolr = 0
      nrowr = 0
      zmaxw = 0.0
      zminw = 0.0

      return
      end
      subroutine setgra ( igrain, lgrain, ngrain )

c*********************************************************************72
c
cc SETGRA sets the gray scale.
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
c  Parameters:
c
c    Input, integer IGRAIN, the gray scale choice.
c
c    Input, character*1 LGRAIN(NGRAIN), the new gray scale, if IGRAIN = 0.
c
c    Input, integer NGRAIN, the number of gray scale entries, if IGRAIN = 0.
c
      implicit none

      integer i
      integer igrain
      integer igray
      character*1 lgrain(1)
      character*1 lgray(30)
      integer margl
      integer ncolr
      integer ngrain
      integer ngray
      integer nrowr
      real zmaxw
      real zminw

      common /gracom/ margl,ncolr,ngray,nrowr,zmaxw,zminw
      common /chrgra/ lgray

      igray = igrain
      if(igray.eq.1)go to 20
      if(igray.eq.2)go to 30
      if(igray.eq.3)go to 40
c
c  USER SUPPLIES GRAY VALUES
c
      if(ngrain.lt.1)then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'SETGRA - Warning!'
        write ( *, '(a)' ) '  The value of NGRAIN was not recognized.'
        write ( *, '(a)' ) '  The input was ignored.'
        return
      end if

      ngray = min0(30,ngrain)

      do i = 1,ngray
        lgray(i) = lgrain(i)
      end do

      return
c
c  GRAY SCALE
c
   20 continue
      ngray = 16
      lgray(1) = ' '
      lgray(2) = '.'
      lgray(3) = ','
      lgray(4) = ':'
      lgray(5) = ';'
      lgray(6) = '-'
      lgray(7) = '='
      lgray(8) = '+'
      lgray(9) = '/'
      lgray(10) = 'Y'
      lgray(11) = 'X'
      lgray(12) = '&'
      lgray(13) = '@'
      lgray(14) = '$'
      lgray(15) = '#'
      lgray(16) = ' '
      return
c
c  NUMERIC SCALE
c
   30 continue
      ngray = 12
      lgray(1) = ' '
      lgray(2) = '0'
      lgray(3) = '1'
      lgray(4) = '2'
      lgray(5) = '3'
      lgray(6) = '4'
      lgray(7) = '5'
      lgray(8) = '6'
      lgray(9) = '7'
      lgray(10) = '8'
      lgray(11) = '9'
      lgray(12) = ' '
      return
c
c  ALPHABETICAL SCALE
c
   40 continue
      ngray = 28
      lgray(1) = ' '
      do i = 2,ngray-1
        lgray(i) = char(i-2+ichar('a'))
      end do
      lgray(ngray) = ' '

      return
      end
      subroutine setigr ( igrid_new )

c*********************************************************************72
c
cc SETIGR sets the value of IGRID.
c
c  Discussion:
c
c    IGRID, controls the drawing of the grid.
c    0 or less, no grid is drawn.
c    1 or more, grid nodes are drawn, and numeric tick marks.
c    2 or more, grid lines are drawn.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    19 May 2008
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer IGRID_NEW, the new value for IGRID.
c
      implicit none

      integer igrid
      integer igrid_new

      common /indata/ igrid

      igrid = igrid_new

      return
      end
      subroutine setlet ( ihnew, isnew, ivnew )

c*********************************************************************72
c
cc SETLET sets the letters used to mark horizontal, vertical and slanted lines.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    19 May 2008
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
      implicit none

      character*1 ibak
      character*1 ihnew
      character*1 ihoriz
      character*1 islant
      character*1 isnew
      character*1 iverti
      character*1 ivnew
      character*1 node

      common /chdata/ ibak,ihoriz,islant,iverti,node

      ihoriz = ihnew
      islant = isnew
      iverti = ivnew

      return
      end
      subroutine setmar ( margin )

c*********************************************************************72
c
cc SETMAR resets the left plot margin.
c
c  Discussion:
c
c    By default, plots begin in column 1, and this is interpreted as
c    though the left margin of the plot, MARGL, was set at 0.
c
c    This routine allows the user to reset the value of MARGL to any
c    positive value.
c
c    Setting MARGL to 10 means that plots will be drawn with a left
c    margin of 10 blank spaces, for instance.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    19 May 2008
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer MARGIN, the new value for the left margin.
c
      implicit none

      integer margin
      integer margl
      integer ncolr
      integer ngray
      integer nrowr
      real zmaxw
      real zminw

      common /gracom/ margl,ncolr,ngray,nrowr,zmaxw,zminw

      if ( margin .lt. 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'SETMAR - Warning:'
        write ( *, '(a)' ) '  0 <= MARGL is required.'
        write ( *, '(a,i8)' ) '  Ignoring request to set it to ', margin
        return
      end if

      margl = margin

      return
      end
      subroutine setnod ( node_new )

c*********************************************************************72
c
cc SETNOD sets the character used to mark data values in dot plots.
c
c  Discussion:
c
c    For dot plots, the location of each "dot" is marked by a particular
c    character stored in NODE.  The default value is 'N'.  The user
c    can reset this character by calling this routine.
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
c  Parameters:
c
c    Input, integer NODE_NEW, the new value for NODE.
c
      implicit none

      character*1 ibak
      character*1 ihoriz
      character*1 islant
      character*1 iverti
      character*1 node
      character*1 node_new

      common /chdata/ ibak,ihoriz,islant,iverti,node

      node = node_new

      return
      end
      subroutine setreg ( ncolin, nrowin )

c*********************************************************************72
c
cc SETREG
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    19 May 2008
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
      implicit none

      integer margl
      integer ncolin
      integer ncolr
      integer ngray
      integer nrowin
      integer nrowr
      real zmaxw
      real zminw

      common /gracom/ margl,ncolr,ngray,nrowr,zmaxw,zminw

      ncolr = ncolin
      nrowr = nrowin

      return
      end
      subroutine setwin ( winmax, winmin )

c*********************************************************************72
c
cc SETWIN
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    19 May 2008
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
      implicit none

      integer margl
      integer ncolr
      integer ngray
      integer nrowr
      real winmax
      real winmin
      real zmaxw
      real zminw

      common /gracom/ margl,ncolr,ngray,nrowr,zmaxw,zminw

      zmaxw = winmax
      zminw = winmin

      return
      end
      subroutine timestamp ( )

c*********************************************************************72
c
cc TIMESTAMP prints out the current YMDHMS date as a timestamp.
c
c  Discussion:
c
c    This FORTRAN77 version is made available for cases where the
c    FORTRAN90 version cannot be used.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    12 January 2007
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    None
c
      implicit none

      character * ( 8 ) ampm
      integer d
      character * ( 8 ) date
      integer h
      integer m
      integer mm
      character * ( 9 ) month(12)
      integer n
      integer s
      character * ( 10 ) time
      integer y

      save month

      data month /
     &  'January  ', 'February ', 'March    ', 'April    ', 
     &  'May      ', 'June     ', 'July     ', 'August   ', 
     &  'September', 'October  ', 'November ', 'December ' /

      call date_and_time ( date, time )

      read ( date, '(i4,i2,i2)' ) y, m, d
      read ( time, '(i2,i2,i2,1x,i3)' ) h, n, s, mm

      if ( h .lt. 12 ) then
        ampm = 'AM'
      else if ( h .eq. 12 ) then
        if ( n .eq. 0 .and. s .eq. 0 ) then
          ampm = 'Noon'
        else
          ampm = 'PM'
        end if
      else
        h = h - 12
        if ( h .lt. 12 ) then
          ampm = 'PM'
        else if ( h .eq. 12 ) then
          if ( n .eq. 0 .and. s .eq. 0 ) then
            ampm = 'Midnight'
          else
            ampm = 'AM'
          end if
        end if
      end if

      write ( *, 
     &  '(i2,1x,a,1x,i4,2x,i2,a1,i2.2,a1,i2.2,a1,i3.3,1x,a)' ) 
     &  d, month(m), y, h, ':', n, ':', s, '.', mm, ampm

      return
      end
