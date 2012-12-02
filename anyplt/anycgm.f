c  anycgm.f   13 April 1998
c
      subroutine anyplt ( icom )

c***********************************************************************
c
cc ANYPLT is an interface routine to a variety of graphics packages.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    24 January 2009
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer ICOM, specifies the graphics request being made.
c    0, enable graphics.
c    1, disable graphics.
c    2, begin plot.
c    3, define plot size.
c    4, move to a point.
c    5, draw to a point.
c    6, clear screen.
c    7, write string at position.
c    8, use virtual cursor.
c    9, end plot.
c    10, ring bell.
c    11, mark data.
c    12, return screen data.
c    13, return version.
c
      real angle
      character*80 carray
      real csize
      real cwide
      character*10 dev
      logical filled
      character*1 flush
      integer icmax
      integer icmin
      integer icom
      integer iplt1
      integer iplt2
      integer itable
      integer ixplt1
      integer ixplt2
      integer iyplt1
      integer iyplt2
      integer lent
      logical leqi
      integer marray
      integer nplot
      integer nval
      real pwide
      real xplt1
      real xplt2
      real xpmax
      real xpmin
      real xval(2)
      real yplt1
      real yplt2
      real ypmax
      real ypmin
      real yval(2)
c
      external leqi
c
      save dev
      save nplot
      save xpmax
      save xpmin
      save ypmax
      save ypmin
c
      common /anycom/ iplt1,iplt2,ixplt1,ixplt2,iyplt1,
     &                iyplt2,marray,xplt1,xplt2,yplt1,yplt2
      common /anychr/ carray
c
c  ICOM = 0  Enable graphics
c
      if ( icom .eq. 0 ) then
 
10      continue
 
        if ( leqi(carray,'cgm') ) then
          dev = 'cgmb'
        else if ( leqi(carray,'cgmb') ) then
          dev = 'cgmb'
        else if ( leqi(carray,'ps') ) then
          dev = 'ps'
        else if ( leqi(carray,'xws') ) then
          dev = 'xws'
        else
          write(*,*)'Enter graphics device.'
          write(*,*)'Choices are cgm, ps, xws'
          read(*,'(a)')carray
          go to 10
        end if
 
        write(*,*)'Using graphics device '//dev
 
        call device(dev)
 
        if ( dev .eq. 'cgmb' ) then
          call outfil('anyplt.cgm')
        else if ( dev .eq. 'ps' ) then
          call outfil('anyplt.ps')
        end if
 
        icmax = 200
        icmin = 2
        itable = 1
        call settab(icmax,icmin,itable)
 
        nplot = 0
c
c  ICOM = 1  Disable graphics
c
      else if ( icom .eq. 1 ) then
        call grfcls
c
c  ICOM = 2  Begin plot
c
      else if ( icom .eq. 2 ) then
 
        if ( nplot .eq. 0 ) then
          call grfini
 
          xval(1) = 0.0
          xval(2) = 1.0
          yval(1) = 0.0
          yval(2) = 1.0
          nval = 2
          call setscl(xval,yval,nval)
          nplot = 1
        else
          call newfrm
          nplot = nplot+1
        end if
c
c  ICOM = 3  Define plot size
c
      else if ( icom .eq. 3 ) then
        xpmin = xplt1
        xpmax = xplt1+xplt2
        ypmin = yplt1
        ypmax = yplt1+yplt2
        xval(1) = xpmin
        xval(2) = xpmax
        yval(1) = ypmin
        yval(2) = ypmax
        nval = 2
        call setscl(xval,yval,nval)
c
c  ICOM = 4  Move to point
c
      else if ( icom .eq. 4 ) then
        call movcgm(xplt1,yplt1)
c
c  ICOM = 5  Draw to point
c
      else if ( icom .eq. 5 ) then
        call drwcgm(xplt1,yplt1)
c
c  ICOM = 6  Clear screen
c
      else if ( icom .eq. 6 ) then
c       call newfrm
c
c  ICOM = 7,  Write string at position
c
      else if ( icom .eq. 7 ) then
 
        angle = 0.0
        cwide = 0.025
        pwide = 1.0
        lent = lenchr ( carray )
        flush = 'c'
 
        call chrplt(angle,cwide,pwide,carray(1:lent),xplt1,yplt1,flush)
c
c  ICOM = 8  Use virtual cursor
c
      else if ( icom .eq. 8 ) then
c
c  ICOM = 9  End plot
c
      else if ( icom .eq. 9 ) then
 
        if ( dev .eq. 'xws' ) then
 
          call movcgm(xpmin,ypmin)
 
          do i = 1,100
            call drwcgm(xpmax,ypmin)
            call drwcgm(xpmax,ypmax)
            call drwcgm(xpmin,ypmax)
            call drwcgm(xpmin,ypmin)
          enddo
 
        end if
c
c  ICOM = 10  Ring bell
c
      else if ( icom .eq. 10 ) then
c
c  ICOM = 11  Mark data
c
      else if ( icom .eq. 11 ) then
 
        filled = .false.
        csize = 0.005
        call circle(xplt1,yplt1,csize,filled)
c
c  ICOM = 12  Return screen data
c
      else if ( icom .eq. 12 ) then
c
c  ICOM = 13  Return version
c
      else if ( icom .eq. 13 ) then
        carray = 'AnyPlt - Version 1.02  13 April 1998  CgmPlt'
c
c  Unknown value of ICOM.
c
      else
        write(*,*)'AnyPlt - Fatal errorc'
        write(*,*)'  Unknown value of ICOM = ',icom
        stop
      end if
 
      return
      end
      subroutine chrplt ( angle, cwide, pwide, string, x, y, flush )
c
c***********************************************************************
c
cc CHRPLT will put a character string onto a graphics image, at any
c  angle and at any size, by plotting the characters.  The plot is
c  assumed to be of size PWIDE by PHITE, although PHITE itself is
c  not input.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    24 January 2009
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c  ANGLE  Input, REAL ANGLE, the angle in degrees at which the
c         string is to be drawn.  0 is typical.  90 degrees would
c         cause the string to be written from top to bottom.
c
c  CWIDE  Input, REAL CWIDE, the width of the characters. This
c         is measured in the same units as the plot width PWIDE.
c         For PWIDE = 1, a plot size of 0.025 would be reasonable,
c         since 40 characters would fit, but 2.0 would be nonsense.
c
c  PWIDE  Input, REAL PWIDE, the width of the plot, in the same
c         units as CWIDE.
c
c  STRING Input, CHARACTER*(*) STRING, contains the information to
c         be plotted.
c
c  X,
c  Y      Input, REAL X, Y.  The coordinates of a point which
c         determines where the string is drawn.  The string will
c         be drawn starting at, centered or, or ending at (X,Y),
c         depending on the value of FLUSH.
c
c  FLUSH  Input, CHARACTER*(*) FLUSH, a string which tells
c         CHRPLT how to place the string.  Only the first
c         character of FLUSH is examined, and the case of
c         the character is not important.
c
c         'L' - the string will be drawn flush left.
c         'C' - the string will be centered.
c         'R' - the string will be drawn flush right.
c
      real pi
      parameter (pi = 3.14159265)
c
      real degrad
      parameter (degrad = pi/180.0)
c
      real angle
      real ca
      character*1 chrtmp
      real cwide
      character*(*) flush
      integer i
      integer iascii
      integer icr
      integer ifont(1617)
      integer ip
      integer ipen
      integer ipoint(95)
      integer iv
      integer lenchr
      integer nchar
      integer nmax
      integer nvec
      real pwide
      logical rotate
      real sa
      real scl2
      character*(*) string
      real x
      real xb
      real xc
      real xcopy
      real xnew
      real xold
      real xrot
      real xt
      real y
      real yb
      real yc
      real ycopy
      real ynew
      real yold
      real yrot
      real yt
c
      external lenchr
c
      save ifont
      save ipoint
c
c  IPOINT is a pointer array into IFONT.
c
c  IPOINT(I) records where the "strokes" for character I begin
c  in the IFONT array.
c
      data (ipoint(i),i = 1,95) /
     &   1,   3,  26,  45,  66, 102, 130, 156, 166, 186, 206, 222, 233,
     & 249, 255, 267, 273, 293, 306, 328, 353, 363, 383, 411, 423, 457,
     & 483, 506, 533, 541, 552, 560, 587, 625, 638, 665, 683, 699, 714,
     & 727, 754, 770, 786, 805, 818, 826, 838, 848, 868, 884, 909, 930,
     & 956, 967, 981, 989,1001,1012,1025,1035,1045,1051,1061,1069,1075,
     &1081,1108,1131,1149,1172,1194,1214,1243,1260,1284,1307,1323,1336,
     &1364,1381,1401,1424,1447,1464,1486,1499,1516,1524,1536,1547,1560,
     &1570,1586,1592,1608/
c
c  IFONT contains the strokes defining the various symbols.
c
      data (ifont(i),i =    1, 396)/
     & 1, 0, 2,10,11, 9,22,10,23,11,22,10,11, 0, 9, 7, 9, 9,11, 9,11, 7,
     & 9, 7, 0, 2, 8,17, 7,23, 9,23, 8,17, 0,14,17,13,23,15,23,14,17, 0,
     & 4, 9,23, 7, 7, 0,13,23,11, 7, 0, 5,17,15,17, 0, 5,13,15,13, 0, 3,
     &15,19,13,21, 9,21, 7,19, 7,17, 9,15,13,15,15,13,15,11,13, 9, 9, 9,
     & 7,11, 0, 9,23, 9, 7, 0,13,23,13, 7, 0, 3, 5,23, 9,23, 9,19, 5,19,
     & 5,23, 0,17,23, 5, 7, 0,13, 7,13,11,17,11,17, 7,13, 7, 0, 1,17, 7,
     & 7,17, 7,19, 9,21,13,21,15,19,15,17, 5,13, 5,11, 9, 7,13, 7,17,15,
     & 0, 1,10,17, 9,23,11,23,10,17, 0, 1,12,23,11,21,10,19, 9,17, 9,15,
     & 9,13,10,11,11, 9,12, 7, 0, 1,12,23,13,21,14,19,15,17,15,15,15,13,
     &14,11,13, 9,12, 7, 0, 3, 7,15,15,15, 0,13,19, 9,11, 0, 9,19,13,11,
     & 0, 2, 7,15,15,15, 0,11,19,11,11, 0, 1,11, 7, 9, 7, 9, 9,11, 9,11,
     & 7,11, 6,10, 4, 0, 1, 7,15,15,15, 0, 1, 9, 7, 9, 9,11, 9,11, 7, 9,
     & 7, 0, 1,15,23, 7, 7, 0, 1, 9,23,13,23,15,19,15,11,13, 7, 9, 7, 7,
     &11, 7,19, 9,23, 0, 2, 7,21, 9,23, 9, 7, 0, 7, 7,11, 7, 0, 1, 5,21,
     & 9,23,15,23,17,21,17,19,15,17, 7,13, 5,10, 5, 7,17, 7, 0, 2, 5,23,
     &17,23,15,17,13,15, 9,15, 0,13,15,17,13,17,10,14, 7, 8, 7, 5,10, 0,
     & 1,13, 7,13,23, 5,13,17,13, 0, 1,17,23, 5,23, 5,17,13,17,17,15,17,
     &11,13, 7, 9, 7, 5,11, 0, 1,17,19,13,23, 9,23, 5,19, 5,13, 9,15,13/
      data (ifont(i),i =  397, 792)/
     &15,17,13,17,11,13, 7, 9, 7, 5,11, 5,13, 0, 1, 5,19, 5,23,17,23,11,
     &15,11, 7, 0, 1, 8,15, 6,17, 6,21, 8,23,14,23,16,21,16,17,14,15, 8,
     &15, 5,13, 5, 9, 8, 7,14, 7,17, 9,17,13,14,15, 0, 1,17,17,15,15, 7,
     &15, 5,17, 5,21, 7,23,15,23,17,21,17,11,15, 7, 7, 7, 5,11, 0, 2, 9,
     &13, 9,15,11,15,11,13, 9,13, 0, 9, 7, 9, 9,11, 9,11, 7, 9, 7, 0, 2,
     & 9,13, 9,15,11,15,11,13, 9,13, 0,11, 7, 9, 7, 9, 9,11, 9,11, 7,11,
     & 6,10, 4, 0, 1,17,21, 5,15,17, 9, 0, 2, 7,15,15,15, 0, 7, 9,15, 9,
     & 0, 1, 5,21,17,15, 5, 9, 0, 2, 7,21, 9,23,13,23,15,21,15,19,11,15,
     &11,11, 0,10, 7,10, 9,12, 9,12, 7,10, 7, 0, 1,13, 7, 9, 7, 5,11, 5,
     &19, 9,23,13,23,17,19,17,11,15, 9,13,11,12,10,10,10, 9,11, 9,15,10,
     &16,12,16,13,15,13,11, 0, 2, 5, 7,11,23,17, 7, 0, 8,15,14,15, 0, 2,
     & 5, 7, 5,23,15,23,17,21,17,17,15,15, 5,15, 0,15,15,17,13,17, 9,15,
     & 7, 5, 7, 0, 1,17,19,13,23, 9,23, 5,19, 5,11, 9, 7,13, 7,17,11, 0,
     & 1, 5, 7, 5,23,13,23,17,19,17,11,13, 7, 5, 7, 0, 2,17,23, 5,23, 5,
     & 7,17, 7, 0, 5,15,12,15, 0, 2, 5, 7, 5,23,17,23, 0, 5,15,12,15, 0,
     & 2,17,19,13,23, 9,23, 5,19, 5,11, 9, 7,13, 7,17,11,17,15,13,15, 0,
     &17,11,17, 7, 0, 3, 5, 7, 5,23, 0, 5,15,17,15, 0,17,23,17, 7, 0, 3,
     & 9,23,13,23, 0,11,23,11, 7, 0, 9, 7,13, 7, 0, 2,15,23,15,11,12, 7/
      data (ifont(i),i =  793,1188)/
     & 8, 7, 5,11, 5,13, 0,13,23,17,23, 0, 2, 5, 7, 5,23, 0,17,23, 5,15,
     &17, 7, 0, 1, 5,23, 5, 7,17, 7, 0, 1, 5, 7, 5,23,11,11,17,23,17, 7,
     & 0, 1, 5, 7, 5,23,17, 7,17,23, 0, 1,17,19,13,23, 9,23, 5,19, 5,11,
     & 9, 7,13, 7,17,11,17,19, 0, 1, 5, 7, 5,23,13,23,17,21,17,17,13,15,
     & 5,15, 0, 2,17,19,13,23, 9,23, 5,19, 5,11, 9, 7,13, 7,17,11,17,19,
     & 0,13,11,17, 7, 0, 2, 5, 7, 5,23,13,23,17,21,17,17,13,15, 5,15, 0,
     &13,15,17, 7, 0, 1,17,19,13,23, 9,23, 5,20, 5,18, 9,15,13,15,17,12,
     &17,10,13, 7, 9, 7, 5,10, 0, 2, 5,23,17,23, 0,11,23,11, 7, 0, 1, 5,
     &23, 5,10, 8, 7,14, 7,17,10,17,23, 0, 1, 5,23,11, 7,17,23, 0, 1, 5,
     &23, 8, 7,11,17,14, 7,17,23, 0, 2, 5,23,17, 7, 0,17,23, 5, 7, 0, 2,
     & 5,23,11,13,17,23, 0,11,13,11, 7, 0, 1, 5,23,17,23, 5, 7,17, 7, 0,
     & 1,11,23, 7,23, 7, 7,11, 7, 0, 1, 7,23,15, 7, 0, 1, 7,23,11,23,11,
     & 7, 7, 7, 0, 1, 7,21,11,23,15,21, 0, 1, 5, 3,17, 3, 0, 1, 9,23,13,
     &19, 0, 2, 7,14, 9,15,13,15,15,14,15, 7, 0,15,12, 9,12, 7,11, 7, 8,
     & 9, 7,13, 7,15, 8, 0, 2, 7,23, 7, 7, 0, 7,13, 9,15,13,15,15,13,15,
     & 9,13, 7, 9, 7, 7, 9, 0, 1,15,13,13,15, 9,15, 7,13, 7, 9, 9, 7,13,
     & 7,15, 9, 0, 2,15,13,13,15, 9,15, 7,13, 7, 9, 9, 7,13, 7,15, 9, 0,
     &15,23,15, 7, 0, 1, 7,11,15,11,15,13,13,15, 9,15, 7,13, 7, 9, 9, 7/
      data (ifont(i),i = 1189,1584)/
     &13, 7,15, 9, 0, 3, 9, 7, 9,23,13,23,13,22, 0, 8,15,12,15, 0, 8, 7,
     &11, 7, 0, 2,15,13,13,15, 9,15, 7,13, 7, 9, 9, 7,13, 7,15, 9, 0,15,
     &13,15, 3,13, 1, 9, 1, 7, 3, 0, 2, 7, 7, 7,23, 0, 7,14, 9,15,13,15,
     &15,14,15, 7, 0, 3, 9,15,11,15,11, 7, 0, 9, 7,13, 7, 0, 9,17, 9,19,
     &11,19,11,17, 9,17, 0, 2, 9,15,11,15,11, 1, 7, 1, 7, 3, 0, 9,17,11,
     &17,11,19, 9,19, 9,17, 0, 3, 7, 7, 7,23, 0,15,15, 7,10, 0, 9,11,15,
     & 7, 0, 2, 9,23,11,23,11, 7, 0, 9, 7,13, 7, 0, 3, 7,15, 7, 7, 0, 7,
     &14, 8,15,10,15,11,14,11, 7, 0,11,14,12,15,14,15,15,14,15, 7, 0, 2,
     & 7, 7, 7,15, 0, 7,14, 9,15,13,15,15,14,15, 7, 0, 1, 7,13, 9,15,13,
     &15,15,13,15, 9,13, 7, 9, 7, 7, 9, 7,13, 0, 2, 7,13, 9,15,13,15,15,
     &13,15, 9,13, 7, 9, 7, 7, 9, 0, 7,14, 7, 1, 0, 2,15,13,13,15, 9,15,
     & 7,13, 7, 9, 9, 7,13, 7,15, 9, 0,15,14,15, 1, 0, 2, 7,15, 9,15, 9,
     & 7, 0, 9,13,11,15,13,15,15,13, 0, 1,15,13,13,15, 9,15, 7,13, 9,11,
     &13,11,15, 9,13, 7, 9, 7, 7, 9, 0, 2, 9,23, 9, 7,11, 7, 0, 7,17,11,
     &17, 0, 2, 7,15, 7, 9, 9, 7,13, 7,15, 9, 0,15,15,15, 7, 0, 1, 7,15,
     &11, 7,15,15, 0, 1, 7,15, 9, 7,11,11,13, 7,15,15, 0, 2, 7,15,15, 7,
     & 0, 7, 7,15,15, 0, 2, 7,15,11, 7, 0,15,15,10, 5, 7, 1, 0, 1, 7,15,
     &15,15, 7, 7,15, 7, 0, 1,11,23, 7,23, 9,17, 7,15, 9,13, 7, 7,11, 7/
      data (ifont(i),i = 1585,1617)/
     & 0, 1, 9,23, 9, 7, 0, 1, 7,23,11,23, 9,17,11,15, 9,13,11, 7, 7, 7,
     & 0, 1, 5,21, 7,23,15,21,17,23, 0/
c
      nchar = lenchr(string)
 
      if ( pwide .le. 0 ) then
        write(*,*)' '
        write(*,*)'ChrPlt - Serious errorc'
        write(*,*)'  The plot width PWIDE is negativec'
        write(*,*)'  PWIDE = ',pwide
        return
      end if
c
c  Chop titles that are too long.  To do this, we need to know the
c  width of the plot (PWIDE) in same units as CWIDE.
c
      nmax = ifix(pwide/cwide)
      if ( nchar .gt. nmax ) then
        nchar = nmax
      end if
c
c  Shift string if centering or right flush option used.
c
      if ( flush(1:1) .eq. 'l'.or.flush(1:1).eq.'L' ) then
        xcopy = x
        ycopy = y
      else if ( flush(1:1) .eq. 'c'.or.flush(1:1).eq.'C' ) then
        xcopy = x-0.5*nchar*cwide*cos(angle*degrad)
        ycopy = y-0.5*nchar*cwide*sin(angle*degrad)
      else if ( flush(1:1) .eq. 'r'.or.flush(1:1).eq.'R' ) then
        xcopy = x-nchar*cwide*cos(angle*degrad)
        ycopy = y-nchar*cwide*sin(angle*degrad)
      else
        xcopy = x
        ycopy = y
      end if
c
c  Note that screen coordinates are used.
c  Thus a width of 0.1 is intended to mean 1/10 of screen size.
c
c  Set the scale factor for character height.
c
      scl2 = cwide/16.0
c
c  Set the starting point for the line of text, the lower left
c  corner of the first character.
c
c  Set the origin about which rotation is performed.
c
      xb = xcopy
      xrot = xcopy
      yb = ycopy
      yrot = ycopy
c
c  Get trig functions if rotation required, converting from
c  degrees to radians.
c
      if ( angle .eq. 0.0 ) then
        rotate = .false.
      else
        ca = cos(angle*degrad)
        sa = sin(angle*degrad)
        rotate = .true.
      end if
c
c  Loop over all characters in the string.
c
      do icr = 1,nchar
 
        xold = x
        yold = y
        xnew = x
        ynew = y
c
c  Get the ASCII code for the character and shift by 31 so that
c  the first printable character becomes code 1.
c
        chrtmp = string(icr:icr)
        iascii = ichar(chrtmp)-31
c
c  Replace any nonprintable characters with blanks.
c
        if ( (iascii.lt.1).or.(iascii.gt.95))iascii = 1
c
c  Get the pointer to this character in font table.
c
        ip = ipoint(iascii)
c
c  Get the number of "vectors" required to draw the character.
c  Here "vectors" means the number of times the pen is lowered, not
c  the number of pen strokes.
c
c  For blanks, this number is 1, due to the way the
c  algorithm is coded.
c
        nvec = ifont(ip)
c
c  Loop over all required pen movements.
c
        do iv = 1,nvec
          ipen = 3
          ip = ip+1
10        continue
          if ( ifont(ip) .eq. 0)go to 20
          xc = xb+(ifont(ip)-1)*scl2
          yc = yb+(ifont(ip+1)-7)*scl2
c
c  Apply rotation if necessary.
c
          if ( rotate ) then
            xt = xc-xrot
            yt = yc-yrot
            xc = ca*xt-sa*yt+xrot
            yc = sa*xt+ca*yt+yrot
          end if
c
c  Plot the pen stroke.
c
          if ( ipen .eq. 3 ) then
            xnew = xc
            ynew = yc
          else
            xold = xnew
            yold = ynew
            xnew = xc
            ynew = yc
c
c  Call the user supplied routine to draw a line from
c  (XOLD,YOLD) to (XNEW,YNEW).
c
            call movcgm(xold,yold)
            call drwcgm(xnew,ynew)
 
          end if
 
          ipen = 2
          ip = ip+2
          go to 10
 
20        continue
 
        enddo
c
c  Advance the base to compensate for character just drawn.
c
        xb = xb+cwide
 
      enddo
 
      return
      end
      subroutine settab ( icmax, icmin, itable )
c
c***********************************************************************
c
cc SETTAB replaces SETCTB, the DRAWCGM routine for setting up
c  the color tables.
c
c  For some reason, SETCTB does not work properly, at least on
c  the IRIS.  Colors beyond color number 201 are typically
c  mangled.
c
c  So SETTAB sets the colors between ICMIN and ICMAX, which
c  should typically be 2 and 201.
c
c  SETTAB will also set the values of color 0 to white, and
c  color 1 to black.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    24 January 2009
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c  ICMAX  Input, INTEGER ICMAX, the maximum color index to be
c         set.
c
c  ICMIN  Input, INTEGER ICMIN, the minimum color index to be
c         set.
c
c  ITABLE Input, INTEGER ITABLE, the desired table.
c         ITABLE should be between 1 and 6.
c
c         1: Gray scale, black to white
c         2: Blue to yellow
c         3: Waves
c         4: Pseudospectral
c         5: Inverted pseudospectral
c         6: Blue to red
c
      real pi
      parameter (pi = 3.14159265)
c
      real bval
      real gauss
      real gval
      integer i
      integer icmax
      integer icmin
      integer itable
      real rval
      real theta
c
      external gauss
c
c  1: Gray scale, ICMIN = black, ICMAX=white
c
      if ( itable .eq. 1 ) then
        do i = icmin,icmax
          bval = real(i-icmin)/real(icmax-icmin)
          gval = real(i-icmin)/real(icmax-icmin)
          rval = real(i-icmin)/real(icmax-icmin)
          call setclr(i,rval,gval,bval)
        enddo
c
c  2: Blue to yellow
c
      else if ( itable .eq. 2 ) then
        do i = icmin,icmax
          rval = real(i-icmin)/real(icmax-icmin)
          gval = real(i-icmin)/real(icmax-icmin)
          bval = (icmax-i)/real(icmax-icmin)
          call setclr(i,rval,gval,bval)
        enddo
c
c  3: Waves
c
      else if ( itable .eq. 3 ) then
        do i = icmin,icmax
          theta = 0.5*pi*real(i-icmin)/real(icmax-icmin)
          rval = cos(theta)**2
          bval = sin(theta)**2
          gval = 0.8*sin(10.0*theta)**6
          call setclr(i,rval,gval,bval)
        enddo
c
c  4: Pseudospectral
c
      else if ( itable .eq. 4 ) then
        do i = icmin,icmax
          theta = 4.0*real(i-icmin)/real(icmax-icmin)
          rval = gauss(theta-1.0)+gauss(theta-4.0)
          gval = gauss(theta-2.0)+gauss(theta-4.0)
          bval = gauss(theta-3.0)+gauss(theta-4.0)
          if ( rval.gt.1.0)rval = 1.0
          if ( gval.gt.1.0)gval = 1.0
          if ( bval.gt.1.0)bval = 1.0
          call setclr(i,rval,gval,bval)
        enddo
c
c  5: Inverted pseudospectral
c
      else if ( itable .eq. 5 ) then
        do i = icmin,icmax
          theta = 4.0*real(icmax-i)/real(icmax-icmin)
          rval = gauss(theta-1.0)+gauss(theta-4.0)
          gval = gauss(theta-2.0)+gauss(theta-4.0)
          bval = gauss(theta-3.0)+gauss(theta-4.0)
          if ( rval.gt.1.0)rval = 1.0
          if ( gval.gt.1.0)gval = 1.0
          if ( bval.gt.1.0)bval = 1.0
          call setclr(i,rval,gval,bval)
        enddo
c
c  6: Blue to red
c
      else if ( itable .eq. 6 ) then
        do i = icmin,icmax
          rval = real(i-icmin)/real(icmax-icmin)
          gval = 0.0
          bval = (icmax-i)/real(icmax-icmin)
          call setclr(i,rval,gval,bval)
        enddo
c
c  Unknown table.
c
      else
        write(*,*)' '
        write(*,*)'SetTab - Fatal errorc'
        write(*,*)'  Legal color table indices are '
        write(*,*)'  between 1 and 6.  Your value was ',itable
      end if
c
c  Background color 0 is to be white.
c
      i = 0
      rval = 1.0
      gval = 1.0
      bval = 1.0
      call setclr(i,rval,gval,bval)
c
c  Foreground color 1 is to be black.
c
      i = 1
      rval = 0.0
      gval = 0.0
      bval = 0.0
      call setclr(i,rval,gval,bval)
 
      return
      end
      function gauss(ratio)
c
c***********************************************************************
c
cc GAUSS is a simple function used by SETTAB.
c
c    GAUSS = EXP(-RATIO**2)
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    24 January 2009
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c  RATIO  Input, REAL RATIO, the parameter used in the formula.
c
c  GAUSS  Output, REAL GAUSS, the value of the gaussian variable.
c
      real gauss
      real ratio

      gauss = exp(-ratio**2)
 
      return
      end
