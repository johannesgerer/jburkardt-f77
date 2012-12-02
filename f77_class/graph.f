      program circ

c*********************************************************************72
c
      integer n
      parameter (n=21)

      integer i
      integer nx
      integer ny
      real x(0:n)
      real y(0:n)

      call gstart ( )

      do i=1,n
        x(i)=real(i-1)/10.0
        y(i)=x(i)**2
      end do
c
c  Plot a set of reference points.
c
      call refpnt ( )
c
c  Choose the marker.
c
      call marker('*')
c
c  Make the plot just large enough to hold the data.
c
      call size(n,x,y)
c
c  Draw a curve through the points.
c
      call curve(n,x,y)
c
c  Wait for the user to hit RETURN.
c
      call wait ( )
c
c  Draw a simple grid.
c
      nx=11
      ny=11
      call cgrid(nx,ny)

      call gstop
      stop
      end
      subroutine capchr(string)

c*********************************************************************72
c
c  CAPCHR accepts a STRING of characters and replaces any lowercase
c  letters by uppercase ones.
c
c  Compare LOWCHR which lowercases a string.
c
c  STRING Input/output, CHARACTER*(*) STRING, the string of
c         characters to be transformed.
c
      integer i
      integer itemp
      integer nchar
      character*(*) string

      intrinsic char
      intrinsic ichar
      intrinsic len

      nchar=len(string)

      do i=1,nchar
        itemp=ichar(string(i:i))
        if(97.le.itemp.and.itemp.le.122)then
          string(i:i)=char(itemp-32)
        endif

      end do

      return
      end
      subroutine cgrid(nx,ny)

c*********************************************************************72
c
c  CGRID is the simplest routine for drawing a grid of horizontal
c  and vertical lines.
c
c  NX     Input, INTEGER NX, the number of grid lines to draw along the
c         X direction.
c
c  NY     Input, INTEGER NY, the number of grid lines to draw along the
c         Y direction.
c
      integer n
      parameter (n=2)

      integer i
      integer nx
      integer ny
      real s(n)
      real smax
      real smin
      real sval
      real t(n)
      real tmax
      real tmin
      real tval

      call rldata('get','smin',smin)
      call rldata('get','smax',smax)
      call rldata('get','tmin',tmin)
      call rldata('get','tmax',tmax)

      do i=1,nx

        if(nx.ne.1)then
          sval=real(i-1)/real(nx-1)
        else
          sval=0.5
        endif

        s(1)=sval
        t(1)=tmin
        s(2)=sval
        t(2)=tmax
        call curve3(n,s,t)

      end do

      do i=1,ny

        if(ny.ne.1)then
          tval=(i-1)/real(ny-1)
        else
          tval=0.5
        endif

        s(1)=smin
        t(1)=tval
        s(2)=smax
        t(2)=tval

        call curve3(n,s,t)

      end do

      return
      end
      subroutine cgrid2(xmin,xmax,nx,ymin,ymax,ny)

c*********************************************************************72
c
c  CGRID2 draws a grid of horizontal and vertical lines.
c
c  XMIN,
c  XMAX   Input, REAL XMIN, XMAX, are the horizontal limits of the location
c         of the grid.
c
c         XMIN and XMAX should be "within" your picture.  If you are using
c         the simple DRAWCGM coordinate system, which uses 0 as the minimum
c         X value and 1 as the maximum X value, then XMIN should be 0 or
c         greater, and XMAX should be 1 or less.
c
c         However, if you have used SETWCD or SETSCL to allow a different
c         range of X values, then XMIN and XMAX may be any values in that
c         range.
c
c         Similar remarks apply to the YMIN and YMAX values.
c
c  NX     Input, INTEGER NX, the number of grid lines to draw along the
c         X direction.
c
c  YMIN,
c  YMAX   Input, REAL YMIN, YMAX, are the vertical limits of the location
c         of the grid.
c
c  NY     Input, INTEGER NY, the number of grid lines to draw along the
c         Y direction.
c
      integer n
      parameter (n=2)

      integer i
      integer nx
      integer ny
      real x(n)
      real xmax
      real xmin
      real xval
      real y(n)
      real ymax
      real ymin
      real yval

      do i=1,nx

        if(nx.ne.1)then
          xval=((nx-i)*xmin+(i-1)*xmax)/real(nx-1)
        else
          xval=0.5*(xmin+xmax)
        endif

        x(1)=xval
        y(1)=ymin
        x(2)=xval
        y(2)=ymax

        call curve2(n,x,xmax,xmin,y,ymax,ymin)

      end do

      do i=1,ny

        if(ny.ne.1)then
          yval=((ny-i)*ymin+(i-1)*ymax)/real(ny-1)
        else
          yval=0.5*(ymin+ymax)
        endif

        x(1)=xmin
        y(1)=yval
        x(2)=xmax
        y(2)=yval

        call curve2(n,x,xmax,xmin,y,ymax,ymin)

      end do

      return
      end
      subroutine curve(n,x,y)

c*********************************************************************72
c
      integer nmax
      parameter (nmax=100)

      integer n

      integer i
      real s(nmax)
      real smax
      real smin
      real t(nmax)
      real tmax
      real tmin
      real x(n)
      real xmax
      real xmin
      real y(n)
      real ymax
      real ymin
c
c  Get X and Y ranges
c
      call rldata('get','xmax',xmax)
      call rldata('get','xmin',xmin)
      call rldata('get','ymax',ymax)
      call rldata('get','ymin',ymin)
c
c  Create copy of this array which lies in (smin,smax) and (tmin,tmax)
c
      call rldata('get','smin',smin)
      call rldata('get','smax',smax)
      call rldata('get','tmin',tmin)
      call rldata('get','tmax',tmax)

      do i=1,n
        s(i)=( (xmax-x(i))*smax + (x(i)-xmin)*smin) / (smax-smin)
        t(i)=( (ymax-y(i))*tmax + (y(i)-ymin)*tmin) / (tmax-tmin)
      end do
c
c  Draw the line through the points.
c
      call gpl(n,s,t)
 
      return
      end
      subroutine size(n,x,y)

c*********************************************************************72
c
      integer n

      real x(n)
      real xmax
      real xmin
      real y(n)
      real ymax
      real ymin
c
c  Get X and Y ranges
c
      call rrange(n,x,xmax,xmin)
      call rrange(n,y,ymax,ymin)
c
c  Force these values to be used to define the size of the next plot.
c
      call rldata('set','xmax',xmax)
      call rldata('set','xmin',xmin)
      call rldata('set','ymax',ymax)
      call rldata('set','ymin',ymin)

      return
      end
      subroutine curve2(n,x,xmax,xmin,y,ymax,ymin)

c*********************************************************************72
c
      integer nmax
      parameter (nmax=100)

      integer n

      integer i
      real s(nmax)
      real smax
      real smin
      real t(nmax)
      real tmax
      real tmin
      real x(n)
      real xmax
      real xmin
      real y(n)
      real ymax
      real ymin
c
c  Create copy of this array which lies in (smin,smax) and (tmin,tmax)
c
      call rldata('get','smin',smin)
      call rldata('get','smax',smax)
      call rldata('get','tmin',tmin)
      call rldata('get','tmax',tmax)

      do i=1,n
        s(i)=( (xmax-x(i))*smax + (x(i)-xmin)*smin) / (smax-smin)
        t(i)=( (ymax-y(i))*tmax + (y(i)-ymin)*tmin) / (tmax-tmin)
      end do

      call gpm(n,s,t)
      call gpl(n,s,t)
      return
      end
      subroutine curve3(n,s,t)

c*********************************************************************72
c
      integer n

      real s(n)
      real t(n)

      call gpm(n,s,t)
      call gpl(n,s,t)

      return
      end
      subroutine gstart ( )

c*********************************************************************72
c
      integer idwork
      integer ires

      idwork=14
c
c  IRES chooses the resolution of the graphics screen.
c    A value of 2 is medium resolution.
c    A value of 3 is the highest resolution available.
c
      ires=3

      call gopks(0)
      call gopwk(idwork,0,ires)
      call gacwk(idwork)

      return
      end
      subroutine gstop ( )
c
c*********************************************************************72
c
      integer idwork

      idwork=14

      call gdawk(idwork)
      call gclwk(idwork)
      call gclks

      return
      end
      subroutine rldata(op,var,value)

c*********************************************************************72
c
c  INDATA works like a sort of COMMON block.  It stores or returns
c  the values of certain variables.  Thus, it allows two routines
c  to "communicate" without having to have data passed up and
c  down the calling tree in argument lists.
c
      character*(*) op
      logical leqi
      real smax
      real smin
      real tmax
      real tmin
      real value
      character*(*) var
      real xmax
      real xmin
      real ymax
      real ymin

      external leqi

      save smax
      save smin
      save tmax
      save tmin
      save xmax
      save xmin
      save ymax
      save ymin

      data smax /0.9/
      data smin /0.1/
      data tmax /0.9/
      data tmin /0.1/
      data xmax /1.0/
      data xmin /0.0/
      data ymax /1.0/
      data ymin /0.0/

      if(leqi(op,'set').and.leqi(var,'smax'))then
        smax=value
      elseif(leqi(op,'set').and.leqi(var,'smin'))then
        smin=value
      elseif(leqi(op,'set').and.leqi(var,'tmax'))then
        tmax=value
      elseif(leqi(op,'set').and.leqi(var,'tmin'))then
        tmin=value
      elseif(leqi(op,'set').and.leqi(var,'xmax'))then
        xmax=value
      elseif(leqi(op,'set').and.leqi(var,'xmin'))then
        xmin=value
      elseif(leqi(op,'set').and.leqi(var,'ymax'))then
        ymax=value
      elseif(leqi(op,'set').and.leqi(var,'ymin'))then
        ymin=value
c
c  Get a value.
c
      elseif(leqi(op,'get').and.leqi(var,'smax'))then
        value=smax
      elseif(leqi(op,'get').and.leqi(var,'smin'))then
        value=smin
      elseif(leqi(op,'get').and.leqi(var,'tmax'))then
        value=tmax
      elseif(leqi(op,'get').and.leqi(var,'tmin'))then
        value=tmin
      elseif(leqi(op,'get').and.leqi(var,'xmax'))then
        value=xmax
      elseif(leqi(op,'get').and.leqi(var,'xmin'))then
        value=xmin
      elseif(leqi(op,'get').and.leqi(var,'ymax'))then
        value=ymax
      elseif(leqi(op,'get').and.leqi(var,'ymin'))then
        value=ymin
      endif

      return
      end
      function leqi(strng1,strng2)

c*********************************************************************72
c
c  LEQI is a case insensitive comparison of two strings for
c  equality.  Thus, LEQI('Anjana','ANJANA') is .TRUE.
c
c  STRNG1,
c  STRNG2 Input, CHARACTER*(*) STRNG1, STRNG2, the strings to
c         compare.
c
c  LEQI   Output, LOGICAL LEQI, the result of the comparison.
c
      integer i
      integer len1
      integer len2
      integer lenc
      logical leqi
      character*1 null
      character*1 s1
      character*1 s2
      character*(*) strng1
      character*(*) strng2

      len1=len(strng1)
      len2=len(strng2)
      lenc=min(len1,len2)

      leqi=.false.
      do i=1,lenc
        s1=strng1(i:i)
        s2=strng2(i:i)
        call capchr(s1)
        call capchr(s2)
        if(s1.ne.s2)return
      end do

      null=char(0)

      do i=lenc+1,len1
        if(strng1(i:i).ne.' '.and.
     &     strng1(i:i).ne.null)return
      end do

      do i=lenc+1,len2
        if(strng2(i:i).ne.' '.and.
     &     strng2(i:i).ne.null)return
      end do

      leqi=.true.

      return
      end
      subroutine marker(mark)

c*********************************************************************72
c
c  MARKER allows the user to set the current marker type.
c
c  MARK   Input, CHARACTER*1 MARK.
c
c         MARK is a CHARACTER variable or a single quoted character,
c         which is the character to use for markers.  There are
c         five legal choices:
c
c         '.', the dot.
c         '+', the plus sign.
c         '*', the star or asterisk.
c         'o', the circle.
c         'x', the diagonal cross.
c
c         If MARK is not one of these choices, the marker type
c         will be set to the dot.
c
      integer idmark
      character*1 mark

      if(mark.eq.'.')then
        idmark=1
      elseif(mark.eq.'+')then
        idmark=2
      elseif(mark.eq.'*')then
        idmark=3
      elseif(mark.eq.'o'.or.
     &       mark.eq.'O'.or.
     &       mark.eq.'0')then
        idmark=4
      elseif(mark.eq.'x'.or.
     &       mark.eq.'X')then
        idmark=5
      else
        idmark=1
      endif

      call gsmk(idmark)

      return
      end
      subroutine points(n,x,y)

c*********************************************************************72
c
      integer nmax
      parameter (nmax=100)

      integer n

      integer i
      real s(nmax)
      real smax
      real smin
      real t(nmax)
      real tmax
      real tmin
      real x(n)
      real xmax
      real xmin
      real y(n)
      real ymax
      real ymin
c
c  Get X and Y ranges
c
      call rrange(n,x,xmax,xmin)
      call rrange(n,y,ymax,ymin)
c
c  Create copy of this array which lies in (smin,smax) and (tmin,tmax)
c
      call rldata('get','smin',smin)
      call rldata('get','smax',smax)
      call rldata('get','tmin',tmin)
      call rldata('get','tmax',tmax)

      do i=1,n
        s(i)=( (xmax-x(i))*smax + (x(i)-xmin)*smin) / (smax-smin)
        t(i)=( (ymax-y(i))*tmax + (y(i)-ymin)*tmin) / (tmax-tmin)
      end do

      call gpm(n,s,t)
      return
      end
      subroutine range(xmin,xmax,ymin,ymax)

c*********************************************************************72

      real xmax
      real xmin
      real ymax
      real ymin

      call rldata('set','xmax',xmax)
      call rldata('set','xmin',xmin)
      call rldata('set','ymax',ymax)
      call rldata('set','ymin',ymin)

      return
      end
      subroutine refpnt ( )

c*********************************************************************72
c
      integer n
      parameter (n=6)

      integer m
      parameter (m=n*n)

      integer i
      integer j
      integer k
      real x(m)
      real y(m)

      k=0
      do i=0,n-1
        do j=0,n-1
          k=k+1
          x(k)=real(i)/real(n-1)
          y(k)=real(j)/real(n-1)
        end do
      end do
c
c  Set the marker to a dot.
c
      call marker('.')
c
c  Draw a marker at each point.
c
      call points(m,x,y)

      return
      end
      subroutine rrange(nval,x,xmax,xmin)

c*********************************************************************72
c
c  RRANGE computes the range of a real array.
c
c  NVAL   Input, INTEGER NVAL, the number of entries in the array.
c
c  X      Input, REAL X(NVAL), the array.
c
c  XMAX,
c  XMIN   Output, REAL XMAX, XMIN, the largest and smallest entries
c         in the array.
c
      integer nval

      integer i
      real x(nval)
      real xmax
      real xmin

      xmax=x(1)
      xmin=x(1)
      do i=2,nval
        xmax=max(xmax,x(i))
        xmin=min(xmin,x(i))
      end do

      return
      end
      subroutine wait ( )

c*********************************************************************72
c
      integer icont
      integer idwork

      read(*,*)
      idwork=14
      icont=0
      call gclrwk(idwork,icont)
      return
      end
