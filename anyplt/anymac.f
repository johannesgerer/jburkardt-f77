c  ANYMAC.F   Version 1.08  09 October 1990
c  ANYPLT/Macintosh interface
c
c  ANYPLT is a subroutine which provides a simple, standard interface
c  between FORTRAN programs and various output devices.  To run a
c  program which calls ANYPLT on a different machine, the program
c  is not modified in any way, but a different version of the ANYPLT
c  program is provided.  Currently, the following versions are available:
c
c  ANYATT - AT&T PC6300 graphics (640 by 400).  Requires ATTPLT.ASM.
c  ANYBUG - Simple debugging output to a file.
c  ANYCAL - CALCOMP file output.  Available on many mainframes.
c  ANYIBM - IBM PC hi resolution (640 by 200).  Requires IBMPLT.ASM.
c  ANYMAC - Macintosh graphics.  (512 by 342) Requires auxilliary TOOLBX.SUB.
c  ANYNCR - NCAR graphics package.
c  ANYNUL - Does nothing.
c  ANYP10 - PLOT10 interactive graphics. (1024 by 768)
c  ANYTTY - Simple 'typewriter' graphics (80 by 24)
c
c  The personal computer routines generally require an auxilliary
c  set of routines written in assembly language.
c
c  ANYMAC was written by John Burkardt
c
c  The symbolic output of characters, numbers and other printable
c  characters was made possible by adaptation of a routine written
c  by Bill Furey of the University of Pittsburgh, Department of
c  Crystallography.
c
      subroutine anyplt(icom)

c***********************************************************************
c
cc ANYPLT is a generic graphic interface.
c
c  This version is for Macintosh graphics.
c
c  When linking, the compiled file 'toolbx.sub' must also be included.
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
      integer*4 lineto,moveto
c     parameter (lineto=z'89109000')
c     parameter (moveto=z'89309000')
c
      save ifont
      save ipoint
      save ixmn,ixmx,iymn,iymx
      save ixmin,ixmax,iymin,iymax
      save rotate
      save xmin,xmax,ymin,ymax
      save xold
      save xsmin,xsmax,ysmin,ysmax
      save xsmn,xsmx,ysmn,ysmx
      save xstart,ystart
      save yold
c
      character carray*80
      character cong*1
      integer*4 grafptr
      integer   ifont(1617)
      integer   ipoint(95)
      logical   rotate
c
      common /anycom/ iplt1,iplt2,ixplt1,ixplt2,iyplt1,
     &                iyplt2,marray,xplt1,xplt2,yplt1,yplt2
      common /anychr/ carray
c
c  Pointer array into IFONT
c
      data (ipoint(i),i=1,95) /
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
      data (ifont(i),i=   1, 396)/
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
      data (ifont(i),i= 397, 792)/
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
      data (ifont(i),i= 793,1188)/
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
      data (ifont(i),i=1189,1584)/
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
      data (ifont(i),i=1585,1617)/
     & 0, 1, 9,23, 9, 7, 0, 1, 7,23,11,23, 9,17,11,15, 9,13,11, 7, 7, 7,
     & 0, 1, 5,21, 7,23,15,21,17,23, 0/
c
c  ICOM=0  Enable graphics.  For interactive graphics, clear the screen.
c  For graphics which create a file, open the file.
c  Also, user will input the portion of the (0,1) by (0,1) output
c  screen that is to be used.  At this point, the effective screen
c  (IXMIN,IXMAX) by (IYMIN,IYMAX) must be computed.
c
c  Note that the whole Macintosh screen is NOT available to Microsoft FORTRAN.
c  Note also the inversion that causes YSMN to be greater than YSMX!
c
      if(icom.eq.0)then
        xsmn=116
        xsmx=396
        ysmn=311
        ysmx=31
        ixmn=xsmn
        ixmx=xsmx
        iymn=ysmn
        iymx=ysmx
        if(xplt1.eq.xplt2)then
          write(*,*)'anymac - warning!  icom=0 command is asking for'
          write(*,*)'         zero screen width (xplt1=xplt2).'
          xplt1=0.0
          xplt2=1.0
          endif
        xsmin=xsmn+xplt1*(xsmx-xsmn)
        xsmax=xsmn+xplt2*(xsmx-xsmn)
        if(yplt1.eq.yplt2)then
          write(*,*)'anymac - warning!  icom=0 command is asking for'
          write(*,*)'         zero screen height (yplt1=yplt2).'
          yplt1=0.0
          yplt2=1.0
          endif
        ysmin=ysmn+yplt1*(ysmx-ysmn)
        ysmax=ysmn+yplt2*(ysmx-ysmn)
        ixmin=xsmin
        ixmax=xsmax
        iymin=ysmin
        iymax=ysmax
        xmin=0.0
        xmax=1.0
        ymin=0.0
        ymax=1.0
        endif
c
c  ICOM=1  Disable graphics.  For interactive graphics, clear the screen.
c  For graphics which create a file, close the file.
c  (Doesn't really apply for Macintosh)
c
c  ICOM=2  Begin plot
c
      if(icom.eq.2)then
        xmin=0.0
        xmax=1.0
        ymin=0.0
        ymax=1.0
        do 33 i=1,48
          write(*,'(1x)')
33        continue
        endif
c
c  ICOM=3  Define user coordinate system (XMIN,XMAX), (YMIN,YMAX).
c
      if(icom.eq.3)then
        xmin=xplt1
        xmax=xmin+xplt2
        ymin=yplt1
        ymax=ymin+yplt2
        endif
c
c  ICOM=4  Move to point.  Only purpose is to begin a line.
c  Input is X, Y in user coordinates.
c
      if(icom.eq.4)then
        xold=xplt1
        yold=yplt1
        endif
c
c  ICOM=5  Draw to point.  Only called after a previous move or draw.
c  Input is X, Y in user coordinates.
c
      if(icom.eq.5)then
        call drwlin(xold,yold,xplt1,yplt1,xmin,xmax,ymin,ymax,
     &  ixmin,ixmax,iymin,iymax)
        xold=xplt1
        yold=yplt1
        endif
c
c  ICOM=6  Clear screen.  Must be a better way on the Macintosh, but I
c  don't know it yet.
c
      if(icom.eq.6)then
        do 330 i=1,48
          write(*,'(1x)')
330       continue
          endif
c
c  ICOM=7,  Write string at position.
c  Variable height and angle should be supported.
c  Note that for this command, screen coordinates are used.
c  Thus a width of 0.1 is intended to mean 1/10 of screen size.
c
      if(icom.eq.7)then
c
c  Set scale factor for character height
c
        csize=xplt2
        angle=yplt2
        scl2=csize*.0625
c
c  Set starting point for line of text (lower left corner of first
c  character) and origin about which rotation is performed.
c
        xb=xplt1
        xrot=xplt1
        yb=yplt1
        yrot=yplt1
        rotate=.false.
c
c  Get trig functions if rotation required, converting from
c  degrees to radians.
c
        if(angle.ne.0.0)then
          ca=cos(angle*.017453)
          sa=sin(angle*.017453)
          rotate=.true.
          endif
c
c  Loop over all characters in string
c
        do 30 icr=1,marray
c
c  Get ASCII code for character and shift by 31 so first printable
c  character becomes code 1
c
          iascii=ichar(carray(icr:icr))-31
c
c  Replace any non-printable characters with blanks
c
          if((iascii.lt.1).or.(iascii.gt.95))iascii=1
c
c  Get pointer to this character in font table
c
          ip=ipoint(iascii)
c
c  Get number of "vectors" required to draw character.
c  Here "vectors" means number of times pen is lowered, not the
c  number of pen strokes. (=1 for blanks, due to way algorithm is coded)
c
          nvec=ifont(ip)
c
c  Loop over all required pen movements
c
          do 20 iv=1,nvec
            ipen=3
            ip=ip+1
15          continue
            if(ifont(ip).eq.0)go to 20
            x=xb+(ifont(ip)-1)*scl2
            y=yb+(ifont(ip+1)-7)*scl2
c
c  Apply rotation if necessary
c
            if(rotate)then
              xt=x-xrot
              yt=y-yrot
              x=ca*xt-sa*yt+xrot
              y=sa*xt+ca*yt+yrot
              endif
c
c  Plot the pen stroke
c
           if(ipen.eq.3)then
              call xtoix(ix2,ixmax,ixmin,x,1.0,0.0)
              call xtoix(iy2,iymax,iymin,y,1.0,0.0)
              ix=ix2
              iy=iy2
              call toolbx(moveto,ix,iy)
            else
              ix1=ix2
              iy1=iy2
              call xtoix(ix2,ixmax,ixmin,x,1.0,0.0)
              call xtoix(iy2,iymax,iymin,y,1.0,0.0)
              ix=ix2
              iy=iy2
              call toolbx(lineto,ix,iy)
              endif
            ipen=2
            ip=ip+2
            go to 15
20          continue
c
c  Advance base to compensate for character just drawn
c
          xb=xb+csize
30        continue
        endif
c
c  ICOM=8  Use virtual cursor.  Not implemented.
c
c  ICOM=9  End plot
c
      if(icom.eq.9)then
        write(*,'(1x)')
        read(*,'(1x)')
        do 70 i=1,48
          write(*,'(1x)')
70        continue
        endif
c
c  ICOM=10  Ring bell
c
      if(icom.eq.10)then
        cong=char(7)
        write(*,'(1x,a1)')cong
        endif
c
c  ICOM=11  Mark data with a set of strokes that are like
c  a *, a +, or an X.  If a '.' is requested, actually try to
c  draw a point (pixel) if possible.
c
      if(icom.eq.11)then
        if(carray(1:1).eq.' ')then
          return
        elseif(carray(1:1).eq.'.')then
          call xtoix(ix1,ixmax,ixmin,xplt1,xmax,xmin)
          call xtoix(iy1,iymax,iymin,yplt1,ymax,ymin)
          ix=ix1
          iy=iy1
          call toolbx(moveto,ix,iy)
          call toolbx(lineto,ix,iy)
        else
          delt=0.5*xlen*14.0/real(ixlen)
          x1=xplt1+delt
          x2=xplt1-delt
          y1=yplt1+delt
          y2=yplt1-delt
          call xtoix(ix1,ixmax,ixmin,x1,xmax,xmin)
          call xtoix(ix2,ixmax,ixmin,x2,xmax,xmin)
          call xtoix(ix3,ixmax,ixmin,x3,xmax,xmin)
          call xtoix(iy1,iymax,iymin,y1,ymax,ymin)
          call xtoix(iy2,iymax,iymin,y2,ymax,ymin)
          call xtoix(iy3,iymax,iymin,y3,ymax,ymin)
          if(carray(1:1).ne.'+')then
            ix=ix1
            iy=iy2
            call toolbx(moveto,ix,iy)
            ix=ix2
            iy=iy1
            call toolbx(lineto,ix,iy)
            ix=ix1
            iy=iy1
            call toolbx(moveto,ix,iy)
            ix=ix2
            iy=iy2
            call toolbx(lineto,ix,iy)
            endif
          if(carray(1:1).ne.'x'.and.carray(1:1).ne.'x')then
            ix=ix1
            iy=iy3
            call toolbx(moveto,ix,iy)
            ix=ix2
            iy=iy3
            call toolbx(lineto,ix,iy)
            ix=ix3
            iy=iy1
            call toolbx(moveto,ix,iy)
            ix=ix3
            iy=iy2
            call toolbx(lineto,ix,iy)
            ix2=ix3
            endif
          endif
        endif
c
c  ICOM=12, Return screen X, Y maximum and minimum in pixels
c  or other 'natural' coordinates.
c
      if(icom.eq.12)then
        xplt1=xsmn
        xplt2=xsmx
        yplt1=ysmn
        yplt2=ysmx
        endif
c
c  ICOM=13  Return version number and date of code.
c
      if(icom.eq.13)then
        carray='anyplt - version 1.08  09 october 1990  macintosh'
        endif
c
c  ICOM=14, draw an arrow of given screen length and screen angle.
c  Even if the X and Y axes are of vastly different scales, this
c  allows one to draw vectors of a predictable direction.
c
      if(icom.eq.14)then
        x1=xplt1
        y1=yplt1
        angle=xplt2
        alen=yplt2
        x2=x1+(xmax-xmin)*alen*cos(angle)
        y2=y1+(ymax-ymin)*alen*sin(angle)
        call drwlin(x1,y1,x2,y2,xmin,xmax,ymin,ymax,
     &  ixmin,ixmax,iymin,iymax)
        endif
      return
      end
      subroutine drwlin(xold,yold,xplt1,yplt1,xmin,xmax,ymin,ymax,
     &ixmin,ixmax,iymin,iymax)
c
      integer*4 ix
      integer*4 iy
      integer*4 lineto,moveto
c
c     parameter (lineto=z'89109000')
c     parameter (moveto=z'89309000')
c
      call clip(xold,yold,xplt1,yplt1,xc,yc,xd,yd,idraw,
     &xmin,xmax,ymin,ymax)
c
      if(idraw.eq.1)then
        call xtoix(ix1,ixmax,ixmin,xc,xmax,xmin)
        call xtoix(iy1,iymax,iymin,yc,ymax,ymin)
        call xtoix(ix2,ixmax,ixmin,xd,xmax,xmin)
        call xtoix(iy2,iymax,iymin,yd,ymax,ymin)
c
        ix=ix1
        iy=iy1
        call toolbx(moveto,ix,iy)
        ix=ix2
        iy=iy2
        call toolbx(lineto,ix,iy)
        endif
      return
      end
      subroutine xtoix(ix,ixmax,ixmin,x,xmax,xmin)
c
c***********************************************************************
c
      ix=int(real(ixmax-ixmin)*(x-xmin)/(xmax-xmin))+ixmin
      if(ix.lt.ixmin.and.ix.lt.ixmax)ix=min(ixmin,ixmax)
      if(ix.gt.ixmax.and.ix.gt.ixmax)ix=max(ixmin,ixmax)
      return
      end
      subroutine clip(xa,ya,xb,yb,xc,yc,xd,yd,idraw,x0,x1,y0,y1)
c
      dimension xval(2),yval(2)
c
      call clip1(xa,ya,xb,yb,xval,yval,next,x0,x1,y0,y1)
      if(next.lt.2)then
        idraw=0
      else
        idraw=1
        endif
      xc=xval(1)
      yc=yval(1)
      xd=xval(2)
      yd=yval(2)
      return
      end
      subroutine clip1(xa,ya,xb,yb,xval,yval,next,x0,x1,y0,y1)
c
      dimension xval(2),yval(2)
c
c  Check to see if both points are interior to the box, and hence
c  nothing need be done.
c
      x00=min(x0,x1)
      x11=max(x0,x1)
      y00=min(y0,y1)
      y11=max(y0,y1)
c
      next=0
      if(
     &  (x00.le.xa.and.xa.le.x11).and.
     &  (y00.le.ya.and.ya.le.y11))then
        next=next+1
        xval(next)=xa
        yval(next)=ya
        endif
      if(
     &  (x00.le.xb.and.xb.le.x11).and.
     &  (y00.le.yb.and.yb.le.y11))then
        next=next+1
        xval(next)=xb
        yval(next)=yb
        endif
      if(next.eq.2)return
c
      call clip2(xa,ya,xb,yb,x00,y00,y11,y,xval,yval,next)
      if(next.eq.2)return
      call clip2(xa,ya,xb,yb,x11,y00,y11,y,xval,yval,next)
      if(next.eq.2)return
      call clip2(ya,xa,yb,xb,y00,x00,x11,x,yval,xval,next)
      if(next.eq.2)return
      call clip2(ya,xa,yb,xb,y11,x00,x11,x,yval,xval,next)
      return
      end
      subroutine clip2(xa,ya,xb,yb,x00,y00,y11,y,xval,yval,next)
c
      dimension xval(*),yval(*)
c
      if(xb.eq.xa)return
c
      t=(x00-xa)/(xb-xa)
      if(t.ge.0.0.and.t.le.1.0)then
        y=(1.0-t)*ya+t*yb
        if(y00.le.y.and.y.le.y11)then
          next=next+1
          xval(next)=x00
          yval(next)=y
          endif
        endif
      return
      end
