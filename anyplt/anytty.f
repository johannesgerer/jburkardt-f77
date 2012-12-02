      SUBROUTINE ANYPLT(ICOM)

C***********************************************************************
C
CC TTYPLT graphics is a crude TTY plotter.
C
C
C  ANYTTY.FOR   Version 1.02  23 October 1987
C  ANYPLT/TTYPLT interface
C
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
      SAVE IXMN,IXMX,IYMN,IYMX
      SAVE IXMIN,IXMAX,IYMIN,IYMAX
      SAVE SCREEN
      SAVE XMIN,XMAX,YMIN,YMAX
      SAVE XSMIN,XSMAX,YSMIN,YSMAX
      SAVE XSMN,XSMX,YSMN,YSMX
      SAVE XSTART,YSTART
      CHARACTER CARRAY*80
      CHARACTER ISAY*1
      COMMON /ANYCOM/ IPLT1,IPLT2,IXPLT1,IXPLT2,IYPLT1,
     1                IYPLT2,MARRAY,XPLT1,XPLT2,YPLT1,YPLT2
      COMMON /ANYCHR/ CARRAY
      CHARACTER CONG*1
      CHARACTER SCREEN(24,80)*1
      XSMN=1
      XSMX=80
      YSMN=24
      YSMX=1
      IXMN=XSMN
      IXMX=XSMX
      IYMN=YSMN
      IYMX=YSMX
C
C  ICOM=0  Enable graphics
C
      IF(ICOM.EQ.0)THEN
        XSMIN=XSMN+XPLT1*(XSMX-XSMN)
        XSMAX=XSMN+XPLT2*(XSMX-XSMN)
        YSMIN=YSMN+YPLT1*(YSMX-YSMN)
        YSMAX=YSMN+YPLT2*(YSMX-YSMN)
        IXMIN=XSMIN
        IXMAX=XSMAX
        IYMIN=YSMIN
        IYMAX=YSMAX
        ENDIF
C
C  ICOM=1  Disable graphics
C
C
C  ICOM=2  Begin plot
C
      IF(ICOM.EQ.2)THEN
        DO 20 I=IYMX,IYMN
          DO 10 J=IXMN,IXMX
            SCREEN(I,J)=' '
10          CONTINUE
20        CONTINUE
        ENDIF
C
C  ICOM=3  Define plot size
C
      IF(ICOM.EQ.3)THEN
        XMIN=XPLT1
        XMAX=XMIN+XPLT2
        YMIN=YPLT1
        YMAX=YMIN+YPLT2
        ENDIF
C
C  ICOM=4  Move to point
C
      IF(ICOM.EQ.4)THEN
        XSTART=XPLT1
        YSTART=YPLT1
        ENDIF
C
C  ICOM=5  Draw to point
C
      IF(ICOM.EQ.5)THEN
        XEND=XPLT1
        YEND=YPLT1
        CALL XTOIX(IXEND,IXMAX,IXMIN,XEND,XMAX,XMIN)
        CALL XTOIX(IXSTR,IXMAX,IXMIN,XSTART,XMAX,XMIN)
        CALL XTOIX(IYEND,IYMAX,IYMIN,YEND,YMAX,YMIN)
        CALL XTOIX(IYSTR,IYMAX,IYMIN,YSTART,YMAX,YMIN)
        IF(IABS(IXEND-IXSTR).GT.IABS(IYEND-IYSTR))THEN
          IHI=IABS(IXEND-IXSTR)+1
          ISIG=+1
          IF(IXEND.LT.IXSTR)ISIG=-1
          DO 90 ICOUNT=1,IHI
            J=IXSTR+ISIG*(ICOUNT-1)
            I=IYSTR+INT(REAL((IYEND-IYSTR)*(J-IXSTR)/(IXEND-IXSTR)))
            IF(SCREEN(I,J).EQ.' ')SCREEN(I,J)='*'
90          CONTINUE
        ELSEIF(IABS(IYEND-IYSTR).GT.0)THEN
          JHI=IABS(IYEND-IYSTR)+1
          ISIG=+1
          IF(IYEND.LT.IYSTR)ISIG=-1
          DO 100 ICOUNT=1,JHI
            I=IYSTR+ISIG*(ICOUNT-1)
            J=IXSTR+INT(REAL((IXEND-IXSTR)*(I-IYSTR)/(IYEND-IYSTR)))
            IF(SCREEN(I,J).EQ.' ')SCREEN(I,J)='*'
100         CONTINUE
        ELSE
          I=IYSTR
          J=IXSTR
          IF(SCREEN(I,J).EQ.' ')SCREEN(I,J)='*'   
          ENDIF
        XSTART=XEND
        YSTART=YEND
        ENDIF
C
C  ICOM=6  Clear screen
C
      IF(ICOM.EQ.6)THEN
        DO 30 I=IYMX,IYMN
          WRITE(*,'(1X)')
30        CONTINUE
        ENDIF
C
C  ICOM=7,  Write string at position
C
      IF(ICOM.EQ.7)THEN
        CALL XTOIX(IX,IXMAX,IXMIN,XPLT1,XMAX,XMIN)
        CALL XTOIX(IY,IYMAX,IYMIN,YPLT1,YMAX,YMIN)
        DO 50 I=1,MARRAY
          IF((IX+I-1).LE.IXMAX)SCREEN(IY,IX+I-1)=CARRAY(I:I)
50        CONTINUE
        ENDIF
C
C  ICOM=8  Use virtual cursor
C
C
C  ICOM=9  End plot
C
      IF(ICOM.EQ.9)THEN
        DO 70 I=IYMAX,IYMIN
          WRITE(*,1000)(SCREEN(I,J),J=1,IXMAX)
70        CONTINUE
        WRITE(*,1040)
        READ(*,1010)
        ENDIF
C
C  ICOM=10  Ring bell
C
      IF(ICOM.EQ.10)THEN
        CONG=CHAR(7)
        WRITE(*,1030)CONG
        ENDIF
C
C  ICOM=11  Mark data
C
      IF(ICOM.EQ.11)THEN
        CALL XTOIX(IX,IXMAX,IXMIN,XPLT1,XMAX,XMIN)
        CALL XTOIX(IY,IYMAX,IYMIN,YPLT1,YMAX,YMIN)
        SCREEN(IY,IX)=CARRAY(1:1)
        ENDIF
C
C  ICOM=12  Return screen data
C
      IF(ICOM.EQ.12)THEN
        XPLT1=XSMN
        XPLT2=XSMX
        YPLT1=YSMN
        YPLT2=YSMX
        ENDIF
C
C  ICOM=13  Return version
C
      IF(ICOM.EQ.13)THEN
        CARRAY='ANYPLT - Version 1.02  23 October 1987  TTYPLT'
        ENDIF
      RETURN
1000  FORMAT(1X,80A1)
1010  FORMAT(1X)
1030  FORMAT(1X,A1)
1040  FORMAT(' Type RETURN')
      END
      SUBROUTINE XTOIX(IX,IXMAX,IXMIN,X,XMAX,XMIN)

C***********************************************************************
C
CC XTOIX converts a real value to an integer value.
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
      IX=INT(REAL(IXMAX-IXMIN)*(X-XMIN)/(XMAX-XMIN))+IXMIN
      RETURN
      END
