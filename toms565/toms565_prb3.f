      program main

c*********************************************************************72
c
c  Reference:
c
c    David Melgaard, Richard Sincovec,
c    Algorithm 565:
c    PDETWO/PSETM/GEARB: Solution of systems of two-dimensional
c    nonlinear partial differential equations,
c    Volume 7, Number 1, March 1981, pages 126-135.
c
      REAL  H,S,X,ERR1,DX,Y,TOUT,REALN1,DY,U1,HUSED,
     *      EXP,R,T0,EPS,ERMAX,ABS,WORK,DL,DLI,DEV,
     *      U2,U3,ERR2,ERR3,REALN2,REALN3
      INTEGER IX,NPDE,NSTEP,NX,NFE,NODE,MF,NY,
     *        NJE,NQUSED,IY,INDEX,I,IWORK,KODE,IK
      COMMON /GEAR3/ HUSED,NQUSED,NSTEP,NFE,NJE
      COMMON /PROB/ DL,DLI,KODE
      DIMENSION U1(10,10),ERR1(10,10),REALN1(10,10)
      DIMENSION U2(5,31),ERR2(5,31),REALN2(5,31)
      DIMENSION U3(2,5,5),ERR3(2,5,5),REALN3(2,5,5)
      DIMENSION WORK(4186),IWORK(155),X(10),Y(31)
      EQUIVALENCE (U1(1,1),U3(1,1,1)),(ERR1(1,1),ERR3(1,1,1)),
     *            (REALN1(1,1),REALN3(1,1,1))
      EQUIVALENCE (U2(1,1),U3(1,1,1)),(ERR2(1,1),ERR3(1,1,1)),
     *            (REALN2(1,1),REALN3(1,1,1))

      save

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TOMS565_PRB3:'
      write ( *, '(a)' ) '  Test the TOMS565 library.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Example 3, Coupled PDE System'
      write ( *, '(a)' ) ' '
C
C  DEFINE THE PROBLEM PARAMETERS.  NPDE,NX, AND NY PRIMARILY DETERMINE
C  THE DIMENSIONS FOR THE ARRAYS IN PDETWO AND THE MODIFIED GEARB.
C
      NX=5
      NY=5
      NPDE=2
      NODE=NX*NY*NPDE
      MF=22
      INDEX=1
      T0=0.0
      H=0.1E-06
      EPS=0.1E-03
      DX=1.0/(FLOAT(NX)-1.0)
      DO 520 IX=1,NX
        X(IX)=FLOAT(IX)*DX-DX
  520   Y(IX)=X(IX)
      IWORK(1) = NPDE
      IWORK(2) = NX
      IWORK(3) = NY
      IWORK(4) = 5
      IWORK(5) = 2308
      IWORK(6) = 50
C
C  DEFINE THE INITIAL CONDITIONS
C
      DO 560 IY=1,NY
        DO 560 IX=1,NX
          U3(1,IX,IY)=X(IX)+Y(IY)
  560     U3(2,IX,IY)=2.0*X(IX)+3.0*Y(IY)
      WRITE (6,570) NODE,T0,H,EPS,MF,((U3(1,IX,IY),IX=1,NX),IY=1,NY),
     *             ((U3(2,IX,IY),IX=1,NX),IY=1,NY)
  570 FORMAT (6H NODE ,I3,4H T0 ,F9.2,3H H ,E8.1,5H EPS ,E8.1,4H MF ,I2
     *        //19H INITIAL U1 VALUES  / 5(/3H    ,5E11.4)
     *        //19H INITIAL U2 VALUES  / 5(/3H    ,5E11.4))
C
C  SET UP THE LOOP FOR CALLING THE INTEGRATOR AT DIFFERENT TOUT VALUES
C
      TOUT=.25
      DO 660 I=1,4
      WRITE (6,600) TOUT
  600 FORMAT (//5H TOUT ,E15.6)
C
C  CALCULATE THE ACTUAL SOLUTION
C
      DO 610 IY=1,NY
        DO 610 IX=1,NX
          REALN3(2,IX,IY)=TOUT+2.0*X(IX)+3.0*Y(IY)
  610     REALN3(1,IX,IY)=TOUT +X(IX)+Y(IY)
C
C  CALL THE INTEGRATOR
C
      CALL DRIVEP (NODE,T0,H,U3,TOUT,EPS,MF,INDEX,WORK,IWORK,X,Y)
C
C  CHECK ERROR RETURN
C
      IF (INDEX .NE. 0) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'TOMS565_PRB3 - Fatal error!'
        write ( *, '(a,i4)' ) '  DRIVEP error return INDEX = ', index
        stop
      end if
C
C  OUTPUT THE RESULTS
C
      WRITE (6,620) HUSED,NQUSED,NSTEP,NFE,NJE
  620 FORMAT (/7H HUSED ,E10.4,7H ORDER ,I3,7H NSTEP ,I6,5H NFE ,I5,
     *        5H NJE ,I5)
      WRITE (6,630) ((U3(1,IX,IY),IX=1,NX),IY=1,NY)
  630 FORMAT (/11H U1 VALUES //(3H    ,5E11.4))
      WRITE (6,635) ((U3(2,IX,IY),IX=1,NX),IY=1,NY)
  635 FORMAT (/11H U2 VALUES //(3H    ,5E11.4))
C
C  DETERMINE THE ERROR
C
      DO 640 IY=1,NY
         DO 640 IX=1,NX
          ERR3(1,IX,IY)=U3(1,IX,IY)-REALN3(1,IX,IY)
          ERR3(2,IX,IY)=U3(2,IX,IY)-REALN3(2,IX,IY)
  640     CONTINUE
      WRITE (6,650) ((ERR3(1,IX,IY),IX=1,NX),IY=1,NY)
  650 FORMAT (/9H ERROR U1 //(3H    ,5E11.4))
      WRITE (6,655) ((ERR3(2,IX,IY),IX=1,NX),IY=1,NY)
  655 FORMAT (/9H ERROR U2 //(3H    ,5E11.4))
        TOUT=TOUT+.25
  660 CONTINUE
c
c  Terminate.
c
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TOMS565_PRB3'
      write ( *, '(a)' ) '  Normal end of execution.'

      STOP
      END
      SUBROUTINE BNDRYV (T,X,Y,U,AV,BV,CV,NPDE)

c*********************************************************************72
C
C  DEFINE THE VERTICAL BOUNDARY CONDITIONS
C
      REAL  T,U,X,Y,BV,AV,CV
      INTEGER NPDE
      COMMON /PROB/ DL,DLI,KODE
      DIMENSION U(NPDE),AV(NPDE),BV(NPDE),CV(NPDE)

      save

      IF (X .NE. 0.0) GO TO 310
      AV(1)=1.0
      BV(1)=0.0
      CV(1)=T+Y
      AV(2)=1.0
      BV(2)=0.0
      CV(2)=T+3.0*Y
      GO TO 400
  310 CONTINUE
      AV(1)=1.0
      BV(1)=0.0
      CV(1)=T+1.0+Y
      AV(2)=0.0
      BV(2)=1.0
      CV(2)=2.0
  400 CONTINUE
      RETURN
      END
      SUBROUTINE BNDRYH (T,X,Y,U,AH,BH,CH,NPDE)

c*********************************************************************72
C
C  DEFINE THE HORIZONTAL BOUNDARY CONDITIONS
C
      REAL  T,U,X,Y,BH,AH,CH
      INTEGER NPDE
      COMMON /PROB/ DL,DLI,KODE
      DIMENSION U(NPDE),AH(NPDE),BH(NPDE),CH(NPDE)

      save

      IF (Y .NE. 0.0) GO TO 310
      AH(1)=0.0
      BH(1)=1.0
      CH(1)=1.0
      AH(2)=1.0
      BH(2)=1.0
      CH(2)=T+2.0*X+3.0
      GO TO 400
  310 CONTINUE
      AH(1)=0.0
      BH(1)=U(2)
      CH(1)=T+2.0*X+3.0
      AH(2)=1.0
      BH(2)=0.0
      CH(2)=T+2.0*X+3.0
  400 CONTINUE

      RETURN
      END
      SUBROUTINE DIFFH (T,X,Y,U,DH,NPDE)

c*********************************************************************72
C
C  DEFINE THE HORIZONTAL DIFFUSION COEFFICIENTS
C
      REAL  T,U,X,Y,DH
      INTEGER NPDE
      COMMON /PROB/ DL,DLI,KODE
      DIMENSION U(NPDE),DH(NPDE,NPDE)

      save

      DH(1,1)=0.0
      DH(1,2)=0.0
      DH(2,1)=U(2)
      DH(2,2)=U(1)

      RETURN
      END
      SUBROUTINE DIFFV (T,X,Y,U,DV,NPDE)

c*********************************************************************72
C
C  DEFINE THE VERTICAL DIFFUSION COEFFICIENTS
C
      REAL  T,U,X,Y,DV
      INTEGER NPDE
      COMMON /PROB/ DL,DLI,KODE
      DIMENSION U(NPDE),DV(NPDE,NPDE)

      save

      DV(1,1)=U(2)
      DV(1,2)=U(1)
      DV(2,1)=0.0
      DV(2,2)=0.0

      RETURN
      END
      SUBROUTINE F(T,X,Y,U,UX,UY,DUXX,DUYY,DUDT,NPDE)

c*********************************************************************72
C
C  DEFINE THE PDE
C
      REAL  T,U,X,Y,UX,UY,DUXX,DUYY,DUDT,EXP
      INTEGER NPDE
      COMMON /PROB/ DL,DLI,KODE
      DIMENSION U(NPDE),UX(NPDE),UY(NPDE),DUXX(NPDE,NPDE),
     * DUYY(NPDE,NPDE),DUDT(NPDE)

      save

      DUDT(1)=U(2)*(DUYY(1,1)+DUYY(1,2))-3.0*U(2)*UX(2)+4.0*UY(1)-UY(2)
      DUDT(2)=DUXX(2,1)+DUXX(2,2)-UY(2)

      RETURN
      END
