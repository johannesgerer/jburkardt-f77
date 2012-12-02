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
      write ( *, '(a)' ) 'TOMS565_PRB1:'
      write ( *, '(a)' ) '  Test the TOMS565 library.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Example 1, Elliptic PDE'
      write ( *, '(a)' ) ' '
C
C  DEFINE THE PROBLEM PARAMETERS.
C
      NX=10
      NY=10
      NPDE=1
      NODE=NPDE*NX*NY
      MF=22
      INDEX=1
      T0=0.0
      H=0.1E-06
      EPS=0.1E-06
      DX = 1.0/(FLOAT(NX)-1.0)
      DY = 1.0/(FLOAT(NY)-1.0)
      DO IX=1,NX
        Y(IX)=FLOAT(IX)*DY-DY
        X(IX)=FLOAT(IX)*DX-DX
      end do
      IWORK(1) = NPDE
      IWORK(2) = NX
      IWORK(3) = NY
      IWORK(4) = 5
      IWORK(5) = 4186
      IWORK(6) = 100
C
C  CALCULATE THE ACTUAL SOLUTION
C
      DO 130 IY = 1,NY
        DO 130 IX = 1,NX
          R = X(IX)
          S = Y(IY)
  130     REALN1(IX,IY)=3.0*EXP(R)*EXP(S)*(R-R*R)*(S-S*S)
      WRITE (6,140) ((REALN1(IX,IY),IX=1,NX),IY=1,NY)
  140 FORMAT (/23H STEADY STATE SOLUTION /10(/3H    ,10E11.4))
C
C  DEFINE THE INITIAL CONDITION
C
      DO 160 IY = 1,NY
        DO 160 IX = 1,NX
  160    U1(IX,IY)=0.0
      WRITE (6,170) NODE,T0,H,EPS,MF,U1
  170 FORMAT (//6H NODE ,I3,4H T0 ,F9.2,3H H ,E8.1,5H EPS ,E8.1,4H MF ,
     *        I2 //16H INITIAL VALUES / 10(/3H    ,10E11.2))
C
C  SET UP THE LOOP FOR CALLING THE INTEGRATOR AT DIFFERENT TOUT VALUES
C
      TOUT=.001

      DO 260 I=1,5

        WRITE (6,200) TOUT
  200   FORMAT (// 5H TOUT,E15.6)
C
C  CALL THE INTEGRATOR
C
      CALL DRIVEP (NODE,T0,H,U1,TOUT,EPS,MF,INDEX,WORK,IWORK,X,Y)
C
C  CHECK ERROR RETURN
C
      IF (INDEX .NE. 0) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'TOMS565_PRB1 - Fatal error!'
        write ( *, '(a,i4)' ) '  DRIVEP error return INDEX = ', index
        stop
      end if
C
C  OUTPUT THE RESULTS
C
      WRITE (6,220) HUSED,NQUSED,NSTEP,NFE,NJE
  220 FORMAT (/7H HUSED ,E10.4,7H ORDER ,I3,7H NSTEP ,I6,5H NFE ,I5,
     *        5H NJE ,I5)
      WRITE (6,230) ((U1(IX,IY),IX=1,NX),IY=1,NY)
  230 FORMAT ( /9H U VALUES  //(3H    ,10E11.4))
C
C  DETERMINE THE ERROR
C
      ERMAX = 0.0
      DO IY = 1,NY
         DO IX = 1,NX
          ERR1(IX,IY) = U1(IX,IY) - REALN1(IX,IY)
          IF (ABS(ERMAX) .LT. ABS(ERR1(IX,IY))) ERMAX=ERR1(IX,IY)
        end do
      end do

      WRITE (6,250) ERMAX
  250 FORMAT (/51H MAXIMUM DIFFERENCE FROM THE STEADY STATE SOLUTION //
     *        (3H    ,E11.4))
        TOUT=TOUT*10.

  260 CONTINUE
c
c  Terminate.
c
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TOMS565_PRB1'
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

      AV(1) = 1.0
      BV(1) = 0.0
      CV(1) = 0.0

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

      AH(1) = 1.0
      BH(1) = 0.0
      CH(1) = 0.0

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

      DH(1,1) = 1.0

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

      DV(1,1) = 1.0

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

      DUDT(1) = DUXX(1,1) + DUYY(1,1) - 6.0*X*Y*EXP(X)*EXP(Y)*
     * (X*Y+X+Y-3.0)

      RETURN
      END
