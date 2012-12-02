      program main

c*********************************************************************72
c
cc MAIN is the main program for TOMS581_PRB.
c
c  Discussion:
c
C
C       DRIVER FOR HYBSVD, A ROUTINE FOR COMPUTING THE SVD OF A MATRIX.
C       RESULTS ARE COMPARED TO THOSE PRODUCED BY SSVDC FROM LINPACK,
C       SVD AND MINFIT FROM EISPACK.
C
C       MATRIX TESTED:  A(I,J) = 1/(I+J),  1ST RHS(I) = 1
C                                          2ND RHS(I) = I
C       (A AND A' )
C
C       PROGRAMMED BY
C           TONY CHAN
C           COMPUTER SCIENCE DEPT., YALE UNIV.,
C           BOX 2158, YALE STATION,
C           NEW HAVEN, CT 06520.
C
C       LAST REVISION.. JANUARY, 1982.
C
      REAL AE(15,15), AH(10,10), AL(20,20), WORK(10), E(10)
      REAL UH(11,10), VH(12,20), Z(13,20), BH(14,5), WH(10)
      REAL UE(15,10), VE(15,10), BE(15,5), WE(10)
      REAL UL(18,10), VL(19,10), WL(10)
      LOGICAL MATU, MATV
      INTEGER NUH, NVH, NZ, NBH
      INTEGER NUE, NVE, NBE
      INTEGER NUL, NVL
      INTEGER NAE, NAH, NAL, M, N, IRHS, IERR, JOB
      DATA NAH /10/, NUH /11/, NVH /12/, NZ /13/, NBH /14/
      DATA NUE /15/, NVE /15/, NBE /15/, NAE /15/
      DATA NUL /18/, NVL /19/, NAL /20/

      call timestamp ( )
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TOMS581_PRB'
      write ( *, '(a)' ) '  FORTRAN77 version'
      write ( *, '(a)' ) '  Test the TOMS581 library.'
C
C  FIRST CASE... A AND RHS B.
C 
      M = 8
      N = 4
      IRHS = 2

      DO IU=1,2
        DO IV=1,2

          IF (IU.EQ.2) MATU = .TRUE.
          IF (IU.EQ.1) MATU = .FALSE.
          IF (IV.EQ.2) MATV = .TRUE.
          IF (IV.EQ.1) MATV = .FALSE.
C
C  JOB PARAMETER NEEDED BY SSVDC.
c
          JOB = (IU-1)*20 + IV - 1
C
C  SET UP MATRIX
C
          DO I=1,M
            BH(I,1) = 1.0
            BH(I,2) = FLOAT(I)
            BE(I,1) = 1.0
            BE(I,2) = FLOAT(I)
            DO J=1,N
              AH(I,J) = 1.0 / FLOAT(I+J)
              AE(I,J) = AH(I,J)
              AL(I,J) = AH(I,J)
            end do
          end do

          write ( *,99999)
99999     FORMAT (/4H A ./)
          DO I=1,M
            write ( *,99998) (AH(I,J),J=1,N)
99998       FORMAT (5E15.7)
          end do

          write ( *,99997)
99997     FORMAT (//45H ******************* HYBSVD *****************//)
          CALL HYBSVD(NAH, NUH, NVH, NZ, NBH, M, N, AH, WH, MATU, UH,
     &     MATV, VH, Z, BH, IRHS, IERR, WORK)
          CALL PRINTS(WH, UH, VH, IERR, BH, NUH, NVH, NBH, MATU, MATV,
     &     M, N, IRHS)

          write ( *,99996)
99996     FORMAT (//45H ******************* EISPAK *****************//)
          CALL SVD(NUE, M, N, AE, WE, MATU, UE, MATV, VE, IERR, WORK)
          CALL MINFIT(NUE, M, N, AE, WE, IRHS, BE, IERR, WORK)
          CALL PRINTS(WE, UE, VE, IERR, BE, NUE, NVE, NBE, MATU, MATV,
     &     M, N, IRHS)

          write ( *,99995)
99995     FORMAT (//45H ******************* SSVDC  *****************//)
          CALL SSVDC(AL, NAL, M, N, WL, E, UL, NUL, VL, NVL, WORK, JOB,
     &     IERR)
          CALL PRINTS(WL, UL, VL, IERR, E, NUL, NVL, NBH, MATU, MATV,
     &     M, N, IRHS)

          CALL DIFF(WE, WH, WL, UE, UH, UL, VE, VH, VL, BH, BE, NUE,
     &     NUH, NUL, NVE, NVH, NVL, NBH, NBE, M, N, IRHS, MATU,
     &     MATV)

        end do
      end do
C
C  SECOND CASE...  A'.
C
      M = 4
      N = 8
      IRHS = 2

      DO IU=1,2
        DO IV=1,2

          IF (IU.EQ.2) MATU = .TRUE.
          IF (IU.EQ.1) MATU = .FALSE.
          IF (IV.EQ.2) MATV = .TRUE.
          IF (IV.EQ.1) MATV = .FALSE.
C
C  JOB PARAMETER NEEDED BY SSVDC.
c
          JOB = (IU-1)*20 + IV - 1
C
C  SET UP MATRIX
C
          DO I=1,M
            BH(I,1) = 1.0
            BH(I,2) = FLOAT(I)
            BE(I,1) = 1.0
            BE(I,2) = FLOAT(I)
            DO J=1,N
              AH(I,J) = 1.0/FLOAT(I+J)
              AE(I,J) = AH(I,J)
              AL(I,J) = AH(I,J)
            end do
          end do

          write ( *,99999)
          DO I=1,M
            write ( *,99998) (AH(I,J),J=1,N)
          end do

          write ( *,99997)
          CALL HYBSVD(NAH, NUH, NVH, NZ, NBH, M, N, AH, WH, MATU, UH,
     &     MATV, VH, Z, BH, IRHS, IERR, WORK)
          CALL PRINTS(WH, UH, VH, IERR, BH, NUH, NVH, NBH, MATU, MATV,
     &     M, N, IRHS)

          write ( *,99996)
          CALL SVD(NUE, M, N, AE, WE, MATU, UE, MATV, VE, IERR, WORK)
          CALL MINFIT(NUE, M, N, AE, WE, IRHS, BE, IERR, WORK)
          CALL PRINTS(WE, UE, VE, IERR, BE, NUE, NVE, NBE, MATU, MATV,
     &     M, N, IRHS)

          write ( *,99995)
          CALL SSVDC(AL, NAL, M, N, WL, E, UL, NUL, VL, NVL, WORK, JOB,
     &     IERR)
          CALL PRINTS(WL, UL, VL, IERR, E, NUL, NVL, NBH, MATU, MATV,
     &     M, N, IRHS)
          CALL DIFF(WE, WH, WL, UE, UH, UL, VE, VH, VL, BH, BE, NUE,
     &     NUH, NUL, NVE, NVH, NVL, NBH, NBE, M, N, IRHS, MATU,
     &     MATV)

        end do
      end do
c
c  Terminate.
c
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TOMS581_PRB'
      write ( *, '(a)' ) '  Normal end of execution.'
      write ( *, '(a)' ) ' '
      call timestamp ( )

      STOP
      END
      SUBROUTINE PRINTS ( W, U, V, IERR, B, NAU, NV, NB, MATU, MATV, M,
     & N, IRHS)

c*********************************************************************72
c
cc PRINTS prints out the singular values and the matrices U, V, B.
C
      implicit none

      integer nb
      integer i 
      integer j
      integer minmn
      INTEGER NAU, NV, IERR, M, N, IRHS
      REAL U(NAU,1), V(NV,1), W(1), B(NB,1)
      LOGICAL MATU, MATV

      MINMN = MIN0(M,N)

      write ( *,99999) MATU, MATV
99999 FORMAT (7H MATU =, L3, 7H MATV =, L3/)
      write ( *,99998) IERR
99998 FORMAT (7H IERR =, I5/)

      write ( *,99997)
99997 FORMAT (17H SINGULAR VALUES./)
      DO I=1,MINMN
        write ( *,99995) W(I)
      end do

      IF (.NOT.MATU) GO TO 30
      write ( *,99996)
99996 FORMAT (/4H U ./)
99995 FORMAT (E15.8)
      DO I=1,M
        write ( *,99993) (U(I,J),J=1,MINMN)
      end do

   30 IF (.NOT.MATV) GO TO 50
      write ( *,99994)
99994 FORMAT (/4H V ./)
99993 FORMAT (5E15.7)
      DO I=1,N
        write ( *,99993) (V(I,J),J=1,MINMN)
      end do

   50 IF (M.LT.N .OR. IRHS.LT.1) GO TO 70
      write ( *,99992)
99992 FORMAT (/4H B ./)
      DO I=1,MINMN
        write ( *,99993) (B(I,J),J=1,IRHS)
      end do

   70 CONTINUE
      RETURN
      END
      SUBROUTINE DIFF(WE, WH, WL, UE, UH, UL, VE, VH, VL, BH, BE, NUE,
     * NUH, NUL, NVE, NVH, NVL, NBH, NBE, M, N, IRHS, MATU, MATV)

c*********************************************************************72
c
cc DIFF computes differences between data computed by different functions.
c
      implicit none

      integer nbe
      integer nbh
      integer nue
      integer nuh
      integer nul
      integer nve
      integer nvh
      integer nvl

      REAL BE(NBE,1), BH(NBH,1)
      real dbhe
      real duel
      real duhe
      real duhl
      real dvel
      real dvhe
      real dvhl
      real dwel
      real dwhe
      real dwhl
      integer i
      integer id
      integer ip1
      INTEGER IRHS
      integer j
      integer krhs
      integer M
      LOGICAL MATU, MATV
      integer minmn
      integer n
      integer nm1
      real t
      REAL WE(1), WH(1), WL(1)
      REAL UE(NUE,1), UH(NUH,1), UL(NUL,1)
      REAL VE(NVE,1), VH(NVH,1), VL(NVL,1)
C
C  SORT SINGULAR VALUES AND EXCHANGE COLUMNS OF U AND V FOR EISPACK.
C  SELECTION SORT MINIMIZES SWAPPING OF U AND V.
C
      DO I = 1, NM1
c
C  FIND INDEX OF MAXIMUM SINGULAR VALUE
c
        ID = I
        IP1 = I + 1
        DO J=IP1,N
          IF (WE(J).GT.WE(ID)) ID = J
        end do
c
C  SWAP SINGULAR VALUES AND VECTORS
c
        IF (ID.NE.I) then

          T = WE(I)
          WE(I) = WE(ID)
          WE(ID) = T
          IF (MATV) CALL SSWAP(N, VE(1,I), 1, VE(1,ID), 1)
          IF (MATU) CALL SSWAP(M, UE(1,I), 1, UE(1,ID), 1)

          DO KRHS=1,IRHS
            T = BE(I,KRHS)
            BE(I,KRHS) = BE(ID,KRHS)
            BE(ID,KRHS) = T
          end do

        end if

      end do

      DWHE = 0.0
      DWHL = 0.0
      DWEL = 0.0
      DUHE = 0.0
      DUHL = 0.0
      DUEL = 0.0
      DVHE = 0.0
      DVHL = 0.0
      DVEL = 0.0
      DBHE = 0.0
      MINMN = MIN0(M,N)

      DO I=1,MINMN
        DWHE = DWHE + (WH(I)-WE(I))**2
        DWHL = DWHL + (WH(I)-WL(I))**2
        DWEL = DWEL + (WE(I)-WL(I))**2
      end do

      DWHE = SQRT(DWHE)
      DWHL = SQRT(DWHL)
      DWEL = SQRT(DWEL)
      write ( *,99999) DWHE, DWHL, DWEL
99999 FORMAT (/35H  2-NORM ( HYBSVD W - EISPAK W ) = , 1PE15.7/6H  2-NO,
     & 29HRM ( HYBSVD W - LINPAK W ) = , 1PE15.7/19H  2-NORM ( EISPAK W,
     & 16H - LINPAK W ) = , 1PE15.7/)
C
C  DIFFERENCE OF ABSOLUTE VALUES BECAUSE OF POSSIBLE SIGN DIFFERENCE.
C
      IF (.NOT.MATU) GO TO 80

      DO J=1,MINMN
        DO I=1,M
          DUHE = DUHE + (ABS(UH(I,J))-ABS(UE(I,J)))**2
          DUHL = DUHL + (ABS(UH(I,J))-ABS(UL(I,J)))**2
          DUEL = DUEL + (ABS(UE(I,J))-ABS(UL(I,J)))**2
        end do
      end do

      DUHE = SQRT(DUHE)
      DUHL = SQRT(DUHL)
      DUEL = SQRT(DUEL)
      write ( *,99998) DUHE, DUHL, DUEL
99998 FORMAT (/35H  2-NORM ( HYBSVD U - EISPAK U ) = , 1PE15.7/6H  2-NO,
     & 29HRM ( HYBSVD U - LINPAK U ) = , 1PE15.7/19H  2-NORM ( EISPAK U,
     & 16H - LINPAK U ) = , 1PE15.7/)

   80 IF (.NOT.MATV) GO TO 110
      DO J=1,MINMN
        DO I=1,N
          DVHE = DVHE + (ABS(VH(I,J))-ABS(VE(I,J)))**2
          DVHL = DVHL + (ABS(VH(I,J))-ABS(VL(I,J)))**2
          DVEL = DVEL + (ABS(VE(I,J))-ABS(VL(I,J)))**2
        end do
      end do
      DVHE = SQRT(DVHE)
      DVHL = SQRT(DVHL)
      DVEL = SQRT(DVEL)
      write ( *,99997) DVHE, DVHL, DVEL
99997 FORMAT (/35H  2-NORM ( HYBSVD V - EISPAK V ) = , 1PE15.7/6H  2-NO,
     & 29HRM ( HYBSVD V - LINPAK V ) = , 1PE15.7/19H  2-NORM ( EISPAK V,
     & 16H - LINPAK V ) = , 1PE15.7/)

  110 IF (M.LT.N .OR. IRHS.LT.1) GO TO 140
      DO J=1,IRHS
        DO I=1,MINMN
          DBHE = DBHE + (ABS(BE(I,J))-ABS(BH(I,J)))**2
        end do
      end do
      DBHE = SQRT(DBHE)
      write ( *,99996) DBHE
99996 FORMAT (/35H  2-NORM ( HYBSVD B - MINFIT B ) = , 1PE15.7/)

  140 CONTINUE
      RETURN
      END
