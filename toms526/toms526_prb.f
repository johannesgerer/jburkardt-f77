      PROGRAM main

c*********************************************************************72
c
cc TTIDBS
c
C THIS PROGRAM IS A TEST PROGRAM FOR THE IDBVIP/IDSFFT SUBPRO-
C GRAM PACKAGE.  ALL ELEMENTS OF RESULTING DZI1 AND DZI2 ARRAYS
C ARE EXPECTED TO BE ZERO.
C THE LUN CONSTANT IN THE LAST DATA INITIALIZATION STATEMENT IS
C THE LOGICAL UNIT NUMBER OF THE STANDARD OUTPUT UNIT AND IS,
C THEREFORE, SYSTEM DEPENDENT.
C DECLARATION STATEMENTS
C
      DIMENSION   XD(30),YD(30),ZD(30),
     1            XI(6),YI(5),ZI(6,5),
     2            ZI1(6,5),ZI2(6,5),DZI1(6,5),DZI2(6,5),
     3            IWK(1030),WK(240)
      DATA  NCP/4/
      DATA  NDP/30/
      DATA  XD(1), XD(2), XD(3), XD(4), XD(5), XD(6),
     1      XD(7), XD(8), XD(9), XD(10),XD(11),XD(12),
     2      XD(13),XD(14),XD(15),XD(16),XD(17),XD(18),
     3      XD(19),XD(20),XD(21),XD(22),XD(23),XD(24),
     4      XD(25),XD(26),XD(27),XD(28),XD(29),XD(30)/
     5      11.16, 24.20, 19.85, 10.35, 19.72,  0.00,
     6      20.87, 19.99, 10.28,  4.51,  0.00, 16.70,
     7       6.08, 25.00, 14.90,  0.00,  9.66,  5.22,
     8      11.77, 15.10, 25.00, 25.00, 14.59, 15.20,
     9       5.23,  2.14,  0.51, 25.00, 21.67,  3.31/
      DATA  YD(1), YD(2), YD(3), YD(4), YD(5), YD(6),
     1      YD(7), YD(8), YD(9), YD(10),YD(11),YD(12),
     2      YD(13),YD(14),YD(15),YD(16),YD(17),YD(18),
     3      YD(19),YD(20),YD(21),YD(22),YD(23),YD(24),
     4      YD(25),YD(26),YD(27),YD(28),YD(29),YD(30)/
     5       1.24, 16.23, 10.72,  4.11,  1.39, 20.00,
     6      20.00,  4.62, 15.16, 20.00,  4.48, 19.65,
     7       4.58, 11.87,  3.12,  0.00, 20.00, 14.66,
     8      10.47, 17.19,  3.87,  0.00,  8.71,  0.00,
     9      10.72, 15.03,  8.37, 20.00, 14.36,  0.13/
      DATA  ZD(1), ZD(2), ZD(3), ZD(4), ZD(5), ZD(6),
     1      ZD(7), ZD(8), ZD(9), ZD(10),ZD(11),ZD(12),
     2      ZD(13),ZD(14),ZD(15),ZD(16),ZD(17),ZD(18),
     3      ZD(19),ZD(20),ZD(21),ZD(22),ZD(23),ZD(24),
     4      ZD(25),ZD(26),ZD(27),ZD(28),ZD(29),ZD(30)/
     5      22.15,  2.83,  7.97, 22.33, 16.83, 34.60,
     6       5.74, 14.72, 21.59, 15.61, 61.77,  6.31,
     7      35.74,  4.40, 21.70, 58.20,  4.73, 40.36,
     8      13.62, 12.57,  8.74, 12.00, 14.81, 21.60,
     9      26.50, 53.10, 49.43,  0.60,  5.52, 44.08/
      DATA  NXI/6/, NYI/5/
      DATA  XI(1), XI(2), XI(3), XI(4), XI(5), XI(6)/
     1       0.00,  5.00, 10.00, 15.00, 20.00, 25.00/
      DATA  YI(1), YI(2), YI(3), YI(4), YI(5)/
     1       0.00,  5.00, 10.00, 15.00, 20.00/
      DATA  ZI(1,1),ZI(2,1),ZI(3,1),ZI(4,1),ZI(5,1),ZI(6,1),
     1      ZI(1,2),ZI(2,2),ZI(3,2),ZI(4,2),ZI(5,2),ZI(6,2),
     2      ZI(1,3),ZI(2,3),ZI(3,3),ZI(4,3),ZI(5,3),ZI(6,3),
     3      ZI(1,4),ZI(2,4),ZI(3,4),ZI(4,4),ZI(5,4),ZI(6,4),
     4      ZI(1,5),ZI(2,5),ZI(3,5),ZI(4,5),ZI(5,5),ZI(6,5)/
     5      58.20, 39.55, 26.90, 21.71, 17.68, 12.00,
     6      61.58, 39.39, 22.04, 21.29, 14.36,  8.04,
     7      59.18, 27.39, 16.78, 13.25,  8.59,  5.36,
     8      52.82, 40.27, 22.76, 16.61,  7.40,  2.88,
     9      34.60, 14.05,  4.12,  3.17,  6.31,  0.60/
      DATA  LUN/6/
C CALCULATION
   10 MD=1
      DO 12  IYI=1,NYI
        DO 11  IXI=1,NXI
          IF(IXI.NE.1.OR.IYI.NE.1)  MD=2
          CALL IDBVIP(MD,NCP,NDP,XD,YD,ZD,1,XI(IXI),YI(IYI),
     1                ZI1(IXI,IYI),IWK,WK)
   11   CONTINUE
   12 CONTINUE
   15 CALL IDSFFT(1,NCP,NDP,XD,YD,ZD,NXI,NYI,XI,YI,ZI2,IWK,WK)
      DO 17  IYI=1,NYI
        DO 16  IXI=1,NXI
          DZI1(IXI,IYI)=ABS(ZI1(IXI,IYI)-ZI(IXI,IYI))
          DZI2(IXI,IYI)=ABS(ZI2(IXI,IYI)-ZI(IXI,IYI))
   16   CONTINUE
   17 CONTINUE
C PRINTING OF INPUT DATA
   20 WRITE (LUN,2020)  NDP
      DO 23  IDP=1,NDP
        IF(MOD(IDP,5).EQ.1)    WRITE (LUN,2021)
        WRITE (LUN,2022)  IDP,XD(IDP),YD(IDP),ZD(IDP)
   23 CONTINUE
C PRINTING OF OUTPUT RESULTS
   30 WRITE (LUN,2030)
      WRITE (LUN,2031)  YI
      DO 33  IXI=1,NXI
        WRITE (LUN,2032)  XI(IXI),(ZI1(IXI,IYI),IYI=1,NYI)
   33 CONTINUE
   40 WRITE (LUN,2040)
      WRITE (LUN,2031)  YI
      DO 43  IXI=1,NXI
        WRITE (LUN,2032)  XI(IXI),(DZI1(IXI,IYI),IYI=1,NYI)
   43 CONTINUE
   50 WRITE (LUN,2050)
      WRITE (LUN,2031)  YI
      DO 53  IXI=1,NXI
        WRITE (LUN,2032)  XI(IXI),(ZI2(IXI,IYI),IYI=1,NYI)
   53 CONTINUE
   60 WRITE (LUN,2060)
      WRITE (LUN,2031)  YI
      DO 63  IXI=1,NXI
        WRITE (LUN,2032)  XI(IXI),(DZI2(IXI,IYI),IYI=1,NYI)
   63 CONTINUE
      STOP
C FORMAT STATEMENTS
 2020 FORMAT(1H1,6HTTIDBS/////3X,10HINPUT DATA,8X,5HNDP =,I3///
     1   30H      I      XD     YD     ZD /)
 2021 FORMAT(1X)
 2022 FORMAT(5X,I2,2X,3F7.2)
 2030 FORMAT(1H1,6HTTIDBS/////3X,17HIDBVIP SUBROUTINE///
     1   26X,10HZI1(XI,YI))
 2031 FORMAT(7X,2HXI,4X,3HYI=/12X,5F7.2/)
 2032 FORMAT(1X/1X,F9.2,2X,5F7.2)
 2040 FORMAT(1X/////3X,10HDIFFERENCE///
     1   25X,11HDZI1(XI,YI))
 2050 FORMAT(1H1,6HTTIDBS/////3X,17HIDSFFT SUBROUTINE///
     1   26X,10HZI2(XI,YI))
 2060 FORMAT(1X/////3X,10HDIFFERENCE///
     1   25X,11HDZI2(XI,YI))
      END
