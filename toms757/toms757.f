      DOUBLE PRECISION FUNCTION ABRAM0(XVALUE)

c*********************************************************************72
c
cc ABRAM0 evaluates the Abramowitz function of order 0.
c
C   DESCRIPTION:
C      This function calculates the Abramowitz function of order 0,
C      defined as
C
C       ABRAM0(x) = integral{ 0 to infinity } exp( -t*t - x/t ) dt
C
C      The code uses Chebyshev expansions with the coefficients
C      given to an accuracy of 20 decimal places.
C
C      This subroutine is set up to work on IEEE machines.
C      For other machines, you should retrieve the code
C      from the general MISCFUN archive.
C
C   ERROR RETURNS:
C      If XVALUE < 0.0, the function prints a message and returns the 
C      value 0.0.
C
C
C   MACHINE-DEPENDENT CONSTANTS:
C
C      NTERMF - INTEGER - No. of terms needed for the array AB0F.
C               Recommended value such that 
C                     ABS( AB0F(NTERMF) ) < EPS/100
C
C      NTERMG - INTEGER - No. of terms needed for array AB0G.
C               Recommended value such that
C                     ABS( AB0G(NTERMG) ) < EPS/100
C
C      NTERMH - INTEGER - No. of terms needed for array AB0H.
C               Recommended value such that
C                     ABS( AB0H(NTERMH) ) < EPS/100
C
C      NTERMA - INTEGER - No. of terms needed for array AB0AS.
C               Recommended value such that
C                     ABS( AB0AS(NTERMA) ) < EPS/100 
C
C     XLOW1 - DOUBLE PRECISION - The value below which 
C              ABRAM0 = root(pi)/2 + X ( ln X - GVAL0 )
C             Recommended value is SQRT(2*EPSNEG)
C
C     LNXMIN - DOUBLE PRECISION - The value of ln XMIN. Used to prevent
C              exponential underflow for large X.
C
C     For values of EPS, EPSNEG, XMIN refer to the file MACHCON.TXT.
C
C     The machine-arithmetic constants are given in DATA
C     statements.
C
C   INTRINSIC FUNCTIONS USED:
C
C     LOG, EXP, SQRT
C
C
C   OTHER MISCFUN SUBROUTINES USED:
C
C          CHEVAL , ERRPRN
C
C
C   AUTHOR:
C
C      DR. ALLAN J. MACLEOD,
C      DEPT. OF MATHEMATICS AND STATISTICS,
C      UNIVERSITY OF PAISLEY ,
C      HIGH ST.,
C      PAISLEY,
C      SCOTLAND.
C      PA1 2BE
C
C      ( e-mail: macl_ms0@paisley.ac.uk ) 
C
C
C   LATEST REVISION:   22 NOVEMBER, 1995
C
C
      INTEGER NTERMA,NTERMF,NTERMG,NTERMH
      DOUBLE PRECISION AB0F(0:8),AB0G(0:8),AB0H(0:8),AB0AS(0:27),
     &     ASLN,ASVAL,CHEVAL,FVAL,GVAL,GVAL0,HALF,HVAL,
     &     LNXMIN,ONEHUN,ONERPI,RTPIB2,RT3BPI,SIX,T,
     &     THREE,TWO,V,X,XLOW1,XVALUE,ZERO
      CHARACTER FNNAME*6,ERRMSG*33
      DATA FNNAME/'ABRAM0'/
      DATA ERRMSG/'FUNCTION CALLED WITH ARGUMENT < 0'/
      DATA AB0F/-0.68121 92709 35494 69816  D    0,
     1          -0.78867 91981 61492 52495  D    0,
     2           0.51215 81776 81881 9543   D   -1,
     3          -0.71092 35289 45412 96     D   -3,
     4           0.36868 18085 04287        D   -5,
     5          -0.91783 23372 37           D   -8,
     6           0.12702 02563              D  -10,
     7          -0.10768 88                 D  -13,
     8           0.599                      D  -17/
      DATA AB0G/-0.60506 03943 08682 73190  D    0,
     1          -0.41950 39816 32017 79803  D    0,
     2           0.17032 65125 19037 0333   D   -1,
     3          -0.16938 91784 24913 97     D   -3,
     4           0.67638 08951 9710         D   -6,
     5          -0.13572 36362 55           D   -8,
     6           0.15629 7065               D  -11,
     7          -0.11288 7                  D  -14,
     8           0.55                       D  -18/
      DATA AB0H/1.38202 65523 05749 89705  D    0,
     1         -0.30097 92907 39749 04355  D    0,
     2          0.79428 88093 64887 241    D   -2,
     3         -0.64319 10276 84756 3      D   -4,
     4          0.22549 83068 4374         D   -6,
     5         -0.41220 96619 5            D   -9,
     6          0.44185 282                D  -12,
     7         -0.30123                    D  -15,
     8          0.14                       D  -18/
      DATA AB0AS(0)/  1.97755 49972 36930 67407  D    0/
      DATA AB0AS(1)/ -0.10460 24792 00481 9485   D   -1/
      DATA AB0AS(2)/  0.69680 79025 36253 66     D   -3/
      DATA AB0AS(3)/ -0.58982 98299 99659 9      D   -4/
      DATA AB0AS(4)/  0.57716 44553 05320        D   -5/
      DATA AB0AS(5)/ -0.61523 01336 5756         D   -6/
      DATA AB0AS(6)/  0.67853 96884 767          D   -7/
      DATA AB0AS(7)/ -0.72306 25379 07           D   -8/
      DATA AB0AS(8)/  0.63306 62736 5            D   -9/
      DATA AB0AS(9)/ -0.98945 3793               D  -11/
      DATA AB0AS(10)/-0.16819 80530              D  -10/
      DATA AB0AS(11)/ 0.67379 9551               D  -11/
      DATA AB0AS(12)/-0.20099 7939               D  -11/
      DATA AB0AS(13)/ 0.54055 903                D  -12/
      DATA AB0AS(14)/-0.13816 679                D  -12/
      DATA AB0AS(15)/ 0.34222 05                 D  -13/
      DATA AB0AS(16)/-0.82668 6                  D  -14/
      DATA AB0AS(17)/ 0.19456 6                  D  -14/
      DATA AB0AS(18)/-0.44268                    D  -15/
      DATA AB0AS(19)/ 0.9562                     D  -16/
      DATA AB0AS(20)/-0.1883                     D  -16/
      DATA AB0AS(21)/ 0.301                      D  -17/
      DATA AB0AS(22)/-0.19                       D  -18/
      DATA AB0AS(23)/-0.14                       D  -18/
      DATA AB0AS(24)/ 0.11                       D  -18/
      DATA AB0AS(25)/-0.4                        D  -19/
      DATA AB0AS(26)/ 0.2                        D  -19/
      DATA AB0AS(27)/-0.1                        D  -19/
      DATA ZERO,HALF,TWO/ 0.0 D 0 , 0.5 D 0, 2.0 D 0/
      DATA THREE,SIX,ONEHUN/ 3.0 D 0, 6.0 D 0 , 100.0 D 0/
      DATA RT3BPI/0.97720 50238 05839 84317 D 0/
      DATA RTPIB2/0.88622 69254 52758 01365 D 0/
      DATA GVAL0/0.13417 65026 47700 70909 D 0/
      DATA ONERPI/0.56418 95835 47756 28695 D 0/
C
C   Machine-dependent constants (suitable for IEEE machines)
C
      DATA NTERMF,NTERMG,NTERMH,NTERMA/8,8,8,22/
      DATA XLOW1,LNXMIN/1.490116D-8,-708.3964D0/
C
C   Start computation
C
      X = XVALUE
C
C   Error test
C
      IF ( X .LT. ZERO ) THEN
         CALL ERRPRN(FNNAME,ERRMSG)
         ABRAM0 = ZERO
         RETURN
      ENDIF
C
C   Code for 0 <= XVALUE <= 2
C
      IF ( X .LE. TWO ) THEN
         IF ( X .EQ. ZERO ) THEN
            ABRAM0 = RTPIB2
            RETURN
         ENDIF
         IF ( X .LT. XLOW1 ) THEN
            ABRAM0 = RTPIB2 + X * ( LOG( X ) - GVAL0 ) 
            RETURN
         ELSE
            T =  ( X * X / TWO - HALF ) - HALF
            FVAL = CHEVAL( NTERMF,AB0F,T ) 
            GVAL = CHEVAL( NTERMG,AB0G,T ) 
            HVAL = CHEVAL( NTERMH,AB0H,T ) 
            ABRAM0 = FVAL/ONERPI + X * ( LOG( X ) * HVAL- GVAL ) 
            RETURN
         ENDIF
      ELSE
C
C   Code for XVALUE > 2
C
         V = THREE *  ( (X / TWO) ** ( TWO / THREE ) ) 
         T =  ( SIX/V - HALF ) - HALF         
         ASVAL = CHEVAL( NTERMA,AB0AS,T ) 
         ASLN = LOG( ASVAL / RT3BPI ) - V
         IF ( ASLN .LT. LNXMIN ) THEN
            ABRAM0 = ZERO
         ELSE
            ABRAM0 = EXP( ASLN ) 
         ENDIF
         RETURN
      ENDIF
      END

      DOUBLE PRECISION FUNCTION ABRAM1(XVALUE)

c*********************************************************************72
c
cc ABRAM1 evaluates the Abramowitz function of order 1.
c
C   DESCRIPTION:
C      This function calculates the Abramowitz function of order 1,
C      defined as
C
C       ABRAM1(x) = integral{ 0 to infinity } t * exp( -t*t - x/t ) dt
C
C      The code uses Chebyshev expansions with the coefficients
C      given to an accuracy of 20 decimal places.
C
C      This subroutine is set up to work on IEEE machines.
C      For other machines, you should retrieve the code
C      from the general MISCFUN archive.
C
C
C   ERROR RETURNS:
C      If XVALUE < 0.0, the function prints a message and returns the 
C      value 0.0.
C
C
C   MACHINE-DEPENDENT CONSTANTS:
C
C      NTERMF - INTEGER - No. of terms needed for the array AB1F.
C               Recommended value such that 
C                     ABS( AB1F(NTERMF) ) < EPS/100
C
C      NTERMG - INTEGER - No. of terms needed for array AB1G.
C               Recommended value such that
C                     ABS( AB1G(NTERMG) ) < EPS/100
C
C      NTERMH - INTEGER - No. of terms needed for array AB1H.
C               Recommended value such that
C                     ABS( AB1H(NTERMH) ) < EPS/100
C
C      NTERMA - INTEGER - No. of terms needed for array AB1AS.
C               Recommended value such that
C                     ABS( AB1AS(NTERMA) ) < EPS/100 
C
C      XLOW - DOUBLE PRECISION - The value below which
C                ABRAM1(x) = 0.5 to machine precision.
C             The recommended value is EPSNEG/2
C
C      XLOW1 - DOUBLE PRECISION - The value below which 
C                ABRAM1(x) = (1 - x ( sqrt(pi) + xln(x) ) / 2 
C              Recommended value is SQRT(2*EPSNEG)
C
C      LNXMIN - DOUBLE PRECISION - The value of ln XMIN. Used to prevent
C              exponential underflow for large X.
C
C      For values of EPS, EPSNEG, XMIN refer to the file MACHCON.TXT
C
C      The machine-arithmetic constants are given in DATA
C      statements.
C
C
C   INTRINSIC FUNCTIONS USED:
C
C     LOG, EXP, SQRT
C
C
C   OTHER MISCFUN SUBROUTINES USED:
C
C          CHEVAL , ERRPRN
C
C
C   AUTHOR:
C
C      DR. ALLAN J. MACLEOD,
C      DEPT. OF MATHEMATICS AND STATISTICS,
C      UNIVERSITY OF PAISLEY,
C      HIGH ST.,
C      PAISLEY,
C      SCOTLAND.
C      PA1 2BE
C
C      ( e-mail: macl_ms0@paisley.ac.uk ) 
C
C
C   LATEST REVISION:   22 NOVEMBER, 1995
C
C
      INTEGER NTERMA,NTERMF,NTERMG,NTERMH
      DOUBLE PRECISION AB1F(0:9),AB1G(0:8),AB1H(0:8),AB1AS(0:27),
     &     ASLN,ASVAL,CHEVAL,FVAL,GVAL,HALF,HVAL,
     &     LNXMIN,ONE,ONEHUN,ONERPI,RT3BPI,SIX,T,THREE,TWO,
     &     V,X,XLOW,XLOW1,XVALUE,ZERO
      CHARACTER FNNAME*6,ERRMSG*33
      DATA FNNAME/'ABRAM1'/
      DATA ERRMSG/'FUNCTION CALLED WITH ARGUMENT < 0'/
      DATA AB1F/1.47285 19257 79788 07369  D    0,
     1          0.10903 49757 01689 56257  D    0,
     2         -0.12430 67536 00565 69753  D    0,
     3          0.30619 79468 53493 315    D   -2,
     4         -0.22184 10323 07651 1      D   -4,
     5          0.69899 78834 451          D   -7,
     6         -0.11597 07644 4            D   -9,
     7          0.11389 776                D  -12,
     8         -0.7173                     D  -16,
     9          0.3                        D  -19/
      DATA AB1G/0.39791 27794 90545 03528  D    0,
     1         -0.29045 28522 64547 20849  D    0,
     2          0.10487 84695 46536 3504   D   -1,
     3         -0.10249 86952 26913 36     D   -3,
     4          0.41150 27939 9110         D   -6,
     5         -0.83652 63894 0            D   -9,
     6          0.97862 595                D  -12,
     7         -0.71868                    D  -15,
     8          0.35                       D  -18/
      DATA AB1H/0.84150 29215 22749 47030  D    0,
     1         -0.77900 50698 77414 3395   D   -1,
     2          0.13399 24558 78390 993    D   -2,
     3         -0.80850 39071 52788        D   -5,
     4          0.22618 58281 728          D   -7,
     5         -0.34413 95838              D  -10,
     6          0.31598 58                 D  -13,
     7         -0.1884                     D  -16,
     8          0.1                        D  -19/
      DATA AB1AS(0)/  2.13013 64342 90655 49448  D    0/
      DATA AB1AS(1)/  0.63715 26795 21853 9933   D   -1/
      DATA AB1AS(2)/ -0.12933 49174 77510 647    D   -2/
      DATA AB1AS(3)/  0.56783 28753 22826 5      D   -4/
      DATA AB1AS(4)/ -0.27943 49391 77646        D   -5/
      DATA AB1AS(5)/  0.56002 14736 787          D   -7/
      DATA AB1AS(6)/  0.23920 09242 798          D   -7/
      DATA AB1AS(7)/ -0.75098 48650 09           D   -8/
      DATA AB1AS(8)/  0.17301 53307 76           D   -8/
      DATA AB1AS(9)/ -0.36648 87795 5            D   -9/
      DATA AB1AS(10)/ 0.75207 58307              D  -10/
      DATA AB1AS(11)/-0.15179 90208              D  -10/
      DATA AB1AS(12)/ 0.30171 3710               D  -11/
      DATA AB1AS(13)/-0.58596 718                D  -12/
      DATA AB1AS(14)/ 0.10914 455                D  -12/
      DATA AB1AS(15)/-0.18705 36                 D  -13/
      DATA AB1AS(16)/ 0.26254 2                  D  -14/
      DATA AB1AS(17)/-0.14627                    D  -15/
      DATA AB1AS(18)/-0.9500                     D  -16/
      DATA AB1AS(19)/ 0.5873                     D  -16/
      DATA AB1AS(20)/-0.2420                     D  -16/
      DATA AB1AS(21)/ 0.868                      D  -17/
      DATA AB1AS(22)/-0.290                      D  -17/
      DATA AB1AS(23)/ 0.93                       D  -18/
      DATA AB1AS(24)/-0.29                       D  -18/
      DATA AB1AS(25)/ 0.9                        D  -19/
      DATA AB1AS(26)/-0.3                        D  -19/
      DATA AB1AS(27)/ 0.1                        D  -19/
      DATA ZERO,HALF,ONE/ 0.0 D 0, 0.5 D 0, 1.0 D 0/
      DATA TWO,THREE,SIX/ 2.0 D 0, 3.0 D 0, 6.0 D 0/
      DATA ONEHUN/100.0 D 0/
      DATA RT3BPI/ 0.97720 50238 05839 84317 D 0/
      DATA ONERPI/ 0.56418 95835 47756 28695 D 0/
C
C   Machine-dependent constants (suitable for IEEE machines)
C
      DATA NTERMF,NTERMG,NTERMH,NTERMA/9,8,8,23/
      DATA XLOW,XLOW1,LNXMIN/1.11023D-16,1.490116D-8,-708.3964D0/
C
C   Start calculation
C
      X = XVALUE
C
C   Error test
C
      IF ( X .LT. ZERO ) THEN
         CALL ERRPRN(FNNAME,ERRMSG)
         ABRAM1 = ZERO
         RETURN
      ENDIF   
C
C   Code for 0 <= XVALUE <= 2
C
      IF ( X .LE. TWO ) THEN
         IF ( X .EQ. ZERO ) THEN
            ABRAM1 = HALF
            RETURN
         ENDIF
         IF ( X .LT. XLOW1 ) THEN
            IF ( X .LT. XLOW ) THEN
               ABRAM1 = HALF
            ELSE
               ABRAM1 = ( ONE - X / ONERPI - X * X * LOG( X ) ) * HALF
            ENDIF
            RETURN
         ELSE
            T =  ( X * X / TWO - HALF ) - HALF
            FVAL = CHEVAL( NTERMF,AB1F,T ) 
            GVAL = CHEVAL( NTERMG,AB1G,T ) 
            HVAL = CHEVAL( NTERMH,AB1H,T ) 
            ABRAM1 = FVAL - X * ( GVAL / ONERPI + X * LOG( X ) * HVAL )
            RETURN
         ENDIF
      ELSE
C
C   Code for XVALUE > 2
C
         V = THREE *  ( (X / TWO) ** ( TWO / THREE ) ) 
         T =  ( SIX / V - HALF ) - HALF         
         ASVAL = CHEVAL( NTERMA,AB1AS,T ) 
         ASLN = LOG( ASVAL * SQRT ( V / THREE ) / RT3BPI ) - V
         IF ( ASLN .LT. LNXMIN ) THEN
            ABRAM1 = ZERO
         ELSE
            ABRAM1 = EXP( ASLN ) 
         ENDIF
         RETURN
      ENDIF
      END

      DOUBLE PRECISION FUNCTION ABRAM2(XVALUE)

c*********************************************************************72
c
cc ABRAM2 evaluates the Abramowitz function of order 2.
c
C   DESCRIPTION:
C      This function calculates the Abramowitz function of order 2,
C      defined as
C
C       ABRAM2(x) = integral{ 0 to infinity } (t**2) * exp( -t*t - x/t ) dt
C
C      The code uses Chebyshev expansions with the coefficients
C      given to an accuracy of 20 decimal places. 
C
C      This subroutine is set up to work on IEEE machines.
C      For other machines, you should retrieve the code
C      from the general MISCFUN archive.
C
C   ERROR RETURNS:
C      If XVALUE < 0.0, the function prints a message and returns the 
C      value 0.0.
C
C
C   MACHINE-DEPENDENT CONSTANTS:
C
C      NTERMF - INTEGER - No. of terms needed for the array AB2F.
C               Recommended value such that 
C                     ABS( AB2F(NTERMF) ) < EPS/100
C
C      NTERMG - INTEGER - No. of terms needed for array AB2G.
C               Recommended value such that
C                     ABS( AB2G(NTERMG) ) < EPS/100
C
C      NTERMH - INTEGER - No. of terms needed for array AB2H.
C               Recommended value such that
C                     ABS( AB2H(NTERMH) ) < EPS/100
C
C      NTERMA - INTEGER - No. of terms needed for array AB2AS.
C               Recommended value such that
C                     ABS( AB2AS(NTERMA) ) < EPS/100 
C
C      XLOW - DOUBLE PRECISION - The value below which 
C               ABRAM2 = root(pi)/4 to machine precision.
C             The recommended value is EPSNEG
C
C      XLOW1 - DOUBLE PRECISION - The value below which 
C                ABRAM2 = root(pi)/4 - x/2 + x**3ln(x)/6
C              Recommended value is SQRT(2*EPSNEG)
C
C      LNXMIN - DOUBLE PRECISION - The value of ln XMIN. Used to prevent
C               exponential underflow for large X.
C
C     For values of EPS, EPSNEG, XMIN refer to the file MACHCON.TXT
C
C     The machine-arithmetic constants are given in DATA
C     statements.
C
C
C   INTRINSIC FUNCTIONS USED:
C
C     LOG, EXP
C
C
C   OTHER MISCFUN SUBROUTINES USED:
C
C          CHEVAL , ERRPRN
C
C
C   AUTHOR:
C
C      DR. ALLAN J. MACLEOD,
C      DEPT. OF MATHEMATICS AND STATISTICS,
C      UNIVERSITY OF PAISLEY,
C      HIGH ST.,
C      PAISLEY,
C      SCOTLAND.
C      PA1 2BE
C
C      ( e-mail: macl_ms0@paisley.ac.uk ) 
C
C
C   LATEST REVISION:   22 NOVEMBER , 1995
C
C
      INTEGER NTERMA,NTERMF,NTERMG,NTERMH
      DOUBLE PRECISION AB2F(0:9),AB2G(0:8),AB2H(0:7),AB2AS(0:26),
     &     ASLN,ASVAL,CHEVAL,FVAL,GVAL,HALF,HVAL,LNXMIN,
     &     ONEHUN,ONERPI,RTPIB4,RT3BPI,SIX,T,THREE,TWO,
     &     V,X,XLOW,XLOW1,XVALUE,ZERO
      CHARACTER FNNAME*6,ERRMSG*33
      DATA FNNAME/'ABRAM2'/
      DATA ERRMSG/'FUNCTION CALLED WITH ARGUMENT < 0'/
      DATA AB2F/1.03612 16280 42437 13846  D    0,
     1          0.19371 24662 67945 70012  D    0,
     2         -0.72587 58839 23300 7378   D   -1,
     3          0.17479 05908 64327 399    D   -2,
     4         -0.12812 23233 75654 9      D   -4,
     5          0.41150 18153 651          D   -7,
     6         -0.69710 47256              D  -10,
     7          0.69901 83                 D  -13,
     8         -0.4492                     D  -16,
     9          0.2                        D  -19/
      DATA AB2G/1.46290 15719 86307 41150  D    0,
     1          0.20189 46688 31540 14317  D    0,
     2         -0.29082 92087 99712 9022   D   -1,
     3          0.47061 04903 52700 50     D   -3,
     4         -0.25792 20803 59333        D   -5,
     5          0.65613 37129 46           D   -8,
     6         -0.91411 0203               D  -11,
     7          0.77427 6                  D  -14,
     8         -0.429                      D  -17/
      DATA AB2H/0.30117 22501 09104 88881  D    0,
     1         -0.15886 67818 31762 3783   D   -1,
     2          0.19295 93693 55845 26     D   -3,
     3         -0.90199 58784 9300         D   -6,
     4          0.20610 50418 37           D   -8,
     5         -0.26511 1806               D  -11,
     6          0.21086 4                  D  -14,
     7         -0.111                      D  -17/
      DATA AB2AS(0)/  2.46492 32530 43348 56893  D    0/
      DATA AB2AS(1)/  0.23142 79742 22489 05432  D    0/
      DATA AB2AS(2)/ -0.94068 17301 00857 73     D   -3/
      DATA AB2AS(3)/  0.82902 70038 08973 3      D   -4/
      DATA AB2AS(4)/ -0.88389 47042 45866        D   -5/
      DATA AB2AS(5)/  0.10663 85435 67985        D   -5/
      DATA AB2AS(6)/ -0.13991 12853 8529         D   -6/
      DATA AB2AS(7)/  0.19397 93208 445          D   -7/
      DATA AB2AS(8)/ -0.27704 99383 75           D   -8/
      DATA AB2AS(9)/  0.39590 68718 6            D   -9/
      DATA AB2AS(10)/-0.54083 54342              D  -10/
      DATA AB2AS(11)/ 0.63554 6076               D  -11/
      DATA AB2AS(12)/-0.38461 613                D  -12/
      DATA AB2AS(13)/-0.11696 067                D  -12/
      DATA AB2AS(14)/ 0.68966 71                 D  -13/
      DATA AB2AS(15)/-0.25031 13                 D  -13/
      DATA AB2AS(16)/ 0.78558 6                  D  -14/
      DATA AB2AS(17)/-0.23033 4                  D  -14/
      DATA AB2AS(18)/ 0.64914                    D  -15/
      DATA AB2AS(19)/-0.17797                    D  -15/
      DATA AB2AS(20)/ 0.4766                     D  -16/
      DATA AB2AS(21)/-0.1246                     D  -16/
      DATA AB2AS(22)/ 0.316                      D  -17/
      DATA AB2AS(23)/-0.77                       D  -18/
      DATA AB2AS(24)/ 0.18                       D  -18/
      DATA AB2AS(25)/-0.4                        D  -19/
      DATA AB2AS(26)/ 0.1                        D  -19/
      DATA ZERO,HALF,TWO/ 0.0 D 0 , 0.5 D 0, 2.0 D 0/
      DATA THREE,SIX,ONEHUN/ 3.0 D 0, 6.0 D 0 , 100.0 D 0/
      DATA RT3BPI/ 0.97720 50238 05839 84317 D 0/
      DATA RTPIB4/ 0.44311 34627 26379 00682 D 0/
      DATA ONERPI/ 0.56418 95835 47756 28695 D 0/
C
C   Machine-dependent constants (suitable for IEEE machines)
C
      DATA NTERMF,NTERMG,NTERMH,NTERMA/9,8,7,23/
      DATA XLOW,XLOW1,LNXMIN/2.22045D-16,1.490116D-8,-708.3964D0/
C
C   Start calculation
C
      X = XVALUE
C
C   Error test
C
      IF ( X .LT. ZERO ) THEN
         CALL ERRPRN(FNNAME,ERRMSG)
         ABRAM2 = ZERO
         RETURN
      ENDIF
C
C   Code for 0 <= XVALUE <= 2
C
      IF ( X .LE. TWO ) THEN
         IF ( X .EQ. ZERO ) THEN
            ABRAM2 = RTPIB4
            RETURN
         ENDIF
         IF ( X .LT. XLOW1 ) THEN
            IF ( X .LT. XLOW ) THEN
               ABRAM2 = RTPIB4
            ELSE
               ABRAM2 = RTPIB4 - HALF * X + X * X * X * LOG( X ) / SIX
            ENDIF 
            RETURN
         ELSE
            T =  ( X * X / TWO - HALF ) - HALF
            FVAL = CHEVAL( NTERMF,AB2F,T ) 
            GVAL = CHEVAL( NTERMG,AB2G,T ) 
            HVAL = CHEVAL( NTERMH,AB2H,T ) 
            ABRAM2 = FVAL/ONERPI + X * ( X * X * LOG(X) * HVAL- GVAL ) 
            RETURN
         ENDIF
      ELSE
C
C   Code for XVALUE > 2
C
         V = THREE *  ( (X / TWO) ** ( TWO / THREE ) ) 
         T =  ( SIX / V - HALF ) - HALF         
         ASVAL = CHEVAL( NTERMA,AB2AS,T ) 
         ASLN = LOG( ASVAL / RT3BPI ) + LOG( V / THREE ) - V
         IF ( ASLN .LT. LNXMIN ) THEN
            ABRAM2 = ZERO
         ELSE
            ABRAM2 = EXP( ASLN ) 
         ENDIF
         RETURN
      ENDIF
      END

      DOUBLE PRECISION FUNCTION AIRINT(XVALUE)

c*********************************************************************72
c
cc AIRINT calculates the integral of the Airy function Ai.
c
C   DESCRIPTION:
C
C      This function calculates the integral of the Airy function Ai,
C      defined as
C
C         AIRINT(x) = {integral 0 to x} Ai(t) dt
C
C      The program uses Chebyshev expansions, the coefficients of which
C      are given to 20 decimal places.
C
C      This subroutine is set up to work on IEEE machines.
C      For other machines, you should retrieve the code
C      from the general MISCFUN archive.
C
C
C   ERROR RETURNS:
C
C      If the argument is too large and negative, it is impossible
C      to accurately compute the necessary SIN and COS functions.
C      An error message is printed, and the program returns the
C      value -2/3 (the value at -infinity).
C
C
C   MACHINE-DEPENDENT CONSTANTS:
C
C      NTERM1 - INTEGER - The no. of terms to be used from the array
C                          AAINT1. The recommended value is such that
C                             ABS(AAINT1(NTERM1)) < EPS/100,
C                          subject to 1 <= NTERM1 <= 25.
C
C      NTERM2 - INTEGER - The no. of terms to be used from the array
C                          AAINT2. The recommended value is such that
C                             ABS(AAINT2(NTERM2)) < EPS/100,
C                          subject to 1 <= NTERM2 <= 21.
C
C      NTERM3 - INTEGER - The no. of terms to be used from the array
C                          AAINT3. The recommended value is such that
C                             ABS(AAINT3(NTERM3)) < EPS/100,
C                          subject to 1 <= NTERM3 <= 40.
C 
C      NTERM4 - INTEGER - The no. of terms to be used from the array
C                          AAINT4. The recommended value is such that
C                             ABS(AAINT4(NTERM4)) < EPS/100,
C                          subject to 1 <= NTERM4 <= 17.
C
C      NTERM5 - INTEGER - The no. of terms to be used from the array
C                          AAINT5. The recommended value is such that
C                             ABS(AAINT5(NTERM5)) < EPS/100,
C                          subject to 1 <= NTERM5 <= 17.
C
C      XLOW1 - DOUBLE PRECISION - The value such that, if |x| < XLOW1,
C                          AIRINT(x) = x * Ai(0)
C                     to machine precision. The recommended value is
C                          2 * EPSNEG.
C
C      XHIGH1 - DOUBLE PRECISION - The value such that, if x > XHIGH1,
C                          AIRINT(x) = 1/3,
C                      to machine precision. The recommended value is
C                          (-1.5*LOG(EPSNEG)) ** (2/3).
C
C      XNEG1 - DOUBLE PRECISION - The value such that, if x < XNEG1,
C                     the trigonometric functions in the asymptotic
C                     expansion cannot be calculated accurately.
C                     The recommended value is
C                          -(1/((EPS)**2/3))
C
C      For values of EPS and EPSNEG, refer to the file MACHCON.TXT.
C
C      The machine-arithmetic constants are given in DATA
C      statements.
C
C 
C   INTRINSIC FUNCTIONS USED:
C                            COS, EXP, SIN, SQRT
C
C
C   OTHER MISCFUN SUBROUTINES USED:
C
C          CHEVAL , ERRPRN
C
C
C   AUTHOR: Dr. Allan J. MacLeod,
C           Dept. of Mathematics and Statistics,
C           Univ. of Paisley,
C           High St.,
C           Paisley,
C           SCOTLAND.
C           PA1 2BE
C 
C           (e-mail:macl_ms0@paisley.ac.uk)
C
C
C   LATEST REVISION: 22 NOVEMBER, 1995.
C
      INTEGER NTERM1,NTERM2,NTERM3,NTERM4,NTERM5
      DOUBLE PRECISION AAINT1(0:25),AAINT2(0:21),AAINT3(0:40),
     1     AAINT4(0:17),AAINT5(0:17),
     2     AIRZER,ARG,CHEVAL,EIGHT,FORTY1,FOUR,FR996,GVAL,
     3     HVAL,NINE,NINHUN,ONE,ONEHUN,PIBY4,PITIM6,RT2B3P,T,TEMP,
     4     THREE,TWO,X,XHIGH1,XLOW1,XNEG1,XVALUE,Z,ZERO
      CHARACTER FNNAME*6,ERRMSG*46
      DATA FNNAME/'AIRINT'/
      DATA ERRMSG/'FUNCTION TOO NEGATIVE FOR ACCURATE COMPUTATION'/
      DATA AAINT1(0)/  0.37713 51769 46836 95526  D    0/
      DATA AAINT1(1)/ -0.13318 86843 24079 47431  D    0/
      DATA AAINT1(2)/  0.31524 97374 78288 4809   D   -1/
      DATA AAINT1(3)/ -0.31854 30764 36574 077    D   -2/
      DATA AAINT1(4)/ -0.87398 76469 86219 15     D   -3/
      DATA AAINT1(5)/  0.46699 49765 53969 71     D   -3/
      DATA AAINT1(6)/ -0.95449 36738 98369 2      D   -4/
      DATA AAINT1(7)/  0.54270 56871 56716        D   -5/
      DATA AAINT1(8)/  0.23949 64062 52188        D   -5/
      DATA AAINT1(9)/ -0.75690 27020 5649         D   -6/
      DATA AAINT1(10)/ 0.90501 38584 518          D   -7/
      DATA AAINT1(11)/ 0.32052 94560 43           D   -8/
      DATA AAINT1(12)/-0.30382 55364 44           D   -8/
      DATA AAINT1(13)/ 0.48900 11859 6            D   -9/
      DATA AAINT1(14)/-0.18398 20572              D  -10/
      DATA AAINT1(15)/-0.71124 7519               D  -11/
      DATA AAINT1(16)/ 0.15177 4419               D  -11/
      DATA AAINT1(17)/-0.10801 922                D  -12/
      DATA AAINT1(18)/-0.96354 2                  D  -14/
      DATA AAINT1(19)/ 0.31342 5                  D  -14/
      DATA AAINT1(20)/-0.29446                    D  -15/
      DATA AAINT1(21)/-0.477                      D  -17/
      DATA AAINT1(22)/ 0.461                      D  -17/
      DATA AAINT1(23)/-0.53                       D  -18/
      DATA AAINT1(24)/ 0.1                        D  -19/
      DATA AAINT1(25)/ 0.1                        D  -19/
      DATA AAINT2(0)/  1.92002 52408 19840 09769  D    0/
      DATA AAINT2(1)/ -0.42200 49417 25628 7021   D   -1/
      DATA AAINT2(2)/ -0.23945 77229 65939 223    D   -2/
      DATA AAINT2(3)/ -0.19564 07048 33529 71     D   -3/
      DATA AAINT2(4)/ -0.15472 52891 05611 2      D   -4/
      DATA AAINT2(5)/ -0.14049 01861 37889        D   -5/
      DATA AAINT2(6)/ -0.12128 01427 1367         D   -6/
      DATA AAINT2(7)/ -0.11791 86050 192          D   -7/
      DATA AAINT2(8)/ -0.10431 55787 88           D   -8/
      DATA AAINT2(9)/ -0.10908 20929 3            D   -9/
      DATA AAINT2(10)/-0.92963 3045               D  -11/
      DATA AAINT2(11)/-0.11094 6520               D  -11/
      DATA AAINT2(12)/-0.78164 83                 D  -13/
      DATA AAINT2(13)/-0.13196 61                 D  -13/
      DATA AAINT2(14)/-0.36823                    D  -15/
      DATA AAINT2(15)/-0.21505                    D  -15/
      DATA AAINT2(16)/ 0.1238                     D  -16/
      DATA AAINT2(17)/-0.557                      D  -17/
      DATA AAINT2(18)/ 0.84                       D  -18/
      DATA AAINT2(19)/-0.21                       D  -18/
      DATA AAINT2(20)/ 0.4                        D  -19/
      DATA AAINT2(21)/-0.1                        D  -19/
      DATA AAINT3(0)/  0.47985 89326 47910 52053  D    0/
      DATA AAINT3(1)/ -0.19272 37512 61696 08863  D    0/
      DATA AAINT3(2)/  0.20511 54129 52542 8189   D   -1/
      DATA AAINT3(3)/  0.63320 00070 73248 8786   D   -1/
      DATA AAINT3(4)/ -0.50933 22261 84575 4082   D   -1/
      DATA AAINT3(5)/  0.12844 24078 66166 3016   D   -1/
      DATA AAINT3(6)/  0.27601 37088 98947 9413   D   -1/
      DATA AAINT3(7)/ -0.15470 66673 86664 9507   D   -1/
      DATA AAINT3(8)/ -0.14968 64655 38931 6026   D   -1/
      DATA AAINT3(9)/  0.33661 76141 73574 541    D   -2/
      DATA AAINT3(10)/ 0.53085 11635 18892 985    D   -2/
      DATA AAINT3(11)/ 0.41371 22645 85550 81     D   -3/
      DATA AAINT3(12)/-0.10249 05799 26726 266    D   -2/
      DATA AAINT3(13)/-0.32508 22167 20258 53     D   -3/
      DATA AAINT3(14)/ 0.86086 60957 16921 3      D   -4/
      DATA AAINT3(15)/ 0.66713 67298 12077 5      D   -4/
      DATA AAINT3(16)/ 0.44920 59993 18095        D   -5/
      DATA AAINT3(17)/-0.67042 72309 58249        D   -5/
      DATA AAINT3(18)/-0.19663 65700 85009        D   -5/
      DATA AAINT3(19)/ 0.22229 67740 7226         D   -6/
      DATA AAINT3(20)/ 0.22332 22294 9137         D   -6/
      DATA AAINT3(21)/ 0.28033 13766 457          D   -7/
      DATA AAINT3(22)/-0.11556 51663 619          D   -7/
      DATA AAINT3(23)/-0.43306 98217 36           D   -8/
      DATA AAINT3(24)/-0.62277 77938              D  -10/
      DATA AAINT3(25)/ 0.26432 66490 3            D   -9/
      DATA AAINT3(26)/ 0.53338 81114              D  -10/
      DATA AAINT3(27)/-0.52295 7269               D  -11/
      DATA AAINT3(28)/-0.38222 9283               D  -11/
      DATA AAINT3(29)/-0.40958 233                D  -12/
      DATA AAINT3(30)/ 0.11515 622                D  -12/
      DATA AAINT3(31)/ 0.38757 66                 D  -13/
      DATA AAINT3(32)/ 0.14028 3                  D  -14/
      DATA AAINT3(33)/-0.14152 6                  D  -14/
      DATA AAINT3(34)/-0.28746                    D  -15/
      DATA AAINT3(35)/ 0.923                      D  -17/
      DATA AAINT3(36)/ 0.1224                     D  -16/
      DATA AAINT3(37)/ 0.157                      D  -17/
      DATA AAINT3(38)/-0.19                       D  -18/
      DATA AAINT3(39)/-0.8                        D  -19/
      DATA AAINT3(40)/-0.1                        D  -19/
      DATA AAINT4/1.99653 30582 85227 30048  D    0,
     1           -0.18754 11776 05417 759    D   -2,
     2           -0.15377 53628 03057 50     D   -3,
     3           -0.12831 12967 68234 9      D   -4,
     4           -0.10812 84819 64162        D   -5,
     5           -0.91821 31174 057          D   -7,
     6           -0.78416 05909 60           D   -8,
     7           -0.67292 45387 8            D   -9,
     8           -0.57963 25198              D  -10,
     9           -0.50104 0991               D  -11,
     X           -0.43420 222                D  -12,
     1           -0.37743 05                 D  -13,
     2           -0.32847 3                  D  -14,
     3           -0.28700                    D  -15,
     4           -0.2502                     D  -16,
     5           -0.220                      D  -17,
     6           -0.19                       D  -18,
     7           -0.2                        D  -19/
      DATA AAINT5/1.13024 60203 44657 16133  D    0,
     1           -0.46471 80646 39872 334    D   -2,
     2           -0.35137 41338 26932 03     D   -3,
     3           -0.27681 17872 54518 5      D   -4,
     4           -0.22205 74525 58107        D   -5,
     5           -0.18089 14236 5974         D   -6,
     6           -0.14876 13383 373          D   -7,
     7           -0.12351 53881 68           D   -8,
     8           -0.10310 10425 7            D   -9,
     9           -0.86749 3013               D  -11,
     X           -0.73080 054                D  -12,
     1           -0.62235 61                 D  -13,
     2           -0.52512 8                  D  -14,
     3           -0.45677                    D  -15,
     4           -0.3748                     D  -16,
     5           -0.356                      D  -17,
     6           -0.23                       D  -18,
     7           -0.4                        D  -19/
      DATA ZERO,ONE,TWO/ 0.0 D 0 , 1.0 D 0 , 2.0 D 0 /
      DATA THREE,FOUR,EIGHT/ 3.0 D 0 , 4.0 D 0 , 8.0 D 0 /
      DATA NINE,FORTY1,ONEHUN/ 9.0 D 0 , 41.0 D 0 , 100.0 D 0/
      DATA NINHUN,FR996/ 900.0 D 0 , 4996.0 D 0 /
      DATA PIBY4/0.78539 81633 97448 30962 D 0/
      DATA PITIM6/18.84955 59215 38759 43078 D 0/
      DATA RT2B3P/0.46065 88659 61780 63902 D 0/
      DATA AIRZER/0.35502 80538 87817 23926 D 0/
C
C   Machine-dependant constants (suitable for IEEE machines)
C
      DATA NTERM1,NTERM2,NTERM3,NTERM4,NTERM5/22,17,37,15,15/
      DATA XLOW1,XHIGH1,XNEG1/2.22045D-16,14.480884D0,-2.727134D10/
C
C   Start computation
C
      X = XVALUE
C
C   Error test
C
      IF ( X .LT. XNEG1 ) THEN
         CALL ERRPRN(FNNAME,ERRMSG)
         AIRINT = -TWO / THREE
         RETURN
      ENDIF
C
C   Code for x >= 0
C
      IF ( X .GE. ZERO ) THEN
         IF ( X .LE. FOUR ) THEN
            IF ( X .LT. XLOW1 ) THEN
               AIRINT = AIRZER * X
            ELSE
               T = X / TWO - ONE
               AIRINT = CHEVAL(NTERM1,AAINT1,T) * X
            ENDIF
         ELSE
            IF ( X .GT. XHIGH1 ) THEN
               TEMP = ZERO
            ELSE 
               Z = ( X + X ) * SQRT(X) / THREE
               TEMP = THREE * Z
               T = ( FORTY1 - TEMP ) / ( NINE + TEMP )
               TEMP = EXP(-Z) * CHEVAL(NTERM2,AAINT2,T) / SQRT(PITIM6*Z)
            ENDIF
            AIRINT = ONE / THREE - TEMP
         ENDIF
      ELSE
C
C   Code for x < 0
C
         IF ( X .GE. -EIGHT ) THEN
            IF ( X .GT. -XLOW1 ) THEN
               AIRINT = AIRZER * X
            ELSE
               T = -X / FOUR - ONE
               AIRINT = X * CHEVAL(NTERM3,AAINT3,T)
            ENDIF
         ELSE
            Z = - ( X + X ) * SQRT(-X) / THREE
            ARG = Z + PIBY4
            TEMP = NINE * Z * Z
            T = ( FR996 - TEMP ) / ( NINHUN + TEMP)
            GVAL = CHEVAL(NTERM4,AAINT4,T)
            HVAL = CHEVAL(NTERM5,AAINT5,T)
            TEMP = GVAL * COS(ARG) + HVAL * SIN(ARG) / Z
            AIRINT = RT2B3P * TEMP / SQRT(Z) - TWO / THREE
         ENDIF
      ENDIF
      RETURN
      END

      DOUBLE PRECISION FUNCTION AIRYGI(XVALUE)

c*********************************************************************72
c
cc AIRYGI computes the modified Airy function Gi(x).
c
C   DESCRIPTION:
C
C      This subroutine computes the modified Airy function Gi(x),
C      defined as
C
C        AIRYGI(x) = [ Integral{0 to infinity} sin(x*t+t^3/3) dt ] / pi
C
C      The approximation uses Chebyshev expansions with the coefficients 
C      given to 20 decimal places.
C
C      This subroutine is set up to work on IEEE machines.
C      For other machines, you should retrieve the code
C      from the general MISCFUN archive.
C
C
C   ERROR RETURNS:
C
C      If x < -XHIGH1*XHIGH1 (see below for definition of XHIGH1), then
C      the trig. functions needed for the asymptotic expansion of Bi(x)
C      cannot be computed to any accuracy. An error message is printed
C      and the code returns the value 0.0.
C 
C
C   MACHINE-DEPENDENT CONSTANTS:
C
C      NTERM1 - INTEGER - The no. of terms to be used from the array
C                         ARGIP1. The recommended value is such that
C                                ABS(ARGIP1(NTERM1)) < EPS/100
C                         subject to 1 <= NTERM1 <= 30.
C
C      NTERM2 - INTEGER - The no. of terms to be used from the array
C                         ARGIP2. The recommended value is such that
C                                ABS(ARGIP2(NTERM2)) < EPS/100
C                         subject to 1 <= NTERM2 <= 29.
C
C      NTERM3 - INTEGER - The no. of terms to be used from the array
C                         ARGIN1. The recommended value is such that
C                                ABS(ARGIN1(NTERM3)) < EPS/100
C                         subject to 1 <= NTERM3 <= 42.
C
C      NTERM4 - INTEGER - The no. of terms to be used from the array
C                         ARBIN1. The recommended value is such that
C                                ABS(ARBIN1(NTERM4)) < EPS/100
C                         subject to 1 <= NTERM4 <= 10.
C
C      NTERM5 - INTEGER - The no. of terms to be used from the array
C                         ARBIN2. The recommended value is such that
C                                ABS(ARBIN2(NTERM5)) < EPS/100
C                         subject to 1 <= NTERM5 <= 11.
C
C      NTERM6 - INTEGER - The no. of terms to be used from the array
C                         ARGH2. The recommended value is such that
C                                ABS(ARHIN1(NTERM6)) < EPS/100
C                         subject to 1 <= NTERM6 <= 15.
C
C      XLOW1 - DOUBLE PRECISION - The value such that, if -XLOW1 < x < XLOW1,
C                     then AIRYGI = Gi(0) to machine precision.
C                     The recommended value is   EPS.
C
C      XHIGH1 - DOUBLE PRECISION - The value such that, if x > XHIGH1, then
C                      AIRYGI = 1/(Pi*x) to machine precision.
C                      Also used for error test - see above.
C                      The recommended value is
C                          cube root( 2/EPS ).
C
C      XHIGH2 - DOUBLE PRECISION - The value above which AIRYGI = 0.0.
C                      The recommended value is 
C                          1/(Pi*XMIN).
C
C      XHIGH3 - DOUBLE PRECISION - The value such that, if x < XHIGH3,
C                      then the Chebyshev expansions for the
C                      asymptotic form of Bi(x) are not needed.
C                      The recommended value is
C                          -8 * cube root( 2/EPSNEG ).
C
C      For values of EPS, EPSNEG, and XMIN refer to the file
C      MACHCON.TXT.
C
C      The machine-arithmetic constants are given in DATA
C      statements.
C
C
C   INTRINSIC FUNCTIONS USED:
C                             COS , SIN , SQRT
C
C
C   OTHER MISCFUN SUBROUTINES USED:
C
C          CHEVAL , ERRPRN
C
C
C   AUTHOR:
C          Dr. Allan J. Macleod,
C          Dept. of Mathematics and Statistics,
C          University of Paisley,
C          High St.,
C          Paisley,
C          SCOTLAND.
C
C          (e-mail: macl_ms0@paisley.ac.uk)
C
C
C   LATEST UPDATE:
C                 22 NOVEMBER, 1995.
C                          
      INTEGER NTERM1,NTERM2,NTERM3,NTERM4,NTERM5,NTERM6
      DOUBLE PRECISION ARGIP1(0:30),ARGIP2(0:29),ARGIN1(0:42),
     1     ARBIN1(0:10),ARBIN2(0:11),ARHIN1(0:15),
     2     BI,CHEB1,CHEB2,CHEVAL,COSZ,FIVE,FIVE14,FOUR,
     3     GIZERO,MINATE,NINE,ONE,ONEBPI,ONEHUN,ONE76,ONE024,PIBY4,
     4     RTPIIN,SEVEN,SEVEN2,SINZ,T,TEMP,THREE,TWELHU,TWENT8,
     5     X,XCUBE,XHIGH1,XHIGH2,XHIGH3,XLOW1,XMINUS,
     6     XVALUE,ZERO,ZETA
      CHARACTER FNNAME*6,ERRMSG*46
      DATA FNNAME/'AIRYGI'/
      DATA ERRMSG/'ARGUMENT TOO NEGATIVE FOR ACCURATE COMPUTATION'/
      DATA ARGIP1(0)/  0.26585 77079 50227 45082  D    0/
      DATA ARGIP1(1)/ -0.10500 33309 75019 22907  D    0/
      DATA ARGIP1(2)/  0.84134 74753 28454 492    D   -2/
      DATA ARGIP1(3)/  0.20210 67387 81343 9541   D   -1/
      DATA ARGIP1(4)/ -0.15595 76113 86355 2234   D   -1/
      DATA ARGIP1(5)/  0.56434 29390 43256 481    D   -2/
      DATA ARGIP1(6)/ -0.59776 84482 66558 09     D   -3/
      DATA ARGIP1(7)/ -0.42833 85026 48677 28     D   -3/
      DATA ARGIP1(8)/  0.22605 66238 09090 27     D   -3/
      DATA ARGIP1(9)/ -0.36083 32945 59226 0      D   -4/
      DATA ARGIP1(10)/-0.78551 89887 88901        D   -5/
      DATA ARGIP1(11)/ 0.47325 24807 46370        D   -5/
      DATA ARGIP1(12)/-0.59743 51397 7694         D   -6/
      DATA ARGIP1(13)/-0.15917 60916 5602         D   -6/
      DATA ARGIP1(14)/ 0.63361 29065 570          D   -7/
      DATA ARGIP1(15)/-0.27609 02326 48           D   -8/
      DATA ARGIP1(16)/-0.25606 41540 85           D   -8/
      DATA ARGIP1(17)/ 0.47798 67685 6            D   -9/
      DATA ARGIP1(18)/ 0.44881 31863              D  -10/
      DATA ARGIP1(19)/-0.23465 08882              D  -10/
      DATA ARGIP1(20)/ 0.76839 085                D  -12/
      DATA ARGIP1(21)/ 0.73227 985                D  -12/
      DATA ARGIP1(22)/-0.85136 87                 D  -13/
      DATA ARGIP1(23)/-0.16302 01                 D  -13/
      DATA ARGIP1(24)/ 0.35676 9                  D  -14/
      DATA ARGIP1(25)/ 0.25001                    D  -15/
      DATA ARGIP1(26)/-0.10859                    D  -15/
      DATA ARGIP1(27)/-0.158                      D  -17/
      DATA ARGIP1(28)/ 0.275                      D  -17/
      DATA ARGIP1(29)/-0.5                        D  -19/
      DATA ARGIP1(30)/-0.6                        D  -19/
      DATA ARGIP2(0)/  2.00473 71227 58014 86391  D    0/
      DATA ARGIP2(1)/  0.29418 41393 64406 724    D   -2/
      DATA ARGIP2(2)/  0.71369 24900 63401 67     D   -3/
      DATA ARGIP2(3)/  0.17526 56343 05022 67     D   -3/
      DATA ARGIP2(4)/  0.43591 82094 02988 2      D   -4/
      DATA ARGIP2(5)/  0.10926 26947 60430 7      D   -4/
      DATA ARGIP2(6)/  0.27238 24183 99029        D   -5/
      DATA ARGIP2(7)/  0.66230 90094 7687         D   -6/
      DATA ARGIP2(8)/  0.15425 32337 0315         D   -6/
      DATA ARGIP2(9)/  0.34184 65242 306          D   -7/
      DATA ARGIP2(10)/ 0.72815 77248 94           D   -8/
      DATA ARGIP2(11)/ 0.15158 85254 52           D   -8/
      DATA ARGIP2(12)/ 0.30940 04803 9            D   -9/
      DATA ARGIP2(13)/ 0.61496 72614              D  -10/
      DATA ARGIP2(14)/ 0.12028 77045              D  -10/
      DATA ARGIP2(15)/ 0.23369 0586               D  -11/
      DATA ARGIP2(16)/ 0.43778 068                D  -12/
      DATA ARGIP2(17)/ 0.79964 47                 D  -13/
      DATA ARGIP2(18)/ 0.14940 75                 D  -13/
      DATA ARGIP2(19)/ 0.24679 0                  D  -14/
      DATA ARGIP2(20)/ 0.37672                    D  -15/
      DATA ARGIP2(21)/ 0.7701                     D  -16/
      DATA ARGIP2(22)/ 0.354                      D  -17/
      DATA ARGIP2(23)/-0.49                       D  -18/
      DATA ARGIP2(24)/ 0.62                       D  -18/
      DATA ARGIP2(25)/-0.40                       D  -18/
      DATA ARGIP2(26)/-0.1                        D  -19/
      DATA ARGIP2(27)/ 0.2                        D  -19/
      DATA ARGIP2(28)/-0.3                        D  -19/
      DATA ARGIP2(29)/ 0.1                        D  -19/
      DATA ARGIN1(0)/ -0.20118 96505 67320 89130  D    0/
      DATA ARGIN1(1)/ -0.72441 75303 32453 0499   D   -1/
      DATA ARGIN1(2)/  0.45050 18923 89478 0120   D   -1/
      DATA ARGIN1(3)/ -0.24221 37112 20787 91099  D    0/
      DATA ARGIN1(4)/  0.27178 84964 36167 8294   D   -1/
      DATA ARGIN1(5)/ -0.57293 21004 81817 9697   D   -1/
      DATA ARGIN1(6)/ -0.18382 10786 03377 63587  D    0/
      DATA ARGIN1(7)/  0.77515 46082 14947 5511   D   -1/
      DATA ARGIN1(8)/  0.18386 56473 39275 60387  D    0/
      DATA ARGIN1(9)/  0.29215 04250 18556 7173   D   -1/
      DATA ARGIN1(10)/-0.61422 94846 78801 8811   D   -1/
      DATA ARGIN1(11)/-0.29993 12505 79461 6238   D   -1/
      DATA ARGIN1(12)/ 0.58593 71183 27706 636    D   -2/
      DATA ARGIN1(13)/ 0.82222 16584 97402 529    D   -2/
      DATA ARGIN1(14)/ 0.13257 98171 66846 893    D   -2/
      DATA ARGIN1(15)/-0.96248 31076 65651 26     D   -3/
      DATA ARGIN1(16)/-0.45065 51599 82118 07     D   -3/
      DATA ARGIN1(17)/ 0.77242 34743 25474        D   -5/
      DATA ARGIN1(18)/ 0.54818 74134 75805 2      D   -4/
      DATA ARGIN1(19)/ 0.12458 98039 74287 6      D   -4/
      DATA ARGIN1(20)/-0.24619 68910 92083        D   -5/
      DATA ARGIN1(21)/-0.16915 41835 45285        D   -5/
      DATA ARGIN1(22)/-0.16769 15316 9442         D   -6/
      DATA ARGIN1(23)/ 0.96365 09337 672          D   -7/
      DATA ARGIN1(24)/ 0.32533 14928 030          D   -7/
      DATA ARGIN1(25)/ 0.50918 04231              D  -10/
      DATA ARGIN1(26)/-0.20918 04535 53           D   -8/
      DATA ARGIN1(27)/-0.41237 38787 0            D   -9/
      DATA ARGIN1(28)/ 0.41633 38253              D  -10/
      DATA ARGIN1(29)/ 0.30325 32117              D  -10/
      DATA ARGIN1(30)/ 0.34058 0529               D  -11/
      DATA ARGIN1(31)/-0.88444 592                D  -12/
      DATA ARGIN1(32)/-0.31639 612                D  -12/
      DATA ARGIN1(33)/-0.15050 76                 D  -13/
      DATA ARGIN1(34)/ 0.11041 48                 D  -13/
      DATA ARGIN1(35)/ 0.24650 8                  D  -14/
      DATA ARGIN1(36)/-0.3107                     D  -16/
      DATA ARGIN1(37)/-0.9851                     D  -16/
      DATA ARGIN1(38)/-0.1453                     D  -16/
      DATA ARGIN1(39)/ 0.118                      D  -17/
      DATA ARGIN1(40)/ 0.67                       D  -18/
      DATA ARGIN1(41)/ 0.6                        D  -19/
      DATA ARGIN1(42)/-0.1                        D  -19/
      DATA ARBIN1/1.99983 76358 35861 55980  D    0,
     1           -0.81046 60923 66941 8      D   -4,
     2            0.13475 66598 4689         D   -6,
     3           -0.70855 84714 3            D   -9,
     4            0.74818 4187               D  -11,
     5           -0.12902 774                D  -12,
     6            0.32250 4                  D  -14,
     7           -0.10809                    D  -15,
     8            0.460                      D  -17,
     9           -0.24                       D  -18,
     X            0.1                        D  -19/
      DATA ARBIN2/0.13872 35645 38791 20276  D    0,
     1           -0.82392 86225 55822 8      D   -4,
     2            0.26720 91950 9866         D   -6,
     3           -0.20742 36853 68           D   -8,
     4            0.28733 92593              D  -10,
     5           -0.60873 521                D  -12,
     6            0.17924 89                 D  -13,
     7           -0.68760                    D  -15,
     8            0.3280                     D  -16,
     9           -0.188                      D  -17,
     X            0.13                       D  -18,
     1           -0.1                        D  -19/
      DATA ARHIN1/1.99647 72039 97796 50525  D    0,
     1           -0.18756 37794 07173 213    D   -2,
     2           -0.12186 47089 77873 39     D   -3,
     3           -0.81402 16096 59287        D   -5,
     4           -0.55050 92595 3537         D   -6,
     5           -0.37630 08043 303          D   -7,
     6           -0.25885 83623 65           D   -8,
     7           -0.17931 82926 5            D   -9,
     8           -0.12459 16873              D  -10,
     9           -0.87171 247                D  -12,
     X           -0.60849 43                 D  -13,
     1           -0.43117 8                  D  -14,
     2           -0.29787                    D  -15,
     3           -0.2210                     D  -16,
     4           -0.136                      D  -17,
     5           -0.14                       D  -18/
      DATA ZERO,ONE,THREE,FOUR/ 0.0 D 0 , 1.0 D 0 , 3.0 D 0 , 4.0 D 0 /
      DATA FIVE,SEVEN,MINATE/ 5.0 D 0 , 7.0 D 0 , -8.0 D 0 /
      DATA NINE,TWENT8,SEVEN2/ 9.0 D 0 , 28.0 D 0 , 72.0 D 0 /
      DATA ONEHUN,ONE76,FIVE14/ 100.0 D 0 , 176.0 D 0 , 514.0 D 0 /
      DATA ONE024,TWELHU/ 1024.0 D 0 , 1200.0 D 0 /
      DATA GIZERO/0.20497 55424 82000 24505 D 0/
      DATA ONEBPI/0.31830 98861 83790 67154 D 0/
      DATA PIBY4/0.78539 81633 97448 30962 D 0/
      DATA RTPIIN/0.56418 95835 47756 28695 D 0/
C
C   Machine-dependent constants (suitable for IEEE machines)
C
      DATA NTERM1,NTERM2,NTERM3/28,23,39/
      DATA NTERM4,NTERM5,NTERM6/9,10,14/
      DATA XLOW1,XHIGH1/2.22045D-16,208063.8307D0/
      DATA XHIGH2,XHIGH3/0.14274D308,-2097152.0D0/
C
C   Start computation
C
      X = XVALUE
C
C   Error test
C
      IF ( X .LT. -XHIGH1*XHIGH1 ) THEN
         CALL ERRPRN(FNNAME,ERRMSG)
         AIRYGI = ZERO
         RETURN
      ENDIF
C
C   Code for x >= 0.0
C
      IF ( X .GE. ZERO ) THEN
         IF ( X .LE. SEVEN ) THEN
            IF ( X .LT. XLOW1 ) THEN
               AIRYGI = GIZERO
            ELSE
               T = ( NINE * X - TWENT8 ) / ( X + TWENT8 )
               AIRYGI = CHEVAL ( NTERM1 , ARGIP1 , T )
            ENDIF
         ELSE
            IF ( X .GT. XHIGH1 ) THEN
               IF ( X .GT. XHIGH2 ) THEN
                  AIRYGI = ZERO
               ELSE
                  AIRYGI = ONEBPI/X
               ENDIF
            ELSE
               XCUBE = X * X * X
               T = ( TWELHU - XCUBE ) / ( FIVE14 + XCUBE )
               AIRYGI = ONEBPI * CHEVAL(NTERM2,ARGIP2,T) / X
            ENDIF
         ENDIF
      ELSE
C
C   Code for x < 0.0
C
         IF ( X .GE. MINATE ) THEN
            IF ( X .GT. -XLOW1 ) THEN
               AIRYGI = GIZERO
            ELSE
               T = -( X + FOUR ) / FOUR
               AIRYGI = CHEVAL(NTERM3,ARGIN1,T)
            ENDIF
         ELSE
            XMINUS = -X
            T = XMINUS * SQRT(XMINUS)
            ZETA = ( T + T ) / THREE
            TEMP = RTPIIN / SQRT(SQRT(XMINUS))
            COSZ = COS ( ZETA + PIBY4 )
            SINZ = SIN ( ZETA + PIBY4 ) / ZETA
            XCUBE = X * X * X
            IF ( X .GT. XHIGH3 ) THEN
               T = - ( ONE024 / ( XCUBE ) + ONE )
               CHEB1 = CHEVAL(NTERM4,ARBIN1,T)
               CHEB2 = CHEVAL(NTERM5,ARBIN2,T)
               BI = ( COSZ * CHEB1 + SINZ * CHEB2 ) * TEMP
            ELSE
               BI = ( COSZ + SINZ * FIVE / SEVEN2 ) * TEMP
            ENDIF
            T = ( XCUBE + TWELHU ) / ( ONE76 - XCUBE )
            AIRYGI = BI + CHEVAL(NTERM6,ARHIN1,T) * ONEBPI / X
         ENDIF
      ENDIF
      RETURN
      END

      DOUBLE PRECISION FUNCTION AIRYHI(XVALUE)

c*********************************************************************72
c
cc AIRYHI computes the modified Airy function Hi(x).
c
C   DESCRIPTION:
C
C      This subroutine computes the modified Airy function Hi(x),
C      defined as
C
C         AIRYHI(x) = [ Integral{0 to infinity} exp(x*t-t^3/3) dt ] / pi
C
C      The approximation uses Chebyshev expansions with the coefficients 
C      given to 20 decimal places.
C
C      This subroutine is set up to work on IEEE machines.
C      For other machines, you should retrieve the code
C      from the general MISCFUN archive.
C
C
C   ERROR RETURNS:
C
C      If x > XHIGH1 (see below for definition of XHIGH1), then
C      the asymptotic expansion of Hi(x) will cause an overflow.
C      An error message is printed and the code returns the largest
C      floating-pt number as the result.
C 
C
C   MACHINE-DEPENDENT CONSTANTS:
C
C      NTERM1 - INTEGER - The no. of terms to be used from the array
C                         ARHIP. The recommended value is such that
C                                ABS(ARHIP(NTERM1)) < EPS/100
C                         subject to 1 <= NTERM1 <= 31.
C
C      NTERM2 - INTEGER - The no. of terms to be used from the array
C                         ARBIP. The recommended value is such that
C                                ABS(ARBIP(NTERM2)) < EPS/100
C                         subject to 1 <= NTERM2 <= 23.
C
C      NTERM3 - INTEGER - The no. of terms to be used from the array
C                         ARGIP. The recommended value is such that
C                                ABS(ARGIP1(NTERM3)) < EPS/100
C                         subject to 1 <= NTERM3 <= 29.
C
C      NTERM4 - INTEGER - The no. of terms to be used from the array
C                         ARHIN1. The recommended value is such that
C                                ABS(ARHIN1(NTERM4)) < EPS/100
C                         subject to 1 <= NTERM4 <= 21.
C
C      NTERM5 - INTEGER - The no. of terms to be used from the array
C                         ARHIN2. The recommended value is such that
C                                ABS(ARHIN2(NTERM5)) < EPS/100
C                         subject to 1 <= NTERM5 <= 15.
C
C      XLOW1 - DOUBLE PRECISION - The value such that, if -XLOW1 < x < XLOW1,
C                     then AIRYGI = Hi(0) to machine precision.
C                     The recommended value is   EPS.
C
C      XHIGH1 - DOUBLE PRECISION - The value such that, if x > XHIGH1, then
C                      overflow might occur. The recommended value is
C                      computed as follows:
C                           compute Z = 1.5*LOG(XMAX)
C                        XHIGH1 = ( Z + LOG(Z)/4 + LOG(PI)/2 )**(2/3)
C
C      XNEG1 - DOUBLE PRECISION - The value below which AIRYHI = 0.0.
C                     The recommended value is 
C                          -1/(Pi*XMIN).
C
C      XNEG2 - DOUBLE PRECISION - The value such that, if x < XNEG2, then
C                      AIRYHI = -1/(Pi*x) to machine precision.
C                      The recommended value is
C                          -cube root( 2/EPS ).
C
C      XMAX - DOUBLE PRECISION - The largest possible floating-pt. number.
C                    This is the value given to the function
C                    if x > XHIGH1.
C
C      For values of EPS, EPSNEG, XMIN  and XMAX refer to the file
C      MACHCON.TXT.
C
C      The machine-arithmetic constants are given in DATA
C      statements.
C
C
C   INTRINSIC FUNCTIONS USED:
C                            EXP , LOG , SQRT
C
C   OTHER MISCFUN SUBROUTINES USED:
C
C          CHEVAL , ERRPRN
C
C
C   AUTHOR:
C          Dr. Allan J. Macleod,
C          Dept. of Mathematics and Statistics,
C          University of Paisley,
C          High St.,
C          Paisley,
C          SCOTLAND.
C
C          (e-mail: macl_ms0@paisley.ac.uk)
C
C
C   LATEST UPDATE:
C                 23 NOVEMBER, 1995.
C                         
      INTEGER NTERM1,NTERM2,NTERM3,NTERM4,NTERM5
      DOUBLE PRECISION ARHIP(0:31),ARBIP(0:23),ARGIP1(0:29),
     1     ARHIN1(0:21),ARHIN2(0:15),
     2     BI,CHEVAL,FIVE14,FOUR,GI,HIZERO,LNRTPI,
     3     MINATE,ONE,ONEBPI,ONEHUN,ONE76,SEVEN,T,TEMP,
     4     THREE,THRE43,TWELHU,TWELVE,TWO,X,XCUBE,
     5     XHIGH1,XLOW1,XMAX,XNEG1,XNEG2,XVALUE,
     6     ZERO,ZETA
      CHARACTER FNNAME*6,ERRMSG*30
      DATA FNNAME/'AIRYHI'/
      DATA ERRMSG/'ARGUMENT TO FUNCTION TOO LARGE'/
      DATA ARHIP(0)/ 1.24013 56256 17628 31114  D    0/
      DATA ARHIP(1)/ 0.64856 34197 39265 35804  D    0/
      DATA ARHIP(2)/ 0.55236 25259 21149 03246  D    0/
      DATA ARHIP(3)/ 0.20975 12207 38575 66794  D    0/
      DATA ARHIP(4)/ 0.12025 66911 80523 73568  D    0/
      DATA ARHIP(5)/ 0.37682 24931 09539 3785   D   -1/
      DATA ARHIP(6)/ 0.16510 88671 54807 1651   D   -1/
      DATA ARHIP(7)/ 0.45592 27552 11570 993    D   -2/
      DATA ARHIP(8)/ 0.16182 84804 77635 013    D   -2/
      DATA ARHIP(9)/ 0.40841 28250 81266 63     D   -3/
      DATA ARHIP(10)/0.12196 47972 13940 51     D   -3/
      DATA ARHIP(11)/0.28650 64098 65761 0      D   -4/
      DATA ARHIP(12)/0.74222 15564 24344        D   -5/
      DATA ARHIP(13)/0.16353 62319 32831        D   -5/
      DATA ARHIP(14)/0.37713 90818 8749         D   -6/
      DATA ARHIP(15)/0.78158 00336 008          D   -7/
      DATA ARHIP(16)/0.16384 47121 370          D   -7/
      DATA ARHIP(17)/0.31985 76659 92           D   -8/
      DATA ARHIP(18)/0.61933 90530 7            D   -9/
      DATA ARHIP(19)/0.11411 16119 1            D   -9/
      DATA ARHIP(20)/0.20649 23454              D  -10/
      DATA ARHIP(21)/0.36001 8664               D  -11/
      DATA ARHIP(22)/0.61401 849                D  -12/
      DATA ARHIP(23)/0.10162 125                D  -12/
      DATA ARHIP(24)/0.16437 01                 D  -13/
      DATA ARHIP(25)/0.25908 4                  D  -14/
      DATA ARHIP(26)/0.39931                    D  -15/
      DATA ARHIP(27)/0.6014                     D  -16/
      DATA ARHIP(28)/0.886                      D  -17/
      DATA ARHIP(29)/0.128                      D  -17/
      DATA ARHIP(30)/0.18                       D  -18/
      DATA ARHIP(31)/0.3                        D  -19/
      DATA ARBIP(0)/  2.00582 13820 97590 64905  D    0/
      DATA ARBIP(1)/  0.29447 84491 70441 549    D   -2/
      DATA ARBIP(2)/  0.34897 54514 77535 5      D   -4/
      DATA ARBIP(3)/  0.83389 73337 4343         D   -6/
      DATA ARBIP(4)/  0.31362 15471 813          D   -7/
      DATA ARBIP(5)/  0.16786 53060 15           D   -8/
      DATA ARBIP(6)/  0.12217 93405 9            D   -9/
      DATA ARBIP(7)/  0.11915 84139              D  -10/
      DATA ARBIP(8)/  0.15414 2553               D  -11/
      DATA ARBIP(9)/  0.24844 455                D  -12/
      DATA ARBIP(10)/ 0.42130 12                 D  -13/
      DATA ARBIP(11)/ 0.50529 3                  D  -14/
      DATA ARBIP(12)/-0.60032                    D  -15/
      DATA ARBIP(13)/-0.65474                    D  -15/
      DATA ARBIP(14)/-0.22364                    D  -15/
      DATA ARBIP(15)/-0.3015                     D  -16/
      DATA ARBIP(16)/ 0.959                      D  -17/
      DATA ARBIP(17)/ 0.616                      D  -17/
      DATA ARBIP(18)/ 0.97                       D  -18/
      DATA ARBIP(19)/-0.37                       D  -18/
      DATA ARBIP(20)/-0.21                       D  -18/
      DATA ARBIP(21)/-0.1                        D  -19/
      DATA ARBIP(22)/ 0.2                        D  -19/
      DATA ARBIP(23)/ 0.1                        D  -19/
      DATA ARGIP1(0)/  2.00473 71227 58014 86391  D    0/
      DATA ARGIP1(1)/  0.29418 41393 64406 724    D   -2/
      DATA ARGIP1(2)/  0.71369 24900 63401 67     D   -3/
      DATA ARGIP1(3)/  0.17526 56343 05022 67     D   -3/
      DATA ARGIP1(4)/  0.43591 82094 02988 2      D   -4/
      DATA ARGIP1(5)/  0.10926 26947 60430 7      D   -4/
      DATA ARGIP1(6)/  0.27238 24183 99029        D   -5/
      DATA ARGIP1(7)/  0.66230 90094 7687         D   -6/
      DATA ARGIP1(8)/  0.15425 32337 0315         D   -6/
      DATA ARGIP1(9)/  0.34184 65242 306          D   -7/
      DATA ARGIP1(10)/ 0.72815 77248 94           D   -8/
      DATA ARGIP1(11)/ 0.15158 85254 52           D   -8/
      DATA ARGIP1(12)/ 0.30940 04803 9            D   -9/
      DATA ARGIP1(13)/ 0.61496 72614              D  -10/
      DATA ARGIP1(14)/ 0.12028 77045              D  -10/
      DATA ARGIP1(15)/ 0.23369 0586               D  -11/
      DATA ARGIP1(16)/ 0.43778 068                D  -12/
      DATA ARGIP1(17)/ 0.79964 47                 D  -13/
      DATA ARGIP1(18)/ 0.14940 75                 D  -13/
      DATA ARGIP1(19)/ 0.24679 0                  D  -14/
      DATA ARGIP1(20)/ 0.37672                    D  -15/
      DATA ARGIP1(21)/ 0.7701                     D  -16/
      DATA ARGIP1(22)/ 0.354                      D  -17/
      DATA ARGIP1(23)/-0.49                       D  -18/
      DATA ARGIP1(24)/ 0.62                       D  -18/
      DATA ARGIP1(25)/-0.40                       D  -18/
      DATA ARGIP1(26)/-0.1                        D  -19/
      DATA ARGIP1(27)/ 0.2                        D  -19/
      DATA ARGIP1(28)/-0.3                        D  -19/
      DATA ARGIP1(29)/ 0.1                        D  -19/
      DATA ARHIN1(0)/  0.31481 01720 64234 04116  D    0/
      DATA ARHIN1(1)/ -0.16414 49921 65889 64341  D    0/
      DATA ARHIN1(2)/  0.61766 51597 73091 3071   D   -1/
      DATA ARHIN1(3)/ -0.19718 81185 93593 3028   D   -1/
      DATA ARHIN1(4)/  0.53690 28300 23331 343    D   -2/
      DATA ARHIN1(5)/ -0.12497 70684 39663 038    D   -2/
      DATA ARHIN1(6)/  0.24835 51559 69949 33     D   -3/
      DATA ARHIN1(7)/ -0.41870 24096 74663 0      D   -4/
      DATA ARHIN1(8)/  0.59094 54379 79124        D   -5/
      DATA ARHIN1(9)/ -0.68063 54118 4345         D   -6/
      DATA ARHIN1(10)/ 0.60728 97629 164          D   -7/
      DATA ARHIN1(11)/-0.36713 03492 42           D   -8/
      DATA ARHIN1(12)/ 0.70780 17552              D  -10/
      DATA ARHIN1(13)/ 0.11878 94334              D  -10/
      DATA ARHIN1(14)/-0.12089 8723               D  -11/
      DATA ARHIN1(15)/ 0.11896 56                 D  -13/
      DATA ARHIN1(16)/ 0.59412 8                  D  -14/
      DATA ARHIN1(17)/-0.32257                    D  -15/
      DATA ARHIN1(18)/-0.2290                     D  -16/
      DATA ARHIN1(19)/ 0.253                      D  -17/
      DATA ARHIN1(20)/ 0.9                        D  -19/
      DATA ARHIN1(21)/-0.2                        D  -19/
      DATA ARHIN2/1.99647 72039 97796 50525  D    0,
     1           -0.18756 37794 07173 213    D   -2,
     2           -0.12186 47089 77873 39     D   -3,
     3           -0.81402 16096 59287        D   -5,
     4           -0.55050 92595 3537         D   -6,
     5           -0.37630 08043 303          D   -7,
     6           -0.25885 83623 65           D   -8,
     7           -0.17931 82926 5            D   -9,
     8           -0.12459 16873              D  -10,
     9           -0.87171 247                D  -12,
     X           -0.60849 43                 D  -13,
     1           -0.43117 8                  D  -14,
     2           -0.29787                    D  -15,
     3           -0.2210                     D  -16,
     4           -0.136                      D  -17,
     5           -0.14                       D  -18/
      DATA ZERO,ONE,TWO/ 0.0 D 0 , 1.0 D 0 , 2.0 D 0/
      DATA THREE,FOUR,SEVEN/ 3.0 D 0 , 4.0 D 0 , 7.0 D 0 /
      DATA MINATE,TWELVE,ONE76/ -8.0 D 0 , 12.0 D 0 , 176.0 D 0 /
      DATA THRE43,FIVE14,TWELHU/ 343.0 D 0 , 514.0 D 0 , 1200.0 D 0 /
      DATA ONEHUN/100.0 D 0/
      DATA HIZERO/0.40995 10849 64000 49010 D 0/
      DATA LNRTPI/0.57236 49429 24700 08707 D 0/
      DATA ONEBPI/0.31830 98861 83790 67154 D 0/
C
C   Machine-dependent constants (suitable for IEEE machines)
C
      DATA NTERM1,NTERM2,NTERM3,NTERM4,NTERM5/29,17,22,19,14/
      DATA XLOW1,XHIGH1/2.220446D-16,104.4175D0/
      DATA XNEG1,XNEG2,XMAX/-0.14274D308,-208063.831D0,1.79D308/
C
C   Start computation
C
      X = XVALUE
C
C   Error test
C
      IF ( X .GT. XHIGH1 ) THEN
         CALL ERRPRN(FNNAME,ERRMSG)
         AIRYHI = XMAX
         RETURN
      ENDIF
C
C   Code for x >= 0.0
C
      IF ( X .GE. ZERO ) THEN
         IF ( X .LE. SEVEN ) THEN
            IF ( X .LT. XLOW1 ) THEN
               AIRYHI = HIZERO
            ELSE
               T = ( X + X ) / SEVEN - ONE
               TEMP = ( X + X + X ) / TWO
               AIRYHI = EXP(TEMP) * CHEVAL(NTERM1,ARHIP,T)
            ENDIF
         ELSE
            XCUBE = X * X * X
            TEMP = SQRT(XCUBE)
            ZETA = ( TEMP + TEMP ) / THREE
            T = TWO * ( SQRT(THRE43/XCUBE) ) - ONE
            TEMP = CHEVAL(NTERM2,ARBIP,T)
            TEMP = ZETA + LOG(TEMP) - LOG(X) / FOUR - LNRTPI
            BI = EXP(TEMP)
            T = ( TWELHU - XCUBE ) / ( XCUBE + FIVE14 )
            GI = CHEVAL(NTERM3,ARGIP1,T) * ONEBPI / X
            AIRYHI = BI - GI
         ENDIF
      ELSE
C
C   Code for x < 0.0
C
         IF ( X .GE. MINATE ) THEN
            IF ( X .GT. -XLOW1 ) THEN
               AIRYHI = HIZERO
            ELSE
               T = ( FOUR * X + TWELVE ) / ( X - TWELVE )
               AIRYHI = CHEVAL(NTERM4,ARHIN1,T)
            ENDIF
         ELSE
            IF ( X .LT. XNEG1 ) THEN
               AIRYHI = ZERO
            ELSE
               IF ( X .LT. XNEG2 ) THEN
                  TEMP = ONE
               ELSE
                  XCUBE = X * X * X
                  T = ( XCUBE + TWELHU ) / ( ONE76 - XCUBE )
                  TEMP = CHEVAL(NTERM5,ARHIN2,T)
               ENDIF
               AIRYHI = - TEMP * ONEBPI / X
            ENDIF
         ENDIF
      ENDIF
      RETURN
      END 

      DOUBLE PRECISION FUNCTION ATNINT(XVALUE)

c*********************************************************************72
c
cc ATNINT calculates the inverse tangent integral.
c
C DESCRIPTION:
C   
C      The function ATNINT calculates the value of the
C      inverse-tangent integral defined by
C
C           ATNINT(x) = integral 0 to x ( (arctan t)/t ) dt
C
C      The approximation uses Chebyshev series with the coefficients
C      given to an accuracy of 20D.
C
C      This subroutine is set up to work on IEEE machines.
C      For other machines, you should retrieve the code
C      from the general MISCFUN archive.
C
C
C ERROR RETURNS:
C
C   There are no error returns from this program.
C
C
C MACHINE-DEPENDENT CONSTANTS:
C
C   NTERMS - INTEGER - The no. of terms of the array ATNINTT.
C                      The recommended value is such that
C                          ATNINA(NTERMS) < EPS/100   
C
C   XLOW - DOUBLE PRECISION - A bound below which ATNINT(x) = x to machine
C                 precision. The recommended value is
C                     sqrt(EPSNEG/2).
C 
C   XUPPER - DOUBLE PRECISION - A bound on x, above which, to machine precision 
C                   ATNINT(x) = (pi/2)ln x
C                   The recommended value is 1/EPS.
C
C     For values of EPSNEG and EPS for various machine/compiler
C     combinations refer to the text file MACHCON.TXT
C
C     The machine-arithmetic constants are given in DATA
C     statements.
C
C
C INTRINSIC FUNCTIONS USED:
C
C    ABS , LOG
C
C
C   OTHER MISCFUN SUBROUTINES USED:
C
C          CHEVAL
C
C
C AUTHOR: Dr. Allan J. MacLeod,
C         Dept. of Mathematics and Statistics,
C         University of Paisley,
C         High St.,
C         PAISLEY
C         SCOTLAND
C
C         (e-mail  macl_ms0@paisley.ac.uk )
C
C
C LATEST MODIFICATION:  23 NOVEMBER , 1995
C
C
C
      INTEGER IND,NTERMS
      DOUBLE PRECISION ATNINA(0:22),CHEVAL,HALF,ONE,ONEHUN,T,TWOBPI,
     &     X,XLOW,XUPPER,XVALUE,ZERO
      DATA ZERO,HALF,ONE/0.0 D 0 , 0.5 D 0 , 1.0 D 0/
      DATA ONEHUN/100.0 D 0/
      DATA TWOBPI/0.63661 97723 67581 34308 D 0/
      DATA ATNINA(0)/  1.91040 36129 62359 37512  D    0/
      DATA ATNINA(1)/ -0.41763 51437 65674 6940   D   -1/
      DATA ATNINA(2)/  0.27539 25507 86367 434    D   -2/
      DATA ATNINA(3)/ -0.25051 80952 62488 81     D   -3/
      DATA ATNINA(4)/  0.26669 81285 12117 1      D   -4/
      DATA ATNINA(5)/ -0.31189 05141 07001        D   -5/
      DATA ATNINA(6)/  0.38833 85313 2249         D   -6/
      DATA ATNINA(7)/ -0.50572 74584 964          D   -7/
      DATA ATNINA(8)/  0.68122 52829 49           D   -8/
      DATA ATNINA(9)/ -0.94212 56165 4            D   -9/
      DATA ATNINA(10)/ 0.13307 87881 6            D   -9/
      DATA ATNINA(11)/-0.19126 78075              D  -10/
      DATA ATNINA(12)/ 0.27891 2620               D  -11/
      DATA ATNINA(13)/-0.41174 820                D  -12/
      DATA ATNINA(14)/ 0.61429 87                 D  -13/
      DATA ATNINA(15)/-0.92492 9                  D  -14/
      DATA ATNINA(16)/ 0.14038 7                  D  -14/
      DATA ATNINA(17)/-0.21460                    D  -15/
      DATA ATNINA(18)/ 0.3301                     D  -16/
      DATA ATNINA(19)/-0.511                      D  -17/
      DATA ATNINA(20)/ 0.79                       D  -18/
      DATA ATNINA(21)/-0.12                       D  -18/
      DATA ATNINA(22)/ 0.2                        D  -19/
C
C   Machine-dependent constants (suitable for IEEE machines)
C
      DATA NTERMS/19/
      DATA XLOW,XUPPER/7.4505806D-9,4.5036D15/
C
C   Start calculation
C
      IND = 1
      X = XVALUE
      IF ( X .LT. ZERO ) THEN
         X = -X
         IND = -1
      ENDIF
C
C   Code for X < =  1.0
C
      IF ( X .LE. ONE ) THEN
         IF ( X .LT. XLOW ) THEN
            ATNINT = X
         ELSE
            T = X * X
            T =  ( T - HALF ) + ( T - HALF ) 
            ATNINT = X * CHEVAL( NTERMS , ATNINA , T ) 
         ENDIF
      ELSE
C
C   Code for X > 1.0
C
         IF ( X .GT. XUPPER ) THEN
            ATNINT = LOG( X ) / TWOBPI
         ELSE
            T = ONE / ( X * X ) 
            T =  ( T - HALF ) + ( T - HALF ) 
            ATNINT = LOG( X ) / TWOBPI + CHEVAL( NTERMS,ATNINA,T ) / X
         ENDIF
      ENDIF
      IF ( IND .LT. 0 ) ATNINT = - ATNINT
      RETURN
      END

      DOUBLE PRECISION FUNCTION BIRINT(XVALUE)

c*********************************************************************72
c
cc BIRINT calculates the integral of the Airy function Bi.
c
C   DESCRIPTION:
C      This function calculates the integral of the Airy function Bi, defined
C
C          BIRINT(x) = integral{0 to x} Bi(t) dt
C
C      The program uses Chebyshev expansions, the coefficients of which
C      are given to 20 decimal places.
C
C      This subroutine is set up to work on IEEE machines.
C      For other machines, you should retrieve the code
C      from the general MISCFUN archive.
C
C
C   ERROR RETURNS:
C
C      If the function is too large and positive the correct
C      value would overflow. An error message is printed and the
C      program returns the value XMAX.
C
C      If the argument is too large and negative, it is impossible
C      to accurately compute the necessary SIN and COS functions,
C      for the asymptotic expansion.
C      An error message is printed, and the program returns the
C      value 0 (the value at -infinity).
C
C
C   MACHINE-DEPENDENT CONSTANTS:
C
C      NTERM1 - INTEGER - The no. of terms to be used from the array
C                          ABINT1. The recommended value is such that
C                             ABS(ABINT1(NTERM1)) < EPS/100,
C                          subject to 1 <= NTERM1 <= 36.
C
C      NTERM2 - INTEGER - The no. of terms to be used from the array
C                          ABINT2. The recommended value is such that
C                             ABS(ABINT2(NTERM2)) < EPS/100,
C                          subject to 1 <= NTERM2 <= 37.
C
C      NTERM3 - INTEGER - The no. of terms to be used from the array
C                          ABINT3. The recommended value is such that
C                             ABS(ABINT3(NTERM3)) < EPS/100,
C                          subject to 1 <= NTERM3 <= 37.
C 
C      NTERM4 - INTEGER - The no. of terms to be used from the array
C                          ABINT4. The recommended value is such that
C                             ABS(ABINT4(NTERM4)) < EPS/100,
C                          subject to 1 <= NTERM4 <= 20.
C
C      NTERM5 - INTEGER - The no. of terms to be used from the array
C                          ABINT5. The recommended value is such that
C                             ABS(ABINT5(NTERM5)) < EPS/100,
C                          subject to 1 <= NTERM5 <= 20.
C
C      XLOW1 - DOUBLE PRECISION - The value such that, if |x| < XLOW1,
C                          BIRINT(x) = x * Bi(0) 
C                     to machine precision. The recommended value is
C                          2 * EPSNEG.
C
C      XHIGH1 - DOUBLE PRECISION - The value such that, if x > XHIGH1,
C                      the function value would overflow.
C                      The recommended value is computed as
C                          z = ln(XMAX) + 0.5ln(ln(XMAX)),
C                          XHIGH1 = (3z/2)^(2/3)
C
C      XNEG1 - DOUBLE PRECISION - The value such that, if x < XNEG1,
C                     the trigonometric functions in the asymptotic
C                     expansion cannot be calculated accurately.
C                     The recommended value is
C                          -(1/((EPS)**2/3))
C
C      XMAX - DOUBLE PRECISION - The value of the largest positive floating-pt
C                    number. Used in giving a value to the function
C                    if x > XHIGH1.
C
C      For values of EPS, EPSNEG, and XMAX see the file MACHCON.TXT.
C
C
C      The machine-arithmetic constants are given in DATA
C      statements.
C
C
C   INTRINSIC FUNCTIONS USED:
C                            COS, EXP, LOG, SIN, SQRT
C
C
C   OTHER MISCFUN SUBROUTINES USED:
C
C          CHEVAL , ERRPRN
C
C
C   AUTHOR: Dr. Allan J. MacLeod,
C           Dept. of Mathematics and Statistics,
C           Univ. of Paisley,
C           High St.,
C           Paisley,
C           SCOTLAND.
C           PA1 2BE
C 
C           (e-mail: macl_ms0@paisley.ac.uk )
C
C
C   LATEST REVISION: 23 NOVEMBER, 1995.
C
      INTEGER NTERM1,NTERM2,NTERM3,NTERM4,NTERM5
      DOUBLE PRECISION ABINT1(0:36),ABINT2(0:37),ABINT3(0:37),
     1     ABINT4(0:20),ABINT5(0:20),
     2     ARG,BIRZER,CHEVAL,EIGHT,FOUR,F1,F2,NINE,NINHUN,
     3     ONE,ONEHUN,ONEPT5,PIBY4,RT2B3P,SIXTEN,SEVEN,T,TEMP,
     4     THREE,THR644,X,XLOW1,XHIGH1,XMAX,XNEG1,XVALUE,
     5     Z,ZERO
      CHARACTER FNNAME*6,ERMSG1*31,ERMSG2*31
      DATA FNNAME/'BIRINT'/
      DATA ERMSG1/'ARGUMENT TOO LARGE AND POSITIVE'/
      DATA ERMSG2/'ARGUMENT TOO LARGE AND NEGATIVE'/
      DATA ABINT1(0)/  0.38683 35244 50385 43350  D    0/
      DATA ABINT1(1)/ -0.88232 13550 88890 8821   D   -1/
      DATA ABINT1(2)/  0.21463 93744 03554 29239  D    0/
      DATA ABINT1(3)/ -0.42053 47375 89131 5126   D   -1/
      DATA ABINT1(4)/  0.59324 22547 49608 6771   D   -1/
      DATA ABINT1(5)/ -0.84078 70811 24270 210    D   -2/
      DATA ABINT1(6)/  0.87182 47727 78487 955    D   -2/
      DATA ABINT1(7)/ -0.12191 60019 96134 55     D   -3/
      DATA ABINT1(8)/  0.44024 82178 60232 34     D   -3/
      DATA ABINT1(9)/  0.27894 68666 63866 78     D   -3/
      DATA ABINT1(10)/-0.70528 04689 78553 7      D   -4/
      DATA ABINT1(11)/ 0.59010 80066 77010 0      D   -4/
      DATA ABINT1(12)/-0.13708 62587 98214 2      D   -4/
      DATA ABINT1(13)/ 0.50596 25737 49073        D   -5/
      DATA ABINT1(14)/-0.51598 83776 6735         D   -6/
      DATA ABINT1(15)/ 0.39751 13123 49           D   -8/
      DATA ABINT1(16)/ 0.95249 85978 055          D   -7/
      DATA ABINT1(17)/-0.36814 35887 321          D   -7/
      DATA ABINT1(18)/ 0.12483 91688 136          D   -7/
      DATA ABINT1(19)/-0.24909 76191 37           D   -8/
      DATA ABINT1(20)/ 0.31775 24555 1            D   -9/
      DATA ABINT1(21)/ 0.54343 65270              D  -10/
      DATA ABINT1(22)/-0.40245 66915              D  -10/
      DATA ABINT1(23)/ 0.13938 55527              D  -10/
      DATA ABINT1(24)/-0.30381 7509               D  -11/
      DATA ABINT1(25)/ 0.40809 511                D  -12/
      DATA ABINT1(26)/ 0.16341 16                 D  -13/
      DATA ABINT1(27)/-0.26838 09                 D  -13/
      DATA ABINT1(28)/ 0.89664 1                  D  -14/
      DATA ABINT1(29)/-0.18308 9                  D  -14/
      DATA ABINT1(30)/ 0.21333                    D  -15/
      DATA ABINT1(31)/ 0.1108                     D  -16/
      DATA ABINT1(32)/-0.1276                     D  -16/
      DATA ABINT1(33)/ 0.363                      D  -17/
      DATA ABINT1(34)/-0.62                       D  -18/
      DATA ABINT1(35)/ 0.5                        D  -19/
      DATA ABINT1(36)/ 0.1                        D  -19/
      DATA ABINT2(0)/  2.04122 07860 25161 35181  D    0/
      DATA ABINT2(1)/  0.21241 33918 62122 1230   D   -1/
      DATA ABINT2(2)/  0.66617 59976 67062 76     D   -3/
      DATA ABINT2(3)/  0.38420 47982 80825 4      D   -4/
      DATA ABINT2(4)/  0.36231 03660 20439        D   -5/
      DATA ABINT2(5)/  0.50351 99011 5074         D   -6/
      DATA ABINT2(6)/  0.79616 48702 253          D   -7/
      DATA ABINT2(7)/  0.71780 84423 36           D   -8/
      DATA ABINT2(8)/ -0.26777 01591 04           D   -8/
      DATA ABINT2(9)/ -0.16848 95146 99           D   -8/
      DATA ABINT2(10)/-0.36811 75725 5            D   -9/
      DATA ABINT2(11)/ 0.47571 28727              D  -10/
      DATA ABINT2(12)/ 0.52636 21945              D  -10/
      DATA ABINT2(13)/ 0.77897 3500               D  -11/
      DATA ABINT2(14)/-0.46054 6143               D  -11/
      DATA ABINT2(15)/-0.18343 3736               D  -11/
      DATA ABINT2(16)/ 0.32191 249                D  -12/
      DATA ABINT2(17)/ 0.29352 060                D  -12/
      DATA ABINT2(18)/-0.16579 35                 D  -13/
      DATA ABINT2(19)/-0.44838 08                 D  -13/
      DATA ABINT2(20)/ 0.27907                    D  -15/
      DATA ABINT2(21)/ 0.71192 1                  D  -14/
      DATA ABINT2(22)/-0.1042                     D  -16/
      DATA ABINT2(23)/-0.11959 1                  D  -14/
      DATA ABINT2(24)/ 0.4606                     D  -16/
      DATA ABINT2(25)/ 0.20884                    D  -15/
      DATA ABINT2(26)/-0.2416                     D  -16/
      DATA ABINT2(27)/-0.3638                     D  -16/
      DATA ABINT2(28)/ 0.863                      D  -17/
      DATA ABINT2(29)/ 0.591                      D  -17/
      DATA ABINT2(30)/-0.256                      D  -17/
      DATA ABINT2(31)/-0.77                       D  -18/
      DATA ABINT2(32)/ 0.66                       D  -18/
      DATA ABINT2(33)/ 0.3                        D  -19/
      DATA ABINT2(34)/-0.15                       D  -18/
      DATA ABINT2(35)/ 0.2                        D  -19/
      DATA ABINT2(36)/ 0.3                        D  -19/
      DATA ABINT2(37)/-0.1                        D  -19/
      DATA ABINT3(0)/  0.31076 96159 86403 49251  D    0/
      DATA ABINT3(1)/ -0.27528 84588 74525 42718  D    0/
      DATA ABINT3(2)/  0.17355 96570 61365 43928  D    0/
      DATA ABINT3(3)/ -0.55440 17909 49284 3130   D   -1/
      DATA ABINT3(4)/ -0.22512 65478 29595 0941   D   -1/
      DATA ABINT3(5)/  0.41073 47447 81252 1894   D   -1/
      DATA ABINT3(6)/  0.98476 12754 64262 480    D   -2/
      DATA ABINT3(7)/ -0.15556 18141 66604 1932   D   -1/
      DATA ABINT3(8)/ -0.56087 18707 30279 234    D   -2/
      DATA ABINT3(9)/  0.24601 77833 22230 475    D   -2/
      DATA ABINT3(10)/ 0.16574 03922 92336 978    D   -2/
      DATA ABINT3(11)/-0.32775 87501 43540 2      D   -4/
      DATA ABINT3(12)/-0.24434 68086 05149 25     D   -3/
      DATA ABINT3(13)/-0.50353 05196 15232 1      D   -4/
      DATA ABINT3(14)/ 0.16302 64722 24785 4      D   -4/
      DATA ABINT3(15)/ 0.85191 40577 80934        D   -5/
      DATA ABINT3(16)/ 0.29790 36300 4664         D   -6/
      DATA ABINT3(17)/-0.64389 70789 6401         D   -6/
      DATA ABINT3(18)/-0.15046 98814 5803         D   -6/
      DATA ABINT3(19)/ 0.15870 13535 823          D   -7/
      DATA ABINT3(20)/ 0.12767 66299 622          D   -7/
      DATA ABINT3(21)/ 0.14057 85341 99           D   -8/
      DATA ABINT3(22)/-0.46564 73974 1            D   -9/
      DATA ABINT3(23)/-0.15682 74879 1            D   -9/
      DATA ABINT3(24)/-0.40389 3560               D  -11/
      DATA ABINT3(25)/ 0.66670 8192               D  -11/
      DATA ABINT3(26)/ 0.12886 9380               D  -11/
      DATA ABINT3(27)/-0.69686 63                 D  -13/
      DATA ABINT3(28)/-0.62543 19                 D  -13/
      DATA ABINT3(29)/-0.71839 2                  D  -14/
      DATA ABINT3(30)/ 0.11529 6                  D  -14/
      DATA ABINT3(31)/ 0.42276                    D  -15/
      DATA ABINT3(32)/ 0.2493                     D  -16/
      DATA ABINT3(33)/-0.971                      D  -17/
      DATA ABINT3(34)/-0.216                      D  -17/
      DATA ABINT3(35)/-0.2                        D  -19/
      DATA ABINT3(36)/ 0.6                        D  -19/
      DATA ABINT3(37)/ 0.1                        D  -19/
      DATA ABINT4(0)/  1.99507 95931 33520 47614  D    0/
      DATA ABINT4(1)/ -0.27373 63759 70692 738    D   -2/
      DATA ABINT4(2)/ -0.30897 11308 12858 50     D   -3/
      DATA ABINT4(3)/ -0.35501 01982 79857 7      D   -4/
      DATA ABINT4(4)/ -0.41217 92715 20133        D   -5/
      DATA ABINT4(5)/ -0.48235 89231 6833         D   -6/
      DATA ABINT4(6)/ -0.56787 30727 927          D   -7/
      DATA ABINT4(7)/ -0.67187 48103 65           D   -8/
      DATA ABINT4(8)/ -0.79811 64985 7            D   -9/
      DATA ABINT4(9)/ -0.95142 71478              D  -10/
      DATA ABINT4(10)/-0.11374 68966              D  -10/
      DATA ABINT4(11)/-0.13635 9969               D  -11/
      DATA ABINT4(12)/-0.16381 418                D  -12/
      DATA ABINT4(13)/-0.19725 75                 D  -13/
      DATA ABINT4(14)/-0.23784 4                  D  -14/
      DATA ABINT4(15)/-0.28752                    D  -15/
      DATA ABINT4(16)/-0.3475                     D  -16/
      DATA ABINT4(17)/-0.422                      D  -17/
      DATA ABINT4(18)/-0.51                       D  -18/
      DATA ABINT4(19)/-0.6                        D  -19/
      DATA ABINT4(20)/-0.1                        D  -19/
      DATA ABINT5(0)/  1.12672 08196 17825 66017  D    0/
      DATA ABINT5(1)/ -0.67140 55675 25561 198    D   -2/
      DATA ABINT5(2)/ -0.69812 91801 78329 69     D   -3/
      DATA ABINT5(3)/ -0.75616 89886 42527 6      D   -4/
      DATA ABINT5(4)/ -0.83498 55745 10207        D   -5/
      DATA ABINT5(5)/ -0.93630 29823 2480         D   -6/
      DATA ABINT5(6)/ -0.10608 55629 6250         D   -6/
      DATA ABINT5(7)/ -0.12131 28916 741          D   -7/
      DATA ABINT5(8)/ -0.13963 11297 65           D   -8/
      DATA ABINT5(9)/ -0.16178 91805 4            D   -9/
      DATA ABINT5(10)/-0.18823 07907              D  -10/
      DATA ABINT5(11)/-0.22027 2985               D  -11/
      DATA ABINT5(12)/-0.25816 189                D  -12/
      DATA ABINT5(13)/-0.30479 64                 D  -13/
      DATA ABINT5(14)/-0.35837 0                  D  -14/
      DATA ABINT5(15)/-0.42831                    D  -15/
      DATA ABINT5(16)/-0.4993                     D  -16/
      DATA ABINT5(17)/-0.617                      D  -17/
      DATA ABINT5(18)/-0.68                       D  -18/
      DATA ABINT5(19)/-0.10                       D  -18/
      DATA ABINT5(20)/-0.1                        D  -19/
      DATA ZERO,ONE,ONEPT5/ 0.0 D 0 , 1.0 D 0 , 1.5 D 0 /
      DATA THREE,FOUR,SEVEN/ 3.0 D 0 , 4.0 D 0 , 7.0 D 0 /
      DATA EIGHT,NINE,SIXTEN/ 8.0 D 0 , 9.0 D 0 , 16.0 D 0 /
      DATA ONEHUN,NINHUN,THR644/100.0 D 0 , 900.0 D 0 , 3644.0 D 0 /
      DATA PIBY4/0.78539 81633 97448 30962 D 0/
      DATA RT2B3P/0.46065 88659 61780 63902 D 0/
      DATA BIRZER/0.61492 66274 46000 73515 D 0/
C
C   Machine-dependent parameters (suitable for IEEE machines)
C
      DATA NTERM1,NTERM2,NTERM3,NTERM4,NTERM5/33,30,34,17,17/
      DATA XLOW1,XHIGH1/2.22044604D-16,104.587632D0/
      DATA XNEG1,XMAX/-2.727134D10,1.79D308/
C
C   Start computation
C
      X = XVALUE
C
C   Error test
C
      IF ( X .GT. XHIGH1 ) THEN
         CALL ERRPRN(FNNAME,ERMSG1)
         BIRINT = XMAX
         RETURN
      ENDIF
      IF ( X .LT. XNEG1 ) THEN
         CALL ERRPRN(FNNAME,ERMSG2)
         BIRINT = ZERO
         RETURN
      ENDIF
C
C   Code for x >= 0.0
C
      IF ( X .GE. ZERO ) THEN
         IF ( X .LT. XLOW1 ) THEN
            BIRINT = BIRZER * X
         ELSE
            IF ( X .LE. EIGHT ) THEN
               T = X / FOUR - ONE
               BIRINT = X * EXP(ONEPT5*X) * CHEVAL(NTERM1,ABINT1,T)
            ELSE
               T = SIXTEN * SQRT(EIGHT/X) / X - ONE
               Z = ( X + X ) * SQRT(X) / THREE
               TEMP = RT2B3P * CHEVAL(NTERM2,ABINT2,T) / SQRT(Z)
               TEMP = Z + LOG(TEMP)
               BIRINT = EXP(TEMP)
            ENDIF
         ENDIF
      ELSE
C
C   Code for x < 0.0
C
         IF ( X .GE. -SEVEN ) THEN
            IF ( X .GT. -XLOW1 ) THEN
               BIRINT = BIRZER * X
            ELSE
               T = - ( X + X ) / SEVEN - ONE
               BIRINT = X * CHEVAL(NTERM3,ABINT3,T)
            ENDIF 
         ELSE
            Z = - ( X + X ) * SQRT(-X) / THREE
            ARG = Z + PIBY4
            TEMP = NINE * Z * Z
            T = (THR644 - TEMP ) / ( NINHUN + TEMP )
            F1 = CHEVAL(NTERM4,ABINT4,T) * SIN(ARG)
            F2 = CHEVAL(NTERM5,ABINT5,T) * COS(ARG) / Z
            BIRINT = ( F2 - F1 ) * RT2B3P / SQRT(Z)
         ENDIF
      ENDIF
      RETURN
      END

      DOUBLE PRECISION FUNCTION CLAUSN(XVALUE)

c*********************************************************************72
c
cc CLAUSN calculates Clausen's integral.
c
C DESCRIPTION:
C
C      This program calculates Clausen's integral defined by
C
C          CLAUSN(x) = integral 0 to x of (-ln(2*sin(t/2))) dt
C
C      The code uses Chebyshev expansions with the coefficients
C      given to 20 decimal places.
C
C      This subroutine is set up to work on IEEE machines.
C      For other machines, you should retrieve the code
C      from the general MISCFUN archive.
C
C
C ERROR RETURNS:
C
C   If |x| is too large it is impossible to reduce the argument
C   to the range [0,2*pi] with any precision. An error message
C   is printed and the program returns the value 0.0
C
C
C MACHINE-DEPENDENT CONSTANTS:
C
C   NTERMS - INTEGER - the no. of terms of the array ACLAUS
C                      to be used. The recommended value is
C                      such that ABS(ACLAUS(NTERMS)) < EPS/100
C                      subject to 1 <= NTERMS <= 15
C  
C   XSMALL - DOUBLE PRECISION - the value below which Cl(x) can be 
C                   approximated by x (1-ln x). The recommended
C                   value is pi*sqrt(EPSNEG/2).
C
C   XHIGH - DOUBLE PRECISION - The value of |x| above which we cannot
C                  reliably reduce the argument to [0,2*pi].
C                  The recommended value is   1/EPS.
C
C     For values of EPS and EPSNEG refer to the file MACHCON.TXT
C
C     The machine-arithmetic constants are given in DATA
C     statements.
C
C
C INTRINSIC FUNCTIONS USED:
C
C   AINT , LOG , SQRT
C
C
C   OTHER MISCFUN SUBROUTINES USED:
C
C          CHEVAL , ERRPRN
C
C
C AUTHOR:  Dr. Allan J. MacLeod,
C          Dept. of Mathematics and Statistics,
C          University of Paisley,
C          High St.
C          PAISLEY
C          SCOTLAND
C
C          ( e-mail: macl_ms0@paisley.ac.uk )
C
C
C LATEST MODIFICATION:   23 NOVEMBER , 1995
C
      INTEGER INDX,NTERMS
      DOUBLE PRECISION ACLAUS(0:15),CHEVAL,HALF,ONE,ONEHUN,PI,PISQ,T,
     &     TWOPI,TWOPIA,TWOPIB,X,XHIGH,XSMALL,XVALUE,ZERO
      CHARACTER FNNAME*6,ERRMSG*26
      DATA FNNAME/'CLAUSN'/
      DATA ERRMSG/'ARGUMENT TOO LARGE IN SIZE'/
      DATA ZERO,HALF,ONE/0.0 D 0 , 0.5 D 0 , 1.0 D 0/
      DATA ONEHUN/100.0 D 0/
      DATA PI/3.14159 26535 89793 2385 D 0/
      DATA PISQ/9.86960 44010 89358 6188 D 0/
      DATA TWOPI/6.28318 53071 79586 4769 D 0/
      DATA TWOPIA,TWOPIB/6.28125 D 0 , 0.19353 07179 58647 69253 D -2/
      DATA ACLAUS/2.14269 43637 66688 44709  D    0,
     1            0.72332 42812 21257 9245   D   -1,
     2            0.10164 24750 21151 164    D   -2,
     3            0.32452 50328 53164 5      D   -4,
     4            0.13331 51875 71472        D   -5,
     5            0.62132 40591 653          D   -7,
     6            0.31300 41353 37           D   -8,
     7            0.16635 72305 6            D   -9,
     8            0.91965 9293               D  -11,
     9            0.52400 462                D  -12,
     X            0.30580 40                 D  -13,
     1            0.18196 9                  D  -14,
     2            0.11004                    D  -15,
     3            0.675                      D  -17,
     4            0.42                       D  -18,
     5            0.3                        D  -19/
C
C  Set machine-dependent constants (suitable for IEEE machines)
C
      DATA NTERMS/13/
      DATA XSMALL,XHIGH/2.3406689D-8,4.5036D15/
C
C  Start execution
C
      X = XVALUE
C
C   Error test
C
      IF ( ABS(X) .GT. XHIGH ) THEN
         CALL ERRPRN(FNNAME,ERRMSG)
         CLAUSN = ZERO
         RETURN
      ENDIF
      INDX = 1
      IF ( X .LT. ZERO ) THEN
         X = -X
         INDX = -1
      ENDIF
C
C  Argument reduced using simulated extra precision
C
      IF ( X .GT. TWOPI ) THEN
         T = AINT( X / TWOPI ) 
         X =  ( X - T * TWOPIA ) - T * TWOPIB
      ENDIF
      IF ( X .GT. PI ) THEN
         X = ( TWOPIA - X ) + TWOPIB
         INDX = -INDX
      ENDIF
C
C  Set result to zero if X multiple of PI
C
      IF ( X .EQ. ZERO ) THEN
         CLAUSN = ZERO
         RETURN
      ENDIF
C
C  Code for X < XSMALL
C
      IF ( X .LT. XSMALL ) THEN
         CLAUSN = X * ( ONE - LOG( X ) ) 
      ELSE
C
C  Code for XSMALL < =  X < =  PI
C
         T =  ( X * X ) / PISQ - HALF
         T = T + T
         IF ( T .GT. ONE ) T = ONE
         CLAUSN = X * CHEVAL( NTERMS,ACLAUS,T ) - X * LOG( X ) 
      ENDIF
      IF ( INDX .LT. 0 ) CLAUSN = -CLAUSN
      RETURN
      END

      DOUBLE PRECISION FUNCTION DEBYE1(XVALUE)

c*********************************************************************72
c
cc DEBYE1 calculates the Debye function of order 1.
c
C   DEFINITION:
C
C      This program calculates the Debye function of order 1, defined as
C
C            DEBYE1(x) = [Integral {0 to x} t/(exp(t)-1) dt] / x
C
C      The code uses Chebyshev series whose coefficients
C      are given to 20 decimal places.
C
C      This subroutine is set up to work on IEEE machines.
C      For other machines, you should retrieve the code
C      from the general MISCFUN archive.
C
C
C   ERROR RETURNS:
C
C      If XVALUE < 0.0 an error message is printed and the
C      function returns the value 0.0
C
C
C   MACHINE-DEPENDENT PARAMETERS:
C
C      NTERMS - INTEGER - The no. of elements of the array ADEB1.
C                         The recommended value is such that
C                             ABS(ADEB1(NTERMS)) < EPS/100 , with
C                                   1 <= NTERMS <= 18
C
C      XLOW - DOUBLE PRECISION - The value below which
C                    DEBYE1 = 1 - x/4 + x*x/36 to machine precision.
C                    The recommended value is 
C                        SQRT(8*EPSNEG)
C
C      XUPPER - DOUBLE PRECISION - The value above which 
C                      DEBYE1 = (pi*pi/(6*x)) - exp(-x)(x+1)/x.
C                      The recommended value is
C                          -LOG(2*EPS)
C
C      XLIM - DOUBLE PRECISION - The value above which DEBYE1 = pi*pi/(6*x)
C                    The recommended value is 
C                          -LOG(XMIN)
C
C      For values of EPS, EPSNEG, and XMIN see the file MACHCON.TXT
C
C      The machine-arithmetic constants are given in DATA
C      statements.
C
C
C   INTRINSIC FUNCTIONS USED:
C
C      AINT , EXP , INT , LOG , SQRT
C      
C
C   OTHER MISCFUN SUBROUTINES USED:
C
C          CHEVAL , ERRPRN
C
C
C   AUTHOR:
C          Dr. Allan J. MacLeod,
C          Dept. of Mathematics and Statistics,
C          University of Paisley
C          High St.
C          PAISLEY
C          SCOTLAND
C          PA1 2BE
C
C          (e-mail:  macl_ms0@paisley.ac.uk )
C
C
C   LATEST UPDATE:  29 NOVEMBER, 1995
C
      INTEGER I,NEXP,NTERMS
      DOUBLE PRECISION ADEB1(0:18),CHEVAL,DEBINF,EIGHT,EXPMX,FOUR,HALF,
     &     NINE,ONE,ONEHUN,QUART,RK,SUM,T,THIRT6,X,XK,XLIM,XLOW,
     &     XUPPER,XVALUE,ZERO
      CHARACTER FNNAME*6,ERRMSG*17
      DATA FNNAME/'DEBYE1'/
      DATA ERRMSG/'ARGUMENT NEGATIVE'/
      DATA ZERO,QUART/0.0 D 0 , 0.25 D 0/
      DATA HALF,ONE/0.5 D 0 , 1.0 D 0/
      DATA FOUR,EIGHT/4.0 D 0 , 8.0 D 0/
      DATA NINE,THIRT6,ONEHUN/9.0 D 0 , 36.0 D 0 , 100.0 D 0/
      DATA DEBINF/0.60792 71018 54026 62866 D 0/
      DATA ADEB1/2.40065 97190 38141 01941  D    0,
     1           0.19372 13042 18936 00885  D    0,
     2          -0.62329 12455 48957 703    D   -2,
     3           0.35111 74770 20648 00     D   -3,
     4          -0.22822 24667 01231 0      D   -4,
     5           0.15805 46787 50300        D   -5,
     6          -0.11353 78197 0719         D   -6,
     7           0.83583 36118 75           D   -8,
     8          -0.62644 24787 2            D   -9,
     9           0.47603 34890              D  -10,
     X          -0.36574 1540               D  -11,
     1           0.28354 310                D  -12,
     2          -0.22147 29                 D  -13,
     3           0.17409 2                  D  -14,
     4          -0.13759                    D  -15,
     5           0.1093                     D  -16,
     6          -0.87                       D  -18,
     7           0.7                        D  -19,
     8          -0.1                        D  -19/
C
C   Machine-dependent constants (suitable for IEEE machines)
C
      DATA XLOW,XUPPER,XLIM/0.298023D-7,35.35051D0,708.39642D0/
      DATA NTERMS/15/
C
C   Start computation
C
      X = XVALUE
C
C   Check XVALUE >= 0.0
C 
      IF ( X .LT. ZERO ) THEN
         CALL ERRPRN(FNNAME,ERRMSG)
         DEBYE1 = ZERO
         RETURN
      ENDIF
C
C   Code for x <= 4.0
C 
      IF ( X .LE. FOUR ) THEN
         IF ( X .LT. XLOW ) THEN
            DEBYE1 = ( ( X - NINE ) * X + THIRT6 ) / THIRT6
         ELSE
            T = ( ( X * X / EIGHT ) - HALF ) - HALF
            DEBYE1 = CHEVAL( NTERMS , ADEB1 , T ) - QUART * X
         ENDIF
      ELSE
C
C   Code for x > 4.0
C 
         DEBYE1 = ONE / ( X * DEBINF )
         IF ( X .LT. XLIM ) THEN
            EXPMX = EXP( -X )
            IF ( X .GT. XUPPER ) THEN
               DEBYE1 = DEBYE1 - EXPMX * ( ONE + ONE / X )
            ELSE
               SUM = ZERO
               RK = AINT( XLIM / X )
               NEXP = INT( RK )
               XK = RK * X
               DO 100 I = NEXP,1,-1
                  T =  ( ONE + ONE / XK ) / RK
                  SUM = SUM * EXPMX + T
                  RK = RK - ONE
                  XK = XK - X
 100           CONTINUE
               DEBYE1 = DEBYE1 - SUM * EXPMX
            ENDIF
         ENDIF
      ENDIF
      RETURN
      END

      DOUBLE PRECISION FUNCTION DEBYE2(XVALUE)

c*********************************************************************72
c
cc DEBYE2 calculates the Debye function of order 2.
c
C   DEFINITION:
C
C      This program calculates the Debye function of order 1, defined as
C
C            DEBYE2(x) = 2*[Integral {0 to x} t*t/(exp(t)-1) dt] / (x*x)
C
C      The code uses Chebyshev series whose coefficients
C      are given to 20 decimal places.
C
C      This subroutine is set up to work on IEEE machines.
C      For other machines, you should retrieve the code
C      from the general MISCFUN archive.
C
C
C   ERROR RETURNS:
C
C      If XVALUE < 0.0 an error message is printed and the
C      function returns the value 0.0
C
C
C   MACHINE-DEPENDENT PARAMETERS:
C
C      NTERMS - INTEGER - The no. of elements of the array ADEB2.
C                         The recommended value is such that
C                             ABS(ADEB2(NTERMS)) < EPS/100,
C                         subject to 1 <= NTERMS <= 18.
C
C      XLOW - DOUBLE PRECISION - The value below which
C                    DEBYE2 = 1 - x/3 + x*x/24 to machine precision.
C                    The recommended value is 
C                        SQRT(8*EPSNEG)
C
C      XUPPER - DOUBLE PRECISION - The value above which 
C                      DEBYE2 = (4*zeta(3)/x^2) - 2*exp(-x)(x^2+2x+1)/x^2.
C                      The recommended value is
C                          -LOG(2*EPS)
C
C      XLIM1 - DOUBLE PRECISION - The value above which DEBYE2 = 4*zeta(3)/x^2
C                     The recommended value is 
C                          -LOG(XMIN)
C
C      XLIM2 - DOUBLE PRECISION - The value above which DEBYE2 = 0.0 to machine
C                     precision. The recommended value is
C                           SQRT(4.8/XMIN)
C
C      For values of EPS, EPSNEG, and XMIN see the file MACHCON.TXT
C
C      The machine-arithmetic constants are given in DATA
C      statements.
C
C
C   INTRINSIC FUNCTIONS USED:
C
C      AINT , EXP , INT , LOG , SQRT
C      
C
C   OTHER MISCFUN SUBROUTINES USED:
C
C          CHEVAL , ERRPRN
C
C
C   AUTHOR:
C          Dr. Allan J. MacLeod,
C          Dept. of Mathematics and Statistics,
C          University of Paisley
C          High St.
C          PAISLEY
C          SCOTLAND
C          PA1 2BE
C
C          (e-mail:  macl_ms0@paisley.ac.uk )
C
C
C   LATEST UPDATE:  29 NOVEMBER, 1995
C
      INTEGER I,NEXP,NTERMS
      DOUBLE PRECISION ADEB2(0:18),CHEVAL,DEBINF,EIGHT,EXPMX,FOUR,
     &     HALF,ONE,ONEHUN,RK,SUM,T,THREE,TWENT4,TWO,X,XK,XLIM1,
     &     XLIM2,XLOW,XUPPER,XVALUE,ZERO
      CHARACTER FNNAME*6,ERRMSG*17
      DATA FNNAME/'DEBYE2'/
      DATA ERRMSG/'ARGUMENT NEGATIVE'/
      DATA ZERO,HALF/0.0 D 0 , 0.5 D 0/
      DATA ONE,TWO,THREE/1.0 D 0 , 2.0 D 0 , 3.0 D 0/
      DATA FOUR,EIGHT,TWENT4/4.0 D 0 , 8.0 D 0 , 24.0 D 0/
      DATA ONEHUN/100.0 D 0/
      DATA DEBINF/4.80822 76126 38377 14160 D 0/
      DATA ADEB2/2.59438 10232 57077 02826  D    0,
     1           0.28633 57204 53071 98337  D    0,
     2          -0.10206 26561 58046 7129   D   -1,
     3           0.60491 09775 34684 35     D   -3,
     4          -0.40525 76589 50210 4      D   -4,
     5           0.28633 82632 88107        D   -5,
     6          -0.20863 94303 0651         D   -6,
     7           0.15523 78758 264          D   -7,
     8          -0.11731 28008 66           D   -8,
     9           0.89735 85888              D  -10,
     X          -0.69317 6137               D  -11,
     1           0.53980 568                D  -12,
     2          -0.42324 05                 D  -13,
     3           0.33377 8                  D  -14,
     4          -0.26455                    D  -15,
     5           0.2106                     D  -16,
     6          -0.168                      D  -17,
     7           0.13                       D  -18,
     8          -0.1                        D  -19/
C
C   Machine-dependent constants
C
      DATA XLOW,XUPPER/0.298023D-7,35.35051D0/
      DATA XLIM1,XLIM2/708.39642D0,2.1572317D154/
      DATA NTERMS/17/
C
C   Start computation
C
      X = XVALUE
C
C   Check XVALUE >= 0.0
C 
      IF ( X .LT. ZERO ) THEN
         CALL ERRPRN(FNNAME,ERRMSG)
         DEBYE2 = ZERO
         RETURN
      ENDIF
C
C   Code for x <= 4.0
C 
      IF ( X .LE. FOUR ) THEN
         IF ( X .LT. XLOW ) THEN
            DEBYE2 = ( ( X - EIGHT ) * X + TWENT4 ) / TWENT4
         ELSE
            T = ( ( X * X / EIGHT ) - HALF ) - HALF
            DEBYE2 = CHEVAL ( NTERMS , ADEB2 , T ) - X / THREE
         ENDIF
      ELSE
C
C   Code for x > 4.0
C 
         IF ( X .GT. XLIM2 ) THEN
            DEBYE2 = ZERO
         ELSE
            DEBYE2 = DEBINF / ( X * X )
            IF ( X .LT. XLIM1 ) THEN
               EXPMX = EXP ( -X )
               IF ( X .GT. XUPPER ) THEN
                  SUM = ( ( X + TWO ) * X + TWO ) / ( X * X )
               ELSE
                  SUM = ZERO
                  RK = AINT ( XLIM1 / X )
                  NEXP = INT ( RK )
                  XK = RK * X
                  DO 100 I = NEXP,1,-1
                     T =  ( ONE + TWO / XK + TWO / ( XK*XK ) ) / RK
                     SUM = SUM * EXPMX + T
                     RK = RK - ONE
                     XK = XK - X
 100              CONTINUE
               ENDIF
               DEBYE2 = DEBYE2 - TWO * SUM * EXPMX
            ENDIF
         ENDIF
      ENDIF
      RETURN
      END

      DOUBLE PRECISION FUNCTION DEBYE3(XVALUE)

c*********************************************************************72
c
cc DEBYE3 calculates the Debye function of order 3.
c
C   DEFINITION:
C
C      This program calculates the Debye function of order 3, defined as
C
C            DEBYE3(x) = 3*[Integral {0 to x} t^3/(exp(t)-1) dt] / (x^3)
C
C      The code uses Chebyshev series whose coefficients
C      are given to 20 decimal places.
C
C      This subroutine is set up to work on IEEE machines.
C      For other machines, you should retrieve the code
C      from the general MISCFUN archive.
C
C
C   ERROR RETURNS:
C
C      If XVALUE < 0.0 an error message is printed and the
C      function returns the value 0.0
C
C
C   MACHINE-DEPENDENT PARAMETERS:
C
C      NTERMS - INTEGER - The no. of elements of the array ADEB3.
C                         The recommended value is such that
C                             ABS(ADEB3(NTERMS)) < EPS/100,
C                         subject to 1 <= NTERMS <= 18
C
C      XLOW - DOUBLE PRECISION - The value below which
C                    DEBYE3 = 1 - 3x/8 + x*x/20 to machine precision.
C                    The recommended value is 
C                        SQRT(8*EPSNEG)
C
C      XUPPER - DOUBLE PRECISION - The value above which 
C               DEBYE3 = (18*zeta(4)/x^3) - 3*exp(-x)(x^3+3x^2+6x+6)/x^3.
C                      The recommended value is
C                          -LOG(2*EPS)
C
C      XLIM1 - DOUBLE PRECISION - The value above which DEBYE3 = 18*zeta(4)/x^3
C                     The recommended value is 
C                          -LOG(XMIN)
C
C      XLIM2 - DOUBLE PRECISION - The value above which DEBYE3 = 0.0 to machine
C                     precision. The recommended value is
C                          CUBE ROOT(19/XMIN)
C
C      For values of EPS, EPSNEG, and XMIN see the file MACHCON.TXT
C
C      The machine-arithmetic constants are given in DATA
C      statements.
C
C
C   INTRINSIC FUNCTIONS USED:
C
C      AINT , EXP , INT , LOG , SQRT
C      
C
C   OTHER MISCFUN SUBROUTINES USED:
C
C          CHEVAL , ERRPRN
C
C
C   AUTHOR:
C          Dr. Allan J. MacLeod,
C          Dept. of Mathematics and Statistics,
C          University of Paisley
C          High St.
C          PAISLEY
C          SCOTLAND
C          PA1 2BE
C
C          (e-mail:  macl_ms0@paisley.ac.uk )
C
C
C   LATEST UPDATE:  29 NOVEMBER, 1995
C
      INTEGER I,NEXP,NTERMS
      DOUBLE PRECISION ADEB3(0:18),CHEVAL,DEBINF,EIGHT,EXPMX,FOUR,
     &     HALF,ONE,ONEHUN,PT375,RK,SEVP5,SIX,SUM,T,THREE,TWENTY,X,
     &     XK,XKI,XLIM1,XLIM2,XLOW,XUPPER,XVALUE,ZERO
      CHARACTER FNNAME*6,ERRMSG*17
      DATA FNNAME/'DEBYE3'/
      DATA ERRMSG/'ARGUMENT NEGATIVE'/
      DATA ZERO,PT375/0.0 D 0 , 0.375 D 0/
      DATA HALF,ONE/0.5 D 0 , 1.0 D 0/
      DATA THREE,FOUR,SIX/3.0 D 0 , 4.0 D 0 , 6.0 D 0/
      DATA SEVP5,EIGHT,TWENTY/7.5 D 0 , 8.0 D 0 , 20.0 D 0/
      DATA ONEHUN/100.0 D 0/
      DATA DEBINF/0.51329 91127 34216 75946 D -1/
      DATA ADEB3/2.70773 70683 27440 94526  D    0,
     1           0.34006 81352 11091 75100  D    0,
     2          -0.12945 15018 44408 6863   D   -1,
     3           0.79637 55380 17381 64     D   -3,
     4          -0.54636 00095 90823 8      D   -4,
     5           0.39243 01959 88049        D   -5,
     6          -0.28940 32823 5386         D   -6,
     7           0.21731 76139 625          D   -7,
     8          -0.16542 09994 98           D   -8,
     9           0.12727 96189 2            D   -9,
     X          -0.98796 3459               D  -11,
     1           0.77250 740                D  -12,
     2          -0.60779 72                 D  -13,
     3           0.48075 9                  D  -14,
     4          -0.38204                    D  -15,
     5           0.3048                     D  -16,
     6          -0.244                      D  -17,
     7           0.20                       D  -18,
     8          -0.2                        D  -19/
C
C   Machine-dependent constants
C
      DATA XLOW,XUPPER/0.298023D-7,35.35051D0/
      DATA XLIM1,XLIM2/708.39642D0,0.9487163D103/
      DATA NTERMS/16/
C
C   Start computation
C
      X = XVALUE
C
C   Error test
C 
      IF ( X .LT. ZERO ) THEN
         CALL ERRPRN(FNNAME,ERRMSG)
         DEBYE3 = ZERO
         RETURN
      ENDIF
C
C   Code for x <= 4.0
C 
      IF ( X .LE. FOUR ) THEN
         IF ( X .LT. XLOW ) THEN
            DEBYE3 = ( ( X - SEVP5 ) * X + TWENTY ) / TWENTY
         ELSE
            T = ( ( X * X / EIGHT ) - HALF ) - HALF
            DEBYE3 = CHEVAL ( NTERMS , ADEB3 , T ) - PT375 * X
         ENDIF
      ELSE
C
C   Code for x > 4.0
C 
         IF ( X .GT. XLIM2 ) THEN
            DEBYE3 = ZERO
         ELSE
            DEBYE3 = ONE / ( DEBINF * X * X * X )
            IF ( X .LT. XLIM1 ) THEN
               EXPMX = EXP ( -X )
               IF ( X .GT. XUPPER ) THEN
                  SUM = (((X+THREE)*X+SIX)*X+SIX) / (X*X*X)
               ELSE
                  SUM = ZERO
                  RK = AINT ( XLIM1 / X )
                  NEXP = INT ( RK )
                  XK = RK * X
                  DO 100 I = NEXP,1,-1
                     XKI = ONE / XK
                     T =  (((SIX*XKI+SIX)*XKI+THREE)*XKI+ONE) / RK
                     SUM = SUM * EXPMX + T
                     RK = RK - ONE
                     XK = XK - X
 100              CONTINUE
               ENDIF
               DEBYE3 = DEBYE3 - THREE * SUM * EXPMX
            ENDIF
         ENDIF
      ENDIF
      RETURN
      END

      DOUBLE PRECISION FUNCTION DEBYE4(XVALUE)

c*********************************************************************72
c
cc DEBYE4 calculates the Debye function of order 4.
c
C   DEFINITION:
C
C      This program calculates the Debye function of order 4, defined as
C
C            DEBYE4(x) = 4*[Integral {0 to x} t^4/(exp(t)-1) dt] / (x^4)
C
C      The code uses Chebyshev series whose coefficients
C      are given to 20 decimal places.
C
C      This subroutine is set up to work on IEEE machines.
C      For other machines, you should retrieve the code
C      from the general MISCFUN archive.
C
C
C   ERROR RETURNS:
C
C      If XVALUE < 0.0 an error message is printed and the
C      function returns the value 0.0
C
C
C   MACHINE-DEPENDENT PARAMETERS:
C
C      NTERMS - INTEGER - The no. of elements of the array ADEB4.
C                         The recommended value is such that
C                             ABS(ADEB4(NTERMS)) < EPS/100,
C                         subject to 1 <= NTERMS <= 18
C
C      XLOW - DOUBLE PRECISION - The value below which
C                    DEBYE4 = 1 - 4x/10 + x*x/18 to machine precision.
C                    The recommended value is 
C                        SQRT(8*EPSNEG)
C
C      XUPPER - DOUBLE PRECISION - The value above which 
C               DEBYE4=(96*zeta(5)/x^4)-4*exp(-x)(x^4+4x^2+12x^2+24x+24)/x^4.
C                      The recommended value is
C                          -LOG(2*EPS)
C
C      XLIM1 - DOUBLE PRECISION - The value above which DEBYE4 = 96*zeta(5)/x^4
C                     The recommended value is 
C                          -LOG(XMIN)
C
C      XLIM2 - DOUBLE PRECISION - The value above which DEBYE4 = 0.0 to machine
C                     precision. The recommended value is
C                          FOURTH ROOT(99/XMIN)
C
C      For values of EPS, EPSNEG, and XMIN see the file MACHCON.TXT
C
C      The machine-arithmetic constants are given in DATA
C      statements.
C
C
C   INTRINSIC FUNCTIONS USED:
C
C      AINT , EXP , INT , LOG , SQRT
C      
C
C   OTHER MISCFUN SUBROUTINES USED:
C
C          CHEVAL , ERRPRN
C
C
C   AUTHOR:
C          Dr. Allan J. MacLeod,
C          Dept. of Mathematics and Statistics,
C          University of Paisley
C          High St.
C          PAISLEY
C          SCOTLAND
C          PA1 2BE
C
C          (e-mail:  macl_ms0@paisley.ac.uk )
C
C
C   LATEST UPDATE:  29 NOVEMBER, 1995
C
      INTEGER I,NEXP,NTERMS
      DOUBLE PRECISION ADEB4(0:18),CHEVAL,DEBINF,EIGHT,EIGHTN,EXPMX,
     1     FIVE,FOUR,FORTY5,HALF,ONE,ONEHUN,RK,SUM,T,TWELVE,TWENT4,
     2     TWOPT5,X,XK,XKI,XLIM1,XLIM2,XLOW,XUPPER,XVALUE,ZERO
      CHARACTER FNNAME*6,ERRMSG*17
      DATA FNNAME/'DEBYE4'/
      DATA ERRMSG/'ARGUMENT NEGATIVE'/
      DATA ZERO,HALF,ONE/0.0 D 0 , 0.5 D 0 , 1.0 D 0/
      DATA TWOPT5,FOUR,FIVE/2.5 D 0 , 4.0 D 0 , 5.0 D 0/
      DATA EIGHT,TWELVE,EIGHTN/8.0 D 0 , 12.0 D 0 , 18.0 D 0/
      DATA TWENT4,FORTY5,ONEHUN/24.0 D 0 , 45.0 D 0 , 100.0 D 0/
      DATA DEBINF/99.54506 44937 63512 92781 D 0/
      DATA ADEB4/2.78186 94150 20523 46008  D    0,
     1           0.37497 67835 26892 86364  D    0,
     2          -0.14940 90739 90315 8326   D   -1,
     3           0.94567 98114 37042 74     D   -3,
     4          -0.66132 91613 89325 5      D   -4,
     5           0.48156 32982 14449        D   -5,
     6          -0.35880 83958 7593         D   -6,
     7           0.27160 11874 160          D   -7,
     8          -0.20807 09912 23           D   -8,
     9           0.16093 83869 2            D   -9,
     X          -0.12547 09791              D  -10,
     1           0.98472 647                D  -12,
     2          -0.77723 69                 D  -13,
     3           0.61648 3                  D  -14,
     4          -0.49107                    D  -15,
     5           0.3927                     D  -16,
     6          -0.315                      D  -17,
     7           0.25                       D  -18,
     8          -0.2                        D  -19/
C
C   Machine-dependent constants
C
      DATA XLOW,XUPPER/0.298023D-7,35.35051D0/
      DATA XLIM1,XLIM2/708.39642D0,2.5826924D77/
      DATA NTERMS/16/
C
C   Start computation
C
      X = XVALUE
C
C   Check XVALUE >= 0.0
C 
      IF ( X .LT. ZERO ) THEN
         CALL ERRPRN(FNNAME,ERRMSG)
         DEBYE4 = ZERO
         RETURN
      ENDIF
C
C   Code for x <= 4.0
C 
      IF ( X .LE. FOUR ) THEN
         IF ( X .LT. XLOW ) THEN
            DEBYE4 = ( ( TWOPT5 * X - EIGHTN ) * X + FORTY5 ) / FORTY5
         ELSE
            T = ( ( X * X / EIGHT ) - HALF ) - HALF
            DEBYE4 = CHEVAL ( NTERMS , ADEB4 , T ) - ( X + X ) / FIVE
         ENDIF
      ELSE
C
C   Code for x > 4.0
C 
         IF ( X .GT. XLIM2 ) THEN
            DEBYE4 = ZERO
         ELSE
            T = X * X
            DEBYE4 = ( DEBINF / T ) / T
            IF ( X .LT. XLIM1 ) THEN
               EXPMX = EXP ( -X )
               IF ( X .GT. XUPPER ) THEN
                  SUM = ( ( ( ( X + FOUR ) * X + TWELVE ) * X +
     &                  TWENT4 ) * X + TWENT4 ) / ( X * X * X * X )
               ELSE
                  SUM = ZERO
                  RK = AINT ( XLIM1 / X )
                  NEXP = INT ( RK )
                  XK = RK * X
                  DO 100 I = NEXP,1,-1
                     XKI = ONE / XK
                     T =  ( ( ( ( TWENT4 * XKI + TWENT4 ) * XKI +
     &                    TWELVE ) * XKI + FOUR ) * XKI + ONE ) / RK
                     SUM = SUM * EXPMX + T
                     RK = RK - ONE
                     XK = XK - X
 100              CONTINUE
               ENDIF
               DEBYE4 = DEBYE4 - FOUR * SUM * EXPMX
            ENDIF
         ENDIF
      ENDIF
      RETURN
      END

      DOUBLE PRECISION FUNCTION EXP3(XVALUE)

c*********************************************************************72
c
cc EXP3 calculates the integral of exp(-t^3).
c
C   DESCRIPTION
C
C      This function calculates 
C          
C           EXP3(X) = integral 0 to X  (exp(-t*t*t)) dt
C
C      The code uses Chebyshev expansions, whose coefficients are
C      given to 20 decimal places.
C
C      This subroutine is set up to work on IEEE machines.
C      For other machines, you should retrieve the code
C      from the general MISCFUN archive.
C
C
C   ERROR RETURNS
C     
C      If XVALUE < 0, an error message is printed and the function 
C      returns the value 0.
C
C
C   MACHINE-DEPENDENT CONSTANTS
C
C      NTERM1 - INTEGER - The no. of terms of the array AEXP3,
C                         The recommended value is such that
C                               AEXP3(NTERM1) < EPS/100. 
C
C      NTERM2 - INTEGER - The no. of terms of the array AEXP3A.
C                         The recommended value is such that
C                               AEXP3A(NTERM2) < EPS/100.
C
C      XLOW - DOUBLE PRECISION - The value below which EXP3(X) = X to machine
C                    precision. The recommended value is
C                          cube root(4*EPSNEG)
C
C      XUPPER - DOUBLE PRECISION - The value above which EXP3(X) = 0.89297...
C                      to machine precision. The recommended value is
C                           cube root(-ln(EPSNEG))
C
C      For values of EPS and EPSNEG for various machine/compiler
C      combinations refer to the file MACHCON.TXT.
C
C      The machine-arithmetic constants are given in DATA
C      statements.
C
C
C   INTRINSIC FUNCTIONS USED
C
C      EXP, LOG
C
C
C   OTHER MISCFUN SUBROUTINES USED:
C
C          CHEVAL , ERRPRN
C
C
C   AUTHOR
C
C      DR. ALLAN J. MACLEOD,
C      DEPARTMENT OF MATHEMATICS AND STATISTICS,
C      UNIVERSITY OF PAISLEY,
C      HIGH ST.,
C      PAISLEY
C      SCOTLAND.
C
C      (e-mail  macl_ms0@paisley.ac.uk )
C
C
C   LATEST MODIFICATION:  29 NOVEMBER, 1995
C
C
      INTEGER NTERM1,NTERM2
      DOUBLE PRECISION AEXP3(0:24),AEXP3A(0:24),CHEVAL,FOUR,
     1       FUNINF,HALF,ONE,ONEHUN,SIXTEN,T,THREE,TWO,X,
     2       XLOW,XUPPER,XVALUE,ZERO
      CHARACTER FNNAME*6,ERRMSG*14
      DATA FNNAME/'EXP3  '/
      DATA ERRMSG/'ARGUMENT < 0.0'/
      DATA ZERO,HALF,ONE/0.0 D 0 , 0.5 D 0 , 1.0 D 0/
      DATA TWO,THREE,FOUR/2.0 D 0 , 3.0 D 0 , 4.0 D 0 /
      DATA SIXTEN,ONEHUN/16.0 D 0 , 100.0 D 0/
      DATA FUNINF/0.89297 95115 69249 21122 D 0/
      DATA AEXP3(0)/  1.26919 84142 21126 01434  D    0/
      DATA AEXP3(1)/ -0.24884 64463 84140 98226  D    0/
      DATA AEXP3(2)/  0.80526 22071 72310 4125   D   -1/
      DATA AEXP3(3)/ -0.25772 73325 19683 2934   D   -1/
      DATA AEXP3(4)/  0.75998 78873 07377 429    D   -2/
      DATA AEXP3(5)/ -0.20306 95581 94040 510    D   -2/
      DATA AEXP3(6)/  0.49083 45866 99329 17     D   -3/
      DATA AEXP3(7)/ -0.10768 22391 42020 77     D   -3/
      DATA AEXP3(8)/  0.21551 72626 42898 4      D   -4/
      DATA AEXP3(9)/ -0.39567 05137 38429        D   -5/
      DATA AEXP3(10)/ 0.66992 40933 8956         D   -6/
      DATA AEXP3(11)/-0.10513 21808 0703         D   -6/
      DATA AEXP3(12)/ 0.15362 58019 825          D   -7/
      DATA AEXP3(13)/-0.20990 96036 36           D   -8/
      DATA AEXP3(14)/ 0.26921 09538 1            D   -9/
      DATA AEXP3(15)/-0.32519 52422              D  -10/
      DATA AEXP3(16)/ 0.37114 8157               D  -11/
      DATA AEXP3(17)/-0.40136 518                D  -12/
      DATA AEXP3(18)/ 0.41233 46                 D  -13/
      DATA AEXP3(19)/-0.40337 5                  D  -14/
      DATA AEXP3(20)/ 0.37658                    D  -15/
      DATA AEXP3(21)/-0.3362                     D  -16/
      DATA AEXP3(22)/ 0.288                      D  -17/
      DATA AEXP3(23)/-0.24                       D  -18/
      DATA AEXP3(24)/ 0.2                        D  -19/
      DATA AEXP3A(0)/  1.92704 64955 06827 37293  D    0/
      DATA AEXP3A(1)/ -0.34929 35652 04813 8054   D   -1/
      DATA AEXP3A(2)/  0.14503 38371 89830 093    D   -2/
      DATA AEXP3A(3)/ -0.89253 36718 32790 3      D   -4/
      DATA AEXP3A(4)/  0.70542 39219 11838        D   -5/
      DATA AEXP3A(5)/ -0.66717 27454 7611         D   -6/
      DATA AEXP3A(6)/  0.72426 75899 824          D   -7/
      DATA AEXP3A(7)/ -0.87825 82560 56           D   -8/
      DATA AEXP3A(8)/  0.11672 23442 78           D   -8/
      DATA AEXP3A(9)/ -0.16766 31281 2            D   -9/
      DATA AEXP3A(10)/ 0.25755 01577              D  -10/
      DATA AEXP3A(11)/-0.41957 8881               D  -11/
      DATA AEXP3A(12)/ 0.72010 412                D  -12/
      DATA AEXP3A(13)/-0.12949 055                D  -12/
      DATA AEXP3A(14)/ 0.24287 03                 D  -13/
      DATA AEXP3A(15)/-0.47331 1                  D  -14/
      DATA AEXP3A(16)/ 0.95531                    D  -15/
      DATA AEXP3A(17)/-0.19914                    D  -15/
      DATA AEXP3A(18)/ 0.4277                     D  -16/
      DATA AEXP3A(19)/-0.944                      D  -17/
      DATA AEXP3A(20)/ 0.214                      D  -17/
      DATA AEXP3A(21)/-0.50                       D  -18/
      DATA AEXP3A(22)/ 0.12                       D  -18/
      DATA AEXP3A(23)/-0.3                        D  -19/
      DATA AEXP3A(24)/ 0.1                        D  -19/
C
C   Machine-dependent constants (suitable for IEEE machines)
C
      DATA XLOW,XUPPER/0.762939D-5,3.3243018D0/
      DATA NTERM1,NTERM2/22,20/
C
C   Start calculation
C
      X = XVALUE
C
C   Error test
C
      IF ( X .LT. ZERO ) THEN
         CALL ERRPRN(FNNAME,ERRMSG)
         EXP3 = ZERO
         RETURN
      ENDIF
C
C   Code for XVALUE < =  2
C
      IF ( X .LE. TWO ) THEN
         IF ( X .LT. XLOW ) THEN
            EXP3 = X
         ELSE
            T =  (  ( X * X * X / FOUR ) - HALF ) - HALF
            EXP3 = X * CHEVAL ( NTERM1,AEXP3,T ) 
         ENDIF
      ELSE
C
C   Code for XVALUE > 2
C
         IF ( X .GT. XUPPER ) THEN
            EXP3 = FUNINF
         ELSE
            T = ( ( SIXTEN/ ( X * X * X ) ) - HALF ) - HALF
            T = CHEVAL ( NTERM2,AEXP3A,T ) 
            T = T * EXP ( -X * X * X ) / ( THREE * X * X ) 
            EXP3 = FUNINF - T
         ENDIF
      ENDIF
      RETURN
      END            

      DOUBLE PRECISION FUNCTION GOODST(XVALUE)

c*********************************************************************72
c
cc GOODST calculates the integral of exp(-t^2/(t+x)).
c
C   DESCRIPTION:
C
C      This function calculates the function defined as
C
C        GOODST(x) = {integral 0 to inf} ( exp(-u*u)/(u+x) ) du
C
C      The code uses Chebyshev expansions whose coefficients are
C      given to 20 decimal places.
C
C      This subroutine is set up to work on IEEE machines.
C      For other machines, you should retrieve the code
C      from the general MISCFUN archive.
C
C
C   ERROR RETURNS:
C
C      If XVALUE <= 0.0, an error message is printed, and the
C      code returns the value 0.0.
C
C
C   MACHINE-DEPENDENT CONSTANTS:
C
C      NTERM1 - The no. of terms to be used in the array AGOST.
C                The recommended value is such that
C                    AGOST(NTERM1) < EPS/100,
C
C      NTERM2 - The no. of terms to be used in the array AGOSTA.
C                The recommended value is such that
C                    AGOSTA(NTERM2) < EPS/100,
C
C      XLOW - The value below which f(x) = -(gamma/2) - ln(x)
C             to machine precision. The recommended value is
C                EPSNEG
C
C      XHIGH - The value above which f(x) = sqrt(pi)/(2x) to
C              machine precision. The recommended value is
C                 2 / EPSNEG
C
C      For values of EPS and EPSNEG refer to the file MACHCON.TXT
C
C      The machine-arithmetic constants are given in DATA
C      statements.
C
C
C   INTRINSIC FUNCTIONS USED:
C 
C       EXP , LOG
C
C
C   OTHER MISCFUN SUBROUTINES USED:
C
C          CHEVAL , ERRPRN
C
C
C   AUTHOR:
C
C      Dr. Allan J. MacLeod,
C      Dept. of Mathematics and Statistics,
C      University of Paisley,
C      High St.,
C      Paisley.
C      SCOTLAND.
C
C      (e-mail: macl_ms0@paisley.ac.uk )
C
C
C   LATEST REVISION:
C                   29 NOVEMBER, 1995 
C      
C
      INTEGER NTERM1,NTERM2
      DOUBLE PRECISION AGOST(0:28),AGOSTA(0:23),
     1     CHEVAL,FVAL,GAMBY2,HALF,ONE,ONEHUN,RTPIB2,SIX,
     2     T,TWO,X,XHIGH,XLOW,XVALUE,ZERO
      CHARACTER FNNAME*6,ERRMSG*15
      DATA FNNAME/'GOODST'/
      DATA ERRMSG/'ARGUMENT <= 0.0'/
      DATA ZERO,HALF,ONE/ 0.0 D 0 , 0.5 D 0 , 1.0 D 0 /
      DATA TWO,SIX/ 2.0 D 0 , 6.0 D 0 /
      DATA ONEHUN/100.0 D 0/
      DATA GAMBY2/0.28860 78324 50766 43030 D 0/
      DATA RTPIB2/0.88622 69254 52758 01365 D 0/
      DATA AGOST(0)/  0.63106 56056 03984 46247  D    0/
      DATA AGOST(1)/  0.25051 73779 32167 08827  D    0/
      DATA AGOST(2)/ -0.28466 20597 90189 40757  D    0/
      DATA AGOST(3)/  0.87615 87523 94862 3552   D   -1/
      DATA AGOST(4)/  0.68260 22672 21252 724    D   -2/
      DATA AGOST(5)/ -0.10811 29544 19225 4677   D   -1/
      DATA AGOST(6)/  0.16910 12441 17152 176    D   -2/
      DATA AGOST(7)/  0.50272 98462 26151 86     D   -3/
      DATA AGOST(8)/ -0.18576 68720 41000 84     D   -3/
      DATA AGOST(9)/ -0.42870 36741 68474        D   -5/
      DATA AGOST(10)/ 0.10095 98903 20290 5      D   -4/
      DATA AGOST(11)/-0.86529 91351 7382         D   -6/
      DATA AGOST(12)/-0.34983 87432 0734         D   -6/
      DATA AGOST(13)/ 0.64832 78683 494          D   -7/
      DATA AGOST(14)/ 0.75759 24985 83           D   -8/
      DATA AGOST(15)/-0.27793 54243 62           D   -8/
      DATA AGOST(16)/-0.48302 35135              D  -10/
      DATA AGOST(17)/ 0.86632 21283              D  -10/
      DATA AGOST(18)/-0.39433 9687               D  -11/
      DATA AGOST(19)/-0.20952 9625               D  -11/
      DATA AGOST(20)/ 0.21501 759                D  -12/
      DATA AGOST(21)/ 0.39590 15                 D  -13/
      DATA AGOST(22)/-0.69227 9                  D  -14/
      DATA AGOST(23)/-0.54829                    D  -15/
      DATA AGOST(24)/ 0.17108                    D  -15/
      DATA AGOST(25)/ 0.376                      D  -17/
      DATA AGOST(26)/-0.349                      D  -17/
      DATA AGOST(27)/ 0.7                        D  -19/
      DATA AGOST(28)/ 0.6                        D  -19/
      DATA AGOSTA(0)/  1.81775 46798 47187 58767  D    0/
      DATA AGOSTA(1)/ -0.99211 46570 74409 7467   D   -1/
      DATA AGOSTA(2)/ -0.89405 86452 54819 243    D   -2/
      DATA AGOSTA(3)/ -0.94955 33127 77267 85     D   -3/
      DATA AGOSTA(4)/ -0.10971 37996 67596 65     D   -3/
      DATA AGOSTA(5)/ -0.13466 94539 57859 0      D   -4/
      DATA AGOSTA(6)/ -0.17274 92743 08265        D   -5/
      DATA AGOSTA(7)/ -0.22931 38019 9498         D   -6/
      DATA AGOSTA(8)/ -0.31278 44178 918          D   -7/
      DATA AGOSTA(9)/ -0.43619 79736 71           D   -8/
      DATA AGOSTA(10)/-0.61958 46474 3            D   -9/
      DATA AGOSTA(11)/-0.89379 91276              D  -10/
      DATA AGOSTA(12)/-0.13065 11094              D  -10/
      DATA AGOSTA(13)/-0.19316 6876               D  -11/
      DATA AGOSTA(14)/-0.28844 270                D  -12/
      DATA AGOSTA(15)/-0.43447 96                 D  -13/
      DATA AGOSTA(16)/-0.65951 8                  D  -14/
      DATA AGOSTA(17)/-0.10080 1                  D  -14/
      DATA AGOSTA(18)/-0.15502                    D  -15/
      DATA AGOSTA(19)/-0.2397                     D  -16/
      DATA AGOSTA(20)/-0.373                      D  -17/
      DATA AGOSTA(21)/-0.58                       D  -18/
      DATA AGOSTA(22)/-0.9                        D  -19/
      DATA AGOSTA(23)/-0.1                        D  -19/
C
C   Machine-dependent constants (suitable for IEEE machines)
C
      DATA NTERM1,NTERM2/26,20/
      DATA XLOW,XHIGH/1.11022303D-16,1.80144D16/
C
C   Start computation
C
      X = XVALUE
C
C   Error test
C
      IF ( X .LE. ZERO ) THEN
         CALL ERRPRN(FNNAME,ERRMSG)
         GOODST = ZERO
         RETURN
      ENDIF
C
C   Computation for 0 < x <= 2
C
      IF ( X .LE. TWO ) THEN
         IF ( X .LT. XLOW ) THEN
            GOODST = - GAMBY2 - LOG(X)
         ELSE
            T = ( X - HALF ) - HALF
            GOODST = CHEVAL(NTERM1,AGOST,T) - EXP(-X*X) * LOG(X)   
         ENDIF
      ELSE
C
C   Computation for x > 2
C
         FVAL = RTPIB2 / X
         IF ( X .GT. XHIGH ) THEN
            GOODST = FVAL
         ELSE
            T = ( SIX - X ) / ( TWO + X )
            GOODST = FVAL * CHEVAL(NTERM2,AGOSTA,T)
         ENDIF
      ENDIF
      RETURN
      END

      DOUBLE PRECISION FUNCTION I0INT(XVALUE)

c*********************************************************************72
c
cc I0INT computes the integral of the modified Bessel function I0(X).
c
C   DESCRIPTION:
C      This program computes the integral of the modified Bessel
C      function I0(x) using the definition
C
C         I0INT(x) = {integral 0 to x} I0(t) dt
C
C      The program uses Chebyshev expansions, the coefficients of
C      which are given to 20 decimal places.
C
C      This subroutine is set up to work on IEEE machines.
C      For other machines, you should retrieve the code
C      from the general MISCFUN archive.
C
C
C   ERROR RETURNS:
C
C      If |XVALUE| larger than a certain limit, the value of 
C      I0INT would cause an overflow. If such a situation occurs
C      the programs prints an error message, and returns the 
C      value sign(XVALUE)*XMAX, where XMAX is the largest
C      acceptable floating-pt. value.
C
C
C   MACHINE-DEPENDENT CONSTANTS:
C
C      NTERM1 - The no. of terms to be used from the array ARI01.
C                The recommended value is such that
C                    ABS(ARI01(NTERM1)) < EPS/100
C
C      NTERM2 - The no. of terms to be used from the array ARI0A.
C                The recommended value is such that
C                    ABS(ARI0A(NTERM2)) < EPS/100
C
C      XLOW - The value below which I0INT(x) = x, to machine precision.
C             The recommended value is 
C                  sqrt(12*EPS).
C
C      XHIGH - The value above which overflow will occur. The
C              recommended value is
C                  ln(XMAX) + 0.5*ln(ln(XMAX)) + ln(2).
C
C      For values of EPS and XMAX refer to the file MACHCON.TXT.
C
C      The machine-arithmetic constants are given in DATA
C      statements.
C
C
C   INTRINSIC FUNCTIONS USED:
C
C      EXP , LOG , SQRT
C
C
C   OTHER MISCFUN SUBROUTINES USED:
C
C          CHEVAL , ERRPRN
C
C
C   AUTHOR:
C
C      Dr. Allan J. MacLeod,
C      Dept. of Mathematics and Statistics,
C      University of Paisley,
C      High St.,
C      Paisley,
C      SCOTLAND
C      PA1 2BE
C
C      (e-mail :   macl_ms0@paisley.ac.uk )
C
C
C   LATEST REVISION:
C                   29 NOVEMBER, 1995
C
      INTEGER IND,NTERM1,NTERM2
      DOUBLE PRECISION ARI01(0:28),ARI0A(0:33),
     1     ATEEN,CHEVAL,HALF,LNR2PI,ONEHUN,T,TEMP,THREE,THIRT6,
     2     X,XHIGH,XLOW,XVALUE,ZERO
      CHARACTER FNNAME*6,ERRMSG*26
      DATA FNNAME/'I0INT '/
      DATA ERRMSG/'SIZE OF ARGUMENT TOO LARGE'/
      DATA ZERO,HALF,THREE/ 0.0 D 0 , 0.5 D 0 , 3.0 D 0 /
      DATA ATEEN,THIRT6,ONEHUN/ 18.0 D 0 , 36.0 D 0 , 100.0 D 0/
      DATA LNR2PI/0.91893 85332 04672 74178 D 0/
      DATA ARI01(0)/  0.41227 90692 67815 16801  D    0/
      DATA ARI01(1)/ -0.34336 34515 00815 19562  D    0/
      DATA ARI01(2)/  0.22667 58871 57512 42585  D    0/
      DATA ARI01(3)/ -0.12608 16471 87422 60032  D    0/
      DATA ARI01(4)/  0.60124 84628 77799 0271   D   -1/
      DATA ARI01(5)/ -0.24801 20462 91335 8248   D   -1/
      DATA ARI01(6)/  0.89277 33895 65563 897    D   -2/
      DATA ARI01(7)/ -0.28325 37299 36696 605    D   -2/
      DATA ARI01(8)/  0.79891 33904 17129 94     D   -3/
      DATA ARI01(9)/ -0.20053 93366 09648 90     D   -3/
      DATA ARI01(10)/ 0.44168 16783 01431 3      D   -4/
      DATA ARI01(11)/-0.82237 70422 46068        D   -5/
      DATA ARI01(12)/ 0.12005 97942 19015        D   -5/
      DATA ARI01(13)/-0.11350 86500 4889         D   -6/
      DATA ARI01(14)/ 0.69606 01446 6            D   -9/
      DATA ARI01(15)/ 0.18062 27728 36           D   -8/
      DATA ARI01(16)/-0.26039 48137 0            D   -9/
      DATA ARI01(17)/-0.16618 8103               D  -11/
      DATA ARI01(18)/ 0.51050 0232               D  -11/
      DATA ARI01(19)/-0.41515 879                D  -12/
      DATA ARI01(20)/-0.73681 38                 D  -13/
      DATA ARI01(21)/ 0.12793 23                 D  -13/
      DATA ARI01(22)/ 0.10324 7                  D  -14/
      DATA ARI01(23)/-0.30379                    D  -15/
      DATA ARI01(24)/-0.1789                     D  -16/
      DATA ARI01(25)/ 0.673                      D  -17/
      DATA ARI01(26)/ 0.44                       D  -18/
      DATA ARI01(27)/-0.14                       D  -18/
      DATA ARI01(28)/-0.1                        D  -19/
      DATA ARI0A(0)/  2.03739 65457 11432 87070  D    0/
      DATA ARI0A(1)/  0.19176 31647 50331 0248   D   -1/
      DATA ARI0A(2)/  0.49923 33451 92881 47     D   -3/
      DATA ARI0A(3)/  0.22631 87103 65981 5      D   -4/
      DATA ARI0A(4)/  0.15868 21082 85561        D   -5/
      DATA ARI0A(5)/  0.16507 85563 6318         D   -6/
      DATA ARI0A(6)/  0.23850 58373 640          D   -7/
      DATA ARI0A(7)/  0.39298 51823 04           D   -8/
      DATA ARI0A(8)/  0.46042 71419 9            D   -9/
      DATA ARI0A(9)/ -0.70725 58172              D  -10/
      DATA ARI0A(10)/-0.67471 83961              D  -10/
      DATA ARI0A(11)/-0.20269 62001              D  -10/
      DATA ARI0A(12)/-0.87320 338                D  -12/
      DATA ARI0A(13)/ 0.17552 0014               D  -11/
      DATA ARI0A(14)/ 0.60383 944                D  -12/
      DATA ARI0A(15)/-0.39779 83                 D  -13/
      DATA ARI0A(16)/-0.80490 48                 D  -13/
      DATA ARI0A(17)/-0.11589 55                 D  -13/
      DATA ARI0A(18)/ 0.82731 8                  D  -14/
      DATA ARI0A(19)/ 0.28229 0                  D  -14/
      DATA ARI0A(20)/-0.77667                    D  -15/
      DATA ARI0A(21)/-0.48731                    D  -15/
      DATA ARI0A(22)/ 0.7279                     D  -16/
      DATA ARI0A(23)/ 0.7873                     D  -16/
      DATA ARI0A(24)/-0.785                      D  -17/
      DATA ARI0A(25)/-0.1281                     D  -16/
      DATA ARI0A(26)/ 0.121                      D  -17/
      DATA ARI0A(27)/ 0.214                      D  -17/
      DATA ARI0A(28)/-0.27                       D  -18/
      DATA ARI0A(29)/-0.36                       D  -18/
      DATA ARI0A(30)/ 0.7                        D  -19/
      DATA ARI0A(31)/ 0.6                        D  -19/
      DATA ARI0A(32)/-0.2                        D  -19/
      DATA ARI0A(33)/-0.1                        D  -19/
C
C   Machine-dependent constants (suitable for IEEE machines)
C
      DATA NTERM1,NTERM2/25,27/
      DATA XLOW,XHIGH/0.5161914D-7,713.758339D0/
C
C   Start computation
C
      IND = 1
      X = XVALUE
      IF ( XVALUE .LT. ZERO ) THEN
         IND = -1
         X = -X
      ENDIF
C
C   Error test
C
      IF ( X .GT. XHIGH ) THEN
         CALL ERRPRN(FNNAME,ERRMSG)
         I0INT = EXP ( XHIGH - LNR2PI - HALF * LOG(XHIGH) )
         IF ( IND .EQ. -1 ) I0INT = -I0INT
         RETURN
      ENDIF
C
C   Code for 0 <= !x! <= 18
C
      IF ( X .LE. ATEEN ) THEN
         IF ( X .LT. XLOW ) THEN
            I0INT = X
         ELSE
            T = ( THREE * X - ATEEN ) / ( X + ATEEN )
            I0INT = X * EXP(X) * CHEVAL(NTERM1,ARI01,T)
         ENDIF
      ELSE
C
C   Code for !x! > 18
C
         T = ( THIRT6 / X - HALF ) - HALF
         TEMP = X - HALF*LOG(X) - LNR2PI + LOG(CHEVAL(NTERM2,ARI0A,T))
         I0INT = EXP(TEMP)
      ENDIF
      IF ( IND .EQ. -1 ) I0INT = -I0INT
      RETURN
      END

      DOUBLE PRECISION FUNCTION I0ML0(XVALUE)

c*********************************************************************72
c
cc I0ML0 calculates the difference between the Bessel I0 and Struve L0 functions.
c
C   DESCRIPTION:
C
C      This program calculates the function I0ML0 defined as
C
C                I0ML0(x) = I0(x) - L0(x)
C
C      where I0(x) is the modified Bessel function of the first kind of
C      order 0, and L0(x) is the modified Struve function of order 0.
C
C      The code uses Chebyshev expansions with the coefficients 
C      given to an accuracy of 20D.
C
C      This subroutine is set up to work on IEEE machines.
C      For other machines, you should retrieve the code
C      from the general MISCFUN archive.
C
C
C   ERROR RETURNS:
C
C      The coefficients are only suitable for XVALUE >= 0.0. If
C      XVALUE < 0.0, an error message is printed and the function
C      returns the value 0.0
C
C
C   MACHINE-DEPENDENT PARAMETERS:
C
C      NTERM1 - INTEGER - The number of terms required for the array
C                         AI0L0. The recommended value is such that
C                              ABS(AI0L0(NTERM1)) < EPS/100
C
C      NTERM2 - INTEGER - The number of terms required for the array
C                         AI0L0A. The recommended value is such that
C                              ABS(AI0L0A(NTERM2)) < EPS/100
C
C      XLOW - DOUBLE PRECISION - The value below which I0ML0(x) = 1 to machine
C                    precision. The recommended value is
C                               EPSNEG
C
C      XHIGH - DOUBLE PRECISION - The value above which I0ML0(x) = 2/(pi*x) to 
C                     machine precision. The recommended value is
C                               SQRT(800/EPS) 
C
C      For values of EPS, and EPSNEG see the file MACHCON.TXT
C
C      The machine-arithmetic constants are given in DATA
C      statements.
C
C
C   INTRINSIC FUNCTIONS USED:
C
C      SQRT  
C
C
C   OTHER MISCFUN SUBROUTINES USED:
C
C          CHEVAL , ERRPRN
C
C
C   AUTHOR:
C          Dr. Allan J. MacLeod
C          Dept. of Mathematics and Statistics
C          University of Paisley
C          High St.
C          Paisley
C          SCOTLAND
C          PA1 2BE
C
C          ( e-mail: macl_ms0@paisley.ac.uk ) 
C
C
C   LATEST REVISION:
C                  29 NOVEMBER, 1995
C
      INTEGER NTERM1,NTERM2
      DOUBLE PRECISION AI0L0(0:23),AI0L0A(0:23),ATEHUN,CHEVAL,
     1     FORTY,ONE,ONEHUN,SIX,SIXTEN,T,TWOBPI,TWO88,X,XHIGH,
     2     XLOW,XSQ,XVALUE,ZERO
      CHARACTER FNNAME*6,ERRMSG*14
      DATA FNNAME/'I0ML0 '/
      DATA ERRMSG/'ARGUMENT < 0.0'/
      DATA ZERO,ONE/ 0.0 D 0 , 1.0 D 0 /
      DATA SIX,SIXTEN/ 6.0 D 0 , 16.0 D 0 /
      DATA FORTY,ONEHUN/ 40.0 D 0 , 100.0 D 0 /
      DATA TWO88,ATEHUN/ 288.0 D 0 , 800.0 D 0 /
      DATA TWOBPI/0.63661 97723 67581 34308 D 0/
      DATA AI0L0(0)/  0.52468 73679 14855 99138  D    0/
      DATA AI0L0(1)/ -0.35612 46069 96505 86196  D    0/
      DATA AI0L0(2)/  0.20487 20286 40099 27687  D    0/
      DATA AI0L0(3)/ -0.10418 64052 04026 93629  D    0/
      DATA AI0L0(4)/  0.46342 11095 54842 9228   D   -1/
      DATA AI0L0(5)/ -0.17905 87192 40349 8630   D   -1/
      DATA AI0L0(6)/  0.59796 86954 81143 177    D   -2/
      DATA AI0L0(7)/ -0.17177 75476 93565 429    D   -2/
      DATA AI0L0(8)/  0.42204 65446 91714 22     D   -3/
      DATA AI0L0(9)/ -0.87961 78522 09412 5      D   -4/
      DATA AI0L0(10)/ 0.15354 34234 86922 3      D   -4/
      DATA AI0L0(11)/-0.21978 07695 84743        D   -5/
      DATA AI0L0(12)/ 0.24820 68393 6666         D   -6/
      DATA AI0L0(13)/-0.20327 06035 607          D   -7/
      DATA AI0L0(14)/ 0.90984 19842 1            D   -9/
      DATA AI0L0(15)/ 0.25617 93929              D  -10/
      DATA AI0L0(16)/-0.71060 9790               D  -11/
      DATA AI0L0(17)/ 0.32716 960                D  -12/
      DATA AI0L0(18)/ 0.23002 15                 D  -13/
      DATA AI0L0(19)/-0.29210 9                  D  -14/
      DATA AI0L0(20)/-0.3566                     D  -16/
      DATA AI0L0(21)/ 0.1832                     D  -16/
      DATA AI0L0(22)/-0.10                       D  -18/
      DATA AI0L0(23)/-0.11                       D  -18/
      DATA AI0L0A(0)/ 2.00326 51024 11606 43125  D    0/
      DATA AI0L0A(1)/ 0.19520 68515 76492 081    D   -2/
      DATA AI0L0A(2)/ 0.38239 52356 99083 28     D   -3/
      DATA AI0L0A(3)/ 0.75342 80817 05443 6      D   -4/
      DATA AI0L0A(4)/ 0.14959 57655 89707 8      D   -4/
      DATA AI0L0A(5)/ 0.29994 05312 10557        D   -5/
      DATA AI0L0A(6)/ 0.60769 60482 2459         D   -6/
      DATA AI0L0A(7)/ 0.12399 49554 4506         D   -6/
      DATA AI0L0A(8)/ 0.25232 62552 649          D   -7/
      DATA AI0L0A(9)/ 0.50463 48573 32           D   -8/
      DATA AI0L0A(10)/0.97913 23623 0            D   -9/
      DATA AI0L0A(11)/0.18389 11524 1            D   -9/
      DATA AI0L0A(12)/0.33763 09278              D  -10/
      DATA AI0L0A(13)/0.61117 9703               D  -11/
      DATA AI0L0A(14)/0.10847 2972               D  -11/
      DATA AI0L0A(15)/0.18861 271                D  -12/
      DATA AI0L0A(16)/0.32803 45                 D  -13/
      DATA AI0L0A(17)/0.56564 7                  D  -14/
      DATA AI0L0A(18)/0.93300                    D  -15/
      DATA AI0L0A(19)/0.15881                    D  -15/
      DATA AI0L0A(20)/0.2791                     D  -16/
      DATA AI0L0A(21)/0.389                      D  -17/
      DATA AI0L0A(22)/0.70                       D  -18/
      DATA AI0L0A(23)/0.16                       D  -18/
C
C   MACHINE-DEPENDENT CONSTANTS (suitable for IEEE-arithmetic machines)
C
      DATA NTERM1,NTERM2/21,21/
      DATA XLOW,XHIGH/1.11022303D-16,1.8981253D9/
C
C   Start computation
C
      X = XVALUE
C
C   Error test
C
      IF ( X .LT. ZERO ) THEN
         CALL ERRPRN(FNNAME,ERRMSG)
         I0ML0 = ZERO
         RETURN
      ENDIF
C
C   Code for x <= 16
C
      IF ( X .LE. SIXTEN ) THEN
         IF ( X .LT. XLOW ) THEN
            I0ML0 = ONE 
            RETURN
         ELSE
            T = ( SIX * X - FORTY ) / ( X + FORTY )
            I0ML0 = CHEVAL(NTERM1,AI0L0,T)
            RETURN
         ENDIF
      ELSE
C
C   Code for x > 16
C
         IF ( X .GT. XHIGH ) THEN
            I0ML0 = TWOBPI / X
         ELSE
            XSQ = X * X
            T = ( ATEHUN - XSQ ) / ( TWO88 + XSQ )
            I0ML0 = CHEVAL(NTERM2,AI0L0A,T) * TWOBPI / X
         ENDIF
      ENDIF
      RETURN
      END

      DOUBLE PRECISION FUNCTION I1ML1(XVALUE)

c*********************************************************************72
c
cc I1ML1 calculates the difference between the Bessel I1 and Struve L1 functions.
c
C   DESCRIPTION:
C
C      This program calculates the function I1ML1 defined as
C
C                I1ML1(x) = I1(x) - L1(x)
C
C      where I1(x) is the modified Bessel function of the first kind of
C      order 1, and L1(x) is the modified Struve function of order 1.
C
C      The code uses Chebyshev expansions with the coefficients 
C      given to an accuracy of 20D.
C
C      This subroutine is set up to work on IEEE machines.
C      For other machines, you should retrieve the code
C      from the general MISCFUN archive.
C
C
C   ERROR RETURNS:
C
C      The coefficients are only suitable for XVALUE >= 0.0. If
C      XVALUE < 0.0, an error message is printed and the function
C      returns the value 0.0
C
C
C   MACHINE-DEPENDENT PARAMETERS:
C
C      NTERM1 - INTEGER - The number of terms required for the array
C                         AI1L1. The recommended value is such that
C                              ABS(AI1L1(NTERM1)) < EPS/100
C
C      NTERM2 - INTEGER - The number of terms required for the array
C                         AI1L1A. The recommended value is such that
C                              ABS(AI1L1A(NTERM2)) < EPS/100
C
C      XLOW - DOUBLE PRECISION - The value below which I1ML1(x) = x/2 to machine
C                    precision. The recommended value is
C                               2*EPSNEG
C
C      XHIGH - DOUBLE PRECISION - The value above which I1ML1(x) = 2/pi to 
C                     machine precision. The recommended value is
C                               SQRT(800/EPS) 
C
C      For values of EPS, and EPSNEG see the file MACHCON.TXT
C
C      The machine-arithmetic constants are given in DATA
C      statements.
C
C
C   INTRINSIC FUNCTIONS USED:
C
C      SQRT  
C
C
C   OTHER MISCFUN SUBROUTINES USED:
C
C          CHEVAL , ERRPRN
C
C
C   AUTHOR:
C          Dr. Allan J. MacLeod
C          Dept. of Mathematics and Statistics
C          University of Paisley
C          High St.
C          Paisley
C          SCOTLAND
C          PA1 2BE
C
C          (e-mail: macl_ms0@paisley.ac.uk )
C
C
C   LATEST REVISION:
C                  29 NOVEMBER 1995
C
      INTEGER NTERM1,NTERM2
      DOUBLE PRECISION AI1L1(0:23),AI1L1A(0:25),ATEHUN,CHEVAL,
     1     FORTY,ONE,ONEHUN,SIX,SIXTEN,T,TWO,TWOBPI,TWO88,
     2     X,XHIGH,XLOW,XSQ,XVALUE,ZERO
      CHARACTER FNNAME*6,ERRMSG*14
      DATA FNNAME/'I1ML1 '/
      DATA ERRMSG/'ARGUMENT < 0.0'/
      DATA ZERO,ONE,TWO/ 0.0 D 0 , 1.0 D 0 , 2.0 D 0 /
      DATA SIX,SIXTEN,FORTY/ 6.0 D 0 , 16.0 D 0 , 40.0 D 0 /
      DATA ONEHUN,TWO88,ATEHUN/ 100.0 D 0 , 288.0 D 0 , 800.0 D 0 /
      DATA TWOBPI/0.63661 97723 67581 34308 D 0/
      DATA AI1L1(0)/  0.67536 36906 23505 76137  D    0/
      DATA AI1L1(1)/ -0.38134 97109 72665 59040  D    0/
      DATA AI1L1(2)/  0.17452 17077 51339 43559  D    0/
      DATA AI1L1(3)/ -0.70621 05887 23502 5061   D   -1/
      DATA AI1L1(4)/  0.25173 41413 55880 3702   D   -1/
      DATA AI1L1(5)/ -0.78709 85616 06423 321    D   -2/
      DATA AI1L1(6)/  0.21481 43686 51922 006    D   -2/
      DATA AI1L1(7)/ -0.50862 19971 79062 36     D   -3/
      DATA AI1L1(8)/  0.10362 60828 04423 30     D   -3/
      DATA AI1L1(9)/ -0.17954 47212 05724 7      D   -4/
      DATA AI1L1(10)/ 0.25978 82745 15414        D   -5/
      DATA AI1L1(11)/-0.30442 40632 4667         D   -6/
      DATA AI1L1(12)/ 0.27202 39894 766          D   -7/
      DATA AI1L1(13)/-0.15812 61441 90           D   -8/
      DATA AI1L1(14)/ 0.18162 09172              D  -10/
      DATA AI1L1(15)/ 0.64796 7659               D  -11/
      DATA AI1L1(16)/-0.54113 290                D  -12/
      DATA AI1L1(17)/-0.30831 1                  D  -14/
      DATA AI1L1(18)/ 0.30563 8                  D  -14/
      DATA AI1L1(19)/-0.9717                     D  -16/
      DATA AI1L1(20)/-0.1422                     D  -16/
      DATA AI1L1(21)/ 0.84                       D  -18/
      DATA AI1L1(22)/ 0.7                        D  -19/
      DATA AI1L1(23)/-0.1                        D  -19/
      DATA AI1L1A(0)/  1.99679 36189 67891 36501  D    0/
      DATA AI1L1A(1)/ -0.19066 32614 09686 132    D   -2/
      DATA AI1L1A(2)/ -0.36094 62241 01744 81     D   -3/
      DATA AI1L1A(3)/ -0.68418 47304 59982 0      D   -4/
      DATA AI1L1A(4)/ -0.12990 08228 50942 6      D   -4/
      DATA AI1L1A(5)/ -0.24715 21887 05765        D   -5/
      DATA AI1L1A(6)/ -0.47147 83969 1972         D   -6/
      DATA AI1L1A(7)/ -0.90208 19982 592          D   -7/
      DATA AI1L1A(8)/ -0.17304 58637 504          D   -7/
      DATA AI1L1A(9)/ -0.33232 36701 59           D   -8/
      DATA AI1L1A(10)/-0.63736 42173 5            D   -9/
      DATA AI1L1A(11)/-0.12180 23975 6            D   -9/
      DATA AI1L1A(12)/-0.23173 46832              D  -10/
      DATA AI1L1A(13)/-0.43906 8833               D  -11/
      DATA AI1L1A(14)/-0.82847 110                D  -12/
      DATA AI1L1A(15)/-0.15562 249                D  -12/
      DATA AI1L1A(16)/-0.29131 12                 D  -13/
      DATA AI1L1A(17)/-0.54396 5                  D  -14/
      DATA AI1L1A(18)/-0.10117 7                  D  -14/
      DATA AI1L1A(19)/-0.18767                    D  -15/
      DATA AI1L1A(20)/-0.3484                     D  -16/
      DATA AI1L1A(21)/-0.643                      D  -17/
      DATA AI1L1A(22)/-0.118                      D  -17/
      DATA AI1L1A(23)/-0.22                       D  -18/
      DATA AI1L1A(24)/-0.4                        D  -19/
      DATA AI1L1A(25)/-0.1                        D  -19/
C
C   MACHINE-DEPENDENT CONSTANTS (suitable for IEEE machines)
C
      DATA NTERM1,NTERM2/20,22/
      DATA XLOW,XHIGH/2.22044605D-16,1.8981253D9/
C
C   Start computation
C
      X = XVALUE
C
C   Error test
C
      IF ( X .LT. ZERO ) THEN
         CALL ERRPRN(FNNAME,ERRMSG)
         I1ML1 = ZERO
         RETURN
      ENDIF
C
C   Code for x <= 16
C
      IF ( X .LE. SIXTEN ) THEN
         IF ( X .LT. XLOW ) THEN
            I1ML1 = X / TWO 
            RETURN
         ELSE
            T = ( SIX * X - FORTY ) / ( X + FORTY )
            I1ML1 = CHEVAL(NTERM1,AI1L1,T) * X / TWO
            RETURN
         ENDIF
      ELSE
C
C   Code for x > 16
C
         IF ( X .GT. XHIGH ) THEN
            I1ML1 = TWOBPI 
         ELSE
            XSQ = X * X
            T = ( ATEHUN - XSQ ) / ( TWO88 + XSQ )
            I1ML1 = CHEVAL(NTERM2,AI1L1A,T) * TWOBPI 
         ENDIF
      ENDIF
      RETURN
      END

      DOUBLE PRECISION FUNCTION J0INT(XVALUE)

c*********************************************************************72
c
cc J0INT calculates the integral of the Bessel function J0.
c
C   DESCRIPTION:
C
C      This function calculates the integral of the Bessel
C      function J0, defined as
C
C        J0INT(x) = {integral 0 to x} J0(t) dt
C
C      The code uses Chebyshev expansions whose coefficients are
C      given to 20 decimal places.
C
C      This subroutine is set up to work on IEEE machines.
C      For other machines, you should retrieve the code
C      from the general MISCFUN archive.
C
C
C   ERROR RETURNS:
C
C      If the value of |x| is too large, it is impossible to 
C      accurately compute the trigonometric functions used. An
C      error message is printed, and the function returns the
C      value 1.0.
C
C
C   MACHINE-DEPENDENT CONSTANTS:
C
C      NTERM1 - The no. of terms to be used from the array
C                ARJ01. The recommended value is such that
C                   ABS(ARJ01(NTERM1)) < EPS/100, provided that
C
C      NTERM2 - The no. of terms to be used from the array
C                ARJ0A1. The recommended value is such that
C                   ABS(ARJ0A1(NTERM2)) < EPS/100, provided that
C
C      NTERM3 - The no. of terms to be used from the array
C                ARJ0A2. The recommended value is such that
C                   ABS(ARJ0A2(NTERM3)) < EPS/100, provided that
C
C      XLOW - The value of |x| below which J0INT(x) = x to
C             machine-precision. The recommended value is
C                 sqrt(12*EPSNEG)
C
C      XHIGH - The value of |x| above which it is impossible
C              to calculate (x-pi/4) accurately. The recommended
C              value is      1/EPSNEG
C
C      For values of EPS and EPSNEG for various machine/compiler
C      combinations refer to the file MACHCON.TXT.
C
C      The machine-arithmetic constants are given in DATA
C      statements.
C
C
C   INTRINSIC FUNCTIONS USED:
C
C      COS , SIN , SQRT
C
C
C   OTHER MISCFUN SUBROUTINES USED:
C
C          CHEVAL , ERRPRN
C
C
C   AUTHOR:
C          Dr. Allan J. MacLeod,
C          Dept. of Mathematics and Statistics,
C          University of Paisley,
C          Paisley,
C          SCOTLAND
C          PA1 2BE
C
C          (e-mail:   macl_ms0@paisley.ac.uk )
C
C
C   LATEST REVISION:
C                   11 JANUARY, 1996
C
      INTEGER IND,NTERM1,NTERM2,NTERM3
      DOUBLE PRECISION ARJ01(0:23),ARJ0A1(0:21),ARJ0A2(0:18),
     1     CHEVAL,FIVE12,ONE,ONEHUN,ONE28,PIB41,PIB411,PIB412,
     2     PIB42,RT2BPI,SIXTEN,T,TEMP,TWELVE,X,XHIGH,XLOW,
     3     XMPI4,XVALUE,ZERO
      CHARACTER FNNAME*6,ERRMSG*26
      DATA FNNAME/'J0INT '/
      DATA ERRMSG/'ARGUMENT TOO LARGE IN SIZE'/
      DATA ZERO,ONE/ 0.0 D 0 , 1.0 D 0 /
      DATA TWELVE,SIXTEN/ 12.0 D 0 , 16.0 D 0 /
      DATA ONEHUN,ONE28,FIVE12/ 100.0 D 0 , 128.0 D 0 , 512 D 0 /
      DATA RT2BPI/0.79788 45608 02865 35588 D 0/
      DATA PIB411,PIB412/ 201.0 D 0 , 256.0 D 0/
      DATA PIB42/0.24191 33974 48309 61566 D -3/      
      DATA ARJ01(0)/  0.38179 27932 16901 73518  D    0/
      DATA ARJ01(1)/ -0.21275 63635 05053 21870  D    0/
      DATA ARJ01(2)/  0.16754 21340 72157 94187  D    0/
      DATA ARJ01(3)/ -0.12853 20977 21963 98954  D    0/
      DATA ARJ01(4)/  0.10114 40545 57788 47013  D    0/
      DATA ARJ01(5)/ -0.91007 95343 20156 8859   D   -1/
      DATA ARJ01(6)/  0.64013 45264 65687 3103   D   -1/
      DATA ARJ01(7)/ -0.30669 63029 92675 4312   D   -1/
      DATA ARJ01(8)/  0.10308 36525 32506 4201   D   -1/
      DATA ARJ01(9)/ -0.25567 06503 99956 918    D   -2/
      DATA ARJ01(10)/ 0.48832 75580 57983 04     D   -3/
      DATA ARJ01(11)/-0.74249 35126 03607 7      D   -4/
      DATA ARJ01(12)/ 0.92226 05637 30861        D   -5/
      DATA ARJ01(13)/-0.95522 82830 7083         D   -6/
      DATA ARJ01(14)/ 0.83883 55845 986          D   -7/
      DATA ARJ01(15)/-0.63318 44888 58           D   -8/
      DATA ARJ01(16)/ 0.41560 50422 1            D   -9/
      DATA ARJ01(17)/-0.23955 29307              D  -10/
      DATA ARJ01(18)/ 0.12228 6885               D  -11/
      DATA ARJ01(19)/-0.55697 11                 D  -13/
      DATA ARJ01(20)/ 0.22782 0                  D  -14/
      DATA ARJ01(21)/-0.8417                     D  -16/
      DATA ARJ01(22)/ 0.282                      D  -17/
      DATA ARJ01(23)/-0.9                        D  -19/
      DATA ARJ0A1(0)/  1.24030 13303 75189 70827  D    0/
      DATA ARJ0A1(1)/ -0.47812 53536 32280 693    D   -2/
      DATA ARJ0A1(2)/  0.66131 48891 70667 8      D   -4/
      DATA ARJ0A1(3)/ -0.18604 27404 86349        D   -5/
      DATA ARJ0A1(4)/  0.83627 35565 080          D   -7/
      DATA ARJ0A1(5)/ -0.52585 70367 31           D   -8/
      DATA ARJ0A1(6)/  0.42606 36325 1            D   -9/
      DATA ARJ0A1(7)/ -0.42117 61024              D  -10/
      DATA ARJ0A1(8)/  0.48894 6426               D  -11/
      DATA ARJ0A1(9)/ -0.64834 929                D  -12/
      DATA ARJ0A1(10)/ 0.96172 34                 D  -13/
      DATA ARJ0A1(11)/-0.15703 67                 D  -13/
      DATA ARJ0A1(12)/ 0.27871 2                  D  -14/
      DATA ARJ0A1(13)/-0.53222                    D  -15/
      DATA ARJ0A1(14)/ 0.10844                    D  -15/
      DATA ARJ0A1(15)/-0.2342                     D  -16/
      DATA ARJ0A1(16)/ 0.533                      D  -17/
      DATA ARJ0A1(17)/-0.127                      D  -17/
      DATA ARJ0A1(18)/ 0.32                       D  -18/
      DATA ARJ0A1(19)/-0.8                        D  -19/
      DATA ARJ0A1(20)/ 0.2                        D  -19/
      DATA ARJ0A1(21)/-0.1                        D  -19/
      DATA ARJ0A2(0)/  1.99616 09630 13416 75339  D    0/
      DATA ARJ0A2(1)/ -0.19037 98192 46668 161    D   -2/
      DATA ARJ0A2(2)/  0.15397 10927 04422 6      D   -4/
      DATA ARJ0A2(3)/ -0.31145 08832 8103         D   -6/
      DATA ARJ0A2(4)/  0.11108 50971 321          D   -7/
      DATA ARJ0A2(5)/ -0.58666 78712 3            D   -9/
      DATA ARJ0A2(6)/  0.41399 26949              D  -10/
      DATA ARJ0A2(7)/ -0.36539 8763               D  -11/
      DATA ARJ0A2(8)/  0.38557 568                D  -12/
      DATA ARJ0A2(9)/ -0.47098 00                 D  -13/
      DATA ARJ0A2(10)/ 0.65022 0                  D  -14/
      DATA ARJ0A2(11)/-0.99624                    D  -15/
      DATA ARJ0A2(12)/ 0.16700                    D  -15/
      DATA ARJ0A2(13)/-0.3028                     D  -16/
      DATA ARJ0A2(14)/ 0.589                      D  -17/
      DATA ARJ0A2(15)/-0.122                      D  -17/
      DATA ARJ0A2(16)/ 0.27                       D  -18/
      DATA ARJ0A2(17)/-0.6                        D  -19/
      DATA ARJ0A2(18)/ 0.1                        D  -19/
C
C   Machine-dependent constants (suitable for IEEE machines)
C
      DATA NTERM1,NTERM2,NTERM3/22,18,16/
      DATA XLOW,XHIGH/3.650024D-8,9.0072D15/
C
C   Start computation
C
      X = XVALUE
      IND = 1
      IF ( X .LT. ZERO ) THEN
         X = -X
         IND = -1
      ENDIF
C
C   Error test 
C
      IF ( X .GT. XHIGH ) THEN
         CALL ERRPRN(FNNAME,ERRMSG)
         J0INT = ONE
         IF ( IND .EQ. -1 ) J0INT = -J0INT
         RETURN
      ENDIF
C
C   Code for 0 <= |x| <= 16
C
      IF ( X .LE. SIXTEN ) THEN
         IF ( X .LT. XLOW ) THEN
            J0INT = X
         ELSE
            T = X * X / ONE28 - ONE
            J0INT = X * CHEVAL(NTERM1,ARJ01,T)
         ENDIF
      ELSE
C
C   Code for |x| > 16
C
         T = FIVE12 / ( X * X ) - ONE
         PIB41 = PIB411 / PIB412
         XMPI4 = ( X - PIB41 ) - PIB42
         TEMP = COS(XMPI4) * CHEVAL(NTERM2,ARJ0A1,T) / X
         TEMP = TEMP - SIN(XMPI4) * CHEVAL(NTERM3,ARJ0A2,T)
         J0INT = ONE - RT2BPI * TEMP / SQRT(X)
      ENDIF
      IF ( IND .EQ. -1 ) J0INT = -J0INT
      RETURN
      END

      DOUBLE PRECISION FUNCTION K0INT(XVALUE)

c*********************************************************************72
c
cc K0INT calculates the integral of the modified Bessel function K0(X).
c
C   DESCRIPTION:
C
C      This function calculates the integral of the modified Bessel function
C      defined by
C
C         K0INT(x) = {integral 0 to x} K0(t) dt
C
C      The code uses Chebyshev expansions, whose coefficients are
C      given to 20 decimal places.
C
C      This subroutine is set up to work on IEEE machines.
C      For other machines, you should retrieve the code
C      from the general MISCFUN archive.
C
C
C   ERROR RETURNS:
C
C      If XVALUE < 0.0, the function is undefined. An error message is
C      printed and the function returns the value 0.0.
C
C
C   MACHINE-DEPENDENT CONSTANTS:
C
C      NTERM1 - The no. of terms to be used in the array AK0IN1. The 
C                recommended value is such that
C                   ABS(AK0IN1(NTERM1)) < EPS/100, 
C
C      NTERM2 - The no. of terms to be used in the array AK0IN2. The 
C                recommended value is such that
C                   ABS(AK0IN2(NTERM2)) < EPS/100,
C
C      NTERM3 - The no. of terms to be used in the array AK0INA. The 
C                recommended value is such that
C                   ABS(AK0INA(NTERM3)) < EPS/100, 
C
C      XLOW - The value below which K0INT = x * ( const - ln(x) ) to
C             machine precision. The recommended value is
C                   sqrt (18*EPSNEG).
C
C      XHIGH - The value above which K0INT = pi/2 to machine precision.
C              The recommended value is
C                   - log (2*EPSNEG)
C
C      For values of EPS and EPSNEG refer to the file MACHCON.TXT.
C
C      The machine-arithmetic constants are given in DATA
C      statements.
C
C
C   INTRINSIC FUNCTIONS USED:
C
C      EXP , LOG , SQRT
C
C
C   OTHER MISCFUN SUBROUTINES USED:
C
C          CHEVAL , ERRPRN
C
C
C   AUTHOR:
C         Dr. Allan J. MacLeod,
C         Dept. of Mathematics and Statistics,
C         University of Paisley,
C         High St.,
C         Paisley,
C         SCOTLAND
C
C         (e-mail: macl_ms0@paisley.ac.uk )
C
C
C   LATEST REVISION:
C                   11 JANUARY, 1996
C
      INTEGER NTERM1,NTERM2,NTERM3
      DOUBLE PRECISION AK0IN1(0:15),AK0IN2(0:15),AK0INA(0:27),
     1     CHEVAL,CONST1,CONST2,EIGHTN,FVAL,HALF,
     2     ONEHUN,PIBY2,RT2BPI,SIX,T,TEMP,TWELVE,X,
     3     XHIGH,XLOW,XVALUE,ZERO
      CHARACTER FNNAME*8,ERRMSG*14
      DATA FNNAME/'K0INT '/
      DATA ERRMSG/'ARGUMENT < 0.0'/
      DATA ZERO,HALF,SIX/ 0.0 D 0 , 0.5 D 0 , 6.0 D 0 /
      DATA TWELVE,EIGHTN,ONEHUN/ 12.0 D 0 , 18.0 D 0 , 100.0 D 0 /
      DATA CONST1/1.11593 15156 58412 44881 D 0/
      DATA CONST2/-0.11593 15156 58412 44881 D 0/
      DATA PIBY2/1.57079 63267 94896 61923 D 0/
      DATA RT2BPI/0.79788 45608 02865 35588 D 0/
      DATA AK0IN1/16.79702 71446 47109 59477  D    0,
     1             9.79134 68767 68894 07070  D    0,
     2             2.80501 31604 43379 39300  D    0,
     3             0.45615 62053 18885 02068  D    0,
     4             0.47162 24457 07476 0784   D   -1,
     5             0.33526 51482 69698 289    D   -2,
     6             0.17335 18119 38747 27     D   -3,
     7             0.67995 18893 64702        D   -5,
     8             0.20900 26835 9924         D   -6,
     9             0.51660 38469 76           D   -8,
     X             0.10485 70833 1            D   -9,
     1             0.17782 9320               D  -11,
     2             0.25568 44                 D  -13,
     3             0.31557                    D  -15,
     4             0.338                      D  -17,
     5             0.3                        D  -19/
      DATA AK0IN2/10.76266 55822 78091 74077  D    0,
     1             5.62333 47984 99975 11550  D    0,
     2             1.43543 66487 92908 67158  D    0,
     3             0.21250 41014 37438 96043  D    0,
     4             0.20365 37393 10000 9554   D   -1,
     5             0.13602 35840 95623 632    D   -2,
     6             0.66753 88699 20909 3      D   -4,
     7             0.25043 00357 07337        D   -5,
     8             0.74064 23741 728          D   -7,
     9             0.17697 47043 14           D   -8,
     X             0.34857 75254              D  -10,
     1             0.57544 785                D  -12,
     2             0.80748 1                  D  -14,
     3             0.9747                     D  -16,
     4             0.102                      D  -17,
     5             0.1                        D  -19/
      DATA AK0INA(0)/  1.91172 06544 50604 53895  D    0/
      DATA AK0INA(1)/ -0.41830 64565 76958 1085   D   -1/
      DATA AK0INA(2)/  0.21335 25080 68147 486    D   -2/
      DATA AK0INA(3)/ -0.15859 49728 45041 81     D   -3/
      DATA AK0INA(4)/  0.14976 24699 85835 1      D   -4/
      DATA AK0INA(5)/ -0.16795 59553 22241        D   -5/
      DATA AK0INA(6)/  0.21495 47247 8804         D   -6/
      DATA AK0INA(7)/ -0.30583 56654 790          D   -7/
      DATA AK0INA(8)/  0.47494 64133 43           D   -8/
      DATA AK0INA(9)/ -0.79424 66043 2            D   -9/
      DATA AK0INA(10)/ 0.14156 55532 5            D   -9/
      DATA AK0INA(11)/-0.26678 25359              D  -10/
      DATA AK0INA(12)/ 0.52814 9717               D  -11/
      DATA AK0INA(13)/-0.10926 3199               D  -11/
      DATA AK0INA(14)/ 0.23518 838                D  -12/
      DATA AK0INA(15)/-0.52479 91                 D  -13/
      DATA AK0INA(16)/ 0.12101 91                 D  -13/
      DATA AK0INA(17)/-0.28763 2                  D  -14/
      DATA AK0INA(18)/ 0.70297                    D  -15/
      DATA AK0INA(19)/-0.17631                    D  -15/
      DATA AK0INA(20)/ 0.4530                     D  -16/
      DATA AK0INA(21)/-0.1190                     D  -16/
      DATA AK0INA(22)/ 0.319                      D  -17/
      DATA AK0INA(23)/-0.87                       D  -18/
      DATA AK0INA(24)/ 0.24                       D  -18/
      DATA AK0INA(25)/-0.7                        D  -19/
      DATA AK0INA(26)/ 0.2                        D  -19/
      DATA AK0INA(27)/-0.1                        D  -19/
C
C   Machine-dependent values (suitable for IEEE machines)
C
      DATA NTERM1,NTERM2,NTERM3/14,14,23/
      DATA XLOW,XHIGH/4.47034836D-8,36.0436534D0/
C
C   Start computation
C
      X = XVALUE
C
C   Error test
C
      IF ( X .LT. ZERO ) THEN
         CALL ERRPRN(FNNAME,ERRMSG)
         K0INT = ZERO
         RETURN
      ENDIF
C
C   Code for 0 <= XVALUE <= 6
C
      IF ( X .LE. SIX ) THEN
         IF ( X .LT. XLOW ) THEN
            FVAL = X
            IF ( X .GT. ZERO ) THEN
               FVAL = FVAL * ( CONST1 - LOG(X) )
            ENDIF
            K0INT = FVAL
         ELSE
            T = ( ( X * X ) / EIGHTN - HALF ) - HALF
            FVAL = ( CONST2 + LOG(X) ) * CHEVAL(NTERM2,AK0IN2,T)
            K0INT = X * ( CHEVAL(NTERM1,AK0IN1,T) - FVAL )
         ENDIF
C
C   Code for x > 6
C
      ELSE
         FVAL = PIBY2
         IF ( X .LT. XHIGH ) THEN
            T = ( TWELVE / X - HALF ) - HALF
            TEMP = EXP(-X) * CHEVAL(NTERM3,AK0INA,T)
            FVAL = FVAL - TEMP / ( SQRT(X) * RT2BPI)
         ENDIF
         K0INT = FVAL
      ENDIF
      RETURN
      END

      DOUBLE PRECISION FUNCTION LOBACH(XVALUE)

c*********************************************************************72
c
cc LOBACH calculates the Lobachevsky function.
c
C   DESCRIPTION:
C
C      This function calculates the Lobachewsky function L(x), defined as
C
C         LOBACH(x) = {integral 0 to x} ( -ln ( | cos t | ) dt
C
C      The code uses Chebyshev expansions whose coefficients are given
C      to 20 decimal places.
C
C      This subroutine is set up to work on IEEE machines.
C      For other machines, you should retrieve the code
C      from the general MISCFUN archive.
C
C
C   ERROR RETURNS:
C
C      If |x| too large, it is impossible to accurately reduce the
C      argument to the range [0,pi]. An error message is printed
C      and the program returns the value 0.0
C
C
C   MACHINE-DEPENDENT CONSTANTS:
C
C      NTERM1 - INTEGER - The no. of terms to be used of the array ARLOB1.
C                          The recommended value is such that
C                          ABS(ARLOB1(NTERM1)) < EPS/100
C
C      NTERM2 - INTEGER - The no. of terms to be used of the array ARLOB2.
C                          The recommended value is such that
C                          ABS(ARLOB2(NTERM2)) < EPS/100
C
C      XLOW1 - DOUBLE PRECISION - The value below which L(x) = 0.0 to machine-precision.
C                     The recommended value is
C                              cube-root ( 6*XMIN )
C
C      XLOW2 - DOUBLE PRECISION - The value below which L(x) = x**3/6 to
C                     machine-precision. The recommended value is
C                              sqrt ( 10*EPS ) 
C
C      XLOW3 - DOUBLE PRECISION - The value below which
C                         L(pi/2) - L(pi/2-x) = x ( 1 - log(x) )
C                     to machine-precision. The recommended value is
C                               sqrt ( 18*EPS )
C
C      XHIGH - DOUBLE PRECISION - The value of |x| above which it is impossible
C                     to accurately reduce the argument. The
C                     recommended value is   1 / EPS.
C
C      For values of EPS, and XMIN, refer to the file MACHCON.TXT
C
C      The machine-arithmetic constants are given in DATA
C      statements.
C
C
C   INTRINSIC FUNCTIONS USED:
C
C      INT , LOG , SQRT
C
C
C   OTHER MISCFUN SUBROUTINES USED:
C
C          CHEVAL , ERRPRN
C
C
C   AUTHOR:
C
C      Dr. Allan J. MacLeod,
C      Dept. of Mathematics and Statistics,
C      University of Paisley,
C      High St.,
C      Paisley,
C      SCOTLAND
C
C      ( e-mail: macl_ms0@paisley.ac.uk )
C
C
C   LATEST UPDATE:
C                 11 JANUARY, 1996
C
      INTEGER INDPI2,INDSGN,NPI,NTERM1,NTERM2
      DOUBLE PRECISION ARLOB1(0:15),ARLOB2(0:10),
     1     CHEVAL,FVAL,FVAL1,HALF,LBPB21,LBPB22,LOBPIA,LOBPIB,
     2     LOBPI1,LOBPI2,ONE,ONEHUN,PI,PIBY2,PIBY21,PIBY22,PIBY4,PI1,
     3     PI11,PI12,PI2,SIX,T,TCON,TEN,TWO,X,XCUB,XHIGH,XLOW1,
     4     XLOW2,XLOW3,XR,XVALUE,ZERO
      CHARACTER FNNAME*6,ERRMSG*26
      DATA FNNAME/'LOBACH'/
      DATA ERRMSG/'ARGUMENT TOO LARGE IN SIZE'/
      DATA ZERO,HALF/ 0.0 D 0 , 0.5 D 0 /
      DATA ONE,TWO,SIX/ 1.0 D 0 , 2.0 D 0 , 6.0 D 0 /
      DATA TEN,ONEHUN/ 10.0 D 0 , 100.0 D 0 /
      DATA LOBPIA,LOBPIB/ 1115.0 D 0 , 512.0 D 0 /
      DATA LOBPI2/-1.48284 69639 78694 99311 D -4/
      DATA LBPB22/-7.41423 48198 93474 96556 D -5/
      DATA PI11,PI12/ 201.0 D 0 , 64.0 D 0 /
      DATA PI2/9.67653 58979 32384 62643 D -4/
      DATA PIBY22/4.83826 79489 66192 31322 D -4/
      DATA TCON/3.24227 78765 54808 68620 D 0/
      DATA ARLOB1/0.34464 88495 34813 00507  D    0,
     1            0.58419 83571 90277 669    D   -2,
     2            0.19175 02969 46003 30     D   -3,
     3            0.78725 16064 56769        D   -5,
     4            0.36507 47741 5804         D   -6,
     5            0.18302 87272 680          D   -7,
     6            0.96890 33300 5            D   -9,
     7            0.53390 55444              D  -10,
     8            0.30340 8025               D  -11,
     9            0.17667 875                D  -12,
     X            0.10493 93                 D  -13,
     1            0.63359                    D  -15,
     2            0.3878                     D  -16,
     3            0.240                      D  -17,
     4            0.15                       D  -18,
     5            0.1                        D  -19/
      DATA ARLOB2/2.03459 41803 61328 51087  D    0,
     1            0.17351 85882 02740 7681   D   -1,
     2            0.55162 80426 09052 1      D   -4,
     3            0.39781 64627 6598         D   -6,
     4            0.36901 80289 18           D   -8,
     5            0.38804 09214              D  -10,
     6            0.44069 698                D  -12,
     7            0.52767 4                  D  -14,
     8            0.6568                     D  -16,
     9            0.84                       D  -18,
     X            0.1                        D  -19/
C
C   Machine-dependent constants (suitable for IEEE machines)
C
      DATA NTERM1,NTERM2/13,9/
      DATA XLOW1,XLOW2/5.11091385D-103,4.71216091D-8/
      DATA XLOW3,XHIGH/6.32202727D-8,4.5035996D15/
C
C   Start computation
C
      X = ABS ( XVALUE )
      INDSGN = 1
      IF ( XVALUE .LT. ZERO ) THEN
         INDSGN = -1
      ENDIF
C
C   Error test
C
      IF ( X .GT. XHIGH ) THEN
         CALL ERRPRN(FNNAME,ERRMSG)
         LOBACH = ZERO
         RETURN
      ENDIF
C
C   Reduce argument to [0,pi]
C
      PI1 = PI11/PI12
      PI = PI1 + PI2
      PIBY2 = PI/TWO
      PIBY21 = PI1/TWO
      PIBY4 = PIBY2/TWO
      NPI = INT ( X / PI )
      XR = ( X - NPI * PI1 ) - NPI * PI2
C
C   Reduce argument to [0,pi/2]
C
      INDPI2 = 0
      IF ( XR .GT. PIBY2 ) THEN
         INDPI2 = 1
         XR = ( PI1 - XR ) + PI2
      ENDIF
C
C   Code for argument in [0,pi/4]
C
      IF ( XR .LE. PIBY4 ) THEN 
         IF ( XR .LT. XLOW1 ) THEN
            FVAL = ZERO
         ELSE
            XCUB = XR * XR * XR
            IF ( XR .LT. XLOW2 ) THEN
               FVAL = XCUB / SIX
            ELSE
               T = ( TCON * XR * XR - HALF ) - HALF
               FVAL = XCUB * CHEVAL(NTERM1,ARLOB1,T)
            ENDIF
         ENDIF
      ELSE
C
C   Code for argument in [pi/4,pi/2]
C
         XR = ( PIBY21 - XR ) + PIBY22
         IF ( XR .EQ. ZERO ) THEN
            FVAL1 = ZERO
         ELSE
            IF ( XR .LT. XLOW3 ) THEN
               FVAL1 = XR * ( ONE - LOG( XR ) )
            ELSE
               T = ( TCON * XR * XR - HALF ) - HALF
               FVAL1 = XR * ( CHEVAL(NTERM2,ARLOB2,T) - LOG( XR ) )
            ENDIF
         ENDIF
         LBPB21 = LOBPIA / ( LOBPIB + LOBPIB )
         FVAL = ( LBPB21 - FVAL1 ) + LBPB22
      ENDIF 
      LOBPI1 = LOBPIA / LOBPIB
C
C   Compute value for argument in [pi/2,pi]
C
      IF ( INDPI2 .EQ. 1 ) THEN
         FVAL = ( LOBPI1 - FVAL ) + LOBPI2
      ENDIF
      LOBACH = FVAL
C
C   Scale up for arguments > pi
C
      IF ( NPI .GT. 0 ) THEN
         LOBACH = ( FVAL + NPI * LOBPI2 ) + NPI * LOBPI1
      ENDIF
      IF ( INDSGN .EQ. -1 ) THEN
         LOBACH = - LOBACH
      ENDIF
      RETURN
      END

      DOUBLE PRECISION FUNCTION STROM(XVALUE)

c*********************************************************************72
c
cc STROM calculates Stromgen's integral.
c
C  DESCRIPTION:
C
C      This program calculates Stromgren's integral, defined as
C 
C          STROM(X) = integral 0 to X { t**7 exp(2t)/[exp(t)-1]**3 } dt
C
C      The code uses a Chebyshev series, the coefficients of which are
C      given to an accuracy of 20 decimal places.
C
C      This subroutine is set up to work on IEEE machines.
C      For other machines, you should retrieve the code
C      from the general MISCFUN archive.
C
C
C  ERROR RETURNS:
C
C    If XVALUE < 0.0, an error message is printed, and the program 
C    returns the value 0.0. 
C
C
C  MACHINE-DEPENDENT CONSTANTS:
C
C    NTERMS - INTEGER - The number of terms of the array ASTROM to be used.
C                       The recommended value is such that
C                             ASTROM(NTERMS) < EPS/100
C
C    XLOW0 - DOUBLE PRECISION - The value below which STROM = 0.0 to machine
C                    precision. The recommended value is
C                          5th root of (130*XMIN)
C
C    XLOW1 - DOUBLE PRECISION - The value below which STROM = 3*(X**5)/(4*(pi**4))
C                   to machine precision. The recommended value is
C                             2*EPSNEG
C
C    EPSLN - DOUBLE PRECISION - The value of ln(EPS). Used to determine the no.
C                   of exponential terms for large X.
C
C    EPNGLN - DOUBLE PRECISION - The value of ln(EPSNEG). Used to prevent
C                    overflow for large X.
C
C    XHIGH - DOUBLE PRECISION - The value above which
C                           STROM = 196.52 - 15*(x**7)*exp(-x)/(4pi**4)
C                   to machine precision. The recommended value is
C                             7 / EPS
C
C     For values of EPS, EPSNEG, and XMIN refer to the file MACHCON.TXT
C
C     The machine-arithmetic constants are given in DATA
C     statements.
C
C
C  INTRINSIC FUNCTIONS USED:
C
C    EXP, INT, LOG
C
C
C   OTHER MISCFUN SUBROUTINES USED:
C
C          CHEVAL , ERRPRN
C
C
C  AUTHOR:
C
C     DR. ALLAN J. MACLEOD,
C     DEPT. OF MATHEMATICS AND STATISTICS,
C     UNIVERSITY OF PAISLEY,
C     HIGH ST.,
C     PAISLEY,
C     SCOTLAND.
C     PA1 2BE.
C
C     (e-mail: macl_ms0@paisley.ac.uk )
C
C
C  LATEST REVISION: 
C                  11 JANUARY, 1996
C
C
      INTEGER K1,K2,NTERMS,NUMEXP
      DOUBLE PRECISION ASTROM(0:26),CHEVAL,EPNGLN,EPSLN,FOUR,
     1     F15BP4,HALF,ONE,ONEHUN,ONE30,ONE5LN,PI4B3,RK,
     2     SEVEN,SUMEXP,SUM2,T,TWO,VALINF,X,XHIGH,
     3     XK,XK1,XLOW0,XLOW1,XVALUE,ZERO
      CHARACTER FNNAME*6,ERRMSG*14
      DATA FNNAME/'STROM '/
      DATA ERRMSG/'ARGUMENT < 0.0'/
      DATA ZERO,HALF,ONE/ 0.0 D 0 , 0.5 D 0 , 1.0 D 0/
      DATA TWO,FOUR,SEVEN/ 2.0 D 0 , 4.0 D 0 , 7.0 D 0 /
      DATA ONEHUN,ONE30,ONE5LN/ 100.0 D 0 , 130.0 D 0 , 0.4055 D 0 /
      DATA F15BP4/0.38497 43345 50662 56959 D -1 /
      DATA PI4B3/1.29878 78804 53365 82982 D 2 /
      DATA VALINF/196.51956 92086 89882 61257 D 0/
      DATA ASTROM(0)/  0.56556 12087 25391 55290  D    0/
      DATA ASTROM(1)/  0.45557 31969 10178 5525   D   -1/
      DATA ASTROM(2)/ -0.40395 35875 93686 9170   D   -1/
      DATA ASTROM(3)/ -0.13339 05720 21486 815    D   -2/
      DATA ASTROM(4)/  0.18586 25062 50538 030    D   -2/
      DATA ASTROM(5)/ -0.46855 55868 05365 9      D   -4/
      DATA ASTROM(6)/ -0.63434 75643 42294 9      D   -4/
      DATA ASTROM(7)/  0.57254 87081 43200        D   -5/
      DATA ASTROM(8)/  0.15935 28122 16822        D   -5/
      DATA ASTROM(9)/ -0.28884 32843 1036         D   -6/
      DATA ASTROM(10)/-0.24466 33604 801          D   -7/
      DATA ASTROM(11)/ 0.10072 50382 374          D   -7/
      DATA ASTROM(12)/-0.12482 98610 4            D   -9/
      DATA ASTROM(13)/-0.26300 62528 3            D   -9/
      DATA ASTROM(14)/ 0.24904 07578              D  -10/
      DATA ASTROM(15)/ 0.48545 4902               D  -11/
      DATA ASTROM(16)/-0.10537 8913               D  -11/
      DATA ASTROM(17)/-0.36044 17                 D  -13/
      DATA ASTROM(18)/ 0.29920 78                 D  -13/
      DATA ASTROM(19)/-0.16397 1                  D  -14/
      DATA ASTROM(20)/-0.61061                    D  -15/
      DATA ASTROM(21)/ 0.9335                     D  -16/
      DATA ASTROM(22)/ 0.709                      D  -17/
      DATA ASTROM(23)/-0.291                      D  -17/
      DATA ASTROM(24)/ 0.8                        D  -19/
      DATA ASTROM(25)/ 0.6                        D  -19/
      DATA ASTROM(26)/-0.1                        D  -19/
C
C  Machine-dependent constants
C
      DATA NTERMS/23/
      DATA XLOW0,XLOW1/7.80293D-62,2.22045D-16/
      DATA EPSLN,EPNGLN/-36.0436534D0,-36.7368006D0/
      DATA XHIGH/3.1525197D16/
C
C  Start execution
C
      X = XVALUE
C
C  Error test
C
      IF ( X .LT. ZERO ) THEN
         CALL ERRPRN(FNNAME,ERRMSG)
         STROM = ZERO
         RETURN
      ENDIF
C
C   Code for x < =  4.0
C
      IF ( X .LE. FOUR ) THEN
         IF ( X .LT. XLOW0 ) THEN
            STROM = ZERO
         ELSE
            IF ( X .LT. XLOW1 ) THEN
               STROM = (X**5) / PI4B3
            ELSE
               T = ( ( X / TWO ) - HALF ) - HALF
               STROM = (X**5) * CHEVAL(NTERMS,ASTROM,T) * F15BP4
            ENDIF
         ENDIF
      ELSE
C 
C  Code for x > 4.0
C
         IF ( X .GT. XHIGH ) THEN
            SUMEXP = ONE
         ELSE
            NUMEXP = INT( EPSLN / (ONE5LN - X ) ) + 1
            IF ( NUMEXP .GT. 1 ) THEN
               T = EXP( -X )
            ELSE
               T = ONE
            ENDIF
            RK = ZERO
            DO 100 K1 = 1 , NUMEXP
               RK = RK + ONE
  100       CONTINUE
            SUMEXP = ZERO
            DO 300 K1 = 1 , NUMEXP
               SUM2 = ONE
               XK = ONE / ( RK * X )
               XK1 = ONE
               DO 200 K2 = 1 , 7
                  SUM2 = SUM2 * XK1 * XK + ONE
                  XK1 = XK1 + ONE
  200          CONTINUE
               SUM2 = SUM2 * ( RK + ONE ) / TWO
               SUMEXP = SUMEXP * T + SUM2
               RK = RK - ONE
  300       CONTINUE
         ENDIF
         T = SEVEN * LOG(X) - X + LOG(SUMEXP)
         IF ( T .LT. EPNGLN ) THEN
            STROM = VALINF
         ELSE
            STROM = VALINF - EXP(T) * F15BP4
         ENDIF
      ENDIF
      RETURN
      END

      DOUBLE PRECISION FUNCTION STRVH0(XVALUE)

c*********************************************************************72
c
cc STRVH0 calculates the Struve function of order 0.
c
C   DESCRIPTION:
C
C      This function calculates the value of the Struve function
C      of order 0, denoted H0(x), for the argument XVALUE, defined
C
C         STRVHO(x) = (2/pi) integral{0 to pi/2} sin(x cos(t)) dt
C
C      H0 also satisfies the second-order equation
C
C                 x*D(Df) + Df + x*f = 2x/pi
C
C      The code uses Chebyshev expansions whose coefficients are
C      given to 20D.
C
C      This subroutine is set up to work on IEEE machines.
C      For other machines, you should retrieve the code
C      from the general MISCFUN archive.
C
C
C   ERROR RETURNS:
C
C      As the asymptotic expansion of H0 involves the Bessel function
C      of the second kind Y0, there is a problem for large x, since
C      we cannot accurately calculate the value of Y0. An error message 
C      is printed and STRVH0 returns the value 0.0.
C
C
C   MACHINE-DEPENDENT CONSTANTS:
C
C      NTERM1 - The no. of terms to be used in the array ARRH0. The
C               recommended value is such that
C                      ABS(ARRH0(NTERM1)) < EPS/100.
C
C      NTERM2 - The no. of terms to be used in the array ARRH0A. The
C               recommended value is such that
C                      ABS(ARRH0A(NTERM2)) < EPS/100.
C
C      NTERM3 - The no. of terms to be used in the array AY0ASP. The
C               recommended value is such that
C                      ABS(AY0ASP(NTERM3)) < EPS/100.
C
C      NTERM4 - The no. of terms to be used in the array AY0ASQ. The
C               recommended value is such that
C                      ABS(AY0ASQ(NTERM4)) < EPS/100.
C
C      XLOW - The value for which H0(x) = 2*x/pi to machine precision, if
C             abs(x) < XLOW. The recommended value is
C                      XLOW = 3 * SQRT(EPSNEG)
C
C      XHIGH - The value above which we are unable to calculate Y0 with
C              any reasonable accuracy. An error message is printed and
C              STRVH0 returns the value 0.0. The recommended value is
C                      XHIGH = 1/EPS.
C
C      For values of EPS and EPSNEG refer to the file MACHCON.TXT.
C
C      The machine-arithmetic constants are given in DATA
C      statements.
C
C
C   INTRINSIC FUNCTIONS USED:
C
C      ABS, COS, SIN, SQRT.
C
C
C   OTHER MISCFUN SUBROUTINES USED:
C
C          CHEVAL , ERRPRN
C
C
C   AUTHOR:
C          ALLAN J. MACLEOD
C          DEPT. OF MATHEMATICS AND STATISTICS
C          UNIVERSITY OF PAISLEY
C          HIGH ST.
C          PAISLEY
C          SCOTLAND
C          PA1 2BE
C
C          (e-mail: macl_ms0@paisley.ac.uk )
C
C
C   LATEST REVISION:
C                   11 JANUARY, 1996
C
C
      INTEGER INDSGN,NTERM1,NTERM2,NTERM3,NTERM4
      DOUBLE PRECISION ARRH0(0:19),ARRH0A(0:20),AY0ASP(0:12),
     1     AY0ASQ(0:13),CHEVAL,EIGHT,ELEVEN,HALF,H0AS,
     2     ONEHUN,ONE,PIBY4,RT2BPI,SIXTP5,T,THR2P5,TWENTY,
     3     TWOBPI,TWO62,X,XHIGH,XLOW,XMP4,XSQ,XVALUE,
     4     Y0P,Y0Q,Y0VAL,ZERO
      CHARACTER FNNAME*6,ERRMSG*26
      DATA FNNAME/'STRVH0'/
      DATA ERRMSG/'ARGUMENT TOO LARGE IN SIZE'/
      DATA ZERO,HALF,ONE/0.0 D 0 , 0.5 D 0 , 1.0 D 0/ 
      DATA EIGHT,ELEVEN/8.0 D 0 , 11.0 D 0/
      DATA TWENTY,ONEHUN/20.0 D 0 , 100.0 D 0/
      DATA SIXTP5,TWO62,THR2P5/60.5 D 0 , 262.0 D 0 , 302.5 D 0/
      DATA PIBY4/0.78539 81633 97448 30962 D 0/
      DATA RT2BPI/0.79788 45608 02865 35588 D 0/
      DATA TWOBPI/0.63661 97723 67581 34308 D 0/
      DATA ARRH0(0)/  0.28696 48739 90132 25740  D    0/
      DATA ARRH0(1)/ -0.25405 33268 16183 52305  D    0/
      DATA ARRH0(2)/  0.20774 02673 93238 94439  D    0/
      DATA ARRH0(3)/ -0.20364 02956 03865 85140  D    0/
      DATA ARRH0(4)/  0.12888 46908 68661 86016  D    0/
      DATA ARRH0(5)/ -0.48256 32815 62226 1202   D   -1/
      DATA ARRH0(6)/  0.11686 29347 56900 1242   D   -1/
      DATA ARRH0(7)/ -0.19811 81356 42418 416    D   -2/
      DATA ARRH0(8)/  0.24899 13851 24212 86     D   -3/
      DATA ARRH0(9)/ -0.24188 27913 78595 0      D   -4/
      DATA ARRH0(10)/ 0.18743 75479 93431        D   -5/
      DATA ARRH0(11)/-0.11873 34607 4362         D   -6/
      DATA ARRH0(12)/ 0.62698 49433 46           D   -8/
      DATA ARRH0(13)/-0.28045 54679 3            D   -9/
      DATA ARRH0(14)/ 0.10769 41205              D  -10/
      DATA ARRH0(15)/-0.35904 793                D  -12/
      DATA ARRH0(16)/ 0.10494 47                 D  -13/
      DATA ARRH0(17)/-0.27119                    D  -15/
      DATA ARRH0(18)/ 0.624                      D  -17/
      DATA ARRH0(19)/-0.13                       D  -18/
      DATA ARRH0A(0)/  1.99291 88575 19923 05515  D    0/
      DATA ARRH0A(1)/ -0.38423 26687 01456 887    D   -2/
      DATA ARRH0A(2)/ -0.32871 99371 23530 50     D   -3/
      DATA ARRH0A(3)/ -0.29411 81203 70340 9      D   -4/
      DATA ARRH0A(4)/ -0.26731 53519 87066        D   -5/
      DATA ARRH0A(5)/ -0.24681 03107 5013         D   -6/
      DATA ARRH0A(6)/ -0.22950 14861 143          D   -7/
      DATA ARRH0A(7)/ -0.21568 22318 33           D   -8/
      DATA ARRH0A(8)/ -0.20303 50648 3            D   -9/
      DATA ARRH0A(9)/ -0.19345 75509              D  -10/
      DATA ARRH0A(10)/-0.18277 3144               D  -11/
      DATA ARRH0A(11)/-0.17768 424                D  -12/
      DATA ARRH0A(12)/-0.16432 96                 D  -13/
      DATA ARRH0A(13)/-0.17156 9                  D  -14/
      DATA ARRH0A(14)/-0.13368                    D  -15/
      DATA ARRH0A(15)/-0.2077                     D  -16/
      DATA ARRH0A(16)/ 0.2                        D  -19/
      DATA ARRH0A(17)/-0.55                       D  -18/
      DATA ARRH0A(18)/ 0.10                       D  -18/
      DATA ARRH0A(19)/-0.4                        D  -19/
      DATA ARRH0A(20)/ 0.1                        D  -19/
      DATA AY0ASP/1.99944 63940 23982 71568  D    0,
     1           -0.28650 77864 70319 58     D   -3,
     2           -0.10050 72797 43762 0      D   -4,
     3           -0.35835 94100 2463         D   -6,
     4           -0.12879 65120 531          D   -7,
     5           -0.46609 48663 6            D   -9,
     6           -0.16937 69454              D  -10,
     7           -0.61852 269                D  -12,
     8           -0.22618 41                 D  -13,
     9           -0.83268                    D  -15,
     X           -0.3042                     D  -16,
     1           -0.115                      D  -17,
     2           -0.4                        D  -19/
      DATA AY0ASQ/1.99542 68138 68286 04092  D    0,
     1           -0.23601 31928 67514 472    D   -2,
     2           -0.76015 38908 50296 6      D   -4,
     3           -0.25610 88714 56343        D   -5,
     4           -0.87502 92185 106          D   -7,
     5           -0.30430 42121 59           D   -8,
     6           -0.10621 42831 4            D   -9,
     7           -0.37737 1479               D  -11,
     8           -0.13213 687                D  -12,
     9           -0.48862 1                  D  -14,
     X           -0.15809                    D  -15,
     1           -0.762                      D  -17,
     2           -0.3                        D  -19,
     3           -0.3                        D  -19/
C
C   MACHINE-DEPENDENT CONSTANTS (Suitable for IEEE-arithmetic machines)
C
      DATA NTERM1,NTERM2,NTERM3,NTERM4/18,18,11,11/
      DATA XLOW,XHIGH/3.1610136D-8,4.50359963D15/
C
C   Start computation
C
      X = XVALUE
      INDSGN = 1
      IF ( X .LT. ZERO ) THEN
         X = -X
         INDSGN = -1
      ENDIF
C
C   Error test
C
      IF ( ABS(XVALUE) .GT. XHIGH ) THEN
         CALL ERRPRN(FNNAME,ERRMSG)
         STRVH0 = ZERO
         RETURN
      ENDIF
C
C   Code for abs(x) <= 11
C
      IF ( X .LE. ELEVEN ) THEN
         IF ( X .LT. XLOW ) THEN
            STRVH0 = TWOBPI * X
         ELSE
            T = ( ( X * X ) / SIXTP5 - HALF ) - HALF
            STRVH0 = TWOBPI * X * CHEVAL ( NTERM1 , ARRH0 , T )
         ENDIF
      ELSE      
C
C   Code for abs(x) > 11
C
         XSQ = X * X
         T = ( TWO62 - XSQ ) / ( TWENTY + XSQ )
         Y0P = CHEVAL ( NTERM3 , AY0ASP , T )
         Y0Q = CHEVAL ( NTERM4 , AY0ASQ , T ) / ( EIGHT * X )
         XMP4 = X - PIBY4
         Y0VAL = Y0P * SIN ( XMP4 ) - Y0Q * COS ( XMP4 )
         Y0VAL = Y0VAL * RT2BPI / SQRT ( X )
         T = ( THR2P5 - XSQ ) / ( SIXTP5 + XSQ )
         H0AS = TWOBPI * CHEVAL ( NTERM2 , ARRH0A , T ) / X
         STRVH0 = Y0VAL + H0AS
      ENDIF
      IF ( INDSGN .EQ. -1 ) STRVH0 = -STRVH0
      RETURN
      END

      DOUBLE PRECISION FUNCTION STRVH1(XVALUE)

c*********************************************************************72
c
cc STRVH1 calculates the Struve function of order 1.
c
C   DESCRIPTION:
C      This function calculates the value of the Struve function
C      of order 1, denoted H1(x), for the argument XVALUE, defined as
C
C                                                                  2
C        STRVH1(x) = (2x/pi) integral{0 to pi/2} sin( x cos(t))*sin t dt
C
C      H1 also satisfies the second-order differential equation
C
C                    2   2                   2            2
C                   x * D f  +  x * Df  +  (x - 1)f  =  2x / pi
C
C      The code uses Chebyshev expansions with the coefficients
C      given to 20D.
C
C      This subroutine is set up to work on IEEE machines.
C      For other machines, you should retrieve the code
C      from the general MISCFUN archive.
C
C
C   ERROR RETURNS:
C
C      As the asymptotic expansion of H1 involves the Bessel function
C      of the second kind Y1, there is a problem for large x, since
C      we cannot accurately calculate the value of Y1. An error message 
C      is printed and STRVH1 returns the value 0.0.
C
C
C   MACHINE-DEPENDENT CONSTANTS:
C
C      NTERM1 - The no. of terms to be used in the array ARRH1. The
C               recommended value is such that
C                      ABS(ARRH1(NTERM1)) < EPS/100.
C
C      NTERM2 - The no. of terms to be used in the array ARRH1A. The
C               recommended value is such that
C                      ABS(ARRH1A(NTERM2)) < EPS/100.
C
C      NTERM3 - The no. of terms to be used in the array AY1ASP. The
C               recommended value is such that
C                      ABS(AY1ASP(NTERM3)) < EPS/100.
C
C      NTERM4 - The no. of terms to be used in the array AY1ASQ. The
C               recommended value is such that
C                      ABS(AY1ASQ(NTERM4)) < EPS/100.
C
C      XLOW1 - The value of x, below which H1(x) set to zero, if
C              abs(x)<XLOW1. The recommended value is 
C                      XLOW1 = 1.5 * SQRT(XMIN)
C
C      XLOW2 - The value for which H1(x) = 2*x*x/pi to machine precision, if
C              abs(x) < XLOW2. The recommended value is
C                      XLOW2 = SQRT(15*EPSNEG)
C
C      XHIGH - The value above which we are unable to calculate Y1 with
C              any reasonable accuracy. An error message is printed and
C              STRVH1 returns the value 0.0. The recommended value is
C                      XHIGH = 1/EPS.
C
C      For values of EPS, EPSNEG and XMIN refer to the file MACHCON.TXT.
C
C      The machine-arithmetic constants are given in DATA
C      statements.
C
C
C   INTRINSIC FUNCTIONS USED:
C
C      ABS, COS, SIN, SQRT.
C
C
C   OTHER MISCFUN SUBROUTINES USED:
C
C          CHEVAL , ERRPRN
C
C
C   AUTHOR:
C          ALLAN J. MACLEOD
C          DEPT. OF MATHEMATICS AND STATISTICS
C          UNIVERSITY OF PAISLEY
C          HIGH ST.
C          PAISLEY
C          SCOTLAND
C          PA1 2BE
C
C          (e-mail: macl_ms0@paisley.ac.uk)
C
C
C   LATEST REVISION:
C                   12 JANUARY, 1996
C
C
      INTEGER NTERM1,NTERM2,NTERM3,NTERM4
      DOUBLE PRECISION ARRH1(0:17),ARRH1A(0:21),AY1ASP(0:14),
     1     AY1ASQ(0:15),CHEVAL,EIGHT,FIFTEN,FORTP5,HALF,
     2     H1AS,NINE,ONEHUN,ONE82,RT2BPI,T,THPBY4,
     3     TWENTY,TWOBPI,TW02P5,X,XHIGH,XLOW1,XLOW2,
     4     XM3P4,XSQ,XVALUE,Y1P,Y1Q,Y1VAL,ZERO
      CHARACTER FNNAME*6,ERRMSG*26
      DATA FNNAME/'STRVH1'/
      DATA ERRMSG/'ARGUMENT TOO LARGE IN SIZE'/
      DATA ZERO,HALF,EIGHT/0.0 D 0 , 0.5 D 0 , 8.0 D 0/
      DATA NINE,FIFTEN/ 9.0 D 0 , 15.0 D 0 / 
      DATA TWENTY,ONEHUN/ 20.0 D 0 , 100.0 D 0/
      DATA FORTP5,ONE82,TW02P5/40.5 D 0 , 182.0 D 0 , 202.5 D 0/
      DATA RT2BPI/0.79788 45608 02865 35588 D 0/
      DATA THPBY4/2.35619 44901 92344 92885 D 0/
      DATA TWOBPI/0.63661 97723 67581 34308 D 0/
      DATA ARRH1/0.17319 06108 36754 39319  D    0,
     1          -0.12606 91759 13526 72005  D    0,
     2           0.79085 76160 49535 7500   D   -1,
     3          -0.31964 93222 32187 0820   D   -1,
     4           0.80804 05814 04918 834    D   -2,
     5          -0.13600 08206 93074 148    D   -2,
     6           0.16227 14861 98894 71     D   -3,
     7          -0.14423 52451 48592 9      D   -4,
     8           0.99219 52573 4072         D   -6,
     9          -0.54416 28049 180          D   -7,
     X           0.24363 16625 63           D   -8,
     1          -0.90770 71338              D  -10,
     2           0.28592 6585               D  -11,
     3          -0.77169 75                 D  -13,
     4           0.18048 9                  D  -14,
     5          -0.3694                     D  -16,
     6           0.67                       D  -18,
     7          -0.1                        D  -19/
      DATA ARRH1A(0)/  2.01083 50495 14733 79407  D    0/
      DATA ARRH1A(1)/  0.59221 86100 36099 903    D   -2/
      DATA ARRH1A(2)/  0.55274 32269 84141 30     D   -3/
      DATA ARRH1A(3)/  0.52698 73856 31103 6      D   -4/
      DATA ARRH1A(4)/  0.50637 45221 40969        D   -5/
      DATA ARRH1A(5)/  0.49028 73642 0678         D   -6/
      DATA ARRH1A(6)/  0.47635 40023 525          D   -7/
      DATA ARRH1A(7)/  0.46525 86522 83           D   -8/
      DATA ARRH1A(8)/  0.45465 16608 1            D   -9/
      DATA ARRH1A(9)/  0.44724 62193              D  -10/
      DATA ARRH1A(10)/ 0.43730 8292               D  -11/
      DATA ARRH1A(11)/ 0.43568 368                D  -12/
      DATA ARRH1A(12)/ 0.41821 90                 D  -13/
      DATA ARRH1A(13)/ 0.44104 4                  D  -14/
      DATA ARRH1A(14)/ 0.36391                    D  -15/
      DATA ARRH1A(15)/ 0.5558                     D  -16/
      DATA ARRH1A(16)/-0.4                        D  -19/
      DATA ARRH1A(17)/ 0.163                      D  -17/
      DATA ARRH1A(18)/-0.34                       D  -18/
      DATA ARRH1A(19)/ 0.13                       D  -18/
      DATA ARRH1A(20)/-0.4                        D  -19/
      DATA ARRH1A(21)/ 0.1                        D  -19/
      DATA AY1ASP/2.00135 24004 58893 96402 D   0,
     1            0.71104 24159 64619 38    D  -3,
     2            0.36659 77028 23244 9     D  -4,
     3            0.19130 15686 57728       D  -5,
     4            0.10046 91138 9777        D  -6,
     5            0.53040 17425 38          D  -8,
     6            0.28100 88617 6           D  -9,
     7            0.14938 86051             D -10,
     8            0.79578 420               D -12,
     9            0.42523 63                D -13,
     X            0.22719 5                 D -14,
     1            0.12216                   D -15,
     2            0.650                     D -17,
     3            0.36                      D -18,
     4            0.2                       D -19/
      DATA AY1ASQ/5.99065 10947 78881 89116 D   0,
     1           -0.48959 32623 36579 635   D  -2,
     2           -0.23238 32130 70706 26    D  -3,
     3           -0.11447 34723 85767 9     D  -4,
     4           -0.57169 92618 9106        D  -6,
     5           -0.28955 16716 917         D  -7,
     6           -0.14751 33456 36          D  -8,
     7           -0.75965 37378             D -10,
     8           -0.39065 8184              D -11,
     9           -0.20464 654               D -12,
     X           -0.10426 36                D -13,
     1           -0.57702                   D -15,
     2           -0.2550                    D -16,
     3           -0.210                     D -17,
     4            0.2                       D -19,
     5           -0.2                       D -19/
C
C   MACHINE-DEPENDENT CONSTANTS (Suitable for IEEE-arithmetic machines)
C
      DATA NTERM1,NTERM2,NTERM3,NTERM4/15,17,12,13/
      DATA XLOW1,XLOW2/2.23750222D-154,4.08085106D-8/
      DATA XHIGH/4.50359963D15/
C
C   Start computation
C
      X = ABS ( XVALUE )
C
C   Error test
C
      IF ( X .GT. XHIGH ) THEN
         CALL ERRPRN(FNNAME,ERRMSG)
         STRVH1 = ZERO
         RETURN
      ENDIF
C
C   Code for abs(x) <= 9
C
      IF ( X .LE. NINE ) THEN
         IF ( X .LT. XLOW1 ) THEN
            STRVH1 = ZERO
         ELSE
            XSQ = X * X
            IF ( X .LT. XLOW2 ) THEN
               STRVH1 = TWOBPI * XSQ
            ELSE
               T = ( XSQ / FORTP5 - HALF ) - HALF
               STRVH1 = TWOBPI * XSQ * CHEVAL ( NTERM1 , ARRH1 , T )
            ENDIF
         ENDIF
      ELSE      
C
C   Code for abs(x) > 9
C
         XSQ = X * X
         T = ( ONE82 - XSQ ) / ( TWENTY + XSQ )
         Y1P = CHEVAL ( NTERM3 , AY1ASP , T )
         Y1Q = CHEVAL ( NTERM4 , AY1ASQ , T ) / ( EIGHT * X)
         XM3P4 = X - THPBY4
         Y1VAL = Y1P * SIN ( XM3P4 ) + Y1Q * COS ( XM3P4 )
         Y1VAL = Y1VAL * RT2BPI / SQRT ( X )
         T = ( TW02P5 - XSQ ) / ( FORTP5 + XSQ )
         H1AS = TWOBPI * CHEVAL ( NTERM2 , ARRH1A , T )
         STRVH1 = Y1VAL + H1AS
      ENDIF
      RETURN
      END

      DOUBLE PRECISION FUNCTION STRVL0(XVALUE)

c*********************************************************************72
c
cc STRVL0 calculates the modified Struve function of order 0.
c
C   DESCRIPTION:
C
C      This function calculates the modified Struve function of
C      order 0, denoted L0(x), defined as the solution of the
C      second-order equation
C
C                  x*D(Df) + Df - x*f  =  2x/pi
C
C      This subroutine is set up to work on IEEE machines.
C      For other machines, you should retrieve the code
C      from the general MISCFUN archive.
C
C
C   ERROR RETURNS:
C
C      If the value of |XVALUE| is too large, the result
C      would cause an floating-pt overflow. An error message
C      is printed and the function returns the value of
C      sign(XVALUE)*XMAX where XMAX is the largest possible
C      floating-pt argument.
C
C
C   MACHINE-DEPENDENT PARAMETERS:
C
C      NTERM1 - INTEGER - The no. of terms for the array ARL0.
C                         The recommended value is such that
C                             ABS(ARL0(NTERM1)) < EPS/100
C
C      NTERM2 - INTEGER - The no. of terms for the array ARL0AS.
C                         The recommended value is such that
C                             ABS(ARL0AS(NTERM2)) < EPS/100 
C
C      NTERM3 - INTEGER - The no. of terms for the array AI0ML0.
C                         The recommended value is such that
C                             ABS(AI0ML0(NTERM3)) < EPS/100
C
C      XLOW - DOUBLE PRECISION - The value of x below which L0(x) = 2*x/pi
C                    to machine precision. The recommended value is
C                             3*SQRT(EPS)
C
C      XHIGH1 - DOUBLE PRECISION - The value beyond which the Chebyshev series
C                      in the asymptotic expansion of I0 - L0 gives
C                      1.0 to machine precision. The recommended value
C                      is   SQRT( 30/EPSNEG )
C
C      XHIGH2 - DOUBLE PRECISION - The value beyond which the Chebyshev series
C                      in the asymptotic expansion of I0 gives 1.0
C                      to machine precision. The recommended value
C                      is   28 / EPSNEG
C
C      XMAX - DOUBLE PRECISION - The value of XMAX, where XMAX is the
C                    largest possible floating-pt argument.
C                    This is used to prevent overflow.
C
C      For values of EPS, EPSNEG and XMAX the user should refer
C      to the file MACHCON.TXT
C
C      The machine-arithmetic constants are given in DATA
C      statements.
C
C
C   INTRINSIC FUNCTIONS USED:
C
C      EXP , LOG , SQRT
C
C
C   OTHER MISCFUN SUBROUTINES USED:
C
C          CHEVAL , ERRPRN
C
C
C   AUTHOR:
C          DR. ALLAN J. MACLEOD
C          DEPT. OF MATHEMATICS AND STATISTICS
C          UNIVERSITY OF PAISLEY
C          HIGH ST.
C          PAISLEY
C          SCOTLAND
C          PA1 2BE
C
C      (e-mail: macl_ms0@paisley.ac.uk )
C
C
C   LATEST REVISION:
C                   12 JANUARY, 1996
C
C
      INTEGER INDSGN,NTERM1,NTERM2,NTERM3
      DOUBLE PRECISION ARL0(0:27),ARL0AS(0:15),AI0ML0(0:23),
     1     ATEHUN,CHEVAL,CH1,CH2,FOUR,LNR2PI,ONE,ONEHUN,
     2     SIXTEN,T,TEST,TWENT4,TWENT8,TWO,TWOBPI,TWO88,
     3     X,XHIGH1,XHIGH2,XLOW,XMAX,XVALUE,XSQ,ZERO
      CHARACTER FNNAME*6,ERRMSG*24
      DATA FNNAME/'STRVL0'/
      DATA ERRMSG/'ARGUMENT CAUSES OVERFLOW'/
      DATA ZERO,ONE,TWO/0.0 D 0 , 1.0 D 0 , 2.0 D 0/
      DATA FOUR,SIXTEN/4.0 D 0 , 16.0 D 0/
      DATA TWENT4,TWENT8,ONEHUN/24.0 D 0 , 28.0 D 0 , 100.0 D 0/
      DATA TWO88,ATEHUN/288.0 D 0 , 800.0 D 0/
      DATA LNR2PI/0.91893 85332 04672 74178 D 0/
      DATA TWOBPI/0.63661 97723 67581 34308 D 0/
      DATA ARL0(0)/  0.42127 45834 99799 24863  D    0/
      DATA ARL0(1)/ -0.33859 53639 12206 12188  D    0/
      DATA ARL0(2)/  0.21898 99481 27107 16064  D    0/
      DATA ARL0(3)/ -0.12349 48282 07131 85712  D    0/
      DATA ARL0(4)/  0.62142 09793 86695 8440   D   -1/
      DATA ARL0(5)/ -0.28178 06028 10954 7545   D   -1/
      DATA ARL0(6)/  0.11574 19676 63809 1209   D   -1/
      DATA ARL0(7)/ -0.43165 85743 06921 179    D   -2/
      DATA ARL0(8)/  0.14614 23499 07298 329    D   -2/
      DATA ARL0(9)/ -0.44794 21180 54614 78     D   -3/
      DATA ARL0(10)/ 0.12364 74610 59437 61     D   -3/
      DATA ARL0(11)/-0.30490 28334 79704 4      D   -4/
      DATA ARL0(12)/ 0.66394 14015 21146        D   -5/
      DATA ARL0(13)/-0.12553 83577 03889        D   -5/
      DATA ARL0(14)/ 0.20073 44645 1228         D   -6/
      DATA ARL0(15)/-0.25882 60170 637          D   -7/
      DATA ARL0(16)/ 0.24114 37427 58           D   -8/
      DATA ARL0(17)/-0.10159 67435 2            D   -9/
      DATA ARL0(18)/-0.12024 30736              D  -10/
      DATA ARL0(19)/ 0.26290 6137               D  -11/
      DATA ARL0(20)/-0.15313 190                D  -12/
      DATA ARL0(21)/-0.15747 60                 D  -13/
      DATA ARL0(22)/ 0.31563 5                  D  -14/
      DATA ARL0(23)/-0.4096                     D  -16/
      DATA ARL0(24)/-0.3620                     D  -16/
      DATA ARL0(25)/ 0.239                      D  -17/
      DATA ARL0(26)/ 0.36                       D  -18/
      DATA ARL0(27)/-0.4                        D  -19/
      DATA ARL0AS(0)/  2.00861 30823 56058 88600  D   0/
      DATA ARL0AS(1)/  0.40373 79665 00438 470    D  -2/
      DATA ARL0AS(2)/ -0.25199 48028 65802 67     D  -3/
      DATA ARL0AS(3)/  0.16057 36682 81117 6      D  -4/
      DATA ARL0AS(4)/ -0.10369 21824 73444        D  -5/
      DATA ARL0AS(5)/  0.67655 78876 305          D  -7/
      DATA ARL0AS(6)/ -0.44499 99067 56           D  -8/
      DATA ARL0AS(7)/  0.29468 88922 8            D  -9/
      DATA ARL0AS(8)/ -0.19621 80522              D -10/
      DATA ARL0AS(9)/  0.13133 0306               D -11/
      DATA ARL0AS(10)/-0.88191 90                 D -13/
      DATA ARL0AS(11)/ 0.59537 6                  D -14/
      DATA ARL0AS(12)/-0.40389                    D -15/
      DATA ARL0AS(13)/ 0.2651                     D -16/
      DATA ARL0AS(14)/-0.208                      D -17/
      DATA ARL0AS(15)/ 0.11                       D -18/
      DATA AI0ML0(0)/ 2.00326 51024 11606 43125  D    0/
      DATA AI0ML0(1)/ 0.19520 68515 76492 081    D   -2/
      DATA AI0ML0(2)/ 0.38239 52356 99083 28     D   -3/
      DATA AI0ML0(3)/ 0.75342 80817 05443 6      D   -4/
      DATA AI0ML0(4)/ 0.14959 57655 89707 8      D   -4/
      DATA AI0ML0(5)/ 0.29994 05312 10557        D   -5/
      DATA AI0ML0(6)/ 0.60769 60482 2459         D   -6/
      DATA AI0ML0(7)/ 0.12399 49554 4506         D   -6/
      DATA AI0ML0(8)/ 0.25232 62552 649          D   -7/
      DATA AI0ML0(9)/ 0.50463 48573 32           D   -8/
      DATA AI0ML0(10)/0.97913 23623 0            D   -9/
      DATA AI0ML0(11)/0.18389 11524 1            D   -9/
      DATA AI0ML0(12)/0.33763 09278              D  -10/
      DATA AI0ML0(13)/0.61117 9703               D  -11/
      DATA AI0ML0(14)/0.10847 2972               D  -11/
      DATA AI0ML0(15)/0.18861 271                D  -12/
      DATA AI0ML0(16)/0.32803 45                 D  -13/
      DATA AI0ML0(17)/0.56564 7                  D  -14/
      DATA AI0ML0(18)/0.93300                    D  -15/
      DATA AI0ML0(19)/0.15881                    D  -15/
      DATA AI0ML0(20)/0.2791                     D  -16/
      DATA AI0ML0(21)/0.389                      D  -17/
      DATA AI0ML0(22)/0.70                       D  -18/
      DATA AI0ML0(23)/0.16                       D  -18/
C
C   MACHINE-DEPENDENT VALUES (Suitable for IEEE-arithmetic machines)
C
      DATA NTERM1,NTERM2,NTERM3/25,14,21/
      DATA XLOW,XMAX/4.4703484D-8,1.797693D308/
      DATA XHIGH1,XHIGH2/5.1982303D8,2.5220158D17/
C
C   Start computation
C
      X = XVALUE
      INDSGN = 1
      IF ( X .LT. ZERO ) THEN
         X = -X
         INDSGN = -1
      ENDIF
C
C   Code for |xvalue| <= 16
C
      IF ( X .LE. SIXTEN ) THEN
         IF ( X .LT. XLOW ) THEN
            STRVL0 = TWOBPI * X
         ELSE
            T = ( FOUR * X - TWENT4 ) / ( X + TWENT4 )
            STRVL0 = TWOBPI * X * CHEVAL(NTERM1,ARL0,T) * EXP(X)
         ENDIF
      ELSE
C
C   Code for |xvalue| > 16
C
         IF ( X .GT. XHIGH2 ) THEN
            CH1 = ONE
         ELSE
            T = ( X - TWENT8 ) / ( FOUR - X )
            CH1 = CHEVAL(NTERM2,ARL0AS,T)
         ENDIF
         IF ( X .GT. XHIGH1 ) THEN
            CH2 = ONE
         ELSE
            XSQ = X * X
            T = ( ATEHUN - XSQ ) / ( TWO88 + XSQ )
            CH2 = CHEVAL(NTERM3,AI0ML0,T)
         ENDIF
         TEST = LOG(CH1) - LNR2PI - LOG(X)/TWO + X
         IF ( TEST .GT. LOG(XMAX) ) THEN
            CALL ERRPRN(FNNAME,ERRMSG)
            STRVL0 = XMAX
         ELSE
            STRVL0 = EXP(TEST) - TWOBPI * CH2 / X
         ENDIF
      ENDIF
      IF ( INDSGN .EQ. -1 ) STRVL0 = -STRVL0
      RETURN
      END

      DOUBLE PRECISION FUNCTION STRVL1(XVALUE)

c*********************************************************************72
c
cc STRVL1 calculates the modified Struve function of order 1.
c
C   DESCRIPTION:
C
C      This function calculates the modified Struve function of
C      order 1, denoted L1(x), defined as the solution of
C
C               x*x*D(Df) + x*Df - (x*x+1)f = 2*x*x/pi
C
C      This subroutine is set up to work on IEEE machines.
C      For other machines, you should retrieve the code
C      from the general MISCFUN archive.
C
C
C   ERROR RETURNS:
C
C      If the value of |XVALUE| is too large, the result
C      would cause an floating-pt overflow. An error message
C      is printed and the function returns the value of
C      sign(XVALUE)*XMAX where XMAX is the largest possible
C      floating-pt argument.
C
C
C   MACHINE-DEPENDENT PARAMETERS:
C
C      NTERM1 - INTEGER - The no. of terms for the array ARL1.
C                         The recommended value is such that
C                             ABS(ARL1(NTERM1)) < EPS/100
C
C      NTERM2 - INTEGER - The no. of terms for the array ARL1AS.
C                         The recommended value is such that
C                             ABS(ARL1AS(NTERM2)) < EPS/100 
C
C      NTERM3 - INTEGER - The no. of terms for the array AI1ML1.
C                         The recommended value is such that
C                             ABS(AI1ML1(NTERM3)) < EPS/100
C
C      XLOW1 - DOUBLE PRECISION - The value of x below which 
C                                     L1(x) = 2*x*x/(3*pi)
C                                 to machine precision. The recommended 
C                                 value is     SQRT(15*EPS)
C
C      XLOW2 - DOUBLE PRECISION - The value of x below which L1(x) set to 0.0.
C                     This is used to prevent underflow. The
C                     recommended value is
C                              SQRT(5*XMIN)
C
C      XHIGH1 - DOUBLE PRECISION - The value of |x| above which the Chebyshev
C                      series in the asymptotic expansion of I1
C                      equals 1.0 to machine precision. The
C                      recommended value is  SQRT( 30 / EPSNEG ).
C
C      XHIGH2 - DOUBLE PRECISION - The value of |x| above which the Chebyshev
C                      series in the asymptotic expansion of I1 - L1
C                      equals 1.0 to machine precision. The recommended
C                      value is   30 / EPSNEG.
C 
C      XMAX - DOUBLE PRECISION - The value of XMAX, where XMAX is the
C                    largest possible floating-pt argument.
C                    This is used to prevent overflow.
C
C      For values of EPS, EPSNEG, XMIN, and XMAX the user should refer 
C      to the file MACHCON.TXT
C
C      The machine-arithmetic constants are given in DATA
C      statements.
C
C
C   INTRINSIC FUNCTIONS USED:
C
C      EXP , LOG , SQRT
C
C
C   OTHER MISCFUN SUBROUTINES USED:
C
C          CHEVAL , ERRPRN
C
C
C   AUTHOR:
C          DR. ALLAN J. MACLEOD
C          DEPT. OF MATHEMATICS AND STATISTICS
C          UNIVERSITY OF PAISLEY
C          HIGH ST.
C          PAISLEY
C          SCOTLAND
C          PA1 2BE
C
C          (e-mail: macl_ms0@paisley.ac.uk )
C
C
C   LATEST UPDATE:
C                 12 JANUARY, 1996
C
C
      INTEGER NTERM1,NTERM2,NTERM3
      DOUBLE PRECISION ARL1(0:26),ARL1AS(0:16),AI1ML1(0:25),
     1     ATEHUN,CHEVAL,CH1,CH2,FOUR,LNR2PI,
     2     ONE,ONEHUN,PI3BY2,SIXTEN,T,TEST,THIRTY,TWENT4,
     3     TWO,TWOBPI,TWO88,X,XHIGH1,XHIGH2,XLOW1,XLOW2,
     4     XMAX,XVALUE,XSQ,ZERO
      CHARACTER FNNAME*6,ERRMSG*24
      DATA FNNAME/'STRVL1'/
      DATA ERRMSG/'ARGUMENT CAUSES OVERFLOW'/
      DATA ZERO,ONE,TWO/0.0 D 0 , 1.0 D 0 , 2.0 D 0/
      DATA FOUR,SIXTEN/4.0 D 0 , 16.0 D 0/
      DATA TWENT4,THIRTY/24.0 D 0 , 30.0 D 0/
      DATA ONEHUN/100.0 D 0/
      DATA TWO88,ATEHUN/288.0 D 0 , 800.0 D 0/
      DATA LNR2PI/0.91893 85332 04672 74178 D 0/
      DATA PI3BY2/4.71238 89803 84689 85769 D 0/
      DATA TWOBPI/0.63661 97723 67581 34308 D 0/
      DATA ARL1(0)/  0.38996 02735 12295 38208  D    0/
      DATA ARL1(1)/ -0.33658 09610 19757 49366  D    0/
      DATA ARL1(2)/  0.23012 46791 25016 45616  D    0/
      DATA ARL1(3)/ -0.13121 59400 79608 32327  D    0/
      DATA ARL1(4)/  0.64259 22289 91284 6518   D   -1/
      DATA ARL1(5)/ -0.27500 32950 61663 5833   D   -1/
      DATA ARL1(6)/  0.10402 34148 63720 8871   D   -1/
      DATA ARL1(7)/ -0.35053 22949 36388 080    D   -2/
      DATA ARL1(8)/  0.10574 84984 21439 717    D   -2/
      DATA ARL1(9)/ -0.28609 42640 36665 58     D   -3/
      DATA ARL1(10)/ 0.69257 08785 94220 8      D   -4/
      DATA ARL1(11)/-0.14896 93951 12271 7      D   -4/
      DATA ARL1(12)/ 0.28103 55825 97128        D   -5/
      DATA ARL1(13)/-0.45503 87929 7776         D   -6/
      DATA ARL1(14)/ 0.60901 71561 770          D   -7/
      DATA ARL1(15)/-0.62354 37248 08           D   -8/
      DATA ARL1(16)/ 0.38430 01206 7            D   -9/
      DATA ARL1(17)/ 0.79054 3916               D  -11/
      DATA ARL1(18)/-0.48982 4083               D  -11/
      DATA ARL1(19)/ 0.46356 884                D  -12/
      DATA ARL1(20)/ 0.68420 5                  D  -14/
      DATA ARL1(21)/-0.56974 8                  D  -14/
      DATA ARL1(22)/ 0.35324                    D  -15/
      DATA ARL1(23)/ 0.4244                     D  -16/
      DATA ARL1(24)/-0.644                      D  -17/
      DATA ARL1(25)/-0.21                       D  -18/
      DATA ARL1(26)/ 0.9                        D  -19/
      DATA ARL1AS(0)/  1.97540 37844 16523 56868  D    0/
      DATA ARL1AS(1)/ -0.11951 30555 08829 4181   D   -1/
      DATA ARL1AS(2)/  0.33639 48526 91960 46     D   -3/
      DATA ARL1AS(3)/ -0.10091 15655 48154 9      D   -4/
      DATA ARL1AS(4)/  0.30638 95132 1998         D   -6/
      DATA ARL1AS(5)/ -0.95370 43703 96           D   -8/
      DATA ARL1AS(6)/  0.29524 73555 8            D   -9/
      DATA ARL1AS(7)/ -0.95107 8318               D  -11/
      DATA ARL1AS(8)/  0.28203 667                D  -12/
      DATA ARL1AS(9)/ -0.11341 75                 D  -13/
      DATA ARL1AS(10)/ 0.147                      D  -17/
      DATA ARL1AS(11)/-0.6232                     D  -16/
      DATA ARL1AS(12)/-0.751                      D  -17/
      DATA ARL1AS(13)/-0.17                       D  -18/
      DATA ARL1AS(14)/ 0.51                       D  -18/
      DATA ARL1AS(15)/ 0.23                       D  -18/
      DATA ARL1AS(16)/ 0.5                        D  -19/
      DATA AI1ML1(0)/  1.99679 36189 67891 36501  D    0/
      DATA AI1ML1(1)/ -0.19066 32614 09686 132    D   -2/
      DATA AI1ML1(2)/ -0.36094 62241 01744 81     D   -3/
      DATA AI1ML1(3)/ -0.68418 47304 59982 0      D   -4/
      DATA AI1ML1(4)/ -0.12990 08228 50942 6      D   -4/
      DATA AI1ML1(5)/ -0.24715 21887 05765        D   -5/
      DATA AI1ML1(6)/ -0.47147 83969 1972         D   -6/
      DATA AI1ML1(7)/ -0.90208 19982 592          D   -7/
      DATA AI1ML1(8)/ -0.17304 58637 504          D   -7/
      DATA AI1ML1(9)/ -0.33232 36701 59           D   -8/
      DATA AI1ML1(10)/-0.63736 42173 5            D   -9/
      DATA AI1ML1(11)/-0.12180 23975 6            D   -9/
      DATA AI1ML1(12)/-0.23173 46832              D  -10/
      DATA AI1ML1(13)/-0.43906 8833               D  -11/
      DATA AI1ML1(14)/-0.82847 110                D  -12/
      DATA AI1ML1(15)/-0.15562 249                D  -12/
      DATA AI1ML1(16)/-0.29131 12                 D  -13/
      DATA AI1ML1(17)/-0.54396 5                  D  -14/
      DATA AI1ML1(18)/-0.10117 7                  D  -14/
      DATA AI1ML1(19)/-0.18767                    D  -15/
      DATA AI1ML1(20)/-0.3484                     D  -16/
      DATA AI1ML1(21)/-0.643                      D  -17/
      DATA AI1ML1(22)/-0.118                      D  -17/
      DATA AI1ML1(23)/-0.22                       D  -18/
      DATA AI1ML1(24)/-0.4                        D  -19/
      DATA AI1ML1(25)/-0.1                        D  -19/
C
C   MACHINE-DEPENDENT VALUES (Suitable for IEEE-arithmetic machines)
C
      DATA NTERM1,NTERM2,NTERM3/24,13,22/
      DATA XLOW1,XLOW2,XMAX/5.7711949D-8,3.3354714D-154,1.797693D308/
      DATA XHIGH1,XHIGH2/5.19823025D8,2.7021597D17/
C
C   START CALCULATION
C
      X = ABS ( XVALUE )
C
C   CODE FOR |XVALUE| <= 16
C
      IF ( X .LE. SIXTEN ) THEN
         IF ( X .LE. XLOW2 ) THEN
            STRVL1 = ZERO
         ELSE 
            XSQ = X * X
            IF ( X .LT. XLOW1 ) THEN
               STRVL1 = XSQ / PI3BY2
            ELSE
               T = ( FOUR * X - TWENT4 ) / ( X + TWENT4 )
               STRVL1 = XSQ * CHEVAL(NTERM1,ARL1,T) * EXP(X) / PI3BY2
            ENDIF
         ENDIF
      ELSE
C
C   CODE FOR |XVALUE| > 16
C
         IF ( X .GT. XHIGH2 ) THEN
            CH1 = ONE
         ELSE
            T = ( X - THIRTY ) / ( TWO - X )
            CH1 = CHEVAL(NTERM2,ARL1AS,T)
         ENDIF
         IF ( X .GT. XHIGH1 ) THEN
            CH2 = ONE
         ELSE
            XSQ = X * X
            T = ( ATEHUN - XSQ ) / ( TWO88 + XSQ )
            CH2 = CHEVAL(NTERM3,AI1ML1,T)
         ENDIF
         TEST = LOG(CH1) - LNR2PI - LOG(X)/TWO + X
         IF ( TEST .GT. LOG(XMAX) ) THEN
            CALL ERRPRN(FNNAME,ERRMSG)
            STRVL1 = XMAX
         ELSE
            STRVL1 = EXP(TEST) - TWOBPI * CH2
         ENDIF
      ENDIF
      RETURN
      END

      DOUBLE PRECISION FUNCTION SYNCH1(XVALUE)

c*********************************************************************72
c
cc SYNCH1 calculates the synchrotron radiation function.
c
C   DESCRIPTION:
C
C      This function calculates the synchrotron radiation function
C      defined as
C
C         SYNCH1(x) = x * Integral{x to inf} K(5/3)(t) dt,
C
C      where K(5/3) is a modified Bessel function of order 5/3.
C
C      The code uses Chebyshev expansions, the coefficients of which 
C      are given to 20 decimal places.
C
C      This subroutine is set up to work on IEEE machines.
C      For other machines, you should retrieve the code
C      from the general MISCFUN archive.
C
C
C   ERROR RETURNS:
C
C      The function is undefined if x < 0.0. If XVALUE < 0.0,
C      an error message is printed and the function returns
C      the value 0.0.
C
C
C   MACHINE-DEPENDENT CONSTANTS:
C
C      NTERM1 - INTEGER - The no. of terms needed from the array
C                         ASYNC1. The recommended value is such that
C                            ABS(ASYNC1(NTERM1)) < EPS/100.
C
C      NTERM2 - INTEGER - The no. of terms needed from the array
C                         ASYNC2. The recommended value is such that
C                            ABS(ASYNC2(NTERM2)) < EPS/100.
C
C      NTERM3 - INTEGER - The no. of terms needed from the array
C                         ASYNCA. The recommended value is such that
C                            ABS(ASYNCA(NTERM3)) < EPS/100.
C
C      XLOW - DOUBLE PRECISION - The value below which
C                        SYNCH1(x) = 2.14952.. * (x**(1/3))
C                    to machine precision. The recommended value
C                    is     sqrt (8*EPSNEG)
C
C      XHIGH1 - DOUBLE PRECISION - The value above which
C                          SYNCH1(x) = 0.0
C                      to machine precision. The recommended value
C                      is     -8*LN(XMIN)/7
C
C      XHIGH2 - DOUBLE PRECISION - The value of LN(XMIN). This is used
C                      to prevent underflow in calculations
C                      for large x.
C
C     For values of EPS, EPSNEG, and XMIN refer to the file MACHCON.TXT
C
C     The machine-arithmetic constants are given in DATA
C     statements.
C
C
C   INTRINSIC FUNCTIONS USED:
C      
C      EXP , LOG , SQRT
C
C
C   OTHER MISCFUN SUBROUTINES USED:
C
C          CHEVAL , ERRPRN
C
C
C   AUTHOR:
C          Dr. Allan J. MacLeod,
C          Dept. of Mathematics and Statistics,
C          University of Paisley,
C          Paisley,
C          SCOTLAND
C          PA1 2BE
C
C          ( e-mail:  macl_ms0@paisley.ac.uk )
C
C
C   LATEST UPDATE:
C                 12 JANUARY, 1996
C
      INTEGER NTERM1,NTERM2,NTERM3
      DOUBLE PRECISION ASYNC1(0:13),ASYNC2(0:11),ASYNCA(0:24),
     1     CHEB1,CHEB2,CHEVAL,CONLOW,EIGHT,FOUR,HALF,
     2     LNRTP2,ONE,ONEHUN,PIBRT3,T,THREE,TWELVE,X,XHIGH1,
     3     XHIGH2,XLOW,XPOWTH,XVALUE,ZERO
      CHARACTER FNNAME*6,ERRMSG*14
      DATA FNNAME/'SYNCH1'/
      DATA ERRMSG/'ARGUMENT < 0.0'/
      DATA ZERO,HALF,ONE/ 0.0 D 0 , 0.5 D 0 , 1.0 D 0 /
      DATA THREE,FOUR/ 3.0 D 0 , 4.0 D 0 /
      DATA EIGHT,TWELVE/ 8.0 D 0 , 12.0 D 0 /
      DATA ONEHUN/ 100.0 D 0 /
      DATA CONLOW/2.14952 82415 34478 63671 D 0/
      DATA PIBRT3/1.81379 93642 34217 85059 D 0/
      DATA LNRTP2/0.22579 13526 44727 43236 D 0/
      DATA ASYNC1/30.36468 29825 01076 27340  D    0,
     1            17.07939 52774 08394 57449  D    0,
     2             4.56013 21335 45072 88887  D    0,
     3             0.54928 12467 30419 97963  D    0,
     4             0.37297 60750 69301 1724   D   -1,
     5             0.16136 24302 01041 242    D   -2,
     6             0.48191 67721 20370 7      D   -4,
     7             0.10512 42528 89384        D   -5,
     8             0.17463 85046 697          D   -7,
     9             0.22815 48654 4            D   -9,
     X             0.24044 3082               D  -11,
     1             0.20865 88                 D  -13,
     2             0.15167                    D  -15,
     3             0.94                       D  -18/
      DATA ASYNC2/0.44907 21623 53266 08443  D    0,
     1            0.89835 36779 94187 2179   D   -1,
     2            0.81044 57377 21512 894    D   -2,
     3            0.42617 16991 08916 19     D   -3,
     4            0.14760 96312 70746 0      D   -4,
     5            0.36286 33615 3998         D   -6,
     6            0.66634 80749 84           D   -8,
     7            0.94907 71655              D  -10,
     8            0.10791 2491               D  -11,
     9            0.10022 01                 D  -13,
     X            0.7745                     D  -16,
     1            0.51                       D  -18/
      DATA ASYNCA(0)/ 2.13293 05161 35500 09848  D    0/
      DATA ASYNCA(1)/ 0.74135 28649 54200 2401   D   -1/
      DATA ASYNCA(2)/ 0.86968 09990 99641 978    D   -2/
      DATA ASYNCA(3)/ 0.11703 82624 87756 921    D   -2/
      DATA ASYNCA(4)/ 0.16451 05798 61919 15     D   -3/
      DATA ASYNCA(5)/ 0.24020 10214 20640 3      D   -4/
      DATA ASYNCA(6)/ 0.35827 75638 93885        D   -5/
      DATA ASYNCA(7)/ 0.54477 47626 9837         D   -6/
      DATA ASYNCA(8)/ 0.83880 28561 957          D   -7/
      DATA ASYNCA(9)/ 0.13069 88268 416          D   -7/
      DATA ASYNCA(10)/0.20530 99071 44           D   -8/
      DATA ASYNCA(11)/0.32518 75368 8            D   -9/
      DATA ASYNCA(12)/0.51791 40412              D  -10/
      DATA ASYNCA(13)/0.83002 9881               D  -11/
      DATA ASYNCA(14)/0.13352 7277               D  -11/
      DATA ASYNCA(15)/0.21591 498                D  -12/
      DATA ASYNCA(16)/0.34996 73                 D  -13/
      DATA ASYNCA(17)/0.56994 2                  D  -14/
      DATA ASYNCA(18)/0.92906                    D  -15/
      DATA ASYNCA(19)/0.15222                    D  -15/
      DATA ASYNCA(20)/0.2491                     D  -16/
      DATA ASYNCA(21)/0.411                      D  -17/
      DATA ASYNCA(22)/0.67                       D  -18/
      DATA ASYNCA(23)/0.11                       D  -18/
      DATA ASYNCA(24)/0.2                        D  -19/
C
C   Machine-dependent constants (suitable for IEEE machines)
C
      DATA NTERM1,NTERM2,NTERM3/12,10,21/
      DATA XLOW/2.98023224D-8/
      DATA XHIGH1,XHIGH2/809.595907D0,-708.396418D0/
C
C   Start calculation
C 
      X = XVALUE
      IF ( X .LT. ZERO ) THEN
         CALL ERRPRN(FNNAME,ERRMSG)
         SYNCH1 = ZERO
         RETURN
      ENDIF
C
C   Code for 0 <= x <= 4
C
      IF ( X .LE. FOUR ) THEN
         XPOWTH = X ** ( ONE / THREE )
         IF ( X .LT. XLOW ) THEN
            SYNCH1 = CONLOW * XPOWTH
         ELSE
            T = ( X * X / EIGHT - HALF ) - HALF
            CHEB1 = CHEVAL(NTERM1,ASYNC1,T)
            CHEB2 = CHEVAL(NTERM2,ASYNC2,T)
            T = XPOWTH * CHEB1 - ( XPOWTH**11 ) * CHEB2
            SYNCH1 = T - PIBRT3 * X
         ENDIF
      ELSE
         IF ( X .GT. XHIGH1 ) THEN
            SYNCH1 = ZERO
         ELSE
            T = ( TWELVE - X ) / ( X + FOUR )
            CHEB1 = CHEVAL(NTERM3,ASYNCA,T)
            T = LNRTP2 - X + LOG( SQRT(X) * CHEB1 )
            IF ( T .LT. XHIGH2 ) THEN
               SYNCH1 = ZERO
            ELSE
               SYNCH1 = EXP(T)
            ENDIF
         ENDIF
      ENDIF
      RETURN
      END

      DOUBLE PRECISION FUNCTION SYNCH2(XVALUE)

c*********************************************************************72
c
cc SYNCH2 calculates the synchrotron radiation function.
c
C   DESCRIPTION:
C
C      This function calculates the synchrotron radiation function
C      defined as
C
C         SYNCH2(x) = x * K(2/3)(x)
C
C      where K(2/3) is a modified Bessel function of order 2/3.
C
C      The code uses Chebyshev expansions, the coefficients of which 
C      are given to 20 decimal places.
C
C      This subroutine is set up to work on IEEE machines.
C      For other machines, you should retrieve the code
C      from the general MISCFUN archive.
C
C
C   ERROR RETURNS:
C
C      The function is undefined if x < 0.0. If XVALUE < 0.0,
C      an error message is printed and the function returns
C      the value 0.0.
C
C
C   MACHINE-DEPENDENT CONSTANTS:
C
C      NTERM1 - INTEGER - The no. of terms needed from the array
C                         ASYNC1. The recommended value is such that
C                            ABS(ASYN21(NTERM1)) < EPS/100.
C
C      NTERM2 - INTEGER - The no. of terms needed from the array
C                         ASYNC2. The recommended value is such that
C                            ABS(ASYN22(NTERM2)) < EPS/100.
C
C      NTERM3 - INTEGER - The no. of terms needed from the array
C                         ASYNCA. The recommended value is such that
C                            ABS(ASYN2A(NTERM3)) < EPS/100.
C
C      XLOW - DOUBLE PRECISION - The value below which
C                        SYNCH2(x) = 1.074764... * (x**(1/3))
C                    to machine precision. The recommended value
C                    is     sqrt (8*EPSNEG)
C
C      XHIGH1 - DOUBLE PRECISION - The value above which
C                          SYNCH2(x) = 0.0
C                      to machine precision. The recommended value
C                      is     -8*LN(XMIN)/7
C
C      XHIGH2 - DOUBLE PRECISION - The value of LN(XMIN). This is used
C                      to prevent underflow in calculations
C                      for large x.
C
C     For values of EPS, EPSNEG, and XMIN refer to the file MACHCON.TXT
C
C     The machine-arithmetic constants are given in DATA
C     statements.
C
C
C   INTRINSIC FUNCTIONS USED:
C
C      EXP , LOG , SQRT
C
C
C   OTHER MISCFUN SUBROUTINES USED:
C
C          CHEVAL , ERRPRN
C
C
C   AUTHOR:
C          Dr. Allan J. MacLeod,
C          Dept. of Mathematics and Statistics,
C          University of Paisley,
C          Paisley,
C          SCOTLAND
C          PA1 2BE
C
C          (e-mail:  macl_ms0@paisley.ac.uk)
C
C
C   LATEST UPDATE:
C                 12 JANUARY, 1996
C
      INTEGER NTERM1,NTERM2,NTERM3
      DOUBLE PRECISION ASYN21(0:14),ASYN22(0:13),ASYN2A(0:18),
     1     CHEB1,CHEB2,CHEVAL,CONLOW,EIGHT,FOUR,HALF,
     2     LNRTP2,ONE,ONEHUN,T,TEN,THREE,TWO,X,XHIGH1,
     3     XHIGH2,XLOW,XPOWTH,XVALUE,ZERO
      CHARACTER FNNAME*6,ERRMSG*14
      DATA FNNAME/'SYNCH2'/
      DATA ERRMSG/'ARGUMENT < 0.0'/
      DATA ZERO,HALF,ONE/ 0.0 D 0 , 0.5 D 0 , 1.0 D 0 /
      DATA TWO,THREE,FOUR/ 2.0 D 0 , 3.0 D 0 , 4.0 D 0 /
      DATA EIGHT,TEN,ONEHUN/ 8.0 D 0 , 10.0 D 0 , 100.0 D 0/
      DATA CONLOW/1.07476 41207 67239 31836 D 0/
      DATA LNRTP2/0.22579 13526 44727 43236 D 0/
      DATA ASYN21/38.61783 99238 43085 48014  D    0,
     1            23.03771 55949 63734 59697  D    0,
     2             5.38024 99868 33570 59676  D    0,
     3             0.61567 93806 99571 07760  D    0,
     4             0.40668 80046 68895 5843   D   -1,
     5             0.17296 27455 26484 141    D   -2,
     6             0.51061 25883 65769 9      D   -4,
     7             0.11045 95950 22012        D   -5,
     8             0.18235 53020 649          D   -7,
     9             0.23707 69803 4            D   -9,
     X             0.24887 2963               D  -11,
     1             0.21528 68                 D  -13,
     2             0.15607                    D  -15,
     3             0.96                       D  -18,
     4             0.1                        D  -19/
      DATA ASYN22/7.90631 48270 66080 42875  D    0,
     1            3.13534 63612 85342 56841  D    0,
     2            0.48548 79477 45371 45380  D    0,
     3            0.39481 66758 27237 2337   D   -1,
     4            0.19661 62233 48088 022    D   -2,
     5            0.65907 89322 93042 0      D   -4,
     6            0.15857 56134 98559        D   -5,
     7            0.28686 53011 233          D   -7,
     8            0.40412 02359 5            D   -9,
     9            0.45568 4443               D  -11,
     X            0.42045 90                 D  -13,
     1            0.32326                    D  -15,
     2            0.210                      D  -17,
     3            0.1                        D  -19/
      DATA ASYN2A/2.02033 70941 70713 60032  D    0,
     1            0.10956 23712 18074 0443   D   -1,
     2            0.85423 84730 11467 55     D   -3,
     3            0.72343 02421 32822 2      D   -4,
     4            0.63124 42796 26992        D   -5,
     5            0.56481 93141 1744         D   -6,
     6            0.51283 24801 375          D   -7,
     7            0.47196 53291 45           D   -8,
     8            0.43807 44214 3            D   -9,
     9            0.41026 81493              D  -10,
     X            0.38623 0721               D  -11,
     1            0.36613 228                D  -12,
     2            0.34802 32                 D  -13,
     3            0.33301 0                  D  -14,
     4            0.31856                    D  -15,
     5            0.3074                     D  -16,
     6            0.295                      D  -17,
     7            0.29                       D  -18,
     8            0.3                        D  -19/
C
C   Machine-dependent constants (suitable for IEEE machines)
C
      DATA NTERM1,NTERM2,NTERM3/13,12,16/
      DATA XLOW/2.98023224D-8/
      DATA XHIGH1,XHIGH2/809.595907D0,-708.396418D0/
C
C   Start calculation
C 
      X = XVALUE
      IF ( X .LT. ZERO ) THEN
         CALL ERRPRN(FNNAME,ERRMSG)
         SYNCH2 = ZERO
         RETURN
      ENDIF
C
C   Code for 0 <= x <= 4
C
      IF ( X .LE. FOUR ) THEN
         XPOWTH = X ** ( ONE / THREE )
         IF ( X .LT. XLOW ) THEN
            SYNCH2 = CONLOW * XPOWTH
         ELSE
            T = ( X * X / EIGHT - HALF ) - HALF
            CHEB1 = CHEVAL(NTERM1,ASYN21,T)
            CHEB2 = CHEVAL(NTERM2,ASYN22,T)
            SYNCH2 = XPOWTH * CHEB1 - ( XPOWTH**5 ) * CHEB2
         ENDIF
      ELSE
         IF ( X .GT. XHIGH1 ) THEN
            SYNCH2 = ZERO
         ELSE
            T = ( TEN - X ) / ( X + TWO )
            CHEB1 = CHEVAL(NTERM3,ASYN2A,T)
            T = LNRTP2 - X + LOG( SQRT(X) * CHEB1 )
            IF ( T .LT. XHIGH2 ) THEN
               SYNCH2 = ZERO
            ELSE
               SYNCH2 = EXP(T)
            ENDIF
         ENDIF
      ENDIF
      RETURN
      END

      DOUBLE PRECISION FUNCTION TRAN02(XVALUE)

c*********************************************************************72
c
cc TRAN02 calculates the transport integral of order 2.
c
C  DESCRIPTION:
C
C      This program calculates the transport integral of order 2, defined as
C
C        TRAN02(X) = integral 0 to X { t**2 exp(t)/[exp(t)-1]**2 } dt
C
C      The program uses a Chebyshev series, the coefficients of which are
C      given to an accuracy of 20 decimal places.
C
C      This subroutine is set up to work on IEEE machines.
C      For other machines, you should retrieve the code
C      from the general MISCFUN archive.
C
C
C  ERROR RETURNS:
C
C    If XVALUE < 0.0, an error message is printed, and the program 
C    returns the value 0.0. 
C
C
C  MACHINE-DEPENDENT CONSTANTS:
C
C    NTERMS - INTEGER - The number of terms of the array ATRAN to be used.
C                       The recommended value is such that
C                             ATRAN(NTERMS) < EPS/100
C
C    XLOW1 - DOUBLE PRECISION - The value below which TRAN02 = x to
C                   machine precision. The recommended value is
C                             sqrt(8*EPSNEG)
C
C    XHIGH1 - DOUBLE PRECISION - The value above which the exponential series for
C                    large x contains only one term. The recommended value
C                    is        - ln(EPS).
C
C    XHIGH2 - DOUBLE PRECISION - The value above which 
C                       TRAN02 = VALINF  -  x**2 exp(-x)
C                    The recommended value is 2/EPS
C
C    XHIGH3 - DOUBLE PRECISION - The value of ln(EPSNEG). Used to prevent overflow
C                    for large x.
C
C     For values of EPS, EPSNEG, and XMIN refer to the file MACHCON.TXT
C
C     The machine-arithmetic constants are given in DATA
C     statements.
C
C
C  INTRINSIC FUNCTIONS USED:
C
C     EXP, INT, LOG, SQRT
C
C
C   OTHER MISCFUN SUBROUTINES USED:
C
C          CHEVAL , ERRPRN
C
C
C  AUTHOR:
C
C     DR. ALLAN J. MACLEOD,
C     DEPT. OF MATHEMATICS AND STATISTICS,
C     UNIVERSITY OF PAISLEY ,
C     HIGH ST.,
C     PAISLEY,
C     SCOTLAND.
C     PA1 2BE.
C
C     (e-mail: macl_ms0@paisley.ac.uk )
C
C
C  LATEST REVISION:     12 JANUARY, 1996
C
C
      INTEGER K1,K2,NTERMS,NUMEXP,NUMJN
      DOUBLE PRECISION ATRAN(0:19),CHEVAL,EIGHT,FOUR,HALF,ONE,ONEHUN,RK,
     &     RNUMJN,SUMEXP,SUM2,T,VALINF,X,XHIGH1,XHIGH2,
     &     XHIGH3,XK,XK1,XLOW1,XVALUE,ZERO
      CHARACTER FNNAME*6,ERRMSG*14
      DATA FNNAME/'TRAN02'/
      DATA ERRMSG/'ARGUMENT < 0.0'/
      DATA ZERO,HALF,ONE/ 0.0 D 0 , 0.5 D 0 , 1.0 D 0 /
      DATA FOUR,EIGHT,ONEHUN/ 4.0 D 0 , 8.0 D 0 , 100.0 D 0 /
      DATA NUMJN,RNUMJN/ 2 , 2.0 D 0 /
      DATA VALINF/0.32898 68133 69645 28729 D 1/
      DATA ATRAN/1.67176 04464 34538 50301  D    0,
     1          -0.14773 53599 46794 48986  D    0,
     2           0.14821 38199 46936 3384   D   -1,
     3          -0.14195 33032 63056 126    D   -2,
     4           0.13065 41324 41570 83     D   -3,
     5          -0.11715 57958 67579 0      D   -4,
     6           0.10333 49844 57557        D   -5,
     7          -0.90191 13042 227          D   -7,
     8           0.78177 16983 31           D   -8,
     9          -0.67445 65684 0            D   -9,
     X           0.57994 63945              D  -10,
     1          -0.49747 6185               D  -11,
     2           0.42596 097                D  -12,
     3          -0.36421 89                 D  -13,
     4           0.31108 6                  D  -14,
     5          -0.26547                    D  -15,
     6           0.2264                     D  -16,
     7          -0.193                      D  -17,
     8           0.16                       D  -18,
     9          -0.1                        D  -19/
C
C  Machine-dependent constants
C
      DATA NTERMS/17/
      DATA XLOW1/2.98023224D-8/
      DATA XHIGH1,XHIGH3/36.04365668D0,-36.73680056D0/
      DATA XHIGH2/9.00719925D15/
C
C  Start execution
C
      X = XVALUE
C
C  Error test
C
      IF ( X .LT. ZERO ) THEN
         CALL ERRPRN(FNNAME,ERRMSG)
         TRAN02 = ZERO
         RETURN
      ENDIF
C
C   Code for x < =  4.0
C
      IF ( X .LE. FOUR ) THEN
         IF ( X .LT. XLOW1 ) THEN
            TRAN02 =  ( X ** ( NUMJN - 1 ) ) / ( RNUMJN - ONE ) 
         ELSE
            T = ( ( ( X * X ) / EIGHT ) - HALF ) - HALF
            TRAN02 = ( X ** ( NUMJN - 1 ) ) * CHEVAL(NTERMS,ATRAN,T) 
         ENDIF
      ELSE
C 
C  Code for x > 4.0
C
         IF ( X .GT. XHIGH2 ) THEN
            SUMEXP = ONE
         ELSE
            IF ( X .LE. XHIGH1 ) THEN
               NUMEXP = INT ( XHIGH1 / X ) + 1
               T = EXP(-X) 
            ELSE
               NUMEXP = 1
               T = ONE
            ENDIF
            RK = ZERO
            DO 100 K1 = 1 , NUMEXP
               RK = RK + ONE
  100       CONTINUE
            SUMEXP = ZERO
            DO 300 K1 = 1 , NUMEXP
               SUM2 = ONE
               XK = ONE / ( RK * X ) 
               XK1 = ONE
               DO 200 K2 = 1 , NUMJN
                  SUM2 = SUM2 * XK1 * XK + ONE
                  XK1 = XK1 + ONE
  200          CONTINUE
               SUMEXP = SUMEXP * T + SUM2
               RK = RK - ONE
  300       CONTINUE
         ENDIF
         T = RNUMJN * LOG(X) - X + LOG(SUMEXP) 
         IF ( T .LT. XHIGH3 ) THEN
            TRAN02 = VALINF
         ELSE
            TRAN02 = VALINF - EXP(T) 
         ENDIF
      ENDIF
      RETURN
      END

      DOUBLE PRECISION FUNCTION TRAN03(XVALUE)

c*********************************************************************72
c
cc TRAN03 calculates the transport integral of order 3.
c
C  DESCRIPTION:
C
C      This program calculates the transport integral of order 3, defined as
C 
C        TRAN03(X) = integral 0 to X { t**3 exp(t)/[exp(t)-1]**2 } dt
C 
C      The program uses a Chebyshev series, the coefficients of which are
C      given to an accuracy of 20 decimal places.
C
C      This subroutine is set up to work on IEEE machines.
C      For other machines, you should retrieve the code
C      from the general MISCFUN archive.
C
C
C  ERROR RETURNS:
C
C    If XVALUE < 0.0, an error message is printed, and the program 
C    returns the value 0.0. 
C
C
C  MACHINE-DEPENDENT CONSTANTS:
C
C    NTERMS - INTEGER - The number of terms of the array ATRAN to be used.
C                       The recommended value is such that
C                             ATRAN(NTERMS) < EPS/100
C
C    XLOW2 - DOUBLE PRECISION - The value below which TRAN03 = 0.0 to machine
C                    precision. The recommended value is
C                          square root of (2*XMIN)
C
C    XLOW1 - DOUBLE PRECISION - The value below which TRAN03 = X**2/2 to
C                   machine precision. The recommended value is
C                             sqrt(8*EPSNEG)
C
C    XHIGH1 - DOUBLE PRECISION - The value above which the exponential series for
C                    large X contains only one term. The recommended value
C                    is        - ln(EPS).
C
C    XHIGH2 - DOUBLE PRECISION - The value above which 
C                       TRAN03 = VALINF  -  X**3 exp(-X)
C                    The recommended value is 3/EPS
C
C    XHIGH3 - DOUBLE PRECISION - The value of ln(EPSNEG). Used to prevent overflow
C                    for large x.
C
C     For values of EPS, EPSNEG, and XMIN refer to the file MACHCON.TXT
C
C     The machine-arithmetic constants are given in DATA
C     statements.
C
C
C  INTRINSIC FUNCTIONS USED:
C
C     EXP, INT, LOG, SQRT
C
C
C   OTHER MISCFUN SUBROUTINES USED:
C
C          CHEVAL , ERRPRN
C
C
C  AUTHOR:
C
C     DR. ALLAN J. MACLEOD,
C     DEPT. OF MATHEMATICS AND STATISTICS,
C     UNIVERSITY OF PAISLEY ,
C     HIGH ST.,
C     PAISLEY,
C     SCOTLAND.
C     PA1 2BE.
C
C     (e-mail: macl_ms0@paisley.ac.uk )
C
C
C  LATEST REVISION:     12 JANUARY, 1996
C
C
      INTEGER K1,K2,NTERMS,NUMEXP,NUMJN
      DOUBLE PRECISION ATRAN(0:19),CHEVAL,EIGHT,FOUR,HALF,ONE,ONEHUN,RK,
     &     RNUMJN,SUMEXP,SUM2,T,VALINF,X,XHIGH1,XHIGH2,
     &     XHIGH3,XK,XK1,XLOW1,XLOW2,XVALUE,ZERO
      CHARACTER FNNAME*6,ERRMSG*14
      DATA FNNAME/'TRAN03'/
      DATA ERRMSG/'ARGUMENT < 0.0'/
      DATA ZERO,HALF,ONE/ 0.0 D 0 , 0.5 D 0 , 1.0 D 0/
      DATA FOUR,EIGHT,ONEHUN/ 4.0 D 0 , 8.0 D 0 , 100.0 D 0 /
      DATA NUMJN,RNUMJN/ 3 , 3.0 D 0 /
      DATA VALINF/0.72123 41418 95756 57124 D 1/
      DATA ATRAN/0.76201 25432 43872 00657  D    0,
     1          -0.10567 43877 05058 53250  D    0,
     2           0.11977 80848 19657 8097   D   -1,
     3          -0.12144 01520 36983 073    D   -2,
     4           0.11550 99769 39285 47     D   -3,
     5          -0.10581 59921 24422 9      D   -4,
     6           0.94746 63385 3018         D   -6,
     7          -0.83622 12128 581          D   -7,
     8           0.73109 09927 75           D   -8,
     9          -0.63505 94778 8            D   -9,
     X           0.54911 82819              D  -10,
     1          -0.47321 3954               D  -11,
     2           0.40676 948                D  -12,
     3          -0.34897 06                 D  -13,
     4           0.29892 3                  D  -14,
     5          -0.25574                    D  -15,
     6           0.2186                     D  -16,
     7          -0.187                      D  -17,
     8           0.16                       D  -18,
     9          -0.1                        D  -19/
C
C  Machine-dependent constants
C
      DATA NTERMS/17/
      DATA XLOW1,XLOW2/2.98023224D-8,2.10953733D-154/
      DATA XHIGH1,XHIGH3/36.04365668D0,-36.73680056D0/
      DATA XHIGH2/1.35107988D16/
C
C  Start execution
C
      X = XVALUE
C
C  Error test
C
      IF ( X .LT. ZERO ) THEN
         CALL ERRPRN(FNNAME,ERRMSG)
         TRAN03 = ZERO
         RETURN
      ENDIF
C
C   Code for x < =  4.0
C
      IF ( X .LE. FOUR ) THEN
         IF ( X .LT. XLOW2 ) THEN
            TRAN03 = ZERO
         ELSE
            IF ( X .LT. XLOW1 ) THEN
               TRAN03 = ( X**(NUMJN-1) ) / ( RNUMJN - ONE )
            ELSE
               T = ( ( ( X*X ) / EIGHT ) - HALF ) - HALF
               TRAN03 = ( X**(NUMJN-1) ) * CHEVAL(NTERMS,ATRAN,T)
            ENDIF
         ENDIF
      ELSE
C 
C  Code for x > 4.0
C
         IF ( X .GT. XHIGH2 ) THEN
            SUMEXP = ONE
         ELSE
            IF ( X .LE. XHIGH1 ) THEN
               NUMEXP = INT(XHIGH1/X) + 1
               T = EXP(-X)
            ELSE
               NUMEXP = 1
               T = ONE
            ENDIF
            RK = ZERO
            DO 100 K1 = 1 , NUMEXP
               RK = RK + ONE
  100       CONTINUE
            SUMEXP = ZERO
            DO 300 K1 = 1 , NUMEXP
               SUM2 = ONE
               XK = ONE / ( RK * X )
               XK1 = ONE
               DO 200 K2 = 1 , NUMJN
                  SUM2 = SUM2 * XK1 * XK + ONE
                  XK1 = XK1 + ONE
  200          CONTINUE
               SUMEXP = SUMEXP * T + SUM2
               RK = RK - ONE
  300       CONTINUE
         ENDIF
         T = RNUMJN * LOG(X) - X + LOG(SUMEXP)
         IF ( T .LT. XHIGH3 ) THEN
            TRAN03 = VALINF
         ELSE
            TRAN03 = VALINF - EXP(T)
         ENDIF
      ENDIF
      RETURN
      END

      DOUBLE PRECISION FUNCTION TRAN04(XVALUE)

c*********************************************************************72
c
cc TRAN04 calculates the transport integral of order 4.
c
C  DESCRIPTION:
C
C      This program calculates the transport integral of order 4, defined as
C
C        TRAN04(X) = integral 0 to X { t**4 exp(t)/[exp(t)-1]**2 } dt
C
C      The program uses a Chebyshev series, the coefficients of which are
C      given to an accuracy of 20 decimal places.
C
C      This subroutine is set up to work on IEEE machines.
C      For other machines, you should retrieve the code
C      from the general MISCFUN archive.
C
C
C  ERROR RETURNS:
C
C    If XVALUE < 0.0, an error message is printed, and the program 
C    returns the value 0.0. 
C
C
C  MACHINE-DEPENDENT CONSTANTS:
C
C    NTERMS - INTEGER - The number of terms of the array ATRAN to be used.
C                       The recommended value is such that
C                             ATRAN(NTERMS) < EPS/100
C
C    XLOW2 - DOUBLE PRECISION - The value below which TRAN04 = 0.0 to machine
C                   precision. The recommended value is
C                          cube root of (3*XMIN)
C
C    XLOW1 - DOUBLE PRECISION - The value below which TRAN04 = X**3/3 to
C                   machine precision. The recommended value is
C                             sqrt(8*EPSNEG)
C
C    XHIGH1 - DOUBLE PRECISION - The value above which the exponential series for
C                    large X contains only one term. The recommended value
C                    is        - ln(EPS).
C
C    XHIGH2 - DOUBLE PRECISION - The value above which 
C                       TRAN04 = VALINF  -  X**4 exp(-X)
C                    The recommended value is 4/EPS
C
C    XHIGH3 - DOUBLE PRECISION - The value of ln(EPSNEG). Used to prevent overflow
C                    for large x.
C
C
C     For values of EPS, EPSNEG, and XMIN refer to the file MACHCON.TXT
C
C     The machine-arithmetic constants are given in DATA
C     statements.
C
C
C  INTRINSIC FUNCTIONS USED:
C
C     EXP, INT, LOG, SQRT
C
C
C   OTHER MISCFUN SUBROUTINES USED:
C
C          CHEVAL , ERRPRN
C
C
C  AUTHOR:
C
C     DR. ALLAN J. MACLEOD,
C     DEPT. OF MATHEMATICS AND STATISTICS,
C     UNIVERSITY OF PAISLEY ,
C     HIGH ST.,
C     PAISLEY,
C     SCOTLAND.
C     PA1 2BE.
C
C     (e-mail: macl_ms0@paisley.ac.uk )
C
C
C  LATEST REVISION:     12 JANUARY, 1996
C
      INTEGER K1,K2,NTERMS,NUMEXP,NUMJN
      DOUBLE PRECISION ATRAN(0:19),CHEVAL,EIGHT,FOUR,HALF,ONE,ONEHUN,RK,
     &     RNUMJN,SUMEXP,SUM2,T,VALINF,X,XHIGH1,XHIGH2,
     &     XHIGH3,XK,XK1,XLOW1,XLOW2,XVALUE,ZERO
      CHARACTER FNNAME*6,ERRMSG*14
      DATA FNNAME/'TRAN04'/
      DATA ERRMSG/'ARGUMENT < 0.0'/
      DATA ZERO,HALF,ONE/ 0.0 D 0 , 0.5 D 0 , 1.0 D 0 /
      DATA FOUR,EIGHT,ONEHUN/ 4.0 D 0 , 8.0 D 0 , 100.0 D 0 /
      DATA NUMJN,RNUMJN/ 4 , 4.0 D 0 / 
      DATA VALINF/0.25975 75760 90673 16596 D 2/
      DATA ATRAN/0.48075 70994 61511 05786  D    0,
     1          -0.81753 78810 32108 3956   D   -1,
     2           0.10027 00665 97516 2973   D   -1,
     3          -0.10599 33935 98201 507    D   -2,
     4           0.10345 06245 03040 53     D   -3,
     5          -0.96442 70548 58991        D   -5,
     6           0.87455 44408 5147         D   -6,
     7          -0.77932 12079 811          D   -7,
     8           0.68649 88614 10           D   -8,
     9          -0.59995 71076 4            D   -9,
     X           0.52136 62413              D  -10,
     1          -0.45118 3819               D  -11,
     2           0.38921 592                D  -12,
     3          -0.33493 60                 D  -13,
     4           0.28766 7                  D  -14,
     5          -0.24668                    D  -15,
     6           0.2113                     D  -16,
     7          -0.181                      D  -17,
     8           0.15                       D  -18,
     9          -0.1                        D  -19/
C
C  Machine-dependent constants
C
      DATA NTERMS/17/
      DATA XLOW1,XLOW2/2.98023224D-8,4.05653502D-103/
      DATA XHIGH1,XHIGH3/36.04365668D0,-36.73680056D0/
      DATA XHIGH2/1.80143985D16/
C
C  Start execution
C
      X = XVALUE
C
C  Error test
C
      IF ( X .LT. ZERO ) THEN
         CALL ERRPRN(FNNAME,ERRMSG)
         TRAN04 = ZERO
         RETURN
      ENDIF
C
C   Code for x < =  4.0
C
      IF ( X .LE. FOUR ) THEN
         IF ( X .LT. XLOW2 ) THEN
            TRAN04 = ZERO
         ELSE
            IF ( X .LT. XLOW1 ) THEN
               TRAN04 =  ( X ** ( NUMJN-1 ) ) / ( RNUMJN - ONE ) 
            ELSE
               T = ( ( ( X * X ) / EIGHT ) - HALF ) - HALF
               TRAN04 = ( X ** ( NUMJN-1 ) ) * CHEVAL(NTERMS,ATRAN,T) 
            ENDIF
         ENDIF
      ELSE
C 
C  Code for x > 4.0
C
         IF ( X .GT. XHIGH2 ) THEN
            SUMEXP = ONE
         ELSE
            IF ( X .LE. XHIGH1 ) THEN
               NUMEXP = INT ( XHIGH1 / X ) + 1
               T = EXP ( -X ) 
            ELSE
               NUMEXP = 1
               T = ONE
            ENDIF
            RK = ZERO
            DO 100 K1 = 1 , NUMEXP
               RK = RK + ONE
  100       CONTINUE
            SUMEXP = ZERO
            DO 300 K1 = 1 , NUMEXP
               SUM2 = ONE
               XK = ONE/ ( RK * X ) 
               XK1 = ONE
               DO 200 K2 = 1 , NUMJN
                  SUM2 = SUM2 * XK1 * XK + ONE
                  XK1 = XK1 + ONE
  200          CONTINUE
               SUMEXP = SUMEXP * T + SUM2
               RK = RK - ONE
  300       CONTINUE
         ENDIF
         T = RNUMJN * LOG( X ) - X + LOG( SUMEXP ) 
         IF ( T .LT. XHIGH3 ) THEN
            TRAN04 = VALINF
         ELSE
            TRAN04 = VALINF - EXP( T ) 
         ENDIF
      ENDIF
      RETURN
      END

      DOUBLE PRECISION FUNCTION TRAN05(XVALUE)

c*********************************************************************72
c
cc TRAN05 calculates the transport integral of order 5.
c
C  DESCRIPTION:
C
C      This program calculates the transport integral of order n, defined as
C
C        TRAN05(X) = integral 0 to X { t**5 exp(t)/[exp(t)-1]**2 } dt
C
C      The program uses a Chebyshev series, the coefficients of which are
C      given to an accuracy of 20 decimal places.
C
C      This subroutine is set up to work on IEEE machines.
C      For other machines, you should retrieve the code
C      from the general MISCFUN archive.
C
C
C  ERROR RETURNS:
C
C    If XVALUE < 0.0, an error message is printed, and the program 
C    returns the value 0.0. 
C
C
C  MACHINE-DEPENDENT CONSTANTS:
C
C    NTERMS - INTEGER - The number of terms of the array ATRAN to be used.
C                       The recommended value is such that
C                             ATRAN(NTERMS) < EPS/100
C
C    XLOW2 - DOUBLE PRECISION - The value below which TRAN05 = 0.0 to machine
C                   precision. The recommended value is
C                          4th root of (4*XMIN)
C
C    XLOW1 - DOUBLE PRECISION - The value below which TRAN05 = X**4/4 to
C                   machine precision. The recommended value is
C                             sqrt(8*EPSNEG)
C
C    XHIGH1 - DOUBLE PRECISION - The value above which the exponential series for
C                    large X contains only one term. The recommended value
C                    is        - ln(EPS).
C
C    XHIGH2 - DOUBLE PRECISION - The value above which 
C                       TRAN05 = VALINF  -  X**5 exp(-X)
C                    The recommended value is 5/EPS
C
C    XHIGH3 - DOUBLE PRECISION - The value of ln(EPSNEG). Used to prevent overflow
C                    for large x.
C
C     For values of EPS, EPSNEG, and XMIN refer to the file MACHCON.TXT
C
C     The machine-arithmetic constants are given in DATA
C     statements.
C
C
C  INTRINSIC FUNCTIONS USED:
C
C     EXP, INT, LOG, SQRT
C
C
C   OTHER MISCFUN SUBROUTINES USED:
C
C          CHEVAL , ERRPRN
C
C
C  AUTHOR:
C
C     DR. ALLAN J. MACLEOD,
C     DEPT. OF MATHEMATICS AND STATISTICS,
C     UNIVERSITY OF PAISLEY ,
C     HIGH ST.,
C     PAISLEY,
C     SCOTLAND.
C     PA1 2BE.
C
C     (e-mail: macl_ms0@paisley.ac.uk )
C
C
C  LATEST REVISION:     12 JANUARY, 1996
C
C
      INTEGER K1,K2,NTERMS,NUMEXP,NUMJN
      DOUBLE PRECISION ATRAN(0:19),CHEVAL,EIGHT,FOUR,HALF,ONE,ONEHUN,RK,
     &     RNUMJN,SUMEXP,SUM2,T,VALINF,X,XHIGH1,XHIGH2,
     &     XHIGH3,XK,XK1,XLOW1,XLOW2,XVALUE,ZERO
      CHARACTER FNNAME*6,ERRMSG*14
      DATA FNNAME/'TRAN05'/
      DATA ERRMSG/'ARGUMENT < 0.0'/
      DATA ZERO,HALF,ONE/ 0.0 D 0 , 0.5 D 0 , 1.0 D 0 /
      DATA FOUR,EIGHT,ONEHUN/ 4.0 D 0 , 8.0 D 0 , 100.0 D 0 /
      DATA NUMJN,RNUMJN/ 5 , 5.0 D 0 /
      DATA VALINF/0.12443 13306 17204 39116 D 3/
      DATA ATRAN/0.34777 77771 33910 78928  D    0,
     1          -0.66456 98897 60504 2801   D   -1,
     2           0.86110 72656 88330 882    D   -2,
     3          -0.93966 82223 75553 84     D   -3,
     4           0.93632 48060 81513 4      D   -4,
     5          -0.88571 31934 08328        D   -5,
     6           0.81191 49891 4503         D   -6,
     7          -0.72957 65423 277          D   -7,
     8           0.64697 14550 45           D   -8,
     9          -0.56849 02825 5            D   -9,
     X           0.49625 59787              D  -10,
     1          -0.43109 3996               D  -11,
     2           0.37310 094                D  -12,
     3          -0.32197 69                 D  -13,
     4           0.27722 0                  D  -14,
     5          -0.23824                    D  -15,
     6           0.2044                     D  -16,
     7          -0.175                      D  -17,
     8           0.15                       D  -18,
     9          -0.1                        D  -19/
C
C  Machine-dependent constants
C
      DATA NTERMS/17/
      DATA XLOW1,XLOW2/2.98023224D-8,1.72723372D-77/
      DATA XHIGH1,XHIGH3/36.04365668D0,-36.73680056D0/
      DATA XHIGH2/2.25179981D16/
C
C  Start execution
C
      X = XVALUE
C
C  Error test
C
      IF ( X .LT. ZERO ) THEN
         CALL ERRPRN(FNNAME,ERRMSG)
         TRAN05 = ZERO
         RETURN
      ENDIF
C
C   Code for x < =  4.0
C
      IF ( X .LE. FOUR ) THEN
         IF ( X .LT. XLOW2 ) THEN
            TRAN05 = ZERO
         ELSE
            IF ( X .LT. XLOW1 ) THEN
               TRAN05 =  ( X ** ( NUMJN - 1 ) ) / ( RNUMJN - ONE ) 
            ELSE
               T = ( ( ( X * X ) / EIGHT ) - HALF ) - HALF
               TRAN05 = ( X ** ( NUMJN-1 ) ) * CHEVAL(NTERMS,ATRAN,T) 
            ENDIF
         ENDIF
      ELSE
C 
C  Code for x > 4.0
C
         IF ( X .GT. XHIGH2 ) THEN
            SUMEXP = ONE
         ELSE
            IF ( X .LE. XHIGH1 ) THEN
               NUMEXP = INT ( XHIGH1 / X )  + 1
               T = EXP ( -X ) 
            ELSE
               NUMEXP = 1
               T = ONE
            ENDIF
            RK = ZERO
            DO 100 K1 = 1 , NUMEXP
               RK = RK + ONE
  100       CONTINUE
            SUMEXP = ZERO
            DO 300 K1 = 1 , NUMEXP
               SUM2 = ONE
               XK = ONE / ( RK * X ) 
               XK1 = ONE
               DO 200 K2 = 1 , NUMJN
                  SUM2 = SUM2 * XK1 * XK + ONE
                  XK1 = XK1 + ONE
  200          CONTINUE
               SUMEXP = SUMEXP * T + SUM2
               RK = RK - ONE
  300       CONTINUE
         ENDIF
         T = RNUMJN * LOG ( X ) - X + LOG( SUMEXP ) 
         IF ( T .LT. XHIGH3 ) THEN
            TRAN05 = VALINF
         ELSE
            TRAN05 = VALINF - EXP( T ) 
         ENDIF
      ENDIF
      RETURN
      END

      DOUBLE PRECISION FUNCTION TRAN06(XVALUE)

c*********************************************************************72
c
cc TRAN06 calculates the transport integral of order 6.
c
C  DESCRIPTION:
C
C      This program calculates the transport integral of order 6, defined as
C
C        TRAN06(X) = integral 0 to X { t**6 exp(t)/[exp(t)-1]**2 } dt
C
C      The program uses a Chebyshev series, the coefficients of which are
C      given to an accuracy of 20 decimal places.
C
C      This subroutine is set up to work on IEEE machines.
C      For other machines, you should retrieve the code
C      from the general MISCFUN archive.
C
C
C  ERROR RETURNS:
C
C    If XVALUE < 0.0, an error message is printed, and the program 
C    returns the value 0.0. 
C
C
C  MACHINE-DEPENDENT CONSTANTS:
C
C    NTERMS - INTEGER - The number of terms of the array ATRAN to be used.
C                       The recommended value is such that
C                             ATRAN(NTERMS) < EPS/100
C
C    XLOW2 - DOUBLE PRECISION - The value below which TRAN06 = 0.0 to machine
C                   precision. The recommended value is
C                          5th root of (5*XMIN)
C
C    XLOW1 - DOUBLE PRECISION - The value below which TRAN06 = X**5/5 to
C                   machine precision. The recommended value is
C                             sqrt(8*EPSNEG)
C
C    XHIGH1 - DOUBLE PRECISION - The value above which the exponential series for
C                    large X contains only one term. The recommended value
C                    is        - ln(EPS).
C
C    XHIGH2 - DOUBLE PRECISION - The value above which 
C                       TRAN06 = VALINF  -  X**6 exp(-X)
C                    The recommended value is 6/EPS
C
C    XHIGH3 - DOUBLE PRECISION - The value of ln(EPSNEG). Used to prevent overflow
C                    for large x.
C
C     For values of EPS, EPSNEG, and XMIN refer to the file MACHCON.TXT
C
C     The machine-arithmetic constants are given in DATA
C     statements.
C
C
C  INTRINSIC FUNCTIONS USED:
C
C     EXP, INT, LOG, SQRT
C
C
C   OTHER MISCFUN SUBROUTINES USED:
C
C          CHEVAL , ERRPRN
C
C
C  AUTHOR:
C
C     DR. ALLAN J. MACLEOD,
C     DEPT. OF MATHEMATICS AND STATISTICS,
C     UNIVERSITY OF PAISLEY ,
C     HIGH ST.,
C     PAISLEY,
C     SCOTLAND.
C     PA1 2BE.
C
C     (e-mail: macl_ms0@paisley.ac.uk )
C
C
C  LATEST REVISION:     12 JANUARY, 1996
C
      INTEGER K1,K2,NTERMS,NUMEXP,NUMJN
      DOUBLE PRECISION ATRAN(0:19),CHEVAL,EIGHT,FOUR,HALF,ONE,ONEHUN,RK,
     &     RNUMJN,SUMEXP,SUM2,T,VALINF,X,XHIGH1,XHIGH2,
     &     XHIGH3,XK,XK1,XLOW1,XLOW2,XVALUE,ZERO
      CHARACTER FNNAME*6,ERRMSG*14
      DATA FNNAME/'TRAN06'/
      DATA ERRMSG/'ARGUMENT < 0.0'/
      DATA ZERO,HALF,ONE/ 0.0 D 0 , 0.5 D 0 , 1.0 D 0 /
      DATA FOUR,EIGHT,ONEHUN/ 4.0 D 0 , 8.0 D 0 , 100.0 D 0 /
      DATA NUMJN,RNUMJN/ 6 , 6.0 D 0 /
      DATA VALINF/0.73248 70046 28803 38059 D 3/
      DATA ATRAN/0.27127 33539 78400 08227  D    0,
     1          -0.55886 10553 19145 3393   D   -1,
     2           0.75391 95132 90083 056    D   -2,
     3          -0.84351 13857 92112 19     D   -3,
     4           0.85490 98079 67670 2      D   -4,
     5          -0.81871 54932 93098        D   -5,
     6           0.75754 24042 7986         D   -6,
     7          -0.68573 06541 831          D   -7,
     8           0.61170 03760 31           D   -8,
     9          -0.54012 70702 4            D   -9,
     X           0.47343 06435              D  -10,
     1          -0.41270 1055               D  -11,
     2           0.35825 603                D  -12,
     3          -0.30997 52                 D  -13,
     4           0.26750 1                  D  -14,
     5          -0.23036                    D  -15,
     6           0.1980                     D  -16,
     7          -0.170                      D  -17,
     8           0.15                       D  -18,
     9          -0.1                        D  -19/
C
C  Machine-dependent constants
C
      DATA NTERMS/17/
      DATA XLOW1,XLOW2/2.98023224D-8,4.06689432D-62/
      DATA XHIGH1,XHIGH3/36.04365668D0,-36.73680056D0/
      DATA XHIGH2/2.70215977D16/
C
C  Start execution
C
      X = XVALUE
C
C  Error test
C
      IF ( X .LT. ZERO ) THEN
         CALL ERRPRN(FNNAME,ERRMSG)
         TRAN06 = ZERO
         RETURN
      ENDIF
C
C   Code for x < =  4 .0
C
      IF ( X .LE. FOUR ) THEN
         IF ( X .LT. XLOW2 ) THEN
            TRAN06 = ZERO
         ELSE
            IF ( X .LT. XLOW1 ) THEN
               TRAN06 =  ( X ** ( NUMJN-1 ) ) / ( RNUMJN - ONE ) 
            ELSE
               T =  ( ( ( X * X ) / EIGHT ) - HALF ) - HALF
               TRAN06 = ( X ** ( NUMJN-1 )  ) * CHEVAL(NTERMS,ATRAN,T) 
            ENDIF
         ENDIF
      ELSE
C 
C  Code for x > 4 .0
C
         IF ( X .GT. XHIGH2 ) THEN
            SUMEXP = ONE
         ELSE
            IF ( X .LE. XHIGH1 ) THEN
               NUMEXP = INT ( XHIGH1 / X ) + 1
               T = EXP( - X ) 
            ELSE
               NUMEXP = 1
               T = ONE
            ENDIF
            RK = ZERO
            DO 100 K1 = 1 , NUMEXP
               RK = RK + ONE
  100       CONTINUE
            SUMEXP = ZERO
            DO 300 K1 = 1 , NUMEXP
               SUM2 = ONE
               XK = ONE / ( RK * X ) 
               XK1 = ONE
               DO 200 K2 = 1 , NUMJN
                  SUM2 = SUM2 * XK1 * XK + ONE
                  XK1 = XK1 + ONE
  200          CONTINUE
               SUMEXP = SUMEXP * T + SUM2
               RK = RK - ONE
  300       CONTINUE
         ENDIF
         T = RNUMJN * LOG( X ) - X + LOG( SUMEXP ) 
         IF ( T .LT. XHIGH3 ) THEN
            TRAN06 = VALINF
         ELSE
            TRAN06 = VALINF - EXP( T ) 
         ENDIF
      ENDIF
      RETURN
      END

      DOUBLE PRECISION FUNCTION TRAN07(XVALUE)

c*********************************************************************72
c
cc TRAN07 calculates the transport integral of order 7.
c
C  DESCRIPTION:
C
C      This program calculates the transport integral of order 7, defined as
C
C        TRAN07(X) = integral 0 to X { t**7 exp(t)/[exp(t)-1]**2 } dt
C
C      The program uses a Chebyshev series, the coefficients of which are
C      given to an accuracy of 20 decimal places.
C
C      This subroutine is set up to work on IEEE machines.
C      For other machines, you should retrieve the code
C      from the general MISCFUN archive.
C
C
C  ERROR RETURNS:
C
C    If XVALUE < 0.0, an error message is printed, and the program 
C    returns the value 0.0. 
C
C
C  MACHINE-DEPENDENT CONSTANTS:
C
C    NTERMS - INTEGER - The number of terms of the array ATRAN to be used.
C                       The recommended value is such that
C                             ATRAN(NTERMS) < EPS/100
C
C    XLOW2 - DOUBLE PRECISION - The value below which TRAN07 = 0.0 to machine
C                   precision. The recommended value is
C                          6th root of (6*XMIN)
C
C    XLOW1 - DOUBLE PRECISION - The value below which TRAN07 = X**6/6 to
C                   machine precision. The recommended value is
C                             sqrt(8*EPSNEG)
C
C    XHIGH1 - DOUBLE PRECISION - The value above which the exponential series for
C                    large X contains only one term. The recommended value
C                    is        - ln(EPS).
C
C    XHIGH2 - DOUBLE PRECISION - The value above which 
C                       TRAN07 = VALINF  -  X**7 exp(-X)
C                    The recommended value is 7/EPS
C
C    XHIGH3 - DOUBLE PRECISION - The value of ln(EPSNEG). Used to prevent overflow
C                    for large x.
C
C     For values of EPS, EPSNEG, and XMIN refer to the file MACHCON.TXT
C
C     The machine-arithmetic constants are given in DATA
C     statements.
C
C
C  INTRINSIC FUNCTIONS USED:
C
C     EXP, INT, LOG, SQRT
C
C
C   OTHER MISCFUN SUBROUTINES USED:
C
C          CHEVAL , ERRPRN
C
C
C  AUTHOR:
C
C     DR. ALLAN J. MACLEOD,
C     DEPT. OF MATHEMATICS AND STATISTICS,
C     UNIVERSITY OF PAISLEY ,
C     HIGH ST.,
C     PAISLEY,
C     SCOTLAND.
C     PA1 2BE.
C
C     (e-mail: macl_ms0@paisley.ac.uk )
C
C
C  LATEST REVISION:     12 JANUARY, 1996
C
      INTEGER K1,K2,NTERMS,NUMEXP,NUMJN
      DOUBLE PRECISION ATRAN(0:19),CHEVAL,EIGHT,FOUR,HALF,ONE,ONEHUN,RK,
     &     RNUMJN,SUMEXP,SUM2,T,VALINF,X,XHIGH1,XHIGH2,
     &     XHIGH3,XK,XK1,XLOW1,XLOW2,XVALUE,ZERO
      CHARACTER FNNAME*6,ERRMSG*14
      DATA FNNAME/'TRAN07'/
      DATA ERRMSG/'ARGUMENT < 0.0'/
      DATA ZERO,HALF,ONE/ 0.0 D 0 , 0.5 D 0 , 1.0 D 0/
      DATA FOUR,EIGHT,ONEHUN/ 4.0 D 0 , 8.0 D 0 , 100.0 D 0 /
      DATA NUMJN,RNUMJN/ 7 , 7.0 D 0/
      DATA VALINF/0.50820 80358 00489 10473 D 4/
      DATA ATRAN/0.22189 25073 40104 04423  D    0,
     1          -0.48167 51061 17799 3694   D   -1,
     2           0.67009 24481 03153 629    D   -2,
     3          -0.76495 18344 30825 57     D   -3,
     4           0.78634 85592 34869 0      D   -4,
     5          -0.76102 51808 87504        D   -5,
     6           0.70991 69629 9917         D   -6,
     7          -0.64680 25624 903          D   -7,
     8           0.58003 92339 60           D   -8,
     9          -0.51443 37014 9            D   -9,
     X           0.45259 44183              D  -10,
     1          -0.39580 0363               D  -11,
     2           0.34453 785                D  -12,
     3          -0.29882 92                 D  -13,
     4           0.25843 4                  D  -14,
     5          -0.22297                    D  -15,
     6           0.1920                     D  -16,
     7          -0.165                      D  -17,
     8           0.14                       D  -18,
     9          -0.1                        D  -19/
C
C  Machine-dependent constants
C
      DATA NTERMS/17/
      DATA XLOW1,XLOW2/2.98023224D-8,7.14906557D-52/
      DATA XHIGH1,XHIGH3/36.04365668D0,-36.73680056D0/
      DATA XHIGH2/3.15251973D16/
C
C  Start execution
C
      X = XVALUE
C
C  Error test
C
      IF ( X .LT. ZERO ) THEN
         CALL ERRPRN(FNNAME,ERRMSG)
         TRAN07 = ZERO
         RETURN
      ENDIF
C
C   Code for x <= 4.0
C
      IF ( X .LE. FOUR ) THEN
         IF ( X .LT. XLOW2 ) THEN
            TRAN07 = ZERO
         ELSE
            IF ( X .LT. XLOW1 ) THEN
               TRAN07 = ( X**(NUMJN-1) ) / ( RNUMJN - ONE )
            ELSE
               T = ( ( ( X*X ) / EIGHT ) - HALF ) - HALF
               TRAN07 = ( X**(NUMJN-1) ) * CHEVAL(NTERMS,ATRAN,T)
            ENDIF
         ENDIF
      ELSE
C 
C  Code for x > 4.0
C
         IF ( X .GT. XHIGH2 ) THEN
            SUMEXP = ONE
         ELSE
            IF ( X .LE. XHIGH1 ) THEN
               NUMEXP = INT ( XHIGH1/X ) + 1
               T = EXP( -X )
            ELSE
               NUMEXP = 1
               T = ONE
            ENDIF
            RK = ZERO
            DO 100 K1 = 1 , NUMEXP
               RK = RK + ONE
  100       CONTINUE
            SUMEXP = ZERO
            DO 300 K1 = 1 , NUMEXP
               SUM2 = ONE
               XK = ONE / ( RK * X )
               XK1 = ONE
               DO 200 K2 = 1 , NUMJN
                  SUM2 = SUM2 * XK1 * XK + ONE
                  XK1 = XK1 + ONE
  200          CONTINUE
               SUMEXP = SUMEXP * T + SUM2
               RK = RK - ONE
  300       CONTINUE
         ENDIF
         T = RNUMJN * LOG(X) - X + LOG(SUMEXP)
         IF ( T .LT. XHIGH3 ) THEN
            TRAN07 = VALINF
         ELSE
            TRAN07 = VALINF - EXP(T)
         ENDIF
      ENDIF
      RETURN
      END

      DOUBLE PRECISION FUNCTION TRAN08(XVALUE)

c*********************************************************************72
c
cc TRAN08 calculates the transport integral of order 8.
c
C  DESCRIPTION:
C
C      This program calculates the transport integral of order 8, defined as
C
C        TRAN08(X) = integral 0 to X { t**8 exp(t)/[exp(t)-1]**2 } dt
C
C      The program uses a Chebyshev series, the coefficients of which are
C      given to an accuracy of 20 decimal places.
C
C      This subroutine is set up to work on IEEE machines.
C      For other machines, you should retrieve the code
C      from the general MISCFUN archive.
C
C
C  ERROR RETURNS:
C
C    If XVALUE < 0.0, an error message is printed, and the program 
C    returns the value 0.0. 
C
C
C  MACHINE-DEPENDENT CONSTANTS:
C
C    NTERMS - INTEGER - The number of terms of the array ATRAN to be used.
C                       The recommended value is such that
C                             ATRAN(NTERMS) < EPS/100
C
C    XLOW2 - DOUBLE PRECISION - The value below which TRAN08 = 0.0 to machine
C                   precision. The recommended value is
C                          7th root of (7*XMIN)
C
C    XLOW1 - DOUBLE PRECISION - The value below which TRAN08 = X**7/7 to
C                   machine precision. The recommended value is
C                             sqrt(8*EPSNEG)
C
C    XHIGH1 - DOUBLE PRECISION - The value above which the exponential series for
C                    large X contains only one term. The recommended value
C                    is        - ln(EPS).
C
C    XHIGH2 - DOUBLE PRECISION - The value above which 
C                       TRAN08 = VALINF  -  X**8 exp(-X)
C                    The recommended value is 8/EPS
C
C    XHIGH3 - DOUBLE PRECISION - The value of ln(EPSNEG). Used to prevent overflow
C                    for large x.
C
C     For values of EPS, EPSNEG, and XMIN refer to the file MACHCON.TXT
C
C     The machine-arithmetic constants are given in DATA
C     statements.
C
C
C  INTRINSIC FUNCTIONS USED:
C
C     EXP, INT, LOG, SQRT
C
C
C   OTHER MISCFUN SUBROUTINES USED:
C
C          CHEVAL , ERRPRN
C
C
C  AUTHOR:
C
C     DR. ALLAN J. MACLEOD,
C     DEPT. OF MATHEMATICS AND STATISTICS,
C     UNIVERSITY OF PAISLEY ,
C     HIGH ST.,
C     PAISLEY,
C     SCOTLAND.
C     PA1 2BE.
C
C     (e-mail: macl_ms0@paisley.ac.uk )
C
C
C  LATEST REVISION:     12 JANUARY, 1996
C
      INTEGER K1,K2,NTERMS,NUMEXP,NUMJN
      DOUBLE PRECISION ATRAN(0:19),CHEVAL,EIGHT,FOUR,HALF,ONE,ONEHUN,RK,
     &     RNUMJN,SUMEXP,SUM2,T,VALINF,X,XHIGH1,XHIGH2,
     &     XHIGH3,XK,XK1,XLOW1,XLOW2,XVALUE,ZERO
      CHARACTER FNNAME*6,ERRMSG*14
      DATA FNNAME/'TRAN08'/
      DATA ERRMSG/'ARGUMENT < 0.0'/
      DATA ZERO,HALF,ONE/ 0.0 D 0 , 0.5 D 0 , 1.0 D 0 /
      DATA FOUR,EIGHT,ONEHUN/ 4.0 D 0 , 8.0 D 0 , 100.0 D 0 /
      DATA NUMJN,RNUMJN/ 8 , 8.0 D 0 /
      DATA VALINF/0.40484 39900 19011 15764 D 5/
      DATA ATRAN/0.18750 69577 40437 19233  D    0,
     1          -0.42295 27646 09367 3337   D   -1,
     2           0.60281 48569 29065 592    D   -2,
     3          -0.69961 05481 18147 76     D   -3,
     4           0.72784 82421 29878 9      D   -4,
     5          -0.71084 62500 50067        D   -5,
     6           0.66786 70689 0115         D   -6,
     7          -0.61201 57501 844          D   -7,
     8           0.55146 52644 74           D   -8,
     9          -0.49105 30705 2            D   -9,
     X           0.43350 00869              D  -10,
     1          -0.38021 8700               D  -11,
     2           0.33182 369                D  -12,
     3          -0.28845 12                 D  -13,
     4           0.24995 8                  D  -14,
     5          -0.21605                    D  -15,
     6           0.1863                     D  -16,
     7          -0.160                      D  -17,
     8           0.14                       D  -18,
     9          -0.1                        D  -19/
C
C  Machine-dependent constants
C
      DATA NTERMS/17/
      DATA XLOW1,XLOW2/2.98023224D-8,1.48029723D-44/
      DATA XHIGH1,XHIGH3/36.04365668D0,-36.73680056D0/
      DATA XHIGH2/3.6028797D16/
C
C  Start execution
C
      X = XVALUE
C
C  Error test
C
      IF ( X .LT. ZERO ) THEN
         CALL ERRPRN(FNNAME,ERRMSG)
         TRAN08 = ZERO
         RETURN
      ENDIF
C
C   Code for x < =  4.0
C
      IF ( X .LE. FOUR ) THEN
         IF ( X .LT. XLOW2 ) THEN
            TRAN08 = ZERO
         ELSE
            IF ( X .LT. XLOW1 ) THEN
               TRAN08 = ( X ** ( NUMJN - 1 ) ) / ( RNUMJN - ONE ) 
            ELSE
               T = ( ( ( X * X ) / EIGHT ) - HALF )  - HALF
               TRAN08 = ( X ** ( NUMJN - 1 ) ) * CHEVAL(NTERMS,ATRAN,T) 
            ENDIF
         ENDIF
      ELSE
C 
C  Code for x > 4.0
C
         IF ( X .GT. XHIGH2 ) THEN
            SUMEXP = ONE
         ELSE
            IF ( X .LE. XHIGH1 ) THEN
               NUMEXP = INT ( XHIGH1 / X ) + 1
               T = EXP ( - X ) 
            ELSE
               NUMEXP = 1
               T = ONE
            ENDIF
            RK = ZERO
            DO 100 K1 = 1 , NUMEXP
               RK = RK + ONE
  100       CONTINUE
            SUMEXP = ZERO
            DO 300 K1 = 1 , NUMEXP
               SUM2 = ONE
               XK = ONE / ( RK * X ) 
               XK1 = ONE
               DO 200 K2 = 1 , NUMJN
                  SUM2 = SUM2 * XK1 * XK + ONE
                  XK1 = XK1 + ONE
  200          CONTINUE
               SUMEXP = SUMEXP * T + SUM2
               RK = RK - ONE
  300       CONTINUE
         ENDIF
         T = RNUMJN * LOG( X ) - X + LOG( SUMEXP ) 
         IF ( T .LT. XHIGH3 ) THEN
            TRAN08 = VALINF
         ELSE
            TRAN08 = VALINF - EXP( T ) 
         ENDIF
      ENDIF
      RETURN
      END

      DOUBLE PRECISION FUNCTION TRAN09(XVALUE)

c*********************************************************************72
c
cc TRAN09 calculates the transport integral of order 9.
c
C  DESCRIPTION:
C
C      This program calculates the transport integral of order 9, defined as
C
C        TRAN09(X) = integral 0 to X { t**9 exp(t)/[exp(t)-1]**2 } dt
C
C      The program uses a Chebyshev series, the coefficients of which are
C      given to an accuracy of 20 decimal places.
C
C      This subroutine is set up to work on IEEE machines.
C      For other machines, you should retrieve the code
C      from the general MISCFUN archive.
C
C
C  ERROR RETURNS:
C
C    If XVALUE < 0.0, an error message is printed, and the program 
C    returns the value 0.0. 
C
C
C  MACHINE-DEPENDENT CONSTANTS:
C
C    NTERMS - INTEGER - The number of terms of the array ATRAN to be used.
C                       The recommended value is such that
C                             ATRAN(NTERMS) < EPS/100
C
C    XLOW2 - DOUBLE PRECISION - The value below which TRAN09 = 0.0 to machine
C                   precision. The recommended value is
C                          8th root of (8*XMIN)
C
C    XLOW1 - DOUBLE PRECISION - The value below which TRAN09 = X**8/8 to
C                   machine precision. The recommended value is
C                             sqrt(8*EPSNEG)
C
C    XHIGH1 - DOUBLE PRECISION - The value above which the exponential series for
C                    large X contains only one term. The recommended value
C                    is        - ln(EPS).
C
C    XHIGH2 - DOUBLE PRECISION - The value above which 
C                       TRAN09 = VALINF  -  X**9 exp(-X)
C                    The recommended value is 9/EPS
C
C    XHIGH3 - DOUBLE PRECISION - The value of ln(EPSNEG). Used to prevent overflow
C                    for large x.
C
C     For values of EPS, EPSNEG, and XMIN refer to the file MACHCON.TXT
C
C     The machine-arithmetic constants are given in DATA
C     statements.
C
C
C  INTRINSIC FUNCTIONS USED:
C
C     EXP, INT, LOG, SQRT
C
C
C   OTHER MISCFUN SUBROUTINES USED:
C
C          CHEVAL , ERRPRN
C
C
C  AUTHOR:
C
C     DR. ALLAN J. MACLEOD,
C     DEPT. OF MATHEMATICS AND STATISTICS,
C     UNIVERSITY OF PAISLEY ,
C     HIGH ST.,
C     PAISLEY,
C     SCOTLAND.
C     PA1 2BE.
C
C     (e-mail: macl_ms0@paisley.ac.uk )
C
C
C  LATEST REVISION:     12 JANUARY, 1996
C
      INTEGER K1,K2,NTERMS,NUMEXP,NUMJN
      DOUBLE PRECISION ATRAN(0:19),CHEVAL,EIGHT,FOUR,HALF,ONE,ONEHUN,RK,
     &     RNUMJN,SUMEXP,SUM2,T,VALINF,X,XHIGH1,XHIGH2,
     &     XHIGH3,XK,XK1,XLOW1,XLOW2,XVALUE,ZERO
      CHARACTER FNNAME*6,ERRMSG*14
      DATA FNNAME/'TRAN09'/
      DATA ERRMSG/'ARGUMENT < 0.0'/
      DATA ZERO,HALF,ONE/ 0.0 D 0 , 0.5 D 0 , 1.0 D 0 /
      DATA FOUR,EIGHT,ONEHUN/ 4.0 D 0 , 8.0 D 0 , 100.0 D 0 /
      DATA NUMJN,RNUMJN/ 9 , 9.0 D 0 /
      DATA VALINF/0.36360 88055 88728 71397 D 6/
      DATA ATRAN/0.16224 04999 19498 46835  D    0,
     1          -0.37683 51452 19593 7773   D   -1,
     2           0.54766 97159 17719 770    D   -2,
     3          -0.64443 94500 94495 21     D   -3,
     4           0.67736 45285 28098 3      D   -4,
     5          -0.66681 34975 82042        D   -5,
     6           0.63047 56001 9047         D   -6,
     7          -0.58074 78663 611          D   -7,
     8           0.52555 13051 23           D   -8,
     9          -0.46968 86176 1            D   -9,
     X           0.41593 95065              D  -10,
     1          -0.36580 8491               D  -11,
     2           0.32000 794                D  -12,
     3          -0.27876 51                 D  -13,
     4           0.24201 7                  D  -14,
     5          -0.20953                    D  -15,
     6           0.1810                     D  -16,
     7          -0.156                      D  -17,
     8           0.13                       D  -18,
     9          -0.1                        D  -19/
C
C  Machine-dependent constants (for IEEE machines) 
C
      DATA NTERMS/17/
      DATA XLOW1,XLOW2/2.98023224D-8,4.5321503D-39/
      DATA XHIGH1,XHIGH3/36.04365668D0,-36.73680056D0/
      DATA XHIGH2/4.05323966D16/
C
C  Start execution
C
      X = XVALUE
C
C  Error test
C
      IF ( X .LT. ZERO ) THEN
         CALL ERRPRN(FNNAME,ERRMSG)
         TRAN09 = ZERO
         RETURN
      ENDIF
C
C   Code for x < =  4.0
C
      IF ( X .LE. FOUR ) THEN
         IF ( X .LT. XLOW2 ) THEN
            TRAN09 = ZERO
         ELSE
            IF ( X .LT. XLOW1 ) THEN
               TRAN09 = ( X ** ( NUMJN - 1 ) ) / ( RNUMJN - ONE ) 
            ELSE
               T = ( ( ( X * X ) / EIGHT ) - HALF ) - HALF
               TRAN09 = ( X ** ( NUMJN-1 ) ) * CHEVAL(NTERMS,ATRAN,T) 
            ENDIF
         ENDIF
      ELSE
C 
C  Code for x > 4.0
C
         IF ( X .GT. XHIGH2 ) THEN
            SUMEXP = ONE
         ELSE
            IF ( X .LE. XHIGH1 ) THEN
               NUMEXP = INT ( XHIGH1 / X ) + 1
               T = EXP( -X ) 
            ELSE
               NUMEXP = 1
               T = ONE
            ENDIF
            RK = ZERO
            DO 100 K1 = 1 , NUMEXP
               RK = RK + ONE
  100       CONTINUE
            SUMEXP = ZERO
            DO 300 K1 = 1 , NUMEXP
               SUM2 = ONE
               XK = ONE / ( RK * X ) 
               XK1 = ONE
               DO 200 K2 = 1 , NUMJN
                  SUM2 = SUM2 * XK1 * XK + ONE
                  XK1 = XK1 + ONE
  200          CONTINUE
               SUMEXP = SUMEXP * T + SUM2
               RK = RK - ONE
  300       CONTINUE
         ENDIF
         T = RNUMJN * LOG( X ) - X + LOG( SUMEXP ) 
         IF ( T.LT.XHIGH3 ) THEN
            TRAN09 = VALINF
         ELSE
            TRAN09 = VALINF - EXP( T ) 
         ENDIF
      ENDIF
      RETURN
      END

      DOUBLE PRECISION FUNCTION Y0INT(XVALUE)

c*********************************************************************72
c
cc BESSEL_Y0_INT calculates the integral of the Bessel function Y0.
c
C   DESCRIPTION:
C
C      This function calculates the integral of the Bessel
C      function Y0, defined as
C
C        Y0INT(x) = {integral 0 to x} Y0(t) dt
C
C      The code uses Chebyshev expansions whose coefficients are
C      given to 20 decimal places.
C
C      This subroutine is set up to work on IEEE machines.
C      For other machines, you should retrieve the code
C      from the general MISCFUN archive.
C
C
C   ERROR RETURNS:
C
C      If x < 0.0, the function is undefined. An error message
C      is printed and the function returns the value 0.0.
C
C      If the value of x is too large, it is impossible to 
C      accurately compute the trigonometric functions used. An
C      error message is printed, and the function returns the
C      value 1.0.
C
C
C   MACHINE-DEPENDENT CONSTANTS:
C
C      NTERM1 - The no. of terms to be used from the array
C                ARJ01. The recommended value is such that
C                   ABS(ARJ01(NTERM1)) < EPS/100
C
C      NTERM2 - The no. of terms to be used from the array
C                ARY01. The recommended value is such that
C                   ABS(ARY01(NTERM2)) < EPS/100
C
C      NTERM3 - The no. of terms to be used from the array
C                ARY0A1. The recommended value is such that
C                   ABS(ARY0A1(NTERM3)) < EPS/100
C
C      NTERM4 - The no. of terms to be used from the array
C                ARY0A2. The recommended value is such that
C                   ABS(ARY0A2(NTERM4)) < EPS/100
C
C      XLOW - The value of x below which
C                  Y0INT(x) = x*(ln(x) - 0.11593)*2/pi
C             to machine-precision. The recommended value is
C                 sqrt(9*EPSNEG)
C
C      XHIGH - The value of x above which it is impossible
C              to calculate (x-pi/4) accurately. The recommended
C              value is      1/EPSNEG
C
C      For values of EPS and EPSNEG, refer to the file MACHCON.TXT
C
C      The machine-arithmetic constants are given in DATA
C      statements.
C
C
C   INTRINSIC FUNCTIONS USED:
C
C      COS , LOG , SIN , SQRT
C
C
C   OTHER MISCFUN SUBROUTINES USED:
C
C          CHEVAL , ERRPRN
C
C
C   AUTHOR:
C          Dr. Allan J. MacLeod,
C          Dept. of Mathematics and Statistics,
C          University of Paisley,
C          Paisley,
C          SCOTLAND
C          PA1 2BE
C
C          (e-mail: macl_ms0@paisley.ac.uk)
C
C
C   LATEST REVISION:
C                   12 JANUARY, 1996
C
      INTEGER NTERM1,NTERM2,NTERM3,NTERM4
      DOUBLE PRECISION ARJ01(0:23),ARY01(0:24),ARY0A1(0:21),
     1     ARY0A2(0:18),CHEVAL,FIVE12,GAL2M1,GAMLN2,
     2     NINE,ONE,ONEHUN,ONE28,PIB41,PIB411,PIB412,
     3     PIB42,RT2BPI,SIXTEN,T,TEMP,TWOBPI,X,XHIGH,
     4     XLOW,XMPI4,XVALUE,ZERO
      CHARACTER FNNAME*6,ERMSG1*14,ERMSG2*18
      DATA FNNAME/'Y0INT '/
      DATA ERMSG1/'ARGUMENT < 0.0'/
      DATA ERMSG2/'ARGUMENT TOO LARGE'/
      DATA ZERO,ONE/ 0.0 D 0 , 1.0 D 0 /
      DATA NINE,SIXTEN/ 9.0 D 0 , 16.0 D 0 /
      DATA ONEHUN,ONE28,FIVE12/ 100.0 D 0 , 128.0 D 0 , 512.0 D 0 /
      DATA RT2BPI/0.79788 45608 02865 35588 D 0/
      DATA PIB411,PIB412/ 201.0 D 0 , 256.0 D 0/
      DATA PIB42/0.24191 33974 48309 61566 D -3/
      DATA TWOBPI/0.63661 97723 67581 34308 D 0/
      DATA GAL2M1/-1.11593 15156 58412 44881 D 0/      
      DATA GAMLN2/-0.11593 15156 58412 44881 D 0/      
      DATA ARJ01(0)/  0.38179 27932 16901 73518  D    0/
      DATA ARJ01(1)/ -0.21275 63635 05053 21870  D    0/
      DATA ARJ01(2)/  0.16754 21340 72157 94187  D    0/
      DATA ARJ01(3)/ -0.12853 20977 21963 98954  D    0/
      DATA ARJ01(4)/  0.10114 40545 57788 47013  D    0/
      DATA ARJ01(5)/ -0.91007 95343 20156 8859   D   -1/
      DATA ARJ01(6)/  0.64013 45264 65687 3103   D   -1/
      DATA ARJ01(7)/ -0.30669 63029 92675 4312   D   -1/
      DATA ARJ01(8)/  0.10308 36525 32506 4201   D   -1/
      DATA ARJ01(9)/ -0.25567 06503 99956 918    D   -2/
      DATA ARJ01(10)/ 0.48832 75580 57983 04     D   -3/
      DATA ARJ01(11)/-0.74249 35126 03607 7      D   -4/
      DATA ARJ01(12)/ 0.92226 05637 30861        D   -5/
      DATA ARJ01(13)/-0.95522 82830 7083         D   -6/
      DATA ARJ01(14)/ 0.83883 55845 986          D   -7/
      DATA ARJ01(15)/-0.63318 44888 58           D   -8/
      DATA ARJ01(16)/ 0.41560 50422 1            D   -9/
      DATA ARJ01(17)/-0.23955 29307              D  -10/
      DATA ARJ01(18)/ 0.12228 6885               D  -11/
      DATA ARJ01(19)/-0.55697 11                 D  -13/
      DATA ARJ01(20)/ 0.22782 0                  D  -14/
      DATA ARJ01(21)/-0.8417                     D  -16/
      DATA ARJ01(22)/ 0.282                      D  -17/
      DATA ARJ01(23)/-0.9                        D  -19/
      DATA ARY01(0)/  0.54492 69630 27243 65490  D    0/
      DATA ARY01(1)/ -0.14957 32358 86847 82157  D    0/
      DATA ARY01(2)/  0.11085 63448 62548 42337  D    0/
      DATA ARY01(3)/ -0.94953 30018 68377 7109   D   -1/
      DATA ARY01(4)/  0.68208 17786 99145 6963   D   -1/
      DATA ARY01(5)/ -0.10324 65338 33682 00408  D    0/
      DATA ARY01(6)/  0.10625 70328 75344 25491  D    0/
      DATA ARY01(7)/ -0.62583 67679 96168 1990   D   -1/
      DATA ARY01(8)/  0.23856 45760 33829 3285   D   -1/
      DATA ARY01(9)/ -0.64486 49130 15404 481    D   -2/
      DATA ARY01(10)/ 0.13128 70828 91002 331    D   -2/
      DATA ARY01(11)/-0.20988 08817 49896 40     D   -3/
      DATA ARY01(12)/ 0.27160 42484 13834 7      D   -4/
      DATA ARY01(13)/-0.29119 91140 14694        D   -5/
      DATA ARY01(14)/ 0.26344 33309 3795         D   -6/
      DATA ARY01(15)/-0.20411 72069 780          D   -7/
      DATA ARY01(16)/ 0.13712 47813 17           D   -8/
      DATA ARY01(17)/-0.80706 80792              D  -10/
      DATA ARY01(18)/ 0.41988 3057               D  -11/
      DATA ARY01(19)/-0.19459 104                D  -12/
      DATA ARY01(20)/ 0.80878 2                  D  -14/
      DATA ARY01(21)/-0.30329                    D  -15/
      DATA ARY01(22)/ 0.1032                     D  -16/
      DATA ARY01(23)/-0.32                       D  -18/
      DATA ARY01(24)/ 0.1                        D  -19/
      DATA ARY0A1(0)/  1.24030 13303 75189 70827  D    0/
      DATA ARY0A1(1)/ -0.47812 53536 32280 693    D   -2/
      DATA ARY0A1(2)/  0.66131 48891 70667 8      D   -4/
      DATA ARY0A1(3)/ -0.18604 27404 86349        D   -5/
      DATA ARY0A1(4)/  0.83627 35565 080          D   -7/
      DATA ARY0A1(5)/ -0.52585 70367 31           D   -8/
      DATA ARY0A1(6)/  0.42606 36325 1            D   -9/
      DATA ARY0A1(7)/ -0.42117 61024              D  -10/
      DATA ARY0A1(8)/  0.48894 6426               D  -11/
      DATA ARY0A1(9)/ -0.64834 929                D  -12/
      DATA ARY0A1(10)/ 0.96172 34                 D  -13/
      DATA ARY0A1(11)/-0.15703 67                 D  -13/
      DATA ARY0A1(12)/ 0.27871 2                  D  -14/
      DATA ARY0A1(13)/-0.53222                    D  -15/
      DATA ARY0A1(14)/ 0.10844                    D  -15/
      DATA ARY0A1(15)/-0.2342                     D  -16/
      DATA ARY0A1(16)/ 0.533                      D  -17/
      DATA ARY0A1(17)/-0.127                      D  -17/
      DATA ARY0A1(18)/ 0.32                       D  -18/
      DATA ARY0A1(19)/-0.8                        D  -19/
      DATA ARY0A1(20)/ 0.2                        D  -19/
      DATA ARY0A1(21)/-0.1                        D  -19/
      DATA ARY0A2(0)/  1.99616 09630 13416 75339  D    0/
      DATA ARY0A2(1)/ -0.19037 98192 46668 161    D   -2/
      DATA ARY0A2(2)/  0.15397 10927 04422 6      D   -4/
      DATA ARY0A2(3)/ -0.31145 08832 8103         D   -6/
      DATA ARY0A2(4)/  0.11108 50971 321          D   -7/
      DATA ARY0A2(5)/ -0.58666 78712 3            D   -9/
      DATA ARY0A2(6)/  0.41399 26949              D  -10/
      DATA ARY0A2(7)/ -0.36539 8763               D  -11/
      DATA ARY0A2(8)/  0.38557 568                D  -12/
      DATA ARY0A2(9)/ -0.47098 00                 D  -13/
      DATA ARY0A2(10)/ 0.65022 0                  D  -14/
      DATA ARY0A2(11)/-0.99624                    D  -15/
      DATA ARY0A2(12)/ 0.16700                    D  -15/
      DATA ARY0A2(13)/-0.3028                     D  -16/
      DATA ARY0A2(14)/ 0.589                      D  -17/
      DATA ARY0A2(15)/-0.122                      D  -17/
      DATA ARY0A2(16)/ 0.27                       D  -18/
      DATA ARY0A2(17)/-0.6                        D  -19/
      DATA ARY0A2(18)/ 0.1                        D  -19/
C
C   Machine-dependent constants (suitable for IEEE machines)
C
      DATA NTERM1,NTERM2,NTERM3,NTERM4/22,22,17,15/
      DATA XLOW,XHIGH/3.16101364D-8,9.007199256D15/
C
C   Start computation
C
      X = XVALUE
C
C   First error test
C
      IF ( X .LT. ZERO ) THEN
         CALL ERRPRN(FNNAME,ERMSG1)
         Y0INT = ZERO
         RETURN
      ENDIF
C
C   Second error test
C
      IF ( X .GT. XHIGH ) THEN
         CALL ERRPRN(FNNAME,ERMSG2)
         Y0INT = ZERO
         RETURN
      ENDIF
C
C   Code for 0 <= x <= 16
C
      IF ( X .LE. SIXTEN ) THEN
         IF ( X .LT. XLOW ) THEN
            IF ( X .EQ. ZERO ) THEN
               Y0INT = ZERO
            ELSE 
               Y0INT = ( LOG(X) + GAL2M1 ) * TWOBPI * X
            ENDIF
         ELSE
            T = X * X / ONE28 - ONE
            TEMP = ( LOG(X) + GAMLN2 ) * CHEVAL(NTERM1,ARJ01,T)
            TEMP = TEMP - CHEVAL(NTERM2,ARY01,T)
            Y0INT = TWOBPI * X * TEMP
         ENDIF
      ELSE
C
C   Code for x > 16
C
         T = FIVE12 / ( X * X ) - ONE
         PIB41 = PIB411 / PIB412
         XMPI4 = ( X - PIB41 ) - PIB42
         TEMP = SIN(XMPI4) * CHEVAL(NTERM3,ARY0A1,T) / X
         TEMP = TEMP + COS(XMPI4) * CHEVAL(NTERM4,ARY0A2,T)
         Y0INT = - RT2BPI * TEMP / SQRT(X)
      ENDIF
      RETURN
      END

      DOUBLE PRECISION FUNCTION CHEVAL(N,A,T)

c*********************************************************************72
c
cc CHEVAL evaluates a Chebyshev series.
c
C   This function evaluates a Chebyshev series, using the
C   Clenshaw method with Reinsch modification, as analysed
C   in the paper by Oliver.
C
C   INPUT PARAMETERS
C
C       N - INTEGER - The no. of terms in the sequence
C
C       A - DOUBLE PRECISION ARRAY, dimension 0 to N - The coefficients of
C           the Chebyshev series
C
C       T - DOUBLE PRECISION - The value at which the series is to be
C           evaluated
C
C
C   REFERENCES
C
C        "An error analysis of the modified Clenshaw method for
C         evaluating Chebyshev and Fourier series" J. Oliver,
C         J.I.M.A., vol. 20, 1977, pp379-391
C
C
C MACHINE-DEPENDENT CONSTANTS: NONE
C
C
C INTRINSIC FUNCTIONS USED;
C
C    ABS
C
C
C AUTHOR:  Dr. Allan J. MacLeod,
C          Dept. of Mathematics and Statistics,
C          University of Paisley ,
C          High St.,
C          PAISLEY,
C          SCOTLAND
C
C
C LATEST MODIFICATION:   21 December , 1992
C
C
      INTEGER I,N
      DOUBLE PRECISION A(0:N),D1,D2,HALF,T,TEST,TT,TWO,U0,U1,U2,ZERO
      DATA ZERO,HALF/ 0.0 D 0 , 0.5 D 0 /
      DATA TEST,TWO/ 0.6 D 0 , 2.0 D 0 /
      U1 = ZERO
C
C   If ABS ( T )  < 0.6 use the standard Clenshaw method
C
      IF ( ABS( T ) .LT. TEST ) THEN
         U0 = ZERO
         TT = T + T
         DO 100 I = N , 0 , -1
            U2 = U1
            U1 = U0
            U0 = TT * U1 + A( I ) - U2
 100     CONTINUE
         CHEVAL =  ( U0 - U2 ) / TWO
      ELSE
C
C   If ABS ( T )  > =  0.6 use the Reinsch modification
C
         D1 = ZERO
C
C   T > =  0.6 code
C
         IF ( T .GT. ZERO ) THEN
            TT =  ( T - HALF ) - HALF
            TT = TT + TT
            DO 200 I = N , 0 , -1
               D2 = D1
               U2 = U1
               D1 = TT * U2 + A( I ) + D2
               U1 = D1 + U2
 200        CONTINUE
            CHEVAL =  ( D1 + D2 ) / TWO
         ELSE
C
C   T < =  -0.6 code
C
            TT =  ( T + HALF ) + HALF
            TT = TT + TT
            DO 300 I = N , 0 , -1
               D2 = D1
               U2 = U1
               D1 = TT * U2 + A( I ) - D2
               U1 = D1 - U2
 300        CONTINUE
            CHEVAL =  ( D1 - D2 ) / TWO
         ENDIF
      ENDIF
      RETURN
      END
      SUBROUTINE ERRPRN(FNNAME,ERRMSG)

c*********************************************************************72
c
cc ERRPRN prints an error message.
c
C   DESCRIPTION:
C      This subroutine prints out an error message if
C      an error has occurred in one of the MISCFUN 
C      functions. 
C
C
C   INPUT PARAMETERS:
C
C      FNNAME - CHARACTER - The name of the function with the error.
C
C      ERRMSG - CHARACTER - The message to be printed out.
C
C
C   MACHINE-DEPENDENT PARAMETER:
C
C      OUTSTR - INTEGER - The numerical value of the output
C                         stream to be used for printing the
C                         error message. The subroutine has the
C                         default value   OUTSTR = 6.
C
C
C   AUTHOR:
C         DR. ALLAN J. MACLEOD,
C         DEPT. OF MATHEMATICS AND STATISTICS,
C         UNIVERSITY OF PAISLEY,
C         HIGH ST.,
C         PAISLEY
C         SCOTLAND
C         PA1 2BE
C
C         (e-mail: macl_ms0@paisley.ac.uk )
C
C
C   LATEST REVISION: 
C                   2 JUNE, 1995
C
      INTEGER OUTSTR
      CHARACTER FNNAME*6,ERRMSG*(*)
      DATA OUTSTR/6/
      WRITE(OUTSTR,1000)FNNAME
      WRITE(OUTSTR,2000)ERRMSG
 1000 FORMAT(/5X,'ERROR IN MISCFUN FUNCTION  ',A6)
 2000 FORMAT(/5X,A50)
      RETURN
      END

