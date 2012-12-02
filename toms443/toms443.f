      FUNCTION WEW_A ( X, EN )

c*********************************************************************72
C
C  ITERATIVE SOLUTION OF X = W * EXP ( W ) WHERE X IS GIVEN.  (NOVEMBER 1970)
C                                                  (REVISED - SEPTEMBER 1971)
C  VERSION A -- CDC 6600 MACHINE ACCURACY.
C
C  INPUT PARAMETER:
C    X  ARGUMENT OF W(X)
C
C  OUTPUT PARAMETERS:
C    WEW  THE DESIRED SOLUTION.
C    EN   THE LAST RELATIVE CORRECTION TO W(X).
C
C  SET CONSTANTS...
C
      SAVE C1, C2, C3, C4, NEWE

      DATA NEWE / 1 /

      IF ( NEWE ) 10, 20, 10
10    NEWE = 0
      C1 = 4. / 3.
      C2 = 7. / 3.
      C3 = 5. / 6.
      C4 = 2. / 3.
C
C  COMPUTE INITIAL GUESS...
C
20    FLOGX = ALOG ( X )
      IF ( X - 6.46 ) 30, 30, 40
30    WN = X * ( 1. + C1 * X ) / ( 1. + X * ( C2 + C3 * X ) )
      ZN = FLOGX - WN - ALOG ( WN )
      GO TO 50
40    WN = FLOGX
      ZN = -ALOG ( WN )
50    CONTINUE
C
C  ITERATION ONE...
C
      TEMP = 1. + WN
      Y = 2. * TEMP * ( TEMP + C4 * ZN ) - ZN
      WN = WN * ( 1. + ZN * Y / ( TEMP * ( Y - ZN ) ) )
C
C  ITERATION TWO...
C
      ZN = FLOGX - WN - ALOG ( WN )
      TEMP = 1. + WN
      TEMP2 = TEMP + C4 * ZN
      EN = ZN * TEMP2 / ( TEMP * TEMP2 - .5 * ZN )
      WN = WN * ( 1. + EN )
C
C  RETURN...
C
      WEW_A = WN
      RETURN
      END
      FUNCTION WEW_B ( X, EN )

c*********************************************************************72
C
C  ITERATIVE SOLUTION OF X = W * EXP ( W ) WHERE X IS GIVEN.  (NOVEMBER 1970)
C                                                  (REVISED - SEPTEMBER 1971)
C  VERSION B -- MAXIMUM RELATIVE ERROR 3.E-7.
C
C  INPUT PARAMETER:
C    X  ARGUMENT OF W(X)
C
C  OUTPUT PARAMETERS:
C    WEW  THE DESIRED SOLUTION.
C    EN   THE LAST RELATIVE CORRECTION TO W(X).
C
C  SET CONSTANTS...
C
      EQUIVALENCE ( F, FLOGX )

      SAVE C1, C2, C3, C4, NEWE

      DATA NEWE / 1 /

      IF ( NEWE ) 10, 20, 10
10    NEWE = 0
      C1 = 4. / 3.
      C2 = 7. / 3.
      C3 = 5. / 6.
      C4 = 2. / 3.
C
C  COMPUTE INITIAL GUESS...
C
20    FLOGX = ALOG ( X )
      IF ( X - .7385 ) 30, 30, 40
30    WN = X * ( 1. + C1 * X ) / ( 1. + X * ( C2 + C3 * X ) )
      GO TO 50
40    WN = F - 24. * ((F+2.) * F-3.) / (( .7 * F + 58. ) * F + 127. )
50    CONTINUE
C
C  ITERATION ONE...
C
      ZN = FLOGX - WN - ALOG ( WN )
      TEMP = 1. + WN
      Y = 2. * TEMP * ( TEMP + C4 * ZN ) - ZN
      EN = ZN * Y / ( TEMP * ( Y - ZN ) )
      WN = WN * ( 1. + EN )
C
C  RETURN...
C
      WEW_B = WN

      RETURN
      END
