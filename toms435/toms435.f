      FUNCTION DGAMMA ( X )

c*********************************************************************72
C
cc DGAMMA evaluates the Gamma function.
c
c  Discussion:
c
C    This routine calculates the GAMMA function for a real argument X.
c
C    Computation is based on an algorithm outlined in reference 1.
C    The program uses rational functions that approximate the GAMMA
C    function to at least 20 significant decimal digits.  Coefficients
C    for the approximation over the interval (1,2) are unpublished.
C    Those for the approximation for X .GE. 12 are from reference 2.
C    The accuracy achieved depends on the arithmetic system, the
C    compiler, the intrinsic functions, and proper selection of the
C    machine-dependent constants.
C
C
C Explanation of machine-dependent constants
C
C beta   - radix for the floating-point representation
C maxexp - the smallest positive power of beta that overflows
C XBIG   - the largest argument for which GAMMA(X) is representable
C          in the machine, i.e., the solution to the equation
C                  GAMMA(XBIG) = beta**maxexp
C XINF   - the largest machine representable floating-point number;
C          approximately beta**maxexp
C EPS    - the smallest positive floating-point number such that
C          1.0+EPS .GT. 1.0
C XMININ - the smallest positive floating-point number such that
C          1/XMININ is machine representable
C
C     Approximate values for some important machines are:
C
C                            beta       maxexp        XBIG
C
C CRAY-1         (S.P.)        2         8191        966.961
C Cyber 180/855
C   under NOS    (S.P.)        2         1070        177.803
C IEEE (IBM/XT,
C   SUN, etc.)   (S.P.)        2          128        35.040
C IEEE (IBM/XT,
C   SUN, etc.)   (D.P.)        2         1024        171.624
C IBM 3033       (D.P.)       16           63        57.574
C VAX D-Format   (D.P.)        2          127        34.844
C VAX G-Format   (D.P.)        2         1023        171.489
C
C                            XINF         EPS        XMININ
C
C CRAY-1         (S.P.)   5.45E+2465   7.11E-15    1.84E-2466
C Cyber 180/855
C   under NOS    (S.P.)   1.26E+322    3.55E-15    3.14E-294
C IEEE (IBM/XT,
C   SUN, etc.)   (S.P.)   3.40E+38     1.19E-7     1.18E-38
C IEEE (IBM/XT,
C   SUN, etc.)   (D.P.)   1.79D+308    2.22D-16    2.23D-308
C IBM 3033       (D.P.)   7.23D+75     2.22D-16    1.39D-76
C VAX D-Format   (D.P.)   1.70D+38     1.39D-17    5.88D-39
C VAX G-Format   (D.P.)   8.98D+307    1.11D-16    1.12D-308
C
C
C Error returns
C
C  The program returns the value XINF for singularities or
C     when overflow would occur.  The computation is believed
C     to be free of underflow and overflow.
C
C References: "An Overview of Software Development for Special
C              Functions", W. J. Cody, Lecture Notes in Mathematics,
C              506, Numerical Analysis Dundee, 1975, G. A. Watson
C              (ed.), Springer Verlag, Berlin, 1976.
C
C              Computer Approximations, Hart, Et. Al., Wiley and
C              sons, New York, 1968.
C
C  Latest modification: October 12, 1989
C
C  Authors: W. J. Cody and L. Stoltz
C           Applied Mathematics Division
C           Argonne National Laboratory
C           Argonne, IL 60439
C
C
      double precision dgamma
      INTEGER I
      integer N
      LOGICAL PARITY
      DOUBLE PRECISION 
     1    C,CONV,EPS,FACT,HALF,ONE,P,PI,Q,RES,SQRTPI,SUM,TWELVE,
     2    TWO,X,XBIG,XDEN,XINF,XMININ,XNUM,Y,Y1,YSQ,Z,ZERO
      DIMENSION C(7),P(8),Q(8)
C
C  Mathematical constants
C
      DATA ONE,HALF,TWELVE,TWO,ZERO/1.0D0,0.5D0,12.0D0,2.0D0,0.0D0/,
     1     SQRTPI/0.9189385332046727417803297D0/,
     2     PI/3.1415926535897932384626434D0/
C
C  Machine dependent parameters
C
      DATA XBIG,XMININ,EPS/34.844D0,5.88D-39,1.39D-17/,
     1     XINF/1.70D+38/
C
C  Numerator and denominator coefficients for rational minimax
C     approximation over (1,2).
C
      DATA P/-1.71618513886549492533811D+0,2.47656508055759199108314D+1,
     1       -3.79804256470945635097577D+2,6.29331155312818442661052D+2,
     2       8.66966202790413211295064D+2,-3.14512729688483675254357D+4,
     3       -3.61444134186911729807069D+4,6.64561438202405440627855D+4/
      DATA Q/-3.08402300119738975254353D+1,3.15350626979604161529144D+2,
     1      -1.01515636749021914166146D+3,-3.10777167157231109440444D+3,
     2        2.25381184209801510330112D+4,4.75584627752788110767815D+3,
     3      -1.34659959864969306392456D+5,-1.15132259675553483497211D+5/
C
C  Coefficients for minimax approximation over (12, INF).
C
      DATA C/-1.910444077728D-03,8.4171387781295D-04,
     1     -5.952379913043012D-04,7.93650793500350248D-04,
     2     -2.777777777777681622553D-03,8.333333333333333331554247D-02,
     3      5.7083835261D-03/
C
C  Statement functions for conversion between integer and float
C
      CONV(I) = DBLE(I)
      PARITY = .FALSE.
      FACT = ONE
      N = 0
      Y = X
      IF (Y .LE. ZERO) THEN
C
C  Argument is negative
C
            Y = -X
            Y1 = AINT(Y)
            RES = Y - Y1
            IF (RES .NE. ZERO) THEN
                  IF (Y1 .NE. AINT(Y1*HALF)*TWO) PARITY = .TRUE.
                  FACT = -PI / SIN(PI*RES)
                  Y = Y + ONE
               ELSE
                  RES = XINF
                  DGAMMA = res
                  return
            END IF
      END IF
C
C  Argument is positive
C
      IF (Y .LT. EPS) THEN
C
C  Argument .LT. EPS
C
            IF (Y .GE. XMININ) THEN
                  RES = ONE / Y
               ELSE
                  RES = XINF
                  dgamma = res
                  return
            END IF
         ELSE IF (Y .LT. TWELVE) THEN
            Y1 = Y
            IF (Y .LT. ONE) THEN
C
C  0.0 .LT. argument .LT. 1.0
C
                  Z = Y
                  Y = Y + ONE
               ELSE
C
C  1.0 .LT. argument .LT. 12.0, reduce argument if necessary
C
                  N = INT(Y) - 1
                  Y = Y - CONV(N)
                  Z = Y - ONE
            END IF
C
C  Evaluate approximation for 1.0 .LT. argument .LT. 2.0
C
            XNUM = ZERO
            XDEN = ONE
            DO I = 1, 8
               XNUM = (XNUM + P(I)) * Z
               XDEN = XDEN * Z + Q(I)
            end do
            RES = XNUM / XDEN + ONE
            IF (Y1 .LT. Y) THEN
C
C  Adjust result for case  0.0 .LT. argument .LT. 1.0
C
                  RES = RES / Y1
               ELSE IF (Y1 .GT. Y) THEN
C
C  Adjust result for case  2.0 .LT. argument .LT. 12.0
C
                  DO I = 1, N
                     RES = RES * Y
                     Y = Y + ONE
                  end do
            END IF
         ELSE
C
C  Evaluate for argument .GE. 12.0,
C
            IF (Y .LE. XBIG) THEN
                  YSQ = Y * Y
                  SUM = C(7)
                  DO I = 1, 6
                     SUM = SUM / YSQ + C(I)
                  end do
                  SUM = SUM/Y - Y + SQRTPI
                  SUM = SUM + (Y-HALF) * LOG(Y)
                  RES = EXP ( SUM )
               ELSE
                  RES = XINF
                  dgamma = res
                  return
            END IF
      END IF
C
C  Final adjustments and return
C
      IF ( PARITY ) then
        RES = -RES
      end if

      IF (FACT .NE. ONE) then
        RES = FACT / RES
      end if

      DGAMMA = RES

      RETURN
      END
      FUNCTION GAMINC ( A, X1, X2, GAM )

c*********************************************************************72
C
cc GAMINC evaluates the modified incomplete Gamma function.
c
C  COMPUTE THE DIFFERENCE BETWEEN TWO MODIFIED INCOMPLETE
C  GAMMA FUNCTIONS FOR (A,X1) AND (A,X2), THEN MULTIPLY BY
C  EXP(X1).  THAT IS, COMPUTE THE INTEGRAL OF ABS(X)**(A-1.)
C  * EXP(X1-X) FROM X1 TO X2.  IF X1 .GT. X2, THEN X1 - X2 MUST
C  BE .LE. EXPLIM.
C  EXPLIM CAN BE A MACHINE DEPENDENT CONSTANT WHICH PREVENTS
C  EXPONENTIATION OVER- AND UNDERFLOWS.  IT IS USED HERE TO
C  SUPPRESS THE CALCULATION OF MIGAM(A,X2) WHEN THE VALUE OF
C  MIGAM(A,X2) IS INSIGNIFICANT.  THIS USAGE REQUIRES X2 +
C  EXPLIM .GE. X1.  (MIGAM IS AN ABBREVIATION FOR MODIFIED IN-
C  COMPLETE GAMMA FUNCTION.)
C  GAM IS THE COMPLETE GAMMA FUNCTION OF A SUPPLIED BY THE
C  CALLING PROGRAM.
C
C  FOR X .GT. 5., GAM - MIGAM(A,X) IS COMPUTED WITH A CONTINUED
C  FRACTION APPROXIMATION.  FOR ABS(X) .LE. 1.0, THE INTEGRAL
C  IS TRANSFORMED AND EXP(Q) IS APPROXIMATED WITH A CHEBYSHEV
C  SERIES SO THAT THE NEW INTEGRAL MAY BE DONE ANALYTICALLY.
C  FOR X .GT. -12, AND X .LT. 5. (ABS (X) .GT. 1.0), A CONTIN-
C  UED FRACTION APPROXIMATION IS USED.  FINALLY FOR X .LE.
C  -12., THE ASYMPTOTIC EXPANSION IS USED.
C
C  SGN IS A SWITCH WHICH, IF NONZERO, INDICATES WHETHER GAM
C  SHOULD BE ADDED OR SUBTRACTED FROM AN INTERMEDIATE RESULT.
C
      DATA EXPLIM / 20.0 /

      Z = X1
      SGN = 0.0
      TIM = -1.0
      EXPDIF = 1.0

5     continue

      IF ( Z .NE. 0. ) GO TO 10
      GAM1 = 0.
      SGN = SGN + TIM
      GO TO 40

10    continue

      IF ( Z .LE. 5. ) GO TO 20
C
C  USE EQUATION 10.
C
      GAM1 = - EXPDIF * Z**A / ( Z + (1.-A)/(1.+1./(Z+(2.-A)/(1.+2.
     &  /(Z+(3.-A)/(1.+3./(Z+1.7) ))))))
      GO TO 40

20    continue

      AZ = ABS ( Z )
      IF ( Z .LE. -12. ) GO TO 30
      SGN = SGN + TIM
C
C  USE EQUATION 17.
C
      IF ( AZ .LE. 1. ) GAM1 = EXPDIF * Z / A *AZ**(A-1.)
     &  *(1.       +Z/(A+1.) *(.9999999+Z/(A+2.)
     &  *(.9999999 +Z/(A+3.) *(1.000008+Z/(A+4.)
     &  *(1.000005 +Z/(A+5.) *(.9994316+Z/(A+6.)
     &  *(.9995587 +Z/(A+7.) *(1.031684+Z/(A+8.)
     &  *1.028125))))))))
C
C  USE EQUATIONS 11 AND 12.  EVALUATION MUST BE DONE
C  IN DOUBLE PRECISION IF COMPUTER HAS 32 OR FEWER BITS
C  PER WORD.
C
      IF ( AZ .GT. 1. ) GAM1 = EXPDIF * Z / A * AZ**(A-1.)
     &  /(1.- A    *Z/( A     *(A+ 1.+   Z/((A+ 2.)
     &  *(1.-(A+1.)*Z/((A+ 2.)*(A+ 3.+2.*Z/((A+ 4.)
     &  *(1.-(A+2.)*Z/((A+ 4.)*(A+ 5.+3.*Z/((A+ 6.)
     &  *(1.-(A+3.)*Z/((A+ 6.)*(A+ 7.+4.*Z/((A+ 8.)
     &  *(1.-(A+4.)*Z/((A+ 8.)*(A+ 9.+5.*Z/((A+10.)
     &  *(1.-(A+5.)*Z/((A+10.)*(A+11.+6.*Z/((A+12.)
     &  *(1.-(A+6.)*Z/((A+12.)*(A+13.+7.*Z/((A+14.)
     &  *(1.00150-A*8.95D-5 + Z*(-.0337062+A*.0004182
     &  +Z*(000999294-A*.000104103))) )))) )))) ))))
     &  )))) )))) )))) ))))
      GO TO 40
C
C  USE EQUATION 10 AND SHANK'S E1 PROCESS ONCE.
C
30    continue

      GAM1 = - EXPDIF * AZ**(A-1.)*(1.+(A-1.)*(1.+(A-2.)*
     &  (1.+(A-3.)*(1.*(A-4.)*(1.*(A-5.)/(Z-A+6.))
     &  /Z)/Z)/Z)/Z)

40    continue

      IF ( TIM .le. 0.0 ) then

        GAMINC = GAM1
C
C  IF TRUE, CONTRIBUTION AT X2 IS .LT. 1.E-7 * (CONTR AT X1),
C  PROVIDED X2 .GT. X1.
C
        IF ( ABS ( X1 - X2 ) .LE. EXPLIM ) THEN
          Z = X2
          EXPDIF = EXP ( X1 - X2 )
          TIM = 1.0
          GO TO 5
        END IF

        GAM1 = 0.0

      end if

      GAMINC = GAM1 - GAMINC

      IF ( SGN .NE. 0.0 ) THEN
        GAMINC = GAMINC - SIGN ( GAM * EXP ( X1 ), SGN )
      end if

      RETURN
      END
      subroutine gamma_inc_values ( n_data, a, x, fx )

c*********************************************************************72
c
cc GAMMA_INC_VALUES returns some values of the incomplete Gamma function.
c
c  Discussion:
c
c    The (normalized) incomplete Gamma function P(A,X) is defined as:
c
c      PN(A,X) = 1/Gamma(A) * Integral ( 0 <= T <= X ) T**(A-1) * exp(-T) dT.
c
c    With this definition, for all A and X,
c
c      0 <= PN(A,X) <= 1
c
c    and
c
c      PN(A,INFINITY) = 1.0
c
c    In Mathematica, the function can be evaluated by:
c
c      1 - GammaRegularized[A,X]
c
c  Modified:
c
c    20 November 2004
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Milton Abramowitz and Irene Stegun,
c    Handbook of Mathematical Functions,
c    US Department of Commerce, 1964.
c
c    Stephen Wolfram,
c    The Mathematica Book,
c    Fourth Edition,
c    Wolfram Media / Cambridge University Press, 1999.
c
c  Parameters:
c
c    Input/output, integer N_DATA.  The user sets N_DATA to 0 before the
c    first call.  On each call, the routine increments N_DATA by 1, and
c    returns the corresponding data; when there is no more data, the
c    output value of N_DATA will be 0 again.
c
c    Output, real ( kind = 8 ) A, the parameter of the function.
c
c    Output, real ( kind = 8 ) X, the argument of the function.
c
c    Output, real ( kind = 8 ) FX, the value of the function.
c
      implicit none

      integer n_max
      parameter ( n_max = 20 )

      real a
      real a_vec(n_max)
      real fx
      real fx_vec(n_max)
      integer n_data
      real x
      real x_vec(n_max)

      save a_vec
      save fx_vec
      save x_vec

      data a_vec /
     &  0.10E+00, 
     &  0.10E+00, 
     &  0.10E+00, 
     &  0.50E+00, 
     &  0.50E+00, 
     &  0.50E+00, 
     &  0.10E+01, 
     &  0.10E+01, 
     &  0.10E+01, 
     &  0.11E+01, 
     &  0.11E+01, 
     &  0.11E+01, 
     &  0.20E+01, 
     &  0.20E+01, 
     &  0.20E+01, 
     &  0.60E+01, 
     &  0.60E+01, 
     &  0.11E+02, 
     &  0.26E+02, 
     &  0.41E+02 /
      data fx_vec /
     &  0.7382350532339351E+00, 
     &  0.9083579897300343E+00, 
     &  0.9886559833621947E+00, 
     &  0.3014646416966613E+00, 
     &  0.7793286380801532E+00, 
     &  0.9918490284064973E+00, 
     &  0.9516258196404043E-01, 
     &  0.6321205588285577E+00, 
     &  0.9932620530009145E+00, 
     &  0.7205974576054322E-01, 
     &  0.5891809618706485E+00, 
     &  0.9915368159845525E+00, 
     &  0.1018582711118352E-01, 
     &  0.4421745996289254E+00, 
     &  0.9927049442755639E+00, 
     &  0.4202103819530612E-01, 
     &  0.9796589705830716E+00, 
     &  0.9226039842296429E+00, 
     &  0.4470785799755852E+00, 
     &  0.7444549220718699E+00 /
      data x_vec /
     &  0.30E-01, 
     &  0.30E+00, 
     &  0.15E+01, 
     &  0.75E-01, 
     &  0.75E+00, 
     &  0.35E+01, 
     &  0.10E+00, 
     &  0.10E+01, 
     &  0.50E+01, 
     &  0.10E+00, 
     &  0.10E+01, 
     &  0.50E+01, 
     &  0.15E+00, 
     &  0.15E+01, 
     &  0.70E+01, 
     &  0.25E+01, 
     &  0.12E+02, 
     &  0.16E+02, 
     &  0.25E+02, 
     &  0.45E+02 /

      if ( n_data < 0 ) then
        n_data = 0
      end if

      n_data = n_data + 1

      if ( n_max < n_data ) then
        n_data = 0
        a = 0.0E+00
        x = 0.0E+00
        fx = 0.0E+00
      else
        a = a_vec(n_data)
        x = x_vec(n_data)
        fx = fx_vec(n_data)
      end if

      return
      end
      REAL FUNCTION GAMMA(X)

c*********************************************************************72
C
cc GAMMA evaluates the Gamma function.
c
C This routine calculates the GAMMA function for a real argument X.
C   Computation is based on an algorithm outlined in reference 1.
C   The program uses rational functions that approximate the GAMMA
C   function to at least 20 significant decimal digits.  Coefficients
C   for the approximation over the interval (1,2) are unpublished.
C   Those for the approximation for X .GE. 12 are from reference 2.
C   The accuracy achieved depends on the arithmetic system, the
C   compiler, the intrinsic functions, and proper selection of the
C   machine-dependent constants.
C
C
C Explanation of machine-dependent constants
C
C beta   - radix for the floating-point representation
C maxexp - the smallest positive power of beta that overflows
C XBIG   - the largest argument for which GAMMA(X) is representable
C          in the machine, i.e., the solution to the equation
C                  GAMMA(XBIG) = beta**maxexp
C XINF   - the largest machine representable floating-point number;
C          approximately beta**maxexp
C EPS    - the smallest positive floating-point number such that
C          1.0+EPS .GT. 1.0
C XMININ - the smallest positive floating-point number such that
C          1/XMININ is machine representable
C
C     Approximate values for some important machines are:
C
C                            beta       maxexp        XBIG
C
C CRAY-1         (S.P.)        2         8191        966.961
C Cyber 180/855
C   under NOS    (S.P.)        2         1070        177.803
C IEEE (IBM/XT,
C   SUN, etc.)   (S.P.)        2          128        35.040
C IEEE (IBM/XT,
C   SUN, etc.)   (D.P.)        2         1024        171.624
C IBM 3033       (D.P.)       16           63        57.574
C VAX D-Format   (D.P.)        2          127        34.844
C VAX G-Format   (D.P.)        2         1023        171.489
C
C                            XINF         EPS        XMININ
C
C CRAY-1         (S.P.)   5.45E+2465   7.11E-15    1.84E-2466
C Cyber 180/855
C   under NOS    (S.P.)   1.26E+322    3.55E-15    3.14E-294
C IEEE (IBM/XT,
C   SUN, etc.)   (S.P.)   3.40E+38     1.19E-7     1.18E-38
C IEEE (IBM/XT,
C   SUN, etc.)   (D.P.)   1.79D+308    2.22D-16    2.23D-308
C IBM 3033       (D.P.)   7.23D+75     2.22D-16    1.39D-76
C VAX D-Format   (D.P.)   1.70D+38     1.39D-17    5.88D-39
C VAX G-Format   (D.P.)   8.98D+307    1.11D-16    1.12D-308
C
C
C Error returns
C
C  The program returns the value XINF for singularities or
C     when overflow would occur.  The computation is believed
C     to be free of underflow and overflow.
C
C
C  Intrinsic functions required are:
C
C     INT, DBLE, EXP, LOG, REAL, SIN
C
C
C References: "An Overview of Software Development for Special
C              Functions", W. J. Cody, Lecture Notes in Mathematics,
C              506, Numerical Analysis Dundee, 1975, G. A. Watson
C              (ed.), Springer Verlag, Berlin, 1976.
C
C              Computer Approximations, Hart, Et. Al., Wiley and
C              sons, New York, 1968.
C
C  Latest modification: October 12, 1989
C
C  Authors: W. J. Cody and L. Stoltz
C           Applied Mathematics Division
C           Argonne National Laboratory
C           Argonne, IL 60439
C
C
      INTEGER I,N
      LOGICAL PARITY
      REAL 
     1    C,CONV,EPS,FACT,HALF,ONE,P,PI,Q,RES,SQRTPI,SUM,TWELVE,
     2    TWO,X,XBIG,XDEN,XINF,XMININ,XNUM,Y,Y1,YSQ,Z,ZERO
      DIMENSION C(7),P(8),Q(8)
C
C  Mathematical constants
C
      DATA ONE,HALF,TWELVE,TWO,ZERO/1.0E0,0.5E0,12.0E0,2.0E0,0.0E0/,
     1     SQRTPI/0.9189385332046727417803297E0/,
     2     PI/3.1415926535897932384626434E0/
C
C  Machine dependent parameters
C
      DATA XBIG,XMININ,EPS/35.040E0,1.18E-38,1.19E-7/,
     1     XINF/3.4E38/
C
C  Numerator and denominator coefficients for rational minimax
C     approximation over (1,2).
C
      DATA P/-1.71618513886549492533811E+0,2.47656508055759199108314E+1,
     1       -3.79804256470945635097577E+2,6.29331155312818442661052E+2,
     2       8.66966202790413211295064E+2,-3.14512729688483675254357E+4,
     3       -3.61444134186911729807069E+4,6.64561438202405440627855E+4/
      DATA Q/-3.08402300119738975254353E+1,3.15350626979604161529144E+2,
     1      -1.01515636749021914166146E+3,-3.10777167157231109440444E+3,
     2        2.25381184209801510330112E+4,4.75584627752788110767815E+3,
     3      -1.34659959864969306392456E+5,-1.15132259675553483497211E+5/
C
C  Coefficients for minimax approximation over (12, INF).
C
      DATA C/-1.910444077728E-03,8.4171387781295E-04,
     1     -5.952379913043012E-04,7.93650793500350248E-04,
     2     -2.777777777777681622553E-03,8.333333333333333331554247E-02,
     3      5.7083835261E-03/
C
C  Statement functions for conversion between integer and float
C
      CONV(I) = REAL(I)
      PARITY = .FALSE.
      FACT = ONE
      N = 0
      Y = X
      IF (Y .LE. ZERO) THEN
C
C  Argument is negative
C
            Y = -X
            Y1 = AINT(Y)
            RES = Y - Y1
            IF (RES .NE. ZERO) THEN
                  IF (Y1 .NE. AINT(Y1*HALF)*TWO) PARITY = .TRUE.
                  FACT = -PI / SIN(PI*RES)
                  Y = Y + ONE
               ELSE
                  RES = XINF
                  GO TO 900
            END IF
      END IF
C
C  Argument is positive
C
      IF (Y .LT. EPS) THEN
C
C  Argument .LT. EPS
C
            IF (Y .GE. XMININ) THEN
                  RES = ONE / Y
               ELSE
                  RES = XINF
                  GO TO 900
            END IF
         ELSE IF (Y .LT. TWELVE) THEN
            Y1 = Y
            IF (Y .LT. ONE) THEN
C
C  0.0 .LT. argument .LT. 1.0
C
                  Z = Y
                  Y = Y + ONE
               ELSE
C
C  1.0 .LT. argument .LT. 12.0, reduce argument if necessary
C
                  N = INT(Y) - 1
                  Y = Y - CONV(N)
                  Z = Y - ONE
            END IF
C
C  Evaluate approximation for 1.0 .LT. argument .LT. 2.0
C
            XNUM = ZERO
            XDEN = ONE
            DO I = 1, 8
               XNUM = (XNUM + P(I)) * Z
               XDEN = XDEN * Z + Q(I)
            end do
            RES = XNUM / XDEN + ONE
            IF (Y1 .LT. Y) THEN
C
C  Adjust result for case  0.0 .LT. argument .LT. 1.0
C
                  RES = RES / Y1
               ELSE IF (Y1 .GT. Y) THEN
C
C  Adjust result for case  2.0 .LT. argument .LT. 12.0
C
                  DO I = 1, N
                     RES = RES * Y
                     Y = Y + ONE
                  end do
            END IF
         ELSE
C
C  Evaluate for argument .GE. 12.0,
C
            IF (Y .LE. XBIG) THEN
                  YSQ = Y * Y
                  SUM = C(7)
                  DO I = 1, 6
                     SUM = SUM / YSQ + C(I)
                  end do
                  SUM = SUM/Y - Y + SQRTPI
                  SUM = SUM + (Y-HALF)*LOG(Y)
                  RES = EXP(SUM)
               ELSE
                  RES = XINF
                  GO TO 900
            END IF
      END IF
C
C  Final adjustments and return
C
      IF (PARITY) RES = -RES
      IF (FACT .NE. ONE) RES = FACT / RES
  900 GAMMA = RES
      RETURN
      END
      subroutine timestamp ( )

c*********************************************************************72
c
cc TIMESTAMP prints out the current YMDHMS date as a timestamp.
c
c  Discussion:
c
c    This FORTRAN77 version is made available for cases where the
c    FORTRAN90 version cannot be used.
c
c  Modified:
c
c    16 September 2005
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    None
c
      implicit none

      character ( len = 8 ) date
      character ( len = 10 ) time

      call date_and_time ( date, time )

      write ( *, '(a8,2x,a10)' ) date, time

      return
      end
