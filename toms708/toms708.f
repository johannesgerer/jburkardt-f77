      function algdiv ( a, b )

c*********************************************************************72
C
cc ALGDIV computes LN(GAMMA(B)/GAMMA(A+B)) when 8 <= B.
C
C
C     IN THIS ALGORITHM, DEL(X) IS THE FUNCTION DEFINED BY
C     LN(GAMMA(X)) = (X - 0.5)*LN(X) - X + 0.5*LN(2*PI) + DEL(X).
C
c  Reference:
c
c    Armido Didonato, Alfred Morris,
c    Algorithm 708: 
c    Significant Digit Computation of the Incomplete Beta Function Ratios,
c    ACM Transactions on Mathematical Software,
c    Volume 18, Number 3, 1992, pages 360-373.
c
      real algdiv

      DATA C0/.833333333333333E-01/, C1/-.277777777760991E-02/,
     *     C2/.793650666825390E-03/, C3/-.595202931351870E-03/,
     *     C4/.837308034031215E-03/, C5/-.165322962780713E-02/

      IF (A .LE. B) GO TO 10
         H = B/A
         C = 1.0/(1.0 + H)
         X = H/(1.0 + H)
         D = A + (B - 0.5)
         GO TO 20
   10 H = A/B
      C = H/(1.0 + H)
      X = 1.0/(1.0 + H)
      D = B + (A - 0.5)
C
C  SET SN = (1 - X**N)/(1 - X)
C
   20 X2 = X*X
      S3 = 1.0 + (X + X2)
      S5 = 1.0 + (X + X2*S3)
      S7 = 1.0 + (X + X2*S5)
      S9 = 1.0 + (X + X2*S7)
      S11 = 1.0 + (X + X2*S9) 
C
C  SET W = DEL(B) - DEL(A + B)
C
      T = (1.0/B)**2
      W = ((((C5*S11*T + C4*S9)*T + C3*S7)*T + C2*S5)*T + C1*S3)*T + C0
      W = W*(C/B)
C
C  Combine the results
c
      u = d*alnrel(a/b)
      v = a*(alog(b) - 1.0)
      if (u .le. v) go to 30
         algdiv = (w - v) - u 
         return
   30 algdiv = (w - u) - v
      return
      end 
      function alnrel ( a ) 

c*********************************************************************72
c
cc ALNREL evaluates the function LN(1 + A).
c
c  Reference:
c
c    Armido Didonato, Alfred Morris,
c    Algorithm 708: 
c    Significant Digit Computation of the Incomplete Beta Function Ratios,
c    ACM Transactions on Mathematical Software,
c    Volume 18, Number 3, 1992, pages 360-373.
c
      real alnrel

      data p1/-.129418923021993e+01/, p2/.405303492862024e+00/,
     *     p3/-.178874546012214e-01/
      data q1/-.162752256355323e+01/, q2/.747811014037616e+00/,
     *     q3/-.845104217945565e-01/

      if (abs(a) .gt. 0.375) go to 10
      t = a/(a + 2.0)
      t2 = t*t
      w = (((p3*t2 + p2)*t2 + p1)*t2 + 1.0)/
     *    (((q3*t2 + q2)*t2 + q1)*t2 + 1.0)
      alnrel = 2.0*t*w
      return

   10 x = 1.d0 + dble(a)
      alnrel = alog(x)
      return
      end 
      function apser ( a, b, x, eps )

c*********************************************************************72
c
cc APSER evaluates I(1-X)(B,A) for A very small.
c
C     APSER YIELDS THE INCOMPLETE BETA RATIO I(SUB(1-X))(B,A) FOR
C     A .LE. MIN(EPS,EPS*B), B*X .LE. 1, AND X .LE. 0.5. USED WHEN
C     A IS VERY SMALL. USE ONLY IF ABOVE INEQUALITIES ARE SATISFIED.
C
c  Reference:
c
c    Armido Didonato, Alfred Morris,
c    Algorithm 708: 
c    Significant Digit Computation of the Incomplete Beta Function Ratios,
c    ACM Transactions on Mathematical Software,
c    Volume 18, Number 3, 1992, pages 360-373.
c
      real apser
      real g
      real j

      data g/.577215664901533/

      bx = b*x
      t = x - bx
      if (b*eps .gt. 2.e-2) go to 10
         c = alog(x) + psi(b) + g + t
         go to 20
   10 c = alog(bx) + g + t

   20 tol = 5.0*eps*abs(c)
      j = 1.0
      s = 0.0
   30    j = j + 1.0
         t = t*(x - bx/j)
         aj = t/j
         s = s + aj 
         if (abs(aj) .gt. tol) go to 30 

      apser = -a*(c + s)

      return
      end 
      function basym ( a, b, lambda, eps )

c*********************************************************************72
c
cc BASYM performs an asymptotic expansion for IX(A,B) for large A and B.
c
C     LAMBDA = (A + B)*Y - B  AND EPS IS THE TOLERANCE USED.
C     IT IS ASSUMED THAT LAMBDA IS NONNEGATIVE AND THAT
C     A AND B ARE GREATER THAN OR EQUAL TO 15.
C
c  Reference:
c
c    Armido Didonato, Alfred Morris,
c    Algorithm 708: 
c    Significant Digit Computation of the Incomplete Beta Function Ratios,
c    ACM Transactions on Mathematical Software,
c    Volume 18, Number 3, 1992, pages 360-373.
c
      real basym
      REAL J0, J1, LAMBDA
      REAL A0(21), B0(21), C(21), D(21) 
C
C  NUM IS THE MAXIMUM VALUE THAT N CAN TAKE IN THE DO LOOP
C  ENDING AT STATEMENT 50. IT IS REQUIRED THAT NUM BE EVEN. 
C  THE ARRAYS A0, B0, C, D HAVE DIMENSION NUM + 1.
C
      DATA NUM/20/
C
C  E0 = 2/SQRT(PI)
C  E1 = 2**(-3/2)
C
      data e0/1.12837916709551/, e1/.353553390593274/

      basym = 0.0
      if (a .ge. b) go to 10
         h = a/b
         r0 = 1.0/(1.0 + h)
         r1 = (b - a)/b
         w0 = 1.0/sqrt(a*(1.0 + h))
         go to 20
   10 h = b/a
      r0 = 1.0/(1.0 + h)
      r1 = (b - a)/a
      w0 = 1.0/sqrt(b*(1.0 + h))

   20 f = a*rlog1(-lambda/a) + b*rlog1(lambda/b)
      t = exp(-f)
      if (t .eq. 0.0) return
      z0 = sqrt(f)
      z = 0.5*(z0/e1)
      z2 = f + f

      a0(1) = (2.0/3.0)*r1
      c(1) = - 0.5*a0(1)
      d(1) = - c(1) 
      j0 = (0.5/e0)*erfc1(1,z0)
      j1 = e1
      sum = j0 + d(1)*w0*j1

      s = 1.0
      h2 = h*h
      hn = 1.0
      w = w0
      znm1 = z
      zn = z2
      do 50 n = 2, num, 2
         hn = h2*hn 
         a0(n) = 2.0*r0*(1.0 + h*hn)/(n + 2.0)
         np1 = n + 1
         s = s + hn 
         a0(np1) = 2.0*r1*s/(n + 3.0)

         do 41 i = n, np1
         r = -0.5*(i + 1.0)
         b0(1) = r*a0(1)
         do 31 m = 2, i
            bsum = 0.0
            mm1 = m - 1
            do 30 j = 1, mm1
               mmj = m - j
   30          bsum = bsum + (j*r - mmj)*a0(j)*b0(mmj)
   31       b0(m) = r*a0(m) + bsum/m
         c(i) = b0(i)/(i + 1.0)

         dsum = 0.0 
         im1 = i - 1
         do j = 1, im1
            imj = i - j
            dsum = dsum + d(imj)*c(j)
         end do

   41    d(i) = -(dsum + c(i))

         j0 = e1*znm1 + (n - 1.0)*j0
         j1 = e1*zn + n*j1
         znm1 = z2*znm1
         zn = z2*zn 
         w = w0*w
         t0 = d(n)*w*j0
         w = w0*w
         t1 = d(np1)*w*j1
         sum = sum + (t0 + t1)
         if ((abs(t0) + abs(t1)) .le. eps*sum) go to 60
   50    continue

   60 u = exp(-bcorr(a,b))
      basym = e0*t*u*sum
      return
      end 
      function bcorr ( a0, b0 )

c*********************************************************************72
C
cc BCORR evaluates a correction term for LN(GAMMA(A)).
c
C     EVALUATION OF  DEL(A0) + DEL(B0) - DEL(A0 + B0)  WHERE
C     LN(GAMMA(A)) = (A - 0.5)*LN(A) - A + 0.5*LN(2*PI) + DEL(A).
C     IT IS ASSUMED THAT A0 .GE. 8 AND B0 .GE. 8. 
C
c  Reference:
c
c    Armido Didonato, Alfred Morris,
c    Algorithm 708: 
c    Significant Digit Computation of the Incomplete Beta Function Ratios,
c    ACM Transactions on Mathematical Software,
c    Volume 18, Number 3, 1992, pages 360-373.
c
      real bcorr

      DATA C0/.833333333333333E-01/, C1/-.277777777760991E-02/,
     *     C2/.793650666825390E-03/, C3/-.595202931351870E-03/,
     *     C4/.837308034031215E-03/, C5/-.165322962780713E-02/

      A = AMIN1(A0, B0)
      B = AMAX1(A0, B0)

      H = A/B
      C = H/(1.0 + H)
      X = 1.0/(1.0 + H)
      X2 = X*X
C
C  SET SN = (1 - X**N)/(1 - X)
C
      S3 = 1.0 + (X + X2)
      S5 = 1.0 + (X + X2*S3)
      S7 = 1.0 + (X + X2*S5)
      S9 = 1.0 + (X + X2*S7)
      S11 = 1.0 + (X + X2*S9) 
C
C  SET W = DEL(B) - DEL(A + B)
C
      T = (1.0/B)**2
      W = ((((C5*S11*T + C4*S9)*T + C3*S7)*T + C2*S5)*T + C1*S3)*T + C0
      W = W*(C/B)
C
C  COMPUTE  DEL(A) + W 
C
      t = (1.0/a)**2
      bcorr = (((((c5*t + c4)*t + c3)*t + c2)*t + c1)*t + c0)/a + w

      return
      end
      subroutine beta_cdf_values ( n_data, a, b, x, fx )

c*********************************************************************72
c
cc BETA_CDF_VALUES returns some values of the Beta CDF.
c
c  Discussion:
c
c    The incomplete Beta function may be written
c
c      BETA_INC(A,B,X) = Integral (0 to X) T**(A-1) * (1-T)**(B-1) dT
c                      / Integral (0 to 1) T**(A-1) * (1-T)**(B-1) dT
c
c    Thus,
c
c      BETA_INC(A,B,0.0) = 0.0
c      BETA_INC(A,B,1.0) = 1.0
c
c    The incomplete Beta function is also sometimes called the
c    "modified" Beta function, or the "normalized" Beta function
c    or the Beta CDF (cumulative density function.
c
c    In Mathematica, the function can be evaluated by:
c
c      BETA[X,A,B] / BETA[A,B]
c
c    The function can also be evaluated by using the Statistics package:
c
c      Needs["Statistics`ContinuousDistributions`"]
c      dist = BetaDistribution [ a, b ]
c      CDF [ dist, x ]
c
c  Modified:
c
c    04 January 2006
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
c    Karl Pearson,
c    Tables of the Incomplete Beta Function,
c    Cambridge University Press, 1968.
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
c    Output, real A, B, the parameters of the function.
c
c    Output, real X, the argument of the function.
c
c    Output, real FX, the value of the function.
c
      implicit none

      integer n_max
      parameter ( n_max = 42 )

      real a
      real a_vec(n_max)
      real b
      real b_vec(n_max)
      real fx
      real fx_vec(n_max)
      integer n_data
      real x
      real x_vec(n_max)

      save a_vec
      save b_vec
      save fx_vec
      save x_vec

      data a_vec /
     &   0.5E+00, 
     &   0.5E+00, 
     &   0.5E+00, 
     &   1.0E+00, 
     &   1.0E+00, 
     &   1.0E+00, 
     &   1.0E+00,
     &   1.0E+00, 
     &   2.0E+00, 
     &   2.0E+00, 
     &   2.0E+00, 
     &   2.0E+00, 
     &   2.0E+00, 
     &   2.0E+00, 
     &   2.0E+00, 
     &   2.0E+00, 
     &   2.0E+00, 
     &   5.5E+00, 
     &  10.0E+00, 
     &  10.0E+00, 
     &  10.0E+00, 
     &  10.0E+00, 
     &  20.0E+00, 
     &  20.0E+00, 
     &  20.0E+00, 
     &  20.0E+00, 
     &  20.0E+00, 
     &  30.0E+00, 
     &  30.0E+00, 
     &  40.0E+00, 
     &  0.1E+01, 
     &   0.1E+01, 
     &   0.1E+01, 
     &   0.1E+01, 
     &   0.1E+01, 
     &   0.1E+01, 
     &   0.1E+01, 
     &   0.1E+01, 
     &   0.2E+01, 
     &   0.3E+01, 
     &   0.4E+01, 
     &   0.5E+01 /
      data b_vec /
     &   0.5E+00, 
     &   0.5E+00, 
     &   0.5E+00, 
     &   0.5E+00, 
     &   0.5E+00, 
     &   0.5E+00, 
     &   0.5E+00, 
     &   1.0E+00, 
     &   2.0E+00, 
     &   2.0E+00, 
     &   2.0E+00, 
     &   2.0E+00, 
     &   2.0E+00, 
     &   2.0E+00, 
     &   2.0E+00, 
     &   2.0E+00, 
     &   2.0E+00, 
     &   5.0E+00, 
     &   0.5E+00, 
     &   5.0E+00, 
     &   5.0E+00, 
     &  10.0E+00, 
     &   5.0E+00, 
     &  10.0E+00, 
     &  10.0E+00, 
     &  20.0E+00, 
     &  20.0E+00, 
     &  10.0E+00, 
     &  10.0E+00, 
     &  20.0E+00, 
     &   0.5E+00, 
     &   0.5E+00, 
     &   0.5E+00, 
     &   0.5E+00, 
     &   0.2E+01, 
     &   0.3E+01, 
     &   0.4E+01, 
     &   0.5E+01, 
     &   0.2E+01, 
     &   0.2E+01, 
     &   0.2E+01, 
     &   0.2E+01 /
      data fx_vec /
     &  0.6376856085851985E-01, 
     &  0.2048327646991335E+00, 
     &  0.1000000000000000E+01, 
     &  0.0000000000000000E+00, 
     &  0.5012562893380045E-02, 
     &  0.5131670194948620E-01, 
     &  0.2928932188134525E+00, 
     &  0.5000000000000000E+00, 
     &  0.2800000000000000E-01, 
     &  0.1040000000000000E+00, 
     &  0.2160000000000000E+00, 
     &  0.3520000000000000E+00, 
     &  0.5000000000000000E+00, 
     &  0.6480000000000000E+00, 
     &  0.7840000000000000E+00, 
     &  0.8960000000000000E+00, 
     &  0.9720000000000000E+00, 
     &  0.4361908850559777E+00, 
     &  0.1516409096347099E+00, 
     &  0.8978271484375000E-01, 
     &  0.1000000000000000E+01, 
     &  0.5000000000000000E+00, 
     &  0.4598773297575791E+00, 
     &  0.2146816102371739E+00, 
     &  0.9507364826957875E+00, 
     &  0.5000000000000000E+00, 
     &  0.8979413687105918E+00, 
     &  0.2241297491808366E+00, 
     &  0.7586405487192086E+00, 
     &  0.7001783247477069E+00, 
     &  0.5131670194948620E-01, 
     &  0.1055728090000841E+00, 
     &  0.1633399734659245E+00, 
     &  0.2254033307585166E+00, 
     &  0.3600000000000000E+00, 
     &  0.4880000000000000E+00, 
     &  0.5904000000000000E+00, 
     &  0.6723200000000000E+00, 
     &  0.2160000000000000E+00, 
     &  0.8370000000000000E-01, 
     &  0.3078000000000000E-01, 
     &  0.1093500000000000E-01 /
      data x_vec /
     &  0.01E+00, 
     &  0.10E+00, 
     &  1.00E+00, 
     &  0.00E+00, 
     &  0.01E+00, 
     &  0.10E+00, 
     &  0.50E+00, 
     &  0.50E+00, 
     &  0.10E+00, 
     &  0.20E+00, 
     &  0.30E+00, 
     &  0.40E+00, 
     &  0.50E+00, 
     &  0.60E+00, 
     &  0.70E+00, 
     &  0.80E+00, 
     &  0.90E+00, 
     &  0.50E+00, 
     &  0.90E+00, 
     &  0.50E+00, 
     &  1.00E+00, 
     &  0.50E+00, 
     &  0.80E+00, 
     &  0.60E+00, 
     &  0.80E+00, 
     &  0.50E+00, 
     &  0.60E+00, 
     &  0.70E+00, 
     &  0.80E+00, 
     &  0.70E+00, 
     &  0.10E+00, 
     &  0.20E+00, 
     &  0.30E+00, 
     &  0.40E+00, 
     &  0.20E+00, 
     &  0.20E+00, 
     &  0.20E+00, 
     &  0.20E+00, 
     &  0.30E+00, 
     &  0.30E+00, 
     &  0.30E+00, 
     &  0.30E+00 /

      if ( n_data < 0 ) then
        n_data = 0
      end if

      n_data = n_data + 1

      if ( n_max < n_data ) then
        n_data = 0
        a = 0.0E+00
        b = 0.0E+00
        x = 0.0E+00
        fx = 0.0E+00
      else
        a = a_vec(n_data)
        b = b_vec(n_data)
        x = x_vec(n_data)
        fx = fx_vec(n_data)
      end if

      return
      end
      subroutine beta_log_values ( n, x, y, fxy )

c*********************************************************************72
c
cc BETA_LOG_VALUES returns some values of the Beta function for testing.
c
c  Modified:
c
c    15 February 2004
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
c  Parameters:
c
c    Input/output, integer N.
c    On input, if N is 0, the first test data is returned, and N is set
c    to the index of the test data.  On each subsequent call, N is
c    incremented and that test data is returned.  When there is no more
c    test data, N is set to 0.
c
c    Output, real X, Y, the arguments of the function.
c
c    Output, real FXY, the value of the function.
c
      implicit none

      integer nmax
      parameter ( nmax = 17 )

      real fxvec ( nmax )
      real fxy
      integer n
      real x
      real xvec ( nmax )
      real y
      real yvec ( nmax )

      data fxvec /
     &  1.609437912E+00,  
     &  0.9162907319E+00,  
     &  0.5108256238E+00,  
     &  0.2231435513E+00, 
     &  1.609437912E+00,  
     &  0.9162907319E+00,  
     &  0.0000000000E+00, 
     &  -1.791759469E+00, 
     &  -3.401197382E+00, 
     &  -4.941642423E+00,  
     &  -6.445719819E+00,  
     &  -3.737669618E+00, 
     &  -5.123963979E+00, 
     &  -6.222576268E+00,  
     &  -7.138867000E+00, 
     &  -7.927324360E+00, 
     &  -9.393661429E+00 /
      data xvec /
     &  0.2E+00, 
     &  0.4E+00, 
     &  0.6E+00, 
     &  0.8E+00, 
     &  1.0E+00, 
     &  1.0E+00, 
     &  1.0E+00, 
     &  2.0E+00, 
     &  3.0E+00, 
     &  4.0E+00, 
     &  5.0E+00, 
     &  6.0E+00, 
     &  6.0E+00, 
     &  6.0E+00, 
     &  6.0E+00, 
     &  6.0E+00, 
     &  7.0E+00 /
      data yvec /
     &  1.0E+00, 
     &  1.0E+00, 
     &  1.0E+00, 
     &  1.0E+00, 
     &  0.2E+00, 
     &  0.4E+00, 
     &  1.0E+00, 
     &  2.0E+00, 
     &  3.0E+00, 
     &  4.0E+00, 
     &  5.0E+00, 
     &  2.0E+00, 
     &  3.0E+00, 
     &  4.0E+00, 
     &  5.0E+00,
     &  6.0E+00, 
     &  7.0E+00 /

      if ( n < 0 ) then
        n = 0
      end if

      n = n + 1

      if ( nmax < n ) then
        n = 0
        x = 0.0E+00
        y = 0.0E+00
        fxy = 0.0E+00
      else
        x = xvec(n)
        y = yvec(n)
        fxy = fxvec(n)
      end if

      return
      end
      function betaln ( a0, b0 )

c*********************************************************************72
c
cc BETALN evaluates the logarithm of the Beta function.
C
C     E = 0.5*LN(2*PI)
C
c  Reference:
c
c    Armido Didonato, Alfred Morris,
c    Algorithm 708: 
c    Significant Digit Computation of the Incomplete Beta Function Ratios,
c    ACM Transactions on Mathematical Software,
c    Volume 18, Number 3, 1992, pages 360-373.
c
      real betaln

      DATA E /.918938533204673/

      A = AMIN1(A0,B0)
      B = AMAX1(A0,B0)
      IF (A .GE. 8.0) GO TO 60
      IF (A .GE. 1.0) GO TO 20
C
C  PROCEDURE WHEN A .LT. 1
C
      IF (B .GE. 8.0) GO TO 10
         BETALN = GAMLN(A) + (GAMLN(B) - GAMLN(A + B))
         RETURN
   10 BETALN = GAMLN(A) + ALGDIV(A,B)
      RETURN
C
C  PROCEDURE WHEN 1 .LE. A .LT. 8
C
   20 IF (A .GT. 2.0) GO TO 30
      IF (B .GT. 2.0) GO TO 21
         BETALN = GAMLN(A) + GAMLN(B) - GSUMLN(A,B)
         RETURN
   21 W = 0.0
      IF (B .LT. 8.0) GO TO 40
         BETALN = GAMLN(A) + ALGDIV(A,B)
         RETURN
C
C  REDUCTION OF A WHEN B .LE. 1000
C
   30 IF (B .GT. 1000.0) GO TO 50
      N = A - 1.0
      W = 1.0
      DO I = 1,N 
         A = A - 1.0
         H = A/B
         W = W * (H/(1.0 + H))
      end do
      W = ALOG(W)
      IF (B .LT. 8.0) GO TO 40
      BETALN = W + GAMLN(A) + ALGDIV(A,B)
      RETURN
C
C  REDUCTION OF B WHEN B .LT. 8
C
   40 N = B - 1.0
      Z = 1.0
      DO I = 1,N 
         B = B - 1.0
         Z = Z * (B/(A + B))
      end do
      BETALN = W + ALOG(Z) + (GAMLN(A) + (GAMLN(B) - GSUMLN(A,B)))
      RETURN
C
C  REDUCTION OF A WHEN B .GT. 1000
C
   50 N = A - 1.0
      W = 1.0
      DO I = 1,N 
         A = A - 1.0
         W = W * (A/(1.0 + A/B))
      end do
      BETALN = (ALOG(W) - N*ALOG(B)) + (GAMLN(A) + ALGDIV(A,B))
      RETURN
C
C  PROCEDURE WHEN A .GE. 8
C
   60 w = bcorr(a,b)
      h = a/b
      c = h/(1.0 + h)
      u = -(a - 0.5)*alog(c)
      v = b*alnrel(h)
      if (u .le. v) go to 61
         betaln = (((-0.5*alog(b) + e) + w) - v) - u
         return
   61 betaln = (((-0.5*alog(b) + e) + w) - u) - v 

      return
      end 
      function bfrac ( a, b, x, y, lambda, eps )

c*********************************************************************72
c
cc BFRAC: continued fraction expansion for IX(A,B) when A and B are greater than 1.
c
C     IT IS ASSUMED THAT  LAMBDA = (A + B)*Y - B. 
C
c  Reference:
c
c    Armido Didonato, Alfred Morris,
c    Algorithm 708: 
c    Significant Digit Computation of the Incomplete Beta Function Ratios,
c    ACM Transactions on Mathematical Software,
c    Volume 18, Number 3, 1992, pages 360-373.
c
      real bfrac
      real lambda
      real n

      BFRAC = BRCOMP(A,B,X,Y) 
      IF (BFRAC .EQ. 0.0) RETURN

      C = 1.0 + LAMBDA
      C0 = B/A
      C1 = 1.0 + 1.0/A
      YP1 = Y + 1.0 

      N = 0.0
      P = 1.0
      S = A + 1.0
      AN = 0.0
      BN = 1.0
      ANP1 = 1.0
      BNP1 = C/C1
      R = C1/C
C
C  CONTINUED FRACTION CALCULATION 
C
   10    N = N + 1.0
         T = N/A
         W = N*(B - N)*X
         E = A/S
         ALPHA = (P*(P + C0)*E*E)*(W*X) 
         E = (1.0 + T)/(C1 + T + T)
         BETA = N + W/S + E*(C + N*YP1) 
         P = 1.0 + T
         S = S + 2.0
C
C  UPDATE AN, BN, ANP1, AND BNP1
C
         T = ALPHA*AN + BETA*ANP1
         AN = ANP1
         ANP1 = T
         T = ALPHA*BN + BETA*BNP1
         BN = BNP1
         BNP1 = T

         R0 = R
         R = ANP1/BNP1
         IF (ABS(R - R0) .LE. EPS*R) GO TO 20
C
C  RESCALE AN, BN, ANP1, AND BNP1 
C
         AN = AN/BNP1
         BN = BN/BNP1
         ANP1 = R
         BNP1 = 1.0 
         GO TO 10
C
C  Termination 
c
   20 bfrac = bfrac*r

      return
      end 
      subroutine bgrat ( a, b, x, y, w, eps, ierr )

c*********************************************************************72
c
cc BGRAT uses asymptotic expansion for IX(A,B) when B < A.
c
C     ASYMPTOTIC EXPANSION FOR IX(A,B) WHEN A IS LARGER THAN B.
C     THE RESULT OF THE EXPANSION IS ADDED TO W. IT IS ASSUMED
C     THAT A .GE. 15 AND B .LE. 1.  EPS IS THE TOLERANCE USED.
C     IERR IS A VARIABLE THAT REPORTS THE STATUS OF THE RESULTS.
C
c  Reference:
c
c    Armido Didonato, Alfred Morris,
c    Algorithm 708: 
c    Significant Digit Computation of the Incomplete Beta Function Ratios,
c    ACM Transactions on Mathematical Software,
c    Volume 18, Number 3, 1992, pages 360-373.
c
      REAL J, L, LNX, NU, N2
      REAL C(30), D(30)

      BM1 = (B - 0.5) - 0.5
      NU = A + 0.5*BM1
      IF (Y .GT. 0.375) GO TO 10
         LNX = ALNREL(-Y)
         GO TO 11
   10 LNX = ALOG(X) 
   11 Z = -NU*LNX
      IF (B*Z .EQ. 0.0) GO TO 100
C
C  COMPUTATION OF THE EXPANSION
C  SET R = EXP(-Z)*Z**B/GAMMA(B)
C
      R = B*(1.0 + GAM1(B))*EXP(B*ALOG(Z))
      R = R*EXP(A*LNX)*EXP(0.5*BM1*LNX) 
      U = ALGDIV(B,A) + B*ALOG(NU)
      U = R*EXP(-U) 
      IF (U .EQ. 0.0) GO TO 100
      CALL GRAT1(B,Z,R,P,Q,EPS)

      V = 0.25*(1.0/NU)**2
      T2 = 0.25*LNX*LNX
      L = W/U
      J = Q/R
      SUM = J
      T = 1.0
      CN = 1.0
      N2 = 0.0
      DO 22 N = 1,30
         BP2N = B + N2
         J = (BP2N*(BP2N + 1.0)*J + (Z + BP2N + 1.0)*T)*V
         N2 = N2 + 2.0
         T = T*T2
         CN = CN/(N2*(N2 + 1.0))
         C(N) = CN
         S = 0.0
         IF (N .EQ. 1) GO TO 21
            NM1 = N - 1
            COEF = B - N
            DO 20 I = 1,NM1
               S = S + COEF*C(I)*D(N-I) 
   20          COEF = COEF + B
   21    D(N) = BM1*CN + S/N
         DJ = D(N)*J
         SUM = SUM + DJ
         IF (SUM .LE. 0.0) GO TO 100
         IF (ABS(DJ) .LE. EPS*(SUM + L)) GO TO 30 
   22 CONTINUE
C
C  ADD THE RESULTS TO W
C
   30 IERR = 0
      W = W + U*SUM 
      RETURN
C
C  The expansion cannot be computed
c
  100 ierr = 1

      return
      end 
      function bpser ( a, b, x, eps ) 

c*********************************************************************72
c
cc BPSER evaluates IX(A,B) when B <= 1 or B*X <= 0.7.
c
C     POWER SERIES EXPANSION FOR EVALUATING IX(A,B) WHEN B .LE. 1
C     OR B*X .LE. 0.7.  EPS IS THE TOLERANCE USED.
C
c  Reference:
c
c    Armido Didonato, Alfred Morris,
c    Algorithm 708: 
c    Significant Digit Computation of the Incomplete Beta Function Ratios,
c    ACM Transactions on Mathematical Software,
c    Volume 18, Number 3, 1992, pages 360-373.
c
      real bpser
      REAL N

      BPSER = 0.0
      IF (X .EQ. 0.0) RETURN
C
C  COMPUTE THE FACTOR X**A/(A*BETA(A,B))
C
      A0 = AMIN1(A,B)
      IF (A0 .LT. 1.0) GO TO 10
         Z = A*ALOG(X) - BETALN(A,B)
         BPSER = EXP(Z)/A
         GO TO 70
   10 B0 = AMAX1(A,B)
      IF (B0 .GE. 8.0) GO TO 60
      IF (B0 .GT. 1.0) GO TO 40
C
C  PROCEDURE FOR A0 .LT. 1 AND B0 .LE. 1
C
      BPSER = X**A
      IF (BPSER .EQ. 0.0) RETURN

      APB = A + B
      IF (APB .GT. 1.0) GO TO 20
         Z = 1.0 + GAM1(APB)
         GO TO 30
   20 U = DBLE(A) + DBLE(B) - 1.D0
      Z = (1.0 + GAM1(U))/APB 

   30 C = (1.0 + GAM1(A))*(1.0 + GAM1(B))/Z
      BPSER = BPSER*C*(B/APB) 
      GO TO 70
C
C  PROCEDURE FOR A0 .LT. 1 AND 1 .LT. B0 .LT. 8
C
   40 U = GAMLN1(A0)
      M = B0 - 1.0
      IF (M .LT. 1) GO TO 50
      C = 1.0
      DO 41 I = 1,M 
         B0 = B0 - 1.0
   41    C = C*(B0/(A0 + B0)) 
      U = ALOG(C) + U

   50 Z = A*ALOG(X) - U
      B0 = B0 - 1.0 
      APB = A0 + B0 
      IF (APB .GT. 1.0) GO TO 51
         T = 1.0 + GAM1(APB)
         GO TO 52
   51 U = DBLE(A0) + DBLE(B0) - 1.D0
      T = (1.0 + GAM1(U))/APB 
   52 BPSER = EXP(Z)*(A0/A)*(1.0 + GAM1(B0))/T
      GO TO 70
C
C  PROCEDURE FOR A0 .LT. 1 AND B0 .GE. 8
C
   60 U = GAMLN1(A0) + ALGDIV(A0,B0)
      Z = A*ALOG(X) - U
      BPSER = (A0/A)*EXP(Z)
   70 IF (BPSER .EQ. 0.0 .OR. A .LE. 0.1*EPS) RETURN
C
C  Compute the series
c
      sum = 0.0
      n = 0.0
      c = 1.0
      tol = eps/a
  100    n = n + 1.0
         c = c*(0.5 + (0.5 - b/n))*x
         w = c/(a + n)
         sum = sum + w
         if (abs(w) .gt. tol) go to 100 
      bpser = bpser*(1.0 + a*sum)

      return
      end 
      subroutine bratio ( a, b, x, y, w, w1, ierr ) 

c*********************************************************************72
C
cc BRATIO evaluates the incomplete Beta function IX(A,B).
c
c  Discussion:
c
C     IT IS ASSUMED THAT A AND B ARE NONNEGATIVE, AND THAT X .LE. 1
C     AND Y = 1 - X.  BRATIO ASSIGNS W AND W1 THE VALUES
C
C                      W  = IX(A,B)
C                      W1 = 1 - IX(A,B) 
C
C     IERR IS A VARIABLE THAT REPORTS THE STATUS OF THE RESULTS.
C     IF NO INPUT ERRORS ARE DETECTED THEN IERR IS SET TO 0 AND
C     W AND W1 ARE COMPUTED. OTHERWISE, IF AN ERROR IS DETECTED,
C     THEN W AND W1 ARE ASSIGNED THE VALUE 0 AND IERR IS SET TO
C     ONE OF THE FOLLOWING VALUES ...
C
C        IERR = 1  IF A OR B IS NEGATIVE
C        IERR = 2  IF A = B = 0
C        IERR = 3  IF X .LT. 0 OR X .GT. 1
C        IERR = 4  IF Y .LT. 0 OR Y .GT. 1
C        IERR = 5  IF X + Y .NE. 1
C        IERR = 6  IF X = A = 0
C        IERR = 7  IF Y = B = 0
C
C     WRITTEN BY ALFRED H. MORRIS, JR.
C        NAVAL SURFACE WARFARE CENTER
C        DAHLGREN, VIRGINIA
C     REVISED ... NOV 1991
c
c  Reference:
c
c    Armido Didonato, Alfred Morris,
c    Algorithm 708: 
c    Significant Digit Computation of the Incomplete Beta Function Ratios,
c    ACM Transactions on Mathematical Software,
c    Volume 18, Number 3, 1992, pages 360-373.
c
      REAL LAMBDA
c
C  EPS IS A MACHINE DEPENDENT CONSTANT. EPS IS THE SMALLEST 
C  FLOATING POINT NUMBER FOR WHICH 1.0 + EPS .GT. 1.0
C
      EPS = SPMPAR(1)
      W = 0.0
      W1 = 0.0
      IF (A .LT. 0.0 .OR. B .LT. 0.0) GO TO 300
      IF (A .EQ. 0.0 .AND. B .EQ. 0.0) GO TO 310
      IF (X .LT. 0.0 .OR. X .GT. 1.0) GO TO 320
      IF (Y .LT. 0.0 .OR. Y .GT. 1.0) GO TO 330
      Z = ((X + Y) - 0.5) - 0.5
      IF (ABS(Z) .GT. 3.0*EPS) GO TO 340

      IERR = 0
      IF (X .EQ. 0.0) GO TO 200
      IF (Y .EQ. 0.0) GO TO 210
      IF (A .EQ. 0.0) GO TO 211
      IF (B .EQ. 0.0) GO TO 201

      EPS = AMAX1(EPS, 1.E-15)
      IF (AMAX1(A,B) .LT. 1.E-3*EPS) GO TO 230

      IND = 0
      A0 = A
      B0 = B
      X0 = X
      Y0 = Y
      IF (AMIN1(A0, B0) .GT. 1.0) GO TO 30
C
C  PROCEDURE FOR A0 .LE. 1 OR B0 .LE. 1
C
      IF (X .LE. 0.5) GO TO 10
      IND = 1
      A0 = B
      B0 = A
      X0 = Y
      Y0 = X

   10 IF (B0 .LT. AMIN1(EPS,EPS*A0)) GO TO 80
      IF (A0 .LT. AMIN1(EPS,EPS*B0) .AND. B0*X0 .LE. 1.0) GO TO 90
      IF (AMAX1(A0, B0) .GT. 1.0) GO TO 20
      IF (A0 .GE. AMIN1(0.2, B0)) GO TO 100
      IF (X0**A0 .LE. 0.9) GO TO 100
      IF (X0 .GE. 0.3) GO TO 110
      N = 20
      GO TO 130

   20 IF (B0 .LE. 1.0) GO TO 100
      IF (X0 .GE. 0.3) GO TO 110
      IF (X0 .GE. 0.1) GO TO 21
      IF ((X0*B0)**A0 .LE. 0.7) GO TO 100
   21 IF (B0 .GT. 15.0) GO TO 131
      N = 20
      GO TO 130
C
C  PROCEDURE FOR A0 .GT. 1 AND B0 .GT. 1
C
   30 IF (A .GT. B) GO TO 31
         LAMBDA = A - (A + B)*X
         GO TO 32
   31 LAMBDA = (A + B)*Y - B
   32 IF (LAMBDA .GE. 0.0) GO TO 40
      IND = 1
      A0 = B
      B0 = A
      X0 = Y
      Y0 = X
      LAMBDA = ABS(LAMBDA)

   40 IF (B0 .LT. 40.0 .AND. B0*X0 .LE. 0.7) GO TO 100
      IF (B0 .LT. 40.0) GO TO 140
      IF (A0 .GT. B0) GO TO 50
         IF (A0 .LE. 100.0) GO TO 120
         IF (LAMBDA .GT. 0.03*A0) GO TO 120
         GO TO 180
   50 IF (B0 .LE. 100.0) GO TO 120
      IF (LAMBDA .GT. 0.03*B0) GO TO 120
      GO TO 180
C
C  EVALUATION OF THE APPROPRIATE ALGORITHM
C
   80 W = FPSER(A0, B0, X0, EPS)
      W1 = 0.5 + (0.5 - W)
      GO TO 220

   90 W1 = APSER(A0, B0, X0, EPS)
      W = 0.5 + (0.5 - W1)
      GO TO 220

  100 W = BPSER(A0, B0, X0, EPS)
      W1 = 0.5 + (0.5 - W)
      GO TO 220

  110 W1 = BPSER(B0, A0, Y0, EPS)
      W = 0.5 + (0.5 - W1)
      GO TO 220

  120 W = BFRAC(A0, B0, X0, Y0, LAMBDA, 15.0*EPS) 
      W1 = 0.5 + (0.5 - W)
      GO TO 220

  130 W1 = BUP(B0, A0, Y0, X0, N, EPS)
      B0 = B0 + N
  131 CALL BGRAT(B0, A0, Y0, X0, W1, 15.0*EPS, IERR1)
      W = 0.5 + (0.5 - W1)
      GO TO 220

  140 N = B0
      B0 = B0 - N
      IF (B0 .NE. 0.0) GO TO 141
         N = N - 1
         B0 = 1.0
  141 W = BUP(B0, A0, Y0, X0, N, EPS)
      IF (X0 .GT. 0.7) GO TO 150
      W = W + BPSER(A0, B0, X0, EPS)
      W1 = 0.5 + (0.5 - W)
      GO TO 220

  150 IF (A0 .GT. 15.0) GO TO 151
         N = 20
         W = W + BUP(A0, B0, X0, Y0, N, EPS)
         A0 = A0 + N
  151 CALL BGRAT(A0, B0, X0, Y0, W, 15.0*EPS, IERR1)
      W1 = 0.5 + (0.5 - W)
      GO TO 220

  180 W = BASYM(A0, B0, LAMBDA, 100.0*EPS)
      W1 = 0.5 + (0.5 - W)
      GO TO 220
C
C  TERMINATION OF THE PROCEDURE
C
  200 IF (A .EQ. 0.0) GO TO 350
  201 W = 0.0
      W1 = 1.0
      RETURN

  210 IF (B .EQ. 0.0) GO TO 360
  211 W = 1.0
      W1 = 0.0
      RETURN

  220 IF (IND .EQ. 0) RETURN
      T = W
      W = W1
      W1 = T
      RETURN
C
C  PROCEDURE FOR A AND B .LT. 1.E-3*EPS
C
  230 W = B/(A + B) 
      W1 = A/(A + B)
      RETURN
C
C  ERROR RETURN
C
  300 IERR = 1
      RETURN
  310 IERR = 2
      RETURN
  320 IERR = 3
      RETURN
  330 IERR = 4
      RETURN
  340 IERR = 5
      RETURN
  350 IERR = 6
      RETURN
  360 IERR = 7

      return
      end 
      function brcmp1 ( mu, a, b, x, y )

c*********************************************************************72
c
cc BRCMP1 evaluates EXP(MU) * (X**A*Y**B/BETA(A,B)).
C
c  Reference:
c
c    Armido Didonato, Alfred Morris,
c    Algorithm 708: 
c    Significant Digit Computation of the Incomplete Beta Function Ratios,
c    ACM Transactions on Mathematical Software,
c    Volume 18, Number 3, 1992, pages 360-373.
c
      real brcmp1
      REAL LAMBDA
      real LNX
      real LNY
C
C  CONST = 1/SQRT(2*PI)
C
      DATA CONST/.398942280401433/

      A0 = AMIN1(A,B)
      IF (A0 .GE. 8.0) GO TO 100

      IF (X .GT. 0.375) GO TO 10
         LNX = ALOG(X)
         LNY = ALNREL(-X)
         GO TO 20
   10 IF (Y .GT. 0.375) GO TO 11
         LNX = ALNREL(-Y)
         LNY = ALOG(Y)
         GO TO 20
   11 LNX = ALOG(X) 
      LNY = ALOG(Y) 

   20 Z = A*LNX + B*LNY
      IF (A0 .LT. 1.0) GO TO 30
      Z = Z - BETALN(A,B)
      BRCMP1 = ESUM(MU,Z)
      RETURN
C
C  PROCEDURE FOR A .LT. 1 OR B .LT. 1 
C
   30 B0 = AMAX1(A,B)
      IF (B0 .GE. 8.0) GO TO 80
      IF (B0 .GT. 1.0) GO TO 60
C
C  ALGORITHM FOR B0 .LE. 1
C
      BRCMP1 = ESUM(MU,Z)
      IF (BRCMP1 .EQ. 0.0) RETURN

      APB = A + B
      IF (APB .GT. 1.0) GO TO 40
         Z = 1.0 + GAM1(APB)
         GO TO 50
   40 U = DBLE(A) + DBLE(B) - 1.D0
      Z = (1.0 + GAM1(U))/APB 

   50 C = (1.0 + GAM1(A))*(1.0 + GAM1(B))/Z
      BRCMP1 = BRCMP1*(A0*C)/(1.0 + A0/B0)
      RETURN
C
C  ALGORITHM FOR 1 .LT. B0 .LT. 8
C
   60 U = GAMLN1(A0)
      N = B0 - 1.0
      IF (N .LT. 1) GO TO 70
      C = 1.0
      DO I = 1,N 
         B0 = B0 - 1.0
         C = C*(B0/(A0 + B0)) 
      end do
      U = ALOG(C) + U

   70 Z = Z - U
      B0 = B0 - 1.0 
      APB = A0 + B0 
      IF (APB .GT. 1.0) GO TO 71
         T = 1.0 + GAM1(APB)
         GO TO 72
   71 U = DBLE(A0) + DBLE(B0) - 1.D0
      T = (1.0 + GAM1(U))/APB 
   72 BRCMP1 = A0*ESUM(MU,Z)*(1.0 + GAM1(B0))/T
      RETURN
C
C  ALGORITHM FOR B0 .GE. 8
C
   80 U = GAMLN1(A0) + ALGDIV(A0,B0)
      BRCMP1 = A0*ESUM(MU,Z - U)
      RETURN
C
C  PROCEDURE FOR A .GE. 8 AND B .GE. 8
C
  100 if (a .gt. b) go to 101 
         h = a/b
         x0 = h/(1.0 + h)
         y0 = 1.0/(1.0 + h)
         lambda = a - (a + b)*x
         go to 110
  101 h = b/a
      x0 = 1.0/(1.0 + h)
      y0 = h/(1.0 + h)
      lambda = (a + b)*y - b

  110 e = -lambda/a 
      if (abs(e) .gt. 0.6) go to 111
         u = rlog1(e)
         go to 120
  111 u = e - alog(x/x0)

  120 e = lambda/b
      if (abs(e) .gt. 0.6) go to 121
         v = rlog1(e)
         go to 130
  121 v = e - alog(y/y0)

  130 z = esum(mu,-(a*u + b*v))
      brcmp1 = const*sqrt(b*x0)*z*exp(-bcorr(a,b))

      return
      end 
      function brcomp ( a, b, x, y ) 

c*********************************************************************72
c
cc BRCOMP evaluates X**A*Y**B/BETA(A,B).
C
c  Reference:
c
c    Armido Didonato, Alfred Morris,
c    Algorithm 708: 
c    Significant Digit Computation of the Incomplete Beta Function Ratios,
c    ACM Transactions on Mathematical Software,
c    Volume 18, Number 3, 1992, pages 360-373.
c
      real brcomp
      REAL LAMBDA, LNX, LNY
C
C  CONST = 1/SQRT(2*PI)
C
      DATA CONST/.398942280401433/

      BRCOMP = 0.0
      IF (X .EQ. 0.0 .OR. Y .EQ. 0.0) RETURN
      A0 = AMIN1(A,B)
      IF (A0 .GE. 8.0) GO TO 100

      IF (X .GT. 0.375) GO TO 10
         LNX = ALOG(X)
         LNY = ALNREL(-X)
         GO TO 20
   10 IF (Y .GT. 0.375) GO TO 11
         LNX = ALNREL(-Y)
         LNY = ALOG(Y)
         GO TO 20
   11 LNX = ALOG(X) 
      LNY = ALOG(Y) 

   20 Z = A*LNX + B*LNY
      IF (A0 .LT. 1.0) GO TO 30
      Z = Z - BETALN(A,B)
      BRCOMP = EXP(Z)
      RETURN
C
C  PROCEDURE FOR A .LT. 1 OR B .LT. 1 
C
   30 B0 = AMAX1(A,B)
      IF (B0 .GE. 8.0) GO TO 80
      IF (B0 .GT. 1.0) GO TO 60
C
C  ALGORITHM FOR B0 .LE. 1
C
      BRCOMP = EXP(Z)
      IF (BRCOMP .EQ. 0.0) RETURN

      APB = A + B
      IF (APB .GT. 1.0) GO TO 40
         Z = 1.0 + GAM1(APB)
         GO TO 50
   40 U = DBLE(A) + DBLE(B) - 1.D0
      Z = (1.0 + GAM1(U))/APB 

   50 C = (1.0 + GAM1(A))*(1.0 + GAM1(B))/Z
      BRCOMP = BRCOMP*(A0*C)/(1.0 + A0/B0)
      RETURN
C
C  ALGORITHM FOR 1 .LT. B0 .LT. 8
C
   60 U = GAMLN1(A0)
      N = B0 - 1.0
      IF (N .LT. 1) GO TO 70
      C = 1.0
      DO I = 1,N 
         B0 = B0 - 1.0
         C = C*(B0/(A0 + B0)) 
      end do
      U = ALOG(C) + U

   70 Z = Z - U
      B0 = B0 - 1.0 
      APB = A0 + B0 
      IF (APB .GT. 1.0) GO TO 71
         T = 1.0 + GAM1(APB)
         GO TO 72
   71 U = DBLE(A0) + DBLE(B0) - 1.D0
      T = (1.0 + GAM1(U))/APB 
   72 BRCOMP = A0*EXP(Z)*(1.0 + GAM1(B0))/T
      RETURN
C
C  ALGORITHM FOR B0 .GE. 8
C
   80 U = GAMLN1(A0) + ALGDIV(A0,B0)
      BRCOMP = A0*EXP(Z - U)
      RETURN
C
C  PROCEDURE FOR A .GE. 8 AND B .GE. 8
C
  100 if (a .gt. b) go to 101 
         h = a/b
         x0 = h/(1.0 + h)
         y0 = 1.0/(1.0 + h)
         lambda = a - (a + b)*x
         go to 110
  101 h = b/a
      x0 = 1.0/(1.0 + h)
      y0 = h/(1.0 + h)
      lambda = (a + b)*y - b

  110 e = -lambda/a 
      if (abs(e) .gt. 0.6) go to 111
         u = rlog1(e)
         go to 120
  111 u = e - alog(x/x0)

  120 e = lambda/b
      if (abs(e) .gt. 0.6) go to 121
         v = rlog1(e)
         go to 130
  121 v = e - alog(y/y0)

  130 z = exp(-(a*u + b*v))
      brcomp = const*sqrt(b*x0)*z*exp(-bcorr(a,b))
      return
      end 
      function bup ( a, b, x, y, n, eps )

c*********************************************************************72
c
cc BUP evaluates IX(A,B) - IX(A+N,B), where N is a positive integer.
c
c  Reference:
c
c    Armido Didonato, Alfred Morris,
c    Algorithm 708: 
c    Significant Digit Computation of the Incomplete Beta Function Ratios,
c    ACM Transactions on Mathematical Software,
c    Volume 18, Number 3, 1992, pages 360-373.
c
C     EPS IS THE TOLERANCE USED.
C
      real bup
      REAL L
C
C  OBTAIN THE SCALING FACTOR EXP(-MU) AND EXP(MU)*(X**A*Y**B/BETA(A,B))/A
C
      APB = A + B
      AP1 = A + 1.0 
      MU = 0
      D = 1.0
      IF (N .EQ. 1 .OR. A .LT. 1.0) GO TO 10
      IF (APB .LT. 1.1*AP1) GO TO 10
         MU = ABS(EXPARG(1))
         K = EXPARG(0)
         IF (K .LT. MU) MU = K
         T = MU
         D = EXP(-T)

   10 BUP = BRCMP1(MU,A,B,X,Y)/A
      IF (N .EQ. 1 .OR. BUP .EQ. 0.0) RETURN
      NM1 = N - 1
      W = D
C
C  LET K BE THE INDEX OF THE MAXIMUM TERM 
C
      K = 0
      IF (B .LE. 1.0) GO TO 40
      IF (Y .GT. 1.E-4) GO TO 20
         K = NM1
         GO TO 30
   20 R = (B - 1.0)*X/Y - A
      IF (R .LT. 1.0) GO TO 40
      K = NM1
      T = NM1
      IF (R .LT. T) K = R
C
C  ADD THE INCREASING TERMS OF THE SERIES 
C
   30 DO I = 1,K 
         L = I - 1
         D = ((APB + L)/(AP1 + L))*X*D
         W = W + D
      end do
      IF (K .EQ. NM1) GO TO 50
C
C  ADD THE REMAINING TERMS OF THE SERIES
C
   40 KP1 = K + 1
      DO I = KP1,NM1
         L = I - 1
         D = ((APB + L)/(AP1 + L))*X*D
         W = W + D
         IF (D .LE. EPS*W) GO TO 50
      end do
C
C  Terminate the procedure 
c
   50 bup = bup*w

      return
      end 
      function erf ( x )

c*********************************************************************72
c
cc ERF evaluates the error function.
C
c  Reference:
c
c    Armido Didonato, Alfred Morris,
c    Algorithm 708: 
c    Significant Digit Computation of the Incomplete Beta Function Ratios,
c    ACM Transactions on Mathematical Software,
c    Volume 18, Number 3, 1992, pages 360-373.
c
      real erf

      real a(5),b(3),p(8),q(8),r(5),s(4)

      data c /.564189583547756/

      data a(1) /.771058495001320e-04/, a(2)/-.133733772997339e-02/,
     *     a(3) /.323076579225834e-01/, a(4) /.479137145607681e-01/,
     *     a(5) /.128379167095513e+00/
      data b(1) /.301048631703895e-02/, b(2) /.538971687740286e-01/,
     *     b(3) /.375795757275549e+00/

      data p(1)/-1.36864857382717e-07/, p(2) /5.64195517478974e-01/,
     *     p(3) /7.21175825088309e+00/, p(4) /4.31622272220567e+01/,
     *     p(5) /1.52989285046940e+02/, p(6) /3.39320816734344e+02/,
     *     p(7) /4.51918953711873e+02/, p(8) /3.00459261020162e+02/
      data q(1) /1.00000000000000e+00/, q(2) /1.27827273196294e+01/,
     *     q(3) /7.70001529352295e+01/, q(4) /2.77585444743988e+02/,
     *     q(5) /6.38980264465631e+02/, q(6) /9.31354094850610e+02/,
     *     q(7) /7.90950925327898e+02/, q(8) /3.00459260956983e+02/

      data r(1) /2.10144126479064e+00/, r(2) /2.62370141675169e+01/,
     *     r(3) /2.13688200555087e+01/, r(4) /4.65807828718470e+00/,
     *     r(5) /2.82094791773523e-01/
      data s(1) /9.41537750555460e+01/, s(2) /1.87114811799590e+02/,
     *     s(3) /9.90191814623914e+01/, s(4) /1.80124575948747e+01/

      ax = abs(x)
      if (ax .gt. 0.5) go to 10
      t = x*x
      top = ((((a(1)*t + a(2))*t + a(3))*t + a(4))*t + a(5)) + 1.0
      bot = ((b(1)*t + b(2))*t + b(3))*t + 1.0
      erf = x*(top/bot)
      return

   10 if (ax .gt. 4.0) go to 20
      top = ((((((p(1)*ax + p(2))*ax + p(3))*ax + p(4))*ax + p(5))*ax 
     *                    + p(6))*ax + p(7))*ax + p(8)
      bot = ((((((q(1)*ax + q(2))*ax + q(3))*ax + q(4))*ax + q(5))*ax 
     *                    + q(6))*ax + q(7))*ax + q(8)
      erf = 0.5 + (0.5 - exp(-x*x)*top/bot)
      if (x .lt. 0.0) erf = -erf
      return

   20 if (ax .ge. 5.8) go to 30
      x2 = x*x
      t = 1.0/x2
      top = (((r(1)*t + r(2))*t + r(3))*t + r(4))*t + r(5)
      bot = (((s(1)*t + s(2))*t + s(3))*t + s(4))*t + 1.0
      erf = (c - top/(x2*bot)) / ax
      erf = 0.5 + (0.5 - exp(-x2)*erf)
      if (x .lt. 0.0) erf = -erf
      return

   30 erf = sign(1.0,x)

      return
      end
      subroutine erf_values ( n_data, x, fx )

c*********************************************************************72
c
cc ERF_VALUES returns some values of the ERF or "error" function for testing.
c
c  Modified:
c
c    17 April 2001
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input/output, integer N_DATA.
c    On input, if N_DATA is 0, the first test data is returned, and
c    N_DATA is set to the index of the test data.  On each subsequent
c    call, N_DATA is incremented and that test data is returned.  When
c    there is no more test data, N_DATA is set to 0.
c
c    Output, real X, the argument of the function.
c
c    Output, real FX, the value of the function.
c
      implicit none

      integer nmax
      parameter ( nmax = 21 )

      real bvec ( nmax )
      real fx
      integer n_data
      real x
      real xvec ( nmax )

      data bvec /
     &  0.0000000000E+00, 
     &  0.1124629160E+00, 
     &  0.2227025892E+00, 
     &  0.3286267595E+00, 
     &  0.4283923550E+00, 
     &  0.5204998778E+00, 
     &  0.6038560908E+00, 
     &  0.6778011938E+00, 
     &  0.7421009647E+00, 
     &  0.7969082124E+00, 
     &  0.8427007929E+00, 
     &  0.8802050696E+00, 
     &  0.9103139782E+00, 
     &  0.9340079449E+00, 
     &  0.9522851198E+00, 
     &  0.9661051465E+00, 
     &  0.9763483833E+00, 
     &  0.9837904586E+00, 
     &  0.9890905016E+00, 
     &  0.9927904292E+00, 
     &  0.9953222650E+00 /
      data xvec /
     &  0.0E+00, 
     &  0.1E+00, 
     &  0.2E+00, 
     &  0.3E+00, 
     &  0.4E+00, 
     &  0.5E+00, 
     &  0.6E+00, 
     &  0.7E+00, 
     &  0.8E+00, 
     &  0.9E+00, 
     &  1.0E+00, 
     &  1.1E+00, 
     &  1.2E+00, 
     &  1.3E+00, 
     &  1.4E+00, 
     &  1.5E+00, 
     &  1.6E+00, 
     &  1.7E+00, 
     &  1.8E+00, 
     &  1.9E+00, 
     &  2.0E+00 /

      if ( n_data < 0 ) then
        n_data = 0
      end if

      n_data = n_data + 1

      if ( nmax < n_data ) then
        n_data = 0
        x = 0.0E+00
        fx = 0.0E+00
      else
        x = xvec(n_data)
        fx = bvec(n_data)
      end if

      return
      end
      function erfc1 ( ind, x )

c*********************************************************************72
c
cc ERFC1 evaluates the complementary error function.
C
C          ERFC1(IND,X) = ERFC(X)            IF IND = 0
C          ERFC1(IND,X) = EXP(X*X)*ERFC(X)   OTHERWISE
C
c  Reference:
c
c    Armido Didonato, Alfred Morris,
c    Algorithm 708: 
c    Significant Digit Computation of the Incomplete Beta Function Ratios,
c    ACM Transactions on Mathematical Software,
c    Volume 18, Number 3, 1992, pages 360-373.
c
      real erfc1

      real a(5),b(3),p(8),q(8),r(5),s(4)
      double precision w

      data c /.564189583547756/

      data a(1) /.771058495001320e-04/, a(2)/-.133733772997339e-02/,
     *     a(3) /.323076579225834e-01/, a(4) /.479137145607681e-01/,
     *     a(5) /.128379167095513e+00/
      data b(1) /.301048631703895e-02/, b(2) /.538971687740286e-01/,
     *     b(3) /.375795757275549e+00/

      data p(1)/-1.36864857382717e-07/, p(2) /5.64195517478974e-01/,
     *     p(3) /7.21175825088309e+00/, p(4) /4.31622272220567e+01/,
     *     p(5) /1.52989285046940e+02/, p(6) /3.39320816734344e+02/,
     *     p(7) /4.51918953711873e+02/, p(8) /3.00459261020162e+02/
      data q(1) /1.00000000000000e+00/, q(2) /1.27827273196294e+01/,
     *     q(3) /7.70001529352295e+01/, q(4) /2.77585444743988e+02/,
     *     q(5) /6.38980264465631e+02/, q(6) /9.31354094850610e+02/,
     *     q(7) /7.90950925327898e+02/, q(8) /3.00459260956983e+02/

      data r(1) /2.10144126479064e+00/, r(2) /2.62370141675169e+01/,
     *     r(3) /2.13688200555087e+01/, r(4) /4.65807828718470e+00/,
     *     r(5) /2.82094791773523e-01/
      data s(1) /9.41537750555460e+01/, s(2) /1.87114811799590e+02/,
     *     s(3) /9.90191814623914e+01/, s(4) /1.80124575948747e+01/
C
C  ABS(X) .LE. 0.5
C
      AX = ABS(X)
      IF (AX .GT. 0.5) GO TO 10
      T = X*X
      TOP = ((((A(1)*T + A(2))*T + A(3))*T + A(4))*T + A(5)) + 1.0
      BOT = ((B(1)*T + B(2))*T + B(3))*T + 1.0
      ERFC1 = 0.5 + (0.5 - X*(TOP/BOT)) 
      IF (IND .NE. 0) ERFC1 = EXP(T) * ERFC1
      RETURN
C
C  0.5 .LT. ABS(X) .LE. 4
C
   10 IF (AX .GT. 4.0) GO TO 20
      TOP = ((((((P(1)*AX + P(2))*AX + P(3))*AX + P(4))*AX + P(5))*AX 
     *                    + P(6))*AX + P(7))*AX + P(8)
      BOT = ((((((Q(1)*AX + Q(2))*AX + Q(3))*AX + Q(4))*AX + Q(5))*AX 
     *                    + Q(6))*AX + Q(7))*AX + Q(8)
      ERFC1 = TOP/BOT
      GO TO 40
C
C  ABS(X) .GT. 4
C
   20 IF (X .LE. -5.6) GO TO 50
      IF (IND .NE. 0) GO TO 30
      IF (X .GT. 100.0) GO TO 60
      IF (X*X .GT. -EXPARG(1)) GO TO 60 

   30 T = (1.0/X)**2
      TOP = (((R(1)*T + R(2))*T + R(3))*T + R(4))*T + R(5)
      BOT = (((S(1)*T + S(2))*T + S(3))*T + S(4))*T + 1.0
      ERFC1 = (C - T*TOP/BOT)/AX
C
C  FINAL ASSEMBLY
C
   40 IF (IND .EQ. 0) GO TO 41
         IF (X .LT. 0.0) ERFC1 = 2.0*EXP(X*X) - ERFC1
         RETURN
   41 W = DBLE(X)*DBLE(X)
      T = W
      E = W - DBLE(T)
      ERFC1 = ((0.5 + (0.5 - E)) * EXP(-T)) * ERFC1
      IF (X .LT. 0.0) ERFC1 = 2.0 - ERFC1
      RETURN
C
C  LIMIT VALUE FOR LARGE NEGATIVE X
C
   50 ERFC1 = 2.0
      IF (IND .NE. 0) ERFC1 = 2.0*EXP(X*X)
      RETURN
C
C  LIMIT VALUE FOR LARGE POSITIVE X WHEN IND = 0
C
   60 erfc1 = 0.0

      return
      end 
      function esum ( mu, x )

c*********************************************************************72
c
cc ESUM evaluates EXP(MU+X).
C
c  Reference:
c
c    Armido Didonato, Alfred Morris,
c    Algorithm 708: 
c    Significant Digit Computation of the Incomplete Beta Function Ratios,
c    ACM Transactions on Mathematical Software,
c    Volume 18, Number 3, 1992, pages 360-373.
c
      real esum

      if (x .gt. 0.0) go to 10

      if (mu .lt. 0) go to 20 
         w = mu + x 
         if (w .gt. 0.0) go to 20
         esum = exp(w)
         return

   10 if (mu .gt. 0) go to 20 
         w = mu + x 
         if (w .lt. 0.0) go to 20
         esum = exp(w)
         return

   20 w = mu
      esum = exp(w)*exp(x)
      return
      end 
      function exparg ( l )

c*********************************************************************72
c
cc EXPARG returns the largest value for which EXP can be computed.
c
C     IF L = 0 THEN  EXPARG(L) = THE LARGEST POSITIVE W FOR WHICH
C     EXP(W) CAN BE COMPUTED. 
C
C     IF L IS NONZERO THEN  EXPARG(L) = THE LARGEST NEGATIVE W FOR
C     WHICH THE COMPUTED VALUE OF EXP(W) IS NONZERO.
C
C     NOTE... ONLY AN APPROXIMATE VALUE FOR EXPARG(L) IS NEEDED.
C
c  Reference:
c
c    Armido Didonato, Alfred Morris,
c    Algorithm 708: 
c    Significant Digit Computation of the Incomplete Beta Function Ratios,
c    ACM Transactions on Mathematical Software,
c    Volume 18, Number 3, 1992, pages 360-373.
c
      integer b
      real exparg
      real lnb

      b = ipmpar(4) 
      if (b .ne. 2) go to 10
         lnb = .69314718055995
         go to 50
   10 if (b .ne. 8) go to 20
         lnb = 2.0794415416798
         go to 50
   20 if (b .ne. 16) go to 30 
         lnb = 2.7725887222398
         go to 50
   30 lnb = alog(float(b))

   50 if (l .eq. 0) go to 60
         m = ipmpar(6) - 1
         exparg = 0.99999 * (m * lnb)
         return
   60 m = ipmpar(7) 
      exparg = 0.99999 * (m * lnb)
      return
      end 
      function fpser ( a, b, x, eps )

c*********************************************************************72
c
cc FPSER uses a series for IX(A,B) with B < min(eps,eps*A) and X <= 0.5.
c
c  Discussion:
c
C                 EVALUATION OF IX(A,B) 
C
C          FOR B .LT. MIN(EPS,EPS*A) AND X .LE. 0.5.
C
C
C                  SET  FPSER = X**A
C
c  Reference:
c
c    Armido Didonato, Alfred Morris,
c    Algorithm 708: 
c    Significant Digit Computation of the Incomplete Beta Function Ratios,
c    ACM Transactions on Mathematical Software,
c    Volume 18, Number 3, 1992, pages 360-373.
c
      real fpser

      FPSER = 1.0
      IF (A .LE. 1.E-3*EPS) GO TO 10
      FPSER = 0.0
      T = A*ALOG(X) 
      IF (T .LT. EXPARG(1)) RETURN
      FPSER = EXP(T)
C
C  NOTE THAT 1/B(A,B) = B 
C
   10 fpser = (b/a)*fpser
      tol = eps/a
      an = a + 1.0
      t = x
      s = t/an
   20    an = an + 1.0
         t = x*t
         c = t/an
         s = s + c
         if (abs(c) .gt. tol) go to 20

      fpser = fpser*(1.0 + a*s)

      return
      end 
      function gam1 ( a )

c*********************************************************************72
c
cc GAM1 evaluates 1/GAMMA(A+1) - 1  for -0.5 <= A <= 1.5
C
c  Reference:
c
c    Armido Didonato, Alfred Morris,
c    Algorithm 708: 
c    Significant Digit Computation of the Incomplete Beta Function Ratios,
c    ACM Transactions on Mathematical Software,
c    Volume 18, Number 3, 1992, pages 360-373.
c
      real gam1
      real p(7), q(5), r(9)

      data p(1)/ .577215664901533e+00/, p(2)/-.409078193005776e+00/,
     *     p(3)/-.230975380857675e+00/, p(4)/ .597275330452234e-01/,
     *     p(5)/ .766968181649490e-02/, p(6)/-.514889771323592e-02/,
     *     p(7)/ .589597428611429e-03/

      data q(1)/ .100000000000000e+01/, q(2)/ .427569613095214e+00/,
     *     q(3)/ .158451672430138e+00/, q(4)/ .261132021441447e-01/,
     *     q(5)/ .423244297896961e-02/

      data r(1)/-.422784335098468e+00/, r(2)/-.771330383816272e+00/,
     *     r(3)/-.244757765222226e+00/, r(4)/ .118378989872749e+00/,
     *     r(5)/ .930357293360349e-03/, r(6)/-.118290993445146e-01/,
     *     r(7)/ .223047661158249e-02/, r(8)/ .266505979058923e-03/,
     *     r(9)/-.132674909766242e-03/

      data s1  / .273076135303957e+00/, s2  / .559398236957378e-01/

      t = a
      d = a - 0.5
      if (d .gt. 0.0) t = d - 0.5
      if (t) 30,10,20

   10 gam1 = 0.0
      return

   20 top = (((((p(7)*t + p(6))*t + p(5))*t + p(4))*t + p(3))*t
     *                  + p(2))*t + p(1)
      bot = (((q(5)*t + q(4))*t + q(3))*t + q(2))*t + 1.0
      w = top/bot
      if (d .gt. 0.0) go to 21
         gam1 = a*w 
         return
   21 gam1 = (t/a)*((w - 0.5) - 0.5)
      return

   30 top = (((((((r(9)*t + r(8))*t + r(7))*t + r(6))*t + r(5))*t
     *                    + r(4))*t + r(3))*t + r(2))*t + r(1)
      bot = (s2*t + s1)*t + 1.0
      w = top/bot
      if (d .gt. 0.0) go to 31
         gam1 = a*((w + 0.5) + 0.5)
         return
   31 gam1 = t*w/a

      return
      end 
      function gamln ( a ) 

c*********************************************************************72
c
cc GAMLN evaluates LN(GAMMA(A)) for positive A.
C
C     WRITTEN BY ALFRED H. MORRIS
C          NAVAL SURFACE WARFARE CENTER 
C          DAHLGREN, VIRGINIA 
C
C     D = 0.5*(LN(2*PI) - 1)
C
c  Reference:
c
c    Armido Didonato, Alfred Morris,
c    Algorithm 708: 
c    Significant Digit Computation of the Incomplete Beta Function Ratios,
c    ACM Transactions on Mathematical Software,
c    Volume 18, Number 3, 1992, pages 360-373.
c
      real gamln

      data d/.418938533204673/

      data c0/.833333333333333e-01/, c1/-.277777777760991e-02/,
     *     c2/.793650666825390e-03/, c3/-.595202931351870e-03/,
     *     c4/.837308034031215e-03/, c5/-.165322962780713e-02/

      if (a .gt. 0.8) go to 10
         gamln = gamln1(a) - alog(a)
         return
   10 if (a .gt. 2.25) go to 20
         t = (a - 0.5) - 0.5
         gamln = gamln1(t)
         return

   20 if (a .ge. 10.0) go to 30
      n = a - 1.25
      t = a
      w = 1.0
      do 21 i = 1,n 
         t = t - 1.0
   21    w = t*w
      gamln = gamln1(t - 1.0) + alog(w) 
      return

   30 t = (1.0/a)**2
      w = (((((c5*t + c4)*t + c3)*t + c2)*t + c1)*t + c0)/a 
      gamln = (d + w) + (a - 0.5)*(alog(a) - 1.0) 

      return
      end 
      function gamln1 ( a )

c*********************************************************************72
c
cc GAMLN1 evaluates LN(GAMMA(1 + A)) for -0.2 <= A <= 1.25
c
c  Reference:
c
c    Armido Didonato, Alfred Morris,
c    Algorithm 708: 
c    Significant Digit Computation of the Incomplete Beta Function Ratios,
c    ACM Transactions on Mathematical Software,
c    Volume 18, Number 3, 1992, pages 360-373.
c
      real gamln1

      data p0/ .577215664901533e+00/, p1/ .844203922187225e+00/,
     *     p2/-.168860593646662e+00/, p3/-.780427615533591e+00/,
     *     p4/-.402055799310489e+00/, p5/-.673562214325671e-01/,
     *     p6/-.271935708322958e-02/
      data q1/ .288743195473681e+01/, q2/ .312755088914843e+01/,
     *     q3/ .156875193295039e+01/, q4/ .361951990101499e+00/,
     *     q5/ .325038868253937e-01/, q6/ .667465618796164e-03/
c
      data r0/.422784335098467e+00/,  r1/.848044614534529e+00/,
     *     r2/.565221050691933e+00/,  r3/.156513060486551e+00/,
     *     r4/.170502484022650e-01/,  r5/.497958207639485e-03/
      data s1/.124313399877507e+01/,  s2/.548042109832463e+00/,
     *     s3/.101552187439830e+00/,  s4/.713309612391000e-02/,
     *     s5/.116165475989616e-03/

      if (a .ge. 0.6) go to 10
      w = ((((((p6*a + p5)*a + p4)*a + p3)*a + p2)*a + p1)*a + p0)/
     *    ((((((q6*a + q5)*a + q4)*a + q3)*a + q2)*a + q1)*a + 1.0)
      gamln1 = -a*w 
      return

   10 x = (a - 0.5) - 0.5
      w = (((((r5*x + r4)*x + r3)*x + r2)*x + r1)*x + r0)/
     *    (((((s5*x + s4)*x + s3)*x + s2)*x + s1)*x + 1.0)
      gamln1 = x*w
      return
      end
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
c    28 August 2004
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
c    Output, real ( kind = 4 ) A, the parameter of the function.
c
c    Output, real ( kind = 4 ) X, the argument of the function.
c
c    Output, real ( kind = 4 ) FX, the value of the function.
c
      implicit none

      integer n_max
      parameter ( n_max = 20 )

      real a
      real a_vec ( n_max )
      real fx
      real fx_vec ( n_max )
      integer n_data
      real x
      real x_vec ( n_max )

      data a_vec /
     &   0.10E+00, 
     &   0.10E+00, 
     &   0.10E+00, 
     &   0.50E+00, 
     &   0.50E+00, 
     &   0.50E+00, 
     &   0.10E+01, 
     &   0.10E+01, 
     &   0.10E+01, 
     &   0.11E+01, 
     &   0.11E+01, 
     &   0.11E+01, 
     &   0.20E+01, 
     &   0.20E+01, 
     &   0.20E+01, 
     &   0.60E+01, 
     &   0.60E+01, 
     &   0.11E+02, 
     &   0.26E+02, 
     &   0.41E+02 /
      data fx_vec /
     &   0.7382350532339351E+00, 
     &   0.9083579897300343E+00, 
     &   0.9886559833621947E+00, 
     &   0.3014646416966613E+00, 
     &   0.7793286380801532E+00, 
     &   0.9918490284064973E+00, 
     &   0.9516258196404043E-01, 
     &   0.6321205588285577E+00, 
     &   0.9932620530009145E+00, 
     &   0.7205974576054322E-01, 
     &   0.5891809618706485E+00, 
     &   0.9915368159845525E+00, 
     &   0.01018582711118352E+00, 
     &   0.4421745996289254E+00, 
     &   0.9927049442755639E+00, 
     &   0.4202103819530612E-01, 
     &   0.9796589705830716E+00, 
     &   0.9226039842296429E+00, 
     &   0.4470785799755852E+00, 
     &   0.7444549220718699E+00 /
      data x_vec /
     &   0.30E-01, 
     &   0.30E+00, 
     &   0.15E+01, 
     &   0.75E-01, 
     &   0.75E+00, 
     &   0.35E+01, 
     &   0.10E+00, 
     &   0.10E+01, 
     &   0.50E+01, 
     &   0.10E+00, 
     &   0.10E+01, 
     &   0.50E+01, 
     &   0.15E+00, 
     &   0.15E+01, 
     &   0.70E+01, 
     &   0.25E+01, 
     &   0.12E+02, 
     &   0.16E+02, 
     &   0.25E+02, 
     &   0.45E+02 /

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
      subroutine gamma_log_values ( n_data, x, fx )

c*********************************************************************72
c
cc GAMMA_LOG_VALUES returns some values of the Log Gamma function for testing.
c
c  Modified:
c
c    17 April 2001
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
c  Parameters:
c
c    Input/output, integer N_DATA.
c    On input, if N_DATA is 0, the first test data is returned, and
c    N_DATA is set to the index of the test data.  On each subsequent
c    call, N_DATA is incremented and that test data is returned.  When
c    there is no more test data, N_DATA is set to 0.
c
c    Output, real X, the argument of the function.
c
c    Output, real FX, the value of the function.
c
      implicit none

      integer nmax
      parameter ( nmax = 18 )

      real bvec ( nmax )
      real fx
      integer n_data
      real x
      real xvec ( nmax )

      data bvec /
     &   1.524064183E+00,    0.7966780066E+00,   0.3982337117E+00,  
     &   0.1520599127E+00,   0.000000000E+00,   -0.04987246543E+00, 
     &  -0.08537410945E+00, -0.1081747934E+00,  -0.1196128950E+00,  
     &  -0.1207822040E+00,  -0.1125917658E+00,  -0.09580771625E+00, 
     &  -0.07108385116E+00, -0.03898428380E+00,  0.000000000E+00,   
     &  12.80182743E+00,    39.33988571E+00,    71.25704193E+00 /
      data xvec /
     &  0.2E+00,  0.4E+00,  0.6E+00,  0.8E+00, 
     &  1.0E+00,  1.1E+00,  1.2E+00,  1.3E+00, 
     &  1.4E+00,  1.5E+00,  1.6E+00,  1.7E+00, 
     &  1.8E+00,  1.9E+00,  2.0E+00, 10.0E+00, 
     & 20.0E+00, 30.0E+00 /

      if ( n_data < 0 ) then
        n_data = 0
      end if

      n_data = n_data + 1

      if ( nmax < n_data ) then
        n_data = 0
        x = 0.0E+00
        fx = 0.0E+00
      else
        x = xvec(n_data)
        fx = bvec(n_data)
      end if

      return
      end
      subroutine grat1 ( a, x, r, p, q, eps )

c*********************************************************************72
c
cc GRAT1 evaluates P(A,X) and Q(A,X) when A <= 1.
c
C        EVALUATION OF THE INCOMPLETE GAMMA RATIO FUNCTIONS 
C                      P(A,X) AND Q(A,X)
C
C     IT IS ASSUMED THAT A .LE. 1.  EPS IS THE TOLERANCE TO BE USED.
C     THE INPUT ARGUMENT R HAS THE VALUE E**(-X)*X**A/GAMMA(A).
C
c  Reference:
c
c    Armido Didonato, Alfred Morris,
c    Algorithm 708: 
c    Significant Digit Computation of the Incomplete Beta Function Ratios,
c    ACM Transactions on Mathematical Software,
c    Volume 18, Number 3, 1992, pages 360-373.
c
      real j
      real l

      IF (A*X .EQ. 0.0) GO TO 130
      IF (A .EQ. 0.5) GO TO 120
      IF (X .LT. 1.1) GO TO 10
      GO TO 50
C
C  TAYLOR SERIES FOR P(A,X)/X**A
C
   10 AN = 3.0
      C = X
      SUM = X/(A + 3.0)
      TOL = 0.1*EPS/(A + 1.0) 
   11    AN = AN + 1.0
         C = -C*(X/AN)
         T = C/(A + AN)
         SUM = SUM + T
         IF (ABS(T) .GT. TOL) GO TO 11
      J = A*X*((SUM/6.0 - 0.5/(A + 2.0))*X + 1.0/(A + 1.0)) 

      Z = A*ALOG(X) 
      H = GAM1(A)
      G = 1.0 + H
      IF (X .LT. 0.25) GO TO 20
         IF (A .LT. X/2.59) GO TO 40
         GO TO 30
   20 IF (Z .GT. -.13394) GO TO 40

   30 W = EXP(Z)
      P = W*G*(0.5 + (0.5 - J))
      Q = 0.5 + (0.5 - P)
      RETURN

   40 L = REXP(Z)
      W = 0.5 + (0.5 + L)
      Q = (W*J - L)*G - H
      IF (Q .LT. 0.0) GO TO 110
      P = 0.5 + (0.5 - Q)
      RETURN
C
C  CONTINUED FRACTION EXPANSION
C
   50 A2NM1 = 1.0
      A2N = 1.0
      B2NM1 = X
      B2N = X + (1.0 - A)
      C = 1.0
   51    A2NM1 = X*A2N + C*A2NM1
         B2NM1 = X*B2N + C*B2NM1
         AM0 = A2NM1/B2NM1
         C = C + 1.0
         CMA = C - A
         A2N = A2NM1 + CMA*A2N
         B2N = B2NM1 + CMA*B2N
         AN0 = A2N/B2N
         IF (ABS(AN0 - AM0) .GE. EPS*AN0) GO TO 51
      Q = R*AN0
      P = 0.5 + (0.5 - Q)
      RETURN
C
C  Special cases
c
  100 p = 0.0
      q = 1.0
      return

  110 p = 1.0
      q = 0.0
      return

  120 if (x .ge. 0.25) go to 121
      p = erf(sqrt(x))
      q = 0.5 + (0.5 - p)
      return
  121 q = erfc1(0,sqrt(x))
      p = 0.5 + (0.5 - q)
      return

  130 if (x .le. a) go to 100 
      go to 110
      end 
      function gsumln ( a, b )

c*********************************************************************72
c
cc GSUMLN evaluates LN(GAMMA(A + B)) for 1 <= A <= 2 and 1 <= B <= 2.
C
c  Reference:
c
c    Armido Didonato, Alfred Morris,
c    Algorithm 708: 
c    Significant Digit Computation of the Incomplete Beta Function Ratios,
c    ACM Transactions on Mathematical Software,
c    Volume 18, Number 3, 1992, pages 360-373.
c
      real gsumln

      x = dble(a) + dble(b) - 2.d0
      if (x .gt. 0.25) go to 10
         gsumln = gamln1(1.0 + x)
         return
   10 if (x .gt. 1.25) go to 20
         gsumln = gamln1(x) + alnrel(x) 
         return
   20 gsumln = gamln1(x - 1.0) + alog(x*(1.0 + x))

      return
      end 
      function ipmpar ( i )

c*********************************************************************72
C
cc IPMPAR sets integer machine constants.
c
C     IPMPAR PROVIDES THE INTEGER MACHINE CONSTANTS FOR THE COMPUTER
C     THAT IS USED. IT IS ASSUMED THAT THE ARGUMENT I IS AN INTEGER
C     HAVING ONE OF THE VALUES 1-10. IPMPAR(I) HAS THE VALUE ...
C
C  INTEGERS.
C
C     ASSUME INTEGERS ARE REPRESENTED IN THE N-DIGIT, BASE-A FORM
C
C               SIGN ( X(N-1)*A**(N-1) + ... + X(1)*A + X(0) )
C
C               WHERE 0 .LE. X(I) .LT. A FOR I=0,...,N-1.
C
C     IPMPAR(1) = A, THE BASE.
C
C     IPMPAR(2) = N, THE NUMBER OF BASE-A DIGITS. 
C
C     IPMPAR(3) = A**N - 1, THE LARGEST MAGNITUDE.
C
C  FLOATING-POINT NUMBERS.
C
C     IT IS ASSUMED THAT THE SINGLE AND DOUBLE PRECISION FLOATING
C     POINT ARITHMETICS HAVE THE SAME BASE, SAY B, AND THAT THE
C     NONZERO NUMBERS ARE REPRESENTED IN THE FORM 
C
C               SIGN (B**E) * (X(1)/B + ... + X(M)/B**M)
C
C               WHERE X(I) = 0,1,...,B-1 FOR I=1,...,M,
C               X(1) .GE. 1, AND EMIN .LE. E .LE. EMAX.
C
C     IPMPAR(4) = B, THE BASE.
C
C  SINGLE-PRECISION 
C
C     IPMPAR(5) = M, THE NUMBER OF BASE-B DIGITS. 
C
C     IPMPAR(6) = EMIN, THE SMALLEST EXPONENT E.
C
C     IPMPAR(7) = EMAX, THE LARGEST EXPONENT E.
C
C  DOUBLE-PRECISION 
C
C     IPMPAR(8) = M, THE NUMBER OF BASE-B DIGITS. 
C
C     IPMPAR(9) = EMIN, THE SMALLEST EXPONENT E.
C
C     IPMPAR(10) = EMAX, THE LARGEST EXPONENT E.
C
C
C
C     TO DEFINE THIS FUNCTION FOR THE COMPUTER BEING USED, ACTIVATE
C     THE DATA STATMENTS FOR THE COMPUTER BY REMOVING THE C FROM
C     COLUMN 1. (ALL THE OTHER DATA STATEMENTS SHOULD HAVE C IN
C     COLUMN 1.)
C
C
C
C     IPMPAR IS AN ADAPTATION OF THE FUNCTION I1MACH, WRITTEN BY
C     P.A. FOX, A.D. HALL, AND N.L. SCHRYER (BELL LABORATORIES).
C     IPMPAR WAS FORMED BY A.H. MORRIS (NSWC). THE CONSTANTS ARE
C     FROM BELL LABORATORIES, NSWC, AND OTHER SOURCES.
C
      INTEGER IMACH(10)
      integer ipmpar
C
C     MACHINE CONSTANTS FOR AMDAHL MACHINES.
C
C     DATA IMACH( 1) /   2 /
C     DATA IMACH( 2) /  31 /
C     DATA IMACH( 3) / 2147483647 /
C     DATA IMACH( 4) /  16 /
C     DATA IMACH( 5) /   6 /
C     DATA IMACH( 6) / -64 /
C     DATA IMACH( 7) /  63 /
C     DATA IMACH( 8) /  14 /
C     DATA IMACH( 9) / -64 /
C     DATA IMACH(10) /  63 /
C
C     MACHINE CONSTANTS FOR THE AT&T 3B SERIES, AT&T
C     PC 7300, AND AT&T 6300. 
C
C     DATA IMACH( 1) /     2 /
C     DATA IMACH( 2) /    31 /
C     DATA IMACH( 3) / 2147483647 /
C     DATA IMACH( 4) /     2 /
C     DATA IMACH( 5) /    24 /
C     DATA IMACH( 6) /  -125 /
C     DATA IMACH( 7) /   128 /
C     DATA IMACH( 8) /    53 /
C     DATA IMACH( 9) / -1021 /
C     DATA IMACH(10) /  1024 /
C
C     MACHINE CONSTANTS FOR THE BURROUGHS 1700 SYSTEM.
C
C     DATA IMACH( 1) /    2 / 
C     DATA IMACH( 2) /   33 / 
C     DATA IMACH( 3) / 8589934591 /
C     DATA IMACH( 4) /    2 / 
C     DATA IMACH( 5) /   24 / 
C     DATA IMACH( 6) / -256 / 
C     DATA IMACH( 7) /  255 / 
C     DATA IMACH( 8) /   60 / 
C     DATA IMACH( 9) / -256 / 
C     DATA IMACH(10) /  255 / 
C
C     MACHINE CONSTANTS FOR THE BURROUGHS 5700 SYSTEM.
C
C     DATA IMACH( 1) /    2 / 
C     DATA IMACH( 2) /   39 / 
C     DATA IMACH( 3) / 549755813887 /
C     DATA IMACH( 4) /    8 / 
C     DATA IMACH( 5) /   13 / 
C     DATA IMACH( 6) /  -50 / 
C     DATA IMACH( 7) /   76 / 
C     DATA IMACH( 8) /   26 / 
C     DATA IMACH( 9) /  -50 / 
C     DATA IMACH(10) /   76 / 
C
C     MACHINE CONSTANTS FOR THE BURROUGHS 6700/7700 SYSTEMS.
C
C     DATA IMACH( 1) /      2 /
C     DATA IMACH( 2) /     39 /
C     DATA IMACH( 3) / 549755813887 /
C     DATA IMACH( 4) /      8 /
C     DATA IMACH( 5) /     13 /
C     DATA IMACH( 6) /    -50 /
C     DATA IMACH( 7) /     76 /
C     DATA IMACH( 8) /     26 /
C     DATA IMACH( 9) / -32754 /
C     DATA IMACH(10) /  32780 /
C
C     MACHINE CONSTANTS FOR THE CDC 6000/7000 SERIES
C     60 BIT ARITHMETIC, AND THE CDC CYBER 995 64 BIT
C     ARITHMETIC (NOS OPERATING SYSTEM).
C
C     DATA IMACH( 1) /    2 / 
C     DATA IMACH( 2) /   48 / 
C     DATA IMACH( 3) / 281474976710655 /
C     DATA IMACH( 4) /    2 / 
C     DATA IMACH( 5) /   48 / 
C     DATA IMACH( 6) / -974 / 
C     DATA IMACH( 7) / 1070 / 
C     DATA IMACH( 8) /   95 / 
c     DATA IMACH( 9) / -926 / 
C     DATA IMACH(10) / 1070 / 
C
C     MACHINE CONSTANTS FOR THE CDC CYBER 995 64 BIT
C     ARITHMETIC (NOS/VE OPERATING SYSTEM).
C
C     DATA IMACH( 1) /     2 /
C     DATA IMACH( 2) /    63 /
C     DATA IMACH( 3) / 9223372036854775807 /
C     DATA IMACH( 4) /     2 /
C     DATA IMACH( 5) /    48 /
C     DATA IMACH( 6) / -4096 /
C     DATA IMACH( 7) /  4095 /
C     DATA IMACH( 8) /    96 /
C     DATA IMACH( 9) / -4096 /
C     DATA IMACH(10) /  4095 /
C
C     MACHINE CONSTANTS FOR THE CRAY 1, XMP, 2, AND 3.
C
C     DATA IMACH( 1) /     2 /
C     DATA IMACH( 2) /    63 /
C     DATA IMACH( 3) / 9223372036854775807 /
C     DATA IMACH( 4) /     2 /
C     DATA IMACH( 5) /    47 /
C     DATA IMACH( 6) / -8189 /
C     DATA IMACH( 7) /  8190 /
C     DATA IMACH( 8) /    94 /
C     DATA IMACH( 9) / -8099 /
C     DATA IMACH(10) /  8190 /
C
C     MACHINE CONSTANTS FOR THE DATA GENERAL ECLIPSE S/200. 
C
C     DATA IMACH( 1) /    2 / 
C     DATA IMACH( 2) /   15 / 
C     DATA IMACH( 3) / 32767 /
C     DATA IMACH( 4) /   16 / 
C     DATA IMACH( 5) /    6 / 
C     DATA IMACH( 6) /  -64 / 
C     DATA IMACH( 7) /   63 / 
C     DATA IMACH( 8) /   14 / 
C     DATA IMACH( 9) /  -64 / 
C     DATA IMACH(10) /   63 / 
C
C     MACHINE CONSTANTS FOR THE HARRIS 220.
C
C     DATA IMACH( 1) /    2 / 
C     DATA IMACH( 2) /   23 / 
C     DATA IMACH( 3) / 8388607 /
C     DATA IMACH( 4) /    2 / 
C     DATA IMACH( 5) /   23 / 
C     DATA IMACH( 6) / -127 / 
C     DATA IMACH( 7) /  127 / 
C     DATA IMACH( 8) /   38 / 
C     DATA IMACH( 9) / -127 / 
C     DATA IMACH(10) /  127 / 
C
C     MACHINE CONSTANTS FOR THE HONEYWELL 600/6000
C     AND DPS 8/70 SERIES.
C
C     DATA IMACH( 1) /    2 / 
C     DATA IMACH( 2) /   35 / 
C     DATA IMACH( 3) / 34359738367 /
C     DATA IMACH( 4) /    2 / 
C     DATA IMACH( 5) /   27 / 
C     DATA IMACH( 6) / -127 / 
C     DATA IMACH( 7) /  127 / 
C     DATA IMACH( 8) /   63 / 
C     DATA IMACH( 9) / -127 / 
C     DATA IMACH(10) /  127 / 
C
C     MACHINE CONSTANTS FOR THE HP 2100 
C     3 WORD DOUBLE PRECISION OPTION WITH FTN4
C
C     DATA IMACH( 1) /    2 / 
C     DATA IMACH( 2) /   15 / 
C     DATA IMACH( 3) / 32767 /
C     DATA IMACH( 4) /    2 / 
C     DATA IMACH( 5) /   23 / 
C     DATA IMACH( 6) / -128 / 
C     DATA IMACH( 7) /  127 / 
C     DATA IMACH( 8) /   39 / 
C     DATA IMACH( 9) / -128 / 
C     DATA IMACH(10) /  127 / 
C
C     MACHINE CONSTANTS FOR THE HP 2100 
C     4 WORD DOUBLE PRECISION OPTION WITH FTN4
C
C     DATA IMACH( 1) /    2 / 
C     DATA IMACH( 2) /   15 / 
C     DATA IMACH( 3) / 32767 /
C     DATA IMACH( 4) /    2 / 
C     DATA IMACH( 5) /   23 / 
C     DATA IMACH( 6) / -128 / 
C     DATA IMACH( 7) /  127 / 
C     DATA IMACH( 8) /   55 / 
C     DATA IMACH( 9) / -128 / 
C     DATA IMACH(10) /  127 / 
C
C     MACHINE CONSTANTS FOR THE HP 9000.
C
C     DATA IMACH( 1) /     2 /
C     DATA IMACH( 2) /    31 /
C     DATA IMACH( 3) / 2147483647 /
C     DATA IMACH( 4) /     2 /
C     DATA IMACH( 5) /    24 /
C     DATA IMACH( 6) /  -126 /
C     DATA IMACH( 7) /   128 /
C     DATA IMACH( 8) /    53 /
C     DATA IMACH( 9) / -1021 /
C     DATA IMACH(10) /  1024 /
C
C     MACHINE CONSTANTS FOR THE IBM 360/370 SERIES,
C     THE ICL 2900, THE ITEL AS/6, THE XEROX SIGMA
C     5/7/9 AND THE SEL SYSTEMS 85/86.
C
C     DATA IMACH( 1) /    2 / 
C     DATA IMACH( 2) /   31 / 
C     DATA IMACH( 3) / 2147483647 /
C     DATA IMACH( 4) /   16 / 
C     DATA IMACH( 5) /    6 / 
C     DATA IMACH( 6) /  -64 / 
C     DATA IMACH( 7) /   63 / 
C     DATA IMACH( 8) /   14 / 
C     DATA IMACH( 9) /  -64 / 
C     DATA IMACH(10) /   63 / 
C
C     MACHINE CONSTANTS FOR THE IBM PC. 
C
      DATA IMACH( 1) /     2 /
      DATA IMACH( 2) /    31 /
      DATA IMACH( 3) / 2147483647 /
      DATA IMACH( 4) /     2 /
      DATA IMACH( 5) /    24 /
      DATA IMACH( 6) /  -125 /
      DATA IMACH( 7) /   128 /
      DATA IMACH( 8) /    53 /
      DATA IMACH( 9) / -1021 /
      DATA IMACH(10) /  1024 /
C
C     MACHINE CONSTANTS FOR THE MACINTOSH II - ABSOFT
C     MACFORTRAN II.
C
C     DATA IMACH( 1) /     2 /
C     DATA IMACH( 2) /    31 /
C     DATA IMACH( 3) / 2147483647 /
C     DATA IMACH( 4) /     2 /
C     DATA IMACH( 5) /    24 /
C     DATA IMACH( 6) /  -125 /
C     DATA IMACH( 7) /   128 /
C     DATA IMACH( 8) /    53 /
C     DATA IMACH( 9) / -1021 /
C     DATA IMACH(10) /  1024 /
C
C     MACHINE CONSTANTS FOR THE MICROVAX - VMS FORTRAN.
C
C     DATA IMACH( 1) /    2 / 
C     DATA IMACH( 2) /   31 / 
C     DATA IMACH( 3) / 2147483647 /
C     DATA IMACH( 4) /    2 / 
C     DATA IMACH( 5) /   24 / 
C     DATA IMACH( 6) / -127 / 
C     DATA IMACH( 7) /  127 / 
C     DATA IMACH( 8) /   56 / 
C     DATA IMACH( 9) / -127 / 
C     DATA IMACH(10) /  127 / 
C
C     MACHINE CONSTANTS FOR THE PDP-10 (KA PROCESSOR).
C
C     DATA IMACH( 1) /    2 / 
C     DATA IMACH( 2) /   35 / 
C     DATA IMACH( 3) / 34359738367 /
C     DATA IMACH( 4) /    2 / 
C     DATA IMACH( 5) /   27 / 
C     DATA IMACH( 6) / -128 / 
C     DATA IMACH( 7) /  127 / 
C     DATA IMACH( 8) /   54 / 
C     DATA IMACH( 9) / -101 / 
C     DATA IMACH(10) /  127 / 
C
C     MACHINE CONSTANTS FOR THE PDP-10 (KI PROCESSOR).
C
C     DATA IMACH( 1) /    2 / 
C     DATA IMACH( 2) /   35 / 
C     DATA IMACH( 3) / 34359738367 /
C     DATA IMACH( 4) /    2 / 
C     DATA IMACH( 5) /   27 / 
C     DATA IMACH( 6) / -128 / 
C     DATA IMACH( 7) /  127 / 
C     DATA IMACH( 8) /   62 / 
C     DATA IMACH( 9) / -128 / 
C     DATA IMACH(10) /  127 / 
C
C     MACHINE CONSTANTS FOR THE PDP-11 FORTRAN SUPPORTING
C     32-BIT INTEGER ARITHMETIC.
C
C     DATA IMACH( 1) /    2 / 
C     DATA IMACH( 2) /   31 / 
C     DATA IMACH( 3) / 2147483647 /
C     DATA IMACH( 4) /    2 / 
C     DATA IMACH( 5) /   24 / 
C     DATA IMACH( 6) / -127 / 
C     DATA IMACH( 7) /  127 / 
C     DATA IMACH( 8) /   56 / 
C     DATA IMACH( 9) / -127 / 
C     DATA IMACH(10) /  127 / 
C
C     MACHINE CONSTANTS FOR THE SEQUENT BALANCE 8000.
C
C     DATA IMACH( 1) /     2 /
C     DATA IMACH( 2) /    31 /
C     DATA IMACH( 3) / 2147483647 /
C     DATA IMACH( 4) /     2 /
C     DATA IMACH( 5) /    24 /
C     DATA IMACH( 6) /  -125 /
C     DATA IMACH( 7) /   128 /
C     DATA IMACH( 8) /    53 /
C     DATA IMACH( 9) / -1021 /
C     DATA IMACH(10) /  1024 /
C
C     MACHINE CONSTANTS FOR THE SILICON GRAPHICS IRIS-4D
C     SERIES (MIPS R3000 PROCESSOR).
C
C     DATA IMACH( 1) /     2 /
C     DATA IMACH( 2) /    31 /
C     DATA IMACH( 3) / 2147483647 /
C     DATA IMACH( 4) /     2 /
C     DATA IMACH( 5) /    24 /
C     DATA IMACH( 6) /  -125 /
C     DATA IMACH( 7) /   128 /
C     DATA IMACH( 8) /    53 /
C     DATA IMACH( 9) / -1021 /
C     DATA IMACH(10) /  1024 /
C
C     MACHINE CONSTANTS FOR THE SUN 3.
C
C     DATA IMACH( 1) /     2 /
C     DATA IMACH( 2) /    31 /
C     DATA IMACH( 3) / 2147483647 /
C     DATA IMACH( 4) /     2 /
C     DATA IMACH( 5) /    24 /
C     DATA IMACH( 6) /  -125 /
C     DATA IMACH( 7) /   128 /
C     DATA IMACH( 8) /    53 /
C     DATA IMACH( 9) / -1021 /
C     DATA IMACH(10) /  1024 /
C
C     MACHINE CONSTANTS FOR THE UNIVAC 1100 SERIES.
C
C     DATA IMACH( 1) /    2 / 
C     DATA IMACH( 2) /   35 / 
C     DATA IMACH( 3) / 34359738367 /
C     DATA IMACH( 4) /    2 / 
C     DATA IMACH( 5) /   27 / 
C     DATA IMACH( 6) / -128 / 
C     DATA IMACH( 7) /  127 / 
C     DATA IMACH( 8) /   60 / 
C     DATA IMACH( 9) /-1024 / 
C     DATA IMACH(10) / 1023 / 
C
C     MACHINE CONSTANTS FOR THE VAX 11/780.
C
C     DATA IMACH( 1) /    2 / 
C     DATA IMACH( 2) /   31 / 
C     DATA IMACH( 3) / 2147483647 /
C     DATA IMACH( 4) /    2 / 
C     DATA IMACH( 5) /   24 / 
C     DATA IMACH( 6) / -127 / 
C     DATA IMACH( 7) /  127 / 
C     DATA IMACH( 8) /   56 / 
C     DATA IMACH( 9) / -127 / 
C     DATA IMACH(10) /  127 / 
C
      IPMPAR = IMACH(I)

      RETURN
      END 
      function psi ( xx )

c*********************************************************************72
C
cc PSI evaluates the Digamma function.
c
C
C     PSI(XX) IS ASSIGNED THE VALUE 0 WHEN THE DIGAMMA FUNCTION CANNOT
C     BE COMPUTED.
C
C     THE MAIN COMPUTATION INVOLVES EVALUATION OF RATIONAL CHEBYSHEV
C     APPROXIMATIONS PUBLISHED IN MATH. COMP. 27, 123-127(1973) BY
C     CODY, STRECOK AND THACHER.
C
C
C     PSI WAS WRITTEN AT ARGONNE NATIONAL LABORATORY FOR THE FUNPACK
C     PACKAGE OF SPECIAL FUNCTION SUBROUTINES. PSI WAS MODIFIED BY
C     A.H. MORRIS (NSWC).
C
      real psi
      REAL P1(7), P2(4), Q1(6), Q2(4)
      DOUBLE PRECISION DX0
C
C  PIOV4 = PI/4
C  DX0 = ZERO OF PSI TO EXTENDED PRECISION
C
      DATA PIOV4/.785398163397448E0/
      DATA DX0/1.461632144968362341262659542325721325D0/
C
C  COEFFICIENTS FOR RATIONAL APPROXIMATION OF
C  PSI(X) / (X - X0),  0.5 .LE. X .LE. 3.0
C
      DATA P1(1)/.895385022981970E-02/,  P1(2)/.477762828042627E+01/, 
     *     P1(3)/.142441585084029E+03/,  P1(4)/.118645200713425E+04/, 
     *     P1(5)/.363351846806499E+04/,  P1(6)/.413810161269013E+04/, 
     *     P1(7)/.130560269827897E+04/
      DATA Q1(1)/.448452573429826E+02/,  Q1(2)/.520752771467162E+03/, 
     *     Q1(3)/.221000799247830E+04/,  Q1(4)/.364127349079381E+04/, 
     *     Q1(5)/.190831076596300E+04/,  Q1(6)/.691091682714533E-05/
C
C  COEFFICIENTS FOR RATIONAL APPROXIMATION OF
C  PSI(X) - LN(X) + 1 / (2*X),  X .GT. 3.0
C
      DATA P2(1)/-.212940445131011E+01/, P2(2)/-.701677227766759E+01/,
     *     P2(3)/-.448616543918019E+01/, P2(4)/-.648157123766197E+00/ 
      DATA Q2(1)/ .322703493791143E+02/, Q2(2)/ .892920700481861E+02/,
     *     Q2(3)/ .546117738103215E+02/, Q2(4)/ .777788548522962E+01/ 
C
C  MACHINE DEPENDENT CONSTANTS ...
C
C        XMAX1  = THE SMALLEST POSITIVE FLOATING POINT CONSTANT
C                 WITH ENTIRELY INTEGER REPRESENTATION.  ALSO USED
C                 AS NEGATIVE OF LOWER BOUND ON ACCEPTABLE NEGATIVE
C                 ARGUMENTS AND AS THE POSITIVE ARGUMENT BEYOND WHICH 
C                 PSI MAY BE REPRESENTED AS ALOG(X).
C
C        XSMALL = ABSOLUTE ARGUMENT BELOW WHICH PI*COTAN(PI*X)
C                 MAY BE REPRESENTED BY 1/X.
C
C
      XMAX1 = IPMPAR(3)
      XMAX1 = AMIN1(XMAX1, 1.0/SPMPAR(1))
      XSMALL = 1.E-9

      X = XX
      AUG = 0.0E0
      IF (X .GE. 0.5E0) GO TO 200
C
C  X .LT. 0.5,  USE REFLECTION FORMULA
C  PSI(1-X) = PSI(X) + PI * COTAN(PI*X)
C
      IF (ABS(X) .GT. XSMALL) GO TO 100 
      IF (X .EQ. 0.0E0) GO TO 400
C
C  0 .LT. ABS(X) .LE. XSMALL.  USE 1/X AS A SUBSTITUTE
C  FOR  PI*COTAN(PI*X)
C
      AUG = -1.0E0 / X
      GO TO 150
C
C  REDUCTION OF ARGUMENT FOR COTAN
C
  100 W = - X
      SGN = PIOV4
      IF (W .GT. 0.0E0) GO TO 120
      W = - W
      SGN = -SGN
C
C  MAKE AN ERROR EXIT IF X .LE. -XMAX1
C
  120 IF (W .GE. XMAX1) GO TO 400
      NQ = INT(W)
      W = W - FLOAT(NQ)
      NQ = INT(W*4.0E0)
      W = 4.0E0 * (W - FLOAT(NQ) * .25E0)
C
C  W IS NOW RELATED TO THE FRACTIONAL PART OF  4.0 * X.
C  ADJUST ARGUMENT TO CORRESPOND TO VALUES IN FIRST
C  QUADRANT AND DETERMINE SIGN
C
      N = NQ / 2
      IF ((N+N) .NE. NQ) W = 1.0E0 - W
      Z = PIOV4 * W 
      M = N / 2
      IF ((M+M) .NE. N) SGN = - SGN
C
C  DETERMINE FINAL VALUE FOR  -PI*COTAN(PI*X)
C
      N = (NQ + 1) / 2
      M = N / 2
      M = M + M
      IF (M .NE. N) GO TO 140 
C
C  CHECK FOR SINGULARITY
C
      IF (Z .EQ. 0.0E0) GO TO 400
C
C  USE COS/SIN AS A SUBSTITUTE FOR COTAN, AND
C  SIN/COS AS A SUBSTITUTE FOR TAN
C
      AUG = SGN * ((COS(Z) / SIN(Z)) * 4.0E0)
      GO TO 150
  140 AUG = SGN * ((SIN(Z) / COS(Z)) * 4.0E0)
  150 X = 1.0E0 - X 
  200 IF (X .GT. 3.0E0) GO TO 300
C
C  0.5 .LE. X .LE. 3.0
C
      DEN = X
      UPPER = P1(1) * X

      DO I = 1, 5
         DEN = (DEN + Q1(I)) * X
         UPPER = (UPPER + P1(I+1)) * X
      end do

      DEN = (UPPER + P1(7)) / (DEN + Q1(6))
      XMX0 = DBLE(X) - DX0
      PSI = DEN * XMX0 + AUG
      RETURN
C
C  IF X .GE. XMAX1, PSI = LN(X)
C
  300 IF (X .GE. XMAX1) GO TO 350
C
C  3.0 .LT. X .LT. XMAX1
C
      W = 1.0E0 / (X * X)
      DEN = W
      UPPER = P2(1) * W

      DO I = 1, 3
         DEN = (DEN + Q2(I)) * W
         UPPER = (UPPER + P2(I+1)) * W
  310 end do

      AUG = UPPER / (DEN + Q2(4)) - 0.5E0 / X + AUG
  350 PSI = AUG + ALOG(X)
      RETURN
C
C  Error return
c
  400 psi = 0.0e0

      return
      end
      subroutine psi_values ( n, x, fx )

c*********************************************************************72
c
cc PSI_VALUES returns some values of the Psi or Digamma function for testing.
c
c  Discussion:
c
c    PSI(X) = d LN ( GAMMA ( X ) ) / d X = GAMMA'(X) / GAMMA(X)
c
c    PSI(1) = - Euler's constant.
c
c    PSI(X+1) = PSI(X) + 1 / X.
c
c  Modified:
c
c    17 May 2001
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
c  Parameters:
c
c    Input/output, integer N.
c    On input, if N is 0, the first test data is returned, and N is set
c    to the index of the test data.  On each subsequent call, N is
c    incremented and that test data is returned.  When there is no more
c    test data, N is set to 0.
c
c    Output, real X, the argument of the function.
c
c    Output, real FX, the value of the function.
c
      implicit none

      integer nmax
      parameter ( nmax = 11 )

      real fx
      real fxvec ( nmax )
      integer n
      real x
      real xvec ( nmax )

      data fxvec /
     &  -0.5772156649E+00,
     &  -0.4237549404E+00,
     &  -0.2890398966E+00, 
     &  -0.1691908889E+00, 
     &  -0.0613845446E+00, 
     &  -0.0364899740E+00, 
     &   0.1260474528E+00,  
     &   0.2085478749E+00,  
     &   0.2849914333E+00,   
     &   0.3561841612E+00,  
     &   0.4227843351E+00 /

      data xvec /
     &  1.0E+00,  
     &  1.1E+00,  
     &  1.2E+00,  
     &  1.3E+00,  
     &  1.4E+00,  
     &  1.5E+00,  
     &  1.6E+00,  
     &  1.7E+00,  
     &  1.8E+00,  
     &  1.9E+00,  
     &  2.0E+00 /

      if ( n < 0 ) then
        n = 0
      end if

      n = n + 1

      if ( nmax < n ) then
        n = 0
        x = 0.0E+00
        fx = 0.0E+00
      else
        x = xvec(n)
        fx = fxvec(n)
      end if

      return
      end
      function r4_epsilon ( )

c*********************************************************************72
c
cc R4_EPSILON returns the round off unit for floating arithmetic.
c
c  Discussion:
c
c    R4_EPSILON is a number R which is a power of 2 with the property that,
c    to the precision of the computer's arithmetic,
c      1 < 1 + R
c    but
c      1 = ( 1 + R / 2 )
c
c    FORTRAN90 provides the superior library routine
c
c      EPSILON ( X )
c
c  Modified:
c
c    28 November 2005
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, real R4_EPSILON, the floating point round-off unit.
c
      implicit none

      real r
      real r4_epsilon
      real r_test
    
      r = 1.0E+00
      r_test = 1.0E+00 + r / 2.0E+00

10    continue
      
      if ( 1.0E+00 < r_test ) then
        r = r / 2.0E+00
        r_test = 1.0E+00 + r / 2.0E+00
        go to 10
      end if

      r4_epsilon = r

      return
      end
      function rexp ( x )

c*********************************************************************72
c
cc REXP evaluates the function EXP(X) - 1.
C
c  Reference:
c
c    Armido Didonato, Alfred Morris,
c    Algorithm 708: 
c    Significant Digit Computation of the Incomplete Beta Function Ratios,
c    ACM Transactions on Mathematical Software,
c    Volume 18, Number 3, 1992, pages 360-373.
c
      real rexp

      data p1/ .914041914819518e-09/, p2/ .238082361044469e-01/,
     *     q1/-.499999999085958e+00/, q2/ .107141568980644e+00/,
     *     q3/-.119041179760821e-01/, q4/ .595130811860248e-03/

      if (abs(x) .gt. 0.15) go to 10
      rexp = x*(((p2*x + p1)*x + 1.0)/((((q4*x + q3)*x + q2)*x
     *                 + q1)*x + 1.0))
      return

   10 w = exp(x)
      if (x .gt. 0.0) go to 20
         rexp = (w - 0.5) - 0.5
         return
   20 rexp = w*(0.5 + (0.5 - 1.0/w))

      return
      end 
      function rlog1 ( x )

c*********************************************************************72
c
cc RLOG1 evaluates the function X - LN(1 + X).
C
c  Reference:
c
c    Armido Didonato, Alfred Morris,
c    Algorithm 708: 
c    Significant Digit Computation of the Incomplete Beta Function Ratios,
c    ACM Transactions on Mathematical Software,
c    Volume 18, Number 3, 1992, pages 360-373.
c
      real rlog1

      DATA A/.566749439387324E-01/
      DATA B/.456512608815524E-01/

      DATA P0/ .333333333333333E+00/, P1/-.224696413112536E+00/,
     *     P2/ .620886815375787E-02/
      DATA Q1/-.127408923933623E+01/, Q2/ .354508718369557E+00/

      IF (X .LT. -0.39 .OR. X .GT. 0.57) GO TO 100
      IF (X .LT. -0.18) GO TO 10
      IF (X .GT.  0.18) GO TO 20
C
C  ARGUMENT REDUCTION
C
      H = X
      W1 = 0.0
      GO TO 30

   10 H = DBLE(X) + 0.3D0
      H = H/0.7
      W1 = A - H*0.3
      GO TO 30

   20 H = 0.75D0*DBLE(X) - 0.25D0
      W1 = B + H/3.0
C
C  Series expansion
c
   30 r = h/(h + 2.0)
      t = r*r
      w = ((p2*t + p1)*t + p0)/((q2*t + q1)*t + 1.0)
      rlog1 = 2.0*t*(1.0/(1.0 - r) - r*w) + w1
      return

  100 w = (x + 0.5) + 0.5
      rlog1 = x - alog(w)
      return
      end 
      function spmpar ( i )

c*********************************************************************72
c
cc SPMPAR returns single precision real machine constants.
c
C     SPMPAR PROVIDES THE SINGLE PRECISION MACHINE CONSTANTS FOR
C     THE COMPUTER BEING USED. IT IS ASSUMED THAT THE ARGUMENT
C     I IS AN INTEGER HAVING ONE OF THE VALUES 1, 2, OR 3. IF THE
C     SINGLE PRECISION ARITHMETIC BEING USED HAS M BASE B DIGITS AND
C     ITS SMALLEST AND LARGEST EXPONENTS ARE EMIN AND EMAX, THEN
C
C        SPMPAR(1) = B**(1 - M), THE MACHINE PRECISION,
C
C        SPMPAR(2) = B**(EMIN - 1), THE SMALLEST MAGNITUDE, 
C
C        SPMPAR(3) = B**EMAX*(1 - B**(-M)), THE LARGEST MAGNITUDE.
C
C
C     WRITTEN BY
C        ALFRED H. MORRIS, JR.
C        NAVAL SURFACE WARFARE CENTER
C        DAHLGREN VIRGINIA
C
      integer emin, emax
      real spmpar

      if (i .gt. 1) go to 10
         b = ipmpar(4)
         m = ipmpar(5)
         spmpar = b**(1 - m)
         return

   10 if (i .gt. 2) go to 20
         b = ipmpar(4)
         emin = ipmpar(6)
         one = float(1)
         binv = one/b
         w = b**(emin + 2)
         spmpar = ((w * binv) * binv) * binv
         return

   20 ibeta = ipmpar(4)
      m = ipmpar(5) 
      emax = ipmpar(7)

      b = ibeta
      bm1 = ibeta - 1
      one = float(1)
      z = b**(m - 1)
      w = ((z - one)*b + bm1)/(b*z)

      z = b**(emax - 2)
      spmpar = ((w * z) * b) * b

      return
      end 
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
c    12 January 2007
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

      character * ( 8 ) ampm
      integer d
      character * ( 8 ) date
      integer h
      integer m
      integer mm
      character * ( 9 ) month(12)
      integer n
      integer s
      character * ( 10 ) time
      integer y

      save month

      data month /
     &  'January  ', 'February ', 'March    ', 'April    ', 
     &  'May      ', 'June     ', 'July     ', 'August   ', 
     &  'September', 'October  ', 'November ', 'December ' /

      call date_and_time ( date, time )

      read ( date, '(i4,i2,i2)' ) y, m, d
      read ( time, '(i2,i2,i2,1x,i3)' ) h, n, s, mm

      if ( h .lt. 12 ) then
        ampm = 'AM'
      else if ( h .eq. 12 ) then
        if ( n .eq. 0 .and. s .eq. 0 ) then
          ampm = 'Noon'
        else
          ampm = 'PM'
        end if
      else
        h = h - 12
        if ( h .lt. 12 ) then
          ampm = 'PM'
        else if ( h .eq. 12 ) then
          if ( n .eq. 0 .and. s .eq. 0 ) then
            ampm = 'Midnight'
          else
            ampm = 'AM'
          end if
        end if
      end if

      write ( *, 
     &  '(i2,1x,a,1x,i4,2x,i2,a1,i2.2,a1,i2.2,a1,i3.3,1x,a)' ) 
     &  d, month(m), y, h, ':', n, ':', s, '.', mm, ampm

      return
      end

