      subroutine bbl ( mew, theta, rl, mrl, lm, rnl )

c*********************************************************************72
C
cc BBL calculates the Beta Binomial log likelihood.
c
c  Modified:
c
c    27 January 2008
c
c  Author:
c
c    D Smith
c    Modifications by John Burkardt
c
c  Reference:
c
c    D Smith,
c    Algorithm AS 189,
c    Maximum Likelihood Estimation of the Parameters of the Beta
c    Binomial Distribution,
c    Applied Statistics, 
c    Volume 32, Number 2, 1983, pages 196-204.
c
c  Parameters:
c
c    Input, double precision MEW, the estimated value of MU.
c
c    Input, double precision THETA, the estimated value of THETA.
c
c    Input, integer RL(MRL,3), array of coefficients of 
c    (MU + R * THETA), ( 1 - MU + R * THETA) and ( 1 + R * THETA ) terms.
c
c    Input, integer MRL, the first dimension of the RL array, 
c    which must be at least the maximum of the values in IN(*).
c
c    Input, integer LM(3), contain the values Max ( IX(J) - 1 ),
c    Max ( IN(J) - IX(J) - 1 ), and Max ( IN(J) - 1 ).
c
c    Output, double precision RNL, the log likelihood. 
c
      implicit none

      integer mrl

      double precision a
      integer i
      integer lm(3)
      double precision mew
      integer mlm
      integer rl(mrl,3)
      double precision rnl
      double precision theta

      rnl = 0.0D+00
      mlm = lm(3)

      do i = 1, mlm

        a = dble ( i - 1 ) * theta

        if(i.le.lm(1)) then
          rnl = rnl + dble ( rl(i,1)) * dlog(mew+a)
        end if

        if(i.le.lm(2)) then
          rnl = rnl + dble ( rl(i,2)) * dlog(1.0D+00-mew+a)
        end if

        rnl = rnl - dble ( rl(i,3)) * dlog(1.0D+00+a)

      end do

      return
      end
      subroutine bbme ( n, ix, in, w, p, inf, mew, theta)

c*********************************************************************72
C
cc BBME estimates MEW and THETA parameters of a Beta Binomial distribution.
c
c  Discussion:
c
C        SUBROUTINE TO ESTIMATE MEW AND THETA OF THE BETA BINOMIAL
C        DISTRIBUTION BY THE METHOD OF MOMENTS
C
c  Modified:
c
c    27 January 2008
c
c  Author:
c
c    D Smith
c    Modifications by John Burkardt
c
c  Reference:
c
c    D Smith,
c    Algorithm AS 189,
c    Maximum Likelihood Estimation of the Parameters of the Beta
c    Binomial Distribution,
c    Applied Statistics, 
c    Volume 32, Number 2, 1983, pages 196-204.
c
c  Parameters:
c
c    Input, integer N, the number of observations or trials.
c
c    Input, integer IX(N), contains the number of successes for 
c    each trial.
c
c    Input, integer IN(N), contains the number tested on each trial.
c
c    Workspace, double precision P(N).
c
c    Input, double precision INF, the maximum acceptable value for THETA.
c
c    Output, double precision MEW, the estimate for MU.
c
c    Output, double precision THETA, the estimate for THETA.
c
      implicit none

      integer n

      double precision d1
      double precision d2
      integer i
      integer in(n)
      double precision inf
      integer ix(n)
      logical j
      double precision mew
      double precision p(n)
      double precision r
      double precision s
      double precision theta
      double precision tp
      double precision w(n)
      double precision wt

      j = .false.
      do i = 1, n
        w(i) = dble ( in(i) )
        p(i) = dble ( ix(i) ) / w(i)
      end do

   10 continue

      wt = 0.0D+00
      tp = 0.0D+00
      do i = 1, n
        wt = wt+w(i)
        tp = tp+w(i)*p(i)
      end do
      tp = tp/wt

      s = 0.0D+00
      d1 = 0.0D+00
      d2 = 0.0D+00
      do i = 1, n
        r = p(i) - tp
        s = s + w(i) * r * r
        r = w(i) * ( 1.0D+00-w(i) / wt )
        d1 = d1 + r/dble ( in(i) )
        d2 = d2 + r
      end do

      s = dble ( n-1)*s/dble ( n)
      r = tp*(1.0D+00-tp)
      if(r.eq.0.0D+00) goto 30
      r = (s-r*d1)/(r*(d2-d1))
      if(r.lt.0.0D+00) r = 0.0D+00
      if(j) goto 30

      do i = 1,n
        w(i) = w(i)/(1.0D+00+r*(w(i)-1.0D+00))
      end do

      j = .true.
      goto 10
   30 mew = tp
      if(r.ge.1.0D+00) goto 35
      theta = r/(1.0D+00-r)
      if(theta.le.inf) return
   35 theta = inf

      return
      end
      subroutine bbml ( n, ix, in, w, p, rl, mrl, iter, ccrit, mew, 
     &  theta, sem, seth, rnl, ifault )

c*********************************************************************72
C
CC BBML estimates parameters of a Beta Binomial distribution.
c
c  Discussion:
c     
C        SUBROUTINE FOR CALCULATING THE MAXIMUM LIKELIHOOD ESTIMATES
C        OF THE PARAMETERS OF THE BETA BINOMIAL DISTRIBUTION
C
c  Modified:
c
c    27 January 2008
c
c  Author:
c
c    D Smith
c    Modifications by John Burkardt
c
c  Reference:
c
c    D Smith,
c    Algorithm AS 189,
c    Maximum Likelihood Estimation of the Parameters of the Beta
c    Binomial Distribution,
c    Applied Statistics, 
c    Volume 32, Number 2, 1983, pages 196-204.
c
c  Parameters:
c
c    Input, integer N, the number of observations or trials.
c
c    Input, integer IX(N), contains the number of successes for 
c    each trial.
c
c    Input, integer IN(N), contains the number tested on each trial.
c
c    Workspace, double precision W(N).
c
c    Workspace, double precision P(N).
c
c    Workspace, integer RL(MRL,3), array of coefficients of 
c    (MU + R * THETA), ( 1 - MU + R * THETA) and ( 1 + R * THETA ) terms.
c
c    Input, integer MRL, the first dimension of the RL array, 
c    which must be at least the maximum of the values in IN(*).
c
c    Input/output, integer ITER;
c    On input, the maximum number of iterations allowed.
c    On output, the number of iterations taken.
c    
c    Input, double precision CCRIT, the convergence criterion.  The iteration is 
c    judged to have converged when abs ( delta MU ) and 
c    abs ( delta THETA) are less than or equal to CCRIT.
c
c    Output, double precision MEW, the maximum likelihood estimate of MU, the mean
c    of the beta binomial distribution.
c
c    Output, double precision THETA, the maximum likelihood estimate of THETA, the
c    shape parameter of the beta binomial distribution.
c
c    Output, double precision SEM, the standard error of the estimate of MU; returned
c    as -1.0 if it cannot be calculated.
c
c    Output, double precision SETH, the standard error of the estimate of THETA;
c    returned as -1.0 if it cannot be calculated.
c
c    Output, double precision RNL, the log likelihood for the maximum likelihood
c    estimates.
c
c    Output, integer IFAULT, error flag.
c    0, no error.
c    1, N <= 1;
c    2, IX(I) = 0 for all I;
c    3, IX(I) = IN(I) for all I;
c    4, max ( IN(I) ) > MRL;
c    5, either IX(I) < 0 or IN(I) < IX(I) for some I;
c    6, MU went outside the range of [0,1], or THETA went outside the
c       range [0,INF], where INF represents Infinity;
c    7, if the maximum number of iterations was exceeded;
c    8, if the damped Newton-Raphson iteration failed.
c
      implicit none

      integer mrl
      integer n

      double precision a
      double precision b
      double precision c
      double precision ccrit
      double precision d
      double precision del
      double precision dum
      double precision e
      double precision eps
      double precision f
      double precision fd(2)
      integer i
      integer ifault
      integer in(n)
      double precision inf
      parameter (inf = 1.0E+06 )
      integer iter
      integer ix(n)
      integer lm(3)
      logical mc
      double precision mew
      integer nnd
      double precision p(n)
      integer rd1(2,2)
      integer rd2(2,3)
      integer rd3(2,4)
      integer rl(mrl,3)
      double precision rnl
      double precision sd(3)
      double precision sem
      double precision seth
      double precision td(4)
      double precision theta
      double precision ub(2)
      double precision w(n)

      data
     *  rd1(1,1), rd1(2,1), rd1(1,2), rd1(2,2)/1,-1,1,1/,
     *  rd2(1,1), rd2(2,1), rd2(1,2), rd2(2,2), 
     *  rd2(1,3), rd2(2,3)/-1,-1,-1,1,-1,-1/,
     *  rd3(1,1), rd3(2,1), rd3(1,2), rd3(2,2), rd3(1,3), 
     *  rd3(2,3), rd3(1,4), rd3(2,4)/2,-2,2,2,2,-2,2,2/

      i = iter
      iter = 0
      mc = .true.
      ub(1) = 0.01D+00
      ub(2) = 0.01D+00
c
c  set the arrays rl and lm
c
      call set ( n, ix, in, rl, mrl, lm, ifault )

      if(ifault.ne.0) then
        return
      end if

      sem = -1.0D+00
      seth = -1.0D+00
      nnd = 0
c
c  calculation of initial estimates (by moments)
c
      call bbme ( n, ix, in, w, p, inf, mew, theta )
      if(theta.eq.inf) goto 50
c
c  newton-raphson iteration on first derivatives
c
    5 if(iter.le.i) goto 10
      ifault = 7
      goto 60
c
c  calculate first derivatives of log likelihood
c
   10 continue

      call gder(mew, theta, rl, mrl, lm, 2, rd1, fd)
c
c  calculate second derivatives of log_likelihood
c
      call gder(mew, theta, rl, mrl, lm, 3, rd2, sd)
c
c  calculate third derivatives of log likelihood
c
      call gder(mew, theta, rl, mrl, lm, 4, rd3, td)
c
c  calculate increments
c
      dum = sd(1)*sd(3) - sd(2)*sd(2)

      if(sd(1).lt.0.0D+00 .and. dum .gt. 0.0D+00 ) goto 15
c
c  non negative definite matrix
c
      nnd = nnd+1
c
c  sd(1) is always negative so a gradient step is made on mew
c
      a = mew - fd(1)/sd(1)
      b = theta
      if(fd(2).ne.0.0D+00) b = b + dsign(ub(2),fd(2))
      if(a.le.0.0D+00) a = 0.0001D+00
      if(a.ge.1.0D+00) a = 0.9999D+00
      if(b.lt.0.0D+00) b = 0.0D+00
      if(b.gt.inf) b = inf
      call bbl(mew, theta, rl, mrl, lm, c)
      call bbl(a, b, rl, mrl, lm, d)
      if(nnd.gt.10.or.c.ge.d) goto 40
      iter = iter+1
      mew = a
      theta = b
      goto 5
   15 del = (fd(2)*sd(2) - fd(1)*sd(3))/dum
      eps = (fd(1)*sd(2) - fd(2)*sd(1))/dum
c
c  check lipschitz condition satisfied
c
      a = sd(2)*td(2) - td(1)*sd(3)
      b = sd(2)*td(3) - td(2)*sd(3)
      c = td(1)*sd(2) - td(2)*sd(1)
      d = sd(2)*td(2) - sd(1)*td(3)
      e = sd(2)*td(4) - td(3)*sd(3)
      f = td(3)*sd(2) - td(4)*sd(1)
      a = del*a + eps*b
      c = del*c + eps*d
      e = del*b + eps*e
      f = del*d + eps*f
      dum = (a*a + c*c + e*e + f*f)/(dum*dum)
      if(dum.ge.1.0) goto 20
      if( dabs(del).le.ccrit.and. dabs(eps).le.ccrit) mc = .false.
      goto 45
c
c  Failure of lipschitz condition. a step in the direction of the
c  gradient is made.
c
   20 a = fd(1)*fd(1)
      b = fd(2)*fd(2)
      c = a*sd(1) + 2.0D+00 * sd(2)*fd(1)*fd(2) + b*sd(3)
      if(c.ne.0.0D+00) goto 25
      del = 0.0D+00
      if(fd(1).ne.0.0D+00) del = dsign(ub(1),fd(1))
      eps = 0.0D+00
      if(fd(2).ne.0.0D+00) eps = dsign(ub(2),fd(2))
      goto 30
   25 c = -(a+b)/c
      del = c*fd(1)
      eps = c*fd(2)
      if( dabs(del).gt.ub(1)) del = dsign(ub(1),del)
      ub(1) = 2.0D+00 * dabs(del)
      if ( dabs(eps).gt.ub(2)) eps = dsign(ub(2),eps)
      ub(2) = 2.0D+00 * dabs(eps)
   30 call bbl(mew, theta, rl, mrl, lm, c)
   35 a = mew + del
      b = theta + eps
      if(a.le.0.0D+00) a = 0.0001D+00
      if(a.ge.1.0D+00) a = 0.9999D+00
      del = a - mew
      if(b.lt.0.0D+00) b = 0.0D+00
      if(b.gt.inf) b = inf
      eps = b - theta
      call bbl(a, b, rl, mrl, lm, d)
c
c  check to see if gradient step has increased log likelihood
c
      if(d.gt.c) goto 45
      del = del/2.0D+00
      eps = eps/2.0D+00
      if( dabs(del).gt.ccrit.or. dabs(eps).gt.ccrit) goto 35
   40 ifault = 8
      goto 60
   45 iter = iter + 1
      a = mew + del
      b = theta + eps
      if ( a.gt.0.0D+00 .and. 
     &     a.lt.1.0D+00 .and.
     &     b.ge.0.0D+00 .and. 
     &     b.le.inf) goto 55
      if(a.le.0.0D+00) mew = 0.0D+00
      if(a.ge.1.0D+00) mew = 1.0D+00
      if(b.lt.0.0D+00) theta = 0.0D+00
      if(b.gt.inf) theta = inf
   50 ifault = 6
      goto 60
   55 mew = a
      theta = b
      if(mc) goto 5
c
c  calculate log likelihood and s.e.s
c
      if(sd(1).lt.0.0D+00) sem = dsqrt(-1.0D+00/sd(1))
      if(sd(3).lt.0.0D+00) seth = dsqrt(-1.0D+00/sd(3))
   60 call bbl(mew, theta, rl, mrl, lm, rnl)

      return
      end
      function beta ( a, b )

c*********************************************************************72
c
cc BETA returns the value of the Beta function.
c
c  Formula:
c
c    BETA(A,B) = ( GAMMA ( A ) * GAMMA ( B ) ) / GAMMA ( A + B )
c              = Integral ( 0 <= T <= 1 ) T**(A-1) (1-T)**(B-1) dT.
c
c  Modified:
c
c    10 July 1998
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, double precision A, B, the parameters of the function.
c    0.0 < A,
c    0.0 < B.
c
c    Output, double precision BETA, the value of the function.
c
      implicit none

      double precision a
      double precision b
      double precision beta
      double precision dlgama
c
c  Check.
c
      if ( a <= 0.0D+00 .or. b <= 0.0D+00 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'BETA - Fatal error!'
        write ( *, '(a)' ) '  Both A and B must be greater than 0.'
        stop
      end if

      beta = exp ( dlgama ( a ) 
     &           + dlgama ( b ) 
     &           - dlgama ( a + b ) )

      return
      end
      subroutine beta_binomial_cdf_inv ( cdf, a, b, c, x )

c*********************************************************************72
c
cc BETA_BINOMIAL_CDF_INV inverts the Beta Binomial CDF.
c
c  Discussion:
c
c    A simple-minded discrete approach is used.
c
c  Modified:
c
c    07 December 1999
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, double precision CDF, the value of the CDF.
c
c    Input, double precision A, B, parameters of the PDF.
c    0.0 < A,
c    0.0 < B.
c
c    Input, integer C, a parameter of the PDF.
c    0 <= C.
c
c    Output, integer X, the smallest X whose cumulative density function
c    is greater than or equal to CDF.
c
      implicit none

      double precision a
      double precision b
      double precision beta
      integer c
      double precision cdf
      double precision cum
      double precision pdf
      integer x
      integer y

      if ( cdf .le. 0.0D+00 ) then

        x = 0

      else

        cum = 0.0D+00

        do y = 0, c

          pdf = beta ( a + dble ( y ), b + dble ( c - y ) ) 
     &    / ( dble ( c + 1 ) 
     &    * beta ( dble ( y + 1), dble ( c - y + 1 ) ) * beta ( a, b ) )

          cum = cum + pdf

          if ( cdf .le. cum ) then
            x = y
            return
          end if

        end do

        x = c

      end if

      return
      end
      subroutine beta_binomial_check ( a, b, c )

c*********************************************************************72
c
cc BETA_BINOMIAL_CHECK checks the parameters of the Beta Binomial PDF.
c
c  Modified:
c
c    07 December 1999
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, double precision A, B, parameters of the PDF.
c    0.0 < A,
c    0.0 < B.
c
c    Input, integer C, a parameter of the PDF.
c    0 <= C.
c
      implicit none

      double precision a
      double precision b
      integer c

      if ( a .le. 0.0D+00 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'BETA_BINOMIAL_CHECK - Fatal error!'
        write ( *, '(a)' ) '  A <= 0.'
        stop
      end if

      if ( b .le. 0.0D+00 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'BETA_BINOMIAL_CHECK - Fatal error!'
        write ( *, '(a)' ) '  B <= 0.'
        stop
      end if

      if ( c .lt. 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'BETA_BINOMIAL_CHECK - Fatal error!'
        write ( *, '(a)' ) '  C < 0.'
        stop
      end if

      return
      end
      subroutine beta_binomial_sample ( a, b, c, seed, x )

c*********************************************************************72
c
c! BETA_BINOMIAL_SAMPLE samples the Beta Binomial CDF.
c
c  Modified:
c
c    07 December 1999
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, double precision A, B, parameters of the PDF.
c    0.0 < A,
c    0.0 < B.
c
c    Input, integer C, a parameter of the PDF.
c    0 <= C.
c
c    Input/output, integer SEED, a seed for the random number generator.
c
c    Output, integer X, a sample of the PDF.
c
      implicit none

      double precision a
      double precision b
      integer c
      double precision cdf
      double precision r8_uniform_01
      integer seed
      integer x

      cdf = r8_uniform_01 ( seed )

      call beta_binomial_cdf_inv ( cdf, a, b, c, x )

      return
      end
      function dlgama ( x )

c*********************************************************************72
c
cc DLGAMA evaluates log ( Gamma ( X ) ) for a real argument.
c
c  Discussion:
c
c    This routine calculates the LOG(GAMMA) function for a positive real
c    argument X.  Computation is based on an algorithm outlined in
c    references 1 and 2.  The program uses rational functions that
c    theoretically approximate LOG(GAMMA) to at least 18 significant
c    decimal digits.  The approximation for X > 12 is from reference
c    3, while approximations for X < 12.0 are similar to those in
c    reference 1, but are unpublished.
c
c  Modified:
c
c    03 April 2007
c
c  Author:
c
c    William Cody, Laura Stoltz
c
c  Reference:
c
c    William Cody, Kenneth Hillstrom,
c    Chebyshev Approximations for the Natural Logarithm of the 
c    Gamma Function,
c    Mathematics of Computation,
c    Volume 21, Number 98, April 1967, pages 198-203.
c
c    Kenneth Hillstrom,
c    ANL/AMD Program ANLC366S, DGAMMA/DLGAMA,
c    May 1969.
c
c    John Hart, Ward Cheney, Charles Lawson, Hans Maehly, 
c    Charles Mesztenyi, John Rice, Henry Thatcher, 
c    Christoph Witzgall,
c    Computer Approximations,
c    Wiley, 1968,
c    LC: QA297.C64.
c
c  Parameters:
c
c    Input, double precision X, the argument of the function.
c
c    Output, double precision DLGAMA, the value of the function.
c
      implicit none

      double precision c(7)
      double precision corr
      double precision d1
      double precision d2
      double precision d4
      double precision dlgama
      double precision eps
      double precision frtbig
      double precision four
      double precision half
      integer i
      double precision one
      double precision pnt68
      double precision p1(8)
      double precision p2(8)
      double precision p4(8)
      double precision q1(8)
      double precision q2(8)
      double precision q4(8)
      double precision res
      double precision sqrtpi
      double precision thrhal
      double precision twelve
      double precision two
      double precision x
      double precision xbig
      double precision xden
      double precision xinf
      double precision xm1
      double precision xm2
      double precision xm4
      double precision xnum
      double precision y
      double precision ysq
      double precision zero
c
c  Mathematical constants
c
      data one /1.0D+00/
      data half /0.5D+00/
      data twelve /12.0D+00/
      data zero /0.0D+00/
      data four /4.0D+00/
      data thrhal /1.5D+00/
      data two /2.0D+00/
      data pnt68 /0.6796875D+00/
      data sqrtpi /0.9189385332046727417803297D+00/
c
c  Machine dependent parameters
c
      data xbig /2.55D+305/
      data xinf /1.79D+308/
      data eps /2.22D-16/
      data frtbig /2.25D+76/
c
c  Numerator and denominator coefficients for rational minimax
c  approximation over (0.5,1.5).
c
      data d1/-5.772156649015328605195174D-01/
      data p1/
     &   4.945235359296727046734888D+00,
     &   2.018112620856775083915565D+02,
     &   2.290838373831346393026739D+03,
     &   1.131967205903380828685045D+04,
     &   2.855724635671635335736389D+04,
     &   3.848496228443793359990269D+04,
     &   2.637748787624195437963534D+04,
     &   7.225813979700288197698961D+03/
      data q1/
     &   6.748212550303777196073036D+01,
     &   1.113332393857199323513008D+03,
     &   7.738757056935398733233834D+03,
     &   2.763987074403340708898585D+04,
     &   5.499310206226157329794414D+04,
     &   6.161122180066002127833352D+04,
     &   3.635127591501940507276287D+04,
     &   8.785536302431013170870835D+03/
c
c  Numerator and denominator coefficients for rational minimax
c  Approximation over (1.5,4.0).
c
      data d2/4.227843350984671393993777D-01/
      data p2/
     &   4.974607845568932035012064D+00,
     &   5.424138599891070494101986D+02,
     &   1.550693864978364947665077D+04,
     &   1.847932904445632425417223D+05,
     &   1.088204769468828767498470D+06,
     &   3.338152967987029735917223D+06,
     &   5.106661678927352456275255D+06,
     &   3.074109054850539556250927D+06/
      data q2/
     &   1.830328399370592604055942D+02,
     &   7.765049321445005871323047D+03,
     &   1.331903827966074194402448D+05,
     &   1.136705821321969608938755D+06,
     &   5.267964117437946917577538D+06,
     &   1.346701454311101692290052D+07,
     &   1.782736530353274213975932D+07,
     &   9.533095591844353613395747D+06/
c
c  Numerator and denominator coefficients for rational minimax
c  Approximation over (4.0,12.0).
c
      data d4/1.791759469228055000094023D+00/
      data p4/
     &   1.474502166059939948905062D+04,
     &   2.426813369486704502836312D+06,
     &   1.214755574045093227939592D+08,
     &   2.663432449630976949898078D+09,
     &   2.940378956634553899906876D+10,
     &   1.702665737765398868392998D+11,
     &   4.926125793377430887588120D+11,
     &   5.606251856223951465078242D+11/
      data q4/
     &   2.690530175870899333379843D+03,
     &   6.393885654300092398984238D+05,
     &   4.135599930241388052042842D+07,
     &   1.120872109616147941376570D+09,
     &   1.488613728678813811542398D+10,
     &   1.016803586272438228077304D+11,
     &   3.417476345507377132798597D+11,
     &   4.463158187419713286462081D+11/
c
c  Coefficients for minimax approximation over (12, INF).
c
      data c/
     &  -1.910444077728D-03,
     &   8.4171387781295D-04,
     &  -5.952379913043012D-04,
     &   7.93650793500350248D-04,
     &  -2.777777777777681622553D-03,
     &   8.333333333333333331554247D-02,
     &   5.7083835261D-03/

      y = x

      if ( zero .lt. y .and. y .le. xbig ) then

        if ( y .le. eps ) then

          res = - dlog ( y )
c
c  EPS < X <= 1.5.
c
        else if ( y .le. thrhal ) then

          if ( y .lt. pnt68 ) then
            corr = - dlog ( y )
            xm1 = y
          else
            corr = zero
            xm1 = ( y - half ) - half
          end if

          if ( y .le. half .or. pnt68 .le. y ) then

            xden = one
            xnum = zero
            do i = 1, 8
              xnum = xnum * xm1 + p1(i)
              xden = xden * xm1 + q1(i)
            end do

            res = corr + ( xm1 * ( d1 + xm1 * ( xnum / xden ) ) )

          else

            xm2 = ( y - half ) - half
            xden = one
            xnum = zero
            do i = 1, 8
              xnum = xnum * xm2 + p2(i)
              xden = xden * xm2 + q2(i)
            end do

            res = corr + xm2 * ( d2 + xm2 * ( xnum / xden ) )

          end if
c
c  1.5 < X <= 4.0.
c
        else if ( y .le. four ) then

          xm2 = y - two
          xden = one
          xnum = zero
          do i = 1, 8
            xnum = xnum * xm2 + p2(i)
            xden = xden * xm2 + q2(i)
          end do

          res = xm2 * ( d2 + xm2 * ( xnum / xden ) )
c
c  4.0 < X <= 12.0.
c
        else if ( y .le. twelve ) then

          xm4 = y - four
          xden = -one
          xnum = zero
          do i = 1, 8
            xnum = xnum * xm4 + p4(i)
            xden = xden * xm4 + q4(i)
          end do

          res = d4 + xm4 * ( xnum / xden )
c
c  Evaluate for 12 <= argument.
c
        else

          res = zero

          if ( y .le. frtbig ) then

            res = c(7)
            ysq = y * y

            do i = 1, 6
              res = res / ysq + c(i)
            end do

          end if

          res = res / y
          corr = dlog ( y )
          res = res + sqrtpi - half * corr
          res = res + y * ( corr - one )

        end if
c
c  Return for bad arguments.
c
      else

        res = xinf

      end if
c
c  Final adjustments and return.
c
      dlgama = res

      return
      end
      subroutine gder ( mew, theta, rl, mrl, lm, ider, rd, pd )

c*********************************************************************72
C
cc GDER is the general derivative subroutine.
c
c  Modified:
c
c    27 January 2008
c
c  Author:
c
c    D Smith
c    Modifications by John Burkardt
c
c  Reference:
c
c    D Smith,
c    Algorithm AS 189,
c    Maximum Likelihood Estimation of the Parameters of the Beta
c    Binomial Distribution,
c    Applied Statistics, 
c    Volume 32, Number 2, 1983, pages 196-204.
c
c  Parameters:
c
c    Input, double precision MEW, the estimated value of MU.
c
c    Input, double precision THETA, the estimated value of THETA.
c
c    Input, integer RL(MRL,3), array of coefficients of 
c    (MU + R * THETA), ( 1 - MU + R * THETA) and ( 1 + R * THETA ) terms.
c
c    Input, integer MRL, the first dimension of the RL array, 
c    which must be at least the maximum of the values in IN(*).
c
c    Input, integer LM(3), contain the values Max ( IX(J) - 1 ),
c    Max ( IN(J) - IX(J) - 1 ), and Max ( IN(J) - 1 ).
c
c    Input, integer IDER, 1 plus the order of the derivative
c    desired.  IDER can be 2, 3 or 4.
c
c    Input, integer RD(2,IDER), an array of coefficients.
c
c    Output, double precision PD(IDER), the derivatives of the log likelihood. 
c
      implicit none

      integer ider
      integer mrl

      double precision a
      double precision b
      double precision c
      double precision d
      integer i
      integer j
      integer k
      integer kk
      integer lm(3)
      double precision mew
      integer mlm
      double precision pd(ider)
      integer rd(2,ider)
      integer rl(mrl,3)
      double precision theta

      mlm = lm(3)
      kk = ider-1

      do i = 1,ider
        pd(i) = 0.0D+00
      end do

      do i = 1,mlm
        c = dble ( i-1)
        a = c*theta
        do 40 j = 1,3
          if(i.gt.lm(j)) goto 40
          goto (10,15,20) j
   10     d = mew+a
          goto 25
   15     d = 1.0D+00 - mew+a
          goto 25
   20     d = 1.0+a
   25     b = dble ( rl(i,j))/d**kk
          if(j.eq.3) goto 35
          do k = 1,ider
            pd(k) = pd(k)+dble ( rd(j,k))*b
            b = b*c
          end do
          goto 40
   35     d = -dble ( rd(1,1))*b*c**kk
          pd(ider) = pd(ider)+d
   40   continue

      end do

      return
      end
      function r8_uniform_01 ( seed )

c*********************************************************************72
c
cc R8_UNIFORM_01 returns a unit pseudorandom R8.
c
c  Discussion:
c
c    This routine implements the recursion
c
c      seed = 16807 * seed mod ( 2**31 - 1 )
c      r8_uniform_01 = seed / ( 2**31 - 1 )
c
c    The integer arithmetic never requires more than 32 bits,
c    including a sign bit.
c
c    If the initial seed is 12345, then the first three computations are
c
c      Input     Output      R8_UNIFORM_01
c      SEED      SEED
c
c         12345   207482415  0.096616
c     207482415  1790989824  0.833995
c    1790989824  2035175616  0.947702
c
c  Modified:
c
c    11 August 2004
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Paul Bratley, Bennett Fox, Linus Schrage,
c    A Guide to Simulation,
c    Second Edition,
c    Springer, 1987,
c    ISBN: 0387964673,
c    LC: QA76.9.C65.B73.
c
c    Bennett Fox,
c    Algorithm 647:
c    Implementation and Relative Efficiency of Quasirandom
c    Sequence Generators,
c    ACM Transactions on Mathematical Software,
c    Volume 12, Number 4, December 1986, pages 362-376.
c
c    Pierre L'Ecuyer,
c    Random Number Generation,
c    in Handbook of Simulation,
c    edited by Jerry Banks,
c    Wiley, 1998,
c    ISBN: 0471134031,
c    LC: T57.62.H37.
c
c    Peter Lewis, Allen Goodman, James Miller,
c    A Pseudo-Random Number Generator for the System/360,
c    IBM Systems Journal,
c    Volume 8, Number 2, 1969, pages 136-143.
c
c  Parameters:
c
c    Input/output, integer SEED, the "seed" value, which should NOT be 0.
c    On output, SEED has been updated.
c
c    Output, double precision R8_UNIFORM_01, a new pseudorandom variate,
c    strictly between 0 and 1.
c
      implicit none

      integer i4_huge
      parameter ( i4_huge = 2147483647 )
      integer k
      double precision r8_uniform_01
      integer seed

      if ( seed .eq. 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'R8_UNIFORM_01 - Fatal error!'
        write ( *, '(a)' ) '  Input value of SEED = 0.'
        stop
      end if

      k = seed / 127773

      seed = 16807 * ( seed - k * 127773 ) - k * 2836

      if ( seed .lt. 0 ) then
        seed = seed + i4_huge
      end if
c
c  Although SEED can be represented exactly as a 32 bit integer,
c  it generally cannot be represented exactly as a 32 bit real number!
c
      r8_uniform_01 = dble ( seed ) * 4.656612875D-10

      return
      end
      subroutine set ( n, ix, in, rl, mrl, lm, ifault )

c*********************************************************************72
C
cc SET sets up an array for the Beta Binomial log likelihood calculation.
c
c  Discussion:
c
C        SUBROUTINE FOR SETTING UP ARRAY FOR CALCULATION OF
C        THE BETA BINOMIAL LOG LIKELIHOOD AND ITS DERIVATIVES
C
c  Modified:
c
c    27 January 2008
c
c  Author:
c
c    D Smith
c    Modifications by John Burkardt
c
c  Reference:
c
c    D Smith,
c    Algorithm AS 189,
c    Maximum Likelihood Estimation of the Parameters of the Beta
c    Binomial Distribution,
c    Applied Statistics, 
c    Volume 32, Number 2, 1983, pages 196-204.
c
c  Parameters:
c
c    Input, integer N, the number of observations or trials.
c
c    Input, integer IX(N), contains the number of successes for 
c    each trial.
c
c    Input, integer IN(N), contains the number tested on each trial.
c
c    Output, integer RL(MRL,3), array of coefficients of 
c    (MU + R * THETA), ( 1 - MU + R * THETA) and ( 1 + R * THETA ) terms.
c
c    Input, integer MRL, the first dimension of the RL array, 
c    which must be at least the maximum of the values in IN(*).
c
c    Output, integer LM(3), contain the values Max ( IX(J) - 1 ),
c    Max ( IN(J) - IX(J) - 1 ), and Max ( IN(J) - 1 ).
c
c    Output, integer IFAULT, error flag.
c    0, no error.
c    1, N <= 1;
c    2, IX(I) = 0 for all I;
c    3, IX(I) = IN(I) for all I;
c    4, max ( IN(I) ) > MRL;
c    5, either IX(I) < 0 or IN(I) < IX(I) for some I.
c
      implicit none

      integer mrl
      integer n

      integer i
      integer ifault
      integer in(n)
      integer ix(n)
      integer j
      integer jj
      integer k
      integer lm(3)
      integer mar
      integer rl(mrl,3)
c
c     test admissibility of data
c
      if(n.gt.1) goto 5
      ifault = 1
      return
    5 continue

      do i = 1,n
        if(ix(i).gt.0) goto 15
      end do

      ifault = 2
      return
   15 continue

      do i = 1,n
        if(ix(i).lt.in(i)) goto 25
      end do

      ifault = 3
      return
c
c  Form matrix of counts.
c
   25 continue

      ifault = 4

      do i = 1,3
        lm(i) = 0
        do j = 1,mrl
          rl(j,i) = 0
        end do
      end do

      do 65 i = 1,n
        jj = ix(i)
        mar = 1
        goto 45
   35   jj = in(i)-ix(i)
        mar = 2
        goto 45
   40   jj = in(i)
        mar = 3
   45   if(jj) 50,60,55
   50   ifault = 5
        return
   55   if(jj.gt.mrl) return
        if(jj.gt.lm(mar)) lm(mar) = jj
        rl(jj,mar) = rl(jj,mar)+1
   60   goto(35,40,65) mar
   65 continue
      ifault = 0
c
c  Evaluate number of calls to different terms of likelihood function.
c
      do i = 1,3
        jj = lm(i)-1
        k = jj
        do j = 1,jj
          rl(k,i) = rl(k,i)+rl(k+1,i)
          k = k-1
        end do
      end do

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

