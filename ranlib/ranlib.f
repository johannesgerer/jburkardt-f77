      subroutine advnst ( k )

c*********************************************************************72
c
cc ADVNST advances the state of the current generator.
c
c  Discussion:
c
c    This routine advances the state  of the current  generator by 2^K values and
c    resets the initial seed to that value.
c
c    This is a transcription from Pascal to Fortran of routine "Advance_State".
c
c  Modified:
c
c    10 December 2007
c
c  Reference:
c
c    Pierre LEcuyer, Serge Cote,
c    Implementing a Random Number Package with Splitting Facilities,
c    ACM Transactions on Mathematical Software,
c    Volume 17, 1991, pages 98-111.
c
c  Parmeters:
c
c    Input, integer K, indicates that the generator is to be advanced by
c    2^K values.
c
      integer numg
      parameter (numg=32)
c     ..
c     .. Scalar Arguments ..
      INTEGER k
c     ..
c     .. Scalars in Common ..
      INTEGER a1,a1vw,a1w,a2,a2vw,a2w,m1,m2
c     ..
c     .. Arrays in Common ..
      INTEGER cg1(numg),cg2(numg),ig1(numg),ig2(numg),lg1(numg),
     +        lg2(numg)
      LOGICAL qanti(numg)
c     ..
c     .. Local Scalars ..
      INTEGER g,i,ib1,ib2
c     ..
c     .. External Functions ..
      INTEGER mltmod
      LOGICAL qrgnin
      EXTERNAL mltmod,qrgnin
c     ..
c     .. External Subroutines ..
      EXTERNAL getcgn,setsd
c     ..
c     .. Common blocks ..
      COMMON /globe/m1,m2,a1,a2,a1w,a2w,a1vw,a2vw,ig1,ig2,lg1,lg2,cg1,
     +       cg2,qanti
c     ..
c     .. Save statement ..
      SAVE /globe/
c
c  Abort unless random number generator initialized
c
      IF ( .not. qrgnin() ) then
        WRITE (*,*) ' ADVNST called before random number generator ',
     +  ' initialized -- abort!'
        STOP ' ADVNST called before random number generator initialized'
      end if
c
c  Get the index of the current random number generator.
c
      call getcgn ( g )

      ib1 = a1
      ib2 = a2

      DO i = 1,k
        ib1 = mltmod ( ib1, ib1, m1 )
        ib2 = mltmod ( ib2, ib2, m2 )
      end do

      CALL setsd ( mltmod(ib1,cg1(g),m1), mltmod(ib2,cg2(g),m2) )
c
c  NOW IB1 = A1**K AND IB2 = A2**K
c
      return
      end
      REAL FUNCTION genbet(aa,bb)

c*********************************************************************72
c
cc GENBET generates a BETA random deviate.
c
c  Discussion:
c
c     Returns a single random deviate from the beta distribution with
c     parameters A and B.  The density of the beta is
c               x^(a-1) * (1-x)^(b-1) / B(a,b) for 0 < x < 1
c
c  Modified:
c
c    10 December 2007
c
c  Reference:
c
c    Russell Cheng,
c    Generating Beta Variates with Nonintegral Shape Parameters,
c    Communications of the ACM,
c    Volume 21, 1978, pages 317-322.
c
c  Parameters:
c
c     A --> First parameter of the beta distribution
c                         REAL A
c
c     B --> Second parameter of the beta distribution
c                         REAL B
c
c

c
c     Close to the largest number that can be exponentiated
c
      REAL expmax
      PARAMETER (expmax=89.0)
c     Close to the largest representable single precision number
      REAL infnty
      PARAMETER (infnty=1.0E38)
c     ..
c     .. Scalar Arguments ..
      REAL aa,bb
c     ..
c     .. Local Scalars ..
      REAL a,alpha,b,beta,delta,gamma,k1,k2,olda,oldb,r,s,t,u1,u2,v,w,y,
     +     z
      LOGICAL qsame
c     ..
c     .. External Functions ..
      REAL ranf
      EXTERNAL ranf
c     ..
c     .. Intrinsic Functions ..
      INTRINSIC exp,log,max,min,sqrt
c     ..
c     .. Save statement ..
      SAVE olda,oldb,alpha,beta,gamma,k1,k2
c     ..
c     .. Data statements ..
      DATA olda,oldb/-1,-1/
c     ..
c     .. Executable Statements ..
      qsame = (olda.EQ.aa) .AND. (oldb.EQ.bb)
      IF (qsame) GO TO 20
      IF (.NOT. (aa.LE.0.0.OR.bb.LE.0.0)) GO TO 10
      WRITE (*,*) ' AA or BB <= 0 in GENBET - Abort!'
      WRITE (*,*) ' AA: ',aa,' BB ',bb
      STOP ' AA or BB <= 0 in GENBET - Abort!'

   10 olda = aa
      oldb = bb
   20 IF (.NOT. (min(aa,bb).GT.1.0)) GO TO 100


c     Alborithm BB

c
c     Initialize
c
      IF (qsame) GO TO 30
      a = min(aa,bb)
      b = max(aa,bb)
      alpha = a + b
      beta = sqrt((alpha-2.0)/ (2.0*a*b-alpha))
      gamma = a + 1.0/beta
   30 CONTINUE
   40 u1 = ranf()
c
c     Step 1
c
      u2 = ranf()
      v = beta*log(u1/ (1.0-u1))
      IF (.NOT. (v.GT.expmax)) GO TO 50
      w = infnty
      GO TO 60

   50 w = a*exp(v)
   60 z = u1**2*u2
      r = gamma*v - 1.3862944
      s = a + r - w
c
c     Step 2
c
      IF ((s+2.609438).GE. (5.0*z)) GO TO 70
c
c     Step 3
c
      t = log(z)
      IF (s.GT.t) GO TO 70
c
c     Step 4
c
      IF ((r+alpha*log(alpha/ (b+w))).LT.t) GO TO 40
c
c     Step 5
c
   70 IF (.NOT. (aa.EQ.a)) GO TO 80
      genbet = w/ (b+w)
      GO TO 90

   80 genbet = b/ (b+w)
   90 GO TO 230


c     Algorithm BC

c
c     Initialize
c
  100 IF (qsame) GO TO 110
      a = max(aa,bb)
      b = min(aa,bb)
      alpha = a + b
      beta = 1.0/b
      delta = 1.0 + a - b
      k1 = delta* (0.0138889+0.0416667*b)/ (a*beta-0.777778)
      k2 = 0.25 + (0.5+0.25/delta)*b
  110 CONTINUE
  120 u1 = ranf()
c
c     Step 1
c
      u2 = ranf()
      IF (u1.GE.0.5) GO TO 130
c
c     Step 2
c
      y = u1*u2
      z = u1*y
      IF ((0.25*u2+z-y).GE.k1) GO TO 120
      GO TO 170
c
c     Step 3
c
  130 z = u1**2*u2
      IF (.NOT. (z.LE.0.25)) GO TO 160
      v = beta*log(u1/ (1.0-u1))
      IF (.NOT. (v.GT.expmax)) GO TO 140
      w = infnty
      GO TO 150

  140 w = a*exp(v)
  150 GO TO 200

  160 IF (z.GE.k2) GO TO 120
c
c     Step 4
c
c
c     Step 5
c
  170 v = beta*log(u1/ (1.0-u1))
      IF (.NOT. (v.GT.expmax)) GO TO 180
      w = infnty
      GO TO 190

  180 w = a*exp(v)
  190 IF ((alpha* (log(alpha/ (b+w))+v)-1.3862944).LT.log(z)) GO TO 120
c
c     Step 6
c
  200 IF (.NOT. (a.EQ.aa)) GO TO 210
      genbet = w/ (b+w)
      GO TO 220

  210 genbet = b/ (b+w)
  220 CONTINUE
  230 RETURN
      END
      REAL FUNCTION genchi(df)

c*********************************************************************72
c
cc GENCHI generates a Chi-Square random deviate.
c
c  Discussion:
c
c    Generates random deviate from the distribution of a chisquare
c    with DF degrees of freedom random variable.
c
c    The algorithm exploits the relation between chisquare and gamma.
c
c  Modified:
c
c    10 December 2007
c
c  Parameters:
c
c     DF --> Degrees of freedom of the chisquare
c            (Must be positive)
c                         REAL DF
c
      REAL df
c     ..
c     .. External Functions ..
      REAL gengam
      EXTERNAL gengam
c     ..
c     .. Executable Statements ..
      IF (.NOT. (df.LE.0.0)) GO TO 10
      WRITE (*,*) 'DF <= 0 in GENCHI - ABORT'
      WRITE (*,*) 'Value of DF: ',df
      STOP 'DF <= 0 in GENCHI - ABORT'

   10 genchi = 2.0*gengam(1.0,df/2.0)

      RETURN
      END
      REAL FUNCTION genexp(av)

c*********************************************************************72
c
cc GENEXP generates an exponential random deviate.
c
c  Discussion:
c
c    Generates a single random deviate from an exponential
c    distribution with mean AV.
c
c  Modified:
c
c    10 December 2007
c
c  Reference:
c
c    Joachim Ahrens, Ulrich Dieter,
c    Computer Methods for Sampling From the
c    Exponential and Normal Distributions,
c    Communications of the ACM,
c    Volume 15, Number 10, October 1972, pages 873-882.
c
c  Parameters:
c
c     AV --> The mean of the exponential distribution from which
c            a random deviate is to be generated.
c                              REAL AV
c
c     GENEXP <-- The random deviate.
c                              REAL GENEXP
c
c

c     .. Scalar Arguments ..
      REAL av
c     ..
c     .. External Functions ..
      REAL sexpo
      EXTERNAL sexpo
c     ..
c     .. Executable Statements ..
      genexp = sexpo()*av
      RETURN
      END
      REAL FUNCTION genf(dfn,dfd)

c*********************************************************************72
c
cc GENF generates an F random deviate.
c
c  Discussion:
c
c    Generates a random deviate from the F (variance ratio)
c    distribution with DFN degrees of freedom in the numerator
c    and DFD degrees of freedom in the denominator.
c
c  Modified:
c
c    10 December 2007
c
c  Parameters:
c
c     DFN --> Numerator degrees of freedom
c             (Must be positive)
c                              REAL DFN
c      DFD --> Denominator degrees of freedom
c             (Must be positive)
c                              REAL DFD
c
c
c                              Method
c
c
c     Directly generates ratio of chisquare variates
c
c     .. Scalar Arguments ..
      REAL dfd,dfn
c     ..
c     .. Local Scalars ..
      REAL xden,xnum
c     ..
c     .. External Functions ..
      REAL genchi
      EXTERNAL genchi
c     ..
c     .. Executable Statements ..
      IF (.NOT. (dfn.LE.0.0.OR.dfd.LE.0.0)) GO TO 10
      WRITE (*,*) 'Degrees of freedom nonpositive in GENF - abort!'
      WRITE (*,*) 'DFN value: ',dfn,'DFD value: ',dfd
      STOP 'Degrees of freedom nonpositive in GENF - abort!'

   10 xnum = genchi(dfn)/dfn
c      GENF = ( GENCHI( DFN ) / DFN ) / ( GENCHI( DFD ) / DFD )
      xden = genchi(dfd)/dfd
      IF (.NOT. (xden.LE. (1.0E-38*xnum))) GO TO 20
      WRITE (*,*) ' GENF - generated numbers would cause overflow'
      WRITE (*,*) ' Numerator ',xnum,' Denominator ',xden
      WRITE (*,*) ' GENF returning 1.0E38'
      genf = 1.0E38
      GO TO 30

   20 genf = xnum/xden
   30 RETURN
      END
      REAL FUNCTION gengam(a,r)

c*********************************************************************72
c
cc GENGAM generates a Gamma random deviate.
c
c  Discussion:
c
c    Generates random deviates from the gamma distribution whose
c    density is (A**R)/Gamma(R) * X**(R-1) * Exp(-A*X)
c
c  Modified:
c
c    10 December 2007
c
c  Reference:
c
c    Joachim Ahrens, Ulrich Dieter,
c    Generating Gamma Variates by a Modified Rejection Technique,
c    Communications of the ACM,
c    Volume 25, Number 1, January 1982, pages 47-54.
c
c    Joachim Ahrens, Ulrich Dieter,
c    Computer Methods for Sampling from Gamma, Beta, Poisson and
c    Binomial Distributions,
c    Computing,
c    Volume 12, Number 3, September 1974, pages 223-246.
c
c  Parameters:
c
c     A --> Location parameter of Gamma distribution
c                              REAL A
c
c     R --> Shape parameter of Gamma distribution
c                              REAL R
c
c
c                              Method
c
c
c     Renames SGAMMA from TOMS as slightly modified by BWB to use RANF
c     instead of SUNIF.
c

c
c     .. Scalar Arguments ..
      REAL a,r
c     ..
c     .. External Functions ..
      REAL sgamma
      EXTERNAL sgamma
c     ..
c     .. Executable Statements ..

      gengam = sgamma(r)
      gengam = gengam/a

      RETURN
      END
      SUBROUTINE genmn(parm,x,work)

c*********************************************************************72
c
cc GENMN generates a multivariate normal deviate.
c
c  Discussion:
c
c  Modified:
c
c    10 December 2007
c
c  Parameters:
c
c     PARM --> Parameters needed to generate multivariate normal
c               deviates (MEANV and Cholesky decomposition of
c               COVM). Set by a previous call to SETGMN.
c               1 : 1                - size of deviate, P
c               2 : P + 1            - mean vector
c               P+2 : P*(P+3)/2 + 1  - upper half of cholesky
c                                       decomposition of cov matrix
c                                             REAL PARM(*)
c
c     X    <-- Vector deviate generated.
c                                             REAL X(P)
c
c     WORK <--> Scratch array
c                                             REAL WORK(P)
c
c
c                              Method
c
c
c     1) Generate P independent standard normal deviates - Ei ~ N(0,1)
c
c     2) Using Cholesky decomposition find A s.t. trans(A)*A = COVM
c
c     3) trans(A)E + MEANV ~ N(MEANV,COVM)
c
c     .. Array Arguments ..
      REAL parm(*),work(*),x(*)
c     ..
c     .. Local Scalars ..
      REAL ae
      INTEGER i,icount,j,p
c     ..
c     .. External Functions ..
      REAL snorm
      EXTERNAL snorm
c     ..
c     .. Intrinsic Functions ..
      INTRINSIC int
c     ..
c     .. Executable Statements ..
      p = int(parm(1))
c
c     Generate P independent normal deviates - WORK ~ N(0,1)
c
      DO i = 1,p
          work(i) = snorm()
      end do

      DO 30,i = 1,p
c
c     PARM (P+2 : P*(P+3)/2 + 1) contains A, the Cholesky
c      decomposition of the desired covariance matrix.
c          trans(A)(1,1) = PARM(P+2)
c          trans(A)(2,1) = PARM(P+3)
c          trans(A)(2,2) = PARM(P+2+P)
c          trans(A)(3,1) = PARM(P+4)
c          trans(A)(3,2) = PARM(P+3+P)
c          trans(A)(3,3) = PARM(P+2-1+2P)  ...
c
c     trans(A)*WORK + MEANV ~ N(MEANV,COVM)
c
          icount = 0
          ae = 0.0
          DO j = 1,i
              icount = icount + j - 1
              ae = ae + parm(i+ (j-1)*p-icount+p+1)*work(j)
          end do
          x(i) = ae + parm(i+1)
   30 CONTINUE

      RETURN
      END
      SUBROUTINE genmul(n,p,ncat,ix)

c*********************************************************************72
c
cc GENMUL generates a multinomial random deviate.
c
c  Discussion:
c
c  Modified:
c
c    10 December 2007
c
c  Reference:
c
c     Luc Devroye,
c     Non-Uniform Random Variate Generation,
c     Springer, 1986,
c     ISBN: 0387963057,
c     LC: QA274.D48.
c
c  Parameters:
c
c     N --> Number of events that will be classified into one of
c           the categories 1..NCAT
c                         INTEGER N
c
c     P --> Vector of probabilities.  P(i) is the probability that
c           an event will be classified into category i.  Thus, P(i)
c           must be [0,1]. Only the first NCAT-1 P(i) must be defined
c           since P(NCAT) is 1.0 minus the sum of the first
c           NCAT-1 P(i).
c                         REAL P(NCAT-1)
c
c     NCAT --> Number of categories.  Length of P and IX.
c                         INTEGER NCAT
c
c     IX <-- Observation from multinomial distribution.  All IX(i)
c            will be nonnegative and their sum will be N.
c                         INTEGER IX(NCAT)
c
c
c     .. Scalar Arguments ..
      INTEGER n,ncat
c     ..
c     .. Array Arguments ..
      REAL p(*)
      INTEGER ix(*)
c     ..
c     .. Local Scalars ..
      REAL prob,ptot,sum
      INTEGER i,icat,ntot
c     ..
c     .. External Functions ..
      INTEGER ignbin
      EXTERNAL ignbin
c     ..
c     .. Intrinsic Functions ..
      INTRINSIC abs
c     ..
c     .. Executable Statements ..

c     Check Arguments
      IF (n.LT.0) STOP 'N < 0 in GENMUL'
      IF (ncat.LE.1) STOP 'NCAT <= 1 in GENMUL'
      ptot = 0.0
      DO 10,i = 1,ncat - 1
          IF (p(i).LT.0.0) STOP 'Some P(i) < 0 in GENMUL'
          IF (p(i).GT.1.0) STOP 'Some P(i) > 1 in GENMUL'
          ptot = ptot + p(i)
   10 CONTINUE
      IF (ptot.GT.0.99999) STOP 'Sum of P(i) > 1 in GENMUL'

c     Initialize variables
      ntot = n
      sum = 1.0
      DO i = 1,ncat
          ix(i) = 0
      end do
c
c     Generate the observation
c
      DO 30,icat = 1,ncat - 1
          prob = p(icat)/sum
          ix(icat) = ignbin(ntot,prob)
          ntot = ntot - ix(icat)
          IF (ntot.LE.0) RETURN
          sum = sum - p(icat)
   30 CONTINUE
      ix(ncat) = ntot

c     Finished
      RETURN
      END
      REAL FUNCTION gennch(df,xnonc)

c*********************************************************************72
c
cc GENNCH generates a noncentral Chi-Square random deviate.
c
c  Discussion:
c
c    Generates random deviate  from the  distribution  of a  noncentral
c    chisquare with DF degrees  of freedom and noncentrality  parameter
c    XNONC.
c
c    Uses fact that  noncentral chisquare  is  the  sum of a  chisquare
c    deviate with DF-1  degrees of freedom plus the  square of a normal
c    deviate with mean XNONC and standard deviation 1.
c
c  Modified:
c
c    10 December 2007
c
c  Parameters:
c
c     DF --> Degrees of freedom of the chisquare
c            (Must be > 1.0)
c                         REAL DF
c
c     XNONC --> Noncentrality parameter of the chisquare
c               (Must be >= 0.0)
c                         REAL XNONC
c
      REAL df,xnonc
c     ..
c     .. External Functions ..
      REAL genchi,gennor
      EXTERNAL genchi,gennor
c     ..
c     .. Intrinsic Functions ..
      INTRINSIC sqrt
c     ..
c     .. Executable Statements ..
      IF (.NOT. (df.LE.1.0.OR.xnonc.LT.0.0)) GO TO 10
      WRITE (*,*) 'DF <= 1 or XNONC < 0 in GENNCH - ABORT'
      WRITE (*,*) 'Value of DF: ',df,' Value of XNONC',xnonc
      STOP 'DF <= 1 or XNONC < 0 in GENNCH - ABORT'

   10 gennch = genchi(df-1.0) + gennor(sqrt(xnonc),1.0)**2
      RETURN
      END
      REAL FUNCTION gennf(dfn,dfd,xnonc)

c*********************************************************************72
c
cc GENNF generates a noncentral F random deviate.
c
c  Discussion:
c
c    Generates a random deviate from the  noncentral F (variance ratio)
c    distribution with DFN degrees of freedom in the numerator, and DFD
c    degrees of freedom in the denominator, and noncentrality parameter
c    XNONC.
c
c  Modified:
c
c    10 December 2007
c
c  Parameters:
c
c     DFN --> Numerator degrees of freedom
c             (Must be >= 1.0)
c                              REAL DFN
c      DFD --> Denominator degrees of freedom
c             (Must be positive)
c                              REAL DFD
c
c     XNONC --> Noncentrality parameter
c               (Must be nonnegative)
c                              REAL XNONC
c
c
c                              Method
c
c
c     Directly generates ratio of noncentral numerator chisquare variate
c     to central denominator chisquare variate.
c
c     .. Scalar Arguments ..
      REAL dfd,dfn,xnonc
c     ..
c     .. Local Scalars ..
      REAL xden,xnum
      LOGICAL qcond
c     ..
c     .. External Functions ..
      REAL genchi,gennch
      EXTERNAL genchi,gennch
c     ..
c     .. Executable Statements ..
      qcond = dfn .LE. 1.0 .OR. dfd .LE. 0.0 .OR. xnonc .LT. 0.0
      IF (.NOT. (qcond)) GO TO 10
      WRITE (*,*) 'In GENNF - Either (1) Numerator DF <= 1.0 or'
      WRITE (*,*) '(2) Denominator DF < 0.0 or '
      WRITE (*,*) '(3) Noncentrality parameter < 0.0'
      WRITE (*,*) 'DFN value: ',dfn,'DFD value: ',dfd,'XNONC value: ',
     +  xnonc
      STOP 'Degrees of freedom or noncent param our of range in GENNF'

   10 xnum = gennch(dfn,xnonc)/dfn
c      GENNF = ( GENNCH( DFN, XNONC ) / DFN ) / ( GENCHI( DFD ) / DFD )
      xden = genchi(dfd)/dfd
      IF (.NOT. (xden.LE. (1.0E-38*xnum))) GO TO 20
      WRITE (*,*) ' GENNF - generated numbers would cause overflow'
      WRITE (*,*) ' Numerator ',xnum,' Denominator ',xden
      WRITE (*,*) ' GENNF returning 1.0E38'
      gennf = 1.0E38
      GO TO 30

   20 gennf = xnum/xden
   30 RETURN
      END
      REAL FUNCTION gennor(av,sd)

c*********************************************************************72
c
cc GENNOR generates a normal random deviate.
c
c  Discussion:
c
c     Generates a single random deviate from a normal distribution
c     with mean, AV, and standard deviation, SD.
c
c  Modified:
c
c    10 December 2007
c
c  Reference:
c
c    Joachim Ahrens, Ulrich Dieter,
c    Extensions of Forsythe's Method for Random
c    Sampling from the Normal Distribution,
c    Mathematics of Computation,
c    Volume 27, Number 124, October 1973, page 927-937.
c
c  Parameters:
c
c     AV --> Mean of the normal distribution.
c                              REAL AV
c
c     SD --> Standard deviation of the normal distribution.
c                              REAL SD
c
c     GENNOR <-- Generated normal deviate.
c                              REAL GENNOR
c
c
      REAL av,sd
c     ..
c     .. External Functions ..
      REAL snorm
      EXTERNAL snorm
c     ..
c     .. Executable Statements ..
      gennor = sd*snorm() + av

      RETURN
      END
      SUBROUTINE genprm(iarray,larray)

c*********************************************************************72
c
cc GENPRM generates and applies a random permutation to an array.
c
c  Discussion:
c
c  Modified:
c
c    10 December 2007
c
c  Parameters:
c
c     IARRAY <--> On output IARRAY is a random permutation of its
c                 value on input
c                         INTEGER IARRAY( LARRAY )
c
c     LARRAY <--> Length of IARRAY
c                         INTEGER LARRAY
c
c     .. Scalar Arguments ..
      INTEGER larray
c     ..
c     .. Array Arguments ..
      INTEGER iarray(larray)
c     ..
c     .. Local Scalars ..
      INTEGER i,itmp,iwhich
c     ..
c     .. External Functions ..
      INTEGER ignuin
      EXTERNAL ignuin
c     ..
c     .. Executable Statements ..
      DO i = 1,larray
          iwhich = ignuin(i,larray)
          itmp = iarray(iwhich)
          iarray(iwhich) = iarray(i)
          iarray(i) = itmp
      end do

      RETURN
      END
      REAL FUNCTION genunf(low,high)

c*********************************************************************72
c
cc GENUNF generates a uniform random deviate.
c
c  Discussion:
c
c     Generates a real uniformly distributed between LOW and HIGH.
c
c  Modified:
c
c    10 December 2007
c
c  Parameters:
c
c     LOW --> Low bound (exclusive) on real value to be generated
c                         REAL LOW
c
c     HIGH --> High bound (exclusive) on real value to be generated
c                         REAL HIGH
c
c     .. Scalar Arguments ..
      REAL high,low
c     ..
c     .. External Functions ..
      REAL ranf
      EXTERNAL ranf
c     ..
c     .. Executable Statements ..
      IF (.NOT. (low.GT.high)) GO TO 10
      WRITE (*,*) 'LOW > HIGH in GENUNF: LOW ',low,' HIGH: ',high
      WRITE (*,*) 'Abort'
      STOP 'LOW > High in GENUNF - Abort'

   10 genunf = low + (high-low)*ranf()

      RETURN
      END
      SUBROUTINE getcgn ( g )

c*********************************************************************72
c
cc GETCGN gets the index of the current random number generator.
c
c  Discussion:
c
c  Modified:
c
c    10 December 2007
c
c  Parameters:
c
c    Output, integer G, the index of the current random number generator,
c    between 1 and 32.
c
      INTEGER curntg
      integer g
      integer numg

      SAVE curntg
      PARAMETER (numg=32)
      DATA curntg/1/

      g = curntg

      RETURN

      ENTRY setcgn(g)

c*********************************************************************72
c
cc SETCGN sets the current random number generator.
c
c  Discussion:
c
c     Sets  the  current  generator to G.    All references to a generat
c     are to the current generator.
c
c  Modified:
c
c    10 December 2007
c
c  Parameters:
c
c     G --> Number of the current random number generator (1..32)
c                    INTEGER G
c
c
c     Abort if generator number out of range
c
      IF (.NOT. (g.LT.0.OR.g.GT.numg)) GO TO 10
      WRITE (*,*) ' Generator number out of range in SETCGN:',
     +  ' Legal range is 1 to ',numg,' -- ABORT!'
      STOP ' Generator number out of range in SETCGN'

   10 curntg = g
      RETURN
      END
      SUBROUTINE getsd(iseed1,iseed2)

c*********************************************************************72
c
cc GETSD returns the value of the random number generator seeds.
c
c  Discussion:
c
c    This  is   a  transcription from  Pascal   to  Fortran  of routine
c    Get_State from the paper
c
c  Modified:
c
c    10 December 2007
c
c  Reference:
c
c    Pierre LEcuyer, Serge Cote,
c    Implementing a Random Number Package with Splitting Facilities,
c    ACM Transactions on Mathematical Software,
c    Volume 17, 1991, pages 98-111.
c
c  Parameters:
c
c     ISEED1 <- First integer seed of generator G
c                                   INTEGER ISEED1
c
c     ISEED2 <- Second integer seed of generator G
c                                   INTEGER ISEED1
c
c     .. Parameters ..
      INTEGER numg
      PARAMETER (numg=32)
c     ..
c     .. Scalar Arguments ..
      INTEGER iseed1,iseed2
c     ..
c     .. Scalars in Common ..
      INTEGER a1,a1vw,a1w,a2,a2vw,a2w,m1,m2
c     ..
c     .. Arrays in Common ..
      INTEGER cg1(numg),cg2(numg),ig1(numg),ig2(numg),lg1(numg),
     +        lg2(numg)
      LOGICAL qanti(numg)
c     ..
c     .. Local Scalars ..
      INTEGER g
c     ..
c     .. External Functions ..
      LOGICAL qrgnin
      EXTERNAL qrgnin
c     ..
c     .. External Subroutines ..
      EXTERNAL getcgn
c     ..
c     .. Common blocks ..
      COMMON /globe/m1,m2,a1,a2,a1w,a2w,a1vw,a2vw,ig1,ig2,lg1,lg2,cg1,
     +       cg2,qanti
c     ..
c     .. Save statement ..
      SAVE /globe/
c     ..
c     .. Executable Statements ..
c     Abort unless random number generator initialized
      IF (qrgnin()) GO TO 10
      WRITE (*,*) ' GETSD called before random number generator ',
     +  ' initialized -- abort!'
      STOP ' GETSD called before random number generator initialized'

   10 CALL getcgn(g)
      iseed1 = cg1(g)
      iseed2 = cg2(g)
      RETURN
      END
      INTEGER FUNCTION ignbin(n,pp)

c*********************************************************************72
c
cc IGNBIN generates a binomial random deviate.
c
c  Discussion:
c
c     Generates a single random deviate from a binomial
c     distribution whose number of trials is N and whose
c     probability of an event in each trial is P.
c
c  Modified:
c
c    10 December 2007
c
c  Reference:
c
c    Voratas Kachitvichyanukul, Bruce Schmeiser,
c    Binomial Random Variate Generation,
c    Communications of the ACM,
c    Volume 31, Number 2, February 1988, page 216-222.
c
c  Parameters:
c
c     N  --> The number of trials in the binomial distribution
c            from which a random deviate is to be generated.
c                              INTEGER N
c
c     P  --> The probability of an event in each trial of the
c            binomial distribution from which a random deviate
c            is to be generated.
c                              REAL P
c
c     IGNBIN <-- A random deviate yielding the number of events
c                from N independent trials, each of which has
c                a probability of event P.
c                              INTEGER IGNBIN
c
c
c                              Note
c
c
c     Uses RANF so the value of the seeds, ISEED1 and ISEED2 must be set
c     by a call similar to the following
c          DUM = RANSET( ISEED1, ISEED2 )
c
c
c                              Method
c
c     SUBROUTINE BTPEC(N,PP,ISEED,JX)
c
c     BINOMIAL RANDOM VARIATE GENERATOR
c     MEAN .LT. 30 -- INVERSE CDF
c       MEAN .GE. 30 -- ALGORITHM BTPE:  ACCEPTANCE-REJECTION VIA
c       FOUR REGION COMPOSITION.  THE FOUR REGIONS ARE A TRIANGLE
c       (SYMMETRIC IN THE CENTER), A PAIR OF PARALLELOGRAMS (ABOVE
c       THE TRIANGLE), AND EXPONENTIAL LEFT AND RIGHT TAILS.
c
c     BTPE REFERS TO BINOMIAL-TRIANGLE-PARALLELOGRAM-EXPONENTIAL.
c     BTPEC REFERS TO BTPE AND "COMBINED."  THUS BTPE IS THE
c       RESEARCH AND BTPEC IS THE IMPLEMENTATION OF A COMPLETE
c       USABLE ALGORITHM.
c     REFERENCE:  VORATAS KACHITVICHYANUKUL AND BRUCE SCHMEISER,
c       "BINOMIAL RANDOM VARIATE GENERATION,"
c       COMMUNICATIONS OF THE ACM, FORTHCOMING
c     WRITTEN:  SEPTEMBER 1980.
c       LAST REVISED:  MAY 1985, JULY 1987
c     REQUIRED SUBPROGRAM:  RAND() -- A UNIFORM (0,1) RANDOM NUMBER
c                           GENERATOR
c     ARGUMENTS
c
c       N : NUMBER OF BERNOULLI TRIALS            (INPUT)
c       PP : PROBABILITY OF SUCCESS IN EACH TRIAL (INPUT)
c       ISEED:  RANDOM NUMBER SEED                (INPUT AND OUTPUT)
c       JX:  RANDOMLY GENERATED OBSERVATION       (OUTPUT)
c
c     VARIABLES
c       PSAVE: VALUE OF PP FROM THE LAST CALL TO BTPEC
c       NSAVE: VALUE OF N FROM THE LAST CALL TO BTPEC
c       XNP:  VALUE OF THE MEAN FROM THE LAST CALL TO BTPEC
c
c       P: PROBABILITY USED IN THE GENERATION PHASE OF BTPEC
c       FFM: TEMPORARY VARIABLE EQUAL TO XNP + P
c       M:  INTEGER VALUE OF THE CURRENT MODE
c       FM:  FLOATING POINT VALUE OF THE CURRENT MODE
c       XNPQ: TEMPORARY VARIABLE USED IN SETUP AND SQUEEZING STEPS
c       P1:  AREA OF THE TRIANGLE
c       C:  HEIGHT OF THE PARALLELOGRAMS
c       XM:  CENTER OF THE TRIANGLE
c       XL:  LEFT END OF THE TRIANGLE
c       XR:  RIGHT END OF THE TRIANGLE
c       AL:  TEMPORARY VARIABLE
c       XLL:  RATE FOR THE LEFT EXPONENTIAL TAIL
c       XLR:  RATE FOR THE RIGHT EXPONENTIAL TAIL
c       P2:  AREA OF THE PARALLELOGRAMS
c       P3:  AREA OF THE LEFT EXPONENTIAL TAIL
c       P4:  AREA OF THE RIGHT EXPONENTIAL TAIL
c       U:  A U(0,P4) RANDOM VARIATE USED FIRST TO SELECT ONE OF THE
c           FOUR REGIONS AND THEN CONDITIONALLY TO GENERATE A VALUE
c           FROM THE REGION
c       V:  A U(0,1) RANDOM NUMBER USED TO GENERATE THE RANDOM VALUE
c           (REGION 1) OR TRANSFORMED INTO THE VARIATE TO ACCEPT OR
c           REJECT THE CANDIDATE VALUE
c       IX:  INTEGER CANDIDATE VALUE
c       X:  PRELIMINARY CONTINUOUS CANDIDATE VALUE IN REGION 2 LOGIC
c           AND A FLOATING POINT IX IN THE ACCEPT/REJECT LOGIC
c       K:  ABSOLUTE VALUE OF (IX-M)
c       F:  THE HEIGHT OF THE SCALED DENSITY FUNCTION USED IN THE
c           ACCEPT/REJECT DECISION WHEN BOTH M AND IX ARE SMALL
c           ALSO USED IN THE INVERSE TRANSFORMATION
c       R: THE RATIO P/Q
c       G: CONSTANT USED IN CALCULATION OF PROBABILITY
c       MP:  MODE PLUS ONE, THE LOWER INDEX FOR EXPLICIT CALCULATION
c            OF F WHEN IX IS GREATER THAN M
c       IX1:  CANDIDATE VALUE PLUS ONE, THE LOWER INDEX FOR EXPLICIT
c             CALCULATION OF F WHEN IX IS LESS THAN M
c       I:  INDEX FOR EXPLICIT CALCULATION OF F FOR BTPE
c       AMAXP: MAXIMUM ERROR OF THE LOGARITHM OF NORMAL BOUND
c       YNORM: LOGARITHM OF NORMAL BOUND
c       ALV:  NATURAL LOGARITHM OF THE ACCEPT/REJECT VARIATE V
c
c       X1,F1,Z,W,Z2,X2,F2, AND W2 ARE TEMPORARY VARIABLES TO BE
c       USED IN THE FINAL ACCEPT/REJECT TEST
c
c       QN: PROBABILITY OF NO SUCCESS IN N TRIALS
c
c     REMARK
c       IX AND JX COULD LOGICALLY BE THE SAME VARIABLE, WHICH WOULD
c       SAVE A MEMORY POSITION AND A LINE OF CODE.  HOWEVER, SOME
c       COMPILERS (E.G.,CDC MNF) OPTIMIZE BETTER WHEN THE ARGUMENTS
c       ARE NOT INVOLVED.
c
c     ISEED NEEDS TO BE DOUBLE PRECISION IF THE IMSL ROUTINE
c     GGUBFS IS USED TO GENERATE UNIFORM RANDOM NUMBER, OTHERWISE
c     TYPE OF ISEED SHOULD BE DICTATED BY THE UNIFORM GENERATOR
c
c  DETERMINE APPROPRIATE ALGORITHM AND WHETHER SETUP IS NECESSARY
c
c     ..
c     .. Scalar Arguments ..
      REAL pp
      INTEGER n
c     ..
c     .. Local Scalars ..
      REAL al,alv,amaxp,c,f,f1,f2,ffm,fm,g,p,p1,p2,p3,p4,psave,q,qn,r,u,
     +     v,w,w2,x,x1,x2,xl,xll,xlr,xm,xnp,xnpq,xr,ynorm,z,z2
      INTEGER i,ix,ix1,k,m,mp,nsave
c     ..
c     .. External Functions ..
      REAL ranf
      EXTERNAL ranf
c     ..
c     .. Intrinsic Functions ..
      INTRINSIC abs,alog,amin1,iabs,int,sqrt
c     ..
c     .. Data statements ..
      DATA psave,nsave/-1.,-1/
c     ..
c     .. Executable Statements ..
      IF (pp.NE.psave) GO TO 10
      IF (n.NE.nsave) GO TO 20
      IF (xnp-30.) 150,30,30
c
c  SETUP, PERFORM ONLY WHEN PARAMETERS CHANGE
c
   10 psave = pp
      p = amin1(psave,1.-psave)
      q = 1. - p
   20 xnp = n*p
      nsave = n
      IF (xnp.LT.30.) GO TO 140
      ffm = xnp + p
      m = ffm
      fm = m
      xnpq = xnp*q
      p1 = int(2.195*sqrt(xnpq)-4.6*q) + 0.5
      xm = fm + 0.5
      xl = xm - p1
      xr = xm + p1
      c = 0.134 + 20.5/ (15.3+fm)
      al = (ffm-xl)/ (ffm-xl*p)
      xll = al* (1.+.5*al)
      al = (xr-ffm)/ (xr*q)
      xlr = al* (1.+.5*al)
      p2 = p1* (1.+c+c)
      p3 = p2 + c/xll
      p4 = p3 + c/xlr
c      WRITE(6,100) N,P,P1,P2,P3,P4,XL,XR,XM,FM
c  100 FORMAT(I15,4F18.7/5F18.7)
c
c  GENERATE VARIATE
c
   30 u = ranf()*p4
      v = ranf()
c
c     TRIANGULAR REGION
c
      IF (u.GT.p1) GO TO 40
      ix = xm - p1*v + u
      GO TO 170
c
c     PARALLELOGRAM REGION
c
   40 IF (u.GT.p2) GO TO 50
      x = xl + (u-p1)/c
      v = v*c + 1. - abs(xm-x)/p1
      IF (v.GT.1. .OR. v.LE.0.) GO TO 30
      ix = x
      GO TO 70
c
c     LEFT TAIL
c
   50 IF (u.GT.p3) GO TO 60
      ix = xl + alog(v)/xll
      IF (ix.LT.0) GO TO 30
      v = v* (u-p2)*xll
      GO TO 70
c
c     RIGHT TAIL
c
   60 ix = xr - alog(v)/xlr
      IF (ix.GT.n) GO TO 30
      v = v* (u-p3)*xlr
c
c  DETERMINE APPROPRIATE WAY TO PERFORM ACCEPT/REJECT TEST
c
   70 k = iabs(ix-m)
      IF (k.GT.20 .AND. k.LT.xnpq/2-1) GO TO 130
c
c     EXPLICIT EVALUATION
c
      f = 1.0
      r = p/q
      g = (n+1)*r
      IF (m-ix) 80,120,100
   80 mp = m + 1
      DO 90 i = mp,ix
          f = f* (g/i-r)
   90 CONTINUE
      GO TO 120

  100 ix1 = ix + 1
      DO 110 i = ix1,m
          f = f/ (g/i-r)
  110 CONTINUE
  120 IF (v-f) 170,170,30
c
c     SQUEEZING USING UPPER AND LOWER BOUNDS ON ALOG(F(X))
c
  130 amaxp = (k/xnpq)* ((k* (k/3.+.625)+.1666666666666)/xnpq+.5)
      ynorm = -k*k/ (2.*xnpq)
      alv = alog(v)
      IF (alv.LT.ynorm-amaxp) GO TO 170
      IF (alv.GT.ynorm+amaxp) GO TO 30
c
c     STIRLING'S FORMULA TO MACHINE ACCURACY FOR
c     THE FINAL ACCEPTANCE/REJECTION TEST
c
      x1 = ix + 1
      f1 = fm + 1.
      z = n + 1 - fm
      w = n - ix + 1.
      z2 = z*z
      x2 = x1*x1
      f2 = f1*f1
      w2 = w*w
      IF (alv- (xm*alog(f1/x1)+ (n-m+.5)*alog(z/w)+ (ix-
     +    m)*alog(w*p/ (x1*q))+ (13860.- (462.- (132.- (99.-
     +    140./f2)/f2)/f2)/f2)/f1/166320.+ (13860.- (462.- (132.- (99.-
     +    140./z2)/z2)/z2)/z2)/z/166320.+ (13860.- (462.- (132.- (99.-
     +    140./x2)/x2)/x2)/x2)/x1/166320.+ (13860.- (462.- (132.- (99.-
     +    140./w2)/w2)/w2)/w2)/w/166320.)) 170,170,30
c
c     INVERSE CDF LOGIC FOR MEAN LESS THAN 30
c
  140 qn = q**n
      r = p/q
      g = r* (n+1)
  150 ix = 0
      f = qn
      u = ranf()
  160 IF (u.LT.f) GO TO 170
      IF (ix.GT.110) GO TO 150
      u = u - f
      ix = ix + 1
      f = f* (g/ix-r)
      GO TO 160

  170 IF (psave.GT.0.5) ix = n - ix
      ignbin = ix
      RETURN
      END
      INTEGER FUNCTION ignlgi()

c*********************************************************************72
c
cc IGNLGI generates a random positive integer.
c
c  Discussion:
c
c    Returns a random integer following a uniform distribution over
c    (1, 2147483562) using the current generator.
c
c    This is a transcription from Pascal to Fortran of routine
c    Random from the paper
c
c  Modified:
c
c    10 December 2007
c
c  Reference:
c
c    Pierre LEcuyer, Serge Cote,
c    Implementing a Random Number Package with Splitting Facilities,
c    ACM Transactions on Mathematical Software,
c    Volume 17, 1991, pages 98-111.
c
c  Parameters:
c
      INTEGER numg
      PARAMETER (numg=32)
c     ..
c     .. Scalars in Common ..
      INTEGER a1,a1vw,a1w,a2,a2vw,a2w,m1,m2
c     ..
c     .. Arrays in Common ..
      INTEGER cg1(numg),cg2(numg),ig1(numg),ig2(numg),lg1(numg),
     +        lg2(numg)
      LOGICAL qanti(numg)
c     ..
c     .. Local Scalars ..
      INTEGER curntg,k,s1,s2,z
      LOGICAL qqssd
c     ..
c     .. External Functions ..
      LOGICAL qrgnin
      EXTERNAL qrgnin
c     ..
c     .. External Subroutines ..
      EXTERNAL getcgn,inrgcm,rgnqsd,setall
c     ..
c     .. Common blocks ..
      COMMON /globe/m1,m2,a1,a2,a1w,a2w,a1vw,a2vw,ig1,ig2,lg1,lg2,cg1,
     +       cg2,qanti
c     ..
c     .. Save statement ..
      SAVE /globe/
c     ..
c     .. Executable Statements ..
c
c     IF THE RANDOM NUMBER PACKAGE HAS NOT BEEN INITIALIZED YET, DO SO.
c     IT CAN BE INITIALIZED IN ONE OF TWO WAYS : 1) THE FIRST CALL TO
c     THIS ROUTINE  2) A CALL TO SETALL.
c
      IF (.NOT. (qrgnin())) CALL inrgcm()
      CALL rgnqsd(qqssd)
      IF (.NOT. (qqssd)) CALL setall(1234567890,123456789)
c
c     Get Current Generator
c
      CALL getcgn(curntg)
      s1 = cg1(curntg)
      s2 = cg2(curntg)
      k = s1/53668
      s1 = a1* (s1-k*53668) - k*12211
      IF (s1.LT.0) s1 = s1 + m1
      k = s2/52774
      s2 = a2* (s2-k*52774) - k*3791
      IF (s2.LT.0) s2 = s2 + m2
      cg1(curntg) = s1
      cg2(curntg) = s2
      z = s1 - s2
      IF (z.LT.1) z = z + m1 - 1
      IF (qanti(curntg)) z = m1 - z
      ignlgi = z
      RETURN
      END
      INTEGER FUNCTION ignnbn(n,p)

c*********************************************************************72
c
cc IGNNBN generates a negative binomial random deviate.
c
c  Discussion:
c
c    Generates a single random deviate from a negative binomial
c    distribution.
c
c  Modified:
c
c    10 December 2007
c
c  Reference:
c
c    Devroye, Luc
c    Non-Uniform Random Variate Generation,
c    Springer, 1986,
c    ISBN: 0387963057,
c    LC: QA274.D48.
c
c  Parameters:
c
c     N  --> Required number of events.
c                              INTEGER N
c
c     P  --> The probability of an event during a Bernoulli trial.
c                              REAL P
c
      REAL p
      INTEGER n
c     ..
c     .. Local Scalars ..
      REAL y,a,r
c     ..
c     .. External Functions ..
      REAL gengam
      INTEGER ignpoi
      EXTERNAL gengam,ignpoi
c     ..
c     .. Intrinsic Functions ..
      INTRINSIC real
c     ..
c     .. Executable Statements ..
c     Check Arguments
      IF (n.LT.0) STOP 'N < 0 in IGNNBN'
      IF (p.LE.0.0) STOP 'P <= 0 in IGNNBN'
      IF (p.GE.1.0) STOP 'P >= 1 in IGNNBN'
c
c     Generate Y, a random gamma (n,(1-p)/p) variable
c
      r = real(n)
      a = p/ (1.0-p)
      y = gengam(a,r)
c
c     Generate a random Poisson(y) variable.
c
      ignnbn = ignpoi(y)
      RETURN
      END
      INTEGER FUNCTION ignpoi(mu)

c*********************************************************************72
c
cc IGNPOI generates a Poisson random deviate.
c
c  Discussion:
c
c    Generates a single random deviate from a Poisson
c    distribution with mean AV.
c
c  Modified:
c
c    10 December 2007
c
c  Reference:
c
c    Joachim Ahrens, Ulrich Dieter,
c    Computer Generation of Poisson Deviates
c    From Modified Normal Distributions,
c    ACM Transactions on Mathematical Software,
c    Volume 8, Number 2, June 1982, pages 163-179.
c
c  Parameters:
c
c     AV --> The mean of the Poisson distribution from which
c            a random deviate is to be generated.
c                              REAL AV
c
c     GENEXP <-- The random deviate.
c                              REAL GENEXP
c
c
c                              Method
c
c
c     Renames KPOIS from TOMS as slightly modified by BWB to use RANF
c     instead of SUNIF.
c                                                                  
c      INTEGER FUNCTION IGNPOI(IR,MU)
c
c     INPUT:  IR=CURRENT STATE OF BASIC RANDOM NUMBER GENERATOR
c             MU=MEAN MU OF THE POISSON DISTRIBUTION
c     OUTPUT: IGNPOI=SAMPLE FROM THE POISSON-(MU)-DISTRIBUTION
c
c
c
c     MUPREV=PREVIOUS MU, MUOLD=MU AT LAST EXECUTION OF STEP P OR B.
c     TABLES: COEFFICIENTS A0-A7 FOR STEP F. FACTORIALS FACT
c     COEFFICIENTS A(K) - FOR PX = FK*V*V*SUM(A(K)*V**K)-DEL
c
c
c
c     SEPARATION OF CASES A AND B
c
c     .. Scalar Arguments ..
      REAL mu
c     ..
c     .. Local Scalars ..
      REAL a0,a1,a2,a3,a4,a5,a6,a7,b1,b2,c,c0,c1,c2,c3,d,del,difmuk,e,
     +     fk,fx,fy,g,muold,muprev,omega,p,p0,px,py,q,s,t,u,v,x,xx
      INTEGER j,k,kflag,l,m
c     ..
c     .. Local Arrays ..
      REAL fact(10),pp(35)
c     ..
c     .. External Functions ..
      REAL ranf,sexpo,snorm
      EXTERNAL ranf,sexpo,snorm
c     ..
c     .. Intrinsic Functions ..
      INTRINSIC abs,alog,exp,float,ifix,max0,min0,sign,sqrt
c     ..
c     .. Data statements ..
      DATA muprev,muold/0.,0./
      DATA a0,a1,a2,a3,a4,a5,a6,a7/-.5,.3333333,-.2500068,.2000118,
     +     -.1661269,.1421878,-.1384794,.1250060/
      DATA fact/1.,1.,2.,6.,24.,120.,720.,5040.,40320.,362880./
c     ..
c     .. Executable Statements ..
      IF (mu.EQ.muprev) GO TO 10
      IF (mu.LT.10.0) GO TO 120
c
c     C A S E  A. (RECALCULATION OF S,D,L IF MU HAS CHANGED)
c
      muprev = mu
      s = sqrt(mu)
      d = 6.0*mu*mu
c
c             THE POISSON PROBABILITIES PK EXCEED THE DISCRETE NORMAL
c             PROBABILITIES FK WHENEVER K >= M(MU). L=IFIX(MU-1.1484)
c             IS AN UPPER BOUND TO M(MU) FOR ALL MU >= 10 .
c
      l = ifix(mu-1.1484)
c
c     STEP N. NORMAL SAMPLE - SNORM(IR) FOR STANDARD NORMAL DEVIATE
c
   10 g = mu + s*snorm()
      IF (g.LT.0.0) GO TO 20
      ignpoi = ifix(g)
c
c     STEP I. IMMEDIATE ACCEPTANCE IF IGNPOI IS LARGE ENOUGH
c
      IF (ignpoi.GE.l) RETURN
c
c     STEP S. SQUEEZE ACCEPTANCE - SUNIF(IR) FOR (0,1)-SAMPLE U
c
      fk = float(ignpoi)
      difmuk = mu - fk
      u = ranf()
      IF (d*u.GE.difmuk*difmuk*difmuk) RETURN
c
c     STEP P. PREPARATIONS FOR STEPS Q AND H.
c             (RECALCULATIONS OF PARAMETERS IF NECESSARY)
c             .3989423=(2*PI)**(-.5)  .416667E-1=1./24.  .1428571=1./7.
c             THE QUANTITIES B1, B2, C3, C2, C1, C0 ARE FOR THE HERMITE
c             APPROXIMATIONS TO THE DISCRETE NORMAL PROBABILITIES FK.
c             C=.1069/MU GUARANTEES MAJORIZATION BY THE 'HAT'-FUNCTION.
c
   20 IF (mu.EQ.muold) GO TO 30
      muold = mu
      omega = .3989423/s
      b1 = .4166667E-1/mu
      b2 = .3*b1*b1
      c3 = .1428571*b1*b2
      c2 = b2 - 15.*c3
      c1 = b1 - 6.*b2 + 45.*c3
      c0 = 1. - b1 + 3.*b2 - 15.*c3
      c = .1069/mu
   30 IF (g.LT.0.0) GO TO 50
c
c             'SUBROUTINE' F IS CALLED (KFLAG=0 FOR CORRECT RETURN)
c
      kflag = 0
      GO TO 70
c
c     STEP Q. QUOTIENT ACCEPTANCE (RARE CASE)
c
   40 IF (fy-u*fy.LE.py*exp(px-fx)) RETURN
c
c     STEP E. EXPONENTIAL SAMPLE - SEXPO(IR) FOR STANDARD EXPONENTIAL
c             DEVIATE E AND SAMPLE T FROM THE LAPLACE 'HAT'
c             (IF T <= -.6744 THEN PK < FK FOR ALL MU >= 10.)
c
   50 e = sexpo()
      u = ranf()
      u = u + u - 1.0
      t = 1.8 + sign(e,u)
      IF (t.LE. (-.6744)) GO TO 50
      ignpoi = ifix(mu+s*t)
      fk = float(ignpoi)
      difmuk = mu - fk
c
c             'SUBROUTINE' F IS CALLED (KFLAG=1 FOR CORRECT RETURN)
c
      kflag = 1
      GO TO 70
c
c     STEP H. HAT ACCEPTANCE (E IS REPEATED ON REJECTION)
c
   60 IF (c*abs(u).GT.py*exp(px+e)-fy*exp(fx+e)) GO TO 50
      RETURN
c
c     STEP F. 'SUBROUTINE' F. CALCULATION OF PX,PY,FX,FY.
c             CASE IGNPOI .LT. 10 USES FACTORIALS FROM TABLE FACT
c
   70 IF (ignpoi.GE.10) GO TO 80
      px = -mu
      py = mu**ignpoi/fact(ignpoi+1)
      GO TO 110
c
c             CASE IGNPOI .GE. 10 USES POLYNOMIAL APPROXIMATION
c             A0-A7 FOR ACCURACY WHEN ADVISABLE
c             .8333333E-1=1./12.  .3989423=(2*PI)**(-.5)
c
   80 del = .8333333E-1/fk
      del = del - 4.8*del*del*del
      v = difmuk/fk
      IF (abs(v).LE.0.25) GO TO 90
      px = fk*alog(1.0+v) - difmuk - del
      GO TO 100

   90 px = fk*v*v* (((((((a7*v+a6)*v+a5)*v+a4)*v+a3)*v+a2)*v+a1)*v+a0) -
     +     del
  100 py = .3989423/sqrt(fk)
  110 x = (0.5-difmuk)/s
      xx = x*x
      fx = -0.5*xx
      fy = omega* (((c3*xx+c2)*xx+c1)*xx+c0)
      IF (kflag) 40,40,60
c
c     C A S E  B. (START NEW TABLE AND CALCULATE P0 IF NECESSARY)
c
  120 muprev = 0.0
      IF (mu.EQ.muold) GO TO 130
      muold = mu
      m = max0(1,ifix(mu))
      l = 0
      p = exp(-mu)
      q = p
      p0 = p
c
c     STEP U. UNIFORM SAMPLE FOR INVERSION METHOD
c
  130 u = ranf()
      ignpoi = 0
      IF (u.LE.p0) RETURN
c
c     STEP T. TABLE COMPARISON UNTIL THE END PP(L) OF THE
c             PP-TABLE OF CUMULATIVE POISSON PROBABILITIES
c             (0.458=PP(9) FOR MU=10)
c
      IF (l.EQ.0) GO TO 150
      j = 1
      IF (u.GT.0.458) j = min0(l,m)
      DO 140 k = j,l
          IF (u.LE.pp(k)) GO TO 180
  140 CONTINUE
      IF (l.EQ.35) GO TO 130
c
c     STEP C. CREATION OF NEW POISSON PROBABILITIES P
c             AND THEIR CUMULATIVES Q=PP(K)
c
  150 l = l + 1
      DO 160 k = l,35
          p = p*mu/float(k)
          q = q + p
          pp(k) = q
          IF (u.LE.q) GO TO 170
  160 CONTINUE
      l = 35
      GO TO 130

  170 l = k
  180 ignpoi = k
      RETURN
      END
      INTEGER FUNCTION ignuin(low,high)

c*********************************************************************72
c
cc IGNUIN generates a random integer in a given range.
c
c  Discussion:
c
c  Modified:
c
c    10 December 2007
c
c  Parameters:
c
c     LOW --> Low bound (inclusive) on integer value to be generated
c                         INTEGER LOW
c
c     HIGH --> High bound (inclusive) on integer value to be generated
c                         INTEGER HIGH
c
c
c                              Note
c
c
c     If (HIGH-LOW) > 2,147,483,561 prints error message on * unit and
c     stops the program.
c

c     IGNLGI generates integers between 1 and 2147483562
c     MAXNUM is 1 less than maximum generable value
c     .. Parameters ..
      INTEGER maxnum
      PARAMETER (maxnum=2147483561)
      CHARACTER*(*) err1,err2
      PARAMETER (err1='LOW > HIGH in IGNUIN',
     +          err2=' ( HIGH - LOW ) > 2,147,483,561 in IGNUIN')
c     ..
c     .. Scalar Arguments ..
      INTEGER high,low
c     ..
c     .. Local Scalars ..
      INTEGER err,ign,maxnow,range,ranp1
c     ..
c     .. External Functions ..
      INTEGER ignlgi
      EXTERNAL ignlgi
c     ..
c     .. Intrinsic Functions ..
      INTRINSIC mod
c     ..
c     .. Executable Statements ..
      IF (.NOT. (low.GT.high)) GO TO 10
      err = 1
c      ABORT-PROGRAM
      GO TO 80

   10 range = high - low
      IF (.NOT. (range.GT.maxnum)) GO TO 20
      err = 2
c      ABORT-PROGRAM
      GO TO 80

   20 IF (.NOT. (low.EQ.high)) GO TO 30
      ignuin = low
      RETURN

      GO TO 70

c     Number to be generated should be in range 0..RANGE
c     Set MAXNOW so that the number of integers in 0..MAXNOW is an
c     integral multiple of the number in 0..RANGE

   30 ranp1 = range + 1
      maxnow = (maxnum/ranp1)*ranp1
   40 ign = ignlgi() - 1
      IF (.NOT. (ign.LE.maxnow)) GO TO 50
      ignuin = low + mod(ign,ranp1)
      RETURN

   50 GO TO 40

   60 CONTINUE
   70 CONTINUE
   80 IF (.NOT. (err.EQ.1)) GO TO 90
      WRITE (*,*) err1
      GO TO 100

c     TO ABORT-PROGRAM
   90 WRITE (*,*) err2
  100 WRITE (*,*) ' LOW: ',low,' HIGH: ',high
      WRITE (*,*) ' Abort on Fatal ERROR'
      IF (.NOT. (err.EQ.1)) GO TO 110
      STOP 'LOW > HIGH in IGNUIN'

      GO TO 120

  110 STOP ' ( HIGH - LOW ) > 2,147,483,561 in IGNUIN'

  120 END
      SUBROUTINE initgn(isdtyp)

c*********************************************************************72
c
cc INITGN initializes the current random number generator.
c
c  Discussion:
c
c    Reinitializes the state of the current generator
c
c    This is a transcription from Pascal to Fortran of routine
c    Init_Generator from the paper
c
c  Modified:
c
c    10 December 2007
c
c  Reference:
c
c    Pierre LEcuyer, Serge Cote,
c    Implementing a Random Number Package with Splitting Facilities,
c    ACM Transactions on Mathematical Software,
c    Volume 17, 1991, pages 98-111.
c
c  Parameters:
c
c     ISDTYP -> The state to which the generator is to be set
c
c          ISDTYP = -1  => sets the seeds to their initial value
c          ISDTYP =  0  => sets the seeds to the first value of
c                          the current block
c          ISDTYP =  1  => sets the seeds to the first value of
c                          the next block
c
c                                   INTEGER ISDTYP
c
c     .. Parameters ..
      INTEGER numg
      PARAMETER (numg=32)
c     ..
c     .. Scalar Arguments ..
      INTEGER isdtyp
c     ..
c     .. Scalars in Common ..
      INTEGER a1,a1vw,a1w,a2,a2vw,a2w,m1,m2
c     ..
c     .. Arrays in Common ..
      INTEGER cg1(numg),cg2(numg),ig1(numg),ig2(numg),lg1(numg),
     +        lg2(numg)
      LOGICAL qanti(numg)
c     ..
c     .. Local Scalars ..
      INTEGER g
c     ..
c     .. External Functions ..
      LOGICAL qrgnin
      INTEGER mltmod
      EXTERNAL qrgnin,mltmod
c     ..
c     .. External Subroutines ..
      EXTERNAL getcgn
c     ..
c     .. Common blocks ..
      COMMON /globe/m1,m2,a1,a2,a1w,a2w,a1vw,a2vw,ig1,ig2,lg1,lg2,cg1,
     +       cg2,qanti
c     ..
c     .. Save statement ..
      SAVE /globe/
c     ..
c     .. Executable Statements ..
c     Abort unless random number generator initialized
      IF (qrgnin()) GO TO 10
      WRITE (*,*) ' INITGN called before random number generator ',
     +  ' initialized -- abort!'
      STOP ' INITGN called before random number generator initialized'

   10 CALL getcgn(g)
      IF ((-1).NE. (isdtyp)) GO TO 20
      lg1(g) = ig1(g)
      lg2(g) = ig2(g)
      GO TO 50

   20 IF ((0).NE. (isdtyp)) GO TO 30
      CONTINUE
      GO TO 50
c     do nothing
   30 IF ((1).NE. (isdtyp)) GO TO 40
      lg1(g) = mltmod(a1w,lg1(g),m1)
      lg2(g) = mltmod(a2w,lg2(g),m2)
      GO TO 50

   40 STOP 'ISDTYP NOT IN RANGE'

   50 cg1(g) = lg1(g)
      cg2(g) = lg2(g)
      RETURN
      END
      SUBROUTINE inrgcm()

c*********************************************************************72
c
cc INRGCM initializes the random number generator common memory.
c
c  Discussion:
c
c    Initializes common area  for random number  generator.  This saves
c    the  nuisance  of  a  BLOCK DATA  routine  and the  difficulty  of
c    assuring that the routine is loaded with the other routines.
c
c  Modified:
c
c    10 December 2007
c
      INTEGER numg
      PARAMETER (numg=32)
c     ..
c     .. Scalars in Common ..
      INTEGER a1,a1vw,a1w,a2,a2vw,a2w,m1,m2
c     ..
c     .. Arrays in Common ..
      INTEGER cg1(numg),cg2(numg),ig1(numg),ig2(numg),lg1(numg),
     +        lg2(numg)
      LOGICAL qanti(numg)
c     ..
c     .. Local Scalars ..
      INTEGER i
      LOGICAL qdum
c     ..
c     .. External Functions ..
      LOGICAL qrgnsn
      EXTERNAL qrgnsn
c     ..
c     .. Common blocks ..
      COMMON /globe/m1,m2,a1,a2,a1w,a2w,a1vw,a2vw,ig1,ig2,lg1,lg2,cg1,
     +       cg2,qanti
c     ..
c     .. Save statement ..
      SAVE /globe/
c     ..
c     .. Executable Statements ..
c     V=20;                            W=30;
c
c     A1W = MOD(A1**(2**W),M1)         A2W = MOD(A2**(2**W),M2)
c     A1VW = MOD(A1**(2**(V+W)),M1)    A2VW = MOD(A2**(2**(V+W)),M2)
c
c   If V or W is changed A1W, A2W, A1VW, and A2VW need to be recomputed.
c    An efficient way to precompute a**(2*j) MOD m is to start with
c    a and square it j times modulo m using the function MLTMOD.
c
      m1 = 2147483563
      m2 = 2147483399
      a1 = 40014
      a2 = 40692
      a1w = 1033780774
      a2w = 1494757890
      a1vw = 2082007225
      a2vw = 784306273
      DO 10,i = 1,numg
          qanti(i) = .FALSE.
   10 CONTINUE
c
c     Tell the world that common has been initialized
c
      qdum = qrgnsn(.TRUE.)
      RETURN
      END
      INTEGER FUNCTION lennob(string)
      IMPLICIT INTEGER (a-p,r-z),LOGICAL (q)

c*********************************************************************72
c
cc LENNOB counts the length of a string, ignoring trailing blanks.
c
c  Discussion:
c
c    Returns the length of STRING up to and including the last
c    non-blank character.
c
c  Modified:
c
c    10 December 2007
c
c  Parameters:
c
c     STRING --> String whose length not counting trailing blanks
c                is returned.
c
      CHARACTER*(*) string

      end = len(string)
      DO 20,i = end,1,-1
          IF (.NOT. (string(i:i).NE.' ')) GO TO 10
          lennob = i
          RETURN

   10     CONTINUE
   20 CONTINUE
      lennob = 0
      RETURN
      END
      INTEGER FUNCTION mltmod(a,s,m)

c*********************************************************************72
c
cc MLTMOD carries out modular multiplication.
c
c  Discussion:
c
c    Returns (A*S) MOD M
c
c    This is a transcription from Pascal to Fortran of routine
c    MULtMod_Decompos from the paper
c
c  Modified:
c
c    10 December 2007
c
c  Reference:
c
c    Pierre LEcuyer, Serge Cote,
c    Implementing a Random Number Package with Splitting Facilities,
c    ACM Transactions on Mathematical Software,
c    Volume 17, 1991, pages 98-111.
c
c  Parameters:
c
c     A, S, M  -->
c                         INTEGER A,S,M
c
c     .. Parameters ..
      INTEGER h
      PARAMETER (h=32768)
c     ..
c     .. Scalar Arguments ..
      INTEGER a,m,s
c     ..
c     .. Local Scalars ..
      INTEGER a0,a1,k,p,q,qh,rh
c     ..
c     .. Executable Statements ..
c
c     H = 2**((b-2)/2) where b = 32 because we are using a 32 bit
c      machine. On a different machine recompute H
c
      IF (.NOT. (a.LE.0.OR.a.GE.m.OR.s.LE.0.OR.s.GE.m)) GO TO 10
      WRITE (*,*) ' A, M, S out of order in MLTMOD - ABORT!'
      WRITE (*,*) ' A = ',a,' S = ',s,' M = ',m
      WRITE (*,*) ' MLTMOD requires: 0 < A < M; 0 < S < M'
      STOP ' A, M, S out of order in MLTMOD - ABORT!'

   10 IF (.NOT. (a.LT.h)) GO TO 20
      a0 = a
      p = 0
      GO TO 120

   20 a1 = a/h
      a0 = a - h*a1
      qh = m/h
      rh = m - h*qh
      IF (.NOT. (a1.GE.h)) GO TO 50
      a1 = a1 - h
      k = s/qh
      p = h* (s-k*qh) - k*rh
   30 IF (.NOT. (p.LT.0)) GO TO 40
      p = p + m
      GO TO 30

   40 GO TO 60

   50 p = 0
c
c     P = (A2*S*H)MOD M
c
   60 IF (.NOT. (a1.NE.0)) GO TO 90
      q = m/a1
      k = s/q
      p = p - k* (m-a1*q)
      IF (p.GT.0) p = p - m
      p = p + a1* (s-k*q)
   70 IF (.NOT. (p.LT.0)) GO TO 80
      p = p + m
      GO TO 70

   80 CONTINUE
   90 k = p/qh
c
c     P = ((A2*H + A1)*S)MOD M
c
      p = h* (p-k*qh) - k*rh
  100 IF (.NOT. (p.LT.0)) GO TO 110
      p = p + m
      GO TO 100

  110 CONTINUE
  120 IF (.NOT. (a0.NE.0)) GO TO 150
c
c     P = ((A2*H + A1)*H*S)MOD M
c
      q = m/a0
      k = s/q
      p = p - k* (m-a0*q)
      IF (p.GT.0) p = p - m
      p = p + a0* (s-k*q)
  130 IF (.NOT. (p.LT.0)) GO TO 140
      p = p + m
      GO TO 130

  140 CONTINUE
  150 mltmod = p
c
      RETURN
      END
      SUBROUTINE phrtsd(phrase,seed1,seed2)

c*********************************************************************72
c
cc PHRTST converts a phrase to a pair of random number generator seeds.
c
c  Discussion:
c
c    Uses a phrase (character string) to generate two seeds for the RGN
c    random number generator.
c
c  Modified:
c
c    10 December 2007
c
c  Parameters:
c
c     PHRASE --> Phrase to be used for random number generation
c                         CHARACTER*(*) PHRASE
c
c     SEED1 <-- First seed for RGN generator
c                         INTEGER SEED1
c
c     SEED2 <-- Second seed for RGN generator
c                         INTEGER SEED2
c
c
c                              Note
c
c
c     Trailing blanks are eliminated before the seeds are generated.
c
c     Generated seed values will fall in the range 1..2^30
c     (1..1,073,741,824)
c
c     .. Parameters ..
      CHARACTER*(*) table
      PARAMETER (table='abcdefghijklmnopqrstuvwxyz'//
     +          'ABCDEFGHIJKLMNOPQRSTUVWXYZ'//'0123456789'//
     +          '!@#$%^&*()_+[];:''"<>?,./')
      INTEGER twop30
      PARAMETER (twop30=1073741824)
c     ..
c     .. Scalar Arguments ..
      INTEGER seed1,seed2
      CHARACTER phrase* (*)
c     ..
c     .. Local Scalars ..
      INTEGER i,ichr,j,lphr
c     ..
c     .. Local Arrays ..
      INTEGER shift(0:4),values(5)
c     ..
c     .. External Functions ..
      INTEGER lennob
      EXTERNAL lennob
c     ..
c     .. Intrinsic Functions ..
      INTRINSIC index,mod
c     ..
c     .. Data statements ..
      DATA shift/1,64,4096,262144,16777216/
c     ..
c     .. Executable Statements ..
      seed1 = 1234567890
      seed2 = 123456789
      lphr = lennob(phrase)
      IF (lphr.LT.1) RETURN
      DO 30,i = 1,lphr
          ichr = mod(index(table,phrase(i:i)),64)
          IF (ichr.EQ.0) ichr = 63
          DO 10,j = 1,5
              values(j) = ichr - j
              IF (values(j).LT.1) values(j) = values(j) + 63
   10     CONTINUE
          DO 20,j = 1,5
              seed1 = mod(seed1+shift(j-1)*values(j),twop30)
              seed2 = mod(seed2+shift(j-1)*values(6-j),twop30)
   20     CONTINUE
   30 CONTINUE
      RETURN

      END
      LOGICAL FUNCTION qrgnin()

c*********************************************************************72
c
cc QRGNIN determines whether the random number generator was initialized.
c
c  Discussion:
c
c    A trivial routine to determine whether or not the random
c    number generator has been initialized.  Returns .TRUE. if
c    it has, else .FALSE.
c
c  Modified:
c
c    10 December 2007
c
      LOGICAL qvalue
c     ..
c     .. Local Scalars ..
      LOGICAL qinit
c     ..
c     .. Entry Points ..
      LOGICAL qrgnsn
c     ..
c     .. Save statement ..
      SAVE qinit
c     ..
c     .. Data statements ..
      DATA qinit/.FALSE./
c     ..
c     .. Executable Statements ..
      qrgnin = qinit
      RETURN

      ENTRY qrgnsn(qvalue)

c*********************************************************************72
c
cc QRGNSN records whether the random number generator was initialized.
c
c  Discussion:
c
c    Sets state of whether random number generator is initialized
c    to QVALUE.
c
c    This routine is actually an entry in QRGNIN, hence it is a
c    logical function.  It returns the (meaningless) value .TRUE.
c
c  Modified:
c
c    10 December 2007
c
      qinit = qvalue
      qrgnsn = .TRUE.
      RETURN
      END
      REAL FUNCTION ranf()

c*********************************************************************72
c
cc RANF returns a uniform random number.
c
c  Discussion:
c
c    Returns a random floating point number from a uniform distribution
c    over 0 - 1 (endpoints of this interval are not returned) using the
c    current generator
c
c    This is a transcription from Pascal to Fortran of routine
c    Uniform_01 from the paper
c
c  Modified:
c
c    10 December 2007
c
c  Reference:
c
c    Pierre LEcuyer, Serge Cote,
c    Implementing a Random Number Package with Splitting Facilities,
c    ACM Transactions on Mathematical Software,
c    Volume 17, 1991, pages 98-111.
c
      INTEGER ignlgi
      EXTERNAL ignlgi
c     ..
c     .. Executable Statements ..
c
c     4.656613057E-10 is 1/M1  M1 is set in a data statement in IGNLGI
c      and is currently 2147483563. If M1 changes, change this also.
c
      ranf = ignlgi()*4.656613057E-10
      RETURN

      END
      REAL FUNCTION sdot(n,sx,incx,sy,incy)
      REAL sx(1),sy(1),stemp
      INTEGER i,incx,incy,ix,iy,m,mp1,n

      stemp = 0.0E0
      sdot = 0.0E0
      IF (n.LE.0) RETURN
      IF (incx.EQ.1 .AND. incy.EQ.1) GO TO 20
      ix = 1
      iy = 1
      IF (incx.LT.0) ix = (-n+1)*incx + 1
      IF (incy.LT.0) iy = (-n+1)*incy + 1
      DO 10 i = 1,n
          stemp = stemp + sx(ix)*sy(iy)
          ix = ix + incx
          iy = iy + incy
   10 CONTINUE
      sdot = stemp
      RETURN

   20 m = mod(n,5)
      IF (m.EQ.0) GO TO 40
      DO 30 i = 1,m
          stemp = stemp + sx(i)*sy(i)
   30 CONTINUE
      IF (n.LT.5) GO TO 60
   40 mp1 = m + 1
      DO 50 i = mp1,n,5
          stemp = stemp + sx(i)*sy(i) + sx(i+1)*sy(i+1) +
     +            sx(i+2)*sy(i+2) + sx(i+3)*sy(i+3) + sx(i+4)*sy(i+4)
   50 CONTINUE
   60 sdot = stemp
      RETURN

      END
      SUBROUTINE setall(iseed1,iseed2)

c*********************************************************************72
c
cc SETALL sets the initial seeds of all the generators.
c
c  Discussion:
c
c    Sets the initial seed of generator 1 to ISEED1 and ISEED2. The
c    initial seeds of the other generators are set accordingly, and
c    all generators states are set to these seeds.
c
c    This is a transcription from Pascal to Fortran of routine
c    Set_Initial_Seed from the paper
c
c  Modified:
c
c    10 December 2007
c
c  Reference:
c
c    Pierre LEcuyer, Serge Cote,
c    Implementing a Random Number Package with Splitting Facilities,
c    ACM Transactions on Mathematical Software,
c    Volume 17, 1991, pages 98-111.
c
c  Parameters:
c
c     ISEED1 -> First of two integer seeds
c                                   INTEGER ISEED1
c
c     ISEED2 -> Second of two integer seeds
c                                   INTEGER ISEED1
c
c     .. Parameters ..
      INTEGER numg
      PARAMETER (numg=32)
c     ..
c     .. Scalar Arguments ..
      INTEGER iseed1,iseed2
      LOGICAL qssd
c     ..
c     .. Scalars in Common ..
      INTEGER a1,a1vw,a1w,a2,a2vw,a2w,m1,m2
c     ..
c     .. Arrays in Common ..
      INTEGER cg1(numg),cg2(numg),ig1(numg),ig2(numg),lg1(numg),
     +        lg2(numg)
      LOGICAL qanti(numg)
c     ..
c     .. Local Scalars ..
      INTEGER g,ocgn
      LOGICAL qqssd
c     ..
c     .. External Functions ..
      INTEGER mltmod
      LOGICAL qrgnin
      EXTERNAL mltmod,qrgnin
c     ..
c     .. External Subroutines ..
      EXTERNAL getcgn,initgn,inrgcm,setcgn
c     ..
c     .. Common blocks ..
      COMMON /globe/m1,m2,a1,a2,a1w,a2w,a1vw,a2vw,ig1,ig2,lg1,lg2,cg1,
     +       cg2,qanti
c     ..
c     .. Save statement ..
      SAVE /globe/,qqssd
c     ..
c     .. Data statements ..
      DATA qqssd/.FALSE./
c     ..
c     .. Executable Statements ..
c
c     TELL IGNLGI, THE ACTUAL NUMBER GENERATOR, THAT THIS ROUTINE
c      HAS BEEN CALLED.
c
      qqssd = .TRUE.
      CALL getcgn(ocgn)
c
c     Initialize Common Block if Necessary
c
      IF (.NOT. (qrgnin())) CALL inrgcm()
      ig1(1) = iseed1
      ig2(1) = iseed2
      CALL initgn(-1)
      DO 10,g = 2,numg
          ig1(g) = mltmod(a1vw,ig1(g-1),m1)
          ig2(g) = mltmod(a2vw,ig2(g-1),m2)
          CALL setcgn(g)
          CALL initgn(-1)
   10 CONTINUE
      CALL setcgn(ocgn)
      RETURN

      ENTRY rgnqsd(qssd)

c*********************************************************************72
c
cc RGNQSD queries whether the random number seed was set.
c
c  Discussion:
c
c    Returns (LOGICAL) QSSD as .TRUE. if SETALL has been invoked,
c    otherwise returns .FALSE.
c
      qssd = qqssd

      RETURN
      END
      SUBROUTINE setant(qvalue)

c*********************************************************************72
c
cc SETANT sets the antithetic switch.
c
c  Discussion:
c
c    Sets whether the current generator produces antithetic values.  If
c    X   is  the value  normally returned  from  a uniform [0,1] random
c    number generator then 1  - X is the antithetic  value. If X is the
c    value  normally  returned  from a   uniform  [0,N]  random  number
c    generator then N - 1 - X is the antithetic value.
c
c    All generators are initialized to NOT generate antithetic values.
c
c  Modified:
c
c    10 December 2007
c
c     This is a transcription from Pascal to Fortran of routine
c     Set_Antithetic from the paper
c
c  Reference:
c
c    Pierre LEcuyer, Serge Cote,
c    Implementing a Random Number Package with Splitting Facilities,
c    ACM Transactions on Mathematical Software,
c    Volume 17, 1991, pages 98-111.
c
c  Parameters:
c
c     QVALUE -> .TRUE. if generator G is to generating antithetic
c                    values, otherwise .FALSE.
c                                   LOGICAL QVALUE
c
c     .. Parameters ..
      INTEGER numg
      PARAMETER (numg=32)
c     ..
c     .. Scalar Arguments ..
      LOGICAL qvalue
c     ..
c     .. Scalars in Common ..
      INTEGER a1,a1vw,a1w,a2,a2vw,a2w,m1,m2
c     ..
c     .. Arrays in Common ..
      INTEGER cg1(numg),cg2(numg),ig1(numg),ig2(numg),lg1(numg),
     +        lg2(numg)
      LOGICAL qanti(numg)
c     ..
c     .. Local Scalars ..
      INTEGER g
c     ..
c     .. External Functions ..
      LOGICAL qrgnin
      EXTERNAL qrgnin
c     ..
c     .. External Subroutines ..
      EXTERNAL getcgn
c     ..
c     .. Common blocks ..
      COMMON /globe/m1,m2,a1,a2,a1w,a2w,a1vw,a2vw,ig1,ig2,lg1,lg2,cg1,
     +       cg2,qanti
c     ..
c     .. Save statement ..
      SAVE /globe/
c     ..
c     .. Executable Statements ..
c     Abort unless random number generator initialized
      IF (qrgnin()) GO TO 10
      WRITE (*,*) ' SETANT called before random number generator ',
     +  ' initialized -- abort!'
      STOP ' SETANT called before random number generator initialized'

   10 CALL getcgn(g)
      qanti(g) = qvalue
      RETURN

      END
      SUBROUTINE setgmn(meanv,covm,p,parm)

c*********************************************************************72
c
cc SETGMN sets the routine that generates multivariate normal deviates.
c
c  Discussion:
c
c    Places P, MEANV, and the Cholesky factoriztion of COVM
c    in GENMN.
c
c  Modified:
c
c    10 December 2007
c
c  Parameters:
c
c     MEANV --> Mean vector of multivariate normal distribution.
c                                        REAL MEANV(P)
c
c     COVM   <--> (Input) Covariance   matrix    of  the  multivariate
c                 normal distribution
c                 (Output) Destroyed on output
c                                        REAL COVM(P,P)
c
c     P     --> Dimension of the normal, or length of MEANV.
c                                        INTEGER P
c
c     PARM <-- Array of parameters needed to generate multivariate norma
c                deviates (P, MEANV and Cholesky decomposition of
c                COVM).
c                1 : 1                - P
c                2 : P + 1            - MEANV
c                P+2 : P*(P+3)/2 + 1  - Cholesky decomposition of COVM
c                                             REAL PARM(P*(P+3)/2 + 1)
c
c     .. Scalar Arguments ..
      INTEGER p
c     ..
c     .. Array Arguments ..
      REAL covm(p,p),meanv(p),parm(p* (p+3)/2+1)
c     ..
c     .. Local Scalars ..
      INTEGER i,icount,info,j
c     ..
c     .. External Subroutines ..
      EXTERNAL spofa
c     ..
c     .. Executable Statements ..
c
c
c     TEST THE INPUT
c
      IF (.NOT. (p.LE.0)) GO TO 10
      WRITE (*,*) 'P nonpositive in SETGMN'
      WRITE (*,*) 'Value of P: ',p
      STOP 'P nonpositive in SETGMN'

   10 parm(1) = p
c
c     PUT P AND MEANV INTO PARM
c
      DO 20,i = 2,p + 1
          parm(i) = meanv(i-1)
   20 CONTINUE
c
c      Cholesky decomposition to find A s.t. trans(A)*(A) = COVM
c
      CALL spofa(covm,p,p,info)
      IF (.NOT. (info.NE.0)) GO TO 30
      WRITE (*,*) ' COVM not positive definite in SETGMN'
      STOP ' COVM not positive definite in SETGMN'

   30 icount = p + 1
c
c     PUT UPPER HALF OF A, WHICH IS NOW THE CHOLESKY FACTOR, INTO PARM
c          COVM(1,1) = PARM(P+2)
c          COVM(1,2) = PARM(P+3)
c                    :
c          COVM(1,P) = PARM(2P+1)
c          COVM(2,2) = PARM(2P+2)  ...
c
      DO 50,i = 1,p
          DO 40,j = i,p
              icount = icount + 1
              parm(icount) = covm(i,j)
   40     CONTINUE
   50 CONTINUE

      RETURN
      END
      SUBROUTINE setsd(iseed1,iseed2)

c*********************************************************************72
c
cc SETSD sets the seed of the current random number generator.
c
c  Discussion:
c
c    Resets the initial  seed of  the current  generator to  ISEED1 and
c    ISEED2. The seeds of the other generators remain unchanged.
c
c    This is a transcription from Pascal to Fortran of routine
c    Set_Seed from the paper
c
c  Modified:
c
c    10 December 2007
c
c  Reference:
c
c    Pierre LEcuyer, Serge Cote,
c    Implementing a Random Number Package with Splitting Facilities,
c    ACM Transactions on Mathematical Software,
c    Volume 17, 1991, pages 98-111.
c
c  Parameters:
c
c     ISEED1 -> First integer seed
c                                   INTEGER ISEED1
c
c     ISEED2 -> Second integer seed
c                                   INTEGER ISEED1
c
c     .. Parameters ..
      INTEGER numg
      PARAMETER (numg=32)
c     ..
c     .. Scalar Arguments ..
      INTEGER iseed1,iseed2
c     ..
c     .. Scalars in Common ..
      INTEGER a1,a1vw,a1w,a2,a2vw,a2w,m1,m2
c     ..
c     .. Arrays in Common ..
      INTEGER cg1(numg),cg2(numg),ig1(numg),ig2(numg),lg1(numg),
     +        lg2(numg)
      LOGICAL qanti(numg)
c     ..
c     .. Local Scalars ..
      INTEGER g
c     ..
c     .. External Functions ..
      LOGICAL qrgnin
      EXTERNAL qrgnin
c     ..
c     .. External Subroutines ..
      EXTERNAL getcgn,initgn
c     ..
c     .. Common blocks ..
      COMMON /globe/m1,m2,a1,a2,a1w,a2w,a1vw,a2vw,ig1,ig2,lg1,lg2,cg1,
     +       cg2,qanti
c     ..
c     .. Save statement ..
      SAVE /globe/
c     ..
c     .. Executable Statements ..
c     Abort unless random number generator initialized
      IF (qrgnin()) GO TO 10
      WRITE (*,*) ' SETSD called before random number generator ',
     +  ' initialized -- abort!'
      STOP ' SETSD called before random number generator initialized'

   10 CALL getcgn(g)
      ig1(g) = iseed1
      ig2(g) = iseed2
      CALL initgn(-1)
      RETURN

      END
      REAL FUNCTION sexpo()

c*********************************************************************72
c
cc SEXPO evaluates the standard exponential distribution.
c                                                                      
c  Discussion:                                                                    
c        
c  Modified:
c
c    10 December 2007
c              
c  Reference:
c                                                
c    Joachim Ahrens, Ulrich Dieter,
c    Computer Methods for Sampling From the
c    Exponential and Normal Distributions,
c    Communications of the ACM,
c    Volume 15, Number 10, October 1972, pages 873-882.           
c                                                                      
c     ALL STATEMENT NUMBERS CORRESPOND TO THE STEPS OF ALGORITHM       
c     'SA' IN THE ABOVE PAPER (SLIGHTLY MODIFIED IMPLEMENTATION)       
c                                                                      
c     Modified by Barry W. Brown, Feb 3, 1988 to use RANF instead of   
c     SUNIF.  The argument IR thus goes away.                          
c                                                                      
c
      DIMENSION q(8)
      EQUIVALENCE (q(1),q1)
c
c     Q(N) = SUM(ALOG(2.0)**K/K!)    K=1,..,N ,      THE HIGHEST N
c     (HERE 8) IS DETERMINED BY Q(N)=1.0 WITHIN STANDARD PRECISION
c
      DATA q/.6931472,.9333737,.9888778,.9984959,.9998293,.9999833,
     +     .9999986,.9999999/
c
   10 a = 0.0
      u = ranf()
      GO TO 30

   20 a = a + q1
   30 u = u + u
      IF (u.LE.1.0) GO TO 20
   40 u = u - 1.0
      IF (u.GT.q1) GO TO 60
   50 sexpo = a + u
      RETURN

   60 i = 1
      ustar = ranf()
      umin = ustar
   70 ustar = ranf()
      IF (ustar.LT.umin) umin = ustar
   80 i = i + 1
      IF (u.GT.q(i)) GO TO 70
   90 sexpo = a + umin*q1
      RETURN

      END
      REAL FUNCTION sgamma(a)

c*********************************************************************72
c
cc SGAMMA evaluates the standard Gamma distribution.
c
c  Discussion:
c
c  Modified:
c
c    10 December 2007
c
c  Reference:
c                                                                                                                                           C
c    Joachim Ahrens, Ulrich Dieter,
c    Generating Gamma Variates by a Modified Rejection Technique,<br>
c    Communications of the ACM,<br>
c    Volume 25, Number 1, January 1982, pages 47-54.               
c                                                                      
c     STEP NUMBERS CORRESPOND TO ALGORITHM 'GD' IN THE ABOVE PAPER     
c                                 (STRAIGHTFORWARD IMPLEMENTATION)     
c                                                                      
c     Modified by Barry W. Brown, Feb 3, 1988 to use RANF instead of   
c     SUNIF.  The argument IR thus goes away.                          
c
c  Parameters:
c                                                                      
c               PARAMETER  0.0 < A < 1.0  !                            
c
c     INPUT: A =PARAMETER (MEAN) OF THE STANDARD GAMMA DISTRIBUTION
c     OUTPUT: SGAMMA = SAMPLE FROM THE GAMMA-(A)-DISTRIBUTION
c
c     COEFFICIENTS Q(K) - FOR Q0 = SUM(Q(K)*A**(-K))
c     COEFFICIENTS A(K) - FOR Q = Q0+(T*T/2)*SUM(A(K)*V**K)
c     COEFFICIENTS E(K) - FOR EXP(Q)-1 = SUM(E(K)*Q**K)
c
      DATA q1,q2,q3,q4,q5,q6,q7/.04166669,.02083148,.00801191,.00144121,
     +     -.00007388,.00024511,.00024240/
      DATA a1,a2,a3,a4,a5,a6,a7/.3333333,-.2500030,.2000062,-.1662921,
     +     .1423657,-.1367177,.1233795/
      DATA e1,e2,e3,e4,e5/1.,.4999897,.1668290,.0407753,.0102930/
c
c     PREVIOUS A PRE-SET TO ZERO - AA IS A', AAA IS A"
c     SQRT32 IS THE SQUAREROOT OF 32 = 5.656854249492380
c
      DATA aa/0.0/,aaa/0.0/,sqrt32/5.656854/
c
c     SAVE STATEMENTS
c
      SAVE aa,aaa,s2,s,d,q0,b,si,c
c
      IF (a.EQ.aa) GO TO 10
      IF (a.LT.1.0) GO TO 140
c
c     STEP  1:  RECALCULATIONS OF S2,S,D IF A HAS CHANGED
c
      aa = a
      s2 = a - 0.5
      s = sqrt(s2)
      d = sqrt32 - 12.0*s
c
c     STEP  2:  T=STANDARD NORMAL DEVIATE,
c               X=(S,1/2)-NORMAL DEVIATE.
c               IMMEDIATE ACCEPTANCE (I)
c
   10 t = snorm()
      x = s + 0.5*t
      sgamma = x*x
      IF (t.GE.0.0) RETURN
c
c     STEP  3:  U= 0,1 -UNIFORM SAMPLE. SQUEEZE ACCEPTANCE (S)
c
      u = ranf()
      IF (d*u.LE.t*t*t) RETURN
c
c     STEP  4:  RECALCULATIONS OF Q0,B,SI,C IF NECESSARY
c
      IF (a.EQ.aaa) GO TO 40
      aaa = a
      r = 1.0/a
      q0 = ((((((q7*r+q6)*r+q5)*r+q4)*r+q3)*r+q2)*r+q1)*r
c
c               APPROXIMATION DEPENDING ON SIZE OF PARAMETER A
c               THE CONSTANTS IN THE EXPRESSIONS FOR B, SI AND
c               C WERE ESTABLISHED BY NUMERICAL EXPERIMENTS
c
      IF (a.LE.3.686) GO TO 30
      IF (a.LE.13.022) GO TO 20
c
c               CASE 3:  A .GT. 13.022
c
      b = 1.77
      si = .75
      c = .1515/s
      GO TO 40
c
c               CASE 2:  3.686 .LT. A .LE. 13.022
c
   20 b = 1.654 + .0076*s2
      si = 1.68/s + .275
      c = .062/s + .024
      GO TO 40
c
c               CASE 1:  A .LE. 3.686
c
   30 b = .463 + s + .178*s2
      si = 1.235
      c = .195/s - .079 + .16*s
c
c     STEP  5:  NO QUOTIENT TEST IF X NOT POSITIVE
c
   40 IF (x.LE.0.0) GO TO 70
c
c     STEP  6:  CALCULATION OF V AND QUOTIENT Q
c
      v = t/ (s+s)
      IF (abs(v).LE.0.25) GO TO 50
      q = q0 - s*t + 0.25*t*t + (s2+s2)*alog(1.0+v)
      GO TO 60

   50 q = q0 + 0.5*t*t* ((((((a7*v+a6)*v+a5)*v+a4)*v+a3)*v+a2)*v+a1)*v
c
c     STEP  7:  QUOTIENT ACCEPTANCE (Q)
c
   60 IF (alog(1.0-u).LE.q) RETURN
c
c     STEP  8:  E=STANDARD EXPONENTIAL DEVIATE
c               U= 0,1 -UNIFORM DEVIATE
c               T=(B,SI)-DOUBLE EXPONENTIAL (LAPLACE) SAMPLE
c
   70 e = sexpo()
      u = ranf()
      u = u + u - 1.0
      t = b + sign(si*e,u)
      IF (.NOT. (u.GE.0.0)) GO TO 80
      t = b + si*e
      GO TO 90

   80 t = b - si*e

c
c     STEP  9:  REJECTION IF T .LT. TAU(1) = -.71874483771719
c
   90 IF (t.LT. (-.7187449)) GO TO 70
c
c     STEP 10:  CALCULATION OF V AND QUOTIENT Q
c
      v = t/ (s+s)
      IF (abs(v).LE.0.25) GO TO 100
      q = q0 - s*t + 0.25*t*t + (s2+s2)*alog(1.0+v)
      GO TO 110

  100 q = q0 + 0.5*t*t* ((((((a7*v+a6)*v+a5)*v+a4)*v+a3)*v+a2)*v+a1)*v
c
c     STEP 11:  HAT ACCEPTANCE (H) (IF Q NOT POSITIVE GO TO STEP 8)
c
  110 IF (q.LE.0.0) GO TO 70
      IF (q.LE.0.5) GO TO 120
      w = exp(q) - 1.0
      GO TO 130

  120 w = ((((e5*q+e4)*q+e3)*q+e2)*q+e1)*q
c
c  IF T IS REJECTED, SAMPLE AGAIN AT STEP 8
c
  130 IF (c*abs(u).GT.w*exp(e-0.5*t*t)) GO TO 70
      x = s + 0.5*t
      sgamma = x*x
      RETURN
c
c  ALTERNATE METHOD FOR PARAMETERS A BELOW 1  (.3678794=EXP(-1.))
c
  140 aa = 0.0
      b = 1.0 + .3678794*a
  150 p = b*ranf()
      IF (p.GE.1.0) GO TO 160
      sgamma = exp(alog(p)/a)
      IF (sexpo().LT.sgamma) GO TO 150
      RETURN

  160 sgamma = -alog((b-p)/a)
      IF (sexpo().LT. (1.0-a)*alog(sgamma)) GO TO 150
      RETURN
      END
      REAL FUNCTION snorm()

c*********************************************************************72
c        
cc SNORM evaluates the standard normal distribution.
c
c  Discussion:
c          
c  Modified:
c
c    10 December 2007
c
c  Reference:
c
c    Joachim Ahrens, Ulrich Dieter,
c    Extensions of Forsythe's Method for Random
c    Sampling from the Normal Distribution,
c    Mathematics of Computation,
c    Volume 27, Number 124, October 1973, page 927-937.      
c                                                                      
c     ALL STATEMENT NUMBERS CORRESPOND TO THE STEPS OF ALGORITHM 'FL'  
c     (M=5) IN THE ABOVE PAPER     (SLIGHTLY MODIFIED IMPLEMENTATION)  
c                                                                      
c     Modified by Barry W. Brown, Feb 3, 1988 to use RANF instead of   
c     SUNIF.  The argument IR thus goes away.                          
c                                                                      
      DIMENSION a(32),d(31),t(31),h(31)
c
c     THE DEFINITIONS OF THE CONSTANTS A(K), D(K), T(K) AND
c     H(K) ARE ACCORDING TO THE ABOVEMENTIONED ARTICLE
c
      DATA a/0.0,.3917609E-1,.7841241E-1,.1177699,.1573107,.1970991,
     +     .2372021,.2776904,.3186394,.3601299,.4022501,.4450965,
     +     .4887764,.5334097,.5791322,.6260990,.6744898,.7245144,
     +     .7764218,.8305109,.8871466,.9467818,1.009990,1.077516,
     +     1.150349,1.229859,1.318011,1.417797,1.534121,1.675940,
     +     1.862732,2.153875/
      DATA d/5*0.0,.2636843,.2425085,.2255674,.2116342,.1999243,
     +     .1899108,.1812252,.1736014,.1668419,.1607967,.1553497,
     +     .1504094,.1459026,.1417700,.1379632,.1344418,.1311722,
     +     .1281260,.1252791,.1226109,.1201036,.1177417,.1155119,
     +     .1134023,.1114027,.1095039/
      DATA t/.7673828E-3,.2306870E-2,.3860618E-2,.5438454E-2,
     +     .7050699E-2,.8708396E-2,.1042357E-1,.1220953E-1,.1408125E-1,
     +     .1605579E-1,.1815290E-1,.2039573E-1,.2281177E-1,.2543407E-1,
     +     .2830296E-1,.3146822E-1,.3499233E-1,.3895483E-1,.4345878E-1,
     +     .4864035E-1,.5468334E-1,.6184222E-1,.7047983E-1,.8113195E-1,
     +     .9462444E-1,.1123001,.1364980,.1716886,.2276241,.3304980,
     +     .5847031/
      DATA h/.3920617E-1,.3932705E-1,.3950999E-1,.3975703E-1,
     +     .4007093E-1,.4045533E-1,.4091481E-1,.4145507E-1,.4208311E-1,
     +     .4280748E-1,.4363863E-1,.4458932E-1,.4567523E-1,.4691571E-1,
     +     .4833487E-1,.4996298E-1,.5183859E-1,.5401138E-1,.5654656E-1,
     +     .5953130E-1,.6308489E-1,.6737503E-1,.7264544E-1,.7926471E-1,
     +     .8781922E-1,.9930398E-1,.1155599,.1404344,.1836142,.2790016,
     +     .7010474/
c
   10 u = ranf()
      s = 0.0
      IF (u.GT.0.5) s = 1.0
      u = u + u - s
   20 u = 32.0*u
      i = int(u)
      IF (i.EQ.32) i = 31
      IF (i.EQ.0) GO TO 100
c
c  START CENTER
c
   30 ustar = u - float(i)
      aa = a(i)
   40 IF (ustar.LE.t(i)) GO TO 60
      w = (ustar-t(i))*h(i)
c
c  EXIT   (BOTH CASES)
c
   50 y = aa + w
      snorm = y
      IF (s.EQ.1.0) snorm = -y
      RETURN
c
c                                CENTER CONTINUED
c
   60 u = ranf()
      w = u* (a(i+1)-aa)
      tt = (0.5*w+aa)*w
      GO TO 80

   70 tt = u
      ustar = ranf()
   80 IF (ustar.GT.tt) GO TO 50
   90 u = ranf()
      IF (ustar.GE.u) GO TO 70
      ustar = ranf()
      GO TO 40
c
c                                START TAIL
c
  100 i = 6
      aa = a(32)
      GO TO 120

  110 aa = aa + d(i)
      i = i + 1
  120 u = u + u
      IF (u.LT.1.0) GO TO 110
  130 u = u - 1.0
  140 w = u*d(i)
      tt = (0.5*w+aa)*w
      GO TO 160

  150 tt = u
  160 ustar = ranf()
      IF (ustar.GT.tt) GO TO 50
  170 u = ranf()
      IF (ustar.GE.u) GO TO 150
      u = ranf()
      GO TO 140

      END
      SUBROUTINE spofa(a,lda,n,info)

c*********************************************************************72
c
cc SPOFA factors a real symmetric positive definite matrix.
c
c  Discussion:
c
c     SPOFA IS USUALLY CALLED BY SPOCO, BUT IT CAN BE CALLED
c     DIRECTLY WITH A SAVING IN TIME IF  RCOND  IS NOT NEEDED.
c     (TIME FOR SPOCO) = (1 + 18/N)*(TIME FOR SPOFA) .
c
c  Modified:
c
c    10 December 2007
c
c  Parameters:
c
c     ON ENTRY
c
c        A       REAL(LDA, N)
c                THE SYMMETRIC MATRIX TO BE FACTORED.  ONLY THE
c                DIAGONAL AND UPPER TRIANGLE ARE USED.
c
c        LDA     INTEGER
c                THE LEADING DIMENSION OF THE ARRAY  A .
c
c        N       INTEGER
c                THE ORDER OF THE MATRIX  A .
c
c     ON RETURN
c
c        A       AN UPPER TRIANGULAR MATRIX  R  SO THAT  A = TRANS(R)*R
c                WHERE  TRANS(R)  IS THE TRANSPOSE.
c                THE STRICT LOWER TRIANGLE IS UNALTERED.
c                IF  INFO .NE. 0 , THE FACTORIZATION IS NOT COMPLETE.
c
c        INFO    INTEGER
c                = 0  FOR NORMAL RETURN.
c                = K  SIGNALS AN ERROR CONDITION.  THE LEADING MINOR
c                     OF ORDER  K  IS NOT POSITIVE DEFINITE.
c
c     LINPACK.  THIS VERSION DATED 08/14/78 .
c     CLEVE MOLER, UNIVERSITY OF NEW MEXICO, ARGONNE NATIONAL LAB.
c
      INTEGER lda,n,info
      REAL a(lda,1)
      REAL sdot,t
      REAL s
      INTEGER j,jm1,k
c     BEGIN BLOCK WITH ...EXITS TO 40
c
c
      DO 30 j = 1,n
          info = j
          s = 0.0E0
          jm1 = j - 1
          IF (jm1.LT.1) GO TO 20
          DO 10 k = 1,jm1
              t = a(k,j) - sdot(k-1,a(1,k),1,a(1,j),1)
              t = t/a(k,k)
              a(k,j) = t
              s = s + t*t
   10     CONTINUE
   20     CONTINUE
          s = a(j,j) - s
c     ......EXIT
          IF (s.LE.0.0E0) GO TO 40
          a(j,j) = sqrt(s)
   30 CONTINUE
      info = 0
   40 CONTINUE
      RETURN

      END
      REAL FUNCTION covar(x,y,n)

c*********************************************************************72
c
cc COVAR computes the covariance of two vectors.
c
c  Discussion:
c
c  Modified:
c
c    10 December 2007
c
c  Parameters:
c
c    Input, real X(N), Y(N), the two vectors.
c
c    Input, integer N, the dimension of the two vectors.
c
c    Output, real COVAR, the covariance of the two vectors.
c
      INTEGER n
c     ..
c     .. Array Arguments ..
      REAL x(n),y(n)
c     ..
c     .. Local Scalars ..
      REAL avx,avy,varx,vary,xmax,xmin
      INTEGER i
c     ..
c     .. External Subroutines ..
      EXTERNAL stat
c     ..
c     .. Intrinsic Functions ..
      INTRINSIC real
c     ..
c     .. Executable Statements ..
      CALL stats(x,n,avx,varx,xmin,xmax)
      CALL stats(y,n,avy,vary,xmin,xmax)

      covar = 0.0
      DO i = 1,n
          covar = covar + (x(i)-avx)* (y(i)-avy)
      end do

      covar = covar/real(n-1)

      RETURN
      END
      SUBROUTINE prcomp(p,mean,xcovar,answer)

c*********************************************************************72
c
cc PRCOMP prints covariance information.
c
c  Discussion:
c
c  Modified:
c
c    10 December 2007
c
c  Parameters:
c
      INTEGER p,maxp
      PARAMETER (maxp=10)
      REAL mean(p),xcovar(p,p),rcovar(maxp,maxp)
      REAL answer(1000,maxp)
      REAL rmean(maxp),rvar(maxp)
      INTEGER maxobs
      PARAMETER (maxobs=1000)

      DO 10,i = 1,p
          CALL stats(answer(1,i),maxobs,rmean(i),rvar(i),dum1,dum2)
          WRITE (*,*) ' Variable Number',i
          WRITE (*,*) ' Mean ',mean(i),' Generated ',rmean(i)
          WRITE (*,*) ' Variance ',xcovar(i,i),' Generated',rvar(i)
   10 CONTINUE
      WRITE (*,*) '                   Covariances'
      DO 30,i = 1,p
          DO 20,j = 1,i - 1
              WRITE (*,*) ' I = ',i,' J = ',j
              rcovar(i,j) = covar(answer(1,i),answer(1,j),maxobs)
              WRITE (*,*) ' Covariance ',xcovar(i,j),' Generated ',
     +          rcovar(i,j)
   20     CONTINUE
   30 CONTINUE
      RETURN

      END
      SUBROUTINE setcov(p,var,corr,covar)

c*********************************************************************72
c
cc SETCOV sets a covariance matrix from variance and common correlation.
c
c  Discussion:
c
c     Set covariance matrix from variance and common correlation
c 
c  Modified:
c
c    10 December 2007
c
c  Parameters:
c
      integer p

      real corr
      real covar(p,p)
      integer i
      integer j
      real var(p)

      do i = 1, p
        do  j = 1, p
          if ( i .eq. j ) then
            covar(i,j) = var(i)
          else
            covar(i,j) = corr * sqrt ( var(i) * var(j) )
          end if
        end do
      end do

      return
      end
      subroutine stats(x,n,av,var,xmin,xmax)

c*********************************************************************72
c
cc STATS computes statistics.
c
c  Discussion:
c
c     Computes AVerage and VARiance of array X(N).
c
c  Modified:
c
c    10 December 2007
c
c  Parameters:
c
      REAL av,var,xmax,xmin
      INTEGER n
c     ..
c     .. Array Arguments ..
      REAL x(n)
c     ..
c     .. Local Scalars ..
      REAL sum
      INTEGER i
c     ..
c     .. Intrinsic Functions ..
      INTRINSIC real
c     ..
c     .. Executable Statements ..
      xmin = x(1)
      xmax = x(1)
      sum = 0.0
      DO 10,i = 1,n
          sum = sum + x(i)
          IF (x(i).LT.xmin) xmin = x(i)
          IF (x(i).GT.xmax) xmax = x(i)
   10 CONTINUE
      av = sum/real(n)
      sum = 0.0
      DO 20,i = 1,n
          sum = sum + (x(i)-av)**2
   20 CONTINUE
      var = sum/real(n-1)
      RETURN

      END
      SUBROUTINE trstat(type,parin,av,var)

c*********************************************************************72
c
cc TRSTAT returns the mean and variance for distributions.
c
c  Discussion:
c
c    Returns mean and variance for a number of statistical distribution
c    as a function of their parameters.
c
c  Modified:
c
c    10 December 2007
c
c  Parameters:
c
c     TYPE --> Character string indicating type of distribution
c             'chis' chisquare
c             'ncch' noncentral chisquare
c             'f'    F (variance ratio)
c             'ncf'  noncentral f
c             'unif' uniform
c             'beta' beta distribution
c                         CHARACTER*(4) TYPE
c
c     PARIN --> Array containing parameters of distribution
c              chisquare
c               PARIN(1) is df
c              noncentral chisquare
c               PARIN(1) is df
c               PARIN(2) is noncentrality parameter
c              F (variance ratio)
c               PARIN(1) is df numerator
c               PARIN(2) is df denominator
c              noncentral F
c               PARIN(1) is df numerator
c               PARIN(2) is df denominator
c               PARIN(3) is noncentrality parameter
c              uniform
c               PARIN(1) is LOW bound
c               PARIN(2) is HIGH bound
c              beta
c               PARIN(1) is A
c               PARIN(2) is B
c                         REAL PARIN(*)
c              binonial
c               PARIN(1) is Number of trials
c               PARIN(2) is Prob Event at Each Trial
c              poisson
c               PARIN(1) is Mean
c
c     AV <-- Mean of specified distribution with specified parameters
c                         REAL AV
c
c     VAR <-- Variance of specified distribution with specified paramete
c                         REAL VAR
c
c
c                              Note
c
c
c     AV and Var will be returned -1 if mean or variance is infinite
c
      IMPLICIT INTEGER (i-n),REAL (a-h,o-p,r-z),LOGICAL (q)

      REAL av,var
      CHARACTER type* (4)
c     ..
c     .. Array Arguments ..
      REAL parin(*)
c     ..
c     .. Local Scalars ..
      REAL a,b,range

      IF ('chis' .eq. type ) then

        av = parin(1)
        var = 2.0*parin(1)
        go to 170

      else if ( 'ncch' .eq. type ) then
  
        a = parin(1) + parin(2)
        b = parin(2)/a
        av = a
        var = 2.0*a* (1.0+b)
        go to 170

      end if

      IF (('f').NE. (type)) GO TO 70
      IF (.NOT. (parin(2).LE.2.0001)) GO TO 30
      av = -1.0
      GO TO 40

   30 av = parin(2)/ (parin(2)-2.0)
   40 IF (.NOT. (parin(2).LE.4.0001)) GO TO 50
      var = -1.0
      GO TO 60

   50 var = (2.0*parin(2)**2* (parin(1)+parin(2)-2.0))/
     +      (parin(1)* (parin(2)-2.0)**2* (parin(2)-4.0))
   60 GO TO 170

   70 IF (('ncf').NE. (type)) GO TO 120
      IF (.NOT. (parin(2).LE.2.0001)) GO TO 80
      av = -1.0
      GO TO 90

   80 av = (parin(2)* (parin(1)+parin(3)))/ ((parin(2)-2.0)*parin(1))
   90 IF (.NOT. (parin(2).LE.4.0001)) GO TO 100
      var = -1.0
      GO TO 110

  100 a = (parin(1)+parin(3))**2 + (parin(1)+2.0*parin(3))*
     +    (parin(2)-2.0)
      b = (parin(2)-2.0)**2* (parin(2)-4.0)
      var = 2.0* (parin(2)/parin(1))**2* (a/b)
  110 GO TO 170

  120 IF (('unif').NE. (type)) GO TO 130
      range = parin(2) - parin(1)
      av = parin(1) + range/2.0
      var = range**2/12.0
      GO TO 170

  130 IF (('beta').NE. (type)) GO TO 140
      av = parin(1)/ (parin(1)+parin(2))
      var = (av*parin(2))/ ((parin(1)+parin(2))*
     +      (parin(1)+parin(2)+1.0))
      WRITE (*,*) ' A, B, AV, VAR ',parin(1),parin(2),av,var
      GO TO 170

  140 IF (('bin').NE. (type)) GO TO 150
      av = parin(1)*parin(2)
      var = av* (1.0-parin(2))
      GO TO 170

  150 IF (('pois').NE. (type)) GO TO 160
      av = parin(1)
      var = parin(1)
      GO TO 170

  160 WRITE (*,*) 'Unimplemented type ',type
      STOP 'Unimplemented type in TRSTAT'

  170 RETURN
      END
