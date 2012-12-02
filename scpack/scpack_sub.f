      subroutine gaussj(n, alpha, beta, b, t, w)

c*********************************************************************72
c
c           this routine computes the nodes t(j) and weights
c        w(j) for gauss-jacobi quadrature formulas.
c        these are used when one wishes to approximate
c
c                 integral (from a to b)  f(x) w(x) dx
c
c                              n
c        by                   sum w  f(t )
c                             j=1  j    j
c
c        (w(x) and w(j) have no connection with each other)
c        where w(x) is the weight function
c
c                   w(x) = (1-x)**alpha * (1+x)**beta
c
c        on (-1, 1), alpha, beta .gt. -1.
c
c     input:
c
c        n        the number of points used for the quadrature rule
c        alpha    see above
c        beta     see above
c        b        real scratch array of length n
c
c     output parameters (both double precision arrays of length n)
c
c        t        will contain the desired nodes.
c        w        will contain the desired weights w(j).
c
c     subroutines required: class, imtql2
c
c     reference:
c
c        the routine has been adapted from the more general
c        routine gaussq by golub and welsch.  see
c        golub, g. h., and welsch, j. h., "calculation of gaussian
c        quadrature rules," mathematics of computation 23 (april,
c        1969), pp. 221-230.
c
      double precision b(n), t(n), w(n), muzero, alpha, beta
c
      call class (n, alpha, beta, b, t, muzero)
      w(1) = 1.0d0
      do 105 i = 2, n
  105    w(i) = 0.0d0
      call imtql2 (n, t, b, w, ierr)
      do 110 i = 1, n
  110    w(i) = muzero * w(i) * w(i)
      return
      end
      subroutine class(n, alpha, beta, b, a, muzero)

c*********************************************************************72
c
c           this procedure supplies the coefficients a(j), b(j) of the
c        recurrence relation
c
c             b p (x) = (x - a ) p   (x) - b   p   (x)
c              j j            j   j-1       j-1 j-2
c
c        for the various classical (normalized) orthogonal polynomials,
c        and the zero-th moment
c
c             muzero = integral w(x) dx
c
c        of the given polynomial's weight function w(x).  since the
c        polynomials are orthonormalized, the tridiagonal matrix is
c        guaranteed to be symmetric.
c
      double precision a(n), b(n), muzero, alpha, beta
      double precision abi, a2b2, dgamma, dsqrt, ab
c
      nm1 = n - 1
c
   50 ab = alpha + beta
      abi = 2.0d0 + ab
      muzero = 2.0d0 ** (ab + 1.0d0) * dgamma(alpha + 1.0d0) * dgamma(
     x beta + 1.0d0) / dgamma(abi)
      a(1) = (beta - alpha)/abi
      b(1) = dsqrt(4.0d0*(1.0d0 + alpha)*(1.0d0 + beta)/((abi + 1.0d0)*
     1  abi*abi))
      a2b2 = beta*beta - alpha*alpha
      do 51 i = 2, nm1
         abi = 2.0d0*i + ab
         a(i) = a2b2/((abi - 2.0d0)*abi)
   51    b(i) = dsqrt (4.0d0*i*(i + alpha)*(i + beta)*(i + ab)/
     1   ((abi*abi - 1)*abi*abi))
      abi = 2.0d0*n + ab
      a(n) = a2b2/((abi - 2.0d0)*abi)
      return
      end
      subroutine imtql2(n, d, e, z, ierr)

c*********************************************************************72
c
c     this is a modified version of the eispack routine imtql2.
c     it finds the eigenvalues and first components of the
c     eigenvectors of a symmetric tridiagonal matrix by the implicit ql
c     method.
c
      integer i, j, k, l, m, n, ii, mml, ierr
      real*8 d(n), e(n), z(n), b, c, f, g, p, r, s, machep
      real*8 dsqrt, dabs, dsign
c
c     :::::::::: machep is a machine dependent parameter specifying
c                the relative precision of floating point arithmetic.
c                machep = 16.0d0**(-13) for long form arithmetic
c                on ibm s360 ::::::::::
      data machep/1.d-15/
c     data machep/z3410000000000000/
c
      ierr = 0
      if (n .eq. 1) go to 1001
c
      e(n) = 0.0d0
      do 240 l = 1, n
         j = 0
c     :::::::::: look for small sub-diagonal element ::::::::::
  105    do 110 m = l, n
            if (m .eq. n) go to 120
            if (dabs(e(m)) .le. machep * (dabs(d(m)) + dabs(d(m+1))))
     x         go to 120
  110    continue
c
  120    p = d(l)
         if (m .eq. l) go to 240
         if (j .eq. 30) go to 1000
         j = j + 1
c     :::::::::: form shift ::::::::::
         g = (d(l+1) - p) / (2.0d0 * e(l))
         r = dsqrt(g*g+1.0d0)
         g = d(m) - p + e(l) / (g + dsign(r, g))
         s = 1.0d0
         c = 1.0d0
         p = 0.0d0
         mml = m - l
c
c     :::::::::: for i=m-1 step -1 until l do -- ::::::::::
         do 200 ii = 1, mml
            i = m - ii
            f = s * e(i)
            b = c * e(i)
            if (dabs(f) .lt. dabs(g)) go to 150
            c = g / f
            r = dsqrt(c*c+1.0d0)
            e(i+1) = f * r
            s = 1.0d0 / r
            c = c * s
            go to 160
  150       s = f / g
            r = dsqrt(s*s+1.0d0)
            e(i+1) = g * r
            c = 1.0d0 / r
            s = s * c
  160       g = d(i+1) - p
            r = (d(i) - g) * s + 2.0d0 * c * b
            p = s * r
            d(i+1) = g + p
            g = c * r - b
c     :::::::::: form first component of vector ::::::::::
            f = z(i+1)
            z(i+1) = s * z(i) + c * f
  200       z(i) = c * z(i) - s * f
c
         d(l) = d(l) - p
         e(l) = g
         e(m) = 0.0d0
         go to 105
  240 continue
c
c     :::::::::: order eigenvalues and eigenvectors ::::::::::
      do 300 ii = 2, n
         i = ii - 1
         k = i
         p = d(i)
c
         do 260 j = ii, n
            if (d(j) .ge. p) go to 260
            k = j
            p = d(j)
  260    continue
c
         if (k .eq. i) go to 300
         d(k) = d(i)
         d(i) = p
         p = z(i)
         z(i) = z(k)
         z(k) = p
  300 continue
c
      go to 1001
c     :::::::::: set error -- no convergence to an
c                eigenvalue after 30 iterations ::::::::::
 1000 ierr = l
 1001 return
      end
      function dgamma(x)

c*********************************************************************72
*
* computes the gamma function via a series expansion from
* abramowitz & stegun
*
* l. n. trefethen, 1/13/84
*
      implicit double precision (a-h,o-z)
      dimension c(26)
      data c /
     1   1.00000 00000 000000,  .57721 56649 015329,
     3   -.65587 80715 202538, -.04200 26350 340952,
     5    .16653 86113 822915, -.04219 77345 555443,
     7   -.00962 19715 278770,  .00721 89432 466630,
     9   -.00116 51675 918591, -.00021 52416 741149,
     1    .00012 80502 823882, -.00002 01348 547807,
     3   -.00000 12504 934821,  .00000 11330 272320,
     5   -.00000 02056 338417,  .00000 00061 160950,
     7    .00000 00050 020075, -.00000 00011 812746,
     9    .00000 00001 043427,  .00000 00000 077823,
     1   -.00000 00000 036968,  .00000 00000 005100,
     3   -.00000 00000 000206, -.00000 00000 000054,
     5    .00000 00000 000014,  .00000 00000 000001/
*
* argument reduction:
      fac = 1.d0
      y = x
    1 if (y.le.1.5d0) goto 2
        y = y - 1.d0
        fac = fac*y
        goto 1
    2 if (y.ge.0.5d0) goto 3
        fac = fac/y
        y = y + 1.d0
        goto 2
*
* series:
    3 g = c(26)
      do 4 i = 1,25
        ii = 26-i
    4   g = y*g + c(ii)
      g = y*g
      dgamma = fac/g
      return
      end
      subroutine ode(f,neqn,y,t,tout,relerr,abserr,iflag,work,iwork)

c*********************************************************************72
      implicit real*8(a-h,o-z)
c
c   double precision subroutine ode integrates a system of neqn
c   first order ordinary differential equations of the form:
c             dy(i)/dt = f(t,y(1),y(2),...,y(neqn))
c             y(i) given at  t .
c   the subroutine integrates from  t  to  tout .  on return the
c   parameters in the call list are set for continuing the integration.
c   the user has only to define a new value  tout  and call  ode  again.
c
c   the differential equations are actually solved by a suite of codes
c   de ,  step , and  intrp .  ode  allocates virtual storage in the
c   arrays  work  and  iwork  and calls  de .  de  is a supervisor which
c   directs the solution.  it calls on the routines  step  and  intrp
c   to advance the integration and to interpolate at output points.
c   step  uses a modified divided difference form of the adams pece
c   formulas and local extrapolation.  it adjusts the order and step
c   size to control the local error per unit step in a generalized
c   sense.  normally each call to  step  advances the solution one step
c   in the direction of  tout .  for reasons of efficiency  de
c   integrates beyond  tout  internally, though never beyond
c   t+10*(tout-t), and calls  intrp  to interpolate the solution at
c   tout .  an option is provided to stop the integration at  tout  but
c   it should be used only if it is impossible to continue the
c   integration beyond  tout .
c
c   this code is completely explained and documented in the text,
c   computer solution of ordinary differential equations:  the initial
c   value problem  by l. f. shampine and m. k. gordon.
c
c   the parameters represent:
c      f -- double precision subroutine f(t,y,yp) to evaluate
c                derivatives yp(i)=dy(i)/dt
c      neqn -- number of equations to be integrated (integer*4)
c      y(*) -- solution vector at t                 (real*8)
c      t -- independent variable                    (real*8)
c      tout -- point at which solution is desired   (real*8)
c      relerr,abserr -- relative and absolute error tolerances for local
c           error test (real*8).  at each step the code requires
c             dabs(local error) .le. dabs(y)*relerr + abserr
c           for each component of the local error and solution vectors
c      iflag -- indicates status of integration     (integer*4)
c      work(*)  (real*8)  -- arrays to hold information internal to
c      iwork(*) (integer*4)    which is necessary for subsequent calls
c
c   first call to ode --
c
c   the user must provide storage in his calling program for the arrays
c   in the call list,
c      y(neqn), work(100+21*neqn), iwork(5),
c   declare  f  in an external statement, supply the double precision
c   subroutine f(t,y,yp)  to evaluate
c      dy(i)/dt = yp(i) = f(t,y(1),y(2),...,y(neqn))
c   and initialize the parameters:
c      neqn -- number of equations to be integrated
c      y(*) -- vector of initial conditions
c      t -- starting point of integration
c      tout -- point at which solution is desired
c      relerr,abserr -- relative and absolute local error tolerances
c      iflag -- +1,-1.  indicator to initialize the code.  normal input
c           is +1.  the user should set iflag=-1 only if it is
c           impossible to continue the integration beyond  tout .
c   all parameters except  f ,  neqn  and  tout  may be altered by the
c   code on output so must be variables in the calling program.
c
c   output from  ode  --
c
c      neqn -- unchanged
c      y(*) -- solution at  t
c      t -- last point reached in integration.  normal return has
c           t = tout .
c      tout -- unchanged
c      relerr,abserr -- normal return has tolerances unchanged.  iflag=
c           signals tolerances increased
c      iflag = 2 -- normal return.  integration reached  tout
c            = 3 -- integration did not reach  tout  because error
c                   tolerances too small.  relerr ,  abserr  increased
c                   appropriately for continuing
c            = 4 -- integration did not reach  tout  because more than
c                   500 steps needed
c            = 5 -- integration did not reach  tout  because equations
c                   appear to be stiff
c            = 6 -- invalid input parameters (fatal error)
c           the value of  iflag  is returned negative when the input
c           value is negative and the integration does not reach  tout ,
c           i.e., -3, -4, -5.
c      work(*),iwork(*) -- information generally of no interest to the
c           user but necessary for subsequent calls.
c
c   subsequent calls to  ode --
c
c   subroutine  ode  returns with all information needed to continue
c   the integration.  if the integration reached  tout , the user need
c   only define a new  tout  and call again.  if the integration did not
c   reach  tout  and the user wants to continue, he just calls again.
c   the output value of  iflag  is the appropriate input value for
c   subsequent calls.  the only situation in which it should be altered
c   is to stop the integration internally at the new  tout , i.e.,
c   change output  iflag=2  to input  iflag=-2 .  error tolerances may
c   be changed by the user before continuing.  all other parameters must
c   remain unchanged.
c
c
c*  subroutines  de  and  step  contain machine dependent constants.   *
c*  be sure they are set before using  ode .                           *
c
c
      logical start,phase1,nornd
      dimension y(neqn),work(1),iwork(5)
      external f
      data ialpha,ibeta,isig,iv,iw,ig,iphase,ipsi,ix,ih,ihold,istart,
     1  itold,idelsn/1,13,25,38,50,62,75,76,88,89,90,91,92,93/
      iyy = 100
      iwt = iyy + neqn
      ip = iwt + neqn
      iyp = ip + neqn
      iypout = iyp + neqn
      iphi = iypout + neqn
      if(iabs(iflag) .eq. 1) go to 1
      start = work(istart) .gt. 0.0d0
      phase1 = work(iphase) .gt. 0.0d0
      nornd = iwork(2) .ne. -1
    1 call de(f,neqn,y,t,tout,relerr,abserr,iflag,work(iyy),
     1  work(iwt),work(ip),work(iyp),work(iypout),work(iphi),
     2  work(ialpha),work(ibeta),work(isig),work(iv),work(iw),work(ig),
     3  phase1,work(ipsi),work(ix),work(ih),work(ihold),start,
     4  work(itold),work(idelsn),iwork(1),nornd,iwork(3),iwork(4),
     5  iwork(5))
      work(istart) = -1.0d0
      if(start) work(istart) = 1.0d0
      work(iphase) = -1.0d0
      if(phase1) work(iphase) = 1.0d0
      iwork(2) = -1
      if(nornd) iwork(2) = 1
      return
      end
      subroutine de(f,neqn,y,t,tout,relerr,abserr,iflag,
     1  yy,wt,p,yp,ypout,phi,alpha,beta,sig,v,w,g,phase1,psi,x,h,hold,
     2  start,told,delsgn,ns,nornd,k,kold,isnold)

c*********************************************************************72

      implicit real*8(a-h,o-z)
c
c   ode  merely allocates storage for  de  to relieve the user of the
c   inconvenience of a long call list.  consequently  de  is used as
c   described in the comments for  ode .
c
c   this code is completely explained and documented in the text,
c   computer solution of ordinary differential equations:  the initial
c   value problem  by l. f. shampine and m. k. gordon.
c
      logical stiff,crash,start,phase1,nornd
      dimension y(neqn),yy(neqn),wt(neqn),phi(neqn,16),p(neqn),yp(neqn),
     1  ypout(neqn),psi(12),alpha(12),beta(12),sig(13),v(12),w(12),g(13)
      external f
c
c
c*  the only machine dependent constant is based on the machine unit   *
c*  roundoff error  u  which is the smallest positive number such that *
c*  1.0+u .gt. 1.0 .  u  must be calculated and  fouru=4.0*u  inserted *
c*  in the following data statement before using  de .  the routine    *
c*  machin  calculates  u .  fouru  and  twou=2.0*u  must also be      *
c*  inserted in subroutine  step  before calling  de .                 *
      data fouru/.888d-15/
c
c
c   the constant  maxnum  is the maximum number of steps allowed in one
c   call to  de .  the user may change this limit by altering the
c   following statement
      data maxnum/500/
c
c            ***            ***            ***
c   test for improper parameters
c
      if(neqn .lt. 1) go to 10
      if(t .eq. tout) go to 10
      if(relerr .lt. 0.0d0  .or.  abserr .lt. 0.0d0) go to 10
      eps = dmax1(relerr,abserr)
      if(eps .le. 0.0d0) go to 10
      if(iflag .eq. 0) go to 10
      isn = isign(1,iflag)
      iflag = iabs(iflag)
      if(iflag .eq. 1) go to 20
      if(t .ne. told) go to 10
      if(iflag .ge. 2  .and.  iflag .le. 5) go to 20
   10 iflag = 6
      return
c
c   on each call set interval of integration and counter for number of
c   steps.  adjust input error tolerances to define weight vector for
c   subroutine  step
c
   20 del = tout - t
      absdel = dabs(del)
      tend = t + 10.0d0*del
      if(isn .lt. 0) tend = tout
      nostep = 0
      kle4 = 0
      stiff = .false.
      releps = relerr/eps
      abseps = abserr/eps
      if(iflag .eq. 1) go to 30
      if(isnold .lt. 0) go to 30
      if(delsgn*del .gt. 0.0d0) go to 50
c
c   on start and restart also set work variables x and yy(*), store the
c   direction of integration and initialize the step size
c
   30 start = .true.
      x = t
      do 40 l = 1,neqn
   40   yy(l) = y(l)
      delsgn = dsign(1.0d0,del)
      h = dsign(dmax1(dabs(tout-x),fouru*dabs(x)),tout-x)
c
c   if already past output point, interpolate and return
c
   50 if(dabs(x-t) .lt. absdel) go to 60
      call intrp(x,yy,tout,y,ypout,neqn,kold,phi,psi)
      iflag = 2
      t = tout
      told = t
      isnold = isn
      return
c
c   if cannot go past output point and sufficiently close,
c   extrapolate and return
c
   60 if(isn .gt. 0  .or.  dabs(tout-x) .ge. fouru*dabs(x)) go to 80
      h = tout - x
      call f(x,yy,yp)
      do 70 l = 1,neqn
   70   y(l) = yy(l) + h*yp(l)
      iflag = 2
      t = tout
      told = t
      isnold = isn
      return
c
c   test for too many steps
c
   80 if(nostep .lt. maxnum) go to 100
      iflag = isn*4
      if(stiff) iflag = isn*5
      do 90 l = 1,neqn
   90   y(l) = yy(l)
      t = x
      told = t
      isnold = 1
      return
c
c   limit step size, set weight vector and take a step
c
  100 h = dsign(dmin1(dabs(h),dabs(tend-x)),h)
      do 110 l = 1,neqn
  110   wt(l) = releps*dabs(yy(l)) + abseps
      call step(x,yy,f,neqn,h,eps,wt,start,
     1  hold,k,kold,crash,phi,p,yp,psi,
     2  alpha,beta,sig,v,w,g,phase1,ns,nornd)
c
c   test for tolerances too small
c
      if(.not.crash) go to 130
      iflag = isn*3
      relerr = eps*releps
      abserr = eps*abseps
      do 120 l = 1,neqn
  120   y(l) = yy(l)
      t = x
      told = t
      isnold = 1
      return
c
c   augment counter on number of steps and test for stiffness
c
  130 nostep = nostep + 1
      kle4 = kle4 + 1
      if(kold .gt. 4) kle4 = 0
      if(kle4 .ge. 50) stiff = .true.
      go to 50
      end
      subroutine step(x,y,f,neqn,h,eps,wt,start,
     1  hold,k,kold,crash,phi,p,yp,psi,
     2  alpha,beta,sig,v,w,g,phase1,ns,nornd)

c*********************************************************************72

      implicit real*8(a-h,o-z)
c
c   double precision subroutine step
c   integrates a system of first order ordinary
c   differential equations one step, normally from x to x+h, using a
c   modified divided difference form of the adams pece formulas.  local
c   extrapolation is used to improve absolute stability and accuracy.
c   the code adjusts its order and step size to control the local error
c   per unit step in a generalized sense.  special devices are included
c   to control roundoff error and to detect when the user is requesting
c   too much accuracy.
c
c   this code is completely explained and documented in the text,
c   computer solution of ordinary differential equations:  the initial
c   value problem  by l. f. shampine and m. k. gordon.
c
c
c   the parameters represent:
c      x -- independent variable             (real*8)
c      y(*) -- solution vector at x          (real*8)
c      yp(*) -- derivative of solution vector at  x  after successful
c           step                             (real*8)
c      neqn -- number of equations to be integrated (integer*4)
c      h -- appropriate step size for next step.  normally determined by
c           code                             (real*8)
c      eps -- local error tolerance.  must be variable  (real*8)
c      wt(*) -- vector of weights for error criterion   (real*8)
c      start -- logical variable set .true. for first step,  .false.
c           otherwise                        (logical*4)
c      hold -- step size used for last successful step  (real*8)
c      k -- appropriate order for next step (determined by code)
c      kold -- order used for last successful step
c      crash -- logical variable set .true. when no step can be taken,
c           .false. otherwise.
c   the arrays  phi, psi  are required for the interpolation subroutine
c   intrp.  the array p is internal to the code.  all are real*8
c
c   input to  step
c
c      first call --
c
c   the user must provide storage in his driver program for all arrays
c   in the call list, namely
c
c     dimension y(neqn),wt(neqn),phi(neqn,16),p(neqn),yp(neqn),psi(12)
c
c   the user must also declare  start  and  crash  logical variables
c   and  f  an external subroutine, supply the subroutine  f(x,y,yp)
c   to evaluate
c      dy(i)/dx = yp(i) = f(x,y(1),y(2),...,y(neqn))
c   and initialize only the following parameters:
c      x -- initial value of the independent variable
c      y(*) -- vector of initial values of dependent variables
c      neqn -- number of equations to be integrated
c      h -- nominal step size indicating direction of integration
c           and maximum size of step.  must be variable
c      eps -- local error tolerance per step.  must be variable
c      wt(*) -- vector of non-zero weights for error criterion
c      start -- .true.
c
c   step  requires the l2 norm of the vector with components
c   local error(l)/wt(l)  be less than  eps  for a successful step.  the
c   array  wt  allows the user to specify an error test appropriate
c   for his problem.  for example,
c      wt(l) = 1.0  specifies absolute error,
c            = dabs(y(l))  error relative to the most recent value of
c                 the l-th component of the solution,
c            = dabs(yp(l))  error relative to the most recent value of
c                 the l-th component of the derivative,
c            = dmax1(wt(l),dabs(y(l)))  error relative to the largest
c                 magnitude of l-th component obtained so far,
c            = dabs(y(l))*relerr/eps + abserr/eps  specifies a mixed
c                 relative-absolute test where  relerr  is relative
c                 error,  abserr  is absolute error and  eps =
c                 dmax1(relerr,abserr) .
c
c      subsequent calls --
c
c   subroutine  step  is designed so that all information needed to
c   continue the integration, including the step size  h  and the order
c   k , is returned with each step.  with the exception of the step
c   size, the error tolerance, and the weights, none of the parameters
c   should be altered.  the array  wt  must be updated after each step
c   to maintain relative error tests like those above.  normally the
c   integration is continued just beyond the desired endpoint and the
c   solution interpolated there with subroutine  intrp .  if it is
c   impossible to integrate beyond the endpoint, the step size may be
c   reduced to hit the endpoint since the code will not take a step
c   larger than the  h  input.  changing the direction of integration,
c   i.e., the sign of  h , requires the user set  start = .true. before
c   calling  step  again.  this is the only situation in which  start
c   should be altered.
c
c   output from  step
c
c      successful step --
c
c   the subroutine returns after each successful step with  start  and
c   crash  set .false. .  x  represents the independent variable
c   advanced one step of length  hold  from its value on input and  y
c   the solution vector at the new value of  x .  all other parameters
c   represent information corresponding to the new  x  needed to
c   continue the integration.
c
c      unsuccessful step --
c
c   when the error tolerance is too small for the machine precision,
c   the subroutine returns without taking a step and  crash = .true. .
c   an appropriate step size and error tolerance for continuing are
c   estimated and all other information is restored as upon input
c   before returning.  to continue with the larger tolerance, the user
c   just calls the code again.  a restart is neither required nor
c   desirable.
c
      logical start,crash,phase1,nornd
      dimension y(neqn),wt(neqn),phi(neqn,16),p(neqn),yp(neqn),psi(12)
      dimension alpha(12),beta(12),sig(13),w(12),v(12),g(13),
     1  gstr(13),two(13)
      external f
c
c*  the only machine dependent constants are based on the machine unit *
c*  roundoff error  u  which is the smallest positive number such that *
c*  1.0+u .gt. 1.0  .  the user must calculate  u  and insert          *
c*  twou=2.0*u  and  fouru=4.0*u  in the data statement before calling *
c*  the code.  the routine  machin  calculates  u .                    *
c*
      data twou,fouru/.444d-15,.888d-15/

      data two/2.0d0,4.0d0,8.0d0,16.0d0,32.0d0,64.0d0,128.0d0,256.0d0,
     1  512.0d0,1024.0d0,2048.0d0,4096.0d0,8192.0d0/
      data gstr/0.500d0,0.0833d0,0.0417d0,0.0264d0,0.0188d0,0.0143d0,
     1  0.0114d0,0.00936d0,0.00789d0,0.00679d0,0.00592d0,0.00524d0,
     2  0.00468d0/
c     data g(1),g(2)/1.0,0.5/,sig(1)/1.0/
c
c
c       ***     begin block 0     ***
c   check if step size or error tolerance is too small for machine
c   precision.  if first step, initialize phi array and estimate a
c   starting step size.
c                   ***
c
c   if step size is too small, determine an acceptable one
c
      crash = .true.
      if(dabs(h) .ge. fouru*dabs(x)) go to 5
      h = dsign(fouru*dabs(x),h)
      return
    5 p5eps = 0.5d0*eps
c
c   if error tolerance is too small, increase it to an acceptable value
c
      round = 0.0d0
      do 10 l = 1,neqn
   10   round = round + (y(l)/wt(l))**2
      round = twou*dsqrt(round)
      if(p5eps .ge. round) go to 15
      eps = 2.0*round*(1.0d0 + fouru)
      return
   15 crash = .false.
      g(1)=1.0d0
      g(2)=0.5d0
      sig(1)=1.0d0
      if(.not.start) go to 99
c
c   initialize.  compute appropriate step size for first step
c
      call f(x,y,yp)
      sum = 0.0d0
      do 20 l = 1,neqn
        phi(l,1) = yp(l)
        phi(l,2) = 0.0d0
   20   sum = sum + (yp(l)/wt(l))**2
      sum = dsqrt(sum)
      absh = dabs(h)
      if(eps .lt. 16.0d0*sum*h*h) absh = 0.25d0*dsqrt(eps/sum)
      h = dsign(dmax1(absh,fouru*dabs(x)),h)
      hold = 0.0d0
      k = 1
      kold = 0
      start = .false.
      phase1 = .true.
      nornd = .true.
      if(p5eps .gt. 100.0d0*round) go to 99
      nornd = .false.
      do 25 l = 1,neqn
   25   phi(l,15) = 0.0d0
   99 ifail = 0
c       ***     end block 0     ***
c
c       ***     begin block 1     ***
c   compute coefficients of formulas for this step.  avoid computing
c   those quantities not changed when step size is not changed.
c                   ***
c
  100 kp1 = k+1
      kp2 = k+2
      km1 = k-1
      km2 = k-2
c
c   ns is the number of steps taken with size h, including the current
c   one.  when k.lt.ns, no coefficients change
c
      if(h .ne. hold) ns = 0
      if(ns.le.kold)   ns=ns+1
      nsp1 = ns+1
      if (k .lt. ns) go to 199
c
c   compute those components of alpha(*),beta(*),psi(*),sig(*) which
c   are changed
c
      beta(ns) = 1.0d0
      realns = ns
      alpha(ns) = 1.0d0/realns
      temp1 = h*realns
      sig(nsp1) = 1.0d0
      if(k .lt. nsp1) go to 110
      do 105 i = nsp1,k
        im1 = i-1
        temp2 = psi(im1)
        psi(im1) = temp1
        beta(i) = beta(im1)*psi(im1)/temp2
        temp1 = temp2 + h
        alpha(i) = h/temp1
        reali = i
  105   sig(i+1) = reali*alpha(i)*sig(i)
  110 psi(k) = temp1
c
c   compute coefficients g(*)
c
c   initialize v(*) and set w(*).  g(2) is set in data statement
c
      if(ns .gt. 1) go to 120
      do 115 iq = 1,k
        temp3 = iq*(iq+1)
        v(iq) = 1.0d0/temp3
  115   w(iq) = v(iq)
      go to 140
c
c   if order was raised, update diagonal part of v(*)
c
  120 if(k .le. kold) go to 130
      temp4 = k*kp1
      v(k) = 1.0d0/temp4
      nsm2 = ns-2
      if(nsm2 .lt. 1) go to 130
      do 125 j = 1,nsm2
        i = k-j
  125   v(i) = v(i) - alpha(j+1)*v(i+1)
c
c   update v(*) and set w(*)
c
  130 limit1 = kp1 - ns
      temp5 = alpha(ns)
      do 135 iq = 1,limit1
        v(iq) = v(iq) - temp5*v(iq+1)
  135   w(iq) = v(iq)
      g(nsp1) = w(1)
c
c   compute the g(*) in the work vector w(*)
c
  140 nsp2 = ns + 2
      if(kp1 .lt. nsp2) go to 199
      do 150 i = nsp2,kp1
        limit2 = kp2 - i
        temp6 = alpha(i-1)
        do 145 iq = 1,limit2
  145     w(iq) = w(iq) - temp6*w(iq+1)
  150   g(i) = w(1)
  199   continue
c       ***     end block 1     ***
c
c       ***     begin block 2     ***
c   predict a solution p(*), evaluate derivatives using predicted
c   solution, estimate local error at order k and errors at orders k,
c   k-1, k-2 as if constant step size were used.
c                   ***
c
c   change phi to phi star
c
      if(k .lt. nsp1) go to 215
      do 210 i = nsp1,k
        temp1 = beta(i)
        do 205 l = 1,neqn
  205     phi(l,i) = temp1*phi(l,i)
  210   continue
c
c   predict solution and differences
c
  215 do 220 l = 1,neqn
        phi(l,kp2) = phi(l,kp1)
        phi(l,kp1) = 0.0d0
  220   p(l) = 0.0d0
      do 230 j = 1,k
        i = kp1 - j
        ip1 = i+1
        temp2 = g(i)
        do 225 l = 1,neqn
          p(l) = p(l) + temp2*phi(l,i)
  225     phi(l,i) = phi(l,i) + phi(l,ip1)
  230   continue
      if(nornd) go to 240
      do 235 l = 1,neqn
        tau = h*p(l) - phi(l,15)
        p(l) = y(l) + tau
  235   phi(l,16) = (p(l) - y(l)) - tau
      go to 250
  240 do 245 l = 1,neqn
  245   p(l) = y(l) + h*p(l)
  250 xold = x
      x = x + h
      absh = dabs(h)
      call f(x,p,yp)
c
c   estimate errors at orders k,k-1,k-2
c
      erkm2 = 0.0d0
      erkm1 = 0.0d0
      erk = 0.0d0
      do 265 l = 1,neqn
        temp3 = 1.0d0/wt(l)
        temp4 = yp(l) - phi(l,1)
        if(km2)265,260,255
  255   erkm2 = erkm2 + ((phi(l,km1)+temp4)*temp3)**2
  260   erkm1 = erkm1 + ((phi(l,k)+temp4)*temp3)**2
  265   erk = erk + (temp4*temp3)**2
      if(km2)280,275,270
  270 erkm2 = absh*sig(km1)*gstr(km2)*dsqrt(erkm2)
  275 erkm1 = absh*sig(k)*gstr(km1)*dsqrt(erkm1)
  280 temp5 = absh*dsqrt(erk)
      err = temp5*(g(k)-g(kp1))
      erk = temp5*sig(kp1)*gstr(k)
      knew = k
c
c   test if order should be lowered
c
      if(km2)299,290,285
  285 if(dmax1(erkm1,erkm2) .le. erk) knew = km1
      go to 299
  290 if(erkm1 .le. 0.5d0*erk) knew = km1
c
c   test if step successful
c
  299 if(err .le. eps) go to 400
c       ***     end block 2     ***
c
c       ***     begin block 3     ***
c   the step is unsuccessful.  restore  x, phi(*,*), psi(*) .
c   if third consecutive failure, set order to one.  if step fails more
c   than three times, consider an optimal step size.  double error
c   tolerance and return if estimated step size is too small for machine
c   precision.
c                   ***
c
c   restore x, phi(*,*) and psi(*)
c
      phase1 = .false.
      x = xold
      do 310 i = 1,k
        temp1 = 1.0d0/beta(i)
        ip1 = i+1
        do 305 l = 1,neqn
  305     phi(l,i) = temp1*(phi(l,i) - phi(l,ip1))
  310   continue
      if(k .lt. 2) go to 320
      do 315 i = 2,k
  315   psi(i-1) = psi(i) - h
c
c   on third failure, set order to one.  thereafter, use optimal step
c   size
c
  320 ifail = ifail + 1
      temp2 = 0.5d0
      if(ifail - 3) 335,330,325
  325 if(p5eps .lt. 0.25d0*erk) temp2 = dsqrt(p5eps/erk)
  330 knew = 1
  335 h = temp2*h
      k = knew
      if(dabs(h) .ge. fouru*dabs(x)) go to 340
      crash = .true.
      h = dsign(fouru*dabs(x),h)
      eps = eps + eps
      return
  340 go to 100
c       ***     end block 3     ***
c
c       ***     begin block 4     ***
c   the step is successful.  correct the predicted solution, evaluate
c   the derivatives using the corrected solution and update the
c   differences.  determine best order and step size for next step.
c                   ***
  400 kold = k
      hold = h
c
c   correct and evaluate
c
      temp1 = h*g(kp1)
      if(nornd) go to 410
      do 405 l = 1,neqn
        rho = temp1*(yp(l) - phi(l,1)) - phi(l,16)
        y(l) = p(l) + rho
  405   phi(l,15) = (y(l) - p(l)) - rho
      go to 420
  410 do 415 l = 1,neqn
  415   y(l) = p(l) + temp1*(yp(l) - phi(l,1))
  420 call f(x,y,yp)
c
c   update differences for next step
c
      do 425 l = 1,neqn
        phi(l,kp1) = yp(l) - phi(l,1)
  425   phi(l,kp2) = phi(l,kp1) - phi(l,kp2)
      do 435 i = 1,k
        do 430 l = 1,neqn
  430     phi(l,i) = phi(l,i) + phi(l,kp1)
  435   continue
c
c   estimate error at order k+1 unless:
c     in first phase when always raise order,
c     already decided to lower order,
c     step size not constant so estimate unreliable
c
      erkp1 = 0.0d0
      if(knew .eq. km1  .or.  k .eq. 12) phase1 = .false.
      if(phase1) go to 450
      if(knew .eq. km1) go to 455
      if(kp1 .gt. ns) go to 460
      do 440 l = 1,neqn
  440   erkp1 = erkp1 + (phi(l,kp2)/wt(l))**2
      erkp1 = absh*gstr(kp1)*dsqrt(erkp1)
c
c   using estimated error at order k+1, determine appropriate order
c   for next step
c
      if(k .gt. 1) go to 445
      if(erkp1 .ge. 0.5d0*erk) go to 460
      go to 450
  445 if(erkm1 .le. dmin1(erk,erkp1)) go to 455
      if(erkp1 .ge. erk  .or.  k .eq. 12) go to 460
c
c   here erkp1 .lt. erk .lt. dmax1(erkm1,erkm2) else order would have
c   been lowered in block 2.  thus order is to be raised
c
c   raise order
c
  450 k = kp1
      erk = erkp1
      go to 460
c
c   lower order
c
  455 k = km1
      erk = erkm1
c
c   with new order determine appropriate step size for next step
c
  460 hnew = h + h
      if(phase1) go to 465
      if(p5eps .ge. erk*two(k+1)) go to 465
      hnew = h
      if(p5eps .ge. erk) go to 465
      temp2 = k+1
      r = (p5eps/erk)**(1.0d0/temp2)
      hnew = absh*dmax1(0.5d0,dmin1(0.9d0,r))
      hnew = dsign(dmax1(hnew,fouru*dabs(x)),h)
  465 h = hnew
      return
c       ***     end block 4     ***
      end
      subroutine intrp(x,y,xout,yout,ypout,neqn,kold,phi,psi)

c*********************************************************************72

      implicit real*8(a-h,o-z)
c
c   the methods in subroutine  step  approximate the solution near  x
c   by a polynomial.  subroutine  intrp  approximates the solution at
c   xout  by evaluating the polynomial there.  information defining this
c   polynomial is passed from  step  so  intrp  cannot be used alone.
c
c   this code is completely explained and documented in the text,
c   computer solution of ordinary differential equations:  the initial
c   value problem  by l. f. shampine and m. k. gordon.
c
c   input to intrp --
c
c   all floating point variables are double precision
c   the user provides storage in the calling program for the arrays in
c   the call list
       dimension y(neqn),yout(neqn),ypout(neqn),phi(neqn,16),psi(12)
c   and defines
c      xout -- point at which solution is desired.
c   the remaining parameters are defined in  step  and passed to  intrp
c   from that subroutine
c
c   output from  intrp --
c
c      yout(*) -- solution at  xout
c      ypout(*) -- derivative of solution at  xout
c   the remaining parameters are returned unaltered from their input
c   values.  integration with  step  may be continued.
c
      dimension g(13),w(13),rho(13)
      data g(1)/1.0d0/,rho(1)/1.0d0/
c
      hi = xout - x
      ki = kold + 1
      kip1 = ki + 1
c
c   initialize w(*) for computing g(*)
c
      do 5 i = 1,ki
        temp1 = i
    5   w(i) = 1.0d0/temp1
      term = 0.0d0
c
c   compute g(*)
c
      do 15 j = 2,ki
        jm1 = j - 1
        psijm1 = psi(jm1)
        gamma = (hi + term)/psijm1
        eta = hi/psijm1
        limit1 = kip1 - j
        do 10 i = 1,limit1
   10     w(i) = gamma*w(i) - eta*w(i+1)
        g(j) = w(1)
        rho(j) = gamma*rho(jm1)
   15   term = psijm1
c
c   interpolate
c
      do 20 l = 1,neqn
        ypout(l) = 0.0d0
   20   yout(l) = 0.0d0
      do 30 j = 1,ki
        i = kip1 - j
        temp2 = g(i)
        temp3 = rho(i)
        do 25 l = 1,neqn
          yout(l) = yout(l) + temp2*phi(l,i)
   25     ypout(l) = ypout(l) + temp3*phi(l,i)
   30   continue
      do 35 l = 1,neqn
   35   yout(l) = y(l) + hi*yout(l)
      return
      end
      subroutine ns01a (n,x,f,ajinv,dstep,dmax,acc,maxfun,iprint,w,
     2                 calfun)

c*********************************************************************72
c
c    this subroutine solves a system of nonlinear algebraic equations o
c    the form f(x) = 0, where f and x are vectors of length n.  the
c    function f should have continuous first derivatives, but there is
c    no need to calculate these, because the algorithm automatically
c    computes finite difference approximations.  the method used in
c    seeking a solution is a compromise between newton's method and
c    steepest descent.  the iterative nature of the algorithm requires
c    that an initial estimate of the solution be supplied, but this can
c    be fairly poor and the process still will usually converge.  this
c    is because the algorithm tends to  take steepest descent steps whe
c    it is far from the solution, switching to the newton direction for
c    rapid convergence as it nears the solution.
c
c    program author: m. j. d. powell
c
c    documentation : technical report aere-r.5947
c                    harwell, england
c                    november, 1968
c
c    the variables in the calling sequence are defined as follows:
c
c    n        number of variables (equations).  must have n>1.
c
c    x        vector of length n containing the variables of the
c             equations.  on entry this should contain the initial
c             estimate of the solution.  on return it will contain the
c             best estimate of the solution.
c
c    f        vector of length n used as working space.  on return it
c             contains the values of the functions at the solution x.
c
c    ajinv    two dimensional array (20 by 20) used as working space.
c             on return it contains an approximation to the inverse of
c             the jacobian matrix of the system at the solution x.
c
c    dstep    step length for use in approximating first derivatives of
c             the functions by differences between function values.
c
c    dmax     generous estimate of the euclidean distance of the
c             solution from the initial estimate.  used to limit maximu
c             size of a single correction step.  must have dmax > dstep
c
c    acc      accuracy requirement.  the iterative process is assumed t
c             have converged and a normal return takes place when an x
c             is found for which the euclidean norm of f(x) is less tha
c             acc.
c
c    maxfun   maximum number of function evaluations (calls to calfun)
c             to be allowed.  (generally 10*n is a sufficient number of
c             evaluations, although for some problems more may be
c             needed.)
c
c    iprint   flag which controls printing.  if iprint=0 there is little
c             printout except error messages.  if iprint=1 then each
c             time calfun is called the values of x and f are printed.
c
c    w        vector of length n*(2*n+5) which is used as working space
c
c    calfun   the subroutine defining the system of nonlinear equations
c             it must have the form
c
c                     subroutine calfun(n,x,f)
c
c             where n, x, and f have the same meaning as given above.
c             calfun should set the components of f to be the function
c             values at x.  must be declared external in the calling
c             program.
c
c    notes:
c
c    there are four possible error exits, each of which is accompanied
c    by an appropriate diagnostic message.  these occur when
c
c          1. the gradient of f at x becomes so small that it is
c             predicted that steps much larger than dmax are needed.
c
c          2. the number of calls to calfun exceeds maxfun.
c
c          3. f(x) fails to decrease for n+4 consecutive iterations whe
c             a decrease is predicted, and x is within dstep of the mos
c             successful vector of variables.
c
c          4. f(x) fails to decrease even though a new (therefore,
c             presumably reliable) jacobian matrix has just been
c             obtained by finite differences, and x is within distance
c             dstep from the point where the jacobian was calculated.
c
c    these error conditions can be caused by a number of factors, such
c    as too small values for maxfun or dmax, or very poor values for
c    dstep or the initial guess for x.  in addition, the problem itself
c    may be so badly conditioned that rounding errors in the
c    computations make finding a solution impossible, or the equations
c    may be inconsistent so that no solution exists.
c
c    another important fact to note is that only one step size, dstep,
c    is used for all the variables.  this implies that all the variable
c    are expected to be of roughly the same order of magnitude.  scalin
c    may be required to accomplish this, but this enhances the numerica
c    stability of the problem anyway and is definitely advisable.
c
c    note that if the order of the system is greater than 20, then ajin
c    (as well as several other arrays not appearing in the calling
c    sequence) must have their dimensions increased.
c
c    matrix inversion is performed by the linpack routines
c    dgefa followed by dgedi.  these require in addition the
c    linpack basic linear algebra subroutines daxpy, dscal,
c    idamax, and dswap.
c
      implicit real*8 (a-h,o-z)
      double precision  x(1),f(1),lawork(20),ajinv(20,20),w(1)
      double precision det(2)
      integer ipvt(20)
      external calfun
c     set various parameters
      maxc=0
c     'maxc' counts the number of calls of calfun
      nt=n+4
      ntest=nt
c     'nt' and 'ntest' cause an error return if f(x) does not decrease
      dtest=dfloat(n+n)-0.5d0
c     'dtest' is used to maintain linear independence
      nx=n*n
      nf=nx+n
      nw=nf+n
      mw=nw+n
      ndc=mw+n
      nd=ndc+n
c     these parameters separate the working space array w
      fmin=0.0d0
c     usually 'fmin' is the least calculated value of f(x),
c     and the best x is in w(nx+1) to w(nx+n)
      dd=0.0d0
c     usually dd is the square of the current step length
      dss=dstep*dstep
      dm=dmax*dmax
      dmm=4.0d0*dm
      is=5
c     'is' controls a 'go to' statement following a call of calfun
      tinc=1.0d0
c     'tinc' is used in the criterion to increase the step length
c     start a new page for printing
      if (iprint) 1,1,85
   85 print 86
   86 format (1h1)
c     call the subroutine calfun
    1 maxc=maxc+1
      call calfun (n,x,f)
c     test for convergence
      fsq=0.0d0
      do 2 i=1,n
      fsq=fsq+f(i)*f(i)
    2 continue
      if (fsq .gt. acc**2) go to 4
c     provide printing of final solution if requested
    3 if (iprint) 5,999,6
    6 print 7,maxc
    7 format (/5x,'the final ns01a solution required',
     1i5,' calls of calfun, and is')
      print 8,(i,x(i),f(i),i=1,n)
    8 format (4x,'i',8x,'x(i)',13x,'f(i)'/(i5,2e17.8))
  999 print 9,maxc,fsq
    9 format ('     the sum-of-squares error at step',i4,' is',e11.4)
    5 return
c     test for error return because f(x) does not decrease
    4 go to (10,11,11,10,11),is
   10 if (fsq-fmin) 15,20,20
   20 if (dd-dss) 12,12,11
   12 ntest=ntest-1
      if (ntest) 13,14,11
   14 print 16,nt
   16 format (///5x,'error return from ns01a because',i5,
     1' calls of calfun failed to improve the residuals')
   17 do 18 i=1,n
      x(i)=w(nx+i)
      f(i)=w(nf+i)
   18 continue
      fsq=fmin
      go to 3
c     error return because a new jacobian is unsuccessful
   13 print 19
   19 format (///5x,'error return from ns01a because f(x) ',
     1'failed to decrease using a new jacobian')
      go to 17
   15 ntest=nt
c     test whether there have been maxfun calls of calfun
   11 if (maxfun-maxc) 21,21,22
   21 print 23,maxc
   23 format (///5x,'error return from ns01a because there have been',
     1i5,' calls of calfun')
      if (fsq-fmin) 3,17,17
c     provide printing if requested
   22 if (iprint) 24,998,25
   25 print 26,maxc
   26 format (/5x,'at the',i5,' th call of calfun we have')
      print 8,(i,x(i),f(i),i=1,n)
  998 print 9,maxc,fsq
   24 go to (27,28,29,87,30),is
c     store the result of the initial call of calfun
   30 fmin=fsq
      do 31 i=1,n
      w(nx+i)=x(i)
      w(nf+i)=f(i)
   31 continue
c     calculate a new jacobian approximation
   32 ic=0
      is=3
   33 ic=ic+1
      x(ic)=x(ic)+dstep
      go to 1
   29 k=ic
      do 34 i=1,n
      w(k)=(f(i)-w(nf+i))/dstep
      k=k+n
   34 continue
      x(ic)=w(nx+ic)
      if (ic-n) 33,35,35
c     calculate the inverse of the jacobian and set the direction matri
   35 k=0
      do 36 i=1,n
      do 37 j=1,n
      k=k+1
      ajinv(i,j)=w(k)
      w(nd+k)=0.0d0
   37 continue
      w(ndc+k+i)=1.0d0
      w(ndc+i)=1.0d0+dfloat(n-i)
   36 continue
c
c invert matrix with linpack routines:
      lda = 20
      info = 0
      call dgefa(ajinv,lda,n,ipvt,info)
      if (info.ne.0) goto 44
      ijob = 1
      call dgedi(ajinv,lda,n,ipvt,det,lawork,ijob)
c
c     start iteration by predicting the descent and newton minima
   38 ds=0.0d0
      dn=0.0d0
      sp=0.0d0
      do 39 i=1,n
      x(i)=0.0d0
      f(i)=0.0d0
      k=i
      do 40 j=1,n
      x(i)=x(i)-w(k)*w(nf+j)
      f(i)=f(i)-ajinv(i,j)*w(nf+j)
      k=k+n
   40 continue
      ds=ds+x(i)*x(i)
      dn=dn+f(i)*f(i)
      sp=sp+x(i)*f(i)
   39 continue
c     test whether a nearby stationary point is predicted
      if (fmin*fmin-dmm*ds) 41,41,42
c     if so then return or revise jacobian
   42 go to (43,43,44),is
   44 print 45
   45 format (///5x,'error return from ns01a because a nearby ',
     1'stationary point of f(x) is predicted')
      go to 17
   43 ntest=0
      do 46 i=1,n
      x(i)=w(nx+i)
   46 continue
      go to 32
c     test whether to apply the full newton correction
   41 is=2
      if (dn-dd) 47,47,48
   47 dd=dmax1(dn,dss)
      ds=0.25d0*dn
      tinc=1.0d0
      if (dn-dss) 49,58,58
   49 is=4
      go to 80
c     calculate the length of the steepest descent step
   48 k=0
      dmult=0.0d0
      do 51 i=1,n
      dw=0.0d0
      do 52 j=1,n
      k=k+1
      dw=dw+w(k)*x(j)
   52 continue
      dmult=dmult+dw*dw
   51 continue
      dmult=ds/dmult
      ds=ds*dmult*dmult
c     test whether to use the steepest descent direction
      if (ds-dd) 53,54,54
c     test whether the initial value of dd has been set
   54 if (dd) 55,55,56
   55 dd=dmax1(dss,dmin1(dm,ds))
      ds=ds/(dmult*dmult)
      go to 41
c     set the multiplier of the steepest descent direction
   56 anmult=0.0d0
      dmult=dmult*dsqrt(dd/ds)
      go to 98
c     interpolate between the steepest descent and the newton direction
   53 sp=sp*dmult
      anmult=(dd-ds)/((sp-ds)+dsqrt((sp-dd)**2+(dn-dd)*(dd-ds)))
      dmult=dmult*(1.0d0-anmult)
c     calculate the change in x and its angle with the first direction
   98 dn=0.0d0
      sp=0.0d0
      do 57 i=1,n
      f(i)=dmult*x(i)+anmult*f(i)
      dn=dn+f(i)*f(i)
      sp=sp+f(i)*w(nd+i)
   57 continue
      ds=0.25d0*dn
c     test whether an extra step is needed for independence
      if (w(ndc+1)-dtest) 58,58,59
   59 if (sp*sp-ds) 60,58,58
c     take the extra step and update the direction matrix
   50 is=2
   60 do 61 i=1,n
      x(i)=w(nx+i)+dstep*w(nd+i)
      w(ndc+i)=w(ndc+i+1)+1.0d0
   61 continue
      w(nd)=1.0d0
      do 62 i=1,n
      k=nd+i
      sp=w(k)
      do 63 j=2,n
      w(k)=w(k+n)
      k=k+n
   63 continue
      w(k)=sp
   62 continue
      go to 1
c     express the new direction in terms of those of the direction
c     matrix, and update the counts in w(ndc+1) etc.
   58 sp=0.0d0
      k=nd
      do 64 i=1,n
      x(i)=dw
      dw=0.0d0
      do 65 j=1,n
      k=k+1
      dw=dw+f(j)*w(k)
   65 continue
      go to (68,66),is
   66 w(ndc+i)=w(ndc+i)+1.0d0
      sp=sp+dw*dw
      if (sp-ds) 64,64,67
   67 is=1
      kk=i
      x(1)=dw
      go to 69
   68 x(i)=dw
   69 w(ndc+i)=w(ndc+i+1)+1.0d0
   64 continue
      w(nd)=1.0d0
c     reorder the directions so that kk is first
      if (kk-1) 70,70,71
   71 ks=ndc+kk*n
      do 72 i=1,n
      k=ks+i
      sp=w(k)
      do 73 j=2,kk
      w(k)=w(k-n)
      k=k-n
   73 continue
      w(k)=sp
   72 continue
c     generate the new orthogonal direction matrix
   70 do 74 i=1,n
      w(nw+i)=0.0d0
   74 continue
      sp=x(1)*x(1)
      k=nd
      do 75 i=2,n
      ds=dsqrt(sp*(sp+x(i)*x(i)))
      dw=sp/ds
      ds=x(i)/ds
      sp=sp+x(i)*x(i)
      do 76 j=1,n
      k=k+1
      w(nw+j)=w(nw+j)+x(i-1)*w(k)
      w(k)=dw*w(k+n)-ds*w(nw+j)
   76 continue
   75 continue
      sp=1.0d0/dsqrt(dn)
      do 77 i=1,n
      k=k+1
      w(k)=sp*f(i)
   77 continue
c     calculate the next vector x, and predict the right hand sides
   80 fnp=0.0d0
      k=0
      do 78 i=1,n
      x(i)=w(nx+i)+f(i)
      w(nw+i)=w(nf+i)
      do 79 j=1,n
      k=k+1
      w(nw+i)=w(nw+i)+w(k)*f(j)
   79 continue
      fnp=fnp+w(nw+i)**2
   78 continue
c     call calfun using the new vector of variables
      go to 1
c     update the step size
   27 dmult=0.9d0*fmin+0.1d0*fnp-fsq
      if (dmult) 82,81,81
   82 dd=dmax1(dss,0.25d0*dd)
      tinc=1.0d0
      if (fsq-fmin) 83,28,28
c     try the test to decide whether to increase the step length
   81 sp=0.0d0
      ss=0.0d0
      do 84 i=1,n
      sp=sp+dabs(f(i)*(f(i)-w(nw+i)))
      ss=ss+(f(i)-w(nw+i))**2
   84 continue
      pj=1.0d0+dmult/(sp+dsqrt(sp*sp+dmult*ss))
      sp=dmin1(4.0d0,tinc,pj)
      tinc=pj/sp
      dd=dmin1(dm,sp*dd)
      go to 83
c     if f(x) improves store the new value of x
   87 if (fsq-fmin) 83,50,50
   83 fmin=fsq
      do 88 i=1,n
      sp=x(i)
      x(i)=w(nx+i)
      w(nx+i)=sp
      sp=f(i)
      f(i)=w(nf+i)
      w(nf+i)=sp
      w(nw+i)=-w(nw+i)
   88 continue
      if (is-1) 28,28,50
c     calculate the changes in f and in x
   28 do 89 i=1,n
      x(i)=x(i)-w(nx+i)
      f(i)=f(i)-w(nf+i)
   89 continue
c     update the approximations to j and to ajinv
      k=0
      do 90 i=1,n
      w(mw+i)=x(i)
      w(nw+i)=f(i)
      do 91 j=1,n
      w(mw+i)=w(mw+i)-ajinv(i,j)*f(j)
      k=k+1
      w(nw+i)=w(nw+i)-w(k)*x(j)
   91 continue
   90 continue
      sp=0.0d0
      ss=0.0d0
      do 92 i=1,n
      ds=0.0d0
      do 93 j=1,n
      ds=ds+ajinv(j,i)*x(j)
   93 continue
      sp=sp+ds*f(i)
      ss=ss+x(i)*x(i)
      f(i)=ds
   92 continue
      dmult=1.0d0
      if (dabs(sp)-0.1d0*ss) 94,95,95
   94 dmult=0.8d0
   95 pj=dmult/ss
      pa=dmult/(dmult*sp+(1.0d0-dmult)*ss)
      k=0
      do 96 i=1,n
      sp=pj*w(nw+i)
      ss=pa*w(mw+i)
      do 97 j=1,n
      k=k+1
      w(k)=w(k)+sp*x(j)
      ajinv(i,j)=ajinv(i,j)+ss*f(j)
   97 continue
   96 continue
      go to 38
      end
      subroutine dgefa(a,lda,n,ipvt,info)

c*********************************************************************72

      integer lda,n,ipvt(1),info
      double precision a(lda,1)
c
c     linpack routine.
c     dgefa factors a double precision matrix by gaussian elimination.
c     requires blas daxpy,dscal,idamax
c
      double precision t
      integer idamax,j,k,kp1,l,nm1
c
c     gaussian elimination with partial pivoting
c
      info = 0
      nm1 = n - 1
      if (nm1 .lt. 1) go to 70
      do 60 k = 1, nm1
         kp1 = k + 1
c
c        find l = pivot index
c
         l = idamax(n-k+1,a(k,k),1) + k - 1
         ipvt(k) = l
c
c        zero pivot implies this column already triangularized
c
         if (a(l,k) .eq. 0.0d0) go to 40
c
c           interchange if necessary
c
            if (l .eq. k) go to 10
               t = a(l,k)
               a(l,k) = a(k,k)
               a(k,k) = t
   10       continue
c
c           compute multipliers
c
            t = -1.0d0/a(k,k)
            call dscal(n-k,t,a(k+1,k),1)
c
c           row elimination with column indexing
c
            do 30 j = kp1, n
               t = a(l,j)
               if (l .eq. k) go to 20
                  a(l,j) = a(k,j)
                  a(k,j) = t
   20          continue
               call daxpy(n-k,t,a(k+1,k),1,a(k+1,j),1)
   30       continue
         go to 50
   40    continue
            info = k
   50    continue
   60 continue
   70 continue
      ipvt(n) = n
      if (a(n,n) .eq. 0.0d0) info = n
      return
      end
      subroutine dgedi(a,lda,n,ipvt,det,work,job)

c*********************************************************************72

      integer lda,n,ipvt(1),job
      double precision a(lda,1),det(2),work(1)
c
c     linpack routine.
c     dgedi computes the determinant and inverse of a matrix
c     using the factors computed by dgeco or dgefa.
c
c     requires blas daxpy,dscal,dswap
c     fortran dabs,mod
c
      double precision t
      double precision ten
      integer i,j,k,kb,kp1,l,nm1
c
c     compute determinant
c
      if (job/10 .eq. 0) go to 70
         det(1) = 1.0d0
         det(2) = 0.0d0
         ten = 10.0d0
         do 50 i = 1, n
            if (ipvt(i) .ne. i) det(1) = -det(1)
            det(1) = a(i,i)*det(1)
c        ...exit
            if (det(1) .eq. 0.0d0) go to 60
   10       if (dabs(det(1)) .ge. 1.0d0) go to 20
               det(1) = ten*det(1)
               det(2) = det(2) - 1.0d0
            go to 10
   20       continue
   30       if (dabs(det(1)) .lt. ten) go to 40
               det(1) = det(1)/ten
               det(2) = det(2) + 1.0d0
            go to 30
   40       continue
   50    continue
   60    continue
   70 continue
c
c     compute inverse(u)
c
      if (mod(job,10) .eq. 0) go to 150
         do 100 k = 1, n
            a(k,k) = 1.0d0/a(k,k)
            t = -a(k,k)
            call dscal(k-1,t,a(1,k),1)
            kp1 = k + 1
            if (n .lt. kp1) go to 90
            do 80 j = kp1, n
               t = a(k,j)
               a(k,j) = 0.0d0
               call daxpy(k,t,a(1,k),1,a(1,j),1)
   80       continue
   90       continue
  100    continue
c
c        form inverse(u)*inverse(l)
c
         nm1 = n - 1
         if (nm1 .lt. 1) go to 140
         do 130 kb = 1, nm1
            k = n - kb
            kp1 = k + 1
            do 110 i = kp1, n
               work(i) = a(i,k)
               a(i,k) = 0.0d0
  110       continue
            do 120 j = kp1, n
               t = work(j)
               call daxpy(n,t,a(1,j),1,a(1,k),1)
  120       continue
            l = ipvt(k)
            if (l .ne. k) call dswap(n,a(1,k),1,a(1,l),1)
  130    continue
  140    continue
  150 continue
      return
      end
      subroutine daxpy(n,da,dx,incx,dy,incy)

c*********************************************************************72
c
c     constant times a vector plus a vector.
c     jack dongarra, linpack, 3/11/78.
c
      double precision dx(1),dy(1),da
      integer i,incx,incy,ixiy,m,mp1,n
      if(n.le.0)return
      if (da .eq. 0.0d0) return
      if(incx.eq.1.and.incy.eq.1)go to 20
      ix = 1
      iy = 1
      if(incx.lt.0)ix = (-n+1)*incx + 1
      if(incy.lt.0)iy = (-n+1)*incy + 1
      do 10 i = 1,n
        dy(iy) = dy(iy) + da*dx(ix)
        ix = ix + incx
        iy = iy + incy
   10 continue
      return
   20 m = mod(n,4)
      if( m .eq. 0 ) go to 40
      do 30 i = 1,m
        dy(i) = dy(i) + da*dx(i)
   30 continue
      if( n .lt. 4 ) return
   40 mp1 = m + 1
      do 50 i = mp1,n,4
        dy(i) = dy(i) + da*dx(i)
        dy(i + 1) = dy(i + 1) + da*dx(i + 1)
        dy(i + 2) = dy(i + 2) + da*dx(i + 2)
        dy(i + 3) = dy(i + 3) + da*dx(i + 3)
   50 continue
      return
      end
      subroutine  dscal(n,da,dx,incx)

c*********************************************************************72
c
c     scales a vector by a constant.
c     jack dongarra, linpack, 3/11/78.
c
      double precision da,dx(1)
      integer i,incx,m,mp1,n,nincx
      if(n.le.0)return
      if(incx.eq.1)go to 20
      nincx = n*incx
      do 10 i = 1,nincx,incx
        dx(i) = da*dx(i)
   10 continue
      return
   20 m = mod(n,5)
      if( m .eq. 0 ) go to 40
      do 30 i = 1,m
        dx(i) = da*dx(i)
   30 continue
      if( n .lt. 5 ) return
   40 mp1 = m + 1
      do 50 i = mp1,n,5
        dx(i) = da*dx(i)
        dx(i + 1) = da*dx(i + 1)
        dx(i + 2) = da*dx(i + 2)
        dx(i + 3) = da*dx(i + 3)
        dx(i + 4) = da*dx(i + 4)
   50 continue
      return
      end
      integer function idamax(n,dx,incx)

c*********************************************************************72
c
c     finds the index of element having max. absolute value.
c     jack dongarra, linpack, 3/11/78.
c
      double precision dx(1),dmax
      integer i,incx,ix,n
      idamax = 0
      if( n .lt. 1 ) return
      idamax = 1
      if(n.eq.1)return
      if(incx.eq.1)go to 20
      ix = 1
      dmax = dabs(dx(1))
      ix = ix + incx
      do 10 i = 2,n
         if(dabs(dx(ix)).le.dmax) go to 5
         idamax = i
         dmax = dabs(dx(ix))
    5    ix = ix + incx
   10 continue
      return
   20 dmax = dabs(dx(1))
      do 30 i = 2,n
         if(dabs(dx(i)).le.dmax) go to 30
         idamax = i
         dmax = dabs(dx(i))
   30 continue
      return
      end
      subroutine  dswap (n,dx,incx,dy,incy)

c*********************************************************************72
c
c     interchanges two vectors.
c     jack dongarra, linpack, 3/11/78.
c
      double precision dx(1),dy(1),dtemp
      integer i,incx,incy,ix,iy,m,mp1,n
      if(n.le.0)return
      if(incx.eq.1.and.incy.eq.1)go to 20
      ix = 1
      iy = 1
      if(incx.lt.0)ix = (-n+1)*incx + 1
      if(incy.lt.0)iy = (-n+1)*incy + 1
      do 10 i = 1,n
        dtemp = dx(ix)
        dx(ix) = dy(iy)
        dy(iy) = dtemp
        ix = ix + incx
        iy = iy + incy
   10 continue
      return
   20 m = mod(n,3)
      if( m .eq. 0 ) go to 40
      do 30 i = 1,m
        dtemp = dx(i)
        dx(i) = dy(i)
        dy(i) = dtemp
   30 continue
      if( n .lt. 3 ) return
   40 mp1 = m + 1
      do 50 i = mp1,n,3
        dtemp = dx(i)
        dx(i) = dy(i)
        dy(i) = dtemp
        dtemp = dx(i + 1)
        dx(i + 1) = dy(i + 1)
        dy(i + 1) = dtemp
        dtemp = dx(i + 2)
        dx(i + 2) = dy(i + 2)
        dy(i + 2) = dtemp
   50 continue
      return
      end
