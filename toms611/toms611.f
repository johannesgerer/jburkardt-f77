      integer function imdcon(k)
c
      integer k
c
c  ***  return integer machine-dependent constants  ***
c
c     ***  k = 1 means return standard output unit number.   ***
c     ***  k = 2 means return alternate output unit number.  ***
c     ***  k = 3 means return  input unit number.            ***
c          (note -- k = 2, 3 are used only by test programs.)
c
c  +++  port version follows...
c
      external i1mach
      integer i1mach
      integer mdperm(3)
      data mdperm(1)/2/, mdperm(2)/4/, mdperm(3)/1/
      imdcon = i1mach(mdperm(k))
c  +++  end of port version  +++
c
c  +++  non-port version follows...
c     integer mdcon(3)
c     data mdcon(1)/6/, mdcon(2)/8/, mdcon(3)/5/
c     imdcon = mdcon(k)
c  +++  end of non-port version  +++
c
 999  return
c  ***  last card of imdcon follows  ***
      end
      double precision function rmdcon(k)
c
c  ***  return machine dependent constants used by nl2sol  ***
c
c +++  comments below contain data statements for various machines.  +++
c +++  to convert to another machine, place a c in column 1 of the   +++
c +++  data statement line(s) that correspond to the current machine +++
c +++  and remove the c from column 1 of the data statement line(s)  +++
c +++  that correspond to the new machine.                           +++
c
      integer k
c
c  ***  the constant returned depends on k...
c
c  ***        k = 1... smallest pos. eta such that -eta exists.
c  ***        k = 2... square root of eta.
c  ***        k = 3... unit roundoff = smallest pos. no. machep such
c  ***                 that 1 + machep .gt. 1 .and. 1 - machep .lt. 1.
c  ***        k = 4... square root of machep.
c  ***        k = 5... square root of big (see k = 6).
c  ***        k = 6... largest machine no. big such that -big exists.
c
      double precision big, eta, machep
      integer bigi(4), etai(4), machei(4)
c/+
      double precision dsqrt
c/
      equivalence (big,bigi(1)), (eta,etai(1)), (machep,machei(1))
c
c  +++  ibm 360, ibm 370, or xerox  +++
c
c     data big/z7fffffffffffffff/, eta/z0010000000000000/,
c    1     machep/z3410000000000000/
c
c  +++  data general  +++
c
c     data big/0.7237005577d+76/, eta/0.5397605347d-78/,
c    1     machep/2.22044605d-16/
c
c  +++  dec 11  +++
c
c     data big/1.7d+38/, eta/2.938735878d-39/, machep/2.775557562d-17/
c
c  +++  hp3000  +++
c
c     data big/1.157920892d+77/, eta/8.636168556d-78/,
c    1     machep/5.551115124d-17/
c
c  +++  honeywell  +++
c
c     data big/1.69d+38/, eta/5.9d-39/, machep/2.1680435d-19/
c
c  +++  dec10  +++
c
c     data big/"377777100000000000000000/,
c    1     eta/"002400400000000000000000/,
c    2     machep/"104400000000000000000000/
c
c  +++  burroughs  +++
c
c     data big/o0777777777777777,o7777777777777777/,
c    1     eta/o1771000000000000,o7770000000000000/,
c    2     machep/o1451000000000000,o0000000000000000/
c
c  +++  control data  +++
c
c     data big/37767777777777777777b,37167777777777777777b/,
c    1     eta/00014000000000000000b,00000000000000000000b/,
c    2     machep/15614000000000000000b,15010000000000000000b/
c
c  +++  prime  +++
c
c     data big/1.0d+9786/, eta/1.0d-9860/, machep/1.4210855d-14/
c
c  +++  univac  +++
c
c     data big/8.988d+307/, eta/1.2d-308/, machep/1.734723476d-18/
c
c  +++  vax  +++
c
c     data big/1.7d+38/, eta/2.939d-39/, machep/1.3877788d-17/
c
c  +++  cray 1  +++
c
c     data bigi(1)/577767777777777777777b/,
c    1     bigi(2)/000007777777777777776b/,
c    2     etai(1)/200004000000000000000b/,
c    3     etai(2)/000000000000000000000b/,
c    4     machei(1)/377224000000000000000b/,
c    5     machei(2)/000000000000000000000b/
c
c  +++  port library -- requires more than just a data statement... +++
c
      external d1mach
      double precision d1mach, zero
      data big/0.d+0/, eta/0.d+0/, machep/0.d+0/, zero/0.d+0/
      if (big .gt. zero) go to 1
         big = d1mach(2)
         eta = d1mach(1)
         machep = d1mach(4)
 1    continue
c
c  +++ end of port +++
c
c-------------------------------  body  --------------------------------
c
      go to (10, 20, 30, 40, 50, 60), k
c
 10   rmdcon = eta
      go to 999
c
 20   rmdcon = dsqrt(256.d+0*eta)/16.d+0
      go to 999
c
 30   rmdcon = machep
      go to 999
c
 40   rmdcon = dsqrt(machep)
      go to 999
c
 50   rmdcon = dsqrt(big/256.d+0)*16.d+0
      go to 999
c
 60   rmdcon = big
c
 999  return
c  ***  last card of rmdcon follows  ***
      end
      subroutine sumsl(n, d, x, calcf, calcg, iv, liv, lv, v,
     1                  uiparm, urparm, ufparm)
c
c  ***  minimize general unconstrained objective function using   ***
c  ***  analytic gradient and hessian approx. from secant update  ***
c
      integer n, liv, lv
      integer iv(liv), uiparm(1)
      double precision d(n), x(n), v(lv), urparm(1)
c     dimension v(71 + n*(n+15)/2), uiparm(*), urparm(*)
      external calcf, calcg, ufparm
c
c  ***  purpose  ***
c
c        this routine interacts with subroutine  sumit  in an attempt
c     to find an n-vector  x*  that minimizes the (unconstrained)
c     objective function computed by  calcf.  (often the  x*  found is
c     a local minimizer rather than a global one.)
c
c--------------------------  parameter usage  --------------------------
c
c n........ (input) the number of variables on which  f  depends, i.e.,
c                  the number of components in  x.
c d........ (input/output) a scale vector such that  d(i)*x(i),
c                  i = 1,2,...,n,  are all in comparable units.
c                  d can strongly affect the behavior of sumsl.
c                  finding the best choice of d is generally a trial-
c                  and-error process.  choosing d so that d(i)*x(i)
c                  has about the same value for all i often works well.
c                  the defaults provided by subroutine deflt (see iv
c                  below) require the caller to supply d.
c x........ (input/output) before (initially) calling sumsl, the call-
c                  er should set  x  to an initial guess at  x*.  when
c                  sumsl returns,  x  contains the best point so far
c                  found, i.e., the one that gives the least value so
c                  far seen for  f(x).
c calcf.... (input) a subroutine that, given x, computes f(x).  calcf
c                  must be declared external in the calling program.
c                  it is invoked by
c                       call calcf(n, x, nf, f, uiparm, urparm, ufparm)
c                  when calcf is called, nf is the invocation
c                  count for calcf.  nf is included for possible use
c                  with calcg.  if x is out of bounds (e.g., if it
c                  would cause overflow in computing f(x)), then calcf
c                  should set nf to 0.  this will cause a shorter step
c                  to be attempted.  (if x is in bounds, then calcf
c                  should not change nf.)  the other parameters are as
c                  described above and below.  calcf should not change
c                  n, p, or x.
c calcg.... (input) a subroutine that, given x, computes g(x), the gra-
c                  dient of f at x.  calcg must be declared external in
c                  the calling program.  it is invoked by
c                       call calcg(n, x, nf, g, uiparm, urparm, ufaprm)
c                  when calcg is called, nf is the invocation
c                  count for calcf at the time f(x) was evaluated.  the
c                  x passed to calcg is usually the one passed to calcf
c                  on either its most recent invocation or the one
c                  prior to it.  if calcf saves intermediate results
c                  for use by calcg, then it is possible to tell from
c                  nf whether they are valid for the current x (or
c                  which copy is valid if two copies are kept).  if g
c                  cannot be computed at x, then calcg should set nf to
c                  0.  in this case, sumsl will return with iv(1) = 65.
c                  (if g can be computed at x, then calcg should not
c                  changed nf.)  the other parameters to calcg are as
c                  described above and below.  calcg should not change
c                  n or x.
c iv....... (input/output) an integer value array of length liv (see
c                  below) that helps control the sumsl algorithm and
c                  that is used to store various intermediate quanti-
c                  ties.  of particular interest are the initialization/
c                  return code iv(1) and the entries in iv that control
c                  printing and limit the number of iterations and func-
c                  tion evaluations.  see the section on iv input
c                  values below.
c liv...... (input) length of iv array.  must be at least 60.  if liv
c                  is too small, then sumsl returns with iv(1) = 15.
c                  when sumsl returns, the smallest allowed value of
c                  liv is stored in iv(lastiv) -- see the section on
c                  iv output values below.  (this is intended for use
c                  with extensions of sumsl that handle constraints.)
c lv....... (input) length of v array.  must be at least 71+n*(n+15)/2.
c                  (at least 77+n*(n+17)/2 for smsno, at least
c                  78+n*(n+12) for humsl).  if lv is too small, then
c                  sumsl returns with iv(1) = 16.  when sumsl returns,
c                  the smallest allowed value of lv is stored in
c                  iv(lastv) -- see the section on iv output values
c                  below.
c v........ (input/output) a floating-point value array of length lv
c                  (see below) that helps control the sumsl algorithm
c                  and that is used to store various intermediate
c                  quantities.  of particular interest are the entries
c                  in v that limit the length of the first step
c                  attempted (lmax0) and specify convergence tolerances
c                  (afctol, lmaxs, rfctol, sctol, xctol, xftol).
c uiparm... (input) user integer parameter array passed without change
c                  to calcf and calcg.
c urparm... (input) user floating-point parameter array passed without
c                  change to calcf and calcg.
c ufparm... (input) user external subroutine or function passed without
c                  change to calcf and calcg.
c
c  ***  iv input values (from subroutine deflt)  ***
c
c iv(1)...  on input, iv(1) should have a value between 0 and 14......
c             0 and 12 mean this is a fresh start.  0 means that
c                  deflt(2, iv, liv, lv, v)
c             is to be called to provide all default values to iv and
c             v.  12 (the value that deflt assigns to iv(1)) means the
c             caller has already called deflt and has possibly changed
c             some iv and/or v entries to non-default values.
c             13 means deflt has been called and that sumsl (and
c             sumit) should only do their storage allocation.  that is,
c             they should set the output components of iv that tell
c             where various subarrays arrays of v begin, such as iv(g)
c             (and, for humsl and humit only, iv(dtol)), and return.
c             14 means that a storage has been allocated (by a call
c             with iv(1) = 13) and that the algorithm should be
c             started.  when called with iv(1) = 13, sumsl returns
c             iv(1) = 14 unless liv or lv is too small (or n is not
c             positive).  default = 12.
c iv(inith).... iv(25) tells whether the hessian approximation h should
c             be initialized.  1 (the default) means sumit should
c             initialize h to the diagonal matrix whose i-th diagonal
c             element is d(i)**2.  0 means the caller has supplied a
c             cholesky factor  l  of the initial hessian approximation
c             h = l*(l**t)  in v, starting at v(iv(lmat)) = v(iv(42))
c             (and stored compactly by rows).  note that iv(lmat) may
c             be initialized by calling sumsl with iv(1) = 13 (see
c             the iv(1) discussion above).  default = 1.
c iv(mxfcal)... iv(17) gives the maximum number of function evaluations
c             (calls on calcf) allowed.  if this number does not suf-
c             fice, then sumsl returns with iv(1) = 9.  default = 200.
c iv(mxiter)... iv(18) gives the maximum number of iterations allowed.
c             it also indirectly limits the number of gradient evalua-
c             tions (calls on calcg) to iv(mxiter) + 1.  if iv(mxiter)
c             iterations do not suffice, then sumsl returns with
c             iv(1) = 10.  default = 150.
c iv(outlev)... iv(19) controls the number and length of iteration sum-
c             mary lines printed (by itsum).  iv(outlev) = 0 means do
c             not print any summary lines.  otherwise, print a summary
c             line after each abs(iv(outlev)) iterations.  if iv(outlev)
c             is positive, then summary lines of length 78 (plus carri-
c             age control) are printed, including the following...  the
c             iteration and function evaluation counts, f = the current
c             function value, relative difference in function values
c             achieved by the latest step (i.e., reldf = (f0-v(f))/f01,
c             where f01 is the maximum of abs(v(f)) and abs(v(f0)) and
c             v(f0) is the function value from the previous itera-
c             tion), the relative function reduction predicted for the
c             step just taken (i.e., preldf = v(preduc) / f01, where
c             v(preduc) is described below), the scaled relative change
c             in x (see v(reldx) below), the step parameter for the
c             step just taken (stppar = 0 means a full newton step,
c             between 0 and 1 means a relaxed newton step, between 1
c             and 2 means a double dogleg step, greater than 2 means
c             a scaled down cauchy step -- see subroutine dbldog), the
c             2-norm of the scale vector d times the step just taken
c             (see v(dstnrm) below), and npreldf, i.e.,
c             v(nreduc)/f01, where v(nreduc) is described below -- if
c             npreldf is positive, then it is the relative function
c             reduction predicted for a newton step (one with
c             stppar = 0).  if npreldf is negative, then it is the
c             negative of the relative function reduction predicted
c             for a step computed with step bound v(lmaxs) for use in
c             testing for singular convergence.
c                  if iv(outlev) is negative, then lines of length 50
c             are printed, including only the first 6 items listed
c             above (through reldx).
c             default = 1.
c iv(parprt)... iv(20) = 1 means print any nondefault v values on a
c             fresh start or any changed v values on a restart.
c             iv(parprt) = 0 means skip this printing.  default = 1.
c iv(prunit)... iv(21) is the output unit number on which all printing
c             is done.  iv(prunit) = 0 means suppress all printing.
c             default = standard output unit (unit 6 on most systems).
c iv(solprt)... iv(22) = 1 means print out the value of x returned (as
c             well as the gradient and the scale vector d).
c             iv(solprt) = 0 means skip this printing.  default = 1.
c iv(statpr)... iv(23) = 1 means print summary statistics upon return-
c             ing.  these consist of the function value, the scaled
c             relative change in x caused by the most recent step (see
c             v(reldx) below), the number of function and gradient
c             evaluations (calls on calcf and calcg), and the relative
c             function reductions predicted for the last step taken and
c             for a newton step (or perhaps a step bounded by v(lmaxs)
c             -- see the descriptions of preldf and npreldf under
c             iv(outlev) above).
c             iv(statpr) = 0 means skip this printing.
c             iv(statpr) = -1 means skip this printing as well as that
c             of the one-line termination reason message.  default = 1.
c iv(x0prt).... iv(24) = 1 means print the initial x and scale vector d
c             (on a fresh start only).  iv(x0prt) = 0 means skip this
c             printing.  default = 1.
c
c  ***  (selected) iv output values  ***
c
c iv(1)........ on output, iv(1) is a return code....
c             3 = x-convergence.  the scaled relative difference (see
c                  v(reldx)) between the current parameter vector x and
c                  a locally optimal parameter vector is very likely at
c                  most v(xctol).
c             4 = relative function convergence.  the relative differ-
c                  ence between the current function value and its lo-
c                  cally optimal value is very likely at most v(rfctol).
c             5 = both x- and relative function convergence (i.e., the
c                  conditions for iv(1) = 3 and iv(1) = 4 both hold).
c             6 = absolute function convergence.  the current function
c                  value is at most v(afctol) in absolute value.
c             7 = singular convergence.  the hessian near the current
c                  iterate appears to be singular or nearly so, and a
c                  step of length at most v(lmaxs) is unlikely to yield
c                  a relative function decrease of more than v(sctol).
c             8 = false convergence.  the iterates appear to be converg-
c                  ing to a noncritical point.  this may mean that the
c                  convergence tolerances (v(afctol), v(rfctol),
c                  v(xctol)) are too small for the accuracy to which
c                  the function and gradient are being computed, that
c                  there is an error in computing the gradient, or that
c                  the function or gradient is discontinuous near x.
c             9 = function evaluation limit reached without other con-
c                  vergence (see iv(mxfcal)).
c            10 = iteration limit reached without other convergence
c                  (see iv(mxiter)).
c            11 = stopx returned .true. (external interrupt).  see the
c                  usage notes below.
c            14 = storage has been allocated (after a call with
c                  iv(1) = 13).
c            17 = restart attempted with n changed.
c            18 = d has a negative component and iv(dtype) .le. 0.
c            19...43 = v(iv(1)) is out of range.
c            63 = f(x) cannot be computed at the initial x.
c            64 = bad parameters passed to assess (which should not
c                  occur).
c            65 = the gradient could not be computed at x (see calcg
c                  above).
c            67 = bad first parameter to deflt.
c            80 = iv(1) was out of range.
c            81 = n is not positive.
c iv(g)........ iv(28) is the starting subscript in v of the current
c             gradient vector (the one corresponding to x).
c iv(lastiv)... iv(44) is the least acceptable value of liv.  (it is
c             only set if liv is at least 44.)
c iv(lastv).... iv(45) is the least acceptable value of lv.  (it is
c             only set if liv is large enough, at least iv(lastiv).)
c iv(nfcall)... iv(6) is the number of calls so far made on calcf (i.e.,
c             function evaluations).
c iv(ngcall)... iv(30) is the number of gradient evaluations (calls on
c             calcg).
c iv(niter).... iv(31) is the number of iterations performed.
c
c  ***  (selected) v input values (from subroutine deflt)  ***
c
c v(bias)..... v(43) is the bias parameter used in subroutine dbldog --
c             see that subroutine for details.  default = 0.8.
c v(afctol)... v(31) is the absolute function convergence tolerance.
c             if sumsl finds a point where the function value is less
c             than v(afctol) in absolute value, and if sumsl does not
c             return with iv(1) = 3, 4, or 5, then it returns with
c             iv(1) = 6.  this test can be turned off by setting
c             v(afctol) to zero.  default = max(10**-20, machep**2),
c             where machep is the unit roundoff.
c v(dinit).... v(38), if nonnegative, is the value to which the scale
c             vector d is initialized.  default = -1.
c v(lmax0).... v(35) gives the maximum 2-norm allowed for d times the
c             very first step that sumsl attempts.  this parameter can
c             markedly affect the performance of sumsl.
c v(lmaxs).... v(36) is used in testing for singular convergence -- if
c             the function reduction predicted for a step of length
c             bounded by v(lmaxs) is at most v(sctol) * abs(f0), where
c             f0  is the function value at the start of the current
c             iteration, and if sumsl does not return with iv(1) = 3,
c             4, 5, or 6, then it returns with iv(1) = 7.  default = 1.
c v(rfctol)... v(32) is the relative function convergence tolerance.
c             if the current model predicts a maximum possible function
c             reduction (see v(nreduc)) of at most v(rfctol)*abs(f0)
c             at the start of the current iteration, where  f0  is the
c             then current function value, and if the last step attempt-
c             ed achieved no more than twice the predicted function
c             decrease, then sumsl returns with iv(1) = 4 (or 5).
c             default = max(10**-10, machep**(2/3)), where machep is
c             the unit roundoff.
c v(sctol).... v(37) is the singular convergence tolerance -- see the
c             description of v(lmaxs) above.
c v(tuner1)... v(26) helps decide when to check for false convergence.
c             this is done if the actual function decrease from the
c             current step is no more than v(tuner1) times its predict-
c             ed value.  default = 0.1.
c v(xctol).... v(33) is the x-convergence tolerance.  if a newton step
c             (see v(nreduc)) is tried that has v(reldx) .le. v(xctol)
c             and if this step yields at most twice the predicted func-
c             tion decrease, then sumsl returns with iv(1) = 3 (or 5).
c             (see the description of v(reldx) below.)
c             default = machep**0.5, where machep is the unit roundoff.
c v(xftol).... v(34) is the false convergence tolerance.  if a step is
c             tried that gives no more than v(tuner1) times the predict-
c             ed function decrease and that has v(reldx) .le. v(xftol),
c             and if sumsl does not return with iv(1) = 3, 4, 5, 6, or
c             7, then it returns with iv(1) = 8.  (see the description
c             of v(reldx) below.)  default = 100*machep, where
c             machep is the unit roundoff.
c v(*)........ deflt supplies to v a number of tuning constants, with
c             which it should ordinarily be unnecessary to tinker.  see
c             section 17 of version 2.2 of the nl2sol usage summary
c             (i.e., the appendix to ref. 1) for details on v(i),
c             i = decfac, incfac, phmnfc, phmxfc, rdfcmn, rdfcmx,
c             tuner2, tuner3, tuner4, tuner5.
c
c  ***  (selected) v output values  ***
c
c v(dgnorm)... v(1) is the 2-norm of (diag(d)**-1)*g, where g is the
c             most recently computed gradient.
c v(dstnrm)... v(2) is the 2-norm of diag(d)*step, where step is the
c             current step.
c v(f)........ v(10) is the current function value.
c v(f0)....... v(13) is the function value at the start of the current
c             iteration.
c v(nreduc)... v(6), if positive, is the maximum function reduction
c             possible according to the current model, i.e., the func-
c             tion reduction predicted for a newton step (i.e.,
c             step = -h**-1 * g,  where  g  is the current gradient and
c             h is the current hessian approximation).
c                  if v(nreduc) is negative, then it is the negative of
c             the function reduction predicted for a step computed with
c             a step bound of v(lmaxs) for use in testing for singular
c             convergence.
c v(preduc)... v(7) is the function reduction predicted (by the current
c             quadratic model) for the current step.  this (divided by
c             v(f0)) is used in testing for relative function
c             convergence.
c v(reldx).... v(17) is the scaled relative change in x caused by the
c             current step, computed as
c                  max(abs(d(i)*(x(i)-x0(i)), 1 .le. i .le. p) /
c                     max(d(i)*(abs(x(i))+abs(x0(i))), 1 .le. i .le. p),
c             where x = x0 + step.
c
c-------------------------------  notes  -------------------------------
c
c  ***  algorithm notes  ***
c
c        this routine uses a hessian approximation computed from the
c     bfgs update (see ref 3).  only a cholesky factor of the hessian
c     approximation is stored, and this is updated using ideas from
c     ref. 4.  steps are computed by the double dogleg scheme described
c     in ref. 2.  the steps are assessed as in ref. 1.
c
c  ***  usage notes  ***
c
c        after a return with iv(1) .le. 11, it is possible to restart,
c     i.e., to change some of the iv and v input values described above
c     and continue the algorithm from the point where it was interrupt-
c     ed.  iv(1) should not be changed, nor should any entries of iv
c     and v other than the input values (those supplied by deflt).
c        those who do not wish to write a calcg which computes the
c     gradient analytically should call smsno rather than sumsl.
c     smsno uses finite differences to compute an approximate gradient.
c        those who would prefer to provide f and g (the function and
c     gradient) by reverse communication rather than by writing subrou-
c     tines calcf and calcg may call on sumit directly.  see the com-
c     ments at the beginning of sumit.
c        those who use sumsl interactively may wish to supply their
c     own stopx function, which should return .true. if the break key
c     has been pressed since stopx was last invoked.  this makes it
c     possible to externally interrupt sumsl (which will return with
c     iv(1) = 11 if stopx returns .true.).
c        storage for g is allocated at the end of v.  thus the caller
c     may make v longer than specified above and may allow calcg to use
c     elements of g beyond the first n as scratch storage.
c
c  ***  portability notes  ***
c
c        the sumsl distribution tape contains both single- and double-
c     precision versions of the sumsl source code, so it should be un-
c     necessary to change precisions.
c        only the functions imdcon and rmdcon contain machine-dependent
c     constants.  to change from one machine to another, it should
c     suffice to change the (few) relevant lines in these functions.
c        intrinsic functions are explicitly declared.  on certain com-
c     puters (e.g. univac), it may be necessary to comment out these
c     declarations.  so that this may be done automatically by a simple
c     program, such declarations are preceded by a comment having c/+
c     in columns 1-3 and blanks in columns 4-72 and are followed by
c     a comment having c/ in columns 1 and 2 and blanks in columns 3-72.
c        the sumsl source code is expressed in 1966 ansi standard
c     fortran.  it may be converted to fortran 77 by commenting out all
c     lines that fall between a line having c/6 in columns 1-3 and a
c     line having c/7 in columns 1-3 and by removing (i.e., replacing
c     by a blank) the c in column 1 of the lines that follow the c/7
c     line and precede a line having c/ in columns 1-2 and blanks in
c     columns 3-72.  these changes convert some data statements into
c     parameter statements, convert some variables from real to
c     character*4, and make the data statements that initialize these
c     variables use character strings delimited by primes instead
c     of hollerith constants.  (such variables and data statements
c     appear only in modules itsum and parck.  parameter statements
c     appear nearly everywhere.)  these changes also add save state-
c     ments for variables given machine-dependent constants by rmdcon.
c
c  ***  references  ***
c
c 1.  dennis, j.e., gay, d.m., and welsch, r.e. (1981), algorithm 573 --
c             an adaptive nonlinear least-squares algorithm, acm trans.
c             math. software 7, pp. 369-383.
c
c 2.  dennis, j.e., and mei, h.h.w. (1979), two new unconstrained opti-
c             mization algorithms which use function and gradient
c             values, j. optim. theory applic. 28, pp. 453-482.
c
c 3.  dennis, j.e., and more, j.j. (1977), quasi-newton methods, motiva-
c             tion and theory, siam rev. 19, pp. 46-89.
c
c 4.  goldfarb, d. (1976), factorized variable metric methods for uncon-
c             strained optimization, math. comput. 30, pp. 796-811.
c
c  ***  general  ***
c
c     coded by david m. gay (winter 1980).  revised summer 1982.
c     this subroutine was written in connection with research
c     supported in part by the national science foundation under
c     grants mcs-7600324, dcr75-10143, 76-14311dss, mcs76-11989,
c     and mcs-7906671.
c.
c
c----------------------------  declarations  ---------------------------
c
      external deflt, sumit
c
c deflt... supplies default iv and v input components.
c sumit... reverse-communication routine that carries out sumsl algo-
c             rithm.
c
      integer g1, iv1, nf
      double precision f
c
c  ***  subscripts for iv   ***
c
      integer nextv, nfcall, nfgcal, g, toobig, vneed
c
c/6
      data nextv/47/, nfcall/6/, nfgcal/7/, g/28/, toobig/2/, vneed/4/
c/7
c     parameter (nextv=47, nfcall=6, nfgcal=7, g=28, toobig=2, vneed=4)
c/
c
c+++++++++++++++++++++++++++++++  body  ++++++++++++++++++++++++++++++++
c
      if (iv(1) .eq. 0) call deflt(2, iv, liv, lv, v)
      iv1 = iv(1)
      if (iv1 .eq. 12 .or. iv1 .eq. 13) iv(vneed) = iv(vneed) + n
      if (iv1 .eq. 14) go to 10
      if (iv1 .gt. 2 .and. iv1 .lt. 12) go to 10
      g1 = 1
      if (iv1 .eq. 12) iv(1) = 13
      go to 20
c
 10   g1 = iv(g)
c
 20   call sumit(d, f, v(g1), iv, liv, lv, n, v, x)
      if (iv(1) - 2) 30, 40, 50
c
 30   nf = iv(nfcall)
      call calcf(n, x, nf, f, uiparm, urparm, ufparm)
      if (nf .le. 0) iv(toobig) = 1
      go to 20
c
 40   call calcg(n, x, iv(nfgcal), v(g1), uiparm, urparm, ufparm)
      go to 20
c
 50   if (iv(1) .ne. 14) go to 999
c
c  ***  storage allocation
c
      iv(g) = iv(nextv)
      iv(nextv) = iv(g) + n
      if (iv1 .ne. 13) go to 10
c
 999  return
c  ***  last card of sumsl follows  ***
      end
      subroutine smsno(n, d, x, calcf, iv, liv, lv, v,
     1                  uiparm, urparm, ufparm)
c
c  ***  minimize general unconstrained objective function using
c  ***  finite-difference gradients and secant hessian approximations.
c
      integer n, liv, lv
      integer iv(liv), uiparm(1)
      double precision d(n), x(n), v(lv), urparm(1)
c     dimension v(77 + n*(n+17)/2), uiparm(*), urparm(*)
      external calcf, ufparm
c
c  ***  purpose  ***
c
c        this routine interacts with subroutine  snoit  in an attempt
c     to find an n-vector  x*  that minimizes the (unconstrained)
c     objective function computed by  calcf.  (often the  x*  found is
c     a local minimizer rather than a global one.)
c
c  ***  parameters  ***
c
c        the parameters for smsno are the same as those for sumsl
c     (which see), except that calcg is omitted.  instead of calling
c     calcg to obtain the gradient of the objective function at x,
c     smsno calls sgrad2, which computes an approximation to the
c     gradient by finite (forward and central) differences using the
c     method of ref. 1.  the following input component is of interest
c     in this regard (and is not described in sumsl).
c
c v(eta0)..... v(42) is an estimated bound on the relative error in the
c             objective function value computed by calcf...
c                  (true value) = (computed value) * (1 + e),
c             where abs(e) .le. v(eta0).  default = machep * 10**3,
c             where machep is the unit roundoff.
c
c        the output values iv(nfcall) and iv(ngcall) have different
c     meanings for smsno than for sumsl...
c
c iv(nfcall)... iv(6) is the number of calls so far made on calcf (i.e.,
c             function evaluations) excluding those made only for
c             computing gradients.  the input value iv(mxfcal) is a
c             limit on iv(nfcall).
c iv(ngcall)... iv(30) is the number of function evaluations made only
c             for computing gradients.  the total number of function
c             evaluations is thus  iv(nfcall) + iv(ngcall).
c
c  ***  reference  ***
c
c 1. stewart, g.w. (1967), a modification of davidon*s minimization
c        method to accept difference approximations of derivatives,
c        j. assoc. comput. mach. 14, pp. 72-83.
c.
c  ***  general  ***
c
c     coded by david m. gay (winter 1980).  revised sept. 1982.
c     this subroutine was written in connection with research
c     supported in part by the national science foundation under
c     grants mcs-7600324, dcr75-10143, 76-14311dss, mcs76-11989,
c     and mcs-7906671.
c
c
c----------------------------  declarations  ---------------------------
c
      external snoit
c
c snoit.... oversees computation of finite-difference gradient and
c         calls sumit to carry out sumsl algorithm.
c
      integer nf
      double precision fx
c
c  ***  subscripts for iv   ***
c
      integer nfcall, toobig
c
c/6
      data nfcall/6/, toobig/2/
c/7
c     parameter (nfcall=6, toobig=2)
c/
c
c+++++++++++++++++++++++++++++++  body  ++++++++++++++++++++++++++++++++
c
 10   call snoit(d, fx, iv, liv, lv, n, v, x)
      if (iv(1) .gt. 2) go to 999
c
c     ***  compute function  ***
c
      nf = iv(nfcall)
      call calcf(n, x, nf, fx, uiparm, urparm, ufparm)
      if (nf .le. 0) iv(toobig) = 1
      go to 10
c
c
 999  return
c  ***  last card of smsno follows  ***
      end
      subroutine sumit(d, fx, g, iv, liv, lv, n, v, x)
c
c  ***  carry out sumsl (unconstrained minimization) iterations, using
c  ***  double-dogleg/bfgs steps.
c
c  ***  parameter declarations  ***
c
      integer liv, lv, n
      integer iv(liv)
      double precision d(n), fx, g(n), v(lv), x(n)
c
c--------------------------  parameter usage  --------------------------
c
c d.... scale vector.
c fx... function value.
c g.... gradient vector.
c iv... integer value array.
c liv.. length of iv (at least 60).
c lv... length of v (at least 71 + n*(n+13)/2).
c n.... number of variables (components in x and g).
c v.... floating-point value array.
c x.... vector of parameters to be optimized.
c
c  ***  discussion  ***
c
c        parameters iv, n, v, and x are the same as the corresponding
c     ones to sumsl (which see), except that v can be shorter (since
c     the part of v that sumsl uses for storing g is not needed).
c     moreover, compared with sumsl, iv(1) may have the two additional
c     output values 1 and 2, which are explained below, as is the use
c     of iv(toobig) and iv(nfgcal).  the value iv(g), which is an
c     output value from sumsl (and smsno), is not referenced by
c     sumit or the subroutines it calls.
c        fx and g need not have been initialized when sumit is called
c     with iv(1) = 12, 13, or 14.
c
c iv(1) = 1 means the caller should set fx to f(x), the function value
c             at x, and call sumit again, having changed none of the
c             other parameters.  an exception occurs if f(x) cannot be
c             (e.g. if overflow would occur), which may happen because
c             of an oversized step.  in this case the caller should set
c             iv(toobig) = iv(2) to 1, which will cause sumit to ig-
c             nore fx and try a smaller step.  the parameter nf that
c             sumsl passes to calcf (for possible use by calcg) is a
c             copy of iv(nfcall) = iv(6).
c iv(1) = 2 means the caller should set g to g(x), the gradient vector
c             of f at x, and call sumit again, having changed none of
c             the other parameters except possibly the scale vector d
c             when iv(dtype) = 0.  the parameter nf that sumsl passes
c             to calcg is iv(nfgcal) = iv(7).  if g(x) cannot be
c             evaluated, then the caller may set iv(nfgcal) to 0, in
c             which case sumit will return with iv(1) = 65.
c.
c  ***  general  ***
c
c     coded by david m. gay (december 1979).  revised sept. 1982.
c     this subroutine was written in connection with research supported
c     in part by the national science foundation under grants
c     mcs-7600324 and mcs-7906671.
c
c        (see sumsl for references.)
c
c+++++++++++++++++++++++++++  declarations  ++++++++++++++++++++++++++++
c
c  ***  local variables  ***
c
      integer dg1, dummy, g01, i, k, l, lstgst, nwtst1, step1,
     1        temp1, w, x01, z
      double precision t
c
c     ***  constants  ***
c
      double precision half, negone, one, onep2, zero
c
c  ***  no intrinsic functions  ***
c
c  ***  external functions and subroutines  ***
c
      external assst, dbdog, deflt, dotprd, itsum, litvmu, livmul,
     1         ltvmul, lupdat, lvmul, parck, reldst, stopx, vaxpy,
     2         vcopy, vscopy, vvmulp, v2norm, wzbfgs
      logical stopx
      double precision dotprd, reldst, v2norm
c
c assst.... assesses candidate step.
c dbdog.... computes double-dogleg (candidate) step.
c deflt.... supplies default iv and v input components.
c dotprd... returns inner product of two vectors.
c itsum.... prints iteration summary and info on initial and final x.
c litvmu... multiplies inverse transpose of lower triangle times vector.
c livmul... multiplies inverse of lower triangle times vector.
c ltvmul... multiplies transpose of lower triangle times vector.
c lupdt.... updates cholesky factor of hessian approximation.
c lvmul.... multiplies lower triangle times vector.
c parck.... checks validity of input iv and v values.
c reldst... computes v(reldx) = relative step size.
c stopx.... returns .true. if the break key has been pressed.
c vaxpy.... computes scalar times one vector plus another.
c vcopy.... copies one vector to another.
c vscopy... sets all elements of a vector to a scalar.
c vvmulp... multiplies vector by vector raised to power (componentwise).
c v2norm... returns the 2-norm of a vector.
c wzbfgs... computes w and z for lupdat corresponding to bfgs update.
c
c  ***  subscripts for iv and v  ***
c
      integer cnvcod, dg, dgnorm, dinit, dstnrm, dst0, f, f0, fdif,
     1        gthg, gtstep, g0, incfac, inith, irc, kagqt, lmat, lmax0,
     2        lmaxs, mode, model, mxfcal, mxiter, nextv, nfcall, nfgcal,
     3        ngcall, niter, nreduc, nwtstp, preduc, radfac, radinc,
     4        radius, rad0, reldx, restor, step, stglim, stlstg, toobig,
     5        tuner4, tuner5, vneed, xirc, x0
c
c  ***  iv subscript values  ***
c
c/6
      data cnvcod/55/, dg/37/, g0/48/, inith/25/, irc/29/, kagqt/33/,
     1     mode/35/, model/5/, mxfcal/17/, mxiter/18/, nfcall/6/,
     2     nfgcal/7/, ngcall/30/, niter/31/, nwtstp/34/, radinc/8/,
     3     restor/9/, step/40/, stglim/11/, stlstg/41/, toobig/2/,
     4     vneed/4/, xirc/13/, x0/43/
c/7
c     parameter (cnvcod=55, dg=37, g0=48, inith=25, irc=29, kagqt=33,
c    1           mode=35, model=5, mxfcal=17, mxiter=18, nfcall=6,
c    2           nfgcal=7, ngcall=30, niter=31, nwtstp=34, radinc=8,
c    3           restor=9, step=40, stglim=11, stlstg=41, toobig=2,
c    4           vneed=4, xirc=13, x0=43)
c/
c
c  ***  v subscript values  ***
c
c/6
      data dgnorm/1/, dinit/38/, dstnrm/2/, dst0/3/, f/10/, f0/13/,
     1     fdif/11/, gthg/44/, gtstep/4/, incfac/23/, lmat/42/,
     2     lmax0/35/, lmaxs/36/, nextv/47/, nreduc/6/, preduc/7/,
     3     radfac/16/, radius/8/, rad0/9/, reldx/17/, tuner4/29/,
     4     tuner5/30/
c/7
c     parameter (dgnorm=1, dinit=38, dstnrm=2, dst0=3, f=10, f0=13,
c    1           fdif=11, gthg=44, gtstep=4, incfac=23, lmat=42,
c    2           lmax0=35, lmaxs=36, nextv=47, nreduc=6, preduc=7,
c    3           radfac=16, radius=8, rad0=9, reldx=17, tuner4=29,
c    4           tuner5=30)
c/
c
c/6
      data half/0.5d+0/, negone/-1.d+0/, one/1.d+0/, onep2/1.2d+0/,
     1     zero/0.d+0/
c/7
c     parameter (half=0.5d+0, negone=-1.d+0, one=1.d+0, onep2=1.2d+0,
c    1           zero=0.d+0)
c/
c
c+++++++++++++++++++++++++++++++  body  ++++++++++++++++++++++++++++++++
c
      i = iv(1)
      if (i .eq. 1) go to 50
      if (i .eq. 2) go to 60
c
c  ***  check validity of iv and v input values  ***
c
      if (iv(1) .eq. 0) call deflt(2, iv, liv, lv, v)
      if (iv(1) .eq. 12 .or. iv(1) .eq. 13)
     1     iv(vneed) = iv(vneed) + n*(n+13)/2
      call parck(2, d, iv, liv, lv, n, v)
      i = iv(1) - 2
      if (i .gt. 12) go to 999
      go to (180, 180, 180, 180, 180, 180, 120, 90, 120, 10, 10, 20), i
c
c  ***  storage allocation  ***
c
10    l = iv(lmat)
      iv(x0) = l + n*(n+1)/2
      iv(step) = iv(x0) + n
      iv(stlstg) = iv(step) + n
      iv(g0) = iv(stlstg) + n
      iv(nwtstp) = iv(g0) + n
      iv(dg) = iv(nwtstp) + n
      iv(nextv) = iv(dg) + n
      if (iv(1) .ne. 13) go to 20
         iv(1) = 14
         go to 999
c
c  ***  initialization  ***
c
 20   iv(niter) = 0
      iv(nfcall) = 1
      iv(ngcall) = 1
      iv(nfgcal) = 1
      iv(mode) = -1
      iv(model) = 1
      iv(stglim) = 1
      iv(toobig) = 0
      iv(cnvcod) = 0
      iv(radinc) = 0
      v(rad0) = zero
      if (v(dinit) .ge. zero) call vscopy(n, d, v(dinit))
      if (iv(inith) .ne. 1) go to 40
c
c     ***  set the initial hessian approximation to diag(d)**-2  ***
c
         l = iv(lmat)
         call vscopy(n*(n+1)/2, v(l), zero)
         k = l - 1
         do 30 i = 1, n
              k = k + i
              t = d(i)
              if (t .le. zero) t = one
              v(k) = t
 30           continue
c
c  ***  compute initial function value  ***
c
 40   iv(1) = 1
      go to 999
c
 50   v(f) = fx
      if (iv(mode) .ge. 0) go to 180
      iv(1) = 2
      if (iv(toobig) .eq. 0) go to 999
         iv(1) = 63
         go to 300
c
c  ***  make sure gradient could be computed  ***
c
 60   if (iv(nfgcal) .ne. 0) go to 70
         iv(1) = 65
         go to 300
c
 70   dg1 = iv(dg)
      call vvmulp(n, v(dg1), g, d, -1)
      v(dgnorm) = v2norm(n, v(dg1))
c
      if (iv(cnvcod) .ne. 0) go to 290
      if (iv(mode) .eq. 0) go to 250
c
c  ***  allow first step to have scaled 2-norm at most v(lmax0)  ***
c
      v(radius) = v(lmax0)
c
      iv(mode) = 0
c
c
c-----------------------------  main loop  -----------------------------
c
c
c  ***  print iteration summary, check iteration limit  ***
c
 80   call itsum(d, g, iv, liv, lv, n, v, x)
 90   k = iv(niter)
      if (k .lt. iv(mxiter)) go to 100
         iv(1) = 10
         go to 300
c
c  ***  update radius  ***
c
 100  iv(niter) = k + 1
      if(k.gt.0)v(radius) = v(radfac) * v(dstnrm)
c
c  ***  initialize for start of next iteration  ***
c
      g01 = iv(g0)
      x01 = iv(x0)
      v(f0) = v(f)
      iv(irc) = 4
      iv(kagqt) = -1
c
c     ***  copy x to x0, g to g0  ***
c
      call vcopy(n, v(x01), x)
      call vcopy(n, v(g01), g)
c
c  ***  check stopx and function evaluation limit  ***
c
 110  if (.not. stopx(dummy)) go to 130
         iv(1) = 11
         go to 140
c
c     ***  come here when restarting after func. eval. limit or stopx.
c
 120  if (v(f) .ge. v(f0)) go to 130
         v(radfac) = one
         k = iv(niter)
         go to 100
c
 130  if (iv(nfcall) .lt. iv(mxfcal)) go to 150
         iv(1) = 9
 140     if (v(f) .ge. v(f0)) go to 300
c
c        ***  in case of stopx or function evaluation limit with
c        ***  improved v(f), evaluate the gradient at x.
c
              iv(cnvcod) = iv(1)
              go to 240
c
c. . . . . . . . . . . . .  compute candidate step  . . . . . . . . . .
c
 150  step1 = iv(step)
      dg1 = iv(dg)
      nwtst1 = iv(nwtstp)
      if (iv(kagqt) .ge. 0) go to 160
         l = iv(lmat)
         call livmul(n, v(nwtst1), v(l), g)
         v(nreduc) = half * dotprd(n, v(nwtst1), v(nwtst1))
         call litvmu(n, v(nwtst1), v(l), v(nwtst1))
         call vvmulp(n, v(step1), v(nwtst1), d, 1)
         v(dst0) = v2norm(n, v(step1))
         call vvmulp(n, v(dg1), v(dg1), d, -1)
         call ltvmul(n, v(step1), v(l), v(dg1))
         v(gthg) = v2norm(n, v(step1))
         iv(kagqt) = 0
 160  call dbdog(v(dg1), lv, n, v(nwtst1), v(step1), v)
      if (iv(irc) .eq. 6) go to 180
c
c  ***  check whether evaluating f(x0 + step) looks worthwhile  ***
c
      if (v(dstnrm) .le. zero) go to 180
      if (iv(irc) .ne. 5) go to 170
      if (v(radfac) .le. one) go to 170
      if (v(preduc) .le. onep2 * v(fdif)) go to 180
c
c  ***  compute f(x0 + step)  ***
c
 170  x01 = iv(x0)
      step1 = iv(step)
      call vaxpy(n, x, one, v(step1), v(x01))
      iv(nfcall) = iv(nfcall) + 1
      iv(1) = 1
      iv(toobig) = 0
      go to 999
c
c. . . . . . . . . . . . .  assess candidate step  . . . . . . . . . . .
c
 180  x01 = iv(x0)
      v(reldx) = reldst(n, d, x, v(x01))
      call assst(iv, liv, lv, v)
      step1 = iv(step)
      lstgst = iv(stlstg)
      if (iv(restor) .eq. 1) call vcopy(n, x, v(x01))
      if (iv(restor) .eq. 2) call vcopy(n, v(lstgst), v(step1))
      if (iv(restor) .ne. 3) go to 190
         call vcopy(n, v(step1), v(lstgst))
         call vaxpy(n, x, one, v(step1), v(x01))
         v(reldx) = reldst(n, d, x, v(x01))
c
 190  k = iv(irc)
      go to (200,230,230,230,200,210,220,220,220,220,220,220,280,250), k
c
c     ***  recompute step with changed radius  ***
c
 200     v(radius) = v(radfac) * v(dstnrm)
         go to 110
c
c  ***  compute step of length v(lmaxs) for singular convergence test.
c
 210  v(radius) = v(lmaxs)
      go to 150
c
c  ***  convergence or false convergence  ***
c
 220  iv(cnvcod) = k - 4
      if (v(f) .ge. v(f0)) go to 290
         if (iv(xirc) .eq. 14) go to 290
              iv(xirc) = 14
c
c. . . . . . . . . . . .  process acceptable step  . . . . . . . . . . .
c
 230  if (iv(irc) .ne. 3) go to 240
         step1 = iv(step)
         temp1 = iv(stlstg)
c
c     ***  set  temp1 = hessian * step  for use in gradient tests  ***
c
         l = iv(lmat)
         call ltvmul(n, v(temp1), v(l), v(step1))
         call lvmul(n, v(temp1), v(l), v(temp1))
c
c  ***  compute gradient  ***
c
 240  iv(ngcall) = iv(ngcall) + 1
      iv(1) = 2
      go to 999
c
c  ***  initializations -- g0 = g - g0, etc.  ***
c
 250  g01 = iv(g0)
      call vaxpy(n, v(g01), negone, v(g01), g)
      step1 = iv(step)
      temp1 = iv(stlstg)
      if (iv(irc) .ne. 3) go to 270
c
c  ***  set v(radfac) by gradient tests  ***
c
c     ***  set  temp1 = diag(d)**-1 * (hessian*step + (g(x0)-g(x)))  ***
c
         call vaxpy(n, v(temp1), negone, v(g01), v(temp1))
         call vvmulp(n, v(temp1), v(temp1), d, -1)
c
c        ***  do gradient tests  ***
c
         if (v2norm(n, v(temp1)) .le. v(dgnorm) * v(tuner4))
     1                  go to 260
              if (dotprd(n, g, v(step1))
     1                  .ge. v(gtstep) * v(tuner5))  go to 270
 260               v(radfac) = v(incfac)
c
c  ***  update h, loop  ***
c
 270  w = iv(nwtstp)
      z = iv(x0)
      l = iv(lmat)
      call wzbfgs(v(l), n, v(step1), v(w), v(g01), v(z))
c
c     ** use the n-vectors starting at v(step1) and v(g01) for scratch..
      call lupdat(v(temp1), v(step1), v(l), v(g01), v(l), n, v(w), v(z))
      iv(1) = 2
      go to 80
c
c. . . . . . . . . . . . . .  misc. details  . . . . . . . . . . . . . .
c
c  ***  bad parameters to assess  ***
c
 280  iv(1) = 64
      go to 300
c
c  ***  print summary of final iteration and other requested items  ***
c
 290  iv(1) = iv(cnvcod)
      iv(cnvcod) = 0
 300  call itsum(d, g, iv, liv, lv, n, v, x)
c
 999  return
c
c  ***  last line of sumit follows  ***
      end
      subroutine snoit(d, fx, iv, liv, lv, n, v, x)
c
c  ***  iteration driver for smsno...
c  ***  minimize general unconstrained objective function using
c  ***  finite-difference gradients and secant hessian approximations.
c
      integer liv, lv, n
      integer iv(liv)
      double precision d(n), fx, x(n), v(lv)
c     dimension v(77 + n*(n+17)/2)
c
c  ***  purpose  ***
c
c        this routine interacts with subroutine  sumit  in an attempt
c     to find an n-vector  x*  that minimizes the (unconstrained)
c     objective function  fx = f(x)  computed by the caller.  (often
c     the  x*  found is a local minimizer rather than a global one.)
c
c  ***  parameters  ***
c
c        the parameters for snoit are the same as those for sumsl
c     (which see), except that calcf, calcg, uiparm, urparm, and ufparm
c     are omitted, and a parameter  fx  for the objective function
c     value at x is added.  instead of calling calcg to obtain the
c     gradient of the objective function at x, snoit calls sgrad2,
c     which computes an approximation to the gradient by finite
c     (forward and central) differences using the method of ref. 1.
c     the following input component is of interest in this regard
c     (and is not described in sumsl).
c
c v(eta0)..... v(42) is an estimated bound on the relative error in the
c             objective function value computed by calcf...
c                  (true value) = (computed value) * (1 + e),
c             where abs(e) .le. v(eta0).  default = machep * 10**3,
c             where machep is the unit roundoff.
c
c        the output values iv(nfcall) and iv(ngcall) have different
c     meanings for smsno than for sumsl...
c
c iv(nfcall)... iv(6) is the number of calls so far made on calcf (i.e.,
c             function evaluations) excluding those made only for
c             computing gradients.  the input value iv(mxfcal) is a
c             limit on iv(nfcall).
c iv(ngcall)... iv(30) is the number of function evaluations made only
c             for computing gradients.  the total number of function
c             evaluations is thus  iv(nfcall) + iv(ngcall).
c
c  ***  references  ***
c
c 1. stewart, g.w. (1967), a modification of davidon*s minimization
c        method to accept difference approximations of derivatives,
c        j. assoc. comput. mach. 14, pp. 72-83.
c.
c  ***  general  ***
c
c     coded by david m. gay (august 1982).
c
c----------------------------  declarations  ---------------------------
c
      external deflt, dotprd, sgrad2, sumit, vscopy
      double precision dotprd
c
c deflt.... supplies default parameter values.
c dotprd... returns inner product of two vectors.
c sgrad2... computes finite-difference gradient approximation.
c sumit.... reverse-communication routine that does sumsl algorithm.
c vscopy... sets all elements of a vector to a scalar.
c
      integer alpha, g1, i, iv1, j, k, w
      double precision one, zero
c
c  ***  subscripts for iv   ***
c
      integer dtype, eta0, f, g, lmat, nextv, nfcall, nfgcal, ngcall,
     1        niter, sgirc, toobig, vneed
c
c/6
      data dtype/16/, eta0/42/, f/10/, g/28/, lmat/42/, nextv/47/,
     1     nfcall/6/, nfgcal/7/, ngcall/30/, niter/31/, sgirc/57/,
     2     toobig/2/, vneed/4/
c/7
c     parameter (dtype=16, eta0=42, f=10, g=28, lmat=42, nextv=47,
c    1           nfcall=6, nfgcal=7, ngcall=30, niter=31, sgirc=57,
c    2           toobig=2, vneed=4)
c/
c/6
      data one/1.d+0/, zero/0.d+0/
c/7
c     parameter (one=1.d+0, zero=0.d+0)
c/
c
c+++++++++++++++++++++++++++++++  body  ++++++++++++++++++++++++++++++++
c
      iv1 = iv(1)
      if (iv1 .eq. 1) go to 10
      if (iv1 .eq. 2) go to 50
      if (iv(1) .eq. 0) call deflt(2, iv, liv, lv, v)
      iv1 = iv(1)
      if (iv1 .eq. 12 .or. iv1 .eq. 13) iv(vneed) = iv(vneed) + 2*n + 6
      if (iv1 .eq. 14) go to 10
      if (iv1 .gt. 2 .and. iv1 .lt. 12) go to 10
      g1 = 1
      if (iv1 .eq. 12) iv(1) = 13
      go to 20
c
 10   g1 = iv(g)
c
 20   call sumit(d, fx, v(g1), iv, liv, lv, n, v, x)
      if (iv(1) - 2) 999, 30, 70
c
c  ***  compute gradient  ***
c
 30   if (iv(niter) .eq. 0) call vscopy(n, v(g1), zero)
      j = iv(lmat)
      k = g1 - n
      do 40 i = 1, n
         v(k) = dotprd(i, v(j), v(j))
         k = k + 1
         j = j + i
 40      continue
c     ***  undo increment of iv(ngcall) done by sumit  ***
      iv(ngcall) = iv(ngcall) - 1
c     ***  store return code from sgrad2 in iv(sgirc)  ***
      iv(sgirc) = 0
c     ***  x may have been restored, so copy back fx... ***
      fx = v(f)
      go to 60
c
c     ***  gradient loop  ***
c
 50   if (iv(toobig) .eq. 0) go to 60
      iv(nfgcal) = 0
      go to 10
c
 60   g1 = iv(g)
      alpha = g1 - n
      w = alpha - 6
      call sgrad2(v(alpha), d, v(eta0), fx, v(g1), iv(sgirc), n, v(w),x)
      if (iv(sgirc) .eq. 0) go to 10
         iv(ngcall) = iv(ngcall) + 1
         go to 999
c
 70   if (iv(1) .ne. 14) go to 999
c
c  ***  storage allocation  ***
c
      iv(g) = iv(nextv) + n + 6
      iv(nextv) = iv(g) + n
      if (iv1 .ne. 13) go to 10
c
 999  return
c  ***  last card of snoit follows  ***
      end
      subroutine dbdog(dig, lv, n, nwtstp, step, v)
c
c  ***  compute double dogleg step  ***
c
c  ***  parameter declarations  ***
c
      integer lv, n
      double precision dig(n), nwtstp(n), step(n), v(lv)
c
c  ***  purpose  ***
c
c        this subroutine computes a candidate step (for use in an uncon-
c     strained minimization code) by the double dogleg algorithm of
c     dennis and mei (ref. 1), which is a variation on powell*s dogleg
c     scheme (ref. 2, p. 95).
c
c--------------------------  parameter usage  --------------------------
c
c    dig (input) diag(d)**-2 * g -- see algorithm notes.
c      g (input) the current gradient vector.
c     lv (input) length of v.
c      n (input) number of components in  dig, g, nwtstp,  and  step.
c nwtstp (input) negative newton step -- see algorithm notes.
c   step (output) the computed step.
c      v (i/o) values array, the following components of which are
c             used here...
c v(bias)   (input) bias for relaxed newton step, which is v(bias) of
c             the way from the full newton to the fully relaxed newton
c             step.  recommended value = 0.8 .
c v(dgnorm) (input) 2-norm of diag(d)**-1 * g -- see algorithm notes.
c v(dstnrm) (output) 2-norm of diag(d) * step, which is v(radius)
c             unless v(stppar) = 0 -- see algorithm notes.
c v(dst0) (input) 2-norm of diag(d) * nwtstp -- see algorithm notes.
c v(grdfac) (output) the coefficient of  dig  in the step returned --
c             step(i) = v(grdfac)*dig(i) + v(nwtfac)*nwtstp(i).
c v(gthg)   (input) square-root of (dig**t) * (hessian) * dig -- see
c             algorithm notes.
c v(gtstep) (output) inner product between g and step.
c v(nreduc) (output) function reduction predicted for the full newton
c             step.
c v(nwtfac) (output) the coefficient of  nwtstp  in the step returned --
c             see v(grdfac) above.
c v(preduc) (output) function reduction predicted for the step returned.
c v(radius) (input) the trust region radius.  d times the step returned
c             has 2-norm v(radius) unless v(stppar) = 0.
c v(stppar) (output) code telling how step was computed... 0 means a
c             full newton step.  between 0 and 1 means v(stppar) of the
c             way from the newton to the relaxed newton step.  between
c             1 and 2 means a true double dogleg step, v(stppar) - 1 of
c             the way from the relaxed newton to the cauchy step.
c             greater than 2 means 1 / (v(stppar) - 1) times the cauchy
c             step.
c
c-------------------------------  notes  -------------------------------
c
c  ***  algorithm notes  ***
c
c        let  g  and  h  be the current gradient and hessian approxima-
c     tion respectively and let d be the current scale vector.  this
c     routine assumes dig = diag(d)**-2 * g  and  nwtstp = h**-1 * g.
c     the step computed is the same one would get by replacing g and h
c     by  diag(d)**-1 * g  and  diag(d)**-1 * h * diag(d)**-1,
c     computing step, and translating step back to the original
c     variables, i.e., premultiplying it by diag(d)**-1.
c
c  ***  references  ***
c
c 1.  dennis, j.e., and mei, h.h.w. (1979), two new unconstrained opti-
c             mization algorithms which use function and gradient
c             values, j. optim. theory applic. 28, pp. 453-482.
c 2. powell, m.j.d. (1970), a hybrid method for non-linear equations,
c             in numerical methods for non-linear equations, edited by
c             p. rabinowitz, gordon and breach, london.
c
c  ***  general  ***
c
c     coded by david m. gay.
c     this subroutine was written in connection with research supported
c     by the national science foundation under grants mcs-7600324 and
c     mcs-7906671.
c
c------------------------  external quantities  ------------------------
c
c  ***  functions and subroutines called  ***
c
      external dotprd, v2norm
      double precision dotprd, v2norm
c
c dotprd... returns inner product of two vectors.
c v2norm... returns 2-norm of a vector.
c
c  ***  intrinsic functions  ***
c/+
      double precision dsqrt
c/
c--------------------------  local variables  --------------------------
c
      integer i
      double precision cfact, cnorm, ctrnwt, ghinvg, femnsq, gnorm,
     1                 nwtnrm, relax, rlambd, t, t1, t2
      double precision half, one, two, zero
c
c  ***  v subscripts  ***
c
      integer bias, dgnorm, dstnrm, dst0, grdfac, gthg, gtstep,
     1        nreduc, nwtfac, preduc, radius, stppar
c
c  ***  data initializations  ***
c
c/6
      data half/0.5d+0/, one/1.d+0/, two/2.d+0/, zero/0.d+0/
c/7
c     parameter (half=0.5d+0, one=1.d+0, two=2.d+0, zero=0.d+0)
c/
c
c/6
      data bias/43/, dgnorm/1/, dstnrm/2/, dst0/3/, grdfac/45/,
     1     gthg/44/, gtstep/4/, nreduc/6/, nwtfac/46/, preduc/7/,
     2     radius/8/, stppar/5/
c/7
c     parameter (bias=43, dgnorm=1, dstnrm=2, dst0=3, grdfac=45,
c    1           gthg=44, gtstep=4, nreduc=6, nwtfac=46, preduc=7,
c    2           radius=8, stppar=5)
c/
c
c+++++++++++++++++++++++++++++++  body  ++++++++++++++++++++++++++++++++
c
      nwtnrm = v(dst0)
      rlambd = one
      if (nwtnrm .gt. zero) rlambd = v(radius) / nwtnrm
      gnorm = v(dgnorm)
      ghinvg = two * v(nreduc)
      v(grdfac) = zero
      v(nwtfac) = zero
      if (rlambd .lt. one) go to 30
c
c        ***  the newton step is inside the trust region  ***
c
         v(stppar) = zero
         v(dstnrm) = nwtnrm
         v(gtstep) = -ghinvg
         v(preduc) = v(nreduc)
         v(nwtfac) = -one
         do 20 i = 1, n
 20           step(i) = -nwtstp(i)
         go to 999
c
 30   v(dstnrm) = v(radius)
      cfact = (gnorm / v(gthg))**2
c     ***  cauchy step = -cfact * g.
      cnorm = gnorm * cfact
      relax = one - v(bias) * (one - gnorm*cnorm/ghinvg)
      if (rlambd .lt. relax) go to 50
c
c        ***  step is between relaxed newton and full newton steps  ***
c
         v(stppar)  =  one  -  (rlambd - relax) / (one - relax)
         t = -rlambd
         v(gtstep) = t * ghinvg
         v(preduc) = rlambd * (one - half*rlambd) * ghinvg
         v(nwtfac) = t
         do 40 i = 1, n
 40           step(i) = t * nwtstp(i)
         go to 999
c
 50   if (cnorm .lt. v(radius)) go to 70
c
c        ***  the cauchy step lies outside the trust region --
c        ***  step = scaled cauchy step  ***
c
         t = -v(radius) / gnorm
         v(grdfac) = t
         v(stppar) = one  +  cnorm / v(radius)
         v(gtstep) = -v(radius) * gnorm
      v(preduc) = v(radius)*(gnorm - half*v(radius)*(v(gthg)/gnorm)**2)
         do 60 i = 1, n
 60           step(i) = t * dig(i)
         go to 999
c
c     ***  compute dogleg step between cauchy and relaxed newton  ***
c     ***  femur = relaxed newton step minus cauchy step  ***
c
 70   ctrnwt = cfact * relax * ghinvg / gnorm
c     *** ctrnwt = inner prod. of cauchy and relaxed newton steps,
c     *** scaled by gnorm**-1.
      t1 = ctrnwt - gnorm*cfact**2
c     ***  t1 = inner prod. of femur and cauchy step, scaled by
c     ***  gnorm**-1.
      t2 = v(radius)*(v(radius)/gnorm) - gnorm*cfact**2
      t = relax * nwtnrm
      femnsq = (t/gnorm)*t - ctrnwt - t1
c     ***  femnsq = square of 2-norm of femur, scaled by gnorm**-1.
      t = t2 / (t1 + dsqrt(t1**2 + femnsq*t2))
c     ***  dogleg step  =  cauchy step  +  t * femur.
      t1 = (t - one) * cfact
      v(grdfac) = t1
      t2 = -t * relax
      v(nwtfac) = t2
      v(stppar) = two - t
      v(gtstep) = t1*gnorm**2 + t2*ghinvg
      v(preduc) = -t1*gnorm * ((t2 + one)*gnorm)
     1                 - t2 * (one + half*t2)*ghinvg
     2                  - half * (v(gthg)*t1)**2
      do 80 i = 1, n
 80      step(i) = t1*dig(i) + t2*nwtstp(i)
c
 999  return
c  ***  last line of dbdog follows  ***
      end
      subroutine ltvmul(n, x, l, y)
c
c  ***  compute  x = (l**t)*y, where  l  is an  n x n  lower
c  ***  triangular matrix stored compactly by rows.  x and y may
c  ***  occupy the same storage.  ***
c
      integer n
      double precision x(n), l(1), y(n)
c     dimension l(n*(n+1)/2)
      integer i, ij, i0, j
      double precision yi, zero
c/6
      data zero/0.d+0/
c/7
c     parameter (zero=0.d+0)
c/
c
      i0 = 0
      do 20 i = 1, n
         yi = y(i)
         x(i) = zero
         do 10 j = 1, i
              ij = i0 + j
              x(j) = x(j) + yi*l(ij)
 10           continue
         i0 = i0 + i
 20      continue
 999  return
c  ***  last card of ltvmul follows  ***
      end
      subroutine lupdat(beta, gamma, l, lambda, lplus, n, w, z)
c
c  ***  compute lplus = secant update of l  ***
c
c  ***  parameter declarations  ***
c
      integer n
      double precision beta(n), gamma(n), l(1), lambda(n), lplus(1),
     1                 w(n), z(n)
c     dimension l(n*(n+1)/2), lplus(n*(n+1)/2)
c
c--------------------------  parameter usage  --------------------------
c
c   beta = scratch vector.
c  gamma = scratch vector.
c      l (input) lower triangular matrix, stored rowwise.
c lambda = scratch vector.
c  lplus (output) lower triangular matrix, stored rowwise, which may
c             occupy the same storage as  l.
c      n (input) length of vector parameters and order of matrices.
c      w (input, destroyed on output) right singular vector of rank 1
c             correction to  l.
c      z (input, destroyed on output) left singular vector of rank 1
c             correction to  l.
c
c-------------------------------  notes  -------------------------------
c
c  ***  application and usage restrictions  ***
c
c        this routine updates the cholesky factor  l  of a symmetric
c     positive definite matrix to which a secant update is being
c     applied -- it computes a cholesky factor  lplus  of
c     l * (i + z*w**t) * (i + w*z**t) * l**t.  it is assumed that  w
c     and  z  have been chosen so that the updated matrix is strictly
c     positive definite.
c
c  ***  algorithm notes  ***
c
c        this code uses recurrence 3 of ref. 1 (with d(j) = 1 for all j)
c     to compute  lplus  of the form  l * (i + z*w**t) * q,  where  q
c     is an orthogonal matrix that makes the result lower triangular.
c        lplus may have some negative diagonal elements.
c
c  ***  references  ***
c
c 1.  goldfarb, d. (1976), factorized variable metric methods for uncon-
c             strained optimization, math. comput. 30, pp. 796-811.
c
c  ***  general  ***
c
c     coded by david m. gay (fall 1979).
c     this subroutine was written in connection with research supported
c     by the national science foundation under grants mcs-7600324 and
c     mcs-7906671.
c
c------------------------  external quantities  ------------------------
c
c  ***  intrinsic functions  ***
c/+
      double precision dsqrt
c/
c--------------------------  local variables  --------------------------
c
      integer i, ij, j, jj, jp1, k, nm1, np1
      double precision a, b, bj, eta, gj, lj, lij, ljj, nu, s, theta,
     1                 wj, zj
      double precision one, zero
c
c  ***  data initializations  ***
c
c/6
      data one/1.d+0/, zero/0.d+0/
c/7
c     parameter (one=1.d+0, zero=0.d+0)
c/
c
c+++++++++++++++++++++++++++++++  body  ++++++++++++++++++++++++++++++++
c
      nu = one
      eta = zero
      if (n .le. 1) go to 30
      nm1 = n - 1
c
c  ***  temporarily store s(j) = sum over k = j+1 to n of w(k)**2 in
c  ***  lambda(j).
c
      s = zero
      do 10 i = 1, nm1
         j = n - i
         s = s + w(j+1)**2
         lambda(j) = s
 10      continue
c
c  ***  compute lambda, gamma, and beta by goldfarb*s recurrence 3.
c
      do 20 j = 1, nm1
         wj = w(j)
         a = nu*z(j) - eta*wj
         theta = one + a*wj
         s = a*lambda(j)
         lj = dsqrt(theta**2 + a*s)
         if (theta .gt. zero) lj = -lj
         lambda(j) = lj
         b = theta*wj + s
         gamma(j) = b * nu / lj
         beta(j) = (a - b*eta) / lj
         nu = -nu / lj
         eta = -(eta + (a**2)/(theta - lj)) / lj
 20      continue
 30   lambda(n) = one + (nu*z(n) - eta*w(n))*w(n)
c
c  ***  update l, gradually overwriting  w  and  z  with  l*w  and  l*z.
c
      np1 = n + 1
      jj = n * (n + 1) / 2
      do 60 k = 1, n
         j = np1 - k
         lj = lambda(j)
         ljj = l(jj)
         lplus(jj) = lj * ljj
         wj = w(j)
         w(j) = ljj * wj
         zj = z(j)
         z(j) = ljj * zj
         if (k .eq. 1) go to 50
         bj = beta(j)
         gj = gamma(j)
         ij = jj + j
         jp1 = j + 1
         do 40 i = jp1, n
              lij = l(ij)
              lplus(ij) = lj*lij + bj*w(i) + gj*z(i)
              w(i) = w(i) + lij*wj
              z(i) = z(i) + lij*zj
              ij = ij + i
 40           continue
 50      jj = jj - j
 60      continue
c
 999  return
c  ***  last card of lupdat follows  ***
      end
      subroutine lvmul(n, x, l, y)
c
c  ***  compute  x = l*y, where  l  is an  n x n  lower triangular
c  ***  matrix stored compactly by rows.  x and y may occupy the same
c  ***  storage.  ***
c
      integer n
      double precision x(n), l(1), y(n)
c     dimension l(n*(n+1)/2)
      integer i, ii, ij, i0, j, np1
      double precision t, zero
c/6
      data zero/0.d+0/
c/7
c     parameter (zero=0.d+0)
c/
c
      np1 = n + 1
      i0 = n*(n+1)/2
      do 20 ii = 1, n
         i = np1 - ii
         i0 = i0 - i
         t = zero
         do 10 j = 1, i
              ij = i0 + j
              t = t + l(ij)*y(j)
 10           continue
         x(i) = t
 20      continue
 999  return
c  ***  last card of lvmul follows  ***
      end
      subroutine sgrad2 (alpha, d, eta0, fx, g, irc, n, w, x)
c
c  ***  compute finite difference gradient by stweart*s scheme  ***
c
c     ***  parameters  ***
c
      integer irc, n
      double precision alpha(n), d(n), eta0, fx, g(n), w(6), x(n)
c
c.......................................................................
c
c     ***  purpose  ***
c
c        this subroutine uses an embellished form of the finite-differ-
c     ence scheme proposed by stewart (ref. 1) to approximate the
c     gradient of the function f(x), whose values are supplied by
c     reverse communication.
c
c     ***  parameter description  ***
c
c  alpha in  (approximate) diagonal elements of the hessian of f(x).
c      d in  scale vector such that d(i)*x(i), i = 1,...,n, are in
c             comparable units.
c   eta0 in  estimated bound on relative error in the function value...
c             (true value) = (computed value)*(1+e),   where
c             abs(e) .le. eta0.
c     fx i/o on input,  fx  must be the computed value of f(x).  on
c             output with irc = 0, fx has been restored to its original
c             value, the one it had when sgrad2 was last called with
c             irc = 0.
c      g i/o on input with irc = 0, g should contain an approximation
c             to the gradient of f near x, e.g., the gradient at the
c             previous iterate.  when sgrad2 returns with irc = 0, g is
c             the desired finite-difference approximation to the
c             gradient at x.
c    irc i/o input/return code... before the very first call on sgrad2,
c             the caller must set irc to 0.  whenever sgrad2 returns a
c             nonzero value for irc, it has perturbed some component of
c             x... the caller should evaluate f(x) and call sgrad2
c             again with fx = f(x).
c      n in  the number of variables (components of x) on which f
c             depends.
c      x i/o on input with irc = 0, x is the point at which the
c             gradient of f is desired.  on output with irc nonzero, x
c             is the point at which f should be evaluated.  on output
c             with irc = 0, x has been restored to its original value
c             (the one it had when sgrad2 was last called with irc = 0)
c             and g contains the desired gradient approximation.
c      w i/o work vector of length 6 in which sgrad2 saves certain
c             quantities while the caller is evaluating f(x) at a
c             perturbed x.
c
c     ***  application and usage restrictions  ***
c
c        this routine is intended for use with quasi-newton routines
c     for unconstrained minimization (in which case  alpha  comes from
c     the diagonal of the quasi-newton hessian approximation).
c
c     ***  algorithm notes  ***
c
c        this code departs from the scheme proposed by stewart (ref. 1)
c     in its guarding against overly large or small step sizes and its
c     handling of special cases (such as zero components of alpha or g).
c
c     ***  references  ***
c
c 1. stewart, g.w. (1967), a modification of davidon*s minimization
c        method to accept difference approximations of derivatives,
c        j. assoc. comput. mach. 14, pp. 72-83.
c
c     ***  history  ***
c
c     designed and coded by david m. gay (summer 1977/summer 1980).
c
c     ***  general  ***
c
c        this routine was prepared in connection with work supported by
c     the national science foundation under grants mcs76-00324 and
c     mcs-7906671.
c
c.......................................................................
c
c     *****  external function  *****
c
      external rmdcon
      double precision rmdcon
c rmdcon... returns machine-dependent constants.
c
c     ***** intrinsic functions *****
c/+
      integer iabs
      double precision dabs, dmax1, dsqrt
c/
c     ***** local variables *****
c
      integer fh, fx0, hsave, i, xisave
      double precision aai, afx, afxeta, agi, alphai, axi, axibar,
     1                 discon, eta, gi, h, hmin
      double precision c2000, four, hmax0, hmin0, h0, machep, one, p002,
     1                 three, two, zero
c
c/6
      data c2000/2.0d+3/, four/4.0d+0/, hmax0/0.02d+0/, hmin0/5.0d+1/,
     1     one/1.0d+0/, p002/0.002d+0/, three/3.0d+0/,
     2     two/2.0d+0/, zero/0.0d+0/
c/7
c     parameter (c2000=2.0d+3, four=4.0d+0, hmax0=0.02d+0, hmin0=5.0d+1,
c    1     one=1.0d+0, p002=0.002d+0, three=3.0d+0,
c    2     two=2.0d+0, zero=0.0d+0)
c/
c/6
      data fh/3/, fx0/4/, hsave/5/, xisave/6/
c/7
c     parameter (fh=3, fx0=4, hsave=5, xisave=6)
c/
c
c---------------------------------  body  ------------------------------
c
      if (irc) 140, 100, 210
c
c     ***  fresh start -- get machine-dependent constants  ***
c
c     store machep in w(1) and h0 in w(2), where machep is the unit
c     roundoff (the smallest positive number such that
c     1 + machep .gt. 1  and  1 - machep .lt. 1),  and  h0 is the
c     square-root of machep.
c
 100  w(1) = rmdcon(3)
      w(2) = dsqrt(w(1))
c
      w(fx0) = fx
c
c     ***  increment  i  and start computing  g(i)  ***
c
 110  i = iabs(irc) + 1
      if (i .gt. n) go to 300
         irc = i
         afx = dabs(w(fx0))
         machep = w(1)
         h0 = w(2)
         hmin = hmin0 * machep
         w(xisave) = x(i)
         axi = dabs(x(i))
         axibar = dmax1(axi, one/d(i))
         gi = g(i)
         agi = dabs(gi)
         eta = dabs(eta0)
         if (afx .gt. zero) eta = dmax1(eta, agi*axi*machep/afx)
         alphai = alpha(i)
         if (alphai .eq. zero) go to 170
         if (gi .eq. zero .or. fx .eq. zero) go to 180
         afxeta = afx*eta
         aai = dabs(alphai)
c
c        *** compute h = stewart*s forward-difference step size.
c
         if (gi**2 .le. afxeta*aai) go to 120
              h = two*dsqrt(afxeta/aai)
              h = h*(one - aai*h/(three*aai*h + four*agi))
              go to 130
 120     h = two*(afxeta*agi/(aai**2))**(one/three)
         h = h*(one - two*agi/(three*aai*h + four*agi))
c
c        ***  ensure that  h  is not insignificantly small  ***
c
 130     h = dmax1(h, hmin*axibar)
c
c        *** use forward difference if bound on truncation error is at
c        *** most 10**-3.
c
         if (aai*h .le. p002*agi) go to 160
c
c        *** compute h = stewart*s step for central difference.
c
         discon = c2000*afxeta
         h = discon/(agi + dsqrt(gi**2 + aai*discon))
c
c        ***  ensure that  h  is neither too small nor too big  ***
c
         h = dmax1(h, hmin*axibar)
         if (h .ge. hmax0*axibar) h = axibar * h0**(two/three)
c
c        ***  compute central difference  ***
c
         irc = -i
         go to 200
c
 140     h = -w(hsave)
         i = iabs(irc)
         if (h .gt. zero) go to 150
         w(fh) = fx
         go to 200
c
 150     g(i) = (w(fh) - fx) / (two * h)
         x(i) = w(xisave)
         go to 110
c
c     ***  compute forward differences in various cases  ***
c
 160     if (h .ge. hmax0*axibar) h = h0 * axibar
         if (alphai*gi .lt. zero) h = -h
         go to 200
 170     h = axibar
         go to 200
 180     h = h0 * axibar
c
 200     x(i) = w(xisave) + h
         w(hsave) = h
         go to 999
c
c     ***  compute actual forward difference  ***
c
 210     g(irc) = (fx - w(fx0)) / w(hsave)
         x(irc) = w(xisave)
         go to 110
c
c  ***  restore fx and indicate that g has been computed  ***
c
 300  fx = w(fx0)
      irc = 0
c
 999  return
c  ***  last card of sgrad2 follows  ***
      end
      subroutine vvmulp(n, x, y, z, k)
c
c ***  set x(i) = y(i) * z(i)**k, 1 .le. i .le. n (for k = 1 or -1)  ***
c
      integer n, k
      double precision x(n), y(n), z(n)
      integer i
c
      if (k .ge. 0) go to 20
      do 10 i = 1, n
 10      x(i) = y(i) / z(i)
      go to 999
c
 20   do 30 i = 1, n
 30      x(i) = y(i) * z(i)
 999  return
c  ***  last card of vvmulp follows  ***
      end
      subroutine wzbfgs (l, n, s, w, y, z)
c
c  ***  compute  y  and  z  for  lupdat  corresponding to bfgs update.
c
      integer n
      double precision l(1), s(n), w(n), y(n), z(n)
c     dimension l(n*(n+1)/2)
c
c--------------------------  parameter usage  --------------------------
c
c l (i/o) cholesky factor of hessian, a lower triang. matrix stored
c             compactly by rows.
c n (input) order of  l  and length of  s,  w,  y,  z.
c s (input) the step just taken.
c w (output) right singular vector of rank 1 correction to l.
c y (input) change in gradients corresponding to s.
c z (output) left singular vector of rank 1 correction to l.
c
c-------------------------------  notes  -------------------------------
c
c  ***  algorithm notes  ***
c
c        when  s  is computed in certain ways, e.g. by  gqtstp  or
c     dbldog,  it is possible to save n**2/2 operations since  (l**t)*s
c     or  l*(l**t)*s is then known.
c        if the bfgs update to l*(l**t) would reduce its determinant to
c     less than eps times its old value, then this routine in effect
c     replaces  y  by  theta*y + (1 - theta)*l*(l**t)*s,  where  theta
c     (between 0 and 1) is chosen to make the reduction factor = eps.
c
c  ***  general  ***
c
c     coded by david m. gay (fall 1979).
c     this subroutine was written in connection with research supported
c     by the national science foundation under grants mcs-7600324 and
c     mcs-7906671.
c
c------------------------  external quantities  ------------------------
c
c  ***  functions and subroutines called  ***
c
      external dotprd, livmul, ltvmul
      double precision dotprd
c dotprd returns inner product of two vectors.
c livmul multiplies l**-1 times a vector.
c ltvmul multiplies l**t times a vector.
c
c  ***  intrinsic functions  ***
c/+
      double precision dsqrt
c/
c--------------------------  local variables  --------------------------
c
      integer i
      double precision cs, cy, eps, epsrt, one, shs, ys, theta
c
c  ***  data initializations  ***
c
c/6
      data eps/0.1d+0/, one/1.d+0/
c/7
c     parameter (eps=0.1d+0, one=1.d+0)
c/
c
c+++++++++++++++++++++++++++++++  body  ++++++++++++++++++++++++++++++++
c
      call ltvmul(n, w, l, s)
      shs = dotprd(n, w, w)
      ys = dotprd(n, y, s)
      if (ys .ge. eps*shs) go to 10
         theta = (one - eps) * shs / (shs - ys)
         epsrt = dsqrt(eps)
         cy = theta / (shs * epsrt)
         cs = (one + (theta-one)/epsrt) / shs
         go to 20
 10   cy = one / (dsqrt(ys) * dsqrt(shs))
      cs = one / shs
 20   call livmul(n, z, l, y)
      do 30 i = 1, n
 30      z(i) = cy * z(i)  -  cs * w(i)
c
 999  return
c  ***  last card of wzbfgs follows  ***
      end
      subroutine assst(iv, liv, lv, v)
c
c  ***  assess candidate step (***sol version 2.3)  ***
c
      integer liv, lv
      integer iv(liv)
      double precision v(lv)
c
c  ***  purpose  ***
c
c        this subroutine is called by an unconstrained minimization
c     routine to assess the next candidate step.  it may recommend one
c     of several courses of action, such as accepting the step, recom-
c     puting it using the same or a new quadratic model, or halting due
c     to convergence or false convergence.  see the return code listing
c     below.
c
c--------------------------  parameter usage  --------------------------
c
c  iv (i/o) integer parameter and scratch vector -- see description
c             below of iv values referenced.
c liv (in)  length of iv array.
c  lv (in)  length of v array.
c   v (i/o) real parameter and scratch vector -- see description
c             below of v values referenced.
c
c  ***  iv values referenced  ***
c
c    iv(irc) (i/o) on input for the first step tried in a new iteration,
c             iv(irc) should be set to 3 or 4 (the value to which it is
c             set when step is definitely to be accepted).  on input
c             after step has been recomputed, iv(irc) should be
c             unchanged since the previous return of assst.
c                on output, iv(irc) is a return code having one of the
c             following values...
c                  1 = switch models or try smaller step.
c                  2 = switch models or accept step.
c                  3 = accept step and determine v(radfac) by gradient
c                       tests.
c                  4 = accept step, v(radfac) has been determined.
c                  5 = recompute step (using the same model).
c                  6 = recompute step with radius = v(lmaxs) but do not
c                       evaulate the objective function.
c                  7 = x-convergence (see v(xctol)).
c                  8 = relative function convergence (see v(rfctol)).
c                  9 = both x- and relative function convergence.
c                 10 = absolute function convergence (see v(afctol)).
c                 11 = singular convergence (see v(lmaxs)).
c                 12 = false convergence (see v(xftol)).
c                 13 = iv(irc) was out of range on input.
c             return code i has precdence over i+1 for i = 9, 10, 11.
c iv(mlstgd) (i/o) saved value of iv(model).
c  iv(model) (i/o) on input, iv(model) should be an integer identifying
c             the current quadratic model of the objective function.
c             if a previous step yielded a better function reduction,
c             then iv(model) will be set to iv(mlstgd) on output.
c iv(nfcall) (in)  invocation count for the objective function.
c iv(nfgcal) (i/o) value of iv(nfcall) at step that gave the biggest
c             function reduction this iteration.  iv(nfgcal) remains
c             unchanged until a function reduction is obtained.
c iv(radinc) (i/o) the number of radius increases (or minus the number
c             of decreases) so far this iteration.
c iv(restor) (out) set to 1 if v(f) has been restored and x should be
c             restored to its initial value, to 2 if x should be saved,
c             to 3 if x should be restored from the saved value, and to
c             0 otherwise.
c  iv(stage) (i/o) count of the number of models tried so far in the
c             current iteration.
c iv(stglim) (in)  maximum number of models to consider.
c iv(switch) (out) set to 0 unless a new model is being tried and it
c             gives a smaller function value than the previous model,
c             in which case assst sets iv(switch) = 1.
c iv(toobig) (in)  is nonzero if step was too big (e.g. if it caused
c             overflow).
c   iv(xirc) (i/o) value that iv(irc) would have in the absence of
c             convergence, false convergence, and oversized steps.
c
c  ***  v values referenced  ***
c
c v(afctol) (in)  absolute function convergence tolerance.  if the
c             absolute value of the current function value v(f) is less
c             than v(afctol), then assst returns with iv(irc) = 10.
c v(decfac) (in)  factor by which to decrease radius when iv(toobig) is
c             nonzero.
c v(dstnrm) (in)  the 2-norm of d*step.
c v(dstsav) (i/o) value of v(dstnrm) on saved step.
c   v(dst0) (in)  the 2-norm of d times the newton step (when defined,
c             i.e., for v(nreduc) .ge. 0).
c      v(f) (i/o) on both input and output, v(f) is the objective func-
c             tion value at x.  if x is restored to a previous value,
c             then v(f) is restored to the corresponding value.
c   v(fdif) (out) the function reduction v(f0) - v(f) (for the output
c             value of v(f) if an earlier step gave a bigger function
c             decrease, and for the input value of v(f) otherwise).
c v(flstgd) (i/o) saved value of v(f).
c     v(f0) (in)  objective function value at start of iteration.
c v(gtslst) (i/o) value of v(gtstep) on saved step.
c v(gtstep) (in)  inner product between step and gradient.
c v(incfac) (in)  minimum factor by which to increase radius.
c  v(lmaxs) (in)  maximum reasonable step size (and initial step bound).
c             if the actual function decrease is no more than twice
c             what was predicted, if a return with iv(irc) = 7, 8, 9,
c             or 10 does not occur, if v(dstnrm) .gt. v(lmaxs), and if
c             v(preduc) .le. v(sctol) * abs(v(f0)), then assst re-
c             turns with iv(irc) = 11.  if so doing appears worthwhile,
c             then assst repeats this test with v(preduc) computed for
c             a step of length v(lmaxs) (by a return with iv(irc) = 6).
c v(nreduc) (i/o)  function reduction predicted by quadratic model for
c             newton step.  if assst is called with iv(irc) = 6, i.e.,
c             if v(preduc) has been computed with radius = v(lmaxs) for
c             use in the singular convervence test, then v(nreduc) is
c             set to -v(preduc) before the latter is restored.
c v(plstgd) (i/o) value of v(preduc) on saved step.
c v(preduc) (i/o) function reduction predicted by quadratic model for
c             current step.
c v(radfac) (out) factor to be used in determining the new radius,
c             which should be v(radfac)*dst, where  dst  is either the
c             output value of v(dstnrm) or the 2-norm of
c             diag(newd)*step  for the output value of step and the
c             updated version, newd, of the scale vector d.  for
c             iv(irc) = 3, v(radfac) = 1.0 is returned.
c v(rdfcmn) (in)  minimum value for v(radfac) in terms of the input
c             value of v(dstnrm) -- suggested value = 0.1.
c v(rdfcmx) (in)  maximum value for v(radfac) -- suggested value = 4.0.
c  v(reldx) (in) scaled relative change in x caused by step, computed
c             (e.g.) by function  reldst  as
c                 max (d(i)*abs(x(i)-x0(i)), 1 .le. i .le. p) /
c                    max (d(i)*(abs(x(i))+abs(x0(i))), 1 .le. i .le. p).
c v(rfctol) (in)  relative function convergence tolerance.  if the
c             actual function reduction is at most twice what was pre-
c             dicted and  v(nreduc) .le. v(rfctol)*abs(v(f0)),  then
c             assst returns with iv(irc) = 8 or 9.
c v(stppar) (in)  marquardt parameter -- 0 means full newton step.
c v(tuner1) (in)  tuning constant used to decide if the function
c             reduction was much less than expected.  suggested
c             value = 0.1.
c v(tuner2) (in)  tuning constant used to decide if the function
c             reduction was large enough to accept step.  suggested
c             value = 10**-4.
c v(tuner3) (in)  tuning constant used to decide if the radius
c             should be increased.  suggested value = 0.75.
c  v(xctol) (in)  x-convergence criterion.  if step is a newton step
c             (v(stppar) = 0) having v(reldx) .le. v(xctol) and giving
c             at most twice the predicted function decrease, then
c             assst returns iv(irc) = 7 or 9.
c  v(xftol) (in)  false convergence tolerance.  if step gave no or only
c             a small function decrease and v(reldx) .le. v(xftol),
c             then assst returns with iv(irc) = 12.
c
c-------------------------------  notes  -------------------------------
c
c  ***  application and usage restrictions  ***
c
c        this routine is called as part of the nl2sol (nonlinear
c     least-squares) package.  it may be used in any unconstrained
c     minimization solver that uses dogleg, goldfeld-quandt-trotter,
c     or levenberg-marquardt steps.
c
c  ***  algorithm notes  ***
c
c        see (1) for further discussion of the assessing and model
c     switching strategies.  while nl2sol considers only two models,
c     assst is designed to handle any number of models.
c
c  ***  usage notes  ***
c
c        on the first call of an iteration, only the i/o variables
c     step, x, iv(irc), iv(model), v(f), v(dstnrm), v(gtstep), and
c     v(preduc) need have been initialized.  between calls, no i/o
c     values execpt step, x, iv(model), v(f) and the stopping toler-
c     ances should be changed.
c        after a return for convergence or false convergence, one can
c     change the stopping tolerances and call assst again, in which
c     case the stopping tests will be repeated.
c
c  ***  references  ***
c
c     (1) dennis, j.e., jr., gay, d.m., and welsch, r.e. (1981),
c        an adaptive nonlinear least-squares algorithm,
c        acm trans. math. software, vol. 7, no. 3.
c
c     (2) powell, m.j.d. (1970)  a fortran subroutine for solving
c        systems of nonlinear algebraic equations, in numerical
c        methods for nonlinear algebraic equations, edited by
c        p. rabinowitz, gordon and breach, london.
c
c  ***  history  ***
c
c        john dennis designed much of this routine, starting with
c     ideas in (2). roy welsch suggested the model switching strategy.
c        david gay and stephen peters cast this subroutine into a more
c     portable form (winter 1977), and david gay cast it into its
c     present form (fall 1978).
c
c  ***  general  ***
c
c     this subroutine was written in connection with research
c     supported by the national science foundation under grants
c     mcs-7600324, dcr75-10143, 76-14311dss, mcs76-11989, and
c     mcs-7906671.
c
c------------------------  external quantities  ------------------------
c
c  ***  no external functions and subroutines  ***
c
c  ***  intrinsic functions  ***
c/+
      double precision dabs, dmax1
c/
c  ***  no common blocks  ***
c
c--------------------------  local variables  --------------------------
c
      logical goodx
      integer i, nfc
      double precision emax, emaxs, gts, rfac1, xmax
      double precision half, one, onep2, two, zero
c
c  ***  subscripts for iv and v  ***
c
      integer afctol, decfac, dstnrm, dstsav, dst0, f, fdif, flstgd, f0,
     1        gtslst, gtstep, incfac, irc, lmaxs, mlstgd, model, nfcall,
     2        nfgcal, nreduc, plstgd, preduc, radfac, radinc, rdfcmn,
     3        rdfcmx, reldx, restor, rfctol, sctol, stage, stglim,
     4        stppar, switch, toobig, tuner1, tuner2, tuner3, xctol,
     5        xftol, xirc
c
c  ***  data initializations  ***
c
c/6
      data half/0.5d+0/, one/1.d+0/, onep2/1.2d+0/, two/2.d+0/,
     1     zero/0.d+0/
c/7
c     parameter (half=0.5d+0, one=1.d+0, onep2=1.2d+0, two=2.d+0,
c    1           zero=0.d+0)
c/
c
c/6
      data irc/29/, mlstgd/32/, model/5/, nfcall/6/, nfgcal/7/,
     1     radinc/8/, restor/9/, stage/10/, stglim/11/, switch/12/,
     2     toobig/2/, xirc/13/
c/7
c     parameter (irc=29, mlstgd=32, model=5, nfcall=6, nfgcal=7,
c    1           radinc=8, restor=9, stage=10, stglim=11, switch=12,
c    2           toobig=2, xirc=13)
c/
c/6
      data afctol/31/, decfac/22/, dstnrm/2/, dst0/3/, dstsav/18/,
     1     f/10/, fdif/11/, flstgd/12/, f0/13/, gtslst/14/, gtstep/4/,
     2     incfac/23/, lmaxs/36/, nreduc/6/, plstgd/15/, preduc/7/,
     3     radfac/16/, rdfcmn/24/, rdfcmx/25/, reldx/17/, rfctol/32/,
     4     sctol/37/, stppar/5/, tuner1/26/, tuner2/27/, tuner3/28/,
     5     xctol/33/, xftol/34/
c/7
c     parameter (afctol=31, decfac=22, dstnrm=2, dst0=3, dstsav=18,
c    1           f=10, fdif=11, flstgd=12, f0=13, gtslst=14, gtstep=4,
c    2           incfac=23, lmaxs=36, nreduc=6, plstgd=15, preduc=7,
c    3           radfac=16, rdfcmn=24, rdfcmx=25, reldx=17, rfctol=32,
c    4           sctol=37, stppar=5, tuner1=26, tuner2=27, tuner3=28,
c    5           xctol=33, xftol=34)
c/
c
c+++++++++++++++++++++++++++++++  body  ++++++++++++++++++++++++++++++++
c
      nfc = iv(nfcall)
      iv(switch) = 0
      iv(restor) = 0
      rfac1 = one
      goodx = .true.
      i = iv(irc)
      if (i .ge. 1 .and. i .le. 12)
     1             go to (20,30,10,10,40,280,220,220,220,220,220,170), i
         iv(irc) = 13
         go to 999
c
c  ***  initialize for new iteration  ***
c
 10   iv(stage) = 1
      iv(radinc) = 0
      v(flstgd) = v(f0)
      if (iv(toobig) .eq. 0) go to 110
         iv(stage) = -1
         iv(xirc) = i
         go to 60
c
c  ***  step was recomputed with new model or smaller radius  ***
c  ***  first decide which  ***
c
 20   if (iv(model) .ne. iv(mlstgd)) go to 30
c        ***  old model retained, smaller radius tried  ***
c        ***  do not consider any more new models this iteration  ***
         iv(stage) = iv(stglim)
         iv(radinc) = -1
         go to 110
c
c  ***  a new model is being tried.  decide whether to keep it.  ***
c
 30   iv(stage) = iv(stage) + 1
c
c     ***  now we add the possibiltiy that step was recomputed with  ***
c     ***  the same model, perhaps because of an oversized step.     ***
c
 40   if (iv(stage) .gt. 0) go to 50
c
c        ***  step was recomputed because it was too big.  ***
c
         if (iv(toobig) .ne. 0) go to 60
c
c        ***  restore iv(stage) and pick up where we left off.  ***
c
         iv(stage) = -iv(stage)
         i = iv(xirc)
         go to (20, 30, 110, 110, 70), i
c
 50   if (iv(toobig) .eq. 0) go to 70
c
c  ***  handle oversize step  ***
c
      if (iv(radinc) .gt. 0) go to 80
         iv(stage) = -iv(stage)
         iv(xirc) = iv(irc)
c
 60      v(radfac) = v(decfac)
         iv(radinc) = iv(radinc) - 1
         iv(irc) = 5
         iv(restor) = 1
         go to 999
c
 70   if (v(f) .lt. v(flstgd)) go to 110
c
c     *** the new step is a loser.  restore old model.  ***
c
      if (iv(model) .eq. iv(mlstgd)) go to 80
         iv(model) = iv(mlstgd)
         iv(switch) = 1
c
c     ***  restore step, etc. only if a previous step decreased v(f).
c
 80   if (v(flstgd) .ge. v(f0)) go to 110
         iv(restor) = 1
         v(f) = v(flstgd)
         v(preduc) = v(plstgd)
         v(gtstep) = v(gtslst)
         if (iv(switch) .eq. 0) rfac1 = v(dstnrm) / v(dstsav)
         v(dstnrm) = v(dstsav)
         nfc = iv(nfgcal)
         goodx = .false.
c
 110  v(fdif) = v(f0) - v(f)
      if (v(fdif) .gt. v(tuner2) * v(preduc)) go to 140
      if(iv(radinc).gt.0) go to 140
c
c        ***  no (or only a trivial) function decrease
c        ***  -- so try new model or smaller radius
c
         if (v(f) .lt. v(f0)) go to 120
              iv(mlstgd) = iv(model)
              v(flstgd) = v(f)
              v(f) = v(f0)
              iv(restor) = 1
              go to 130
 120     iv(nfgcal) = nfc
 130     iv(irc) = 1
         if (iv(stage) .lt. iv(stglim)) go to 160
              iv(irc) = 5
              iv(radinc) = iv(radinc) - 1
              go to 160
c
c  ***  nontrivial function decrease achieved  ***
c
 140  iv(nfgcal) = nfc
      rfac1 = one
      v(dstsav) = v(dstnrm)
      if (v(fdif) .gt. v(preduc)*v(tuner1)) go to 190
c
c  ***  decrease was much less than predicted -- either change models
c  ***  or accept step with decreased radius.
c
      if (iv(stage) .ge. iv(stglim)) go to 150
c        ***  consider switching models  ***
         iv(irc) = 2
         go to 160
c
c     ***  accept step with decreased radius  ***
c
 150  iv(irc) = 4
c
c  ***  set v(radfac) to fletcher*s decrease factor  ***
c
 160  iv(xirc) = iv(irc)
      emax = v(gtstep) + v(fdif)
      v(radfac) = half * rfac1
      if (emax .lt. v(gtstep)) v(radfac) = rfac1 * dmax1(v(rdfcmn),
     1                                           half * v(gtstep)/emax)
c
c  ***  do false convergence test  ***
c
 170  if (v(reldx) .le. v(xftol)) go to 180
         iv(irc) = iv(xirc)
         if (v(f) .lt. v(f0)) go to 200
              go to 230
c
 180  iv(irc) = 12
      go to 240
c
c  ***  handle good function decrease  ***
c
 190  if (v(fdif) .lt. (-v(tuner3) * v(gtstep))) go to 210
c
c     ***  increasing radius looks worthwhile.  see if we just
c     ***  recomputed step with a decreased radius or restored step
c     ***  after recomputing it with a larger radius.
c
      if (iv(radinc) .lt. 0) go to 210
      if (iv(restor) .eq. 1) go to 210
c
c        ***  we did not.  try a longer step unless this was a newton
c        ***  step.
c
         v(radfac) = v(rdfcmx)
         gts = v(gtstep)
         if (v(fdif) .lt. (half/v(radfac) - one) * gts)
     1            v(radfac) = dmax1(v(incfac), half*gts/(gts + v(fdif)))
         iv(irc) = 4
         if (v(stppar) .eq. zero) go to 230
         if (v(dst0) .ge. zero .and. (v(dst0) .lt. two*v(dstnrm)
     1             .or. v(nreduc) .lt. onep2*v(fdif)))  go to 230
c             ***  step was not a newton step.  recompute it with
c             ***  a larger radius.
              iv(irc) = 5
              iv(radinc) = iv(radinc) + 1
c
c  ***  save values corresponding to good step  ***
c
 200  v(flstgd) = v(f)
      iv(mlstgd) = iv(model)
      if (iv(restor) .ne. 1) iv(restor) = 2
      v(dstsav) = v(dstnrm)
      iv(nfgcal) = nfc
      v(plstgd) = v(preduc)
      v(gtslst) = v(gtstep)
      go to 230
c
c  ***  accept step with radius unchanged  ***
c
 210  v(radfac) = one
      iv(irc) = 3
      go to 230
c
c  ***  come here for a restart after convergence  ***
c
 220  iv(irc) = iv(xirc)
      if (v(dstsav) .ge. zero) go to 240
         iv(irc) = 12
         go to 240
c
c  ***  perform convergence tests  ***
c
 230  iv(xirc) = iv(irc)
 240  if (iv(restor) .eq. 1 .and. v(flstgd) .lt. v(f0)) iv(restor) = 3
      if (dabs(v(f)) .lt. v(afctol)) iv(irc) = 10
      if (half * v(fdif) .gt. v(preduc)) go to 999
      emax = v(rfctol) * dabs(v(f0))
      emaxs = v(sctol) * dabs(v(f0))
      if (v(dstnrm) .gt. v(lmaxs) .and. v(preduc) .le. emaxs)
     1                       iv(irc) = 11
      if (v(dst0) .lt. zero) go to 250
      i = 0
      if ((v(nreduc) .gt. zero .and. v(nreduc) .le. emax) .or.
     1    (v(nreduc) .eq. zero. and. v(preduc) .eq. zero))  i = 2
      if (v(stppar) .eq. zero .and. v(reldx) .le. v(xctol)
     1                        .and. goodx)                  i = i + 1
      if (i .gt. 0) iv(irc) = i + 6
c
c  ***  consider recomputing step of length v(lmaxs) for singular
c  ***  convergence test.
c
 250  if (iv(irc) .gt. 5 .and. iv(irc) .ne. 12) go to 999
      if (v(dstnrm) .gt. v(lmaxs)) go to 260
         if (v(preduc) .ge. emaxs) go to 999
              if (v(dst0) .le. zero) go to 270
                   if (half * v(dst0) .le. v(lmaxs)) go to 999
                        go to 270
 260  if (half * v(dstnrm) .le. v(lmaxs)) go to 999
      xmax = v(lmaxs) / v(dstnrm)
      if (xmax * (two - xmax) * v(preduc) .ge. emaxs) go to 999
 270  if (v(nreduc) .lt. zero) go to 290
c
c  ***  recompute v(preduc) for use in singular convergence test  ***
c
      v(gtslst) = v(gtstep)
      v(dstsav) = v(dstnrm)
      if (iv(irc) .eq. 12) v(dstsav) = -v(dstsav)
      v(plstgd) = v(preduc)
      i = iv(restor)
      iv(restor) = 2
      if (i .eq. 3) iv(restor) = 0
      iv(irc) = 6
      go to 999
c
c  ***  perform singular convergence test with recomputed v(preduc)  ***
c
 280  v(gtstep) = v(gtslst)
      v(dstnrm) = dabs(v(dstsav))
      iv(irc) = iv(xirc)
      if (v(dstsav) .le. zero) iv(irc) = 12
      v(nreduc) = -v(preduc)
      v(preduc) = v(plstgd)
      iv(restor) = 3
 290  if (-v(nreduc) .le. v(rfctol) * dabs(v(f0))) iv(irc) = 11
c
 999  return
c
c  ***  last card of assst follows  ***
      end
      subroutine deflt(alg, iv, liv, lv, v)
c
c  ***  supply ***sol (version 2.3) default values to iv and v  ***
c
c  ***  alg = 1 means regression constants.
c  ***  alg = 2 means general unconstrained optimization constants.
c
      integer liv, lv
      integer alg, iv(liv)
      double precision v(lv)
c
      external imdcon, vdflt
      integer imdcon
c imdcon... returns machine-dependent integer constants.
c vdflt.... provides default values to v.
c
      integer miv, mv
      integer miniv(2), minv(2)
c
c  ***  subscripts for iv  ***
c
      integer algsav, covprt, covreq, dtype, hc, ierr, inith, inits,
     1        ipivot, ivneed, lastiv, lastv, lmat, mxfcal, mxiter,
     2        nfcov, ngcov, nvdflt, outlev, parprt, parsav, perm,
     3        prunit, qrtyp, rdreq, rmat, solprt, statpr, vneed,
     4        vsave, x0prt
c
c  ***  iv subscript values  ***
c
c/6
      data algsav/51/, covprt/14/, covreq/15/, dtype/16/, hc/71/,
     1     ierr/75/, inith/25/, inits/25/, ipivot/76/, ivneed/3/,
     2     lastiv/44/, lastv/45/, lmat/42/, mxfcal/17/, mxiter/18/,
     3     nfcov/52/, ngcov/53/, nvdflt/50/, outlev/19/, parprt/20/,
     4     parsav/49/, perm/58/, prunit/21/, qrtyp/80/, rdreq/57/,
     5     rmat/78/, solprt/22/, statpr/23/, vneed/4/, vsave/60/,
     6     x0prt/24/
c/7
c     parameter (algsav=51, covprt=14, covreq=15, dtype=16, hc=71,
c    1           ierr=75, inith=25, inits=25, ipivot=76, ivneed=3,
c    2           lastiv=44, lastv=45, lmat=42, mxfcal=17, mxiter=18,
c    3           nfcov=52, ngcov=53, nvdflt=50, outlev=19, parprt=20,
c    4           parsav=49, perm=58, prunit=21, qrtyp=80, rdreq=57,
c    5           rmat=78, solprt=22, statpr=23, vneed=4, vsave=60,
c    6           x0prt=24)
c/
      data miniv(1)/80/, miniv(2)/59/, minv(1)/98/, minv(2)/71/
c
c-------------------------------  body  --------------------------------
c
      if (alg .lt. 1 .or. alg .gt. 2) go to 40
      miv = miniv(alg)
      if (liv .lt. miv) go to 20
      mv = minv(alg)
      if (lv .lt. mv) go to 30
      call vdflt(alg, lv, v)
      iv(1) = 12
      iv(algsav) = alg
      iv(ivneed) = 0
      iv(lastiv) = miv
      iv(lastv) = mv
      iv(lmat) = mv + 1
      iv(mxfcal) = 200
      iv(mxiter) = 150
      iv(outlev) = 1
      iv(parprt) = 1
      iv(perm) = miv + 1
      iv(prunit) = imdcon(1)
      iv(solprt) = 1
      iv(statpr) = 1
      iv(vneed) = 0
      iv(x0prt) = 1
c
      if (alg .ge. 2) go to 10
c
c  ***  regression  values
c
      iv(covprt) = 3
      iv(covreq) = 1
      iv(dtype) = 1
      iv(hc) = 0
      iv(ierr) = 0
      iv(inits) = 0
      iv(ipivot) = 0
      iv(nvdflt) = 32
      iv(parsav) = 67
      iv(qrtyp) = 1
      iv(rdreq) = 3
      iv(rmat) = 0
      iv(vsave) = 58
      go to 999
c
c  ***  general optimization values
c
 10   iv(dtype) = 0
      iv(inith) = 1
      iv(nfcov) = 0
      iv(ngcov) = 0
      iv(nvdflt) = 25
      iv(parsav) = 47
      go to 999
c
 20   iv(1) = 15
      go to 999
c
 30   iv(1) = 16
      go to 999
c
 40   iv(1) = 67
c
 999  return
c  ***  last card of deflt follows  ***
      end
      double precision function dotprd(p, x, y)
c
c  ***  return the inner product of the p-vectors x and y.  ***
c
      integer p
      double precision x(p), y(p)
c
      integer i
      double precision one, sqteta, t, zero
c/+
      double precision dmax1, dabs
c/
      external rmdcon
      double precision rmdcon
c
c  ***  rmdcon(2) returns a machine-dependent constant, sqteta, which
c  ***  is slightly larger than the smallest positive number that
c  ***  can be squared without underflowing.
c
c/6
      data one/1.d+0/, sqteta/0.d+0/, zero/0.d+0/
c/7
c     parameter (one=1.d+0, zero=0.d+0)
c     data sqteta/0.d+0/
c/
c
      dotprd = zero
      if (p .le. 0) go to 999
      if (sqteta .eq. zero) sqteta = rmdcon(2)
      do 20 i = 1, p
         t = dmax1(dabs(x(i)), dabs(y(i)))
         if (t .gt. one) go to 10
         if (t .lt. sqteta) go to 20
         t = (x(i)/sqteta)*y(i)
         if (dabs(t) .lt. sqteta) go to 20
 10      dotprd = dotprd + x(i)*y(i)
 20   continue
c
 999  return
c  ***  last card of dotprd follows  ***
      end
      subroutine itsum(d, g, iv, liv, lv, p, v, x)
c
c  ***  print iteration summary for ***sol (version 2.3)  ***
c
c  ***  parameter declarations  ***
c
      integer liv, lv, p
      integer iv(liv)
      double precision d(p), g(p), v(lv), x(p)
c
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
c  ***  local variables  ***
c
      integer alg, i, iv1, m, nf, ng, ol, pu
c/6
      real model1(6), model2(6)
c/7
c     character*4 model1(6), model2(6)
c/
      double precision nreldf, oldf, preldf, reldf, zero
c
c  ***  intrinsic functions  ***
c/+
      integer iabs
      double precision dabs, dmax1
c/
c  ***  no external functions or subroutines  ***
c
c  ***  subscripts for iv and v  ***
c
      integer algsav, dstnrm, f, fdif, f0, needhd, nfcall, nfcov, ngcov,
     1        ngcall, niter, nreduc, outlev, preduc, prntit, prunit,
     2        reldx, solprt, statpr, stppar, sused, x0prt
c
c  ***  iv subscript values  ***
c
c/6
      data algsav/51/, needhd/36/, nfcall/6/, nfcov/52/, ngcall/30/,
     1     ngcov/53/, niter/31/, outlev/19/, prntit/39/, prunit/21/,
     2     solprt/22/, statpr/23/, sused/64/, x0prt/24/
c/7
c     parameter (algsav=51, needhd=36, nfcall=6, nfcov=52, ngcall=30,
c    1           ngcov=53, niter=31, outlev=19, prntit=39, prunit=21,
c    2           solprt=22, statpr=23, sused=64, x0prt=24)
c/
c
c  ***  v subscript values  ***
c
c/6
      data dstnrm/2/, f/10/, f0/13/, fdif/11/, nreduc/6/, preduc/7/,
     1     reldx/17/, stppar/5/
c/7
c     parameter (dstnrm=2, f=10, f0=13, fdif=11, nreduc=6, preduc=7,
c    1           reldx=17, stppar=5)
c/
c
c/6
      data zero/0.d+0/
c/7
c     parameter (zero=0.d+0)
c/
c/6
      data model1(1)/4h    /, model1(2)/4h    /, model1(3)/4h    /,
     1     model1(4)/4h    /, model1(5)/4h  g /, model1(6)/4h  s /,
     2     model2(1)/4h g  /, model2(2)/4h s  /, model2(3)/4hg-s /,
     3     model2(4)/4hs-g /, model2(5)/4h-s-g/, model2(6)/4h-g-s/
c/7
c     data model1/'    ','    ','    ','    ','  g ','  s '/,
c    1     model2/' g  ',' s  ','g-s ','s-g ','-s-g','-g-s'/
c/
c
c-------------------------------  body  --------------------------------
c
      pu = iv(prunit)
      if (pu .eq. 0) go to 999
      iv1 = iv(1)
      if (iv1 .gt. 62) iv1 = iv1 - 51
      ol = iv(outlev)
      alg = iv(algsav)
      if (iv1 .lt. 2 .or. iv1 .gt. 15) go to 370
      if (iv1 .ge. 12) go to 120
      if (iv1 .eq. 2 .and. iv(niter) .eq. 0) go to 390
      if (ol .eq. 0) go to 120
      if (iv1 .ge. 10 .and. iv(prntit) .eq. 0) go to 120
      if (iv1 .gt. 2) go to 10
         iv(prntit) = iv(prntit) + 1
         if (iv(prntit) .lt. iabs(ol)) go to 999
 10   nf = iv(nfcall) - iabs(iv(nfcov))
      iv(prntit) = 0
      reldf = zero
      preldf = zero
      oldf = dmax1(dabs(v(f0)), dabs(v(f)))
      if (oldf .le. zero) go to 20
         reldf = v(fdif) / oldf
         preldf = v(preduc) / oldf
 20   if (ol .gt. 0) go to 60
c
c        ***  print short summary line  ***
c
         if (iv(needhd) .eq. 1 .and. alg .eq. 1) write(pu,30)
 30   format(/10h   it   nf,6x,1hf,7x,5hreldf,3x,6hpreldf,3x,5hreldx,
     1       2x,13hmodel  stppar)
         if (iv(needhd) .eq. 1 .and. alg .eq. 2) write(pu,40)
 40   format(/11h    it   nf,7x,1hf,8x,5hreldf,4x,6hpreldf,4x,5hreldx,
     1       3x,6hstppar)
         iv(needhd) = 0
         if (alg .eq. 2) go to 50
         m = iv(sused)
         write(pu,100) iv(niter), nf, v(f), reldf, preldf, v(reldx),
     1                 model1(m), model2(m), v(stppar)
         go to 120
c
 50      write(pu,110) iv(niter), nf, v(f), reldf, preldf, v(reldx),
     1                 v(stppar)
         go to 120
c
c     ***  print long summary line  ***
c
 60   if (iv(needhd) .eq. 1 .and. alg .eq. 1) write(pu,70)
 70   format(/11h    it   nf,6x,1hf,7x,5hreldf,3x,6hpreldf,3x,5hreldx,
     1       2x,13hmodel  stppar,2x,6hd*step,2x,7hnpreldf)
      if (iv(needhd) .eq. 1 .and. alg .eq. 2) write(pu,80)
 80   format(/11h    it   nf,7x,1hf,8x,5hreldf,4x,6hpreldf,4x,5hreldx,
     1       3x,6hstppar,3x,6hd*step,3x,7hnpreldf)
      iv(needhd) = 0
      nreldf = zero
      if (oldf .gt. zero) nreldf = v(nreduc) / oldf
      if (alg .eq. 2) go to 90
      m = iv(sused)
      write(pu,100) iv(niter), nf, v(f), reldf, preldf, v(reldx),
     1             model1(m), model2(m), v(stppar), v(dstnrm), nreldf
      go to 120
c
 90   write(pu,110) iv(niter), nf, v(f), reldf, preldf,
     1             v(reldx), v(stppar), v(dstnrm), nreldf
 100  format(i6,i5,d10.3,2d9.2,d8.1,a3,a4,2d8.1,d9.2)
 110  format(i6,i5,d11.3,2d10.2,3d9.1,d10.2)
c
 120  if (iv(statpr) .lt. 0) go to 430
      go to (999, 999, 130, 150, 170, 190, 210, 230, 250, 270, 290, 310,
     1       330, 350, 520), iv1
c
 130  write(pu,140)
 140  format(/26h ***** x-convergence *****)
      go to 430
c
 150  write(pu,160)
 160  format(/42h ***** relative function convergence *****)
      go to 430
c
 170  write(pu,180)
 180  format(/49h ***** x- and relative function convergence *****)
      go to 430
c
 190  write(pu,200)
 200  format(/42h ***** absolute function convergence *****)
      go to 430
c
 210  write(pu,220)
 220  format(/33h ***** singular convergence *****)
      go to 430
c
 230  write(pu,240)
 240  format(/30h ***** false convergence *****)
      go to 430
c
 250  write(pu,260)
 260  format(/38h ***** function evaluation limit *****)
      go to 430
c
 270  write(pu,280)
 280  format(/28h ***** iteration limit *****)
      go to 430
c
 290  write(pu,300)
 300  format(/18h ***** stopx *****)
      go to 430
c
 310  write(pu,320)
 320  format(/44h ***** initial f(x) cannot be computed *****)
c
      go to 390
c
 330  write(pu,340)
 340  format(/37h ***** bad parameters to assess *****)
      go to 999
c
 350  write(pu,360)
 360  format(/43h ***** gradient could not be computed *****)
      if (iv(niter) .gt. 0) go to 480
      go to 390
c
 370  write(pu,380) iv(1)
 380  format(/14h ***** iv(1) =,i5,6h *****)
      go to 999
c
c  ***  initial call on itsum  ***
c
 390  if (iv(x0prt) .ne. 0) write(pu,400) (i, x(i), d(i), i = 1, p)
 400  format(/23h     i     initial x(i),8x,4hd(i)//(1x,i5,d17.6,d14.3))
c     *** the following are to avoid undefined variables when the
c     *** function evaluation limit is 1...
      v(dstnrm) = zero
      v(fdif) = zero
      v(nreduc) = zero
      v(preduc) = zero
      v(reldx)  = zero
      if (iv1 .ge. 12) go to 999
      iv(needhd) = 0
      iv(prntit) = 0
      if (ol .eq. 0) go to 999
      if (ol .lt. 0 .and. alg .eq. 1) write(pu,30)
      if (ol .lt. 0 .and. alg .eq. 2) write(pu,40)
      if (ol .gt. 0 .and. alg .eq. 1) write(pu,70)
      if (ol .gt. 0 .and. alg .eq. 2) write(pu,80)
      if (alg .eq. 1) write(pu,410) v(f)
      if (alg .eq. 2) write(pu,420) v(f)
 410  format(/11h     0    1,d10.3)
c365  format(/11h     0    1,e11.3)
 420  format(/11h     0    1,d11.3)
      go to 999
c
c  ***  print various information requested on solution  ***
c
 430  iv(needhd) = 1
      if (iv(statpr) .eq. 0) go to 480
         oldf = dmax1(dabs(v(f0)), dabs(v(f)))
         preldf = zero
         nreldf = zero
         if (oldf .le. zero) go to 440
              preldf = v(preduc) / oldf
              nreldf = v(nreduc) / oldf
 440     nf = iv(nfcall) - iv(nfcov)
         ng = iv(ngcall) - iv(ngcov)
         write(pu,450) v(f), v(reldx), nf, ng, preldf, nreldf
 450  format(/9h function,d17.6,8h   reldx,d17.3/12h func. evals,
     1   i8,9x,11hgrad. evals,i8/7h preldf,d16.3,6x,7hnpreldf,d15.3)
c
         if (iv(nfcov) .gt. 0) write(pu,460) iv(nfcov)
 460     format(/1x,i4,50h extra func. evals for covariance and diagnost
     1ics.)
         if (iv(ngcov) .gt. 0) write(pu,470) iv(ngcov)
 470     format(1x,i4,50h extra grad. evals for covariance and diagnosti
     1cs.)
c
 480  if (iv(solprt) .eq. 0) go to 999
         iv(needhd) = 1
         write(pu,490)
 490  format(/22h     i      final x(i),8x,4hd(i),10x,4hg(i)/)
         do 500 i = 1, p
              write(pu,510) i, x(i), d(i), g(i)
 500          continue
 510     format(1x,i5,d16.6,2d14.3)
      go to 999
c
 520  write(pu,530)
 530  format(/24h inconsistent dimensions)
 999  return
c  ***  last card of itsum follows  ***
      end
      subroutine litvmu(n, x, l, y)
c
c  ***  solve  (l**t)*x = y,  where  l  is an  n x n  lower triangular
c  ***  matrix stored compactly by rows.  x and y may occupy the same
c  ***  storage.  ***
c
      integer n
      double precision x(n), l(1), y(n)
      integer i, ii, ij, im1, i0, j, np1
      double precision xi, zero
c/6
      data zero/0.d+0/
c/7
c     parameter (zero=0.d+0)
c/
c
      do 10 i = 1, n
 10      x(i) = y(i)
      np1 = n + 1
      i0 = n*(n+1)/2
      do 30 ii = 1, n
         i = np1 - ii
         xi = x(i)/l(i0)
         x(i) = xi
         if (i .le. 1) go to 999
         i0 = i0 - i
         if (xi .eq. zero) go to 30
         im1 = i - 1
         do 20 j = 1, im1
              ij = i0 + j
              x(j) = x(j) - xi*l(ij)
 20           continue
 30      continue
 999  return
c  ***  last card of litvmu follows  ***
      end
      subroutine livmul(n, x, l, y)
c
c  ***  solve  l*x = y, where  l  is an  n x n  lower triangular
c  ***  matrix stored compactly by rows.  x and y may occupy the same
c  ***  storage.  ***
c
      integer n
      double precision x(n), l(1), y(n)
      external dotprd
      double precision dotprd
      integer i, j, k
      double precision t, zero
c/6
      data zero/0.d+0/
c/7
c     parameter (zero=0.d+0)
c/
c
      do 10 k = 1, n
         if (y(k) .ne. zero) go to 20
         x(k) = zero
 10      continue
      go to 999
 20   j = k*(k+1)/2
      x(k) = y(k) / l(j)
      if (k .ge. n) go to 999
      k = k + 1
      do 30 i = k, n
         t = dotprd(i-1, l(j+1), x)
         j = j + i
         x(i) = (y(i) - t)/l(j)
 30      continue
 999  return
c  ***  last card of livmul follows  ***
      end
      subroutine parck(alg, d, iv, liv, lv, n, v)
c
c  ***  check ***sol (version 2.3) parameters, print changed values  ***
c
c  ***  alg = 1 for regression, alg = 2 for general unconstrained opt.
c
      integer alg, liv, lv, n
      integer iv(liv)
      double precision d(n), v(lv)
c
      external rmdcon, vcopy, vdflt
      double precision rmdcon
c rmdcon -- returns machine-dependent constants.
c vcopy  -- copies one vector to another.
c vdflt  -- supplies default parameter values to v alone.
c/+
      integer max0
c/
c
c  ***  local variables  ***
c
      integer i, ii, iv1, j, k, l, m, miv1, miv2, ndfalt, parsv1, pu
      integer ijmp, jlim(2), miniv(2), ndflt(2)
c/6
      integer varnm(2), sh(2)
      real cngd(3), dflt(3), vn(2,34), which(3)
c/7
c     character*1 varnm(2), sh(2)
c     character*4 cngd(3), dflt(3), vn(2,34), which(3)
c/
      double precision big, machep, tiny, vk, vm(34), vx(34), zero
c
c  ***  iv and v subscripts  ***
c
      integer algsav, dinit, dtype, dtype0, epslon, inits, ivneed,
     1        lastiv, lastv, lmat, nextiv, nextv, nvdflt, oldn,
     2        parprt, parsav, perm, prunit, vneed
c
c
c/6
      data algsav/51/, dinit/38/, dtype/16/, dtype0/54/, epslon/19/,
     1     inits/25/, ivneed/3/, lastiv/44/, lastv/45/, lmat/42/,
     2     nextiv/46/, nextv/47/, nvdflt/50/, oldn/38/, parprt/20/,
     3     parsav/49/, perm/58/, prunit/21/, vneed/4/
c/7
c     parameter (algsav=51, dinit=38, dtype=16, dtype0=54, epslon=19,
c    1           inits=25, ivneed=3, lastiv=44, lastv=45, lmat=42,
c    2           nextiv=46, nextv=47, nvdflt=50, oldn=38, parprt=20,
c    3           parsav=49, perm=58, prunit=21, vneed=4)
c     save big, machep, tiny
c/
c
      data big/0.d+0/, machep/-1.d+0/, tiny/1.d+0/, zero/0.d+0/
c/6
      data vn(1,1),vn(2,1)/4hepsl,4hon../
      data vn(1,2),vn(2,2)/4hphmn,4hfc../
      data vn(1,3),vn(2,3)/4hphmx,4hfc../
      data vn(1,4),vn(2,4)/4hdecf,4hac../
      data vn(1,5),vn(2,5)/4hincf,4hac../
      data vn(1,6),vn(2,6)/4hrdfc,4hmn../
      data vn(1,7),vn(2,7)/4hrdfc,4hmx../
      data vn(1,8),vn(2,8)/4htune,4hr1../
      data vn(1,9),vn(2,9)/4htune,4hr2../
      data vn(1,10),vn(2,10)/4htune,4hr3../
      data vn(1,11),vn(2,11)/4htune,4hr4../
      data vn(1,12),vn(2,12)/4htune,4hr5../
      data vn(1,13),vn(2,13)/4hafct,4hol../
      data vn(1,14),vn(2,14)/4hrfct,4hol../
      data vn(1,15),vn(2,15)/4hxcto,4hl.../
      data vn(1,16),vn(2,16)/4hxfto,4hl.../
      data vn(1,17),vn(2,17)/4hlmax,4h0.../
      data vn(1,18),vn(2,18)/4hlmax,4hs.../
      data vn(1,19),vn(2,19)/4hscto,4hl.../
      data vn(1,20),vn(2,20)/4hdini,4ht.../
      data vn(1,21),vn(2,21)/4hdtin,4hit../
      data vn(1,22),vn(2,22)/4hd0in,4hit../
      data vn(1,23),vn(2,23)/4hdfac,4h..../
      data vn(1,24),vn(2,24)/4hdltf,4hdc../
      data vn(1,25),vn(2,25)/4hdltf,4hdj../
      data vn(1,26),vn(2,26)/4hdelt,4ha0../
      data vn(1,27),vn(2,27)/4hfuzz,4h..../
      data vn(1,28),vn(2,28)/4hrlim,4hit../
      data vn(1,29),vn(2,29)/4hcosm,4hin../
      data vn(1,30),vn(2,30)/4hhube,4hrc../
      data vn(1,31),vn(2,31)/4hrspt,4hol../
      data vn(1,32),vn(2,32)/4hsigm,4hin../
      data vn(1,33),vn(2,33)/4heta0,4h..../
      data vn(1,34),vn(2,34)/4hbias,4h..../
c/7
c     data vn(1,1),vn(2,1)/'epsl','on..'/
c     data vn(1,2),vn(2,2)/'phmn','fc..'/
c     data vn(1,3),vn(2,3)/'phmx','fc..'/
c     data vn(1,4),vn(2,4)/'decf','ac..'/
c     data vn(1,5),vn(2,5)/'incf','ac..'/
c     data vn(1,6),vn(2,6)/'rdfc','mn..'/
c     data vn(1,7),vn(2,7)/'rdfc','mx..'/
c     data vn(1,8),vn(2,8)/'tune','r1..'/
c     data vn(1,9),vn(2,9)/'tune','r2..'/
c     data vn(1,10),vn(2,10)/'tune','r3..'/
c     data vn(1,11),vn(2,11)/'tune','r4..'/
c     data vn(1,12),vn(2,12)/'tune','r5..'/
c     data vn(1,13),vn(2,13)/'afct','ol..'/
c     data vn(1,14),vn(2,14)/'rfct','ol..'/
c     data vn(1,15),vn(2,15)/'xcto','l...'/
c     data vn(1,16),vn(2,16)/'xfto','l...'/
c     data vn(1,17),vn(2,17)/'lmax','0...'/
c     data vn(1,18),vn(2,18)/'lmax','s...'/
c     data vn(1,19),vn(2,19)/'scto','l...'/
c     data vn(1,20),vn(2,20)/'dini','t...'/
c     data vn(1,21),vn(2,21)/'dtin','it..'/
c     data vn(1,22),vn(2,22)/'d0in','it..'/
c     data vn(1,23),vn(2,23)/'dfac','....'/
c     data vn(1,24),vn(2,24)/'dltf','dc..'/
c     data vn(1,25),vn(2,25)/'dltf','dj..'/
c     data vn(1,26),vn(2,26)/'delt','a0..'/
c     data vn(1,27),vn(2,27)/'fuzz','....'/
c     data vn(1,28),vn(2,28)/'rlim','it..'/
c     data vn(1,29),vn(2,29)/'cosm','in..'/
c     data vn(1,30),vn(2,30)/'hube','rc..'/
c     data vn(1,31),vn(2,31)/'rspt','ol..'/
c     data vn(1,32),vn(2,32)/'sigm','in..'/
c     data vn(1,33),vn(2,33)/'eta0','....'/
c     data vn(1,34),vn(2,34)/'bias','....'/
c/
c
      data vm(1)/1.0d-3/, vm(2)/-0.99d+0/, vm(3)/1.0d-3/, vm(4)/1.0d-2/,
     1     vm(5)/1.2d+0/, vm(6)/1.d-2/, vm(7)/1.2d+0/, vm(8)/0.d+0/,
     2     vm(9)/0.d+0/, vm(10)/1.d-3/, vm(11)/-1.d+0/, vm(13)/0.d+0/,
     3     vm(15)/0.d+0/, vm(16)/0.d+0/, vm(19)/0.d+0/, vm(20)/-10.d+0/,
     4     vm(21)/0.d+0/, vm(22)/0.d+0/, vm(23)/0.d+0/, vm(27)/1.01d+0/,
     5     vm(28)/1.d+10/, vm(30)/0.d+0/, vm(31)/0.d+0/, vm(32)/0.d+0/,
     6     vm(34)/0.d+0/
      data vx(1)/0.9d+0/, vx(2)/-1.d-3/, vx(3)/1.d+1/, vx(4)/0.8d+0/,
     1     vx(5)/1.d+2/, vx(6)/0.8d+0/, vx(7)/1.d+2/, vx(8)/0.5d+0/,
     2     vx(9)/0.5d+0/, vx(10)/1.d+0/, vx(11)/1.d+0/, vx(14)/0.1d+0/,
     3     vx(15)/1.d+0/, vx(16)/1.d+0/, vx(19)/1.d+0/, vx(23)/1.d+0/,
     4     vx(24)/1.d+0/, vx(25)/1.d+0/, vx(26)/1.d+0/, vx(27)/1.d+10/,
     5     vx(29)/1.d+0/, vx(31)/1.d+0/, vx(32)/1.d+0/, vx(33)/1.d+0/,
     6     vx(34)/1.d+0/
c
c/6
      data varnm(1)/1hp/, varnm(2)/1hn/, sh(1)/1hs/, sh(2)/1hh/
      data cngd(1),cngd(2),cngd(3)/4h---c,4hhang,4hed v/,
     1     dflt(1),dflt(2),dflt(3)/4hnond,4hefau,4hlt v/
c/7
c     data varnm(1)/'p'/, varnm(2)/'n'/, sh(1)/'s'/, sh(2)/'h'/
c     data cngd(1),cngd(2),cngd(3)/'---c','hang','ed v'/,
c    1     dflt(1),dflt(2),dflt(3)/'nond','efau','lt v'/
c/
      data ijmp/33/, jlim(1)/0/, jlim(2)/24/, ndflt(1)/32/, ndflt(2)/25/
      data miniv(1)/80/, miniv(2)/59/
c
c...............................  body  ................................
c
      pu = 0
      if (prunit .le. liv) pu = iv(prunit)
      if (alg .lt. 1 .or. alg .gt. 2) go to 340
      if (iv(1) .eq. 0) call deflt(alg, iv, liv, lv, v)
      iv1 = iv(1)
      if (iv1 .ne. 13 .and. iv1 .ne. 12) go to 10
      miv1 = miniv(alg)
      if (perm .le. liv) miv1 = max0(miv1, iv(perm) - 1)
      if (ivneed .le. liv) miv2 = miv1 + max0(iv(ivneed), 0)
      if (lastiv .le. liv) iv(lastiv) = miv2
      if (liv .lt. miv1) go to 300
      iv(ivneed) = 0
      iv(lastv) = max0(iv(vneed), 0) + iv(lmat) - 1
      iv(vneed) = 0
      if (liv .lt. miv2) go to 300
      if (lv .lt. iv(lastv)) go to 320
 10   if (alg .eq. iv(algsav)) go to 30
         if (pu .ne. 0) write(pu,20) alg, iv(algsav)
 20      format(/39h the first parameter to deflt should be,i3,
     1          12h rather than,i3)
         iv(1) = 82
         go to 999
 30   if (iv1 .lt. 12 .or. iv1 .gt. 14) go to 60
         if (n .ge. 1) go to 50
              iv(1) = 81
              if (pu .eq. 0) go to 999
              write(pu,40) varnm(alg), n
 40           format(/8h /// bad,a1,2h =,i5)
              go to 999
 50      if (iv1 .ne. 14) iv(nextiv) = iv(perm)
         if (iv1 .ne. 14) iv(nextv) = iv(lmat)
         if (iv1 .eq. 13) go to 999
         k = iv(parsav) - epslon
         call vdflt(alg, lv-k, v(k+1))
         iv(dtype0) = 2 - alg
         iv(oldn) = n
         which(1) = dflt(1)
         which(2) = dflt(2)
         which(3) = dflt(3)
         go to 110
 60   if (n .eq. iv(oldn)) go to 80
         iv(1) = 17
         if (pu .eq. 0) go to 999
         write(pu,70) varnm(alg), iv(oldn), n
 70      format(/5h /// ,1a1,14h changed from ,i5,4h to ,i5)
         go to 999
c
 80   if (iv1 .le. 11 .and. iv1 .ge. 1) go to 100
         iv(1) = 80
         if (pu .ne. 0) write(pu,90) iv1
 90      format(/13h ///  iv(1) =,i5,28h should be between 0 and 14.)
         go to 999
c
 100  which(1) = cngd(1)
      which(2) = cngd(2)
      which(3) = cngd(3)
c
 110  if (iv1 .eq. 14) iv1 = 12
      if (big .gt. tiny) go to 120
         tiny = rmdcon(1)
         machep = rmdcon(3)
         big = rmdcon(6)
         vm(12) = machep
         vx(12) = big
         vx(13) = big
         vm(14) = machep
         vm(17) = tiny
         vx(17) = big
         vm(18) = tiny
         vx(18) = big
         vx(20) = big
         vx(21) = big
         vx(22) = big
         vm(24) = machep
         vm(25) = machep
         vm(26) = machep
         vx(28) = rmdcon(5)
         vm(29) = machep
         vx(30) = big
         vm(33) = machep
 120  m = 0
      i = 1
      j = jlim(alg)
      k = epslon
      ndfalt = ndflt(alg)
      do 150 l = 1, ndfalt
         vk = v(k)
         if (vk .ge. vm(i) .and. vk .le. vx(i)) go to 140
              m = k
              if (pu .ne. 0) write(pu,130) vn(1,i), vn(2,i), k, vk,
     1                                    vm(i), vx(i)
 130          format(/6h ///  ,2a4,5h.. v(,i2,3h) =,d11.3,7h should,
     1               11h be between,d11.3,4h and,d11.3)
 140     k = k + 1
         i = i + 1
         if (i .eq. j) i = ijmp
 150     continue
c
      if (iv(nvdflt) .eq. ndfalt) go to 170
         iv(1) = 51
         if (pu .eq. 0) go to 999
         write(pu,160) iv(nvdflt), ndfalt
 160     format(/13h iv(nvdflt) =,i5,13h rather than ,i5)
         go to 999
 170  if ((iv(dtype) .gt. 0 .or. v(dinit) .gt. zero) .and. iv1 .eq. 12)
     1                  go to 200
      do 190 i = 1, n
         if (d(i) .gt. zero) go to 190
              m = 18
              if (pu .ne. 0) write(pu,180) i, d(i)
 180     format(/8h ///  d(,i3,3h) =,d11.3,19h should be positive)
 190     continue
 200  if (m .eq. 0) go to 210
         iv(1) = m
         go to 999
c
 210  if (pu .eq. 0 .or. iv(parprt) .eq. 0) go to 999
      if (iv1 .ne. 12 .or. iv(inits) .eq. alg-1) go to 230
         m = 1
         write(pu,220) sh(alg), iv(inits)
 220     format(/22h nondefault values..../5h init,a1,14h..... iv(25) =,
     1          i3)
 230  if (iv(dtype) .eq. iv(dtype0)) go to 250
         if (m .eq. 0) write(pu,260) which
         m = 1
         write(pu,240) iv(dtype)
 240     format(20h dtype..... iv(16) =,i3)
 250  i = 1
      j = jlim(alg)
      k = epslon
      l = iv(parsav)
      ndfalt = ndflt(alg)
      do 290 ii = 1, ndfalt
         if (v(k) .eq. v(l)) go to 280
              if (m .eq. 0) write(pu,260) which
 260          format(/1h ,3a4,9halues..../)
              m = 1
              write(pu,270) vn(1,i), vn(2,i), k, v(k)
 270          format(1x,2a4,5h.. v(,i2,3h) =,d15.7)
 280     k = k + 1
         l = l + 1
         i = i + 1
         if (i .eq. j) i = ijmp
 290     continue
c
      iv(dtype0) = iv(dtype)
      parsv1 = iv(parsav)
      call vcopy(iv(nvdflt), v(parsv1), v(epslon))
      go to 999
c
 300  iv(1) = 15
      if (pu .eq. 0) go to 999
      write(pu,310) liv, miv2
 310  format(/10h /// liv =,i5,17h must be at least,i5)
      if (liv .lt. miv1) go to 999
      if (lv .lt. iv(lastv)) go to 320
      go to 999
c
 320  iv(1) = 16
      if (pu .eq. 0) go to 999
      write(pu,330) lv, iv(lastv)
 330  format(/9h /// lv =,i5,17h must be at least,i5)
      go to 999
c
 340  iv(1) = 67
      if (pu .eq. 0) go to 999
      write(pu,350) alg
 350  format(/10h /// alg =,i5,15h must be 1 or 2)
c
 999  return
c  ***  last card of parck follows  ***
      end
      double precision function reldst(p, d, x, x0)
c
c  ***  compute and return relative difference between x and x0  ***
c  ***  nl2sol version 2.2  ***
c
      integer p
      double precision d(p), x(p), x0(p)
c/+
      double precision dabs
c/
      integer i
      double precision emax, t, xmax, zero
c/6
      data zero/0.d+0/
c/7
c     parameter (zero=0.d+0)
c/
c
      emax = zero
      xmax = zero
      do 10 i = 1, p
         t = dabs(d(i) * (x(i) - x0(i)))
         if (emax .lt. t) emax = t
         t = d(i) * (dabs(x(i)) + dabs(x0(i)))
         if (xmax .lt. t) xmax = t
 10      continue
      reldst = zero
      if (xmax .gt. zero) reldst = emax / xmax
 999  return
c  ***  last card of reldst follows  ***
      end
      logical function stopx(idummy)
c     *****parameters...
      integer idummy
c
c     ..................................................................
c
c     *****purpose...
c     this function may serve as the stopx (asynchronous interruption)
c     function for the nl2sol (nonlinear least-squares) package at
c     those installations which do not wish to implement a
c     dynamic stopx.
c
c     *****algorithm notes...
c     at installations where the nl2sol system is used
c     interactively, this dummy stopx should be replaced by a
c     function that returns .true. if and only if the interrupt
c     (break) key has been pressed since the last call on stopx.
c
c     ..................................................................
c
      stopx = .false.
      return
      end
      subroutine vaxpy(p, w, a, x, y)
c
c  ***  set w = a*x + y  --  w, x, y = p-vectors, a = scalar  ***
c
      integer p
      double precision a, w(p), x(p), y(p)
c
      integer i
c
      do 10 i = 1, p
 10      w(i) = a*x(i) + y(i)
      return
      end
      subroutine vcopy(p, y, x)
c
c  ***  set y = x, where x and y are p-vectors  ***
c
      integer p
      double precision x(p), y(p)
c
      integer i
c
      do 10 i = 1, p
 10      y(i) = x(i)
      return
      end
      subroutine vdflt(alg, lv, v)
c
c  ***  supply ***sol (version 2.3) default values to v  ***
c
c  ***  alg = 1 means regression constants.
c  ***  alg = 2 means general unconstrained optimization constants.
c
      integer alg, lv
      double precision v(lv)
c/+
      double precision dmax1
c/
      external rmdcon
      double precision rmdcon
c rmdcon... returns machine-dependent constants
c
      double precision machep, mepcrt, one, sqteps, three
c
c  ***  subscripts for v  ***
c
      integer afctol, bias, cosmin, decfac, delta0, dfac, dinit, dltfdc,
     1        dltfdj, dtinit, d0init, epslon, eta0, fuzz, huberc,
     2        incfac, lmax0, lmaxs, phmnfc, phmxfc, rdfcmn, rdfcmx,
     3        rfctol, rlimit, rsptol, sctol, sigmin, tuner1, tuner2,
     4        tuner3, tuner4, tuner5, xctol, xftol
c
c/6
      data one/1.d+0/, three/3.d+0/
c/7
c     parameter (one=1.d+0, three=3.d+0)
c/
c
c  ***  v subscript values  ***
c
c/6
      data afctol/31/, bias/43/, cosmin/47/, decfac/22/, delta0/44/,
     1     dfac/41/, dinit/38/, dltfdc/42/, dltfdj/43/, dtinit/39/,
     2     d0init/40/, epslon/19/, eta0/42/, fuzz/45/, huberc/48/,
     3     incfac/23/, lmax0/35/, lmaxs/36/, phmnfc/20/, phmxfc/21/,
     4     rdfcmn/24/, rdfcmx/25/, rfctol/32/, rlimit/46/, rsptol/49/,
     5     sctol/37/, sigmin/50/, tuner1/26/, tuner2/27/, tuner3/28/,
     6     tuner4/29/, tuner5/30/, xctol/33/, xftol/34/
c/7
c     parameter (afctol=31, bias=43, cosmin=47, decfac=22, delta0=44,
c    1           dfac=41, dinit=38, dltfdc=42, dltfdj=43, dtinit=39,
c    2           d0init=40, epslon=19, eta0=42, fuzz=45, huberc=48,
c    3           incfac=23, lmax0=35, lmaxs=36, phmnfc=20, phmxfc=21,
c    4           rdfcmn=24, rdfcmx=25, rfctol=32, rlimit=46, rsptol=49,
c    5           sctol=37, sigmin=50, tuner1=26, tuner2=27, tuner3=28,
c    6           tuner4=29, tuner5=30, xctol=33, xftol=34)
c/
c
c-------------------------------  body  --------------------------------
c
      machep = rmdcon(3)
      v(afctol) = 1.d-20
      if (machep .gt. 1.d-10) v(afctol) = machep**2
      v(decfac) = 0.5d+0
      sqteps = rmdcon(4)
      v(dfac) = 0.6d+0
      v(delta0) = sqteps
      v(dtinit) = 1.d-6
      mepcrt = machep ** (one/three)
      v(d0init) = 1.d+0
      v(epslon) = 0.1d+0
      v(incfac) = 2.d+0
      v(lmax0) = 1.d+0
      v(lmaxs) = 1.d+0
      v(phmnfc) = -0.1d+0
      v(phmxfc) = 0.1d+0
      v(rdfcmn) = 0.1d+0
      v(rdfcmx) = 4.d+0
      v(rfctol) = dmax1(1.d-10, mepcrt**2)
      v(sctol) = v(rfctol)
      v(tuner1) = 0.1d+0
      v(tuner2) = 1.d-4
      v(tuner3) = 0.75d+0
      v(tuner4) = 0.5d+0
      v(tuner5) = 0.75d+0
      v(xctol) = sqteps
      v(xftol) = 1.d+2 * machep
c
      if (alg .ge. 2) go to 10
c
c  ***  regression  values
c
      v(cosmin) = dmax1(1.d-6, 1.d+2 * machep)
      v(dinit) = 0.d+0
      v(dltfdc) = mepcrt
      v(dltfdj) = sqteps
      v(fuzz) = 1.5d+0
      v(huberc) = 0.7d+0
      v(rlimit) = rmdcon(5)
      v(rsptol) = 1.d-3
      v(sigmin) = 1.d-4
      go to 999
c
c  ***  general optimization values
c
 10   v(bias) = 0.8d+0
      v(dinit) = -1.0d+0
      v(eta0) = 1.0d+3 * machep
c
 999  return
c  ***  last card of vdflt follows  ***
      end
      subroutine vscopy(p, y, s)
c
c  ***  set p-vector y to scalar s  ***
c
      integer p
      double precision s, y(p)
c
      integer i
c
      do 10 i = 1, p
 10      y(i) = s
      return
      end
      double precision function v2norm(p, x)
c
c  ***  return the 2-norm of the p-vector x, taking  ***
c  ***  care to avoid the most likely underflows.    ***
c
      integer p
      double precision x(p)
c
      integer i, j
      double precision one, r, scale, sqteta, t, xi, zero
c/+
      double precision dabs, dsqrt
c/
      external rmdcon
      double precision rmdcon
c
c/6
      data one/1.d+0/, zero/0.d+0/
c/7
c     parameter (one=1.d+0, zero=0.d+0)
c     save sqteta
c/
      data sqteta/0.d+0/
c
      if (p .gt. 0) go to 10
         v2norm = zero
         go to 999
 10   do 20 i = 1, p
         if (x(i) .ne. zero) go to 30
 20      continue
      v2norm = zero
      go to 999
c
 30   scale = dabs(x(i))
      if (i .lt. p) go to 40
         v2norm = scale
         go to 999
 40   t = one
      if (sqteta .eq. zero) sqteta = rmdcon(2)
c
c     ***  sqteta is (slightly larger than) the square root of the
c     ***  smallest positive floating point number on the machine.
c     ***  the tests involving sqteta are done to prevent underflows.
c
      j = i + 1
      do 60 i = j, p
         xi = dabs(x(i))
         if (xi .gt. scale) go to 50
              r = xi / scale
              if (r .gt. sqteta) t = t + r*r
              go to 60
 50           r = scale / xi
              if (r .le. sqteta) r = zero
              t = one  +  t * r*r
              scale = xi
 60      continue
c
      v2norm = scale * dsqrt(t)
 999  return
c  ***  last card of v2norm follows  ***
      end
      subroutine humsl(n, d, x, calcf, calcgh, iv, liv, lv, v,
     1                  uiparm, urparm, ufparm)
c
c  ***  minimize general unconstrained objective function using   ***
c  ***  (analytic) gradient and hessian provided by the caller.   ***
c
      integer liv, lv, n
      integer iv(liv), uiparm(1)
      double precision d(n), x(n), v(lv), urparm(1)
c     dimension v(78 + n*(n+12)), uiparm(*), urparm(*)
      external calcf, calcgh, ufparm
c
c------------------------------  discussion  ---------------------------
c
c        this routine is like sumsl, except that the subroutine para-
c     meter calcg of sumsl (which computes the gradient of the objec-
c     tive function) is replaced by the subroutine parameter calcgh,
c     which computes both the gradient and (lower triangle of the)
c     hessian of the objective function.  the calling sequence is...
c             call calcgh(n, x, nf, g, h, uiparm, urparm, ufparm)
c     parameters n, x, nf, g, uiparm, urparm, and ufparm are the same
c     as for sumsl, while h is an array of length n*(n+1)/2 in which
c     calcgh must store the lower triangle of the hessian at x.  start-
c     ing at h(1), calcgh must store the hessian entries in the order
c     (1,1), (2,1), (2,2), (3,1), (3,2), (3,3), ...
c        the value printed (by itsum) in the column labelled stppar
c     is the levenberg-marquardt used in computing the current step.
c     zero means a full newton step.  if the special case described in
c     ref. 1 is detected, then stppar is negated.  the value printed
c     in the column labelled npreldf is zero if the current hessian
c     is not positive definite.
c        it sometimes proves worthwhile to let d be determined from the
c     diagonal of the hessian matrix by setting iv(dtype) = 1 and
c     v(dinit) = 0.  the following iv and v components are relevant...
c
c iv(dtol)..... iv(59) gives the starting subscript in v of the dtol
c             array used when d is updated.  (iv(dtol) can be
c             initialized by calling humsl with iv(1) = 13.)
c iv(dtype).... iv(16) tells how the scale vector d should be chosen.
c             iv(dtype) .le. 0 means that d should not be updated, and
c             iv(dtype) .ge. 1 means that d should be updated as
c             described below with v(dfac).  default = 0.
c v(dfac)..... v(41) and the dtol and d0 arrays (see v(dtinit) and
c             v(d0init)) are used in updating the scale vector d when
c             iv(dtype) .gt. 0.  (d is initialized according to
c             v(dinit), described in sumsl.)  let
c                  d1(i) = max(sqrt(abs(h(i,i))), v(dfac)*d(i)),
c             where h(i,i) is the i-th diagonal element of the current
c             hessian.  if iv(dtype) = 1, then d(i) is set to d1(i)
c             unless d1(i) .lt. dtol(i), in which case d(i) is set to
c                  max(d0(i), dtol(i)).
c             if iv(dtype) .ge. 2, then d is updated during the first
c             iteration as for iv(dtype) = 1 (after any initialization
c             due to v(dinit)) and is left unchanged thereafter.
c             default = 0.6.
c v(dtinit)... v(39), if positive, is the value to which all components
c             of the dtol array (see v(dfac)) are initialized.  if
c             v(dtinit) = 0, then it is assumed that the caller has
c             stored dtol in v starting at v(iv(dtol)).
c             default = 10**-6.
c v(d0init)... v(40), if positive, is the value to which all components
c             of the d0 vector (see v(dfac)) are initialized.  if
c             v(dfac) = 0, then it is assumed that the caller has
c             stored d0 in v starting at v(iv(dtol)+n).  default = 1.0.
c
c  ***  reference  ***
c
c 1. gay, d.m. (1981), computing optimal locally constrained steps,
c         siam j. sci. statist. comput. 2, pp. 186-197.
c.
c  ***  general  ***
c
c     coded by david m. gay (winter 1980).  revised sept. 1982.
c     this subroutine was written in connection with research supported
c     in part by the national science foundation under grants
c     mcs-7600324 and mcs-7906671.
c
c----------------------------  declarations  ---------------------------
c
      external deflt, humit
c
c deflt... provides default input values for iv and v.
c humit... reverse-communication routine that does humsl algorithm.
c
      integer g1, h1, iv1, lh, nf
      double precision f
c
c  ***  subscripts for iv   ***
c
      integer g, h, nextv, nfcall, nfgcal, toobig, vneed
c
c/6
      data nextv/47/, nfcall/6/, nfgcal/7/, g/28/, h/56/, toobig/2/,
     1     vneed/4/
c/7
c     parameter (nextv=47, nfcall=6, nfgcal=7, g=28, h=56, toobig=2,
c    1           vneed=4)
c/
c
c+++++++++++++++++++++++++++++++  body  ++++++++++++++++++++++++++++++++
c
      lh = n * (n + 1) / 2
      if (iv(1) .eq. 0) call deflt(2, iv, liv, lv, v)
      if (iv(1) .eq. 12 .or. iv(1) .eq. 13)
     1     iv(vneed) = iv(vneed) + n*(n+3)/2
      iv1 = iv(1)
      if (iv1 .eq. 14) go to 10
      if (iv1 .gt. 2 .and. iv1 .lt. 12) go to 10
      g1 = 1
      h1 = 1
      if (iv1 .eq. 12) iv(1) = 13
      go to 20
c
 10   g1 = iv(g)
      h1 = iv(h)
c
 20   call humit(d, f, v(g1), v(h1), iv, lh, liv, lv, n, v, x)
      if (iv(1) - 2) 30, 40, 50
c
 30   nf = iv(nfcall)
      call calcf(n, x, nf, f, uiparm, urparm, ufparm)
      if (nf .le. 0) iv(toobig) = 1
      go to 20
c
 40   call calcgh(n, x, iv(nfgcal), v(g1), v(h1), uiparm, urparm,
     1            ufparm)
      go to 20
c
 50   if (iv(1) .ne. 14) go to 999
c
c  ***  storage allocation
c
      iv(g) = iv(nextv)
      iv(h) = iv(g) + n
      iv(nextv) = iv(h) + n*(n+1)/2
      if (iv1 .ne. 13) go to 10
c
 999  return
c  ***  last card of humsl follows  ***
      end
      subroutine humit(d, fx, g, h, iv, lh, liv, lv, n, v, x)
c
c  ***  carry out humsl (unconstrained minimization) iterations, using
c  ***  hessian matrix provided by the caller.
c
c  ***  parameter declarations  ***
c
      integer lh, liv, lv, n
      integer iv(liv)
      double precision d(n), fx, g(n), h(lh), v(lv), x(n)
c
c--------------------------  parameter usage  --------------------------
c
c d.... scale vector.
c fx... function value.
c g.... gradient vector.
c h.... lower triangle of the hessian, stored rowwise.
c iv... integer value array.
c lh... length of h = p*(p+1)/2.
c liv.. length of iv (at least 60).
c lv... length of v (at least 78 + n*(n+21)/2).
c n.... number of variables (components in x and g).
c v.... floating-point value array.
c x.... parameter vector.
c
c  ***  discussion  ***
c
c        parameters iv, n, v, and x are the same as the corresponding
c     ones to humsl (which see), except that v can be shorter (since
c     the part of v that humsl uses for storing g and h is not needed).
c     moreover, compared with humsl, iv(1) may have the two additional
c     output values 1 and 2, which are explained below, as is the use
c     of iv(toobig) and iv(nfgcal).  the value iv(g), which is an
c     output value from humsl, is not referenced by humit or the
c     subroutines it calls.
c
c iv(1) = 1 means the caller should set fx to f(x), the function value
c             at x, and call humit again, having changed none of the
c             other parameters.  an exception occurs if f(x) cannot be
c             computed (e.g. if overflow would occur), which may happen
c             because of an oversized step.  in this case the caller
c             should set iv(toobig) = iv(2) to 1, which will cause
c             humit to ignore fx and try a smaller step.  the para-
c             meter nf that humsl passes to calcf (for possible use by
c             calcgh) is a copy of iv(nfcall) = iv(6).
c iv(1) = 2 means the caller should set g to g(x), the gradient of f at
c             x, and h to the lower triangle of h(x), the hessian of f
c             at x, and call humit again, having changed none of the
c             other parameters except perhaps the scale vector d.
c                  the parameter nf that humsl passes to calcg is
c             iv(nfgcal) = iv(7).  if g(x) and h(x) cannot be evaluated,
c             then the caller may set iv(nfgcal) to 0, in which case
c             humit will return with iv(1) = 65.
c                  note -- humit overwrites h with the lower triangle
c             of  diag(d)**-1 * h(x) * diag(d)**-1.
c.
c  ***  general  ***
c
c     coded by david m. gay (winter 1980).  revised sept. 1982.
c     this subroutine was written in connection with research supported
c     in part by the national science foundation under grants
c     mcs-7600324 and mcs-7906671.
c
c        (see sumsl and humsl for references.)
c
c+++++++++++++++++++++++++++  declarations  ++++++++++++++++++++++++++++
c
c  ***  local variables  ***
c
      integer dg1, dummy, i, j, k, l, lstgst, nn1o2, step1,
     1        temp1, w1, x01
      double precision t
c
c     ***  constants  ***
c
      double precision one, onep2, zero
c
c  ***  no intrinsic functions  ***
c
c  ***  external functions and subroutines  ***
c
      external assst, deflt, dotprd, dupdu, gqtst, itsum, parck,
     1         reldst, slvmul, stopx, vaxpy, vcopy, vscopy, v2norm
      logical stopx
      double precision dotprd, reldst, v2norm
c
c assst.... assesses candidate step.
c deflt.... provides default iv and v input values.
c dotprd... returns inner product of two vectors.
c dupdu.... updates scale vector d.
c gqtst.... computes optimally locally constrained step.
c itsum.... prints iteration summary and info on initial and final x.
c parck.... checks validity of input iv and v values.
c reldst... computes v(reldx) = relative step size.
c slvmul... multiplies symmetric matrix times vector, given the lower
c             triangle of the matrix.
c stopx.... returns .true. if the break key has been pressed.
c vaxpy.... computes scalar times one vector plus another.
c vcopy.... copies one vector to another.
c vscopy... sets all elements of a vector to a scalar.
c v2norm... returns the 2-norm of a vector.
c
c  ***  subscripts for iv and v  ***
c
      integer cnvcod, dg, dgnorm, dinit, dstnrm, dtinit, dtol,
     1        dtype, d0init, f, f0, fdif, gtstep, incfac, irc, kagqt,
     2        lmat, lmax0, lmaxs, mode, model, mxfcal, mxiter, nextv,
     3        nfcall, nfgcal, ngcall, niter, preduc, radfac, radinc,
     4        radius, rad0, reldx, restor, step, stglim, stlstg, stppar,
     5        toobig, tuner4, tuner5, vneed, w, xirc, x0
c
c  ***  iv subscript values  ***
c
c/6
      data cnvcod/55/, dg/37/, dtol/59/, dtype/16/, irc/29/, kagqt/33/,
     1     lmat/42/, mode/35/, model/5/, mxfcal/17/, mxiter/18/,
     2     nextv/47/, nfcall/6/, nfgcal/7/, ngcall/30/, niter/31/,
     3     radinc/8/, restor/9/, step/40/, stglim/11/, stlstg/41/,
     4     toobig/2/, vneed/4/, w/34/, xirc/13/, x0/43/
c/7
c     parameter (cnvcod=55, dg=37, dtol=59, dtype=16, irc=29, kagqt=33,
c    1           lmat=42, mode=35, model=5, mxfcal=17, mxiter=18,
c    2           nextv=47, nfcall=6, nfgcal=7, ngcall=30, niter=31,
c    3           radinc=8, restor=9, step=40, stglim=11, stlstg=41,
c    4           toobig=2, vneed=4, w=34, xirc=13, x0=43)
c/
c
c  ***  v subscript values  ***
c
c/6
      data dgnorm/1/, dinit/38/, dstnrm/2/, dtinit/39/, d0init/40/,
     1     f/10/, f0/13/, fdif/11/, gtstep/4/, incfac/23/, lmax0/35/,
     2     lmaxs/36/, preduc/7/, radfac/16/, radius/8/, rad0/9/,
     3     reldx/17/, stppar/5/, tuner4/29/, tuner5/30/
c/7
c     parameter (dgnorm=1, dinit=38, dstnrm=2, dtinit=39, d0init=40,
c    1           f=10, f0=13, fdif=11, gtstep=4, incfac=23, lmax0=35,
c    2           lmaxs=36, preduc=7, radfac=16, radius=8, rad0=9,
c    3           reldx=17, stppar=5, tuner4=29, tuner5=30)
c/
c
c/6
      data one/1.d+0/, onep2/1.2d+0/, zero/0.d+0/
c/7
c     parameter (one=1.d+0, onep2=1.2d+0, zero=0.d+0)
c/
c
c+++++++++++++++++++++++++++++++  body  ++++++++++++++++++++++++++++++++
c
      i = iv(1)
      if (i .eq. 1) go to 30
      if (i .eq. 2) go to 40
c
c  ***  check validity of iv and v input values  ***
c
      if (iv(1) .eq. 0) call deflt(2, iv, liv, lv, v)
      if (iv(1) .eq. 12 .or. iv(1) .eq. 13)
     1     iv(vneed) = iv(vneed) + n*(n+21)/2 + 7
      call parck(2, d, iv, liv, lv, n, v)
      i = iv(1) - 2
      if (i .gt. 12) go to 999
      nn1o2 = n * (n + 1) / 2
      if (lh .ge. nn1o2) go to (210,210,210,210,210,210,160,120,160,
     1                          10,10,20), i
         iv(1) = 66
         go to 350
c
c  ***  storage allocation  ***
c
 10   iv(dtol) = iv(lmat) + nn1o2
      iv(x0) = iv(dtol) + 2*n
      iv(step) = iv(x0) + n
      iv(stlstg) = iv(step) + n
      iv(dg) = iv(stlstg) + n
      iv(w) = iv(dg) + n
      iv(nextv) = iv(w) + 4*n + 7
      if (iv(1) .ne. 13) go to 20
         iv(1) = 14
         go to 999
c
c  ***  initialization  ***
c
 20   iv(niter) = 0
      iv(nfcall) = 1
      iv(ngcall) = 1
      iv(nfgcal) = 1
      iv(mode) = -1
      iv(model) = 1
      iv(stglim) = 1
      iv(toobig) = 0
      iv(cnvcod) = 0
      iv(radinc) = 0
      v(rad0) = zero
      v(stppar) = zero
      if (v(dinit) .ge. zero) call vscopy(n, d, v(dinit))
      k = iv(dtol)
      if (v(dtinit) .gt. zero) call vscopy(n, v(k), v(dtinit))
      k = k + n
      if (v(d0init) .gt. zero) call vscopy(n, v(k), v(d0init))
      iv(1) = 1
      go to 999
c
 30   v(f) = fx
      if (iv(mode) .ge. 0) go to 210
      iv(1) = 2
      if (iv(toobig) .eq. 0) go to 999
         iv(1) = 63
         go to 350
c
c  ***  make sure gradient could be computed  ***
c
 40   if (iv(nfgcal) .ne. 0) go to 50
         iv(1) = 65
         go to 350
c
c  ***  update the scale vector d  ***
c
 50   dg1 = iv(dg)
      if (iv(dtype) .le. 0) go to 70
      k = dg1
      j = 0
      do 60 i = 1, n
         j = j + i
         v(k) = h(j)
         k = k + 1
 60      continue
      call dupdu(d, v(dg1), iv, liv, lv, n, v)
c
c  ***  compute scaled gradient and its norm  ***
c
 70   dg1 = iv(dg)
      k = dg1
      do 80 i = 1, n
         v(k) = g(i) / d(i)
         k = k + 1
 80      continue
      v(dgnorm) = v2norm(n, v(dg1))
c
c  ***  compute scaled hessian  ***
c
      k = 1
      do 100 i = 1, n
         t = one / d(i)
         do 90 j = 1, i
              h(k) = t * h(k) / d(j)
              k = k + 1
 90           continue
 100     continue
c
      if (iv(cnvcod) .ne. 0) go to 340
      if (iv(mode) .eq. 0) go to 300
c
c  ***  allow first step to have scaled 2-norm at most v(lmax0)  ***
c
      v(radius) = v(lmax0)
c
      iv(mode) = 0
c
c
c-----------------------------  main loop  -----------------------------
c
c
c  ***  print iteration summary, check iteration limit  ***
c
 110  call itsum(d, g, iv, liv, lv, n, v, x)
 120  k = iv(niter)
      if (k .lt. iv(mxiter)) go to 130
         iv(1) = 10
         go to 350
c
 130  iv(niter) = k + 1
c
c  ***  initialize for start of next iteration  ***
c
      dg1 = iv(dg)
      x01 = iv(x0)
      v(f0) = v(f)
      iv(irc) = 4
      iv(kagqt) = -1
c
c     ***  copy x to x0  ***
c
      call vcopy(n, v(x01), x)
c
c  ***  update radius  ***
c
      if (k .eq. 0) go to 150
      step1 = iv(step)
      k = step1
      do 140 i = 1, n
         v(k) = d(i) * v(k)
         k = k + 1
 140     continue
      v(radius) = v(radfac) * v2norm(n, v(step1))
c
c  ***  check stopx and function evaluation limit  ***
c
 150  if (.not. stopx(dummy)) go to 170
         iv(1) = 11
         go to 180
c
c     ***  come here when restarting after func. eval. limit or stopx.
c
 160  if (v(f) .ge. v(f0)) go to 170
         v(radfac) = one
         k = iv(niter)
         go to 130
c
 170  if (iv(nfcall) .lt. iv(mxfcal)) go to 190
         iv(1) = 9
 180     if (v(f) .ge. v(f0)) go to 350
c
c        ***  in case of stopx or function evaluation limit with
c        ***  improved v(f), evaluate the gradient at x.
c
              iv(cnvcod) = iv(1)
              go to 290
c
c. . . . . . . . . . . . .  compute candidate step  . . . . . . . . . .
c
 190  step1 = iv(step)
      dg1 = iv(dg)
      l = iv(lmat)
      w1 = iv(w)
      call gqtst(d, v(dg1), h, iv(kagqt), v(l), n, v(step1), v, v(w1))
      if (iv(irc) .eq. 6) go to 210
c
c  ***  check whether evaluating f(x0 + step) looks worthwhile  ***
c
      if (v(dstnrm) .le. zero) go to 210
      if (iv(irc) .ne. 5) go to 200
      if (v(radfac) .le. one) go to 200
      if (v(preduc) .le. onep2 * v(fdif)) go to 210
c
c  ***  compute f(x0 + step)  ***
c
 200  x01 = iv(x0)
      step1 = iv(step)
      call vaxpy(n, x, one, v(step1), v(x01))
      iv(nfcall) = iv(nfcall) + 1
      iv(1) = 1
      iv(toobig) = 0
      go to 999
c
c. . . . . . . . . . . . .  assess candidate step  . . . . . . . . . . .
c
 210  x01 = iv(x0)
      v(reldx) = reldst(n, d, x, v(x01))
      call assst(iv, liv, lv, v)
      step1 = iv(step)
      lstgst = iv(stlstg)
      if (iv(restor) .eq. 1) call vcopy(n, x, v(x01))
      if (iv(restor) .eq. 2) call vcopy(n, v(lstgst), v(step1))
      if (iv(restor) .ne. 3) go to 220
         call vcopy(n, v(step1), v(lstgst))
         call vaxpy(n, x, one, v(step1), v(x01))
         v(reldx) = reldst(n, d, x, v(x01))
c
 220  k = iv(irc)
      go to (230,260,260,260,230,240,250,250,250,250,250,250,330,300), k
c
c     ***  recompute step with new radius  ***
c
 230     v(radius) = v(radfac) * v(dstnrm)
         go to 150
c
c  ***  compute step of length v(lmaxs) for singular convergence test.
c
 240  v(radius) = v(lmaxs)
      go to 190
c
c  ***  convergence or false convergence  ***
c
 250  iv(cnvcod) = k - 4
      if (v(f) .ge. v(f0)) go to 340
         if (iv(xirc) .eq. 14) go to 340
              iv(xirc) = 14
c
c. . . . . . . . . . . .  process acceptable step  . . . . . . . . . . .
c
 260  if (iv(irc) .ne. 3) go to 290
         temp1 = lstgst
c
c     ***  prepare for gradient tests  ***
c     ***  set  temp1 = hessian * step + g(x0)
c     ***             = diag(d) * (h * step + g(x0))
c
c        use x0 vector as temporary.
         k = x01
         do 270 i = 1, n
              v(k) = d(i) * v(step1)
              k = k + 1
              step1 = step1 + 1
 270          continue
         call slvmul(n, v(temp1), h, v(x01))
         do 280 i = 1, n
              v(temp1) = d(i) * v(temp1) + g(i)
              temp1 = temp1 + 1
 280          continue
c
c  ***  compute gradient and hessian  ***
c
 290  iv(ngcall) = iv(ngcall) + 1
      iv(1) = 2
      go to 999
c
 300  iv(1) = 2
      if (iv(irc) .ne. 3) go to 110
c
c  ***  set v(radfac) by gradient tests  ***
c
      temp1 = iv(stlstg)
      step1 = iv(step)
c
c     ***  set  temp1 = diag(d)**-1 * (hessian*step + (g(x0)-g(x)))  ***
c
      k = temp1
      do 310 i = 1, n
         v(k) = (v(k) - g(i)) / d(i)
         k = k + 1
 310     continue
c
c     ***  do gradient tests  ***
c
      if (v2norm(n, v(temp1)) .le. v(dgnorm) * v(tuner4)) go to 320
           if (dotprd(n, g, v(step1))
     1               .ge. v(gtstep) * v(tuner5))  go to 110
 320            v(radfac) = v(incfac)
                go to 110
c
c. . . . . . . . . . . . . .  misc. details  . . . . . . . . . . . . . .
c
c  ***  bad parameters to assess  ***
c
 330  iv(1) = 64
      go to 350
c
c  ***  print summary of final iteration and other requested items  ***
c
 340  iv(1) = iv(cnvcod)
      iv(cnvcod) = 0
 350  call itsum(d, g, iv, liv, lv, n, v, x)
c
 999  return
c
c  ***  last card of humit follows  ***
      end
      subroutine dupdu(d, hdiag, iv, liv, lv, n, v)
c
c  ***  update scale vector d for humsl  ***
c
c  ***  parameter declarations  ***
c
      integer liv, lv, n
      integer iv(liv)
      double precision d(n), hdiag(n), v(lv)
c
c  ***  local variables  ***
c
      integer dtoli, d0i, i
      double precision t, vdfac
c
c  ***  intrinsic functions  ***
c/+
      double precision dabs, dmax1, dsqrt
c/
c  ***  subscripts for iv and v  ***
c
      integer dfac, dtol, dtype, niter
c/6
      data dfac/41/, dtol/59/, dtype/16/, niter/31/
c/7
c     parameter (dfac=41, dtol=59, dtype=16, niter=31)
c/
c
c-------------------------------  body  --------------------------------
c
      i = iv(dtype)
      if (i .eq. 1) go to 10
         if (iv(niter) .gt. 0) go to 999
c
 10   dtoli = iv(dtol)
      d0i = dtoli + n
      vdfac = v(dfac)
      do 20 i = 1, n
         t = dmax1(dsqrt(dabs(hdiag(i))), vdfac*d(i))
         if (t .lt. v(dtoli)) t = dmax1(v(dtoli), v(d0i))
         d(i) = t
         dtoli = dtoli + 1
         d0i = d0i + 1
 20      continue
c
 999  return
c  ***  last card of dupdu follows  ***
      end
      subroutine gqtst(d, dig, dihdi, ka, l, p, step, v, w)
c
c  *** compute goldfeld-quandt-trotter step by more-hebden technique ***
c  ***  (nl2sol version 2.2), modified a la more and sorensen  ***
c
c  ***  parameter declarations  ***
c
      integer ka, p
      double precision d(p), dig(p), dihdi(1), l(1), v(21), step(p),
     1                 w(1)
c     dimension dihdi(p*(p+1)/2), l(p*(p+1)/2), w(4*p+7)
c
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
c  ***  purpose  ***
c
c        given the (compactly stored) lower triangle of a scaled
c     hessian (approximation) and a nonzero scaled gradient vector,
c     this subroutine computes a goldfeld-quandt-trotter step of
c     approximate length v(radius) by the more-hebden technique.  in
c     other words, step is computed to (approximately) minimize
c     psi(step) = (g**t)*step + 0.5*(step**t)*h*step  such that the
c     2-norm of d*step is at most (approximately) v(radius), where
c     g  is the gradient,  h  is the hessian, and  d  is a diagonal
c     scale matrix whose diagonal is stored in the parameter d.
c     (gqtst assumes  dig = d**-1 * g  and  dihdi = d**-1 * h * d**-1.)
c
c  ***  parameter description  ***
c
c     d (in)  = the scale vector, i.e. the diagonal of the scale
c              matrix  d  mentioned above under purpose.
c   dig (in)  = the scaled gradient vector, d**-1 * g.  if g = 0, then
c              step = 0  and  v(stppar) = 0  are returned.
c dihdi (in)  = lower triangle of the scaled hessian (approximation),
c              i.e., d**-1 * h * d**-1, stored compactly by rows., i.e.,
c              in the order (1,1), (2,1), (2,2), (3,1), (3,2), etc.
c    ka (i/o) = the number of hebden iterations (so far) taken to deter-
c              mine step.  ka .lt. 0 on input means this is the first
c              attempt to determine step (for the present dig and dihdi)
c              -- ka is initialized to 0 in this case.  output with
c              ka = 0  (or v(stppar) = 0)  means  step = -(h**-1)*g.
c     l (i/o) = workspace of length p*(p+1)/2 for cholesky factors.
c     p (in)  = number of parameters -- the hessian is a  p x p  matrix.
c  step (i/o) = the step computed.
c     v (i/o) contains various constants and variables described below.
c     w (i/o) = workspace of length 4*p + 6.
c
c  ***  entries in v  ***
c
c v(dgnorm) (i/o) = 2-norm of (d**-1)*g.
c v(dstnrm) (output) = 2-norm of d*step.
c v(dst0)   (i/o) = 2-norm of d*(h**-1)*g (for pos. def. h only), or
c             overestimate of smallest eigenvalue of (d**-1)*h*(d**-1).
c v(epslon) (in)  = max. rel. error allowed for psi(step).  for the
c             step returned, psi(step) will exceed its optimal value
c             by less than -v(epslon)*psi(step).  suggested value = 0.1.
c v(gtstep) (out) = inner product between g and step.
c v(nreduc) (out) = psi(-(h**-1)*g) = psi(newton step)  (for pos. def.
c             h only -- v(nreduc) is set to zero otherwise).
c v(phmnfc) (in)  = tol. (together with v(phmxfc)) for accepting step
c             (more*s sigma).  the error v(dstnrm) - v(radius) must lie
c             between v(phmnfc)*v(radius) and v(phmxfc)*v(radius).
c v(phmxfc) (in)  (see v(phmnfc).)
c             suggested values -- v(phmnfc) = -0.25, v(phmxfc) = 0.5.
c v(preduc) (out) = psi(step) = predicted obj. func. reduction for step.
c v(radius) (in)  = radius of current (scaled) trust region.
c v(rad0)   (i/o) = value of v(radius) from previous call.
c v(stppar) (i/o) is normally the marquardt parameter, i.e. the alpha
c             described below under algorithm notes.  if h + alpha*d**2
c             (see algorithm notes) is (nearly) singular, however,
c             then v(stppar) = -alpha.
c
c  ***  usage notes  ***
c
c     if it is desired to recompute step using a different value of
c     v(radius), then this routine may be restarted by calling it
c     with all parameters unchanged except v(radius).  (this explains
c     why step and w are listed as i/o).  on an initial call (one with
c     ka .lt. 0), step and w need not be initialized and only compo-
c     nents v(epslon), v(stppar), v(phmnfc), v(phmxfc), v(radius), and
c     v(rad0) of v must be initialized.
c
c  ***  algorithm notes  ***
c
c        the desired g-q-t step (ref. 2, 3, 4, 6) satisfies
c     (h + alpha*d**2)*step = -g  for some nonnegative alpha such that
c     h + alpha*d**2 is positive semidefinite.  alpha and step are
c     computed by a scheme analogous to the one described in ref. 5.
c     estimates of the smallest and largest eigenvalues of the hessian
c     are obtained from the gerschgorin circle theorem enhanced by a
c     simple form of the scaling described in ref. 7.  cases in which
c     h + alpha*d**2 is nearly (or exactly) singular are handled by
c     the technique discussed in ref. 2.  in these cases, a step of
c     (exact) length v(radius) is returned for which psi(step) exceeds
c     its optimal value by less than -v(epslon)*psi(step).  the test
c     suggested in ref. 6 for detecting the special case is performed
c     once two matrix factorizations have been done -- doing so sooner
c     seems to degrade the performance of optimization routines that
c     call this routine.
c
c  ***  functions and subroutines called  ***
c
c dotprd - returns inner product of two vectors.
c litvmu - applies inverse-transpose of compact lower triang. matrix.
c livmul - applies inverse of compact lower triang. matrix.
c lsqrt  - finds cholesky factor (of compactly stored lower triang.).
c lsvmin - returns approx. to min. sing. value of lower triang. matrix.
c rmdcon - returns machine-dependent constants.
c v2norm - returns 2-norm of a vector.
c
c  ***  references  ***
c
c 1.  dennis, j.e., gay, d.m., and welsch, r.e. (1981), an adaptive
c             nonlinear least-squares algorithm, acm trans. math.
c             software, vol. 7, no. 3.
c 2.  gay, d.m. (1981), computing optimal locally constrained steps,
c             siam j. sci. statist. computing, vol. 2, no. 2, pp.
c             186-197.
c 3.  goldfeld, s.m., quandt, r.e., and trotter, h.f. (1966),
c             maximization by quadratic hill-climbing, econometrica 34,
c             pp. 541-551.
c 4.  hebden, m.d. (1973), an algorithm for minimization using exact
c             second derivatives, report t.p. 515, theoretical physics
c             div., a.e.r.e. harwell, oxon., england.
c 5.  more, j.j. (1978), the levenberg-marquardt algorithm, implemen-
c             tation and theory, pp.105-116 of springer lecture notes
c             in mathematics no. 630, edited by g.a. watson, springer-
c             verlag, berlin and new york.
c 6.  more, j.j., and sorensen, d.c. (1981), computing a trust region
c             step, technical report anl-81-83, argonne national lab.
c 7.  varga, r.s. (1965), minimal gerschgorin sets, pacific j. math. 15,
c             pp. 719-729.
c
c  ***  general  ***
c
c     coded by david m. gay.
c     this subroutine was written in connection with research
c     supported by the national science foundation under grants
c     mcs-7600324, dcr75-10143, 76-14311dss, mcs76-11989, and
c     mcs-7906671.
c
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
c  ***  local variables  ***
c
      logical restrt
      integer dggdmx, diag, diag0, dstsav, emax, emin, i, im1, inc, irc,
     1        j, k, kalim, kamin, k1, lk0, phipin, q, q0, uk0, x
      double precision alphak, aki, akk, delta, dst, eps, gtsta, lk,
     1                 oldphi, phi, phimax, phimin, psifac, rad, radsq,
     2                 root, si, sk, sw, t, twopsi, t1, t2, uk, wi
c
c     ***  constants  ***
      double precision big, dgxfac, epsfac, four, half, kappa, negone,
     1                 one, p001, six, three, two, zero
c
c  ***  intrinsic functions  ***
c/+
      double precision dabs, dmax1, dmin1, dsqrt
c/
c  ***  external functions and subroutines  ***
c
      external dotprd, litvmu, livmul, lsqrt, lsvmin, rmdcon, v2norm
      double precision dotprd, lsvmin, rmdcon, v2norm
c
c  ***  subscripts for v  ***
c
      integer dgnorm, dstnrm, dst0, epslon, gtstep, stppar, nreduc,
     1        phmnfc, phmxfc, preduc, radius, rad0
c/6
      data dgnorm/1/, dstnrm/2/, dst0/3/, epslon/19/, gtstep/4/,
     1     nreduc/6/, phmnfc/20/, phmxfc/21/, preduc/7/, radius/8/,
     2     rad0/9/, stppar/5/
c/7
c     parameter (dgnorm=1, dstnrm=2, dst0=3, epslon=19, gtstep=4,
c    1           nreduc=6, phmnfc=20, phmxfc=21, preduc=7, radius=8,
c    2           rad0=9, stppar=5)
c/
c
c/6
      data epsfac/50.0d+0/, four/4.0d+0/, half/0.5d+0/,
     1     kappa/2.0d+0/, negone/-1.0d+0/, one/1.0d+0/, p001/1.0d-3/,
     2     six/6.0d+0/, three/3.0d+0/, two/2.0d+0/, zero/0.0d+0/
c/7
c     parameter (epsfac=50.0d+0, four=4.0d+0, half=0.5d+0,
c    1     kappa=2.0d+0, negone=-1.0d+0, one=1.0d+0, p001=1.0d-3,
c    2     six=6.0d+0, three=3.0d+0, two=2.0d+0, zero=0.0d+0)
c     save dgxfac
c/
      data big/0.d+0/, dgxfac/0.d+0/
c
c  ***  body  ***
c
c     ***  store largest abs. entry in (d**-1)*h*(d**-1) at w(dggdmx).
      dggdmx = p + 1
c     ***  store gerschgorin over- and underestimates of the largest
c     ***  and smallest eigenvalues of (d**-1)*h*(d**-1) at w(emax)
c     ***  and w(emin) respectively.
      emax = dggdmx + 1
      emin = emax + 1
c     ***  for use in recomputing step, the final values of lk, uk, dst,
c     ***  and the inverse derivative of more*s phi at 0 (for pos. def.
c     ***  h) are stored in w(lk0), w(uk0), w(dstsav), and w(phipin)
c     ***  respectively.
      lk0 = emin + 1
      phipin = lk0 + 1
      uk0 = phipin + 1
      dstsav = uk0 + 1
c     ***  store diag of (d**-1)*h*(d**-1) in w(diag),...,w(diag0+p).
      diag0 = dstsav
      diag = diag0 + 1
c     ***  store -d*step in w(q),...,w(q0+p).
      q0 = diag0 + p
      q = q0 + 1
c     ***  allocate storage for scratch vector x  ***
      x = q + p
      rad = v(radius)
      radsq = rad**2
c     ***  phitol = max. error allowed in dst = v(dstnrm) = 2-norm of
c     ***  d*step.
      phimax = v(phmxfc) * rad
      phimin = v(phmnfc) * rad
      psifac = two * v(epslon) / (three * (four * (v(phmnfc) + one) *
     1                       (kappa + one)  +  kappa  +  two) * rad**2)
c     ***  oldphi is used to detect limits of numerical accuracy.  if
c     ***  we recompute step and it does not change, then we accept it.
      oldphi = zero
      eps = v(epslon)
      irc = 0
      restrt = .false.
      kalim = ka + 50
c
c  ***  start or restart, depending on ka  ***
c
      if (ka .ge. 0) go to 290
c
c  ***  fresh start  ***
c
      k = 0
      uk = negone
      ka = 0
      kalim = 50
      v(dgnorm) = v2norm(p, dig)
      v(nreduc) = zero
      v(dst0) = zero
      kamin = 3
      if (v(dgnorm) .eq. zero) kamin = 0
c
c     ***  store diag(dihdi) in w(diag0+1),...,w(diag0+p)  ***
c
      j = 0
      do 10 i = 1, p
         j = j + i
         k1 = diag0 + i
         w(k1) = dihdi(j)
 10      continue
c
c     ***  determine w(dggdmx), the largest element of dihdi  ***
c
      t1 = zero
      j = p * (p + 1) / 2
      do 20 i = 1, j
         t = dabs(dihdi(i))
         if (t1 .lt. t) t1 = t
 20      continue
      w(dggdmx) = t1
c
c  ***  try alpha = 0  ***
c
 30   call lsqrt(1, p, l, dihdi, irc)
      if (irc .eq. 0) go to 50
c        ***  indef. h -- underestimate smallest eigenvalue, use this
c        ***  estimate to initialize lower bound lk on alpha.
         j = irc*(irc+1)/2
         t = l(j)
         l(j) = one
         do 40 i = 1, irc
 40           w(i) = zero
         w(irc) = one
         call litvmu(irc, w, l, w)
         t1 = v2norm(irc, w)
         lk = -t / t1 / t1
         v(dst0) = -lk
         if (restrt) go to 210
         go to 70
c
c     ***  positive definite h -- compute unmodified newton step.  ***
 50   lk = zero
      t = lsvmin(p, l, w(q), w(q))
      if (t .ge. one) go to 60
         if (big .le. zero) big = rmdcon(6)
         if (v(dgnorm) .ge. t*t*big) go to 70
 60   call livmul(p, w(q), l, dig)
      gtsta = dotprd(p, w(q), w(q))
      v(nreduc) = half * gtsta
      call litvmu(p, w(q), l, w(q))
      dst = v2norm(p, w(q))
      v(dst0) = dst
      phi = dst - rad
      if (phi .le. phimax) go to 260
      if (restrt) go to 210
c
c  ***  prepare to compute gerschgorin estimates of largest (and
c  ***  smallest) eigenvalues.  ***
c
 70   k = 0
      do 100 i = 1, p
         wi = zero
         if (i .eq. 1) go to 90
         im1 = i - 1
         do 80 j = 1, im1
              k = k + 1
              t = dabs(dihdi(k))
              wi = wi + t
              w(j) = w(j) + t
 80           continue
 90      w(i) = wi
         k = k + 1
 100     continue
c
c  ***  (under-)estimate smallest eigenvalue of (d**-1)*h*(d**-1)  ***
c
      k = 1
      t1 = w(diag) - w(1)
      if (p .le. 1) go to 120
      do 110 i = 2, p
         j = diag0 + i
         t = w(j) - w(i)
         if (t .ge. t1) go to 110
              t1 = t
              k = i
 110     continue
c
 120  sk = w(k)
      j = diag0 + k
      akk = w(j)
      k1 = k*(k-1)/2 + 1
      inc = 1
      t = zero
      do 150 i = 1, p
         if (i .eq. k) go to 130
         aki = dabs(dihdi(k1))
         si = w(i)
         j = diag0 + i
         t1 = half * (akk - w(j) + si - aki)
         t1 = t1 + dsqrt(t1*t1 + sk*aki)
         if (t .lt. t1) t = t1
         if (i .lt. k) go to 140
 130     inc = i
 140     k1 = k1 + inc
 150     continue
c
      w(emin) = akk - t
      uk = v(dgnorm)/rad - w(emin)
      if (v(dgnorm) .eq. zero) uk = uk + p001 + p001*uk
      if (uk .le. zero) uk = p001
c
c  ***  compute gerschgorin (over-)estimate of largest eigenvalue  ***
c
      k = 1
      t1 = w(diag) + w(1)
      if (p .le. 1) go to 170
      do 160 i = 2, p
         j = diag0 + i
         t = w(j) + w(i)
         if (t .le. t1) go to 160
              t1 = t
              k = i
 160     continue
c
 170  sk = w(k)
      j = diag0 + k
      akk = w(j)
      k1 = k*(k-1)/2 + 1
      inc = 1
      t = zero
      do 200 i = 1, p
         if (i .eq. k) go to 180
         aki = dabs(dihdi(k1))
         si = w(i)
         j = diag0 + i
         t1 = half * (w(j) + si - aki - akk)
         t1 = t1 + dsqrt(t1*t1 + sk*aki)
         if (t .lt. t1) t = t1
         if (i .lt. k) go to 190
 180     inc = i
 190     k1 = k1 + inc
 200     continue
c
      w(emax) = akk + t
      lk = dmax1(lk, v(dgnorm)/rad - w(emax))
c
c     ***  alphak = current value of alpha (see alg. notes above).  we
c     ***  use more*s scheme for initializing it.
      alphak = dabs(v(stppar)) * v(rad0)/rad
c
      if (irc .ne. 0) go to 210
c
c  ***  compute l0 for positive definite h  ***
c
      call livmul(p, w, l, w(q))
      t = v2norm(p, w)
      w(phipin) = dst / t / t
      lk = dmax1(lk, phi*w(phipin))
c
c  ***  safeguard alphak and add alphak*i to (d**-1)*h*(d**-1)  ***
c
 210  ka = ka + 1
      if (-v(dst0) .ge. alphak .or. alphak .lt. lk .or. alphak .ge. uk)
     1                      alphak = uk * dmax1(p001, dsqrt(lk/uk))
      if (alphak .le. zero) alphak = half * uk
      if (alphak .le. zero) alphak = uk
      k = 0
      do 220 i = 1, p
         k = k + i
         j = diag0 + i
         dihdi(k) = w(j) + alphak
 220     continue
c
c  ***  try computing cholesky decomposition  ***
c
      call lsqrt(1, p, l, dihdi, irc)
      if (irc .eq. 0) go to 240
c
c  ***  (d**-1)*h*(d**-1) + alphak*i  is indefinite -- overestimate
c  ***  smallest eigenvalue for use in updating lk  ***
c
      j = (irc*(irc+1))/2
      t = l(j)
      l(j) = one
      do 230 i = 1, irc
 230     w(i) = zero
      w(irc) = one
      call litvmu(irc, w, l, w)
      t1 = v2norm(irc, w)
      lk = alphak - t/t1/t1
      v(dst0) = -lk
      go to 210
c
c  ***  alphak makes (d**-1)*h*(d**-1) positive definite.
c  ***  compute q = -d*step, check for convergence.  ***
c
 240  call livmul(p, w(q), l, dig)
      gtsta = dotprd(p, w(q), w(q))
      call litvmu(p, w(q), l, w(q))
      dst = v2norm(p, w(q))
      phi = dst - rad
      if (phi .le. phimax .and. phi .ge. phimin) go to 270
      if (phi .eq. oldphi) go to 270
      oldphi = phi
      if (phi .lt. zero) go to 330
c
c  ***  unacceptable alphak -- update lk, uk, alphak  ***
c
 250  if (ka .ge. kalim) go to 270
c     ***  the following dmin1 is necessary because of restarts  ***
      if (phi .lt. zero) uk = dmin1(uk, alphak)
c     *** kamin = 0 only iff the gradient vanishes  ***
      if (kamin .eq. 0) go to 210
      call livmul(p, w, l, w(q))
      t1 = v2norm(p, w)
      alphak = alphak  +  (phi/t1) * (dst/t1) * (dst/rad)
      lk = dmax1(lk, alphak)
      go to 210
c
c  ***  acceptable step on first try  ***
c
 260  alphak = zero
c
c  ***  successful step in general.  compute step = -(d**-1)*q  ***
c
 270  do 280 i = 1, p
         j = q0 + i
         step(i) = -w(j)/d(i)
 280     continue
      v(gtstep) = -gtsta
      v(preduc) = half * (dabs(alphak)*dst*dst + gtsta)
      go to 410
c
c
c  ***  restart with new radius  ***
c
 290  if (v(dst0) .le. zero .or. v(dst0) - rad .gt. phimax) go to 310
c
c     ***  prepare to return newton step  ***
c
         restrt = .true.
         ka = ka + 1
         k = 0
         do 300 i = 1, p
              k = k + i
              j = diag0 + i
              dihdi(k) = w(j)
 300          continue
         uk = negone
         go to 30
c
 310  kamin = ka + 3
      if (v(dgnorm) .eq. zero) kamin = 0
      if (ka .eq. 0) go to 50
c
      dst = w(dstsav)
      alphak = dabs(v(stppar))
      phi = dst - rad
      t = v(dgnorm)/rad
      uk = t - w(emin)
      if (v(dgnorm) .eq. zero) uk = uk + p001 + p001*uk
      if (uk .le. zero) uk = p001
      if (rad .gt. v(rad0)) go to 320
c
c        ***  smaller radius  ***
         lk = zero
         if (alphak .gt. zero) lk = w(lk0)
         lk = dmax1(lk, t - w(emax))
         if (v(dst0) .gt. zero) lk = dmax1(lk, (v(dst0)-rad)*w(phipin))
         go to 250
c
c     ***  bigger radius  ***
 320  if (alphak .gt. zero) uk = dmin1(uk, w(uk0))
      lk = dmax1(zero, -v(dst0), t - w(emax))
      if (v(dst0) .gt. zero) lk = dmax1(lk, (v(dst0)-rad)*w(phipin))
      go to 250
c
c  ***  decide whether to check for special case... in practice (from
c  ***  the standpoint of the calling optimization code) it seems best
c  ***  not to check until a few iterations have failed -- hence the
c  ***  test on kamin below.
c
 330  delta = alphak + dmin1(zero, v(dst0))
      twopsi = alphak*dst*dst + gtsta
      if (ka .ge. kamin) go to 340
c     *** if the test in ref. 2 is satisfied, fall through to handle
c     *** the special case (as soon as the more-sorensen test detects
c     *** it).
      if (delta .ge. psifac*twopsi) go to 370
c
c  ***  check for the special case of  h + alpha*d**2  (nearly)
c  ***  singular.  use one step of inverse power method with start
c  ***  from lsvmin to obtain approximate eigenvector corresponding
c  ***  to smallest eigenvalue of (d**-1)*h*(d**-1).  lsvmin returns
c  ***  x and w with  l*w = x.
c
 340  t = lsvmin(p, l, w(x), w)
c
c     ***  normalize w  ***
      do 350 i = 1, p
 350     w(i) = t*w(i)
c     ***  complete current inv. power iter. -- replace w by (l**-t)*w.
      call litvmu(p, w, l, w)
      t2 = one/v2norm(p, w)
      do 360 i = 1, p
 360     w(i) = t2*w(i)
      t = t2 * t
c
c  ***  now w is the desired approximate (unit) eigenvector and
c  ***  t*x = ((d**-1)*h*(d**-1) + alphak*i)*w.
c
      sw = dotprd(p, w(q), w)
      t1 = (rad + dst) * (rad - dst)
      root = dsqrt(sw*sw + t1)
      if (sw .lt. zero) root = -root
      si = t1 / (sw + root)
c
c  ***  the actual test for the special case...
c
      if ((t2*si)**2 .le. eps*(dst**2 + alphak*radsq)) go to 380
c
c  ***  update upper bound on smallest eigenvalue (when not positive)
c  ***  (as recommended by more and sorensen) and continue...
c
      if (v(dst0) .le. zero) v(dst0) = dmin1(v(dst0), t2**2 - alphak)
      lk = dmax1(lk, -v(dst0))
c
c  ***  check whether we can hope to detect the special case in
c  ***  the available arithmetic.  accept step as it is if not.
c
c     ***  if not yet available, obtain machine dependent value dgxfac.
 370  if (dgxfac .eq. zero) dgxfac = epsfac * rmdcon(3)
c
      if (delta .gt. dgxfac*w(dggdmx)) go to 250
         go to 270
c
c  ***  special case detected... negate alphak to indicate special case
c
 380  alphak = -alphak
      v(preduc) = half * twopsi
c
c  ***  accept current step if adding si*w would lead to a
c  ***  further relative reduction in psi of less than v(epslon)/3.
c
      t1 = zero
      t = si*(alphak*sw - half*si*(alphak + t*dotprd(p,w(x),w)))
      if (t .lt. eps*twopsi/six) go to 390
         v(preduc) = v(preduc) + t
         dst = rad
         t1 = -si
 390  do 400 i = 1, p
         j = q0 + i
         w(j) = t1*w(i) - w(j)
         step(i) = w(j) / d(i)
 400     continue
      v(gtstep) = dotprd(p, dig, w(q))
c
c  ***  save values for use in a possible restart  ***
c
 410  v(dstnrm) = dst
      v(stppar) = alphak
      w(lk0) = lk
      w(uk0) = uk
      v(rad0) = rad
      w(dstsav) = dst
c
c     ***  restore diagonal of dihdi  ***
c
      j = 0
      do 420 i = 1, p
         j = j + i
         k = diag0 + i
         dihdi(j) = w(k)
 420     continue
c
 999  return
c
c  ***  last card of gqtst follows  ***
      end
      subroutine lsqrt(n1, n, l, a, irc)
c
c  ***  compute rows n1 through n of the cholesky factor  l  of
c  ***  a = l*(l**t),  where  l  and the lower triangle of  a  are both
c  ***  stored compactly by rows (and may occupy the same storage).
c  ***  irc = 0 means all went well.  irc = j means the leading
c  ***  principal  j x j  submatrix of  a  is not positive definite --
c  ***  and  l(j*(j+1)/2)  contains the (nonpos.) reduced j-th diagonal.
c
c  ***  parameters  ***
c
      integer n1, n, irc
      double precision l(1), a(1)
c     dimension l(n*(n+1)/2), a(n*(n+1)/2)
c
c  ***  local variables  ***
c
      integer i, ij, ik, im1, i0, j, jk, jm1, j0, k
      double precision t, td, zero
c
c  ***  intrinsic functions  ***
c/+
      double precision dsqrt
c/
c/6
      data zero/0.d+0/
c/7
c     parameter (zero=0.d+0)
c/
c
c  ***  body  ***
c
      i0 = n1 * (n1 - 1) / 2
      do 50 i = n1, n
         td = zero
         if (i .eq. 1) go to 40
         j0 = 0
         im1 = i - 1
         do 30 j = 1, im1
              t = zero
              if (j .eq. 1) go to 20
              jm1 = j - 1
              do 10 k = 1, jm1
                   ik = i0 + k
                   jk = j0 + k
                   t = t + l(ik)*l(jk)
 10                continue
 20           ij = i0 + j
              j0 = j0 + j
              t = (a(ij) - t) / l(j0)
              l(ij) = t
              td = td + t*t
 30           continue
 40      i0 = i0 + i
         t = a(i0) - td
         if (t .le. zero) go to 60
         l(i0) = dsqrt(t)
 50      continue
c
      irc = 0
      go to 999
c
 60   l(i0) = t
      irc = i
c
 999  return
c
c  ***  last card of lsqrt  ***
      end
      double precision function lsvmin(p, l, x, y)
c
c  ***  estimate smallest sing. value of packed lower triang. matrix l
c
c  ***  parameter declarations  ***
c
      integer p
      double precision l(1), x(p), y(p)
c     dimension l(p*(p+1)/2)
c
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
c  ***  purpose  ***
c
c     this function returns a good over-estimate of the smallest
c     singular value of the packed lower triangular matrix l.
c
c  ***  parameter description  ***
c
c  p (in)  = the order of l.  l is a  p x p  lower triangular matrix.
c  l (in)  = array holding the elements of  l  in row order, i.e.
c             l(1,1), l(2,1), l(2,2), l(3,1), l(3,2), l(3,3), etc.
c  x (out) if lsvmin returns a positive value, then x is a normalized
c             approximate left singular vector corresponding to the
c             smallest singular value.  this approximation may be very
c             crude.  if lsvmin returns zero, then some components of x
c             are zero and the rest retain their input values.
c  y (out) if lsvmin returns a positive value, then y = (l**-1)*x is an
c             unnormalized approximate right singular vector correspond-
c             ing to the smallest singular value.  this approximation
c             may be crude.  if lsvmin returns zero, then y retains its
c             input value.  the caller may pass the same vector for x
c             and y (nonstandard fortran usage), in which case y over-
c             writes x (for nonzero lsvmin returns).
c
c  ***  algorithm notes  ***
c
c     the algorithm is based on (1), with the additional provision that
c     lsvmin = 0 is returned if the smallest diagonal element of l
c     (in magnitude) is not more than the unit roundoff times the
c     largest.  the algorithm uses a random number generator proposed
c     in (4), which passes the spectral test with flying colors -- see
c     (2) and (3).
c
c  ***  subroutines and functions called  ***
c
c        v2norm - function, returns the 2-norm of a vector.
c
c  ***  references  ***
c
c     (1) cline, a., moler, c., stewart, g., and wilkinson, j.h.(1977),
c         an estimate for the condition number of a matrix, report
c         tm-310, applied math. div., argonne national laboratory.
c
c     (2) hoaglin, d.c. (1976), theoretical properties of congruential
c         random-number generators --  an empirical view,
c         memorandum ns-340, dept. of statistics, harvard univ.
c
c     (3) knuth, d.e. (1969), the art of computer programming, vol. 2
c         (seminumerical algorithms), addison-wesley, reading, mass.
c
c     (4) smith, c.s. (1971), multiplicative pseudo-random number
c         generators with prime modulus, j. assoc. comput. mach. 18,
c         pp. 586-593.
c
c  ***  history  ***
c
c     designed and coded by david m. gay (winter 1977/summer 1978).
c
c  ***  general  ***
c
c     this subroutine was written in connection with research
c     supported by the national science foundation under grants
c     mcs-7600324, dcr75-10143, 76-14311dss, and mcs76-11989.
c
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
c  ***  local variables  ***
c
      integer i, ii, ix, j, ji, jj, jjj, jm1, j0, pm1
      double precision b, sminus, splus, t, xminus, xplus
c
c  ***  constants  ***
c
      double precision half, one, r9973, zero
c
c  ***  intrinsic functions  ***
c/+
      integer mod
      real float
      double precision dabs
c/
c  ***  external functions and subroutines  ***
c
      external dotprd, v2norm, vaxpy
      double precision dotprd, v2norm
c
c/6
      data half/0.5d+0/, one/1.d+0/, r9973/9973.d+0/, zero/0.d+0/
c/7
c     parameter (half=0.5d+0, one=1.d+0, r9973=9973.d+0, zero=0.d+0)
c/
c
c  ***  body  ***
c
      ix = 2
      pm1 = p - 1
c
c  ***  first check whether to return lsvmin = 0 and initialize x  ***
c
      ii = 0
      j0 = p*pm1/2
      jj = j0 + p
      if (l(jj) .eq. zero) go to 110
      ix = mod(3432*ix, 9973)
      b = half*(one + float(ix)/r9973)
      xplus = b / l(jj)
      x(p) = xplus
      if (p .le. 1) go to 60
      do 10 i = 1, pm1
         ii = ii + i
         if (l(ii) .eq. zero) go to 110
         ji = j0 + i
         x(i) = xplus * l(ji)
 10      continue
c
c  ***  solve (l**t)*x = b, where the components of b have randomly
c  ***  chosen magnitudes in (.5,1) with signs chosen to make x large.
c
c     do j = p-1 to 1 by -1...
      do 50 jjj = 1, pm1
         j = p - jjj
c       ***  determine x(j) in this iteration. note for i = 1,2,...,j
c       ***  that x(i) holds the current partial sum for row i.
         ix = mod(3432*ix, 9973)
         b = half*(one + float(ix)/r9973)
         xplus = (b - x(j))
         xminus = (-b - x(j))
         splus = dabs(xplus)
         sminus = dabs(xminus)
         jm1 = j - 1
         j0 = j*jm1/2
         jj = j0 + j
         xplus = xplus/l(jj)
         xminus = xminus/l(jj)
         if (jm1 .eq. 0) go to 30
         do 20 i = 1, jm1
              ji = j0 + i
              splus = splus + dabs(x(i) + l(ji)*xplus)
              sminus = sminus + dabs(x(i) + l(ji)*xminus)
 20           continue
 30      if (sminus .gt. splus) xplus = xminus
         x(j) = xplus
c       ***  update partial sums  ***
         if (jm1 .gt. 0) call vaxpy(jm1, x, xplus, l(j0+1), x)
 50      continue
c
c  ***  normalize x  ***
c
 60   t = one/v2norm(p, x)
      do 70 i = 1, p
 70      x(i) = t*x(i)
c
c  ***  solve l*y = x and return lsvmin = 1/twonorm(y)  ***
c
      do 100 j = 1, p
         jm1 = j - 1
         j0 = j*jm1/2
         jj = j0 + j
         t = zero
         if (jm1 .gt. 0) t = dotprd(jm1, l(j0+1), y)
         y(j) = (x(j) - t) / l(jj)
 100     continue
c
      lsvmin = one/v2norm(p, y)
      go to 999
c
 110  lsvmin = zero
 999  return
c  ***  last card of lsvmin follows  ***
      end
      subroutine slvmul(p, y, s, x)
c
c  ***  set  y = s * x,  s = p x p symmetric matrix.  ***
c  ***  lower triangle of  s  stored rowwise.         ***
c
c  ***  parameter declarations  ***
c
      integer p
      double precision s(1), x(p), y(p)
c     dimension s(p*(p+1)/2)
c
c  ***  local variables  ***
c
      integer i, im1, j, k
      double precision xi
c
c  ***  no intrinsic functions  ***
c
c  ***  external function  ***
c
      external dotprd
      double precision dotprd
c
c-----------------------------------------------------------------------
c
      j = 1
      do 10 i = 1, p
         y(i) = dotprd(i, s(j), x)
         j = j + i
 10      continue
c
      if (p .le. 1) go to 999
      j = 1
      do 40 i = 2, p
         xi = x(i)
         im1 = i - 1
         j = j + 1
         do 30 k = 1, im1
              y(k) = y(k) + s(j)*xi
              j = j + 1
 30           continue
 40      continue
c
 999  return
c  ***  last card of slvmul follows  ***
      end
      function d1mach ( i )

c*********************************************************************72
c
cc D1MACH returns double precision real machine-dependent constants.
c
c  Discussion:
c
c    D1MACH can be used to obtain machine-dependent parameters
c    for the local machine environment.  It is a function
c    with one input argument, and can be called as follows:
c
c      D = D1MACH ( I )
c
c    where I=1,...,5.  The output value of D above is
c    determined by the input value of I:.
c
c    D1MACH ( 1) = B**(EMIN-1), the smallest positive magnitude.
c    D1MACH ( 2) = B**EMAX*(1 - B**(-T)), the largest magnitude.
c    D1MACH ( 3) = B**(-T), the smallest relative spacing.
c    D1MACH ( 4) = B**(1-T), the largest relative spacing.
c    D1MACH ( 5) = LOG10(B)
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    25 April 2007
c
c  Author:
c
c    Original FORTRAN77 version by Phyllis Fox, Andrew Hall, Norman Schryer.
c    This FORTRAN77 version by John Burkardt.
c
c  Reference:
c
c    Phyllis Fox, Andrew Hall, Norman Schryer,
c    Algorithm 528:
c    Framework for a Portable Library,
c    ACM Transactions on Mathematical Software,
c    Volume 4, Number 2, June 1978, page 176-188.
c
c  Parameters:
c
c    Input, integer I, the index of the desired constant.
c
c    Output, double precision D1MACH, the value of the constant.
c
      implicit none

      double precision d1mach
      integer i

      if ( i < 1 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'D1MACH - Fatal error!'
        write ( *, '(a)' ) '  The input argument I is out of bounds.'
        write ( *, '(a)' ) '  Legal values satisfy 1 <= I <= 5.'
        write ( *, '(a,i12)' ) '  I = ', i
        d1mach = 0.0D+00
        stop
      else if ( i == 1 ) then
        d1mach = 4.450147717014403D-308
      else if ( i == 2 ) then
        d1mach = 8.988465674311579D+307
      else if ( i == 3 ) then
        d1mach = 1.110223024625157D-016
      else if ( i == 4 ) then
        d1mach = 2.220446049250313D-016
      else if ( i == 5 ) then
        d1mach = 0.301029995663981D+000
      else if ( 5 < i ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'D1MACH - Fatal error!'
        write ( *, '(a)' ) '  The input argument I is out of bounds.'
        write ( *, '(a)' ) '  Legal values satisfy 1 <= I <= 5.'
        write ( *, '(a,i12)' ) '  I = ', i
        d1mach = 0.0D+00
        stop
      end if

      return
      end
      function i1mach ( i )

c*********************************************************************72
c
cc I1MACH returns integer machine dependent constants.
c
c  Discussion:
c
c    Input/output unit numbers.
c
c      I1MACH(1) = the standard input unit.
c      I1MACH(2) = the standard output unit.
c      I1MACH(3) = the standard punch unit.
c      I1MACH(4) = the standard error message unit.
c
c    Words.
c
c      I1MACH(5) = the number of bits per integer storage unit.
c      I1MACH(6) = the number of characters per integer storage unit.
c
c    Integers.
c
c    Assume integers are represented in the S digit base A form:
c
c      Sign * (X(S-1)*A**(S-1) + ... + X(1)*A + X(0))
c
c    where 0 <= X(1:S-1) < A.
c
c      I1MACH(7) = A, the base.
c      I1MACH(8) = S, the number of base A digits.
c      I1MACH(9) = A**S-1, the largest integer.
c
c    Floating point numbers
c
c    Assume floating point numbers are represented in the T digit 
c    base B form:
c
c      Sign * (B**E) * ((X(1)/B) + ... + (X(T)/B**T) )
c
c    where 0 <= X(I) < B for I=1 to T, 0 < X(1) and EMIN <= E <= EMAX.
c
c      I1MACH(10) = B, the base.
c
c    Single precision
c
c      I1MACH(11) = T, the number of base B digits.
c      I1MACH(12) = EMIN, the smallest exponent E.
c      I1MACH(13) = EMAX, the largest exponent E.
c
c    Double precision
c
c      I1MACH(14) = T, the number of base B digits.
c      I1MACH(15) = EMIN, the smallest exponent E.
c      I1MACH(16) = EMAX, the largest exponent E.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    25 April 2007
c
c  Author:
c
c    Original FORTRAN77 version by Phyllis Fox, Andrew Hall, Norman Schryer.
c    This FORTRAN77 version by John Burkardt.
c
c  Reference:
c
c    Phyllis Fox, Andrew Hall, Norman Schryer,
c    Algorithm 528,
c    Framework for a Portable Library,
c    ACM Transactions on Mathematical Software,
c    Volume 4, Number 2, June 1978, page 176-188.
c
c  Parameters:
c
c    Input, integer I, chooses the parameter to be returned.
c    1 <= I <= 16.
c
c    Output, integer I1MACH, the value of the chosen parameter.
c
      implicit none

      integer i
      integer i1mach

      if ( i < 1 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'I1MACH - Fatal error!'
        write ( *, '(a)' ) '  The input argument I is out of bounds.'
        write ( *, '(a)' ) '  Legal values satisfy 1 <= I <= 16.'
        write ( *, '(a,i12)' ) '  I = ', i
        i1mach = 0
        stop
      else if ( i == 1 ) then
        i1mach = 5
      else if ( i == 2 ) then
        i1mach = 6
      else if ( i == 3 ) then
        i1mach = 7
      else if ( i == 4 ) then
        i1mach = 6
      else if ( i == 5 ) then
        i1mach = 32
      else if ( i == 6 ) then
        i1mach = 4
      else if ( i == 7 ) then
        i1mach = 2
      else if ( i == 8 ) then
        i1mach = 31
      else if ( i == 9 ) then
        i1mach = 2147483647
      else if ( i == 10 ) then
        i1mach = 2
      else if ( i == 11 ) then
        i1mach = 24
      else if ( i == 12 ) then
        i1mach = -125
      else if ( i == 13 ) then
        i1mach = 128
      else if ( i == 14 ) then
        i1mach = 53
      else if ( i == 15 ) then
        i1mach = -1021
      else if ( i == 16 ) then
        i1mach = 1024
      else if ( 16 < i ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'I1MACH - Fatal error!'
        write ( *, '(a)' ) '  The input argument I is out of bounds.'
        write ( *, '(a)' ) '  Legal values satisfy 1 <= I <= 16.'
        write ( *, '(a,i12)' ) '  I = ', i
        i1mach = 0
        stop
      end if

      return
      end
      function r1mach ( i )

c*********************************************************************72
c
cc R1MACH returns single precision real machine constants.
c
c  Discussion:
c
c    Assume that single precision real numbers are stored with a mantissa 
c    of T digits in base B, with an exponent whose value must lie 
c    between EMIN and EMAX.  Then for values of I between 1 and 5, 
c    R1MACH will return the following values:
c
c      R1MACH(1) = B**(EMIN-1), the smallest positive magnitude.
c      R1MACH(2) = B**EMAX*(1-B**(-T)), the largest magnitude.
c      R1MACH(3) = B**(-T), the smallest relative spacing.
c      R1MACH(4) = B**(1-T), the largest relative spacing.
c      R1MACH(5) = log10(B)
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    25 April 2007
c
c  Author:
c
c    Original FORTRAN77 version by Phyllis Fox, Andrew Hall, Norman Schryer.
c    This FORTRAN77 version by John Burkardt.
c
c  Reference:
c
c    Phyllis Fox, Andrew Hall, Norman Schryer,
c    Algorithm 528,
c    Framework for a Portable Library,
c    ACM Transactions on Mathematical Software,
c    Volume 4, Number 2, June 1978, page 176-188.
c
c  Parameters:
c
c    Input, integer I, chooses the parameter to be returned.
c    1 <= I <= 5.
c
c    Output, real R1MACH, the value of the chosen parameter.
c
      implicit none

      integer i
      real r1mach

      if ( i < 1 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'R1MACH - Fatal error!'
        write ( *, '(a)' ) '  The input argument I is out of bounds.'
        write ( *, '(a)' ) '  Legal values satisfy 1 <= I <= 5.'
        write ( *, '(a,i12)' ) '  I = ', i
        r1mach = 0.0E+00
        stop
      else if ( i == 1 ) then
        r1mach = 1.1754944E-38
      else if ( i == 2 ) then
        r1mach = 3.4028235E+38
      else if ( i == 3 ) then
        r1mach = 5.9604645E-08
      else if ( i == 4 ) then
        r1mach = 1.1920929E-07
      else if ( i == 5 ) then
        r1mach = 0.3010300E+00
      else if ( 5 < i ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'R1MACH - Fatal error!'
        write ( *, '(a)' ) '  The input argument I is out of bounds.'
        write ( *, '(a)' ) '  Legal values satisfy 1 <= I <= 5.'
        write ( *, '(a,i12)' ) '  I = ', i
        r1mach = 0.0E+00
        stop
      end if

      return
      end
