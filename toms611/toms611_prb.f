      program main

c     ***  simple test program for smsno, sumsl, and humsl  ***
c
c  in these examples,  n = 4  and  f(x) = (1.0 + 0.5*(x1**t)*a*x1)**0.5,
c  where  x1(i) = d1(i)*x(i) - i,  **t denotes transpose, and  a  is a
c  matrix having fives on the main diagonal and ones everywhere else.
c  the scale vector  d1  is passed to qdrtf, the subroutine that
c  evaluates  f,  as part of urparm.  specifically, the matrix  urp
c  declared below is passed for ufparm, and  d1  is urp(*,1), the first
c  column of urp.  this main program repeatedly minimizes f, starting
c  from  x = 0,  by calling smsno, sumsl, and humsl.  we actually
c  use two different objective functions, since we change  d1  after
c  the first call on sumsl.  all runs but the last use  d = d1.
c
c     f(x) is minimized at  x1 = 0  (a vector of zeros), i.e., at
c  x(i) = i/d1(i).
c
      external deflt, humsl, imdcon, qdrtf, qdrtg, qdrtgh, smsno,
     1         sumsl
      integer imdcon
c
c imdcon - supplies nout, the output unit number.
c qdrtf  - passed for calcf to smsno, sumsl, and humsl.
c qdrtg  - passed for calcg to sumsl.
c qdrtgh - passed for calcgh to humsl.
c
      integer i, iv(60), liv, lv, nout, uip(1)
      double precision d(4), urp(4,3), v(150), x(4)
c
c  the length of  v  is dictated by humsl.
c
      data liv/60/, lv/150/
c
c  initialize nout, d, d1, and x.
c
      nout = imdcon(1)
c
      do 10 i = 1, 4
         d(i) = 1.d+0
         urp(i,1) = d(i)
         x(i) = 0.d+0
 10      continue
      write(nout,20)
 20   format(/16h smsno on qdrtf )
c
c  before this first call, we set iv(1) to 0 so that all input
c  components of iv and v will be given default values.  before
c  subsequent calls, we set iv(1) to 12 so that the old input values
c  of iv and v are used.
c
c  qdrtf does not make use of ufparm.  in calling smsno, we
c  arbitrarily pass qdrtf for ufparm to satifsy the calling sequence.
c
      iv(1) = 0
      call smsno(4, d, x, qdrtf, iv, liv, lv, v, uip, urp, qdrtf)
c
c  we reinitialize x and minimize  f  again, this time using sumsl.
c  qdrtg, the subroutine passed for calcg, assumes that ufparm is qdrtf.
c
      do 30 i = 1, 4
 30      x(i) = 0.d+0
      write(nout,40)
 40   format(/16h sumsl on qdrtf )
c
      iv(1) = 12
      call sumsl(4, d, x, qdrtf, qdrtg, iv, liv, lv, v, uip,urp,qdrtf)
c
c  now we modify  f  by using a different choice of d1.  we still use
c  d = d1, so the performance of sumsl should stay the same -- only d
c  and the final x and gradient should be affected.
c
      do 50 i = 1, 4
         x(i) = 0.d+0
         d(i) = 1.d2 ** i
         urp(i,1) = d(i)
 50      continue
      write(nout,40)
c
      iv(1) = 12
      call sumsl(4, d, x, qdrtf, qdrtg, iv, liv, lv, v, uip,urp,qdrtf)
c
c  using the last choice of d and d1, we now use humsl to minimize f.
c  like qdrtg, qdrtgh assumes that ufparm is qdrtf.
c
      do 60 i = 1, 4
 60      x(i) = 0.d+0
      write(nout,70)
 70   format(/16h humsl on qdrtf )
c
      iv(1) = 12
      call humsl(4, d, x, qdrtf, qdrtgh, iv, liv, lv, v, uip,urp,qdrtf)
c
c  we repeat the last run with iv(dtype) = 1 and v(dinit) = 0.0, so
c  that humsl will determine  d  from the diagonal of the hessian.
c  this run also demonstrates the use of subroutine deflt and the
c  passing of nondefault parameters.  (since the iv and v input
c  components still have their default values at his point, it is not
c  really necessary to call deflt.  it is necessary to reset iv(1) to
c  12, however, and deflt does this for us.)
c
      write(nout,80)
 80   format(/18h humsl updating d )
c
      do 90 i = 1, 4
 90      x(i) = 0.d+0
c
      call deflt(2, iv, liv, lv, v)
      iv(16) = 1
      v(38) = 0.d+0
      call humsl(4, d, x, qdrtf, qdrtgh, iv, liv, lv, v, uip,urp,qdrtf)
c
 999  stop
      end
c
c***********************************************************************
c
c     q d r t f
c
c***********************************************************************
c
      subroutine qdrtf(n, x, nf, f, uip, urp, ufp)
c
c  this routine evaluates the objective function f(x) described in the
c  main program above.  it stores in urp(*,2) and urp(*,3) some
c  information useful in evaluating the gradient and hessian of f, and
c  it stores  nf  in uip(1) to identify the x corresponding to this
c  information.  f(x) has the form  f(x) = phi(q(x)),  where  q  is a
c  quadratic form and  phi(y) = y**0.5.  the gradient of  f  is
c  g(x) = phiprm(q(x))*gq(x),  where  phiprm  is the derivative of phi
c  and  gq  is the gradient of q.  this routine stores phiprm(q(x)) in
c  urp(1,3) and gq(x) in urp(*,2).  the hessian of f is
c  h(x) = phi2prm(q(x))*gq(x)*gq(x)**t + phiprm(q(x))*hq(x),  where
c  phi2prm  is the second derivative of phi, **t denotes transpose,
c  and  hq  is the hessian of q.  this routine stores phi2prm(q(x)) in
c  urp(2,3).  the subroutines qdrtg and qdrtgh given below would work
c  without change on any other choice of phi.  qdrtg would also work
c  with any other differentiable function q.  qdrtgh, on the other
c  hand, assumes that  hq(x)  is the matrix  a  described in the main
c  program above.
c
      integer n, nf, uip(1)
      double precision x(n), f, urp(n,3)
      external ufp
c/+
      real float
      double precision dsqrt
c/
      integer i
      double precision dn, f2, t, t1
c
c
      uip(1) = nf
      dn = n
      t = 0.d+0
      do 10 i = 1, n
         urp(i,2) = urp(i,1)*x(i) - float(i)
         t = t + urp(i,2)
 10      continue
      f2 = 0.d+0
      do 20 i = 1, n
         t1 = dn*urp(i,2) + t
         f2 = f2 + t1*urp(i,2)
         urp(i,2) = urp(i,1) * t1
 20      continue
      f2 = 1.d+0  +   0.5d+0 * f2
      f = dsqrt(f2)
      urp(1,3) = 0.5d+0 / f
      urp(2,3) = -0.5d+0 / (f * f2)
 999  return
      end
c
c***********************************************************************
c
c     q d r t g
c
c***********************************************************************
c
      subroutine qdrtg(n, x, nf, g, uip, urp, qdrtf)
c
c  this routine evaluates the gradient of the objective function f(x)
c  described in the main program above.  see the comments there and in
c  subroutine qdrtf above.
c
      integer n, nf, uip(1)
      double precision x(n), g(n), urp(n,3)
      external qdrtf
c
      integer i
      double precision f
c
      if (nf .ne. uip(1)) call qdrtf(n, x, nf, f, uip, urp, qdrtf)
      do 10 i = 1, n
 10      g(i) = urp(1,3) * urp(i,2)
 999  return
      end
c
c***********************************************************************
c
c     q d r t g h
c
c***********************************************************************
c
      subroutine qdrtgh(n, x, nf, g, h, uip, urp, qdrtf)
c
c  this routine evaluates the gradient and hessian of the objective
c  function f(x) described in the main program above.  see the comments
c  there and in subroutine qdrtf above.  note that the  h  returned is
c  the lower triangle of the hessian, stored row-wise.
c
      integer n, nf, uip(1)
      double precision x(n), g(n), h(1), urp(n,3)
c     dimension h(n*(n+1)/2)
      external qdrtf
c
      integer i, j, k
      double precision dn, f, t1, t2
c
      if (nf .ne. uip(1)) call qdrtf(n, x, nf, f, uip, urp, qdrtf)
      k = 0
      dn = n
      do 20 i = 1, n
         g(i) = urp(1,3) * urp(i,2)
         t1 = urp(1,3) * urp(i,1)
         t2 = urp(i,2) * urp(2,3)
         do 10 j = 1, i
              k = k + 1
              h(k) = t2*urp(j,2) + t1*urp(j,1)
 10           continue
         h(k) = h(k) + dn*urp(i,1)*t1
 20      continue
 999  return
      end
