      subroutine qinit ( n, betam, nptsq, qwork )

c*********************************************************************72
c
cc QINIT computes nodes and weights for Gauss-Jacobi quadrature.
c
c  License:
c
c    The package may be freely copied and used, but reference to
c    the paper by Nick Trefethen should be given in any publications arising
c    from its use.
c 
c  Modified:
c
c    07 December 2007
c
c  Author:
c
c    Nick Trefethen
c
c  Reference:
c
c    Nick Trefethen,
c    Numerical Computation of the Schwarz-Christoffel Transformation,
c    SIAM Journal of Scientific and Statistical Computing,
c    Volume 1, 1980, pages 82-102.
c
c calling sequence parameters: see comments in scsolv
c
c the work array qwork must be dimensioned at least nptsq * (2n+3).
c it is divided up into 2n+3 vectors of length nptsq: the first
c n+1 contain quadrature nodes on output, the next n+1 contain
c quadrature weights on output, and the final one is a
c scratch vector needed by gaussj.
c
      implicit double precision (a-b,d-h,o-v,x-y)
      implicit complex*16(c,w,z)
      dimension qwork(1),betam(n)
c
c for each finite vertex w(k), compute nodes and weights for
c one-sided gauss-jacobi quadrature along a curve beginning at z(k):
      iscr = nptsq*(2*n+2) + 1
      do 1 k = 1,n
        inodes = nptsq*(k-1) + 1
        iwts = nptsq*(n+k) + 1
    1   if (betam(k).gt.-1.d0) call gaussj(nptsq,0.d0,betam(k),
     &    qwork(iscr),qwork(inodes),qwork(iwts))
c
c compute nodes and weights for pure gaussian quadrature:
      inodes = nptsq*n + 1
      iwts = nptsq*(2*n+1) + 1
      call gaussj(nptsq,0.d0,0.d0,
     &  qwork(iscr),qwork(inodes),qwork(iwts))

      return
      end
      subroutine scsolv(iprint,iguess,tol,errest,n,c,z,wc,
     &   w,betam,nptsq,qwork)

c*********************************************************************72
c
cc SCSOLV computes parameters for the Schwarz-Christoffel mapping.
c
c  License:
c
c    The package may be freely copied and used, but reference to
c    the paper by Nick Trefethen should be given in any publications arising
c    from its use.
c 
c  Modified:
c
c    07 December 2007
c
c  Author:
c
c    Nick Trefethen
c
c  Reference:
c
c    Nick Trefethen,
c    Numerical Computation of the Schwarz-Christoffel Transformation,
c    SIAM Journal of Scientific and Statistical Computing,
c    Volume 1, 1980, pages 82-102.
c
c this subroutine computes the accessory parameters c and
c z(k) for the schwarz-christoffel transformation
c which sends the unit disk to the interior of the polygon
c w(1),...,w(n).  this mapping is of the form:
c
c                         z    n
c    w   =   wc  +   c * int prod (1-z/z(k))**betam(k) dz .    (1.2)
c                         0   k=1
c
c the image polygon may be unbounded; permitted angles lie in the
c range -3.le.betam(k).le.1.  w(n) and w(1) must be finite.
c we normalize by the conditions:
c
c    z(n) = 1                                                  (2.1)
c    w(0.0) = wc  (a point in the interior of the polygon)     (2.1)
c
c calling sequence:
c
c   iprint  -2,-1,0, or 1 for increasing amounts of output (input)
c
c   iguess  1 if an initial guess for z is supplied, otherwise 0
c           (input)
c
c   tol     desired accuracy in solution of nonlinear system
c           (input).  recommended value: 10.**(-nptsq-1) * typical
c           size of vertices w(k)
c
c   errest  estimtated error in solution (output).  more
c           precisely, errest is an approximate bound for how far
c           the true vertices of the image polygon may be from those
c           computed by numerical integration using the
c           numerically determined prevertices z(k).
c
c   n       number of vertices of the image polygon (input).
c           must be .le. 20
c
c   c       complex scale factor in formula above (output)
c
c   z       complex array of prevertices on the unit circle.
c           dimension at least n.  if an initial guess is
c           being supplied it should be in z on input, with z(n)=1.
c           in any case the correct prevertices will be in z on output.
c
c   wc      complex image of 0 in the polygon, as in above formula
c           (input).  it is safest to pick wc to be as central as
c           possible in the polygon in the sense that as few parts
c           of the polygon as possible are shielded from wc by
c           reentrant edges.
c
c   w       complex array of vertices of the image polygon
c           (input).  dimension at least n.  it is a good idea
c           to keep the w(k) roughly on the scale of unity.
c           w(k) will be ignored when the vertex lies at infinity;
c           see betam, below.  each connected boundary component
c           must include at least one vertex w(k), even if it
c           has to be a degenerate vertex with betam(k) = 0.
c           w(n) and w(1) must be finite.
c
c   betam   real array with betam(k) the external angle in the
c           polygon at vertex k divided by minus pi (input).
c           dimension at least n.  permitted values lie in
c           the range -3.le.betam(k).le.1.  (examples: each
c           betam(k) is -1/2 for a rectangle, -2/3 for an equi-
c           lateral triangle, +1 at the end of a slit.)  the
c           sum of the betam(k) will be -2 if they have been
c           set correctly.  betam(n-1) should not be 0 or 1.
c           w(k) lies at infinity if and only if betam(k).le.-1.
c
c   nptsq   the number of points to be used per subinterval
c           in gauss-jacobi quadrature (input).  recommended
c           value: equal to one more than the number of digits
c           of accuracy desired in the answer.  must be the same
c           as in the call to qinit which filled the vector qwork.
c
c   qwork   real quadrature work array (input).  dimension
c           at least nptsq * (2n+3) but no greater than 460.
c           before calling scsolv qwork must have been filled
c           by subroutine qinit.
c
c the problem is solved by finding the
c solution to a system of n-1 nonlinear equations in the n-1
c unknowns y(1),...,y(n-1), which are related to the points
c z(k) by the formula:
c
c      y(k) = log ((th(k)-th(k-1))/(th(k+1)-th(k)))            (2.7)
c
c where th(k) denotes the argument of z(k).
c subroutine scfun defines this system of equations.
c the original problem is subject to the contraints th(k) < th(k+1),
c but these vanish in the transformation from z to y.
c
c reference:  l. n. trefethen, "numerical computation of the
c   schwarz-christoffel transformation," siam j. sci. stat. comp. 1
c   (1980), 82-102.  equation nos. above are taken from this paper.
c
c lloyd n. trefethen
c courant institute of mathematical sciences, new york university
c 251 mercer st.,   new york, ny 10012
c (212) 460-7224
c october 1979 (version 1); july 1983 (version 2)
c
      implicit double precision (a-b,d-h,o-v,x-y)
      implicit complex*16(c,w,z)

      common /param1/ kfix(20),krat(20),ncomp,nptsq2,c2,
     &  qwork2(460),betam2(20),z2(20),wc2,w2(20)
      dimension z(n),w(n),betam(n),qwork(1)
      dimension ajinv(20,20),scr(900),fval(19),y(19)
      external scfun
      nm = n-1
c
c check input data:
      call check(tol,n,w,betam)
c
c determine number of boundary components, etc.:
c   pass 1: one fixed point for each infinite vertex:
      ncomp = 0
      do 1 k = 2,nm
        if (betam(k).gt.-1.d0) goto 1
        ncomp = ncomp + 1
        kfix(ncomp) = k - 1
        if (ncomp.eq.1) kfix(ncomp) = 1
    1 continue
      if (ncomp.gt.0) goto 2
      ncomp = 1
      kfix(ncomp) = 1
c   pass 2: one ratio for each line segment:
    2 continue
      neq = 2*ncomp
      do 3 k = 1,nm
        if (neq.eq.n-1) goto 4
        if (betam(k).le.-1.d0.or.betam(k+1).le.-1.d0) goto 3
        neq = neq + 1
        krat(neq) = k
    3 continue
    4 z(n) = (1.d0,0.d0)
c
c initial guess, case iguess.eq.0:
c (vertices equally spaced around circle):
      if (iguess.ne.0) goto 11
      do 5 k = 1,nm
    5   y(k) = 0.d0
      goto 12
c
c initial guess, case iguess.ne.0:
c (vertices supplied on input):
   11 continue
      do 9 k = 1,nm
        km = k-1
        if (km.eq.0) km = n
        tmp1 = dimag(log(z(k+1)/z(k)))
        if (tmp1.lt.0.d0) tmp1 = tmp1 + 2.d0 * acos(-1.d0)
        tmp2 = dimag(log(z(k)/z(km)))
        if (tmp2.lt.0.d0) tmp2 = tmp2 + 2.d0 * acos(-1.d0)
    9   y(k) = log(tmp2) - log(tmp1)
   12 continue
c
c ns01a control parameters:
      dstep = 1.d-6
      dmax = 20.d0
      maxfun = (n-1) * 15
c
c copy input data to /param1/ for scfun:
c (this is necessary because ns01a requires a fixed calling
c sequence in subroutine scfun.)
      nptsq2 = nptsq
      wc2 = wc
      do 6 k = 1,n
        z2(k) = z(k)
        betam2(k) = betam(k)
    6   w2(k) = w(k)
      nwdim = nptsq * (2*n+3)
      do 7 i = 1,nwdim
    7   qwork2(i) = qwork(i)
c
c solve nonlinear system with ns01a:
      call ns01a(nm,y,fval,ajinv,dstep,dmax,tol,maxfun,
     &  iprint,scr,scfun)
c
c copy output data from /param1/:
      c = c2
      do 8 k = 1,nm
    8   z(k) = z2(k)
c
c print results and test accuracy:
      if (iprint.ge.0) call scoutp(n,c,z,wc,w,betam,nptsq)
      call sctest(errest,n,c,z,wc,w,betam,nptsq,qwork)
      if (iprint.ge.-1) write (6,201) errest
  201 format (' errest:',e12.4/)
      return
      end
      function wsc(zz,kzz,z0,w0,k0,n,c,z,betam,nptsq,qwork)

c*********************************************************************72
c
cc WSC computes the forward map w(zz).
c
c  License:
c
c    The package may be freely copied and used, but reference to
c    the paper by Nick Trefethen should be given in any publications arising
c    from its use.
c 
c  Modified:
c
c    07 December 2007
c
c  Author:
c
c    Nick Trefethen
c
c  Reference:
c
c    Nick Trefethen,
c    Numerical Computation of the Schwarz-Christoffel Transformation,
c    SIAM Journal of Scientific and Statistical Computing,
c    Volume 1, 1980, pages 82-102.
c
c calling sequence:
c
c   zz     point in the disk at which w(zz) is desired (input)
c
c   kzz    k if zz = z(k) for some k, otherwise 0 (input)
c
c   z0     nearby point in the disk at which w(z0) is known and
c          finite (input)
c
c   w0     w(z0)  (input)
c
c   k0     k if z0 = z(k) for some k, otherwise 0 (input)
c
c   n,c,z,betam,nptsq,qwork     as in scsolv (input)
c
c convenient values of z0, w0, and k0 for most applications can be
c supplied by subroutine nearz.
c
      implicit double precision (a-b,d-h,o-v,x-y)
      implicit complex*16(c,w,z)

      dimension z(n),betam(n),qwork(1)

      wsc = w0 + c * zquad(z0,k0,zz,kzz,n,z,betam,nptsq,qwork)

      return
      end
      function zsc(ww,iguess,zinit,z0,w0,k0,eps,ier,n,c,
     &  z,wc,w,betam,nptsq,qwork)

c*********************************************************************72
c
cc ZSC computes the inverse map Z(WW) by Newton iteration.
c
c  License:
c
c    The package may be freely copied and used, but reference to
c    the paper by Nick Trefethen should be given in any publications arising
c    from its use.
c 
c  Modified:
c
c    07 December 2007
c
c  Author:
c
c    Nick Trefethen
c
c  Reference:
c
c    Nick Trefethen,
c    Numerical Computation of the Schwarz-Christoffel Transformation,
c    SIAM Journal of Scientific and Statistical Computing,
c    Volume 1, 1980, pages 82-102.
c
c calling sequence:
c
c   ww     point in the polygon at which z(ww) is desired (input)
c
c   iguess (input)
c          .eq.1 - initial guess is supplied as parameter zinit
c          .ne.1 - get initial guess from program ode (slower).
c                  for this the segment wc-ww must lie within
c                  the polygon.
c
c   zinit  initial guess if iguess.eq.1, otherwise ignored (input).
c          may not be a prevertex z(k)
c
c   z0     point in the disk near z(ww) at which w(z0) is known and
c          finite (input).
c
c   w0     w(z0)  (input).  the line segment from w0 to ww must
c          lie entirely within the closed polygon.
c
c   k0     k if z0 = z(k) for some k, otherwise 0 (input)
c
c   eps    desired accuracy in answer z(ww)  (input)
c
c   ier    error flag (input and output).
c          on input, give ier.ne.0 to suppress error messages.
c          on output, ier.ne.0 indicates unsuccessful computation --
c          try again with a better initial guess.
c
c   n,c,z,wc,w,betam,nptsq,qwork     as in scsolv (input)
c
c convenient values of z0, w0, and k0 for some applications can be
c supplied by subroutine nearw.
c
      implicit double precision (a-b,d-h,o-v,x-y)
      implicit complex*16(c,w,z)

      dimension scr(142),iscr(5)
      dimension z(n),w(n),betam(n),qwork(1)
      external zfode
      logical odecal
      common /param2/ cdwdt,z2(20),betam2(20),n2
c
      odecal = .false.
      if (iguess.ne.1) goto 1
      zi = zinit
      goto 3
c
c get initial guess zi from program ode:
    1 n2 = n
      do 2 k = 1,n
        z2(k) = z(k)
    2   betam2(k) = betam(k)
      zi = (0.d0,0.d0)
      t = 0.d0
      iflag = -1
      relerr = 0.d0
      abserr = 1.d-4
      cdwdt = (ww-wc)/c
      call ode(zfode,2,zi,t,1.d0,relerr,abserr,iflag,scr,iscr)
      if (iflag.ne.2.and.ier.eq.0) write (6,201) iflag
      odecal = .true.
c
c refine answer by newton iteration:
    3 continue
      do 4 iter = 1,10
        zfnwt = ww - wsc(zi,0,z0,w0,k0,n,c,z,betam,nptsq,qwork)
        zi = zi + zfnwt/(c*zprod(zi,0,n,z,betam))
        if (abs(zi).ge.1.1d0) zi = .5d0 * zi/abs(zi)
        if (abs(zfnwt).lt.eps) goto 5
    4   continue
      if (.not.odecal) goto 1
      if (ier.eq.0) write (6,202)
      ier = 1
    5 zsc = zi
c
  201 format (/' *** nonstandard return from ode in zsc: iflag =',i2/)
  202 format (/' *** possible error in zsc: no convergence in 10'/
     &      '     iterations.  may need a better initial guess zinit')
      return
      end
      subroutine zfode(t,zz,zdzdt)

c*********************************************************************72
c
cc ZFODE computes the function ZDZDT needed by ODE in ZSC.
c
c  License:
c
c    The package may be freely copied and used, but reference to
c    the paper by Nick Trefethen should be given in any publications arising
c    from its use.
c 
c  Modified:
c
c    07 December 2007
c
c  Author:
c
c    Nick Trefethen
c
c  Reference:
c
c    Nick Trefethen,
c    Numerical Computation of the Schwarz-Christoffel Transformation,
c    SIAM Journal of Scientific and Statistical Computing,
c    Volume 1, 1980, pages 82-102.
c
      implicit double precision (a-b,d-h,o-v,x-y)
      implicit complex*16(c,w,z)

      common /param2/ cdwdt,z(20),betam(20),n

      zdzdt = cdwdt / zprod(zz,0,n,z,betam)

      return
      end
      subroutine check(eps,n,w,betam)

c*********************************************************************72
c
cc CHECK checks geometry of the problem.
c
c  License:
c
c    The package may be freely copied and used, but reference to
c    the paper by Nick Trefethen should be given in any publications arising
c    from its use.
c 
c  Modified:
c
c    07 December 2007
c
c  Author:
c
c    Nick Trefethen
c
c  Reference:
c
c    Nick Trefethen,
c    Numerical Computation of the Schwarz-Christoffel Transformation,
c    SIAM Journal of Scientific and Statistical Computing,
c    Volume 1, 1980, pages 82-102.
c
      implicit double precision (a-b,d-h,o-v,x-y)
      implicit complex*16(c,w,z)

      dimension w(n),betam(n)

      sum = 0.d0
      do 1 k = 1,n
    1   sum = sum + betam(k)
      if (abs(sum+2.d0).lt.eps) goto 2
      write (6,301)
    2 if (betam(1).gt.-1.d0) goto 3
      write (6,302)
      stop
    3 if (betam(n).gt.-1.d0) goto 4
      write (6,303)
      stop
    4 if (abs(betam(n-1)).gt.eps) goto 5
      write (6,304)
      write (6,306)
    5 if (abs(betam(n-1)-1.d0).gt.eps) goto 6
      write (6,305)
      write (6,306)
      stop
    6 do 7 k = 2,n
        if (betam(k).le.-1.d0.or.betam(k-1).le.-1.d0) goto 7
        if (abs(w(k)-w(k-1)).le.eps) goto 8
    7   continue
      if (abs(w(1)-w(n)).gt.eps) goto 9
    8 write (6,307)
      stop
    9 if (n.ge.3) goto 10
      write (6,309)
      stop
   10 if (n.le.20) goto 11
      write (6,310)
      stop
   11 continue
      return
c
  301 format (/' *** error in check: angles do not add up to 2'/)
  302 format (/' *** error in check: w(1) must be finite'/)
  303 format (/' *** error in check: w(n) must be finite'/)
  304 format (/' *** warning in check: w(n-1) not determined'/)
  305 format (/' *** error in check: w(n-1) not determined')
  306 format (/'   renumber vertices so that betam(n-1) is not 0 or 1')
  307 format (/' *** error in check: two adjacent vertices are equal'/)
  309 format (/' *** error in check: n must be no less than 3'/)
  310 format (/' *** error in check: n must be no more than 20'/)
      end
      subroutine yztran(n,y,z)

c*********************************************************************72
c
cc YZTRAN transforms y(k) to z(k).
c
c  License:
c
c    The package may be freely copied and used, but reference to
c    the paper by Nick Trefethen should be given in any publications arising
c    from its use.
c 
c  Modified:
c
c    07 December 2007
c
c  Author:
c
c    Nick Trefethen
c
c  Reference:
c
c    Nick Trefethen,
c    Numerical Computation of the Schwarz-Christoffel Transformation,
c    SIAM Journal of Scientific and Statistical Computing,
c    Volume 1, 1980, pages 82-102.
c
c  Discussion:
c
c    See comments in subroutine scsolv.
c
      implicit double precision (a-b,d-h,o-v,x-y)
      implicit complex*16(c,w,z)

      dimension y(n),z(n)
      nm = n - 1
      pi = acos(-1.d0)
c
      dth = 1.d0
      thsum = dth
      do 1 k = 1,nm
        dth = dth / exp(y(k))
    1   thsum = thsum + dth
c
      dth = 2.d0 * pi / thsum
      thsum = dth
      z(1) = dcmplx(cos(dth),sin(dth))
      do 2 k = 2,nm
        dth = dth / exp(y(k-1))
        thsum = thsum + dth
    2   z(k) = dcmplx(cos(thsum),sin(thsum))

      return
      end
      subroutine scfun(ndim,y,fval)

c*********************************************************************72
c
cc SCFUN is the function whose zero must be found in SCSOLV.
c
c  License:
c
c    The package may be freely copied and used, but reference to
c    the paper by Nick Trefethen should be given in any publications arising
c    from its use.
c 
c  Modified:
c
c    07 December 2007
c
c  Author:
c
c    Nick Trefethen
c
c  Reference:
c
c    Nick Trefethen,
c    Numerical Computation of the Schwarz-Christoffel Transformation,
c    SIAM Journal of Scientific and Statistical Computing,
c    Volume 1, 1980, pages 82-102.
c
      implicit double precision (a-b,d-h,o-v,x-y)
      implicit complex*16(c,w,z)

      dimension fval(ndim),y(ndim)
      common /param1/ kfix(20),krat(20),ncomp,nptsq,c,
     &  qwork(460),betam(20),z(20),wc,w(20)
      n = ndim+1
c
c transform y(k) to z(k):
      call yztran(n,y,z)
c
c set up: compute integral from 0 to z(n):
      wdenom = zquad((0.d0,0.d0),0,z(n),n,n,z,betam,nptsq,qwork)
      c = (w(n)-wc) / wdenom
c
c case 1: w(k) and w(k+1) finite:
c (compute integral along chord z(k)-z(k+1)):
      nfirst = 2*ncomp + 1
      if (nfirst.gt.n-1) goto 11
      do 10 neq = nfirst,ndim
        kl = krat(neq)
        kr = kl+1
        zint = zquad(z(kl),kl,z(kr),kr,n,z,betam,nptsq,qwork)
        fval(neq) = abs(w(kr)-w(kl)) - abs(c*zint)
   10 continue
c
c case 2: w(k+1) infinite:
c (compute contour integral along radius 0-z(k)):
   11 do 20 nvert = 1,ncomp
        kr = kfix(nvert)
        zint = zquad((0.d0,0.d0),0,z(kr),kr,n,z,betam,nptsq,qwork)
        zfval = w(kr) - wc - c*zint
        fval(2*nvert-1) = dreal(zfval)
        fval(2*nvert) = dimag(zfval)
   20 continue
      return
      end
      subroutine scoutp(n,c,z,wc,w,betam,nptsq)

c*********************************************************************72
c
cc SCOUTP prints K, W(K), TH(K), BETAM(K) and Z(K).
c
c  Discussion:
c
c    The routine also prints the constants n, nptsq, wc, c.
c
c  License:
c
c    The package may be freely copied and used, but reference to
c    the paper by Nick Trefethen should be given in any publications arising
c    from its use.
c 
c  Modified:
c
c    07 December 2007
c
c  Author:
c
c    Nick Trefethen
c
c  Reference:
c
c    Nick Trefethen,
c    Numerical Computation of the Schwarz-Christoffel Transformation,
c    SIAM Journal of Scientific and Statistical Computing,
c    Volume 1, 1980, pages 82-102.
c
      implicit double precision (a-b,d-h,o-v,x-y)
      implicit complex*16(c,w,z)

      dimension z(n),w(n),betam(n)
c
      write (6,102) n, nptsq
      pi = acos(-1.d0)
      do 1 k = 1,n
        thdpi = dimag(log(z(k))) / pi
        if (thdpi.le.0.d0) thdpi = thdpi + 2.d0
        if (betam(k).gt.-1.d0) write (6,103) k,w(k),thdpi,betam(k),z(k)
    1   if (betam(k).le.-1.d0) write (6,104) k,thdpi,betam(k),z(k)
      write (6,105) wc,c
      return
c
  102 format (/' parameters defining map:',15x,'(n =',
     &  i3,')',6x,'(nptsq =',i3,')'//
     &  '  k',10x,'w(k)',10x,'th(k)/pi',5x,'betam(k)',
     &  13x,'z(k)'/
     &  ' ---',9x,'----',10x,'--------',5x,'--------',
     &  13x,'----'/)
  103 format (i3,'    (',f6.3,',',f6.3,')',f14.8,f12.5,
     &  3x,'(',f11.8,',',f11.8,')')
  104 format (i3,'        infinity   ',f14.8,f12.5,
     &  3x,'(',f11.8,',',f11.8,')')
  105 format (/' wc = (',e15.8,',',e15.8,')'/
     &          '  c = (',e15.8,',',e15.8,')'/)
      end
      subroutine sctest(errest,n,c,z,wc,w,betam,nptsq,qwork)

c*********************************************************************72
c
cc SCTEST tests the computed map for accuracy.
c
c  Modified:
c  License:
c
c    The package may be freely copied and used, but reference to
c    the paper by Nick Trefethen should be given in any publications arising
c    from its use.
c 
c  Modified:
c
c    07 December 2007
c
c  Author:
c
c    Nick Trefethen
c
c  Reference:
c
c    Nick Trefethen,
c    Numerical Computation of the Schwarz-Christoffel Transformation,
c    SIAM Journal of Scientific and Statistical Computing,
c    Volume 1, 1980, pages 82-102.
c
      implicit double precision (a-b,d-h,o-v,x-y)
      implicit complex*16(c,w,z)

      dimension z(n),w(n),betam(n),qwork(1)
c
c test length of radii:
      errest = 0.d0
      do 10 k = 2,n
        if (betam(k).gt.-1.d0) rade = abs(wc -
     &    wsc((0.d0,0.d0),0,z(k),w(k),k,n,c,z,betam,nptsq,qwork))
        if (betam(k).le.-1.d0) rade = abs(wsc((.1d0,.1d0),0,
     &    z(k-1),w(k-1),k-1,n,c,z,betam,nptsq,qwork)
     &    - wsc((.1d0,.1d0),0,z(k+1),w(k+1),k+1,
     &    n,c,z,betam,nptsq,qwork))
        errest = max(errest,rade)
   10   continue
      return
      end
      function zquad(za,ka,zb,kb,n,z,betam,nptsq,qwork)

c*********************************************************************72
c
cc ZQUAD computes the line integral of ZPROD.
c
c  Discussion:
c
c    computes the complex line integral of zprod from za to zb along a
c    straight line segment within the unit disk.  function zquad1 is
c    called twice, once for each half of this integral.
c
c  License:
c
c    The package may be freely copied and used, but reference to
c    the paper by Nick Trefethen should be given in any publications arising
c    from its use.
c 
c  Modified:
c
c    07 December 2007
c
c  Author:
c
c    Nick Trefethen
c
c  Reference:
c
c    Nick Trefethen,
c    Numerical Computation of the Schwarz-Christoffel Transformation,
c    SIAM Journal of Scientific and Statistical Computing,
c    Volume 1, 1980, pages 82-102.
c
      implicit double precision (a-b,d-h,o-v,x-y)
      implicit complex*16(c,w,z)

      dimension z(n),betam(n),qwork(1)
c
      if (abs(za).gt.1.1d0.or.abs(zb).gt.1.1d0) write (6,301)
  301 format (/' *** warning in zquad: z outside the disk')
c
      zmid = (za + zb) / 2.d0
      zquad = zquad1(za,zmid,ka,n,z,betam,nptsq,qwork)
     &  - zquad1(zb,zmid,kb,n,z,betam,nptsq,qwork)
      return
      end
      function zquad1(za,zb,ka,n,z,betam,nptsq,qwork)

c*********************************************************************72
c
cc ZQUAD1 computes half of the line integral of ZPROD.
c
c  Discussion:
c
c    The routine computes the complex line integral of zprod from za to zb along a
c    straight line segment within the unit disk.  compound one-sided
c    gauss-jacobi quadrature is used, using function dist to determine
c    the distance to the nearest singularity z(k).
c
c  License:
c
c    The package may be freely copied and used, but reference to
c    the paper by Nick Trefethen should be given in any publications arising
c    from its use.
c 
c  Modified:
c
c    07 December 2007
c
c  Author:
c
c    Nick Trefethen
c
c  Reference:
c
c    Nick Trefethen,
c    Numerical Computation of the Schwarz-Christoffel Transformation,
c    SIAM Journal of Scientific and Statistical Computing,
c    Volume 1, 1980, pages 82-102.
c
      implicit double precision (a-b,d-h,o-v,x-y)
      implicit complex*16(c,w,z)

      dimension z(n),betam(n),qwork(1)

      data resprm /1.d0/
c
c check for zero-length integrand:
      if (abs(za-zb).gt.0.d0) goto 1
      zquad1 = (0.d0,0.d0)
      return
c
c step 1: one-sided gauss-jacobi quadrature for left endpoint:
    1 r = min(1.d0,dist(za,ka,n,z)*resprm/abs(za-zb))
      zaa = za + r*(zb-za)
      zquad1 = zqsum(za,zaa,ka,n,z,betam,nptsq,qwork)
c
c step 2: adjoin intervals of pure gaussian quadrature if necessary:
   10 if (r.eq.1.d0) return
      r = min(1.d0,dist(zaa,0,n,z)*resprm/abs(zaa-zb))
      zbb = zaa + r*(zb-zaa)
      zquad1 = zquad1 + zqsum(zaa,zbb,0,n,z,betam,nptsq,qwork)
      zaa = zbb
      goto 10
      end
      function dist(zz,ks,n,z)

c*********************************************************************72
c
cc DIST determines the distance from zz to the nearest singularity z(k).
c
c  Discussion:
c
c    The singularity Z(KS) is to be ignored.
c
c  License:
c
c    The package may be freely copied and used, but reference to
c    the paper by Nick Trefethen should be given in any publications arising
c    from its use.
c 
c  Modified:
c
c    07 December 2007
c
c  Author:
c
c    Nick Trefethen
c
c  Reference:
c
c    Nick Trefethen,
c    Numerical Computation of the Schwarz-Christoffel Transformation,
c    SIAM Journal of Scientific and Statistical Computing,
c    Volume 1, 1980, pages 82-102.
c
      implicit double precision (a-b,d-h,o-v,x-y)
      implicit complex*16(c,w,z)
      dimension z(n)

      dist = 99.d0
      do 1 k = 1,n
        if (k.eq.ks) goto 1
        dist = min(dist,abs(zz-z(k)))
    1   continue
      return
      end
      function zqsum(za,zb,ka,n,z,betam,nptsq,qwork)

c*********************************************************************72
c
cc ZQSUM computes the integral of ZPROD.
c
c computes the integral of zprod from za to zb by applying a
c one-sided gauss-jacobi formula with possible singularity at za.
c
c  License:
c
c    The package may be freely copied and used, but reference to
c    the paper by Nick Trefethen should be given in any publications arising
c    from its use.
c 
c  Modified:
c
c    07 December 2007
c
c  Author:
c
c    Nick Trefethen
c
c  Reference:
c
c    Nick Trefethen,
c    Numerical Computation of the Schwarz-Christoffel Transformation,
c    SIAM Journal of Scientific and Statistical Computing,
c    Volume 1, 1980, pages 82-102.
c
      implicit double precision (a-b,d-h,o-v,x-y)
      implicit complex*16(c,w,z)
      dimension z(n),betam(n),qwork(1)
c
      zs = (0.d0,0.d0)
      zh = (zb-za) / 2.d0
      zc = (za+zb) / 2.d0
      k = ka
      if (k.eq.0) k = n+1
      iwt1 = nptsq*(k-1) + 1
      iwt2 = iwt1 + nptsq - 1
      ioffst = nptsq*(n+1)
      do 1 i = iwt1,iwt2
    1   zs = zs + qwork(ioffst+i)*zprod(zc+zh*qwork(i),ka,n,z,betam)
      zqsum = zs*zh
      if (abs(zh).ne.0.d0.and.k.ne.n+1)
     &  zqsum = zqsum*abs(zh)**betam(k)
      return
      end
      function zprod(zz,ks,n,z,betam)

c*********************************************************************72
c
cc ZPROD computes a Schwarz-Christoffel integrand.
c
c  The integrand is
c
c           n
c         prod  (1-zz/z(k))**betam(k)  ,
c          k=1
c
c taking argument only (not modulus) for term k = ks.
c
c  note -- in practice this is the innermost subroutine
c  in scpack calculations.  the complex log calculation below
c  may account for as much as half of the total execution time.
c
c  License:
c
c    The package may be freely copied and used, but reference to
c    the paper by Nick Trefethen should be given in any publications arising
c    from its use.
c 
c  Modified:
c
c    07 December 2007
c
c  Author:
c
c    Nick Trefethen
c
c  Reference:
c
c    Nick Trefethen,
c    Numerical Computation of the Schwarz-Christoffel Transformation,
c    SIAM Journal of Scientific and Statistical Computing,
c    Volume 1, 1980, pages 82-102.
c
      implicit double precision (a-b,d-h,o-v,x-y)
      implicit complex*16(c,w,z)
      dimension z(n),betam(n)
*/*/* new line:
      common /logcnt/ ncount,ncnt1
c
      zsum = (0.d0,0.d0)
      do 1 k = 1,n
        ztmp = (1.d0,0.d0) - zz/z(k)
        if (k.eq.ks) ztmp = ztmp / abs(ztmp)
    1   zsum = zsum + betam(k)*log(ztmp)
      zprod = exp(zsum)
*/*/* new line:
      ncount = ncount + n
      return
      end
      function rprod(zz,n,z,betam)

c*********************************************************************72
c
cc RPROD computes the absolute value of a Schwarz-Christoffel integrand.
c
c  The integrand is
c
c           n
c         prod  (1-zz/z(k))**betam(k)
c          k=1
c
c  License:
c
c    The package may be freely copied and used, but reference to
c    the paper by Nick Trefethen should be given in any publications arising
c    from its use.
c 
c  Modified:
c
c    07 December 2007
c
c  Author:
c
c    Nick Trefethen
c
c  Reference:
c
c    Nick Trefethen,
c    Numerical Computation of the Schwarz-Christoffel Transformation,
c    SIAM Journal of Scientific and Statistical Computing,
c    Volume 1, 1980, pages 82-102.
c
      implicit double precision (a-b,d-h,o-v,x-y)
      implicit complex*16(c,w,z)
      dimension z(n),betam(n)

      sum = 0.d0
      do 1 k = 1,n
        ztmp = (1.d0,0.d0) - zz/z(k)
    1   sum = sum + betam(k)*log(abs(ztmp))
      rprod = exp(sum)
      return
      end
      subroutine nearz(zz,zn,wn,kn,n,z,wc,w,betam)

c*********************************************************************72
c
cc NEARZ returns information about the nearest prevertex.
c
c returns information associated with the nearest prevertex z(k)
c to the point zz, or with 0 if 0 is closer than any z(k).
c zn = prevertex position, wn = w(zn), kn = prevertex no. (0 to n)
c
c  License:
c
c    The package may be freely copied and used, but reference to
c    the paper by Nick Trefethen should be given in any publications arising
c    from its use.
c 
c  Modified:
c
c    07 December 2007
c
c  Author:
c
c    Nick Trefethen
c
c  Reference:
c
c    Nick Trefethen,
c    Numerical Computation of the Schwarz-Christoffel Transformation,
c    SIAM Journal of Scientific and Statistical Computing,
c    Volume 1, 1980, pages 82-102.
c
      implicit double precision (a-b,d-h,o-v,x-y)
      implicit complex*16(c,w,z)
      dimension z(n),w(n),betam(n)

      dist = abs(zz)
      kn = 0
      zn = (0.d0,0.d0)
      wn = wc
      if (dist.le..5d0) return

      do k = 1,n
        if (betam(k).gt.-1.d0) then
          distzk = abs(zz-z(k))
          if (distzk.lt.dist) then
            dist = distzk
            kn = k
          end if
        end if
      end do

      if (kn.eq.0) return
      zn = z(kn)
      wn = w(kn)
      return
      end
      subroutine nearw(ww,zn,wn,kn,n,z,wc,w,betam)

c*********************************************************************72
c
c  NEARW returns information about the nearest vertex.
c
c  Discussion:
c
c    The routine returns information associated with the nearest vertex w(k)
c    to the point ww, or with wc if wc is closer than any w(k).
c    zn = prevertex position, wn = w(zn), kn = vertex no. (0 to n)
c
c  License:
c
c    The package may be freely copied and used, but reference to
c    the paper by Nick Trefethen should be given in any publications arising
c    from its use.
c 
c  Modified:
c
c    07 December 2007
c
c  Author:
c
c    Nick Trefethen
c
c  Reference:
c
c    Nick Trefethen,
c    Numerical Computation of the Schwarz-Christoffel Transformation,
c    SIAM Journal of Scientific and Statistical Computing,
c    Volume 1, 1980, pages 82-102.
c
      implicit double precision (a-b,d-h,o-v,x-y)
      implicit complex*16(c,w,z)
      dimension z(n),w(n),betam(n)

      dist = abs(ww-wc)
      kn = 0
      zn = (0.d0,0.d0)
      wn = wc

      do k = 1,n
        if (betam(k).le.-1.d0) goto 1
        distwk = abs(ww-w(k))
        if (distwk.ge.dist) goto 1
        dist = distwk
        kn = k
    1   continue
      end do

      if (kn.eq.0) return
      zn = z(kn)
      wn = w(kn)
      return
      end
      subroutine angles ( n, w, betam )

c*********************************************************************72
c
cc ANGLES computes external angles.
c
c  Discussion:
c
c    computes external angles -pi*betam(k) from knowledge of
c    the vertices w(k).  an angle betam(k) is computed for each
c    k for which w(k-1), w(k), and w(k+1) are finite.
c    to get this information across any vertices at infinity
c    should be signaled by the value w(k) = (99.,99.) on input.
c
c  License:
c
c    The package may be freely copied and used, but reference to
c    the paper by Nick Trefethen should be given in any publications arising
c    from its use.
c 
c  Modified:
c
c    07 December 2007
c
c  Author:
c
c    Nick Trefethen
c
c  Reference:
c
c    Nick Trefethen,
c    Numerical Computation of the Schwarz-Christoffel Transformation,
c    SIAM Journal of Scientific and Statistical Computing,
c    Volume 1, 1980, pages 82-102.
c
      implicit double precision (a-b,d-h,o-v,x-y)
      implicit complex*16(c,w,z)

      dimension w(n),betam(n)
      c9 = (99.d0,99.d0)

      pi = acos(-1.d0)
      do 1 k = 1,n
        km = mod(k+n-2,n)+1
        kp = mod(k,n)+1
        if (w(km).eq.c9.or.w(k).eq.c9.or.w(kp).eq.c9) goto 1
        betam(k) = dimag(log((w(km)-w(k))/(w(kp)-w(k))))/pi - 1.d0
        if (betam(k).le.-1.d0) betam(k) = betam(k) + 2.d0
    1   continue
      return
      end
      subroutine count0

c*********************************************************************72
c
cc COUNT0 initializes internal counters.
c
c  License:
c
c    The package may be freely copied and used, but reference to
c    the paper by Nick Trefethen should be given in any publications arising
c    from its use.
c 
c  Modified:
c
c    07 December 2007
c
c  Author:
c
c    Nick Trefethen
c
c  Reference:
c
c    Nick Trefethen,
c    Numerical Computation of the Schwarz-Christoffel Transformation,
c    SIAM Journal of Scientific and Statistical Computing,
c    Volume 1, 1980, pages 82-102.
c
      common /logcnt/ ncount, ncnt1

      ncount = 0
      ncnt1 = 0

      write (6,1)
    1 format (' ------- log counter set to zero')

      return
      end
      subroutine count

c*********************************************************************72
c
cc COUNT prints the internal counter.
c
c  License:
c
c    The package may be freely copied and used, but reference to
c    the paper by Nick Trefethen should be given in any publications arising
c    from its use.
c 
c  Modified:
c
c    07 December 2007
c
c  Author:
c
c    Nick Trefethen
c
c  Reference:
c
c    Nick Trefethen,
c    Numerical Computation of the Schwarz-Christoffel Transformation,
c    SIAM Journal of Scientific and Statistical Computing,
c    Volume 1, 1980, pages 82-102.
c
      common /logcnt/ ncount,ncnt1

      ncdiff = ncount - ncnt1

      write (6,2) ncdiff,ncount
    2 format (' ------- no. logs: since last count',i7,',   total',i8)

      ncnt1 = ncount

      return
      end
