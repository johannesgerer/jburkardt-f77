      subroutine cheb_r4 (n,a,b,fnu,alpha,beta,s,ierr,s0,s1,s2)

c*********************************************************************72
c
cc CHEB_R4 generates recursion coefficients ALPHA and BETA.
c
c Given a set of polynomials  p(0),p(1),...,p(2*n-1)  satisfying
c
c        p(k+1)(x)=(x-a(k))*p(k)(x)-b(k)*p(k-1)(x),
c                        k=0,1,...,2*n-2,
c
c        p(-1)(x)=0,  p(0)(x)=1,
c
c and associated modified moments
c
c           fnu(k)=integral of p(k)(x)*dlambda(x),
c                        k=0,1,...,2*n-1,
c
c this subroutine uses the modified Chebyshev algorithm to generate the
c recursion coefficients  alpha(k),beta(k), k=0,1,...,n-1, for the
c polynomials  pi(k)  orthogonal with respect to the integration
c measure  dlambda(x), i.e.,
c
c        pi(k+1)(x)=(x-alpha(k))*pi(k)(x)-beta(k)*pi(k-1)(x),
c                          k=0,1,...,n-1,
c
c        pi(-1)(x)=0,  pi(0)(x)=1.
c
c  Reference:
c
c    Walter Gautschi,
c    On Generating Orthogonal Polynomials,
c    SIAM Journal on Scientific and Statistical Computing,
c    Volume 3, Number 3, 1982, pages 289-317.
c
c     Input:    n - - the number of recursion coefficients desired
c               a,b-- arrays of dimension 2*n-1 to be filled with the
c                     values of  a(k-1),b(k-1), k=1,2,...,2*n-1
c               fnu-- array of dimension  2*n  to be filled with the
c                     values of the modified moments  fnu(k-1), k=1,2,
c                     ...,2*n
c     Output:   alpha,beta-- arrays containing, respectively, the
c                     recursion coefficients  alpha(k-1),beta(k-1),
c                     k=1,2,...,n, where  beta(0)  is the total mass.
c               s - - array containing the normalization factors
c                     s(k)=integral [pi(k)(x)]**2 dlambda(x), k=0,1,
c                     2,...,n-1.
c               ierr- an error flag, equal to  0  on normal return,
c                     equal to  1  if  abs(fnu(0))  is less than the
c                     machine zero, equal to  2  if  n  is out of range,
c                     equal to  -k  if  s(k), k=0,1,2,...,n-1, is about
c                     to underflow, and equal to  +k  if it is about to
c                     overflow.
c
c The arrays  s0,s1,s2  are needed for working space.
c
c On machines with limited exponent range, the occurrence of underflow
c [overflow] in the computation of the  alpha's  and  beta's  can often
c be avoided by multiplying all modified moments by a sufficiently large
c [small] scaling factor and dividing the new  beta(0)  by the same
c scaling factor.
c
c The routine uses the function subroutine  r1mach.
c
      dimension a(*),b(*),fnu(*),alpha(n),beta(n),s(n),s0(*),s1(*),
     *s2(*)
c
c The arrays  a,b  are assumed to have dimension  2*n-1, the arrays
c fnu,s0,s1,s2  dimension  2*n.
c
      nd=2*n
      tiny=10.0*r1mach(1)
      huge=0.1*r1mach(2)
      ierr=0
      if(abs(fnu(1)).lt.tiny) then
        ierr=1
        return
      end if
      if(n.lt.1) then
        ierr=2
        return
      end if
c
c Initialization
c
      alpha(1)=a(1)+fnu(2)/fnu(1)
      beta(1)=fnu(1)
      if(n.eq.1) return
      s(1)=fnu(1)
      do l=1,nd
        s0(l)=0.0
        s1(l)=fnu(l)
      end do
c
c Continuation
c
      do k=2,n
        lk=nd-k+1
        do l=k,lk
c
c The quantities  s2(l)  for l > k are auxiliary quantities which may
c be zero or may become so small as to underflow, without however
c causing any harm.
c
          s2(l)=s1(l+1)-(alpha(k-1)-a(l))*s1(l)-beta(k-1)*s0(l)
     *      +b(l)*s1(l-1)
          if(l.eq.k) s(k)=s2(k)
c
c Check impending underflow or overflow
c
          if(abs(s(k)).lt.tiny) then
            ierr=-(k-1)
            return
          else if(abs(s(k)).gt.huge) then
            ierr=k-1
            return
          end if
        end do
c
c Compute the alpha- and beta-coefficient
c
        alpha(k)=a(k)+(s2(k+1)/s2(k))-(s1(k)/s1(k-1))
        beta(k)=s2(k)/s1(k-1)
        do l=k,lk
          s0(l)=s1(l)
          s1(l)=s2(l)
        end do
      end do

      return
      end
      subroutine chri(n,iopt,a,b,x,y,hr,hi,alpha,beta,ierr)

c*********************************************************************72
c
cc CHRI_R4 implements the Christoffel or generalized Christoffel theorem.
c
c
c This subroutine implements the Christoffel or generalized Christoffel
c theorem. In all cases except  iopt=7, it uses nonlinear recurrence
c algorithms described in W. Gautschi,An algorithmic implementation
c of the generalized Christoffel theorem'', Numerical Integration
c (G. Haemmerlin, ed.), Birkhaeuser, Basel, 1982, pp. 89-106. The case
c iopt=7  incorporates a QR step with shift  x  in the manner of
c J. Kautsky and G.H. Golub, On the calculation of Jacobi matrices'',
c Linear Algebra Appl. 52/53, 1983, 439-455, using the algorithm of
c Eq. (67.11) on p. 567 in J.H. Wilkinson,The Algebraic Eigenvalue
c Problem'', Clarendon Press, Oxford, 1965. Given the recursion
c coefficients  a(k),b(k), k=0,1,...,n, for the (monic) orthogonal
c polynomials with respect to some measure  dlambda(t), it generates
c the recursion coefficients  alpha(k),beta(k), k=0,1,...,n-1, for the
c measure
c
c              (t-x)dlambda(t)               if  iopt=1
c              [(t-x)**2+y**2]dlambda(t)     if  iopt=2
c              (t**2+y**2)dlambda(t) with    if  iopt=3
c                dlambda(t) and supp(dlambda)
c                symmetric  with respect to
c                the origin
c              dlambda(t)/(t-x)              if  iopt=4
c              dlambda(t)/[(t-x)**2+y**2]    if  iopt=5
c              dlambda(t)/(t**2+y**2) with   if  iopt=6
c                dlambda(t) and supp(dlambda)
c                symmetric with respect to
c                the origin
c              [(t-x)**2]dlambda(t)          if  iopt=7
c
c
c      Input:   n  - - - the number of recurrence coefficients
c                        desired; type integer
c               iopt - - an integer selecting the desired weight
c                        distribution
c               a,b  - - arrays of dimension  n+1  containing the
c                        recursion coefficients a(k-1),b(k-1),k=1,2,
c                        ...,n+1, of the polynomials orthogonal with
c                        respect to the given measure  dlambda(t)
c               x,y  - - real parameters defining the linear and
c                        quadratic factors, or divisors, of  dlambda(t)
c               hr,hi  - the real and imaginary part, respectively, of
c                        the integral of dlambda(t)/(z-t), where z=x+iy;
c                        the parameter  hr  is used only if  iopt=4 or
c                        5, the parameter  hi  only if  iopt=5 or 6
c
c      Output:  alpha,beta - - arrays of dimension  n  containing the
c                         desired recursion coefficients  alpha(k-1),
c                         beta(k-1), k=1,2,...,n
c
c It is assumed that  n  is larger than or equal to 2. Otherwise, the
c routine exits immediately with the error flag  ierr  set equal to 1.
c If  iopt  is not between 1 and 7, the routine exits with  ierr=2.
c
c The routine uses the function subroutine  r1mach  to evaluate the
c constant  eps, which is used only if  iopt=7.
c
      dimension a(*),b(*),alpha(n),beta(n)
c
c The arrays  a,b  are assumed to have dimension  n+1.
c
      eps=5.*r1mach(3)
c
c The quantity  eps  is a constant slightly larger than the machine
c precision.
c
      ierr=0
      if(n.lt.2) then
        ierr=1
        return
      end if
c
c What follows implements Eq. (3.7) of W. Gautschi, op. cit.
c
      if (iopt.eq.1) then
        e=0.
        do 10 k=1,n
          q=a(k)-e-x
          beta(k)=q*e
          e=b(k+1)/q
          alpha(k)=x+q+e
   10   continue
c
c Set the first beta-coefficient as discussed in Section 5.1 of the
c companion paper.
c
        beta(1)=b(1)*(a(1)-x)
        return
c
c What follows implements Eq. (4.7) of W. Gautschi, op. cit.
c
      else if(iopt.eq.2) then
        s=x-a(1)
        t=y
        eio=0.
        do 20 k=1,n
          d=s*s+t*t
          er=-b(k+1)*s/d
          ei=b(k+1)*t/d
          s=x+er-a(k+1)
          t=y+ei
          alpha(k)=x+t*er/ei-s*ei/t
          beta(k)=t*eio*(1.+(er/ei)**2)
          eio=ei
   20   continue
c
c Set the first beta-coefficient.
c
        beta(1)=b(1)*(b(2)+(a(1)-x)**2+y*y)
        return
c
c What follows implements Eq. (4.8) of W. Gautschi, op. cit.
c
      else if(iopt.eq.3) then
        t=y
        eio=0.
        do 30 k=1,n
          ei=b(k+1)/t
          t=y+ei
          alpha(k)=0.
          beta(k)=t*eio
          eio=ei
   30   continue
c
c Set the first beta-coefficient.
c
        beta(1)=b(1)*(b(2)+y*y)
        return
c
c What follows implements Eqs. (5.1),(5.2) of W. Gautschi, op. cit.
c
      else if(iopt.eq.4) then
        alpha(1)=x-b(1)/hr
        beta(1)=-hr
        q=-b(1)/hr
        do 40 k=2,n
          e=a(k-1)-x-q
          beta(k)=q*e
          q=b(k)/e
          alpha(k)=q+e+x
   40   continue
        return
c
c What follows implements Eq. (5.8) of W. Gautschi, op. cit.
c
      else if(iopt.eq.5) then
        nm1=n-1
        d=hr*hr+hi*hi
        eroo=a(1)-x+b(1)*hr/d
        eioo=-b(1)*hi/d-y
        alpha(1)=x+hr*y/hi
        beta(1)=-hi/y
        alpha(2)=x-b(1)*hi*eroo/(d*eioo)+hr*eioo/hi
        beta(2)=y*eioo*(1.+(hr/hi)**2)
        if(n.eq.2) return
        so=b(2)/(eroo**2+eioo**2)
        ero=a(2)-x-so*eroo
        eio=so*eioo-y
        alpha(3)=x+eroo*eio/eioo+so*eioo*ero/eio
        beta(3)=-b(1)*hi*eio*(1.+(eroo/eioo)**2)/d
        if(n.eq.3) return
        do 50 k=3,nm1
          s=b(k)/(ero**2+eio**2)
          er=a(k)-x-s*ero
          ei=s*eio-y
          alpha(k+1)=x+ero*ei/eio+s*eio*er/ei
          beta(k+1)=so*eioo*ei*(1.+(ero/eio)**2)
          eroo=ero
          eioo=eio
          ero=er
          eio=ei
          so=s
   50   continue
        return
c
c What follows implements Eq. (5.9) of W. Gautschi, op. cit.
c
      else if(iopt.eq.6) then
        nm1=n-1
        eoo=-b(1)/hi-y
        eo=b(2)/eoo-y
        alpha(1)=0.
        beta(1)=-hi/y
        alpha(2)=0.
        beta(2)=y*eoo
        if(n.eq.2) return
        alpha(3)=0.
        beta(3)=-b(1)*eo/hi
        if(n.eq.3) return
        do 60 k=3,nm1
          e=b(k)/eo-y
          beta(k+1)=b(k-1)*e/eoo
          alpha(k+1)=0.
          eoo=eo
          eo=e
   60   continue
        return
c
c What follows implements a QR step with shift  x.
c
      else if(iopt.eq.7) then
        u=0.
        c=1.
        c0=0.
        do 70 k=1,n
          gamma=a(k)-x-u
          cm1=c0
          c0=c
          if(abs(c0).gt.eps) then
            p2=(gamma**2)/c0
          else
            p2=cm1*b(k)
          end if
          if(k.gt.1) beta(k)=s*(p2+b(k+1))
          s=b(k+1)/(p2+b(k+1))
          c=p2/(p2+b(k+1))
          u=s*(gamma+a(k+1)-x)
          alpha(k)=gamma+u+x
   70   continue
        beta(1)=b(1)*(b(2)+(x-a(1))**2)
        return
      else
        ierr=2
        return
      end if
      end
      subroutine cheb_r8 ( n, da, db, dnu, dalpha, dbeta, ds, iderr,
     &  ds0, ds1, ds2 )

c*********************************************************************72
c
cc CHEB_R8 generates recursion coefficients ALPHA and BETA.
c
c  Reference:
c
c    Walter Gautschi,
c    On Generating Orthogonal Polynomials,
c    SIAM Journal on Scientific and Statistical Computing,
c    Volume 3, Number 3, 1982, pages 289-317.
c
      double precision da,db,dnu,dalpha,dbeta,ds,ds0,ds1,ds2,dtiny,
     *d1mach,dhuge
      dimension da(*),db(*),dnu(*),dalpha(n),dbeta(n),ds(n),
     *ds0(*),ds1(*),ds2(*)
c
c The arrays  da,db  are assumed to have dimension  2*n-1, the arrays
c dnu,ds0,ds1,ds2  dimension  2*n.
c
      nd=2*n
      dtiny=10.d0*d1mach(1)
      dhuge=.1d0*d1mach(2)
      iderr=0
      if(dabs(dnu(1)).lt.dtiny) then
        iderr=1
        return
      end if
      if(n.lt.1) then
        iderr=2
        return
      end if
      dalpha(1)=da(1)+dnu(2)/dnu(1)
      dbeta(1)=dnu(1)
      if(n.eq.1) return
      ds(1)=dnu(1)
      do 10 l=1,nd
        ds0(l)=0.d0
        ds1(l)=dnu(l)
   10 continue
      do 40 k=2,n
        lk=nd-k+1
        do 20 l=k,lk
          ds2(l)=ds1(l+1)-(dalpha(k-1)-da(l))*ds1(l)-dbeta(k-1)*ds0(l)
     *      +db(l)*ds1(l-1)
        if(l.eq.k) ds(k)=ds2(k)
   20   continue
        if(dabs(ds(k)).lt.dtiny) then
          iderr=-(k-1)
          return
        else if(dabs(ds(k)).gt.dhuge) then
          iderr=k-1
          return
        end if
        dalpha(k)=da(k)+(ds2(k+1)/ds2(k))-(ds1(k)/ds1(k-1))
        dbeta(k)=ds2(k)/ds1(k-1)
        do 30 l=k,lk
          ds0(l)=ds1(l)
          ds1(l)=ds2(l)
   30   continue
   40 continue
      return
      end
      subroutine chri_r8(n,iopt,da,db,dx,dy,dhr,dhi,dalpha,dbeta,ierr)

c*********************************************************************72
c
cc CHRI_R8 implements the Christoffel or generalized Christoffel theorem.
c
      double precision da,db,dx,dy,dhr,dhi,dalpha,dbeta,deps,d1mach,
     *de,dq,ds,dt,deio,dd,der,dei,deroo,deioo,dso,dero,deoo,deo,du,dc,
     *dc0,dgam,dcm1,dp2
      dimension da(*),db(*),dalpha(n),dbeta(n)
c
c The arrays  da,db  are assumed to have dimension  n+1.
c
      deps=5.d0*d1mach(3)
      ierr=0
      if(n.lt.2) then
        ierr=1
        return
      end if
      if(iopt.eq.1) then
        de=0.d0
        do 10 k=1,n
          dq=da(k)-de-dx
          dbeta(k)=dq*de
          de=db(k+1)/dq
          dalpha(k)=dx+dq+de
   10   continue
        dbeta(1)=db(1)*(da(1)-dx)
        return
      else if(iopt.eq.2) then
        ds=dx-da(1)
        dt=dy
        deio=0.d0
        do 20 k=1,n
          dd=ds*ds+dt*dt
          der=-db(k+1)*ds/dd
          dei=db(k+1)*dt/dd
          ds=dx+der-da(k+1)
          dt=dy+dei
          dalpha(k)=dx+dt*der/dei-ds*dei/dt
          dbeta(k)=dt*deio*(1.d0+(der/dei)**2)
          deio=dei
   20   continue
        dbeta(1)=db(1)*(db(2)+(da(1)-dx)**2+dy*dy)
        return
      else if(iopt.eq.3) then
        dt=dy
        deio=0.d0
        do 30 k=1,n
          dei=db(k+1)/dt
          dt=dy+dei
          dalpha(k)=0.d0
          dbeta(k)=dt*deio
          deio=dei
   30   continue
        dbeta(1)=db(1)*(db(2)+dy*dy)
        return
      else if(iopt.eq.4) then
        dalpha(1)=dx-db(1)/dhr
        dbeta(1)=-dhr
        dq=-db(1)/dhr
        do 40 k=2,n
          de=da(k-1)-dx-dq
          dbeta(k)=dq*de
          dq=db(k)/de
          dalpha(k)=dq+de+dx
   40   continue
        return
      else if(iopt.eq.5) then
        nm1=n-1
        dd=dhr*dhr+dhi*dhi
        deroo=da(1)-dx+db(1)*dhr/dd
        deioo=-db(1)*dhi/dd-dy
        dalpha(1)=dx+dhr*dy/dhi
        dbeta(1)=-dhi/dy
        dalpha(2)=dx-db(1)*dhi*deroo/(dd*deioo)+dhr*deioo/dhi
        dbeta(2)=dy*deioo*(1.d0+(dhr/dhi)**2)
        if(n.eq.2) return
        dso=db(2)/(deroo**2+deioo**2)
        dero=da(2)-dx-dso*deroo
        deio=dso*deioo-dy
        dalpha(3)=dx+deroo*deio/deioo+dso*deioo*dero/deio
        dbeta(3)=-db(1)*dhi*deio*(1.d0+(deroo/deioo)**2)/dd
        if(n.eq.3) return
        do 50 k=3,nm1
          ds=db(k)/(dero**2+deio**2)
          der=da(k)-dx-ds*dero
          dei=ds*deio-dy
          dalpha(k+1)=dx+dero*dei/deio+ds*deio*der/dei
          dbeta(k+1)=dso*deioo*dei*(1.d0+(dero/deio)**2)
          deroo=dero
          deioo=deio
          dero=der
          deio=dei
          dso=ds
   50   continue
        return
      else if(iopt.eq.6) then
        nm1=n-1
        deoo=-db(1)/dhi-dy
        deo=db(2)/deoo-dy
        dalpha(1)=0.d0
        dbeta(1)=-dhi/dy
        dalpha(2)=0.d0
        dbeta(2)=dy*deoo
        if(n.eq.2) return
        dalpha(3)=0.d0
        dbeta(3)=-db(1)*deo/dhi
        if(n.eq.3) return
        do 60 k=3,nm1
          de=db(k)/deo-dy
          dbeta(k+1)=db(k-1)*de/deoo
          dalpha(k+1)=0.d0
          deoo=deo
          deo=de
   60   continue
        return
      else if(iopt.eq.7) then
        du=0.d0
        dc=1.d0
        dc0=0.d0
        do 70 k=1,n
          dgam=da(k)-dx-du
          dcm1=dc0
          dc0=dc
          if(dabs(dc0).gt.deps) then
            dp2=(dgam**2)/dc0
          else
            dp2=dcm1*db(k)
          end if
          if(k.gt.1) dbeta(k)=ds*(dp2+db(k+1))
          ds=db(k+1)/(dp2+db(k+1))
          dc=dp2/(dp2+db(k+1))
          du=ds*(dgam+da(k+1)-dx)
          dalpha(k)=dgam+du+dx
   70   continue
        dbeta(1)=db(1)*(db(2)+(dx-da(1))**2)
        return
      else
        ierr=2
        return
      end if
      end
      subroutine dgauss(n,dalpha,dbeta,deps,dzero,dweigh,ierr,de)

c*********************************************************************72
c
cc DGAUSS is a double-precision version of the routine  gauss.
c
      double precision dalpha,dbeta,deps,dzero,dweigh,de,dp,dg,dr,
     *ds,dc,df,db
      dimension dalpha(n),dbeta(n),dzero(n),dweigh(n),de(n)
      if(n.lt.1) then
        ierr=-1
        return
      end if
      ierr=0
      dzero(1)=dalpha(1)
      if(dbeta(1).lt.0.d0) then
        ierr=-2
        return
      end if
      dweigh(1)=dbeta(1)
      if (n.eq.1) return
      dweigh(1)=1.d0
      de(n)=0.d0
      do 100 k=2,n
        dzero(k)=dalpha(k)
        if(dbeta(k).lt.0.d0) then
          ierr=-2
          return
        end if
        de(k-1)=dsqrt(dbeta(k))
        dweigh(k)=0.d0
  100 continue
      do 240 l=1,n
        j=0
  105   do 110 m=l,n
          if(m.eq.n) goto 120
          if(dabs(de(m)).le.deps*(dabs(dzero(m))+dabs(dzero(m+1))))
     *      goto 120
  110   continue
  120   dp=dzero(l)
        if(m.eq.l) goto 240
        if(j.eq.30) goto 400
        j=j+1
        dg=(dzero(l+1)-dp)/(2.d0*de(l))
        dr=dsqrt(dg*dg+1.d0)
        dg=dzero(m)-dp+de(l)/(dg+dsign(dr,dg))
        ds=1.d0
        dc=1.d0
        dp=0.d0
        mml=m-l
        do 200 ii=1,mml
          i=m-ii
          df=ds*de(i)
          db=dc*de(i)
          if(dabs(df).lt.dabs(dg)) goto 150
          dc=dg/df
          dr=dsqrt(dc*dc+1.d0)
          de(i+1)=df*dr
          ds=1.d0/dr
          dc=dc*ds
          goto 160
  150     ds=df/dg
          dr=dsqrt(ds*ds+1.d0)
          de(i+1)=dg*dr
          dc=1.d0/dr
          ds=ds*dc
  160     dg=dzero(i+1)-dp
          dr=(dzero(i)-dg)*ds+2.d0*dc*db
          dp=ds*dr
          dzero(i+1)=dg+dp
          dg=dc*dr-db
          df=dweigh(i+1)
          dweigh(i+1)=ds*dweigh(i)+dc*df
          dweigh(i)=dc*dweigh(i)-ds*df
  200   continue
        dzero(l)=dzero(l)-dp
        de(l)=dg
        de(m)=0.d0
        goto 105
  240 continue
      do 300 ii=2,n
        i=ii-1
        k=i
        dp=dzero(i)
        do 260 j=ii,n
          if(dzero(j).ge.dp) goto 260
          k=j
          dp=dzero(j)
  260   continue
        if(k.eq.i) goto 300
        dzero(k)=dzero(i)
        dzero(i)=dp
        dp=dweigh(i)
        dweigh(i)=dweigh(k)
        dweigh(k)=dp
  300 continue
      do 310 k=1,n
        dweigh(k)=dbeta(1)*dweigh(k)*dweigh(k)
  310 continue
      return
  400 ierr=l
      return
      end
      subroutine dgchri(n,iopt,nu0,numax,deps,da,db,dx,dy,dalpha,dbeta,
     *nu,ierr,ierrc,dnu,drhor,drhoi,droldr,droldi,ds,ds0,ds1,ds2)

c*********************************************************************72
c
cc DGCHRI is a double-precision version of the routine  gchri.
c
      double precision deps,da(numax),db(numax),dx,dy,dalpha(n),
     *dbeta(n),dnu(*),drhor(*),drhoi(*),droldr(*),droldi(*),
     *ds(n),ds0(*),ds1(*),ds2(*)
c
c The arrays  dnu,drhor,drhoi,droldr,droldi,ds0,ds1,ds2  are assumed
c to have dimension  2*n.
c
      if(n.lt.1) then
        ierr=-1
        return
      end if
      ierr=0
      nd=2*n
      ndm1=nd-1
      if(iopt.eq.1) then
        call dknum(ndm1,nu0,numax,dx,dy,deps,da,db,drhor,drhoi,nu,
     *    ierr,droldr,droldi)
        do 10 k=1,nd
          dnu(k)=-drhor(k)
   10   continue
        call cheb_r8 ( n,da,db,dnu,dalpha,dbeta,ds,ierrc,ds0,ds1,ds2)
        return
      else if(iopt.eq.2) then
        dy=dabs(dy)
        call dknum(ndm1,nu0,numax,dx,dy,deps,da,db,drhor,drhoi,nu,
     *    ierr,droldr,droldi)
        do 20 k=1,nd
          dnu(k)=-drhoi(k)/dy
   20   continue
        call cheb_r8 (n,da,db,dnu,dalpha,dbeta,ds,ierrc,ds0,ds1,ds2)
        return
      else
        ierr=1
        return
      end if
      end
      subroutine dkern(n,nu0,numax,dx,dy,deps,da,db,dkerr,dkeri,
     *  nu,ierr,droldr,droldi)

c*********************************************************************72
c
cc DKERN is a double-precision version of the routine  kern.
c
      double precision dx,dy,deps,da(numax),db(numax),dkerr(*),
     *  dkeri(*),droldr(*),droldi(*),dp0r,dp0i,dpr,dpi,dpm1r,
     *  dpm1i,dden,dt
c
c The arrays  dkerr,dkeri,droldr,droldi  are assumed to have
c dimension  n+1.
c
      call dknum(n,nu0,numax,dx,dy,deps,da,db,dkerr,dkeri,nu,ierr,
     *  droldr,droldi)
      if(ierr.ne.0) return
      dp0r=0.d0
      dp0i=0.d0
      dpr=1.d0
      dpi=0.d0
      do 10 k=1,n
        dpm1r=dp0r
        dpm1i=dp0i
        dp0r=dpr
        dp0i=dpi
        dpr=(dx-da(k))*dp0r-dy*dp0i-db(k)*dpm1r
        dpi=(dx-da(k))*dp0i+dy*dp0r-db(k)*dpm1i
        dden=dpr**2+dpi**2
        dt=(dkerr(k+1)*dpr+dkeri(k+1)*dpi)/dden
        dkeri(k+1)=(dkeri(k+1)*dpr-dkerr(k+1)*dpi)/dden
        dkerr(k+1)=dt
   10 continue
      return
      end
      subroutine dknum(n,nu0,numax,dx,dy,deps,da,db,drhor,drhoi,nu,
     *ierr,droldr,droldi)

c*********************************************************************72
c
cc DKNUM is a double-precision version of the routine  knum.
c
      double precision dx,dy,deps,da(numax),db(numax),drhor(*),
     *drhoi(*),droldr(*),droldi(*),drr,dri,dden,dt
c
c The arrays  drhor,drhoi,droldr,droldi  are assumed to have
c dimension  n+1.
c
      ierr=0
      np1=n+1
      if(nu0.gt.numax) then
        ierr=nu0
        return
      end if
      if(nu0.lt.np1) nu0=np1
      nu=nu0-5
      do 10 k=1,np1
        drhor(k)=0.d0
        drhoi(k)=0.d0
   10 continue
   20 nu=nu+5
      if(nu.gt.numax) then
        ierr=numax
        goto 60
      end if
      do 30 k=1,np1
        droldr(k)=drhor(k)
        droldi(k)=drhoi(k)
   30 continue
      drr=0.d0
      dri=0.d0
      do 40 j=1,nu
        j1=nu-j+1
        dden=(dx-da(j1)-drr)**2+(dy-dri)**2
        drr=db(j1)*(dx-da(j1)-drr)/dden
        dri=-db(j1)*(dy-dri)/dden
        if(j1.le.np1) then
          drhor(j1)=drr
          drhoi(j1)=dri
        end if
   40 continue
      do 50 k=1,np1
C Following statement replaced -- authors remark
c
c       if((drhor(k)-droldr(k))**2+(drhoi(k)-droldi(k))**2.gt.
c    *    deps*(drhor(k)**2+drhoi(k)**2)) goto 20
        if((drhor(k)-droldr(k))**2+(drhoi(k)-droldi(k))**2.gt.
     *    (deps**2)*(drhor(k)**2+drhoi(k)**2)) goto 20
   50 continue
   60 if(n.eq.0) return
      do 70 k=2,np1
        dt=drhor(k)*drhor(k-1)-drhoi(k)*drhoi(k-1)
        drhoi(k)=drhor(k)*drhoi(k-1)+drhoi(k)*drhor(k-1)
        drhor(k)=dt
   70 continue
      return
      end
      subroutine dlancz(n,ncap,dx,dw,dalpha,dbeta,ierr,dp0,dp1)

c*********************************************************************72
c
cc DLANCZ is a double-precision version of the routine  lancz.
c
      double precision dx(ncap),dw(ncap),dalpha(n),dbeta(n),
     *dp0(ncap),dp1(ncap),dpi,dgam,dsig,dt,dxlam,drho,dtmp,
     *dtsig,dtk
      if(n.le.0 .or. n.gt.ncap) then
        ierr=1
        return
      else
        ierr=0
      end if
      do 10 i=1,ncap
        dp0(i)=dx(i)
        dp1(i)=0.d0
   10 continue
      dp1(1)=dw(1)
      do 30 i=1,ncap-1
        dpi=dw(i+1)
        dgam=1.d0
        dsig=0.d0
        dt=0.d0
        dxlam=dx(i+1)
        do 20 k=1,i+1
          drho=dp1(k)+dpi
          dtmp=dgam*drho
          dtsig=dsig
          if(drho.le.0.d0) then
            dgam=1.d0
            dsig=0.d0
          else
            dgam=dp1(k)/drho
            dsig=dpi/drho
          end if
          dtk=dsig*(dp0(k)-dxlam)-dgam*dt
          dp0(k)=dp0(k)-(dtk-dt)
          dt=dtk
          if(dsig.le.0.d0) then
            dpi=dtsig*dp1(k)
          else
            dpi=(dt**2)/dsig
          end if
          dtsig=dsig
          dp1(k)=dtmp
   20   continue
   30 continue
      do 40 k=1,n
        dalpha(k)=dp0(k)
        dbeta(k)=dp1(k)
   40 continue
      return
      end
      subroutine dlob(n,dalpha,dbeta,dleft,dright,dzero,dweigh,
     *ierr,de,da,db)

c*********************************************************************72
c
cc DLOB is a double-precision version of the routine  lob.
c
      double precision dleft,dright,depsma,dp0l,dp0r,dp1l,dp1r,dpm1l,
     *dpm1r,ddet,dalpha(*),dbeta(*),dzero(*),dweigh(*),de(*),da(*),
     *db(*),d1mach
c
c The arrays  dalpha,dbeta,dzero,dweigh,de,da,db  are assumed to have
c dimension  n+2.
c
      depsma=d1mach(3)
c
c depsma is the machine double precision.
c
      np1=n+1
      np2=n+2
      do 10 k=1,np2
        da(k)=dalpha(k)
        db(k)=dbeta(k)
   10 continue
      dp0l=0.d0
      dp0r=0.d0
      dp1l=1.d0
      dp1r=1.d0
      do 20 k=1,np1
        dpm1l=dp0l
        dp0l=dp1l
        dpm1r=dp0r
        dp0r=dp1r
        dp1l=(dleft-da(k))*dp0l-db(k)*dpm1l
        dp1r=(dright-da(k))*dp0r-db(k)*dpm1r
   20 continue
      ddet=dp1l*dp0r-dp1r*dp0l
      da(np2)=(dleft*dp1l*dp0r-dright*dp1r*dp0l)/ddet
      db(np2)=(dright-dleft)*dp1l*dp1r/ddet
      call dgauss(np2,da,db,depsma,dzero,dweigh,ierr,de)
      return
      end
      subroutine dmcdis(n,ncapm,mc,mp,dxp,dyp,dquad,deps,iq,idelta,
     *irout,finld,finrd,dendl,dendr,dxfer,dwfer,dalpha,dbeta,ncap,
     *kount,ierrd,ied,dbe,dx,dw,dxm,dwm,dp0,dp1,dp2)

c*********************************************************************72
c
cc DMCDIS is a double-precision version of the routine  mcdis.
c
      double precision dxp,dyp,deps,dendl,dendr,dxfer,dwfer,dalpha,
     *dbeta,dbe,dx,dw,dxm,dwm,dp0,dp1,dp2
      dimension dxp(*),dyp(*),dendl(mc),dendr(mc),dxfer(ncapm),
     *dwfer(ncapm),dalpha(n),dbeta(n),dbe(n),dx(ncapm),dw(ncapm),
     *dxm(*),dwm(*),dp0(*),dp1(*),dp2(*)
      logical finld,finrd
c
c The arrays  dxp,dyp  are assumed to have dimension  mp  if mp > 0,
c the arrays  dxm,dwm,dp0,dp1,dp2  dimension  mc*ncapm+mp.
c
      if(idelta.le.0) idelta=1
      if(n.lt.1) then
        ierrd=-1
        return
      end if
      incap=1
      kount=-1
      ierr=0
      do 10 k=1,n
        dbeta(k)=0.d0
   10 continue
      ncap=(2*n-1)/idelta
   20 do 30 k=1,n
        dbe(k)=dbeta(k)
   30 continue
      kount=kount+1
      if(kount.gt.1) incap=2**(kount/5)*n
      ncap=ncap+incap
      if(ncap.gt.ncapm) then
        ierrd=ncapm
        return
      end if
      mtncap=mc*ncap
      do 50 i=1,mc
        im1tn=(i-1)*ncap
        if(iq.eq.1) then
          call dquad(ncap,dx,dw,i,ierr)
        else
          call dqgp(ncap,dx,dw,i,ierr,mc,finld,finrd,dendl,dendr,dxfer,
     *      dwfer)
        end if
        if(ierr.ne.0) then
          ierrd=i
          return
        end if
        do 40 k=1,ncap
          dxm(im1tn+k)=dx(k)
          dwm(im1tn+k)=dw(k)
   40   continue
   50 continue
      if(mp.ne.0) then
        do 60 k=1,mp
          dxm(mtncap+k)=dxp(k)
          dwm(mtncap+k)=dyp(k)
   60   continue
      end if
      if(irout.eq.1) then
        call dsti(n,mtncap+mp,dxm,dwm,dalpha,dbeta,ied,dp0,dp1,dp2)
      else
        call dlancz(n,mtncap+mp,dxm,dwm,dalpha,dbeta,ied,dp0,dp1)
      end if
      do 70 k=1,n
        if(dabs(dbeta(k)-dbe(k)).gt.deps*dabs(dbeta(k))) goto 20
   70 continue
      return
      end
      subroutine dmcheb(n,ncapm,mcd,mp,dxp,dyp,dquad,deps,iq,
     *idelta,finld,finrd,dendl,dendr,dxfer,dwfer,da,db,dnu,dalpha,
     *dbeta,ncap,kount,ierrd,dbe,dx,dw,dxm,dwm,ds,ds0,ds1,ds2)

c*********************************************************************72
c
cc DMCHEB is a double-precision version of the routine  mccheb.
c
      double precision dxp,dyp,deps,dendl,dendr,dxfer,dwfer,da,db,
     *dnu,dalpha,dbeta,dbe,dx,dw,dxm,dwm,ds,ds0,ds1,ds2,dsum,dp1,
     *dp,dpm1
      dimension dxp(*),dyp(*),dendl(mcd),dendr(mcd),dxfer(ncapm),
     *dwfer(ncapm),da(*),db(*),dnu(*),dalpha(n),dbeta(n),dbe(n),
     *dx(ncapm),dw(ncapm),dxm(*),dwm(*),ds(n),ds0(*),ds1(*),ds2(*)
      logical finld,finrd
c
c The arrays  dxp,dyp  are assumed to have dimension  mp  if mp > 0,
c the arrays  da,db  dimension 2*n-1, the arrays  dnu,ds0,ds1,ds2
c dimension  2*n, and the arrays  dxm,dwm  dimension  mc*ncapm+mp.
c
      nd=2*n
      if(idelta.le.0) idelta=1
      if(n.lt.1) then
        ierrd=-1
        return
      end if
      incap=1
      kount=-1
      ierrd=0
      do 10 k=1,n
        dbeta(k)=0.d0
   10 continue
      ncap=(nd-1)/idelta
   20 do 30 k=1,n
        dbe(k)=dbeta(k)
   30 continue
      kount=kount+1
      if(kount.gt.1) incap=2**(kount/5)*n
      ncap=ncap+incap
      if(ncap.gt.ncapm) then
        ierrd=ncapm
        return
      end if
      mtncap=mcd*ncap
      do 50 i=1,mcd
        im1tn=(i-1)*ncap
        if(iq.eq.1) then
          call dquad(ncap,dx,dw,i,ierr)
        else
          call dqgp(ncap,dx,dw,i,ierr,mcd,finld,finrd,dendl,dendr,
     *      dxfer,dwfer)
        end if
        if(ierr.ne.0) then
          ierrd=i
          return
        end if
        do 40 k=1,ncap
          dxm(im1tn+k)=dx(k)
          dwm(im1tn+k)=dw(k)
   40   continue
   50 continue
      if(mp.ne.0) then
        do 60 k=1,mp
          dxm(mtncap+k)=dxp(k)
          dwm(mtncap+k)=dyp(k)
   60   continue
      end if
      mtnpmp=mtncap+mp
      do 90 k=1,nd
        km1=k-1
        dsum=0.d0
        do 80 i=1,mtnpmp
          dp1=0.d0
          dp=1.d0
          if(k.gt.1) then
            do 70 l=1,km1
              dpm1=dp1
              dp1=dp
              dp=(dxm(i)-da(l))*dp1-db(l)*dpm1
   70       continue
          end if
          dsum=dsum+dwm(i)*dp
   80   continue
        dnu(k)=dsum
   90 continue
      call cheb_r8 (n,da,db,dnu,dalpha,dbeta,ds,ierr,ds0,ds1,ds2)
      do 100 k=1,n
        if(dabs(dbeta(k)-dbe(k)).gt.deps*dabs(dbeta(k))) goto 20
  100 continue
      return
      end
      subroutine dqgp(n,dx,dw,i,ierr,mcd,finld,finrd,dendl,dendr,
     *  dxfer,dwfer)

c*********************************************************************72
c
cc DQGP is a double-precision version of the routine  qgp. 
c
c  The user has to supply the routine
c
c              double precision function dwf(dx,i),
c
c which evaluates the weight function in double precision at the
c point  dx  on the i-th component interval.
c
      double precision dx,dw,dendl,dendr,dxfer,dwfer,dphi,dphi1,dwf
      dimension dx(n),dw(n),dendl(mcd),dendr(mcd),dxfer(*),dwfer(*)
      logical finld,finrd
c
c The arrays  dxfer,dwfer  are dimensioned in the routine  dmcdis.
c
      ierr=0
      if(i.eq.1) call dfejer(n,dxfer,dwfer)
      if(i.gt.1 .and. i.lt.mcd) goto 60
      if(mcd.eq.1) then
        if(finld.and.finrd) goto 60
        if(finld) goto 20
        if(finrd) goto 40
        do 10 k=1,n
          call dsymtr(dxfer(k),dphi,dphi1)
          dx(k)=dphi
          dw(k)=dwfer(k)*dwf(dphi,i)*dphi1
   10   continue
        return
      else
        if((i.eq.1.and.finld).or.(i.eq.mcd.and.finrd)) goto 60
        if(i.eq.1) goto 40
      end if
   20 do 30 k=1,n
        call dtr(dxfer(k),dphi,dphi1)
        dx(k)=dendl(mcd)+dphi
        dw(k)=dwfer(k)*dwf(dx(k),mcd)*dphi1
   30 continue
      return
   40 do 50 k=1,n
        call dtr(-dxfer(k),dphi,dphi1)
        dx(k)=dendr(1)-dphi
        dw(k)=dwfer(k)*dwf(dx(k),1)*dphi1
   50 continue
      return
   60 do 70 k=1,n
        dx(k)=.5d0*((dendr(i)-dendl(i))*dxfer(k)+dendr(i)+dendl(i))
        dw(k)=.5d0*(dendr(i)-dendl(i))*dwfer(k)*dwf(dx(k),i)
   70 continue
      return
      end
      subroutine dsymtr(dt,dphi,dphi1)

c*********************************************************************72
c
cc DSYMTR is a double-precision version of  symtr.
c
      double precision dt,dphi,dphi1,dt2
      dt2=dt*dt
      dphi=dt/(1.-dt2)
      dphi1=(dt2+1.d0)/(dt2-1.d0)**2
      return
      end
      subroutine dtr(dt,dphi,dphi1)

c*********************************************************************72
c
cc DTR is a double-precision version of  tr.
c
      double precision dt,dphi,dphi1
      dphi=(1.d0+dt)/(1.d0-dt)
      dphi1=2.d0/(dt-1.d0)**2
      return
      end
      subroutine dfejer(n,dx,dw)

c*********************************************************************72
c
cc DFEJER is a double-precision version of  fejer.
c
      double precision dx,dw,dpi,dn,dc1,dc0,dt,dsum,dc2
      dimension dx(n),dw(n)
      dpi=4.d0*datan(1.d0)
      nh=n/2
      np1h=(n+1)/2
      dn=dble(n)
      do 10 k=1,nh
        dx(n+1-k)=dcos(.5d0*dble(2*k-1)*dpi/dn)
        dx(k)=-dx(n+1-k)
   10 continue
      if(2*nh.ne.n) dx(np1h)=0.d0
      do 30 k=1,np1h
        dc1=1.d0
        dc0=2.d0*dx(k)*dx(k)-1.d0
        dt=2.d0*dc0
        dsum=dc0/3.d0
        do 20 m=2,nh
          dc2=dc1
          dc1=dc0
          dc0=dt*dc1-dc2
          dsum=dsum+dc0/dble(4*m*m-1)
   20   continue
        dw(k)=2.d0*(1.d0-2.d0*dsum)/dn
        dw(n+1-k)=dw(k)
   30 continue
      return
      end
      subroutine dradau(n,dalpha,dbeta,dend,dzero,dweigh,ierr,de,
     *da,db)

c*********************************************************************72
c
c This is a double-precision version of the routine  radau.
c
      double precision dend,depsma,dp0,dp1,dpm1,dalpha(*),dbeta(*),
     *dzero(*),dweigh(*),de(*),da(*),db(*),d1mach
c
c The arrays  dalpha,dbeta,dzero,dweigh,de,da,db  are assumed to have
c dimension  n+1.
c
      depsma=d1mach(3)
c
c depsma is the machine double precision.
c
      np1=n+1
      do 10 k=1,np1
        da(k)=dalpha(k)
        db(k)=dbeta(k)
   10 continue
      dp0=0.d0
      dp1=1.d0
      do 20 k=1,n
        dpm1=dp0
        dp0=dp1
        dp1=(dend-da(k))*dp0-db(k)*dpm1
   20 continue
      da(np1)=dend-db(np1)*dp0/dp1
      call dgauss(np1,da,db,depsma,dzero,dweigh,ierr,de)
      return
      end
      subroutine drecur(n,ipoly,dal,dbe,da,db,iderr)

c*********************************************************************72
c
c This is a double-precision version of the routine  recur.
c
      external dgamma
      double precision dal,dbe,da,db,dlmach,d1mach,dkm1,dalpbe,dt,
     *dlga,dal2,dbe2,dgamma
      dimension da(n),db(n)
      if(n.lt.1) then
        iderr=3
        return
      end if
      dlmach=dlog(d1mach(2))
      iderr=0
      do 10 k=1,n
        da(k)=0.d0
   10 continue
      if(ipoly.eq.1) then
        db(1)=2.d0
        if (n.eq.1) return
        do 20 k=2,n
          dkm1=dble(k-1)
          db(k)=1.d0/(4.d0-1.d0/(dkm1*dkm1))
   20   continue
        return
      else if(ipoly.eq.2) then
        da(1)=.5d0
        db(1)=1.d0
        if(n.eq.1) return
        do 30 k=2,n
          da(k)=.5d0
          dkm1=dble(k-1)
          db(k)=.25d0/(4.d0-1.d0/(dkm1*dkm1))
   30   continue
        return
      else if(ipoly.eq.3) then
        db(1)=4.d0*datan(1.d0)
        if(n.eq.1) return
        db(2)=.5d0
        if(n.eq.2) return
        do 40 k=3,n
          db(k)=.25d0
   40   continue
        return
      else if(ipoly.eq.4) then
        db(1)=2.d0*datan(1.d0)
        if(n.eq.1) return
        do 50 k=2,n
          db(k)=.25d0
   50   continue
        return
      else if(ipoly.eq.5) then
        db(1)=4.d0*datan(1.d0)
        da(1)=.5d0
        if(n.eq.1) return
        do 60 k=2,n
          db(k)=.25d0
   60   continue
        return
      else if(ipoly.eq.6) then
        if(dal.le.-1.d0 .or. dbe.le.-1.d0) then
          iderr=1
          return
        else
          dalpbe=dal+dbe
          da(1)=(dbe-dal)/(dalpbe+2.d0)
          dt=(dalpbe+1.d0)*dlog(2.d0)+dlga(dal+1.d0)+dlga(dbe+1.d0)-
     *      dlga(dalpbe+2.d0)
          if(dt.gt.dlmach) then
            iderr=2
            db(1)=d1mach(2)
          else
            db(1)=dexp(dt)
          end if
          if(n.eq.1) return
          dal2=dal*dal
          dbe2=dbe*dbe
          da(2)=(dbe2-dal2)/((dalpbe+2.d0)*(dalpbe+4.d0))
          db(2)=4.d0*(dal+1.d0)*(dbe+1.d0)/((dalpbe+3.d0)*(dalpbe+
     *      2.d0)**2)
          if(n.eq.2) return
          do 70 k=3,n
            dkm1=dble(k-1)
            da(k)=.25d0*(dbe2-dal2)/(dkm1*dkm1*(1.d0+.5d0*dalpbe/dkm1)
     *        *(1.d0+.5d0*(dalpbe+2.d0)/dkm1))
            db(k)=.25d0*(1.d0+dal/dkm1)*(1.d0+dbe/dkm1)*(1.d0+dalpbe/
     *        dkm1)/((1.d0+.5d0*(dalpbe+1.d0)/dkm1)*(1.d0+.5d0*(dalpbe
     *      -1.d0)/dkm1)*(1.d0+.5d0*dalpbe/dkm1)**2)
   70     continue
          return
        end if
      else if(ipoly.eq.7) then
        if(dal.le.-1.d0) then
          iderr=1
          return
        else
          da(1)=dal+1.d0
          db(1)=dgamma(dal+1.d0,iderr)
          if(iderr.eq.2) db(1)=d1mach(2)
          if(n.eq.1) return
          do 80 k=2,n
            dkm1=dble(k-1)
            da(k)=2.d0*dkm1+dal+1.d0
            db(k)=dkm1*(dkm1+dal)
   80     continue
          return
        end if
      else if(ipoly.eq.8) then
        db(1)=dsqrt(4.d0*datan(1.d0))
        if(n.eq.1) return
        do 90 k=2,n
          db(k)=.5d0*dble(k-1)
   90   continue
        return
      else
        iderr=4
      end if
      return
      end
      double precision function dlga(dx)

c*********************************************************************72
c
c This routine evaluates the logarithm of the gamma function by a
c combination of recurrence and asymptotic approximation.
c
c  Reference:
c
c    William Cody, Kenneth Hillstrom,
c    Chebyshev Approximations for the Natural Logarithm of the
c    Gamma Function,
c    Mathematics of Computation,
c    Volume 21, Number 98, April 1967, pages 198-203.
c
c
c The entries in the next data statement are the numerators and
c denominators, respectively, of the quantities B[16]/(16*15),
c B[14]/(14*13),..., B[2]/(2*1), where B[2n] are the Bernoulli
c numbers.
c
      double precision dbnum,dbden,dx,d1mach,dc,dp,dy,dt,ds
      dimension dbnum(8),dbden(8)

      data dbnum/-3.617d3,1.d0,-6.91d2,1.d0,-1.d0,1.d0,-1.d0,1.d0/,
     *     dbden/1.224d5,1.56d2,3.6036d5,1.188d3,1.68d3,1.26d3,3.6d2,
     *1.2d1/
c
c The quantity  dprec  in the next statement is the number of decimal
c digits carried in double-precision floating-point arithmetic.
c
      dprec=-alog10(sngl(d1mach(3)))
      dc=.5d0*dlog(8.d0*datan(1.d0))
      dp=1.d0
      dy=dx
      y=sngl(dy)
c
c The quantity  y0  below is the threshold value beyond which asymptotic
c evaluation gives sufficient accuracy; see Eq. 6.1.42 in M. Abramowitz
c and I.A. Stegun,Handbook of Mathematical Functions''. The constants
c are .12118868... = ln(10)/19 and .05390522... = ln(|B[20]|/190)/19.
c
      y0=exp(.121189*dprec+.053905)
   10 if(y.gt.y0) goto 20
      dp=dy*dp
      dy=dy+1.d0
      y=sngl(dy)
      goto 10
   20 dt=1.d0/(dy*dy)
c
c The right-hand side of the next assignment statement is B[18]/(18*17).
c
      ds=4.3867d4/2.44188d5
      do 30 i=1,8
        ds=dt*ds+dbnum(i)/dbden(i)
   30 continue
      dlga=(dy-.5d0)*dlog(dy)-dy+dc+ds/dy-dlog(dp)
      return
      end
      double precision function dgamma(dx,iderr)

c*********************************************************************72
c
c This evaluates the gamma function for real positive  dx, using the
c function subroutine  dlga.
c
      double precision dx,dlmach,d1mach,dt,dlga
      dlmach=dlog(d1mach(2))
      iderr=0
      dt=dlga(dx)
      if(dt.ge.dlmach) then
        iderr=2
        dgamma=d1mach(2)
        return
      else
        dgamma=dexp(dt)
        return
      end if
      end
      subroutine dsti(n,ncap,dx,dw,dalpha,dbeta,ierr,dp0,dp1,dp2)

c*********************************************************************72
c
c This is a double-precision version of the routine  sti.
c
c  Reference:
c
c    Walter Gautschi,
c    On Generating Orthogonal Polynomials,
c    SIAM Journal on Scientific and Statistical Computing,
c    Volume 3, Number 3, 1982, pages 289-317.
c
      double precision dx,dw,dalpha,dbeta,dp0,dp1,dp2,dtiny,d1mach,
     *dhuge,dsum0,dsum1,dsum2,dt
      dimension dx(ncap),dw(ncap),dalpha(n),dbeta(n),dp0(ncap),
     *dp1(ncap),dp2(ncap)
      dtiny=10.d0*d1mach(1)
      dhuge=.1d0*d1mach(2)
      ierr=0
      if(n.le.0 .or. n.gt.ncap) then
        ierr=1
        return
      end if
      nm1=n-1
      dsum0=0.d0
      dsum1=0.d0
      do 10 m=1,ncap
        dsum0=dsum0+dw(m)
        dsum1=dsum1+dw(m)*dx(m)
   10 continue
      dalpha(1)=dsum1/dsum0
      dbeta(1)=dsum0
      if(n.eq.1) return
      do 20 m=1,ncap
        dp1(m)=0.d0
        dp2(m)=1.d0
   20 continue
      do 40 k=1,nm1
        dsum1=0.d0
        dsum2=0.d0
        do 30 m=1,ncap
          if(dw(m).eq.0.d0) goto 30
          dp0(m)=dp1(m)
          dp1(m)=dp2(m)
          dp2(m)=(dx(m)-dalpha(k))*dp1(m)-dbeta(k)*dp0(m)
          if(dabs(dp2(m)).gt.dhuge .or. dabs(dsum2).gt.dhuge) then
            ierr=k
            return
          end if
          dt=dw(m)*dp2(m)*dp2(m)
          dsum1=dsum1+dt
          dsum2=dsum2+dt*dx(m)
   30   continue
        if(dabs(dsum1).lt.dtiny) then
          ierr=-k
          return
        end if
        dalpha(k+1)=dsum2/dsum1
        dbeta(k+1)=dsum1/dsum0
        dsum0=dsum1
   40 continue
      return
      end
      subroutine gauss(n,alpha,beta,eps,zero,weight,ierr,e)

c*********************************************************************72
c
c Given  n  and a measure  dlambda, this routine generates the n-point
c Gaussian quadrature formula
c
c     integral over supp(dlambda) of f(x)dlambda(x)
c
c        = sum from k=1 to k=n of w(k)f(x(k)) + R(n;f).
c
c The nodes are returned as  zero(k)=x(k) and the weights as
c weight(k)=w(k), k=1,2,...,n. The user has to supply the recursion
c coefficients  alpha(k), beta(k), k=0,1,2,...,n-1, for the measure
c dlambda. The routine computes the nodes as eigenvalues, and the
c weights in term of the first component of the respective normalized
c eigenvectors of the n-th order Jacobi matrix associated with  dlambda.
c It uses a translation and adaptation of the algol procedure  imtql2,
c Numer. Math. 12, 1968, 377-383, by Martin and Wilkinson, as modified
c by Dubrulle, Numer. Math. 15, 1970, 450. See also Handbook for
c Autom. Comput., vol. 2 - Linear Algebra, pp.241-248, and the eispack
c routine  imtql2.
c
c        Input:  n - - the number of points in the Gaussian quadrature
c                      formula; type integer
c                alpha,beta - - arrays of dimension  n  to be filled
c                      with the values of  alpha(k-1), beta(k-1), k=1,2,
c                      ...,n
c                eps - the relative accuracy desired in the nodes
c                      and weights
c
c        Output: zero- array of dimension  n  containing the Gaussian
c                      nodes (in increasing order)  zero(k)=x(k), k=1,2,
c                      ...,n
c                weight - array of dimension  n  containing the
c                      Gaussian weights  weight(k)=w(k), k=1,2,...,n
c                ierr- an error flag equal to  0  on normal return,
c                      equal to  i  if the QR algorithm does not
c                      converge within 30 iterations on evaluating the
c                      i-th eigenvalue, equal to  -1  if  n  is not in
c                      range, and equal to  -2  if one of the beta's is
c                      negative.
c
c The array  e  is needed for working space.
c
      dimension alpha(n),beta(n),zero(n),weight(n),e(n)
      if(n.lt.1) then
        ierr=-1
        return
      end if
      ierr=0
      zero(1)=alpha(1)
      if(beta(1).lt.0.) then
        ierr=-2
        return
      end if
      weight(1)=beta(1)
      if (n.eq.1) return
      weight(1)=1.
      e(n)=0.
      do 100 k=2,n
        zero(k)=alpha(k)
        if(beta(k).lt.0.) then
          ierr=-2
          return
        end if
        e(k-1)=sqrt(beta(k))
        weight(k)=0.
  100 continue
      do 240 l=1,n
        j=0
c
c Look for a small subdiagonal element.
c
  105   do 110 m=l,n
          if(m.eq.n) goto 120
          if(abs(e(m)).le.eps*(abs(zero(m))+abs(zero(m+1)))) goto 120
  110   continue
  120   p=zero(l)
        if(m.eq.l) goto 240
        if(j.eq.30) goto 400
        j=j+1
c
c Form shift.
c
        g=(zero(l+1)-p)/(2.*e(l))
        r=sqrt(g*g+1.)
        g=zero(m)-p+e(l)/(g+sign(r,g))
        s=1.
        c=1.
        p=0.
        mml=m-l
c
c For i=m-1 step -1 until l do ...
c
        do 200 ii=1,mml
          i=m-ii
          f=s*e(i)
          b=c*e(i)
          if(abs(f).lt.abs(g)) goto 150
          c=g/f
          r=sqrt(c*c+1.)
          e(i+1)=f*r
          s=1./r
          c=c*s
          goto 160
  150     s=f/g
          r=sqrt(s*s+1.)
          e(i+1)=g*r
          c=1./r
          s=s*c
  160     g=zero(i+1)-p
          r=(zero(i)-g)*s +2.*c*b
          p=s*r
          zero(i+1)=g+p
          g=c*r-b
c
c Form first component of vector.
c
          f=weight(i+1)
          weight(i+1)=s*weight(i)+c*f
          weight(i)=c*weight(i)-s*f
  200   continue
        zero(l)=zero(l)-p
        e(l)=g
        e(m)=0.
        goto 105
  240 continue
c
c Order eigenvalues and eigenvectors.
c
      do 300 ii=2,n
        i=ii-1
        k=i
        p=zero(i)
        do 260 j=ii,n
          if(zero(j).ge.p) goto 260
          k=j
          p=zero(j)
  260   continue
        if(k.eq.i) goto 300
        zero(k)=zero(i)
        zero(i)=p
        p=weight(i)
        weight(i)=weight(k)
        weight(k)=p
  300 continue
      do 310 k=1,n
        weight(k)=beta(1)*weight(k)*weight(k)
  310 continue
      return
c
c Set error - no convergence to an eigenvalue after 30 iterations.
c
  400 ierr=l
      return
      end
      subroutine gchri(n,iopt,nu0,numax,eps,a,b,x,y,alpha,beta,nu,
     *  ierr,ierrc,fnu,rho,rold,s,s0,s1,s2)

c*********************************************************************72
c
c This routine implements the generalized Christoffel theorem, using
c the method of modified moments (cf. Section 4 of W. Gautschi,
c Minimal solutions of three-term recurrence relations and orthogonal
c polynomials'', Math. Comp. 36, 1981, 547-554). Given the recursion
c coefficients  a(k), b(k), k=0,1,...n, for the (monic) orthogonal
c polynomials with respect to some measure  dlambda(t), it generates
c the recursion coefficients  alpha(k), beta(k), k=0,1,2,...,n-1 for
c the measure
c
c         dlambda(t)/(t-x)        if iopt=1
c         dlambda(t)/{(t-x)**2+y**2} if iopt=2
c
c   Input:  n   - - the number of recurrence coefficients desired;
c                   type integer
c           iopt  - an integer selecting the desired weight distribution
c           nu0   - an integer estimating the starting backward
c                   recurrence index; in the absence of any better
c                   choice, take  nu0 = 3*n
c           numax - an integer controlling termination of backward
c                   recursion in case of nonconvergence; a conservative
c                   choice is  numax = 500
c           eps - - a relative error tolerance; type real
c           a,b - - arrays of dimension numax to be supplied with the
c                   recursion coefficients a(k)=alpha(k-1),b(k)=beta(k),
c                   k=1,2,...,numax, for the measure  dlambda
c           x,y - - real parameters defining the linear and quadratic
c                   divisors of  dlambda
c
c   Output: alpha,beta - arrays of dimension  n  containing the desired
c                   recursion coefficients  alpha(k-1), beta(k-1), k=1,
c                   2,...,n
c           nu  - - the backward recurrence index yielding convergence;
c                   in case of nonconvergence,  nu  will have the value
c                   numax
c           ierr  - an error flag, where
c                   ierr=0     on normal return
c                   ierr=1     if  iopt  is neither 1 nor 2
c                   ierr=nu0   if  nu0 > numax
c                   ierr=numax if the backward recurrence algorithm does
c                              not converge
c                   ierr=-1    if  n  is not in range
c           ierrc - an error flag inherited from the routine  cheb
c
c The arrays  fnu,s,s0,s1,s2  are working space. The routine calls
c upon the routines  knum  and  cheb.
c
      complex rho,rold,z
      dimension a(numax),b(numax),alpha(n),beta(n),fnu(*),rho(*),
     *rold(*),s(n),s0(*),s1(*),s2(*)
c
c The arrays  fnu,rho,rold,s0,s1,s2  are assumed to have dimension  2*n.
c
      if(n.lt.1) then
        ierr=-1
        return
      end if
      ierr=0
      nd=2*n
      ndm1=nd-1
c
c Linear divisor
c
      if(iopt.eq.1) then
c
c Generate the modified moments of  dlambda.
c
        z=cmplx(x,0.)
        call knum(ndm1,nu0,numax,z,eps,a,b,rho,nu,ierr,rold)
        do 10 k=1,nd
          fnu(k)=-real(rho(k))
   10   continue
c
c Compute the desired recursion coefficients by means of the modified
c Chebyshev algorithm.
c
        call cheb_r4 (n,a,b,fnu,alpha,beta,s,ierrc,s0,s1,s2)
        return
c
c Quadratic divisor
c
      else if(iopt.eq.2) then
c
c Generate the modified moments of  dlambda.
c
        y=abs(y)
        z=cmplx(x,y)
        call knum(ndm1,nu0,numax,z,eps,a,b,rho,nu,ierr,rold)
        do 20 k=1,nd
          fnu(k)=-aimag(rho(k))/y
   20   continue
c
c Compute the desired recursion coefficients by means of the modified
c Chebyshev algorithm.
c
        call cheb_r4 (n,a,b,fnu,alpha,beta,s,ierrc,s0,s1,s2)
        return
      else
        ierr=1
        return
      end if
      end
      subroutine kern(n,nu0,numax,z,eps,a,b,ker,nu,ierr,rold)

c*********************************************************************72
c
c This routine generates the kernels in the Gauss quadrature remainder
c term, namely
c
c           K(k)(z)=rho(k)(z)/pi(k)(z), k=0,1,2,...,n,
c
c where  rho(k)  are the output quantities of the routine  knum, and
c pi(k)  the (monic) orthogonal polynomials. The results are returned
c in the array  ker  as ker(k)=K(k-1)(z), k=1,2,...,n+1. All the other
c input and output parameters have the same meaning as in the routine
c knum.
c
      complex z,ker,rold,p0,p,pm1
      dimension a(numax),b(numax),ker(*),rold(*)
c
c The arrays  ker,rold  are assumed to have dimension  n+1.
c
      call knum(n,nu0,numax,z,eps,a,b,ker,nu,ierr,rold)
      p0=(0.,0.)
      p=(1.,0.)
      do 10 k=1,n
        pm1=p0
        p0=p
        p=(z-a(k))*p0-b(k)*pm1
        ker(k+1)=ker(k+1)/p
   10 continue
      return
      end
      subroutine knum(n,nu0,numax,z,eps,a,b,rho,nu,ierr,rold)

c*********************************************************************72
c
c This routine generates
c
c   rho(k)(z)=integral pi(k)(t)dlambda(t)/(z-t), k=0,1,2,...,n,
c
c where  pi(k)(t)  is the (monic) k-th degree orthogonal polynomial
c with respect to the measure  dlambda(t), and the integral is extended
c over the support of  dlambda. It is assumed that  z  is a complex
c number outside the smallest interval containing the support of
c dlambda. The quantities  rho(k)(z)  are computed as the first  n+1
c members of the minimal solution of the basic three-term recurrence
c relation
c
c      y(k+1)(z)=(z-a(k))y(k)(z)-b(k)y(k-1)(z), k=0,1,2,...,
c
c satisfied by the orthogonal polynomials  pi(k)(z).
c
c   Input:  n  - -  the largest integer  k  for which  rho(k)  is
c                   desired
c           nu0  -  an estimate of the starting backward recurrence
c                   index; if no better estimate is known, set
c                   nu0 = 3*n/2; for Jacobi, Laguerre and Hermite
c                   weight functions, estimates of  nu0  are generated
c                   respectively by the routines  nu0jac,nu0lag  and
c                   nu0her
c           numax - an integer larger than  n  cutting off backward
c                   recursion in case of nonconvergence; if  nu0
c                   exceeds  numax, then the routine aborts with the
c                   error flag  ierr  set equal to  nu0
c           z - - - the variable in  rho(k)(z); type complex
c           eps - - the relative accuracy to which the  rho(k)  are
c                   desired
c           a,b - - arrays of dimension  numax  to be supplied with the
c                   recurrence coefficients  a(k-1), b(k-1), k=1,2,...,
c                   numax.
c
c   Output: rho - - an array of dimension  n+1  containing the results
c                   rho(k)=rho(k-1)(z), k=1,2,...,n+1; type complex
c           nu  - - the starting backward recurrence index that yields
c                   convergence
c           ierr  - an error flag equal to zero on normal return, equal
c                   to  nu0  if  nu0 > numax, and equal to  numax in
c                   case of nonconvergence.
c
c The complex array  rold  of dimension  n+1  is used for working space.
c
      complex z,rho,rold,r
      dimension a(numax),b(numax),rho(*),rold(*)
c
c The arrays  rho,rold  are assumed to have dimension  n+1.
c
      ierr=0
      np1=n+1
      if(nu0.gt.numax) then
        ierr=nu0
        return
      end if
      if(nu0.lt.np1) nu0=np1
      nu=nu0-5
      do 10 k=1,np1
        rho(k)=(0.,0.)
   10 continue
   20 nu=nu+5
      if(nu.gt.numax) then
        ierr=numax
        goto 60
      end if
      do 30 k=1,np1
        rold(k)=rho(k)
   30 continue
      r=(0.,0.)
      do 40 j=1,nu
        j1=nu-j+1
        r=cmplx(b(j1),0.)/(z-cmplx(a(j1),0.)-r)
        if(j1.le.np1) rho(j1)=r
   40 continue
      do 50 k=1,np1
        if(cabs(rho(k)-rold(k)).gt.eps*cabs(rho(k))) goto 20
   50 continue
   60 if(n.eq.0) return
      do 70 k=2,np1
        rho(k)=rho(k)*rho(k-1)
   70 continue
      return
      end
      subroutine lancz(n,ncap,x,w,alpha,beta,ierr,p0,p1)

c*********************************************************************72
c
c This routine carries out the same task as the routine  sti, but
c uses the more stable Lanczos method. The meaning of the input
c and output parameters is the same as in the routine  sti. (This
c routine is adapted from the routine RKPW in W.B. Gragg and
c W.J. Harrod,The numerically stable reconstruction of Jacobi
c matrices from spectral data'', Numer. Math. 44, 1984, 317-335.)
c
      dimension x(ncap),w(ncap),alpha(n),beta(n),p0(ncap),p1(ncap)
      if(n.le.0 .or. n.gt.ncap) then
        ierr=1
        return
      else
        ierr=0
      end if
      do 10 i=1,ncap
        p0(i)=x(i)
        p1(i)=0.
   10 continue
      p1(1)=w(1)
      do 30 i=1,ncap-1
        pi=w(i+1)
        gam=1.
        sig=0.
        t=0.
        xlam=x(i+1)
        do 20 k=1,i+1
          rho=p1(k)+pi
          tmp=gam*rho
          tsig=sig
          if(rho.le.0.) then
            gam=1.
            sig=0.
          else
            gam=p1(k)/rho
            sig=pi/rho
          end if
          tk=sig*(p0(k)-xlam)-gam*t
          p0(k)=p0(k)-(tk-t)
          t=tk
          if(sig.le.0.) then
            pi=tsig*p1(k)
          else
            pi=(t**2)/sig
          end if
          tsig=sig
          p1(k)=tmp
   20   continue
   30 continue
      do 40 k=1,n
        alpha(k)=p0(k)
        beta(k)=p1(k)
   40 continue
      return
      end
      subroutine lob(n,alpha,beta,aleft,right,zero,weight,ierr,e,a,b)

c*********************************************************************72
c
c Given  n  and a measure  dlambda, this routine generates the
c (n+2)-point Gauss-Lobatto quadrature formula
c
c   integral over supp(dlambda) of f(x)dlambda(x)
c
c      = w(0)f(x(0)) + sum from k=1 to k=n of w(k)f(x(k))
c
c              + w(n+1)f(x(n+1)) + R(n;f).
c
c The nodes are returned as  zero(k)=x(k), the weights as  weight(k)
c =w(k), k=0,1,...,n,n+1. The user has to supply the recursion
c coefficients  alpha(k), beta(k), k=0,1,...,n,n+1, for the measure
c dlambda. The nodes and weights are computed in terms of the
c eigenvalues and first component of the normalized eigenvectors of
c a slightly modified Jacobi matrix of order  n+2. The routine calls
c upon the subroutine  gauss  and the function subroutine  r1mach.
c
c   Input:  n - -  the number of interior points in the Gauss-Lobatto
c                  formula; type integer
c           alpha,beta - arrays of dimension  n+2  to be supplied with
c                  the recursion coefficients  alpha(k-1), beta(k-1),
c                  k=1,2,...,n+2, of the underlying measure; the
c                  routine does not use  alpha(n+2), beta(n+2)
c           aleft,right - the prescribed left and right endpoints
c                  x(0)  and  x(n+1)  of the Gauss-Lobatto formula
c
c   Output: zero - an array of dimension  n+2  containing the nodes (in
c                  increasing order)  zero(k)=x(k), k=0,1,...,n,n+1
c           weight-an array of dimension  n+2  containing the weights
c                  weight(k)=w(k), k=0,1,...,n,n+1
c           ierr - an error flag inherited from the routine  gauss
c
c The arrays  e,a,b  are needed for working space.
c
      dimension alpha(*),beta(*),zero(*),weight(*),e(*),a(*),b(*)
c
c The arrays  alpha,beta,zero,weight,e,a,b  are assumed to have
c dimension  n+2.
c
      epsma=r1mach(3)
c
c epsma is the machine single precision.
c
      np1=n+1
      np2=n+2

      do k=1,np2
        a(k)=alpha(k)
        b(k)=beta(k)
      end do

      p0l=0.0
      p0r=0.0
      p1l=1.0
      p1r=1.0

      do k=1,np1
        pm1l=p0l
        p0l=p1l
        pm1r=p0r
        p0r=p1r
        p1l=(aleft-a(k))*p0l-b(k)*pm1l
        p1r=(right-a(k))*p0r-b(k)*pm1r
      end do

      det=p1l*p0r-p1r*p0l
      a(np2)=(aleft*p1l*p0r-right*p1r*p0l)/det
      b(np2)=(right-aleft)*p1l*p1r/det

      call gauss(np2,a,b,epsma,zero,weight,ierr,e)

      return
      end
      subroutine mccheb(n,ncapm,mc,mp,xp,yp,quad,eps,iq,idelta,
     *finl,finr,endl,endr,xfer,wfer,a,b,fnu,alpha,beta,ncap,kount,
     *ierr,be,x,w,xm,wm,s,s0,s1,s2)

c*********************************************************************72
c
c This is a multiple-component discretized modified Chebyshev
c algorithm, basically a modified Chebyshev algorithm in which the
c modified moments are discretized in the same manner as the inner
c product in the discretization procedure  mcdis. The input and
c output parameters are as in  mcdis. In addition, the arrays  a,b
c must be filled with the recursion coefficients  a(k-1),b(k-1),
c k=1,2,...,2*n-1, defining the modified moments. The arrays
c be,x,w,xm,wm,s,s0,s1,s2  are used for working space. The routine
c calls upon the subroutine  cheb. The routine exits immediately with
c ierr=-1  if  n  is not in range.
c
      dimension xp(*),yp(*),endl(mc),endr(mc),xfer(ncapm),
     *wfer(ncapm),a(*),b(*),fnu(*),alpha(n),beta(n),be(n),x(ncapm),
     *w(ncapm),xm(*),wm(*),s(n),s0(*),s1(*),s2(*)
      logical finl,finr
c
c The arrays  xp,yp  are assumed to have dimension  mp  if mp > 0,
c the arrays  a,b  dimension 2*n-1, the arrays  fnu,s0,s1,s2  dimension
c 2*n, and the arrays  xm,wm  dimension  mc*ncapm+mp.
c
      nd=2*n
      if(idelta.le.0) idelta=1
      if(n.lt.1) then
        ierr=-1
        return
      end if
c
c Initialization
c
      incap=1
      kount=-1
      ierr=0
      do 10 k=1,n
        beta(k)=0.
   10 continue
      ncap=(nd-1)/idelta
   20 do 30 k=1,n
        be(k)=beta(k)
   30 continue
      kount=kount+1
      if(kount.gt.1) incap=2**(kount/5)*n
      ncap=ncap+incap
      if(ncap.gt.ncapm) then
        ierr=ncapm
        return
      end if
c
c Discretization of the modified moments
c
      mtncap=mc*ncap
      do 50 i=1,mc
        im1tn=(i-1)*ncap
        if(iq.eq.1) then
          call quad(ncap,x,w,i,ierr)
        else
          call qgp(ncap,x,w,i,ierr,mc,finl,finr,endl,endr,xfer,wfer)
        end if
        if(ierr.ne.0) then
          ierr=i
          return
        end if
        do 40 k=1,ncap
          xm(im1tn+k)=x(k)
          wm(im1tn+k)=w(k)
   40   continue
   50 continue
      if(mp.ne.0) then
        do 60 k=1,mp
          xm(mtncap+k)=xp(k)
          wm(mtncap+k)=yp(k)
   60   continue
      end if
      mtnpmp=mtncap+mp
      do 90 k=1,nd
        km1=k-1
        sum=0.
        do 80 i=1,mtnpmp
          p1=0.
          p=1.
          if(k.gt.1) then
            do 70 l=1,km1
              pm1=p1
              p1=p
              p=(xm(i)-a(l))*p1-b(l)*pm1
   70       continue
          end if
          sum=sum+wm(i)*p
   80   continue
        fnu(k)=sum
   90 continue
c
c Computation of the desired recursion coefficients
c
      call cheb_r4 (n,a,b,fnu,alpha,beta,s,ierr,s0,s1,s2)
c
c In the following statement, the absolute value of the beta's is
c used to guard against failure in cases where the routine is applied
c to variable-sign weight functions and hence the positivity of the
c beta's is not guaranteed.
c
      do 100 k=1,n
        if(abs(beta(k)-be(k)).gt.eps*abs(beta(k))) goto 20
  100 continue
      return
      end
      subroutine mcdis(n,ncapm,mc,mp,xp,yp,quad,eps,iq,idelta,irout,
     *finl,finr,endl,endr,xfer,wfer,alpha,beta,ncap,kount,ierr,ie,
     *be,x,w,xm,wm,p0,p1,p2)

c*********************************************************************72
c
c This is a multiple-component discretization procedure as described in
c Section 4.3 of the companion paper. It generates to a relative
c accuracy of  eps  the recursion coefficients  alpha(k), beta(k),
c k=0,1,...,n-1, for the polynomials orthogonal with respect to a
c weight distribution consisting of the sum of  mc  continuous
c components and a discrete component with  mp  points. The continuous
c part of the spectrum is made up of  mc  weight functions, each
c supported on its own interval. These intervals may or may not be
c disjoint. The discretization of the inner product on the i-th
c interval is furnished either by a user-supplied subroutine  quad,
c or by the general-purpose subroutine  qgp  provided in this package,
c depending on whether  iq  is equal, or not equal, to  1, respectively.
c The user-supplied routine must have the form  quad(n,x,w,i,ierr)  and
c is assumed to supply the abscissas  x(k)  and weights  w(k), k=1,2,
c ...,n, to be used in approximating the i-th inner product
c
c               integral of p(x)*q(x)*wf(x,i)dx
c
c by the
c
c       sum over k from 1 to n of w(k)*p(x(k))*q(x(k)),
c
c                                        i=1,2,...,mc.
c
c The desired recurrence coefficients are then approximated by the
c recursion coefficients of the discrete orthogonal polynomials
c belonging to the discretized inner product, which in turn are
c computed by either the Stieltjes procedure or the Lanczos algorithm
c according as  irout  is equal to, or not equal to  1, respectively.
c Two error flags  ierr,ie  are provided which signal the occurrence
c of an error condition in the quadrature process, or in the routine
c sti  or  lancz  (whichever is used), respectively. The point spectrum
c is given through its abscissas  xp  and jumps  yp.
c
c If the quadrature routine  quad  has polynomial degree of exactness
c at least  id(n)  for each i, and if  id(n)/n = idelta + O(1/n)  as
c n  goes to infinity, then the procedure is designed to converge after
c one iteration, provided  idelta  is set with the appropriate
c integer. Normally,  idelta=1 (for interpolatory rules) or  idelta=2
c (for Gaussian rules). The default value is  idelta=1.
c
c    Input:  n    - - the number of recursion coefficients desired;
c                     type integer
c            ncapm  - a discretization parameter indicating an upper
c                     limit of the fineness of the discretization;
c                     ncapm=500  will usually be satisfactory; type
c                     integer
c            mc  - -  the number of disjoint intervals in the
c                     continuous part of the spectrum; type integer
c            mp  - -  the number of points in the discrete part of
c                     the spectrum; type integer. If there is no
c                     point spectrum, set  mp=0.
c            xp  - -  an array of dimension  mp  containing the
c                     abscissas of the point spectrum
c            yp  - -  an array of dimension  mp  containing the jumps
c                     of the point spectrum
c            quad  -  a subroutine determining the discretization of
c                     the inner product on each component interval,
c                     or a dummy routine if  iq  is not equal to  1
c                     (see below)
c            eps  - - the desired relative accuracy of the nonzero
c                     recursion coefficients; type real
c            iq   - - an integer selecting a user-supplied quadrature
c                     routine  quad  if  iq=1  or the ORTHPOL routine
c                     qgp  otherwise
c            idelta - a nonzero integer, typically  1  or  2, inducing
c                     fast convergence in the case of special quadrature
c                     routines
c            irout  - an integer selecting the routine for generating
c                     the recursion coefficients from the discrete
c                     inner product. Specifically,  irout=1  selects the
c                     routine  sti, whereas any other value selects the
c                     routine  lancz
c
c The logical variables  finl,finr, the arrays  endl,endr  of
c dimension  mc, and the arrays  xfer,wfer  of dimension  ncapm  are
c input variables to the subroutine  qgp  and are used (and hence need
c to be properly dimensioned) only if  iq  is not equal to  1.
c
c    Output:  alpha,beta - arrays of dimension n, holding as k-th
c                     element  alpha(k-1), beta(k-1), k=1,2,...,n,
c                     respectively
c             ncap  - an integer indicating the fineness of the
c                     discretization that yields convergence within
c                     the eps-tolerance
c             kount - the number of iterations used
c             ierr  - an error flag, equal to  0  on normal return,
c                     equal to  -1  if  n  is not in the proper range,
c                     equal to  i  if there is an error condition in
c                     the discretization of the i-th interval,
c                     and equal to  ncapm  if the discretized
c                     Stieltjes procedure does not converge within the
c                     discretization resolution specified by  ncapm
c             ie - -  an error flag inherited from the routine  sti
c                     or  lancz  (whichever is used)
c
c The array  be  of dimension  n, the arrays  x,w  of dimension  ncapm,
c and the arrays  xm,wm,p0,p1,p2  of dimension mc*ncapm+mp  are used
c for working space. The routine calls upon the subroutine  sti  or
c lancz, depending on the choice of  irout.
c
      dimension xp(*),yp(*),endl(mc),endr(mc),xfer(ncapm),wfer(ncapm),
     *alpha(n),beta(n),be(n),x(ncapm),w(ncapm),xm(*),wm(*),p0(*),p1(*),
     *p2(*)
      logical finl,finr
c
c The arrays  xp,yp  are assumed to have dimension  mp  if mp > 0, and
c the arrays  xm,wm,p0,p1,p2  dimension  mc*ncapm+mp.
c
      if(idelta.le.0) idelta=1
      if(n.lt.1) then
        ierr=-1
        return
      end if
c
c Initialization
c
      incap=1
      kount=-1
      ierr=0
      do 10 k=1,n
        beta(k)=0.
   10 continue
      ncap=(2*n-1)/idelta
   20 do 30 k=1,n
        be(k)=beta(k)
   30 continue
      kount=kount+1
      if(kount.gt.1) incap=2**(kount/5)*n
      ncap=ncap+incap
      if(ncap.gt.ncapm) then
        ierr=ncapm
        return
      end if
c
c Discretization of the inner product
c
      mtncap=mc*ncap
      do 50 i=1,mc
        im1tn=(i-1)*ncap
        if(iq.eq.1) then
          call quad(ncap,x,w,i,ierr)
        else
          call qgp(ncap,x,w,i,ierr,mc,finl,finr,endl,endr,xfer,
     *      wfer)
        end if
        if(ierr.ne.0) then
          ierr=i
          return
        end if
        do 40 k=1,ncap
          xm(im1tn+k)=x(k)
          wm(im1tn+k)=w(k)
   40   continue
   50 continue
      if(mp.ne.0) then
        do 60 k=1,mp
          xm(mtncap+k)=xp(k)
          wm(mtncap+k)=yp(k)
   60   continue
      end if
c
c Computation of the desired recursion coefficients
c
      if(irout.eq.1) then
        call sti(n,mtncap+mp,xm,wm,alpha,beta,ie,p0,p1,p2)
      else
        call lancz(n,mtncap+mp,xm,wm,alpha,beta,ie,p0,p1)
      end if
c
c In the following statement, the absolute value of the beta's is
c used to guard against failure in cases where the routine is applied
c to variable-sign weight functions and hence the positivity of the
c beta's is not guaranteed.
c
      do 70 k=1,n
        if(abs(beta(k)-be(k)).gt.eps*abs(beta(k))) goto 20
   70 continue
      return
      end
      function nu0her(n,z,eps)

c*********************************************************************72
c
c This is an auxiliary function routine providing a starting backward
c recurrence index for the Hermite measure that can be used in place
c of  nu0  in the routines  knum  and  dknum.
c
      complex z
      nu0her=2.*(sqrt(.5*real(n+1))+.25*alog(1./eps)/
     *  abs(aimag(z)))**2
      return
      end
      function nu0jac(n,z,eps)

c*********************************************************************72
c
c This is an auxiliary function routine providing a starting backward
c recurrence index for the Jacobi measure that can be used in place
c of  nu0  in the routines  knum  and  dknum.
c
      complex z
      pi=4.*atan(1.)
      x=real(z)
      y=abs(aimag(z))
      if(x.lt.1.) then
        if(x.lt.-1.) angle=.5*(2.*pi+atan(y/(x-1.))+atan(y/(x+1.)))
        if(x.eq.-1.) angle=.5*(1.5*pi-atan(.5*y))
        if(x.gt.-1.) angle=.5*(pi+atan(y/(x-1.))+atan(y/(x+1.)))
      else
        if(x.eq.1.) angle=.5*(.5*pi+atan(.5*y))
        if(x.gt.1.) angle=.5*(atan(y/(x-1.))+atan(y/(x+1.)))
      end if
      x2=x*x
      y2=y*y
      r=((x2-y2-1.)**2+4.*x2*y2)**.25
      r=sqrt((x+r*cos(angle))**2+(y+r*sin(angle))**2)
      nu0jac=real(n+1)+.5*alog(1./eps)/alog(r)
      return
      end
      function nu0lag(n,z,al,eps)

c*********************************************************************72
c
c This is an auxiliary function routine providing a starting backward
c recurrence index for the Laguerre measure that can be used in place
c of  nu0  in the routines  knum  and  dknum.
c
      complex z
      pi=4.*atan(1.)
      x=real(z)
      y=aimag(z)
      phi=.5*pi
      if(y.lt.0.) phi=1.5*pi
      if(x.eq.0.) goto 10
      phi=atan(y/x)
      if(y.gt.0. .and. x.gt.0.) goto 10
      phi=phi+pi
      if(x.lt.0.) goto 10
      phi=phi+pi
   10 nu0lag=(sqrt(real(n+1)+.5*(al+1.))+alog(1./eps)/(4.*(x*x+
     *  y*y)**.25*cos(.5*(phi-pi))))**2-.5*(al+1.)
      return
      end
      subroutine qgp(n,x,w,i,ierr,mc,finl,finr,endl,endr,xfer,wfer)

c*********************************************************************72
c
c This is a general-purpose discretization routine that can be used
c as an alternative to the routine  quad  in the multiple-component
c discretization procedure  mcdis. It takes no account of the special
c nature of the weight function involved and hence may result in slow
c convergence of the discretization procedure. This routine, therefore,
c should be used only as a last resort, when no better, more natural
c discretization can be found.
c
c It is assumed that there are  mc.ge.1  disjoint component intervals.
c The discretization is effected by the Fejer quadrature rule,
c suitably transformed to the respective interval. An interval that
c extends to minus infinity has to be indexed by  1; one that extends
c to plus infinity has to be indexed by  mc.
c
c The output variable  ierr  is given the value  0. Additional input
c parameters and working space used by this routine are as follows:
c
c          mc      - the number of component intervals; type integer
c          finl    - a logical variable to be set .true. if the
c                    extreme left interval is finite, and .false.
c                    otherwise
c          finr    - a logical variable to be set .true. if the
c                    extreme right interval is finite, and .false.
c                    otherwise
c          endl    - an array of dimension  mc  containing the left
c                    endpoints of the component intervals; if the
c                    first of these extends to minus infinity,  endl(1)
c                    can be set to an arbitrary value
c          endr    - an array of dimension  mc  containing the right
c                    endpoints of the component intervals; if the
c                    last of these extends to plus infinity,  endr(mc)
c                    can be set to an arbitrary value
c          xfer,wfer-working arrays holding the Fejer nodes and
c                    weights, respectively, for the interval [-1,1].
c
c The user has to supply the routine
c
c                     function wf(x,i),
c
c which evaluates the weight function at the point  x  on the i-th
c component interval. The routine also uses the subroutines  fejer,
c symtr  and  tr, which are appended.
c
      dimension x(n),w(n),endl(mc),endr(mc),xfer(*),wfer(*)
      logical finl,finr
c
c The arrays  xfer,wfer  are dimensioned in the routine  mcdis.
c
      ierr=0
      if(i.eq.1) call fejer(n,xfer,wfer)
      if(i.gt.1 .and. i.lt.mc) goto 60
      if(mc.eq.1) then
        if(finl.and.finr) goto 60
        if(finl) goto 20
        if(finr) goto 40
        do 10 k=1,n
          call symtr(xfer(k),phi,phi1)
          x(k)=phi
          w(k)=wfer(k)*wf(phi,i)*phi1
   10   continue
        return
      else
        if((i.eq.1.and.finl).or.(i.eq.mc.and.finr)) goto 60
        if(i.eq.1) goto 40
      end if
   20 do 30 k=1,n
        call tr(xfer(k),phi,phi1)
        x(k)=endl(mc)+phi
        w(k)=wfer(k)*wf(x(k),mc)*phi1
   30 continue
      return
   40 do 50 k=1,n
        call tr(-xfer(k),phi,phi1)
        x(k)=endr(1)-phi
        w(k)=wfer(k)*wf(x(k),1)*phi1
   50 continue
      return
   60 do 70 k=1,n
        x(k)=.5*((endr(i)-endl(i))*xfer(k)+endr(i)+endl(i))
        w(k)=.5*(endr(i)-endl(i))*wfer(k)*wf(x(k),i)
   70 continue
      return
      end
      subroutine symtr(t,phi,phi1)

c*********************************************************************72
c
c This implements a particular transformation  x=phi(t)  mapping
c the t-interval [-1,1] to the x-interval [-oo,oo].
c
c        input:   t
c        output:  phi=phi(t)
c                 phi1=derivative of phi(t)
c
      t2=t*t
      phi=t/(1.-t2)
      phi1=(t2+1.)/(t2-1.)**2
      return
      end
      subroutine tr(t,phi,phi1)

c*********************************************************************72
c
c This implements a particular transformation  x=phi(t)  mapping
c the t-interval [-1,1] to the x-interval [0,oo].
c
c         input:   t
c         output:  phi=phi(t)
c                  phi1=derivative of phi(t)
c
      phi=(1.+t)/(1.-t)
      phi1=2./(t-1.)**2
      return
      end
      subroutine fejer(n,x,w)

c*********************************************************************72
c
c This routine generates the n-point Fejer quadrature rule.
c
c         input:   n   - the number of quadrature nodes
c         output:  x,w - arrays of dimension  n  holding the quadrature
c                        nodes and weights, respectively; the nodes
c                        are ordered increasingly
c
      dimension x(n),w(n)
      pi=4.*atan(1.)
      nh=n/2
      np1h=(n+1)/2
      fn=real(n)
      do 10 k=1,nh
        x(n+1-k)=cos(.5*real(2*k-1)*pi/fn)
        x(k)=-x(n+1-k)
   10 continue
      if(2*nh.ne.n) x(np1h)=0.
      do 30 k=1,np1h
        c1=1.
        c0=2.*x(k)*x(k)-1.
        t=2.*c0
        sum=c0/3.
        do 20 m=2,nh
          c2=c1
          c1=c0
          c0=t*c1-c2
          sum=sum+c0/real(4*m*m-1)
   20   continue
        w(k)=2.*(1.-2.*sum)/fn
        w(n+1-k)=w(k)
   30 continue
      return
      end
      subroutine radau(n,alpha,beta,end,zero,weight,ierr,e,a,b)

c*********************************************************************72
c
c Given  n  and a measure  dlambda, this routine generates the
c (n+1)-point Gauss-Radau quadrature formula
c
c   integral over supp(dlambda) of f(t)dlambda(t)
c
c     = w(0)f(x(0)) + sum from k=1 to k=n of w(k)f(x(k)) + R(n;f).
c
c The nodes are returned as  zero(k)=x(k), the weights as  weight(k)
c =w(k), k=0,1,2,...,n. The user has to supply the recursion
c coefficients  alpha(k), beta(k), k=0,1,2,...,n, for the measure
c dlambda. The nodes and weights are computed as eigenvalues and
c in terms of the first component of the respective normalized
c eigenvectors of a slightly modified Jacobi matrix of order  n+1.
c To do this, the routine calls upon the subroutine  gauss. It also
c uses the function subroutine  r1mach.
c
c    Input:  n - -  the number of interior points in the Gauss-Radau
c                   formula; type integer
c            alpha,beta - arrays of dimension  n+1  to be supplied with
c                   the recursion coefficients  alpha(k-1), beta(k-1),
c                   k=1,2,...,n+1; the coefficient  alpha(n+1)  is not
c                   used by the routine
c            end -  the prescribed endpoint  x(0)  of the Gauss-Radau
c                   formula; type real
c
c    Output: zero - array of dimension  n+1  containing the nodes (in
c                   increasing order)  zero(k)=x(k), k=0,1,2,...,n
c            weight-array of dimension  n+1  containing the weights
c                   weight(k)=w(k), k=0,1,2,...,n
c            ierr - an error flag inherited from the routine  gauss
c
c The arrays  e,a,b  are needed for working space.
c
      dimension alpha(*),beta(*),zero(*),weight(*),e(*),a(*),b(*)
c
c The arrays  alpha,beta,zero,weight,e,a,b  are assumed to have
c dimension  n+1.
c
      epsma=r1mach(3)
c
c epsma is the machine single precision.
c
      np1=n+1
      do 10 k=1,np1
        a(k)=alpha(k)
        b(k)=beta(k)
   10 continue
      p0=0.
      p1=1.
      do 20 k=1,n
        pm1=p0
        p0=p1
        p1=(end-a(k))*p0-b(k)*pm1
   20 continue
      a(np1)=end-b(np1)*p0/p1
      call gauss(np1,a,b,epsma,zero,weight,ierr,e)
      return
      end
      subroutine recur(n,ipoly,al,be,a,b,ierr)

c*********************************************************************72
c
c This subroutine generates the coefficients  a(k),b(k), k=0,1,...,n-1,
c in the recurrence relation
c
c       p(k+1)(x)=(x-a(k))*p(k)(x)-b(k)*p(k-1)(x),
c                            k=0,1,...,n-1,
c
c       p(-1)(x)=0,  p(0)(x)=1,
c
c for some classical (monic) orthogonal polynomials, and sets  b(0)
c equal to the total mass of the weight distribution. The results are
c stored in the arrays  a,b,  which hold, respectively, the coefficients
c a(k-1),b(k-1), k=1,2,...,n.
c
c       Input:  n - - the number of recursion coefficients desired
c               ipoly-integer identifying the polynomial as follows:
c                     1=Legendre polynomial on (-1,1)
c                     2=Legendre polynomial on (0,1)
c                     3=Chebyshev polynomial of the first kind
c                     4=Chebyshev polynomial of the second kind
c                     5=Jacobi polynomial with parameters  al=-.5,be=.5
c                     6=Jacobi polynomial with parameters  al,be
c                     7=generalized Laguerre polynomial with
c                       parameter  al
c                     8=Hermite polynomial
c               al,be-input parameters for Jacobi and generalized
c                     Laguerre polynomials
c
c       Output: a,b - arrays containing, respectively, the recursion
c                     coefficients  a(k-1),b(k-1), k=1,2,...,n.
c               ierr -an error flag, equal to  0  on normal return,
c                     equal to  1  if  al  or  be  are out of range
c                     when  ipoly=6  or  ipoly=7, equal to  2  if  b(0)
c                     overflows when  ipoly=6  or  ipoly=7, equal to  3
c                     if  n  is out of range, and equal to  4  if  ipoly
c                     is not an admissible integer. In the case  ierr=2,
c                     the coefficient  b(0)  is set equal to the largest
c                     machine-representable number.
c
c The subroutine calls for the function subroutines  r1mach,gamma  and
c alga. The routines  gamma  and  alga, which are included in this file,
c evaluate respectively the gamma function and its logarithm for
c positive arguments. They are used only in the cases  ipoly=6  and
c ipoly=7.
c
      external gamma
      dimension a(n),b(n)
      if(n.lt.1) then
        ierr=3
        return
      end if
      almach=alog(r1mach(2))
      ierr=0
      do 10 k=1,n
        a(k)=0.
   10 continue
      if(ipoly.eq.1) then
        b(1)=2.
        if (n.eq.1) return
        do 20 k=2,n
          fkm1=real(k-1)
          b(k)=1./(4.-1./(fkm1*fkm1))
   20   continue
        return
      else if (ipoly.eq.2) then
        a(1)=.5
        b(1)=1.
        if(n.eq.1) return
        do 30 k=2,n
          a(k)=.5
          fkm1=real(k-1)
          b(k)=.25/(4.-1./(fkm1*fkm1))
   30   continue
        return
      else if(ipoly.eq.3) then
        b(1)=4.*atan(1.)
        if(n.eq.1) return
        b(2)=.5
        if(n.eq.2) return
        do 40 k=3,n
          b(k)=.25
   40   continue
        return
      else if(ipoly.eq.4) then
        b(1)=2.*atan(1.)
        if(n.eq.1) return
        do 50 k=2,n
          b(k)=.25
   50   continue
        return
      else if(ipoly.eq.5) then
        b(1)=4.*atan(1.)
        a(1)=.5
        if(n.eq.1) return
        do 60 k=2,n
          b(k)=.25
   60   continue
        return
      else if(ipoly.eq.6) then
        if(al.le.-1. .or. be.le.-1.) then
          ierr=1
          return
        else
          alpbe=al+be
          a(1)=(be-al)/(alpbe+2.)
          t=(alpbe+1.)*alog(2.)+alga(al+1.)+alga(be+1.)-
     *      alga(alpbe+2.)
          if(t.gt.almach) then
            ierr=2
            b(1)=r1mach(2)
          else
            b(1)=exp(t)
          end if
          if(n.eq.1) return
          al2=al*al
          be2=be*be
          a(2)=(be2-al2)/((alpbe+2.)*(alpbe+4.))
          b(2)=4.*(al+1.)*(be+1.)/((alpbe+3.)*(alpbe+2.)**2)
          if(n.eq.2) return
          do 70 k=3,n
            fkm1=real(k-1)
            a(k)=.25*(be2-al2)/(fkm1*fkm1*(1.+.5*alpbe/fkm1)*
     *        (1.+.5*(alpbe+2.)/fkm1))
            b(k)=.25*(1.+al/fkm1)*(1.+be/fkm1)*(1.+alpbe/fkm1)/
     *      ((1.+.5*(alpbe+1.)/fkm1)*(1.+.5*(alpbe-1.)/fkm1)
     *      *(1.+.5*alpbe/fkm1)**2)
   70     continue
          return
        end if
      else if(ipoly.eq.7) then
        if(al.le.-1.) then
          ierr=1
          return
        else
          a(1)=al+1.
          b(1)=gamma(al+1.,ierr)
          if(ierr.eq.2) b(1)=r1mach(2)
          if(n.eq.1) return
          do 80 k=2,n
            fkm1=real(k-1)
            a(k)=2.*fkm1+al+1.
            b(k)=fkm1*(fkm1+al)
   80     continue
          return
        end if
      else if(ipoly.eq.8) then
        b(1)=sqrt(4.*atan(1.))
        if(n.eq.1) return
        do 90 k=2,n
          b(k)=.5*real(k-1)
   90   continue
        return
      else
        ierr=4
      end if
      end
      function alga(x)

c*********************************************************************72
c
c This is an auxiliary function subroutine (not optimized in any
c sense) evaluating the logarithm of the gamma function for positive
c arguments  x. It is called by the subroutine  gamma. The integer  m0
c in the first executable statement is the smallest integer  m  such
c that  1*3*5* ... *(2*m+1)/(2**m)  is greater than or equal to the
c largest machine-representable number. The routine is based on a
c rational approximation valid on [.5,1.5] due to W.J. Cody and
c K.E. Hillstrom; see Math. Comp. 21, 1967, 198-203, in particular the
c case  n=7  in Table II. For the computation of  m0  it calls upon the
c function subroutines  t  and  r1mach. The former, appended below,
c evaluates the inverse function  t = t(y)  of  y = t ln t.
c
c  Reference:
c
c    William Cody, Kenneth Hillstrom,
c    Chebyshev Approximations for the Natural Logarithm of the
c    Gamma Function,
c    Mathematics of Computation,
c    Volume 21, Number 98, April 1967, pages 198-203.
c
      dimension cnum(8),cden(8)
      data cnum/4.120843185847770,85.68982062831317,243.175243524421,
     *-261.7218583856145,-922.2613728801522,-517.6383498023218,
     *-77.41064071332953,-2.208843997216182/,
     *cden/1.,45.64677187585908,377.8372484823942,951.323597679706,
     *846.0755362020782,262.3083470269460,24.43519662506312,
     *.4097792921092615/
c
c The constants in the statement below are  exp(1.)  and  .5*alog(8.).
c
      m0=2.71828*t((alog(r1mach(2))-1.03972)/2.71828)
      xi=aint(x)
      if(x-xi.gt..5) xi=xi+1.
      m=ifix(xi)-1
c
c Computation of log gamma on the standard interval (1/2,3/2]
c
      xe=x-real(m)
      snum=cnum(1)
      sden=cden(1)
      do 10 k=2,8
        snum=xe*snum+cnum(k)
        sden=xe*sden+cden(k)
   10 continue
      alga=(xe-1.)*snum/sden
c
c Computation of log gamma on (0,1/2]
c
      if(m.eq.-1) then
        alga=alga-alog(x)
        return
      else if(m.eq.0) then
        return
      else
c
c Computation of log gamma on (3/2,5/2]
c
        p=xe
        if(m.eq.1) then
          alga=alga+alog(p)
          return
        else
c
c Computation of log gamma for arguments larger than 5/2
c
          mm1=m-1
c
c The else-clause in the next statement is designed to avoid possible
c overflow in the computation of  p  in the if-clause, at the expense
c of computing many logarithms.
c
          if(m.lt.m0) then
            do 20 k=1,mm1
              p=(xe+real(k))*p
   20       continue
            alga=alga+alog(p)
            return
          else
            alga=alga+alog(xe)
            do 30 k=1,mm1
              alga=alga+alog(xe+real(k))
   30       continue
            return
          end if
        end if
      end if
      end
      function gamma(x,ierr)

c*********************************************************************72
c
c This evaluates the gamma function for real positive  x, using the
c function subroutines  alga  and  r1mach. In case of overflow, the
c routine returns the largest machine-representable number and the
c error flag  ierr=2.
c
      almach=alog(r1mach(2))
      ierr=0
      t=alga(x)
      if(t.ge.almach) then
        ierr=2
        gamma=r1mach(2)
        return
      else
        gamma=exp(t)
        return
      end if
      end
      function t(y)

c*********************************************************************72
c
c This evaluates the inverse function  t = t(y)  of y = t ln t  for
c nonnegative  y  to an accuracy of about one percent. For the
c approximation used, see pp. 51-52 in W. Gautschi,Computational
c aspects of three-term recurrence relations'', SIAM Rev. 9, 1967,
c 24-82.
c
      if(y.le.10.) then
        p=.000057941*y-.00176148
        p=y*p+.0208645
        p=y*p-.129013
        p=y*p+.85777
        t=y*p+1.0125
      else
        z=alog(y)-.775
        p=(.775-alog(z))/(1.+z)
        p=1./(1.+p)
        t=y*p/z
      end if
      return
      end
      subroutine sti(n,ncap,x,w,alpha,beta,ierr,p0,p1,p2)

c*********************************************************************72
c
c This routine applies Stieltjes's procedure'' to generate the recursion
c coefficients  alpha(k), beta(k) , k=0,1,...,n-1, for the discrete
c (monic) orthogonal polynomials associated with the inner product
c
c     (f,g)=sum over k from 1 to ncap of w(k)*f(x(k))*g(x(k)).
c
c The integer  n  must be between  1  and  ncap, inclusive; otherwise,
c there is an error exit with  ierr=1. The results are stored in the
c arrays  alpha, beta; the arrays  p0, p1, p2  are working arrays.
c
c If there is a threat of underflow or overflow in the calculation
c of the coefficients  alpha(k)  and  beta(k), the routine exits with
c the error flag  ierr  set equal to  -k  (in the case of underflow)
c or  +k  (in the case of overflow), where  k  is the recursion index
c for which the problem occurs. The former [latter] can often be avoided
c by multiplying all weights  w(k)  by a sufficiently large [small]
c scaling factor prior to entering the routine, and, upon exit, divide
c the coefficient  beta(0)  by the same factor.
c
c This routine should be used with caution if  n  is relatively close
c to  ncap, since there is a distinct possibility of numerical
c instability developing. (See W. Gautschi,Is the recurrence relation
c for orthogonal polynomials always stable?'', BIT, 1993, to appear.)
c In that case, the routine  lancz  should be used.
c
c  Reference:
c
c    Walter Gautschi,
c    On Generating Orthogonal Polynomials,
c    SIAM Journal on Scientific and Statistical Computing,
c    Volume 3, Number 3, 1982, pages 289-317.
c
      dimension x(ncap),w(ncap),alpha(n),beta(n),p0(ncap),p1(ncap),
     *p2(ncap)
      tiny=10.*r1mach(1)
      huge=.1*r1mach(2)
      ierr=0
      if(n.le.0 .or. n.gt.ncap) then
        ierr=1
        return
      end if
      nm1=n-1
c
c Compute the first alpha- and beta-coefficient.
c
      sum0=0.
      sum1=0.
      do 10 m=1,ncap
        sum0=sum0+w(m)
        sum1=sum1+w(m)*x(m)
   10 continue
      alpha(1)=sum1/sum0
      beta(1)=sum0
      if(n.eq.1) return
c
c Compute the remaining alpha- and beta-coefficients.
c
      do 20 m=1,ncap
        p1(m)=0.
        p2(m)=1.
   20 continue
      do 40 k=1,nm1
        sum1=0.
        sum2=0.
        do 30 m=1,ncap
c
c The following statement is designed to avoid an overflow condition
c in the computation of  p2(m)  when the weights  w(m)  go to zero
c faster (and underflow) than the  p2(m)  grow.
c
          if(w(m).eq.0.) goto 30
          p0(m)=p1(m)
          p1(m)=p2(m)
          p2(m)=(x(m)-alpha(k))*p1(m)-beta(k)*p0(m)
c
c Check for impending overflow.
c
          if(abs(p2(m)).gt.huge .or. abs(sum2).gt.huge) then
            ierr=k
            return
          end if
          t=w(m)*p2(m)*p2(m)
          sum1=sum1+t
          sum2=sum2+t*x(m)
   30   continue
c
c Check for impending underflow.
c
        if(abs(sum1).lt.tiny) then
          ierr=-k
          return
        end if
        alpha(k+1)=sum2/sum1
        beta(k+1)=sum1/sum0
        sum0=sum1
   40 continue
      return
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
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
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
