      subroutine bakslv(nr,n,a,x,b)

c*********************************************************************72
c
cc BAKSLV solves A'*X=B, A lower triangular.
c
c  Discussion:
c
c    solve  ax=b  where a is upper triangular matrix.
c    note that a is input as a lower triangular matrix and
c    that this routine takes its transpose implicitly.
c
c  Modified:
c
c    09 January 2009
c
c  Author:
c
c    Robert Schnabel, John Koontz, Barry Weiss
c
c  Reference:
c
c    John Dennis, Robert Schnabel,
c    Numerical Methods for Unconstrained Optimization 
c    and Nonlinear Equations,
c    SIAM, 1996,
c    ISBN13: 978-0-898713-64-0,
c    LC: QA402.5.D44.
c
c    Robert Schnabel, John Koontz, Barry Weiss,
c    A modular system of algorithms for unconstrained minimization,
c    Technical Report CU-CS-240-82, 
c    Computer Science Department,
c    University of Colorado at Boulder, 1982.
c
c  Parameters:
c
c nr           --> row dimension of matrix
c n            --> dimension of problem
c a(n,n)       --> lower triangular matrix (preserved)
c x(n)        <--  solution vector
c b(n)         --> right-hand side vector
c
c note
c ----
c if b is no longer required by calling routine,
c then vectors b and x may share the same storage.
c
      implicit double precision ( a-h, o-z )

      dimension a(nr,1)
      dimension b(n)
      dimension x(n)
c
c solve (l-transpose)x=b. (back solve)
c
      i=n
      x(i)=b(i)/a(i,i)
      if(n.eq. 1 ) return
   10 ip1=i
      i=i-1
      sum = 0.0D+00
      do j=ip1,n
        sum=sum+a(j,i)*x(j)
      end do
      x(i)=(b(i)-sum)/a(i,i)
      if(i.gt. 1 ) go to 10

      return
      end
      subroutine chlhsn(nr,n,a,epsm,sx,udiag)

c*********************************************************************72
c
cc CHLHSN Cholesky decomposes the perturbed model Hessian matrix.
c
c  Discussion:
c
c    find the l(l-transpose), written ll+, decomposition of the perturbed
c    model hessian matrix a+mu*i(where mu.ne.0 and i is identity matrix)
c    which is safely positive definite.  if a is safely positive definite
c    upon entry, then mu=0.
c
c  Modified:
c
c    09 January 2009
c
c  Author:
c
c    Robert Schnabel, John Koontz, Barry Weiss
c
c  Reference:
c
c    John Dennis, Robert Schnabel,
c    Numerical Methods for Unconstrained Optimization 
c    and Nonlinear Equations,
c    SIAM, 1996,
c    ISBN13: 978-0-898713-64-0,
c    LC: QA402.5.D44.
c
c    Robert Schnabel, John Koontz, Barry Weiss,
c    A modular system of algorithms for unconstrained minimization,
c    Technical Report CU-CS-240-82, 
c    Computer Science Department,
c    University of Colorado at Boulder, 1982.
c
c  Parameters:
c
c nr           --> row dimension of matrix
c n            --> dimension of problem
c a(n,n)      <--> on entry; "a" is model hessian (only lower
c                  triangular part and diagonal stored)
c                  on exit:  a contains l of ll+ decomposition of
c                  perturbed model hessian in lower triangular
c                  part and diagonal and contains hessian in upper
c                  triangular part and udiag
c epsm         --> machine epsilon
c sx(n)        --> diagonal scaling matrix for x
c udiag(n)    <--  on exit: contains diagonal of hessian
c
c internal variables
c ------------------
c tol              tolerance
c diagmn           minimum element on diagonal of a
c diagmx           maximum element on diagonal of a
c offmax           maximum off-diagonal element of a
c offrow           sum of off-diagonal elements in a row of a
c evmin            minimum eigenvalue of a
c evmax            maximum eigenvalue of a
c
c description
c -----------
c 1. if "a" has any negative diagonal elements, then choose mu>0
c such that the diagonal of a:=a+mu*i is all positive
c with the ratio of its smallest to largest element on the
c order of sqrt(epsm).
c
c 2. "a" undergoes a perturbed cholesky decomposition which
c results in an ll+ decomposition of a+d, where d is a
c non-negative diagonal matrix which is implicitly added to
c "a" during the decomposition if "a" is not positive definite.
c "a" is retained and not changed during this process by
c copying l into the upper triangular part of "a" and the
c diagonal into udiag.  then the cholesky decomposition routine
c is called.  on return, addmax contains maximum element of d.
c
c 3. if addmax=0, "a" was positive definite going into step 2
c and return is made to calling program.  otherwise,
c the minimum number sdd which must be added to the
c diagonal of a to make it safely strictly diagonally dominant
c is calculated.  since a+addmax*i and a+sdd*i are safely
c positive definite, choose mu=min(addmax,sdd) and decompose
c a+mu*i to obtain l.
c
      implicit double precision ( a-h, o-z )

      dimension a(nr,1),sx(n),udiag(n)
c
c scale hessian
c pre- and post- multiply "a" by inv(sx)
c
      do j=1,n
        do i=j,n
          a(i,j)=a(i,j)/(sx(i)*sx(j))
        end do
      end do
c
c step1
c -----
c note:  if a different tolerance is desired throughout this
c algorithm, change tolerance here:
c
      tol=sqrt(epsm)

      diagmx=a(1,1)
      diagmn=a(1,1)
      if(n.eq. 1 ) go to 40
      do i=2,n
        if(a(i,i).lt.diagmn) diagmn=a(i,i)
        if(a(i,i).gt.diagmx) diagmx=a(i,i)
      end do

   40 posmax=max(diagmx, 0.0D+00 )
c
c diagmn .le. 0
c
      if(diagmn.gt.posmax*tol) go to 110
c     if(diagmn.le.posmax*tol)
c     then
        amu=tol*(posmax-diagmn)-diagmn
        if(amu.ne. 0.0D+00 ) go to 90
c       if(amu.eq. 0.0D+00 )
c       then
c
c find largest off-diagonal element of a
c
          offmax=0.0D+00
          if(n.eq. 1 ) go to 70
          do 60 i=2,n
            im1=i-1
            do j=1,im1
              if(abs(a(i,j)).gt.offmax) offmax=abs(a(i,j))
            end do
   60     continue
   70     amu=offmax
          if(amu.ne. 0.0D+00 ) go to 80
c         if(amu.eq. 0.0D+00 )
c         then
            amu = 1.00D+00
            go to 90
c         else
   80       amu=amu* ( 1.00D+00 + tol )
c         end if
c       end if
c
c a=a + mu*i
c
   90   do 100 i=1,n
          a(i,i)=a(i,i)+amu
  100   continue
        diagmx=diagmx+amu
c     end if
c
c step2
c -----
c copy lower triangular part of "a" to upper triangular part
c and diagonal of "a" to udiag
c
  110 continue
      do 130 j=1,n
        udiag(j)=a(j,j)
        if(j.eq.n) go to 130
        jp1=j+1
        do 120 i=jp1,n
          a(j,i)=a(i,j)
  120   continue
  130 continue
c
      call choldc(nr,n,a,diagmx,tol,addmax)
c
c
c step3
c -----
c if addmax=0, "a" was positive definite going into step 2,
c the ll+ decomposition has been done, and we return.
c otherwise, addmax>0.  perturb "a" so that it is safely
c diagonally dominant and find ll+ decomposition
c
      if(addmax.le. 0.0D+00 ) go to 220
c     if(addmax.gt. 0.0D+00 )
c     then
c
c restore original "a" (lower triangular part and diagonal)
c
        do 150 j=1,n
          a(j,j)=udiag(j)
          if(j.eq.n) go to 150
          jp1=j+1
          do 140 i=jp1,n
            a(i,j)=a(j,i)
  140     continue
  150   continue
c
c find sdd such that a+sdd*i is safely positive definite
c note:  evmin<0 since a is not positive definite;
c
        evmin = 0.0D+00
        evmax=a(1,1)
        do 200 i=1,n
          offrow = 0.0D+00
          if(i.eq. 1 ) go to 170
          im1=i-1
          do 160 j=1,im1
            offrow=offrow+abs(a(i,j))
  160     continue
  170     if(i.eq.n) go to 190
          ip1=i+1
          do 180 j=ip1,n
            offrow=offrow+abs(a(j,i))
  180     continue
  190     evmin=min(evmin,a(i,i)-offrow)
          evmax=max(evmax,a(i,i)+offrow)
  200   continue
        sdd=tol*(evmax-evmin)-evmin
c
c perturb "a" and decompose again
c
        amu=min(sdd,addmax)
        do i=1,n
          a(i,i)=a(i,i)+amu
          udiag(i)=a(i,i)
        end do
c
c "a" now guaranteed safely positive definite
c
        call choldc(nr,n,a, 0.0D+00, tol, addmax )
c     end if
c
c unscale hessian and cholesky decomposition matrix
c
  220 do 260 j=1,n
        do i=j,n
          a(i,j)=sx(i)*a(i,j)
        end do
        if(j.eq. 1 ) go to 250
        jm1=j-1
        do 240 i=1,jm1
          a(i,j)=sx(i)*sx(j)*a(i,j)
  240   continue
  250   udiag(j)=udiag(j)*sx(j)*sx(j)
  260 continue
      return
      end
      subroutine choldc(nr,n,a,diagmx,tol,addmax)

c*********************************************************************72
c
cc CHOLDC finds the perturbed Cholesky decomposition of A+D.
c
c  Discussion:
c
c    find the perturbed l(l-transpose), written ll+, decomposition
c    of a+d, where d is a non-negative diagonal matrix added to a if
c    necessary to allow the cholesky decomposition to continue.
c
c  Modified:
c
c    09 January 2009
c
c  Author:
c
c    Robert Schnabel, John Koontz, Barry Weiss
c
c  Reference:
c
c    John Dennis, Robert Schnabel,
c    Numerical Methods for Unconstrained Optimization 
c    and Nonlinear Equations,
c    SIAM, 1996,
c    ISBN13: 978-0-898713-64-0,
c    LC: QA402.5.D44.
c
c    Robert Schnabel, John Koontz, Barry Weiss,
c    A modular system of algorithms for unconstrained minimization,
c    Technical Report CU-CS-240-82, 
c    Computer Science Department,
c    University of Colorado at Boulder, 1982.
c
c  Parameters:
c
c nr           --> row dimension of matrix
c n            --> dimension of problem
c a(n,n)      <--> on entry: matrix for which to find perturbed
c                       cholesky decomposition
c                  on exit:  contains l of ll+ decomposition
c                  in lower triangular part and diagonal of "a"
c diagmx       --> maximum diagonal element of "a"
c tol          --> tolerance
c addmax      <--  maximum amount implicitly added to diagonal of "a"
c                  in forming the cholesky decomposition of a+d
c internal variables
c ------------------
c aminl    smallest element allowed on diagonal of l
c amnlsq   =aminl**2
c offmax   maximum off-diagonal element in column of a
c
c
c description
c -----------
c the normal cholesky decomposition is performed.  however, if at any
c point the algorithm would attempt to set l(i,i)=sqrt(temp)
c with temp < tol*diagmx, then l(i,i) is set to sqrt(tol*diagmx)
c instead.  this is equivalent to adding tol*diagmx-temp to a(i,i)
c
      implicit double precision ( a-h, o-z )

      dimension a(nr,1)

      addmax = 0.0D+00
      aminl=sqrt(diagmx*tol)
      amnlsq=aminl*aminl
c
c form column j of l
c
      do 100 j=1,n
c
c find diagonal elements of l
c
        sum = 0.0D+00
        if(j.eq. 1 ) go to 20
        jm1=j-1
        do k=1,jm1
          sum=sum + a(j,k)*a(j,k)
        end do
   20   temp=a(j,j)-sum
        if(temp.lt.amnlsq) go to 30
c       if(temp.ge.aminl**2)
c       then
          a(j,j)=sqrt(temp)
          go to 60
c       else
c
c find maximum off-diagonal element in column
c
   30     offmax = 0.0D+00
          if(j.eq.n) go to 50
          jp1=j+1
          do i=jp1,n
            if(abs(a(i,j)).gt.offmax) offmax=abs(a(i,j))
          end do
   50     if(offmax.le.amnlsq) offmax=amnlsq
c
c add to diagonal element  to allow cholesky decomposition to continue
c
          a(j,j)=sqrt(offmax)
          addmax=max(addmax,offmax-temp)
c       end if
c
c find i,j element of lower triangular matrix
c
   60   if(j.eq.n) go to 100
        jp1=j+1
        do i=jp1,n
          sum = 0.0D+00
          if(j.eq. 1 ) go to 80
          jm1=j-1
          do k=1,jm1
            sum=sum+a(i,k)*a(j,k)
          end do
   80     a(i,j)=(a(i,j)-sum)/a(j,j)
        end do
  100 continue
      return
      end
      subroutine d1fcn(n,x,g)

c*********************************************************************72
c
cc D1FCN is a dummy routine for the analytic gradient.
c
c  Modified:
c
c    09 January 2009
c
c  Author:
c
c    Robert Schnabel, John Koontz, Barry Weiss
c
c  Reference:
c
c    John Dennis, Robert Schnabel,
c    Numerical Methods for Unconstrained Optimization 
c    and Nonlinear Equations,
c    SIAM, 1996,
c    ISBN13: 978-0-898713-64-0,
c    LC: QA402.5.D44.
c
c    Robert Schnabel, John Koontz, Barry Weiss,
c    A modular system of algorithms for unconstrained minimization,
c    Technical Report CU-CS-240-82, 
c    Computer Science Department,
c    University of Colorado at Boulder, 1982.
c
c  Parameters:
c
c dummy routine to prevent unsatisfied external diagnostic
c when specific analytic gradient function not supplied.
c
      implicit double precision ( a-h, o-z )

      dimension x(n),g(n)

      write ( *, 1000 )
      write ( *, 1010 )
      stop
 1000 format(' d1fcn - dummy gradient routine called.')
 1010 format(' d1fcn - program forced to stop.')
      end
      subroutine d2fcn(nr,n,x,h)

c*********************************************************************72
c
cc D2FCN is a dummy routine for the analytic Hessian.
c
c  Modified:
c
c    09 January 2009
c
c  Author:
c
c    Robert Schnabel, John Koontz, Barry Weiss
c
c  Reference:
c
c    John Dennis, Robert Schnabel,
c    Numerical Methods for Unconstrained Optimization 
c    and Nonlinear Equations,
c    SIAM, 1996,
c    ISBN13: 978-0-898713-64-0,
c    LC: QA402.5.D44.
c
c    Robert Schnabel, John Koontz, Barry Weiss,
c    A modular system of algorithms for unconstrained minimization,
c    Technical Report CU-CS-240-82, 
c    Computer Science Department,
c    University of Colorado at Boulder, 1982.
c
c  Parameters:
c
c dummy routine to prevent unsatisfied external diagnostic
c when specific analytic hessian function not supplied.
c
      implicit double precision ( a-h, o-z )

      dimension x(n),h(nr,1)

      write ( *, 1000 )
      write ( *, 1010 )
      stop
 1000 format(' d2fcn - dummy hessian routine called')
 1010 format(' d2fcn - program forced to stop')
      end
      subroutine dfault(n,typsiz,fscale,method,iexp,msg,ndigit,
     &     itnlim,iagflg,iahflg,iprint,lounit,dlt,gradtl,stepmx,steptl)

c*********************************************************************72
c
cc DFAULT sets default values for each input variable to the minimization package.
c
c  Discussion:
c
c    set default values for each input variable to
c    minimization algorithm.
c
c  Modified:
c
c    09 January 2009
c
c  Author:
c
c    Robert Schnabel, John Koontz, Barry Weiss
c
c  Reference:
c
c    John Dennis, Robert Schnabel,
c    Numerical Methods for Unconstrained Optimization 
c    and Nonlinear Equations,
c    SIAM, 1996,
c    ISBN13: 978-0-898713-64-0,
c    LC: QA402.5.D44.
c
c    Robert Schnabel, John Koontz, Barry Weiss,
c    A modular system of algorithms for unconstrained minimization,
c    Technical Report CU-CS-240-82, 
c    Computer Science Department,
c    University of Colorado at Boulder, 1982.
c
c  Parameters:
c
c n            --> dimension of problem
c typsiz(n)   <--  typical size for each component of x
c                  set to 1.
c fscale      <--  estimate of scale of minimization function
c                  set to 1.
c method      <--  algorithm to use to solve minimization problem
c                  set to 1
c iexp        <--  =0 if minimization function not expensive to evaluate
c                  set to 1
c msg         <--  message to inhibit certain automatic checks + output
c                  set to 0
c ndigit      <--  number of good digits in minimization function
c                  use negative number if function assumed to have
c                  maximum possible number of good digits.
c                  set to -1
c itnlim      <--  maximum number of allowable iterations
c                  set to 150.
c iagflg      <--  =0 if analytic gradient not supplied
c                  set to 0.
c iahflg      <--  =0 if analytic hessian not supplied
c                  set to 0.
c iprint       --> amount of output desired.
c                  set to 1
c lounit      <--  device to which to send output
c                  0=no output.
c                  set to 6.
c dlt         <--  trust region radius
c                  if unknown, set to -1.
c                  set to -1.
c gradtl      <--  tolerance at which gradient considered close enough
c                  to zero to terminate algorithm
c                  set to cube root of relative machine precision.
c stepmx      <--  value of zero to trip default maximum in optchk
c                  set to 0.0
c steptl      <--  tolerance at which successive iterates considered
c                  close enough to terminate algorithm
c                  set to square root of relative machine precision.
c
      implicit double precision ( a-h, o-z )

      double precision r8_epsilon
      dimension typsiz(n)
c
c set typical size of x and minimization function
c
      do i=1,n
        typsiz(i) = 1.0D+00
      end do
      fscale = 1.0D+00
c
c set tolerances
c
      dlt = -1.0D+00

      epsm = r8_epsilon ( )

      gradtl=epsm**( 1.0D+00 / 3.0D+00 )
      stepmx = 0.0D+00
      steptl=sqrt(epsm)
c
c set flags
c
      method=1
      iexp=1
      msg=0
      ndigit=-1
      itnlim=150
      iagflg=0
      iahflg=0
      iprint=1
      lounit=6

      return
      end
      subroutine dogdrv(nr,n,x,f,g,a,p,xpls,fpls,fcn,sx,stepmx,
     &     steptl,dlt,iretcd,mxtake,sc,wrk1,wrk2,wrk3,iprint,lounit)

c*********************************************************************72
c
cc DOGDRV finds the next Newton iterate by the double dogleg method.
c
c  Discussion:
c
c    find a next newton iterate (xpls) by the double dogleg method
c
c  Modified:
c
c    09 January 2009
c
c  Author:
c
c    Robert Schnabel, John Koontz, Barry Weiss
c
c  Reference:
c
c    John Dennis, Robert Schnabel,
c    Numerical Methods for Unconstrained Optimization 
c    and Nonlinear Equations,
c    SIAM, 1996,
c    ISBN13: 978-0-898713-64-0,
c    LC: QA402.5.D44.
c
c    Robert Schnabel, John Koontz, Barry Weiss,
c    A modular system of algorithms for unconstrained minimization,
c    Technical Report CU-CS-240-82, 
c    Computer Science Department,
c    University of Colorado at Boulder, 1982.
c
c  Parameters:
c
c nr           --> row dimension of matrix
c n            --> dimension of problem
c x(n)         --> old iterate x(k-1)
c f            --> function value at old iterate, f(x)
c g(n)         --> gradient  at old iterate, g(x), or approximate
c a(n,n)       --> cholesky decomposition of hessian
c                  in lower triangular part and diagonal
c p(n)         --> newton step
c xpls(n)     <--  new iterate x(k)
c fpls        <--  function value at new iterate, f(xpls)
c fcn          --> name of subroutine to evaluate function
c sx(n)        --> diagonal scaling matrix for x
c stepmx       --> maximum allowable step size
c steptl       --> relative step size at which successive iterates
c                  considered close enough to terminate algorithm
c dlt         <--> trust region radius
c                  (retain value between successive calls)
c iretcd      <--  return code
c                    =0 satisfactory xpls found
c                    =1 failed to find satisfactory xpls sufficiently
c                       distinct from x
c mxtake      <--  boolean flag indicating step of maximum length used
c sc(n)        --> workspace (current step)
c wrk1(n)      --> workspace (and place holding argument to tregup)
c wrk2(n)      --> workspace
c wrk3(n)      --> workspace
c iprint       --> amount of output desired.
c lounit       --> device to which to send output
c                  0=no output.
c
      implicit double precision ( a-h, o-z )

      dimension x(n),xpls(n),g(n),p(n)
      dimension sx(n)
      dimension sc(n),wrk1(n),wrk2(n),wrk3(n)
      dimension a(nr,1)
      logical fstdog,nwtake,mxtake
      external fcn

      iretcd=4
      fstdog=.true.
      tmp = 0.0D+00
      do i=1,n
        tmp=tmp+sx(i)*sx(i)*p(i)*p(i)
      end do
      rnwtln=sqrt(tmp)
      if(iprint.gt. 1 )write(lounit,1030) rnwtln

   20 continue
c
c find new step by double dogleg algorithm
c
      call dogstp(nr,n,g,a,p,sx,rnwtln,dlt,nwtake,fstdog,
     &     wrk1,wrk2,cln,eta,sc,iprint,lounit,stepmx)
c
c check new point and update trust region
c
      call tregup(nr,n,x,f,g,a,fcn,sc,sx,nwtake,stepmx,steptl,dlt,
     &     iretcd,wrk3,fplsp,xpls,fpls,mxtake,iprint,lounit,2,wrk1)

      if(iretcd.le. 1 ) return
      go to 20
 1000 format(42h dogdrv    initial trust region not given.,
     &       22h  compute cauchy step.)
 1010 format(18h dogdrv    alpha =,g16.8/
     &       18h dogdrv    beta  =,g16.8/
     &       18h dogdrv    dlt   =,g16.8/
     &       18h dogdrv    nwtake=,l1    )
 1020 format(28h dogdrv    current step (sc))
 1030 format(18h0dogdrv    rnwtln=,g16.8)
 1040 format(14h dogdrv       ,5(g16.8,3x))
      end
      subroutine dogstp(nr,n,g,a,p,sx,rnwtln,dlt,nwtake,fstdog,
     &     ssd,v,cln,eta,sc,iprint,lounit,stepmx)

c*********************************************************************72
c
cc DOGSTP finds the next double dogleg stepsize.
c
c  Discussion:
c
c    find new step by double dogleg algorithm
c
c  Modified:
c
c    09 January 2009
c
c  Author:
c
c    Robert Schnabel, John Koontz, Barry Weiss
c
c  Reference:
c
c    John Dennis, Robert Schnabel,
c    Numerical Methods for Unconstrained Optimization 
c    and Nonlinear Equations,
c    SIAM, 1996,
c    ISBN13: 978-0-898713-64-0,
c    LC: QA402.5.D44.
c
c    Robert Schnabel, John Koontz, Barry Weiss,
c    A modular system of algorithms for unconstrained minimization,
c    Technical Report CU-CS-240-82, 
c    Computer Science Department,
c    University of Colorado at Boulder, 1982.
c
c  Parameters:
c
c nr           --> row dimension of matrix
c n            --> dimension of problem
c g(n)         --> gradient at current iterate, g(x)
c a(n,n)       --> cholesky decomposition of hessian in
c                  lower part and diagonal
c p(n)         --> newton step
c sx(n)        --> diagonal scaling matrix for x
c rnwtln       --> newton step length
c dlt         <--> trust region radius
c nwtake      <--> boolean, =.true. if newton step taken
c fstdog      <--> boolean, =.true. if on first leg of dogleg
c ssd(n)      <--> workspace (cauchy step to the minimum of the
c                  quadratic model in the scaled steepest descent
c                  direction) (retain value between successive calls)
c v(n)        <--> workspace  (retain value between successive calls)
c cln         <--> cauchy length
c                  (retain value between successive calls)
c eta              (retain value between successive calls)
c sc(n)       <--  current step
c iprint       --> amount of output desired.
c lounit       --> device to which to send output
c                   0=no output.
c stepmx       --> maximum allowable step size
c
c internal variables
c ------------------
c cln              length of cauchy step
c
      implicit double precision ( a-h, o-z )

      double precision r8vec_dot
      dimension g(n),p(n)
      dimension sx(n)
      dimension sc(n),ssd(n),v(n)
      dimension a(nr,1)
      logical nwtake,fstdog
c
c can we take newton step
c
      if(rnwtln.gt.dlt) go to 20
c     if(rnwtln.le.dlt)
c     then
        nwtake=.true.
        do i=1,n
          sc(i)=p(i)
        end do
        dlt=rnwtln
        if(iprint.gt. 1 )write(lounit,1000)
        go to 140
c     else
c
c newton step too long
c cauchy step is on double dogleg curve
c
   20   nwtake=.false.
        if(.not.fstdog) go to 80
c       if(fstdog)
c       then
c
c         calculate double dogleg curve (ssd)
          fstdog=.false.

          alpha = 0.0D+00
          do i=1,n
            alpha=alpha + (g(i)*g(i))/(sx(i)*sx(i))
          end do

          beta = 0.0D+00
          do i=1,n
            tmp = 0.0D+00
            do j=i,n
              tmp=tmp + (a(j,i)*g(j))/(sx(j)*sx(j))
            end do
            beta=beta+tmp*tmp
          end do

          do i=1,n
            ssd(i)=-(alpha/beta)*g(i)/sx(i)
          end do

          cln=alpha*sqrt(alpha)/beta

          eta= 0.2D+00 + ( 0.8D+00 *alpha*alpha)
     &      /(-beta * r8vec_dot ( n, g, p ) )

          do i=1,n
            v(i)=eta*sx(i)*p(i) - ssd(i)
          end do

          if (dlt .eq. ( -1.00D+00 )) dlt = min(cln, stepmx)
          if(iprint.gt. 1)write(lounit,1030) alpha,beta,cln,eta
          if(iprint.gt. 1)write(lounit,1040)
          if(iprint.gt. 1)write(lounit,1090) (ssd(i),i=1,n)
          if(iprint.gt. 1)write(lounit,1050)
          if(iprint.gt. 1)write(lounit,1090) (v(i),i=1,n)
c       end if
   80   if(eta*rnwtln.gt.dlt) go to 100
c       if(eta*rnwtln .le. dlt)
c       then
c
c  take partial step in newton direction
c
          do i=1,n
            sc(i)=(dlt/rnwtln)*p(i)
          end do
          if(iprint.gt. 1 )write(lounit,1060)
          go to 140
c       else
  100     if(cln.lt.dlt) go to 120
c
c  if(cln.ge.dlt) then
c  take step in steepest descent direction
c
            do i=1,n
              sc(i)=(dlt/cln)*ssd(i)/sx(i)
            end do
            if(iprint.gt. 1 )write(lounit,1070)
            go to 140
c         else
c
c           calculate convex combination of ssd and eta*p
c           which has scaled length dlt
c
  120       dot1 = r8vec_dot ( n, v, ssd )
            dot2 = r8vec_dot ( n, v, v )
            alam=(-dot1+sqrt((dot1*dot1)-dot2*(cln*cln-dlt*dlt)))/dot2
            do i=1,n
              sc(i)=(ssd(i) + alam*v(i))/sx(i)
            end do
            if(iprint.gt. 1 )write(lounit,1080)
c         end if
c       end if
c     end if
  140 continue
      if(iprint.gt. 1)write(lounit,1010) fstdog,nwtake,rnwtln,dlt
      if(iprint.gt. 1)write(lounit,1020)
      if(iprint.gt. 1)write(lounit,1090) (sc(i),i=1,n)
      return
 1000 format(27h0dogstp    take newton step)
 1010 format(18h dogstp    fstdog=,l1/
     &       18h dogstp    nwtake=,l1/
     &       18h dogstp    rnwtln=,g16.8/
     &       18h dogstp    dlt   =,g16.8)
 1020 format(28h dogstp    current step (sc))
 1030 format(18h dogstp    alpha =,g16.8/
     &       18h dogstp    beta  =,g16.8/
     &       18h dogstp    cln   =,g16.8/
     &       18h dogstp    eta   =,g16.8)
 1040 format(28h dogstp    cauchy step (ssd))
 1050 format(12h dogstp    v)
 1060 format(48h0dogstp    take partial step in newton direction)
 1070 format(50h0dogstp    take step in steepest descent direction)
 1080 format(39h0dogstp    take convex combination step)
 1090 format(14h dogstp       ,5(g16.8,3x))
      end
      subroutine explain ( itrmcd )

c*********************************************************************72
c
cc EXPLAIN prints an explanation for the UNCMIN termination code.
c
c  Modified:
c
c    09 January 2009
c
c  Author:
c
c    Robert Schnabel, John Koontz, Barry Weiss
c
c  Reference:
c
c    John Dennis, Robert Schnabel,
c    Numerical Methods for Unconstrained Optimization 
c    and Nonlinear Equations,
c    SIAM, 1996,
c    ISBN13: 978-0-898713-64-0,
c    LC: QA402.5.D44.
c
c    Robert Schnabel, John Koontz, Barry Weiss,
c    A modular system of algorithms for unconstrained minimization,
c    Technical Report CU-CS-240-82, 
c    Computer Science Department,
c    University of Colorado at Boulder, 1982.
c
c  Parameters:
c
c    Input, integer ITRMCD, the UNCMIN termination code.
c
      implicit double precision ( a-h, o-z )

      integer itrmcd

      write(*,*)' '
      write(*,*)'EXPLAIN:'
      write(*,*)'  UNCMIN returned a termination code of ', itrmcd
      write(*,*)' '
      write(*,*)'This code has the following meaning:'
      write(*,*)' '

      if(itrmcd.eq.0)then

        write(*,*)'This termination code means UNCMIN received'
        write(*,*)'illegal or unacceptable input.'
      elseif(itrmcd.eq. 1)then
        write(*,*)'The gradient is relatively close to zero.'
        write(*,*)'The current iterate is probably a good solution.'
      elseif(itrmcd.eq. 2)then
        write(*,*)'Successive iterates are not changing much.'
        write(*,*)'The current iterate is probably a good solution.'
      elseif(itrmcd.eq. 3)then
        write(*,*)'The last global step could not find a point with '
        write(*,*)'lower function value than that achieved by XPLS.'
        write(*,*)' '
        write(*,*)'XPLS may be an approximate local minimum, or '
        write(*,*)' '
        write(*,*)'The function is too nonlinear for this program, or'
        write(*,*)'the stepsize tolerance is too large.'
      elseif(itrmcd.eq. 4)then
        write(*,*)'The iteration limit was exceeded.'
        write(*,*)'The current iterate does not satisfy the '
        write (*,*)'requirements.'
      elseif(itrmcd.eq. 5)then
        write(*,*)'maximum stepsize was exceeded 5 consecutive times.'
        write(*,*)' '
        write(*,*)'Either the function is unbounded below, or '
        write(*,*)'becomes asymptotic to a finite value from above in'
        write(*,*)'some direction, or'
        write(*,*)'the value of STEPMX is simply too small.'
      else
        write(*,*)'This is an unknown error code!'
      end if

      write(*,*)' '

      return
      end
      subroutine forslv(nr,n,a,x,b)

c*********************************************************************72
c
cc FORSLV solves A*X=B, where A is lower triangular.
c
c  Discussion:
c
c    solve  ax=b  where a is lower triangular matrix
c
c  Modified:
c
c    09 January 2009
c
c  Author:
c
c    Robert Schnabel, John Koontz, Barry Weiss
c
c  Reference:
c
c    John Dennis, Robert Schnabel,
c    Numerical Methods for Unconstrained Optimization 
c    and Nonlinear Equations,
c    SIAM, 1996,
c    ISBN13: 978-0-898713-64-0,
c    LC: QA402.5.D44.
c
c    Robert Schnabel, John Koontz, Barry Weiss,
c    A modular system of algorithms for unconstrained minimization,
c    Technical Report CU-CS-240-82, 
c    Computer Science Department,
c    University of Colorado at Boulder, 1982.
c
c  Parameters:
c
c nr           --> row dimension of matrix
c n            --> dimension of problem
c a(n,n)       --> lower triangular matrix (preserved)
c x(n)        <--  solution vector
c b(n)         --> right-hand side vector
c
c note
c ----
c if b is no longer required by calling routine,
c then vectors b and x may share the same storage.
c
      implicit double precision ( a-h, o-z )

      dimension a(nr,1),x(n),b(n)
c
c solve lx=b. (foreward solve)
c
      x(1)=b(1)/a(1,1)
      if(n.eq. 1) return
      do i=2,n
        sum = 0.0D+00
        im1=i-1
        do j=1,im1
          sum=sum+a(i,j)*x(j)
        end do
        x(i)=(b(i)-sum)/a(i,i)
      end do

      return
      end
      subroutine fstocd(n,x,fcn,sx,rnoise,g)

c*********************************************************************72
c
cc FSTOCD finds central difference approximation to the gradient.
c
c  Discussion:
c
c    find central difference approximation g to the first derivative
c    (gradient) of the function defined by fcn at the point x.
c
c  Modified:
c
c    09 January 2009
c
c  Author:
c
c    Robert Schnabel, John Koontz, Barry Weiss
c
c  Reference:
c
c    John Dennis, Robert Schnabel,
c    Numerical Methods for Unconstrained Optimization 
c    and Nonlinear Equations,
c    SIAM, 1996,
c    ISBN13: 978-0-898713-64-0,
c    LC: QA402.5.D44.
c
c    Robert Schnabel, John Koontz, Barry Weiss,
c    A modular system of algorithms for unconstrained minimization,
c    Technical Report CU-CS-240-82, 
c    Computer Science Department,
c    University of Colorado at Boulder, 1982.
c
c  Parameters:
c
c n            --> dimension of problem
c x            --> point at which gradient is to be approximated.
c fcn          --> name of subroutine to evaluate function.
c sx           --> diagonal scaling matrix for x.
c rnoise       --> relative noise in fcn (f(x)).
c g           <--  central difference approximation to gradient.
c
      implicit double precision ( a-h, o-z )

      dimension x(n)
      dimension sx(n)
      dimension g(n)
      external fcn
c
c find i th  stepsize, evaluate two neighbors in direction of i th
c unit vector, and evaluate i th  component of gradient.
c
      third = 1.0D+00 / 3.0D+00
      do i = 1, n
         stepi = rnoise**third * max(abs(x(i)), 1.0D+00 / sx(i))
         xtempi = x(i)
         x(i) = xtempi + stepi
         call fcn (n, x, fplus)
         x(i) = xtempi - stepi
         call fcn (n, x, fminus)
         x(i) = xtempi
         g(i) = (fplus - fminus) / ( 2.0D+00 * stepi)
      end do

      return
      end
      subroutine fstofd(nr,m,n,xpls,fcn,fpls,a,sx,rnoise,fhat,icase)

c*********************************************************************72
c
cc FSTOFD finds forward difference approximation to the gradient.
c
c  Discussion:
c
c    find first order forward finite difference approximation "a" to the
c    first derivative of the function defined by the subprogram "fname"
c    evaluated at the new iterate "xpls".
c
c
c    for optimization use this routine to estimate:
c    1) the first derivative (gradient) of the optimization function "fcn
c       analytic user routine has been supplied;
c    2) the second derivative (hessian) of the optimization function
c       if no analytic user routine has been supplied for the hessian but
c       one has been supplied for the gradient ("fcn") and if the
c       optimization function is inexpensive to evaluate
c
c    _m=1 (optimization) algorithm estimates the gradient of the function
c         (fcn).   fcn(x) # f: r(n)-->r(1)
c    _m=n (systems) algorithm estimates the jacobian of the function
c         fcn(x) # f: r(n)-->r(n).
c    _m=n (optimization) algorithm estimates the hessian of the optimizatio
c         function, where the hessian is the first derivative of "fcn"
c
c  Modified:
c
c    09 January 2009
c
c  Author:
c
c    Robert Schnabel, John Koontz, Barry Weiss
c
c  Reference:
c
c    John Dennis, Robert Schnabel,
c    Numerical Methods for Unconstrained Optimization 
c    and Nonlinear Equations,
c    SIAM, 1996,
c    ISBN13: 978-0-898713-64-0,
c    LC: QA402.5.D44.
c
c    Robert Schnabel, John Koontz, Barry Weiss,
c    A modular system of algorithms for unconstrained minimization,
c    Technical Report CU-CS-240-82, 
c    Computer Science Department,
c    University of Colorado at Boulder, 1982.
c
c  Parameters:
c
c nr           --> row dimension of matrix
c m            --> number of rows in a
c n            --> number of columns in a; dimension of problem
c xpls(n)      --> new iterate:  x(k)
c fcn          --> name of subroutine to evaluate function
c fpls(m)      --> _m=1 (optimization) function value at new iterate:
c                       fcn(xpls)
c                  _m=n (optimization) value of first derivative
c                       (gradient) given by user function fcn
c                  _m=n (systems)  function value of associated
c                       minimization function
c a(nr,n)     <--  finite difference approximation (see note).  only
c                  lower triangular matrix and diagonal are returned
c sx(n)        --> diagonal scaling matrix for x
c rnoise       --> relative noise in fcn (f(x))
c fhat(m)      --> workspace
c icase        --> =1 optimization (gradient)
c                  =2 systems
c                  =3 optimization (hessian)
c
c internal variables
c ------------------
c stepsz - stepsize in the j-th variable direction
c
      implicit double precision ( a-h, o-z )

      dimension xpls(n),fpls(m)
      dimension fhat(m)
      dimension sx(n)
      dimension a(nr,1)
      external fcn
c
c find j-th column of a
c each column is derivative of f(fcn) with respect to xpls(j)
c
      do j=1,n
        stepsz=sqrt(rnoise)*max(abs(xpls(j)), 1.0D+00 /sx(j))
        xtmpj=xpls(j)
        xpls(j)=xtmpj+stepsz
        call fcn(n,xpls,fhat)
        xpls(j)=xtmpj
        do i=1,m
          a(i,j)=(fhat(i)-fpls(i))/stepsz
        end do
      end do

      if(icase.ne. 3) return
c
c if computing hessian, a must be symmetric
c
      if(n.eq. 1) return
      nm1=n-1
      do j=1,nm1
        jp1=j+1
        do i=jp1,m
          a(i,j)=(a(i,j)+a(j,i)) / 2.0D+00
        end do
      end do

      return
      end
      subroutine grdchk(n,x,fcn,f,g,typsiz,sx,fscale,rnf,
     &     analtl,wrk1,msg,iprint,lounit)

c*********************************************************************72
c
cc GRDCHK compares the analytic gradient against a finite difference estimate.
c
c  Discussion:
c
c    check analytic gradient against estimated gradient
c
c  Modified:
c
c    09 January 2009
c
c  Author:
c
c    Robert Schnabel, John Koontz, Barry Weiss
c
c  Reference:
c
c    John Dennis, Robert Schnabel,
c    Numerical Methods for Unconstrained Optimization 
c    and Nonlinear Equations,
c    SIAM, 1996,
c    ISBN13: 978-0-898713-64-0,
c    LC: QA402.5.D44.
c
c    Robert Schnabel, John Koontz, Barry Weiss,
c    A modular system of algorithms for unconstrained minimization,
c    Technical Report CU-CS-240-82, 
c    Computer Science Department,
c    University of Colorado at Boulder, 1982.
c
c  Parameters:
c
c n            --> dimension of problem
c x(n)         --> estimate to a root of fcn
c fcn          --> name of subroutine to evaluate optimization function
c                  must be declared external in calling routine
c                       fcn:  r(n) --> r(1)
c f            --> function value:  fcn(x)
c g(n)         --> gradient:  g(x)
c typsiz(n)    --> typical size for each component of x
c sx(n)        --> diagonal scaling matrix:  sx(i)=1./typsiz(i)
c fscale       --> estimate of scale of objective function fcn
c rnf          --> relative noise in optimization function fcn
c analtl       --> tolerance for comparison of estimated and
c                  analytical gradients
c wrk1(n)      --> workspace
c msg         <--  message or error code
c                    on output: =-21, probable coding error of gradient
c iprint       --> amount of output desired.
c lounit       --> device to which to send output
c                  0=no output.
c
      implicit double precision ( a-h, o-z )

      dimension x(n),g(n)
      dimension sx(n),typsiz(n)
      dimension wrk1(n)
      external fcn
c
c compute first order finite difference gradient and compare to
c analytic gradient.
c
      call fstofd(1,1,n,x,fcn,f,wrk1,sx,rnf,wrk,1)
c
      ker=0
      do i=1,n
        gs=max(abs(f),fscale)/max(abs(x(i)),typsiz(i))
        if(abs(g(i)-wrk1(i)).gt.max(abs(g(i)),gs)*analtl) then
		  ker=1
		  write(*,*)'grdchk suspects an error in gradient ',i
		  end if
        end do

      if(ker.ne.0)then
        if(iprint.gt. 0)write(lounit,1000)
        if(iprint.gt. 1)write(lounit,1010)
        if(iprint.gt. 1)write(lounit,1020) (i,g(i),wrk1(i),i=1,n)
        msg=-21
		end if

      if(ker.eq.0.and.iprint.gt. 1)write(lounit,1030)
      return
 1000 format(47h grdchk    probable error in coding of analytic,
     &       19h gradient function.)
 1010 format(16h grdchk     comp,12x,8hanalytic,12x,8hestimate)
 1020 format(11h grdchk    ,i5,3x,g16.8,3x,g16.8)
 1030 format(' grdchk    user gradient seems correct')
      end
      subroutine heschk(nr,n,x,fcn,d1fcn,d2fcn,f,g,a,typsiz,sx,rnf,
     &     analtl,iagflg,udiag,wrk1,wrk2,msg,iprint,lounit)

c*********************************************************************72
c
cc HESCHK compares the analytic Hessian against a finite difference estimate.
c
c  Discussion:
c
c    check analytic hessian against estimated hessian
c    (this may be done only if the user supplied analytic hessian
c    d2fcn fills only the lower triangular part and diagonal of a)
c
c  Modified:
c
c    09 January 2009
c
c  Author:
c
c    Robert Schnabel, John Koontz, Barry Weiss
c
c  Reference:
c
c    John Dennis, Robert Schnabel,
c    Numerical Methods for Unconstrained Optimization 
c    and Nonlinear Equations,
c    SIAM, 1996,
c    ISBN13: 978-0-898713-64-0,
c    LC: QA402.5.D44.
c
c    Robert Schnabel, John Koontz, Barry Weiss,
c    A modular system of algorithms for unconstrained minimization,
c    Technical Report CU-CS-240-82, 
c    Computer Science Department,
c    University of Colorado at Boulder, 1982.
c
c  Parameters:
c
c nr           --> row dimension of matrix
c n            --> dimension of problem
c x(n)         --> estimate to a root of fcn
c fcn          --> name of subroutine to evaluate optimization function
c                  must be declared external in calling routine
c                       fcn:  r(n) --> r(1)
c d1fcn        --> name of subroutine to evaluate gradient of fcn.
c                  must be declared external in calling routine
c d2fcn        --> name of subroutine to evaluate hessian of fcn.
c                  must be declared external in calling routine
c f            --> function value:  fcn(x)
c g(n)        <--  gradient:  g(x)
c a(n,n)      <--  on exit:  hessian in lower triangular part and diag
c typsiz(n)    --> typical size for each component of x
c sx(n)        --> diagonal scaling matrix:  sx(i)=1./typsiz(i)
c rnf          --> relative noise in optimization function fcn
c analtl       --> tolerance for comparison of estimated and
c                  analytical gradients
c iagflg       --> =1 if analytic gradient supplied
c udiag(n)     --> workspace
c wrk1(n)      --> workspace
c wrk2(n)      --> workspace
c msg         <--> message or error code
c                    on input : if =1xx do not compare anal + est hess
c                    on output: =-22, probable coding error of hessian
c iprint       --> amount of output desired.
c lounit       --> device to which to send output
c                  0=no output.
c
      implicit double precision ( a-h, o-z )

      dimension x(n),g(n),a(nr,1)
      dimension typsiz(n),sx(n)
      dimension udiag(n),wrk1(n),wrk2(n)
      external fcn
      external d1fcn
      external d2fcn
c
c compute finite difference approximation a to the hessian.
c
      if(iagflg.eq. 1) call fstofd(nr,n,n,x,d1fcn,g,a,sx,rnf,wrk1,3)
      if(iagflg.ne. 1) call sndofd(nr,n,x,fcn,f,a,sx,rnf,wrk1,wrk2)
c
      ker=0
c
c copy lower triangular part of "a" to upper triangular part
c and diagonal of "a" to udiag
c
      do 20 j=1,n
        udiag(j)=a(j,j)
        if(j.eq.n) go to 20
        jp1=j+1
        do i=jp1,n
          a(j,i)=a(i,j)
        end do
   20 continue
c
c compute analytic hessian and compare to finite difference
c approximation.
c
      call d2fcn(nr,n,x,a)
      do 40 j=1,n
        hs=max(abs(g(j)), 1.0D+00 )/max(abs(x(j)),typsiz(j))
        if(abs(a(j,j)-udiag(j)).gt.max(abs(udiag(j)),hs)*analtl)
     &       ker=1
        if(j.eq.n) go to 40
        jp1=j+1
        do i=jp1,n
          if(abs(a(i,j)-a(j,i)).gt.max(abs(a(i,j)),hs)*analtl) ker=1
        end do
   40 continue

      if(ker.eq.0) go to 80
        if(iprint.gt. 0)write(lounit,1000)
        if(iprint.gt. 1)write(lounit,1010)
        do i=1,n
          if(i.eq. 1) go to 60
          im1=i-1
          do j=1,im1
            if(iprint.gt. 1)write(lounit,1020) i,j,a(i,j),a(j,i)
          end do
   60     if(iprint.gt. 1)write(lounit,1020) i,i,a(i,i),udiag(i)
        end do
        msg=-22

   80 continue
      if(ker.eq.0.and.iprint.gt. 1)write(lounit,1030)
      return
 1000 format(47h heschk    probable error in coding of analytic,
     &       18h hessian function.)
 1010 format(21h heschk      row  col,14x,8hanalytic,14x,10h(estimate))
 1020 format(11h heschk    ,2i5,2x,g16.8,2x,1h(,g16.8,1h))
 1030 format('0heschk    user hessian seems correct')
      end
      subroutine hookdr(nr,n,x,f,g,a,udiag,p,xpls,fpls,fcn,sx,stepmx,
     &     steptl,dlt,iretcd,mxtake,amu,dltp,phi,phip0,
     &     sc,xplsp,wrk0,epsm,itncnt,iprint,lounit)

c*********************************************************************72
c
cc HOOKDR finds the next Newton iterate by the More-Hebdon method.
c
c  Discussion:
c
c    find a next newton iterate (xpls) by the more-hebdon method
c
c  Modified:
c
c    09 January 2009
c
c  Author:
c
c    Robert Schnabel, John Koontz, Barry Weiss
c
c  Reference:
c
c    John Dennis, Robert Schnabel,
c    Numerical Methods for Unconstrained Optimization 
c    and Nonlinear Equations,
c    SIAM, 1996,
c    ISBN13: 978-0-898713-64-0,
c    LC: QA402.5.D44.
c
c    Robert Schnabel, John Koontz, Barry Weiss,
c    A modular system of algorithms for unconstrained minimization,
c    Technical Report CU-CS-240-82, 
c    Computer Science Department,
c    University of Colorado at Boulder, 1982.
c
c  Parameters:
c
c nr           --> row dimension of matrix
c n            --> dimension of problem
c x(n)         --> old iterate x(k-1)
c f            --> function value at old iterate, f(x)
c g(n)         --> gradient at old iterate, g(x), or approximate
c a(n,n)       --> cholesky decomposition of hessian in lower
c                  triangular part and diagonal.
c                  hessian in upper triangular part and udiag.
c udiag(n)     --> diagonal of hessian in a(.,.)
c p(n)         --> newton step
c xpls(n)     <--  new iterate x(k)
c fpls        <--  function value at new iterate, f(xpls)
c fcn          --> name of subroutine to evaluate function
c sx(n)        --> diagonal scaling matrix for x
c stepmx       --> maximum allowable step size
c steptl       --> relative step size at which successive iterates
c                  considered close enough to terminate algorithm
c dlt         <--> trust region radius
c iretcd      <--  return code
c                    =0 satisfactory xpls found
c                    =1 failed to find satisfactory xpls sufficiently
c                       distinct from x
c mxtake      <--  boolean flag indicating step of maximum length used
c amu         <--> (retain value between successive calls)
c dltp        <--> (retain value between successive calls)
c phi         <--> (retain value between successive calls)
c phip0       <--> (retain value between successive calls)
c sc(n)        --> workspace
c xplsp(n)     --> workspace
c wrk0(n)      --> workspace
c epsm         --> machine epsilon
c itncnt       --> iteration count
c iprint       --> amount of output desired.
c lounit       --> device to which to send output
c                  0=no output.
c
      implicit double precision ( a-h, o-z )

      dimension x(n),g(n),p(n),xpls(n),sx(n)
      dimension a(nr,1),udiag(n)
      dimension sc(n),xplsp(n),wrk0(n)
      logical mxtake,nwtake
      logical fstime
      external fcn

      iretcd=4
      fstime=.true.
      tmp = 0.0D+00
      do i=1,n
        tmp=tmp+sx(i)*sx(i)*p(i)*p(i)
      end do
      rnwtln=sqrt(tmp)
      if(iprint.gt. 1)write(lounit,1030) rnwtln
c
      if(itncnt.gt. 1) go to 50
c     if(itncnt.eq. 1)
c     then
        amu=0.0D+00
c
c       if first iteration and trust region not provided by user,
c       compute initial trust region.
c
        if(dlt.ne. ( -1.0D+00 )) go to 50
c       if(dlt.eq. ( -1.0D+00 ))
c       then
          alpha = 0.0D+00
          do i=1,n
            alpha=alpha+(g(i)*g(i))/(sx(i)*sx(i))
          end do

          beta=0.0D+00
          do i=1,n
            tmp = 0.0D+00
            do j=i,n
              tmp=tmp + (a(j,i)*g(j))/(sx(j)*sx(j))
            end do
            beta=beta+tmp*tmp
          end do
          dlt=alpha*sqrt(alpha)/beta
          dlt = min(dlt, stepmx)
          if(iprint.gt. 1)write(lounit,1000)
          if(iprint.gt. 1)write(lounit,1010) alpha,beta,dlt
c       end if
c     end if
c
   50 continue
c
c find new step by more-hebdon algorithm
c
      call hookst(nr,n,g,a,udiag,p,sx,rnwtln,dlt,amu,
     &     dltp,phi,phip0,fstime,sc,nwtake,wrk0,epsm,iprint,lounit)
      dltp=dlt
c
c check new point and update trust region
c
      call tregup(nr,n,x,f,g,a,fcn,sc,sx,nwtake,stepmx,steptl,
     &   dlt,iretcd,xplsp,fplsp,xpls,fpls,mxtake,iprint,lounit,3,udiag)
      if(iretcd.le. 1) return
      go to 50
 1000 format(43h hookdr    initial trust region not given. ,
     &       21h compute cauchy step.)
 1010 format(18h hookdr    alpha =,g16.8/
     &       18h hookdr    beta  =,g16.8/
     &       18h hookdr    dlt   =,g16.8)
 1020 format(28h hookdr    current step (sc))
 1030 format(18h0hookdr    rnwtln=,g16.8)
 1040 format(14h hookdr       ,5(g16.8,3x))
      end
      subroutine hookst(nr,n,g,a,udiag,p,sx,rnwtln,dlt,amu,
     &     dltp,phi,phip0,fstime,sc,nwtake,wrk0,epsm,iprint,lounit)

c*********************************************************************72
c
cc HOOKST finds the More-Hebdon stepsize.
c
c  Discussion:
c
c    find new step by more-hebdon algorithm
c
c  Modified:
c
c    09 January 2009
c
c  Author:
c
c    Robert Schnabel, John Koontz, Barry Weiss
c
c  Reference:
c
c    John Dennis, Robert Schnabel,
c    Numerical Methods for Unconstrained Optimization 
c    and Nonlinear Equations,
c    SIAM, 1996,
c    ISBN13: 978-0-898713-64-0,
c    LC: QA402.5.D44.
c
c    Robert Schnabel, John Koontz, Barry Weiss,
c    A modular system of algorithms for unconstrained minimization,
c    Technical Report CU-CS-240-82, 
c    Computer Science Department,
c    University of Colorado at Boulder, 1982.
c
c  Parameters:
c
c nr           --> row dimension of matrix
c n            --> dimension of problem
c g(n)         --> gradient at current iterate, g(x)
c a(n,n)       --> cholesky decomposition of hessian in
c                  lower triangular part and diagonal.
c                  hessian or approx in upper triangular part
c udiag(n)     --> diagonal of hessian in a(.,.)
c p(n)         --> newton step
c sx(n)        --> diagonal scaling matrix for n
c rnwtln       --> newton step length
c dlt         <--> trust region radius
c amu         <--> (retain value between successive calls)
c dltp         --> trust region radius at last exit from this routine
c phi         <--> (retain value between successive calls)
c phip0       <--> (retain value between successive calls)
c fstime      <--> boolean. =.true. if first entry to this routine
c                  during k-th iteration
c sc(n)       <--  current step
c nwtake      <--  boolean, =.true. if newton step taken
c wrk0         --> workspace
c epsm         --> machine epsilon
c iprint       --> amount of output desired.
c lounit       --> device to which to send output
c                  0=no output.
c
      implicit double precision ( a-h, o-z )

      double precision r8vec_norm_l2
      dimension g(n),p(n),sx(n),sc(n),wrk0(n)
      dimension a(nr,1),udiag(n)
      logical nwtake,done
      logical fstime
c
c hi and alo are constants used in this routine.
c change here if other values are to be substituted.
c
      hi = 1.5D+00
      alo = 0.75D+00
      if(rnwtln.gt.hi*dlt) go to 20
c
c  take newton step
c
        nwtake=.true.
        do i=1,n
          sc(i)=p(i)
        end do
        dlt=min(dlt,rnwtln)
        amu = 0.0D+00
        if(iprint.gt. 1)write(lounit,1000)
        return
c
c  newton step not taken
c
   20   continue
        if(iprint.gt. 1)write(lounit,1010)
        nwtake=.false.
        if(amu .le. 0.0D+00 ) go to 30
          amu=amu- (phi+dltp) *((dltp-dlt)+phi)/(dlt*phip)
          if(iprint.gt. 1)write(lounit,1050) amu
   30   continue
        phi=rnwtln-dlt
        if(.not.fstime) go to 50
          do i=1,n
            wrk0(i)=sx(i)*sx(i)*p(i)
          end do
c
c  solve l*y = (sx**2)*p
c
          call forslv(nr,n,a,wrk0,wrk0)

          phip0 = - ( r8vec_norm_l2 ( n, wrk0 ) )**2 / rnwtln

          fstime=.false.
   50   phip=phip0
        amulo=-phi/phip
        amuup = 0.0D+00
        do i=1,n
          amuup=amuup+(g(i)*g(i))/(sx(i)*sx(i))
        end do
        amuup=sqrt(amuup)/dlt
        done=.false.
        if(iprint.gt. 1)write(lounit,1050) amu
        if(iprint.gt. 1)write(lounit,1080) phi
        if(iprint.gt. 1)write(lounit,1090) phip
        if(iprint.gt. 1)write(lounit,1060) amulo
        if(iprint.gt. 1)write(lounit,1070) amuup
c
c  test value of amu; generate next amu if necessary
c
   70   continue
        if(done) return
        if(iprint.gt. 1)write(lounit,1110)
        if(amu.ge.amulo .and. amu.le.amuup) go to 80
          amu=max(sqrt(amulo*amuup),amuup * 1.0D-03 )
          if(iprint.gt. 1)write(lounit,1050) amu
   80   continue
c
c  copy (h,udiag) to l
c  where h <-- h+amu*(sx**2) (do not actually change (h,udiag))
c
        do 100 j=1,n
          a(j,j)=udiag(j) + amu*sx(j)*sx(j)
          if(j.eq.n) go to 100
          jp1=j+1
          do i=jp1,n
            a(i,j)=a(j,i)
          end do
  100   continue
c
c       factor h=l(l+)
c
        call choldc(nr,n,a, 0.0D+00, sqrt(epsm),addmax)
c
c       solve h*p = l(l+)*sc = -g
c
        do i=1,n
          wrk0(i)=-g(i)
        end do

        call lltslv(nr,n,a,sc,wrk0)
        if(iprint.gt. 1)write(lounit,1040)
        if(iprint.gt. 1)write(lounit,1120) (sc(i),i=1,n)
c
c  reset h.  note since udiag has not been destroyed we need do
c  nothing here.  h is in the upper part and in udiag, still intact
c
        stepln = 0.0D+00
        do i=1,n
          stepln=stepln + sx(i)*sx(i)*sc(i)*sc(i)
        end do

        stepln=sqrt(stepln)
        phi=stepln-dlt
        do i=1,n
          wrk0(i)=sx(i)*sx(i)*sc(i)
        end do
        call forslv(nr,n,a,wrk0,wrk0)

        phip = - ( r8vec_norm_l2 ( n, wrk0 ) )**2 / stepln

        if(iprint.gt. 1)write(lounit,1100) dlt,stepln
        if(iprint.gt. 1)write(lounit,1080) phi
        if(iprint.gt. 1)write(lounit,1090) phip
        if((alo*dlt.gt.stepln .or. stepln.gt.hi*dlt) .and.
     &       (amuup-amulo.gt. 0.0D+00 )) go to 140
c
c  sc is acceptable hookstep
c
          if(iprint.gt. 1)write(lounit,1030)
          done=.true.
          go to 70
c       else
c
c  sc not acceptable hookstep.  select new amu
c
  140     continue
          if(iprint.gt. 1)write(lounit,1020)
          amulo=max(amulo,amu-(phi/phip))
          if(phi.lt. 0.0D+00 ) amuup=min(amuup,amu)
          amu=amu-(stepln*phi)/(dlt*phip)
          if(iprint.gt. 1)write(lounit,1050) amu
          if(iprint.gt. 1)write(lounit,1060) amulo
          if(iprint.gt. 1)write(lounit,1070) amuup
          go to 70
c       end if
c     end if
c
 1000 format(27h0hookst    take newton step)
 1010 format(32h0hookst    newton step not taken)
 1020 format(31h hookst    sc is not acceptable)
 1030 format(27h hookst    sc is acceptable)
 1040 format(28h hookst    current step (sc))
 1050 format(18h hookst    amu   =,g16.8)
 1060 format(18h hookst    amulo =,g16.8)
 1070 format(18h hookst    amuup =,g16.8)
 1080 format(18h hookst    phi   =,g16.8)
 1090 format(18h hookst    phip  =,g16.8)
 1100 format(18h hookst    dlt   =,g16.8/
     &       18h hookst    stepln=,g16.8)
 1110 format(23h0hookst    find new amu)
 1120 format(14h hookst       ,5(g16.8,3x))
      end
      subroutine hsnint(nr,n,a,sx,method)

c*********************************************************************72
c
cc HSNINT approximates the initial Hessian by secant updates.
c
c  Discussion:
c
c    provide initial hessian when using secant updates
c
c  Modified:
c
c    09 January 2009
c
c  Author:
c
c    Robert Schnabel, John Koontz, Barry Weiss
c
c  Reference:
c
c    John Dennis, Robert Schnabel,
c    Numerical Methods for Unconstrained Optimization 
c    and Nonlinear Equations,
c    SIAM, 1996,
c    ISBN13: 978-0-898713-64-0,
c    LC: QA402.5.D44.
c
c    Robert Schnabel, John Koontz, Barry Weiss,
c    A modular system of algorithms for unconstrained minimization,
c    Technical Report CU-CS-240-82, 
c    Computer Science Department,
c    University of Colorado at Boulder, 1982.
c
c  Parameters:
c
c nr           --> row dimension of matrix
c n            --> dimension of problem
c a(n,n)      <--  initial hessian (lower triangular matrix)
c sx(n)        --> diagonal scaling matrix for x
c method       --> algorithm to use to solve minimization problem
c                    =1,2 factored secant method used
c                    =3   unfactored secant method used
c
      implicit double precision ( a-h, o-z )

      dimension a(nr,1),sx(n)

      do 20 j=1,n
        if(method.eq. 3) a(j,j)=sx(j)*sx(j)
        if(method.ne. 3) a(j,j)=sx(j)
        if(j.eq.n) go to 20
        jp1=j+1
        do i=jp1,n
          a(i,j) = 0.0D+00
        end do
   20 continue

      return
      end
      subroutine lltslv(nr,n,a,x,b)

c*********************************************************************72
c
cc LLTSLV solves A*X=B when A has been Cholesky factored.
c
c  Discussion:
c
c    solve ax=b where a has the form l(l-transpose)
c    but only the lower triangular part, l, is stored.
c
c  Modified:
c
c    09 January 2009
c
c  Author:
c
c    Robert Schnabel, John Koontz, Barry Weiss
c
c  Reference:
c
c    John Dennis, Robert Schnabel,
c    Numerical Methods for Unconstrained Optimization 
c    and Nonlinear Equations,
c    SIAM, 1996,
c    ISBN13: 978-0-898713-64-0,
c    LC: QA402.5.D44.
c
c    Robert Schnabel, John Koontz, Barry Weiss,
c    A modular system of algorithms for unconstrained minimization,
c    Technical Report CU-CS-240-82, 
c    Computer Science Department,
c    University of Colorado at Boulder, 1982.
c
c  Parameters:
c
c nr           --> row dimension of matrix
c n            --> dimension of problem
c a(n,n)       --> matrix of form l(l-transpose).
c                  on return a is unchanged.
c x(n)        <--  solution vector
c b(n)         --> right-hand side vector
c
c note
c ----
c if b is not required by calling program, then
c b and x may share the same storage.
c
      implicit double precision ( a-h, o-z )

      dimension a(nr,1),x(n),b(n)
c
c forward solve, result in x
c
      call forslv(nr,n,a,x,b)
c
c back solve, result in x
c
      call bakslv(nr,n,a,x,x)
      return
      end
      subroutine lnsrch(n,x,f,g,p,xpls,fpls,fcn,mxtake,
     &   iretcd,stepmx,steptl,sx,iprint,lounit)

c*********************************************************************72
c
cc LNSRCH finds the next Newton iterate by line search.
c
c  Discussion:
c
c    find a next newton iterate by line search.
c
c  Modified:
c
c    09 January 2009
c
c  Author:
c
c    Robert Schnabel, John Koontz, Barry Weiss
c
c  Reference:
c
c    John Dennis, Robert Schnabel,
c    Numerical Methods for Unconstrained Optimization 
c    and Nonlinear Equations,
c    SIAM, 1996,
c    ISBN13: 978-0-898713-64-0,
c    LC: QA402.5.D44.
c
c    Robert Schnabel, John Koontz, Barry Weiss,
c    A modular system of algorithms for unconstrained minimization,
c    Technical Report CU-CS-240-82, 
c    Computer Science Department,
c    University of Colorado at Boulder, 1982.
c
c  Parameters:
c
c n            --> dimension of problem
c x(n)         --> old iterate:   x(k-1)
c f            --> function value at old iterate, f(x)
c g(n)         --> gradient at old iterate, g(x), or approximate
c p(n)         --> non-zero newton step
c xpls(n)     <--  new iterate x(k)
c fpls        <--  function value at new iterate, f(xpls)
c fcn          --> name of subroutine to evaluate function
c iretcd      <--  return code
c mxtake      <--  boolean flag indicating step of maximum length used
c stepmx       --> maximum allowable step size
c steptl       --> relative step size at which successive iterates
c                  considered close enough to terminate algorithm
c sx(n)        --> diagonal scaling matrix for x
c iprint       --> amount of output desired.
c lounit       --> device to which to send output
c
c internal variables
c ------------------
c sln              newton length
c rln              relative length of newton step
c
      implicit double precision ( a-h, o-z )

      double precision r8vec_dot
      integer n,iretcd
      dimension sx(n)
      dimension x(n),g(n),p(n)
      dimension xpls(n)
      logical mxtake
      external fcn

      mxtake=.false.
      iretcd=2
      if(iprint.gt. 1)write(lounit,1040)
      if(iprint.gt. 1)write(lounit,1060) (p(i),i=1,n)

      tmp = 0.0D+00
      do i=1,n
        tmp=tmp+sx(i)*sx(i)*p(i)*p(i)
      end do

      sln=sqrt(tmp)
      if(sln.le.stepmx) go to 20
c
c newton step longer than maximum allowed
c
        scl=stepmx/sln
        call sclmul(n,scl,p,p)
        sln=stepmx
        if(iprint.gt. 1)write(lounit,1050)
        if(iprint.gt. 1)write(lounit,1060) (p(i),i=1,n)
   20 continue

      slp = r8vec_dot ( n, g, p )

      rln = 0.0D+00
      do i=1,n
        rln=max(rln,abs(p(i))/max(abs(x(i)), 1.0D+00 /sx(i)))
      end do

      rmnlmb=steptl/rln
      almbda = 1.0D+00
      if(iprint.gt. 1)write(lounit,1020) sln,slp,rmnlmb
c
c loop
c check if new iterate satisfactory.  generate new lambda if necessary.
c
   40 continue
      if(iretcd.lt. 2) return
      do i=1,n
        xpls(i)=x(i) + almbda*p(i)
      end do

      call fcn(n,xpls,fpls)
      if(iprint.gt. 1)write(lounit,1000) almbda
      if(iprint.gt. 1)write(lounit,1010)
      if(iprint.gt. 1)write(lounit,1060) (xpls(i),i=1,n)
      if(iprint.gt. 1)write(lounit,1030) fpls
      if(fpls.gt. f+slp * 1.0D-04 * almbda) go to 60
c
c solution found
c
        iretcd=0
        if(almbda.eq. 1.0D+00 .and. sln.gt. 0.99D+00 * stepmx) then
          mxtake=.true.
        end if
        go to 40
c
c solution not (yet) found
c
   60   if(almbda .ge. rmnlmb) go to 70
c
c no satisfactory xpls found sufficiently distinct from x
c
          iretcd=1
          go to 40
c
c calculate new lambda
c
   70     if(almbda.ne. 1.0D+00 ) go to 80
c
c first backtrack: quadratic fit
c
            tlmbda=-slp/( 2.0D+00 *(fpls-f-slp))
            go to 110
c
c all subsequent backtracks: cubic fit
c
   80       t1=fpls-f-almbda*slp
            t2=pfpls-f-plmbda*slp
            t3 = 1.0D+00 /(almbda-plmbda)
            a=t3*(t1/(almbda*almbda) - t2/(plmbda*plmbda))
            b=t3*(t2*almbda/(plmbda*plmbda)
     &           - t1*plmbda/(almbda*almbda) )
            disc=b*b- 3.0D+00 * a*slp
            if(disc.le. b*b) go to 90
c
c only one positive critical point, must be minimum
c
              tlmbda=(-b+sign ( 1.0D+00, a)*sqrt(disc)) 
     &          / ( 3.0D+00 * a )
              go to 100
c
c both critical points positive, first is minimum
c
   90         tlmbda=(-b-sign ( 1.0D+00, a )*sqrt(disc)) 
     &          / ( 3.0D+00 * a )
  100       if(tlmbda.gt. 0.5D+00 * almbda) tlmbda= 0.5D+00 *almbda
  110     plmbda=almbda
          pfpls=fpls

          if(tlmbda.ge. almbda * 0.1D+00 ) go to 120
            almbda=almbda * 0.1D+00
            go to 130
  120       almbda=tlmbda
  130 go to 40
 
 1000 format(18h lnsrch    almbda=,g16.8)
 1010 format(29h lnsrch    new iterate (xpls))
 1020 format(' lnsrch    newton length =',g16.8/
     &       ' lnsrch    slp =          ',g16.8/
     &       ' lnsrch    rmnlmb   =     ',g16.8)
 1030 format(19h lnsrch    f(xpls)=,g16.8)
 1040 format(26h0lnsrch    newton step (p))
 1050 format(' lnsrch      reduced newton step')
 1060 format(14h lnsrch       ,5(g16.8,3x))
      end
      subroutine mvmltl(nr,n,a,x,y)

c*********************************************************************72
c
cc MVMLTL computes Y = L * X, where L is lower triangular.
c
c  Discussion:
c
c    compute y=lx where l is a lower triangular matrix stored in a
c
c  Modified:
c
c    09 January 2009
c
c  Author:
c
c    Robert Schnabel, John Koontz, Barry Weiss
c
c  Reference:
c
c    John Dennis, Robert Schnabel,
c    Numerical Methods for Unconstrained Optimization 
c    and Nonlinear Equations,
c    SIAM, 1996,
c    ISBN13: 978-0-898713-64-0,
c    LC: QA402.5.D44.
c
c    Robert Schnabel, John Koontz, Barry Weiss,
c    A modular system of algorithms for unconstrained minimization,
c    Technical Report CU-CS-240-82, 
c    Computer Science Department,
c    University of Colorado at Boulder, 1982.
c
c  Parameters:
c
c nr           --> row dimension of matrix
c n            --> dimension of problem
c a(n,n)       --> lower triangular (n*n) matrix
c x(n)         --> operand vector
c y(n)        <--  result vector
c
c note
c ----
c x and y cannot share storage
c
      implicit double precision ( a-h, o-z )

      dimension a(nr,1),x(n),y(n)

      do i=1,n
        sum = 0.0D+00
        do j=1,i
          sum=sum+a(i,j)*x(j)
        end do
        y(i)=sum
      end do

      return
      end
      subroutine mvmlts(nr,n,a,x,y)

c*********************************************************************72
c
cc MVMLTS computes Y = A * X where A is in symmetric storage.
c
c  Discussion:
c
c    compute y=ax where "a" is a symmetric (n*n) matrix stored in its lower
c    triangular part and x,y are n-vectors
c
c  Modified:
c
c    09 January 2009
c
c  Author:
c
c    Robert Schnabel, John Koontz, Barry Weiss
c
c  Reference:
c
c    John Dennis, Robert Schnabel,
c    Numerical Methods for Unconstrained Optimization 
c    and Nonlinear Equations,
c    SIAM, 1996,
c    ISBN13: 978-0-898713-64-0,
c    LC: QA402.5.D44.
c
c    Robert Schnabel, John Koontz, Barry Weiss,
c    A modular system of algorithms for unconstrained minimization,
c    Technical Report CU-CS-240-82, 
c    Computer Science Department,
c    University of Colorado at Boulder, 1982.
c
c  Parameters:
c
c nr           --> row dimension of matrix
c n            --> dimension of problem
c a(n,n)       --> symmetric (n*n) matrix stored in
c                  lower triangular part and diagonal
c x(n)         --> operand vector
c y(n)        <--  result vector
c
c note
c ----
c x and y cannot share storage.
c
      implicit double precision ( a-h, o-z )

      dimension a(nr,1),x(n),y(n)

      do i=1,n
        sum = 0.0D+00
        do j=1,i
          sum=sum+a(i,j)*x(j)
        end do
        if(i.eq.n) go to 30
        ip1=i+1
        do j=ip1,n
          sum=sum+a(j,i)*x(j)
        end do
   30   y(i)=sum
      end do

      return
      end
      subroutine mvmltu(nr,n,a,x,y)

c*********************************************************************72
c
cc MVMLTU computes Y = L' * X, where L is lower triangular.
c
c  Discussion:
c
c    compute y=(l+)x where l is a lower triangular matrix stored in a
c    (l-transpose (l+) is taken implicitly)
c
c  Modified:
c
c    09 January 2009
c
c  Author:
c
c    Robert Schnabel, John Koontz, Barry Weiss
c
c  Reference:
c
c    John Dennis, Robert Schnabel,
c    Numerical Methods for Unconstrained Optimization 
c    and Nonlinear Equations,
c    SIAM, 1996,
c    ISBN13: 978-0-898713-64-0,
c    LC: QA402.5.D44.
c
c    Robert Schnabel, John Koontz, Barry Weiss,
c    A modular system of algorithms for unconstrained minimization,
c    Technical Report CU-CS-240-82, 
c    Computer Science Department,
c    University of Colorado at Boulder, 1982.
c
c  Parameters:
c
c nr           --> row dimension of matrix
c n            --> dimension of problem
c a(nr,1)       --> lower triangular (n*n) matrix
c x(n)         --> operand vector
c y(n)        <--  result vector
c
c note
c ----
c x and y cannot share storage
c
      implicit double precision ( a-h, o-z )

      dimension a(nr,1),x(n),y(n)

      do i=1,n
        sum = 0.0D+00
        do j=i,n
          sum=sum+a(j,i)*x(j)
        end do
        y(i)=sum
      end do

      return
      end
      subroutine optchk(n,x,typsiz,sx,fscale,gradtl,itnlim,ndigit,epsm,
     &     dlt,method,iexp,iagflg,iahflg,stepmx,msg,iprint,lounit)

c*********************************************************************72
c
cc OPTCHK checks the input for reasonableness.
c
c  Discussion:
c
c    check input for reasonableness
c
c  Modified:
c
c    09 January 2009
c
c  Author:
c
c    Robert Schnabel, John Koontz, Barry Weiss
c
c  Reference:
c
c    John Dennis, Robert Schnabel,
c    Numerical Methods for Unconstrained Optimization 
c    and Nonlinear Equations,
c    SIAM, 1996,
c    ISBN13: 978-0-898713-64-0,
c    LC: QA402.5.D44.
c
c    Robert Schnabel, John Koontz, Barry Weiss,
c    A modular system of algorithms for unconstrained minimization,
c    Technical Report CU-CS-240-82, 
c    Computer Science Department,
c    University of Colorado at Boulder, 1982.
c
c  Parameters:
c
c n            --> dimension of problem
c x(n)         --> on entry, estimate to root of fcn
c typsiz(n)   <--> typical size of each component of x
c sx(n)       <--  diagonal scaling matrix for x
c fscale      <--> estimate of scale of objective function fcn
c gradtl       --> tolerance at which gradient considered close
c                  enough to zero to terminate algorithm
c itnlim      <--> maximum number of allowable iterations
c ndigit      <--> number of good digits in optimization function fcn
c epsm         --> machine epsilon
c dlt         <--> trust region radius
c method      <--> algorithm indicator
c iexp        <--> expense flag
c iagflg      <--> =1 if analytic gradient supplied
c iahflg      <--> =1 if analytic hessian supplied
c stepmx      <--> maximum step size
c msg         <--> message and error code
c iprint       --> amount of output desired.
c lounit       --> device to which to send output
c
      implicit double precision ( a-h, o-z )

      dimension x(n),typsiz(n),sx(n)
c
c  check output unit first
c
      if(lounit.le.0)iprint=0
c
c  check that parameters only take on acceptable values.
c  if not, set them to default values.
c
      if(method.lt. 1 .or. method.gt. 3) method=1
      if(iagflg.ne. 1 ) iagflg=0
      if(iahflg.ne. 1 ) iahflg=0
      if(iexp.ne.0 ) iexp=1
c
c  check dimension of problem
c
      if(n.le.0) go to 40
c
c  compute scale matrix
c
      do i=1,n
        if(typsiz(i).eq. 0.0D+00 ) typsiz(i) = 1.0D+00
        if(typsiz(i).lt. 0.0D+00 ) typsiz(i)=-typsiz(i)
        sx(i)=1.0D+00 / typsiz(i)
      end do
c
c  check maximum step size
c
      if ( stepmx .le. 0.0D+00 ) then

        stpsiz = 0.0D+00
        do i = 1, n
          stpsiz = stpsiz + x(i)*x(i)*sx(i)*sx(i)
        end do
        stpsiz = sqrt(stpsiz)
        stepmx = max ( 1.0D+03 * stpsiz, 1.0D+03 )

      end if
c
c  check function scale
c
      if ( fscale .eq. 0.0D+00 ) then
        fscale = 1.0D+00
      else if ( fscale .lt. 0.0D+00 ) then
        fscale = - fscale
      end if
c
c  check gradient tolerance
c
      if(gradtl.lt. 0.0D+00 ) go to 50
c
c  check iteration limit
c
      if(itnlim.le.0) go to 60
c
c  check number of digits of accuracy in function fcn
c
      if(ndigit.eq.0) go to 70
      if(ndigit.lt.0) ndigit=-log10(epsm)
c
c  check trust region radius
c
      if(dlt.le. 0.0D+00 ) dlt = -1.0D+00
      if (dlt .gt. stepmx) dlt = stepmx
      return
c
c  error exits
c
   40 if(iprint.gt.0)write(lounit,1000) n
      msg=-1
      go to 80
   50 if(iprint.gt.0)write(lounit,1010) gradtl
      msg=-3
      go to 80
   60 if(iprint.gt.0)write(lounit,1020) itnlim
      msg=-4
      go to 80
   70 if(iprint.gt.0)write(lounit,1030) ndigit
      msg=-5
      go to 80
   80 return

 1000 format(32h0optchk    illegal dimension, n=,i5)
 1010 format(38h0optchk    illegal tolerance.  gradtl=,g16.8)
 1020 format(44h0optchk    illegal iteration limit.  itnlim=,i5)
 1030 format(52h0optchk    minimization function has no good digits.,
     &        9h  ndigit=,i5)
      end
      subroutine optdrv(nr,n,x,fcn,d1fcn,d2fcn,typsiz,fscale,
     &     method,iexp,msg,ndigit,itnlim,iagflg,iahflg,iprint,lounit,
     &     dlt,gradtl,stepmx,steptl,xpls,fpls,gpls,itrmcd,
     &     a,udiag,g,p,sx,wrk0,wrk1,wrk2,wrk3)

c*********************************************************************72
c
cc OPTDRV is a driver for the nonlinear optimization code.
c
c  Discussion:
c
c    driver for non-linear optimization problem
c
c  Modified:
c
c    09 January 2009
c
c  Author:
c
c    Robert Schnabel, John Koontz, Barry Weiss
c
c  Reference:
c
c    John Dennis, Robert Schnabel,
c    Numerical Methods for Unconstrained Optimization 
c    and Nonlinear Equations,
c    SIAM, 1996,
c    ISBN13: 978-0-898713-64-0,
c    LC: QA402.5.D44.
c
c    Robert Schnabel, John Koontz, Barry Weiss,
c    A modular system of algorithms for unconstrained minimization,
c    Technical Report CU-CS-240-82, 
c    Computer Science Department,
c    University of Colorado at Boulder, 1982.
c
c  Parameters:
c
c nr           --> row dimension of matrix
c n            --> dimension of problem
c x(n)         --> on entry: estimate to a root of fcn
c fcn          --> name of subroutine to evaluate optimization function
c                  must be declared external in calling routine
c                            fcn: r(n) --> r(1)
c d1fcn        --> (optional) name of subroutine to evaluate gradient
c                  of fcn.  must be declared external in calling routine
c d2fcn        --> (optional) name of subroutine to evaluate hessian of
c                  of fcn.  must be declared external in calling routine
c typsiz(n)    --> typical size for each component of x
c fscale       --> estimate of scale of objective function
c method       --> algorithm to use to solve minimization problem
c                    =1 line search
c                    =2 double dogleg
c                    =3 more-hebdon
c iexp         --> =1 if optimization function fcn is expensive to
c                  evaluate, =0 otherwise.  if set then hessian will
c                  be evaluated by secant update instead of
c                  analytically or by finite differences
c msg         <--> on input:  (.gt.0) message to inhibit certain
c                    automatic checks
c                  on output: (.lt.0) error code; =0 no error
c ndigit       --> number of good digits in optimization function fcn
c itnlim       --> maximum number of allowable iterations
c iagflg       --> =1 if analytic gradient supplied
c iahflg       --> =1 if analytic hessian supplied
c iprint       --> amount of output desired.
c lounit       --> device to which to send output
c dlt          --> trust region radius
c gradtl       --> tolerance at which gradient considered close
c                  enough to zero to terminate algorithm
c stepmx       --> maximum allowable step size
c steptl       --> relative step size at which successive iterates
c                  considered close enough to terminate algorithm
c xpls(n)     <--> on exit:  xpls is local minimum
c fpls        <--> on exit:  function value at solution, xpls
c gpls(n)     <--> on exit:  gradient at solution xpls
c itrmcd      <--  termination code
c a(n,n)       --> workspace for hessian (or estimate)
c                  and its cholesky decomposition
c udiag(n)     --> workspace (for diagonal of hessian)
c g(n)         --> workspace (for gradient at current iterate)
c p(n)         --> workspace for step
c sx(n)        --> workspace (for diagonal scaling matrix)
c wrk0(n)      --> workspace
c wrk1(n)      --> workspace
c wrk2(n)      --> workspace
c wrk3(n)      --> workspace
c
c
c internal variables
c ------------------
c analtl           tolerance for comparison of estimated and
c                  analytical gradients and hessians
c epsm             machine epsilon
c f                function value: fcn(x)
c itncnt           current iteration, k
c rnf              relative noise in optimization function fcn.
c                       noise=10.**(-ndigit)
c
      implicit double precision ( a-h, o-z )

      double precision r8_epsilon
      dimension x(n),xpls(n),g(n),gpls(n),p(n)
      dimension typsiz(n),sx(n)
      dimension a(nr,1),udiag(n)
      dimension wrk0(n),wrk1(n),wrk2(n),wrk3(n)
      logical mxtake,noupdt
      external fcn,d1fcn,d2fcn
c
c initialization
c
      do i=1,n
        p(i) = 0.0D+00
      end do

      itncnt=0
      iretcd=-1

      epsm = r8_epsilon ( )

      call optchk(n,x,typsiz,sx,fscale,gradtl,itnlim,ndigit,epsm,
     &     dlt,method,iexp,iagflg,iahflg,stepmx,msg,iprint,lounit)
      if(msg.lt.0) return
      rnf=max( 10.0D+00 **(-ndigit),epsm)
      analtl=max( 1.0D-02, sqrt(rnf))

      if(mod(msg/8,2).eq. 1) go to 20
      if(iprint.gt. 1)write(lounit,1010)
      if(iprint.gt. 1)write(lounit,1000) (typsiz(i),i=1,n)
      if(iprint.gt. 1)write(lounit,1020)
      if(iprint.gt. 1)write(lounit,1000) (sx(i),i=1,n)
      if(iprint.gt. 1)write(lounit,1030) fscale
      if(iprint.gt. 1)write(lounit,1040)
     & ndigit,iagflg,iahflg,iexp,method,itnlim,epsm
      if(iprint.gt. 1)write(lounit,1050) stepmx,steptl,gradtl,dlt,rnf,
     & analtl
   20 continue
c
c evaluate fcn(x)
c
      call fcn(n,x,f)
c
c evaluate analytic or finite difference gradient and check analytic
c gradient, if requested.
c
      if (iagflg .eq. 1) go to 30
c     if (iagflg .eq. 0)
c     then
        call fstofd (1, 1, n, x, fcn, f, g, sx, rnf, wrk, 1)
        go to 40
c
   30 call d1fcn (n, x, g)
      if (mod(msg/2,2) .eq. 1) go to 40
c     if (mod(msg/2,2).eq.0)
c     then
        call grdchk(n,x,fcn,f,g,typsiz,sx,fscale,
     &  rnf,analtl,wrk1,msg,iprint,lounit)
        if (msg .lt. 0) return
   40 continue

      call optstp(n,x,f,g,wrk1,itncnt,icscmx,
     &            itrmcd,gradtl,steptl,sx,fscale,itnlim,iretcd,mxtake,
     &            iprint,lounit,msg)
      if(itrmcd.ne.0) go to 250

      if(iexp.ne. 1) go to 50
c
c if optimization function expensive to evaluate (iexp=1), then
c hessian will be obtained by secant updates.  get initial hessian.
c
      call hsnint(nr,n,a,sx,method)
      go to 90
   50 continue
c
c evaluate analytic or finite difference hessian and check analytic
c hessian if requested (only if user-supplied analytic hessian
c routine d2fcn fills only lower triangular part and diagonal of a).
c
      if (iahflg .eq. 1) go to 60
c     if (iahflg .eq. 0)
c     then
         if (iagflg .eq. 1) call fstofd (nr, n, n, x, d1fcn, g, a, sx,
     &      rnf, wrk1, 3)
         if (iagflg .ne. 1) call sndofd (nr, n, x, fcn, f, a, sx, rnf,
     &      wrk1, wrk2)
         go to 80
c
c     else
   60    if (mod(msg/4,2).eq.0) go to 70
c        if (mod(msg/4, 2) .eq. 1)
c        then
            call d2fcn (nr, n, x, a)
            go to 80
c
c        else
   70       call heschk (nr, n, x, fcn, d1fcn, d2fcn, f, g, a, typsiz,
     &         sx,rnf,analtl,iagflg,udiag,wrk1,wrk2,msg,iprint,lounit)
c
c           heschk evaluates d2fcn and checks it against the finite
c           difference hessian which it calculates by calling fstofd
c           (if iagflg .eq. 1) or sndofd (otherwise).
c
            if (msg .lt. 0) return
   80 continue

   90 if(mod(msg/8,2).eq.0 .and. iprint.gt.0)then
        call result(nr,n,x,f,g,a,p,itncnt,1)
       end if
c
c
c iteration
c 
  100 itncnt=itncnt+1
c
c find perturbed local model hessian and its ll+ decomposition
c (skip this step if line search or dogstep techniques being used with
c secant updates.  cholesky decomposition l already obtained from
c secfac.)
c
      if(iexp.eq. 1 .and. method.ne. 3) go to 120
  110   call chlhsn(nr,n,a,epsm,sx,udiag)
  120 continue
c
c solve for newton step:  ap=-g
c
      do i=1,n
        wrk1(i)=-g(i)
      end do

      call lltslv(nr,n,a,p,wrk1)
c
c decide whether to accept newton step  xpls=x + p
c or to choose xpls by a global strategy.
c
      if (iagflg .ne. 0 .or. method .eq. 1) go to 140
      dltsav = dlt
      if (method .eq. 2) go to 140
      amusav = amu
      dlpsav = dltp
      phisav = phi
      phpsav = phip0
  140 if(method.eq. 1)
     &     call lnsrch(n,x,f,g,p,xpls,fpls,fcn,mxtake,iretcd,
     &     stepmx,steptl,sx,iprint,lounit)
      if(method.eq. 2)
     &     call dogdrv(nr,n,x,f,g,a,p,xpls,fpls,fcn,sx,stepmx,
     &     steptl,dlt,iretcd,mxtake,wrk0,wrk1,wrk2,wrk3,iprint,lounit)
      if(method.eq. 3)
     &     call hookdr(nr,n,x,f,g,a,udiag,p,xpls,fpls,fcn,sx,stepmx,
     &     steptl,dlt,iretcd,mxtake,amu,dltp,phi,phip0,wrk0,
     &     wrk1,wrk2,epsm,itncnt,iprint,lounit)
c
c if could not find satisfactory step and forward difference
c gradient was used, retry using central difference gradient.
c
      if (iretcd .ne. 1 .or. iagflg .ne. 0) go to 150
c     if (iretcd .eq. 1 .and. iagflg .eq. 0)
c     then
c
c        set iagflg for central differences
c
         iagflg = -1
         if(iprint.gt. 1)write(lounit,1060) itncnt

         call fstocd (n, x, fcn, sx, rnf, g)
         if (method .eq. 1) go to 120
         dlt = dltsav
         if (method .eq. 2) go to 120
         amu = amusav
         dltp = dlpsav
         phi = phisav
         phip0 = phpsav
         go to 110
c     end if
c
c calculate step for output
c
  150 continue

      do i = 1, n
         p(i) = xpls(i) - x(i)
      end do
c
c calculate gradient at xpls
c
      if (iagflg .eq. (-1)) go to 170
      if (iagflg .eq. 0) go to 180
c
c analytic gradient
c
      call d1fcn (n, xpls, gpls)
      go to 190
c
c central difference gradient
c
  170 call fstocd (n, xpls, fcn, sx, rnf, gpls)
      go to 190
c
c forward difference gradient
c
  180 call fstofd (1, 1, n, xpls, fcn, fpls, gpls, sx, rnf, wrk, 1)
  190 continue
c
c check whether stopping criteria satisfied
c
      call optstp(n,xpls,fpls,gpls,x,itncnt,icscmx,
     &            itrmcd,gradtl,steptl,sx,fscale,itnlim,iretcd,mxtake,
     &            iprint,lounit,msg)
      if(itrmcd.ne.0) go to 240
c
c evaluate hessian at xpls
c
      if(iexp.eq.0) go to 200

      if(method.eq. 3)
     &     call secunf(nr,n,x,g,a,udiag,xpls,gpls,epsm,itncnt,rnf,
     &     iagflg,noupdt,wrk1,wrk2,wrk3)
      if(method.ne. 3)
     &     call secfac(nr,n,x,g,a,xpls,gpls,epsm,itncnt,rnf,iagflg,
     &     noupdt,wrk0,wrk1,wrk2,wrk3)
      go to 220
  200 if(iahflg.eq. 1) go to 210
      if(iagflg.eq. 1)
     &     call fstofd(nr,n,n,xpls,d1fcn,gpls,a,sx,rnf,wrk1,3)

      if(iagflg.ne. 1) then
        call sndofd(nr,n,xpls,fcn,fpls,a,sx,rnf,wrk1,wrk2)
      end if

      go to 220
  210 call d2fcn(nr,n,xpls,a)
  220 continue
      if(mod(msg/16,2).eq. 1 .and. iprint.gt.0)then
        call result(nr,n,xpls,fpls,gpls,a,p,itncnt,1)
      end if
c
c x <-- xpls  and  g <-- gpls  and  f <-- fpls
c
      f=fpls
      do i=1,n
        x(i)=xpls(i)
        g(i)=gpls(i)
      end do
      go to 100
c
c termination
c reset xpls,fpls,gpls,  if previous iterate solution
c
  240 if(itrmcd.ne. 3) go to 270
  250 continue
      fpls=f
      do i=1,n
        xpls(i)=x(i)
        gpls(i)=g(i)
      end do
c
c print results
c
  270 continue
      if(mod(msg/8,2).eq.0 .and. iprint.gt.0)then
        call result(nr,n,xpls,fpls,gpls,a,p,itncnt,0)
      end if
      msg=0
      return

 1000 format(14h optdrv       ,5(g16.8,3x))
 1010 format(20h0optdrv    typical x)
 1020 format(40h optdrv    diagonal scaling matrix for x)
 1030 format(22h optdrv    typical f =,g16.8)
 1040 format(40h0optdrv    number of good digits in fcn=,i5/
     &       27h optdrv    gradient flag  =,i5,18h   (=1 if gradient,
     &       ' supplied)'/
     &       27h optdrv    hessian flag   =,i5,'   (=1 if hessian',
     &       ' supplied)'/
     &       27h optdrv    expense flag   =,i5, 9h   (=1 if,
     &       ' function expensive)'/
     &       ' optdrv    method to use=',i5/
     &       '           1=line search'/
     &       '           2=double dogleg'/
     &       '           3=more-hebdon'/
     &       27h optdrv    iteration limit=,i5/
     &       27h optdrv    machine epsilon=,g16.8)
 1050 format(33h0optdrv    maximum step size    =,g16.8/
     &       33h optdrv    step tolerance       =,g16.8/
     &       33h optdrv    gradient tolerance   =,g16.8/
     &       33h optdrv    trust region radius  =,g16.8/
     &       33h optdrv    rel noise in fcn     =,g16.8/
     &       33h optdrv    anal-fd tolerance    =,g16.8)
 1060 format(52h optdrv    shift from forward to central differences,
     &   14h in iteration , i5)
      end
      subroutine optif0 ( nr, n, x, fcn, xpls, fpls, gpls, itrmcd, 
     &  a, wrk )

c*********************************************************************72
c
cc OPTIF0 provides the simplest interface to the optimization package.
c
c  Discussion:
c
c    provide simplest interface to minimization package.
c    user has no control over options.
c
c  Modified:
c
c    09 January 2009
c
c  Author:
c
c    Robert Schnabel, John Koontz, Barry Weiss
c
c  Reference:
c
c    John Dennis, Robert Schnabel,
c    Numerical Methods for Unconstrained Optimization 
c    and Nonlinear Equations,
c    SIAM, 1996,
c    ISBN13: 978-0-898713-64-0,
c    LC: QA402.5.D44.
c
c    Robert Schnabel, John Koontz, Barry Weiss,
c    A modular system of algorithms for unconstrained minimization,
c    Technical Report CU-CS-240-82, 
c    Computer Science Department,
c    University of Colorado at Boulder, 1982.
c
c  Parameters:
c
c nr           --> row dimension of matrix
c                  a positive integer specifying the row dimension
c                  of the matrices a and wrk in the user's calling
c                  program.  nr must be greater or equal to n.
c
c n            --> dimension of problem
c                  a positive integer specifying the order or dimension
c                  of the problem.  the program will abort if
c                  n.le.0.  the program is inefficient for n.eq.1.
c                  for n.eq.1, do not call this routine, but rather
c                  optif9.
c
c x(n)        <--> initial estimate of minimum on input.
c                  on return, x contains x(k-1), the next to last
c                  iterate.
c
c fcn          --> name of routine to evaluate minimization function.
c                  fcn must be declared external in the calling
c                  routine, and must have the form
c                  subroutine fcn(n,x,f) where x is a vector of
c                  length n.  the subroutine must not alter x or n.
c
c xpls(n)     <--  local minimum
c                  an n dimensional array containing the best
c                  approximation to the local minimum on return
c                  if the algorithm has converged.
c
c fpls        <--  function value at local minimum xpls
c
c gpls(n)     <--  gradient at local minimum xpls
c
c itrmcd      <--  termination code
c                  0=bad input (see msg).
c                  1=relative gradient is close to zero.  current
c                  iterate is probably solution.
c                  2=successive iterates within tolerance. current
c                  iterate is probably solution.
c                  3=last global step failed to locate a point
c                  lower than xpls.  either xpls is an
c                  approximate local minimum of the function,
c                  the function is too nonlinear for this program,
c                  or steptl is too large.
c                  4=iteration limit exceeded.
c                  5=maximum step size stepmx exceeded 5 consecutive
c                  times.  either the function is unbounded below,
c                  becomes asymptotic to a finite value from above
c                  in some direction, or stepmx is too small.
c
c a(n,n)       --> workspace
c                  a matrix used to store the hessian and its
c                  cholesky decomposition.  the user must declare
c                  this array of dimensions nr by nc where nr
c                  and nc are both at least n.  the hessian is
c                  not evaluated at the last iterate so on return
c                  the information in a is old.
c
c wrk(n,9)     --> workspace
c                  workspace required by the subroutine.
c
      implicit double precision ( a-h, o-z )

      dimension x(n),xpls(n),gpls(n)
      dimension a(nr,1),wrk(nr,1)
      external fcn,d1fcn,d2fcn
c
c equivalence wrk(n,1) = udiag(n)
c             wrk(n,2) = g(n)
c             wrk(n,3) = p(n)
c             wrk(n,4) = typsiz(n)
c             wrk(n,5) = sx(n)
c             wrk(n,6) = wrk0(n)
c             wrk(n,7) = wrk1(n)
c             wrk(n,8) = wrk2(n)
c             wrk(n,9) = wrk3(n)
c
      call dfault(n,wrk(1,4),fscale,method,iexp,msg,ndigit,
     &     itnlim,iagflg,iahflg,iprint,lounit,dlt,gradtl,stepmx,steptl)

      call optdrv(nr,n,x,fcn,d1fcn,d2fcn,wrk(1,4),fscale,
     &     method,iexp,msg,ndigit,itnlim,iagflg,iahflg,iprint,lounit,
     &     dlt,gradtl,stepmx,steptl,
     &     xpls,fpls,gpls,itrmcd,
     &     a,wrk(1,1),wrk(1,2),wrk(1,3),wrk(1,5),wrk(1,6),
     &     wrk(1,7),wrk(1,8),wrk(1,9))

      return
      end
      subroutine optif9(nr,n,x,fcn,d1fn,d2fn,typsiz,fscale,
     &     method,iexp,msg,ndigit,itnlim,iagflg,iahflg,iprint,lounit,
     &     dlt,gradtl,stepmx,steptl,xpls,fpls,gpls,itrmcd,a,wrk)

c*********************************************************************72
c
cc OPTIF9 provides a complete interface to the minimization package.
c
c  Discussion:
c
c    provide complete interface to minimization package.
c    user has full control over options.
c
c  Modified:
c
c    09 January 2009
c
c  Author:
c
c    Robert Schnabel, John Koontz, Barry Weiss
c
c  Reference:
c
c    John Dennis, Robert Schnabel,
c    Numerical Methods for Unconstrained Optimization 
c    and Nonlinear Equations,
c    SIAM, 1996,
c    ISBN13: 978-0-898713-64-0,
c    LC: QA402.5.D44.
c
c    Robert Schnabel, John Koontz, Barry Weiss,
c    A modular system of algorithms for unconstrained minimization,
c    Technical Report CU-CS-240-82, 
c    Computer Science Department,
c    University of Colorado at Boulder, 1982.
c
c  Parameters:
c
c nr           --> row dimension of matrix
c                  the row dimension of the matrices a and wrk
c                  as declared in the user calling routine.
c                  nr must be greater than or equal to n.
c
c n            --> dimension of problem, the number of variables.
c                  the program will abort if n.le.0.
c
c x(n)        <--> on entry: estimate to a root of fcn
c                  on return, the next to the last iterate x(k-1).
c
c fcn          --> name of subroutine to evaluate optimization function
c                  must be declared external in calling routine
c                  the calling statement must be
c                  subroutine fcn(n,x,f) where x is an n vector.
c                  fcn must not change the value of x or n.
c
c d1fn         --> name of subroutine to evaluate gradient
c                  of fcn.  must be declared external in calling routine
c                  the calling statement must be
c                  subroutine d1fn(n,x,g) where x and g are n vectors.
c                  the subroutine must not alter x or n.  on return,
c                  g must be the gradient of f, that is,
c                  g(i)=df/dx(i).  if no analytic gradient routine is
c                  used, the user must use the dummy name d1fcn.
c                  if an analytic gradient is supplied, the program
c                  will automatically compare the gradient with one
c                  generated by finite differences unless msg is
c                  appropriately set.
c
c d2fn         --> name of subroutine to evaluate hessian of
c                  fcn.  must be declared external in calling routine
c                  the subroutine statement must be
c                  subroutine d2fn(nr,n,x,h) where x is a vector
c                  of length n, and h is an n by n matrix stored
c                  with a row dimension of nr.  the subroutine must
c                  not change n, nr, or x.  if no analytic hessian
c                  is to be supplied, the user must use the name
c                  d2fcn.  otherwise, d2fn should store the hessian
c                  in h, that is, h(i,j)=d d(f)/d x(i) d x(j).
c                  if an analytic hessian is supplied, the program
c                  will automatically compare it to one computed by
c                  finite differences, unless msg is appropriately set
c
c typsiz(n)    --> typical size for each component of x
c                  each entry should be positive.  this information
c                  is used to scale the information.  if variables
c                  x(i) and x(j) are typically of much different
c                  magnitudes, the scaling of typsiz may help to
c                  solve the problem.
c
c fscale       --> estimate of scale of objective function
c                  a positive scalar estimating the magnitude of
c                  the optimizing function near the minimum x*.
c                  if too large a value is supplied, the program
c                  may terminate prematurely.
c
c method       --> algorithm to use to solve minimization problem
c                  =1 line search
c                  =2 double dogleg
c                  =3 more-hebdon
c
c iexp         --> =1 if optimization function fcn is expensive to
c                  evaluate, =0 otherwise.  if set then hessian will
c                  be evaluated by secant update instead of
c                  analytically or by finite differences
c
c msg         <--> on input:  (.gt.0) message to inhibit certain
c                  automatic checks.  set msg to a value between
c                  0 and 31 by adding the value of each option
c                  below that is to be indicated.
c                  0=no options.
c                  2=do not check analytic gradient routine d1fn
c                  against a finite difference estimate.
c                  4=do not check analytic hessian routine d2fn
c                  against a finite difference estimate.
c                  8=suppress printing of input state and final
c                  results.
c                  16=print intermediate results.
c
c                  on output: (.lt.0) error code; =0 no error
c
c                  -1 = illegal dimension, n.le.0.
c                  -3 = illegal tolerance on gradient, gradtl.lt.0
c                  -4 = iteration limit itnlim.le.0
c                  -5 = no good digits in optimization function,
c                  ndigit=0.
c                  -21 = probable coding error in user gradient
c                  routine d1fn.  analytic and finite difference
c                  gradients do not agree.
c                  -22 = probable coding error in user hessian
c                  routine d2fn.  analytic and finite difference
c                  hessians do not agree.
c
c ndigit       --> number of good digits in optimization function fcn
c                  negative means function has maximum possible
c                  good digits.  0 is not allowed.  positive
c                  value is appropriate when the function is only
c                  computed approximately via quadrature, iteration
c                  or so on.
c
c itnlim       --> maximum number of allowable iterations
c                  nonpositive values are illegal.  this number
c                  limits the number of iterations to be taken.
c
c iagflg       --> =0 if no analytic gradient supplied, d1fcn to be
c                  used.
c                  =1 if analytic gradient routine supplied.
c
c iahflg       --> =0 if no analytic hessian supplied, d2fcn to be
c                  used.
c                  =1 if analytic hessian routine supplied.
c
c  iprint      --> amount of output desired
c                  0=none,
c                  1=starting, stopping, and errors
c                  2=all that, plus intermediate information
c
c lounit       --> device to which to send output
c                  0 means no output.
c                  6 means output to tty or lineprinter.
c
c dlt          --> trust region radius
c                  not needed when line search used (method=1).
c                  otherwise, if no estimate is known, set dlt=-1.
c                  if an estimate of the maximum reasonable scaled
c                  step is known, set dlt to this value.
c
c gradtl       --> tolerance at which gradient considered close
c                  enough to zero to terminate algorithm
c                  the scaled gradient is a measure of the
c                  relative change in f in each direction x(i),
c                  divided by the relative change in x(i).
c
c stepmx       --> maximum allowable step size
c                  stepmx is used to keep the function from
c                  overflowing, or the arguments from leaving the
c                  area of interest, and to detect divergence.
c
c steptl       --> relative step size at which successive iterates
c                  considered close enough to terminate algorithm
c                  steptl can be set to 10**(-p) where p is the
c                  number of digits of agreement between successive
c                  iterates which the user considers to mean
c                  satisfactory convergence.
c
c xpls(n)     <--  on exit:  xpls is local minimum
c                  if the algorithm has converged.
c
c fpls        <--  on exit:  function value at solution, xpls
c
c gpls(n)     <--  on exit:  gradient at solution xpls
c
c itrmcd      <--  termination code
c                  0=erroneous input detected.  see msg.
c                  1=relative gradient is close to zero.
c                  current iterate is probably solution.
c                  2=successive iterates within tolerance.
c                  current iterate is probably solution.
c                  3=last global step failed to locate a
c                  point lower than xpls.  either xpls is a
c                  local minimum, or the function is too nonlinear,
c                  or steptl is too large.
c                  4=iteration limit exceeded.
c                  5=stepmx maximum stepsize exceeded 5 times
c                  in a row.  either the function is unbounded below,
c                  or becomes asymptotic to a finite value from
c                  above, or stepmx is too small.
c
c a(n,n)       --> workspace for hessian (or estimate)
c                  and its cholesky decomposition
c                  dimensioned (nr,n) in calling program.  the hessian
c                  is not evaluated at the last iterate, so information
c                  in a is slightly out of date on return.
c
c wrk(n,8)     --> workspace
c                  workspace needed by the program.
c
      implicit double precision ( a-h, o-z )

      dimension x(n),xpls(n),gpls(n),typsiz(n)
      dimension a(nr,1),wrk(nr,1)
      external fcn,d1fn,d2fn
c
c equivalence wrk(n,1) = udiag(n)
c             wrk(n,2) = g(n)
c             wrk(n,3) = p(n)
c             wrk(n,4) = sx(n)
c             wrk(n,5) = wrk0(n)
c             wrk(n,6) = wrk1(n)
c             wrk(n,7) = wrk2(n)
c             wrk(n,8) = wrk3(n)
c
      call optdrv(nr,n,x,fcn,d1fn,d2fn,typsiz,fscale,
     & method,iexp,msg,ndigit,itnlim,iagflg,iahflg,iprint,lounit,
     & dlt,gradtl,stepmx,steptl,xpls,fpls,gpls,itrmcd,
     & a,wrk(1,1),wrk(1,2),wrk(1,3),wrk(1,4),wrk(1,5),
     & wrk(1,6),wrk(1,7),wrk(1,8))

      return
      end
      subroutine optstp(n,xpls,fpls,gpls,x,itncnt,icscmx,
     & itrmcd,gradtl,steptl,sx,fscale,itnlim,iretcd,mxtake,iprint,
     & lounit,msg)

c*********************************************************************72
c
cc OPTSTP checks the optimization stopping criteria.
c
c  Discussion:
c
c    find whether the algorithm should terminate, due to any
c    of the following:
c    1) problem solved within user tolerance
c    2) convergence within user tolerance
c    3) iteration limit reached
c    4) divergence or too restrictive maximum step (stepmx) suspected
c
c  Modified:
c
c    09 January 2009
c
c  Author:
c
c    Robert Schnabel, John Koontz, Barry Weiss
c
c  Reference:
c
c    John Dennis, Robert Schnabel,
c    Numerical Methods for Unconstrained Optimization 
c    and Nonlinear Equations,
c    SIAM, 1996,
c    ISBN13: 978-0-898713-64-0,
c    LC: QA402.5.D44.
c
c    Robert Schnabel, John Koontz, Barry Weiss,
c    A modular system of algorithms for unconstrained minimization,
c    Technical Report CU-CS-240-82, 
c    Computer Science Department,
c    University of Colorado at Boulder, 1982.
c
c  Parameters:
c
c n            --> dimension of problem
c xpls(n)      --> new iterate x(k)
c fpls         --> function value at new iterate f(xpls)
c gpls(n)      --> gradient at new iterate, g(xpls), or approximate
c x(n)         --> old iterate x(k-1)
c itncnt       --> current iteration k
c icscmx      <--> number consecutive steps .ge. stepmx
c                  (retain value between successive calls)
c itrmcd      <--  termination code
c gradtl       --> tolerance at which relative gradient considered close
c                  enough to zero to terminate algorithm
c steptl       --> relative step size at which successive iterates
c                  considered close enough to terminate algorithm
c sx(n)        --> diagonal scaling matrix for x
c fscale       --> estimate of scale of objective function
c itnlim       --> maximum number of allowable iterations
c iretcd       --> return code
c mxtake       --> boolean flag indicating step of maximum length used
c iprint       --> amount of output desired.
c lounit       --> device to which to send output
c msg          --> if msg includes a term 8, suppress output
c
      implicit double precision ( a-h, o-z )

      integer n,itncnt,icscmx,itrmcd,itnlim
      dimension sx(n)
      dimension xpls(n),gpls(n),x(n)
      logical mxtake

      itrmcd=0
c
c last global step failed to locate a point lower than x
      if(iretcd.ne. 1) go to 10
c     if(iretcd.eq. 1)
c     then
        jtrmcd=3
        go to 50
c     end if
   10 continue
c
c find direction in which relative gradient maximum.
c check whether within tolerance
c
      d=max(abs(fpls),fscale)
      rgx = 0.0D+00
      do i=1,n
        relgrd=abs(gpls(i))*max(abs(xpls(i)), 1.0D+00 /sx(i))/d
        rgx=max(rgx,relgrd)
      end do

      jtrmcd=1
      if(rgx.le.gradtl) go to 50

      if(itncnt.eq.0) return
c
c find direction in which relative stepsize maximum
c check whether within tolerance.
c
      rsx=0.0D+00
      do i=1,n
        relstp=abs(xpls(i)-x(i))/max(abs(xpls(i)), 1.0D+00 /sx(i))
        rsx=max(rsx,relstp)
      end do
      jtrmcd=2
      if(rsx.le.steptl) go to 50
c
c check iteration limit
c
      jtrmcd=4
      if(itncnt.ge.itnlim) go to 50
c
c check number of consecutive steps.ne.stepmx
c
      if(mxtake) go to 40
c     if(.not.mxtake)
c     then
        icscmx=0
        return
c     else
   40   continue
        if (mod(msg/8,2).eq. 0 .and. iprint.gt. 1)write(lounit,1000)
        icscmx=icscmx+1
        if(icscmx.lt. 5) return
        jtrmcd=5
c     end if
c
c
c print termination code
c
   50 itrmcd=jtrmcd
      if (mod(msg/8,2) .eq. 0) go to(60,70,80,90,100), itrmcd
      go to 110
   60 if(iprint.gt. 1)write(lounit,1010)
      go to 110
   70 if(iprint.gt. 1)write(lounit,1020)
      go to 110
   80 if(iprint.gt. 1)write(lounit,1030)
      go to 110
   90 if(iprint.gt. 1)write(lounit,1040)
      go to 110
  100 if(iprint.gt. 1)write(lounit,1050)

  110 return

 1000 format(48h0optstp    step of maximum length (stepmx) taken)
 1010 format(43h0optstp    relative gradient close to zero./
     &       48h optstp    current iterate is probably solution.)
 1020 format(48h0optstp    successive iterates within tolerance./
     &       48h optstp    current iterate is probably solution.)
 1030 format(52h0optstp    last global step failed to locate a point,
     &       14h lower than x./
     &       51h optstp    either x is an approximate local minimum,
     &       17h of the function,/
     &       50h optstp    the function is too non-linear for this,
     &       11h algorithm,/
     &       34h optstp    or steptl is too large.)
 1040 format(36h0optstp    iteration limit exceeded./
     &       28h optstp    algorithm failed.)
 1050 format(39h0optstp    maximum step size exceeded 5,
     &       19h consecutive times./
     &       50h optstp    either the function is unbounded below,/
     &       47h optstp    becomes asymptotic to a finite value,
     &       30h from above in some direction,/
     &       33h optstp    or stepmx is too small)
      end
      subroutine qraux1(nr,n,r,i)

c*********************************************************************72
c
cc QRAUX1 interchanges two rows of an upper Hessenberg matrix.
c
c  Discussion:
c
c    interchange rows i,i+1 of the upper hessenberg matrix r,
c    columns i to n
c
c  Modified:
c
c    09 January 2009
c
c  Author:
c
c    Robert Schnabel, John Koontz, Barry Weiss
c
c  Reference:
c
c    John Dennis, Robert Schnabel,
c    Numerical Methods for Unconstrained Optimization 
c    and Nonlinear Equations,
c    SIAM, 1996,
c    ISBN13: 978-0-898713-64-0,
c    LC: QA402.5.D44.
c
c    Robert Schnabel, John Koontz, Barry Weiss,
c    A modular system of algorithms for unconstrained minimization,
c    Technical Report CU-CS-240-82, 
c    Computer Science Department,
c    University of Colorado at Boulder, 1982.
c
c  Parameters:
c
c nr           --> row dimension of matrix
c n            --> dimension of matrix
c r(n,n)      <--> upper hessenberg matrix
c i            --> index of row to interchange (i.lt.n)
c
      implicit double precision ( a-h, o-z )

      dimension r(nr,1)

      do j=i,n
        tmp=r(i,j)
        r(i,j)=r(i+1,j)
        r(i+1,j)=tmp
      end do

      return
      end
      subroutine qraux2(nr,n,r,i,a,b)

c*********************************************************************72
c
cc QRAUX2 premultiplies a matrix by a Jacobi rotation.
c
c  Discussion:
c
c    pre-multiply r by the jacobi rotation j(i,i+1,a,b)
c
c  Modified:
c
c    09 January 2009
c
c  Author:
c
c    Robert Schnabel, John Koontz, Barry Weiss
c
c  Reference:
c
c    John Dennis, Robert Schnabel,
c    Numerical Methods for Unconstrained Optimization 
c    and Nonlinear Equations,
c    SIAM, 1996,
c    ISBN13: 978-0-898713-64-0,
c    LC: QA402.5.D44.
c
c    Robert Schnabel, John Koontz, Barry Weiss,
c    A modular system of algorithms for unconstrained minimization,
c    Technical Report CU-CS-240-82, 
c    Computer Science Department,
c    University of Colorado at Boulder, 1982.
c
c  Parameters:
c
c nr           --> row dimension of matrix
c n            --> dimension of matrix
c r(n,n)      <--> upper hessenberg matrix
c i            --> index of row
c a            --> scalar
c b            --> scalar
c
      implicit double precision ( a-h, o-z )

      dimension r(nr,1)

      den=sqrt(a*a + b*b)
      c=a/den
      s=b/den
      do j=i,n
        y=r(i,j)
        z=r(i+1,j)
        r(i,j)=c*y - s*z
        r(i+1,j)=s*y + c*z
      end do

      return
      end
      subroutine qrupdt(nr,n,a,u,v)

c*********************************************************************72
c
cc QRUPDT updates a QR factorization.
c
c  Discussion:
c
c    find an orthogonal (n*n) matrix (q*) and an upper triangular (n*n)
c    matrix (r*) such that (q*)(r*)=r+u(v+)
c
c  Modified:
c
c    09 January 2009
c
c  Author:
c
c    Robert Schnabel, John Koontz, Barry Weiss
c
c  Reference:
c
c    John Dennis, Robert Schnabel,
c    Numerical Methods for Unconstrained Optimization 
c    and Nonlinear Equations,
c    SIAM, 1996,
c    ISBN13: 978-0-898713-64-0,
c    LC: QA402.5.D44.
c
c    Robert Schnabel, John Koontz, Barry Weiss,
c    A modular system of algorithms for unconstrained minimization,
c    Technical Report CU-CS-240-82, 
c    Computer Science Department,
c    University of Colorado at Boulder, 1982.
c
c  Parameters:
c
c nr           --> row dimension of matrix
c n            --> dimension of problem
c a(n,n)      <--> on input:  contains r
c                  on output: contains (r*)
c u(n)         --> vector
c v(n)         --> vector
c
      implicit double precision ( a-h, o-z )

      dimension a(nr,1)
      dimension u(n),v(n)
c
c determine last non-zero in u(.)
c
      k=n
   10 continue

      if(u(k).ne. 0.0D+00 .or. k.eq. 1) go to 20
c     if(u(k).eq. 0.0D+00 .and. k.gt. 1)
c     then
        k=k-1
        go to 10
c     end if
c
c (k-1) jacobi rotations transform
c     r + u(v+) --> (r*) + (u(1)*e1)(v+)
c which is upper hessenberg
c
   20 if(k.le. 1) go to 50
        km1=k-1
        do 40 ii=1,km1
          i=km1-ii+1
          if(u(i).ne. 0.0D+00 ) go to 30
            call qraux1(nr,n,a,i)
            u(i)=u(i+1)
            go to 40
   30       call qraux2(nr,n,a,i,u(i),-u(i+1))
            u(i)=sqrt(u(i)*u(i) + u(i+1)*u(i+1))
   40   continue
c     end if
c
c  r <-- r + (u(1)*e1)(v+)
c
   50 do j=1,n
        a(1,j)=a(1,j) +u(1)*v(j)
      end do
c
c  (k-1) jacobi rotations transform upper hessenberg r
c  to upper triangular (r*)
c
      if(k.le. 1) go to 90
        km1=k-1
        do 80 i=1,km1
          if(a(i,i).ne. 0.0D+00 ) go to 70
            icopy=i
            call qraux1(nr,n,a,icopy)
            go to 80
   70       t1=a(i,i)
            t2=-a(i+1,i)
            icopy=i
            call qraux2(nr,n,a,icopy,t1,t2)
   80   continue
   90 return
      end
      function r8_epsilon ( )

c*********************************************************************72
c
cc R8_EPSILON returns the R8 roundoff unit.
c
c  Discussion:
c
c    The roundoff unit is a number R which is a power of 2 with the
c    property that, to the precision of the computer's arithmetic,
c      1 .lt. 1 + R
c    but
c      1 = ( 1 + R / 2 )
c
c    FORTRAN90 provides the superior library routine
c
c      EPSILON ( X )
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    06 March 2006
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, double precision R8_EPSILON, the R8 roundoff unit.
c
      implicit none

      double precision r8
      double precision r8_epsilon
      double precision r8_test

      r8 = 1.0D+00
      r8_test = 1.0D+00 + ( r8 / 2.0D+00 )

10    continue

      if ( 1.0D+00 .lt. r8_test ) then
        r8 = r8 / 2.0D+00
        r8_test = 1.0D+00 + ( r8 / 2.0D+00 )
        go to 10
      end if

      r8_epsilon = r8

      return
      end
      function r8vec_dot ( n, v1, v2 )

c*********************************************************************72
c
cc R8VEC_DOT finds the dot product of a pair of R8VEC's.
c
c  Discussion:
c
c    An R8VEC is a vector of R8 values.
c
c    In FORTRAN90, the system routine DOT_PRODUCT should be called
c    directly.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    27 May 2008
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the dimension of the vectors.
c
c    Input, double precision V1(N), V2(N), the vectors.
c
c    Output, double precision R8VEC_DOT, the dot product.
c
      implicit none

      integer n

      integer i
      double precision r8vec_dot
      double precision v1(n)
      double precision v2(n)
      double precision value

      value = 0.0D+00
      do i = 1, n
        value = value + v1(i) * v2(i)
      end do

      r8vec_dot = value

      return
      end
      function r8vec_norm_l2 ( n, a )

c*********************************************************************72
c
cc R8VEC_NORM_L2 returns the L2 norm of an R8VEC.
c
c  Discussion:
c
c    An R8VEC is a vector of R8 values.
c
c    The vector L2 norm is defined as:
c
c      R8VEC_NORM_L2 = sqrt ( sum ( 1 <= I <= N ) A(I)**2 ).
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    27 May 2008
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the number of entries in A.
c
c    Input, double precision A(N), the vector whose L2 norm is desired.
c
c    Output, double precision R8VEC_NORM_L2, the L2 norm of A.
c
      implicit none

      integer n

      double precision a(n)
      integer i
      double precision r8vec_norm_l2
      double precision value

      value = 0.0D+00
      do i = 1, n
        value = value + a(i) * a(i)
      end do
      value = sqrt ( value )

      r8vec_norm_l2 = value

      return
      end
      subroutine result(nr,n,x,f,g,a,p,itncnt,iflg)

c*********************************************************************72
c
cc RESULT prints information from the optimization procedure.
c
c  Discussion:
c
c    print information
c
c  Modified:
c
c    09 January 2009
c
c  Author:
c
c    Robert Schnabel, John Koontz, Barry Weiss
c
c  Reference:
c
c    John Dennis, Robert Schnabel,
c    Numerical Methods for Unconstrained Optimization 
c    and Nonlinear Equations,
c    SIAM, 1996,
c    ISBN13: 978-0-898713-64-0,
c    LC: QA402.5.D44.
c
c    Robert Schnabel, John Koontz, Barry Weiss,
c    A modular system of algorithms for unconstrained minimization,
c    Technical Report CU-CS-240-82, 
c    Computer Science Department,
c    University of Colorado at Boulder, 1982.
c
c  Parameters:
c
c nr           --> row dimension of matrix
c n            --> dimension of problem
c x(n)         --> iterate x(k)
c f            --> function value at x(k)
c g(n)         --> gradient at x(k)
c a(n,n)       --> hessian at x(k)
c p(n)         --> step taken
c itncnt       --> iteration number k
c iflg         --> flag controlling info to print
c
      implicit double precision ( a-h, o-z )

      dimension a(nr,1)
      dimension g(n)
      dimension p(n)
      dimension x(n)

      write(*,*)' '
      write(*,*)'RESULT - Iterate ',itncnt
      if(iflg.ne.0)then
        write(*,*)' '
        write(*,*)'Step'
        write(*,*)' '
        write(*,'(5g14.6)') (p(i),i=1,n)
      end if
      write(*,*)' '
      write(*,*)'X'
      write(*,*)' '
      write(*,'(5g14.6)') (x(i),i=1,n)
      write(*,*)'Function value = ',f
      write(*,*)' '
      write(*,*)'Gradient vector'
      write(*,*)' '
      write(*,'(5g14.6)') (g(i),i=1,n)
      if(iflg.ne.0)then
        write(*,*)' '
        write(*,*)'Hessian matrix'
        write(*,*)' '
        do i=1,n
          write(*,'(5g14.6)') (a(i,j),j=1,i),(a(k,i),k=i+1,n)
        end do
      end if

      return
      end
      subroutine sclmul(n,s,v,z)

c*********************************************************************72
c
cc SCLMUL multiplies a vector by a scalar.
c
c  Discussion:
c
c    multiply vector by scalar
c    result vector may be operand vector
c
c  Modified:
c
c    09 January 2009
c
c  Author:
c
c    Robert Schnabel, John Koontz, Barry Weiss
c
c  Reference:
c
c    John Dennis, Robert Schnabel,
c    Numerical Methods for Unconstrained Optimization 
c    and Nonlinear Equations,
c    SIAM, 1996,
c    ISBN13: 978-0-898713-64-0,
c    LC: QA402.5.D44.
c
c    Robert Schnabel, John Koontz, Barry Weiss,
c    A modular system of algorithms for unconstrained minimization,
c    Technical Report CU-CS-240-82, 
c    Computer Science Department,
c    University of Colorado at Boulder, 1982.
c
c  Parameters:
c
c n            --> dimension of vectors
c s            --> scalar
c v(n)         --> operand vector
c z(n)        <--  result vector
c
      implicit double precision ( a-h, o-z )

      dimension v(n),z(n)

      do i=1,n
        z(i)=s*v(i)
      end do

      return
      end
      subroutine secfac(nr,n,x,g,a,xpls,gpls,epsm,itncnt,rnf,
     &     iagflg,noupdt,s,y,u,w)

c*********************************************************************72
c
cc SECFAC updates the Hessian matrix by the BFGS factored method.
c
c  Discussion:
c
c    update hessian by the bfgs factored method
c
c  Modified:
c
c    09 January 2009
c
c  Author:
c
c    Robert Schnabel, John Koontz, Barry Weiss
c
c  Reference:
c
c    John Dennis, Robert Schnabel,
c    Numerical Methods for Unconstrained Optimization 
c    and Nonlinear Equations,
c    SIAM, 1996,
c    ISBN13: 978-0-898713-64-0,
c    LC: QA402.5.D44.
c
c    Robert Schnabel, John Koontz, Barry Weiss,
c    A modular system of algorithms for unconstrained minimization,
c    Technical Report CU-CS-240-82, 
c    Computer Science Department,
c    University of Colorado at Boulder, 1982.
c
c  Parameters:
c
c nr           --> row dimension of matrix
c n            --> dimension of problem
c x(n)         --> old iterate, x(k-1)
c g(n)         --> gradient or approximate at old iterate
c a(n,n)      <--> on entry: cholesky decomposition of hessian in
c                    lower part and diagonal.
c                  on exit:  updated cholesky decomposition of hessian
c                    in lower triangular part and diagonal
c xpls(n)      --> new iterate, x(k)
c gpls(n)      --> gradient or approximate at new iterate
c epsm         --> machine epsilon
c itncnt       --> iteration count
c rnf          --> relative noise in optimization function fcn
c iagflg       --> =1 if analytic gradient supplied, =0 itherwise
c noupdt      <--> boolean: no update yet
c                  (retain value between successive calls)
c s(n)         --> workspace
c y(n)         --> workspace
c u(n)         --> workspace
c w(n)         --> workspace
c
      implicit double precision ( a-h, o-z )

      double precision r8vec_dot
      double precision r8vec_norm_l2
      dimension x(n),xpls(n),g(n),gpls(n)
      dimension a(nr,1)
      dimension s(n),y(n),u(n),w(n)
      logical noupdt,skpupd

      if(itncnt.eq. 1) noupdt=.true.
      do i=1,n
        s(i)=xpls(i)-x(i)
        y(i)=gpls(i)-g(i)
      end do

      den1 = r8vec_dot ( n, x, y )

      snorm2 = r8vec_norm_l2 ( n, s )
      ynrm2 = r8vec_norm_l2 ( n, y )

      if ( den1 .lt. sqrt(epsm)*snorm2*ynrm2) go to 160
        call mvmltu(nr,n,a,s,u)
        den2 = r8vec_dot ( n, u, u )
c
c  l <-- sqrt(den1/den2)*l
c
        alp=sqrt(den1/den2)
        if(.not.noupdt) go to 40
          do j=1,n
            u(j)=alp*u(j)
            do i=j,n
              a(i,j)=alp*a(i,j)
            end do
          end do
          noupdt=.false.
          den2=den1
          alp = 1.0D+00
   40   skpupd=.true.
c
c  w = l(l+)s = hs
c
        call mvmltl(nr,n,a,u,w)
        i=1
        if(iagflg.ne.0) go to 50
          reltol=sqrt(rnf)
          go to 60
   50     reltol=rnf
   60   if(i.gt.n .or. .not.skpupd) go to 80
          if(abs(y(i)-w(i)) .lt. reltol*max(abs(g(i)),abs(gpls(i))))
     &         go to 70
            skpupd=.false.
            go to 60
   70       i=i+1
            go to 60
   80   if(skpupd) go to 160
c
c  w=y-alp*l(l+)s
c
          do i=1,n
            w(i)=y(i)-alp*w(i)
          end do
c
c  alp=1/sqrt(den1*den2)
c
          alp=alp/den1
c
c  u=(l+)/sqrt(den1*den2) = (l+)s/sqrt((y+)s * (s+)l(l+)s)
c
          do i=1,n
            u(i)=alp*u(i)
          end do
c
c  copy l into upper triangular part.  zero l.
c
          if(n.eq. 1) go to 130
          do i=2,n
            im1=i-1
            do j=1,im1
              a(j,i)=a(i,j)
              a(i,j)= 0.0D+00
            end do
          end do
c
c  find q, (l+) such that  q(l+) = (l+) + u(w+)
c
  130     call qrupdt(nr,n,a,u,w)
c
c  upper triangular part and diagonal of a now contain updated
c  cholesky decomposition of hessian.  copy back to lower
c  triangular part.
c
          if(n.eq. 1) go to 160
          do i=2,n
            im1=i-1
            do j=1,im1
              a(i,j)=a(j,i)
            end do
          end do
c       end if
c     end if
  160 return
      end
      subroutine secunf(nr,n,x,g,a,udiag,xpls,gpls,epsm,itncnt,
     &     rnf,iagflg,noupdt,s,y,t)

c*********************************************************************72
c
cc SECUNF updates the Hessian by the BFGS unfactored method.
c
c  Discussion:
c
c    update hessian by the bfgs unfactored method
c
c  Modified:
c
c    09 January 2009
c
c  Author:
c
c    Robert Schnabel, John Koontz, Barry Weiss
c
c  Reference:
c
c    John Dennis, Robert Schnabel,
c    Numerical Methods for Unconstrained Optimization 
c    and Nonlinear Equations,
c    SIAM, 1996,
c    ISBN13: 978-0-898713-64-0,
c    LC: QA402.5.D44.
c
c    Robert Schnabel, John Koontz, Barry Weiss,
c    A modular system of algorithms for unconstrained minimization,
c    Technical Report CU-CS-240-82, 
c    Computer Science Department,
c    University of Colorado at Boulder, 1982.
c
c  Parameters:
c
c nr           --> row dimension of matrix
c n            --> dimension of problem
c x(n)         --> old iterate, x(k-1)
c g(n)         --> gradient or approximate at old iterate
c a(n,n)      <--> on entry: approximate hessian at old iterate
c                    in upper triangular part (and udiag)
c                  on exit:  updated approx hessian at new iterate
c                    in lower triangular part and diagonal
c                  (lower triangular part of symmetric matrix)
c udiag        --> on entry: diagonal of hessian
c xpls(n)      --> new iterate, x(k)
c gpls(n)      --> gradient or approximate at new iterate
c epsm         --> machine epsilon
c itncnt       --> iteration count
c rnf          --> relative noise in optimization function fcn
c iagflg       --> =1 if analytic gradient supplied, =0 otherwise
c noupdt      <--> boolean: no update yet
c                  (retain value between successive calls)
c s(n)         --> workspace
c y(n)         --> workspace
c t(n)         --> workspace
c
      implicit double precision ( a-h, o-z )

      double precision r8vec_dot
      double precision r8vec_norm_l2
      dimension x(n),g(n),xpls(n),gpls(n)
      dimension a(nr,1)
      dimension udiag(n)
      dimension s(n),y(n),t(n)
      logical noupdt,skpupd
c
c  copy hessian in upper triangular part and udiag to
c  lower triangular part and diagonal
c
      do 20 j=1,n
        a(j,j)=udiag(j)
        if(j.eq.n) go to 20
        jp1=j+1
        do i=jp1,n
          a(i,j)=a(j,i)
        end do
   20 continue

      if(itncnt.eq. 1) noupdt=.true.

      do i=1,n
        s(i)=xpls(i)-x(i)
        y(i)=gpls(i)-g(i)
      end do

      den1 = r8vec_dot ( n, s, y )
      snorm2 = r8vec_norm_l2 ( n, s )
      ynrm2 = r8vec_norm_l2 ( n, y )
      if(den1 .lt.sqrt(epsm)*snorm2*ynrm2) go to 110
c     if(den1 .ge.sqrt(epsm)*snorm2*ynrm2)
c     then
        call mvmlts(nr,n,a,s,t)
        den2 = r8vec_dot ( n, s, t )
        if(.not. noupdt) go to 60
c       if(noupdt)
c       then
c
c         h <-- ((s+)y/(s+)hs)h
c
          gam=den1/den2
          den2=gam*den2
          do j=1,n
            t(j)=gam*t(j)
            do i=j,n
              a(i,j)=gam*a(i,j)
            end do
          end do
          noupdt=.false.
c       end if
   60   skpupd=.true.
c
c  check update condition on row i
c
        do 70 i=1,n
          tol=rnf*max(abs(g(i)),abs(gpls(i)))
          if(iagflg.eq.0) tol=tol/sqrt(rnf)
          if(abs(y(i)-t(i)).lt.tol) go to 70
c         if(abs(y(i)-t(i)).ge.tol)
c         then
            skpupd=.false.
            go to 80
c         end if
   70   continue
   80   if(skpupd) go to 110
c
c  bfgs update
c
          do j=1,n
            do i=j,n
              a(i,j)=a(i,j)+y(i)*y(j)/den1-t(i)*t(j)/den2
            end do
          end do
c       end if
c     end if
  110 return
      end
      subroutine sndofd(nr,n,xpls,fcn,fpls,a,sx,rnoise,stepsz,anbr)

c*********************************************************************72
c
cc SNDOFD finds a second order approximation to the Hessian.
c
c  Discussion:
c
c    find second order forward finite difference approximation "a"
c    to the second derivative (hessian) of the function defined by the subp
c    "fcn" evaluated at the new iterate "xpls"
c
c    for optimization use this routine to estimate
c    1) the second derivative (hessian) of the optimization function
c    if no analytical user function has been supplied for either
c    the gradient or the hessian and if the optimization function
c    "fcn" is inexpensive to evaluate.
c
c  Modified:
c
c    09 January 2009
c
c  Author:
c
c    Robert Schnabel, John Koontz, Barry Weiss
c
c  Reference:
c
c    John Dennis, Robert Schnabel,
c    Numerical Methods for Unconstrained Optimization 
c    and Nonlinear Equations,
c    SIAM, 1996,
c    ISBN13: 978-0-898713-64-0,
c    LC: QA402.5.D44.
c
c    Robert Schnabel, John Koontz, Barry Weiss,
c    A modular system of algorithms for unconstrained minimization,
c    Technical Report CU-CS-240-82, 
c    Computer Science Department,
c    University of Colorado at Boulder, 1982.
c
c  Parameters:
c
c nr           --> row dimension of matrix
c n            --> dimension of problem
c xpls(n)      --> new iterate:   x(k)
c fcn          --> name of subroutine to evaluate function
c fpls         --> function value at new iterate, f(xpls)
c a(n,n)      <--  finite difference approximation to hessian
c                  only lower triangular matrix and diagonal
c                  are returned
c sx(n)        --> diagonal scaling matrix for x
c rnoise       --> relative noise in fname (f(x))
c stepsz(n)    --> workspace (stepsize in i-th component direction)
c anbr(n)      --> workspace (neighbor in i-th direction)
c
      implicit double precision ( a-h, o-z )

      dimension xpls(n)
      dimension sx(n)
      dimension stepsz(n),anbr(n)
      dimension a(nr,1)
      external fcn
c
c  find i-th stepsize and evaluate neighbor in direction
c  of i-th unit vector.
c
      ov3 = 1.00D+00 / 3.00D+00
      do i=1,n
        stepsz(i)=rnoise**ov3 * max(abs(xpls(i)), 1.0D+00 /sx(i))
        xtmpi=xpls(i)
        xpls(i)=xtmpi+stepsz(i)
        call fcn(n,xpls,anbr(i))
        xpls(i)=xtmpi
      end do
c
c  calculate column i of a
c
      do i=1,n

        xtmpi=xpls(i)
        xpls(i)=xtmpi+ 2.0D+00 * stepsz(i)
        call fcn(n,xpls,fhat)
        a(i,i)=((fpls-anbr(i))+(fhat-anbr(i)))/(stepsz(i)*stepsz(i))
c
c  calculate sub-diagonal elements of column
c
        if(i.eq.n) go to 30
        xpls(i)=xtmpi+stepsz(i)
        ip1=i+1
        do j=ip1,n
          xtmpj=xpls(j)
          xpls(j)=xtmpj+stepsz(j)
          call fcn(n,xpls,fhat)
          a(j,i)=((fpls-anbr(i))+(fhat-anbr(j)))/(stepsz(i)*stepsz(j))
          xpls(j)=xtmpj
        end do
   30   xpls(i)=xtmpi

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
      subroutine tregup(nr,n,x,f,g,a,fcn,sc,sx,nwtake,stepmx,steptl,
     & dlt,iretcd,xplsp,fplsp,xpls,fpls,mxtake,iprint,lounit,method,
     & udiag)

c*********************************************************************72
c
cc TREGUP accepts the next iterate and updates the trust region.
c
c  Discussion:
c
c    decide whether to accept xpls=x+sc as the next iterate and update the
c    trust region dlt.
c
c  Modified:
c
c    09 January 2009
c
c  Author:
c
c    Robert Schnabel, John Koontz, Barry Weiss
c
c  Reference:
c
c    John Dennis, Robert Schnabel,
c    Numerical Methods for Unconstrained Optimization 
c    and Nonlinear Equations,
c    SIAM, 1996,
c    ISBN13: 978-0-898713-64-0,
c    LC: QA402.5.D44.
c
c    Robert Schnabel, John Koontz, Barry Weiss,
c    A modular system of algorithms for unconstrained minimization,
c    Technical Report CU-CS-240-82, 
c    Computer Science Department,
c    University of Colorado at Boulder, 1982.
c
c  Parameters:
c
c nr           --> row dimension of matrix
c n            --> dimension of problem
c x(n)         --> old iterate x(k-1)
c f            --> function value at old iterate, f(x)
c g(n)         --> gradient at old iterate, g(x), or approximate
c a(n,n)       --> cholesky decomposition of hessian in
c                  lower triangular part and diagonal.
c                  hessian or approx in upper triangular part
c fcn          --> name of subroutine to evaluate function
c sc(n)        --> current step
c sx(n)        --> diagonal scaling matrix for x
c nwtake       --> boolean, =.true. if newton step taken
c stepmx       --> maximum allowable step size
c steptl       --> relative step size at which successive iterates
c                  considered close enough to terminate algorithm
c dlt         <--> trust region radius
c iretcd      <--> return code
c                    =0 xpls accepted as next iterate;
c                       dlt trust region for next iteration.
c                    =1 xpls unsatisfactory but accepted as next iterate
c                       because xpls-x .lt. smallest allowable
c                       step length.
c                    =2 f(xpls) too large.  continue current iteration
c                       with new reduced dlt.
c                    =3 f(xpls) sufficiently small, but quadratic model
c                       predicts f(xpls) sufficiently well to continue
c                       current iteration with new doubled dlt.
c xplsp(n)    <--> workspace (value needs to be retained between
c                  succesive calls of k-th global step)
c fplsp       <--> (retain value between successive calls)
c xpls(n)     <--  new iterate x(k)
c fpls        <--  function value at new iterate, f(xpls)
c mxtake      <--  boolean flag indicating step of maximum length used
c iprint       --> amount of output desired.
c lounit       --> device to which to send output
c method       --> algorithm to use to solve minimization problem
c                    =1 line search
c                    =2 double dogleg
c                    =3 more-hebdon
c udiag(n)     --> diagonal of hessian in a(.,.)
c
      implicit double precision ( a-h, o-z )

      double precision r8vec_dot
      dimension x(n),xpls(n),g(n)
      dimension sx(n),sc(n),xplsp(n)
      dimension a(nr,1)
      logical nwtake,mxtake
      dimension udiag(n)
      external fcn

      mxtake=.false.
      do i=1,n
        xpls(i)=x(i)+sc(i)
      end do

      call fcn(n,xpls,fpls)
      dltf=fpls-f
      slp = r8vec_dot ( n, g, sc )
c
c  next statement added for case of compilers which do not optimize
c  evaluation of next "if" statement (in which case fplsp could be
c  undefined).
c
      if(iretcd.eq. 4) fplsp = 0.0D+00
      if(iprint.gt. 1)write(lounit,1100) iretcd,fpls,fplsp,dltf,slp
      if(iretcd.ne. 3 .or. (fpls.lt.fplsp.and.dltf.le. 1.0D-04*slp))
     &                                                     go to 30
c
c  reset xpls to xplsp and terminate global step
c
        iretcd=0

        do i=1,n
          xpls(i)=xplsp(i)
        end do

        fpls=fplsp
        dlt = 0.5D+00 * dlt
        if(iprint.gt. 1)write(lounit,1000)
        go to 180
c     else
c
c  fpls too large
c
   30   if(dltf.le. 1.0D-04 *slp) go to 80
          if(iprint.gt. 1)write(lounit,1010)
          rln = 0.0D+00
          do i=1,n
            rln=max(rln,abs(sc(i))/max(abs(xpls(i)), 1.0D+00 /sx(i)))
          end do
          if(iprint.gt. 1)write(lounit,1110) rln
          if(rln.ge.steptl) go to 50
c         if(rln.lt.steptl)
c         then
c
c  cannot find satisfactory xpls sufficiently distinct from x
c
            iretcd=1
            if(iprint.gt. 1)write(lounit,1030)
            go to 180
c         else
c
c  reduce trust region and continue global step
c
   50       iretcd=2
            dltmp=-slp*dlt/( 2.0D+00 *(dltf-slp))
            if(iprint.gt. 1)write(lounit,1120) dltmp
            if(dltmp.ge. 0.1D+00 * dlt) go to 60
              dlt= 0.1D+00 * dlt
              go to 70
   60         dlt=dltmp
   70       continue
            if(iprint.gt. 1)write(lounit,1040)
            go to 180
c
c  fpls sufficiently small
c
   80     continue
          if(iprint.gt. 1)write(lounit,1070)
          dltfp = 0.0D+00
          if (method .eq. 3) go to 110
          do i = 1, n
            temp = 0.0D+00
            do j = i, n
              temp = temp + (a(j, i)*sc(j))
            end do
            dltfp = dltfp + temp*temp
          end do
          go to 140
c
c         else
c
  110     do 130 i = 1, n
             dltfp = dltfp + udiag(i)*sc(i)*sc(i)
             if (i .eq. n) go to 130
             temp = 0
             ip1 = i + 1
             do j = ip1, n
                temp = temp + a(i, j)*sc(i)*sc(j)
             end do
             dltfp = dltfp + 2.0D+00 * temp
  130     continue
c
c         end if
c
  140     dltfp = slp + dltfp / 2.0D+00
          if(iprint.gt. 1)write(lounit,1130) dltfp,nwtake
          if(iretcd.eq. 2 .or. (abs(dltfp-dltf).gt. 0.1D+00 *abs(dltf))
     &         .or. nwtake .or. (dlt.gt. 0.99D+00 *stepmx)) go to 160
c         if(iretcd.ne. 2 .and. (abs(dltfp-dltf) .le. 0.1D+00 *abs(dltf))
c    &         .and. (.not.nwtake) .and. (dlt.le. 0.99D+00 *stepmx))
c         then
c
c  double trust region and continue global step
c
            iretcd=3
            do i=1,n
              xplsp(i)=xpls(i)
            end do
            fplsp=fpls
            dlt=min ( 2.0D+00 * dlt,stepmx)
            if(iprint.gt. 1)write(lounit,1080)
            go to 180
c         else
c
c  accept xpls as next iterate.  choose new trust region.
c
  160       continue
            if(iprint.gt. 1)write(lounit,1090)
            iretcd=0
            if(dlt.gt. 0.99D+00 *stepmx) mxtake=.true.
            if(dltf.lt. 0.1D+00 *dltfp) go to 170
c           if(dltf.ge. 0.1D+00 *dltfp)
c           then
c
c  decrease trust region for next iteration
c
              dlt= 0.5D+00 * dlt
              go to 180
c           else
c
c  check whether to increase trust region for next iteration
c
  170         if(dltf.le. 0.75D+00 *dltfp) 
     &          dlt=min ( 2.0D+00 *dlt,stepmx)

c           end if
c         end if
c       end if
c     end if
  180 continue
      if(iprint.gt. 1)write(lounit,1020)
      if(iprint.gt. 1)write(lounit,1050) iretcd,mxtake,dlt,fpls
      if(iprint.gt. 1)write(lounit,1060)
      if(iprint.gt. 1)write(lounit,1140) (xpls(i),i=1,n)
      return
 1000 format('0tregup    reset xpls to xplsp. terminate global step')
 1010 format(26h0tregup    fpls too large.)
 1020 format(38h0tregup    values after call to tregup)
 1030 format(54h tregup    cannot find satisfactory xpls distinct from,
     &       27h x.  terminate global step.)
 1040 format(53h tregup    reduce trust region. continue global step.)
 1050 format(21h tregup       iretcd=,i3/
     &       21h tregup       mxtake=,l1/
     &       21h tregup       dlt   =,g16.8/
     &       21h tregup       fpls  =,g16.8)
 1060 format(32h tregup       new iterate (xpls))
 1070 format(35h tregup    fpls sufficiently small.)
 1080 format(54h tregup    double trust region.  continue global step.)
 1090 format(' tregup    accept xpls.  choose new trust region.',
     &       ' terminate global step.')
 1100 format(18h tregup    iretcd=,i5/
     &       18h tregup    fpls  =,g16.8/
     &       18h tregup    fplsp =,g16.8/
     &       18h tregup    dltf  =,g16.8/
     &       18h tregup    slp   =,g16.8)
 1110 format(18h tregup    rln   =,g16.8)
 1120 format(18h tregup    dltmp =,g16.8)
 1130 format(18h tregup    dltfp =,g16.8/
     &       18h tregup    nwtake=,l1)
 1140 format(14h tregup       ,5(g16.8,3x))
      end
