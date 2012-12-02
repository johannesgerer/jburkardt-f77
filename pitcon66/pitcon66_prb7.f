      program main

c*********************************************************************72
c
cc MAIN is the main program for PITCON66_PRB7.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    22 August 2011
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Ivo Babuska, Werner Rheinboldt,
c    Reliable Error Estimations and Mesh Adaptation for the Finite 
c    Element Method,
c    in International Conference on Computational Methods 
c    in Nonlinear Mechanics,
c    edited by John Oden,
c    Elsevier, 1980,
c    ISBN: 0444853820,
c    LC: QA808.I57.
c
c    John Oden,
c    Finite Elements of Nonlinear Continua,
c    Dover, 2006,
c    ISBN: 0486449734,
c    LC: QA808.2.O33.
c
c
c  The general problem to be solved is the two point boundary value problem:
c
c    -D/DX PHI(X,U,DU/DX,T)  +  PSI(X,U,DU/DX,T) = 0,
c
c    U(X=0)=U0, U(X=1)=U1.
c
c  The particular problem solved here has the form:
c
c    U'' + T*SIN(U+U**2+U**3)=0
c
c    U(0)=U(1)=0.0
c
c  U is approximated by dividing the interval into NINT subintervals,
c  and using the NPOLYS Legendre functions of orders 0 through NPOLYS-1
c  on each subinterval.  The NINT*NPOLYS coefficients of these basis
c  functions form the representation of U.
c
c  Moreover, there are NINT*NPOLYS finite element equations set up,
c  of the form
c
c  INTEGRAL(Subinterval I) -PHI(X,U,U',T)*PL'(J) + PSI(X,U,U',T)*PL(J) = 0
c
c  There is a complication however, in that Lagrange multipliers are
c  used to enforce continuity of the solution U at the interface
c  nodes between neighboring subintervals.  These conditions create
c  new variables, and new equations, and even modify the above standard
c  finite element equations.
c
c  Options:
c
c  ICHOOZ=1  Piecewise linear, 1 continuity condition
c  ICHOOZ=2  Piecewise cubic, 1 continuity condition
c  ICHOOZ=3  Piecewise cubic, 2 continuity conditions
c  ICHOOZ=4  Piecewise quintic, 1 continuity condition
c  ICHOOZ=5  Piecewise quintic, 2 continuity conditions
c  ICHOOZ=6  Piecewise quintic, 3 continuity conditions
c
c  All options use 8 intervals.
c
      implicit none

      integer liw
      integer lrw
      integer mint
      integer mvar

      integer ngauss
      parameter (ngauss=8)

      integer maxpol
      parameter (maxpol=10)

      parameter (mint=8)
      parameter (mvar=mint*6+2+(mint-1)*3)
      parameter (liw=mvar+29)
      parameter (lrw=29+(6+mvar)*mvar)

      external bval
      external denslv
      external fp0011
      external fx0011
      external init
      external pitcon
      external uval

      intrinsic mod

      double precision bcone
      double precision bczero
      double precision dbcodt
      double precision dbczdt
      double precision fpar(1)
      double precision gcoef
      double precision gpoint
      integer i
      integer ichooz
      integer ierror
      integer ihold
      integer ii
      integer ipar(1)
      integer iskip
      integer iwork(liw)
      character name*12
      integer nbcz
      integer nbco
      integer ndrv
      integer nint
      integer npolys
      integer nvar
      integer nvary
      integer nvarz
      double precision pl
      double precision pld
      double precision rwork(lrw)
      double precision theta
      double precision u
      double precision uprym
      double precision xg
      double precision xl
      double precision xr(mvar)
      double precision xrr

      save /intmem/
      save /relmem/

      common /intmem/ npolys,ndrv,nvary,nvarz,nbcz,nbco,nint

      common /relmem/ bczero(8),dbczdt(8),bcone(8),dbcodt(8),
     &  pl(8),pld(8),gcoef(8),gpoint(8),theta(10,10)
c
c  Zero out work arrays
c
      do i=1,liw
        iwork(i)=0
      enddo

      do i=1,lrw
        rwork(i)=0.0
      enddo
c
c  Initialize COMMON arrays
c
      call init(gcoef,gpoint,maxpol,ngauss,theta)
c
c  Set number of variables
c
      ichooz=1
      nint=8
      if(ichooz.eq.1)npolys=2
      if(ichooz.eq.2.or.ichooz.eq.3)npolys=4
      if(ichooz.eq.4.or.ichooz.eq.5.or.ichooz.eq.6)npolys=6
      nvary=nint*npolys
      nbcz=1
      if(ichooz.eq.1.or.ichooz.eq.2.or.ichooz.eq.4)ndrv=1
      if(ichooz.eq.3.or.ichooz.eq.5)ndrv=2
      if(ichooz.eq.6)ndrv=3
      nbco=1
      nvarz=nbcz+(nint-1)*ndrv+nbco
      nvar=nvary+nvarz+1
c
c  Starting point
c
      do i=1,nvar
        xr(i)=0.0
      enddo
c
c  We perturb the zero solution, by setting the 5-th coefficient
c  to some nonzero value.  We will then also require that the program
c  find an exact solution by correcting this starting point, but
c  without altering the 5-th coefficient, thus guaranteeing that the
c  solution is not all zero.
c
      ihold=5
      xr(ihold)=0.001
c
c  IWORK(1)=0     ; This is a startup
c  IWORK(2)=IHOLD ; Use index IHOLD for first parameter
c  IWORK(3)=0     ; Program may choose index
c  IWORK(4)=0     ; Update jacobian every newton step
c  IWORK(5)=0     ; Seek no target values.
c  IWORK(6)=NVAR  ; Seek limit points in index NVAR
c  IWORK(7)=1     ; Moderate amount of output
c  IWORK(9)=0     ; Use user's jacobian routine
c
c  IWORK(17)=20   ; (Unusual setting)  Increase the allowed number of Newton
c                   iterations to 20.
c
      iwork(1)=0
      iwork(2)=ihold
      iwork(3)=0
      iwork(4)=0
      iwork(5)=0
      iwork(6)=nvar
      iwork(7)=2
      iwork(9)=0
      iwork(17)=20
c
c  RWORK(1)=0.0001 ; Absolute error tolerance
c  RWORK(2)=0.0001 ; Relative error tolerance
c  RWORK(3)=0.5    ; Minimum stepsize
c  RWORK(4)=0.5    ; Maximum stepsize
c  RWORK(5)=0.25   ; Starting stepsize
c  RWORK(6)=1.0    ; Starting direction
c  RWORK(7)=0.0    ; Target value
c
      rwork(1)=0.0001
      rwork(2)=0.0001
      rwork(3)=0.001
      rwork(4)=0.5
      rwork(5)=0.25
      rwork(6)=1.0
      rwork(7)=0.0

      write ( *, '(a)' ) ' '
      call timestamp ( )
      write(*,*)' '
      write(*,*)'PCPRB7'
      write ( *, '(a)' ) '  FORTRAN77 version'
      write(*,*)'  PITCON sample program'
      write(*,*)'  The materially nonlinear problem.'
      write(*,*)' '

      if(ichooz.eq.1)then
        write(*,*)'Piecewise linears'
      elseif(ichooz.eq.2.or.ichooz.eq.3)then
        write(*,*)'Piecewise cubics'
      else
        write(*,*)'Piecewise quintics.'
      endif

      if(ichooz.eq.6)then
        write(*,*)'3 continuity conditions at breakpoints'
      elseif(ichooz.eq.5.or.ichooz.eq.3)then
        write(*,*)'2 continuity conditions at breakpoints'
      else
        write(*,*)'1 continuity condition at breakpoints.'
      endif

      write(*,*)' '
      write(*,*)'Step  Type of point     Lambda'
      write(*,*)' '
      i=0
      name='Start point  '
      write(*,'(1x,i3,2x,a12,2x,3g14.6)')i,name,xr(nvar)
c
c  Take another step.
c
      do i=1,30

        call pitcon(fp0011,fpar,fx0011,ierror,ipar,iwork,liw,
     &    nvar,rwork,lrw,xr,denslv)

        if(iwork(1).eq.1)then
          name='corrected    '
        elseif(iwork(1).eq.2)then
          name='continuation '
        elseif(iwork(1).eq.3)then
          name='target point '
        elseif(iwork(1).eq.4)then
          name='limit point  '
        elseif(iwork(1).lt.0)then
          name='jacobian   '
        endif

        write(*,'(1x,i3,2x,a12,2x,3g14.6)')i,name,xr(nvar)
c
c  Print out shape of rod every 5 steps
c
        if(mod(i,5).eq.0)then
          write(*,*)' '
          write(*,*)'Current rod coordinates:'
          write(*,*)' '
          write(*,*)'X,  U(X),  U''(X)'
          write(*,*)' '

          do ii=1,nint

            xl=(ii-1)/real(nint)
            xrr=(ii)/real(nint)
            iskip=(ii-1)*npolys
            xg=0.5*(xl+xrr)
            call bval(npolys,pl,pld,xg,xl,xrr)
            call uval(npolys,pl,pld,xg,u,uprym,xr,iskip)
            write(*,'(1x,3g14.6)')xg,u,uprym
          enddo

          write(*,*)' '

        endif

        if(iwork(1).eq.3)go to 60

        if(ierror.ne.0)then
          write(*,*)'Iteration halted, IERROR=',ierror
          go to 60
        endif

      enddo

60    continue
      write(*,*)' '
      write(*,*)'Jacobians      ',iwork(19)
      write(*,*)'Factorizations ',iwork(20)
      write(*,*)'Solves         ',iwork(21)
      write(*,*)'Functions      ',iwork(22)
      write ( *, '(a)' ) ' '
      call timestamp ( )
      stop
      end
      subroutine init(gcoef,gpoint,maxpol,ngauss,theta)

c*********************************************************************72
c
cc INIT initializes the values of the COMMON block variables.
c
c  Discussion:
c
c    GCOEF, GPOINT, and THETA are the 8 point Gauss integration
c    coefficients and abscissas, and the value and derivatives of the Legendre
c    polynomials at 1.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
      implicit none

      integer maxpol
      integer ngauss

      double precision gcoef(ngauss)
      double precision gpoint(ngauss)
      integer i
      integer j
      double precision theta(maxpol,maxpol)

      gcoef(1)=.1012285363D+00
      gcoef(2)=.2223810345D+00
      gcoef(3)=.3137066459D+00
      gcoef(4)=.3626837834D+00
      gcoef(5)=.3626837834D+00
      gcoef(6)=.3137066459D+00
      gcoef(7)=.2223810345D+00
      gcoef(8)=.1012285363D+00

      gpoint(1)= .9602898565D+00
      gpoint(2)= .7966664774D+00
      gpoint(3)= .5255324099D+00
      gpoint(4)= .1834346425D+00
      gpoint(5)=-.1834346425D+00
      gpoint(6)=-.5255324099D+00
      gpoint(7)=-.7966664774D+00
      gpoint(8)=-.9602898565D+00
c
c  THETA(I,J) is derivative (I-1) of Legendre polynomial (J-1) at X=1.
c
      do j=1,maxpol

        theta(1,j)=1.0
        do i=2,maxpol
          theta(i,j)=(j+i-2)*(j-i+1)*theta(i-1,j)/(2*i-2)
        enddo

      enddo

      return
      end
      subroutine fx0011(nvar,fpar,ipar,x,fx,ierror)

c*********************************************************************72
c
cc FX0011 evaluates the nonlinear function which constrains the shape of the rod.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
      implicit none

      external bounds
      external bval
      external fval
      external uval

      integer nvar

      double precision bcone
      double precision bczero
      double precision dbcodt
      double precision dbczdt
      double precision dtdx
      double precision dtdxl
      double precision dtdxr
      double precision fpar(*)
      double precision fx(nvar)
      double precision gcoef
      double precision gpoint
      double precision h2i
      double precision h2il
      double precision h2ir
      integer i
      integer ieqn
      integer ierror
      integer indx
      integer ipar(*)
      integer iskip
      integer ivar
      integer j
      integer k
      integer khi
      integer l
      integer lhil
      integer lhir
      integer nbcz
      integer nbco
      integer ncl
      integer ncr
      integer ndrv
      integer ndsum
      integer neqn
      integer nint
      integer nodes
      integer npolys
      integer npsum
      integer nvary
      integer nvarz
      double precision phi
      double precision phipt
      double precision phipu
      double precision phipup
      double precision pl
      double precision pld
      double precision psi
      double precision psipt
      double precision psipu
      double precision psipup
      double precision s
      double precision tcon
      double precision term
      double precision terml
      double precision termr
      double precision theta
      double precision u
      double precision uprym
      double precision x(nvar)
      double precision xc
      double precision xg
      double precision xl
      double precision xr

      save /intmem/
      save /relmem/

      common /intmem/ npolys,ndrv,nvary,nvarz,nbcz,nbco,nint

      common /relmem/ bczero(8),dbczdt(8),bcone(8),dbcodt(8),
     &  pl(8),pld(8),gcoef(8),gpoint(8),theta(10,10)
c
c  Zero out FX
c
      tcon=x(nvar)
      call bounds
      neqn=nvar-1

      do j=1,neqn
        fx(j)=0.0
      enddo
c
c  1. Set up the terms (A*Y) involving the bivariate form.
c
c  For each interval I:
c
      do 40 i=1,nint

        xl=(i-1)/real(nint)
        xr=(i)/real(nint)
        iskip=(i-1)*npolys
c
c  For each Gauss point J, evaluate the integrand.
c
        do 30 j=1,8

          xg=xl+(gpoint(j)+1.0)*(xr-xl)/2.0
          call bval(npolys,pl,pld,xg,xl,xr)
          call uval(npolys,pl,pld,xg,u,uprym,x,iskip)
          call fval(xg,u,uprym,tcon,phi,phipu,phipup,phipt,
     *              psi,psipu,psipup,psipt)
c
c  L   PROJECT ONTO EACH TEST FUNCTION PL(L) AND PLD(L)
c
          ieqn=iskip

          do l=1,npolys
            ieqn=ieqn+1
            term=gcoef(j)*(xr-xl)*(psi*pl(l)+phi*pld(l))/2.0
            fx(ieqn)=fx(ieqn)+term
          enddo

   30     continue
   40   continue
c
c  2. ADD THE TERMS INVOLVING THE CONTINUITY OF THE TEST FUNCTIONS
c  WHICH ARE THE TERMS B*Z IN  F=A*Y + B*Z
c
      ncl=nvary
      ncr=nvary+nbcz
c
c  I  FOR EACH INTERVAL
c
      do 110 i=1,nint

        xl=(i-1)/real(nint)
        xr=(i)/real(nint)
        dtdx=2.0/(xr-xl)
c
c  J   FOR THE POLYNOMIALS USED IN APPROXIMATING EACH U,
c  COUNT CONDITIONS AT LEFT ENDPOINT, LHIL, AND AT RIGHT, LHIR
c  IF WE ARE IN THE FIRST OR LAST INTERVAL, ONE OF
c  THESE WILL BE BOUNDARY CONDITIONS
c
        lhil=ndrv
        lhir=ndrv
        if (i.eq.1) lhil=nbcz
        if (i.eq.nint) lhir=nbco
        if (lhil.le.0.and.lhir.le.0) go to 100
c
c  K  FOR EACH TEST FUNCTION PL(K)
c
        do 90 k=1,npolys
          terml=0.0
          termr=0.0
          s=(-1.0)**(k+1)
          ieqn=(i-1)*npolys+k
c
c  L  FOR DERIVATIVES 0 THRU NBCZ SPECIFIED AT 0.0 TO BE ZERO
c  OR DERIVATIVES 0 THRU NDRV SPECIFIED TO BE CONTINUOUS
c  AT INTERIOR NODES
c  OR DERIVATIVES 0 THRU NBCO SPECIFIED AT 1.0 TO BE ZERO
c
          h2i=1.0
          do l=1,lhil
            s=-s
            indx=ncl+l
            terml=terml+s*h2i*theta(l,k)*x(indx)
            h2i=h2i*dtdx
          enddo

          h2i=1.0
          do l=1,lhir
            indx=ncr+l
            termr=termr+h2i*theta(l,k)*x(indx)
            h2i=h2i*dtdx
          enddo

          fx(ieqn)=fx(ieqn)+terml+termr
   90     continue
        if (lhil.gt.0) ncl=ncl+lhil
        if (lhir.gt.0) ncr=ncr+lhir
  100   continue
  110   continue
c
c  3. CREATE THE TERMS FOR THE U FUNCTIONS AND THEIR DERIVATIVES
c  THE MATRIX TERMS  ( C*Y )
c  ONE EQUATION IS GENERATED FOR COMPONENT AND CONDITION
c
      npsum=0
      dtdxr=0.0
      dtdxl=0.0
c
c  I  FOR EACH NODE
c
      ndsum=nvary
      nodes=nint+1
      do 210 i=1,nodes

        if (i.gt.1) xl=(i-2)/real(nint)
        xc=(i-1)/real(nint)
        if (i.lt.nodes) xr=(i)/real(nint)
        if (xc.ne.xl) dtdxl=2.0/(xc-xl)
        if (xr.ne.xc) dtdxr=2.0/(xr-xc)
        h2il=1.0
        h2ir=1.0
c
c  K   FOR DERIVATIVES 0 THRU NBCZ SPECIFIED AT 0.0
c  OR DERIVATIVES 0 THRU NDRV SPECIFIED TO BE CONTINUOUS
c  AT INTERIOR NODES
c  OR DERIVATIVES 0 THRU NBCO SPECIFIED AT 1.0
c
        khi=ndrv
        if (i.eq.1) khi=nbcz
        if (i.eq.nodes) khi=nbco
        if (khi.le.0) go to 190
        do 180 k=1,khi
          s=(-1.0)**(k+1)
c
c  L   SET UP THE TERM FROM THE LEFT HAND INTERVAL
c
          terml=0.0
          if (i.gt.1) go to 120
          terml=bczero(k)
          go to 140
  120     do 130 l=1,npolys
            ivar=npsum+l-npolys
            terml=terml+x(ivar)*h2il*theta(k,l)
  130       continue
c
c  L   SET UP THE TERM FROM THE RIGHT HAND INTERVAL
c
  140     if (i.lt.nodes) go to 150
          termr=-bcone(k)
          go to 170
  150     termr=0.0
          do 160 l=1,npolys
            ivar=npsum+l
            s=-s
            termr=termr+s*x(ivar)*h2ir*theta(k,l)
  160       continue
  170     ieqn=ndsum+k
          fx(ieqn)=terml+termr
          h2il=h2il*dtdxl
          h2ir=h2ir*dtdxr
  180     continue
        ndsum=ndsum+khi
  190   npsum=npsum+npolys
  200   continue
  210   continue
      return
      end
      subroutine fp0011(nvar,fpar,ipar,x,fprime,ierror)

c*********************************************************************72
c
cc FP0011 evaluates the jacobian matrix.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
      implicit none

      external bounds
      external bval
      external fval
      external uval

      integer nvar

      double precision bcone
      double precision bczero
      double precision dbcodt
      double precision dbczdt
      double precision dtdx
      double precision fpar(*)
      double precision fprime(nvar,nvar)
      double precision gcoef
      double precision gpoint
      double precision h2i
      integer i
      integer ieqn
      integer ierror
      integer ipar(*)
      integer iskip
      integer ivar
      integer j
      integer k
      integer khil
      integer khir
      integer l
      integer lhil
      integer lhir
      integer n
      integer nbcz
      integer nbco
      integer ncl
      integer ncr
      integer ndrv
      integer nint
      integer npolys
      integer npsum
      integer nvary
      integer nvarz
      double precision phi
      double precision phipt
      double precision phipu
      double precision phipup
      double precision pl
      double precision pld
      double precision psi
      double precision psipt
      double precision psipu
      double precision psipup
      double precision s
      double precision sum
      double precision tcon
      double precision term
      double precision term1
      double precision term2
      double precision term3
      double precision term4
      double precision theta
      double precision u
      double precision uprym
      double precision x(nvar)
      double precision xg
      double precision xl
      double precision xr

      save /intmem/
      save /relmem/

      common /intmem/ npolys,ndrv,nvary,nvarz,nbcz,nbco,nint

      common /relmem/ bczero(8),dbczdt(8),bcone(8),dbcodt(8),
     &  pl(8),pld(8),gcoef(8),gpoint(8),theta(10,10)
c
c  ZERO OUT THE MATRIX
c
      tcon=x(nvar)
      call bounds
      do 20 i=1,nvar
        do 10 j=1,nvar
          fprime(i,j)=0.0
   10     continue
   20   continue
c
c  1.  SET UP THE TERMS FROM THE BIVARIATE FORM (A*Y )
c
c      I  FOR EACH INTERVAL (XL,XR)
c
      do 80 i=1,nint

        xl=(i-1)/real(nint)
        xr=(i)/real(nint)
        iskip=(i-1)*npolys
c
c  J  FOR EACH GAUSS POINT XG IN THE INTERVAL
c  NOTE THAT THE LEGENDRE POLYNOMIALS AND DERIVATIVES ARE SET UP
c  BY THE CALL TO VALUE
c
        do 70 j=1,8
          xg=xl+0.5*(xr-xl)*(gpoint(j)+1.0)
          call bval(npolys,pl,pld,xg,xl,xr)
          call uval(npolys,pl,pld,xg,u,uprym,x,iskip)
          call fval(xg,u,uprym,tcon,phi,phipu,phipup,phipt,
     *              psi,psipu,psipup,psipt)
c
c  L  FOR EACH LEGENDRE POLYNOMIAL COEFFICIENT
c
          ieqn=iskip
          do 50 l=1,npolys
            ieqn=ieqn+1
            term=gcoef(j)*0.5*(xr-xl)*(psipt*pl(l)+phipt*pld(l))
            fprime(ieqn,nvar)=fprime(ieqn,nvar)+term
c
c  N   FOR EACH Y-COEFFICIENT OF A U
c
            do 30 n=1,npolys
              ivar=npolys*(i-1)+n
              term1=psipu*pl(n)*pl(l)
              term2=psipup*pld(n)*pl(l)
              term3=phipu*pl(n)*pld(l)
              term4=phipup*pld(n)*pld(l)
              sum=term1+term2+term3+term4
              fprime(ieqn,ivar)=fprime(ieqn,ivar)+
     *        gcoef(j)*0.5*(xr-xl)*sum
   30         continue
   50       continue
   70     continue
   80   continue
c
c  2. ADD THE TERMS INVOLVING THE CONTINUITY OF THE TEST FUNCTIONS
c  WHICH ARE THE TERMS B*Z IN  F=A*Y+B*Z
c
c
c     I  FOR EACH INTERVAL,
c
      do 150 i=1,nint
        ncl=nvary
        if (i.gt.1) ncl=nvary+nbcz+(i-2)*ndrv
        ncr=nvary+nbcz+(i-1)*ndrv

        xl=(i-1)/real(nint)
        xr=(i)/real(nint)
        dtdx=2.0/(xr-xl)
c
c  J   FOR THE POLYNOMIALS USED IN APPROXIMATING EACH U,
c      COUNT CONDITIONS AT LEFT ENDPOINT, LHIL, AND AT RIGHT, LHIR
c      IF WE ARE IN THE FIRST OR LAST INTERVAL, ONE OF
c      THESE WILL BE BOUNDARY CONDITIONS
c
          lhil=ndrv
          lhir=ndrv
          if (i.eq.1) lhil=nbcz
          if (i.eq.nint) lhir=nbco
          if (lhil.le.0.and.lhir.le.0) go to 130
c
c  K  FOR EACH TEST FUNCTION PL(K)
c
          do 120 k=1,npolys
            s=(-1.0)**(k+1)
            ieqn=(i-1)*npolys+k
c
c  L  FOR DERIVATIVES 0 THRU NBCZ SPECIFIED AT 0.0 TO BE ZERO
c     OR DERIVATIVES 0 THRU NDRV SPECIFIED TO BE CONTINUOUS
c     AT INTERIOR NODES
c     OR DERIVATIVES 0 THRU NBCO SPECIFIED AT 1.0 TO BE ZERO
c
c     EVALUATE CONTRIBUTION FROM LEFT ENDPOINT
c
            h2i=1.0
            do 90 l=1,lhil
              s=-s
              ivar=ncl+l
              fprime(ieqn,ivar)=s*h2i*theta(l,k)
              h2i=h2i*dtdx
   90         continue
c
c  EVALUATE CONTRIBUTION FROM RIGHT ENDPOINT
c
            h2i=1.0
            do 110 l=1,lhir
              ivar=ncr+l
              fprime(ieqn,ivar)=h2i*theta(l,k)
              h2i=h2i*dtdx
  110         continue
  120       continue
          if (lhil.gt.0) ncl=ncl+lhil
          if (lhir.gt.0) ncr=ncr+lhir
  130     continue
  140     continue
  150   continue
c
c  3. CREATE THE TERMS FOR THE U FUNCTIONS AND THEIR DERIVATIVES
c  THE MATRIX TERMS  ( C*Y )
c  ONE EQUATION IS GENERATED FOR COMPONENT AND CONDITION
c
c
c  I  FOR EACH INTERVAL
c
      do 230 i=1,nint
        ncl=nvary
        if (i.gt.1) ncl=nvary+nbcz+(i-2)*ndrv
        ncr=nvary+nbcz+(i-1)*ndrv
        npsum=(i-1)*npolys

        xl=(i-1)/real(nint)
        xr=(i)/real(nint)
        dtdx=2.0/(xr-xl)
        h2i=1.0
c
c  K   FOR DERIVATIVES 0 THRU NBCZ SPECIFIED AT 0.0
c      OR DERIVATIVES 0 THRU NDRV SPECIFIED TO BE CONTINUOUS
c      AT INTERIOR NODES
c      OR DERIVATIVES 0 THRU NBCO SPECIFIED AT 1.0
c
        khil=ndrv
        if (i.eq.1) khil=nbcz
        do 170 k=1,khil
          ieqn=ncl+k
c
c  L   SET UP THE TERM FROM THE LEFT HAND ENDPOINT
c
          if (i.eq.1) fprime(ieqn,nvar)=dbczdt(k)
          s=(-1.0)**(k+1)
          do 160 l=1,npolys
            ivar=npsum+l
            s=-s
            fprime(ieqn,ivar)=s*h2i*theta(k,l)
  160       continue
          h2i=h2i*dtdx
  170     continue
        ncl=ncl+khil
        h2i=1.0
        khir=ndrv
        if (i.eq.nint) khir=nbco
        do 200 k=1,khir
          ieqn=ncr+k
          fprime(ieqn,nvar)=0.0
          if (i.eq.nint) fprime(ieqn,nvar)=-dbcodt(k)
          do 190 l=1,npolys
            ivar=npsum+l
            fprime(ieqn,ivar)=h2i*theta(k,l)
  190       continue
          h2i=h2i*dtdx
  200     continue
        ncr=ncr+khir
        npsum=npsum+npolys
  230   continue
      return
      end
      subroutine bval(npolys,pl,pld,xg,xl,xr)

c*********************************************************************72
c
cc BVAL evaluates Legendre polynomials and derivatives.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
      implicit none

      double precision a
      integer i
      integer npolys
      double precision pl(8)
      double precision pld(8)
      double precision t
      double precision xg
      double precision xl
      double precision xr

      t=-1.0+2.0*(xg-xl)/(xr-xl)
      pl(1)=1.0
      pld(1)=0.0
      pl(2)=t
      pld(2)=1.0
      a=0.0
      do 10 i=3,8
        a=a+1.0
        pl(i)=((a+a+1.0)*t*pl(i-1)-a*pl(i-2))/(a+1.0)
        pld(i)=((a+a+1.0)*(t*pld(i-1)+pl(i-1))-a*pld(i-2))/(a+1.0)
10      continue

      do 30 i=1,8
        pld(i)=2.0*pld(i)/(xr-xl)
30      continue
      return
      end
      subroutine uval(npolys,pl,pld,xg,u,uprym,x,iskip)

c*********************************************************************72
c
cc UVAL evaluates the solution U and its spatial derivative UPRYM.
c
c  The routine sums the coefficients times the polynomial values.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
      implicit none

      integer npolys

      integer indx
      integer iskip
      integer j
      double precision pl(npolys)
      double precision pld(npolys)
      double precision u
      double precision uprym
      double precision x(*)
      double precision xg

      u=0.0
      uprym=0.0
      do j=1,npolys
        indx=iskip+j
        u=u+x(indx)*pl(j)
        uprym=uprym+x(indx)*pld(j)
      enddo

      return
      end
      subroutine fval(xg,u,uprym,tcon,phi,phipu,phipup,phipt,
     * psi,psipu,psipup,psipt)

c*********************************************************************72
c
cc FVAL is an auxilliary function.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
      implicit none

      intrinsic cos
      intrinsic sin

      double precision phi
      double precision phipu
      double precision phipup
      double precision phipt
      double precision psi
      double precision psipu
      double precision psipup
      double precision psipt
      double precision tcon
      double precision u
      double precision uprym
      double precision xg

      phi=-uprym
      phipu=0.0
      phipup=-1.0
      phipt=0.0

      psi=tcon*sin(u*(1.0+u*(1.0+u)))
      psipu=tcon*(1.0+u*(2.0+3.0*u))*cos(u*(1.0+u*(1.0+u)))
      psipup=0.0
      psipt=sin(u*(1.0+u*(1.0+u)))

      return
      end
      subroutine bounds

c*********************************************************************72
c
cc BOUNDS sets boundary values for the materially nonlinear problem.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
      implicit none

      double precision bcone
      double precision bczero
      double precision dbcodt
      double precision dbczdt
      double precision gcoef
      double precision gpoint
      double precision pl
      double precision pld
      double precision theta

      save /relmem/

      common /relmem/ bczero(8),dbczdt(8),bcone(8),dbcodt(8),
     &  pl(8),pld(8),gcoef(8),gpoint(8),theta(10,10)

      bcone(1)=0.0
      bczero(1)=0.0
      dbcodt(1)=0.0
      dbczdt(1)=0.0

      return
      end
