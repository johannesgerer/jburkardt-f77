      program main

c*********************************************************************72
c
cc MAIN is the main program for TOMS726_PRB.
c
c  Discussion:
c
c    TOMS726_PRB calls sample problems for the TOMS726 library.
c
c  Modified:
c
c    26 January 2007
c
c  Author:
c
c    John Burkardt
c
      implicit none

      call timestamp ( )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TOMS726_PRB'
      write ( *, '(a)' ) '  FORTRAN77 version'
      write ( *, '(a)' ) '  Test the TOMS726 library.'

      call test01 ( )
      call test02 ( )
      call test03 ( )
c
c  This test is causing floating point exceptions.
c
      if ( .false. ) then
        call test04 ( )
      end if
      call test05 ( )
      call test06 ( )
      call test07 ( )
      call test08 ( )
      call test09 ( )

      call test10 ( )
      call test11 ( )
c
c  Terminate.
c
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TOMS726_PRB'
      write ( *, '(a)' ) '  Normal end of execution.'

      write ( *, '(a)' ) ' '
      call timestamp ( )

      stop
      end
      subroutine test01 ( )

c*********************************************************************72
c
cc TEST01
c
c This test generates the first n beta-coefficients in the recurrence
c relation for the orthogonal polynomials relative to the weight
c function
c
c         ((1-om2*x**2)*(1-x**2))**(-1/2)  on (-1,1)
c
c for om2=.1(.2).9,.99,.999, both in single and double precision,
c using modified moments if  modmom=.true.  and ordinary moments
c otherwise. In the former case, n=80, in the latter, n=20. Printed
c are the double-precision values of the coefficients along with the
c relative errors of the single-precision values.
c
c  Modified:
c
c    17 March 2008
c
c  Reference:
c
c    Walter Gautschi,
c    On Generating Orthogonal Polynomials,
c    SIAM Journal on Scientific and Statistical Computing,
c    Volume 3, Number 3, 1982, pages 289-317.
c
      dimension fnu(160),f(80),f0(80),rr(80),a(159),b(159),
     *alpha(80),beta(80),s(80),s0(160),s1(160),s2(160)
      double precision doom2(7),deps,d1mach,dom2,dnu(160),d(80),
     *d0(80),drr(80),da(159),db(159),dalpha(80),dbeta(80),ds(80),
     *ds0(160),ds1(160),ds2(160)
      logical modmom

      data doom2/.1d0,.3d0,.5d0,.7d0,.9d0,.99d0,.999d0/

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST01'

      write(*,1)
    1 format(/)
      modmom=.true.
      eps=r1mach(3)
      deps=d1mach(3)
      if(modmom) then
        n=80
      else
        n=20
      end if
      ndm1=2*n-1
      do 30 iom=1,7
        dom2=doom2(iom)
        om2=sngl(dom2)
c
c Compute the modified resp. ordinary moments using Eqs. (3.7) and (3.9)
c of the companion paper. On machines with limited exponent range, some
c of the high-order modified moments may underflow, without this having
c any deteriorating effect on the accuracy.
c
        call mm_01_r4(n,eps,modmom,om2,fnu,ierr,f,f0,rr)
        call mm_01_r8(n,deps,modmom,dom2,dnu,iderr,d,d0,drr)
        if(ierr.ne.0 .or. iderr.ne.0) then
          write(*,2) ierr,iderr,om2
    2     format(/5x,'ierr in mm_01_r4 = ',i1,
     &    '  iderr in mm_01_r8 = ',i1,
     *      '  for om2 = ',f8.4/)
          goto 30
        end if
c
c Generate the recursion coefficients for the polynomials defining the
c modified resp. ordinary moments.
c
        if(modmom) then
          call recur(ndm1,3,0.,0.,a,b,ierr)
          call drecur(ndm1,3,0.d0,0.d0,da,db,iderr)
        else
          do k=1,ndm1
            a(k)=0.
            b(k)=0.
            da(k)=0.d0
            db(k)=0.d0
          end do
        end if
c
c Compute the desired recursion coefficients by means of the modified
c Chebyshev algorithm.
c
        call cheb_r4(n,a,b,fnu,alpha,beta,s,ierr,s0,s1,s2)
c
c On machines with limited single-precision exponent range, the routine
c CHEB_R4  may generate an underflow exception, which however is harmless
c and can be ignored.
c
        call cheb_r8(n,da,db,dnu,dalpha,dbeta,ds,iderr,ds0,ds1,ds2)
        write(*,3) ierr,iderr
    3   format(/5x,'ierr in cheb_r4 = ',i3,'  iderr in cheb_r8 = ',i3/)
        write(*,4)
    4   format(/5x,'k',14x,'dbeta(k)'/)
        do 20 k=1,n
          km1=k-1
          if(iderr.eq.0 .or. km1.lt.abs(iderr)) then
            if(ierr.eq.0 .or. km1.lt.abs(ierr)) then
              errb=sngl(dabs(dble(beta(k))-dbeta(k))/dbeta(k))
              if(k.eq.1) then
                write(*,5) km1,dbeta(k),errb,om2
    5           format(1x,i5,d36.28,e12.4,'   om2 =',f6.3)
              else
                write(*,6) km1,dbeta(k),errb
    6           format(1x,i5,d36.28,e12.4)
              end if
            else
              write(*,7) km1,dbeta(k)
    7         format(1x,i5,d36.28)
            end if
          end if
   20   continue
        write(*,1)
   30 continue
      return
      end
      subroutine mm_01_r4(n,eps,modmom,om2,fnu,ierr,f,f0,rr)

c*********************************************************************72
c
cc MM_01_R4 generates the modified (Chebyshev) resp. ordinary
c moments of the weight function
c
c          ((1-om2*x**2)*(1-x**2))**(-1/2)  on (-1,1)
c
c using Eqs. (3.7) resp. (3.9) of the companion paper.
c
      dimension fnu(*),f(n),f0(n),rr(n)
      logical modmom
c
c The array  fnu  is assumed to have dimension  2*n.
c
      ierr=0
      nd=2*n
      ndm1=nd-1
      pi=4.*atan(1.)
c
c Compute the Fourier coefficients of ((1-om2*sin(theta)**2))**(-1/2)
c as minimal solution of a three-term recurrence relation as described
c on pp.310-311 of W. Gautschi,On generating orthogonal polynomials'',
c SIAM J. Sci. Statist. Comput. 3, 1982, 289-317.
c
      q=om2/(2.-om2+2.*sqrt(1.-om2))
      q1=(1.+q*q)/q
      do k=1,n
        f(k)=0.
      end do
      nu=nd
   20 nu=nu+10
      do 30 k=1,n
        f0(k)=f(k)
   30 continue
      if(nu.gt.500) then
        ierr=1
        return
      end if
      r=0.
      s=0.
      do 40 k=1,nu
        n1=nu-k+1
        fn1=real(n1)
        r=-(fn1-.5)/(fn1*q1+(fn1+.5)*r)
        s=r*(2.+s)
        if(n1.le.n) rr(n1)=r
   40 continue
      c0=1./(1.+s)
      f(1)=rr(1)*c0
      if(n.gt.1) then
        do 50 k=2,n
          f(k)=rr(k)*f(k-1)
   50   continue
      end if
      do 60 k=1,n
        if(abs(f(k)-f0(k)).gt.eps*abs(f(k))) goto 20
   60 continue
c
c Compute the desired modified resp. ordinary moments in term of
c the above Fourier coefficients.
c
      fnu(1)=pi*c0
      if(n.eq.1) return
      fnu(2)=0.
      if(n.eq.2) return
      if(modmom) then
        c=2.*pi
        do 70 k=3,ndm1,2
          k1=(k-1)/2
          c=-.25*c
          fnu(k)=c*f(k1)
          fnu(k+1)=0.
   70   continue
      else
        c=.5*pi
        fnu(3)=c*(c0-f(1))
        fnu(4)=0.
        c=-c
        do 90 k=5,ndm1,2
          k1=(k-1)/2
          k1m1=k1-1
          c=-.25*c
          c1=1.
          sum=f(k1)
          do 80 i=1,k1m1
            c1=-c1*real(2*k1-i+1)/real(i)
            sum=sum+c1*f(k1-i)
   80     continue
          c1=-c1*real(k1+1)/real(2*k1)
          sum=sum+c1*c0
          fnu(k)=c*sum
          fnu(k+1)=0.
   90   continue
      end if

      return
      end
      subroutine mm_01_r8(n,deps,modmom,dom2,dnu,ierrd,d,d0,drr)

c*********************************************************************72
c
cc MM_01_R8 is a double-precision version of the routine MM_01_R4.
c
      double precision deps,dom2,dnu(*),d(n),d0(n),drr(n),dpi,dq,
     *dq1,dr,ds,dn1,dc0,dc,dc1,dsum
      logical modmom
c
c The array  dnu  is assumed to have dimension  2*n.
c
      ierrd=0
      nd=2*n
      ndm1=nd-1
      dpi=4.d0*datan(1.d0)
      dq=dom2/(2.d0-dom2+2.d0*dsqrt(1.d0-dom2))
      dq1=(1.d0+dq*dq)/dq
      do k=1,n
        d(k)=0.d0
      end do
      nud=nd
   20 nud=nud+10
      do 30 k=1,n
        d0(k)=d(k)
   30 continue
      if(nud.gt.1000) then
        ierrd=1
        return
      end if
      dr=0.d0
      ds=0.d0
      do 40 k=1,nud
        n1=nud-k+1
        dn1=dble(n1)
        dr=-(dn1-.5d0)/(dn1*dq1+(dn1+.5d0)*dr)
        ds=dr*(2.d0+ds)
        if(n1.le.n) drr(n1)=dr
   40 continue
      dc0=1.d0/(1.d0+ds)
      d(1)=drr(1)*dc0
      if(n.gt.1) then
        do 50 k=2,n
          d(k)=drr(k)*d(k-1)
   50   continue
      end if
      do 60 k=1,n
        if(dabs(d(k)-d0(k)).gt.deps*dabs(d(k))) goto 20
   60 continue
      dnu(1)=dpi*dc0
      if(n.eq.1) return
      dnu(2)=0.d0
      if(n.eq.2) return
      if(modmom) then
        dc=2.d0*dpi
        do 70 k=3,ndm1,2
          k1=(k-1)/2
          dc=-.25d0*dc
          dnu(k)=dc*d(k1)
          dnu(k+1)=0.d0
   70   continue
      else
        dc=.5d0*dpi
        dnu(3)=dc*(dc0-d(1))
        dnu(4)=0.d0
        dc=-dc
        do 90 k=5,ndm1,2
          k1=(k-1)/2
          k1m1=k1-1
          dc=-.25d0*dc
          dc1=1.d0
          dsum=d(k1)
          do 80 i=1,k1m1
            dc1=-dc1*dble(2*k1-i+1)/dble(i)
            dsum=dsum+dc1*d(k1-i)
   80     continue
          dc1=-dc1*dble(k1+1)/dble(2*k1)
          dsum=dsum+dc1*dc0
          dnu(k)=dc*dsum
          dnu(k+1)=0.d0
   90   continue
      end if

      return
      end
      subroutine test02 ( )

c*********************************************************************72
c
c This test generates the first n recursion coefficients for the
c orthogonal polynomials relative to the weight function
c
c        (x**sigma)*ln(1/x)  on (0,1],  sigma = -.5, 0, .5,
c
c where n=100 when using modified (Legendre) moments, and n=12 when
c using ordinary moments. It prints the double-precision values of the
c coefficients as well as the relative errors of the respective single-
c precision values and the maximum relative errors.
c
      dimension a(199),b(199),fnu(200),alpha(100),beta(100),s(100),
     *s0(200),s1(200),s2(200)
      double precision dsigma,da(199),db(199),dnu(200),dalpha(100),
     *dbeta(100),ds(100),ds0(200),ds1(200),ds2(200)
      logical modmom,intexp

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST02'

      modmom=.true.
c
c Generate the recursion coefficients for the polynomials defining the
c modified resp. ordinary moments.
c
      if(modmom) then
        n=100
        ndm1=2*n-1
        call recur(ndm1,2,0.,0.,a,b,ierr)
        call drecur(ndm1,2,0.d0,0.d0,da,db,iderr)
      else
        n=12
        ndm1=2*n-1
        do k=1,ndm1
          a(k)=0.
          b(k)=0.
          da(k)=0.d0
          db(k)=0.d0
        end do
      end if
      do 30 is=1,3
        dsigma=-.5d0+.5d0*dble(is-1)
        sigma=sngl(dsigma)
        if(is.eq.2) then
          intexp=.true.
        else
          intexp=.false.
        end if
c
c Compute the modified resp. ordinary moments using Eqs. (3.12) and
c (3.11) of the companion paper. On machines with limited exponent
c range, some of the high-order modified moments may underflow, without
c this having any deteriorating effect on the accuracy.
c
        call mm_02_r4(n,modmom,intexp,sigma,fnu)
        call mm_02_r8(n,modmom,intexp,dsigma,dnu)
c
c Compute the desired recursion coefficients by means of the modified
c Chebyshev algorithm; for the latter, see, e.g., Section 2.4 of
c W. Gautschi, On generating orthogonal polynomials'', SIAM J. Sci.
c Statist. Comput. 3, 1982, 289-317.
c
        call cheb_r4(n,a,b,fnu,alpha,beta,s,ierr,s0,s1,s2)
c
c On machines with limited single-precision exponent range, the routine
c CHEB_R4  may generate an underflow exception, which however is harmless
c and can be ignored.
c
        call cheb_r8(n,da,db,dnu,dalpha,dbeta,ds,iderr,ds0,ds1,ds2)
        write(*,1) ierr,iderr
    1   format(/6x,'ierr in cheb_r4 = ',i4,' iderr in cheb_r8 = ',i4/)
c
c Compute and print the relative errors and their maxima.
c
        eamax=0.
        ebmax=0.
        write(*,2) sigma
    2   format(/32x,'sigma =',f5.1)
        write(*,3)
    3   format(/3x,'k',13x,'dalpha(k)',25x,'dbeta(k)'/)
        do 20 k=1,n
          km1=k-1
          if(iderr.eq.0 .or. km1.lt.abs(iderr)) then
            write(*,4) km1,dalpha(k),dbeta(k)
    4       format(1x,i3,2d33.25)
            if(ierr.eq.0 .or. km1.lt.abs(ierr)) then
              erra=sngl(dabs(dble(alpha(k))-dalpha(k))/dalpha(k))
              errb=sngl(dabs(dble(beta(k))-dbeta(k))/dbeta(k))
              write(*,5) erra,errb
    5         format(4x,e12.4,21x,e12.4)
              if(erra.gt.eamax) then
                eamax=erra
                kamax=km1
              end if
              if(errb.gt.ebmax) then
                ebmax=errb
                kbmax=km1
              end if
            end if
          end if
   20   continue
        write(*,6) eamax,kamax,ebmax,kbmax
    6   format(/6x,'eamax =',e11.4,' at',i3,9x,'ebmax =',e11.4,
     *    ' at',i3//)
   30 continue

      return
      end
      subroutine mm_02_r4(n,modmom,intexp,sigma,fnu)

c*********************************************************************72
c
c This generates the first  2*n  modified moments (if modmom=.true.)
c relative to shifted monic Legendre polynomials, using Eq. (3.12) of
c the companion paper, and the first  2*n  ordinary moments (if modmom
c =.false.) by Eq. (3.11), of the weight function
c
c          (x**sigma)*ln(1/x)  on (0,1],   sigma > -1,
c
c for sigma an integer (if intexp=.true.) or a real number (if intexp
c =.false.). In either case, the input variable  sigma  is of type real.
c
      dimension fnu(*)
      logical modmom,intexp
c
c The array  fnu  is assumed to have dimension  2*n.
c
      nd=2*n
      sigp1=sigma+1.
      if(modmom) then
        isigma=int(sigma)
        isigp1=isigma+1
        isigp2=isigma+2
        isigp3=isigma+3
        if(intexp .and. isigp1.lt.nd) then
          kmax=isigp1
        else
          kmax=nd
        end if
        c=1.
        do 20 k=1,kmax
          km1=k-1
          fk=real(k)
          p=1.
          s=1./sigp1
          if(kmax.gt.1) then
            do i=1,km1
              fi=real(i)
              p=(sigp1-fi)*p/(sigp1+fi)
              s=s+1./(sigp1+fi)-1./(sigp1-fi)
            end do
          end if
          fnu(k)=c*s*p/sigp1
          c=fk*c/(4.*fk-2.)
   20   continue
        if(.not.intexp .or. isigp1.ge.nd) return
        q=-.5
        if(isigma.gt.0) then
          do 30 iq=1,isigma
            fiq=real(iq)
            q=fiq*fiq*q/((2.*fiq+1.)*(2.*fiq+2.))
   30     continue
        end if
        fnu(isigp2)=c*q
        if(isigp2.eq.nd) return
        do 40 k=isigp3,nd
          km1=k-1
          fkm1=real(km1)
          fnu(k)=-fkm1*(fkm1-sigp1)*fnu(km1)/((4.*fkm1-2.)*
     *      (fkm1+sigp1))
   40   continue
        return
      else
        do 50 k=1,nd
          fkm1=real(k-1)
          fnu(k)=(1./(sigp1+fkm1))**2
   50   continue
      end if
      end
      subroutine mm_02_r8(n,modmom,intexp,dsigma,dnu)

c*********************************************************************72
c
cc MM_02_R8 is a double-precision version of the routine MM_02_R4.
c
      double precision dsigma,dnu(*),dsigp1,dc,dk,dp,ds,di,dq,diq,dkm1
      logical modmom,intexp
c
c The array  dnu  is assumed to have dimension  2*n.
c
      nd=2*n
      dsigp1=dsigma+1.d0
      if(modmom) then
        isigma=idint(dsigma)
        isigp1=isigma+1
        isigp2=isigma+2
        isigp3=isigma+3
        if(intexp .and. isigp1.lt.nd) then
          kmax=isigp1
        else
          kmax=nd
        end if
        dc=1.d0
        do 20 k=1,kmax
          km1=k-1
          dk=dble(k)
          dp=1.d0
          ds=1.d0/dsigp1
          if(kmax.gt.1) then
            do i=1,km1
              di=dble(i)
              dp=(dsigp1-di)*dp/(dsigp1+di)
              ds=ds+1.d0/(dsigp1+di)-1.d0/(dsigp1-di)
            end do
          end if
          dnu(k)=dc*ds*dp/dsigp1
          dc=dk*dc/(4.d0*dk-2.d0)
   20   continue
        if(.not.intexp .or. isigp1.ge.nd) return
        dq=-.5d0
        if(isigma.gt.0) then
          do 30 iq=1,isigma
            diq=dble(iq)
            dq=diq*diq*dq/((2.d0*diq+1.d0)*(2.d0*diq+2.d0))
   30     continue
        end if
        dnu(isigp2)=dc*dq
        if(isigp2.eq.nd) return
        do 40 k=isigp3,nd
          km1=k-1
          dkm1=dble(km1)
          dnu(k)=-dkm1*(dkm1-dsigp1)*dnu(km1)/((4.d0*dkm1-2.d0)*
     *      (dkm1+dsigp1))
   40   continue
        return
      else
        do 50 k=1,nd
          dkm1=dble(k-1)
          dnu(k)=(1.d0/(dsigp1+dkm1))**2
   50   continue
      end if
      return
      end
      subroutine test03 ( )

c*********************************************************************72
c
      dimension x(320),w(320),bexact(320),p0(320),p1(320),p2(320),
     *als(320),bes(320),all(320),bel(320)
c
c This test applies both the Stieltjes procedure (cf. Section 2.1 of
c W. Gautschi, On generating orthogonal polynomials'', SIAM J. Sci.
c Statist. Comput. 3, 1982, 289-317) and the Lanczos algorithm (cf.
c W.B. Gragg and W.J. Harrod, The numerically stable reconstruction of
c Jacobi matrices from spectral data'', Numer. Math. 44, 1984, 317-335)
c to generate the  N  recursion coefficients, N = 40, 80, 160, 320,
c for the (monic) orthogonal polynomials relative to the discrete inner
c product supported on  N  equally spaced points on [-1,1] (including
c the end points) and having equal weights 2/N. The routine prints the
c absolute errors in the alpha-coefficients and the relative errors in
c the beta-coefficients, these being computed using known formulae for
c the coefficients. The maxima of these errors are also printed.
c
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST03'

      ncap=20
      do 30 incap=1,4
        ncap=2*ncap
        fncap=real(ncap)
        fncm1=real(ncap-1)
c
c Generate the abscissae and weights of the discrete inner product and
c the exact beta-coefficients.
c
        x(1)=-1.
        w(1)=2./fncap
        bexact(1)=2.
        do k=2,ncap
          fkm1=real(k-1)
          x(k)=-1.+2.*fkm1/fncm1
          w(k)=w(1)
          bexact(k)=(1.+1./fncm1)**2*(1.-(fkm1/fncap)**2)/
     *      (4.-(1./fkm1)**2)
        end do
c
c Compute the desired coefficients, first by the Stieltjes procedure,
c and then by the Lanczos algorithm. Indicate via the error flag  ierrs
c whether a critical underflow condition has arisen in Stieltjes's
c procedure. (There may, in addition, occur harmless underflow, which
c the routine  sti  does not test for.)
c
        call sti(ncap,ncap,x,w,als,bes,ierrs,p0,p1,p2)
        call lancz(ncap,ncap,x,w,all,bel,ierrl,p0,p1)
        write(*,1) ierrs,ierrl
    1   format(/5x,'ierr in sti = ',i4,11x,'ierr in lancz = ',i3/)
c
c Compute and print the absolute errors of the alpha-coefficients and
c the relative errors of the beta-coefficients as well as the maximum
c respective errors.
c
        erralm=0.
        errblm=0.
        write(*,2)
    2   format(5x,'k',4x,'erra',8x,'errb',10x,'erra',8x,'errb'/)
        do 20 in=1,ncap
          inm1=in-1
          erras=abs(als(in))
          errbs=abs((bes(in)-bexact(in))/bexact(in))
          erral=abs(all(in))
          errbl=abs((bel(in)-bexact(in))/bexact(in))
          if(erral.gt.erralm) erralm=erral
          if(errbl.gt.errblm) errblm=errbl
          if(ierrs.eq.0 .or. inm1.lt.abs(ierrs)) then
            if(in.eq.1) then
              write(*,3) inm1,erras,errbs,erral,errbl,ncap
    3         format(1x,i5,2e12.4,2x,2e12.4,'   N =',i4)
            else
              write(*,4) inm1,erras,errbs,erral,errbl
    4         format(1x,i5,2e12.4,2x,2e12.4)
            end if
          else
            if(in.eq.1) then
              write(*,5) inm1,erral,errbl,ncap
    5         format(1x,i5,26x,2e12.4,'   N =',i4)
            else
              write(*,6) inm1,erral,errbl
    6         format(1x,i5,26x,2e12.4)
            end if
          end if
   20   continue
        write(*,7) erralm,errblm
    7   format(/32x,2e12.4//)
   30 continue
      return
      end
      subroutine test04 ( )

c*********************************************************************72
c
cc TEST04 is a test of the routine  mcdis, which is applied to generate
c the first 80 recurrence coefficients for the orthogonal polynomials
c belonging to the weight function
c
c       (1-x**2)**(-1/2) + c   on (-1,1),  c = 1, 10, 100.
c
c The corresponding inner product is discretized by applying the
c Gauss-Chebyshev quadrature rule to the first term of the weight
c function, and the Gauss-Legendre rule to the second term. In
c addition to the beta-coefficients (all alpha's are zero), the
c routine prints the variables  ncap  and  kount  to confirm
c convergence after one iteration.
c
      external qchle_r4

      dimension a(81),b(81),e(81),xp(1),yp(1),endl(1),endr(1),
     *xfer(1),wfer(1),alpha(80),beta(80),be(80),x(81),w(81),
     *xm(162),wm(162)
      real p0(162)
      real p1(162)
      real p2(162)
      real betap(80,3)
      logical finl,finr

      common /common04/ c,a,b,e,epsma

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST04'
c
c epsma is the machine single precision.
c
      epsma=r1mach(3)
      iq=1
      idelta=2
      irout=1
      n=80
      mc=2
      mp=0
      ncapm=81
      eps=5000.*epsma
      iem=0
      c=.1

      do ic=1,3

        c=10.*c
c
c Compute the desired recursion coefficients. On machines with limited
c exponent range, harmless underflow may occur in the routine  gauss
c used in  QCHLE_R4  to generate the Gauss-Legendre quadrature rule.
c
        call mcdis(n,ncapm,mc,mp,xp,yp,qchle_r4,eps,iq,idelta,irout,finl,
     *    finr,endl,endr,xfer,wfer,alpha,beta,ncap,kount,ierr,ie,be,
     *    x,w,xm,wm,p0,p1,p2)

        write(*,2) ncap,kount,ierr,ie,c
    2   format(3x,'ncap = ',i3,' kount = ',i2,' ierr = ',i3,
     *    ' ie = ',i3,' for c = ',f5.1)

        iem = max ( iem, abs ( ie ) )

        if(ie.ne.0 .and. abs(ie).le.n) then
          call mcdis(abs(ie)-1,ncapm,mc,mp,xp,yp,qchle_r4,eps,iq,
     &    idelta, irout,finl,finr,endl,endr,xfer,wfer,alpha,beta,
     &    ncap,kount,ierr,ie,be,x,w,xm,wm,p0,p1,p2)
          write(*,2) ncap,kount,ierr,ie,c
          write ( * , '(a)' ) ' '
        end if
c
c Assemble the results in an array.
c
        do k=1,n
          if(ie.eq.0 .or. k-1.lt.abs(ie)) then
            betap(k,ic)=beta(k)
          else
            betap(k,ic)=0.
          end if
        end do

      end do
c
c Print the results.
c
      write(*,3)
    3 format(//3x,'k',4x,'beta(k), c=1',6x,'beta(k), c=10',4x,
     *  'beta(k), c=100'/)
      do 30 k=1,n
        km1=k-1
        if(iem.eq.0 .or. km1.lt.iem) then
          write(*,4) km1,betap(k,1),betap(k,2),betap(k,3)
        end if
    4   format(1x,i3,3e18.10)
   30 continue

      return
      end
      subroutine qchle_r4(n,x,w,i,ierr)

c*********************************************************************72
c
cc QCHLE_R4 returns Gauss-Chebyshev or Gauss-Legendre rules.
c
      dimension x(n),w(n),a(81),b(81),e(81)

      common /common04/ c,a,b,e,epsma

      fn=real(n)
      pi=4.*atan(1.)

      if(i.eq.1) then

        do 10 k=1,n
          fk=real(k)
          x(k)=cos((2.*fk-1.)*pi/(2.*fn))
          w(k)=pi/fn
   10   continue

      else

        call recur(n,1,0.,0.,a,b,ierr)
        call gauss(n,a,b,epsma,x,w,ierr,e)
        do 20 k=1,n
          w(k)=c*w(k)
   20   continue

      end if

      return
      end
      subroutine test05 ( )

c*********************************************************************72
c
      external qjac_r4
      dimension a(41),b(41),e(41),xp(1),yp(1),endl(1),endr(1),xfer(1),
     *wfer(1),alpha(40), beta(40),be(40),x(41),w(41),xm(42),wm(42),
     *p0(42),p1(42),p2(42)
      double precision dy,dalj,dbej,dalpbe,da(41),db(41),daex(40),
     *dbex(40),dnum,dd,dkm1,dden
      logical finl,finr

      common /common05/ a,b,e,epsma
c
c This test generates the first 40 recursion coefficients for
c polynomials orthogonal with respect to the Jacobi weight function
c with parameters  alj = -.8(.2)1., bej = -.8(.2)1.  and an added mass
c point of strength  y = .5, 1, 2, 4, 8  at the left end point. It
c also computes the maximum relative errors (absolute errors for alpha-
c coefficients near zero) of the computed coefficients by comparing
c them against the exact coefficients known analytically.
c
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST05'

      write(*,1)
    1 format(/)
c
c epsma is the machine single precision
c
      epsma=r1mach(3)
      iq=1
      idelta=2
      irout=1
      n=40
      mc=1
      mp=1
      xp(1)=-1.
      ncapm=41
      eps=5000.*epsma
      y=.25
      do 40 iy=1,5
        y=2.*y
        dy=dble(y)
        yp(1)=y
        write(*,2) y
    2   format(/1x,'y = ',f6.2/)
        write(*,3)
    3   format(2x,'alj',2x,'bej',5x,'erra',7x,'errb',6x,'alpha',
     *    4x,'beta',4x,'ka',2x,'kb',1x,'ierr',1x,'ie',1x,'it'/)
        do 30 ia=1,10
          alj=-1.+.2*real(ia)
          dalj=dble(alj)
          do 20 ib=1,10
            bej=-1.+.2*real(ib)
            alpbe=alj+bej
            dbej=dble(bej)
            dalpbe=dalj+dbej
c
c Generate the Jacobi recurrence coefficients.
c
            call recur(ncapm,6,alj,bej,a,b,ierr)
            call drecur(ncapm,6,dalj,dbej,da,db,iderr)
c
c Compute the desired recursion coefficients.
c
            call mcdis(n,ncapm,mc,mp,xp,yp,qjac_r4,eps,iq,idelta,irout,
     *        finl,finr,endl,endr,xfer,wfer,alpha,beta,ncap,kount,
     *        ierr,ie,be,x,w,xm,wm,p0,p1,p2)
c
c Compute the exact coefficients by Eqs. (4.19)-(4.21) of the companion
c paper along with the relative errors (absolute errors for alpha-
c coefficients close to zero).
c
            daex(1)=(da(1)-dy)/(1.d0+dy)
            dbex(1)=1.d0+dy
            erra=abs(sngl(dble(alpha(1))-daex(1)))
            if(abs(sngl(daex(1))).gt.eps) erra=erra/abs(sngl(daex(1)))
            errb=abs(sngl((dble(beta(1))-dbex(1))/dbex(1)))
            erram=erra
            alpham=alpha(1)
            errbm=errb
            betam=beta(1)
            kam=0
            kbm=0
            dnum=1.d0+dy
            dd=1.d0
            do 10 k=2,n
              km1=k-1
              dkm1=dble(km1)
              dden=dnum
              if(k.gt.2) dd=(dbej+dkm1)*(dalpbe+dkm1)*dd/((dalj+dkm1
     *          -1.d0)*(dkm1-1.d0))
              dnum=(1.d0+(dbej+dkm1+1.d0)*(dalpbe+dkm1+1.d0)*dy*dd/
     *          (dkm1*(dalj+dkm1)))/(1.d0+dy*dd)
              daex(k)=da(k)+2.d0*dkm1*(dkm1+dalj)*(dnum-1.d0)/
     *          ((dalpbe+2.d0*dkm1)*(dalpbe+2.d0*dkm1+1.d0))
     *          +2.d0*(dbej+dkm1+1.d0)*(dalpbe+dkm1+1.d0)*((1.d0/dnum)
     *          -1.d0)/((dalpbe+2.d0*dkm1+1.d0)*(dalpbe+2.d0*dkm1+2.d0))
              dbex(k)=dnum*db(k)/dden
              erra=abs(sngl(dble(alpha(k))-daex(k)))
              if(abs(sngl(daex(k))).gt.eps) erra=erra/abs(sngl(daex(k)))
              errb=abs(sngl((dble(beta(k))-dbex(k))/dbex(k)))
              if(erra.gt.erram) then
                erram=erra
                alpham=alpha(k)
                kam=km1
              end if
              if(errb.gt.errbm) then
                errbm=errb
                betam=beta(k)
                kbm=km1
              end if
   10       continue
c
c Print the results.
c
            write(*,4) alj,bej,erram,errbm,alpham,betam,kam,kbm,ierr,
     *        ie,kount
    4       format(1x,2f5.2,2e11.4,2f9.6,4i4,i2)
   20     continue
          write(*,1)
   30   continue
   40 continue
      return
      end
      subroutine qjac_r4(n,x,w,i,ierr)

c*********************************************************************72
c
cc QJAC_R4 returns a Gauss-Jacobi rule.
c
      dimension x(n),w(n),a(41),b(41),e(41)
      common /common05/ a,b,e,epsma

      call gauss(n,a,b,epsma,x,w,ierr,e)
      do 10 k=1,n
        w(k)=w(k)/b(1)
   10 continue
      return
      end
      subroutine test06 ( )

c*********************************************************************72
c
      external qlag,dqlag
      dimension a(100),b(100),e(100),xp(1),yp(1),endl(1),endr(1),
     *xfer(1),wfer(1),alpha(40),beta(40),be(40),x(100),w(100),
     *xm(200),wm(200),p0(200),p1(200),p2(200)
      double precision d1mach,da(300),db(300),de(300),depsma,dxp(1),
     *dyp(1),dendl(1),dendr(1),dxfer(1),dwfer(1),dalpha(40),dbeta(40),
     *dbe(40),dx(300),dw(300),dxm(600),dwm(600),dp0(600),dp1(600),
     *dp2(600),deps
      logical finl,finr,finld,finrd
      common /common06/ a,b,e,epsma
      common/common06_2/da,db,de,depsma
c
c This test generates in single and double precision the first 40
c recursion coefficients of the orthogonal polynomials belonging to
c the logistic density function
c
c              exp(-x)/((1+exp(-x))**2)  on (-oo, oo).
c
c It prints the double-precision beta-coefficients (the alpha's being
c all zero) along with the absolute and relative errors of the alpha-
c coefficients resp. beta-coefficients.
c
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST06'

      write(*,1)
    1 format(/)
c
c epsma and depsma are the machine single and double precision.
c
      iq=1
      idelta=1
      irout=1
      epsma=r1mach(3)
      depsma=d1mach(3)
      n=40
      mc=2
      mp=0
      ncapm=100
      ncapmm=mc*ncapm+mp
      ncpmd=300
      ncpmmd=mc*ncpmd+mp
      eps=5000.*epsma
      deps=1000.*depsma
c
c Compute the desired coefficients. On machines with limited exponent
c range, some of the weights in the Gauss-Laguerre quadrature rule may
c underflow.
c
      call mcdis(n,ncapm,mc,mp,xp,yp,qlag,eps,iq,idelta,irout,finl,
     *  finr,endl,endr,xfer,wfer,alpha,beta,ncap,kount,ierr,ie,be,x,w,
     *  xm,wm,p0,p1,p2)
      call dmcdis(n,ncpmd,mc,mp,dxp,dyp,dqlag,deps,iq,idelta,irout,
     *  finld,finrd,dendl,dendr,dxfer,dwfer,dalpha,dbeta,ncapd,kountd,
     *  ierrd,ied,dbe,dx,dw,dxm,dwm,dp0,dp1,dp2)
      write(*,2) ncap,kount,ierr,ie
    2 format(/1x,'ncap = ',i3,' kount = ',i2,' ierr = ',i3,' ie = ',i5)
      write(*,3) ncapd,kountd,ierrd,ied
    3 format(1x,'ncapd= ',i3,' kountd= ',i2,' ierrd= ',i3,' ied= ',i5/)
c
c Print the results.
c
      write(*,4)
    4 format(/5x,'k',13x,'dbeta(k)',17x,'erra',8x,'errb'/)
      do 10 k=1,n
        km1=k-1
        erra=abs(alpha(k))
        errb=abs(sngl((dble(beta(k))-dbeta(k))/dbeta(k)))
        if(ied.eq.0 .or. km1.lt.ied) then
          if(ie.eq.0 .or. km1.lt.ie) then
            write(*,5) km1,dbeta(k),erra,errb
    5       format(1x,i5,d33.25,2e12.4)
          else
            write(*,6) km1,dbeta(k)
    6       format(1x,i5,d33.25)
          end if
        end if
   10 continue
      return
      end
      subroutine qlag(n,x,w,i,ierr)

c*********************************************************************72
c
      dimension x(n),w(n),a(100),b(100),e(100)
      common/common06/a,b,e,epsma
      call recur(n,7,0.,0.,a,b,ierr)
      call gauss(n,a,b,epsma,x,w,ierr,e)
      do 10 k=1,n
        w(k)=w(k)/((1.+exp(-x(k)))**2)
        if(i.eq.1) x(k)=-x(k)
   10 continue
      return
      end
      subroutine dqlag(n,dx,dw,i,ierr)

c*********************************************************************72
c
      double precision dx(n),dw(n),da(300),db(300),de(300),depsma
      common/common06_2/da,db,de,depsma
      call drecur(n,7,0.d0,0.d0,da,db,ierr)
      call dgauss(n,da,db,depsma,dx,dw,ierr,de)
      do 10 k=1,n
        dw(k)=dw(k)/((1.d0+dexp(-dx(k)))**2)
        if(i.eq.1) dx(k)=-dx(k)
   10 continue
      return
      end
      subroutine test07 ( )

c*********************************************************************72
c
      external quad,dquad
      dimension endl(4),endr(4),xp(1),yp(1),xfer(100),wfer(100),
     *alpha(40),beta(40),be(40),x(100),w(100),xm(400),wm(400),p0(400),
     *p1(400),p2(400)
      double precision d1mach,depsma,dendl(4),dendr(4),dxp(1),dyp(1),
     *dxfer(250),dwfer(250),dalpha(40),dbeta(40),dbe(40),dx(250),
     *dw(250),dxm(1000),dwm(1000),dp0(1000),dp1(1000),dp2(1000),di,deps
      logical finl,finr,finld,finrd
c
c This test generates in single and double precision the first 40
c recurrence coefficients of the orthogonal polynomials for the half-
c range Hermite weight function
c
c                   exp(-x**2)  on  (0,oo).
c
c Printed are the double-precision values of the alpha- and beta-
c coefficients along with the respective relative errors of the single-
c precision values.
c
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST07'

      finl=.true.
      finr=.false.
      finld=.true.
      finrd=.false.
      epsma=r1mach(3)
      depsma=d1mach(3)
c
c epsma and depsma are the machine single and double precision.
c
      iq=2
      idelta=1
      irout=1
      n=40
      mc=4
      mcd=4
      mp=0
      ncapm=100
      ncpmd=250
c
c Set up the partition for the discretization of the inner product.
c
      do 10 i=1,4
        fi=real(i)
        di=dble(fi)
        endl(i)=3.*(fi-1.)
        endr(i)=3.*fi
        dendl(i)=3.d0*(di-1.d0)
        dendr(i)=3.d0*di
   10 continue
      eps=50.*epsma
      deps=1000.*depsma
c
c Compute the desired recursion coefficients by the multiple-component
c discretization procedure. On the third and fourth subinterval of the
c partition, the quadrature weights produced by  qgp  resp.  dqgp  may
c underflow, on the fourth subinterval even on machines with large
c exponent range.
c
      call mcdis(n,ncapm,mc,mp,xp,yp,quad,eps,iq,idelta,irout,finl,
     *  finr,endl,endr,xfer,wfer,alpha,beta,ncap,kount,ierr,ie,be,x,w,
     *  xm,wm,p0,p1,p2)
      call dmcdis(n,ncpmd,mcd,mp,dxp,dyp,dquad,deps,iq,idelta,irout,
     *  finld,finrd,dendl,dendr,dxfer,dwfer,dalpha,dbeta,ncapd,kountd,
     *  ierrd,ied,dbe,dx,dw,dxm,dwm,dp0,dp1,dp2)
      write(*,1) ncap,kount,ierr,ie
    1 format(/1x,'ncap = ',i3,' kount = ',i2,' ierr = ',i3,' ie = ',i3)
      write(*,2) ncapd,kountd,ierrd,ied
    2 format(1x,'ncapd =',i3,' kountd =',i2,' ierrd =',i3,' ied =',i3/)
c
c Print the results.
c
      write(*,3)
    3 format(/5x,'k',13x,'dalpha(k)',25x,'dbeta(k)')
      write(*,4)
    4 format(11x,'erra',29x,'errb')
      do 20 k=1,n
        km1=k-1
        erra=abs(sngl((dble(alpha(k))-dalpha(k))/dalpha(k)))
        errb=abs(sngl((dble(beta(k))-dbeta(k))/dbeta(k)))
        write(*,5) km1,dalpha(k),dbeta(k)
    5   format(1x,i5,2d33.25)
        write(*,6) erra,errb
    6   format(6x,e12.4,21x,e12.4)
   20 continue
      return
      end
      subroutine quad(n,x,w,i,ierr)

c*********************************************************************72
c
      dimension x(n),w(n)
      print *,' User has selected the wrong SP-quadrature routine.'
      stop
      end
      subroutine dquad(n,dx,dw,i,ierr)

c*********************************************************************72
c
      double precision dx(n),dw(n)
      print *,' User has selected the wrong DP-quadrature routine.'
      stop
      end
      function wf(x,i)

c*********************************************************************72
c
      wf=exp(-x*x)
      return
      end
      double precision function dwf(dx,i)

c*********************************************************************72
c
      double precision dx
      dwf=dexp(-dx*dx)
      return
      end
      subroutine test08 ( )

c*********************************************************************72
c
      external qcheb
      dimension oom2(7),xp(1),yp(1),endl(1),endr(1),xfer(1),wfer(1),
     *a(79),b(79),fnu(80),alpha(40),beta(40),be(40),x(500),w(500),
     *xm(500),wm(500),s(40),s0(80),s1(80),s2(80)
      logical finl,finr
      common /common08/ om2
      data oom2/.1,.3,.5,.7,.9,.99,.999/
c
c This test reproduces the results of  test1  in single precision,
c for n=40, using the routine  mccheb  in place of CHEB_R4  and
c Gauss-Chebyshev quadrature to discretize the modified moments.
c
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST08'

      write(*,1)
    1 format(/)
      epsma=r1mach(3)
c
c epsma is the machine single precision.
c
      iq=1
      idelta=1
      n=40
      ndm1=2*n-1
      mc=1
      mp=0
      ncapm=500
      eps=100.*epsma
c
c Generate the recurrence coefficients for the (Chebyshev) polynomials
c defining the modified moments.
c
      call recur(ndm1,3,0.,0.,a,b,ierr)
c
c Compute the desired recursion coefficients by the discretized
c Chebyshev algorithm.
c
      do 20 iom=1,7
        om2=oom2(iom)
        call mccheb(n,ncapm,mc,mp,xp,yp,qcheb,eps,iq,idelta,finl,finr,
     *    endl,endr,xfer,wfer,a,b,fnu,alpha,beta,ncap,kount,ierr,be,x,
     *    w,xm,wm,s,s0,s1,s2)
c
c On machines with limited single-precision exponent range, the routine
c CHEB_R4 may have generated an underflow exception, which however is
c harmless and can be ignored.
c
        write(*,2) ncap,kount,ierr
    2   format(/'ncap=',i3,'  kount=',i3,'  ierr=',i3/)
c
c Print the results.
c
        write(*,3)
    3   format(5x,'k',6x,'beta(k)'/)
        do 10 k=1,n
          km1=k-1
          if(k.eq.1) then
            write(*,4) km1,beta(k),om2
    4       format(1x,i5,e18.10,'   om2 =',f6.3)
          else
            write(*,5) km1,beta(k)
    5     format(1x,i5,e18.10)
          end if
   10   continue
        write(*,1)
   20 continue
      return
      end
      subroutine qcheb(n,x,w,i,ierr)

c*********************************************************************72
c
      dimension x(n),w(n)
      common/common08/om2
      fn=real(n)
      pi=4.*atan(1.)
      do 10 k=1,n
        fk=real(k)
        x(k)=cos((2.*fk-1.)*pi/(2.*fn))
        w(k)=pi/(fn*sqrt(1.-om2*x(k)**2))
   10 continue
      return
      end
      subroutine test09 ( )

c*********************************************************************72
c
      dimension a(199),b(199),fnu(200),alpha(100),beta(100),s(100),
     *s0(200),s1(200),s2(200),alphc(100),betc(100)
      double precision dsigma,da(199),db(199),dnu(200),dalpha(100),
     *dbeta(100),ds(100),ds0(200),ds1(200),ds2(200),dalphc(100),
     *dbetc(100)
c
c This test recomputes the results of  test2  for  sigma=.5  by applying
c the routine  chri  with  iopt=1, x=0  to the weight function with
c parameter  sigma=-.5. Printed are the relative discrepancies in both
c single and double precision between these results and those obtained
c in  test 2  by the modified Chebyshev algorithm. The test is embedded
c in the routine  test2, from which all print statements have been
c removed.
c
      logical modmom,intexp

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST09'

      modmom=.true.
c
c Generate the recursion coefficients for the polynomials defining the
c modified resp. ordinary moments.
c
      if(modmom) then
        n=100
        ndm1=2*n-1
        call recur(ndm1,2,0.,0.,a,b,ierr)
        call drecur(ndm1,2,0.d0,0.d0,da,db,iderr)
      else
        n=12
        ndm1=2*n-1
        do 10 k=1,ndm1
          a(k)=0.
          b(k)=0.
          da(k)=0.d0
          db(k)=0.d0
   10   continue
      end if
      do 30 is=1,3
        dsigma=-.5d0+.5d0*dble(is-1)
        sigma=sngl(dsigma)
        if(is.eq.2) then
          intexp=.true.
        else
          intexp=.false.
        end if
c
c Compute the modified resp. ordinary moments using Eqs. (3.12) and
c (3.11) of the companion paper. On machines with limited exponent
c range, some of the high-order modified moments may underflow, without
c this having any deteriorating effect on the accuracy.
c
        call mm_09_r4(n,modmom,intexp,sigma,fnu)
        call mm_09_r8(n,modmom,intexp,dsigma,dnu)
c
c Compute the desired recursion coefficients by means of the modified
c Chebyshev algorithm; for the latter, see, e.g., Section 2.4 of
c W. Gautschi, On generating orthogonal polynomials'', SIAM J. Sci.
c Statist. Comput. 3, 1982, 289-317.
c
        call cheb_r4(n,a,b,fnu,alpha,beta,s,ierr,s0,s1,s2)
c
c On machines with limited single-precision exponent range, the routine
c CHEB_R4 may generate an underflow exception, which however is harmless
c and can be ignored.
c
        call cheb_r8(n,da,db,dnu,dalpha,dbeta,ds,iderr,ds0,ds1,ds2)
c
c Up to this point the code is identical with the one of  test2.
c
        if(is.eq.1) then
          write(*,1) ierr,iderr
    1     format(/1x,'ierr in cheb_r4 = ',i4,' iderr in cheb_r8 = ',i4/)
          if(ierr.ne.0) then
            nc=abs(ierr)
          else
            nc=n
          end if
          if(iderr.ne.0) then
            ncd=abs(iderr)
          else
            ncd=n
          end if
c
c Compute the desired recursion coefficients by a modification
c algorithm.
c
          nm1=nc-1
          call chri(nm1,1,alpha,beta,0.,0.,0.,0.,alphc,betc,ierr)
          nm1=ncd-1
          call chri_r8(nm1,1,dalpha,dbeta,0.d0,0.d0,0.d0,0.d0,dalphc,
     *      dbetc,iderr)
        end if
        if(is.eq.3) then
          write(*,2)
    2     format(/1x,'test of the results for sigma=1/2'/)
          np=nc
          if(ncd.lt.nc) np=ncd
          nm1=np-1
c
c Compute and print the relative discrepancies between the results of
c the modified Chebyshev algorithm and the modification algorithm.
c
          write(*,3)
    3     format(3x,'k',2x,'err alpha',4x,'err beta',12x,'err dalpha',
     *      2x,'err dbeta'/)
          do 15 k=1,nm1
            km1=k-1
            errac=abs(alpha(k)-alphc(k))/alpha(k)
            errbc=abs(beta(k)-betc(k))/beta(k)
            errdac=sngl(dabs(dalpha(k)-dalphc(k))/dalpha(k))
            errdbc=sngl(dabs(dbeta(k)-dbetc(k))/dbetc(k))
            write(*,4) km1,errac,errbc,errdac,errdbc
    4       format(1x,i3,2e12.4,9x,2e12.4)
   15     continue
          write(*,5)
    5     format(/1x,'end of test'/)
        end if
c
c The rest of the code is essentially the same as the corresponding
c piece of code in  test2  with all print statements removed.
c
        eamax=0.
        ebmax=0.
        do 20 k=1,n
          km1=k-1
          erra=sngl(dabs(dble(alpha(k))-dalpha(k))/dalpha(k))
          errb=sngl(dabs(dble(beta(k))-dbeta(k))/dbeta(k))
          if(erra.gt.eamax) then
            eamax=erra
            kamax=km1
          end if
          if(errb.gt.ebmax) then
            ebmax=errb
            kbmax=km1
          end if
   20   continue
   30 continue
      return
      end
      subroutine mm_09_r4(n,modmom,intexp,sigma,fnu)

c*********************************************************************72
c
c
cc MM_09_R4 generates the first  2*n  modified moments (if modmom=.true.)
c relative to shifted monic Legendre polynomials, using Eq. (3.12) of
c the companion paper, and the first  2*n  ordinary moments (if modmom
c =.false.) by Eq. (3.11), of the weight function
c
c          (x**sigma)*ln(1/x)  on (0,1],   sigma > -1,
c
c for sigma an integer (if intexp=.true.) or a real number (if intexp
c =.false.). In either case, the input variable sigma is of type real.
c
      dimension fnu(*)
      logical modmom,intexp
c
c The array  fnu  is assumed to have dimension  2*n.
c
      nd=2*n
      sigp1=sigma+1.
      if(modmom) then
        isigma=int(sigma)
        isigp1=isigma+1
        isigp2=isigma+2
        isigp3=isigma+3
        if(intexp .and. isigp1.lt.nd) then
          kmax=isigp1
        else
          kmax=nd
        end if
        c=1.
        do 20 k=1,kmax
          km1=k-1
          fk=real(k)
          p=1.
          s=1./sigp1
          if(kmax.gt.1) then
            do 10 i=1,km1
              fi=real(i)
              p=(sigp1-fi)*p/(sigp1+fi)
              s=s+1./(sigp1+fi)-1./(sigp1-fi)
   10       continue
          end if
          fnu(k)=c*s*p/sigp1
          c=fk*c/(4.*fk-2.)
   20   continue
        if(.not.intexp .or. isigp1.ge.nd) return
        q=-.5
        if(isigma.gt.0) then
          do 30 iq=1,isigma
            fiq=real(iq)
            q=fiq*fiq*q/((2.*fiq+1.)*(2.*fiq+2.))
   30     continue
        end if
        fnu(isigp2)=c*q
        if(isigp2.eq.nd) return
        do 40 k=isigp3,nd
          km1=k-1
          fkm1=real(km1)
          fnu(k)=-fkm1*(fkm1-sigp1)*fnu(km1)/((4.*fkm1-2.)*
     *      (fkm1+sigp1))
   40   continue
        return
      else
        do 50 k=1,nd
          fkm1=real(k-1)
          fnu(k)=(1./(sigp1+fkm1))**2
   50   continue
      end if
      end
      subroutine mm_09_r8(n,modmom,intexp,dsigma,dnu)

c*********************************************************************72
c
c
cc MM_09_R8 is a double-precision version of the routine MM_09_R4.
c
      double precision dsigma,dnu(*),dsigp1,dc,dk,dp,ds,di,dq,diq,dkm1
      logical modmom,intexp
c
c The array  dnu  is assumed to have dimension  2*n.
c
      nd=2*n
      dsigp1=dsigma+1.d0
      if(modmom) then
        isigma=idint(dsigma)
        isigp1=isigma+1
        isigp2=isigma+2
        isigp3=isigma+3
        if(intexp .and. isigp1.lt.nd) then
          kmax=isigp1
        else
          kmax=nd
        end if
        dc=1.d0
        do 20 k=1,kmax
          km1=k-1
          dk=dble(k)
          dp=1.d0
          ds=1.d0/dsigp1
          if(kmax.gt.1) then
            do 10 i=1,km1
              di=dble(i)
              dp=(dsigp1-di)*dp/(dsigp1+di)
              ds=ds+1.d0/(dsigp1+di)-1.d0/(dsigp1-di)
   10       continue
          end if
          dnu(k)=dc*ds*dp/dsigp1
          dc=dk*dc/(4.d0*dk-2.d0)
   20   continue
        if(.not.intexp .or. isigp1.ge.nd) return
        dq=-.5d0
        if(isigma.gt.0) then
          do 30 iq=1,isigma
            diq=dble(iq)
            dq=diq*diq*dq/((2.d0*diq+1.d0)*(2.d0*diq+2.d0))
   30     continue
        end if
        dnu(isigp2)=dc*dq
        if(isigp2.eq.nd) return
        do 40 k=isigp3,nd
          km1=k-1
          dkm1=dble(km1)
          dnu(k)=-dkm1*(dkm1-dsigp1)*dnu(km1)/((4.d0*dkm1-2.d0)*
     *      (dkm1+dsigp1))
   40   continue
        return
      else
        do 50 k=1,nd
          dkm1=dble(k-1)
          dnu(k)=(1.d0/(dsigp1+dkm1))**2
   50   continue
      end if
      end
      subroutine test10 ( )

c*********************************************************************72
c
      dimension a(31),b(31),alpha(31),beta(31),z(11),w(11),e(11),
     *a1(31),b1(31),betap(20,12),erram(12),errbm(12)
      double precision depsma,d1mach,da(31),db(31),dalpha(31),dbeta(31),
     *dz(11),dw(11),de(11),da1(31),db1(31)
c
c epsma and depsma are the machine single and double precision.
c

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST10'

      epsma=r1mach(3)
      depsma=d1mach(3)
c
c This test applies the routines  indp  and  dindp  to generate the
c first 20 recursion coefficients of the induced Legendre polynomials
c pind(k,m)(.), m=0,1,2,...,11, that is, of the polynomials orthogonal
c relative to the weight function
c
c                 [p(m)(x)]**2    on [-1,1],
c
c where  p(m)(.)  is the (monic) Legendre polynomial of degree m.
c (When m=0, then  pind(k,0)(.)=p(k)(.).) The routine also prints the
c absolute and relative errors, respectively, of the alpha- and beta-
c coefficients.
c
      n=20
      do 20 im=1,12
        m=im-1
        npm=n+m
c
c Generate the Legendre recurrence coefficients required in the
c routines  indp  and dindp.
c
        call recur(npm,1,0.,0.,a,b,ierr)
        call drecur(npm,1,0.d0,0.d0,da,db,ierr)
c
c Compute the desired recursion coefficients.
c
        call indp(n,m,a,b,epsma,alpha,beta,ierr,z,w,e,a1,b1)
        call dindp(n,m,da,db,depsma,dalpha,dbeta,ierr,dz,dw,
     *    de,da1,db1)
c
c Compute and print the respective errors.
c
        erram(im)=0.
        errbm(im)=0.
        do 10 k=1,n
          erra=sngl(dabs(dble(alpha(k))-dalpha(k)))
          errb=sngl(dabs((dble(beta(k))-dbeta(k))/dbeta(k)))
          if(erra.gt.erram(im)) erram(im)=erra
          if(errb.gt.errbm(im)) errbm(im)=errb
          betap(k,im)=sngl(dbeta(k))
   10   continue
   20 continue
      do 40 ip=1,3
        ip4=1+4*(ip-1)
        write(*,1) ip4-1,ip4,ip4+1,ip4+2
    1   format(5x,'k',2x,'m=',i1,'  beta(k)',2x,'m=',i1,'  beta(k)',2x,
     *    'm=',i2,' beta(k)',2x,'m=',i2,' beta(k)'/)
        do 30 k=1,n
          km1=k-1
          write(*,2) km1,betap(k,ip4),betap(k,ip4+1),betap(k,ip4+2),
     *      betap(k,ip4+3)
    2     format(1x,i5,4f14.10)
   30   continue
        write(*,3) erram(ip4),erram(ip4+1),erram(ip4+2),erram(ip4+3)
    3   format(/4x,'erra',e12.4,3e14.4)
        write(*,4) errbm(ip4),errbm(ip4+1),errbm(ip4+2),errbm(ip4+3)
    4   format(4x,'errb',e12.4,3e14.4//)
   40 continue
      return
      end
      subroutine indp(n,m,a,b,eps,alpha,beta,ierr,z,w,e,a1,b1)

c*********************************************************************72
c
c
c If  p(m)(.)  denotes the (monic) orthogonal polynomial of degree  m
c relative to the weight function  w(x), then the corresponding m-th
c induced orthogonal polynomials  pind(k,m)(.), k=0,1,2,..., are those
c orthogonal with respect to the weight function
c
c                   (p(m)(x)**2)*w(x).
c
c (For background on induced orthogonal polynomials, including an
c algorithm for generating their recursion coefficients, see W. Gautschi
c and S. Li,A set of orthogonal polynomials induced by a given
c orthogonal polynomial'', Aequationes Math., to appear.) This routine
c obtains the first n recurrence coefficients of the m-th induced
c orthogonal polynomials by an m-fold application of the routine  chri
c with  iopt=7, the shifts taken being, in succession, the zeros of
c p(m)(.).
c
      dimension a(*),b(*),alpha(*),beta(*),z(m),w(m),e(m),
     *a1(*),b1(*)
c
c The arrays  a,b,alpha,beta,a1,b1  are assumed to have dimension  n+m.
c
      npm=n+m
      do 10 k=1,npm
        alpha(k)=a(k)
        beta(k)=b(k)
   10 continue
      if(m.eq.0) return
      call gauss(m,a,b,eps,z,w,ierr,e)
      do 30 imu=1,m
        mi=npm-imu
        do 20 k=1,mi+1
          a1(k)=alpha(k)
          b1(k)=beta(k)
   20   continue
        x=z(imu)
        call chri(mi,7,a1,b1,x,0.,0.,0.,alpha,beta,ierrc)
   30 continue
      return
      end
      subroutine dindp(n,m,da,db,deps,dalpha,dbeta,ierr,dz,dw,de,
     *da1,db1)

c*********************************************************************72
c
c
c This is a double-precision version of the routine  indp.
c
      double precision da(*),db(*),deps,dalpha(*),dbeta(*),
     *dz(m),dw(m),de(m),da1(*),db1(*),dx
c
c The arrays  da,db,dalpha,dbeta,da1,db1  are assumed to have
c dimension  n+m.
c
      npm=n+m
      do 10 k=1,npm
        dalpha(k)=da(k)
        dbeta(k)=db(k)
   10 continue
      if(m.eq.0) return
      call dgauss(m,da,db,deps,dz,dw,ierr,de)
      do 30 imu=1,m
        mi=npm-imu
        do 20 k=1,mi+1
          da1(k)=dalpha(k)
          db1(k)=dbeta(k)
   20   continue
        dx=dz(imu)
        call chri_r8(mi,7,da1,db1,dx,0.d0,0.d0,0.d0,dalpha,
     *    dbeta,ierrc)
   30 continue
      return
      end
      subroutine test11 ( )

c*********************************************************************72
c
      complex rho,rold,z,e
      dimension xx(5),rr(5),a(500),b(500),alpha(40),beta(40),alphc(40),
     *betc(40),fnu(80),rho(80),rold(80),s(40),s0(80),s1(80),s2(80),
     *alphr(40),betr(40),alphcr(40),betcr(40)
      double precision d1mach,depsma,deps,dal,dbe,dhi,da(800),db(800),
     *dx,dy,dalpha(40),dbeta(40),dnu(80),drhor(80),drhoi(80),
     *droldr(80),droldi(80),ds(40),ds0(80),ds1(80),ds2(80),dhr,
     *dalphc(40),dbetc(40),dalphr(40),dbetr(40),dalcr(40),dbecr(40)
      data xx/1.001,1.01,1.04,1.07,1.1/
      data rr/1.05,1.1625,1.275,1.3875,1.5/
c
c This test is to illustrate the dissimilar performance of the routines
c chri  and  gchri  in the case of division of the Jacobi weight
c function  w(t;alj,bej)  with parameters  alj,bej  by either a linear
c divisor  t-x  or a quadratic divisor  (t-x)**2 + y**2 . In either
c case, the parameters selected are  alj=-.8(.4).8, bej=alj(.4).8.  In
c the former case, x = -1.001, -1.01, -1.04, -1.07 and -1.1, whereas in
c the latter case,  x and  y  are chosen to lie, regularly spaced, on
c the upper half of an ellipse with foci at +1 and -1 and sum of the
c semiaxes equal to  rho = 1.05, 1.1625, 1.275, 1.3875 and 1.5. The
c routines are run in both single and double precision with n=40, the
c results of the latter being used to calculate, and print, the maximum
c absolute and relative error of the single-precision alpha- and beta-
c coefficients, respectively. Also printed are the starting recurrence
c indexes required in the backward recurrence schemes of  gchri,dgchri
c to achieve single- resp. double-precision accuracy. This information
c is contained in the first line of each 3-line block of the output,
c where in the case of quadratic divisors only average values (averaged
c over the upper half of the respective ellipse) are shown. The second
c and third line of each 3-line block display the maximum
c reconstruction error'', that is, the maximum errors in the alpha's
c and beta's if the coefficients produced by  gchri,chri  and  dgchri,
c chri_r8  are fed back to the routines  chri  and  chri_r8  with  iopt=1
c to recover the original recursion coefficients in single and double
c precision.
c
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST11'

      write(*,1)
    1 format(/)
      epsma=r1mach(3)
      depsma=d1mach(3)
c
c epsma and depsma are the machine single and double precision.
c
      n=40
      np1=n+1
      nm1=n-1
      nd=2*n
      ndm1=nd-1
      numax=500
      numaxd=800
      eps=10.*epsma
      deps=100.d0*depsma
      epsd=sngl(deps)
      ipoly=6
      do 70 ial=1,5
        al=-.8+.4*real(ial-1)
        dal=dble(al)
        ibemax=6-ial
        do 60 ibe=1,ibemax
          be=al+.4*real(ibe-1)
          dbe=dble(be)
          write(*,2) al,be
    2     format(///1x,'al = ',f6.2,'  be = ',f6.2//)
          hi=0.
          dhi=0.d0
c
c Generate the Jacobi recurrence coefficients to be used in the
c backward recurrence algorithm of the routines  gchri  and  dgchri.
c
          call recur(numax,ipoly,al,be,a,b,ierr)
          call drecur(numaxd,ipoly,dal,dbe,da,db,ierrd)
          write(*,3)
    3     format(30x,'gchri',20x,'chri')
          write(*,4)
    4     format(5x,'x',5x,'nu0',2x,'nud0',4x,'erra',8x,'errb',10x,
     *      'erra',7x,'errb'/)
          do 20 ix=1,5
            x=-xx(ix)
            dx=dble(x)
            y=0.
            dy=0.d0
            z=cmplx(x,y)
c
c Compute the starting index for backward recurrence.
c
            nu0=nu0jac(ndm1,z,eps)
            nu0d=nu0jac(ndm1,z,epsd)
c
c Generate the recurrence coefficients for the Jacobi weight function
c divided by a linear divisor, using the routines  gchri,dgchri.
c
            call gchri(n,1,nu0,numax,eps,a,b,x,y,alpha,beta,nu,ierrg,
     *        ierrc,fnu,rho,rold,s,s0,s1,s2)
c
c On machines with limited single-precision exponent range, the routine
c CHEB_R4 used in  gchri  may have generated an underflow exception,
c which however is harmless and can be ignored.
c
            call dgchri(n,1,nu0d,numaxd,deps,da,db,dx,dy,dalpha,dbeta,
     *        nud,ierrgd,ierrcd,dnu,drhor,drhoi,droldr,droldi,ds,ds0,
     *        ds1,ds2)
            if(ierrg.ne.0 .or. ierrc.ne.0 .or. ierrgd.ne.0 .or.ierrcd
     *        .ne.0) then
              write(*,5) ierrg,ierrgd,al,be,x
    5         format(/1x,'ierrg in gchri = ',i4,' ierrg in dgchri = ',
     *          i4,' for al = ',f6.2,' be = ',f6.2,' x = ',f7.4)
              write(*,6) ierrc,ierrcd,al,be,x
    6         format(1x,'ierrc in gchri = ',i4,' ierrc in dgchri = ',
     *          i4,' for al = ',f6.2,' be = ',f6.2,' x = ',f7.4/)
              goto 20
            end if
c
c Generate the recurrence coefficients for the Jacobi weight function
c divided by a linear divisor, using the routines  chri, chri_r8.
c
            hr=real(rho(1))
            dhr=drhor(1)
            call chri(n,4,a,b,x,y,hr,hi,alphc,betc,ierr)
            call chri_r8(n,4,da,db,dx,dy,dhr,dhi,dalphc,dbetc,ierr)
c
c Do the reconstruction.
c
            call chri(nm1,1,alpha,beta,x,y,0.,0.,alphr,betr,ierr)
            call chri_r8(nm1,1,dalpha,dbeta,dx,dy,0.d0,0.d0,dalphr,
     &        dbetr, ierr)
            call chri(nm1,1,alphc,betc,x,y,0.,0.,alphcr,betcr,ierr)
            call chri_r8(nm1,1,dalphc,dbetc,dx,dy,0.d0,0.d0,dalcr,dbecr,
     *        ierr)
c
c Compute and print the maximum errors.
c
            erragm=0.
            errbgm=0.
            erracm=0.
            errbcm=0.
            erram=0.
            errbm=0.
            errdam=0.
            errdbm=0.
            eracrm=0.
            erbcrm=0.
            edacrm=0.
            edbcrm=0.
            do 10 k=1,n
              km1=k-1
              errag=abs(sngl(dble(alpha(k))-dalpha(k)))
              errbg=abs(sngl((dble(beta(k))-dbeta(k))/dbeta(k)))
              errac=abs(sngl(dble(alphc(k))-dalphc(k)))
              errbc=abs(sngl((dble(betc(k))-dbetc(k))/dbetc(k)))
              if(k.lt.n) then
                erra=abs(alphr(k)-a(k))
                errb=abs((betr(k)-b(k))/b(k))
                errda=abs(sngl(dalphr(k)-da(k)))
                errdb=abs(sngl((dbetr(k)-db(k))/db(k)))
                eracr=abs(alphcr(k)-a(k))
                erbcr=abs((betcr(k)-b(k))/b(k))
                edacr=abs(sngl(dalcr(k)-da(k)))
                edbcr=abs(sngl((dbecr(k)-db(k))/db(k)))
                if(erra.gt.erram) erram=erra
                if(errb.gt.errbm) errbm=errb
                if(errda.gt.errdam) errdam=errda
                if(errdb.gt.errdbm) errdbm=errdb
                if(eracr.gt.eracrm) eracrm=eracr
                if(erbcr.gt.erbcrm) erbcrm=erbcr
                if(edacr.gt.edacrm) edacrm=edacr
                if(edbcr.gt.edbcrm) edbcrm=edbcr
              end if
              if(errag.gt.erragm) erragm=errag
              if(errbg.gt.errbgm) errbgm=errbg
              if(errac.gt.erracm) erracm=errac
              if(errbc.gt.errbcm) errbcm=errbc
   10       continue
            write(*,7) x,nu0,nu0d,erragm,errbgm,erracm,errbcm
    7       format(/1x,f7.4,2i6,2e12.4,2x,2e12.4)
            if(ix.eq.1) then
              write(*,8) erram,errbm,errdam,errdbm
    8         format(11x,'reconstr.',2e12.4,2x,2e12.4)
              write(*,9) eracrm,erbcrm,edacrm,edbcrm
    9         format(12x,'errors',2x,2e12.4,2x,2e12.4)
            else
              write(*,11) erram,errbm,errdam,errdbm
   11         format(20x,2e12.4,2x,2e12.4)
              write(*,11) eracrm,erbcrm,edacrm,edbcrm
            end if
   20     continue
          write(*,1)
          ndiv=20
          ndivm1=ndiv-1
          fndiv=real(ndiv)
          fndm1=real(ndivm1)
          pi=4.*atan(1.)
          write(*,3)
          write(*,12)
   12     format(4x,'rho',4x,'nu0',2x,'nud0',4x,'erra',8x,'errb',10x,
     *      'erra',7x,'errb'/)
          do 50 ir=1,5
            r=rr(ir)
            agmv=0.
            bgmv=0.
            acmv=0.
            bcmv=0.
            amv=0.
            bmv=0.
            damv=0.
            dbmv=0.
            acrmv=0.
            bcrmv=0.
            dacrmv=0.
            dbcrmv=0.
            nu0v=0
            nu0dv=0
            do 40 ith=1,ndivm1
c
c Generate the points on the ellipse.
c
              theta=pi*real(ith)/fndiv
              e=cmplx(cos(theta),sin(theta))
              z=.5*(r*e+1./(r*e))
              x=real(z)
              y=aimag(z)
              dx=dble(x)
              dy=dble(y)
c
c Compute the starting index for backward recurrence.
c
              nu0=nu0jac(ndm1,z,eps)
              nu0d=nu0jac(ndm1,z,epsd)
c
c Generate the recurrence coefficients for the Jacobi weight function
c divided by a quadratic divisor, using the routines  gchri,dgchri.
c
              call gchri(n,2,nu0,numax,eps,a,b,x,y,alpha,beta,nu,ierrg,
     *          ierrc,fnu,rho,rold,s,s0,s1,s2)
              call dgchri(n,2,nu0d,numaxd,deps,da,db,dx,dy,dalpha,dbeta,
     *          nud,ierrgd,ierrcd,dnu,drhor,drhoi,droldr,droldi,ds,ds0,
     *          ds1,ds2)
              if(ierrg.ne.0 .or.ierrc.ne.0 .or. ierrgd.ne.0 .or. ierrcd
     *          .ne.0) then
                write(*,5) ierrg,ierrgd,al,be,x
                write(*,6) ierrc,ierrcd,al,be,x
                goto 40
              end if
              nu0v=nu0v+nu0
              nu0dv=nu0dv+nu0d
c
c Generate the recurrence coefficients for the Jacobi weight function
c divided by a quadratic divisor, using the routines  chri, chri_r8.
c
              hr=real(rho(1))
              hi=aimag(rho(1))
              dhr=drhor(1)
              dhi=drhoi(1)
              call chri(n,5,a,b,x,y,hr,hi,alphc,betc,ierr)
              call chri_r8(n,5,da,db,dx,dy,dhr,dhi,dalphc,dbetc,ierr)
c
c Do the reconstruction.
c
              call chri(nm1,2,alpha,beta,x,y,0.,0.,alphr,betr,ierr)
              call chri_r8(nm1,2,dalpha,dbeta,dx,dy,0.d0,0.d0,dalphr,
     *          dbetr,ierr)
              call chri(nm1,2,alphc,betc,x,y,0.,0.,alphcr,betcr,ierr)
              call chri_r8(nm1,2,dalphc,dbetc,dx,dy,0.d0,0.d0,dalcr,
     *          dbecr,ierr)
c
c Compute and print the maximum average errors.
c
              erragm=0.
              errbgm=0.
              erracm=0.
              errbcm=0.
              erram=0.
              errbm=0.
              errdam=0.
              errdbm=0.
              eracrm=0.
              erbcrm=0.
              edacrm=0.
              edbcrm=0.
              do 30 k=1,n
                km1=k-1
                errag=abs(sngl(dble(alpha(k))-dalpha(k)))
                errbg=abs(sngl((dble(beta(k))-dbeta(k))/dbeta(k)))
                errac=abs(sngl(dble(alphc(k))-dalphc(k)))
                errbc=abs(sngl((dble(betc(k))-dbetc(k))/dbetc(k)))
                if(k.lt.n) then
                  erra=abs(alphr(k)-a(k))
                  errb=abs((betr(k)-b(k))/b(k))
                  errda=abs(sngl(dalphr(k)-da(k)))
                  errdb=abs(sngl((dbetr(k)-db(k))/db(k)))
                  eracr=abs(alphcr(k)-a(k))
                  erbcr=abs((betcr(k)-b(k))/b(k))
                  edacr=abs(sngl(dalcr(k)-da(k)))
                  edbcr=abs(sngl((dbecr(k)-db(k))/db(k)))
                  if(erra.gt.erram) erram=erra
                  if(errb.gt.errbm) errbm=errb
                  if(errda.gt.errdam) errdam=errda
                  if(errdb.gt.errdbm) errdbm=errdb
                  if(eracr.gt.eracrm) eracrm=eracr
                  if(erbcr.gt.erbcrm) erbcrm=erbcr
                  if(edacr.gt.edacrm) edacrm=edacr
                  if(edbcr.gt.edbcrm) edbcrm=edbcr
                end if
                if(errag.gt.erragm) erragm=errag
                if(errbg.gt.errbgm) errbgm=errbg
                if(errac.gt.erracm) erracm=errac
                if(errbc.gt.errbcm) errbcm=errbc
   30         continue
              agmv=agmv+erragm
              bgmv=bgmv+errbgm
              acmv=acmv+erracm
              bcmv=bcmv+errbcm
              amv=amv+erram
              bmv=bmv+errbm
              damv=damv+errdam
              dbmv=dbmv+errdbm
              acrmv=acrmv+eracrm
              bcrmv=bcrmv+erbcrm
              dacrmv=dacrmv+edacrm
              dbcrmv=dbcrmv+edbcrm
   40       continue
            nu0=real(nu0v)/fndm1
            nu0d=real(nu0dv)/fndm1
            erragm=agmv/fndm1
            errbgm=bgmv/fndm1
            erracm=acmv/fndm1
            errbcm=bcmv/fndm1
            erram=amv/fndm1
            errbm=bmv/fndm1
            errdam=damv/fndm1
            errdbm=dbmv/fndm1
            eracrm=acrmv/fndm1
            erbcrm=bcrmv/fndm1
            edacrm=dacrmv/fndm1
            edbcrm=dbcrmv/fndm1
            if(ir.eq.1) then
              write(*,7) r,nu0,nu0d,erragm,errbgm,erracm,errbcm
              write(*,8) erram,errbm,errdam,errdbm
              write(*,9) eracrm,erbcrm,edacrm,edbcrm
            else
              write(*,7) r,nu0,nu0d,erragm,errbgm,erracm,errbcm
              write(*,11) erram,errbm,errdam,errdbm
              write(*,11) eracrm,erbcrm,edacrm,edbcrm
            end if
   50     continue
   60   continue
   70 continue

      return
      end
