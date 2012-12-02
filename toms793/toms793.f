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
c    D1MACH ( 1) = B^(EMIN-1), the smallest positive magnitude.
c    D1MACH ( 2) = B^EMAX*(1 - B^(-T)), the largest magnitude.
c    D1MACH ( 3) = B^(-T), the smallest relative spacing.
c    D1MACH ( 4) = B^(1-T), the largest relative spacing.
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
      subroutine dabmod(n,ncapmd,mcap,deps,irout,da,db,dxir,dxii,is,
     +                  dalpha,dbeta,ncapd,kountd,ierrd,dbe,dx,dw,de,
     +                  dp0,dp1,dp2)

c*********************************************************************72
c
cc DABMOD ???
c
c  Modified:
c
c    08 November 2012
c
c  Author:
c
c    Walter Gautschi
c
c  Reference:
c
c    Walter Gautschi,
c    Algorithm 793: GQRAT, Gauss Quadrature for Rational Functions,
c    ACM Transactions on Mathematical Software,
c    Volume 25, Number 2, pages 213-239, June 1999.
c
c  Parameters:
c
      implicit none

      double precision deps
      integer ierrd,irout,kountd,mcap,n,ncapd,ncapmd
c     ..
c     .. array arguments ..  the upper bounds dxii, dxir and is are all
c                            mcap, but this doesn't work in fortran 77
c                            if mcap = 0.
      double precision da(ncapmd),dalpha(n),db(ncapmd),dbe(n),dbeta(n),
     +                 de(ncapmd),dp0(ncapmd),dp1(ncapmd),dp2(ncapmd),
     +                 dw(ncapmd),dx(ncapmd),dxii(*),dxir(*)
      integer is(*)
      double precision depsm,dp
      integer ied,ierrgd,imu,incap,k,mu
c     ..
c     .. external functions ..
      double precision d1mach
      external d1mach
c     ..
c     .. external subroutines ..
      external dgauss,dlancz,dsti
c     ..
c     .. intrinsic functions ..
      intrinsic abs
c     ..
      integer i1mach, nout
      nout = i1mach(2)
      ierrd = 0
      depsm = d1mach(3)
      incap = 1

      do k = 1,n
        dalpha(k) = da(k)
        dbeta(k) = db(k)
      end do

      if (mcap.gt.0) then

          kountd = -1

          do k = 1,n
            dbeta(k) = 0.d0
          end do

          ncapd = (2*n-1)/2

   30     continue

          do k = 1,n
            dbe(k) = dbeta(k)
          end do

          kountd = kountd + 1
          if (kountd.gt.1) incap = 2** (kountd/5)*n
          ncapd = ncapd + incap
          if (ncapd.gt.ncapmd) then
              ierrd = ncapmd
              return

          end if

          call dgauss(ncapd,da,db,depsm,dx,dw,ierrgd,de)
          if (ierrgd.ne.0) then
              write (nout,fmt=9000) ierrgd
              ierrd = 1
              return

          end if

          do 60 k = 1,ncapd
              dp = 1.d0
              imu = 0
              do 50 mu = 1,mcap
                  if (imu.eq.0) then
                      if (dxii(mu).eq.0.d0) then
                          dp = ((1.d0+dxir(mu)*dx(k))**is(mu))*dp

                      else
                          dp = (((1.d0+dxir(mu)*dx(k))**2+
     +                         (dxii(mu)*dx(k))**2)**is(mu))*dp
                          imu = 1
                      end if

                  else
                      imu = 0
                  end if

   50         continue
              dw(k) = dw(k)/dp
c          dw(k)=dsqrt(1.d0+.5d0*theta*dx(k))*dw(k)/dp
   60     continue
          if (irout.eq.1) then
              call dsti(n,ncapd,dx,dw,dalpha,dbeta,ied,dp0,dp1,dp2)
              if (ied.ne.0) then
                  ierrd = 2
                  return

              end if

          else
              call dlancz(n,ncapd,dx,dw,dalpha,dbeta,ied,dp0,dp1)
              if (ied.ne.0) then
                  ierrd = 2
                  return

              end if

          end if

          do k = 1,n
            if (abs(dbeta(k)-dbe(k)).gt.deps*dbeta(k)) then
              go to 30
            end if
          end do

      end if

      return

 9000 format (1x,'ierrgd in dgauss=',i5)
      end
      subroutine dgauss(n,dalpha,dbeta,deps,dzero,dweigh,ierr,de)

c*********************************************************************72
c
cc DGAUSS ???
c
c  Modified:
c
c    08 November 2012
c
c  Author:
c
c    Walter Gautschi
c
c  Reference:
c
c    Walter Gautschi,
c    Algorithm 793: GQRAT, Gauss Quadrature for Rational Functions,
c    ACM Transactions on Mathematical Software,
c    Volume 25, Number 2, pages 213-239, June 1999.
c
c  Parameters:
c
      implicit none

      double precision deps
      integer ierr,n

      double precision dalpha(n),dbeta(n),de(n),dweigh(n),dzero(n)
      double precision db,dc,df,dg,dp,dr,ds
      integer i,ii,j,k,l,m,mml
      intrinsic abs,sign,sqrt

      if (n.lt.1) then
        ierr = -1
        return
      end if

      ierr = 0
      dzero(1) = dalpha(1)

      if (dbeta(1).lt.0.d0) then
        ierr = -2
        return
      end if

      dweigh(1) = dbeta(1)
      if (n.eq.1) return
      dweigh(1) = 1.d0
      de(n) = 0.d0

      do k = 2,n
        dzero(k) = dalpha(k)
        if (dbeta(k).lt.0.d0) then
          ierr = -2
          return
        end if
        de(k-1) = sqrt(dbeta(k))
        dweigh(k) = 0.d0
      end do

      do 80 l = 1,n
          j = 0
   20     do 30 m = l,n
              if (m.eq.n) go to 40
              if (abs(de(m)).le.deps* (abs(dzero(m))+abs(dzero(m+
     +            1)))) go to 40
   30     continue
   40     dp = dzero(l)
          if (m.eq.l) go to 80
          if (j.eq.30) go to 120
          j = j + 1
          dg = (dzero(l+1)-dp)/ (2.d0*de(l))
          dr = sqrt(dg*dg+1.d0)
          dg = dzero(m) - dp + de(l)/ (dg+sign(dr,dg))
          ds = 1.d0
          dc = 1.d0
          dp = 0.d0
          mml = m - l
          do 70 ii = 1,mml
              i = m - ii
              df = ds*de(i)
              db = dc*de(i)
              if (abs(df).lt.abs(dg)) go to 50
              dc = dg/df
              dr = sqrt(dc*dc+1.d0)
              de(i+1) = df*dr
              ds = 1.d0/dr
              dc = dc*ds
              go to 60

   50         ds = df/dg
              dr = sqrt(ds*ds+1.d0)
              de(i+1) = dg*dr
              dc = 1.d0/dr
              ds = ds*dc
   60         dg = dzero(i+1) - dp
              dr = (dzero(i)-dg)*ds + 2.d0*dc*db
              dp = ds*dr
              dzero(i+1) = dg + dp
              dg = dc*dr - db
              df = dweigh(i+1)
              dweigh(i+1) = ds*dweigh(i) + dc*df
              dweigh(i) = dc*dweigh(i) - ds*df
   70     continue
          dzero(l) = dzero(l) - dp
          de(l) = dg
          de(m) = 0.d0
          go to 20

   80 continue
      do 100 ii = 2,n
          i = ii - 1
          k = i
          dp = dzero(i)
          do 90 j = ii,n
              if (dzero(j).ge.dp) go to 90
              k = j
              dp = dzero(j)
   90     continue
          if (k.eq.i) go to 100
          dzero(k) = dzero(i)
          dzero(i) = dp
          dp = dweigh(i)
          dweigh(i) = dweigh(k)
          dweigh(k) = dp
  100 continue
      do 110 k = 1,n
          dweigh(k) = dbeta(1)*dweigh(k)*dweigh(k)
  110 continue
      return

  120 ierr = l

      return
      end
      subroutine dgchrs(n,iopt,da,db,dx,dhp,dhn,dalpha,dbeta,ierrd)

c*********************************************************************72
c
cc DGCHRS ???
c
c  Modified:
c
c    08 November 2012
c
c  Author:
c
c    Walter Gautschi
c
c  Reference:
c
c    Walter Gautschi,
c    Algorithm 793: GQRAT, Gauss Quadrature for Rational Functions,
c    ACM Transactions on Mathematical Software,
c    Volume 25, Number 2, pages 213-239, June 1999.
c
c  Parameters:
c
      implicit none

      double precision dhn,dhp,dx
      integer ierrd,iopt,n
c     ..
c     .. array arguments ..
      double precision da(n),dalpha(n),db(n),dbeta(n)
c     ..
c     .. local scalars ..
      double precision de,deh,dq,dqh
      integer k
c     ..
      ierrd = 0

      if (iopt.eq.1) then

          dalpha(1) = dx - db(1)/dhp
          dbeta(1) = -dhp
          dq = -db(1)/dhp
          do k = 2,n
              de = da(k-1) - dx - dq
              dbeta(k) = dq*de
              dq = db(k)/de
              dalpha(k) = dq + de + dx
          end do

      else if (iopt.eq.2) then

          dalpha(1) = dx* (dhp+dhn)/ (dhp-dhn)
          dbeta(1) = - (dhp-dhn)/ (2.d0*dx)
          dq = -db(1)/dhp
          dqh = -dhp/dbeta(1)
          de = 0.d0
          do 20 k = 2,n
              deh = dq + de + 2.d0*dx - dqh
              dbeta(k) = dqh*deh
              de = da(k-1) - dx - dq
              dqh = dq*de/deh
              dalpha(k) = dqh + deh - dx
              dq = db(k)/de
   20     continue

      else
          ierrd = 1
      end if

      return
      end
      subroutine dgqrat(n,mcap,dalpha,dbeta,dxir,dxii,is,dzg,dwg,dconst,
     +                  ierrd,de)

c*********************************************************************72
c
cc DGQRAT ???
c
c  Modified:
c
c    08 November 2012
c
c  Author:
c
c    Walter Gautschi
c
c  Reference:
c
c    Walter Gautschi,
c    Algorithm 793: GQRAT, Gauss Quadrature for Rational Functions,
c    ACM Transactions on Mathematical Software,
c    Volume 25, Number 2, pages 213-239, June 1999.
c
c  Parameters:
c
      implicit none

      double precision dconst
      integer ierrd,mcap,n
c     ..
c     .. array arguments ..  the upper bounds dxii, dxir and is are all
c                            mcap, but this doesn't work in fortran 77
c                            if mcap = 0.
      double precision dalpha(n+1),dbeta(n+1),de(n),dwg(n),dxii(*),
     +                 dxir(*),dzg(n)
      integer is(*)
c     ..
c     .. local scalars ..
      double precision depsm,dp
      integer ierrgd,imu,k,mu
c     ..
c     .. external functions ..
      double precision d1mach
      external d1mach
c     ..
c     .. external subroutines ..
      external dgauss
c     ..
c     .. intrinsic functions ..
      intrinsic dble
c     ..
      integer i1mach, nout
      nout = i1mach(2)
      ierrd = 0
      depsm = d1mach(3)
      call dgauss(n,dalpha,dbeta,depsm,dzg,dwg,ierrgd,de)
      if (ierrgd.ne.0) then
          write (nout,fmt=9000) ierrgd
          ierrd = 1
          return

      end if

      dconst = dbeta(1)
      do 20 k = 1,n
          dconst = dbeta(k+1)*dconst/dble(2*k* (2*k-1))
          if (mcap.gt.0) then
              dp = 1.d0
              imu = 0
              do 10 mu = 1,mcap
                  if (imu.eq.0) then
                      if (dxii(mu).eq.0.d0) then
                          dp = ((1.d0+dxir(mu)*dzg(k))**is(mu))*dp

                      else
                          dp = (((1.d0+dxir(mu)*dzg(k))**2+
     +                         (dxii(mu)*dzg(k))**2)**is(mu))*dp
                          imu = 1
                      end if

                  else
                      imu = 0
                  end if

   10         continue
              dwg(k) = dp*dwg(k)
          end if

   20 continue
      return

 9000 format (1x,'ierrgd in dgauss=',i5)
      end
      subroutine dknum(n,nu0,numax,dx,dy,deps,da,db,drhor,drhoi,nu,ierr,
     +                 droldr,droldi)

c*********************************************************************72
c
cc DKNUM ???
c
c  Modified:
c
c    08 November 2012
c
c  Author:
c
c    Walter Gautschi
c
c  Reference:
c
c    Walter Gautschi,
c    Algorithm 793: GQRAT, Gauss Quadrature for Rational Functions,
c    ACM Transactions on Mathematical Software,
c    Volume 25, Number 2, pages 213-239, June 1999.
c
c  Parameters:
c
      implicit none
c
c The arrays  drhor,drhoi,droldr,droldi  are assumed to have
c dimension  n+1.
c
      double precision deps,dx,dy
      integer ierr,n,nu,nu0,numax
c     ..
c     .. array arguments ..
      double precision da(numax),db(numax),drhoi(*),drhor(*),droldi(*),
     +                 droldr(*)
c     ..
c     .. local scalars ..
      double precision dden,dri,drr,dt
      integer j,j1,k,np1
c     ..
      ierr = 0
      np1 = n + 1

      if (nu0.gt.numax) then
          ierr = nu0
          return
      end if

      if (nu0.lt.np1) nu0 = np1
      nu = nu0 - 5
      do k = 1,np1
          drhor(k) = 0.d0
          drhoi(k) = 0.d0
      end do

   20 continue

      nu = nu + 5
      if (nu.gt.numax) then
          ierr = numax
          go to 60
      end if

      do 30 k = 1,np1
          droldr(k) = drhor(k)
          droldi(k) = drhoi(k)
   30 continue
      drr = 0.d0
      dri = 0.d0
      do 40 j = 1,nu
          j1 = nu - j + 1
          dden = (dx-da(j1)-drr)**2 + (dy-dri)**2
          drr = db(j1)* (dx-da(j1)-drr)/dden
          dri = -db(j1)* (dy-dri)/dden
          if (j1.le.np1) then
              drhor(j1) = drr
              drhoi(j1) = dri
          end if

   40 continue
      do 50 k = 1,np1
          if ((drhor(k)-droldr(k))**2+ (drhoi(k)-droldi(k))**2.gt.
     +        deps* (drhor(k)**2+drhoi(k)**2)) go to 20
   50 continue

   60 continue

      do k = 2, np1
        dt = drhor(k)*drhor(k-1) - drhoi(k)*drhoi(k-1)
        drhoi(k) = drhor(k)*drhoi(k-1) + drhoi(k)*drhor(k-1)
        drhor(k) = dt
      end do

      return
      end
      subroutine dlancz(n,ncap,dx,dw,dalpha,dbeta,ierr,dp0,dp1)

c*********************************************************************72
c
cc DLANCZ ???
c
c  Modified:
c
c    08 November 2012
c
c  Author:
c
c    Walter Gautschi
c
c  Reference:
c
c    Walter Gautschi,
c    Algorithm 793: GQRAT, Gauss Quadrature for Rational Functions,
c    ACM Transactions on Mathematical Software,
c    Volume 25, Number 2, pages 213-239, June 1999.
c
c  Parameters:
c
      implicit none

      integer ierr,n,ncap
c     ..
c     .. array arguments ..
      double precision dalpha(n),dbeta(n),dp0(ncap),dp1(ncap),dw(ncap),
     +                 dx(ncap)
c     ..
c     .. local scalars ..
      double precision dgam,dpi,drho,dsig,dt,dtk,dtmp,dtsig,dxlam
      integer i,k
c     ..
      if (n.le.0 .or. n.gt.ncap) then
          ierr = 1
          return

      else
          ierr = 0
      end if

      do 10 i = 1,ncap
          dp0(i) = dx(i)
          dp1(i) = 0.d0
   10 continue

      dp1(1) = dw(1)
      do 30 i = 1,ncap - 1
          dpi = dw(i+1)
          dgam = 1.d0
          dsig = 0.d0
          dt = 0.d0
          dxlam = dx(i+1)
          do 20 k = 1,i + 1
              drho = dp1(k) + dpi
              dtmp = dgam*drho
              dtsig = dsig
              if (drho.le.0.d0) then
                  dgam = 1.d0
                  dsig = 0.d0

              else
                  dgam = dp1(k)/drho
                  dsig = dpi/drho
              end if

              dtk = dsig* (dp0(k)-dxlam) - dgam*dt
              dp0(k) = dp0(k) - (dtk-dt)
              dt = dtk
              if (dsig.le.0.d0) then
                  dpi = dtsig*dp1(k)

              else
                  dpi = (dt**2)/dsig
              end if

              dtsig = dsig
              dp1(k) = dtmp
   20     continue
   30 continue

      do k = 1,n
        dalpha(k) = dp0(k)
        dbeta(k) = dp1(k)
      end do

      return
      end
      subroutine drecur(n,ipoly,dal,dbe,da,db,iderr)

c*********************************************************************72
c
cc DRECUR ???
c
c  Modified:
c
c    08 November 2012
c
c  Author:
c
c    Walter Gautschi
c
c  Reference:
c
c    Walter Gautschi,
c    Algorithm 793: GQRAT, Gauss Quadrature for Rational Functions,
c    ACM Transactions on Mathematical Software,
c    Volume 25, Number 2, pages 213-239, June 1999.
c
c  Parameters:
c
      implicit none

      double precision dal,dbe
      integer iderr,ipoly,n
c     ..
c     .. array arguments ..
      double precision da(n),db(n)
c     ..
c     .. local scalars ..
      double precision dal2,dalpbe,dbe2,dkm1,dlmach,dt
      integer k
c     ..
c     .. external functions ..
      double precision d1mach,dgamma,dlga
      external d1mach,dgamma,dlga
c     ..
c     .. intrinsic functions ..
      intrinsic atan,dble,exp,log,sqrt

      if (n.lt.1) then
          iderr = 3
          return
      end if

      dlmach = log(d1mach(2))
      iderr = 0

      do k = 1,n
          da(k) = 0.d0
      end do

      if (ipoly.eq.1) then
          db(1) = 2.d0
          if (n.eq.1) return
          do 20 k = 2,n
              dkm1 = dble(k-1)
              db(k) = 1.d0/ (4.d0-1.d0/ (dkm1*dkm1))
   20     continue
          return

      else if (ipoly.eq.2) then
          da(1) = .5d0
          db(1) = 1.d0
          if (n.eq.1) return
          do 30 k = 2,n
              da(k) = .5d0
              dkm1 = dble(k-1)
              db(k) = .25d0/ (4.d0-1.d0/ (dkm1*dkm1))
   30     continue
          return

      else if (ipoly.eq.3) then
          db(1) = 4.d0*atan(1.d0)
          if (n.eq.1) return
          db(2) = .5d0
          if (n.eq.2) return
          do 40 k = 3,n
              db(k) = .25d0
   40     continue
          return

      else if (ipoly.eq.4) then
          db(1) = 2.d0*atan(1.d0)
          if (n.eq.1) return
          do 50 k = 2,n
              db(k) = .25d0
   50     continue
          return

      else if (ipoly.eq.5) then
          db(1) = 4.d0*atan(1.d0)
          da(1) = .5d0
          if (n.eq.1) return
          do 60 k = 2,n
              db(k) = .25d0
   60     continue
          return

      else if (ipoly.eq.6) then
          if (dal.le.-1.d0 .or. dbe.le.-1.d0) then
              iderr = 1
              return

          else
              dalpbe = dal + dbe
              da(1) = (dbe-dal)/ (dalpbe+2.d0)
              dt = (dalpbe+1.d0)*log(2.d0) + dlga(dal+1.d0) +
     +             dlga(dbe+1.d0) - dlga(dalpbe+2.d0)
              if (dt.gt.dlmach) then
                  iderr = 2
                  db(1) = d1mach(2)

              else
                  db(1) = exp(dt)
              end if

              if (n.eq.1) return
              dal2 = dal*dal
              dbe2 = dbe*dbe
              da(2) = (dbe2-dal2)/ ((dalpbe+2.d0)* (dalpbe+4.d0))
              db(2) = 4.d0* (dal+1.d0)* (dbe+1.d0)/
     +                ((dalpbe+3.d0)* (dalpbe+2.d0)**2)
              if (n.eq.2) return
              do 70 k = 3,n
                  dkm1 = dble(k-1)
                  da(k) = .25d0* (dbe2-dal2)/
     +                    (dkm1*dkm1* (1.d0+.5d0*dalpbe/dkm1)*
     +                    (1.d0+.5d0* (dalpbe+2.d0)/dkm1))
                  db(k) = .25d0* (1.d0+dal/dkm1)* (1.d0+dbe/dkm1)*
     +                    (1.d0+dalpbe/dkm1)/ ((1.d0+.5d0* (dalpbe+
     +                    1.d0)/dkm1)* (1.d0+.5d0* (dalpbe-1.d0)/dkm1)*
     +                    (1.d0+.5d0*dalpbe/dkm1)**2)
   70         continue
              return

          end if

      else if (ipoly.eq.7) then
          if (dal.le.-1.d0) then
              iderr = 1
              return

          else
              da(1) = dal + 1.d0
              db(1) = dgamma(dal+1.d0,iderr)
              if (iderr.eq.2) db(1) = d1mach(2)
              if (n.eq.1) return
              do 80 k = 2,n
                  dkm1 = dble(k-1)
                  da(k) = 2.d0*dkm1 + dal + 1.d0
                  db(k) = dkm1* (dkm1+dal)
   80         continue
              return

          end if

      else if (ipoly.eq.8) then
          db(1) = sqrt(4.d0*atan(1.d0))
          if (n.eq.1) return
          do 90 k = 2,n
              db(k) = .5d0*dble(k-1)
   90     continue
          return

      else
          iderr = 4
      end if

      return
      end
      function dlga(dx)

c*********************************************************************72
c
cc DLGA evaluates the logarithm of the Gamma function.
c
c  Modified:
c
c    08 November 2012
c
c  Author:
c
c    Walter Gautschi
c
c  Reference:
c
c    Walter Gautschi,
c    Algorithm 793: GQRAT, Gauss Quadrature for Rational Functions,
c    ACM Transactions on Mathematical Software,
c    Volume 25, Number 2, pages 213-239, June 1999.
c
c  Parameters:
c
      implicit none

      double precision dlga
      double precision dx
      double precision dc,dp,ds,dt,dy
      real dprec,y,y0
      integer i
c     ..
c     .. local arrays ..
      double precision dbden(8),dbnum(8)
c     ..
c     .. external functions ..
      double precision d1mach
      external d1mach
c     ..
c     .. intrinsic functions ..
      intrinsic atan,exp,log,log10,real
C     ..
C     .. Data statements ..
c
c This routine evaluates the logarithm of the gamma function by a
c combination of recurrence and asymptotic approximation.
c
c The entries in the next data statement are the numerators and
c denominators, respectively, of the quantities B[16]/(16*15),
c B[14]/(14*13),..., B[2]/(2*1), where B[2n] are the Bernoulli
c numbers.
c
      data dbnum/-3.617d3,1.d0,-6.91d2,1.d0,-1.d0,1.d0,-1.d0,1.d0/,
     +     dbden/1.224d5,1.56d2,3.6036d5,1.188d3,1.68d3,1.26d3,3.6d2,
     +     1.2d1/
C     ..
c
c The quantity  dprec  in the next statement is the number of decimal
c digits carried in double-precision floating-point arithmetic.
c
      dprec = -log10(real(d1mach(3)))
      dc = .5d0*log(8.d0*atan(1.d0))
      dp = 1.d0
      dy = dx
      y = real(dy)
c
c The quantity  y0  below is the threshold value beyond which asymptotic
c evaluation gives sufficient accuracy; see Eq. 6.1.42 in M. Abramowitz
c and I.A. Stegun,Handbook of Mathematical Functions''. The constants
c are .12118868... = ln(10)/19 and .05390522... = ln(|B[20]|/190)/19.
c
      y0 = exp(.121189*dprec+.053905)
   10 if (y.gt.y0) go to 20
      dp = dy*dp
      dy = dy + 1.d0
      y = real(dy)
      go to 10

   20 dt = 1.d0/ (dy*dy)
c
c the right-hand side of the next assignment statement is b[18]/(18*17).
c
      ds = 4.3867d4/2.44188d5
      do i = 1,8
        ds = dt*ds + dbnum(i)/dbden(i)
      end do
      dlga = (dy-.5d0)*log(dy) - dy + dc + ds/dy - log(dp)

      return
      end
      function dgamma(dx,iderr)

c*********************************************************************72
c
cc DGAMMA evaluates the gamma function.
c
c  Modified:
c
c    08 November 2012
c
c  Author:
c
c    Walter Gautschi
c
c  Reference:
c
c    Walter Gautschi,
c    Algorithm 793: GQRAT, Gauss Quadrature for Rational Functions,
c    ACM Transactions on Mathematical Software,
c    Volume 25, Number 2, pages 213-239, June 1999.
c
c  Parameters:
c
c This evaluates the gamma function for real positive  dx, using the
c function subroutine  dlga.
c
      implicit none

      double precision dgamma
      double precision dx
      integer iderr
c     ..
c     .. local scalars ..
      double precision dlmach,dt
c     ..
c     .. external functions ..
      double precision d1mach,dlga
      external d1mach,dlga
c     ..
c     .. intrinsic functions ..
      intrinsic exp,log
c     ..
      dlmach = log(d1mach(2))
      iderr = 0
      dt = dlga(dx)
      if (dt.ge.dlmach) then
          iderr = 2
          dgamma = d1mach(2)
          return

      else
          dgamma = exp(dt)
          return

      end if

      end
      subroutine dsti(n,ncap,dx,dw,dalpha,dbeta,ierr,dp0,dp1,dp2)

c*********************************************************************72
c
cc DSTI ???
c
c  Modified:
c
c    08 November 2012
c
c  Author:
c
c    Walter Gautschi
c
c  Reference:
c
c    Walter Gautschi,
c    Algorithm 793: GQRAT, Gauss Quadrature for Rational Functions,
c    ACM Transactions on Mathematical Software,
c    Volume 25, Number 2, pages 213-239, June 1999.
c
c  Parameters:
c
      implicit none

      integer ierr,n,ncap
c     ..
c     .. array arguments ..
      double precision dalpha(n),dbeta(n),dp0(ncap),dp1(ncap),dp2(ncap),
     +                 dw(ncap),dx(ncap)
c     ..
c     .. local scalars ..
      double precision dhuge,dsum0,dsum1,dsum2,dt,dtiny
      integer k,m,nm1
c     ..
c     .. external functions ..
      double precision d1mach
      external d1mach
c     ..
c     .. intrinsic functions ..
      intrinsic abs
c     ..
      dtiny = 10.d0*d1mach(1)
      dhuge = .1d0*d1mach(2)
      ierr = 0
      if (n.le.0 .or. n.gt.ncap) then
          ierr = 1
          return

      end if

      nm1 = n - 1
      dsum0 = 0.d0
      dsum1 = 0.d0
      do 10 m = 1,ncap
          dsum0 = dsum0 + dw(m)
          dsum1 = dsum1 + dw(m)*dx(m)
   10 continue
      dalpha(1) = dsum1/dsum0
      dbeta(1) = dsum0
      if (n.eq.1) return
      do m = 1,ncap
          dp1(m) = 0.d0
          dp2(m) = 1.d0
      end do

      do 40 k = 1,nm1
          dsum1 = 0.d0
          dsum2 = 0.d0
          do 30 m = 1,ncap
              if (dw(m).eq.0.d0) go to 30
              dp0(m) = dp1(m)
              dp1(m) = dp2(m)
              dp2(m) = (dx(m)-dalpha(k))*dp1(m) - dbeta(k)*dp0(m)
              if (abs(dp2(m)).gt.dhuge .or. abs(dsum2).gt.dhuge) then
                  ierr = k
                  return

              end if

              dt = dw(m)*dp2(m)*dp2(m)
              dsum1 = dsum1 + dt
              dsum2 = dsum2 + dt*dx(m)
   30     continue
          if (abs(dsum1).lt.dtiny) then
              ierr = -k
              return

          end if

          dalpha(k+1) = dsum2/dsum1
          dbeta(k+1) = dsum1/dsum0
          dsum0 = dsum1
   40 continue

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
c      Sign * (X(S-1)*A^(S-1) + ... + X(1)*A + X(0))
c
c    where 0 <= X(1:S-1) < A.
c
c      I1MACH(7) = A, the base.
c      I1MACH(8) = S, the number of base A digits.
c      I1MACH(9) = A^S-1, the largest integer.
c
c    Floating point numbers
c
c    Assume floating point numbers are represented in the T digit 
c    base B form:
c
c      Sign * (B^E) * ((X(1)/B) + ... + (X(T)/B^T) )
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
      FUNCTION NU0HER(N,Z,EPS)

c*********************************************************************72
c
cc NU0HER ???
c
c  Modified:
c
c    08 November 2012
c
c  Author:
c
c    Walter Gautschi
c
c  Reference:
c
c    Walter Gautschi,
c    Algorithm 793: GQRAT, Gauss Quadrature for Rational Functions,
c    ACM Transactions on Mathematical Software,
c    Volume 25, Number 2, pages 213-239, June 1999.
c
c  Parameters:
c
c This is an auxiliary function routine providing a starting backward
c recurrence index for the Hermite measure that can be used in place
c of  nu0  in the routines  knum  and  dknum.
c
      implicit none

      integer nu0her
      complex z
      real eps
      integer n
      intrinsic abs,aimag,log,real,sqrt

      nu0her = 2.* (sqrt(.5*real(n+1))+.25*log(1./eps)/abs(aimag(z)))**2

      return
      end
      function nu0jac(n,z,eps)

c*********************************************************************72
c
cc NU0JAC ???
c
c  Modified:
c
c    08 November 2012
c
c  Author:
c
c    Walter Gautschi
c
c  Reference:
c
c    Walter Gautschi,
c    Algorithm 793: GQRAT, Gauss Quadrature for Rational Functions,
c    ACM Transactions on Mathematical Software,
c    Volume 25, Number 2, pages 213-239, June 1999.
c
c  Parameters:
c
c This is an auxiliary function routine providing a starting backward
c recurrence index for the Jacobi measure that can be used in place
c of  nu0  in the routines  knum  and  dknum.
c
      implicit none

      integer nu0jac
      complex z
      real eps
      integer n
c     ..
c     .. local scalars ..
      real angle,pi,r,x,x2,y,y2
c     ..
c     .. intrinsic functions ..
      intrinsic abs,aimag,atan,cos,log,real,sin,sqrt
c     ..
      pi = 4.*atan(1.)
      x = real(z)
      y = abs(aimag(z))
      if (x.lt.1.) then
          if (x.lt.-1.) angle = .5* (2.*pi+atan(y/ (x-1.))+
     +                          atan(y/ (x+1.)))
          if (x.eq.-1.) angle = .5* (1.5*pi-atan(.5*y))
          if (x.gt.-1.) angle = .5* (pi+atan(y/ (x-1.))+atan(y/ (x+1.)))

      else
          if (x.eq.1.) angle = .5* (.5*pi+atan(.5*y))
          if (x.gt.1.) angle = .5* (atan(y/ (x-1.))+atan(y/ (x+1.)))
      end if

      x2 = x*x
      y2 = y*y
      r = ((x2-y2-1.)**2+4.*x2*y2)**.25
      r = sqrt((x+r*cos(angle))**2+ (y+r*sin(angle))**2)
      nu0jac = real(n+1) + .5*log(1./eps)/log(r)

      return
      end
      function nu0lag(n,z,al,eps)

c*********************************************************************72
c
cc NU0LAG ???
c
c  Modified:
c
c    08 November 2012
c
c  Author:
c
c    Walter Gautschi
c
c  Reference:
c
c    Walter Gautschi,
c    Algorithm 793: GQRAT, Gauss Quadrature for Rational Functions,
c    ACM Transactions on Mathematical Software,
c    Volume 25, Number 2, pages 213-239, June 1999.
c
c  Parameters:
c
c This is an auxiliary function routine providing a starting backward
c recurrence index for the Laguerre measure that can be used in place
c of  nu0  in the routines  knum  and  dknum.
c
      implicit none

      integer nu0lag
      complex z
      real al,eps
      integer n
c     ..
c     .. local scalars ..
      real phi,pi,x,y
c     ..
c     .. intrinsic functions ..
      intrinsic aimag,atan,cos,log,real,sqrt
c     ..
      pi = 4.*atan(1.)
      x = real(z)
      y = aimag(z)
      phi = .5*pi
      if (y.lt.0.) phi = 1.5*pi
      if (x.eq.0.) go to 10
      phi = atan(y/x)
      if (y.gt.0. .and. x.gt.0.) go to 10
      phi = phi + pi
      if (x.lt.0.) go to 10
      phi = phi + pi
   10 nu0lag = (sqrt(real(n+1)+.5* (al+1.))+
     +         log(1./eps)/ (4.* (x*x+y*y)**.25*cos(.5* (phi-pi))))**2 -
     +         .5* (al+1.)
      return
      end
      subroutine timestamp ( )

c*********************************************************************72
c
cc TIMESTAMP prints out the current YMDHMS date as a timestamp.
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
