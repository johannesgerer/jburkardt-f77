      program main

c*********************************************************************72
c
cc MAIN is the main program of TOMS612_PRB.
c
c  Discussion:
c
c    This program demonstrates the use of TRIEX for the integration
c    of a function over a triangle.  Four integrand functions are tried.
c
c    The form in which the test integrals are supplied to the
c    integrator corresponds to
c      I(PA,PB) I(PC*X+PD,PG*X+PH) F(X,Y) DXDY
c    with
c      F(X,Y) = G(U(X),V(X,Y))*V1(X)/D
c    where
c      D = PB-PA
c      V1(X) = (1-(X-PA)/D)/((PG-PC)*X+PH-PD)
c      U(X) = (X-PA)/D
c      V(X,Y) = V1(X)*(Y-PC*X-PD)
c
c    Hence this is the transformed version of the integral over
c    the unit triangle given by
c      I(0,1) I(0,1-X) G(X,Y) DXDY
c
c    The following list provides, for each test integral, the
c    parameters PA, PB, PC, PD, PG, PH, and the vertices of the
c    triangle over which is thus integrated: VER(1,K) = ABSCISSAE
c    and VER(2,K) = ordinates of the vertices indexed K.
c
c    Integral number 1
c
c    G(X,Y) = X**(-0.2)
c    PA=1  PB=5  PC=0  PD=0  PG=0.25  PH=-0.25
c    VER(1,1)=1  VER(1,2)=5  VER(1,3)=5
c    VER(2,1)=0  VER(2,2)=0  VER(2,3)=1
c
c    Integral number 2
c
c    G(X,Y) = (X+Y)**(-0.2)
c    PA=0  PB=1  PC=0  PD=0  PG=-1  PH=1
c    VER(1,1)=0  VER(1,2)=1  VER(1,3)=0
c    VER(2,1)=0  VER(2,2)=0  VER(2,3)=1
c
c    Integral number 3
c
c    G(X,Y) = (1-X-Y)**(-0.2)
c    PA=-1  PB=3  PC=0.25  PD=-2.75  PG=-1  PH=1
c    VER(1,1)=-1  VER(1,2)=3  VER(1,3)=-1
c    VER(2,1)=-3  VER(2,2)=-2  VER(2,3)=2
c
c    Integral number 4
c
c    G(X,Y) = (X*Y)**(-0.2)
c    PA=0  PB=-7  PC=0  PD=0  PG=-3/7  PH=-3
c    VER(1,1)=0  VER(1,2)=-7  VER(1,3)=0
c    VER(2,1)=0  VER(2,2)=0  VER(2,3)=-3
c
c  Modified:
c
c    07 December 2006
c
c  Author:
c
c    Elise deDoncker, Ian Robinson,
c
c  Reference:
c
c    Elise deDoncker, Ian Robinson,
c    Algorithm 612:
c    Integration over a Triangle Using Nonlinear Extrapolation,
c    ACM Transactions on Mathematical Software,
c    Volume 10, Number 1, March 1984, pages 17-22.
c
      implicit none

      double precision abserr
      double precision c(4)
      double precision d(4)
      double precision epsabs
      double precision epsrel
      double precision er
      double precision es
      double precision exa
      double precision exact(4)
      double precision f
      external f
      double precision g(4)
      double precision h(4)
      integer ier
      integer indacc
      integer inum
      integer kount
      integer limev
      integer n2
      integer n2m1
      integer nout
      integer num
      double precision p(8)
      double precision pa
      double precision pb
      double precision pc
      double precision pd
      double precision pg
      double precision ph
      double precision q(8)
      double precision r(8)
      double precision result
      double precision ver(2,3)

      common / cnt / kount, num
      common / par / pa, pb, pc, pd, pg, ph
c
c  Exact values
c
      data exact(1), exact(2), exact(3), exact(4) /
     &  0.6944444444444444D+00,
     &  0.5555555555555556D+00,
     &  0.6944444444444444D+00,
     &  0.9481026454955768D+00 /
c
c   (P(2*NUM-1),P(2*NUM)) is 1st vertex for integral NUM=1,2,3,4
c
      data p(1), p(2), p(3), p(4), p(5), p(6), p(7), p(8)
     & /1.0D+00,0.0D+00,0.0D+00,0.0D+00,-1.0D+00,-3.0D+00,0.0D+00,
     & 0.0D+00/
c
c   (Q(2*NUM-1),Q(2*NUM)) is 2nd vertex for integral NUM=1,2,3,4
c
      data q(1), q(2), q(3), q(4), q(5), q(6), q(7), q(8)
     & /5.0D+00,0.0D+00,1.0D+00,0.0D+00,3.0D+00,-2.0D+00,-7.0D+00,
     & 0.0D+00/
c
c   (R(2*NUM-1),R(2*NUM)) is 3rd vertex for integral NUM=1,2,3,4
c
      data r(1), r(2), r(3), r(4), r(5), r(6), r(7), r(8)
     & /5.0D+00,1.0D+00,0.0D+00,1.0D+00,-1.0D+00,2.0D+00,0.0D+00,
     & -3.0D+00/
c
c   Parameters
c
      data C(1), C(2), C(3), C(4) /0.0D+00,0.0D+00,0.25D+00,0.0D+00/
      data D(1), D(2), D(3), D(4) /0.0D+00,0.0D+00,-2.75D+00,0.0D+00/
      data G(1), G(2), G(3), G(4) /0.25D+00,-1.0D+00,-1.0D+00,
     & -0.428571428571428571D+00/
      data h(1), h(2), h(3), h(4) /-0.25D+00,1.0D+00,1.0D+00,-3.0D+00/

      call timestamp ( )
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TOMS612_PRB'
      write ( *, '(a)' ) '  FORTRAN77 version'
      write ( *, '(a)' ) '  Test the TOMS612 library.'

      WRITE ( *,99998)
99998 FORMAT (6X, 8HRELATIVE, 9X, 8HCOMPUTED, 10X, 8HRELATIVE, 2X,
     * 6HESTIM., 3X, 6HNUMBER/1X, 3HNR., 2X, 8HACCURACY, 9X, 8HINTEGRAL,
     * 10X, 5HERROR, 5X, 8HRELATIVE, 1X, 6HOF  F , 4H IER/6X, 7HREQUEST,
     * 2HED, 36X, 5HERROR, 4X, 5HEVAL.)

      epsabs = 0.0d+00
      limev = 50000

      do inum = 1, 4

        exa = exact(inum)
        epsrel = 1.0d-02

        do indacc = 1, 3

          kount = 0
          num = inum
          n2 = num + num
          n2m1 = n2 - 1
          ver(1,1) = p(n2m1)
          ver(2,1) = p(n2)
          ver(1,2) = q(n2m1)
          ver(2,2) = q(n2)
          ver(1,3) = r(n2m1)
          ver(2,3) = r(n2)
          pa = p(n2m1)
          pb = q(n2m1)
          pc = c(num)
          pd = d(num)
          pg = g(num)
          ph = h(num)

          call triex ( f, ver, epsabs, epsrel, limev, result,
     &      abserr, ier )

          er = dabs ( ( result - exa ) / exa )
          es = abserr / dabs ( result )

          if ( indacc .le. 1 ) then
            write ( *, '(a)' ) ' '
            write ( *,99997) num, epsrel, result, er, es, kount, ier
99997       format (1x, i2, d11.2, d25.16, 2d10.2, i7, i4)
          else
            write ( *,99996) epsrel, result, er, es, kount, ier
99996       format (3X, D11.2, D25.16, 2D10.2, I7, I4)
          end if

          epsrel = epsrel * 1.0D-02

        end do

        write ( *, '(a,d25.16)' ) '      Exact = ', exa

      end do
c
c  Terminate.
c
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TOMS612_PRB'
      write ( *, '(a)' ) '  Normal end of execution.'

      write ( *, '(a)' ) ' '
      call timestamp ( )

      stop
      end
      function f ( x, y )

c*********************************************************************72
c
cc F evaluates the test integrand functions.
c
c  Modified:
c
c    07 December 2006
c
c  Author:
c
c    Elise deDoncker, Ian Robinson,
c
c  Reference:
c
c    Elise deDoncker, Ian Robinson,
c    Algorithm 612:
c    Integration over a Triangle Using Nonlinear Extrapolation,
c    ACM Transactions on Mathematical Software,
c    Volume 10, Number 1, March 1984, pages 17-22.
c
c  Parameters:
c
c    Input, double precision X, Y, the evaluation point.
c
c    Output, double precision F, the value of the integrand at (X,Y).
c
      implicit none

      double precision a
      double precision d
      double precision d1mach
      double precision drelpr
      double precision f
      integer kount
      integer num
      double precision p
      double precision pa
      double precision pb
      double precision pc
      double precision pd
      double precision pg
      double precision ph
      double precision u
      double precision v
      double precision v1
      double precision x
      double precision y

      common / cnt / kount, num
      common / par / pa, pb, pc, pd, pg, ph
c
c   DRELPR  is the largest relative spacing
c
c   CALL FUNCTION SUBPROGRAM D1MACH FROM  P O R T  MATHEMATICAL
c   SUBROUTINE LIBRARY (BELL LABORATORIES)
c
      drelpr = d1mach ( 4 )
      p = -0.2D+00
      kount = kount + 1
      f = 0.0D+00
      d = pb - pa
      u = ( x - pa ) / d
      v1 = ( 1.0D+00 - u ) / ( ( pg - pc ) * x + ph - pd )
      v = v1 * ( y - pc * x - pd )
c
c  Singularity at U=0, I.E. X=PA
c  (Vertex of 1st sample triangle)
c
      if ( num .eq. 1 ) then

        if ( 0.0D+00 .lt. u ) then
          f = u**p
        else
          f = 0.0D+00
        end if
c
c  Singularity at U=V=0, I.E. X=PA and Y=PC*X+PD
c  (Vertex of 2nd sample triangle)
c
      else if ( num .eq. 2 ) then

        a = u + v

        if ( 0.0D+00 .lt. a ) then
          f = a**p
        else
          f = 0.0D+00
        end if
c
c  Singularity along U+V=1
c  (Side of 3rd sample triangle)
c
      else if ( num .eq. 3 ) then

        a = 1.0D+00 - u - v

        if ( drelpr .le. a ) then
          f = a**p
        else
          f = 0.0D+00
        end if
c
c  Singularity along U=0 and V=0
c  (Two sides of 4th sample triangle)
c
      else if ( num .eq. 4 ) then

        a = u * v

        if ( 0.0D+00 .lt. a ) then
          f = a**p
        else
          f = 0.0D+00
        end if

      end if

      f = f * v1 / d

      return
      end
