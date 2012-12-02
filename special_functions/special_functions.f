      subroutine airya ( x, ai, bi, ad, bd )

c*********************************************************************72
c
cc AIRYA computes Airy functions and their derivatives.
c
c  Licensing:
c
c    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
c    they give permission to incorporate this routine into a user program 
c    provided that the copyright is acknowledged.
c
c  Modified:
c
c    14 March 2012
c
c  Author:
c
c    Shanjie Zhang, Jianming Jin
c
c  Reference:
c
c    Shanjie Zhang, Jianming Jin,
c    Computation of Special Functions,
c    Wiley, 1996,
c    ISBN: 0-471-11963-6,
c    LC: QA351.C45.
c       
c  Parameters:
c
c    Input, double precision X, the argument of the Airy function.
c
c    Output, double precision AI, BI, AD, BD, the values of Ai(x), Bi(x),
c    Ai'(x), Bi'(x).
c
      implicit none

      double precision ad
      double precision ai
      double precision bd
      double precision bi
      double precision c1
      double precision c2
      double precision pir
      double precision sr3
      double precision vi1
      double precision vi2
      double precision vj1
      double precision vj2
      double precision vk1
      double precision vk2
      double precision vy1
      double precision vy2
      double precision x
      double precision xa
      double precision xq
      double precision z

      xa = abs ( x )
      pir = 0.318309886183891D+00
      c1 = 0.355028053887817D+00
      c2 = 0.258819403792807D+00
      sr3 = 1.732050807568877D+00
      z = xa ** 1.5D+00 / 1.5D+00
      xq = sqrt ( xa )

      call ajyik ( z, vj1, vj2, vy1, vy2, vi1, vi2, vk1, vk2 )

      if ( x .eq. 0.0D+00 ) then
        ai = c1
        bi = sr3 * c1
        ad = - c2
        bd = sr3 * c2
      else if ( 0.0D+00 .lt. x ) then
        ai = pir * xq / sr3 * vk1
        bi = xq * ( pir * vk1 + 2.0D+00 / sr3 * vi1 )
        ad = - xa / sr3 * pir * vk2
        bd = xa * ( pir * vk2 + 2.0D+00 / sr3 * vi2 )
      else
        ai = 0.5D+00 * xq * ( vj1 - vy1 / sr3 )
        bi = - 0.5D+00 * xq * ( vj1 / sr3 + vy1 )
        ad = 0.5D+00 * xa * ( vj2 + vy2 / sr3 )
        bd = 0.5D+00 * xa * ( vj2 / sr3 - vy2 )
      end if

      return
      end
      subroutine airyb ( x, ai, bi, ad, bd )

c*********************************************************************72
c
cc AIRYB computes Airy functions and their derivatives.
c
c  Licensing:
c
c    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
c    they give permission to incorporate this routine into a user program 
c    provided that the copyright is acknowledged.
c
c  Modified:
c
c    02 June 2012
c
c  Author:
c
c    Shanjie Zhang, Jianming Jin
c
c  Reference:
c
c    Shanjie Zhang, Jianming Jin,
c    Computation of Special Functions,
c    Wiley, 1996,
c    ISBN: 0-471-11963-6,
c    LC: QA351.C45.
c
c  Parameters:
c
c    Input, double precision X, argument of Airy function.
c
c    Output, double precision AI, Ai(x).
c
c    Output, double precision BI, Bi(x).
c
c    Output, double precision AD, Ai'(x).
c
c    Output, double precision BD, Bi'(x).
c
      implicit none

      double precision ad
      double precision ai
      double precision bd
      double precision bi
      double precision c1
      double precision c2
      double precision ck(41)
      double precision df
      double precision dg
      double precision dk(41)
      double precision eps
      double precision fx
      double precision gx
      integer k
      integer km
      double precision pi
      double precision r
      double precision rp
      double precision sad
      double precision sai
      double precision sbd
      double precision sbi
      double precision sda
      double precision sdb
      double precision sr3
      double precision ssa
      double precision ssb
      double precision x
      double precision xa
      double precision xar
      double precision xcs
      double precision xe
      double precision xf
      double precision xm
      double precision xp1
      double precision xq
      double precision xr1
      double precision xr2
      double precision xss

      eps = 1.0D-15
      pi = 3.141592653589793D+00
      c1 = 0.355028053887817D+00
      c2 = 0.258819403792807D+00
      sr3 = 1.732050807568877D+00
      xa = abs ( x )
      xq = sqrt ( xa )

      if ( x .le. 0.0D+00 ) then
        xm = 8.0D+00
      else
        xm = 5.0D+00
      end if

      if ( x .eq. 0.0D+00 ) then
        ai = c1
        bi = sr3 * c1
        ad = -c2
        bd = sr3 * c2
        return
      end if

      if ( xa .le. xm ) then

        fx = 1.0D+00
        r = 1.0D+00
        do k = 1, 40
          r = r * x / ( 3.0D+00 * k ) * x 
     &      / ( 3.0D+00 * k - 1.0D+00 ) * x
          fx = fx + r
          if ( abs ( r ) .lt. abs ( fx ) * eps ) then
            go to 10
          end if
        end do

10      continue

        gx = x
        r = x
        do k = 1, 40
          r = r * x / ( 3.0D+00 * k ) * x 
     &      / ( 3.0D+00 * k + 1.0D+00 ) * x
          gx = gx + r
          if ( abs ( r ) .lt. abs ( gx ) * eps ) then
            go to 20
          end if
        end do

20      continue

        ai = c1 * fx - c2 * gx
        bi = sr3 * ( c1 * fx + c2 * gx )
        df = 0.5D+00 * x * x
        r = df
        do k = 1, 40
          r = r * x / ( 3.0D+00 * k ) * x 
     &      / ( 3.0D+00 * k + 2.0D+00 ) * x
          df = df + r
          if ( abs ( r ) .lt. abs ( df ) * eps ) then
            go to 30
          end if
        end do

30      continue

        dg = 1.0D+00
        r = 1.0D+00
        do k = 1,40
          r = r * x / ( 3.0D+00 * k ) * x 
     &      / ( 3.0D+00 * k - 2.0D+00 ) * x
          dg = dg + r
          if ( abs ( r ) .lt. abs ( dg ) * eps ) then
            go to 40
          end if
        end do

40      continue

        ad = c1 * df - c2 * dg
        bd = sr3 * ( c1 * df + c2 * dg )

      else

        xe = xa * xq / 1.5D+00
        xr1 = 1.0D+00 / xe
        xar = 1.0D+00 / xq
        xf = sqrt ( xar )
        rp = 0.5641895835477563D+00
        r = 1.0D+00
        do k = 1, 40
          r = r * ( 6.0D+00 * k - 1.0D+00 ) / 216.0D+00 
     &      * ( 6.0D+00 * k - 3.0D+00 ) / k
     &      * ( 6.0D+00 * k - 5.0D+00 ) / ( 2.0D+00 * k - 1.0D+00 )
          ck(k) = r
          dk(k) = - ( 6.0D+00 * k + 1.0D+00 ) 
     &      / ( 6.0D+00 * k - 1.0D+00 ) * ck(k)
        end do

        km = int ( 24.5D+00 - xa )

        if ( xa .lt. 6.0D+00 ) then
          km = 14
        end if

        if ( 15.0D+00 .lt. xa ) then
          km = 10
        end if

        if ( 0.0D+00 .lt. x ) then
          sai = 1.0D+00
          sad = 1.0D+00
          r = 1.0D+00
          do k = 1, km
            r = -r * xr1
            sai = sai + ck(k) * r
            sad = sad + dk(k) * r
          end do
          sbi = 1.0D+00
          sbd = 1.0D+00
          r = 1.0D+00
          do k = 1, km
            r = r * xr1
            sbi = sbi + ck(k) * r
            sbd = sbd + dk(k) * r
          end do
          xp1 = exp ( - xe )
          ai = 0.5D+00 * rp * xf * xp1 * sai
          bi = rp * xf / xp1 * sbi
          ad = -0.5D+00 * rp / xf * xp1 * sad
          bd = rp / xf / xp1 * sbd
        else
          xcs = cos ( xe + pi / 4.0D+00 )
          xss = sin ( xe + pi / 4.0D+00 )
          ssa = 1.0D+00
          sda = 1.0D+00
          r = 1.0D+00
          xr2 = 1.0D+00 / ( xe * xe )
          do k = 1, km
            r = - r * xr2
            ssa = ssa + ck(2*k) * r
            sda = sda + dk(2*k) * r
          end do
          ssb = ck(1) * xr1
          sdb = dk(1) * xr1
          r = xr1
          do k = 1, km
            r = -r * xr2
            ssb = ssb + ck(2*k+1) * r
            sdb = sdb + dk(2*k+1) * r
          end do

          ai =   rp * xf * ( xss * ssa - xcs * ssb )
          bi =   rp * xf * ( xcs * ssa + xss * ssb )
          ad = - rp / xf * ( xcs * sda + xss * sdb )
          bd =   rp / xf * ( xss * sda - xcs * sdb )

        end if

      end if

      return
      end
      subroutine airyzo ( nt, kf, xa, xb, xc, xd )

c*********************************************************************72
c
cc AIRYZO computes the first NT zeros of Ai(x) and Ai'(x).
c
c   Discussion:
c
c    Compute the first NT zeros of Airy functions Ai(x) and Ai'(x), 
c    a and a', and the associated values of Ai(a') and Ai'(a); and 
c    the first NT zeros of Airy functions Bi(x) and Bi'(x), b and
c    b', and the associated values of Bi(b') and Bi'(b).
c
c  Licensing:
c
c    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
c    they give permission to incorporate this routine into a user program 
c    provided that the copyright is acknowledged.
c
c  Modified:
c
c    14 March 2012
c
c  Author:
c
c    Shanjie Zhang, Jianming Jin
c
c  Reference:
c
c    Shanjie Zhang, Jianming Jin,
c    Computation of Special Functions,
c    Wiley, 1996,
c    ISBN: 0-471-11963-6,
c    LC: QA351.C45.
c
c  Parameters:
c
c    Input, integer NT, the number of zeros.
c
c    Input, integer KF, the function code.
c    1 for Ai(x) and Ai'(x);
c    2 for Bi(x) and Bi'(x).
c
c    Output, double precision XA(m), a, the m-th zero of Ai(x) or
c    b, the m-th zero of Bi(x).
c
c    Output, double precision XB(m), a', the m-th zero of Ai'(x) or
c    b', the m-th zero of Bi'(x).
c
c    Output, double precision XC(m), Ai(a') or Bi(b').
c
c    Output, double precision XD(m), Ai'(a) or Bi'(b)
c
      implicit none

      integer nt

      double precision ad
      double precision ai
      double precision bd
      double precision bi
      integer i
      integer kf
      double precision pi
      double precision rt
      double precision rt0
      double precision u
      double precision u1
      double precision x
      double precision xa(nt)
      double precision xb(nt)
      double precision xc(nt)
      double precision xd(nt)

      pi = 3.141592653589793D+00

      do i = 1, nt

        if ( kf .eq. 1 ) then

          u = 3.0D+00 * pi * ( 4.0D+00 * i - 1 ) / 8.0D+00
          u1 = 1.0D+00 / ( u * u )
          rt0 = - ( u * u ) ** ( 1.0D+00 / 3.0D+00 ) * ((((
     &      - 15.5902D+00 * u1
     &      + 0.929844D+00 ) * u1
     &      - 0.138889D+00 ) * u1
     &      + 0.10416667D+00 ) * u1
     &      + 1.0D+00 )

        else if ( kf .eq. 2 ) then

          if ( i .eq. 1 ) then
            rt0 = -1.17371D+00
          else
            u = 3.0D+00 * pi * ( 4.0D+00 * i - 3.0D+00 ) / 8.0D+00
            u1 = 1.0D+00 / ( u * u )
            rt0 = - ( u * u ) ** ( 1.0D+00 / 3.0D+00 ) * ((((
     &        - 15.5902D+00 * u1
     &        + 0.929844D+00 ) * u1
     &        - 0.138889D+00 ) * u1
     &        + 0.10416667D+00 ) * u1
     &        + 1.0D+00 )
          end if

        end if

10      continue

        x = rt0
        call airyb ( x, ai, bi, ad, bd )

        if ( kf .eq. 1 ) then
          rt = rt0 - ai / ad
        else
          rt = rt0 - bi / bd
        end if

        if ( 1.0D-09 .lt. abs ( ( rt - rt0 ) / rt ) ) then
          rt0 = rt
          go to 10
        else
          xa(i) = rt
          if ( kf .eq. 1 ) then
            xd(i) = ad
          else
            xd(i) = bd
          end if
        end if

      end do

      do i = 1, nt

        if ( kf .eq. 1 ) then
          if ( i .eq. 1 ) then
            rt0 = -1.01879D+00
          else
            u = 3.0D+00 * pi * ( 4.0D+00 * i - 3.0D+00 ) / 8.0D+00
            u1 = 1.0D+00 / ( u * u )
            rt0 = - ( u * u ) ** ( 1.0D+00 / 3.0D+00 ) * ((((
     &          15.0168D+00 * u1
     &        - 0.873954D+00 ) * u1
     &        + 0.121528D+00 ) * u1
     &        - 0.145833D+00 ) * u1
     &        + 1.0D+00 )
          end if
        else if ( kf .eq. 2 ) then
          if ( i .eq. 1 ) then
            rt0 = -2.29444D+00
          else
            u = 3.0D+00 * pi * ( 4.0D+00 * i - 1.0D+00 ) / 8.0D+00
            u1 = 1.0D+00 / ( u * u )
            rt0 = - ( u * u ) ** ( 1.0D+00 / 3.0D+00 ) * ((((
     &          15.0168D+00 * u1
     &        - 0.873954D+00 ) * u1
     &        + 0.121528D+00 ) * u1
     &        - 0.145833D+00 ) * u1 + 1.0D+00 )
          end if
        end if

20      continue

        x = rt0
        call airyb ( x, ai, bi, ad, bd )

        if ( kf .eq. 1 ) then
          rt = rt0 - ad / ( ai * x )
        else
          rt = rt0 - bd / ( bi * x )
        end if

        if ( 1.0D-09 .lt. abs ( ( rt - rt0 ) / rt ) ) then
          rt0 = rt
          go to 20
        else
          xb(i) = rt
          if ( kf .eq. 1 ) then
            xc(i) = ai
          else
            xc(i) = bi
          end if
        end if

      end do

      return
      end
      subroutine ajyik ( x, vj1, vj2, vy1, vy2, vi1, vi2, vk1, vk2 )

c*********************************************************************72
c
cc AJYIK computes Bessel functions Jv(x), Yv(x), Iv(x), Kv(x).
c
c  Discussion: 
c
c    Compute Bessel functions Jv(x) and Yv(x), and modified Bessel functions 
c    Iv(x) and Kv(x), and their derivatives with v = 1/3, 2/3.
c
c  Licensing:
c
c    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
c    they give permission to incorporate this routine into a user program 
c    provided that the copyright is acknowledged.
c
c  Modified:
c
c    31 July 2012
c
c  Author:
c
c    Shanjie Zhang, Jianming Jin
c
c  Reference:
c
c    Shanjie Zhang, Jianming Jin,
c    Computation of Special Functions,
c    Wiley, 1996,
c    ISBN: 0-471-11963-6,
c    LC: QA351.C45.
c
c  Parameters:
c
c    Input, double precision X, the argument.  X should not be zero.
c
c    Output, double precision VJ1, VJ2, VY1, VY2, VI1, VI2, VK1, VK2,
c    the values of J1/3(x), J2/3(x), Y1/3(x), Y2/3(x), I1/3(x), I2/3(x),
c    K1/3(x), K2/3(x).
c
      implicit none

      double precision a0
      double precision b0
      double precision c0
      double precision ck
      double precision gn
      double precision gn1
      double precision gn2
      double precision gp1
      double precision gp2
      integer k
      integer k0
      integer l
      double precision pi
      double precision pv1
      double precision pv2
      double precision px
      double precision qx
      double precision r
      double precision rp
      double precision rp2
      double precision rq
      double precision sk
      double precision sum
      double precision uj1
      double precision uj2
      double precision uu0
      double precision vi1
      double precision vi2
      double precision vil
      double precision vj1
      double precision vj2
      double precision vjl
      double precision vk1
      double precision vk2
      double precision vl
      double precision vsl
      double precision vv
      double precision vv0
      double precision vy1
      double precision vy2
      double precision x
      double precision x2
      double precision xk

      if ( x .eq. 0.0D+00 ) then
        vj1 = 0.0D+00
        vj2 = 0.0D+00
        vy1 = -1.0D+300
        vy2 = 1.0D+300
        vi1 = 0.0D+00
        vi2 = 0.0D+00
        vk1 = -1.0D+300
        vk2 = -1.0D+300
        return
      end if

      pi = 3.141592653589793D+00
      rp2 = 0.63661977236758D+00
      gp1 = 0.892979511569249D+00
      gp2 = 0.902745292950934D+00
      gn1 = 1.3541179394264D+00
      gn2 = 2.678938534707747D+00
      vv0 = 0.444444444444444D+00
      uu0 = 1.1547005383793D+00
      x2 = x * x

      if ( x .lt. 35.0D+00 ) then
        k0 = 12
      else if ( x .lt. 50.0D+00 ) then
        k0 = 10
      else
        k0 = 8
      end if

      if ( x .le. 12.0D+00 ) then

        do l = 1, 2
          vl = l / 3.0D+00
          vjl = 1.0D+00
          r = 1.0D+00
          do k = 1, 40
            r = -0.25D+00 * r * x2 / ( k * ( k + vl ) )
            vjl = vjl + r
            if ( abs ( r ) .lt. 1.0D-15 ) then
              go to 10
            end if
          end do

10        continue

          a0 = ( 0.5D+00 * x ) ** vl
          if ( l .eq. 1 ) then
            vj1 = a0 / gp1 * vjl
          else
            vj2 = a0 / gp2 * vjl
          end if

        end do

      else

        do l = 1, 2

          vv = vv0 * l * l
          px = 1.0D+00
          rp = 1.0D+00

          do k = 1, k0
            rp = - 0.78125D-02 * rp 
     &        * ( vv - ( 4.0D+00 * k - 3.0D+00 ) ** 2 )
     &        * ( vv - ( 4.0D+00 * k - 1.0D+00 ) ** 2 ) 
     &        / ( k * ( 2.0D+00 * k - 1.0D+00 ) * x2 )
            px = px + rp
          end do

          qx = 1.0D+00
          rq = 1.0D+00
          do k = 1, k0
            rq = - 0.78125D-02 * rq 
     &        * ( vv - ( 4.0D+00 * k - 1.0D+00 ) ** 2 )
     &        * ( vv - ( 4.0D+00 * k + 1.0D+00 ) ** 2 )
     &        / ( k * ( 2.0D+00 * k + 1.0D+00 ) * x2 )
            qx = qx + rq
          end do

          qx = 0.125D+00 * ( vv - 1.0D+00 ) * qx / x
          xk = x - ( 0.5D+00 * l / 3.0D+00 + 0.25D+00 ) * pi
          a0 = sqrt ( rp2 / x )
          ck = cos ( xk )
          sk = sin ( xk )
          if ( l .eq. 1) then
            vj1 = a0 * ( px * ck - qx * sk )
            vy1 = a0 * ( px * sk + qx * ck )
          else
            vj2 = a0 * ( px * ck - qx * sk )
            vy2 = a0 * ( px * sk + qx * ck )
          end if

        end do

      end if

      if ( x .le. 12.0D+00 ) then

        do l = 1, 2

          vl = l / 3.0D+00
          vjl = 1.0D+00
          r = 1.0D+00
          do k = 1, 40
            r = -0.25D+00 * r * x2 / ( k * ( k - vl ) )
            vjl = vjl + r
            if ( abs ( r ) .lt. 1.0D-15 ) then
              go to 20
            end if
          end do

20        continue

          b0 = ( 2.0D+00 / x ) ** vl
          if ( l .eq. 1 ) then
            uj1 = b0 * vjl / gn1
          else
             uj2 = b0 * vjl / gn2
          end if

        end do

        pv1 = pi / 3.0D+00
        pv2 = pi / 1.5D+00
        vy1 = uu0 * ( vj1 * cos ( pv1 ) - uj1 )
        vy2 = uu0 * ( vj2 * cos ( pv2 ) - uj2 )

      end if

      if ( x .le. 18.0D+00 ) then

        do l = 1, 2
          vl = l / 3.0D+00
          vil = 1.0D+00
          r = 1.0D+00
          do k = 1, 40
            r = 0.25D+00 * r * x2 / ( k * ( k + vl ) )
            vil = vil + r
            if ( abs ( r ) .lt. 1.0D-15 ) then
              go to 30
            end if
          end do

30        continue

          a0 = ( 0.5D+00 * x ) ** vl

          if ( l .eq. 1 ) then
            vi1 = a0 / gp1 * vil
          else
            vi2 = a0 / gp2 * vil
          end if

        end do

      else

        c0 = exp ( x ) / sqrt ( 2.0D+00 * pi * x )

        do l = 1, 2
          vv = vv0 * l * l
          vsl = 1.0D+00
          r = 1.0D+00
          do k = 1, k0
            r = - 0.125D+00 * r 
     &        * ( vv - ( 2.0D+00 * k - 1.0D+00 ) ** 2 ) / ( k * x )
            vsl = vsl + r
          end do
          if ( l .eq. 1 ) then
            vi1 = c0 * vsl
          else
            vi2 = c0 * vsl
          end if
        end do

      end if

      if ( x .le. 9.0D+00 ) then

        do l = 1, 2
          vl = l / 3.0D+00
          if ( l .eq. 1 ) then
            gn = gn1
          else
            gn = gn2
          end if
          a0 = ( 2.0D+00 / x ) ** vl / gn
          sum = 1.0D+00
          r = 1.0D+00
          do k = 1, 60
            r = 0.25D+00 * r * x2 / ( k * ( k - vl ) )
            sum = sum + r
            if ( abs ( r ) .lt. 1.0D-15 ) then
              go to 40
            end if
          end do

40        continue

          if ( l .eq. 1 ) then
            vk1 = 0.5D+00 * uu0 * pi * ( sum * a0 - vi1 )
          else
            vk2 = 0.5D+00 * uu0 * pi * ( sum * a0 - vi2 )
          end if

        end do

      else

        c0 = exp ( - x ) * sqrt ( 0.5D+00 * pi / x )

        do l = 1, 2
          vv = vv0 * l * l
          sum = 1.0D+00
          r = 1.0D+00
          do k = 1, k0
            r = 0.125D+00 * r * ( vv - ( 2.0D+00 * k - 1.0D+00 ) ** 2 ) 
     &        / ( k * x )
            sum = sum + r
          end do
          if ( l .eq. 1 ) then
            vk1 = c0 * sum
          else
            vk2 = c0 * sum
          end if
        end do

      end if

      return
      end
      subroutine aswfa ( m, n, c, x, kd, cv, s1f, s1d )

c*********************************************************************72
c
cc ASWFA computes prolate and oblate spheroidal angular functions of the first kind.
c
c  Licensing:
c
c    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
c    they give permission to incorporate this routine into a user program 
c    provided that the copyright is acknowledged.
c
c  Modified:
c
c    14 March 2012
c
c  Author:
c
c    Shanjie Zhang, Jianming Jin
c
c  Reference:
c
c    Shanjie Zhang, Jianming Jin,
c    Computation of Special Functions,
c    Wiley, 1996,
c    ISBN: 0-471-11963-6,
c    LC: QA351.C45.
c
c  Parameters:
c
c    Input, integer M, the mode parameter.
c
c    Input, integer N, the mode parameter, with N = M, M+1, ...
c
c    Input, double precision C, the spheroidal parameter.
c
c    Input, double precision X, the argument of the angular function.
c    |X| < 1.0.
c
c    Input, integer KD, the function code.
c    1, the prolate function.
c    -1, the oblate function.
c
c    Input, double precision CV, the characteristic value.
c
c    Output, double precision S1F, S1D, the angular function of the first
c    kind and its derivative.
c
      implicit none

      double precision a0
      double precision c
      double precision ck(200)
      double precision cv
      double precision d0
      double precision d1
      double precision df(200)
      double precision eps
      integer ip
      integer k
      integer kd
      integer m
      integer n
      integer nm
      integer nm2
      double precision r
      double precision s1d
      double precision s1f
      double precision su1
      double precision su2
      double precision x
      double precision x0
      double precision x1

      eps = 1.0D-14
      x0 = x
      x = abs ( x )

      if ( n - m .eq. 2 * int ( ( n - m ) / 2 ) ) then
        ip = 0
      else
        ip = 1
      end if

      nm = 10 + int ( ( n - m ) / 2 + c )
      nm2 = nm / 2 - 2 
      call sdmn ( m, n, c, cv, kd, df )
      call sckb ( m, n, c, df, ck )
      x1 = 1.0D+00 - x * x

      if ( m .eq. 0 .and. x1 .eq. 0.0D+00 ) then
        a0 = 1.0D+00
      else
        a0 = x1 ** ( 0.5D+00 * m )
      end if

      su1 = ck(1)
      do k = 1, nm2
        r = ck(k+1) * x1 ** k
        su1 = su1 + r
        if ( k .ge. 10 .and. abs ( r / su1 ) .lt. eps ) then
          go to 10
        end if
      end do

10    continue

      s1f = a0 * x ** ip * su1

      if ( x .eq. 1.0D+00 ) then

        if ( m .eq. 0 ) then
          s1d = ip * ck(1) - 2.0D+00 * ck(2)
        else if ( m .eq. 1 ) then
          s1d = -1.0D+100
        else if ( m .eq. 2 ) then
          s1d = -2.0D+00 * ck(1)
        else if ( 3 .le. m ) then
          s1d = 0.0D+00
        end if

      else
        d0 = ip - m / x1 * x ** ( ip + 1.0D+00 )
        d1 = -2.0D+00 * a0 * x ** ( ip + 1.0D+00 )
        su2 = ck(2)
        do k = 2, nm2
          r = k * ck(k+1) * x1 ** ( k - 1.0D+00 )
          su2 = su2 + r
          if ( 10 .le. k .and. abs ( r / su2 ) .lt. eps ) then
            go to 20
          end if
        end do

20      continue

        s1d = d0 * a0 * su1 + d1 * su2

      end if

      if ( x0 .lt. 0.0D+00 ) then
        if ( ip .eq. 0 ) then
          s1d = -s1d
        else if ( ip .eq. 1 ) then
          s1f = -s1f
        end if
      end if

      x = x0

      return
      end
      subroutine aswfb ( m, n, c, x, kd, cv, s1f, s1d )

c*********************************************************************72
c
cc ASWFB: prolate and oblate spheroidal angular functions of the first kind.
c
c  Licensing:
c
c    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
c    they give permission to incorporate this routine into a user program 
c    provided that the copyright is acknowledged.
c
c  Modified:
c
c    20 July 2012
c
c  Author:
c
c    Shanjie Zhang, Jianming Jin
c
c  Reference:
c
c    Shanjie Zhang, Jianming Jin,
c    Computation of Special Functions,
c    Wiley, 1996,
c    ISBN: 0-471-11963-6,
c    LC: QA351.C45.
c
c  Parameters:
c
c    Input, integer M, the mode parameter, m = 0, 1, 2, ...
c
c    Input, integer N, mode parameter, N = M, M+1, M+2, ...
c
c    Input, double precision C, the spheroidal parameter.
c
c    Input, double precision X, the argument, with |X| < 1.0.
c
c    Input, integer KD, the function code.
c    1, the prolate function.
c    -1, the oblate function.
c
c    Input, double precision CV, the characteristic value.
c
c    Output, double precision S1F, S1D, the angular function of the first
c    kind and its derivative.
c
      implicit none

      double precision c
      double precision cv
      double precision df(200)
      double precision eps
      integer ip
      integer k
      integer kd
      integer m
      integer mk
      integer n
      integer nm
      integer nm2
      double precision pd(0:251)
      double precision pm(0:251)
      double precision s1d
      double precision s1f
      double precision su1
      double precision sw
      double precision x

      eps = 1.0D-14

      if ( n - m .eq. 2 * int ( ( n - m ) / 2 ) ) then
        ip = 0
      else
        ip = 1
      end if

      nm = 25 + int ( ( n - m ) / 2 + c )
      nm2 = 2 * nm + m
      call sdmn ( m, n, c, cv, kd, df )
      call lpmns ( m, nm2, x, pm, pd )
      su1 = 0.0D+00
      do k = 1, nm
        mk = m + 2 * ( k - 1 ) + ip
        su1 = su1 + df(k) * pm(mk)
        if ( abs ( sw - su1 ) .lt. abs ( su1 ) * eps ) then
          go to 10
        end if
        sw = su1
      end do

10    continue

      s1f = ( -1.0D+00 ) ** m * su1

      su1 = 0.0D+00
      do k = 1, nm
        mk = m + 2 * ( k - 1 ) + ip
        su1 = su1 + df(k) * pd(mk)
        if ( abs ( sw - su1 ) .lt. abs ( su1 ) * eps ) then
          go to 20
        end if
        sw = su1
      end do

20    continue

      s1d = ( -1.0D+00 ) ** m * su1

      return
      end
      subroutine bernoa ( n, bn )

c*********************************************************************72
c
cc BERNOA computes the Bernoulli number Bn.
c
c  Licensing:
c
c    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
c    they give permission to incorporate this routine into a user program 
c    provided that the copyright is acknowledged.
c
c  Modified:
c
c    11 July 2012
c
c  Author:
c
c    Shanjie Zhang, Jianming Jin
c
c  Reference:
c
c    Shanjie Zhang, Jianming Jin,
c    Computation of Special Functions,
c    Wiley, 1996,
c    ISBN: 0-471-11963-6,
c    LC: QA351.C45.
c
c  Parameters:
c
c    Input, integer N, the index.
c
c    Output, double precision BN, the value of the N-th Bernoulli number.
c
      implicit none

      integer n

      double precision bn(0:n)
      integer j
      integer k
      integer m
      double precision r
      double precision s

      bn(0) = 1.0D+00
      bn(1) = -0.5D+00

      do m = 2, n
        s = - ( 1.0D+00 / ( m + 1.0D+00 ) - 0.5D+00 )
        do k = 2, m - 1
          r = 1.0D+00
          do j = 2, k
            r = r * ( j + m - k ) / j
          end do
        s = s - r * bn(k)
       end do
       bn(m) = s
      end do

      do m = 3, n, 2
        bn(m) = 0.0D+00
      end do

      return
      end
      subroutine bernob ( n, bn )

c*********************************************************************72
c
cc BERNOB computes the Bernoulli number Bn.
c
c  Licensing:
c
c    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
c    they give permission to incorporate this routine into a user program 
c    provided that the copyright is acknowledged.
c
c  Modified:
c
c    11 July 2012
c
c  Author:
c
c    Shanjie Zhang, Jianming Jin
c
c  Reference:
c
c    Shanjie Zhang, Jianming Jin,
c    Computation of Special Functions,
c    Wiley, 1996,
c    ISBN: 0-471-11963-6,
c    LC: QA351.C45.
c
c  Parameters:
c
c    Input, integer N, the index.
c
c    Output, double precision BN, the value of the N-th Bernoulli number.
c
      implicit none

      integer n

      double precision bn(0:n)
      integer k
      integer m
      double precision r1
      double precision r2
      double precision s
      double precision tpi

      tpi = 6.283185307179586D+00
      bn(0) = 1.0D+00
      bn(1) = -0.5D+00
      bn(2) = 1.0D+00 / 6.0D+00
      r1 = ( 2.0D+00 / tpi ) ** 2

      do m = 4, n, 2

        r1 = - r1 * ( m - 1 ) * m / ( tpi * tpi )
        r2 = 1.0D+00
 
        do k = 2, 10000
          s = ( 1.0D+00 / k ) ** m
          r2 = r2 + s
          if ( s .lt. 1.0D-15 ) then
            go to 10
          end if
        end do

10      continue

        bn(m) = r1 * r2

      end do

      return
      end
      subroutine beta ( p, q, bt )

c*********************************************************************72
c
cc BETA computes the Beta function B(p,q).
c
c  Licensing:
c
c    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
c    they give permission to incorporate this routine into a user program 
c    provided that the copyright is acknowledged.
c
c  Modified:
c
c    12 March 2012
c
c  Author:
c
c    Shanjie Zhang, Jianming Jin.
c
c  Reference:
c
c    Shanjie Zhang, Jianming Jin,
c    Computation of Special Functions,
c    Wiley, 1996,
c    ISBN: 0-471-11963-6,
c    LC: QA351.C45
c
c  Parameters:
c
c    Input, double precision P, Q, the parameters.
c    0 < P, 0 < Q.
c
c    Output, double precision BT, the value of B(P,Q).
c
      implicit none

      double precision bt
      double precision gp
      double precision gpq
      double precision gq
      double precision p
      double precision ppq
      double precision q

      call gamma ( p, gp )
      call gamma ( q, gq )
      ppq = p + q
      call gamma ( ppq, gpq )
      bt = gp * gq / gpq

      return
      end
      subroutine bjndd ( n, x, bj, dj, fj )

c*********************************************************************72
c
cc BJNDD computes Bessel functions Jn(x) and first and second derivatives.
c
c  Licensing:
c
c    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
c    they give permission to incorporate this routine into a user program 
c    provided that the copyright is acknowledged.
c
c  Modified:
c
c    11 July 2012
c
c  Author:
c
c    Shanjie Zhang, Jianming Jin
c
c  Reference:
c
c    Shanjie Zhang, Jianming Jin,
c    Computation of Special Functions,
c    Wiley, 1996,
c    ISBN: 0-471-11963-6,
c    LC: QA351.C45.
c
c  Parameters:
c
c    Input, integer N, the order.
c
c    Input, double precision X, the argument.
c
c    Output, double precision BJ(N+1), DJ(N+1), FJ(N+1), the values of 
c    Jn(x), Jn'(x) and Jn''(x) in the last entries.
c
      implicit none

      integer n

      double precision bj(n+1)
      double precision bs
      double precision dj(n+1)
      double precision f
      double precision f0
      double precision f1
      double precision fj(n+1)
      integer k
      integer m
      integer mt
      integer nt
      double precision x

      do nt = 1, 900
        mt = int ( 0.5D+00 * log10 ( 6.28D+00 * nt ) 
     &    - nt * log10 ( 1.36D+00 * abs ( x ) / nt ) )
        if ( 20 .lt. mt ) then
          go to 10
        end if
      end do

10    continue

      m = nt
      bs = 0.0D+00
      f0 = 0.0D+00
      f1 = 1.0D-35
      do k = m, 0, -1
        f = 2.0D+00 * ( k + 1.0D+00 ) * f1 / x - f0
        if ( k .le. n ) then
          bj(k+1) = f
        end if
        if ( k .eq. 2 * int ( k / 2 ) ) then
          bs = bs + 2.0D+00 * f
        end if
        f0 = f1
        f1 = f
      end do

      do k = 0, n
        bj(k+1) = bj(k+1) / ( bs - f )
      end do

      dj(1) = -bj(2)
      fj(1) = -1.0D+00 * bj(1) - dj(1) / x
      do k = 1, n
        dj(k+1) = bj(k) - k * bj(k+1) / x
        fj(k+1) = ( k * k / ( x * x ) - 1.0D+00 ) * bj(k+1) 
     &    - dj(k+1) / x
      end do

      return
      end
      subroutine cbk ( m, n, c, cv, qt, ck, bk )

c*********************************************************************72
c
cc CBK computes coefficients for oblate radial functions with small argument.
c
c  Licensing:
c
c    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
c    they give permission to incorporate this routine into a user program 
c    provided that the copyright is acknowledged.
c
c  Modified:
c
c    20 July 2012
c
c  Author:
c
c    Shanjie Zhang, Jianming Jin
c
c  Reference:
c
c    Shanjie Zhang, Jianming Jin,
c    Computation of Special Functions,
c    Wiley, 1996,
c    ISBN: 0-471-11963-6,
c    LC: QA351.C45.
c
c  Parameters:
c
c    Input, integer M, the mode parameter;  M = 0, 1, 2, ...
c
c    Input, integer N, mode parameter, N = M, M + 1, M + 2, ...
c
c    Input, double precision C, spheroidal parameter.
c
c    Input, double precision CV, the characteristic value.
c
c    Input, double precision QT, ?
c
c    Input, double precision CK(*), ?
c
c    Output, double precision BK(*), the coefficients.
c
      implicit none

      double precision bk(200)
      double precision c
      double precision ck(200)
      double precision cv
      double precision eps
      integer i
      integer i1
      integer ip
      integer j
      integer k
      integer m
      integer n
      integer n2
      integer nm
      double precision qt
      double precision r1
      double precision s1
      double precision sw
      double precision t
      double precision u(200)
      double precision v(200)
      double precision w(200)

      eps = 1.0D-14
      if ( n - m .eq. 2 * int ( ( n - m ) / 2 ) ) then
        ip = 0
      else
        ip = 1
      end if
      nm = 25 + int ( 0.5D+00 * ( n - m ) + c )
      u(1) = 0.0D+00
      n2 = nm - 2
      do j = 2, n2
        u(j) = c * c
      end do

      do j = 1, n2
        v(j) = ( 2.0D+00 * j - 1.0D+00 - ip ) 
     &    * ( 2.0D+00 * ( j - m ) - ip ) + m * ( m - 1.0D+00 ) - cv
      end do

      do j = 1, nm - 1
        w(j) = ( 2.0D+00 * j - ip ) * ( 2.0D+00 * j + 1.0D+00 - ip )
      end do

      if ( ip .eq. 0 ) then

        do k = 0, n2 - 1

          s1 = 0.0D+00
          i1 = k - m + 1

          do i = i1, nm
            if ( 0 .le. i ) then
              r1 = 1.0D+00
              do j = 1, k
                r1 = r1 * ( i + m - j ) / j
              end do
              s1 = s1 + ck(i+1) * ( 2.0D+00 * i + m ) * r1
              if ( abs ( s1 - sw ) .lt. abs ( s1 ) * eps ) then
                go to 10
              end if
              sw = s1
            end if
          end do

10        continue

          bk(k+1) = qt * s1

        end do

      else if ( ip .eq. 1 ) then

        do k = 0, n2 - 1

          s1 = 0.0D+00
          i1 = k - m + 1

          do i = i1, nm

            if ( 0 .le. i ) then

              r1 = 1.0D+00
              do j = 1, k
                r1 = r1 * ( i + m - j ) / j
              end do

              if ( 0 .lt. i ) then
                s1 = s1 + ck(i) * ( 2.0D+00 * i + m - 1 ) * r1
              end if
              s1 = s1 - ck(i+1) * ( 2.0D+00 * i + m ) * r1
              if ( abs ( s1 - sw ) .lt. abs ( s1 ) * eps ) then
                go to 20
              end if
              sw = s1

            end if

          end do

20        continue

          bk(k+1) = qt * s1

        end do

      end if

      w(1) = w(1) / v(1)
      bk(1) = bk(1) / v(1)
      do k = 2, n2
        t = v(k) - w(k-1) * u(k)
        w(k) = w(k) / t
        bk(k) = ( bk(k) - bk(k-1) * u(k) ) / t
      end do

      do k = n2 - 1, 1, -1
        bk(k) = bk(k) - w(k) * bk(k+1)
      end do

      return
      end
      subroutine cchg ( a, b, z, chg )

c*********************************************************************72
c
cc CCHG computes the confluent hypergeometric function.
c
c  Discussion:
c
c    This function computes the confluent hypergeometric function
c    M(a,b,z) with real parameters a, b and complex argument z.
c
c  Licensing:
c
c    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
c    they give permission to incorporate this routine into a user program 
c    provided that the copyright is acknowledged.
c
c  Modified:
c
c    26 July 2012
c
c  Author:
c
c    Shanjie Zhang, Jianming Jin
c
c  Reference:
c
c    Shanjie Zhang, Jianming Jin,
c    Computation of Special Functions,
c    Wiley, 1996,
c    ISBN: 0-471-11963-6,
c    LC: QA351.C45.
c
c  Parameters:
c
c    Input, double precision A, B, parameter values.
c
c    Input, complex*16 Z, the argument.
c
c    Output, complex*16 CHG, the value of M(a,b,z).
c
      implicit none

      double precision a
      double precision a0
      double precision a1
      double precision b
      double precision ba
      complex*16 cfac
      complex*16 chg
      complex*16 chg1
      complex*16 chg2
      complex*16 chw
      complex*16 ci
      complex*16 cr
      complex*16 cr1
      complex*16 cr2
      complex*16 crg
      complex*16 cs1
      complex*16 cs2
      complex*16 cy0
      complex*16 cy1
      double precision g1
      double precision g2
      double precision g3
      integer i
      integer j
      integer k
      integer la
      integer m
      integer n
      integer nl
      integer ns
      double precision phi
      double precision pi
      double precision x
      double precision x0
      double precision y
      complex*16 z
      complex*16 z0

      pi = 3.141592653589793D+00
      ci = cmplx ( 0.0D+00, 1.0D+00 )
      a0 = a
      a1 = a
      z0 = z

      if ( b .eq. 0.0D+00 .or. b .eq. - int ( abs ( b ) ) ) then
        chg = cmplx ( 1.0D+30, 0.0D+00 )
      else if ( a .eq. 0.0D+00 .or. z .eq. 0.0D+00 ) then
        chg = cmplx ( 1.0D+00, 0.0D+00 )
      else if ( a .eq. -1.0D+00 ) then
        chg = 1.0D+00 - z / b
      else if ( a .eq. b ) then
        chg = cdexp ( z )
      else if ( a - b .eq. 1.0D+00 ) then
        chg = ( 1.0D+00 + z / b ) * cdexp ( z )
      else if ( a .eq. 1.0D+00 .and. b .eq. 2.0D+00 ) then
        chg = ( cdexp ( z ) - 1.0D+00 ) / z
      else if ( a .eq. int ( a ) .and. a .lt. 0.0D+00 ) then
        m = int ( - a )
        cr = cmplx ( 1.0D+00, 0.0D+00 )
        chg = cmplx ( 1.0D+00, 0.0D+00 )
        do k = 1, m
          cr = cr * ( a + k - 1.0D+00 ) / k / ( b + k - 1.0D+00 ) * z
          chg = chg + cr
        end do
      else

        x0 = real ( z )
        if ( x0 .lt. 0.0D+00 ) then
          a = b - a
          a0 = a
          z = - z
        end if

        if ( a .lt. 2.0D+00 ) then
          nl = 0
        else
          nl = 1
          la = int ( a )
          a = a - la - 1.0D+00
        end if

        do n = 0, nl

          if ( 2.0D+00 .le. a0 ) then
            a = a + 1.0D+00
          end if

          if ( cdabs ( z ) .lt. 20.0D+00 + abs ( b ) .or. 
     &      a .lt. 0.0D+00 ) then

            chg = cmplx ( 1.0D+00, 0.0D+00 )
            crg = cmplx ( 1.0D+00, 0.0D+00 )
            do j = 1, 500
              crg = crg * ( a + j - 1.0D+00 ) 
     &          / ( j * ( b + j - 1.0D+00 ) ) * z
              chg = chg + crg
              if ( cdabs ( ( chg - chw ) / chg ) .lt. 1.0D-15 ) then
                go to 10
              end if
              chw = chg
            end do

10          continue

          else

            call gamma ( a, g1 )
            call gamma ( b, g2 )
            ba = b - a
            call gamma ( ba, g3 )
            cs1 = cmplx ( 1.0D+00, 0.0D+00 )
            cs2 = cmplx ( 1.0D+00, 0.0D+00 )
            cr1 = cmplx ( 1.0D+00, 0.0D+00 )
            cr2 = cmplx ( 1.0D+00, 0.0D+00 )

            do i = 1, 8
              cr1 = - cr1 * (     a + i - 1.0D+00 ) 
     &          * ( a - b + i ) / ( z * i )
              cr2 =   cr2 * ( b - a + i - 1.0D+00 ) 
     &          * ( i - a ) / ( z * i )
              cs1 = cs1 + cr1
              cs2 = cs2 + cr2
            end do

            x = real ( z )
            y = dimag ( z )

            if ( x .eq. 0.0D+00 .and. 0.0D+00 .le. y ) then
              phi = 0.5D+00 * pi
            else if ( x .eq. 0.0D+00 .and. y .le. 0.0D+00 ) then
              phi = -0.5D+00 * pi
            else
              phi = atan ( y / x )
            end if

            if ( -1.5D+00 * pi .lt. phi .and.
     &        phi .le. -0.5 * pi ) then
              ns = -1
            else if ( -0.5D+00 * pi .lt. phi .and. 
     &        phi .lt. 1.5D+00 * pi ) then
              ns = 1
            end if

            if ( y .eq. 0.0D+00 ) then
              cfac = cos ( pi * a )
            else
              cfac = cdexp ( ns * ci * pi * a )
            end if

            chg1 = g2 / g3 * z ** ( - a ) * cfac * cs1
            chg2 = g2 / g1 * cdexp ( z ) * z ** ( a - b ) * cs2
            chg = chg1 + chg2

          end if

          if ( n .eq. 0 ) then
            cy0 = chg
          else if ( n .eq. 1 ) then
            cy1 = chg
          end if

        end do

        if ( 2.0D+00 .le. a0 ) then
          do i = 1, la - 1
            chg = ( ( 2.0D+00 * a - b + z ) * cy1 
     &        + ( b - a ) * cy0 ) / a
            cy0 = cy1
            cy1 = chg
            a = a + 1.0D+00
          end do
        end if

        if ( x0 .lt. 0.0D+00 ) then
          chg = chg * cdexp ( - z )
        end if

      end if

      a = a1
      z = z0

      return
      end
      subroutine cerf ( z, cer, cder )

c*********************************************************************72
c
cc CERF computes the error function and derivative for a complex argument.
c
c  Licensing:
c
c    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
c    they give permission to incorporate this routine into a user program 
c    provided that the copyright is acknowledged.
c
c  Modified:
c
c    25 July 2012
c
c  Author:
c
c    Shanjie Zhang, Jianming Jin
c
c  Reference:
c
c    Shanjie Zhang, Jianming Jin,
c    Computation of Special Functions,
c    Wiley, 1996,
c    ISBN: 0-471-11963-6,
c    LC: QA351.C45.
c 
c  Parameters:
c
c    Input, complex*16, the argument.
c
c    Output, complex*16 CER, CDER, the values of erf(z) and erf'(z).
c
      implicit none

      complex*16 c0
      complex*16 cder
      complex*16 cer
      complex*16 cs
      double precision ei1
      double precision ei2
      double precision eps
      double precision er
      double precision er0
      double precision er1
      double precision er2 
      double precision eri
      double precision err
      integer k
      integer n
      double precision pi
      double precision r
      double precision ss
      double precision w
      double precision w1
      double precision w2
      double precision x
      double precision x2
      double precision y
      complex*16 z

      eps = 1.0D-12
      pi = 3.141592653589793D+00
      x = real ( z )
      y = dimag ( z )
      x2 = x * x

      if ( x .le. 3.5D+00 ) then

        er = 1.0D+00
        r = 1.0D+00
        do k = 1, 100
          r = r * x2 / ( k + 0.5D+00 )
          er = er + r
          if ( abs ( er - w ) .le. eps * abs ( er ) ) then
            go to 10
          end if
          w = er
        end do

10      continue

        c0 = 2.0D+00 / sqrt ( pi ) * x * exp ( - x2 )
        er0 = c0 * er

      else

        er = 1.0D+00
        r = 1.0D+00
        do k = 1, 12
          r = - r * ( k - 0.5D+00 ) / x2
          er = er + r
        end do
        c0 = exp ( - x2 ) / ( x * sqrt ( pi ) )
        er0 = 1.0D+00 - c0 * er

      end if

      if ( y .eq. 0.0D+00 ) then
        err = er0
        eri = 0.0D+00
      else
        cs = cos ( 2.0D+00 * x * y )
        ss = sin ( 2.0D+00 * x * y )
        er1 = exp ( - x2 ) * ( 1.0D+00 - cs ) / ( 2.0D+00 * pi * x )
        ei1 = exp ( - x2 ) * ss / ( 2.0D+00 * pi * x )
        er2 = 0.0D+00
        do n = 1, 100
          er2 = er2 + exp ( - 0.25D+00 * n * n ) 
     &      / ( n * n + 4.0D+00 * x2 ) * ( 2.0D+00 * x
     &      - 2.0D+00 * x * cosh ( n * y ) * cs 
     &      + n * sinh ( n * y ) * ss )
          if ( abs ( ( er2 - w1 ) / er2 ) .lt. eps ) then
            go to 20
          end if
          w1 = er2
        end do

20      continue

        c0 = 2.0D+00 * exp ( - x2 ) / pi
        err = er0 + er1 + c0 * er2
        ei2 = 0.0D+00
        do n = 1, 100
          ei2 = ei2 + exp ( - 0.25D+00 * n * n ) 
     &      / ( n * n + 4.0D+00 * x2 ) * ( 2.0D+00 * x
     &      * cosh ( n * y ) * ss + n * sinh ( n * y ) * cs )
          if ( abs ( ( ei2 - w2 ) / ei2 ) .lt. eps ) then
            go to 30
          end if
          w2 = ei2
        end do

30      continue

        eri = ei1 + c0 * ei2

      end if

      cer = cmplx ( err, eri )
      cder = 2.0D+00 / sqrt ( pi ) * cdexp ( - z * z )

      return
      end
      subroutine cerror ( z, cer )

c*********************************************************************72
c
cc CERROR computes the error function for a complex argument.
c
c  Licensing:
c
c    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
c    they give permission to incorporate this routine into a user program 
c    provided that the copyright is acknowledged.
c
c  Modified:
c
c    15 July 2012
c
c  Author:
c
c    Shanjie Zhang, Jianming Jin
c
c  Reference:
c
c    Shanjie Zhang, Jianming Jin,
c    Computation of Special Functions,
c    Wiley, 1996,
c    ISBN: 0-471-11963-6,
c    LC: QA351.C45.
c 
c  Parameters:
c
c    Input, complex*16 Z, the argument.
c
c    Output, complex*16 CER, the function value.
c
      implicit none

      double precision a0
      complex*16 c0
      complex*16 cer
      complex*16 cl
      complex*16 cr
      complex*16 cs
      integer k
      double precision pi
      complex*16 z
      complex*16 z1

      a0 = cdabs ( z )
      c0 = cdexp ( - z * z )
      pi = 3.141592653589793D+00
      z1 = z

      if ( real ( z ) .lt. 0.0D+00 ) then
        z1 = - z
      end if

      if ( a0 .le. 5.8D+00 ) then    

        cs = z1
        cr = z1
        do k = 1, 120
          cr = cr * z1 * z1 / ( k + 0.5D+00 )
          cs = cs + cr
          if ( cdabs ( cr / cs ) .lt. 1.0D-15 ) then
            go to 10
          end if
        end do

10      continue

        cer = 2.0D+00 * c0 * cs / dsqrt ( pi )

      else

        cl = 1.0D+00 / z1              
        cr = cl
        do k = 1, 13
          cr = -cr * ( k - 0.5D+00 ) / ( z1 * z1 )
          cl = cl + cr
          if ( cdabs ( cr / cl ) .lt. 1.0D-15 ) then
            go to 20
          end if
        end do

20      continue

        cer = 1.0D+00 - c0 * cl / sqrt ( pi )

      end if

      if ( real ( z ) .lt. 0.0D+00 ) then
        cer = -cer
      end if

      return
      end
      subroutine cerzo ( nt, zo )

c*********************************************************************72
c
cc CERZO evaluates the complex zeros of the error function.
c
c  Discussion:
c
c    The modified Newton method is used.
c
c  Licensing:
c
c    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
c    they give permission to incorporate this routine into a user program 
c    provided that the copyright is acknowledged.
c
c  Modified:
c
c    15 July 2012
c
c  Author:
c
c    Shanjie Zhang, Jianming Jin
c
c  Reference:
c
c    Shanjie Zhang, Jianming Jin,
c    Computation of Special Functions,
c    Wiley, 1996,
c    ISBN: 0-471-11963-6,
c    LC: QA351.C45.
c
c  Parameters:
c
c    Input, integer NT, the number of zeros.
c
c    Output, complex*16 ZO(NT), the zeros.
c
      implicit none

      integer nt

      integer i
      integer it
      integer j
      integer nr
      double precision pi
      double precision pu
      double precision pv
      double precision px
      double precision py
      double precision w
      double precision w0
      complex*16 z
      complex*16 zd
      complex*16 zf
      complex*16 zfd
      complex*16 zgd
      complex*16 zo(nt)
      complex*16 zp
      complex*16 zq
      complex*16 zw

      pi = 3.141592653589793D+00

      do nr = 1, nt

        pu = sqrt ( pi * ( 4.0D+00 * nr - 0.5D+00 ) )
        pv = pi * sqrt ( 2.0D+00 * nr - 0.25D+00 )
        px = 0.5D+00 * pu - 0.5D+00 * log ( pv ) / pu
        py = 0.5D+00 * pu + 0.5D+00 * log ( pv ) / pu
        z = cmplx ( px, py )
        it = 0

10      continue

        it = it + 1
        call cerf ( z, zf, zd )
        zp = cmplx ( 1.0D+00, 0.0D+00 )
        do i = 1, nr - 1
          zp = zp * ( z - zo(i) )
        end do
        zfd = zf / zp
        zq = cmplx ( 0.0D+00, 0.0D+00 )
        do i = 1, nr - 1
          zw = cmplx ( 1.0D+00, 0.0D+00 )
          do j = 1, nr - 1
            if ( j .ne. i ) then
              zw = zw * ( z - zo(j) )
            end if
          end do
          zq = zq + zw
        end do
        zgd = ( zd - zq * zfd ) / zp
        z = z - zfd / zgd
        w0 = w
        w = cdabs ( z )

        if ( it .le. 50 .and. 1.0D-11 .lt. abs ( ( w - w0 ) / w) ) then
          go to 10
        end if

        zo(nr) = z

      end do

      return
      end
      subroutine cfc ( z, zf, zd )

c*********************************************************************72
c
cc CFC computes the complex Fresnel integral C(z) and C'(z).
c
c  Licensing:
c
c    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
c    they give permission to incorporate this routine into a user program 
c    provided that the copyright is acknowledged.
c
c  Modified:
c
c    26 July 2012
c
c  Author:
c
c    Shanjie Zhang, Jianming Jin
c
c  Reference:
c
c    Shanjie Zhang, Jianming Jin,
c    Computation of Special Functions,
c    Wiley, 1996,
c    ISBN: 0-471-11963-6,
c    LC: QA351.C45.
c
c  Parameters:
c
c    Input, complex*16 Z, the argument.
c
c    Output, complex*16 ZF, ZD, the values of C(z) and C'(z).
c
      implicit none

      complex*16 c
      complex*16 cf
      complex*16 cf0
      complex*16 cf1
      complex*16 cg
      complex*16 cr
      double precision eps
      integer k
      integer m
      double precision pi
      double precision w0
      double precision wa
      double precision wa0
      complex*16 z
      complex*16 z0
      complex*16 zd
      complex*16 zf
      complex*16 zp
      complex*16 zp2

      eps = 1.0D-14
      pi = 3.141592653589793D+00
      w0 = cdabs ( z )
      zp = 0.5D+00 * pi * z * z
      zp2 = zp * zp
      z0 = cmplx ( 0.0D+00, 0.0D+00 )

      if ( z .eq. z0 ) then

        c = z0

      else if ( w0 .le. 2.5D+00 ) then

        cr = z
        c = cr
        do k = 1, 80
          cr = -0.5D+00 * cr * ( 4.0D+00 * k - 3.0D+00 )
     &      / k / ( 2.0D+00 * k - 1.0D+00 )
     &      / ( 4.0D+00 * k + 1.0D+00 ) * zp2
          c = c + cr
          wa = cdabs ( c )
          if ( abs ( ( wa - wa0 ) / wa ) .lt. eps .and. 10 .lt. k ) then
            go to 10
          end if
          wa0 = wa
        end do

 10     continue

      else if ( 2.5D+00 .lt. w0 .and. w0 .lt. 4.5D+00 ) then

        m = 85
        c = z0
        cf1 = z0
        cf0 = cmplx ( 1.0D-30, 0.0D+00 )
        do k = m, 0, -1
          cf = ( 2.0D+00 * k + 3.0D+00 ) * cf0 / zp - cf1
          if ( k .eq. int ( k / 2 ) * 2 ) then
            c = c + cf
          end if
          cf1 = cf0
          cf0 = cf
        end do
        c = cdsqrt ( 2.0D+00 / ( pi * zp ) ) * cdsin ( zp ) / cf * c

      else

        cr = cmplx ( 1.0D+00, 0.0D+00 )
        cf = cmplx ( 1.0D+00, 0.0D+00 )
        do k = 1, 20
          cr = - 0.25D+00 * cr * ( 4.0D+00 * k - 1.0D+00 )
     &      * ( 4.0D+00 * k - 3.0D+00 ) / zp2
          cf = cf + cr
        end do
        cr = 1.0D+00 / ( pi * z * z )
        cg = cr
        do k = 1, 12
          cr = - 0.25D+00 * cr * ( 4.0D+00 * k + 1.0D+00 )
     &      * ( 4.0D+00 * k - 1.0D+00 ) / zp2
          cg = cg + cr
        end do
        c = 0.5D+00 + ( cf * cdsin ( zp ) - cg * cdcos ( zp ) ) 
     &    / ( pi * z )

      end if

      zf = c
      zd = cdcos ( 0.5D+00 * pi * z * z )

      return
      end
      subroutine cfs ( z, zf, zd )

c*********************************************************************72
c
cc CFS computes the complex Fresnel integral S(z) and S'(z).
c
c  Licensing:
c
c    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
c    they give permission to incorporate this routine into a user program 
c    provided that the copyright is acknowledged.
c
c  Modified:
c
c    24 July 2012
c
c  Author:
c
c    Shanjie Zhang, Jianming Jin
c
c  Reference:
c
c    Shanjie Zhang, Jianming Jin,
c    Computation of Special Functions,
c    Wiley, 1996,
c    ISBN: 0-471-11963-6,
c    LC: QA351.C45.
c
c  Parameters:
c
c    Input, complex*16 Z, the argument.
c
c    Output, complex*16 ZF, ZD, the values of S(z) and S'(z).
c
      implicit none

      complex*16 cf
      complex*16 cf0
      complex*16 cf1
      complex*16 cg
      complex*16 cr
      double precision eps
      integer k
      integer m
      double precision pi
      complex*16 s
      double precision w0
      double precision wb
      double precision wb0
      complex*16 z
      complex*16 z0
      complex*16 zd
      complex*16 zf
      complex*16 zp
      complex*16 zp2

      eps = 1.0D-14
      pi = 3.141592653589793D+00
      w0 = cdabs ( z )
      zp = 0.5D+00 * pi * z * z
      zp2 = zp * zp
      z0 = cmplx ( 0.0D+00, 0.0D+00 )

      if ( z .eq. z0 ) then

        s = z0

      else if ( w0 .le. 2.5D+00 ) then

        s = z * zp / 3.0D+00
        cr = s
        do k = 1, 80
          cr = -0.5D+00 * cr * ( 4.0D+00 * k - 1.0D+00 ) / k 
     &      / ( 2.0D+00 * k + 1.0D+00 )
     &      / ( 4.0D+00 * k + 3.0D+00 ) * zp2
          s = s + cr
          wb = cdabs ( s )
          if ( abs ( wb - wb0 ) .lt. eps .and. 10 .lt. k ) then
            go to 10
          end if
          wb0 = wb
        end do

10      continue

      else if ( 2.5D+00 .lt. w0 .and. w0 .lt. 4.5D+00 ) then

        m = 85
        s = z0
        cf1 = z0
        cf0 = cmplx ( 1.0D-30, 0.0D+00 )
        do k = m, 0, -1
          cf = ( 2.0D+00 * k + 3.0D+00 ) * cf0 / zp - cf1
          if ( k .ne. int ( k / 2 ) * 2 ) then
            s = s + cf
          end if
          cf1 = cf0
          cf0 = cf
        end do
        s = cdsqrt ( 2.0D+00 / ( pi * zp ) ) * cdsin ( zp ) / cf * s

      else

        cr = cmplx ( 1.0D+00, 0.0D+00 )
        cf = cmplx ( 1.0D+00, 0.0D+00 )
        do k = 1, 20
          cr = -0.25D+00 * cr * ( 4.0D+00 * k - 1.0D+00 )
     &      * ( 4.0D+00 * k - 3.0D+00 ) / zp2
          cf = cf + cr
        end do
        cr = 1.0D+00 / ( pi * z * z )
        cg = cr
        do k = 1, 12
          cr = -0.25D+00 * cr * ( 4.0D+00 * k + 1.0D+00 )
     &      * ( 4.0D+00 * k - 1.0D+00 ) / zp2
          cg = cg + cr
        end do
        s = 0.5D+00 - ( cf * cdcos ( zp ) + cg * cdsin ( zp ) ) 
     &    / ( pi * z )

      end if

      zf = s
      zd = cdsin ( 0.5D+00 * pi * z * z )

      return
      end
      subroutine cgama ( x, y, kf, gr, gi )

c*********************************************************************72
c
cc CGAMA computes the Gamma function for complex argument.
c
c  Discussion:
c
c    This procedcure computes the gamma function ‚(z) or ln[‚(z)]
c    for a complex argument
c
c  Licensing:
c
c    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
c    they give permission to incorporate this routine into a user program 
c    provided that the copyright is acknowledged.
c
c  Modified:
c
c    26 July 2012
c
c  Author:
c
c    Shanjie Zhang, Jianming Jin
c
c  Reference:
c
c    Shanjie Zhang, Jianming Jin,
c    Computation of Special Functions,
c    Wiley, 1996,
c    ISBN: 0-471-11963-6,
c    LC: QA351.C45.
c
c  Parameters:
c
c    Input, double precision X, Y, the real and imaginary parts of 
c    the argument Z.
c
c    Input, integer KF, the function code.
c    0 for ln[‚(z)]
c    1 for ‚(z)
c
c    Output, double precision GR, GI, the real and imaginary parts of
c    the selected function.
c
      implicit none

      double precision a(10)
      double precision g0
      double precision gi
      double precision gi1
      double precision gr
      double precision gr1
      integer j
      integer k
      integer kf
      integer na
      double precision pi
      double precision si
      double precision sr
      double precision t
      double precision th
      double precision th1
      double precision th2
      double precision x
      double precision x0
      double precision x1
      double precision y
      double precision y1
      double precision z1
      double precision z2

      save a

      data a / 8.333333333333333D-02,-2.777777777777778D-03,
     &         7.936507936507937D-04,-5.952380952380952D-04,
     &         8.417508417508418D-04,-1.917526917526918D-03,
     &         6.410256410256410D-03,-2.955065359477124D-02,
     &         1.796443723688307D-01,-1.39243221690590D+00/

      pi = 3.141592653589793D+00

      if ( y .eq. 0.0D+00 .and. x .eq. int ( x ) .and. 
     &  x .le. 0.0D+00 ) then
        gr = 1.0D+300
        gi = 0.0D+00
        return
      else if ( x .lt. 0.0D+00 ) then
        x1 = x
        y1 = y
        x = -x
        y = -y
      end if

      x0 = x

      if ( x .le. 7.0D+00 ) then
        na = int ( 7 - x )
        x0 = x + na
      end if

      z1 = sqrt ( x0 * x0 + y * y )
      th = atan ( y / x0 )
      gr = ( x0 - 0.5D+00 ) * log ( z1 ) - th * y - x0 
     &  + 0.5D+00 * log ( 2.0D+00 * pi )
      gi = th * ( x0 - 0.5D+00 ) + y * log ( z1 ) - y

      do k = 1, 10
        t = z1 ** ( 1 - 2 * k )
        gr = gr + a(k) * t * cos ( ( 2.0D+00 * k - 1.0D+00 ) * th )
        gi = gi - a(k) * t * sin ( ( 2.0D+00 * k - 1.0D+00 ) * th )
      end do

      if ( x .le. 7.0D+00 ) then
        gr1 = 0.0D+00
        gi1 = 0.0D+00
        do j = 0, na - 1
          gr1 = gr1 + 0.5D+00 * log ( ( x + j ) ** 2 + y * y )
          gi1 = gi1 + atan ( y / ( x + j ) )
        end do
        gr = gr - gr1
        gi = gi - gi1
      end if

      if ( x1 .lt. 0.0D+00 ) then
        z1 = sqrt ( x * x + y * y )
        th1 = atan ( y / x )
        sr = - sin ( pi * x ) * cosh ( pi * y )
        si = - cos ( pi * x ) * sinh ( pi * y )
        z2 = sqrt ( sr * sr + si * si )
        th2 = atan ( si / sr )
        if ( sr .lt. 0.0D+00 ) then
          th2 = pi + th2
        end if
        gr = log ( pi / ( z1 * z2 ) ) - gr
        gi = - th1 - th2 - gi
        x = x1
        y = y1
      end if

      if ( kf .eq. 1 ) then
        g0 = exp ( gr )
        gr = g0 * cos ( gi )
        gi = g0 * sin ( gi )
      end if

      return
      end
      subroutine ch12n ( n, z, nm, chf1, chd1, chf2, chd2 )

c*********************************************************************72
c
cc CH12N computes Hankel functions of first and second kinds, complex argument.
c
c  Discussion:
c
c    Both the Hankel functions and their derivatives are computed.
c
c  Licensing:
c
c    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
c    they give permission to incorporate this routine into a user program 
c    provided that the copyright is acknowledged.
c
c  Modified:
c
c    26 July 2012
c
c  Author:
c
c    Shanjie Zhang, Jianming Jin
c
c  Reference:
c
c    Shanjie Zhang, Jianming Jin,
c    Computation of Special Functions,
c    Wiley, 1996,
c    ISBN: 0-471-11963-6,
c    LC: QA351.C45.
c 
c  Parameters:
c
c    Input, integer N, the order of the functions.
c
c    Input, complex*16 Z, the argument.
c
c    Output, integer NM, the highest order computed.
c
c    Output, complex*16 CHF1(0:n), CHD1(0:n), CHF2(0:n), CHD2(0:n), the
c    values of Hn(1)(z), Hn(1)'(z), Hn(2)(z), Hn(2)'(z).
c
      implicit none

      integer n

      complex*16 cbi(0:250)
      complex*16 cbj(0:250)
      complex*16 cbk(0:250)
      complex*16 cby(0:250)
      complex*16 cdi(0:250)
      complex*16 cdj(0:250)
      complex*16 cdk(0:250)
      complex*16 cdy(0:250)
      complex*16 chd1(0:n)
      complex*16 chd2(0:n)
      complex*16 cf1
      complex*16 cfac
      complex*16 chf1(0:n)
      complex*16 chf2(0:n)
      complex*16 ci
      integer k
      integer nm
      double precision pi
      complex*16 z
      complex*16 zi

      ci = cmplx ( 0.0D+00, 1.0D+00 )
      pi = 3.141592653589793D+00

      if ( dimag ( z ) .lt. 0.0D+00 ) then

        call cjynb ( n, z, nm, cbj, cdj, cby, cdy )

        do k = 0, nm
          chf1(k) = cbj(k) + ci * cby(k)
          chd1(k) = cdj(k) + ci * cdy(k)
        end do

        zi = ci * z
        call ciknb ( n, zi, nm, cbi, cdi, cbk, cdk )
        cfac = -2.0D+00 / ( pi * ci )

        do k = 0, nm
          chf2(k) = cfac * cbk(k)
          chd2(k) = cfac * ci * cdk(k)
          cfac = cfac * ci
        end do

      else if ( 0.0D+00 .lt. dimag ( z ) ) then

        zi = - ci * z
        call ciknb ( n, zi, nm, cbi, cdi, cbk, cdk )
        cf1 = -ci
        cfac = 2.0D+00 / ( pi * ci )

        do k = 0, nm
          chf1(k) = cfac * cbk(k)
          chd1(k) = -cfac * ci * cdk(k)
          cfac = cfac * cf1
        end do

        call cjynb ( n, z, nm, cbj, cdj, cby, cdy )

        do k = 0, nm
          chf2(k) = cbj(k) - ci * cby(k)
          chd2(k) = cdj(k) - ci * cdy(k)
        end do

      else

        call cjynb ( n, z, nm, cbj, cdj, cby, cdy )

        do k = 0, nm
          chf1(k) = cbj(k) + ci * cby(k)
          chd1(k) = cdj(k) + ci * cdy(k)
          chf2(k) = cbj(k) - ci * cby(k)
          chd2(k) = cdj(k) - ci * cdy(k)
        end do

      end if

      return
      end
      subroutine chgm ( a, b, x, hg )

c*********************************************************************72
c
cc CHGM computes the confluent hypergeometric function M(a,b,x).
c
c  Licensing:
c
c    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
c    they give permission to incorporate this routine into a user program 
c    provided that the copyright is acknowledged.
c
c  Modified:
c
c    27 July 2012
c
c  Author:
c
c    Shanjie Zhang, Jianming Jin
c
c  Reference:
c
c    Shanjie Zhang, Jianming Jin,
c    Computation of Special Functions,
c    Wiley, 1996,
c    ISBN: 0-471-11963-6,
c    LC: QA351.C45.
c
c  Parameters:
c
c    Input, double precision A, B, parameters.
c
c    Input, double precision X, the argument.
c
c    Output, double precision HG, the value of M(a,b,x).
c
      implicit none

      double precision a
      double precision a0
      double precision a1
      double precision aa
      double precision b
      double precision hg
      double precision hg1
      double precision hg2
      integer i
      integer j
      integer k
      integer la
      integer m
      integer n
      integer nl
      double precision pi
      double precision r
      double precision r1
      double precision r2
      double precision rg
      double precision sum1
      double precision sum2
      double precision ta
      double precision tb
      double precision tba
      double precision x
      double precision x0
      double precision xg
      double precision y0
      double precision y1

      pi = 3.141592653589793D+00
      a0 = a
      a1 = a
      x0 = x
      hg = 0.0D+00

      if ( b .eq. 0.0D+00 .or. b .eq. - abs ( int ( b ) ) ) then
        hg = 1.0D+300
      else if ( a .eq. 0.0D+00 .or. x .eq. 0.0D+00 ) then
        hg = 1.0D+00
      else if ( a .eq. -1.0D+00 ) then
        hg = 1.0D+00 - x / b
      else if ( a .eq. b ) then
        hg = exp ( x )
      else if ( a - b .eq. 1.0D+00 ) then
        hg = ( 1.0D+00 + x / b ) * exp ( x )
      else if ( a .eq. 1.0D+00 .and. b .eq. 2.0D+00 ) then
        hg = ( exp ( x ) - 1.0D+00 ) / x
      else if ( a .eq. int ( a ) .and. a .lt. 0.0D+00 ) then
        m = int ( - a )
        r = 1.0D+00
        hg = 1.0D+00
        do k = 1, m
          r = r * ( a + k - 1.0D+00 ) / k / ( b + k - 1.0D+00 ) * x
          hg = hg + r
        end do
      end if

      if ( hg .ne. 0.0D+00 ) then
        return
      end if

      if ( x .lt. 0.0D+00 ) then
        a = b - a
        a0 = a
        x = abs ( x )
      end if

      if ( a .lt. 2.0D+00 ) then
        nl = 0
      end if

      if ( 2.0D+00 .le. a ) then
        nl = 1
        la = int ( a )
        a = a - la - 1.0D+00
      end if

      do n = 0, nl

        if ( 2.0D+00 .le. a0 ) then
          a = a + 1.0D+00
        end if

        if ( x .le. 30.0D+00 + abs ( b ) .or. a .lt. 0.0D+00 ) then

          hg = 1.0D+00
          rg = 1.0D+00
          do j = 1, 500
            rg = rg * ( a + j - 1.0D+00 ) 
     &        / ( j * ( b + j - 1.0D+00 ) ) * x
            hg = hg + rg
            if ( abs ( rg / hg ) .lt. 1.0D-15 ) then
              go to 10
            end if
          end do

10        continue

        else

          call gamma ( a, ta )
          call gamma ( b, tb )
          xg = b - a
          call gamma ( xg, tba )
          sum1 = 1.0D+00
          sum2 = 1.0D+00
          r1 = 1.0D+00
          r2 = 1.0D+00
          do i = 1, 8
            r1 = - r1 * ( a + i - 1.0D+00 ) * ( a - b + i ) / ( x * i )
            r2 = - r2 * ( b - a + i - 1.0D+00 ) * ( a - i ) / ( x * i )
            sum1 = sum1 + r1
            sum2 = sum2 + r2
          end do
          hg1 = tb / tba * x ** ( - a ) * cos ( pi * a ) * sum1
          hg2 = tb / ta * exp ( x ) * x ** ( a - b ) * sum2
          hg = hg1 + hg2

        end if

        if ( n .eq. 0 ) then
          y0 = hg
        else if ( n .eq. 1 ) then
          y1 = hg
        end if

      end do

      if ( 2.0D+00 .le. a0 ) then
        do i = 1, la - 1
          hg = ( ( 2.0D+00 * a - b + x ) * y1 + ( b - a ) * y0 ) / a
          y0 = y1
          y1 = hg
          a = a + 1.0D+00
        end do
      end if

      if ( x0 .lt. 0.0D+00 ) then
        hg = hg * exp ( x0 )
      end if

      a = a1
      x = x0

      return
      end
      subroutine chgu ( a, b, x, hu, md )

c*********************************************************************72
c
cc CHGU computes the confluent hypergeometric function U(a,b,x).
c
c  Licensing:
c
c    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
c    they give permission to incorporate this routine into a user program 
c    provided that the copyright is acknowledged.
c
c  Modified:
c
c    27 July 2012
c
c  Author:
c
c    Shanjie Zhang, Jianming Jin
c
c  Reference:
c
c    Shanjie Zhang, Jianming Jin,
c    Computation of Special Functions,
c    Wiley, 1996,
c    ISBN: 0-471-11963-6,
c    LC: QA351.C45.
c 
c  Parameters:
c
c    Input, double precision A, B, parameters.
c
c    Input, double precision X, the argument.
c
c    Output, double precision HU, U(a,b,x).
c
c    Output, integer MD, the method code.
c
      implicit none

      double precision a
      double precision a00
      double precision aa
      double precision b
      double precision b00
      logical bl1
      logical bl2
      logical bl3
      logical bn
      double precision hu
      double precision hu1
      integer id
      integer id1
      logical il1
      logical il2
      logical il3
      integer md
      double precision x

      aa = a - b + 1.0D+00
      il1 = a .eq. int ( a ) .and. a .le. 0.0D+00
      il2 = aa .eq. int ( aa ) .and. aa .le. 0.0D+00
      il3 = abs ( a * ( a - b + 1.0D+00 ) ) / x .le. 2.0D+00
      bl1 = x .le. 5.0D+00 .or. ( x .le. 10.0D+00 .and. a .le. 2.0D+00 )
      bl2 = ( 5.0D+00 .lt. x .and. x .le. 12.5D+00 ) .and. 
     &  ( 1.0D+00 .le. a .and. a + 4.0D+00 .le. b )
      bl3 = 12.5D+00 .lt. x .and. 5.0D+00 .le. a .and. 
     &  a + 5.0D+00 .le. b
      bn = b .eq. int ( b ) .and. b .ne. 0.0D+00
      id1 = -100

      if ( b .ne. int ( b ) ) then
        call chgus ( a, b, x, hu, id1 )
        md = 1
        if ( 6 .le. id1 ) then
          return
        end if
        hu1 = hu
      end if

      if ( il1 .or. il2 .or. il3 ) then
        call chgul ( a, b, x, hu, id )
        md = 2
        if ( 6 .le. id ) then
          return
        end if
        if ( id .lt. id1 ) then
          md = 1
          id = id1
          hu = hu1
        end if
      end if

      if ( 0.0D+00 .le. a ) then
        if ( bn .and. ( bl1 .or. bl2 .or. bl3 ) ) then
          call chgubi ( a, b, x, hu, id )
          md = 3
        else
          call chguit ( a, b, x, hu, id )
          md = 4
        end if
      else
        if ( b .le. a ) then
          a00 = a
          b00 = b
          a = a - b + 1.0D+00
          b = 2.0D+00 - b
          call chguit ( a, b, x, hu, id )
          hu = x ** ( 1.0D+00 - b00 ) * hu
          a = a00
          b = b00
          md = 4
        else if ( bn .and. ( .not. il1 ) ) then
          call chgubi ( a, b, x, hu, id )
          md = 3
        end if
      end if

      if ( id .lt. 6 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'CHGU - Warning!'
        write ( *, '(a)' ) '  Accurate results were not obtained.'
      end if

      return
      end
      subroutine chgubi ( a, b, x, hu, id )

c*********************************************************************72
c
cc CHGUBI computes the confluent hypergeometric function with integer argument B.
c
c  Discussion:
c
c    This procedure computes the confluent hypergeometric function
c    U(a,b,x) with integer b ( b = Ò1,Ò2,... )
c
c  Licensing:
c
c    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
c    they give permission to incorporate this routine into a user program 
c    provided that the copyright is acknowledged.
c
c  Modified:
c
c    31 July 2012
c
c  Author:
c
c    Shanjie Zhang, Jianming Jin
c
c  Reference:
c
c    Shanjie Zhang, Jianming Jin,
c    Computation of Special Functions,
c    Wiley, 1996,
c    ISBN: 0-471-11963-6,
c    LC: QA351.C45.
c
c  Parameters:
c
c    Input, double precision A, B, parameters.
c
c    Input, double precision X, the argument.
c
c    Output, double precision HU, the value of U(a,b,x).
c
c    Output, integer ID, the estimated number of significant digits.
c
      implicit none

      double precision a
      double precision a0
      double precision a1
      double precision a2
      double precision b
      double precision da1
      double precision da2
      double precision db1
      double precision db2
      double precision el
      double precision ga
      double precision ga1
      double precision h0
      double precision hm1
      double precision hm2
      double precision hm3
      double precision hmax
      double precision hmin
      double precision hu
      double precision hu1
      double precision hu2
      double precision hw
      integer id
      integer id1
      integer id2
      integer j 
      integer k
      integer m
      integer n
      double precision ps
      double precision r
      double precision rn
      double precision rn1
      double precision s0
      double precision s1
      double precision s2
      double precision sa
      double precision sb
      double precision ua
      double precision ub
      double precision x

      id = -100
      el = 0.5772156649015329D+00
      n = abs ( b - 1 )
      rn1 = 1.0D+00
      rn = 1.0D+00
      do j = 1, n
        rn = rn * j
        if ( j .eq. n - 1 ) then
          rn1 = rn
        end if
      end do

      call psi ( a, ps )
      call gamma ( a, ga )

      if ( 0.0D+00 .lt. b ) then
        a0 = a
        a1 = a - n
        a2 = a1
        call gamma ( a1, ga1 )
        ua = ( - 1 ) ** ( n - 1 ) / ( rn * ga1 )
        ub = rn1 / ga * x ** ( - n )
      else
        a0 = a + n
        a1 = a0
        a2 = a
        call gamma ( a1, ga1 )
        ua = ( - 1 ) ** ( n - 1 ) / ( rn * ga ) * x ** n
        ub = rn1 / ga1
      end if

      hm1 = 1.0D+00
      r = 1.0D+00
      hmax = 0.0D+00
      hmin = 1.0D+300

      do k = 1, 150
        r = r * ( a0 + k - 1.0D+00 ) * x / ( ( n + k ) * k )
        hm1 = hm1 + r
        hu1 = abs ( hm1 )
        hmax = max ( hmax, hu1 )
        hmin = min ( hmin, hu1 )
        if ( abs ( hm1 - h0 ) .lt. abs ( hm1 ) * 1.0D-15 ) then
          go to 10
        end if
        h0 = hm1
      end do

10    continue

      da1 = log10 ( hmax )
      if ( hmin .ne. 0.0D+00 ) then
        da2 = log10 ( hmin )
      end if
      id = 15 - abs ( da1 - da2 )
      hm1 = hm1 * log ( x )
      s0 = 0.0D+00
      do m = 1, n
        if ( 0.0D+00 .le. b ) then
          s0 = s0 - 1.0D+00 / m
        else
          s0 = s0 + ( 1.0D+00 - a ) / ( m * ( a + m - 1.0D+00 ) )
        end if
      end do
      hm2 = ps + 2.0D+00 * el + s0
      r = 1.0D+00
      hmax = 0.0D+00
      hmin = 1.0D+300
      do k = 1, 150
        s1 = 0.0D+00
        s2 = 0.0D+00
        if ( 0.0D+00 .lt. b ) then
          do m = 1, k
            s1 = s1 - ( m + 2.0D+00 * a - 2.0D+00 ) 
     &        / ( m * ( m + a - 1.0D+00 ) )
          end do
          do m = 1, n
            s2 = s2 + 1.0D+00 / ( k + m )
          end do
        else
          do m = 1, k + n
            s1 = s1 + ( 1.0D+00 - a ) / ( m * ( m + a - 1.0D+00 ) )
          end do
          do m = 1, k
            s2 = s2 + 1.0D+00 / m
          end do
        end if
        hw = 2.0D+00 * el + ps + s1 - s2
        r = r * ( a0 + k - 1.0D+00 ) * x / ( ( n + k ) * k )
        hm2 = hm2 + r * hw
        hu2 = abs ( hm2 )
        hmax = max ( hmax, hu2 )
        hmin = min ( hmin, hu2 )

        if ( abs ( ( hm2 - h0 ) / hm2 ) .lt. 1.0D-15 ) then
          go to 20
        end if
        h0 = hm2
      end do

20    continue

      db1 = log10 ( hmax )
      if ( hmin .ne. 0.0D+00 ) then
        db2 = log10 ( hmin )
      end if
      id1 = 15 - abs ( db1 - db2 )
      id = min ( id, id1 )

      if ( n .eq. 0 ) then
        hm3 = 0.0D+00
      else
        hm3 = 1.0D+00
      end if

      r = 1.0D+00
      do k = 1, n - 1
        r = r * ( a2 + k - 1.0D+00 ) / ( ( k - n ) * k ) * x
        hm3 = hm3 + r
      end do

      sa = ua * ( hm1 + hm2 )
      sb = ub * hm3
      hu = sa + sb

      if ( sa .ne. 0.0D+00 ) then
        id1 = int ( log10 ( abs ( sa ) ) )
      end if

      if ( hu .ne. 0.0D+00 ) then
        id2 = int ( log10 ( abs ( hu ) ) )
      end if

      if ( sa * sb .lt. 0.0D+00 ) then
        id = id - abs ( id1 - id2 )
      end if

      return
      end
      subroutine chguit ( a, b, x, hu, id )

c*********************************************************************72
c
cc CHGUIT computes the hypergeometric function using Gauss-Legendre integration.
c
c  Discussion:
c
c    This procedure computes the hypergeometric function U(a,b,x) by
c    using Gaussian-Legendre integration (n = 60)
c
c  Licensing:
c
c    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
c    they give permission to incorporate this routine into a user program 
c    provided that the copyright is acknowledged.
c
c  Modified:
c
c    22 July 2012
c
c  Author:
c
c    Shanjie Zhang, Jianming Jin
c
c  Reference:
c
c    Shanjie Zhang, Jianming Jin,
c    Computation of Special Functions,
c    Wiley, 1996,
c    ISBN: 0-471-11963-6,
c    LC: QA351.C45.
c
c  Parameters:
c
c    Input, double precision A, B, parameters.
c
c    Input, double precision X, the argument.
c
c    Output, double precision HU, U(a,b,z).
c
c    Output, integer ID, the estimated number of significant digits.
c
      implicit none

      double precision a
      double precision a1
      double precision b
      double precision b1
      double precision c
      double precision d
      double precision f1
      double precision f2
      double precision g
      double precision ga
      double precision hu
      double precision hu0
      double precision hu1
      double precision hu2
      integer id
      integer j
      integer k
      integer m
      double precision s
      double precision t(30)
      double precision t1
      double precision t2
      double precision t3
      double precision t4
      double precision w(30)
      double precision x

      save t
      save w

      data t /  0.259597723012478D-01, 0.778093339495366D-01,
     &          0.129449135396945D+00, 0.180739964873425D+00,
     &          0.231543551376029D+00, 0.281722937423262D+00,
     &          0.331142848268448D+00, 0.379670056576798D+00,
     &          0.427173741583078D+00, 0.473525841761707D+00,
     &          0.518601400058570D+00, 0.562278900753945D+00,
     &          0.604440597048510D+00, 0.644972828489477D+00,
     &          0.683766327381356D+00, 0.720716513355730D+00,
     &          0.755723775306586D+00, 0.788693739932264D+00,
     &          0.819537526162146D+00, 0.848171984785930D+00,
     &          0.874519922646898D+00, 0.898510310810046D+00,
     &          0.920078476177628D+00, 0.939166276116423D+00,
     &          0.955722255839996D+00, 0.969701788765053D+00,
     &          0.981067201752598D+00, 0.989787895222222D+00,
     &          0.995840525118838D+00, 0.999210123227436D+00 /

      data w /  0.519078776312206D-01, 0.517679431749102D-01,
     &          0.514884515009810D-01, 0.510701560698557D-01,
     &          0.505141845325094D-01, 0.498220356905502D-01,
     &          0.489955754557568D-01, 0.480370318199712D-01,
     &          0.469489888489122D-01, 0.457343797161145D-01,
     &          0.443964787957872D-01, 0.429388928359356D-01,
     &          0.413655512355848D-01, 0.396806954523808D-01,
     &          0.378888675692434D-01, 0.359948980510845D-01,
     &          0.340038927249464D-01, 0.319212190192963D-01,
     &          0.297524915007890D-01, 0.275035567499248D-01,
     &          0.251804776215213D-01, 0.227895169439978D-01,
     &          0.203371207294572D-01, 0.178299010142074D-01,
     &          0.152746185967848D-01, 0.126781664768159D-01,
     &          0.100475571822880D-01, 0.738993116334531D-02,
     &          0.471272992695363D-02, 0.202681196887362D-02 /

      id = 7
      a1 = a - 1.0D+00
      b1 = b - a - 1.0D+00
      c = 12.0D+00 / x

      do m = 10, 100, 5

        hu1 = 0.0D+00
        g = 0.5D+00 * c / m
        d = g
        do j = 1, m
          s = 0.0D+00
          do k = 1, 30
            t1 = d + g * t(k)
            t2 = d - g * t(k)
            f1 = exp ( - x * t1 ) * t1 ** a1 * ( 1.0D+00 + t1 ) ** b1
            f2 = exp ( - x * t2 ) * t2 ** a1 * ( 1.0D+00 + t2 ) ** b1
            s = s + w(k) * ( f1 + f2 )
          end do
          hu1 = hu1 + s * g
          d = d + 2.0D+00 * g
        end do

        if ( abs ( 1.0D+00 - hu0 / hu1 ) .lt. 1.0D-07 ) then
          go to 10
        end if

        hu0 = hu1

      end do

10    continue

      call gamma ( a, ga )
      hu1 = hu1 / ga

      do m = 2, 10, 2
        hu2 = 0.0D+00
        g = 0.5D+00 / m
        d = g
        do j = 1, m
          s = 0.0D+00
          do k = 1, 30
            t1 = d + g * t(k)
            t2 = d - g * t(k)
            t3 = c / ( 1.0D+00 - t1 )
            t4 = c / ( 1.0D+00 - t2 ) 
            f1 = t3 * t3 / c * exp ( - x * t3 ) * t3 ** a1 
     &        * ( 1.0D+00 + t3 ) ** b1
            f2 = t4 * t4 / c * exp ( - x * t4 ) * t4 ** a1 
     &        * ( 1.0D+00 + t4 ) ** b1
            s = s + w(k) * ( f1 + f2 )
          end do
          hu2 = hu2 + s * g
          d = d + 2.0D+00 * g
        end do

        if ( abs ( 1.0D+00 - hu0 / hu2 ) .lt. 1.0D-07 ) then
          go to 20
        end if

        hu0 = hu2

      end do

20    continue

      call gamma ( a, ga )
      hu2 = hu2 / ga
      hu = hu1 + hu2

      return
      end
      subroutine chgul ( a, b, x, hu, id )

c*********************************************************************72
c
cc CHGUL: confluent hypergeometric function U(a,b,x) for large argument X.
c
c  Licensing:
c
c    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
c    they give permission to incorporate this routine into a user program 
c    provided that the copyright is acknowledged.
c
c  Modified:
c
c    22 July 2012
c
c  Author:
c
c    Shanjie Zhang, Jianming Jin
c
c  Reference:
c
c    Shanjie Zhang, Jianming Jin,
c    Computation of Special Functions,
c    Wiley, 1996,
c    ISBN: 0-471-11963-6,
c    LC: QA351.C45.
c 
c  Parameters:
c
c    Input, double precision A, B, parameters.
c
c    Input, double precision X, the argument.
c
c    Output, double precision HU, the value of U(a,b,x).
c
c    Output, integer ID, the estimated number of significant digits.
c
      implicit none

      double precision a
      double precision aa
      double precision b
      double precision hu
      integer id
      logical il1
      logical il2
      integer k
      integer nm
      double precision r
      double precision ra
      double precision r0
      double precision x

      id = -100
      aa = a - b + 1.0D+00
      il1 = ( a .eq. int ( a ) ) .and. ( a .le. 0.0D+00 )
      il2 = ( aa .eq. int ( aa ) ) .and. ( aa .le. 0.0D+00 )

      if ( il1 .or. il2 ) then

        if ( il1 ) then
          nm = abs ( a )
        end if

        if ( il2 ) then
          nm = abs ( aa )
        end if

        hu = 1.0D+00
        r = 1.0D+00
        do k = 1, nm
          r = - r * ( a + k - 1.0D+00 ) * ( a - b + k ) / ( k * x )
          hu = hu + r
        end do
        hu = x ** ( - a ) * hu
        id = 10

      else

        hu = 1.0D+00
        r = 1.0D+00
        do k = 1, 25
          r = - r * ( a + k - 1.0D+00 ) * ( a - b + k ) / ( k * x )
          ra = abs ( r )
          if ( ( 5 .lt. k .and. r0 .le. ra ) .or. ra .lt. 1.0D-15 ) then
            go to 10
          end if
          r0 = ra
          hu = hu + r
        end do

10      continue

        id = abs ( log10 ( ra ) )
        hu = x ** ( - a ) * hu

      end if

      return
      end
      subroutine chgus ( a, b, x, hu, id )

c*********************************************************************72
c
cc CHGUS: confluent hypergeometric function U(a,b,x) for small argument X.
c
c  Licensing:
c
c    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
c    they give permission to incorporate this routine into a user program 
c    provided that the copyright is acknowledged.
c
c  Modified:
c
c    27 July 2012
c
c  Author:
c
c    Shanjie Zhang, Jianming Jin
c
c  Reference:
c
c    Shanjie Zhang, Jianming Jin,
c    Computation of Special Functions,
c    Wiley, 1996,
c    ISBN: 0-471-11963-6,
c    LC: QA351.C45.
c
c  Parameters:
c
c    Input, double precision A, B, parameters.
c
c    Input, double precision X, the argument.
c
c    Output, double precision HU, U(a,b,x).
c
c    Output, integer ID, the estimated number of significant digits.
c
      implicit none

      double precision a
      double precision b
      double precision d1
      double precision d2
      double precision ga
      double precision gab
      double precision gb
      double precision gb2
      double precision h0
      double precision hmax
      double precision hmin
      double precision hu
      double precision hu0
      double precision hua
      integer id
      integer j
      double precision pi
      double precision r1
      double precision r2
      double precision x
      double precision xg1
      double precision xg2

      id = -100
      pi = 3.141592653589793D+00
      call gamma ( a, ga )
      call gamma ( b, gb )
      xg1 = 1.0D+00 + a - b
      call gamma ( xg1, gab )
      xg2 = 2.0D+00 - b
      call gamma ( xg2, gb2 )
      hu0 = pi / sin ( pi * b )
      r1 = hu0 / ( gab * gb )
      r2 = hu0 * x ** ( 1.0D+00 - b ) / ( ga * gb2 )
      hu = r1 - r2
      hmax = 0.0D+00
      hmin = 1.0D+300
      do j = 1, 150
        r1 = r1 * ( a + j - 1.0D+00 ) / ( j * ( b + j - 1.0D+00 ) ) * x
        r2 = r2 * ( a - b + j ) / ( j * ( 1.0D+00 - b + j ) ) * x
        hu = hu + r1 - r2
        hua = abs ( hu )
        hmax = max ( hmax, hua )
        hmin = min ( hmin, hua )
        if ( abs ( hu - h0 ) .lt. abs ( hu ) * 1.0D-15 ) then
          go to 10
        end if
        h0 = hu
      end do

10    continue

      d1 = log10 ( hmax )
      if ( hmin .ne. 0.0D+00 ) then
        d2 = log10 ( hmin )
      end if
      id = 15 - abs ( d1 - d2 )

      return
      end
      subroutine cik01 ( z, cbi0, cdi0, cbi1, cdi1, cbk0, cdk0, 
     &  cbk1, cdk1 )

c*********************************************************************72
c
cc CIK01: modified Bessel I0(z), I1(z), K0(z) and K1(z) for complex argument.
c
c  Discussion:
c
c    This procedure computes the modified Bessel functions I0(z), I1(z), 
c    K0(z), K1(z), and their derivatives for a complex argument.
c
c  Licensing:
c
c    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
c    they give permission to incorporate this routine into a user program 
c    provided that the copyright is acknowledged.
c
c  Modified:
c
c    31 July 2012
c
c  Author:
c
c    Shanjie Zhang, Jianming Jin
c
c  Reference:
c
c    Shanjie Zhang, Jianming Jin,
c    Computation of Special Functions,
c    Wiley, 1996,
c    ISBN: 0-471-11963-6,
c    LC: QA351.C45.
c 
c  Parameters:
c
c    Input, complex*16 Z, the argument.
c
c    Output, complex*16 CBI0, CDI0, CBI1, CDI1, CBK0, CDK0, CBK1, CDK1,
c    the values of I0(z), I0'(z), I1(z), I1'(z), K0(z), K0'(z), K1(z), K1'(z).
c
      implicit none

      double precision a(12)
      double precision a0
      double precision a1(10)
      double precision b(12)
      complex*16 ca
      complex*16 cb
      complex*16 cbi0
      complex*16 cbi1
      complex*16 cbk0
      complex*16 cbk1
      complex*16 cdi0
      complex*16 cdi1
      complex*16 cdk0
      complex*16 cdk1
      complex*16 ci
      complex*16 cr
      complex*16 cs
      complex*16 ct
      complex*16 cw
      integer k
      integer k0
      double precision pi
      double precision w0
      complex*16 z
      complex*16 z1
      complex*16 z2
      complex*16 zr
      complex*16 zr2

      save a
      save a1
      save b

      data a / 0.125D+00,7.03125D-02,
     &            7.32421875D-02,1.1215209960938D-01,
     &            2.2710800170898D-01,5.7250142097473D-01,
     &            1.7277275025845D+00,6.0740420012735D+00,
     &            2.4380529699556D+01,1.1001714026925D+02,
     &            5.5133589612202D+02,3.0380905109224D+03/

      data a1 / 0.125D+00,0.2109375D+00,
     &             1.0986328125D+00,1.1775970458984D+01,
     &             2.1461706161499D+002,5.9511522710323D+03,
     &             2.3347645606175D+05,1.2312234987631D+07,
     &             8.401390346421D+08,7.2031420482627D+10/

      data b / -0.375D+00,-1.171875D-01,
     &            -1.025390625D-01,-1.4419555664063D-01,
     &            -2.7757644653320D-01,-6.7659258842468D-01,
     &            -1.9935317337513D+00,-6.8839142681099D+00,
     &            -2.7248827311269D+01,-1.2159789187654D+02,
     &            -6.0384407670507D+02,-3.3022722944809D+03/

      pi = 3.141592653589793D+00
      ci = cmplx ( 0.0D+00, 1.0D+00 )
      a0 = cdabs ( z )
      z2 = z * z
      z1 = z

      if ( a0 .eq. 0.0D+00 ) then
        cbi0 = cmplx ( 1.0D+00, 0.0D+00 )
        cbi1 = cmplx ( 0.0D+00, 0.0D+00 )
        cdi0 = cmplx ( 0.0D+00, 0.0D+00 )
        cdi1 = cmplx ( 0.5D+00, 0.0D+00 )
        cbk0 = cmplx ( 1.0D+30, 0.0D+00 )
        cbk1 = cmplx ( 1.0D+30, 0.0D+00 )
        cdk0 = - cmplx ( 1.0D+30, 0.0D+00 )
        cdk1 = - cmplx ( 1.0D+30, 0.0D+00 )
        return
      end if

      if ( real ( z ) .lt. 0.0D+00 ) then
        z1 = -z
      end if

      if ( a0 .le. 18.0D+00 ) then

        cbi0 = cmplx ( 1.0D+00, 0.0D+00 )
        cr = cmplx ( 1.0D+00, 0.0D+00 )
        do k = 1, 50
          cr = 0.25D+00 * cr * z2 / ( k * k )
          cbi0 = cbi0 + cr
          if ( cdabs ( cr / cbi0 ) .lt. 1.0D-15 ) then
            go to 10
          end if
        end do

10      continue

        cbi1 = cmplx ( 1.0D+00, 0.0D+00 )
        cr = cmplx ( 1.0D+00, 0.0D+00 )
        do k = 1, 50
          cr = 0.25D+00 * cr * z2 / ( k * ( k + 1 ) )
          cbi1 = cbi1 + cr
          if ( cdabs ( cr / cbi1 ) .lt. 1.0D-15 ) then
            go to 20
          end if
        end do

20      continue

        cbi1 = 0.5D+00 * z1 * cbi1

      else

        if ( a0 .lt. 35.0D+00 ) then
          k0 = 12
        else if ( a0 .lt. 50.0D+00 ) then
          k0 = 9
        else
          k0 = 7
        end if

        ca = cdexp ( z1 ) / cdsqrt ( 2.0D+00 * pi * z1 )
        cbi0 = cmplx ( 1.0D+00, 0.0D+00 )
        zr = 1.0D+00 / z1
        do k = 1, k0
          cbi0 = cbi0 + a(k) * zr ** k
        end do
        cbi0 = ca * cbi0
        cbi1 = cmplx ( 1.0D+00, 0.0D+00 )
        do k = 1, k0
          cbi1 = cbi1 + b(k) * zr ** k
        end do
        cbi1 = ca * cbi1

      end if

      if ( a0 .le. 9.0D+00 ) then

        cs = cmplx ( 0.0D+00, 0.0D+00 )
        ct = - cdlog ( 0.5D+00 * z1 ) - 0.5772156649015329D+00
        w0 = 0.0D+00
        cr = cmplx ( 1.0D+00, 0.0D+00 )
        do k = 1, 50
          w0 = w0 + 1.0D+00 / k
          cr = 0.25D+00 * cr / ( k * k ) * z2
          cs = cs + cr * ( w0 + ct )
          if ( cdabs ( ( cs - cw ) / cs ) .lt. 1.0D-15 ) then
            go to 30
          end if
          cw = cs
        end do

30      continue

        cbk0 = ct + cs

      else

        cb = 0.5D+00 / z1
        zr2 = 1.0D+00 / z2
        cbk0 = cmplx ( 1.0D+00, 0.0D+00 )
        do k = 1, 10
          cbk0 = cbk0 + a1(k) * zr2 ** k
        end do
        cbk0 = cb * cbk0 / cbi0

      end if

      cbk1 = ( 1.0D+00 / z1 - cbi1 * cbk0 ) / cbi0

      if ( real ( z ) .lt. 0.0D+00 ) then

        if ( dimag ( z ) .lt. 0.0D+00 ) then
          cbk0 = cbk0 + ci * pi * cbi0
          cbk1 = - cbk1 + ci * pi * cbi1
        else
          cbk0 = cbk0 - ci * pi * cbi0
          cbk1 = - cbk1 - ci * pi * cbi1
        end if

        cbi1 = - cbi1

      end if

      cdi0 = cbi1
      cdi1 = cbi0 - 1.0D+00 / z * cbi1
      cdk0 = - cbk1
      cdk1 = - cbk0 - 1.0D+00 / z * cbk1

      return
      end
      subroutine ciklv ( v, z, cbiv, cdiv, cbkv, cdkv )

c*********************************************************************72
c
cc CIKLV: modified Bessel functions Iv(z), Kv(z), complex argument, large order.
c
c  Discussion:
c
c    This procedure computes modified Bessel functions Iv(z) and
c    Kv(z) and their derivatives with a complex argument and a large order.
c
c  Licensing:
c
c    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
c    they give permission to incorporate this routine into a user program 
c    provided that the copyright is acknowledged.
c
c  Modified:
c
c    31 July 2012
c
c  Author:
c
c    Shanjie Zhang, Jianming Jin
c
c  Reference:
c
c    Shanjie Zhang, Jianming Jin,
c    Computation of Special Functions,
c    Wiley, 1996,
c    ISBN: 0-471-11963-6,
c    LC: QA351.C45.
c
c  Parameters:
c
c    Input, double precision V, the order of Iv(z) and Kv(z).
c
c    Input, complex*16 Z, the argument.
c
c    Output, double precision CBIV, CDIV, CBKV, CDKV, the values of
c    Iv(z), Iv'(z), Kv(z), Kv'(z).
c
      implicit none

      double precision a(91)
      complex*16 cbiv
      complex*16 cbkv
      complex*16 cdiv
      complex*16 cdkv
      complex*16 ceta
      complex*16 cf(12)
      complex*16 cfi
      complex*16 cfk
      complex*16 csi
      complex*16 csk
      complex*16 ct
      complex*16 ct2
      complex*16 cws
      integer i
      integer k
      integer km
      integer l
      integer l0
      integer lf
      double precision pi
      double precision v
      double precision v0
      double precision vr
      complex*16 z

      pi = 3.141592653589793D+00
      km = 12
      call cjk ( km, a )

      do l = 1, 0, -1

        v0 = v - l
        cws = cdsqrt ( 1.0D+00 + ( z / v0 ) * ( z / v0 ) )
        ceta = cws + cdlog ( z / v0 / ( 1.0D+00 + cws ) )
        ct = 1.0D+00 / cws
        ct2 = ct * ct
        do k = 1, km
          l0 = k * ( k + 1 ) / 2 + 1
          lf = l0 + k
          cf(k) = a(lf)
          do i = lf - 1, l0, -1
            cf(k) = cf(k) * ct2 + a(i)
          end do
          cf(k) = cf(k) * ct ** k
        end do
        vr = 1.0D+00 / v0
        csi = cmplx ( 1.0D+00, 0.0D+00 )
        do k = 1, km
          csi = csi + cf(k) * vr ** k
        end do
        cbiv = cdsqrt ( ct / ( 2.0D+00 * pi * v0 ) ) 
     &    * cdexp ( v0 * ceta ) * csi
        if ( l .eq. 1 ) then
          cfi = cbiv
        end if
        csk = cmplx ( 1.0D+00, 0.0D+00 )
        do k = 1, km
          csk = csk + ( - 1 ) ** k * cf(k) * vr ** k
        end do
        cbkv = cdsqrt ( pi * ct / ( 2.0D+00 * v0 ) ) 
     &    * cdexp ( - v0 * ceta ) * csk

        if ( l .eq. 1 ) then
          cfk = cbkv
        end if

      end do

      cdiv =   cfi - v / z * cbiv
      cdkv = - cfk - v / z * cbkv

      return
      end
      subroutine cikna ( n, z, nm, cbi, cdi, cbk, cdk )

c*********************************************************************72
c
cc CIKNA: modified Bessel functions In(z), Kn(z), derivatives, complex argument.
c
c  Licensing:
c
c    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
c    they give permission to incorporate this routine into a user program 
c    provided that the copyright is acknowledged.
c
c  Modified:
c
c    30 July 2012
c
c  Author:
c
c    Shanjie Zhang, Jianming Jin
c
c  Reference:
c
c    Shanjie Zhang, Jianming Jin,
c    Computation of Special Functions,
c    Wiley, 1996,
c    ISBN: 0-471-11963-6,
c    LC: QA351.C45.
c
c  Parameters:
c
c    Input, integer N, the order of In(z) and Kn(z).
c
c    Input, complex*16 Z, the argument.
c
c    Output, integer NM, the highest order computed.
c
c    Output, complex*16 CBI((0:N), CDI(0:N), CBK(0:N), CDK(0:N), the values of
c    In(z), In'(z), Kn(z), Kn'(z).
c
      implicit none

      integer n

      double precision a0
      complex*16 cbi(0:n)
      complex*16 cbi0
      complex*16 cbi1
      complex*16 cbk(0:n)
      complex*16 cbk0
      complex*16 cbk1
      complex*16 cdi(0:n)
      complex*16 cdi0
      complex*16 cdi1
      complex*16 cdk(0:n)
      complex*16 cdk0
      complex*16 cdk1
      complex*16 cf
      complex*16 cf1
      complex*16 cf2
      complex*16 ckk
      complex*16 cs
      integer k
      integer m
      integer msta1
      integer msta2
      integer nm
      complex*16 z

      a0 = cdabs ( z )
      nm = n

      if ( a0 .lt. 1.0D-100 ) then
        do k = 0, n
          cbi(k) = cmplx ( 0.0D+00, 0.0D+00 )
          cdi(k) = cmplx ( 0.0D+00, 0.0D+00 )
          cbk(k) = - cmplx ( 1.0D+30, 0.0D+00 )
          cdk(k) = cmplx ( 1.0D+30, 0.0D+00 )
        end do
        cbi(0) = cmplx ( 1.0D+00, 0.0D+00 )
        cdi(1) = cmplx ( 0.5D+00, 0.0D+00 )
        return
      end if

      call cik01 ( z, cbi0, cdi0, cbi1, cdi1, cbk0, cdk0, cbk1, cdk1 )

      cbi(0) = cbi0
      cbi(1) = cbi1
      cbk(0) = cbk0
      cbk(1) = cbk1
      cdi(0) = cdi0
      cdi(1) = cdi1
      cdk(0) = cdk0
      cdk(1) = cdk1

      if ( n .le. 1 ) then
        return
      end if

      m = msta1 ( a0, 200 )

      if ( m .lt. n ) then
        nm = m
      else
        m = msta2 ( a0, n, 15 )
      end if

      cf2 = cmplx ( 0.0D+00, 0.0D+00 )
      cf1 = cmplx ( 1.0D-30, 0.0D+00 )
      do k = m, 0, -1
        cf = 2.0D+00 * ( k + 1.0D+00 ) / z * cf1 + cf2
        if ( k .le. nm ) then
          cbi(k) = cf
        end if
        cf2 = cf1
        cf1 = cf
      end do

      cs = cbi0 / cf
      do k = 0, nm
        cbi(k) = cs * cbi(k)
      end do

      do k = 2, nm
        if ( cdabs ( cbi(k-2) ) .lt. cdabs ( cbi(k-1) ) ) then
          ckk = ( 1.0D+00 / z - cbi(k) * cbk(k-1) ) / cbi(k-1)
        else
          ckk = ( cbi(k) * cbk(k-2) + 2.0D+00 * ( k - 1.0D+00 ) 
     &      / ( z * z ) ) / cbi(k-2)
        end if
        cbk(k) = ckk
      end do

      do k = 2, nm
        cdi(k) =   cbi(k-1) - k / z * cbi(k)
        cdk(k) = - cbk(k-1) - k / z * cbk(k)
      end do

      return
      end
      subroutine ciknb ( n, z, nm, cbi, cdi, cbk, cdk )

c*********************************************************************72
c
cc CIKNB computes modified Bessel functions In(z) and Kn(z) for complex argument.
c
c  Discussion:
c
c    This procedure also evaluates the derivatives.
c
c  Licensing:
c
c    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
c    they give permission to incorporate this routine into a user program 
c    provided that the copyright is acknowledged.
c
c  Modified:
c
c    30 July 2012
c
c  Author:
c
c    Shanjie Zhang, Jianming Jin
c
c  Reference:
c
c    Shanjie Zhang, Jianming Jin,
c    Computation of Special Functions,
c    Wiley, 1996,
c    ISBN: 0-471-11963-6,
c    LC: QA351.C45.
c
c  Parameters:
c
c    Input, integer N, the order of In(z) and Kn(z).
c
c    Input, complex*16 Z, the argument.
c
c    Output, integer NM, the highest order computed.
c
c    Output, complex*16 CBI((0:N), CDI(0:N), CBK(0:N), CDK(0:N), the values of
c    In(z), In'(z), Kn(z), Kn'(z).
c
      implicit none

      integer n

      double precision a0
      complex*16 c
      complex*16 ca0
      complex*16 cbi(0:n)
      complex*16 cbkl
      complex*16 cbs
      complex*16 cdi(0:n)
      complex*16 cbk(0:n)
      complex*16 cdk(0:n)
      complex*16 cf
      complex*16 cf0
      complex*16 cf1
      complex*16 cg
      complex*16 cg0
      complex*16 cg1
      complex*16 ci
      complex*16 cr
      complex*16 cs0
      complex*16 csk0
      double precision el
      double precision fac
      integer k
      integer k0
      integer l
      integer m
      integer msta1
      integer msta2
      integer nm
      double precision pi
      double precision vt
      complex*16 z
      complex*16 z1

      pi = 3.141592653589793D+00
      el = 0.57721566490153D+00
      a0 = cdabs ( z )
      nm = n

      if ( a0 .lt. 1.0D-100 ) then
        do k = 0, n
          cbi(k) = cmplx ( 0.0D+00, 0.0D+00 )
          cbk(k) = cmplx ( 1.0D+30, 0.0D+00 )
          cdi(k) = cmplx ( 0.0D+00, 0.0D+00 )
          cdk(k) = - cmplx ( 1.0D+30, 0.0D+00 )
        end do
        cbi(0) = cmplx ( 1.0D+00, 0.0D+00 )
        cdi(1) = cmplx ( 0.5D+00, 0.0D+00 ) 
        return
      end if

      ci = cmplx ( 0.0D+00, 1.0D+00 )

      if ( real ( z ) .lt. 0.0D+00 ) then
        z1 = -z
      else
        z1 = z
      end if

      if ( n .eq. 0 ) then
        nm = 1
      end if

      m = msta1 ( a0, 200 )

      if ( m .lt. nm ) then
        nm = m
      else
        m = msta2 ( a0, nm, 15 )
      end if

      cbs = 0.0D+00
      csk0 = 0.0D+00
      cf0 = 0.0D+00
      cf1 = 1.0D-100

      do k = m, 0, -1
        cf = 2.0D+00 * ( k + 1.0D+00 ) * cf1 / z1 + cf0
        if ( k .le. nm ) then
          cbi(k) = cf
        end if
        if ( k .ne. 0 .and. k .eq. 2 * int ( k / 2 ) ) then
          csk0 = csk0 + 4.0D+00 * cf / k
        end if
        cbs = cbs + 2.0D+00 * cf
        cf0 = cf1
        cf1 = cf
      end do

      cs0 = cdexp ( z1 ) / ( cbs - cf )

      do k = 0, nm
        cbi(k) = cs0 * cbi(k)
      end do

      if ( a0 .le. 9.0D+00 ) then

        cbk(0) = - ( cdlog ( 0.5D+00 * z1 ) + el ) * cbi(0) + cs0 * csk0
        cbk(1) = ( 1.0D+00 / z1 - cbi(1) * cbk(0) ) / cbi(0)

      else

        ca0 = cdsqrt ( pi / ( 2.0D+00 * z1 ) ) * cdexp ( -z1 )

        if ( a0 .lt. 25.0D+00 ) then
          k0 = 16
        else if ( a0 .lt. 80.0D+00 ) then
          k0 = 10
        else if ( a0 .lt. 200.0D+00 ) then
          k0 = 8
        else
          k0 = 6
        end if

        do l = 0, 1
          cbkl = 1.0D+00
          vt = 4.0D+00 * l
          cr = cmplx ( 1.0D+00, 0.0D+00 )
          do k = 1, k0
            cr = 0.125D+00 * cr 
     &        * ( vt - ( 2.0D+00 * k - 1.0D+00 ) ** 2 ) / ( k * z1 )
            cbkl = cbkl + cr
          end do
          cbk(l) = ca0 * cbkl
        end do
      end if

      cg0 = cbk(0)
      cg1 = cbk(1)
      do k = 2, nm
        cg = 2.0D+00 * ( k - 1.0D+00 ) / z1 * cg1 + cg0
        cbk(k) = cg
        cg0 = cg1
        cg1 = cg
      end do

      if ( real ( z ) .lt. 0.0D+00 ) then
        fac = 1.0D+00
        do k = 0, nm
          if ( dimag ( z ) .lt. 0.0D+00 ) then
            cbk(k) = fac * cbk(k) + ci * pi * cbi(k)
          else
            cbk(k) = fac * cbk(k) - ci * pi * cbi(k)
          end if
          cbi(k) = fac * cbi(k)
          fac = - fac
        end do
      end if

      cdi(0) = cbi(1)
      cdk(0) = -cbk(1)
      do k = 1, nm
        cdi(k) = cbi(k-1) - k / z * cbi(k)
        cdk(k) = - cbk(k-1) - k / z * cbk(k)
      end do

      return
      end
      subroutine cikva ( v, z, vm, cbi, cdi, cbk, cdk )

c*********************************************************************72
c
cc CIKVA: modified Bessel functions Iv(z), Kv(z), arbitrary order, complex.
c
c  Discussion:
c
c    Compute the modified Bessel functions Iv(z), Kv(z)
c    and their derivatives for an arbitrary order and
c    complex argument
c
c  Licensing:
c
c    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
c    they give permission to incorporate this routine into a user program 
c    provided that the copyright is acknowledged.
c
c  Modified:
c
c    31 July 2012
c
c  Author:
c
c    Shanjie Zhang, Jianming Jin
c
c  Reference:
c
c    Shanjie Zhang, Jianming Jin,
c    Computation of Special Functions,
c    Wiley, 1996,
c    ISBN: 0-471-11963-6,
c    LC: QA351.C45.
c
c  Parameters:       
c
c    Input, double precision V, the order of the functions.
c
c    Input, complex*16 Z, the argument.
c
c    Output, double precision VM, the highest order computed.
c
c    Output, double precision CBI(0:N), CDI(0:N), CBK(0:N), CDK(0:N),
c    the values of In+v0(z), In+v0'(z), Kn+v0(z), Kn+v0'(z).
c
      implicit none

      double precision a0
      complex*16 ca
      complex*16 ca1
      complex*16 ca2
      complex*16 cb
      complex*16 cbi(0:*)
      complex*16 cbi0
      complex*16 cdi(0:*)
      complex*16 cbk(0:*)
      complex*16 cbk0
      complex*16 cbk1
      complex*16 cdk(0:*)
      complex*16 cf
      complex*16 cf1
      complex*16 cf2
      complex*16 cg0
      complex*16 cg1
      complex*16 cgk
      complex*16 ci
      complex*16 ci0
      complex*16 cp
      complex*16 cr
      complex*16 cr1
      complex*16 cr2
      complex*16 cs
      complex*16 csu
      complex*16 ct
      complex*16 cvk
      double precision gan
      double precision gap
      integer k
      integer k0
      integer m
      integer msta1
      integer msta2
      integer n
      double precision pi
      double precision piv
      double precision v
      double precision v0
      double precision v0n
      double precision v0p
      double precision vm
      double precision vt
      double precision w0
      double precision ws
      double precision ws0
      complex*16 z
      complex*16 z1
      complex*16 z2

      pi = 3.141592653589793D+00
      ci = cmplx ( 0.0D+00, 1.0D+00 )
      a0 = cdabs ( z )
      z1 = z
      z2 = z * z
      n = int ( v )
      v0 = v - n
      piv = pi * v0
      vt = 4.0D+00 * v0 * v0

      if ( n .eq. 0 ) then
        n = 1
      end if

      if ( a0 .lt. 1.0D-100 ) then

        do k = 0, n
          cbi(k) = 0.0D+00
          cdi(k) = 0.0D+00
          cbk(k) = -1.0D+300
          cdk(k) = 1.0D+300
        end do

        if ( v0 .eq. 0.0D+00 ) then
          cbi(0) = cmplx ( 1.0D+00, 0.0D+00 )
          cdi(1) = cmplx ( 0.5D+00, 0.0D+00 )
        end if

        vm = v
        return

      end if

      if ( a0 .lt. 35.0D+00 ) then
        k0 = 14
      else if ( a0 .lt. 50.0D+00 ) then
        k0 = 10
      else
        k0 = 8
      end if

      if ( real ( z ) .lt. 0.0D+00 ) then
        z1 = -z
      end if

      if ( a0 .lt. 18.0D+00 ) then

        if ( v0 .eq. 0.0D+00 ) then
          ca1 = cmplx (1.0D+00, 0.0D+00 )
        else
          v0p = 1.0D+00 + v0
          call gamma ( v0p, gap )
          ca1 = ( 0.5D+00 * z1 ) ** v0 / gap
        end if

        ci0 = cmplx ( 1.0D+00, 0.0D+00 )
        cr = cmplx ( 1.0D+00, 0.0D+00 )
        do k = 1, 50
          cr = 0.25D+00 * cr * z2 / ( k * ( k + v0 ) )
          ci0 = ci0 + cr
          if ( cdabs ( cr ) .lt. cdabs ( ci0 ) * 1.0D-15 ) then
            go to 10
          end if
        end do

10      continue

        cbi0 = ci0 * ca1

      else

        ca = cdexp ( z1 ) / cdsqrt ( 2.0D+00 * pi * z1 )
        cs = cmplx ( 1.0D+00, 0.0D+00 )
        cr = cmplx ( 1.0D+00, 0.0D+00 )
        do k = 1, k0
          cr = - 0.125D+00 * cr 
     &      * ( vt - ( 2.0D+00 * k - 1.0D+00 ) ** 2 ) / ( k * z1 )
          cs = cs + cr
        end do
        cbi0 = ca * cs

      end if

      m = msta1 ( a0, 200 )

      if ( m .lt. n ) then
         n = m
      else
         m = msta2 ( a0, n, 15 )
      end if

      cf2 = cmplx ( 0.0D+00, 0.0D+00 )
      cf1 = cmplx ( 1.0D-30, 0.0D+00 )
      do k = m, 0, -1
        cf = 2.0D+00 * ( v0 + k + 1.0D+00 ) / z1 * cf1 + cf2
        if ( k .le. n ) then
          cbi(k) = cf
        end if
        cf2 = cf1
        cf1 = cf
      end do

      cs = cbi0 / cf
      do k = 0, n
        cbi(k) = cs * cbi(k)
      end do

      if ( a0 .le. 9.0D+00 ) then

        if ( v0 .eq. 0.0D+00 ) then
          ct = -cdlog ( 0.5D+00 * z1 ) - 0.5772156649015329D+00
          cs = cmplx ( 0.0D+00, 0.0D+00 )
          w0 = 0.0D+00
          cr = cmplx ( 1.0D+00, 0.0D+00 )
          do k = 1, 50
            w0 = w0 + 1.0D+00 / k
            cr = 0.25D+00 * cr / ( k * k ) * z2
            cp = cr * ( w0 + ct )
            cs = cs + cp
            if ( 10 .le. k .and. cdabs ( cp / cs ) .lt. 1.0D-15 ) then
              go to 20
            end if
          end do

20        continue

          cbk0 = ct + cs

        else

          v0n = 1.0D+00 - v0
          call gamma ( v0n, gan )
          ca2 = 1.0D+00 / ( gan * ( 0.5D+00 * z1 ) ** v0 )
          ca1 = ( 0.5D+00 * z1 ) ** v0 / gap
          csu = ca2 - ca1
          cr1 = cmplx ( 1.0D+00, 0.0D+00 )
          cr2 = cmplx ( 1.0D+00, 0.0D+00 )
          do k = 1, 50
            cr1 = 0.25D+00 * cr1 * z2 / ( k * ( k - v0 ) )
            cr2 = 0.25D+00 * cr2 * z2 / ( k * ( k + v0 ) )
            csu = csu + ca2 * cr1 - ca1 * cr2
            ws = cdabs ( csu )
            if ( 10 .le. k .and. 
     &        abs ( ws - ws0 ) / ws .lt. 1.0D-15 ) then
              go to 30
            end if
            ws0 = ws
          end do

30        continue

          cbk0 = 0.5D+00 * pi * csu / sin ( piv )

        end if

      else

        cb = cdexp ( - z1 ) * cdsqrt ( 0.5D+00 * pi / z1 )
        cs = cmplx ( 1.0D+00, 0.0D+00 )
        cr = cmplx ( 1.0D+00, 0.0D+00 )
        do k = 1, k0
          cr = 0.125D+00 * cr 
     &      * ( vt - ( 2.0D+00 * k - 1.0D+00 ) ** 2 ) / ( k * z1 )
          cs = cs + cr
        end do
        cbk0 = cb * cs

      end if

      cbk1 = ( 1.0D+00 / z1 - cbi(1) * cbk0 ) / cbi(0)
      cbk(0) = cbk0
      cbk(1) = cbk1
      cg0 = cbk0
      cg1 = cbk1

      do k = 2, n
        cgk = 2.0D+00 * ( v0 + k - 1.0D+00 ) / z1 * cg1 + cg0
        cbk(k) = cgk
        cg0 = cg1
        cg1 = cgk
      end do

      if ( real ( z ) .lt. 0.0D+00 ) then
        do k = 0, n
          cvk = cdexp ( ( k + v0 ) * pi * ci )
          if ( dimag ( z ) .lt. 0.0D+00 ) then
            cbk(k) = cvk * cbk(k) + pi * ci * cbi(k)
            cbi(k) = cbi(k) / cvk
          else if ( 0.0D+00 .lt. dimag ( z ) ) then
            cbk(k) = cbk(k) / cvk - pi * ci * cbi(k)
            cbi(k) = cvk * cbi(k)
          end if
        end do
      end if

      cdi(0) = v0 / z * cbi(0) + cbi(1)
      cdk(0) = v0 / z * cbk(0) - cbk(1)
      do k = 1, n
        cdi(k) = - ( k + v0 ) / z * cbi(k) + cbi(k-1)
        cdk(k) = - ( k + v0 ) / z * cbk(k) - cbk(k-1)
      end do

      vm = n + v0

      return
      end
      subroutine cikvb ( v, z, vm, cbi, cdi, cbk, cdk )

c*********************************************************************72
c
cc CIKVB: modified Bessel functions,Iv(z), Kv(z), arbitrary order, complex.
c
c  Licensing:
c
c    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
c    they give permission to incorporate this routine into a user program 
c    provided that the copyright is acknowledged.
c
c  Modified:
c
c    02 August 2012
c
c  Author:
c
c    Shanjie Zhang, Jianming Jin
c
c  Reference:
c
c    Shanjie Zhang, Jianming Jin,
c    Computation of Special Functions,
c    Wiley, 1996,
c    ISBN: 0-471-11963-6,
c    LC: QA351.C45.
c
c  Parameters:
c
c    Input, double precision V, the order of the functions.
c
c    Input, complex*16 Z, the argument.
c
c    Output, double precision VM, the highest order computed.
c
c    Output, double precision CBI(0:N), CDI(0:N), CBK(0:N), CDK(0:N),
c    the values of In+v0(z), In+v0'(z), Kn+v0(z), Kn+v0'(z).
c
      implicit none

      double precision a0
      complex*16 ca
      complex*16 ca1
      complex*16 ca2
      complex*16 cb
      complex*16 cbi(0:*)
      complex*16 cbi0
      complex*16 cdi(0:*)
      complex*16 cbk(0:*)
      complex*16 cbk0
      complex*16 cdk(0:*)
      complex*16 cf
      complex*16 cf1
      complex*16 cf2
      complex*16 ci
      complex*16 ci0
      complex*16 ckk
      complex*16 cp
      complex*16 cr
      complex*16 cr1
      complex*16 cr2
      complex*16 cs
      complex*16 csu
      complex*16 ct
      complex*16 cvk
      double precision gan
      double precision gap
      integer k
      integer k0
      integer m
      integer msta1
      integer msta2
      integer n
      double precision pi
      double precision piv
      double precision v
      double precision v0
      double precision v0n
      double precision v0p
      double precision vm
      double precision vt
      double precision w0
      complex*16 z
      complex*16 z1
      complex*16 z2

      z1 = z
      z2 = z * z
      a0 = cdabs ( z )
      pi = 3.141592653589793D+00
      ci = cmplx ( 0.0D+00, 1.0D+00 )
      n = int ( v )
      v0 = v - n
      piv = pi * v0
      vt = 4.0D+00 * v0 * v0

      if ( n .eq. 0 ) then
        n = 1
      end if

      if ( a0 .lt. 1.0D-100 ) then
        do k = 0, n
          cbi(k) = 0.0D+00
          cdi(k) = 0.0D+00
          cbk(k) = -1.0D+300
          cdk(k) = 1.0D+300
        end do
        if ( v0 .eq. 0.0D+00 ) then
          cbi(0) = cmplx ( 1.0D+00, 0.0D+00 )
          cdi(1) = cmplx ( 0.5D+00, 0.0D+00 )
        end if
        vm = v
        return
      end if

      if ( a0 .lt. 35.0D+00 ) then
        k0 = 14
      else if ( a0 .lt. 50.0D+00 ) then
        k0 = 10
      else
        k0 = 8
      end if

      if ( real ( z ) .lt. 0.0D+00 ) then
        z1 = -z
      end if

      if ( a0 .lt. 18.0D+00 ) then

        if ( v0 .eq. 0.0D+00 ) then
          ca1 = cmplx ( 1.0D+00, 0.0D+00 )
        else
          v0p = 1.0D+00 + v0
          call gamma ( v0p, gap )
          ca1 = ( 0.5D+00 * z1 ) ** v0 / gap
        end if

        ci0 = cmplx ( 1.0D+00, 0.0D+00 )
        cr = cmplx ( 1.0D+00, 0.0D+00 )
        do k = 1, 50
          cr = 0.25D+00 * cr * z2 / ( k * ( k + v0 ) )
          ci0 = ci0 + cr
          if ( cdabs ( cr / ci0 ) .lt. 1.0D-15 ) then
            go to 10
          end if
        end do

10      continue

        cbi0 = ci0 * ca1

      else

        ca = cdexp ( z1 ) / cdsqrt ( 2.0D+00 * pi * z1 )
        cs = cmplx ( 1.0D+00, 0.0D+00 )
        cr = cmplx ( 1.0D+00, 0.0D+00 )
        do k = 1, k0
          cr = -0.125D+00 * cr 
     &      * ( vt - ( 2.0D+00 * k - 1.0D+00 ) ** 2 ) / ( k * z1 )
          cs = cs + cr
        end do
        cbi0 = ca * cs

      end if

      m = msta1 ( a0, 200 )
      if ( m .lt. n ) then
        n = m
      else
        m = msta2 ( a0, n, 15 )
      end if

      cf2 = cmplx ( 0.0D+00, 0.0D+00 )
      cf1 = cmplx ( 1.0D-30, 0.0D+00 )
      do k = m, 0, -1
        cf = 2.0D+00 * ( v0 + k + 1.0D+00 ) / z1 * cf1 + cf2
        if ( k .le. n ) then
          cbi(k) = cf
        end if
        cf2 = cf1
        cf1 = cf
      end do
      cs = cbi0 / cf

      do k = 0, n
        cbi(k) = cs * cbi(k)
      end do

      if ( a0 .le. 9.0D+00 ) then

        if ( v0 .eq. 0.0D+00 ) then

          ct = - cdlog ( 0.5D+00 * z1 ) - 0.5772156649015329D+00
          cs = cmplx ( 0.0D+00, 0.0D+00 )
          w0 = 0.0D+00
          cr = cmplx ( 1.0D+00, 0.0D+00 )
          do k = 1, 50
            w0 = w0 + 1.0D+00 / k
            cr = 0.25D+00 * cr / ( k * k ) * z2
            cp = cr * ( w0 + ct )
            cs = cs + cp
            if ( 10 .le. k .and. cdabs ( cp / cs ) .lt. 1.0D-15 ) then
              go to 20
            end if
          end do

20        continue

          cbk0 = ct + cs

        else

          v0n = 1.0D+00 - v0
          call gamma ( v0n, gan )
          ca2 = 1.0D+00 / ( gan * ( 0.5D+00 * z1 ) ** v0 )
          ca1 = ( 0.5D+00 * z1 ) ** v0 / gap
          csu = ca2 - ca1
          cr1 = cmplx ( 1.0D+00, 0.0D+00 )
          cr2 = cmplx ( 1.0D+00, 0.0D+00 )
          do k = 1, 50
            cr1 = 0.25D+00 * cr1 * z2 / ( k * ( k - v0 ) )
            cr2 = 0.25D+00 * cr2 * z2 / ( k * ( k + v0 ) )
            cp = ca2 * cr1 - ca1 * cr2
            csu = csu + cp
            if ( 10 .le. k .and. cdabs ( cp / csu ) .lt. 1.0D-15 ) then
              go to 30
            end if
          end do

30        continue

          cbk0 = 0.5D+00 * pi * csu / sin ( piv )

        end if

      else

        cb = cdexp ( -z1 ) * cdsqrt ( 0.5D+00 * pi / z1 )
        cs = cmplx ( 1.0D+00, 0.0D+00 )
        cr = cmplx ( 1.0D+00, 0.0D+00 )
        do k = 1, k0
          cr = 0.125D+00 * cr * ( vt - ( 2.0D+00 * k - 1.0D+00 ) ** 2 ) 
     &      / ( k * z1 )
          cs = cs + cr
        end do

        cbk0 = cb * cs

      end if

      cbk(0) = cbk0

      if ( real ( z ) .lt. 0.0D+00 ) then
        do k = 0, n
          cvk = cdexp ( ( k + v0 ) * pi * ci )
          if ( dimag ( z ) .lt. 0.0D+00 ) then
            cbk(k) = cvk * cbk(k) + pi * ci * cbi(k)
            cbi(k) = cbi(k) / cvk
          else if ( 0.0D+00 .lt. dimag ( z ) ) then
            cbk(k) = cbk(k) / cvk - pi * ci * cbi(k)
            cbi(k) = cvk * cbi(k)
          end if
        end do
      end if

      do k = 1, n
        ckk = ( 1.0D+00 / z - cbi(k) * cbk(k-1) ) / cbi(k-1)
        cbk(k) = ckk
      end do

      cdi(0) = v0 / z * cbi(0) + cbi(1)
      cdk(0) = v0 / z * cbk(0) - cbk(1)
      do k = 1, n
        cdi(k) = - ( k + v0 ) / z * cbi(k) + cbi(k-1)
        cdk(k) = - ( k + v0 ) / z * cbk(k) - cbk(k-1)
      end do 
 
      vm = n + v0

      return
      end
      subroutine cisia ( x, ci, si )

c*********************************************************************72
c
cc CISIA computes cosine Ci(x) and sine integrals Si(x).
c
c  Licensing:
c
c    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
c    they give permission to incorporate this routine into a user program 
c    provided that the copyright is acknowledged.
c
c  Modified:
c
c    02 July 2012
c
c  Author:
c
c    Shanjie Zhang, Jianming Jin
c
c  Reference:
c
c    Shanjie Zhang, Jianming Jin,
c    Computation of Special Functions,
c    Wiley, 1996,
c    ISBN: 0-471-11963-6,
c    LC: QA351.C45.
c
c  Parameters:
c
c    Input, double precision X, the argument of Ci(x) and Si(x).
c
c    Output, double precision CI, SI, the values of Ci(x) and Si(x).
c
      implicit none

      double precision bj(101)
      double precision ci
      double precision el
      double precision eps
      integer k
      integer m
      double precision p2
      double precision si
      double precision x
      double precision x2
      double precision xa
      double precision xa0
      double precision xa1
      double precision xcs
      double precision xf
      double precision xg
      double precision xg1
      double precision xg2
      double precision xr
      double precision xs
      double precision xss

      p2 = 1.570796326794897D+00
      el = 0.5772156649015329D+00
      eps = 1.0D-15
      x2 = x * x

      if ( x .eq. 0.0D+00 ) then

        ci = -1.0D+300
        si = 0.0D+00

      else if ( x .le. 16.0D+00 ) then

        xr = -0.25D+00 * x2
        ci = el + log ( x ) + xr
        do k = 2, 40
          xr = -0.5D+00 * xr * ( k - 1 ) 
     &      / ( k * k * ( 2 * k - 1 ) ) * x2
          ci = ci + xr
          if ( abs ( xr ) .lt. abs ( ci ) * eps ) then
            go to 10
          end if
        end do

10      continue

        xr = x
        si = x
        do k = 1, 40
          xr = -0.5D+00 * xr * ( 2 * k - 1 ) / k 
     &      / ( 4 * k * k + 4 * k + 1 ) * x2
          si = si + xr
          if ( abs ( xr ) .lt. abs ( si ) * eps ) then
            return
          end if
        end do

      else if ( x .le. 32.0D+00 ) then

        m = int ( 47.2D+00 + 0.82D+00 * x )
        xa1 = 0.0D+00
        xa0 = 1.0D-100
        do k = m, 1, -1
          xa = 4.0D+00 * k * xa0 / x - xa1
          bj(k) = xa
          xa1 = xa0
          xa0 = xa
        end do
        xs = bj(1)
        do k = 3, m, 2
          xs = xs + 2.0D+00 * bj(k)
        end do
        bj(1) = bj(1) / xs
        do k = 2, m
          bj(k) = bj(k) / xs
        end do
        xr = 1.0D+00
        xg1 = bj(1)
        do k = 2, m
          xr = 0.25D+00 * xr * ( 2.0D+00 * k - 3.0D+00 ) ** 2
     &      / ( ( k - 1.0D+00 ) * ( 2.0D+00 * k - 1.0D+00 ) ** 2 ) * x
          xg1 = xg1 + bj(k) * xr
        end do
        xr = 1.0D+00
        xg2 = bj(1)
        do k = 2, m
          xr = 0.25D+00 * xr * ( 2.0D+00 * k - 5.0D+00 ) ** 2
     &      / ( ( k - 1.0D+00 ) * ( 2.0D+00 * k - 3.0D+00 ) ** 2 ) * x
          xg2 = xg2 + bj(k) * xr
        end do
        xcs = cos ( x / 2.0D+00 )
        xss = sin ( x / 2.0D+00 )
        ci = el + log ( x ) - x * xss * xg1 + 2 * xcs * xg2 
     &    - 2 * xcs * xcs
        si = x * xcs * xg1 + 2 * xss * xg2 - sin ( x )

      else

        xr = 1.0D+00
        xf = 1.0D+00
        do k = 1, 9
          xr = -2.0D+00 * xr * k * ( 2 * k - 1 ) / x2
          xf = xf + xr
        end do
        xr = 1.0D+00 / x
        xg = xr
        do k = 1, 8
          xr = -2.0D+00 * xr * ( 2 * k + 1 ) * k / x2
          xg = xg + xr
        end do
        ci = xf * sin ( x ) / x - xg * cos ( x ) / x
        si = p2 - xf * cos ( x ) / x - xg * sin ( x ) / x

      end if

      return
      end
      subroutine cisib ( x, ci, si )

c*********************************************************************72
c
cc CISIB computes cosine and sine integrals.
c
c  Licensing:
c
c    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
c    they give permission to incorporate this routine into a user program 
c    provided that the copyright is acknowledged.
c
c  Modified:
c
c    20 March 2012
c
c  Author:
c
c    Shanjie Zhang, Jianming Jin
c
c  Reference:
c
c    Shanjie Zhang, Jianming Jin,
c    Computation of Special Functions,
c    Wiley, 1996,
c    ISBN: 0-471-11963-6,
c    LC: QA351.C45.
c
c  Parameters:
c
c    Input, double precision X, the argument of Ci(x) and Si(x).
c
c    Output, double precision CI, SI, the values of Ci(x) and Si(x).
c
      implicit none

      double precision ci
      double precision fx
      double precision gx
      double precision si
      double precision x
      double precision x2

      x2 = x * x

      if ( x .eq. 0.0D+00 ) then

        ci = -1.0D+300
        si = 0.0D+00

      else if ( x .le. 1.0D+00 ) then

        ci = (((( -3.0D-08        * x2 
     &           + 3.10D-06     ) * x2 
     &           - 2.3148D-04   ) * x2 
     &           + 1.041667D-02 ) * x2 
     &           - 0.25D+00     ) * x2 + 0.577215665D+00 + log ( x )

         si = (((( 3.1D-07        * x2
     &           - 2.834D-05    ) * x2
     &           + 1.66667D-03  ) * x2
     &           - 5.555556D-02 ) * x2 + 1.0D+00 ) * x

      else

        fx = (((( x2
     &    + 38.027264D+00  ) * x2
     &    + 265.187033D+00 ) * x2
     &    + 335.67732D+00  ) * x2
     &    + 38.102495D+00  ) /
     &    (((( x2
     &    + 40.021433D+00  ) * x2
     &    + 322.624911D+00 ) * x2
     &    + 570.23628D+00  ) * x2
     &    + 157.105423D+00 )

        gx = (((( x2
     &    + 42.242855D+00  ) * x2
     &    + 302.757865D+00 ) * x2
     &    + 352.018498D+00 ) * x2
     &    + 21.821899D+00 ) /
     &    (((( x2
     &    + 48.196927D+00   ) * x2
     &    + 482.485984D+00  ) * x2
     &    + 1114.978885D+00 ) * x2
     &    + 449.690326D+00  ) / x

        ci = fx * sin ( x ) / x - gx * cos ( x ) / x

        si = 1.570796327D+00 - fx * cos ( x ) / x - gx * sin ( x ) / x

      end if

      return
      end
      subroutine cjk ( km, a )

c*********************************************************************72
c
cc CJK computes asymptotic expansion coefficients for Bessel functions of large order.
c
c  Licensing:
c
c    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
c    they give permission to incorporate this routine into a user program 
c    provided that the copyright is acknowledged.
c
c  Modified:
c
c    13 July 2012
c
c  Author:
c
c    Shanjie Zhang, Jianming Jin
c
c  Reference:
c
c    Shanjie Zhang, Jianming Jin,
c    Computation of Special Functions,
c    Wiley, 1996,
c    ISBN: 0-471-11963-6,
c    LC: QA351.C45.
c
c  Parameters:
c
c    Input, integer KM, the maximum value of K.
c
c    Output, double precision A(L), the value of Cj(k) where j and k are 
c    related to L by L = j+1+[k*(k+1)]/2; j,k = 0,1,...,Km.
c
      implicit none

      double precision a(*)
      double precision f
      double precision f0
      double precision g
      double precision g0
      integer j
      integer k
      integer km
      integer l1
      integer l2
      integer l3
      integer l4

      a(1) = 1.0D+00
      f0 = 1.0D+00
      g0 = 1.0D+00
      do k = 0, km - 1
        l1 = ( k + 1 ) * ( k + 2 ) / 2 + 1
        l2 = ( k + 1 ) * ( k + 2 ) / 2 + k + 2
        f = ( 0.5D+00 * k + 0.125D+00 / ( k + 1 ) ) * f0
        g = - ( 1.5D+00 * k + 0.625D+00 
     &    / ( 3.0D+00 * ( k + 1.0D+00 ) ) ) * g0
        a(l1) = f
        a(l2) = g
        f0 = f
        g0 = g
      end do

      do k = 1, km - 1
        do j = 1, k
          l3 = k * ( k + 1 ) / 2 + j + 1
          l4 = ( k + 1 ) * ( k + 2 ) / 2 + j + 1
          a(l4) = ( j + 0.5D+00 * k + 0.125D+00 
     &      / ( 2.0D+00 * j + k + 1.0D+00 ) ) * a(l3)
     &      - ( j + 0.5D+00 * k - 1.0D+00 + 0.625D+00 
     &      / ( 2.0D+00 * j + k + 1.0D+00 ) ) * a(l3-1)
        end do
      end do

      return
      end
      subroutine cjy01 ( z, cbj0, cdj0, cbj1, cdj1, cby0, cdy0, 
     & cby1, cdy1 )

c*********************************************************************72
c
cc CJY01: Bessel functions, derivatives, J0(z), J1(z), Y0(z), Y1(z), complex argument.
c
c  Licensing:
c
c    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
c    they give permission to incorporate this routine into a user program 
c    provided that the copyright is acknowledged.
c
c  Modified:
c
c    02 August 2012
c
c  Author:
c
c    Shanjie Zhang, Jianming Jin
c
c  Reference:
c
c    Shanjie Zhang, Jianming Jin,
c    Computation of Special Functions,
c    Wiley, 1996,
c    ISBN: 0-471-11963-6,
c    LC: QA351.C45.
c 
c  Parameters:
c
c    Input, complex*16 Z, the argument.
c
c    Output, complex*16 CBJ0, CDJ0, CBJ1, CDJ1, CBY0, CDY0, CBY1, CDY1,
c    the values of J0(z), J0'(z), J1(z), J1'(z), Y0(z), Y0'(z), Y1(z), Y1'(z).
c
      implicit none

      double precision a(12)
      double precision a0
      double precision a1(12)
      double precision b(12)
      double precision b1(12)
      complex*16 cbj0
      complex*16 cbj1
      complex*16 cby0
      complex*16 cby1
      complex*16 cdj0
      complex*16 cdj1
      complex*16 cdy0
      complex*16 cdy1
      complex*16 ci
      complex*16 cp
      complex*16 cp0
      complex*16 cp1
      complex*16 cq0
      complex*16 cq1
      complex*16 cr
      complex*16 cs
      complex*16 ct1
      complex*16 ct2
      complex*16 cu
      double precision el
      integer k
      integer k0
      double precision pi
      double precision rp2
      double precision w0
      double precision w1
      complex*16 z
      complex*16 z1
      complex*16 z2

      save a
      save a1
      save b
      save b1

      data a     /-0.703125D-01,0.112152099609375D+00,
     &            -0.5725014209747314D+00,0.6074042001273483D+01,
     &            -0.1100171402692467D+03,0.3038090510922384D+04,
     &            -0.1188384262567832D+06,0.6252951493434797D+07,
     &            -0.4259392165047669D+09,0.3646840080706556D+11,
     &            -0.3833534661393944D+13,0.4854014686852901D+15/

      data a1     /0.1171875D+00,-0.144195556640625D+00,
     &             0.6765925884246826D+00,-0.6883914268109947D+01,
     &             0.1215978918765359D+03,-0.3302272294480852D+04,
     &             0.1276412726461746D+06,-0.6656367718817688D+07,
     &             0.4502786003050393D+09,-0.3833857520742790D+11,
     &             0.4011838599133198D+13,-0.5060568503314727D+15/

      data b     / 0.732421875D-01,-0.2271080017089844D+00,
     &             0.1727727502584457D+01,-0.2438052969955606D+02,
     &             0.5513358961220206D+03,-0.1825775547429318D+05,
     &             0.8328593040162893D+06,-0.5006958953198893D+08,
     &             0.3836255180230433D+10,-0.3649010818849833D+12,
     &             0.4218971570284096D+14,-0.5827244631566907D+16/

      data b1     /-0.1025390625D+00,0.2775764465332031D+00,
     &             -0.1993531733751297D+01,0.2724882731126854D+02,
     &             -0.6038440767050702D+03,0.1971837591223663D+05,
     &             -0.8902978767070678D+06,0.5310411010968522D+08,
     &             -0.4043620325107754D+10,0.3827011346598605D+12,
     &             -0.4406481417852278D+14,0.6065091351222699D+16/

      pi = 3.141592653589793D+00
      el = 0.5772156649015329D+00
      rp2 = 2.0D+00 / pi
      ci = cmplx ( 0.0D+00, 1.0D+00 )
      a0 = cdabs ( z )
      z2 = z * z
      z1 = z

      if ( a0 .eq. 0.0D+00 ) then
        cbj0 = cmplx ( 1.0D+00, 0.0D+00 )
        cbj1 = cmplx ( 0.0D+00, 0.0D+00 )
        cdj0 = cmplx ( 0.0D+00, 0.0D+00 )
        cdj1 = cmplx ( 0.5D+00, 0.0D+00 )
        cby0 = - cmplx ( 1.0D+30, 0.0D+00 )
        cby1 = - cmplx ( 1.0D+30, 0.0D+00 )
        cdy0 = cmplx ( 1.0D+30, 0.0D+00 )
        cdy1 = cmplx ( 1.0D+30, 0.0D+00 )
        return
      end if

      if ( real ( z ) .lt. 0.0D+00 ) then
        z1 = -z
      end if

      if ( a0 .le. 12.0D+00 ) then

        cbj0 = cmplx ( 1.0D+00, 0.0D+00 )
        cr = cmplx ( 1.0D+00, 0.0D+00 )
        do k = 1, 40
          cr = -0.25D+00 * cr * z2 / ( k * k )
          cbj0 = cbj0 + cr
          if ( cdabs ( cr ) .lt. cdabs ( cbj0 ) * 1.0D-15 ) then
            go to 10
          end if
        end do

10      continue

        cbj1 = cmplx ( 1.0D+00, 0.0D+00 )
        cr = cmplx ( 1.0D+00, 0.0D+00 )
        do k = 1, 40
          cr = -0.25D+00 * cr * z2 / ( k * ( k + 1.0D+00 ) )
          cbj1 = cbj1 + cr
          if ( cdabs ( cr ) .lt. cdabs ( cbj1 ) * 1.0D-15 ) then
            go to 20
          end if
        end do

20      continue

        cbj1 = 0.5D+00 * z1 * cbj1
        w0 = 0.0D+00
        cr = cmplx ( 1.0D+00, 0.0D+00 )
        cs = cmplx ( 0.0D+00, 0.0D+00 )
        do k = 1, 40
          w0 = w0 + 1.0D+00 / k
          cr = -0.25D+00 * cr / ( k * k ) * z2
          cp = cr * w0
          cs = cs + cp
          if ( cdabs ( cp ) .lt. cdabs ( cs ) * 1.0D-15 ) then
            go to 30
          end if
        end do

30      continue

        cby0 = rp2 * ( cdlog ( z1 / 2.0D+00 ) + el ) * cbj0 - rp2 * cs
        w1 = 0.0D+00
        cr = cmplx ( 1.0D+00, 0.0D+00 )
        cs = cmplx ( 1.0D+00, 0.0D+00 )
        do k = 1, 40
          w1 = w1 + 1.0D+00 / k
          cr = -0.25D+00 * cr / ( k * ( k + 1 ) ) * z2
          cp = cr * ( 2.0D+00 * w1 + 1.0D+00 / ( k + 1.0D+00 ) )
          cs = cs + cp
          if ( cdabs ( cp ) .lt. cdabs ( cs ) * 1.0D-15 ) then
            go to 40
          end if
        end do

40      continue

        cby1 = rp2 * ( ( cdlog ( z1 / 2.0D+00 ) + el ) * cbj1
     &    - 1.0D+00 / z1 - 0.25D+00 * z1 * cs )

      else

        if ( a0 .lt. 35.0D+00 ) then
          k0 = 12
        else if ( a0 .lt. 50.0D+00 ) then
          k0 = 10
        else
          k0 = 8
        end if

        ct1 = z1 - 0.25D+00 * pi
        cp0 = cmplx ( 1.0D+00, 0.0D+00 )
        do k = 1, k0
          cp0 = cp0 + a(k) * z1 ** ( - 2 * k )
        end do
        cq0 = -0.125D+00 / z1
        do k = 1, k0
          cq0 = cq0 + b(k) * z1 ** ( - 2 * k - 1 )
        end do
        cu = cdsqrt ( rp2 / z1 )
        cbj0 = cu * ( cp0 * cdcos ( ct1 ) - cq0 * cdsin ( ct1 ) )
        cby0 = cu * ( cp0 * cdsin ( ct1 ) + cq0 * cdcos ( ct1 ) )
        ct2 = z1 - 0.75D+00 * pi
        cp1 = cmplx ( 1.0D+00, 0.0D+00 )
        do k = 1, k0
          cp1 = cp1 + a1(k) * z1 ** ( - 2 * k )
        end do
        cq1 = 0.375D+00 / z1
        do k = 1, k0
          cq1 = cq1 + b1(k) * z1 ** ( - 2 * k - 1 )
        end do
        cbj1 = cu * ( cp1 * cdcos ( ct2 ) - cq1 * cdsin ( ct2 ) )
        cby1 = cu * ( cp1 * cdsin ( ct2 ) + cq1 * cdcos ( ct2 ) )

      end if

      if ( real ( z ) .lt. 0.0D+00 ) then
        if ( dimag ( z ) .lt. 0.0D+00 ) then
          cby0 = cby0 - 2.0D+00 * ci * cbj0
          cby1 = - ( cby1 - 2.0D+00 * ci * cbj1 )
        else
          cby0 = cby0 + 2.0D+00 * ci * cbj0
          cby1 = - ( cby1 + 2.0D+00 * ci * cbj1 )
        end if
        cbj1 = -cbj1
      end if

      cdj0 = -cbj1
      cdj1 = cbj0 - 1.0D+00 / z * cbj1
      cdy0 = -cby1
      cdy1 = cby0 - 1.0D+00 / z * cby1

      return
      end
      subroutine cjylv ( v, z, cbjv, cdjv, cbyv, cdyv )

c*********************************************************************72
c
cc CJYLV computes Bessel functions Jv(z), Yv(z) of complex argument and large order v.
c
c  Licensing:
c
c    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
c    they give permission to incorporate this routine into a user program 
c    provided that the copyright is acknowledged.
c
c  Modified:
c
c    25 July 2012
c
c  Author:
c
c    Shanjie Zhang, Jianming Jin
c
c  Reference:
c
c    Shanjie Zhang, Jianming Jin,
c    Computation of Special Functions,
c    Wiley, 1996,
c    ISBN: 0-471-11963-6,
c    LC: QA351.C45.
c
c  Parameters:
c
c    Input, double precision V, the order of Jv(z) and Yv(z).
c
c    Input, complex*16 Z, the argument.
c
c    Output, complex*16 CBJV, CDJV, CBYV, CDYV, the values of Jv(z), Jv'(z),
c    Yv(z), Yv'(z).
c
      implicit none

      double precision a(91)
      complex*16 cbjv
      complex*16 cbyv
      complex*16 cdjv
      complex*16 cdyv
      complex*16 ceta
      complex*16 cf(12)
      complex*16 cfj
      complex*16 cfy
      complex*16 csj
      complex*16 csy
      complex*16 ct
      complex*16 ct2
      complex*16 cws
      integer i
      integer k
      integer km
      integer l
      integer l0
      integer lf
      double precision pi
      double precision v
      double precision v0
      double precision vr
      complex*16 z

      km = 12
      call cjk ( km, a )
      pi = 3.141592653589793D+00

      do l = 1, 0, -1

        v0 = v - l
        cws = cdsqrt ( 1.0D+00 - ( z / v0 ) * ( z / v0 ) )
        ceta = cws + cdlog ( z / v0 / ( 1.0D+00 + cws ) )
        ct = 1.0D+00 / cws
        ct2 = ct * ct

        do k = 1, km
          l0 = k * ( k + 1 ) / 2 + 1
          lf = l0 + k
          cf(k) = a(lf)
          do i = lf - 1, l0, -1
            cf(k) = cf(k) * ct2 + a(i)
          end do
          cf(k) = cf(k) * ct ** k
        end do

        vr = 1.0D+00 / v0
        csj = cmplx ( 1.0D+00, 0.0D+00 )
        do k = 1, km
          csj = csj + cf(k) * vr ** k
        end do
        cbjv = cdsqrt ( ct / ( 2.0D+00 * pi * v0 ) ) 
     &    * cdexp ( v0 * ceta ) * csj
        if ( l .eq. 1 ) then
          cfj = cbjv
        end if
        csy = cmplx ( 1.0D+00, 0.0D+00 )
        do k = 1, km
          csy = csy + ( -1.0D+00 ) ** k * cf(k) * vr ** k
        end do
        cbyv = - cdsqrt ( 2.0D+00 * ct / ( pi * v0 ) ) 
     &    * cdexp ( - v0 * ceta ) * csy
        if ( l .eq. 1 ) then
          cfy = cbyv
        end if

      end do

      cdjv = - v / z * cbjv + cfj
      cdyv = - v / z * cbyv + cfy

      return
      end
      subroutine cjyna ( n, z, nm, cbj, cdj, cby, cdy )

c*********************************************************************72
c
cc CJYNA: Bessel functions and derivatives, Jn(z) and Yn(z) of complex argument.
c
c  Licensing:
c
c    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
c    they give permission to incorporate this routine into a user program 
c    provided that the copyright is acknowledged.
c
c  Modified:
c
c    02 August 2012
c
c  Author:
c
c    Shanjie Zhang, Jianming Jin
c
c  Reference:
c
c    Shanjie Zhang, Jianming Jin,
c    Computation of Special Functions,
c    Wiley, 1996,
c    ISBN: 0-471-11963-6,
c    LC: QA351.C45.
c
c  Parameters:
c
c    Input, integer N, the order of Jn(z) and Yn(z).
c
c    Input, complex*16 Z, the argument of Jn(z) and Yn(z).
c
c    Output, integer NM, the highest order computed.
c
c    Output, complex*16, CBJ(0:N), CDJ(0:N), CBY(0:N), CDY(0:N),
c    the values of Jn(z), Jn'(z), Yn(z), Yn'(z).
c
      implicit none

      integer n

      double precision a0
      complex*16 cbj(0:n)
      complex*16 cbj0
      complex*16 cbj1
      complex*16 cby(0:n)
      complex*16 cby0
      complex*16 cby1
      complex*16 cdj(0:n)
      complex*16 cdj0
      complex*16 cdj1
      complex*16 cdy(0:n)
      complex*16 cdy0
      complex*16 cdy1
      complex*16 cf
      complex*16 cf1
      complex*16 cf2
      complex*16 cg0
      complex*16 cg1
      complex*16 ch0
      complex*16 ch1
      complex*16 ch2
      complex*16 cj0
      complex*16 cj1
      complex*16 cjk
      complex*16 cp11
      complex*16 cp12
      complex*16 cp21
      complex*16 cp22
      complex*16 cs
      complex*16 cyk
      complex*16 cyl1
      complex*16 cyl2
      complex*16 cylk
      integer k
      integer lb
      integer lb0
      integer m
      integer msta1
      integer msta2
      integer nm
      double precision pi
      double precision wa
      double precision ya0
      double precision ya1
      double precision yak
      complex*16 z

      pi = 3.141592653589793D+00
      a0 = cdabs ( z )
      nm = n

      if ( a0 .lt. 1.0D-100 ) then
        do k = 0, n
          cbj(k) = cmplx ( 0.0D+00, 0.0D+00 )
          cdj(k) = cmplx ( 0.0D+00, 0.0D+00 )
          cby(k) = - cmplx ( 1.0D+30, 0.0D+00 )
          cdy(k) = cmplx ( 1.0D+30, 0.0D+00 )
        end do
        cbj(0) = cmplx ( 1.0D+00, 0.0D+00 )
        cdj(1) = cmplx ( 0.5D+00, 0.0D+00 )
        return
      end if

      call cjy01 ( z, cbj0, cdj0, cbj1, cdj1, cby0, cdy0, cby1, cdy1 )
      cbj(0) = cbj0
      cbj(1) = cbj1
      cby(0) = cby0
      cby(1) = cby1
      cdj(0) = cdj0
      cdj(1) = cdj1
      cdy(0) = cdy0
      cdy(1) = cdy1

      if ( n .le. 1 ) then
        return
      end if

      if ( n .lt. int ( 0.25D+00 * a0 ) ) then

        cj0 = cbj0
        cj1 = cbj1
        do k = 2, n
          cjk = 2.0D+00 * ( k - 1.0D+00 ) / z * cj1 - cj0
          cbj(k) = cjk
          cj0 = cj1
          cj1 = cjk
        end do

      else

        m = msta1 ( a0, 200 )

        if ( m .lt. n ) then
          nm = m
        else
          m = msta2 ( a0, n, 15 )
        end if

        cf2 = cmplx ( 0.0D+00, 0.0D+00 ) 
        cf1 = cmplx ( 1.0D-30, 0.0D+00 )
        do k = m, 0, -1
          cf = 2.0D+00 * ( k + 1.0D+00 ) / z * cf1 - cf2
          if ( k .le. nm ) then
            cbj(k) = cf
          end if
          cf2 = cf1
          cf1 = cf
        end do

        if ( cdabs ( cbj1 ) .lt. cdabs ( cbj0 ) ) then
          cs = cbj0 / cf
        else
          cs = cbj1 / cf2
        end if

        do k = 0, nm
          cbj(k) = cs * cbj(k)
        end do

      end if

      do k = 2, nm
        cdj(k) = cbj(k-1) - k / z * cbj(k)
      end do
      ya0 = cdabs ( cby0 )
      lb = 0
      cg0 = cby0
      cg1 = cby1
      do k = 2, nm
        cyk = 2.0D+00 * ( k - 1.0D+00 ) / z * cg1 - cg0
        if ( cdabs ( cyk ) .le. 1.0D+290 ) then         
          yak = cdabs ( cyk )
          ya1 = cdabs ( cg0 )
          if ( yak .lt. ya0 .and. yak .lt. ya1 ) then
            lb = k
          end if
          cby(k) = cyk
          cg0 = cg1
          cg1 = cyk
        end if
      end do

      if ( lb .le. 4 .or. dimag ( z ) .eq. 0.0D+00 ) then
        go to 20
      end if

10    continue

      if (lb .eq. lb0) then
        go to 20
      end if

      ch2 = cmplx ( 1.0D+00, 0.0D+00 )
      ch1 = cmplx ( 0.0D+00, 0.0D+00 )
      lb0 = lb
      do k = lb, 1, -1
        ch0 = 2.0D+00 * k / z * ch1 - ch2
        ch2 = ch1
        ch1 = ch0
      end do
      cp12 = ch0
      cp22 = ch2
      ch2 = cmplx ( 0.0D+00, 0.0D+00 )
      ch1 = cmplx ( 1.0D+00, 0.0D+00 )
      do k = lb, 1, -1
        ch0 = 2.0D+00 * k / z * ch1 - ch2
        ch2 = ch1
        ch1 = ch0
      end do
      cp11 = ch0
      cp21 = ch2

      if ( lb .eq. nm ) then
        cbj(lb+1) = 2.0D+00 * lb / z * cbj(lb) - cbj(lb-1)
      end if

      if ( cdabs ( cbj(1) ) .lt. cdabs ( cbj(0) ) ) then
        cby(lb+1) = ( cbj(lb+1) * cby0 - 2.0D+00 * cp11 
     &    / ( pi * z ) ) / cbj(0)
        cby(lb) = ( cbj(lb) * cby0 + 2.0D+00 * cp12 / ( pi * z ) ) 
     &    / cbj(0)
      else
        cby(lb+1) = ( cbj(lb+1) * cby1 - 2.0D+00 * cp21 
     &    / ( pi * z ) ) / cbj(1)
        cby(lb) = ( cbj(lb) * cby1 + 2.0D+00 * cp22 / ( pi * z ) ) 
     &    / cbj(1)
      end if

      cyl2 = cby(lb+1)
      cyl1 = cby(lb)
      do k = lb - 1, 0, -1
        cylk = 2.0D+00 * ( k + 1.0D+00 ) / z * cyl1 - cyl2
        cby(k) = cylk
        cyl2 = cyl1
        cyl1 = cylk
      end do

      cyl1 = cby(lb)
      cyl2 = cby(lb+1)
      do k = lb + 1, nm - 1
        cylk = 2.0D+00 * k / z * cyl2 - cyl1
        cby(k+1) = cylk
        cyl1 = cyl2
        cyl2 = cylk
      end do

      do k = 2, nm
        wa = cdabs ( cby(k) )
        if ( wa .lt. cdabs ( cby(k-1) ) ) then
          lb = k
        end if
      end do

      go to 10

20    continue

      do k = 2, nm
        cdy(k) = cby(k-1) - k / z * cby(k)
      end do

      return
      end
      subroutine cjynb ( n, z, nm, cbj, cdj, cby, cdy )

c*********************************************************************72
c
cc CJYNB: Bessel functions, derivatives, Jn(z) and Yn(z) of complex argument.
c
c  Licensing:
c
c    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
c    they give permission to incorporate this routine into a user program 
c    provided that the copyright is acknowledged.
c
c  Modified:
c
c    03 August 2012
c
c  Author:
c
c    Shanjie Zhang, Jianming Jin
c
c  Reference:
c
c    Shanjie Zhang, Jianming Jin,
c    Computation of Special Functions,
c    Wiley, 1996,
c    ISBN: 0-471-11963-6,
c    LC: QA351.C45.
c 
c  Parameters:
c
c    Input, integer N, the order of Jn(z) and Yn(z).
c
c    Input, complex*16 Z, the argument of Jn(z) and Yn(z).
c
c    Output, integer NM, the highest order computed.
c
c    Output, complex*16 CBJ(0:N), CDJ(0:N), CBY(0:N), CDY(0:N), 
c    the values of Jn(z), Jn'(z), Yn(z), Yn'(z).
c
      implicit none

      integer n

      double precision a(4)
      double precision a0
      double precision a1(4)
      double precision b(4)
      double precision b1(4)
      complex*16 cbj(0:n)
      complex*16 cbj0
      complex*16 cbj1
      complex*16 cbjk
      complex*16 cbs
      complex*16 cby(0:n)
      complex*16 cby0
      complex*16 cby1
      complex*16 cdj(0:n)
      complex*16 cdy(0:n)
      complex*16 ce
      complex*16 cf
      complex*16 cf1
      complex*16 cf2
      complex*16 cp0
      complex*16 cp1
      complex*16 cq0
      complex*16 cq1
      complex*16 cs0
      complex*16 csu
      complex*16 csv
      complex*16 ct1
      complex*16 ct2
      complex*16 cu
      complex*16 cyy
      double precision el
      integer k
      integer m
      integer msta1
      integer msta2
      integer nm
      double precision pi
      double precision r2p
      double precision y0
      complex*16 z

      save a
      save a1
      save b
      save b1

      data a    / -0.7031250000000000D-01, 0.1121520996093750D+00,
     &            -0.5725014209747314D+00, 0.6074042001273483D+01/

      data a1    / 0.1171875000000000D+00,-0.1441955566406250D+00,
     &             0.6765925884246826D+00,-0.6883914268109947D+01/

      data b     / 0.7324218750000000D-01,-0.2271080017089844D+00,
     &             0.1727727502584457D+01,-0.2438052969955606D+02/

      data b1    / -0.1025390625000000D+00,0.2775764465332031D+00,
     &             -0.1993531733751297D+01,0.2724882731126854D+02/

      el = 0.5772156649015329D+00
      pi = 3.141592653589793D+00
      r2p = 0.63661977236758D+00
      y0 = abs ( dimag ( z ) )
      a0 = cdabs ( z )
      nm = n

      if ( a0 .lt. 1.0D-100 ) then
        do k = 0, n
          cbj(k) = cmplx ( 0.0D+00, 0.0D+00 )
          cdj(k) = cmplx ( 0.0D+00, 0.0D+00 )
          cby(k) = - cmplx ( 1.0D+30, 0.0D+00 )
          cdy(k) = cmplx ( 1.0D+30, 0.0D+00 )
        end do
        cbj(0) = cmplx ( 1.0D+00, 0.0D+00 )
        cdj(1) = cmplx ( 0.5D+00, 0.0D+00 )
        return
      end if

      if ( a0 .le. 300.0D+00 .or. 80 .lt. n ) then

        if ( n .eq. 0 ) then
          nm = 1
        end if
        m = msta1 ( a0, 200 )
        if ( m .lt. nm ) then
          nm = m
        else
          m = msta2 ( a0, nm, 15 )
        end if

        cbs = cmplx ( 0.0D+00, 0.0D+00 )
        csu = cmplx ( 0.0D+00, 0.0D+00 )
        csv = cmplx ( 0.0D+00, 0.0D+00 )
        cf2 = cmplx ( 0.0D+00, 0.0D+00 )
        cf1 = cmplx ( 1.0D-30, 0.0D+00 )

        do k = m, 0, -1
          cf = 2.0D+00 * ( k + 1.0D+00 ) / z * cf1 - cf2
          if ( k .le. nm ) then
            cbj(k) = cf
          end if
          if ( k .eq. 2 * int ( k / 2 ) .and. k .ne. 0 ) then
            if ( y0 .le. 1.0D+00 ) then
              cbs = cbs + 2.0D+00 * cf
            else
              cbs = cbs + ( -1.0D+00 ) ** ( k / 2 ) * 2.0D+00 * cf
            end if
            csu = csu + ( -1.0D+00 ) ** ( k / 2 ) * cf / k
          else if ( 1 .lt. k ) then
            csv = csv + ( -1.0D+00 ) ** ( k / 2 ) * k 
     &        / ( k * k - 1.0D+00 ) * cf
          end if
          cf2 = cf1
          cf1 = cf
        end do

        if ( y0 .le. 1.0D+00 ) then
          cs0 = cbs + cf
        else
          cs0 = ( cbs + cf ) / cdcos ( z )
        end if

        do k = 0, nm
          cbj(k) = cbj(k) / cs0
        end do

        ce = cdlog ( z / 2.0D+00 ) + el
        cby(0) = r2p * ( ce * cbj(0) - 4.0D+00 * csu / cs0 )
        cby(1) = r2p * ( - cbj(0) / z + ( ce - 1.0D+00 ) * cbj(1) 
     &    - 4.0D+00 * csv / cs0 )

      else

        ct1 = z - 0.25D+00 * pi
        cp0 = cmplx ( 1.0D+00, 0.0D+00 )
        do k = 1, 4
          cp0 = cp0 + a(k) * z ** ( - 2 * k )
        end do
        cq0 = -0.125D+00 / z
        do k = 1, 4
          cq0 = cq0 + b(k) * z ** ( - 2 * k - 1 )
        end do
        cu = cdsqrt ( r2p / z )
        cbj0 = cu * ( cp0 * cdcos ( ct1 ) - cq0 * cdsin ( ct1 ) )
        cby0 = cu * ( cp0 * cdsin ( ct1 ) + cq0 * cdcos ( ct1 ) )
        cbj(0) = cbj0
        cby(0) = cby0
        ct2 = z - 0.75D+00 * pi
        cp1 = cmplx ( 1.0D+00, 0.0D+00 )
        do k = 1, 4
          cp1 = cp1 + a1(k) * z ** ( - 2 * k )
        end do
        cq1 = 0.375D+00 / z
        do k = 1, 4
          cq1 = cq1 + b1(k) * z ** ( - 2 * k - 1 )
        end do
        cbj1 = cu * ( cp1 * cdcos ( ct2 ) - cq1 * cdsin ( ct2 ) )
        cby1 = cu * ( cp1 * cdsin ( ct2 ) + cq1 * cdcos ( ct2 ) )
        cbj(1) = cbj1
        cby(1) = cby1
        do k = 2, nm
          cbjk = 2.0D+00 * ( k - 1.0D+00 ) / z * cbj1 - cbj0
          cbj(k) = cbjk
          cbj0 = cbj1
          cbj1 = cbjk
        end do
      end if

      cdj(0) = -cbj(1)
      do k = 1, nm
        cdj(k) = cbj(k-1) - k / z * cbj(k)
      end do

      if ( 1.0D+00 .lt. cdabs ( cbj(0) ) ) then
        cby(1) = ( cbj(1) * cby(0) - 2.0D+00 / ( pi * z ) ) / cbj(0)
      end if

      do k = 2, nm
        if ( cdabs ( cbj(k-2) ) .le. cdabs ( cbj(k-1) ) ) then
          cyy = ( cbj(k) * cby(k-1) - 2.0D+00 / ( pi * z ) ) / cbj(k-1)
        else
          cyy = ( cbj(k) * cby(k-2) - 4.0D+00 * ( k - 1.0D+00 )
     &      / ( pi * z * z ) ) / cbj(k-2)
        end if
        cby(k) = cyy
      end do

      cdy(0) = -cby(1)
      do k = 1, nm
        cdy(k) = cby(k-1) - k / z * cby(k)
      end do

      return
      end
      subroutine cjyva ( v, z, vm, cbj, cdj, cby, cdy )

c*********************************************************************72
c
cc CJYVA: Bessel functions and derivatives, Jv(z) and Yv(z) of complex argument.
c
c  Licensing:
c
c    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
c    they give permission to incorporate this routine into a user program 
c    provided that the copyright is acknowledged.
c
c  Modified:
c
c    03 August 2012
c
c  Author:
c
c    Shanjie Zhang, Jianming Jin
c
c  Reference:
c
c    Shanjie Zhang, Jianming Jin,
c    Computation of Special Functions,
c    Wiley, 1996,
c    ISBN: 0-471-11963-6,
c    LC: QA351.C45.
c
c  Parameters:
c
c    Input, double precision V, the order of Jv(z) and Yv(z).
c
c    Input, complex*16, the argument.
c
c    Output, double precision VM, the highest order computed.
c
c    Output, double precision CBJ(0:*), CDJ(0:*), CBY(0:*), CDY(0:*), 
c    the values of Jn+v0(z), Jn+v0'(z), Yn+v0(z), Yn+v0'(z).
c
      implicit none

      double precision a0
      complex*16 ca
      complex*16 ca0
      complex*16 cb
      complex*16 cbj(0:*)
      complex*16 cby(0:*)
      complex*16 cck
      complex*16 cdj(0:*)
      complex*16 cdy(0:*)
      complex*16 cec
      complex*16 cf
      complex*16 cf0
      complex*16 cf1
      complex*16 cf2
      complex*16 cfac0
      complex*16 cfac1
      complex*16 cg0
      complex*16 cg1
      complex*16 ch0
      complex*16 ch1
      complex*16 ch2
      complex*16 ci
      complex*16 cju0
      complex*16 cju1
      complex*16 cjv0
      complex*16 cjv1
      complex*16 cjvl
      complex*16 cp11
      complex*16 cp12
      complex*16 cp21
      complex*16 cp22
      complex*16 cpz
      complex*16 cqz
      complex*16 cr
      complex*16 cr0
      complex*16 cr1
      complex*16 crp
      complex*16 crq
      complex*16 cs
      complex*16 cs0
      complex*16 cs1
      complex*16 csk
      complex*16 cyk
      complex*16 cyl1
      complex*16 cyl2
      complex*16 cylk
      complex*16 cyv0
      complex*16 cyv1
      double precision ga
      double precision gb
      integer j
      integer k
      integer k0
      integer l
      integer lb
      integer lb0
      integer m
      integer msta1
      integer msta2
      integer n
      double precision pi
      double precision pv0
      double precision pv1
      double precision rp2
      double precision v
      double precision v0
      double precision vg
      double precision vl
      double precision vm
      double precision vv
      double precision w0
      double precision w1
      double precision wa
      double precision ya0
      double precision ya1
      double precision yak
      complex*16 z
      complex*16 z1
      complex*16 z2
      complex*16 zk

      pi = 3.141592653589793D+00
      rp2 = 0.63661977236758D+00
      ci = cmplx ( 0.0D+00, 1.0D+00 )
      a0 = cdabs ( z )
      z1 = z
      z2 = z * z
      n = int ( v )
      v0 = v - n
      pv0 = pi * v0
      pv1 = pi * ( 1.0D+00 + v0 )

      if ( a0 .lt. 1.0D-100 ) then

        do k = 0, n
          cbj(k) = cmplx ( 0.0D+00, 0.0D+00 )
          cdj(k) = cmplx ( 0.0D+00, 0.0D+00 )
          cby(k) = - cmplx ( 1.0D+30, 0.0D+00 )
          cdy(k) = cmplx ( 1.0D+30, 0.0D+00 )
        end do

        if ( v0 .eq. 0.0D+00 ) then
          cbj(0) = cmplx ( 1.0D+00, 0.0D+00 )
          cdj(1) = cmplx ( 0.5D+00, 0.0D+00 )
        else
          cdj(0) = cmplx ( 1.0D+30, 0.0D+00 )
        end if

        vm = v                     
        return

      end if

      if ( real ( z ) .lt. 0.0D+00 ) then
        z1 = -z
      end if

      if ( a0 .le. 12.0D+00 ) then

        do l = 0, 1
          vl = v0 + l
          cjvl = cmplx ( 1.0D+00, 0.0D+00 )
          cr = cmplx ( 1.0D+00, 0.0D+00 )
          do k = 1, 40
            cr = -0.25D+00 * cr * z2 / ( k * ( k + vl ) )
            cjvl = cjvl + cr
            if ( cdabs ( cr ) .lt. cdabs ( cjvl ) * 1.0D-15 ) then
              go to 10
            end if
          end do

10        continue

          vg = 1.0D+00 + vl
          call gamma ( vg, ga )
          ca = ( 0.5D+00 * z1 ) ** vl / ga

          if ( l .eq. 0 ) then
            cjv0 = cjvl * ca
          else
            cjv1 = cjvl * ca
          end if

        end do

      else

        if ( a0 .lt. 35.0D+00 ) then
          k0 = 11
        else if ( a0 .lt.50.0D+00 ) then
          k0 = 10
        else
          k0 = 8
        end if

        do j = 0, 1
          vv = 4.0D+00 * ( j + v0 ) * ( j + v0 )
          cpz = cmplx ( 1.0D+00, 0.0D+00 )
          crp = cmplx ( 1.0D+00, 0.0D+00 )
          do k = 1, k0
            crp = - 0.78125D-02 * crp 
     &        * ( vv - ( 4.0D+00 * k - 3.0D+00 ) ** 2 )
     &        * ( vv - ( 4.0D+00 * k - 1.0D+00 ) ** 2 ) 
     &        / ( k * ( 2.0D+00 * k - 1.0D+00 ) * z2 )
            cpz = cpz + crp
          end do
          cqz = cmplx ( 1.0D+00, 0.0D+00 )
          crq = cmplx ( 1.0D+00, 0.0D+00 )
          do k = 1, k0
            crq = -0.78125D-02 * crq 
     &        * ( vv - ( 4.0D+00 * k - 1.0D+00 ) ** 2 ) 
     &        * ( vv - ( 4.0D+00 * k + 1.0D+00 ) ** 2 )
     &        / ( k * ( 2.0D+00 * k + 1.0D+00 ) * z2 )
            cqz = cqz + crq
          end do
          cqz = 0.125D+00 * ( vv - 1.0D+00 ) * cqz / z1
          zk = z1 - ( 0.5D+00 * ( j + v0 ) + 0.25D+00 ) * pi
          ca0 = cdsqrt ( rp2 / z1 )
          cck = cdcos ( zk )
          csk = cdsin ( zk )
          if ( j .eq. 0 ) then
            cjv0 = ca0 * ( cpz * cck - cqz * csk )
            cyv0 = ca0 * ( cpz * csk + cqz * cck )
          else if ( j .eq. 1 ) then
            cjv1 = ca0 * ( cpz * cck - cqz * csk )
            cyv1 = ca0 * ( cpz * csk + cqz * cck )
          end if
        end do

      end if

      if ( a0 .le. 12.0D+00 ) then

        if ( v0 .ne. 0.0D+00 ) then

          do l = 0, 1
            vl = v0 + l
            cjvl = cmplx ( 1.0D+00, 0.0D+00 )
            cr = cmplx ( 1.0D+00, 0.0D+00 )
            do k = 1, 40
              cr = -0.25D+00 * cr * z2 / ( k * ( k - vl ) )
              cjvl = cjvl + cr
              if ( cdabs ( cr ) .lt. cdabs ( cjvl ) * 1.0D-15 ) then
                go to 20
              end if
            end do

20          continue

            vg = 1.0D+00 - vl
            call gamma ( vg, gb )
            cb = ( 2.0D+00 / z1 ) ** vl / gb
            if ( l .eq. 0 ) then
              cju0 = cjvl * cb
            else
              cju1 = cjvl * cb
            end if
          end do
          cyv0 = ( cjv0 * cos ( pv0 ) - cju0 ) / sin ( pv0 )
          cyv1 = ( cjv1 * cos ( pv1 ) - cju1 ) / sin ( pv1 )

        else

          cec = cdlog ( z1 / 2.0D+00 ) + 0.5772156649015329D+00
          cs0 = cmplx ( 0.0D+00, 0.0D+00)
          w0 = 0.0D+00
          cr0 = cmplx ( 1.0D+00, 0.0D+00 )
          do k = 1, 30
            w0 = w0 + 1.0D+00 / k
            cr0 = -0.25D+00 * cr0 / ( k * k ) * z2
            cs0 = cs0 + cr0 * w0
          end do
          cyv0 = rp2 * ( cec * cjv0 - cs0 )
          cs1 = cmplx ( 1.0D+00, 0.0D+00 )
          w1 = 0.0D+00
          cr1 = cmplx ( 1.0D+00, 0.0D+00 )
          do k = 1, 30
            w1 = w1 + 1.0D+00 / k
            cr1 = -0.25D+00 * cr1 / ( k * ( k + 1 ) ) * z2
            cs1 = cs1 + cr1 * ( 2.0D+00 * w1 + 1.0D+00 
     &        / ( k + 1.0D+00 ) )
          end do
          cyv1 = rp2 * ( cec * cjv1 - 1.0D+00 / z1 
     &      - 0.25D+00 * z1 * cs1 )

        end if

      end if

      if ( real ( z ) .lt. 0.0D+00 ) then

        cfac0 = cdexp ( pv0 * ci )
        cfac1 = cdexp ( pv1 * ci )

        if ( dimag ( z ) .lt. 0.0D+00 ) then
          cyv0 = cfac0 * cyv0 - 2.0D+00 * ci * cos ( pv0 ) * cjv0
          cyv1 = cfac1 * cyv1 - 2.0D+00 * ci * cos ( pv1 ) * cjv1
          cjv0 = cjv0 / cfac0
          cjv1 = cjv1 / cfac1
        else if ( 0.0D+00 .lt. dimag ( z ) ) then
          cyv0 = cyv0 / cfac0 + 2.0D+00 * ci * cos ( pv0 ) * cjv0
          cyv1 = cyv1 / cfac1 + 2.0D+00 * ci * cos ( pv1 ) * cjv1
          cjv0 = cfac0 * cjv0
          cjv1 = cfac1 * cjv1
        end if

      end if

      cbj(0) = cjv0
      cbj(1) = cjv1

      if ( 2 .le. n .and. n .le. int ( 0.25D+00 * a0 ) ) then

        cf0 = cjv0
        cf1 = cjv1
        do k = 2, n
          cf = 2.0D+00 * ( k + v0 - 1.0D+00 ) / z * cf1 - cf0
          cbj(k) = cf
          cf0 = cf1
          cf1 = cf
        end do

      else if ( 2 .le. n ) then

        m = msta1 ( a0, 200 )
        if ( m .lt. n ) then
          n = m
        else
          m = msta2 ( a0, n, 15 )
        end if
        cf2 = cmplx ( 0.0D+00, 0.0D+00 )
        cf1 = cmplx ( 1.0D-30, 0.0D+00 )
        do k = m, 0, -1
          cf = 2.0D+00 * ( v0 + k + 1.0D+00 ) / z * cf1 - cf2
          if ( k .le. n ) then
            cbj(k) = cf
          end if
          cf2 = cf1
          cf1 = cf
        end do
        if ( cdabs ( cjv1 ) .lt. cdabs ( cjv0 ) ) then
          cs = cjv0 / cf
        else
          cs = cjv1 / cf2
        end if

        do k = 0, n
          cbj(k) = cs * cbj(k)
        end do

      end if

        cdj(0) = v0 / z * cbj(0) - cbj(1)
        do k = 1, n
          cdj(k) = - ( k + v0 ) / z * cbj(k) + cbj(k-1)
        end do

        cby(0) = cyv0
        cby(1) = cyv1
        ya0 = cdabs ( cyv0 )
        lb = 0
        cg0 = cyv0
        cg1 = cyv1
        do k = 2, n
          cyk = 2.0D+00 * ( v0 + k - 1.0D+00 ) / z * cg1 - cg0
          if ( cdabs ( cyk ) .le. 1.0D+290 ) then
            yak = cdabs ( cyk )
            ya1 = cdabs ( cg0 )
            if ( yak .lt. ya0 .and. yak .lt. ya1 ) then
              lb = k
            end if
            cby(k) = cyk
            cg0 = cg1
            cg1 = cyk
          end if
        end do

        if ( lb .le. 4 .or. dimag ( z ) .eq. 0.0D+00 ) then
          go to 40
        end if

30      continue

        if ( lb .eq. lb0 ) then
          go to 40
        end if

        ch2 = cmplx ( 1.0D+00, 0.0D+00 )
        ch1 = cmplx ( 0.0D+00, 0.0D+00 )
        lb0 = lb
        do k = lb, 1, -1
          ch0 = 2.0D+00 * ( k + v0 ) / z * ch1 - ch2
          ch2 = ch1
          ch1 = ch0
        end do
        cp12 = ch0
        cp22 = ch2
        ch2 = cmplx ( 0.0D+00, 0.0D+00 )
        ch1 = cmplx ( 1.0D+00, 0.0D+00 )
        do k = lb, 1, -1
          ch0 = 2.0D+00 * ( k + v0 ) / z * ch1 - ch2
          ch2 = ch1
          ch1 = ch0
        end do
        cp11 = ch0
        cp21 = ch2

        if ( lb .eq. n ) then
          cbj(lb+1) = 2.0D+00 * ( lb + v0 ) / z * cbj(lb) - cbj(lb-1)
        end if

        if ( cdabs ( cbj(1) ) .lt. cdabs ( cbj(0) ) ) then
          cby(lb+1) = ( cbj(lb+1) * cyv0 - 2.0D+00 * cp11 
     &      / ( pi * z ) ) / cbj(0)
          cby(lb) = ( cbj(lb) * cyv0 + 2.0D+00 * cp12 
     &      / ( pi * z ) ) / cbj(0)
        else
          cby(lb+1) = ( cbj(lb+1) * cyv1 - 2.0D+00 * cp21 
     &      / ( pi * z ) ) / cbj(1)
          cby(lb) = ( cbj(lb) * cyv1 + 2.0D+00 * cp22 
     &      / ( pi * z ) ) / cbj(1)
        end if

        cyl2 = cby(lb+1)
        cyl1 = cby(lb)
        do k = lb - 1, 0, -1
          cylk = 2.0D+00 * ( k + v0 + 1.0D+00 ) / z * cyl1 - cyl2
          cby(k) = cylk
          cyl2 = cyl1
          cyl1 = cylk
        end do

      cyl1 = cby(lb)
      cyl2 = cby(lb+1)
      do k = lb + 1, n - 1
        cylk = 2.0D+00 * ( k + v0 ) / z * cyl2 - cyl1
        cby(k+1) = cylk
        cyl1 = cyl2
        cyl2 = cylk
      end do

      do k = 2, n
        wa = cdabs ( cby(k) )
        if ( wa .lt. cdabs ( cby(k-1) ) ) then
          lb = k
        end if
      end do

      go to 30

40    continue

      cdy(0) = v0 / z * cby(0) - cby(1)
      do k = 1, n
        cdy(k) = cby(k-1) - ( k + v0 ) / z * cby(k)
      end do
      vm = n + v0

      return
      end
      subroutine cjyvb ( v, z, vm, cbj, cdj, cby, cdy )

c*********************************************************************72
c
cc CJYVB: Bessel functions and derivatives, Jv(z) and Yv(z) of complex argument.
c
c  Licensing:
c
c    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
c    they give permission to incorporate this routine into a user program 
c    provided that the copyright is acknowledged.
c
c  Modified:
c
c    03 August 2012
c
c  Author:
c
c    Shanjie Zhang, Jianming Jin
c
c  Reference:
c
c    Shanjie Zhang, Jianming Jin,
c    Computation of Special Functions,
c    Wiley, 1996,
c    ISBN: 0-471-11963-6,
c    LC: QA351.C45.
c
c  Parameters:
c
c    Input, double precision V, the order of Jv(z) and Yv(z).
c
c    Input, complex*16 Z, the argument.
c
c    Output, double precision VM, the highest order computed.
c
c    Output, double precision CBJ(0:*), CDJ(0:*), CBY(0:*), CDY(0:*), 
c    the values of Jn+v0(z), Jn+v0'(z), Yn+v0(z), Yn+v0'(z).
c
      implicit none

      double precision a0
      complex*16 ca
      complex*16 ca0
      complex*16 cb
      complex*16 cbj(0:*)
      complex*16 cby(0:*)
      complex*16 cck
      complex*16 cdj(0:*)
      complex*16 cdy(0:*)
      complex*16 cec
      complex*16 cf
      complex*16 cf1
      complex*16 cf2
      complex*16 cfac0
      complex*16 ci
      complex*16 cju0
      complex*16 cjv0
      complex*16 cjvn
      complex*16 cpz
      complex*16 cqz
      complex*16 cr
      complex*16 cr0
      complex*16 crp
      complex*16 crq
      complex*16 cs
      complex*16 cs0
      complex*16 csk
      complex*16 cyv0
      complex*16 cyy
      double precision ga
      double precision gb
      integer k
      integer k0
      integer m
      integer msta1
      integer msta2
      integer n
      double precision pi
      double precision pv0
      double precision rp2
      double precision v
      double precision v0
      double precision vg
      double precision vm
      double precision vv
      double precision w0
      complex*16 z
      complex*16 z1
      complex*16 z2
      complex*16 zk

      pi = 3.141592653589793D+00
      rp2 = 0.63661977236758D+00
      ci = cmplx ( 0.0D+00, 1.0D+00 )
      a0 = cdabs ( z )
      z1 = z
      z2 = z * z
      n = int ( v )
      v0 = v - n
      pv0 = pi * v0
  
      if ( a0 .lt. 1.0D-100 ) then

        do k = 0, n
          cbj(k) = cmplx ( 0.0D+00, 0.0D+00 )
          cdj(k) = cmplx ( 0.0D+00, 0.0D+00 )
          cby(k) = - cmplx ( 1.0D+30, 0.0D+00 )
          cdy(k) = cmplx ( 1.0D+30, 0.0D+00 )
        end do

        if ( v0 .eq. 0.0D+00 ) then
          cbj(0) = cmplx ( 1.0D+00, 0.0D+00 )
          cdj(1) = cmplx ( 0.5D+00, 0.0D+00 )
        else
          cdj(0) = cmplx ( 1.0D+30, 0.0D+00 )
        end if

        vm = v
        return

      end if

      if ( real ( z ) .lt. 0.0D+00 ) then
        z1 = -z
      end if

      if ( a0 .le. 12.0D+00 ) then

        cjv0 = cmplx ( 1.0D+00, 0.0D+00 )
        cr = cmplx ( 1.0D+00, 0.0D+00 )
        do k = 1, 40
          cr = -0.25D+00 * cr * z2 / ( k * ( k + v0 ) )
          cjv0 = cjv0 + cr
          if ( cdabs ( cr ) .lt. cdabs ( cjv0 ) * 1.0D-15 ) then
            go to 10
          end if
        end do

10      continue

        vg = 1.0D+00 + v0
        call gamma ( vg, ga )
        ca = ( 0.5D+00 * z1 ) ** v0 / ga
        cjv0 = cjv0 * ca

      else

        if ( a0 .lt. 35.0D+00 ) then
          k0 = 11
        else if ( a0 .lt. 50.0D+00 ) then
          k0 = 10
        else
          k0 = 8
        end if

        vv = 4.0D+00 * v0 * v0
        cpz = cmplx ( 1.0D+00, 0.0D+00 )
        crp = cmplx ( 1.0D+00, 0.0D+00 )
        do k = 1, k0
          crp = -0.78125D-02 * crp 
     &      * ( vv - ( 4.0D+00 * k - 3.0D+00 ) ** 2 )
     &      * ( vv - ( 4.0D+00 * k - 1.0D+00 ) **2 ) 
     &      / ( k * ( 2.0D+00 * k - 1.0D+00 ) * z2 )
          cpz = cpz + crp
        end do
        cqz = cmplx ( 1.0D+00, 0.0D+00 )
        crq = cmplx ( 1.0D+00, 0.0D+00 )
        do k = 1, k0
          crq = -0.78125D-02 * crq
     &      * ( vv - ( 4.0D+00 * k - 1.0D+00 ) ** 2 )
     &      * ( vv - ( 4.0D+00 * k + 1.0D+00 ) ** 2 ) 
     &      / ( k * ( 2.0D+00 * k + 1.0D+00 ) * z2 )
          cqz = cqz + crq
        end do
        cqz = 0.125D+00 * ( vv - 1.0D+00 ) * cqz / z1
        zk = z1 - ( 0.5D+00 * v0 + 0.25D+00 ) * pi
        ca0 = cdsqrt ( rp2 / z1 )
        cck = cdcos ( zk )
        csk = cdsin ( zk )
        cjv0 = ca0 * ( cpz * cck - cqz * csk )
        cyv0 = ca0 * ( cpz * csk + cqz * cck )

      end if

      if ( a0 .le. 12.0D+00 ) then

        if ( v0 .ne. 0.0D+00 ) then

          cjvn = cmplx ( 1.0D+00, 0.0D+00 )
          cr = cmplx ( 1.0D+00, 0.0D+00 )
          do k = 1, 40
            cr = -0.25D+00 * cr * z2 / ( k * ( k - v0 ) )
            cjvn = cjvn + cr
            if ( cdabs ( cr ) .lt. cdabs ( cjvn ) * 1.0D-15 ) then
              go to 20
            end if
          end do

20        continue

          vg = 1.0D+00 - v0
          call gamma ( vg, gb )
          cb = ( 2.0D+00 / z1 ) ** v0 / gb
          cju0 = cjvn * cb
          cyv0 = ( cjv0 * cos ( pv0 ) - cju0 ) / sin ( pv0 )

        else

          cec = cdlog ( z1 / 2.0D+00 ) + 0.5772156649015329D+00
          cs0 = cmplx ( 0.0D+00, 0.0D+00 )
          w0 = 0.0D+00
          cr0 = cmplx ( 1.0D+00, 0.0D+00 )
          do k = 1, 30
            w0 = w0 + 1.0D+00 / k
            cr0 = -0.25D+00 * cr0 / ( k * k ) * z2
            cs0 = cs0 + cr0 * w0
          end do
          cyv0 = rp2 * ( cec * cjv0 - cs0 )

        end if

      end if

      if ( n .eq. 0 ) then
        n = 1
      end if

      m = msta1 ( a0, 200 )
      if ( m .lt. n ) then
        n = m
      else
        m = msta2 ( a0, n, 15 )
      end if

      cf2 = cmplx ( 0.0D+00, 0.0D+00 )
      cf1 = cmplx ( 1.0D-30, 0.0D+00 )
      do k = m, 0, -1
        cf = 2.0D+00 * ( v0 + k + 1.0D+00 ) / z1 * cf1 - cf2
        if ( k .le. n ) then
          cbj(k) = cf
        end if
        cf2 = cf1
        cf1 = cf
      end do

      cs = cjv0 / cf
      do k = 0, n
        cbj(k) = cs * cbj(k)
      end do

      if ( real ( z ) .lt. 0.0D+00) then

        cfac0 = cdexp ( pv0 * ci )
        if ( dimag ( z ) .lt. 0.0D+00 ) then
          cyv0 = cfac0 * cyv0 - 2.0D+00 * ci * cos ( pv0 ) * cjv0
        else if ( 0.0D+00 .lt. dimag ( z ) ) then
          cyv0 = cyv0 / cfac0 + 2.0D+00 * ci * cos ( pv0 ) * cjv0
        end if

        do k = 0, n
          if ( dimag ( z ) .lt. 0.0D+00) then
            cbj(k) = cdexp ( - pi * ( k + v0 ) * ci ) * cbj(k)
          else if ( 0.0D+00 .lt. dimag ( z ) ) then
            cbj(k) = cdexp ( pi * ( k + v0 ) * ci ) * cbj(k)
          end if
        end do

        z1 = z1

      end if

      cby(0) = cyv0
      do k = 1, n
        cyy = ( cbj(k) * cby(k-1) - 2.0D+00 / ( pi * z ) ) / cbj(k-1)
        cby(k) = cyy
      end do

      cdj(0) = v0 / z * cbj(0) - cbj(1)
      do k = 1, n
        cdj(k) = - ( k + v0 ) / z * cbj(k) + cbj(k-1)
      end do

      cdy(0) = v0 / z * cby(0) - cby(1)
      do k = 1, n
        cdy(k) = cby(k-1) - ( k + v0 ) / z * cby(k)
      end do

      vm = n + v0

      return
      end
      subroutine clpmn ( mm, m, n, x, y, cpm, cpd )

c*********************************************************************72
c
cc CLPMN: associated Legendre functions and derivatives for complex argument.
c
c  Discussion:
c
c    Compute the associated Legendre functions Pmn(z)   
c    and their derivatives Pmn'(z) for a complex argument
c
c  Licensing:
c
c    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
c    they give permission to incorporate this routine into a user program 
c    provided that the copyright is acknowledged.
c
c  Modified:
c
c    01 August 2012
c
c  Author:
c
c    Shanjie Zhang, Jianming Jin
c
c  Reference:
c
c    Shanjie Zhang, Jianming Jin,
c    Computation of Special Functions,
c    Wiley, 1996,
c    ISBN: 0-471-11963-6,
c    LC: QA351.C45.
c
c  Parameters:
c
c    Input, integer MM, the physical dimension of CPM and CPD.
c
c    Input, integer M, N, the order and degree of Pmn(z).
c
c    Input, double precision X, Y, the real and imaginary parts of the argument Z.
c
c    Output, complex*16 CPM(0:MM,0:N), CPD(0:MM,0:N), the values of
c    Pmn(z) and Pmn'(z).
c
      implicit none

      integer mm

      complex*16 cpd(0:mm,0:n)
      complex*16 cpm(0:mm,0:n)
      integer i
      integer j
      integer ls
      integer m
      integer n
      double precision x
      double precision y
      complex*16 z
      complex*16 zq
      complex*16 zs

      z = cmplx ( x, y )

      do i = 0, n
        do j = 0, m
          cpm(j,i) = cmplx ( 0.0D+00, 0.0D+00 )
          cpd(j,i) = cmplx ( 0.0D+00, 0.0D+00 )
        end do
      end do

      cpm(0,0) = cmplx ( 1.0D+00, 0.0D+00 )

      if ( abs ( x ) .eq. 1.0D+00 .and. y .eq. 0.0D+00 ) then

        do i = 1, n
          cpm(0,i) = x ** i
          cpd(0,i) = 0.5D+00 * i * ( i + 1 ) * x ** ( i + 1 )
        end do

        do j = 1, n
          do i = 1, m
            if ( i .eq. 1 ) then
              cpd(i,j) = cmplx ( 1.0D+30, 0.0D+00 )
            else if ( i .eq. 2 ) then
              cpd(i,j) = -0.25D+00 
     &          * ( j + 2 ) * ( j + 1 ) * j * ( j - 1 ) * x ** ( j + 1 )
            end if
          end do
        end do

        return

      end if

      if ( 1.0D+00 .lt. cdabs ( z ) ) then
        ls = -1
      else
        ls = 1
      end if

      zq = cdsqrt ( ls * ( 1.0D+00 - z * z ) )
      zs = ls * ( 1.0D+00 - z * z )
      do i = 1, m
        cpm(i,i) = -ls * ( 2.0D+00 * i - 1.0D+00 ) * zq * cpm(i-1,i-1)
      end do
      do i = 0, m
        cpm(i,i+1) = ( 2.0D+00 * i + 1.0D+00 ) * z * cpm(i,i)
      end do

      do i = 0, m
        do j = i + 2, n
          cpm(i,j) = ( ( 2.0D+00 * j - 1.0D+00 ) * z * cpm(i,j-1)
     &      - ( i + j - 1.0D+00 ) * cpm(i,j-2) ) / ( j - i )
        end do
      end do

      cpd(0,0) = cmplx ( 0.0D+00, 0.0D+00 )
      do j = 1, n
        cpd(0,j) = ls * j * ( cpm(0,j-1) - z * cpm(0,j) ) / zs
      end do 

      do i = 1, m
        do j = i, n
          cpd(i,j) = ls * i * z * cpm(i,j) / zs 
     &      + ( j + i ) * ( j - i + 1.0D+00 ) / zq * cpm(i-1,j)
        end do
      end do

      return
      end
      subroutine clpn ( n, x, y, cpn, cpd )

c*********************************************************************72
c
cc CLPN computes Legendre functions and derivatives for complex argument.
c
c  Discussion:
c
c    Compute Legendre polynomials Pn(z) and their derivatives Pn'(z) for a complex
c    argument.
c
c  Licensing:
c
c    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
c    they give permission to incorporate this routine into a user program 
c    provided that the copyright is acknowledged.
c
c  Modified:
c
c    15 July 2012
c
c  Author:
c
c    Shanjie Zhang, Jianming Jin
c
c  Reference:
c
c    Shanjie Zhang, Jianming Jin,
c    Computation of Special Functions,
c    Wiley, 1996,
c    ISBN: 0-471-11963-6,
c    LC: QA351.C45.
c
c  Parameters:
c
c    Input, integer N, the degree.
c
c    Input, double precision X, Y, the real and imaginary parts of the argument.
c
c    Output, complex*16 CPN(0:N), CPD(0:N), the values of Pn(z), Pn'(z).
c
      implicit none

      integer n

      complex*16 cp0
      complex*16 cp1
      complex*16 cpd(0:n)
      complex*16 cpf
      complex*16 cpn(0:n)
      integer k
      double precision x
      double precision y
      complex*16 z

      z = cmplx ( x, y )

      cpn(0) = cmplx ( 1.0D+00, 0.0D+00 )
      cpn(1) = z
      cpd(0) = cmplx ( 0.0D+00, 0.0D+00 )
      cpd(1) = cmplx ( 1.0D+00, 0.0D+00 )

      cp0 = cmplx ( 1.0D+00, 0.0D+00 )
      cp1 = z
      do k = 2, n
        cpf = ( 2.0D+00 * k - 1.0D+00 ) / k * z * cp1 
     &    - ( k - 1.0D+00 ) / k * cp0
        cpn(k) = cpf
        if ( abs ( x ) .eq. 1.0D+00 .and. y .eq. 0.0D+00 ) then
          cpd(k) = 0.5D+00 * x ** ( k + 1 ) * k * ( k + 1.0D+00 )
        else
          cpd(k) = k * ( cp1 - z * cpf ) / ( 1.0D+00 - z * z )
        end if
        cp0 = cp1
        cp1 = cpf
      end do

      return
      end
      subroutine clqmn ( mm, m, n, x, y, cqm, cqd )

c*********************************************************************72
c
cc CLQMN computes associated Legendre functions and derivatives for complex argument.
c
c  Discussion:
c
c    This procedure computes the associated Legendre functions of the second 
c    kind, Qmn(z) and Qmn'(z), for a complex argument.
c
c  Licensing:
c
c    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
c    they give permission to incorporate this routine into a user program 
c    provided that the copyright is acknowledged.
c
c  Modified:
c
c    02 August 2012
c
c  Author:
c
c    Shanjie Zhang, Jianming Jin
c
c  Reference:
c
c    Shanjie Zhang, Jianming Jin,
c    Computation of Special Functions,
c    Wiley, 1996,
c    ISBN: 0-471-11963-6,
c    LC: QA351.C45.
c 
c  Parameters:
c
c    Input, integer MM, the physical dimension of CQM and CQD.
c
c    Input, integer M, N, the order and degree of Qmn(z).
c
c    Input, double precision X, Y, the real and imaginary parts of the 
c    argument Z.
c
c    Output, complex*16 CQM(0:MM,0:N), CQD(0:MM,0:N), the values of
c    Qmn(z) and Qmn'(z).
c
      implicit none

      integer mm
      integer n 

      complex*16 cq0
      complex*16 cq1
      complex*16 cq10
      complex*16 cqf
      complex*16 cqf0
      complex*16 cqf1
      complex*16 cqf2
      complex*16 cqm(0:mm,0:n)
      complex*16 cqd(0:mm,0:n)
      integer i
      integer j
      integer k
      integer km
      integer ls
      integer m
      double precision x
      double precision xc
      double precision y
      complex*16 z
      complex*16 zq
      complex*16 zs

      z = cmplx ( x, y )

      if ( abs ( x ) .eq. 1.0D+00 .and. y .eq. 0.0D+00 ) then
        do i = 0, m
          do j = 0, n
            cqm(i,j) = cmplx ( 1.0D+30, 0.0D+00 )
            cqd(i,j) = cmplx ( 1.0D+30, 0.0D+00 )
          end do
        end do
        return
      end if

      xc = cdabs ( z )

      if ( dimag ( z ) .eq. 0.0D+00 .or. xc .lt. 1.0D+00 ) then
        ls = 1
      end if

      if ( 1.0D+00 .lt. xc ) then
        ls = -1
      end if

      zq = cdsqrt ( ls * ( 1.0D+00 - z * z ) )
      zs = ls * ( 1.0D+00 - z * z )
      cq0 = 0.5D+00 * cdlog ( ls * ( 1.0D+00 + z ) / ( 1.0D+00 - z ) )

      if ( xc .lt. 1.0001D+00 ) then

        cqm(0,0) = cq0
        cqm(0,1) = z * cq0 - 1.0D+00
        cqm(1,0) = -1.0D+00 / zq
        cqm(1,1) = - zq * ( cq0 + z / ( 1.0D+00 - z * z ) )
        do i = 0, 1
          do j = 2, n
            cqm(i,j) = ( ( 2.0D+00 * j - 1.0D+00 ) * z * cqm(i,j-1)
     &        - ( j + i - 1.0D+00 ) * cqm(i,j-2) ) / ( j - i )
          end do
        end do

        do j = 0, n
          do i = 2, m
            cqm(i,j) = -2.0D+00 * ( i - 1.0D+00 ) * z / zq * cqm(i-1,j)
     &        - ls * ( j + i - 1.0D+00 ) 
     &        * ( j - i + 2.0D+00 ) * cqm(i-2,j)
          end do
        end do

      else

        if ( 1.1D+00 .lt. xc ) then
          km = 40 + m + n
        else
          km = ( 40 + m + n ) 
     &      * int ( - 1.0D+00 - 1.8D+00 * log ( xc - 1.0D+00 ) )
        end if

        cqf2 = cmplx ( 0.0D+00, 0.0D+00 )
        cqf1 = cmplx ( 1.0D+00, 0.0D+00 )
        do k = km, 0, -1
          cqf0 = ( ( 2 * k + 3.0D+00 ) * z * cqf1
     &      - ( k + 2.0D+00 ) * cqf2 ) / ( k + 1.0D+00 )
          if ( k .le. n ) then
            cqm(0,k) = cqf0
          end if
          cqf2 = cqf1
          cqf1 = cqf0
        end do

        do k = 0, n
          cqm(0,k) = cq0 * cqm(0,k) / cqf0
        end do
 
        cqf2 = 0.0D+00
        cqf1 = 1.0D+00
        do k = km, 0, -1
          cqf0 = ( ( 2 * k + 3.0D+00 ) * z * cqf1
     &      - ( k + 1.0D+00 ) * cqf2 ) / ( k + 2.0D+00 )
          if ( k .le. n ) then
            cqm(1,k) = cqf0
          end if
          cqf2 = cqf1
          cqf1 = cqf0
        end do

        cq10 = -1.0D+00 / zq
        do k = 0, n 
          cqm(1,k) = cq10 * cqm(1,k) / cqf0
        end do

        do j = 0, n
          cq0 = cqm(0,j)
          cq1 = cqm(1,j)
          do i = 0, m - 2
            cqf = -2.0D+00 * ( i + 1 ) * z / zq * cq1
     &        + ( j - i ) * ( j + i + 1.0D+00 ) * cq0
            cqm(i+2,j) = cqf
            cq0 = cq1
            cq1 = cqf
          end do
        end do

      end if

      cqd(0,0) = ls / zs
      do j = 1, n
        cqd(0,j) = ls * j * ( cqm(0,j-1) - z * cqm(0,j) ) / zs
      end do

      do j = 0, n
        do i = 1, m
          cqd(i,j) = ls * i * z / zs * cqm(i,j)
     &      + ( i + j ) * ( j - i + 1.0D+00 ) / zq * cqm(i-1,j)
        end do
      end do

      return
      end
      subroutine clqn ( n, x, y, cqn, cqd )

c*********************************************************************72
c
cc CLQN: Legendre function Qn(z) and derivative Wn'(z) for complex argument.
c
c  Licensing:
c
c    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
c    they give permission to incorporate this routine into a user program 
c    provided that the copyright is acknowledged.
c
c  Modified:
c
c    01 August 2012
c
c  Author:
c
c    Shanjie Zhang, Jianming Jin
c
c  Reference:
c
c    Shanjie Zhang, Jianming Jin,
c    Computation of Special Functions,
c    Wiley, 1996,
c    ISBN: 0-471-11963-6,
c    LC: QA351.C45.
c
c  Parameters:
c
c    Input, integer N, the degree of Qn(z).
c
c    Input, double precision X, Y, the real and imaginary parts of the argument Z.
c
c    Output, complex*16 CQN(0:N), CQD(0:N), the values of Qn(z) and Qn'(z.
c
      implicit none

      integer n

      complex*16 cq0
      complex*16 cq1
      complex*16 cqf0
      complex*16 cqf1
      complex*16 cqf2
      complex*16 cqn(0:n)
      complex*16 cqd(0:n)
      integer k
      integer km
      integer ls
      double precision x
      double precision y
      complex*16 z

      z = cmplx ( x, y )

      if ( z .eq. 1.0D+00 ) then
        do k = 0, n
          cqn(k) = cmplx ( 1.0D+30, 0.0D+00 )
          cqd(k) = cmplx ( 1.0D+30, 0.0D+00 )
        end do
        return
      end if

      if ( 1.0D+00 .lt. cdabs ( z ) ) then
        ls = -1
      else
        ls = +1
      end if

      cq0 = 0.5D+00 * cdlog ( ls * ( 1.0D+00 + z ) / ( 1.0D+00 - z ) )
      cq1 = z * cq0 - 1.0D+00
      cqn(0) = cq0
      cqn(1) = cq1

      if ( cdabs ( z ) .lt. 1.0001D+00 ) then

        cqf0 = cq0
        cqf1 = cq1
        do k = 2, n
          cqf2 = ( ( 2.0D+00 * k - 1.0D+00 ) * z * cqf1 
     &      - ( k - 1.0D+00 ) * cqf0 ) / k
          cqn(k) = cqf2
          cqf0 = cqf1
          cqf1 = cqf2
        end do

      else

        if ( 1.1D+00 .lt. cdabs ( z ) ) then
          km = 40 + n
        else
          km = ( 40 + n ) * int ( - 1.0D+00 
     &      - 1.8D+00 * log ( cdabs ( z - 1.0D+00 ) ) )
        end if

        cqf2 = 0.0D+00
        cqf1 = 1.0D+00
        do k = km, 0, -1
          cqf0 = ( ( 2 * k + 3.0D+00 ) * z * cqf1 
     &      - ( k + 2.0D+00 ) * cqf2 ) / ( k + 1.0D+00 )
          if ( k .le. n ) then
            cqn(k) = cqf0
          end if
          cqf2 = cqf1
          cqf1 = cqf0
        end do
        do k = 0, n
          cqn(k) = cqn(k) * cq0 / cqf0
        end do
      end if

      cqd(0) = ( cqn(1) - z * cqn(0) ) / ( z * z - 1.0D+00 )
      do k = 1, n
        cqd(k) = ( k * z * cqn(k) - k * cqn(k-1) ) / ( z * z - 1.0D+00 )
      end do

      return
      end
      subroutine comelp ( hk, ck, ce )

c*********************************************************************72
c
cc COMELP computes complete elliptic integrals K(k) and E(k).
c
c  Licensing:
c
c    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
c    they give permission to incorporate this routine into a user program 
c    provided that the copyright is acknowledged.
c
c  Modified:
c
c    07 July 2012
c
c  Author:
c
c    Shanjie Zhang, Jianming Jin
c
c  Reference:
c
c    Shanjie Zhang, Jianming Jin,
c    Computation of Special Functions,
c    Wiley, 1996,
c    ISBN: 0-471-11963-6,
c    LC: QA351.C45.
c
c  Parameters:
c
c    Input, double precision HK, the modulus.  0 <= HK <= 1.
c
c    Output, double precision CK, CE, the values of K(HK) and E(HK).
c
      implicit none

      double precision ae
      double precision ak
      double precision be
      double precision bk
      double precision ce
      double precision ck
      double precision hk
      double precision pk

      pk = 1.0D+00 - hk * hk

      if ( hk .eq. 1.0D+00 ) then

        ck = 1.0D+300
        ce = 1.0D+00

      else

        ak = (((
     &      0.01451196212D+00   * pk
     &    + 0.03742563713D+00 ) * pk
     &    + 0.03590092383D+00 ) * pk
     &    + 0.09666344259D+00 ) * pk
     &    + 1.38629436112D+00

        bk = (((
     &      0.00441787012D+00   * pk
     &    + 0.03328355346D+00 ) * pk
     &    + 0.06880248576D+00 ) * pk
     &    + 0.12498593597D+00 ) * pk
     &    + 0.5D+00

        ck = ak - bk * log ( pk )

        ae = (((
     &      0.01736506451D+00   * pk
     &    + 0.04757383546D+00 ) * pk
     &    + 0.0626060122D+00  ) * pk
     &    + 0.44325141463D+00 ) * pk
     &    + 1.0D+00

        be = (((
     &      0.00526449639D+00   * pk
     &    + 0.04069697526D+00 ) * pk
     &    + 0.09200180037D+00 ) * pk
     &    + 0.2499836831D+00  ) * pk

        ce = ae - be * log ( pk )

      end if

      return
      end
      subroutine cpbdn ( n, z, cpb, cpd )

c*********************************************************************72
c
cc CPBDN: parabolic cylinder function Dn(z) and Dn'(z) for complex argument.
c
c  Licensing:
c
c    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
c    they give permission to incorporate this routine into a user program 
c    provided that the copyright is acknowledged.
c
c  Modified:
c
c    29 July 2012
c
c  Author:
c
c    Shanjie Zhang, Jianming Jin
c
c  Reference:
c
c    Shanjie Zhang, Jianming Jin,
c    Computation of Special Functions,
c    Wiley, 1996,
c    ISBN: 0-471-11963-6,
c    LC: QA351.C45.
c
c  Parameters:
c
c    Input, integer N, the order.
c
c    Input, complex*16 Z, the argument.
c
c    Output, complex*16 CPB(0:N), CPD(0:N), the values of Dn(z) and Dn'(z).
c
      implicit none

      integer n

      double precision a0
      complex*16 c0
      complex*16 ca0
      complex*16 cf
      complex*16 cf0
      complex*16 cf1
      complex*16 cfa
      complex*16 cfb
      complex*16 cpb(0:n)
      complex*16 cpd(0:n)
      complex*16 cs0
      integer k
      integer m
      integer n0
      integer n1
      integer nm1
      double precision pi
      double precision x
      complex*16 z
      complex*16 z1

      pi = 3.141592653589793D+00
      x = real ( z )
      a0 = cdabs ( z )
      c0 = cmplx ( 0.0D+00, 0.0D+00 )
      ca0 = cdexp ( -0.25D+00 * z * z )

      if ( 0 .le. n ) then

        cf0 = ca0
        cf1 = z * ca0
        cpb(0) = cf0
        cpb(1) = cf1
        do k = 2, n
          cf = z * cf1 - ( k - 1.0D+00 ) * cf0
          cpb(k) = cf
          cf0 = cf1
          cf1 = cf
        end do

      else

        n0 = -n

        if ( x .le. 0.0D+00 .or. cdabs ( z ) .eq. 0.0D+00 ) then

          cf0 = ca0
          cpb(0) = cf0
          z1 = - z
          if ( a0 .le. 7.0D+00 ) then
            call cpdsa ( -1, z1, cf1 )
          else
            call cpdla ( -1, z1, cf1 )
          end if
          cf1 = dsqrt ( 2.0D+00 * pi ) / ca0 - cf1
          cpb(1) = cf1
          do k = 2, n0
            cf = ( - z * cf1 + cf0 ) / ( k - 1.0D+00 )
            cpb(k) = cf
            cf0 = cf1
            cf1 = cf
          end do

        else

          if ( a0 .le. 3.0D+00 ) then

            call cpdsa ( -n0, z, cfa )
            cpb(n0) = cfa
            n1 = n0 + 1
            call cpdsa ( -n1, z, cfb )
            cpb(n1) = cfb
            nm1 = n0 - 1
            do k = nm1, 0, -1
              cf = z * cfa + ( k + 1.0D+00 ) * cfb
              cpb(k) = cf
              cfb = cfa
              cfa = cf
            end do

          else

            m = 100 + abs ( n )
            cfa = c0
            cfb = cmplx ( 1.0D-30, 0.0D+00 )
            do k = m, 0, -1
              cf = z * cfb + ( k + 1.0D+00 ) * cfa
              if ( k .le. n0 ) then
                cpb(k) = cf
              end if
              cfa = cfb
              cfb = cf
            end do
            cs0 = ca0 / cf
            do k = 0, n0
              cpb(k) = cs0 * cpb(k)
            end do

          end if

        end if

      end if

      cpd(0) = -0.5D+00 * z * cpb(0)

      if ( 0 .le. n ) then
        do k = 1, n
          cpd(k) = -0.5D+00 * z * cpb(k) + k * cpb(k-1)
        end do
      else
        do k = 1, n0
          cpd(k) = 0.5D+00 * z * cpb(k) - cpb(k-1)
        end do
      end if

      return
      end
      subroutine cpdla ( n, z, cdn )

c*********************************************************************72
c
cc CPDLA computes complex parabolic cylinder function Dn(z) for large argument.
c
c  Licensing:
c
c    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
c    they give permission to incorporate this routine into a user program 
c    provided that the copyright is acknowledged.
c
c  Modified:
c
c    07 July 2012
c
c  Author:
c
c    Shanjie Zhang, Jianming Jin
c
c  Reference:
c
c    Shanjie Zhang, Jianming Jin,
c    Computation of Special Functions,
c    Wiley, 1996,
c    ISBN: 0-471-11963-6,
c    LC: QA351.C45.
c
c  Parameters:
c
c    Input, integer N, the order.
c
c    Input, complex*16 Z, the argument.
c
c    Output, complex*16 CDN, the function value.
c
      implicit none

      complex*16 cb0
      complex*16 cdn
      complex*16 cr
      integer k
      integer n
      complex*16 z

      cb0 = z ** n * cdexp ( -0.25D+00 * z * z )
      cr = cmplx ( 1.0D+00, 0.0D+00 )
      cdn = cmplx ( 1.0D+00, 0.0D+00 )

      do k = 1, 16

        cr = -0.5D+00 * cr * ( 2.0D+00 * k - n - 1.0D+00 ) 
     &    * ( 2.0D+00 * k - n - 2.0D+00 ) / ( k * z * z )

        cdn = cdn + cr

        if ( cdabs ( cr ) .lt. cdabs ( cdn ) * 1.0D-12 ) then
          go to 10
        end if

      end do

10    continue

      cdn = cb0 * cdn

      return
      end
      subroutine cpdsa ( n, z, cdn )

c*********************************************************************72
c
cc CPDSA computes complex parabolic cylinder function Dn(z) for small argument.
c
c  Licensing:
c
c    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
c    they give permission to incorporate this routine into a user program 
c    provided that the copyright is acknowledged.
c
c  Modified:
c
c    29 July 2012
c
c  Author:
c
c    Shanjie Zhang, Jianming Jin
c
c  Reference:
c
c    Shanjie Zhang, Jianming Jin,
c    Computation of Special Functions,
c    Wiley, 1996,
c    ISBN: 0-471-11963-6,
c    LC: QA351.C45.
c
c  Parameters:
c
c    Input, integer N, the order.
c
c    Input, complex*16 Z, the argument.
c
c    Output, complex*16 CDN, the value of DN(z).
c
      implicit none

      complex*16 ca0
      complex*16 cb0
      complex*16 cdn
      complex*16 cdw
      complex*16 cr
      double precision eps
      double precision g0
      double precision g1
      double precision ga0
      double precision gm
      integer m
      integer n
      double precision pd
      double precision pi
      double precision sq2
      double precision va0
      double precision vm
      double precision vt
      double precision xn
      complex*16 z

      eps = 1.0D-15
      pi = 3.141592653589793D+00
      sq2 = sqrt ( 2.0D+00 )
      ca0 = cdexp ( - 0.25D+00 * z * z )
      va0 = 0.5D+00 * ( 1.0D+00 - n )

      if ( n .eq. 0 ) then

        cdn = ca0

      else

        if ( cdabs ( z ) .eq. 0.0D+00 ) then

          if ( va0 .le. 0.0D+00 .and. va0 .eq. int ( va0 ) ) then
            cdn = 0.0D+00
          else
            call gaih ( va0, ga0 )
            pd = sqrt ( pi ) / ( 2.0D+00 ** ( -0.5D+00 * n ) * ga0 )
            cdn = cmplx ( pd, 0.0D+00 )
          end if

        else

          xn = - n
          call gaih ( xn, g1 )
          cb0 = 2.0D+00 ** ( -0.5D+00 * n - 1.0D+00 ) * ca0 / g1
          vt = -0.5D+00 * n
          call gaih ( vt, g0 )
          cdn = cmplx ( g0, 0.0D+00 )
          cr = cmplx ( 1.0D+00, 0.0D+00 )

          do m = 1, 250
            vm = 0.5D+00 * ( m - n )
            call gaih ( vm, gm )
            cr = - cr * sq2 * z / m
            cdw = gm * cr
            cdn = cdn + cdw
            if ( cdabs ( cdw ) .lt. cdabs ( cdn ) * eps ) then
              go to 10
            end if
          end do

10        continue

          cdn = cb0 * cdn

        end if

      end if

      return
      end
      subroutine cpsi ( x, y, psr, psi )

c*********************************************************************72
c
cc CPSI computes the psi function for a complex argument.
c
c  Licensing:
c
c    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
c    they give permission to incorporate this routine into a user program 
c    provided that the copyright is acknowledged.
c
c  Modified:
c
c    16 July 2012
c
c  Author:
c
c    Shanjie Zhang, Jianming Jin
c
c  Reference:
c
c    Shanjie Zhang, Jianming Jin,
c    Computation of Special Functions,
c    Wiley, 1996,
c    ISBN: 0-471-11963-6,
c    LC: QA351.C45.
c
c  Parameters:
c
c    Input, double precision X, Y, the real and imaginary parts of the argument.
c
c    Output, double precision PSR, PSI, the real and imaginary parts
c    of the function value.
c
      implicit none

      double precision a(8)
      double precision ct2
      integer k
      integer n
      double precision pi
      double precision psi
      double precision psr
      double precision ri
      double precision rr
      double precision th
      double precision tm
      double precision tn
      double precision x
      double precision x0
      double precision x1
      double precision y
      double precision y1
      double precision z0
      double precision z2

      save a

      data a / -0.8333333333333D-01, 0.83333333333333333D-02,
     &         -0.39682539682539683D-02, 0.41666666666666667D-02,
     &         -0.75757575757575758D-02, 0.21092796092796093D-01,
     &         -0.83333333333333333D-01, 0.4432598039215686D+00 /

      pi = 3.141592653589793D+00

      if ( y .eq. 0.0D+00 .and. x .eq. int ( x ) .and. 
     &  x .le. 0.0D+00 ) then

        psr = 1.0D+300
        psi = 0.0D+00

      else

        if ( x .lt. 0.0D+00 ) then
          x1 = x
          y1 = y
          x = -x
          y = -y
        end if

        x0 = x

        if ( x .lt. 8.0D+00 ) then
          n = 8 - int ( x )
          x0 = x + n
        end if

        if ( x0 .eq. 0.0D+00 ) then
          if ( y .ne. 0.0D+00 ) then
            th = 0.5D+00 * pi
          else
            th = 0.0D+00
          end if
        else
          th = atan ( y / x0 )
        end if

        z2 = x0 * x0 + y * y
        z0 = sqrt ( z2 )
        psr = log ( z0 ) - 0.5D+00 * x0 / z2
        psi = th + 0.5D+00 * y / z2
        do k = 1, 8
          psr = psr + a(k) * z2 ** ( - k ) * cos ( 2.0D+00 * k * th )
          psi = psi - a(k) * z2 ** ( - k ) * sin ( 2.0D+00 * k * th )
        end do

        if ( x .lt. 8.0D+00 ) then
          rr = 0.0D+00
          ri = 0.0D+00
          do k = 1, n
            rr = rr + ( x0 - k ) / ( ( x0 - k ) ** 2.0D+00 + y * y )
            ri = ri + y / ( ( x0 - k ) ** 2.0D+00 + y * y )
          end do
          psr = psr - rr
          psi = psi + ri
        end if

        if ( x1 .lt. 0.0D+00 ) then
          tn = tan ( pi * x )
          tm = tanh ( pi * y )
          ct2 = tn * tn + tm * tm
          psr = psr + x / ( x * x + y * y ) 
     &      + pi * ( tn - tn * tm * tm ) / ct2
          psi = psi - y / ( x * x + y * y ) 
     &      - pi * tm * ( 1.0D+00 + tn * tn ) / ct2
          x = x1
          y = y1
        end if

      end if

      return
      end
      subroutine csphik ( n, z, nm, csi, cdi, csk, cdk )

c*********************************************************************72
c
cc CSPHIK: modified spherical Bessel functions and derivatives for complex argument.
c
c  Licensing:
c
c    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
c    they give permission to incorporate this routine into a user program 
c    provided that the copyright is acknowledged.
c
c  Modified:
c
c    29 July 2012
c
c  Author:
c
c    Shanjie Zhang, Jianming Jin
c
c  Reference:
c
c    Shanjie Zhang, Jianming Jin,
c    Computation of Special Functions,
c    Wiley, 1996,
c    ISBN: 0-471-11963-6,
c    LC: QA351.C45.
c 
c  Parameters:
c
c    Input, integer N, the order of in(z) and kn(z).
c
c    Input, complex*16 Z, the argument.
c
c    Output, integer NM, the highest order computed.
c
c    Output, complex*16 CSI(0:N), CDI(0:N), CSK(0:N), CDK(0:N),
c    the values of in(z), in'(z), kn(z), kn'(z).
c
      implicit none

      integer n

      double precision a0
      complex*16 ccosh1
      complex*16 cdi(0:n)
      complex*16 cdk(0:n)
      complex*16 cf
      complex*16 cf0
      complex*16 cf1
      complex*16 ci
      complex*16 cs
      complex*16 csi(0:n)
      complex*16 csi0
      complex*16 csi1
      complex*16 csinh1
      complex*16 csk(0:n)
      integer k
      integer m
      integer msta1
      integer msta2
      integer nm
      double precision pi
      complex*16 z

      pi = 3.141592653589793D+00
      a0 = cdabs ( z )    
      nm = n

      if ( a0 .lt. 1.0D-60 ) then
        do k = 0, n
          csi(k) = 0.0D+00
          cdi(k) = 0.0D+00
          csk(k) = 1.0D+300
          cdk(k) = -1.0D+300
        end do
        csi(0) = 1.0D+00
        cdi(1) = 0.3333333333333333D+00
        return
      end if

      ci = cmplx ( 0.0D+00, 1.0D+00 )
      csinh1 = cdsin ( ci * z ) / ci
      ccosh1 = cdcos ( ci * z )
      csi0 = csinh1 / z
      csi1 = ( - csinh1 / z + ccosh1 ) / z
      csi(0) = csi0
      csi(1) = csi1

      if ( 2 .le. n ) then

        m = msta1 ( a0, 200 )
        if ( m .lt. n ) then
          nm = m
        else
          m = msta2 ( a0, n, 15 )
        end if

        cf0 = 0.0D+00
        cf1 = 1.0D+00-100
        do k = m, 0, -1
          cf = ( 2.0D+00 * k + 3.0D+00 ) * cf1 / z + cf0
          if ( k .le. nm ) then
            csi(k) = cf
          end if
          cf0 = cf1
          cf1 = cf
        end do

        if ( cdabs ( csi0 ) .le. cdabs ( csi1 ) ) then
          cs = csi1 / cf0
        else
          cs = csi0 / cf
        end if

        do k = 0, nm
          csi(k) = cs * csi(k)
        end do

      end if

      cdi(0) = csi(1)
      do k = 1, nm
        cdi(k) = csi(k-1) - ( k + 1.0D+00 ) * csi(k) / z
      end do

      csk(0) = 0.5D+00 * pi / z * cdexp ( - z )
      csk(1) = csk(0) * ( 1.0D+00 + 1.0D+00 / z )
      do k = 2, nm
        if ( cdabs ( csi(k-2) ) .lt. cdabs ( csi(k-1) ) ) then
          csk(k) = ( 0.5D+00 * pi / ( z * z ) - csi(k) * csk(k-1) ) 
     &      / csi(k-1)
        else
          csk(k) = ( csi(k) * csk(k-2) + ( k - 0.5D+00 ) 
     &      * pi / z ** 3 ) / csi(k-2)
        end if
      end do

      cdk(0) = -csk(1)
      do k = 1, nm
        cdk(k) = - csk(k-1) - ( k + 1.0D+00 ) * csk(k) / z
      end do

      return
      end
      subroutine csphjy ( n, z, nm, csj, cdj, csy, cdy )

c*********************************************************************72
c
cc CSPHJY computes spherical Bessel functions jn(z) and yn(z) for complex argument.
c
c  Discussion:
c
c    This procedure computes spherical Bessel functions jn(z) and yn(z)
c    and their derivatives for a complex argument.
c
c  Licensing:
c
c    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
c    they give permission to incorporate this routine into a user program 
c    provided that the copyright is acknowledged.
c
c  Modified:
c
c    01 August 2012
c
c  Author:
c
c    Shanjie Zhang, Jianming Jin
c
c  Reference:
c
c    Shanjie Zhang, Jianming Jin,
c    Computation of Special Functions,
c    Wiley, 1996,
c    ISBN: 0-471-11963-6,
c    LC: QA351.C45.
c 
c  Parameters:
c
c    Input, integer N, the order of jn(z) and yn(z).
c
c    Input, complex*16 Z, the argument.
c
c    Output, integer NM, the highest order computed.
c
c    Output, complex*16 CSJ(0:N0, CDJ(0:N), CSY(0:N), CDY(0:N),
c    the values of jn(z), jn'(z), yn(z), yn'(z).
c
      implicit none

      integer n

      double precision a0
      complex*16 csj(0:n)
      complex*16 cdj(0:n)
      complex*16 csy(0:n)
      complex*16 cdy(0:n)
      complex*16 cf
      complex*16 cf0
      complex*16 cf1
      complex*16 cs
      complex*16 csa
      complex*16 csb
      integer k
      integer m
      integer msta1
      integer msta2
      integer nm
      complex*16 z

      a0 = cdabs ( z )
      nm = n

      if ( a0 .lt. 1.0D-60 ) then
        do k = 0, n
          csj(k) = 0.0D+00
          cdj(k) = 0.0D+00
          csy(k) = -1.0D+300
          cdy(k) = 1.0D+300
        end do
        csj(0) = cmplx ( 1.0D+00, 0.0D+00 )
        cdj(1) = cmplx ( 0.333333333333333D+00, 0.0D+00 )
        return
      end if

      csj(0) = cdsin ( z ) / z
      csj(1) = ( csj(0) - cdcos ( z ) ) / z

      if ( 2 .le. n ) then
        csa = csj(0)
        csb = csj(1)
        m = msta1 ( a0, 200 )
        if ( m .lt. n ) then
          nm = m
        else
          m = msta2 ( a0, n, 15 )
        end if
        cf0 = 0.0D+00
        cf1 = 1.0D+00-100
        do k = m, 0, -1
          cf = ( 2.0D+00 * k + 3.0D+00 ) * cf1 / z - cf0
          if ( k .le. nm ) then
            csj(k) = cf
          end if
          cf0 = cf1
          cf1 = cf
        end do

        if ( cdabs ( csa ) .le. cdabs ( csb ) ) then
          cs = csb / cf0
        else
          cs = csa / cf
        end if

        do k = 0, nm
          csj(k) = cs * csj(k)
        end do

      end if

      cdj(0) = ( cdcos ( z ) - cdsin ( z ) / z ) / z
      do k = 1, nm
        cdj(k) = csj(k-1) - ( k + 1.0D+00 ) * csj(k) / z
      end do
      csy(0) = - cdcos ( z ) / z
      csy(1) = ( csy(0) - cdsin ( z ) ) / z
      cdy(0) = ( cdsin ( z ) + cdcos ( z ) / z ) / z
      cdy(1) = ( 2.0D+00 * cdy(0) - cdcos ( z ) )  / z

      do k = 2, nm
        if ( cdabs ( csj(k-2) ) .lt. cdabs ( csj(k-1) ) ) then
          csy(k) = ( csj(k) * csy(k-1) - 1.0D+00 / ( z * z ) ) 
     &      / csj(k-1)
        else
          csy(k) = ( csj(k) * csy(k-2)
     &      - ( 2.0D+00 * k - 1.0D+00 ) / z ** 3 ) / csj(k-2)
        end if
      end do

      do k = 2, nm
        cdy(k) = csy(k-1) - ( k + 1.0D+00 ) * csy(k) / z
      end do

      return
      end
      subroutine cv0 ( kd, m, q, a0 )

c*********************************************************************72
c
cc CV0 computes the initial characteristic value of Mathieu functions.
c
c  Discussion:
c
c    This procedure computes the initial characteristic value of Mathieu 
c    functions for m <= 12 or q <= 300 or q <= m*m.
c
c  Licensing:
c
c    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
c    they give permission to incorporate this routine into a user program 
c    provided that the copyright is acknowledged.
c
c  Modified:
c
c   03 August 2012
c
c  Author:
c
c    Shanjie Zhang, Jianming Jin
c
c  Reference:
c
c    Shanjie Zhang, Jianming Jin,
c    Computation of Special Functions,
c    Wiley, 1996,
c    ISBN: 0-471-11963-6,
c    LC: QA351.C45.
c
c  Parameters:
c
c    Input, integer KD, the case code:
c    1, for cem(x,q)  ( m = 0,2,4,...)
c    2, for cem(x,q)  ( m = 1,3,5,...)
c    3, for sem(x,q)  ( m = 1,3,5,...)
c    4, for sem(x,q)  ( m = 2,4,6,...)
c
c    Input, integer M, the order of the functions.
c
c    Input, double precision Q, the parameter of the functions.
c
c    Output, double precision A0, the characteristic value.
c
      implicit none

      double precision a0
      integer kd
      integer m
      double precision q
      double precision q2

      q2 = q * q

      if ( m .eq. 0 ) then

        if ( q .le. 1.0D+00 ) then

          a0 = (((
     &        0.0036392D+00   * q2
     &      - 0.0125868D+00 ) * q2
     &      + 0.0546875D+00 ) * q2
     &      - 0.5D+00 )       * q2

        else if ( q .le. 10.0D+00 ) then

          a0 = ((
     &        3.999267D-03   * q
     &      - 9.638957D-02 ) * q
     &      - 0.88297D+00 )  * q
     &      + 0.5542818D+00 

        else

          call cvql ( kd, m, q, a0 )

        end if

      else if ( m .eq. 1 ) then

        if ( q .le. 1.0D+00 .and. kd .eq. 2 ) then

          a0 = (((
     &      - 6.51D-04 * q
     &      - 0.015625D+00 ) * q
     &      - 0.125D+00 ) * q
     &      + 1.0D+00 ) * q
     &      + 1.0D+00 

        else if ( q .le. 1.0D+00 .and. kd .eq. 3 ) then

          a0 = (((
     &      - 6.51D-04 * q
     &      + 0.015625D+00 ) * q
     &      - 0.125D+00 ) * q
     &      - 1.0D+00 ) * q
     &      + 1.0D+00
 
        else if ( q .le. 10.0D+00 .and. kd .eq. 2 ) then

          a0 = (((
     &      - 4.94603D-04 * q
     &      + 1.92917D-02 ) * q
     &      - 0.3089229D+00 ) * q
     &      + 1.33372D+00 ) * q
     &      + 0.811752D+00 

        else if ( q .le. 10.0D+00 .and. kd .eq. 3 ) then

          a0 = ((
     &        1.971096D-03 * q
     &      - 5.482465D-02 ) * q
     &      - 1.152218D+00 ) * q
     &      + 1.10427D+00 

        else

          call cvql ( kd, m, q, a0 )

        end if

      else if ( m .eq. 2 ) then

        if ( q .le. 1.0D+00 .and. kd .eq. 1 ) then

          a0 = (((
     &      - 0.0036391D+00   * q2
     &      + 0.0125888D+00 ) * q2
     &      - 0.0551939D+00 ) * q2
     &      + 0.416667D+00 )  * q2 + 4.0D+00 

        else if ( q .le. 1.0D+00 .and. kd .eq. 4 ) then

          a0 = ( 
     &        0.0003617D+00 * q2 
     &      - 0.0833333D+00 ) * q2 + 4.0D+00 

        else if ( q .le. 15.0D+00 .and. kd .eq. 1 ) then

          a0 = (((
     &        3.200972D-04    * q
     &      - 8.667445D-03 )  * q
     &      - 1.829032D-04 )  * q
     &      + 0.9919999D+00 ) * q
     &      + 3.3290504D+00 

        else if ( q .le. 10.0D+00 .and. kd .eq. 4 ) then

          a0 = ((
     &        2.38446D-03 * q
     &      - 0.08725329D+00 ) * q
     &      - 4.732542D-03 ) * q
     &      + 4.00909D+00 

        else

          call cvql ( kd, m, q, a0 )

        end if

      else if ( m .eq. 3 ) then

        if ( q .le. 1.0D+00 .and. kd .eq. 2 ) then
          a0 = ((
     &        6.348D-04 * q
     &      + 0.015625D+00 ) * q
     &      + 0.0625 ) * q2 
     &      + 9.0D+00 
        else if ( q .le. 1.0D+00 .and. kd .eq. 3 ) then
          a0 = ((
     &        6.348D-04 * q
     &      - 0.015625D+00 ) * q
     &      + 0.0625D+00 ) * q2
     &      + 9.0D+00 
        else if ( q .le. 20.0D+00 .and. kd .eq. 2 ) then
          a0 = (((
     &        3.035731D-04 * q
     &      - 1.453021D-02 ) * q
     &      + 0.19069602D+00 ) * q
     &      - 0.1039356D+00 ) * q
     &      + 8.9449274D+00 
        else if ( q .le. 15.0D+00 .and. kd .eq. 3 ) then
          a0 = ((
     &        9.369364D-05 * q
     &      - 0.03569325D+00 ) * q
     &      + 0.2689874D+00 ) * q
     &      + 8.771735D+00 
        else
          call cvql ( kd, m, q, a0 )
        end if

      else if ( m .eq. 4 ) then

        if ( q .le. 1.0D+00 .and. kd .eq. 1 ) then
          a0 = ((
     &      - 2.1D-06 * q2
     &      + 5.012D-04 ) * q2
     &      + 0.0333333 ) * q2
     &      + 16.0D+00
        else if ( q .le. 1.0D+00 .and. kd .eq. 4 ) then
          a0 = ((
     &        3.7D-06 * q2
     &      - 3.669D-04 ) * q2
     &      + 0.0333333D+00 ) * q2
     &      + 16.0D+00
        else if ( q .le. 25.0D+00 .and. kd .eq. 1 ) then
          a0 = (((
     &        1.076676D-04 * q
     &      - 7.9684875D-03 ) * q
     &      + 0.17344854D+00 ) * q
     &      - 0.5924058D+00 ) * q
     &      + 16.620847D+00
        else if ( q .le. 20.0D+00 .and. kd .eq. 4 ) then
          a0 = ((
     &      - 7.08719D-04 * q
     &      + 3.8216144D-03 ) * q
     &      + 0.1907493D+00 ) * q
     &      + 15.744D+00
        else
          call cvql ( kd, m, q, a0 )
        end if

      else if ( m .eq. 5 ) then

        if ( q .le. 1.0D+00 .and. kd .eq. 2 ) then
          a0 = ((
     &        6.8D-6 * q
     &      + 1.42D-05 ) * q2
     &      + 0.0208333D+00 ) * q2
     &      + 25.0D+00
        else if ( q .le. 1.0D+00 .and. kd .eq. 3 ) then
          a0 = ((
     &      - 6.8D-06 * q
     &      + 1.42D-05 ) * q2
     &      + 0.0208333D+00 ) * q2
     &      + 25.0D+00
        else if ( q .le. 35.0D+00 .and. kd .eq. 2 ) then
          a0 = (((
     &        2.238231D-05 * q
     &      - 2.983416D-03 ) * q
     &      + 0.10706975D+00 ) * q
     &      - 0.600205D+00 ) * q
     &      + 25.93515D+00
        else if ( q .le. 25.0D+00 .and. kd .eq. 3 ) then
          a0 = ((
     &      - 7.425364D-04 * q
     &      + 2.18225D-02 ) * q
     &      + 4.16399D-02 ) * q
     &      + 24.897D+00
        else
          call cvql ( kd, m, q, a0 )
        end if

      else if ( m .eq. 6 ) then

        if ( q .le. 1.0D+00 ) then
          a0 = ( 0.4D-06 * q2 + 0.0142857 ) * q2 + 36.0D+00
        else if ( q .le. 40.0D+00 .and. kd .eq. 1 ) then
          a0 = (((
     &      - 1.66846D-05 * q
     &      + 4.80263D-04 ) * q
     &      + 2.53998D-02 ) * q
     &      - 0.181233D+00 ) * q 
     &      + 36.423D+00
        else if ( q .le. 35.0D+00 .and. kd .eq. 4 ) then
          a0 = ((
     &      - 4.57146D-04 * q
     &      + 2.16609D-02 ) * q
     &      - 2.349616D-02 ) * q
     &      + 35.99251D+00
        else
          call cvql ( kd, m, q, a0 )
        end if

      else if ( m .eq. 7 ) then

        if ( q .le. 10.0D+00 ) then
          call cvqm ( m, q, a0 )
        else if ( q .le. 50.0D+00 .and. kd .eq. 2 ) then
          a0 = (((
     &      - 1.411114D-05 * q
     &      + 9.730514D-04 ) * q
     &      - 3.097887D-03 ) * q
     &      + 3.533597D-02 ) * q
     &      + 49.0547D+00
        else if ( q .le. 40.0D+00 .and. kd .eq. 3 ) then
          a0 = ((
     &      - 3.043872D-04 * q
     &      + 2.05511D-02 ) * q
     &      - 9.16292D-02 ) * q
     &      + 49.19035D+00
        else
          call cvql ( kd, m, q, a0 )
        end if

      else if ( 8 .le. m ) then

        if ( q .le. 3.0D+00 * m ) then
          call cvqm ( m, q, a0 )
        else if ( m * m .lt. q ) then
          call cvql ( kd, m, q, a0 )
        else if ( m .eq. 8 .and. kd .eq. 1 ) then
          a0 = (((
     &        8.634308D-06 * q
     &      - 2.100289D-03 ) * q
     &      + 0.169072D+00 ) * q
     &      - 4.64336D+00 ) * q
     &      + 109.4211D+00
        else if ( m .eq. 8 .and. kd .eq. 4 ) then
          a0 = ((
     &      - 6.7842D-05 * q
     &      + 2.2057D-03 ) * q
     &      + 0.48296D+00 ) * q
     &      + 56.59D+00
        else if ( m .eq. 9 .and. kd .eq. 2 ) then
          a0 = (((
     &        2.906435D-06 * q
     &      - 1.019893D-03 ) * q
     &      + 0.1101965D+00 ) * q
     &      - 3.821851D+00 ) * q
     &      + 127.6098D+00
        else if ( m .eq. 9 .and. kd .eq. 3 ) then
          a0 = ((
     &      - 9.577289D-05 * q
     &      + 0.01043839D+00 ) * q
     &      + 0.06588934D+00 ) * q
     &      + 78.0198D+00
        else if ( m .eq. 10 .and. kd .eq. 1 ) then
          a0 = (((
     &        5.44927D-07 * q
     &      - 3.926119D-04 ) * q
     &      + 0.0612099D+00 ) * q
     &      - 2.600805D+00 ) * q
     &      + 138.1923D+00
        else if ( m .eq. 10 .and. kd .eq. 4 ) then
          a0 = ((
     &      - 7.660143D-05 * q
     &      + 0.01132506D+00 ) * q
     &      - 0.09746023D+00 ) * q
     &      + 99.29494D+00
        else if ( m .eq. 11 .and. kd .eq. 2 ) then
          a0 = (((
     &      - 5.67615D-07 * q
     &      + 7.152722D-06 ) * q
     &      + 0.01920291D+00 ) * q
     &      - 1.081583D+00 ) * q
     &      + 140.88D+00
        else if ( m .eq. 11 .and. kd .eq. 3 ) then
          a0 = ((
     &      - 6.310551D-05 * q
     &      + 0.0119247D+00 ) * q
     &      - 0.2681195D+00 ) * q
     &      + 123.667D+00
        else if ( m .eq. 12 .and. kd .eq. 1 ) then
          a0 = (((
     &      - 2.38351D-07 * q
     &      - 2.90139D-05 ) * q
     &      + 0.02023088D+00 ) * q
     &      - 1.289D+00 ) * q
     &      + 171.2723D+00
        else if ( m .eq. 12 .and. kd .eq. 4 ) then
          a0 = (((
     &        3.08902D-07 * q
     &      - 1.577869D-04 ) * q
     &      + 0.0247911D+00 ) * q
     &      - 1.05454D+00 ) * q 
     &      + 161.471D+00

        end if

      end if

      return
      end
      subroutine cva1 ( kd, m, q, cv )

c*********************************************************************72
c
cc CVA1 computes a sequence of characteristic values of Mathieu functions.
c
c  Licensing:
c
c    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
c    they give permission to incorporate this routine into a user program 
c    provided that the copyright is acknowledged.
c
c  Modified:
c
c    25 July 2012
c
c  Author:
c
c    Shanjie Zhang, Jianming Jin
c
c  Reference:
c
c    Shanjie Zhang, Jianming Jin,
c    Computation of Special Functions,
c    Wiley, 1996,
c    ISBN: 0-471-11963-6,
c    LC: QA351.C45.
c
c  Parameters:
c
c    Input, integer KD, the case code.
c    1, for cem(x,q)  ( m = 0,2,4,˙˙˙ )
c    2, for cem(x,q)  ( m = 1,3,5,˙˙˙ )
c    3, for sem(x,q)  ( m = 1,3,5,˙˙˙ )
c    4, for sem(x,q)  ( m = 2,4,6,˙˙˙ )
c
c    Input, integer M, the maximum order of the Mathieu functions.
c
c    Input, double precision Q, the parameter of the Mathieu functions.
c
c    Output, double precision CV(*), characteristic values.
c    For KD = 1, CV(1), CV(2), CV(3),..., correspond to
c    the characteristic values of cem for m = 0,2,4,...
c    For KD = 2, CV(1), CV(2), CV(3),..., correspond to
c    the characteristic values of cem for m = 1,3,5,...
c    For KD = 3, CV(1), CV(2), CV(3),..., correspond to
c    the characteristic values of sem for m = 1,3,5,...
c    For KD = 4, CV(1), CV(2), CV(3),..., correspond to
c    the characteristic values of sem for m = 0,2,4,...
c       
      implicit none

      double precision cv(200)
      double precision d(500)
      double precision e(500)
      double precision eps
      double precision f(500)
      double precision g(200)
      double precision h(200)
      integer i
      integer ic
      integer icm
      integer j
      integer k
      integer k1
      integer kd
      integer m
      integer nm
      integer nm1
      double precision q
      double precision s
      double precision t
      double precision t1
      double precision x1
      double precision xa
      double precision xb

      eps = 1.0D-14

      if ( kd .eq. 4 ) then
        icm = m / 2
      else
        icm = int ( m / 2 ) + 1
      end if

      if ( q .eq. 0.0D+00 ) then

        if ( kd .eq. 1 ) then
          do ic = 1, icm
            cv(ic) = 4.0D+00 * ( ic - 1.0D+00 ) ** 2
          end do
        else if ( kd .ne. 4 ) then
          do ic = 1, icm
            cv(ic) = ( 2.0D+00 * ic - 1.0D+00 ) ** 2
          end do
        else
          do ic = 1, icm
            cv(ic) = 4.0D+00 * ic * ic
          end do
        end if

      else

        nm = int ( 10D+00 + 1.5D+00 * m + 0.5D+00 * q )
        e(1) = 0.0D+00
        f(1) = 0.0D+00

        if ( kd .eq. 1 ) then

          d(1) = 0.0D+00
          do i = 2, nm
            d(i) = 4.0D+00 * ( i - 1.0D+00 ) ** 2
            e(i) = q
            f(i) = q * q
          end do
          e(2) = sqrt ( 2.0D+00 ) * q
          f(2) = 2.0D+00 * q * q

        else if ( kd .ne. 4 ) then

          d(1) = 1.0D+00 + ( -1.0D+00 ) ** kd * q
          do i = 2, nm
            d(i) = ( 2.0D+00 * i - 1.0D+00 ) ** 2
            e(i) = q
            f(i) = q * q
          end do

        else

          d(1) = 4.0D+00
          do i = 2, nm
            d(i) = 4.0D+00 * i * i
            e(i) = q
            f(i) = q * q
          end do

        end if

        xa = d(nm) + abs ( e(nm) )
        xb = d(nm) - abs ( e(nm) )

        nm1 = nm - 1
        do i = 1, nm1
          t = abs ( e(i) ) + abs ( e(i+1) )
          t1 = d(i) + t
          xa = max ( xa, t1 )
          t1 = d(i) - t
          xb = min ( xb, t1 )
        end do

        do i = 1, icm
          g(i) = xa
          h(i) = xb
        end do

        do k = 1, icm

          do k1 = k, icm
            if ( g(k1) .lt. g(k) ) then
              g(k) = g(k1)
              go to 10
            end if
          end do

10        continue

          if ( k .ne. 1 .and. h(k) .lt. h(k-1) ) then
            h(k) = h(k-1)
          end if

20        continue

          x1 = ( g(k) + h(k) ) /2.0D+00
          cv(k) = x1

          if ( abs ( ( g(k) - h(k) ) / x1 ) .lt. eps ) then
            go to 30
          end if

          j = 0
          s = 1.0D+00
          do i = 1, nm
            if ( s .eq. 0.0D+00 ) then
              s = s + 1.0D-30
            end if
            t = f(i) / s
            s = d(i) - t - x1
            if ( s .lt. 0.0D+00 ) then
              j = j + 1
            end if
          end do

          if ( j .lt. k ) then
            h(k) = x1
          else
            g(k) = x1
            if ( icm .le. j ) then
              g(icm) = x1
            else
              h(j+1) = max ( h(j+1), x1 )
              g(j) = min ( g(j), x1 )
            end if
          end if

          go to 20

30        continue

          cv(k) = x1

        end do

      end if

      return
      end
      subroutine cva2 ( kd, m, q, a )

c*********************************************************************72
c
cc CVA2 computes a specific characteristic value of Mathieu functions.
c
c  Licensing:
c
c    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
c    they give permission to incorporate this routine into a user program 
c    provided that the copyright is acknowledged.
c
c  Modified:
c
c    22 July 2012
c
c  Author:
c
c    Shanjie Zhang, Jianming Jin
c
c  Reference:
c
c    Shanjie Zhang, Jianming Jin,
c    Computation of Special Functions,
c    Wiley, 1996,
c    ISBN: 0-471-11963-6,
c    LC: QA351.C45.
c
c  Parameters:
c
c    Input, integer KD, the case code:
c    1, for cem(x,q)  ( m = 0,2,4,...)
c    2, for cem(x,q)  ( m = 1,3,5,...)
c    3, for sem(x,q)  ( m = 1,3,5,...)
c    4, for sem(x,q)  ( m = 2,4,6,...)
c
c    Input, integer M, the order of the Mathieu functions.
c
c    Input, double precision Q, the parameter of the Mathieu functions.
c
c    Output, double precision A, the characteristic value.
c
      implicit none

      double precision a
      double precision a1
      double precision a2
      double precision delta
      integer i
      integer iflag
      integer kd
      integer m
      integer ndiv
      integer nn
      double precision q
      double precision q1
      double precision q2
      double precision qq

      if ( m .le. 12 .or. q .le. 3.0D+00 * m .or. m * m .lt. q ) then

        call cv0 ( kd, m, q, a )

        if ( q .ne. 0.0D+00 ) then
          call refine ( kd, m, q, a, 1 )
        end if

      else

        ndiv = 10
        delta = ( m - 3.0D+00 ) * m / dble ( ndiv )

        if ( ( q - 3.0D+00 * m ) .le. ( m * m - q ) ) then

10        continue

          nn = int ( ( q - 3.0D+00 * m ) / delta ) + 1
          delta = ( q - 3.0D+00 * m ) / nn
          q1 = 2.0D+00 * m
          call cvqm ( m, q1, a1 )
          q2 = 3.0D+00 * m
          call cvqm ( m, q2, a2 )
          qq = 3.0D+00 * m
 
          do i = 1, nn

            qq = qq + delta
            a = ( a1 * q2 - a2 * q1 + ( a2 - a1 ) * qq ) / ( q2 - q1 )
 
            if ( i .eq. nn ) then
              iflag = -1
            else
              iflag = 1
            end if

            call refine ( kd, m, qq, a, iflag )
            q1 = q2
            q2 = qq
            a1 = a2
            a2 = a

          end do

          if ( iflag .eq. -10 ) then
            ndiv = ndiv * 2
            delta = ( m - 3.0D+00 ) * m / dble ( ndiv )
            go to 10
          end if

        else

20        continue

          nn = int ( ( m * m - q ) / delta ) + 1
          delta = ( m * m - q ) / nn
          q1 = m * ( m - 1.0D+00 )
          call cvql ( kd, m, q1, a1 )
          q2 = m * m
          call cvql ( kd, m, q2, a2 )
          qq = m * m

          do i = 1, nn

            qq = qq - delta
            a = ( a1 * q2 - a2 * q1 + ( a2 - a1 ) * qq ) / ( q2 - q1 )

            if ( i .eq. nn ) then
              iflag = -1
            else
              iflag = 1
            end if

            call refine ( kd, m, qq, a, iflag )
            q1 = q2
            q2 = qq
            a1 = a2
            a2 = a

          end do

          if ( iflag .eq. -10 ) then
            ndiv = ndiv * 2
            delta = ( m - 3.0D+00 ) * m / dble ( ndiv )
            go to 20
          end if

        end if

      end if

      return
      end
      subroutine cvf ( kd, m, q, a, mj, f )

c*********************************************************************72
c
cc CVF computes the value of F for the characteristic equation of Mathieu functions.
c
c  Licensing:
c
c    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
c    they give permission to incorporate this routine into a user program 
c    provided that the copyright is acknowledged.
c
c  Modified:
c
c    16 July 2012
c
c  Author:
c
c    Shanjie Zhang, Jianming Jin
c
c  Reference:
c
c    Shanjie Zhang, Jianming Jin,
c    Computation of Special Functions,
c    Wiley, 1996,
c    ISBN: 0-471-11963-6,
c    LC: QA351.C45.
c
c  Parameters:
c
c    Input, integer KD, the case code:
c    1, for cem(x,q)  ( m = 0,2,4,...)
c    2, for cem(x,q)  ( m = 1,3,5,...)
c    3, for sem(x,q)  ( m = 1,3,5,...)
c    4, for sem(x,q)  ( m = 2,4,6,...)
c
c    Input, integer M, the order of the Mathieu functions.
c
c    Input, double precision Q, the parameter of the Mathieu functions.
c
c    Input, double precision A, the characteristic value.
c
c    Input, integer MJ, ?
c
c    Output, double precision F, the value of the function for the
c    characteristic equation.
c
      implicit none

      double precision a
      double precision b
      double precision f
      integer ic
      integer j
      integer j0
      integer jf
      integer kd
      integer l
      integer l0
      integer m
      integer mj
      double precision q
      double precision t0
      double precision t1
      double precision t2

      b = a
      ic = int ( m / 2 )
      l = 0
      l0 = 0
      j0 = 2
      jf = ic

      if ( kd .eq. 1 ) then
        l0 = 2
        j0 = 3
      else if ( kd .eq. 2 .or. kd .eq. 3 ) then
        l = 1
      else if ( kd .eq. 4 ) then
        jf = ic - 1
      end if

      t1 = 0.0D+00
      do j = mj, ic + 1, -1
        t1 = - q * q / ( ( 2.0D+00 * j + l ) ** 2 - b + t1 )
      end do

      if ( m .le. 2 ) then

        t2 = 0.0D+00

        if ( kd .eq. 1 ) then
          if ( m .eq. 0 ) then
            t1 = t1 + t1
          else if ( m .eq. 2 ) then
            t1 = - 2.0D+00 * q * q / ( 4.0D+00 - b + t1 ) - 4.0D+00
          end if
        else if ( kd .eq. 2 ) then
          if ( m .eq. 1 ) then
            t1 = t1 + q
          end if
        else if ( kd .eq. 3 ) then
          if ( m .eq. 1 ) then
            t1 = t1 - q
          end if
        end if

      else

        if ( kd .eq. 1 ) then
          t0 = 4.0D+00 - b + 2.0D+00 * q * q / b
        else if ( kd .eq. 2 ) then
          t0 = 1.0D+00 - b + q
        else if ( kd .eq. 3 ) then
          t0 = 1.0D+00 - b - q
        else if ( kd .eq. 4 ) then
          t0 = 4.0D+00 - b
        end if

        t2 = - q * q / t0
        do j = j0, jf
          t2 = - q * q / ( ( 2.0D+00 * j - l - l0 ) ** 2 - b + t2 )
        end do

      end if

      f = ( 2.0D+00 * ic + l ) ** 2 + t1 + t2 - b

      return
      end
      subroutine cvql ( kd, m, q, a0 )

c*********************************************************************72
c
cc CVQL computes the characteristic value of Mathieu functions for q <= 3*m.
c
c  Licensing:
c
c    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
c    they give permission to incorporate this routine into a user program 
c    provided that the copyright is acknowledged.
c
c  Modified:
c
c    10 July 2012
c
c  Author:
c
c    Shanjie Zhang, Jianming Jin
c
c  Reference:
c
c    Shanjie Zhang, Jianming Jin,
c    Computation of Special Functions,
c    Wiley, 1996,
c    ISBN: 0-471-11963-6,
c    LC: QA351.C45.
c
c  Parameters:
c
c    Input, integer KD, the case code:
c    1, for cem(x,q)  ( m = 0,2,4,...)
c    2, for cem(x,q)  ( m = 1,3,5,...)
c    3, for sem(x,q)  ( m = 1,3,5,...)
c    4, for sem(x,q)  ( m = 2,4,6,...)
c
c    Input, integer M, the order of the Mathieu functions.
c
c    Input, double precision Q, the parameter value.
c
c    Output, double precision A0, the initial characteristic value.
c
      implicit none

      double precision a0
      double precision c1
      double precision cv1
      double precision cv2
      double precision d1
      double precision d2
      double precision d3
      double precision d4
      integer kd
      integer m
      double precision p1
      double precision p2
      double precision q
      double precision w
      double precision w2
      double precision w3
      double precision w4
      double precision w6

      if ( kd .eq. 1 .or. kd .eq. 2 ) then
        w = 2.0D+00 * m + 1.0D+00
      else
        w = 2.0D+00 * m - 1.0D+00
      end if

      w2 = w * w
      w3 = w * w2
      w4 = w2 * w2
      w6 = w2 * w4
      d1 = 5.0D+00 + 34.0D+00 / w2 + 9.0D+00 / w4
      d2 = ( 33.0D+00 + 410.0D+00 / w2 + 405.0D+00 / w4 ) / w
      d3 = ( 63.0D+00 + 1260.0D+00 / w2 + 2943.0D+00 / w4 
     &  + 486.0D+00 / w6 ) / w2
      d4 = ( 527.0D+00 + 15617.0D+00 / w2 + 69001.0D+00 / w4 
     &  + 41607.0D+00 / w6 ) / w3
      c1 = 128.0D+00
      p2 = q / w4
      p1 = sqrt ( p2 )
      cv1 = - 2.0D+00 * q + 2.0D+00 * w * sqrt ( q ) 
     &  - ( w2 + 1.0D+00 ) / 8.0D+00
      cv2 = ( w + 3.0D+00 / w ) + d1 / ( 32.0D+00 * p1 ) + d2 
     &  / ( 8.0D+00 * c1 * p2 )
      cv2 = cv2 + d3 / ( 64.0D+00 * c1 * p1 * p2 ) + d4 
     &  / ( 16.0D+00 * c1 * c1 * p2 * p2 )
      a0 = cv1 - cv2 / ( c1 * p1 )

      return
      end
      subroutine cvqm ( m, q, a0 )

c*********************************************************************72
c
cc CVQM computes the characteristic value of Mathieu functions for q <= m*m.
c
c  Licensing:
c
c    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
c    they give permission to incorporate this routine into a user program 
c    provided that the copyright is acknowledged.
c
c  Modified:
c
c    27 March 2012
c
c  Author:
c
c    Shanjie Zhang, Jianming Jin
c
c  Reference:
c
c    Shanjie Zhang, Jianming Jin,
c    Computation of Special Functions,
c    Wiley, 1996,
c    ISBN: 0-471-11963-6,
c    LC: QA351.C45.
c
c  Parameters:
c
c    Input, integer M, the order of the Mathieu functions.
c
c    Input, double precision Q, the parameter value.
c
c    Output, double precision A0, the initial characteristic value.
c
      implicit none

      double precision a0
      double precision hm1
      double precision hm3
      double precision hm5
      integer m
      double precision q

      hm1 = 0.5D+00 * q / ( m * m - 1.0D+00 )
      hm3 = 0.25D+00 * hm1 ** 3 / ( m * m - 4.0D+00 )
      hm5 = hm1 * hm3 * q 
     &  / ( ( m * m - 1.0D+00 ) * ( m * m - 9.0D+00 ) )
      a0 = m * m + q * ( hm1 + ( 5.0D+00 * m * m + 7.0D+00 ) * hm3
     &  + ( 9.0D+00 * m ** 4 + 58.0D+00 * m * m + 29.0D+00 ) * hm5 )

      return
      end
      subroutine cy01 ( kf, z, zf, zd )

c*********************************************************************72
c
cc CY01 computes complex Bessel functions Y0(z) and Y1(z) and derivatives.
c
c  Licensing:
c
c    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
c    they give permission to incorporate this routine into a user program 
c    provided that the copyright is acknowledged.
c
c  Modified:
c
c    01 August 2012
c
c  Author:
c
c    Shanjie Zhang, Jianming Jin
c
c  Reference:
c
c    Shanjie Zhang, Jianming Jin,
c    Computation of Special Functions,
c    Wiley, 1996,
c    ISBN: 0-471-11963-6,
c    LC: QA351.C45.
c
c  Parameters:
c
c    Input, integer KF, the function choice.
c    0 for ZF = Y0(z) and ZD = Y0'(z);
c    1 for ZF = Y1(z) and ZD = Y1'(z);
c    2 for ZF = Y1'(z) and ZD = Y1''(z).
c
c    Input, complex*16 Z, the argument.
c
c    Output, complex*16 ZF, ZD, the values of the requested function 
c    and derivative.
c
      implicit none

      double precision a(12)
      double precision a0
      double precision b(12)
      double precision a1(12)
      double precision b1(12)
      complex*16 cbj0
      complex*16 cbj1
      complex*16 cby0
      complex*16 cby1
      complex*16 cdy0
      complex*16 cdy1
      complex*16 ci
      complex*16 cp
      complex*16 cp0
      complex*16 cp1
      complex*16 cq0
      complex*16 cq1
      complex*16 cr
      complex*16 cs
      complex*16 ct1
      complex*16 ct2
      complex*16 cu
      double precision el
      integer k
      integer k0
      integer kf
      double precision pi
      double precision rp2
      double precision w0
      double precision w1
      complex*16 z
      complex*16 z1
      complex*16 z2
      complex*16 zd
      complex*16 zf

      save a
      save a1
      save b
      save b1

      data a / -0.703125D-01, 0.112152099609375D+00,
     &            -0.5725014209747314D+00, 0.6074042001273483D+01,
     &            -0.1100171402692467D+03, 0.3038090510922384D+04,
     &            -0.1188384262567832D+06, 0.6252951493434797D+07,
     &            -0.4259392165047669D+09, 0.3646840080706556D+11,
     &            -0.3833534661393944D+13, 0.4854014686852901D+15/

      data a1 / 0.1171875D+00, -0.144195556640625D+00,
     &             0.6765925884246826D+00, -0.6883914268109947D+01,
     &             0.1215978918765359D+03, -0.3302272294480852D+04,
     &             0.1276412726461746D+06, -0.6656367718817688D+07,
     &             0.4502786003050393D+09, -0.3833857520742790D+11,
     &             0.4011838599133198D+13, -0.5060568503314727D+15/

      data b / 0.732421875D-01, -0.2271080017089844D+00,
     &             0.1727727502584457D+01, -0.2438052969955606D+02,
     &             0.5513358961220206D+03, -0.1825775547429318D+05,
     &             0.8328593040162893D+06, -0.5006958953198893D+08,
     &             0.3836255180230433D+10, -0.3649010818849833D+12,
     &             0.4218971570284096D+14, -0.5827244631566907D+16/

      data b1 / -0.1025390625D+00, 0.2775764465332031D+00,
     &             -0.1993531733751297D+01, 0.2724882731126854D+02,
     &             -0.6038440767050702D+03, 0.1971837591223663D+05,
     &             -0.8902978767070678D+06, 0.5310411010968522D+08,
     &             -0.4043620325107754D+10, 0.3827011346598605D+12,
     &             -0.4406481417852278D+14, 0.6065091351222699D+16/

      pi = 3.141592653589793D+00
      el = 0.5772156649015329D+00
      rp2 = 2.0D+00 / pi
      ci = cmplx ( 0.0D+00, 1.0D+00 )
      a0 = cdabs ( z )
      z2 = z * z
      z1 = z

      if ( a0 .eq. 0.0D+00 ) then
        cbj0 = cmplx ( 1.0D+00, 0.0D+00 )
        cbj1 = cmplx ( 0.0D+00, 0.0D+00 )
        cby0 = cmplx ( -1.0D+30, 0.0D+00 )
        cby1 = cmplx ( -1.0D+30, 0.0D+00 )
        cdy0 = cmplx ( 1.0D+30, 0.0D+00 )
        cdy1 = cmplx ( 1.0D+30, 0.0D+00 )
        go to 50
      end if

      if ( real ( z ) .lt. 0.0D+00) then
        z1 = -z
      end if

      if ( a0 .le. 12.0D+00 ) then

        cbj0 = cmplx ( 1.0D+00, 0.0D+00 )
        cr = cmplx ( 1.0D+00, 0.0D+00 )
        do k = 1, 40
          cr = - 0.25D+00 * cr * z2 / ( k * k )
          cbj0 = cbj0 + cr
          if ( cdabs ( cr ) .lt. cdabs ( cbj0 ) * 1.0D-15 ) then
            go to 10
          end if
        end do

10      continue

        cbj1 = cmplx ( 1.0D+00, 0.0D+00 )
        cr = cmplx ( 1.0D+00, 0.0D+00 )
        do k = 1, 40
          cr = -0.25D+00 * cr * z2 / ( k * ( k + 1.0D+00 ) )
          cbj1 = cbj1 + cr
          if ( cdabs ( cr ) .lt. cdabs ( cbj1 ) * 1.0D-15 ) then
            go to 20
          end if
        end do

20      continue

        cbj1 = 0.5D+00 * z1 * cbj1
        w0 = 0.0D+00
        cr = cmplx ( 1.0D+00, 0.0D+00 )
        cs = cmplx ( 0.0D+00, 0.0D+00 )
        do k = 1, 40
          w0 = w0 + 1.0D+00 / k
          cr = -0.25D+00 * cr / ( k * k ) * z2
          cp = cr * w0
          cs = cs + cp
          if ( cdabs ( cp ) .lt. cdabs ( cs ) * 1.0D-15 ) then
            go to 30
          end if
        end do

30      continue

        cby0 = rp2 * ( cdlog ( z1 / 2.0D+00 ) + el ) * cbj0 - rp2 * cs
        w1 = 0.0D+00
        cr = cmplx ( 1.0D+00, 0.0D+00 )
        cs = cmplx ( 1.0D+00, 0.0D+00 )
        do k = 1, 40
          w1 = w1 + 1.0D+00 / k
          cr = - 0.25D+00 * cr / ( k * ( k + 1 ) ) * z2
          cp = cr * ( 2.0D+00 * w1 + 1.0D+00 / ( k + 1.0D+00 ) )
          cs = cs + cp
          if ( cdabs ( cp ) .lt. cdabs ( cs ) * 1.0D-15 ) then
            go to 40
          end if
        end do

40      continue

        cby1 = rp2 * ( ( cdlog ( z1 / 2.0D+00 ) + el ) * cbj1
     &    - 1.0D+00 / z1 - 0.25D+00 * z1 * cs )

      else

        if ( a0 .lt. 35.0D+00 ) then
          k0 = 12
        else if ( a0 .lt. 50.0D+00 ) then
          k0 = 10
        else
          k0 = 8
        end if

        ct1 = z1 - 0.25D+00 * pi
        cp0 = cmplx ( 1.0D+00, 0.0D+00 )
        do k = 1, k0
          cp0 = cp0 + a(k) * z1 ** ( - 2 * k )
        end do
        cq0 = -0.125D+00 / z1
        do k = 1, k0
          cq0 = cq0 + b(k) * z1 ** ( - 2 * k - 1 )
        end do
        cu = cdsqrt ( rp2 / z1 )
        cbj0 = cu * ( cp0 * cdcos ( ct1 ) - cq0 * cdsin ( ct1 ) )
        cby0 = cu * ( cp0 * cdsin ( ct1 ) + cq0 * cdcos ( ct1 ) )
        ct2 = z1 - 0.75D+00 * pi
        cp1 = cmplx ( 1.0D+00, 0.0D+00 )
        do k = 1, k0
          cp1 = cp1 + a1(k) * z1 ** ( - 2 * k )
        end do
        cq1 = 0.375D+00 / z1
        do k = 1, k0
          cq1 = cq1 + b1(k) * z1 ** ( - 2 * k - 1 )
        end do
        cbj1 = cu * ( cp1 * cdcos ( ct2 ) - cq1 * cdsin ( ct2 ) )
        cby1 = cu * ( cp1 * cdsin ( ct2 ) + cq1 * cdcos ( ct2 ) )

      end if

      if ( real ( z ) .lt. 0.0D+00 ) then

        if ( dimag ( z ) .lt. 0.0D+00 ) then
          cby0 = cby0 - 2.0D+00 * ci * cbj0
        else
          cby0 = cby0 + 2.0D+00 * ci * cbj0
        end if

        if ( dimag ( z ) .lt. 0.0D+00 ) then
          cby1 = - ( cby1 - 2.0D+00 * ci * cbj1 )
        else
          cby1 = - ( cby1 + 2.0D+00 * ci * cbj1 )
        end if
        cbj1 = - cbj1

      end if

      cdy0 = - cby1
      cdy1 = cby0 - 1.0D+00 / z * cby1

50    continue

      if ( kf .eq. 0 ) then
        zf = cby0
        zd = cdy0
      else if ( kf .eq. 1 ) then
        zf = cby1
        zd = cdy1
      else if ( kf .eq. 2 ) then
        zf = cdy1
        zd = - cdy1 / z - ( 1.0D+00 - 1.0D+00 / ( z * z ) ) * cby1
      end if

      return
      end
      subroutine cyzo ( nt, kf, kc, zo, zv )

c*********************************************************************72
c
cc CYZO computes zeros of complex Bessel functions Y0(z) and Y1(z) and Y1'(z).
c
c  Parameters:
c
c    Ths procedure computes the complex zeros of Y0(z), Y1(z) and Y1'(z), 
c    and their associated values at the zeros using the modified Newton's 
c    iteration method.
c
c  Licensing:
c
c    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
c    they give permission to incorporate this routine into a user program 
c    provided that the copyright is acknowledged.
c
c  Modified:
c
c    22 July 2012
c
c  Author:
c
c    Shanjie Zhang, Jianming Jin
c
c  Reference:
c
c    Shanjie Zhang, Jianming Jin,
c    Computation of Special Functions,
c    Wiley, 1996,
c    ISBN: 0-471-11963-6,
c    LC: QA351.C45.
c
c  Parameters:
c
c    Input, integer NT, the number of zeros.
c
c    Input, integer KF, the function choice.
c    0 for Y0(z) and Y1(z0);
c    1 for Y1(z) and Y0(z1);
c    2 for Y1'(z) and Y1(z1').
c
c    Input, integer KC, complex/real choice.
c    0, for complex roots;
c    1, for real roots.
c
c    Output, double precision ZO(NT), ZV(NT), the zeros of Y0(z) or Y1(z) 
c    or Y1'(z), and the value of Y0'(z) or Y1'(z) or Y1(z) at the L-th zero.
c
      implicit none

      integer nt

      double precision h
      integer i
      integer it
      integer j
      integer kc
      integer kf
      integer nr
      double precision w
      double precision w0
      double precision x
      double precision y
      complex*16 z
      complex*16 zd
      complex*16 zero
      complex*16 zf
      complex*16 zfd
      complex*16 zgd
      complex*16 zo(nt)
      complex*16 zp
      complex*16 zq
      complex*16 zv(nt)
      complex*16 zw

      if ( kc .eq. 0 ) then
        x = -2.4D+00
        y = 0.54D+00
        h = 3.14D+00
      else if ( kc .eq. 1 ) then
        x = 0.89D+00
        y = 0.0D+00
        h = -3.14D+00
      end if

      if ( kf .eq. 1 ) then
        x = -0.503D+00
      else if ( kf .eq. 2 ) then
        x = 0.577D+00
      end if

      zero = cmplx ( x, y )

      do nr = 1, nt

        if ( nr .eq. 1 ) then
          z = zero
        else
          z = zo(nr-1) - h
        end if

        it = 0

10      continue

        it = it + 1
        call cy01 ( kf, z, zf, zd )

        zp = cmplx ( 1.0D+00, 0.0D+00 )
        do i = 1, nr - 1
          zp = zp * ( z - zo(i) )
        end do

        zfd = zf / zp

        zq = cmplx ( 0.0D+00, 0.0D+00 )
        do i = 1, nr - 1
          zw = cmplx ( 1.0D+00, 0.0D+00 )
          do j = 1, nr - 1
            if ( j .ne. i ) then
              zw = zw * ( z - zo(j) )
            end if
          end do
          zq = zq + zw
        end do

        zgd = ( zd - zq * zfd ) / zp
        z = z - zfd / zgd
        w0 = w
        w = cdabs ( z )

        if ( it .le. 50 .and. 1.0D-12 .lt. abs ( ( w - w0 ) / w ) ) then
          go to 10
        end if

        zo(nr) = z

      end do

      do i = 1, nt
        z = zo(i)
        if ( kf .eq. 0 .or. kf .eq. 2 ) then
          call cy01 ( 1, z, zf, zd )
          zv(i) = zf
        else if ( kf .eq. 1 ) then
          call cy01 ( 0, z, zf, zd )
          zv(i) = zf
        end if
      end do

      return
      end
      subroutine dvla ( va, x, pd )

c*********************************************************************72
c
cc DVLA computes parabolic cylinder functions Dv(x) for large argument.
c
c  Licensing:
c
c    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
c    they give permission to incorporate this routine into a user program 
c    provided that the copyright is acknowledged.
c
c  Modified:
c
c    06 July 2012
c
c  Author:
c
c    Shanjie Zhang, Jianming Jin
c
c  Reference:
c
c    Shanjie Zhang, Jianming Jin,
c    Computation of Special Functions,
c    Wiley, 1996,
c    ISBN: 0-471-11963-6,
c    LC: QA351.C45.
c 
c  Parameters:
c
c    Input, double precision VA, the order.
c
c    Input, double precision X, the argument.
c
c    Output, double precision PD, the function value.
c
      implicit none

      double precision a0
      double precision ep
      double precision eps
      double precision gl
      integer k
      double precision pd
      double precision pi
      double precision r
      double precision va
      double precision vl
      double precision x
      double precision x1

      pi = 3.141592653589793D+00
      eps = 1.0D-12
      ep = exp ( -0.25D+00 * x * x )
      a0 = abs ( x ) ** va * ep
      r = 1.0D+00
      pd = 1.0D+00
      do k = 1, 16
        r = -0.5D+00 * r * ( 2.0D+00 * k - va - 1.0D+00 ) 
     &    * ( 2.0D+00 * k - va - 2.0D+00 ) / ( k * x * x )
        pd = pd + r
        if ( abs ( r / pd ) .lt. eps ) then
          go to 10
        end if
      end do

10    continue

      pd = a0 * pd

      if ( x .lt. 0.0D+00 ) then
        x1 = - x
        call vvla ( va, x1, vl )
        call gamma ( -va, gl )
        pd = pi * vl / gl + cos ( pi * va ) * pd
      end if

      return
      end
      subroutine dvsa ( va, x, pd )

c*********************************************************************72
c
cc DVSA computes parabolic cylinder functions Dv(x) for small argument.
c
c  Licensing:
c
c    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
c    they give permission to incorporate this routine into a user program 
c    provided that the copyright is acknowledged.
c
c  Modified:
c
c    07 July 2012
c
c  Author:
c
c    Shanjie Zhang, Jianming Jin
c
c  Reference:
c
c    Shanjie Zhang, Jianming Jin,
c    Computation of Special Functions,
c    Wiley, 1996,
c    ISBN: 0-471-11963-6,
c    LC: QA351.C45.
c
c  Parameters:
c
c    Input, double precision VA, the order.
c
c    Input, double precision X, the argument.
c
c    Output, double precision PD, the function value.
c
      implicit none

      double precision a0
      double precision ep
      double precision eps
      double precision g0
      double precision g1
      double precision ga0
      double precision gm
      integer m
      double precision pd
      double precision pi
      double precision r
      double precision r1
      double precision sq2
      double precision va
      double precision va0
      double precision vm
      double precision vt
      double precision x

      eps = 1.0D-15
      pi = 3.141592653589793D+00
      sq2 = sqrt ( 2.0D+00 )
      ep = exp ( -0.25D+00 * x * x )
      va0 = 0.5D+00 * ( 1.0D+00 - va )

      if ( va .eq. 0.0D+00 ) then

        pd = ep

      else

        if ( x .eq. 0.0D+00 ) then
          if ( va0 .le. 0.0D+00 .and. va0 .eq. int ( va0 ) ) then
            pd = 0.0D+00
          else
            call gamma ( va0, ga0 )
            pd = sqrt ( pi ) / ( 2.0D+00 ** ( - 0.5D+00 * va ) * ga0 )
          end if

        else

          call gamma ( -va, g1 )
          a0 = 2.0D+00 ** ( -0.5D+00 * va - 1.0D+00 ) * ep / g1
          vt = -0.5D+00 * va
          call gamma ( vt, g0 )
          pd = g0
          r = 1.0D+00
          do m = 1, 250
            vm = 0.5D+00 * ( m - va )
            call gamma ( vm, gm )
            r = - r * sq2 * x / m
            r1 = gm * r
            pd = pd + r1
            if ( abs ( r1 ) .lt. abs ( pd ) * eps ) then
              go to 10
            end if
          end do

10        continue

          pd = a0 * pd

        end if

      end if

      return
      end
      subroutine e1xa ( x, e1 )

c*********************************************************************72
c
cc E1XA computes the exponential integral E1(x).
c
c  Licensing:
c
c    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
c    they give permission to incorporate this routine into a user program 
c    provided that the copyright is acknowledged.
c
c  Modified:
c
c    06 July 2012
c
c  Author:
c
c    Shanjie Zhang, Jianming Jin
c
c  Reference:
c
c    Shanjie Zhang, Jianming Jin,
c    Computation of Special Functions,
c    Wiley, 1996,
c    ISBN: 0-471-11963-6,
c    LC: QA351.C45.
c
c  Parameters:
c
c    Input, double precision X, the argument.
c
c    Output, double precision E1, the function value.
c
      implicit none

      double precision e1
      double precision es1
      double precision es2
      double precision x

      if ( x .eq. 0.0D+00 ) then

        e1 = 1.0D+300

      else if ( x .le. 1.0D+00 ) then

        e1 = -log ( x ) + ((((
     &      1.07857D-03 * x
     &    - 9.76004D-03 ) * x
     &    + 5.519968D-02 ) * x
     &    - 0.24991055D+00 ) * x
     &    + 0.99999193D+00 ) * x
     &    - 0.57721566D+00

      else

        es1 = ((( x
     &    +  8.5733287401D+00 ) * x 
     &    + 18.059016973D+00 )  * x
     &    +  8.6347608925D+00 ) * x 
     &    +  0.2677737343D+00

        es2 = ((( x
     &    +  9.5733223454D+00 ) * x
     &    + 25.6329561486D+00 ) * x
     &    + 21.0996530827D+00 ) * x
     &    +  3.9584969228D+00

        e1 = exp ( - x ) / x * es1 / es2

      end if

      return
      end
      subroutine e1xb ( x, e1 )

c*********************************************************************72
c
cc E1XB computes the exponential integral E1(x).
c
c  Licensing:
c
c    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
c    they give permission to incorporate this routine into a user program 
c    provided that the copyright is acknowledged.
c
c  Modified:
c
c    06 July 2012
c
c  Author:
c
c    Shanjie Zhang, Jianming Jin
c
c  Reference:
c
c    Shanjie Zhang, Jianming Jin,
c    Computation of Special Functions,
c    Wiley, 1996,
c    ISBN: 0-471-11963-6,
c    LC: QA351.C45.
c
c  Parameters:
c
c    Input, double precision X, the argument.
c
c    Output, double precision E1, the function value.
c
      implicit none

      double precision e1
      double precision ga
      integer k
      integer m
      double precision r
      double precision t
      double precision t0
      double precision x

      if ( x .eq. 0.0D+00 ) then

        e1 = 1.0D+300

      else if ( x .le. 1.0D+00 ) then

        e1 = 1.0D+00
        r = 1.0D+00

        do k = 1, 25
          r = -r * k * x / ( k + 1.0D+00 )**2
          e1 = e1 + r
          if ( abs ( r ) .le. abs ( e1 ) * 1.0D-15 ) then
            go to 10
          end if
        end do

10      continue

        ga = 0.5772156649015328D+00
        e1 = - ga - log ( x ) + x * e1

      else

        m = 20 + int ( 80.0D+00 / x )
        t0 = 0.0D+00
        do k = m, 1, -1
          t0 = k / ( 1.0D+00 + k / ( x + t0 ) )
        end do
        t = 1.0D+00 / ( x + t0 )
        e1 = exp ( -x ) * t

      end if

      return
      end
      subroutine e1z ( z, ce1 )

c*********************************************************************72
c
cc E1Z computes the complex exponential integral E1(z).
c
c  Licensing:
c
c    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
c    they give permission to incorporate this routine into a user program 
c    provided that the copyright is acknowledged.
c
c  Modified:
c
c    16 July 2012
c
c  Author:
c
c    Shanjie Zhang, Jianming Jin
c
c  Reference:
c
c    Shanjie Zhang, Jianming Jin,
c    Computation of Special Functions,
c    Wiley, 1996,
c    ISBN: 0-471-11963-6,
c    LC: QA351.C45.
c 
c  Parameters:
c
c    Input, complex*16 Z, the argument.
c
c    Output, complex*16 CE1, the function value.
c
      implicit none

      double precision a0
      complex*16 ce1
      complex*16 cr
      complex*16 ct
      complex*16 ct0
      double precision el
      integer k
      double precision pi
      double precision x
      complex*16 z

      pi = 3.141592653589793D+00
      el = 0.5772156649015328D+00
      x = real ( z )
      a0 = cdabs ( z )

      if ( a0 .eq. 0.0D+00 ) then
        ce1 = cmplx ( 1.0D+30, 0.0D+00 )
      else if ( a0 .le. 10.0D+00 .or. 
     &  ( x .lt. 0.0D+00 .and. a0 .lt. 20.0D+00 ) ) then
        ce1 = cmplx ( 1.0D+00, 0.0D+00 )
        cr = cmplx ( 1.0D+00, 0.0D+00 )
        do k = 1, 150
          cr = - cr * k * z / ( k + 1.0D+00 )**2
          ce1 = ce1 + cr
          if ( cdabs ( cr ) .le. cdabs ( ce1 ) * 1.0D-15 ) then
            go to 10
          end if
        end do

10      continue

        ce1 = - el - cdlog ( z ) + z * ce1

      else

        ct0 = cmplx ( 0.0D+00, 0.0D+00 )
        do k = 120, 1, -1
          ct0 = k / ( 1.0D+00 + k / ( z + ct0 ) )
        end do
        ct = 1.0D+00 / ( z + ct0 )

        ce1 = cdexp ( - z ) * ct
        if ( x .le. 0.0D+00 .and. dimag ( z ) .eq. 0.0D+00 ) then
          ce1 = ce1 - pi * cmplx ( 0.0D+00, 1.0D+00 )
        end if

      end if

      return
      end
      subroutine eix ( x, ei )

c*********************************************************************72
c
cc EIX computes the exponential integral Ei(x).
c
c  Licensing:
c
c    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
c    they give permission to incorporate this routine into a user program 
c    provided that the copyright is acknowledged.
c
c  Modified:
c
c    10 July 2012
c
c  Author:
c
c    Shanjie Zhang, Jianming Jin
c
c  Reference:
c
c    Shanjie Zhang, Jianming Jin,
c    Computation of Special Functions,
c    Wiley, 1996,
c    ISBN: 0-471-11963-6,
c    LC: QA351.C45.
c
c  Parameters:
c
c    Input, double precision X, the argument.
c
c    Output, double precision EI, the function value.
c
      implicit none

      double precision ei
      double precision ga
      integer k
      double precision r
      double precision x

      if ( x .eq. 0.0D+00 ) then

        ei = -1.0D+300

      else if ( x .le. 40.0D+00 ) then

        ei = 1.0D+00
        r = 1.0D+00
        do k = 1, 100
          r = r * k * x / ( k + 1.0D+00 )**2
          ei = ei + r
          if ( abs ( r / ei ) .le. 1.0D-15 ) then
            go to 10
          end if
        end do

10      continue

        ga = 0.5772156649015328D+00
        ei = ga + log ( x ) + x * ei

      else

        ei = 1.0D+00
        r = 1.0D+00
        do k = 1, 20
          r = r * k / x
          ei = ei + r
        end do
        ei = exp ( x ) / x * ei

      end if

      return
      end
      subroutine elit ( hk, phi, fe, ee )

c*********************************************************************72
c
cc ELIT computes complete and incomplete elliptic integrals F(k,phi) and E(k,phi).
c
c  Licensing:
c
c    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
c    they give permission to incorporate this routine into a user program 
c    provided that the copyright is acknowledged.
c
c  Modified:
c
c    12 July 2012
c
c  Author:
c
c    Shanjie Zhang, Jianming Jin
c
c  Reference:
c
c    Shanjie Zhang, Jianming Jin,
c    Computation of Special Functions,
c    Wiley, 1996,
c    ISBN: 0-471-11963-6,
c    LC: QA351.C45.
c
c  Parameters:
c
c    Input, double precision HK, the modulus, between 0 and 1.
c
c    Input, double precision PHI, the argument in degrees.
c
c    Output, double precision FE, EE, the values of F(k,phi) and E(k,phi).
c
      implicit none

      double precision a
      double precision a0
      double precision b
      double precision b0
      double precision c
      double precision ce
      double precision ck
      double precision d
      double precision d0
      double precision ee
      double precision fac
      double precision fe
      double precision g
      double precision hk
      integer n
      double precision phi
      double precision pi
      double precision r

      g = 0.0D+00
      pi = 3.14159265358979D+00
      a0 = 1.0D+00
      b0 = sqrt ( 1.0D+00 - hk * hk )
      d0 = ( pi / 180.0D+00 ) * phi
      r = hk * hk

      if ( hk .eq. 1.0D+00 .and. phi .eq. 90.0D+00 ) then

        fe = 1.0D+300
        ee = 1.0D+00

      else if ( hk .eq. 1.0D+00 ) then

        fe = log ( ( 1.0D+00 + sin ( d0 ) ) / cos ( d0 ) )
        ee = sin ( d0 )

      else

        fac = 1.0D+00
        do n = 1, 40
          a = ( a0 + b0 ) /2.0D+00
          b = sqrt ( a0 * b0 )
          c = ( a0 - b0 ) / 2.0D+00
          fac = 2.0D+00 * fac
          r = r + fac * c * c
          if ( phi .ne. 90.0D+00 ) then
            d = d0 + atan ( ( b0 / a0 ) * tan ( d0 ) )
            g = g + c * sin( d )
            d0 = d + pi * int ( d / pi + 0.5D+00 )
          end if
          a0 = a
          b0 = b
          if ( c .lt. 1.0D-07 ) then
            go to 10
          end if
        end do

10      continue

        ck = pi / ( 2.0D+00 * a )
        ce = pi * ( 2.0D+00 - r ) / ( 4.0D+00 * a )
        if ( phi .eq. 90.0D+00 ) then
          fe = ck
          ee = ce
        else
          fe = d / ( fac * a )
          ee = fe * ce / ck + g
        end if

      end if

      return
      end
      subroutine elit3 ( phi, hk, c, el3 )

c*********************************************************************72
c
cc ELIT3 computes the elliptic integral of the third kind.
c
c  Discussion:
c
c    Gauss-Legendre quadrature is used.
c
c  Licensing:
c
c    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
c    they give permission to incorporate this routine into a user program 
c    provided that the copyright is acknowledged.
c
c  Modified:
c
c    14 July 2012
c
c  Author:
c
c    Shanjie Zhang, Jianming Jin
c
c  Reference:
c
c    Shanjie Zhang, Jianming Jin,
c    Computation of Special Functions,
c    Wiley, 1996,
c    ISBN: 0-471-11963-6,
c    LC: QA351.C45.
c
c  Parameters:
c
c    Input, double precision PHI, the argument in degrees.
c
c    Input, double precision HK, the modulus, between 0 and 1.
c
c    Input, double precision C, the parameter, between 0 and 1.
c
c    Output, double precision EL3, the value of the elliptic integral
c    of the third kind.
c
      implicit none

      double precision c
      double precision c0
      double precision c1
      double precision c2
      double precision el3
      double precision f1
      double precision f2
      double precision hk
      integer i 
      logical lb1
      logical lb2
      double precision phi
      double precision t(10)
      double precision t1
      double precision t2
      double precision w(10)

      save t
      save w

      data t / 0.9931285991850949D+00, 0.9639719272779138D+00,
     &         0.9122344282513259D+00, 0.8391169718222188D+00,
     &         0.7463319064601508D+00, 0.6360536807265150D+00,
     &         0.5108670019508271D+00, 0.3737060887154195D+00,
     &         0.2277858511416451D+00, 0.7652652113349734D-01/

      data w / 0.1761400713915212D-01, 0.4060142980038694D-01,
     &         0.6267204833410907D-01, 0.8327674157670475D-01,
     &         0.1019301198172404D+00, 0.1181945319615184D+00,
     &         0.1316886384491766D+00, 0.1420961093183820D+00,
     &         0.1491729864726037D+00, 0.1527533871307258D+00 /

      lb1 = ( hk .eq. 1.0D+00 ) .and.
     &  ( abs ( phi - 90.0D+00 ) .le. 1.0D-08 )

      lb2 = c .eq. 1.0D+00 .and. abs ( phi - 90.0D+00 ) .le. 1.0D-08

      if ( lb1 .or. lb2 ) then
        el3 = 1.0D+300
        return
      end if

      c1 = 0.87266462599716D-02 * phi
      c2 = c1

      el3 = 0.0D+00
      do i = 1, 10
        c0 = c2 * t(i)
        t1 = c1 + c0
        t2 = c1 - c0
        f1 = 1.0D+00 / ( ( 1.0D+00 - c * sin(t1) * sin(t1) ) 
     &    * sqrt ( 1.0D+00 - hk * hk * sin ( t1 ) * sin ( t1 ) ) )
        f2 = 1.0D+00 / ( ( 1.0D+00 - c * sin ( t2 ) * sin ( t2 ) ) *
     &    sqrt( 1.0D+00 - hk * hk * sin ( t2 ) * sin ( t2 ) ) )
        el3 = el3 + w(i) * ( f1 + f2 )
      end do

      el3 = c1 * el3

      return
      end
      function envj ( n, x )

c*********************************************************************72
c
cc ENVJ is a utility function used by MSTA1 and MSTA2.
c
c  Licensing:
c
c    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
c    they give permission to incorporate this routine into a user program 
c    provided that the copyright is acknowledged.
c
c  Modified:
c
c    14 March 2012
c
c  Author:
c
c    Shanjie Zhang, Jianming Jin
c
c  Reference:
c
c    Shanjie Zhang, Jianming Jin,
c    Computation of Special Functions,
c    Wiley, 1996,
c    ISBN: 0-471-11963-6,
c    LC: QA351.C45.
c
c  Parameters:
c
c    Input, integer N, ?
c
c    Input, double precision X, ?
c
c    Output, double precision ENVJ, ?
c
      implicit none

      double precision envj
      integer n
      double precision x

      envj = 0.5D+00 * log10 ( 6.28D+00 * n ) 
     &  - n * log10 ( 1.36D+00 * x / n )

      return
      end
      subroutine enxa ( n, x, en )

c*********************************************************************72
c
cc ENXA computes the exponential integral En(x).
c
c  Licensing:
c
c    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
c    they give permission to incorporate this routine into a user program 
c    provided that the copyright is acknowledged.
c
c  Modified:
c
c    07 July 2012
c
c  Author:
c
c    Shanjie Zhang, Jianming Jin
c
c  Reference:
c
c    Shanjie Zhang, Jianming Jin,
c    Computation of Special Functions,
c    Wiley, 1996,
c    ISBN: 0-471-11963-6,
c    LC: QA351.C45.
c
c  Parameters:
c
c    Input, integer N, the order.
c
c    Input, double precision X, the argument.
c
c    Output, double precision EN(0:N), the function values.
c
      implicit none

      integer n

      double precision e1
      double precision ek
      double precision en(0:n)
      integer k
      double precision x

      en(0) = exp ( - x ) / x 
      call e1xb ( x, e1 )

      en(1) = e1
      do k = 2, n
        ek = ( exp ( - x ) - x * e1 ) / ( k - 1.0D+00 )
        en(k) = ek
        e1 = ek
      end do

      return
      end
      subroutine enxb ( n, x, en )

c*********************************************************************72
c
cc ENXB computes the exponential integral En(x).
c
c  Licensing:
c
c    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
c    they give permission to incorporate this routine into a user program 
c    provided that the copyright is acknowledged.
c
c  Modified:
c
c    10 July 2012
c
c  Author:
c
c    Shanjie Zhang, Jianming Jin
c
c  Reference:
c
c    Shanjie Zhang, Jianming Jin,
c    Computation of Special Functions,
c    Wiley, 1996,
c    ISBN: 0-471-11963-6,
c    LC: QA351.C45.
c
c  Parameters:
c
c    Input, integer N, the order.
c
c    Input, double precision X, the argument.
c
c    Output, double precision EN(0:N), the function values.
c
      implicit none

      integer n

      double precision en(0:n)
      double precision ens
      integer j
      integer k
      integer l
      integer m
      double precision ps
      double precision r
      double precision rp
      double precision s
      double precision s0
      double precision t
      double precision t0
      double precision x

      if ( x .eq. 0.0D+00 ) then

        en(0) = 1.0D+300
        en(1) = 1.0D+300
        do k = 2, n
          en(k) = 1.0D+00 / ( k - 1.0D+00 )
        end do
        return

      else if ( x .le. 1.0D+00 ) then

        en(0) = exp ( - x ) / x
        do l = 1, n
          rp = 1.0D+00
          do j = 1, l - 1
            rp = - rp * x / j
          end do
          ps = -0.5772156649015328D+00
          do m = 1, l - 1
            ps = ps + 1.0D+00 / m
          end do
          ens = rp * ( - log ( x ) + ps )
          s = 0.0D+00
          do m = 0, 20
            if ( m .ne. l - 1 ) then
              r = 1.0D+00
              do j = 1, m
                r = - r * x / j
              end do
              s = s + r / ( m - l + 1.0D+00 )
              if ( abs ( s - s0 ) .lt. abs ( s ) * 1.0D-15 ) then
                go to 10
              end if
              s0 = s
            end if
          end do

10        continue

          en(l) = ens - s

        end do

      else

        en(0) = exp ( - x ) / x
        m = 15 + int ( 100.0D+00 / x )
        do l = 1, n
          t0 = 0.0D+00
          do k = m, 1, -1
            t0 = ( l + k - 1.0D+00 ) / ( 1.0D+00 + k / ( x + t0 ) )
          end do
          t = 1.0D+00 / ( x + t0 )
          en(l) = exp ( - x ) * t
        end do

      end if

      return
      end
      subroutine error ( x, err )

c*********************************************************************72
c
cc ERROR evaluates the error function.
c
c  Licensing:
c
c    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
c    they give permission to incorporate this routine into a user program 
c    provided that the copyright is acknowledged.
c
c  Modified:
c
c    10 July 2012
c
c  Author:
c
c    Shanjie Zhang, Jianming Jin
c
c  Reference:
c
c    Shanjie Zhang, Jianming Jin,
c    Computation of Special Functions,
c    Wiley, 1996,
c    ISBN: 0-471-11963-6,
c    LC: QA351.C45.
c
c  Parameters:
c
c    Input, double precision X, the argument.
c
c    Output, double precision ERR, the function value.
c
      implicit none

      double precision c0
      double precision eps
      double precision er
      double precision err
      integer k
      double precision pi
      double precision r
      double precision x
      double precision x2

      eps = 1.0D-15
      pi = 3.141592653589793D+00
      x2 = x * x

      if ( abs ( x ) .lt. 3.5D+00 ) then

        er = 1.0D+00
        r = 1.0D+00

        do k = 1, 50
          r = r * x2 / ( k + 0.5D+00 )
          er = er+r
          if ( abs ( r ) .le. abs ( er ) * eps ) then
            go to 10
          end if
        end do

10      continue

        c0 = 2.0D+00 / sqrt ( pi ) * x * exp ( - x2 )
        err = c0 * er

      else

        er = 1.0D+00
        r = 1.0D+00
        do k = 1, 12
          r = - r * ( k - 0.5D+00 ) / x2
          er = er + r
        end do

        c0 = exp ( - x2 ) / ( abs ( x ) * sqrt ( pi ) )

        err = 1.0D+00 - c0 * er
        if ( x .lt. 0.0D+00 ) then
          err = -err
        end if

      end if

      return
      end
      subroutine eulera ( n, en )

c*********************************************************************72
c
cc EULERA computes the Euler number En.
c
c  Licensing:
c
c    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
c    they give permission to incorporate this routine into a user program 
c    provided that the copyright is acknowledged.
c
c  Modified:
c
c    10 July 2012
c
c  Author:
c
c    Shanjie Zhang, Jianming Jin
c
c  Reference:
c
c    Shanjie Zhang, Jianming Jin,
c    Computation of Special Functions,
c    Wiley, 1996,
c    ISBN: 0-471-11963-6,
c    LC: QA351.C45.
c
c  Parameters:
c
c    Input, integer N, the index of the highest value to compute.
c
c    Output, double precision EN(0:N), the Euler numbers up to the N-th value.
c
      implicit none

      integer n

      double precision en(0:n)
      integer j
      integer k
      integer m
      double precision r
      double precision s

      en(0) = 1.0D+00

      do m = 1, n / 2
        s = 1.0D+00
        do k = 1, m - 1
          r = 1.0D+00
          do j = 1, 2 * k
            r = r * ( 2.0D+00 * m - 2.0D+00 * k + j ) / j
          end do
          s = s + r * en(2*k)
        end do
        en(2*m) = -s
      end do

      return
      end
      subroutine eulerb ( n, en )

c*********************************************************************72
c
cc EULERB computes the Euler number En.
c
c  Licensing:
c
c    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
c    they give permission to incorporate this routine into a user program 
c    provided that the copyright is acknowledged.
c
c  Modified:
c
c    09 July 2012
c
c  Author:
c
c    Shanjie Zhang, Jianming Jin
c
c  Reference:
c
c    Shanjie Zhang, Jianming Jin,
c    Computation of Special Functions,
c    Wiley, 1996,
c    ISBN: 0-471-11963-6,
c    LC: QA351.C45.
c
c  Parameters:
c
c    Input, integer N, the index of the highest value to compute.
c
c    Output, double precision EN(0:N), the Euler numbers up to the N-th value.
c
      implicit none

      integer n

      double precision en(0:n)
      double precision hpi
      double precision isgn
      integer k
      integer m
      double precision r1
      double precision r2
      double precision s

      hpi = 2.0D+00 / 3.141592653589793D+00
      en(0) = 1.0D+00
      en(2) = -1.0D+00
      r1 = -4.0D+00 * hpi ** 3

      do m = 4, n, 2
        r1 = - r1 * ( m - 1 ) * m * hpi * hpi
        r2 = 1.0D+00
        isgn = 1.0D+00
        do k = 3, 1000, 2
          isgn = - isgn
          s = ( 1.0D+00 / k ) ** ( m + 1 )
          r2 = r2 + isgn * s
          if ( s .lt. 1.0D-15 ) then
            go to 10
          end if
        end do

10      continue

        en(m) = r1 * r2

        end do

      return
      end
      subroutine fcoef ( kd, m, q, a, fc )

c*********************************************************************72
c
cc FCOEF: expansion coefficients for Mathieu and modified Mathieu functions.
c
c  Licensing:
c
c    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
c    they give permission to incorporate this routine into a user program 
c    provided that the copyright is acknowledged.
c
c  Modified:
c
c    01 August 2012
c
c  Author:
c
c    Shanjie Zhang, Jianming Jin
c
c  Reference:
c
c    Shanjie Zhang, Jianming Jin,
c    Computation of Special Functions,
c    Wiley, 1996,
c    ISBN: 0-471-11963-6,
c    LC: QA351.C45.
c
c  Parameters:
c
c    Input, integer KD, the case code.
c    1, for cem(x,q)  ( m = 0,2,4,...)
c    2, for cem(x,q)  ( m = 1,3,5,...)
c    3, for sem(x,q)  ( m = 1,3,5,...)
c    4, for sem(x,q)  ( m = 2,4,6,...)
c
c    Input, integer M, the order of the Mathieu function.
c
c    Input, double precision Q, the parameter of the Mathieu functions.
c
c    Input, double precision A, the characteristic value of the Mathieu
c    functions for given m and q.
c
c    Output, double precision FC(*), the expansion coefficients of Mathieu
c    functions ( k =  1,2,...,KM ).  FC(1),FC(2),FC(3),... correspond to
c    A0,A2,A4,... for KD = 1 case, 
c    A1,A3,A5,... for KD = 2 case,
c    B1,B3,B5,... for KD = 3 case,
c    B2,B4,B6,... for KD = 4 case.
c
      implicit none

      double precision a
      double precision f
      double precision f1
      double precision f2
      double precision f3
      double precision fc(251)
      integer i
      integer j
      integer k
      integer kb
      integer kd
      integer km
      integer m
      double precision q
      double precision qm
      double precision s
      double precision s0
      double precision sp
      double precision ss
      double precision u
      double precision v

      if ( q .le. 1.0D+00 ) then
        qm = 7.5D+00 + 56.1D+00 * sqrt ( q ) - 134.7D+00 * q 
     &    + 90.7D+00 * sqrt ( q ) * q
      else
        qm = 17.0D+00 + 3.1D+00 * sqrt ( q ) - 0.126D+00 * q 
     &    + 0.0037D+00 * sqrt ( q ) * q
      end if

      km = int ( qm + 0.5D+00 * m )

      if ( q .eq. 0.0D+00 ) then

        do k = 1, km
          fc(k) = 0.0D+00
        end do

        if ( kd .eq. 1 ) then
          fc((m+2)/2) = 1.0D+00
          if (m .eq. 0 ) then
            fc(1) = 1.0D+00 / sqrt ( 2.0D+00 )
          end if
        else if ( kd .eq. 4 ) then
          fc(m/2) = 1.0D+00
        else
          fc((m+1)/2) = 1.0D+00
        end if

        return

      end if

      kb = 0
      s = 0.0D+00
      f = 1.0D-100
      u = 0.0D+00
      fc(km) = 0.0D+00

      if ( kd .eq. 1 ) then

        do k = km, 3, -1

          v = u
          u = f
          f = ( a - 4.0D+00 * k * k ) * u / q - v

          if ( abs ( f ) .lt. abs ( fc(k+1) ) ) then

            kb = k
            fc(1) = 1.0D-100
            sp = 0.0D+00
            f3 = fc(k+1)
            fc(2) = a / q * fc(1)
            fc(3) = ( a - 4.0D+00 ) * fc(2) / q - 2.0D+00 * fc(1)
            u = fc(2)
            f1 = fc(3)

            do i = 3, kb
              v = u
              u = f1
              f1 = ( a - 4.0D+00 * ( i - 1.0D+00 ) ** 2 ) * u / q - v
              fc(i+1) = f1
              if ( i .eq. kb ) then
                f2 = f1
              else
                sp = sp + f1 * f1
              end if
            end do

            sp = sp + 2.0D+00 * fc(1) ** 2 + fc(2) ** 2 + fc(3) ** 2
            ss = s + sp * ( f3 / f2 ) ** 2
            s0 = sqrt ( 1.0D+00 / ss )
            do j = 1, km
              if ( j .le. kb + 1 ) then
                fc(j) = s0 * fc(j) * f3 / f2
              else
                fc(j) = s0 * fc(j)
              end if
            end do
            go to 30
          else
            fc(k) = f
            s = s + f * f
          end if

        end do

        fc(2) = q * fc(3) / ( a - 4.0D+00 - 2.0D+00 * q * q / a )
        fc(1) = q / a * fc(2)
        s = s + 2.0D+00 * fc(1) ** 2 + fc(2) ** 2
        s0 = sqrt ( 1.0D+00 / s )
        do k = 1, km
          fc(k) = s0 * fc(k)
        end do

      else if ( kd .eq. 2 .or. kd .eq. 3 ) then

        do k = km, 3, -1

          v = u
          u = f
          f = ( a - ( 2.0D+00 * k - 1 ) ** 2 ) * u / q - v

          if ( abs ( fc(k) ) .le. abs ( f ) ) then
            fc(k-1) = f
            s = s + f * f
          else
            kb = k
            f3 = fc(k)
            go to 10
          end if

        end do

        fc(1) = q / ( a - 1.0D+00 - ( - 1 ) ** kd * q ) * fc(2)
        s = s + fc(1) * fc(1)
        s0 = sqrt ( 1.0D+00 / s )
        do k = 1, km
          fc(k) = s0 * fc(k)
        end do

        go to 30

10      continue

        fc(1) = 1.0D-100
        fc(2) = ( a - 1.0D+00 - ( - 1 ) ** kd * q ) / q * fc(1)
        sp = 0.0D+00
        u = fc(1)
        f1 = fc(2)
        do i = 2, kb - 1
          v = u
          u = f1
          f1 = ( a - ( 2.0D+00 * i - 1.0D+00 ) ** 2 ) * u / q - v
          if ( i .ne. kb - 1 ) then
            fc(i+1) = f1
            sp = sp + f1 * f1
          else
            f2 = f1
          end if
        end do

        sp = sp + fc(1) ** 2 + fc(2) ** 2
        ss = s + sp * ( f3 / f2 ) ** 2
        s0 = 1.0D+00 / sqrt ( ss )
        do j = 1, km
          if ( j .lt. kb ) then
            fc(j) = s0 * fc(j) * f3 / f2
          else
            fc(j) = s0 * fc(j)
          end if
        end do

      else if ( kd .eq. 4 ) then

        do k = km, 3, -1
          v = u
          u = f
          f = ( a - 4.0D+00 * k * k ) * u / q - v
          if ( abs ( fc(k) ) .le. abs ( f ) ) then
            fc(k-1) = f
            s = s + f * f
          else
            kb = k
            f3 = fc(k)
            go to 20
          end if
        end do

        fc(1) = q / ( a - 4.0D+00 ) * fc(2)
        s = s + fc(1) * fc(1)
        s0 = sqrt ( 1.0D+00 / s )
        do k = 1, km
          fc(k) = s0 * fc(k)
        end do

        go to 30

20      continue

        fc(1) = 1.0D-100
        fc(2) = ( a - 4.0D+00 ) / q * fc(1)
        sp = 0.0D+00
        u = fc(1)
        f1 = fc(2)

        do i = 2, kb - 1
          v = u
          u = f1
          f1 = ( a - 4.0D+00 * i * i ) * u / q - v
          if ( i .ne. kb - 1 ) then
            fc(i+1) = f1
            sp = sp + f1 * f1
          else
            f2 = f1
          end if
        end do

        sp = sp + fc(1) ** 2 + fc(2) ** 2
        ss = s + sp * ( f3 / f2 ) ** 2
        s0 = 1.0D+00 / sqrt ( ss )

        do j = 1, km
          if ( j .lt. kb ) then
            fc(j) = s0 * fc(j) * f3 / f2
          else
            fc(j) = s0 * fc(j)
          end if
        end do

      end if

30    continue

      if ( fc(1) .lt. 0.0D+00 ) then
        do j = 1, km
          fc(j) = -fc(j)
        end do
      end if

      return
      end
      subroutine fcs ( x, c, s )

c*********************************************************************72
c
cc FCS computes Fresnel integrals C(x) and S(x).
c
c  Licensing:
c
c    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
c    they give permission to incorporate this routine into a user program 
c    provided that the copyright is acknowledged.
c
c  Modified:
c
c    17 July 2012
c
c  Author:
c
c    Shanjie Zhang, Jianming Jin
c
c  Reference:
c
c    Shanjie Zhang, Jianming Jin,
c    Computation of Special Functions,
c    Wiley, 1996,
c    ISBN: 0-471-11963-6,
c    LC: QA351.C45.
c 
c  Parameters:
c
c    Input, double precision X, the argument.
c
c    Output, double precision C, S, the function values.
c
      implicit none

      double precision c
      double precision eps
      double precision f
      double precision f0
      double precision f1
      double precision g
      integer k
      integer m
      double precision pi
      double precision px
      double precision q
      double precision r
      double precision s
      double precision su
      double precision t
      double precision t0
      double precision t2
      double precision x
      double precision xa

      eps = 1.0D-15
      pi = 3.141592653589793D+00
      xa = abs ( x )
      px = pi * xa
      t = 0.5D+00 * px * xa
      t2 = t * t

      if ( xa .eq. 0.0D+00 ) then

        c = 0.0D+00
        s = 0.0D+00

      else if ( xa .lt. 2.5D+00 ) then

        r = xa
        c = r
        do k = 1, 50
          r = -0.5D+00 * r * ( 4.0D+00 * k - 3.0D+00 ) / k 
     &      / ( 2.0D+00 * k - 1.0D+00 ) / ( 4.0D+00 * k + 1.0D+00 ) * t2
          c = c + r
          if ( abs ( r ) .lt. abs ( c ) * eps ) then
            go to 10
          end if
        end do

10      continue

        s = xa * t / 3.0D+00
        r = s
        do k = 1, 50
          r = - 0.5D+00 * r * ( 4.0D+00 * k - 1.0D+00 ) / k 
     &      / ( 2.0D+00 * k + 1.0D+00 ) / ( 4.0D+00 * k + 3.0D+00 ) * t2
          s = s + r
          if ( abs ( r ) .lt. abs ( s ) * eps ) then
            if ( x .lt. 0.0D+00 ) then
              c = -c
              s = -s
            end if
            return
          end if
        end do

      else if ( xa .lt. 4.5D+00 ) then

        m = int ( 42.0D+00 + 1.75D+00 * t )
        su = 0.0D+00
        c = 0.0D+00
        s = 0.0D+00
        f1 = 0.0D+00
        f0 = 1.0D-100

        do k = m, 0, -1
          f = ( 2.0D+00 * k + 3.0D+00 ) * f0 / t - f1
          if ( k .eq. int ( k / 2 ) * 2 ) then
            c = c + f
          else
            s = s + f
          end if
          su = su + ( 2.0D+00 * k + 1.0D+00 ) * f * f
          f1 = f0
          f0 = f
        end do

        q = sqrt ( su )
        c = c * xa / q
        s = s * xa / q

      else

        r = 1.0D+00
        f = 1.0D+00
        do k = 1, 20
          r = -0.25D+00 * r * ( 4.0D+00 * k - 1.0D+00 ) 
     &      * ( 4.0D+00 * k - 3.0D+00 ) / t2
          f = f + r
        end do
        r = 1.0D+00 / ( px * xa )
        g = r
        do k = 1, 12
          r = -0.25D+00 * r * ( 4.0D+00 * k + 1.0D+00 ) 
     &      * ( 4.0D+00 * k - 1.0D+00 ) / t2
          g = g + r
        end do

        t0 = t - int ( t / ( 2.0D+00 * pi ) ) * 2.0D+00 * pi
        c = 0.5D+00 + ( f * sin ( t0 ) - g * cos ( t0 ) ) / px
        s = 0.5D+00 - ( f * cos ( t0 ) + g * sin ( t0 ) ) / px

      end if

      if ( x .lt. 0.0D+00 ) then
        c = -c
        s = -s
      end if

      return
      end
      subroutine fcszo ( kf, nt, zo )

c*********************************************************************72
c
cc FCSZO computes complex zeros of Fresnel integrals C(x) or S(x).
c
c  Licensing:
c
c    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
c    they give permission to incorporate this routine into a user program 
c    provided that the copyright is acknowledged.
c
c  Modified:
c
c    17 July 2012
c
c  Author:
c
c    Shanjie Zhang, Jianming Jin
c
c  Reference:
c
c    Shanjie Zhang, Jianming Jin,
c    Computation of Special Functions,
c    Wiley, 1996,
c    ISBN: 0-471-11963-6,
c    LC: QA351.C45.
c
c  Parameters:
c
c    Input, integer KF, the function code.
c    1 for C(z);
c    2 for S(z)
c
c    Input, integer NT, the total number of zeros desired.
c
c    Output, complex*16 Z0(NT), the zeros.
c
      implicit none

      integer nt

      integer i
      integer it
      integer j
      integer kf
      integer nr
      double precision pi
      double precision psq
      double precision px
      double precision py
      double precision w
      double precision w0
      complex*16 z
      complex*16 zd
      complex*16 zf
      complex*16 zfd
      complex*16 zgd
      complex*16 zo(nt)
      complex*16 zp
      complex*16 zq
      complex*16 zw

      pi = 3.141592653589793D+00

      do nr = 1, nt

        if ( kf .eq. 1 ) then
          psq = sqrt ( 4.0D+00 * nr - 1.0D+00 )
        else
          psq = 2.0D+00 * sqrt ( dble ( nr ) )
        end if

        px = psq - log ( pi * psq ) / ( pi * pi * psq ** 3.0D+00 )
        py = log ( pi * psq ) / ( pi * psq )
        z = cmplx ( px, py )

        if ( kf .eq. 2 ) then
          if ( nr .eq. 2 ) then
            z = cmplx ( 2.8334D+00, 0.2443D+00 )
          else if ( nr .eq. 3 ) then
            z = cmplx ( 3.4674D+00, 0.2185D+00 )
          else if ( nr .eq. 4 ) then
            z = cmplx ( 4.0025D+00, 0.2008D+00 )
          end if
        end if

        it = 0

10      continue

        it = it + 1

        if ( kf .eq. 1 ) then
          call cfc ( z, zf, zd )
        else
          call cfs ( z, zf, zd )
        end if

        zp = cmplx ( 1.0D+00, 0.0D+00 )
        do i = 1, nr - 1
          zp = zp * ( z - zo(i) )
        end do
        zfd = zf / zp
        zq = cmplx ( 0.0D+00, 0.0D+00 )

        do i = 1, nr - 1
          zw = cmplx ( 1.0D+00, 0.0D+00 )
          do j = 1, nr - 1
            if ( j .ne. i ) then
              zw = zw * ( z - zo(j) )
            end if
          end do
          zq = zq + zw
        end do

        zgd = ( zd - zq * zfd ) / zp
        z = z - zfd / zgd
        w0 = w
        w = cdabs ( z )

        if ( it .le. 50 .and. abs ( ( w - w0 ) / w ) .gt. 1.0D-12 ) then
          go to 10
        end if

        zo(nr) = z

      end do

      return
      end
      subroutine ffk ( ks, x, fr, fi, fm, fa, gr, gi, gm, ga )

c*********************************************************************72
c
cc FFK computes modified Fresnel integrals F+/-(x) and K+/-(x).
c
c  Licensing:
c
c    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
c    they give permission to incorporate this routine into a user program 
c    provided that the copyright is acknowledged.
c
c  Modified:
c
c    23 July 2012
c
c  Author:
c
c    Shanjie Zhang, Jianming Jin
c
c  Reference:
c
c    Shanjie Zhang, Jianming Jin,
c    Computation of Special Functions,
c    Wiley, 1996,
c    ISBN: 0-471-11963-6,
c    LC: QA351.C45.
c
c  Parameters:
c
c    Input, integer KS, the sign code.
c    0, to calculate F+(x) and K+(x);
c    1, to calculate F_(x) and K_(x).
c
c    Input, double precision X, the argument.
c
c    Output, double precision FR, FI, FM, FA, the values of
c    Re[F+/-(x)], Im[F+/-(x)], |F+/-(x)|, Arg[F+/-(x)]  (Degs.).
c
c    Output, double precision GR, GI, GM, GA, the values of
c    Re[K+/-(x)], Im[K+/-(x)], |K+/-(x)|, Arg[K+/-(x)]  (Degs.).
c       
      implicit none

      double precision c1
      double precision cs
      double precision eps
      double precision fa
      double precision fi
      double precision fi0
      double precision fm
      double precision fr
      double precision ga
      double precision gi
      double precision gm
      double precision gr
      integer k
      integer ks
      integer m
      double precision p2p
      double precision pi
      double precision pp2
      double precision s1
      double precision srd
      double precision ss
      double precision x
      double precision x2
      double precision x4
      double precision xa
      double precision xc
      double precision xf
      double precision xf0
      double precision xf1
      double precision xg
      double precision xp
      double precision xq
      double precision xq2
      double precision xr
      double precision xs
      double precision xsu
      double precision xw

      srd = 57.29577951308233D+00
      eps = 1.0D-15
      pi = 3.141592653589793D+00
      pp2 = 1.2533141373155D+00
      p2p = 0.7978845608028654D+00
      xa = abs ( x )
      x2 = x * x
      x4 = x2 * x2

      if ( x .eq. 0.0D+00 ) then

        fr = 0.5D+00 * sqrt ( 0.5D+00 * pi )
        fi = ( -1.0D+00 ) ** ks * fr
        fm = sqrt ( 0.25D+00 * pi )
        fa = ( -1.0D+00 ) ** ks * 45.0D+00
        gr = 0.5D+00
        gi = 0.0D+00
        gm = 0.5D+00
        ga = 0.0D+00

      else

        if ( xa .le. 2.5D+00 ) then

          xr = p2p * xa
          c1 = xr
          do k = 1, 50
            xr = -0.5D+00 * xr * ( 4.0D+00 * k - 3.0D+00 ) / k 
     &        / ( 2.0D+00 * k - 1.0D+00 )
     &        / ( 4.0D+00 * k + 1.0D+00 ) * x4
            c1 = c1 + xr
            if ( abs ( xr / c1 ) .lt. eps ) then
              go to 10
            end if
          end do

10        continue

          s1 = p2p * xa * xa * xa / 3.0D+00
          xr = s1
          do k = 1, 50
            xr = -0.5D+00 * xr * ( 4.0D+00 * k - 1.0D+00 )
     &        / k / ( 2.0D+00 * k + 1.0D+00 )
     &        / ( 4.0D+00 * k + 3.0D+00 ) * x4
            s1 = s1 + xr
            if ( abs ( xr / s1 ) .lt. eps ) then
              go to 20
            end if
          end do

20        continue

        else if ( xa .lt. 5.5D+00 ) then

          m = int ( 42.0D+00 + 1.75D+00 * x2 )
          xsu = 0.0D+00
          xc = 0.0D+00
          xs = 0.0D+00
          xf1 = 0.0D+00
          xf0 = 1.0D-100
          do k = m, 0, -1
            xf = ( 2.0D+00 * k + 3.0D+00 ) * xf0 / x2 - xf1
            if ( k .eq. 2 * int ( k / 2 ) )  then
              xc = xc + xf
            else
              xs = xs + xf
            end if
            xsu = xsu + ( 2.0D+00 * k + 1.0D+00 ) * xf * xf
            xf1 = xf0
            xf0 = xf
          end do
          xq = sqrt ( xsu )
          xw = p2p * xa / xq
          c1 = xc * xw
          s1 = xs * xw

        else

          xr = 1.0D+00
          xf = 1.0D+00
          do k = 1, 12
            xr = -0.25D+00 * xr * ( 4.0D+00 * k - 1.0D+00 )
     &        * ( 4.0D+00 * k - 3.0D+00 ) / x4
            xf = xf + xr
          end do
          xr = 1.0D+00 / ( 2.0D+00 * xa * xa )
          xg = xr
          do k = 1, 12
            xr = -0.25D+00 * xr * ( 4.0D+00 * k + 1.0D+00 )
     &        * ( 4.0D+00 * k - 1.0D+00 ) / x4
            xg = xg + xr
          end do
          c1 = 0.5D+00 + ( xf * sin ( x2 ) - xg * cos ( x2 ) )
     &      / sqrt ( 2.0D+00 * pi ) / xa
          s1 = 0.5D+00 - ( xf * cos ( x2 ) + xg * sin ( x2 ) )
     &      / sqrt ( 2.0D+00 * pi ) / xa

        end if
 
        fr = pp2 * ( 0.5D+00 - c1 )
        fi0 = pp2 * ( 0.5D+00 - s1 )
        fi = ( -1.0D+00 ) ** ks * fi0
        fm = sqrt ( fr * fr + fi * fi )

        if ( 0.0D+00 .le. fr ) then
          fa = srd * atan ( fi / fr )
        else if ( 0.0D+00 .lt. fi ) then
          fa = srd * ( atan ( fi / fr ) + pi )
        else if ( fi .lt. 0.0D+00 ) then
          fa = srd * ( atan ( fi / fr ) - pi )
        end if

        xp = x * x + pi / 4.0D+00
        cs = cos ( xp )
        ss = sin ( xp )
        xq2 = 1.0D+00 / sqrt ( pi )
        gr = xq2 * ( fr * cs + fi0 * ss )
        gi = ( -1.0D+00 ) ** ks * xq2 * ( fi0 * cs - fr * ss )
        gm = sqrt ( gr * gr + gi * gi )

        if ( 0.0D+00 .le. gr ) then
          ga = srd * atan ( gi / gr )
        else if ( 0.0D+00 .lt. gi ) then
          ga = srd * ( atan ( gi / gr ) + pi )
        else if ( gi .lt. 0.0D+00 ) then
          ga = srd * ( atan ( gi / gr ) - pi )
        end if

        if ( x .lt. 0.0D+00 ) then
          fr = pp2 - fr
          fi = ( -1.0D+00 ) ** ks * pp2 - fi
          fm = sqrt ( fr * fr + fi * fi )
          fa = srd * atan ( fi / fr )
          gr = cos ( x * x ) - gr
          gi = - ( -1.0D+00 ) ** ks * sin ( x * x ) - gi
          gm = sqrt ( gr * gr + gi * gi )
          ga = srd * atan ( gi / gr )
        end if

      end if

      return
      end
      subroutine gaih ( x, ga )

c*********************************************************************72
c
cc GAIH computes the GammaH function.
c
c  Licensing:
c
c    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
c    they give permission to incorporate this routine into a user program 
c    provided that the copyright is acknowledged.
c
c  Modified:
c
c    09 July 2012
c
c  Author:
c
c    Shanjie Zhang, Jianming Jin
c
c  Reference:
c
c    Shanjie Zhang, Jianming Jin,
c    Computation of Special Functions,
c    Wiley, 1996,
c    ISBN: 0-471-11963-6,
c    LC: QA351.C45.
c
c  Parameters:
c
c    Input, double precision X, the argument.
c
c    Output, double precision GA, the function value.
c
      implicit none

      double precision ga
      integer k
      integer m
      integer m1
      double precision pi
      double precision x

      pi = 3.141592653589793D+00

      if ( x .eq. int ( x ) .and. 0.0 .lt. x ) then
        ga = 1.0D+00
        m1 = int ( x - 1.0D+00 )
        do k = 2, m1
          ga = ga * k
        end do
      else if ( x + 0.5D+00 .eq. int ( x + 0.5D+00) .and. 
     &  0.0D+00 .lt. x ) then
        m = int ( x )
        ga = sqrt ( pi )
        do k = 1, m
          ga = 0.5D+00 * ga * ( 2.0D+00 * k - 1.0D+00 )
        end do
      end if

      return
      end
      subroutine gam0 ( x, ga )

c*********************************************************************72
c
cc GAM0 computes the Gamma function for the LAMV function.
c
c  Licensing:
c
c    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
c    they give permission to incorporate this routine into a user program 
c    provided that the copyright is acknowledged.
c
c  Modified:
c
c    09 July 2012
c
c  Author:
c
c    Shanjie Zhang, Jianming Jin
c
c  Reference:
c
c    Shanjie Zhang, Jianming Jin,
c    Computation of Special Functions,
c    Wiley, 1996,
c    ISBN: 0-471-11963-6,
c    LC: QA351.C45.
c
c  Parameters:
c
c    Input, double precision X, the argument.
c
c    Output, double precision GA, the function value.
c   
      implicit none

      double precision g(25)
      double precision ga
      double precision gr
      integer k
      double precision x

      data g /
     &   1.0D+00,
     &   0.5772156649015329D+00,
     &  -0.6558780715202538D+00,
     &  -0.420026350340952D-01,
     &   0.1665386113822915D+00,
     &  -0.421977345555443D-01,
     &  -0.96219715278770D-02,
     &   0.72189432466630D-02,
     &  -0.11651675918591D-02,
     &  -0.2152416741149D-03,
     &   0.1280502823882D-03,
     &  -0.201348547807D-04,
     &  -0.12504934821D-05,
     &   0.11330272320D-05,
     &  -0.2056338417D-06,
     &   0.61160950D-08,
     &   0.50020075D-08,
     &  -0.11812746D-08,
     &   0.1043427D-09,
     &   0.77823D-11,
     &  -0.36968D-11,
     &   0.51D-12,
     &  -0.206D-13,
     &  -0.54D-14,
     &   0.14D-14/

      gr = g(25)
      do k = 24, 1, -1
        gr = gr * x + g(k)
      end do

      ga = 1.0D+00 / ( gr * x )

      return
      end
      subroutine gamma ( x, ga )

c*********************************************************************72
c
cc GAMMA computes the Gamma function.
c
c  Licensing:
c
c    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
c    they give permission to incorporate this routine into a user program 
c    provided that the copyright is acknowledged.
c
c  Modified:
c
c    14 July 2012
c
c  Author:
c
c    Shanjie Zhang, Jianming Jin
c
c  Reference:
c
c    Shanjie Zhang, Jianming Jin,
c    Computation of Special Functions,
c    Wiley, 1996,
c    ISBN: 0-471-11963-6,
c    LC: QA351.C45.
c
c  Parameters:
c
c    Input, double precision X, the argument.
c
c    Output, double precision GA, the function value.
c
      implicit none

      double precision g(26)
      double precision ga
      double precision gr
      integer k
      integer m
      integer m1
      double precision pi
      double precision r
      double precision x
      double precision z

      save g

      data g /
     &   1.0D+00, 
     &   0.5772156649015329D+00,
     &  -0.6558780715202538D+00, 
     &  -0.420026350340952D-01,
     &   0.1665386113822915D+00,
     &  -0.421977345555443D-01,
     &  -0.96219715278770D-02, 
     &   0.72189432466630D-02,
     &  -0.11651675918591D-02, 
     &  -0.2152416741149D-03,
     &   0.1280502823882D-03,
     &  -0.201348547807D-04,
     &  -0.12504934821D-05, 
     &   0.11330272320D-05,
     &  -0.2056338417D-06, 
     &   0.61160950D-08,
     &   0.50020075D-08, 
     &  -0.11812746D-08,
     &   0.1043427D-09, 
     &   0.77823D-11,
     &  -0.36968D-11, 
     &   0.51D-12,
     &  -0.206D-13, 
     &  -0.54D-14, 
     &   0.14D-14, 
     &   0.1D-15 /

      pi = 3.141592653589793D+00

      if ( x .eq. int ( x ) ) then

        if ( 0.0D+00 .lt. x ) then
          ga = 1.0D+00
          m1 = x - 1
          do k = 2, m1
            ga = ga * k
          end do
        else
          ga = 1.0D+300
        end if

      else

        if ( 1.0D+00 .lt. abs ( x ) ) then
          z = abs ( x )
          m = int ( z )
          r = 1.0D+00
          do k = 1, m
            r = r * ( z - k )
          end do
          z = z - m
        else
          z = x
        end if

        gr = g(26)
        do k = 25, 1, -1
          gr = gr * z + g(k)
        end do
        ga = 1.0D+00 / ( gr * z )

        if ( 1.0D+00 .lt. abs ( x ) ) then
          ga = ga * r
          if ( x .lt. 0.0D+00 ) then
            ga = - pi / ( x * ga * sin ( pi * x ) )
          end if
        end if

      end if

      return
      end
      subroutine gmn ( m, n, c, x, bk, gf, gd )

c*********************************************************************72
c
cc GMN computes quantities needed for oblate radial functions with small argument.
c
c  Discussion:
c
c    This procedure computes Gmn(-ic,ix) and its derivative for oblate
c    radial functions with a small argument.
c
c  Licensing:
c
c    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
c    they give permission to incorporate this routine into a user program 
c    provided that the copyright is acknowledged.
c
c  Modified:
c
c    15 July 2012
c
c  Author:
c
c    Shanjie Zhang, Jianming Jin
c
c  Reference:
c
c    Shanjie Zhang, Jianming Jin,
c    Computation of Special Functions,
c    Wiley, 1996,
c    ISBN: 0-471-11963-6,
c    LC: QA351.C45.
c
c  Parameters:
c
c    Input, integer M, the mode parameter;  M = 0, 1, 2, ...
c
c    Input, integer N, mode parameter, N = M, M + 1, M + 2, ...
c
c    Input, double precision C, spheroidal parameter.
c
c    Input, double precision X, the argument.
c
c    Input, double precision BK(*), coefficients.
c
c    Output, double precision GF, GD, the value of Gmn(-C,X) and Gmn'(-C,X).
c
      implicit none

      double precision bk(200)
      double precision c
      double precision eps
      double precision gd
      double precision gd0
      double precision gd1
      double precision gf
      double precision gf0
      double precision gw
      integer ip
      integer k
      integer m
      integer n
      integer nm
      double precision x
      double precision xm

      eps = 1.0D-14

      if ( n - m .eq. 2 * int ( ( n - m ) / 2 ) ) then
        ip = 0
      else
        ip = 1
      end if

      nm = 25 + int ( 0.5D+00 * ( n - m ) + c )
      xm = ( 1.0D+00 + x * x ) ** ( -0.5D+00 * m )
      gf0 = 0.0D+00
      do k = 1, nm
        gf0 = gf0 + bk(k) * x ** ( 2.0D+00 * k - 2.0D+00 )
        if ( abs ( ( gf0 - gw ) / gf0 ) .lt. eps .and. 10 .le. k ) then
          go to 10
        end if
        gw = gf0
      end do

10    continue

      gf = xm * gf0 * x ** ( 1 - ip )

      gd1 = - m * x / ( 1.0D+00 + x * x ) * gf
      gd0 = 0.0D+00

      do k = 1, nm

        if ( ip .eq. 0 ) then
          gd0 = gd0 + ( 2.0D+00 * k - 1.0D+00 ) * bk(k) 
     &      * x ** ( 2.0D+00 * k - 2.0D+00 )
        else
          gd0 = gd0 + 2.0D+00 * k * bk(k+1) 
     &      * x ** ( 2.0D+00 * k - 1.0D+00 )
        end if

        if ( abs ( ( gd0 - gw ) / gd0 ) .lt. eps .and. 10 .le. k ) then
          go to 20
        end if

        gw = gd0

      end do

20    continue

      gd = gd1 + xm * gd0

      return
      end
      subroutine herzo ( n, x, w )

c*********************************************************************72
c
cc HERZO computes the zeros the Hermite polynomial Hn(x).
c
c  Discussion:
c
c    This procedure computes the zeros of Hermite polynomial Ln(x)
c    in the interval [-1,+1], and the corresponding
c    weighting coefficients for Gauss-Hermite integration.
c
c  Licensing:
c
c    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
c    they give permission to incorporate this routine into a user program 
c    provided that the copyright is acknowledged.
c
c  Modified:
c
c    15 July 2012
c
c  Author:
c
c    Shanjie Zhang, Jianming Jin
c
c  Reference:
c
c    Shanjie Zhang, Jianming Jin,
c    Computation of Special Functions,
c    Wiley, 1996,
c    ISBN: 0-471-11963-6,
c    LC: QA351.C45.
c
c  Parameters:
c
c    Input, integer N, the order of the polynomial.
c
c    Output, double precision X(N), the zeros.
c
c    Output, double precision W(N), the corresponding weights.
c
      implicit none

      integer n

      double precision f0
      double precision f1
      double precision fd
      double precision gd
      double precision hd
      double precision hf
      double precision hn
      integer i
      integer it
      integer j
      integer k
      integer nr
      double precision p
      double precision q
      double precision r
      double precision r1
      double precision r2
      double precision w(n)
      double precision wp
      double precision x(n)
      double precision x0
      double precision z
      double precision z0
      double precision zl

      hn = 1.0D+00 / n
      zl = -1.1611D+00 + 1.46D+00 * sqrt ( dble ( n ) )

      do nr = 1, n / 2

        if ( nr .eq. 1 ) then
          z = zl
        else
          z = z - hn * ( n / 2 + 1 - nr )
        end if

        it = 0

10      continue

        it = it + 1
        z0 = z
        f0 = 1.0D+00
        f1 = 2.0D+00 * z
        do k = 2, n
          hf = 2.0D+00 * z * f1 - 2.0D+00 * ( k - 1.0D+00 ) * f0
          hd = 2.0D+00 * k * f1
          f0 = f1
          f1 = hf
        end do

        p = 1.0D+00
        do i = 1, nr - 1
          p = p * ( z - x(i) )
        end do
        fd = hf / p

        q = 0.0D+00
        do i = 1, nr - 1
          wp = 1.0D+00
          do j = 1, nr - 1
            if ( j .ne. i ) then
              wp = wp * ( z - x(j) )
            end if
          end do
          q = q + wp
        end do

        gd = ( hd - q * fd ) / p
        z = z - fd / gd

        if ( it .le. 40 .and. 1.0D-15 .lt. abs ( ( z - z0 ) / z ) ) then
          go to 10
        end if

        x(nr) = z
        x(n+1-nr) = -z
        r = 1.0D+00
        do k = 1, n
          r = 2.0D+00 * r * k
        end do
        w(nr) = 3.544907701811D+00 * r / ( hd * hd )
        w(n+1-nr) = w(nr)

      end do

      if ( n .ne. 2 * int ( n / 2 ) ) then
        r1 = 1.0D+00
        r2 = 1.0D+00
        do j = 1, n
          r1 = 2.0D+00 * r1 * j
          if ( ( n + 1 ) / 2 .le. j ) then
            r2 = r2 * j
          end if
        end do
        w(n/2+1) = 0.88622692545276D+00 * r1 / ( r2 * r2 )
        x(n/2+1) = 0.0D+00
      end if

      return
      end
      subroutine hygfx ( a, b, c, x, hf )

c*********************************************************************72
c
cc HYGFX computes the hypergeometric function F(a,b,c,x).
c
c  Licensing:
c
c    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
c    they give permission to incorporate this routine into a user program 
c    provided that the copyright is acknowledged.
c
c  Modified:
c
c    04 April 2012
c
c  Author:
c
c    Shanjie Zhang, Jianming Jin
c
c  Reference:
c
c    Shanjie Zhang, Jianming Jin,
c    Computation of Special Functions,
c    Wiley, 1996,
c    ISBN: 0-471-11963-6,
c    LC: QA351.C45.
c 
c  Parameters:
c
c    Input, double precision A, B, C, parameters.
c
c    Input, double precision X, the argument.
c
c    Output, double precision HF, the function value.
c
      implicit double precision (a-h,o-z)

      double precision a
      double precision b
      double precision c
      double precision c0
      double precision el
      double precision eps
      double precision g0
      double precision g1
      double precision g2
      double precision g3
      double precision gc
      double precision hf
      integer k
      logical l0
      logical l1
      logical l2
      logical l3
      logical l4
      logical l5
      double precision pi
      double precision r
      double precision x

      pi = 3.141592653589793D+00
      el = 0.5772156649015329D+00
      l0 = c .eq. int ( c ) .and. c .lt. 0.0D+00
      l1 = 1.0D+00 - x .lt. 1.0D-15 .and. c - a - b .le. 0.0D+00
      l2 = a .eq. int ( a ) .and. a .lt. 0.0D+00
      l3 = b .eq. int ( b ) .and. b .lt. 0.0D+00
      l4 = c - a .eq. int ( c - a ) .and. c - a .le. 0.0D+00
      l5 = c - b .eq. int ( c - b ) .and. c - b .le. 0.0D+00

      if ( l0 .or. l1 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'HYGFX - Fatal error!'
        write ( *, '(a)' ) '  The hypergeometric series is divergent.'
        stop
      end if

      if ( 0.95D+00 .lt. x ) then
        eps = 1.0D-08
      else
        eps = 1.0D-15
      end if

      if ( x .eq. 0.0D+00 .or.a .eq. 0.0D+00 .or.b .eq. 0.0D+00 ) then
        hf = 1.0D+00
        return
      else if ( 1.0D+00 - x .eq. eps .and. 0.0D+00 .lt. c - a - b ) then
        call gamma ( c, gc )
        call gamma ( c - a - b, gcab )
        call gamma ( c - a, gca )
        call gamma ( c - b, gcb )
        hf = gc * gcab / ( gca * gcb )
        return
      else if ( 1.0D+00 + x .le. eps .and. 
     &  abs ( c - a + b - 1.0D+00 ) .le. eps ) then
        g0 = dsqrt ( pi ) * 2.0D+00 ** ( - a )
        call gamma ( c, g1 )
        call gamma ( 1.0D+00 + a / 2.0D+00 - b, g2 )
        call gamma ( 0.5D+00 + 0.5D+00 * a, g3 )
        hf = g0 * g1 / ( g2 * g3 )
        return
      else if ( l2 .or. l3 ) then
        if ( l2 ) then
          nm = int ( abs ( a ) )
        end if
        if ( l3 ) then
          nm = int ( abs ( b ) )
        end if
        hf = 1.0D+00
        r = 1.0D+00
        do k = 1, nm
          r = r * ( a + k - 1.0D+00 ) * ( b + k - 1.0D+00) 
     &      / ( k * ( c + k - 1.0D+00 ) ) * x
          hf = hf + r
        end do
        return
      else if ( l4 .or. l5 ) then
        if ( l4 ) then
          nm = int ( abs ( c - a ) )
        end if
        if ( l5 ) then
          nm = int ( abs ( c - b ) )
        end if
        hf = 1.0D+00
        r = 1.0D+00
        do k = 1, nm
          r = r * ( c - a + k - 1.0D+00 ) * ( c - b + k - 1.0D+00 )
     &      / ( k * ( c + k - 1.0D+00 ) ) * x
          hf = hf + r
        end do
        hf = ( 1.0D+00 - x ) ** ( c - a - b ) * hf
        return
      end if

      aa = a
      bb = b
      x1 = x

      if ( x .lt. 0.0D+00 ) then
        x = x / ( x - 1.0D+00 )
        if ( a .lt. c .and. b .lt. a .and. 0.0D+00 .lt. b ) then
          a = bb
          b = aa
        end if
        b = c - b
      end if

      if ( 0.75D+00 .le. x ) then
        gm = 0.0D+00
        if ( abs ( c - a - b - int ( c - a - b ) ) .lt. 1.0D-15 ) then
          m = int ( c - a - b )
          call gamma ( a, ga )
          call gamma ( b, gb )
          call gamma ( c, gc )
          call gamma ( a + m, gam )
          call gamma ( b + m, gbm )
          call psi ( a, pa )
          call psi ( b, pb )
          if ( m .ne. 0 ) then
            gm = 1.0D+00
          end if
          do j = 1, abs ( m ) - 1
            gm = gm * j
          end do
          rm = 1.0D+00
          do j = 1, abs ( m )
            rm = rm * j
          end do
          f0 = 1.0D+00
          r0 = 1.0D+00
          r1 = 1.0D+00
          sp0 = 0.0D+00
          sp = 0.0D+00
          if ( 0 .le. m ) then
            c0 = gm * gc / ( gam * gbm )
            c1 = - gc * ( x - 1.0D+00 ) ** m / ( ga * gb * rm )
            do k = 1, m - 1
              r0 = r0 * ( a + k - 1.0D+00 ) * ( b + k - 1.0D+00 )
     &          / ( k * ( k - m ) ) * ( 1.0D+00 - x )
              f0 = f0 + r0
            end do

            do k = 1, m
              sp0 = sp0 + 1.0D+00 / ( a + k - 1.0D+00) 
     &          + 1.0D+00 / ( b + k - 1.0D+00 ) - 1.0D+00 / k
            end do

            f1 = pa + pb + sp0 + 2.0D+00 * el + log ( 1.0D+00 - x )

            do k = 1, 250
              sp = sp + (1.0D+00-a) / (k*(a+k-1.0))
     &          +(1.0-b)/(k*(b+k-1.0))
              sm = 0.0D+00
              do j = 1,m
                sm = sm+(1.0D+00-a)/((j+k)*(a+j+k-1.0))+1.0/
     &            (b+j+k-1.0D+00 )
              end do
              rp = pa+pb+2.0D+00*el+sp+sm+log(1.0D+00-x)
              r1 = r1*(a+m+k-1.0D+00)*(b+m+k-1.0)
     &          /(k*(m+k))*(1.0D+00 -x)
              f1 = f1+r1*rp
              if (abs(f1-hw) .lt. abs(f1)*eps) then
                go to 10
              end if
              hw = f1
            end do

10          continue

            hf = f0*c0+f1*c1

          else if (m .lt. 0) then

            m = -m
            c0 = gm*gc/(ga*gb*(1.0D+00-x)**m)
            c1 = -(-1)**m*gc/(gam*gbm*rm)

                 do k = 1,m-1
                    r0 = r0*(a-m+k-1.0D+00)*(b-m+k-1.0D+00)
     &                /(k*(k-m))*(1.0D+00-x)
                    f0 = f0+r0
                 end do

                 do k = 1,m
                   sp0 = sp0+1.0D+00/k
                 end do
                 f1 = pa+pb-sp0+2.0D+00*el+log(1.0D+00-x)
                 do k = 1,250
                    sp = sp + (1.0D+00-a)/(k*(a+k-1.0D+00))
     &                +(1.0D+00-b)/(k*(b+k-1.0D+00))
                    sm = 0.0D+00
                    do j = 1,m
                      sm = sm+1.0D+00/(j+k)
                    end do
                    rp = pa+pb+2.0D+00*el+sp-sm+log(1.0D+00-x)
                    r1 = r1*(a+k-1.0D+00)*(b+k-1.0D+00)
     &                /(k*(m+k))*(1.0D+00-x)
                    f1 = f1+r1*rp
                    if (abs(f1-hw) .lt. abs(f1)*eps) then
                      go to 20
                    end if
                    hw = f1
                 end do

20               continue

                 hf = f0*c0+f1*c1
              end if
           else
              call gamma(a,ga)
              call gamma(b,gb)
              call gamma(c,gc)
              call gamma(c-a,gca)
              call gamma(c-b,gcb)
              call gamma(c-a-b,gcab)
              call gamma(a+b-c,gabc)
              c0 = gc*gcab/(gca*gcb)
              c1 = gc*gabc/(ga*gb)*(1.0D+00-x)**(c-a-b)
              hf = 0.0D+00
              r0 = c0
              r1 = c1
              do k = 1,250
                 r0 = r0*(a+k-1.0D+00)*(b+k-1.0D+00)
     &             /(k*(a+b-c+k))*(1.0D+00-x)
                 r1 = r1*(c-a+k-1.0D+00)*(c-b+k-1.0D+00)/(k*(c-a-b+k))
     &              *(1.0D+00-x)
                 hf = hf+r0+r1
                 if (abs(hf-hw) .lt. abs(hf)*eps) then
                   go to 30
                 end if
                 hw = hf
              end do

30            continue

              hf = hf+c0+c1

           end if
        else
           a0 = 1.0D+00
           if (c.gt.a.and.c .lt. 2.0D+00*a.and.
     &         c.gt.b.and.c .lt. 2.0D+00*b) then
              a0 = ( 1.0D+00 - x ) ** ( c - a - b )
              a = c-a
              b = c-b
           end if
           hf = 1.0D+00
           r = 1.0D+00
           do k = 1,250
              r = r*(a+k-1.0D+00)*(b+k-1.0D+00)/(k*(c+k-1.0D+00))*x
              hf = hf+r
              if (abs(hf-hw).le.abs(hf)*eps) then
                go to 40
              end if
              hw = hf
           end do

40         continue

           hf = a0*hf

      end if

      if ( x1 .lt. 0.0D+00 ) then
        x = x1
        c0 = 1.0D+00 / ( 1.0D+00 - x ) ** aa
        hf = c0 * hf
      end if

      a = aa
      b = bb

      if ( 120 .lt. k ) then
        write(*,115)
115     format('Warning! You should check the accuracy')
      end if

      return
      end
      subroutine hygfz ( a, b, c, z, zhf )

c*********************************************************************72
c
cc HYGFZ computes the hypergeometric function F(a,b,c,x) for complex argument.
c
c  Licensing:
c
c    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
c    they give permission to incorporate this routine into a user program 
c    provided that the copyright is acknowledged.
c
c  Modified:
c
c    03 August 2012
c
c  Author:
c
c    Shanjie Zhang, Jianming Jin
c
c  Reference:
c
c    Shanjie Zhang, Jianming Jin,
c    Computation of Special Functions,
c    Wiley, 1996,
c    ISBN: 0-471-11963-6,
c    LC: QA351.C45.
c
c  Parameters:
c
c    Input, double precision A, B, C, parameters.
c
c    Input, complex*16 Z, the argument.
c
c    Output, complex*16 ZHF, the value of F(a,b,c,z).
c
      implicit none

      double precision a
      double precision a0
      double precision aa
      double precision b
      double precision bb
      double precision c
      double precision ca
      double precision cb
      double precision el
      double precision eps
      double precision g0
      double precision g1
      double precision g2
      double precision g3
      double precision ga
      double precision gab
      double precision gabc
      double precision gam
      double precision gb
      double precision gba
      double precision gbm
      double precision gc
      double precision gca
      double precision gcab
      double precision gcb
      double precision gcbk
      double precision gm
      integer j
      integer k
      logical l0
      logical l1
      logical l2
      logical l3
      logical l4
      logical l5
      logical l6
      integer m
      integer mab
      integer mcab
      integer nca
      integer ncb
      integer nm
      double precision pa
      double precision pac
      double precision pb
      double precision pca
      double precision pi
      double precision rk1
      double precision rk2
      double precision rm
      double precision sj1
      double precision sj2
      double precision sm
      double precision sp
      double precision sp0
      double precision sq
      double precision t0
      double precision w0
      double precision ws
      double precision x
      double precision y
      complex*16 z
      complex*16 z00
      complex*16 z1
      complex*16 zc0
      complex*16 zc1
      complex*16 zf0
      complex*16 zf1
      complex*16 zhf
      complex*16 zp
      complex*16 zp0
      complex*16 zr
      complex*16 zr0
      complex*16 zr1
      complex*16 zw

      x = real ( z )
      y = dimag ( z )
      eps = 1.0D-15
      l0 = c .eq. int ( c ) .and. c .lt. 0.0D+00
      l1 = abs ( 1.0D+00 - x ) .lt. eps .and. y .eq. 0.0D+00 .and.
     &  c - a - b .le. 0.0D+00
      l2 = cdabs ( z + 1.0D+00 ) .lt. eps .and. 
     &  abs ( c - a + b - 1.0D+00 ) .lt. eps
      l3 = a .eq. int ( a ) .and. a .lt. 0.0D+00
      l4 = b .eq. int ( b ) .and. b .lt. 0.0D+00
      l5 = c - a .eq. int ( c - a ) .and. c - a .le. 0.0D+00
      l6 = c - b .eq. int ( c - b ) .and. c - b .le. 0.0D+00
      aa = a
      bb = b
      a0 = cdabs ( z )
      if ( 0.95D+00 .lt. a0 ) then
        eps = 1.0D-08
      end if
      pi = 3.141592653589793D+00
      el = 0.5772156649015329D+00

      if ( l0 .or. l1 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'HYGFZ - Fatal error!'
        write ( *, '(a)' ) '  The hypergeometric series is divergent.'
        stop
      end if

      if ( a0 .eq. 0.0D+00 .or. a .eq. 0.0D+00 .or. 
     &  b .eq. 0.0D+00 ) then

        zhf = cmplx ( 1.0D+00, 0.0D+00 )

      else if ( z .eq. 1.0D+00.and. 0.0D+00 .lt. c - a - b ) then

        call gamma ( c, gc )
        call gamma ( c - a - b, gcab )
        call gamma ( c - a, gca )
        call gamma ( c - b, gcb )
        zhf = gc * gcab / ( gca * gcb )

      else if ( l2 ) then

        g0 = dsqrt ( pi ) * 2.0D+00 ** ( - a )
        call gamma ( c, g1 )
        call gamma ( 1.0D+00 + a / 2.0D+00 - b, g2 )
        call gamma ( 0.5D+00 + 0.5D+00 * a, g3 )
        zhf = g0 * g1 / ( g2 * g3 )

      else if ( l3 .or. l4 ) then

        if ( l3 ) then
          nm = int ( abs ( a ) )
        end if

        if ( l4 ) then
          nm = int ( abs ( b ) )
        end if

        zhf = cmplx ( 1.0D+00, 0.0D+00 )
        zr = cmplx ( 1.0D+00, 0.0D+00 )
        do k = 1, nm
          zr = zr * ( a + k - 1.0D+00 ) * ( b + k - 1.0D+00 )
     &      / ( k * ( c + k - 1.0D+00 ) ) * z
          zhf = zhf+zr
        end do

      else if ( l5 .or. l6 ) then

        if ( l5 ) then
          nm = int ( abs ( c - a ) )
        end if

        if ( l6 ) then
          nm = int ( abs ( c - b ) )
        end if

        zhf = cmplx ( 1.0D+00, 0.0D+00 )
        zr = cmplx ( 1.0D+00, 0.0D+00 )
        do k = 1, nm
          zr = zr * ( c - a + k - 1.0D+00 ) * ( c - b + k - 1.0D+00 )
     &      / ( k * ( c + k - 1.0D+00 ) ) * z
          zhf = zhf + zr
        end do
        zhf = ( 1.0D+00 - z ) ** ( c - a - b ) * zhf

      else if ( a0 .le. 1.0D+00 ) then

        if ( x .lt. 0.0D+00 ) then

          z1 = z / ( z - 1.0D+00 )
          if ( a .lt. c .and. b .lt. a .and. 0.0D+00 .lt. b ) then  
            a = bb
            b = aa
          end if
          zc0 = 1.0D+00 / ( ( 1.0D+00 - z ) ** a )
          zhf = cmplx ( 1.0D+00, 0.0D+00 )
          zr0 = cmplx ( 1.0D+00, 0.0D+00 )
          do k = 1, 500
            zr0 = zr0 * ( a + k - 1.0D+00 ) * ( c - b + k - 1.0D+00 )
     &        / ( k * ( c + k - 1.0D+00 ) ) * z1
            zhf = zhf + zr0
            if ( cdabs ( zhf - zw ) .lt. cdabs ( zhf ) * eps ) then
              go to 10
            end if
            zw = zhf
          end do

10        continue

          zhf = zc0 * zhf

        else if ( 0.90D+00 .le. a0 ) then

          gm = 0.0D+00
          mcab = int ( c - a - b + eps * dsign ( 1.0D+00, c - a - b ) )

          if ( abs ( c - a - b - mcab ) .lt. eps ) then

            m = int ( c - a - b )
            call gamma ( a, ga )
            call gamma ( b, gb )
            call gamma ( c, gc )
            call gamma ( a + m, gam )
            call gamma ( b + m, gbm ) 
            call psi ( a, pa )
            call psi ( b, pb )
            if ( m .ne. 0 ) then
              gm = 1.0D+00
            end if
            do j = 1, abs ( m ) - 1
              gm = gm * j
            end do
            rm = 1.0D+00
            do j = 1, abs ( m )
              rm = rm * j
            end do
            zf0 = cmplx ( 1.0D+00, 0.0D+00 )
            zr0 = cmplx ( 1.0D+00, 0.0D+00 )
            zr1 = cmplx ( 1.0D+00, 0.0D+00 )
            sp0 = 0.0D+00
            sp = 0.0D+00

            if ( 0 .le. m ) then

              zc0 = gm * gc / ( gam * gbm )
              zc1 = - gc * ( z - 1.0D+00 ) ** m / ( ga * gb * rm )
              do k = 1, m - 1
                zr0 = zr0 * ( a + k - 1.0D+00 ) * ( b + k - 1.0D+00 )
     &            / ( k * ( k - m ) ) * ( 1.0D+00 - z )
                zf0 = zf0 + zr0
              end do
              do k = 1, m
                sp0 = sp0 + 1.0D+00 / ( a + k - 1.0D+00 ) 
     &            + 1.0D+00 / ( b + k - 1.0D+00 ) - 1.0D+00 / k
              end do
              zf1 = pa + pb + sp0 + 2.0D+00 * el 
     &          + cdlog ( 1.0D+00 - z )
              do k = 1, 500
                sp = sp + ( 1.0D+00 - a )
     &            / ( k * ( a + k - 1.0D+00 ) ) + ( 1.0D+00 - b )
     &            / ( k * ( b + k - 1.0D+00 ) )
                sm = 0.0D+00
                do j = 1, m
                  sm = sm + ( 1.0D+00 - a ) / ( ( j + k ) 
     &              * ( a + j + k - 1.0D+00 ) )
     &              + 1.0D+00 / ( b + j + k - 1.0D+00 )
                end do
                zp = pa + pb + 2.0D+00 * el + sp + sm 
     &            + cdlog ( 1.0D+00 - z )
                zr1 = zr1 * ( a + m + k - 1.0D+00 )
     &            * ( b + m + k - 1.0D+00 ) / ( k * ( m + k ) )
     &            * ( 1.0D+00 - z )
                zf1 = zf1 + zr1 * zp
                if ( cdabs ( zf1 - zw ) .lt. cdabs ( zf1 ) * eps ) then
                  go to 20
                end if
                zw = zf1
              end do

20            continue

              zhf = zf0 * zc0 + zf1 * zc1

            else if ( m .lt. 0 ) then

              m = - m
              zc0 = gm * gc / ( ga * gb * ( 1.0D+00 - z ) ** m )
              zc1 = - ( - 1.0D+00 ) ** m * gc / ( gam * gbm * rm )
              do k = 1, m - 1
                zr0 = zr0 * ( a - m + k - 1.0D+00 ) 
     &            * ( b - m + k - 1.0D+00 ) / ( k * ( k - m ) )
     &            * ( 1.0D+00 - z )
                zf0 = zf0 + zr0
              end do

              do k = 1, m
                sp0 = sp0 + 1.0D+00 / k
              end do

              zf1 = pa + pb - sp0 + 2.0D+00 * el + cdlog ( 1.0D+00 - z )

              do k = 1, 500
                sp = sp + ( 1.0D+00 - a ) / ( k * ( a + k - 1.0D+00 ) )
     &            + ( 1.0D+00 - b ) / ( k * ( b + k - 1.0D+00 ) )
                sm = 0.0D+00
                do j = 1, m
                  sm = sm + 1.0D+00 / ( j + k )
                end do
                zp = pa + pb + 2.0D+00 * el + sp - sm 
     &            + cdlog ( 1.0D+00 - z )
                zr1 = zr1 * ( a + k - 1.0D+00 ) * ( b + k - 1.0D+00 )
     &            / ( k * ( m + k ) ) * ( 1.0D+00 - z )
                zf1 = zf1 + zr1 * zp
                if ( cdabs ( zf1 - zw ) .lt. cdabs ( zf1 ) * eps ) then
                  go to 30
                end if
                zw = zf1

              end do

30            continue

              zhf = zf0 * zc0 + zf1 * zc1

            end if

          else

            call gamma ( a, ga )
            call gamma ( b, gb )
            call gamma ( c, gc )
            call gamma ( c - a, gca )
            call gamma ( c - b, gcb )
            call gamma ( c - a - b, gcab )
            call gamma ( a + b - c, gabc )
            zc0 = gc * gcab / ( gca * gcb )
            zc1 = gc * gabc / ( ga * gb ) 
     &        * ( 1.0D+00 - z ) ** ( c - a - b )
            zhf = cmplx ( 0.0D+00, 0.0D+00 )
            zr0 = zc0
            zr1 = zc1
            do k = 1, 500
              zr0 = zr0 * ( a + k - 1.0D+00 ) * ( b + k - 1.0D+00 )
     &          / ( k * ( a + b - c + k ) ) * ( 1.0D+00 - z )
              zr1 = zr1 * ( c - a + k - 1.0D+00 ) 
     &          * ( c - b + k - 1.0D+00 ) / ( k * ( c - a - b + k ) )
     &          * ( 1.0D+00 - z )
              zhf = zhf + zr0 + zr1
              if ( cdabs ( zhf - zw ) .lt. cdabs ( zhf ) * eps ) then
                go to 40
              end if
              zw = zhf
            end do

40          continue

            zhf = zhf + zc0 + zc1

          end if

        else

          z00 = cmplx ( 1.0D+00, 0.0D+00 )

          if ( c - a .lt. a .and. c - b .lt. b ) then
            z00 = ( 1.0D+00 - z ) ** ( c - a - b )
            a = c - a
            b = c - b
          end if

          zhf = cmplx ( 1.0D+00, 0.0D+00 )
          zr = cmplx ( 1.0D+00, 0.0D+00 )

          do k = 1, 1500
            zr = zr * ( a + k - 1.0D+00 ) * ( b + k - 1.0D+00 ) 
     &        / ( k * ( c + k - 1.0D+00 ) ) * z
            zhf = zhf + zr
            if ( cdabs ( zhf - zw ) .le. cdabs ( zhf ) * eps ) then
              go to 50
            end if
            zw = zhf
          end do

50        continue

          zhf = z00 * zhf

        end if

      else if ( 1.0D+00 .lt. a0 ) then

        mab = int ( a - b + eps * dsign ( 1.0D+00, a - b ) )

        if ( abs ( a - b - mab ) .lt. eps .and. a0 .le. 1.1D+00 ) then
          b = b + eps
        end if

        if ( eps .lt. abs ( a - b - mab ) ) then

          call gamma ( a, ga )
          call gamma ( b, gb )
          call gamma ( c, gc )
          call gamma ( a - b, gab )
          call gamma ( b - a, gba )
          call gamma ( c - a, gca )
          call gamma ( c - b, gcb )
          zc0 = gc * gba / ( gca * gb * ( - z ) ** a )
          zc1 = gc * gab / ( gcb * ga * ( - z ) ** b )
          zr0 = zc0
          zr1 = zc1
          zhf = cmplx ( 0.0D+00, 0.0D+00 )

          do k = 1, 500
            zr0 = zr0 * ( a + k - 1.0D+00 ) * ( a - c + k ) 
     &        / ( ( a - b + k ) * k * z )
            zr1 = zr1 * ( b + k - 1.0D+00 ) * ( b - c + k ) 
     &        / ( ( b - a + k ) * k * z )
            zhf = zhf + zr0 + zr1
            if ( cdabs ( ( zhf - zw ) / zhf ) .le. eps ) then
              go to 60
            end if
            zw = zhf
          end do

60        continue

          zhf = zhf + zc0 + zc1

        else

          if ( a - b .lt. 0.0D+00 ) then
            a = bb
            b = aa
          end if

          ca = c - a
          cb = c - b
          nca = int ( ca + eps * dsign ( 1.0D+00, ca ) )
          ncb = int ( cb + eps * dsign ( 1.0D+00, cb ) )

          if ( abs ( ca - nca ) .lt. eps .or. 
     &      abs ( cb - ncb ) .lt. eps ) then
            c = c + eps
          end if

          call gamma ( a, ga )
          call gamma ( c, gc )
          call gamma ( c - b, gcb )
          call psi ( a, pa )
          call psi ( c - a, pca )
          call psi ( a - c, pac )
          mab = int ( a - b + eps )
          zc0 = gc / ( ga * ( - z ) ** b )
          call gamma ( a - b, gm )
          zf0 = gm / gcb * zc0
          zr = zc0
          do k = 1, mab - 1
            zr = zr * ( b + k - 1.0D+00 ) / ( k * z )
            t0 = a - b - k
            call gamma ( t0, g0 )
            call gamma ( c - b - k, gcbk )
            zf0 = zf0 + zr * g0 / gcbk
          end do

          if ( mab .eq. 0 ) then
            zf0 = cmplx ( 0.0D+00, 0.0D+00 )
          end if

          zc1 = gc / ( ga * gcb * ( - z ) ** a )
          sp = -2.0D+00 * el - pa - pca
          do j = 1, mab
            sp = sp + 1.0D+00 / j
          end do
          zp0 = sp + cdlog ( - z )
          sq = 1.0D+00
          do j = 1, mab
            sq = sq * ( b + j - 1.0D+00 ) * ( b - c + j ) / j
          end do
          zf1 = ( sq * zp0 ) * zc1
          zr = zc1
          rk1 = 1.0D+00
          sj1 = 0.0D+00

          do k = 1, 10000
            zr = zr / z
            rk1 = rk1 * ( b + k - 1.0D+00 ) * ( b - c + k ) / ( k * k )
            rk2 = rk1
            do j = k + 1, k + mab
              rk2 = rk2 * ( b + j - 1.0D+00 ) * ( b - c + j ) / j
            end do
            sj1 = sj1 + ( a - 1.0D+00 ) / ( k * ( a + k - 1.0D+00 ) )
     &        + ( a - c - 1.0D+00 ) / ( k * ( a - c + k - 1.0D+00 ) )
            sj2 = sj1
            do j = k + 1, k + mab
              sj2 = sj2 + 1.0D+00 / j
            end do
            zp = -2.0D+00 * el - pa - pac + sj2 
     &        - 1.0D+00 / ( k + a - c )
     &        - pi / dtan ( pi * ( k + a - c ) ) + cdlog ( - z ) 
            zf1 = zf1 + rk2 * zr * zp
            ws = cdabs ( zf1 )
            if ( abs ( ( ws - w0 ) / ws ) .lt. eps ) then
              go to 70
            end if
            w0 = ws
          end do

70        continue

          zhf = zf0 + zf1

        end if

      end if

      a = aa
      b = bb
      if ( 150 .lt. k ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'HYGFZ - Warning!'
        write ( *, '(a)' ) 
     &    '  The solution returned may have low accuracy.'
      end if

      return
      end
      subroutine ik01a ( x, bi0, di0, bi1, di1, bk0, dk0, bk1, dk1 )

c*********************************************************************72
c
cc IK01A compute Bessel function I0(x), I1(x), K0(x), and K1(x).
c
c  Discussion:
c
c    This procedure computes modified Bessel functions I0(x), I1(x),
c    K0(x) and K1(x), and their derivatives.
c
c  Licensing:
c
c    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
c    they give permission to incorporate this routine into a user program 
c    provided that the copyright is acknowledged.
c
c  Modified:
c
c    16 July 2012
c
c  Author:
c
c    Shanjie Zhang, Jianming Jin
c
c  Reference:
c
c    Shanjie Zhang, Jianming Jin,
c    Computation of Special Functions,
c    Wiley, 1996,
c    ISBN: 0-471-11963-6,
c    LC: QA351.C45.
c
c  Parameters:
c
c    Input, double precision X, the argument.
c
c    Output, double precision BI0, DI0, BI1, DI1, BK0, DK0, BK1, DK1, the
c    values of I0(x), I0'(x), I1(x), I1'(x), K0(x), K0'(x), K1(x), K1'(x).
c
      implicit none

      double precision a(12)
      double precision a1(8)
      double precision b(12)
      double precision bi0
      double precision bi1
      double precision bk0
      double precision bk1
      double precision ca
      double precision cb
      double precision ct
      double precision di0
      double precision di1
      double precision dk0
      double precision dk1
      double precision el
      integer k
      integer k0
      double precision pi
      double precision r
      double precision w0
      double precision ww
      double precision x
      double precision x2
      double precision xr
      double precision xr2

      save a
      save a1
      save b

      data a / 0.125D+00, 7.03125D-02,
     &            7.32421875D-02, 1.1215209960938D-01,
     &            2.2710800170898D-01, 5.7250142097473D-01,
     &            1.7277275025845D+00, 6.0740420012735D+00,
     &            2.4380529699556D+01, 1.1001714026925D+002,
     &            5.5133589612202D+02, 3.0380905109224D+03/

      data a1 / 0.125D+00, 0.2109375D+00,
     &             1.0986328125D+00, 1.1775970458984D+01,
     &             2.1461706161499D+02, 5.9511522710323D+03,
     &             2.3347645606175D+05, 1.2312234987631D+07/

      data b / -0.375D+00, -1.171875D-01,
     &            -1.025390625D-01, -1.4419555664063D-01,
     &            -2.7757644653320D-01, -6.7659258842468D-01,
     &            -1.9935317337513D+00, -6.8839142681099D+00,
     &            -2.7248827311269D+01, -1.2159789187654D+02,
     &            -6.0384407670507D+02, -3.3022722944809D+03/

      pi = 3.141592653589793D+00
      el = 0.5772156649015329D+00
      x2 = x * x

      if ( x .eq. 0.0D+00 ) then

        bi0 = 1.0D+00
        bi1 = 0.0D+00
        bk0 = 1.0D+300
        bk1 = 1.0D+300
        di0 = 0.0D+00
        di1 = 0.5D+00
        dk0 = -1.0D+300
        dk1 = -1.0D+300
        return

      else if ( x .le. 18.0D+00 ) then

        bi0 = 1.0D+00
        r = 1.0D+00
        do k = 1, 50
          r = 0.25D+00 * r * x2 / ( k * k )
          bi0 = bi0 + r
          if ( abs ( r / bi0 ) .lt. 1.0D-15 ) then
            go to 10
          end if
        end do

10      continue

        bi1 = 1.0D+00
        r = 1.0D+00
        do k = 1, 50
          r = 0.25D+00 * r * x2 / ( k * ( k + 1 ) )
          bi1 = bi1 + r
          if ( abs ( r / bi1 ) .lt. 1.0D-15 ) then
            go to 20
          end if
        end do

20      continue

        bi1 = 0.5D+00 * x * bi1

      else

        if ( x .lt. 35.0D+00 ) then
          k0 = 12
        else if ( x .lt. 50.0D+00 ) then
          k0 = 9
        else
          k0 = 7
        end if

        ca = exp ( x ) / sqrt ( 2.0D+00 * pi * x )
        bi0 = 1.0D+00
        xr = 1.0D+00 / x
        do k = 1, k0
          bi0 = bi0 + a(k) * xr ** k
        end do
        bi0 = ca * bi0
        bi1 = 1.0D+00
        do k = 1, k0
          bi1 = bi1 + b(k) * xr ** k
        end do
        bi1 = ca * bi1

      end if

      if ( x .le. 9.0D+00 ) then

        ct = - ( log ( x / 2.0D+00 ) + el )
        bk0 = 0.0D+00
        w0 = 0.0D+00
        r = 1.0D+00
        do k = 1, 50
          w0 = w0 + 1.0D+00 / k
          r = 0.25D+00 * r / ( k * k ) * x2
          bk0 = bk0 + r * ( w0 + ct )
          if ( abs ( ( bk0 - ww ) / bk0 ) .lt. 1.0D-15 ) then
            go to 30
          end if
          ww = bk0
        end do

30      continue

        bk0 = bk0 + ct

      else

        cb = 0.5D+00 / x
        xr2 = 1.0D+00 / x2
        bk0 = 1.0D+00
        do k = 1, 8
          bk0 = bk0 + a1(k) * xr2 ** k
        end do
        bk0 = cb * bk0 / bi0

      end if

      bk1 = ( 1.0D+00 / x - bi1 * bk0 ) / bi0
      di0 = bi1
      di1 = bi0 - bi1 / x
      dk0 = - bk1
      dk1 = - bk0 - bk1 / x

      return
      end
      subroutine ik01b ( x, bi0, di0, bi1, di1, bk0, dk0, bk1, dk1 )

c*********************************************************************72
c
cc IK01B: Bessel functions I0(x), I1(x), K0(x), and K1(x) and derivatives.
c
c  Licensing:
c
c    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
c    they give permission to incorporate this routine into a user program 
c    provided that the copyright is acknowledged.
c
c  Modified:
c
c    17 July 2012
c
c  Author:
c
c    Shanjie Zhang, Jianming Jin
c
c  Reference:
c
c    Shanjie Zhang, Jianming Jin,
c    Computation of Special Functions,
c    Wiley, 1996,
c    ISBN: 0-471-11963-6,
c    LC: QA351.C45.
c
c  Parameters:
c
c    Input, double precision X, the argument.
c
c    Output, double precision BI0, DI0, BI1, DI1, BK0, DK0, BK1, DK1, the
c    values of I0(x), I0'(x), I1(x), I1'(x), K0(x), K0'(x), K1(x), K1'(x).
c
      implicit none

      double precision bi0
      double precision bi1
      double precision bk0
      double precision bk1
      double precision di0
      double precision di1
      double precision dk0
      double precision dk1
      double precision t
      double precision t2
      double precision x

      if ( x .eq. 0.0D+00 ) then

        bi0 = 1.0D+00
        bi1 = 0.0D+00
        bk0 = 1.0D+300
        bk1 = 1.0D+300
        di0 = 0.0D+00
        di1 = 0.5D+00
        dk0 = -1.0D+300
        dk1 = -1.0D+300
        return

      else if ( x .le. 3.75D+00 ) then

        t = x / 3.75D+00
        t2 = t * t

        bi0 = ((((( 
     &      0.0045813D+00   * t2
     &    + 0.0360768D+00 ) * t2
     &    + 0.2659732D+00 ) * t2
     &    + 1.2067492D+00 ) * t2
     &    + 3.0899424D+00 ) * t2
     &    + 3.5156229D+00 ) * t2
     &    + 1.0D+00

        bi1 = x * ((((((
     &      0.00032411D+00   * t2
     &    + 0.00301532D+00 ) * t2
     &    + 0.02658733D+00 ) * t2
     &    + 0.15084934D+00 ) * t2
     &    + 0.51498869D+00 ) * t2
     &    + 0.87890594D+00 ) * t2
     &    + 0.5D+00 )

      else

        t = 3.75D+00 / x

        bi0 = ((((((((
     &      0.00392377D+00   * t
     &    - 0.01647633D+00 ) * t
     &    + 0.02635537D+00 ) * t
     &    - 0.02057706D+00 ) * t
     &    + 0.916281D-02 ) * t
     &    - 0.157565D-02 ) * t
     &    + 0.225319D-02 ) * t
     &    + 0.01328592D+00 ) * t
     &    + 0.39894228D+00 ) * exp ( x ) / sqrt ( x )

        bi1 = ((((((((
     &    - 0.420059D-02     * t
     &    + 0.01787654D+00 ) * t
     &    - 0.02895312D+00 ) * t
     &    + 0.02282967D+00 ) * t
     &    - 0.01031555D+00 ) * t
     &    + 0.163801D-02 ) * t
     &    - 0.00362018D+00 ) * t
     &    - 0.03988024D+00 ) * t
     &    + 0.39894228D+00 ) * exp ( x ) / sqrt ( x )

      end if

      if ( x .le. 2.0D+00 ) then

        t = x / 2.0D+00
        t2 = t * t

        bk0 = (((((
     &      0.0000074D+00   * t2
     &    + 0.0001075D+00 ) * t2
     &    + 0.00262698D+00 ) * t2
     &    + 0.0348859D+00 ) * t2
     &    + 0.23069756D+00 ) * t2
     &    + 0.4227842D+00 ) * t2
     &    - 0.57721566D+00 - bi0 * log ( t )

        bk1 = ((((((
     &    - 0.00004686D+00   * t2
     &    - 0.00110404D+00 ) * t2
     &    - 0.01919402D+00 ) * t2
     &    - 0.18156897D+00 ) * t2
     &    - 0.67278579D+00 ) * t2
     &    + 0.15443144D+00 ) * t2
     &    + 1.0D+00 ) / x + bi1 * log ( t )

      else

        t = 2.0D+00 / x
        t2 = t * t

        bk0 = ((((((
     &      0.00053208D+00   * t
     &    - 0.0025154D+00 )  * t
     &    + 0.00587872D+00 ) * t
     &    - 0.01062446D+00 ) * t
     &    + 0.02189568D+00 ) * t
     &    - 0.07832358D+00 ) * t
     &    + 1.25331414D+00 ) * exp ( - x ) / sqrt ( x )

        bk1 = ((((((
     &    - 0.00068245D+00   * t 
     &    + 0.00325614D+00 ) * t
     &    - 0.00780353D+00 ) * t
     &    + 0.01504268D+00 ) * t
     &    - 0.0365562D+00  ) * t 
     &    + 0.23498619D+00 ) * t
     &    + 1.25331414D+00 ) * exp ( - x ) / sqrt ( x )

      end if

      di0 = bi1
      di1 = bi0 - bi1 / x
      dk0 = -bk1
      dk1 = -bk0 - bk1 / x

      return
      end
      subroutine ikna ( n, x, nm, bi, di, bk, dk )

c*********************************************************************72
c
cc IKNA compute Bessel function In(x) and Kn(x), and derivatives.
c
c  Licensing:
c
c    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
c    they give permission to incorporate this routine into a user program 
c    provided that the copyright is acknowledged.
c
c  Modified:
c
c    16 July 2012
c
c  Author:
c
c    Shanjie Zhang, Jianming Jin
c
c  Reference:
c
c    Shanjie Zhang, Jianming Jin,
c    Computation of Special Functions,
c    Wiley, 1996,
c    ISBN: 0-471-11963-6,
c    LC: QA351.C45.
c
c  Parameters:
c
c    Input, integer N, the order of In(x) and Kn(x).
c
c    Input, double precision X, the argument.
c
c    Output, integer NM, the highest order computed.
c
c    Output, double precision BI(0:N), DI(0:N), BK(0:N), DK(0:N),
c    the values of In(x), In'(x), Kn(x), Kn'(x).
c
      implicit none

      integer n

      double precision bi(0:n)
      double precision bi0
      double precision bi1
      double precision bk(0:n)
      double precision bk0
      double precision bk1
      double precision di(0:n)
      double precision di0
      double precision di1
      double precision dk(0:n)
      double precision dk0
      double precision dk1
      double precision f
      double precision f0
      double precision f1
      double precision g
      double precision g0
      double precision g1
      double precision h
      double precision h0
      double precision h1
      integer k
      integer m
      integer msta1
      integer msta2
      integer nm
      double precision s0
      double precision x

      nm = n

      if ( x .le. 1.0D-100 ) then
        do k = 0, n
          bi(k) = 0.0D+00
          di(k) = 0.0D+00
          bk(k) = 1.0D+300
          dk(k) = -1.0D+300
        end do
        bi(0) = 1.0D+00
        di(1) = 0.5D+00
        return
      end if

      call ik01a ( x, bi0, di0, bi1, di1, bk0, dk0, bk1, dk1 )
      bi(0) = bi0
      bi(1) = bi1
      bk(0) = bk0
      bk(1) = bk1
      di(0) = di0
      di(1) = di1
      dk(0) = dk0
      dk(1) = dk1

      if ( n .le. 1 ) then
        return
      end if

      if ( 40.0D+00 .lt. x .and. n .lt. int ( 0.25D+00 * x ) ) then

        h0 = bi0
        h1 = bi1
        do k = 2, n
          h = -2.0D+00 * ( k - 1.0D+00 ) / x * h1 + h0
          bi(k) = h
          h0 = h1
          h1 = h
        end do

      else

        m = msta1 ( x, 200 )

        if ( m .lt. n ) then
          nm = m
        else
          m = msta2 ( x, n, 15 )
        end if

        f0 = 0.0D+00
        f1 = 1.0D-100
        do k = m, 0, -1
          f = 2.0D+00 * ( k + 1.0D+00 ) * f1 / x + f0
          if ( k .le. nm ) then
            bi(k) = f
          end if
          f0 = f1
          f1 = f
        end do
        s0 = bi0 / f
        do k = 0, nm
          bi(k) = s0 * bi(k)
        end do
      end if

      g0 = bk0
      g1 = bk1
      do k = 2, nm
        g = 2.0D+00 * ( k - 1.0D+00 ) / x * g1 + g0
        bk(k) = g
        g0 = g1
        g1 = g
      end do

      do k = 2, nm
        di(k) = bi(k-1) - k / x * bi(k)
        dk(k) = - bk(k-1) - k / x * bk(k)
      end do

      return
      end
      subroutine iknb ( n, x, nm, bi, di, bk, dk )

c*********************************************************************72
c
cc IKNB compute Bessel function In(x) and Kn(x).
c
c  Discussion:
c
c    Compute modified Bessel functions In(x) and Kn(x),
c    and their derivatives.
c
c  Licensing:
c
c    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
c    they give permission to incorporate this routine into a user program 
c    provided that the copyright is acknowledged.
c
c  Modified:
c
c    17 July 2012
c
c  Author:
c
c    Shanjie Zhang, Jianming Jin
c
c  Reference:
c
c    Shanjie Zhang, Jianming Jin,
c    Computation of Special Functions,
c    Wiley, 1996,
c    ISBN: 0-471-11963-6,
c    LC: QA351.C45.
c
c  Parameters:
c
c    Input, integer N, the order of In(x) and Kn(x).
c
c    Input, double precision X, the argument.
c
c    Output, integer NM, the highest order computed.
c
c    Output, double precision BI(0:N), DI(0:N), BK(0:N), DK(0:N),
c    the values of In(x), In'(x), Kn(x), Kn'(x).
c
      implicit none

      integer n

      double precision a0
      double precision bi(0:n)
      double precision bk(0:n)
      double precision bkl
      double precision bs
      double precision di(0:n)
      double precision dk(0:n)
      double precision el
      double precision f
      double precision f0
      double precision f1
      double precision g
      double precision g0
      double precision g1
      integer k
      integer k0
      integer l
      integer m
      integer msta1
      integer msta2
      integer nm
      double precision pi
      double precision r
      double precision s0
      double precision sk0
      double precision vt
      double precision x

      pi = 3.141592653589793D+00
      el = 0.5772156649015329D+00
      nm = n

      if ( x .le. 1.0D-100 ) then
        do k = 0, n
          bi(k) = 0.0D+00
          di(k) = 0.0D+00
          bk(k) = 1.0D+300
          dk(k) = -1.0D+300
        end do
        bi(0) = 1.0D+00
        di(1) = 0.5D+00
        return
      end if

      if ( n .eq. 0 ) then
        nm = 1
      end if

      m = msta1 ( x, 200 )
      if ( m .lt. nm ) then
        nm = m
      else
        m = msta2 ( x, nm, 15 )
      end if

      bs = 0.0D+00
      sk0 = 0.0D+00
      f0 = 0.0D+00
      f1 = 1.0D-100
      do k = m, 0, -1
        f = 2.0D+00 * ( k + 1.0D+00 ) / x * f1 + f0
        if ( k .le. nm ) then
          bi(k) = f
        end if
        if ( k .ne. 0 .and. k .eq. 2 * int ( k / 2 ) ) then
          sk0 = sk0 + 4.0D+00 * f / k
        end if
        bs = bs + 2.0D+00 * f
        f0 = f1
        f1 = f
      end do

      s0 = exp ( x ) / ( bs - f )
      do k = 0, nm
        bi(k) = s0 * bi(k)
      end do

      if ( x .le. 8.0D+00 ) then
        bk(0) = - ( log ( 0.5D+00 * x ) + el ) * bi(0) + s0 * sk0
        bk(1) = ( 1.0D+00 / x - bi(1) * bk(0) ) / bi(0)
      else
        a0 = sqrt ( pi / ( 2.0D+00 * x ) ) * exp ( - x ) 

        if ( x .lt. 25.0D+00 ) then
          k0 = 16
        else if ( x .lt. 80.0D+00 ) then
          k0 = 10
        else if ( x .lt. 200.0D+00 ) then
          k0 = 8
        else
          k0 = 6
        end if

        do l = 0, 1
          bkl = 1.0D+00
          vt = 4.0D+00 * l
          r = 1.0D+00
          do k = 1, k0
            r = 0.125D+00 * r * ( vt - ( 2.0D+00 * k - 1.0D+00 ) ** 2 ) 
     &        / ( k * x )
            bkl = bkl + r
          end do
          bk(l) = a0 * bkl
        end do
      end if

      g0 = bk(0)
      g1 = bk(1)
      do k = 2, nm
        g = 2.0D+00 * ( k - 1.0D+00 ) / x * g1 + g0
        bk(k) = g
        g0 = g1
        g1 = g
      end do

      di(0) = bi(1)
      dk(0) = -bk(1)
      do k = 1, nm
        di(k) = bi(k-1) - k / x * bi(k)
        dk(k) = -bk(k-1) - k / x * bk(k)
      end do

      return
      end
      subroutine ikv ( v, x, vm, bi, di, bk, dk )

c*********************************************************************72
c
cc IKV compute modified Bessel function Iv(x) and Kv(x) and their derivatives.
c
c  Licensing:
c
c    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
c    they give permission to incorporate this routine into a user program 
c    provided that the copyright is acknowledged.
c
c  Modified:
c
c    17 July 2012
c
c  Author:
c
c    Shanjie Zhang, Jianming Jin
c
c  Reference:
c
c    Shanjie Zhang, Jianming Jin,
c    Computation of Special Functions,
c    Wiley, 1996,
c    ISBN: 0-471-11963-6,
c    LC: QA351.C45.
c
c  Parameters:
c
c    Input, double precision V, the order of Iv(x) and Kv(x).
c    V = N + V0.
c
c    Input, double precision X, the argument.
c
c    Output, double precision VM, the highest order computed.
c
c    Output, double precision BI(0:N), DI(0:N), BK(0:N), DK(0:N), the
c    values of In+v0(x), In+v0'(x), Kn+v0(x), Kn+v0'(x).
c
      implicit none

      double precision a1
      double precision a2
      double precision bi(0:*)
      double precision bi0
      double precision bk(0:*)
      double precision bk0
      double precision bk1
      double precision bk2
      double precision ca
      double precision cb
      double precision cs
      double precision ct
      double precision di(0:*)
      double precision dk(0:*)
      double precision f
      double precision f1
      double precision f2
      double precision gan
      double precision gap
      integer k
      integer k0
      integer m
      integer msta1
      integer msta2
      integer n
      double precision pi
      double precision piv
      double precision r
      double precision r1
      double precision r2
      double precision sum
      double precision v
      double precision v0
      double precision v0n
      double precision v0p
      double precision vm
      double precision vt
      double precision w0
      double precision wa
      double precision ww
      double precision x
      double precision x2

      pi = 3.141592653589793D+00
      x2 = x * x
      n = int ( v )
      v0 = v - n
      if ( n .eq. 0 ) then
        n = 1
      end if

      if ( x .lt. 1.0D-100 ) then

        do k = 0, n
          bi(k) = 0.0D+00
          di(k) = 0.0D+00
          bk(k) = -1.0D+300
          dk(k) = 1.0D+300
        end do

        if ( v .eq. 0.0D+00 ) then
          bi(0) = 1.0D+00
          di(1) = 0.5D+00
        end if

        vm = v
        return

      end if

      piv = pi * v0
      vt = 4.0D+00 * v0 * v0

      if ( v0 .eq. 0.0D+00 ) then
        a1 = 1.0D+00
      else
        v0p = 1.0D+00 + v0
        call gamma ( v0p, gap )
        a1 = ( 0.5D+00 * x ) ** v0 / gap
      end if

      if ( x .lt. 35.0D+00 ) then
        k0 = 14
      else if ( x .lt. 50.0D+00 ) then
        k0 = 10
      else
        k0 = 8
      end if
 
      if ( x .le. 18.0D+00 ) then

        bi0 = 1.0D+00
        r = 1.0D+00
        do k = 1, 30
          r = 0.25D+00 * r * x2 / ( k * ( k + v0 ) )
          bi0 = bi0 + r
          if ( abs ( r / bi0 ) .lt. 1.0D-15 ) then
            go to 10
          end if
        end do

10      continue

        bi0 = bi0 * a1

      else

        ca = exp ( x ) / sqrt ( 2.0D+00 * pi * x )
        sum = 1.0D+00
        r = 1.0D+00
        do k = 1, k0
          r = -0.125D+00 * r * ( vt - ( 2.0D+00 * k - 1.0D+00 ) ** 2 ) 
     &      / ( k * x )
          sum = sum + r
        end do
        bi0 = ca * sum

      end if

      m = msta1 ( x, 200 )

      if ( m .lt. n ) then
        n = m
      else
        m = msta2 ( x, n, 15 )
      end if

      f2 = 0.0D+00
      f1 = 1.0D-100
      do k = m, 0, -1
        f = 2.0D+00 * ( v0 + k + 1.0D+00 ) / x * f1 + f2
        if ( k .le. n ) then
          bi(k) = f
        end if
        f2 = f1
        f1 = f
      end do

      cs = bi0 / f
      do k = 0, n
        bi(k) = cs * bi(k)
      end do

      di(0) = v0 / x * bi(0) + bi(1)
      do k = 1, n
        di(k) = - ( k + v0 ) / x * bi(k) + bi(k-1)
      end do

      if ( x .le. 9.0D+00 ) then

        if ( v0 .eq. 0.0D+00 ) then

          ct = - log ( 0.5D+00 * x ) - 0.5772156649015329D+00
          cs = 0.0D+00
          w0 = 0.0D+00
          r = 1.0D+00
          do k = 1, 50
            w0 = w0 + 1.0D+00 / k
            r = 0.25D+00 * r / ( k * k ) * x2
            cs = cs + r * ( w0 + ct )
            wa = abs ( cs )
            if ( abs ( ( wa - ww ) / wa ) .lt. 1.0D-15 ) then
              go to 20
            end if
            ww = wa
          end do

20        continue

          bk0 = ct + cs

        else

          v0n = 1.0D+00 - v0
          call gamma ( v0n, gan )
          a2 = 1.0D+00 / ( gan * ( 0.5D+00 * x ) ** v0 )
          a1 = ( 0.5D+00 * x ) ** v0 / gap
          sum = a2 - a1
          r1 = 1.0D+00
          r2 = 1.0D+00
          do k = 1, 120
            r1 = 0.25D+00 * r1 * x2 / ( k * ( k - v0 ) )
            r2 = 0.25D+00 * r2 * x2 / ( k * ( k + v0 ) )
            sum = sum + a2 * r1 - a1 * r2
            wa = abs ( sum )
            if ( abs ( ( wa - ww ) / wa ) .lt. 1.0D-15 ) then
              go to 30
            end if
            ww = wa
          end do

30        continue

          bk0 = 0.5D+00 * pi * sum / sin ( piv )

        end if

      else

        cb = exp ( - x ) * sqrt ( 0.5D+00 * pi / x )
        sum = 1.0D+00
        r = 1.0D+00
        do k = 1, k0
          r = 0.125D+00 * r * ( vt - ( 2.0D+00 * k - 1.0D+00 ) ** 2 ) 
     &      / ( k * x )
          sum = sum + r
        end do
        bk0 = cb * sum

      end if

      bk1 = ( 1.0D+00 / x - bi(1) * bk0 ) / bi(0)
      bk(0) = bk0
      bk(1) = bk1
      do k = 2, n
        bk2 = 2.0D+00 * ( v0 + k - 1.0D+00 ) / x * bk1 + bk0
        bk(k) = bk2
        bk0 = bk1
        bk1 = bk2
      end do

      dk(0) = v0 / x * bk(0) - bk(1)
      do k = 1, n
        dk(k) = - ( k + v0 ) / x * bk(k) - bk(k-1)
      end do

      vm = n + v0

      return
      end
      subroutine incob ( a, b, x, bix )

c*********************************************************************72
c
cc INCOB computes the incomplete beta function Ix(a,b).
c
c  Licensing:
c
c    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
c    they give permission to incorporate this routine into a user program 
c    provided that the copyright is acknowledged.
c
c  Modified:
c
c    22 July 2012
c
c  Author:
c
c    Shanjie Zhang, Jianming Jin
c
c  Reference:
c
c    Shanjie Zhang, Jianming Jin,
c    Computation of Special Functions,
c    Wiley, 1996,
c    ISBN: 0-471-11963-6,
c    LC: QA351.C45.
c
c  Parameters:
c
c    Input, double precision A, B, parameters.
c
c    Input, double precision X, the argument.
c
c    Output, double precision BIX, the function value.
c
      implicit none

      double precision a
      double precision b
      double precision bix
      double precision bt
      double precision dk(51)
      double precision fk(51)
      integer k
      double precision s0
      double precision t1
      double precision t2
      double precision ta
      double precision tb
      double precision x

      s0 = ( a + 1.0D+00 ) / ( a + b + 2.0D+00 )
      call beta ( a, b, bt )

      if ( x .le. s0 ) then

        do k = 1, 20
          dk(2*k) = k * ( b - k ) * x / 
     &      ( a + 2.0D+00 * k - 1.0D+00 ) / ( a + 2.0D+00 * k )
        end do

        do k = 0, 20
          dk(2*k+1) = - ( a + k ) * ( a + b + k ) * x 
     &      / ( a + 2.0D+00 * k ) / ( a + 2.0D+00 * k + 1.0D+00 )
        end do

        t1 = 0.0D+00
        do k = 20, 1, -1
          t1 = dk(k) / ( 1.0D+00 + t1 )
        end do
        ta = 1.0D+00 / ( 1.0D+00 + t1 )
        bix = x ** a * ( 1.0D+00 - x ) ** b / ( a * bt ) * ta

      else

        do k = 1, 20
          fk(2*k) = k * ( a - k ) * ( 1.0D+00 - x ) 
     &      / ( b + 2.0D+00 * k - 1.0D+00 ) / ( b + 2.0D+00 * k )
        end do

        do k = 0,20
          fk(2*k+1) = - ( b + k ) * ( a + b + k ) * ( 1.0D+00 - x )
     &      / ( b + 2.0D+00 * k ) / ( b + 2.0D+00 * k + 1.0D+00 )
        end do

        t2 = 0.0D+00
        do k = 20, 1, -1
          t2 = fk(k) / ( 1.0D+00 + t2 )
        end do
        tb = 1.0D+00 / ( 1.0D+00 + t2 )
        bix = 1.0D+00 - x ** a * ( 1.0D+00 - x ) ** b / ( b * bt ) * tb

      end if

      return
      end
      subroutine incog ( a, x, gin, gim, gip )

c*********************************************************************72
c
cc INCOG computes the incomplete gamma function r(a,x), ,(a,x), P(a,x).
c
c  Licensing:
c
c    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
c    they give permission to incorporate this routine into a user program 
c    provided that the copyright is acknowledged.
c
c  Modified:
c
c    22 July 2012
c
c  Author:
c
c    Shanjie Zhang, Jianming Jin
c
c  Reference:
c
c    Shanjie Zhang, Jianming Jin,
c    Computation of Special Functions,
c    Wiley, 1996,
c    ISBN: 0-471-11963-6,
c    LC: QA351.C45.
c
c  Parameters:
c
c    Input, double precision A, the parameter.
c
c    Input, double precision X, the argument.
c
c    Output, double precision GIN, GIM, GIP, the values of
c    r(a,x), ‚(a,x), P(a,x).
c
      implicit none

      double precision a
      double precision ga
      double precision gim
      double precision gin
      double precision gip
      integer k
      double precision r
      double precision s
      double precision t0
      double precision x
      double precision xam

      xam = -  x + a * log ( x )

      if ( 700.0D+00 .lt. xam .or. 170.0D+00 .lt. a ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'INCOG - Fatal error!'
        write ( *, '(a)' ) '  A and/or X is too large!'
        stop
      end if

      if ( x .eq. 0.0D+00 ) then

        gin = 0.0D+00
        call gamma ( a, ga )
        gim = ga
        gip = 0.0D+00

      else if ( x .le. 1.0D+00 + a ) then

        s = 1.0D+00 / a
        r = s
        do k = 1, 60
          r = r * x / ( a + k )
          s = s + r
          if ( abs ( r / s ) .lt. 1.0D-15 ) then
            go to 10
          end if
        end do

10      continue

        gin = exp ( xam ) * s
        call gamma ( a, ga )
        gip = gin / ga
        gim = ga - gin

      else if ( 1.0D+00 + a .lt. x ) then

        t0 = 0.0D+00
        do k = 60, 1, -1
          t0 = ( k - a ) / ( 1.0D+00 + k / ( x + t0 ) )
        end do
        gim = exp ( xam ) / ( x + t0 )
        call gamma ( a, ga )
        gin = ga - gim
        gip = 1.0D+00 - gim / ga

      end if

      return
      end
      subroutine itairy ( x, apt, bpt, ant, bnt )

c*********************************************************************72
c
cc ITAIRY computes the integrals of Airy functions.
c
c  Discussion:
c
c    Compute the integrals of Airy functions with respect to t,
c    from 0 and x.
c
c  Licensing:
c
c    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
c    they give permission to incorporate this routine into a user program 
c    provided that the copyright is acknowledged.
c
c  Modified:
c
c    19 July 2012
c
c  Author:
c
c    Shanjie Zhang, Jianming Jin
c
c  Reference:
c
c    Shanjie Zhang, Jianming Jin,
c    Computation of Special Functions,
c    Wiley, 1996,
c    ISBN: 0-471-11963-6,
c    LC: QA351.C45.
c
c  Parameters:
c
c    Input, double precision X, the upper limit of the integral.
c
c    Output, double precision APT, BPT, ANT, BNT, the integrals, from 0 to x,
c    of Ai(t), Bi(t), Ai(-t), and Bi(-t).
c
      implicit none

      double precision a(16)
      double precision ant
      double precision apt
      double precision bnt
      double precision bpt
      double precision c1
      double precision c2
      double precision eps
      double precision fx
      double precision gx
      integer k
      integer l
      double precision pi
      double precision q0
      double precision q1
      double precision q2
      double precision r
      double precision sr3
      double precision su1
      double precision su2
      double precision su3
      double precision su4
      double precision su5
      double precision su6
      double precision x
      double precision xe
      double precision xp6
      double precision xr1
      double precision xr2

      save a

      data a / 0.569444444444444D+00, 0.891300154320988D+00,
     &             0.226624344493027D+01, 0.798950124766861D+01,
     &             0.360688546785343D+02, 0.198670292131169D+03,
     &             0.129223456582211D+04, 0.969483869669600D+04,
     &             0.824184704952483D+05, 0.783031092490225D+06,
     &             0.822210493622814D+07, 0.945557399360556D+08,
     &             0.118195595640730D+10, 0.159564653040121D+11,
     &             0.231369166433050D+12, 0.358622522796969D+13/

      eps = 1.0D-15
      pi = 3.141592653589793D+00
      c1 = 0.355028053887817D+00
      c2 = 0.258819403792807D+00
      sr3 = 1.732050807568877D+00

      if ( x .eq. 0.0D+00 ) then

        apt = 0.0D+00
        bpt = 0.0D+00
        ant = 0.0D+00
        bnt = 0.0D+00

      else

        if ( abs ( x ) .le. 9.25D+00 ) then

          do l = 0, 1

            x = ( -1.0D+00 ) ** l * x
            fx = x
            r = x

            do k = 1, 40
              r = r * ( 3.0D+00 * k - 2.0D+00 ) 
     &          / ( 3.0D+00 * k + 1.0D+00 ) * x / ( 3.0D+00 * k )
     &          * x / ( 3.0D+00 * k - 1.0D+00 ) * x 
              fx = fx + r
              if ( abs ( r ) .lt. abs ( fx ) * eps ) then
                go to 10
              end if
            end do

10          continue

            gx = 0.5D+00 * x * x
            r = gx

            do k = 1, 40
              r = r * ( 3.0D+00 * k - 1.0D+00 ) 
     &          / ( 3.0D+00 * k + 2.0D+00 ) * x / ( 3.0D+00 * k ) * x 
     &          / ( 3.0D+00 * k + 1.0D+00 ) * x
              gx = gx + r
              if ( abs ( r ) .lt. abs ( gx ) * eps ) then
                go to 20
              end if
            end do

20          continue

            ant = c1 * fx - c2 * gx
            bnt = sr3 * ( c1 * fx + c2 * gx )

            if ( l .eq. 0 ) then
              apt = ant
              bpt = bnt
            else
              ant = -ant
              bnt = -bnt
              x = -x
            end if

          end do

        else

          q2 = 1.414213562373095D+00
          q0 = 0.3333333333333333D+00
          q1 = 0.6666666666666667D+00
          xe = x * sqrt ( x ) / 1.5D+00
          xp6 = 1.0D+00 / sqrt ( 6.0D+00 * pi * xe )
          su1 = 1.0D+00
          r = 1.0D+00
          xr1 = 1.0D+00 / xe
          do k = 1, 16
            r = - r * xr1
            su1 = su1 + a(k) * r
          end do
          su2 = 1.0D+00
          r = 1.0D+00
          do k = 1, 16
            r = r * xr1
            su2 = su2 + a(k) * r
          end do

          apt = q0 - exp ( - xe ) * xp6 * su1
          bpt = 2.0D+00 * exp ( xe ) * xp6 * su2
          su3 = 1.0D+00
          r = 1.0D+00
          xr2 = 1.0D+00 / ( xe * xe )
          do k = 1, 8
            r = - r * xr2
            su3 = su3 + a(2*k) * r
          end do
          su4 = a(1) * xr1
          r = xr1
          do k = 1, 7
            r = -r * xr2
            su4 = su4 + a(2*k+1) * r
          end do
          su5 = su3 + su4
          su6 = su3 - su4
          ant = q1 - q2 * xp6 * ( su5 * cos ( xe ) - su6 * sin ( xe ) )
          bnt = q2 * xp6 * ( su5 * sin ( xe ) + su6 * cos ( xe ) )

        end if

      end if

      return
      end
      subroutine itika ( x, ti, tk )

c*********************************************************************72
c
cc ITIKA computes the integral of the modified Bessel functions I0(t) and K0(t).
c
c  Discussion:
c
c    This procedure integrates modified Bessel functions I0(t) and
c    K0(t) with respect to t from 0 to x.
c
c  Licensing:
c
c    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
c    they give permission to incorporate this routine into a user program 
c    provided that the copyright is acknowledged.
c
c  Modified:
c
c    18 July 2012
c
c  Author:
c
c    Shanjie Zhang, Jianming Jin
c
c  Reference:
c
c    Shanjie Zhang, Jianming Jin,
c    Computation of Special Functions,
c    Wiley, 1996,
c    ISBN: 0-471-11963-6,
c    LC: QA351.C45.
c
c  Parameters:
c
c    Input, double precision X, the upper limit of the integral.
c
c    Output, double precision TI, TK, the integrals of I0(t) and K0(t)
c    from 0 to X.
c
      implicit none

      double precision a(10)
      double precision b1
      double precision b2
      double precision e0
      double precision el
      integer k
      double precision pi
      double precision r
      double precision rc1
      double precision rc2
      double precision rs
      double precision ti
      double precision tk
      double precision tw
      double precision x
      double precision x2

      save a

      data a / 0.625D+00, 1.0078125D+00,
     &       2.5927734375D+00, 9.1868591308594D+00,
     &       4.1567974090576D+01, 2.2919635891914D+02,
     &       1.491504060477D+03, 1.1192354495579D+04,
     &       9.515939374212D+04, 9.0412425769041D+05 /

      pi = 3.141592653589793D+00
      el = 0.5772156649015329D+00

      if ( x .eq. 0.0D+00 ) then

        ti = 0.0D+00
        tk = 0.0D+00
        return

      else if ( x .lt. 20.0D+00 ) then

        x2 = x * x
        ti = 1.0D+00
        r = 1.0D+00
        do k = 1, 50
          r = 0.25D+00 * r * ( 2 * k - 1.0D+00 ) / ( 2 * k + 1.0D+00 ) 
     &      / ( k * k ) * x2
          ti = ti + r
          if ( abs ( r / ti ) .lt. 1.0D-12 ) then
            go to 10
          end if
        end do

10      continue

        ti = ti * x

      else

        ti = 1.0D+00
        r = 1.0D+00
        do k = 1, 10
          r = r / x
          ti = ti + a(k) * r
        end do
        rc1 = 1.0D+00 / sqrt ( 2.0D+00 * pi * x )
        ti = rc1 * exp ( x ) * ti

      end if

      if ( x .lt. 12.0D+00 ) then

        e0 = el + log ( x / 2.0D+00 )
        b1 = 1.0D+00 - e0
        b2 = 0.0D+00
        rs = 0.0D+00
        r = 1.0D+00
        do k = 1, 50
          r = 0.25D+00 * r * ( 2 * k - 1.0D+00 ) 
     &      / ( 2 * k + 1.0D+00 ) / ( k * k ) * x2
          b1 = b1 + r * ( 1.0D+00 / ( 2 * k + 1 ) - e0 )
          rs = rs + 1.0D+00 / k
          b2 = b2 + r * rs
          tk = b1 + b2
          if ( abs ( ( tk - tw ) / tk ) .lt. 1.0D-12 ) then
            go to 20
          end if
          tw = tk
        end do

20      continue

        tk = tk * x

      else

        tk = 1.0D+00
        r = 1.0D+00
        do k = 1, 10
          r = -r / x
          tk = tk + a(k) * r
        end do 
        rc2 = sqrt ( pi / ( 2.0D+00 * x ) )
        tk = pi / 2.0D+00 - rc2 * tk * exp ( - x )

      end if

      return
      end
      subroutine itikb ( x, ti, tk )

c*********************************************************************72
c
cc ITIKB computes the integral of the Bessel functions I0(t) and K0(t).
c
c  Discussion:
c
c    This procedure integrates Bessel functions I0(t) and K0(t)
c    with respect to t from 0 to x.
c
c  Licensing:
c
c    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
c    they give permission to incorporate this routine into a user program 
c    provided that the copyright is acknowledged.
c
c  Modified:
c
c    24 July 2012
c
c  Author:
c
c    Shanjie Zhang, Jianming Jin
c
c  Reference:
c
c    Shanjie Zhang, Jianming Jin,
c    Computation of Special Functions,
c    Wiley, 1996,
c    ISBN: 0-471-11963-6,
c    LC: QA351.C45.
c 
c  Parameters:
c
c    Input, double precision X, the upper limit of the integral.
c
c    Output, double precision TI, TK, the integral of I0(t) and K0(t)
c    from 0 to X.
c
      implicit none

      double precision pi
      double precision t
      double precision t1
      double precision ti
      double precision tk
      double precision x

      pi = 3.141592653589793D+00

      if ( x .eq. 0.0D+00 ) then

        ti = 0.0D+00

      else if ( x .lt. 5.0D+00 ) then

        t1 = x / 5.0D+00
        t = t1 * t1
        ti = ((((((((
     &      0.59434D-03 * t
     &    + 0.4500642D-02 ) * t
     &    + 0.044686921D+00 ) * t
     &    + 0.300704878D+00 ) * t
     &    + 1.471860153D+00 ) * t
     &    + 4.844024624D+00 ) * t
     &    + 9.765629849D+00 ) * t
     &    +10.416666367D+00 ) * t
     &    + 5.0D+00 ) * t1

      else if ( 5.0D+00 .le. x .and. x .le. 8.0D+00 ) then

        t = 5.0D+00 / x
        ti = (((
     &    - 0.015166D+00 * t
     &    - 0.0202292D+00 ) * t
     &    + 0.1294122D+00 ) * t
     &    - 0.0302912D+00 ) * t
     &    + 0.4161224D+00
        ti = ti * exp ( x ) / sqrt ( x )

      else

        t = 8.0D+00 / x
        ti = (((((
     &    - 0.0073995D+00 * t
     &    + 0.017744D+00 ) * t
     &    - 0.0114858D+00 ) * t
     &    + 0.55956D-02 ) * t
     &    + 0.59191D-02 ) * t
     &    + 0.0311734D+00 ) * t
     &    + 0.3989423D+00
        ti = ti * exp ( x ) / sqrt ( x )

      end if

      if ( x .eq. 0.0D+00 ) then

        tk = 0.0D+00

      else if ( x .le. 2.0D+00 ) then

        t1 = x / 2.0D+00
        t = t1 * t1
        tk = ((((((
     &      0.116D-05        * t
     &    + 0.2069D-04 )     * t
     &    + 0.62664D-03 )    * t
     &    + 0.01110118D+00 ) * t
     &    + 0.11227902D+00 ) * t
     &    + 0.50407836D+00 ) * t
     &    + 0.84556868D+00 ) * t1
        tk = tk - log ( x / 2.0D+00 ) * ti

      else if ( 2.0D+00 .lt. x .and. x .le. 4.0D+00 ) then

        t = 2.0D+00 / x
        tk = (((
     &      0.0160395D+00   * t
     &    - 0.0781715D+00 ) * t
     &    + 0.185984D+00 )  * t
     &    - 0.3584641D+00 ) * t
     &    + 1.2494934D+00
        tk = pi / 2.0D+00 - tk * exp ( - x ) / sqrt ( x )

      else if ( 4.0D+00 .lt. x .and. x .le. 7.0D+00 ) then

        t = 4.0D+00 / x
        tk = (((((
     &      0.37128D-02 * t
     &    - 0.0158449D+00 ) * t
     &    + 0.0320504D+00 ) * t
     &    - 0.0481455D+00 ) * t
     &    + 0.0787284D+00 ) * t
     &    - 0.1958273D+00 ) * t
     &    + 1.2533141D+00
        tk = pi / 2.0D+00 - tk * exp ( - x ) / sqrt ( x )

      else

        t = 7.0D+00 / x
        tk = (((((
     &      0.33934D-03      * t
     &    - 0.163271D-02 )   * t
     &    + 0.417454D-02 )   * t
     &    - 0.933944D-02 )   * t
     &    + 0.02576646D+00 ) * t
     &    - 0.11190289D+00 ) * t
     &    + 1.25331414D+00
        tk = pi / 2.0D+00 - tk * exp ( - x ) / sqrt ( x )

      end if

      return
      end
      subroutine itjya ( x, tj, ty )

c*********************************************************************72
c
cc ITJYA computes integrals of Bessel functions J0(t) and Y0(t).
c
c  Discussion:
c
c    This procedure integrates Bessel functions J0(t) and Y0(t) with
c    respect to t from 0 to x.
c
c  Licensing:
c
c    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
c    they give permission to incorporate this routine into a user program 
c    provided that the copyright is acknowledged.
c
c  Modified:
c
c    25 July 2012
c
c  Author:
c
c    Shanjie Zhang, Jianming Jin
c
c  Reference:
c
c    Shanjie Zhang, Jianming Jin,
c    Computation of Special Functions,
c    Wiley, 1996,
c    ISBN: 0-471-11963-6,
c    LC: QA351.C45.
c
c  Parameters:
c
c    Input, double precision X, the upper limit of the integral.
c
c    Output, double precision TJ, TY, the integrals of J0(t) and Y0(t) 
c    from 0 to x.
c
      implicit none

      double precision a(18)
      double precision a0
      double precision a1
      double precision af
      double precision bf
      double precision bg
      double precision el
      double precision eps
      integer k
      double precision pi
      double precision r
      double precision r2
      double precision rc
      double precision rs
      double precision tj
      double precision ty
      double precision ty1
      double precision ty2
      double precision x
      double precision x2
      double precision xp

      pi = 3.141592653589793D+00
      el = 0.5772156649015329D+00
      eps = 1.0D-12

      if ( x .eq. 0.0D+00 ) then

        tj = 0.0D+00
        ty = 0.0D+00

      else if ( x .le. 20.0D+00 ) then

        x2 = x * x
        tj = x
        r = x
        do k = 1, 60
          r = -0.25D+00 * r * ( 2 * k - 1.0D+00 ) / ( 2 * k + 1.0D+00 ) 
     &      / ( k * k ) * x2
          tj = tj + r
          if ( abs ( r ) .lt. abs ( tj ) * eps ) then
            go to 10
          end if
        end do

10      continue

        ty1 = ( el + log ( x / 2.0D+00 ) ) * tj
        rs = 0.0D+00
        ty2 = 1.0D+00
        r = 1.0D+00

        do k = 1, 60
          r = -0.25D+00 * r * ( 2 * k - 1.0D+00 ) / ( 2 * k + 1.0D+00 ) 
     &      / ( k * k ) * x2
          rs = rs + 1.0D+00 / k
          r2 = r * ( rs + 1.0D+00 / ( 2.0D+00 * k + 1.0D+00 ) )
          ty2 = ty2 + r2
          if ( abs ( r2 ) .lt. abs ( ty2 ) * eps ) then
            go to 20
          end if
        end do

20      continue

        ty = ( ty1 - x * ty2 ) * 2.0D+00 / pi

      else

        a0 = 1.0D+00
        a1 = 5.0D+00 / 8.0D+00
        a(1) = a1

        do k = 1, 16
          af = ( ( 1.5D+00 * ( k + 0.5D+00 ) * ( k + 5.0D+00 / 6.0D+00 )
     &      * a1 - 0.5D+00 * ( k + 0.5D+00 ) * ( k + 0.5D+00 ) 
     &      * ( k - 0.5D+00 ) * a0 ) ) / ( k + 1.0D+00 )
          a(k+1) = af
          a0 = a1
          a1 = af
        end do

        bf = 1.0D+00
        r = 1.0D+00
        do k = 1, 8
          r = -r / ( x * x )
          bf = bf + a(2*k) * r
        end do
        bg = a(1) / x
        r = 1.0D+00 / x
        do k = 1, 8
          r = -r / ( x * x )
          bg = bg + a(2*k+1) * r
        end do
        xp = x + 0.25D+00 * pi
        rc = sqrt ( 2.0D+00 / ( pi * x ) )
        tj = 1.0D+00 - rc * ( bf * cos ( xp ) + bg * sin ( xp ) )
        ty = rc * ( bg * cos ( xp ) - bf * sin ( xp ) )

      end if

      return
      end
      subroutine itjyb ( x, tj, ty )

c*********************************************************************72
c
cc ITJYB computes integrals of Bessel functions J0(t) and Y0(t).
c
c  Discussion:
c
c    This procedure integrates Bessel functions J0(t) and Y0(t)
c    with respect to t from 0 to x.
c
c  Licensing:
c
c    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
c    they give permission to incorporate this routine into a user program 
c    provided that the copyright is acknowledged.
c
c  Modified:
c
c    25 July 2012
c
c  Author:
c
c    Shanjie Zhang, Jianming Jin
c
c  Reference:
c
c    Shanjie Zhang, Jianming Jin,
c    Computation of Special Functions,
c    Wiley, 1996,
c    ISBN: 0-471-11963-6,
c    LC: QA351.C45.
c
c  Parameters:
c
c    Input, double precision X, the upper limit of the integral.
c
c    Output, double precision TJ, TY, the integrals of J0(t) and Y0(t) 
c    from 0 to x.
c
      implicit none

      double precision f0
      double precision g0
      double precision pi
      double precision t
      double precision tj
      double precision ty
      double precision x
      double precision x1
      double precision xt

      pi = 3.141592653589793D+00

      if ( x .eq. 0.0D+00 ) then

        tj = 0.0D+00
        ty = 0.0D+00

      else if ( x .le. 4.0D+00 ) then

        x1 = x / 4.0D+00
        t = x1 * x1

        tj = (((((((
     &    - 0.133718D-03      * t
     &    + 0.2362211D-02 )   * t
     &    - 0.025791036D+00 ) * t
     &    + 0.197492634D+00 ) * t
     &    - 1.015860606D+00 ) * t
     &    + 3.199997842D+00 ) * t
     &    - 5.333333161D+00 ) * t
     &    + 4.0D+00 ) * x1

        ty = ((((((((
     &      0.13351D-04       * t
     &    - 0.235002D-03 )    * t
     &    + 0.3034322d-02 )   * t
     &    - 0.029600855D+00 ) * t
     &    + 0.203380298D+00 ) * t
     &    - 0.904755062D+00 ) * t
     &    + 2.287317974D+00 ) * t
     &    - 2.567250468D+00 ) * t
     &    + 1.076611469D+00 ) * x1

        ty = 2.0D+00 / pi * log ( x / 2.0D+00 ) * tj - ty

      else if ( x .le. 8.0D+00 ) then

        xt = x - 0.25D+00 * pi
        t = 16.0D+00 / ( x * x )

        f0 = ((((((
     &      0.1496119D-02     * t
     &    - 0.739083D-02 )    * t
     &    + 0.016236617D+00 ) * t
     &    - 0.022007499D+00 ) * t
     &    + 0.023644978D+00 ) * t
     &    - 0.031280848D+00 ) * t
     &    + 0.124611058D+00 ) * 4.0D+00 / x

        g0 = (((((
     &      0.1076103D-02     * t
     &    - 0.5434851D-02 )   * t
     &    + 0.01242264D+00 )  * t
     &    - 0.018255209D+00 ) * t
     &    + 0.023664841D+00 ) * t
     &    - 0.049635633D+00 ) * t
     &    + 0.79784879D+00

        tj = 1.0D+00 
     &    - ( f0 * cos ( xt ) - g0 * sin ( xt ) ) / sqrt ( x )

        ty = - ( f0 * sin ( xt ) + g0 * cos ( xt ) ) / sqrt ( x )

      else

        t = 64.0D+00 / ( x * x )
        xt = x-0.25D+00 * pi

        f0 = (((((((
     &    - 0.268482D-04     * t
     &    + 0.1270039D-03 )  * t
     &    - 0.2755037D-03 )  * t
     &    + 0.3992825D-03 )  * t
     &    - 0.5366169D-03 )  * t
     &    + 0.10089872D-02 ) * t
     &    - 0.40403539D-02 ) * t
     &    + 0.0623347304D+00 ) * 8.0D+00 / x

        g0 = ((((((
     &    - 0.226238D-04        * t
     &    + 0.1107299D-03 )     * t
     &    - 0.2543955D-03 )     * t
     &    + 0.4100676D-03 )     * t
     &    - 0.6740148D-03 )     * t
     &    + 0.17870944D-02 )    * t
     &    - 0.01256424405D+00 ) * t
     &    + 0.79788456D+00

        tj = 1.0D+00 
     &    - ( f0 * cos ( xt ) - g0 * sin ( xt ) ) / sqrt ( x )

        ty = - ( f0 * sin ( xt ) + g0 * cos ( xt ) ) / sqrt ( x )

      end if

      return
      end
      subroutine itsh0 ( x, th0 )

c*********************************************************************72
c
cc ITSH0 integrates the Struve function H0(t) from 0 to x.
c
c  Discussion:
c
c    This procedure evaluates the integral of Struve function
c    H0(t) with respect to t from 0 and x.
c
c  Licensing:
c
c    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
c    they give permission to incorporate this routine into a user program 
c    provided that the copyright is acknowledged.
c
c  Modified:
c
c    25 July 2012
c
c  Author:
c
c    Shanjie Zhang, Jianming Jin
c
c  Reference:
c
c    Shanjie Zhang, Jianming Jin,
c    Computation of Special Functions,
c    Wiley, 1996,
c    ISBN: 0-471-11963-6,
c    LC: QA351.C45.
c
c  Parameters:
c
c    Input, double precision X, the upper limit of the integral.
c
c    Output, double precision TH0, the integral of H0(t) from 0 to x.
c
      implicit none

      double precision a(25)
      double precision a0
      double precision a1
      double precision af
      double precision bf
      double precision bg
      double precision el
      integer k
      double precision pi
      double precision r
      double precision rd
      double precision s
      double precision s0
      double precision th0
      double precision ty
      double precision x
      double precision xp

      pi = 3.141592653589793D+00
      r = 1.0D+00            

      if ( x .le. 30.0D+00 ) then

        s = 0.5D+00

        do k = 1, 100

          if ( k .eq. 1 ) then
            rd = 0.5D+00
          else
            rd = 1.0D+00
          end if

          r = - r * rd * k / ( k + 1.0D+00 ) 
     &      * ( x / ( 2.0D+00 * k + 1.0D+00 ) ) ** 2
          s = s + r

          if ( abs ( r ) .lt. abs ( s ) * 1.0D-12 ) then
            go to 10
          end if

        end do

10      continue

        th0 = 2.0D+00 / pi * x * x * s

      else

        s = 1.0D+00
        do k = 1, 12
          r = - r * k / ( k + 1.0D+00 ) 
     &      * ( ( 2.0D+00 * k + 1.0D+00 ) / x ) ** 2
          s = s + r
          if ( abs ( r ) .lt. abs ( s ) * 1.0D-12 ) then
            go to 20
          end if
        end do

20      continue

        el = 0.57721566490153D+00
        s0 = s / ( pi * x * x ) + 2.0D+00 / pi 
     &    * ( log ( 2.0D+00 * x ) + el )
        a0 = 1.0D+00
        a1 = 5.0D+00 / 8.0D+00
        a(1) = a1
        do k = 1, 20
          af = ( ( 1.5D+00 * ( k + 0.5D+00 ) 
     &      * ( k + 5.0D+00 / 6.0D+00 ) * a1 - 0.5D+00
     &      * ( k + 0.5D+00 ) * ( k + 0.5D+00 ) 
     &      * ( k - 0.5D+00 ) * a0 ) ) / ( k + 1.0D+00 )
          a(k+1) = af
          a0 = a1
          a1 = af
        end do

        bf = 1.0D+00
        r = 1.0D+00
        do k = 1, 10
          r = - r / ( x * x )
          bf = bf + a(2*k) * r
        end do
        bg = a(1) / x
        r = 1.0D+00 / x
        do k = 1, 10
          r = - r / ( x * x ) 
          bg = bg + a(2*k+1) * r
        end do
        xp = x + 0.25D+00 * pi
        ty = sqrt ( 2.0D+00 / ( pi * x ) ) 
     &    * ( bg * cos ( xp ) - bf * sin ( xp ) )
        th0 = ty + s0

      end if

      return
      end
      subroutine itsl0 ( x, tl0 )

c*********************************************************************72
c
cc ITSL0 integrates the Struve function L0(t) from 0 to x.
c
c  Discussion:
c
c    This procedure evaluates the integral of modified Struve function
c    L0(t) with respect to t from 0 to x.
c
c  Licensing:
c
c    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
c    they give permission to incorporate this routine into a user program 
c    provided that the copyright is acknowledged.
c
c  Modified:
c
c    31 July 2012
c
c  Author:
c
c    Shanjie Zhang, Jianming Jin
c
c  Reference:
c
c    Shanjie Zhang, Jianming Jin,
c    Computation of Special Functions,
c    Wiley, 1996,
c    ISBN: 0-471-11963-6,
c    LC: QA351.C45.
c
c  Parameters:
c
c    Input, double precision X, the upper limit of the integral.
c
c    Output, double precision TL0, the integral of L0(t) from 0 to x.
c
      implicit none

      double precision a(18)
      double precision a0
      double precision a1
      double precision af
      double precision el
      integer k
      double precision pi
      double precision r
      double precision rd
      double precision s
      double precision s0
      double precision ti
      double precision tl0
      double precision x

      pi = 3.141592653589793D+00
      r = 1.0D+00

      if ( x .le. 20.0D+00 ) then

        s = 0.5D+00
        do k = 1, 100
 
          if ( k .eq. 1 ) then
            rd = 0.5D+00
          else
            rd = 1.0D+00
          end if
          r = r * rd * k / ( k + 1.0D+00 ) 
     &      * ( x / ( 2.0D+00 * k + 1.0D+00 ) ) ** 2
          s = s + r
          if ( abs ( r / s ) .lt. 1.0D-12 ) then
            go to 10
          end if
        end do

10      continue

        tl0 = 2.0D+00 / pi * x * x * s

      else

        s = 1.0D+00
        do k = 1, 10
          r = r * k / ( k + 1.0D+00 ) 
     &      * ( ( 2.0D+00 * k + 1.0D+00 ) / x ) ** 2
          s = s + r
          if ( abs ( r / s ) .lt. 1.0D-12 ) then
            go to 20
          end if
        end do

20      continue

        el = 0.57721566490153D+00
        s0 = - s / ( pi * x * x ) + 2.0D+00 / pi 
     &    * ( log ( 2.0D+00 * x ) + el )
        a0 = 1.0D+00
        a1 = 5.0D+00 / 8.0D+00
        a(1) = a1
        do k = 1, 10
          af = ( ( 1.5D+00 * ( k + 0.50D+00 ) 
     &      * ( k + 5.0D+00 / 6.0D+00 ) * a1 - 0.5D+00 
     &      * ( k + 0.5D+00 ) ** 2 * ( k -0.5D+00 ) * a0 ) ) 
     &      / ( k + 1.0D+00 )
          a(k+1) = af
          a0 = a1
          a1 = af
        end do

        ti = 1.0D+00
        r = 1.0D+00
        do k = 1, 11
          r = r / x
          ti = ti + a(k) * r
        end do
        tl0 = ti / sqrt ( 2.0D+00 * pi * x ) * exp ( x ) + s0

      end if

      return
      end
      subroutine itth0 ( x, tth )

c*********************************************************************72
c
cc ITTH0 integrates H0(t)/t from x to oo.
c
c  Licensing:
c
c    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
c    they give permission to incorporate this routine into a user program 
c    provided that the copyright is acknowledged.
c
c  Modified:
c
c    23 July 2012
c
c  Author:
c
c    Shanjie Zhang, Jianming Jin
c
c  Reference:
c
c    Shanjie Zhang, Jianming Jin,
c    Computation of Special Functions,
c    Wiley, 1996,
c    ISBN: 0-471-11963-6,
c    LC: QA351.C45.
c
c  Parameters:
c
c    Input, double precision X, the lower limit of the integral.
c
c    Output, double precision TTH, the integral of H0(t)/t from x to oo.
c
      implicit none

      double precision f0
      double precision g0
      integer k
      double precision pi
      double precision r
      double precision s
      double precision t
      double precision tth
      double precision tty
      double precision x
      double precision xt

      pi = 3.141592653589793D+00
      s = 1.0D+00
      r = 1.0D+00

      if ( x .lt. 24.5D+00 ) then

        do k = 1, 60
          r = - r * x * x * ( 2.0D+00 * k - 1.0D+00 )
     &      / ( 2.0D+00 * k + 1.0D+00 ) ** 3
          s = s + r
          if ( abs ( r ) .lt. abs ( s ) * 1.0D-12 ) then
            go to 10
          end if
        end do

10      continue

        tth = pi / 2.0D+00 - 2.0D+00 / pi * x * s

      else

        do k = 1, 10
          r = - r * ( 2.0D+00 * k - 1.0D+00 ) ** 3 
     &      / ( ( 2.0D+00 * k + 1.0D+00 ) * x * x )
          s = s + r
          if ( abs ( r ) .lt. abs ( s ) * 1.0D-12 ) then
            go to 20
          end if
        end do

20      continue

        tth = 2.0D+00 / ( pi * x ) * s
        t = 8.0D+00 / x
        xt = x + 0.25D+00 * pi
        f0 = (((((
     &      0.18118D-02 * t
     &    - 0.91909D-02 ) * t
     &    + 0.017033D+00 ) * t
     &    - 0.9394D-03 ) * t
     &    - 0.051445D+00 ) * t
     &    - 0.11D-05 ) * t
     &    + 0.7978846D+00
        g0 = (((((
     &    - 0.23731D-02 * t
     &    + 0.59842D-02 ) * t
     &    + 0.24437D-02 ) * t
     &    - 0.0233178D+00 ) * t
     &    + 0.595D-04 ) * t
     &    + 0.1620695D+00 ) * t
        tty = ( f0 * sin ( xt ) - g0 * cos ( xt ) ) / ( sqrt ( x ) * x )
        tth = tth + tty

        end if

      return
      end
      subroutine ittika ( x, tti, ttk )

c*********************************************************************72
c
cc ITTIKA integrates (I0(t)-1)/t from 0 to x, K0(t)/t from x to infinity.
c
c  Licensing:
c
c    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
c    they give permission to incorporate this routine into a user program 
c    provided that the copyright is acknowledged.
c
c  Modified:
c
c    23 July 2012
c
c  Author:
c
c    Shanjie Zhang, Jianming Jin
c
c  Reference:
c
c    Shanjie Zhang, Jianming Jin,
c    Computation of Special Functions,
c    Wiley, 1996,
c    ISBN: 0-471-11963-6,
c    LC: QA351.C45.
c
c  Parameters:
c
c    Input, double precision X, the integral limit.
c
c    Output, double precision TTI, TTK, the integrals of [I0(t)-1]/t from 0 to x,
c    and of K0(t)/t from x to oo.
c
      implicit none

      double precision b1
      double precision c(8)
      double precision e0
      double precision el
      integer k
      double precision pi
      double precision r
      double precision r2
      double precision rc
      double precision rs
      double precision tti
      double precision ttk
      double precision x

      save c

      data c / 1.625D+00, 4.1328125D+00,
     &       1.45380859375D+01, 6.553353881835D+01,
     &       3.6066157150269D+02, 2.3448727161884D+03,
     &       1.7588273098916D+04, 1.4950639538279D+05 /

      pi = 3.141592653589793D+00
      el = 0.5772156649015329D+00

      if ( x .eq. 0.0D+00 ) then
        tti = 0.0D+00
        ttk = 1.0D+300
        return
      end if

      if ( x .lt. 40.0D+00 ) then
        tti = 1.0D+00
        r = 1.0D+00
        do k = 2, 50
          r = 0.25D+00 * r * ( k - 1.0D+00 ) / ( k * k * k ) * x * x
          tti = tti + r
          if ( abs ( r / tti ) .lt. 1.0D-12 ) then
            go to 10
          end if
        end do

10      continue

        tti = tti * 0.125D+00 * x * x

      else

        tti = 1.0D+00
        r = 1.0D+00
        do k = 1, 8
          r = r / x
          tti = tti + c(k) * r
        end do
        rc = x * sqrt ( 2.0D+00 * pi * x )
        tti = tti * exp ( x ) / rc

      end if

      if ( x .le. 12.0D+00 ) then

        e0 = ( 0.5D+00 * log ( x / 2.0D+00 ) + el ) 
     &    * log ( x / 2.0D+00 ) + pi * pi / 24.0D+00 + 0.5D+00 * el * el
        b1 = 1.5D+00 - ( el + log ( x / 2.0D+00 ) )
        rs = 1.0D+00
        r = 1.0D+00
        do k = 2, 50
          r = 0.25D+00 * r * ( k - 1.0D+00 ) / ( k * k * k ) * x * x
          rs = rs + 1.0D+00 / k
          r2 = r * ( rs + 1.0D+00 / ( 2.0D+00 * k ) 
     &      - ( el + log ( x / 2.0D+00 ) ) )
          b1 = b1 + r2
          if ( abs ( r2 / b1 ) .lt. 1.0D-12 ) then
            go to 20
          end if
        end do

20      continue

        ttk = e0 - 0.125D+00 * x * x * b1

      else

        ttk = 1.0D+00
        r = 1.0D+00
        do k = 1, 8
          r = - r / x
          ttk = ttk + c(k) * r
        end do
        rc = x * sqrt ( 2.0D+00 / pi * x )
        ttk = ttk * exp ( - x ) / rc

      end if

      return
      end
      subroutine ittikb ( x, tti, ttk )

c*********************************************************************72
c
cc ITTIKB integrates (I0(t)-1)/t from 0 to x, K0(t)/t from x to infinity.
c
c  Licensing:
c
c    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
c    they give permission to incorporate this routine into a user program 
c    provided that the copyright is acknowledged.
c
c  Modified:
c
c    28 July 2012
c
c  Author:
c
c    Shanjie Zhang, Jianming Jin
c
c  Reference:
c
c    Shanjie Zhang, Jianming Jin,
c    Computation of Special Functions,
c    Wiley, 1996,
c    ISBN: 0-471-11963-6,
c    LC: QA351.C45.
c
c  Parameters:
c
c    Input, double precision X, the integral limit.
c
c    Output, double precision TTI, TTK, the integrals of
c    [I0(t)-1]/t from 0 to x, and K0(t)/t from x to oo.
c
      implicit none

      double precision e0
      double precision el
      double precision pi
      double precision t
      double precision t1
      double precision tti
      double precision ttk
      double precision x
      double precision x1

      pi = 3.141592653589793D+00
      el = 0.5772156649015329D+00

      if ( x .eq. 0.0D+00 ) then

        tti = 0.0D+00

      else if ( x .le. 5.0D+00 ) then

        x1 = x / 5.0D+00
        t = x1 * x1
        tti = (((((((
     &      0.1263D-03       * t
     &    + 0.96442D-03 )    * t
     &    + 0.968217D-02 )   * t
     &    + 0.06615507D+00 ) * t
     &    + 0.33116853D+00 ) * t
     &    + 1.13027241D+00 ) * t
     &    + 2.44140746D+00 ) * t
     &    + 3.12499991D+00 ) * t

      else

        t = 5.0D+00 / x
        tti = (((((((((
     &       2.1945464D+00   * t
     &    -  3.5195009D+00 ) * t
     &    - 11.9094395D+00 ) * t
     &    + 40.394734D+00  ) * t
     &    - 48.0524115D+00 ) * t
     &    + 28.1221478D+00 ) * t
     &    -  8.6556013D+00 ) * t
     &    +  1.4780044D+00 ) * t
     &    -  0.0493843D+00 ) * t
     &    +  0.1332055D+00 ) * t
     &    +  0.3989314D+00
        tti = tti * exp ( x ) / ( sqrt ( x ) * x )

      end if

      if ( x .eq. 0.0D+00 ) then

        ttk = 1.0D+300

      else if ( x .le. 2.0D+00 ) then

        t1 = x / 2.0D+00
        t = t1 * t1
        ttk = (((((
     &      0.77D-06         * t
     &    + 0.1544D-04 )     * t
     &    + 0.48077D-03 )    * t
     &    + 0.925821D-02 )   * t
     &    + 0.10937537D+00 ) * t
     &    + 0.74999993D+00 ) * t
        e0 = el + log ( x / 2.0D+00 )
        ttk = pi * pi / 24.0D+00 + e0 * ( 0.5D+00 * e0 + tti ) - ttk

      else if ( x .le. 4.0D+00 ) then

        t = 2.0D+00 / x
        ttk = (((
     &      0.06084D+00    * t
     &    - 0.280367D+00 ) * t
     &    + 0.590944D+00 ) * t
     &    - 0.850013D+00 ) * t
     &    + 1.234684D+00
        ttk = ttk * exp ( - x ) / ( sqrt ( x ) * x )

      else

        t = 4.0D+00 / x
        ttk = (((((
     &      0.02724D+00     * t
     &    - 0.1110396D+00 ) * t
     &    + 0.2060126D+00 ) * t
     &    - 0.2621446D+00 ) * t
     &    + 0.3219184D+00 ) * t
     &    - 0.5091339D+00 ) * t
     &    + 1.2533141D+00
        ttk = ttk * exp ( - x ) / ( sqrt ( x ) * x )

      end if

      return
      end
      subroutine ittjya ( x, ttj, tty )

c*********************************************************************72
c
cc ITTJYA integrates (1-J0(t))/t from 0 to x, and Y0(t)/t from x to infinity.
c
c  Licensing:
c
c    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
c    they give permission to incorporate this routine into a user program 
c    provided that the copyright is acknowledged.
c
c  Modified:
c
c    28 July 2012
c
c  Author:
c
c    Shanjie Zhang, Jianming Jin
c
c  Reference:
c
c    Shanjie Zhang, Jianming Jin,
c    Computation of Special Functions,
c    Wiley, 1996,
c    ISBN: 0-471-11963-6,
c    LC: QA351.C45.
c
c  Parameters:
c
c    Input, double precision X, the integral limit.
c
c    Output, double precision TTJ, TTY, the integrals of [1-J0(t)]/t 
c    from 0 to x and of Y0(t)/t from x to oo.
c
      implicit none

      double precision a0
      double precision b1
      double precision bj0
      double precision bj1
      double precision by0
      double precision by1
      double precision e0
      double precision el
      double precision g0
      double precision g1
      integer k
      integer l
      double precision pi
      double precision px
      double precision qx
      double precision r
      double precision r0
      double precision r1
      double precision r2
      double precision rs
      double precision t
      double precision ttj
      double precision tty
      double precision vt
      double precision x
      double precision xk

      pi = 3.141592653589793D+00
      el = 0.5772156649015329D+00

      if ( x .eq. 0.0D+00 ) then

        ttj = 0.0D+00
        tty = -1.0D+300

      else if ( x .le. 20.0D+00 ) then

        ttj = 1.0D+00
        r = 1.0D+00
        do k = 2, 100
          r = - 0.25D+00 * r * ( k - 1.0D+00 ) / ( k * k * k ) * x * x
          ttj = ttj + r
          if ( abs ( r ) .lt. abs ( ttj ) * 1.0D-12 ) then
            go to 10
          end if
        end do

10      continue

        ttj = ttj * 0.125D+00 * x * x
        e0 = 0.5D+00 * ( pi * pi / 6.0D+00 - el * el ) 
     &    - ( 0.5D+00 * log ( x / 2.0D+00 ) + el )
     &    * log ( x / 2.0D+00 )
        b1 = el + log ( x / 2.0D+00 ) - 1.5D+00
        rs = 1.0D+00
        r = -1.0D+00
        do k = 2, 100
          r = - 0.25D+00 * r * ( k - 1.0D+00 ) / ( k * k * k ) * x * x
          rs = rs + 1.0D+00 / k
          r2 = r * ( rs + 1.0D+00 / ( 2.0D+00 * k ) 
     &      - ( el + log ( x / 2.0D+00 ) ) ) 
          b1 = b1 + r2
          if ( abs ( r2 ) .lt. abs ( b1 ) * 1.0D-12 ) then
            go to 20
          end if
        end do

20      continue

        tty = 2.0D+00 / pi * ( e0 + 0.125D+00 * x * x * b1 )

      else

        a0 = sqrt ( 2.0D+00 / ( pi * x ) )

        do l = 0, 1

          vt = 4.0D+00 * l * l
          px = 1.0D+00
          r = 1.0D+00
          do k = 1, 14
            r = - 0.0078125D+00 * r 
     &        * ( vt - ( 4.0D+00 * k - 3.0D+00 ) ** 2 )
     &        / ( x * k ) * ( vt - ( 4.0D+00 * k - 1.0D+00 ) ** 2 )
     &        / ( ( 2.0D+00 * k - 1.0D+00 ) * x )
            px = px + r
            if ( abs ( r ) .lt. abs ( px ) * 1.0D-12 ) then
              go to 30
            end if
          end do

30        continue

          qx = 1.0D+00
          r = 1.0D+00
          do k = 1, 14
            r = -0.0078125D+00 * r 
     &        * ( vt - ( 4.0D+00 * k - 1.0D+00 ) ** 2 )
     &        / ( x * k ) * ( vt - ( 4.0D+00 * k + 1.0D+00 ) ** 2 )
     &        / ( 2.0D+00 * k + 1.0D+00 ) / x
            qx = qx + r
            if ( abs ( r ) .lt. abs ( qx ) * 1.0D-12 ) then
              go to 40
            end if
          end do

40        continue

          qx = 0.125D+00 * ( vt - 1.0D+00 ) / x * qx
          xk = x - ( 0.25D+00 + 0.5D+00 * l ) * pi
          bj1 = a0 * ( px * cos ( xk ) - qx * sin ( xk ) )
          by1 = a0 * ( px * sin ( xk ) + qx * cos ( xk ) )
          if ( l .eq. 0 ) then
            bj0 = bj1
            by0 = by1
          end if

        end do

        t = 2.0D+00 / x
        g0 = 1.0D+00
        r0 = 1.0D+00
        do k = 1, 10
          r0 = - k * k * t * t *r0
          g0 = g0 + r0
        end do

        g1 = 1.0D+00
        r1 = 1.0D+00
        do k = 1, 10
          r1 = - k * ( k + 1.0D+00 ) * t * t * r1
          g1 = g1 + r1
        end do

        ttj = 2.0D+00 * g1 * bj0 / ( x * x ) - g0 * bj1 / x 
     &    + el + log ( x / 2.0D+00 )
        tty = 2.0D+00 * g1 * by0 / ( x * x ) - g0 * by1 / x

      end if

      return
      end
      subroutine ittjyb ( x, ttj, tty )

c*********************************************************************72
c
cc ITTJYB integrates (1-J0(t))/t from 0 to x, and Y0(t)/t from x to infinity.
c
c  Licensing:
c
c    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
c    they give permission to incorporate this routine into a user program 
c    provided that the copyright is acknowledged.
c
c  Modified:
c
c    01 August 2012
c
c  Author:
c
c    Shanjie Zhang, Jianming Jin
c
c  Reference:
c
c    Shanjie Zhang, Jianming Jin,
c    Computation of Special Functions,
c    Wiley, 1996,
c    ISBN: 0-471-11963-6,
c    LC: QA351.C45.
c
c  Parameters:
c
c    Input, double precision X, the integral limit.
c
c    Output, double precision TTJ, TTY, the integrals of [1-J0(t)]/t 
c    from 0 to x and of Y0(t)/t from x to oo.
c
      implicit none

      double precision e0
      double precision el
      double precision f0
      double precision g0
      double precision pi
      double precision t
      double precision t1
      double precision ttj
      double precision tty
      double precision x
      double precision x1
      double precision xt

      pi = 3.141592653589793D+00
      el = 0.5772156649015329D+00

      if ( x .eq. 0.0D+00 ) then

        ttj = 0.0D+00
        tty = -1.0D+300

      else if ( x .le. 4.0D+00 ) then

        x1 = x / 4.0D+00
        t = x1 * x1

        ttj = ((((((
     &      0.35817D-04 * t
     &    - 0.639765D-03 ) * t
     &    + 0.7092535D-02 ) * t
     &    - 0.055544803D+00 ) * t
     &    + 0.296292677D+00 ) * t
     &    - 0.999999326D+00 ) * t
     &    + 1.999999936D+00 ) * t

        tty = (((((((
     &    - 0.3546D-05        * t
     &    + 0.76217D-04 )     * t
     &    - 0.1059499D-02 )   * t
     &    + 0.010787555D+00 ) * t
     &    - 0.07810271D+00 )  * t
     &    + 0.377255736D+00 ) * t
     &    - 1.114084491D+00 ) * t
     &    + 1.909859297D+00 ) * t

        e0 = el + log ( x / 2.0D+00 )
        tty = pi / 6.0D+00 + e0 / pi * ( 2.0D+00 * ttj - e0 ) - tty

      else if ( x .le. 8.0D+00 ) then

        xt = x + 0.25D+00 * pi
        t1 = 4.0D+00 / x
        t = t1 * t1

        f0 = (((((
     &      0.0145369D+00 * t
     &    - 0.0666297D+00 ) * t
     &    + 0.1341551D+00 ) * t
     &    - 0.1647797D+00 ) * t
     &    + 0.1608874D+00 ) * t
     &    - 0.2021547D+00 ) * t
     &    + 0.7977506D+00

        g0 = ((((((
     &      0.0160672D+00   * t
     &    - 0.0759339D+00 ) * t
     &    + 0.1576116D+00 ) * t
     &    - 0.1960154D+00 ) * t
     &    + 0.1797457D+00 ) * t
     &    - 0.1702778D+00 ) * t
     &    + 0.3235819D+00 ) * t1

        ttj = ( f0 * cos ( xt ) + g0 * sin ( xt ) ) / ( sqrt ( x ) * x )
        ttj = ttj + el + log ( x / 2.0D+00 )
        tty = ( f0 * sin ( xt ) - g0 * cos ( xt ) ) / ( sqrt ( x ) * x )

      else

        t = 8.0D+00 / x
        xt = x + 0.25D+00 * pi

        f0 = (((((
     &      0.18118D-02    * t
     &    - 0.91909D-02 )  * t
     &    + 0.017033D+00 ) * t
     &    - 0.9394D-03 )   * t
     &    - 0.051445D+00 ) * t
     &    - 0.11D-05 )     * t
     &    + 0.7978846D+00

        g0 = (((((
     &    - 0.23731D-02     * t
     &    + 0.59842D-02 )   * t
     &    + 0.24437D-02 )   * t
     &    - 0.0233178D+00 ) * t
     &    + 0.595D-04 )     * t
     &    + 0.1620695D+00 ) * t

        ttj = ( f0 * cos ( xt ) + g0 * sin ( xt ) ) 
     &    / ( sqrt ( x ) * x ) + el + log ( x / 2.0D+00 )
        tty = ( f0 * sin ( xt ) - g0 * cos ( xt ) ) 
     &    / ( sqrt ( x ) * x )

      end if

      return
      end
      subroutine jdzo ( nt, n, m, p, zo )

c*********************************************************************72
c
cc JDZO computes the zeros of Bessel functions Jn(x) and Jn'(x).
c
c  Discussion:
c
c    This procedure computes the zeros of Bessel functions Jn(x) and
c    Jn'(x), and arrange them in the order of their magnitudes.
c
c  Licensing:
c
c    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
c    they give permission to incorporate this routine into a user program 
c    provided that the copyright is acknowledged.
c
c  Modified:
c
c    01 August 2012
c
c  Author:
c
c    Shanjie Zhang, Jianming Jin
c
c  Reference:
c
c    Shanjie Zhang, Jianming Jin,
c    Computation of Special Functions,
c    Wiley, 1996,
c    ISBN: 0-471-11963-6,
c    LC: QA351.C45.
c
c  Parameters:
c
c    Input, integer NT, the number of zeros.
c
c    Output, integer N(*), the  order of Jn(x) or Jn'(x) associated
c    with the L-th zero.
c
c    Output, integer M(*), the serial number of the zeros of Jn(x)
c    or Jn'(x) associated with the L-th zero ( L is the serial number of all the
c    zeros of Jn(x) and Jn'(x) ).
c
c    Output, character*4 P(L), 'TM' or 'TE', a code for designating the
c    zeros of Jn(x)  or Jn'(x).  In the waveguide applications, the zeros
c    of Jn(x) correspond to TM modes and those of Jn'(x) correspond to TE modes.
c
c    Output, double precision ZO(*), the zeros of Jn(x) and Jn'(x).
c
      implicit none

      double precision bj(101)
      double precision dj(101)
      double precision fj(101)
      integer i
      integer j
      integer k
      integer l
      integer l0
      integer l1
      integer l2
      integer m(1400)
      integer m1(70)
      integer mm
      integer n(1400)
      integer n1(70)
      integer nm
      integer nt
      character p(1400)*4
      character p1(70)*4
      double precision x
      double precision x0
      double precision x1
      double precision x2
      double precision xm
      double precision zo(1400)
      double precision zoc(70)

      if ( nt .lt. 600 ) then
        xm = -1.0D+00 + 2.248485D+00 * dble ( nt ) ** 0.5D+00 
     &    - 0.0159382D+00 * nt + 3.208775D-04 * dble ( nt ) ** 1.5D+00
        nm = int ( 14.5D+00 + 0.05875D+00 * nt )
        mm = int ( 0.02D+00 * nt ) + 6
      else
        xm = 5.0D+00 + 1.445389D+00 * ( dble ( nt ) ) ** 0.5D+00
     &    + 0.01889876D+00 * nt 
     &    - 2.147763D-04 * ( dble ( nt ) ) ** 1.5D+00
        nm = int ( 27.8D+00 + 0.0327D+00 * nt )
        mm = int ( 0.01088D+00 * nt ) + 10
      end if

      l0 = 0

      do i = 1,nm

        x1 = 0.407658D+00 + 0.4795504D+00 
     &    * ( dble ( i - 1 ) ) ** 0.5D+00 + 0.983618D+00 * ( i - 1 )
        x2 = 1.99535D+00 + 0.8333883 * ( dble ( i - 1 ) ) ** 0.5D+00 
     &    + 0.984584D+00 * ( i - 1 )
        l1 = 0

        do j = 1, mm

          if ( i .eq. 1 .and. j .eq. 1 ) then
            go to 20
          end if

          x = x1

10        continue

          call bjndd ( i, x, bj, dj, fj )
          x0 = x
          x = x - dj(i) / fj(i)

          if ( xm .lt. x1 ) then
            go to 30
          end if

          if ( 1.0D-10 .lt. abs ( x - x0 ) ) then
            go to 10
          end if

20        continue

          l1 = l1 + 1
          n1(l1) = i - 1
          m1(l1) = j
          if ( i .eq. 1 ) then
            m1(l1) = j - 1
          end if
          p1(l1) = 'TE'
          zoc(l1) = x

          if ( i .le. 15 ) then
            x1 = x + 3.057D+00 + 0.0122D+00 * ( i - 1 )
     &        + ( 1.555D+00 + 0.41575D+00 * ( i - 1 ) ) / ( j + 1 ) ** 2
          else
            x1 = x + 2.918D+00 + 0.01924D+00 * ( i - 1 )
     &        + ( 6.26D+00 + 0.13205D+00 * ( i - 1 ) ) / ( j + 1 ) ** 2
          end if

30        continue

          x = x2

40        continue

          call bjndd ( i, x, bj, dj, fj )
          x0 = x
          x = x - bj(i) / dj(i)

          if ( x .le. xm ) then

            if ( 1.0D-10 .lt. abs ( x - x0 ) ) then
              go to 40
            end if
            l1 = l1 + 1
            n1(l1) = i - 1
            m1(l1) = j
            p1(l1) = 'TM'
            zoc(l1) = x
            if ( i .le. 15 ) then
              x2 = x + 3.11D+00 + 0.0138D+00 * ( i - 1 )
     &          + ( 0.04832D+00 + 0.2804D+00 * ( i - 1 ) ) 
     &          / ( j + 1 ) ** 2
            else
              x2 = x + 3.001D+00 + 0.0105D+00 * ( i - 1 )
     &          + ( 11.52D+00 + 0.48525D+00 * ( i - 1 ) ) 
     &          / ( j + 3 ) ** 2
            end if

          end if

        end do

        l = l0 + l1
        l2 = l

50      continue

        if ( l0 .eq. 0 ) then
          do k = 1, l
            zo(k) = zoc(k)
            n(k) = n1(k)
            m(k) = m1(k)
            p(k) = p1(k)
          end do
          l1 = 0
        else if ( l0 .ne. 0 ) then
          if ( zoc(l1) .le. zo(l0) ) then
            zo(l0+l1) = zo(l0)
            n(l0+l1) = n(l0)
            m(l0+l1) = m(l0)
            p(l0+l1) = p(l0)
            l0 = l0 - 1
          else
            zo(l0+l1) = zoc(l1)
            n(l0+l1) = n1(l1)
            m(l0+l1) = m1(l1)
            p(l0+l1) = p1(l1)
            l1 = l1 - 1
          end if
        end if

        if ( l1 .ne. 0 ) then
          go to 50
        end if

        l0 = l2

      end do

      return
      end
      subroutine jelp ( u, hk, esn, ecn, edn, eph )

c*********************************************************************72
c
cc JELP computes Jacobian elliptic functions SN(u), CN(u), DN(u).
c
c  Licensing:
c
c    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
c    they give permission to incorporate this routine into a user program 
c    provided that the copyright is acknowledged.
c
c  Modified:
c
c    08 July 2012
c
c  Author:
c
c    Shanjie Zhang, Jianming Jin
c
c  Reference:
c
c    Shanjie Zhang, Jianming Jin,
c    Computation of Special Functions,
c    Wiley, 1996,
c    ISBN: 0-471-11963-6,
c    LC: QA351.C45.
c
c  Parameters:
c
c    Input, double precision U, the argument.
c
c    Input, double precision HK, the modulus, between 0 and 1.
c
c    Output, double precision ESN, ECN, EDN, EPH, the values of
c    sn(u), cn(u), dn(u), and phi (in degrees).
c
      implicit none

      double precision a
      double precision a0
      double precision b
      double precision b0
      double precision c
      double precision d
      double precision dn
      double precision ecn
      double precision edn
      double precision eph
      double precision esn
      double precision hk
      integer j
      integer n 
      double precision pi
      double precision r(40)
      double precision sa
      double precision t
      double precision u

      pi = 3.14159265358979D+00
      a0 = 1.0D+00
      b0 = sqrt ( 1.0D+00 - hk * hk )

      do n = 1, 40

        a = ( a0 + b0 ) / 2.0D+00
        b = sqrt ( a0 * b0 )
        c = ( a0 - b0 ) / 2.0D+00
        r(n) = c / a

        if ( c .lt. 1.0D-07 ) then
          go to 10
        end if

        a0 = a
        b0 = b

      end do

10    continue

      dn = 2.0D+00 ** n * a * u

      do j = n, 1, -1
        t = r(j) * sin ( dn )
        sa = atan ( t / sqrt ( abs ( 1.0D+00 - t * t )))
        d = 0.5D+00 * ( dn + sa )
        dn = d
      end do

      eph = d * 180.0D+00 / pi
      esn = sin ( d )
      ecn = cos ( d )
      edn = sqrt ( 1.0D+00 - hk * hk * esn * esn )

      return
      end
      subroutine jy01a ( x, bj0, dj0, bj1, dj1, by0, dy0, by1, dy1 )

c*********************************************************************72
c
cc JY01A computes Bessel functions J0(x), J1(x), Y0(x), Y1(x) and derivatives.
c
c  Licensing:
c
c    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
c    they give permission to incorporate this routine into a user program 
c    provided that the copyright is acknowledged.
c
c  Modified:
c
c    01 August 2012
c
c  Author:
c
c    Shanjie Zhang, Jianming Jin
c
c  Reference:
c
c    Shanjie Zhang, Jianming Jin,
c    Computation of Special Functions,
c    Wiley, 1996,
c    ISBN: 0-471-11963-6,
c    LC: QA351.C45.
c
c  Parameters:
c
c    Input, double precision X, the argument.
c
c    Output, double precision BJ0, DJ0, BJ1, DJ1, BY0, DY0, BY1, DY1,
c    the values of J0(x), J0'(x), J1(x), J1'(x), Y0(x), Y0'(x), Y1(x), Y1'(x).
c
      implicit none

      double precision a(12)
      double precision b(12)
      double precision a1(12)
      double precision b1(12)
      double precision bj0
      double precision bj1
      double precision by0
      double precision by1
      double precision cs0
      double precision cs1
      double precision cu
      double precision dj0
      double precision dj1
      double precision dy0
      double precision dy1
      double precision ec
      integer k
      integer k0
      double precision p0
      double precision p1
      double precision pi
      double precision q0
      double precision q1
      double precision r
      double precision r0
      double precision r1
      double precision rp2
      double precision t1
      double precision t2
      double precision w0
      double precision w1
      double precision x
      double precision x2

      save a
      save a1
      save b
      save b1

      data a     /-0.7031250000000000D-01, 0.1121520996093750D+00,
     &            -0.5725014209747314D+00, 0.6074042001273483D+01,
     &            -0.1100171402692467D+03, 0.3038090510922384D+04,
     &            -0.1188384262567832D+06, 0.6252951493434797D+07,
     &            -0.4259392165047669D+09, 0.3646840080706556D+11,
     &            -0.3833534661393944D+13, 0.4854014686852901D+15/

      data a1     /0.1171875000000000D+00, -0.1441955566406250D+00,
     &             0.6765925884246826D+00, -0.6883914268109947D+01,
     &             0.1215978918765359D+03, -0.3302272294480852D+04,
     &             0.1276412726461746D+06, -0.6656367718817688D+07,
     &             0.4502786003050393D+09, -0.3833857520742790D+11,
     &             0.4011838599133198D+13, -0.5060568503314727D+15/

      data b     / 0.7324218750000000D-01, -0.2271080017089844D+00,
     &             0.1727727502584457D+01, -0.2438052969955606D+02,
     &             0.5513358961220206D+03, -0.1825775547429318D+05,
     &             0.8328593040162893D+06, -0.5006958953198893D+08,
     &             0.3836255180230433D+10, -0.3649010818849833D+12,
     &             0.4218971570284096D+14, -0.5827244631566907D+16/

      data b1     /-0.1025390625000000D+00, 0.2775764465332031D+00,
     &             -0.1993531733751297D+01, 0.2724882731126854D+02,
     &             -0.6038440767050702D+03, 0.1971837591223663D+05,
     &             -0.8902978767070678D+06, 0.5310411010968522D+08,
     &             -0.4043620325107754D+10, 0.3827011346598605D+12,
     &             -0.4406481417852278D+14, 0.6065091351222699D+16/

      pi = 3.141592653589793D+00
      rp2 = 0.63661977236758D+00
      x2 = x * x

      if ( x .eq. 0.0D+00 ) then
        bj0 = 1.0D+00
        bj1 = 0.0D+00
        dj0 = 0.0D+00
        dj1 = 0.5D+00
        by0 = -1.0D+300
        by1 = -1.0D+300
        dy0 = 1.0D+300
        dy1 = 1.0D+300
        return
      end if

      if ( x .le. 12.0D+00 ) then

        bj0 = 1.0D+00
        r = 1.0D+00
        do k = 1,30
          r = -0.25D+00 * r * x2 / ( k * k )
          bj0 = bj0 + r
          if ( abs ( r ) .lt. abs ( bj0 ) * 1.0D-15 ) then
            go to 10
          end if
        end do

10      continue

        bj1 = 1.0D+00
        r = 1.0D+00
        do k = 1, 30
          r = -0.25D+00 * r * x2 / ( k * ( k + 1.0D+00 ) )
          bj1 = bj1 + r
          if ( abs ( r ) .lt. abs ( bj1 ) * 1.0D-15 ) then
            go to 20
          end if
        end do

20      continue

        bj1 = 0.5D+00 * x * bj1
        ec = log ( x / 2.0D+00 ) + 0.5772156649015329D+00
        cs0 = 0.0D+00
        w0 = 0.0D+00
        r0 = 1.0D+00
        do k = 1, 30
          w0 = w0 + 1.0D+00 / k
          r0 = -0.25D+00 * r0 / ( k * k ) * x2
          r = r0 * w0
          cs0 = cs0 + r
          if ( abs ( r ) .lt. abs ( cs0 ) * 1.0D-15 ) then
            go to 30
          end if
        end do

30      continue

        by0 = rp2 * ( ec * bj0 - cs0 )
        cs1 = 1.0D+00
        w1 = 0.0D+00
        r1 = 1.0D+00
        do k = 1, 30
          w1 = w1 + 1.0D+00 / k
          r1 = -0.25D+00 * r1 / ( k * ( k + 1 ) ) * x2
          r = r1 * ( 2.0D+00 * w1 + 1.0D+00 / ( k + 1.0D+00 ) )
          cs1 = cs1 + r
          if ( abs ( r ) .lt. abs ( cs1 ) * 1.0D-15 ) then
            go to 40
          end if
        end do

40      continue

        by1 = rp2 * ( ec * bj1 - 1.0D+00 / x - 0.25D+00 * x * cs1 )

      else

        if ( x .lt. 35.0D+00 ) then
          k0 = 12
        else if ( x .lt. 50.0D+00 ) then
          k0 = 10
        else
          k0 = 8
        end if

        t1 = x - 0.25D+00 * pi
        p0 = 1.0D+00
        q0 = -0.125D+00 / x
        do k = 1, k0
          p0 = p0 + a(k) * x ** ( - 2 * k )
          q0 = q0 + b(k) * x ** ( - 2 * k - 1 )
        end do
        cu = dsqrt ( rp2 / x )
        bj0 = cu * ( p0 * cos ( t1 ) - q0 * sin ( t1 ) )
        by0 = cu * ( p0 * sin ( t1 ) + q0 * cos ( t1 ) )
        t2 = x - 0.75D+00 * pi
        p1 = 1.0D+00
        q1 = 0.375D+00 / x
        do k = 1, k0
          p1 = p1 + a1(k) * x ** ( - 2 * k )
          q1 = q1 + b1(k) * x ** ( - 2 * k - 1 )
        end do
        cu = dsqrt ( rp2 / x )
        bj1 = cu * ( p1 * cos ( t2 ) - q1 * sin ( t2 ) )
        by1 = cu * ( p1 * sin ( t2 ) + q1 * cos ( t2 ) )

      end if

      dj0 = - bj1
      dj1 = bj0 - bj1 / x
      dy0 = - by1
      dy1 = by0 - by1 / x

      return
      end
      subroutine jy01b ( x, bj0, dj0, bj1, dj1, by0, dy0, by1, dy1 )

c*********************************************************************72
c
cc JY01B computes Bessel functions J0(x), J1(x), Y0(x), Y1(x) and derivatives.
c
c  Licensing:
c
c    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
c    they give permission to incorporate this routine into a user program 
c    provided that the copyright is acknowledged.
c
c  Modified:
c
c    02 August 2012
c
c  Author:
c
c    Shanjie Zhang, Jianming Jin
c
c  Reference:
c
c    Shanjie Zhang, Jianming Jin,
c    Computation of Special Functions,
c    Wiley, 1996,
c    ISBN: 0-471-11963-6,
c    LC: QA351.C45.
c
c  Parameters:
c
c    Input, double precision X, the argument.
c
c    Output, double precision BJ0, DJ0, BJ1, DJ1, BY0, DY0, BY1, DY1,
c    the values of J0(x), J0'(x), J1(x), J1'(x), Y0(x), Y0'(x), Y1(x), Y1'(x).
c
      implicit none

      double precision a0
      double precision bj0
      double precision bj1
      double precision by0
      double precision by1
      double precision dj0
      double precision dj1
      double precision dy0
      double precision dy1
      double precision p0
      double precision p1
      double precision pi
      double precision q0
      double precision q1
      double precision t
      double precision t2
      double precision ta0
      double precision ta1
      double precision x

      pi = 3.141592653589793D+00

      if ( x .eq. 0.0D+00 ) then

        bj0 = 1.0D+00
        bj1 = 0.0D+00
        dj0 = 0.0D+00
        dj1 = 0.5D+00
        by0 = -1.0D+300
        by1 = -1.0D+300
        dy0 = 1.0D+300
        dy1 = 1.0D+300
        return

      else if ( x .le. 4.0D+00 ) then

        t = x / 4.0D+00
        t2 = t * t

        bj0 = ((((((
     &    - 0.5014415D-03 * t2
     &    + 0.76771853D-02 ) * t2
     &    - 0.0709253492D+00 ) * t2
     &    + 0.4443584263D+00 ) * t2
     &    - 1.7777560599D+00 ) * t2
     &    + 3.9999973021D+00 ) * t2
     &    - 3.9999998721D+00 ) * t2
     &    + 1.0D+00

        bj1 = t * (((((((
     &    - 0.1289769D-03 * t2
     &    + 0.22069155D-02 ) * t2
     &    - 0.0236616773D+00 ) * t2
     &    + 0.1777582922D+00 ) * t2
     &    - 0.8888839649D+00 ) * t2
     &    + 2.6666660544D+00 ) * t2
     &    - 3.9999999710D+00 ) * t2
     &    + 1.9999999998D+00 )

        by0 = (((((((
     &    - 0.567433D-04 * t2
     &    + 0.859977D-03 ) * t2
     &    - 0.94855882D-02 ) * t2
     &    + 0.0772975809D+00 ) * t2
     &    - 0.4261737419D+00 ) * t2
     &    + 1.4216421221D+00 ) * t2
     &    - 2.3498519931D+00 ) * t2
     &    + 1.0766115157D+00 ) * t2
     &    + 0.3674669052D+00

        by0 = 2.0D+00 / pi * log ( x / 2.0D+00 ) * bj0 + by0

        by1 = ((((((((
     &      0.6535773D-03 * t2
     &    - 0.0108175626D+00 ) * t2
     &    + 0.107657606D+00 ) * t2
     &    - 0.7268945577D+00 ) * t2
     &    + 3.1261399273D+00 ) * t2
     &    - 7.3980241381D+00 ) * t2
     &    + 6.8529236342D+00 ) * t2
     &    + 0.3932562018D+00 ) * t2
     &    - 0.6366197726D+00 ) / x

        by1 = 2.0D+00 / pi * log ( x / 2.0D+00 ) * bj1 + by1

      else

        t = 4.0D+00 / x
        t2 = t * t
        a0 = dsqrt ( 2.0D+00 / ( pi * x ) )

        p0 = ((((
     &    - 0.9285D-05 * t2
     &    + 0.43506D-04 ) * t2
     &    - 0.122226D-03 ) * t2
     &    + 0.434725D-03 ) * t2
     &    - 0.4394275D-02 ) * t2
     &    + 0.999999997D+00

        q0 = t * (((((
     &      0.8099D-05 * t2
     &    - 0.35614D-04 ) * t2
     &    + 0.85844D-04 ) * t2
     &    - 0.218024D-03 ) * t2
     &    + 0.1144106D-02 ) * t2
     &    - 0.031249995D+00 )

        ta0 = x - 0.25D+00 * pi
        bj0 = a0 * ( p0 * cos ( ta0 ) - q0 * sin ( ta0 ) )
        by0 = a0 * ( p0 * sin ( ta0 ) + q0 * cos ( ta0 ) )

        p1 = ((((
     &      0.10632D-04 * t2
     &    - 0.50363D-04 ) * t2
     &    + 0.145575D-03 ) * t2
     &    - 0.559487D-03 ) * t2
     &    + 0.7323931D-02 ) * t2
     &    + 1.000000004D+00

        q1 = t * (((((
     &    - 0.9173D-05      * t2
     &    + 0.40658D-04 )   * t2
     &    - 0.99941D-04 )   * t2
     &    + 0.266891D-03 )  * t2
     &    - 0.1601836D-02 ) * t2
     &    + 0.093749994D+00 )

        ta1 = x - 0.75D+00 * pi
        bj1 = a0 * ( p1 * cos ( ta1 ) - q1 * sin ( ta1 ) )
        by1 = a0 * ( p1 * sin ( ta1 ) + q1 * cos ( ta1 ) )

      end if

      dj0 = - bj1
      dj1 = bj0 - bj1 / x
      dy0 = - by1
      dy1 = by0 - by1 / x

      return
      end
      subroutine jyna ( n, x, nm, bj, dj, by, dy )

c*********************************************************************72
c
cc JYNA computes Bessel functions Jn(x) and Yn(x) and derivatives.
c
c  Licensing:
c
c    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
c    they give permission to incorporate this routine into a user program 
c    provided that the copyright is acknowledged.
c
c  Modified:
c
c    29 April 2012
c
c  Author:
c
c    Shanjie Zhang, Jianming Jin
c
c  Reference:
c
c    Shanjie Zhang, Jianming Jin,
c    Computation of Special Functions,
c    Wiley, 1996,
c    ISBN: 0-471-11963-6,
c    LC: QA351.C45.
c
c  Parameters:
c
c    Input, integer N, the order.
c
c    Input, double precision X, the argument.
c
c    Output, integer NM, the highest order computed.
c
c    Output, double precision BJ(0:N), DJ(0:N), BY(0:N), DY(0:N), the values
c    of Jn(x), Jn'(x), Yn(x), Yn'(x).
c
      implicit none

      integer n

      double precision bj(0:n)
      double precision bj0
      double precision bj1
      double precision bjk
      double precision by(0:n)
      double precision by0
      double precision by1
      double precision cs
      double precision dj(0:n)
      double precision dj0
      double precision dj1
      double precision dy(0:n)
      double precision dy0
      double precision dy1
      double precision f
      double precision f0
      double precision f1
      double precision f2
      integer k
      integer m
      integer msta1
      integer msta2
      integer nm
      double precision x

      nm = n

      if ( x .lt. 1.0D-100 ) then

        do k = 0, n
          bj(k) = 0.0D+00
          dj(k) = 0.0D+00
          by(k) = -1.0D+300
          dy(k) = 1.0D+300
        end do
        bj(0) = 1.0D+00
        dj(1) = 0.5D+00
        return

      end if

      call jy01b ( x, bj0, dj0, bj1, dj1, by0, dy0, by1, dy1 )
      bj(0) = bj0
      bj(1) = bj1
      by(0) = by0
      by(1) = by1
      dj(0) = dj0
      dj(1) = dj1
      dy(0) = dy0
      dy(1) = dy1

      if ( n .le. 1 ) then
        return
      end if

      if ( n .lt. int ( 0.9D+00 * x) ) then

        do k = 2, n
          bjk = 2.0D+00 * ( k - 1.0D+00 ) / x * bj1 - bj0
          bj(k) = bjk
          bj0 = bj1
          bj1 = bjk
        end do

      else

        m = msta1 ( x, 200 )

        if ( m .lt. n ) then
          nm = m
        else
          m = msta2 ( x, n, 15 )
        end if

        f2 = 0.0D+00
        f1 = 1.0D-100
        do k = m, 0, -1
          f = 2.0D+00 * ( k + 1.0D+00 ) / x * f1 - f2
          if ( k .le. nm ) then
            bj(k) = f
          end if
          f2 = f1
          f1 = f
        end do

        if ( abs ( bj1 ) .lt. abs ( bj0 ) ) then
          cs = bj0 / f
        else
          cs = bj1 / f2
        end if

        do k = 0, nm
          bj(k) = cs * bj(k)
        end do

      end if

      do k = 2, nm
        dj(k) = bj(k-1) - k / x * bj(k)
      end do

      f0 = by(0)
      f1 = by(1)
      do k = 2, nm
        f = 2.0D+00 * ( k - 1.0D+00 ) / x * f1 - f0
        by(k) = f
        f0 = f1
        f1 = f
      end do

      do k = 2, nm
        dy(k) = by(k-1) - k * by(k) / x
      end do

      return
      end
      subroutine jynb ( n, x, nm, bj, dj, by, dy )

c*********************************************************************72
c
cc JYNB computes Bessel functions Jn(x) and Yn(x) and derivatives.
c
c  Licensing:
c
c    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
c    they give permission to incorporate this routine into a user program 
c    provided that the copyright is acknowledged.
c
c  Modified:
c
c    02 August 2012
c
c  Author:
c
c    Shanjie Zhang, Jianming Jin
c
c  Reference:
c
c    Shanjie Zhang, Jianming Jin,
c    Computation of Special Functions,
c    Wiley, 1996,
c    ISBN: 0-471-11963-6,
c    LC: QA351.C45.
c
c  Parameters:
c
c    Input, integer N, the order.
c
c    Input, double precision X, the argument.
c
c    Output, integer NM, the highest order computed.
c
c    Output, double precision BJ(0:N), DJ(0:N), BY(0:N), DY(0:N), the values
c    of Jn(x), Jn'(x), Yn(x), Yn'(x).
c
      implicit none

      integer n

      double precision a(4)
      double precision a1(4)
      double precision b(4)
      double precision b1(4)
      double precision bj(0:n)
      double precision bj0
      double precision bj1
      double precision bjk
      double precision bs
      double precision by(0:n)
      double precision by0
      double precision by1
      double precision byk
      double precision cu
      double precision dj(0:n)
      double precision dy(0:n)
      double precision ec
      double precision f
      double precision f1
      double precision f2
      integer k
      integer m
      integer msta1
      integer msta2
      integer nm
      double precision p0
      double precision p1
      double precision pi
      double precision q0
      double precision q1
      double precision r2p
      double precision s0
      double precision su
      double precision sv
      double precision t1
      double precision t2
      double precision x

      save a
      save a1
      save b
      save b1

      data a / -0.7031250000000000D-01, 0.1121520996093750D+00,
     &         -0.5725014209747314D+00, 0.6074042001273483D+01/

      data a1 / 0.1171875000000000D+00, -0.1441955566406250D+00,
     &          0.6765925884246826D+00, -0.6883914268109947D+01/

      data b / 0.7324218750000000D-01, -0.2271080017089844D+00,
     &         0.1727727502584457D+01, -0.2438052969955606D+02/

      data b1 /-0.1025390625000000D+00, 0.2775764465332031D+00,
     &         -0.1993531733751297D+01, 0.2724882731126854D+02/

      pi = 3.141592653589793D+00
      r2p = 0.63661977236758D+00
      nm = n

      if ( x .lt. 1.0D-100 ) then
        do k = 0, n
          bj(k) = 0.0D+00
          dj(k) = 0.0D+00
          by(k) = -1.0D+300
          dy(k) = 1.0D+300
        end do
        bj(0) = 1.0D+00
        dj(1) = 0.5D+00
        return
      end if

      if ( x .le. 300.0D+00 .or. int ( 0.9D+00 * x ) .lt. n ) then

        if ( n .eq. 0 ) then
          nm = 1
        end if

        m = msta1 ( x, 200 )

        if ( m .lt. nm ) then
          nm = m
        else
          m = msta2 ( x, nm, 15 )
        end if

        bs = 0.0D+00
        su = 0.0D+00
        sv = 0.0D+00
        f2 = 0.0D+00
        f1 = 1.0D-100

        do k = m, 0, -1
          f = 2.0D+00 * ( k + 1.0D+00 ) / x * f1 - f2
          if ( k .le. nm ) then
            bj(k) = f
          end if
          if ( k .eq. 2 * int ( k / 2 ) .and. k .ne. 0 ) then
            bs = bs + 2.0D+00 * f
            su = su + ( -1.0D+00 ) ** ( k / 2 ) * f / k
          else if ( 1 .lt. k ) then
            sv = sv + ( -1.0D+00 ) ** ( k / 2 ) * k 
     &        / ( k * k - 1.0D+00 ) * f
          end if
          f2 = f1
          f1 = f
        end do

        s0 = bs + f
        do k = 0, nm
          bj(k) = bj(k) / s0
        end do

        ec = log ( x / 2.0D+00 ) + 0.5772156649015329D+00
        by0 = r2p * ( ec * bj(0) - 4.0D+00 * su / s0 )
        by(0) = by0
        by1 = r2p * ( ( ec - 1.0D+00 ) * bj(1) - bj(0) / x 
     &    - 4.0D+00 * sv / s0 )
        by(1) = by1

      else

        t1 = x - 0.25D+00 * pi
        p0 = 1.0D+00
        q0 = -0.125D+00 / x
        do k = 1, 4
          p0 = p0 + a(k) * x ** ( - 2 * k )
          q0 = q0 + b(k) * x ** ( - 2 * k - 1 )
        end do
        cu = dsqrt ( r2p / x )
        bj0 = cu * ( p0 * cos ( t1 ) - q0 * sin ( t1 ) )
        by0 = cu * ( p0 * sin ( t1 ) + q0 * cos ( t1 ) )
        bj(0) = bj0
        by(0) = by0
        t2 = x - 0.75D+00 * pi
        p1 = 1.0D+00
        q1 = 0.375D+00 / x
        do k = 1, 4
          p1 = p1 + a1(k) * x ** ( - 2 * k )
          q1 = q1 + b1(k) * x ** ( - 2 * k - 1 )
        end do
        bj1 = cu * ( p1 * cos ( t2 ) - q1 * sin ( t2 ) )
        by1 = cu * ( p1 * sin ( t2 ) + q1 * cos ( t2 ) )
        bj(1) = bj1
        by(1) = by1
        do k = 2, nm
          bjk = 2.0D+00 * ( k - 1.0D+00 ) / x * bj1 - bj0
          bj(k) = bjk
          bj0 = bj1
          bj1 = bjk
        end do
      end if

      dj(0) = -bj(1)
      do k = 1, nm
        dj(k) = bj(k-1) - k / x * bj(k)
      end do

      do k = 2, nm
        byk = 2.0D+00 * ( k - 1.0D+00 ) * by1 / x - by0
        by(k) = byk
        by0 = by1
        by1 = byk
      end do

      dy(0) = -by(1)
      do k = 1, nm
        dy(k) = by(k-1) - k * by(k) / x
      end do

      return
      end
      subroutine jyndd ( n, x, bjn, djn, fjn, byn, dyn, fyn )

c*********************************************************************72
c
cc JYNDD computes Bessel functions Jn(x) and Yn(x), first and second derivatives.
c
c  Licensing:
c
c    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
c    they give permission to incorporate this routine into a user program 
c    provided that the copyright is acknowledged.
c
c  Modified:
c
c    07 April 2012
c
c  Author:
c
c    Shanjie Zhang, Jianming Jin
c
c  Reference:
c
c    Shanjie Zhang, Jianming Jin,
c    Computation of Special Functions,
c    Wiley, 1996,
c    ISBN: 0-471-11963-6,
c    LC: QA351.C45.
c
c  Parameters:
c
c    Input, integer N, the order.
c
c    Input, double precision X, the argument.
c
c    Output, double precision BJN, DJN, FJN, BYN, DYN, FYN, the values of
c    Jn(x), Jn'(x), Jn"(x), Yn(x), Yn'(x), Yn"(x).
c
      implicit none

      double precision bj(102)
      double precision bjn
      double precision byn
      double precision bs
      double precision by(102)
      double precision djn
      double precision dyn
      double precision e0
      double precision ec
      double precision f
      double precision f0
      double precision f1
      double precision fjn
      double precision fyn
      integer k
      integer m
      integer mt
      integer n
      integer nt
      double precision s1
      double precision su
      double precision x

      do nt = 1, 900
        mt = int ( 0.5D+00 * log10 ( 6.28D+00 * nt )
     &    - nt * log10 ( 1.36D+00 * abs ( x ) / nt ) )
        if ( 20 .lt. mt ) then
          go to 10
        end if
      end do

10    continue

      m = nt
      bs = 0.0D+00
      f0 = 0.0D+00
      f1 = 1.0D-35
      su = 0.0D+00
      do k = m, 0, -1
        f = 2.0D+00 * ( k + 1.0D+00 ) * f1 / x - f0
        if ( k .le. n + 1 ) then
          bj(k+1) = f
        end if
        if ( k .eq. 2 * int ( k / 2 ) ) then
          bs = bs + 2.0D+00 * f
          if ( k .ne. 0 ) then
            su = su + ( -1.0D+00 ) ** ( k / 2 ) * f / k
          end if
        end if
        f0 = f1
        f1 = f
      end do

      do k = 0, n + 1
        bj(k+1) = bj(k+1) / ( bs - f )
      end do

      bjn = bj(n+1)
      ec = 0.5772156649015329D+00
      e0 = 0.3183098861837907D+00
      s1 = 2.0D+00 * e0 * ( log ( x / 2.0D+00 ) + ec ) * bj(1)
      f0 = s1 - 8.0D+00 * e0 * su / ( bs - f )
      f1 = ( bj(2) * f0 - 2.0D+00 * e0 / x ) / bj(1)

      by(1) = f0
      by(2) = f1
      do k = 2, n + 1 
        f = 2.0D+00 * ( k - 1.0D+00 ) * f1 / x - f0
        by(k+1) = f
        f0 = f1
        f1 = f
      end do

      byn = by(n+1)
      djn = - bj(n+2) + n * bj(n+1) / x
      dyn = - by(n+2) + n * by(n+1) / x
      fjn = ( n * n / ( x * x ) - 1.0D+00 ) * bjn - djn / x
      fyn = ( n * n / ( x * x ) - 1.0D+00 ) * byn - dyn / x

      return
      end
      subroutine jyv ( v, x, vm, bj, dj, by, dy )

c*********************************************************************72
c
cc JYV computes Bessel functions Jv(x) and Yv(x) and their derivatives.
c
c  Licensing:
c
c    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
c    they give permission to incorporate this routine into a user program 
c    provided that the copyright is acknowledged.
c
c  Modified:
c
c    02 August 2012
c
c  Author:
c
c    Shanjie Zhang, Jianming Jin
c
c  Reference:
c
c    Shanjie Zhang, Jianming Jin,
c    Computation of Special Functions,
c    Wiley, 1996,
c    ISBN: 0-471-11963-6,
c    LC: QA351.C45.
c
c  Parameters:
c
c    Input, double precision V, the order of Jv(x) and Yv(x).
c
c    Input, double precision X, the argument of Jv(x) and Yv(x).
c
c    Output, double precision VM, the highest order computed.
c
c    Output, double precision BJ(0:N), DJ(0:N), BY(0:N), DY(0:N),
c    the values of Jn+v0(x), Jn+v0'(x), Yn+v0(x), Yn+v0'(x).
c
      implicit none

      double precision a
      double precision a0
      double precision b
      double precision bj(0:*)
      double precision bju0
      double precision bju1
      double precision bjv0
      double precision bjv1
      double precision bjvl
      double precision by(0:*)
      double precision byv0
      double precision byv1
      double precision byvk
      double precision ck
      double precision cs
      double precision cs0
      double precision cs1
      double precision dj(0:*)
      double precision dy(0:*)
      double precision ec
      double precision el
      double precision f
      double precision f0
      double precision f1
      double precision f2
      double precision ga
      double precision gb
      integer j
      integer k
      integer k0
      integer l
      integer m
      integer msta1
      integer msta2
      integer n
      double precision pi
      double precision pv0
      double precision pv1
      double precision px
      double precision qx
      double precision r
      double precision r0
      double precision r1
      double precision rp
      double precision rp2
      double precision rq
      double precision sk
      double precision v
      double precision v0
      double precision vg
      double precision vl
      double precision vm
      double precision vv
      double precision w0
      double precision w1
      double precision x
      double precision x2
      double precision xk

      el = 0.5772156649015329D+00
      pi = 3.141592653589793D+00
      rp2 = 0.63661977236758D+00
      x2 = x * x
      n = int ( v )
      v0 = v - n

      if ( x .lt. 1.0D-100 ) then

        do k = 0, n
          bj(k) = 0.0D+00
          dj(k) = 0.0D+00
          by(k) = -1.0D+300
          dy(k) = 1.0D+300
        end do

        if ( v0 .eq. 0.0D+00 ) then
          bj(0) = 1.0D+00
          dj(1) = 0.5D+00
        else
          dj(0) = 1.0D+300
        end if
        vm = v  
        return

      end if

      if ( x .le. 12.0D+00 ) then

        do l = 0, 1
          vl = v0 + l
          bjvl = 1.0D+00
          r = 1.0D+00
          do k = 1, 40
            r = -0.25D+00 * r * x2 / ( k * ( k + vl ) )
            bjvl = bjvl + r
            if ( abs ( r ) .lt. abs ( bjvl ) * 1.0D-15 ) then
              go to 10
            end if
          end do

10        continue

          vg = 1.0D+00 + vl
          call gamma ( vg, ga )
          a = ( 0.5D+00 * x ) ** vl / ga

          if ( l .eq. 0 ) then
            bjv0 = bjvl * a
          else
            bjv1 = bjvl * a
          end if

        end do

      else

        if ( x .lt. 35.0D+00 ) then
          k0 = 11
        else if ( x .lt. 50.0D+00 ) then
          k0 = 10
        else
          k0 = 8
        end if

        do j = 0, 1

          vv = 4.0D+00 * ( j + v0 ) * ( j + v0 )
          px = 1.0D+00
          rp = 1.0D+00
          do k = 1, k0
            rp = -0.78125D-02 * rp 
     &        * ( vv - ( 4.0D+00 * k - 3.0D+00 ) ** 2 ) 
     &        * ( vv - ( 4.0D+00 * k - 1.0D+00 ) ** 2 ) 
     &        / ( k * ( 2.0D+00 * k - 1.0D+00 ) * x2 )
            px = px + rp
          end do
          qx = 1.0D+00
          rq = 1.0D+00
          do k = 1, k0
            rq = -0.78125D-02 * rq 
     &        * ( vv - ( 4.0D+00 * k - 1.0D+00 ) ** 2 )
     &        * ( vv - ( 4.0D+00 * k + 1.0D+00 ) ** 2 )
     &        / ( k * ( 2.0D+00 * k + 1.0D+00 ) * x2 )
            qx = qx + rq
          end do
          qx = 0.125D+00 * ( vv - 1.0D+00 ) * qx / x
          xk = x - ( 0.5D+00 * ( j + v0 ) + 0.25D+00 ) * pi
          a0 = dsqrt ( rp2 / x )
          ck = cos ( xk )
          sk = sin ( xk )
          if ( j .eq. 0 ) then
            bjv0 = a0 * ( px * ck - qx * sk )
            byv0 = a0 * ( px * sk + qx * ck )
          else if ( j .eq. 1 ) then
            bjv1 = a0 * ( px * ck - qx * sk )
            byv1 = a0 * ( px * sk + qx * ck )
          end if

        end do

      end if

      bj(0) = bjv0
      bj(1) = bjv1
      dj(0) = v0 / x * bj(0) - bj(1)
      dj(1) = - ( 1.0D+00 + v0 ) / x * bj(1) + bj(0)

      if ( 2 .le. n .and. n .le. int ( 0.9D+00 * x ) ) then
        f0 = bjv0
        f1 = bjv1
        do k = 2, n
          f = 2.0D+00 * ( k + v0 - 1.0D+00 ) / x * f1 - f0
          bj(k) = f
          f0 = f1
          f1 = f
        end do
      else if ( 2 .le. n ) then
        m = msta1 ( x, 200 )
        if ( m .lt. n ) then
          n = m
        else
          m = msta2 ( x, n, 15 )
        end if
        f2 = 0.0D+00
        f1 = 1.0D-100
        do k = m, 0, -1
          f = 2.0D+00 * ( v0 + k + 1.0D+00 ) / x * f1 - f2
          if ( k .le. n ) then
            bj(k) = f
          end if
          f2 = f1
          f1 = f
        end do

        if ( abs ( bjv1 ) .lt. abs ( bjv0 ) ) then
          cs = bjv0 / f
        else
          cs = bjv1 / f2
        end if
        do k = 0, n
          bj(k) = cs * bj(k)
        end do
      end if

      do k = 2, n
        dj(k) = - ( k + v0 ) / x * bj(k) + bj(k-1)
      end do

      if ( x .le. 12.0D+00 ) then

        if ( v0 .ne. 0.0D+00 ) then

          do l = 0, 1

            vl = v0 + l
            bjvl = 1.0D+00
            r = 1.0D+00
            do k = 1, 40
              r = -0.25D+00 * r * x2 / ( k * ( k - vl ) )
              bjvl = bjvl + r
              if ( abs ( r ) .lt. abs ( bjvl ) * 1.0D-15 ) then
                go to 20
              end if
            end do

20          continue

            vg = 1.0D+00 - vl
            call gamma ( vg, gb )
            b = ( 2.0D+00 / x ) ** vl / gb

            if ( l .eq. 0 ) then
              bju0 = bjvl * b
            else
              bju1 = bjvl * b
            end if

          end do

          pv0 = pi * v0
          pv1 = pi * ( 1.0D+00 + v0 )
          byv0 = ( bjv0 * cos ( pv0 ) - bju0 ) / sin ( pv0 )
          byv1 = ( bjv1 * cos ( pv1 ) - bju1 ) / sin ( pv1 )

        else

          ec = log ( x / 2.0D+00 ) + el
          cs0 = 0.0D+00
          w0 = 0.0D+00
          r0 = 1.0D+00
          do k = 1, 30
            w0 = w0 + 1.0D+00 / k
            r0 = -0.25D+00 * r0 / ( k * k ) * x2
            cs0 = cs0 + r0 * w0
          end do
          byv0 = rp2 * ( ec * bjv0 - cs0 )
          cs1 = 1.0D+00
          w1 = 0.0D+00
          r1 = 1.0D+00
          do k = 1, 30
            w1 = w1 + 1.0D+00 / k
            r1 = -0.25D+00 * r1 / ( k * ( k + 1 ) ) * x2
            cs1 = cs1 + r1 * ( 2.0D+00 * w1 + 1.0D+00 
     &        / ( k + 1.0D+00 ) )
          end do
          byv1 = rp2 * ( ec * bjv1 - 1.0D+00 / x - 0.25D+00 * x * cs1 )

        end if

      end if

      by(0) = byv0
      by(1) = byv1
      do k = 2, n
        byvk = 2.0D+00 * ( v0 + k - 1.0D+00 ) / x * byv1 - byv0
        by(k) = byvk
        byv0 = byv1
        byv1 = byvk
      end do

      dy(0) = v0 / x * by(0) - by(1)
      do k = 1, n
        dy(k) = - ( k + v0 ) / x * by(k) + by(k-1)
      end do

      vm = n + v0

      return
      end
      subroutine jyzo ( n, nt, rj0, rj1, ry0, ry1 )

c*********************************************************************72
c
cc JYZO computes the zeros of Bessel functions Jn(x), Yn(x) and derivatives.
c
c  Licensing:
c
c    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
c    they give permission to incorporate this routine into a user program 
c    provided that the copyright is acknowledged.
c
c  Modified:
c
c    28 July 2012
c
c  Author:
c
c    Shanjie Zhang, Jianming Jin
c
c  Reference:
c
c    Shanjie Zhang, Jianming Jin,
c    Computation of Special Functions,
c    Wiley, 1996,
c    ISBN: 0-471-11963-6,
c    LC: QA351.C45.
c
c  Parameters:
c
c    Input, integer N, the order of the Bessel functions.
c
c    Input, integer NT, the number of zeros.
c
c    Output, double precision RJ0(NT), RJ1(NT), RY0(NT), RY1(NT), the zeros 
c    of Jn(x), Jn'(x), Yn(x), Yn'(x).
c
      implicit none

      integer nt

      double precision bjn
      double precision byn
      double precision djn
      double precision dyn
      double precision fjn
      double precision fyn
      integer l
      integer n
      double precision rj0(nt)
      double precision rj1(nt)
      double precision ry0(nt)
      double precision ry1(nt)
      double precision x
      double precision x0

      if ( n .le. 20 ) then
        x = 2.82141D+00 + 1.15859D+00 * dble ( n ) 
      else
        x = n + 1.85576D+00 * dble ( n ) ** 0.33333D+00 
     &    + 1.03315D+00 / dble ( n ) ** 0.33333D+00
      end if

      l = 0

10    continue

      x0 = x
      call jyndd ( n, x, bjn, djn, fjn, byn, dyn, fyn )
      x = x - bjn / djn
      if ( 1.0D-09 .lt. abs ( x - x0 ) ) then
        go to 10
      end if
      l = l + 1
      rj0(l) = x
      x = x + 3.1416D+00 + ( 0.0972D+00 + 0.0679D+00 * dble ( n ) 
     &  - 0.000354D+00 * dble ( n ) ** 2 ) / l

      if ( l .lt. nt ) then
        go to 10
      end if

      if ( n .le. 20 ) then
        x = 0.961587D+00 + 1.07703D+00 * dble ( n ) 
      else
        x =  dble ( n ) + 0.80861D+00 * dble ( n ) ** 0.33333D+00
     &    + 0.07249D+00 / dble ( n ) ** 0.33333D+00
      end if

      if ( n .eq. 0 ) then
        x = 3.8317D+00
      end if

      l = 0

20    continue

      x0 = x
      call jyndd ( n, x, bjn, djn, fjn, byn, dyn, fyn )
      x = x - djn / fjn
      if ( 1.0D-09 .lt. abs ( x - x0 ) ) then
        go to 20
      end if
      l = l + 1
      rj1(l) = x
      x = x + 3.1416D+00 + ( 0.4955D+00 + 0.0915D+00 * dble ( n ) 
     &  - 0.000435D+00 * dble ( n ) ** 2 ) / l

      if ( l .lt. nt ) then
        go to 20
      end if

      if ( n .le. 20 ) then
        x = 1.19477D+00 + 1.08933D+00 * dble ( n ) 
      else
        x =  dble ( n ) + 0.93158D+00 * dble ( n ) ** 0.33333D+00
     &    + 0.26035D+00 / dble ( n ) ** 0.33333D+00
      end if
 
      l = 0

30    continue

      x0 = x
      call jyndd ( n, x, bjn, djn, fjn, byn, dyn, fyn )
      x = x - byn / dyn

      if ( 1.0D-09 .lt. abs ( x - x0 ) ) then
        go to 30
      end if

      l = l + 1
      ry0(l) = x
      x = x + 3.1416D+00 + ( 0.312D+00 + 0.0852D+00 * dble ( n ) 
     &  - 0.000403D+00 * dble ( n ) ** 2 ) / l

      if ( l .lt. nt ) then
        go to 30
      end if

      if ( n .le. 20 ) then
        x = 2.67257D+00 + 1.16099D+00 * dble ( n ) 
      else
        x =  dble ( n ) + 1.8211D+00 * dble ( n ) ** 0.33333D+00
     &    + 0.94001D+00 / dble ( n ) ** 0.33333D+00
      end if
  
      l = 0

40    continue

      x0 = x
      call jyndd ( n, x, bjn, djn, fjn, byn, dyn, fyn )
      x = x - dyn / fyn
      if ( 1.0D-09 .lt. abs ( x - x0 ) ) then
        go to 40
      end if
      l = l + 1
      ry1(l) = x
      x = x + 3.1416D+00 + ( 0.197D+00 + 0.0643D+00 * dble ( n ) 
     &  -0.000286D+00 * dble ( n ) ** 2 ) / l 

      if ( l .lt. nt ) then
        go to 40
      end if

      return
      end
      subroutine klvna ( x, ber, bei, ger, gei, der, dei, her, hei )

c*********************************************************************72
c
cc KLVNA: Kelvin functions ber(x), bei(x), ker(x), and kei(x), and derivatives.
c
c  Licensing:
c
c    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
c    they give permission to incorporate this routine into a user program 
c    provided that the copyright is acknowledged.
c
c  Modified:
c
c    03 August 2012
c
c  Author:
c
c    Shanjie Zhang, Jianming Jin
c
c  Reference:
c
c    Shanjie Zhang, Jianming Jin,
c    Computation of Special Functions,
c    Wiley, 1996,
c    ISBN: 0-471-11963-6,
c    LC: QA351.C45.
c
c  Parameters:
c
c    Input, double precision X, the argument.
c
c    Output, double precision BER, BEI, GER, GEI, DER, DEI, HER, HEI, 
c    the values of ber x, bei x, ker x, kei x, ber'x, bei'x, ker'x, kei'x.
c
      implicit none

      double precision bei
      double precision ber
      double precision cn0
      double precision cp0
      double precision cs
      double precision dei
      double precision der
      double precision el
      double precision eps
      double precision fac
      double precision gei
      double precision ger
      double precision gs
      double precision hei
      double precision her
      integer k
      integer km
      integer m
      double precision pi
      double precision pn0
      double precision pn1
      double precision pp0
      double precision pp1
      double precision qn0
      double precision qn1
      double precision qp0
      double precision qp1
      double precision r
      double precision r0
      double precision r1
      double precision rc
      double precision rs
      double precision sn0
      double precision sp0
      double precision ss
      double precision x
      double precision x2
      double precision x4
      double precision xc1
      double precision xc2
      double precision xd
      double precision xe1
      double precision xe2
      double precision xt

      pi = 3.141592653589793D+00
      el = 0.5772156649015329D+00
      eps = 1.0D-15

      if ( x .eq. 0.0D+00 ) then
        ber = 1.0D+00
        bei = 0.0D+00
        ger = 1.0D+300
        gei = -0.25D+00 * pi
        der = 0.0D+00
        dei = 0.0D+00
        her = -1.0D+300
        hei = 0.0D+00
        return
      end if

      x2 = 0.25D+00 * x * x
      x4 = x2 * x2

      if ( abs ( x ) .lt. 10.0D+00 ) then

        ber = 1.0D+00
        r = 1.0D+00
        do m = 1, 60
          r = -0.25D+00 * r / ( m * m ) 
     &      / ( 2.0D+00 * m - 1.0D+00 ) ** 2 * x4
          ber = ber + r
          if ( abs ( r ) .lt. abs ( ber ) * eps ) then
            go to 10
          end if
        end do

10      continue

        bei = x2
        r = x2
        do m = 1, 60
          r = -0.25D+00 * r / ( m * m ) 
     &      / ( 2.0D+00 * m + 1.0D+00 ) ** 2 * x4
          bei = bei + r
          if ( abs ( r ) .lt. abs ( bei ) * eps ) then
            go to 20
          end if
        end do

20      continue

        ger = - ( log ( x / 2.0D+00 ) + el ) * ber + 0.25D+00 * pi * bei
        r = 1.0D+00
        gs = 0.0D+00
        do m = 1, 60
          r = -0.25D+00 * r / ( m * m ) 
     &      / ( 2.0D+00 * m - 1.0D+00 ) ** 2 * x4
          gs = gs + 1.0D+00 / ( 2.0D+00 * m - 1.0D+00 ) + 1.0D+00 
     &      / ( 2.0D+00 * m )
          ger = ger + r * gs
          if ( abs ( r * gs ) .lt. abs ( ger ) * eps ) then
            go to 30
          end if
        end do

30      continue

        gei = x2 - ( log ( x / 2.0D+00 ) + el ) * bei 
     &    - 0.25D+00 * pi * ber
        r = x2
        gs = 1.0D+00
        do m = 1, 60
          r = -0.25D+00 * r / ( m * m ) 
     &      / ( 2.0D+00 * m + 1.0D+00 ) ** 2 * x4
          gs = gs + 1.0D+00 / ( 2.0D+00 * m ) + 1.0D+00 
     &      / ( 2.0D+00 * m + 1.0D+00 )
          gei = gei + r * gs
          if ( abs ( r * gs ) .lt. abs ( gei ) * eps ) then
            go to 40
          end if
        end do

40      continue

        der = -0.25D+00 * x * x2
        r = der
        do m = 1, 60
          r = -0.25D+00 * r / m / ( m + 1.0D+00 ) 
     &      / ( 2.0D+00 * m + 1.0D+00 ) ** 2 * x4
          der = der + r
          if ( abs ( r ) .lt. abs ( der ) * eps ) then
            go to 50
          end if
        end do

50      continue

        dei = 0.5D+00 * x
        r = dei
        do m = 1, 60
          r = -0.25D+00 * r / ( m * m ) / ( 2.0D+00 * m - 1.0D+00 )
     &      / ( 2.0D+00 * m + 1.0D+00 ) * x4
          dei = dei + r
          if ( abs ( r ) .lt. abs ( dei ) * eps ) then
            go to 60
          end if
        end do

60      continue

        r = -0.25D+00 * x * x2
        gs = 1.5D+00
        her = 1.5D+00 * r - ber / x 
     &    - ( log ( x / 2.0D+00 ) + el ) * der + 0.25D+00 * pi * dei
        do m = 1, 60
          r = -0.25D+00 * r / m / ( m + 1.0D+00 ) 
     &      / ( 2.0D+00 * m + 1.0D+00 ) ** 2 * x4
          gs = gs + 1.0D+00 / ( 2 * m + 1.0D+00 ) + 1.0D+00 
     &      / ( 2 * m + 2.0D+00 )
          her = her + r * gs
          if ( abs ( r * gs ) .lt. abs ( her ) * eps ) then
            go to 70
          end if
        end do

70      continue

        r = 0.5D+00 * x
        gs = 1.0D+00
        hei = 0.5D+00 * x - bei / x 
     &    - ( log ( x / 2.0D+00 ) + el ) * dei - 0.25D+00 * pi * der
        do m = 1, 60
          r = -0.25D+00 * r / ( m * m ) / ( 2 * m - 1.0D+00 )
     &      / ( 2 * m + 1.0D+00 ) * x4
          gs = gs + 1.0D+00 / ( 2.0D+00 * m ) + 1.0D+00 
     &      / ( 2 * m + 1.0D+00 )
          hei = hei + r * gs
          if ( abs ( r * gs ) .lt. abs ( hei ) * eps ) then 
            return
          end if
        end do

      else

        pp0 = 1.0D+00
        pn0 = 1.0D+00
        qp0 = 0.0D+00
        qn0 = 0.0D+00
        r0 = 1.0D+00

        if ( abs ( x ) .lt. 40.0D+00 ) then
          km = 18
        else
          km = 10
        end if

        fac = 1.0D+00
        do k = 1, km
          fac = -fac
          xt = 0.25D+00 * k * pi - int ( 0.125D+00 * k ) * 2.0D+00 * pi
          cs = cos ( xt )
          ss = sin ( xt )
          r0 = 0.125D+00 * r0 * ( 2.0D+00 * k - 1.0D+00 ) ** 2 / k / x
          rc = r0 * cs
          rs = r0 * ss
          pp0 = pp0 + rc
          pn0 = pn0 + fac * rc
          qp0 = qp0 + rs
          qn0 = qn0 + fac * rs
        end do

        xd = x / dsqrt (2.0D+00 )
        xe1 = dexp ( xd )
        xe2 = dexp ( - xd )
        xc1 = 1.0D+00 / dsqrt ( 2.0D+00 * pi * x )
        xc2 = dsqrt ( 0.5D+00 * pi / x )
        cp0 = cos ( xd + 0.125D+00 * pi )
        cn0 = cos ( xd - 0.125D+00 * pi )
        sp0 = sin ( xd + 0.125D+00 * pi )
        sn0 = sin ( xd - 0.125D+00 * pi )
        ger = xc2 * xe2 * (  pn0 * cp0 - qn0 * sp0 )
        gei = xc2 * xe2 * ( -pn0 * sp0 - qn0 * cp0 )
        ber = xc1 * xe1 * (  pp0 * cn0 + qp0 * sn0 ) - gei / pi
        bei = xc1 * xe1 * (  pp0 * sn0 - qp0 * cn0 ) + ger / pi
        pp1 = 1.0D+00
        pn1 = 1.0D+00
        qp1 = 0.0D+00
        qn1 = 0.0D+00
        r1 = 1.0D+00
        fac = 1.0D+00

        do k = 1, km
          fac = -fac
          xt = 0.25D+00 * k * pi - int ( 0.125D+00 * k ) * 2.0D+00 * pi
          cs = cos ( xt )
          ss = sin ( xt )
          r1 = 0.125D+00 * r1 
     &      * ( 4.0D+00 - ( 2.0D+00 * k - 1.0D+00 ) ** 2 ) / k / x
          rc = r1 * cs
          rs = r1 * ss
          pp1 = pp1 + fac * rc
          pn1 = pn1 + rc
          qp1 = qp1 + fac * rs
          qn1 = qn1 + rs
        end do

        her = xc2 * xe2 * ( - pn1 * cn0 + qn1 * sn0 )
        hei = xc2 * xe2 * (   pn1 * sn0 + qn1 * cn0 )
        der = xc1 * xe1 * (   pp1 * cp0 + qp1 * sp0 ) - hei / pi
        dei = xc1 * xe1 * (   pp1 * sp0 - qp1 * cp0 ) + her / pi

      end if

      return
      end
      subroutine klvnb ( x, ber, bei, ger, gei, der, dei, her, hei )

c*********************************************************************72
c
cc KLVNB: Kelvin functions ber(x), bei(x), ker(x), and kei(x), and derivatives.
c
c  Licensing:
c
c    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
c    they give permission to incorporate this routine into a user program 
c    provided that the copyright is acknowledged.
c
c  Modified:
c
c    11 April 2012
c
c  Author:
c
c    Shanjie Zhang, Jianming Jin
c
c  Reference:
c
c    Shanjie Zhang, Jianming Jin,
c    Computation of Special Functions,
c    Wiley, 1996,
c    ISBN: 0-471-11963-6,
c    LC: QA351.C45.
c
c  Parameters:
c
c    Input, double precision X, the argument.
c
c    Output, double precision BER, BEI, GER, GEI, DER, DEI, HER, HEI, 
c    the values of ber x, bei x, ker x, kei x, ber'x, bei'x, ker'x, kei'x.
c
      implicit none

      double precision bei
      double precision ber
      double precision csn
      double precision csp
      double precision dei
      double precision der
      double precision fxi
      double precision fxr
      double precision gei
      double precision ger
      double precision hei
      double precision her
      integer l
      double precision pi
      double precision pni
      double precision pnr
      double precision ppi
      double precision ppr
      double precision ssn
      double precision ssp
      double precision t
      double precision t2
      double precision tni
      double precision tnr
      double precision tpi
      double precision tpr
      double precision u
      double precision v
      double precision x
      double precision yc1
      double precision yc2
      double precision yci
      double precision ye1
      double precision ye2
      double precision yei
      double precision yd

      pi = 3.141592653589793D+00

      if ( x .eq. 0.0D+00 ) then

        ber = 1.0D+00
        bei = 0.0D+00
        ger = 1.0D+300
        gei = -0.25D+00 * pi
        der = 0.0D+00
        dei = 0.0D+00
        her = -1.0D+300
        hei = 0.0D+00

      else if ( x .lt. 8.0D+00 ) then

        t = x / 8.0D+00
        t2 = t * t
        u = t2 * t2

        ber = ((((((
     &    - 0.901D-05 * u
     &    + 0.122552D-02 ) * u
     &    - 0.08349609D+00 ) * u
     &    + 2.64191397D+00 ) * u
     &    - 32.36345652D+00 ) * u
     &    + 113.77777774D+00 ) * u
     &    - 64.0D+00 ) * u
     &    + 1.0D+00

        bei = t * t * ((((((
     &      0.11346D-03 * u
     &    - 0.01103667D+00 ) * u
     &    + 0.52185615D+00 ) * u
     &    - 10.56765779D+00 ) * u
     &    + 72.81777742D+00 ) * u
     &    - 113.77777774D+00 ) * u
     &    + 16.0D+00 )

        ger = ((((((
     &    - 0.2458D-04 * u
     &    + 0.309699D-02 ) * u
     &    - 0.19636347D+00 ) * u
     &    + 5.65539121D+00 ) * u
     &    - 60.60977451D+00 ) * u
     &    + 171.36272133D+00 ) * u
     &    - 59.05819744D+00 ) * u
     &    - 0.57721566D+00

        ger = ger - log ( 0.5D+00 * x ) * ber + 0.25D+00 * pi * bei

        gei = t2 * ((((((
     &      0.29532D-03 * u
     &    - 0.02695875D+00 ) * u
     &    + 1.17509064D+00 ) * u
     &    - 21.30060904D+00 ) * u
     &    + 124.2356965D+00 ) * u
     &    - 142.91827687D+00 ) * u
     &    + 6.76454936D+00 )

        gei = gei - log ( 0.5D+00 * x ) * bei - 0.25D+00 * pi * ber

        der = x * t2 * ((((((
     &    - 0.394D-05 * u
     &    + 0.45957D-03 ) * u
     &    - 0.02609253D+00 ) * u
     &    + 0.66047849D+00 ) * u
     &    - 6.0681481D+00 ) * u
     &    + 14.22222222D+00 ) * u
     &    - 4.0D+00 )

        dei = x * ((((((
     &      0.4609D-04 * u
     &    - 0.379386D-02 ) * u
     &    + 0.14677204D+00 ) * u
     &    - 2.31167514D+00 ) * u
     &    + 11.37777772D+00 ) * u
     &    - 10.66666666D+00 ) * u
     &    + 0.5D+00 ) 

        her = x * t2 * ((((((
     &    - 0.1075D-04 * u
     &    + 0.116137D-02 ) * u
     &    - 0.06136358D+00 ) * u
     &    + 1.4138478D+00 ) * u
     &    - 11.36433272D+00 ) * u
     &    + 21.42034017D+00 ) * u
     &    - 3.69113734D+00 )

        her = her - log ( 0.5D+00 * x ) * der - ber / x 
     &    + 0.25D+00 * pi * dei

        hei = x * ((((((
     &      0.11997D-03 * u
     &    - 0.926707D-02 ) * u
     &    + 0.33049424D+00 ) * u
     &    - 4.65950823D+00 ) * u
     &    + 19.41182758D+00 ) * u
     &    - 13.39858846D+00 ) * u
     &    + 0.21139217D+00 )

        hei = hei - log ( 0.5D+00 * x ) * dei - bei / x 
     &    - 0.25D+00 * pi * der

      else

        t = 8.0D+00 / x

        do l = 1, 2

          v = ( -1.0D+00 ) ** l * t

          tpr = ((((
     &        0.6D-06 * v
     &      - 0.34D-05 ) * v
     &      - 0.252D-04 ) * v
     &      - 0.906D-04 ) * v * v
     &      + 0.0110486D+00 ) * v

          tpi = ((((
     &        0.19D-05 * v
     &      + 0.51D-05 ) * v * v
     &      - 0.901D-04 ) * v
     &      - 0.9765D-03 ) * v
     &      - 0.0110485D+00 ) * v
     &      - 0.3926991D+00

          if ( l .eq. 1 ) then
            tnr = tpr
            tni = tpi
          end if

        end do

        yd = x / dsqrt ( 2.0D+00 )
        ye1 = dexp ( yd + tpr )
        ye2 = dexp ( - yd + tnr )
        yc1 = 1.0D+00 / dsqrt ( 2.0D+00 * pi * x )
        yc2 = dsqrt ( pi / ( 2.0D+00 * x ) )
        csp = cos ( yd + tpi )
        ssp = sin ( yd + tpi )
        csn = cos ( - yd + tni )
        ssn = sin ( - yd + tni )
        ger = yc2 * ye2 * csn
        gei = yc2 * ye2 * ssn
        fxr = yc1 * ye1 * csp
        fxi = yc1 * ye1 * ssp
        ber = fxr - gei / pi
        bei = fxi + ger / pi

        do l = 1, 2

          v = ( -1.0D+00 ) ** l * t

          ppr = (((((
     &        0.16D-05 * v
     &      + 0.117D-04 ) * v
     &      + 0.346D-04 ) * v
     &      + 0.5D-06 ) * v
     &      - 0.13813D-02 ) * v
     &      - 0.0625001D+00 ) * v
     &      + 0.7071068D+00

          ppi = (((((
     &      - 0.32D-05 * v
     &      - 0.24D-05 ) * v
     &      + 0.338D-04 ) * v
     &      + 0.2452D-03 ) * v
     &      + 0.13811D-02 ) * v
     &      - 0.1D-06 ) * v
     &      + 0.7071068D+00

          if ( l .eq. 1 ) then
            pnr = ppr
            pni = ppi
          end if

        end do

        her =     gei * pni - ger * pnr
        hei = - ( gei * pnr + ger * pni )
        der = fxr * ppr - fxi * ppi - hei / pi
        dei = fxi * ppr + fxr * ppi + her / pi

      end if

      return
      end
      subroutine klvnzo ( nt, kd, zo )

c*********************************************************************72
c
cc KLVNZO computes zeros of the Kelvin functions.
c
c  Licensing:
c
c    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
c    they give permission to incorporate this routine into a user program 
c    provided that the copyright is acknowledged.
c
c  Modified:
c
c    15 July 2012
c
c  Author:
c
c    Shanjie Zhang, Jianming Jin
c
c  Reference:
c
c    Shanjie Zhang, Jianming Jin,
c    Computation of Special Functions,
c    Wiley, 1996,
c    ISBN: 0-471-11963-6,
c    LC: QA351.C45.
c
c  Parameters:
c
c    Input, integer NT, the number of zeros.
c
c    Input, integer KD, the function code.
c    1 for ber x, 
c    2 for bei x,
c    3 for ker x, 
c    4 for kei x,
c    5 for ber' x, 
c    6 for bei' x,
c    7 for ker' x, 
c    8 for kei' x.
c
c    Output, double precision ZO(NT), the zeros of the given Kelvin function.
c
      implicit none

      integer nt

      double precision bei
      double precision ber
      double precision ddi
      double precision ddr
      double precision dei
      double precision der
      double precision gdi
      double precision gdr
      double precision gei
      double precision ger
      double precision hei
      double precision her
      integer kd
      integer m
      double precision rt
      double precision rt0(8)
      double precision zo(nt)

      rt0(1) = 2.84891D+00
      rt0(2) = 5.02622D+00
      rt0(3) = 1.71854D+00
      rt0(4) = 3.91467D+00
      rt0(5) = 6.03871D+00
      rt0(6) = 3.77268D+00
      rt0(7) = 2.66584D+00
      rt0(8) = 4.93181D+00

      rt = rt0(kd)

      do m = 1,nt

10      continue

        call klvna ( rt, ber, bei, ger, gei, der, dei, her, hei )

        if ( kd .eq. 1 ) then
          rt = rt - ber / der
        else if ( kd .eq. 2 ) then
          rt = rt - bei / dei
        else if ( kd .eq. 3 ) then
          rt = rt - ger / her
        else if ( kd .eq. 4 ) then
          rt = rt - gei / hei
        else if ( kd .eq. 5 ) then
          ddr = - bei - der / rt
          rt = rt - der / ddr
        else if ( kd .eq. 6 ) then
          ddi = ber - dei / rt
          rt = rt - dei / ddi
        else if ( kd .eq. 7 ) then
          gdr = - gei - her / rt
          rt = rt - her / gdr
        else
          gdi = ger - hei / rt
          rt = rt - hei / gdi
        end if

        if ( 5.0D-10 .lt. abs ( rt - rt0(kd) ) ) then
          rt0(kd) = rt
          go to 10
        end if

        zo(m) = rt
        rt = rt + 4.44D+00

      end do

      return
      end
      subroutine kmn ( m, n, c, cv, kd, df, dn, ck1, ck2 )

c*********************************************************************72
c
cc KMN: expansion coefficients of prolate or oblate spheroidal functions.
c
c  Licensing:
c
c    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
c    they give permission to incorporate this routine into a user program 
c    provided that the copyright is acknowledged.
c
c  Modified:
c
c    02 August 2012
c
c  Author:
c
c    Shanjie Zhang, Jianming Jin
c
c  Reference:
c
c    Shanjie Zhang, Jianming Jin,
c    Computation of Special Functions,
c    Wiley, 1996,
c    ISBN: 0-471-11963-6,
c    LC: QA351.C45.
c
c  Parameters:
c
c    Input, integer M, the mode parameter;  M = 0, 1, 2, ...
c
c    Input, integer N, mode parameter, N = M, M + 1, M + 2, ...
c
c    Input, double precision C, spheroidal parameter.
c
c    Input, double precision CV, the characteristic value.
c
c    Input, integer KD, the function code.
c    1, the prolate function.
c    -1, the oblate function.
c
c    Input, double precision DF(*), the expansion coefficients.
c
      implicit none

      double precision c
      double precision ck1
      double precision ck2
      double precision cs
      double precision cv
      double precision df(200)
      double precision dn(200)
      double precision dnp
      double precision g0
      double precision gk0
      double precision gk1
      double precision gk2
      double precision gk3
      integer i
      integer ip
      integer j
      integer k
      integer kd
      integer l
      integer m
      integer n
      integer nm
      integer nm1
      integer nn
      double precision r
      double precision r1
      double precision r2
      double precision r3
      double precision r4
      double precision r5
      double precision rk(200)
      double precision sa0
      double precision sb0
      double precision su0
      double precision sw
      double precision t
      double precision tp(200)
      double precision u(200)
      double precision v(200)
      double precision w(200)

      nm = 25 + int ( 0.5D+00 * ( n - m ) + c )
      nn = nm + m
      cs = c * c * kd

      if ( n - m .eq. 2 * int ( ( n - m ) / 2 ) ) then
        ip = 0
      else
        ip = 1
      end if

      do i = 1, nn + 3

        if ( ip .eq. 0 ) then
          k = - 2 * ( i - 1 )
        else
          k = - ( 2 * i - 3 )
        end if

        gk0 = 2.0D+00 * m + k
        gk1 = ( m + k ) * ( m + k + 1.0D+00 )
        gk2 = 2.0D+00 * ( m + k ) - 1.0D+00
        gk3 = 2.0D+00 * ( m + k ) + 3.0D+00
        u(i) = gk0 * ( gk0 - 1.0D+00 ) * cs 
     &    / ( gk2 * ( gk2 + 2.0D+00 ) )
        v(i) = gk1 - cv + ( 2.0D+00 * ( gk1 - m * m ) - 1.0D+00 ) * cs 
     &    / ( gk2 * gk3 )
        w(i) = ( k + 1.0D+00 ) * ( k + 2.0D+00 ) * cs 
     &    / ( ( gk2 + 2.0D+00 ) * gk3 )

      end do

      do k = 1, m
        t = v(m+1)
        do l = 0, m - k - 1
          t = v(m-l) - w(m-l+1) * u(m-l) / t
        end do
        rk(k) = -u(k) / t
      end do

      r = 1.0D+00
      do k = 1, m
        r = r * rk(k)
        dn(k) = df(1) * r
      end do

      tp(nn) = v(nn+1)
      do k = nn - 1, m + 1,-1
        tp(k) = v(k+1) - w(k+2) * u(k+1) / tp(k+1)
        if ( m + 1 .lt. k ) then
          rk(k) = -u(k) / tp(k)
        end if
      end do

      if ( m .eq. 0 ) then
        dnp = df(1)
      else
        dnp = dn(m)
      end if

      dn(m+1) = ( - 1.0D+00 ) ** ip * dnp * cs 
     &  / ( ( 2.0D+00 * m - 1.0D+00 ) 
     &  * ( 2.0D+00 * m + 1.0D+00 - 4.0D+00 * ip ) * tp(m+1) )
      do k = m + 2, nn
        dn(k) = rk(k) * dn(k-1)
      end do

      r1 = 1.0D+00
      do j = 1, ( n + m + ip ) / 2
           r1 = r1 * ( j + 0.5D+00 * ( n + m + ip ) )
      end do
      nm1 = ( n - m ) / 2
      r = 1.0D+00
      do j = 1, 2 * m + ip
        r = r * j
      end do
      su0 = r * df(1)

      do k = 2, nm
        r = r * ( m + k - 1.0D+00 ) * ( m + k + ip - 1.5D+00 )
     &    / ( k - 1.0D+00 ) / ( k + ip - 1.5D+00 )
        su0 = su0 + r * df(k)
        if ( nm1 .lt. k .and. 
     &    abs ( ( su0 - sw ) / su0 ) .lt. 1.0D-14 ) then
          go to 10
        end if
        sw = su0
      end do

10    continue

      if ( kd .eq. 1 ) then
        go to 20
      end if

      r2 = 1.0D+00
      do j = 1,m
        r2 = 2.0D+00 * c * r2 * j
      end do
      r3 = 1.0D+00
      do j = 1, ( n - m - ip ) / 2
        r3 = r3 * j
      end do
      sa0 = ( 2.0D+00 * ( m + ip ) + 1.0D+00 ) * r1 
     &  / ( 2.0D+00 ** n * c ** ip * r2 * r3 * df(1) )
      ck1 = sa0 * su0

      if ( kd .eq. -1 ) then
        return
      end if

20    continue

      r4 = 1.0D+00
      do j = 1, ( n - m - ip ) / 2
        r4 = 4.0D+00 * r4 * j
      end do
      r5 = 1.0D+00
      do j = 1, m
        r5 = r5 * ( j + m ) / c
      end do

      if ( m .eq. 0 ) then
        g0 = df(1)
      else
        g0 = dn(m)
      end if

      sb0 = ( ip + 1.0D+00 ) * c ** ( ip + 1 ) 
     &  / ( 2.0D+00 * ip * ( m - 2.0D+00 ) + 1.0D+00 ) 
     &  / ( 2.0D+00 * m - 1.0D+00 )

      ck2 = ( -1 ) ** ip * sb0 * r4 * r5 * g0 / r1 * su0

      return
      end
      subroutine lagzo ( n, x, w )

c*********************************************************************72
c
cc LAGZO computes zeros of the Laguerre polynomial, and integration weights.
c
c  Discussion:
c
c    This procedure computes the zeros of Laguerre polynomial Ln(x) in the 
c    interval [0,Ï], and the corresponding weighting coefficients for Gauss-Laguerre
c    integration.
c
c  Licensing:
c
c    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
c    they give permission to incorporate this routine into a user program 
c    provided that the copyright is acknowledged.
c
c  Modified:
c
c    07 April 2012
c
c  Author:
c
c    Shanjie Zhang, Jianming Jin
c
c  Reference:
c
c    Shanjie Zhang, Jianming Jin,
c    Computation of Special Functions,
c    Wiley, 1996,
c    ISBN: 0-471-11963-6,
c    LC: QA351.C45.
c
c  Parameters:
c
c    Input, integer N, the order of the Laguerre polynomial.
c
c    Output, double precision X(N), the zeros of the Laguerre polynomial.
c
c    Output, double precision W(N), the weighting coefficients.
c
      implicit none

      integer n

      double precision f0
      double precision f1
      double precision fd
      double precision gd
      double precision hn
      integer i
      integer it
      integer j
      integer k
      integer nr
      double precision p
      double precision pd
      double precision pf
      double precision q
      double precision w(n)
      double precision wp
      double precision x(n)
      double precision z
      double precision z0

      hn = 1.0D+00 / dble ( n )

      do nr = 1, n

        if ( nr .eq. 1 ) then
          z = hn
        else
          z = x(nr-1) + hn * nr ** 1.27D+00
        end if

        it = 0

10      continue

        it = it + 1
        z0 = z
        p = 1.0D+00
        do i = 1, nr - 1
          p = p * ( z - x(i) )
        end do

        f0 = 1.0D+00
        f1 = 1.0D+00 - z
        do k = 2, n
          pf = (( 2.0D+00 * k - 1.0D+00 - z ) * f1 
     &      - ( k - 1.0D+00 ) * f0 ) / k
          pd = k / z * ( pf - f1 )
          f0 = f1
          f1 = pf
        end do

        fd = pf / p

        q = 0.0D+00
        do i = 1, nr - 1
          wp = 1.0D+00
          do j = 1, nr - 1
            if ( j .ne. i ) then
              wp = wp * ( z - x(j) )
            end if
          end do
          q = q + wp
        end do

        gd = ( pd - q * fd ) / p
        z = z - fd / gd

        if ( it .le. 40 .and. abs ( ( z - z0 ) / z ) .gt. 1.0D-15 ) then
          go to 10
        end if

        x(nr) = z
        w(nr) = 1.0D+00 / ( z * pd * pd )

      end do

      return
      end
      subroutine lamn ( n, x, nm, bl, dl )

c*********************************************************************72
c
cc LAMN computes lambda functions and derivatives.
c
c  Licensing:
c
c    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
c    they give permission to incorporate this routine into a user program 
c    provided that the copyright is acknowledged.
c
c  Modified:
c
c    14 July 2012
c
c  Author:
c
c    Shanjie Zhang, Jianming Jin
c
c  Reference:
c
c    Shanjie Zhang, Jianming Jin,
c    Computation of Special Functions,
c    Wiley, 1996,
c    ISBN: 0-471-11963-6,
c    LC: QA351.C45.
c
c  Parameters:
c
c    Input, integer N, the order.
c
c    Input, double precision X, the argument.
c
c    Output, integer NM, the highest order computed.
c
c    Output, double precision BL(0:N), DL(0:N), the
c    value of the lambda function and its derivative of orders 0 through N.
c
      implicit none

      integer n

      double precision bg
      double precision bk
      double precision bl(0:n)
      double precision bs
      double precision dl(0:n)
      double precision f
      double precision f0
      double precision f1
      integer i
      integer k
      integer m
      integer msta1
      integer msta2
      integer nm
      double precision r
      double precision r0
      double precision uk
      double precision x
      double precision x2

      nm = n

      if ( abs ( x ) .lt. 1.0D-100 ) then
        do k = 0, n
          bl(k) = 0.0D+00
          dl(k) = 0.0D+00
        end do
        bl(0) = 1.0D+00
        dl(1) = 0.5D+00
        return
      end if

      if ( x .le. 12.0D+00 ) then

        x2 = x * x

        do k = 0, n
          bk = 1.0D+00
          r = 1.0D+00
          do i = 1, 50
            r = -0.25D+00 * r * x2 / ( i * ( i + k ) )
            bk = bk + r
            if ( abs ( r ) .lt. abs ( bk ) * 1.0D-15 ) then
              go to 10
            end if
          end do

10        continue

          bl(k) = bk
          if ( 1 .le. k ) then
            dl(k-1) = - 0.5D+00 * x / k * bk
          end if

        end do

        uk = 1.0D+00
        r = 1.0D+00
        do i = 1, 50
          r = -0.25D+00 * r * x2 / ( i * ( i + n + 1.0D+00 ) )
          uk = uk + r
          if ( abs ( r ) .lt. abs ( uk ) * 1.0D-15 ) then
            go to 20
          end if
        end do

20      continue

        dl(n) = -0.5D+00 * x / ( n + 1.0D+00 ) * uk
        return

      end if

      if ( n .eq. 0 ) then
        nm = 1
      end if

      m = msta1 ( x, 200 )

      if ( m .lt. nm ) then
        nm = m
      else
        m = msta2 ( x, nm, 15 )
      end if

      bs = 0.0D+00
      f0 = 0.0D+00
      f1 = 1.0D-100
      do k = m, 0, -1
        f = 2.0D+00 * ( k + 1.0D+00 ) * f1 / x - f0
        if ( k .le. nm ) then
          bl(k) = f
        end if
        if ( k .eq. 2 * int ( k / 2 ) ) then
          bs = bs + 2.0D+00 * f
        end if
        f0 = f1
        f1 = f
      end do

      bg = bs - f
      do k = 0, nm
        bl(k) = bl(k) / bg
      end do

      r0 = 1.0D+00
      do k = 1, nm
        r0 = 2.0D+00 * r0 * k / x
        bl(k) = r0 * bl(k)
      end do

      dl(0) = -0.5D+00 * x * bl(1)
      do k = 1, nm
        dl(k) = 2.0D+00 * k / x * ( bl(k-1) - bl(k) )
      end do

      return
      end
      subroutine lamv ( v, x, vm, vl, dl )

c*********************************************************************72
c
cc LAMV computes lambda functions and derivatives of arbitrary order.
c
c  Licensing:
c
c    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
c    they give permission to incorporate this routine into a user program 
c    provided that the copyright is acknowledged.
c
c  Modified:
c
c    31 July 2012
c
c  Author:
c
c    Shanjie Zhang, Jianming Jin
c
c  Reference:
c
c    Shanjie Zhang, Jianming Jin,
c    Computation of Special Functions,
c    Wiley, 1996,
c    ISBN: 0-471-11963-6,
c    LC: QA351.C45.
c
c  Parameters:
c
c    Input, double precision V, the order.
c
c    Input, double precision X, the argument.
c
c    Output, double precision VM, the highest order computed.
c
c    Output, double precision VL(0:*), DL(0:*), the Lambda function and 
c    derivative, of orders N+V0.
c
      implicit none

      double precision v

      double precision a0
      double precision bjv0
      double precision bjv1
      double precision bk
      double precision ck
      double precision cs
      double precision dl(0:int(v))
      double precision f
      double precision f0
      double precision f1
      double precision f2
      double precision fac
      double precision ga
      integer i
      integer j
      integer k
      integer k0
      integer m
      integer msta1
      integer msta2
      integer n
      double precision pi
      double precision px
      double precision qx
      double precision r
      double precision r0
      double precision rc
      double precision rp
      double precision rp2
      double precision rq
      double precision sk
      double precision uk
      double precision v0
      double precision vk
      double precision vl(0:int(v))
      double precision vm
      double precision vv
      double precision x
      double precision x2
      double precision xk

      pi = 3.141592653589793D+00
      rp2 = 0.63661977236758D+00
      x = abs ( x )
      x2 = x * x
      n = int ( v )
      v0 = v - n
      vm = v

      if ( x .le. 12.0D+00 ) then

        do k = 0, n

          vk = v0 + k
          bk = 1.0D+00
          r = 1.0D+00

          do i = 1, 50
            r = -0.25D+00 * r * x2 / ( i * ( i + vk ) )
            bk = bk + r
            if ( abs ( r ) .lt. abs ( bk ) * 1.0D-15 ) then
              go to 10
            end if
          end do

10        continue

          vl(k) = bk
          uk = 1.0D+00
          r = 1.0D+00
          do i = 1, 50
            r = -0.25D+00 * r * x2 / ( i * ( i + vk + 1.0D+00 ))
            uk = uk + r
            if ( abs ( r ) .lt. abs ( uk ) * 1.0D-15 ) then
              go to 20
            end if
          end do

20        continue

          dl(k) = - 0.5D+00 * x / ( vk + 1.0D+00 ) * uk

        end do

        return

      end if

      if ( x .lt. 35.0D+00 ) then
        k0 = 11
      else if ( x .lt. 50.0D+00 ) then
        k0 = 10
      else
        k0 = 8
      end if

      do j = 0, 1
        vv = 4.0D+00 * ( j + v0 ) * ( j + v0 )
        px = 1.0D+00
        rp = 1.0D+00
        do k = 1, k0
          rp = - 0.78125D-02 * rp 
     &      * ( vv - ( 4.0D+00 * k - 3.0D+00 ) ** 2 ) 
     &      * ( vv - ( 4.0D+00 * k - 1.0D+00 ) ** 2 ) 
     &      / ( k * ( 2.0 * k - 1.0D+00 ) * x2 )
          px = px + rp
        end do
        qx = 1.0D+00
        rq = 1.0D+00
        do k = 1, k0
          rq = - 0.78125D-02 * rq 
     &      * ( vv - ( 4.0D+00 * k - 1.0D+00 ) ** 2 )
     &      * ( vv - ( 4.0D+00 * k + 1.0D+00 ) ** 2 )
     &      / ( k * ( 2.0D+00 * k + 1.0D+00 ) * x2 )
          qx = qx + rq
        end do
        qx = 0.125D+00 * ( vv - 1.0D+00 ) * qx / x
        xk = x - ( 0.5D+00 * ( j + v0 ) + 0.25D+00 ) * pi
        a0 = sqrt ( rp2 / x )
        ck = cos ( xk )
        sk = sin ( xk )
        if ( j .eq. 0 ) then
          bjv0 = a0 * ( px * ck - qx * sk )
        else
          bjv1 = a0 * ( px * ck - qx * sk )
        end if
      end do

      if ( v0 .eq. 0.0D+00 ) then
        ga = 1.0D+00
      else
        call gam0 ( v0, ga )
        ga = v0 * ga
      end if

      fac = ( 2.0D+00 / x ) ** v0 * ga
      vl(0) = bjv0
      dl(0) = - bjv1 + v0 / x * bjv0
      vl(1) = bjv1
      dl(1) = bjv0 - ( 1.0D+00 + v0 ) / x * bjv1
      r0 = 2.0D+00 * ( 1.0D+00 + v0 ) / x

      if ( n .le. 1 ) then
        vl(0) = fac * vl(0)
        dl(0) = fac * dl(0) - v0 / x * vl(0)
        vl(1) = fac * r0 * vl(1)
        dl(1) = fac * r0 * dl(1) - ( 1.0D+00 + v0 ) / x * vl(1)
        return
      end if

      if ( 2 .le. n .and. n .le. int ( 0.9D+00 * x ) ) then
        f0 = bjv0
        f1 = bjv1
        do k = 2, n
          f = 2.0D+00 * ( k + v0 - 1.0D+00 ) / x * f1 - f0
          f0 = f1
          f1 = f
          vl(k) = f
        end do
      else if ( 2 .le. n ) then
        m = msta1 ( x, 200 )
        if ( m .lt. n ) then
          n = m
        else
          m = msta2 ( x, n, 15 )
        end if
        f2 = 0.0D+00
        f1 = 1.0D-100
        do k = m, 0, -1
          f = 2.0D+00 * ( v0 + k + 1.0D+00 ) / x * f1 - f2
          if ( k .le. n ) then
            vl(k) = f
          end if
          f2 = f1
          f1 = f
        end do

        if ( abs ( bjv0 ) .le. abs ( bjv1 ) ) then
          cs = bjv1 / f2
        else
          cs = bjv0 / f
        end if

        do k = 0, n
          vl(k) = cs * vl(k)
        end do

      end if

      vl(0) = fac * vl(0)
      do j = 1, n
        rc = fac * r0
        vl(j) = rc * vl(j)
        dl(j-1) = - 0.5D+00 * x / ( j + v0 ) * vl(j)
        r0 = 2.0D+00 * ( j + v0 + 1 ) / x * r0
      end do
      dl(n) = 2.0D+00 * ( v0 + n ) * ( vl(n-1) - vl(n) ) / x
      vm = n + v0

      return
      end
      subroutine legzo ( n, x, w )

c*********************************************************************72
c
cc LEGZO computes the zeros of Legendre polynomials, and integration weights.
c
c  Discussion:
c
c    This procedure computes the zeros of Legendre polynomial Pn(x) in the interval 
c    [-1,1], and the corresponding weighting coefficients for Gauss-Legendre
c    integration.
c
c  Licensing:
c
c    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
c    they give permission to incorporate this routine into a user program 
c    provided that the copyright is acknowledged.
c
c  Modified:
c
c    13 July 2012
c
c  Author:
c
c    Shanjie Zhang, Jianming Jin
c
c  Reference:
c
c    Shanjie Zhang, Jianming Jin,
c    Computation of Special Functions,
c    Wiley, 1996,
c    ISBN: 0-471-11963-6,
c    LC: QA351.C45.
c
c  Parameters:
c
c    Input, integer N, the order of the polynomial.
c
c    Output, double precision X(N), W(N), the zeros of the polynomial,
c    and the corresponding weights.
c
      implicit none

      integer n

      double precision f0
      double precision f1
      double precision fd
      double precision gd
      integer i
      integer j
      integer k
      integer n0
      integer nr
      double precision p
      double precision pd
      double precision pf
      double precision q
      double precision w(n)
      double precision wp
      double precision x(n)
      double precision z
      double precision z0

      n0 = ( n + 1 ) / 2

      do nr = 1, n0

        z = cos ( 3.1415926D+00 * ( nr - 0.25D+00 ) / n )

10      continue

        z0 = z
        p = 1.0D+00
        do i = 1, nr - 1
          p = p * ( z - x(i))
        end do
        f0 = 1.0D+00
        if ( nr .eq. n0 .and. n .ne. 2 * int ( n / 2 ) ) then
          z = 0.0D+00
        end if
        f1 = z
        do k = 2, n
          pf = ( 2.0D+00 - 1.0D+00 / k ) * z * f1 
     &      - ( 1.0D+00 - 1.0D+00 / k ) * f0
          pd = k * ( f1 - z * pf ) / ( 1.0D+00 - z * z )
          f0 = f1
          f1 = pf
        end do

        if ( z .ne. 0.0D+00 ) then

          fd = pf / p
          q = 0.0D+00
          do i = 1, nr - 1
            wp = 1.0D+00
            do j = 1, nr - 1
              if ( j .ne. i ) then
                wp = wp * ( z - x(j) )
              end if
            end do
            q = q + wp
          end do
          gd = ( pd - q * fd ) / p
          z = z - fd / gd
          if ( abs ( z - z0 ) .gt. abs ( z ) * 1.0D-15 ) then
            go to 10
          end if

        end if

        x(nr) = z
        x(n+1-nr) = - z
        w(nr) = 2.0D+00 / ( ( 1.0D+00 - z * z ) * pd * pd )
        w(n+1-nr) = w(nr)

      end do

      return
      end
      subroutine lgama ( kf, x, gl )

c*********************************************************************72
c
cc LGAMA computes the gamma function or its logarithm.
c
c  Licensing:
c
c    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
c    they give permission to incorporate this routine into a user program 
c    provided that the copyright is acknowledged.
c
c  Modified:
c
c    15 July 2012
c
c  Author:
c
c    Shanjie Zhang, Jianming Jin
c
c  Reference:
c
c    Shanjie Zhang, Jianming Jin,
c    Computation of Special Functions,
c    Wiley, 1996,
c    ISBN: 0-471-11963-6,
c    LC: QA351.C45.
c
c  Parameters:
c
c    Input, integer KF, the argument code.
c    1, for gamma(x);
c    2, for ln(gamma(x)).
c
c    Input, double precision X, the argument.
c
c    Output, double precision GL, the function value.
c
      implicit none

      double precision a(10)
      double precision gl
      double precision gl0
      integer k
      integer kf
      integer n
      double precision x
      double precision x0
      double precision x2
      double precision xp

      save a

      data a / 8.333333333333333D-02, -2.777777777777778D-03,
     &         7.936507936507937D-04, -5.952380952380952D-04,
     &         8.417508417508418D-04, -1.917526917526918D-03,
     &         6.410256410256410D-03, -2.955065359477124D-02,
     &         1.796443723688307D-01, -1.39243221690590D+00 /

      x0 = x

      if ( x .eq. 1.0D+00 .or. x .eq. 2.0D+00 ) then
        gl = 0.0D+00
        if ( kf .eq. 1 ) then
          gl = 1.0D+00
        end if
        return
      else if ( x .le. 7.0D+00 ) then
        n = int ( 7.0D+00 - x )
        x0 = x + n
      end if

      x2 = 1.0D+00 / ( x0 * x0 )
      xp = 6.283185307179586477D+00
      gl0 = a(10)

      do k = 9, 1, -1
        gl0 = gl0 * x2 + a(k)
      end do

      gl = gl0 / x0 + 0.5D+00 * log ( xp ) 
     &  + ( x0 - 0.5D+00 ) * log ( x0 ) - x0

      if ( x .le. 7.0D+00 ) then
        do k = 1, n
          gl = gl - log ( x0 - 1.0D+00 )
          x0 = x0 - 1.0D+00
        end do
      end if

      if ( kf .eq. 1 ) then
        gl = exp ( gl )
      end if

      return
      end
      subroutine lpmn ( mm, m, n, x, pm, pd )

c*********************************************************************72
c
cc LPMN computes associated Legendre functions Pmn(X) and derivatives P'mn(x).
c
c  Licensing:
c
c    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
c    they give permission to incorporate this routine into a user program 
c    provided that the copyright is acknowledged.
c
c  Modified:
c
c    19 July 2012
c
c  Author:
c
c    Shanjie Zhang, Jianming Jin
c
c  Reference:
c
c    Shanjie Zhang, Jianming Jin,
c    Computation of Special Functions,
c    Wiley, 1996,
c    ISBN: 0-471-11963-6,
c    LC: QA351.C45.
c
c  Parameters:
c
c    Input, integer MM, the leading dimension of PM and PD.
c
c    Input, integer M, the order of Pmn(x).
c
c    Input, integer N, the degree of Pmn(x).
c
c    Input, double precision X, the argument of Pmn(x).
c
c    Output, double precision PM(0:MM,0:N), PD(0:MM,0:N), the
c    values of Pmn(x) and Pmn'(x).
c
      implicit none

      integer mm
      integer n

      integer i
      integer j
      integer ls
      integer m
      double precision pd(0:mm,0:n)
      double precision pm(0:mm,0:n)
      double precision x
      double precision xq
      double precision xs

      do i = 0, n
        do j = 0, m
          pm(j,i) = 0.0D+00
          pd(j,i) = 0.0D+00
        end do
      end do

      pm(0,0) = 1.0D+00

      if ( abs ( x ) .eq. 1.0D+00 ) then

        do i = 1, n
          pm(0,i) = x ** i
          pd(0,i) = 0.5D+00 * i * ( i + 1.0D+00 ) * x ** ( i + 1 )
        end do

        do j = 1, n
          do i = 1, m
            if ( i .eq. 1 ) then
              pd(i,j) = 1.0D+300
            else if ( i .eq. 2 ) then
              pd(i,j) = -0.25D+00 * ( j + 2 ) * ( j + 1 ) * j 
     &          * ( j - 1 ) * x ** ( j + 1 )
            end if
          end do
        end do

        return

      end if

      if ( 1.0D+00 .lt. abs ( x ) ) then
        ls = -1
      else
        ls = +1
      end if

      xq = sqrt ( ls * ( 1.0D+00 - x * x ) )
      xs = ls * ( 1.0D+00 - x * x )
      do i = 1, m
        pm(i,i) = - ls * ( 2.0D+00 * i - 1.0D+00 ) * xq * pm(i-1,i-1)
      end do

      do i = 0, m
        pm(i,i+1) = ( 2.0D+00 * i + 1.0D+00 ) * x * pm(i,i)
      end do

      do i = 0, m
        do j = i + 2, n
          pm(i,j) = ( ( 2.0D+00 * j - 1.0D+00 ) * x * pm(i,j-1) -
     &      ( i + j - 1.0D+00 ) * pm(i,j-2) ) / ( j - i )
        end do
      end do

      pd(0,0) = 0.0D+00
      do j = 1, n
        pd(0,j) = ls * j * ( pm(0,j-1) - x * pm(0,j) ) / xs
      end do

      do i = 1, m
        do j = i, n
          pd(i,j) = ls * i * x * pm(i,j) / xs + ( j + i )
     &      * ( j - i + 1.0D+00 ) / xq * pm(i-1,j)
        end do
      end do

      return
      end
      subroutine lpmns ( m, n, x, pm, pd )

c*********************************************************************72
c
cc LPMNS computes associated Legendre functions Pmn(X) and derivatives P'mn(x).
c
c  Licensing:
c
c    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
c    they give permission to incorporate this routine into a user program 
c    provided that the copyright is acknowledged.
c
c  Modified:
c
c    18 July 2012
c
c  Author:
c
c    Shanjie Zhang, Jianming Jin
c
c  Reference:
c
c    Shanjie Zhang, Jianming Jin,
c    Computation of Special Functions,
c    Wiley, 1996,
c    ISBN: 0-471-11963-6,
c    LC: QA351.C45.
c
c  Parameters:
c
c    Input, integer M, the order of Pmn(x).
c
c    Input, integer N, the degree of Pmn(x).
c
c    Input, double precision X, the argument.
c
c    Output, double precision PM(0:N), PD(0:N), the values and derivatives
c    of the function from degree 0 to N.
c
      implicit none

      integer n

      integer k
      integer m
      double precision pm(0:n)
      double precision pm0
      double precision pm1
      double precision pm2
      double precision pmk
      double precision pd(0:n)
      double precision x
      double precision x0

      do k = 0, n
        pm(k) = 0.0D+00
        pd(k) = 0.0D+00
      end do

      if ( abs ( x ) .eq. 1.0D+00 ) then

        do k = 0, n
          if ( m .eq. 0 ) then
            pm(k) = 1.0D+00
            pd(k) = 0.5D+00 * k * ( k + 1.0D+00 )
            if ( x .lt. 0.0D+00 ) then
              pm(k) = ( -1.0D+00 ) ** k * pm(k)
              pd(k) = ( -1.0D+00 ) ** ( k + 1 ) * pd(k)
            end if
          else if ( m .eq. 1 ) then
            pd(k) = 1.0D+300
          else if ( m .eq. 2 ) then
            pd(k) = -0.25D+00 * ( k + 2.0D+00 ) * ( k + 1.0D+00 ) 
     &        * k * ( k - 1.0D+00 )
            if ( x .lt. 0.0D+00 ) then
              pd(k) = ( -1.0D+00 ) ** ( k + 1 ) * pd(k)
            end if
          end if
        end do
        return
      end if

      x0 = abs ( 1.0D+00 - x * x )
      pm0 = 1.0D+00
      pmk = pm0
      do k = 1, m
        pmk = ( 2.0D+00 * k - 1.0D+00 ) * sqrt ( x0 ) * pm0
        pm0 = pmk
      end do
      pm1 = ( 2.0D+00 * m + 1.0D+00 ) * x * pm0
      pm(m) = pmk
      pm(m+1) = pm1
      do k = m + 2, n
        pm2 = ( ( 2.0D+00 * k - 1.0D+00 ) * x * pm1 
     &    - ( k + m - 1.0D+00 ) * pmk ) / ( k - m )
        pm(k) = pm2
        pmk = pm1
        pm1 = pm2
      end do

      pd(0) = ( ( 1.0D+00 - m ) * pm(1) - x * pm(0) ) 
     &  / ( x * x - 1.0D+00 )  
      do k = 1, n
        pd(k) = ( k * x * pm(k) - ( k + m ) * pm(k-1) ) 
     &    / ( x * x - 1.0D+00 )
      end do

      return
      end
      subroutine lpmv ( v, m, x, pmv )

c*********************************************************************72
c
cc LPMV computes associated Legendre functions Pmv(X) with arbitrary degree.
c
c  Discussion:
c
c    Compute the associated Legendre function Pmv(x) with an integer order 
c    and an arbitrary nonnegative degree v.
c
c  Licensing:
c
c    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
c    they give permission to incorporate this routine into a user program 
c    provided that the copyright is acknowledged.
c
c  Modified:
c
c    19 July 2012
c
c  Author:
c
c    Shanjie Zhang, Jianming Jin
c
c  Reference:
c
c    Shanjie Zhang, Jianming Jin,
c    Computation of Special Functions,
c    Wiley, 1996,
c    ISBN: 0-471-11963-6,
c    LC: QA351.C45.
c
c  Parameters:
c
c    Input, double precision V, the degree of Pmv(x).
c
c    Input, integer M, the order of Pmv(x).
c
c    Input, double precision X, the argument of Pm(x).
c
c    Output, double precision PMV, the value of Pm(x).
c
      implicit none

      double precision c0
      double precision el
      double precision eps
      integer j
      integer k
      integer m
      integer nv
      double precision pa
      double precision pi
      double precision pmv
      double precision pss
      double precision psv
      double precision pv0
      double precision qr
      double precision r
      double precision r0
      double precision r1
      double precision r2
      double precision rg
      double precision s
      double precision s0
      double precision s1
      double precision s2
      double precision v
      double precision v0
      double precision vs
      double precision x
      double precision xq

      pi = 3.141592653589793D+00
      el = 0.5772156649015329D+00
      eps = 1.0D-14
      nv = int ( v )
      v0 = v - nv

      if ( x .eq. -1.0D+00 .and. v .ne. nv ) then
        if ( m .eq. 0 ) then
          pmv = -1.0D+300
        else
          pmv = 1.0D+300
        end if
        return
      end if

      c0 = 1.0D+00

      if ( m .ne. 0 ) then

        rg = v * ( v + m )
        do j = 1, m - 1 
          rg = rg * ( v * v - j * j )
        end do
        xq = sqrt ( 1.0D+00 - x * x )
        r0 = 1.0D+00
        do j = 1, m
          r0 = 0.5D+00 * r0 * xq / j
        end do
        c0 = r0 * rg

      end if

      if ( v0 .eq. 0.0D+00 ) then

        pmv = 1.0D+00
        r = 1.0D+00
        do k = 1, nv - m
          r = 0.5D+00 * r * ( - nv + m + k - 1.0D+00 ) 
     &      * ( nv + m + k ) / ( k * ( k + m ) ) * ( 1.0D+00 + x )
          pmv = pmv + r
        end do
        pmv = ( -1.0D+00 ) ** nv * c0 * pmv

      else

        if ( -0.35D+00 .le. x ) then

          pmv = 1.0D+00
          r = 1.0D+00
          do k = 1, 100
            r = 0.5D+00 * r * ( - v + m + k - 1.0D+00 )
     &        * ( v + m + k ) / ( k * ( m + k ) ) * ( 1.0D+00 - x )
            pmv = pmv + r
            if ( 12 .lt. k .and. abs ( r / pmv ) .lt. eps ) then
              go to 10
            end if
          end do

10        continue

          pmv = ( -1.0D+00 ) ** m * c0 * pmv

        else

          vs = sin ( v * pi ) / pi
          pv0 = 0.0D+00

          if ( m .ne. 0 ) then

            qr = sqrt ( ( 1.0D+00 - x ) / ( 1.0D+00 + x ) )
            r2 = 1.0D+00
            do j = 1, m
              r2 = r2 * qr * j
            end do
            s0 = 1.0D+00
            r1 = 1.0D+00
            do k = 1, m - 1 
              r1 = 0.5D+00 * r1 * ( - v + k - 1 ) * ( v + k ) 
     &          / ( k * ( k - m ) ) * ( 1.0D+00 + x )
              s0 = s0 + r1
            end do
            pv0 = - vs * r2 / m * s0

          end if

          call psi ( v, psv )
          pa = 2.0D+00 * ( psv + el ) + pi / tan ( pi * v ) 
     &      + 1.0D+00 / v

          s1 = 0.0D+00
          do j = 1, m
            s1 = s1 + ( j * j + v * v ) / ( j * ( j * j - v * v ) )
          end do

          pmv = pa + s1 - 1.0D+00 / ( m - v ) 
     &      + log ( 0.5D+00 * ( 1.0D+00 + x ) )
          r = 1.0D+00
          do k = 1, 100
            r = 0.5D+00 * r * ( - v + m + k - 1.0D+00 ) * ( v + m + k )
     &        / ( k * ( k + m ) ) * ( 1.0D+00 + x )
            s = 0.0D+00
            do j = 1, m
              s = s + ( ( k + j ) ** 2 + v * v ) 
     &          / ( ( k + j ) * ( ( k + j ) ** 2 - v * v ) )
            end do
            s2 = 0.0D+00
            do j = 1, k
              s2 = s2 + 1.0D+00 / ( j * ( j * j - v * v ) )
            end do
            pss = pa + s + 2.0D+00 * v * v * s2 
     &        - 1.0D+00 / ( m + k - v ) 
     &        + log ( 0.5D+00 * ( 1.0D+00 + x ) )
            r2 = pss * r
            pmv = pmv + r2
            if ( abs ( r2 / pmv ) .lt. eps ) then
              go to 20
            end if
          end do

20        continue

          pmv = pv0 + pmv * vs * c0

        end if

      end if

      return
      end
      subroutine lpn ( n, x, pn, pd )

c*********************************************************************72
c
cc LPN computes Legendre polynomials Pn(x) and derivatives Pn'(x).
c
c  Licensing:
c
c    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
c    they give permission to incorporate this routine into a user program 
c    provided that the copyright is acknowledged.
c
c  Modified:
c
c    07 July 2012
c
c  Author:
c
c    Shanjie Zhang, Jianming Jin
c
c  Reference:
c
c    Shanjie Zhang, Jianming Jin,
c    Computation of Special Functions,
c    Wiley, 1996,
c    ISBN: 0-471-11963-6,
c    LC: QA351.C45.
c
c  Parameters:
c
c    Input, integer N, the maximum degree.
c
c    Input, double precision X, the argument.
c
c    Output, double precision PN(0:N), PD(0:N), the values and derivatives
c    of the polyomials of degrees 0 to N at X.
c
      implicit none

      integer n

      integer k
      double precision p0
      double precision p1
      double precision pd(0:n)
      double precision pf
      double precision pn(0:n)
      double precision x

      pn(0) = 1.0D+00
      pn(1) = x
      pd(0) = 0.0D+00
      pd(1) = 1.0D+00
      p0 = 1.0D+00
      p1 = x

      do k = 2, n

        pf = ( 2.0D+00 * k - 1.0D+00 ) / k * x * p1 
     &    - ( k - 1.0D+00 ) / k * p0
        pn(k) = pf

        if ( abs ( x ) .eq. 1.0D+00 ) then
          pd(k) = 0.5D+00 * x ** ( k + 1 ) * k * ( k + 1.0D+00 )
        else
          pd(k) = k * ( p1 - x * pf ) / ( 1.0D+00 - x * x )
        end if

        p0 = p1
        p1 = pf

      end do

      return
      end
      subroutine lpni ( n, x, pn, pd, pl )

c*********************************************************************72
c
cc LPNI computes Legendre polynomials Pn(x), derivatives, and integrals.
c
c  Discussion:
c
c    This routine computes Legendre polynomials Pn(x), Pn'(x)
c    and the integral of Pn(t) from 0 to x.
c
c  Licensing:
c
c    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
c    they give permission to incorporate this routine into a user program 
c    provided that the copyright is acknowledged.
c
c  Modified:
c
c    13 July 2012
c
c  Author:
c
c    Shanjie Zhang, Jianming Jin
c
c  Reference:
c
c    Shanjie Zhang, Jianming Jin,
c    Computation of Special Functions,
c    Wiley, 1996,
c    ISBN: 0-471-11963-6,
c    LC: QA351.C45.
c
c  Parameters:
c
c    Input, integer N, the maximum degree.
c
c    Input, double precision X, the argument.
c
c    Output, double precision PN(0:N), PD(0:N), PL(0:N), the values, 
c    derivatives and integrals of the polyomials of degrees 0 to N at X.
c
      implicit none

      integer n

      integer j
      integer k
      integer n1
      double precision p0
      double precision p1
      double precision pd(0:n)
      double precision pf
      double precision pl(0:n)
      double precision pn(0:n)
      double precision r
      double precision x

      pn(0) = 1.0D+00
      pn(1) = x
      pd(0) = 0.0D+00
      pd(1) = 1.0D+00
      pl(0) = x
      pl(1) = 0.5D+00 * x * x
      p0 = 1.0D+00
      p1 = x

      do k = 2, n

        pf = ( 2.0D+00 * k - 1.0D+00 ) / k * x * p1 - ( k - 1.0D+00 ) 
     &    / k * p0
        pn(k) = pf

        if ( abs ( x ) .eq. 1.0D+00 ) then
          pd(k) = 0.5D+00 * x ** ( k + 1 ) * k * ( k + 1.0D+00 )
        else
          pd(k) = k * ( p1 - x * pf ) / ( 1.0D+00 - x * x )
        end if

        pl(k) = ( x * pn(k) - pn(k-1) ) / ( k + 1.0D+00 )
        p0 = p1
        p1 = pf

        if ( k .ne. 2 * int ( k / 2 ) ) then

          r = 1.0D+00 / ( k + 1.0D+00 )
          n1 = ( k - 1 ) / 2
          do j = 1, n1
            r = ( 0.5D+00 / j - 1.0D+00 ) * r
          end do
          pl(k) = pl(k) + r

        end if

      end do

      return
      end
      subroutine lqmn ( mm, m, n, x, qm, qd )

c*********************************************************************72
c
cc LQMN computes associated Legendre functions Qmn(x) and derivatives.
c
c  Discussion:
c
c    This routine computes the associated Legendre functions of the
c    second kind, Qmn(x) and Qmn'(x).
c
c  Licensing:
c
c    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
c    they give permission to incorporate this routine into a user program 
c    provided that the copyright is acknowledged.
c
c  Modified:
c
c    13 July 2012
c
c  Author:
c
c    Shanjie Zhang, Jianming Jin
c
c  Reference:
c
c    Shanjie Zhang, Jianming Jin,
c    Computation of Special Functions,
c    Wiley, 1996,
c    ISBN: 0-471-11963-6,
c    LC: QA351.C45.
c
c  Parameters:
c
c    Input, integer MM, determines the leading dimension of QM and QD.
c
c    Input, integer M, the order of Qmn(x).
c
c    Input, integer N, the degree of Qmn(x).
c
c    Output, double precision QM(0:MM,0:N), QD(0:MM,0:N), contains the values
c    of Qmn(x) and Qmn'(x).
c
      implicit none

      integer mm
      integer n

      integer i
      integer j
      integer k
      integer km
      integer ls
      integer m
      double precision q0
      double precision q1
      double precision q10
      double precision qd(0:mm,0:n)
      double precision qf
      double precision qf0
      double precision qf1
      double precision qf2
      double precision qm(0:mm,0:n)
      double precision x
      double precision xq
      double precision xs

      if ( abs ( x ) .eq. 1.0D+00 ) then
        do i = 0, m
          do j = 0, n
            qm(i,j) = 1.0D+300
            qd(i,j) = 1.0D+300
          end do
        end do
        return
      end if

      if ( 1.0D+00 .lt. abs ( x ) ) then
        ls = -1
      else
        ls = 1
      end if

      xs = ls * ( 1.0D+00 - x * x )
      xq = sqrt ( xs )
      q0 = 0.5D+00 * log ( abs ( ( x + 1.0D+00 ) / ( x - 1.0D+00 ) ) )

      if ( abs ( x ) .lt. 1.0001D+00 ) then
        qm(0,0) = q0
        qm(0,1) = x * q0 - 1.0D+00
        qm(1,0) = -1.0D+00 / xq
        qm(1,1) = -xq * ( q0 + x / ( 1.0D+00 - x * x ) )
        do i = 0, 1
          do j = 2, n
            qm(i,j) = ( ( 2.0D+00 * j - 1.0D+00 ) * x * qm(i,j-1)
     &            - ( j + i - 1.0D+00 ) * qm(i,j-2))/ ( j - i )
          end do
        end do

        do j = 0, n
          do i = 2, m
            qm(i,j) = -2.0D+00 * ( i - 1.0D+00 ) * x / xq * qm(i-1,j) 
     &        - ls * ( j + i - 1.0D+00 ) 
     &        * ( j - i + 2.0D+00 ) * qm(i-2,j)
          end do
        end do

      else

        if ( 1.1D+00 .lt. abs ( x ) ) then
          km = 40 + m + n
        else
          km = ( 40 + m + n ) 
     &      * int ( -1.0D+00 - 1.8D+00 * log ( x - 1.0D+00 ) )
        end if

        qf2 = 0.0D+00
        qf1 = 1.0D+00
        do k = km, 0, -1
          qf0 = ( ( 2 * k + 3.0D+00 ) * x * qf1 
     &      - ( k + 2.0D+00 ) * qf2 ) / ( k + 1.0D+00 )
          if ( k .le. n ) then
            qm(0,k) = qf0
          end if
          qf2 = qf1
          qf1 = qf0
        end do

        do k = 0, n
          qm(0,k) = q0 * qm(0,k) / qf0
        end do

        qf2 = 0.0D+00
        qf1 = 1.0D+00
        do k = km, 0, -1
          qf0 = ( ( 2 * k + 3.0D+00 ) * x * qf1 
     &      - ( k + 1.0D+00 ) * qf2 ) / ( k + 2.0D+00 )
          if ( k .le. n ) then
            qm(1,k) = qf0
          end if
          qf2 = qf1
          qf1 = qf0
        end do

        q10 = -1.0D+00 / xq
        do k = 0, n
          qm(1,k) = q10 * qm(1,k) / qf0
        end do

        do j = 0, n
          q0 = qm(0,j)
          q1 = qm(1,j)
          do i = 0, m - 2
            qf = -2.0D+00 * ( i + 1 ) * x / xq * q1 
     &        + ( j - i ) * ( j + i + 1.0D+00 ) * q0
            qm(i+2,j) = qf
            q0 = q1
            q1 = qf
          end do
        end do

      end if

      qd(0,0) = ls / xs
      do j = 1, n
        qd(0,j) = ls * j * ( qm(0,j-1) - x * qm(0,j) ) / xs
      end do

      do j = 0, n
        do i = 1, m
          qd(i,j) = ls * i * x / xs * qm(i,j) 
     &      + ( i + j ) * ( j - i + 1.0D+00 ) / xq * qm(i-1,j)
        end do
      end do

      return
      end
      subroutine lqmns ( m, n, x, qm, qd )

c*********************************************************************72
c
cc LQMNS computes associated Legendre functions Qmn(x) and derivatives Qmn'(x).
c
c  Licensing:
c
c    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
c    they give permission to incorporate this routine into a user program 
c    provided that the copyright is acknowledged.
c
c  Modified:
c
c    28 July 2012
c
c  Author:
c
c    Shanjie Zhang, Jianming Jin
c
c  Reference:
c
c    Shanjie Zhang, Jianming Jin,
c    Computation of Special Functions,
c    Wiley, 1996,
c    ISBN: 0-471-11963-6,
c    LC: QA351.C45.
c
c  Parameters:
c
c    Input, integer M, the order.
c
c    Input, integer N, the degree.
c
c    Input, double precision X, the argument.
c
c    Output, double precision QM(0:N), QD(0:N), the values of Qmn(x) and Qmn'(x).
c
      implicit none

      integer n

      integer k
      integer km
      integer l
      integer ls
      integer m
      double precision q0
      double precision q00
      double precision q01
      double precision q0l
      double precision q10
      double precision q11
      double precision q1l
      double precision qd(0:n)
      double precision qf0
      double precision qf1
      double precision qf2
      double precision qg0
      double precision qg1
      double precision qh0
      double precision qh1
      double precision qh2
      double precision qm(0:n)
      double precision qm0
      double precision qm1
      double precision qmk
      double precision x
      double precision xq

      do k = 0, n
        qm(k) = 0.0D+00
        qd(k) = 0.0D+00
      end do

      if ( abs ( x ) .eq. 1.0D+00 ) then
         do k = 0, n
           qm(k) = 1.0D+300
           qd(k) = 1.0D+300
         end do
         return
      end if

      if ( 1.0D+00 .lt. abs ( x ) ) then
        ls = -1
      else
        ls = +1
      end if

      xq = sqrt ( ls * ( 1.0D+00 - x * x ) )
      q0 = 0.5D+00 * log ( abs ( ( x + 1.0D+00 ) / ( x - 1.0D+00 ) ) )
      q00 = q0
      q10 = -1.0D+00 / xq
      q01 = x * q0 - 1.0D+00
      q11 = - ls * xq * ( q0 + x / ( 1.0D+00 - x * x ) )
      qf0 = q00
      qf1 = q10
      do k = 2, m
        qm0 = -2.0D+00 * ( k - 1.0D+00 ) / xq * x * qf1
     &    - ls * ( k - 1.0D+00 ) * ( 2.0D+00 - k ) * qf0
        qf0 = qf1
        qf1 = qm0
      end do

      if ( m .eq. 0 ) then
        qm0 = q00
      else if ( m .eq. 1 ) then
        qm0 = q10
      end if

      qm(0) = qm0

      if ( abs ( x ) .lt. 1.0001D+00 ) then

        if ( m .eq. 0 .and. 0 .lt. n ) then

          qf0 = q00
          qf1 = q01
          do k = 2, n
            qf2 = ( ( 2.0D+00 * k - 1.0D+00 ) * x * qf1 
     &        - ( k - 1.0D+00 ) * qf0 ) / k
            qm(k) = qf2
            qf0 = qf1
            qf1 = qf2
          end do

        end if
        qg0 = q01
        qg1 = q11
        do k = 2, m
          qm1 = - 2.0D+00 * ( k - 1.0D+00 ) / xq * x * qg1 
     &      - ls * k * ( 3.0D+00 - k ) * qg0
          qg0 = qg1
          qg1 = qm1
        end do

        if ( m .eq. 0 ) then
          qm1 = q01
        else if ( m .eq. 1 ) then
          qm1 = q11
        end if
        qm(1) = qm1

        if ( m .eq. 1 .and. 1 .lt. n ) then

          qh0 = q10
          qh1 = q11
          do k = 2, n
            qh2 = ( ( 2.0D+00 * k - 1.0D+00 ) * x * qh1 - k * qh0 ) 
     &        / ( k - 1.0D+00 )
            qm(k) = qh2
            qh0 = qh1
            qh1 = qh2
          end do

        else if ( 2 .le. m ) then

          qg0 = q00
          qg1 = q01
          qh0 = q10
          qh1 = q11

          do l = 2, n
            q0l = ( ( 2.0D+00 * l - 1.0D+00 ) * x * qg1 
     &        - ( l - 1.0D+00 ) * qg0 ) / l
            q1l = ( ( 2.0D+00 * l - 1.0D+00 ) * x * qh1 - l * qh0 ) 
     &        / ( l - 1.0D+00 )
            qf0 = q0l
            qf1 = q1l
            do k = 2, m
              qmk = - 2.0D+00 * ( k - 1.0D+00 ) / xq * x * qf1 
     &          - ls * ( k + l - 1.0D+00 ) * ( l + 2.0D+00 - k ) * qf0
              qf0 = qf1
              qf1 = qmk
            end do
            qm(l) = qmk
            qg0 = qg1
            qg1 = q0l
            qh0 = qh1
            qh1 = q1l
          end do

        end if

      else

        if ( 1.1D+00 .lt. abs ( x ) ) then
          km = 40 + m + n
        else
          km = ( 40 + m + n ) 
     &      * int ( - 1.0D+00 - 1.8D+00 * log ( x - 1.0D+00 ) )
        end if

        qf2 = 0.0D+00
        qf1 = 1.0D+00
        do k = km, 0, -1
          qf0 = ( ( 2.0D+00 * k + 3.0D+00 ) * x * qf1 
     &      - ( k + 2.0D+00 - m ) * qf2 ) / ( k + m + 1.0D+00 )
          if ( k .le. n ) then
            qm(k) = qf0
          end if
          qf2 = qf1
          qf1 = qf0
        end do

        do k = 0, n
         qm(k) = qm(k) * qm0 / qf0
        end do

      end if

      if ( abs ( x ) .lt. 1.0D+00 ) then
        do k = 0, n
          qm(k) = ( -1 ) ** m * qm(k)
        end do
      end if

      qd(0) = ( ( 1.0D+00 - m ) * qm(1) - x * qm(0) ) 
     &  / ( x * x - 1.0D+00 )
      do k = 1, n
        qd(k) = ( k * x * qm(k) - ( k + m ) * qm(k-1) ) 
     &    / ( x * x - 1.0D+00 )
      end do

      return
      end
      subroutine lqna ( n, x, qn, qd )

c*********************************************************************72
c
cc LQNA computes Legendre function Qn(x) and derivatives Qn'(x).
c
c  Licensing:
c
c    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
c    they give permission to incorporate this routine into a user program 
c    provided that the copyright is acknowledged.
c
c  Modified:
c
c    19 July 2012
c
c  Author:
c
c    Shanjie Zhang, Jianming Jin
c
c  Reference:
c
c    Shanjie Zhang, Jianming Jin,
c    Computation of Special Functions,
c    Wiley, 1996,
c    ISBN: 0-471-11963-6,
c    LC: QA351.C45.
c
c  Parameters:
c
c    Input, integer N, the degree of Qn(x).
c
c    Input, double precision X, the argument of Qn(x).
c
c    Output, double precision QN(0:N), QD(0:N), the values of
c    Qn(x) and Qn'(x).
c
      implicit none

      integer n

      integer k
      double precision q0
      double precision q1
      double precision qd(0:n)
      double precision qf
      double precision qn(0:n)
      double precision x

      if ( abs ( x ) .eq. 1.0D+00 ) then

        do k = 0, n
          qn(k) = 1.0D+300
          qd(k) = -1.0D+300
        end do

      else if ( abs ( x ) .lt. 1.0D+00 ) then

        q0 = 0.5D+00 * log ( ( 1.0D+00 + x ) / ( 1.0D+00 - x ) )
        q1 = x * q0 - 1.0D+00
        qn(0) = q0
        qn(1) = q1
        qd(0) = 1.0D+00 / ( 1.0D+00 - x * x )
        qd(1) = qn(0) + x * qd(0)
        do k = 2, n
          qf = ( ( 2 * k - 1 ) * x * q1 - ( k - 1 ) * q0 ) / k
          qn(k) = qf
          qd(k) = ( qn(k-1) - x * qf ) * k / ( 1.0D+00 - x * x )
          q0 = q1
          q1 = qf
        end do

      end if

      return
      end
      subroutine lqnb ( n, x, qn, qd )

c*********************************************************************72
c
cc LQNB computes Legendre function Qn(x) and derivatives Qn'(x).
c
c  Licensing:
c
c    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
c    they give permission to incorporate this routine into a user program 
c    provided that the copyright is acknowledged.
c
c  Modified:
c
c    19 July 2012
c
c  Author:
c
c    Shanjie Zhang, Jianming Jin
c
c  Reference:
c
c    Shanjie Zhang, Jianming Jin,
c    Computation of Special Functions,
c    Wiley, 1996,
c    ISBN: 0-471-11963-6,
c    LC: QA351.C45.
c
c  Parameters:
c
c    Input, integer N, the degree of Qn(x).
c
c    Input, double precision X, the argument of Qn(x).
c
c    Output, double precision QN(0:N), QD(0:N), the values of
c    Qn(x) and Qn'(x).
c
      implicit none

      integer n

      double precision eps
      integer j
      integer k
      integer l
      integer nl
      double precision q0
      double precision q1
      double precision qc1
      double precision qc2
      double precision qd(0:n)
      double precision qf
      double precision qf0
      double precision qf1
      double precision qf2
      double precision qn(0:n)
      double precision qr
      double precision x
      double precision x2

      eps = 1.0D-14

      if ( abs ( x ) .eq. 1.0D+00 ) then
        do k = 0, n
          qn(k) = 1.0D+300
          qd(k) = 1.0D+300
        end do
        return
      end if

      if ( x .le. 1.021D+00 ) then

        x2 = abs ( ( 1.0D+00 + x ) / ( 1.0D+00 - x ) )
        q0 = 0.5D+00 * log ( x2 )
        q1 = x * q0 - 1.0D+00
        qn(0) = q0
        qn(1) = q1
        qd(0) = 1.0D+00 / ( 1.0D+00 - x * x )
        qd(1) = qn(0) + x * qd(0)
        do k = 2, n
          qf = ( ( 2.0D+00 * k - 1.0D+00 ) * x * q1 
     &      - ( k - 1.0D+00 ) * q0 ) / k
          qn(k) = qf
          qd(k) = ( qn(k-1) - x * qf ) * k / ( 1.0D+00 - x * x )
          q0 = q1
          q1 = qf
        end do

      else

        qc2 = 1.0D+00 / x
        do j = 1, n
          qc2 = qc2 * j / ( ( 2.0D+00 * j + 1.0D+00 ) * x )
          if ( j .eq. n - 1 ) then
            qc1 = qc2
          end if
        end do

        do l = 0, 1

          nl = n + l
          qf = 1.0D+00
          qr = 1.0D+00
          do k = 1, 500
            qr = qr * ( 0.5D+00 * nl + k - 1.0D+00 ) 
     &        * ( 0.5D+00 * ( nl - 1 ) + k )
     &        / ( ( nl + k - 0.5D+00 ) * k * x * x )
            qf = qf + qr
            if ( abs ( qr / qf ) .lt. eps ) then
              go to 10
            end if
          end do

10        continue

          if ( l .eq. 0 ) then
            qn(n-1) = qf * qc1
          else
            qn(n) = qf * qc2
          end if

        end do

        qf2 = qn(n)
        qf1 = qn(n-1)
        do k = n, 2, -1
          qf0 = ( ( 2.0D+00 * k - 1.0D+00 ) * x * qf1 - k * qf2 ) 
     &      / ( k - 1.0D+00 )
          qn(k-2) = qf0
          qf2 = qf1
          qf1 = qf0
        end do

        qd(0) = 1.0D+00 / ( 1.0D+00 - x * x )
        do k = 1, n
          qd(k) = k * ( qn(k-1) - x * qn(k) ) / ( 1.0D+00 - x * x )
        end do

      end if

      return
      end
      function msta1 ( x, mp )

c*********************************************************************72
c
cc MSTA1 determines a backward recurrence starting point for Jn(x).
c
c  Discussion:
c
c    This procedure determines the starting point for backward  
c    recurrence such that the magnitude of    
c    Jn(x) at that point is about 10^(-MP).
c
c  Licensing:
c
c    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
c    they give permission to incorporate this routine into a user program 
c    provided that the copyright is acknowledged.
c
c  Modified:
c
c    08 July 2012
c
c  Author:
c
c    Shanjie Zhang, Jianming Jin
c
c  Reference:
c
c    Shanjie Zhang, Jianming Jin,
c    Computation of Special Functions,
c    Wiley, 1996,
c    ISBN: 0-471-11963-6,
c    LC: QA351.C45.
c
c  Parameters:
c
c    Input, double precision X, the argument.
c
c    Input, integer MP, the negative logarithm of the desired magnitude.
c
c    Output, integer MSTA1, the starting point.
c
      implicit none

      double precision a0
      double precision envj
      double precision f
      double precision f0
      double precision f1
      integer it
      integer mp
      integer msta1
      integer n0
      integer n1
      integer nn
      double precision x

      a0 = abs ( x )
      n0 = int ( 1.1D+00 * a0 ) + 1
      f0 = envj ( n0, a0 ) - mp
      n1 = n0 + 5
      f1 = envj ( n1, a0 ) - mp
      do it = 1, 20       
        nn = n1 - ( n1 - n0 ) / ( 1.0D+00 - f0 / f1 )                  
        f = envj ( nn, a0 ) - mp
        if ( abs ( nn - n1 ) .lt. 1 ) then
          go to 10
        end if
        n0 = n1
        f0 = f1
        n1 = nn
        f1 = f
      end do

 10   continue

      msta1 = nn

      return
      end
      function msta2 ( x, n, mp )

c*********************************************************************72
c
cc MSTA2 determines a backward recurrence starting point for Jn(x).
c
c  Discussion:
c
c    This procedure determines the starting point for a backward
c    recurrence such that all Jn(x) has MP significant digits.
c
c  Licensing:
c
c    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
c    they give permission to incorporate this routine into a user program 
c    provided that the copyright is acknowledged.
c
c  Modified:
c
c    08 July 2012
c
c  Author:
c
c    Shanjie Zhang, Jianming Jin
c
c  Reference:
c
c    Shanjie Zhang, Jianming Jin,
c    Computation of Special Functions,
c    Wiley, 1996,
c    ISBN: 0-471-11963-6,
c    LC: QA351.C45.
c
c  Parameters:
c
c    Input, double precision X, the argument of Jn(x).
c
c    Input, integer N, the order of Jn(x).
c
c    Input, integer MP, the number of significant digits.
c
c    Output, integer MSTA2, the starting point.
c
      implicit none

      double precision a0
      double precision ejn
      double precision envj
      double precision f
      double precision f0
      double precision f1
      double precision hmp
      integer it
      integer mp
      integer msta2
      integer n
      integer n0
      integer n1
      integer nn
      double precision obj
      double precision x

      a0 = abs ( x )
      hmp = 0.5D+00 * mp
      ejn = envj ( n, a0 )

      if ( ejn .le. hmp ) then
        obj = mp
        n0 = int ( 1.1D+00 * a0 )
      else
        obj = hmp + ejn
        n0 = n
      end if

      f0 = envj ( n0, a0 ) - obj
      n1 = n0 + 5
      f1 = envj ( n1, a0 ) - obj

      do it = 1, 20
        nn = n1 - ( n1 - n0 ) / ( 1.0D+00 - f0 / f1 )
        f = envj ( nn, a0 ) - obj
        if ( abs ( nn - n1 ) .lt. 1 ) then
          go to 10
        end if
        n0 = n1
        f0 = f1
        n1 = nn
        f1 = f
      end do

10    continue

      msta2 = nn + 10

      return
      end
      subroutine mtu0 ( kf, m, q, x, csf, csd )

c*********************************************************************72
c
cc MTU0 computes Mathieu functions CEM(x,q) and SEM(x,q) and derivatives.
c
c  Licensing:
c
c    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
c    they give permission to incorporate this routine into a user program 
c    provided that the copyright is acknowledged.
c
c  Modified:
c
c    20 July 2012
c
c  Author:
c
c    Shanjie Zhang, Jianming Jin
c
c  Reference:
c
c    Shanjie Zhang, Jianming Jin,
c    Computation of Special Functions,
c    Wiley, 1996,
c    ISBN: 0-471-11963-6,
c    LC: QA351.C45.
c
c  Parameters:
c
c    Input, integer KF, the function code.
c    1 for computing cem(x,q) and cem'(x,q)
c    2 for computing sem(x,q) and sem'(x,q).
c
c    Input, integer M, the order of the Mathieu functions.
c
c    Input, double precision Q, the parameter of the Mathieu functions.
c
c    Input, double precision X, the argument of the Mathieu functions in degrees.
c
c    Output, double precision CSF, CSD, the values of cem(x,q) and cem'(x,q),
c    or of sem(x,q) and sem'(x,q).
c
      implicit none

      double precision a
      double precision csd
      double precision csf
      double precision eps
      double precision fg(251)
      integer ic
      integer k
      integer kd
      integer kf
      integer km
      integer m
      double precision q
      double precision qm
      double precision rd
      double precision x
      double precision xr

      eps = 1.0D-14

      if ( kf .eq. 1 ) then

        if ( m .eq. 2 * int ( m / 2 ) ) then
          kd = 1
        else
          kd = 2
        end if

      else

        if ( m .ne. 2 * int ( m / 2 ) ) then
          kd = 3
        else
          kd = 4
        end if

      end if

      call cva2 ( kd, m, q, a )

      if ( q .le. 1.0D+00 ) then
        qm = 7.5D+00 + 56.1D+00 * sqrt ( q ) - 134.7D+00 * q 
     &    + 90.7D+00 * sqrt ( q ) * q
      else
        qm = 17.0D+00 + 3.1D+00 * sqrt ( q ) - 0.126D+00 * q 
     &    + 0.0037D+00 * sqrt ( q ) * q
      end if

      km = int ( qm + 0.5D+00 * m )
      call fcoef ( kd, m, q, a, fg )
      ic = int ( m / 2 ) + 1
      rd = 1.74532925199433D-02
      xr = x * rd

      csf = 0.0D+00

      do k = 1, km

        if ( kd .eq. 1 ) then
          csf = csf + fg(k) * cos ( ( 2.0D+00 * k - 2.0D+00 ) * xr )
        else if ( kd .eq. 2 ) then
          csf = csf + fg(k) * cos ( ( 2.0D+00 * k - 1.0D+00 ) * xr )
        else if ( kd .eq. 3 ) then
          csf = csf + fg(k) * sin ( ( 2.0D+00 * k - 1.0D+00 ) * xr )
        else if ( kd .eq. 4 ) then
          csf = csf + fg(k) * sin ( 2.0D+00 * k * xr )
        end if

        if ( ic .le. k .and. abs ( fg(k) ) .lt. abs ( csf ) * eps ) then
          go to 10
        end if

      end do

10    continue

      csd = 0.0D+00

      do k = 1, km

        if ( kd .eq. 1 ) then
          csd = csd - ( 2 * k - 2 ) * fg(k) * sin ( ( 2 * k - 2 ) * xr )
        else if ( kd .eq. 2 ) then
          csd = csd - ( 2 * k - 1 ) * fg(k) * sin ( ( 2 * k - 1 ) * xr )
        else if ( kd .eq. 3 ) then
          csd = csd + ( 2 * k - 1 ) * fg(k) * cos ( ( 2 * k - 1 ) * xr )
        else if ( kd .eq. 4 ) then
          csd = csd + 2.0D+00 * k * fg(k) * cos ( 2 * k * xr )
        end if

        if ( ic .le. k .and. abs ( fg(k) ) .lt. abs ( csd ) * eps ) then
          go to 20
        end if

      end do

20    continue

      return
      end
      subroutine mtu12 ( kf, kc, m, q, x, f1r, d1r, f2r, d2r )

c*********************************************************************72
c
cc MTU12 computes modified Mathieu functions of the first and second kind.
c
c  Discussion:
c
c    This procedure computes modified Mathieu functions of the first and
c    second kinds, Mcm(1)(2)(x,q) and Msm(1)(2)(x,q),
c    and their derivatives.
c
c  Licensing:
c
c    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
c    they give permission to incorporate this routine into a user program 
c    provided that the copyright is acknowledged.
c
c  Modified:
c
c    10 April 2012
c
c  Author:
c
c    Shanjie Zhang, Jianming Jin
c
c  Reference:
c
c    Shanjie Zhang, Jianming Jin,
c    Computation of Special Functions,
c    Wiley, 1996,
c    ISBN: 0-471-11963-6,
c    LC: QA351.C45.
c
c  Parameters:
c
c    Input, integer KF, the function code.
c    1 for computing Mcm(x,q);
c    2 for computing Msm(x,q).
c
c    Input, integer KC, the function code.
c    1, for computing the first kind
c    2, for computing the second kind or Msm(2)(x,q) and Msm(2)'(x,q)
c    3, for computing both the first and second kinds.
c
c    Input, integer M, the order of the Mathieu functions.
c
c    Input, double precision Q, the parameter of the Mathieu functions.
c
c    Input, double precision X, the argument of the Mathieu functions.
c
c    Output, double precision F1R, D1R, F2R, D2R, the values of 
c    Mcm(1)(x,q) or Msm(1)(x,q), Derivative of Mcm(1)(x,q) or Msm(1)(x,q),
c    Mcm(2)(x,q) or Msm(2)(x,q), Derivative of Mcm(2)(x,q) or Msm(2)(x,q).
c
      implicit none

      double precision a
      double precision bj1(0:251)
      double precision bj2(0:251)
      double precision by1(0:251)
      double precision by2(0:251)
      double precision c1
      double precision c2
      double precision d1r
      double precision d2r
      double precision dj1(0:251)
      double precision dj2(0:251)
      double precision dy1(0:251)
      double precision dy2(0:251)
      double precision eps
      double precision f1r
      double precision f2r
      double precision fg(251)
      integer ic
      integer k
      integer kc
      integer kd
      integer kf
      integer km
      integer m
      integer nm
      double precision q
      double precision qm
      double precision u1
      double precision u2
      double precision w1
      double precision w2
      double precision x

      eps = 1.0D-14

      if ( kf .eq. 1 ) then
        if ( m .eq. 2 * int ( m / 2 ) ) then
          kd = 1
        else
          kd = 2
        end if
      else
        if ( m .ne. 2 * int ( m / 2 ) ) then
          kd = 3
        else
          kd = 4
        end if
      end if

      call cva2 ( kd, m, q, a )

      if ( q .le. 1.0D+00 ) then
        qm = 7.5D+00 + 56.1D+00 * sqrt ( q ) - 134.7D+00 * q 
     &    + 90.7D+00 * sqrt ( q ) * q
      else
        qm = 17.0D+00 + 3.1D+00 * sqrt ( q ) - 0.126D+00 * q 
     &    + 0.0037D+00 * sqrt ( q ) * q
      end if

      km = int ( qm + 0.5D+00 * m )              
      call fcoef ( kd, m, q, a, fg )

      if ( kd .eq. 4 ) then
        ic = m / 2
      else
        ic = int ( m / 2 ) + 1
      end if

      c1 = exp ( - x )
      c2 = exp ( x )
      u1 = sqrt ( q ) * c1
      u2 = sqrt ( q ) * c2

      call jynb ( km, u1, nm, bj1, dj1, by1, dy1 )
      call jynb ( km, u2, nm, bj2, dj2, by2, dy2 )

      if ( kc .eq. 2 ) then
        go to 30
      end if

      f1r = 0.0D+00

      do k = 1, km

        if ( kd .eq. 1 ) then
          f1r = f1r + ( - 1.0D+00 ) ** ( ic + k ) * fg(k) 
     &      * bj1(k-1) * bj2(k-1)
        else if ( kd .eq. 2 .or. kd .eq. 3 ) then
          f1r = f1r + ( -1.0D+00 ) ** ( ic + k ) * fg(k) 
     &      * ( bj1(k-1) * bj2(k) 
     &      + ( - 1.0D+00 ) ** kd * bj1(k) * bj2(k-1) )
        else
          f1r = f1r + ( -1.0D+00 ) ** ( ic + k ) * fg(k)
     &      * ( bj1(k-1) * bj2(k+1) - bj1(k+1) * bj2(k-1) )
        end if

        if ( 5 .le. k .and. 
     &    abs ( f1r - w1 ) .lt. abs ( f1r ) * eps ) then
          go to 10
        end if

        w1 = f1r

      end do

10    continue

      f1r = f1r / fg(1)
      d1r = 0.0D+00
      do k = 1, km
        if ( kd .eq. 1 ) then
          d1r = d1r + ( - 1.0D+00 ) ** ( ic + k ) * fg(k) 
     &      * ( c2 * bj1(k-1) * dj2(k-1) - c1 * dj1(k-1) * bj2(k-1) )
        else if ( kd .eq. 2 .or. kd .eq. 3 ) then
          d1r = d1r + ( -1.0D+00 ) ** ( ic + k ) * fg(k) 
     &      * ( c2 * ( bj1(k-1) * dj2(k)
     &      + ( -1.0D+00 ) ** kd * bj1(k) * dj2(k-1) ) 
     &      - c1 * ( dj1(k-1) * bj2(k)
     &      + ( -1.0D+00 ) ** kd * dj1(k) * bj2(k-1) ) )
        else
          d1r = d1r + ( -1.0D+00 ) ** ( ic + k ) * fg(k) 
     &      * ( c2 * ( bj1(k-1) * dj2(k+1) - bj1(k+1) * dj2(k-1) ) 
     &      - c1 * ( dj1(k-1) * bj2(k+1) - dj1(k+1) * bj2(k-1) ) )
        end if
        if ( 5 .le. k .and. 
     &    abs ( d1r - w2 ) .lt. abs ( d1r ) * eps ) then
          go to 20
        end if
        w2 = d1r
      end do

20    continue

      d1r = d1r * sqrt ( q ) / fg(1)

      if ( kc .eq. 1 ) then
        return
      end if

30    continue

      f2r = 0.0D+00

      do k = 1, km
        if ( kd .eq. 1 ) then
          f2r = f2r + ( -1.0D+00 ) ** ( ic + k ) * fg(k) 
     &      * bj1(k-1) * by2(k-1)
        else if ( kd .eq. 2 .or. kd .eq. 3 ) then
          f2r = f2r + ( -1.0D+00 ) ** ( ic + k ) * fg(k)
     &      * ( bj1(k-1) * by2(k) 
     &      + ( -1.0D+00 ) ** kd * bj1(k) * by2(k-1) )
        else
          f2r = f2r + ( -1.0D+00 ) ** ( ic + k ) * fg(k)
     &      * ( bj1(k-1) * by2(k+1) - bj1(k+1) * by2(k-1) )
        end if
        if ( 5 .le. k .and. 
     &    abs ( f2r - w1 ) .lt. abs ( f2r ) * eps ) then
          go to 40
        end if
        w1 = f2r
      end do

40    continue

      f2r = f2r / fg(1)
      d2r = 0.0D+00

      do k = 1, km
        if ( kd .eq. 1 ) then
          d2r = d2r + ( -1.0D+00 ) ** ( ic + k ) * fg(k) 
     &      * ( c2 * bj1(k-1) * dy2(k-1) - c1 * dj1(k-1) * by2(k-1) )
        else if ( kd .eq. 2 .or. kd .eq. 3 ) then
          d2r = d2r + ( -1.0D+00 ) ** ( ic + k ) * fg(k)
     &      * ( c2 * ( bj1(k-1) * dy2(k)
     &      + ( -1.0D+00 ) ** kd * bj1(k) * dy2(k-1) )
     &      - c1 * ( dj1(k-1) * by2(k)
     &      + ( -1.0D+00 ) ** kd * dj1(k) * by2(k-1) ) )
        else
          d2r = d2r + ( -1.0D+00 ) ** ( ic + k ) * fg(k)
     &      * ( c2 * ( bj1(k-1) * dy2(k+1) - bj1(k+1) * dy2(k-1) )
     &      - c1 * ( dj1(k-1) * by2(k+1) - dj1(k+1) * by2(k-1) ) )
        end if

        if ( 5 .le. k .and. 
     &    abs ( d2r - w2 ) .lt. abs ( d2r ) * eps ) then
          go to 50
        end if

        w2 = d2r

      end do

50    continue

      d2r = d2r * sqrt ( q ) / fg(1)

      return
      end
      subroutine othpl ( kf, n, x, pl, dpl )

c*********************************************************************72
c
cc OTHPL computes orthogonal polynomials Tn(x), Un(x), Ln(x) or Hn(x).
c
c  Discussion:
c
c    This procedure computes orthogonal polynomials: Tn(x) or Un(x),
c    or Ln(x) or Hn(x), and their derivatives.
c
c  Licensing:
c
c    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
c    they give permission to incorporate this routine into a user program 
c    provided that the copyright is acknowledged.
c
c  Modified:
c
c    08 July 2012
c
c  Author:
c
c    Shanjie Zhang, Jianming Jin
c
c  Reference:
c
c    Shanjie Zhang, Jianming Jin,
c    Computation of Special Functions,
c    Wiley, 1996,
c    ISBN: 0-471-11963-6,
c    LC: QA351.C45.
c
c  Parameters:
c
c    Input, integer KT, the function code:
c    1 for Chebyshev polynomial Tn(x)
c    2 for Chebyshev polynomial Un(x)
c    3 for Laguerre polynomial Ln(x)
c    4 for Hermite polynomial Hn(x)
c
c    Input, integer N, the order.
c
c    Input, double precision X, the argument.
c
c    Output, double precision PL(0:N), DPL(0:N), the value and derivative of
c    the polynomials of order 0 through N at X.
c
      implicit none

      integer n

      double precision a
      double precision b
      double precision c
      double precision dpl(0:n)
      double precision dy0
      double precision dy1
      double precision dyn
      integer k
      integer kf
      double precision pl(0:n)
      double precision x
      double precision y0
      double precision y1
      double precision yn

      a = 2.0D+00
      b = 0.0D+00
      c = 1.0D+00
      y0 = 1.0D+00
      y1 = 2.0D+00 * x
      dy0 = 0.0D+00
      dy1 = 2.0D+00
      pl(0) = 1.0D+00
      pl(1) = 2.0D+00 * x
      dpl(0) = 0.0D+00
      dpl(1) = 2.0D+00

      if ( kf .eq. 1 ) then
        y1 = x
        dy1 = 1.0D+00
        pl(1) = x
        dpl(1) = 1.0D+00
      else if ( kf .eq. 3 ) then
        y1 = 1.0D+00 - x
        dy1 = -1.0D+00
        pl(1) = 1.0D+00 - x
        dpl(1) = -1.0D+00
      end if

      do k = 2, n

        if ( kf .eq. 3 ) then
          a = -1.0D+00 / k
          b = 2.0D+00 + a
          c = 1.0D+00 + a
        else if ( kf .eq. 4 ) then
          c = 2.0D+00 * ( k - 1.0D+00 )
        end if

        yn = ( a * x + b ) * y1 - c * y0
        dyn = a * y1 + ( a * x + b ) * dy1 - c * dy0
        pl(k) = yn
        dpl(k) = dyn
        y0 = y1
        y1 = yn
        dy0 = dy1
        dy1 = dyn

      end do

      return
      end
      subroutine pbdv ( v, x, dv, dp, pdf, pdd )

c*********************************************************************72
c
cc PBDV computes parabolic cylinder functions Dv(x) and derivatives.
c
c  Licensing:
c
c    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
c    they give permission to incorporate this routine into a user program 
c    provided that the copyright is acknowledged.
c
c  Modified:
c
c    29 July 2012
c
c  Author:
c
c    Shanjie Zhang, Jianming Jin
c
c  Reference:
c
c    Shanjie Zhang, Jianming Jin,
c    Computation of Special Functions,
c    Wiley, 1996,
c    ISBN: 0-471-11963-6,
c    LC: QA351.C45.
c
c  Parameters:
c
c    Input, double precision V, the order.
c
c    Input, double precision X, the argument.
c
c    Output, double precision DV(0:*), DP(0:*), the values of
c    Dn+v0(x), Dn+v0'(x).
c
c    Output, double precision PDF, PDD, the values of Dv(x) and Dv'(x).
c
      implicit none

      double precision dp(0:*)
      double precision dv(0:*)
      double precision ep
      double precision f
      double precision f0
      double precision f1
      integer ja
      integer k
      integer l
      integer m
      integer na
      integer nk
      integer nv
      double precision pd
      double precision pd0
      double precision pd1
      double precision pdd
      double precision pdf
      double precision s0
      double precision v
      double precision v0
      double precision v1
      double precision v2
      double precision vh
      double precision x
      double precision xa

      xa = abs ( x )
      vh = v
      v = v + sign ( 1.0D+00, v )
      nv = int ( v )
      v0 = v - nv
      na = abs ( nv )
      ep = exp ( -0.25D+00 * x * x )

      if ( 1 .le. na ) then
        ja = 1
      end if

      if ( 0.0D+00 .le. v ) then
        if ( v0 .eq. 0.0D+00 ) then
          pd0 = ep
          pd1 = x * ep
        else
          do l = 0, ja
            v1 = v0 + l
            if ( xa .le. 5.8D+00 ) then
              call dvsa ( v1, x, pd1 )
            else
              call dvla ( v1, x, pd1 )
            end if
            if ( l .eq. 0 ) then
              pd0 = pd1
            end if
          end do
        end if

        dv(0) = pd0
        dv(1) = pd1
        do k = 2, na
          pdf = x * pd1 - ( k + v0 - 1.0D+00 ) * pd0
          dv(k) = pdf
          pd0 = pd1
          pd1 = pdf
        end do

      else

        if ( x .le. 0.0D+00 ) then

          if ( xa .le. 5.8D+00 )  then
            call dvsa ( v0, x, pd0 )
            v1 = v0 - 1.0D+00
            call dvsa ( v1, x, pd1 )
          else
            call dvla ( v0, x, pd0 )
            v1 = v0 - 1.0D+00
            call dvla ( v1, x, pd1 )
          end if

          dv(0) = pd0
          dv(1) = pd1
          do k = 2, na
            pd = ( - x * pd1 + pd0 ) / ( k - 1.0D+00 - v0 )
            dv(k) = pd
            pd0 = pd1
            pd1 = pd
          end do

        else if ( x .le. 2.0D+00 ) then

          v2 = nv + v0
          if ( nv .eq. 0 ) then
            v2 = v2 - 1.0D+00
          end if

          nk = int ( - v2 )
          call dvsa ( v2, x, f1 )
          v1 = v2 + 1.0D+00
          call dvsa ( v1, x, f0 )
          dv(nk) = f1
          dv(nk-1) = f0
          do k = nk - 2, 0, -1
            f = x * f0 + ( k - v0 + 1.0D+00 ) * f1
            dv(k) = f
            f1 = f0
            f0 = f
          end do

        else

          if ( xa .le. 5.8D+00 ) then
            call dvsa ( v0, x, pd0 )
          else
            call dvla ( v0, x, pd0 )
          end if

          dv(0) = pd0
          m = 100 + na
          f1 = 0.0D+00
          f0 = 1.0D-30
          do k = m, 0, -1
            f = x * f0 + ( k - v0 + 1.0D+00 ) * f1
            if ( k .le. na ) then
              dv(k) = f
            end if
            f1 = f0
            f0 = f
          end do
          s0 = pd0 / f
          do k = 0, na
            dv(k) = s0 * dv(k)
          end do

        end if

      end if

      do k = 0, na - 1
        v1 = abs ( v0 ) + k
        if ( 0.0D+00 .le. v ) then
          dp(k) = 0.5D+00 * x * dv(k) - dv(k+1)
        else
          dp(k) = -0.5D+00 * x * dv(k) - v1 * dv(k+1)
        end if
      end do

      pdf = dv(na-1)
      pdd = dp(na-1)
      v = vh

      return
      end
      subroutine pbvv ( v, x, vv, vp, pvf, pvd )

c*********************************************************************72
c
cc PBVV computes parabolic cylinder functions Vv(x) and their derivatives.
c
c  Licensing:
c
c    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
c    they give permission to incorporate this routine into a user program 
c    provided that the copyright is acknowledged.
c
c  Modified:
c
c    29 July 2012
c
c  Author:
c
c    Shanjie Zhang, Jianming Jin
c
c  Reference:
c
c    Shanjie Zhang, Jianming Jin,
c    Computation of Special Functions,
c    Wiley, 1996,
c    ISBN: 0-471-11963-6,
c    LC: QA351.C45.
c
c  Parameters:
c
c    Input, double precision V, the order.
c
c    Input, double precision X, the argument.
c
c    Output, double precision VV(0:*), VP(0:*), the values of Vv(x), Vv'(x).
c7
c    Output, double precision PVF, PVD, the values of Vv(x) and Vv'(x).
c
      implicit none

      double precision f
      double precision f0
      double precision f1
      integer ja
      integer k
      integer kv
      integer l
      integer m
      integer na
      integer nv
      double precision pi
      double precision pv0
      double precision pvd
      double precision pvf
      double precision q2p
      double precision qe
      double precision s0
      double precision v
      double precision v0
      double precision v1
      double precision v2
      double precision vh
      double precision vp(0:*)
      double precision vv(0:*)
      double precision x
      double precision xa

      pi = 3.141592653589793D+00
      xa = abs ( x )
      vh = v
      v = v + sign ( 1.0D+00, v )
      nv = int ( v )
      v0 = v - nv
      na = abs ( nv )
      qe = exp ( 0.25D+00 * x * x )
      q2p = sqrt ( 2.0D+00 / pi )

      if ( 1 .le. na ) then
        ja = 1
      end if

      if ( v .le. 0.0D+00 ) then

        if ( v0 .eq. 0.0D+00 ) then

          if ( xa .le. 7.5D+00 ) then 
            call vvsa ( v0, x, pv0 )
          else
            call vvla ( v0, x, pv0 )
          end if

          f0 = q2p * qe
          f1 = x * f0
          vv(0) = pv0
          vv(1) = f0
          vv(2) = f1

        else

          do l = 0, ja
            v1 = v0 - l
            if ( xa .le. 7.5D+00 ) then
              call vvsa ( v1, x, f1 )
            else
              call vvla ( v1, x, f1 )
            end if
            if ( l .eq. 0 ) then
              f0 = f1
            end if
          end do

          vv(0) = f0
          vv(1) = f1

        end if

        if ( v0 .eq. 0.0D+00 ) then
          kv = 3
        else
          kv = 2
        end if

        do k = kv, na
          f = x * f1 + ( k - v0 - 2.0D+00 ) * f0
          vv(k) = f
          f0 = f1
          f1 = f
        end do

      else

        if ( 0.0D+00 .le. x .and. x .le. 7.5D+00 ) then

          v2 = v
          if ( v2 .lt. 1.0D+00 ) then
            v2 = v2 + 1.0D+00
          end if

          call vvsa ( v2, x, f1 )
          v1 = v2 - 1.0D+00
          kv = int ( v2 )
          call vvsa ( v1, x, f0 )
          vv(kv) = f1
          vv(kv-1) = f0
          do k = kv - 2, 0, - 1
            f = x * f0 - ( k + v0 + 2.0D+00 ) * f1
            if ( k .le. na ) then
              vv(k) = f
            end if
            f1 = f0
            f0 = f
          end do

        else if ( 7.5D+00 .lt. x ) then

          call vvla ( v0, x, pv0 )
          m = 100 + abs ( na )
          vv(1) = pv0
          f1 = 0.0D+00
          f0 = 1.0D-40
          do k = m, 0, -1
            f = x * f0 - ( k + v0 + 2.0D+00 ) * f1
            if ( k .le. na ) then
              vv(k) = f
            end if
            f1 = f0
            f0 = f
          end do
          s0 = pv0 / f
          do k = 0, na
            vv(k) = s0 * vv(k)
          end do

        else

          if ( xa .le. 7.5D+00 ) then
            call vvsa ( v0, x, f0 )
            v1 = v0 + 1.0D+00
            call vvsa ( v1, x, f1 )
          else
            call vvla ( v0, x, f0 )
            v1 = v0 + 1.0D+00
            call vvla ( v1, x, f1 )
          end if

          vv(0) = f0
          vv(1) = f1
          do k = 2, na
            f = ( x * f1 - f0 ) / ( k + v0 )
            vv(k) = f
            f0 = f1
            f1 = f
          end do

        end if

      end if

      do k = 0, na - 1
        v1 = v0 + k
        if ( 0.0D+00 .le. v ) then
          vp(k) = 0.5D+00 * x * vv(k) - ( v1 + 1.0D+00 ) * vv(k+1)
        else
          vp(k) = - 0.5D+00 * x * vv(k) + vv(k+1)
        end if
      end do

      pvf = vv(na-1)
      pvd = vp(na-1)
      v = vh

      return
      end
      subroutine pbwa ( a, x, w1f, w1d, w2f, w2d )

c*********************************************************************72
c
cc PBWA computes parabolic cylinder functions W(a,x) and derivatives.
c
c  Licensing:
c
c    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
c    they give permission to incorporate this routine into a user program 
c    provided that the copyright is acknowledged.
c
c  Modified:
c
c    29 July 2012
c
c  Author:
c
c    Shanjie Zhang, Jianming Jin
c
c  Reference:
c
c    Shanjie Zhang, Jianming Jin,
c    Computation of Special Functions,
c    Wiley, 1996,
c    ISBN: 0-471-11963-6,
c    LC: QA351.C45.
c
c  Parameters:
c
c    Input, double precision A, the parameter.
c
c    Input, double precision X, the argument.
c
c    Output, double precision W1F, W1D, W2F, W2D, the values of
c    W(a,x), W'(a,x), W(a,-x), W'(a,-x).
c
      implicit none

      double precision a
      double precision d(100)
      double precision d1
      double precision d2
      double precision dl
      double precision eps
      double precision f1
      double precision f2
      double precision g1
      double precision g2
      double precision h(100)
      double precision h0
      double precision h1
      double precision hl
      integer k
      integer l1
      integer l2
      integer m
      double precision p0
      double precision r
      double precision r1
      double precision ugi
      double precision ugr
      double precision vgi
      double precision vgr
      double precision w1d
      double precision w1f
      double precision w2d
      double precision w2f
      double precision x
      double precision x1
      double precision x2
      double precision y1
      double precision y1d
      double precision y1f
      double precision y2d
      double precision y2f

      eps = 1.0D-15
      p0 = 0.59460355750136D+00

      if ( a .eq. 0.0D+00 ) then
        g1 = 3.625609908222D+00
        g2 = 1.225416702465D+00
      else
        x1 = 0.25D+00
        y1 = 0.5D+00 * a
        call cgama ( x1, y1, 1, ugr, ugi )
        g1 = sqrt ( ugr * ugr + ugi * ugi )
        x2 = 0.75D+00
        call cgama ( x2, y1, 1, vgr, vgi )
        g2 = sqrt ( vgr * vgr + vgi * vgi )
      end if

      f1 = sqrt ( g1 / g2 )
      f2 = sqrt ( 2.0D+00 * g2 / g1 )
      h0 = 1.0D+00
      h1 = a
      h(1) = a
      do l1 = 4, 200, 2
        m = l1 / 2
        hl = a * h1 - 0.25D+00 * ( l1 - 2.0D+00 ) 
     &    * ( l1 - 3.0D+00 ) * h0
        h(m) = hl
        h0 = h1
        h1 = hl
      end do
      y1f = 1.0D+00
      r = 1.0D+00
      do k = 1, 100
        r = 0.5D+00 * r * x * x / ( k * ( 2.0D+00 * k - 1.0D+00 ) )
        r1 = h(k) * r
        y1f = y1f + r1
        if ( abs ( r1 / y1f ) .le. eps .and. 30 .lt. k ) then
          go to 10
        end if
      end do

10    continue

      y1d = a
      r = 1.0D+00
      do k = 1, 100
        r = 0.5D+00 * r * x * x / ( k * ( 2.0D+00 * k + 1.0D+00 ) )
        r1 = h(k+1) * r
        y1d = y1d + r1
        if ( abs ( r1 / y1d ) .le. eps. and. 30 .lt. k ) then
          go to 20
        end if
      end do

20    continue

      y1d = x * y1d
      d1 = 1.0D+00
      d2 = a
      d(1) = 1.0D+00
      d(2) = a
      do l2 = 5, 160, 2
        m = ( l2 + 1 ) / 2
        dl = a * d2 - 0.25D+00 * ( l2 - 2.0D+00 ) 
     &    * ( l2 - 3.0D+00 ) * d1
        d(m) = dl
        d1 = d2
        d2 = dl
      end do

      y2f = 1.0D+00
      r = 1.0D+00
      do k = 1, 100
        r = 0.5D+00 * r * x * x / ( k * ( 2.0D+00 * k + 1.0D+00 ) )
        r1 = d(k+1) * r
        y2f = y2f + r1
        if ( abs ( r1 / y2f ) .le. eps .and. 30 .lt. k ) then
          go to 30
        end if
      end do

30    continue

      y2f = x * y2f
      y2d = 1.0D+00
      r = 1.0D+00
      do k = 1, 100
        r = 0.5D+00 * r * x * x / ( k * ( 2.0D+00 * k - 1.0D+00 ) )
        r1 = d(k+1) * r
        y2d = y2d + r1
        if ( abs ( r1 / y2d ) .le. eps .and. 30 .lt. k ) then
          go to 40
        end if
      end do

40    continue

      w1f = p0 * ( f1 * y1f - f2 * y2f )
      w2f = p0 * ( f1 * y1f + f2 * y2f )
      w1d = p0 * ( f1 * y1d - f2 * y2d )
      w2d = p0 * ( f1 * y1d + f2 * y2d )

      return
      end
      subroutine psi ( x, ps )

c*********************************************************************72
c
cc PSI computes the Psi function.
c
c  Licensing:
c
c    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
c    they give permission to incorporate this routine into a user program 
c    provided that the copyright is acknowledged.
c
c  Modified:
c
c    22 July 2012
c
c  Author:
c
c    Shanjie Zhang, Jianming Jin
c
c  Reference:
c
c    Shanjie Zhang, Jianming Jin,
c    Computation of Special Functions,
c    Wiley, 1996,
c    ISBN: 0-471-11963-6,
c    LC: QA351.C45.
c
c  Parameters:
c
c    Input, double precision X, the argument.
c
c    Output, double precision PS, the value of the function.
c
      implicit none

      double precision a1
      double precision a2
      double precision a3
      double precision a4
      double precision a5
      double precision a6
      double precision a7
      double precision a8
      double precision el
      integer k
      integer n
      double precision pi
      double precision ps
      double precision s
      double precision x
      double precision x2
      double precision xa

      xa = abs ( x )
      pi = 3.141592653589793D+00
      el = 0.5772156649015329D+00
      s = 0.0D+00

      if ( x .eq. int ( x ) .and. x .le. 0.0D+00 ) then

        ps = 1.0D+300
        return

      else if ( xa .eq. int ( xa ) ) then

        n = int ( xa )
        do k = 1, n - 1
          s = s + 1.0D+00 / k
        end do
        ps = - el + s

      else if ( xa + 0.5D+00 .eq. int ( xa + 0.5D+00 ) ) then

        n = xa - 0.5D+00
        do k = 1, n
          s = s + 1.0D+00 / ( 2.0D+00 * k - 1.0D+00 )
        end do
        ps = - el + 2.0D+00 * s - 1.386294361119891D+00

      else

        if ( xa .lt. 10.0D+00 ) then
          n = 10 - int ( xa )
          do k = 0, n - 1 
            s = s + 1.0D+00 / ( xa + k )
          end do
          xa = xa + n
        end if

        x2 = 1.0D+00 / ( xa * xa )
        a1 = -0.8333333333333D-01
        a2 = 0.83333333333333333D-02
        a3 = -0.39682539682539683D-02
        a4 = 0.41666666666666667D-02
        a5 = -0.75757575757575758D-02
        a6 = 0.21092796092796093D-01
        a7 = -0.83333333333333333D-01
        a8 = 0.4432598039215686D+00
        ps = log ( xa ) - 0.5D+00 / xa + x2 * (((((((
     &      a8   * x2
     &    + a7 ) * x2
     &    + a6 ) * x2
     &    + a5 ) * x2
     &    + a4 ) * x2
     &    + a3 ) * x2
     &    + a2 ) * x2
     &    + a1 )
        ps = ps-s

      end if

      if ( x .lt. 0.0D+00 ) then
        ps = ps - pi * cos ( pi * x ) / sin ( pi * x ) - 1.0D+00 / x
      end if

      return
      end
      subroutine qstar ( m, n, c, ck, ck1, qs, qt )

c*********************************************************************72
c
cc QSTAR computes Q*mn(-ic) for oblate radial functions with a small argument.
c
c  Licensing:
c
c    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
c    they give permission to incorporate this routine into a user program 
c    provided that the copyright is acknowledged.
c
c  Modified:
c
c    18 July 2012
c
c  Author:
c
c    Shanjie Zhang, Jianming Jin
c
c  Reference:
c
c    Shanjie Zhang, Jianming Jin,
c    Computation of Special Functions,
c    Wiley, 1996,
c    ISBN: 0-471-11963-6,
c    LC: QA351.C45.
c
c  Parameters:
c
c    Input, integer M, the mode parameter;  M = 0, 1, 2, ...
c
c    Input, integer N, mode parameter, N = M, M + 1, M + 2, ...
c
c    Input, double precision C, spheroidal parameter.
c
c    Input, double precision CK(*), ?
c
c    Input, double precision CK1, ?
c
c    Output, double precision QS, ?
c
c    Output, double precision QT, ?
c
      implicit none

      double precision ap(200)
      double precision c
      double precision ck(200)
      double precision ck1
      integer i
      integer ip
      integer k
      integer l
      integer m
      integer n
      double precision qs
      double precision qs0
      double precision qt
      double precision r
      double precision s
      double precision sk

      if ( n - m .eq. 2 * int ( ( n - m ) / 2 ) ) then
        ip = 0
      else
        ip = 1
      end if

      r = 1.0D+00 / ck(1) ** 2
      ap(1) = r
      do i = 1, m
        s = 0.0D+00
        do l = 1, i
          sk = 0.0D+00
          do k = 0, l
            sk = sk + ck(k+1) * ck(l-k+1)
          end do
          s = s + sk * ap(i-l+1)
        end do
        ap(i+1) = -r * s
      end do 

      qs0 = ap(m+1)     
      do l = 1, m
        r = 1.0D+00
        do k = 1, l
          r = r * ( 2.0D+00 * k + ip ) 
     &      * ( 2.0D+00 * k - 1.0D+00 + ip ) / ( 2.0D+00 * k ) ** 2
        end do
        qs0 = qs0 + ap(m-l+1) * r
      end do

      qs = ( -1.0D+00 ) ** ip * ck1 * ( ck1 * qs0 ) / c
      qt = - 2.0D+00 / ck1 * qs

      return
      end
      subroutine rctj ( n, x, nm, rj, dj )

c*********************************************************************72
c
cc RCTJ computes Riccati-Bessel function of the first kind, and derivatives.
c
c  Licensing:
c
c    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
c    they give permission to incorporate this routine into a user program 
c    provided that the copyright is acknowledged.
c
c  Modified:
c
c    18 July 2012
c
c  Author:
c
c    Shanjie Zhang, Jianming Jin
c
c  Reference:
c
c    Shanjie Zhang, Jianming Jin,
c    Computation of Special Functions,
c    Wiley, 1996,
c    ISBN: 0-471-11963-6,
c    LC: QA351.C45.
c
c  Parameters:
c
c    Input, integer N, the order of jn(x).
c
c    Input, double precision X, the argument.
c
c    Output, integer NM, the highest order computed.
c
c    Output, double precision RJ(0:N), the values of x jn(x).
c
c    Output, double precision DJ(0:N), the values of [x jn(x)]'.
c
      implicit none

      integer n

      double precision cs
      double precision dj(0:n)
      double precision f
      double precision f0
      double precision f1
      integer k
      integer m
      integer msta1
      integer msta2
      integer nm
      double precision rj(0:n)
      double precision rj0
      double precision rj1
      double precision x

      nm = n

      if ( abs ( x ) .lt. 1.0D-100 ) then
        do k = 0, n
          rj(k) = 0.0D+00
          dj(k) = 0.0D+00
        end do
        dj(0) = 1.0D+00
        return
      end if

      rj(0) = sin ( x )
      rj(1) = rj(0) / x - cos ( x )
      rj0 = rj(0)
      rj1 = rj(1)

      if ( 2 .le. n ) then

        m = msta1 ( x, 200 )

        if ( m .lt. n ) then
          nm = m
        else
          m = msta2 ( x, n, 15 )
        end if

        f0 = 0.0D+00
        f1 = 1.0D-100
        do k = m, 0, -1
          f = ( 2.0D+00 * k + 3.0D+00 ) * f1 / x - f0
          if ( k .le. nm ) then
            rj(k) = f
          end if
          f0 = f1
          f1 = f
        end do

        if ( abs ( rj1 ) .lt. abs ( rj0 ) ) then
          cs = rj0 / f
        else
          cs = rj1 / f0
        end if

        do k = 0, nm
          rj(k) = cs * rj(k)
        end do

      end if

      dj(0) = cos ( x )
      do k = 1, nm
        dj(k) = - k * rj(k) / x + rj(k-1)
      end do

      return
      end
      subroutine rcty ( n, x, nm, ry, dy )

c*********************************************************************72
c
cc RCTY computes Riccati-Bessel function of the second kind, and derivatives.
c
c  Licensing:
c
c    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
c    they give permission to incorporate this routine into a user program 
c    provided that the copyright is acknowledged.
c
c  Modified:
c
c    18 July 2012
c
c  Author:
c
c    Shanjie Zhang, Jianming Jin
c
c  Reference:
c
c    Shanjie Zhang, Jianming Jin,
c    Computation of Special Functions,
c    Wiley, 1996,
c    ISBN: 0-471-11963-6,
c    LC: QA351.C45.
c
c  Parameters:
c
c    Input, integer N, the order of yn(x).
c
c    Input, double precision X, the argument.
c
c    Output, integer NM, the highest order computed.
c
c    Output, double precision RY(0:N), the values of x yn(x).
c
c    Output, double precision DY(0:N), the values of [x yn(x)]'.
c
      implicit none

      integer n

      double precision dy(0:n)
      integer k
      integer nm
      double precision rf0
      double precision rf1
      double precision rf2
      double precision ry(0:n)
      double precision x

      nm = n

      if ( x .lt. 1.0D-60 ) then
        do k = 0, n
          ry(k) = -1.0D+300
          dy(k) = 1.0D+300
        end do
        ry(0) = -1.0D+00
        dy(0) = 0.0D+00
        return
      end if

      ry(0) = - cos ( x )
      ry(1) = ry(0) / x - sin ( x )
      rf0 = ry(0)
      rf1 = ry(1)
      do k = 2, n
        rf2 = ( 2.0D+00 * k - 1.0D+00 ) * rf1 / x - rf0
        if ( 1.0D+300 .lt. abs ( rf2 ) ) then
          go to 10
        end if
        ry(k) = rf2
        rf0 = rf1
        rf1 = rf2
      end do

10    continue

      nm = k - 1
      dy(0) = sin ( x )
      do k = 1, nm
        dy(k) = - k * ry(k) / x + ry(k-1)
      end do

      return
      end
      subroutine refine ( kd, m, q, a, iflag )

c*********************************************************************72
c
cc REFINE refines an estimate of the characteristic value of Mathieu functions.
c
c  Discussion:
c
c    This procedure calculates the accurate characteristic value
c    by the secant method.
c
c  Licensing:
c
c    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
c    they give permission to incorporate this routine into a user program 
c    provided that the copyright is acknowledged.
c
c  Modified:
c
c    20 July 2012
c
c  Author:
c
c    Shanjie Zhang, Jianming Jin
c
c  Reference:
c
c    Shanjie Zhang, Jianming Jin,
c    Computation of Special Functions,
c    Wiley, 1996,
c    ISBN: 0-471-11963-6,
c    LC: QA351.C45.
c
c  Parameters:
c
c    Input, integer KD, the case code:
c    1, for cem(x,q)  ( m = 0,2,4,...)
c    2, for cem(x,q)  ( m = 1,3,5,...)
c    3, for sem(x,q)  ( m = 1,3,5,...)
c    4, for sem(x,q)  ( m = 2,4,6,...)
c
c    Input, integer M, the order of the Mathieu functions.
c
c    Input, double precision Q, the parameter of the Mathieu functions.
c
c    Input/output, double precision A, the characteristic value, which
c    should have been refined on output.
c
      implicit none

      double precision a
      double precision ca
      double precision delta
      double precision eps
      double precision f
      double precision f0
      double precision f1
      integer it
      integer iflag
      integer kd
      integer m
      integer mj
      double precision q
      double precision x
      double precision x0
      double precision x1

      eps = 1.0D-14
      mj = 10 + m
      ca = a
      delta = 0.0D+00
      x0 = a
      call cvf ( kd, m, q, x0, mj, f0 )
      x1 = 1.002D+00 * a
      call cvf ( kd, m, q, x1, mj, f1 )

10    continue
  
      do it = 1, 100
        mj = mj + 1
        x = x1 - ( x1 - x0 ) / ( 1.0D+00 - f0 / f1 )
        call cvf ( kd, m, q, x, mj, f )
        if ( abs ( 1.0D+00 - x1 / x ) .lt. eps .or. 
     &    f .eq. 0.0D+00 ) then
          go to 20
        end if
        x0 = x1
        f0 = f1
        x1 = x
        f1 = f
      end do

20    continue

      a = x

      if ( 0.05D+00 .lt. delta ) then
        a = ca
        if ( iflag .lt. 0 ) then
          iflag = -10
        end if
        return
      end if

      if ( 0.05D+00 .lt. abs ( ( a - ca ) / ca) ) then
        x0 = ca
        delta = delta + 0.005D+00
        call cvf ( kd, m, q, x0, mj, f0 )
        x1 = ( 1.0D+00 + delta ) * ca
        call cvf ( kd, m, q, x1, mj, f1 )
        go to 10
      end if

      return
      end
      subroutine rmn1 ( m, n, c, x, df, kd, r1f, r1d )

c*********************************************************************72
c
cc RMN1 computes prolate and oblate spheroidal functions of the first kind.
c
c  Discussion:
c
c    This procedure computes prolate and oblate spheroidal radial
c    functions of the first kind for given m, n, c and x.
c
c  Licensing:
c
c    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
c    they give permission to incorporate this routine into a user program 
c    provided that the copyright is acknowledged.
c
c  Modified:
c
c    29 July 2012
c
c  Author:
c
c    Shanjie Zhang, Jianming Jin
c
c  Reference:
c
c    Shanjie Zhang, Jianming Jin,
c    Computation of Special Functions,
c    Wiley, 1996,
c    ISBN: 0-471-11963-6,
c    LC: QA351.C45.
c
c  Parameters:
c
c    Input, integer M, the mode parameter;  M = 0, 1, 2, ...
c
c    Input, integer N, mode parameter, N = M, M + 1, M + 2, ...
c
c    Input, double precision C, spheroidal parameter.
c
c    Input, double precision X, the argument.
c
c    Input, double precision DF(*), the expansion coefficients.
c
c    Input, integer KD, the function code.
c    1, the prolate function.
c    -1, the oblate function.
c
c    Output, double precision R1F, R1D, the function and derivative.
c
      implicit none

      double precision a0
      double precision b0
      double precision c
      double precision ck(200)
      double precision cx
      double precision df(200)
      double precision dj(0:251)
      double precision eps
      integer ip
      integer j
      integer k
      integer kd
      integer l
      integer lg
      integer m
      integer n
      integer nm
      integer nm1
      integer nm2
      integer np
      double precision r
      double precision r0
      double precision r1
      double precision r1d
      double precision r1f
      double precision r2
      double precision r3
      double precision reg
      double precision sa0
      double precision sj(0:251)
      double precision suc
      double precision sud
      double precision sum
      double precision sw
      double precision sw1
      double precision x

      eps = 1.0D-14
      nm1 = int ( ( n - m ) / 2 )
      if ( n - m .eq. 2 * nm1 ) then
        ip = 0
      else
        ip = 1
      end if
      nm = 25 + nm1 + int ( c )
      reg = 1.0D+00
      if ( 80 .lt. m + nm ) then
        reg = 1.0D-200
      end if
      r0 = reg
      do j = 1, 2 * m + ip
        r0 = r0 * j
      end do
      r = r0    
      suc = r * df(1)
      do k = 2, nm
        r = r * ( m + k - 1.0D+00 ) * ( m + k + ip - 1.5D+00 )
     &    / ( k - 1.0D+00 ) / ( k + ip - 1.5D+00 )
        suc = suc + r * df(k)

        if ( nm1 .lt. k .and. 
     &    abs ( suc - sw ) .lt. abs ( suc ) * eps ) then
          go to 10
        end if

        sw = suc

      end do

10    continue

      if ( x .eq. 0.0D+00 ) then

        call sckb ( m, n, c, df, ck )
        sum = 0.0D+00
        do j = 1, nm
          sum = sum + ck(j)
          if ( abs ( sum - sw1 ) .lt. abs ( sum ) * eps ) then
            go to 20
          end if
          sw1 = sum
        end do

20      continue

        r1 = 1.0D+00
        do j = 1, ( n + m + ip ) / 2
          r1 = r1 * ( j + 0.5D+00 * ( n + m + ip ) )
        end do

        r2 = 1.0D+00
        do j = 1, m
          r2 = 2.0D+00 * c * r2 * j
        end do

        r3 = 1.0D+00
        do j = 1, ( n - m - ip ) / 2
          r3 = r3 * j
        end do

        sa0 = ( 2.0D+00 * ( m + ip ) + 1.0D+00 ) * r1
     &    / ( 2.0D+00 ** n * c ** ip * r2 * r3 )

        if ( ip .eq. 0 ) then
          r1f = sum / ( sa0 * suc ) * df(1) * reg
          r1d = 0.0D+00
        else if ( ip .eq. 1 ) then
          r1f = 0.0D+00
          r1d = sum / ( sa0 * suc ) * df(1) * reg
        end if

        return

      end if

      cx = c * x
      nm2 = 2 * nm + m
      call sphj ( nm2, cx, nm2, sj, dj )
      a0 = ( 1.0D+00 - kd / ( x * x ) ) ** ( 0.5D+00 * m ) / suc  
      r1f = 0.0D+00
      do k = 1, nm
        l = 2 * k + m - n - 2 + ip
        if ( l .eq. 4 * int ( l / 4 ) ) then
          lg = 1
        else
          lg = -1
        end if
        if ( k .eq. 1 ) then
          r = r0
        else
          r = r * ( m + k - 1.0D+00 ) * ( m + k + ip - 1.5D+00 )
     &      / ( k - 1.0D+00 ) / ( k + ip - 1.5D+00 )
        end if
        np = m + 2 * k - 2 + ip
        r1f = r1f + lg * r * df(k) * sj(np)
        if ( nm1 .lt. k .and. 
     &    abs ( r1f - sw ) .lt. abs ( r1f ) * eps ) then
          go to 30
        end if
        sw = r1f
      end do

30    continue

      r1f = r1f * a0
      b0 = kd * m / x ** 3.0D+00 / ( 1.0D+00 - kd / ( x * x ) ) * r1f    
      sud = 0.0D+00

      do k = 1, nm

        l = 2 * k + m - n - 2 + ip

        if ( l .eq. 4 * int ( l / 4 ) ) then
          lg = 1
        else
          lg = -1
        end if

        if ( k .eq. 1 ) then
          r = r0
        else
          r = r * ( m + k - 1.0D+00 ) * ( m + k + ip - 1.5D+00 )
     &      / ( k - 1.0D+00 ) / ( k + ip - 1.5D+00 )
        end if

        np = m + 2 * k - 2 + ip
        sud = sud + lg * r * df(k) * dj(np)
        if ( nm1 .lt. k .and. 
     &    abs ( sud - sw ) .lt. abs ( sud ) * eps ) then
          go to 40
        end if
        sw = sud
      end do

40    continue

      r1d = b0 + a0 * c * sud

      return
      end
      subroutine rmn2l ( m, n, c, x, df, kd, r2f, r2d, id )

c*********************************************************************72
c
cc RMN2L: prolate and oblate spheroidal functions, second kind, large CX.
c
c  Discussion:
c
c    This procedure computes prolate and oblate spheroidal radial functions 
c    of the second kind for given m, n, c and a large cx.
c
c  Licensing:
c
c    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
c    they give permission to incorporate this routine into a user program 
c    provided that the copyright is acknowledged.
c
c  Modified:
c
c    30 July 2012
c
c  Author:
c
c    Shanjie Zhang, Jianming Jin
c
c  Reference:
c
c    Shanjie Zhang, Jianming Jin,
c    Computation of Special Functions,
c    Wiley, 1996,
c    ISBN: 0-471-11963-6,
c    LC: QA351.C45.
c
c  Parameters:
c
c    Input, integer M, the mode parameter;  M = 0, 1, 2, ...
c
c    Input, integer N, mode parameter, N = M, M + 1, M + 2, ...
c
c    Input, double precision C, spheroidal parameter.
c
c    Input, double precision X, the argument.
c
c    Input, double precision DF(*), the expansion coefficients.
c
c    Input, integer KD, the function code.
c    1, the prolate function.
c    -1, the oblate function.
c
c    Output, double precision R2F, R2D, the function and derivative values.
c
      implicit none

      double precision a0
      double precision b0
      double precision c
      double precision cx
      double precision df(200)
      double precision dy(0:251)
      double precision eps
      double precision eps1
      double precision eps2
      integer id
      integer id1
      integer id2
      integer ip
      integer j
      integer k
      integer kd
      integer l
      integer lg
      integer m
      integer n
      integer nm
      integer nm1
      integer nm2
      integer np
      double precision r
      double precision r0
      double precision r2d
      double precision r2f
      double precision reg
      double precision sw
      double precision suc
      double precision sud
      double precision sy(0:251)
      double precision x

      eps = 1.0D-14

      nm1 = int ( ( n - m ) / 2 )

      if ( n - m .eq. 2 * nm1 ) then
        ip = 0
      else
        ip = 1
      end if
      nm = 25 + nm1 + int ( c )

      if ( 80 .lt. m + nm ) then
        reg = 1.0D-200
      else
        reg = 1.0D+00
      end if
      nm2 = 2 * nm + m
      cx = c * x
      call sphy ( nm2, cx, nm2, sy, dy )
      r0 = reg
      do j = 1, 2 * m + ip
        r0 = r0 * j
      end do
      r = r0    
      suc = r * df(1)
      do k = 2, nm
        r = r * ( m + k - 1.0D+00 ) * ( m + k + ip - 1.5D+00 )
     &    / ( k - 1.0D+00 ) / ( k + ip - 1.5D+00 )
        suc = suc + r * df(k)
        if ( nm1 .lt. k .and. 
     &    abs ( suc - sw ) .lt. abs ( suc ) * eps ) then
          go to 10
        end if
        sw = suc
      end do

10    continue

      a0 = ( 1.0D+00 - kd / ( x * x ) ) ** ( 0.5D+00 * m ) / suc
      r2f = 0.0D+00
      do k = 1, nm
        l = 2 * k + m - n - 2 + ip
        if ( l .eq. 4 * int ( l / 4 ) ) then
          lg = 1
        else
          lg = -1
        end if

        if ( k .eq. 1 ) then
          r = r0
        else
          r = r * ( m + k - 1.0D+00 ) * ( m + k + ip - 1.5D+00 )
     &      / ( k - 1.0D+00 ) / ( k + ip - 1.5D+00 )
        end if

        np = m + 2 * k - 2 + ip
        r2f = r2f + lg * r * ( df(k) * sy(np) )
        eps1 = abs ( r2f - sw )
        if ( nm1 .lt. k .and. eps1 .lt. abs ( r2f ) * eps ) then
          go to 20
        end if
        sw = r2f
      end do

20    continue

      id1 = int ( log10 ( eps1 / abs ( r2f ) + eps ) )
      r2f = r2f * a0

      if ( nm2 .le. np ) then
        id = 10
        return
      end if

      b0 = kd * m / x ** 3.0D+00 / ( 1.0D+00 - kd / ( x * x ) ) * r2f                
      sud = 0.0D+00
      do k = 1, nm
        l = 2 * k + m - n - 2 + ip
        if ( l .eq. 4 * int ( l / 4 ) ) then
          lg = 1
        else
          lg = -1
        end if
        if (k .eq. 1) then
          r = r0
        else
          r = r * ( m + k - 1.0D+00 ) * ( m + k + ip - 1.5D+00 )
     &      / ( k - 1.0D+00 ) / ( k + ip - 1.5D+00 )
        end if
        np = m + 2 * k - 2 + ip
        sud = sud + lg * r * ( df(k) * dy(np) )
        eps2 = abs ( sud - sw )
        if ( nm1 .lt. k .and. eps2 .lt. abs ( sud ) * eps ) then
          go to 30
        end if
        sw = sud
      end do

30    continue

      r2d = b0 + a0 * c * sud
      id2 = int ( log10 ( eps2 / abs ( sud ) + eps ) )
      id = max ( id1, id2 )

      return
      end
      subroutine rmn2so ( m, n, c, x, cv, df, kd, r2f, r2d )

c*********************************************************************72
c
cc RMN2SO computes oblate radial functions of the second kind with small argument.
c
c  Discussion:
c
c    This procedure computes oblate radial functions of the second kind
c    with a small argument, Rmn(-ic,ix) and Rmn'(-ic,ix).
c
c  Licensing:
c
c    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
c    they give permission to incorporate this routine into a user program 
c    provided that the copyright is acknowledged.
c
c  Modified:
c
c    27 July 2012
c
c  Author:
c
c    Shanjie Zhang, Jianming Jin
c
c  Reference:
c
c    Shanjie Zhang, Jianming Jin,
c    Computation of Special Functions,
c    Wiley, 1996,
c    ISBN: 0-471-11963-6,
c    LC: QA351.C45.
c 
c  Parameters:
c
c    Input, integer M, the mode parameter;  M = 0, 1, 2, ...
c
c    Input, integer N, mode parameter, N = M, M + 1, M + 2, ...
c
c    Input, double precision C, spheroidal parameter.
c
c    Input, double precision X, the argument.
c
c    Input, double precision CV, the characteristic value.
c
c    Input, double precision DF(*), the expansion coefficients.
c
c    Input, integer KD, the function code.
c    1, the prolate function.
c    -1, the oblate function.
c
c    Output, double precision R2F, R2D, the values of Rmn(-ic,ix) and Rmn'(-ic,ix).
c
      implicit none

      double precision bk(200)
      double precision c
      double precision ck(200)
      double precision ck1
      double precision ck2
      double precision cv
      double precision df(200)
      double precision dn(200)
      double precision eps
      double precision gd
      double precision gf
      double precision h0
      integer ip
      integer j
      integer kd
      integer m
      integer n
      integer nm
      double precision pi
      double precision qs
      double precision qt
      double precision r1d
      double precision r1f
      double precision r2d
      double precision r2f
      double precision sum
      double precision sw
      double precision x

      if ( abs ( df(1) ) .le. 1.0D-280 ) then
        r2f = 1.0D+300
        r2d = 1.0D+300
        return
      end if

      eps = 1.0D-14
      pi = 3.141592653589793D+00
      nm = 25 + int ( ( n - m ) / 2 + c )
      if ( n - m .eq. 2 * int ( ( n - m ) / 2 ) ) then
        ip = 0
      else
        ip = 1
      end if

      call sckb ( m, n, c, df, ck )
      call kmn ( m, n, c, cv, kd, df, dn, ck1, ck2 )
      call qstar ( m, n, c, ck, ck1, qs, qt )
      call cbk ( m, n, c, cv, qt, ck, bk )

      if ( x .eq. 0.0D+00 ) then

        sum = 0.0D+00
        do j = 1, nm
          sum = sum + ck(j)
          if ( abs ( sum - sw ) .lt. abs ( sum ) * eps ) then
            go to 10
          end if
          sw = sum
        end do

10      continue

        if ( ip .eq. 0 ) then
          r1f = sum / ck1
          r2f = - 0.5D+00 * pi * qs * r1f
          r2d = qs * r1f + bk(1)
        else if ( ip .eq. 1 ) then
           r1d = sum / ck1
           r2f = bk(1)
           r2d = -0.5D+00 * pi * qs * r1d
        end if

        return

      else

        call gmn ( m, n, c, x, bk, gf, gd )
        call rmn1 ( m, n, c, x, df, kd, r1f, r1d )
        h0 = atan ( x ) - 0.5D+00 * pi
        r2f = qs * r1f * h0 + gf
        r2d = qs * ( r1d * h0 + r1f / ( 1.0D+00 + x * x ) ) + gd

      end if

      return
      end
      subroutine rmn2sp ( m, n, c, x, cv, df, kd, r2f, r2d )

c*********************************************************************72
c
cc RMN2SP: prolate, oblate spheroidal radial functions, kind 2, small argument.
c
c  Licensing:
c
c    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
c    they give permission to incorporate this routine into a user program 
c    provided that the copyright is acknowledged.
c
c  Modified:
c
c    28 July 2012
c
c  Author:
c
c    Shanjie Zhang, Jianming Jin
c
c  Reference:
c
c    Shanjie Zhang, Jianming Jin,
c    Computation of Special Functions,
c    Wiley, 1996,
c    ISBN: 0-471-11963-6,
c    LC: QA351.C45.
c
c  Parameters:
c
c    Input, integer M, the mode parameter;  M = 0, 1, 2, ...
c
c    Input, integer N, mode parameter, N = M, M + 1, M + 2, ...
c
c    Input, double precision C, spheroidal parameter.
c
c    Input, double precision X, the argument.
c
c    Input, double precision CV, the characteristic value.
c
c    Input, double precision DF(*), the expansion coefficients.
c
c    Input, integer KD, the function code.
c    1, the prolate function.
c    -1, the oblate function.
c
c    Output, double precision R2F, R2D, the values of the function and 
c    its derivative.
c
      implicit none

      double precision c
      double precision ck1
      double precision ck2
      double precision cv
      double precision df(200)
      double precision dn(200)
      double precision eps
      double precision ga
      double precision gb
      double precision gc
      integer ip
      integer j
      integer j1
      integer j2
      integer k
      integer kd
      integer ki
      integer l1
      integer m
      integer n
      integer nm
      integer nm1
      integer nm2
      integer nm3
      double precision pd(0:251)
      double precision pm(0:251)
      double precision qd(0:251)
      double precision qm(0:251)
      double precision r1
      double precision r2
      double precision r2d
      double precision r2f
      double precision r3
      double precision r4
      double precision sd
      double precision sd0
      double precision sd1
      double precision sd2
      double precision sdm
      double precision sf
      double precision spd1
      double precision spd2
      double precision spl
      double precision su0
      double precision su1
      double precision su2
      double precision sum
      double precision sw
      double precision x

      if ( abs ( df(1) ) .lt. 1.0D-280 ) then
        r2f = 1.0D+300
        r2d = 1.0D+300
        return
      end if

      eps = 1.0D-14

      nm1 = int ( ( n - m ) / 2 )

      if ( n - m .eq. 2 * nm1 ) then
        ip = 0
      else
        ip = 1
      end if

      nm = 25 + nm1 + int ( c )
      nm2 = 2 * nm + m
      call kmn ( m, n, c, cv, kd, df, dn, ck1, ck2 )
      call lpmns ( m, nm2, x, pm, pd )
      call lqmns ( m, nm2, x, qm, qd )

      su0 = 0.0D+00
      do k = 1, nm
        j = 2 * k - 2 + m + ip
        su0 = su0 + df(k) * qm(j)
        if ( nm1 .lt. k .and. 
     &    abs ( su0 - sw ) .lt. abs ( su0 ) * eps ) then
          go to 10                                                                  
        end if
        sw = su0
      end do

10    continue

      sd0 = 0.0D+00

      do k = 1, nm
        j = 2 * k - 2 + m + ip
        sd0 = sd0 + df(k) * qd(j)
        if ( nm1 .lt. k .and. 
     &    abs ( sd0 - sw ) .lt. abs ( sd0 ) * eps ) then
          go to 20
        end if
        sw = sd0
      end do

20    continue

      su1 = 0.0D+00
      sd1 = 0.0D+00
      do k = 1, m
        j = m - 2 * k + ip
        if ( j .lt. 0 ) then
          j = - j - 1
        end if
        su1 = su1 + dn(k) * qm(j)
        sd1 = sd1 + dn(k) * qd(j)
      end do

      ga = ( ( x - 1.0D+00 ) / ( x + 1.0D+00 ) ) ** ( 0.5D+00 * m )

      do k = 1, m

        j = m - 2 * k + ip

        if ( 0 .le. j ) then
          go to 30
        end if
        if ( j .lt. 0 ) then
          j = - j - 1
        end if 
        r1 = 1.0D+00
        do j1 = 1, j
          r1 = ( m + j1 ) * r1
        end do
        r2 = 1.0D+00
        do j2 = 1, m - j - 2
          r2 = j2 * r2
        end do
        r3 = 1.0D+00
        sf = 1.0D+00
        do l1 = 1, j
          r3 = 0.5D+00 * r3 * ( - j + l1 - 1.0D+00 ) * ( j + l1 ) 
     &      / ( ( m + l1 ) * l1 ) * ( 1.0D+00 - x )
          sf = sf + r3
        end do

        if ( m - j .le. 1 ) then
          gb = 1.0D+00
        else
          gb = ( m - j - 1.0D+00 ) * r2
        end if

        spl = r1 * ga * gb * sf
        su1 = su1 + ( -1 ) ** ( j + m ) * dn(k) * spl
        spd1 = m / ( x * x - 1.0D+00 ) * spl
        gc = 0.5D+00 * j * ( j + 1.0 ) / ( m + 1.0D+00 )
        sd = 1.0D+00
        r4 = 1.0D+00
        do l1 = 1, j - 1
          r4 = 0.5D+00 * r4 * ( - j + l1 ) * ( j + l1 + 1.0D+00 ) 
     &      / ( ( m + l1 + 1.0D+00 ) * l1 ) * ( 1.0D+00 - x )
          sd = sd + r4
        end do

        spd2 = r1 * ga * gb * gc * sd
        sd1 = sd1 + ( - 1 ) ** ( j + m ) * dn(k) * ( spd1 + spd2 )

30      continue

      end do

      su2 = 0.0D+00
      ki = ( 2 * m + 1 + ip ) / 2
      nm3 = nm + ki
      do k = ki, nm3
        j = 2 * k - 1 - m - ip
        su2 = su2 + dn(k) * pm(j)
        if ( m .lt. j .and. 
     &    abs ( su2 - sw ) .lt. abs ( su2 ) * eps ) then
          go to 40
        end if
        sw = su2
      end do

40    continue

      sd2 = 0.0D+00

      do k = ki, nm3
        j = 2 * k - 1 - m - ip
        sd2 = sd2 + dn(k) * pd(j)
        if ( m .lt. j .and. 
     &    abs ( sd2 - sw ) .lt. abs ( sd2 ) * eps ) then
          go to 50
        end if
        sw = sd2
      end do

50    continue

      sum = su0 + su1 + su2
      sdm = sd0 + sd1 + sd2
      r2f = sum / ck2
      r2d = sdm / ck2

      return
      end
      subroutine rswfo ( m, n, c, x, cv, kf, r1f, r1d, r2f, r2d )

c*********************************************************************72
c
cc RSWFO computes prolate spheroidal radial function of first and second kinds.
c
c  Discussion:
c
c    This procedure computes oblate radial functions of the first
c    and second kinds, and their derivatives.
c
c  Licensing:
c
c    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
c    they give permission to incorporate this routine into a user program 
c    provided that the copyright is acknowledged.
c
c  Modified:
c
c    07 July 2012
c
c  Author:
c
c    Shanjie Zhang, Jianming Jin
c
c  Reference:
c
c    Shanjie Zhang, Jianming Jin,
c    Computation of Special Functions,
c    Wiley, 1996,
c    ISBN: 0-471-11963-6,
c    LC: QA351.C45.
c 
c  Parameters:
c
c    Input, integer M, the mode parameter;  M = 0, 1, 2, ...
c
c    Input, integer N, mode parameter, N = M, M + 1, M + 2, ...
c
c    Input, double precision C, spheroidal parameter.
c
c    Input, double precision X, the argument.
c
c    Input, double precision CV, the characteristic value.
c
c    Input, integer KF, the function code.
c    1, for the first kind
c    2, for the second kind
c    3, for both the first and second kinds.
c
c    Output, double precision R1F, the radial function of the first kind;
c
c    Output, double precision R1D, the derivative of the radial function of
c    the first kind;
c
c    Output, double precision R2F, the radial function of the second kind;
c
c    Output, double precision R2D, the derivative of the radial function of
c    the second kind;
c
      implicit none

      double precision c
      double precision cv
      double precision df(200)
      integer id
      integer kd
      integer kf
      integer m
      integer n
      double precision r1d
      double precision r1f
      double precision r2d
      double precision r2f
      double precision x

      kd = -1
      call sdmn ( m, n, c, cv, kd, df )

      if ( kf .ne. 2 ) then
        call rmn1 ( m, n, c, x, df, kd, r1f, r1d )
      end if

      if ( 1 .lt. kf ) then
        id = 10
        if ( 1.0D-08 .lt. x ) then
          call rmn2l ( m, n, c, x, df, kd, r2f, r2d, id )
        end if
        if ( -1 .lt. id ) then
          call rmn2so ( m, n, c, x, cv, df, kd, r2f, r2d )
        end if
      end if

      return
      end
      subroutine rswfp ( m, n, c, x, cv, kf, r1f, r1d, r2f, r2d )

c*********************************************************************72
c
cc RSWFP computes prolate spheroidal radial function of first and second kinds.
c
c  Discussion:
c
c    This procedure computes prolate spheriodal radial functions of the
c    first and second kinds, and their derivatives.
c
c  Licensing:
c
c    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
c    they give permission to incorporate this routine into a user program 
c    provided that the copyright is acknowledged.
c
c  Modified:
c
c    07 July 2012
c
c  Author:
c
c    Shanjie Zhang, Jianming Jin
c
c  Reference:
c
c    Shanjie Zhang, Jianming Jin,
c    Computation of Special Functions,
c    Wiley, 1996,
c    ISBN: 0-471-11963-6,
c    LC: QA351.C45.
c
c  Parameters:
c
c    Input, integer M, the mode parameter;  M = 0, 1, 2, ...
c
c    Input, integer N, mode parameter, N = M, M + 1, M + 2, ...
c
c    Input, double precision C, spheroidal parameter.
c
c    Input, double precision X, the argument of the radial function, 1 < X.
c
c    Input, double precision CV, the characteristic value.
c
c    Input, integer KF, the function code.
c    1, for the first kind
c    2, for the second kind
c    3, for both the first and second kinds.
c
c    Output, double precision R1F, the radial function of the first kind;
c
c    Output, double precision R1D, the derivative of the radial function of
c    the first kind;
c
c    Output, double precision R2F, the radial function of the second kind;
c
c    Output, double precision R2D, the derivative of the radial function of
c    the second kind;
c
      implicit none

      double precision c
      double precision cv
      double precision df(200)
      integer id
      integer kd
      integer kf
      integer m
      integer n
      double precision r1d
      double precision r1f
      double precision r2d
      double precision r2f
      double precision x

      kd = 1
      call sdmn ( m, n, c, cv, kd, df )

      if ( kf .ne. 2 ) then
        call rmn1 ( m, n, c, x, df, kd, r1f, r1d )
      end if

      if ( 1 .lt. kf ) then
        call rmn2l ( m, n, c, x, df, kd, r2f, r2d, id )
        if ( -8 .lt. id ) then
          call rmn2sp ( m, n, c, x, cv, df, kd, r2f, r2d )
        end if
      end if

      return
      end
      subroutine scka ( m, n, c, cv, kd, ck )

c*********************************************************************72
c
cc SCKA: expansion coefficients for prolate and oblate spheroidal functions.
c
c  Licensing:
c
c    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
c    they give permission to incorporate this routine into a user program 
c    provided that the copyright is acknowledged.
c
c  Modified:
c
c    22 July 2012
c
c  Author:
c
c    Shanjie Zhang, Jianming Jin
c
c  Reference:
c
c    Shanjie Zhang, Jianming Jin,
c    Computation of Special Functions,
c    Wiley, 1996,
c    ISBN: 0-471-11963-6,
c    LC: QA351.C45.
c
c  Parameters:
c
c    Input, integer M, the mode parameter.
c
c    Input, integer N, the mode parameter.
c
c    Input, double precision C, the spheroidal parameter.
c
c    Input, double precision CV, the characteristic value.
c
c    Input, integer KD, the function code.
c    1, the prolate function.
c    -1, the oblate function.
c
c    Output, double precision CK(*), the expansion coefficients.
c    CK(1), CK(2),... correspond to c0, c2,..., and so on.
c       
      implicit none

      double precision c
      double precision ck(200)
      double precision cs
      double precision cv
      double precision f
      double precision f0
      double precision f1
      double precision f2
      double precision fl
      double precision fs
      integer ip
      integer j
      integer k
      integer k1
      integer kb
      integer kd
      integer m
      integer n
      integer nm
      double precision r1
      double precision r2
      double precision s0
      double precision su1
      double precision su2

      if ( c .le. 1.0D-10 ) then
        c = 1.0D-10
      end if

      nm = 25 + int ( ( n - m ) / 2 + c )
      cs = c * c * kd

      if ( n - m .eq. 2 * int ( ( n - m ) / 2 ) ) then
        ip = 0
      else
        ip = 1
      end if

      fs = 1.0D+00
      f1 = 0.0D+00
      f0 = 1.0D-100
      kb = 0
      ck(nm+1) = 0.0D+00

      do k = nm, 1, -1

        f = ((( 2.0D+00 * k + m + ip ) 
     &    * ( 2.0D+00 * k + m + 1.0D+00 + ip ) - cv + cs ) * f0
     &    - 4.0D+00 * ( k + 1.0D+00 ) * ( k + m + 1.0D+00 ) * f1 ) / cs

        if ( abs ( ck(k+1) ) .lt. abs ( f ) ) then

          ck(k) = f
          f1 = f0
          f0 = f

          if ( 1.0D+100 .lt. abs ( f ) ) then
            do k1 = nm, k, -1
              ck(k1) = ck(k1) * 1.0D-100
            end do
            f1 = f1 * 1.0D-100
            f0 = f0 * 1.0D-100
          end if

        else

          kb = k
          fl = ck(k+1)
          f1 = 1.0D+00
          f2 = 0.25D+00 * ( ( m + ip ) * ( m + ip + 1.0D+00 ) 
     &      - cv + cs ) / ( m + 1.0D+00 ) * f1
          ck(1) = f1

          if ( kb .eq. 1 ) then
            fs = f2
          else if (kb .eq. 2 ) then
            ck(2) = f2
            fs = 0.125D+00 * ( ( ( m + ip + 2.0D+00 ) 
     &        * ( m + ip + 3.0D+00 ) - cv + cs ) * f2
     &        - cs * f1 ) / ( m + 2.0D+00 )
          else
            ck(2) = f2
            do j = 3, kb + 1
              f = 0.25D+00 * ( ( ( 2.0D+00 * j + m + ip - 4.0D+00 )
     &          * ( 2.0D+00 * j + m + ip - 3.0D+00 ) - cv + cs ) * f2
     &          - cs * f1 ) / ( ( j - 1.0D+00 ) * ( j + m - 1.0D+00 ) )
              if ( j .le. kb ) then
                ck(j) = f
              end if
              f1 = f2
              f2 = f
            end do
            fs = f
          end if

          go to 10

        end if

      end do

10    continue

      su1 = 0.0D+00
      do k = 1, kb
        su1 = su1 + ck(k)
      end do

      su2 = 0.0D+00
      do k = kb + 1, nm
        su2 = su2 + ck(k)
      end do

      r1 = 1.0D+00
      do j = 1, ( n + m + ip ) / 2
        r1 = r1 * ( j + 0.5D+00 * ( n + m + ip ) )
      end do

      r2 = 1.0D+00
      do j = 1, ( n - m - ip ) / 2
        r2 = - r2 * j
      end do

      if ( kb .eq. 0 ) then
        s0 = r1 / ( 2.0D+00 ** n * r2 * su2 )
      else
        s0 = r1 / ( 2.0D+00 ** n * r2 * ( fl / fs * su1 + su2 ) )
      end if

      do k = 1, kb
        ck(k) = fl / fs * s0 * ck(k)
      end do

      do k = kb + 1, nm
        ck(k) = s0 * ck(k)
      end do

      return
      end
      subroutine sckb ( m, n, c, df, ck )

c*********************************************************************72
c
cc SCKB: expansion coefficients for prolate and oblate spheroidal functions.
c
c  Licensing:
c
c    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
c    they give permission to incorporate this routine into a user program 
c    provided that the copyright is acknowledged.
c
c  Modified:
c
c    15 July 2012
c
c  Author:
c
c    Shanjie Zhang, Jianming Jin
c
c  Reference:
c
c    Shanjie Zhang, Jianming Jin,
c    Computation of Special Functions,
c    Wiley, 1996,
c    ISBN: 0-471-11963-6,
c    LC: QA351.C45.
c
c  Parameters:
c
c    Input, integer M, the mode parameter.
c
c    Input, integer N, the mode parameter.
c
c    Input, double precision C, the spheroidal parameter.
c
c    Input, double precision DF(*), the expansion coefficients DF.
c
c    Output, double precision CK(*), the expansion coefficients CK.
c
      implicit none

      double precision c
      double precision ck(200)
      double precision d1
      double precision d2
      double precision d3
      double precision df(200)
      double precision fac
      integer i
      integer i1
      integer i2
      integer ip
      integer k
      integer m
      integer n
      integer nm
      double precision r
      double precision r1
      double precision reg
      double precision sum
      double precision sw

      c = max ( c, 1.0D-10 )

      nm = 25 + int ( 0.5D+00 * ( n - m ) + c )
 
      if ( n - m .eq. 2 * int ( ( n - m ) / 2 ) ) then
        ip = 0
      else
        ip = 1
      end if

      if ( 80 .lt. m + nm ) then
        reg = 1.0D-200
      else
        reg = 1.0D+00
      end if

      fac = - 0.5D+00 ** m

      do k = 0, nm - 1

        fac = - fac
        i1 = 2 * k + ip + 1
        r = reg
        do i = i1, i1 + 2 * m - 1
          r = r * i
        end do 

        i2 = k + m + ip
        do i = i2, i2 + k - 1
          r = r * ( i + 0.5D+00 )
        end do

        sum = r * df(k+1)
        do i = k + 1, nm
          d1 = 2.0D+00 * i + ip
          d2 = 2.0D+00 * m + d1
          d3 = i + m + ip - 0.5D+00
          r = r * d2 * ( d2 - 1.0D+00 ) * i * ( d3 + k ) 
     &      / ( d1 * ( d1 - 1.0D+00 ) * ( i - k ) * d3 )
          sum = sum + r * df(i+1)
          if ( abs ( sw - sum ) .lt. abs ( sum ) * 1.0D-14 ) then
            go to 10
          end if
          sw = sum
        end do

10      continue

        r1 = reg
        do i = 2, m + k
          r1 = r1 * i
        end do

        ck(k+1) = fac * sum / r1

      end do

      return
      end
      subroutine sdmn ( m, n, c, cv, kd, df )

c*********************************************************************72
c
cc SDMN: expansion coefficients for prolate and oblate spheroidal functions.
c
c  Licensing:
c
c    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
c    they give permission to incorporate this routine into a user program 
c    provided that the copyright is acknowledged.
c
c  Modified:
c
c    29 July 2012
c
c  Author:
c
c    Shanjie Zhang, Jianming Jin
c
c  Reference:
c
c    Shanjie Zhang, Jianming Jin,
c    Computation of Special Functions,
c    Wiley, 1996,
c    ISBN: 0-471-11963-6,
c    LC: QA351.C45.
c
c  Parameters:
c
c    Input, integer M, the mode parameter.
c
c    Input, integer N, the mode parameter.
c
c    Input, double precision C, the spheroidal parameter.
c
c    Input, double precision CV, the characteristic value.
c
c    Input, integer KD, the function code.
c    1, the prolate function.
c    -1, the oblate function.
c
c    Output, double precision DF(*), expansion coefficients;
c    DF(1), DF(2), ... correspond to d0, d2, ... for even n-m and d1,
c    d3, ... for odd n-m
c
      implicit none

      double precision a(200)
      double precision c
      double precision cs
      double precision cv
      double precision d(200)
      double precision d2k
      double precision df(200)
      double precision dk0
      double precision dk1
      double precision dk2
      double precision f
      double precision f0
      double precision f1
      double precision f2
      double precision fl
      double precision fs
      double precision g(200)
      integer i
      integer ip
      integer j
      integer k
      integer k1
      integer kb
      integer kd
      integer m
      integer n
      integer nm
      double precision r1
      double precision r3
      double precision r4
      double precision s0
      double precision su1
      double precision su2
      double precision sw

      nm = 25 + int ( 0.5D+00 * ( n - m ) + c )

      if ( c .lt. 1.0D-10 ) then
         do i = 1, nm
           df(i) = 0D+00
         end do
         df((n-m)/2+1) = 1.0D+00
         return
      end if   

      cs = c * c * kd

      if ( n - m .eq. 2 * int ( ( n - m ) / 2 ) ) then
        ip = 0
      else
        ip = 1
      end if

      do i = 1, nm + 2
        if ( ip .eq. 0 ) then
          k = 2 * ( i - 1 )
        else
          k = 2 * i - 1
        end if
        dk0 = m + k
        dk1 = m + k + 1
        dk2 = 2 * ( m + k )
        d2k = 2 * m + k
        a(i) = ( d2k + 2.0D+00 ) * ( d2k + 1.0D+00 ) 
     &    / ( ( dk2 + 3.0D+00 ) * ( dk2 + 5.0D+00 ) ) * cs
        d(i) = dk0 * dk1 
     &    + ( 2.0D+00 * dk0 * dk1 - 2.0D+00 * m * m - 1.0D+00 )
     &    / ( ( dk2 - 1.0D+00 ) * ( dk2 + 3.0D+00 ) ) * cs
        g(i) = k * ( k - 1.0D+00 ) / ( ( dk2 - 3.0D+00 )
     &    * ( dk2 - 1.0D+00 ) ) * cs
      end do

      fs = 1.0D+00
      f1 = 0.0D+00
      f0 = 1.0D-100
      kb = 0
      df(nm+1) = 0.0D+00

      do k = nm, 1, -1

        f = - ( ( d(k+1) - cv ) * f0 + a(k+1) * f1 ) / g(k+1)

        if ( abs ( df(k+1) ) .lt. abs ( f ) ) then

          df(k) = f
          f1 = f0
          f0 = f
          if ( 1.0D+100 .lt. abs ( f ) ) then
            do k1 = k, nm
              df(k1) = df(k1) * 1.0D-100
            end do
            f1 = f1 * 1.0D-100
            f0 = f0 * 1.0D-100
          end if  

        else

          kb = k
          fl = df(k+1)
          f1 = 1.0D-100
          f2 = - ( d(1) - cv ) / a(1) * f1
          df(1) = f1

          if ( kb .eq. 1 ) then

            fs = f2

          else if ( kb .eq. 2 ) then

            df(2) = f2
            fs = - ( ( d(2) - cv ) * f2 + g(2) * f1 ) / a(2)

          else 

            df(2) = f2
            do j = 3, kb + 1
              f = - ( ( d(j-1) - cv ) * f2 + g(j-1) * f1 ) / a(j-1)
              if ( j .le. kb ) then
                df(j) = f
              end if
              if ( 1.0D+100 .lt. abs ( f ) ) then
                do k1 = 1, j
                  df(k1) = df(k1) * 1.0D-100
                end do
                f = f * 1.0D-100
                f2 = f2 * 1.0D-100
              end if  
              f1 = f2
              f2 = f
            end do
            fs = f

          end if

          go to 10
        end if

      end do

10    continue

      su1 = 0.0D+00

      r1 = 1.0D+00
      do j = m + ip + 1, 2 * ( m + ip )
        r1 = r1 * j
      end do

      su1 = df(1) * r1
      do k = 2, kb
        r1 = - r1 * ( k + m + ip - 1.5D+00 ) / ( k - 1.0D+00 )
        su1 = su1 + r1 * df(k)
      end do

      su2 = 0.0D+00
      do k = kb + 1, nm
        if ( k .ne. 1 ) then
          r1 = - r1 * ( k + m + ip - 1.5D+00 ) / ( k - 1.0D+00 )
        end if
        su2 = su2 + r1 * df(k)
        if ( abs ( sw - su2 ) .lt. abs ( su2 ) * 1.0D-14 ) then
          go to 20
        end if
        sw = su2
      end do

20    continue

      r3 = 1.0D+00
      do j = 1, ( m + n + ip ) / 2
        r3 = r3 * ( j + 0.5D+00 * ( n + m + ip ) )
      end do

      r4 = 1.0D+00
      do j = 1, ( n - m - ip ) / 2
        r4 = -4.0D+00 * r4 * j
      end do

      s0 = r3 / ( fl * ( su1 / fs ) + su2 ) / r4
      do k = 1, kb
        df(k) = fl / fs * s0 * df(k)
      end do

      do k = kb + 1, nm
        df(k) = s0 * df(k)
      end do

      return
      end
      subroutine segv ( m, n, c, kd, cv, eg )

c*********************************************************************72
c
cc SEGV computes the characteristic values of spheroidal wave functions.
c
c  Licensing:
c
c    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
c    they give permission to incorporate this routine into a user program 
c    provided that the copyright is acknowledged.
c
c  Modified:
c
c    28 July 2012
c
c  Author:
c
c    Shanjie Zhang, Jianming Jin
c
c  Reference:
c
c    Shanjie Zhang, Jianming Jin,
c    Computation of Special Functions,
c    Wiley, 1996,
c    ISBN: 0-471-11963-6,
c    LC: QA351.C45.
c
c  Parameters:
c
c    Input, integer M, the mode parameter.
c
c    Input, integer N, the mode parameter.
c
c    Input, double precision C, the spheroidal parameter.
c
c    Input, integer KD, the function code.
c    1, the prolate function.
c    -1, the oblate function.
c
c    Output, double precision CV, the characteristic value.
c
c    Output, double precision EG(*), the characteristic value for 
c    mode parameters m and n.  ( L = n - m + 1 )
c
      implicit none

      double precision a(300)
      double precision b(100)
      double precision c
      double precision cs
      double precision cv
      double precision cv0(100)
      double precision d(300)
      double precision d2k
      double precision dk0
      double precision dk1
      double precision dk2
      double precision e(300)
      double precision eg(200)
      double precision f(300)
      double precision g(300)
      double precision h(100)
      integer i
      integer icm
      integer j
      integer k
      integer k1
      integer kd
      integer l
      integer m
      integer n
      integer nm
      integer nm1
      double precision s
      double precision t
      double precision t1
      double precision x1
      double precision xa
      double precision xb

      if ( c .lt. 1.0D-10 ) then
        do i = 1, n
          eg(i) = ( i + m ) * ( i + m - 1.0D+00 )
        end do
        cv = eg(n-m+1)
        return
      end if

      icm = ( n - m + 2 ) / 2
      nm = 10 + int ( 0.5D+00 * ( n - m ) + c )
      cs = c * c * kd

      do l = 0, 1

        do i = 1, nm
          if ( l .eq. 0 ) then
            k = 2 * ( i - 1 )
          else
            k = 2 * i - 1
          end if
          dk0 = m + k
          dk1 = m + k + 1
          dk2 = 2 * ( m + k )
          d2k = 2 * m + k
          a(i) = ( d2k + 2.0D+00 ) * ( d2k + 1.0D+00 ) 
     &      / ( ( dk2 + 3.0D+00 ) * ( dk2 + 5.0D+00 ) ) * cs
          d(i) = dk0 * dk1 + ( 2.0D+00 * dk0 * dk1 
     &      - 2.0 * m * m - 1.0D+00 )
     &      / ( ( dk2 - 1.0D+00 ) * ( dk2 + 3.0D+00 ) ) * cs
          g(i) = k * ( k - 1.0D+00 ) / ( ( dk2 - 3.0D+00 ) 
     &      * ( dk2 - 1.0D+00 ) ) * cs
        end do

        do k = 2, nm
          e(k) = sqrt ( a(k-1) * g(k) )
          f(k) = e(k) * e(k)
        end do

        f(1) = 0.0D+00
        e(1) = 0.0D+00
        xa = d(nm) + abs ( e(nm) )
        xb = d(nm) - abs ( e(nm) )
        nm1 = nm - 1
        do i = 1, nm1
          t = abs ( e(i) ) + abs ( e(i+1) )
          t1 = d(i) + t
          if ( xa .lt. t1 ) then
            xa = t1
          end if
          t1 = d(i) - t
          if ( t1 .lt. xb ) then
            xb = t1
          end if
        end do

        do i = 1, icm
          b(i) = xa
          h(i) = xb
        end do

        do k = 1, icm

          do k1 = k, icm
            if ( b(k1) .lt. b(k) ) then
              b(k) = b(k1)
              go to 10
            end if
          end do

10        continue

          if ( k .ne. 1 .and. h(k) .lt. h(k-1) ) then
            h(k) = h(k-1)
          end if

20        continue

          x1 = ( b(k) + h(k) ) /2.0D+00
          cv0(k) = x1

          if ( abs ( ( b(k) - h(k) ) / x1 ) .lt. 1.0D-14 ) then
            go to 30
          end if

          j = 0
          s = 1.0D+00
          do i = 1, nm
            if ( s .eq. 0.0D+00 ) then
              s = s + 1.0D-30
            end if
            t = f(i) / s
            s = d(i) - t - x1
            if ( s .lt. 0.0D+00 ) then
              j = j + 1
            end if
          end do

          if ( j .lt. k ) then

            h(k) = x1

          else

            b(k) = x1
            if ( icm .le. j ) then
              b(icm) = x1
            else
              if ( h(j+1) .lt. x1 ) then
                h(j+1) = x1
              end if
              if ( x1 .lt. b(j) ) then
                b(j) = x1
              end if
            end if

          end if

          go to 20

30        continue

          cv0(k) = x1

          if ( l .eq. 0 ) then
            eg(2*k-1) = cv0(k)
          else
            eg(2*k) = cv0(k)
          end if

        end do

      end do

      cv = eg(n-m+1)

      return
      end
      subroutine sphi ( n, x, nm, si, di )

c*********************************************************************72
c
cc SPHI computes spherical Bessel functions in(x) and their derivatives in'(x).
c
c  Licensing:
c
c    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
c    they give permission to incorporate this routine into a user program 
c    provided that the copyright is acknowledged.
c
c  Modified:
c
c    18 July 2012
c
c  Author:
c
c    Shanjie Zhang, Jianming Jin
c
c  Reference:
c
c    Shanjie Zhang, Jianming Jin,
c    Computation of Special Functions,
c    Wiley, 1996,
c    ISBN: 0-471-11963-6,
c    LC: QA351.C45.
c
c  Parameters:
c
c    Input, integer N, the order of In(X).
c
c    Input, double precision X, the argument.
c
c    Output, integer NM, the highest order computed.
c
c    Output, double precision SI(0:N), DI(0:N), the values and derivatives
c    of the function of orders 0 through N.
c
      implicit none

      integer n

      double precision cs
      double precision di(0:n)
      double precision f
      double precision f0
      double precision f1
      integer k
      integer m
      integer msta1
      integer msta2
      integer nm
      double precision si(0:n)
      double precision si0
      double precision x

      nm = n

      if ( abs ( x ) .lt. 1.0D-100 ) then
        do k = 0, n
          si(k) = 0.0D+00
          di(k) = 0.0D+00
        end do
        si(0) = 1.0D+00
        di(1) = 0.333333333333333D+00
        return
      end if

      si(0) = dsinh ( x ) / x
      si(1) = -( dsinh ( x ) / x - dcosh ( x ) ) / x
      si0 = si(0)

      if ( 2 .le. n ) then

        m = msta1 ( x, 200 )
        if ( m .lt. n ) then
          nm = m
        else
          m = msta2 ( x, n, 15 )
        end if
        f0 = 0.0D+00
        f1 = 1.0D+00-100
        do k = m, 0, -1
          f = ( 2.0D+00 * k + 3.0D+00 ) * f1 / x + f0
          if ( k .le. nm ) then
            si(k) = f
          end if
          f0 = f1
          f1 = f
        end do
        cs = si0 / f
        do k = 0, nm
          si(k) = cs * si(k)
        end do

      end if

      di(0) = si(1)
      do k = 1, nm
        di(k) = si(k-1) - ( k + 1.0D+00 ) / x * si(k)
      end do

      return
      end
      subroutine sphj ( n, x, nm, sj, dj )

c*********************************************************************72
c
cc SPHJ computes spherical Bessel functions jn(x) and their derivatives.
c
c  Licensing:
c
c    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
c    they give permission to incorporate this routine into a user program 
c    provided that the copyright is acknowledged.
c
c  Modified:
c
c    15 July 2012
c
c  Author:
c
c    Shanjie Zhang, Jianming Jin
c
c  Reference:
c
c    Shanjie Zhang, Jianming Jin,
c    Computation of Special Functions,
c    Wiley, 1996,
c    ISBN: 0-471-11963-6,
c    LC: QA351.C45.
c
c  Parameters:
c
c    Input, integer N, the order.
c
c    Input, double precision X, the argument.
c
c    Output, integer NM, the highest order computed.
c
c    Output, double precision SJ(0:N), the values of jn(x).
c
c    Output, double precision DJ(0:N), the values of jn'(x).
c
      implicit none

      integer n

      double precision cs
      double precision dj(0:n)
      double precision f
      double precision f0
      double precision f1
      integer k
      integer m
      integer msta1
      integer msta2
      integer nm
      double precision sa
      double precision sb
      double precision sj(0:n)
      double precision x

      nm = n

      if ( abs ( x ) .le. 1.0D-100 ) then
        do k = 0, n
          sj(k) = 0.0D+00
          dj(k) = 0.0D+00
        end do
        sj(0) = 1.0D+00
        dj(1) = 0.3333333333333333D+00
        return
      end if

      sj(0) = sin ( x ) / x
      sj(1) = ( sj(0) - cos ( x ) ) / x

      if ( 2 .le. n ) then

        sa = sj(0)
        sb = sj(1)
        m = msta1 ( x, 200 )
        if ( m .lt. n ) then
          nm = m
        else
          m = msta2 ( x, n, 15 )
        end if

        f0 = 0.0D+00
        f1 = 1.0D+00-100
        do k = m, 0, -1
          f = ( 2.0D+00 * k + 3.0D+00 ) * f1 / x - f0
          if ( k .le. nm ) then
            sj(k) = f
          end if
          f0 = f1
          f1 = f
        end do

        if ( abs ( sa ) .le. abs ( sb ) ) then
          cs = sb / f0
        else
          cs = sa / f
        end if

        do k = 0, nm
          sj(k) = cs * sj(k)
        end do

      end if      

      dj(0) = ( cos(x) - sin(x) / x ) / x
      do k = 1, nm
        dj(k) = sj(k-1) - ( k + 1.0D+00 ) * sj(k) / x
      end do

      return
      end
      subroutine sphk ( n, x, nm, sk, dk )

c*********************************************************************72
c
cc SPHK computes modified spherical Bessel functions kn(x) and their derivatives.
c
c  Discussion:
c
c    This procedure computes modified spherical Bessel functions
c    of the second kind, kn(x) and kn'(x).
c
c  Licensing:
c
c    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
c    they give permission to incorporate this routine into a user program 
c    provided that the copyright is acknowledged.
c
c  Modified:
c
c    15 July 2012
c
c  Author:
c
c    Shanjie Zhang, Jianming Jin
c
c  Reference:
c
c    Shanjie Zhang, Jianming Jin,
c    Computation of Special Functions,
c    Wiley, 1996,
c    ISBN: 0-471-11963-6,
c    LC: QA351.C45.
c
c  Parameters:
c
c    Input, integer N, the order.
c
c    Input, double precision X, the argument.
c
c    Output, integer NM, the highest order computed.
c
c    Output, double precision SK(0:N), DK(0:N), the values of kn(x) and kn'(x).
c
      implicit none

      integer n

      double precision dk(0:n)
      double precision f
      double precision f0
      double precision f1
      integer k
      integer nm
      double precision sk(0:n)
      double precision pi
      double precision x

      pi = 3.141592653589793D+00
      nm = n
      if ( x .lt. 1.0D-60 ) then
        do k = 0,n
          sk(k) = 1.0D+300
          dk(k) = -1.0D+300
        end do
        return
      end if

      sk(0) = 0.5D+00 * pi / x * exp ( - x )
      sk(1) = sk(0) * ( 1.0D+00 + 1.0D+00 / x )
      f0 = sk(0)
      f1 = sk(1)
      do k = 2, n
        f = ( 2.0D+00 * k - 1.0D+00 ) * f1 / x + f0
        sk(k) = f
        if ( 1.0D+300 .lt. abs ( f ) ) then
          go to 10
        end if
        f0 = f1
        f1 = f
      end do

10    continue

      nm = k - 1

      dk(0) = -sk(1)
      do k = 1, nm
        dk(k) = -sk(k-1) - ( k + 1.0D+00 ) / x * sk(k)
      end do

      return
      end
      subroutine sphy ( n, x, nm, sy, dy )

c*********************************************************************72
c
cc SPHY computes spherical Bessel functions yn(x) and their derivatives.
c
c  Licensing:
c
c    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
c    they give permission to incorporate this routine into a user program 
c    provided that the copyright is acknowledged.
c
c  Modified:
c
c    15 July 2012
c
c  Author:
c
c    Shanjie Zhang, Jianming Jin
c
c  Reference:
c
c    Shanjie Zhang, Jianming Jin,
c    Computation of Special Functions,
c    Wiley, 1996,
c    ISBN: 0-471-11963-6,
c    LC: QA351.C45.
c
c  Parameters:
c
c    Input, integer N, the order.
c
c    Input, double precision X, the argument.
c
c    Output, integer NM, the highest order computed.
c
c    Output, double precision SY(0:N), DY(0:N), the values of yn(x) and yn'(x).
c 
      implicit none

      integer n

      double precision dy(0:n)
      double precision f
      double precision f0
      double precision f1
      integer k
      integer nm
      double precision sy(0:n)
      double precision x

      nm = n

      if ( x .lt. 1.0D-60 ) then
        do k = 0, n
          sy(k) = -1.0D+300
          dy(k) = 1.0D+300
        end do
        return
      end if

      sy(0) = -cos ( x ) / x
      sy(1) = ( sy(0) - sin ( x ) ) / x
      f0 = sy(0)
      f1 = sy(1)
      do k = 2, n
        f = ( 2.0D+00 * k - 1.0D+00 ) * f1 / x - f0
        sy(k) = f
        if ( 1.0D+300 .le. abs ( f ) ) then
          go to 10
        end if              
        f0 = f1
        f1 = f
      end do

10    continue

      nm = k - 1
      dy(0) = ( sin ( x ) + cos ( x ) / x ) / x
      do k = 1, nm
        dy(k) = sy(k-1) - ( k + 1.0D+00 ) * sy(k) / x
      end do

      return
      end
      subroutine stvh0 ( x, sh0 )

c*********************************************************************72
c
cc STVH0 computes the Struve function H0(x).
c
c  Licensing:
c
c    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
c    they give permission to incorporate this routine into a user program 
c    provided that the copyright is acknowledged.
c
c  Modified:
c
c    22 July 2012
c
c  Author:
c
c    Shanjie Zhang, Jianming Jin
c
c  Reference:
c
c    Shanjie Zhang, Jianming Jin,
c    Computation of Special Functions,
c    Wiley, 1996,
c    ISBN: 0-471-11963-6,
c    LC: QA351.C45.
c
c  Parameters:
c
c    Input, double precision X, the argument.
c
c    Output, double precision SH0, the value of H0(x).
c
      implicit none

      double precision a0
      double precision by0
      integer k
      integer km
      double precision p0
      double precision pi
      double precision q0
      double precision r
      double precision s
      double precision sh0
      double precision t
      double precision t2
      double precision ta0
      double precision x

      pi = 3.141592653589793D+00
      s = 1.0D+00
      r = 1.0D+00

      if ( x .le. 20.0D+00 ) then
        a0 = 2.0D+00 * x / pi
        do k = 1, 60
          r = - r * x / ( 2.0D+00 * k + 1.0D+00 ) * x 
     &      / ( 2.0D+00 * k + 1.0D+00 )
          s = s + r
          if ( abs ( r ) .lt. abs ( s ) * 1.0D-12 ) then
            go to 10
          end if
        end do

10      continue

        sh0 = a0 * s

      else

        if ( x .lt. 50.0D+00 ) then
          km = int ( 0.5D+00 * ( x + 1.0D+00 ) )
        else
          km = 25
        end if

        do k = 1, km
          r = - r * ( ( 2.0D+00 * k - 1.0D+00 ) / x ) ** 2
          s = s + r
          if ( abs ( r ) .lt. abs ( s ) * 1.0D-12 ) then
            go to 20
          end if
        end do

20      continue

        t = 4.0D+00 / x
        t2 = t * t

        p0 = ((((
     &    - 0.37043D-05     * t2 
     &    + 0.173565D-04 )  * t2 
     &    - 0.487613D-04 )  * t2 
     &    + 0.17343D-03 )   * t2
     &    - 0.1753062D-02 ) * t2 
     &    + 0.3989422793D+00

        q0 = t * (((((
     &      0.32312D-05     * t2
     &    - 0.142078D-04 )  * t2
     &    + 0.342468D-04 )  * t2 
     &    - 0.869791D-04 )  * t2
     &    + 0.4564324D-03 ) * t2
     &    - 0.0124669441D+00 )

        ta0 = x - 0.25D+00 * pi
        by0 = 2.0D+00 / sqrt ( x ) 
     &    * ( p0 * sin ( ta0 ) + q0 * cos ( ta0 ) )
        sh0 = 2.0D+00 / ( pi * x ) * s + by0

      end if

      return
      end
      subroutine stvh1 ( x, sh1 )

c*********************************************************************72
c
cc STVH1 computes the Struve function H1(x).
c
c  Licensing:
c
c    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
c    they give permission to incorporate this routine into a user program 
c    provided that the copyright is acknowledged.
c
c  Modified:
c
c    22 July 2012
c
c  Author:
c
c    Shanjie Zhang, Jianming Jin
c
c  Reference:
c
c    Shanjie Zhang, Jianming Jin,
c    Computation of Special Functions,
c    Wiley, 1996,
c    ISBN: 0-471-11963-6,
c    LC: QA351.C45.
c
c  Parameters:
c
c    Input, double precision X, the argument.
c
c    Output, double precision SH1, the value of H1(x).
c
      implicit none

      double precision a0
      double precision by1
      integer k
      integer km
      double precision p1
      double precision pi
      double precision q1
      double precision r
      double precision s
      double precision sh1
      double precision t
      double precision t2
      double precision ta1
      double precision x

      pi = 3.141592653589793D+00
      r = 1.0D+00

      if ( x .le. 20.0D+00 ) then

        s = 0.0D+00
        a0 = - 2.0D+00 / pi
        do k = 1, 60
          r = - r * x * x / ( 4.0D+00 * k * k - 1.0D+00 )
          s = s + r
          if ( abs ( r ) .lt. abs ( s ) * 1.0D-12 ) then
            go to 10
          end if
        end do

10      continue

        sh1 = a0 * s

      else

        s = 1.0D+00

        if ( x .le. 50.0D+00 ) then
          km = int ( 0.5D+00 * x )
        else
          km = 25
        end if

        do k = 1, km
          r = - r * ( 4.0D+00 * k * k - 1.0D+00 ) / ( x * x )
          s = s + r
          if ( abs ( r ) .lt. abs ( s ) * 1.0D-12 ) then
            go to 20
          end if
        end do

20      continue

        t = 4.0D+00 / x
        t2 = t * t

        p1 = (((( 
     &      0.42414D-05      * t2
     &    - 0.20092d-04 )    * t2
     &    + 0.580759D-04 )   * t2
     &    - 0.223203D-03 )   * t2
     &    + 0.29218256D-02 ) * t2
     &    + 0.3989422819D+00

        q1 = t * (((((
     &    - 0.36594D-05     * t2 
     &    + 0.1622D-04 )    * t2
     &    - 0.398708D-04 )  * t2
     &    + 0.1064741D-03 ) * t2
     &    - 0.63904D-03 )   * t2
     &    + 0.0374008364D+00 )

        ta1 = x - 0.75D+00 * pi
        by1 = 2.0D+00 / sqrt ( x ) * ( p1 * sin ( ta1 ) 
     &    + q1 * cos ( ta1 ) )
        sh1 = 2.0D+00 / pi * ( 1.0D+00 + s / ( x * x ) ) + by1

      end if

      return
      end
      subroutine stvhv ( v, x, hv )

c*********************************************************************72
c
cc STVHV computes the Struve function Hv(x) with arbitrary order v.
c
c  Licensing:
c
c    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
c    they give permission to incorporate this routine into a user program 
c    provided that the copyright is acknowledged.
c
c  Modified:
c
c    24 July 2012
c
c  Author:
c
c    Shanjie Zhang, Jianming Jin
c
c  Reference:
c
c    Shanjie Zhang, Jianming Jin,
c    Computation of Special Functions,
c    Wiley, 1996,
c    ISBN: 0-471-11963-6,
c    LC: QA351.C45.
c
c  Parameters:
c
c    Input, double precision V, the order of the function.
c
c    Input, double precision X, the argument.
c
c    Output, double precision HV, the value of Hv(x).
c
      implicit none

      double precision bf
      double precision bf0
      double precision bf1
      double precision by0
      double precision by1
      double precision byv
      double precision ga
      double precision gb
      double precision hv
      integer k
      integer l
      integer n
      double precision pi
      double precision pu0
      double precision pu1
      double precision qu0
      double precision qu1
      double precision r1
      double precision r2
      double precision s
      double precision s0
      double precision sa
      double precision sr
      double precision t0
      double precision t1
      double precision u
      double precision u0
      double precision v
      double precision v0
      double precision va
      double precision vb
      double precision vt
      double precision x

      pi = 3.141592653589793D+00

      if ( x .eq. 0.0D+00 ) then
        if ( -1.0D+00 .lt. v .or. int ( v ) - v .eq. 0.5D+00 ) then
          hv = 0.0D+00
        else if ( v .lt. -1.0D+00 ) then
          hv = ( -1.0D+00 ) ** ( int ( 0.5D+00 - v ) - 1 ) * 1.0D+300
        else if ( v .eq. -1.0D+00 ) then
          hv = 2.0D+00 / pi
        end if
        return
      end if

      if ( x .le. 20.0D+00 ) then

        v0 = v + 1.5D+00
        call gamma ( v0, ga )
        s = 2.0D+00 / ( sqrt ( pi ) * ga )
        r1 = 1.0D+00

        do k = 1, 100
          va = k + 1.5D+00
          call gamma ( va, ga )
          vb = v + k + 1.5D+00
          call gamma ( vb, gb )
          r1 = -r1 * ( 0.5D+00 * x ) ** 2
          r2 = r1 / ( ga * gb )
          s = s + r2
          if ( abs ( r2 ) .lt. abs ( s ) * 1.0D-12 ) then
            go to 10
          end if
        end do

10      continue

        hv = ( 0.5D+00 * x ) ** ( v + 1.0D+00 ) * s

      else

        sa = ( 0.5D+00 * x ) ** ( v - 1.0D+00 ) / pi
        v0 = v + 0.5D+00
        call gamma ( v0, ga )
        s = sqrt ( pi ) / ga
        r1 = 1.0D+00

        do k = 1, 12
          va = k + 0.5D+00
          call gamma ( va, ga )
          vb = - k + v + 0.5D+00
          call gamma ( vb, gb )
          r1 = r1 / ( 0.5D+00 * x ) ** 2
          s = s + r1 * ga / gb
        end do

        s0 = sa * s
        u = abs ( v )
        n = int ( u )
        u0 = u - n

        do l = 0, 1

          vt = 4.0D+00 * ( u0 + l ) ** 2
          r1 = 1.0D+00
          pu1 = 1.0D+00
          do k = 1, 12
            r1 = -0.0078125D+00 * r1 
     &        * ( vt - ( 4.0D+00 * k - 3.0D+00 ) ** 2 ) 
     &        * ( vt - ( 4.0D+00 * k - 1.0D+00 ) ** 2 )
     &        / ( ( 2.0D+00 * k - 1.0D+00 ) * k * x * x )
            pu1 = pu1 + r1
          end do

          qu1 = 1.0D+00
          r2 = 1.0D+00
          do k = 1, 12
            r2 = -0.0078125D+00 * r2
     &        * ( vt - ( 4.0D+00 * k - 1.0D+00 ) ** 2 )
     &        * ( vt - ( 4.0D+00 * k + 1.0D+00 ) ** 2 )
     &        / ( ( 2.0D+00 * k + 1.0D+00 ) * k * x * x )
            qu1 = qu1 + r2
          end do
          qu1 = 0.125D+00 * ( vt - 1.0D+00 ) / x * qu1

          if ( l .eq. 0 ) then
            pu0 = pu1
            qu0 = qu1
          end if

        end do

        t0 = x - ( 0.5D+00 * u0 + 0.25D+00 ) * pi
        t1 = x - ( 0.5D+00 * u0 + 0.75D+00 ) * pi
        sr = sqrt ( 2.0D+00 / ( pi * x ) )
        by0 = sr * ( pu0 * sin ( t0 ) + qu0 * cos ( t0 ) )
        by1 = sr * ( pu1 * sin ( t1 ) + qu1 * cos ( t1 ) )
        bf0 = by0
        bf1 = by1
        do k = 2, n
          bf = 2.0D+00 * ( k - 1.0D+00 + u0 ) / x * bf1 - bf0
          bf0 = bf1
          bf1 = bf
        end do

        if ( n .eq. 0 ) then
          byv = by0
        else if ( n .eq. 1 ) then
          byv = by1
        else
          byv = bf
        end if
        hv = byv + s0
      end if

      return
      end
      subroutine stvl0 ( x, sl0 )

c*********************************************************************72
c
cc STVL0 computes the modified Struve function L0(x).
c
c  Licensing:
c
c    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
c    they give permission to incorporate this routine into a user program 
c    provided that the copyright is acknowledged.
c
c  Modified:
c
c    22 July 2012
c
c  Author:
c
c    Shanjie Zhang, Jianming Jin
c
c  Reference:
c
c    Shanjie Zhang, Jianming Jin,
c    Computation of Special Functions,
c    Wiley, 1996,
c    ISBN: 0-471-11963-6,
c    LC: QA351.C45.
c
c  Parameters:
c
c    Input, double precision X, the argument.
c
c    Output, double precision SL0, the function value.
c
      implicit none

      double precision a0
      double precision a1
      double precision bi0
      integer k
      integer km
      double precision pi
      double precision r
      double precision s
      double precision sl0
      double precision x

      pi = 3.141592653589793D+00
      s = 1.0D+00
      r = 1.0D+00

      if ( x .le. 20.0D+00 ) then

        a0 = 2.0D+00 * x / pi

        do k = 1, 60
          r = r * ( x / ( 2.0D+00 * k + 1.0D+00 ) ) ** 2
          s = s + r
          if ( abs ( r / s ) .lt. 1.0D-12 ) then
            go to 10
          end if
        end do

10      continue

        sl0 = a0 * s

      else

        if ( x .lt. 50.0D+00 ) then
          km = int ( 0.5D+00 * ( x + 1.0D+00 ) )
        else
          km = 25
        end if

        do k = 1, km
          r = r * ( ( 2.0D+00 * k - 1.0D+00 ) / x ) ** 2
          s = s + r
          if ( abs ( r / s ) .lt. 1.0D-12 ) then
            go to 20
          end if
        end do

20      continue

        a1 = exp ( x ) / sqrt ( 2.0D+00 * pi * x )
        r = 1.0D+00
        bi0 = 1.0D+00
        do k = 1, 16
          r = 0.125D+00 * r * ( 2.0D+00 * k - 1.0D+00 ) ** 2 / ( k * x )
          bi0 = bi0 + r
          if ( abs ( r / bi0 ) .lt. 1.0D-12 ) then
            go to 30
          end if
        end do

30      continue

        bi0 = a1 * bi0
        sl0 = - 2.0D+00 / ( pi * x ) * s + bi0

      end if

      return
      end
      subroutine stvl1 ( x, sl1 )

c*********************************************************************72
c
cc STVL1 computes the modified Struve function L1(x).
c
c  Licensing:
c
c    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
c    they give permission to incorporate this routine into a user program 
c    provided that the copyright is acknowledged.
c
c  Modified:
c
c    24 July 2012
c
c  Author:
c
c    Shanjie Zhang, Jianming Jin
c
c  Reference:
c
c    Shanjie Zhang, Jianming Jin,
c    Computation of Special Functions,
c    Wiley, 1996,
c    ISBN: 0-471-11963-6,
c    LC: QA351.C45.
c
c  Parameters:
c
c    Input, double precision X, the argument.
c
c    Output, double precision SL1, the function value.
c
      implicit none

      double precision a1
      double precision bi1
      integer k
      integer km
      double precision pi
      double precision r
      double precision s
      double precision sl1
      double precision x

      pi = 3.141592653589793D+00
      r = 1.0D+00

      if ( x .le. 20.0D+00 ) then
        s = 0.0D+00
        do k = 1, 60
          r = r * x * x / ( 4.0D+00 * k * k - 1.0D+00 )
          s = s + r
          if ( abs ( r / s ) .lt. 1.0D-12 ) then
            go to 10
          end if
        end do

10      continue

        sl1 = 2.0D+00 / pi * s

      else

        s = 1.0D+00
        km = int ( 0.50D+00 * x )
        km = min ( km, 25 )

        do k = 1, km
          r = r * ( 2.0D+00 * k + 3.0D+00 ) 
     &      * ( 2.0D+00 * k + 1.0D+00 ) / ( x * x )
          s = s + r
          if ( abs ( r / s ) .lt. 1.0D-12 ) then
            go to 20
          end if
        end do

20      continue

        sl1 = 2.0D+00 / pi * ( -1.0D+00 + 1.0D+00 
     &    / ( x * x ) + 3.0D+00 * s / x**4 )
        a1 = exp ( x ) / sqrt ( 2.0D+00 * pi * x )
        r = 1.0D+00
        bi1 = 1.0D+00
        do k = 1, 16
          r = -0.125D+00 * r 
     &      * ( 4.0D+00 - ( 2.0D+00 * k - 1.0D+00 ) ** 2 ) / ( k * x )
          bi1 = bi1 + r
          if ( abs ( r / bi1 ) .lt. 1.0D-12 ) then
            go to 30
          end if
        end do

30      continue

        sl1 = sl1 + a1 * bi1

      end if

      return
      end
      subroutine stvlv ( v, x, slv )

c*********************************************************************72
c
cc STVLV computes the modified Struve function Lv(x) with arbitary order.
c
c  Licensing:
c
c    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
c    they give permission to incorporate this routine into a user program 
c    provided that the copyright is acknowledged.
c
c  Modified:
c
c    04 July 2012
c
c  Author:
c
c    Shanjie Zhang, Jianming Jin
c
c  Reference:
c
c    Shanjie Zhang, Jianming Jin,
c    Computation of Special Functions,
c    Wiley, 1996,
c    ISBN: 0-471-11963-6,
c    LC: QA351.C45.
c
c  Parameters:
c
c    Input, double precision V, the order of Lv(x).
c
c    Input, double precision X, the argument of Lv(x).
c
c    Output, double precision SLV, the value of Lv(x).
c
      implicit none

      double precision bf
      double precision bf0
      double precision bf1
      double precision biv
      double precision biv0
      double precision ga
      double precision gb
      integer k
      integer l
      integer n
      double precision pi
      double precision r
      double precision r1
      double precision r2
      double precision s
      double precision s0
      double precision sa
      double precision slv
      double precision u
      double precision u0
      double precision v
      double precision v0
      double precision va
      double precision vb
      double precision vt
      double precision x

      pi = 3.141592653589793D+00

      if ( x .eq. 0.0D+00 ) then

        if ( -1.0D+00 .lt. v .or. int ( v ) - v .eq. 0.5D+00 ) then
          slv = 0.0D+00
        else if ( v .lt. -1.0D+00 ) then
          slv = ( -1 ) ** ( int ( 0.5D+00 - v ) - 1 ) * 1.0D+300
        else if ( v .eq. -1.0D+00 ) then
          slv = 2.0D+00 / pi
        end if

      else if ( x .le. 40.0D+00 ) then

        v0 = v + 1.5D+00
        call gamma ( v0, ga )
        s = 2.0D+00 / ( sqrt ( pi ) * ga )
        r1 = 1.0D+00
        do k = 1, 100
          va = k + 1.5D+00
          call gamma ( va, ga )
          vb = v + k + 1.5D+00
          call gamma ( vb, gb )
          r1 = r1 * ( 0.5D+00 * x ) ** 2
          r2 = r1 / ( ga * gb )
          s = s + r2
          if ( abs ( r2 / s ) .lt. 1.0D-12 ) then
            go to 10
          end if
        end do

10      continue

        slv = ( 0.5D+00 * x ) ** ( v + 1.0D+00 ) * s

      else

        sa = -1.0D+00 / pi * ( 0.5D+00 * x ) ** ( v - 1.0D+00 )
        v0 = v + 0.5D+00
        call gamma ( v0, ga )
        s = - sqrt ( pi ) / ga
        r1 = -1.0D+00
        do k = 1, 12
          va = k + 0.5D+00
          call gamma ( va, ga )
          vb = - k + v + 0.5D+00
          call gamma ( vb, gb )
          r1 = - r1 / ( 0.5D+00 * x ) ** 2
          s = s + r1 * ga / gb
        end do
        s0 = sa * s
        u = abs ( v )
        n = int ( u )
        u0 = u - n
        do l = 0, 1
          vt = u0 + l
          r = 1.0D+00
          biv = 1.0D+00
          do k = 1, 16
            r = -0.125D+00 * r * ( 4.0D+00 * vt * vt - 
     &        ( 2.0D+00 * k - 1.0D+00 )**2 ) / ( k * x )
            biv = biv + r
            if ( abs ( r / biv ) .lt. 1.0D-12 ) then
              go to 20
            end if
          end do

20        continue

          if ( l .eq. 0 ) then
            biv0 = biv
          end if

        end do

        bf0 = biv0
        bf1 = biv
        do k = 2, n
          bf = - 2.0D+00 * ( k - 1.0D+00 + u0 ) / x * bf1 + bf0
          bf0 = bf1
          bf1 = bf
        end do

        if ( n .eq. 0 ) then
          biv = biv0
        else if ( 1 .lt. n ) then
          biv = bf
        end if

        slv = exp ( x ) / sqrt ( 2.0D+00 * pi * x ) * biv + s0

      end if

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
      subroutine vvla ( va, x, pv )

c*********************************************************************72
c
cc VVLA computes parabolic cylinder function Vv(x) for large arguments.
c
c  Licensing:
c
c    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
c    they give permission to incorporate this routine into a user program 
c    provided that the copyright is acknowledged.
c
c  Modified:
c
c    06 April 2012
c
c  Author:
c
c    Shanjie Zhang, Jianming Jin
c
c  Reference:
c
c    Shanjie Zhang, Jianming Jin,
c    Computation of Special Functions,
c    Wiley, 1996,
c    ISBN: 0-471-11963-6,
c    LC: QA351.C45.
c
c  Parameters:
c
c    Input, double precision X, the argument.
c
c    Input, double precision VA, the order nu.
c
c    Output, double precision PV, the value of V(nu,x).
c
      implicit none

      double precision a0
      double precision dsl
      double precision eps
      double precision gl
      integer k
      double precision pdl
      double precision pi
      double precision pv
      double precision qe
      double precision r
      double precision va
      double precision x
      double precision x1

      pi = 3.141592653589793D+00
      eps = 1.0D-12
      qe = exp ( 0.25D+00 * x * x )
      a0 = abs ( x ) ** ( -va - 1.0D+00 ) * sqrt ( 2.0D+00 / pi ) * qe

      r = 1.0D+00
      pv = 1.0D+00
      do k = 1, 18
        r = 0.5D+00 * r * ( 2.0D+00 * k + va - 1.0D+00 ) 
     &    * ( 2.0D+00 * k + va ) / ( k * x * x )
        pv = pv + r
        if ( abs ( r / pv ) .lt. eps ) then
          go to 10
        end if
      end do

10    continue

      pv = a0 * pv

      if ( x .lt. 0.0D+00 ) then
        x1 = -x
        call dvla ( va, x1, pdl )
        call gamma ( -va, gl )
        dsl = sin ( pi * va ) * sin ( pi * va )
        pv = dsl * gl / pi * pdl - cos ( pi * va ) * pv
      end if

      return
      end
      subroutine vvsa ( va, x, pv )

c*********************************************************************72
c
cc VVSA computes parabolic cylinder function V(nu,x) for small arguments.
c
c  Licensing:
c
c    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
c    they give permission to incorporate this routine into a user program 
c    provided that the copyright is acknowledged.
c
c  Modified:
c
c    06 April 2012
c
c  Author:
c
c    Shanjie Zhang, Jianming Jin
c
c  Reference:
c
c    Shanjie Zhang, Jianming Jin,
c    Computation of Special Functions,
c    Wiley, 1996,
c    ISBN: 0-471-11963-6,
c    LC: QA351.C45.
c
c  Parameters:
c
c    Input, double precision X, the argument.
c
c    Input, double precision VA, the order nu.
c
c    Output, double precision PV, the value of V(nu,x).
c
      implicit none

      double precision a0
      double precision ep
      double precision eps
      double precision fac
      double precision g1
      double precision ga0
      double precision gm
      double precision gw
      integer m
      double precision pi
      double precision pv
      double precision r
      double precision r1
      double precision sq2
      double precision sv
      double precision sv0
      double precision v1
      double precision va
      double precision va0
      double precision vb0
      double precision vm
      double precision x

      eps = 1.0D-15
      pi = 3.141592653589793D+00
      ep = exp ( -0.25D+00 * x * x )
      va0 = 1.0D+00 + 0.5D+00 * va

      if ( x .eq. 0.0D+00 ) then

        if ( ( va0 .le. 0.0D+00 .and. va0 .eq. int ( va0 ) ) .or.
     &    va .eq. 0.0D+00 ) then
          pv = 0.0D+00
        else
          vb0 = -0.5D+00 * va
          sv0 = sin ( va0 * pi )
          call gamma ( va0, ga0 )
          pv = 2.0D+00 ** vb0 * sv0 / ga0
        end if

      else

        sq2 = sqrt ( 2.0D+00 )
        a0 = 2.0D+00 ** ( -0.5D+00 * va ) * ep / ( 2.0D+00 * pi )
        sv = sin ( - ( va + 0.5D+00 ) * pi )
        v1 = -0.5D+00 * va
        call gamma ( v1, g1 )
        pv = ( sv + 1.0D+00 ) * g1
        r = 1.0D+00
        fac = 1.0D+00

        do m = 1, 250
          vm = 0.5D+00 * ( m - va )
          call gamma ( vm, gm )
          r = r * sq2 * x / m
          fac = - fac
          gw = fac * sv + 1.0D+00
          r1 = gw * r * gm
          pv = pv + r1
          if ( abs ( r1 / pv ) .lt. eps .and. gw .ne. 0.0D+00 ) then
            go to 10
          end if
        end do

10      continue

        pv = a0 * pv

      end if

      return
      end
