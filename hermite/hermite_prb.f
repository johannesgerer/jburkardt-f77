      program main

c*********************************************************************72
c
cc MAIN is the main program for HERMITE_PRB.
c
c  Discussion:
c
c    HERMITE_PRB tests HERMITE.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    01 November 2011
c
c  Author:
c
c    John Burkardt
c
      implicit none

      call timestamp ( )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'HERMITE_PRB'
      write ( *, '(a)' ) '  FORTRAN77 version'
      write ( *, '(a)' ) '  Test the HERMITE library.'

      call test01 ( )
      call test02 ( )
      call test03 ( )
      call test04 ( )
      call test05 ( )
      call test06 ( )
      call test07 ( )
      call test08 ( )
c
c  Terminate.
c
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'HERMITE_PRB'
      write ( *, '(a)' ) '  Normal end of execution.'

      write ( *, '(a)' ) ' '
      call timestamp ( )

      stop
      end
      subroutine test01 ( )

c*********************************************************************72
c
cc TEST01 uses f(x) = 1 + 2x + 3x^2 at x = 0, 1, 2.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    19 May 2011
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n 
      parameter ( n = 3 )

      double precision x(n)
      double precision y(n)
      double precision yp(n)

      save x
      save y
      save yp

      data x / 0.0D+00, 1.0D+00,  2.0D+00 /
      data y / 1.0D+00, 6.0D+00, 17.0D+00 /
      data yp / 2.0D+00, 8.0D+00, 14.0D+00 /

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST01'
      write ( *, '(a)' ) 
     &  '  HERMITE computes the Hermite interpolant to data.'
      write ( *, '(a)' ) '  Here, f(x) = 1 + 2x + 3x^2.'

      call hermite_demo ( n, x, y, yp )

      return
      end
      subroutine test02 ( )

c*********************************************************************72
c
cc TEST02 uses f(x) = 6 + 5x + 4x^2 + 3x^3 + 2x^4 + x^5 at x = 0, 1, 2.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    19 May 2011
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n
      parameter ( n = 3 )

      integer i
      double precision x(n)
      double precision y(n)
      double precision yp(n)

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST02'
      write ( *, '(a)' ) 
     &  '  HERMITE computes the Hermite interpolant to data.'
      write ( *, '(a)' ) 
     &  '  Here, f(x) = 6 + 5x + 4x^2 + 3x^3 + 2x^4 + x^5.'

      do i = 1, 3

        x(i) = dble ( i - 1 )

        y(i) = 6.0D+00 + x(i) * ( 
     &         5.0D+00 + x(i) * ( 
     &         4.0D+00 + x(i) * ( 
     &         3.0D+00 + x(i) * ( 
     &         2.0D+00 + x(i) ) ) ) )

        yp(i) = 5.0D+00 + x(i) * ( 
     &          8.0D+00 + x(i) * ( 
     &          9.0D+00 + x(i) * ( 
     &          8.0D+00 + x(i) *   
     &          5.0D+00 ) ) )

      end do

      call hermite_demo ( n, x, y, yp )

      return
      end
      subroutine test03 ( )

c*********************************************************************72
c
cc TEST03 uses f(x) = r1 + r2x + r3x^2 + r4x^3 + r5x^4 + r6x^5 at x = r7 r8 r9
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    19 May 2011
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n
      parameter ( n = 3 )

      double precision c(0:2*n-1)
      integer i
      integer seed
      double precision x(n)
      double precision y(n)
      double precision yp(n)

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST03'
      write ( *, '(a)' ) 
     &  '  HERMITE computes the Hermite interpolant to data.'
      write ( *, '(a)' ) 
     &  '  Here, f(x) is a fifth order polynomial with random'
      write ( *, '(a)' ) 
     &  '  coefficients, and the abscissas are random.'

      seed = 123456789

      call r8vec_uniform_01 ( n, seed, x )
      call r8vec_print ( n, x, '  Random abscissas' )

      call r8vec_uniform_01 ( 2 * n, seed, c )
      call r8vec_print ( 2 * n, c, '  Random polynomial coefficients.' )

      do i = 1, 3

        y(i) = c(0) + x(i) * ( 
     &         c(1) + x(i) * ( 
     &         c(2) + x(i) * ( 
     &         c(3) + x(i) * ( 
     &         c(4) + x(i) * ( 
     &         c(5) ) ) ) ) )

        yp(i) = c(1)           + x(i) * ( 
     &          c(2) * 2.0D+00 + x(i) * ( 
     &          c(3) * 3.0D+00 + x(i) * ( 
     &          c(4) * 4.0D+00 + x(i) *   
     &          c(5) * 5.0D+00 ) ) )
      end do

      call hermite_demo ( n, x, y, yp )

      return
      end
      subroutine test04 ( )

c*********************************************************************72
c
cc TEST04 interpolates the Runge function with equally spaced data.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    31 October 2011
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer m
      parameter ( m = 15 )
      integer md
      parameter ( md = m * 2 )
      integer ms
      parameter ( ms = m * 10 )

      integer i
      double precision max_dif
      integer n
      integer nd
      integer ns
      double precision x(m)
      double precision xd(md)
      double precision xdp(md-1)
      double precision xs(ms)
      double precision xhi
      double precision xlo
      double precision xt
      double precision y(m)
      double precision yd(md)
      double precision ydp(md-1)
      double precision yp(m)
      double precision ys(ms)
      double precision yt

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST04'
      write ( *, '(a)' ) 
     &  '  HERMITE computes the Hermite interpolant to data.'
      write ( *, '(a)' ) '  Here, f(x) is the Runge function'
      write ( *, '(a)' ) 
     &  '  and the data is evaluated at equally spaced points.'
      write ( *, '(a)' ) '  As N increases, the maximum error grows.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '     N     Max | F(X) - H(F(X)) |'
      write ( *, '(a)' ) ' '

      do n = 3, 15, 2

        nd = 2 * n
        ns = 10 * ( n - 1 ) + 1

        xlo = -5.0D+00
        xhi = +5.0D+00
        call r8vec_linspace ( n, xlo, xhi, x )

        do i = 1, n
          y(i) = 1.0D+00 / ( 1.0D+00 + x(i)**2 )
          yp(i) = - 2.0D+00 * x(i) / ( 1.0D+00 + x(i)**2 )**2
        end do

        call hermite_interpolant ( n, x, y, yp, xd, yd, xdp, ydp )
c
c  Compare exact and interpolant at sample points.
c
        call r8vec_linspace ( ns, xlo, xhi, xs )

        call dif_vals ( nd, xd, yd, ns, xs, ys )

        max_dif = 0.0D+00
        do i = 1, ns
          xt = xs(i)
          yt = 1.0D+00 / ( 1.0D+00 + xt * xt )
          max_dif = max ( max_dif, abs ( ys(i) - yt ) )
        end do

        write ( *, '(2x,i4,2x,g14.6)' ) n, max_dif

      end do

      return
      end
      subroutine test05 ( )

c*********************************************************************72
c
cc TEST05 interpolates the Runge function with Chebyshev spaced data.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    31 October 2011
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer m
      parameter ( m = 15 )
      integer md
      parameter ( md = m * 2 )
      integer ms
      parameter ( ms = m * 10 )

      integer i
      double precision max_dif
      integer n
      integer nd
      integer ns
      double precision x(m)
      double precision xd(md)
      double precision xdp(md-1)
      double precision xs(ms)
      double precision xhi
      double precision xlo
      double precision xt
      double precision y(m)
      double precision yd(md)
      double precision ydp(md-1)
      double precision yp(m)
      double precision ys(ms)
      double precision yt

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST05'
      write ( *, '(a)' ) 
     &  '  HERMITE computes the Hermite interpolant to data.'
      write ( *, '(a)' ) '  Here, f(x) is the Runge function'
      write ( *, '(a)' ) 
     &  '  and the data is evaluated at Chebyshev spaced points.'
      write ( *, '(a)' ) 
     &  '  As N increases, the maximum error decreases.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '     N     Max | F(X) - H(F(X)) |'
      write ( *, '(a)' ) ' '

      do n = 3, 15, 2

        nd = 2 * n
        ns = 10 * ( n - 1 ) + 1

        xlo = -5.0D+00
        xhi = +5.0D+00
        call r8vec_chebyshev ( n, xlo, xhi, x )

        do i = 1, n
          y(i) = 1.0D+00 / ( 1.0D+00 + x(i)**2 )
          yp(i) = - 2.0D+00 * x(i) / ( 1.0D+00 + x(i)**2 )**2
        end do

        call hermite_interpolant ( n, x, y, yp, xd, yd, xdp, ydp )
c
c  Compare exact and interpolant at sample points.
c
        call r8vec_linspace ( ns, xlo, xhi, xs )

        call dif_vals ( nd, xd, yd, ns, xs, ys )

        max_dif = 0.0D+00
        do i = 1, ns
          xt = xs(i)
          yt = 1.0D+00 / ( 1.0D+00 + xt * xt )
          max_dif = max ( max_dif, abs ( ys(i) - yt ) )
        end do

        write ( *, '(2x,i4,2x,g14.6)' ) n, max_dif

      end do

      return
      end
      subroutine test06 ( )

c*********************************************************************72
c
cc TEST06 tests HERMITE_BASIS_0 and HERMITE_BASIS_1.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    28 May 2011
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer nd
      parameter ( nd = 2 )

      double precision f01
      double precision f02
      double precision f11
      double precision f12
      integer i
      double precision xd(nd)
      double precision xv
      double precision yd(nd)
      double precision yh
      double precision ypd(nd)
      double precision yv

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST06:'
      write ( *, '(a)' ) 
     &  '  HERMITE_BASIS_0 and HERMITE_BASIS_1 evaluate the'
      write ( *, '(a)' ) 
     &  '  Hermite global polynomial basis functions'
      write ( *, '(a)' ) 
     &  '  of type 0: associated with function values, and'
      write ( *, '(a)' ) 
     &  '  of type 1: associated with derivative values.'
c
c  Let y = x^3 + x^2 + x + 1,
c  and compute the Hermite global polynomial interpolant based on two 
c  abscissas:
c
      xd(1) = 0.0D+00
      xd(2) = 10.0D+00

      do i = 1, nd
        yd(i) = xd(i)**3 + xd(i)**2 + xd(i) + 1.0D+00
        ypd(i) = 3.0D+00 * xd(i)**2 + 2.0D+00 * xd(i) + 1.0D+00
      end do

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Interpolate y = x^3 + x^2 + x + 1.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '     XD         Y(XD)      Y''(XD)'
      write ( *, '(a)' ) ' '
      do i = 1, nd
        write ( *, '(2x,g10.4,2x,g10.4,2x,g10.4)' ) xd(i), yd(i), ypd(i)
      end do

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '     XV         Y(XV)      H(XV)'
      write ( *, '(a)' ) ' '

      do i = 1, 11

        xv = dble ( i - 1 )

        yv = xv**3 + xv**2 + xv + 1.0D+00

        call hermite_basis_0 ( 2, xd, 1, xv, f01 )
        call hermite_basis_1 ( 2, xd, 1, xv, f11 )
        call hermite_basis_0 ( 2, xd, 2, xv, f02 )
        call hermite_basis_1 ( 2, xd, 2, xv, f12 )

        yh = yd(1) * f01 + ypd(1) * f11 + yd(2) * f02 + ypd(2) * f12

        write ( *, '(2x,g10.4,2x,g10.4,2x,g10.4)' ), xv, yv, yh

      end do

      return
      end
      subroutine test07 ( )

c*********************************************************************72
c
cc TEST07 tests HERMITE_INTERPOLANT_RULE.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    16 June 2011
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n_max
      parameter ( n_max = 11 )

      double precision a
      double precision b
      integer i
      integer k
      integer n
      double precision q
      double precision w(2*n_max)
      double precision x(n_max)

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST07:'
      write ( *, '(a)' ) '  HERMITE_INTERPOLANT_RULE'
      write ( *, '(a)' ) 
     &  '  is given a set of N abscissas for a Hermite interpolant'
      write ( *, '(a)' ) '  and returns N pairs of quadrature weights'
      write ( *, '(a)' ) 
     &  '  for function and derivative values at the abscissas.'

      n = 3
      a = 0.0D+00
      b = 10.0D+00
      call r8vec_linspace ( n, a, b, x )
      call hermite_interpolant_rule ( n, a, b, x, w )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 
     &  '     I       X               W(F(X))        W(F''(X))'
      write ( *, '(a)' ) ' '
      k = 1
      do i = 1, n
        write ( *, '(2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) 
     &    i, x(i), w(k), w(k+1)
        k = k + 2
      end do

      write ( *, '(a)' ) ' '
      write ( *, '(a,g14.6,a,g14.6)' ) 
     &  '  Use the quadrature rule over interval ', a, ' to ' , b
      write ( *, '(a)' ) ' '

      q = 0.0D+00
      k = 1
      do i = 1, n
        q = q + w(k) * 1 + w(k+1) * 0.0D+00
        k = k + 2
      end do
      write ( *, * ) '  Estimate integral of 1 = ', q

      q = 0.0D+00
      k = 1
      do i = 1, n
        q = q + w(k) * x(i) + w(k+1) * 1.0D+00
        k = k + 2
      end do
      write ( *, * ) '  Estimate integral of X = ', q

      q = 0.0D+00
      k = 1
      do i = 1, n
        q = q + w(k) * x(i)**2 + w(k+1) * 2.0D+00 * x(i)
        k = k + 2
      end do
      write ( *, * ) '  Estimate integral of X^2 = ', q

      q = 0.0D+00
      k = 1
      do i = 1, n
        q = q + w(k) * sin ( x(i) ) + w(k+1) * cos ( x(i) )
        k = k + 2
      end do
      write ( *, * ) '  Estimate integral of SIN(X) = ', q

      n = 3
      a = 0.0D+00
      b = 1.0D+00
      call r8vec_linspace ( n, a, b, x )
      call hermite_interpolant_rule ( n, a, b, x, w )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 
     &  '     I       X               W(F(X))        W(F''(X))'
      write ( *, '(a)' ) ' '
      k = 1
      do i = 1, n
        write ( *, '(2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) 
     &    i, x(i), w(k), w(k+1)
        k = k + 2
      end do

      write ( *, '(a)' ) ' '
      write ( *, '(a,g14.6,a,g14.6)' ) 
     &  '  Use the quadrature rule over interval ', a, ' to ' , b
      write ( *, '(a)' ) ' '

      q = 0.0D+00
      k = 1
      do i = 1, n
        q = q + w(k) * 1 + w(k+1) * 0.0D+00
        k = k + 2
      end do
      write ( *, * ) '  Estimate integral of 1 = ', q

      q = 0.0D+00
      k = 1
      do i = 1, n
        q = q + w(k) * x(i) + w(k+1) * 1.0D+00
        k = k + 2
      end do
      write ( *, * ) '  Estimate integral of X = ', q

      q = 0.0D+00
      k = 1
      do i = 1, n
        q = q + w(k) * x(i)**2 + w(k+1) * 2.0D+00 * x(i)
        k = k + 2
      end do
      write ( *, * ) '  Estimate integral of X^2 = ', q

      q = 0.0D+00
      k = 1
      do i = 1, n
        q = q + w(k) * sin ( x(i) ) + w(k+1) * cos ( x(i) )
        k = k + 2
      end do
      write ( *, * ) '  Estimate integral of SIN(X) = ', q

      n = 11
      a = 0.0D+00
      b = 10.0D+00
      call r8vec_linspace ( n, a, b, x )
      call hermite_interpolant_rule ( n, a, b, x, w )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 
     &  '     I       X               W(F(X))        W(F''(X))'
      write ( *, '(a)' ) ' '
      k = 1
      do i = 1, n
        write ( *, '(2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) 
     &    i, x(i), w(k), w(k+1)
        k = k + 2
      end do

      write ( *, '(a)' ) ' '
      write ( *, '(a,g14.6,a,g14.6)' ) 
     &  '  Use the quadrature rule over interval ', a, ' to ' , b
      write ( *, '(a)' ) ' '

      q = 0.0D+00
      k = 1
      do i = 1, n
        q = q + w(k) * 1 + w(k+1) * 0.0D+00
        k = k + 2
      end do
      write ( *, * ) '  Estimate integral of 1 = ', q

      q = 0.0D+00
      k = 1
      do i = 1, n
        q = q + w(k) * x(i) + w(k+1) * 1.0D+00
        k = k + 2
      end do
      write ( *, * ) '  Estimate integral of X = ', q

      q = 0.0D+00
      k = 1
      do i = 1, n
        q = q + w(k) * x(i)**2 + w(k+1) * 2.0D+00 * x(i)
        k = k + 2
      end do
      write ( *, * ) '  Estimate integral of X^2 = ', q

      q = 0.0D+00
      k = 1
      do i = 1, n
        q = q + w(k) * sin ( x(i) ) + w(k+1) * cos ( x(i) )
        k = k + 2
      end do
      write ( *, * ) '  Estimate integral of SIN(X) = ', q

      n = 11
      a = 0.0D+00
      b = 1.0D+00
      call r8vec_linspace ( n, a, b, x )
      call hermite_interpolant_rule ( n, a, b, x, w )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 
     &  '     I       X               W(F(X))        W(F''(X))'
      write ( *, '(a)' ) ' '
      k = 1
      do i = 1, n
        write ( *, '(2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) 
     &    i, x(i), w(k), w(k+1)
        k = k + 2
      end do

      write ( *, '(a)' ) ' '
      write ( *, '(a,g14.6,a,g14.6)' ) 
     &  '  Use the quadrature rule over interval ', a, ' to ' , b
      write ( *, '(a)' ) ' '

      q = 0.0D+00
      k = 1
      do i = 1, n
        q = q + w(k) * 1 + w(k+1) * 0.0D+00
        k = k + 2
      end do
      write ( *, * ) '  Estimate integral of 1 = ', q

      q = 0.0D+00
      k = 1
      do i = 1, n
        q = q + w(k) * x(i) + w(k+1) * 1.0D+00
        k = k + 2
      end do
      write ( *, * ) '  Estimate integral of X = ', q

      q = 0.0D+00
      k = 1
      do i = 1, n
        q = q + w(k) * x(i)**2 + w(k+1) * 2.0D+00 * x(i)
        k = k + 2
      end do
      write ( *, * ) '  Estimate integral of X^2 = ', q

      q = 0.0D+00
      k = 1
      do i = 1, n
        q = q + w(k) * sin ( x(i) ) + w(k+1) * cos ( x(i) )
        k = k + 2
      end do
      write ( *, * ) '  Estimate integral of SIN(X) = ', q

      n = 11
      a = 0.0D+00
      b = 1.0D+00
      call r8vec_chebyshev ( n, a, b, x )
      call hermite_interpolant_rule ( n, a, b, x, w )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 
     &  '     I       X               W(F(X))        W(F''(X))'
      write ( *, '(a)' ) ' '
      k = 1
      do i = 1, n
        write ( *, '(2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) 
     &    i, x(i), w(k), w(k+1)
        k = k + 2
      end do

      write ( *, '(a)' ) ' '
      write ( *, '(a,g14.6,a,g14.6)' ) 
     &  '  Use the quadrature rule over interval ', a, ' to ' , b
      write ( *, '(a)' ) ' '

      q = 0.0D+00
      k = 1
      do i = 1, n
        q = q + w(k) * 1 + w(k+1) * 0.0D+00
        k = k + 2
      end do
      write ( *, * ) '  Estimate integral of 1 = ', q

      q = 0.0D+00
      k = 1
      do i = 1, n
        q = q + w(k) * x(i) + w(k+1) * 1.0D+00
        k = k + 2
      end do
      write ( *, * ) '  Estimate integral of X = ', q

      q = 0.0D+00
      k = 1
      do i = 1, n
        q = q + w(k) * x(i)**2 + w(k+1) * 2.0D+00 * x(i)
        k = k + 2
      end do
      write ( *, * ) '  Estimate integral of X^2 = ', q

      q = 0.0D+00
      k = 1
      do i = 1, n
        q = q + w(k) * sin ( x(i) ) + w(k+1) * cos ( x(i) )
        k = k + 2
      end do
      write ( *, * ) '  Estimate integral of SIN(X) = ', q

      return
      end
      subroutine test08 ( )

c*********************************************************************72
c
cc TEST08 tabulates the interpolant and its derivative. 
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    01 November 2011
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n
      parameter ( n = 5 )

      integer nd
      parameter ( nd = 2 * n )
      integer ndp
      parameter ( ndp = 2 * n - 1 )
      integer ns
      parameter ( ns = 4 * ( n - 1 ) + 1 )

      integer i
      double precision x(n)
      double precision xd(nd)
      double precision xdp(ndp)
      double precision xs(ns)
      double precision y(n)
      double precision yd(nd)
      double precision ydp(ndp)
      double precision yp(n)
      double precision ys(ns)
      double precision ysp(ns)

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST08'
      write ( *, '(a)' ) 
     &  '  HERMITE_INTERPOLANT sets up the Hermite interpolant.'
      write ( *, '(a)' ) '  HERMITE_INTERPOLANT_VALUE evaluates it.'
      write ( *, '(a)' ) '  Consider data for y=sin(x) at x=0,1,2,3,4.'

      call r8vec_linspace ( n, 0.0D+00, 4.0D+00, x )

      do i = 1, n
        y(i) = sin ( x(i) )
        yp(i) = cos ( x(i) )
      end do

      call hermite_interpolant ( n, x, y, yp, xd, yd, xdp, ydp )
c
c  Now sample the interpolant at NS points, which include data values.
c
      call r8vec_linspace ( ns, 0.0D+00, 4.0D+00, xs )

      call hermite_interpolant_value ( nd, xd, yd, xdp, ydp, ns, xs, 
     &  ys, ysp )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 
     &  '  In the following table, there should be perfect'
      write ( *, '(a)' ) '  agreement between F and H, and F'' and H'''
      write ( *, '(a)' ) '  at the data points X = 0, 1, 2, 3, and 4.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 
     &  '  In between, H and H'' approximate F and F''.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 
     &  '     I       X(I)          F(X(I))         H(X(I))' // 
     &  '        F''(X(I))        H''(X(I))'
      write ( *, '(a)' ) ' '
      do i = 1, ns
        write ( *, 
     &    '(2x,i4,2x,g14.6,2x,g14.6,2x,g14.6,2x,g14.6,2x,g14.6)' ) 
     &    i, xs(i), sin ( xs(i) ), ys(i), cos ( xs(i) ), ysp(i)
      end do

      return
      end
