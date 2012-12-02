      function besra2 ( v, k )

c*********************************************************************72
c
cc BESRA2 returns the Bessel function ratio I1(V) / I0(V).
c
c  Discussion:
c
c    k decimal digit version of besrat - c1, c2 as in comment cards of be
c
c  Modified:
c
c    30 December 2006
c
c  Author:
c
c    Geoffrey Hill
c
c  Reference:
c
c    Geoffrey Hill,
c    Algorithm 571,
c    Statistics for von Mises's and Fisher's Distributions of Directions,
c    I1(x)/I0(x), I1.5(x)/I0.5(x), and Their Inverses,
c    ACM Transactions on Mathematical Software,
c    Volume 7, Number 2, June 1981, pages 233-238.
c
c  Parameters:
c
c    Input, real V, the concentration parameter of the von Mises
c    distribution.  The sign of V is ignored.
c
c    Input, integer K, ?
c
c    Output, real BESRA2, the expected moduluse of the mean vector
c    sum of unit vectors sampled from the von Mises distribution
c    with concentration parameter V.
c
      implicit none

      real besra2
      real c1
      real c2
      real cx
      real d
      integer j
      integer k
      integer n
      real s
      real v
      real x
      real xx
      real y

      d = real ( k )
      c1 = ( d + 9.0E+00 - 8.0E+00 / d ) * 0.0351E+00
      c2 = ( ( d - 5.0E+00 )**3 / 180.0E+00 + d - 5.0E+00 ) / 10.0E+00
      cx = d * 0.5E+00 + 11.0E+00
      y = 0.0E+00
      x = v

      if ( x .le. cx ) then

        n = int ( ( x + 16.0E+00
     &    - 16.0E+00 / ( x + c1 + 0.75E+00 ) ) * c1 )
        x = x * 0.5E+00
        s = x * x

        do j = 1, n
          y = s / ( real ( n - j + 2 ) + y )
        end do

        besra2 = x / ( 1.0E+00 + y )

      else

        n = int ( ( 68.0E+00 / x + 1.0E+00 ) * c2 ) + 1
        x = x * 4.0E+00
        xx = real ( n * 2 + 1 )

        do j = 1, n
          y = xx / ( ( -2.0E+00 - y ) * xx + x )
          xx = xx - 2.0E+00
        end do

        besra2 = 1.0E+00 - 2.0E+00 / ( x - 1.0E+00 - y )

      end if

      return
      end
      function besrat ( v )

c*********************************************************************72
c
cc BESRAT returns the Bessel function ratio I1(V) / I0(V).
c
c  Discussion:
c
c    This routine returns 
c
c      BESRAT = A(V) 
c
c    where A(V) is the expected modulus of the mean vector 
c    sum of unit vectors sampled from the von Mises distribution of 
c    directions in 2D with concentration parameter = V.
c
c  Modified:
c
c    30 December 2006
c
c  Author:
c
c    Geoffrey Hill
c
c  Reference:
c
c    Geoffrey Hill,
c    Algorithm 571,
c    Statistics for von Mises's and Fisher's Distributions of Directions,
c    I1(x)/I0(x), I1.5(x)/I0.5(x), and Their Inverses,
c    ACM Transactions on Mathematical Software,
c    Volume 7, Number 2, June 1981, pages 233-238.
c
c  Parameters:
c
c    Input, real V, the concentration parameter of the von Mises
c    distribution.  The sign of V is ignored.
c
c    Output, real BESRAT, the expected moduluse of the mean vector
c    sum of unit vectors sampled from the von Mises distribution
c    with concentration parameter V.
c
      implicit none

      real besrat
      real c1
      real c2
      real cx
      integer j
      integer n
      real v
      real x
      real xx
      real y
c
c  adjust to s decimal digit precision by setting data constants -
c     c1 = (s+9.0-8.0/s)*0.0351
c     c2 = ((s-5.0)**3/180.0+s-5.0)/10.0
c     cx = s*0.5 + 11.0
c  for s in range (5,14).  thus for s = 9.3 :
c
      save c1
      save c2
      save cx

      data c1 / 0.613E+00 /
      data c2 / 0.475E+00 /
      data cx / 15.65E+00 /

      y = 0.0E+00
      x = abs ( v )
c
c  For small X, ratio = x/(2+x*x/(4+x*x/(6+x*x/(8+ ... )))
c
      if ( x .le. cx ) then

        n = int ( ( x + 16.0E+00
     &    - 16.0E+00 / ( x + c1 + 0.75E+00 ) ) * c1 )
        x = x * 0.5E+00
        xx = x * x

        do j = 1, n
          y = xx / ( real ( n - j + 2 ) + y )
        end do

        besrat = x / ( 1.0E+00 + y )
c
c  For large X, ratio = 1-2/(4x-1-1/(4x/3-2-1/(4x/5-2- ... )))
c
      else

        n = int ( ( 68.0E+00 / x + 1.0E+00 ) * c2 ) + 1
        x = x * 4.0E+00
        xx = real ( n * 2 + 1 )

        do j = 1, n
          y = xx / ( ( -2.0E+00 - y ) * xx + x )
          xx = xx - 2.0E+00
        end do

        besrat = 1.0E+00 - 2.0E+00 / ( x - 1.0E+00 - y )

      end if

      return
      end
      function r4_cappa3 ( r )

c*********************************************************************72
c
cc R4_CAPPA3 estimates KAPPA from the sample mean, for the Fisher distribution.
c
c  Discussion:
c
c    This routine returns cappa3 = the maximum likelihood estimate of 
c    'kappa', the concentration parameter of the fisher distribution of 
c    directions in 3 dimensions, corresponding to a sample mean vector 
c    modulus r.
c
c    cappa3 = k(b), the inverse function of  b(k) = ratio of modified
c    bessel functions of the first kind of orders 3/2 and 1/2.
c
c    for precision up to 8s (significant decimal digits) omit three
c    lines following statement labeled 20. for greater precision up
c    to 16s the auxiliary subprogram, function R4_SPHERR(X), is needed
c    with precision one decimal digit greater than for cappa3(r).
c    to avoid exponential overflow set bigx = 1.1513*s + 0.4; for 48b
c
c  Modified:
c
c    30 December 2006
c
c  Author:
c
c    Geoffrey Hill
c
c  Reference:
c
c    Geoffrey Hill,
c    Algorithm 571,
c    Statistics for von Mises's and Fisher's Distributions of Directions,
c    I1(x)/I0(x), I1.5(x)/I0.5(x), and Their Inverses,
c    ACM Transactions on Mathematical Software,
c    Volume 7, Number 2, June 1981, pages 233-238.
c
c  Parameters:
c
c    Input, real R, the modulus of the sample mean vector.
c
c    Output, real R4_CAPPA3, the maximum likelihood estimate for the
c    concentration parameter.
c
      implicit none

      real bigx
      real r
      real r4_cappa3
      real r4_spherr
      real s
      real t
      real w1
      real w2
      real x
      real y

      save bigx
      save w1
      save w2

      data bigx / 17.0E+00 /
      data w1 / 10.0E+00 /
      data w2 / 0.01E+00 /

      y = r
      x = -1.0E+00
c
c  Error signal: value -1.0 returned if argument negative or 1.0 or more
c
      if ( y .lt. 0.0E+00 .or. 1.0E+00 .le. y ) then
        r4_cappa3 = x
        return
      end if
c
c  For small R, use inverse Gauss continued fraction. (yields 8.4s)
c
      if ( y .lt. 0.5E+00 ) then

        x = 3.0E+00 * y
        s = x * x
        t = 12.375E+00

        if ( 0.7E+00 .lt. x ) then
          t = ( ( 5.0E+00
     &      * x - 14.74E+00 )
     &      * x + 16.5198E+00 )
     &      * x + 6.2762E+00
        end if

        x = x
     &    / ( 1.0E+00 - s
     &    / ( 15.0E+00 - 4.0E+00 * s
     &    / ( 7.0E+00 - s
     &    / ( 5.0E+00 - s / t ))))
c
c  For large R, approx. inverse of R = coth(k) - 1/k. (yields 8.1s)
c
      else

        x = 1.0E+00 / ( 1.0E+00 - y )

        if ( bigx .lt. x ) then
          r4_cappa3 = x
          return
        end if

        s = 2.0E+00 * x
        t = exp ( s ) - s * s - 1.0E+00

        if ( y .lt. 0.8E+00 ) then
          t = ((( 0.00254E+00
     &      * s - 0.071042E+00 )
     &      * s + 0.6943388E+00 )
     &      * s - 2.3816184E+00 )
     &      * s + 0.1508478E+00 - 0.14789E+00 * y + t
        end if

        x = t * x / ( t + s )

        if ( w1 .lt. x ) then
          r4_cappa3 = x
          return
        end if

      end if
c
c  One stage Newton-Raphson inversion doubles precision.
c
      if ( x .lt. w2 ) then
        r4_cappa3 = x
        return
      end if

      s = exp ( x )
      s = s * 2.0E+00 / ( s * s - 1.0E+00 )
      x = x + ( y - r4_spherr ( x ) ) / ( 1.0E+00 / ( x * x ) - s * s )

      r4_cappa3 = x

      return
      end
      function r4_precis ( a, b )

c*********************************************************************72
c
cc R4_PRECIS returns the number of matching decimal digits.
c
c  Modified:
c
c    30 December 2006
c
c  Author:
c
c    Geoffrey Hill
c
c  Reference:
c
c    Geoffrey Hill,
c    Algorithm 571,
c    Statistics for von Mises's and Fisher's Distributions of Directions,
c    I1(x)/I0(x), I1.5(x)/I0.5(x), and Their Inverses,
c    ACM Transactions on Mathematical Software,
c    Volume 7, Number 2, June 1981, pages 233-238.
c
c  Parameters:
c
c    Input, real A, B, the numbers to be compared.
c
c    Output, real R4_PRECIS, the number of matching decimal digits.
c
      implicit none

      real a
      real b
      real c
      real r4_precis
      real x

      save c

      data c / -0.4342944819E+00 /

      if ( b .eq. 0.0E+00 ) then
        r4_precis = 20.0E+00
        return
      end if

      x = abs ( 1.0E+00 - a / b )

      if ( x .eq. 0.0E+00 ) then
        r4_precis = 20.0E+00
        return
      end if

      r4_precis = alog ( x ) * c

      return
      end
      function r4_spherr ( cappa )

c*********************************************************************72
c
cc R4_SPHERR computes the expected sample mean for the Fisher distribution.
c
c  Discussion:
c
c    This routine returns b(k) for k = abs(cappa).  b(k) is the expected
c    modulus of the mean vector sum of unit vectors sampled from the
c    fisher distribution of directions in 3d with parameter = cappa.
c    b(k) = the ratio of modified bessel functions of the first kind
c    of orders 3/2 and 1/2, equivalent to coth(k) - 1/k.
c
c    for s decimal digit precision and to avoid exponential overflow
c    set son4 = s/4 and bigx = 1.1513*s+0.4, e.g., for 48b = 14.45s;
c
c  Modified:
c
c    30 December 2006
c
c  Author:
c
c    Geoffrey Hill
c
c  Reference:
c
c    Geoffrey Hill,
c    Algorithm 571,
c    Statistics for von Mises's and Fisher's Distributions of Directions,
c    I1(x)/I0(x), I1.5(x)/I0.5(x), and Their Inverses,
c    ACM Transactions on Mathematical Software,
c    Volume 7, Number 2, June 1981, pages 233-238.
c
c  Parameters:
c
c    Input, real CAPPA, the concentration parameter.
c
c    Output, real R4_SPHERR, ?
c
      implicit none

      real bigx
      real cappa
      integer n
      real r4_spherr
      real son4
      real t
      real x
      real xx

      save bigx
      save son4

      data bigx / 17.0E+00 /
      data son4 / 3.61E+00 /

      t = 0.0E+00
      x = abs ( cappa )
c
c  For small X, evaluate Gauss continued fraction x/(3+x*x/(5+...))
c  up from j-th level, where n=2*j+1 yields s decimals precision.
c
      if ( x .le. 1.0E+00 ) then

        n = int ( ( x + 0.88E+00 ) * son4 ) * 2 + 5
        xx = x * x

   10   continue

        t = xx / ( real ( n ) + t )
        n = n - 2

        if ( 3 .lt. n ) then
          go to 10
        end if

        r4_spherr = x / ( t + 3.0E+00 )
c
c  For large X use exponential form of coth(x) - 1/x
c
      else

        if ( x .lt. bigx ) then
          t = 2.0E+00 / ( exp ( 2.0E+00 * x ) - 1.0E+00 )
        end if

        r4_spherr = t + 1.0E+00 - 1.0E+00 / x

      end if

      return
      end
      function r8_cappa3 ( r, d )

c*********************************************************************72
c
cc R8_CAPPA3 estimates KAPPA from the sample mean, for the Fisher distribution. 
c
c  Modified:
c
c    30 December 2006
c
c  Author:
c
c    Geoffrey Hill
c
c  Reference:
c
c    Geoffrey Hill,
c    Algorithm 571,
c    Statistics for von Mises's and Fisher's Distributions of Directions,
c    I1(x)/I0(x), I1.5(x)/I0.5(x), and Their Inverses,
c    ACM Transactions on Mathematical Software,
c    Volume 7, Number 2, June 1981, pages 233-238.
c
c  Parameters:
c
c    Input, double precision R, the modulus of the sample mean vector.
c
c    Input, double precision D, ?
c
c    Output, double precision R8_CAPPA3, the maximum likelihood estimate
c    for the concentration parameter.
c
      implicit none

      double precision d
      integer j
      integer l
      double precision r
      double precision r8_cappa3
      double precision r8_spherr
      double precision s
      double precision t
      double precision x
      double precision y

      y = r

      if ( y .lt. 0.5D+00 ) then

        x = 3.0D+00 * y
        s = x * x
        t = 12.375D+00

        if ( 0.7D+00 .lt. x ) then
          t = ( ( 5.0D+00
     &      * x - 14.74D+00 )
     &      * x + 16.5198D+00 )
     &      * x + 6.2762D+00
        end if

        x = x / ( 1.0D+00 - s
     &        / ( 15.0D+00 - 4.0D+00 * s
     &        / ( 7.0D+00 - s
     &        / ( 5.0D+00 - s / t ))))

      else

        x = 1.0D+00 / ( 1.0D+00 - y )

        if ( 25.0D+00 .lt. x ) then
          r8_cappa3 = x
          return
        end if

        s = 2.0D+00 * x
        t = dexp ( s ) - s * s - 1.0D+00

        if ( y .lt. 0.8D+00 ) then
          t = ((( 0.00254D+00
     &      * s - 0.071042D+00 )
     &      * s + 0.6943388D+00 )
     &      * s - 2.3816184D+00 )
     &      * s + 0.1508478D+00 - 0.14789D+00 * y + t
        end if

        x = t * x / ( t + s )

        if ( 100.0D+00 .lt. x ) then
          r8_cappa3 = x
          return
        end if

      end if

      j = int ( d / 8.0D+00 )

      if ( j .eq. 0 .or. x .lt. 0.01D+00 ) then
        r8_cappa3 = x
        return
      end if

      s = dexp ( x )
      s = s * 2.0D+00 / ( s * s - 1.0D+00 )
      do l = 1, j
        x = x + ( y - r8_spherr ( x, 0.0D+00 ) )
     &    / ( 1.0D+00 / ( x * x ) - s * s )
      end do

      r8_cappa3 = x

      return
      end
      function r8_precis ( a, b )

c*********************************************************************72
c
cc R8_PRECIS returns the number of matching decimal digits.
c
c  Discussion:
c
c   default value 30.0 if perfect match.
c
c  Modified:
c
c    30 December 2006
c
c  Author:
c
c    Geoffrey Hill
c
c  Reference:
c
c    Geoffrey Hill,
c    Algorithm 571,
c    Statistics for von Mises's and Fisher's Distributions of Directions,
c    I1(x)/I0(x), I1.5(x)/I0.5(x), and Their Inverses,
c    ACM Transactions on Mathematical Software,
c    Volume 7, Number 2, June 1981, pages 233-238.
c
c  Parameters:
c
c    Input, double precision A, B, the numbers to be compared.
c
c    Output, double precision R8_PRECIS, the number of matching decimal digits.
c
      implicit none

      double precision a
      double precision b
      double precision c
      double precision r8_precis
      double precision x

      save c

      data c / -0.4342944819D+00 /

      if ( b .eq. 0.0D+00 ) then
        r8_precis = 30.0D+00
        return
      end if

      x = 1.0D+00 - a / b

      if ( x .eq. 0.0D+00 ) then
        r8_precis = 30.0D+00
        return
      end if

      r8_precis = dlog ( dabs ( x ) ) * c

      return
      end
      function r8_spherr ( dcappa, d )

c*********************************************************************72
c
cc R8_SPHERR computes the expected sample mean for the Fisher distribution.
c
c  Modified:
c
c    30 December 2006
c
c  Author:
c
c    Geoffrey Hill
c
c  Reference:
c
c    Geoffrey Hill,
c    Algorithm 571,
c    Statistics for von Mises's and Fisher's Distributions of Directions,
c    I1(x)/I0(x), I1.5(x)/I0.5(x), and Their Inverses,
c    ACM Transactions on Mathematical Software,
c    Volume 7, Number 2, June 1981, pages 233-238.
c
c  Parameters:
c
c    Input, double precision DCAPPA, the concentration parameter.
c
c    Input, double precision D, ?
c
c    Output, double precision R8_SPHERR, ?
c
      implicit none

      double precision d
      double precision dcappa
      integer n
      double precision r8_spherr
      double precision s
      double precision x
      double precision xx

      x = dcappa
      s = 0.0D+00
      xx = x * x

      if ( d .le. 0.0D+00 ) go to 30

      n = int ( d ) * 2 + 1

   10 continue

      if ( n .le. 3 ) then
        r8_spherr = x / ( s + 3.0D+00 )
        return
      end if

      s = xx / ( dble ( n ) + s )
      n = n - 2
      go to 10

   30 continue

      n = 31

      if ( x .le. 1.0D+00 ) then
        go to 10
      end if

      if ( x .le. 34.0D+00 ) then
        r8_spherr = 2.0D+00 / ( dexp ( 2.0D+00 * x ) - 1.0D+00 )
     &    - 1.0D+00 / x + 1.0D+00
      else
        r8_spherr = 1.0D+00 - 1.0D+00 / x
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
c  Modified:
c
c    28 November 2006
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

      character*(8) date
      character*(10) time

      call date_and_time ( date, time )

      write ( *, '(a8,2x,a10)' ) date, time

      return
      end
      function vkappa ( r )

c*********************************************************************72
c
cc VKAPPA estimates KAPPA from the sample mean, for the von Mises distribution.
c
c  Discussion:
c
c    This routine returns vkappa = the maximum likelihood estimate of 
c    'kappa', the concentration parameter of von mises' distribution of 
c    directions in 2 dimensions, corresponding to a sample mean vector 
c    modulus r.
c
c    vkappa = k(a), the inverse function of  a(k) = ratio of modified
c    bessel functions of the first kind, viz., a(k) = i1(k)/i0(k).
c
c  Modified:
c
c    30 December 2006
c
c  Author:
c
c    Geoffrey Hill
c
c  Reference:
c
c    Geoffrey Hill,
c    Algorithm 571,
c    Statistics for von Mises's and Fisher's Distributions of Directions,
c    I1(x)/I0(x), I1.5(x)/I0.5(x), and Their Inverses,
c    ACM Transactions on Mathematical Software,
c    Volume 7, Number 2, June 1981, pages 233-238.
c
c  Parameters:
c
c    Input, real R, the sample mean vector modulus.
c
c    Output, real VKAPPA, the maximum likelihood estimate for the
c    concentration parameter.
c
      implicit none

      real a
      real besrat
      real r
      real s
      real v1
      real v2
      real vkappa
      real x
      real y

      save v1
      save v2
c
c  For 8s (significant decimal digits) precision auxiliary routine
c  function besrat(v) must be set to at least 9.3s
c
      data v1 / 0.642E+00 /
      data v2 / 0.95E+00 /

      a = r
      s = -1.0E+00
c
c  Error signal: value -1.0 returned if argument negative or 1.0 or more.
c
      if ( a .lt. 0.0E+00 .or. 1.0E+00 .le. a ) then
        vkappa = s
        return
      end if

      y = 2.0E+00 / ( 1.0E+00 - a )
c
c   For R above 0.85 use continued fraction approximation.
c
      if ( 0.85E+00 .lt. a ) then

        if ( 0.95E+00 .le. a ) then
          x = 32.0E+00 / ( 120.0E+00 * a - 131.5E+00 + y )
        else
          x = ( -2326.0E+00
     &      * a + 4317.5526E+00 )
     &      * a - 2001.035224E+00
        end if

        s = ( y + 1.0E+00 + 3.0E+00
     &    / ( y - 5.0E+00 - 12.0E+00
     &    / ( y - 10.0E+00 - x ) ) ) * 0.25E+00

        if ( v2 .le. a ) then
          vkappa = s
          return
        end if
c
c   For R below 0.85, use adjusted inverse Taylor series.
c
      else

        x = a * a
        s = ((( a - 5.6076E+00 )
     &    * a + 5.0797E+00 )
     &    * a - 4.6494E+00 ) * y * x - 1.0E+00

        s = (((( s
     &    * x + 15.0E+00 )
     &    * x + 60.0E+00 )
     &    * x / 360.0E+00 + 1.0E+00 )
     &    * x - 2.0E+00 ) * a / ( x - 1.0E+00 )

        if ( a .lt. v1 ) then
          vkappa = s
          return
        end if

      end if
c
c   For R in (0.642,0.95) apply Newton-Raphson, twice if R in
c   (0.75,0.875), for 8s precision, using approximate derivative -
c
      y = ( ( 0.00048E+00
     &  * y - 0.1589E+00 )
     &  * y + 0.744E+00 )
     &  * y - 4.2932E+00

      if ( a .le. 0.875E+00 ) then
        s = ( besrat ( s ) - a ) * y + s
      end if

      if ( 0.75E+00 .le. a ) then
        s = ( besrat ( s ) - a ) * y + s
      end if

      vkappa = s

      return
      end
