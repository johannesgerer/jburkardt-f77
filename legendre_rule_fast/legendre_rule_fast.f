      program main

c*********************************************************************72
c
cc MAIN is the main program for LEGENDRE_RULE_FAST.
c
c  Discussion:
c
c    This program computes a standard Gauss-Legendre quadrature rule
c    and writes it to a file.
c
c  Usage:
c
c    legendre_rule_fast n a b
c
c    where
c
c    * n is the number of points in the rule;
c    * a is the left endpoint;
c    * b is the right endpoint.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    15 September 2010
c
c  Author:
c
c    John Burkardt
c
      implicit none

      double precision  a
      integer arg_num
      double precision  b
      integer iarg
      integer iargc
      integer ierror
      integer last
      integer n
      character * 255 string

      call timestamp ( )
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'LEGENDRE_RULE_FAST'
      write ( *, '(a)' ) '  FORTRAN77 version'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 
     &  '  Compute a Gauss-Legendre rule for approximating'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '    Integral ( b <= x <= b ) f(x) dx'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  of order N.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  The computed rule is written to 3 files:'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '    leg_oN_w.txt - the weight file'
      write ( *, '(a)' ) '    leg_oN_x.txt - the abscissa file.'
      write ( *, '(a)' ) '    leg_oN_r.txt - the region file.'
c
c  Get the number of command line arguments.
c
      arg_num = iargc ( )
c
c  Get N.
c
      if ( 1 .le. arg_num ) then
      
        iarg = 1
        call getarg ( iarg, string )
        call s_to_i4 ( string, n, ierror, last )
        
      else
      
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) '  Enter the rule order N:'
        read ( *, * ) n
        
      end if

      write ( *, '(a)' ) ' '
      write ( *, '(a,i8)' ) '  The requested order N ', n
c
c  Get A.
c
      if ( 2 .le. arg_num ) then
      
        iarg = 2
        call getarg ( iarg, string )
        call s_to_r8 ( string, a, ierror, last )
        
      else
      
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) '  Enter the left endpoint A:'
        read ( *, * ) a
        
      end if

      write ( *, '(a)' ) ' '
      write ( *, '(a,g14.6)' ) '  The left endpoint A is ', a
c
c  Get B.
c
      if ( 3 .le. arg_num ) then
      
        iarg = 3
        call getarg ( iarg, string )
        call s_to_r8 ( string, b, ierror, last )
        
      else
      
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) '  Enter the right endpoint B:'
        read ( *, * ) b
        
      end if

      write ( *, '(a)' ) ' '
      write ( *, '(a,g14.6)' ) '  The right endpoint B is ', b
c
c  Construct the rule and output it.
c
      call legendre_handle ( n, a, b )
c
c  Terminate.
c
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'LEGENDRE_RULE_FAST:'
      write ( *, '(a)' ) '  Normal end of execution.'
      write ( *, '(a)' ) ' '
      call timestamp ( )

      stop
      end
      subroutine ch_cap ( ch )

c*********************************************************************72
c
cc CH_CAP capitalizes a single character.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    03 January 2007
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input/output, character CH, the character to capitalize.
c
      implicit none

      character ch
      integer itemp

      itemp = ichar ( ch )

      if ( 97 .le. itemp .and. itemp .le. 122 ) then
        ch = char ( itemp - 32 )
      end if

      return
      end
      function ch_eqi ( c1, c2 )

c*********************************************************************72
c
cc CH_EQI is a case insensitive comparison of two characters for equality.
c
c  Example:
c
c    CH_EQI ( 'A', 'a' ) is TRUE.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    03 January 2007
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, character C1, C2, the characters to compare.
c
c    Output, logical CH_EQI, the result of the comparison.
c
      implicit none

      character c1
      character c1_cap
      character c2
      character c2_cap
      logical ch_eqi

      c1_cap = c1
      c2_cap = c2

      call ch_cap ( c1_cap )
      call ch_cap ( c2_cap )

      if ( c1_cap .eq. c2_cap ) then
        ch_eqi = .true.
      else
        ch_eqi = .false.
      end if

      return
      end
      subroutine ch_to_digit ( c, digit )

c*********************************************************************72
c
cc CH_TO_DIGIT returns the integer value of a base 10 digit.
c
c  Example:
c
c     C   DIGIT
c    ---  -----
c    '0'    0
c    '1'    1
c    ...  ...
c    '9'    9
c    ' '    0
c    'X'   -1
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    04 August 1999
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, character C, the decimal digit, '0' through '9' or blank
c    are legal.
c
c    Output, integer DIGIT, the corresponding integer value.  If C was
c    'illegal', then DIGIT is -1.
c
      implicit none

      character c
      integer digit

      if ( lge ( c, '0' ) .and. lle ( c, '9' ) ) then

        digit = ichar ( c ) - 48

      else if ( c .eq. ' ' ) then

        digit = 0

      else

        digit = -1

      end if

      return
      end
      subroutine digit_to_ch ( digit, c )

c*********************************************************************72
c
cc DIGIT_TO_CH returns the character representation of a decimal digit.
c
c  Example:
c
c    DIGIT   C
c    -----  ---
c      0    '0'
c      1    '1'
c    ...    ...
c      9    '9'
c     17    '*'
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    04 August 1999
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer DIGIT, the digit value between 0 and 9.
c
c    Output, character C, the corresponding character, or '*' if DIGIT
c    was illegal.
c
      implicit none

      character c
      integer digit

      if ( 0 .le. digit .and. digit .le. 9 ) then

        c = char ( digit + 48 )

      else

        c = '*'

      end if

      return
      end
      subroutine get_unit ( unit )

c*********************************************************************72
c
cc GET_UNIT returns a free FORTRAN unit number.
c
c  Discussion:
c
c    A "free" FORTRAN unit number is a value between 1 and 99 which
c    is not currently associated with an I/O device.  A free FORTRAN unit
c    number is needed in order to open a file with the OPEN command.
c
c    If UNIT = 0, then no free FORTRAN unit could be found, although
c    all 99 units were checked (except for units 5, 6 and 9, which
c    are commonly reserved for console I/O).
c
c    Otherwise, UNIT is a value between 1 and 99, representing a
c    free FORTRAN unit.  Note that GET_UNIT assumes that units 5 and 6
c    are special, and will never return those values.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    26 October 2008
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, integer UNIT, the free unit number.
c
      implicit none

      integer i
      integer unit

      unit = 0

      do i = 1, 99

        if ( i .ne. 5 .and. i .ne. 6 .and. i .ne. 9 ) then

          open ( unit = i, err = 10, status = 'scratch' )
          close ( unit = i )

          unit = i

          return
        end if

10      continue

      end do

      return
      end
      subroutine i4_to_s_left ( i4, s )

c*********************************************************************72
c
cc I4_TO_S_LEFT converts an I4 to a left-justified string.
c
c  Discussion:
c
c    An I4 is an integer.
c
c  Example:
c
c    Assume that S is 6 characters long:
c
c        I4  S
c
c         1  1
c        -1  -1
c         0  0
c      1952  1952
c    123456  123456
c   1234567  ******  <-- Not enough room!
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    15 September 2010
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer I4, an integer to be converted.
c
c    Output, character * ( * ) S, the representation of the integer.
c    The integer will be left-justified.  If there is not enough space,
c    the string will be filled with stars.
c
      implicit none

      character c
      integer i
      integer i4
      integer idig
      integer ihi
      integer ilo
      integer ipos
      integer ival
      character * ( * ) s

      s = ' '

      ilo = 1
      ihi = len ( s )

      if ( ihi .le. 0 ) then
        return
      end if
c
c  Make a copy of the integer.
c
      ival = i4
c
c  Handle the negative sign.
c
      if ( ival .lt. 0 ) then

        if ( ihi .le. 1 ) then
          s(1:1) = '*'
          return
        end if

        ival = -ival
        s(1:1) = '-'
        ilo = 2

      end if
c
c  The absolute value of the integer goes into S(ILO:IHI).
c
      ipos = ihi
c
c  Find the last digit of IVAL, strip it off, and stick it into the string.
c
10    continue

        idig = mod ( ival, 10 )
        ival = ival / 10

        if ( ipos .lt. ilo ) then
          do i = 1, ihi
            s(i:i) = '*'
          end do
          return
        end if

        call digit_to_ch ( idig, c )

        s(ipos:ipos) = c
        ipos = ipos - 1

        if ( ival .eq. 0 ) then
          go to 20
        end if

      go to 10

20    continue
c
c  Shift the string to the left.
c
      s(ilo:ilo+ihi-ipos-1) = s(ipos+1:ihi)
      s(ilo+ihi-ipos:ihi) = ' '
 
      return
      end
      subroutine legendre_compute_glr ( n, x, w )

c*********************************************************************72
c
cc LEGENDRE_COMPUTE_GLR: Legendre quadrature by the Glaser-Liu-Rokhlin method.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    20 October 2009
c
c  Author:
c
c    Original MATLAB version by Nick Hale.
c    FORTRAN77 version by John Burkardt.
c
c  Reference:
c
c    Andreas Glaser, Xiangtao Liu, Vladimir Rokhlin, 
c    A fast algorithm for the calculation of the roots of special functions, 
c    SIAM Journal on Scientific Computing,
c    Volume 29, Number 4, pages 1420-1438, 2007.
c
c  Parameters:
c
c    Input, integer N, the order.
c
c    Output, double precision X(N), the abscissas.
c
c    Output, double precision W(N), the weights.
c
      implicit none

      integer n

      integer i
      double precision p
      double precision pp
      double precision r8vec_sum
      double precision w(n)
      double precision w_sum
      double precision x(n)
c
c  Get the value and derivative of the N-th Legendre polynomial at 0.0.
c
      call legendre_compute_glr0 ( n, p, pp )
c
c  If N is odd, then zero is a root.
c
      if ( mod ( n, 2 ) .eq. 1 ) then

        x((n+1)/2) = 0.0D+00
        w((n+1)/2) = pp
c
c  If N is even, we have to compute a root.
c
      else

        call legendre_compute_glr2 ( p, n, x((n/2)+1), w((n/2)+1) )

      end if
c
c  Get the complete set of roots and derivatives.
c
      call legendre_compute_glr1 ( n, x, w )
c
c  Compute W.
c
      do i = 1, n
        w(i) = 2.0D+00 / 
     &    ( 1.0D+00 - x(i) ) / ( 1.0D+00 + x(i) ) / w(i) / w(i)
      end do

      w_sum = r8vec_sum ( n, w )

      do i = 1, n
        w(i) = 2.0D+00 * w(i) / w_sum
      end do

      return
      end
      subroutine legendre_compute_glr0 ( n, p, pp )

c*********************************************************************72
c
cc LEGENDRE_COMPUTE_GLR0 gets a starting value for the fast algorithm.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    16 October 2009
c
c  Author:
c
c    Original MATLAB version by Nick Hale.
c    FORTRAN77 version by John Burkardt.
c
c  Reference:
c
c    Andreas Glaser, Xiangtao Liu, Vladimir Rokhlin, 
c    A fast algorithm for the calculation of the roots of special functions, 
c    SIAM Journal on Scientific Computing,
c    Volume 29, Number 4, pages 1420-1438, 2007.
c
c  Parameters:
c
c    Input, integer N, the order of the Legendre polynomial.
c
c    Output, double precision P, PP, the value of the N-th Legendre polynomial
c    and its derivative at 0.
c
      implicit none

      integer n

      integer k
      double precision p
      double precision pm1
      double precision pm2
      double precision pp
      double precision ppm1
      double precision ppm2
      double precision rk
c
c  Compute coefficients of P_m(0), Pm'(0), m = 0,..,N
c
      pm2 = 0.0D+00
      pm1 = 1.0D+00
      ppm2 = 0.0D+00
      ppm1 = 0.0D+00

      do k = 0, n - 1
        rk = dble ( k )
        p = - rk * pm2 / ( rk + 1.0D+00 )
        pp = ( ( 2.0D+00 * rk + 1.0D+00 ) * pm1 
     &                   - rk             * ppm2 ) 
     &       / (           rk + 1.0D+00 )
        pm2 = pm1
        pm1 = p
        ppm2 = ppm1
        ppm1 = pp
      end do

      return
      end
      subroutine legendre_compute_glr1 ( n, x, w )

c*********************************************************************72
c
cc LEGENDRE_COMPUTE_GLR1 gets the complete set of Legendre points and weights.
c
c  Discussion:
c
c    This routine requires that a starting estimate be provided for one
c    root and its derivative.  This information will be stored in entry
c    (N+1)/2 if N is odd, or N/2 if N is even, of X and W.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    15 November 2009
c
c  Author:
c
c    Original C++ version by Nick Hale.
c    FORTRAN77 version by John Burkardt.
c
c  Reference:
c
c    Andreas Glaser, Xiangtao Liu, Vladimir Rokhlin, 
c    A fast algorithm for the calculation of the roots of special functions, 
c    SIAM Journal on Scientific Computing,
c    Volume 29, Number 4, pages 1420-1438, 2007.
c
c  Parameters:
c
c    Input, integer N, the order of the Legendre polynomial.
c
c    Input/output, double precision X(N).  On input, a starting value
c    has been set in one entry.  On output, the roots of the Legendre 
c    polynomial.
c
c    Input/output, double precision W(N).  On input, a starting value
c    has been set in one entry.  On output, the derivatives of the Legendre 
c    polynomial at the zeros.
c
c  Local Parameters:
c
c    Local, integer M, the number of terms in the Taylor expansion.
c
      implicit none

      integer m
      parameter ( m = 30 )
      integer n

      double precision dk
      double precision dn
      double precision h
      integer j
      integer k
      integer l
      integer n2
      double precision pi
      parameter ( pi = 3.141592653589793D+00 )
      double precision rk2_leg
      integer s
      double precision ts_mult
      double precision u(m+2)
      double precision up(m+1)
      double precision w(n)
      double precision x(n)
      double precision xp

      if ( mod ( n, 2 ) .eq. 1 ) then
        n2 = ( n - 1 ) / 2 - 1
        s = 1
      else
        n2 = n / 2 - 1
        s = 0
      end if

      dn = dble ( n )

      do j = n2 + 1, n - 2

        xp = x(j+1)

        h = rk2_leg ( pi/2.0D+00, -pi/2.0D+00, xp, n ) - xp

        u(1) = 0.0D+00
        u(2) = 0.0D+00
        u(3) = w(j+1)

        up(1) = 0.0D+00
        up(2) = u(3)

        do k = 0, m - 2

          dk = dble ( k )

          u(k+4) = 
     &    ( 
     &      2.0D+00 * xp * ( dk + 1.0D+00 ) * u(k+3) 
     &      + ( dk * ( dk + 1.0D+00 ) - dn * ( dn + 1.0D+00 ) ) 
     &      * u(k+2) / ( dk + 1.0D+00 ) 
     &    ) / ( 1.0D+00 - xp ) / ( 1.0D+00 + xp ) / ( dk + 2.0D+00 )

          up(k+3) = ( dk + 2.0D+00 ) * u(k+4)

        end do

        do l = 0, 4
          h = h - ts_mult ( u, h, m ) / ts_mult ( up, h, m-1 )
        end do

        x(j+2) = xp + h
        w(j+2) = ts_mult ( up, h, m - 1 )   

      end do

      do k = 0, n2 + s
        x(k+1) = - x(n-1-k+1)
        w(k+1) = w(n-1-k+1)
      end do

      return
      end
      subroutine legendre_compute_glr2 ( pn0, n, x1, d1 )

c*********************************************************************72
c
cc LEGENDRE_COMPUTE_GLR2 finds the first real root.
c
c  Discussion:
c
c    This routine is only called if N is even.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    20 October 2009
c
c  Author:
c
c    Original MATLAB version by Nick Hale.
c    FORTRAN77 version by John Burkardt.
c
c  Reference:
c
c    Andreas Glaser, Xiangtao Liu, Vladimir Rokhlin, 
c    A fast algorithm for the calculation of the roots of special functions, 
c    SIAM Journal on Scientific Computing,
c    Volume 29, Number 4, pages 1420-1438, 2007.
c
c  Parameters:
c
c    Input, double precision PN0, the value of the N-th Legendre polynomial
c    at 0.
c
c    Input, integer N, the order of the Legendre polynomial.
c
c    Output, double precision X1, the first real root.
c
c    Output, double precision D1, the derivative at X1.
c
c  Local Parameters:
c
c    Local, integer M, the number of terms in the Taylor expansion.
c
      implicit none

      integer m
      parameter ( m = 30 )

      double precision d1
      double precision eps
      integer i
      integer k
      integer kk
      integer l
      integer n
      double precision pi
      parameter ( pi = 3.141592653589793D+00 )
      double precision pn0
      double precision r8_epsilon
      double precision r8_huge
      double precision r8vec_dot_product
      double precision rk
      double precision rn
      double precision scale
      double precision step
      double precision theta
      double precision u(m+1)
      double precision up(m+1)
      double precision x1
      double precision x1k(m+1)

      k = ( n + 1 ) / 2

      theta = pi * dble ( 4 * k - 1 ) / dble ( 4 * n + 2 )

      x1 = ( 1.0D+00 - dble ( n - 1 ) 
     &  / dble ( 8 * n * n * n ) 
     &  - 1.0D+00 / dble ( 384 * n * n * n * n ) 
     &  * ( 39.0D+00 - 28.0D+00 / ( sin ( theta ) * sin ( theta ) ) ) ) 
     &  * cos ( theta )
c
c  Scaling.
c
      scale = 1.0D+00 / x1
c
c  Recurrence relation for Legendre polynomials.
c
      u(1:m+1) = 0.0D+00
      up(1:m+1) = 0.0D+00

      rn = dble ( n )

      u(1) = pn0

      do k = 0, m - 2, 2

        rk = dble ( k )

        u(k+3) = ( rk * ( rk + 1.0D+00 ) 
     &    - rn * ( rn + 1.0D+00 ) ) * u(k+1) 
     &    / ( rk + 1.0D+00 ) / ( rk + 2.0D+00 ) / scale / scale

        up(k+2) = ( rk + 2.0D+00 ) * u(k+3) * scale

      end do
c
c  Flip for more accuracy in inner product calculation
c
      u = u(m+1:1:-1)
      up = up(m+1:1:-1)

      x1k(1:m+1) = 1.0D+00

      step = r8_huge ( )
      l = 0
c
c  Newton iteration.
c
      eps = r8_epsilon ( )

      do while ( eps < abs ( step ) .and. l .lt. 10 )
        l = l + 1
        step = r8vec_dot_product ( m + 1, u,  x1k ) 
     &       / r8vec_dot_product ( m + 1, up, x1k )
        x1 = x1 - step
        x1k(1) = 1.0D+00
        x1k(2) = scale * x1
        do kk = 3, m + 1
          x1k(kk) = x1k(kk-1) * scale * x1
        end do
        x1k(1:m+1) = x1k(m+1:1:-1)
      end do

      d1 = r8vec_dot_product ( m + 1, up, x1k )

      return
      end
      subroutine legendre_handle ( n, a, b )

c*********************************************************************72
c
cc LEGENDRE_HANDLE computes the requested Gauss-Legendre rule and outputs it.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    24 October 2009
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the order of the rule.
c
c    Input, double precision A, B, the left and right endpoints.
c 
      implicit none

      double precision a
      double precision b
      integer i
      integer n
      character ( len = 80 ) output_r
      character ( len = 80 ) output_w
      character ( len = 80 ) output_x
      double precision r(2)
      double precision t1
      double precision t2
      character * 10 tag
      double precision, allocatable, dimension ( : ) :: w
      double precision, allocatable, dimension ( : ) :: x

      r(1) = a
      r(2) = b
c
c  Compute the rule.
c
      allocate ( w(n) )
      allocate ( x(n) )

      call cpu_time ( t1 )
      call legendre_compute_glr ( n, x, w )
      call cpu_time ( t2 )
      
      write ( *, '(a)' ) ' '
      write ( *, '(a,g14.6,a)' ) 
     &  '  Computation required ', t2 - t1, ' seconds.'
c
c  Rescale the data.
c
      call rescale ( n, a, b, x, w )
c
c  Write the data to files.
c
      call i4_to_s_left ( n, tag )

      output_w = 'leg_o' // trim ( tag ) // '_w.txt'
      output_x = 'leg_o' // trim ( tag ) // '_x.txt'
      output_r = 'leg_o' // trim ( tag ) // '_r.txt'

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 
     &  '  Weight file will be   "' // trim ( output_w ) // '".'
      write ( *, '(a)' ) 
     &  '  Abscissa file will be "' // trim ( output_x ) // '".'
      write ( *, '(a)' ) 
     &  '  Region file will be   "' // trim ( output_r ) // '".'
                
      call r8mat_write ( output_w, 1, n, w )
      call r8mat_write ( output_x, 1, n, x )
      call r8mat_write ( output_r, 1, 2, r )
          
      deallocate ( w )
      deallocate ( x )

      return
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

      double precision one
      double precision r8_epsilon
      double precision temp
      double precision test
      double precision value

      one = dble ( 1 )

      value = one
      temp = value / 2.0D+00
      test = one + temp

10    continue

      if ( one .lt. test ) then
        value = temp
        temp = value / 2.0D+00
        test = one + temp
        go to 10
      end if

      r8_epsilon = value

      return
      end
      function r8_huge ( )

c*********************************************************************72
c
cc R8_HUGE returns a "huge" R8.
c
c  Discussion:
c
c    The value returned by this function is NOT required to be the
c    maximum representable R8.  This value varies from machine to machine,
c    from compiler to compiler, and may cause problems when being printed.
c    We simply want a "very large" but non-infinite number.
c
c    FORTRAN90 provides a built-in routine HUGE ( X ) that
c    can return the maximum representable number of the same datatype
c    as X, if that is what is really desired.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    13 April 2004
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, double precision R8_HUGE, a huge number.
c
      implicit none

      double precision r8_huge

      r8_huge = 1.0D+30

      return
      end
      subroutine r8mat_write ( output_filename, m, n, table )

c*********************************************************************72
c
cc R8MAT_WRITE writes a R8MAT file.
c
c  Discussion:
c
c    An R8MAT is an array of R8's.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    22 October 2009
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, character * ( * ) OUTPUT_FILENAME, the output file name.
c
c    Input, integer M, the spatial dimension.
c
c    Input, integer N, the number of points.
c
c    Input, double precision TABLE(M,N), the table data.
c
      implicit none

      integer m
      integer n

      integer j
      character * ( * ) output_filename
      integer output_unit
      character * ( 30 ) string
      double precision table(m,n)
c
c  Open the file.
c
      call get_unit ( output_unit )

      open ( unit = output_unit, file = output_filename,
     &  status = 'replace' )
c
c  Create the format string.
c
      if ( 0 .lt. m .and. 0 .lt. n ) then

        write ( string, '(a1,i8,a1,i8,a1,i8,a1)' )
     &    '(', m, 'g', 24, '.', 16, ')'
c
c  Write the data.
c
        do j = 1, n
          write ( output_unit, string ) table(1:m,j)
        end do

      end if
c
c  Close the file.
c
      close ( unit = output_unit )

      return
      end
      function r8vec_dot_product ( n, v1, v2 )

c*********************************************************************72
c
cc R8VEC_DOT_PRODUCT finds the dot product of a pair of R8VEC's.
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
c    Output, double precision R8VEC_DOT_PRODUCT, the dot product.
c
      implicit none

      integer n

      integer i
      double precision r8vec_dot_product
      double precision v1(n)
      double precision v2(n)
      double precision value

      value = 0.0D+00
      do i = 1, n
        value = value + v1(i) * v2(i)
      end do

      r8vec_dot_product = value

      return
      end
      function r8vec_sum ( n, v1 )

c*********************************************************************72
c
cc R8VEC_SUM sums the entries of an R8VEC.
c
c  Discussion:
c
c    An R8VEC is a vector of R8's.
c
c    In FORTRAN90, the system routine SUM should be called
c    directly.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    22 July 2008
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the dimension of the vectors.
c
c    Input, double precision V1(N), the vector.
c
c    Output, double precision R8VEC_SUM, the sum of the entries.
c
      implicit none

      integer n

      integer i
      double precision r8vec_sum
      double precision v1(n)
      double precision value

      value = 0.0D+00
      do i = 1, n
        value = value + v1(i)
      end do

      r8vec_sum = value

      return
      end
      subroutine rescale ( n, a, b, x, w )

c*********************************************************************72
c
cc RESCALE rescales a Legendre quadrature rule from [-1,+1] to [A,B].
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    18 October 2009
c
c  Author:
c
c    Original MATLAB version by Nick Hale.
c    FORTRAN77 version by John Burkardt.
c
c  Reference:
c
c    Andreas Glaser, Xiangtao Liu, Vladimir Rokhlin, 
c    A fast algorithm for the calculation of the roots of special functions, 
c    SIAM Journal on Scientific Computing,
c    Volume 29, Number 4, pages 1420-1438, 2007.
c
c  Parameters:
c
c    Input, integer N, the order.
c
c    Input, double precision A, B, the endpoints of the new interval.
c
c    Input/output, double precision X(N), on input, the abscissas for [-1,+1].
c    On output, the abscissas for [A,B].
c
c    Input/output, double precision W(N), on input, the weights for [-1,+1].
c    On output, the weights for [A,B].
c
      implicit none

      integer n

      double precision a
      double precision b
      integer i
      double precision w(n)
      double precision x(n)

      do i = 1, n
        x(i) = ( ( a + b ) + ( b - a ) * x(i) ) / 2.0D+00
      end do

      do i = 1, n
        w(i) = ( b - a ) * w(i) / 2.0D+00
      end do

      return
      end
      function rk2_leg ( t1, t2, x, n )

c*********************************************************************72
c
cc RK2_LEG advances the value of X(T) using a Runge-Kutta method.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    15 November 2009
c
c  Author:
c
c    Original C++ version by Nick Hale.
c    FORTRAN77 version by John Burkardt.
c
c  Parameters:
c
c    Input, double precision T1, T2, the range of the integration interval.
c
c    Input, double precision X, the value of X at T1.
c
c    Input, integer N, the number of steps to take.
c
c    Output, double precision RK2_LEG, the value of X at T2.
c
      implicit none

      double precision f
      double precision h
      integer j
      double precision k1
      double precision k2
      integer m
      parameter ( m = 10 )
      integer n
      double precision rk2_leg
      double precision snn1
      double precision t
      double precision t1
      double precision t2
      double precision x
      double precision x2

      x2 = x

      h = ( t2 - t1 ) / dble ( m )
      snn1 = sqrt ( dble ( n * ( n + 1 ) ) )
      t = t1

      do j = 0, m - 1

        f = ( 1.0D+00 - x2 ) * ( 1.0D+00 + x2 )
        k1 = - h * f / ( snn1 * sqrt ( f ) 
     &    - 0.5D+00 * x2 * sin ( 2.0D+00 * t ) )
        x2 = x2 + k1

        t = t + h

        f = ( 1.0D+00 - x2 ) * ( 1.0D+00 + x2 )
        k2 = - h * f / ( snn1 * sqrt ( f ) 
     &    - 0.5D+00 * x2 * sin ( 2.0D+00 * t ) )
        x2 = x2 + 0.5D+00 * ( k2 - k1 )

      end do

      rk2_leg = x2

      return
      end
      function s_len_trim ( s )

c*********************************************************************72
c
cc S_LEN_TRIM returns the length of a string to the last nonblank.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    05 March 2004
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, character*(*) S, a string.
c
c    Output, integer S_LEN_TRIM, the length of the string to the last nonblank.
c
      implicit none

      integer i
      character*(*) s
      integer s_len_trim

      do i = len ( s ), 1, -1

        if ( s(i:i) .ne. ' ' ) then
          s_len_trim = i
          return
        end if

      end do

      s_len_trim = 0

      return
      end
      subroutine s_to_i4 ( s, ival, ierror, length )

c*********************************************************************72
c
cc S_TO_I4 reads an I4 from a string.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    28 April 2008
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, character * ( * ) S, a string to be examined.
c
c    Output, integer IVAL, the integer value read from the string.
c    If the string is blank, then IVAL will be returned 0.
c
c    Output, integer IERROR, an error flag.
c    0, no error.
c    1, an error occurred.
c
c    Output, integer LENGTH, the number of characters of S
c    used to make IVAL.
c
      implicit none

      character c
      integer i
      integer ierror
      integer isgn
      integer istate
      integer ival
      integer length
      character * ( * ) s
      integer s_len_trim

      ierror = 0
      istate = 0
      isgn = 1
      ival = 0

      do i = 1, s_len_trim ( s )

        c = s(i:i)
c
c  Haven't read anything.
c
        if ( istate .eq. 0 ) then

          if ( c .eq. ' ' ) then

          else if ( c .eq. '-' ) then
            istate = 1
            isgn = -1
          else if ( c .eq. '+' ) then
            istate = 1
            isgn = + 1
          else if ( lle ( '0', c ) .and. lle ( c, '9' ) ) then
            istate = 2
            ival = ichar ( c ) - ichar ( '0' )
          else
            ierror = 1
            return
          end if
c
c  Have read the sign, expecting digits.
c
        else if ( istate .eq. 1 ) then

          if ( c .eq. ' ' ) then

          else if ( lle ( '0', c ) .and. lle ( c, '9' ) ) then
            istate = 2
            ival = ichar ( c ) - ichar ( '0' )
          else
            ierror = 1
            return
          end if
c
c  Have read at least one digit, expecting more.
c
        else if ( istate .eq. 2 ) then

          if ( lle ( '0', c ) .and. lle ( c, '9' ) ) then
            ival = 10 * ival + ichar ( c ) - ichar ( '0' )
          else
            ival = isgn * ival
            length = i - 1
            return
          end if

        end if

      end do
c
c  If we read all the characters in the string, see if we're OK.
c
      if ( istate .eq. 2 ) then
        ival = isgn * ival
        length = s_len_trim ( s )
      else
        ierror = 1
        length = 0
      end if

      return
      end
      subroutine s_to_r8 ( s, dval, ierror, length )

c*********************************************************************72
c
cc S_TO_R8 reads an R8 from a string.
c
c  Discussion:
c
c    The routine will read as many characters as possible until it reaches
c    the end of the string, or encounters a character which cannot be
c    part of the number.
c
c    Legal input is:
c
c       1 blanks,
c       2 '+' or '-' sign,
c       2.5 blanks
c       3 integer part,
c       4 decimal point,
c       5 fraction part,
c       6 'E' or 'e' or 'D' or 'd', exponent marker,
c       7 exponent sign,
c       8 exponent integer part,
c       9 exponent decimal point,
c      10 exponent fraction part,
c      11 blanks,
c      12 final comma or semicolon,
c
c    with most quantities optional.
c
c  Example:
c
c    S                 DVAL
c
c    '1'               1.0
c    '     1   '       1.0
c    '1A'              1.0
c    '12,34,56'        12.0
c    '  34 7'          34.0
c    '-1E2ABCD'        -100.0
c    '-1X2ABCD'        -1.0
c    ' 2E-1'           0.2
c    '23.45'           23.45
c    '-4.2E+2'         -420.0
c    '17d2'            1700.0
c    '-14e-2'         -0.14
c    'e2'              100.0
c    '-12.73e-9.23'   -12.73 * 10.0**(-9.23)
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    28 April 2008
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, character * ( * ) S, the string containing the
c    data to be read.  Reading will begin at position 1 and
c    terminate at the end of the string, or when no more
c    characters can be read to form a legal real.  Blanks,
c    commas, or other nonnumeric data will, in particular,
c    cause the conversion to halt.
c
c    Output, double precision DVAL, the value read from the string.
c
c    Output, integer IERROR, error flag.
c    0, no errors occurred.
c    1, 2, 6 or 7, the input number was garbled.  The
c    value of IERROR is the last type of input successfully
c    read.  For instance, 1 means initial blanks, 2 means
c    a plus or minus sign, and so on.
c
c    Output, integer LENGTH, the number of characters read
c    to form the number, including any terminating
c    characters such as a trailing comma or blanks.
c
      implicit none

      logical ch_eqi
      character c
      double precision dval
      integer ierror
      integer ihave
      integer isgn
      integer iterm
      integer jbot
      integer jsgn
      integer jtop
      integer length
      integer nchar
      integer ndig
      double precision rbot
      double precision rexp
      double precision rtop
      character * ( * ) s
      integer s_len_trim

      nchar = s_len_trim ( s )

      ierror = 0
      dval = 0.0D+00
      length = -1
      isgn = 1
      rtop = 0
      rbot = 1
      jsgn = 1
      jtop = 0
      jbot = 1
      ihave = 1
      iterm = 0

10    continue

        length = length + 1

        if ( nchar .lt. length+1 ) then
          go to 20
        end if

        c = s(length+1:length+1)
c
c  Blank character.
c
        if ( c .eq. ' ' ) then

          if ( ihave .eq. 2 ) then

          else if ( ihave .eq. 6 .or. ihave .eq. 7 ) then
            iterm = 1
          else if ( 1 .lt. ihave ) then
            ihave = 11
          end if
c
c  Comma.
c
        else if ( c .eq. ',' .or. c .eq. ';' ) then

          if ( ihave .ne. 1 ) then
            iterm = 1
            ihave = 12
            length = length + 1
          end if
c
c  Minus sign.
c
        else if ( c .eq. '-' ) then

          if ( ihave .eq. 1 ) then
            ihave = 2
            isgn = -1
          else if ( ihave .eq. 6 ) then
            ihave = 7
            jsgn = -1
          else
            iterm = 1
          end if
c
c  Plus sign.
c
        else if ( c .eq. '+' ) then

          if ( ihave .eq. 1 ) then
            ihave = 2
          else if ( ihave .eq. 6 ) then
            ihave = 7
          else
            iterm = 1
          end if
c
c  Decimal point.
c
        else if ( c .eq. '.' ) then

          if ( ihave .lt. 4 ) then
            ihave = 4
          else if ( 6 .le. ihave .and. ihave .le. 8 ) then
            ihave = 9
          else
            iterm = 1
          end if
c
c  Scientific notation exponent marker.
c
        else if ( ch_eqi ( c, 'E' ) .or. ch_eqi ( c, 'D' ) ) then

          if ( ihave .lt. 6 ) then
            ihave = 6
          else
            iterm = 1
          end if
c
c  Digit.
c
        else if ( ihave .lt. 11 .and. lle ( '0', c ) 
     &    .and. lle ( c, '9' ) ) then

          if ( ihave .le. 2 ) then
            ihave = 3
          else if ( ihave .eq. 4 ) then
            ihave = 5
          else if ( ihave .eq. 6 .or. ihave .eq. 7 ) then
            ihave = 8
          else if ( ihave .eq. 9 ) then
            ihave = 10
          end if

          call ch_to_digit ( c, ndig )

          if ( ihave .eq. 3 ) then
            rtop = 10.0D+00 * rtop + dble ( ndig )
          else if ( ihave .eq. 5 ) then
            rtop = 10.0D+00 * rtop + dble ( ndig )
            rbot = 10.0D+00 * rbot
          else if ( ihave .eq. 8 ) then
            jtop = 10 * jtop + ndig
          else if ( ihave .eq. 10 ) then
            jtop = 10 * jtop + ndig
            jbot = 10 * jbot
          end if
c
c  Anything else is regarded as a terminator.
c
        else
          iterm = 1
        end if
c
c  If we haven't seen a terminator, and we haven't examined the
c  entire string, go get the next character.
c
        if ( iterm .eq. 1 ) then
          go to 20
        end if

        go to 10

20    continue
c
c  If we haven't seen a terminator, and we have examined the
c  entire string, then we're done, and LENGTH is equal to NCHAR.
c
      if ( iterm .ne. 1 .and. length+1 .eq. nchar ) then
        length = nchar
      end if
c
c  Number seems to have terminated.  Have we got a legal number?
c  Not if we terminated in states 1, 2, 6 or 7.
c
      if ( ihave .eq. 1 .or. ihave .eq. 2 .or. 
     &     ihave .eq. 6 .or. ihave .eq. 7 ) then
        ierror = ihave
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'S_TO_R8 - Serious error!'
        write ( *, '(a)' ) '  Illegal or nonnumeric input:'
        write ( *, '(a,a)' ) '    ', s
        return
      end if
c
c  Number seems OK.  Form it.
c
      if ( jtop .eq. 0 ) then
        rexp = 1.0D+00
      else
        if ( jbot .eq. 1 ) then
          rexp = 10.0D+00 ** ( jsgn * jtop )
        else
          rexp = 10.0D+00 ** ( dble ( jsgn * jtop ) / dble ( jbot ) )
        end if
      end if

      dval = dble ( isgn ) * rexp * rtop / rbot

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
      function ts_mult ( u, h, n )

c*********************************************************************72
c
cc TS_MULT...
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    15 November 2009
c
c  Author:
c
c    Original C++ version by Nick Hale.
c    FORTRAN77 version by John Burkardt.
c
c  Parameters:
c
c    Input, double precision U(N+1), ...
c
c    Input, double precision H, ...
c
c    Input, integer N, ...
c
c    Output, double precision TS_MULT, ...
c
      implicit none

      integer n

      double precision h
      double precision hk
      integer k
      double precision ts
      double precision ts_mult
      double precision u(n+1)
      
      ts = 0.0D+00
      hk = 1.0D+00
      do k = 1, n
        ts = ts + u(k+1) * hk
        hk = hk * h
      end do

      ts_mult = ts

      return
      end
