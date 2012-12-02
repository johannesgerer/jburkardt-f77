      subroutine angle_shift ( alpha, beta, gamma )

c*********************************************************************72
c
cc ANGLE_SHIFT shifts angle ALPHA to lie between BETA and BETA+2PI.
c
c  Discussion:
c
c    The input angle ALPHA is shifted by multiples of 2 * PI to lie
c    between BETA and BETA+2*PI.
c
c    The resulting angle GAMMA has all the same trigonometric function
c    values as ALPHA.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    24 October 2007
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, double precision ALPHA, the angle to be shifted.
c
c    Input, double precision BETA, defines the lower endpoint of
c    the angle range.
c
c    Output, double precision GAMMA, the shifted angle.
c
      implicit none

      double precision alpha
      double precision beta
      double precision gamma
      double precision pi
      parameter ( pi = 3.141592653589793D+00 )

      if ( alpha .lt. beta ) then
        gamma = 
     &    beta - mod ( beta - alpha, 2.0D+00 * pi ) + 2.0D+00 * pi
      else
        gamma = beta + mod ( alpha - beta, 2.0D+00 * pi )
      end if

      return
      end
      subroutine angle_shift_deg ( alpha, beta, gamma )

c*********************************************************************72
c
cc ANGLE_SHIFT_DEG shifts angle ALPHA to lie between BETA and BETA+360.
c
c  Discussion:
c
c    The input angle ALPHA is shifted by multiples of 360 to lie
c    between BETA and BETA+360.
c
c    The resulting angle GAMMA has all the same trigonometric function
c    values as ALPHA.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    24 October 2007
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, double precision ALPHA, the angle to be shifted, in degrees.
c
c    Input, double precision BETA, defines the lower endpoint of
c    the angle range.
c
c    Output, double precision GAMMA, the shifted angle.
c
      implicit none

      double precision alpha
      double precision beta
      double precision gamma

      if ( alpha .lt. beta ) then
        gamma = beta - mod ( beta - alpha, 360.0D+00 ) + 360.0D+00
      else
        gamma = beta + mod ( alpha - beta, 360.0D+00 )
      end if

      return
      end
      subroutine angle_to_rgb ( angle, r, g, b )

c*********************************************************************72
c
cc ANGLE_TO_RGB returns a color on the perimeter of the color hexagon.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    07 December 2004
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, double precision ANGLE, the angle in the color hexagon.  
c    The sextants are defined by the following points:
c        0 degrees, 1, 0, 0, red;
c       60 degrees, 1, 1, 0, yellow;
c      120 degrees, 0, 1, 0, green;
c      180 degrees, 0, 1, 1, cyan;
c      240 degrees, 0, 0, 1, blue;
c      300 degrees, 1, 0, 1, magenta.
c
c    Output, double precision R, G, B, RGB specifications for the color 
c    that lies at the given angle, on the perimeter of the color hexagon.  One
c    value will be 1, and one value will be 0.
c
      implicit none

      double precision angle
      double precision angle2
      double precision b
      double precision g
      double precision degrees_to_radians
      parameter ( degrees_to_radians 
     &  = 3.141592653589793D+00 / 180.0D+00 )
      double precision r

      angle = mod ( angle, 360.0D+00 )

      if ( angle .lt. 0.0D+00 ) then
        angle = angle + 360.0D+00
      end if

      if ( angle .le. 60.0D00 ) then

        angle2 = 3.0D+00 * angle / 4.0D+00
        angle2 = degrees_to_radians * angle2
        r = 1.0D+00
        g = dtan ( angle2 )
        b = 0.0D+00

      else if ( angle .le. 120.0D+00 ) then

        angle2 = 3.0D+00 * angle / 4.0D+00
        angle2 = degrees_to_radians * angle2
        r = dcos ( angle2 ) / dsin ( angle2 )
        g = 1.0D+00
        b = 0.0D+00

      else if ( angle .le. 180.0D+00 ) then

        angle2 = 3.0D+00 * ( angle - 120.0D+00 ) / 4.0D+00
        angle2 = degrees_to_radians * angle2
        r = 0.0D+00
        g = 1.0D+00
        b = dtan ( angle2 )

      else if ( angle .le. 240.0D+00 ) then

        angle2 = 3.0D+00 * ( angle - 120.0D+00 ) / 4.0D+00
        angle2 = degrees_to_radians * angle2
        r = 0.0D+00
        g = dcos ( angle2 ) / dsin ( angle2 )
        b = 1.0D+00

      else if ( angle .le. 300.0D+00 ) then

        angle2 = 3.0D+00 * ( angle - 240.0D+00 ) / 4.0D+00
        angle2 = degrees_to_radians * angle2
        r = dtan ( angle2 )
        g = 0.0D+00
        b = 1.0D+00

      else if ( angle .le. 360.0D+00 ) then

        angle2 = 3.0D+00 * ( angle - 240.0D+00 ) / 4.0D+00
        angle2 = degrees_to_radians * angle2
        r = 1.0D+00
        g = 0.0D+00
        b = dcos ( angle2 ) / dsin ( angle2 )

      end if

      return
      end
      subroutine axis_limits ( xmin, xmax, ndivs, pxmin, pxmax, pxdiv, 
     &  nticks )

c*********************************************************************72
c
cc AXIS_LIMITS returns "nice" axis limits for a plot.
c
c  Discussion:
c
c    The routine is given information about the range of a variable, and
c    the number of divisions desired.  It returns suggestions for
c    labeling a plotting axis for the variable, including the
c    starting and ending points, the length of a single division,
c    and a suggested tick marking for the axis.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    01 July 2008
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, double precision XMIN, XMAX, the lower and upper values that
c    must be included on the axis.  XMIN must be less than XMAX.
c
c    Input, integer NDIVS, the number of divisions desired along
c    the axis.
c
c    Output, double precision PXMIN, PXMAX, the recommended lower and upper
c    axis bounds.  It will be the case that PXMIN <= XMIN < XMAX <= PXMAX.
c
c    Output, double precision PXDIV, the recommended size of a single division.
c
c    Output, integer NTICKS, a suggested number of ticks to use,
c    if subdividing each of the NDIVS divisions of the axis.
c
      implicit none

      integer nsteps
      parameter ( nsteps = 5 )

      double precision best
      double precision good
      integer i
      integer ihi
      integer ilo
      integer intlog
      integer iticks(nsteps) 
      data iticks / 5, 4, 4, 5, 5 /
      integer ival
      integer j
      integer ndivs
      integer nticks
      double precision pxmax
      double precision pxmax2
      double precision pxmin
      double precision pxmin2
      double precision pxdiv
      double precision pxdiv2
      double precision r8_log_10
      double precision reldif
      double precision steps(nsteps)
      data steps / 1.0D+00,  2.0D+00,  4.0D+00,  5.0D+00, 10.0D+00 /
      double precision temp
      double precision xmax
      double precision xmin

      if ( xmin .eq. xmax ) then
        xmin = xmin - 0.5D+00
        xmax = xmax + 0.5D+00
      else if ( xmax .lt. xmin ) then
        temp = xmin
        xmin = xmax
        xmax = temp
      end if

      if ( ndivs .le. 0 ) then
        ndivs = 5
      end if
c
c  Set RELDIF, the size of the X interval divided by the largest X.
c
      if ( xmax .ne. xmin ) then
        reldif = ( xmax - xmin ) / max ( abs ( xmax ), abs ( xmin ) )
      else
        reldif = 0.0D+00
      end if
c
c  If RELDIF tells us that XMIN and XMAX are extremely close,
c  do some simple things.
c
      if ( reldif .lt. 0.00001D+00 ) then

        if ( xmax .eq. 0.0D+00 ) then

          pxdiv = 1.0D+00

        else

          intlog = int ( r8_log_10 ( xmax ) )

          if ( intlog .lt. 0 ) then
            intlog = intlog - 1
          end if

          pxdiv = 10.0D+00**intlog

          if ( 1.0D+00 .lt. pxdiv ) then
            pxdiv = 1.0D+00
          end if

        end if

        nticks = 5
        pxmin = xmax - dble ( ndivs / 2 ) * pxdiv
        pxmax = xmax + dble ( ndivs - ( ndivs / 2 ) ) * pxdiv
c
c  But now handle the more general case, when XMIN and XMAX
c  are relatively far apart.
c
      else

        best = -999.0D+00
c
c  On second loop, increase INTLOG by 1.
c
        do j = 1, 2
c
c  Compute INTLOG, roughly the logarithm base 10 of the range
c  divided by the number of divisions.
c
          intlog = int ( r8_log_10 ( ( xmax - xmin ) 
     &      / dble ( ndivs ) ) ) + ( j - 1 )

          if ( xmax - xmin .lt. dble ( ndivs ) ) then
            intlog = intlog - 1
          end if
c
c  Now consider taking 1, 2, 4, 5 or 10 steps of size 10**INTLOG:
c
          do i = 1, nsteps
c
c  Compute the size of each step.
c
            pxdiv2 = steps(i) * 10.0D+00**intlog
c
c  Make sure NDIVS steps can reach from XMIN to XMAX, at least.
c
            if ( xmax .le. xmin + ndivs * pxdiv2 ) then
c
c  Now decide where to start the axis.
c  Start the axis at PXMIN2, to the left of XMIN, and
c  representing a whole number of steps of size PXDIV2.
c
              if ( 0.0D+00 .le. xmin ) then
                ival = int ( xmin / pxdiv2 )
              else
                ival = int ( xmin / pxdiv2 ) - 1
              end if

              pxmin2 = ival * pxdiv2
c
c  PXMAX2 is, of course, NDIVS steps above PXMIN2.
c
              pxmax2 = pxmin2 + ndivs * pxdiv2
c
c  Only consider going on if PXMAX2 is at least XMAX.
c
              if ( xmax .le. pxmax2 ) then
c
c  Now judge this grid by the relative amount of wasted axis length.
c
                good = ( xmax - xmin ) / ( pxmax2 - pxmin2 )

                if ( best .lt. good ) then
                  best = good
                  pxmax = pxmax2
                  pxmin = pxmin2
                  pxdiv = pxdiv2
                  nticks = iticks(i)
                end if

              end if

            end if

          end do
    
        end do

      end if
c
c  If necessary, adjust the locations of PXMIN and PXMAX so that the
c  interval is more symmetric in containing XMIN through XMAX.
c
10    continue

        ilo = int ( xmin - pxmin ) / pxdiv
        ihi = int ( pxmax - xmax ) / pxdiv

        if ( ihi .lt. ilo + 2 ) then
          go to 20
        end if

        pxmin = pxmin - pxdiv
        pxmax = pxmax - pxdiv

      go to 10

20    continue

      return
      end
      subroutine bar_check ( digit, check )

c*********************************************************************72
c
cc BAR_CHECK computes the check digit for a barcode.
c
c  Formula:
c
c    CHECK = SUM ( I = 1, 11, by 2's ) DIGIT(I)
c       + 3 * SUM ( I = 2, 10, by 2's ) DIGIT(I)
c
c    CHECK = MOD ( 10 - MOD ( CHECK, 10 ), 10 )
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    01 August 2008
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer DIGIT(12), entries 1 through 11 of DIGIT 
c    contain the digits of the bar code.  Each entry must be between 0 and 9.
c    The 12th digit should be the check digit.
c
c    Output, integer CHECK, the correct check digit.  If the bar
c    code is correct, then DIGIT(12) should equal CHECK.
c
      implicit none

      integer check
      integer digit(12)
      integer i

      check = 0
      do i = 1, 11, 2
        check = check + digit(i)
      end do

      do i = 2, 10, 2
        check = check + 3 * digit(i)
      end do

      check = mod ( 10 - mod ( check, 10 ), 10 )

      return
      end
      subroutine bar_code ( digit, bar )

c*********************************************************************72
c
cc BAR_CODE constructs the 113 character barcode from 11 digits.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    25 July 2010
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input/output, integer DIGIT(12).
c    On input, the first 11 entries of DIGIT contain a code to be
c    turned into a barcode.
c    On output, the 12-th entry of DIGIT is a check digit.
c
c    Output, character*113 BAR, the bar code corresponding to the
c    digit information.
c
      implicit none

      character*113 bar
      integer check
      character*7 codel
      character*7 coder
      integer digit(12)
      integer i
c
c  9 character quiet zone.
c
      bar(1:9) = '000000000'
c
c  3 character guard pattern.
c
      bar(10:12) = '101'
c
c  7 character product category.
c
      call bar_digit_code_left ( digit(1), codel )
      bar(13:19) = codel
c
c  35 characters contain the 5 digit manufacturer code.
c
      do i = 1, 5
        call bar_digit_code_left ( digit(i+1), codel )
        bar(20+(i-1)*7:20+(i-1)*7+6) = codel
      end do
c
c  Center guard pattern.
c
      bar(55:59) = '01010'
c
c  35 characters contain the 5 digit product code.
c
      do i = 1, 5
        call bar_digit_code_right ( digit(i+6), coder )
        bar(60+(i-1)*7:60+(i-1)*7+6) = coder
      end do
c
c  Compute the check digit.
c
      call bar_check ( digit, check )
      digit(12) = check

      call bar_digit_code_right ( digit(12), coder )
      bar(95:101) = coder
c
c  Guard pattern.
c
      bar(102:104) = '101'
c
c  Quiet zone.
c
      bar(105:113) = '000000000'

      return
      end
      subroutine bar_digit_code_left ( digit, codel )

c*********************************************************************72
c
cc BAR_DIGIT_CODE_LEFT returns the 7 character left bar code for a digit.
c
c  Example:
c
c    DIGIT = 3
c    CODEL = '0111101'
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    25 July 2010
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer DIGIT, the digit, between 0 and 9.
c
c    Output, character*7 CODEL, the left code for the digit.
c
      implicit none

      character*7 codel
      integer digit

      if ( digit .eq. 0 ) then
        codel = '0001101'
      else if ( digit .eq. 1 ) then
        codel = '0011001'
      else if ( digit .eq. 2 ) then
        codel = '0010011'
      else if ( digit .eq. 3 ) then
        codel = '0111101'
      else if ( digit .eq. 4 ) then
        codel = '0100011'
      else if ( digit .eq. 5 ) then
        codel = '0110001'
      else if ( digit .eq. 6 ) then
        codel = '0101111'
      else if ( digit .eq. 7 ) then
        codel = '0111011'
      else if ( digit .eq. 8 ) then
        codel = '0110111'
      else if ( digit .eq. 9 ) then
        codel = '0001011'
      else
        codel = '???????'
      end if

      return
      end
      subroutine bar_digit_code_right ( digit, coder )

c*********************************************************************72
c
cc BAR_DIGIT_CODE_RIGHT returns the 7 character right bar code for a digit.
c
c  Example:
c
c    DIGIT = 3
c    CODER = '1000010'
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    26 August 2004
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer DIGIT, the digit, between 0 and 9.
c
c    Output, character*7 CODER, the right code for the digit.
c
      implicit none

      character*7 coder
      integer digit

      if ( digit .eq. 0 ) then
        coder = '1110010'
      else if ( digit .eq. 1 ) then
        coder = '1100110'
      else if ( digit .eq. 2 ) then
        coder = '1101100'
      else if ( digit .eq. 3 ) then
        coder = '1000010'
      else if ( digit .eq. 4 ) then
        coder = '1011100'
      else if ( digit .eq. 5 ) then
        coder = '1001110'
      else if ( digit .eq. 6 ) then
        coder = '1010000'
      else if ( digit .eq. 7 ) then
        coder = '1000100'
      else if ( digit .eq. 8 ) then
        coder = '1001000'
      else if ( digit .eq. 9 ) then
        coder = '1110100'
      else
        coder = '???????'
      end if

      return
      end
      function bmi_english ( w_lb, h_ft, h_in )

c*********************************************************************72
c
cc BMI_ENGLISH computes the body mass index given English measurements.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    07 December 2004
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, double precision W_LB, the body weight in pounds.
c
c    Input, double precision H_FT, H_IN, the body height in feet and inches
c
c    Output, double precision BMI_ENGLISH, the body mass index.
c
      implicit none

      double precision bmi_english
      double precision bmi_metric
      double precision feet_to_meters
      double precision h_ft
      double precision h_in
      double precision h_m
      double precision pounds_to_kilograms
      double precision w_kg
      double precision w_lb

      w_kg = pounds_to_kilograms ( w_lb )

      h_m = feet_to_meters ( h_ft + ( h_in / 12.0D+00 ) )

      bmi_english = bmi_metric ( w_kg, h_m )

      return
      end
      function bmi_metric ( w_kg, h_m )

c*********************************************************************72
c
cc BMI_METRIC computes the body mass index given metric measurements.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    20 May 2008
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, double precision W_KG, the body weight in kilograms.
c
c    Input, double precision H_M, the body height in meters.
c
c    Output, double precision BMI_METRIC, the body mass index.
c
      implicit none

      double precision bmi_metric
      double precision h_m
      double precision w_kg

      bmi_metric = ( w_kg / h_m ) / h_m

      return
      end
      function ch_is_digit ( c )

c*********************************************************************72
c
cc CH_IS_DIGIT returns TRUE if a character is a decimal digit.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    15 January 1999
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, character C, the character to be analyzed.
c
c    Output, logical CH_IS_DIGIT, TRUE if C is a digit, FALSE otherwise.
c
      implicit none

      character c
      logical ch_is_digit

      if ( lge ( c, '0' ) .and. lle ( c, '9' ) ) then
        ch_is_digit = .true.
      else
        ch_is_digit = .false.
      end if

      return
      end
      function degrees_to_radians ( degrees )

c*********************************************************************72
c
cc DEGREES_TO_RADIANS converts an angle measure from degrees to radians.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    01 January 2007
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, double precision DEGREES, the angle measure in degrees.
c
c    Output, double precision DEGREES_TO_RADIANS, the angle measure in radians.
c
      implicit none

      double precision degrees
      double precision degrees_to_radians
      double precision pi
      parameter ( pi = 3.141592653589793D+00 )

      degrees_to_radians = ( degrees / 180.0D+00 ) * pi

      return
      end
      function e_constant ( )

c*********************************************************************72
c
cc E_CONSTANT returns the value of E.
c
c  Discussion:
c
c    "E" was named in honor of Euler, but is known as Napier's constant.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    30 October 2007
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, double precision E_CONSTANT, the base of the natural 
c    logarithm system.
c
      implicit none

      double precision e_constant

      e_constant = 2.718281828459045D+00

      return
      end
      function euler_constant ( )

c*********************************************************************72
c
cc EULER_CONSTANT returns the value of the Euler-Mascheroni constant.
c
c  Discussion:
c
c    The Euler-Mascheroni constant is often denoted by a lower-case
c    Gamma.  Gamma is defined as
c
c      Gamma = limit ( M -> Infinity ) ( Sum ( 1 <= N <= M ) 1 / N ) - Log ( M )
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    30 October 2007
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, double precision EULER_CONSTANT, the value of the 
c    Euler-Mascheroni constant.
c
      implicit none

      double precision euler_constant

      euler_constant = 0.5772156649015328D+00

      return
      end
      subroutine fac_div ( prime_num, npower1, npower2, npower3 )

c*********************************************************************72
c
cc FAC_DIV divides two quantities represented as prime factors.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    30 October 2007
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer PRIME_NUM, the index of the highest prime 
c    number used in the representations.
c
c    Input, integer NPOWER1(PRIME_NUM), the powers of primes
c    in the representation of the first quantity.
c
c    Input, integer NPOWER2(PRIME_NUM), the powers of primes
c    in the representation of the second quantity.
c
c    Output, integer NPOWER3(PRIME_NUM), the powers of primes
c    in the representation of the quotient.
c
      implicit none

      integer prime_num

      integer npower1(prime_num)
      integer npower2(prime_num)
      integer npower3(prime_num)
      integer prime

      do prime = 1, prime_num

        npower3(prime) = npower1(prime) - npower2(prime)

      end do

      return
      end
      subroutine fac_gcd ( prime_num, npower1, npower2, npower3 )

c*********************************************************************72
c
cc FAC_GCD finds the GCD of two products of prime factors.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    01 August 2008
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer PRIME_NUM, the index of the highest prime 
c    number used in the representations.
c
c    Input, integer NPOWER1(PRIME_NUM), the powers of primes
c    in the representation of the first quantity.  All the powers
c    must be nonnegative.
c
c    Input, integer NPOWER2(PRIME_NUM), the powers of primes
c    in the representation of the second quantity.  All the powers
c    must be nonnegative.
c
c    Output, integer NPOWER3(PRIME_NUM), the powers of primes
c    in the representation of the GCD.
c
      implicit none

      integer prime_num

      integer i
      integer npower1(prime_num)
      integer npower2(prime_num)
      integer npower3(prime_num)

      do i = 1, prime_num

        if ( npower1(i) .lt. 0 ) then
          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) 'FAC_GCD - Fatal error!'
          write ( *, '(a)' ) '  One of the powers is negative!'
          stop
        end if

        if ( npower2(i) .lt. 0 ) then
          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) 'FAC_GCD - Fatal error!'
          write ( *, '(a)' ) '  One of the powers is negative!'
          stop
        end if

        npower3(i) = min ( npower1(i), npower2(i) )

      end do

      return
      end
      subroutine fac_lcm ( prime_num, npower1, npower2, npower3 )

c*********************************************************************72
c
cc FAC_LCM finds the LCM of two products of prime factors.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    22 August 2008
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer PRIME_NUM, the index of the highest prime 
c    number used in the representations.
c
c    Input, integer NPOWER1(PRIME_NUM), the powers of primes
c    in the representation of the first quantity.
c
c    Input, integer NPOWER2(PRIME_NUM), the powers of primes
c    in the representation of the second quantity.
c
c    Output, integer NPOWER3(PRIME_NUM), the powers of primes
c    in the representation of the LCM.
c
      implicit none

      integer prime_num

      integer i
      integer npower1(prime_num)
      integer npower2(prime_num)
      integer npower3(prime_num)

      do i = 1, prime_num

        if ( npower1(i) .lt. 0 ) then
          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) 'FAC_LCM - Fatal error!'
          write ( *, '(a)' ) '  One of the powers is negative!'
          stop
        end if

        if ( npower2(i) .lt. 0 ) then
          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) 'FAC_LCM - Fatal error!'
          write ( *, '(a)' ) '  One of the powers is negative!'
          stop
        end if

        npower3(i) = max ( npower1(i), npower2(i) )

      end do

      return
      end
      subroutine fac_mul ( prime_num, npower1, npower2, npower3 )

c*********************************************************************72
c
cc FAC_MUL multiplies two quantities represented as prime factors.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    22 August 2008
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer PRIME_NUM, the index of the highest prime 
c    number used in the representations.
c
c    Input, integer NPOWER1(PRIME_NUM), the powers of primes
c    in the representation of the first quantity.
c
c    Input, integer NPOWER2(PRIME_NUM), the powers of primes
c    in the representation of the second quantity.
c
c    Output, integer NPOWER3(PRIME_NUM), the powers of primes
c    in the representation of the product.
c
      implicit none

      integer prime_num

      integer npower1(prime_num)
      integer npower2(prime_num)
      integer npower3(prime_num)
      integer prime

      do prime = 1, prime_num
        npower3(prime) = npower1(prime) + npower2(prime)
      end do

      return
      end
      subroutine fac_print ( prime_num, npower )

c*********************************************************************72
c
cc FAC_PRINT prints a product of prime factors.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    22 August 2008
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer PRIME_NUM, the index of the highest prime 
c    number used in the representations.
c
c    Input, integer NPOWER(PRIME_NUM), the powers of primes
c    in the representation of the quantity.
c
      implicit none

      integer prime_num

      integer i
      integer npower(prime_num)
      integer prime

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '   Prime     Power'
      write ( *, '(a)' ) ' '
      do i = 1, prime_num
        if ( npower(i) .ne. 0 ) then
          write ( *, '(i8,2x,i8)' ) prime(i), npower(i)
        end if
      end do

      return
      end
      subroutine fac_to_i4 ( prime_num, npower, intval )

c*********************************************************************72
c
cc FAC_TO_I4 converts a product of prime factors into an integer.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    22 August 2008
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer PRIME_NUM, the index of the highest prime 
c    number used in the representations.
c
c    Input, integer NPOWER(PRIME_NUM), the powers of primes
c    in the representation of the quantity.  If any of these powers
c    are negative, then INTVAL will be set to 0.
c
c    Output, integer INTVAL, the integer represented by the 
c    product of the prime factors.
c
      implicit none

      integer prime_num

      integer i
      integer intval
      integer npower(prime_num)
      integer prime

      intval = 1
      do i = 1, prime_num

        if ( npower(i) .lt. 0 ) then
          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) 'FAC_TO_I4 - Fatal error!'
          write ( *, '(a)' ) '  One of the powers is negative!'
          stop
        end if

        intval = intval * prime(i)**npower(i)

      end do

      return
      end
      subroutine fac_to_rat ( prime_num, npower, top, bot )

c*********************************************************************72
c
cc FAC_TO_RAT converts a prime factorization into a rational value.
c
c  Example:
c
c    Start with the prime factorization representation:
c
c      40/9 = 2**3 * 3**(-2) * 5
c
c    Input:
c
c      NPOWER = ( 3, -2, 1 )
c
c    Output:
c
c      TOP = 40 ( = 2**3 * 5**1 = PRIME(1)**3                 * PRIME(3)**1 )
c      BOT = 9  ( = 3**2        =               PRIME(2)**2 )
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    16 August 1999
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer PRIME_NUM, the index of the highest prime 
c    number used in the representations.
c
c    Input, integer NPOWER(PRIME_NUM).  NPOWER(I) is the power of
c    the I-th prime in the prime factorization.  NPOWER(I) may
c    be positive or negative.
c
c    Output, integer TOP, BOT, the top and bottom of a rational 
c    value.
c
      implicit none

      integer prime_num

      integer bot
      integer i
      integer npower(prime_num)
      integer prime
      integer top

      top = 1
      bot = 1
      do i = 1, prime_num
        if ( 0 .lt. npower(i) ) then
          top = top * prime(i)**npower(i)
        else if ( npower(i) .lt. 0 ) then
          bot = bot * prime(i)**(-npower(i))
        end if
      end do

      return
      end
      function feet_to_meters ( ft )

c*********************************************************************72
c
cc FEET_TO_METERS converts a measurement in feet to meters.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    20 May 2008
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, double precision FT, the length in feet.
c
c    Output, double precision FEET_TO_METERS, the corresponding 
c    length in meters.
c
      implicit none

      double precision feet_to_meters
      double precision ft

      feet_to_meters = 0.0254D+00 * 12.0D+00 * ft

      return
      end
      function gauss_sum ( dim_num, n, amplitude, center, width, x )
 
c*********************************************************************72
c
cc GAUSS_SUM evaluates a function that is the sum of Gaussians.
c
c  Discussion:
c
c    Gauss_Sum(X) = Sum ( 1 <= J <= Ngauss ) Amplitude(I) * exp ( -Arg )
c
c    where
c
c      Arg = sum ( 1 <= I <= DIM_NUM ) ( ( ( X(I) - Center(I,J) ) / Width(J) )**2 )
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    22 August 2008
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer DIM_NUM, the spatial dimension.
c
c    Input, integer N, the number of component Gaussian functions.
c
c    Input, double precision AMPLITUDE(N), CENTER(DIM_NUM,N), WIDTH(N),
c    the amplitude, center and width for the component Gaussian functions.
c
c    Input, double precision X(DIM_NUM), the point at which the function 
c    is to be evaluated.
c
c    Output, double precision GAUSS_SUM, the value of the function.
c
      implicit none

      integer dim_num
      integer n

      double precision amplitude(n)
      double precision arg
      double precision center(dim_num,n)
      double precision gauss_sum
      integer i
      integer j
      double precision width(n)
      double precision x(dim_num)

      gauss_sum = 0.0D+00

      do j = 1, n

        arg = 0.0D+00
        do i = 1, dim_num
          arg = arg + ( ( x(i) - center(i,j) ) / width(j) )**2
        end do

        gauss_sum = gauss_sum + amplitude(j) * exp ( - arg )

      end do

      return
      end
      subroutine get_seed ( seed )

c*********************************************************************72
c
cc GET_SEED returns a seed for the random number generator.
c
c  Discussion:
c
c    The seed depends on the current time, and ought to be (slightly)
c    different every millisecond.  Thus, calling this routine several
c    times in succession will probably return the SAME seed, but
c    calling it a few minutes or days apart will turn a suitably
c    "random" seed.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    11 January 2007
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, integer SEED, a pseudorandom seed value.
c
      implicit none

      integer day
      integer hour
      integer i4_huge
      parameter ( i4_huge = 2147483647 )
      integer milli
      integer minute
      integer month
      integer second
      integer seed
      double precision temp
      character * ( 10 ) time
      character * ( 8 ) date
      integer year

      call date_and_time ( date, time )

      read ( date, '(i4,i2,i2)' ) year, month, day
      read ( time, '(i2,i2,i2,1x,i3)' ) hour, minute, second, milli

      temp = 0.0D+00
      temp = temp + dble ( month - 1 ) / 11.0D+00
      temp = temp + dble ( day   - 1 ) / 30.0D+00
      temp = temp + dble ( hour      ) / 23.0D+00
      temp = temp + dble ( minute    ) / 59.0D+00
      temp = temp + dble ( second    ) / 59.0D+00
      temp = temp + dble ( milli     ) / 999.0D+00

      temp = temp / 6.0D+00
c
c  Force 0 < TEMP <= 1.
c
10    continue

      if ( temp .le. 0.0D+00 ) then
        temp = temp + 1.0D+00
        go to 10
      end if

20    continue

      if ( 1.0D+00 .lt. temp ) then
        temp = temp - 1.0D+00
        go to 20
      end if

      seed = int ( dble ( i4_huge ) * temp )
c
c  Never use a seed of 0 or maximum integer.
c
      if ( seed .eq. 0 ) then
        seed = 1
      end if

      if ( seed .eq. i4_huge ) then
        seed = seed - 1
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

        if ( i /= 5 .and. i /= 6 .and. i /= 9 ) then

          open ( unit = i, err = 10, status = 'scratch' )
          close ( unit = i )

          unit = i

          return
        end if

10      continue

      end do

      return
      end
      subroutine grid1 ( dim_num, nstep, x1, x2, x )

c*********************************************************************72
c
cc GRID1 finds grid points between X1 and X2 in N dimensions.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    22 August 2008
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer DIM_NUM, the dimension of the points X1 and X2.
c
c    Input, integer NSTEP, the number of points to be generated.
c    NSTEP must be at least 2.
c
c    Input, double precision X1(DIM_NUM), X2(DIM_NUM), the first and last
c    points, between which the equally spaced points are
c    to be computed.
c
c    Output, double precision X(DIM_NUM,NSTEP), the set of equally spaced
c    points.  Each column of X represents one point, with X(*,1) = X1
c    and X(*,NSTEP) = X2.
c
      implicit none

      integer dim_num
      integer nstep

      integer i
      integer j
      double precision x(dim_num,nstep)
      double precision x1(dim_num)
      double precision x2(dim_num)

      if ( dim_num .lt. 1 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'GRID1 - Fatal error!'
        write ( *, '(a)' ) '  DIM_NUM < 1.'
        write ( *, '(a,i8)' ) '  DIM_NUM = ', dim_num
        stop
      end if

      if ( nstep .lt. 2 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'GRID1 - Fatal error!'
        write ( *, '(a)' ) '  NSTEP < 2.'
        write ( *, '(a,i8)' ) '  NSTEP = ', nstep
        stop
      end if

      do j = 1, nstep
        do i = 1, dim_num
          x(i,j) = 
     &      ( dble ( nstep - j     ) * x1(i)   
     &      + dble (         j - 1 ) * x2(i) ) 
     &      / dble ( nstep     - 1 )
        end do
      end do

      return
      end
      subroutine grid1n ( j, dim_num, nstep, x1, x2, x )

c*********************************************************************72
c
cc GRID1N finds the I-th grid point between X1 and X2 in N dimensions.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    01 September 2008
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer J, the number of the desired point.
c    Normally J would be between 1 and NSTEP, but that is
c    not necessary.  Note that J = 1 returns X1 and J = NSTEP
c    returns X2.
c
c    Input, integer DIM_NUM, the dimension of the points X, 
c    X1 and X2.
c
c    Input, integer NSTEP, this is the number of equally
c    spaced points that are between X1 and X2.  NSTEP must
c    be at least 2, because X1 and X2 are always included
c    in the set of points.
c
c    Input, double precision X1(DIM_NUM), X2(DIM_NUM), the first and last
c    points, between which the equally spaced points lie.
c
c    Output, double precision X(DIM_NUM), the J-th grid point between X1
c    and X2.
c
      implicit none

      integer dim_num
      integer nstep

      integer dim
      integer j
      double precision x(dim_num)
      double precision x1(dim_num)
      double precision x2(dim_num)

      if ( dim_num .lt. 1 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'GRID1N - Fatal error!'
        write ( *, '(a)' ) '  DIM_NUM < 1.'
        write ( *, '(a,i8)' ) '  DIM_NUM = ', dim_num
        stop
      end if

      if ( nstep .le. 1 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'GRID1N - Fatal error!'
        write ( *, '(a)' ) '  NSTEP <= 1.'
        write ( *, '(a,i8)' ) '  NSTEP = ', nstep
        stop
      end if

      do dim = 1, dim_num
        x(dim) = ( dble ( nstep - j     ) * x1(dim) 
     &           + dble (         j - 1 ) * x2(dim) ) 
     &           / dble ( nstep     - 1 )
      end do

      return
      end
      subroutine grid2 ( j1, j2, dim_num, nstep, x1, x2, x )

c*********************************************************************72
c
cc GRID2 computes grid points between X1 and X2 in N dimensions.
c
c  Discussion:
c
c    GRID2 computes grid points between X1 and X2 in N dimensions.
c
c    However, X1 need not be the first point computed, nor X2 the last.
c    The user must specify the steps on which X1 and X2 are passed
c    through.  These steps may even be outside the range of 1 through NSTEP.
c
c    We assume that a set of equally spaced points have
c    been drawn on the line through X1 and X2, and that
c    they have been numbered, with X1 labeled J1 and X2
c    labeled J2.  J1 or J2 may be between 1 and NSTEP,
c    in which case X1 or X2 will actually be returned in the
c    X array, but there is no requirement that J1 or J2
c    satisfy this condition.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    06 December 2008
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer J1, J2.  J1 specifies the step on which
c    X1 would be computed, and similarly for J2.  
c    J1 and J2 must be distinct.
c
c    Input, integer DIM_NUM, the dimension of the points X1 and X2.
c
c    Input, integer NSTEP, this is the number of equally
c    spaced points that are to be generated.
c    NSTEP should be at least 1.
c
c    Input, double precision X1(DIM_NUM), X2(DIM_NUM), the points that define
c    the line along which the equally spaced points are generated, and
c    which may or may not be included in the set of computed points.
c
c    Output, double precision X(DIM_NUM,NSTEP), the set of equally spaced
c    points.  Each column of X represents one point.
c    If 1 <= J1 <= NSTEP, then X(*,J1) = X1, and similarly for J2.
c
      implicit none

      integer dim_num
      integer nstep

      integer i
      integer j
      integer j1
      integer j2
      double precision x(dim_num,nstep)
      double precision x1(dim_num)
      double precision x2(dim_num)

      if ( dim_num .lt. 1 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'GRID2 - Fatal error!'
        write ( *, '(a)' ) '  DIM_NUM < 1.'
        write ( *, '(a,i8)' ) '  DIM_NUM = ', dim_num
        stop
      end if

      if ( j1 .eq. j2 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'GRID2 - Fatal error!'
        write ( *, '(a)' ) '  J1 = J2, leading to zero denominator.'
        write ( *, '(a,i8)' ) '  J1 = ', j1
        write ( *, '(a,i8)' ) '  J2 = ', j2
        stop
      end if

      do j = 1, nstep
        do i = 1, dim_num
          x(i,j) = ( dble ( j2 - j      ) * x1(i)   
     &             + dble (      j - j1 ) * x2(i) ) 
     &             / dble ( j2     - j1 )
        end do
      end do

      return
      end
      subroutine grid2n ( j, j1, j2, dim_num, x1, x2, x )

c*********************************************************************72
c
cc GRID2N computes one grid point between X1 and X2 in N dimensions.
c
c  Discussion:
c
c    However, X1 need not be the first point computed, nor X2 the last.
c    The user must specify the steps on which X1 and X2 are passed through.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    04 June 2010
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer J, the J coordinate of the desired point.
c    Note that if J = J1, X will be returned as X1, and if
c    J = J2, X will be returned as X2.
c
c    Input, integer J1, J2.  J1 specifies the step on which
c    X1 would be computed, and similarly for J2.  That is,
c    we assume that a set of equally spaced points have
c    been drawn on the line through X1 and X2, and that
c    they have been numbered, with X1 labeled J1 and X2
c    labeled J2.  J1 and J2 must be distinct.
c
c    Input, integer DIM_NUM, the dimension of the points X1 and X2.
c
c    Input, double precision X1(DIM_NUM), X2(DIM_NUM), the points that define
c    the line along which the equally spaced points are
c    generated, and which may or may not be included in the
c    set of computed points.
c
c    Output, double precision X(DIM_NUM).  X(I) is the J-th point from the
c    set of equally spaced points.
c
      implicit none

      integer dim_num

      integer i
      integer j
      integer j1
      integer j2
      double precision x(dim_num)
      double precision x1(dim_num)
      double precision x2(dim_num)

      if ( dim_num .lt. 1 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'GRID2N - Fatal error!'
        write ( *, '(a)' ) '  DIM_NUM .lt. 1.'
        write ( *, '(a,i8)' ) '  DIM_NUM = ', dim_num
        stop
      end if

      if ( j1 .eq. j2 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'GRID2N - Fatal error!'
        write ( *, '(a)' ) '  J1 = J2, leading to zero denominator.'
        write ( *, '(a,i8)' ) '  J1 = ', j1
        write ( *, '(a,i8)' ) '  J2 = ', j2
        stop
      end if

      do i = 1, dim_num
        x(i) = ( dble ( j2 - j      ) * x1(i)   
     &         + dble (      j - j1 ) * x2(i) ) 
     &         / dble ( j2     - j1 )
      end do

      return
      end
      subroutine grid3 ( dim_num, nstep1, nstep2, x1, x2, x3, x )

c*********************************************************************72
c
cc GRID3 computes a grid on the parallelogram set by X1, X2 and X3 in N space.
c
c  Discussion:
c
c    The line between X1 and X2 will have NSTEP1 points generated along 
c    it, and the line between X1 and X3 will have NSTEP2 points generated
c    along it.
c
c    Fixing the second and third indices of X represents one point, with
c    the following special values:
c
c      X(*,1,1)      = X1
c      X(*,NSTEP1,1) = X2
c      X(*,1,NSTEP2) = X3.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    04 June 2010
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer DIM_NUM, the dimension of the points X1, 
c    X2 and X3.
c
c    Input, integer NSTEP1, NSTEP2.  These are the number of
c    equally spaced points to generate in the first and second
c    directions.  NSTEP1 and NSTEP2 must be at least 2, because X1, X2 and
c    X3 are always included in the set of points.
c
c    Input, double precision X1(DIM_NUM), X2(DIM_NUM), X3(DIM_NUM), the points
c    which define three corners of the parallelogram on
c    which the grid will be generated.
c
c    Output, double precision X(DIM_NUM,NSTEP1,NSTEP2), the set of equally
c    spaced points.  
c
      implicit none

      integer dim_num
      integer nstep1
      integer nstep2

      integer i
      integer j
      integer k
      double precision psi1
      double precision psi2
      double precision psi3
      double precision x(dim_num,nstep1,nstep2)
      double precision x1(dim_num)
      double precision x2(dim_num)
      double precision x3(dim_num)

      if ( dim_num .lt. 1 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'GRID3 - Fatal error!'
        write ( *, '(a)' ) '  DIM_NUM .lt. 1.'
        write ( *, '(a,i8)' ) '  DIM_NUM = ', dim_num
        stop
      end if

      if ( nstep1 .le. 1 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'GRID3 - Fatal error!'
        write ( *, '(a)' ) '  NSTEP1 .le. 1.'
        write ( *, '(a,i8)' ) '  NSTEP1 = ', nstep1
        stop
      end if

      if ( nstep2 .le. 1 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'GRID3 - Fatal error!'
        write ( *, '(a)' ) '  NSTEP2 .le. 1.'
        write ( *, '(a,i8)' ) '  NSTEP2 = ', nstep2
        stop
      end if

      do j = 1, nstep1

        psi2 = dble ( j - 1 ) / dble ( nstep1 - 1 )

        do k = 1, nstep2

          psi3 = dble ( k - 1 ) / dble ( nstep2 - 1 )

          psi1 = 1.0D+00 - psi2 - psi3

          do i = 1, dim_num
            x(i,j,k) = psi1 * x1(i) + psi2 * x2(i) + psi3 * x3(i)
          end do

        end do
      end do

      return
      end
      subroutine grid3n ( j, k, dim_num, nstep1, nstep2, x1, x2, x3, x )

c*********************************************************************72
c
cc GRID3N computes a parallelogram grid on 3 points in N dimensions.
c
c  Discussion:
c
c    The line between X1 and X2 will have NSTEP1
c    points generated along it, and the line between X1 and
c    X3 will have NSTEP2 points generated along it.
c
c    The following special values are:
c
c      J       K         X
c
c      1       1         X1
c      NSTEP1  1         X2
c      1       NSTEP2    X3
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    04 June 2010
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer J, K, the parallelogram coordinates
c    of the point.  J measures steps from X1 to X2, and
c    K measures steps from X1 to X3.  Normally, J would
c    be between 1 and NSTEP1, K between 1 and NSTEP2,
c    but this is not necessary.
c
c    Input, integer DIM_NUM, the dimension of the points X1, 
c    X2 and X3.
c
c    Input, integer NSTEP1, NSTEP2.  These are the number of
c    equally spaced points to generate in the first and second
c    directions.  NSTEP1 and NSTEP2 must be at least 2, because X1, X2 and
c    X3 are always included in the set of points.
c
c    Input, double precision X1(DIM_NUM), X2(DIM_NUM), X3(DIM_NUM), the points
c    which define three corners of the parallelogram on
c    which the grid will be generated.
c
c    Output, double precision X(DIM_NUM), the point with coordinates (J,K)
c    from the the set of equally  spaced points.  
c
      implicit none

      integer dim_num
      integer nstep1
      integer nstep2

      integer i
      integer j
      integer k
      double precision psi1
      double precision psi2
      double precision psi3
      double precision x(dim_num)
      double precision x1(dim_num)
      double precision x2(dim_num)
      double precision x3(dim_num)

      if ( dim_num .lt. 1 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'GRID3N - Fatal error!'
        write ( *, '(a)' ) '  DIM_NUM .lt. 1.'
        write ( *, '(a,i8)' ) '  DIM_NUM = ', dim_num
        stop
      end if

      if ( nstep1 .le. 1 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'GRID3N - Fatal error!'
        write ( *, '(a)' ) '  NSTEP1 .le. 1.'
        write ( *, '(a,i8)' ) '  NSTEP1 = ', nstep1
        stop
      end if

      if ( nstep2 .le. 1 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'GRID3N - Fatal error!'
        write ( *, '(a)' ) '  NSTEP2 .le. 1.'
        write ( *, '(a,i8)' ) '  NSTEP2 = ', nstep2
        stop
      end if

      psi2 = dble ( j - 1 ) / dble ( nstep1 - 1 )

      psi3 = dble ( k - 1 ) / dble ( nstep2 - 1 )

      psi1 = 1.0D+00 - psi2 - psi3

      do i = 1, dim_num
        x(i) = psi1 * x1(i) + psi2 * x2(i) + psi3 * x3(i)
      end do

      return
      end
      subroutine grid4 ( j1, j2, k1, k2, dim_num, nstep1, nstep2, 
     &  x1, x2, x3, x )

c*********************************************************************72
c
cc GRID4 computes a grid on the parallelogram set by X1, X2 and X3 in N space.
c
c  Discussion:
c
c    Unlike GRID3, GRID4 does not necessarily place X1 at the
c    "origin" of the parallelogram, with X2 and X3 set at the
c    extreme J and K coordinates.  Instead, the user is free
c    to specify the J and K coordinates of the points, although
c    they are required to lie on a subparallelogram of the
c    larger one.
c
c    The line through X1 and X2 will have NSTEP1
c    points generated along it, and the line through X1 and
c    X3 will have NSTEP2 points generated along it.
c
c    If we imagine that the
c    main parallelogram is drawn first, with coordinate
c    ranges 1 .le. J .le. NSTEP1 and 1 .le. K .le. NSTEP2, then
c    these indices determine the (J,K) coordinates of the
c    three points, namely:
c
c      X1 : (J1,K1)
c      X2 : (J2,K1)
c      X3 : (J1,K2)
c
c    Of course, we actually start with the points X1, X2,
c    and X3, and they define a parallelogram and a (J,K)
c    coordinate system over the plane containing them.  We
c    then are free to consider the parallelogram defined
c    by the three points (1,1), (NSTEP1,1) and (1,NSTEP2),
c    which may or may not contain any of the points X1, X2
c    and X3.
c
c    Assuming that the indices J1, J2, K1 and K2 are "within
c    bounds", the following special values will be computed:
c
c      X(*,J1,K1) = X1
c      X(*,J2,K1) = X2
c      X(*,J1,K2) = X3.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    04 June 2010
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer J1, J2, K1, K2, the indices.  
c
c    Input, integer DIM_NUM, the dimension of the points X1, 
c    X2 and X3.
c
c    Input, integer NSTEP1, NSTEP2.  These are the number of
c    equally spaced points to generate in the first and second
c    directions.  NSTEP1 and NSTEP2 should be at least 1.
c
c    Input, double precision X1(DIM_NUM), X2(DIM_NUM), X3(DIM_NUM), the points
c    which define three corners of the parallelogram on
c    which the grid will be generated.
c
c    Output, double precision X(DIM_NUM,NSTEP1,NSTEP2), the set of equally
c    spaced points.  Fixing the second and third indices
c    of X represents one point.
c
      implicit none

      integer dim_num
      integer nstep1
      integer nstep2

      integer i
      integer j
      integer j1
      integer j2
      integer k
      integer k1
      integer k2
      double precision psi1
      double precision psi2
      double precision psi3
      double precision x(dim_num,nstep1,nstep2)
      double precision x1(dim_num)
      double precision x2(dim_num)
      double precision x3(dim_num)

      if ( dim_num .lt. 1 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'GRID4 - Fatal error!'
        write ( *, '(a)' ) '  DIM_NUM .lt. 1.'
        write ( *, '(a,i8)' ) '  DIM_NUM = ', dim_num
        stop
      end if

      if ( nstep1 .le. 1 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'GRID4 - Fatal error!'
        write ( *, '(a)' ) '  NSTEP1 .le. 1.'
        write ( *, '(a,i8)' ) '  NSTEP1 = ', nstep1
        stop
      end if

      if ( nstep2 .le. 1 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'GRID4 - Fatal error!'
        write ( *, '(a)' ) '  NSTEP2 .le. 1.'
        write ( *, '(a,i8)' ) '  NSTEP2 = ', nstep2
        stop
      end if

      if ( j1 .eq. j2 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'GRID4 - Fatal error!'
        write ( *, '(a)' ) '  J1 = J2'
        write ( *, '(a,i8)' ) '  J1 = ', j1
        write ( *, '(a,i8)' ) '  J2 = ', j2
        stop
      end if

      if ( k1 .eq. k2 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'GRID4 - Fatal error!'
        write ( *, '(a)' ) '  K1 = K2'
        write ( *, '(a,i8)' ) '  K1 = ', k1
        write ( *, '(a,i8)' ) '  K2 = ', k2
        stop
      end if

      do j = 1, nstep1

        psi2 = dble (  j - j1 ) / dble ( j2 - j1 )

        do k = 1, nstep2

          psi3 = dble (  k - k1 ) / dble ( k2 - k1 )

          psi1 = 1.0D+00 - psi2 - psi3

          do i = 1, dim_num
            x(i,j,k) = psi1 * x1(i) + psi2 * x2(i) + psi3 * x3(i)
          end do

        end do
      end do

      return
      end
      subroutine grid4n ( j, j1, j2, k, k1, k2, dim_num, nstep1, 
     &  nstep2, x1, x2, x3, x )

c*********************************************************************72
c
cc GRID4N computes a single point on a parallelogram grid in N space.
c
c  Discussion:
c
c    The computation is identical to that of GRID4, except that
c    only one point at a time is computed.
c
c    The line through X1 and X2 will have NSTEP1
c    points generated along it, and the line through X1 and
c    X3 will have NSTEP2 points generated along it.
c
c    The following special values will be computed:
c
c      J  K  X
c
c      J1 K1 X1
c      J2 K2 X2
c      J1 K2 X3
c
c    If we imagine that the main parallelogram is drawn first, with 
c    coordinate ranges 1 .le. J .le. NSTEP1 and 1 .le. K .le. NSTEP2, then
c    the indices J and K determine the (J,K) coordinates of the
c    three points X1, X2, and X3, namely:
c
c      X1 : (J1,K1)
c      X2 : (J2,K1)
c      X3 : (J1,K2)
c
c    Of course, we actually start with the points X1, X2,
c    and X3, and they define a parallelogram and an (J,K)
c    coordinate system over the plane containing them.  We
c    then are free to consider the parallelogram defined
c    by the three points (1,1), (NSTEP1,1) and (1,NSTEP2),
c    which may or may not contain any of the points X1, X2
c    and X3.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    04 June 2010
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer J, the J coordinate of the point X.
c
c    Input, integer J1, J2.  See discussion.
c
c    Input, integer K, the K coordinate of the point X.
c
c    Input, integer K1, K2.  See discussion.
c
c    Input, integer DIM_NUM, the dimension of the points 
c    X, X1, X2 and X3.
c
c    Input, integer NSTEP1, NSTEP2.  These are the number of
c    equally spaced points generated in the first and second
c    directions.
c    NSTEP1 and NSTEP2 should be at least 1.
c
c    Input, double precision X1(DIM_NUM), X2(DIM_NUM), X3(DIM_NUM), the points
c    which define three corners of the parallelogram on
c    which the grid will be generated.
c
c    Output, double precision X(DIM_NUM), the point whose parallelogram
c    coordinates are (J,K).
c
      implicit none

      integer dim_num
      integer nstep1
      integer nstep2

      integer i
      integer j
      integer j1
      integer j2
      integer k
      integer k1
      integer k2
      double precision psi1
      double precision psi2
      double precision psi3
      double precision x(dim_num)
      double precision x1(dim_num)
      double precision x2(dim_num)
      double precision x3(dim_num)

      if ( dim_num .lt. 1 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'GRID4N - Fatal error!'
        write ( *, '(a)' ) '  DIM_NUM .lt. 1.'
        write ( *, '(a,i8)' ) '  DIM_NUM = ', dim_num
        stop
      end if

      if ( nstep1 .le. 1 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'GRID4N - Fatal error!'
        write ( *, '(a)' ) '  NSTEP1 .le. 1.'
        write ( *, '(a,i8)' ) '  NSTEP1 = ', nstep1
        stop
      end if

      if ( nstep2 .le. 1 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'GRID4N - Fatal error!'
        write ( *, '(a)' ) '  NSTEP2 .le. 1.'
        write ( *, '(a,i8)' ) '  NSTEP2 = ', nstep2
        stop
      end if

      if ( j1 .eq. j2 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'GRID4N - Fatal error!'
        write ( *, '(a)' ) '  J1 = J2'
        write ( *, '(a,i8)' ) '  J1 = ', j1
        write ( *, '(a,i8)' ) '  J2 = ', j2
        stop
      end if

      if ( k1 .eq. k2 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'GRID4N - Fatal error!'
        write ( *, '(a)' ) '  K1 = K2'
        write ( *, '(a,i8)' ) '  K1 = ', k1
        write ( *, '(a,i8)' ) '  K2 = ', k2
        stop
      end if

      psi2 = dble ( j  - j1 ) / dble ( j2 - j1 )

      psi3 = dble ( k  - k1 ) / dble ( k2 - k1 )

      psi1 = 1.0D+00 - psi2 - psi3

      do i = 1, dim_num
        x(i) = psi1 * x1(i) + psi2 * x2(i) + psi3 * x3(i)
      end do

      return
      end
      function i4_modp ( i, j )

c*********************************************************************72
c
cc I4_MODP returns the nonnegative remainder of integer division.
c
c  Discussion:
c
c    If
c      NREM = I4_MODP ( I, J )
c      NMULT = ( I - NREM ) / J
c    then
c      I = J * NMULT + NREM
c    where NREM is always nonnegative.
c
c    The MOD function computes a result with the same sign as the
c    quantity being divided.  Thus, suppose you had an angle A,
c    and you wanted to ensure that it was between 0 and 360.
c    Then mod(A,360) would do, if A was positive, but if A
c    was negative, your result would be between -360 and 0.
c
c    On the other hand, I4_MODP(A,360) is between 0 and 360, always.
c
c  Example:
c
c        I     J     MOD I4_MODP    Factorization
c
c      107    50       7       7    107 =  2 *  50 + 7
c      107   -50       7       7    107 = -2 * -50 + 7
c     -107    50      -7      43   -107 = -3 *  50 + 43
c     -107   -50      -7      43   -107 =  3 * -50 + 43
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    30 December 2006
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer I, the number to be divided.
c
c    Input, integer J, the number that divides I.
c
c    Output, integer I4_MODP, the nonnegative remainder when I is
c    divided by J.
c
      implicit none

      integer i
      integer i4_modp
      integer j
      integer value

      if ( j .eq. 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'I4_MODP - Fatal error!'
        write ( *, '(a,i8)' ) '  Illegal divisor J = ', j
        stop
      end if

      value = mod ( i, j )

      if ( value .lt. 0 ) then
        value = value + abs ( j )
      end if

      i4_modp = value

      return
      end
      function i4_sign ( x )

c*********************************************************************72
c
cc I4_SIGN evaluates the sign of an I4.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    31 December 2006
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer X, the number whose sign is desired.
c
c    Output, integer I4_SIGN, the sign of the number.
c
      implicit none

      integer i4_sign
      integer x

      if ( x .lt. 0 ) then
        i4_sign = -1
      else
        i4_sign = +1
      end if

      return
      end
      subroutine i4_swap ( i, j )

c*********************************************************************72
c
cc I4_SWAP switches two I4's.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    06 January 2006
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input/output, integer I, J.  On output, the values of I and
c    J have been interchanged.
c
      implicit none

      integer i
      integer j
      integer k

      k = i
      i = j
      j = k

      return
      end
      subroutine i4_to_digits_decimal ( i, n, digit )

c*********************************************************************72
c
cc I4_TO_DIGITS_DECIMAL determines the last N decimal digits of an I4.
c
c  Discussion:
c
c    An I4 is an integer value.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    30 October 2007
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer I, the integer to be analyzed.
c
c    Input, integer N, the number of digits to determine.
c
c    Output, integer DIGIT(N), the last N decimal digits of I.
c    DIGIT(I) is the "coefficient" of 10**(I-1).
c
      implicit none

      integer n

      integer digit(n)
      integer i
      integer i_copy
      integer i4_ten
      parameter ( i4_ten = 10 )
      integer j

      i_copy = i

      do j = 1, n
        digit(j) = mod ( i_copy, i4_ten )
        i_copy = ( i_copy - digit(j) ) / 10
      end do

      return
      end
      subroutine i4_to_fac ( intval, prime_num, npower )

c*********************************************************************72
c
cc I4_TO_FAC converts an I4 into a product of prime factors.
c
c  Discussion:
c
c    This routine will fail if the input integer is not positive,
c    or if PRIME_NUM is too small to account for the factors of the integer.
c
c    An I4 is an integer value.
c
c    The formula is:
c
c      INTVAL = Product ( 1 <= I <= PRIME_NUM ) PRIME(I)**NPOWER(I).
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    30 October 2007
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer INTVAL, the integer to be factored.
c
c    Input, integer PRIME_NUM, the number of prime factors for
c    which storage has been allocated.
c
c    Output, integer NPOWER(PRIME_NUM), the powers of the primes.
c
      implicit none

      integer prime_num

      integer i
      integer intcopy
      integer intval
      integer npower(prime_num)
      integer p
      integer prime

      if ( intval .le. 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'I4_TO_FAC - Fatal error!'
        write ( *, '(a)' ) '  Input integer is not positive.'
        stop
      end if
c
c  Try dividing the remainder by each prime.
c
      intcopy = intval

      do i = 1, prime_num

        npower(i) = 0

        p = prime ( i )

10      continue

        if ( mod ( intcopy, p ) .eq. 0 ) then
          npower(i) = npower(i) + 1
          intcopy = intcopy / p
          go to 10
        end if

      end do

      return
      end
      function i4_to_isbn ( i )

c*********************************************************************72
c
cc I4_TO_ISBN converts an I4 to an ISBN digit.
c
c  Discussion:
c
c    Only the integers 0 through 10 can be input.  The representation
c    of 10 is 'X'.
c
c    An I4 is an integer value.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    30 October 2007
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Book Industry Study Group,
c    The Evolution in Product Identification:
c    Sunrise 2005 and the ISBN-13,
c    http://www.bisg.org/docs/The_Evolution_in_Product_ID.pdf
c
c  Parameters:
c
c    Input, integer I, an integer between 0 and 10.
c
c    Output, character I4_TO_ISBN, the ISBN character code of the integer.
c    If I is illegal, then I4_TO_ISBN is set to '?'.
c
      implicit none

      integer i
      character i4_to_isbn

           if ( i .eq. 0 ) then
        i4_to_isbn = '0'
      else if ( i .eq. 1 ) then
        i4_to_isbn = '1'
      else if ( i .eq. 2 ) then
        i4_to_isbn = '2'
      else if ( i .eq. 3 ) then
        i4_to_isbn = '3'
      else if ( i .eq. 4 ) then
        i4_to_isbn = '4'
      else if ( i .eq. 5 ) then
        i4_to_isbn = '5'
      else if ( i .eq. 6 ) then
        i4_to_isbn = '6'
      else if ( i .eq. 7 ) then
        i4_to_isbn = '7'
      else if ( i .eq. 8 ) then
        i4_to_isbn = '8'
      else if ( i .eq. 9 ) then
        i4_to_isbn = '9'
      else if ( i .eq. 10 ) then
        i4_to_isbn = 'X'
      else
        i4_to_isbn = '?'
      end if

      return
      end
      function i4_uniform ( a, b, seed )

c*********************************************************************72
c
cc I4_UNIFORM returns a scaled pseudorandom I4.
c
c  Discussion:
c
c    An I4 is an integer value.
c
c    The pseudorandom number should be uniformly distributed
c    between A and B.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    12 November 2006
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Paul Bratley, Bennett Fox, Linus Schrage,
c    A Guide to Simulation,
c    Springer Verlag, pages 201-202, 1983.
c
c    Pierre L'Ecuyer,
c    Random Number Generation,
c    in Handbook of Simulation,
c    edited by Jerry Banks,
c    Wiley Interscience, page 95, 1998.
c
c    Bennett Fox,
c    Algorithm 647:
c    Implementation and Relative Efficiency of Quasirandom
c    Sequence Generators,
c    ACM Transactions on Mathematical Software,
c    Volume 12, Number 4, pages 362-376, 1986.
c
c    Peter Lewis, Allen Goodman, James Miller
c    A Pseudo-Random Number Generator for the System/360,
c    IBM Systems Journal,
c    Volume 8, pages 136-143, 1969.
c
c  Parameters:
c
c    Input, integer A, B, the limits of the interval.
c
c    Input/output, integer SEED, the "seed" value, which should NOT be 0.
c    On output, SEED has been updated.
c
c    Output, integer I4_UNIFORM, a number between A and B.
c
      implicit none

      integer a
      integer b
      integer i4_uniform
      integer k
      real r
      integer seed
      integer value

      if ( seed .eq. 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'I4_UNIFORM - Fatal error!'
        write ( *, '(a)' ) '  Input value of SEED = 0.'
        stop
      end if

      k = seed / 127773

      seed = 16807 * ( seed - k * 127773 ) - k * 2836

      if ( seed .lt. 0 ) then
        seed = seed + 2147483647
      end if

      r = real ( seed ) * 4.656612875E-10
c
c  Scale R to lie between A-0.5 and B+0.5.
c
      r = ( 1.0E+00 - r ) * ( real ( min ( a, b ) ) - 0.5E+00 )
     &  +             r   * ( real ( max ( a, b ) ) + 0.5E+00 )
c
c  Use rounding to convert R to an integer between A and B.
c
      value = nint ( r )

      value = max ( value, min ( a, b ) )
      value = min ( value, max ( a, b ) )

      i4_uniform = value

      return
      end
      subroutine i4vec_indicator ( n, a )

c*********************************************************************72
c
cc I4VEC_INDICATOR sets an I4VEC to the indicator vector.
c
c  Discussion:
c
c    An I4VEC is a vector of I4's.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    02 January 2007
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the number of elements of A.
c
c    Output, integer A(N), the array to be initialized.
c
      implicit none

      integer n

      integer a(n)
      integer i

      do i = 1, n
        a(i) = i
      end do

      return
      end
      subroutine i4vec_min ( n, a, amin )

c*********************************************************************72
c
cc I4VEC_MIN computes the minimum element of an I4VEC.
c
c  Discussion:
c
c    An I4VEC is a vector of I4's.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    19 January 2007
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the number of entries in the array.
c
c    Input, integer A(N), the array.
c
c    Output, integer AMIN, the value of the smallest entry.
c
      implicit none

      integer n

      integer a(n)
      integer amin
      integer i

      amin = a(1)

      do i = 2, n
        amin = min ( amin, a(i) )
      end do

      return
      end
      subroutine i4vec_permute ( n, p, a )

c*********************************************************************72
c
cc I4VEC_PERMUTE permutes an I4VEC in place.
c
c  Discussion:
c
c    An I4VEC is a vector of I4's.
c
c    This routine permutes an array of integer "objects", but the same
c    logic can be used to permute an array of objects of any arithmetic
c    type, or an array of objects of any complexity.  The only temporary
c    storage required is enough to store a single object.  The number
c    of data movements made is N + the number of cycles of order 2 or more,
c    which is never more than N + N/2.
c
c  Example:
c
c    Input:
c
c      N = 5
c      P = (   2,   4,   5,   1,   3 )
c      A = (   1,   2,   3,   4,   5 )
c
c    Output:
c
c      A    = (   2,   4,   5,   1,   3 ).
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    18 July 2010
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the number of objects.
c
c    Input, integer P(N), the permutation.  P(I) = J means
c    that the I-th element of the output array should be the J-th
c    element of the input array.  
c
c    Input/output, integer A(N), the array to be permuted.
c
      implicit none

      integer n

      integer a(n)
      integer a_temp
      integer base
      parameter ( base = 1 )
      integer i
      integer ierror
      integer iget
      integer iput
      integer istart
      integer p(n)

      call perm_check ( n, p, base, ierror )

      if ( ierror .ne. 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'I4VEC_PERMUTE - Fatal error!'
        write ( *, '(a)' ) '  PERM_CHECK rejects this permutation.'
        stop
      end if
c
c  Search for the next element of the permutation that has not been used.
c
      do istart = 1, n

        if ( p(istart) .lt. 0 ) then

          go to 20

        else if ( p(istart) .eq. istart ) then

          p(istart) = - p(istart)
          go to 20

        else

          a_temp = a(istart)
          iget = istart
c
c  Copy the new value into the vacated entry.
c
10        continue

            iput = iget
            iget = p(iget)

            p(iput) = - p(iput)

            if ( iget .lt. 1 .or. n .lt. iget ) then
              write ( *, '(a)' ) ' '
              write ( *, '(a)' ) 'I4VEC_PERMUTE - Fatal error!'
              write ( *, '(a)' ) '  An index is out of range.'
              write ( *, '(a,i8,a,i8)' ) '  P(', iput, ') = ', iget
              stop
            end if

            if ( iget .eq. istart ) then
              a(iput) = a_temp
              go to 20
            end if

            a(iput) = a(iget)

          go to 10

        end if

20      continue

      end do
c
c  Restore the signs of the entries.
c
      do i = 1, n
        p(1:n) = - p(1:n)
      end do

      return
      end
      subroutine i4vec_print ( n, a, title )

c*********************************************************************72
c
cc I4VEC_PRINT prints an I4VEC.
c
c  Discussion:
c
c    An I4VEC is a vector of integer values.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    27 November 2006
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the number of components of the vector.
c
c    Input, integer A(N), the vector to be printed.
c
c    Input, character*(*) TITLE, a title.
c
      implicit none

      integer n

      integer a(n)
      integer i
      character*(*) title

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) trim ( title )
      write ( *, '(a)' ) ' '
      do i = 1, n
        write ( *, '(2x,i8,a,1x,i12)' ) i, ':', a(i)
      end do

      return
      end
      function i4vec_sum ( n, a )

c*********************************************************************72
c
cc I4VEC_SUM returns the sum of the entries of an I4VEC.
c
c  Discussion:
c
c    An I4VEC is a vector of I4's.
c
c    In FORTRAN90, this facility is offered by the built in
c    SUM function:
c
c      I4VEC_SUM ( N, A ) = SUM ( A(1:N) )
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
c    Input, integer N, the number of entries in the array.
c
c    Input, integer A(N), the array.
c
c    Output, integer I4VEC_SUM, the sum of the entries.
c
      implicit none

      integer n

      integer a(n)
      integer i
      integer i4vec_sum

      i4vec_sum = 0

      do i = 1, n
        i4vec_sum = i4vec_sum + a(i)
      end do

      return
      end
      subroutine i4vec_uniform ( n, a, b, seed, x )

c*********************************************************************72
c
cc I4VEC_UNIFORM returns a scaled pseudorandom I4VEC.
c
c  Discussion:
c
c    An I4VEC is a vector of I4's.
c
c    The pseudorandom numbers should be uniformly distributed
c    between A and B.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    12 November 2006
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the dimension of the vector.
c
c    Input, integer A, B, the limits of the interval.
c
c    Input/output, integer SEED, the "seed" value, which should NOT be 0.
c    On output, SEED has been updated.
c
c    Output, integer X(N), a vector of numbers between A and B.
c
      implicit none

      integer n

      integer a
      integer b
      integer i
      integer k
      real r
      integer seed
      integer value
      integer x(n)

      if ( seed .eq. 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'I4VEC_UNIFORM - Fatal error!'
        write ( *, '(a)' ) '  Input value of SEED = 0.'
        stop
      end if

      do i = 1, n

        k = seed / 127773

        seed = 16807 * ( seed - k * 127773 ) - k * 2836

        if ( seed .lt. 0 ) then
          seed = seed + 2147483647
        end if

        r = real ( seed ) * 4.656612875E-10
c
c  Scale R to lie between A-0.5 and B+0.5.
c
        r = ( 1.0E+00 - r ) * ( real ( min ( a, b ) ) - 0.5E+00 )
     &    +             r   * ( real ( max ( a, b ) ) + 0.5E+00 )
c
c  Use rounding to convert R to an integer between A and B.
c
        value = nint ( r )

        value = max ( value, min ( a, b ) )
        value = min ( value, max ( a, b ) )

        x(i) = value

      end do

      return
      end
      subroutine ij_next ( i, j, n )

c*********************************************************************72
c
cc IJ_NEXT returns the next matrix index.
c
c  Discussion:
c
c    For N = 3, the sequence of indices returned is:
c
c      (1,1), (1,2), (1,3), (2,1), (2,2), (2,3), (3,1), (3,2), (3,3), (0,0).
c
c    Note that once the value (N,N) is returned, the next value returned
c    will be (0,0).
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    22 July 2010
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input/output, integer I, J.  On input, the current pair of 
c    indices.  On output, the next pair of indices.  If either index is illegal
c    on input, the output value of (I,J) will be (1,1).
c
c    Input, integer N, the maximum value for I and J.
c
      implicit none

      integer i
      integer j
      integer n

      if ( n .lt. 1 ) then
        i = 0
        j = 0
        return
      end if

      if ( i .lt. 1 .or. n .lt. i .or. j .lt. 1 .or. n .lt. j ) then
        i = 1
        j = 1
        return
      end if

      if ( j .lt. n ) then
        j = j + 1
      else if ( i .lt. n ) then
        i = i + 1
        j = 1
      else
        i = 0
        j = 0
      end if

      return
      end
      subroutine ij_next_gt ( i, j, n )

c*****************************************************************************80
c
cc IJ_NEXT_GT returns the next matrix index, with the constraint that I < J.
c
c  Discussion:
c
c    For N = 5, the sequence of indices returned is:
c
c      (1,2), (1,3), (1,4), (1,5), (2,3), (2,4), (2,5), (3,4), (3,5), (4,5).
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    22 July 2010
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input/output, integer I, J.  On input, the current pair of 
c    indices.  On output, the next pair of indices.  If either index is illegal 
c    on input, the output value of (I,J) will be (1,2).
c
c    Input, integer N, the maximum value for I and J.
c    A value of N less than 2 is nonsense.
c
      implicit none

      integer i
      integer j
      integer n

      if ( n .lt. 2 ) then
        i = 0
        j = 0
        return
      end if

      if ( i .lt. 1 .or. n .lt. i .or. j .lt. 1 .or. n .lt. j .or. 
     &  j .le. i ) then
        i = 1
        j = 2
        return
      end if

      if ( j .lt. n ) then
        j = j + 1
      else if ( i .lt. n - 1 ) then
        i = i + 1
        j = i + 1
      else
        i = 0
        j = 0
      end if

      return
      end
      subroutine index_box2_next_2d ( n1, n2, ic, jc, i, j, more )

c*********************************************************************72
c
cc INDEX_BOX2_NEXT_2D produces indices on the surface of a box in 2D.
c
c  Discussion:
c
c    The box has center at (IC,JC), and has half-widths N1 and N2.
c    The indices are exactly those which are between (IC-N1,JC-N2) and
c    (IC+N1,JC+N2) with the property that at least one of I and J
c    is an "extreme" value.
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
c    Input, integer N1, N2, the half-widths of the box, that is, the
c    maximum distance allowed between (IC,JC) and (I,J).
c
c    Input, integer IC, JC, the central cell of the box.
c
c    Input/output, integer I, J.  On input, the previous index set.
c    On output, the next index set.  On the first call, MORE should
c    be set to FALSE, and the input values of I and J are ignored.
c
c    Input/output, logical MORE.
c    On the first call for a given box, the user should set MORE to FALSE.
c    On return, the routine sets MORE to TRUE.
c    When there are no more indices, the routine sets MORE to FALSE.
c
      implicit none

      integer i
      integer ic
      integer j
      integer jc
      logical more
      integer n1
      integer n2

      if ( .not. more ) then
        more = .true.
        i = ic - n1
        j = jc - n2
        return
      end if

      if ( i .eq. ic + n1 .and. j .eq. jc + n2 ) then
        more = .false.
        return
      end if
c
c  Increment J.
c
      j = j + 1
c
c  Check J.
c
      if ( jc + n2 .lt. j ) then
        j = jc - n2
        i = i + 1
      else if ( j .lt. jc + n2 .and. 
     &  ( i .eq. ic - n1 .or. i .eq. ic + n1 ) ) then
        return
      else
        j = jc + n2
        return
      end if

      return
      end
      subroutine index_box2_next_3d ( n1, n2, n3, ic, jc, kc, i, j, k, 
     &  more )

c*********************************************************************72
c
cc INDEX_BOX2_NEXT_3D produces indices on the surface of a box in 3D.
c
c  Discussion:
c
c    The box has a central cell of (IC,JC,KC), with a half widths of
c    (N1,N2,N3).  The indices are exactly those between (IC-N1,JC-N2,KC-N3)
c    and (IC+N1,JC+N2,KC+N3) with the property that at least one of I, J,
c    and K is an "extreme" value.
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
c    Input, integer N1, N2, N3, the "half widths" of the box, that is, the
c    maximum distances from the central cell allowed for I, J and K.
c
c    Input, integer IC, JC, KC, the central cell of the box.
c
c    Input/output, integer I, J, K.  On input, the previous index set.
c    On output, the next index set.  On the first call, MORE should
c    be set to FALSE, and the input values of I, J, and K are ignored.
c
c    Input/output, logical MORE.
c    On the first call for a given box, the user should set MORE to FALSE.
c    On return, the routine sets MORE to TRUE.
c    When there are no more indices, the routine sets MORE to FALSE.
c
      implicit none

      integer i
      integer ic
      integer j
      integer jc
      integer k
      integer kc
      logical more
      integer n1
      integer n2
      integer n3

      if ( .not. more ) then
        more = .true.
        i = ic - n1
        j = jc - n2
        k = kc - n3
        return
      end if

      if ( i .eq. ic + n1 .and. j .eq. jc + n2 .and. 
     &  k .eq. kc + n3 ) then
        more = .false.
        return
      end if
c
c  Increment K.
c
      k = k + 1
c
c  Check K.
c
      if ( kc + n3 .lt. k ) then
        k = kc - n3
        j = j + 1
      else if ( k .lt. kc + n3 .and. 
     &  ( i .eq. ic - n1 .or. i .eq. ic + n1 .or. 
     &    j .eq. jc - n2 .or. j .eq. jc + n2 ) ) then
        return
      else
        k = kc + n3
        return
      end if
c
c  Check J.
c
      if ( jc + n2 .lt. j ) then
        j = jc - n2
        i = i + 1
      else if ( j .lt. jc + n2 .and. 
     &  ( i .eq. ic - n1 .or. i .eq. ic + n1 .or. 
     &    k .eq. kc - n3 .or. k .eq. kc + n3 ) ) then
        return
      else
        j = jc + n2
        return
      end if

      return
      end
      function index1_col ( i_min, i, i_max, index_min )

c*********************************************************************72
c
cc INDEX1_COL indexes a 1D vector by columns.
c
c  Discussion:
c
c    This 1D routine is provided primarily for analogy.
c    Moreover, in 1D there is no difference between row and column indexing.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    03 April 2010
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer I_MIN, I, I_MAX, for the first index,
c    the minimum, the index, and the maximum.
c
c    Input, integer INDEX_MIN, the index of element I_MIN.
c    Typically, this is 0 or 1.
c
c    Output, integer INDEX1_COL, the index of element I.
c
      implicit none

      integer i
      integer i_max
      integer i_min
      integer index_min
      integer index1_col

      index1_col = index_min + ( i - i_min )

      return
      end
      function index1_row ( i_min, i, i_max, index_min )

c*********************************************************************72
c
cc INDEX1_ROW indexes a 1D vector by rows.
c
c  Discussion:
c
c    This 1D routine is provided primarily for analogy.
c    Moreover, in 1D there is no difference between row and column indexing.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    03 April 2010
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer I_MIN, I, I_MAX, for the first index,
c    the minimum, the index, and the maximum.
c
c    Input, integer INDEX_MIN, the index of element I_MIN.
c    Typically, this is 0 or 1.
c
c    Output, integer INDEX1_ROW, the index of element I.
c
      implicit none

      integer i
      integer i_max
      integer i_min
      integer index_min
      integer index1_row

      index1_row = index_min + ( i - i_min )

      return
      end
      function index2_col ( i_min, i, i_max, j_min, j, j_max, 
     &  index_min )

c*********************************************************************72
c
cc INDEX2_COL indexes a 2D array by columns.
c
c  Discussion:
c
c    Entries of the array are indexed starting at entry (I_MIN,J_MIN), 
c    and increasing the row index first.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    03 April 2010
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer I_MIN, I, I_MAX, for row indices,
c    the minimum, the index, and the maximum.
c
c    Input, integer J_MIN, J, J_MAX, for column indices,
c    the minimum, the index, and the maximum.
c
c    Input, integer INDEX_MIN, the index of element (I_MIN,J_MIN).
c    Typically, this is 0 or 1.
c
c    Output, integer INDEX2_COL, the index of element (I,J).
c
      implicit none

      integer i
      integer i_max
      integer i_min
      integer index_min
      integer index2_col
      integer j
      integer j_max
      integer j_min

      index2_col = index_min + ( i - i_min ) 
     &  + ( j - j_min ) * ( i_max + 1 - i_min )

      return
      end
      function index2_row ( i_min, i, i_max, j_min, j, j_max, 
     &  index_min )

c*********************************************************************72
c
cc INDEX2_ROW indexes a 2D array by row.
c
c  Discussion:
c
c    Entries of the array are indexed starting at entry (I_MIN,J_MIN), 
c    and increasing the column index first.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    03 April 2010
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer I_MIN, I, I_MAX, for row indices,
c    the minimum, the index, and the maximum.
c
c    Input, integer J_MIN, J, J_MAX, for column indices,
c    the minimum, the index, and the maximum.
c
c    Input, integer INDEX_MIN, the index of element (I_MIN,J_MIN).
c    Typically, this is 0 or 1.
c
c    Output, integer INDEX2_ROW, the index of element (I,J).
c
      implicit none

      integer i
      integer i_max
      integer i_min
      integer index_min
      integer index2_row
      integer j
      integer j_max
      integer j_min

      index2_row = index_min + ( j - j_min ) 
     &  + ( i - i_min ) * ( j_max + 1 - j_min )

      return
      end
      function index3_col ( i_min, i, i_max, j_min, j, j_max, 
     &  k_min, k, k_max, index_min )

c*********************************************************************72
c
cc INDEX3_COL indexes a 3D array by columns.
c
c  Discussion:
c
c    Entries of the array are indexed starting at entry (I_MIN,J_MIN,K_MIN), 
c    and increasing the row index first, then the column index.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    03 April 2010
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer I_MIN, I, I_MAX, for row indices,
c    the minimum, the index, and the maximum.
c
c    Input, integer J_MIN, J, J_MAX, for column indices,
c    the minimum, the index, and the maximum.
c
c    Input, integer K_MIN, K, K_MAX, for plane indices,
c    the minimum, the index, and the maximum.
c
c    Input, integer INDEX_MIN, the index of (I_MIN,J_MIN,K_MIN).
c    Typically, this is 0 or 1.
c
c    Output, integer INDEX3_COL, the index of element (I,J,K).
c
      implicit none

      integer i
      integer i_max
      integer i_min
      integer index_min
      integer index3_col
      integer j
      integer j_max
      integer j_min
      integer k
      integer k_max
      integer k_min

      index3_col = index_min 
     &  + ( i - i_min ) 
     &  + ( j - j_min ) * ( i_max + 1 - i_min ) 
     &  + ( k - k_min ) * ( j_max + 1 - j_min ) * ( i_max + 1 - i_min )

      return
      end
      function index3_row ( i_min, i, i_max, j_min, j, j_max, 
     &  k_min, k, k_max, index_min )

c*********************************************************************72
c
cc INDEX3_ROW indexes a 3D array by rows.
c
c  Discussion:
c
c    When we say "by rows", we really just mean that entries of the array are 
c    indexed starting at entry (I_MIN,J_MIN,K_MIN), and the increasing the LAST
c    index first, then the next-to-the-last, and so on.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    03 April 2010
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer I_MIN, I, I_MAX, for row indices,
c    the minimum, the index, and the maximum.
c
c    Input, integer J_MIN, J, J_MAX, for column indices,
c    the minimum, the index, and the maximum.
c
c    Input, integer K_MIN, K, K_MAX, for plane indices,
c    the minimum, the index, and the maximum.
c
c    Input, integer INDEX_MIN, the index of (I_MIN,J_MIN,K_MIN).
c    Typically, this is 0 or 1.
c
c    Output, integer INDEX3_ROW, the index of element (I,J,K).
c
      implicit none

      integer i
      integer i_max
      integer i_min
      integer index_min
      integer index3_row
      integer j
      integer j_max
      integer j_min
      integer k
      integer k_max
      integer k_min

      index3_row = index_min 
     &  + ( k - k_min ) 
     &  + ( j - j_min ) * ( k_max + 1 - k_min ) 
     &  + ( i - i_min ) * ( j_max + 1 - j_min ) * ( k_max + 1 - k_min )

      return
      end
      function index4_col ( i1_min, i1, i1_max, i2_min, i2, i2_max, 
     &  i3_min, i3, i3_max, i4_min, i4, i4_max, index_min )

c*********************************************************************72
c
cc INDEX4_COL indexes a 4D array by columns.
c
c  Discussion:
c
c    Entries of the array are indexed starting at (I1_MIN,I2_MIN,I3_MIN,I4_MIN), 
c    and increasing the initial index first, then the second, third and so on.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    03 April 2010
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer I1_MIN, I1, I1_MAX, for index 1,
c    the minimum, the index, and the maximum.
c
c    Input, integer I2_MIN, I2, I2_MAX, for index 2,
c    the minimum, the index, and the maximum.
c
c    Input, integer I3_MIN, I3, I3_MAX, for index 3,
c    the minimum, the index, and the maximum.
c
c    Input, integer I4_MIN, I4, I4_MAX, for index 4,
c    the minimum, the index, and the maximum.
c
c    Input, integer INDEX_MIN, the index of 
c    (I1_MIN,I2_MIN,I3_MIN,I4_MIN).  Typically, this is 0 or 1.
c
c    Output, integer INDEX4_COL, the index of element (I1,I2,I3,I4).
c
      implicit none

      integer i1
      integer i1_max
      integer i1_min
      integer i2
      integer i2_max
      integer i2_min
      integer i3
      integer i3_max
      integer i3_min
      integer i4
      integer i4_max
      integer i4_min
      integer index_min
      integer index4_col

      index4_col = index_min 
     &  + ( i1 - i1_min ) 
     &  + ( i2 - i2_min ) * ( i1_max + 1 - i1_min ) 
     &  + ( i3 - i3_min ) * ( i2_max + 1 - i2_min ) 
     &  * ( i1_max + 1 - i1_min ) 
     &  + ( i4 - i4_min ) * ( i3_max + 1 - i3_min ) 
     &  * ( i2_max + 1 - i2_min ) * ( i1_max + 1 - i1_min )

      return
      end
      function index4_row ( i1_min, i1, i1_max, i2_min, i2, i2_max, 
     &  i3_min, i3, i3_max, i4_min, i4, i4_max, index_min )

c*********************************************************************72
c
cc INDEX4_ROW indexes a 4D array by rows.
c
c  Discussion:
c
c    Entries of the array are indexed starting at (I1_MIN,I2_MIN,I3_MIN,I4_MIN), 
c    and increasing the last index, then the next to last, and so on.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    03 April 2010
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer I1_MIN, I1, I1_MAX, for index 1,
c    the minimum, the index, and the maximum.
c
c    Input, integer I2_MIN, I2, I2_MAX, for index 2,
c    the minimum, the index, and the maximum.
c
c    Input, integer I3_MIN, I3, I3_MAX, for index 3,
c    the minimum, the index, and the maximum.
c
c    Input, integer I4_MIN, I4, I4_MAX, for index 4,
c    the minimum, the index, and the maximum.
c
c    Input, integer INDEX_MIN, the index of element 
c    (I1_MIN,I2_MIN,I3_MIN,I4_MIN).  Typically, this is 0 or 1.
c
c    Output, integer INDEX4_ROW, the index of element (I1,I2,I3,I4).
c
      implicit none

      integer i1
      integer i1_max
      integer i1_min
      integer i2
      integer i2_max
      integer i2_min
      integer i3
      integer i3_max
      integer i3_min
      integer i4
      integer i4_max
      integer i4_min
      integer index_min
      integer index4_row

      index4_row = index_min 
     &  + ( i4 - i4_min ) 
     &  + ( i3 - i3_min ) * ( i4_max + 1 - i4_min ) 
     &  + ( i2 - i2_min ) * ( i3_max + 1 - i3_min ) 
     &  * ( i4_max + 1 - i4_min ) 
     &  + ( i1 - i1_min ) * ( i2_max + 1 - i2_min ) 
     &  * ( i3_max + 1 - i3_min ) * ( i4_max + 1 - i4_min )

      return
      end
      function indexn_col ( n, i_min, i, i_max, index_min )

c*********************************************************************72
c
cc INDEXN_COL indexes an ND array by columns.
c
c  Discussion:
c
c    Entries of the array are indexed starting at entry 
c      ( I_MIN(1), I_MIN(2),...,I_MIN(N) ), 
c    and increasing the first index up to I_MAX(1), 
c    then the second and so on.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    03 April 2010
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the number of indices.
c
c    Input, integer I_MIN(N), the minimum indices.
c
c    Input, integer I(N), the indices.
c
c    Input, integer I_MAX(N), for maximum indices.
c
c    Input, integer INDEX_MIN, the index of 
c    ( I_MIN(1), I_MIN(2),...,I_MIN(N) ).  Typically, this is 0 or 1.
c
c    Output, integer INDEXN_COL, the index of element I.
c
      implicit none

      integer n

      integer i(n)
      integer i_max(n)
      integer i_min(n)
      integer index_min
      integer indexn_col
      integer j

      indexn_col = ( i(n) - i_min(n) )

      do j = n - 1, 1, - 1
        indexn_col = indexn_col * ( i_max(j) + 1 - i_min(j) ) 
     &    + ( i(j) - i_min(j) )
      end do
      indexn_col = indexn_col + index_min

      return
      end
      function indexn_row ( n, i_min, i, i_max, index_min )

c*********************************************************************72
c
cc INDEXN_ROW indexes an ND array by rows.
c
c  Discussion:
c
c    Entries of the array are indexed starting at entry 
c      ( I_MIN(1), I_MIN(2),...,I_MIN(N) ), 
c    and increasing the last index up to I_MAX(N), 
c    then the next-to-last and so on.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    03 April 2010
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the number of indices.
c
c    Input, integer I_MIN(N), the minimum indices.
c
c    Input, integer I(N), the indices.
c
c    Input, integer I_MAX(N), for maximum indices.
c
c    Input, integer INDEX_MIN, the index of 
c    ( I_MIN(1), I_MIN(2),...,I_MIN(N) ).  Typically, this is 0 or 1.
c
c    Output, integer INDEXN_ROW, the index of element I.
c
      implicit none

      integer n

      integer i(n)
      integer i_max(n)
      integer i_min(n)
      integer index_min
      integer indexn_row
      integer j

      indexn_row = ( i(1) - i_min(1) )

      do j = 2, n
        indexn_row = indexn_row * ( i_max(j) + 1 - i_min(j) ) 
     &    + ( i(j) - i_min(j) )
      end do
      indexn_row = indexn_row + index_min

      return
      end
      subroutine isbn_check ( isbn, check )

c*********************************************************************72
c
cc ISBN_CHECK checks an ISBN code.
c
c  Discussion:
c
c    ISBN stands for International Standard Book Number.  A unique ISBN
c    is assigned to each new book.  The ISBN includes 10 digits.  There is
c    an initial digit, then a dash, then a set of digits which are a
c    code for the publisher, another digit, and then the check digit:
c
c      initial-publisher-book-check
c
c    The number of digits used for the publisher and book codes can vary,
c    but the check digit is always one digit, and the total number of
c    digits is always 10.
c
c    The check digit is interesting, because it is a way of trying to
c    make sure that an ISBN has not been incorrectly copied.  Specifically,
c    if the ISBN is correct, then its ten digits will satisfy
c
c       10 * A + 9 * B + 8 * C + 7 * D + 6 * E
c      + 5 * F * 4 * G * 3 * H + 2 * I +     J  = 0 mod 11.
c
c    Here, we've taken 'A' to represent the first digit and 'J' the
c    last (which is the check digit).  In order for the code to work,
c    the value of J must be allowed to be anything from 0 to 10.  In
c    cases where J works out to be 10, the special digit 'X' is used.
c    An 'X' digit can only occur in this last check-digit position
c    on an ISBN.
c
c  Example:
c
c    0-8493-9640-9
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    25 July 2010
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Book Industry Study Group,
c    The Evolution in Product Identification:
c    Sunrise 2005 and the ISBN-13,
c    http://www.bisg.org/docs/The_Evolution_in_Product_ID.pdf
c
c  Parameters:
c
c    Input, character * ( * ) ISBN, an ISBN code.
c
c    Output, integer CHECK, the value of the ISBN check sum.
c    If CHECK is zero, the ISBN code is legitimate.
c    If CHECK is -1, then the ISBN code is not legitimate because it does
c    not contain exactly 10 digits.  If CHECK is between 1 and 10, then
c    the ISBN code has the right number of digits, but at least one of
c    the digits is incorrect.
c
      implicit none

      character c
      logical ch_is_digit
      integer check
      integer digit(10)
      integer i
      character * ( * ) isbn
      integer isbn_to_i4
      integer lenc
      integer num_digit
c
c  Determine how many digits have been supplied.
c
      lenc = len_trim ( isbn )

      i = 0
      num_digit = 0

10    continue

        i = i + 1

        if ( lenc .lt. i ) then
          go to 20
        end if

        c = isbn(i:i)

        if ( ch_is_digit ( c ) ) then

          num_digit = num_digit + 1
          digit(num_digit) = isbn_to_i4 ( c )

        else if ( ( num_digit .eq. 9 .and. isbn(i:i) .eq. 'X' ) .or. 
     &            ( num_digit .eq. 9 .and. isbn(i:i) .eq. 'x' ) ) then

          num_digit = num_digit + 1
          digit(num_digit) = isbn_to_i4 ( c )

        end if

        if ( 10 .le. num_digit ) then
          go to 20
        end if

      go to 10

20    continue
c
c  If we didn't get exactly 10 digits, return with an error.
c
      if ( num_digit .ne. 10 ) then
        check = -1
        return
      end if
c
c  Compute the checksum.
c
      check = 0
      do i = 1, 10
        check = check + ( 11 - i ) * digit(i)
      end do

      check = mod ( check, 11 )

      return
      end
      subroutine isbn_fill ( isbn )

c*****************************************************************************80
c
cc ISBN_FILL fills in a missing digit in an ISBN code.
c
c  Example:
c
c    Input:
c
c      0-8493-9?40-9
c
c    Output:
c
c      0-8493-9640-9
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    25 July 2010
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Book Industry Study Group,
c    The Evolution in Product Identification:
c    Sunrise 2005 and the ISBN-13,
c    http://www.bisg.org/docs/The_Evolution_in_Product_ID.pdf
c
c  Parameters:
c
c    Input/output, character*(*) ISBN, a partial ISBN code.  On input,
c    a single digit has been replaced by the character '?', signifying
c    that that digit is missing.  The routine replaces the question
c    mark by the correct digit.
c
      implicit none

      character c
      logical ch_is_digit
      integer check
      integer digit(10)
      integer digit_pos
      integer i
      character i4_to_isbn
      character * ( * ) isbn
      integer isbn_pos
      integer isbn_to_i4
      integer j
      integer k
      integer lenc
      integer num_digit

      lenc = len_trim ( isbn )

      i = 0
      isbn_pos = -1
      digit_pos = -1
      num_digit = 0

10    continue

        i = i + 1

        if ( lenc .lt. i ) then
          go to 20
        end if

        c = isbn(i:i)

        if ( ch_is_digit ( c ) ) then

          num_digit = num_digit + 1
          digit(num_digit) = isbn_to_i4 ( c )

        else if ( ( num_digit .eq. 9 .and. isbn(i:i) .eq. 'X' ) .or. 
     &            ( num_digit .eq. 9 .and. isbn(i:i) .eq. 'x' ) ) then

          num_digit = num_digit + 1
          digit(num_digit) = isbn_to_i4 ( c )

        else if ( c .eq. '?' ) then

          if ( isbn_pos .eq. -1 ) then

            num_digit = num_digit + 1
            digit(num_digit) = 0
            digit_pos = num_digit
            isbn_pos = i

          else
            write ( *, '(a)' ) ' '
            write ( *, '(a)' ) 'ISBN_FILL - Fatal error!'
            write ( *, '(a)' ) '  Only one question mark is allowed.'
            return
          end if

        end if

        if ( 10 .le. num_digit ) then
          go to 20
        end if

      go to 10

20    continue

      if ( num_digit .ne. 10 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'ISBN_FILL - Fatal error!'
        write ( *, '(a)' ) 
     &    '  The input ISBN code did not have 10 digits.'
        return
      end if

      if ( isbn_pos .eq. -1 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'ISBN_FILL - Fatal error!'
        write ( *, '(a)' ) '  A question mark is required.'
        return
      end if

      check = 0
      do i = 1, 10
        check = check + ( 11 - i ) * digit(i)
      end do

      check = mod ( check, 11 )

      if ( check .eq. 0 ) then

        k = 0
c
c  Need to solve the modular equation:
c
c    A * X = B mod C
c
c  Below is a stupid way.  One day I will come back and fix this up.
c
      else

        do i = 1, 10
          j = ( 11 - digit_pos ) * i + check
          if ( mod ( j, 11 ) .eq. 0 ) then
            k = i
          end if
        end do

      end if

      isbn(isbn_pos:isbn_pos) = i4_to_isbn ( k )

      return
      end
      function isbn_to_i4 ( c )

c*********************************************************************72
c
cc ISBN_TO_I4 converts an ISBN character into an integer.
c
c  Discussion:
c
c    The characters '0' through '9' stand for themselves, but
c    the character 'X' or 'x' stands for 10.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    22 July 2010
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Book Industry Study Group,
c    The Evolution in Product Identification:
c    Sunrise 2005 and the ISBN-13,
c    http://www.bisg.org/docs/The_Evolution_in_Product_ID.pdf
c
c  Parameters:
c
c    Input, character C, the ISBN character code to be converted.
c
c    Output, integer ISBN_TO_I4, the numeric value of the character
c    code, between 0 and 10.  This value is returned as -1 if C is
c    not a valid character code.
c
      implicit none

      character c
      integer isbn_to_i4

           if ( c .eq. '0' ) then
        isbn_to_i4 = 0
      else if ( c .eq. '1' ) then
        isbn_to_i4 = 1
      else if ( c .eq. '2' ) then
        isbn_to_i4 = 2
      else if ( c .eq. '3' ) then
        isbn_to_i4 = 3
      else if ( c .eq. '4' ) then
        isbn_to_i4 = 4
      else if ( c .eq. '5' ) then
        isbn_to_i4 = 5
      else if ( c .eq. '6' ) then
        isbn_to_i4 = 6
      else if ( c .eq. '7' ) then
        isbn_to_i4 = 7
      else if ( c .eq. '8' ) then
        isbn_to_i4 = 8
      else if ( c .eq. '9' ) then
        isbn_to_i4 = 9
      else if ( c .eq. 'X' .or. c .eq. 'x' ) then
        isbn_to_i4 = 10
      else
        isbn_to_i4 = -1
      end if

      return
      end
      function iset2_compare ( x1, y1, x2, y2 )

!*********************************************************************72
!
!! ISET2_COMPARE compares two I2 sets.
!
!  Discussion:
!
!    The I2 set (X1,Y1) .lt. (X2,Y2) if
!
!      min ( X1, Y1 ) .lt. min ( X2, Y2 ) or
!      min ( X1, Y1 ) = min ( X2, Y2 ) and max ( X1, Y1 ) .lt. max ( X2, Y2 )
!
!    The I2 set (X1,Y1) = (X2,Y2) if
!
!      min ( X1, Y1 ) = min ( X2, Y2 ) and max ( X1, Y1 ) = max ( X2, Y2 )
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    25 July 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer X1, Y1, the first I2 set.
!
!    Input, integer X2, Y2, the second I2 set.
!
!    Output, integer ISET2_COMPARE: 
!    -1, (X1,Y1) .lt. (X2,Y2);
!     0, (X1,Y1) = (X2,Y2);
!    +1, (X1,Y1) > (X2,Y2).
!
      implicit none

      integer a1
      integer a2
      integer b1
      integer b2
      integer iset2_compare
      integer value
      integer x1
      integer x2
      integer y1
      integer y2

      a1 = min ( x1, y1 )
      b1 = max ( x1, y1 )

      a2 = min ( x2, y2 )
      b2 = max ( x2, y2 )

      if ( a1 .lt. a2 ) then
        value = -1
      else if ( a2 .lt. a1 ) then
        value = +1
      else if ( b1 .lt. b2 ) then
        value = -1
      else if ( b2 .lt. b1 ) then
        value = +1
      else
        value = 0
      end if

      iset2_compare = value

      return
      end
      subroutine iset2_index_insert_unique ( n_max, n, x, y, indx, 
     &  xval, yval, ival, ierror )

!*********************************************************************72
!
!! ISET2_INDEX_INSERT_UNIQUE inserts unique values into an indexed sorted list.
!
!  Discussion:
!
!    If the input value does not occur in the list, then N, X, Y and INDX
!    are updated.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    25 July 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N_MAX, the maximum size of the list.
!
!    Input/output, integer N, the size of the list.
!
!    Input/output, integer X(N), Y(N), the list of I2 sets.
!
!    Input/output, integer INDX(N), the sort index of the list.
!
!    Input, integer XVAL, YVAL, the value to be inserted if it is
!    not already in the list.
!
!    Output, integer IVAL, the index in INDX corresponding to the
!    value XVAL, YVAL.
!
!    Output, integer IERROR, 0 for no error, 1 if an error 
!    occurred.
!
      implicit none

      integer n_max

      integer equal
      integer i
      integer ierror
      integer indx(n_max)
      integer ival
      integer less
      integer more
      integer n
      integer x(n_max)
      integer xval
      integer y(n_max)
      integer yval

      ierror = 0

      if ( n .le. 0 ) then

        if ( n_max .le. 0 ) then
          ierror = 1
          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) 'ISET2_INDEX_INSERT_UNIQUE - Fatal error!'
          write ( *, '(a)' ) '  Not enough space to store new data.'
          return
        end if

        n = 1
        x(1) = min ( xval, yval )
        y(1) = max ( xval, yval )
        indx(1) = 1
        ival = 1
        return

      end if
!
!  Does ( XVAL, YVAL ) already occur in the list?
!
      call iset2_index_search ( n_max, n, x, y, indx, xval, yval, 
     &  less, equal, more )

      if ( equal .eq. 0 ) then

        if ( n_max .le. n ) then
          ierror = 1
          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) 'ISET2_INDEX_INSERT_UNIQUE - Fatal error!'
          write ( *, '(a)' ) '  Not enough space to store new data.'
          return
        end if

        x(n+1) = min ( xval, yval )
        y(n+1) = max ( xval, yval )
        ival = more
        do i = n, more, -1
          indx(i+1) = indx(i)
        end do
        indx(more) = n + 1
        n = n + 1

      else

        ival = equal

      end if

      return
      end
      subroutine iset2_index_search ( n_max, n, x, y, indx, xval, yval, 
     &  less, equal, more )

!*********************************************************************72
!
!! ISET2_INDEX_SEARCH searches for an I2 set value in an indexed sorted list.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    25 July 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N_MAX, the maximum size of the list.
!
!    Input, integer N, the size of the current list.
!
!    Input, integer X(N), Y(N), the list.
!
!    Input, integer INDX(N), the sort index of the list.
!
!    Input, integer XVAL, YVAL, the value to be sought.
!
!    Output, integer LESS, EQUAL, MORE, the indexes in INDX of the
!    list entries that are just less than, equal to, and just greater
!    than the test value.  If the test value does not occur in the list,
!    then EQUAL is zero.  If the test value is the minimum entry of the
!    list, then LESS is 0.  If the test value is the greatest entry of
!    the list, then MORE is N+1.
!
      implicit none

      integer n_max

      integer compare
      integer equal
      integer hi
      integer indx(n_max)
      integer less
      integer lo
      integer mid
      integer more
      integer n
      integer iset2_compare
      integer x(n_max)
      integer xhi
      integer xlo
      integer xmid
      integer xval
      integer y(n_max)
      integer yhi
      integer ylo
      integer ymid
      integer yval

      if ( n .le. 0 ) then
        less = 0
        equal = 0
        more = 0
        return
      end if

      lo = 1
      hi = n

      xlo = x(indx(lo))
      ylo = y(indx(lo))

      xhi = x(indx(hi))
      yhi = y(indx(hi))

      compare = iset2_compare ( xval, yval, xlo, ylo )

      if ( compare .eq. -1 ) then
        less = 0
        equal = 0
        more = 1
        return
      else if ( compare .eq. 0 ) then
        less = 0
        equal = 1
        more = 2
        return
      end if

      compare = iset2_compare ( xval, yval, xhi, yhi )

      if ( compare .eq. +1 ) then
        less = n
        equal = 0
        more = n + 1
        return
      else if ( compare .eq. 0 ) then
        less = n - 1
        equal = n
        more = n + 1
        return
      end if

10    continue

        if ( lo + 1 .eq. hi ) then
          less = lo
          equal = 0
          more = hi
          go to 20
        end if

        mid = ( lo + hi ) / 2
        xmid = x(indx(mid))
        ymid = y(indx(mid))

        compare = iset2_compare ( xval, yval, xmid, ymid )

        if ( compare .eq. 0 ) then
          equal = mid
          less = equal - 1
          more = equal + 1
          return
        else if ( compare .eq. -1 ) then
          hi = mid
        else if ( compare .eq. +1 ) then
          lo = mid
        end if

      go to 10

20    continue

      return
      end
      function lcm_12n ( n )

c*********************************************************************72
c
cc LCM_12N computes the least common multiple of the integers 1 through N.
c
c  Example:
c
c    N    LCM_12N
c
c    1          1
c    2          2
c    3          3
c    4         12
c    5         60
c    6         60
c    7        420
c    8        840
c    9       2520
c   10       2520
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    22 July 2010
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the value of N.
c
c    Output, integer LCM_12N, the least common multiple of the 
c    integers 1 to N.
c
      implicit none

      integer i
      integer imult
      integer j
      integer lcm_12n
      integer n

      lcm_12n = 1

      do i = 2, n

        imult = i

        do j = 1, i - 1

          if ( mod ( imult, ( i - j ) ) .eq. 0 ) then
            imult = imult / ( i - j )
          end if

        end do

        lcm_12n = lcm_12n * imult

      end do

      return
      end
      subroutine lmat_print ( m, n, a, title )

c*********************************************************************72
c
cc LMAT_PRINT prints an LMAT.
c
c  Discussion:
c
c    An LMAT is an array of L values.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    03 November 2011
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer M, the number of rows in A.
c
c    Input, integer N, the number of columns in A.
c
c    Input, logical A(M,N), the matrix.
c
c    Input, character * ( * ) TITLE, a title.
c
      implicit none

      integer m
      integer n

      logical a(m,n)
      character * ( * ) title

      call lmat_print_some ( m, n, a, 1, 1, m, n, title )

      return
      end
      subroutine lmat_print_some ( m, n, a, ilo, jlo, ihi, jhi, title )

c*********************************************************************72
c
cc LMAT_PRINT_SOME prints some of an LMAT.
c
c  Discussion:
c
c    An LMAT is an array of L values.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    03 November 2011
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer M, N, the number of rows and columns.
c
c    Input, logical A(M,N), an M by N matrix to be printed.
c
c    Input, integer ILO, JLO, the first row and column to print.
c
c    Input, integer IHI, JHI, the last row and column to print.
c
c    Input, character * ( * ) TITLE, a title.
c
      implicit none

      integer incx
      parameter ( incx = 35 )
      integer m
      integer n

      logical a(m,n)
      character * ( 14 ) ctemp(incx)
      integer i
      integer i2hi
      integer i2lo
      integer ihi
      integer ilo
      integer inc
      integer j
      integer j2
      integer j2hi
      integer j2lo
      integer jhi
      integer jlo
      character * ( * ) title

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) trim ( title )

      do j2lo = max ( jlo, 1 ), min ( jhi, n ), incx

        j2hi = j2lo + incx - 1
        j2hi = min ( j2hi, n )
        j2hi = min ( j2hi, jhi )

        inc = j2hi + 1 - j2lo

        write ( *, '(a)' ) ' '

        if ( 100 .le. j2hi ) then
          do j = j2lo, j2hi
            j2 = j + 1 - j2lo
            write ( ctemp(j2), '(1x,i1)' ) j / 100
          end do
          write ( *, '(''      '',35a2)' ) ctemp(1:inc)
        end if

        if ( 10 .le. j2hi ) then
          do j = j2lo, j2hi
            j2 = j + 1 - j2lo
            write ( ctemp(j2), '(1x,i1)' ) mod ( j / 10, 10 )
          end do
          write ( *, '(''      '',35a2)' ) ctemp(1:inc)
        end if

        do j = j2lo, j2hi
          j2 = j + 1 - j2lo
          write ( ctemp(j2), '(1x,i1)' ) mod ( j, 10 )
        end do
        write ( *, '(''  Col '',35a2)' ) ctemp(1:inc)

        write ( *, '(a)' ) '  Row'
        write ( *, '(a)' ) ' '

        i2lo = max ( ilo, 1 )
        i2hi = min ( ihi, m )

        do i = i2lo, i2hi
          write ( *, '(i5,a1,35(1x,l1))' ) i, ':', a(i,j2lo:j2hi)
        end do

      end do

      return
      end
      subroutine lmat_transpose_print ( m, n, a, title )

c*********************************************************************72
c
cc LMAT_TRANSPOSE_PRINT prints an LMAT, transposed.
c
c  Discussion:
c
c    An LMAT is an array of L values.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    03 November 2011
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer M, N, the number of rows and columns.
c
c    Input, logical A(M,N), an M by N matrix to be printed.
c
c    Input, character * ( * ) TITLE, a title.
c
      implicit none

      integer m
      integer n

      logical a(m,n)
      character * ( * ) title

      call lmat_transpose_print_some ( m, n, a, 1, 1, m, n, title )

      return
      end
      subroutine lmat_transpose_print_some ( m, n, a, ilo, jlo, ihi, 
     &  jhi, title )

c*********************************************************************72
c
cc LMAT_TRANSPOSE_PRINT_SOME prints some of an LMAT, transposed.
c
c  Discussion:
c
c    An LMAT is an array of L values.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    03 November 2011
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer M, N, the number of rows and columns.
c
c    Input, logical A(M,N), an M by N matrix to be printed.
c
c    Input, integer ILO, JLO, the first row and column to print.
c
c    Input, integer IHI, JHI, the last row and column to print.
c
c    Input, character * ( * ) TITLE, a title.
c
      implicit none

      integer incx
      parameter ( incx = 35 )
      integer m
      integer n

      logical a(m,n)
      character * ( 14 ) ctemp(incx)
      integer i
      integer i2
      integer i2hi
      integer i2lo
      integer ihi
      integer ilo
      integer inc
      integer j
      integer j2hi
      integer j2lo
      integer jhi
      integer jlo
      character * ( * ) title

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) trim ( title )

      do i2lo = max ( ilo, 1 ), min ( ihi, m ), incx

        i2hi = i2lo + incx - 1
        i2hi = min ( i2hi, m )
        i2hi = min ( i2hi, ihi )

        inc = i2hi + 1 - i2lo

        write ( *, '(a)' ) ' '

        if ( 100 .le. i2hi ) then
          do i = i2lo, i2hi
            i2 = i + 1 - i2lo
            write ( ctemp(i2), '(1x,i1)' ) i / 100
          end do
          write ( *, '(''      '',35a2)' ) ctemp(1:inc)
        end if

        if ( 10 .le. i2hi ) then
          do i = i2lo, i2hi
            i2 = i + 1 - i2lo
            write ( ctemp(i2), '(1x,i1)' ) mod ( i / 10, 10 )
          end do
          write ( *, '(''      '',35a2)' ) ctemp(1:inc)
        end if

        do i = i2lo, i2hi
          i2 = i + 1 - i2lo
          write ( ctemp(i2), '(1x,i1)' ) mod ( i, 10 )
        end do
        write ( *, '(''  Row '',35a2)' ) ctemp(1:inc)

        write ( *, '(a)' ) '  Col'
        write ( *, '(a)' ) ' '

        j2lo = max ( jlo, 1 )
        j2hi = min ( jhi, n )

        do j = j2lo, j2hi
          write ( *, '(i5,a1,35(1x,l1))' ) j, ':', a(i2lo:i2hi,j)
        end do

      end do

      return
      end
      subroutine luhn_check ( digit_num, digit, check_sum )

c*********************************************************************72
c
cc LUHN_CHECK computes the Luhn checksum for a string of digits.
c
c  Discussion:
c
c    To compute the Luhn checksum, begin at the end of the string, and double
c    every other digit.  If a doubled digit is greater than 9, subtract 9.
c    Then sum the digits to get CHECK_SUM.
c
c    If mod ( CHECK_SUM, 10 ) = 0 the digit sequence is accepted.
c
c    The Luhn check sum will detect any single digit error, as well as
c    most errors involving a transposition of adjacent digits; it cannot
c    detect the transposition of the adjacent digits 0 and 9, however.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    20 July 2010
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer DIGIT_NUM, the number of digits.
c
c    Input, integer DIGIT(DIGIT_NUM), the string of digits.
c    Normally, these are true digits, that is, values between 0 and 9.
c
c    Output, integer CHECK_SUM, the check sum, which
c    should be divisible by 10 if the digit sequence is correct.
c
      implicit none

      integer digit_num

      integer check_sum
      integer digit(digit_num)
      integer digit_copy(digit_num)
      integer i
      integer i4vec_sum

      do i = 1, digit_num
        digit_copy(i) = digit(i)
      end do

      do i = digit_num - 1, 1, -2
        digit_copy(i) = 2 * digit_copy(i)
      end do

      do i = 1, digit_num
        if ( 9 .lt. digit_copy(i) ) then
          digit_copy(i) = digit_copy(i) - 9
        end if
      end do

      check_sum = i4vec_sum ( digit_num, digit_copy )

      return
      end
      subroutine lvec_print ( n, a, title )

c*********************************************************************72
c
cc LVEC_PRINT prints an LVEC.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    02 June 2009
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the number of components of the vector.
c
c    Input, logical A(N), the vector to be printed.
c
c    Input, character * ( * ) TITLE, a title.
c
      implicit none

      integer n

      logical a(n)
      character * ( * ) title

      call lvec_print_some ( n, a, 1, n, title )

      return
      end
      subroutine lvec_print_some ( n, a, i_lo, i_hi, title )

c*********************************************************************72
c
cc LVEC_PRINT_SOME prints "some" of an LVEC.
c
c  Discussion:
c
c    An LVEC is a vector of logical values.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    12 July 2007
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the number of entries of the vector.
c
c    Input, logical A(N), the vector to be printed.
c
c    Input, integer I_LO, I_HI, the first and last indices to 
c    print. The routine expects 1 <= I_LO <= I_HI <= N.
c
c    Input, character * ( * ) TITLE, a title.
c
      implicit none

      integer n

      logical a(n)
      integer i
      integer i_hi
      integer i_lo
      character * ( * ) title

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) trim ( title )
      write ( *, '(a)' ) ' '

      do i = max ( i_lo, 1 ), min ( i_hi, n )
        write ( *, '(2x,i8,a,1x,l1)' ) i, ':', a(i)
      end do

      return
      end
      subroutine perm_check ( n, p, base, ierror )

c*********************************************************************72
c
cc PERM_CHECK checks that a vector represents a permutation.
c
c  Discussion:
c
c    The routine verifies that each of the integers from BASE to
c    to BASE+N-1 occurs among the N entries of the permutation.
c
c    Set the input quantity BASE to 0, if P is a 0-based permutation,
c    or to 1 if P is a 1-based permutation.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    03 June 2009
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the number of entries.
c
c    Input, integer P(N), the array to check.
c
c    Input, integer BASE, the index base.
c
c    Output, integer IERROR, error flag.
c    0, the array represents a permutation.
c    nonzero, the array does not represent a permutation.  The smallest
c    missing value is equal to IERROR.
c
      implicit none

      integer n

      integer base
      integer find
      integer ierror
      integer p(n)
      integer seek

      ierror = 0

      do seek = base, base + n - 1

        ierror = 1

        do find = 1, n
          if ( p(find) .eq. seek ) then
            ierror = 0
            go to 10
          end if
        end do

10      continue

        if ( ierror .ne. 0 ) then
          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) 'PERM_CHECK - Fatal error!'
          write ( *, '(a)' ) '  The input array does not represent'
          write ( *, '(a)' ) '  a proper permutation.'
          stop
        end if

      end do

      return
      end
      subroutine perm_cycle ( n, iopt, p, isgn, ncycle )

c*********************************************************************72
c
cc PERM_CYCLE analyzes a permutation.
c
c  Discussion:
c
c    The routine will count cycles, find the sign of a permutation,
c    and tag a permutation.
c
c  Example:
c
c    Input:
c
c      N = 9
c      IOPT = 1
c      P = 2, 3, 9, 6, 7, 8, 5, 4, 1
c
c    Output:
c
c      NCYCLE = 3
c      ISGN = +1
c      P = -2, 3, 9, -6, -7, 8, 5, 4, 1
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    25 January 2007
c
c  Author:
c
c    Original FORTRAN77 version by Albert Nijenhuis, Herbert Wilf.
c    This FORTRAN77 version by John Burkardt.
c
c  Reference:
c
c    Albert Nijenhuis, Herbert Wilf,
c    Combinatorial Algorithms for Computers and Calculators,
c    Second Edition,
c    Academic Press, 1978,
c    ISBN: 0-12-519260-6,
c    LC: QA164.N54.
c
c  Parameters:
c
c    Input, integer N, the number of objects being permuted.
c
c    Input, integer IOPT, requests tagging.
c    0, the permutation will not be tagged.
c    1, the permutation will be tagged.
c
c    Input/output, integer P(N).  On input, P describes a
c    permutation, in the sense that entry I is to be moved to P(I).
c    If IOPT = 0, then P will not be changed by this routine.
c    If IOPT = 1, then on output, P will be "tagged".  That is,
c    one element of every cycle in P will be negated.  In this way,
c    a user can traverse a cycle by starting at any entry I1 of P
c    which is negative, moving to I2 = ABS(P(I1)), then to
c    P(I2), and so on, until returning to I1.
c
c    Output, integer ISGN, the "sign" of the permutation, which is
c    +1 if the permutation is even, -1 if odd.  Every permutation
c    may be produced by a certain number of pairwise switches.
c    If the number of switches is even, the permutation itself is
c    called even.
c
c    Output, integer NCYCLE, the number of cycles in the permutation.
c
      implicit none

      integer n

      integer base
      integer i
      integer i1
      integer i2
      integer ierror
      integer iopt
      integer is
      integer isgn
      integer ncycle
      integer p(n)

      base = 1
      call perm_check ( n, p, base, ierror )

      is = 1
      ncycle = n

      do i = 1, n

        i1 = p(i)

10      continue

        if ( i .lt. i1 ) then
          ncycle = ncycle - 1
          i2 = p(i1)
          p(i1) = -i2
          i1 = i2
          go to 10
        end if

        if ( iopt .ne. 0 ) then
          is = -sign ( 1, p(i) )
        end if

        p(i) = sign ( p(i), is )

      end do

      isgn = 1 - 2 * mod ( n - ncycle, 2 )

      return
      end
      subroutine perm_free ( npart, ipart, nfree, ifree )

c*********************************************************************72
c
cc PERM_FREE reports the unused items in a partial permutation.
c
c  Discussion:
c
c    It is assumed that the N objects being permuted are the integers
c    from 1 to N, and that IPART contains a "partial" permutation, that
c    is, the NPART entries of IPART represent the beginning of a
c    permutation of all N items.
c
c    The routine returns in IFREE the items that have not been used yet.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    25 January 2007
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer NPART, the number of entries in IPART.  NPART may be 0.
c
c    Input, integer IPART(NPART), the partial permutation, which should
c    contain, at most once, some of the integers between 1 and
c    NPART+NFREE.
c
c    Input, integer NFREE, the number of integers that have not been
c    used in IPART.  This is simply N - NPART.  NFREE may be zero.
c
c    Output, integer IFREE(NFREE), the integers between 1 and NPART+NFREE
c    that were not used in IPART.
c
      implicit none

      integer nfree
      integer npart

      integer i
      integer ifree(nfree)
      integer ipart(npart)
      integer j
      integer k
      integer match
      integer n

      n = npart + nfree

      if ( npart .lt. 0 ) then

        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'PERM_FREE - Fatal error!'
        write ( *, '(a)' ) '  NPART .lt. 0.'
        write ( *, '(a,i8)' ) '  NPART = ', npart
        stop

      else if ( npart .eq. 0 ) then

        call i4vec_indicator ( n, ifree )

      else if ( nfree .lt. 0 ) then

        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'PERM_FREE - Fatal error!'
        write ( *, '(a)' ) '  NFREE .lt. 0.'
        write ( *, '(a,i8)' ) '  NFREE = ', nfree
        stop

      else if ( nfree .eq. 0 ) then

        return

      else

        k = 0

        do i = 1, n

          match = 0

          do j = 1, npart
            if ( ipart(j) .eq. i ) then
              match = j
              exit
            end if
          end do

          if ( match .eq. 0 ) then

            k = k + 1

            if ( nfree .lt. k ) then
              write ( *, '(a)' ) ' '
              write ( *, '(a)' ) 'PERM_FREE - Fatal error!'
              write ( *, '(a)' ) '  The partial permutation is illegal.'
              write ( *, '(a)' ) '  It should contain, at most once,'
              write ( *, '(a,i8)' ) 
     &          '  some of the integers between 1 and ', n
              write ( *, '(a)' ) ' '
              write ( *, '(a)' ) '  Our error is that NFREE .lt. K,'
              write ( *, '(a)' ) '  We have TOO MANY missing values.'
              write ( *, '(a,i8)' ) '  Value of NFREE = ', nfree
              write ( *, '(a,i8)' ) '  Value of K =     ', k
              call i4vec_print ( npart, ipart, 
     &          '  Partial permutation:' )
              stop
            end if

            ifree(k) = i

          end if

        end do

      end if

      return
      end
      subroutine perm_inverse ( n, p )

c*********************************************************************72
c
cc PERM_INVERSE inverts a permutation "in place".
c
c  Discussion:
c
c    This algorithm assumes that the entries in the permutation vector are
c    strictly positive.  In particular, the value 0 must not occur.
c
c    When necessary, this function shifts the data temporarily so that
c    this requirement is satisfied.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    03 June 2009
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the number of objects being permuted.
c
c    Input/output, integer P(N), the permutation, in standard index form.
c    On output, P describes the inverse permutation
c
      implicit none

      integer n

      integer base
      integer i
      integer i0
      integer i1
      integer i2
      integer ierror
      integer is
      integer p(n)
      integer p_min

      if ( n .le. 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'PERM_INVERSE - Fatal error!'
        write ( *, '(a,i8)' ) '  Input value of N = ', n
        stop
      end if
c
c  Find the least value, and shift data so it begins at 1.
c
      call i4vec_min ( n, p, p_min )
      base = 1
      do i = 1, n
        p(i) = p(i) - p_min + base
      end do
c
c  Check the permutation.
c
      call perm_check ( n, p, base, ierror )

      if ( ierror .ne. 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'PERM_INVERSE - Fatal error!'
        write ( *, '(a)' ) '  PERM_CHECK rejects this permutation.'
        stop
      end if
c
c  Invert the permutation.
c
      is = 1

      do i = 1, n

        i1 = p(i)

10      continue

        if ( i .lt. i1 ) then
          i2 = p(i1)
          p(i1) = -i2
          i1 = i2
          go to 10
        end if

        is = -sign ( 1, p(i) )
        p(i) = sign ( p(i), is )

      end do

      do i = 1, n

        i1 = -p(i)

        if ( 0 .le. i1 ) then

          i0 = i

20        continue

            i2 = p(i1)
            p(i1) = i0

            if ( i2 .lt. 0 ) then
              go to 30
            end if

            i0 = i1
            i1 = i2

          go to 20

30        continue

        end if

      end do
c
c  Undo the shift.
c
      do i = 1, n
        p(i) = p(i) + p_min - base
      end do

      return
      end
      subroutine perm_next ( n, p, more, even )

c*********************************************************************72
c
cc PERM_NEXT computes all of the permutations of N objects, one at a time.
c
c  Discussion:
c
c    The routine is initialized by calling with MORE = TRUE, in which case
c    it returns the identity permutation.
c
c    If the routine is called with MORE = FALSE, then the successor of the
c    input permutation is computed.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    25 January 2007
c
c  Author:
c
c    Original FORTRAN77 version by Albert Nijenhuis, Herbert Wilf.
c    This FORTRAN77 version by John Burkardt.
c
c  Reference:
c
c    Albert Nijenhuis, Herbert Wilf,
c    Combinatorial Algorithms for Computers and Calculators,
c    Second Edition,
c    Academic Press, 1978,
c    ISBN: 0-12-519260-6,
c    LC: QA164.N54.
c
c  Parameters:
c
c    Input, integer N, the number of objects being permuted.
c
c    Input/output, integer P(N), the permutation, in standard index form.
c    On the first call, the input value is unimportant.
c    On subsequent calls, the input value should be the same
c    as the output value from the previous call.  In other words, the
c    user should just leave P alone.
c    On output, contains the "next" permutation.
c
c    Input/output, logical MORE.
c    Set MORE = FALSE before the first call.
c    MORE will be reset to TRUE and a permutation will be returned.
c    Each new call produces a new permutation until
c    MORE is returned FALSE.
c
c    Input/output, logical EVEN.
c    The input value of EVEN should simply be its output value from the
c    previous call; (the input value on the first call doesn't matter.)
c    On output, EVEN is TRUE if the output permutation is even, that is,
c    involves an even number of transpositions.
c
      implicit none

      integer n

      logical even
      integer i
      integer i1
      integer ia
      integer id
      integer is
      integer j
      integer l
      integer m
      logical more
      integer p(n)

      if ( .not. more ) then

        call i4vec_indicator ( n, p )
        more = .true.
        even = .true.

        if ( n .eq. 1 ) then
          more = .false.
          return
        end if

        if ( p(n) .ne. 1 .or. p(1) .ne. 2 + mod ( n, 2 ) ) then
          return
        end if

        do i = 1, n-3
          if ( p(i+1) .ne. p(i)+1 ) then
            return
          end if
        end do

        more = .false.

      else

        if ( n .eq. 1 ) then
          p(1) = 0
          more = .false.
          return
        end if

        if ( even ) then

          ia = p(1)
          p(1) = p(2)
          p(2) = ia
          even = .false.

          if ( p(n) .ne. 1 .or. p(1) .ne. 2 + mod ( n, 2 ) ) then
            return
          end if

          do i = 1, n-3
            if ( p(i+1) .ne. p(i)+1 ) then
              return
            end if
          end do

          more = .false.
          return

        else

          more = .false.

          is = 0

          do i1 = 2, n

            ia = p(i1)
            i = i1 - 1
            id = 0

            do j = 1, i
              if ( ia .lt. p(j) ) then
                id = id + 1
              end if
            end do

            is = id + is
            if ( id .ne. i * mod ( is, 2 ) ) then
              more = .true.
              exit
            end if

          end do

          if ( .not. more ) then
            p(1) = 0
            return
          end if

        end if

        m = mod ( is+1, 2 ) * ( n + 1 )

        do j = 1, i

          if ( sign ( 1, p(j)-ia ) .ne. sign ( 1, p(j)-m ) ) then
            m = p(j)
            l = j
          end if

        end do

        p(l) = ia
        p(i1) = m
        even = .true.

      end if

      return
      end
      subroutine perm_print ( n, p, title )

c*********************************************************************72
c
cc PERM_PRINT prints a permutation.
c
c  Example:
c
c    Input:
c
c      P = 7 2 4 1 5 3 6
c
c    Printed output:
c
c      "This is the permutation:"
c
c      1 2 3 4 5 6 7
c      7 2 4 1 5 3 6
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    25 January 2007
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the number of objects permuted.
c
c    Input, integer P(N), the permutation, in standard index form.
c
c    Input, character * ( * ) TITLE, an optional title.
c    If no title is supplied, then only the permutation is printed.
c
      implicit none

      integer n

      integer i
      integer ihi
      integer ilo
      integer inc
      parameter ( inc = 20 )
      integer p(n)
      character * ( * ) title
      integer title_length

      title_length = len_trim ( title )

      if ( 0 .lt. title_length ) then

        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) trim ( title )

        do ilo = 1, n, inc
          ihi = min ( n, ilo + inc - 1 )
          write ( *, '(a)' ) ' '
          write ( *, '(2x,20i4)' ) ( i, i = ilo, ihi )
          write ( *, '(2x,20i4)' ) ( p(i), i = ilo, ihi )
        end do

      else

        do ilo = 1, n, inc
          ihi = min ( n, ilo + inc - 1 )
          write ( *, '(2x,20i4)' ) ( p(i), i = ilo, ihi )
        end do

      end if

      return
      end
      subroutine perm_uniform ( n, seed, p )

c*********************************************************************72
c
cc PERM_UNIFORM selects a random permutation of N objects.
c
c  Discussion:
c
c    The routine assumes the objects are labeled 1, 2, ... N.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    25 January 2007
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Albert Nijenhuis, Herbert Wilf,
c    Combinatorial Algorithms for Computers and Calculators,
c    Second Edition,
c    Academic Press, 1978,
c    ISBN: 0-12-519260-6,
c    LC: QA164.N54.
c
c  Parameters:
c
c    Input, integer N, the number of objects to be permuted.
c
c    Input/output, integer SEED, a seed for the random number generator.
c
c    Output, integer P(N), a permutation of ( 1, 2, ..., N ), in standard 
c    index form.
c
      implicit none

      integer n

      integer i
      integer i4_uniform
      integer j
      integer p(n)
      integer pk
      integer seed

      do i = 1, n
        p(i) = i
      end do

      do i = 1, n
        j = i4_uniform ( i, n, seed )
        pk = p(i)
        p(i) = p(j)
        p(j) = pk
      end do

      return
      end
      function pounds_to_kilograms ( lb )

c*********************************************************************72
c
cc POUNDS_TO_KILOGRAMS converts a measurement in pounds to kilograms.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    20 May 2008
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, double precision LB, the weight in pounds.
c
c    Output, double precision POUNDS_TO_KILOGRAMS, the corresponding 
c    weight in kilograms.
c
      implicit none

      double precision lb
      double precision pounds_to_kilograms

      pounds_to_kilograms = 0.4535924D+00 * lb

      return
      end
      function prime ( n )

c*********************************************************************72
c
cc PRIME returns any of the first PRIME_MAX prime numbers.
c
c  Discussion:
c
c    PRIME_MAX is 1600, and the largest prime stored is 13499.
c
c    Thanks to Bart Vandewoestyne for pointing out a typo, 18 February 2005.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    21 January 2007
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Milton Abramowitz, Irene Stegun,
c    Handbook of Mathematical Functions,
c    National Bureau of Standards, 1964,
c    ISBN: 0-486-61272-4,
c    LC: QA47.A34.
c
c    Daniel Zwillinger,
c    CRC Standard Mathematical Tables and Formulae,
c    30th Edition,
c    CRC Press, 1996, pages 95-98.
c
c  Parameters:
c
c    Input, integer N, the index of the desired prime number.
c    In general, is should be true that 0 <= N <= PRIME_MAX.
c    N = -1 returns PRIME_MAX, the index of the largest prime available.
c    N = 0 is legal, returning PRIME = 1.
c
c    Output, integer PRIME, the N-th prime.  If N is out of range,
c    PRIME is returned as -1.
c
      implicit none

      integer prime_max
      parameter ( prime_max = 1600 )

      integer i
      integer n
      integer npvec(prime_max)
      integer prime

      save npvec

      data ( npvec(i), i = 1, 100 ) /
     &      2,    3,    5,    7,   11,   13,   17,   19,   23,   29,
     &     31,   37,   41,   43,   47,   53,   59,   61,   67,   71,
     &     73,   79,   83,   89,   97,  101,  103,  107,  109,  113,
     &    127,  131,  137,  139,  149,  151,  157,  163,  167,  173,
     &    179,  181,  191,  193,  197,  199,  211,  223,  227,  229,
     &    233,  239,  241,  251,  257,  263,  269,  271,  277,  281,
     &    283,  293,  307,  311,  313,  317,  331,  337,  347,  349,
     &    353,  359,  367,  373,  379,  383,  389,  397,  401,  409,
     &    419,  421,  431,  433,  439,  443,  449,  457,  461,  463,
     &    467,  479,  487,  491,  499,  503,  509,  521,  523,  541 /

      data ( npvec(i), i = 101, 200 ) /
     &    547,  557,  563,  569,  571,  577,  587,  593,  599,  601,
     &    607,  613,  617,  619,  631,  641,  643,  647,  653,  659,
     &    661,  673,  677,  683,  691,  701,  709,  719,  727,  733,
     &    739,  743,  751,  757,  761,  769,  773,  787,  797,  809,
     &    811,  821,  823,  827,  829,  839,  853,  857,  859,  863,
     &    877,  881,  883,  887,  907,  911,  919,  929,  937,  941,
     &    947,  953,  967,  971,  977,  983,  991,  997, 1009, 1013,
     &   1019, 1021, 1031, 1033, 1039, 1049, 1051, 1061, 1063, 1069,
     &   1087, 1091, 1093, 1097, 1103, 1109, 1117, 1123, 1129, 1151,
     &   1153, 1163, 1171, 1181, 1187, 1193, 1201, 1213, 1217, 1223 /

      data ( npvec(i), i = 201, 300 ) /
     &   1229, 1231, 1237, 1249, 1259, 1277, 1279, 1283, 1289, 1291,
     &   1297, 1301, 1303, 1307, 1319, 1321, 1327, 1361, 1367, 1373,
     &   1381, 1399, 1409, 1423, 1427, 1429, 1433, 1439, 1447, 1451,
     &   1453, 1459, 1471, 1481, 1483, 1487, 1489, 1493, 1499, 1511,
     &   1523, 1531, 1543, 1549, 1553, 1559, 1567, 1571, 1579, 1583,
     &   1597, 1601, 1607, 1609, 1613, 1619, 1621, 1627, 1637, 1657,
     &   1663, 1667, 1669, 1693, 1697, 1699, 1709, 1721, 1723, 1733,
     &   1741, 1747, 1753, 1759, 1777, 1783, 1787, 1789, 1801, 1811,
     &   1823, 1831, 1847, 1861, 1867, 1871, 1873, 1877, 1879, 1889,
     &   1901, 1907, 1913, 1931, 1933, 1949, 1951, 1973, 1979, 1987 /

      data ( npvec(i), i = 301, 400 ) /
     &   1993, 1997, 1999, 2003, 2011, 2017, 2027, 2029, 2039, 2053,
     &   2063, 2069, 2081, 2083, 2087, 2089, 2099, 2111, 2113, 2129,
     &   2131, 2137, 2141, 2143, 2153, 2161, 2179, 2203, 2207, 2213,
     &   2221, 2237, 2239, 2243, 2251, 2267, 2269, 2273, 2281, 2287,
     &   2293, 2297, 2309, 2311, 2333, 2339, 2341, 2347, 2351, 2357,
     &   2371, 2377, 2381, 2383, 2389, 2393, 2399, 2411, 2417, 2423,
     &   2437, 2441, 2447, 2459, 2467, 2473, 2477, 2503, 2521, 2531,
     &   2539, 2543, 2549, 2551, 2557, 2579, 2591, 2593, 2609, 2617,
     &   2621, 2633, 2647, 2657, 2659, 2663, 2671, 2677, 2683, 2687,
     &   2689, 2693, 2699, 2707, 2711, 2713, 2719, 2729, 2731, 2741 /

      data ( npvec(i), i = 401, 500 ) /
     &   2749, 2753, 2767, 2777, 2789, 2791, 2797, 2801, 2803, 2819,
     &   2833, 2837, 2843, 2851, 2857, 2861, 2879, 2887, 2897, 2903,
     &   2909, 2917, 2927, 2939, 2953, 2957, 2963, 2969, 2971, 2999,
     &   3001, 3011, 3019, 3023, 3037, 3041, 3049, 3061, 3067, 3079,
     &   3083, 3089, 3109, 3119, 3121, 3137, 3163, 3167, 3169, 3181,
     &   3187, 3191, 3203, 3209, 3217, 3221, 3229, 3251, 3253, 3257,
     &   3259, 3271, 3299, 3301, 3307, 3313, 3319, 3323, 3329, 3331,
     &   3343, 3347, 3359, 3361, 3371, 3373, 3389, 3391, 3407, 3413,
     &   3433, 3449, 3457, 3461, 3463, 3467, 3469, 3491, 3499, 3511,
     &   3517, 3527, 3529, 3533, 3539, 3541, 3547, 3557, 3559, 3571 /

      data ( npvec(i), i = 501, 600 ) /
     &   3581, 3583, 3593, 3607, 3613, 3617, 3623, 3631, 3637, 3643,
     &   3659, 3671, 3673, 3677, 3691, 3697, 3701, 3709, 3719, 3727,
     &   3733, 3739, 3761, 3767, 3769, 3779, 3793, 3797, 3803, 3821,
     &   3823, 3833, 3847, 3851, 3853, 3863, 3877, 3881, 3889, 3907,
     &   3911, 3917, 3919, 3923, 3929, 3931, 3943, 3947, 3967, 3989,
     &   4001, 4003, 4007, 4013, 4019, 4021, 4027, 4049, 4051, 4057,
     &   4073, 4079, 4091, 4093, 4099, 4111, 4127, 4129, 4133, 4139,
     &   4153, 4157, 4159, 4177, 4201, 4211, 4217, 4219, 4229, 4231,
     &   4241, 4243, 4253, 4259, 4261, 4271, 4273, 4283, 4289, 4297,
     &   4327, 4337, 4339, 4349, 4357, 4363, 4373, 4391, 4397, 4409 /

      data ( npvec(i), i = 601, 700 ) /
     &   4421, 4423, 4441, 4447, 4451, 4457, 4463, 4481, 4483, 4493,
     &   4507, 4513, 4517, 4519, 4523, 4547, 4549, 4561, 4567, 4583,
     &   4591, 4597, 4603, 4621, 4637, 4639, 4643, 4649, 4651, 4657,
     &   4663, 4673, 4679, 4691, 4703, 4721, 4723, 4729, 4733, 4751,
     &   4759, 4783, 4787, 4789, 4793, 4799, 4801, 4813, 4817, 4831,
     &   4861, 4871, 4877, 4889, 4903, 4909, 4919, 4931, 4933, 4937,
     &   4943, 4951, 4957, 4967, 4969, 4973, 4987, 4993, 4999, 5003,
     &   5009, 5011, 5021, 5023, 5039, 5051, 5059, 5077, 5081, 5087,
     &   5099, 5101, 5107, 5113, 5119, 5147, 5153, 5167, 5171, 5179,
     &   5189, 5197, 5209, 5227, 5231, 5233, 5237, 5261, 5273, 5279 /

      data ( npvec(i), i = 701, 800 ) /
     &   5281, 5297, 5303, 5309, 5323, 5333, 5347, 5351, 5381, 5387,
     &   5393, 5399, 5407, 5413, 5417, 5419, 5431, 5437, 5441, 5443,
     &   5449, 5471, 5477, 5479, 5483, 5501, 5503, 5507, 5519, 5521,
     &   5527, 5531, 5557, 5563, 5569, 5573, 5581, 5591, 5623, 5639,
     &   5641, 5647, 5651, 5653, 5657, 5659, 5669, 5683, 5689, 5693,
     &   5701, 5711, 5717, 5737, 5741, 5743, 5749, 5779, 5783, 5791,
     &   5801, 5807, 5813, 5821, 5827, 5839, 5843, 5849, 5851, 5857,
     &   5861, 5867, 5869, 5879, 5881, 5897, 5903, 5923, 5927, 5939,
     &   5953, 5981, 5987, 6007, 6011, 6029, 6037, 6043, 6047, 6053,
     &   6067, 6073, 6079, 6089, 6091, 6101, 6113, 6121, 6131, 6133 /

      data ( npvec(i), i = 801, 900 ) /
     &   6143, 6151, 6163, 6173, 6197, 6199, 6203, 6211, 6217, 6221,
     &   6229, 6247, 6257, 6263, 6269, 6271, 6277, 6287, 6299, 6301,
     &   6311, 6317, 6323, 6329, 6337, 6343, 6353, 6359, 6361, 6367,
     &   6373, 6379, 6389, 6397, 6421, 6427, 6449, 6451, 6469, 6473,
     &   6481, 6491, 6521, 6529, 6547, 6551, 6553, 6563, 6569, 6571,
     &   6577, 6581, 6599, 6607, 6619, 6637, 6653, 6659, 6661, 6673,
     &   6679, 6689, 6691, 6701, 6703, 6709, 6719, 6733, 6737, 6761,
     &   6763, 6779, 6781, 6791, 6793, 6803, 6823, 6827, 6829, 6833,
     &   6841, 6857, 6863, 6869, 6871, 6883, 6899, 6907, 6911, 6917,
     &   6947, 6949, 6959, 6961, 6967, 6971, 6977, 6983, 6991, 6997 /

      data ( npvec(i), i = 901, 1000 ) /
     &   7001, 7013, 7019, 7027, 7039, 7043, 7057, 7069, 7079, 7103,
     &   7109, 7121, 7127, 7129, 7151, 7159, 7177, 7187, 7193, 7207,
     &   7211, 7213, 7219, 7229, 7237, 7243, 7247, 7253, 7283, 7297,
     &   7307, 7309, 7321, 7331, 7333, 7349, 7351, 7369, 7393, 7411,
     &   7417, 7433, 7451, 7457, 7459, 7477, 7481, 7487, 7489, 7499,
     &   7507, 7517, 7523, 7529, 7537, 7541, 7547, 7549, 7559, 7561,
     &   7573, 7577, 7583, 7589, 7591, 7603, 7607, 7621, 7639, 7643,
     &   7649, 7669, 7673, 7681, 7687, 7691, 7699, 7703, 7717, 7723,
     &   7727, 7741, 7753, 7757, 7759, 7789, 7793, 7817, 7823, 7829,
     &   7841, 7853, 7867, 7873, 7877, 7879, 7883, 7901, 7907, 7919 /

      data ( npvec(i), i = 1001, 1100 ) /
     &   7927, 7933, 7937, 7949, 7951, 7963, 7993, 8009, 8011, 8017,
     &   8039, 8053, 8059, 8069, 8081, 8087, 8089, 8093, 8101, 8111,
     &   8117, 8123, 8147, 8161, 8167, 8171, 8179, 8191, 8209, 8219,
     &   8221, 8231, 8233, 8237, 8243, 8263, 8269, 8273, 8287, 8291,
     &   8293, 8297, 8311, 8317, 8329, 8353, 8363, 8369, 8377, 8387,
     &   8389, 8419, 8423, 8429, 8431, 8443, 8447, 8461, 8467, 8501,
     &   8513, 8521, 8527, 8537, 8539, 8543, 8563, 8573, 8581, 8597,
     &   8599, 8609, 8623, 8627, 8629, 8641, 8647, 8663, 8669, 8677,
     &   8681, 8689, 8693, 8699, 8707, 8713, 8719, 8731, 8737, 8741,
     &   8747, 8753, 8761, 8779, 8783, 8803, 8807, 8819, 8821, 8831 /

      data ( npvec(i), i = 1101, 1200 ) /
     &   8837, 8839, 8849, 8861, 8863, 8867, 8887, 8893, 8923, 8929,
     &   8933, 8941, 8951, 8963, 8969, 8971, 8999, 9001, 9007, 9011,
     &   9013, 9029, 9041, 9043, 9049, 9059, 9067, 9091, 9103, 9109,
     &   9127, 9133, 9137, 9151, 9157, 9161, 9173, 9181, 9187, 9199,
     &   9203, 9209, 9221, 9227, 9239, 9241, 9257, 9277, 9281, 9283,
     &   9293, 9311, 9319, 9323, 9337, 9341, 9343, 9349, 9371, 9377,
     &   9391, 9397, 9403, 9413, 9419, 9421, 9431, 9433, 9437, 9439,
     &   9461, 9463, 9467, 9473, 9479, 9491, 9497, 9511, 9521, 9533,
     &   9539, 9547, 9551, 9587, 9601, 9613, 9619, 9623, 9629, 9631,
     &   9643, 9649, 9661, 9677, 9679, 9689, 9697, 9719, 9721, 9733 /

      data ( npvec(i), i = 1201, 1300 ) /
     &   9739, 9743, 9749, 9767, 9769, 9781, 9787, 9791, 9803, 9811,
     &   9817, 9829, 9833, 9839, 9851, 9857, 9859, 9871, 9883, 9887,
     &   9901, 9907, 9923, 9929, 9931, 9941, 9949, 9967, 9973,10007,
     &  10009,10037,10039,10061,10067,10069,10079,10091,10093,10099,
     &  10103,10111,10133,10139,10141,10151,10159,10163,10169,10177,
     &  10181,10193,10211,10223,10243,10247,10253,10259,10267,10271,
     &  10273,10289,10301,10303,10313,10321,10331,10333,10337,10343,
     &  10357,10369,10391,10399,10427,10429,10433,10453,10457,10459,
     &  10463,10477,10487,10499,10501,10513,10529,10531,10559,10567,
     &  10589,10597,10601,10607,10613,10627,10631,10639,10651,10657 /

      data ( npvec(i), i = 1301, 1400 ) /
     &  10663,10667,10687,10691,10709,10711,10723,10729,10733,10739,
     &  10753,10771,10781,10789,10799,10831,10837,10847,10853,10859,
     &  10861,10867,10883,10889,10891,10903,10909,10937,10939,10949,
     &  10957,10973,10979,10987,10993,11003,11027,11047,11057,11059,
     &  11069,11071,11083,11087,11093,11113,11117,11119,11131,11149,
     &  11159,11161,11171,11173,11177,11197,11213,11239,11243,11251,
     &  11257,11261,11273,11279,11287,11299,11311,11317,11321,11329,
     &  11351,11353,11369,11383,11393,11399,11411,11423,11437,11443,
     &  11447,11467,11471,11483,11489,11491,11497,11503,11519,11527,
     &  11549,11551,11579,11587,11593,11597,11617,11621,11633,11657 /

      data ( npvec(i), i = 1401, 1500 ) /
     &  11677,11681,11689,11699,11701,11717,11719,11731,11743,11777,
     &  11779,11783,11789,11801,11807,11813,11821,11827,11831,11833,
     &  11839,11863,11867,11887,11897,11903,11909,11923,11927,11933,
     &  11939,11941,11953,11959,11969,11971,11981,11987,12007,12011,
     &  12037,12041,12043,12049,12071,12073,12097,12101,12107,12109,
     &  12113,12119,12143,12149,12157,12161,12163,12197,12203,12211,
     &  12227,12239,12241,12251,12253,12263,12269,12277,12281,12289,
     &  12301,12323,12329,12343,12347,12373,12377,12379,12391,12401,
     &  12409,12413,12421,12433,12437,12451,12457,12473,12479,12487,
     &  12491,12497,12503,12511,12517,12527,12539,12541,12547,12553 /

      data ( npvec(i), i = 1501, 1600 ) /
     &  12569,12577,12583,12589,12601,12611,12613,12619,12637,12641,
     &  12647,12653,12659,12671,12689,12697,12703,12713,12721,12739,
     &  12743,12757,12763,12781,12791,12799,12809,12821,12823,12829,
     &  12841,12853,12889,12893,12899,12907,12911,12917,12919,12923,
     &  12941,12953,12959,12967,12973,12979,12983,13001,13003,13007,
     &  13009,13033,13037,13043,13049,13063,13093,13099,13103,13109,
     &  13121,13127,13147,13151,13159,13163,13171,13177,13183,13187,
     &  13217,13219,13229,13241,13249,13259,13267,13291,13297,13309,
     &  13313,13327,13331,13337,13339,13367,13381,13397,13399,13411,
     &  13417,13421,13441,13451,13457,13463,13469,13477,13487,13499 /

      if ( n .eq. -1 ) then
        prime = prime_max
      else if ( n .eq. 0 ) then
        prime = 1
      else if ( n .le. prime_max ) then
        prime = npvec(n)
      else
        prime = -1
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'PRIME - Fatal error!'
        write ( *, '(a,i8)' ) '  Illegal prime index N = ', n
        write ( *, '(a,i8)' )
     &    '  N should be between 1 and PRIME_MAX =', prime_max
        stop
      end if

      return
      end
      function prime_ge ( n )

c*********************************************************************72
c
cc PRIME_GE returns the smallest prime greater than or equal to N.
c
c  Example:
c
c    N     PRIME_GE
c
c    -10    2
c      1    2
c      2    2
c      3    3
c      4    5
c      5    5
c      6    7
c      7    7
c      8   11
c      9   11
c     10   11
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    23 August 2008
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the number to be bounded.
c
c    Output, integer PRIME_GE, the smallest prime number that is 
c    greater than or equal to N.  However, if N is larger than the largest
c    prime stored, then PRIME_GE is returned as -1.
c
      implicit none

      integer i_hi
      integer i_lo
      integer i_mid
      integer n
      integer p_hi
      integer p_lo
      integer p_mid
      integer prime
      integer prime_ge

      if ( n .le. 2 ) then
        prime_ge = 2
        return
      end if

      i_lo = 1
      p_lo = prime(i_lo)
      i_hi = prime(-1)
      p_hi = prime(i_hi)

      if ( p_hi .lt. n ) then
        prime_ge = -p_hi
        return
      end if

10    continue

        if ( i_lo + 1 .eq. i_hi ) then
          prime_ge = p_hi
          go to 20
        end if

        i_mid = ( i_lo + i_hi ) / 2
        p_mid = prime(i_mid)

        if ( p_mid .lt. n ) then
          i_lo = i_mid
          p_lo = p_mid
        else if ( n .le. p_mid ) then
          i_hi = i_mid
          p_hi = p_mid
        end if

      go to 10

20    continue

      return
      end
      subroutine primer ( n, iprime )

c*********************************************************************72
c
cc PRIMER computes the prime numbers up to a given limit.
c
c  Discussion:
c
c    PRIMER returns the results of its computations in the vector
c    IPRIME.  IPRIME(I) is -1 if the number I is not prime, and
c    1 if I is prime.
c
c    The algorithm is a simple-minded sieve of Eratosthenes, with
c    no attempt at efficiency.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    23 August 2008
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the dimension of IPRIME, and the maximum
c    value that will be considered.
c
c    Output, integer IPRIME(N), records the results for each 
c    integer.  IPRIME(I) = -1 if I is not prime, and IPRIME(I) = 1 if I is
c    prime.  By convention, IPRIME(1) will be set to -1.
c
      implicit none

      integer n

      integer i
      integer iprime(n)
      integer next
c
c  IPRIME(I) = 0 means we don't know if I is prime.
c
      do i = 1, n
        iprime(i) = 0
      end do
c
c  By convention, 1 is not prime.
c
      iprime(1) = -1
      next = 1
c
c  Examine the integers in order.
c
      do next = 2, n

        if ( iprime(next) .eq. 0 ) then
          iprime(next) = 1
          do i = 2 * next, n, next
            iprime(i) = -1
          end do
        end if

      end do

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
      function r8_log_10 ( x )

c*********************************************************************72
c
cc R8_LOG_10 returns the logarithm base 10 of an R8.
c
c  Discussion:
c
c    value = Log10 ( |X| )
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    16 January 2007
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, double precision X, the number whose base 2 logarithm is desired.
c    X should not be 0.
c
c    Output, double precision R8_LOG_10, the logarithm base 10 of the absolute
c    value of X.  It should be true that |X| = 10**R_LOG_10.
c
      implicit none

      double precision r8_huge
      double precision r8_log_10
      double precision x

      if ( x .eq. 0.0D+00 ) then
        r8_log_10 = - r8_huge ( x )
      else
        r8_log_10 = log10 ( abs ( x ) )
      end if

      return
      end
      function r8_uniform ( a, b, seed )

c*********************************************************************72
c
cc R8_UNIFORM returns a scaled pseudorandom R8.
c
c  Discussion:
c
c    The pseudorandom number should be uniformly distributed
c    between A and B.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    06 January 2006
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, double precision A, B, the limits of the interval.
c
c    Input/output, integer SEED, the "seed" value, which should NOT be 0.
c    On output, SEED has been updated.
c
c    Output, double precision R8_UNIFORM, a number strictly between A and B.
c
      implicit none

      double precision a
      double precision b
      integer k
      double precision r8_uniform
      integer seed

      if ( seed .eq. 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'R8_UNIFORM - Fatal error!'
        write ( *, '(a)' ) '  Input value of SEED = 0.'
        stop
      end if

      k = seed / 127773

      seed = 16807 * ( seed - k * 127773 ) - k * 2836

      if ( seed .lt. 0 ) then
        seed = seed + 2147483647
      end if
c
c  Although SEED can be represented exactly as a 32 bit integer,
c  it generally cannot be represented exactly as a 32 bit real number!
c
      r8_uniform = a + ( b - a ) * dble ( seed ) * 4.656612875D-10

      return
      end
      function r8_uniform_01 ( seed )

c*********************************************************************72
c
cc R8_UNIFORM_01 returns a unit pseudorandom R8.
c
c  Discussion:
c
c    This routine implements the recursion
c
c      seed = 16807 * seed mod ( 2**31 - 1 )
c      r8_uniform_01 = seed / ( 2**31 - 1 )
c
c    The integer arithmetic never requires more than 32 bits,
c    including a sign bit.
c
c    If the initial seed is 12345, then the first three computations are
c
c      Input     Output      R8_UNIFORM_01
c      SEED      SEED
c
c         12345   207482415  0.096616
c     207482415  1790989824  0.833995
c    1790989824  2035175616  0.947702
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    11 August 2004
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Paul Bratley, Bennett Fox, Linus Schrage,
c    A Guide to Simulation,
c    Second Edition,
c    Springer, 1987,
c    ISBN: 0387964673,
c    LC: QA76.9.C65.B73.
c
c    Bennett Fox,
c    Algorithm 647:
c    Implementation and Relative Efficiency of Quasirandom
c    Sequence Generators,
c    ACM Transactions on Mathematical Software,
c    Volume 12, Number 4, December 1986, pages 362-376.
c
c    Pierre L'Ecuyer,
c    Random Number Generation,
c    in Handbook of Simulation,
c    edited by Jerry Banks,
c    Wiley, 1998,
c    ISBN: 0471134031,
c    LC: T57.62.H37.
c
c    Peter Lewis, Allen Goodman, James Miller,
c    A Pseudo-Random Number Generator for the System/360,
c    IBM Systems Journal,
c    Volume 8, Number 2, 1969, pages 136-143.
c
c  Parameters:
c
c    Input/output, integer SEED, the "seed" value, which should NOT be 0.
c    On output, SEED has been updated.
c
c    Output, double precision R8_UNIFORM_01, a new pseudorandom variate,
c    strictly between 0 and 1.
c
      implicit none

      integer i4_huge
      parameter ( i4_huge = 2147483647 )
      integer k
      double precision r8_uniform_01
      integer seed

      if ( seed .eq. 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'R8_UNIFORM_01 - Fatal error!'
        write ( *, '(a)' ) '  Input value of SEED = 0.'
        stop
      end if

      k = seed / 127773

      seed = 16807 * ( seed - k * 127773 ) - k * 2836

      if ( seed .lt. 0 ) then
        seed = seed + i4_huge
      end if
c
c  Although SEED can be represented exactly as a 32 bit integer,
c  it generally cannot be represented exactly as a 32 bit real number!
c
      r8_uniform_01 = dble ( seed ) * 4.656612875D-10

      return
      end
      subroutine r8mat_print ( m, n, a, title )

c*********************************************************************72
c
cc R8MAT_PRINT prints an R8MAT.
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
c    20 May 2004
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer M, the number of rows in A.
c
c    Input, integer N, the number of columns in A.
c
c    Input, double precision A(M,N), the matrix.
c
c    Input, character ( len = * ) TITLE, a title.
c
      implicit none

      integer m
      integer n

      double precision a(m,n)
      character ( len = * ) title

      call r8mat_print_some ( m, n, a, 1, 1, m, n, title )

      return
      end
      subroutine r8mat_print_some ( m, n, a, ilo, jlo, ihi, jhi,
     &  title )

c*********************************************************************72
c
cc R8MAT_PRINT_SOME prints some of an R8MAT.
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
c    25 January 2007
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer M, N, the number of rows and columns.
c
c    Input, double precision A(M,N), an M by N matrix to be printed.
c
c    Input, integer ILO, JLO, the first row and column to print.
c
c    Input, integer IHI, JHI, the last row and column to print.
c
c    Input, character ( len = * ) TITLE, a title.
c
      implicit none

      integer incx
      parameter ( incx = 5 )
      integer m
      integer n

      double precision a(m,n)
      character * ( 14 ) ctemp(incx)
      integer i
      integer i2hi
      integer i2lo
      integer ihi
      integer ilo
      integer inc
      integer j
      integer j2
      integer j2hi
      integer j2lo
      integer jhi
      integer jlo
      character * ( * ) title

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) trim ( title )

      do j2lo = max ( jlo, 1 ), min ( jhi, n ), incx

        j2hi = j2lo + incx - 1
        j2hi = min ( j2hi, n )
        j2hi = min ( j2hi, jhi )

        inc = j2hi + 1 - j2lo

        write ( *, '(a)' ) ' '

        do j = j2lo, j2hi
          j2 = j + 1 - j2lo
          write ( ctemp(j2), '(i7,7x)') j
        end do

        write ( *, '(''  Col   '',5a14)' ) ( ctemp(j), j = 1, inc )
        write ( *, '(a)' ) '  Row'
        write ( *, '(a)' ) ' '

        i2lo = max ( ilo, 1 )
        i2hi = min ( ihi, m )

        do i = i2lo, i2hi

          do j2 = 1, inc

            j = j2lo - 1 + j2

            write ( ctemp(j2), '(g14.6)' ) a(i,j)

          end do

          write ( *, '(i5,a,5a14)' ) i, ':', ( ctemp(j), j = 1, inc )

        end do

      end do

      return
      end
      subroutine r8mat_transpose_print ( m, n, a, title )

c*********************************************************************72
c
cc R8MAT_TRANSPOSE_PRINT prints an R8MAT, transposed.
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
c    28 April 2008
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer M, N, the number of rows and columns.
c
c    Input, double precision A(M,N), an M by N matrix to be printed.
c
c    Input, character*(*) TITLE, a title.
c
      implicit none

      integer m
      integer n

      double precision a(m,n)
      character*(*) title

      call r8mat_transpose_print_some ( m, n, a, 1, 1, m, n, title )

      return
      end
      subroutine r8mat_transpose_print_some ( m, n, a, ilo, jlo, ihi, 
     &  jhi, title )

c*********************************************************************72
c
cc R8MAT_TRANSPOSE_PRINT_SOME prints some of an R8MAT transposed.
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
c    28 April 2008
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer M, N, the number of rows and columns.
c
c    Input, double precision A(M,N), an M by N matrix to be printed.
c
c    Input, integer ILO, JLO, the first row and column to print.
c
c    Input, integer IHI, JHI, the last row and column to print.
c
c    Input, character * ( * ) TITLE, a title.
c
      implicit none

      integer incx
      parameter ( incx = 5 )
      integer m
      integer n

      double precision a(m,n)
      character * ( 14 ) ctemp(incx)
      integer i
      integer i2
      integer i2hi
      integer i2lo
      integer ihi
      integer ilo
      integer inc
      integer j
      integer j2hi
      integer j2lo
      integer jhi
      integer jlo
      character * ( * ) title

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) trim ( title )

      do i2lo = max ( ilo, 1 ), min ( ihi, m ), incx

        i2hi = i2lo + incx - 1
        i2hi = min ( i2hi, m )
        i2hi = min ( i2hi, ihi )

        inc = i2hi + 1 - i2lo

        write ( *, '(a)' ) ' '

        do i = i2lo, i2hi
          i2 = i + 1 - i2lo
          write ( ctemp(i2), '(i8,6x)') i
        end do

        write ( *, '(''       Row'',5a14)' ) ctemp(1:inc)
        write ( *, '(a)' ) '       Col'

        j2lo = max ( jlo, 1 )
        j2hi = min ( jhi, n )

        do j = j2lo, j2hi

          do i2 = 1, inc
            i = i2lo - 1 + i2
            write ( ctemp(i2), '(g14.6)' ) a(i,j)
          end do

          write ( *, '(2x,i8,a,5a14)' ) j, ':', ( ctemp(i), i = 1, inc )

        end do

      end do

      return
      end
      subroutine r8poly_degree ( na, a, degree )

c*********************************************************************72
c
cc R8POLY_DEGREE returns the degree of a polynomial.
c
c  Discussion:
c
c    The degree of a polynomial is the index of the highest power
c    of X with a nonzero coefficient.
c
c    The degree of a constant polynomial is 0.  The degree of the
c    zero polynomial is debatable, but this routine returns the
c    degree as 0.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    21 March 2001
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer NA, the dimension of A.
c
c    Input, double precision A(0:NA), the coefficients of the polynomials.
c
c    Output, integer DEGREE, the degree of A.
c
      implicit none

      integer na

      double precision a(0:na)
      integer degree

      degree = na

10    continue

      if ( 0 .lt. degree ) then

        if ( a(degree) .ne. 0.0D+00 ) then
          return
        end if

        degree = degree - 1

        go to 10

      end if

      return
      end
      subroutine r8poly_print ( n, a, title )

c*********************************************************************72
c
cc R8POLY_PRINT prints out a polynomial.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    22 July 2010
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the dimension of A.
c
c    Input, double precision A(0:N), the polynomial coefficients.
c    A(0) is the constant term and
c    A(N) is the coefficient of X**N.
c
c    Input, character * ( * ) TITLE, a title.
c
      implicit none

      integer n

      double precision a(0:n)
      integer i
      double precision mag
      integer n2
      character plus_minus
      character * ( * )  title

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) trim ( title )
      write ( *, '(a)' ) ' '

      call r8poly_degree ( n, a, n2 )

      if ( n2 .le. 0 ) then
        write ( *, '( ''  p(x) = 0'' )' )
        return
      end if

      if ( a(n2) .lt. 0.0D+00 ) then
        plus_minus = '-'
      else
        plus_minus = ' '
      end if

      mag = abs ( a(n2) )

      if ( 2 .le. n2 ) then
        write ( *, '( ''  p(x) = '', a1, g14.6, '' * x ^ '', i3 )' ) 
     &    plus_minus, mag, n2
      else if ( n2 .eq. 1 ) then
        write ( *, '( ''  p(x) = '', a1, g14.6, '' * x'' )' ) 
     &    plus_minus, mag
      else if ( n2 .eq. 0 ) then
        write ( *, '( ''  p(x) = '', a1, g14.6 )' ) plus_minus, mag
      end if

      do i = n2 - 1, 0, -1

        if ( a(i) .lt. 0.0D+00 ) then
          plus_minus = '-'
        else
          plus_minus = '+'
        end if

        mag = abs ( a(i) )

        if ( mag .ne. 0.0D+00 ) then

          if ( 2 .le. i ) then
            write ( *, 
     &        ' ( ''         '', a1, g14.6, '' * x ^ '', i3 )' ) 
     &        plus_minus, mag, i
          else if ( i .eq. 1 ) then
            write ( *, 
     &        ' ( ''         '', a1, g14.6, '' * x'' )' )
     &        plus_minus, mag
          else if ( i .eq. 0 ) then
            write ( *, ' ( ''         '', a1, g14.6 )' ) 
     &        plus_minus, mag
          end if
        end if

      end do

      return
      end
      subroutine r8vec_indicator ( n, a )

c*********************************************************************72
c
cc R8VEC_INDICATOR sets an R8VEC to the indicator vector.
c
c  Discussion:
c
c    An R8VEC is a vector of R8's.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    22 January 2007
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the number of elements of A.
c
c    Output, double precision A(N), the array to be initialized.
c
      implicit none

      integer n

      double precision a(n)
      integer i

      do i = 1, n
        a(i) = dble ( i )
      end do

      return
      end
      subroutine r8vec_mean ( n, a, mean )

c*********************************************************************72
c
cc R8VEC_MEAN returns the mean of an R8VEC.
c
c  Discussion:
c
c    An R8VEC is a vector of R8's.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    29 May 2009
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the number of entries in the vector.
c
c    Input, double precision A(N), the vector whose mean is desired.
c
c    Output, double precision MEAN, the mean of the vector entries.
c
      implicit none

      integer n

      double precision a(n)
      integer i
      double precision mean

      mean = 0.0D+00
      do i = 1, n
        mean = mean + a(i)
      end do
      mean = mean / dble ( n )

      return
      end
      subroutine r8vec_print ( n, a, title )

c*********************************************************************72
c
cc R8VEC_PRINT prints an R8VEC.
c
c  Discussion:
c
c    An R8VEC is a vector of R8's.
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
c    Input, integer N, the number of components of the vector.
c
c    Input, double precision A(N), the vector to be printed.
c
c    Input, character * ( * ) TITLE, a title.
c
      implicit none

      integer n

      double precision a(n)
      integer i
      character ( len = * ) title

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) trim ( title )
      write ( *, '(a)' ) ' '
      do i = 1, n
        write ( *, '(2x,i8,a,1x,g16.8)' ) i, ':', a(i)
      end do

      return
      end
      subroutine r8vec_variance ( n, a, variance )

c*********************************************************************72
c
cc R8VEC_VARIANCE returns the variance of an R8VEC.
c
c  Discussion:
c
c    An R8VEC is a vector of R8's.
c
c    The variance of a vector X of length N is defined as
c
c      mean ( X(1:n) ) = sum ( X(1:n) ) / n
c
c      var ( X(1:n) ) = sum ( ( X(1:n) - mean )**2 ) / ( n - 1 )
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    04 July 2009
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the number of entries in the vector.
c    N should be at least 2.
c
c    Input, double precision A(N), the vector.
c
c    Output, double precision VARIANCE, the variance of the vector.
c
      implicit none

      integer n

      double precision a(n)
      integer i
      double precision mean
      double precision variance

      if ( n .lt. 2 ) then

        variance = 0.0D+00
 
      else

        mean = 0.0D+00
        do i = 1, n
          mean = mean + a(i)
        end do
        mean = mean / dble ( n )

        variance = 0.0D+00
        do i = 1, n
          variance = variance + ( a(i) - mean )**2
        end do
        variance = variance / dble ( n - 1 )

      end if

      return
      end
      function radians_to_degrees ( radians )

c*********************************************************************72
c
cc RADIANS_TO_DEGREES converts an angle measure from radians to degrees.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    11 December 2004
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, double precision RADIANS, the angle measure in radians.
c
c    Output, double precision RADIANS_TO_DEGREES, the angle measure in degrees.
c
      implicit none
 
      double precision pi
      parameter ( pi = 3.141592653589793D+00 )
      double precision radians
      double precision radians_to_degrees
 
      radians_to_degrees = ( radians / pi ) * 180.0D+00

      return
      end
      subroutine rat_factor ( m, n, maxfactor, factor_num, factor, 
     &  power, mleft, nleft )

c*********************************************************************72
c
cc RAT_FACTOR factors a rational value into a product of prime factors.
c
c  Discussion:
c
c    ( M / N ) = ( MLEFT / NLEFT ) * Product ( 1 <= I <= FACTOR_NUM )
c      FACTOR(I)**POWER(I).
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    22 July 2010
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer M, N, the top and bottom of a rational value.
c    The ratio of M and N must be positive.
c
c    Input, integer MAXFACTOR, the maximum number of factors for
c    which storage has been allocated.
c
c    Output, integer FACTOR_NUM, the number of prime factors 
c    of M/N.
c
c    Output, integer FACTOR(MAXFACTOR), the prime factors of M/N.
c
c    Output, integer POWER(MAXFACTOR).  POWER(I) is the power of
c    the FACTOR(I) in the representation of M/N.
c
c    Output, integer MLEFT, NLEFT, the top and bottom of 
c    the factor of M / N that remains.  If ABS ( MLEFT / NLEFT ) is not 1, then
c    the rational value was not completely factored.
c
      implicit none

      integer maxfactor

      integer factor(maxfactor)
      integer factor_num
      integer i
      integer ( kind = 4 ) m
      integer ( kind = 4 ) mleft
      integer ( kind = 4 ) n
      integer ( kind = 4 ) nleft
      integer ( kind = 4 ) p
      integer ( kind = 4 ) power(maxfactor)
      integer ( kind = 4 ) prime
      integer ( kind = 4 ) prime_max

      factor_num = 0

      mleft = m
      nleft = n
c
c  NLEFT should be nonnegative.
c
      if ( nleft .lt. 0 ) then
        mleft = -mleft
        nleft = -nleft
      end if

      if ( m .eq. 0 .or. n .eq. 0 ) then
        return
      end if

      if ( m .eq. n ) then
        factor_num = 1
        factor(1) = 1
        power(1) = 1
        return
      end if
c
c  Find out how many primes we stored.
c
      prime_max = prime ( -1 )

      do i = 1, prime_max

        p = prime ( i )

        if ( mod ( nleft, p ) .eq. 0 .or. 
     &       mod ( abs ( mleft ), p ) .eq. 0 ) then

          if ( factor_num .lt. maxfactor ) then

            factor_num = factor_num + 1
            factor(factor_num) = p
            power(factor_num) = 0
c
c  Divide MLEFT by PRIME(I) as often as you can.
c
            if ( mod ( abs ( mleft ), p ) .eq. 0  ) then

10            continue
 
                power(factor_num) = power(factor_num) + 1
                mleft = mleft / p

                if ( mod ( abs ( mleft ), p ) .ne. 0 ) then
                  go to 20
                end if

              go to 10

20            continue

            end if
c
c  Divide NLEFT by PRIME(I) as often as you can.
c
            if ( mod ( nleft, p ) .eq. 0  ) then

30            continue

                power(factor_num) = power(factor_num) - 1
                nleft = nleft / p

                if ( mod ( nleft, p ) .ne. 0 ) then
                  go to 40
                end if

              go to 30

40            continue

            end if

            if ( power(factor_num) .eq. 0 ) then
              factor_num = factor_num - 1
            end if

          end if

        end if

      end do

      return
      end
      subroutine rickey ( ab, bb, er, f, h, hb, hp, r, so, tb, g )

c*********************************************************************72
c
cc RICKEY evaluates Branch Rickey's baseball index.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    22 July 2010
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Alan Schwarz,
c    Looking Beyond the Batting Average,
c    The New York Times, 
c    Sunday, 1 August 2004.
c    
c    Branch Rickey,
c    Goodby to Some Old Baseball Ideas,
c    Life Magazine,
c    2 August 1954.
c
c  Parameters:
c
c    Input, integer AB, number of at-bats.
c
c    Input, integer BB, base on balls.
c
c    Input, integer ER, earned runs.
c
c    Input, double precision F, the fielding rating.
c
c    Input, integer H, number of hits.
c 
c    Input, integer HB, hit batsmen.
c
c    Input, integer HP, hit by pitcher.
c
c    Input, integer R, runs.
c
c    Input, integer SO, strike outs.
c
c    Input, integer TB, total bases.
c
c    Output, double precision G, the Branch Rickey index, an estimate for the
c    expected winning percentage of a team with the given statistics.
c    (0.5 has already been subtracted from this value.)
c
      implicit none

      integer ab
      integer bb
      integer er
      double precision f
      double precision g
      integer h
      integer hb
      double precision hitting
      integer hp
      double precision pitching
      integer r
      integer so
      integer tb

      hitting = 
     &    dble (  h + bb + hp   ) / dble ( ab + bb + hp ) 
     &  + dble ( 3 * ( tb - h ) ) / dble ( 4 * ab       ) 
     &  + dble (              r ) / dble ( h + bb + hp  )

      pitching = 
     &    dble ( h       ) / dble ( ab           ) 
     &  + dble ( bb + hb ) / dble ( ab + bb + hb ) 
     &  + dble ( er      ) / dble ( h + bb + hb  ) 
     &  - dble ( so      ) / dble ( 8 * ( ab + bb + hb ) )

      g = hitting - pitching - f

      return
      end
      subroutine roots_to_i4poly ( n, x, c )

c*********************************************************************72
c
cc ROOTS_TO_I4POLY converts polynomial roots to polynomial coefficients.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    22 July 2010
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the number of roots specified.
c
c    Input, integer X(N), the roots.
c
c    Output, integer C(0:N), the coefficients of the polynomial.
c
      implicit none

      integer n

      integer c(0:n)
      integer i
      integer j
      integer x(n)
!
!  Initialize C to (0, 0, ..., 0, 1).
!  Essentially, we are setting up a divided difference table.
!
      do i = 0, n - 1
        c(i) = 0
      end do
      c(n) = 1
!
!  Convert to standard polynomial form by shifting the abscissas
!  of the divided difference table to 0.
!
      do j = 1, n
        do i = 1, n + 1 - j
          c(n-i) = c(n-i) - x(n+1-i-j+1) * c(n-i+1)
        end do
      end do

      return
      end
      subroutine roots_to_r8poly ( n, x, c )

c*********************************************************************72
c
cc ROOTS_TO_R8POLY converts polynomial roots to polynomial coefficients.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    17 July 2010
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the number of roots specified.
c
c    Input, double precision X(N), the roots.
c
c    Output, double precision C(0:N), the coefficients of the polynomial.
c
      implicit none

      integer n

      double precision c(0:n)
      integer i
      integer j
      double precision x(n)
c
c  Initialize C to (0, 0, ..., 0, 1).
c  Essentially, we are setting up a divided difference table.
c
      do i = 0, n - 1
        c(i) = 0.0D+00
      end do
      c(n) = 1.0D+00
c
c  Convert to standard polynomial form by shifting the abscissas
c  of the divided difference table to 0.
c
      do j = 1, n
        do i = 1, n + 1 - j
          c(n-i) = c(n-i) - x(n+1-i-j+1) * c(n-i+1)
        end do
      end do

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
      subroutine sort_heap_external ( n, indx, i, j, isgn )

c*********************************************************************72
c
cc SORT_HEAP_EXTERNAL externally sorts a list of items into ascending order.
c
c  Discussion:
c
c    The actual list of data is not passed to the routine.  Hence this
c    routine may be used to sort integers, reals, numbers, names,
c    dates, shoe sizes, and so on.  After each call, the routine asks
c    the user to compare or interchange two items, until a special
c    return value signals that the sorting is completed.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    25 January 2007
c
c  Author:
c
c    Original FORTRAN77 version by Albert Nijenhuis, Herbert Wilf.
c    This FORTRAN77 version by John Burkardt.
c
c  Reference:
c
c    Albert Nijenhuis, Herbert Wilf,
c    Combinatorial Algorithms for Computers and Calculators,
c    Second Edition,
c    Academic Press, 1978,
c    ISBN: 0-12-519260-6,
c    LC: QA164.N54.
c
c  Parameters:
c
c    Input, integer N, the number of items to be sorted.
c
c    Input/output, integer INDX, the main communication signal.
c
c    The user must set INDX to 0 before the first call.
c    Thereafter, the user should not change the value of INDX until
c    the sorting is done.
c
c    On return, if INDX is
c
c      greater than 0,
c      * interchange items I and J;
c      * call again.
c
c      less than 0,
c      * compare items I and J;
c      * set ISGN = -1 if I .lt. J, ISGN = +1 if J .lt. I;
c      * call again.
c
c      equal to 0, the sorting is done.
c
c    Output, integer I, J, the indices of two items.
c    On return with INDX positive, elements I and J should be interchanged.
c    On return with INDX negative, elements I and J should be compared, and
c    the result reported in ISGN on the next call.
c
c    Input, integer ISGN, results of comparison of elements I and J.
c    (Used only when the previous call returned INDX less than 0).
c    ISGN .le. 0 means I is less than or equal to J;
c    0 .le. ISGN means I is greater than or equal to J.
c
      implicit none

      integer i
      integer i_save
      integer indx
      integer isgn
      integer j
      integer j_save
      integer k
      integer k1
      integer n
      integer n1

      save i_save
      save j_save
      save k
      save k1
      save n1

      data i_save / 0 /
      data j_save / 0 /
      data k / 0 /
      data k1 / 0 /
      data n1 / 0 /
c
c  INDX = 0: This is the first call.
c
      if ( indx .eq. 0 ) then

        i_save = 0
        j_save = 0
        k = n / 2
        k1 = k
        n1 = n
c
c  INDX .lt. 0: The user is returning the results of a comparison.
c
      else if ( indx .lt. 0 ) then

        if ( indx .eq. -2 ) then

          if ( isgn .lt. 0 ) then
            i_save = i_save + 1
          end if

          j_save = k1
          k1 = i_save
          indx = -1
          i = i_save
          j = j_save
          return

        end if

        if ( 0 .lt. isgn ) then
          indx = 2
          i = i_save
          j = j_save
          return
        end if

        if ( k .le. 1 ) then

          if ( n1 .eq. 1 ) then
            i_save = 0
            j_save = 0
            indx = 0
          else
            i_save = n1
            n1 = n1 - 1
            j_save = 1
            indx = 1
          end if

          i = i_save
          j = j_save
          return

        end if

        k = k - 1
        k1 = k
c
c  0 .lt. INDX, the user was asked to make an interchange.
c
      else if ( indx .eq. 1 ) then

        k1 = k

      end if

10    continue

        i_save = 2 * k1

        if ( i_save .eq. n1 ) then
          j_save = k1
          k1 = i_save
          indx = -1
          i = i_save
          j = j_save
          return
        else if ( i_save .le. n1 ) then
          j_save = i_save + 1
          indx = -2
          i = i_save
          j = j_save
          return
        end if

        if ( k .le. 1 ) then
          go to 20
        end if

        k = k - 1
        k1 = k

      go to 10

20    continue

      if ( n1 .eq. 1 ) then
        i_save = 0
        j_save = 0
        indx = 0
        i = i_save
        j = j_save
      else
        i_save = n1
        n1 = n1 - 1
        j_save = 1
        indx = 1
        i = i_save
        j = j_save
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
      subroutine tuple_next2 ( n, xmin, xmax, rank, x )

c*********************************************************************72
c
cc TUPLE_NEXT2 computes the next element of an integer tuple space.
c
c  Discussion:
c
c    The elements X are N vectors.
c
c    Each entry X(I) is constrained to lie between XMIN(I) and XMAX(I).
c
c    The elements are produced one at a time.
c
c    The first element is
c      (XMIN(1), XMIN(2), ..., XMIN(N)),
c    the second is (probably)
c      (XMIN(1), XMIN(2), ..., XMIN(N)+1),
c    and the last element is
c      (XMAX(1), XMAX(2), ..., XMAX(N))
c
c    Intermediate elements are produced in a lexicographic order, with
c    the first index more important than the last, and the ordering of
c    values at a fixed index implicitly defined by the sign of
c    XMAX(I) - XMIN(I).
c
c  Example:
c
c    N = 2,
c    XMIN = (/ 1, 10 /)
c    XMAX = (/ 3,  8 /)
c
c    RANK    X
c    ----  -----
c      1   1 10
c      2   1  9
c      3   1  8
c      4   2 10
c      5   2  9
c      6   2  8
c      7   3 10
c      8   3  9
c      9   3  8
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    25 January 2007
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the number of components.
c
c    Input, integer XMIN(N), XMAX(N), the "minimum" and "maximum" entry values.
c    These values are minimum and maximum only in the sense of the lexicographic
c    ordering.  In fact, XMIN(I) may be less than, equal to, or greater
c    than XMAX(I).
c
c    Input/output, integer RANK, the rank of the item.  On first call,
c    set RANK to 0 to start up the sequence.  On return, if RANK is zero,
c    there are no more items in the sequence.
c
c    Input/output, integer X(N), on input the previous tuple.
c    On output, the next tuple.
c
      implicit none

      integer n

      integer i
      integer prod
      integer rank
      integer x(n)
      integer xmin(n)
      integer xmax(n)

      if ( rank .lt. 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'TUPLE_NEXT2 - Fatal error!'
        write ( *, '(a,i8)' ) '  Illegal value of RANK = ', rank
        stop
      end if

      prod = 1
      do i = 1, n
        prod = prod * ( 1 + abs ( xmax(i) - xmin(i) ) )
      end do

      if ( prod .lt. rank ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'TUPLE_NEXT2 - Fatal error!'
        write ( *, '(a,i8)' ) '  Illegal value of RANK = ', rank
        stop
      end if

      if ( rank .eq. 0 ) then
        do i = 1, n
          x(i) = xmin(i)
        end do
        rank = 1
        return
      end if

      rank = rank + 1
      i = n

10    continue

        if ( x(i) .ne. xmax(i) ) then
          x(i) = x(i) + sign ( 1, xmax(i) - xmin(i) )
          go to 20
        end if

        x(i) = xmin(i)

        if ( i .eq. 1 ) then
          rank = 0
          go to 20
        end if

        i = i - 1

      go to 10

20    continue

      return
      end
      subroutine tvec_even ( nt, t )

c*********************************************************************72
c
cc TVEC_EVEN computes evenly spaced angles between 0 and 2*PI.
c
c  Discussion:
c
c    The computation realizes that 0 = 2 * PI.
c
c  Example:
c
c    NT = 4
c
c    T = ( 0, PI/2, PI, 3*PI/2 )
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    07 June 2008
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer NT, the number of values to compute.
c
c    Output, double precision TVEC(NT), the evenly spaced angles, in radians.
c
      implicit none

      integer nt

      integer i
      double precision pi
      parameter ( pi = 3.141592653589793D+00 )
      double precision t(nt)

      do i = 1, nt
        t(i) = dble ( 2 * ( i - 1 ) ) * pi / dble ( nt )
      end do

      return
      end
      subroutine tvec_even2 ( nt, t )

c*********************************************************************72
c
cc TVEC_EVEN2 computes evenly spaced angles between 0 and 2*PI.
c
c  Discussion:
c
c    The computation realizes that 0 = 2 * PI.  The values are equally
c    spaced in the circle, do not include 0, and are symmetric about 0.
c
c  Example:
c
c    NT = 4
c
c    T = ( PI/4, 3*PI/4, 5*PI/4, 7*PI/4 )
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    07 June 2008
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer NT, the number of values to compute.
c
c    Output, double precision TVEC(NT), the evenly spaced angles, in radians.
c
      implicit none

      integer nt

      integer i
      double precision pi
      parameter ( pi = 3.141592653589793D+00 )
      double precision t(nt)

      do i = 1, nt
        t(i) = dble ( 2 * i - 1 ) * pi / dble ( nt )
      end do

      return
      end
      subroutine tvec_even3 ( nt, t )

c*********************************************************************72
c
cc TVEC_EVEN3 computes evenly spaced angles between 0 and 2*PI.
c
c  Discussion:
c
c    The angles begin with 0 and end with 2*PI.
c
c  Example:
c
c    NT = 4
c
c    T = ( 0, 2*PI/3, 4*PI/3 2*PI )
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    07 June 2008
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer NT, the number of values to compute.
c
c    Output, double precision TVEC(NT), the evenly spaced angles, in radians.
c
      implicit none

      integer nt

      integer i
      double precision pi
      parameter ( pi = 3.141592653589793D+00 )
      double precision t(nt)

      if ( nt .eq. 1 ) then
        t(1) = pi
      else
        do i = 1, nt
          t(i) = dble ( 2 * ( i - 1 ) ) * pi / dble ( nt - 1 )
        end do
      end if

      return
      end
      subroutine tvec_even_bracket ( nt, theta1, theta2, t )

c*********************************************************************72
c
cc TVEC_EVEN_BRACKET computes evenly spaced angles between THETA1 and THETA2.
c
c  Discussion:
c
c    The interval between THETA1 and THETA2 is divided into NT-1 subintervals.
c
c    The angles returned are the breakpoints of these subintervals,
c    including THETA1 and THETA2.
c
c  Example:
c
c    NT = 4
c    THETA1 = 30
c    THETA2 = 90
c
c    T = ( 30, 50, 70, 90 )
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    07 June 2008
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer NT, the number of values to compute.
c
c    Input, double precision THETA1, THETA2, the limiting angles.
c
c    Output, double precision TVEC(NT), the evenly spaced angles.
c
      implicit none

      integer nt

      integer i
      double precision t(nt)
      double precision theta1
      double precision theta2

      if ( nt .eq. 1 ) then

        t(1) = ( theta1 + theta2 ) / 2.0D+00

      else

        do i = 1, nt
          t(i) = ( dble ( nt - i     ) * theta1   
     &           + dble (      i - 1 ) * theta2 ) 
     &           / dble ( nt     - 1 )
        end do

      end if

      return
      end
      subroutine tvec_even_bracket2 ( nt, theta1, theta2, t )

c*********************************************************************72
c
cc TVEC_EVEN_BRACKET2 computes evenly spaced angles from THETA1 to THETA2.
c
c  Discussion:
c
c    The interval between THETA1 and THETA2 is divided into NT+1 subintervals.
c
c    The angles returned are the internal NT breakpoints of the subintervals.
c
c  Example:
c
c    NT = 5
c    THETA1 = 30
c    THETA2 = 90
c
c    T = ( 40, 50, 60, 70, 80 )
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    30 October 2007
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer NT, the number of values to compute.
c
c    Input, double precision THETA1, THETA2, the limiting angles.
c
c    Output, double precision TVEC(NT), the evenly spaced angles.
c
      implicit none

      integer nt

      integer i
      double precision t(nt)
      double precision theta1
      double precision theta2

      do i = 1, nt
        t(i) = ( dble ( nt + 1 - i ) * theta1   
     &         + dble (          i ) * theta2 ) 
     &         / dble ( nt + 1     )
      end do

      return
      end
      subroutine tvec_even_bracket3 ( nt, theta1, theta2, t )

c*********************************************************************72
c
cc TVEC_EVEN_BRACKET3 computes evenly spaced angles between THETA1 and THETA2.
c
c  Discussion:
c
c    The interval between THETA1 and THETA2 is divided into NT subintervals.
c
c    The angles returned are the midpoints of each subinterval.
c
c  Example:
c
c    NT = 3
c    THETA1 = 30
c    THETA2 = 90
c
c    T = ( 40, 60, 80 )
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    01 January 2007
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer NT, the number of values to compute.
c
c    Input, double precision THETA1, THETA2, the limiting angles.
c
c    Output, double precision TVEC(NT), the evenly spaced angles.
c
      implicit none

      integer nt

      integer i
      double precision t(nt)
      double precision theta1
      double precision theta2

      do i = 1, nt
        t(i) = ( dble ( 2 * nt - 2 * i + 1 ) * theta1   
     &         + dble (          2 * i - 1 ) * theta2 ) 
     &         / dble ( 2 * nt             )
      end do

      return
      end
      subroutine upc_check_digit ( p, l, r, c )

c*********************************************************************72
c
cc UPC_CHECK_DIGIT returns the check digit of a UPC.
c
c  Discussion:
c
c    UPC stands for Universal Price Code.
c
c    A full UPC is a string of 12 digits, in groups of size 1, 5, 5, and 1,
c    of the form P-LLLLL-RRRRR-C, where:
c
c      P is the one-digit product type code.
c      L is the five-digit manufacturer code.
c      R is the five_digit product code
c      C is the check digit.
c
c  Example:
c
c    0-72890-00011-8
c    0-12345-67890-5
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    01 July 2008
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer P, the one-digit product type code.
c
c    Input, integer L, the five-digit manufacturer code.
c
c    Input, integer R, the five-digit product code.
c
c    Output, integer C, the check digit.
c
      implicit none

      integer c
      integer i4_ten
      parameter ( i4_ten = 10 )
      integer l
      integer lc(5)
      integer p
      integer r
      integer rc(5)

      if ( p .lt. 0 .or. 9 .lt. p ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'UPC_CHECK_DIGIT - Fatal error!'
        write ( *, '(a)' ) '  P .lt. 0 or 9 .lt. P.'
        stop
      end if

      if ( l .lt. 0 .or. 99999 .lt. l ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'UPC_CHECK_DIGIT - Fatal error!'
        write ( *, '(a)' ) '  L .lt. 0 or 99999 .lt. L.'
        stop
      end if

      if ( r .lt. 0 .or. 99999 .lt. r ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'UPC_CHECK_DIGIT - Fatal error!'
        write ( *, '(a)' ) '  R .lt. 0 or 99999 .lt. R.'
        stop
      end if

      call i4_to_digits_decimal ( l, 5, lc )
      call i4_to_digits_decimal ( r, 5, rc )

      c = ( p + lc(2) + lc(4) + rc(1) + rc(3) + rc(5) ) * 3 
     &        + lc(1) + lc(3) + lc(5) + rc(2) + rc(4)

      c = mod ( c, i4_ten )

      c = mod ( i4_ten - c, i4_ten )

      return
      end
      function versine_pulse ( t, ta, tb, v1, amp )

c*********************************************************************72
c
cc VERSINE_PULSE adds a versine pulse to a constant.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    01 March 2010
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, double precision T, the current time.
c
c    Input, double precision TA, the time at which the pulse begins.
c
c    Input, double precision TB, the time at which the pulse finishes.
c
c    Input, double precision V1, the constant value.
c
c    Input, double precision AMP, the amplitude of the pulse.
c
c    Output, double precision VERSINE_PULSE, the value of the signal at time T.
c
      implicit none

      double precision amp
      double precision angle
      double precision pi
      double precision t
      double precision ta
      double precision tb
      double precision v
      double precision v1
      double precision versine_pulse

      parameter ( pi = 3.141592653589793D+00 )

      v = v1

      if ( ta .le. t .and. t .le. tb ) then
        angle = 2.0D+00 * pi * ( t - ta ) / ( tb - ta ) 
        v = v + 0.5D+00 * amp * ( 1.0D+00 - cos ( angle ) )
      end if

      versine_pulse = v

      return
      end
