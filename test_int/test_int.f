      function besj0 ( x )

c*********************************************************************72
c
cc BESJ0 evaluates the Bessel J0(X) function.
c
c  Discussion:
c
c    This routine computes approximate values for Bessel functions
c    of the first kind of order zero for arguments  |X| .le. XMAX
c
c    See comments heading CALJY0.
c
c  Modified:
c
c    03 April 2007
c
c  Author:
c
c    William Cody
c
c  Parameters:
c
c    Input, double precision X, the argument of the function.
c
c    Output, double precision BESJ0, the value of the function.
c
      implicit none

      double precision besj0
      integer jint
      double precision result
      double precision x

      jint = 0
      call caljy0 ( x, result, jint )
      besj0 = result

      return
      end
      function besy0 ( x )

c*********************************************************************72
c
cc BESY0 evaluates the Bessel Y0(X) function.
c
c  Discussion:
c
c    This routine computes approximate values for Bessel functions
c    of the second kind of order zero for arguments 0 < X .le. XMAX.
c
c    See comments heading CALJY0.
c
c  Modified:
c
c    03 April 2007
c
c  Author:
c
c    William Cody
c
c  Parameters:
c
c    Input, double precision X, the argument of the function.
c
c    Output, double precision BESY0, the value of the function.
c
      implicit none

      double precision besy0
      integer jint
      double precision result
      double precision x

      jint = 1
      call caljy0 ( x, result, jint )
      besy0 = result

      return
      end
      subroutine caljy0 ( arg, result, jint )

c*********************************************************************72
c
cc CALJY0 computes various J0 and Y0 Bessel functions.
c
c  Discussion:
c
c    This routine computes zero-order Bessel functions of the first and
c    second kind (J0 and Y0), for real arguments X, where 0 < X .le. XMAX
c    for Y0, and |X| .le. XMAX for J0.  
c
c  Modified:
c
c    03 April 2007
c
c  Author:
c
c    William Cody
c
c  Reference:
c
c    John Hart, Ward Cheney, Charles Lawson, Hans Maehly, 
c    Charles Mesztenyi, John Rice, Henry Thatcher, 
c    Christoph Witzgall,
c    Computer Approximations,
c    Wiley, 1968,
c    LC: QA297.C64.
c
c  Parameters:
c
c    Input, double precision ARG, the argument.  If JINT = 0, ARG
c    must satisfy 
c     -XMAX < ARG < XMAX;
c    If JINT = 1, then ARG must satisfy
c      0 < ARG < XMAX.
c
c    Output, double precision RESULT, the value of the function,
c    which depends on the input value of JINT:
c    0, RESULT = J0(x);
c    1, RESULT = Y0(x);
c
c    Input, integer JINT, chooses the function to be computed.
c    0, J0(x);
c    1, Y0(x);
c
      implicit none

      integer i
      integer jint
      double precision arg
      double precision ax
      double precision cons
      double precision down
      double precision eight
      double precision five5
      double precision four
      double precision one
      double precision oneov8
      double precision pi2
      double precision pj0(7)
      double precision pj1(8)
      double precision plg(4)
      double precision prod
      double precision py0(6)
      double precision py1(7)
      double precision py2(8)
      double precision p0(6)
      double precision p1(6)
      double precision p17
      double precision qj0(5)
      double precision qj1(7)
      double precision qlg(4)
      double precision qy0(5)
      double precision qy1(6)
      double precision qy2(7)
      double precision q0(5)
      double precision q1(5)
      double precision resj
      double precision result
      double precision r0
      double precision r1
      double precision sixty4
      double precision three
      double precision twopi
      double precision twopi1
      double precision twopi2
      double precision two56
      double precision up
      double precision w
      double precision wsq
      double precision xden
      double precision xinf
      double precision xmax
      double precision xnum
      double precision xsmall
      double precision xj0
      double precision xj1
      double precision xj01
      double precision xj02
      double precision xj11
      double precision xj12
      double precision xy
      double precision xy0
      double precision xy01
      double precision xy02
      double precision xy1
      double precision xy11
      double precision xy12
      double precision xy2
      double precision xy21
      double precision xy22
      double precision z
      double precision zero
      double precision zsq
c
c  Mathematical constants
c  CONS = ln(.5) + Euler's gamma
c
      data zero / 0.0d0 /
      data one /1.0d0 /
      data three /3.0d0 /
      data four /4.0d0 /
      data eight /8.0d0/
      data five5 / 5.5d0 /
      data sixty4 /64.0d0 /
      data oneov8 /0.125d0 /
      data p17 /1.716d-1/
      data two56 /256.0d0/
      data cons / -1.1593151565841244881d-1/
      data pi2 /6.3661977236758134308d-1/
      data twopi /6.2831853071795864769d0/
      data twopi1 /6.28125d0 /
      data twopi2 / 1.9353071795864769253d-3/
c
c  Machine-dependent constants
c
      data xmax /1.07d+09/
      data xsmall /9.31d-10/
      data xinf /1.7d+38/
c
c  Zeroes of Bessel functions
c
      data xj0 /2.4048255576957727686d+0/
      data xj1 /5.5200781102863106496d+0/
      data xy0 /8.9357696627916752158d-1/
      data xy1 /3.9576784193148578684d+0/
      data xy2 /7.0860510603017726976d+0/
      data xj01 / 616.0d+0/
      data xj02 /-1.4244423042272313784d-03/
      data xj11 /1413.0d+0/
      data xj12 / 5.4686028631064959660d-04/
      data xy01 / 228.0d+0/
      data xy02 / 2.9519662791675215849d-03/
      data xy11 /1013.0d+0/
      data xy12 / 6.4716931485786837568d-04/
      data xy21 /1814.0d+0/
      data xy22 / 1.1356030177269762362d-04/
c
c  Coefficients for rational approximation to ln(x/a)
c
      data plg/-2.4562334077563243311d+01,2.3642701335621505212d+02,
     &         -5.4989956895857911039d+02,3.5687548468071500413d+02/
      data qlg/-3.5553900764052419184d+01,1.9400230218539473193d+02,
     &         -3.3442903192607538956d+02,1.7843774234035750207d+02/
c
c  Coefficients for rational approximation of
c  J0(X) / (X**2 - XJ0**2),  XSMALL < |X| .le. 4.0
c
      data pj0/6.6302997904833794242d+06,-6.2140700423540120665d+08,
     &         2.7282507878605942706d+10,-4.1298668500990866786d+11,
     &        -1.2117036164593528341d-01, 1.0344222815443188943d+02,
     &        -3.6629814655107086448d+04/
      data qj0/4.5612696224219938200d+05, 1.3985097372263433271d+08,
     &         2.6328198300859648632d+10, 2.3883787996332290397d+12,
     &         9.3614022392337710626d+02/
c
c  Coefficients for rational approximation of
c  J0(X) / (X**2 - XJ1**2), 4.0 < |X| .le. 8.0
c
      data pj1/4.4176707025325087628d+03, 1.1725046279757103576d+04,
     &         1.0341910641583726701d+04,-7.2879702464464618998d+03,
     &        -1.2254078161378989535d+04,-1.8319397969392084011d+03,
     &         4.8591703355916499363d+01, 7.4321196680624245801d+02/
      data qj1/3.3307310774649071172d+02,-2.9458766545509337327d+03,
     &         1.8680990008359188352d+04,-8.4055062591169562211d+04,
     &         2.4599102262586308984d+05,-3.5783478026152301072d+05,
     &        -2.5258076240801555057d+01/
c
c  Coefficients for rational approximation of
c  (Y0(X) - 2 LN(X/XY0) J0(X)) / (X**2 - XY0**2),
c  XSMALL < |X| .le. 3.0
c
      data py0/1.0102532948020907590d+04,-2.1287548474401797963d+06,
     &         2.0422274357376619816d+08,-8.3716255451260504098d+09,
     &         1.0723538782003176831d+11,-1.8402381979244993524d+01/
      data qy0/6.6475986689240190091d+02, 2.3889393209447253406d+05,
     &         5.5662956624278251596d+07, 8.1617187777290363573d+09,
     &         5.8873865738997033405d+11/
c
c  Coefficients for rational approximation of
c  (Y0(X) - 2 LN(X/XY1) J0(X)) / (X**2 - XY1**2),
c  3.0 < |X| .le. 5.5
c
      data py1/-1.4566865832663635920d+04, 4.6905288611678631510d+06,
     &         -6.9590439394619619534d+08, 4.3600098638603061642d+10,
     &         -5.5107435206722644429d+11,-2.2213976967566192242d+13,
     &          1.7427031242901594547d+01/
      data qy1/ 8.3030857612070288823d+02, 4.0669982352539552018d+05,
     &          1.3960202770986831075d+08, 3.4015103849971240096d+10,
     &          5.4266824419412347550d+12, 4.3386146580707264428d+14/
c
c  Coefficients for rational approximation of
c  (Y0(X) - 2 LN(X/XY2) J0(X)) / (X**2 - XY2**2),
c  5.5 < |X| .le. 8.0
c
      data py2/ 2.1363534169313901632d+04,-1.0085539923498211426d+07,
     &          2.1958827170518100757d+09,-1.9363051266772083678d+11,
     &         -1.2829912364088687306d+11, 6.7016641869173237784d+14,
     &         -8.0728726905150210443d+15,-1.7439661319197499338d+01/
      data qy2/ 8.7903362168128450017d+02, 5.3924739209768057030d+05,
     &          2.4727219475672302327d+08, 8.6926121104209825246d+10,
     &          2.2598377924042897629d+13, 3.9272425569640309819d+15,
     &          3.4563724628846457519d+17/
c
c  Coefficients for Hart,s approximation, 8.0 < |X|.
c
      data p0/3.4806486443249270347d+03, 2.1170523380864944322d+04,
     &        4.1345386639580765797d+04, 2.2779090197304684302d+04,
     &        8.8961548424210455236d-01, 1.5376201909008354296d+02/
      data q0/3.5028735138235608207d+03, 2.1215350561880115730d+04,
     &        4.1370412495510416640d+04, 2.2779090197304684318d+04,
     &        1.5711159858080893649d+02/
      data p1/-2.2300261666214198472d+01,-1.1183429920482737611d+02,
     &        -1.8591953644342993800d+02,-8.9226600200800094098d+01,
     &        -8.8033303048680751817d-03,-1.2441026745835638459d+00/
      data q1/1.4887231232283756582d+03, 7.2642780169211018836d+03,
     &        1.1951131543434613647d+04, 5.7105024128512061905d+03,
     &        9.0593769594993125859d+01/
c
c  Check for error conditions.
c
      ax = abs ( arg )

      if ( jint .eq. 1 .and. arg .le. zero ) then
        result = -xinf
        return
      else if ( xmax .lt. ax ) then
        result = zero
        return
      end if

      if ( eight .lt. ax ) then
        go to 800
      end if

      if ( ax .le. xsmall ) then
        if ( jint .eq. 0 ) then
          result = one
        else
          result = pi2 * ( log ( ax ) + cons )
        end if
        return
      end if
c
c  Calculate J0 for appropriate interval, preserving
c  accuracy near the zero of J0.
c
      zsq = ax * ax

      if ( ax .le. four ) then
        xnum = ( pj0(5) * zsq + pj0(6) ) * zsq + pj0(7)
        xden = zsq + qj0(5)
        do i = 1, 4
          xnum = xnum * zsq + pj0(i)
          xden = xden * zsq + qj0(i)
        end do
        prod = ( ( ax - xj01 / two56 ) - xj02 ) * ( ax + xj0 )
      else
        wsq = one - zsq / sixty4
        xnum = pj1(7) * wsq + pj1(8)
        xden = wsq + qj1(7)
        do i = 1, 6
          xnum = xnum * wsq + pj1(i)
          xden = xden * wsq + qj1(i)
        end do
        prod = ( ax + xj1 ) * ( ( ax - xj11 / two56 ) - xj12 )
      end if

      result = prod * xnum / xden

      if ( jint .eq. 0 ) then
        return
      end if
c
c  Calculate Y0.  First find  RESJ = pi/2 ln(x/xn) J0(x),
c  where xn is a zero of Y0.
c
      if ( ax .le. three ) then
        up = ( ax - xy01 / two56 ) - xy02
        xy = xy0
      else if ( ax .le. five5 ) then
        up = ( ax - xy11 / two56 ) - xy12
        xy = xy1
      else
        up = ( ax - xy21 / two56 ) - xy22
        xy = xy2
      end if

      down = ax + xy

      if ( abs ( up ) .lt. p17 * down ) then
        w = up / down
        wsq = w * w
        xnum = plg(1)
        xden = wsq + qlg(1)
        do i = 2, 4
          xnum = xnum * wsq + plg(i)
          xden = xden * wsq + qlg(i)
        end do
        resj = pi2 * result * w * xnum / xden
      else
        resj = pi2 * result * log ( ax / xy )
      end if
c
c  Now calculate Y0 for appropriate interval, preserving
c  accuracy near the zero of Y0.
c
      if ( ax .le. three ) then
        xnum = py0(6) * zsq + py0(1)
        xden = zsq + qy0(1)
        do i = 2, 5
          xnum = xnum * zsq + py0(i)
          xden = xden * zsq + qy0(i)
        end do
      else if ( ax .le. five5 ) then
        xnum = py1(7) * zsq + py1(1)
        xden = zsq + qy1(1)
        do i = 2, 6
          xnum = xnum * zsq + py1(i)
          xden = xden * zsq + qy1(i)
        end do
      else
        xnum = py2(8) * zsq + py2(1)
        xden = zsq + qy2(1)
        do i = 2, 7
          xnum = xnum * zsq + py2(i)
          xden = xden * zsq + qy2(i)
        end do
      end if

      result = resj + up * down * xnum / xden

      return
c
c  Calculate J0 or Y0 for 8.0 < |ARG|.
c
  800 continue

      z = eight / ax
      w = ax / twopi
      w = aint ( w ) + oneov8
      w = ( ax - w * twopi1 ) - w * twopi2
      zsq = z * z
      xnum = p0(5) * zsq + p0(6)
      xden = zsq + q0(5)
      up = p1(5) * zsq + p1(6)
      down = zsq + q1(5)

      do i = 1, 4
        xnum = xnum * zsq + p0(i)
        xden = xden * zsq + q0(i)
        up = up * zsq + p1(i)
        down = down * zsq + q1(i)
      end do

      r0 = xnum / xden
      r1 = up / down

      if ( jint .eq. 0 ) then
        result = sqrt ( pi2 / ax ) 
     &    * ( r0 * cos ( w ) - z * r1 * sin ( w ) )
      else
        result = sqrt ( pi2 / ax ) 
     &    * ( r0 * sin ( w ) + z * r1 * cos ( w ) )
      end if

      return
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
c      Gamma = limit ( M -> +oo ) ( Sum ( 1 .le. N .le. M ) 1 / N ) - Log ( M )
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    27 April 1999
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

      euler_constant = 0.577215664901532860606512090082402431042D+00

      return
      end
      function i4_huge ( )

c*********************************************************************72
c
cc I4_HUGE returns a "huge" I4.
c
c  Discussion:
c
c    On an IEEE 32 bit machine, I4_HUGE should be 2**31 - 1, and its
c    bit pattern should be
c
c     01111111111111111111111111111111
c
c    In this case, its numerical value is 2147483647.
c
c    Using the Dec/Compaq/HP Alpha FORTRAN compiler FORT, I could
c    use I4_HUGE() and HUGE interchangeably.
c
c    However, when using the G95, the values returned by HUGE were
c    not equal to 2147483647, apparently, and were causing severe
c    and obscure errors in my random number generator, which needs to
c    add I4_HUGE to the seed whenever the seed is negative.  So I
c    am backing away from invoking HUGE, whereas I4_HUGE is under
c    my control.
c
c    Explanation: because under G95 the default integer type is 64 bits!
c    So HUGE ( 1 ) = a very very huge integer indeed, whereas
c    I4_HUGE ( ) = the same old 32 bit big value.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    26 January 2007
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, integer I4_HUGE, a "huge" I4.
c
      implicit none

      integer i4
      integer i4_huge

      i4_huge = 2147483647

      return
      end
      subroutine i4_to_halton_number_sequence ( seed, base, n, r )

c*********************************************************************72
c
cc I4_TO_HALTON_NUMBER_SEQUENCE: next N elements of a scalar Halton sequence.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    26 February 2001
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    John Halton,
c    On the efficiency of certain quasi-random sequences of points
c    in evaluating multi-dimensional integrals,
c    Numerische Mathematik,
c    Volume 2, pages 84-90, 1960.
c
c  Parameters:
c
c    Input, integer SEED, the index of the desired element.
c    Only the absolute value of SEED is considered.
c    SEED = 0 is allowed, and returns R = 0.
c
c    Input, integer BASE, the Halton base, which should
c    be a prime number.  This routine only checks that BASE is greater
c    than 1.
c
c    Input, integer N, the number of elements desired.
c
c    Output, double precision R(N), the SEED-th through (SEED+N-1)-th
c    elements of the Halton sequence for base BASE.
c
      implicit none

      integer n

      integer base
      double precision base_inv
      integer digit
      integer i
      double precision r(n)
      integer seed
      integer seed2
c
c  Set SEED2 = ( SEED, SEED+1, SEED+2, ..., SEED+N-1 )
c
      if ( base .le. 1 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'I4_TO_HALTON_NUMBER_SEQUENCE - Fatal error!'
        write ( *, '(a)' ) '  The input base BASE is .le. 1!'
        write ( *, '(a,i8)' ) '  BASE = ', base
        stop
      end if

      do i = 1, n

        seed2 = i + abs ( seed ) - 1

        r(i) = 0.0D+00

        base_inv = 1.0D+00 / dble ( base )

10      continue

        if ( seed2 .ne. 0 ) then
          digit = mod ( seed2, base )
          r(i) = r(i) + dble ( digit ) * base_inv
          base_inv = base_inv / dble ( base )
          seed2 = seed2 / base
          go to 10
        end if

      end do

      return
      end
      subroutine i4vec_indicator ( n, a )

c*********************************************************************72
c
cc I4VEC_INDICATOR sets an I4VEC to the indicator vector A(I)=I.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    09 November 2000
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
      subroutine p00_even ( prob, int_num, result )

c*********************************************************************72
c
cc P00_EVEN uses evenly spaced points to integrate a function.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    21 April 2003
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer PROB, the problem index.
c
c    Input, integer INT_NUM, the number of sample points.
c
c    Output, double precision RESULT, the approximate integral.
c
      implicit none

      integer int_num

      double precision a
      double precision b
      double precision fx(int_num)
      integer i
      integer prob
      double precision r8vec_sum
      double precision result
      double precision x(int_num)

      call p00_lim ( prob, a, b )

      if ( int_num .eq. 1 ) then
        x(1) = ( a + b ) / 2.0D+00
      else
        do i = 1, int_num
          x(i) = ( dble ( int_num - i     ) * a    
     &           + dble (           i - 1 ) * b  ) 
     &           / dble ( int_num     - 1 )
        end do
      end if

      call p00_fun ( prob, int_num, x, fx )

      result = ( b - a ) * r8vec_sum ( int_num, fx ) / dble ( int_num )

      return
      end
      subroutine p00_exact ( prob, exact )

c*********************************************************************72
c
cc P00_EXACT returns the exact integral for any problem.
c
c  Discussion:
c
c    This routine provides a "generic" interface to the exact integral
c    routines for the various problems, and allows a problem to be called
c    by number (PROB) rather than by name.
c
c    In some cases, the "exact" value of the integral is in fact
c    merely a respectable approximation.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    20 April 2003
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer PROB, the number of the desired test problem.
c
c    Output, double precision EXACT, the exact value of the integral.
c
      implicit none

      double precision exact
      integer prob

      if ( prob .eq. 1 ) then
        call p01_exact ( exact )
      else if ( prob .eq. 2 ) then
        call p02_exact ( exact )
      else if ( prob .eq. 3 ) then
        call p03_exact ( exact )
      else if ( prob .eq. 4 ) then
        call p04_exact ( exact )
      else if ( prob .eq. 5 ) then
        call p05_exact ( exact )
      else if ( prob .eq. 6 ) then
        call p06_exact ( exact )
      else if ( prob .eq. 7 ) then
        call p07_exact ( exact )
      else if ( prob .eq. 8 ) then
        call p08_exact ( exact )
      else if ( prob .eq. 9 ) then
        call p09_exact ( exact )
      else if ( prob .eq. 10 ) then
        call p10_exact ( exact )
      else if ( prob .eq. 11 ) then
        call p11_exact ( exact )
      else if ( prob .eq. 12 ) then
        call p12_exact ( exact )
      else if ( prob .eq. 13 ) then
        call p13_exact ( exact )
      else if ( prob .eq. 14 ) then
        call p14_exact ( exact )
      else if ( prob .eq. 15 ) then
        call p15_exact ( exact )
      else if ( prob .eq. 16 ) then
        call p16_exact ( exact )
      else if ( prob .eq. 17 ) then
        call p17_exact ( exact )
      else if ( prob .eq. 18 ) then
        call p18_exact ( exact )
      else if ( prob .eq. 19 ) then
        call p19_exact ( exact )
      else if ( prob .eq. 20 ) then
        call p20_exact ( exact )
      else if ( prob .eq. 21 ) then
        call p21_exact ( exact )
      else if ( prob .eq. 22 ) then
        call p22_exact ( exact )
      else if ( prob .eq. 23 ) then
        call p23_exact ( exact )
      else if ( prob .eq. 24 ) then
        call p24_exact ( exact )
      else if ( prob .eq. 25 ) then
        call p25_exact ( exact )
      else if ( prob .eq. 26 ) then
        call p26_exact ( exact )
      else if ( prob .eq. 27 ) then
        call p27_exact ( exact )
      else if ( prob .eq. 28 ) then
        call p28_exact ( exact )
      else if ( prob .eq. 29 ) then
        call p29_exact ( exact )
      else if ( prob .eq. 30 ) then
        call p30_exact ( exact )
      else if ( prob .eq. 31 ) then
        call p31_exact ( exact )
      else if ( prob .eq. 32 ) then
        call p32_exact ( exact )
      else if ( prob .eq. 33 ) then
        call p33_exact ( exact )
      else if ( prob .eq. 34 ) then
        call p34_exact ( exact )
      else if ( prob .eq. 35 ) then
        call p35_exact ( exact )
      else if ( prob .eq. 36 ) then
        call p36_exact ( exact )
      else if ( prob .eq. 37 ) then
        call p37_exact ( exact )
      else if ( prob .eq. 38 ) then
        call p38_exact ( exact )
      else if ( prob .eq. 39 ) then
        call p39_exact ( exact )
      else if ( prob .eq. 40 ) then
        call p40_exact ( exact )
      else if ( prob .eq. 41 ) then
        call p41_exact ( exact )
      else if ( prob .eq. 42 ) then
        call p42_exact ( exact )
      else if ( prob .eq. 43 ) then
        call p43_exact ( exact )
      else if ( prob .eq. 44 ) then
        call p44_exact ( exact )
      else if ( prob .eq. 45 ) then
        call p45_exact ( exact )
      else if ( prob .eq. 46 ) then
        call p46_exact ( exact )
      else if ( prob .eq. 47 ) then
        call p47_exact ( exact )
      else if ( prob .eq. 48 ) then
        call p48_exact ( exact )
      else if ( prob .eq. 49 ) then
        call p49_exact ( exact )
      else if ( prob .eq. 50 ) then
        call p50_exact ( exact )
      else if ( prob .eq. 51 ) then
        call p51_exact ( exact )
      else if ( prob .eq. 52 ) then
        call p52_exact ( exact )
      else if ( prob .eq. 53 ) then
        call p53_exact ( exact )
      else if ( prob .eq. 54 ) then
        call p54_exact ( exact )
      else if ( prob .eq. 55 ) then
        call p55_exact ( exact )
      else if ( prob .eq. 56 ) then
        call p56_exact ( exact )
      else if ( prob .eq. 57 ) then
        call p57_exact ( exact )
      else
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'P00_EXACT - Fatal error!'
        write ( *, '(a,i8)' ) '  Illegal problem number = ', prob
        stop
      end if

      return
      end
      subroutine p00_fun ( prob, n, x, fx )

c*********************************************************************72
c
cc P00_FUN evaluates the integrand for any problem.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    20 April 2003
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer PROB, the number of the desired test problem.
c
c    Input, integer N, the number of evaluation points.
c
c    Input, double precision X(N), the evaluation points.
c
c    Output, double precision FX(N), the integrand values.
c
      implicit none

      integer n

      double precision fx(n)
      integer prob
      double precision x(n)

      if ( prob .eq. 1 ) then
        call p01_fun ( n, x, fx )
      else if ( prob .eq. 2 ) then
        call p02_fun ( n, x, fx )
      else if ( prob .eq. 3 ) then
        call p03_fun ( n, x, fx )
      else if ( prob .eq. 4 ) then
        call p04_fun ( n, x, fx )
      else if ( prob .eq. 5 ) then
        call p05_fun ( n, x, fx )
      else if ( prob .eq. 6 ) then
        call p06_fun ( n, x, fx )
      else if ( prob .eq. 7 ) then
        call p07_fun ( n, x, fx )
      else if ( prob .eq. 8 ) then
        call p08_fun ( n, x, fx )
      else if ( prob .eq. 9 ) then
        call p09_fun ( n, x, fx )
      else if ( prob .eq. 10 ) then
        call p10_fun ( n, x, fx )
      else if ( prob .eq. 11 ) then
        call p11_fun ( n, x, fx )
      else if ( prob .eq. 12 ) then
        call p12_fun ( n, x, fx )
      else if ( prob .eq. 13 ) then
        call p13_fun ( n, x, fx )
      else if ( prob .eq. 14 ) then
        call p14_fun ( n, x, fx )
      else if ( prob .eq. 15 ) then
        call p15_fun ( n, x, fx )
      else if ( prob .eq. 16 ) then
        call p16_fun ( n, x, fx )
      else if ( prob .eq. 17 ) then
        call p17_fun ( n, x, fx )
      else if ( prob .eq. 18 ) then
        call p18_fun ( n, x, fx )
      else if ( prob .eq. 19 ) then
        call p19_fun ( n, x, fx )
      else if ( prob .eq. 20 ) then
        call p20_fun ( n, x, fx )
      else if ( prob .eq. 21 ) then
        call p21_fun ( n, x, fx )
      else if ( prob .eq. 22 ) then
        call p22_fun ( n, x, fx )
      else if ( prob .eq. 23 ) then
        call p23_fun ( n, x, fx )
      else if ( prob .eq. 24 ) then
        call p24_fun ( n, x, fx )
      else if ( prob .eq. 25 ) then
        call p25_fun ( n, x, fx )
      else if ( prob .eq. 26 ) then
        call p26_fun ( n, x, fx )
      else if ( prob .eq. 27 ) then
        call p27_fun ( n, x, fx )
      else if ( prob .eq. 28 ) then
        call p28_fun ( n, x, fx )
      else if ( prob .eq. 29 ) then
        call p29_fun ( n, x, fx )
      else if ( prob .eq. 30 ) then
        call p30_fun ( n, x, fx )
      else if ( prob .eq. 31 ) then
        call p31_fun ( n, x, fx )
      else if ( prob .eq. 32 ) then
        call p32_fun ( n, x, fx )
      else if ( prob .eq. 33 ) then
        call p33_fun ( n, x, fx )
      else if ( prob .eq. 34 ) then
        call p34_fun ( n, x, fx )
      else if ( prob .eq. 35 ) then
        call p35_fun ( n, x, fx )
      else if ( prob .eq. 36 ) then
        call p36_fun ( n, x, fx )
      else if ( prob .eq. 37 ) then
        call p37_fun ( n, x, fx )
      else if ( prob .eq. 38 ) then
        call p38_fun ( n, x, fx )
      else if ( prob .eq. 39 ) then
        call p39_fun ( n, x, fx )
      else if ( prob .eq. 40 ) then
        call p40_fun ( n, x, fx )
      else if ( prob .eq. 41 ) then
        call p41_fun ( n, x, fx )
      else if ( prob .eq. 42 ) then
        call p42_fun ( n, x, fx )
      else if ( prob .eq. 43 ) then
        call p43_fun ( n, x, fx )
      else if ( prob .eq. 44 ) then
        call p44_fun ( n, x, fx )
      else if ( prob .eq. 45 ) then
        call p45_fun ( n, x, fx )
      else if ( prob .eq. 46 ) then
        call p46_fun ( n, x, fx )
      else if ( prob .eq. 47 ) then
        call p47_fun ( n, x, fx )
      else if ( prob .eq. 48 ) then
        call p48_fun ( n, x, fx )
      else if ( prob .eq. 49 ) then
        call p49_fun ( n, x, fx )
      else if ( prob .eq. 50 ) then
        call p50_fun ( n, x, fx )
      else if ( prob .eq. 51 ) then
        call p51_fun ( n, x, fx )
      else if ( prob .eq. 52 ) then
        call p52_fun ( n, x, fx )
      else if ( prob .eq. 53 ) then
        call p53_fun ( n, x, fx )
      else if ( prob .eq. 54 ) then
        call p54_fun ( n, x, fx )
      else if ( prob .eq. 55 ) then
        call p55_fun ( n, x, fx )
      else if ( prob .eq. 56 ) then
        call p56_fun ( n, x, fx )
      else if ( prob .eq. 57 ) then
        call p57_fun ( n, x, fx )
      else
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'P00_FUN - Fatal error!'
        write ( *, '(a,i8)' ) '  Illegal problem number = ', prob
        stop
      end if

      return
      end
      subroutine p00_gauss_legendre ( prob, int_num, result )

c*********************************************************************72
c
cc P00_GAUSS_LEGENDRE applies a composite Gauss-Legendre rule.
c
c  Discussion:
c
c    A 4 point rule is used.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    12 November 2009
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer PROB, the problem index.
c
c    Input, integer INT_NUM, the number of subintervals.
c
c    Output, double precision RESULT, the approximate integral.
c
      implicit none

      integer gauss_num
      parameter ( gauss_num = 4 )

      double precision a
      double precision a_sub
      double precision b
      double precision b_sub
      double precision fx(gauss_num)
      double precision gauss_abs(gauss_num)
      double precision gauss_weight(gauss_num)
      double precision h
      integer i
      integer int_i
      integer int_num
      integer j
      integer prob
      double precision r8vec_dot_product
      double precision result
      double precision x(gauss_num)

      save gauss_abs
      save gauss_weight

      data gauss_abs / 
     &  -0.861136311594052575223946488893D+00, 
     &  -0.339981043584856264802665759103D+00, 
     &   0.339981043584856264802665759103D+00, 
     &   0.861136311594052575223946488893D+00 /
      data gauss_weight / 
     &  0.347854845137453857373063949222D+00, 
     &  0.652145154862546142626936050778D+00, 
     &  0.652145154862546142626936050778D+00, 
     &  0.347854845137453857373063949222D+00 /

      call p00_lim ( prob, a, b )

      h = ( b - a ) / dble ( int_num )

      result = 0.0D+00

      do int_i = 1, int_num

        a_sub = ( dble ( int_num - int_i + 1 ) * a   
     &          + dble (           int_i - 1 ) * b ) 
     &          / dble ( int_num             )

        b_sub = ( dble ( int_num - int_i ) * a   
     &          + dble (           int_i ) * b ) 
     &          / dble ( int_num         )

        do j = 1, gauss_num

          x(j) = ( ( 1.0D+00 - gauss_abs(j) ) * a_sub   
     &           + ( 1.0D+00 + gauss_abs(j) ) * b_sub ) 
     &            /  2.0D+00
        end do

        call p00_fun ( prob, gauss_num, x, fx )

        result = result + 0.5D+00 * h 
     &    * r8vec_dot_product ( gauss_num, gauss_weight, fx )

      end do

      return
      end
      subroutine p00_halton ( prob, int_num, result )

c*********************************************************************72
c
cc P00_HALTON applies a Halton sequence rule to integrate a function.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    19 December 2002
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    John Halton,
c    On the efficiency of certain quasi-random sequences of points
c    in evaluating multi-dimensional integrals,
c    Numerische Mathematik,
c    Volume 2, pages 84-90, 1960.
c
c  Parameters:
c
c    Input, integer PROB, the problem index.
c
c    Input, integer INT_NUM, the number of sample points.
c
c    Output, double precision RESULT, the approximate integral.
c
      implicit none

      integer int_num

      double precision a
      double precision b
      integer base
      double precision fx(int_num)
      integer i
      integer j
      integer prob
      double precision r8vec_sum
      double precision result
      integer seed
      double precision x(int_num)

      call p00_lim ( prob, a, b )

      seed = 1
      base = 2
      call i4_to_halton_number_sequence ( seed, base, int_num, x )

      do j = 1, int_num
        x(j) = a + ( b - a ) * x(j)
      end do

      call p00_fun ( prob, int_num, x, fx )

      result = ( b - a ) * r8vec_sum ( int_num, fx ) / dble ( int_num )

      return
      end
      subroutine p00_lim ( prob, a, b )

c*********************************************************************72
c
cc P00_LIM returns the integration limits for any problem.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    01 October 2002
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer PROB, the number of the desired test problem.
c
c    Output, double precision A, B, the limits of integration.
c
      implicit none

      double precision a
      double precision b
      integer prob

      if ( prob .eq. 1 ) then
        call p01_lim ( a, b )
      else if ( prob .eq. 2 ) then
        call p02_lim ( a, b )
      else if ( prob .eq. 3 ) then
        call p03_lim ( a, b )
      else if ( prob .eq. 4 ) then
        call p04_lim ( a, b )
      else if ( prob .eq. 5 ) then
        call p05_lim ( a, b )
      else if ( prob .eq. 6 ) then
        call p06_lim ( a, b )
      else if ( prob .eq. 7 ) then
        call p07_lim ( a, b )
      else if ( prob .eq. 8 ) then
        call p08_lim ( a, b )
      else if ( prob .eq. 9 ) then
        call p09_lim ( a, b )
      else if ( prob .eq. 10 ) then
        call p10_lim ( a, b )
      else if ( prob .eq. 11 ) then
        call p11_lim ( a, b )
      else if ( prob .eq. 12 ) then
        call p12_lim ( a, b )
      else if ( prob .eq. 13 ) then
        call p13_lim ( a, b )
      else if ( prob .eq. 14 ) then
        call p14_lim ( a, b )
      else if ( prob .eq. 15 ) then
        call p15_lim ( a, b )
      else if ( prob .eq. 16 ) then
        call p16_lim ( a, b )
      else if ( prob .eq. 17 ) then
        call p17_lim ( a, b )
      else if ( prob .eq. 18 ) then
        call p18_lim ( a, b )
      else if ( prob .eq. 19 ) then
        call p19_lim ( a, b )
      else if ( prob .eq. 20 ) then
        call p20_lim ( a, b )
      else if ( prob .eq. 21 ) then
        call p21_lim ( a, b )
      else if ( prob .eq. 22 ) then
        call p22_lim ( a, b )
      else if ( prob .eq. 23 ) then
        call p23_lim ( a, b )
      else if ( prob .eq. 24 ) then
        call p24_lim ( a, b )
      else if ( prob .eq. 25 ) then
        call p25_lim ( a, b )
      else if ( prob .eq. 26 ) then
        call p26_lim ( a, b )
      else if ( prob .eq. 27 ) then
        call p27_lim ( a, b )
      else if ( prob .eq. 28 ) then
        call p28_lim ( a, b )
      else if ( prob .eq. 29 ) then
        call p29_lim ( a, b )
      else if ( prob .eq. 30 ) then
        call p30_lim ( a, b )
      else if ( prob .eq. 31 ) then
        call p31_lim ( a, b )
      else if ( prob .eq. 32 ) then
        call p32_lim ( a, b )
      else if ( prob .eq. 33 ) then
        call p33_lim ( a, b )
      else if ( prob .eq. 34 ) then
        call p34_lim ( a, b )
      else if ( prob .eq. 35 ) then
        call p35_lim ( a, b )
      else if ( prob .eq. 36 ) then
        call p36_lim ( a, b )
      else if ( prob .eq. 37 ) then
        call p37_lim ( a, b )
      else if ( prob .eq. 38 ) then
        call p38_lim ( a, b )
      else if ( prob .eq. 39 ) then
        call p39_lim ( a, b )
      else if ( prob .eq. 40 ) then
        call p40_lim ( a, b )
      else if ( prob .eq. 41 ) then
        call p41_lim ( a, b )
      else if ( prob .eq. 42 ) then
        call p42_lim ( a, b )
      else if ( prob .eq. 43 ) then
        call p43_lim ( a, b )
      else if ( prob .eq. 44 ) then
        call p44_lim ( a, b )
      else if ( prob .eq. 45 ) then
        call p45_lim ( a, b )
      else if ( prob .eq. 46 ) then
        call p46_lim ( a, b )
      else if ( prob .eq. 47 ) then
        call p47_lim ( a, b )
      else if ( prob .eq. 48 ) then
        call p48_lim ( a, b )
      else if ( prob .eq. 49 ) then
        call p49_lim ( a, b )
      else if ( prob .eq. 50 ) then
        call p50_lim ( a, b )
      else if ( prob .eq. 51 ) then
        call p51_lim ( a, b )
      else if ( prob .eq. 52 ) then
        call p52_lim ( a, b )
      else if ( prob .eq. 53 ) then
        call p53_lim ( a, b )
      else if ( prob .eq. 54 ) then
        call p54_lim ( a, b )
      else if ( prob .eq. 55 ) then
        call p55_lim ( a, b )
      else if ( prob .eq. 56 ) then
        call p56_lim ( a, b )
      else if ( prob .eq. 57 ) then
        call p57_lim ( a, b )
      else
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'P00_LIM - Fatal error!'
        write ( *, '(a,i8)' ) '  Illegal problem number = ', prob
        stop
      end if

      return
      end
      subroutine p00_midpoint ( prob, int_num, result )

c*********************************************************************72
c
cc P00_MIDPOINT applies the composite midpoint rule to integrate a function.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    19 March 2002
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer PROB, the problem index.
c
c    Input, integer INT_NUM, the number of subintervals.
c
c    Output, double precision RESULT, the approximate integral.
c
      implicit none

      integer int_num

      double precision a
      double precision a_sub
      double precision b
      double precision b_sub
      double precision fx(int_num)
      integer int_i
      integer prob
      double precision r8vec_sum
      double precision result
      double precision x(int_num)

      call p00_lim ( prob, a, b )

      do int_i = 1, int_num

        a_sub = ( dble ( int_num - int_i + 1 ) * a   
     &          + dble (           int_i - 1 ) * b ) 
     &          / dble ( int_num             )

        b_sub = ( dble ( int_num - int_i ) * a   
     &          + dble (           int_i ) * b ) 
     &          / dble ( int_num         )

        x(int_i) = 0.5D+00 * ( a_sub + b_sub )

      end do

      call p00_fun ( prob, int_num, x, fx )

      result =  ( b - a ) * r8vec_sum ( int_num, fx ) 
     &  / dble ( int_num )

      return
      end
      subroutine p00_montecarlo ( prob, int_num, result )

c*********************************************************************72
c
cc P00_MONTECARLO applies the Monte Carlo rule to integrate a function.
c
c  Discussion:
c
c    This routine originally used an automatic array for X.  However,
c    under the G95 compiler, this was causing bizarre errors.  Replacing
c    the automatic array by an allocatable array made the problems
c    disappear.  Not an entirely satisfactory conclusion!
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    30 December 2011
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer PROB, the problem index.
c
c    Input, integer INT_NUM, the number of sample points.
c
c    Output, double precision RESULT, the approximate integral.
c
      implicit none

      integer int_num

      double precision a
      double precision b
      double precision fx(int_num)
      integer i
      integer prob
      double precision r8vec_sum
      double precision result
      integer seed
      double precision x(int_num)

      seed = 123456789

      call p00_lim ( prob, a, b )

      call r8vec_uniform ( int_num, a, b, seed, x )

      call p00_fun ( prob, int_num, x, fx )

      result = ( b - a ) * r8vec_sum ( int_num, fx ) 
     &  / dble ( int_num )

      return
      end
      subroutine p00_prob_num ( prob_num )

c*********************************************************************72
c
cc P00_PROB_NUM returns the number of test integration problems.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    20 April 2003
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, integer PROB_NUM, the number of test integration
c    problems.
c
      implicit none

      integer prob_num

      prob_num = 57

      return
      end
      subroutine p00_simpson ( prob, int_num, result )

c*********************************************************************72
c
cc P00_SIMPSON applies the composite Simpson rule to integrate a function.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    19 March 2002
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer PROB, the problem index.
c
c    Input, integer INT_NUM, the number of subintervals.
c
c    Output, double precision RESULT, the approximate integral.
c
      implicit none

      double precision a
      double precision a_sub
      double precision b
      double precision b_sub
      double precision fx1(1)
      double precision fx2(1)
      double precision fx3(1)
      double precision h
      integer int_i
      integer int_num
      integer prob
      double precision result
      double precision x1(1)
      double precision x2(1)
      double precision x3(1)

      call p00_lim ( prob, a, b )

      h = ( b - a ) / dble ( int_num )

      result = 0.0D+00

      do int_i = 1, int_num

        a_sub = ( dble ( int_num - int_i + 1 ) * a   
     &          + dble (           int_i - 1 ) * b ) 
     &          / dble ( int_num             )

        b_sub = ( dble ( int_num - int_i ) * a 
     &          + dble (           int_i ) * b ) 
     &          / dble ( int_num         )

        x1(1) = a_sub
        call p00_fun ( prob, 1, x1, fx1 )
        x2(1) = 0.5D+00 * ( a_sub + b_sub )
        call p00_fun ( prob, 1, x2, fx2 )
        x3(1) = b_sub
        call p00_fun ( prob, 1, x3, fx3 )

        result = result + h * ( 
     &                   fx1(1) 
     &       + 4.0D+00 * fx2(1) 
     &       +           fx3(1) ) / 6.0D+00

      end do

      return
      end
      subroutine p00_trapezoid ( prob, int_num, result )

c*********************************************************************72
c
cc P00_TRAPEZOID applies the composite trapezoid rule to integrate a function.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    01 October 2002
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer PROB, the problem index.
c
c    Input, integer INT_NUM, the number of subintervals.
c
c    Output, double precision RESULT, the approximate integral.
c
      implicit none

      double precision a
      double precision a_sub
      double precision b
      double precision b_sub
      double precision fx1(1)
      double precision fx2(1)
      double precision h
      integer int_i
      integer int_num
      integer prob
      double precision result
      double precision x1(1)
      double precision x2(1)

      call p00_lim ( prob, a, b )

      h = ( b - a ) / dble ( int_num )

      result = 0.0D+00

      do int_i = 1, int_num

        a_sub = ( dble ( int_num - int_i + 1 ) * a   
     &          + dble (           int_i - 1 ) * b ) 
     &          / dble ( int_num             )

        b_sub = ( dble ( int_num - int_i ) * a   
     &          + dble (           int_i ) * b ) 
     &          / dble ( int_num         )

        x1(1) = a_sub
        x2(1) = b_sub

        call p00_fun ( prob, 1, x1, fx1 )
        call p00_fun ( prob, 1, x2, fx2 )

        result = result + 0.5D+00 * h * ( fx1(1) + fx2(1) )

      end do

      return
      end
      subroutine p01_exact ( exact )

c*********************************************************************72
c
cc P01_EXACT returns the exact integral for problem 1.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    18 December 1998
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, double precision EXACT, the value of the integral.
c
      implicit none

      double precision exact

      exact = exp ( 1.0D+00 ) - 1.0D+00

      return
      end
      subroutine p01_fun ( n, x, fx )

c*********************************************************************72
c
cc P01_FUN evaluates the integrand for problem 1.
c
c  Interval:
c
c    0 .le. x .le. 1
c
c  Integrand:
c
c    exp ( x )
c
c  Antiderivative:
c
c    exp ( x )
c
c  Exact Integral:
c
c    exp ( 1 ) - 1
c
c  Approximate Integral (25 digits):
c
c    1.718281828459045235360287...
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    18 December 1998
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    David Kahaner,
c    Comparison of Numerical Quadrature Formulas,
c    in Mathematical Software, edited by John R Rice,
c    Academic Press, 1971.
c
c  Parameters:
c
c    Input, integer N, the number of evaluation points.
c
c    Input, double precision X(N), the evaluation points.
c
c    Output, double precision FX(N), the integrand values.
c
      implicit none

      integer n

      double precision fx(n)
      integer i
      double precision x(n)

      do i = 1, n
        fx(i) = exp ( x(i) )
      end do

      return
      end
      subroutine p01_lim ( a, b )

c*********************************************************************72
c
cc P01_LIM returns the integration limits for problem 1.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    22 December 1998
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, double precision A, B, the limits of integration.
c
      implicit none

      double precision a
      double precision b

      a = 0.0D+00
      b = 1.0D+00

      return
      end
      subroutine p02_exact ( exact )

c*********************************************************************72
c
cc P02_EXACT returns the exact integral for problem 2.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    18 December 1998
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, double precision EXACT, the value of the integral.
c
      implicit none

      double precision exact

      exact = 0.7D+00

      return
      end
      subroutine p02_fun ( n, x, fx )

c*********************************************************************72
c
cc P02_FUN evaluates the integrand for problem 2.
c
c  Discussion:
c
c    The integrand is discontinuous at X = 0.3.
c
c  Interval:
c
c    0 .le. x .le. 1
c
c  Integrand:
c
c    if ( x < 0.3 )
c      f(x) = 0
c    else
c      f(x) = 1
c
c  Antiderivative:
c
c    if ( x < 0.3 )
c      g(x) = 0
c    else
c      g(x) = X - 0.3
c
c  Exact Integral:
c
c    0.7
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    18 December 1998
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    David Kahaner,
c    Comparison of Numerical Quadrature Formulas,
c    in Mathematical Software, edited by John R Rice,
c    Academic Press, 1971.
c
c  Parameters:
c
c    Input, integer N, the number of evaluation points.
c
c    Input, double precision X(N), the evaluation points.
c
c    Output, double precision FX(N), the integrand values.
c
      implicit none

      integer n

      double precision fx(n)
      integer i
      double precision x(n)

      do i = 1, n

        if ( x(i) .lt. 0.3D+00 ) then
          fx(i) = 0.0D+00
        else
          fx(i) = 1.0D+00
        end if

      end do

      return
      end
      subroutine p02_lim ( a, b )

c*********************************************************************72
c
cc P02_LIM returns the integration limits for problem 2.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    22 December 1998
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, double precision A, B, the limits of integration.
c
      implicit none

      double precision a
      double precision b

      a = 0.0D+00
      b = 1.0D+00

      return
      end
      subroutine p03_exact ( exact )

c*********************************************************************72
c
cc P03_EXACT returns the exact integral for problem 3.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    18 December 1998
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, double precision EXACT, the value of the integral.
c
      implicit none

      double precision exact

      exact = 2.0D+00 / 3.0D+00

      return
      end
      subroutine p03_fun ( n, x, fx )

c*********************************************************************72
c
cc P03_FUN evaluates the integrand for problem 3.
c
c  Discussion:
c
c    The integrand is not differentiable at X = 0.
c
c  Interval:
c
c    0 .le. x .le. 1
c
c  Integrand:
c
c    sqrt ( x )
c
c  Antiderivative:
c
c    ( 2 / 3 ) * x^(3/2)
c
c  Exact Integral:
c
c    2 / 3
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    31 October 2000
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    David Kahaner,
c    Comparison of Numerical Quadrature Formulas,
c    in Mathematical Software, edited by John R Rice,
c    Academic Press, 1971.
c
c  Parameters:
c
c    Input, integer N, the number of evaluation points.
c
c    Input, double precision X(N), the evaluation points.
c
c    Output, double precision FX(N), the integrand values.
c
      implicit none

      integer n

      double precision fx(n)
      integer i
      double precision x(n)

      do i = 1, n
        fx(i) = sqrt ( x(i) )
      end do

      return
      end
      subroutine p03_lim ( a, b )

c*********************************************************************72
c
cc P03_LIM returns the integration limits for problem 3.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    22 December 1998
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, double precision A, B, the limits of integration.
c
      implicit none

      double precision a
      double precision b

      a = 0.0D+00
      b = 1.0D+00

      return
      end
      subroutine p04_exact ( exact )

c*********************************************************************72
c
cc P04_EXACT returns the estimated integral for problem 4.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    18 December 1998
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, double precision EXACT, the estimated value of the integral.
c
      implicit none

      double precision exact

      exact = 0.47942822668880166736D+00

      return
      end
      subroutine p04_fun ( n, x, fx )

c*********************************************************************72
c
cc P04_FUN evaluates the integrand for problem 4.
c
c  Interval:
c
c    -1 .le. x .le. 1
c
c  Integrand:
c
c    0.92 * cosh ( x ) - cos ( x )
c
c  Antiderivative:
c
c    0.92 * sinh ( x ) - sin ( x )
c
c  Exact Integral:
c
c    1.84 * sinh ( 1 ) - 2 * sin ( 1 )
c
c  Approximate Integral (20 digits):
c
c    0.47942822668880166736...
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    31 October 2000
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Charles Clenshaw, Alan Curtis,
c    A Method for Numerical Integration on an Automatic Computer,
c    Numerische Mathematik,
c    Volume 2, Number 1, December 1960, pages 197-205.
c
c    David Kahaner,
c    Comparison of Numerical Quadrature Formulas,
c    in Mathematical Software, edited by John R Rice,
c    Academic Press, 1971.
c
c  Parameters:
c
c    Input, integer N, the number of evaluation points.
c
c    Input, double precision X(N), the evaluation points.
c
c    Output, double precision FX(N), the integrand values.
c
      implicit none

      integer n

      double precision fx(n)
      integer i
      double precision x(n)

      do i = 1, n
        fx(i) = 0.92D+00 * cosh ( x(i) ) - cos ( x(i) )
      end do

      return
      end
      subroutine p04_lim ( a, b )

c*********************************************************************72
c
cc P04_LIM returns the integration limits for problem 4.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    22 December 1998
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, double precision A, B, the limits of integration.
c
      implicit none

      double precision a
      double precision b

      a = -1.0D+00
      b = 1.0D+00

      return
      end
      subroutine p05_exact ( exact )

c*********************************************************************72
c
cc P05_EXACT returns the estimated integral for problem 5.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    06 October 2006
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, double precision EXACT, the estimated value of the integral.
c
      implicit none

      double precision exact

      exact = 1.5822329637296729331D+00

      return
      end
      subroutine p05_fun ( n, x, fx )

c*********************************************************************72
c
cc P05_FUN evaluates the integrand for problem 5.
c
c  Interval:
c
c    -1 .le. x .le. 1
c
c  Integrand:
c
c    1 / ( x^4 + x^2 + 0.9 )
c
c  Approximate Integral (20 digits):
c
c    1.5822329637296729331...
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    06 October 2006
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Charles Clenshaw, Alan Curtis,
c    A Method for Numerical Integration on an Automatic Computer,
c    Numerische Mathematik,
c    Volume 2, Number 1, December 1960, pages 197-205.
c
c  Parameters:
c
c    Input, integer N, the number of evaluation points.
c
c    Input, double precision X(N), the evaluation points.
c
c    Output, double precision FX(N), the integrand values.
c
      implicit none

      integer n

      double precision fx(n)
      integer i
      double precision x(n)

      do i = 1, n 
        fx(i) = 1.0D+00 / ( x(i)**4 + x(i)**2 + 0.9D+00 )
      end do

      return
      end
      subroutine p05_lim ( a, b )

c*********************************************************************72
c
cc P05_LIM returns the integration limits for problem 5.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    06 October 2006
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, double precision A, B, the limits of integration.
c
      implicit none

      double precision a
      double precision b

      a = -1.0D+00
      b = 1.0D+00

      return
      end
      subroutine p06_exact ( exact )

c*********************************************************************72
c
cc P06_EXACT returns the exact integral for problem 6.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    06 October 2006
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, double precision EXACT, the estimated value of the integral.
c
      implicit none

      double precision exact

      exact = 1.460447131787105D+00

      return
      end
      subroutine p06_fun ( n, x, fx )

c*********************************************************************72
c
cc P06_FUN evaluates the integrand for problem 6.
c
c  Interval:
c
c    -1 .le. x .le. 1
c
c  Integrand:
c
c    sqrt ( abs ( x + 0.5 ) )
c
c  Exact Integral:
c
c    ( sqrt ( 2 ) + 3 * sqrt ( 6 ) ) / 6 = 1.460447131787105
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    06 October 2006
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Charles Clenshaw, Alan Curtis,
c    A Method for Numerical Integration on an Automatic Computer,
c    Numerische Mathematik,
c    Volume 2, Number 1, December 1960, pages 197-205.
c
c  Parameters:
c
c    Input, integer N, the number of evaluation points.
c
c    Input, double precision X(N), the evaluation points.
c
c    Output, double precision FX(N), the integrand values.
c
      implicit none

      integer n

      double precision fx(n)
      integer i
      double precision x(n)

      do i = 1, n
        fx(i) = sqrt ( abs ( x(i) + 0.5D+00 ) )
      end do

      return
      end
      subroutine p06_lim ( a, b )

c*********************************************************************72
c
cc P06_LIM returns the integration limits for problem 6.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    06 October 2006
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, double precision A, B, the limits of integration.
c
      implicit none

      double precision a
      double precision b

      a = -1.0D+00
      b = 1.0D+00

      return
      end
      subroutine p07_exact ( exact )

c*********************************************************************72
c
cc P07_EXACT returns the exact integral for problem 7.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    18 December 1998
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, double precision EXACT, the value of the integral.
c
      implicit none

      double precision exact

      exact = 2.0D+00

      return
      end
      subroutine p07_fun ( n, x, fx )

c*********************************************************************72
c
cc P07_FUN evaluates the integrand for problem 7.
c
c  Discussion:
c
c    The integrand is singular at X = 0.
c
c  Interval:
c
c    0 .le. X .le. 1
c
c  Integrand:
c
c    1 / sqrt ( X )
c
c  Antiderivative:
c
c    2 * sqrt ( X )
c
c  Exact Integral:
c
c    2
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    31 October 2000
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    David Kahaner,
c    Comparison of Numerical Quadrature Formulas,
c    in Mathematical Software, edited by John R Rice,
c    Academic Press, 1971
c
c  Parameters:
c
c    Input, integer N, the number of evaluation points.
c
c    Input, double precision X(N), the evaluation points.
c
c    Output, double precision FX(N), the integrand values.
c
      implicit none

      integer n

      double precision fx(n)
      integer i
      double precision x(n)

      do i = 1, n

        if ( 0.0D+00 .lt. x(i) ) then
          fx(i) = 1.0D+00 / sqrt ( x(i) )
        else
          fx(i) = 0.0D+00
        end if

      end do

      return
      end
      subroutine p07_lim ( a, b )

c*********************************************************************72
c
cc P07_LIM returns the integration limits for problem 7.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    22 December 1998
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, double precision A, B, the limits of integration.
c
      implicit none

      double precision a
      double precision b

      a = 0.0D+00
      b = 1.0D+00

      return
      end
      subroutine p08_exact ( exact )

c*********************************************************************72
c
cc P08_EXACT returns the estimated integral for problem 8.
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
c    Output, double precision EXACT, the estimated value of the integral.
c
      implicit none

      double precision exact

      exact = 0.86697298733991103757D+00

      return
      end
      subroutine p08_fun ( n, x, fx )

c*********************************************************************72
c
cc P08_FUN evaluates the integrand for problem 8.
c
c  Interval:
c
c    0 .le. x .le. 1
c
c  Integrand:
c
c    1 / ( 1 + X^4 )
c
c  Antiderivative:
c
c    (1/8) * sqrt ( 2 )
c    * ln ( ( X^2 + sqrt ( 2 ) * X + 1 ) / ( X^2 - sqrt ( 2 ) * X + 1 ) )
c    + (1/4) * sqrt ( 2 ) * arctan ( sqrt ( 2 ) * X / ( 1 - X^2 ) )
c
c  Approximate Integral (20 digits):
c
c    0.86697298733991103757...
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    31 October 2000
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    David Kahaner,
c    Comparison of Numerical Quadrature Formulas,
c    in Mathematical Software, edited by John R Rice,
c    Academic Press, 1971.
c
c  Parameters:
c
c    Input, integer N, the number of evaluation points.
c
c    Input, double precision X(N), the evaluation points.
c
c    Output, double precision FX(N), the integrand values.
c
      implicit none

      integer n

      double precision fx(n)
      integer i
      double precision x(n)

      do i = 1, n
        fx(i) = 1.0D+00 / ( 1.0D+00 + x(i)**4 )
      end do

      return
      end
      subroutine p08_lim ( a, b )

c*********************************************************************72
c
cc P08_LIM returns the integration limits for problem 8.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    22 December 1998
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, double precision A, B, the limits of integration.
c
      implicit none

      double precision a
      double precision b

      a = 0.0D+00
      b = 1.0D+00

      return
      end
      subroutine p09_exact ( exact )

c*********************************************************************72
c
cc P09_EXACT returns the estimated integral for problem 9.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    18 December 1998
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, double precision EXACT, the estimated value of the integral.
c
      implicit none

      double precision exact

      exact = 1.1547005383792515290D+00

      return
      end
      subroutine p09_fun ( n, x, fx )

c*********************************************************************72
c
cc P09_FUN evaluates the integrand for problem 9.
c
c  Discussion:
c
c    The integrand is oscillatory, going through 5 periods in [0,1].
c
c  Interval:
c
c    0 .le. x .le. 1
c
c  Integrand:
c
c    2 / ( 2 + sin ( 10 * pi * X ) )
c
c  Antiderivative:
c
c    1 / ( 5 * pi * sqrt ( 3 ) ) *
c    arctan ( ( 1 + 2 * tan ( 5 * pi * X ) ) / sqrt ( 3 ) )
c
c  Exact Integral:
c
c    2 / sqrt ( 3 )
c
c  Approximate Integral (20 digits):
c
c    1.1547005383792515290...
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    29 October 2000
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    David Kahaner,
c    Comparison of Numerical Quadrature Formulas,
c    in Mathematical Software, edited by John R Rice,
c    Academic Press, 1971.
c
c  Parameters:
c
c    Input, integer N, the number of evaluation points.
c
c    Input, double precision X(N), the evaluation points.
c
c    Output, double precision FX(N), the integrand values.
c
      implicit none

      integer n

      double precision fx(n)
      integer i
      double precision x(n)
      double precision pi
      parameter ( pi = 3.141592653589793D+00 )

      do i = 1, n
        fx(i) = 2.0D+00 / ( 2.0D+00 + sin ( 10.0D+00 * pi * x(i) ) )
      end do

      return
      end
      subroutine p09_lim ( a, b )

c*********************************************************************72
c
cc P09_LIM returns the integration limits for problem 9.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    22 December 1998
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, double precision A, B, the limits of integration.
c
      implicit none

      double precision a
      double precision b

      a = 0.0D+00
      b = 1.0D+00

      return
      end
      subroutine p10_exact ( exact )

c*********************************************************************72
c
cc P10_EXACT returns the estimated integral for problem 10.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    18 December 1998
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, double precision EXACT, the estimated value of the integral.
c
      implicit none

      double precision exact

      exact = 0.6931471805599453094172321D+00

      return
      end
      subroutine p10_fun ( n, x, fx )

c*********************************************************************72
c
cc P10_FUN evaluates the integrand for problem 10.
c
c  Interval:
c
c    0 .le. x .le. 1
c
c  Integrand:
c
c    1 / ( 1 + X )
c
c  Antiderivative:
c
c    ln ( 1 + X )
c
c  Exact Integral:
c
c    ln ( 2 )
c
c  Approximate Integral (25 digits):
c
c    0.6931471805599453094172321...
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    31 October 2000
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    David Kahaner,
c    Comparison of Numerical Quadrature Formulas,
c    in Mathematical Software, edited by John R Rice,
c    Academic Press, 1971.
c
c  Parameters:
c
c    Input, integer N, the number of evaluation points.
c
c    Input, double precision X(N), the evaluation points.
c
c    Output, double precision FX(N), the integrand values.
c
      implicit none

      integer n

      double precision fx(n)
      integer i
      double precision x(n)

      do i = 1, n
        fx(i) = 1.0D+00 / ( 1.0D+00 + x(i) )
      end do

      return
      end
      subroutine p10_lim ( a, b )

c*********************************************************************72
c
cc P10_LIM returns the integration limits for problem 10.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    22 December 1998
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, double precision A, B, the limits of integration.
c
      implicit none

      double precision a
      double precision b

      a = 0.0D+00
      b = 1.0D+00

      return
      end
      subroutine p11_exact ( exact )

c*********************************************************************72
c
cc P11_EXACT returns the estimated integral for problem 11.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    18 December 1998
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, double precision EXACT, the estimated value of the integral.
c
      implicit none

      double precision exact

      exact = 0.37988549304172247537D+00

      return
      end
      subroutine p11_fun ( n, x, fx )

c*********************************************************************72
c
cc P11_FUN evaluates the integrand for problem 11.
c
c  Interval:
c
c    0 .le. x .le. 1
c
c  Integrand:
c
c    1 / ( 1 + exp ( X ) )
c
c  Antiderivative:
c
c    ln ( exp ( X ) / ( 1 + exp ( X ) ) )
c
c  Exact Integral:
c
c    ln ( 2 * E / ( 1 + E ) )
c
c  Approximate Integral (20 digits):
c
c    0.37988549304172247537...
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    18 December 1998
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    David Kahaner,
c    Comparison of Numerical Quadrature Formulas,
c    in Mathematical Software, edited by John R Rice,
c    Academic Press, 1971.
c
c  Parameters:
c
c    Input, integer N, the number of evaluation points.
c
c    Input, double precision X(N), the evaluation points.
c
c    Output, double precision FX(N), the integrand values.
c
      implicit none

      integer n

      double precision fx(n)
      integer i
      double precision x(n)

      do i = 1, n
        fx(i) = 1.0D+00 / ( 1.0D+00 + exp ( x(i) ) )
      end do

      return
      end
      subroutine p11_lim ( a, b )

c*********************************************************************72
c
cc P11_LIM returns the integration limits for problem 11.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    22 December 1998
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, double precision A, B, the limits of integration.
c
      implicit none

      double precision a
      double precision b

      a = 0.0D+00
      b = 1.0D+00

      return
      end
      subroutine p12_exact ( exact )

c*********************************************************************72
c
cc P12_EXACT returns the estimated integral for problem 12.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    18 December 1998
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, double precision EXACT, the estimated value of the integral.
c
      implicit none

      double precision exact

      exact = 0.77750463411224827642D+00

      return
      end
      subroutine p12_fun ( n, x, fx )

c*********************************************************************72
c
cc P12_FUN evaluates the integrand for problem 12.
c
c  Discussion:
c
c    The integrand has a removable singularity at X = 0.
c
c  Interval:
c
c    0 .le. x .le. 1
c
c  Integrand:
c
c    X / ( exp ( X ) - 1 )
c
c  Antiderivative:
c
c    The Debye function.
c
c  Approximate Integral (20 digits):
c
c    0.77750463411224827642...
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    18 December 1998
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    David Kahaner,
c    Comparison of Numerical Quadrature Formulas,
c    in Mathematical Software, edited by John R Rice,
c    Academic Press, 1971.
c
c  Parameters:
c
c    Input, integer N, the number of evaluation points.
c
c    Input, double precision X(N), the evaluation points.
c
c    Output, double precision FX(N), the integrand values.
c
      implicit none

      integer n

      double precision fx(n)
      integer i
      double precision x(n)

      do i = 1, n

        if ( x(i) .eq. 0.0D+00 ) then
          fx(i) = 1.0D+00
        else
          fx(i) = x(i) / ( exp ( x(i) ) - 1.0D+00 )
        end if

      end do

      return
      end
      subroutine p12_lim ( a, b )

c*********************************************************************72
c
cc P12_LIM returns the integration limits for problem 12.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    22 December 1998
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, double precision A, B, the limits of integration.
c
      implicit none

      double precision a
      double precision b

      a = 0.0D+00
      b = 1.0D+00

      return
      end
      subroutine p13_exact ( exact )

c*********************************************************************72
c
cc P13_EXACT returns the estimated integral for problem 13.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    06 November 2000
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, double precision EXACT, the estimated value of the integral.
c
      implicit none

      double precision a
      double precision b
      double precision exact
      double precision r8_si

      call p13_lim ( a, b )

      exact = r8_si ( b ) - r8_si ( a )

      return
      end
      subroutine p13_fun ( n, x, fx )

c*********************************************************************72
c
cc P13_FUN evaluates the integrand for problem 13.
c
c  Interval:
c
c    0 .le. x .le. 10
c
c  Integrand:
c
c    sin ( X ) / X
c
c  Approximate Integral (20 digits):
c
c    1.6583475942188740493...
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    06 November 2000
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the number of evaluation points.
c
c    Input, double precision X(N), the evaluation points.
c
c    Output, double precision FX(N), the integrand values.
c
      implicit none

      integer n

      double precision fx(n)
      integer i
      double precision x(n)

      do i = 1, n

        if ( x(i) .eq. 0.0D+00 ) then
          fx(i) = 1.0D+00
        else
          fx(i) = sin ( x(i) ) / x(i)
        end if

      end do

      return
      end
      subroutine p13_lim ( a, b )

c*********************************************************************72
c
cc P13_LIM returns the integration limits for problem 13.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    06 November 2000
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, double precision A, B, the limits of integration.
c
      implicit none

      double precision a
      double precision b

      a = 0.0D+00
      b = 10.0D+00

      return
      end
      subroutine p14_exact ( exact )

c*********************************************************************72
c
cc P14_EXACT returns the estimated integral for problem 14.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    18 December 1998
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, double precision EXACT, the estimated value of the integral.
c
      implicit none

      double precision exact

      exact = 0.500000211166D+00

      return
      end
      subroutine p14_fun ( n, x, fx )

c*********************************************************************72
c
cc P14_FUN evaluates the integrand for problem 14.
c
c  Discussion:
c
c    For X's that aren't actually very big, the function becomes very
c    small.  Some compilers may product code that fails in these cases.
c    An attempt has been made to return a value of 0 when the computed
c    value of F(X) would be extremely small.
c
c  Interval:
c
c    0 .le. x .le. 10
c
c  Integrand:
c
c    sqrt ( 50 ) * exp ( - 50 * pi * x * x )
c
c  Exact Integral:
c
c    0.5 * erf ( 50 * sqrt ( 2 * pi ) )
c
c  Approximate Integral:
c
c    0.500000211166...
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    17 May 2007
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    David Kahaner,
c    Comparison of Numerical Quadrature Formulas,
c    in Mathematical Software, edited by John R Rice,
c    Academic Press, 1971.
c
c  Parameters:
c
c    Input, integer N, the number of evaluation points.
c
c    Input, double precision X(N), the evaluation points.
c
c    Output, double precision FX(N), the integrand values.
c
      implicit none

      integer n

      double precision fx(n)
      integer i
      double precision x(n)
      double precision pi
      parameter ( pi = 3.141592653589793D+00 )
      double precision r8_epsilon
      double precision x_max

      save x_max

      data x_max / 0.0D+00 /

      if ( x_max .eq. 0.0D+00 ) then
        x_max = sqrt ( log ( max ( r8_epsilon ( ), 1.0D-10 ) ) 
     &   / ( - 50.0D+00 * pi ) )
      end if

      do i = 1, n

        if ( x_max .lt. abs ( x(i) ) ) then
          fx(i) = 0.0D+00
        else
          fx(i) = sqrt ( 50.0D+00 ) 
     &      * exp ( - 50.0D+00 * pi * x(i) * x(i) )
        end if

      end do

      return
      end
      subroutine p14_lim ( a, b )

c*********************************************************************72
c
cc P14_LIM returns the integration limits for problem 14.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    22 December 1998
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, double precision A, B, the limits of integration.
c
      implicit none

      double precision a
      double precision b

      a = 0.0D+00
      b = 10.0D+00

      return
      end
      subroutine p15_exact ( exact )

c*********************************************************************72
c
cc P15_EXACT returns the exact integral for problem 15.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    18 December 1998
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, double precision EXACT, the value of the integral.
c
      implicit none

      double precision exact

      exact = 1.0D+00

      return
      end
      subroutine p15_fun ( n, x, fx )

c*********************************************************************72
c
cc P15_FUN evaluates the integrand for problem 15.
c
c  Interval:
c
c    0 .le. x .le. 10
c
c  Integrand:
c
c    25 * exp ( - 25 * X )
c
c  Antiderivative:
c
c    - exp ( - 25 * X )
c
c  Exact Integral:
c
c    1 - exp ( - 250 )
c
c  Approximate Integral:
c
c    1.00000000...
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    18 December 1998
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    David Kahaner,
c    Comparison of Numerical Quadrature Formulas,
c    in Mathematical Software, edited by John R Rice,
c    Academic Press, 1971.
c
c  Parameters:
c
c    Input, integer N, the number of evaluation points.
c
c    Input, double precision X(N), the evaluation points.
c
c    Output, double precision FX(N), the integrand values.
c
      implicit none

      integer n

      double precision fx(n)
      integer i
      double precision x(n)

      do i = 1, n
        fx(i) = 25.0D+00 * exp ( - 25.0D+00 * x(i) )
      end do

      return
      end
      subroutine p15_lim ( a, b )

c*********************************************************************72
c
cc P15_LIM returns the integration limits for problem 15.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    22 December 1998
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, double precision A, B, the limits of integration.
c
      implicit none

      double precision a
      double precision b

      a = 0.0D+00
      b = 10.0D+00

      return
      end
      subroutine p16_exact ( exact )

c*********************************************************************72
c
cc P16_EXACT returns the exact integral for problem 16.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    23 July 2008
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, double precision EXACT, the value of the integral.
c
      implicit none

      double precision exact

      exact = 0.49936338107645674464D+00

      return
      end
      subroutine p16_fun ( n, x, fx )

c*********************************************************************72
c
cc P16_FUN evaluates the integrand for problem 16.
c
c  Interval:
c
c    0 .le. x .le. 10
c
c  Integrand:
c
c    50.0 / ( pi * ( 2500.0 * X * X + 1.0 ) )
c
c  Antiderivative:
c
c    ( 1 / pi ) * arctan ( 50 * X )
c
c  Approximate Integral (20 digits):
c
c    0.49936338107645674464...
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    23 July 2008
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    David Kahaner,
c    Comparison of Numerical Quadrature Formulas,
c    in Mathematical Software, edited by John R Rice,
c    Academic Press, 1971.
c
c  Parameters:
c
c    Input, integer N, the number of evaluation points.
c
c    Input, double precision X(N), the evaluation points.
c
c    Output, double precision FX(N), the integrand values.
c
      implicit none

      integer n

      double precision fx(n)
      integer i
      double precision pi
      parameter ( pi = 3.141592653589793D+00 )
      double precision x(n)

      do i = 1, n
        fx(i) = 50.0D+00 / pi / ( 2500.0D+00 * x(i) * x(i) + 1.0D+00 )
      end do

      return
      end
      subroutine p16_lim ( a, b )

c*********************************************************************72
c
cc P16_LIM returns the integration limits for problem 16.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    22 December 1998
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, double precision A, B, the limits of integration.
c
      implicit none

      double precision a
      double precision b

      a = 0.0D+00
      b = 1.0D+00

      return
      end
      subroutine p17_exact ( exact )

c*********************************************************************72
c
cc P17_EXACT returns the estimated integral for problem 17.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    09 May 2001
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, double precision EXACT, the estimated value of the integral.
c
      implicit none

      double precision exact

      exact = 0.5D+00

      return
      end
      subroutine p17_fun ( n, x, fx )

c*********************************************************************72
c
cc P17_FUN evaluates the integrand for problem 17.
c
c  Interval:
c
c    0 .le. x .le. 1
c
c  Integrand:
c
c    ( sin ( 50 * pi * X ) )^2
c
c  Antiderivative:
c
c    1/2 X - sin ( 100 * pi * X ) / ( 200 * pi )
c
c  Approximate Integral:
c
c    0.5
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    09 May 2001
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the number of evaluation points.
c
c    Input, double precision X(N), the evaluation points.
c
c    Output, double precision FX(N), the integrand values.
c
      implicit none

      integer n

      double precision fx(n)
      integer i
      double precision x(n)
      double precision pi
      parameter ( pi = 3.141592653589793D+00 )

      do i = 1, n
        fx(i) = ( sin ( 50.0D+00 * pi * x(i) ) )**2
      end do

      return
      end
      subroutine p17_lim ( a, b )

c*********************************************************************72
c
cc P17_LIM returns the integration limits for problem 17.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    22 December 1998
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, double precision A, B, the limits of integration.
c
      implicit none

      double precision a
      double precision b

      a = 0.0D+00
      b = 1.0D+00

      return
      end
      subroutine p18_exact ( exact )

c*********************************************************************72
c
cc P18_EXACT returns the estimated integral for problem 18.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    01 November 2000
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, double precision EXACT, the estimated value of the integral.
c
      implicit none

      double precision exact

      exact = 0.17055734950243820437D+00

      return
      end
      subroutine p18_fun ( n, x, fx )

c*********************************************************************72
c
cc P18_FUN evaluates the integrand for problem 18.
c
c  Interval:
c
c    0 .le. x .le. 1
c
c  Integrand:
c
c    X / ( exp ( X ) + 1 )
c
c  Approximate Integral (20 digits):
c
c    0.17055734950243820437...
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    01 November 2000
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Hermann Engels,
c    Numerical Quadrature and Cubature,
c    Academic Press, 1980.
c
c  Parameters:
c
c    Input, integer N, the number of evaluation points.
c
c    Input, double precision X(N), the evaluation points.
c
c    Output, double precision FX(N), the integrand values.
c
      implicit none

      integer n

      double precision fx(n)
      integer i
      double precision x(n)

      do i = 1, n
        fx(i) = x(i) / ( exp ( x(i) ) + 1.0D+00 )
      end do

      return
      end
      subroutine p18_lim ( a, b )

c*********************************************************************72
c
cc P18_LIM returns the integration limits for problem 18.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    01 November 2000
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, double precision A, B, the limits of integration.
c
      implicit none

      double precision a
      double precision b

      a = 0.0D+00
      b = 1.0D+00

      return
      end
      subroutine p19_exact ( exact )

c*********************************************************************72
c
cc P19_EXACT returns the exact integral for problem 19.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    18 December 1998
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, double precision EXACT, the value of the integral.
c
      implicit none

      double precision exact

      exact = - 1.0D+00

      return
      end
      subroutine p19_fun ( n, x, fx )

c*********************************************************************72
c
cc P19_FUN evaluates the integrand for problem 19.
c
c  Discussion:
c
c    The integrand is singular at X = 0.
c
c  Interval:
c
c    0 .le. x .le. 1
c
c  Integrand:
c
c    ln ( X )
c
c  Antiderivative:
c
c    X * ln ( X ) - X
c
c  Exact Integral:
c
c    -1
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    18 December 1998
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    David Kahaner,
c    Comparison of Numerical Quadrature Formulas,
c    in Mathematical Software, edited by John R Rice,
c    Academic Press, 1971.
c
c  Parameters:
c
c    Input, integer N, the number of evaluation points.
c
c    Input, double precision X(N), the evaluation points.
c
c    Output, double precision FX(N), the integrand values.
c
      implicit none

      integer n

      double precision fx(n)
      integer i
      double precision x(n)

      do i = 1, n

        if ( x(i) .le. 1.0D-15 ) then
          fx(i) = 0.0D+00
        else
          fx(i) = log ( x(i) )
        end if

      end do

      return
      end
      subroutine p19_lim ( a, b )

c*********************************************************************72
c
cc P19_LIM returns the integration limits for problem 19.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    22 December 1998
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, double precision A, B, the limits of integration.
c
      implicit none

      double precision a
      double precision b

      a = 0.0D+00
      b = 1.0D+00

      return
      end
      subroutine p20_exact ( exact )

c*********************************************************************72
c
cc P20_EXACT returns the estimated integral for problem 20.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    18 December 1998
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, double precision EXACT, the estimated value of the integral.
c
      double precision exact

      exact = 1.5643964440690497731D+00

      return
      end
      subroutine p20_fun ( n, x, fx )

c*********************************************************************72
c
cc P20_FUN evaluates the integrand for problem 20.
c
c  Interval:
c
c    -1 .le. x .le. 1
c
c  Integrand:
c
c    1 / ( X^2 + 1.005 )
c
c  Antiderivative:
c
c    ( 1 / sqrt ( 1.005 ) ) * arctan ( X / sqrt ( 1.005 ) )
c
c  Approximate Integral (20 digits):
c
c    1.5643964440690497731...
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    18 December 1998
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    David Kahaner,
c    Comparison of Numerical Quadrature Formulas,
c    in Mathematical Software, edited by John R Rice,
c    Academic Press, 1971.
c
c  Parameters:
c
c    Input, integer N, the number of evaluation points.
c
c    Input, double precision X(N), the evaluation points.
c
c    Output, double precision FX(N), the integrand values.
c
      implicit none

      integer n

      double precision fx(n)
      integer i
      double precision x(n)

      do i = 1, n
        fx(i) = 1.0D+00 / ( x(i)**2 + 1.005D+00 )
      end do

      return
      end
      subroutine p20_lim ( a, b )

c*********************************************************************72
c
cc P20_LIM returns the integration limits for problem 20.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    22 December 1998
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, double precision A, B, the limits of integration.
c
      implicit none

      double precision a
      double precision b

      a = -1.0D+00
      b = 1.0D+00

      return
      end
      subroutine p21_exact ( exact )

c*********************************************************************72
c
cc P21_EXACT returns the estimated integral for problem 21.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    18 December 1998
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, double precision EXACT, the estimated value of the integral.
c
      implicit none

      double precision exact

      exact = 0.21080273631018169851D+00

      return
      end
      subroutine p21_fun ( n, x, fx )

c*********************************************************************72
c
cc P21_FUN evaluates the integrand for problem 21.
c
c  Interval:
c
c    0 .le. x .le. 1
c
c  Integrand:
c
c       ( sech (   10.0 * ( x - 0.2 ) ) )^2
c     + ( sech (  100.0 * ( x - 0.4 ) ) )^4
c     + ( sech ( 1000.0 * ( x - 0.6 ) ) )^6
c
c  Exact Integral:
c
c    ( 1 + tanh ( 8 ) * tanh ( 2 ) ) / 10.0 + 2 / 150 + 2 / 1875
c
c  Approximate Integral (20 digits):
c
c    0.21080273631018169851...
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    29 December 2011
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    David Kahaner,
c    Comparison of Numerical Quadrature Formulas,
c    in Mathematical Software, edited by John R Rice,
c    Academic Press, 1971.
c
c  Parameters:
c
c    Input, integer N, the number of evaluation points.
c
c    Input, double precision X(N), the evaluation points.
c
c    Output, double precision FX(N), the integrand values.
c
      implicit none

      integer n

      double precision fx(n)
      integer i
      double precision x(n)

      do i = 1, n
        fx(i) = 
     &      ( 1.0D+00 / cosh (   10.0D+00 * ( x(i) - 0.2D+00 ) ) )**2 
     &    + ( 1.0D+00 / cosh (  100.0D+00 * ( x(i) - 0.4D+00 ) ) )**4 
     &    + ( 1.0D+00 / cosh ( 1000.0D+00 * ( x(i) - 0.6D+00 ) ) )**6
      end do

      return
      end
      subroutine p21_lim ( a, b )

c*********************************************************************72
c
cc P21_LIM returns the integration limits for problem 21.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    22 December 1998
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, double precision A, B, the limits of integration.
c
      implicit none

      double precision a
      double precision b

      a = 0.0D+00
      b = 1.0D+00

      return
      end
      subroutine p22_exact ( exact )

c*********************************************************************72
c
cc P22_EXACT returns the estimated integral for problem 22.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    30 October 2000
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, double precision EXACT, the estimated value of the integral.
c
      implicit none

      double precision exact
      double precision pi
      parameter ( pi = 3.141592653589793D+00 )

      exact = 0.125D+00 * log ( 9.0D+00 ) + pi / sqrt ( 48.0D+00 )

      return
      end
      subroutine p22_fun ( n, x, fx )

c*********************************************************************72
c
cc P22_FUN evaluates the integrand for problem 22.
c
c  Interval:
c
c    0 .le. X .le. 1
c
c  Integrand:
c
c    1 / ( X^4 + X^2 + 1 )
c
c  Exact integral:
c
c    ln ( 9 ) / 8 + pi / sqrt ( 48 )
c
c  Approximate Integral (20 digits):
c
c    0.72810291322558188550...
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    30 October 2000
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Philip Davis, Philip Rabinowitz,
c    Methods of Numerical Integration,
c    Second Edition,
c    Dover, 2007,
c    ISBN: 0486453391,
c    LC: QA299.3.D28.
c
c  Parameters:
c
c    Input, integer N, the number of evaluation points.
c
c    Input, double precision X(N), the evaluation points.
c
c    Output, double precision FX(N), the integrand values.
c
      implicit none

      integer n

      double precision fx(n)
      integer i
      double precision x(n)

      do i = 1, n
        fx(i) = 1.0D+00 / ( x(i)**4 + x(i)**2 + 1.0D+00 )
      end do

      return
      end
      subroutine p22_lim ( a, b )

c*********************************************************************72
c
cc P22_LIM returns the integration limits for problem 22.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    22 December 1998
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, double precision A, B, the limits of integration.
c
      implicit none

      double precision a
      double precision b

      a = 0.0D+00
      b = 1.0D+00

      return
      end
      subroutine p23_exact ( exact )

c*********************************************************************72
c
cc P23_EXACT returns the estimated integral for problem 23.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    23 December 1998
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, double precision EXACT, the estimated value of the integral.
c
      implicit none

      double precision exact

      exact = 0.62471325642771360429D+00

      return
      end
      subroutine p23_fun ( n, x, fx )

c*********************************************************************72
c
cc P23_FUN evaluates the integrand for problem 23.
c
c  Discussion:
c
c    The integrand has a singularity at X = 0.
c    The integrand is discontinuous at X = 0.
c    The integrand is arbitrarily oscillatory as X decreases to 0.
c    The integrand becomes unbounded as X decreases to 0.
c
c    Integral ( 0 < X < 1 ) ( 1 / X ) sin ( 1 / X ) dX
c    = Integral ( 1 < X < +oo ) ( 1 / X ) * sin ( X ) dX.
c
c  Interval:
c
c    0 .le. x .le. 1
c
c  Integrand:
c
c    ( 1 / x ) sin ( 1 / x )
c
c  Approximate Integral (20 digits):
c
c    0.62471325642771360429...
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    23 December 1998
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Philip Davis, Philip Rabinowitz,
c    Methods of Numerical Integration,
c    Second Edition,
c    Dover, 2007,
c    ISBN: 0486453391,
c    LC: QA299.3.D28.
c
c  Parameters:
c
c    Input, integer N, the number of evaluation points.
c
c    Input, double precision X(N), the evaluation points.
c
c    Output, double precision FX(N), the integrand values.
c
      implicit none

      integer n

      double precision fx(n)
      integer i
      double precision x(n)

      do i = 1, n

        if ( x(i) .eq. 0.0D+00 ) then
          fx(i) = 0.0D+00
        else
          fx(i) = ( 1.0D+00 / x(i) ) * sin ( 1.0D+00 / x(i) )
        end if

      end do

      return
      end
      subroutine p23_lim ( a, b )

c*********************************************************************72
c
cc P23_LIM returns the integration limits for problem 23.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    23 December 1998
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, double precision A, B, the limits of integration.
c
      implicit none

      double precision a
      double precision b

      a = 0.0D+00
      b = 1.0D+00

      return
      end
      subroutine p24_exact ( exact )

c*********************************************************************72
c
cc P24_EXACT returns the estimated integral for problem 24.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    09 January 1999
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, double precision EXACT, the estimated value of the integral.
c
      implicit none

      double precision exact

      exact = - 0.0067547455D+00

      return
      end
      subroutine p24_fun ( n, x, fx )

c*********************************************************************72
c
cc P24_FUN evaluates the integrand for problem 24.
c
c  Discussion:
c
c    The integrand is continuous, but nowhere differentiable.
c
c  Interval:
c
c    0 .le. X .le. 0.5
c
c  Integrand:
c
c    ( 1 / pi ) * sum ( 1 .le. I < +oo ) 2^(-I) * cos ( 7^I * pi * X )
c
c  Approximate Integral:
c
c    - 0.0067547455
c
c  Antiderivative:
c
c    ( 1 / pi^2 ) * sum ( 1 .le. I < +oo ) 14^(-I) * sin ( 7^I * pi * X )
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    22 November 1999
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Philip Davis, Philip Rabinowitz,
c    Methods of Numerical Integration,
c    Second Edition,
c    Dover, 2007,
c    ISBN: 0486453391,
c    LC: QA299.3.D28.
c
c    Herbert Salzer, Norman Levine,
c    Table of a Weierstrass Continuous Nondifferentiable Function,
c    Mathematics of Computation,
c    Volume 15, pages 120 - 130, 1961.
c
c  Parameters:
c
c    Input, integer N, the number of evaluation points.
c
c    Input, double precision X(N), the evaluation points.
c
c    Output, double precision FX(N), the integrand values.
c
      implicit none

      integer n

      double precision fx(n)
      integer i
      integer j
      double precision x(n)
      integer, parameter :: n_term = 40
      double precision pi
      parameter ( pi = 3.141592653589793D+00 )

      do i = 1, n

        fx(i) = 0.0D+00
        do j = 1, n_term
          fx(i) = fx(i) + cos ( 7.0D+00**j * pi * x(i) ) / 2.0D+00**j
        end do

        fx(i) = fx(i) / pi

      end do

      return
      end
      subroutine p24_lim ( a, b )

c*********************************************************************72
c
cc P24_LIM returns the integration limits for problem 24.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    09 January 1999
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, double precision A, B, the limits of integration.
c
      implicit none

      double precision a
      double precision b

      a = 0.0D+00
      b = 0.5D+00

      return
      end
      subroutine p25_exact ( exact )

c*********************************************************************72
c
cc P25_EXACT returns the estimated integral for problem 25.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    22 November 1999
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, double precision EXACT, the estimated value of the integral.
c
      implicit none

      double precision exact

      exact = 0.3D+00 * log ( 0.3D+00 ) + 0.7D+00 * log ( 0.7D+00 ) - 1.0D+00

      return
      end
      subroutine p25_fun ( n, x, fx )

c*********************************************************************72
c
cc P25_FUN evaluates the integrand for problem 25.
c
c  Interval:
c
c    0 .le. X .le. 1.
c
c  Integrand:
c
c    ln ( abs ( x - 0.7 ) )
c
c  Exact Integral:
c
c    0.3 * ln ( 0.3 ) + 0.7 * ln ( 0.7 ) - 1
c
c  Approximate Integral (20 digits):
c
c    -1.6108643020548934630
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    22 November 1999
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Kendall Atkinson,
c    An Introduction to Numerical Analysis,
c    Prentice Hall, 1984, page 303.
c
c  Parameters:
c
c    Input, integer N, the number of evaluation points.
c
c    Input, double precision X(N), the evaluation points.
c
c    Output, double precision FX(N), the integrand values.
c
      implicit none

      integer n

      double precision fx(n)
      integer i
      double precision x(n)

      do i = 1, n

        if ( x(i) .eq. 0.7D+00 ) then
          fx(i) = 0.0D+00
        else
          fx(i) = log ( abs ( x(i) - 0.7D+00 ) )
        end if

      end do

      return
      end
      subroutine p25_lim ( a, b )

c*********************************************************************72
c
cc P25_LIM returns the integration limits for problem 25.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    22 November 1999
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, double precision A, B, the limits of integration.
c
      implicit none

      double precision a
      double precision b

      a = 0.0D+00
      b = 1.0D+00

      return
      end
      subroutine p26_exact ( exact )

c*********************************************************************72
c
cc P26_EXACT returns the exact integral for problem 26.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    22 November 1999
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, double precision EXACT, the value of the integral.
c
      implicit none

      double precision exact

      exact = 7.9549265210128452745D+00

      return
      end
      subroutine p26_fun ( n, x, fx )

c*********************************************************************72
c
cc P26_FUN evaluates the integrand for problem 26.
c
c  Interval:
c
c    0 .le. x .le. 2 pi
c
c  Integrand:
c
c    exp ( cos ( x ) )
c
c  Exact Integral:
c
c    2 * Pi * I0(1)
c
c  Approximate Integral (20 digits):
c
c    7.9549265210128452745...
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    22 November 1999
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Kendall Atkinson,
c    An Introduction to Numerical Analysis,
c    Prentice Hall, 1984, page 262.
c
c  Parameters:
c
c    Input, integer N, the number of evaluation points.
c
c    Input, double precision X(N), the evaluation points.
c
c    Output, double precision FX(N), the integrand values.
c
      implicit none

      integer n

      double precision fx(n)
      integer i
      double precision x(n)

      do i = 1, n
        fx(i) = exp ( cos ( x(i) ) )
      end do

      return
      end
      subroutine p26_lim ( a, b )

c*********************************************************************72
c
cc P26_LIM returns the integration limits for problem 26.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    22 November 1999
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, double precision A, B, the limits of integration.
c
      implicit none

      double precision a
      double precision b
      double precision pi
      parameter ( pi = 3.141592653589793D+00 )

      a = 0.0D+00
      b = 2.0D+00 * pi

      return
      end
      subroutine p27_exact ( exact )

c*********************************************************************72
c
cc P27_EXACT returns the exact integral for problem 27.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    09 January 1999
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, double precision EXACT, the value of the integral.
c
      implicit none

      double precision exact

      exact = 5.0D+00 - 6.0D+00 * log ( 2.0D+00 )

      return
      end
      subroutine p27_fun ( n, x, fx )

c*********************************************************************72
c
cc P27_FUN evaluates the integrand for problem 27.
c
c  Interval:
c
c    0 .le. X .le. 1
c
c  Integrand:
c
c    1 / ( X^(1/2) + X^(1/3) )
c
c  Exact Integral:
c
c    5 - 6 * ln ( 2 )
c
c  Approximate Integral (20 digits):
c
c    0.84111691664032814350...
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    22 November 1999
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Philip Davis, Philip Rabinowitz,
c    Methods of Numerical Integration,
c    Second Edition,
c    Dover, 2007,
c    ISBN: 0486453391,
c    LC: QA299.3.D28.
c
c  Parameters:
c
c    Input, integer N, the number of evaluation points.
c
c    Input, double precision X(N), the evaluation points.
c
c    Output, double precision FX(N), the integrand values.
c
      implicit none

      integer n

      double precision fx(n)
      integer i
      double precision x(n)

      do i = 1, n

        if ( x(i) .eq. 0.0D+00 ) then
          fx(i) = 0.0D+00
        else
          fx(i) = 1.0D+00 / ( sqrt ( x(i) ) + x(i)**(1.0D+00/3.0D+00) )
        end if

      end do

      return
      end
      subroutine p27_lim ( a, b )

c*********************************************************************72
c
cc P27_LIM returns the integration limits for problem 27.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    09 January 1999
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, double precision A, B, the limits of integration.
c
      implicit none

      double precision a
      double precision b

      a = 0.0D+00
      b = 1.0D+00

      return
      end
      subroutine p28_exact ( exact )

c*********************************************************************72
c
cc P28_EXACT returns the exact integral for problem 28.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    22 November 1999
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, double precision EXACT, the value of the integral.
c
      implicit none

      double precision exact
      double precision pi
      parameter ( pi = 3.141592653589793D+00 )

      exact = ( 50.0D+00 / 2501.0D+00 ) 
     &  * ( 1.0D+00 - exp ( - 2.0D+00 * pi ) )

      return
      end
      subroutine p28_fun ( n, x, fx )

c*********************************************************************72
c
cc P28_FUN evaluates the integrand for problem 28.
c
c  Interval:
c
c    0 .le. X .le. 2 PI
c
c  Integrand:
c
c    exp ( - X ) * sin ( 50 * X )
c
c  Exact Integral:
c
c    50 / ( 2501 ) * ( 1 - exp ( - 2 * PI ) )
c
c  Approximate Integral (20 digits):
c
c    0.019954669277654778312...
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    22 November 1999
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Kendall Atkinson,
c    An Introduction to Numerical Analysis,
c    Prentice Hall, 1984, page 303.
c
c  Parameters:
c
c    Input, integer N, the number of evaluation points.
c
c    Input, double precision X(N), the evaluation points.
c
c    Output, double precision FX(N), the integrand values.
c
      implicit none

      integer n

      double precision fx(n)
      integer i
      double precision x(n)

      do i = 1, n
        fx(i) = exp ( - x(i) ) * sin ( 50.0D+00 * x(i) )
      end do

      return
      end
      subroutine p28_lim ( a, b )

c*********************************************************************72
c
cc P28_LIM returns the integration limits for problem 28.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    22 November 1999
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, double precision A, B, the limits of integration.
c
      implicit none

      double precision a
      double precision b
      double precision pi
      parameter ( pi = 3.141592653589793D+00 )

      a = 0.0D+00
      b = 2.0D+00 * pi

      return
      end
      subroutine p29_exact ( exact )

c*********************************************************************72
c
cc P29_EXACT returns the exact integral for problem 29.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    11 January 1999
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, double precision EXACT, the value of the integral.
c
      implicit none

      double precision exact

      exact = 1.0D+00 - log ( 2.0D+00 )

      return
      end
      subroutine p29_fun ( n, x, fx )

c*********************************************************************72
c
cc P29_FUN evaluates the integrand for problem 29.
c
c  Interval:
c
c    0 .le. X .le. 1
c
c  Integrand:
c
c    F ( X ) = 1 / ( X + 2 )   for 0 < X < E - 2
c            = 0               otherwise
c
c  Exact Integral:
c
c    1 - ln ( 2 )
c
c  Approximate Integral (20 digits):
c
c    0.30685281944005469058...
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    11 January 1999
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Philip Davis, Philip Rabinowitz,
c    Methods of Numerical Integration,
c    Second Edition,
c    Dover, 2007,
c    ISBN: 0486453391,
c    LC: QA299.3.D28.
c
c  Parameters:
c
c    Input, integer N, the number of evaluation points.
c
c    Input, double precision X(N), the evaluation points.
c
c    Output, double precision FX(N), the integrand values.
c
      implicit none

      integer n

      double precision fx(n)
      integer i
      double precision x(n)

      do i = 1, n

        if ( 0.0D+00 .le. x(i) .and. 
     &    x(i) .le. exp ( 1.0D+00 ) - 2.0D+00 ) then
          fx(i) = 1.0D+00 / ( x(i) + 2.0D+00 )
        else
          fx(i) = 0.0D+00
        end if

      end do

      return
      end
      subroutine p29_lim ( a, b )

c*********************************************************************72
c
cc P29_LIM returns the integration limits for problem 29.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    11 January 1999
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, double precision A, B, the limits of integration.
c
      implicit none

      double precision a
      double precision b

      a = 0.0D+00
      b = 1.0D+00

      return
      end
      subroutine p30_exact ( exact )

c*********************************************************************72
c
cc P30_EXACT returns the exact integral for problem 30.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    04 November 2009
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, double precision EXACT, the estimated value of the integral.
c
      implicit none

      double precision exact

      exact = -4.5275696251606720278D+00

      return
      end
      subroutine p30_fun ( n, x, fx )

c*********************************************************************72
c
cc P30_FUN evaluates the integrand for problem 30.
c
c  Interval:
c
c    2 .le. x .le. 7
c
c  Integrand:
c
c          cos (       x )
c    + 5 * cos ( 1.6 * x )
c    - 2 * cos ( 2.0 * x )
c    + 5 * cos ( 4.5 * x )
c    + 7 * cos ( 9.0 * x )
c
c  Antiderivative:
c
c          sin (       x )
c    + 5 * sin ( 1.6 * x ) / 1.6
c    - 2 * sin ( 2.0 * x ) / 2.0
c    + 5 * sin ( 4.5 * x ) / 4.5
c    + 7 * sin ( 9.0 * x ) / 9.0
c
c  Exact Integral:
c
c    -4.5275696251606720278
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    05 November 2009
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Dianne OLeary,
c    Scientific Computing with Case Studies,
c    SIAM, 2008,
c    ISBN13: 978-0-898716-66-5,
c    LC: QA401.O44.
c
c  Parameters:
c
c    Input, integer N, the number of evaluation points.
c
c    Input, double precision X(N), the evaluation points.
c
c    Output, double precision FX(N), the integrand values.
c
      implicit none

      integer n

      double precision fx(n)
      integer i
      double precision x(n)

      do i = 1, n
        fx(i) = 
     &              cos (           x(i) ) 
     &  + 5.0D+00 * cos ( 1.6D+00 * x(i) ) 
     &  - 2.0D+00 * cos ( 2.0D+00 * x(i) ) 
     &  + 5.0D+00 * cos ( 4.5D+00 * x(i) ) 
     &  + 7.0D+00 * cos ( 9.0D+00 * x(i) )
      end do

      return
      end
      subroutine p30_lim ( a, b )

c*********************************************************************72
c
cc P30_LIM returns the integration limits for problem 30.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    04 November 2009
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, double precision A, B, the limits of integration.
c
      implicit none

      double precision a
      double precision b

      a = 2.0D+00
      b = 7.0D+00

      return
      end
      subroutine p31_exact ( exact )

c*********************************************************************72
c
cc P31_EXACT returns the exact integral for problem 31.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    21 November 1999
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, double precision EXACT, the value of the integral.
c
      implicit none

      double precision exact

      exact = 2.0D+00 * atan ( 4.0D+00 )

      return
      end
      subroutine p31_fun ( n, x, fx )

c*********************************************************************72
c
cc P31_FUN evaluates the integrand for problem 31.
c
c  Discussion:
c
c    A simple Newton-Cotes quadrature rule, in which the order of the
c    rule is increased, but the interval is not subdivided, diverges
c    for this integrand.
c
c    This is Runge's function, a standard example of the perils of
c    using high order polynomial interpolation at equally spaced nodes.
c    Since this is exactly what is really going on in a Newton Cotes
c    rule, it is little wonder that the result is so poor.
c
c  Interval:
c
c    -4 .le. x .le. 4
c
c  Integrand:
c
c    1 / ( 1 + x^2 )
c
c  Antiderivative:
c
c    arctan ( x )
c
c  Exact Integral:
c
c    2 * arctan ( 4 )
c
c  Approximate Integral (20 digits):
c
c    2.6516353273360649301...
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    21 November 1999
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Kendall Atkinson,
c    An Introduction to Numerical Analysis,
c    Prentice Hall, 1984, page 266.
c
c  Parameters:
c
c    Input, integer N, the number of evaluation points.
c
c    Input, double precision X(N), the evaluation points.
c
c    Output, double precision FX(N), the integrand values.
c
      implicit none

      integer n

      double precision fx(n)
      integer i
      double precision x(n)

      do i = 1, n
        fx(i) = 1.0D+00 / ( 1.0D+00 + x(i)**2 )
      end do

      return
      end
      subroutine p31_lim ( a, b )

c*********************************************************************72
c
cc P31_LIM returns the integration limits for problem 31.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    21 November 1999
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, double precision A, B, the limits of integration.
c
      implicit none

      double precision a
      double precision b

      a = - 4.0D+00
      b =   4.0D+00

      return
      end
      subroutine p32_exact ( exact )

c*********************************************************************72
c
cc P32_EXACT returns the exact integral for problem 32.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    21 November 1999
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, double precision EXACT, the value of the integral.
c
      implicit none

      double precision exact
      double precision pi
      parameter ( pi = 3.141592653589793D+00 )

      exact = - 0.5D+00 * ( exp ( pi ) + 1.0D+00 )

      return
      end
      subroutine p32_fun ( n, x, fx )

c*********************************************************************72
c
cc P32_FUN evaluates the integrand for problem 32.
c
c  Interval:
c
c    0 .le. X .le. PI
c
c  Integrand:
c
c    exp ( X ) * cos ( X )
c
c  Antiderivative:
c
c    0.5 * exp ( X ) * ( sin ( X ) + cos ( X ) )
c
c  Exact Integral:
c
c    - 0.5 * ( exp ( PI ) + 1 )
c
c  Approximate Integral (20 digits):
c
c    -12.070346316389634503...
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    21 November 1999
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Kendall Atkinson,
c    An Introduction to Numerical Analysis,
c    Prentice Hall, 1984, page 254, 277, 297.
c
c  Parameters:
c
c    Input, integer N, the number of evaluation points.
c
c    Input, double precision X(N), the evaluation points.
c
c    Output, double precision FX(N), the integrand values.
c
      implicit none

      integer n

      double precision fx(n)
      integer i
      double precision x(n)

      do i = 1, n
        fx(i) = exp ( x(i) ) * cos ( x(i) )
      end do

      return
      end
      subroutine p32_lim ( a, b )

c*********************************************************************72
c
cc P32_LIM returns the integration limits for problem 32.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    21 November 1999
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, double precision A, B, the limits of integration.
c
      implicit none

      double precision a
      double precision b
      double precision pi
      parameter ( pi = 3.141592653589793D+00 )

      a = 0.0D+00
      b = pi

      return
      end
      subroutine p33_exact ( exact )

c*********************************************************************72
c
cc P33_EXACT returns the exact integral for problem 33.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    21 November 1999
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, double precision EXACT, the value of the integral.
c
      implicit none

      double precision exact
      double precision pi
      parameter ( pi = 3.141592653589793D+00 )

      exact = 0.5D+00 * sqrt ( pi )

      return
      end
      subroutine p33_fun ( n, x, fx )

c*********************************************************************72
c
cc P33_FUN evaluates the integrand for problem 33.
c
c  Discussion:
c
c    The integrand is singular at both endpoints of the interval.
c
c  Interval:
c
c    0 .le. X .le. 1
c
c  Integrand:
c
c    sqrt ( - ln ( X ) )
c
c  Exact Integral:
c
c    sqrt ( pi ) / 2
c
c  Approximate Integral (20 digits):
c
c    0.88622692545275801365...
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    22 November 1999
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Kendall Atkinson,
c    An Introduction to Numerical Analysis,
c    Prentice Hall, 1984, page 307.
c
c  Parameters:
c
c    Input, integer N, the number of evaluation points.
c
c    Input, double precision X(N), the evaluation points.
c
c    Output, double precision FX(N), the integrand values.
c
      implicit none

      integer n

      double precision fx(n)
      integer i
      double precision x(n)

      do i = 1, n

        if ( x(i) .le. 0.0D+00 ) then
          fx(i) = 0.0D+00
        else
          fx(i) = sqrt ( - log ( x(i) ) )
        end if

      end do

      return
      end
      subroutine p33_lim ( a, b )

c*********************************************************************72
c
cc P33_LIM returns the integration limits for problem 33.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    21 November 1999
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, double precision A, B, the limits of integration.
c
      implicit none

      double precision a
      double precision b

      a = 0.0D+00
      b = 1.0D+00

      return
      end
      subroutine p34_exact ( exact )

c*********************************************************************72
c
cc P34_EXACT returns the exact integral for problem 34.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    29 October 2000
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, double precision EXACT, the value of the integral.
c
      implicit none

      double precision exact

      exact = 1627879.0D+00 / 1500.0D+00

      return
      end
      subroutine p34_fun ( n, x, fx )

c*********************************************************************72
c
cc P34_FUN evaluates the integrand for problem 34.
c
c  Interval:
c
c    0 .le. X .le. 1
c
c  Integrand:
c
c    ( 10 * X - 1 ) * ( 10 * X - 1.1 ) * ( 10 * X - 1.2 ) * ( 10 * X - 1.3 )
c
c  Exact Integral:
c
c    1627879 / 1500
c
c  Approximate Integral (20 digits):
c
c    1085.2526666666666666...
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    29 October 2000
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Hermann Engels,
c    Numerical Quadrature and Cubature,
c    Academic Press, 1980.
c
c  Parameters:
c
c    Input, integer N, the number of evaluation points.
c
c    Input, double precision X(N), the evaluation points.
c
c    Output, double precision FX(N), the integrand values.
c
      implicit none

      integer n

      double precision fx(n)
      integer i
      double precision x(n)

      do i = 1, n
        fx(i) = ( 10.0D+00 * x(i) - 1.0D+00 ) 
     &    * ( 10.0D+00 * x(i) - 1.1D+00 ) 
     &    * ( 10.0D+00 * x(i) - 1.2D+00 ) 
     &    * ( 10.0D+00 * x(i) - 1.3D+00 )
      end do

      return
      end
      subroutine p34_lim ( a, b )

c*********************************************************************72
c
cc P34_LIM returns the integration limits for problem 34.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    29 October 2000
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, double precision A, B, the limits of integration.
c
      implicit none

      double precision a
      double precision b

      a = 0.0D+00
      b = 1.0D+00

      return
      end
      subroutine p35_exact ( exact )

c*********************************************************************72
c
cc P35_EXACT returns the exact integral for problem 35.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    31 October 2000
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, double precision EXACT, the value of the integral.
c
      implicit none

      double precision exact

      exact = 26.0D+00

      return
      end
      subroutine p35_fun ( n, x, fx )

c*********************************************************************72
c
cc P35_FUN evaluates the integrand for problem 35.
c
c  Interval:
c
c    -9 .le. X .le. 100
c
c  Integrand:
c
c    1 / sqrt ( abs ( X ) )
c
c  Exact Integral:
c
c    26
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    31 October 2000
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Hermann Engels,
c    Numerical Quadrature and Cubature,
c    Academic Press, 1980.
c
c  Parameters:
c
c    Input, integer N, the number of evaluation points.
c
c    Input, double precision X(N), the evaluation points.
c
c    Output, double precision FX(N), the integrand values.
c
      implicit none

      integer n

      double precision fx(n)
      integer i
      double precision x(n)

      do i = 1, n

        if ( x(i) .eq. 0.0D+00 ) then
          fx(i) = 0.0D+00
        else
          fx(i) = 1.0D+00 / sqrt ( abs ( x(i) ) )
        end if

      end do

      return
      end
      subroutine p35_lim ( a, b )

c*********************************************************************72
c
cc P35_LIM returns the integration limits for problem 35.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    31 October 2000
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, double precision A, B, the limits of integration.
c
      implicit none

      double precision a
      double precision b

      a = -9.0D+00
      b = 100.0D+00

      return
      end
      subroutine p36_exact ( exact )

c*********************************************************************72
c
cc P36_EXACT returns the exact integral for problem 36.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    09 March 2002
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, double precision EXACT, the value of the integral.
c
      implicit none

      double precision alpha
      double precision exact

      call p36_param_get ( alpha )

      exact = 1.0D+00 / ( alpha + 1.0D+00 )**2

      return
      end
      subroutine p36_fun ( n, x, fx )

c*********************************************************************72
c
cc P36_FUN evaluates the integrand for problem 36.
c
c  Discussion:
c
c    The problem has a parameter ALPHA that can be set by calling
c    P36_PARAM_SET.  It had a default value of -0.9.
c
c    The integrand has an endpoint singularity at X=0.
c
c    Suggested values of ALPHA include -0.9 through 2.6.
c
c  Interval:
c
c    0 .le. X .le. 1
c
c  Integrand:
c
c    X^alpha * ln ( 1 / X )
c
c  Exact Integral:
c
c    1 / ( alpha + 1 )^2
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    09 March 2002
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Robert Piessens, Elise de Doncker-Kapenga,
c    Christian Ueberhuber, David Kahaner,
c    QUADPACK: A Subroutine Package for Automatic Integration,
c    Springer, 1983, page 83.
c
c  Parameters:
c
c    Input, integer N, the number of evaluation points.
c
c    Input, double precision X(N), the evaluation points.
c
c    Output, double precision FX(N), the integrand values.
c
      implicit none

      integer n

      double precision alpha
      double precision fx(n)
      integer i
      double precision x(n)
      double precision value

      call p36_param_get ( alpha )

      do i = 1, n

        if ( x(i) .le. 0.0D+00 ) then
          fx(i) = 0.0D+00
        else
          fx(i) = x(i)**(alpha) * log ( 1.0D+00 / x(i) )
        end if

      end do

      return
      end
      subroutine p36_lim ( a, b )

c*********************************************************************72
c
cc P36_LIM returns the integration limits for problem 36.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    09 March 2002
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, double precision A, B, the limits of integration.
c
      implicit none

      double precision a
      double precision b

      a = 0.0D+00
      b = 1.0D+00

      return
      end
      subroutine p36_param ( action, name, value )

c*********************************************************************72
c
cc P36_PARAM gets or sets the parameter values for problem 36.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    09 March 2002
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, character ( len = * ) ACTION.
c    'get' to get the value.
c    'set' to set the value.
c
c    Input, character ( len = * ) NAME, the name of the parameter.
c    'alpha' is the only option.
c
c    Input/output, double precision VALUE.
c    If the action is 'get', then VALUE returns the current parameter value.
c    If ACTION is 'set', then the parameter value is set to VALUE.
c
      implicit none

      character ( len = * ) action
      double precision alpha
      character ( len = * ) name
      logical s_eqi
      double precision value

      save alpha

      data alpha / - 0.9D+00 /

      if ( s_eqi ( action, 'get' ) ) then
        if ( s_eqi ( name, 'alpha' ) ) then
          value = alpha
        else
          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) 'P36_PARAM - Fatal error!'
          write ( *, '(a)' ) '  Unrecognized name.'
          stop
        end if
      else if ( s_eqi ( action, 'set' ) ) then
        if ( s_eqi ( name, 'alpha' ) ) then
          alpha = value
        else
          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) 'P36_PARAM - Fatal error!'
          write ( *, '(a)' ) '  Unrecognized name.'
          stop
        end if
      else
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'P36_PARAM - Fatal error!'
        write ( *, '(a)' ) '  Unrecognized action.'
        stop
      end if

      return
      end
      subroutine p36_param_get ( alpha )

c*********************************************************************72
c
cc P36_PARAM_GET returns the parameter values for problem 36.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    09 March 2002
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, double precision ALPHA, the current value of the parameter.
c
      implicit none

      double precision alpha

      call p36_param ( 'get', 'alpha', alpha )

      return
      end
      subroutine p36_param_set ( alpha )

c*********************************************************************72
c
cc P36_PARAM_SET sets the parameter values for problem 36.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    09 March 2002
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, double precision ALPHA, the new value of the parameter.
c
      implicit none

      double precision alpha

      call p36_param ( 'set', 'alpha', alpha )

      return
      end
      subroutine p37_exact ( exact )

c*********************************************************************72
c
cc P37_EXACT returns the exact integral for problem 37.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    09 March 2002
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, double precision EXACT, the value of the integral.
c
      implicit none

      double precision alpha
      double precision exact
      double precision pi
      parameter ( pi = 3.141592653589793D+00 )

      call p37_param_get ( alpha )

      exact = atan ( ( 4.0D+00 - pi ) * 4.0D+00**( alpha - 1.0D+00 ) ) 
     &      + atan (             pi   * 4.0D+00**( alpha - 1.0D+00 ) )

      return
      end
      subroutine p37_fun ( n, x, fx )

c*********************************************************************72
c
cc P37_FUN evaluates the integrand for problem 37.
c
c  Discussion:
c
c    The problem has a parameter ALPHA that can be set by calling
c    P37_PARAM_SET.  It had a default value of 5.0.
c
c    The integrand has a peak of height 4^ALPHA at X = PI/4.
c
c    Suggested values of ALPHA include 0 through 20.
c
c  Interval:
c
c    0 .le. X .le. 1
c
c  Integrand:
c
c    4^(-ALPHA) / ( (X-PI/4)^2 + 16^(-ALPHA) )
c
c  Exact Integral:
c
c    atan ( ( 4 - PI ) * 4^(ALPHA-1) ) + atan ( PI * 4^(ALPHA-1) )
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    09 March 2002
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Robert Piessens, Elise de Doncker-Kapenga,
c    Christian Ueberhuber, David Kahaner,
c    QUADPACK: A Subroutine Package for Automatic Integration,
c    Springer, 1983, page 83.
c
c  Parameters:
c
c    Input, integer N, the number of evaluation points.
c
c    Input, double precision X(N), the evaluation points.
c
c    Output, double precision FX(N), the integrand values.
c
      implicit none

      integer n

      double precision alpha
      double precision fx(n)
      integer i
      double precision pi
      parameter ( pi = 3.141592653589793D+00 )
      double precision x(n)

      call p37_param_get ( alpha )

      do i = 1, n

        fx(i) = 4.0D+00 ** ( -alpha ) 
     &    / ( ( x(i) - 0.25D+00 * pi )**2 + 16.0D+00**(-alpha) )

      end do

      return
      end
      subroutine p37_lim ( a, b )

c*********************************************************************72
c
cc P37_LIM returns the integration limits for problem 37.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    09 March 2002
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, double precision A, B, the limits of integration.
c
      implicit none

      double precision a
      double precision b

      a = 0.0D+00
      b = 1.0D+00

      return
      end
      subroutine p37_param ( action, name, value )

c*********************************************************************72
c
cc P37_PARAM gets or sets the parameter values for problem 37.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    09 March 2002
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, character ( len = * ) ACTION.
c    'get' to get the value.
c    'set' to set the value.
c
c    Input, character ( len = * ) NAME, the name of the parameter.
c    'alpha' is the only option.
c
c    Input/output, double precision VALUE.
c    If the action is 'get', then VALUE returns the current parameter value.
c    If ACTION is 'set', then the parameter value is set to VALUE.
c
      implicit none

      character ( len = * ) action
      double precision alpha
      character ( len = * ) name
      logical s_eqi
      double precision value

      save alpha

      data alpha / 5.0D+00 /

      if ( s_eqi ( action, 'get' ) ) then
        if ( s_eqi ( name, 'alpha' ) ) then
          value = alpha
        else
          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) 'P37_PARAM - Fatal error!'
          write ( *, '(a)' ) '  Unrecognized name.'
          stop
        end if
      else if ( s_eqi ( action, 'set' ) ) then
        if ( s_eqi ( name, 'alpha' ) ) then
          alpha = value
        else
          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) 'P37_PARAM - Fatal error!'
          write ( *, '(a)' ) '  Unrecognized name.'
          stop
        end if
      else
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'P37_PARAM - Fatal error!'
        write ( *, '(a)' ) '  Unrecognized action.'
        stop
      end if

      return
      end
      subroutine p37_param_get ( alpha )

c*********************************************************************72
c
cc P37_PARAM_GET returns the parameter values for problem 37.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    09 March 2002
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, double precision ALPHA, the current value of the parameter.
c
      implicit none

      double precision alpha

      call p37_param ( 'get', 'alpha', alpha )

      return
      end
      subroutine p37_param_set ( alpha )

c*********************************************************************72
c
cc P37_PARAM_SET sets the parameter values for problem 37.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    09 March 2002
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, double precision ALPHA, the new value of the parameter.
c
      implicit none

      double precision alpha

      call p37_param ( 'set', 'alpha', alpha )

      return
      end
      subroutine p38_exact ( exact )

c*********************************************************************72
c
cc P38_EXACT returns the exact integral for problem 38.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    09 March 2002
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, double precision EXACT, the value of the integral.
c
      implicit none

      double precision alpha
      double precision besj0
      double precision exact
      double precision pi
      parameter ( pi = 3.141592653589793D+00 )
      double precision x

      call p38_param_get ( alpha )

      x = 2.0D+00**alpha

      exact = pi * besj0 ( x )

      return
      end
      subroutine p38_fun ( n, x, fx )

c*********************************************************************72
c
cc P38_FUN evaluates the integrand for problem 38.
c
c  Discussion:
c
c    The problem has a parameter ALPHA that can be set by calling
c    P38_PARAM_SET.
c
c    The integrand oscillates more strongly as ALPHA is increased.
c
c    The suggested range for ALPHA is 0 to 10.
c
c  Interval:
c
c    0 .le. X .le. PI
c
c  Integrand:
c
c    cos ( 2^ALPHA * sin ( x ) )
c
c  Exact Integral:
c
c    pi * J0 ( 2^ALPHA )
c
c    where J0 ( x ) is the J Bessel function of order 0.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    09 March 2002
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Robert Piessens, Elise de Doncker-Kapenga,
c    Christian Ueberhuber, David Kahaner,
c    QUADPACK: A Subroutine Package for Automatic Integration,
c    Springer, 1983, page 83.
c
c  Parameters:
c
c    Input, integer N, the number of evaluation points.
c
c    Input, double precision X(N), the evaluation points.
c
c    Output, double precision FX(N), the integrand values.
c
      implicit none

      integer n

      double precision fx(n)
      integer i
      double precision x(n)
      double precision alpha
      double precision value

      call p38_param_get ( alpha )

      do i = 1, n
        fx(i) = cos ( 2.0D+00**alpha * sin ( x(i) ) )
      end do

      return
      end
      subroutine p38_lim ( a, b )

c*********************************************************************72
c
cc P38_LIM returns the integration limits for problem 38.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    09 March 2002
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, double precision A, B, the limits of integration.
c
      implicit none

      double precision a
      double precision b
      double precision pi
      parameter ( pi = 3.141592653589793D+00 )

      a = 0.0D+00
      b = pi

      return
      end
      subroutine p38_param ( action, name, value )

c*********************************************************************72
c
cc P38_PARAM gets or sets the parameter values for problem 38.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    09 March 2002
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, character ( len = * ) ACTION.
c    'get' to get the value.
c    'set' to set the value.
c
c    Input, character ( len = * ) NAME, the name of the parameter.
c    'alpha' is the only option.
c
c    Input/output, double precision VALUE.
c    If the action is 'get', then VALUE returns the current parameter value.
c    If ACTION is 'set', then the parameter value is set to VALUE.
c
      implicit none

      character ( len = * ) action
      double precision alpha
      character ( len = * ) name
      logical s_eqi
      double precision value

      save alpha

      data alpha / 3.0D+00 /

      if ( s_eqi ( action, 'get' ) ) then
        if ( s_eqi ( name, 'alpha' ) ) then
          value = alpha
        else
          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) 'P38_PARAM - Fatal error!'
          write ( *, '(a)' ) '  Unrecognized name.'
          stop
        end if
      else if ( s_eqi ( action, 'set' ) ) then
        if ( s_eqi ( name, 'alpha' ) ) then
          alpha = value
        else
          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) 'P38_PARAM - Fatal error!'
          write ( *, '(a)' ) '  Unrecognized name.'
          stop
        end if
      else
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'P38_PARAM - Fatal error!'
        write ( *, '(a)' ) '  Unrecognized action.'
        stop
      end if

      return
      end
      subroutine p38_param_get ( alpha )

c*********************************************************************72
c
cc P38_PARAM_GET returns the parameter values for problem 38.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    09 March 2002
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, double precision ALPHA, the current value of the parameter.
c
      implicit none

      double precision alpha

      call p38_param ( 'get', 'alpha', alpha )

      return
      end
      subroutine p38_param_set ( alpha )

c*********************************************************************72
c
cc P38_PARAM_SET sets the parameter values for problem 38.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    09 March 2002
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, double precision ALPHA, the new value of the parameter.
c
      implicit none

      double precision alpha

      call p38_param ( 'set', 'alpha', alpha )

      return
      end
      subroutine p39_exact ( exact )

c*********************************************************************72
c
cc P39_EXACT returns the exact integral for problem 39.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    09 March 2002
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, double precision EXACT, the value of the integral.
c
      implicit none

      double precision alpha
      double precision exact

      call p39_param_get ( alpha )

      exact = ( ( 2.0D+00 / 3.0D+00 )**( alpha + 1.0D+00 ) 
     &        + ( 1.0D+00 / 3.0D+00 )**( alpha + 1.0D+00 ) ) 
     &        / ( alpha + 1.0D+00 )

      return
      end
      subroutine p39_fun ( n, x, fx )

c*********************************************************************72
c
cc P39_FUN evaluates the integrand for problem 39.
c
c  Discussion:
c
c    The problem has a parameter ALPHA that can be set by calling
c    P39_PARAM_SET.
c
c    The integrand has a singularity at an internal point ( x = 1/3 )
c    with small binary period.
c
c    The suggested range for ALPHA is -0.8 through 2.1.
c
c  Interval:
c
c    0 .le. X .le. 1
c
c  Integrand:
c
c    ( abs ( x - 1/3 ) )^alpha
c
c  Exact Integral:
c
c    ( (2/3)^(alpha+1) + (1/3)^(alpha+1) ) / ( alpha + 1 )
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    14 November 2009
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Robert Piessens, Elise de Doncker-Kapenga,
c    Christian Ueberhuber, David Kahaner,
c    QUADPACK: A Subroutine Package for Automatic Integration,
c    Springer, 1983, page 83.
c
c  Parameters:
c
c    Input, integer N, the number of evaluation points.
c
c    Input, double precision X(N), the evaluation points.
c
c    Output, double precision FX(N), the integrand values.
c
      implicit none

      integer n

      double precision fx(n)
      integer i
      double precision x(n)
      double precision alpha

      call p39_param_get ( alpha )

      do i = 1, n

        if ( x(i) - 1.0D+00 / 3.0D+00 .eq. 0.0D+00 ) then
          fx(i) = 0.0D+00
        else
          fx(i) = ( abs ( x(i) - 1.0D+00 / 3.0D+00 ) )**alpha
        end if

      end do

      return
      end
      subroutine p39_lim ( a, b )

c*********************************************************************72
c
cc P39_LIM returns the integration limits for problem 39.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    09 March 2002
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, double precision A, B, the limits of integration.
c
      implicit none

      double precision a
      double precision b

      a = 0.0D+00
      b = 1.0D+00

      return
      end
      subroutine p39_param ( action, name, value )

c*********************************************************************72
c
cc P39_PARAM gets or sets the parameter values for problem 39.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    09 March 2002
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, character ( len = * ) ACTION.
c    'get' to get the value.
c    'set' to set the value.
c
c    Input, character ( len = * ) NAME, the name of the parameter.
c    'alpha' is the only option.
c
c    Input/output, double precision VALUE.
c    If the action is 'get', then VALUE returns the current parameter value.
c    If ACTION is 'set', then the parameter value is set to VALUE.
c
      implicit none

      character ( len = * ) action
      double precision alpha
      character ( len = * ) name
      logical s_eqi
      double precision value

      save alpha

      data alpha / -0.5D+00 /

      if ( s_eqi ( action, 'get' ) ) then
        if ( s_eqi ( name, 'alpha' ) ) then
          value = alpha
        else
          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) 'P39_PARAM - Fatal error!'
          write ( *, '(a)' ) '  Unrecognized name.'
          stop
        end if
      else if ( s_eqi ( action, 'set' ) ) then
        if ( s_eqi ( name, 'alpha' ) ) then
          alpha = value
        else
          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) 'P39_PARAM - Fatal error!'
          write ( *, '(a)' ) '  Unrecognized name.'
          stop
        end if
      else
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'P39_PARAM - Fatal error!'
        write ( *, '(a)' ) '  Unrecognized action.'
        stop
      end if

      return
      end
      subroutine p39_param_get ( alpha )

c*********************************************************************72
c
cc P39_PARAM_GET returns the parameter values for problem 39.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    09 March 2002
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, double precision ALPHA, the current value of the parameter.
c
      implicit none

      double precision alpha

      call p39_param ( 'get', 'alpha', alpha )

      return
      end
      subroutine p39_param_set ( alpha )

c*********************************************************************72
c
cc P39_PARAM_SET sets the parameter values for problem 39.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    09 March 2002
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, double precision ALPHA, the new value of the parameter.
c
      implicit none

      double precision alpha

      call p39_param ( 'set', 'alpha', alpha )

      return
      end
      subroutine p40_exact ( exact )

c*********************************************************************72
c
cc P40_EXACT returns the exact integral for problem 40.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    09 March 2002
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, double precision EXACT, the value of the integral.
c
      implicit none

      double precision alpha
      double precision exact
      double precision pi
      parameter ( pi = 3.141592653589793D+00 )

      call p40_param_get ( alpha )

      exact = ( ( 1.0 - 0.25D+00 * pi )**( alpha + 1.0D+00 ) 
     &        + (     + 0.25D+00 * pi )**( alpha + 1.0D+00 ) ) 
     &        / ( alpha + 1.0D+00 )

      return
      end
      subroutine p40_fun ( n, x, fx )

c*********************************************************************72
c
cc P40_FUN evaluates the integrand for problem 40.
c
c  Discussion:
c
c    The problem has a parameter ALPHA that can be set by calling
c    P40_PARAM_SET.
c
c    The integrand has a singularity at an internal point ( x = pi/4 ).
c
c    The suggested range for ALPHA is -0.8 through 2.1.
c
c  Interval:
c
c    0 .le. X .le. 1
c
c  Integrand:
c
c    ( abs ( x - pi/4 ) )^alpha
c
c  Exact Integral:
c
c    ( (1-pi/4)^(alpha+1) + (pi/4)^(alpha+1) ) / ( alpha + 1 )
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    09 March 2002
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Robert Piessens, Elise de Doncker-Kapenga,
c    Christian Ueberhuber, David Kahaner,
c    QUADPACK: A Subroutine Package for Automatic Integration,
c    Springer, 1983, page 83.
c
c  Parameters:
c
c    Input, integer N, the number of evaluation points.
c
c    Input, double precision X(N), the evaluation points.
c
c    Output, double precision FX(N), the integrand values.
c
      implicit none

      integer n

      double precision fx(n)
      integer i
      double precision x(n)
      double precision alpha
      double precision pi
      parameter ( pi = 3.141592653589793D+00 )

      call p40_param_get ( alpha )

      do i = 1, n
        fx(i) = ( abs ( x(i) - 0.25D+00 * pi ) )**alpha
      end do

      return
      end
      subroutine p40_lim ( a, b )

c*********************************************************************72
c
cc P40_LIM returns the integration limits for problem 40.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    09 March 2002
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, double precision A, B, the limits of integration.
c
      implicit none

      double precision a
      double precision b

      a = 0.0D+00
      b = 1.0D+00

      return
      end
      subroutine p40_param ( action, name, value )

c*********************************************************************72
c
cc P40_PARAM gets or sets the parameter values for problem 40.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    09 March 2002
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, character ( len = * ) ACTION.
c    'get' to get the value.
c    'set' to set the value.
c
c    Input, character ( len = * ) NAME, the name of the parameter.
c    'alpha' is the only option.
c
c    Input/output, double precision VALUE.
c    If the action is 'get', then VALUE returns the current parameter value.
c    If ACTION is 'set', then the parameter value is set to VALUE.
c
      implicit none

      character ( len = * ) action
      double precision alpha
      character ( len = * ) name
      logical s_eqi
      double precision value

      save alpha

      data alpha / -0.5D+00 /

      if ( s_eqi ( action, 'get' ) ) then
        if ( s_eqi ( name, 'alpha' ) ) then
          value = alpha
        else
          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) 'P40_PARAM - Fatal error!'
          write ( *, '(a)' ) '  Unrecognized name.'
          stop
        end if
      else if ( s_eqi ( action, 'set' ) ) then
        if ( s_eqi ( name, 'alpha' ) ) then
          alpha = value
        else
          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) 'P40_PARAM - Fatal error!'
          write ( *, '(a)' ) '  Unrecognized name.'
          stop
        end if
      else
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'P40_PARAM - Fatal error!'
        write ( *, '(a)' ) '  Unrecognized action.'
        stop
      end if

      return
      end
      subroutine p40_param_get ( alpha )

c*********************************************************************72
c
cc P40_PARAM_GET returns the parameter values for problem 40.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    09 March 2002
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, double precision ALPHA, the current value of the parameter.
c
      implicit none

      double precision alpha

      call p40_param ( 'get', 'alpha', alpha )

      return
      end
      subroutine p40_param_set ( alpha )

c*********************************************************************72
c
cc P40_PARAM_SET sets the parameter values for problem 40.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    09 March 2002
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, double precision ALPHA, the new value of the parameter.
c
      implicit none

      double precision alpha

      call p40_param ( 'set', 'alpha', alpha )

      return
      end
      subroutine p41_exact ( exact )

c*********************************************************************72
c
cc P41_EXACT returns the exact integral for problem 41.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    10 March 2002
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, double precision EXACT, the value of the integral.
c
      implicit none

      double precision alpha
      double precision exact
      double precision pi
      parameter ( pi = 3.141592653589793D+00 )

      call p41_param_get ( alpha )

      exact = pi / sqrt ( ( 1.0D+00 + 2.0D+00**(-alpha) )**2 - 1.0D+00 )

      return
      end
      subroutine p41_fun ( n, x, fx )

c*********************************************************************72
c
cc P41_FUN evaluates the integrand for problem 41.
c
c  Discussion:
c
c    The problem has a parameter ALPHA that can be set by calling
c    P41_PARAM_SET.
c
c    The integrand has a singularity at both endpoints, whose
c    severity increases with ALPHA.
c
c    The suggested range for ALPHA is 1 through 20.
c
c  Interval:
c
c    -1 .le. X .le. 1
c
c  Integrand:
c
c    1 / ( sqrt ( 1 - x^2 ) * ( x + 1 + 2^(-alpha) ) )
c
c  Exact Integral:
c
c    pi / sqrt ( ( 1 + 2^(-alpha) ) - 1 )
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    14 November 2009
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Robert Piessens, Elise de Doncker-Kapenga,
c    Christian Ueberhuber, David Kahaner,
c    QUADPACK: A Subroutine Package for Automatic Integration,
c    Springer, 1983, page 83.
c
c  Parameters:
c
c    Input, integer N, the number of evaluation points.
c
c    Input, double precision X(N), the evaluation points.
c
c    Output, double precision FX(N), the integrand values.
c
      implicit none

      integer n

      double precision fx(n)
      integer i
      double precision x(n)
      double precision alpha

      call p41_param_get ( alpha )

      do i = 1, n

        if ( 1.0D+00 - x(i)**2 .eq. 0.0D+00 .or. 
     &       x(i) + 1.0D+00 + 0.5D+00**alpha .eq. 0.0D+00 ) then
          fx(i) = 0.0D+00
        else
          fx(i) = 1.0D+00 / ( sqrt ( 1.0D+00 - x(i)**2 ) 
     &      * ( x(i) + 1.0D+00 + 0.5D+00**alpha ) )
        end if

      end do

      return
      end
      subroutine p41_lim ( a, b )

c*********************************************************************72
c
cc P41_LIM returns the integration limits for problem 41.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    10 March 2002
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, double precision A, B, the limits of integration.
c
      implicit none

      double precision a
      double precision b

      a = -1.0D+00
      b = 1.0D+00

      return
      end
      subroutine p41_param ( action, name, value )

c*********************************************************************72
c
cc P41_PARAM gets or sets the parameter values for problem 41.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    10 March 2002
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, character ( len = * ) ACTION.
c    'get' to get the value.
c    'set' to set the value.
c
c    Input, character ( len = * ) NAME, the name of the parameter.
c    'alpha' is the only option.
c
c    Input/output, double precision VALUE.
c    If the action is 'get', then VALUE returns the current parameter value.
c    If ACTION is 'set', then the parameter value is set to VALUE.
c
      implicit none

      character ( len = * ) action
      double precision alpha
      character ( len = * ) name
      logical s_eqi
      double precision value

      save alpha

      data alpha / 3.0D+00 /

      if ( s_eqi ( action, 'get' ) ) then
        if ( s_eqi ( name, 'alpha' ) ) then
          value = alpha
        else
          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) 'P41_PARAM - Fatal error!'
          write ( *, '(a)' ) '  Unrecognized name.'
          stop
        end if
      else if ( s_eqi ( action, 'set' ) ) then
        if ( s_eqi ( name, 'alpha' ) ) then
          alpha = value
        else
          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) 'P41_PARAM - Fatal error!'
          write ( *, '(a)' ) '  Unrecognized name.'
          stop
        end if
      else
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'P41_PARAM - Fatal error!'
        write ( *, '(a)' ) '  Unrecognized action.'
        stop
      end if

      return
      end
      subroutine p41_param_get ( alpha )

c*********************************************************************72
c
cc P41_PARAM_GET returns the parameter values for problem 41.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    10 March 2002
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, double precision ALPHA, the current value of the parameter.
c
      implicit none

      double precision alpha

      call p41_param ( 'get', 'alpha', alpha )

      return
      end
      subroutine p41_param_set ( alpha )

c*********************************************************************72
c
cc P41_PARAM_SET sets the parameter values for problem 41.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    10 March 2002
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, double precision ALPHA, the new value of the parameter.
c
      implicit none

      double precision alpha

      call p41_param ( 'set', 'alpha', alpha )

      return
      end
      subroutine p42_exact ( exact )

c*********************************************************************72
c
cc P42_EXACT returns the exact integral for problem 42.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    11 March 2002
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, double precision EXACT, the value of the integral.
c
      implicit none

      double precision alpha
      double precision exact
      double precision r8_gamma

      call p42_param_get ( alpha )

      exact = 2.0D+00**( alpha - 2.0D+00 ) 
     &  * ( r8_gamma ( alpha / 2.0D+00 ) )**2 
     &  / r8_gamma ( alpha )

      return
      end
      subroutine p42_fun ( n, x, fx )

c*********************************************************************72
c
cc P42_FUN evaluates the integrand for problem 42.
c
c  Discussion:
c
c    The problem has a parameter ALPHA that can be set by calling
c    P42_PARAM_SET.
c
c    The integrand has a singularity at X = 0 if ALPHA < 1.
c
c    The suggested range for ALPHA is 0.1 through 2.
c
c  Interval:
c
c    0 .le. X .le. pi/2
c
c  Integrand:
c
c    ( sin(x) )^( alpha - 1 )
c
c  Exact Integral:
c
c    2^( alpha - 2 ) * ( Gamma(alpha/2) )^2 / Gamma(alpha)
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    14 November 2009
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Robert Piessens, Elise de Doncker-Kapenga,
c    Christian Ueberhuber, David Kahaner,
c    QUADPACK: A Subroutine Package for Automatic Integration,
c    Springer, 1983, page 83.
c
c  Parameters:
c
c    Input, integer N, the number of evaluation points.
c
c    Input, double precision X(N), the evaluation points.
c
c    Output, double precision FX(N), the integrand values.
c
      implicit none

      integer n

      double precision fx(n)
      integer i
      double precision x(n)
      double precision alpha
      double precision base

      call p42_param_get ( alpha )

      do i = 1, n

        base = sin ( x(i) )

        if ( base .eq. 0.0D+00 ) then

          if ( 1.0D+00 .lt. alpha ) then
            fx(i) = 0.0D+00
          else if ( alpha .eq. 1.0D+00 ) then
            fx(i) = 1.0D+00
          else
            fx(i) = 0.0D+00
          end if

        else

          fx(i) = base**( alpha - 1.0D+00 )

        end if

      end do

      return
      end
      subroutine p42_lim ( a, b )

c*********************************************************************72
c
cc P42_LIM returns the integration limits for problem 42.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    11 March 2002
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, double precision A, B, the limits of integration.
c
      implicit none

      double precision a
      double precision b
      double precision pi
      parameter ( pi = 3.141592653589793D+00 )

      a = 0.0D+00
      b = pi / 2.0D+00

      return
      end
      subroutine p42_param ( action, name, value )

c*********************************************************************72
c
cc P42_PARAM gets or sets the parameter values for problem 42.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    11 March 2002
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, character ( len = * ) ACTION.
c    'get' to get the value.
c    'set' to set the value.
c
c    Input, character ( len = * ) NAME, the name of the parameter.
c    'alpha' is the only option.
c
c    Input/output, double precision VALUE.
c    If the action is 'get', then VALUE returns the current parameter value.
c    If ACTION is 'set', then the parameter value is set to VALUE.
c
      implicit none

      character ( len = * ) action
      double precision alpha
      character ( len = * ) name
      logical s_eqi
      double precision value

      save alpha

      data alpha / 0.3D+00 /

      if ( s_eqi ( action, 'get' ) ) then
        if ( s_eqi ( name, 'alpha' ) ) then
          value = alpha
        else
          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) 'P42_PARAM - Fatal error!'
          write ( *, '(a)' ) '  Unrecognized name.'
          stop
        end if
      else if ( s_eqi ( action, 'set' ) ) then
        if ( s_eqi ( name, 'alpha' ) ) then
          alpha = value
        else
          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) 'P42_PARAM - Fatal error!'
          write ( *, '(a)' ) '  Unrecognized name.'
          stop
        end if
      else
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'P42_PARAM - Fatal error!'
        write ( *, '(a)' ) '  Unrecognized action.'
        stop
      end if

      return
      end
      subroutine p42_param_get ( alpha )

c*********************************************************************72
c
cc P42_PARAM_GET returns the parameter values for problem 42.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    11 March 2002
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, double precision ALPHA, the current value of the parameter.
c
      implicit none

      double precision alpha

      call p42_param ( 'get', 'alpha', alpha )

      return
      end
      subroutine p42_param_set ( alpha )

c*********************************************************************72
c
cc P42_PARAM_SET sets the parameter values for problem 42.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    11 March 2002
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, double precision ALPHA, the new value of the parameter.
c
      implicit none

      double precision alpha

      call p42_param ( 'set', 'alpha', alpha )

      return
      end
      subroutine p43_exact ( exact )

c*********************************************************************72
c
cc P43_EXACT returns the exact integral for problem 43.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    11 March 2002
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, double precision EXACT, the value of the integral.
c
      implicit none

      double precision alpha
      double precision exact
      double precision r8_gamma

      call p43_param_get ( alpha )

      exact = r8_gamma ( alpha )

      return
      end
      subroutine p43_fun ( n, x, fx )

c*********************************************************************72
c
cc P43_FUN evaluates the integrand for problem 43.
c
c  Discussion:
c
c    The problem has a parameter ALPHA that can be set by calling
c    P43_PARAM_SET.
c
c    The suggested parameter range is 0.1 .le. ALPHA .le. 2.0.
c
c    The integrand has an algebraic endpoint singularity at X = 1
c    times a singular factor.
c
c  Interval:
c
c    0 .le. X .le. 1
c
c  Integrand:
c
c    ( ln ( 1 / x ) )^( alpha - 1 )
c
c  Exact Integral:
c
c    Gamma(alpha)
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    11 March 2002
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Robert Piessens, Elise de Doncker-Kapenga,
c    Christian Ueberhuber, David Kahaner,
c    QUADPACK: A Subroutine Package for Automatic Integration,
c    Springer, 1983, page 84.
c
c  Parameters:
c
c    Input, integer N, the number of evaluation points.
c
c    Input, double precision X(N), the evaluation points.
c
c    Output, double precision FX(N), the integrand values.
c
      implicit none

      integer n

      double precision fx(n)
      integer i
      double precision x(n)
      double precision alpha

      call p43_param_get ( alpha )

      do i = 1, n

        if ( x(i) .le. 0.0D+00 ) then

          fx(i) = 0.0D+00

        else if ( x(i) .eq. 0.0D+00 ) then

          if ( alpha - 1.0D+00 .lt. 0.0D+00 ) then
            fx(i) = 0.0D+00
          else if ( alpha - 1.0D+00 .eq. 0.0D+00 ) then
            fx(i) = 1.0D+00
          else
            fx(i) = 0.0D+00
          end if

        else if ( x(i) .eq. 1.0D+00 ) then

          if ( alpha - 1.0D+00 .lt. 0.0D+00 ) then
            fx(i) = 0.0D+00
          else if ( alpha - 1.0D+00 .eq. 0.0D+00 ) then
            fx(i) = 1.0D+00
          else
            fx(i) = 0.0D+00
          end if

        else

          fx(i) = ( log ( 1.0D+00 / x(i) ) )**( alpha - 1.0D+00 )

        end if

      end do

      return
      end
      subroutine p43_lim ( a, b )

c*********************************************************************72
c
cc P43_LIM returns the integration limits for problem 43.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    11 March 2002
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, double precision A, B, the limits of integration.
c
      implicit none

      double precision a
      double precision b

      a = 0.0D+00
      b = 1.0D+00

      return
      end
      subroutine p43_param ( action, name, value )

c*********************************************************************72
c
cc P43_PARAM gets or sets the parameter values for problem 43.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    11 March 2002
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, character ( len = * ) ACTION.
c    'get' to get the value.
c    'set' to set the value.
c
c    Input, character ( len = * ) NAME, the name of the parameter.
c    'alpha' is the only option.
c
c    Input/output, double precision VALUE.
c    If the action is 'get', then VALUE returns the current parameter value.
c    If ACTION is 'set', then the parameter value is set to VALUE.
c
      implicit none

      character ( len = * ) action
      double precision alpha
      character ( len = * ) name
      logical s_eqi
      double precision value

      save alpha

      data alpha / 0.3D+00 /

      if ( s_eqi ( action, 'get' ) ) then
        if ( s_eqi ( name, 'alpha' ) ) then
          value = alpha
        else
          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) 'P43_PARAM - Fatal error!'
          write ( *, '(a)' ) '  Unrecognized name.'
          stop
        end if
      else if ( s_eqi ( action, 'set' ) ) then
        if ( s_eqi ( name, 'alpha' ) ) then
          alpha = value
        else
          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) 'P43_PARAM - Fatal error!'
          write ( *, '(a)' ) '  Unrecognized name.'
          stop
        end if
      else
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'P43_PARAM - Fatal error!'
        write ( *, '(a)' ) '  Unrecognized action.'
        stop
      end if

      return
      end
      subroutine p43_param_get ( alpha )

c*********************************************************************72
c
cc P43_PARAM_GET returns the parameter values for problem 43.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    11 March 2002
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, double precision ALPHA, the current value of the parameter.
c
      implicit none

      double precision alpha

      call p43_param ( 'get', 'alpha', alpha )

      return
      end
      subroutine p43_param_set ( alpha )

c*********************************************************************72
c
cc P43_PARAM_SET sets the parameter values for problem 43.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    11 March 2002
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, double precision ALPHA, the new value of the parameter.
c
      implicit none

      double precision alpha

      call p43_param ( 'set', 'alpha', alpha )

      return
      end
      subroutine p44_exact ( exact )

c*********************************************************************72
c
cc P44_EXACT returns the exact integral for problem 44.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    12 March 2002
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, double precision EXACT, the value of the integral.
c
      implicit none

      double precision alpha
      double precision exact

      call p44_param_get ( alpha )

      exact = ( 20.0D+00 * sin ( 2.0D+00**alpha ) 
     &  - 2.0D+00**alpha * cos ( 2.0D+00**alpha ) 
     &  + 2.0D+00**alpha * exp ( -20.0D+00 ) ) 
     &  / ( 400.0D+00 + 4.0D+00**alpha )

      return
      end
      subroutine p44_fun ( n, x, fx )

c*********************************************************************72
c
cc P44_FUN evaluates the integrand for problem 44.
c
c  Discussion:
c
c    The problem has a parameter ALPHA that can be set by calling
c    P44_PARAM_SET.
c
c    The suggested parameter range is 0.0 .le. ALPHA .le. 9.0.
c
c    As ALPHA increases, the integrand becomes more oscillatory.
c
c  Interval:
c
c    0 .le. X .le. 1
c
c  Integrand:
c
c    exp ( 20 * ( x - 1 ) ) * sin ( 2^alpha * x )
c
c  Exact Integral:
c
c    ( 20 sin ( 2^alpha ) - 2^alpha cos ( 2^alpha )
c    + 2^alpha exp ( -20 ) ) / ( 400 + 4^alpha )
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    12 March 2002
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Robert Piessens, Elise de Doncker-Kapenga,
c    Christian Ueberhuber, David Kahaner,
c    QUADPACK: A Subroutine Package for Automatic Integration,
c    Springer, 1983, page 84.
c
c  Parameters:
c
c    Input, integer N, the number of evaluation points.
c
c    Input, double precision X(N), the evaluation points.
c
c    Output, double precision FX(N), the integrand values.
c
      implicit none

      integer n

      double precision fx(n)
      integer i
      double precision x(n)
      double precision alpha

      call p44_param_get ( alpha )

      do i = 1, n
        fx(i) = exp ( 20.0D+00 * ( x(i) - 1.0D+00 ) ) 
     &    * sin ( 2.0D+00**alpha * x(i) )
      end do

      return
      end
      subroutine p44_lim ( a, b )

c*********************************************************************72
c
cc P44_LIM returns the integration limits for problem 44.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    12 March 2002
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, double precision A, B, the limits of integration.
c
      implicit none

      double precision a
      double precision b

      a = 0.0D+00
      b = 1.0D+00

      return
      end
      subroutine p44_param ( action, name, value )

c*********************************************************************72
c
cc P44_PARAM gets or sets the parameter values for problem 44.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    12 March 2002
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, character ( len = * ) ACTION.
c    'get' to get the value.
c    'set' to set the value.
c
c    Input, character ( len = * ) NAME, the name of the parameter.
c    'alpha' is the only option.
c
c    Input/output, double precision VALUE.
c    If the action is 'get', then VALUE returns the current parameter value.
c    If ACTION is 'set', then the parameter value is set to VALUE.
c
      implicit none

      character ( len = * ) action
      double precision alpha
      character ( len = * ) name
      logical s_eqi
      double precision value

      save alpha

      data alpha / 2.0D+00 /

      if ( s_eqi ( action, 'get' ) ) then
        if ( s_eqi ( name, 'alpha' ) ) then
          value = alpha
        else
          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) 'P44_PARAM - Fatal error!'
          write ( *, '(a)' ) '  Unrecognized name.'
          stop
        end if
      else if ( s_eqi ( action, 'set' ) ) then
        if ( s_eqi ( name, 'alpha' ) ) then
          alpha = value
        else
          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) 'P44_PARAM - Fatal error!'
          write ( *, '(a)' ) '  Unrecognized name.'
          stop
        end if
      else
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'P44_PARAM - Fatal error!'
        write ( *, '(a)' ) '  Unrecognized action.'
        stop
      end if

      return
      end
      subroutine p44_param_get ( alpha )

c*********************************************************************72
c
cc P44_PARAM_GET returns the parameter values for problem 44.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    12 March 2002
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, double precision ALPHA, the current value of the parameter.
c
      implicit none

      double precision alpha

      call p44_param ( 'get', 'alpha', alpha )

      return
      end
      subroutine p44_param_set ( alpha )

c*********************************************************************72
c
cc P44_PARAM_SET sets the parameter values for problem 44.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    12 March 2002
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, double precision ALPHA, the new value of the parameter.
c
      implicit none

      double precision alpha

      call p44_param ( 'set', 'alpha', alpha )

      return
      end
      subroutine p45_exact ( exact )

c*********************************************************************72
c
cc P45_EXACT returns the exact integral for problem 45.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    16 March 2002
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, double precision EXACT, the value of the integral.
c
      implicit none

      double precision alpha
      double precision besj0
      double precision exact
      double precision pi
      parameter ( pi = 3.141592653589793D+00 )

      call p45_param_get ( alpha )

      exact = pi * cos ( 2.0D+00**( alpha - 1.0D+00 ) ) 
     &  * besj0 ( 2.0D+00**( alpha - 1.0D+00 ) )

      return
      end
      subroutine p45_fun ( n, x, fx )

c*********************************************************************72
c
cc P45_FUN evaluates the integrand for problem 45.
c
c  Discussion:
c
c    The problem has a parameter ALPHA that can be set by calling
c    P45_PARAM_SET.
c
c    The suggested parameter range is 0.0 .le. ALPHA .le. 8.0.
c
c    The function is singular at 0 and 1.
c
c  Interval:
c
c    0 .le. X .le. 1
c
c  Integrand:
c
c    cos ( 2^alpha * x ) / sqrt ( x * ( 1 - x ) )
c
c  Exact Integral:
c
c    pi * cos ( 2^(alpha-1) ) * J0 ( 2^(alpha-1) )
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    16 March 2002
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Robert Piessens, Elise de Doncker-Kapenga,
c    Christian Ueberhuber, David Kahaner,
c    QUADPACK: A Subroutine Package for Automatic Integration,
c    Springer, 1983, page 84.
c
c  Parameters:
c
c    Input, integer N, the number of evaluation points.
c
c    Input, double precision X(N), the evaluation points.
c
c    Output, double precision FX(N), the integrand values.
c
      implicit none

      integer n

      double precision fx(n)
      integer i
      double precision x(n)
      double precision alpha

      call p45_param_get ( alpha )

      do i = 1, n

        if ( x(i) .eq. 0.0D+00 ) then
          fx(i) = 0.0D+00
        else if ( x(i) .eq. 1.0D+00 ) then
          fx(i) = 0.0D+00
        else
          fx(i) = cos ( 2.0D+00**alpha * x(i) ) 
     &      / sqrt ( x(i) * ( 1.0D+00 - x(i) ) )
        end if

      end do

      return
      end
      subroutine p45_lim ( a, b )

c*********************************************************************72
c
cc P45_LIM returns the integration limits for problem 45.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    16 March 2002
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, double precision A, B, the limits of integration.
c
      implicit none

      double precision a
      double precision b

      a = 0.0D+00
      b = 1.0D+00

      return
      end
      subroutine p45_param ( action, name, value )

c*********************************************************************72
c
cc P45_PARAM gets or sets the parameter values for problem 45.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    16 March 2002
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, character ( len = * ) ACTION.
c    'get' to get the value.
c    'set' to set the value.
c
c    Input, character ( len = * ) NAME, the name of the parameter.
c    'alpha' is the only option.
c
c    Input/output, double precision VALUE.
c    If the action is 'get', then VALUE returns the current parameter value.
c    If ACTION is 'set', then the parameter value is set to VALUE.
c
      implicit none

      character ( len = * ) action
      double precision alpha
      character ( len = * ) name
      logical s_eqi
      double precision value

      save alpha

      data alpha / 2.0D+00 /

      if ( s_eqi ( action, 'get' ) ) then
        if ( s_eqi ( name, 'alpha' ) ) then
          value = alpha
        else
          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) 'P45_PARAM - Fatal error!'
          write ( *, '(a)' ) '  Unrecognized name.'
          stop
        end if
      else if ( s_eqi ( action, 'set' ) ) then
        if ( s_eqi ( name, 'alpha' ) ) then
          alpha = value
        else
          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) 'P45_PARAM - Fatal error!'
          write ( *, '(a)' ) '  Unrecognized name.'
          stop
        end if
      else
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'P45_PARAM - Fatal error!'
        write ( *, '(a)' ) '  Unrecognized action.'
        stop
      end if

      return
      end
      subroutine p45_param_get ( alpha )

c*********************************************************************72
c
cc P45_PARAM_GET returns the parameter values for problem 45.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    16 March 2002
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, double precision ALPHA, the current value of the parameter.
c
      implicit none

      double precision alpha

      call p45_param ( 'get', 'alpha', alpha )

      return
      end
      subroutine p45_param_set ( alpha )

c*********************************************************************72
c
cc P45_PARAM_SET sets the parameter values for problem 45.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    16 March 2002
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, double precision ALPHA, the new value of the parameter.
c
      implicit none

      double precision alpha

      call p45_param ( 'set', 'alpha', alpha )

      return
      end
      subroutine p46_exact ( exact )

c*********************************************************************72
c
cc P46_EXACT returns the exact integral for problem 46.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    26 December 2011
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, double precision EXACT, the value of the integral.
c
      implicit none

      double precision exact

      exact = 6.0690909595647754101D+00

      return
      end
      subroutine p46_fun ( n, x, fx )

c*********************************************************************72
c
cc P46_FUN evaluates the integrand for problem 46.
c
c  Discussion:
c
c    The problem has a parameter ALPHA that can be set by calling
c    P63_PARAM_SET.
c
c    The integrand is the radius of an ellipse as a function of angle.
c
c    The integral represents the arc length of the ellipse.
c
c    The suggested parameter range is 0.0 .le. ALPHA < 1.0.  ALPHA is
c    the eccentricity of the ellipse.
c
c  Interval:
c
c    0 .le. theta .le. 2 pi
c
c  Integrand:
c
c    r(theta) = ( 1 - alpha^2 ) / ( 1 - alpha * cos ( theta ) )
c
c  Exact Integral:
c
c    When alpha = sin ( pi / 12 ), then
c
c      6.0690909595647754101
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    26 December 2011
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Richard Crandall,
c    Projects in Scientific Computing,
c    Springer, 2000, page 47.
c
c  Parameters:
c
c    Input, integer N, the number of evaluation points.
c
c    Input, double precision X(N), the evaluation points.
c
c    Output, double precision FX(N), the integrand values.
c
      implicit none

      integer n

      double precision fx(n)
      integer i
      double precision x(n)
      double precision alpha

      call p46_param_get ( alpha )

      do i = 1, n
        fx(i) = ( 1.0D+00 - alpha**2 ) 
     &    / ( 1.0D+00 - alpha * cos ( x(i) ) )
      end do

      return
      end
      subroutine p46_lim ( a, b )

c*********************************************************************72
c
cc P46_LIM returns the integration limits for problem 46.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    20 March 2002
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, double precision A, B, the limits of integration.
c
      implicit none

      double precision a
      double precision b
      double precision pi
      parameter ( pi = 3.141592653589793D+00 )

      a = 0.0D+00
      b = 2.0D+00 * pi

      return
      end
      subroutine p46_param ( action, name, value )

c*********************************************************************72
c
cc P46_PARAM gets or sets the parameter values for problem 46.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    25 December 2011
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, character ( len = * ) ACTION.
c    'get' to get the value.
c    'set' to set the value.
c
c    Input, character ( len = * ) NAME, the name of the parameter.
c    'alpha' is the only option.
c
c    Input/output, double precision VALUE.
c    If the action is 'get', then VALUE returns the current parameter value.
c    If ACTION is 'set', then the parameter value is set to VALUE.
c
      implicit none

      character ( len = * ) action
      double precision alpha
      character ( len = * ) name
      double precision pi
      parameter ( pi = 3.141592653589793D+00 )
      logical s_eqi
      logical set
      double precision value

      save alpha
      save set

      data alpha / 0.0D+00 /
      data set / .false. /

      if ( .not. set ) then
        alpha = sin ( pi / 12.0D+00 )
        set = .true.
      end if

      if ( s_eqi ( action, 'get' ) ) then
        if ( s_eqi ( name, 'alpha' ) ) then
          value = alpha
        else
          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) 'P46_PARAM - Fatal error!'
          write ( *, '(a)' ) '  Unrecognized name.'
          stop
        end if
      else if ( s_eqi ( action, 'set' ) ) then
        if ( s_eqi ( name, 'alpha' ) ) then
          alpha = value
        else
          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) 'P46_PARAM - Fatal error!'
          write ( *, '(a)' ) '  Unrecognized name.'
          stop
        end if
      else
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'P46_PARAM - Fatal error!'
        write ( *, '(a)' ) '  Unrecognized action.'
        stop
      end if

      return
      end
      subroutine p46_param_get ( alpha )

c*********************************************************************72
c
cc P46_PARAM_GET returns the parameter values for problem 46.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    20 March 2002
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, double precision ALPHA, the current value of the parameter.
c
      implicit none

      double precision alpha

      call p46_param ( 'get', 'alpha', alpha )

      return
      end
      subroutine p46_param_set ( alpha )

c*********************************************************************72
c
cc P46_PARAM_SET sets the parameter values for problem 46.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    20 March 2002
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, double precision ALPHA, the new value of the parameter.
c
      implicit none

      double precision alpha

      call p46_param ( 'set', 'alpha', alpha )

      return
      end
      subroutine p47_exact ( exact )

c*********************************************************************72
c
cc P47_EXACT returns the exact integral for problem 47.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    21 March 2002
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, double precision EXACT, the value of the integral.
c
      implicit none

      double precision exact

      exact = - 4.0D+00 / 9.0D+00

      return
      end
      subroutine p47_fun ( n, x, fx )

c*********************************************************************72
c
cc P47_FUN evaluates the integrand for problem 47.
c
c  Discussion:
c
c    The function is singular at the left endpoint.
c
c  Interval:
c
c    0 .le. x .le. 1
c
c  Integrand:
c
c    sqrt ( x ) * ln ( x )
c
c  Exact Integral:
c
c    -4/9 = -0.4444...
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    14 November 2009
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Robert Piessens, Elise de Doncker-Kapenga,
c    Christian Ueberhuber, David Kahaner,
c    QUADPACK: A Subroutine Package for Automatic Integration,
c    Springer, 1983, page 101.
c
c  Parameters:
c
c    Input, integer N, the number of evaluation points.
c
c    Input, double precision X(N), the evaluation points.
c
c    Output, double precision FX(N), the integrand values.
c
      implicit none

      integer n

      double precision fx(n)
      integer i
      double precision x(n)

      do i = 1, n

        if ( x(i) .eq. 0.0D+00 ) then
          fx(i) = 0.0D+00
        else
          fx(i) = sqrt ( x(i) ) * log ( x(i) )
        end if

      end do

      return
      end
      subroutine p47_lim ( a, b )

c*********************************************************************72
c
cc P47_LIM returns the integration limits for problem 47.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    21 March 2002
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, double precision A, B, the limits of integration.
c
      implicit none

      double precision a
      double precision b

      a = 0.0D+00
      b = 1.0D+00

      return
      end
      subroutine p48_exact ( exact )

c*********************************************************************72
c
cc P48_EXACT returns the exact integral for problem 48.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    21 March 2002
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, double precision EXACT, the value of the integral.
c
      implicit none

      double precision exact

      exact = -4.0D+00

      return
      end
      subroutine p48_fun ( n, x, fx )

c*********************************************************************72
c
cc P48_FUN evaluates the integrand for problem 48.
c
c  Discussion:
c
c    The function is singular at the left endpoint.
c
c  Interval:
c
c    0 .le. x .le. 1
c
c  Integrand:
c
c    ln ( x ) / sqrt ( x )
c
c  Exact Integral:
c
c    -4
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    14 November 2009
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Robert Piessens, Elise de Doncker-Kapenga,
c    Christian Ueberhuber, David Kahaner,
c    QUADPACK: A Subroutine Package for Automatic Integration,
c    Springer, 1983, page 103.
c
c  Parameters:
c
c    Input, integer N, the number of evaluation points.
c
c    Input, double precision X(N), the evaluation points.
c
c    Output, double precision FX(N), the integrand values.
c
      implicit none

      integer n

      double precision fx(n)
      integer i
      double precision x(n)

      do i = 1, n

        if ( x(i) .eq. 0.0D+00 ) then
          fx(i) = 0.0D+00
        else
          fx(i) = log ( x(i) ) / sqrt ( x(i) )
        end if

      end do

      return
      end
      subroutine p48_lim ( a, b )

c*********************************************************************72
c
cc P48_LIM returns the integration limits for problem 48.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    21 March 2002
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, double precision A, B, the limits of integration.
c
      implicit none

      double precision a
      double precision b

      a = 0.0D+00
      b = 1.0D+00

      return
      end
      subroutine p49_exact ( exact )

c*********************************************************************72
c
cc P49_EXACT returns the exact integral for problem 49.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    22 March 2002
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, double precision EXACT, the value of the integral.
c
      implicit none

      double precision exact

      exact = 61.0D+00 * log ( 2.0D+00 ) 
     &  + 77.0D+00 * log ( 7.0D+00 ) / 4.0D+00 - 27.0D+00

      return
      end
      subroutine p49_fun ( n, x, fx )

c*********************************************************************72
c
cc P49_FUN evaluates the integrand for problem 49.
c
c  Discussion:
c
c    The function is singular at two internal points, 1 and sqrt(2).
c
c  Interval:
c
c    0 .le. x .le. 3
c
c  Integrand:
c
c    x^3 * log ( abs ( ( x^2 - 1 ) * ( x^2 - 2 ) ) )
c
c  Exact Integral:
c
c    61 log ( 2 ) + (77/4) log ( 7 ) - 27 = 52.7408...
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    14 November 2009
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Robert Piessens, Elise de Doncker-Kapenga,
c    Christian Ueberhuber, David Kahaner,
c    QUADPACK: A Subroutine Package for Automatic Integration,
c    Springer, 1983, page 104.
c
c  Parameters:
c
c    Input, integer N, the number of evaluation points.
c
c    Input, double precision X(N), the evaluation points.
c
c    Output, double precision FX(N), the integrand values.
c
      implicit none

      integer n

      double precision fx(n)
      integer i
      double precision x(n)

      do i = 1, n

        if ( ( x(i)**2 - 1.0D+00 ) 
     &    * ( x(i)**2 - 2.0D+00 ) .eq. 0.0D+00 ) then
          fx(i) = 0.0D+00
        else
          fx(i) = x(i)**3 * log ( abs ( ( x(i)**2 - 1.0D+00 ) 
     &      * ( x(i)**2 - 2.0D+00 ) ) )
        end if

      end do

      return
      end
      subroutine p49_lim ( a, b )

c*********************************************************************72
c
cc P49_LIM returns the integration limits for problem 49.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    22 March 2002
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, double precision A, B, the limits of integration.
c
      implicit none

      double precision a
      double precision b

      a = 0.0D+00
      b = 3.0D+00

      return
      end
      subroutine p50_exact ( exact )

c*********************************************************************72
c
cc P50_EXACT returns the exact integral for problem 50.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    23 March 2002
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, double precision EXACT, the value of the integral.
c
      implicit none

      double precision exact
      double precision euler_constant
      double precision pi
      parameter ( pi = 3.141592653589793D+00 )
      double precision r8_ci
      double precision t

      t = 10.0D+00 * pi

      exact = ( - euler_constant ( ) - log ( t ) + r8_ci ( t ) ) / t

      return
      end
      subroutine p50_fun ( n, x, fx )

c*********************************************************************72
c
cc P50_FUN evaluates the integrand for problem 50.
c
c  Discussion:
c
c    The function has a removable singularity at x = 0.
c
c  Interval:
c
c    0 .le. x .le. 1
c
c  Integrand:
c
c    log ( x ) * sin ( 10 * pi * x )
c
c  Exact Integral:
c
c    ( - gamma - log ( 10 * pi ) + Ci ( 10 * pi ) ) / 10 * pi = -0.1281316...
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    23 March 2002
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Robert Piessens, Elise de Doncker-Kapenga,
c    Christian Ueberhuber, David Kahaner,
c    QUADPACK: A Subroutine Package for Automatic Integration,
c    Springer, 1983, page 106.
c
c  Parameters:
c
c    Input, integer N, the number of evaluation points.
c
c    Input, double precision X(N), the evaluation points.
c
c    Output, double precision FX(N), the integrand values.
c
      implicit none

      integer n

      double precision fx(n)
      integer i
      double precision x(n)
      double precision pi
      parameter ( pi = 3.141592653589793D+00 )

      do i = 1, n

        if ( x(i) .eq. 0.0D+00 ) then
          fx(i) = 0.0D+00
        else
          fx(i) = log ( x(i) ) * sin ( 10.0D+00 * pi * x(i) )
        end if

      end do

      return
      end
      subroutine p50_lim ( a, b )

c*********************************************************************72
c
cc P50_LIM returns the integration limits for problem 50.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    23 March 2002
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, double precision A, B, the limits of integration.
c
      implicit none

      double precision a
      double precision b

      a = 0.0D+00
      b = 1.0D+00

      return
      end
      subroutine p51_exact ( exact )

c*********************************************************************72
c
cc P51_EXACT returns the exact integral for problem 51.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    12 November 2009
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, double precision EXACT, the value of the integral.
c
      implicit none

      double precision exact
      double precision pi
      parameter ( pi = 3.141592653589793D+00 )
      double precision r8_ci
      double precision r8_si

      exact = - ( r8_ci ( 1.0D+00 ) * sin ( 1.0D+00 ) + 
     &  ( 0.5D+00 * pi - r8_si ( 1.0D+00 ) ) * cos ( 1.0D+00 ) ) / pi

      return
      end
      subroutine p51_fun ( n, x, fx )

c*********************************************************************72
c
cc P51_FUN evaluates the integrand for problem 51.
c
c  Discussion:
c
c    The function has a removable singularity at x = 0.
c
c  Interval:
c
c    0 .le. x .le. 1
c
c  Integrand:
c
c    ln ( x ) / ( 1 + ( ln(x) )^2 )^2
c
c  Exact Integral:
c
c    - ( ci(1) * sin(1) + ( pi/2 - si(1) ) * cos(1) ) / pi = - 0.1892752...
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    24 March 2002
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Robert Piessens, Elise de Doncker-Kapenga,
c    Christian Ueberhuber, David Kahaner,
c    QUADPACK: A Subroutine Package for Automatic Integration,
c    Springer, 1983, page 108.
c
c  Parameters:
c
c    Input, integer N, the number of evaluation points.
c
c    Input, double precision X(N), the evaluation points.
c
c    Output, double precision FX(N), the integrand values.
c
      implicit none

      integer n

      double precision fx(n)
      integer i
      double precision x(n)

      do i = 1, n

        if ( x(i) .eq. 0.0D+00 ) then
          fx(i) = 0.0D+00
        else
          fx(i) = log ( x(i) ) / ( 1.0D+00 + ( log ( x(i) ) )**2 )**2
        end if

      end do

      return
      end
      subroutine p51_lim ( a, b )

c*********************************************************************72
c
cc P51_LIM returns the integration limits for problem 51.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    24 March 2002
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, double precision A, B, the limits of integration.
c
      implicit none

      double precision a
      double precision b

      a = 0.0D+00
      b = 1.0D+00

      return
      end
      subroutine p52_exact ( exact )

c*********************************************************************72
c
cc P52_EXACT returns the exact integral for problem 52.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    24 March 2002
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, double precision EXACT, the value of the integral.
c
      implicit none

      double precision exact

      exact = log ( 125.0D+00 / 631.0D+00 ) / 18.0D+00

      return
      end
      subroutine p52_fun ( n, x, fx )

c*********************************************************************72
c
cc P52_FUN evaluates the integrand for problem 52.
c
c  Discussion:
c
c    The function has a singularity at x = 0.
c
c  Interval:
c
c    -1 .le. x .le. 5
c
c  Integrand:
c
c    1 / ( x * ( 5 * x^3 + 6 ) )
c
c  Exact Integral:
c
c    ln ( 125 / 631 ) / 18 = -0.08994401...
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    24 March 2002
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Robert Piessens, Elise de Doncker-Kapenga,
c    Christian Ueberhuber, David Kahaner,
c    QUADPACK: A Subroutine Package for Automatic Integration,
c    Springer, 1983, page 109.
c
c  Parameters:
c
c    Input, integer N, the number of evaluation points.
c
c    Input, double precision X(N), the evaluation points.
c
c    Output, double precision FX(N), the integrand values.
c
      implicit none

      integer n

      double precision fx(n)
      integer i
      double precision x(n)

      do i = 1, n

        if ( x(i) .eq. 0.0D+00 ) then
          fx(i) = 0.0D+00
        else
          fx(i) = 1.0D+00 / ( x(i) * ( 5.0D+00 * x(i)**3 + 6.0D+00 ) )
        end if

      end do

      return
      end
      subroutine p52_lim ( a, b )

c*********************************************************************72
c
cc P52_LIM returns the integration limits for problem 52.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    24 March 2002
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, double precision A, B, the limits of integration.
c
      implicit none

      double precision a
      double precision b

      a = -1.0D+00
      b = 5.0D+00

      return
      end
      subroutine p53_exact ( exact )

c*********************************************************************72
c
cc P53_EXACT returns the exact integral for problem 53.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    26 December 2011
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, double precision EXACT, the value of the integral.
c
      implicit none

      double precision exact
      double precision pi
      parameter ( pi = 3.141592653589793D+00 )

      exact = 0.5D+00 * pi - atan ( 1.0D+00 / sqrt ( 2.0D+00 ) ) 
     &  + log ( 3.0D+00 ) / 2.0D+00

      return
      end
      subroutine p53_fun ( n, x, fx )

c*********************************************************************72
c
cc P53_FUN evaluates the integrand for problem 53.
c
c  Discussion:
c
c    The integrand is singular at x = -1 + sqrt ( 3 ) = 0.732...
c
c  Interval:
c
c    0 .le. x .le. 1
c
c  Integrand:
c
c    1 / sqrt ( abs ( x^2 + 2 * x - 2 ) )
c
c  Exact Integral:
c
c    pi / 2 - arctan ( 1 / sqrt ( 2 ) ) + ln ( 3 ) / 2 = 1.504622...
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    26 December 2011
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Robert Piessens, Elise de Doncker-Kapenga,
c    Christian Ueberhuber, David Kahaner,
c    QUADPACK: A Subroutine Package for Automatic Integration,
c    Springer, 1983, page 110.
c
c  Parameters:
c
c    Input, integer N, the number of evaluation points.
c
c    Input, double precision X(N), the evaluation points.
c
c    Output, double precision FX(N), the integrand values.
c
      implicit none

      integer n

      double precision fx(n)
      integer i
      double precision x(n)

      do i = 1, n
        fx(i) = 1.0D+00 
     &    / sqrt ( abs ( x(i)**2 + 2.0D+00 * x(i) - 2.0D+00 ) )
      end do

      return
      end
      subroutine p53_lim ( a, b )

c*********************************************************************72
c
cc P53_LIM returns the integration limits for problem 53.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    24 March 2002
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, double precision A, B, the limits of integration.
c
      implicit none

      double precision a
      double precision b

      a = 0.0D+00
      b = 1.0D+00

      return
      end
      subroutine p54_exact ( exact )

c*********************************************************************72
c
cc P54_EXACT returns the exact integral for problem 54.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    01 October 2002
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, double precision EXACT, the value of the integral.
c
      implicit none

      double precision exact

      exact = 2.0D+00 / sqrt ( 3.0D+00 )

      return
      end
      subroutine p54_fun ( n, x, fx )

c*********************************************************************72
c
cc P54_FUN evaluates the integrand for problem 54.
c
c  Discussion:
c
c    The reference claims that this integrand is more closely approximated
c    by the trapezoid rule than by Gauss-Legendre quadrature.
c
c    Points  Trapezoid  Gauss-Legendre
c     4      1.91667    2.53883
c    12      2.1594     2.25809
c
c    However, the stated results hardly give one confidence in
c    the convergence of the trapezoid results, and I am unable to
c    confirm them, because my results for 4 points give good results
c    (about 1.14) for BOTH Trapezoid and Gauss-Legendre!
c
c  Interval:
c
c    0 .le. x .le. 1
c
c  Integrand:
c
c    2 / ( 2 + sin ( 10 * PI * X ) )
c
c  Exact Integral:
c
c    2 / sqrt ( 3 ) = 1.1547...
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    01 October 2002
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Prem Kythe, Pratap Puri,
c    Computational Methods for Linear Integral Equations,
c    Birkhaeuser, 2002.
c
c  Parameters:
c
c    Input, integer N, the number of evaluation points.
c
c    Input, double precision X(N), the evaluation points.
c
c    Output, double precision FX(N), the integrand values.
c
      implicit none

      integer n

      double precision fx(n)
      integer i
      double precision x(n)
      double precision pi
      parameter ( pi = 3.141592653589793D+00 )

      do i = 1, n
        fx(i) = 2.0D+00 / ( 2.0D+00 + sin ( 10.0D+00 * pi * x(i) ) )
      end do

      return
      end
      subroutine p54_lim ( a, b )

c*********************************************************************72
c
cc P54_LIM returns the integration limits for problem 54.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    01 October 2002
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, double precision A, B, the limits of integration.
c
      implicit none

      double precision a
      double precision b

      a = 0.0D+00
      b = 1.0D+00

      return
      end
      subroutine p55_exact ( exact )

c*********************************************************************72
c
cc P55_EXACT returns the exact integral for problem 55.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    21 April 2003
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, double precision EXACT, the value of the integral.
c
      implicit none

      double precision a
      double precision b
      double precision c
      double precision exact
      double precision pi
      parameter ( pi = 3.141592653589793D+00 )
      double precision r8_erf
      double precision x0

      call p55_lim ( a, b )
      call p55_param ( 'get', 'c', c )
      call p55_param ( 'get', 'x0', x0 )

      exact = sqrt ( pi ) * 
     &  ( r8_erf ( c * ( b - x0 ) ) - r8_erf ( c * ( a - x0 ) ) ) 
     &  / ( 2.0D+00 * c )

      return
      end
      subroutine p55_fun ( n, x, fx )

c*********************************************************************72
c
cc P55_FUN evaluates the integrand for problem 55.
c
c  Interval:
c
c    0 .le. x .le. 1
c
c  Integrand:
c
c    exp ( - c^2 * ( x - x0 )^2 )
c
c  Exact Integral:
c
c    sqrt ( pi )
c    * ( erf ( c * ( b - x0 ) ) - erf ( c * ( a - x0 ) ) )
c    / ( 2 * c )
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    20 April 2003
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the number of evaluation points.
c
c    Input, double precision X(N), the evaluation points.
c
c    Output, double precision FX(N), the integrand values.
c
      implicit none

      integer n

      double precision fx(n)
      integer i
      double precision x(n)
      double precision c
      double precision x0

      call p55_param ( 'get', 'c', c )
      call p55_param ( 'get', 'x0', x0 )

      do i = 1, n
        fx(i) = exp ( - c**2 * ( x(i) - x0 )**2 )
      end do

      return
      end
      subroutine p55_lim ( a, b )

c*********************************************************************72
c
cc P55_LIM returns the integration limits for problem 55.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    20 April 2003
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, double precision A, B, the limits of integration.
c
      implicit none

      double precision a
      double precision b

      a = 0.0D+00
      b = 1.0D+00

      return
      end
      subroutine p55_param ( action, name, value )

c*********************************************************************72
c
cc P55_PARAM sets or gets real scalar parameters for problem 55.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    20 April 2003
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, character ( len = * ) ACTION,
c    'get' to get a parameter.
c    'set' to set a parameter.
c
c    Input, character ( len = * ) NAME, the name of the variable.
c    'C' is the coefficient.
c    'X0' is the base point.
c
c    Input/output, double precision VALUE.
c    * If ACTION = 'set', then VALUE is an input quantity, and is the
c      new value to be assigned to NAME.
c    * If ACTION = 'get', then VALUE is an output quantity, and is the
c      current value of NAME.
c
      implicit none

      character ( len = * ) action
      double precision c
      character ( len = * ) name
      double precision value
      double precision, save :: x0 = 0.75D+00

      save c

      data c / 3.0D+00 /

      if ( action(1:1) .eq. 'G' .or. action(1:1) .eq. 'g' ) then

        if ( name .eq. 'C' .or. name .eq. 'c' ) then
          value = c
        else if ( name .eq. 'X0' .or. name .eq. 'x0' ) then
          value = x0
        else
          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) 'P55_PARAM - Fatal error!'
          write ( *, '(a)' ) 
     &      '  Unrecognized name = "' // trim ( name ) // '".'
          stop
        end if

      else if ( action(1:1) .eq. 'S' .or. action(1:1) .eq. 's' ) then

        if ( name .eq. 'C' .or. name .eq. 'c' ) then
          c = value
        else if ( name .eq. 'X0' .or. name .eq. 'x0' ) then
          x0 = value
        else
          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) 'P55_PARAM - Fatal error!'
          write ( *, '(a)' ) 
     &      '  Unrecognized name = "' // trim ( name ) // '".'
          stop
        end if

      else

        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'P55_PARAM - Fatal error!'
        write ( *, '(a)' ) 
     &    '  Unrecognized action = "' // trim ( action ) // '".'
        stop

      end if

      return
      end
      subroutine p56_exact ( exact )

c*********************************************************************72
c
cc P56_EXACT returns the estimated integral for problem 56.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    31 October 2000
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, double precision EXACT, the estimated value of the integral.
c
      implicit none

      double precision exact

      exact = 1.9922524079504000171D+00

      return
      end
      subroutine p56_fun ( n, x, fx )

c*********************************************************************72
c
cc P56_FUN evaluates the integrand for problem 56.
c
c  Interval:
c
c    -1 .le. x .le. 1
c
c  Integrand:
c
c    1 / ( x^6 + 0.9 )
c
c  Approximate Integral (20 digits):
c
c    1.9922524079504000171...
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    31 October 2000
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    David Kahaner,
c    Comparison of Numerical Quadrature Formulas,
c    in Mathematical Software,
c    edited by John R Rice,
c    Academic Press, 1971.
c
c  Parameters:
c
c    Input, integer N, the number of evaluation points.
c
c    Input, double precision X(N), the evaluation points.
c
c    Output, double precision FX(N), the integrand values.
c
      implicit none

      integer n

      double precision fx(n)
      integer i
      double precision x(n)

      do i = 1, n
        fx(i) = 1.0D+00 / ( x(i)**6 + 0.9D+00 )
      end do

      return
      end
      subroutine p56_lim ( a, b )

c*********************************************************************72
c
cc P56_LIM returns the integration limits for problem 56.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    22 December 1998
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, double precision A, B, the limits of integration.
c
      implicit none

      double precision a
      double precision b

      a = - 1.0D+00
      b = 1.0D+00

      return
      end
      subroutine p57_exact ( exact )

c*********************************************************************72
c
cc P57_EXACT returns the exact integral for problem 57.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    18 December 1998
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, double precision EXACT, the estimated value of the integral.
c
      implicit none

      double precision exact

      exact = 0.4D+00

      return
      end
      subroutine p57_fun ( n, x, fx )

c*********************************************************************72
c
cc P57_FUN evaluates the integrand for problem 57.
c
c  Interval:
c
c    0 .le. x .le. 1
c
c  Integrand:
c
c    x^(3/2)
c
c  Antiderivative:
c
c    (2/5) * x^(5/2)
c
c  Exact Integral:
c
c    0.4
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    31 October 2000
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    David Kahaner,
c    Comparison of Numerical Quadrature Formulas,
c    in Mathematical Software, edited by John R Rice,
c    Academic Press, 1971.
c
c  Parameters:
c
c    Input, integer N, the number of evaluation points.
c
c    Input, double precision X(N), the evaluation points.
c
c    Output, double precision FX(N), the integrand values.
c
      implicit none

      integer n

      double precision fx(n)
      integer i
      double precision x(n)

      do i = 1, n
        fx(i) = sqrt ( x(i)**3 )
      end do

      return
      end
      subroutine p57_lim ( a, b )

c*********************************************************************72
c
cc P57_LIM returns the integration limits for problem 57.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    22 December 1998
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, double precision A, B, the limits of integration.
c
      implicit none

      double precision a
      double precision b

      a = 0.0D+00
      b = 1.0D+00

      return
      end
      function r8_ci ( x )

c*********************************************************************72
c
cc R8_CI evaluates the cosine integral Ci of an R8 argument.
c
c  Discussion:
c
c    The cosine integral is defined by
c
c      CI(X) = - integral ( X .le. T < +oo ) ( cos ( T ) ) / T  dT
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    25 March 2010
c
c  Author:
c
c    Original FORTRAN77 version by Wayne Fullerton.
c    This FORTRAN77 version by John Burkardt.
c
c  Reference:
c
c    Wayne Fullerton,
c    Portable Special Function Routines,
c    in Portability of Numerical Software,
c    edited by Wayne Cowell,
c    Lecture Notes in Computer Science, Volume 57,
c    Springer 1977,
c    ISBN: 978-3-540-08446-4,
c    LC: QA297.W65.
c
c  Parameters:
c
c    Input, double precision X, the argument.
c
c    Output, double precision R8_CI, the cosine integral Ci evaluated at X.
c
      implicit none

      double precision cics(19)
      double precision f
      double precision g
      integer nci
      double precision r8_ci
      double precision r8_csevl
      integer r8_inits
      double precision r8_mach
      double precision sinx
      double precision x
      double precision xsml
      double precision y

      save cics
      save nci
      save xsml

      data cics(  1) / -0.34004281856055363156281076633129873D+00 /
      data cics(  2) / -1.03302166401177456807159271040163751D+00 /
      data cics(  3) /  0.19388222659917082876715874606081709D+00 /
      data cics(  4) / -0.01918260436019865893946346270175301D+00 /
      data cics(  5) /  0.00110789252584784967184098099266118D+00 /
      data cics(  6) / -0.00004157234558247208803840231814601D+00 /
      data cics(  7) /  0.00000109278524300228715295578966285D+00 /
      data cics(  8) / -0.00000002123285954183465219601280329D+00 /
      data cics(  9) /  0.00000000031733482164348544865129873D+00 /
      data cics( 10) / -0.00000000000376141547987683699381798D+00 /
      data cics( 11) /  0.00000000000003622653488483964336956D+00 /
      data cics( 12) / -0.00000000000000028911528493651852433D+00 /
      data cics( 13) /  0.00000000000000000194327860676494420D+00 /
      data cics( 14) / -0.00000000000000000001115183182650184D+00 /
      data cics( 15) /  0.00000000000000000000005527858887706D+00 /
      data cics( 16) / -0.00000000000000000000000023907013943D+00 /
      data cics( 17) /  0.00000000000000000000000000091001612D+00 /
      data cics( 18) / -0.00000000000000000000000000000307233D+00 /
      data cics( 19) /  0.00000000000000000000000000000000926D+00 /

      data nci / 0 /
      data xsml / 0.0D+00 /

      if ( nci .eq. 0 ) then
        nci = r8_inits ( cics, 19, 0.1D+00 * r8_mach ( 3 ) )
        xsml = dsqrt ( r8_mach ( 3 ) )
      end if

      if ( x .le. 0.0D+00 ) then
        
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'R8_CI - Fatal error!'
        write ( *, '(a)' ) '  X .le. 0.0.'
        stop
      
      else if ( x .le. xsml ) then
        y = - 1.0D+00
        r8_ci = dlog ( x ) - 0.5D+00 + r8_csevl ( y, cics, nci )
      else if ( x .le. 4.0D+00 ) then
        y = ( x * x - 8.0D+00 ) * 0.125D+00
        r8_ci = dlog ( x ) - 0.5D+00 + r8_csevl ( y, cics, nci )
      else
        call r8_sifg ( x, f, g )
        sinx = dsin ( x )
        r8_ci = f * sinx - g * dcos ( x )
      end if

      return
      end
      function r8_csevl ( x, a, n )

c*********************************************************************72
c
cc R8_CSEVL evaluates a Chebyshev series.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    09 March 2010
c
c  Author:
c
c    This FORTRAN77 version by John Burkardt.
c
c  Reference:
c
c    Roger Broucke,
c    Algorithm 446:
c    Ten Subroutines for the Manipulation of Chebyshev Series,
c    Communications of the ACM,
c    Volume 16, Number 4, April 1973, pages 254-256.
c
c  Parameters:
c
c    Input, double precision X, the evaluation point.
c
c    Input, double precision CS(N), the Chebyshev coefficients.
c
c    Input, integer N, the number of Chebyshev coefficients.
c
c    Output, double precision R8_CSEVL, the Chebyshev series evaluated at X.
c
      implicit none

      integer n

      double precision a(n)
      double precision b0
      double precision b1
      double precision b2
      integer i
      double precision r8_csevl
      double precision twox
      double precision x

      if ( n .lt. 1 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'R8_CSEVL - Fatal error!'
        write ( *, '(a)' ) '  Number of terms .le. 0.'
        stop
      end if

      if ( 1000 .lt. n ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'R8_CSEVL - Fatal error!'
        write ( *, '(a)' ) '  Number of terms > 1000.'
        stop
      end if

      if ( x .lt. -1.1D+00 .or. 1.1D+00 .lt. x ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'R8_CSEVL - Fatal error!'
        write ( *, '(a)' ) '  X outside (-1,+1)'
        write ( *, '(a,g14.6)' ) '  X = ', x
        stop
      end if

      twox = 2.0D+00 * x
      b1 = 0.0D+00
      b0 = 0.0D+00

      do i = n, 1, -1
        b2 = b1
        b1 = b0
        b0 = twox * b1 - b2 + a(i)
      end do

      r8_csevl = 0.5D+00 * ( b0 - b2 )

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

      double precision r8_epsilon

      r8_epsilon = 2.220446049250313D-016

      return
      end
      function r8_erf ( x )

c*********************************************************************72
c
cc R8_ERF evaluates the error function of an R8 argument.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    12 March 2010
c
c  Author:
c
c    Original FORTRAN77 version by Wayne Fullerton.
c    This FORTRAN77 version by John Burkardt.
c
c  Reference:
c
c    Wayne Fullerton,
c    Portable Special Function Routines,
c    in Portability of Numerical Software,
c    edited by Wayne Cowell,
c    Lecture Notes in Computer Science, Volume 57,
c    Springer 1977,
c    ISBN: 978-3-540-08446-4,
c    LC: QA297.W65.
c
c  Parameters:
c
c    Input, double precision X, the argument.
c
c    Output, double precision R8_ERF, the error function of X.
c
      implicit none

      double precision erfcs(21)
      integer nterf
      double precision r8_csevl
      double precision r8_erf
      double precision r8_erfc
      integer r8_inits
      double precision r8_mach
      double precision sqeps
      double precision sqrtpi
      double precision value
      double precision x
      double precision xbig
      double precision y

      save erfcs
      save nterf
      save sqeps
      save sqrtpi
      save xbig

      data erfcs(  1) / -0.49046121234691808039984544033376D-01 /
      data erfcs(  2) / -0.14226120510371364237824741899631D+00 /
      data erfcs(  3) / +0.10035582187599795575754676712933D-01 /
      data erfcs(  4) / -0.57687646997674847650827025509167D-03 /
      data erfcs(  5) / +0.27419931252196061034422160791471D-04 /
      data erfcs(  6) / -0.11043175507344507604135381295905D-05 /
      data erfcs(  7) / +0.38488755420345036949961311498174D-07 /
      data erfcs(  8) / -0.11808582533875466969631751801581D-08 /
      data erfcs(  9) / +0.32334215826050909646402930953354D-10 /
      data erfcs( 10) / -0.79910159470045487581607374708595D-12 /
      data erfcs( 11) / +0.17990725113961455611967245486634D-13 /
      data erfcs( 12) / -0.37186354878186926382316828209493D-15 /
      data erfcs( 13) / +0.71035990037142529711689908394666D-17 /
      data erfcs( 14) / -0.12612455119155225832495424853333D-18 /
      data erfcs( 15) / +0.20916406941769294369170500266666D-20 /
      data erfcs( 16) / -0.32539731029314072982364160000000D-22 /
      data erfcs( 17) / +0.47668672097976748332373333333333D-24 /
      data erfcs( 18) / -0.65980120782851343155199999999999D-26 /
      data erfcs( 19) / +0.86550114699637626197333333333333D-28 /
      data erfcs( 20) / -0.10788925177498064213333333333333D-29 /
      data erfcs( 21) / +0.12811883993017002666666666666666D-31 /

      data nterf / 0 /
      data sqeps / 0.0D+00 /
      data sqrtpi / 1.77245385090551602729816748334115D+00 /
      data xbig / 0.0D+00 /

      if ( nterf .eq. 0 ) then
        nterf = r8_inits ( erfcs, 21, 0.1D+00 * r8_mach ( 3 ) )
        xbig = dsqrt ( - dlog ( sqrtpi * r8_mach ( 3 ) ) )
        sqeps = dsqrt ( 2.0D+00 * r8_mach ( 3 ) )
      end if

      y = dabs ( x )

      if ( y .le. sqeps ) then
        value = 2.0D+00 * x / sqrtpi
      else if ( y .le. 1.0D+00 ) then
        value = x * ( 1.0D+00 
     &    + r8_csevl ( 2.0D+00 * x * x - 1.0D+00, erfcs, nterf ) )
      else if ( y .le. xbig ) then
        value = 1.0D+00 - r8_erfc ( y )
        if ( x .lt. 0.0D+00 ) then
          value = - value
        end if
      else
        value = 1.0D+00
        if ( x .lt. 0.0D+00 ) then
          value = - value
        end if
      end if

      r8_erf = value

      return
      end
      function r8_erfc ( x )

c*********************************************************************72
c
cc R8_ERFC evaluates the co-error function of an R8 argument.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    14 March 2010
c
c  Author:
c
c    Original FORTRAN77 version by Wayne Fullerton.
c    This FORTRAN77 version by John Burkardt.
c
c  Reference:
c
c    Wayne Fullerton,
c    Portable Special Function Routines,
c    in Portability of Numerical Software,
c    edited by Wayne Cowell,
c    Lecture Notes in Computer Science, Volume 57,
c    Springer 1977,
c    ISBN: 978-3-540-08446-4,
c    LC: QA297.W65.
c
c  Parameters:
c
c    Input, double precision X, the argument.
c
c    Output, double precision R8_ERFC, the co-error function of X.
c
      implicit none

      double precision erc2cs(49)
      double precision erfccs(59)
      double precision erfcs(21)
      double precision eta
      integer nterc2
      integer nterf
      integer nterfc
      double precision r8_csevl
      double precision r8_erfc
      integer r8_inits
      double precision r8_mach
      double precision sqeps
      double precision sqrtpi
      double precision x
      double precision xmax
      double precision xsml
      double precision y

      save erfccs
      save erfcs
      save erc2cs
      save nterc2
      save nterf
      save nterfc
      save sqeps
      save sqrtpi
      save xmax
      save xsml

      data erfcs(  1) / -0.49046121234691808039984544033376D-01 /
      data erfcs(  2) / -0.14226120510371364237824741899631D+00 /
      data erfcs(  3) / +0.10035582187599795575754676712933D-01 /
      data erfcs(  4) / -0.57687646997674847650827025509167D-03 /
      data erfcs(  5) / +0.27419931252196061034422160791471D-04 /
      data erfcs(  6) / -0.11043175507344507604135381295905D-05 /
      data erfcs(  7) / +0.38488755420345036949961311498174D-07 /
      data erfcs(  8) / -0.11808582533875466969631751801581D-08 /
      data erfcs(  9) / +0.32334215826050909646402930953354D-10 /
      data erfcs( 10) / -0.79910159470045487581607374708595D-12 /
      data erfcs( 11) / +0.17990725113961455611967245486634D-13 /
      data erfcs( 12) / -0.37186354878186926382316828209493D-15 /
      data erfcs( 13) / +0.71035990037142529711689908394666D-17 /
      data erfcs( 14) / -0.12612455119155225832495424853333D-18 /
      data erfcs( 15) / +0.20916406941769294369170500266666D-20 /
      data erfcs( 16) / -0.32539731029314072982364160000000D-22 /
      data erfcs( 17) / +0.47668672097976748332373333333333D-24 /
      data erfcs( 18) / -0.65980120782851343155199999999999D-26 /
      data erfcs( 19) / +0.86550114699637626197333333333333D-28 /
      data erfcs( 20) / -0.10788925177498064213333333333333D-29 /
      data erfcs( 21) / +0.12811883993017002666666666666666D-31 /

      data erc2cs(  1) / -0.6960134660230950112739150826197D-01 /
      data erc2cs(  2) / -0.4110133936262089348982212084666D-01 /
      data erc2cs(  3) / +0.3914495866689626881561143705244D-02 /
      data erc2cs(  4) / -0.4906395650548979161280935450774D-03 /
      data erc2cs(  5) / +0.7157479001377036380760894141825D-04 /
      data erc2cs(  6) / -0.1153071634131232833808232847912D-04 /
      data erc2cs(  7) / +0.1994670590201997635052314867709D-05 /
      data erc2cs(  8) / -0.3642666471599222873936118430711D-06 /
      data erc2cs(  9) / +0.6944372610005012589931277214633D-07 /
      data erc2cs( 10) / -0.1371220902104366019534605141210D-07 /
      data erc2cs( 11) / +0.2788389661007137131963860348087D-08 /
      data erc2cs( 12) / -0.5814164724331161551864791050316D-09 /
      data erc2cs( 13) / +0.1238920491752753181180168817950D-09 /
      data erc2cs( 14) / -0.2690639145306743432390424937889D-10 /
      data erc2cs( 15) / +0.5942614350847910982444709683840D-11 /
      data erc2cs( 16) / -0.1332386735758119579287754420570D-11 /
      data erc2cs( 17) / +0.3028046806177132017173697243304D-12 /
      data erc2cs( 18) / -0.6966648814941032588795867588954D-13 /
      data erc2cs( 19) / +0.1620854541053922969812893227628D-13 /
      data erc2cs( 20) / -0.3809934465250491999876913057729D-14 /
      data erc2cs( 21) / +0.9040487815978831149368971012975D-15 /
      data erc2cs( 22) / -0.2164006195089607347809812047003D-15 /
      data erc2cs( 23) / +0.5222102233995854984607980244172D-16 /
      data erc2cs( 24) / -0.1269729602364555336372415527780D-16 /
      data erc2cs( 25) / +0.3109145504276197583836227412951D-17 /
      data erc2cs( 26) / -0.7663762920320385524009566714811D-18 /
      data erc2cs( 27) / +0.1900819251362745202536929733290D-18 /
      data erc2cs( 28) / -0.4742207279069039545225655999965D-19 /
      data erc2cs( 29) / +0.1189649200076528382880683078451D-19 /
      data erc2cs( 30) / -0.3000035590325780256845271313066D-20 /
      data erc2cs( 31) / +0.7602993453043246173019385277098D-21 /
      data erc2cs( 32) / -0.1935909447606872881569811049130D-21 /
      data erc2cs( 33) / +0.4951399124773337881000042386773D-22 /
      data erc2cs( 34) / -0.1271807481336371879608621989888D-22 /
      data erc2cs( 35) / +0.3280049600469513043315841652053D-23 /
      data erc2cs( 36) / -0.8492320176822896568924792422399D-24 /
      data erc2cs( 37) / +0.2206917892807560223519879987199D-24 /
      data erc2cs( 38) / -0.5755617245696528498312819507199D-25 /
      data erc2cs( 39) / +0.1506191533639234250354144051199D-25 /
      data erc2cs( 40) / -0.3954502959018796953104285695999D-26 /
      data erc2cs( 41) / +0.1041529704151500979984645051733D-26 /
      data erc2cs( 42) / -0.2751487795278765079450178901333D-27 /
      data erc2cs( 43) / +0.7290058205497557408997703680000D-28 /
      data erc2cs( 44) / -0.1936939645915947804077501098666D-28 /
      data erc2cs( 45) / +0.5160357112051487298370054826666D-29 /
      data erc2cs( 46) / -0.1378419322193094099389644800000D-29 /
      data erc2cs( 47) / +0.3691326793107069042251093333333D-30 /
      data erc2cs( 48) / -0.9909389590624365420653226666666D-31 /
      data erc2cs( 49) / +0.2666491705195388413323946666666D-31 /

      data erfccs(  1) / +0.715179310202924774503697709496D-01 /
      data erfccs(  2) / -0.265324343376067157558893386681D-01 /
      data erfccs(  3) / +0.171115397792085588332699194606D-02 /
      data erfccs(  4) / -0.163751663458517884163746404749D-03 /
      data erfccs(  5) / +0.198712935005520364995974806758D-04 /
      data erfccs(  6) / -0.284371241276655508750175183152D-05 /
      data erfccs(  7) / +0.460616130896313036969379968464D-06 /
      data erfccs(  8) / -0.822775302587920842057766536366D-07 /
      data erfccs(  9) / +0.159214187277090112989358340826D-07 /
      data erfccs( 10) / -0.329507136225284321486631665072D-08 /
      data erfccs( 11) / +0.722343976040055546581261153890D-09 /
      data erfccs( 12) / -0.166485581339872959344695966886D-09 /
      data erfccs( 13) / +0.401039258823766482077671768814D-10 /
      data erfccs( 14) / -0.100481621442573113272170176283D-10 /
      data erfccs( 15) / +0.260827591330033380859341009439D-11 /
      data erfccs( 16) / -0.699111056040402486557697812476D-12 /
      data erfccs( 17) / +0.192949233326170708624205749803D-12 /
      data erfccs( 18) / -0.547013118875433106490125085271D-13 /
      data erfccs( 19) / +0.158966330976269744839084032762D-13 /
      data erfccs( 20) / -0.472689398019755483920369584290D-14 /
      data erfccs( 21) / +0.143587337678498478672873997840D-14 /
      data erfccs( 22) / -0.444951056181735839417250062829D-15 /
      data erfccs( 23) / +0.140481088476823343737305537466D-15 /
      data erfccs( 24) / -0.451381838776421089625963281623D-16 /
      data erfccs( 25) / +0.147452154104513307787018713262D-16 /
      data erfccs( 26) / -0.489262140694577615436841552532D-17 /
      data erfccs( 27) / +0.164761214141064673895301522827D-17 /
      data erfccs( 28) / -0.562681717632940809299928521323D-18 /
      data erfccs( 29) / +0.194744338223207851429197867821D-18 /
      data erfccs( 30) / -0.682630564294842072956664144723D-19 /
      data erfccs( 31) / +0.242198888729864924018301125438D-19 /
      data erfccs( 32) / -0.869341413350307042563800861857D-20 /
      data erfccs( 33) / +0.315518034622808557122363401262D-20 /
      data erfccs( 34) / -0.115737232404960874261239486742D-20 /
      data erfccs( 35) / +0.428894716160565394623737097442D-21 /
      data erfccs( 36) / -0.160503074205761685005737770964D-21 /
      data erfccs( 37) / +0.606329875745380264495069923027D-22 /
      data erfccs( 38) / -0.231140425169795849098840801367D-22 /
      data erfccs( 39) / +0.888877854066188552554702955697D-23 /
      data erfccs( 40) / -0.344726057665137652230718495566D-23 /
      data erfccs( 41) / +0.134786546020696506827582774181D-23 /
      data erfccs( 42) / -0.531179407112502173645873201807D-24 /
      data erfccs( 43) / +0.210934105861978316828954734537D-24 /
      data erfccs( 44) / -0.843836558792378911598133256738D-25 /
      data erfccs( 45) / +0.339998252494520890627359576337D-25 /
      data erfccs( 46) / -0.137945238807324209002238377110D-25 /
      data erfccs( 47) / +0.563449031183325261513392634811D-26 /
      data erfccs( 48) / -0.231649043447706544823427752700D-26 /
      data erfccs( 49) / +0.958446284460181015263158381226D-27 /
      data erfccs( 50) / -0.399072288033010972624224850193D-27 /
      data erfccs( 51) / +0.167212922594447736017228709669D-27 /
      data erfccs( 52) / -0.704599152276601385638803782587D-28 /
      data erfccs( 53) / +0.297976840286420635412357989444D-28 /
      data erfccs( 54) / -0.126252246646061929722422632994D-28 /
      data erfccs( 55) / +0.539543870454248793985299653154D-29 /
      data erfccs( 56) / -0.238099288253145918675346190062D-29 /
      data erfccs( 57) / +0.109905283010276157359726683750D-29 /
      data erfccs( 58) / -0.486771374164496572732518677435D-30 /
      data erfccs( 59) / +0.152587726411035756763200828211D-30 /

      data nterc2 / 0 /
      data nterf / 0 /
      data nterfc / 0 /
      data sqeps / 0.0D+00 /
      data sqrtpi / 1.77245385090551602729816748334115D+00 /
      data xmax / 0.0D+00 /
      data xsml / 0.0D+00 /

      if ( nterf .eq. 0 ) then

        eta = 0.1D+00 * r8_mach ( 3 )
        nterf = r8_inits ( erfcs, 21, eta )
        nterfc = r8_inits ( erfccs, 59, eta )
        nterc2 = r8_inits ( erc2cs, 49, eta )

        xsml = - dsqrt ( - dlog ( sqrtpi * r8_mach ( 3 ) ) )
        xmax = dsqrt (- dlog ( sqrtpi * r8_mach ( 1 ) ) )
        xmax = xmax - 0.5D+00 * dlog ( xmax ) / xmax - 0.01D+00
        sqeps = dsqrt ( 2.0D+00 * r8_mach ( 3 ) )

      end if

      if ( x .le. xsml ) then

        r8_erfc = 2.0D+00
        return

      end if

      if ( xmax .lt. x ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'R8_ERFC - Warning!'
        write ( *, '(a)' ) '  X so big that ERFC underflows.'
        r8_erfc = 0.0D+00
        return
      end if

      y = dabs ( x )

      if ( y .lt. sqeps ) then
        r8_erfc = 1.0D+00 - 2.0D+00 * x / sqrtpi
        return
      else if ( y .le. 1.0D+00 ) then
        r8_erfc = 1.0D+00 - x * ( 1.0D+00 
     &    + r8_csevl ( 2.0D+00 * x * x - 1.0D+00, erfcs, nterf ) )
        return
      end if

      y = y * y

      if ( y .le. 4.0D+00 ) then
        r8_erfc = dexp ( - y ) / dabs ( x ) * ( 0.5D+00 
     &    + r8_csevl ( ( 8.0D+00 / y - 5.0D+00 ) / 3.0D+00, erc2cs, 
     &    nterc2 ) )
      else 
        r8_erfc = dexp ( - y ) / dabs ( x ) * ( 0.5D+00 
     &    + r8_csevl ( 8.0D+00 / y - 1.0D+00, erfccs, nterfc ) )
      end if

      if ( x .lt. 0.0D+00 ) then
        r8_erfc = 2.0D+00 - r8_erfc
      end if

      return
      end
      subroutine r8_gaml ( xmin, xmax )

c*********************************************************************72
c
cc R8_GAML evaluates bounds for an R8 argument of the gamma function.
c
c  Discussion:
c
c    This function calculates the minimum and maximum legal bounds 
c    for X in the evaluation of GAMMA ( X ).
c
c    XMIN and XMAX are not the only bounds, but they are the only 
c    non-trivial ones to calculate.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    20 March 2010
c
c  Author:
c
c    Original FORTRAN77 version by Wayne Fullerton.
c    This FORTRAN77 version by John Burkardt.
c
c  Reference:
c
c    Wayne Fullerton,
c    Portable Special Function Routines,
c    in Portability of Numerical Software,
c    edited by Wayne Cowell,
c    Lecture Notes in Computer Science, Volume 57,
c    Springer 1977,
c    ISBN: 978-3-540-08446-4,
c    LC: QA297.W65.
c
c  Parameters:
c
c    Output, double precision XMIN, XMAX, the bounds.
c
      implicit none

      double precision alnbig
      double precision alnsml
      integer i
      integer j
      double precision r8_mach
      double precision xln
      double precision xmax
      double precision xmin
      double precision xold

      alnsml = dlog ( r8_mach ( 1 ) )
      xmin = - alnsml

      do i = 1, 10

        xold = xmin
        xln = dlog ( xmin )
        xmin = xmin - xmin * ( ( xmin + 0.5D+00 ) * xln - xmin 
     &    - 0.2258D+00 + alnsml ) / ( xmin * xln + 0.5D+00 )

        if ( dabs ( xmin - xold ) .lt. 0.005D+00 ) then

          xmin = - xmin + 0.01D+00

          alnbig = dlog ( r8_mach ( 2 ) )
          xmax = alnbig

          do j = 1, 10

            xold = xmax
            xln = dlog ( xmax )
            xmax = xmax - xmax * ( ( xmax - 0.5D+00 ) * xln - xmax 
     &        + 0.9189D+00 - alnbig ) / ( xmax * xln - 0.5D+00 )

            if ( dabs ( xmax - xold ) .lt. 0.005D+00 ) then
              xmax = xmax - 0.01D+00
              xmin = dmax1 ( xmin, - xmax + 1.0D+00 )
              return
            end if

          end do

          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) 'R8_GAML - Fatal error!'
          write ( *, '(a)' ) '  Unable to find XMAX.'
          stop

        end if

      end do

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'R8_GAML - Fatal error!'
      write ( *, '(a)' ) '  Unable to find XMIN.'

      stop
      end
      function r8_gamma ( x )

c*********************************************************************72
c
cc R8_GAMMA evaluates the gamma function of an R8 argument.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    15 March 2010
c
c  Author:
c
c    Original FORTRAN77 version by Wayne Fullerton.
c    This FORTRAN77 version by John Burkardt.
c
c  Reference:
c
c    Wayne Fullerton,
c    Portable Special Function Routines,
c    in Portability of Numerical Software,
c    edited by Wayne Cowell,
c    Lecture Notes in Computer Science, Volume 57,
c    Springer 1977,
c    ISBN: 978-3-540-08446-4,
c    LC: QA297.W65.
c
c  Parameters:
c
c    Input, double precision X, the argument.
c
c    Output, double precision R8_GAMMA, the gamma function of X.
c
      implicit none

      double precision dxrel
      double precision gcs(42)
      integer i
      integer n
      integer ngcs
      double precision pi
      double precision r8_csevl
      double precision r8_gamma
      integer r8_inits
      double precision r8_lgmc
      double precision r8_mach
      double precision sinpiy
      double precision sq2pil
      double precision x
      double precision xmax
      double precision xmin
      double precision xsml
      double precision y

      save dxrel
      save gcs
      save ngcs
      save pi
      save sq2pil
      save xmax
      save xmin
      save xsml

      data gcs(  1) / +0.8571195590989331421920062399942D-02 /
      data gcs(  2) / +0.4415381324841006757191315771652D-02 /
      data gcs(  3) / +0.5685043681599363378632664588789D-01 /
      data gcs(  4) / -0.4219835396418560501012500186624D-02 /
      data gcs(  5) / +0.1326808181212460220584006796352D-02 /
      data gcs(  6) / -0.1893024529798880432523947023886D-03 /
      data gcs(  7) / +0.3606925327441245256578082217225D-04 /
      data gcs(  8) / -0.6056761904460864218485548290365D-05 /
      data gcs(  9) / +0.1055829546302283344731823509093D-05 /
      data gcs( 10) / -0.1811967365542384048291855891166D-06 /
      data gcs( 11) / +0.3117724964715322277790254593169D-07 /
      data gcs( 12) / -0.5354219639019687140874081024347D-08 /
      data gcs( 13) / +0.9193275519859588946887786825940D-09 /
      data gcs( 14) / -0.1577941280288339761767423273953D-09 /
      data gcs( 15) / +0.2707980622934954543266540433089D-10 /
      data gcs( 16) / -0.4646818653825730144081661058933D-11 /
      data gcs( 17) / +0.7973350192007419656460767175359D-12 /
      data gcs( 18) / -0.1368078209830916025799499172309D-12 /
      data gcs( 19) / +0.2347319486563800657233471771688D-13 /
      data gcs( 20) / -0.4027432614949066932766570534699D-14 /
      data gcs( 21) / +0.6910051747372100912138336975257D-15 /
      data gcs( 22) / -0.1185584500221992907052387126192D-15 /
      data gcs( 23) / +0.2034148542496373955201026051932D-16 /
      data gcs( 24) / -0.3490054341717405849274012949108D-17 /
      data gcs( 25) / +0.5987993856485305567135051066026D-18 /
      data gcs( 26) / -0.1027378057872228074490069778431D-18 /
      data gcs( 27) / +0.1762702816060529824942759660748D-19 /
      data gcs( 28) / -0.3024320653735306260958772112042D-20 /
      data gcs( 29) / +0.5188914660218397839717833550506D-21 /
      data gcs( 30) / -0.8902770842456576692449251601066D-22 /
      data gcs( 31) / +0.1527474068493342602274596891306D-22 /
      data gcs( 32) / -0.2620731256187362900257328332799D-23 /
      data gcs( 33) / +0.4496464047830538670331046570666D-24 /
      data gcs( 34) / -0.7714712731336877911703901525333D-25 /
      data gcs( 35) / +0.1323635453126044036486572714666D-25 /
      data gcs( 36) / -0.2270999412942928816702313813333D-26 /
      data gcs( 37) / +0.3896418998003991449320816639999D-27 /
      data gcs( 38) / -0.6685198115125953327792127999999D-28 /
      data gcs( 39) / +0.1146998663140024384347613866666D-28 /
      data gcs( 40) / -0.1967938586345134677295103999999D-29 /
      data gcs( 41) / +0.3376448816585338090334890666666D-30 /
      data gcs( 42) / -0.5793070335782135784625493333333D-31 /

      data dxrel / 0.0D+00 /
      data ngcs / 0 /
      data pi / 3.14159265358979323846264338327950D+00 /
      data sq2pil / 0.91893853320467274178032973640562D+00 /
      data xmax / 0.0D+00 /
      data xmin / 0.0D+00 /
      data xsml / 0.0D+00 /

      if ( ngcs .eq. 0 ) then
        ngcs = r8_inits ( gcs, 42, 0.1D+00 * r8_mach ( 3 ) )
        call r8_gaml ( xmin, xmax )
        xsml = dexp ( dmax1 ( dlog ( r8_mach ( 1 ) ),
     &    - dlog ( r8_mach ( 2 ) ) ) + 0.01D+00 )
        dxrel = dsqrt ( r8_mach ( 4 ) )
      end if

      y = dabs ( x )

      if ( y .le. 10.0D+00 ) then

        n = int ( x )
        if ( x .lt. 0.0D+00 ) then
          n = n - 1
        end if
        y = x - dble ( n )
        n = n - 1
        r8_gamma = 0.9375D+00 + r8_csevl ( 2.0D+00 * y - 1.0D+00, 
     &    gcs, ngcs )

        if ( n .eq. 0 ) then

          return

        else if ( n .lt. 0 ) then

          n = - n

          if ( x .eq. 0.0D+00 ) then
            write ( *, '(a)' ) ' '
            write ( *, '(a)' ) 'R8_GAMMA - Fatal error!'
            write ( *, '(a)' ) '  X is 0.'
            stop
          end if

          if ( x .lt. 0.0D+00 .and. 
     &      x + dble ( n - 2 ) .eq. 0.0D+00 ) then
            write ( *, '(a)' ) ' '
            write ( *, '(a)' ) 'R8_GAMMA - Fatal error!'
            write ( *, '(a)' ) '  X is a negative integer.'
            stop
          end if

          if ( x .lt. - 0.5D+00 .and. 
     &      dabs ( ( x - dint ( x - 0.5D+00 ) ) / x ) .lt. dxrel ) then
            write ( *, '(a)' ) ' '
            write ( *, '(a)' ) 'R8_GAMMA - Warning!'
            write ( *, '(a)' ) '  X too near a negative integer,'
            write ( *, '(a)' ) '  answer is half precision.'
          end if

          if ( y .lt. xsml ) then
            write ( *, '(a)' ) ' '
            write ( *, '(a)' ) 'R8_GAMMA - Fatal error!'
            write ( *, '(a)' ) 
     &        '  X is so close to zero that Gamma overflows.'
            stop
          end if

          do i = 1, n
            r8_gamma = r8_gamma / ( x + dble ( i - 1 ) )
          end do

        else if ( n .eq. 0 ) then

        else

          do i = 1, n
            r8_gamma = ( y + dble ( i ) ) * r8_gamma
          end do

        end if

      else

        if ( xmax .lt. x ) then
          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) 'R8_GAMMA - Fatal error!'
          write ( *, '(a)' ) '  X so big that Gamma overflows.'
          stop
        end if
c
c  Underflow.
c
        if ( x .lt. xmin ) then
          r8_gamma = 0.0D+00
          return
        end if

        r8_gamma = dexp ( ( y - 0.5D+00 ) * dlog ( y ) - y + sq2pil 
     &    + r8_lgmc ( y ) )

        if ( 0.0D+00 .lt. x ) then
          return
        end if

        if ( dabs ( ( x - dint ( x - 0.5D+00 ) ) / x ) .lt. dxrel ) then
          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) 'R8_GAMMA - Warning!'
          write ( *, '(a)' ) '  X too near a negative integer,'
          write ( *, '(a)' ) '  answer is half precision.'
        end if

        sinpiy = dsin ( pi * y )

        if ( sinpiy .eq. 0.0D+00 ) then
          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) 'R8_GAMMA - Fatal error!'
          write ( *, '(a)' ) '  X is a negative integer.'
          stop
        end if

        r8_gamma = - pi / ( y * sinpiy * r8_gamma )

      end if

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
      function r8_inits ( dos, nos, eta )

c*********************************************************************72
c
cc R8_INITS initializes a Chebyshev series.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    09 March 2010
c
c  Author:
c
c    This FORTRAN77 version by John Burkardt.
c
c  Reference:
c
c    Roger Broucke,
c    Algorithm 446:
c    Ten Subroutines for the Manipulation of Chebyshev Series,
c    Communications of the ACM,
c    Volume 16, Number 4, April 1973, pages 254-256.
c
c  Parameters:
c
c    Input, double precision DOS(NOS), the Chebyshev coefficients.
c
c    Input, integer NOS, the number of coefficients.
c
c    Input, double precision ETA, the desired accuracy.
c
c    Output, integer R8_INITS, the number of terms of the series needed
c    to ensure the requested accuracy.
c
      implicit none

      integer nos

      double precision dos(nos)
      double precision err 
      double precision eta
      integer i
      integer r8_inits

      if ( nos .lt. 1 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'R8_INITS - Fatal error!'
        write ( *, '(a)' ) '  Number of coefficients < 1.'
        stop
      end if

      err = 0.0D+00

      do i = nos, 1, -1
        err = err + dabs ( dos(i) )
        if ( eta .lt. err ) then
          r8_inits = i
          return
        end if
      end do

      r8_inits = nos
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'R8_INITS - Warning!'
      write ( *, '(a)' ) '  ETA may be too small.'

      return
      end
      function r8_lgmc ( x )

c*********************************************************************72
c
cc R8_LGMC evaluates the log gamma correction factor for an R8 argument.
c
c  Discussion:
c
c    For 10 <= X, compute the log gamma correction factor so that
c
c      log ( gamma ( x ) ) = log ( sqrt ( 2 * pi ) ) 
c                          + ( x - 0.5 ) * log ( x ) - x 
c                          + r8_lgmc ( x )
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    16 March 2010
c
c  Author:
c
c    Original FORTRAN77 version by Wayne Fullerton.
c    This FORTRAN77 version by John Burkardt.
c
c  Reference:
c
c    Wayne Fullerton,
c    Portable Special Function Routines,
c    in Portability of Numerical Software,
c    edited by Wayne Cowell,
c    Lecture Notes in Computer Science, Volume 57,
c    Springer 1977,
c    ISBN: 978-3-540-08446-4,
c    LC: QA297.W65.
c
c  Parameters:
c
c    Input, double precision X, the argument.
c
c    Output, double precision R8_LGMC, the correction factor.
c
      implicit none

      double precision algmcs(15)
      integer nalgm
      double precision r8_csevl
      integer r8_inits
      double precision r8_lgmc
      double precision r8_mach
      double precision x
      double precision xbig
      double precision xmax

      save algmcs
      save nalgm
      save xbig
      save xmax

      data algmcs(  1) / +0.1666389480451863247205729650822D+00 /
      data algmcs(  2) / -0.1384948176067563840732986059135D-04 /
      data algmcs(  3) / +0.9810825646924729426157171547487D-08 /
      data algmcs(  4) / -0.1809129475572494194263306266719D-10 /
      data algmcs(  5) / +0.6221098041892605227126015543416D-13 /
      data algmcs(  6) / -0.3399615005417721944303330599666D-15 /
      data algmcs(  7) / +0.2683181998482698748957538846666D-17 /
      data algmcs(  8) / -0.2868042435334643284144622399999D-19 /
      data algmcs(  9) / +0.3962837061046434803679306666666D-21 /
      data algmcs( 10) / -0.6831888753985766870111999999999D-23 /
      data algmcs( 11) / +0.1429227355942498147573333333333D-24 /
      data algmcs( 12) / -0.3547598158101070547199999999999D-26 /
      data algmcs( 13) / +0.1025680058010470912000000000000D-27 /
      data algmcs( 14) / -0.3401102254316748799999999999999D-29 /
      data algmcs( 15) / +0.1276642195630062933333333333333D-30 /

      data nalgm / 0 /
      data xbig / 0.0D+00 /
      data xmax / 0.0D+00 /

      if ( nalgm .eq. 0 ) then
        nalgm = r8_inits ( algmcs, 15, r8_mach ( 3 ) )
        xbig = 1.0D+00 / dsqrt ( r8_mach ( 3 ) )
        xmax = dexp ( dmin1 ( dlog ( r8_mach ( 2 ) / 12.0D+00 ), 
     &    - dlog ( 12.0D+00 * r8_mach ( 1 ) ) ) )
      end if

      if ( x .lt. 10.0D+00 ) then

        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'R8_LGMC - Fatal error!'
        write ( *, '(a)' ) '  X must be at least 10.'
        stop

      else if ( x .lt. xbig ) then

        r8_lgmc = r8_csevl ( 2.0D+00 * ( 10.0D+00 / x ) 
     &    * ( 10.0D+00 / x ) - 1.0D+00, algmcs, nalgm ) / x

      else if ( x .lt. xmax ) then

        r8_lgmc = 1.0D+00 / ( 12.0D+00 * x )

      else

        r8_lgmc = 0.0D+00

      end if

      return
      end
      function r8_mach ( i )

c*********************************************************************72
c
cc R8_MACH returns double precision real machine-dependent constants.
c
c  Discussion:
c
c    R8_MACH can be used to obtain machine-dependent parameters
c    for the local machine environment.  It is a function
c    with one input argument, and can be called as follows:
c
c      D = R8_MACH ( I )
c
c    where I=1,...,5.  The output value of D above is
c    determined by the input value of I:.
c
c    R8_MACH ( 1) = B^(EMIN-1), the smallest positive magnitude.
c    R8_MACH ( 2) = B^EMAX*(1 - B**(-T)), the largest magnitude.
c    R8_MACH ( 3) = B^(-T), the smallest relative spacing.
c    R8_MACH ( 4) = B^(1-T), the largest relative spacing.
c    R8_MACH ( 5) = LOG10(B)
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
c    Output, double precision R8_MACH, the value of the constant.
c
      implicit none

      double precision r8_mach
      integer i

      if ( i .lt. 1 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'R8_MACH - Fatal error!'
        write ( *, '(a)' ) '  The input argument I is out of bounds.'
        write ( *, '(a)' ) '  Legal values satisfy 1 .le. I .le. 5.'
        write ( *, '(a,i12)' ) '  I = ', i
        r8_mach = 0.0D+00
        stop
      else if ( i .eq. 1 ) then
        r8_mach = 4.450147717014403D-308
      else if ( i .eq. 2 ) then
        r8_mach = 8.988465674311579D+307
      else if ( i .eq. 3 ) then
        r8_mach = 1.110223024625157D-016
      else if ( i .eq. 4 ) then
        r8_mach = 2.220446049250313D-016
      else if ( i .eq. 5 ) then
        r8_mach = 0.301029995663981D+000
      else if ( 5 .lt. i ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'R8_MACH - Fatal error!'
        write ( *, '(a)' ) '  The input argument I is out of bounds.'
        write ( *, '(a)' ) '  Legal values satisfy 1 .le. I .le. 5.'
        write ( *, '(a,i12)' ) '  I = ', i
        r8_mach = 0.0D+00
        stop
      end if

      return
      end
      function r8_sech ( x )

c*********************************************************************72
c
cc R8_SECH evaluates the hyperbolic secant, while avoiding COSH overflow.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    29 August 2000
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, double precision X, the argument of the function.
c
c    Output, double precision R8_SECH, the value of the function.
c
      implicit none

      double precision log_huge 
      parameter ( log_huge = 80.0D+00 )
      double precision r8_sech
      double precision x

      if ( log_huge .lt. abs ( x ) ) then
        r8_sech = 0.0D+00
      else
        r8_sech = 1.0D+00 / cosh ( x )
      end if

      return
      end
      function r8_si ( x )

c*********************************************************************72
c
cc R8_SI evaluates the sine integral Si of an R8 argument.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    26 March 2010
c
c  Author:
c
c    Original FORTRAN77 version by Wayne Fullerton.
c    This FORTRAN77 version by John Burkardt.
c
c  Reference:
c
c    Wayne Fullerton,
c    Portable Special Function Routines,
c    in Portability of Numerical Software,
c    edited by Wayne Cowell,
c    Lecture Notes in Computer Science, Volume 57,
c    Springer 1977,
c    ISBN: 978-3-540-08446-4,
c    LC: QA297.W65.
c
c  Parameters:
c
c    Input, double precision X, the argument.
c
c    Output, double precision R8_SI, the sine integral Si evaluated at X.
c
      implicit none

      double precision absx
      double precision cosx
      double precision f
      double precision g
      integer nsi
      double precision pi2
      double precision r8_csevl
      integer r8_inits
      double precision r8_mach
      double precision r8_si
      double precision sics(18)
      double precision x
      double precision xsml

      save nsi
      save pi2
      save sics
      save xsml

      data sics(  1) / -0.1315646598184841928904275173000457D+00 /
      data sics(  2) / -0.2776578526973601892048287660157299D+00 /
      data sics(  3) /  0.0354414054866659179749135464710086D+00 /
      data sics(  4) / -0.0025631631447933977658752788361530D+00 /
      data sics(  5) /  0.0001162365390497009281264921482985D+00 /
      data sics(  6) / -0.0000035904327241606042670004347148D+00 /
      data sics(  7) /  0.0000000802342123705710162308652976D+00 /
      data sics(  8) / -0.0000000013562997692540250649931846D+00 /
      data sics(  9) /  0.0000000000179440721599736775567759D+00 /
      data sics( 10) / -0.0000000000001908387343087145490737D+00 /
      data sics( 11) /  0.0000000000000016669989586824330853D+00 /
      data sics( 12) / -0.0000000000000000121730988368503042D+00 /
      data sics( 13) /  0.0000000000000000000754181866993865D+00 /
      data sics( 14) / -0.0000000000000000000004014178842446D+00 /
      data sics( 15) /  0.0000000000000000000000018553690716D+00 /
      data sics( 16) / -0.0000000000000000000000000075166966D+00 /
      data sics( 17) /  0.0000000000000000000000000000269113D+00 /
      data sics( 18) / -0.0000000000000000000000000000000858D+00 /

      data nsi / 0 /
      data pi2 / 1.57079632679489661923132169163975D+00 /
      data xsml / 0.0D+00 /

      if ( nsi .eq. 0 ) then
        nsi = r8_inits ( sics, 18, 0.1D+00 * r8_mach ( 3 ) )
        xsml = dsqrt ( r8_mach ( 3 ) )
      end if

      absx = dabs ( x )

      if ( absx .lt. xsml ) then
        r8_si = x
      else if ( absx .le. 4.0D+00 ) then
        r8_si = x * ( 0.75D+00 
     &    + r8_csevl ( ( x * x - 8.0D+00 ) * 0.125D+00, sics, nsi ) )
      else
        call r8_sifg ( absx, f, g )
        cosx = dcos ( absx )
        r8_si = pi2 - f * cosx - g * dsin ( x )
        if ( x .lt. 0.0D+00 ) then
          r8_si = - r8_si
        end if
      end if

      return
      end
      subroutine r8_sifg ( x, f, g )

c*********************************************************************72
c
cc R8_SIFG is a utility routine.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    25 March 2010
c
c  Author:
c
c    Original FORTRAN77 version by Wayne Fullerton.
c    This FORTRAN77 version by John Burkardt.
c
c  Reference:
c
c    Wayne Fullerton,
c    Portable Special Function Routines,
c    in Portability of Numerical Software,
c    edited by Wayne Cowell,
c    Lecture Notes in Computer Science, Volume 57,
c    Springer 1977,
c    ISBN: 978-3-540-08446-4,
c    LC: QA297.W65.
c
c  Parameters:
c
c    Input, double precision X, the argument.
c
c    Output, double precision F, G.
c
      implicit none

      double precision f
      double precision f1cs(43)
      double precision f2cs(99)
      double precision g
      double precision g1cs(44)
      double precision g2cs(44)
      double precision g3cs(56)
      integer nf1
      integer nf2
      integer ng1
      integer ng2
      integer ng3
      double precision r8_csevl
      integer r8_inits
      double precision r8_mach
      double precision tol
      double precision x
      double precision xbig
      double precision xbnd
      double precision xbndg
      double precision xmaxf
      double precision xmaxg

      save f1cs
      save f2cs
      save g1cs
      save g2cs
      save g3cs
      save nf1
      save nf2
      save ng1
      save ng2
      save ng3
      save xbig
      save xbnd
      save xbndg
      save xmaxf
      save xmaxg

      data f1cs(  1) / -0.1191081969051363610348201965828918D+00 /
      data f1cs(  2) / -0.0247823144996236247590074150823133D+00 /
      data f1cs(  3) /  0.0011910281453357821268120363054457D+00 /
      data f1cs(  4) / -0.0000927027714388561748308600360706D+00 /
      data f1cs(  5) /  0.0000093373141568270996868204582766D+00 /
      data f1cs(  6) / -0.0000011058287820557143938979426306D+00 /
      data f1cs(  7) /  0.0000001464772071460162169336550799D+00 /
      data f1cs(  8) / -0.0000000210694496287689532601227548D+00 /
      data f1cs(  9) /  0.0000000032293492366848236382857374D+00 /
      data f1cs( 10) / -0.0000000005206529617529375828014986D+00 /
      data f1cs( 11) /  0.0000000000874878884570278750268316D+00 /
      data f1cs( 12) / -0.0000000000152176187056123668294574D+00 /
      data f1cs( 13) /  0.0000000000027257192405419573900583D+00 /
      data f1cs( 14) / -0.0000000000005007053075968556290255D+00 /
      data f1cs( 15) /  0.0000000000000940240902726068511779D+00 /
      data f1cs( 16) / -0.0000000000000180014444791803678336D+00 /
      data f1cs( 17) /  0.0000000000000035062621432741785826D+00 /
      data f1cs( 18) / -0.0000000000000006935282926769149709D+00 /
      data f1cs( 19) /  0.0000000000000001390925136454216568D+00 /
      data f1cs( 20) / -0.0000000000000000282486885074170585D+00 /
      data f1cs( 21) /  0.0000000000000000058031305693579081D+00 /
      data f1cs( 22) / -0.0000000000000000012046901573375820D+00 /
      data f1cs( 23) /  0.0000000000000000002525052443655940D+00 /
      data f1cs( 24) / -0.0000000000000000000533980268805594D+00 /
      data f1cs( 25) /  0.0000000000000000000113855786274122D+00 /
      data f1cs( 26) / -0.0000000000000000000024462861505259D+00 /
      data f1cs( 27) /  0.0000000000000000000005293659320439D+00 /
      data f1cs( 28) / -0.0000000000000000000001153184940277D+00 /
      data f1cs( 29) /  0.0000000000000000000000252786568318D+00 /
      data f1cs( 30) / -0.0000000000000000000000055738645378D+00 /
      data f1cs( 31) /  0.0000000000000000000000012358245621D+00 /
      data f1cs( 32) / -0.0000000000000000000000002754350842D+00 /
      data f1cs( 33) /  0.0000000000000000000000000616906808D+00 /
      data f1cs( 34) / -0.0000000000000000000000000138817443D+00 /
      data f1cs( 35) /  0.0000000000000000000000000031375329D+00 /
      data f1cs( 36) / -0.0000000000000000000000000007121249D+00 /
      data f1cs( 37) /  0.0000000000000000000000000001622778D+00 /
      data f1cs( 38) / -0.0000000000000000000000000000371206D+00 /
      data f1cs( 39) /  0.0000000000000000000000000000085221D+00 /
      data f1cs( 40) / -0.0000000000000000000000000000019633D+00 /
      data f1cs( 41) /  0.0000000000000000000000000000004538D+00 /
      data f1cs( 42) / -0.0000000000000000000000000000001052D+00 /
      data f1cs( 43) /  0.0000000000000000000000000000000245D+00 /

      data f2cs(  1) / -0.03484092538970132330836049733745577D+00 /
      data f2cs(  2) / -0.01668422056779596873246786312278676D+00 /
      data f2cs(  3) /  0.00067529012412377385045207859239727D+00 /
      data f2cs(  4) / -0.00005350666225447013628785577557429D+00 /
      data f2cs(  5) /  0.00000626934217790075267050759431626D+00 /
      data f2cs(  6) / -0.00000095266388019916680677790414293D+00 /
      data f2cs(  7) /  0.00000017456292242509880425504427666D+00 /
      data f2cs(  8) / -0.00000003687954030653093307097646628D+00 /
      data f2cs(  9) /  0.00000000872026777051395264075816938D+00 /
      data f2cs( 10) / -0.00000000226019703919738748530423167D+00 /
      data f2cs( 11) /  0.00000000063246249765250612520444877D+00 /
      data f2cs( 12) / -0.00000000018889118884717869240911480D+00 /
      data f2cs( 13) /  0.00000000005967746729997813372620472D+00 /
      data f2cs( 14) / -0.00000000001980443117372239011196007D+00 /
      data f2cs( 15) /  0.00000000000686413954772103383713264D+00 /
      data f2cs( 16) / -0.00000000000247310193070199106074890D+00 /
      data f2cs( 17) /  0.00000000000092263594549941404196042D+00 /
      data f2cs( 18) / -0.00000000000035523634999261784497297D+00 /
      data f2cs( 19) /  0.00000000000014076049625351591461820D+00 /
      data f2cs( 20) / -0.00000000000005726228499747652794311D+00 /
      data f2cs( 21) /  0.00000000000002386537545413171810106D+00 /
      data f2cs( 22) / -0.00000000000001017141890764597142232D+00 /
      data f2cs( 23) /  0.00000000000000442594531078364424968D+00 /
      data f2cs( 24) / -0.00000000000000196344933049189761979D+00 /
      data f2cs( 25) /  0.00000000000000088688748314810461024D+00 /
      data f2cs( 26) / -0.00000000000000040743345027311546948D+00 /
      data f2cs( 27) /  0.00000000000000019016837215675339859D+00 /
      data f2cs( 28) / -0.00000000000000009009707297478042442D+00 /
      data f2cs( 29) /  0.00000000000000004329211274095668667D+00 /
      data f2cs( 30) / -0.00000000000000002108144465322479526D+00 /
      data f2cs( 31) /  0.00000000000000001039637907026452274D+00 /
      data f2cs( 32) / -0.00000000000000000518891007948931936D+00 /
      data f2cs( 33) /  0.00000000000000000261955324869899371D+00 /
      data f2cs( 34) / -0.00000000000000000133690399951301570D+00 /
      data f2cs( 35) /  0.00000000000000000068941057702931664D+00 /
      data f2cs( 36) / -0.00000000000000000035905362610437250D+00 /
      data f2cs( 37) /  0.00000000000000000018878077255791706D+00 /
      data f2cs( 38) / -0.00000000000000000010016125265594380D+00 /
      data f2cs( 39) /  0.00000000000000000005360725691578228D+00 /
      data f2cs( 40) / -0.00000000000000000002893198974944827D+00 /
      data f2cs( 41) /  0.00000000000000000001574065100202625D+00 /
      data f2cs( 42) / -0.00000000000000000000863027106431206D+00 /
      data f2cs( 43) /  0.00000000000000000000476715602862288D+00 /
      data f2cs( 44) / -0.00000000000000000000265222739998504D+00 /
      data f2cs( 45) /  0.00000000000000000000148582865063866D+00 /
      data f2cs( 46) / -0.00000000000000000000083797235923135D+00 /
      data f2cs( 47) /  0.00000000000000000000047565916422711D+00 /
      data f2cs( 48) / -0.00000000000000000000027169073353112D+00 /
      data f2cs( 49) /  0.00000000000000000000015612738881686D+00 /
      data f2cs( 50) / -0.00000000000000000000009024555078347D+00 /
      data f2cs( 51) /  0.00000000000000000000005246097049119D+00 /
      data f2cs( 52) / -0.00000000000000000000003066450818697D+00 /
      data f2cs( 53) /  0.00000000000000000000001801996250957D+00 /
      data f2cs( 54) / -0.00000000000000000000001064443050752D+00 /
      data f2cs( 55) /  0.00000000000000000000000631942158881D+00 /
      data f2cs( 56) / -0.00000000000000000000000377013812246D+00 /
      data f2cs( 57) /  0.00000000000000000000000225997542918D+00 /
      data f2cs( 58) / -0.00000000000000000000000136100844814D+00 /
      data f2cs( 59) /  0.00000000000000000000000082333232003D+00 /
      data f2cs( 60) / -0.00000000000000000000000050025986091D+00 /
      data f2cs( 61) /  0.00000000000000000000000030526245684D+00 /
      data f2cs( 62) / -0.00000000000000000000000018705164021D+00 /
      data f2cs( 63) /  0.00000000000000000000000011508404393D+00 /
      data f2cs( 64) / -0.00000000000000000000000007108714611D+00 /
      data f2cs( 65) /  0.00000000000000000000000004408065533D+00 /
      data f2cs( 66) / -0.00000000000000000000000002743760867D+00 /
      data f2cs( 67) /  0.00000000000000000000000001714144851D+00 /
      data f2cs( 68) / -0.00000000000000000000000001074768860D+00 /
      data f2cs( 69) /  0.00000000000000000000000000676259777D+00 /
      data f2cs( 70) / -0.00000000000000000000000000426981348D+00 /
      data f2cs( 71) /  0.00000000000000000000000000270500637D+00 /
      data f2cs( 72) / -0.00000000000000000000000000171933331D+00 /
      data f2cs( 73) /  0.00000000000000000000000000109636138D+00 /
      data f2cs( 74) / -0.00000000000000000000000000070132573D+00 /
      data f2cs( 75) /  0.00000000000000000000000000045001784D+00 /
      data f2cs( 76) / -0.00000000000000000000000000028963835D+00 /
      data f2cs( 77) /  0.00000000000000000000000000018697009D+00 /
      data f2cs( 78) / -0.00000000000000000000000000012104646D+00 /
      data f2cs( 79) /  0.00000000000000000000000000007859065D+00 /
      data f2cs( 80) / -0.00000000000000000000000000005116867D+00 /
      data f2cs( 81) /  0.00000000000000000000000000003340627D+00 /
      data f2cs( 82) / -0.00000000000000000000000000002186851D+00 /
      data f2cs( 83) /  0.00000000000000000000000000001435340D+00 /
      data f2cs( 84) / -0.00000000000000000000000000000944523D+00 /
      data f2cs( 85) /  0.00000000000000000000000000000623117D+00 /
      data f2cs( 86) / -0.00000000000000000000000000000412101D+00 /
      data f2cs( 87) /  0.00000000000000000000000000000273208D+00 /
      data f2cs( 88) / -0.00000000000000000000000000000181558D+00 /
      data f2cs( 89) /  0.00000000000000000000000000000120934D+00 /
      data f2cs( 90) / -0.00000000000000000000000000000080737D+00 /
      data f2cs( 91) /  0.00000000000000000000000000000054022D+00 /
      data f2cs( 92) / -0.00000000000000000000000000000036227D+00 /
      data f2cs( 93) /  0.00000000000000000000000000000024348D+00 /
      data f2cs( 94) / -0.00000000000000000000000000000016401D+00 /
      data f2cs( 95) /  0.00000000000000000000000000000011074D+00 /
      data f2cs( 96) / -0.00000000000000000000000000000007497D+00 /
      data f2cs( 97) /  0.00000000000000000000000000000005091D+00 /
      data f2cs( 98) / -0.00000000000000000000000000000003470D+00 /
      data f2cs( 99) /  0.00000000000000000000000000000002377D+00 /

      data g1cs(  1) / -0.3040578798253495954499726682091083D+00 /
      data g1cs(  2) / -0.0566890984597120587731339156118269D+00 /
      data g1cs(  3) /  0.0039046158173275643919984071554082D+00 /
      data g1cs(  4) / -0.0003746075959202260618619339867489D+00 /
      data g1cs(  5) /  0.0000435431556559843679552220840065D+00 /
      data g1cs(  6) / -0.0000057417294453025046561970723475D+00 /
      data g1cs(  7) /  0.0000008282552104502629741937616492D+00 /
      data g1cs(  8) / -0.0000001278245892594642727883913223D+00 /
      data g1cs(  9) /  0.0000000207978352948687884439257529D+00 /
      data g1cs( 10) / -0.0000000035313205921990798042032682D+00 /
      data g1cs( 11) /  0.0000000006210824236308951068631449D+00 /
      data g1cs( 12) / -0.0000000001125215474446292649336987D+00 /
      data g1cs( 13) /  0.0000000000209088917684421605267019D+00 /
      data g1cs( 14) / -0.0000000000039715831737681727689158D+00 /
      data g1cs( 15) /  0.0000000000007690431314272089939005D+00 /
      data g1cs( 16) / -0.0000000000001514696742731613519826D+00 /
      data g1cs( 17) /  0.0000000000000302892146552359684119D+00 /
      data g1cs( 18) / -0.0000000000000061399703834708825400D+00 /
      data g1cs( 19) /  0.0000000000000012600605829510933553D+00 /
      data g1cs( 20) / -0.0000000000000002615029250939483683D+00 /
      data g1cs( 21) /  0.0000000000000000548278844891796821D+00 /
      data g1cs( 22) / -0.0000000000000000116038182129526571D+00 /
      data g1cs( 23) /  0.0000000000000000024771654107129795D+00 /
      data g1cs( 24) / -0.0000000000000000005330672753223389D+00 /
      data g1cs( 25) /  0.0000000000000000001155666075598465D+00 /
      data g1cs( 26) / -0.0000000000000000000252280547744957D+00 /
      data g1cs( 27) /  0.0000000000000000000055429038550786D+00 /
      data g1cs( 28) / -0.0000000000000000000012252208421297D+00 /
      data g1cs( 29) /  0.0000000000000000000002723664318684D+00 /
      data g1cs( 30) / -0.0000000000000000000000608707831422D+00 /
      data g1cs( 31) /  0.0000000000000000000000136724874476D+00 /
      data g1cs( 32) / -0.0000000000000000000000030856626806D+00 /
      data g1cs( 33) /  0.0000000000000000000000006995212319D+00 /
      data g1cs( 34) / -0.0000000000000000000000001592587569D+00 /
      data g1cs( 35) /  0.0000000000000000000000000364051056D+00 /
      data g1cs( 36) / -0.0000000000000000000000000083539465D+00 /
      data g1cs( 37) /  0.0000000000000000000000000019240303D+00 /
      data g1cs( 38) / -0.0000000000000000000000000004446816D+00 /
      data g1cs( 39) /  0.0000000000000000000000000001031182D+00 /
      data g1cs( 40) / -0.0000000000000000000000000000239887D+00 /
      data g1cs( 41) /  0.0000000000000000000000000000055976D+00 /
      data g1cs( 42) / -0.0000000000000000000000000000013100D+00 /
      data g1cs( 43) /  0.0000000000000000000000000000003074D+00 /
      data g1cs( 44) / -0.0000000000000000000000000000000723D+00 /

      data g2cs(  1) / -0.1211802894731646263541834046858267D+00 /
      data g2cs(  2) / -0.0316761386394950286701407923505610D+00 /
      data g2cs(  3) /  0.0013383199778862680163819429492182D+00 /
      data g2cs(  4) / -0.0000895511011392252425531905069518D+00 /
      data g2cs(  5) /  0.0000079155562961718213115249467924D+00 /
      data g2cs(  6) / -0.0000008438793322241520181418982080D+00 /
      data g2cs(  7) /  0.0000001029980425677530146647227274D+00 /
      data g2cs(  8) / -0.0000000139295750605183835795834444D+00 /
      data g2cs(  9) /  0.0000000020422703959875980400677594D+00 /
      data g2cs( 10) / -0.0000000003196534694206427035434752D+00 /
      data g2cs( 11) /  0.0000000000528147832657267698615312D+00 /
      data g2cs( 12) / -0.0000000000091339554672671033735289D+00 /
      data g2cs( 13) /  0.0000000000016426251238967760444819D+00 /
      data g2cs( 14) / -0.0000000000003055897039322660002410D+00 /
      data g2cs( 15) /  0.0000000000000585655825785779717892D+00 /
      data g2cs( 16) / -0.0000000000000115229197730940120563D+00 /
      data g2cs( 17) /  0.0000000000000023209469119988537310D+00 /
      data g2cs( 18) / -0.0000000000000004774355834177535025D+00 /
      data g2cs( 19) /  0.0000000000000001000996765800180573D+00 /
      data g2cs( 20) / -0.0000000000000000213533778082256704D+00 /
      data g2cs( 21) /  0.0000000000000000046277190777367671D+00 /
      data g2cs( 22) / -0.0000000000000000010175807410227657D+00 /
      data g2cs( 23) /  0.0000000000000000002267657399884672D+00 /
      data g2cs( 24) / -0.0000000000000000000511630776076426D+00 /
      data g2cs( 25) /  0.0000000000000000000116767014913108D+00 /
      data g2cs( 26) / -0.0000000000000000000026935427672470D+00 /
      data g2cs( 27) /  0.0000000000000000000006275665841146D+00 /
      data g2cs( 28) / -0.0000000000000000000001475880557531D+00 /
      data g2cs( 29) /  0.0000000000000000000000350145314739D+00 /
      data g2cs( 30) / -0.0000000000000000000000083757732152D+00 /
      data g2cs( 31) /  0.0000000000000000000000020191815152D+00 /
      data g2cs( 32) / -0.0000000000000000000000004903567705D+00 /
      data g2cs( 33) /  0.0000000000000000000000001199123348D+00 /
      data g2cs( 34) / -0.0000000000000000000000000295170610D+00 /
      data g2cs( 35) /  0.0000000000000000000000000073113112D+00 /
      data g2cs( 36) / -0.0000000000000000000000000018217843D+00 /
      data g2cs( 37) /  0.0000000000000000000000000004565148D+00 /
      data g2cs( 38) / -0.0000000000000000000000000001150151D+00 /
      data g2cs( 39) /  0.0000000000000000000000000000291267D+00 /
      data g2cs( 40) / -0.0000000000000000000000000000074125D+00 /
      data g2cs( 41) /  0.0000000000000000000000000000018953D+00 /
      data g2cs( 42) / -0.0000000000000000000000000000004868D+00 /
      data g2cs( 43) /  0.0000000000000000000000000000001256D+00 /
      data g2cs( 44) / -0.0000000000000000000000000000000325D+00 /

      data g3cs(  1) / -0.0280574367809472928402815264335299D+00 /
      data g3cs(  2) / -0.0137271597162236975409100508089556D+00 /
      data g3cs(  3) /  0.0002894032638760296027448941273751D+00 /
      data g3cs(  4) / -0.0000114129239391197145908743622517D+00 /
      data g3cs(  5) /  0.0000006813965590726242997720207302D+00 /
      data g3cs(  6) / -0.0000000547952289604652363669058052D+00 /
      data g3cs(  7) /  0.0000000055207429918212529109406521D+00 /
      data g3cs(  8) / -0.0000000006641464199322920022491428D+00 /
      data g3cs(  9) /  0.0000000000922373663487041108564960D+00 /
      data g3cs( 10) / -0.0000000000144299088886682862611718D+00 /
      data g3cs( 11) /  0.0000000000024963904892030710248705D+00 /
      data g3cs( 12) / -0.0000000000004708240675875244722971D+00 /
      data g3cs( 13) /  0.0000000000000957217659216759988140D+00 /
      data g3cs( 14) / -0.0000000000000207889966095809030537D+00 /
      data g3cs( 15) /  0.0000000000000047875099970877431627D+00 /
      data g3cs( 16) / -0.0000000000000011619070583377173759D+00 /
      data g3cs( 17) /  0.0000000000000002956508969267836974D+00 /
      data g3cs( 18) / -0.0000000000000000785294988256492025D+00 /
      data g3cs( 19) /  0.0000000000000000216922264368256612D+00 /
      data g3cs( 20) / -0.0000000000000000062113515831676342D+00 /
      data g3cs( 21) /  0.0000000000000000018384568838450977D+00 /
      data g3cs( 22) / -0.0000000000000000005610887482137276D+00 /
      data g3cs( 23) /  0.0000000000000000001761862805280062D+00 /
      data g3cs( 24) / -0.0000000000000000000568111050541451D+00 /
      data g3cs( 25) /  0.0000000000000000000187786279582313D+00 /
      data g3cs( 26) / -0.0000000000000000000063531694151124D+00 /
      data g3cs( 27) /  0.0000000000000000000021968802368238D+00 /
      data g3cs( 28) / -0.0000000000000000000007754666550395D+00 /
      data g3cs( 29) /  0.0000000000000000000002791018356581D+00 /
      data g3cs( 30) / -0.0000000000000000000001023178525247D+00 /
      data g3cs( 31) /  0.0000000000000000000000381693403919D+00 /
      data g3cs( 32) / -0.0000000000000000000000144767895606D+00 /
      data g3cs( 33) /  0.0000000000000000000000055779512634D+00 /
      data g3cs( 34) / -0.0000000000000000000000021817239071D+00 /
      data g3cs( 35) /  0.0000000000000000000000008656646309D+00 /
      data g3cs( 36) / -0.0000000000000000000000003482157895D+00 /
      data g3cs( 37) /  0.0000000000000000000000001419188130D+00 /
      data g3cs( 38) / -0.0000000000000000000000000585714314D+00 /
      data g3cs( 39) /  0.0000000000000000000000000244660482D+00 /
      data g3cs( 40) / -0.0000000000000000000000000103387099D+00 /
      data g3cs( 41) /  0.0000000000000000000000000044177299D+00 /
      data g3cs( 42) / -0.0000000000000000000000000019080079D+00 /
      data g3cs( 43) /  0.0000000000000000000000000008326038D+00 /
      data g3cs( 44) / -0.0000000000000000000000000003669553D+00 /
      data g3cs( 45) /  0.0000000000000000000000000001632875D+00 /
      data g3cs( 46) / -0.0000000000000000000000000000733357D+00 /
      data g3cs( 47) /  0.0000000000000000000000000000332327D+00 /
      data g3cs( 48) / -0.0000000000000000000000000000151906D+00 /
      data g3cs( 49) /  0.0000000000000000000000000000070020D+00 /
      data g3cs( 50) / -0.0000000000000000000000000000032539D+00 /
      data g3cs( 51) /  0.0000000000000000000000000000015240D+00 /
      data g3cs( 52) / -0.0000000000000000000000000000007193D+00 /
      data g3cs( 53) /  0.0000000000000000000000000000003420D+00 /
      data g3cs( 54) / -0.0000000000000000000000000000001638D+00 /
      data g3cs( 55) /  0.0000000000000000000000000000000790D+00 /
      data g3cs( 56) / -0.0000000000000000000000000000000383D+00 /

      data nf1 / 0 /
      data nf2 / 0 /
      data ng1 / 0 /
      data ng2 / 0 /
      data ng3 / 0 /
      data xbig / 0.0D+00 /
      data xbnd / 0.0D+00 /
      data xbndg / 0.0D+00 /
      data xmaxf / 0.0D+00 /
      data xmaxg / 0.0D+00 /

      if ( nf1 .eq. 0 ) then
        tol = 0.1D+00 * r8_mach ( 3 )
        nf1 = r8_inits ( f1cs, 43, tol )
        nf2 = r8_inits ( f2cs, 99, tol )
        ng1 = r8_inits ( g1cs, 44, tol )
        ng2 = r8_inits ( g2cs, 44, tol )
        ng3 = r8_inits ( g3cs, 56, tol )
        xbig = dsqrt ( 1.0D+00 / r8_mach ( 3 ) )
        xmaxf = dexp ( dmin1 ( - dlog ( r8_mach ( 1 ) ), 
     &    dlog ( r8_mach ( 2 ) ) ) - 0.01D+00 )
        xmaxg = 1.0D+00 / dsqrt ( r8_mach ( 1 ) )
        xbnd = dsqrt ( 50.0D+00 )
        xbndg = dsqrt ( 200.0D+00 )
      end if

      if ( x .lt. 4.0D+00 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'R8_SIFG - Fatal error!'
        write ( *, '(a)' ) '  Approximation invalid for X < 4.'
        stop
      else if ( x .le. xbnd ) then
        f = ( 1.0D+00 
     &    + r8_csevl ( ( 1.0D+00 / x / x - 0.04125D+00 )
     &    / 0.02125D+00, f1cs, nf1 ) ) / x
        g = ( 1.0D+00 
     &    + r8_csevl ( ( 1.0D+00 / x / x - 0.04125D+00 )
     &    / 0.02125D+00, g1cs, ng1 ) ) / x / x
      else if ( x .le. xbig ) then
        f = ( 1.0D+00 
     &    + r8_csevl ( 100.D+00 / x / x - 1.0D+00, f2cs, nf2 ) ) / x
        if ( x .le. xbndg ) then 
          g = ( 1.0D+00
     &      + r8_csevl ( ( 10000.0D+00 / x / x - 125.0D+00 ) 
     &      / 75.0D+00, g2cs, ng2 ) ) / x / x
        else
          g = ( 1.0D+00 
     &      + r8_csevl ( 400.0D+00 / x / x - 1.0D+00, g3cs, ng3 ) ) 
     &      / x / x
        end if
      else
        if ( x .lt. xmaxf ) then
          f = 1.0D+00 / x
        else
          f = 0.0D+00
        end if
        if ( x .lt. xmaxg ) then
          g = 1.0D+00 / x / x
        else
          g = 0.0D+00
        end if
      end if

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
      subroutine r8vec_uniform ( n, a, b, seed, r )

c*********************************************************************72
c
cc R8VEC_UNIFORM returns a scaled pseudorandom R8VEC.
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
c    29 January 2005
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
c    Bennett Fox,
c    Algorithm 647:
c    Implementation and Relative Efficiency of Quasirandom
c    Sequence Generators,
c    ACM Transactions on Mathematical Software,
c    Volume 12, Number 4, pages 362-376, 1986.
c
c    Peter Lewis, Allen Goodman, James Miller,
c    A Pseudo-Random Number Generator for the System/360,
c    IBM Systems Journal,
c    Volume 8, pages 136-143, 1969.
c
c  Parameters:
c
c    Input, integer M, the number of entries in the vector.
c
c    Input, double precision A, B, the lower and upper limits.
c
c    Input/output, integer SEED, the "seed" value, which should NOT be 0.
c    On output, SEED has been updated.
c
c    Output, double precision R(N), the vector of pseudorandom values.
c
      implicit none

      integer n

      double precision a
      double precision b
      integer i
      integer k
      integer seed
      double precision r(n)

      if ( seed .eq. 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'R8VEC_UNIFORM - Fatal error!'
        write ( *, '(a)' ) '  Input value of SEED = 0.'
        stop
      end if

      do i = 1, n

        k = seed / 127773

        seed = 16807 * ( seed - k * 127773 ) - k * 2836

        if ( seed .lt. 0 ) then
          seed = seed + 2147483647
        end if

        r(i) = a + ( b - a ) * dble ( seed ) * 4.656612875D-10

      end do

      return
      end
      function s_eqi ( s1, s2 )

c*********************************************************************72
c
cc S_EQI is a case insensitive comparison of two strings for equality.
c
c  Example:
c
c    S_EQI ( 'Anjana', 'ANJANA' ) is TRUE.
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
c    Input, character*(*) S1, S2, the strings to compare.
c
c    Output, logical S_EQI, the result of the comparison.
c
      implicit none

      character c1
      character c2
      integer i
      integer lenc
      logical s_eqi
      character*(*) s1
      integer s1_length
      character*(*) s2
      integer s2_length

      s1_length = len ( s1 )
      s2_length = len ( s2 )
      lenc = min ( s1_length, s2_length )

      s_eqi = .false.

      do i = 1, lenc

        c1 = s1(i:i)
        c2 = s2(i:i)
        call ch_cap ( c1 )
        call ch_cap ( c2 )

        if ( c1 .ne. c2 ) then
          return
        end if

      end do

      do i = lenc + 1, s1_length
        if ( s1(i:i) .ne. ' ' ) then
          return
        end if
      end do

      do i = lenc + 1, s2_length
        if ( s2(i:i) .ne. ' ' ) then
          return
        end if
      end do

      s_eqi = .true.

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
