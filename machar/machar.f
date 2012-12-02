      subroutine r4_machar ( ibeta, it, irnd, ngrd, machep, negep,
     &  iexp, minexp, maxexp, eps, epsneg, xmin, xmax )

c*********************************************************************72
c
cc R4_MACHAR determines single precision machine constants.
c
c  Discussion:
c
c    This routine determines the parameters of the single precision
c    floating-point arithmetic system.  The determination of the first 
c    three uses an extension of an algorithm due to Malcolm, 
c    incorporating some of the improvements suggested by Gentleman and 
c    Marovich.  
c
c    This routine appeared as ACM algorithm 665.
c
c    An earlier version of this program was published in Cody and Waite.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    05 May 2010
c
c  Author:
c
c    Original FORTRAN77 version by William Cody.
c    This FORTRAN77 version by John Burkardt.
c
c  Reference:
c
c    William Cody,
c    ACM Algorithm 665, MACHAR, a subroutine to dynamically determine 
c    machine parameters,
c    ACM Transactions on Mathematical Software,
c    Volume 14, Number 4, pages 303-311, 1988.
c
c    William Cody, William Waite,
c    Software Manual for the Elementary Functions,
c    Prentice Hall, 1980.
c
c    Morven Gentleman, Scott Marovich,
c    Communications of the ACM,
c    Volume 17, pages 276-277, 1974.
c
c    Michael Malcolm,
c    Communications of the ACM,
c    Volume 15, pages 949-951, 1972.
c
c  Parameters:
c
c    Output, integer IBETA, the radix for the floating-point 
c    representation.
c
c    Output, integer IT, the number of base IBETA digits 
c    in the floating-point significand.
c
c    Output, integer IRND:
c    0, if floating-point addition chops.
c    1, if floating-point addition rounds, but not in the IEEE style.
c    2, if floating-point addition rounds in the IEEE style.
c    3, if floating-point addition chops, and there is partial underflow.
c    4, if floating-point addition rounds, but not in the IEEE style, and 
c      there is partial underflow.
c    5, if floating-point addition rounds in the IEEE style, and there is 
c      partial underflow.
c
c    Output, integer NGRD, the number of guard digits for 
c    multiplication with truncating arithmetic.  It is
c    0, if floating-point arithmetic rounds, or if it truncates and only 
c      IT base IBETA digits participate in the post-normalization shift of the
c      floating-point significand in multiplication;
c    1, if floating-point arithmetic truncates and more than IT base IBETA
c      digits participate in the post-normalization shift of the floating-point
c      significand in multiplication.
c
c    Output, integer MACHEP, the largest negative integer such that
c      1.0 + real ( IBETA ) ^ MACHEP /= 1.0, 
c    except that MACHEP is bounded below by - ( IT + 3 ).
c
c    Output, integer NEGEPS, the largest negative integer such that
c      1.0 - real ( IBETA ) ^ NEGEPS /= 1.0, 
c    except that NEGEPS is bounded below by - ( IT + 3 ).
c
c    Output, integer IEXP, the number of bits (decimal places 
c    if IBETA = 10) reserved for the representation of the exponent (including 
c    the bias or sign) of a floating-point number.
c
c    Output, integer MINEXP, the largest in magnitude negative 
c    integer such that
c      real ( IBETA ) ^ MINEXP 
c    is positive and normalized.
c
c    Output, integer MAXEXP, the smallest positive power of 
c    BETA that overflows.
c
c    Output, real EPS, the smallest positive floating-point 
c    number such that  
c      1.0 + EPS /= 1.0. 
c    in particular, if either IBETA = 2  or IRND = 0, 
c      EPS = real ( IBETA ) ^ MACHEP.
c    Otherwise,  
c      EPS = ( real ( IBETA ) ^ MACHEP ) / 2.
c
c    Output, real EPSNEG, a small positive floating-point number 
c    such that
c      1.0 - EPSNEG /= 1.0. 
c    In particular, if IBETA = 2 or IRND = 0, 
c      EPSNEG = real ( IBETA ) ^ NEGEPS.
c    Otherwise,  
c      EPSNEG = ( real ( IBETA ) ^ NEGEPS ) / 2.  
c    Because NEGEPS is bounded below by - ( IT + 3 ), EPSNEG might not be the
c    smallest number that can alter 1.0 by subtraction.
c
c    Output, real XMIN, the smallest non-vanishing normalized 
c    floating-point power of the radix:
c      XMIN = real ( IBETA ) ^ MINEXP
c
c    Output, real XMAX, the largest finite floating-point number.
c    In particular,
c      XMAX = ( 1.0 - EPSNEG ) * real ( IBETA ) ^ MAXEXP
c    On some machines, the computed value of XMAX will be only the second, 
c    or perhaps third, largest number, being too small by 1 or 2 units in 
c    the last digit of the significand.
c
      implicit none

      real a
      real b
      real beta
      real betah
      real betain
      real eps
      real epsneg
      integer i
      integer ibeta
      integer iexp
      integer irnd
      integer it
      integer itemp
      integer iz
      integer j
      integer k
      integer machep
      integer maxexp
      integer minexp
      integer mx
      integer negep
      integer ngrd
      integer nxres
      real one
      real t
      real temp
      real temp1
      real tempa
      real two
      real xmax
      real xmin
      real y
      real z
      real zero

      one = real ( 1 )
      two = one + one
      zero = one - one
c
c  Determine IBETA, BETA ala Malcolm.
c
      a = one

10    continue

      a = a + a
      temp = a + one
      temp1 = temp - a
      if ( temp1 - one .eq. zero ) then
        go to 10
      end if

      b = one

20    continue

      b = b + b
      temp = a + b
      itemp = int ( temp - a )
      if ( itemp .eq. 0 ) then
        go to 20
      end if

      ibeta = itemp
      beta = real ( ibeta )
c
c  Determine IT, IRND.
c
      it = 0
      b = one

30    continue

      it = it + 1
      b = b * beta
      temp = b + one
      temp1 = temp-b
      if ( temp1 - one .eq. zero ) then
        go to 30
      end if

      irnd = 0
      betah = beta / two
      temp = a + betah
      if ( temp - a .ne. zero ) then
        irnd = 1
      end if
      tempa = a + beta
      temp = tempa + betah
      if ( ( irnd .eq. 0 ) .and. ( temp - tempa .ne. zero ) ) then
        irnd = 2
      end if
c
c  Determine NEGEP, EPSNEG.
c
      negep = it + 3
      betain = one / beta
      a = one
      do i = 1, negep
        a = a * betain
      end do
      b = a

40    continue

      temp = one - a

      if ( temp - one .eq. zero ) then
        a = a * beta
        negep = negep - 1
        go to 40
      end if

      negep = - negep
      epsneg = a
      if ( ( ibeta .ne. 2 ) .and. ( irnd .ne. 0 ) ) then
        a = ( a * ( one + a ) ) / two
        temp = one - a
        if ( temp - one .ne. zero ) then
          epsneg = a
        end if
      end if
c
c  Determine MACHEP, EPS.
c
      machep = - it - 3
      a = b

50    continue

      temp = one + a

      if ( temp - one .eq. zero ) then
        a = a * beta
        machep = machep + 1
        go to 50
      end if

      eps = a
      temp = tempa + beta * ( one + eps )

      if ( ( ibeta .ne. 2 ) .and. ( irnd .ne. 0 ) ) then
        a = ( a * ( one + a ) ) / two
        temp = one + a
        if ( temp - one .ne. zero ) then
          eps = a
        end if
      end if
c
c  Determine NGRD.
c
      ngrd = 0
      temp = one + eps
      if ( ( irnd .eq. 0 ) .and. ( temp * one - one .ne. zero ) ) then
        ngrd = 1
      endif
c
c  Determine IEXP, MINEXP, XMIN.
c
c  Loop to determine largest I and K = 2^I such that
c         (1/BETA) ^ (2^(I))
c  does not underflow.
c  Exit from loop is signaled by an underflow.
c
      i = 0
      k = 1
      z = betain
      t = one + eps
      nxres = 0

60    continue

      y = z
      z = y * y
c
c  Check for underflow here.
c
      a = z * one
      temp = z * t

      if ( ( a + a .ne. zero ) .and. ( abs ( z ) .lt. y ) ) then

        temp1 = temp * betain

        if ( temp1 * beta .ne. z ) then
          i = i + 1
          k = k + k
          go to 60
        end if

      end if
c
c  This segment is for decimal machines only.
c
      if ( ibeta .eq. 10 ) then

        iexp = 2
        iz = ibeta

70      continue

        if ( iz .le. k ) then
          iz = iz * ibeta
          iexp = iexp + 1
          go to 70
        end if

        mx = iz + iz - 1

      else

        iexp = i + 1
        mx = k + k

      end if
c
c  Loop to determine MINEXP, XMIN.
c  Exit from loop is signaled by an underflow.
c
80    continue

      xmin = y
      y = y * betain
c
c  Check for underflow here.
c
      a = y * one
      temp = y * t

      if ( ( ( a + a ) .ne. zero ) .and. ( abs ( y ) .lt. xmin ) ) then

        k = k + 1
        temp1 = temp * betain

        if ( temp1 * beta .ne. y ) then
          go to 80
        end if

        nxres = 3
        xmin = y

      end if

      minexp = - k
c
c  Determine MAXEXP, XMAX.
c
      if ( ( mx .le. k + k - 3 ) .and. ( ibeta .ne. 10 ) ) then
        mx = mx + mx
        iexp = iexp + 1
      end if

      maxexp = mx + minexp
c
c  Adjust IRND to reflect partial underflow.
c
      irnd = irnd + nxres
c
c  Adjust for IEEE-style machines.
c
      if ( ( irnd .eq. 2 ) .or. ( irnd .eq. 5 ) ) then
        maxexp = maxexp - 2
      end if
c
c  Adjust for non-IEEE machines with partial underflow.
c
      if ( ( irnd .eq. 3 ) .or. ( irnd .eq. 4 ) ) then
        maxexp = maxexp - it
      end if
c
c  Adjust for machines with implicit leading bit in binary
c  significand, and machines with radix point at extreme
c  right of significand.
c
      i = maxexp + minexp
      if ( ( ibeta .eq. 2 ) .and. ( i .eq. 0 ) ) then
        maxexp = maxexp - 1
      end if
      if ( 20 .lt. i ) then
        maxexp = maxexp - 1
      end if
      if ( a .ne. y ) then
        maxexp = maxexp - 2
      end if
      xmax = one - epsneg
      if ( xmax * one .ne. xmax ) then
        xmax = one - beta * epsneg
      end if
      xmax = xmax / ( beta * beta * beta * xmin )
      i = maxexp + minexp + 3

      do j = 1, i
        if ( ibeta .eq. 2 ) then
          xmax = xmax + xmax
        else
          xmax = xmax * beta
        end if
      end do

      return
      end
      subroutine r8_machar ( ibeta, it, irnd, ngrd, machep, negep,
     &  iexp, minexp, maxexp, eps, epsneg, xmin, xmax )

c*********************************************************************72
c
cc R8_MACHAR determines double precision machine constants.
c
c  Discussion:
c
c    This routine determines the parameters of the double precision
c    floating-point arithmetic system.  The determination of the first 
c    three uses an extension of an algorithm due to Malcolm, 
c    incorporating some of the improvements suggested by Gentleman and 
c    Marovich.  
c
c    This routine appeared as ACM algorithm 665.
c
c    An earlier version of this program was published in Cody and Waite.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    28 January 2002
c
c  Author:
c
c    Original FORTRAN77 version by William Cody.
c    This FORTRAN77 version by John Burkardt.
c
c  Reference:
c
c    William Cody,
c    ACM Algorithm 665, MACHAR, a subroutine to dynamically determine 
c    machine parameters,
c    ACM Transactions on Mathematical Software,
c    Volume 14, Number 4, pages 303-311, 1988.
c
c    William Cody, William Waite,
c    Software Manual for the Elementary Functions,
c    Prentice Hall, 1980.
c
c    Morven Gentleman, Scott Marovich,
c    Communications of the ACM,
c    Volume 17, pages 276-277, 1974.
c
c    Michael Malcolm,
c    Communications of the ACM,
c    Volume 15, pages 949-951, 1972.
c
c  Parameters:
c
c    Output, integer IBETA, the radix for the floating-point
c    representation.
c
c    Output, integer IT, the number of base IBETA digits in 
c    the floating-point significand.
c
c    Output, integer IRND:
c    0, if floating-point addition chops.
c    1, if floating-point addition rounds, but not in the IEEE style.
c    2, if floating-point addition rounds in the IEEE style.
c    3, if floating-point addition chops, and there is partial underflow.
c    4, if floating-point addition rounds, but not in the IEEE style, and 
c      there is partial underflow.
c    5, if floating-point addition rounds in the IEEE style, and there is 
c      partial underflow.
c
c    Output, integer NGRD, the number of guard digits for 
c    multiplication with truncating arithmetic.  It is
c    0, if floating-point arithmetic rounds, or if it truncates and only 
c      IT base IBETA digits participate in the post-normalization shift of the
c      floating-point significand in multiplication;
c    1, if floating-point arithmetic truncates and more than IT base IBETA
c      digits participate in the post-normalization shift of the floating-point
c      significand in multiplication.
c
c    Output, integer MACHEP, the largest negative integer 
c    such that
c      1.0 < 1.0 + dble ( IBETA ) ^ MACHEP, 
c    except that MACHEP is bounded below by - ( IT + 3 ).
c
c    Output, integer NEGEPS, the largest negative integer 
c    such that
c      1.0 - dble ( IBETA ) ^ NEGEPS < 1.0, 
c    except that NEGEPS is bounded below by - ( IT + 3 ).
c
c    Output, integer IEXP, the number of bits (decimal places 
c    if IBETA = 10) reserved for the representation of the exponent (including
c    the bias or sign) of a floating-point number.
c
c    Output, integer MINEXP, the largest in magnitude negative
c    integer such that
c      dble ( IBETA ) ^ MINEXP 
c    is positive and normalized.
c
c    Output, integer MAXEXP, the smallest positive power of
c    BETA that overflows.
c
c    Output, double precision EPS, the smallest positive floating-point number
c    such that  
c      1.0 + EPS /= 1.0. 
c    in particular, if either IBETA = 2  or IRND = 0, 
c      EPS = dble ( IBETA ) ^ MACHEP.
c    Otherwise,  
c      EPS = ( dble ( IBETA ) ^ MACHEP ) / 2.
c
c    Output, double precision EPSNEG, a small positive floating-point number
c    such that
c      1.0 - EPSNEG < 1.0. 
c    In particular, if IBETA = 2 or IRND = 0, 
c      EPSNEG = dble ( IBETA ) ^ NEGEPS.
c    Otherwise,  
c      EPSNEG = ( dble ( IBETA ) ^ NEGEPS ) / 2.  
c    Because NEGEPS is bounded below by - ( IT + 3 ), EPSNEG might not be the
c    smallest number that can alter 1.0 by subtraction.
c
c    Output, double precision XMIN, the smallest non-vanishing normalized
c    floating-point power of the radix:
c      XMIN = dble ( IBETA ) ^ MINEXP
c
c    Output, double precision XMAX, the largest finite floating-point number.
c    In particular,
c      XMAX = ( 1.0 - EPSNEG ) * dble ( IBETA ) ^ MAXEXP
c    On some machines, the computed value of XMAX will be only the second, 
c    or perhaps third, largest number, being too small by 1 or 2 units in 
c    the last digit of the significand.
c
      implicit none

      double precision a
      double precision b
      double precision beta
      double precision betah
      double precision betain
      double precision eps
      double precision epsneg
      integer i
      integer ibeta
      integer iexp
      integer irnd
      integer it
      integer itemp
      integer iz
      integer j
      integer k
      integer machep
      integer maxexp
      integer minexp
      integer mx
      integer negep
      integer ngrd
      integer nxres
      double precision one
      double precision t
      double precision temp
      double precision temp1
      double precision tempa
      double precision two
      double precision xmax
      double precision xmin
      double precision y
      double precision z
      double precision zero

      one = dble ( 1 )
      two = one + one
      zero = one - one
c
c  Determine IBETA, BETA ala Malcolm.
c
      a = one
   10 continue
      a = a + a
      temp = a+one
      temp1 = temp-a
      if (temp1-one .eq. zero) then
        go to 10
      end if

      b = one
   20 continue
      b = b + b
      temp = a+b
      itemp = int(temp-a)
      if (itemp .eq. 0) then
        go to 20
      end if

      ibeta = itemp
      beta = dble ( ibeta )
c
c  Determine IT, IRND.
c
      it = 0
      b = one
  100 it = it + 1
         b = b * beta
         temp = b+one
         temp1 = temp-b
         if (temp1-one .eq. zero) go to 100
      irnd = 0
      betah = beta / two
      temp = a+betah
      if (temp-a .ne. zero) irnd = 1
      tempa = a + beta
      temp = tempa+betah
      if ((irnd .eq. 0) .and. (temp-tempa .ne. zero)) irnd = 2
c
c  Determine NEGEP, EPSNEG.
c
      negep = it + 3
      betain = one / beta
      a = one
      do i = 1, negep
         a = a * betain
      end do
      b = a
  210 temp = one-a
         if (temp-one .ne. zero) go to 220
         a = a * beta
         negep = negep - 1
      go to 210
  220 negep = -negep
      epsneg = a
      if ((ibeta .eq. 2) .or. (irnd .eq. 0)) go to 300
      a = (a*(one+a)) / two
      temp = one-a
      if (temp-one .ne. zero) epsneg = a
c
c  Determine MACHEP, EPS.
c
  300 machep = -it - 3
      a = b
  310 temp = one+a
         if (temp-one .ne. zero) go to 320
         a = a * beta
         machep = machep + 1
      go to 310
  320 eps = a
      temp = tempa+beta*(one+eps)
      if ((ibeta .eq. 2) .or. (irnd .eq. 0)) go to 350
      a = (a*(one+a)) / two
      temp = one+a
      if (temp-one .ne. zero) eps = a
c
c  Determine NGRD.
c
  350 ngrd = 0
      temp = one+eps
      if ((irnd .eq. 0) .and. (temp*one-one .ne. zero)) ngrd = 1
c
c  Determine IEXP, MINEXP, XMIN.
c
c  Loop to determine largest I and K = 2^I such that
c         (1/BETA) ^ (2^(I))
c  does not underflow.
c  Exit from loop is signaled by an underflow.
c
      i = 0
      k = 1
      z = betain
      t = one + eps
      nxres = 0
  400 y = z
         z = y * y
c
c  Check for underflow here.
c
         a = z * one
         temp = z * t
         if ((a+a .eq. zero) .or. (abs(z) .ge. y)) go to 410
         temp1 = temp * betain
         if (temp1*beta .eq. z) go to 410
         i = i + 1
         k = k + k
      go to 400
  410 if (ibeta .eq. 10) go to 420
      iexp = i + 1
      mx = k + k
      go to 450
c
c  This segment is for decimal machines only.
c
  420 iexp = 2
      iz = ibeta
  430 if (k .lt. iz) go to 440
         iz = iz * ibeta
         iexp = iexp + 1
      go to 430
  440 mx = iz + iz - 1
c
c  Loop to determine MINEXP, XMIN.
c  Exit from loop is signaled by an underflow.
c
  450 xmin = y
         y = y * betain
c
c  Check for underflow here.
c
         a = y * one
         temp = y * t
         if (((a+a) .eq. zero) .or. (abs(y) .ge. xmin)) go to 460
         k = k + 1
         temp1 = temp * betain
         if (temp1*beta .ne. y) go to 450
      nxres = 3
      xmin = y
  460 minexp = -k
c
c  Determine MAXEXP, XMAX.
c
      if ((mx .gt. k+k-3) .or. (ibeta .eq. 10)) go to 500
      mx = mx + mx
      iexp = iexp + 1
  500 maxexp = mx + minexp
c
c  Adjust IRND to reflect partial underflow.
c
      irnd = irnd + nxres
c
c  Adjust for IEEE-style machines.
c
      if ((irnd .eq. 2) .or. (irnd .eq. 5)) maxexp = maxexp - 2
c
c  Adjust for non-IEEE machines with partial underflow.
c
      if ((irnd .eq. 3) .or. (irnd .eq. 4)) maxexp = maxexp - it
c
c  Adjust for machines with implicit leading bit in binary
c  significand, and machines with radix point at extreme
c  right of significand.
c
      i = maxexp + minexp
      if ((ibeta .eq. 2) .and. (i .eq. 0)) maxexp = maxexp - 1
      if (i .gt. 20) maxexp = maxexp - 1
      if (a .ne. y) maxexp = maxexp - 2
      xmax = one - epsneg
      if (xmax*one .ne. xmax) xmax = one - beta * epsneg
      xmax = xmax / (beta * beta * beta * xmin)
      i = maxexp + minexp + 3

      do j = 1, i
          if (ibeta .eq. 2) xmax = xmax + xmax
          if (ibeta .ne. 2) xmax = xmax * beta
      end do

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
