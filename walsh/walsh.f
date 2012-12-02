      subroutine ffwt ( n, x )

c*********************************************************************72
c
cc FFWT performs an in-place fast Walsh transform.
c
c  Discussion:
c
c    This routine performs a fast Walsh transform on an input series X
c    leaving the transformed results in X. 
c    X is dimensioned N, which must be a power of 2.
c    The results of this Walsh transform are in sequency order.
c
c    The output sequence could be normalized by dividing by N.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    11 March 2011
c
c  Author:
c
c    Ken Beauchamp
c
c  Reference:
c
c    Ken Beauchamp,
c    Walsh functions and their applications,
c    Academic Press, 1975,
c    ISBN: 0-12-084050-2,
c    LC: QA404.5.B33.
c
c  Parameters:
c
c    Input, integer N, the number of items in X.
c    N must be a power of 2.
c
c    Input/output, double precision X(N), the data to be transformed.
c
      implicit none

      integer n

      double precision hold
      integer i
      integer i4_log_2
      integer ii
      integer j
      integer j2
      integer js
      integer k
      integer l
      integer m
      integer mw
      integer mw1
      integer nw
      integer nz
      integer nz2
      integer nzi
      integer nzn
      integer two_power(24)
      double precision x(n)
      double precision z

      m = i4_log_2 ( n )

      do i = 1, m
        two_power(i) = 2**( m - i )
      end do

      do l = 1, m

        nz = 2**( l - 1 )
        nzi = 2 * nz
        nzn = n / nzi
        nz2 = nz / 2
        if ( nz2 .eq. 0 ) then
          nz2 = 1
        end if

        do i = 1, nzn

          js = ( i - 1 ) * nzi
          z = 1.0D+00
          do ii = 1, 2
            do j = 1, nz2
              js = js + 1
              j2 = js + nz
              hold = x(js) + z * x(j2)
              z = - z
              x(j2) = x(js) + z * x(j2)
              x(js) = hold
              z = - z
            end do
            if ( l .eq. 1 ) then
              go to 3
            end if
            z = - 1.0D+00
          end do
3         continue
        end do
      end do
c
c  Bit reversal section.
c
      nw = 0
      do k = 1, n
c
c  Choose correct index and switch elements if not already switched.
c
        if ( k .lt. nw + 1 ) then
          hold = x(nw+1)
          x(nw+1) = x(k)
          x(k) = hold
        end if
c
c  Bump up series by 1.
c
        do i = 1, m

          ii = i
          if ( nw .lt. two_power(i) ) then
            go to 80
          end if
          mw = nw / two_power(i)
          mw1 = mw / 2
          if ( mw .le. 2 * mw1 ) then
            go to 80
          end if

          nw = nw - two_power(i)

        end do

80      continue

        nw = nw + two_power(ii)

      end do

      return
      end
      subroutine fwt ( n, x, y )

c*********************************************************************72
c
cc FWT performs a fast Walsh transform.
c
c  Discussion:
c
c    This routine performs a fast Walsh transform on an input series X
c    leaving the transformed results in X. 
c    X is dimensioned N, which must be a power of 2.
c    The results of this Walsh transform are in sequency order.
c
c    The output sequence could be normalized by dividing by N.
c
c    Note that the program text in the reference included the line
c      y(jd) = abs ( x(j) - x(j2) )
c    which has been corrected to:
c      y(jd) = x(j) - x(j2)
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    11 March 2011
c
c  Author:
c
c    Ken Beauchamp
c
c  Reference:
c
c    Ken Beauchamp,
c    Walsh functions and their applications,
c    Academic Press, 1975,
c    ISBN: 0-12-084050-2,
c    LC: QA404.5.B33.
c
c  Parameters:
c
c    Input, integer N, the number of items in X.
c    N must be a power of 2.
c
c    Input/output, double precision X(N), the data to be transformed.
c
c    Workspace, real Y(N).
c
      implicit none

      integer n

      integer i
      integer i4_log_2
      integer j
      integer j2
      integer jd
      integer js
      integer l
      integer m
      integer n2
      integer nx
      integer ny
      integer nz
      integer nzi
      integer nzn
      double precision x(n)
      double precision y(n)

      n2 = n / 2
      m = i4_log_2 ( n )

      do l = 1, m

        ny = 0
        nz = 2**(l-1)
        nzi = 2 * nz
        nzn = n / nzi

        do i = 1, nzn

          nx = ny + 1
          ny = ny + nz
          js = ( i - 1 ) * nzi
          jd = js + nzi + 1

          do j = nx, ny
            js = js + 1
            j2 = j + n2
            y(js) = x(j) + x(j2)
            jd = jd - 1
            y(jd) = x(j) - x(j2)
          end do

        end do

        do j = 1, n
          x(j) = y(j)
        end do

      end do

      return
      end
      subroutine haar ( n, x, y )

c*********************************************************************72
c
cc HAAR performs a Haar transform.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    11 March 2011
c
c  Author:
c
c    Ken Beauchamp
c
c  Reference:
c
c    Ken Beauchamp,
c    Walsh functions and their applications,
c    Academic Press, 1975,
c    ISBN: 0-12-084050-2,
c    LC: QA404.5.B33.
c
c  Parameters:
c
c    Input, integer N, the number of items in X.
c    N must be a power of 2.
c
c    Input/output, double precision X(N), the data to be transformed.
c
c    Workspace, real Y(N).
c
      implicit none

      integer n

      integer i
      integer i1
      integer i4_log_2
      integer j
      integer jj
      integer k
      integer l
      integer l2
      integer l3
      double precision x(n)
      double precision y(n)

      k = i4_log_2 ( n )

      do i = 1, k

        l = k + 1 - i
        l2 = 2**( l - 1 )

        do i1 = 1, 2 * l2
          y(i1) = x(i1)
        end do

        do j = 1, l2
           l3 = l2 + j
           jj = 2 * j - 1
           x(j) = y(jj) + y(jj+1)
           x(l3) = y(jj) - y(jj+1)
        end do

      end do

      return
      end
      subroutine haarin ( n, x, y )

c*********************************************************************72
c
cc HAARIN inverts a Haar transform.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    11 March 2011
c
c  Author:
c
c    Ken Beauchamp
c
c  Reference:
c
c    Ken Beauchamp,
c    Walsh functions and their applications,
c    Academic Press, 1975,
c    ISBN: 0-12-084050-2,
c    LC: QA404.5.B33.
c
c  Parameters:
c
c    Input, integer N, the number of items in X.
c    N must be a power of 2.
c
c    Input/output, double precision X(N), the data to be transformed.
c
c    Workspace, real Y(N).
c
      implicit none

      integer n

      integer i
      integer i1
      integer i4_log_2
      integer j
      integer jj
      integer jj1
      integer k
      integer l
      integer lj
      double precision x(n)
      double precision y(n)

      k = i4_log_2 ( n )

      do i = 1, k

        l = 2**( i - 1 )

        do i1 = 1, 2 * l
          y(i1) = x(i1)
        end do

        do j = 1, l
          lj = l + j
          jj = 2 * j
          jj1 = jj - 1
          x(jj) = y(j) - y(lj)
          x(jj1) = y(j) + y(lj)
        end do

      end do

      return
      end
      subroutine hnorm ( n, x )

c*********************************************************************72
c
cc HNORM computes the normalization factors for a forward or inverse Haar transform.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    11 March 2011
c
c  Author:
c
c    Ken Beauchamp
c
c  Reference:
c
c    Ken Beauchamp,
c    Walsh functions and their applications,
c    Academic Press, 1975,
c    ISBN: 0-12-084050-2,
c    LC: QA404.5.B33.
c
c  Parameters:
c
c    Input, integer N, the number of items in X.
c    N must be a power of 2.
c
c    Input/output, double precision X(N), the data to be transformed.
c
      implicit none

      integer n

      integer i
      integer i4_log_2
      integer ii
      integer j
      integer jmax
      integer jmin
      integer k
      double precision wlk
      double precision x(n)

      k = i4_log_2 ( n )

      x(1) = x(1) / 2.0D+00**k

      if ( 1 .le. k ) then
        x(2) = x(2) / 2.0D+00**k
      end if

      do ii = 2, k

        i = ii - 1
        wlk = 1.0D+00 / 2.0D+00**( k - i )
        jmin = 2**i + 1
        jmax = 2**ii

        do j = jmin, jmax
          x(j) = x(j) * wlk
        end do

      end do

      return
      end
      function i4_log_2 ( i )

c*********************************************************************72
c
cc I4_LOG_2 returns the integer part of the logarithm base 2 of an I4.
c
c  Discussion:
c
c    For positive I4_LOG_2(I), it should be true that
c      2**I4_LOG_2(X) .le. |I| < 2**(I4_LOG_2(I)+1).
c    The special case of I4_LOG_2(0) returns -HUGE().
c
c    An I4 is an integer value.
c
c  Example:
c
c     I  I4_LOG_2
c
c     0  -1
c     1,  0
c     2,  1
c     3,  1
c     4,  2
c     5,  2
c     6,  2
c     7,  2
c     8,  3
c     9,  3
c    10,  3
c   127,  6
c   128,  7
c   129,  7
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    28 October 2007
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer I, the number whose logarithm base 2
c    is desired.
c
c    Output, integer I4_LOG_2, the integer part of the
c    logarithm base 2 of the absolute value of I.
c
      implicit none

      integer i
      integer i_abs
      integer i4_log_2
      integer i4_huge
      parameter ( i4_huge = 2147483647 )

      if ( i .eq. 0 ) then

        i4_log_2 = - i4_huge

      else

        i4_log_2 = 0

        i_abs = abs ( i )

10      continue

        if ( 2 .le. i_abs ) then
          i_abs = i_abs / 2
          i4_log_2 = i4_log_2 + 1
          go to 10
        end if

      end if

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
      function i4_wrap ( ival, ilo, ihi )

c*********************************************************************72
c
cc I4_WRAP forces an I4 to lie between given limits by wrapping.
c
c  Example:
c
c    ILO = 4, IHI = 8
c
c    I  Value
c
c    -2     8
c    -1     4
c     0     5
c     1     6
c     2     7
c     3     8
c     4     4
c     5     5
c     6     6
c     7     7
c     8     8
c     9     4
c    10     5
c    11     6
c    12     7
c    13     8
c    14     4
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
c    Input, integer IVAL, an integer value.
c
c    Input, integer ILO, IHI, the desired bounds for the integer value.
c
c    Output, integer I4_WRAP, a "wrapped" version of IVAL.
c
      implicit none

      integer i4_modp
      integer i4_wrap
      integer ihi
      integer ilo
      integer ival
      integer jhi
      integer jlo
      integer value
      integer wide

      jlo = min ( ilo, ihi )
      jhi = max ( ilo, ihi )

      wide = jhi - jlo + 1

      if ( wide .eq. 1 ) then
        value = jlo
      else
        value = jlo + i4_modp ( ival - jlo, wide )
      end if

      i4_wrap = value

      return
      end
      subroutine r8vec_copy ( n, a1, a2 )

c*********************************************************************72
c
cc R8VEC_COPY copies an R8VEC.
c
c  Discussion:
c
c    An R8VEC is a vector of R8 values.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    13 August 2008
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the length of the vectors.
c
c    Input, double precision A1(N), the vector to be copied.
c
c    Output, double precision A2(N), a copy of A1.
c
      implicit none

      integer n

      double precision a1(n)
      double precision a2(n)
      integer i

      do i = 1, n
        a2(i) = a1(i)
      end do

      return
      end
      subroutine r8vec_shift_circular ( shift, n, x )

c*********************************************************************72
c
cc R8VEC_SHIFT_CIRCULAR performs a circular shift on an R8VEC.
c
c  Discussion:
c
c    An R8VEC is a vector of R8 values.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    11 March 2011
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer SHIFT, the amount by which each entry is to
c    be shifted.
c
c    Input, integer N, the length of the vector.
c
c    Input/output, double precision A1(N), the vector to be shifted.
c
      implicit none

      integer n

      integer i
      integer i4_wrap
      integer j
      integer shift
      double precision x(n)
      double precision y(n)

      do i = 1, n
        y(i) = x(i)
      end do

      do i = 1, n
        j = i4_wrap ( i - shift, 1, n )
        x(i) = y(j)
      end do

      return
      end
      subroutine r8vec_uniform_01 ( n, seed, r )

c*********************************************************************72
c
cc R8VEC_UNIFORM_01 returns a unit pseudorandom R8VEC.
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
c    17 July 2006
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
c    Input, integer N, the number of entries in the vector.
c
c    Input/output, integer SEED, the "seed" value, which should NOT be 0.
c    On output, SEED has been updated.
c
c    Output, double precision R(N), the vector of pseudorandom values.
c
      implicit none

      integer n

      integer i
      integer k
      integer seed
      double precision r(n)

      do i = 1, n

        k = seed / 127773

        seed = 16807 * ( seed - k * 127773 ) - k * 2836

        if ( seed .lt. 0 ) then
          seed = seed + 2147483647
        end if

        r(i) = dble ( seed ) * 4.656612875D-10

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
      subroutine walsh ( n, x, y )

c*********************************************************************72
c
cc WALSH performs a fast Walsh transform.
c
c  Discussion:
c
c    This routine performs a fast Wash transform on an input series X
c    leaving the transformed results in X.  The array Y is used for working space.
c    X and Y are dimensioned N, which must be a power of 2.
c    The results of this Walsh transform are in sequency order.
c
c    The output sequence could be normalized by dividing by N.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    11 March 2011
c
c  Author:
c
c    Ken Beauchamp
c
c  Reference:
c
c    Ken Beauchamp,
c    Walsh functions and their applications,
c    Academic Press, 1975,
c    ISBN: 0-12-084050-2,
c    LC: QA404.5.B33.
c
c  Parameters:
c
c    Input, integer N, the number of items in X.
c    N must be a power of 2.
c
c    Input/output, double precision X(N), the data to be transformed.
c
c    Workspace, real Y(N).
c
      implicit none

      integer n

      double precision a
      integer i
      integer i1
      integer i4_log_2
      integer is
      integer j
      integer j1
      integer l
      integer m
      integer n1
      integer n2
      double precision w
      double precision x(n)
      double precision y(n/2)
      double precision z

      n2 = n / 2
      m = i4_log_2 ( n )
      z = -1.0D+00

      do j = 1, m

        n1 = 2**( m - j + 1 )
        j1 = 2**( j - 1 )

        do l = 1, j1

          is = ( l - 1 ) * n1 + 1
          i1 = 0
          w = z

          do i = is, is + n1 - 1, 2
            a = x(i)
            x(is+i1) = a + x(i+1)
            i1 = i1 + 1
            y(i1) = ( x(i+1) - a ) * w
            w = w * z
          end do

          do i = 1, n1 / 2
            x(n1/2+is-1+i) = y(i)
          end do

        end do

      end do

      return
      end
