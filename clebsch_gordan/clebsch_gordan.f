      function cof3j ( j1, j2, j3, m1, m2, m3 )

c*********************************************************************72
c
cc COF3J evaluates the Wigner 3J coefficients.
c
c  Discussion:
c
c    The J and M parameters should be integer or half integer.
c
c    The M values should sum to 0.
c
c    The J values should be nonnegative, and should satisfy the
c    triangle inequalities; that is, 
c
c      J1 <= J2 + J3
c      J2 <= J1 + J3
c      J3 <= J1 + J2
c
c    In Mathematica, the function can be evaluated by:
c
c      ThreeJSymbol[{j1,m1},{j2,m2},{j3,m3}]
c
c  Modified:
c
c    08 February 2007
c
c  Reference:
c
c    Albert Messiah,
c    Quantum Mechanics,
c    Volume II,
c    North Holland, 1963,
c    ISBN13: 978-0486409245,
c    LC: QC174.1.M413.
c
c  Parameters:
c
c    Input, double precision J1, J2, J3, M1, M2, M3, the parameters.
c
c    Output, double precision COF3J, the value of the Wigner 3J coefficient.
c
      implicit none

      double precision cof3j
      double precision h
      double precision j1
      integer j1_copy
      double precision j2
      integer j2_copy
      double precision j3
      integer j3_copy
      double precision m1
      integer m1_copy
      double precision m2
      integer m2_copy
      double precision m3
      double precision vcc

      h = j1 - j2 - m3

      j1_copy = nint ( 2.0D+00 * j1 )
      j2_copy = nint ( 2.0D+00 * j2 )
      j3_copy = nint ( 2.0D+00 * j3 )

      m1_copy = nint ( 2.0D+00 * m1 )
      m2_copy = nint ( 2.0D+00 * m2 )

      cof3j = vcc ( j1_copy, j2_copy, j3_copy, m1_copy, m2_copy ) 
     &  / sqrt ( dble ( 2 * j3 + 1 ) )

      if ( mod ( nint ( abs ( h ) ), 2 ) .ne. 0 ) then
        cof3j = - cof3j
      end if

      return
      end
      function cof6j ( j1, j2, j3, j4, j5, j6 )

c*********************************************************************72
c
cc COF6J evaluates the Wigner 6J coefficients.
c
c  Discussion:
c
c    In Mathematica, the function can be evaluated by:
c
c      SixJSymbol[{j1,j2,j3},{j4,j5,j6}]
c
c    The triangle inequality should be satisfied by certain triples 
c    of the J values, namely
c
c      (J1,J2,J3),
c      (J1,J5,J6),
c      (J2,J4,J6) and
c      (J3,J4,J5).
c
c    For (J1,J2,J3), this implies:
c
c      J1 <= J2 + J3
c      J2 <= J1 + J3
c      J3 <= J1 + J2
c
c  Modified:
c
c    09 February 2007
c
c  Reference:
c
c    Albert Messiah,
c    Quantum Mechanics,
c    Volume II,
c    North Holland, 1963,
c    ISBN13: 978-0486409245,
c    LC: QC174.1.M413.
c
c  Parameters:
c
c    Input, double precision J1, J2, J3, J4, J5, J6, the parameters.
c
c    Output, double precision COF6J, the value of the Wigner 6J coefficient.
c
      implicit none

      double precision cof6j
      integer ih
      double precision j1
      integer j12
      double precision j2
      integer j22
      double precision j3
      integer j32
      double precision j4
      integer j42
      double precision j5
      integer j52
      double precision j6
      integer j62
      double precision racah

      ih = nint ( j1 + j2 + j4 + j5 )

      j12 = nint ( 2.0D+00 * j1 )
      j22 = nint ( 2.0D+00 * j2 )
      j32 = nint ( 2.0D+00 * j3 )
      j42 = nint ( 2.0D+00 * j4 )
      j52 = nint ( 2.0D+00 * j5 )
      j62 = nint ( 2.0D+00 * j6 )

      cof6j = racah ( j12, j22, j32, j42, j52, j62 )

      if ( mod ( ih, 2 ) .ne. 0 ) then
        cof6j = - cof6j
      end if

      return
      end
      function cof9j ( j1, j2, j3, j4, j5, j6, j7, j8, j9 )

c*********************************************************************72
c
cc COF9J evaluates the Wigner 9J coefficients.
c
c  Discussion:
c
c    The triangle inequality should be satisfied by certain triples 
c    of the J values, namely
c
c      (J1,J2,J3),
c      (J4,J5,J6),
c      (J7,J8,J9),
c      (J1,J4,J7),
c      (J2,J5,J8) and
c      (J3,J6,J9).
c
c    For (J1,J2,J3), this implies:
c
c      J1 <= J2 + J3
c      J2 <= J1 + J3
c      J3 <= J1 + J2
c
c  Modified:
c
c    09 February 2007
c
c  Reference:
c
c    Albert Messiah,
c    Quantum Mechanics,
c    Volume II,
c    North Holland, 1963,
c    ISBN13: 978-0486409245,
c    LC: QC174.1.M413.
c
c  Parameters:
c
c    Input, double precision J1, J2, J3, J4, J5, J6, J7, J8, J9,
c    the parameters.
c
c    Output, double precision COF9J, the value of the Wigner 9J coefficient.
c
      implicit none

      double precision cof9j
      double precision j1
      integer j12
      double precision j2
      integer j22
      double precision j3
      integer j32
      double precision j4
      integer j42
      double precision j5
      integer j52
      double precision j6
      integer j62
      double precision j7
      integer j72
      double precision j8
      integer j82
      double precision j9
      integer j92
      double precision winej

      j12 = nint ( 2.0D+00 * j1 )
      j22 = nint ( 2.0D+00 * j2 )
      j32 = nint ( 2.0D+00 * j3 )
      j42 = nint ( 2.0D+00 * j4 )
      j52 = nint ( 2.0D+00 * j5 )
      j62 = nint ( 2.0D+00 * j6 )
      j72 = nint ( 2.0D+00 * j7 )
      j82 = nint ( 2.0D+00 * j8 )
      j92 = nint ( 2.0D+00 * j9 )

      cof9j = winej ( j12, j22, j32, j42, j52, j62, j72, j82, j92 )

      return
      end
      subroutine facalc ( setup )

c*********************************************************************72
c
cc FACALC calculates log factorials and saves them in common.
c
c  Modified:
c
c    04 February 2007
c
c  Parameters:
c
c    Output, logical SETUP, is returned as FALSE, indicating that
c    the data does not need to be set up again.
c
      double precision fl(322)
      integer n
      logical setup

      save / factrl /

      common / factrl / fl

      fl(1) = 0.0D+00
      fl(2) = 0.0D+00
      do n = 3, 322
        fl(n) = fl(n-1) + log ( dble ( n - 1 ) )
      end do

      setup = .false.

      return
      end
      subroutine i4_swap ( i, j )

c*********************************************************************72
c
cc I4_SWAP switches two I4's.
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
      subroutine nine_j_values ( n_data, j1, j2, j3, j4, j5, j6, j7, 
     &  j8, j9, fx )

c*********************************************************************72
c
cc NINE_J_VALUES returns some values of the Wigner 9J function.
c
c  Modified:
c
c    09 February 2007
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
c  Parameters:
c
c    Input/output, integer N_DATA.  The user sets N_DATA to 0 before the
c    first call.  On each call, the routine increments N_DATA by 1, and
c    returns the corresponding data; when there is no more data, the
c    output value of N_DATA will be 0 again.
c
c    Output, double precision J1, J2, J3, J4, J5, J6, J7, J8, J9, 
c    the arguments of the function.
c
c    Output, double precision FX, the value of the function.
c
      implicit none

      integer n_max
      parameter ( n_max = 9 )

      double precision fx
      double precision fx_vec(n_max)
      integer n_data
      double precision j1
      double precision j1_vec(n_max)
      double precision j2
      double precision j2_vec(n_max)
      double precision j3
      double precision j3_vec(n_max)
      double precision j4
      double precision j4_vec(n_max)
      double precision j5
      double precision j5_vec(n_max)
      double precision j6
      double precision j6_vec(n_max)
      double precision j7
      double precision j7_vec(n_max)
      double precision j8
      double precision j8_vec(n_max)
      double precision j9
      double precision j9_vec(n_max)

      save fx_vec
      save j1_vec
      save j2_vec
      save j3_vec
      save j4_vec
      save j5_vec
      save j6_vec
      save j7_vec
      save j8_vec
      save j9_vec

      data fx_vec / 
     &   0.0004270039294528318D+00, 
     &  -0.001228915451058514D+00, 
     &  -0.0001944260688400887D+00, 
     &   0.003338419923885592D+00, 
     &  -0.0007958936865080434D+00, 
     &  -0.004338208690251972D+00, 
     &   0.05379143536399187D+00, 
     &   0.006211299937499411D+00, 
     &   0.03042903097250921D+00 /
      data j1_vec / 
     &  1.0D+00, 
     &  1.5D+00, 
     &  2.0D+00, 
     &  1.0D+00, 
     &  1.5D+00, 
     &  2.0D+00, 
     &  0.5D+00, 
     &  1.0D+00, 
     &  1.5D+00  /
      data j2_vec / 
     &  8.0D+00, 
     &  8.0D+00, 
     &  8.0D+00, 
     &  3.0D+00, 
     &  3.0D+00, 
     &  3.0D+00, 
     &  0.5D+00, 
     &  0.5D+00, 
     &  0.5D+00  /
      data j3_vec / 
     &  7.0D+00, 
     &  7.0D+00, 
     &  7.0D+00, 
     &  2.0D+00, 
     &  2.0D+00, 
     &  2.0D+00, 
     &  1.0D+00, 
     &  1.0D+00, 
     &  1.0D+00 /
      data j4_vec / 
     &  6.5D+00, 
     &  6.5D+00, 
     &  6.5D+00, 
     &  4.0D+00, 
     &  4.0D+00, 
     &  4.0D+00, 
     &  2.0D+00, 
     &  2.0D+00, 
     &  2.0D+00 /
      data j5_vec / 
     &  7.5D+00, 
     &  7.5D+00, 
     &  7.5D+00, 
     &  1.5D+00, 
     &  1.5D+00, 
     &  1.5D+00, 
     &  1.0D+00, 
     &  1.0D+00, 
     &  1.0D+00 /
      data j6_vec / 
     &  7.5D+00, 
     &  7.5D+00, 
     &  7.5D+00, 
     &  3.0D+00, 
     &  3.0D+00, 
     &  3.0D+00, 
     &  1.5D+00, 
     &  1.5D+00, 
     &  1.5D+00 /
      data j7_vec / 
     &  6.0D+00, 
     &  6.0D+00, 
     &  6.0D+00, 
     &  3.5D+00, 
     &  3.5D+00, 
     &  3.5D+00, 
     &  1.5D+00, 
     &  1.5D+00, 
     &  1.5D+00 /
      data j8_vec / 
     &  10.0D+00, 
     &  10.0D+00, 
     &  10.0D+00, 
     &   2.0D+00, 
     &   2.0D+00, 
     &   2.0D+00, 
     &   0.5D+00, 
     &   0.5D+00, 
     &   0.5D+00 /
      data j9_vec / 
     &  6.0D+00, 
     &  6.0D+00, 
     &  6.0D+00, 
     &  2.0D+00, 
     &  2.0D+00, 
     &  2.0D+00, 
     &  1.5D+00, 
     &  1.5D+00, 
     &  1.5D+00 /

      if ( n_data .lt. 0 ) then
        n_data = 0
      end if

      n_data = n_data + 1

      if ( n_max .lt. n_data ) then
        n_data = 0
        j1 = 0.0D+00
        j2 = 0.0D+00
        j3 = 0.0D+00
        j4 = 0.0D+00
        j5 = 0.0D+00
        j6 = 0.0D+00
        j7 = 0.0D+00
        j8 = 0.0D+00
        j9 = 0.0D+00
        fx = 0.0D+00
      else
        j1 = j1_vec(n_data)
        j2 = j2_vec(n_data)
        j3 = j3_vec(n_data)
        j4 = j4_vec(n_data)
        j5 = j5_vec(n_data)
        j6 = j6_vec(n_data)
        j7 = j7_vec(n_data)
        j8 = j8_vec(n_data)
        j9 = j9_vec(n_data)
        fx = fx_vec(n_data)
      end if

      return
      end
      function phasef ( n )

c*********************************************************************72
c
cc PHASEF returns 1 for even arguments and -1 for odd ones.
c
c  Modified:
c
c    04 February 2007
c
c  Parameters:
c
c    Input, integer N, an integer to be tested.
c
c    Output, double precision PHASEF, is 1 if N is even, and -1 if N is odd.
c
      implicit none

      integer n
      double precision phasef

      if ( mod ( n, 2 ) .eq. 0 ) then
        phasef = 1.0D+00
      else
        phasef = -1.0D+00
      end if

      return
      end
      function racah ( jd1, jd2, jd3, ld1, ld2, ld3 )

c*********************************************************************72
c
cc RACAH ???
c
c  Modified:
c
c    04 February 2007
c
c  Reference:
c
c    Albert Messiah,
c    Quantum Mechanics,
c    Volume II,
c    North Holland, 1963,
c    ISBN13: 978-0486409245,
c    LC: QC174.1.M413.
c
c  Parameters:
c
c    Input, integer JD1, JD2, JD3, LD1, LD2, LD3, the parameters.
c
c    Output, double precision RACAH, ???
c
      implicit none

      double precision f6j
      integer i
      integer i1
      integer j1
      integer j2
      integer j3
      integer jd1
      integer jd2
      integer jd3
      integer l1
      integer l2
      integer l3
      integer ld1
      integer ld2
      integer ld3
      integer med(12)
      integer n
      double precision phasef
      double precision racah
      double precision s6j

      racah = 0.0D+00

      j1 = jd1
      j2 = jd2
      j3 = jd3

      l1 = ld1
      l2 = ld2
      l3 = ld3
c
c  Angular momentum coupling tests.
c
      i = - j1 + j2 + j3
      i1  = i / 2

      if ( i .ne. 2 * i1 ) then
        return
      end if

      med(1) = i1
      i = + j1 - j2 + j3
      i1 = i / 2

      if ( i .ne. 2 * i1 ) then
        return
      end if

      med(2) = i1
      i = + j1 + j2 - j3
      i1 = i / 2

      if ( i .ne. 2 * i1 ) then
        return
      end if

      med(3) = i1
      i = - j1 + l2 + l3
      i1 = i / 2

      if ( i .ne. 2 * i1 ) then
        return
      end if

      med(4) = i1
      i = + j1 - l2 + l3
      i1 = i / 2

      if ( i .ne. 2 * i1 ) then
        return
      end if

      med(5) = i1
      i = + j1 + l2 - l3
      i1 = i / 2

      if ( i .ne. 2 * i1 ) then
        return
      end if

      med(6) = i1
      i = - l1 + j2 + l3
      i1 = i / 2

      if ( i .ne. 2 * i1 ) then
        return
      end if

      med(7) = i1
      i = + l1 - j2 + l3
      i1 = i / 2

      if ( i .ne. 2 * i1 ) then
        return
      end if

      med(8) = i1
      i = + l1 + j2 - l3
      i1 = i / 2

      if ( i .ne. 2 * i1 ) then
        return
      end if

      med(9) = i1
      i = - l1 + l2 + j3
      i1 = i / 2

      if ( i .ne. 2 * i1 ) then
        return
      end if

      med(10) = i1
      i = + l1 - l2 + j3
      i1 = i / 2

      if ( i .ne. 2 * i1 ) then
        return
      end if

      med(11) = i1
      i = + l1 + l2 - j3
      i1 = i / 2

      if ( i .ne. 2 * i1 ) then
        return
      end if

      med(12) = i1

      do n = 1, 12
        if ( med(n) .lt. 0 ) then
          return
        end if
      end do

      f6j = s6j ( j1, j2, j3, l1, l2, l3 )
      
      racah = f6j * phasef ( ( j1 + j2 + l1 + l2 ) / 2 )

      return
      end
      function s6j ( jd1, jd2, jd3, ld1, ld2, ld3 )

c*********************************************************************72
c
cc S6J ???
c
c  Modified:
c
c    04 February 2007
c
c  Parameters:
c
c    Input, integer JD1, JD2, JD3, LD1, LD2, LD3, the parameters.
c
c    Output, double precision S6J, ???
c
      implicit none

      double precision delog
      double precision fl(322)
      integer j1
      integer j2
      integer j3
      integer jd1
      integer jd2
      integer jd3
      integer k
      integer kmax
      integer l1
      integer l2
      integer l3
      integer ld1
      integer ld2
      integer ld3
      integer ma(4)
      integer max_a
      integer mb(3)
      integer med(12)
      integer min_b
      integer min1
      integer min2
      integer min3
      integer min4
      integer min5
      integer min6
      integer min7
      integer min8
      integer minp1
      integer n
      integer num1
      integer num2
      integer num3
      integer num4
      double precision p
      double precision plog
      double precision q
      double precision s
      double precision s6j
      logical setup
      double precision uk
      double precision ulog

      save / factrl /
      save setup

      common / factrl / fl

      data setup / .true. /

      j1 = jd1
      j2 = jd2
      j3 = jd3

      l1 = ld1
      l2 = ld2
      l3 = ld3

      if ( setup ) then
        call facalc ( setup )
      end if

      med(1)  = ( - j1 + j2 + j3 ) / 2
      med(2)  = ( + j1 - j2 + j3 ) / 2
      med(3)  = ( + j1 + j2 - j3 ) / 2
      med(4)  = ( - j1 + l2 + l3 ) / 2
      med(5)  = ( + j1 - l2 + l3 ) / 2
      med(6)  = ( + j1 + l2 - l3 ) / 2
      med(7)  = ( - l1 + j2 + l3 ) / 2
      med(8)  = ( + l1 - j2 + l3 ) / 2
      med(9)  = ( + l1 + j2 - l3 ) / 2
      med(10) = ( - l1 + l2 + j3 ) / 2
      med(11) = ( + l1 - l2 + j3 ) / 2
      med(12) = ( + l1 + l2 - j3 ) / 2

      ma(1) = med(1)  + med(2)  + med(3)
      ma(2) = med(4)  + med(5)  + med(6)
      ma(3) = med(7)  + med(8)  + med(9)
      ma(4) = med(10) + med(11) + med(12)

      mb(1) = ma(1) + med(12)
      mb(2) = ma(1) + med(4)
      mb(3) = ma(1) + med(8)
c
c  Determine maximum of (j1+j2+j3),(j1+l2+l3),(l1+j2+l3),(l1+l2+j3)
c
      max_a = ma(1)
      do n = 2, 4
        if ( max_a .lt. ma(n) ) then
          max_a = ma(n)
        end if
      end do
c
c  Determine minimum of (j1+j2+l1+l2),(j2+j3+l2+l3),(j3+j1+l3+l1)
c
      min_b = mb(1)
      do n = 2, 3
        if ( mb(n) .lt. min_b ) then 
          min_b = mb(n)
        end if
      end do

      kmax = min_b - max_a
      minp1 = min_b + 1
      min1 = min_b + 1 - ma(1)
      min2 = min_b + 1 - ma(2)
      min3 = min_b + 1 - ma(3)
      min4 = min_b + 1 - ma(4)
      min5 = min_b + 2
      min6 = mb(1) - min_b
      min7 = mb(2) - min_b
      min8 = mb(3) - min_b
c
c  Sum the series.
c
      uk = 1.0D-15
      s = 1.0D-15
      do k = 1, kmax

        uk = - uk 
     &    * dble ( ( min1 - k ) * ( min2 - k ) 
     &           * ( min3 - k ) * ( min4 - k ) ) 
     &    / dble ( ( min5 - k ) * ( min6 + k ) 
     &           * ( min7 + k ) * ( min8 + k ) )

        if ( abs ( uk ) .le. 1.0D-25 ) then
          go to 10
        end if

        s = s + uk

      end do

10    continue

      s = s * 1.0D+15
c
c  Calculate delta functions.
c
      delog = 0.0D+00
      do n = 1, 12
        delog = delog + fl ( med(n)+1 )
      end do

      num1 = ma(1) + 2
      num2 = ma(2) + 2
      num3 = ma(3) + 2
      num4 = ma(4) + 2

      delog = delog - fl(num1) - fl(num2) - fl(num3) - fl(num4)
      delog = 0.5D+00 * delog
      ulog = fl(min5) - fl(min1) - fl(min2) - fl(min3) - fl(min4)
     &  - fl(min6+1) - fl(min7+1) - fl(min8+1)
      plog = delog + ulog

      if ( plog .lt. -64.0D+00 ) then
        q = plog + 64.0D+00
        q = exp ( q )
        s6j = q * s
        if ( abs ( s6j ) .le. 1.0D+00 ) then
          s6j = 0.0D+00
          return
        end if
        s6j = s6j * exp ( -64.0D+00 )
      else
        p = exp ( plog )
        s6j = p * s
      end if

      if ( mod ( min_b, 2 ) .eq. 1 ) then
        s6j = - s6j
      end if

      return
      end
      subroutine six_j_values ( n_data, j1, j2, j3, j4, j5, j6, fx )

c*********************************************************************72
c
cc SIX_J_VALUES returns some values of the Wigner 6J function.
c
c  Discussion:
c
c    In Mathematica, the function can be evaluated by:
c
c      SixJSymbol[{j1,j2,j3},{j4,j5,j6}]
c
c  Modified:
c
c    07 February 2007
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
c    Stephen Wolfram,
c    The Mathematica Book,
c    Fourth Edition,
c    Cambridge University Press, 1999,
c    ISBN: 0-521-64314-7,
c    LC: QA76.95.W65.
c
c  Parameters:
c
c    Input/output, integer N_DATA.  The user sets N_DATA to 0 before the
c    first call.  On each call, the routine increments N_DATA by 1, and
c    returns the corresponding data; when there is no more data, the
c    output value of N_DATA will be 0 again.
c
c    Output, double precision J1, J2, J3, J4, J5, J6, the arguments
c    of the function.
c
c    Output, double precision FX, the value of the function.
c
      implicit none

      integer n_max
      parameter ( n_max = 15 )

      double precision fx
      double precision fx_vec(n_max)
      double precision j1
      double precision j1_vec(n_max)
      double precision j2
      double precision j2_vec(n_max)
      double precision j3
      double precision j3_vec(n_max)
      double precision j4
      double precision j4_vec(n_max)
      double precision j5
      double precision j5_vec(n_max)
      double precision j6
      double precision j6_vec(n_max)
      integer n_data

      save fx_vec
      save j1_vec
      save j2_vec
      save j3_vec
      save j4_vec
      save j5_vec
      save j6_vec

      data fx_vec /
     &   0.03490905138373300D+00,  
     &  -0.03743025039659792D+00,  
     &   0.01890866390959560D+00,  
     &   0.007342448254928643D+00, 
     &  -0.02358935185081794D+00,  
     &   0.01913476955215437D+00,  
     &   0.001288017397724172D+00, 
     &  -0.01930018366290527D+00,   
     &   0.01677305949382889D+00,  
     &   0.005501147274850949D+00, 
     &  -0.02135439790896831D+00,  
     &   0.003460364451435387D+00, 
     &   0.02520950054795585D+00,  
     &   0.01483990561221713D+00,  
     &   0.002708577680633186D+00 /
      data j1_vec /
     &  1.0D+00, 
     &  2.0D+00, 
     &  3.0D+00, 
     &  4.0D+00, 
     &  5.0D+00, 
     &  6.0D+00, 
     &  7.0D+00, 
     &  8.0D+00, 
     &  9.0D+00, 
     & 10.0D+00, 
     & 11.0D+00, 
     & 12.0D+00, 
     & 13.0D+00, 
     & 14.0D+00, 
     & 15.0D+00 /
      data j2_vec / 
     &  8.0D+00, 
     &  8.0D+00, 
     &  8.0D+00, 
     &  8.0D+00, 
     &  8.0D+00, 
     &  8.0D+00, 
     &  8.0D+00, 
     &  8.0D+00, 
     &  8.0D+00, 
     &  8.0D+00, 
     &  8.0D+00, 
     &  8.0D+00,
     &  8.0D+00, 
     &  8.0D+00, 
     &  8.0D+00 /
      data j3_vec / 
     &  7.0D+00, 
     &  7.0D+00, 
     &  7.0D+00, 
     &  7.0D+00, 
     &  7.0D+00, 
     &  7.0D+00, 
     &  7.0D+00, 
     &  7.0D+00, 
     &  7.0D+00, 
     &  7.0D+00, 
     &  7.0D+00, 
     &  7.0D+00, 
     &  7.0D+00, 
     &  7.0D+00, 
     &  7.0D+00 /
      data j4_vec / 
     &  6.5D+00, 
     &  6.5D+00, 
     &  6.5D+00, 
     &  6.5D+00, 
     &  6.5D+00, 
     &  6.5D+00, 
     &  6.5D+00, 
     &  6.5D+00, 
     &  6.5D+00, 
     &  6.5D+00, 
     &  6.5D+00, 
     &  6.5D+00, 
     &  6.5D+00, 
     &  6.5D+00, 
     &  6.5D+00 /
      data j5_vec / 
     &  7.5D+00, 
     &  7.5D+00, 
     &  7.5D+00, 
     &  7.5D+00, 
     &  7.5D+00, 
     &  7.5D+00, 
     &  7.5D+00, 
     &  7.5D+00, 
     &  7.5D+00, 
     &  7.5D+00, 
     &  7.5D+00, 
     &  7.5D+00, 
     &  7.5D+00, 
     &  7.5D+00, 
     &  7.5D+00 /
      data j6_vec / 
     &  7.5D+00, 
     &  7.5D+00, 
     &  7.5D+00, 
     &  7.5D+00, 
     &  7.5D+00, 
     &  7.5D+00, 
     &  7.5D+00, 
     &  7.5D+00, 
     &  7.5D+00, 
     &  7.5D+00, 
     &  7.5D+00, 
     &  7.5D+00, 
     &  7.5D+00, 
     &  7.5D+00, 
     &  7.5D+00 /

      if ( n_data .lt. 0 ) then
        n_data = 0
      end if

      n_data = n_data + 1

      if ( n_max .lt. n_data ) then
        n_data = 0
        j1 = 0.0D+00
        j2 = 0.0D+00
        j3 = 0.0D+00
        j4 = 0.0D+00
        j5 = 0.0D+00
        j6 = 0.0D+00
        fx = 0.0D+00
      else
        j1 = j1_vec(n_data)
        j2 = j2_vec(n_data)
        j3 = j3_vec(n_data)
        j4 = j4_vec(n_data)
        j5 = j5_vec(n_data)
        j6 = j6_vec(n_data)
        fx = fx_vec(n_data)
      end if

      return
      end
      subroutine three_j_values ( n_data, j1, j2, j3, m1, m2, m3, fx )

c*********************************************************************72
c
cc THREE_J_VALUES returns some values of the Wigner 3J function.
c
c  Discussion:
c
c    In Mathematica, the function can be evaluated by:
c
c      ThreeJSymbol[{j1,m1},{j2,m2},{j3,m3}]
c
c  Modified:
c
c    07 February 2007
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
c    Stephen Wolfram,
c    The Mathematica Book,
c    Fourth Edition,
c    Cambridge University Press, 1999,
c    ISBN: 0-521-64314-7,
c    LC: QA76.95.W65.
c
c  Parameters:
c
c    Input/output, integer N_DATA.  The user sets N_DATA to 0 before the
c    first call.  On each call, the routine increments N_DATA by 1, and
c    returns the corresponding data; when there is no more data, the
c    output value of N_DATA will be 0 again.
c
c    Output, double precision J1, J2, J3, M1, M2, M3, the arguments
c    of the function.
c
c    Output, double precision FX, the value of the function.
c
      implicit none

      integer n_max
      parameter ( n_max = 8 )

      double precision fx
      double precision fx_vec(n_max)
      integer n_data
      double precision j1
      double precision j1_vec(n_max)
      double precision j2
      double precision j2_vec(n_max)
      double precision j3
      double precision j3_vec(n_max)
      double precision m1
      double precision m1_vec(n_max)
      double precision m2
      double precision m2_vec(n_max)
      double precision m3
      double precision m3_vec(n_max)

      save fx_vec
      save j1_vec
      save j2_vec
      save j3_vec
      save m1_vec
      save m2_vec
      save m3_vec

      data fx_vec /
     &   0.2788866755113585D+00,
     &  -0.09534625892455923D+00,
     &  -0.06741998624632421D+00,
     &   0.1533110351679666D+00,
     &  -0.1564465546936860D+00,
     &   0.1099450412156551D+00,
     &  -0.05536235693131719D+00,
     &   0.01799835451137786D+00 /
      data j1_vec /
     &  1.0D+00,
     &  2.0D+00,
     &  3.0D+00,
     &  4.0D+00,
     &  5.0D+00,
     &  6.0D+00,
     &  7.0D+00,
     &  8.0D+00 /
      data j2_vec /
     &  4.5D+00,
     &  4.5D+00,
     &  4.5D+00,
     &  4.5D+00,
     &  4.5D+00,
     &  4.5D+00,
     &  4.5D+00,
     &  4.5D+00 /
      data j3_vec /
     &  3.5D+00,
     &  3.5D+00,
     &  3.5D+00,
     &  3.5D+00,
     &  3.5D+00,
     &  3.5D+00,
     &  3.5D+00,
     &  3.5D+00 /
      data m1_vec /
     &  1.0D+00,
     &  1.0D+00,
     &  1.0D+00,
     &  1.0D+00,
     &  1.0D+00,
     &  1.0D+00,
     &  1.0D+00,
     &  1.0D+00 /
      data m2_vec /
     &  -3.5D+00,
     &  -3.5D+00,
     &  -3.5D+00,
     &  -3.5D+00,
     &  -3.5D+00,
     &  -3.5D+00,
     &  -3.5D+00,
     &  -3.5D+00 /
      data m3_vec /
     &  2.5D+00,
     &  2.5D+00,
     &  2.5D+00,
     &  2.5D+00,
     &  2.5D+00,
     &  2.5D+00,
     &  2.5D+00,
     &  2.5D+00 /

      if ( n_data .lt. 0 ) then
        n_data = 0
      end if

      n_data = n_data + 1

      if ( n_max .lt. n_data ) then
        n_data = 0
        j1 = 0.0D+00
        j2 = 0.0D+00
        j3 = 0.0D+00
        m1 = 0.0D+00
        m2 = 0.0D+00
        m3 = 0.0D+00
        fx = 0.0D+00
      else
        j1 = j1_vec(n_data)
        j2 = j2_vec(n_data)
        j3 = j3_vec(n_data)
        m1 = m1_vec(n_data)
        m2 = m2_vec(n_data)
        m3 = m3_vec(n_data)
        fx = fx_vec(n_data)
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
      function vcc ( jd1, jd2, jd3, md1, md2 )

c*********************************************************************72
c
cc VCC evaluates the Racah V coefficient.
c
c  Discussion:
c
c    The Racah V coefficient is written V(j1,j2,j; m1, m2, m );
c    Generally, it is required that m = - m1 - m2.
c
c  Modified:
c
c    04 February 2007
c
c  Parameters:
c
c    Input, integer JD1, JD2, JD3, MD1, MD2, the parameters.
c
c    Output, double precision VCC, the value of the Racah V coefficient.
c 
      implicit none

      double precision delog
      double precision f3j
      double precision fl(322)
      integer i
      integer i1
      integer j1
      integer j2
      integer j3
      integer jd1
      integer jd2
      integer jd3
      integer k
      integer kmax
      integer kmin
      integer m1
      integer m2
      integer m3
      integer md1
      integer md2
      integer min1
      integer min2
      integer min3
      integer min4
      integer min5
      integer mtri(9)
      integer n
      integer ncut
      integer num
      double precision p
      double precision phasef
      double precision plog
      double precision s
      logical setup
      double precision sig
      double precision slog
      double precision uk
      double precision ulog
      double precision vcc

      save / factrl /
      save setup

      common / factrl / fl

      data setup / .true. /

      vcc = 0.0D+00

      j1 = jd1
      j2 = jd2
      j3 = jd3

      m1 = md1
      m2 = md2
      m3 = - m1 - m2

      if ( setup ) then
        call facalc ( setup )
      end if

      i = j1 + j2 - j3
      i1 = i / 2

      if ( i .ne. 2 * i1 ) then
        return
      end if

      mtri(1) = i1
      i = j1 - j2 + j3
      i1 = i / 2

      if ( i .ne. 2 * i1 ) then
        return
      end if

      mtri(2) = i1
      i = - j1 + j2 + j3
      i1 = i / 2

      if ( i .ne. 2 * i1 ) then
        return
      end if

      mtri(3) = i1
      i = j1 + m1
      i1 = i / 2

      if ( i .ne. 2 * i1 ) then
        return
      end if

      mtri(4) = i1
      mtri(5) = ( j1 - m1 ) / 2
      i = j2 + m2
      i1 = i / 2

      if ( i .ne. 2 * i1 ) then
        return
      end if

      mtri(6) = i1
      mtri(7) = ( j2 - m2 ) / 2
      i = j3 + m3
      i1 = i / 2

      if ( i .ne. 2 * i1 ) then
        return
      end if

      mtri(8) = i1
      mtri(9) = ( j3 - m3 ) / 2

      do n = 1, 9
        if ( mtri(n) .lt. 0 ) then
          return
        end if
      end do

      min4 = j3 - j2 + m1
      if ( min4 .lt. 0 ) then
        kmin = - min4
      else
        kmin = 0
      end if

      min5 = j3 - j1 - m2
      if ( kmin .lt. - min5 ) then
        kmin = - min5
      end if
      kmin = kmin / 2

      if ( j2 - j3 + m1 .lt. 0 )then
        kmax = j1 + j2 - j3
      else
        kmax = j1 - m1
      end if

      if ( j2 + m2 .lt. kmax ) then
        kmax = j2 + m2
      else
        kmax = kmax / 2
      end if

      min1 = mtri(1) - kmin + 1
      min2 = mtri(5) - kmin + 1
      min3 = mtri(6) - kmin + 1
      min4 = min4 / 2 + kmin
      min5 = min5 / 2 + kmin
c
c  Sum the series.
c
      uk = 1.0D-10
      s = 1.0D-10
      ncut = 0
      kmax = kmax - kmin

      do k = 1, kmax

        uk = - uk 
     &    * dble ( ( min1 - k ) * ( min2 - k ) * ( min3 - k ) )
     &    / dble ( ( kmin + k ) * ( min4 + k ) * ( min5 + k ) )

        if ( 1.0D+30 .le. abs ( uk ) ) then
          uk = 1.0D-10 * uk
          s = 1.0D-10 * s
          ncut = ncut + 1
        end if

        if ( abs ( uk ) .lt. 1.0D-20 ) then
          go to 10
        end if

        s = s + uk

      end do
c
c  Calculate delta functions.
c
10    continue

      delog = 0.0D+00
      do n = 1, 9
        delog = delog + fl ( mtri(n)+1 )
      end do

      num = ( j1 + j2 + j3 ) / 2 + 2
      delog = 0.5D+00 * ( delog - fl(num) )
      ulog = -fl(kmin+1) - fl(min1) - fl(min2) - fl(min3) 
     &  - fl(min4+1) - fl(min5+1)
      plog = delog + ulog

      if ( plog .ge. -80.0D+00 .and. ncut .le. 0 ) then
        s = s * 1.0D+10
        p = exp ( plog )
        f3j = p * s
      else
        sig = sign ( 1.0D+00, s )
        s = abs ( s )
        slog = dlog ( s ) + dble ( ncut + 1 ) * dlog ( 1.0D+10 )
        f3j = sig * exp ( slog + plog )
      end if

      vcc = sqrt ( dble ( j3 + 1 ) ) * f3j * phasef ( kmin )

      return
      end
      function winej ( jd1, jd2, jd3, jd4, jd5, jd6, jd7, jd8, jd9 )

c*********************************************************************72
c
cc WINEJ ???
c
c  Modified:
c
c    04 February 2007
c
c  Parameters:
c
c    Input, integer JD1, JD2, JD3, JD4, JD5, JD6, JD7, JD8, JD9, 
c    the parameters.
c
c    Output, double precision WINEJ, ???
c
      implicit none

      double precision flk
      integer i
      integer i1
      integer j1
      integer j2
      integer j3
      integer j4
      integer j5
      integer j6
      integer j7
      integer j8
      integer j9
      integer jd1
      integer jd2
      integer jd3
      integer jd4
      integer jd5
      integer jd6
      integer jd7
      integer jd8
      integer jd9
      integer k
      integer kmax
      integer kmin
      integer kn(6)
      integer ksign
      integer kx(6)
      integer mtria(18)
      integer n
      integer nn(6)
      double precision phasef
      double precision s6j
      double precision sig
      double precision sum
      double precision term
      double precision winej

      winej = 0.0D+00

      j1 = jd1
      j2 = jd2
      j3 = jd3
      j4 = jd4
      j5 = jd5
      j6 = jd6
      j7 = jd7
      j8 = jd8
      j9 = jd9
c
c  Angular momentum coupling tests.
c
      i = - j1 + j2 + j3
      i1 = i / 2

      if ( i .ne. 2 * i1 ) then
c       write ( *, * ) '#1'
c       return
      end if

      mtria(1) = i1
      i = + j1 - j2 + j3
      i1 = i / 2

      if ( i .ne. 2 * i1 ) then
c       write ( *, * ) '#2'
c       return
      end if

      mtria(2) = i1
      i = + j1 + j2 - j3
      i1 = i / 2

      if ( i .ne. 2 * i1 ) then
c       write ( *, * ) '#3'
c       return
      end if

      mtria(3) = i1
      i = - j4 + j5 + j6
      i1 = i / 2

      if ( i .ne. 2 * i1 ) then
c       write ( *, * ) '#4'
c       return
      end if

      mtria(4) = i1
      i = + j4 - j5 + j6
      i1 = i / 2

      if ( i .ne. 2 * i1 ) then
c       write ( *, * ) '#5'
c       return
      end if

      mtria(5) = i1
      i = + j4 + j5 - j6
      i1 = i / 2

      if ( i .ne. 2 * i1 ) then
c       write ( *, * ) '#6'
c       return
      end if

      mtria(6) = i1
      i = - j7 + j8 + j9
      i1 = i / 2

      if ( i .ne. 2 * i1 ) then
c       write ( *, * ) '#7'
c       return
      end if

      mtria(7) = i1
      i = + j7 - j8 + j9
      i1 = i / 2

      if ( i .ne. 2 * i1 ) then
c       write ( *, * ) '#7.5'
c       return
      end if

      mtria(8) = i1
      i = + j7 + j8 - j9
      i1 = i / 2

      if ( i .ne. 2 * i1 ) then
c       write ( *, * ) '#8'
c       return
      end if
c
c  Now check columns.
c
      mtria(9) = i1
      i = - j1 + j4 + j7
      i1 = i / 2

      if ( i .ne. 2 * i1 ) then
c       write ( *, * ) '#9'
c       return
      end if

      mtria(10) = i1
      i = + j1 - j4 + j7
      i1 = i / 2

      if ( i .ne. 2 * i1 ) then
c       write ( *, * ) '#10'
c       return
      end if

      mtria(11) = i1
      i = + j1 + j4 - j7
      i1 = i / 2

      if ( i .ne. 2 * i1 ) then
c       write ( *, * ) '#11'
c       return
      end if

      mtria(12) = i1
      i = - j2 + j5 + j8
      i1 = i / 2

      if ( i .ne. 2 * i1 ) then
c       write ( *, * ) '#12'
c       return
      end if

      mtria(13) = i1
      i = + j2 - j5 + j8
      i1 = i / 2

      if ( i .ne. 2 * i1 ) then
c       write ( *, * ) '#13'
c       return
      end if

      mtria(14) = i1
      i = + j2 + j5 - j8
      i1 = i / 2

      if ( i .ne. 2 * i1 ) then
c       write ( *, * ) '#14'
c       return
      end if

      mtria(15) = i1
      i = - j3 + j6 + j9
      i1 = i / 2

      if ( i .ne. 2 * i1 ) then
c       write ( *, * ) '#14'
c       return
      end if

      mtria(16) = i1
      i = + j3 - j6 + j9
      i1 = i / 2

      if ( i .ne. 2 * i1 ) then
c       write ( *, * ) '#15'
c       return
      end if

      mtria(17) = i1
      i = + j3 + j6 - j9
      i1 = i / 2

      if ( i .ne. 2 * i1 ) then
c       write ( *, * ) '#16'
c       return
      end if

      mtria(18) = i1

      do n = 1, 18
        if ( mtria(n) .lt. 0 ) then
c         write ( *, * ) '#17'
c         return
        end if
      end do

      kn(1) = max ( iabs ( j2 - j6 ),
     &              iabs ( j1 - j9 ),
     &              iabs ( j4 - j8 ) )

      kn(2) = max ( iabs ( j2 - j7 ), 
     &              iabs ( j5 - j9 ),
     &              iabs ( j4 - j3 ) )

      kn(3) = max ( iabs ( j6 - j7 ), 
     &              iabs ( j5 - j1 ),
     &              iabs ( j8 - j3 ) )

      kn(4) = max ( iabs ( j6 - j1 ), 
     &              iabs ( j2 - j9 ),
     &              iabs ( j5 - j7 ) )

      kn(5) = max ( iabs ( j2 - j4 ), 
     &              iabs ( j3 - j7 ),
     &              iabs ( j6 - j8 ) )

      kn(6) = max ( iabs ( j3 - j5 ), 
     &              iabs ( j1 - j8 ),
     &              iabs ( j4 - j9 ) )

      kx(1) = min ( j2 + j6, j1 + j9, j4 + j8 )
      kx(2) = min ( j2 + j7, j5 + j9, j4 + j3 )
      kx(3) = min ( j6 + j7, j5 + j1, j8 + j3 )
      kx(4) = min ( j1 + j6, j2 + j9, j5 + j7 )
      kx(5) = min ( j2 + j4, j3 + j7, j6 + j8 )
      kx(6) = min ( j3 + j5, j1 + j8, j4 + j9 )

      do k = 1, 6
        nn(k) = kx(k) - kn(k)
      end do

      ksign = 1

      i = min ( nn(1), nn(2), nn(3), nn(4), nn(5), nn(6) )

      do k = 1, 6
        if ( i .eq. nn(k) ) then
          go to 10
        end if
      end do

      k = 6

10    continue

      kmin = kn(k) + 1
      kmax = kx(k) + 1

      if ( k .eq. 1 ) then

        go to 30

      else if ( k .eq. 2 ) then

        call i4_swap ( j1, j5 )
        call i4_swap ( j3, j8 )
        call i4_swap ( j6, j7 )
        go to 30

      else if ( k .eq. 3 ) then

        call i4_swap ( j2, j7 )
        call i4_swap ( j3, j4 )
        call i4_swap ( j5, j9 )
        go to 30

      else if ( k .eq. 4 ) then

        call i4_swap ( j1, j2 )
        call i4_swap ( j4, j5 )
        call i4_swap ( j7, j8 )
        go to 20

      else if ( k .eq. 5 ) then

        call i4_swap ( j1, j3 )
        call i4_swap ( j4, j6 )
        call i4_swap ( j7, j9 )
        go to 20

      else if ( k .eq. 6 ) then

        call i4_swap ( j2, j3 )
        call i4_swap ( j5, j6 )
        call i4_swap ( j8, j9 )

      end if

20    continue

      ksign = j1 + j2 + j3 + j4 + j5 + j6 + j7 + j8 + j9
      ksign = 1 - mod ( ksign, 4 )
c
c  Sum the series.
c
30    continue

      sum = 0.0D+00
      sig = phasef ( kmin - 1 ) * ksign
      flk = dble ( kmin )
      write ( *, * ) 'DEBUG: Kmin, Kmax = ', kmin, kmax
      do k = kmin, kmax, 2

        term = flk 
     &    * s6j ( j1, j4, j7, j8,  j9,  k-1 ) 
     &    * s6j ( j2, j5, j8, j4,  k-1, j6  )
     &    * s6j ( j3, j6, j9, k-1, j1,  j2  )

        flk = flk + 2.0D+00
        sum = sum + term

      end do

      winej = sum * sig

      return
      end
