      subroutine i4vec_copy ( n, a1, a2 )

c*********************************************************************72
c
cc I4VEC_COPY copies an I4VEC.
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
c    Input, integer N, the length of the vectors.
c
c    Input, integer A1(N), the vector to be copied.
c
c    Output, integer A2(N), a copy of A1.
c
      implicit none

      integer n

      integer a1(n)
      integer a2(n)
      integer i

      do i = 1, n
        a2(i) = a1(i)
      end do

      return
      end
      function i4vec_dot_product ( n, x, y )

c*********************************************************************72
c
cc I4VEC_DOT_PRODUCT computes the dot product of two I4VEC's.
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
c    19 December 2011
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the size of the array.
c
c    Input, integer X(N), Y(N), the arrays.
c
c    Output, integer I4VEC_DOT_PRODUCT, the dot product of X and Y.
c
      implicit none

      integer n

      integer i
      integer i4vec_dot_product
      integer value
      integer x(n)
      integer y(n)

      value = 0
      do i = 1, n
        value = value + x(i) * y(i)
      end do

      i4vec_dot_product = value

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
      subroutine partition_brute ( n, w, c, discrepancy )

c*********************************************************************72
c
cc PARTITION_BRUTE approaches the partition problem using brute force.
c
c  Discussion:
c
c    We are given a set of N integers W.
c
c    We seek to partition W into subsets W0 and W1, such that the subsets
c    have equal sums.
c
c    The "discrepancy" is the absolute value of the difference between the
c    two sums, and will be zero if we have solved the problem.
c
c    For a given set of integers, there may be zero, one, or many solutions.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    09 May 2012
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the size of the set.
c
c    Input, integer W(N), the integers.
c
c    Output, integer C(N), indicates the proposed solution.
c    C(I) is 0 for items in set W0 and 1 for items in set W1.
c
c    Output, integer DISCREPANCY, the discrepancy.
c
      implicit none

      integer n

      integer c(n)
      integer d(n)
      integer d_discrepancy
      integer discrepancy
      integer i4vec_dot_product
      integer i4vec_sum
      integer rank
      integer w(n)
      integer w_sum

      w_sum = i4vec_sum ( n, w )
      discrepancy = w_sum

      rank = -1

10    continue

        call subset_next ( n, d, rank )

        if ( rank .eq. -1 ) then
          go to 20
        end if

        d_discrepancy = 
     &    abs ( w_sum - 2 * i4vec_dot_product ( n, d, w ) )

        if ( d_discrepancy .lt. discrepancy ) then
          discrepancy = d_discrepancy
          call i4vec_copy ( n, d, c )
        end if

        if ( discrepancy .eq. 0 ) then
          go to 20
        end if

      go to 10

20    continue

      return
      end
      subroutine partition_count ( n, w, count )

c*********************************************************************72
c
cc PARTITION_COUNT counts the solutions to a partition problem.
c
c  Discussion:
c
c    We are given a set of N integers W.
c
c    We seek to partition W into subsets W0 and W1, such that the subsets
c    have equal sums.
c
c    The "discrepancy" is the absolute value of the difference between the
c    two sums, and will be zero if we have solved the problem.
c
c    For a given set of integers, there may be zero, one, or many solutions.
c
c    In the case where the weights are distinct, the count returned by this
c    function may be regarded as twice as big as it should be, since the
c    partition (W0,W1) is counted a second time as (W1,W0).  A more serious
c    overcount can occur if the set W contains duplicate elements - in the
c    extreme case, W might be entirely 1's, in which case there is really
c    only one (interesting) solution, but this function will count many.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    10 May 2012
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the size of the set.
c
c    Input, integer W(N), the integers.
c
c    Output, integer COUNT, the number of solutions.
c
      implicit none

      integer n

      integer c(n)
      integer count
      integer discrepancy
      integer i4vec_dot_product
      integer i4vec_sum
      integer rank
      integer w(n)
      integer w_sum

      w_sum = i4vec_sum ( n, w )

      rank = -1
      count = 0

10    continue

        call subset_next ( n, c, rank )

        if ( rank .eq. -1 ) then
          go to 20
        end if

        discrepancy = 
     &    abs ( w_sum - 2 * i4vec_dot_product ( n, c, w ) )

        if ( discrepancy .eq. 0 ) then
          count = count + 1
        end if

      go to 10

20    continue

      return
      end
      subroutine subset_next ( n, t, rank )

c*********************************************************************72
c
cc SUBSET_NEXT computes the subset lexicographic successor.
c
c  Discussion:
c
c    This is a lightly modified version of "subset_lex_successor()" from COMBO.
c
c  Example:
c
c    On initial call, N is 5 and the input value of RANK is -1.
c    Then here are the successive outputs from the program:
c
c   Rank   T1   T2   T3   T4   T5
c   ----   --   --   --   --   --
c      0    0    0    0    0    0
c      1    0    0    0    0    1
c      2    0    0    0    1    0
c      3    0    0    0    1    1
c     ..   ..   ..   ..   ..   ..
c     30    1    1    1    1    0
c     31    1    1    1    1    1
c     -1    0    0    0    0    0  <-- Reached end of cycle.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    10 May 2012
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Donald Kreher, Douglas Simpson,
c    Combinatorial Algorithms,
c    CRC Press, 1998,
c    ISBN: 0-8493-3988-X,
c    LC: QA164.K73.
c
c  Parameters:
c
c    Input, integer N, the number of elements in the master set.
c    N must be positive.
c
c    Input/output, integer T(N), describes a subset.  T(I) is 0 if
c    the I-th element of the master set is not in the subset, and is
c    1 if the I-th element is part of the subset.
c    On input, T describes a subset.
c    On output, T describes the next subset in the ordering.
c
c    Input/output, integer RANK, the rank.
c    If RANK = -1 on input, then the routine understands that this is
c    the first call, and that the user wishes the routine to supply
c    the first element in the ordering, which has RANK = 0.
c    In general, the input value of RANK is increased by 1 for output,
c    unless the very last element of the ordering was input, in which
c    case the output value of RANK is -1.
c
      implicit none

      integer n

      integer i
      integer rank
      integer t(n)
c
c  Return the first element.
c
      if ( rank == -1 ) then
        do i = 1, n
          t(i) = 0
        end do
        rank = 0
        return
      end if

      do i = n, 1, -1

        if ( t(i) .eq. 0 ) then
          t(i) = 1
          rank = rank + 1
          return
        else
          t(i) = 0
        end if

      end do

      rank = -1

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
