      subroutine i4_to_digits_binary ( i, n, c )

c*********************************************************************72
c
cc I4_TO_DIGITS_BINARY produces the binary digits of an I4.
c
c  Discussion:
c
c    An I4 is an integer.
c
c  Example:
c
c     I    N     C               Binary
c    --  ---   ---         ------------
c     0    1   0                      0
c     0    2   0, 0                  00
c     1    3   1, 0, 0              100
c     2    3   0, 1, 0              010
c     3    3   1, 1, 0              011
c     4    3   0, 0, 1              100
c     8    3   0, 0, 0           (1)000
c     8    5   0, 0, 0, 1, 0      01000
c    -8    5   0, 0, 0, 1, 0  (-) 01000
c
c     0    3   0, 0, 0
c     1    3   1, 0, 0
c     2    3   0, 1, 0
c     3    3   1, 1, 0
c     4    3   0, 0, 1
c     5    3   1, 0, 1
c     6    3   0, 1, 1
c     7    3   1, 1, 1
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
c    Input, integer I, an integer to be represented.
c
c    Input, integer N, the number of binary digits to produce.
c
c    Output, integer C(N), the first N binary digits of I,
c    with C(1) being the units digit.
c
      implicit none

      integer n

      integer c(n)
      integer i
      integer i_copy
      integer j

      i_copy = abs ( i )

      do j = 1, n

        c(j) = mod ( i_copy, 2 )
        i_copy = i_copy / 2

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
      subroutine subset_sum_count ( n, w, t, ind_min, ind_max, 
     &  solution_num )

c*********************************************************************72
c
cc SUBSET_SUM_COUNT counts solutions to the subset sum problem in a range.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    02 February 2012
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the size of the set.
c
c    Input, integer W(N), a set of weights.  The length of this
c    array must be no more than 31.
c
c    Input, integer T, the target value.
c
c    Input, integer IND_MIN, IND_MAX, the lower and upper
c    limits to be searched.  0 <= IND_MIN <= IND_MAX <= (2^N)-1.
c
c    Output, integer SOLUTION_NUM, the number of distinct
c    solutions of the subset sum problem found within the given range.
c
      implicit none

      integer n

      integer c(n)
      integer i4vec_dot_product
      integer ind
      integer ind_max
      integer ind_max2
      integer ind_min
      integer ind_min2
      integer solution_num
      integer t
      integer w(n)
c
c  Check the data.
c
      if ( n .lt. 1 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'SUBSET_SUM_COUNT - Fatal error!'
        write ( *, '(a)' ) '  N < 1.'
        stop
      end if

      if ( 31 .lt. n ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'SUBSET_SUM_COUNT - Fatal error!'
        write ( *, '(a)' ) '  31 < N.'
        stop
      end if

      ind_min2 = max ( ind_min, 0 )
      ind_max2 = min ( ind_max, ( 2 ** n ) - 1 )
c
c  Run through the range.
c
      write ( *, '(a)' ) ' '
      write ( *, '(a,i8)' ) '  Searching from IND_MIN = ', ind_min2
      write ( *, '(a,i8)' ) '  through IND_MAX = ', ind_max2

      solution_num = 0

      do ind = ind_min2, ind_max2
c
c  Convert INDEX into vector of indices in W.
c
        call i4_to_digits_binary ( ind, n, c )
c
c  If the sum of those weights matches the target, return combination.
c
        if ( i4vec_dot_product ( n, c, w ) .eq. t ) then
          solution_num = solution_num + 1
        end if

      end do

      return
      end
      subroutine subset_sum_find ( n, w, t, ind_min, ind_max, ind, c )

c*********************************************************************72
c
cc SUBSET_SUM seeks a subset of a set that has a given sum.
c
c  Discussion:
c
c    This function tries to compute a target value as the sum of
c    a selected subset of a given set of weights.
c
c    This function works by brute force, that is, it tries every
c    possible subset to see if it sums to the desired value.
c
c    Given N weights, every possible selection can be described by
c    one of the N-digit binary numbers from 0 to 2^N-1.
c
c    This function includes a range, which allows the user to
c    control which subsets are to be checked.  Thus, if there are
c    N weights, specifying a range of [ 0, 2^N-1] indicates that
c    all subsets should be checked.  On the other hand, this full
c    range could be broken down into smaller subranges, each of
c    which could be checked independently.
c
c    It is possible that, in the given range, there may be multiple
c    solutions of the problem.  This function will only return
c    one such solution, if found.  However, the function may be called
c    again, with an appropriate restriction of the range, to continue
c    the search for other solutions.
c
c  Example:
c
c    w = [ 1, 2, 4, 8, 16, 32 ];
c    t = 22;
c    r = [ 0, 2^6 - 1 ];
c
c    call subset_sum ( w, t, r, c, ind )
c
c    c = [ 2, 3, 5 ]
c    index = 22
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    02 February 2012
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the size of the set.
c
c    Input, integer W(N), a set of weights.  The length of this
c    array must be no more than 31.
c
c    Input, integer T, the target value.
c
c    Input, integer IND_MIN, IND_MAX, the lower and upper
c    limits to be searched.  0 <= IND_MIN <= IND_MAX <= (2^N)-1.
c
c    Output, integer IND, the index of the solution.
c    If IND is -1, no solution was found in the range.
c
c    Output, integer C(N), indicates the solution, assuming
c    that IND is not -1.  In that case, the sum T is made by selecting
c    those weights W(I) for which C(I) is 1.  In fact,
c    T = sum ( 1 <= I <= N ) C(I) * W(I).
c
      implicit none

      integer n

      integer c(n)
      integer i4vec_dot_product
      integer ind
      integer ind_max
      integer ind_max2
      integer ind_min
      integer ind_min2
      integer t
      integer w(n)
c
c  Check the data.
c
      if ( n .lt. 1 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'SUBSET_SUM - Fatal error!'
        write ( *, '(a)' ) '  N < 1.'
        stop
      end if

      if ( 31 .lt. n ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'SUBSET_SUM - Fatal error!'
        write ( *, '(a)' ) '  31 < N.'
        stop
      end if

      ind_min2 = max ( ind_min, 0 )
      ind_max2 = min ( ind_max, ( 2 ** n ) - 1 )
c
c  Run through the range.
c
      write ( *, '(a)' ) ' '
      write ( *, '(a,i8)' ) '  Searching from IND_MIN = ', ind_min2
      write ( *, '(a,i8)' ) '  through IND_MAX = ', ind_max2

      do ind = ind_min2, ind_max2
c
c  Convert INDEX into vector of indices in W.
c
        call i4_to_digits_binary ( ind, n, c )
c
c  If the sum of those weights matches the target, return combination.
c
        if ( i4vec_dot_product ( n, c, w ) .eq. t ) then
          return
        end if

      end do

      ind = - 1

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
