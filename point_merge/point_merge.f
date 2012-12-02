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
      subroutine point_radial_tol_unique_count ( m, n, a, tol, seed, 
     &  unique_num )

c*********************************************************************72
c
cc POINT_RADIAL_TOL_UNIQUE_COUNT counts the tolerably unique points.
c
c  Discussion:
c
c    The input data is an M x N array A, representing the M-dimensional
c    coordinates of N points.
c
c    The output is the number of tolerably unique points in the list.
c
c    This program performs the same task as POINT_TOL_UNIQUE_COUNT.
c    But that program is guaranteed to use N^2 comparisons.
c
c    It is hoped that this function, on the other hand, will tend
c    to use O(N) comparisons after an O(NLog(N)) sort.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    24 July 2010
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer M, the number of rows.
c
c    Input, integer N, the number of columns.
c
c    Input, double precision A(M,N), the array of N columns of data.
c
c    Input, double precision TOL, a tolerance for equality.
c
c    Input/output, integer SEED, a seed for the random
c    number generator.
c
c    Output, integer UNIQUE_NUM, the number of tolerably
c    unique points.
c
      implicit none

      integer m
      integer n

      double precision a(m,n)
      double precision dist
      integer hi
      integer i
      integer indx(n)
      integer j
      integer k
      double precision r(n)
      double precision r8vec_sum
      integer seed
      double precision tol
      logical unique(n)
      integer unique_num
      double precision w(n)
      double precision w_sum
      double precision z(m)

      if ( n .le. 0 ) then
        unique_num = 0
        return
      end if
c
c  Assign a base point Z randomly in the convex hull.
c
      call r8vec_uniform_01 ( n, seed, w )

      w_sum = r8vec_sum ( n, w )
      do i = 1, n
        w(i) = w(i) / w_sum
      end do

      do i = 1, m
        z(i) = 0.0D+00
        do j = 1, n
          z(i) = z(i) + a(i,j) * w(j)
        end do
      end do
c
c  Compute the radial distance R of each point to Z.
c
      do j = 1, n
        r(j) = 0.0D+00
        do i = 1, m
          r(j) = r(j) + ( a(i,j) - z(i) )**2
        end do
        r(j) = sqrt ( r(j) )
      end do
c
c  Implicitly sort the R array.
c
      call r8vec_sort_heap_index_a ( n, r, indx )
c
c  To determine if a point I is tolerably unique, we only have to check
c  whether it is distinct from all points J such that R(I) <= R(J) <= R(J)+TOL.
c
      unique_num = 0

      do i = 1, n
        unique(1:n) = .true.
      end do

      do i = 1, n

        if ( unique(indx(i)) ) then
c
c  Point INDX(I) is unique, in that no earlier point is near it.
c
          unique_num = unique_num + 1
c
c  Look for later points which are close to point INDX(I)
c  in terms of R.
c
          hi = i

10        continue

          if ( hi .lt. n ) then
            if ( r(indx(hi+1)) .le. r(indx(i)) + tol ) then
              hi = hi + 1
              go to 10
            end if
          end if
c
c  Points INDX(I+1) through INDX(HI) have an R value close to 
c  point INDX(I).  Are they truly close to point INDEX(I)?
c
          do j = i + 1, hi
            if ( unique(indx(j)) ) then
              dist = 0.0D+00
              do k = 1, m
                dist = dist + ( a(k,indx(i)) - a(k,indx(j)) )**2
              end do
              dist = sqrt ( dist )
              if ( dist .le. tol ) then
                unique(indx(j)) = .false.
              end if
            end if
          end do

        end if

      end do

      return
      end
      subroutine point_radial_tol_unique_index ( m, n, a, tol, seed, 
     &  unique_num, undx, xdnu )

c*********************************************************************72
c
cc POINT_RADIAL_TOL_UNIQUE_INDEX indexes the tolerably unique points.
c
c  Discussion:
c
c    The input data is an M x N array A, representing the M-dimensional
c    coordinates of N points.
c
c    The output is:
c    * the number of tolerably unique points in the list;
c    * the index, in the list of unique items, of the representatives 
c      of each point;
c    * the index, in A, of the tolerably unique representatives.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    28 July 2010
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer M, the number of rows.
c
c    Input, integer N, the number of columns.
c
c    Input, double precision A(M,N), the array of N columns of data.
c
c    Input, double precision TOL, a tolerance for equality.
c
c    Input/output, integer SEED, a seed for the random
c    number generator.
c
c    Output, integer UNIQUE_NUM, the number of tolerably
c    unique points.
c
c    Output, integer UNDX(UNIQUE_NUM), the index, in A, of the 
c    tolerably unique points.
c
c    Output, integer XDNU(N), the index, in UNDX, of the 
c    tolerably unique point that "represents" this point.
c
      implicit none

      integer m
      integer n

      double precision a(m,n)
      double precision dist
      integer hi
      integer i
      integer indx(n)
      integer j
      integer k
      double precision r(n)
      double precision r8vec_sum
      integer seed
      double precision tol
      integer undx(n)
      logical unique(n)
      integer unique_num
      double precision w(n)
      double precision w_sum
      integer xdnu(n)
      double precision z(m)

      if ( n .le. 0 ) then
        unique_num = 0
        return
      end if
c
c  Assign a base point Z randomly in the convex hull.
c
      call r8vec_uniform_01 ( n, seed, w )

      w_sum = r8vec_sum ( n, w )
      do i = 1, n
        w(i) = w(i) / w_sum
      end do

      do i = 1, m
        z(i) = 0.0D+00
        do j = 1, n
          z(i) = z(i) + a(i,j) * w(j)
        end do
      end do
c
c  Compute the radial distance R of each point to Z.
c
      do j = 1, n
        r(j) = 0.0D+00
        do i = 1, m
          r(j) = r(j) + ( a(i,j) - z(i) )**2
        end do
        r(j) = sqrt ( r(j) )
      end do
c
c  Implicitly sort the R array.
c
      call r8vec_sort_heap_index_a ( n, r, indx )
c
c  To determine if a point I is tolerably unique, we only have to check
c  whether it is distinct from all points J such that R(I) <= R(J) <= R(J)+TOL.
c
      unique_num = 0
      do i = 1, n
        unique(i) = .true.
      end do

      do i = 1, n

        if ( unique(indx(i)) ) then
c
c  Point INDX(I) is unique, in that no earlier point is near it.
c
          unique_num = unique_num + 1
          xdnu(indx(i)) = unique_num
          undx(unique_num) = indx(i)
c
c  Look for later points which are close to point INDX(I)
c  in terms of R.
c
          hi = i

10        continue

          if ( hi .lt. n ) then
            if ( r(indx(hi+1)) .le. r(indx(i)) + tol ) then
              hi = hi + 1
              go to 10
            end if
          end if
c
c  Points INDX(I+1) through INDX(HI) have an R value close to 
c  point INDX(I).  Are they truly close to point INDEX(I)?
c
          do j = i + 1, hi
            if ( unique(indx(j)) ) then
              dist = 0.0D+00
              do k = 1, m
                dist = dist + ( a(k,indx(i)) - a(k,indx(j)) )**2
              end do
              dist = sqrt ( dist )
              if ( dist .le. tol ) then
                unique(indx(j)) = .false.
                xdnu(indx(j)) = xdnu(indx(i))
              end if
            end if
          end do

        end if

      end do

      return
      end
      subroutine point_radial_unique_count ( m, n, a, seed, unique_num )

c*********************************************************************72
c
cc POINT_RADIAL_UNIQUE_COUNT counts the unique points.
c
c  Discussion:
c
c    The input data is an M x N array A, representing the M-dimensional
c    coordinates of N points.
c
c    The output is the number of unique points in the list.
c
c    This program performs the same task as POINT_UNIQUE_COUNT, and
c    carries out more work.  Hence, it is not a substitute for
c    POINT_UNIQUE_COUNT.  Instead, it is intended to be a starting point
c    for a similar program which includes a tolerance.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    24 July 2010
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer M, the number of rows.
c
c    Input, integer N, the number of columns.
c
c    Input, double precision A(M,N), the array of N columns of data.
c
c    Input/output, integer SEED, a seed for the random
c    number generator.
c
c    Output, integer UNIQUE_NUM, the number of unique points.
c
      implicit none

      integer m
      integer n

      double precision a(m,n)
      logical equal
      integer hi
      integer i
      integer indx(n)
      integer j
      integer j1
      integer j2
      integer lo
      double precision r(n)
      double precision r8vec_sum
      integer seed
      integer unique_index
      integer unique_num
      double precision w(n)
      double precision w_sum
      double precision z(m)

      if ( n .le. 0 ) then
        unique_num = 0
        return
      end if
c
c  Assign a base point Z randomly in the convex hull.
c
      call r8vec_uniform_01 ( n, seed, w )

      w_sum = r8vec_sum ( n, w )
      do i = 1, n
        w(i) = w(i) / w_sum
      end do

      do i = 1, m
        z(i) = 0.0D+00
        do j = 1, n
          z(i) = z(i) + a(i,j) * w(j)
        end do
      end do
c
c  Compute the radial distance R of each point to Z.
c
      do j = 1, n
        r(j) = 0.0D+00
        do i = 1, m
          r(j) = r(j) + ( a(i,j) - z(i) )**2
        end do
        r(j) = sqrt ( r(j) )
      end do
c
c  Implicitly sort the R array.
c
      call r8vec_sort_heap_index_a ( n, r, indx )
c
c  To determine if a point is unique, we only have to check
c  whether it is distinct from all points with the same
c  R value and lower ordering.
c
      unique_num = 0
      hi = 0

10    continue

      if ( hi .lt. n ) then
c
c  Advance LO.
c
        lo = hi + 1
c
c  Extend HI.
c
        hi = lo

20      continue

        if ( hi .lt. n ) then
          if ( r(indx(hi+1)) .eq. r(indx(lo)) ) then
            hi = hi + 1
            go to 20
          end if
        end if
c
c  Points INDX(LO) through INDX(HI) have same R value.
c
c  Find the unique ones.
c
        unique_num = unique_num + 1

        do j1 = lo + 1, hi

          do j2 = lo, j1 - 1

            equal = .true.

            do i = 1, m
              if ( a(i,indx(j2)) .ne. a(i,indx(j1)) ) then
                equal = .false.
                go to 30
              end if
            end do

30          continue

            if ( equal ) then
              go to 40
            end if

          end do

          unique_num = unique_num + 1

40        continue

        end do

        go to 10

      end if

      return
      end
      subroutine point_tol_unique_count ( m, n, a, tol, unique_num )

c*********************************************************************72
c
cc POINT_TOL_UNIQUE_COUNT counts the tolerably unique points.
c
c  Discussion:
c
c    The input data is an M x N array A, representing the M-dimensional
c    coordinates of N points.
c
c    This function uses a simple but expensive approach.  The first point
c    is accepted as unique.  Each subsequent point is accepted as unique
c    only if it is at least a tolerance away from all accepted unique points.
c    This means the expected amount of work is O(N^2).
c
c    The output is the number of unique points in the list.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    24 July 2010
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer M, the number of rows.
c
c    Input, integer N, the number of columns.
c
c    Input, double precision A(M,N), the array of N columns of data.
c
c    Input, double precision TOL, a tolerance.
c
c    Output, integer UNIQUE_NUM, the number of unique points.
c
      implicit none

      integer m
      integer n

      double precision a(m,n)
      double precision dist
      integer i
      integer j
      integer k
      double precision tol
      logical unique(n)
      integer unique_num

      do i = 1, n
        unique(i) = .true.
      end do

      unique_num = n

      do i = 2, n

        do j = 1, i - 1

          if ( unique(j) ) then

            dist = 0.0D+00
            do k = 1, m
              dist = dist + ( a(k,i) - a(k,j) )**2
            end do
            dist = sqrt ( dist )

            if ( dist .le. tol ) then
              unique(i) = .false.
              unique_num = unique_num - 1
              go to 10
            end if

          end if

        end do

10      continue

      end do

      return
      end
      subroutine point_tol_unique_index ( m, n, a, tol, unique_num, 
     &  xdnu )

c*********************************************************************72
c
cc POINT_TOL_UNIQUE_INDEX indexes the tolerably unique points.
c
c  Discussion:
c
c    This routine uses an algorithm that is O(N^2).
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    24 July 2010
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer M, the dimension of the data values.
c
c    Input, integer N, the number of data values.
c
c    Input, double precision A(M,N), the data values.
c
c    Input, double precision TOL, a tolerance for equality.
c
c    Output, integer UNIQUE_NUM, the number of tolerably
c    unique points.
c
c    Output, integer XDNU(N), the index, in A, of the tolerably unique 
c    point that "represents" this point.
c
      implicit none

      integer m
      integer n

      double precision a(m,n)
      double precision dist
      integer i
      integer j
      integer k
      double precision tol
      logical unique(n)
      integer unique_num
      integer xdnu(n)

      do i = 1, n
        unique(i) = .true.
      end do

      do i = 1, n
        xdnu(i) = i
      end do
      unique_num = n

      i = 1
      xdnu(i) = 1

      do i = 2, n

        do j = 1, i - 1

          if ( unique(j) ) then

            dist = 0.0D+00
            do k = 1, m
              dist = dist + ( a(k,i) - a(k,j) )**2
            end do
            dist = sqrt ( dist )

            if ( dist .le. tol ) then
              unique(i) = .false.
              unique_num = unique_num - 1
              xdnu(i) = j
              go to 10
            end if
          end if

        end do

10      continue

      end do

      return
      end
      subroutine point_unique_count ( m, n, a, unique_num )

c*********************************************************************72
c
cc POINT_UNIQUE_COUNT counts the unique points.
c
c  Discussion:
c
c    The input data is an M x N array A, representing the M-dimensional
c    coordinates of N points.
c
c    The algorithm relies on the fact that, in a sorted list, points that
c    are exactly equal must occur consecutively.
c
c    The output is the number of unique points in the list.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    24 July 2010
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer M, the number of rows.
c
c    Input, integer N, the number of columns.
c
c    Input, double precision A(M,N), the array of N columns of data.
c
c    Output, integer UNIQUE_NUM, the number of unique points.
c
      implicit none

      integer m
      integer n

      double precision a(m,n)
      integer i
      integer indx(n)
      integer j
      integer unique_index
      integer unique_num

      if ( n .le. 0 ) then
        unique_num = 0
        return
      end if
c
c  Implicitly sort the array.
c
      call r8col_sort_heap_index_a ( m, n, a, indx )
c
c  Two points are considered equal only if they exactly match.
c  In that case, equal points can only occur as consecutive items
c  in the sorted list.   This makes counting easy.
c
      unique_num = 1
      unique_index = indx(1)

      do j = 2, n
        do i = 1, m
          if ( a(i,unique_index) .ne. a(i,indx(j)) ) then
            unique_num = unique_num + 1
            unique_index = indx(j)
            go to 10
          end if
        end do

10      continue

      end do

      return
      end
      subroutine point_unique_index ( m, n, a, unique_num, undx, xdnu )

c*********************************************************************72
c
cc POINT_UNIQUE_INDEX indexes unique points.
c
c  Discussion:
c
c    An R8COL is an M by N array of R8's, regarded as an array of N columns,
c    each of length M.
c
c    The goal of this routine is to determine a vector UNDX,
c    which points to the unique elements of A, in sorted order,
c    and a vector XDNU, which identifies, for each entry of A, the index of
c    the unique sorted element of A.
c
c    This is all done with index vectors, so that the elements of
c    A are never moved.
c
c    The first step of the algorithm requires the indexed sorting
c    of A, which creates arrays INDX and XDNI.  (If all the entries
c    of A are unique, then these arrays are the same as UNDX and XDNU.)
c
c    We then use INDX to examine the entries of A in sorted order,
c    noting the unique entries, creating the entries of XDNU and
c    UNDX as we go.
c
c    Once this process has been completed, the object A could be
c    replaced by a compressed object XU, containing the unique entries
c    of X in sorted order, using the formula
c
c      XU(*) = A(UNDX(*)).
c
c    We could then, if we wished, reconstruct the entire vector A, or
c    any element of it, by index, as follows:
c
c      A(I) = XU(XDNU(I)).
c
c    We could then replace A by the combination of XU and XDNU.
c
c    Later, when we need the I-th entry of A, we can locate it as
c    the XDNU(I)-th entry of XU.
c
c    Here is an example of a vector A, the sort and inverse sort
c    index vectors, and the unique sort and inverse unique sort vectors
c    and the compressed unique sorted vector.
c
c      I    A   Indx  Xdni      XU   Undx  Xdnu
c    ----+-----+-----+-----+--------+-----+-----+
c      1 | 11.     1     1 |    11.     1     1
c      2 | 22.     3     5 |    22.     2     2
c      3 | 11.     6     2 |    33.     4     1
c      4 | 33.     9     8 |    55.     5     3
c      5 | 55.     2     9 |                  4
c      6 | 11.     7     3 |                  1
c      7 | 22.     8     6 |                  2
c      8 | 22.     4     7 |                  2
c      9 | 11.     5     4 |                  1
c
c    INDX(2) = 3 means that sorted item(2) is A(3).
c    XDNI(2) = 5 means that A(2) is sorted item(5).
c
c    UNDX(3) = 4 means that unique sorted item(3) is at A(4).
c    XDNU(8) = 2 means that A(8) is at unique sorted item(2).
c
c    XU(XDNU(I))) = A(I).
c    XU(I)        = A(UNDX(I)).
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
c    Input, integer M, the dimension of the data values.
c
c    Input, integer N, the number of data values.
c
c    Input, double precision A(M,N), the data values.
c
c    Input, integer UNIQUE_NUM, the number of unique values 
c    in A.  This value is only required for languages in which the size of
c    UNDX must be known in advance.
c
c    Output, integer UNDX(UNIQUE_NUM), the UNDX vector.
c
c    Output, integer XDNU(N), the XDNU vector.
c
      implicit none

      integer m
      integer n
      integer unique_num

      double precision a(m,n)
      double precision diff
      integer i
      integer indx(n)
      integer j
      integer k
      integer undx(unique_num)
      integer xdnu(n)
c
c  Implicitly sort the array.
c
      call r8col_sort_heap_index_a ( m, n, a, indx )
c
c  Walk through the implicitly sorted array.
c
      i = 1
      j = 1
      undx(j) = indx(i)
      xdnu(indx(i)) = j

      do i = 2, n

        diff = 0.0D+00
        do k = 1, m
          diff = max ( diff, abs ( a(k,indx(i)) - a(k,undx(j)) ) )
        end do

        if ( 0.0D+00 .lt. diff ) then
          j = j + 1
          undx(j) = indx(i)
        end if

        xdnu(indx(i)) = j

      end do

      return
      end
      subroutine r8col_duplicates ( m, n, n_unique, seed, a )

c*********************************************************************72
c
cc R8COL_DUPLICATES generates an R8COL with some duplicate columns.
c
c  Discussion:
c
c    An R8COL is an M by N array of R8's, regarded as an array of N columns,
c    each of length M.
c
c    This routine generates a random R8COL with a specified number of
c    duplicate columns.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    21 July 2010
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer M, the number of rows in each column of A.
c
c    Input, integer N, the number of columns in A.
c
c    Input, integer N_UNIQUE, the number of unique columns in A.
c    1 <= N_UNIQUE <= N.
c
c    Input/output, integer SEED, a seed for the random
c    number generator.
c
c    Output, double precision A(M,N), the array.
c
      implicit none

      integer m
      integer n

      double precision a(m,n)
      integer i
      integer i4_uniform
      integer j1
      integer j2
      integer n_unique
      integer seed
      double precision temp(m)

      if ( n_unique .lt. 1 .or. n .lt. n_unique ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'R8COL_DUPLICATES - Fatal error!'
        write ( *, '(a)' ) '  1 <= N_UNIQUE <= N is required.'
        stop
      end if

      call r8mat_uniform_01 ( m, n_unique, seed, a )
c
c  Randomly copy unique columns.
c
      do j1 = n_unique + 1, n
        j2 = i4_uniform ( 1, n_unique, seed )
        do i = 1, m
          a(i,j1) = a(i,j2)
        end do
      end do
c
c  Permute the columns.
c
      do j1 = 1, n
        j2 = i4_uniform ( j1, n, seed )
        do i = 1, m
          temp(i) = a(i,j1)
          a(i,j1) = a(i,j2)
          a(i,j2) = temp(i)
        end do
      end do

      return
      end
      subroutine r8col_sort_heap_index_a ( m, n, a, indx )

c*********************************************************************72
c
cc R8COL_SORT_HEAP_INDEX_A does an indexed heap ascending sort of an R8COL.
c
c  Discussion:
c
c    An R8COL is an M by N array of R8's, regarded as an array of N columns,
c    each of length M.
c
c    The sorting is not actually carried out.  Rather an index array is
c    created which defines the sorting.  This array may be used to sort
c    or index the array, or to sort or index related arrays keyed on the
c    original array.
c
c    A(*,J1) < A(*,J2) if the first nonzero entry of A(*,J1)-A(*,J2) 
c    is negative.
c
c    Once the index array is computed, the sorting can be carried out
c    "implicitly:
c
c      A(*,INDX(*)) is sorted,
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
c    Input, integer M, the number of rows in each column of A.
c
c    Input, integer N, the number of columns in A.
c
c    Input, double precision A(M,N), the array.
c
c    Output, integer INDX(N), the sort index.  The I-th element
c    of the sorted array is column INDX(I).
c
      implicit none

      integer m
      integer n

      double precision a(m,n)
      double precision column(m)
      integer i
      integer indx(n)
      integer indxt
      integer ir
      integer isgn
      integer j
      integer l

      if ( n .lt. 1 ) then
        return
      end if

      do i = 1, n
        indx(i) = i
      end do

      if ( n .eq. 1 ) then
        return
      end if

      l = ( n / 2 ) + 1
      ir = n

10    continue

        if ( 1 .lt. l ) then

          l = l - 1
          indxt = indx(l)
          do i = 1, m
            column(i) = a(i,indxt)
          end do

        else

          indxt = indx(ir)
          do i = 1, m
            column(i) = a(i,indxt)
          end do
          indx(ir) = indx(1)
          ir = ir - 1

          if ( ir .eq. 1 ) then
            indx(1) = indxt
            go to 30
          end if

        end if

        i = l
        j = l + l

20      continue

        if ( j .le. ir ) then

          if ( j .lt. ir ) then

            call r8vec_compare ( m, a(1,indx(j)), a(1,indx(j+1)), isgn )

            if ( isgn .lt. 0 ) then
              j = j + 1
            end if

          end if

          call r8vec_compare ( m, column, a(1,indx(j)), isgn )

          if ( isgn .lt. 0 ) then
            indx(i) = indx(j)
            i = j
            j = j + j
          else
            j = ir + 1
          end if

          go to 20

        end if

        indx(i) = indxt

      go to 10

30    continue

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
      subroutine r8mat_uniform_01 ( m, n, seed, r )

c*********************************************************************72
c
cc R8MAT_UNIFORM_01 returns a unit pseudorandom R8MAT.
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
c    Input, integer M, N, the number of rows and columns in the array.
c
c    Input/output, integer SEED, the "seed" value, which should NOT be 0.
c    On output, SEED has been updated.
c
c    Output, double precision R(M,N), the array of pseudorandom values.
c
      implicit none

      integer m
      integer n

      integer i
      integer j
      integer k
      integer seed
      double precision r(m,n)

      if ( seed .eq. 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'R8MAT_UNIFORM_01 - Fatal error!'
        write ( *, '(a)' ) '  Input value of SEED = 0.'
        stop
      end if

      do j = 1, n

        do i = 1, m

          k = seed / 127773

          seed = 16807 * ( seed - k * 127773 ) - k * 2836

          if ( seed .lt. 0 ) then
            seed = seed + 2147483647
          end if

          r(i,j) = dble ( seed ) * 4.656612875D-10

        end do
      end do

      return
      end
      subroutine r8vec_compare ( n, a1, a2, isgn )

c*********************************************************************72
c
cc R8VEC_COMPARE compares two R8VEC's.
c
c  Discussion:
c
c    An R8VEC is a vector of R8 values.
c
c    The lexicographic ordering is used.
c
c  Example:
c
c    Input:
c
c      A1 = ( 2.0, 6.0, 2.0 )
c      A2 = ( 2.0, 8.0, 12.0 )
c
c    Output:
c
c      ISGN = -1
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    23 February 1999
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the number of entries in the vectors.
c
c    Input, double precision A1(N), A2(N), the vectors to be compared.
c
c    Output, integer ISGN, the results of the comparison:
c    -1, A1 < A2,
c     0, A1 = A2,
c    +1, A1 > A2.
c
      implicit none

      integer n

      double precision a1(n)
      double precision a2(n)
      integer isgn
      integer k

      isgn = 0

      k = 1

10    continue

      if ( k .le. n ) then

        if ( a1(k) .lt. a2(k) ) then
          isgn = -1
          return
        else if ( a2(k) .lt. a1(k) ) then
          isgn = + 1
          return
        end if

        k = k + 1

        go to 10

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
      subroutine r8vec_sort_heap_index_a ( n, a, indx )

c*********************************************************************72
c
cc R8VEC_SORT_HEAP_INDEX_A does an indexed heap ascending sort of an R8VEC.
c
c  Discussion:
c
c    An R8VEC is a vector of R8's.
c
c    The sorting is not actually carried out.  Rather an index array is
c    created which defines the sorting.  This array may be used to sort
c    or index the array, or to sort or index related arrays keyed on the
c    original array.
c
c    Once the index array is computed, the sorting can be carried out
c    "implicitly:
c
c      A(INDX(I:N)) is sorted,
c
c    or explicitly, by the call
c
c      call r8vec_permute ( n, indx, a )
c
c    after which A(1:N) is sorted.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    01 June 2009
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the number of entries in the array.
c
c    Input, double precision A(N), an array to be index-sorted.
c
c    Output, integer INDX(N), the sort index.  The
c    I-th element of the sorted array is A(INDX(I)).
c
      implicit none

      integer n

      double precision a(n)
      double precision aval
      integer i
      integer indx(n)
      integer indxt
      integer ir
      integer j
      integer l

      if ( n .lt. 1 ) then
        return
      end if

      do i = 1, n
        indx(i) = i
      end do

      if ( n .eq. 1 ) then
        return
      end if

      l = n / 2 + 1
      ir = n

10    continue

        if ( 1 .lt. l ) then

          l = l - 1
          indxt = indx(l)
          aval = a(indxt)

        else

          indxt = indx(ir)
          aval = a(indxt)
          indx(ir) = indx(1)
          ir = ir - 1

          if ( ir .eq. 1 ) then
            indx(1) = indxt
            go to 30
          end if

        end if

        i = l
        j = l + l

20      continue

        if ( j .le. ir ) then

          if ( j .lt. ir ) then
            if ( a(indx(j)) .lt. a(indx(j+1)) ) then
              j = j + 1
            end if
          end if

          if ( aval .lt. a(indx(j)) ) then
            indx(i) = indx(j)
            i = j
            j = j + j
          else
            j = ir + 1
          end if

          go to 20

        end if

        indx(i) = indxt

      go to 10

30    continue

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
