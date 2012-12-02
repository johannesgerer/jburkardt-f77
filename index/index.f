      function index0 ( i_min, i, i_max )

c*********************************************************************72
c
cc INDEX0 indexes a 1D vector using a zero base.
c
c  Discussion:
c
c    Index       Element
c    ---------   --------
c    0           I_MIN
c    INDEX0      I
c   (INDEX_MAX)  I_MAX
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    27 November 2012
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
c    Output, integer INDEX0, the index of element I.
c
      implicit none

      integer i
      integer i_max
      integer i_min
      integer index_min
      parameter ( index_min = 0 )
      integer index0

      index0 = index_min + ( i - i_min )

      return
      end
      function index01 ( i_min, i, i_max, j_min, j, j_max )

c*********************************************************************72
c
cc INDEX01 indexes a 2D array by columns, with a zero base.
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
c    27 November 2012
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
c    Output, integer INDEX01, the index of element (I,J).
c
      implicit none

      integer i
      integer i_max
      integer i_min
      integer index_min
      parameter ( index_min = 0 )
      integer index01
      integer j
      integer j_max
      integer j_min

      index01 = 
     &  index_min 
     &  + (         i - i_min ) 
     &  + ( i_max + 1 - i_min ) * ( j - j_min )

      return
      end
      function index012 ( i_min, i, i_max, j_min, j, j_max, k_min, k, 
     &  k_max )

c*********************************************************************72
c
cc INDEX012 indexes a 3D array by columns with zero base.
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
c    27 November 2012
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
c    Output, integer INDEX012, the index of element (I,J,K).
c
      implicit none

      integer i
      integer i_max
      integer i_min
      integer index_min
      parameter ( index_min = 0 )
      integer index012
      integer j
      integer j_max
      integer j_min
      integer k
      integer k_max
      integer k_min

      index012 = 
     &    index_min 
     &  + (         i - i_min )
     &  + ( i_max + 1 - i_min ) * (         j - j_min ) *  
     &  + ( i_max + 1 - i_min ) * ( j_max + 1 - j_min ) * ( k - k_min )

      return
      end
      function index0123 ( i1_min, i1, i1_max, i2_min, i2, i2_max, 
     &  i3_min, i3, i3_max, i4_min, i4, i4_max )

c*********************************************************************72
c
cc INDEX0123 indexes a 4D array by columns, with a zero base.
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
c    27 November 2012
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
c    Output, integer INDEX0123, the index of (I1,I2,I3,I4).
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
      parameter ( index_min = 0 )
      integer index0123

      index0123 = 
     &    index_min 
     &  + (         i1 - i1_min ) 
     &  + ( i1_max + 1 - i1_min ) * (         i2 - i2_min ) 
     &  + ( i1_max + 1 - i1_min ) * ( i2_max + 1 - i2_min ) 
     &  * (         i3 - i3_min ) 
     &  + ( i1_max + 1 - i1_min ) * ( i2_max + 1 - i2_min ) 
     &  * ( i3_max + 1 - i3_min ) * (         i4 - i4_min )

      return
      end
      function index0n ( n, i_min, i, i_max )

c*********************************************************************72
c
cc INDEX0N indexes an N-dimensional array by columns, with zero base.
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
c    27 November 2012
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
c    Output, integer INDEX0N, the index of element I.
c
      implicit none

      integer n

      integer i(n)
      integer i_max(n)
      integer i_min(n)
      integer index_min
      parameter ( index_min = 0 )
      integer index0n
      integer j
      integer value

      value = ( i(n) - i_min(n) )

      do j = n - 1, 1, - 1
        value = value * ( i_max(j) + 1 - i_min(j) ) 
     &    + ( i(j) - i_min(j) )
      end do
      value = value + index_min

      index0n = value

      return
      end
      function index1 ( i_min, i, i_max )

c*********************************************************************72
c
cc INDEX1 indexes a 1D vector using a unit base.
c
c  Discussion:
c
c    Index       Element
c    ---------   --------
c    1           I_MIN
c    INDEX1      I
c   (INDEX_MAX)  I_MAX
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    27 November 2012
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
c    Output, integer INDEX1, the index of element I.
c
      implicit none

      integer i
      integer i_max
      integer i_min
      integer index_min
      parameter ( index_min = 1 )
      integer index1

      index1 = index_min + ( i - i_min )

      return
      end
      function index10 ( i_min, i, i_max, j_min, j, j_max )

c*********************************************************************72
c
cc INDEX10 indexes a 2D array by rows, with a zero base.
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
c    27 November 2012
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
c    Output, integer INDEX10, the index of element (I,J).
c
      implicit none

      integer i
      integer i_max
      integer i_min
      integer index_min
      parameter ( index_min = 0 )
      integer index10
      integer j
      integer j_max
      integer j_min

      index10 = index_min 
     &           +                         ( j - j_min ) 
     &           + ( i - i_min ) * ( j_max + 1 - j_min )

      return
      end
      function index12 ( i_min, i, i_max, j_min, j, j_max )

c*********************************************************************72
c
cc INDEX12 indexes a 2D array by columns, with a unit base.
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
c    27 November 2012
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
c    Output, integer INDEX12, the index of element (I,J).
c
      implicit none

      integer i
      integer i_max
      integer i_min
      integer index_min
      parameter ( index_min = 1 )
      integer index12
      integer j
      integer j_max
      integer j_min

      index12 = 
     &  index_min 
     &  + (         i - i_min ) 
     &  + ( i_max + 1 - i_min ) * ( j - j_min )

      return
      end
      function index123 ( i_min, i, i_max, j_min, j, j_max, k_min, k, 
     &  k_max )

c*********************************************************************72
c
cc INDEX123 indexes a 3D array by columns with unit base.
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
c    27 November 2012
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
c    Output, integer INDEX123, the index of element (I,J,K).
c
      implicit none

      integer i
      integer i_max
      integer i_min
      integer index_min
      parameter ( index_min = 1 )
      integer index123
      integer j
      integer j_max
      integer j_min
      integer k
      integer k_max
      integer k_min

      index123 = 
     &    index_min 
     &  + (         i - i_min ) 
     &  + ( i_max + 1 - i_min ) * (         j - j_min ) *  
     &  + ( i_max + 1 - i_min ) * ( j_max + 1 - j_min ) * ( k - k_min )

      return
      end
      function index1234 ( i1_min, i1, i1_max, i2_min, i2, i2_max, 
     &  i3_min, i3, i3_max, i4_min, i4, i4_max )

c*********************************************************************72
c
cc INDEX1234 indexes a 4D array by columns, with a unit base.
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
c    16 November 2012
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
c    Output, integer INDEX1234, the index of (I1,I2,I3,I4).
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
      parameter ( index_min = 1 )
      integer index1234

      index1234 = 
     &    index_min 
     &  + (         i1 - i1_min ) 
     &  + ( i1_max + 1 - i1_min ) * (         i2 - i2_min ) 
     &  + ( i1_max + 1 - i1_min ) * ( i2_max + 1 - i2_min ) 
     &  * (         i3 - i3_min ) 
     &  + ( i1_max + 1 - i1_min ) * ( i2_max + 1 - i2_min ) 
     &  * ( i3_max + 1 - i3_min ) * (         i4 - i4_min )

      return
      end
      function index1n ( n, i_min, i, i_max )

c*********************************************************************72
c
cc INDEX1N indexes an N-dimensional array by columns, with unit base.
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
c    27 November 2012
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
c    Output, integer INDEX1N, the index of element I.
c
      implicit none

      integer n

      integer i(n)
      integer i_max(n)
      integer i_min(n)
      integer index_min
      parameter ( index_min = 1 )
      integer index1n
      integer j
      integer value

      value = ( i(n) - i_min(n) )

      do j = n - 1, 1, - 1
        value = value * ( i_max(j) + 1 - i_min(j) ) 
     &    + ( i(j) - i_min(j) )
      end do
      value = value + index_min

      index1n = value

      return
      end
      function index21 ( i_min, i, i_max, j_min, j, j_max )

c*********************************************************************72
c
cc INDEX21 indexes a 2D array by rows, with a unit base.
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
c    27 November 2012
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
c    Output, integer INDEX21, the index of element (I,J).
c
      implicit none

      integer i
      integer i_max
      integer i_min
      integer index_min
      parameter ( index_min = 1 )
      integer index21
      integer j
      integer j_max
      integer j_min

      index21 = index_min 
     &           +                         ( j - j_min ) 
     &           + ( i - i_min ) * ( j_max + 1 - j_min )

      return
      end
      function index210 ( i_min, i, i_max, j_min, j, j_max, k_min, 
     &  k, k_max )

c*********************************************************************72
c
cc INDEX210 indexes a 3D array by rows, with zero base.
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
c    27 November 2012
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
c    Output, integer INDEX210, the index of element (I,J,K).
c
      implicit none

      integer i
      integer i_max
      integer i_min
      integer index_min
      parameter ( index_min = 0 )
      integer index210
      integer j
      integer j_max
      integer j_min
      integer k
      integer k_max
      integer k_min

      index210 = 
     &    index_min 
     &  +                                                 ( k - k_min ) 
     &  +                         ( j - j_min ) * ( k_max + 1 - k_min ) 
     &  + ( i - i_min ) * ( j_max + 1 - j_min ) * ( k_max + 1 - k_min )

      return
      end
      function index321 ( i_min, i, i_max, j_min, j, j_max, k_min, k, 
     &  k_max )

c*********************************************************************72
c
cc INDEX321 indexes a 3D array by rows, with zero base.
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
c    27 November 2012
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
c    Output, integer INDEX321, the index of element (I,J,K).
c
      implicit none

      integer i
      integer i_max
      integer i_min
      integer index_min
      parameter ( index_min = 1 )
      integer index321
      integer j
      integer j_max
      integer j_min
      integer k
      integer k_max
      integer k_min

      index321 = 
     &    index_min 
     &  +                                                 ( k - k_min ) 
     &  +                         ( j - j_min ) * ( k_max + 1 - k_min ) 
     &  + ( i - i_min ) * ( j_max + 1 - j_min ) * ( k_max + 1 - k_min )

      return
      end
      function index3210 ( i1_min, i1, i1_max, i2_min, i2, i2_max, 
     &  i3_min, i3, i3_max, i4_min, i4, i4_max )

c*********************************************************************72
c
cc INDEX3210 indexes a 4D array by rows, with zero base.
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
c    27 November 2012
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
c    Output, integer INDEX3210, the index of (I1,I2,I3,I4).
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
      parameter ( index_min = 0 )
      integer index3210

      index3210 = 
     &  index_min 
     &  + ( i4 - i4_min ) 
     &  + ( i3 - i3_min ) 
     &  * ( i4_max + 1 - i4_min ) 
     &  + ( i2 - i2_min ) * ( i3_max + 1 - i3_min ) 
     &  * ( i4_max + 1 - i4_min ) 
     &  + ( i1 - i1_min ) * ( i2_max + 1 - i2_min ) 
     &  * ( i3_max + 1 - i3_min ) * ( i4_max + 1 - i4_min )

      return
      end
      function index4321 ( i1_min, i1, i1_max, i2_min, i2, i2_max, 
     &  i3_min, i3, i3_max, i4_min, i4, i4_max )

c*********************************************************************72
c
cc INDEX4321 indexes a 4D array by rows, with unit base.
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
c    16 November 2012
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
c    Output, integer INDEX4321, the index of (I1,I2,I3,I4).
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
      parameter ( index_min = 1 )
      integer index4321

      index4321 = 
     &  index_min 
     &  + ( i4 - i4_min ) 
     &  + ( i3 - i3_min ) 
     &  * ( i4_max + 1 - i4_min ) 
     &  + ( i2 - i2_min ) * ( i3_max + 1 - i3_min ) 
     &  * ( i4_max + 1 - i4_min ) 
     &  + ( i1 - i1_min ) * ( i2_max + 1 - i2_min ) 
     &  * ( i3_max + 1 - i3_min ) * ( i4_max + 1 - i4_min )

      return
      end
      function indexn0 ( n, i_min, i, i_max )

c*********************************************************************72
c
cc INDEXN0 indexes an N-dimensional array by rows, with zero base.
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
c    16 November 2012
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
c    Output, integer INDEXN0, the index of element I.
c
      implicit none

      integer n

      integer i(n)
      integer i_max(n)
      integer i_min(n)
      integer index_min
      parameter ( index_min = 0 )
      integer indexn0
      integer j
      integer value

      value = ( i(1) - i_min(1) )

      do j = 2, n
        value = value * ( i_max(j) + 1 - i_min(j) ) 
     &    + ( i(j) - i_min(j) )
      end do
      value = value + index_min

      indexn0 = value

      return
      end
      function indexn1 ( n, i_min, i, i_max )

c*********************************************************************72
c
cc INDEXN1 indexes an N-dimensional array by rows, with unit base.
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
c    16 November 2012
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
c    Output, integer INDEXN1, the index of element I.
c
      implicit none

      integer n

      integer i(n)
      integer i_max(n)
      integer i_min(n)
      integer index_min
      parameter ( index_min = 1 )
      integer indexn1
      integer j
      integer value

      value = ( i(1) - i_min(1) )

      do j = 2, n
        value = value * ( i_max(j) + 1 - i_min(j) ) 
     &    + ( i(j) - i_min(j) )
      end do
      value = value + index_min

      indexn1 = value

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
