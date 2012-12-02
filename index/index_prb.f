      program main

c*********************************************************************72
c
cc INDEX_PRB tests INDEX.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    28 November 2012
c
c  Author:
c
c    John Burkardt
c
      implicit none

      call timestamp ( )
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'INDEX_PRB:'
      write ( *, '(a)' ) '  FORTRAN77 version'
      write ( *, '(a)' ) '  Test the INDEX library.'

      call test01 ( )
      call test02 ( )
      call test03 ( )
      call test04 ( )
      call test05 ( )
c
c  Terminate.
c
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'INDEX_PRB:'
      write ( *, '(a)' ) '  Normal end of execution.'
      write ( *, '(a)' ) ' '
      call timestamp ( )

      stop
      end
      subroutine test01 ( )

c*********************************************************************72
c
cc TEST01 tests INDEX0 and INDEX1.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    28 November 2012
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer i
      integer i_max
      integer i_min
      integer index0
      integer index1
      integer index_max
      integer index_min
      integer value

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST01'
      write ( *, '(a)' ) '  INDEX0 indexes a 1D array with zero base,'
      write ( *, '(a)' ) '  INDEX1 indexes a 1D array with  unit base.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '             Min Index   Max'
      write ( *, '(a)' ) ' '

      i_min = 1
      i = 3
      i_max = 5
      write ( *, '(2x,a,2x,i4,2x,i4,2x,i4)' ) 
     &  '1D Index', i_min,     i,     i_max

      value = index0 ( i_min, i, i_max )
      index_min = 0
      index_max = index_min + i_max - i_min
      write ( *, '(2x,a,2x,i4,2x,i4,2x,i4)' ) 
     &  'Index0: ', index_min, value, index_max

      value = index1 ( i_min, i, i_max )
      index_min = 1
      index_max = index_min + i_max - i_min
      write ( *, '(2x,a,2x,i4,2x,i4,2x,i4)' ) 
     &  'Index1: ', index_min, value, index_max

      return
      end
      subroutine test02 ( )

c*********************************************************************72
c
cc TEST02 tests INDEX01, INDEX10, INDEX12 and INDEX21.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    28 November 2012
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer i
      integer i_max
      integer i_min
      integer index01
      integer index10
      integer index12
      integer index21
      integer index_max
      integer index_min
      integer j
      integer j_max
      integer j_min
      integer value

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST02'
      write ( *, '(a)' ) '  For a 2D array,'
      write ( *, '(a)' ) '  INDEX01 column indexes with zero base,'
      write ( *, '(a)' ) '  INDEX10 row indexes with zero base,'
      write ( *, '(a)' ) '  INDEX12 column indexes with unit base,'
      write ( *, '(a)' ) '  INDEX21 row indexes with unit base.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '                Min   Index     Max'
      write ( *, '(a)' ) ' '

      i_min = 1
      i = 3
      i_max = 5
      j_min = 1
      j = 2
      j_max = 4
      write ( *, '(2x,a,2x,i3,i3,2x,i3,i3,2x,i3,i3)' ) '2D Index:', 
     &  i_min, j_min, i, j, i_max, j_max

      value = index01 ( i_min, i, i_max, j_min, j, j_max )
      index_min = 0
      index_max = index_min + ( i_max - i_min + 1 ) 
     &  * ( j_max - j_min + 1 ) - 1
      write ( *, '(2x,a,2x,i6,2x,i6,2x,i6)' ) 
     &  'INDEX01: ', index_min, value, index_max

      value = index10 ( i_min, i, i_max, j_min, j, j_max )
      index_min = 0
      index_max = index_min + ( i_max - i_min + 1 ) 
     &  * ( j_max - j_min + 1 ) - 1
      write ( *, '(2x,a,2x,i6,2x,i6,2x,i6)' ) 
     &  'INDEX10: ', index_min, value, index_max

      value = index12 ( i_min, i, i_max, j_min, j, j_max )
      index_min = 1
      index_max = index_min + ( i_max - i_min + 1 ) 
     &  * ( j_max - j_min + 1 ) - 1
      write ( *, '(2x,a,2x,i6,2x,i6,2x,i6)' ) 
     &  'INDEX12: ', index_min, value, index_max

      value = index21 ( i_min, i, i_max, j_min, j, j_max )
      index_min = 1
      index_max = index_min + ( i_max - i_min + 1 ) 
     &  * ( j_max - j_min + 1 ) - 1
      write ( *, '(2x,a,2x,i6,2x,i6,2x,i6)' ) 
     &  'INDEX21: ', index_min, value, index_max

      return
      end
      subroutine test03 ( )

c*********************************************************************72
c
cc TEST03 tests INDEX012, INDEX123, INDEX210, and INDEX321.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    28 November 2012
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer i
      integer i_max
      integer i_min
      integer index012
      integer index123
      integer index210
      integer index321
      integer index_max
      integer index_min
      integer j
      integer j_max
      integer j_min
      integer k
      integer k_max
      integer k_min
      integer m
      integer value

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST03'
      write ( *, '(a)' ) '  For a 3D array,'
      write ( *, '(a)' ) '  INDEX012 column indexes with zero base,'
      write ( *, '(a)' ) '  INDEX123 column indexes with unit base,'
      write ( *, '(a)' ) '  INDEX210 row indexes with zero base,'
      write ( *, '(a)' ) '  INDEX321 row indexes with unit base.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '                   Min      Index        Max'
      write ( *, '(a)' ) ' '

      i_min = 1
      i = 3
      i_max = 5
      j_min = 1
      j = 2
      j_max = 4
      k_min = 1
      k = 1
      k_max = 3

      m = ( i_max - i_min + 1 ) 
     &  * ( j_max - j_min + 1 ) 
     &  * ( k_max - k_min + 1 )

      write ( *, '(2x,a,2x,i3,i3,i3,2x,i3,i3,i3,2x,i3,i3,i3)' ) 
     &  '3D Index:', i_min, j_min, k_min, i, j, k, i_max, j_max, k_max

      value = index012 ( i_min, i, i_max, j_min, j, j_max, k_min, 
     &  k, k_max )
      index_min = 0
      index_max = index_min + m - 1
      write ( *, '(2x,a,2x,i9,2x,i9,2x,i9)' ) 
     &  'INDEX012:', index_min, value, index_max

      value = index123 ( i_min, i, i_max, j_min, j, j_max, k_min, 
     &  k, k_max )
      index_min = 1
      index_max = index_min + m - 1
      write ( *, '(2x,a,2x,i9,2x,i9,2x,i9)' ) 
     &  'INDEX123:', index_min, value, index_max

      value = index210 ( i_min, i, i_max, j_min, j, j_max, k_min, 
     &  k, k_max )
      index_min = 0
      index_max = index_min + m - 1
      write ( *, '(2x,a,2x,i9,2x,i9,2x,i9)' ) 
     &  'INDEX210:', index_min, value, index_max

      value = index321 ( i_min, i, i_max, j_min, j, j_max, k_min, 
     &  k, k_max )
      index_min = 1
      index_max = index_min + m - 1
      write ( *, '(2x,a,2x,i9,2x,i9,2x,i9)' ) 
     &  'INDEX321:', index_min, value, index_max

      return
      end
      subroutine test04 ( )

c*********************************************************************72
c
cc TEST04 tests INDEX0123, INDEX1234, INDEX3210, and INDEX4321.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    28 November 2012
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer i
      integer i_max
      integer i_min
      integer index0123
      integer index1234
      integer index3210
      integer index4321
      integer index_max
      integer index_min
      integer j
      integer j_max
      integer j_min
      integer k
      integer k_max
      integer k_min
      integer l
      integer l_max
      integer l_min
      integer m
      integer value

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST04'
      write ( *, '(a)' ) '  For a 4D array,'
      write ( *, '(a)' ) '  INDEX0123 column indexes with zero base,'
      write ( *, '(a)' ) '  INDEX1234 column indexes with unit base,'
      write ( *, '(a)' ) '  INDEX3210 row indexes with zero base,'
      write ( *, '(a)' ) '  INDEX4321 row indexes with unit base.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 
     &  '                       Min         Index           Max'
      write ( *, '(a)' ) ' '

      i_min = 1
      i = 3
      i_max = 5
      j_min = 1
      j = 2
      j_max = 4
      k_min = 1
      k = 1
      k_max = 3
      l_min = 1
      l = 2
      l_max = 2

      m = ( i_max - i_min + 1 ) 
     &  * ( j_max - j_min + 1 )
     &  * ( k_max - k_min + 1 ) 
     &  * ( l_max - l_min + 1 )

      write ( *, '(2x,a,2x,i3,i3,i3,i3,2x,i3,i3,i3,i3,2x,i3,i3,i3,i3)' )
     &  '4D Index: ', i_min, j_min, k_min, l_min, i, j, k, l, 
     &  i_max, j_max, k_max, l_max

      value = index0123 ( i_min, i, i_max, j_min, j, j_max, 
     &  k_min, k, k_max, l_min, l, l_max )
      index_min = 0
      index_max = index_min + m - 1
      write ( *, '(2x,a,2x,i12,2x,i12,2x,i12)' ) 
     &  'INDEX0123:', index_min, value, index_max

      value = index1234 ( i_min, i, i_max, j_min, j, j_max, 
     &  k_min, k, k_max, l_min, l, l_max )
      index_min = 1
      index_max = index_min + m - 1
      write ( *, '(2x,a,2x,i12,2x,i12,2x,i12)' ) 
     &  'INDEX1234:', index_min, value, index_max

      value = index3210 ( i_min, i, i_max, j_min, j, j_max, 
     &  k_min, k, k_max, l_min, l, l_max )
      index_min = 0
      index_max = index_min + m - 1
      write ( *, '(2x,a,2x,i12,2x,i12,2x,i12)' ) 
     &  'INDEX3210:', index_min, value, index_max

      value = index4321 ( i_min, i, i_max, j_min, j, j_max, 
     &  k_min, k, k_max, l_min, l, l_max )
      index_min = 1
      index_max = index_min + m - 1
      write ( *, '(2x,a,2x,i12,2x,i12,2x,i12)' ) 
     &  'INDEX4321:', index_min, value, index_max

      return
      end
      subroutine test05 ( )

c*********************************************************************72
c
cc TEST05 tests INDEX0N, INDEX1N, INDEXN0 and INDEXN1.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    28 November 2012
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n
      parameter ( n = 4 )

      integer i(n)
      integer i_max(n)
      integer i_min(n)
      integer index0n
      integer index1n
      integer indexn0
      integer indexn1
      integer index_max
      integer index_min
      integer j
      integer m
      integer value

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST05'
      write ( *, '(a)' ) '  For an N-dimensional array,'
      write ( *, '(a)' ) '  INDEX0N column indexes with zero base,'
      write ( *, '(a)' ) '  INDEX1N column indexes with unit base,'
      write ( *, '(a)' ) '  INDEXN0 row indexes with zero base,'
      write ( *, '(a)' ) '  INDEXN1 row indexes with unit base.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 
     &  '                       Min         Index           Max'

      i_min(1) = 1
      i(1) = 3
      i_max(1) = 5
      i_min(2) = 1
      i(2) = 2
      i_max(2) = 4
      i_min(3) = 1
      i(3) = 1
      i_max(3) = 3
      i_min(4) = 1
      i(4) = 2
      i_max(4) = 2

      m = 1
      do j = 1, n
        m = m * ( i_max(j) - i_min(j) + 1 )
      end do

      write ( *, '(2x,a,2x,i3,i3,i3,i3,2x,i3,i3,i3,i3,2x,i3,i3,i3,i3)' ) 
     &  'ND Index: ', i_min(1:n), i(1:n), i_max(1:n)
      value = index0n ( n, i_min, i, i_max )
      index_min = 0
      index_max = index_min + m - 1
      write ( *, '(2x,a,2x,i12,2x,i12,2x,i12)' ) 
     &  'INDEX0N:  ', index_min, value, index_max

      value = index1n ( n, i_min, i, i_max )
      index_min = 1
      index_max = index_min + m - 1
      write ( *, '(2x,a,2x,i12,2x,i12,2x,i12)' ) 
     &  'INDEX1N:  ', index_min, value, index_max

      value = indexn0 ( n, i_min, i, i_max )
      index_min = 0
      index_max = index_min + m - 1
      write ( *, '(2x,a,2x,i12,2x,i12,2x,i12)' ) 
     &  'INDEXN0:  ', index_min, value, index_max

      value = indexn1 ( n, i_min, i, i_max )
      index_min = 1
      index_max = index_min + m - 1
      write ( *, '(2x,a,2x,i12,2x,i12,2x,i12)' ) 
     &  'INDEXN1:  ', index_min, value, index_max

      return
      end


