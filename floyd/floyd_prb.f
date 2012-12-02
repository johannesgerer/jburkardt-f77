      program main

c*********************************************************************72
c
cc MAIN is the main program for FLOYD_PRB.
c
c  Discussion:
c
c    FLOYD_PRB calls a set of problems for FLOYD.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    21 July 2011
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n
      double precision wtime

      call timestamp ( )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'FLOYD_PRB'
      write ( *, '(a)' ) '  FORTRAN77 version'
      write ( *, '(a)' ) '  Test the FLOYD library.'

      call test01 ( )
      call test02 ( )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'FLOYD_TEST03'
      write ( *, '(a)' ) '  Test I4MAT_FLOYD on the MOD(I,J) matrix.'
      write ( *, '(a)' ) '  The work is roughly N^3.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '         N   Time (seconds)  Time/N^3'
      write ( *, '(a)' ) ' '

      n = 1

10    continue

      if ( n .le. 1024 ) then
        call test03 ( n, wtime )
        write ( *, '(2x,i8,2x,g14.6,2x,g14.6)' ) 
     &    n, wtime, 1000000.0D+00 * wtime / dble ( n**3 )
        n = n * 2
        go to 10

      end if
c
c  Terminate.
c
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'FLOYD_PRB'
      write ( *, '(a)' ) '  Normal end of execution.'

      write ( *, '(a)' ) ' '
      call timestamp ( )

      stop
      end
      subroutine test01 ( )

c*********************************************************************72
c
cc TEST01 tests I4MAT_FLOYD.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    21 July 2011
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n
      parameter ( n = 6 )

      integer a(6,6)
      integer huge
      integer i
      integer i4_huge
      integer j

      save a

      data a /
     &   0, -1, -1, -1, -1, -1, 
     &   2,  0, -1, -1, -1,  5, 
     &   5,  7,  0, -1,  2, -1, 
     &  -1,  1,  4,  0, -1,  2, 
     &  -1, -1, -1,  3,  0,  4, 
     &  -1,  8, -1, -1,  3,  0  /

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST01'
      write ( *, '(a)' ) 
     &  '  I4MAT_FLOYO uses Floyd''s algorithm to find the'
      write ( *, '(a)' ) 
     &  '  shortest distance between all pairs of nodes'
      write ( *, '(a)' ) 
     &  '  in a directed graph, starting from the initial array'
      write ( *, '(a)' ) '  of direct node-to-node distances.'

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  In the initial direct distance array, if'
      write ( *, '(a)' ) '    A(I,J) = -1,'
      write ( *, '(a)' ) 
     &  '  this indicates there is NO directed link from'
      write ( *, '(a)' ) 
     &  '  node I to node J.  In that case, the value of'
      write ( *, '(a)' ) '  of A(I,J) is essentially "infinity".'

      call i4mat_print ( n, n, a, '  Initial direct distance array:' )

      huge = i4_huge ( ) / 2

      do i = 1, n
        do j = 1, n
          if ( a(i,j) .eq. - 1 ) then
            a(i,j) = huge
          end if
        end do
      end do

      call i4mat_floyd ( n, a )

      do i = 1, n
        do j = 1, n
          if ( huge .le. a(i,j) ) then
            a(i,j) = - 1
          end if
        end do
      end do

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  In the final shortest distance array, if'
      write ( *, '(a)' ) '    A(I,J) = -1,'
      write ( *, '(a)' ) 
     &  '  this indicates there is NO directed path from'
      write ( *, '(a)' ) '  node I to node J.'

      call i4mat_print ( n, n, a, '  Final shortest distance array:' )

      return
      end
      subroutine test02 ( )

c*********************************************************************72
c
cc TEST02 tests R8MAT_FLOYD.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    21 July 2011
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n
      parameter ( n = 6 )

      double precision a(6,6)
      integer i
      integer j
      double precision r8_huge

      save a

      data a /
     &   0.0D+00, -1.0D+00, -1.0D+00, -1.0D+00, -1.0D+00, -1.0D+00, 
     &   2.0D+00,  0.0D+00, -1.0D+00, -1.0D+00, -1.0D+00,  5.0D+00, 
     &   5.0D+00,  7.0D+00,  0.0D+00, -1.0D+00,  2.0D+00, -1.0D+00, 
     &  -1.0D+00,  1.0D+00,  4.0D+00,  0.0D+00, -1.0D+00,  2.0D+00, 
     &  -1.0D+00, -1.0D+00, -1.0D+00,  3.0D+00,  0.0D+00,  4.0D+00, 
     &  -1.0D+00,  8.0D+00, -1.0D+00, -1.0D+00,  3.0D+00,  0.0D+00  /

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST02'
      write ( *, '(a)' ) 
     &  '  R8MAT_FLOYO uses Floyd''s algorithm to find the'
      write ( *, '(a)' ) 
     &  '  shortest distance between all pairs of nodes'
      write ( *, '(a)' ) 
     &  '  in a directed graph, starting from the initial array'
      write ( *, '(a)' ) '  of direct node-to-node distances.'

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  In the initial direct distance array, if'
      write ( *, '(a)' ) '    A(I,J) = -1,'
      write ( *, '(a)' ) 
     &  '  this indicates there is NO directed link from'
      write ( *, '(a)' ) 
     &  '  node I to node J.  In that case, the value of'
      write ( *, '(a)' ) '  of A(I,J) is essentially "infinity".'

      call r8mat_print ( n, n, a, '  Initial direct distance array:' )

      do i = 1, n
        do j = 1, n
          if ( a(i,j) .eq. - 1.0D+00 ) then
            a(i,j) = r8_huge ( )
          end if
        end do
      end do

      call r8mat_floyd ( n, a )

      do i = 1, n
        do j = 1, n
          if ( r8_huge ( ) .le. a(i,j) ) then
            a(i,j) = - 1.0D+00
          end if
        end do
      end do

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  In the final shortest distance array, if'
      write ( *, '(a)' ) '    A(I,J) = -1,'
      write ( *, '(a)' ) 
     &  '  this indicates there is NO directed path from'
      write ( *, '(a)' ) '  node I to node J.'

      call r8mat_print ( n, n, a, '  Final shortest distance array:' )

      return
      end
      subroutine test03 ( n, wtime )

c*********************************************************************72
c
cc TEST03 tests I4MAT_FLOYD.
c
c  Discussion:
c
c    The matrix size is input by the user.
c
c    The matrix A has the property that
c
c      A(I,J) = 1 if I is divisible by J.
c
c  Example:
c
c    N = 6
c
c    1 0 0 0 0 0
c    1 1 0 0 0 0
c    1 0 1 0 0 0
c    1 1 0 1 0 0
c    1 0 0 0 1 0
c    1 1 1 0 0 1
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    20 July 2011
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the size of the matrix.
c
c    Output, double precision WTIME, the CPU  time required by I4MAT_FLOYD.
c
      implicit none

      integer n

      integer a(n,n)
      integer huge
      integer i
      integer i4_huge
      integer j
      double precision time1
      double precision time2
      double precision wtime

      huge = i4_huge ( ) / 2

      do j = 1, n
        do i = 1, n
          if ( mod ( i, j ) .eq. 0 ) then
            a(i,j) = 1
          else
            a(i,j) = huge
          end if
        end do
      end do

      call cpu_time ( time1 )

      call i4mat_floyd ( n, a )

      call cpu_time ( time2 )

      wtime = time2 - time1

      return
      end
