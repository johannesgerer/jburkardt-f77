      program main

c*********************************************************************72
c
cc SIMPLEX_COORDINATES_PRB tests SIMPLEX_COORDINATES.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    19 September 2010
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n

      call timestamp ( )
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'SIMPLEX_COORDINATES_PRB'
      write ( *, '(a)' ) '  FORTRAN77 version'
      write ( *, '(a)' ) '  Test the SIMPLEX_COORDINATES library.'

      n = 3
      call test01 ( n )
      call test02 ( n )

      n = 4
      call test01 ( n )
      call test02 ( n )
c
c  Terminate.
c
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'SIMPLEX_COORDINATES_PRB'
      write ( *, '(a)' ) '  Normal end of execution.'
      write ( *, '(a)' ) ' '

      call timestamp ( )

      stop
      end
      subroutine test01 ( n )

c*********************************************************************72
c
cc TEST01 calls SIMPLEX_COORDINATES1.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    19 September 2010
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the spatial dimension.
c
      implicit none

      integer n

      integer i
      integer j
      integer k
      double precision r8_factorial
      double precision side
      double precision volume
      double precision volume2
      double precision x(n,n+1)
      double precision xtx(n+1,n+1)

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST01'
      write ( *, '(a)' ) '  Call SIMPLEX_COORDINATES1'

      call simplex_coordinates1 ( n, x )

      call r8mat_transpose_print ( n, n + 1, x, 
     &  '  Simplex vertex coordinates:' )

      side = 0.0D+00
      do i = 1, n
        side = side + ( x(i,1) - x(i,2) )**2
      end do
      side = sqrt ( side )

      call simplex_volume ( n, x, volume )

      volume2 = sqrt ( real ( n + 1, kind = 8 ) ) / r8_factorial ( n ) 
     &  / sqrt ( 2.0D+00**n ) * side**n

      write ( *, '(a)' ) ' '
      write ( *, '(a,g14.6)' ) '  Side length =     ', side
      write ( *, '(a,g14.6)' ) '  Volume =          ', volume
      write ( *, '(a,g14.6)' ) '  Expected volume = ', volume2

      do i = 1, n + 1
        do j = 1, n + 1
          xtx(i,j) = 0.0D+00
          do k = 1, n
            xtx(i,j) = xtx(i,j) + x(k,i) * x(k,j)
          end do
        end do
      end do

      call r8mat_transpose_print ( n + 1, n + 1, xtx, 
     &  '  Dot product matrix:' )

      return
      end
      subroutine test02 ( n )

c*********************************************************************72
c
cc TEST02 calls SIMPLEX_COORDINATES2.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    19 September 2010
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the spatial dimension.
c
      implicit none

      integer n

      integer i
      integer j
      integer k
      double precision r8_factorial
      double precision side
      double precision volume
      double precision volume2
      double precision x(n,n+1)
      double precision xtx(n+1,n+1)

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST02'
      write ( *, '(a)' ) '  Call SIMPLEX_COORDINATES2'

      call simplex_coordinates2 ( n, x )

      call r8mat_transpose_print ( n, n + 1, x, 
     &  '  Simplex vertex coordinates:' )

      side = 0.0D+00
      do i = 1, n
        side = side + ( x(i,1) - x(i,2) )**2
      end do
      side = sqrt ( side )

      call simplex_volume ( n, x, volume )

      volume2 = sqrt ( dble ( n + 1 ) ) / r8_factorial ( n ) 
     &  / sqrt ( 2.0D+00**n ) * side**n

      write ( *, '(a)' ) ' '
      write ( *, '(a,g14.6)' ) '  Side length =     ', side
      write ( *, '(a,g14.6)' ) '  Volume =          ', volume
      write ( *, '(a,g14.6)' ) '  Expected volume = ', volume2

      do i = 1, n + 1
        do j = 1, n + 1
          xtx(i,j) = 0.0D+00
          do k = 1, n
            xtx(i,j) = xtx(i,j) + x(k,i) * x(k,j)
          end do
        end do
      end do

      call r8mat_transpose_print ( n + 1, n + 1, xtx, 
     &  '  Dot product matrix:' )

      return
      end
