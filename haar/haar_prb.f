      program main

c*********************************************************************72
c
cc MAIN is the main program for HAAR_PRB.
c
c  Discussion:
c
c    HAAR_PRB tests the HAAR library.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    17 March 2011
c
c  Author:
c
c    John Burkardt
c
      implicit none

      call timestamp ( )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'HAAR_PRB'
      write ( *, '(a)' ) '  FORTRAN77 version'
      write ( *, '(a)' ) '  Test the HAAR library.'

      call test01 ( )
      call test02 ( )
c
c  Terminate.
c
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'HAAR_PRB'
      write ( *, '(a)' ) '  Normal end of execution.'

      write ( *, '(a)' ) ' '
      call timestamp ( )

      return
      end
      subroutine test01 ( )

c*********************************************************************72
c
cc TEST01 tests HAAR_1D.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    14 March 2011
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n_max
      parameter ( n_max = 16 )

      double precision a_first
      double precision a_last
      integer i
      integer n
      integer seed
      double precision u(n_max)
      double precision v(n_max)
      double precision w(n_max)

      write ( *, '(a)' ) ''
      write ( *, '(a)' ) 'TEST01'
      write ( *, '(a)' ) 
     &  '  HAAR_1D computes the Haar transform of a vector.'
c
c  Random data.
c
      n = 16
      seed = 123456789
      call r8vec_uniform_01 ( n, seed, u )
      call r8vec_copy ( n, u, v )

      call haar_1d ( n, v )

      call r8vec_copy ( n, v, w )
      call haar_1d_inverse ( n, w )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '   i      U(i)        H(U)(i)  Hinv(H(U))(i)'
      write ( *, '(a)' ) ' '
      do i = 1, n
        write ( *, '(2x,i2,2x,f10.4,2x,f10.4,2x,f10.4)' ) 
     &    i, u(i), v(i), w(i)
      end do
c
c  Constant signal.
c
      n = 8
      call r8vec_ones ( n, u )
      call r8vec_copy ( n, u, v )

      call haar_1d ( n, v )

      call r8vec_copy ( n, v, w )
      call haar_1d_inverse ( n, w )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '   i      U(i)        H(U)(i)  Hinv(H(U))(i)'
      write ( *, '(a)' ) ' '
      do i = 1, n
        write ( *, '(2x,i2,2x,f10.4,2x,f10.4,2x,f10.4)' ) 
     &    i, u(i), v(i), w(i)
      end do
c
c  Linear signal.
c
      n = 16
      a_first = 1.0D+00
      a_last = dble ( n )
      call r8vec_linspace ( n, a_first, a_last, u )
      call r8vec_copy ( n, u, v )

      call haar_1d ( n, v )

      call r8vec_copy ( n, v, w )
      call haar_1d_inverse ( n, w )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '   i      U(i)        H(U)(i)  Hinv(H(U))(i)'
      write ( *, '(a)' ) ' '
      do i = 1, n
        write ( *, '(2x,i2,2x,f10.4,2x,f10.4,2x,f10.4)' ) 
     &    i, u(i), v(i), w(i)
      end do
c
c  Quadratic data.
c
      n = 8
      u(1) = 25.0D+00
      u(2) = 16.0D+00
      u(3) = 9.0D+00
      u(4) = 4.0D+00
      u(5) = 1.0D+00
      u(6) = 0.0D+00
      u(7) = 1.0D+00
      u(8) = 4.0D+00
      call r8vec_copy ( n, u, v )

      call haar_1d ( n, v )

      call r8vec_copy ( n, v, w )
      call haar_1d_inverse ( n, w )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '   i      U(i)        H(U)(i)  Hinv(H(U))(i)'
      write ( *, '(a)' ) ' '
      do i = 1, n
        write ( *, '(2x,i2,2x,f10.4,2x,f10.4,2x,f10.4)' ) 
     &    i, u(i), v(i), w(i)
      end do

      return
      end
      subroutine test02 ( )

c*********************************************************************72
c
cc TEST02 tests HAAR_2D and HAAR_2D_INVERSE.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    17 March 2011
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer m
      parameter ( m = 16 )
      integer n
      parameter ( n = 4 )

      integer i
      integer j
      integer seed
      real ( kind = 8 ) u(m,n)
      real ( kind = 8 ) v(m,n)
      real ( kind = 8 ) w(m,n)

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST02'
      write ( *, '(a)' ) 
     &  '  HAAR_2D computes the Haar transform of an array.'
      write ( *, '(a)' ) '  HAAR_2D_INVERSE inverts the transform.'
!
!  Demonstrate successful inversion.
!
      seed = 123456789
      call r8mat_uniform_01 ( m, n, seed, u )

      call r8mat_print ( m, n, u, '  Input array U:' )

      do j = 1, n
        do i = 1, m
          v(i,j) = u(i,j)
        end do
      end do

      call haar_2d ( m, n, v )

      call r8mat_print ( m, n, v, '  Transformed array V:' )

      do j = 1, n
        do i = 1, m
          w(i,j) = v(i,j)
        end do
      end do

      call haar_2d_inverse ( m, n, w )

      call r8mat_print ( m, n, w, '  Recovered array W:' )

      return
      end
