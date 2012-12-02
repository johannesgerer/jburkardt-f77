      program main

c*********************************************************************72
c
cc MATRIX_EXPONENTIAL_TEST tests some matrix exponential algorithms.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    17 October 2012
c
c  Author:
c
c    John Burkardt
c
      implicit none

      call timestamp ( )
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'MATRIX_EXPONENTIAL_TEST:'
      write ( *, '(a)' ) '  FORTRAN77 version'
      write ( *, '(a)' ) '  Test the MATRIX_EXPONENTIAL library.'
      write ( *, '(a)' ) '  The R8LIB library is needed.'
      write ( *, '(a)' ) 
     &  '  This test needs the TEST_MATRIX_EXPONENTIAL library.'

      call matrix_exponential_test01 ( )
c
c  Terminate.
c
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'MATRIX_EXPONENTIAL_TEST:'
      write ( *, '(a)' ) '  Normal end of execution.'
      write ( *, '(a)' ) ' '
      call timestamp ( )

      return
      end
      subroutine matrix_exponential_test01 ( )

c*********************************************************************72
c
cc MATRIX_EXPONENTIAL_TEST01 compares matrix exponential algorithms.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    01 December 2011
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n_max
      parameter ( n_max = 10 )

      double precision a(n_max,n_max)
      double precision a_exp(n_max,n_max)
      integer n
      integer test
      integer test_num

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'MATRIX_EXPONENTIAL_TEST01:'
      write ( *, '(a)' ) 
     &  '  EXPM1 is an M-file equivalent to MATLAB''s EXPM'
      write ( *, '(a)' ) '  EXPM2 uses a Taylor series approach'
      write ( *, '(a)' ) '  EXPM3 relies on an eigenvalue calculation.'

      call mexp_test_num ( test_num )

      do test = 1, test_num

        write ( *, '(a)' ) ' '
        write ( *, '(a,i4)' ) '  Test #', test

        call mexp_story ( test )

        call mexp_n ( test, n )

        write ( *, '(a,i4)' ) '  Matrix order N = ', n

        call mexp_a ( test, n, a )

        call r8mat_print ( n, n, a, '  Matrix:' )

        call expm1 ( n, a, a_exp )
        call r8mat_print ( n, n, a_exp, '  EXPM1(A):' )

        call expm2 ( n, a, a_exp )
        call r8mat_print ( n, n, a_exp, '  EXPM2(A):' )

c   call expm3 ( n, a, a_exp )
c   call r8mat_print ( n, n, a_exp, '  EXPM3(A):' )

        call mexp_expa ( test, n, a_exp )
        call r8mat_print ( n, n, a_exp, '  Exact Exponential:' )

      end do

      return
      end
