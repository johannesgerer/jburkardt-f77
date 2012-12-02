      program main

c*********************************************************************72
c
cc TEST_MATRIX_EXPONENTIAL_TEST tests some matrix exponential algorithms.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    29 November 2011
c
c  Author:
c
c    John Burkardt
c
      implicit none

      call timestamp ( );
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST_MATRIX_EXPONENTIAL_TEST:'
      write ( *, '(a)' ) '  FORTRAN77 version'
      write ( *, '(a)' ) '  Test the TEST_MATRIX_EXPONENTIAL library.'
      write ( *, '(a)' ) '  The R8LIB library is needed.'

      call test_matrix_exponential_test01 ( );
c
c  Terminate.
c
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST_MATRIX_EXPONENTIAL_TEST:'
      write ( *, '(a)' ) '  Normal end of execution.'
      write ( *, '(a)' ) ' '
      call timestamp ( );

      return
      end
      subroutine test_matrix_exponential_test01 ( )

c*********************************************************************72
c
cc TEST_MATRIX_EXPONENTIAL_TEST01 retrieves the test data.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    29 November 2011
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n_max
      parameter ( n_max = 10 )

      double precision a(n_max,n_max)
      double precision expa(n_max,n_max)
      integer n
      integer test
      integer test_num

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST_MATRIX_EXPONENTIAL_TEST01:'
      write ( *, '(a)' ) 
     &  '  Retrieve the data for each matrix exponential test.'

      call mexp_test_num ( test_num )

      do test = 1, test_num

        write ( *, '(a)' ) ' '
        write ( *, '(a,i4)' ) '  Test #', test

        call mexp_n ( test, n )

        call mexp_story ( test );

        write ( *, '(a)' ) ' '
        write ( *, '(a,i4)' ) '  Matrix order N = ', n

        call mexp_a ( test, n, a )

        call r8mat_print ( n, n, a, '  Matrix A:' );

        call mexp_expa ( test, n, expa )
        call r8mat_print ( n, n, expa, '  Exact Exponential exp(A):' )

      end do

      return
      end
