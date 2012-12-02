      program main

c*********************************************************************72
c
cc MAIN is the main program for FEM2D_PACK_PRB.
c
c  Discussion:
c
c    FEM2D_PACK_PRB calls the various FEM2D_PACK tests.
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
      implicit none

      call timestamp ( )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'FEM2D_PACK_PRB'
      write ( *, '(a)' ) '  FORTRAN77 version'
      write ( *, '(a)' ) '  Test the FEM2D_PACK library.'

      call test03 ( )
c
c  Terminate.
c
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'FEM2D_PACK_PRB'
      write ( *, '(a)' ) '  Normal end of execution.'

      write ( *, '(a)' ) ' '
      call timestamp ( )

      stop
      end
      subroutine test03 ( )

c*********************************************************************72
c
cc TEST03 tests BASIS_11_**_TEST.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    29 June 2009
c
c  Author:
c
c    John Burkardt
c
      implicit none

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST03'
      write ( *, '(a)' ) '  Test the computation of ONE basis function'
      write ( *, '(a)' ) '  at ONE point in a given element:'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  BASIS_11_Q4_TEST : Q4 element.'
      write ( *, '(a)' ) '  BASIS_11_T3_TEST : T3 element.'
      write ( *, '(a)' ) '  BASIS_11_T4_TEST : T4 element.'
      write ( *, '(a)' ) '  BASIS_11_T6_TEST : T6 element.'

      call basis_11_q4_test ( )

      call basis_11_t3_test ( )

      call basis_11_t4_test ( )

      call basis_11_t6_test ( )

      return
      end
