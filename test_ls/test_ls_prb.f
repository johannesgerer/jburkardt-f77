      program main

c*********************************************************************72
c
cc MAIN is the main program for TEST_LS_PRB.
c
c  Discussion:
c
c    TEST_LS_PRB tests the TEST_LS library.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    07 September 2012
c
c  Author:
c
c    John Burkardt
c
      implicit none

      call timestamp ( )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST_LS_PRB'
      write ( *, '(a)' ) '  FORTRAN77 version'
      write ( *, '(a)' ) '  Test the TEST_LS library.'

      call test01 ( )
c
c  Terminate.
c
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST_LS_PRB'
      write ( *, '(a)' ) '  Normal end of execution.'

      write ( *, '(a)' ) ' '
      call timestamp ( )

      stop
      end
      subroutine test01 ( )

c*********************************************************************72
c
cc TEST01 summarizes the test data.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    07 September 2012
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer m_max
      parameter ( m_max = 10 )
      integer n_max
      parameter ( n_max = 10 )

      double precision a(m_max,n_max)
      double precision b(m_max)
      double precision b_norm
      integer i
      integer m
      integer n
      integer prob
      integer prob_num
      double precision r(m_max)
      double precision r_norm
      double precision r8vec_norm
      double precision x(n_max)
      double precision x_norm

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST01'
      write ( *, '(a)' ) 
     &  '  Get each least squares test, compute the maximum residual.'
      write ( *, '(a)' ) 
     &  '  The L2 norm of the residual MUST be no greater than the'
      write ( *, '(a)' ) 
     &  '  L2 norm of the right hand side, else 0 is a better solution.'

      call p00_prob_num ( prob_num )

      write ( *, '(a)' ) ' '
      write ( *, '(a,i4)' ) '  Number of problems = ', prob_num
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 
     &  '  Index     M     N     ||B||         ||X||         ||R||'
      write ( *, '(a)' ) ' '

      do prob = 1, prob_num

        call p00_m ( prob, m )
        call p00_n ( prob, n )

        call p00_a ( prob, m, n, a )
        call p00_b ( prob, m, b )
        call p00_x ( prob, n, x )

        call r8mat_mv ( m, n, a, x, r )
        do i = 1, m
          r(i) = r(i) - b(i)
        end do

        b_norm = r8vec_norm ( m, b )
        x_norm = r8vec_norm ( n, x )
        r_norm = r8vec_norm ( m, r )

        write ( *, '(2x,i5,2x,i4,2x,i4,2x,g12.4,2x,g12.4,2x,g12.4)' ) 
     &  prob, m, n, b_norm, x_norm, r_norm

      end do

      return
      end
