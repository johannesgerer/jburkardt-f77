      program main

c*********************************************************************72
c
cc CONDITION_PRB tests the CONDITION library.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    22 August 2012
c
c  Author:
c
c    John Burkardt
c
      implicit none

      call timestamp ( )
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'CONDITION_PRB:'
      write ( *, '(a)' ) '  FORTRAN77 version'
      write ( *, '(a)' ) '  Test the CONDITION library.'

      call test01 ( )
      call test02 ( )
      call test03 ( )
c
c  Terminate.
c
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'CONDITION_PRB:'
      write ( *, '(a)' ) '  Normal end of execution.'
      write ( *, '(a)' ) ' '
      call timestamp ( )

      return
      end
      subroutine test01 ( )

c*********************************************************************72
c
cc TEST01 tests CONDITION_LINPACK.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    22 August 2012
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n_max
      parameter ( n_max = 5 )

      double precision a(n_max,n_max)
      double precision a_inverse(n_max,n_max)
      double precision a_inverse_norm_l1
      double precision a_norm_l1
      double precision alpha
      double precision beta
      double precision cond
      double precision cond_l1
      integer i
      integer i1
      integer info
      integer j
      integer n
      character * ( 80 ) name
      integer pivot(n_max)
      double precision r8mat_norm_l1
      integer seed
      double precision z(n_max)

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST01'
      write ( *, '(a)' ) '  For a matrix in general storage,'
      write ( *, '(a)' ) 
     &  '  CONDITION_LINPACK estimates the L1 condition number.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 
     &  '  Matrix               Order   Condition         Linpack'
      write ( *, '(a)' ) ' '
c
c  Combinatorial matrix.
c
      name = 'Combinatorial'
      n = 4
      alpha = 2.0D+00
      beta = 3.0D+00
      call combin ( alpha, beta, n, a )
      call combin_inverse ( alpha, beta, n, a_inverse )
      a_norm_l1         = r8mat_norm_l1 ( n, n, a )
      a_inverse_norm_l1 = r8mat_norm_l1 ( n, n, a_inverse )
      cond_l1           = a_norm_l1 * a_inverse_norm_l1
      call condition_linpack ( n, a, pivot, cond, z )
      write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6)' ) 
     &  name, n, cond_l1, cond
c
c  CONEX1
c
      name = 'CONEX1'
      n = 4
      alpha = 100.0D+00
      call conex1 ( alpha, a )
      call conex1_inverse ( alpha, a_inverse )
      a_norm_l1         = r8mat_norm_l1 ( n, n, a )
      a_inverse_norm_l1 = r8mat_norm_l1 ( n, n, a_inverse )
      cond_l1           = a_norm_l1 * a_inverse_norm_l1
      call condition_linpack ( n, a, pivot, cond, z )
      write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6)' ) 
     &  name, n, cond_l1, cond
c
c  CONEX2
c
      name = 'CONEX2'
      n = 3
      alpha = 100.0D+00
      call conex2 ( alpha, a )
      call conex2_inverse ( alpha, a_inverse )
      a_norm_l1         = r8mat_norm_l1 ( n, n, a )
      a_inverse_norm_l1 = r8mat_norm_l1 ( n, n, a_inverse )
      cond_l1           = a_norm_l1 * a_inverse_norm_l1
      call condition_linpack ( n, a, pivot, cond, z )
      write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6)' ) 
     &  name, n, cond_l1, cond
c
c  CONEX3
c
      name = 'CONEX3'
      n = 5
      call conex3 ( n, a )
      call conex3_inverse ( n, a_inverse )
      a_norm_l1         = r8mat_norm_l1 ( n, n, a )
      a_inverse_norm_l1 = r8mat_norm_l1 ( n, n, a_inverse )
      cond_l1           = a_norm_l1 * a_inverse_norm_l1
      call condition_linpack ( n, a, pivot, cond, z )
      write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6)' ) 
     &  name, n, cond_l1, cond
c
c  CONEX4
c
      name = 'CONEX4'
      n = 4
      call conex4 ( a )
      call conex4_inverse ( a_inverse )
      a_norm_l1         = r8mat_norm_l1 ( n, n, a )
      a_inverse_norm_l1 = r8mat_norm_l1 ( n, n, a_inverse )
      cond_l1           = a_norm_l1 * a_inverse_norm_l1
      call condition_linpack ( n, a, pivot, cond, z )
      write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6)' ) 
     &  name, n, cond_l1, cond
c
c  KAHAN
c
      name = 'KAHAN'
      n = 4
      alpha = 0.25D+00
      call kahan ( alpha, n, n, a )
      call kahan_inverse ( alpha, n, a_inverse )
      a_norm_l1         = r8mat_norm_l1 ( n, n, a )
      a_inverse_norm_l1 = r8mat_norm_l1 ( n, n, a_inverse )
      cond_l1           = a_norm_l1 * a_inverse_norm_l1
      call condition_linpack ( n, a, pivot, cond, z )
      write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6)' ) 
     &  name, n, cond_l1, cond
c
c  Random
c
      seed = 123456789

      do i = 1, 5
        name = 'RANDOM'
        n = 4
        call r8mat_uniform_01 ( n, n, seed, a )
        call r8mat_copy ( n, n, a, a_inverse )
        call r8ge_fa ( n, a_inverse, pivot, info )
        call r8ge_inverse ( n, a_inverse, pivot )
        a_norm_l1         = r8mat_norm_l1 ( n, n, a )
        a_inverse_norm_l1 = r8mat_norm_l1 ( n, n, a_inverse )
        cond_l1           = a_norm_l1 * a_inverse_norm_l1
        call condition_linpack ( n, a, pivot, cond, z )
        write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6)' ) 
     &    name, n, cond_l1, cond
      end do

      return
      end
      subroutine test02 ( )

c*********************************************************************72
c
cc TEST02 tests CONDITION_SAMPLE1.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    22 August 2012
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n_max
      parameter ( n_max = 5 )

      double precision a(n_max,n_max)
      double precision a_inverse(n_max,n_max)
      double precision a_inverse_norm_l1
      double precision a_norm_l1
      double precision alpha
      double precision beta
      double precision cond
      double precision cond_l1
      integer i
      integer i1
      integer info
      integer j
      integer j1
      integer m
      integer m_test(3)
      integer n
      character * ( 80 ) name
      integer pivot(n_max)
      double precision r8mat_norm_l1
      integer seed

      save m_test

      data m_test / 10, 1000, 100000 /

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST02'
      write ( *, '(a)' ) '  For a matrix in general storage,'
      write ( *, '(a)' ) 
     &  '  CONDITION_SAMPLE1 estimates the L1 condition'
     &  // ' number using sampling.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 
     &  '  Matrix                 Samples Order   '
     &  // 'Condition        Estimate'
c
c  Combinatorial matrix.
c
      name = 'Combinatorial'
      n = 4
      alpha = 2.0D+00
      beta = 3.0D+00
      call combin ( alpha, beta, n, a )
      call combin_inverse ( alpha, beta, n, a_inverse )
      a_norm_l1         = r8mat_norm_l1 ( n, n, a )
      a_inverse_norm_l1 = r8mat_norm_l1 ( n, n, a_inverse )
      cond_l1           = a_norm_l1 * a_inverse_norm_l1
      write ( *, '(a)' ) ' '
      do i = 1, 3
        m = m_test(i)
        call condition_sample1 ( n, a, m, cond )
        write ( *, '(2x,a20,2x,i8,2x,i4,2x,g14.6,2x,g14.6)' ) 
     &    name, m, n, cond_l1, cond
      end do
c
c  CONEX1
c
      name = 'CONEX1'
      n = 4
      alpha = 100.0D+00
      call conex1 ( alpha, a )
      call conex1_inverse ( alpha, a_inverse )
      a_norm_l1         = r8mat_norm_l1 ( n, n, a )
      a_inverse_norm_l1 = r8mat_norm_l1 ( n, n, a_inverse )
      cond_l1           = a_norm_l1 * a_inverse_norm_l1
      write ( *, '(a)' ) ' '
      do i = 1, 3
        m = m_test(i)
        call condition_sample1 ( n, a, m, cond )
        write ( *, '(2x,a20,2x,i8,2x,i4,2x,g14.6,2x,g14.6)' ) 
     &    name, m, n, cond_l1, cond
      end do
c
c  CONEX2
c
      name = 'CONEX2'
      n = 3
      alpha = 100.0D+00
      call conex2 ( alpha, a )
      call conex2_inverse ( alpha, a_inverse )
      a_norm_l1         = r8mat_norm_l1 ( n, n, a )
      a_inverse_norm_l1 = r8mat_norm_l1 ( n, n, a_inverse )
      cond_l1           = a_norm_l1 * a_inverse_norm_l1
      write ( *, '(a)' ) ' '
      do i = 1, 3
        m = m_test(i)
        call condition_sample1 ( n, a, m, cond )
        write ( *, '(2x,a20,2x,i8,2x,i4,2x,g14.6,2x,g14.6)' ) 
     &    name, m, n, cond_l1, cond
      end do
c
c  CONEX3
c
      name = 'CONEX3'
      n = 5
      call conex3 ( n, a )
      call conex3_inverse ( n, a_inverse )
      a_norm_l1         = r8mat_norm_l1 ( n, n, a )
      a_inverse_norm_l1 = r8mat_norm_l1 ( n, n, a_inverse )
      cond_l1           = a_norm_l1 * a_inverse_norm_l1
      write ( *, '(a)' ) ' '
      do i = 1, 3
        m = m_test(i)
        call condition_sample1 ( n, a, m, cond )
        write ( *, '(2x,a20,2x,i8,2x,i4,2x,g14.6,2x,g14.6)' ) 
     &    name, m, n, cond_l1, cond
      end do
c
c  CONEX4
c
      name = 'CONEX4'
      n = 4
      call conex4 ( a )
      call conex4_inverse ( a_inverse )
      a_norm_l1         = r8mat_norm_l1 ( n, n, a )
      a_inverse_norm_l1 = r8mat_norm_l1 ( n, n, a_inverse )
      cond_l1           = a_norm_l1 * a_inverse_norm_l1
      write ( *, '(a)' ) ' '
      do i = 1, 3
        m = m_test(i)
        call condition_sample1 ( n, a, m, cond )
        write ( *, '(2x,a20,2x,i8,2x,i4,2x,g14.6,2x,g14.6)' ) 
     &    name, m, n, cond_l1, cond
      end do
c
c  KAHAN
c
      name = 'KAHAN'
      n = 4
      alpha = 0.25D+00
      call kahan ( alpha, n, n, a )
      call kahan_inverse ( alpha, n, a_inverse )
      a_norm_l1         = r8mat_norm_l1 ( n, n, a )
      a_inverse_norm_l1 = r8mat_norm_l1 ( n, n, a_inverse )
      cond_l1           = a_norm_l1 * a_inverse_norm_l1
      write ( *, '(a)' ) ' '
      do i = 1, 3
        m = m_test(i)
        call condition_sample1 ( n, a, m, cond )
        write ( *, '(2x,a20,2x,i8,2x,i4,2x,g14.6,2x,g14.6)' ) 
     &    name, m, n, cond_l1, cond
      end do
c
c  Random
c
      seed = 123456789

      do j = 1, 5
        name = 'RANDOM'
        n = 4
        call r8mat_uniform_01 ( n, n, seed, a )
        call r8mat_copy ( n, n, a, a_inverse )
        call r8ge_fa ( n, a_inverse, pivot, info )
        call r8ge_inverse ( n, a_inverse, pivot )
        a_norm_l1         = r8mat_norm_l1 ( n, n, a )
        a_inverse_norm_l1 = r8mat_norm_l1 ( n, n, a_inverse )
        cond_l1           = a_norm_l1 * a_inverse_norm_l1
        write ( *, '(a)' ) ' '
        do i = 1, 3
          m = m_test(i)
          call condition_sample1 ( n, a, m, cond )
          write ( *, '(2x,a20,2x,i8,2x,i4,2x,g14.6,2x,g14.6)' ) 
     &      name, m, n, cond_l1, cond
        end do
      end do

      return
      end
      subroutine test03 ( )

c*********************************************************************72
c
cc TEST03 tests CONDITION_HAGER.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    22 August 2012
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n_max
      parameter ( n_max = 5 )

      double precision a(n_max,n_max)
      double precision a_inverse(n_max,n_max)
      double precision a_inverse_norm_l1
      double precision a_norm_l1
      double precision alpha
      double precision beta
      double precision cond
      double precision cond_l1
      integer i
      integer i1
      integer info
      integer j
      integer n
      character * ( 80 ) name
      integer pivot(n_max)
      double precision r8mat_norm_l1
      integer seed

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST03'
      write ( *, '(a)' ) '  For a matrix in general storage,'
      write ( *, '(a)' ) 
     &  '  CONDITION_HAGER estimates the L1 condition number.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 
     &  '  Matrix               Order   Condition         Hager'
      write ( *, '(a)' ) ' '
c
c  Combinatorial matrix.
c
      name = 'Combinatorial'
      n = 4
      alpha = 2.0D+00
      beta = 3.0D+00
      call combin ( alpha, beta, n, a )
      call combin_inverse ( alpha, beta, n, a_inverse )
      a_norm_l1         = r8mat_norm_l1 ( n, n, a )
      a_inverse_norm_l1 = r8mat_norm_l1 ( n, n, a_inverse )
      cond_l1           = a_norm_l1 * a_inverse_norm_l1
      call condition_hager ( n, a, cond )
      write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6)' ) 
     &  name, n, cond_l1, cond
c
c  CONEX1
c
      name = 'CONEX1'
      n = 4
      alpha = 100.0D+00
      call conex1 ( alpha, a )
      call conex1_inverse ( alpha, a_inverse )
      a_norm_l1         = r8mat_norm_l1 ( n, n, a )
      a_inverse_norm_l1 = r8mat_norm_l1 ( n, n, a_inverse )
      cond_l1           = a_norm_l1 * a_inverse_norm_l1
      call condition_hager ( n, a, cond )
      write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6)' ) 
     &  name, n, cond_l1, cond
c
c  CONEX2
c
      name = 'CONEX2'
      n = 3
      alpha = 100.0D+00
      call conex2 ( alpha, a )
      call conex2_inverse ( alpha, a_inverse )
      a_norm_l1         = r8mat_norm_l1 ( n, n, a )
      a_inverse_norm_l1 = r8mat_norm_l1 ( n, n, a_inverse )
      cond_l1           = a_norm_l1 * a_inverse_norm_l1
      call condition_hager ( n, a, cond )
      write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6)' ) 
     &  name, n, cond_l1, cond
c
c  CONEX3
c
      name = 'CONEX3'
      n = 5
      call conex3 ( n, a )
      call conex3_inverse ( n, a_inverse )
      a_norm_l1         = r8mat_norm_l1 ( n, n, a )
      a_inverse_norm_l1 = r8mat_norm_l1 ( n, n, a_inverse )
      cond_l1           = a_norm_l1 * a_inverse_norm_l1
      call condition_hager ( n, a, cond )
      write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6)' ) 
     &  name, n, cond_l1, cond
c
c  CONEX4
c
      name = 'CONEX4'
      n = 4
      call conex4 ( a )
      call conex4_inverse ( a_inverse )
      a_norm_l1         = r8mat_norm_l1 ( n, n, a )
      a_inverse_norm_l1 = r8mat_norm_l1 ( n, n, a_inverse )
      cond_l1           = a_norm_l1 * a_inverse_norm_l1
      call condition_hager ( n, a, cond )
      write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6)' ) 
     &  name, n, cond_l1, cond
c
c  KAHAN
c
      name = 'KAHAN'
      n = 4
      alpha = 0.25D+00
      call kahan ( alpha, n, n, a )
      call kahan_inverse ( alpha, n, a_inverse )
      a_norm_l1         = r8mat_norm_l1 ( n, n, a )
      a_inverse_norm_l1 = r8mat_norm_l1 ( n, n, a_inverse )
      cond_l1           = a_norm_l1 * a_inverse_norm_l1
      call condition_hager ( n, a, cond )
      write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6)' ) 
     &  name, n, cond_l1, cond
c
c  Random
c
      seed = 123456789

      do i = 1, 5
        name = 'RANDOM'
        n = 4
        call r8mat_uniform_01 ( n, n, seed, a )
        call r8mat_copy ( n, n, a, a_inverse )
        call r8ge_fa ( n, a_inverse, pivot, info )
        call r8ge_inverse ( n, a_inverse, pivot )
        a_norm_l1         = r8mat_norm_l1 ( n, n, a )
        a_inverse_norm_l1 = r8mat_norm_l1 ( n, n, a_inverse )
        cond_l1           = a_norm_l1 * a_inverse_norm_l1
        call condition_hager ( n, a, cond )
        write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6)' ) 
     &    name, n, cond_l1, cond
      end do

      return
      end
