      program main

c*********************************************************************72
c
cc MAIN is the main program for TEST_MAT_PRB.
c
c  Discussion:
c
c    TEST_MAT_PRB calls the TEST_MAT test routines.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    07 April 2012
c
c  Author:
c
c    John Burkardt
c
      implicit none

      call timestamp ( )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST_MAT_PRB'
      write ( *, '(a)' ) '  FORTRAN77 version'
      write ( *, '(a)' ) '  Test the TEST_MAT library.'

      call test_cond ( )
      call test_determinant ( )
      call test_eigen ( )
      call test_inverse ( )
      call test_null ( )
      call test_plu ( )
      call test_solution ( )
c
c  Terminate.
c
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST_MAT_PRB'
      write ( *, '(a)' ) '  Normal end of execution.'

      write ( *, '(a)' ) ' '
      call timestamp ( )

      stop
      end
      subroutine test_cond ( )

c*********************************************************************72
c
cc TEST_COND tests the condition number computations.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    10 April 2012
c
c  Author:
c
c    John Burkardt
c
      implicit none

      double precision alpha
      double precision beta
      double precision cond
      integer n
      double precision r8_uniform_01

      integer seed
      character * ( 20 ) title

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST_COND'
      write ( *, '(a)' ) '  Compute the condition number of an example'
      write ( *, '(a)' ) '  of each test matrix'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Matrix title             N      COND'
      write ( *, '(a)' ) ' '
c
c  AEGERTER matrix.
c
      title = 'AEGERTER'
      n = 5
      call aegerter_condition ( n, cond )
      write ( *, '(2x,a20,2x,i4,2x,g14.6)' ) title, n, cond
c
c  BAB matrix.
c
      title = 'BAB'
      seed = 123456789
      alpha = anint ( 50.0D+00 * r8_uniform_01 ( seed ) - 25.0D+00 ) 
     &  / 5.0D+00
      beta = anint ( 50.0D+00 * r8_uniform_01 ( seed ) - 25.0D+00 ) 
     &  / 5.0D+00
      call bab_condition ( n, alpha, beta, cond )
      write ( *, '(2x,a20,2x,i4,2x,g14.6)' ) title, n, cond
c
c  BODEWIG matrix.
c
      title = 'BODEWIG'
      n = 4
      call bodewig_condition ( cond )
      write ( *, '(2x,a20,2x,i4,2x,g14.6)' ) title, n, cond
c
c  COMBIN matrix.
c
      title = 'COMBIN'
      n = 3
      seed = 123456789
      alpha = r8_uniform_01 ( seed )
      alpha = anint ( 50.0D+00 * alpha ) / 5.0D+00
      beta = r8_uniform_01 ( seed )
      beta = anint ( 50.0D+00 * beta ) / 5.0D+00
      call combin_condition ( alpha, beta, n, cond )
      write ( *, '(2x,a20,2x,i4,2x,g14.6)' ) title, n, cond
c
c  CONEX3 matrix.
c
      title = 'CONEX3'
      n = 5
      call conex3_condition ( n, cond )
      write ( *, '(2x,a20,2x,i4,2x,g14.6)' ) title, n, cond
c
c  RUTIS5 matrix.
c
      title = 'RUTIS5'
      n = 4
      call rutis5_condition ( cond )
      write ( *, '(2x,a20,2x,i4,2x,g14.6)' ) title, n, cond
c
c  SUMMATION matrix.
c
      title = 'SUMMATION'
      n = 5
      call summation_condition ( n, cond )
      write ( *, '(2x,a20,2x,i4,2x,g14.6)' ) title, n, cond
c
c  TRI_UPPER matrix.
c
      title = 'TRI_UPPER'
      n = 5
      seed = 123456789
      alpha = r8_uniform_01 ( seed )
      alpha = anint ( 50.0D+00 * alpha - 25.0D+00 ) / 5.0D+00
      call tri_upper_condition ( alpha, n, cond )
      write ( *, '(2x,a20,2x,i4,2x,g14.6)' ) title, n, cond
c
c  WILK03 matrix.
c
      title = 'WILK03'
      n = 3
      call wilk03_condition ( cond )
      write ( *, '(2x,a20,2x,i4,2x,g14.6)' ) title, n, cond
c
c  WILSON matrix.
c
      title = 'WILSON'
      n = 4
      call wilson_condition ( cond )
      write ( *, '(2x,a20,2x,i4,2x,g14.6)' ) title, n, cond
 
      return
      end
      subroutine test_determinant ( )

c*********************************************************************72
c
cc TEST_DETERMINANT tests the determinant computations.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    07 July 2011
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n_max
      parameter ( n_max = 100 )

      double precision a(n_max,n_max)
      double precision alpha
      double precision b
      double precision beta
      integer col_num
      double precision d1
      double precision d2
      double precision d3
      double precision d4
      double precision d5
      double precision da
      double precision determ1
      double precision determ2
      double precision di
      double precision gamma
      integer i
      integer i4_uniform
      integer ii
      integer jj
      integer k
      double precision l(n_max,n_max)
      integer m
      integer n
      double precision norm_frobenius
      double precision p(n_max,n_max)
      double precision perturb
      integer pivot(n_max)
      double precision prob
      double precision r8_uniform_01
      double precision r8mat_norm_fro
      double precision r8vec_sum
      integer rank
      integer row_num
      integer seed
      integer seed_save
      character*80 title
      double precision u(n_max,n_max)
      double precision v1(n_max)
      double precision v2(n_max)
      double precision v3(n_max)
      double precision w(n_max)
      double precision x(2*n_max-1)
      integer x_n
      double precision y(n_max)
      integer y_n
      double precision y_sum
      double precision z(n_max)

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST_DETERMINANT'
      write ( *, '(a)' ) 
     &  '  Compute the determinants of an example of each'
      write ( *, '(a)' ) 
     &  '  test matrix; compare with the determinant routine,'
      write ( *, '(a)' ) 
     &  '  if available.  Print the matrix Frobenius norm'
      write ( *, '(a)' ) '  for an estimate of magnitude.'
      write ( *, '(a)' ) ' '
      write ( *, '(a,a)' ) '  Matrix title             N      ',
     &               'Determ          Determ           ||A||'
      write ( *, '(a)' ) ' '
c
c  AEGERTER matrix.
c
      title = 'AEGERTER'
      n = 5
      call aegerter ( n, a )
      call aegerter_determinant ( n, determ1 )
      call r8mat_determinant ( n, a, determ2 )
      norm_frobenius = r8mat_norm_fro ( n, n, a )
      write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) 
     &  title, n, determ1, determ2, norm_frobenius
c
c  ANTICIRCULANT matrix.
c
      title = 'ANTICIRCULANT'
      n = 3
      seed = 123456789
      call r8vec_uniform_01 ( n, seed, x )
      do i = 1, n
        x(i) = ( anint ( 50.0D+00 * x(i) - 25.0D+00 ) ) / 5.0D+00
      end do
      call anticirculant ( n, n, x, a )
      call anticirculant_determinant ( n, x, determ1 )
      call r8mat_determinant ( n, a, determ2 )
      norm_frobenius = r8mat_norm_fro ( n, n, a )
      write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) 
     &  title, n, determ1, determ2, norm_frobenius
c
c  ANTICIRCULANT matrix.
c
      title = 'ANTICIRCULANT'
      n = 4
      seed = 123456789
      call r8vec_uniform_01 ( n, seed, x )
      do i = 1, n
        x(i) = ( anint ( 50.0D+00 * x(i) - 25.0D+00 ) ) / 5.0D+00
      end do
      call anticirculant ( n, n, x, a )
      call anticirculant_determinant ( n, x, determ1 )
      call r8mat_determinant ( n, a, determ2 )
      norm_frobenius = r8mat_norm_fro ( n, n, a )
      write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) 
     &  title, n, determ1, determ2, norm_frobenius
c
c  ANTICIRCULANT matrix.
c
      title = 'ANTICIRCULANT'
      n = 5
      seed = 123456789
      call r8vec_uniform_01 ( n, seed, x )
      do i = 1, n
        x(i) = ( anint ( 50.0D+00 * x(i) - 25.0D+00 ) ) / 5.0D+00
      end do
      call anticirculant ( n, n, x, a )
      call anticirculant_determinant ( n, x, determ1 )
      call r8mat_determinant ( n, a, determ2 )
      norm_frobenius = r8mat_norm_fro ( n, n, a )
      write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) 
     &  title, n, determ1, determ2, norm_frobenius
c
c  ANTIHADAMARD matrix.
c
      title = 'ANTIHADAMARD'
      n = 5
      call antihadamard ( n, a )
      call antihadamard_determinant ( n, determ1 )
      call r8mat_determinant ( n, a, determ2 )
      norm_frobenius = r8mat_norm_fro ( n, n, a )
      write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) 
     &  title, n, determ1, determ2, norm_frobenius
c
c  ANTISYMM_RANDOM matrix.
c
      title = 'ANTISYMM_RANDOM'
      n = 5
      seed = 123456789
      call antisymm_random ( n, seed, a )
      call r8mat_determinant ( n, a, determ2 )
      norm_frobenius = r8mat_norm_fro ( n, n, a )
      write ( *, '(2x,a20,2x,i4,2x,14x,2x,g14.6,2x,g14.6)' )
     &  title, n,          determ2, norm_frobenius
c
c  ANTISYMM_RANDOM matrix.
c
      title = 'ANTISYMM_RANDOM'
      n = 6
      seed = 123456789
      call antisymm_random ( n, seed, a )
      call r8mat_determinant ( n, a, determ2 )
      norm_frobenius = r8mat_norm_fro ( n, n, a )
      write ( *, '(2x,a20,2x,i4,2x,14x,2x,g14.6,2x,g14.6)' )
     &  title, n,          determ2, norm_frobenius
c
c  BAB matrix.
c
      title = 'BAB'
      n = 5
      seed = 123456789
      alpha = anint ( 50.0D+00 * r8_uniform_01 ( seed ) - 25.0D+00 ) / 5.0D+00
      beta = anint ( 50.0D+00 * r8_uniform_01 ( seed ) - 25.0D+00 ) / 5.0D+00
      call bab ( n, alpha, beta, a )
      call bab_determinant ( n, alpha, beta, determ1 )
      call r8mat_determinant ( n, a, determ2 )
      norm_frobenius = r8mat_norm_fro ( n, n, a )
      write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) 
     &  title, n, determ1, determ2, norm_frobenius
c
c  BIMARKOV_RANDOM matrix.
c
      title = 'BIMARKOV_RANDOM'
      n = 5
      seed = 123456789
      call bimarkov_random ( n, seed, a )
      call r8mat_determinant ( n, a, determ2 )
      norm_frobenius = r8mat_norm_fro ( n, n, a )
      write ( *, '(2x,a20,2x,i4,2x,14x,2x,g14.6,2x,g14.6)' )
     &  title, n,          determ2, norm_frobenius
c
c  BIS matrix.
c
      title = 'BIS'
      n = 5
      seed = 123456789
      alpha = anint ( 50.0D+00 * r8_uniform_01 ( seed ) ) / 5.0D+00
      beta = anint ( 50.0D+00 * r8_uniform_01 ( seed ) ) / 5.0D+00
      call bis ( alpha, beta, n, n, a )
      call bis_determinant ( alpha, beta, n, determ1 )
      call r8mat_determinant ( n, a, determ2 )
      norm_frobenius = r8mat_norm_fro ( n, n, a )
      write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) 
     &  title, n, determ1, determ2, norm_frobenius
c
c  BODEWIG matrix.
c
      title = 'BODEWIG'
      n = 4
      call bodewig ( a )
      call bodewig_determinant ( determ1 )
      call r8mat_determinant ( n, a, determ2 )
      norm_frobenius = r8mat_norm_fro ( n, n, a )
      write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) 
     &  title, n, determ1, determ2, norm_frobenius
c
c  BOOTHROYD matrix.
c
      title = 'BOOTHROYD'
      n = 5
      call boothroyd ( n, a )
      call boothroyd_determinant ( n, determ1 )
      call r8mat_determinant ( n, a, determ2 )
      norm_frobenius = r8mat_norm_fro ( n, n, a )
      write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) 
     &  title, n, determ1, determ2, norm_frobenius
c
c  BORDERBAND matrix.
c
      title = 'BORDERBAND'
      n = 5
      call borderband ( n, a )
      call borderband_determinant ( n, determ1 )
      call r8mat_determinant ( n, a, determ2 )
      norm_frobenius = r8mat_norm_fro ( n, n, a )
      write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) 
     &  title, n, determ1, determ2, norm_frobenius
c
c  CARRY matrix.
c
      title = 'CARRY'
      n = 5
      seed = 123456789
      k = i4_uniform ( 2, 20, seed )
      call carry ( k, n, a )
      call carry_determinant ( k, n, determ1 )
      call r8mat_determinant ( n, a, determ2 )
      norm_frobenius = r8mat_norm_fro ( n, n, a )
      write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) 
     &  title, n, determ1, determ2, norm_frobenius
c
c  CAUCHY matrix.
c
      title = 'CAUCHY'
      n = 5
      seed = 123456789
      call r8vec_uniform_01 ( n, seed, x )
      call r8vec_uniform_01 ( n, seed, y )
      call cauchy ( n, x, y, a )
      call cauchy_determinant ( n, x, y, determ1 )
      call r8mat_determinant ( n, a, determ2 )
      norm_frobenius = r8mat_norm_fro ( n, n, a )
      write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) 
     &  title, n, determ1, determ2, norm_frobenius
c
c  CHEBY_DIFF1 matrix.
c
      title = 'CHEBY_DIFF1'
      n = 5
      call cheby_diff1 ( n, a )
      call r8mat_determinant ( n, a, determ2 )
      norm_frobenius = r8mat_norm_fro ( n, n, a )
      write ( *, '(2x,a20,2x,i4,2x,14x,2x,g14.6,2x,g14.6)' )
     &  title, n,          determ2, norm_frobenius
c
c  CHEBY_DIFF1 matrix.
c
      title = 'CHEBY_DIFF1'
      n = 6
      call cheby_diff1 ( n, a )
      call r8mat_determinant ( n, a, determ2 )
      norm_frobenius = r8mat_norm_fro ( n, n, a )
      write ( *, '(2x,a20,2x,i4,2x,14x,2x,g14.6,2x,g14.6)' )
     &  title, n,          determ2, norm_frobenius
c
c  CHEBY_T matrix.
c
      title = 'CHEBY_T'
      n = 5
      call cheby_t ( n, a )
      call cheby_t_determinant ( n, determ1 )
      call r8mat_determinant ( n, a, determ2 )
      norm_frobenius = r8mat_norm_fro ( n, n, a )
      write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) 
     &  title, n, determ1, determ2, norm_frobenius
c
c  CHEBY_U matrix.
c
      title = 'CHEBY_U'
      n = 5
      call cheby_u ( n, a )
      call cheby_u_determinant ( n, determ1 )
      call r8mat_determinant ( n, a, determ2 )
      norm_frobenius = r8mat_norm_fro ( n, n, a )
      write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) 
     &  title, n, determ1, determ2, norm_frobenius
c
c  CHEBY_VAN1 matrix.
c
      title = 'CHEBY_VAN1'
      n = 5
      seed = 123456789
      call r8vec_uniform_01 ( n, seed, x )
      do i = 1, n
        x(i) = anint ( 10.0D+00 * x(i) - 5.0D+00 )
      end do
      call cheby_van1 ( n, x, a )
      call r8mat_determinant ( n, a, determ2 )
      norm_frobenius = r8mat_norm_fro ( n, n, a )
      write ( *, '(2x,a20,2x,i4,2x,14x,2x,g14.6,2x,g14.6)' ) 
     &  title, n,          determ2, norm_frobenius
c
c  CHEBY_VAN2 matrix.
c
      do n = 2, 10
        title = 'CHEBY_VAN2'
        call cheby_van2 ( n, a )
        call cheby_van2_determinant ( n, determ1 )
        call r8mat_determinant ( n, a, determ2 )
        norm_frobenius = r8mat_norm_fro ( n, n, a )
        write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) 
     &    title, n, determ1, determ2, norm_frobenius
      end do
c
c  CHEBY_VAN3 matrix.
c
      title = 'CHEBY_VAN3'
      n = 5
      call cheby_van3 ( n, a )
      call cheby_van3_determinant ( n, determ1 )
      call r8mat_determinant ( n, a, determ2 )
      norm_frobenius = r8mat_norm_fro ( n, n, a )
      write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) 
     &  title, n, determ1, determ2, norm_frobenius
c
c  CHOW matrix.
c
      title = 'CHOW'
      n = 5
      seed = 123456789
      alpha = r8_uniform_01 ( seed );
      alpha = anint ( 50.0D+00 * alpha ) / 5.0D+00
      beta = r8_uniform_01 ( seed )
      beta = anint ( 50.0D+00 * beta ) / 5.0D+00
      call chow ( alpha, beta, n, n, a )
      call chow_determinant ( alpha, beta, n, determ1 )
      call r8mat_determinant ( n, a, determ2 )
      norm_frobenius = r8mat_norm_fro ( n, n, a )
      write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) 
     &  title, n, determ1, determ2, norm_frobenius
c
c  CIRCULANT matrix.
c
      title = 'CIRCULANT'
      n = 5
      seed = 123456789;
      call r8vec_uniform_01 ( n, seed, x )
      do i = 1, n
        x(i) = anint ( 50.0D+00 * x(i) - 25.0D+00 ) / 5.0D+00
      end do
      call circulant ( n, n, x, a )
      call circulant_determinant ( n, x, determ1 )
      call r8mat_determinant ( n, a, determ2 )
      norm_frobenius = r8mat_norm_fro ( n, n, a )
      write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) 
     &  title, n, determ1, determ2, norm_frobenius
c
c  CIRCULANT2 matrix.
c
      title = 'CIRCULANT2'
      n = 3
      call circulant2 ( n, a )
      call circulant2_determinant ( n, determ1 )
      call r8mat_determinant ( n, a, determ2 )
      norm_frobenius = r8mat_norm_fro ( n, n, a )
      write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) 
     &  title, n, determ1, determ2, norm_frobenius
c
c  CIRCULANT2 matrix.
c
      title = 'CIRCULANT2'
      n = 4
      call circulant2 ( n, a )
      call circulant2_determinant ( n, determ1 )
      call r8mat_determinant ( n, a, determ2 )
      norm_frobenius = r8mat_norm_fro ( n, n, a )
      write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) 
     &  title, n, determ1, determ2, norm_frobenius
c
c  CIRCULANT2 matrix.
c
      title = 'CIRCULANT2'
      n = 5
      call circulant2 ( n, a )
      call circulant2_determinant ( n, determ1 )
      call r8mat_determinant ( n, a, determ2 )
      norm_frobenius = r8mat_norm_fro ( n, n, a )
      write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) 
     &  title, n, determ1, determ2, norm_frobenius
c
c  CLEMENT1 matrix.
c
      title = 'CLEMENT1'
      n = 5
      call clement1 ( n, a )
      call clement1_determinant ( n, determ1 )
      call r8mat_determinant ( n, a, determ2 )
      norm_frobenius = r8mat_norm_fro ( n, n, a )
      write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) 
     &  title, n, determ1, determ2, norm_frobenius
c
c  CLEMENT1 matrix.
c
      title = 'CLEMENT1'
      n = 6
      call clement1 ( n, a )
      call clement1_determinant ( n, determ1 )
      call r8mat_determinant ( n, a, determ2 )
      norm_frobenius = r8mat_norm_fro ( n, n, a )
      write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) 
     &  title, n, determ1, determ2, norm_frobenius
c
c  CLEMENT2 matrix.
c
      title = 'CLEMENT2'
      n = 5
      call clement2 ( n, a )
      call clement2_determinant ( n, determ1 )
      call r8mat_determinant ( n, a, determ2 )
      norm_frobenius = r8mat_norm_fro ( n, n, a )
      write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) 
     &  title, n, determ1, determ2, norm_frobenius
c
c  CLEMENT2 matrix.
c
      title = 'CLEMENT2'
      n = 6
      call clement2 ( n, a )
      call clement2_determinant ( n, determ1 )
      call r8mat_determinant ( n, a, determ2 )
      norm_frobenius = r8mat_norm_fro ( n, n, a )
      write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) 
     &  title, n, determ1, determ2, norm_frobenius
c
c  CLEMENT3.
c
      title = 'CLEMENT3'
      n = 5
      seed = 123456789
      call r8vec_uniform_01 ( n-1, seed, x )
      do i = 1, n - 1
        x(i) = anint ( 10.0D+00 * x(i) - 5.0D+00 )
      end do
      call r8vec_uniform_01 ( n-1, seed, y )
      do i = 1, n - 1
        y(i) = anint ( 10.0D+00 * y(i) - 5.0D+00 )
      end do
      call clement3 ( n, x, y, a )
      call clement3_determinant ( n, x, y, determ1 )
      call r8mat_determinant ( n, a, determ2 )
      norm_frobenius = r8mat_norm_fro ( n, n, a )
      write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) 
     &  title, n, determ1, determ2, norm_frobenius
c
c  CLEMENT3.
c
      title = 'CLEMENT3'
      n = 6
      seed = 123456789
      call r8vec_uniform_01 ( n-1, seed, x )
      do i = 1, n - 1
        x(i) = anint ( 10.0D+00 * x(i) - 5.0D+00 )
      end do
      call r8vec_uniform_01 ( n-1, seed, y )
      do i = 1, n - 1
        y(i) = anint ( 10.0D+00 * y(i) - 5.0D+00 )
      end do
      call clement3 ( n, x, y, a )
      call clement3_determinant ( n, x, y, determ1 )
      call r8mat_determinant ( n, a, determ2 )
      norm_frobenius = r8mat_norm_fro ( n, n, a )
      write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) 
     &  title, n, determ1, determ2, norm_frobenius
c
c  COMBIN matrix.
c
      title = 'COMBIN'
      n = 5
      seed = 123456789
      alpha = r8_uniform_01 ( seed )
      alpha = anint ( 50.0D+00 * alpha ) / 5.0D+00
      beta = r8_uniform_01 ( seed )
      beta = anint ( 50.0D+00 * beta ) / 5.0D+00
      call combin ( alpha, beta, n, a )
      call combin_determinant ( alpha, beta, n, determ1 )
      call r8mat_determinant ( n, a, determ2 )
      norm_frobenius = r8mat_norm_fro ( n, n, a )
      write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) 
     &  title, n, determ1, determ2, norm_frobenius
c
c  COMPANION.
c
      title = 'COMPANION'
      n = 5
      seed = 123456789
      call r8vec_uniform_01 ( n, seed, x )
      do i = 1, n
        x(i) = anint ( 10.0D+00 * x(i) - 5.0D+00 )
      end do
      call companion ( n, x, a )
      call companion_determinant ( n, x, determ1 )
      call r8mat_determinant ( n, a, determ2 )
      norm_frobenius = r8mat_norm_fro ( n, n, a )
      write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) 
     &  title, n, determ1, determ2, norm_frobenius
c
c  COMPLEX_I matrix.
c
      title = 'COMPLEX_I'
      n = 2
      call complex_i ( a )
      call complex_i_determinant ( determ1 )
      call r8mat_determinant ( n, a, determ2 )
      norm_frobenius = r8mat_norm_fro ( n, n, a )
      write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) 
     &  title, n, determ1, determ2, norm_frobenius
c
c  CONEX1 matrix.
c
      title = 'CONEX1'
      n = 4
      seed = 123456789
      alpha = r8_uniform_01 ( seed )
      alpha = anint ( 50.0D+00 * alpha ) / 5.0D+00
      call conex1 ( alpha, a )
      call conex1_determinant ( alpha, determ1 )
      call r8mat_determinant ( n, a, determ2 )
      norm_frobenius = r8mat_norm_fro ( n, n, a )
      write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) 
     &  title, n, determ1, determ2, norm_frobenius
c
c  CONEX2 matrix.
c
      title = 'CONEX2'
      n = 3
      seed = 123456789
      alpha = r8_uniform_01 ( seed )
      alpha = anint ( 50.0D+00 * alpha ) / 5.0D+00
      call conex2 ( alpha, a )
      call conex2_determinant ( alpha, determ1 )
      call r8mat_determinant ( n, a, determ2 )
      norm_frobenius = r8mat_norm_fro ( n, n, a )
      write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) 
     &  title, n, determ1, determ2, norm_frobenius
c
c  CONEX3 matrix.
c
      title = 'CONEX3'
      n = 5
      call conex3 ( n, a )
      call conex3_determinant ( n, determ1 )
      call r8mat_determinant ( n, a, determ2 )
      norm_frobenius = r8mat_norm_fro ( n, n, a )
      write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) 
     &  title, n, determ1, determ2, norm_frobenius
c
c  CONFERENCE matrix.
c
      title = 'CONFERENCE'
      n = 6
      call conference ( n, a )
      call conference_determinant ( n, determ1 )
      call r8mat_determinant ( n, a, determ2 )
      norm_frobenius = r8mat_norm_fro ( n, n, a )
      write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) 
     &  title, n, determ1, determ2, norm_frobenius
c
c  CREATION matrix.
c
      title = 'CREATION'
      n = 5
      call creation ( n, n, a )
      call creation_determinant ( n, determ1 )
      call r8mat_determinant ( n, a, determ2 )
      norm_frobenius = r8mat_norm_fro ( n, n, a )
      write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) 
     &  title, n, determ1, determ2, norm_frobenius
c
c  DAUB2 matrix.
c
      title = 'DAUB2'
      n = 4
      call daub2 ( n, a )
      call daub2_determinant ( n, determ1 )
      call r8mat_determinant ( n, a, determ2 )
      norm_frobenius = r8mat_norm_fro ( n, n, a )
      write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) 
     &  title, n, determ1, determ2, norm_frobenius
c
c  DAUB4 matrix.
c
      title = 'DAUB4'
      n = 8
      call daub4 ( n, a )
      call daub4_determinant ( n, determ1 )
      call r8mat_determinant ( n, a, determ2 )
      norm_frobenius = r8mat_norm_fro ( n, n, a )
      write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) 
     &  title, n, determ1, determ2, norm_frobenius
c
c  DAUB6 matrix.
c
      title = 'DAUB6'
      n = 12
      call daub6 ( n, a )
      call daub6_determinant ( n, determ1 )
      call r8mat_determinant ( n, a, determ2 )
      norm_frobenius = r8mat_norm_fro ( n, n, a )
      write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) 
     &  title, n, determ1, determ2, norm_frobenius
c
c  DAUB8 matrix.
c
      title = 'DAUB8'
      n = 16
      call daub8 ( n, a )
      call daub8_determinant ( n, determ1 )
      call r8mat_determinant ( n, a, determ2 )
      norm_frobenius = r8mat_norm_fro ( n, n, a )
      write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) 
     &  title, n, determ1, determ2, norm_frobenius
c
c  DAUB10 matrix.
c
      title = 'DAUB10'
      n = 20
      call daub10 ( n, a )
      call daub10_determinant ( n, determ1 )
      call r8mat_determinant ( n, a, determ2 )
      norm_frobenius = r8mat_norm_fro ( n, n, a )
      write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) 
     &  title, n, determ1, determ2, norm_frobenius
c
c  DAUB12 matrix.
c
      title = 'DAUB12'
      n = 24
      call daub12 ( n, a )
      call daub12_determinant ( n, determ1 )
      call r8mat_determinant ( n, a, determ2 )
      norm_frobenius = r8mat_norm_fro ( n, n, a )
      write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) 
     &  title, n, determ1, determ2, norm_frobenius
c
c  DIAGONAL.
c
      title = 'DIAGONAL'
      n = 5
      seed = 123456789
      call r8vec_uniform_01 ( n, seed, x )
      do i = 1, n
        x(i) = anint ( 10.0D+00 * x(i) - 5.0D+00 )
      end do
      call diagonal ( n, n, x, a )
      call diagonal_determinant ( n, x, determ1 )
      call r8mat_determinant ( n, a, determ2 )
      norm_frobenius = r8mat_norm_fro ( n, n, a )
      write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) 
     &  title, n, determ1, determ2, norm_frobenius
c
c  DIF1 matrix.
c
      title = 'DIF1'
      n = 5
      call dif1 ( n, n, a )
      call dif1_determinant ( n, determ1 )
      call r8mat_determinant ( n, a, determ2 )
      norm_frobenius = r8mat_norm_fro ( n, n, a )
      write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) 
     &  title, n, determ1, determ2, norm_frobenius
c
c  DIF1CYCLIC matrix.
c
      title = 'DIF1CYCLIC'
      n = 5
      call dif1cyclic ( n, a )
      call dif1cyclic_determinant ( n, determ1 )
      call r8mat_determinant ( n, a, determ2 )
      norm_frobenius = r8mat_norm_fro ( n, n, a )
      write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) 
     &  title, n, determ1, determ2, norm_frobenius
c
c  DIF2 matrix.
c
      title = 'DIF2'
      n = 5
      call dif2 ( n, n, a )
      call dif2_determinant ( n, determ1 )
      call r8mat_determinant ( n, a, determ2 )
      norm_frobenius = r8mat_norm_fro ( n, n, a )
      write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) 
     &  title, n, determ1, determ2, norm_frobenius
c
c  DIF2CYCLIC matrix.
c
      title = 'DIF2CYCLIC'
      n = 5
      call dif2cyclic ( n, a )
      call dif2cyclic_determinant ( n, determ1 )
      call r8mat_determinant ( n, a, determ2 )
      norm_frobenius = r8mat_norm_fro ( n, n, a )
      write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) 
     &  title, n, determ1, determ2, norm_frobenius
c
c  DORR matrix.
c
      title = 'DORR'
      n = 5
      seed = 123456789
      alpha = r8_uniform_01 ( seed )
      alpha = anint ( 50.0D+00 * alpha ) / 5.0D+00
      call dorr ( alpha, n, a )
      call r8mat_determinant ( n, a, determ2 )
      norm_frobenius = r8mat_norm_fro ( n, n, a )
      write ( *, '(2x,a20,2x,i4,2x,14x,2x,g14.6,2x,g14.6)' ) 
     &  title, n,          determ2, norm_frobenius
c
c  DOWNSHIFT matrix.
c
      title = 'DOWNSHIFT'
      n = 5
      call downshift ( n, a )
      call downshift_determinant ( n, determ1 )
      call r8mat_determinant ( n, a, determ2 )
      norm_frobenius = r8mat_norm_fro ( n, n, a )
      write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) 
     &  title, n, determ1, determ2, norm_frobenius
c
c  EBERLEIN matrix.
c
      title = 'EBERLEIN'
      n = 5
      seed = 123456789
      alpha = r8_uniform_01 ( seed )
      alpha = anint ( 50.0D+00 * alpha ) / 5.0D+00
      call eberlein ( alpha, n, a )
      call eberlein_determinant ( alpha, n, determ1 )
      call r8mat_determinant ( n, a, determ2 )
      norm_frobenius = r8mat_norm_fro ( n, n, a )
      write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) 
     &  title, n, determ1, determ2, norm_frobenius
c
c  EULERIAN matrix.
c
      title = 'EULERIAN'
      n = 5
      call eulerian ( n, n, a )
      call eulerian_determinant ( n, determ1 )
      call r8mat_determinant ( n, a, determ2 )
      norm_frobenius = r8mat_norm_fro ( n, n, a )
      write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) 
     &  title, n, determ1, determ2, norm_frobenius
c
c  EXCHANGE matrix.
c
      title = 'EXCHANGE'
      n = 5
      call exchange ( n, n, a )
      call exchange_determinant ( n, determ1 )
      call r8mat_determinant ( n, a, determ2 )
      norm_frobenius = r8mat_norm_fro ( n, n, a )
      write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) 
     &  title, n, determ1, determ2, norm_frobenius
c
c  FIBONACCI1 matrix.
c
      title = 'FIBONACCI1'
      n = 5
      seed = 123456789
      alpha = r8_uniform_01 ( seed )
      alpha = anint ( 50.0D+00 * alpha ) / 5.0D+00
      beta = r8_uniform_01 ( seed )
      beta = anint ( 50.0D+00 * beta ) / 5.0D+00
      call fibonacci1 ( n, alpha, beta, a )
      call fibonacci1_determinant ( n, alpha, beta, determ1 )
      call r8mat_determinant ( n, a, determ2 )
      norm_frobenius = r8mat_norm_fro ( n, n, a )
      write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) 
     &  title, n, determ1, determ2, norm_frobenius
c
c  FIBONACCI2 matrix.
c
      title = 'FIBONACCI2'
      n = 5
      call fibonacci2 ( n, a )
      call fibonacci2_determinant ( n, determ1 )
      call r8mat_determinant ( n, a, determ2 )
      norm_frobenius = r8mat_norm_fro ( n, n, a )
      write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) 
     &  title, n, determ1, determ2, norm_frobenius
c
c  FIBONACCI3 matrix.
c
      title = 'FIBONACCI3'
      n = 5
      call fibonacci3 ( n, a )
      call fibonacci3_determinant ( n, determ1 )
      call r8mat_determinant ( n, a, determ2 )
      norm_frobenius = r8mat_norm_fro ( n, n, a )
      write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) 
     &  title, n, determ1, determ2, norm_frobenius
c
c  FIEDLER.
c
      title = 'FIEDLER'
      n = 7
      seed = 123456789
      call r8vec_uniform_01 ( n, seed, x )
      x(1:n) = anint ( 50.0D+00 * x(1:n) - 25.0D+00 ) / 5.0D+00
      call fiedler ( n, n, x, a )
      call fiedler_determinant ( n, x, determ1 )
      call r8mat_determinant ( n, a, determ2 )
      norm_frobenius = r8mat_norm_fro ( n, n, a )
      write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) 
     &  title, n, determ1, determ2, norm_frobenius
c
c  FORSYTHE matrix.
c
      title = 'FORSYTHE'
      n = 5
      seed = 123456789
      alpha = r8_uniform_01 ( seed )
      alpha = anint ( 50.0D+00 * alpha ) / 5.0D+00
      beta = r8_uniform_01 ( seed )
      beta = anint ( 50.0D+00 * beta ) / 5.0D+00
      call forsythe ( alpha, beta, n, a )
      call forsythe_determinant ( alpha, beta, n, determ1 )
      call r8mat_determinant ( n, a, determ2 )
      norm_frobenius = r8mat_norm_fro ( n, n, a )
      write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) 
     &  title, n, determ1, determ2, norm_frobenius
c
c  FOURIER_COSINE matrix.
c
      title = 'FOURIER_COSINE'
      n = 5
      call fourier_cosine ( n, a )
      call fourier_cosine_determinant ( n, determ1 )
      call r8mat_determinant ( n, a, determ2 )
      norm_frobenius = r8mat_norm_fro ( n, n, a )
      write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) 
     &  title, n, determ1, determ2, norm_frobenius
c
c  FOURIER_SINE matrix.
c
      title = 'FOURIER_SINE'
      n = 5
      call fourier_sine ( n, a )
      call fourier_sine_determinant ( n, determ1 )
      call r8mat_determinant ( n, a, determ2 )
      norm_frobenius = r8mat_norm_fro ( n, n, a )
      write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) 
     &  title, n, determ1, determ2, norm_frobenius
c
c  FRANK matrix.
c
      title = 'FRANK'
      n = 5
      call frank ( n, a )
      call frank_determinant ( n, determ1 )
      call r8mat_determinant ( n, a, determ2 )
      norm_frobenius = r8mat_norm_fro ( n, n, a )
      write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) 
     &  title, n, determ1, determ2, norm_frobenius
c
c  GEAR matrix.
c
      do n = 4, 8
        title = 'GEAR'
        seed = 123456789
        ii = i4_uniform ( -n, n, seed )
        jj = i4_uniform ( -n, n, seed )
        call gear ( ii, jj, n, a )
        call gear_determinant ( ii, jj, n, determ1 )
        call r8mat_determinant ( n, a, determ2 )
        norm_frobenius = r8mat_norm_fro ( n, n, a )
        write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) 
     &    title, n, determ1, determ2, norm_frobenius
      end do
c
c  GFPP matrix.
c
      title = 'GFPP'
      n = 5
      seed = 123456789
      alpha = r8_uniform_01 ( seed )
      alpha = anint ( 50.0D+00 * alpha ) / 5.0D+00
      call gfpp ( n, alpha, a )
      call gfpp_determinant ( n, alpha, determ1 )
      call r8mat_determinant ( n, a, determ2 )
      norm_frobenius = r8mat_norm_fro ( n, n, a )
      write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) 
     &  title, n, determ1, determ2, norm_frobenius
c
c  GIVENS matrix.
c
      title = 'GIVENS'
      n = 5
      call givens ( n, n, a )
      call givens_determinant ( n, determ1 )
      call r8mat_determinant ( n, a, determ2 )
      norm_frobenius = r8mat_norm_fro ( n, n, a )
      write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) 
     &  title, n, determ1, determ2, norm_frobenius
c
c  GK316 matrix.
c
      title = 'GK316'
      n = 5
      call gk316 ( n, a )
      call gk316_determinant ( n, determ1 )
      call r8mat_determinant ( n, a, determ2 )
      norm_frobenius = r8mat_norm_fro ( n, n, a )
      write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) 
     &  title, n, determ1, determ2, norm_frobenius
c
c  GK323 matrix.
c
      title = 'GK323'
      n = 5
      call gk323 ( n, n, a )
      call gk323_determinant ( n, determ1 )
      call r8mat_determinant ( n, a, determ2 )
      norm_frobenius = r8mat_norm_fro ( n, n, a )
      write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) 
     &  title, n, determ1, determ2, norm_frobenius
c
c  GK324 matrix.
c
      title = 'GK324'
      n = 5
      seed = 123456789
      call r8vec_uniform_01 ( n, seed, x )
      do i = 1, n
        x(i) = anint ( 50.0D+00 * x(i) - 25.0D+00 ) / 5.0D+00
      end do
      call gk324 ( n, n, x, a )
      call gk324_determinant ( n, x, determ1 )
      call r8mat_determinant ( n, a, determ2 )
      norm_frobenius = r8mat_norm_fro ( n, n, a )
      write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) 
     &  title, n, determ1, determ2, norm_frobenius
c
c  GRCAR matrix.
c
      title = 'GRCAR'
      n = 5
      seed = 123456789
      k = i4_uniform ( 1, n-1, seed )
      call grcar ( n, n, k, a )
      call r8mat_determinant ( n, a, determ2 )
      norm_frobenius = r8mat_norm_fro ( n, n, a )
      write ( *, '(2x,a20,2x,i4,2x,14x,2x,g14.6,2x,g14.6)' ) 
     &  title, n,          determ2, norm_frobenius
c
c  HADAMARD matrix.
c
      title = 'HADAMARD'
      n = 5
      call hadamard ( n, n, a )
      call r8mat_determinant ( n, a, determ2 )
      norm_frobenius = r8mat_norm_fro ( n, n, a )
      write ( *, '(2x,a20,2x,i4,2x,14x,2x,g14.6,2x,g14.6)' ) 
     &  title, n,          determ2, norm_frobenius
c
c  HANKEL matrix.
c
      title = 'HANKEL'
      n = 5
      seed = 123456789
      call r8vec_uniform_01 ( 2 * n - 1, seed, x )
      do i = 1, 2 * n - 1
        x(i) = anint ( 50.0D+00 * x(i) - 25.0D+00 ) / 5.0D+00
      end do
      call hankel ( n, x, a )
      call r8mat_determinant ( n, a, determ2 )
      norm_frobenius = r8mat_norm_fro ( n, n, a )
      write ( *, '(2x,a20,2x,i4,2x,14x,2x,g14.6,2x,g14.6)' ) 
     &  title, n,          determ2, norm_frobenius
c
c  HANOWA matrix.
c
      title = 'HANOWA'
      n = 6
      seed = 123456789
      alpha = r8_uniform_01 ( seed )
      alpha = anint ( 50.0D+00 * alpha ) / 5.0D+00
      call hanowa ( alpha, n, a )
      call hanowa_determinant ( alpha, n, determ1 )
      call r8mat_determinant ( n, a, determ2 )
      norm_frobenius = r8mat_norm_fro ( n, n, a )
      write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) 
     &  title, n, determ1, determ2, norm_frobenius
c
c  HARMAN matrix.
c
      title = 'HARMAN'
      n = 8
      call harman ( a )
      call harman_determinant ( determ1 )
      call r8mat_determinant ( n, a, determ2 )
      norm_frobenius = r8mat_norm_fro ( n, n, a )
      write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) 
     &  title, n, determ1, determ2, norm_frobenius
c
c  HARTLEY matrix.
c
      title = 'HARTLEY'
      do n = 5, 8
        call hartley ( n, a )
        call hartley_determinant ( n, determ1 )
        call r8mat_determinant ( n, a, determ2 )
        norm_frobenius = r8mat_norm_fro ( n, n, a )
        write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) 
     &    title, n, determ1, determ2, norm_frobenius
      end do
c
c  HELMERT matrix.
c
      title = 'HELMERT'
      n = 5
      call helmert ( n, a )
      call helmert_determinant ( n, determ1 )
      call r8mat_determinant ( n, a, determ2 )
      norm_frobenius = r8mat_norm_fro ( n, n, a )
      write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) 
     &  title, n, determ1, determ2, norm_frobenius
c
c  HELMERT2 matrix.
c
      title = 'HELMERT2'
      n = 5
      seed = 123456789
      call r8vec_uniform_01 ( n, seed, x )
      do i = 1, n
        x(i) = anint ( 50.0D+00 * x(i) - 25.0D+00 ) / 5.0D+00
      end do
      call helmert2 ( n, x, a )
      call r8mat_determinant ( n, a, determ2 )
      norm_frobenius = r8mat_norm_fro ( n, n, a )
      write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) 
     &  title, n, determ1, determ2, norm_frobenius
c
c  HERMITE matrix.
c
      title = 'HERMITE'
      n = 5
      call hermite ( n, a )
      call hermite_determinant ( n, determ1 )
      call r8mat_determinant ( n, a, determ2 )
      norm_frobenius = r8mat_norm_fro ( n, n, a )
      write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) 
     &  title, n, determ1, determ2, norm_frobenius
c
c  HERNDON matrix.
c
      title = 'HERNDON'
      n = 5
      call herndon ( n, a )
      call herndon_determinant ( n, determ1 )
      call r8mat_determinant ( n, a, determ2 )
      norm_frobenius = r8mat_norm_fro ( n, n, a )
      write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) 
     &  title, n, determ1, determ2, norm_frobenius
c
c  HILBERT matrix.
c
      title = 'HILBERT'
      n = 5
      call hilbert ( n, n, a )
      call hilbert_determinant ( n, determ1 )
      call r8mat_determinant ( n, a, determ2 )
      norm_frobenius = r8mat_norm_fro ( n, n, a )
      write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) 
     &  title, n, determ1, determ2, norm_frobenius
c
c  HOUSEHOLDER matrix.
c
      title = 'HOUSEHOLDER'
      n = 5
      seed = 123456789
      call r8vec_uniform_01 ( n, seed, x )
      do i = 1, n
        x(i) = anint ( 50.0D+00 * x(i) - 25.0D+00 ) / 5.0D+00
      end do
      call householder ( n, x, a )
      call householder_determinant ( n, x, determ1 )
      call r8mat_determinant ( n, a, determ2 )
      norm_frobenius = r8mat_norm_fro ( n, n, a )
      write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) 
     &  title, n, determ1, determ2, norm_frobenius
c
c  IDEM_RANDOM matrix.
c
      title = 'IDEM_RANDOM'
      n = 5
      seed = 123456789
      rank = i4_uniform ( 0, n, seed )
      call idem_random ( n, rank, seed, a )
      call idem_random_determinant ( n, rank, determ1 )
      call r8mat_determinant ( n, a, determ2 )
      norm_frobenius = r8mat_norm_fro ( n, n, a )
      write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) 
     &  title, n, determ1, determ2, norm_frobenius
c
c  IDENTITY matrix.
c
      title = 'IDENTITY'
      n = 5
      call identity ( n, n, a )
      call identity_determinant ( n, determ1 )
      call r8mat_determinant ( n, a, determ2 )
      norm_frobenius = r8mat_norm_fro ( n, n, a )
      write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) 
     &  title, n, determ1, determ2, norm_frobenius
c
c  IJFACT1 matrix.
c
      title = 'IJFACT1'
      n = 5
      call ijfact1 ( n, a )
      call ijfact1_determinant ( n, determ1 )
      call r8mat_determinant ( n, a, determ2 )
      norm_frobenius = r8mat_norm_fro ( n, n, a )
      write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) 
     &  title, n, determ1, determ2, norm_frobenius
c
c  IJFACT2 matrix.
c
      title = 'IJFACT2'
      n = 5
      call ijfact2 ( n, a )
      call ijfact2_determinant ( n, determ1 )
      call r8mat_determinant ( n, a, determ2 )
      norm_frobenius = r8mat_norm_fro ( n, n, a )
      write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) 
     &  title, n, determ1, determ2, norm_frobenius
c
c  ILL3 matrix.
c
      title = 'ILL3'
      n = 3
      call ill3 ( a )
      call ill3_determinant ( determ1 )
      call r8mat_determinant ( n, a, determ2 )
      norm_frobenius = r8mat_norm_fro ( n, n, a )
      write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) 
     &  title, n, determ1, determ2, norm_frobenius
c
c  INTEGRATION matrix.
c
      title = 'INTEGRATION'
      n = 6
      seed = 123456789
      alpha = r8_uniform_01 ( seed )
      alpha = anint ( 50.0D+00 * alpha ) / 5.0D+00
      call integration ( alpha, n, a )
      call integration_determinant ( alpha, n, determ1 )
      call r8mat_determinant ( n, a, determ2 )
      norm_frobenius = r8mat_norm_fro ( n, n, a )
      write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) 
     &  title, n, determ1, determ2, norm_frobenius
c
c  INVOL matrix.
c
      title = 'INVOL'
      n = 5
      call invol ( n, a )
      call invol_determinant ( n, determ1 )
      call r8mat_determinant ( n, a, determ2 )
      norm_frobenius = r8mat_norm_fro ( n, n, a )
      write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) 
     &  title, n, determ1, determ2, norm_frobenius
c
c  INVOL_RANDOM matrix.
c
      title = 'INVOL_RANDOM'
      n = 5
      seed = 123456789
      rank = i4_uniform ( 0, n, seed )
      call invol_random ( n, rank, seed, a )
      call r8mat_determinant ( n, a, determ2 )
      norm_frobenius = r8mat_norm_fro ( n, n, a )
      write ( *, '(2x,a20,2x,i4,2x,14x,2x,g14.6,2x,g14.6)' ) 
     &  title, n,          determ2, norm_frobenius
c
c  JACOBI matrix.
c
      title = 'JACOBI'
      n = 5
      call jacobi ( n, n, a )
      call jacobi_determinant ( n, determ1 )
      call r8mat_determinant ( n, a, determ2 )
      norm_frobenius = r8mat_norm_fro ( n, n, a )
      write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) 
     &  title, n, determ1, determ2, norm_frobenius
c
c  JACOBI matrix.
c
      title = 'JACOBI'
      n = 6
      call jacobi ( n, n, a )
      call jacobi_determinant ( n, determ1 )
      call r8mat_determinant ( n, a, determ2 )
      norm_frobenius = r8mat_norm_fro ( n, n, a )
      write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) 
     &  title, n, determ1, determ2, norm_frobenius
c
c  JORDAN matrix.
c
      title = 'JORDAN'
      n = 6
      seed = 123456789
      alpha = r8_uniform_01 ( seed )
      alpha = anint ( 50.0D+00 * alpha ) / 5.0D+00
      call jordan ( alpha, n, n, a )
      call jordan_determinant ( alpha, n, determ1 )
      call r8mat_determinant ( n, a, determ2 )
      norm_frobenius = r8mat_norm_fro ( n, n, a )
      write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) 
     &  title, n, determ1, determ2, norm_frobenius
c
c  KAHAN matrix.
c
      title = 'KAHAN'
      n = 5
      seed = 123456789
      alpha = r8_uniform_01 ( seed )
      alpha = anint ( 50.0D+00 * alpha ) / 5.0D+00
      call kahan ( alpha, n, n, a )
      call kahan_determinant ( alpha, n, determ1 )
      call r8mat_determinant ( n, a, determ2 )
      norm_frobenius = r8mat_norm_fro ( n, n, a )
      write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) 
     &  title, n, determ1, determ2, norm_frobenius
c
c  KERSHAW matrix.
c
      title = 'KERSHAW'
      n = 4
      call kershaw ( a )
      call kershaw_determinant ( determ1 )
      call r8mat_determinant ( n, a, determ2 )
      norm_frobenius = r8mat_norm_fro ( n, n, a )
      write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) 
     &  title, n, determ1, determ2, norm_frobenius
c
c  KERSHAWTRI matrix.
c
      title = 'KERSHAWTRI'
      n = 5
      x_n = ( n + 1 ) / 2
      seed = 123456789
      call r8vec_uniform_01 ( x_n, seed, x )
      do i = 1, x_n
        x(i) = anint ( 50.0D+00 * x(i) - 25.0D+00 ) / 5.0D+00
      end do
      call kershawtri ( n, x, a )
      call kershawtri_determinant ( n, x, determ1 )
      call r8mat_determinant ( n, a, determ2 )
      norm_frobenius = r8mat_norm_fro ( n, n, a )
      write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) 
     &  title, n, determ1, determ2, norm_frobenius
c
c  KMS matrix.
c
      title = 'KMS'
      n = 5
      seed = 123456789
      alpha = r8_uniform_01 ( seed )
      alpha = anint ( 50.0D+00 * alpha ) / 5.0D+00
      call kms ( alpha, n, n, a )
      call kms_determinant ( alpha, n, determ1 )
      call r8mat_determinant ( n, a, determ2 )
      norm_frobenius = r8mat_norm_fro ( n, n, a )
      write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) 
     &  title, n, determ1, determ2, norm_frobenius
c
c  LAGUERRE matrix.
c
      title = 'LAGUERRE'
      n = 5
      call laguerre ( n, a )
      call laguerre_determinant ( n, determ1 )
      call r8mat_determinant ( n, a, determ2 )
      norm_frobenius = r8mat_norm_fro ( n, n, a )
      write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) 
     &  title, n, determ1, determ2, norm_frobenius
c
c  LEHMER matrix.
c
      title = 'LEHMER'
      n = 5
      call lehmer ( n, n, a )
      call r8mat_determinant ( n, a, determ2 )
      norm_frobenius = r8mat_norm_fro ( n, n, a )
      write ( *, '(2x,a20,2x,i4,2x,14x,2x,g14.6,2x,g14.6)' ) 
     &  title, n,          determ2, norm_frobenius
c
c  LESLIE matrix.
c
      title = 'LESLIE'
      n = 4
      b =  0.025D+00
      di = 0.010D+00
      da = 0.100D+00
      call leslie ( b, di, da, a )
      call leslie_determinant ( b, di, da, determ1 )
      call r8mat_determinant ( n, a, determ2 )
      norm_frobenius = r8mat_norm_fro ( n, n, a )
      write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) 
     &  title, n, determ1, determ2, norm_frobenius
c
c  LESP matrix.
c
      title = 'LESP'
      n = 5
      call lesp ( n, n, a )
      call lesp_determinant ( n, determ1 )
      call r8mat_determinant ( n, a, determ2 )
      norm_frobenius = r8mat_norm_fro ( n, n, a )
      write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) 
     &  title, n, determ1, determ2, norm_frobenius
c
c  LIETZKE matrix.
c
      title = 'LIETZKE'
      n = 5
      call lietzke ( n, a )
      call lietzke_determinant ( n, determ1 )
      call r8mat_determinant ( n, a, determ2 )
      norm_frobenius = r8mat_norm_fro ( n, n, a )
      write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) 
     &  title, n, determ1, determ2, norm_frobenius
c
c  LIGHTS_OUT matrix.
c
      title = 'LIGHTS_OUT'
      row_num = 5
      col_num = 5
      n = row_num * col_num;
      call lights_out ( row_num, col_num, n, a )
      call r8mat_determinant ( n, a, determ2 )
      norm_frobenius = r8mat_norm_fro ( n, n, a )
      write ( *, '(2x,a20,2x,i4,2x,14x,2x,g14.6,2x,g14.6)' ) 
     &  title, n,          determ2, norm_frobenius
c
c  LINE_ADJ matrix.
c
      title = 'LINE_ADJ'
      n = 5
      call line_adj ( n, a )
      call line_adj_determinant ( n, determ1 )
      call r8mat_determinant ( n, a, determ2 )
      norm_frobenius = r8mat_norm_fro ( n, n, a )
      write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) 
     &  title, n, determ1, determ2, norm_frobenius
c
c  LINE_LOOP_ADJ matrix.
c
      title = 'LINE_LOOP_ADJ'
      n = 5
      call line_loop_adj ( n, a )
      call line_loop_adj_determinant ( n, determ1 )
      call r8mat_determinant ( n, a, determ2 )
      norm_frobenius = r8mat_norm_fro ( n, n, a )
      write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) 
     &  title, n, determ1, determ2, norm_frobenius
c
c  LOEWNER matrix.
c
      title = 'LOEWNER'
      n = 5
      seed = 123456789
      call r8vec_uniform_01 ( n, seed, w )
      call r8vec_uniform_01 ( n, seed, x )
      call r8vec_uniform_01 ( n, seed, y )
      call r8vec_uniform_01 ( n, seed, z )
      call loewner ( w, x, y, z, n, a )
      call r8mat_determinant ( n, a, determ2 )
      norm_frobenius = r8mat_norm_fro ( n, n, a )
      write ( *, '(2x,a20,2x,i4,2x,14x,2x,g14.6,2x,g14.6)' ) 
     &  title, n,          determ2, norm_frobenius
c
c  LOTKIN matrix.
c
      title = 'LOTKIN'
      n = 5
      call lotkin ( n, n, a )
      call lotkin_determinant ( n, determ1 )
      call r8mat_determinant ( n, a, determ2 )
      norm_frobenius = r8mat_norm_fro ( n, n, a )
      write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) 
     &  title, n, determ1, determ2, norm_frobenius
c
c  MARKOV_RANDOM matrix.
c
      title = 'MARKOV_RANDOM'
      n = 5
      seed = 123456789
      call markov_random ( n, seed, a )
      call r8mat_determinant ( n, a, determ2 )
      norm_frobenius = r8mat_norm_fro ( n, n, a )
      write ( *, '(2x,a20,2x,i4,2x,14x,2x,g14.6,2x,g14.6)' ) 
     &  title, n,          determ2, norm_frobenius
c
c  MAXIJ matrix.
c
      title = 'MAXIJ'
      n = 5
      call maxij ( n, n, a )
      call maxij_determinant ( n, determ1 )
      call r8mat_determinant ( n, a, determ2 )
      norm_frobenius = r8mat_norm_fro ( n, n, a )
      write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) 
     &  title, n, determ1, determ2, norm_frobenius
c
c  MILNES matrix.
c
      title = 'MILNES'
      n = 5
      seed = 123456789
      call r8vec_uniform_01 ( n, seed, x )
      do i = 1, n
        x(i) = anint ( 50.0D+00 * x(i) - 25.0D+00 ) / 5.0D+00
      end do
      call milnes ( n, n, x, a )
      call milnes_determinant ( n, x, determ1 )
      call r8mat_determinant ( n, a, determ2 )
      norm_frobenius = r8mat_norm_fro ( n, n, a )
      write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) 
     &  title, n, determ1, determ2, norm_frobenius
c
c  MINIJ matrix.
c
      title = 'MINIJ'
      n = 5
      call minij ( n, n, a )
      call minij_determinant ( n, determ1 )
      call r8mat_determinant ( n, a, determ2 )
      norm_frobenius = r8mat_norm_fro ( n, n, a )
      write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) 
     &  title, n, determ1, determ2, norm_frobenius
c
c  MOLER1 matrix.
c
      title = 'MOLER1'
      n = 5
      seed = 123456789
      alpha = r8_uniform_01 ( seed )
      alpha = anint ( 50.0D+00 * alpha ) / 5.0D+00
      call moler1 ( alpha, n, n, a )
      call moler1_determinant ( alpha, n, determ1 )
      call r8mat_determinant ( n, a, determ2 )
      norm_frobenius = r8mat_norm_fro ( n, n, a )
      write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) 
     &  title, n, determ1, determ2, norm_frobenius
c
c  MOLER2 matrix.
c
      title = 'MOLER2'
      n = 5
      call moler2 ( a )
      call moler2_determinant ( determ1 )
      call r8mat_determinant ( n, a, determ2 )
      norm_frobenius = r8mat_norm_fro ( n, n, a )
      write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) 
     &  title, n, determ1, determ2, norm_frobenius
c
c  MOLER3 matrix.
c
      title = 'MOLER3'
      n = 5
      call moler3 ( n, n, a )
      call moler3_determinant ( n, determ1 )
      call r8mat_determinant ( n, a, determ2 )
      norm_frobenius = r8mat_norm_fro ( n, n, a )
      write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) 
     &  title, n, determ1, determ2, norm_frobenius
c
c  NEUMANN matrix.
c
      title = 'NEUMANN'
      row_num = 5
      col_num = 5
      n = row_num * col_num
      call neumann ( row_num, col_num, n, a )
      call neumann_determinant ( n, determ1 )
      call r8mat_determinant ( n, a, determ2 )
      norm_frobenius = r8mat_norm_fro ( n, n, a )
      write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) 
     &  title, n, determ1, determ2, norm_frobenius
c
c  ONE matrix.
c
      title = 'ONE'
      n = 5
      call one ( n, n, a )
      call one_determinant ( n, determ1 )
      call r8mat_determinant ( n, a, determ2 )
      norm_frobenius = r8mat_norm_fro ( n, n, a )
      write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) 
     &  title, n, determ1, determ2, norm_frobenius
c
c  ORTEGA matrix.
c
      title = 'ORTEGA'
      n = 5
      seed = 123456789
      call r8vec_uniform_01 ( n, seed, v1 )
      do i = 1, n
        v1(i) = anint ( 50.0D+00 * v1(i) - 25.0D+00  ) / 5.0D+00
      end do
      call r8vec_uniform_01 ( n, seed, v2 )
      do i = 1, n
        v2(i) = anint ( 50.0D+00 * v2(i) - 25.0D+00  ) / 5.0D+00
      end do
      call r8vec_uniform_01 ( n, seed, v3 )
      do i = 1, n
        v3(i) = anint ( 50.0D+00 * v3(i) - 25.0D+00  ) / 5.0D+00
      end do
      call ortega ( n, v1, v2, v3, a )
      call ortega_determinant ( n, v1, v2, v3, determ1 )
      call r8mat_determinant ( n, a, determ2 )
      norm_frobenius = r8mat_norm_fro ( n, n, a )
      write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) 
     &  title, n, determ1, determ2, norm_frobenius
c
c  ORTH_RANDOM matrix.
c
      title = 'ORTH_RANDOM'
      n = 5
      seed = 123456789
      call orth_random ( n, seed, a )
      call orth_random_determinant ( n, seed, determ1 )
      call r8mat_determinant ( n, a, determ2 )
      norm_frobenius = r8mat_norm_fro ( n, n, a )
      write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) 
     &  title, n, determ1, determ2, norm_frobenius
c
c  ORTH_SYMM matrix.
c
      title = 'ORTH_SYMM'
      n = 5
      call orth_symm ( n, a )
      call orth_symm_determinant ( n, determ1 )
      call r8mat_determinant ( n, a, determ2 )
      norm_frobenius = r8mat_norm_fro ( n, n, a )
      write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) 
     &  title, n, determ1, determ2, norm_frobenius
c
c  OTO matrix.
c
      title = 'OTO'
      n = 5
      call oto ( n, n, a )
      call oto_determinant ( n, determ1 )
      call r8mat_determinant ( n, a, determ2 )
      norm_frobenius = r8mat_norm_fro ( n, n, a )
      write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) 
     &  title, n, determ1, determ2, norm_frobenius
c
c  PARTER matrix.
c
      title = 'PARTER'
      n = 5
      call parter ( n, n, a )
      call parter_determinant ( n, determ1 )
      call r8mat_determinant ( n, a, determ2 )
      norm_frobenius = r8mat_norm_fro ( n, n, a )
      write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) 
     &  title, n, determ1, determ2, norm_frobenius
c
c  PASCAL1 matrix.
c
      title = 'PASCAL1'
      n = 5
      call pascal1 ( n, a )
      call pascal1_determinant ( n, determ1 )
      call r8mat_determinant ( n, a, determ2 )
      norm_frobenius = r8mat_norm_fro ( n, n, a )
      write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) 
     &  title, n, determ1, determ2, norm_frobenius
c
c  PASCAL2 matrix.
c
      title = 'PASCAL2'
      n = 5
      call pascal2 ( n, a )
      call pascal2_determinant ( n, determ1 )
      call r8mat_determinant ( n, a, determ2 )
      norm_frobenius = r8mat_norm_fro ( n, n, a )
      write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) 
     &  title, n, determ1, determ2, norm_frobenius
c
c  PASCAL3 matrix.
c
      title = 'PASCAL3'
      n = 5
      seed = 123456789
      alpha = r8_uniform_01 ( seed )
      alpha = anint ( 50.0D+00 * alpha - 25.0D+00 ) / 5.0D+00
      call pascal3 ( n, alpha, a )
      call pascal3_determinant ( n, alpha, determ1 )
      call r8mat_determinant ( n, a, determ2 )
      norm_frobenius = r8mat_norm_fro ( n, n, a )
      write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) 
     &  title, n, determ1, determ2, norm_frobenius
c
c  PDS_RANDOM matrix.
c
      title = 'PDS_RANDOM'
      n = 5
      seed = 123456789
      seed_save = seed
      call pds_random ( n, seed, a )
      seed = seed_save
      call pds_random_determinant ( n, seed, determ1 )
      call r8mat_determinant ( n, a, determ2 )
      norm_frobenius = r8mat_norm_fro ( n, n, a )
      write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) 
     &  title, n, determ1, determ2, norm_frobenius
c
c  PEI matrix.
c
      title = 'PEI'
      n = 5
      seed = 123456789
      alpha = r8_uniform_01 ( seed )
      alpha = anint ( 50.0D+00 * alpha ) / 5.0D+00
      call pei ( alpha, n, a )
      call pei_determinant ( alpha, n, determ1 )
      call r8mat_determinant ( n, a, determ2 )
      norm_frobenius = r8mat_norm_fro ( n, n, a )
      write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) 
     &  title, n, determ1, determ2, norm_frobenius
c
c  PERMUTATION_RANDOM matrix.
c
      title = 'PERMUTATION_RANDOM'
      n = 5
      seed = 123456789
      seed_save = seed
      call permutation_random ( n, seed, a )
      seed = seed_save
      call permutation_random_determinant ( n, seed, determ1 )
      call r8mat_determinant ( n, a, determ2 )
      norm_frobenius = r8mat_norm_fro ( n, n, a )
      write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) 
     &  title, n, determ1, determ2, norm_frobenius
c
c  PLU matrix.
c
      title = 'PLU'
      n = 5
      do i = 1, n
        pivot(i) = i
      end do
      call plu ( n, pivot, p, l, u, a )
      call plu_determinant ( n, p, l, u, determ1 )
      call r8mat_determinant ( n, a, determ2 )
      norm_frobenius = r8mat_norm_fro ( n, n, a )
      write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) 
     &  title, n, determ1, determ2, norm_frobenius
c
c  POISSON matrix.
c
      title = 'POISSON'
      row_num = 5
      col_num = 5
      n = row_num * col_num
      call poisson ( row_num, col_num, n, a )
      call poisson_determinant ( row_num, col_num, n, determ1 )
      call r8mat_determinant ( n, a, determ2 )
      norm_frobenius = r8mat_norm_fro ( n, n, a )
      write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) 
     &  title, n, determ1, determ2, norm_frobenius
c
c  PROLATE matrix.
c
      title = 'PROLATE'
      n = 5
      seed = 123456789
      alpha = r8_uniform_01 ( seed )
      alpha = anint ( 50.0D+00 * alpha ) / 5.0D+00
      call prolate ( alpha, n, a )
      call r8mat_determinant ( n, a, determ2 )
      norm_frobenius = r8mat_norm_fro ( n, n, a )
      write ( *, '(2x,a20,2x,i4,2x,14x,2x,g14.6,2x,g14.6)' ) 
     &  title, n,          determ2, norm_frobenius
c
c  RECTANGLE_ADJ matrix.
c
      title = 'RECTANGLE_ADJ'
      row_num = 5
      col_num = 5
      n = row_num * col_num
      call rectangle_adj ( row_num, col_num, n, a )
      call rectangle_adj_determinant ( row_num, col_num, determ1 )
      call r8mat_determinant ( n, a, determ2 )
      norm_frobenius = r8mat_norm_fro ( n, n, a )
      write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) 
     &  title, n, determ1, determ2, norm_frobenius
c
c  REDHEFFER matrix.
c
      title = 'REDHEFFER'
      n = 5
      call redheffer ( n, a )
      call redheffer_determinant ( n, determ1 )
      call r8mat_determinant ( n, a, determ2 )
      norm_frobenius = r8mat_norm_fro ( n, n, a )
      write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) 
     &  title, n, determ1, determ2, norm_frobenius
c
c  REF_RANDOM matrix.
c
      title = 'REF_RANDOM'
      n = 5
      prob = 0.65D+00
      seed_save = 123456789
      seed = seed_save
      call ref_random ( n, n, prob, seed, a )
      seed = seed_save
      call ref_random_determinant ( n, prob, seed, determ1 )
      call r8mat_determinant ( n, a, determ2 )
      norm_frobenius = r8mat_norm_fro ( n, n, a )
      write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) 
     &  title, n, determ1, determ2, norm_frobenius
c
c  REF_RANDOM matrix.
c
      title = 'REF_RANDOM'
      n = 5
      prob = 0.85D+00
      seed_save = 123456789
      seed = seed_save
      call ref_random ( n, n, prob, seed, a )
      seed = seed_save
      call ref_random_determinant ( n, prob, seed, determ1 )
      call r8mat_determinant ( n, a, determ2 )
      norm_frobenius = r8mat_norm_fro ( n, n, a )
      write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) 
     &  title, n, determ1, determ2, norm_frobenius
c
c  RIEMANN matrix.
c
      title = 'RIEMANN'
      n = 5
      call riemann ( n, n, a )
      call r8mat_determinant ( n, a, determ2 )
      norm_frobenius = r8mat_norm_fro ( n, n, a )
      write ( *, '(2x,a20,2x,i4,2x,14x,2x,g14.6,2x,g14.6)' ) 
     &  title, n,          determ2, norm_frobenius
c
c  RING_ADJ matrix.
c
      do n = 1, 8
        title = 'RING_ADJ'
        call ring_adj ( n, a )
        call ring_adj_determinant ( n, determ1 )
        call r8mat_determinant ( n, a, determ2 )
        norm_frobenius = r8mat_norm_fro ( n, n, a )
        write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) 
     &    title, n, determ1, determ2, norm_frobenius
      end do
c
c  RIS matrix.
c
      title = 'RIS'
      n = 5
      call ris ( n, a )
      call ris_determinant ( n, determ1 )
      call r8mat_determinant ( n, a, determ2 )
      norm_frobenius = r8mat_norm_fro ( n, n, a )
      write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) 
     &  title, n, determ1, determ2, norm_frobenius
c
c  RODMAN matrix.
c
      title = 'RODMAN'
      n = 5
      seed = 123456789
      alpha = r8_uniform_01 ( seed )
      alpha = anint ( 50.0D+00 * alpha - 25.0D+00 ) / 5.0D+00
      call rodman ( alpha, n, n, a )
      call rodman_determinant ( alpha, n, determ1 )
      call r8mat_determinant ( n, a, determ2 )
      norm_frobenius = r8mat_norm_fro ( n, n, a )
      write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) 
     &  title, n, determ1, determ2, norm_frobenius
c
c  ROSSER1 matrix.
c
c  Note that while the correct determinant of this matrix is 0,
c  that value is very difficult to calculate correctly.  MATLAB
c  returns det ( A ) = -10611, for instance.
c
      title = 'ROSSER1'
      n = 8
      call rosser1 ( a )
      call rosser1_determinant ( determ1 )
      call r8mat_determinant ( n, a, determ2 )
      norm_frobenius = r8mat_norm_fro ( n, n, a )
      write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) 
     &  title, n, determ1, determ2, norm_frobenius
c
c  ROUTH matrix.
c
      title = 'ROUTH'
      n = 5
      seed = 123456789
      call r8vec_uniform_01 ( n, seed, x )
      do i = 1, n
        x(i) = anint ( 50.0D+00 * x(i) ) / 5.0D+00
      end do
      call routh ( n, x, a )
      call routh_determinant ( n, x, determ1 )
      call r8mat_determinant ( n, a, determ2 )
      norm_frobenius = r8mat_norm_fro ( n, n, a )
      write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) 
     &  title, n, determ1, determ2, norm_frobenius
c
c  RUTIS1 matrix.
c
      title = 'RUTIS1'
      n = 4
      call rutis1 ( a )
      call rutis1_determinant ( determ1 )
      call r8mat_determinant ( n, a, determ2 )
      norm_frobenius = r8mat_norm_fro ( n, n, a )
      write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) 
     &  title, n, determ1, determ2, norm_frobenius
c
c  RUTIS2 matrix.
c
      title = 'RUTIS2'
      n = 4
      call rutis2 ( a )
      call rutis2_determinant ( determ1 )
      call r8mat_determinant ( n, a, determ2 )
      norm_frobenius = r8mat_norm_fro ( n, n, a )
      write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) 
     &  title, n, determ1, determ2, norm_frobenius
c
c  RUTIS3 matrix.
c
      title = 'RUTIS3'
      n = 4
      call rutis3 ( a )
      call rutis3_determinant ( determ1 )
      call r8mat_determinant ( n, a, determ2 )
      norm_frobenius = r8mat_norm_fro ( n, n, a )
      write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) 
     &  title, n, determ1, determ2, norm_frobenius
c
c  RUTIS4 matrix.
c
      title = 'RUTIS4'
      n = 5
      call rutis4 ( n, a )
      call rutis4_determinant ( n, determ1 )
      call r8mat_determinant ( n, a, determ2 )
      norm_frobenius = r8mat_norm_fro ( n, n, a )
      write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) 
     &  title, n, determ1, determ2, norm_frobenius
c
c  RUTIS5 matrix.
c
      title = 'RUTIS5'
      n = 4
      call rutis5 ( a )
      call rutis5_determinant ( determ1 )
      call r8mat_determinant ( n, a, determ2 )
      norm_frobenius = r8mat_norm_fro ( n, n, a )
      write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) 
     &  title, n, determ1, determ2, norm_frobenius
c
c  SCHUR_BLOCK matrix.
c
      title = 'SCHUR_BLOCK'
      n = 5
      x_n = ( n + 1 ) / 2
      y_n = n / 2
      seed = 123456789
      call r8vec_uniform_01 ( x_n, seed, x )
      x(1:x_n) = anint ( 50.0D+00 * x(1:x_n) - 25.0D+00 ) / 5.0D+00
      call r8vec_uniform_01 ( y_n, seed, y )
      y(1:y_n) = anint ( 50.0D+00 * y(1:y_n) - 25.0D+00 ) / 5.0D+00
      call schur_block ( n, x, y, a )
      call schur_block_determinant ( n, x, y, determ1 )
      call r8mat_determinant ( n, a, determ2 )
      norm_frobenius = r8mat_norm_fro ( n, n, a )
      write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) 
     &  title, n, determ1, determ2, norm_frobenius
c
c  SKEW_CIRCULANT matrix.
c
      title = 'SKEW_CIRCULANT'
      n = 5
      seed = 123456789
      call r8vec_uniform_01 ( n, seed, x )
      do i = 1, n
        x(i) = anint ( 50.0D+00 * x(i) - 25.0D+00 ) / 5.0D+00
      end do
      call skew_circulant ( n, n, x, a )
      call skew_circulant_determinant ( n, x, determ1 )
      call r8mat_determinant ( n, a, determ2 )
      norm_frobenius = r8mat_norm_fro ( n, n, a )
      write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) 
     &  title, n, determ1, determ2, norm_frobenius
c
c  SPLINE matrix.
c
      title = 'SPLINE'
      n = 5
      seed = 123456789
      call r8vec_uniform_01 ( n, seed, x )
      do i = 1, n
        x(i) = anint ( 50.0D+00 * x(i) - 25.0D+00 ) / 5.0D+00
      end do
      call spline ( n, x, a )
      call spline_determinant ( n, x, determ1 )
      call r8mat_determinant ( n, a, determ2 )
      norm_frobenius = r8mat_norm_fro ( n, n, a )
      write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) 
     &  title, n, determ1, determ2, norm_frobenius
c
c  STIRLING matrix.
c
      title = 'STIRLING'
      n = 5
      call stirling ( n, n, a )
      call stirling_determinant ( n, determ1 )
      call r8mat_determinant ( n, a, determ2 )
      norm_frobenius = r8mat_norm_fro ( n, n, a )
      write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) 
     &  title, n, determ1, determ2, norm_frobenius
c
c  STRIPE matrix.
c
      title = 'STRIPE'
      n = 5
      call stripe ( n, a )
      call r8mat_determinant ( n, a, determ2 )
      norm_frobenius = r8mat_norm_fro ( n, n, a )
      write ( *, '(2x,a20,2x,i4,2x,14x,2x,g14.6,2x,g14.6)' ) 
     &  title, n,          determ2, norm_frobenius
c
c  SUMMATION matrix.
c
      title = 'SUMMATION'
      n = 5
      call summation ( n, n, a )
      call summation_determinant ( n, determ1 )
      call r8mat_determinant ( n, a, determ2 )
      norm_frobenius = r8mat_norm_fro ( n, n, a )
      write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) 
     &  title, n, determ1, determ2, norm_frobenius
c
c  SWEET1 matrix.
c
      title = 'SWEET1'
      n = 6
      perturb = 0.0D+00
      call sweet1 ( perturb, a )
      call r8mat_determinant ( n, a, determ2 )
      norm_frobenius = r8mat_norm_fro ( n, n, a )
      write ( *, '(2x,a20,2x,i4,2x,14x,2x,g14.6,2x,g14.6)' ) 
     &  title, n,          determ2, norm_frobenius
c
c  SWEET2 matrix.
c
      title = 'SWEET2'
      n = 6
      perturb = 0.0D+00
      call sweet2 ( perturb, a )
      call r8mat_determinant ( n, a, determ2 )
      norm_frobenius = r8mat_norm_fro ( n, n, a )
      write ( *, '(2x,a20,2x,i4,2x,14x,2x,g14.6,2x,g14.6)' ) 
     &  title, n,          determ2, norm_frobenius
c
c  SWEET3 matrix.
c
      title = 'SWEET3'
      n = 6
      perturb = 0.0D+00
      call sweet3 ( perturb, a )
      call r8mat_determinant ( n, a, determ2 )
      norm_frobenius = r8mat_norm_fro ( n, n, a )
      write ( *, '(2x,a20,2x,i4,2x,14x,2x,g14.6,2x,g14.6)' ) 
     &  title, n,          determ2, norm_frobenius
c
c  SWEET4 matrix.
c
      title = 'SWEET4'
      n = 13
      call sweet4 ( perturb, a )
      call r8mat_determinant ( n, a, determ2 )
      norm_frobenius = r8mat_norm_fro ( n, n, a )
      write ( *, '(2x,a20,2x,i4,2x,14x,2x,g14.6,2x,g14.6)' ) 
     &  title, n,          determ2, norm_frobenius
c
c  SYLVESTER matrix.
c
      title = 'SYLVESTER'
      n = 5
      x_n = 3 + 1
      y_n = 2 + 1
      seed = 123456789
      call r8vec_uniform_01 ( x_n, seed, x )
      do i = 1, x_n
        x(i) = anint ( 50.0D+00 * x(i) - 25.0D+00 ) / 5.0D+00
      end do
      call r8vec_uniform_01 ( y_n, seed, y )
      do i = 1, y_n
        y(i) = anint ( 50.0D+00 * y(i) - 25.0D+00 ) / 5.0D+00
      end do
      call sylvester ( n, x_n - 1, x, y_n - 1, y, a )
      call r8mat_determinant ( n, a, determ2 )
      norm_frobenius = r8mat_norm_fro ( n, n, a )
      write ( *, '(2x,a20,2x,i4,2x,14x,2x,g14.6,2x,g14.6)' ) 
     &  title, n,          determ2, norm_frobenius
c
c  SYMM_RANDOM matrix.
c
      title = 'SYMM_RANDOM'
      n = 5
      seed = 123456789
      call r8vec_uniform_01 ( n, seed, x )
      do i = 1, n
        x(i) = anint ( 50.0D+00 * x(i) - 25.0D+00 ) / 5.0D+00
      end do
      call symm_random ( n, x, seed, a )
      call symm_random_determinant ( n, x, determ1 )
      call r8mat_determinant ( n, a, determ2 )
      norm_frobenius = r8mat_norm_fro ( n, n, a )
      write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) 
     &  title, n, determ1, determ2, norm_frobenius
c
c  TOEPLITZ matrix.
c
      title = 'TOEPLITZ'
      n = 5
      seed = 123456789
      call r8vec_uniform_01 ( 2 * n - 1, seed, x )
      do i = 1, 2 * n - 1
        x(i) = anint ( 50.0D+00 * x(i) - 25.0D+00 ) / 5.0D+00
      end do
      call toeplitz ( n, x, a )
      call r8mat_determinant ( n, a, determ2 )
      norm_frobenius = r8mat_norm_fro ( n, n, a )
      write ( *, '(2x,a20,2x,i4,2x,14x,2x,g14.6,2x,g14.6)' ) 
     &  title, n,          determ2, norm_frobenius
c
c  TOEPLITZ_5DIAG matrix.
c
      title = 'TOEPLITZ_5DIAG'
      n = 5
      seed = 123456789
      d1 = r8_uniform_01 ( seed )
      d1 = anint ( 50.0D+00 * d1 - 25.0D+00 ) / 5.0D+00
      d2 = r8_uniform_01 ( seed )
      d2 = anint ( 50.0D+00 * d2 - 25.0D+00 ) / 5.0D+00
      d3 = r8_uniform_01 ( seed )
      d3 = anint ( 50.0D+00 * d3 - 25.0D+00 ) / 5.0D+00
      d4 = r8_uniform_01 ( seed )
      d4 = anint ( 50.0D+00 * d4 - 25.0D+00 ) / 5.0D+00
      d5 = r8_uniform_01 ( seed )
      d5 = anint ( 50.0D+00 * d5 - 25.0D+00 ) / 5.0D+00
      call toeplitz_5diag ( n, d1, d2, d3, d4, d5, a )
      call r8mat_determinant ( n, a, determ2 )
      norm_frobenius = r8mat_norm_fro ( n, n, a )
      write ( *, '(2x,a20,2x,i4,2x,14x,2x,g14.6,2x,g14.6)' ) 
     &  title, n,          determ2, norm_frobenius
c
c  TOEPLITZ_5S matrix.
c
      title = 'TOEPLITZ_5S'
      row_num = 5
      col_num = 5
      n = row_num * col_num
      seed = 123456789
      alpha = r8_uniform_01 ( seed )
      alpha = anint ( 50.0D+00 * alpha - 25.0D+00 ) / 5.0D+00
      beta = r8_uniform_01 ( seed )
      beta = anint ( 50.0D+00 * beta - 25.0D+00 ) / 5.0D+00
      gamma = r8_uniform_01 ( seed )
      gamma = anint ( 50.0D+00 * gamma - 25.0D+00 ) / 5.0D+00
      call toeplitz_5s ( row_num, col_num, alpha, beta, gamma, n, a )
      call r8mat_determinant ( n, a, determ2 )
      norm_frobenius = r8mat_norm_fro ( n, n, a )
      write ( *, '(2x,a20,2x,i4,2x,14x,2x,g14.6,2x,g14.6)' ) 
     &  title, n,          determ2, norm_frobenius
c
c  TOEPLITZ_PDS matrix.
c
      title = 'TOEPLITZ_PDS'
      m = 3
      n = 5
      seed = 123456789
      call r8vec_uniform_01 ( m, seed, x )
      call r8vec_uniform_01 ( m, seed, y )
      y_sum = r8vec_sum ( m, y )
      do i = 1, m
        y(i) = y(i) / y_sum
      end do
      call toeplitz_pds ( m, n, x, y, a )
      call r8mat_determinant ( n, a, determ2 )
      norm_frobenius = r8mat_norm_fro ( n, n, a )
      write ( *, '(2x,a20,2x,i4,2x,14x,2x,g14.6,2x,g14.6)' ) 
     &  title, n,          determ2, norm_frobenius
c
c  TOURNAMENT_RANDOM matrix.
c
      title = 'TOURNAMENT_RANDOM'
      n = 5
      seed_save = 123456789
      seed = seed_save
      call tournament_random ( n, seed, a )
      seed = seed_save
      call tournament_random_determinant ( n, seed, determ1 )
      call r8mat_determinant ( n, a, determ2 )
      norm_frobenius = r8mat_norm_fro ( n, n, a )
      write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) 
     &  title, n, determ1, determ2, norm_frobenius
c
c  TRANSITION_RANDOM matrix.
c
      title = 'TRANSITION_RANDOM'
      n = 5
      seed = 123456789
      call transition_random ( n, seed, a )
      call r8mat_determinant ( n, a, determ2 )
      norm_frobenius = r8mat_norm_fro ( n, n, a )
      write ( *, '(2x,a20,2x,i4,2x,14x,2x,g14.6,2x,g14.6)' ) 
     &  title, n,          determ2, norm_frobenius
c
c  TRENCH matrix.
c
      title = 'TRENCH'
      n = 5
      seed = 123456789
      alpha = r8_uniform_01 ( seed )
      alpha = anint ( 50.0D+00 * alpha - 25.0D+00 ) / 5.0D+00
      call trench ( alpha, n, n, a )
      call r8mat_determinant ( n, a, determ2 )
      norm_frobenius = r8mat_norm_fro ( n, n, a )
      write ( *, '(2x,a20,2x,i4,2x,14x,2x,g14.6,2x,g14.6)' ) 
     &  title, n,          determ2, norm_frobenius
c
c  TRI_UPPER matrix.
c
      title = 'TRI_UPPER'
      n = 5
      seed = 123456789
      alpha = r8_uniform_01 ( seed )
      alpha = anint ( 50.0D+00 * alpha - 25.0D+00 ) / 5.0D+00
      call tri_upper ( alpha, n, a )
      call tri_upper_determinant ( alpha, n, determ1 )
      call r8mat_determinant ( n, a, determ2 )
      norm_frobenius = r8mat_norm_fro ( n, n, a )
      write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) 
     &  title, n, determ1, determ2, norm_frobenius
c
c  TRIS matrix.
c
      title = 'TRIS'
      n = 5
      seed = 123456789
      alpha = r8_uniform_01 ( seed )
      alpha = anint ( 50.0D+00 * alpha - 25.0D+00 ) / 5.0D+00
      beta = r8_uniform_01 ( seed )
      beta = anint ( 50.0D+00 * beta - 25.0D+00 ) / 5.0D+00
      gamma = r8_uniform_01 ( seed )
      gamma = anint ( 50.0D+00 * gamma - 25.0D+00 ) / 5.0D+00
      call tris ( n, n, alpha, beta, gamma, a )
      call tris_determinant ( n, alpha, beta, gamma, determ1 )
      call r8mat_determinant ( n, a, determ2 )
      norm_frobenius = r8mat_norm_fro ( n, n, a )
      write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) 
     &  title, n, determ1, determ2, norm_frobenius
c
c  TRIV matrix.
c
      title = 'TRIV'
      n = 5
      seed = 123456789
      call r8vec_uniform_01 ( n - 1, seed, x )
      do i = 1, n - 1
        x(i) = anint ( 50.0D+00 * x(i) - 25.0D+00 ) / 5.0D+00
      end do
      call r8vec_uniform_01 ( n, seed, y )
      do i = 1, n
        y(i) = anint ( 50.0D+00 * y(i) - 25.0D+00 ) / 5.0D+00
      end do
      call r8vec_uniform_01 ( n - 1, seed, z )
      do i = 1, n - 1
        z(i) = anint ( 50.0D+00 * z(i) - 25.0D+00 ) / 5.0D+00
      end do
      call triv ( n, x, y, z, a )
      call triv_determinant ( n, x, y, z, determ1 )
      call r8mat_determinant ( n, a, determ2 )
      norm_frobenius = r8mat_norm_fro ( n, n, a )
      write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) 
     &  title, n, determ1, determ2, norm_frobenius

c
c  TRIW matrix.
c
      title = 'TRIW'
      n = 5
      seed = 123456789
      k = i4_uniform ( 0, n-1, seed )
      alpha = r8_uniform_01 ( seed )
      alpha = anint ( 50.0D+00 * alpha - 25.0D+00 ) / 5.0D+00
      call triw ( alpha, k, n, a )
      call triw_determinant ( alpha, k, n, determ1 )
      call r8mat_determinant ( n, a, determ2 )
      norm_frobenius = r8mat_norm_fro ( n, n, a )
      write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) 
     &  title, n, determ1, determ2, norm_frobenius
c
c  UPSHIFT matrix.
c
      title = 'UPSHIFT'
      n = 5
      call upshift ( n, a )
      call upshift_determinant ( n, determ1 )
      call r8mat_determinant ( n, a, determ2 )
      norm_frobenius = r8mat_norm_fro ( n, n, a )
      write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) 
     &  title, n, determ1, determ2, norm_frobenius
c
c  VAND1 matrix.
c
      title = 'VAND1'
      n = 5
      seed = 123456789
      call r8vec_uniform_01 ( n, seed, x )
      do i = 1, n
        x(i) = anint ( 50.0D+00 * x(i) - 25.0D+00 ) / 5.0D+00
      end do
      call vand1 ( n, x, a )
      call vand1_determinant ( n, x, determ1 )
      call r8mat_determinant ( n, a, determ2 )
      norm_frobenius = r8mat_norm_fro ( n, n, a )
      write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) 
     &  title, n, determ1, determ2, norm_frobenius
c
c  VAND2 matrix.
c
      title = 'VAND2'
      n = 5
      seed = 123456789
      call r8vec_uniform_01 ( n, seed, x )
      do i = 1, n
        x(i) = anint ( 50.0D+00 * x(i) - 25.0D+00 ) / 5.0D+00
      end do
      call vand2 ( n, x, a )
      call vand2_determinant ( n, x, determ1 )
      call r8mat_determinant ( n, a, determ2 )
      norm_frobenius = r8mat_norm_fro ( n, n, a )
      write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) 
     &  title, n, determ1, determ2, norm_frobenius
c
c  WATHEN matrix.
c
      title = 'WATHEN'
      row_num = 5
      col_num = 5
      call wathen_order ( row_num, col_num, n )
      call wathen ( row_num, col_num, n, a )
      call r8mat_determinant ( n, a, determ2 )
      call r8mat_determinant ( n, a, determ2 )
      norm_frobenius = r8mat_norm_fro ( n, n, a )
      write ( *, '(2x,a20,2x,i4,2x,14x,2x,g14.6,2x,g14.6)' ) 
     &  title, n,          determ2, norm_frobenius
c
c  WILK03 matrix.
c
      title = 'WILK03'
      n = 3
      call wilk03 ( a )
      call wilk03_determinant ( determ1 )
      call r8mat_determinant ( n, a, determ2 )
      norm_frobenius = r8mat_norm_fro ( n, n, a )
      write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) 
     &  title, n, determ1, determ2, norm_frobenius
c
c  WILK04 matrix.
c
      title = 'WILK04'
      n = 4
      call wilk04 ( a )
      call wilk04_determinant ( determ1 )
      call r8mat_determinant ( n, a, determ2 )
      norm_frobenius = r8mat_norm_fro ( n, n, a )
      write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) 
     &  title, n, determ1, determ2, norm_frobenius
c
c  WILK05 matrix.
c
      title = 'WILK05'
      n = 5
      call wilk05 ( a )
      call wilk05_determinant ( determ1 )
      call r8mat_determinant ( n, a, determ2 )
      norm_frobenius = r8mat_norm_fro ( n, n, a )
      write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) 
     &  title, n, determ1, determ2, norm_frobenius
c
c  WILK12 matrix.
c
      title = 'WILK12'
      n = 12
      call wilk12 ( a )
      call wilk12_determinant ( determ1 )
      call r8mat_determinant ( n, a, determ2 )
      norm_frobenius = r8mat_norm_fro ( n, n, a )
      write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) 
     &  title, n, determ1, determ2, norm_frobenius
c
c  WILK20 matrix.
c
      title = 'WILK20'
      n = 20
      seed = 123456789
      alpha = r8_uniform_01 ( seed )
      alpha = anint ( 50.0D+00 * alpha - 25.0D+00 ) / 5.0D+00
      call wilk20 ( alpha, a )
      call r8mat_determinant ( n, a, determ2 )
      norm_frobenius = r8mat_norm_fro ( n, n, a )
      write ( *, '(2x,a20,2x,i4,2x,14x,2x,g14.6,2x,g14.6)' )
     &  title, n,          determ2, norm_frobenius
c
c  WILK21 matrix.
c
      title = 'WILK21'
      n = 21
      call wilk21 ( n, a )
      call wilk21_determinant ( n, determ1 )
      call r8mat_determinant ( n, a, determ2 )
      norm_frobenius = r8mat_norm_fro ( n, n, a )
      write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) 
     &  title, n, determ1, determ2, norm_frobenius
c
c  WILSON matrix.
c
      title = 'WILSON'
      n = 4
      call wilson ( a )
      call wilson_determinant ( n, determ1 )
      call r8mat_determinant ( n, a, determ2 )
      norm_frobenius = r8mat_norm_fro ( n, n, a )
      write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) 
     &  title, n, determ1, determ2, norm_frobenius
c
c  ZERO matrix.
c
      title = 'ZERO'
      n = 5
      call zero ( n, n, a )
      call zero_determinant ( n, determ1 )
      call r8mat_determinant ( n, a, determ2 )
      norm_frobenius = r8mat_norm_fro ( n, n, a )
      write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) 
     &  title, n, determ1, determ2, norm_frobenius
c
c  ZIELKE matrix.
c
      title = 'ZIELKE'
      n = 5
      seed = 123456789
      d1 = r8_uniform_01 ( seed )
      d1 = anint ( 50.0D+00 * d1 - 25.0D+00 ) / 5.0D+00
      d2 = r8_uniform_01 ( seed )
      d2 = anint ( 50.0D+00 * d2 - 25.0D+00 ) / 5.0D+00
      d3 = r8_uniform_01 ( seed )
      d3 = anint ( 50.0D+00 * d3 - 25.0D+00 ) / 5.0D+00
      call zielke ( n, d1, d2, d3, a )
      call r8mat_determinant ( n, a, determ2 )
      norm_frobenius = r8mat_norm_fro ( n, n, a )
      write ( *, '(2x,a20,2x,i4,2x,14x,2x,g14.6,2x,g14.6)' )
     &  title, n,          determ2, norm_frobenius

      return
      end
      subroutine test_eigen ( )

c*********************************************************************72
c
cc TEST_EIGEN tests the eigenvalue computations.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    03 July 2011
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n_max
      parameter ( n_max = 25 )

      double precision a(n_max,n_max)
      double precision alpha
      double precision beta
      double precision error_frobenius
      double precision gamma
      integer i
      integer i1
      integer i4_uniform
      integer k
      double precision lambda(n_max)
      integer n
      double precision norm_frobenius
      double precision r8_uniform_01
      double precision r8mat_norm_fro
      integer rank
      integer seed
      integer seed_save
      character * ( 20 ) title
      double precision v1(n_max)
      double precision v2(n_max)
      double precision v3(n_max)
      double precision x(n_max,n_max)

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST_EIGEN'
      write ( *, '(a)' ) 
     &  '  Compute the Frobenius norm of the eigenvalue error:'
      write ( *, '(a)' ) '    A * X - X * LAMBDA'
      write ( *, '(a)' ) 
     &  '  given a set of K eigenvectors X and eigenvalues LAMBDA.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 
     &  '  Matrix title             N     K      ||A||' 
     &  // '          ||(A-Lambda*I)*X||'
      write ( *, '(a)' ) ' '
c
c  BODEWIG matrix.
c
      title = 'BODEWIG'
      n = 4
      k = 4
      call bodewig ( a )
      call bodewig_eigenvalues ( lambda )
      call bodewig_right ( x )
      call r8mat_is_eigen_right ( n, k, a, x, lambda, error_frobenius )
      norm_frobenius = r8mat_norm_fro ( n, n, a )
      write ( *, '(2x,a20,2x,i4,2x,i4,2x,g14.6,2x,g14.6)' ) 
     &  title, n, k, norm_frobenius, error_frobenius
c
c  CARRY matrix.
c
      title = 'CARRY'
      n = 5
      k = 5
      seed = 123456789
      i1 = i4_uniform ( 2, 20, seed )
      call carry ( i1, n, a )
      call carry_eigenvalues ( i1, n, lambda )
      call carry_right ( n, x )
      call r8mat_is_eigen_right ( n, k, a, x, lambda, error_frobenius )
      norm_frobenius = r8mat_norm_fro ( n, n, a )
      write ( *, '(2x,a20,2x,i4,2x,i4,2x,g14.6,2x,g14.6)' ) 
     &  title, n, k, norm_frobenius, error_frobenius
c
c  CHOW matrix.
c
      title = 'CHOW'
      n = 5
      k = 5
      seed = 123456789
      alpha = r8_uniform_01 ( seed )
      alpha = anint ( 50.0D+00 * alpha ) / 5.0D+00
      beta = r8_uniform_01 ( seed )
      beta = anint ( 50.0D+00 * beta ) / 5.0D+00
      call chow ( alpha, beta, n, n, a )
      call chow_eigenvalues ( alpha, beta, n, lambda )
      call chow_right ( alpha, beta, n, x )
      call r8mat_is_eigen_right ( n, k, a, x, lambda, error_frobenius )
      norm_frobenius = r8mat_norm_fro ( n, n, a )
      write ( *, '(2x,a20,2x,i4,2x,i4,2x,g14.6,2x,g14.6)' ) 
     &  title, n, k, norm_frobenius, error_frobenius
c
c  COMBIN matrix.
c
      title = 'COMBIN'
      n = 5
      k = 5
      seed = 123456789
      alpha = r8_uniform_01 ( seed )
      alpha = anint ( 50.0D+00 * alpha ) / 5.0D+00
      beta = r8_uniform_01 ( seed )
      beta = anint ( 50.0D+00 * beta ) / 5.0D+00
      call combin ( alpha, beta, n, a )
      call combin_eigenvalues ( alpha, beta, n, lambda )
      call combin_right ( alpha, beta, n, x )
      call r8mat_is_eigen_right ( n, k, a, x, lambda, error_frobenius )
      norm_frobenius = r8mat_norm_fro ( n, n, a )
      write ( *, '(2x,a20,2x,i4,2x,i4,2x,g14.6,2x,g14.6)' ) 
     &  title, n, k, norm_frobenius, error_frobenius
c
c  DIF2 matrix.
c
      title = 'DIF2'
      n = 5
      k = 5
      call dif2 ( n, n, a )
      call dif2_eigenvalues ( n, lambda )
      call dif2_right ( n, x )
      call r8mat_is_eigen_right ( n, k, a, x, lambda, error_frobenius )
      norm_frobenius = r8mat_norm_fro ( n, n, a )
      write ( *, '(2x,a20,2x,i4,2x,i4,2x,g14.6,2x,g14.6)' ) 
     &  title, n, k, norm_frobenius, error_frobenius
c
c  EXCHANGE matrix.
c
      title = 'EXCHANGE'
      n = 5
      k = 5
      call exchange ( n, n, a )
      call exchange_eigenvalues ( n, lambda )
      call exchange_right ( n, x )
      call r8mat_is_eigen_right ( n, k, a, x, lambda, error_frobenius )
      norm_frobenius = r8mat_norm_fro ( n, n, a )
      write ( *, '(2x,a20,2x,i4,2x,i4,2x,g14.6,2x,g14.6)' ) 
     &  title, n, k, norm_frobenius, error_frobenius
c
c  IDEM_RANDOM matrix.
c
      title = 'IDEM_RANDOM'
      n = 5
      k = 5
      rank = 3
      seed_save = 987654321
      seed = seed_save
      call idem_random ( n, rank, seed, a )
      call idem_random_eigenvalues ( n, rank, lambda )
      seed = seed_save
      call idem_random_right ( n, rank, seed, x )
      call r8mat_is_eigen_right ( n, k, a, x, lambda, error_frobenius )
      norm_frobenius = r8mat_norm_fro ( n, n, a )
      write ( *, '(2x,a20,2x,i4,2x,i4,2x,g14.6,2x,g14.6)' ) 
     &  title, n, k, norm_frobenius, error_frobenius
c
c  IDENTITY matrix.
c
      title = 'IDENTITY'
      n = 5
      k = 5
      call identity ( n, n, a )
      call identity_eigenvalues ( n, lambda )
      call identity_right ( n, x )
      call r8mat_is_eigen_right ( n, k, a, x, lambda, error_frobenius )
      norm_frobenius = r8mat_norm_fro ( n, n, a )
      write ( *, '(2x,a20,2x,i4,2x,i4,2x,g14.6,2x,g14.6)' ) 
     &  title, n, k, norm_frobenius, error_frobenius
c
c  ILL3 matrix.
c
      title = 'ILL3'
      n = 3
      k = 3
      call ill3 ( a )
      call ill3_eigenvalues ( lambda )
      call ill3_right ( x )
      call r8mat_is_eigen_right ( n, k, a, x, lambda, error_frobenius )
      norm_frobenius = r8mat_norm_fro ( n, n, a )
      write ( *, '(2x,a20,2x,i4,2x,i4,2x,g14.6,2x,g14.6)' ) 
     &  title, n, k, norm_frobenius, error_frobenius
c
c  KMS matrix.
c
      title = 'KMS'
      n = 5
      k = 5
      seed = 123456789
      alpha = r8_uniform_01 ( seed )
      call kms ( alpha, n, n, a )
      call kms_eigenvalues ( alpha, n, lambda )
      call kms_right ( alpha, n, x )
      call r8mat_is_eigen_right ( n, k, a, x, lambda, error_frobenius )
      norm_frobenius = r8mat_norm_fro ( n, n, a )
      write ( *, '(2x,a20,2x,i4,2x,i4,2x,g14.6,2x,g14.6)' ) 
     &  title, n, k, norm_frobenius, error_frobenius
c
c  ONE matrix.
c
      title = 'ONE'
      n = 5
      k = 5
      call one ( n, n, a )
      call one_eigenvalues ( n, lambda )
      call one_right ( n, x )
      call r8mat_is_eigen_right ( n, k, a, x, lambda, error_frobenius )
      norm_frobenius = r8mat_norm_fro ( n, n, a )
      write ( *, '(2x,a20,2x,i4,2x,i4,2x,g14.6,2x,g14.6)' ) 
     &  title, n, k, norm_frobenius, error_frobenius
c
c  ORTEGA matrix.
c
      title = 'ORTEGA'
      n = 5
      k = n
      seed = 123456789
      call r8vec_uniform_01 ( n, seed, v1 )
      do i = 1, n
        v1(i) = anint ( 50.0D+00 * v1(i) - 25.0D+00  ) / 5.0D+00
      end do
      call r8vec_uniform_01 ( n, seed, v2 )
      do i = 1, n
        v2(i) = anint ( 50.0D+00 * v2(i) - 25.0D+00  ) / 5.0D+00
      end do
      call r8vec_uniform_01 ( n, seed, v3 )
      do i = 1, n
        v3(i) = anint ( 50.0D+00 * v3(i) - 25.0D+00  ) / 5.0D+00
      end do
      call ortega ( n, v1, v2, v3, a )
      call ortega_eigenvalues ( n, v1, v2, v3, lambda )
      call ortega_right ( n, v1, v2, v3, x )
      call r8mat_is_eigen_right ( n, k, a, x, lambda, error_frobenius )
      norm_frobenius = r8mat_norm_fro ( n, n, a )
      write ( *, '(2x,a20,2x,i4,2x,i4,2x,g14.6,2x,g14.6)' ) 
     &  title, n, k, norm_frobenius, error_frobenius
c
c  OTO matrix.
c
      title = 'OTO'
      n = 5
      k = 5
      call oto ( n, n, a )
      call oto_eigenvalues ( n, lambda )
      call oto_right ( n, x )
      call r8mat_is_eigen_right ( n, k, a, x, lambda, error_frobenius )
      norm_frobenius = r8mat_norm_fro ( n, n, a )
      write ( *, '(2x,a20,2x,i4,2x,i4,2x,g14.6,2x,g14.6)' ) 
     &  title, n, k, norm_frobenius, error_frobenius
c
c  PDS_RANDOM matrix.
c
      title = 'PDS_RANDOM'
      n = 5
      k = 5
      seed_save = 123456789
      seed = seed_save
      call pds_random ( n, seed, a )
      seed = seed_save
      call pds_random_eigenvalues ( n, seed, lambda )
      seed = seed_save
      call pds_random_right ( n, seed, x )
      call r8mat_is_eigen_right ( n, k, a, x, lambda, error_frobenius )
      norm_frobenius = r8mat_norm_fro ( n, n, a )
      write ( *, '(2x,a20,2x,i4,2x,i4,2x,g14.6,2x,g14.6)' ) 
     &  title, n, k, norm_frobenius, error_frobenius
c
c  PEI matrix.
c
      title = 'PEI'
      n = 5
      k = 5
      seed = 123456789
      alpha = r8_uniform_01 ( seed )
      alpha = anint ( 50.0D+00 * alpha ) / 5.0D+00
      call pei ( alpha, n, a )
      call pei_eigenvalues ( alpha, n, lambda )
      call pei_right ( alpha, n, x )
      call r8mat_is_eigen_right ( n, k, a, x, lambda, error_frobenius )
      norm_frobenius = r8mat_norm_fro ( n, n, a )
      write ( *, '(2x,a20,2x,i4,2x,i4,2x,g14.6,2x,g14.6)' ) 
     &  title, n, k, norm_frobenius, error_frobenius
c
c  RODMAN matrix.
c
      title = 'RODMAN'
      n = 5
      k = 5
      seed = 123456789
      alpha = r8_uniform_01 ( seed )
      alpha = anint ( 50.0D+00 * alpha ) / 5.0D+00
      call rodman ( alpha, n, n, a )
      call rodman_eigenvalues ( alpha, n, lambda )
      call rodman_right ( alpha, n, x )
      call r8mat_is_eigen_right ( n, k, a, x, lambda, error_frobenius )
      norm_frobenius = r8mat_norm_fro ( n, n, a )
      write ( *, '(2x,a20,2x,i4,2x,i4,2x,g14.6,2x,g14.6)' ) 
     &  title, n, k, norm_frobenius, error_frobenius
c
c  ROSSER1 matrix.
c
      title = 'ROSSER1'
      n = 8
      k = 8
      call rosser1 ( a )
      call rosser1_eigenvalues ( lambda )
      call rosser1_right ( x )
      call r8mat_is_eigen_right ( n, k, a, x, lambda, error_frobenius )
      norm_frobenius = r8mat_norm_fro ( n, n, a )
      write ( *, '(2x,a20,2x,i4,2x,i4,2x,g14.6,2x,g14.6)' ) 
     &  title, n, k, norm_frobenius, error_frobenius
c
c  RUTIS1 matrix.
c
      title = 'RUTIS1'
      n = 4
      k = 4
      call rutis1 ( a )
      call rutis1_eigenvalues ( lambda )
      call rutis1_right ( x )
      call r8mat_is_eigen_right ( n, k, a, x, lambda, error_frobenius )
      norm_frobenius = r8mat_norm_fro ( n, n, a )
      write ( *, '(2x,a20,2x,i4,2x,i4,2x,g14.6,2x,g14.6)' ) 
     &  title, n, k, norm_frobenius, error_frobenius
c
c  RUTIS2 matrix.
c
      title = 'RUTIS2'
      n = 4
      k = 4
      call rutis2 ( a )
      call rutis2_eigenvalues ( lambda )
      call rutis2_right ( x )
      call r8mat_is_eigen_right ( n, k, a, x, lambda, error_frobenius )
      norm_frobenius = r8mat_norm_fro ( n, n, a )
      write ( *, '(2x,a20,2x,i4,2x,i4,2x,g14.6,2x,g14.6)' ) 
     &  title, n, k, norm_frobenius, error_frobenius
c
c  RUTIS3 matrix.
c  COMPLEX eigenvalues cannot be handled yetc
c
      if ( .false. ) then
        title = 'RUTIS3'
        n = 4
        k = 4
        call rutis3 ( a )
        call rutis3_eigenvalues ( lambda )
        call rutis3_right ( x )
c       call r8mat_is_eigen_right ( n, k, a, x, lambda, error_frobenius )
        norm_frobenius = r8mat_norm_fro ( n, n, a )
        write ( *, '(2x,a20,2x,i4,2x,i4,2x,g14.6,2x,g14.6)' ) 
     &    title, n, k, norm_frobenius, error_frobenius

      end if
c
c  RUTIS5 matrix.
c
      title = 'RUTIS5'
      n = 4
      k = 4
      call rutis5 ( a )
      call rutis5_eigenvalues ( lambda )
      call rutis5_right ( x )
      call r8mat_is_eigen_right ( n, k, a, x, lambda, error_frobenius )
      norm_frobenius = r8mat_norm_fro ( n, n, a )
      write ( *, '(2x,a20,2x,i4,2x,i4,2x,g14.6,2x,g14.6)' ) 
     &  title, n, k, norm_frobenius, error_frobenius
c
c  WILK12 matrix.
c
      title = 'WILK12'
      n = 12
      k = 12
      call wilk12 ( a )
      call wilk12_eigenvalues ( lambda )
      call wilk12_right ( x )
      call r8mat_is_eigen_right ( n, k, a, x, lambda, error_frobenius )
      norm_frobenius = r8mat_norm_fro ( n, n, a )
      write ( *, '(2x,a20,2x,i4,2x,i4,2x,g14.6,2x,g14.6)' ) 
     &  title, n, k, norm_frobenius, error_frobenius
c
c  WILSON matrix.
c
      title = 'WILSON'
      n = 4
      k = 4
      call wilson ( a )
      call wilson_eigenvalues ( lambda )
      call wilson_right ( x )
      call r8mat_is_eigen_right ( n, k, a, x, lambda, error_frobenius )
      norm_frobenius = r8mat_norm_fro ( n, n, a )
      write ( *, '(2x,a20,2x,i4,2x,i4,2x,g14.6,2x,g14.6)' ) 
     &  title, n, k, norm_frobenius, error_frobenius
c
c  ZERO matrix.
c
      title = 'ZERO'
      n = 5
      k = 5
      call zero ( n, n, a )
      call zero_eigenvalues ( n, lambda )
      call zero_right ( n, x )
      call r8mat_is_eigen_right ( n, k, a, x, lambda, error_frobenius )
      norm_frobenius = r8mat_norm_fro ( n, n, a )
      write ( *, '(2x,a20,2x,i4,2x,i4,2x,g14.6,2x,g14.6)' ) 
     &  title, n, k, norm_frobenius, error_frobenius

      return
      end
      subroutine test_inverse ( )

c*********************************************************************72
c
cc TEST_INVERSE tests the inverse computations.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    10 July 2011
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n_max
      parameter ( n_max = 25 )

      double precision a(n_max,n_max)
      double precision alpha
      double precision b(n_max,n_max)
      double precision beta
      double precision c(n_max,n_max)
      double precision error_ab
      double precision error_ac
      double precision gamma
      integer i
      integer i4_uniform
      integer k
      double precision l(n_max,n_max)
      integer n
      double precision norm_frobenius
      double precision p(n_max,n_max)
      integer pivot(n_max)
      double precision r8_uniform_01
      double precision r8mat_norm_fro
      integer seed
      integer seed_save
      character * ( 20 ) title
      double precision u(n_max,n_max)
      double precision v1(n_max)
      double precision v2(n_max)
      double precision v3(n_max)
      double precision x(n_max)
      integer x_n
      double precision y(n_max)
      integer y_n
      double precision z(n_max)

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST_INVERSE'
      write ( *, '(a)' ) '  A = a test matrix of order N;'
      write ( *, '(a)' ) '  B = inverse as computed by a routine.'
      write ( *, '(a)' ) '  C = inverse as computed by R8MAT_INVERSE.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  ||A||    = Frobenius norm of A.'
      write ( *, '(a)' ) '  ||I-AB|| = Frobenius norm of I-A*B.'
      write ( *, '(a)' ) '  ||I-AC|| = Frobenius norm of I-A*C.'

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Matrix title             N        ' // 
     &               '||A||        ||I-AB||        ||I-AC||'
      write ( *, '(a)' ) ' '
c
c  AEGERTER matrix.
c
      title = 'AEGERTER'
      n = 5
      call aegerter ( n, a )
      call aegerter_inverse ( n, b )
      call r8mat_inverse ( n, a, c )
      call r8mat_is_inverse ( n, a, b, error_ab )
      call r8mat_is_inverse ( n, a, c, error_ac )
      norm_frobenius = r8mat_norm_fro ( n, n, a )
      write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) 
     &  title, n, norm_frobenius, error_ab, error_ac
c
c  BAB matrix.
c
      title = 'BAB'
      n = 5
      seed = 123456789
      alpha = r8_uniform_01 ( seed )
      alpha = anint ( 50.0D+00 * alpha - 25.0D+00 ) / 5.0D+00
      beta = r8_uniform_01 ( seed )
      beta = anint ( 50.0D+00 * beta - 25.0D+00 ) / 5.0D+00
      call bab ( n, alpha, beta, a )
      call bab_inverse ( n, alpha, beta, b )
      call r8mat_inverse ( n, a, c )
      call r8mat_is_inverse ( n, a, b, error_ab )
      call r8mat_is_inverse ( n, a, c, error_ac )
      norm_frobenius = r8mat_norm_fro ( n, n, a )
      write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) 
     &  title, n, norm_frobenius, error_ab, error_ac
c
c  BERNSTEIN matrix.
c
      title = 'BERNSTEIN'
      n = 5
      call bernstein ( n, a )
      call bernstein_inverse ( n, b )
      call r8mat_inverse ( n, a, c )
      call r8mat_is_inverse ( n, a, b, error_ab )
      call r8mat_is_inverse ( n, a, c, error_ac )
      norm_frobenius = r8mat_norm_fro ( n, n, a )
      write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) 
     &  title, n, norm_frobenius, error_ab, error_ac
c
c  BIS matrix.
c
      title = 'BIS'
      n = 5
      seed = 123456789
      alpha = r8_uniform_01 ( seed )
      alpha = anint ( 50.0D+00 * alpha ) / 5.0D+00
      beta = r8_uniform_01 ( seed )
      beta = anint ( 50.0D+00 * beta ) / 5.0D+00
      call bis ( alpha, beta, n, n, a );
      call bis_inverse ( alpha, beta, n, b )
      call r8mat_inverse ( n, a, c )
      call r8mat_is_inverse ( n, a, b, error_ab )
      call r8mat_is_inverse ( n, a, c, error_ac )
      norm_frobenius = r8mat_norm_fro ( n, n, a )
      write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) 
     &  title, n, norm_frobenius, error_ab, error_ac
c
c  BODEWIG matrix.
c
      title = 'BODEWIG'
      n = 4
      call bodewig ( a )
      call bodewig_inverse ( b )
      call r8mat_inverse ( n, a, c )
      call r8mat_is_inverse ( n, a, b, error_ab )
      call r8mat_is_inverse ( n, a, c, error_ac )
      norm_frobenius = r8mat_norm_fro ( n, n, a )
      write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) 
     &  title, n, norm_frobenius, error_ab, error_ac
c
c  BOOTHROYD matrix.
c
      title = 'BOOTHROYD'
      n = 5
      call boothroyd ( n, a )
      call boothroyd_inverse ( n, b )
      call r8mat_inverse ( n, a, c )
      call r8mat_is_inverse ( n, a, b, error_ab )
      call r8mat_is_inverse ( n, a, c, error_ac )
      norm_frobenius = r8mat_norm_fro ( n, n, a )
      write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) 
     &  title, n, norm_frobenius, error_ab, error_ac
c
c  BORDERBAND matrix.
c
      title = 'BORDERBAND'
      n = 5
      call borderband ( n, a )
      call borderband_inverse ( n, b )
      call r8mat_inverse ( n, a, c )
      call r8mat_is_inverse ( n, a, b, error_ab )
      call r8mat_is_inverse ( n, a, c, error_ac )
      norm_frobenius = r8mat_norm_fro ( n, n, a )
      write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) 
     &  title, n, norm_frobenius, error_ab, error_ac
c
c  CARRY matrix.
c
      title = 'CARRY'
      n = 5
      seed = 123456789
      k = i4_uniform ( 2, 20, seed )
      call carry ( k, n, a )
      call carry_inverse ( k, n, b )
      call r8mat_inverse ( n, a, c )
      call r8mat_is_inverse ( n, a, b, error_ab )
      call r8mat_is_inverse ( n, a, c, error_ac )
      norm_frobenius = r8mat_norm_fro ( n, n, a )
      write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) 
     &  title, n, norm_frobenius, error_ab, error_ac
c
c  CAUCHY matrix.
c
      title = 'CAUCHY'
      n = 5
      seed = 123456789
      call r8vec_uniform_01 ( n, seed, x )
      call r8vec_uniform_01 ( n, seed, y )
      call cauchy ( n, x, y, a )
      call cauchy_inverse ( n, x, y, b )
      call r8mat_inverse ( n, a, c )
      call r8mat_is_inverse ( n, a, b, error_ab )
      call r8mat_is_inverse ( n, a, c, error_ac )
      norm_frobenius = r8mat_norm_fro ( n, n, a )
      write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) 
     &  title, n, norm_frobenius, error_ab, error_ac
c
c  CHEBY_T matrix.
c
      title = 'CHEBY_T'
      n = 5
      call cheby_t ( n, a )
      call cheby_t_inverse ( n, b )
      call r8mat_inverse ( n, a, c )
      call r8mat_is_inverse ( n, a, b, error_ab )
      call r8mat_is_inverse ( n, a, c, error_ac )
      norm_frobenius = r8mat_norm_fro ( n, n, a )
      write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) 
     &  title, n, norm_frobenius, error_ab, error_ac
c
c  CHEBY_U matrix.
c
      title = 'CHEBY_U'
      n = 5
      call cheby_u ( n, a )
      call cheby_u_inverse ( n, b )
      call r8mat_inverse ( n, a, c )
      call r8mat_is_inverse ( n, a, b, error_ab )
      call r8mat_is_inverse ( n, a, c, error_ac )
      norm_frobenius = r8mat_norm_fro ( n, n, a )
      write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) 
     &  title, n, norm_frobenius, error_ab, error_ac
c
c  CHEBY_VAN2 matrix.
c
      title = 'CHEBY_VAN2'
      n = 5
      call cheby_van2 ( n, a )
      call cheby_van2_inverse ( n, b )
      call r8mat_inverse ( n, a, c )
      call r8mat_is_inverse ( n, a, b, error_ab )
      call r8mat_is_inverse ( n, a, c, error_ac )
      norm_frobenius = r8mat_norm_fro ( n, n, a )
      write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) 
     &  title, n, norm_frobenius, error_ab, error_ac
c
c  CHEBY_VAN3 matrix.
c
      title = 'CHEBY_VAN3'
      n = 5
      call cheby_van3 ( n, a )
      call cheby_van3_inverse ( n, b )
      call r8mat_inverse ( n, a, c )
      call r8mat_is_inverse ( n, a, b, error_ab )
      call r8mat_is_inverse ( n, a, c, error_ac )
      norm_frobenius = r8mat_norm_fro ( n, n, a )
      write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) 
     &  title, n, norm_frobenius, error_ab, error_ac
c
c  CHOW matrix.
c
      title = 'CHOW'
      n = 5
      seed = 123456789
      alpha = r8_uniform_01 ( seed )
      alpha = anint ( 50.0D+00 * alpha ) / 5.0D+00
      beta = r8_uniform_01 ( seed )
      beta = anint ( 50.0D+00 * beta ) / 5.0D+00
      call chow ( alpha, beta, n, n, a )
      call chow_inverse ( alpha, beta, n, b )
      call r8mat_inverse ( n, a, c )
      call r8mat_is_inverse ( n, a, b, error_ab )
      call r8mat_is_inverse ( n, a, c, error_ac )
      norm_frobenius = r8mat_norm_fro ( n, n, a )
      write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) 
     &  title, n, norm_frobenius, error_ab, error_ac
c
c  CIRCULANT matrix.
c
      if ( .false. ) then
      title = 'CIRCULANT'
      n = 5
      seed = 123456789
      call r8vec_uniform_01 ( n, seed, x )
      do i = 1, n
        x(i) = anint ( 50.0D+00 * x(i) - 25.0D+00 ) / 5.0D+00
      end do
      call circulant ( n, n, x, a )
c     call circulant_inverse ( n, x, b )
      call r8mat_inverse ( n, a, c )
      call r8mat_is_inverse ( n, a, b, error_ab )
      call r8mat_is_inverse ( n, a, c, error_ac )
      norm_frobenius = r8mat_norm_fro ( n, n, a )
      write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) 
     &  title, n, norm_frobenius, error_ab, error_ac
      else
      write ( *, '(2x,a)' ) 'CIRCULANT --- NOT READY'
      end if
c
c  CIRCULANT2 matrix.
c
      title = 'CIRCULANT2'
      n = 5
      call circulant2 ( n, a )
      call circulant2_inverse ( n, b )
      call r8mat_inverse ( n, a, c )
      call r8mat_is_inverse ( n, a, b, error_ab )
      call r8mat_is_inverse ( n, a, c, error_ac )
      norm_frobenius = r8mat_norm_fro ( n, n, a )
      write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) 
     &  title, n, norm_frobenius, error_ab, error_ac
c
c  CLEMENT1 matrix.
c
      title = 'CLEMENT1'
      n = 6
      call clement1 ( n, a )
      call clement1_inverse ( n, b )
      call r8mat_inverse ( n, a, c )
      call r8mat_is_inverse ( n, a, b, error_ab )
      call r8mat_is_inverse ( n, a, c, error_ac )
      norm_frobenius = r8mat_norm_fro ( n, n, a )
      write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) 
     &  title, n, norm_frobenius, error_ab, error_ac
c
c  CLEMENT2 matrix.
c
      title = 'CLEMENT2'
      n = 6
      call clement2 ( n, a )
      call clement2_inverse ( n, b )
      call r8mat_inverse ( n, a, c )
      call r8mat_is_inverse ( n, a, b, error_ab )
      call r8mat_is_inverse ( n, a, c, error_ac )
      norm_frobenius = r8mat_norm_fro ( n, n, a )
      write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) 
     &  title, n, norm_frobenius, error_ab, error_ac
c
c  CLEMENT3.
c
      title = 'CLEMENT3'
      n = 6
      seed = 123456789
      call r8vec_uniform_01 ( n - 1, seed, x )
      do i = 1, n - 1
        x(i) = anint ( 50.0D+00 * x(i) - 25.0D+00 ) / 5.0D+00
      end do
      call r8vec_uniform_01 ( n - 1, seed, y )
      do i = 1, n - 1
        y(i) = anint ( 50.0D+00 * y(i) - 25.0D+00 ) / 5.0D+00
      end do
      call clement3 ( n, x, y, a )
      call clement3_inverse ( n, x, y, b )
      call r8mat_inverse ( n, a, c )
      call r8mat_is_inverse ( n, a, b, error_ab )
      call r8mat_is_inverse ( n, a, c, error_ac )
      norm_frobenius = r8mat_norm_fro ( n, n, a )
      write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) 
     &  title, n, norm_frobenius, error_ab, error_ac
c
c  COMBIN matrix.
c
      title = 'COMBIN'
      n = 5
      seed = 123456789
      alpha = r8_uniform_01 ( seed )
      alpha = anint ( 50.0D+00 * alpha ) / 5.0D+00
      beta = r8_uniform_01 ( seed )
      beta = anint ( 50.0D+00 * beta ) / 5.0D+00
      call combin ( alpha, beta, n, a )
      call combin_inverse ( alpha, beta, n, b )
      call r8mat_inverse ( n, a, c )
      call r8mat_is_inverse ( n, a, b, error_ab )
      call r8mat_is_inverse ( n, a, c, error_ac )
      norm_frobenius = r8mat_norm_fro ( n, n, a )
      write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) 
     &  title, n, norm_frobenius, error_ab, error_ac
c
c  COMPANION.
c
      title = 'COMPANION'
      n = 5
      seed = 123456789
      call r8vec_uniform_01 ( n, seed, x )
      do i = 1, n
        x(i) = anint ( 10.0D+00 * x(i) - 5.0D+00 )
      end do
      call companion ( n, x, a )
      call companion_inverse ( n, x, b )
      call r8mat_inverse ( n, a, c )
      call r8mat_is_inverse ( n, a, b, error_ab )
      call r8mat_is_inverse ( n, a, c, error_ac )
      norm_frobenius = r8mat_norm_fro ( n, n, a )
      write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) 
     &  title, n, norm_frobenius, error_ab, error_ac
c
c  COMPLEX_I
c
      title = 'COMPLEX_I'
      n = 2
      call complex_i ( a )
      call complex_i_inverse ( b )
      call r8mat_inverse ( n, a, c )
      call r8mat_is_inverse ( n, a, b, error_ab )
      call r8mat_is_inverse ( n, a, c, error_ac )
      norm_frobenius = r8mat_norm_fro ( n, n, a )
      write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) 
     &  title, n, norm_frobenius, error_ab, error_ac
c
c  CONEX1 matrix.
c
      title = 'CONEX1'
      n = 4
      seed = 123456789
      alpha = r8_uniform_01 ( seed )
      alpha = anint ( 50.0D+00 * alpha ) / 5.0D+00
      call conex1 ( alpha, a )
      call conex1_inverse ( alpha, b )
      call r8mat_inverse ( n, a, c )
      call r8mat_is_inverse ( n, a, b, error_ab )
      call r8mat_is_inverse ( n, a, c, error_ac )
      norm_frobenius = r8mat_norm_fro ( n, n, a )
      write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) 
     &  title, n, norm_frobenius, error_ab, error_ac
c
c  CONEX2 matrix.
c
      title = 'CONEX2'
      n = 3
      seed = 123456789
      alpha = r8_uniform_01 ( seed )
      alpha = anint ( 50.0D+00 * alpha ) / 5.0D+00
      call conex2 ( alpha, a )
      call conex2_inverse ( alpha, b )
      call r8mat_inverse ( n, a, c )
      call r8mat_is_inverse ( n, a, b, error_ab )
      call r8mat_is_inverse ( n, a, c, error_ac )
      norm_frobenius = r8mat_norm_fro ( n, n, a )
      write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) 
     &  title, n, norm_frobenius, error_ab, error_ac
c
c  CONEX3 matrix.
c
      title = 'CONEX3'
      n = 5
      call conex3 ( n, a )
      call conex3_inverse ( n, b )
      call r8mat_inverse ( n, a, c )
      call r8mat_is_inverse ( n, a, b, error_ab )
      call r8mat_is_inverse ( n, a, c, error_ac )
      norm_frobenius = r8mat_norm_fro ( n, n, a )
      write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) 
     &  title, n, norm_frobenius, error_ab, error_ac
c
c  CONFERENCE matrix.
c
      title = 'CONFERENCE'
      n = 6
      call conference ( n, a )
      call conference_inverse ( n, b )
      call r8mat_inverse ( n, a, c )
      call r8mat_is_inverse ( n, a, b, error_ab )
      call r8mat_is_inverse ( n, a, c, error_ac )
      norm_frobenius = r8mat_norm_fro ( n, n, a )
      write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) 
     &  title, n, norm_frobenius, error_ab, error_ac
c
c  DAUB2 matrix.
c
      title = 'DAUB2'
      n = 4
      call daub2 ( n, a )
      call daub2_inverse ( n, b )
      call r8mat_inverse ( n, a, c )
      call r8mat_is_inverse ( n, a, b, error_ab )
      call r8mat_is_inverse ( n, a, c, error_ac )
      norm_frobenius = r8mat_norm_fro ( n, n, a )
      write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) 
     &  title, n, norm_frobenius, error_ab, error_ac
c
c  DAUB4 matrix.
c
      title = 'DAUB4'
      n = 8
      call daub4 ( n, a )
      call daub4_inverse ( n, b )
      call r8mat_inverse ( n, a, c )
      call r8mat_is_inverse ( n, a, b, error_ab )
      call r8mat_is_inverse ( n, a, c, error_ac )
      norm_frobenius = r8mat_norm_fro ( n, n, a )
      write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) 
     &  title, n, norm_frobenius, error_ab, error_ac
c
c  DAUB6.
c
      title = 'DAUB6'
      n = 12
      call daub6 ( n, a )
      call daub6_inverse ( n, b )
      call r8mat_inverse ( n, a, c )
      call r8mat_is_inverse ( n, a, b, error_ab )
      call r8mat_is_inverse ( n, a, c, error_ac )
      norm_frobenius = r8mat_norm_fro ( n, n, a )
      write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) 
     &  title, n, norm_frobenius, error_ab, error_ac
c
c  DAUB8.
c
      title = 'DAUB8'
      n = 16
      call daub8 ( n, a )
      call daub8_inverse ( n, b )
      call r8mat_inverse ( n, a, c )
      call r8mat_is_inverse ( n, a, b, error_ab )
      call r8mat_is_inverse ( n, a, c, error_ac )
      norm_frobenius = r8mat_norm_fro ( n, n, a )
      write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) 
     &  title, n, norm_frobenius, error_ab, error_ac
c
c  DAUB10 matrix.
c
      title = 'DAUB10'
      n = 20
      call daub10 ( n, a )
      call daub10_inverse ( n, b )
      call r8mat_inverse ( n, a, c )
      call r8mat_is_inverse ( n, a, b, error_ab )
      call r8mat_is_inverse ( n, a, c, error_ac )
      norm_frobenius = r8mat_norm_fro ( n, n, a )
      write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) 
     &  title, n, norm_frobenius, error_ab, error_ac
c
c  DAUB12 matrix.
c
      title = 'DAUB12'
      n = 24
      call daub12 ( n, a )
      call daub12_inverse ( n, b )
      call r8mat_inverse ( n, a, c )
      call r8mat_is_inverse ( n, a, b, error_ab )
      call r8mat_is_inverse ( n, a, c, error_ac )
      norm_frobenius = r8mat_norm_fro ( n, n, a )
      write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) 
     &  title, n, norm_frobenius, error_ab, error_ac
c
c  DIAGONAL.
c
      title = 'DIAGONAL'
      n = 5
      seed = 123456789
      call r8vec_uniform_01 ( n, seed, x )
      do i = 1, n
        x(i) = anint ( 50.0D+00 * x(i) - 25.0D+00 ) / 5.0D+00
      end do
      call diagonal ( n, n, x, a )
      call diagonal_inverse ( n, x, b )
      call r8mat_inverse ( n, a, c )
      call r8mat_is_inverse ( n, a, b, error_ab )
      call r8mat_is_inverse ( n, a, c, error_ac )
      norm_frobenius = r8mat_norm_fro ( n, n, a )
      write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) 
     &  title, n, norm_frobenius, error_ab, error_ac
c
c  DIF2 matrix.
c
      title = 'DIF2'
      n = 5
      call dif2 ( n, n, a )
      call dif2_inverse ( n, b )
      call r8mat_inverse ( n, a, c )
      call r8mat_is_inverse ( n, a, b, error_ab )
      call r8mat_is_inverse ( n, a, c, error_ac )
      norm_frobenius = r8mat_norm_fro ( n, n, a )
      write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) 
     &  title, n, norm_frobenius, error_ab, error_ac
c
c  DOWNSHIFT matrix.
c
      title = 'DOWNSHIFT'
      n = 5
      call downshift ( n, a )
      call downshift_inverse ( n, b )
      call r8mat_inverse ( n, a, c )
      call r8mat_is_inverse ( n, a, b, error_ab )
      call r8mat_is_inverse ( n, a, c, error_ac )
      norm_frobenius = r8mat_norm_fro ( n, n, a )
      write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) 
     &  title, n, norm_frobenius, error_ab, error_ac
c
c  DRMAC
c

c
c  EULERIAN matrix.
c
      title = 'EULERIAN'
      n = 5
      call eulerian ( n, n, a )
      call eulerian_inverse ( n, b )
      call r8mat_inverse ( n, a, c )
      call r8mat_is_inverse ( n, a, b, error_ab )
      call r8mat_is_inverse ( n, a, c, error_ac )
      norm_frobenius = r8mat_norm_fro ( n, n, a )
      write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) 
     &  title, n, norm_frobenius, error_ab, error_ac
c
c  EXCHANGE matrix.
c
      title = 'EXCHANGE'
      n = 5
      call exchange ( n, n, a )
      call exchange_inverse ( n, b )
      call r8mat_inverse ( n, a, c )
      call r8mat_is_inverse ( n, a, b, error_ab )
      call r8mat_is_inverse ( n, a, c, error_ac )
      norm_frobenius = r8mat_norm_fro ( n, n, a )
      write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) 
     &  title, n, norm_frobenius, error_ab, error_ac
c
c  FIBONACCI2 matrix.
c
      title = 'FIBONACCI2'
      n = 5
      call fibonacci2 ( n, a )
      call fibonacci2_inverse ( n, b )
      call r8mat_inverse ( n, a, c )
      call r8mat_is_inverse ( n, a, b, error_ab )
      call r8mat_is_inverse ( n, a, c, error_ac )
      norm_frobenius = r8mat_norm_fro ( n, n, a )
      write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) 
     &  title, n, norm_frobenius, error_ab, error_ac
c
c  FIBONACCI3 matrix.
c
      title = 'FIBONACCI3'
      n = 5
      call fibonacci3 ( n, a )
      call fibonacci3_inverse ( n, b )
      call r8mat_inverse ( n, a, c )
      call r8mat_is_inverse ( n, a, b, error_ab )
      call r8mat_is_inverse ( n, a, c, error_ac )
      norm_frobenius = r8mat_norm_fro ( n, n, a )
      write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) 
     &  title, n, norm_frobenius, error_ab, error_ac
c
c  FIEDLER.
c  The FIEDLER_INVERSE routine assumes the X vector is sorted.
c
      title = 'FIEDLER'
      n = 7
      seed = 123456789
      call r8vec_uniform_01 ( n, seed, x )
      do i = 1, n
        x(i) = anint ( 50.0D+00 * x(i) - 25.0D+00 ) / 5.0D+00
      end do
      call r8vec_sort_bubble_a ( n, x )
      call fiedler ( n, n, x, a )
      call fiedler_inverse ( n, x, b )
      call r8mat_inverse ( n, a, c )
      call r8mat_is_inverse ( n, a, b, error_ab )
      call r8mat_is_inverse ( n, a, c, error_ac )
      norm_frobenius = r8mat_norm_fro ( n, n, a )
      write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) 
     &  title, n, norm_frobenius, error_ab, error_ac
c
c  FORSYTHE matrix.
c
      title = 'FORSYTHE'
      n = 5
      seed = 123456789
      alpha = r8_uniform_01 ( seed )
      alpha = anint ( 50.0D+00 * alpha ) / 5.0D+00
      beta = r8_uniform_01 ( seed )
      beta = anint ( 50.0D+00 * beta ) / 5.0D+00
      call forsythe ( alpha, beta, n, a )
      call forsythe_inverse ( alpha, beta, n, b )
      call r8mat_inverse ( n, a, c )
      call r8mat_is_inverse ( n, a, b, error_ab )
      call r8mat_is_inverse ( n, a, c, error_ac )
      norm_frobenius = r8mat_norm_fro ( n, n, a )
      write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) 
     &  title, n, norm_frobenius, error_ab, error_ac
c
c  FOURIER_COSINE matrix.
c
      title = 'FOURIER_COSINE'
      n = 5
      call fourier_cosine ( n, a )
      call fourier_cosine_inverse ( n, b )
      call r8mat_inverse ( n, a, c )
      call r8mat_is_inverse ( n, a, b, error_ab )
      call r8mat_is_inverse ( n, a, c, error_ac )
      norm_frobenius = r8mat_norm_fro ( n, n, a )
      write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) 
     &  title, n, norm_frobenius, error_ab, error_ac
c
c  FOURIER_SINE matrix.
c
      title = 'FOURIER_SINE'
      n = 5
      call fourier_sine ( n, a )
      call fourier_sine_inverse ( n, b )
      call r8mat_inverse ( n, a, c )
      call r8mat_is_inverse ( n, a, b, error_ab )
      call r8mat_is_inverse ( n, a, c, error_ac )
      norm_frobenius = r8mat_norm_fro ( n, n, a )
      write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) 
     &  title, n, norm_frobenius, error_ab, error_ac
c
c  FRANK matrix.
c
      title = 'FRANK'
      n = 5
      call frank ( n, a )
      call frank_inverse ( n, b )
      call r8mat_inverse ( n, a, c )
      call r8mat_is_inverse ( n, a, b, error_ab )
      call r8mat_is_inverse ( n, a, c, error_ac )
      norm_frobenius = r8mat_norm_fro ( n, n, a )
      write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) 
     &  title, n, norm_frobenius, error_ab, error_ac
c
c  GFPP matrix.
c
      title = 'GFPP'
      n = 5
      seed = 123456789
      alpha = r8_uniform_01 ( seed )
      call gfpp ( n, alpha, a )
      call gfpp_inverse ( n, alpha, b )
      call r8mat_inverse ( n, a, c )
      call r8mat_is_inverse ( n, a, b, error_ab )
      call r8mat_is_inverse ( n, a, c, error_ac )
      norm_frobenius = r8mat_norm_fro ( n, n, a )
      write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) 
     &  title, n, norm_frobenius, error_ab, error_ac
c
c  GIVENS matrix.
c
      title = 'GIVENS'
      n = 5
      call givens ( n, n, a )
      call givens_inverse ( n, b )
      call r8mat_inverse ( n, a, c )
      call r8mat_is_inverse ( n, a, b, error_ab )
      call r8mat_is_inverse ( n, a, c, error_ac )
      norm_frobenius = r8mat_norm_fro ( n, n, a )
      write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) 
     &  title, n, norm_frobenius, error_ab, error_ac
c
c  GK316 matrix.
c
      title = 'GK316'
      n = 5
      call gk316 ( n, a )
      call gk316_inverse ( n, b )
      call r8mat_inverse ( n, a, c )
      call r8mat_is_inverse ( n, a, b, error_ab )
      call r8mat_is_inverse ( n, a, c, error_ac )
      norm_frobenius = r8mat_norm_fro ( n, n, a )
      write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) 
     &  title, n, norm_frobenius, error_ab, error_ac
c
c  GK323 matrix.
c
      title = 'GK323'
      n = 5
      call gk323 ( n, n, a )
      call gk323_inverse ( n, b )
      call r8mat_inverse ( n, a, c )
      call r8mat_is_inverse ( n, a, b, error_ab )
      call r8mat_is_inverse ( n, a, c, error_ac )
      norm_frobenius = r8mat_norm_fro ( n, n, a )
      write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) 
     &  title, n, norm_frobenius, error_ab, error_ac
c
c  GK324 matrix.
c
      title = 'GK324'
      n = 5
      seed = 123456789
      call r8vec_uniform_01 ( n, seed, x )
      do i = 1, n
        x(i) = anint ( 50.0D+00 * x(i) - 25.0D+00 ) / 5.0D+00
      end do
      call gk324 ( n, n, x, a )
      call gk324_inverse ( n, x, b )
      call r8mat_inverse ( n, a, c )
      call r8mat_is_inverse ( n, a, b, error_ab )
      call r8mat_is_inverse ( n, a, c, error_ac )
      norm_frobenius = r8mat_norm_fro ( n, n, a )
      write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) 
     &  title, n, norm_frobenius, error_ab, error_ac
c
c  HANOWA matrix.
c
      title = 'HANOWA'
      n = 8
      seed = 123456789
      alpha = r8_uniform_01 ( seed )
      call hanowa ( alpha, n, a )
      call hanowa_inverse ( alpha, n, b )
      call r8mat_inverse ( n, a, c )
      call r8mat_is_inverse ( n, a, b, error_ab )
      call r8mat_is_inverse ( n, a, c, error_ac )
      norm_frobenius = r8mat_norm_fro ( n, n, a )
      write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) 
     &  title, n, norm_frobenius, error_ab, error_ac
c
c  HARMAN matrix.
c
      title = 'HARMAN'
      n = 8
      call harman ( a )
      call harman_inverse ( b )
      call r8mat_inverse ( n, a, c )
      call r8mat_is_inverse ( n, a, b, error_ab )
      call r8mat_is_inverse ( n, a, c, error_ac )
      norm_frobenius = r8mat_norm_fro ( n, n, a )
      write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) 
     &  title, n, norm_frobenius, error_ab, error_ac
c
c  HARTLEY matrix.
c
      title = 'HARTLEY'
      n = 5
      call hartley ( n, a )
      call hartley_inverse ( n, b )
      call r8mat_inverse ( n, a, c )
      call r8mat_is_inverse ( n, a, b, error_ab )
      call r8mat_is_inverse ( n, a, c, error_ac )
      norm_frobenius = r8mat_norm_fro ( n, n, a )
      write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) 
     &  title, n, norm_frobenius, error_ab, error_ac
c
c  HELMERT matrix.
c
      title = 'HELMERT'
      n = 5
      call helmert ( n, a )
      call helmert_inverse ( n, b )
      call r8mat_inverse ( n, a, c )
      call r8mat_is_inverse ( n, a, b, error_ab )
      call r8mat_is_inverse ( n, a, c, error_ac )
      norm_frobenius = r8mat_norm_fro ( n, n, a )
      write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) 
     &  title, n, norm_frobenius, error_ab, error_ac
c
c  HELMERT2 matrix.
c
      title = 'HELMERT2'
      n = 5
      seed = 123456789
      call r8vec_uniform_01 ( n, seed, x )
      do i = 1, n
        x(i) = anint ( 50.0D+00 * x(i) - 25.0D+00 ) / 5.0D+00
      end do
      call helmert2 ( n, x, a )
      call helmert2_inverse ( n, x, b )
      call r8mat_inverse ( n, a, c )
      call r8mat_is_inverse ( n, a, b, error_ab )
      call r8mat_is_inverse ( n, a, c, error_ac )
      norm_frobenius = r8mat_norm_fro ( n, n, a )
      write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) 
     &  title, n, norm_frobenius, error_ab, error_ac
c
c  HERMITE matrix.
c
      title = 'HERMITE'
      n = 5
      call hermite ( n, a )
      call hermite_inverse ( n, b )
      call r8mat_inverse ( n, a, c )
      call r8mat_is_inverse ( n, a, b, error_ab )
      call r8mat_is_inverse ( n, a, c, error_ac )
      norm_frobenius = r8mat_norm_fro ( n, n, a )
      write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) 
     &  title, n, norm_frobenius, error_ab, error_ac
c
c  HERNDON matrix.
c
      title = 'HERNDON'
      n = 5
      call herndon ( n, a )
      call herndon_inverse ( n, b )
      call r8mat_inverse ( n, a, c )
      call r8mat_is_inverse ( n, a, b, error_ab )
      call r8mat_is_inverse ( n, a, c, error_ac )
      norm_frobenius = r8mat_norm_fro ( n, n, a )
      write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) 
     &  title, n, norm_frobenius, error_ab, error_ac
c
c  HILBERT matrix.
c
      title = 'HILBERT'
      n = 5
      call hilbert ( n, n, a )
      call hilbert_inverse ( n, b )
      call r8mat_inverse ( n, a, c )
      call r8mat_is_inverse ( n, a, b, error_ab )
      call r8mat_is_inverse ( n, a, c, error_ac )
      norm_frobenius = r8mat_norm_fro ( n, n, a )
      write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) 
     &  title, n, norm_frobenius, error_ab, error_ac
c
c  HOUSEHOLDER matrix.
c
      title = 'HOUSEHOLDER'
      n = 5
      seed = 123456789
      call r8vec_uniform_01 ( n, seed, x )
      do i = 1, n
        x(i) = anint ( 50.0D+00 * x(i) - 25.0D+00 ) / 5.0D+00
      end do
      call householder ( n, x, a )
      call householder_inverse ( n, x, b )
      call r8mat_inverse ( n, a, c )
      call r8mat_is_inverse ( n, a, b, error_ab )
      call r8mat_is_inverse ( n, a, c, error_ac )
      norm_frobenius = r8mat_norm_fro ( n, n, a )
      write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) 
     &  title, n, norm_frobenius, error_ab, error_ac
c
c  IDENTITY matrix.
c
      title = 'IDENTITY'
      n = 5
      call identity ( n, n, a )
      call identity_inverse ( n, b )
      call r8mat_inverse ( n, a, c )
      call r8mat_is_inverse ( n, a, b, error_ab )
      call r8mat_is_inverse ( n, a, c, error_ac )
      norm_frobenius = r8mat_norm_fro ( n, n, a )
      write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) 
     &  title, n, norm_frobenius, error_ab, error_ac
c
c  ILL3 matrix.
c
      title = 'ILL3'
      n = 3
      call ill3 ( a )
      call ill3_inverse ( b )
      call r8mat_inverse ( n, a, c )
      call r8mat_is_inverse ( n, a, b, error_ab )
      call r8mat_is_inverse ( n, a, c, error_ac )
      norm_frobenius = r8mat_norm_fro ( n, n, a )
      write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) 
     &  title, n, norm_frobenius, error_ab, error_ac
c
c  INTEGRATION matrix.
c
      title = 'INTEGRATION'
      n = 6
      seed = 123456789
      alpha = r8_uniform_01 ( seed )
      alpha = anint ( 50.0D+00 * alpha ) / 5.0D+00
      call integration ( alpha, n, a )
      call integration_inverse ( alpha, n, b )
      call r8mat_inverse ( n, a, c )
      call r8mat_is_inverse ( n, a, b, error_ab )
      call r8mat_is_inverse ( n, a, c, error_ac )
      norm_frobenius = r8mat_norm_fro ( n, n, a )
      write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) 
     &  title, n, norm_frobenius, error_ab, error_ac
c
c  INVOL matrix.
c
      title = 'INVOL'
      n = 5
      call invol ( n, a )
      call invol_inverse ( n, b )
      call r8mat_inverse ( n, a, c )
      call r8mat_is_inverse ( n, a, b, error_ab )
      call r8mat_is_inverse ( n, a, c, error_ac )
      norm_frobenius = r8mat_norm_fro ( n, n, a )
      write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) 
     &  title, n, norm_frobenius, error_ab, error_ac
c
c  JORDAN matrix.
c
      title = 'JORDAN'
      n = 6
      seed = 123456789
      alpha = r8_uniform_01 ( seed )
      alpha = anint ( 50.0D+00 * alpha ) / 5.0D+00
      call jordan ( alpha, n, n, a )
      call jordan_inverse ( alpha, n, b )
      call r8mat_inverse ( n, a, c )
      call r8mat_is_inverse ( n, a, b, error_ab )
      call r8mat_is_inverse ( n, a, c, error_ac )
      norm_frobenius = r8mat_norm_fro ( n, n, a )
      write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) 
     &  title, n, norm_frobenius, error_ab, error_ac
c
c  KAHAN matrix.
c
      title = 'KAHAN'
      n = 5
      seed = 123456789
      alpha = r8_uniform_01 ( seed )
      alpha = anint ( 50.0D+00 * alpha ) / 5.0D+00
      call kahan ( alpha, n, n, a )
      call kahan_inverse ( alpha, n, b )
      call r8mat_inverse ( n, a, c )
      call r8mat_is_inverse ( n, a, b, error_ab )
      call r8mat_is_inverse ( n, a, c, error_ac )
      norm_frobenius = r8mat_norm_fro ( n, n, a )
      write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) 
     &  title, n, norm_frobenius, error_ab, error_ac
c
c  KERSHAW matrix.
c
      title = 'KERSHAW'
      n = 4
      call kershaw ( a )
      call kershaw_inverse ( b )
      call r8mat_inverse ( n, a, c )
      call r8mat_is_inverse ( n, a, b, error_ab )
      call r8mat_is_inverse ( n, a, c, error_ac )
      norm_frobenius = r8mat_norm_fro ( n, n, a )
      write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) 
     &  title, n, norm_frobenius, error_ab, error_ac
c
c  KERSHAWTRI matrix.
c
      title = 'KERSHAWTRI'
      n = 5
      x_n = ( n + 1 ) / 2
      seed = 123456789
      call r8vec_uniform_01 ( x_n, seed, x )
      do i = 1, x_n
        x(i) = anint ( 50.0D+00 * x(i) - 25.0D+00 ) / 5.0D+00
      end do
      call kershawtri ( n, x, a )
      call kershawtri_inverse ( n, x, b )
      call r8mat_inverse ( n, a, c )
      call r8mat_is_inverse ( n, a, b, error_ab )
      call r8mat_is_inverse ( n, a, c, error_ac )
      norm_frobenius = r8mat_norm_fro ( n, n, a )
      write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) 
     &  title, n, norm_frobenius, error_ab, error_ac
c
c  KMS matrix.
c
      title = 'KMS'
      n = 5
      seed = 123456789
      alpha = r8_uniform_01 ( seed )
      alpha = anint ( 50.0D+00 * alpha ) / 5.0D+00
      call kms ( alpha, n, n, a )
      call kms_inverse ( alpha, n, b )
      call r8mat_inverse ( n, a, c )
      call r8mat_is_inverse ( n, a, b, error_ab )
      call r8mat_is_inverse ( n, a, c, error_ac )
      norm_frobenius = r8mat_norm_fro ( n, n, a )
      write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) 
     &  title, n, norm_frobenius, error_ab, error_ac
c
c  LAGUERRE matrix.
c
      title = 'LAGUERRE'
      n = 5
      call laguerre ( n, a )
      call laguerre_inverse ( n, b )
      call r8mat_inverse ( n, a, c )
      call r8mat_is_inverse ( n, a, b, error_ab )
      call r8mat_is_inverse ( n, a, c, error_ac )
      norm_frobenius = r8mat_norm_fro ( n, n, a )
      write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) 
     &  title, n, norm_frobenius, error_ab, error_ac
c
c  LEGENDRE matrix.
c
      title = 'LEGENDRE'
      n = 5
      call legendre ( n, a )
      call legendre_inverse ( n, b )
      call r8mat_inverse ( n, a, c )
      call r8mat_is_inverse ( n, a, b, error_ab )
      call r8mat_is_inverse ( n, a, c, error_ac )
      norm_frobenius = r8mat_norm_fro ( n, n, a )
      write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) 
     &  title, n, norm_frobenius, error_ab, error_ac
c
c  LEHMER matrix.
c
      title = 'LEHMER'
      n = 5
      call lehmer ( n, n, a )
      call lehmer_inverse ( n, b )
      call r8mat_inverse ( n, a, c )
      call r8mat_is_inverse ( n, a, b, error_ab )
      call r8mat_is_inverse ( n, a, c, error_ac )
      norm_frobenius = r8mat_norm_fro ( n, n, a )
      write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) 
     &  title, n, norm_frobenius, error_ab, error_ac
c
c  LIETZKE matrix.
c
      title = 'LIETZKE'
      n = 5
      call lietzke ( n, a )
      call lietzke_inverse ( n, b )
      call r8mat_inverse ( n, a, c )
      call r8mat_is_inverse ( n, a, b, error_ab )
      call r8mat_is_inverse ( n, a, c, error_ac )
      norm_frobenius = r8mat_norm_fro ( n, n, a )
      write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) 
     &  title, n, norm_frobenius, error_ab, error_ac
c
c  LOTKIN matrix.
c
      title = 'LOTKIN'
      n = 5
      call lotkin ( n, n, a )
      call lotkin_inverse ( n, b )
      call r8mat_inverse ( n, a, c )
      call r8mat_is_inverse ( n, a, b, error_ab )
      call r8mat_is_inverse ( n, a, c, error_ac )
      norm_frobenius = r8mat_norm_fro ( n, n, a )
      write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) 
     &  title, n, norm_frobenius, error_ab, error_ac
c
c  MAXIJ matrix.
c
      title = 'MAXIJ'
      n = 5
      call maxij ( n, n, a )
      call maxij_inverse ( n, b )
      call r8mat_inverse ( n, a, c )
      call r8mat_is_inverse ( n, a, b, error_ab )
      call r8mat_is_inverse ( n, a, c, error_ac )
      norm_frobenius = r8mat_norm_fro ( n, n, a )
      write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) 
     &  title, n, norm_frobenius, error_ab, error_ac
c
c  MILNES matrix.
c
      title = 'MILNES'
      n = 5
      seed = 123456789
      call r8vec_uniform_01 ( n, seed, x )
      do i = 1, n
        x(i) = anint ( 50.0D+00 * x(i) - 25.0D+00 ) / 5.0D+00
      end do
      call milnes ( n, n, x, a )
      call milnes_inverse ( n, x, b )
      call r8mat_inverse ( n, a, c )
      call r8mat_is_inverse ( n, a, b, error_ab )
      call r8mat_is_inverse ( n, a, c, error_ac )
      norm_frobenius = r8mat_norm_fro ( n, n, a )
      write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) 
     &  title, n, norm_frobenius, error_ab, error_ac
c
c  MINIJ matrix.
c
      title = 'MINIJ'
      n = 5
      call minij ( n, n, a )
      call minij_inverse ( n, b )
      call r8mat_inverse ( n, a, c )
      call r8mat_is_inverse ( n, a, b, error_ab )
      call r8mat_is_inverse ( n, a, c, error_ac )
      norm_frobenius = r8mat_norm_fro ( n, n, a )
      write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) 
     &  title, n, norm_frobenius, error_ab, error_ac
c
c  MOLER1 matrix.
c
      title = 'MOLER1'
      n = 5
      seed = 123456789
      alpha = r8_uniform_01 ( seed )
      alpha = anint ( 50.0D+00 * alpha ) / 5.0D+00
      call moler1 ( alpha, n, n, a )
      call moler1_inverse ( alpha, n, b )
      call r8mat_inverse ( n, a, c )
      call r8mat_is_inverse ( n, a, b, error_ab )
      call r8mat_is_inverse ( n, a, c, error_ac )
      norm_frobenius = r8mat_norm_fro ( n, n, a )
      write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) 
     &  title, n, norm_frobenius, error_ab, error_ac
c
c  MOLER3 matrix.
c
      title = 'MOLER3'
      n = 5
      call moler3 ( n, n, a )
      call moler3_inverse ( n, b )
      call r8mat_inverse ( n, a, c )
      call r8mat_is_inverse ( n, a, b, error_ab )
      call r8mat_is_inverse ( n, a, c, error_ac )
      norm_frobenius = r8mat_norm_fro ( n, n, a )
      write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) 
     &  title, n, norm_frobenius, error_ab, error_ac
c
c  ORTEGA matrix.
c
      title = 'ORTEGA'
      n = 5
      seed = 123456789
      call r8vec_uniform_01 ( n, seed, v1 )
      do i = 1, n
        v1(i) = anint ( 50.0D+00 * v1(i) - 25.0D+00  ) / 5.0D+00
      end do
      call r8vec_uniform_01 ( n, seed, v2 )
      do i = 1, n
        v2(i) = anint ( 50.0D+00 * v2(i) - 25.0D+00  ) / 5.0D+00
      end do
      call r8vec_uniform_01 ( n, seed, v3 )
      do i = 1, n
        v3(i) = anint ( 50.0D+00 * v3(i) - 25.0D+00  ) / 5.0D+00
      end do
      call ortega ( n, v1, v2, v3, a )
      call ortega_inverse ( n, v1, v2, v3, b )
      call r8mat_inverse ( n, a, c )
      call r8mat_is_inverse ( n, a, b, error_ab )
      call r8mat_is_inverse ( n, a, c, error_ac )
      norm_frobenius = r8mat_norm_fro ( n, n, a )
      write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) 
     &  title, n, norm_frobenius, error_ab, error_ac
c
c  ORTH_SYMM matrix.
c
      title = 'ORTH_SYMM'
      n = 5
      call orth_symm ( n, a )
      call orth_symm_inverse ( n, b )
      call r8mat_inverse ( n, a, c )
      call r8mat_is_inverse ( n, a, b, error_ab )
      call r8mat_is_inverse ( n, a, c, error_ac )
      norm_frobenius = r8mat_norm_fro ( n, n, a )
      write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) 
     &  title, n, norm_frobenius, error_ab, error_ac
c
c  OTO matrix.
c
      title = 'OTO'
      n = 5
      call oto ( n, n, a )
      call oto_inverse ( n, b )
      call r8mat_inverse ( n, a, c )
      call r8mat_is_inverse ( n, a, b, error_ab )
      call r8mat_is_inverse ( n, a, c, error_ac )
      norm_frobenius = r8mat_norm_fro ( n, n, a )
      write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) 
     &  title, n, norm_frobenius, error_ab, error_ac
c
c  PARTER matrix.
c
      title = 'PARTER'
      n = 5
      call parter ( n, n, a )
      call parter_inverse ( n, b )
      call r8mat_inverse ( n, a, c )
      call r8mat_is_inverse ( n, a, b, error_ab )
      call r8mat_is_inverse ( n, a, c, error_ac )
      norm_frobenius = r8mat_norm_fro ( n, n, a )
      write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) 
     &  title, n, norm_frobenius, error_ab, error_ac
c
c  PASCAL1 matrix.
c
      title = 'PASCAL1'
      n = 5
      call pascal1 ( n, a )
      call pascal1_inverse ( n, b )
      call r8mat_inverse ( n, a, c )
      call r8mat_is_inverse ( n, a, b, error_ab )
      call r8mat_is_inverse ( n, a, c, error_ac )
      norm_frobenius = r8mat_norm_fro ( n, n, a )
      write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) 
     &  title, n, norm_frobenius, error_ab, error_ac
c
c  PASCAL2 matrix.
c
      title = 'PASCAL2'
      n = 5
      call pascal2 ( n, a )
      call pascal2_inverse ( n, b )
      call r8mat_inverse ( n, a, c )
      call r8mat_is_inverse ( n, a, b, error_ab )
      call r8mat_is_inverse ( n, a, c, error_ac )
      norm_frobenius = r8mat_norm_fro ( n, n, a )
      write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) 
     &  title, n, norm_frobenius, error_ab, error_ac
c
c  PASCAL3 matrix.
c
      title = 'PASCAL3'
      n = 5
      seed = 123456789
      alpha = r8_uniform_01 ( seed )
      alpha = anint ( 50.0D+00 * alpha - 25.0D+00 ) / 5.0D+00
      call pascal3 ( n, alpha, a )
      call pascal3_inverse ( n, alpha, b )
      call r8mat_inverse ( n, a, c )
      call r8mat_is_inverse ( n, a, b, error_ab )
      call r8mat_is_inverse ( n, a, c, error_ac )
      norm_frobenius = r8mat_norm_fro ( n, n, a )
      write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) 
     &  title, n, norm_frobenius, error_ab, error_ac
c
c  PDS_RANDOM matrix.
c
      title = 'PDS_RANDOM'
      n = 5
      seed_save = 123456789
      seed = seed_save
      call pds_random ( n, seed, a )
      seed = seed_save
      call pds_random_inverse ( n, seed, b )
      call r8mat_inverse ( n, a, c )
      call r8mat_is_inverse ( n, a, b, error_ab )
      call r8mat_is_inverse ( n, a, c, error_ac )
      norm_frobenius = r8mat_norm_fro ( n, n, a )
      write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) 
     &  title, n, norm_frobenius, error_ab, error_ac
c
c  PEI matrix.
c
      title = 'PEI'
      n = 5
      seed = 123456789
      alpha = r8_uniform_01 ( seed )
      alpha = anint ( 50.0D+00 * alpha ) / 5.0D+00
      call pei ( alpha, n, a )
      call pei_inverse ( alpha, n, b )
      call r8mat_inverse ( n, a, c )
      call r8mat_is_inverse ( n, a, b, error_ab )
      call r8mat_is_inverse ( n, a, c, error_ac )
      norm_frobenius = r8mat_norm_fro ( n, n, a )
      write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) 
     &  title, n, norm_frobenius, error_ab, error_ac
c
c  PERMUTATION_RANDOM matrix.
c
      title = 'PERMUTATION_RANDOM'
      n = 5
      seed = 123456789
      seed_save = seed
      call permutation_random ( n, seed, a )
      seed = seed_save
      call permutation_random_inverse ( n, seed, b )
      call r8mat_inverse ( n, a, c )
      call r8mat_is_inverse ( n, a, b, error_ab )
      call r8mat_is_inverse ( n, a, c, error_ac )
      norm_frobenius = r8mat_norm_fro ( n, n, a )
      write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) 
     &  title, n, norm_frobenius, error_ab, error_ac
c
c  PLU matrix.
c
      title = 'PLU'
      n = 5
      do i = 1, n
        pivot(i) = i
      end do
      call plu ( n, pivot, p, l, u, a )
      call plu_inverse ( n, p, l, u, b )
      call r8mat_inverse ( n, a, c )
      call r8mat_is_inverse ( n, a, b, error_ab )
      call r8mat_is_inverse ( n, a, c, error_ac )
      norm_frobenius = r8mat_norm_fro ( n, n, a )
      write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) 
     &  title, n, norm_frobenius, error_ab, error_ac
c
c  RIS matrix.
c
      title = 'RIS'
      n = 5
      call ris ( n, a )
      call ris_inverse ( n, b )
      call r8mat_inverse ( n, a, c )
      call r8mat_is_inverse ( n, a, b, error_ab )
      call r8mat_is_inverse ( n, a, c, error_ac )
      norm_frobenius = r8mat_norm_fro ( n, n, a )
      write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) 
     &  title, n, norm_frobenius, error_ab, error_ac
c
c  RODMAN matrix.
c
      title = 'RODMAN'
      n = 5
      seed = 123456789
      alpha = r8_uniform_01 ( seed )
      alpha = anint ( 50.0D+00 * alpha - 25.0D+00 ) / 5.0D+00
      call rodman ( alpha, n, n, a )
      call rodman_inverse ( alpha, n, b )
      call r8mat_inverse ( n, a, c )
      call r8mat_is_inverse ( n, a, b, error_ab )
      call r8mat_is_inverse ( n, a, c, error_ac )
      norm_frobenius = r8mat_norm_fro ( n, n, a )
      write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) 
     &  title, n, norm_frobenius, error_ab, error_ac
c
c  RUTIS1 matrix.
c
      title = 'RUTIS1'
      n = 4
      call rutis1 ( a )
      call rutis1_inverse ( b )
      call r8mat_inverse ( n, a, c )
      call r8mat_is_inverse ( n, a, b, error_ab )
      call r8mat_is_inverse ( n, a, c, error_ac )
      norm_frobenius = r8mat_norm_fro ( n, n, a )
      write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) 
     &  title, n, norm_frobenius, error_ab, error_ac
c
c  RUTIS2 matrix.
c
      title = 'RUTIS2'
      n = 4
      call rutis2 ( a )
      call rutis2_inverse ( b )
      call r8mat_inverse ( n, a, c )
      call r8mat_is_inverse ( n, a, b, error_ab )
      call r8mat_is_inverse ( n, a, c, error_ac )
      norm_frobenius = r8mat_norm_fro ( n, n, a )
      write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) 
     &  title, n, norm_frobenius, error_ab, error_ac
c
c  RUTIS3 matrix.
c
      title = 'RUTIS3'
      n = 4
      call rutis3 ( a )
      call rutis3_inverse ( b )
      call r8mat_inverse ( n, a, c )
      call r8mat_is_inverse ( n, a, b, error_ab )
      call r8mat_is_inverse ( n, a, c, error_ac )
      norm_frobenius = r8mat_norm_fro ( n, n, a )
      write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) 
     &  title, n, norm_frobenius, error_ab, error_ac
c
c  RUTIS4 matrix.
c
      title = 'RUTIS4'
      n = 5
      call rutis4 ( n, a )
      call rutis4_inverse ( n, b )
      call r8mat_inverse ( n, a, c )
      call r8mat_is_inverse ( n, a, b, error_ab )
      call r8mat_is_inverse ( n, a, c, error_ac )
      norm_frobenius = r8mat_norm_fro ( n, n, a )
      write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) 
     &  title, n, norm_frobenius, error_ab, error_ac
c
c  RUTIS5 matrix.
c
      title = 'RUTIS5'
      n = 4
      call rutis5 ( a )
      call rutis5_inverse ( b )
      call r8mat_inverse ( n, a, c )
      call r8mat_is_inverse ( n, a, b, error_ab )
      call r8mat_is_inverse ( n, a, c, error_ac )
      norm_frobenius = r8mat_norm_fro ( n, n, a )
      write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) 
     &  title, n, norm_frobenius, error_ab, error_ac
c
c  SCHUR_BLOCK matrix.
c
      title = 'SCHUR_BLOCK'
      n = 5
      x_n = ( n + 1 ) / 2
      y_n = n / 2
      seed = 123456789
      call r8vec_uniform_01 ( x_n, seed, x )
      do i = 1, x_n
        x(i) = anint ( 50.0D+00 * x(i) - 25.0D+00 ) / 5.0D+00
      end do
      call r8vec_uniform_01 ( y_n, seed, y )
      do i = 1, y_n
        y(i) = anint ( 50.0D+00 * y(i) - 25.0D+00 ) / 5.0D+00
      end do
      call schur_block ( n, x, y, a )
      call schur_block_inverse ( n, x, y, b )
      call r8mat_inverse ( n, a, c )
      call r8mat_is_inverse ( n, a, b, error_ab )
      call r8mat_is_inverse ( n, a, c, error_ac )
      norm_frobenius = r8mat_norm_fro ( n, n, a )
      write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) 
     &  title, n, norm_frobenius, error_ab, error_ac
c
c  SPLINE matrix.
c
      title = 'SPLINE'
      n = 5
      seed = 123456789
      call r8vec_uniform_01 ( n - 1, seed, x )
      do i = 1, n - 1
        x(i) = anint ( 50.0D+00 * x(i) - 25.0D+00 ) / 5.0D+00
      end do
      call spline ( n, x, a )
      call spline_inverse ( n, x, b )
      call r8mat_inverse ( n, a, c )
      call r8mat_is_inverse ( n, a, b, error_ab )
      call r8mat_is_inverse ( n, a, c, error_ac )
      norm_frobenius = r8mat_norm_fro ( n, n, a )
      write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) 
     &  title, n, norm_frobenius, error_ab, error_ac
c
c  STIRLING matrix.
c
      title = 'STIRLING'
      n = 5
      call stirling ( n, n, a )
      call stirling_inverse ( n, b )
      call r8mat_inverse ( n, a, c )
      call r8mat_is_inverse ( n, a, b, error_ab )
      call r8mat_is_inverse ( n, a, c, error_ac )
      norm_frobenius = r8mat_norm_fro ( n, n, a )
      write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) 
     &  title, n, norm_frobenius, error_ab, error_ac
c
c  SUMMATION matrix.
c
      title = 'SUMMATION'
      n = 5
      call summation ( n, n, a )
      call summation_inverse ( n, b )
      call r8mat_inverse ( n, a, c )
      call r8mat_is_inverse ( n, a, b, error_ab )
      call r8mat_is_inverse ( n, a, c, error_ac )
      norm_frobenius = r8mat_norm_fro ( n, n, a )
      write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) 
     &  title, n, norm_frobenius, error_ab, error_ac
c
c  TRI_UPPER matrix.
c
      title = 'TRI_UPPER'
      n = 5
      seed = 123456789
      alpha = r8_uniform_01 ( seed )
      alpha = anint ( 50.0D+00 * alpha - 25.0D+00 ) / 5.0D+00
      call tri_upper ( alpha, n, a )
      call tri_upper_inverse ( alpha, n, b )
      call r8mat_inverse ( n, a, c )
      call r8mat_is_inverse ( n, a, b, error_ab )
      call r8mat_is_inverse ( n, a, c, error_ac )
      norm_frobenius = r8mat_norm_fro ( n, n, a )
      write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) 
     &  title, n, norm_frobenius, error_ab, error_ac
c
c  TRIS matrix.
c
      title = 'TRIS'
      n = 5
      seed = 123456789
      alpha = r8_uniform_01 ( seed )
      alpha = anint ( 50.0D+00 * alpha - 25.0D+00 ) / 5.0D+00
      beta = r8_uniform_01 ( seed )
      beta = anint ( 50.0D+00 * beta - 25.0D+00 ) / 5.0D+00
      gamma = r8_uniform_01 ( seed )
      gamma = anint ( 50.0D+00 * gamma - 25.0D+00 ) / 5.0D+00
      call tris ( n, n, alpha, beta, gamma, a )
      call tris_inverse ( n, alpha, beta, gamma, b )
      call r8mat_inverse ( n, a, c )
      call r8mat_is_inverse ( n, a, b, error_ab )
      call r8mat_is_inverse ( n, a, c, error_ac )
      norm_frobenius = r8mat_norm_fro ( n, n, a )
      write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) 
     &  title, n, norm_frobenius, error_ab, error_ac
c
c  TRIV matrix.
c
      title = 'TRIV'
      n = 5
      seed = 123456789
      call r8vec_uniform_01 ( n - 1, seed, x )
      do i = 1, n - 1
        x(i) = anint ( 50.0D+00 * x(i) - 25.0D+00 ) / 5.0D+00
      end do
      call r8vec_uniform_01 ( n, seed, y )
      do i = 1, n
        y(i) = anint ( 50.0D+00 * y(i) - 25.0D+00 ) / 5.0D+00
      end do
      call r8vec_uniform_01 ( n - 1, seed, z )
      do i = 1, n - 1
        z(i) = anint ( 50.0D+00 * z(i) - 25.0D+00 ) / 5.0D+00
      end do
      call triv ( n, x, y, z, a )
      call triv_inverse ( n, x, y, z, b )
      call r8mat_inverse ( n, a, c )
      call r8mat_is_inverse ( n, a, b, error_ab )
      call r8mat_is_inverse ( n, a, c, error_ac )
      norm_frobenius = r8mat_norm_fro ( n, n, a )
      write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) 
     &  title, n, norm_frobenius, error_ab, error_ac
c
c  TRIW matrix.
c
      title = 'TRIW'
      n = 5
      seed = 123456789
      k = i4_uniform ( 0, n - 1, seed )
      alpha = r8_uniform_01 ( seed )
      alpha = anint ( 50.0D+00 * alpha - 25.0D+00 ) / 5.0D+00
      call triw ( alpha, k, n, a )
      call triw_inverse ( alpha, k, n, b )
      call r8mat_inverse ( n, a, c )
      call r8mat_is_inverse ( n, a, b, error_ab )
      call r8mat_is_inverse ( n, a, c, error_ac )
      norm_frobenius = r8mat_norm_fro ( n, n, a )
      write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) 
     &  title, n, norm_frobenius, error_ab, error_ac
c
c  UPSHIFT matrix.
c
      title = 'UPSHIFT'
      n = 5
      call upshift ( n, a )
      call upshift_inverse ( n, b )
      call r8mat_inverse ( n, a, c )
      call r8mat_is_inverse ( n, a, b, error_ab )
      call r8mat_is_inverse ( n, a, c, error_ac )
      norm_frobenius = r8mat_norm_fro ( n, n, a )
      write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) 
     &  title, n, norm_frobenius, error_ab, error_ac
c
c  VAND1 matrix.
c
      title = 'VAND1'
      n = 5
      seed = 123456789
      call r8vec_uniform_01 ( n, seed, x )
      do i = 1, n
        x(i) = anint ( 50.0D+00 * x(i) - 25.0D+00 ) / 5.0D+00
      end do
      call vand1 ( n, x, a )
      call vand1_inverse ( n, x, b )
      call r8mat_inverse ( n, a, c )
      call r8mat_is_inverse ( n, a, b, error_ab )
      call r8mat_is_inverse ( n, a, c, error_ac )
      norm_frobenius = r8mat_norm_fro ( n, n, a )
      write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) 
     &  title, n, norm_frobenius, error_ab, error_ac
c
c  VAND2 matrix.
c
      title = 'VAND2'
      n = 5
      seed = 123456789
      call r8vec_uniform_01 ( n, seed, x )
      do i = 1, n
        x(i) = anint ( 50.0D+00 * x(i) - 25.0D+00 ) / 5.0D+00
      end do
      call vand2 ( n, x, a )
      call vand2_inverse ( n, x, b )
      call r8mat_inverse ( n, a, c )
      call r8mat_is_inverse ( n, a, b, error_ab )
      call r8mat_is_inverse ( n, a, c, error_ac )
      norm_frobenius = r8mat_norm_fro ( n, n, a )
      write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) 
     &  title, n, norm_frobenius, error_ab, error_ac
c
c  WILK03 matrix.
c
      title = 'WILK03'
      n = 3
      call wilk03 ( a )
      call wilk03_inverse ( b )
      call r8mat_inverse ( n, a, c )
      call r8mat_is_inverse ( n, a, b, error_ab )
      call r8mat_is_inverse ( n, a, c, error_ac )
      norm_frobenius = r8mat_norm_fro ( n, n, a )
      write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) 
     &  title, n, norm_frobenius, error_ab, error_ac
c
c  WILK04 matrix.
c
      title = 'WILK04'
      n = 4
      call wilk04 ( a )
      call wilk04_inverse ( b )
      call r8mat_inverse ( n, a, c )
      call r8mat_is_inverse ( n, a, b, error_ab )
      call r8mat_is_inverse ( n, a, c, error_ac )
      norm_frobenius = r8mat_norm_fro ( n, n, a )
      write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) 
     &  title, n, norm_frobenius, error_ab, error_ac
c
c  WILK05 matrix.
c
      title = 'WILK05'
      n = 5
      call wilk05 ( a )
      call wilk05_inverse ( b )
      call r8mat_inverse ( n, a, c )
      call r8mat_is_inverse ( n, a, b, error_ab )
      call r8mat_is_inverse ( n, a, c, error_ac )
      norm_frobenius = r8mat_norm_fro ( n, n, a )
      write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) 
     &  title, n, norm_frobenius, error_ab, error_ac
c
c  WILK21 matrix.
c
      title = 'WILK21'
      n = 21
      call wilk21 ( n, a )
      call wilk21_inverse ( n, b )
      call r8mat_inverse ( n, a, c )
      call r8mat_is_inverse ( n, a, b, error_ab )
      call r8mat_is_inverse ( n, a, c, error_ac )
      norm_frobenius = r8mat_norm_fro ( n, n, a )
      write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) 
     &  title, n, norm_frobenius, error_ab, error_ac
c
c  WILSON matrix.
c
      title = 'WILSON'
      n = 4
      call wilson ( a )
      call wilson_inverse ( b )
      call r8mat_inverse ( n, a, c )
      call r8mat_is_inverse ( n, a, b, error_ab )
      call r8mat_is_inverse ( n, a, c, error_ac )
      norm_frobenius = r8mat_norm_fro ( n, n, a )
      write ( *, '(2x,a20,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) 
     &  title, n, norm_frobenius, error_ab, error_ac

      return
      end
      subroutine test_null ( )

c*********************************************************************72
c
cc TEST_NULL tests the null vectors.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    13 June 2011
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n_max
      parameter ( n_max = 25 )

      double precision a(n_max,n_max)
      double precision at(n_max,n_max)
      double precision alpha
      integer col_num
      double precision error_l2
      double precision f1
      double precision f2
      integer m
      integer mt
      integer n
      integer nt
      double precision norm_a_frobenius
      double precision norm_x_l2
      double precision r8_uniform_01
      double precision r8mat_norm_fro
      double precision r8vec_norm_l2
      integer row_num
      integer seed
      character*20 title
      double precision x(n_max)

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST_NULL'
      write ( *, '(a)' ) '  A = a test matrix of order M by N'
      write ( *, '(a)' ) 
     &  '  x = an N vector, candidate for a null vector.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  ||A|| = Frobenius norm of A.'
      write ( *, '(a)' ) '  ||x|| = L2 norm of x.'
      write ( *, '(a)' ) 
     &  '  ||A*x||/||x|| = L2 norm of A*x over L2 norm of x.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Matrix title	           M     N      ' //
     &         '||A||            ||x||        ||A*x||/||x||'
      write ( *, '(a)' ) ' '
c
c  ARCHIMEDES matrix.
c
      title = 'ARCHIMEDES'
      m = 7
      n = 8
      call archimedes ( a )
      call archimedes_null ( x )
      call r8mat_is_null_vector ( m, n, a, x, error_l2 )
      norm_a_frobenius = r8mat_norm_fro ( m, n, a )
      norm_x_l2 = r8vec_norm_l2 ( n, x )
      write ( *, '(2x,a20,2x,i4,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) 
     &  title, m, n, norm_a_frobenius, norm_x_l2, error_l2 
c
c  CHEBY_DIFf1 matrix.
c
      title = 'CHEBY_DIFF1'
      m = 5
      n = m
      call cheby_diff1 ( n, a )
      call cheby_diff1_null ( n, x )
      call r8mat_is_null_vector ( m, n, a, x, error_l2 )
      norm_a_frobenius = r8mat_norm_fro ( m, n, a )
      norm_x_l2 = r8vec_norm_l2 ( n, x )
      write ( *, '(2x,a20,2x,i4,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) 
     &  title, m, n, norm_a_frobenius, norm_x_l2, error_l2
c
c  CREATION matrix.
c
      title = 'CREATION'
      m = 5
      n = m
      call creation ( m, n, a )
      call creation_null ( n, x )
      call r8mat_is_null_vector ( m, n, a, x, error_l2 )
      norm_a_frobenius = r8mat_norm_fro ( m, n, a )
      norm_x_l2 = r8vec_norm_l2 ( n, x )
      write ( *, '(2x,a20,2x,i4,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) 
     &  title, m, n, norm_a_frobenius, norm_x_l2, error_l2 
c
c  DIF1 matrix.
c  Only has null vectors for N odd.
c
      title = 'DIF1'
      m = 5
      n = m
      call dif1 ( m, n, a )
      call dif1_null ( n, x )
      call r8mat_is_null_vector ( m, n, a, x, error_l2 )
      norm_a_frobenius = r8mat_norm_fro ( m, n, a )
      norm_x_l2 = r8vec_norm_l2 ( n, x )
      write ( *, '(2x,a20,2x,i4,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) 
     &  title, m, n, norm_a_frobenius, norm_x_l2, error_l2 
c
c  DIF1CYCLIC matrix.
c
      title = 'DIF1CYCLIC'
      m = 5
      n = m
      call dif1cyclic ( n, a )
      call dif1cyclic_null ( n, x )
      call r8mat_is_null_vector ( m, n, a, x, error_l2 )
      norm_a_frobenius = r8mat_norm_fro ( m, n, a )
      norm_x_l2 = r8vec_norm_l2 ( n, x )
      write ( *, '(2x,a20,2x,i4,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) 
     &  title, m, n, norm_a_frobenius, norm_x_l2, error_l2 
c
c  DIF2CYCLIC matrix.
c
      title = 'DIF2CYCLIC'
      m = 5
      n = m
      call dif2cyclic ( n, a )
      call dif2cyclic_null ( n, x )
      call r8mat_is_null_vector ( m, n, a, x, error_l2 )
      norm_a_frobenius = r8mat_norm_fro ( m, n, a )
      norm_x_l2 = r8vec_norm_l2 ( n, x )
      write ( *, '(2x,a20,2x,i4,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) 
     &  title, m, n, norm_a_frobenius, norm_x_l2, error_l2 
c
c  EBERLEIN matrix.
c  We have a LEFT null vector.
c
      title = 'EBERLEIN (left)'
      m = 5
      n = 5
      seed = 123456789
      alpha = r8_uniform_01 ( seed )
      alpha = anint ( 50.0D+00 * alpha - 25.0D+00 ) / 5.0D+00
      call eberlein ( alpha, n, a )
      mt = n
      nt = m
      call r8mat_transpose ( m, n, a, at )
      call eberlein_null_left ( n, x )
      call r8mat_is_null_vector ( mt, nt, at, x, error_l2 )
      norm_a_frobenius = r8mat_norm_fro ( mt, nt, at )
      norm_x_l2 = r8vec_norm_l2 ( nt, x )
      write ( *, '(2x,a20,2x,i4,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) 
     &  title, m, n, norm_a_frobenius, norm_x_l2, error_l2 
c
c  FIBONACCI1 matrix.
c
      title = 'FIBONACCI1'
      m = 5
      n = m
      seed = 123456789
      f1 = r8_uniform_01 ( seed )
      f1 = anint ( 50.0D+00 * f1 - 25.0D+00 ) / 5.0D+00
      f2 = r8_uniform_01 ( seed )
      f2 = anint ( 50.0D+00 * f2 - 25.0D+00 ) / 5.0D+00
      call fibonacci1 ( n, f1, f2, a )
      call fibonacci1_null ( n, f1, f2, x )
      call r8mat_is_null_vector ( m, n, a, x, error_l2 )
      norm_a_frobenius = r8mat_norm_fro ( m, n, a )
      norm_x_l2 = r8vec_norm_l2 ( n, x )
      write ( *, '(2x,a20,2x,i4,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) 
     &  title, m, n, norm_a_frobenius, norm_x_l2, error_l2 
c
c  LAUCHLI matrix.
c  We have a LEFT null vector of a RECTANGULAR matrix.
c
      title = 'LAUCHLI (left)'
      m = 6
      n = m - 1
      seed = 123456789
      alpha = r8_uniform_01 ( seed )
      alpha = anint ( 50.0D+00 * alpha - 25.0D+00 ) / 5.0D+00
      call lauchli ( alpha, m, n, a )
      mt = n
      nt = m
      call r8mat_transpose ( m, n, a, at )
      call lauchli_null_left ( alpha, m, n, x )
      call r8mat_is_null_vector ( mt, nt, at, x, error_l2 )
      norm_a_frobenius = r8mat_norm_fro ( mt, nt, at )
      norm_x_l2 = r8vec_norm_l2 ( nt, x )
      write ( *, '(2x,a20,2x,i4,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) 
     &  title, m, n, norm_a_frobenius, norm_x_l2, error_l2 
c
c  LINE_ADJ matrix.
c
      title = 'LINE_ADJ'
      m = 7
      n = m
      call line_adj ( n, a )
      call line_adj_null ( n, x )
      call r8mat_is_null_vector ( m, n, a, x, error_l2 )
      norm_a_frobenius = r8mat_norm_fro ( m, n, a )
      norm_x_l2 = r8vec_norm_l2 ( n, x )
      write ( *, '(2x,a20,2x,i4,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) 
     &  title, m, n, norm_a_frobenius, norm_x_l2, error_l2 
c
c  MOLER2 matrix.
c
      title = 'MOLER2'
      m = 5
      n = 5
      call moler2 ( a )
      call moler2_null ( x )
      call r8mat_is_null_vector ( m, n, a, x, error_l2 )
      norm_a_frobenius = r8mat_norm_fro ( m, n, a )
      norm_x_l2 = r8vec_norm_l2 ( n, x )
      write ( *, '(2x,a20,2x,i4,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) 
     &  title, m, n, norm_a_frobenius, norm_x_l2, error_l2 
c
c  NEUMANN matrix.
c
      title = 'NEUMANN'
      row_num = 5
      col_num = 5
      m = row_num * col_num
      n = row_num * col_num
      call neumann ( row_num, col_num, n, a )
      call neumann_null ( row_num, col_num, x )
      call r8mat_is_null_vector ( m, n, a, x, error_l2 )
      norm_a_frobenius = r8mat_norm_fro ( m, n, a )
      norm_x_l2 = r8vec_norm_l2 ( n, x )
      write ( *, '(2x,a20,2x,i4,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) 
     &  title, m, n, norm_a_frobenius, norm_x_l2, error_l2 
c
c  ONE matrix.
c
      title = 'ONE'
      m = 5
      n = 5
      call one ( m, n, a )
      call one_null ( n, x )
      call r8mat_is_null_vector ( m, n, a, x, error_l2 )
      norm_a_frobenius = r8mat_norm_fro ( m, n, a )
      norm_x_l2 = r8vec_norm_l2 ( n, x )
      write ( *, '(2x,a20,2x,i4,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) 
     &  title, m, n, norm_a_frobenius, norm_x_l2, error_l2 
c
c  RING_ADJ matrix.
c  N must be a multiple of 4 for there to be a null vector.
c
      title = 'RING_ADJ'
      m = 12
      n = 12
      call ring_adj ( n, a )
      call ring_adj_null ( n, x )
      call r8mat_is_null_vector ( m, n, a, x, error_l2 )
      norm_a_frobenius = r8mat_norm_fro ( m, n, a )
      norm_x_l2 = r8vec_norm_l2 ( n, x )
      write ( *, '(2x,a20,2x,i4,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) 
     &  title, m, n, norm_a_frobenius, norm_x_l2, error_l2 
c
c  ROSSER1 matrix.
c
      title = 'ROSSER1'
      m = 8
      n = 8
      call rosser1 ( a )
      call rosser1_null ( x )
      call r8mat_is_null_vector ( m, n, a, x, error_l2 )
      norm_a_frobenius = r8mat_norm_fro ( m, n, a )
      norm_x_l2 = r8vec_norm_l2 ( n, x )
      write ( *, '(2x,a20,2x,i4,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) 
     &  title, m, n, norm_a_frobenius, norm_x_l2, error_l2 
c
c  ZERO matrix.
c
      title = 'ZERO'
      m = 5
      n = 5
      call zero ( m, n, a )
      call zero_null ( n, x )
      call r8mat_is_null_vector ( m, n, a, x, error_l2 )
      norm_a_frobenius = r8mat_norm_fro ( m, n, a )
      norm_x_l2 = r8vec_norm_l2 ( n, x )
      write ( *, '(2x,a20,2x,i4,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) 
     &  title, m, n, norm_a_frobenius, norm_x_l2, error_l2 

      return
      end
      subroutine test_plu ( )

c*********************************************************************72
c
cc TEST_PLU tests the PLU factors.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    26 June 2011
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer m_max
      parameter ( m_max = 25 )
      integer n_max
      parameter ( n_max = 25 )

      double precision a(m_max,n_max)
      double precision alpha
      double precision error_frobenius
      double precision l(m_max,n_max)
      integer m
      integer n
      double precision norm_a_frobenius
      double precision p(m_max,n_max)
      double precision r8_uniform_01
      double precision r8mat_norm_fro
      integer seed
      character * ( 20 ) title
      double precision u(m_max,n_max)

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST_PLU'
      write ( *, '(a)' ) '  A = a test matrix of order M by N'
      write ( *, '(a)' ) '  P, L, U are the PLU factors.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  ||A|| = Frobenius norm of A.'
      write ( *, '(a)' ) '  ||A-PLU|| = Frobenius norm of A-P*L*U.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Matrix title	           M     N      ' 
     &  // '||A||            ||A-PLU||'
      write ( *, '(a)' ) ' '
c
c  BODEWIG matrix.
c
      title = 'BODEWIG'
      m = 4
      n = 4
      call bodewig ( a )
      call bodewig_plu ( p, l, u )
      call r8mat_is_plu ( m, n, a, p, l, u, error_frobenius )
      norm_a_frobenius = r8mat_norm_fro ( m, n, a )
      write ( *, '(2x,a20,2x,i4,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) 
     &  title, m, n, norm_a_frobenius, error_frobenius
c
c  BORDERBAND matrix.
c
      title = 'BORDERBAND'
      m = 5
      n = 5
      call borderband ( n, a )
      call borderband_plu ( n, p, l, u )
      call r8mat_is_plu ( m, n, a, p, l, u, error_frobenius )
      norm_a_frobenius = r8mat_norm_fro ( m, n, a )
      write ( *, '(2x,a20,2x,i4,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) 
     &  title, m, n, norm_a_frobenius, error_frobenius
c
c  DIF2 matrix.
c
      title = 'DIF2'
      m = 5
      n = 5
      call dif2 ( m, n, a )
      call dif2_plu ( n, p, l, u )
      call r8mat_is_plu ( m, n, a, p, l, u, error_frobenius )
      norm_a_frobenius = r8mat_norm_fro ( m, n, a )
      write ( *, '(2x,a20,2x,i4,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) 
     &  title, m, n, norm_a_frobenius, error_frobenius
c
c  GFPP matrix.
c
      title = 'GFPP'
      m = 5
      n = 5
      seed = 123456789
      alpha = r8_uniform_01 ( seed )
      alpha = anint ( 50.0D+00 * alpha ) / 5.0D+00
      call gfpp ( n, alpha, a )
      call gfpp_plu ( n, alpha, p, l, u )
      call r8mat_is_plu ( m, n, a, p, l, u, error_frobenius )
      norm_a_frobenius = r8mat_norm_fro ( m, n, a )
      write ( *, '(2x,a20,2x,i4,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) 
     &  title, m, n, norm_a_frobenius, error_frobenius
c
c  GIVENS matrix.
c
      title = 'GIVENS'
      m = 5
      n = 5
      call givens ( n, n, a )
      call givens_plu ( n, p, l, u  )
      call r8mat_is_plu ( m, n, a, p, l, u, error_frobenius )
      norm_a_frobenius = r8mat_norm_fro ( m, n, a )
      write ( *, '(2x,a20,2x,i4,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) 
     &  title, m, n, norm_a_frobenius, error_frobenius
c
c  KMS matrix.
c
      title = 'KMS'
      m = 5
      n = 5
      seed = 123456789
      alpha = r8_uniform_01 ( seed )
      alpha = anint ( 50.0D+00 * alpha ) / 5.0D+00
      call kms ( alpha, m, n, a )
      call kms_plu ( alpha, n, p, l, u  )
      call r8mat_is_plu ( m, n, a, p, l, u, error_frobenius )
      norm_a_frobenius = r8mat_norm_fro ( m, n, a )
      write ( *, '(2x,a20,2x,i4,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) 
     &  title, m, n, norm_a_frobenius, error_frobenius
c
c  MAXIJ matrix.
c
      title = 'MAXIJ'
      m = 5
      n = 5
      call maxij ( n, n, a )
      call maxij_plu ( n, p, l, u  )
      call r8mat_is_plu ( m, n, a, p, l, u, error_frobenius )
      norm_a_frobenius = r8mat_norm_fro ( m, n, a )
      write ( *, '(2x,a20,2x,i4,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) 
     &  title, m, n, norm_a_frobenius, error_frobenius
c
c  MINIJ matrix.
c
      title = 'MINIJ'
      m = 5
      n = 5
      call minij ( n, n, a )
      call minij_plu ( n, p, l, u  )
      call r8mat_is_plu ( m, n, a, p, l, u, error_frobenius )
      norm_a_frobenius = r8mat_norm_fro ( m, n, a )
      write ( *, '(2x,a20,2x,i4,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) 
     &  title, m, n, norm_a_frobenius, error_frobenius
c
c  MOLER1 matrix.
c
      title = 'MOLER1'
      m = 5
      n = 5
      seed = 123456789
      alpha = r8_uniform_01 ( seed )
      alpha = anint ( 50.0D+00 * alpha ) / 5.0D+00
      call moler1 ( alpha, n, n, a )
      call moler1_plu ( alpha, n, p, l, u  )
      call r8mat_is_plu ( m, n, a, p, l, u, error_frobenius )
      norm_a_frobenius = r8mat_norm_fro ( m, n, a )
      write ( *, '(2x,a20,2x,i4,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) 
     &  title, m, n, norm_a_frobenius, error_frobenius
c
c  MOLER3 matrix.
c
      title = 'MOLER3'
      m = 5
      n = 5
      call moler3 ( m, n, a )
      call moler3_plu ( n, p, l, u  )
      call r8mat_is_plu ( m, n, a, p, l, u, error_frobenius )
      norm_a_frobenius = r8mat_norm_fro ( m, n, a )
      write ( *, '(2x,a20,2x,i4,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) 
     &  title, m, n, norm_a_frobenius, error_frobenius
c
c  OTO matrix.
c
      title = 'OTO'
      m = 5
      n = 5
      call oto ( m, n, a )
      call oto_plu ( n, p, l, u  )
      call r8mat_is_plu ( m, n, a, p, l, u, error_frobenius )
      norm_a_frobenius = r8mat_norm_fro ( m, n, a )
      write ( *, '(2x,a20,2x,i4,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) 
     &  title, m, n, norm_a_frobenius, error_frobenius
c
c  PASCAL2 matrix.
c
      title = 'PASCAL2'
      m = 5
      n = 5
      call pascal2 ( n, a )
      call pascal2_plu ( n, p, l, u  )
      call r8mat_is_plu ( m, n, a, p, l, u, error_frobenius )
      norm_a_frobenius = r8mat_norm_fro ( m, n, a )
      write ( *, '(2x,a20,2x,i4,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) 
     &  title, m, n, norm_a_frobenius, error_frobenius
c
c  WILSON matrix.
c
      title = 'WILSON'
      m = 4
      n = 4
      call wilson ( a )
      call wilson_plu ( p, l, u  )
      call r8mat_is_plu ( m, n, a, p, l, u, error_frobenius )
      norm_a_frobenius = r8mat_norm_fro ( m, n, a )
      write ( *, '(2x,a20,2x,i4,2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) 
     &  title, m, n, norm_a_frobenius, error_frobenius

      return
      end
      subroutine test_solution ( )

c*********************************************************************72
c
cc TEST_SOLUTION tests the linear solution computations.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    12 June 2011
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n_max
      parameter ( n_max = 20 )

      double precision a(n_max,n_max)
      double precision alpha
      double precision b(n_max,n_max)
      double precision beta
      double precision error_frobenius
      double precision gamma
      integer i1
      integer i4_uniform
      integer k
      integer m
      integer n
      integer ncol
      double precision norm_frobenius
      integer nrow
      double precision r8_uniform_01
      double precision r8mat_norm_fro
      integer seed
      integer seed_save
      character*20 title
      double precision x(n_max)

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST_SOLUTION'
      write ( *, '(a)' ) 
     &  '  Compute the Frobenius norm of the solution error:'
      write ( *, '(a)' ) '    A * X - B'
      write ( *, '(a)' ) 
     &  '  given MxN matrix A, NxK solution X, MxK right hand side B.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 
     &  '  Matrix title             M     N' // 
     &  '     K      ||A||         ||A*X-B||'
      write ( *, '(a)' ) ' '
c
c  BODEWIG matrix.
c
      title = 'BODEWIG'
      m = 4
      n = 4
      k = 1
      call bodewig ( a )
      call bodewig_rhs ( b )
      call bodewig_solution ( x )
      call r8mat_is_solution ( m, n, k, a, x, b, error_frobenius )
      norm_frobenius = r8mat_norm_fro ( n, n, a )
      write ( *, '(2x,a20,2x,i4,2x,i4,2x,i4,2x,g14.6,2x,g14.6)' ) 
     &  title, m, n, k, norm_frobenius, error_frobenius
c
c  DIF2 matrix.
c
      title = 'DIF2'
      m = 10
      n = 10
      k = 2
      call dif2 ( m, n, a )
      call dif2_rhs ( m, k, b )
      call dif2_solution ( n, k, x )
      call r8mat_is_solution ( m, n, k, a, x, b, error_frobenius )
      norm_frobenius = r8mat_norm_fro ( n, n, a )
      write ( *, '(2x,a20,2x,i4,2x,i4,2x,i4,2x,g14.6,2x,g14.6)' ) 
     &  title, m, n, k, norm_frobenius, error_frobenius
c
c  FRANK matrix.
c
      title = 'FRANK'
      m = 10
      n = 10
      k = 2
      call frank ( n, a )
      call frank_rhs ( m, k, b )
      call frank_solution ( n, k, x )
      call r8mat_is_solution ( m, n, k, a, x, b, error_frobenius )
      norm_frobenius = r8mat_norm_fro ( n, n, a )
      write ( *, '(2x,a20,2x,i4,2x,i4,2x,i4,2x,g14.6,2x,g14.6)' ) 
     &  title, m, n, k, norm_frobenius, error_frobenius
c
c  POISSON matrix.
c
      title = 'POISSON'
      nrow = 4
      ncol = 5
      m = nrow * ncol
      n = nrow * ncol
      k = 1
      call poisson ( nrow, ncol, n, a )
      call poisson_rhs ( nrow, ncol, n, b )
      call poisson_solution ( nrow, ncol, n, x )
      call r8mat_is_solution ( m, n, k, a, x, b, error_frobenius )
      norm_frobenius = r8mat_norm_fro ( n, n, a )
      write ( *, '(2x,a20,2x,i4,2x,i4,2x,i4,2x,g14.6,2x,g14.6)' ) 
     &  title, m, n, k, norm_frobenius, error_frobenius
c
c  WILK03 matrix.
c
      title = 'WILK03'
      m = 3
      n = 3
      k = 1
      call wilk03 ( a )
      call wilk03_rhs ( b )
      call wilk03_solution ( x )
      call r8mat_is_solution ( m, n, k, a, x, b, error_frobenius )
      norm_frobenius = r8mat_norm_fro ( n, n, a )
      write ( *, '(2x,a20,2x,i4,2x,i4,2x,i4,2x,g14.6,2x,g14.6)' ) 
     &  title, m, n, k, norm_frobenius, error_frobenius
c
c  WILK04 matrix.
c
      title = 'WILK04'
      m = 4
      n = 4
      k = 1
      call wilk04 ( a )
      call wilk04_rhs ( b )
      call wilk04_solution ( x )
      call r8mat_is_solution ( m, n, k, a, x, b, error_frobenius )
      norm_frobenius = r8mat_norm_fro ( n, n, a )
      write ( *, '(2x,a20,2x,i4,2x,i4,2x,i4,2x,g14.6,2x,g14.6)' ) 
     &  title, m, n, k, norm_frobenius, error_frobenius
c
c  WILSON matrix.
c
      title = 'WILSON'
      m = 4
      n = 4
      k = 1
      call wilson ( a )
      call wilson_rhs ( b )
      call wilson_solution ( x )
      call r8mat_is_solution ( m, n, k, a, x, b, error_frobenius )
      norm_frobenius = r8mat_norm_fro ( n, n, a )
      write ( *, '(2x,a20,2x,i4,2x,i4,2x,i4,2x,g14.6,2x,g14.6)' ) 
     &  title, m, n, k, norm_frobenius, error_frobenius

      return
      end
