      subroutine mexp_a ( test, n, a )

c*********************************************************************72
c
cc MEXP_A returns the matrix for a given test.
c
c  Discussion:
c
c     1) Diagonal example
c     2) Symmetric example
c     3) Laub
c     4) Moler and Van Loan
c     5) Moler and Van Loan
c     6) Moler and Van Loan
c     7) Moler and Van Loan
c     8) Wikipedia example
c     9) NAG F01ECF
c    10) Ward #1
c    11) Ward #2
c    12) Ward #3
c    13) Ward #4
c    14) Moler example
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
c  Reference:
c
c    Alan Laub,
c    Review of "Linear System Theory" by Joao Hespanha,
c    SIAM Review,
c    Volume 52, Number 4, December 2010, page 779-781.
c
c    Cleve Moler, Charles VanLoan,
c    Nineteen Dubious Ways to Compute the Exponential of a Matrix,
c    Twenty-Five Years Later,
c    SIAM Review,
c    Volume 45, Number 1, March 2003, pages 3-49.
c
c    Cleve Moler,
c    Cleve's Corner: A Balancing Act for the Matrix Exponential,
c    July 23rd, 2012.
c
c    Robert Ward,
c    Numerical computation of the matrix exponential with accuracy estimate,
c    SIAM Journal on Numerical Analysis,
c    Volume 14, Number 4, September 1977, pages 600-610.
c
c  Parameters:
c
c    Input, integer TEST, the index of the test case.
c
c    Input, integer N, the order of the matrix.
c
c    Output, double precision A(N,N), the matrix.
c
      implicit none

      integer n

      double precision a(n,n)
      double precision a01(2,2)
      double precision a02(2,2)
      double precision a03(2,2)
      double precision a04(2,2)
      double precision a05(4,4)
      double precision a06(2,2)
      double precision a08(3,3)
      double precision a09(4,4)
      double precision a10(3,3)
      double precision a11(3,3)
      double precision a12(3,3)
      integer i
      integer j
      integer test
      double precision r8_epsilon

      save a01
      save a02
      save a03
      save a04
      save a05
      save a06
      save a08
      save a09
      save a10
      save a11
      save a12

      data a01 /
     &    1.0D+00, 0.0D+00, 
     &     0.0D+00, 2.0D+00 /
      data a02 /
     &     1.0D+00, 3.0D+00, 
     &     3.0D+00, 2.0D+00 /
      data a03 /
     &     0.0D+00, -39.0D+00, 
     &     1.0D+00, -40.0D+00 /
      data a04 /
     &     -49.0D+00, -64.0D+00, 
     &      24.0D+00,  31.0D+00 /
      data a05 /
     &     0.0D+00, 0.0D+00, 0.0D+00, 0.0D+00, 
     &     6.0D+00, 0.0D+00, 0.0D+00, 0.0D+00, 
     &     0.0D+00, 6.0D+00, 0.0D+00, 0.0D+00, 
     &     0.0D+00, 0.0D+00, 6.0D+00, 0.0D+00 /
      data a06 /
     &     1.0D+00, 0.0D+00, 
     &     1.0D+00, 1.0D+00 /
      data a08 /
     &     21.0D+00,  -5.0D+00,   4.0D+00, 
     &     17.0D+00,  -1.0D+00,   4.0D+00, 
     &      6.0D+00,  -6.0D+00,  16.0D+00 /
      data a09 /
     &     1.0D+00, 3.0D+00, 3.0D+00, 3.0D+00, 
     &     2.0D+00, 1.0D+00, 2.0D+00, 3.0D+00, 
     &     2.0D+00, 1.0D+00, 1.0D+00, 3.0D+00, 
     &     2.0D+00, 2.0D+00, 2.0D+00, 1.0D+00 /
      data a10 /
     &     4.0D+00, 1.0D+00, 1.0D+00, 
     &     2.0D+00, 4.0D+00, 1.0D+00, 
     &     0.0D+00, 1.0D+00, 4.0D+00 /
      data a11 /
     &     29.87942128909879D+00, 
     &      0.7815750847907159D+00, 
     &     -2.289519314033932D+00, 
     &      0.7815750847907159D+00, 
     &     25.72656945571064D+00, 
     &      8.680737820540137D+00, 
     &     -2.289519314033932D+00, 
     &      8.680737820540137D+00, 
     &     34.39400925519054D+00 /
      data a12 /
     &    -131.0D+00, -390.0D+00, -387.0D+00, 
     &      19.0D+00,   56.0D+00,   57.0D+00, 
     &      18.0D+00,   54.0D+00,   52.0D+00 /

      if ( test .eq. 1 ) then
        call r8mat_copy ( n, n, a01, a )
      else if ( test .eq. 2 ) then
        call r8mat_copy ( n, n, a02, a )
      else if ( test .eq. 3 ) then
        call r8mat_copy ( n, n, a03, a )
      else if ( test .eq. 4 ) then
        call r8mat_copy ( n, n, a04, a )
      else if ( test .eq. 5 ) then
        call r8mat_copy ( n, n, a05, a )
      else if ( test .eq. 6 ) then
        call r8mat_copy ( n, n, a06, a )
      else if ( test .eq. 7 ) then
        a(1,1) = 1.0D+00 + r8_epsilon ( )
        a(2,1) = 0.0D+00
        a(1,2) = 1.0D+00
        a(2,2) = 1.0D+00 - r8_epsilon ( )
      else if ( test .eq. 8 ) then
        call r8mat_copy ( n, n, a08, a )
      else if ( test .eq. 9 ) then
        call r8mat_copy ( n, n, a09, a )
      elseif ( test .eq. 10 ) then
        call r8mat_copy ( n, n, a10, a )
      elseif ( test .eq. 11 ) then
        call r8mat_copy ( n, n, a11, a )
      elseif ( test .eq. 12 ) then
        call r8mat_copy ( n, n, a12, a )
      elseif ( test .eq. 13 ) then
        do j = 1, n
          do i = 1, n
            if ( j .eq. i + 1 ) then
              a(i,j) = 1.0D+00
            else if ( i .eq. n .and. j .eq. 1 ) then
              a(i,j) = 1.0D-10
            else
              a(i,j) = 0.0D+00
            end if
          end do
        end do
      elseif ( test .eq. 14 ) then
        a(1,1) = 0.0D+00
        a(1,2) = 1.0D-08
        a(1,3) = 0.0D+00
        a(2,1) = - 2.0D+10 - 2.0D+08 / 3.0D+00
        a(2,2) = - 3.0D+00
        a(2,3) = 2.0D+10
        a(3,1) = 200.0D+00 / 3.0D+00
        a(3,2) = 0.0D+00
        a(3,3) = - 200.0D+00 / 3.0D+00
      else
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'MEXP_A - Fatal error!'
        write ( *, '(a,i8)' ) '  Illegal value of TEST = ', test
      end if

      return
      end
      subroutine mexp_expa ( test, n, expa )

c*********************************************************************72
c
cc MEXP_EXPA returns the "exact" exponential matrix for a given test.
c
c  Discussion:
c
c    In some cases, the "exact" value is given to six significant digits.
c
c     1) Diagonal example
c     2) Symmetric example
c     3) Laub
c     4) Moler and Van Loan
c     5) Moler and Van Loan
c     6) Moler and Van Loan
c     7) Moler and Van Loan
c     8) Wikipedia example
c     9) NAG F01ECF
c    10) Ward #1
c    11) Ward #2
c    12) Ward #3
c    13) Ward #4
c    14) Moler example
c
c    Thanks to Alex Griffing for correcting the value of matrix 3,
c    17 October 2012.
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
c  Reference:
c
c    Alan Laub,
c    Review of "Linear System Theory" by Joao Hespanha,
c    SIAM Review,
c    Volume 52, Number 4, December 2010, page 779-781.
c
c    Cleve Moler, Charles VanLoan,
c    Nineteen Dubious Ways to Compute the Exponential of a Matrix,
c    Twenty-Five Years Later,
c    SIAM Review,
c    Volume 45, Number 1, March 2003, pages 3-49.
c
c    Cleve Moler,
c    Cleve's Corner: A Balancing Act for the Matrix Exponential,
c    July 23rd, 2012.
c
c    Robert Ward,
c    Numerical computation of the matrix exponential with accuracy estimate,
c    SIAM Journal on Numerical Analysis,
c    Volume 14, Number 4, September 1977, pages 600-610.
c
c  Parameters:
c
c    Input, integer TEST, the index of the test case.
c
c    Input, integer N, the order of the matrix.
c
c    Output, double precision EXPA(N,N), the exponential of the test matrix.
c
      implicit none

      integer n

      double precision exp16
      double precision exp4
      double precision expa(n,n)
      double precision expa01(2,2)
      double precision expa02(2,2)
      double precision expa03(2,2)
      double precision expa04(2,2)
      double precision expa05(4,4)
      double precision expa06(2,2)
      double precision expa07(2,2)
      double precision expa09(4,4)
      double precision expa10(3,3)
      double precision expa11(3,3)
      double precision expa12(3,3)
      double precision expa14(3,3)
      integer i
      integer j
      integer test

      save expa01
      save expa02
      save expa03
      save expa04
      save expa05
      save expa06
      save expa07
      save expa09
      save expa10
      save expa11
      save expa12
      save expa14

      data expa01 /
     &   2.718281828459046D+00, 0.0D+00, 
     &    0.0D+00,               7.389056098930650D+00 /
      data expa02 /
     &     39.322809708033859D+00,  46.166301438885753D+00, 
     &     46.166301438885768D+00,  54.711576854329110D+00 /
      data expa03 /
     &  0.37756048D+00,  0.00968104D+00,
     & -0.37756048D+00, -0.00968104D+00 /
      data expa04 /
     &     -0.735759D+00, -1.471518D+00, 
     &      0.551819D+00,  1.103638D+00 /
      data expa05 /
     &     1.0D+00,  0.0D+00, 0.0D+00, 0.0D+00, 
     &     6.0D+00,  1.0D+00, 0.0D+00, 0.0D+00, 
     &    18.0D+00,  6.0D+00, 1.0D+00, 0.0D+00, 
     &    36.0D+00, 18.0D+00, 6.0D+00, 1.0D+00 /
      data expa06 /
     &     2.718281828459046D+00, 0.0D+00, 
     &     2.718281828459046D+00, 2.718281828459046D+00 /
      data expa07 /
     &     2.718309D+00, 0.0D+00, 
     &     2.718282D+00, 2.718255D+00 /
      data expa09 /
     &     740.7038D+00, 731.2510D+00, 823.7630D+00, 998.4355D+00, 
     &     610.8500D+00, 603.5524D+00, 679.4257D+00, 823.7630D+00, 
     &     542.2743D+00, 535.0884D+00, 603.5524D+00, 731.2510D+00, 
     &     549.1753D+00, 542.2743D+00, 610.8500D+00, 740.7038D+00 /
      data expa10 /
     &     147.8666224463699D+00, 
     &     127.7810855231823D+00, 
     &     127.7810855231824D+00, 
     &     183.7651386463682D+00, 
     &     183.7651386463682D+00, 
     &     163.6796017231806D+00, 
     &     71.79703239999647D+00, 
     &     91.88256932318415D+00, 
     &    111.9681062463718D+00 /
      data expa11 /
     &    5.496313853692378D+15, 
     &   -1.823188097200899D+16, 
     &   -3.047577080858001D+16, 
     &   -1.823188097200898D+16, 
     &    6.060522870222108D+16, 
     &    1.012918429302482D+17, 
     &   -3.047577080858001D+16, 
     &    1.012918429302482D+17, 
     &    1.692944112408493D+17 /
      data expa12 /
     &   -1.509644158793135D+00, 
     &   -5.632570799891469D+00, 
     &   -4.934938326088363D+00, 
     &    0.3678794391096522D+00, 
     &    1.471517758499875D+00, 
     &    1.103638317328798D+00, 
     &    0.1353352811751005D+00, 
     &    0.4060058435250609D+00, 
     &    0.5413411267617766D+00 /
      data expa14 /
     &    4.468494682831735D-01, 
     &   -5.743067779479621D+06, 
     &    4.477229778494929D-01, 
     &    1.540441573839520D-09, 
     &   -1.528300386868247D-02, 
     &    1.542704845195912D-09, 
     &    4.628114535587735D-01, 
     &   -4.526542712784168D+06, 
     &    4.634806488376499D-01 /

      if ( test .eq. 1 ) then
        call r8mat_copy ( n, n, expa01, expa )
      else if ( test .eq. 2 ) then
        call r8mat_copy ( n, n, expa02, expa )
      else if ( test .eq. 3 ) then
        call r8mat_copy ( n, n, expa03, expa )
      else if ( test .eq. 4 ) then
        call r8mat_copy ( n, n, expa04, expa )
      else if ( test .eq. 5 ) then
        call r8mat_copy ( n, n, expa05, expa )
      else if ( test .eq. 6 ) then
        call r8mat_copy ( n, n, expa06, expa )
      else if ( test .eq. 7 ) then
        call r8mat_copy ( n, n, expa07, expa )
      else if ( test .eq. 8 ) then
        exp16 = exp ( 16.0D+00 )
        exp4 = exp ( 4.0D+00 )
        expa(1,1) = 0.25D+00 * ( 13.0D+00 * exp16           - exp4 )
        expa(2,1) = 0.25D+00 * ( -9.0D+00 * exp16           + exp4 )
        expa(3,1) = 0.25D+00 *   16.0D+00 * exp16
        expa(1,2) = 0.25D+00 * ( 13.0D+00 * exp16 - 5.0D+00 * exp4 )
        expa(2,2) = 0.25D+00 * ( -9.0D+00 * exp16 + 5.0D+00 * exp4 )
        expa(2,3) = 0.25D+00 *   16.0D+00 * exp16
        expa(1,3) = 0.25D+00 * (  2.0D+00 * exp16 - 2.0D+00 * exp4 )
        expa(2,3) = 0.25D+00 * ( -2.0D+00 * exp16 + 2.0D+00 * exp4 )
        expa(3,3) = 0.25D+00 *    4.0D+00 * exp16 
      else if ( test .eq. 9 ) then
        call r8mat_copy ( n, n, expa09, expa )
      elseif ( test .eq. 10 ) then
        call r8mat_copy ( n, n, expa10, expa )
      elseif ( test .eq. 11 ) then
        call r8mat_copy ( n, n, expa11, expa )
      elseif ( test .eq. 12 ) then
        call r8mat_copy ( n, n, expa12, expa )
      elseif ( test .eq. 13 ) then
        do j = 1, n
          do i = 1, n
            expa(i,j) = 0.0D+00
          end do
        end do
      elseif ( test .eq. 14 ) then
        call r8mat_copy ( n, n, expa14, expa )
      else
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'MEXP_EXPA - Fatal error!'
        write ( *, '(a,i8)' ) '  Illegal value of TEST = ', test
      end if

      return
      end
      subroutine mexp_n ( test, n )

c*********************************************************************72
c
cc MEXP_N returns the matrix order for a given test.
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
c  Parameters:
c
c    Input, integer TEST, the index of the test case.
c
c    Output, integer N, the order of the matrix.
c
      implicit none

      integer n
      integer test

      if ( test .eq. 1 ) then
        n = 2
      else if ( test .eq. 2 ) then
        n = 2
      else if ( test .eq. 3 ) then
        n = 2
      else if ( test .eq. 4 ) then
        n = 2
      else if ( test .eq. 5 ) then
        n = 4
      else if ( test .eq. 6 ) then
        n = 2
      else if ( test .eq. 7 ) then
        n = 2
      else if ( test .eq. 8 ) then
        n = 3
      else if ( test .eq. 9 ) then
        n = 4
      else if ( test .eq. 10 ) then
        n = 3
      else if ( test .eq. 11 ) then
        n = 3
      else if ( test .eq. 12 ) then
        n = 3
      else if ( test .eq. 13 ) then
        n = 10
      else if ( test .eq. 14 ) then
        n = 3
      else
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'MEXP_N - Fatal error!'
        write ( *, '(a,i8)' ) '  Illegal value of TEST = ', test
      end if

      return
      end
      subroutine mexp_story ( test )

c*********************************************************************72
c
cc MEXP_STORY prints explanatory text for each problem.
c
c  Discussion:
c
c     1) Diagonal example
c     2) Symmetric example
c     3) Laub
c     4) Moler and Van Loan
c     5) Moler and Van Loan
c     6) Moler and Van Loan
c     7) Moler and Van Loan
c     8) Wikipedia example
c     9) NAG F01ECF
c    10) Ward #1
c    11) Ward #2
c    12) Ward #3
c    13) Ward #4
c    14) Moler example
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
c  Reference:
c
c    Alan Laub,
c    Review of "Linear System Theory" by Joao Hespanha,
c    SIAM Review,
c    Volume 52, Number 4, December 2010, page 779-781.
c
c    Cleve Moler, Charles VanLoan,
c    Nineteen Dubious Ways to Compute the Exponential of a Matrix,
c    Twenty-Five Years Later,
c    SIAM Review,
c    Volume 45, Number 1, March 2003, pages 3-49.
c
c    Cleve Moler,
c    Cleve's Corner: A Balancing Act for the Matrix Exponential,
c    July 23rd, 2012.
c
c    Robert Ward,
c    Numerical computation of the matrix exponential with accuracy estimate,
c    SIAM Journal on Numerical Analysis,
c    Volume 14, Number 4, September 1977, pages 600-610.
c
c  Parameters:
c
c    Input, integer TEST, the index of the test case.
c
      implicit none

      integer test

      if ( test .eq. 1 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) '  This matrix is diagonal.'
        write ( *, '(a)' ) 
     &    '  The calculation of the matrix exponential is simple.'
      else if ( test .eq. 2 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) '  This matrix is symmetric.'
        write ( *, '(a)' ) 
     &   '  The calculation is straightforward.'
      else if ( test .eq. 3 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) '  This example is due to Laub.'
        write ( *, '(a)' ) 
     &    '  This matrix is ill-suited for the Taylor series approach.'
        write ( *, '(a)' ) 
     &    '  As powers of A are computed, the entries blow up quickly.'
      else if ( test .eq. 4 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 
     &    '  This example is due to Moler and Van Loan.'
        write ( *, '(a)' ) 
     &    '  The example causes problems for series summation,'
        write ( *, '(a)' ) '  and for diagonal Pade approximations.'
      else if ( test .eq. 5 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 
     &    '  This example is due to Moler and Van Loan.'
        write ( *, '(a)' ) '  This matrix is strictly upper triangular'
        write ( *, '(a)' ) 
     &    '  All powers of A are zero beyond some (low) limit.' 
        write ( *, '(a)' ) 
     &    '  This example will cause problems for Pade approximations.'
      else if ( test .eq. 6 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 
     &    '  This example is due to Moler and Van Loan.'
        write ( *, '(a)' ) 
     &    '  This matrix does not have a complete set of eigenvectors.'
        write ( *, '(a)' ) 
     &    '  That means the eigenvector approach will fail.'
      else if ( test .eq. 7 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 
     &    '  This example is due to Moler and Van Loan.'
        write ( *, '(a)' ) '  This matrix is very close to example 5.'
        write ( *, '(a)' ) 
     &    '  Mathematically, it has a complete set of eigenvectors.'
        write ( *, '(a)' ) 
     &    '  Numerically, however, the calculation will be suspect.'
      else if ( test .eq. 8 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) '  This matrix was an example in Wikipedia.'
      else if ( test .eq. 9 ) then
        write ( *, '(a)' ) ' ' 
        write ( *, '(a)' ) '  This matrix is due to the NAG Library.'
        write ( *, '(a)' ) '  It is an example for function F01ECF.'
      else if ( test .eq. 10 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) '  This is Ward''s example #1.'
        write ( *, '(a)' ) '  It is defective and nonderogatory.'
        write ( *, '(a)' ) '  The eigenvalues are 3, 3 and 6.'
      else if ( test .eq. 11 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) '  This is Ward''s example #2.'
        write ( *, '(a)' ) '  It is a symmetric matrix.'
        write ( *, '(a)' ) '  The eigenvalues are 20, 30, 40.'
      else if ( test .eq. 12 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) '  This is Ward''s example #3.'
        write ( *, '(a)' ) 
     &    '  Ward''s algorithm has difficulty estimating the accuracy'
        write ( *, '(a)' ) 
     &    '  of its results.  The eigenvalues are -1, -2, -20.'
      else if ( test .eq. 13 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) '  This is Ward''s example #4.'
        write ( *, '(a)' ) '  This is a version of the Forsythe matrix.'
        write ( *, '(a)' ) 
     &  '  The eigenvector problem is badly conditioned.'
        write ( *, '(a)' ) 
     &    '  Ward''s algorithm has difficulty estimating the accuracy'
        write ( *, '(a)' ) '  of its results for this problem.'
      else if ( test == 14 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) '  This is Moler''s example.'
        write ( *, '(a)' ) '  This badly scaled matrix caused ' //
     &    'problems for MATLAB''s expm().'
      else
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'MEXP_STORY - Fatal error!'
        write ( *, '(a,i8)' ) '  Illegal value of TEST = ', test
      end if

      return
      end
      subroutine mexp_test_num ( test_num )

c*********************************************************************72
c
cc MEXP_TEST_NUM returns the number of matrix exponential tests.
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
c  Parameters:
c
c    Output, integer TEST_NUM, the number of tests.
c
      implicit none

      integer test_num

      test_num = 14

      return
      end
