      program main

c*********************************************************************72
c
cc MAIN is the main program for COMBO_PRB.
c
c  Discussion:
c
c    COMBO_PRB calls the COMBO tests.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    21 January 2007
c
c  Author:
c
c    John Burkardt
c
      implicit none

      call timestamp ( )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'COMBO_PRB'
      write ( *, '(a)' ) '  FORTRAN77 version'
      write ( *, '(a)' ) '  Test the COMBO library.'

      call test01 ( )
      call test02 ( )
      call test03 ( )
      call test04 ( )
      call test05 ( )
      call test06 ( )
      call test07 ( )
      call test08 ( )
      call test09 ( )

      call test10 ( )
      call test11 ( )
      call test12 ( )
      call test13 ( )
      call test14 ( )
      call test15 ( )
      call test16 ( )
      call test17 ( )
      call test18 ( )
      call test19 ( )

      call test20 ( )
      call test21 ( )
      call test22 ( )
      call test23 ( )
      call test24 ( )
      call test25 ( )
      call test26 ( )
      call test27 ( )
      call test28 ( )
      call test29 ( )

      call test30 ( )
      call test31 ( )
      call test32 ( )
      call test33 ( )
      call test34 ( )
      call test35 ( )
      call test36 ( )
      call test37 ( )
      call test38 ( )
      call test39 ( )

      call test40 ( )
      call test41 ( )
      call test42 ( )
c
c  Terminate.
c
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'COMBO_PRB'
      write ( *, '(a)' ) '  Normal end of execution.'

      write ( *, '(a)' )  ' '
      call timestamp ( )

      stop
      end 
      subroutine test01 ( )

c*********************************************************************72
c
cc TEST01 tests BAL_SEQ_ENUM, BAL_SEQ_RANK, BAL_SEQ_SUCCESSOR, BAL_SEQ_UNRANK.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    06 July 2007
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n
      parameter ( n = 5 )

      integer nseq
      integer rank
      integer rank_old
      integer t(2*n)

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST01'
      write ( *, '(a)' ) '  Balanced sequences:'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  BAL_SEQ_ENUM enumerates,'
      write ( *, '(a)' ) '  BAL_SEQ_RANK ranks,'
      write ( *, '(a)' ) '  BAL_SEQ_SUCCESSOR lists,'
      write ( *, '(a)' ) '  BAL_SEQ_UNRANK unranks.'
c
c  Enumerate.
c
      call bal_seq_enum ( n, nseq )

      write ( *, '(a)' ) ' '
      write ( *, '(a,i8)' ) '  For N = ', n
      write ( *, '(a,i8)' ) 
     &  '  the number of balanced sequences is ', nseq
      write ( *, '(a)' ) ' '
c
c  List.
c
      rank = -1

10    continue

        rank_old = rank

        call bal_seq_successor ( n, t, rank )

        if ( rank .le. rank_old ) then
          go to 20
        end if

        write ( *, '(4x,i3,2x,10i2)' ) rank, t(1:2*n)

      go to 10

20    continue
c
c  Unrank.
c
      rank = nseq / 2

      call bal_seq_unrank ( rank, n, t )

      write ( *, '(a)' ) ' '
      write ( *, '(a,i8)' ) '  The element of rank ', rank
      write ( *, '(a)' ) ' '
      write ( *, '(4x,10i2)' ) t(1:2*n)
c
c  Rank.
c
      call bal_seq_rank ( n, t, rank )

      call i4vec_print ( 2*n, t, '  Element to be ranked:' )

      write ( *, '(a)' ) ' '
      write ( *, '(a,i8)' ) '  Computed rank: ', rank

      return
      end 
      subroutine test02 ( )

c*********************************************************************72
c
cc TEST02 tests BAL_SEQ_TO_TABLEAU, TABLEAU_TO_BAL_SEQ.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    06 July 2007
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n
      parameter ( n = 4 )

      integer rank
      integer t(2*n)
      integer tab(2,n)

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST02'
      write ( *, '(a)' ) '  BAL_SEQ_TO_TABLEAU converts a balanced'
      write ( *, '(a)' ) '  sequence to a tableau;'
      write ( *, '(a)' ) '  TABLEAU_TO_BAL_SEQ converts a tableau'
      write ( *, '(a)' ) '  to a balanced sequence.'
c
c  Pick a "random" balanced sequence.
c
      rank = 7

      call bal_seq_unrank ( rank, n, t )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  "Random" balanced sequence:'
      write ( *, '(a)' ) ' '
      write ( *, '(4x,8i2)' ) t(1:2*n)
c
c  Convert to a tableau.
c
      call bal_seq_to_tableau ( n, t, tab )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Corresponding tableau'
      write ( *, '(a)' ) ' '
      write ( *, '(4x,4i2)' ) tab(1,1:n)
      write ( *, '(4x,4i2)' ) tab(2,1:n)
c
c  Convert to a balanced sequence.
c
      call tableau_to_bal_seq ( n, tab, t )

      call i4vec_print ( 2*n, t, '  Corresponding balanced sequence:' )

      return
      end 
      subroutine test03 ( )

c*********************************************************************72
c
cc TEST03 tests BELL_NUMBERS.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    25 July 2011
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n_max
      parameter ( n_max = 11 )

      integer b(0:n_max)
      integer bn
      integer n
      integer n_data

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST03'
      write ( *, '(a)' ) '  BELL_NUMBERS computes Bell numbers.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '         N        BELL(N)    BELL_NUMBERS(N)'
      write ( *, '(a)' ) ' '

      n_data = 0

10    continue

        call bell_values ( n_data, n, bn )

        if ( n_data .eq. 0 ) then
          go to 20
        end if

        call bell_numbers ( n, b )

        write ( *, '(2x,i8,2x,i12,2x,i12)' ) n, bn, b(n)

      go to 10

20    continue

      return
      end 
      subroutine test04 ( )

c*********************************************************************72
c
cc TEST04 tests BINOMIAL.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    06 July 2007
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer binomial
      integer i
      integer j

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST04'
      write ( *, '(a)' ) '  BINOMIAL computes binomial coefficients.'

      do i = -1, 5
        do j = -1, 5
          write ( *, '(3i8)' ) i, j, binomial ( i, j )
        end do
      end do

      return
      end 
      subroutine test05 ( )

c*********************************************************************72
c
cc TEST05 tests CYCLE_TO_PERM, PERM_TO_CYCLE.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    06 July 2007
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n
      parameter ( n = 7 )

      integer i
      integer jlo
      integer index(n)
      integer ncycle
      integer nperm
      integer p(n)
      integer rank
      integer t(n)

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST05'
      write ( *, '(a)' ) '  CYCLE_TO_PERM converts a permutation from'
      write ( *, '(a)' ) '  cycle to array form;'
      write ( *, '(a)' ) '  PERM_TO_CYCLE converts a permutation from'
      write ( *, '(a)' ) '  array to cycle form.'
c
c  Enumerate.
c
      call perm_enum ( n, nperm )
c
c  Choose a "random" permutation.
c
      rank = nperm / 2

      call perm_lex_unrank ( rank, n, p )

      call perm_print ( n, p, '  "Random" permutation:' )
c
c  Convert the permutation to cycle form.
c
      call perm_to_cycle ( n, p, ncycle, t, index )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Corresponding cycle form:'
      write ( *, '(a,i8)' ) '  Number of cycles is ', ncycle
      write ( *, '(a)' ) ' '
      jlo = 0
      do i = 1, ncycle
        write ( *, '(4x,20i4)' ) t(jlo+1:jlo+index(i))
        jlo = jlo + index(i)
      end do
c
c  Convert the set partition back to an RGF.
c
      call cycle_to_perm ( n, ncycle, t, index, p )

      call perm_print ( n, p, '  Corresponding permutation:' )

      return
      end 
      subroutine test06 ( )

c*********************************************************************72
c
cc TEST06 tests DIST_ENUM and DIST_NEXT.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    06 July 2007
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer k
      parameter ( k = 3 )

      integer idist
      integer m
      logical more
      integer num_dist
      integer q(k)

      m = 5
      more = .false.

      call dist_enum ( k, m, num_dist )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST06'
      write ( *, '(a)' ) '  For a distribution of M indistinguishable'
      write ( *, '(a)' ) '  objects among K distinguishable slots:'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  DIST_ENUM enumerates them;'
      write ( *, '(a)' ) '  DIST_NEXT produces the "next" one.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Number of:'
      write ( *, '(a,i8)' ) '    indistinguishable objects = ', m
      write ( *, '(a,i8)' ) '    distinguishable slots =     ', k
      write ( *, '(a,i8)' ) '    distributions is            ', num_dist
      write ( *, '(a)' ) ' '

      idist = 0

10    continue

        call dist_next ( k, m, q, more )

        if ( .not. more ) then
          go to 20
        end if

        idist = idist + 1
        write ( *, '(4x,6i5)' ) idist, q(1:k)

      go to 10

20    continue

      return
      end 
      subroutine test07 ( )

c*********************************************************************72
c
cc TEST07 tests I4_FACTORIAL.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    06 July 2007
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer i4_factorial
      integer fx
      integer fx2
      integer n
      integer x

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST07:'
      write ( *, '(a)' ) 
     &  '  I4_FACTORIAL evaluates the factorial function.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '     X       Exact F       FACTORIAL(X)'
      write ( *, '(a)' ) ' '

      n = 0

10    continue

        call i4_factorial_values ( n, x, fx )

        if ( n .eq. 0 ) then
          go to 20
        end if

        if ( x .le. 0.0D+00 ) then
          go to 10
        end if

        fx2 = i4_factorial ( x )

        write ( *, '(i4,2i12)' ) x, fx, fx2

      go to 10

20    continue

      return
      end 
      subroutine test08 ( )

c*********************************************************************72
c
cc TEST08 tests GRAY_CODE_*.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    06 July 2007
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n
      parameter ( n = 5 )

      integer ngray
      integer rank
      integer rank_old
      integer t(n)

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST08'
      write ( *, '(a)' ) '  Gray codes:'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  GRAY_CODE_ENUM enumerates,'
      write ( *, '(a)' ) '  GRAY_CODE_RANK ranks,'
      write ( *, '(a)' ) '  GRAY_CODE_SUCCESSOR lists,'
      write ( *, '(a)' ) '  GRAY_CODE_UNRANK unranks.'
c
c  Enumerate.
c
      call gray_code_enum ( n, ngray )

      write ( *, '(a)' ) ' '
      write ( *, '(a,i8)' ) '  For N = ', n
      write ( *, '(a,i8)' ) 
     &  '  the number of Gray code elements is ', ngray
      write ( *, '(a)' ) ' '
c
c  List
c
      rank = -1

10    continue

        rank_old = rank

        call gray_code_successor ( n, t, rank )

        if ( rank .le. rank_old ) then
          go to 20
        end if

        write ( *, '(4x,6i5)' ) rank, t(1:n)

      go to 10

20    continue
c
c  Unrank.
c
      rank = ngray / 2

      call gray_code_unrank ( rank, n, t )

      write ( *, '(a)' ) ' '
      write ( *, '(a,i8)' ) '  The element of rank ', rank
      write ( *, '(a)' ) ' '
      write ( *, '(4x,6i5)' ) t(1:n)
c
c  Rank.
c
      call gray_code_rank ( n, t, rank )

      call i4vec_print ( n, t, '  Element to be ranked:' )

      write ( *, '(a)' ) ' '
      write ( *, '(a,i8)' ) '  Computed rank: ', rank

      return
      end 
      subroutine test09 ( )

c*********************************************************************72
c
cc TEST09 tests I4VEC_SEARCH_BINARY_A and I4VEC_SORT_INSERT_A.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    06 July 2007
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n
      parameter ( n = 10 )

      integer a(n)
      integer b
      integer index

      save a

      data a / 6, 7, 1, 0, 4, 3, 2, 1, 5, 8 /

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST09'
      write ( *, '(a)' ) '  Integer vectors:'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  I4VEC_SORT_INSERT_A ascending sorts;'
      write ( *, '(a)' ) 
     &  '  I4VEC_SEARCH_BINARY_A searches a ascending sorted vector.'

      call i4vec_print ( n, a, '  Before ascending sort:' )

      call i4vec_sort_insert_a ( n, a )

      call i4vec_print ( n, a, '  After ascending sort:' )

      b = 5

      write ( *, '(a)' ) ' '
      write ( *, '(a,i8)' ) 
     &  '  Now search for an instance of the value ', b

      call i4vec_search_binary_a ( n, a, b, index )

      write ( *, '(a)' ) ' '
      if ( index .eq. 0 ) then
        write ( *, '(a)' ) '  The value does not occur.'
      else
        write ( *, '(a,i8)' ) '  The value occurs at index = ', index
      end if

      return
      end 
      subroutine test10 ( )

c*********************************************************************72
c
cc TEST10 tests I4VEC_SEARCH_BINARY_D and I4VEC_SORT_INSERT_D.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    06 July 2007
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n
      parameter ( n = 10 )

      integer a(n)
      integer b
      integer index

      save a

      data a / 6, 7, 1, 0, 4, 3, 2, 1, 5, 8 /

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST10'
      write ( *, '(a)' ) '  Integer vectors:'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  I4VEC_SORT_INSERT_D descending sorts;'
      write ( *, '(a)' ) '  I4VEC_SEARCH_BINARY_D searches a descending'
      write ( *, '(a)' ) '  sorted vector.'

      call i4vec_print ( n, a, '  Before descending sort:' )

      call i4vec_sort_insert_d ( n, a )

      call i4vec_print ( n, a, '  After descending sort:' )

      b = 5

      write ( *, '(a)' ) ' '
      write ( *, '(a,i8)' ) 
     &  '  Now search for an instance of the value ', b

      call i4vec_search_binary_d ( n, a, b, index )

      write ( *, '(a)' ) ' '
      if ( index .eq. 0 ) then
        write ( *, '(a)' ) '  The value does not occur.'
      else
        write ( *, '(a,i8)' ) '  The value occurs at index = ', index
      end if

      return
      end 
      subroutine test11 ( )

c*********************************************************************72
c
cc TEST11 tests KNAPSACK_REORDER and KNAPSACK_01.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    06 July 2007
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n
      parameter ( n = 5 )

      integer i
      double precision mass
      double precision mass_limit
      double precision p(n)
      double precision  profit
      double precision w(n)
      double precision x(n)

      save p
      save w

      data p / 24.0, 13.0, 23.0, 15.0, 16.0 /
      data w / 12.0,  7.0, 11.0,  8.0,  9.0 /

      mass_limit = 26.D+00

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST11'
      write ( *, '(a)' ) 
     &  '  KNAPSACK_REORDER reorders the knapsack data.'
      write ( *, '(a)' ) 
     &  '  KNAPSACK_01 solves the 0/1 knapsack problem.'

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Object, Profit, Mass, "Profit Density"'
      write ( *, '(a)' ) ' '
      do i = 1, n
        write ( *, '(2x,i6,2x,f7.3,2x,f7.3,2x,f7.3)' ) 
     &    i, p(i), w(i), p(i)/w(i)
      end do

      call knapsack_reorder ( n, p, w )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  After reordering by Profit Density:'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Object, Profit, Mass, "Profit Density"'
      write ( *, '(a)' ) ' '
      do i = 1, n
        write ( *, '(2x,i6,2x,f7.3,2x,f7.3,2x,f7.3)' ) 
     &    i, p(i), w(i), p(i) / w(i)
      end do

      write ( *, '(a)' ) ' '
      write ( *, '(a,f7.3)' ) '  Total mass restriction is ', mass_limit

      call knapsack_01 ( n, mass_limit, p, w, x, mass, profit )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Object, Density, Choice, Profit, Mass'
      write ( *, '(a)' ) ' '
      do i = 1, n
        write ( *, '(i6,f7.3,f7.3,2f7.3)' ) i, p(i)/w(i), x(i), 
     &    x(i) * p(i), x(i) * w(i)
      end do

      write ( *, '(a)' ) ' '
      write ( *, '(a,2f7.3)' ) '  Total:            ', profit, mass

      return
      end 
      subroutine test12 ( )

c*********************************************************************72
c
cc TEST12 tests KNAPSACK_REORDER and KNAPSACK_RATIONAL.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    06 July 2007
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n
      parameter ( n = 5 )

      integer i
      double precision mass
      double precision mass_limit
      double precision p(n)
      double precision  profit
      double precision w(n)
      double precision x(n)

      save p
      save w

      data p / 24.0, 13.0, 23.0, 15.0, 16.0 /
      data w / 12.0,  7.0, 11.0,  8.0,  9.0 /

      mass_limit = 26.0D+00

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST12'
      write ( *, '(a)' ) 
     &  '  KNAPSACK_REORDER reorders the knapsack data.'
      write ( *, '(a)' ) 
     &  '  KNAPSACK_RATIONAL solves the rational knapsack problem.'

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Object, Profit, Mass, "Profit Density"'
      write ( *, '(a)' ) ' '
      do i = 1, n
        write ( *, '(i6,3f7.3)' ) i, p(i), w(i), p(i) / w(i)
      end do

      call knapsack_reorder ( n, p, w )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  After reordering by Profit Density:'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Object, Profit, Mass, "Profit Density"'
      write ( *, '(a)' ) ' '
      do i = 1, n
        write ( *, '(i6,3f7.3)' ) i, p(i), w(i), p(i) / w(i)
      end do

      write ( *, '(a)' ) ' '
      write ( *, '(a,f7.3)' ) 
     &  '  Total mass restriction is ', mass_limit

      call knapsack_rational ( n, mass_limit, p, w, x, mass, profit )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Object, Density, Choice, Profit, Mass'
      write ( *, '(a)' ) ' '
      do i = 1, n
        write ( *, '(i6,f7.3,f7.3,2f7.3)' ) i, p(i) / w(i), x(i), 
     &    x(i) * p(i), x(i) * w(i)
      end do

      write ( *, '(a)' ) ' '
      write ( *, '(a,2f7.3)' ) '  Total:            ', profit, mass

      return
      end 
      subroutine test13 ( )

c*********************************************************************72
c
cc TEST13 tests KSUBSET_COLEX_RANK, _SUCCESSOR, _UNRANK, _ENUM.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    06 July 2007
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer k
      parameter ( k = 3 )

      integer n
      integer nksub
      integer rank
      integer rank_old
      integer t(k)

      n = 5

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST13'
      write ( *, '(a)' ) '  K-subsets of an N set,'
      write ( *, '(a)' ) '  using the colexicographic ordering:'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  KSUBSET_COLEX_RANK ranks,'
      write ( *, '(a)' ) '  KSUBSET_COLEX_SUCCESSOR lists,'
      write ( *, '(a)' ) '  KSUBSET_COLEX_UNRANK unranks.'
      write ( *, '(a)' ) '  KSUBSET_ENUM enumerates,'
c
c  Enumerate.
c
      call ksubset_enum ( k, n, nksub )

      write ( *, '(a)' ) ' '
      write ( *, '(a,i8)' ) '  For N = ', n
      write ( *, '(a,i8)' ) '  the number of K subsets is ', nksub
      write ( *, '(a)' ) ' '
c
c  List
c
      rank = -1

10    continue

        rank_old = rank

        call ksubset_colex_successor ( k, n, t, rank )

        if ( rank .le. rank_old ) then
          go to 20
        end if

        write ( *, '(4x,6i5)' ) rank, t(1:k)

      go to 10

20    continue
c
c  Unrank.
c
      rank = nksub / 2

      call ksubset_colex_unrank ( rank, k, n, t )

      write ( *, '(a)' ) ' '
      write ( *, '(a,i8)' ) '  The element of rank ', rank
      write ( *, '(a)' ) ' '
      write ( *, '(4x,6i5)' ) t(1:k)
c
c  Rank.
c
      call ksubset_colex_rank ( k, n, t, rank )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  The rank of the element:'
      write ( *, '(a)' ) ' '
      write ( *, '(4x,6i5)' ) t(1:k)
      write ( *, '(a)' ) ' '
      write ( *, '(a,i8)' ) '  is computed as ', rank

      return
      end 
      subroutine test14 ( )

c*********************************************************************72
c
cc TEST14 tests KSUBSET_ENUM, _LEX_RANK, _LEX_SUCCESSOR, _LEX_UNRANK.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    06 July 2007
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer k
      parameter ( k = 3 )

      integer n
      integer nksub
      integer rank
      integer rank_old
      integer t(k)

      n = 5

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST14'
      write ( *, '(a)' ) '  K-subsets of an N set,'
      write ( *, '(a)' ) '  using the lexicographic ordering:'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  KSUBSET_ENUM enumerates,'
      write ( *, '(a)' ) '  KSUBSET_LEX_RANK ranks,'
      write ( *, '(a)' ) '  KSUBSET_LEX_SUCCESSOR lists,'
      write ( *, '(a)' ) '  KSUBSET_LEX_UNRANK unranks.'
c
c  Enumerate.
c
      call ksubset_enum ( k, n, nksub )

      write ( *, '(a)' ) ' '
      write ( *, '(a,i8)' ) '  For N = ', n
      write ( *, '(a,i8)' ) '  the number of K subsets is ', nksub
      write ( *, '(a)' ) ' '
c
c  List
c
      rank = -1

10    continue

        rank_old = rank

        call ksubset_lex_successor ( k, n, t, rank )

        if ( rank .le. rank_old ) then
          go to 20
        end if

        write ( *, '(6i5)' ) rank, t(1:k)

      go to 10

20    continue
c
c  Unrank.
c
      rank = nksub / 2

      call ksubset_lex_unrank ( rank, k, n, t )

      write ( *, '(a)' ) ' '
      write ( *, '(a,i8)' ) '  The element of rank ', rank
      write ( *, '(a)' ) ' '
      write ( *, '(6i5)' ) t(1:k)
c
c  Rank.
c
      call ksubset_lex_rank ( k, n, t, rank )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  The rank of the element:'
      write ( *, '(a)' ) ' '
      write ( *, '(6i5)' ) t(1:k)
      write ( *, '(a)' ) ' '
      write ( *, '(a,i8)' ) '  is computed as ', rank

      return
      end 
      subroutine test15 ( )

c*********************************************************************72
c
cc TEST15 tests KSUBSET_ENUM, _REVDOOR_RANK, _REVDOOR_SUCCESSOR, _REVDOOR_UNRANK.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    06 July 2007
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer k
      parameter ( k = 3 )

      integer n
      integer nksub
      integer rank
      integer rank_old
      integer t(k)

      n = 5

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST15'
      write ( *, '(a)' ) '  K-subsets of an N set,'
      write ( *, '(a)' ) '  using the revolving door ordering:'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  KSUBSET_ENUM enumerates,'
      write ( *, '(a)' ) '  KSUBSET_REVDOOR_RANK ranks,'
      write ( *, '(a)' ) '  KSUBSET_REVDOOR_SUCCESSOR lists,'
      write ( *, '(a)' ) '  KSUBSET_REVDOOR_UNRANK unranks.'
c
c  Enumerate.
c
      call ksubset_enum ( k, n, nksub )

      write ( *, '(a)' ) ' '
      write ( *, '(a,i8)' ) '  For N = ', n
      write ( *, '(a,i8)' ) '  the number of K subsets is ', nksub
      write ( *, '(a)' ) ' '
c
c  List
c
      rank = -1

10    continue

        rank_old = rank

        call ksubset_revdoor_successor ( k, n, t, rank )

        if ( rank .le. rank_old ) then
          go to 20
        end if

        write ( *, '(6i5)' ) rank, t(1:k)

      go to 10

20    continue
c
c  Unrank.
c
      rank = nksub / 2

      call ksubset_revdoor_unrank ( rank, k, n, t )

      write ( *, '(a)' ) ' '
      write ( *, '(a,i8)' ) '  The element of rank ', rank
      write ( *, '(a)' ) ' '
      write ( *, '(6i5)' ) t(1:k)
c
c  Rank.
c
      call ksubset_revdoor_rank ( k, n, t, rank )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  The rank of the element:'
      write ( *, '(a)' ) ' '
      write ( *, '(6i5)' ) t(1:k)
      write ( *, '(a)' ) ' '
      write ( *, '(a,i8)' ) '  is computed as ', rank

      return
      end 
      subroutine test16 ( )

c*********************************************************************72
c
cc TEST16 tests MARRIAGE.
c
c  Discussion:
c
c    PREFER(M,W) is the index of women W on man M's list.
c    RANK(W,M) is the index of man M on woman W's list.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    06 July 2007
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n
      parameter ( n = 5 )

      integer fiancee(n)
      integer i
      integer next(n)
      integer prefer(n,n)
      integer rank(n,n)

      save prefer
      save rank

      data prefer /
     &  2, 1, 2, 1, 5, 
     &  5, 2, 3, 3, 3, 
     &  1, 3, 5, 2, 2, 
     &  3, 4, 4, 4, 1, 
     &  4, 5, 1, 5, 4  /

      data rank /
     &  2, 4, 1, 4, 5, 
     &  4, 3, 3, 2, 2, 
     &  5, 5, 4, 1, 3, 
     &  3, 1, 2, 3, 1, 
     &  1, 2, 5, 5, 4 /

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST16'
      write ( *, '(a)' ) '  MARRIAGE arranges a set of stable marriages'
      write ( *, '(a)' ) '  given a set of preferences.'

      call marriage ( n, prefer, rank, fiancee, next )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Man, Wife''s rank, Wife'
      write ( *, '(a)' ) ' '
      do i = 1, n
        write ( *, '(3i8)' ) i, next(i), prefer(i,next(i))
      end do

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Woman, Husband''s rank, Husband'
      write ( *, '(a)' ) ' '
      do i = 1, n
        write ( *, '(3i8)' ) i, rank(i,fiancee(i)), fiancee(i)
      end do

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Correct result:'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  M:W 1  2  3  4  5'
      write ( *, '(a)' ) '   1  +  .  .  .  .'
      write ( *, '(a)' ) '   2  .  .  .  +  .'
      write ( *, '(a)' ) '   3  .  .  .  .  +'
      write ( *, '(a)' ) '   4  .  .  +  .  .'
      write ( *, '(a)' ) '   5  .  +  .  .  .'

      return
      end 
      subroutine test17 ( )

c*********************************************************************72
c
cc TEST17 tests MOUNTAIN.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    06 July 2007
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n
      parameter ( n = 5 )

      integer mxy
      integer row(0:2*n)
      integer x
      integer y

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST17'
      write ( *, '(a)' ) '  MOUNTAIN computes mountain numbers.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '   Y    MXY'
      write ( *, '(a)' ) ' '

      do y = 0, n
        do x = 0, 2*n
          call mountain ( n, x, y, mxy )
          row(x) = mxy
        end do
        write ( *, '(2x,i2,3x,11i4)' ) y, row(0:2*n)
      end do

      return
      end 
      subroutine test18 ( )

c*********************************************************************72
c
cc TEST18 tests NPART_ENUM, _RSF_LEX_RANK, _RSF_LEX_SUCCESSOR, _RSF_LEX_UNRANK.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    06 July 2007
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer npart
      parameter ( npart = 3 )

      integer n
      integer npartitions
      integer rank
      integer rank_old
      integer t(npart)

      n = 12

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST18'
      write ( *, '(a)' ) '  Partitions of N with NPART parts'
      write ( *, '(a)' ) '  in reverse standard form:'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  NPART_ENUM enumerates,'
      write ( *, '(a)' ) '  NPART_RSF_LEX_RANK ranks,'
      write ( *, '(a)' ) '  NPART_RSF_LEX_SUCCESSOR lists;'
      write ( *, '(a)' ) '  NPART_RSF_LEX_UNRANK unranks.'
c
c  Enumerate.
c
      call npart_enum ( n, npart, npartitions )

      write ( *, '(a)' ) ' '
      write ( *, '(a,i8)' ) '  For N = ', n
      write ( *, '(a,i8)' ) '  and NPART = ', npart
      write ( *, '(a,i8)' ) 
     &  '  the number of partitions is ', npartitions
      write ( *, '(a)' ) ' '
c
c  List.
c
      rank = -1

10    continue

        rank_old = rank

        call npart_rsf_lex_successor ( n, npart, t, rank )

        if ( rank .le. rank_old ) then
          go to 20
        end if

        write ( *, '(6i5)' ) rank, t(1:npart)

      go to 10

20    continue
c
c  Unrank.
c
      rank = npartitions / 3

      call npart_rsf_lex_unrank ( rank, n, npart, t )

      write ( *, '(a)' ) ' '
      write ( *, '(a,i8)' ) '  The element of rank ', rank
      write ( *, '(a)' ) ' '
      write ( *, '(6i5)' ) t(1:npart)
c
c  Rank.
c
      call npart_rsf_lex_rank ( n, npart, t, rank )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  The rank of the element:'
      write ( *, '(a)' ) ' '
      write ( *, '(6i5)' ) t(1:npart)
      write ( *, '(a)' ) ' '
      write ( *, '(a,i8)' ) '  is computed as ', rank

      return
      end 
      subroutine test19 ( )

c*********************************************************************72
c
cc TEST19 tests NPART_RSF_LEX_RANDOM;
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    06 July 2007
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer npart
      parameter ( npart = 3 )

      integer i
      integer n
      integer seed
      integer t(npart)

      n = 12
      seed = 123456789

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST19'
      write ( *, '(a)' ) '  Partitions of N with NPART parts'
      write ( *, '(a)' ) '  in reverse standard form:'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 
     &  '  NPART_RSF_LEX_RANDOM produces random examples.'

      do i = 1, 10

        call npart_rsf_lex_random ( n, npart, seed, t )

        write ( *, '(6i5)' ) t(1:npart)

      end do

      return
      end 
      subroutine test20 ( )

c*********************************************************************72
c
cc TEST20 tests NPART_ENUM and NPART_SF_SUCCESSOR;
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    06 July 2007
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer npart
      parameter ( npart = 3 )

      integer n
      integer npartitions
      integer rank
      integer rank_old
      integer t(npart)

      n = 12

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST20'
      write ( *, '(a)' ) '  Partitions of N with NPART parts'
      write ( *, '(a)' ) '  in standard form:'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  NPART_ENUM enumerates,'
      write ( *, '(a)' ) '  NPART_SF_LEX_SUCCESSOR lists.'
c
c  Enumerate.
c
      call npart_enum ( n, npart, npartitions )

      write ( *, '(a)' ) ' '
      write ( *, '(a,i8)' ) '  For N = ', n
      write ( *, '(a,i8)' ) '  and NPART = ', npart
      write ( *, '(a,i8)' ) 
     &  '  the number of partitions is ', npartitions
      write ( *, '(a)' ) ' '
c
c  List.
c
      rank = -1

10    continue

        rank_old = rank

        call npart_sf_lex_successor ( n, npart, t, rank )

        if ( rank .le. rank_old ) then
          go to 20
        end if

        write ( *, '(6i5)' ) rank, t(1:npart)

      go to 10

20    continue

      return
      end 
      subroutine test21 ( )

c*********************************************************************72
c
cc TEST21 tests NPART_TABLE and PART_TABLE.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    06 July 2007
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer maxn
      parameter ( maxn = 10 )
      integer maxpart
      parameter ( maxpart = 5 )

      integer i
      integer p(0:maxn,0:maxpart)
      integer p2(0:maxn)

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST21'
      write ( *, '(a)' ) '  NPART_TABLE tabulates partitions'
      write ( *, '(a)' ) '  of N with NPART parts;'
      write ( *, '(a)' ) '  PART_TABLE tabulates partitions of N.'

      call npart_table ( maxn, maxpart, maxn, p )

      call part_table ( maxn, p2 )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 
     &  '    I P(I)  P(I,0) P(I,1) P(I,2) P(I,3) P(I,4) P(I,5)'
      write ( *, '(a)' ) ' '

      do i = 0, maxn
        write ( *, '(11i5)' ) i, p2(i), p(i,0:maxpart)
      end do

      return
      end 
      subroutine test22 ( )

c*********************************************************************72
c
cc TEST22 tests PART_ENUM and PART_SUCCESSOR.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    06 July 2007
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n
      parameter ( n = 8 )

      integer npart
      integer npartitions
      integer rank
      integer rank_old
      integer t(n)

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST22'
      write ( *, '(a)' ) '  PART_SUCCESSOR produces partitions of N,'
      write ( *, '(a)' ) '  PART_ENUM enumerates.'
      write ( *, '(a)' ) ' '
      write ( *, '(a,i8)' ) '  Partitions of N = ', n
c
c  Enumerate.
c
      call part_enum ( n, npartitions )

      write ( *, '(a)' ) ' '
      write ( *, '(a,i8)' ) '  For N = ', n
      write ( *, '(a,i8)' ) 
     &  '  the number of partitions is ', npartitions
      write ( *, '(a)' ) ' '
c
c  List.
c
      rank = -1

10    continue

        rank_old = rank

        call part_successor ( n, npart, t, rank )

        if ( rank .le. rank_old ) then
          go to 20
        end if

        write ( *, '(2x,i2,3x,10i3)' ) rank, t(1:npart)

      go to 10

20    continue

      return
      end 
      subroutine test23 ( )

c*********************************************************************72
c
cc TEST23 tests PART_SUCCESSOR and PART_SF_CONJUGATE.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    06 July 2007
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n
      parameter ( n = 8 )

      integer b(n)
      integer npart
      integer npartb
      integer rank
      integer rank_old
      integer t(n)

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST23'
      write ( *, '(a)' ) '  PART_SUCCESSOR produces partitions of N,'
      write ( *, '(a)' ) 
     &  '  PART_SF_CONJUGATE produces the conjugate of a partition.'
      write ( *, '(a)' ) ' '
      write ( *, '(a,i8)' ) '  Partitions of N = ', n
c
c  List.
c
      rank = -1

10    continue

        rank_old = rank

        call part_successor ( n, npart, t, rank )

        if ( rank .le. rank_old ) then
          go to 20
        end if

        write ( *, '(2x,i2,4x,10i3)' ) rank, t(1:npart)
        call part_sf_conjugate ( n, npart, t, npartb, b )
        write ( *, '(2x,a4,2x,10i3)' ) 'Con:', b(1:npartb)

      go to 10

20    continue

      return
      end 
      subroutine test24 ( )

c*********************************************************************72
c
cc TEST24 tests PART_SF_MAJORIZE.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    06 July 2007
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n
      parameter ( n = 8 )

      integer a(n)
      integer b(n)
      integer c(n)
      integer i
      integer nparta
      parameter ( nparta = 5 )
      integer npartb
      parameter ( npartb = 6 )
      integer npartc
      parameter ( npartc = 6 )
      integer result

      save a
      save b
      save c

      data a / 2, 2, 2, 1, 1, 0, 0, 0 /
      data b / 3, 1, 1, 1, 1, 1, 0, 0 /
      data c / 2, 2, 1, 1, 1, 1, 0, 0 /

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST24'
      write ( *, '(a)' ) 
     &  '  PART_SF_MAJORIZE determines if one partition'
      write ( *, '(a)' ) '  majorizes another.'
      write ( *, '(a)' ) ' '
      write ( *, '(a,i8)' ) '  Partitions of N = ', n
      write ( *, '(a)' ) ' '
      write ( *, '(2x,a,2x,10i3)' ) 'A:', ( a(i), i = 1, nparta )
      write ( *, '(2x,a,2x,10i3)' ) 'B:', ( b(i), i = 1, npartb )
      write ( *, '(2x,a,2x,10i3)' ) 'C:', ( c(i), i = 1, npartc )

      call part_sf_majorize ( n, nparta, a, npartb, b, result )
      write ( *, '(a)' ) ' '
      write ( *, '(a,i8)' ) '  A compare B: ', result
      call part_sf_majorize ( n, npartb, b, npartc, c, result )
      write ( *, '(a,i8)' ) '  B compare C: ', result
      call part_sf_majorize ( n, npartc, c, nparta, a, result )
      write ( *, '(a,i8)' ) '  C compare A: ', result
      call part_sf_majorize ( n, npartc, c, npartc, c, result )
      write ( *, '(a,i8)' ) '  C compare C: ', result

      return
      end 
      subroutine test25 ( )

c*********************************************************************72
c
cc TEST25 tests PARTITION_GREEDY.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    06 July 2007
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n
      parameter ( n = 10 )

      integer a(n)
      integer i
      integer indx(n)
      integer sums(2)

      save a

      data a / 2, 10, 3, 8, 5, 7, 9, 5, 3, 2 /
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST25'
      write ( *, '(a)' ) 
     &  '  PARTITION_GREEDY partitions an integer vector into'
      write ( *, '(a)' ) '  two subsets with nearly equal sum.'

      call partition_greedy ( n, a, indx )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'Data set #1 partitioned:'
      write ( *, '(a)' ) ' '

      sums(1) = 0
      sums(2) = 0

      do i = 1, n

        if ( indx(i) .eq. 1 ) then
          sums(1) = sums(1) + a(i)
          write ( *, '(2x,i6)' ) a(i)
        else
          write ( *, '(2x,6x,i6)' ) a(i)
          sums(2) = sums(2) + a(i)
        end if

      end do

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Sums:'
      write ( *, '(a)' ) ' '
      write ( *, '(2i6)' ) sums(1), sums(2)

      a = (/ 771, 121, 281, 854, 885, 734, 486, 1003, 83, 62 /)

      call partition_greedy ( n, a, indx )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Data set #2 partitioned:'
      write ( *, '(a)' ) ' '

      sums(1:2) = 0

      do i = 1, n

        if ( indx(i) .eq. 1 ) then
          sums(1) = sums(1) + a(i)
          write ( *, '(i6)' ) a(i)
        else
          write ( *, '(6x,i6)' ) a(i)
          sums(2) = sums(2) + a(i)
        end if

      end do

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Sums:'
      write ( *, '(a)' ) ' '
      write ( *, '(2i6)' ) sums(1), sums(2)

      return
      end 
      subroutine test26 ( )

c*********************************************************************72
c
cc TEST26 tests PARTN_ENUM, PARTN_SUCCESSOR and PART_SF_CONJUGATE.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    06 July 2007
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n
      parameter ( n = 11 )

      integer b(n)
      integer nmax
      integer npart
      integer npart2
      integer npartitions
      integer rank
      integer rank_old
      integer t(n)

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST26'
      write ( *, '(a)' ) '  Partitions of N with maximum element NMAX:'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  PARTN_SUCCESSOR lists;'
      write ( *, '(a)' ) '  PARTN_ENUM enumerates.'

      nmax = 4
c
c  Enumerate.
c
      call partn_enum ( n, nmax, npartitions )

      write ( *, '(a)' ) ' '
      write ( *, '(a,i8)' ) '  For N = ', n
      write ( *, '(a,i8)' ) '  and NMAX = ', nmax
      write ( *, '(a,i8)' ) 
     &  '  the number of partitions is ', npartitions
      write ( *, '(a)' ) ' '
c
c  List.
c
      rank = -1

10    continue

        rank_old = rank

        call partn_successor ( n, nmax, npart, t, rank )

        if ( rank .le. rank_old ) then
          go to 20
        end if

        write ( *, '(2x,i2,3x,15i3)' ) rank, t(1:npart)

      go to 10

20    continue
c
c  List conjugates.
c
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Repeat, but list RSF conjugated partitions.'
      write ( *, '(a)' ) ' '
      rank = -1

30    continue

        rank_old = rank

        call partn_successor ( n, nmax, npart, t, rank )

        if ( rank .le. rank_old ) then
          go to 40
        end if

        call part_sf_conjugate ( n, npart, t, npart2, b )
        call i4vec_reverse ( npart2, b )
        write ( *, '(2x,i2,3x,15i3)' ) rank, b(1:npart2)

      go to 30

40    continue

      return
      end 
      subroutine test27 ( )

c*********************************************************************72
c
cc TEST27 tests PERM_INV and PERM_MUL.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    06 July 2007
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n
      parameter ( n = 4 )

      integer nperm
      integer p(n)
      integer q(n)
      integer r(n)
      integer rank

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST27'
      write ( *, '(a)' ) '  Permutations of the integers:'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  PERM_INV computes an inverse permutation,'
      write ( *, '(a)' ) '  PERM_MUL multiplies two permutations.'
c
c  Enumerate.
c
      call perm_enum ( n, nperm )
c
c  Unrank.
c
      rank = nperm / 2

      call perm_lex_unrank ( rank, n, p )

      call perm_print ( n, p, '  The permutation P is ' )
c
c  Invert.
c
      call perm_inv ( n, p, q )

      call perm_print ( n, q, '  The inverse permutation Q is ' )
c
c  Multiply.
c
      call perm_mul ( n, p, q, r )

      call perm_print ( n, r, '  The product R = P * Q is ' )

      return
      end 
      subroutine test28 ( )

c*********************************************************************72
c
cc TEST28 tests PERM_ENUM, _LEX_RANK, _LEX_SUCCESSOR, _LEX_UNRANK.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    06 July 2007
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n
      parameter ( n = 4 )

      integer nperm
      integer pi(n)
      integer rank
      integer rank_old

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST28'
      write ( *, '(a)' ) '  Permutations of the integers,'
      write ( *, '(a)' ) '  using the lexicographic ordering:'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  PERM_ENUM enumerates,'
      write ( *, '(a)' ) '  PERM_LEX_RANK ranks,'
      write ( *, '(a)' ) '  PERM_LEX_SUCCESSOR lists,'
      write ( *, '(a)' ) '  PERM_LEX_UNRANK unranks.'
c
c  Enumerate.
c
      call perm_enum ( n, nperm )

      write ( *, '(a)' ) ' '
      write ( *, '(a,i8)' ) '  For N = ', n
      write ( *, '(a,i8)' ) '  the number of permutations is ', nperm
      write ( *, '(a)' ) ' '
c
c  List
c
      rank = -1

10    continue

        rank_old = rank

        call perm_lex_successor ( n, pi, rank )

        if ( rank .le. rank_old ) then
          go to 20
        end if

        write ( *, '(6i5)' ) rank, pi(1:n)

      go to 10

20    continue
c
c  Unrank.
c
      rank = nperm / 2

      call perm_lex_unrank ( rank, n, pi )

      write ( *, '(a)' ) ' '
      write ( *, '(a,i8)' ) '  The element of rank ', rank

      call perm_print ( n, pi, ' ' )
c
c  Rank.
c
      call perm_lex_rank ( n, pi, rank )

      call perm_print ( n, pi, '  The rank of the element:' )
      write ( *, '(a)' ) ' '
      write ( *, '(a,i8)' ) 'is computed as ', rank

      return
      end 
      subroutine test29 ( )

c*********************************************************************72
c
cc TEST29 tests PERM_TJ_ENUM, _TJ_RANK, _TJ_SUCCESSOR, _TJ_UNRANK.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    06 July 2007
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n
      parameter ( n = 4 )

      integer nperm
      integer pi(n)
      integer rank
      integer rank_old

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST29'
      write ( *, '(a)' ) '  Permutations of the integers'
      write ( *, '(a)' ) '  using the Trotter-Johnson ordering:'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  PERM_ENUM enumerates,'
      write ( *, '(a)' ) '  PERM_TJ_RANK ranks,'
      write ( *, '(a)' ) '  PERM_TJ_SUCCESSOR lists,'
      write ( *, '(a)' ) '  PERM_TJ_UNRANK unranks.'
c
c  Enumerate.
c
      call perm_enum ( n, nperm )

      write ( *, '(a)' ) ' '
      write ( *, '(a,i8)' ) '  For N = ', n
      write ( *, '(a,i8)' ) '  the number of permutations is ', nperm
      write ( *, '(a)' ) ' '
c
c  List
c
      rank = -1

10    continue

        rank_old = rank

        call perm_tj_successor ( n, pi, rank )

        if ( rank .le. rank_old ) then
          go to 20
        end if

        write ( *, '(6i5)' ) rank, pi(1:n)

      go to 10

20    continue
c
c  Unrank.
c
      rank = nperm / 2

      call perm_tj_unrank ( rank, n, pi )

      write ( *, '(a)' ) ' '
      write ( *, '(a,i4)' ) '  The element of rank ', rank

      call perm_print ( n, pi, ' ' )
c
c  Rank.
c
      call perm_tj_rank ( n, pi, rank )

      call perm_print ( n, pi, '  The rank of the element:' )

      write ( *, '(a)' ) ' '
      write ( *, '(a,i8)' ) '  is computed as ', rank

      return
      end 
      subroutine test30 ( )

c*********************************************************************72
c
cc TEST30 tests PRUEFER_ENUM, PRUEFER_RANK, PRUEFER_SUCCESSOR, PRUEFER_UNRANK.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    06 July 2007
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n
      parameter ( n = 4 )

      integer ncode
      integer p(n-2)
      integer rank
      integer rank_old

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST30'
      write ( *, '(a)' ) '  Pruefer codes:'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  PRUEFER_ENUM enumerates,'
      write ( *, '(a)' ) '  PRUEFER_RANK ranks,'
      write ( *, '(a)' ) '  PRUEFER_SUCCESSOR lists,'
      write ( *, '(a)' ) '  PRUEFER_UNRANK unranks.'
c
c  Enumerate.
c
      call pruefer_enum ( n, ncode )

      write ( *, '(a)' ) ' '
      write ( *, '(a,i8)' ) '  For N = ', n
      write ( *, '(a,i8)' ) '  the number of Pruefer codes is ', ncode
      write ( *, '(a)' ) ' '
c
c  List
c
      rank = -1

10    continue

        rank_old = rank

        call pruefer_successor ( n, p, rank )

        if ( rank .le. rank_old ) then
          go to 20
        end if

        write ( *, '(6i5)' ) rank, p(1:n-2)

      go to 10

20    continue
c
c  Unrank.
c
      rank = ncode / 2

      call pruefer_unrank ( rank, n, p )

      write ( *, '(a)' ) ' '
      write ( *, '(a,i8)' ) '  The element of rank ', rank
      write ( *, '(a)' ) ' '
      write ( *, '(6i5)' ) p(1:n-2)
c
c  Rank.
c
      call pruefer_rank ( n, p, rank )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  The rank of the element:'
      write ( *, '(a)' ) ' '
      write ( *, '(6i5)' ) p(1:n-2)
      write ( *, '(a)' ) ' '
      write ( *, '(a,i8)' ) '  is computed as ', rank

      return
      end 
      subroutine test31 ( )

c*********************************************************************72
c
cc TEST31 tests PRUEFER_TO_TREE and TREE_TO_PRUEFER.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    06 July 2007
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n
      parameter ( n = 5 )

      integer i4_hi
      integer i4_lo
      integer i4_uniform
      integer j
      integer p(n-2)
      integer pruefer_num
      integer rank
      integer :: seed = 123456789
      integer t(2,n-1)
      integer test
      integer, parameter :: test_num = 5

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST31'
      write ( *, '(a)' ) '  PRUEFER_TO_TREE converts a Pruefer code to'
      write ( *, '(a)' ) '  a tree;'
      write ( *, '(a)' ) 
     &  '  TREE_TO_PRUEFER converts a tree to a Pruefer'
      write ( *, '(a)' ) '  code.'

      call pruefer_enum ( n, pruefer_num )

      i4_lo = 0
      i4_hi = pruefer_num - 1

      do test = 1, test_num
c
c  Pick a "random" Pruefer code.
c
        rank = i4_uniform ( i4_lo, i4_hi, seed )

        call pruefer_unrank ( rank, n, p )

        write ( *, '(a)' ) ' '
        write ( *, '(a,i8)' ) '  Random Pruefer code of rank ', rank
        write ( *, '(6i5)' ) p(1:n-2)
c
c  Convert the Pruefer code to a tree.
c
        call pruefer_to_tree ( n, p, t )

        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) '  Edge list for the corresponding tree:'
        write ( *, '(a)' ) ' '
        do j = 1, n - 1
          write ( *, '(6i5)' ) j, t(1:2,j)
        end do
c
c  Convert the tree to a Pruefer code.
c
        call tree_to_pruefer ( n, t, p )

        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) '  Corresponding Pruefer code:'
        write ( *, '(6i5)' ) p(1:n-2)

      end do

      return
      end 
      subroutine test32 ( )

c*********************************************************************72
c
cc TEST32 tests QUEENS and BACKTRACK.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    06 July 2007
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n
      parameter ( n = 8 )

      integer maxstack
      parameter ( maxstack = n * n )

      integer iarray(n)
      integer indx
      integer istack(maxstack)
      integer k
      integer nstack

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST32'
      write ( *, '(a)' ) '  QUEENS produces nonattacking queens'
      write ( *, '(a)' ) '  on a chessboard.'
      write ( *, '(a)' ) '  BACKTRACK supervises a backtrack search.'
      write ( *, '(a)' ) ' '

      indx = 0

10    continue

        call backtrack ( n, iarray, indx, k, nstack, istack, maxstack )

        if ( indx .eq. 1 ) then

          write ( *, '(19i4)' ) iarray(1:n)

        else if ( indx .eq. 2 ) then

          call queens ( n, iarray, k, nstack, istack, maxstack )

        else

          go to 20

        end if

      go to 10

20    continue

      return
      end 
      subroutine test33 ( )

c*********************************************************************72
c
cc TEST33 tests RGF_G_TABLE.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    25 July 2011
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer mmax
      parameter ( mmax = 8 )

      integer d(0:mmax,0:mmax)
      integer i
      integer m

      m = 6

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST33'
      write ( *, '(a)' ) 
     &  '  RGF_G_TABLE tabulates generalized restricted'
      write ( *, '(a)' ) '  growth functions.'
      write ( *, '(a)' ) ' '

      call rgf_g_table ( m, mmax, d )

      do i = 0, m
        write ( *, '(7i6)' ) d(i,0:m-i)
      end do

      return
      end 
      subroutine test34 ( )

c*********************************************************************72
c
cc TEST34 tests RGF_ENUM, RGF_RANK, RGF_SUCCESSOR, RGF_UNRANK.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    06 July 2007
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer m
      parameter ( m = 4 )

      integer f(m)
      integer nrgf
      integer rank
      integer rank_old

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST34'
      write ( *, '(a)' ) '  Restricted growth functions:'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  RGF_ENUM enumerates,'
      write ( *, '(a)' ) '  RGF_RANK ranks,'
      write ( *, '(a)' ) '  RGF_SUCCESSOR lists;'
      write ( *, '(a)' ) '  RGF_UNRANK unranks.'
c
c  Enumerate.
c
      call rgf_enum ( m, nrgf )

      write ( *, '(a)' ) ' '
      write ( *, '(a,i8)' ) '  For M = ', m
      write ( *, '(a,i8)' ) '  the number of RGF''s is ', nrgf
      write ( *, '(a)' ) ' '
c
c  List.
c
      rank = -1

10    continue

        rank_old = rank

        call rgf_successor ( m, f, rank )

        if ( rank .le. rank_old ) then
          go to 20
        end if

        write ( *, '(6i5)' ) rank, f(1:m)

      go to 10

20    continue
c
c  Unrank.
c
      rank = nrgf / 2

      call rgf_unrank ( rank, m, f )

      write ( *, '(a)' ) ' '
      write ( *, '(a,i8)' ) '  The element of rank ', rank
      write ( *, '(a)' ) ' '
      write ( *, '(6i5)' ) f(1:m)
c
c  Rank.
c
      call rgf_rank ( m, f, rank )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  The rank of the element:'
      write ( *, '(a)' ) ' '
      write ( *, '(6i5)' ) f(1:m)
      write ( *, '(a)' ) ' '
      write ( *, '(a,i8)' ) '  is computed as ', rank
      return
      end 
      subroutine test35 ( )

c*********************************************************************72
c
cc TEST35 tests RGF_TO_SETPART and SETPART_TO_RGF.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    06 July 2007
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer m
      parameter ( m = 8 )

      integer i
      integer jlo
      integer f(m)
      integer index(m)
      integer nsub
      integer rank
      integer s(m)

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST35'
      write ( *, '(a)' ) '  RGF_TO_SETPART converts a balanced'
      write ( *, '(a)' ) '  sequence to a restricted growth function;'
      write ( *, '(a)' ) '  SETPART_TO_RGF converts a restricted growth'
      write ( *, '(a)' ) '  function to a balanced sequence.'
c
c  Choose a "random" RGF.
c
      rank = 7
      call rgf_unrank ( rank, m, f )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  "Random" restricted growth function:'
      write ( *, '(a)' ) ' '
      write ( *, '(8i2)' ) f(1:m)
c
c  Convert the RGF to a set partition.
c
      call rgf_to_setpart ( m, f, nsub, s, index )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Corresponding set partition'
      write ( *, '(a)' ) ' '
      jlo = 1
      do i = 1, nsub
        write ( *, '(8i4)' ) s(jlo:index(i))
        jlo = index(i) + 1
      end do
c
c  Convert the set partition back to an RGF.
c
      call setpart_to_rgf ( m, nsub, s, index, f )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'Corresponding restricted growth function:'
      write ( *, '(a)' ) ' '
      write ( *, '(8i2)' ) f(1:m)

      return
      end 
      subroutine test36 ( )

c*********************************************************************72
c
cc TEST36 tests SETPART_ENUM.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    06 July 2007
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n
      integer npart

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST36'
      write ( *, '(a)' ) '  Set partitions:'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  SETPART_ENUM enumerates,'
      write ( *, '(a)' ) ' '
c
c  Enumerate.
c
      do n = 1, 6
        call setpart_enum ( n, npart )
        write ( *, '(i6,i6)' ) n, npart
      end do

      return
      end 
      subroutine test37 ( )

c*********************************************************************72
c
cc TEST37 tests STIRLING_NUMBERS1.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    06 July 2007
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer maxm
      parameter ( maxm = 6 )
      integer maxn
      parameter ( maxn = 6 )

      integer i
      integer s(0:maxm,0:maxn)

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST37'
      write ( *, '(a)' ) 
     &  '  STIRLING_NUMBERS1 computes a table of Stirling'
      write ( *, '(a)' ) '  numbers of the first kind.'

      call stirling_numbers1 ( maxm, maxn, s )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 
     &  '    I S(I,0) S(I,1) S(I,2) S(I,3) S(I,4) S(I,5)'
      write ( *, '(a)' ) ' '

      do i = 0, maxm
        write ( *, '(11i5)' ) i, s(i,0:maxn)
      end do

      return
      end 
      subroutine test38 ( )

c*********************************************************************72
c
cc TEST38 tests STIRLING_NUMBERS2.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    06 July 2007
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer maxm
      parameter ( maxm = 6 )
      integer maxn
      parameter ( maxn = 6 )

      integer i
      integer s(0:maxm,0:maxn)

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST38'
      write ( *, '(a)' ) 
     &  '  STIRLING_NUMBERS2 computes a table of Stirling'
      write ( *, '(a)' ) '  numbers of the second kind.'

      call stirling_numbers2 ( maxm, maxn, s )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 
     &  '    I S(I,0) S(I,1) S(I,2) S(I,3) S(I,4) S(I,5)'
      write ( *, '(a)' ) ' '

      do i = 0, maxm
        write ( *, '(11i5)' ) i, s(i,0:maxn)
      end do

      return
      end 
      subroutine test39 ( )

c*********************************************************************72
c
cc TEST39 tests SUBSET_COLEX_RANK, _COLEX_SUCCESSOR, _COLEX_UNRANK, _ENUM.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    06 July 2007
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n
      parameter ( n = 5 )

      integer nsub
      integer rank
      integer rank_old
      integer t(n)

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST39'
      write ( *, '(a)' ) '  All subsets of a set,'
      write ( *, '(a)' ) '  using the colexicographic ordering:'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  SUBSET_COLEX_RANK ranks,'
      write ( *, '(a)' ) '  SUBSET_COLEX_SUCCESSOR lists,'
      write ( *, '(a)' ) '  SUBSET_COLEX_UNRANK unranks.'
      write ( *, '(a)' ) '  SUBSET_ENUM enumerates.'
c
c  Enumerate.
c
      call subset_enum ( n, nsub )

      write ( *, '(a)' ) ' '
      write ( *, '(a,i8)' ) '  For N = ', n
      write ( *, '(a,i8)' ) '  the number of subsets is ', nsub
      write ( *, '(a)' ) ' '
c
c  List
c
      rank = -1

10    continue

        rank_old = rank

        call subset_colex_successor ( n, t, rank )

        if ( rank .le. rank_old ) then
          go to 20
        end if

        write ( *, '(6i5)' ) rank, t(1:n)

      go to 10

20    continue
c
c  Unrank.
c
      rank = nsub / 3

      call subset_colex_unrank ( rank, n, t )

      write ( *, '(a)' ) ' '
      write ( *, '(a,i8)' ) '  The element of rank ', rank
      write ( *, '(a)' ) ' '
      write ( *, '(6i5)' ) t(1:n)
c
c  Rank.
c
      call subset_colex_rank ( n, t, rank )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  The rank of the element:'
      write ( *, '(a)' ) ' '
      write ( *, '(6i5)' ) t(1:n)
      write ( *, '(a)' ) ' '
      write ( *, '(a,i8)' ) '  is computed as ', rank

      return
      end 
      subroutine test40 ( )

c*********************************************************************72
c
cc TEST40 tests SUBSET_ENUM, _LEX_RANK, _LEX_SUCCESSOR, _LEX_UNRANK.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    06 July 2007
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n
      parameter ( n = 5 )

      integer nsub
      integer rank
      integer rank_old
      integer t(n)

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST40'
      write ( *, '(a)' ) '  All subsets of a set,'
      write ( *, '(a)' ) '  using the lexicographic ordering:'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  SUBSET_ENUM enumerates,'
      write ( *, '(a)' ) '  SUBSET_LEX_RANK ranks,'
      write ( *, '(a)' ) '  SUBSET_LEX_SUCCESSOR lists,'
      write ( *, '(a)' ) '  SUBSET_LEX_UNRANK unranks.'
c
c  Enumerate.
c
      call subset_enum ( n, nsub )

      write ( *, '(a)' ) ' '
      write ( *, '(a,i8)' ) '  For N = ', n
      write ( *, '(a,i8)' ) '  the number of subsets is ', nsub
      write ( *, '(a)' ) ' '
c
c  List
c
      rank = -1

10    continue

        rank_old = rank

        call subset_lex_successor ( n, t, rank )

        if ( rank .le. rank_old ) then
          go to 20
        end if

        write ( *, '(6i5)' ) rank, t(1:n)

      go to 10

20    continue
c
c  Unrank.
c
      rank = nsub / 3

      call subset_lex_unrank ( rank, n, t )

      write ( *, '(a)' ) ' '
      write ( *, '(a,i8)' ) '  The element of rank ', rank
      write ( *, '(a)' ) ' '
      write ( *, '(6i5)' ) t(1:n)
c
c  Rank.
c
      call subset_lex_rank ( n, t, rank )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  The rank of the element:'
      write ( *, '(a)' ) ' '
      write ( *, '(6i5)' ) t(1:n)
      write ( *, '(a)' ) ' '
      write ( *, '(a,i8)' ) '  is computed as ', rank

      return
      end 
      subroutine test41 ( )

c*********************************************************************72
c
cc TEST41 tests SUBSETSUM_SWAP.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    06 July 2007
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n
      parameter ( n = 7 )

      integer a(n)
      integer i
      integer index(n)
      integer sum_achieved
      integer sum_desired

      save a

      data a / 12, 8, 11, 30, 8, 3, 7 /

      sum_desired = 17

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST41'
      write ( *, '(a)' ) 
     &  '  SUBSETSUM_SWAP seeks a solution of the subset'
      write ( *, '(a)' ) '  sum problem using pair swapping.'
      write ( *, '(a)' ) ' '
      write ( *, '(a,i4)' ) '  The desired sum is ', sum_desired

      call subsetsum_swap ( n, a, sum_desired, index, sum_achieved )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '    A(I), INDEX(I)'
      write ( *, '(a)' ) ' '

      do i = 1, n
        write ( *, '(2x,2i5)' ) a(i), index(i)
      end do

      write ( *, '(a)' ) ' '
      write ( *, '(a,i8)' ) '  The achieved sum is ', sum_achieved

      return
      end 
      subroutine test42 ( )

c*********************************************************************72
c
cc TEST42 tests TREE_ENUM, TREE_RANK, TREE_SUCCESSOR, TREE_UNRANK.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    06 July 2007
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n
      parameter ( n = 4 )

      integer rank
      integer rank_old
      integer t(2,n-1)
      integer tree_num

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST42'
      write ( *, '(a)' ) '  Trees:'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  TREE_ENUM enumerates,'
      write ( *, '(a)' ) '  TREE_RANK ranks,'
      write ( *, '(a)' ) '  TREE_SUCCESSOR lists,'
      write ( *, '(a)' ) '  TREE_UNRANK unranks.'
c
c  Enumerate.
c
      call tree_enum ( n, tree_num )

      write ( *, '(a)' ) ' '
      write ( *, '(a,i8)' ) '  For N = ', n
      write ( *, '(a,i8)' ) '  the number of trees is ', tree_num
      write ( *, '(a)' ) ' '
c
c  List
c
      rank = -1

10    continue

        rank_old = rank

        call tree_successor ( n, t, rank )

        if ( rank .le. rank_old ) then
          go to 20
        end if

        write ( *, '(i5,2x,5i5)' ) rank, t(1,1:n-1)
        write ( *, '(5x,2x,5i5)' )       t(2,1:n-1)

      go to 10

20    continue
c
c  Unrank.
c
      rank = tree_num / 2

      call tree_unrank ( rank, n, t )

      write ( *, '(a)' ) ' '
      write ( *, '(a,i8)' ) '  The element of rank ', rank
      write ( *, '(a)' ) ' '
      write ( *, '(2x,5i5)' ) t(1,1:n-1)
      write ( *, '(2x,5i5)' ) t(2,1:n-1)
c
c  Rank.
c
      call tree_rank ( n, t, rank )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  The rank of the element:'
      write ( *, '(a)' ) ' '
      write ( *, '(2x,5i5)' ) t(1,1:n-1)
      write ( *, '(2x,5i5)' ) t(2,1:n-1)
      write ( *, '(a)' ) ' '
      write ( *, '(a,i8)' ) '  is computed as ', rank

      return
      end 
