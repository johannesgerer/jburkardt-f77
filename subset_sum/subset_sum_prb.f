      program main

c*********************************************************************72
c
cc MAIN is the main program for SUBSET_SUM_PRB.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    18 May 2012
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n_max
      parameter ( n_max = 21 )
      integer ind
      integer ind_max
      integer ind_min
      integer n
      integer t
      integer test
      integer test_num
      parameter ( test_num = 9 )
      integer w(n_max)

      call timestamp ( )
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'SUBSET_SUM_TEST:'
      write ( *, '(a)' ) '  FORTRAN77 version'
      write ( *, '(a)' ) '  Test the SUBSET_SUM library.'
c
c  Find individual solutions.
c
      do test = 1, test_num

        if ( test .eq. 1 ) then
          n = 8
          w( 1) = 15
          w( 2) = 22
          w( 3) = 14
          w( 4) = 26
          w( 5) = 32
          w( 6) =  9
          w( 7) = 16
          w( 8) =  8
          t = 53
          ind_min = 0
          ind_max = 2**n - 1
        else if ( test .eq. 2 ) then
          n = 8
          w( 1) = 15
          w( 2) = 22
          w( 3) = 14
          w( 4) = 26
          w( 5) = 32
          w( 6) =  9
          w( 7) = 16
          w( 8) =  8
          t = 53
          ind_min = ind + 1
          ind_max = 2**n - 1
        else if ( test .eq. 3 ) then
          n = 8
          w( 1) = 15
          w( 2) = 22
          w( 3) = 14
          w( 4) = 26
          w( 5) = 32
          w( 6) =  9
          w( 7) = 16
          w( 8) =  8
          t = 53
          ind_min = ind + 1
          ind_max = 2**n - 1
        else if ( test .eq. 4 ) then
          n = 10
          w( 1) =  267
          w( 2) =  493
          w( 3) =  869
          w( 4) =  961
          w( 5) = 1000
          w( 6) = 1153
          w( 7) = 1246
          w( 8) = 1598
          w( 9) = 1766
          w(10) = 1922
          t = 5842
          ind_min = 0
          ind_max = 2**n - 1
        else if ( test .eq. 5 ) then
          n = 21
          w( 1) =  518533
          w( 2) = 1037066
          w( 3) = 2074132
          w( 4) = 1648264
          w( 5) =  796528
          w( 6) = 1593056
          w( 7) =  686112
          w( 8) = 1372224
          w( 9) =  244448
          w(10) =  488896
          w(11) =  977792
          w(12) = 1955584
          w(13) = 1411168
          w(14) =  322336
          w(15) =  644672
          w(16) = 1289344
          w(17) =   78688
          w(18) =  157376
          w(19) =  314752
          w(20) =  629504
          w(21) = 1259008
          t = 2463098
          ind_min = 0
          ind_max = 2**n - 1
        else if ( test .eq. 6 ) then
          n = 10
          w( 1) = 41
          w( 2) = 34
          w( 3) = 21
          w( 4) = 20
          w( 5) =  8
          w( 6) =  7
          w( 7) =  7
          w( 8) =  4
          w( 9) =  3
          w(10) =  3
          t = 50
          ind_min = 0
          ind_max = 2**n - 1
        else if ( test .eq. 7 ) then
          n = 9
          w( 1) = 81
          w( 2) = 80
          w( 3) = 43
          w( 4) = 40
          w( 5) = 30
          w( 6) = 26
          w( 7) = 12
          w( 8) = 11
          w( 9) =  9
          t = 100
          ind_min = 0
          ind_max = 2**n - 1
        else if ( test .eq. 8 ) then
          n = 6
          w( 1) =  1
          w( 2) =  2
          w( 3) =  4
          w( 4) =  8
          w( 5) = 16
          w( 6) = 32
          t = 22
          ind_min = 0
          ind_max = 2**n - 1
        else if ( test .eq. 9 ) then
          n = 10
          w( 1) = 25
          w( 2) = 27
          w( 3) =  3
          w( 4) = 12
          w( 5) =  6
          w( 6) = 15
          w( 7) =  9
          w( 8) = 30
          w( 9) = 21
          w(10) = 19
          t = 50
          ind_min = 0
          ind_max = 2**n - 1
        end if

        call test01 ( n, w, t, ind_min, ind_max, ind )

      end do
c
c  Simply count solutions.
c
      do test = 1, test_num

        if ( test .eq. 1 ) then
          n = 8
          w( 1) = 15
          w( 2) = 22
          w( 3) = 14
          w( 4) = 26
          w( 5) = 32
          w( 6) =  9
          w( 7) = 16
          w( 8) =  8
          t = 53
          ind_min = 0
          ind_max = 2**n - 1
        else if ( test .eq. 2 ) then
          n = 8
          w( 1) = 15
          w( 2) = 22
          w( 3) = 14
          w( 4) = 26
          w( 5) = 32
          w( 6) =  9
          w( 7) = 16
          w( 8) =  8
          t = 53
          ind_min = 68
          ind_max = 2**n - 1
        else if ( test .eq. 3 ) then
          n = 8
          w( 1) = 15
          w( 2) = 22
          w( 3) = 14
          w( 4) = 26
          w( 5) = 32
          w( 6) =  9
          w( 7) = 16
          w( 8) =  8
          t = 53
          ind_min = 167
          ind_max = 2**n - 1
        else if ( test .eq. 4 ) then
          n = 10
          w( 1) =  267
          w( 2) =  493
          w( 3) =  869
          w( 4) =  961
          w( 5) = 1000
          w( 6) = 1153
          w( 7) = 1246
          w( 8) = 1598
          w( 9) = 1766
          w(10) = 1922
          t = 5842
          ind_min = 0
          ind_max = 2**n - 1
        else if ( test .eq. 5 ) then
          n = 21
          w( 1) =  518533
          w( 2) = 1037066
          w( 3) = 2074132
          w( 4) = 1648264
          w( 5) =  796528
          w( 6) = 1593056
          w( 7) =  686112
          w( 8) = 1372224
          w( 9) =  244448
          w(10) =  488896
          w(11) =  977792
          w(12) = 1955584
          w(13) = 1411168
          w(14) =  322336
          w(15) =  644672
          w(16) = 1289344
          w(17) =   78688
          w(18) =  157376
          w(19) =  314752
          w(20) =  629504
          w(21) = 1259008
          t = 2463098
          ind_min = 0
          ind_max = 2**n - 1
        else if ( test .eq. 6 ) then
          n = 10
          w( 1) = 41
          w( 2) = 34
          w( 3) = 21
          w( 4) = 20
          w( 5) =  8
          w( 6) =  7
          w( 7) =  7
          w( 8) =  4
          w( 9) =  3
          w(10) =  3
          t = 50
          ind_min = 0
          ind_max = 2**n - 1
        else if ( test .eq. 7 ) then
          n = 9
          w( 1) = 81
          w( 2) = 80
          w( 3) = 43
          w( 4) = 40
          w( 5) = 30
          w( 6) = 26
          w( 7) = 12
          w( 8) = 11
          w( 9) =  9
          t = 100
          ind_min = 0
          ind_max = 2**n - 1
        else if ( test .eq. 8 ) then
          n = 6
          w( 1) =  1
          w( 2) =  2
          w( 3) =  4
          w( 4) =  8
          w( 5) = 16
          w( 6) = 32
          t = 22
          ind_min = 0
          ind_max = 2**n - 1
        else if ( test .eq. 9 ) then
          n = 10
          w( 1) = 25
          w( 2) = 27
          w( 3) =  3
          w( 4) = 12
          w( 5) =  6
          w( 6) = 15
          w( 7) =  9
          w( 8) = 30
          w( 9) = 31
          w(10) = 19
          t = 50
          ind_min = 0
          ind_max = 2**n - 1
        end if

        call test02 ( n, w, t, ind_min, ind_max )

      end do
c
c  Terminate.
c
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'SUBSET_SUM_TEST'
      write ( *, '(a)' ) '  Normal end of execution.'
      write ( *, '(a)' ) ' '
      call timestamp ( )

      return
      end
      subroutine test01 ( n, w, t, ind_min, ind_max, ind )

c*********************************************************************72
c
cc TEST01 seeks a subset of a set that has a given sum.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    02 February 2012
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the number of weights.
c
c    Input, integer W(N), a set of weights.  The length of this
c    array must be no more than 31.
c
c    Input, integer T, the target value.
c
c    Input, integer IND_MIN, IND_MAX, the lower and upper
c    limits to be searched.
c
c    Output, integer IND, the index of a solution, if found,
c    or the value -1 otherwise.
c
      implicit none

      integer n

      integer c(n)
      integer i
      integer ind
      integer ind_max
      integer ind_min
      integer t
      integer w(n)

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST01:'
      write ( *, '(a)' ) '  Seek a subset of W that sums to T.'
      write ( *, '(a)' ) ' '
      write ( *, '(a,i8)' ) '  Target value T = ', t
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '   I       W(I)'
      write ( *, '(a)' ) ' '
      do i = 1, n
        write ( *, '(2x,i2,2x,i8)' ) i, w(i)
      end do

      call subset_sum_find ( n, w, t, ind_min, ind_max, ind, c )

      if ( ind .eq. -1 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) '  No solution was found.'
        return
      end if

      write ( *, '(a)' ) ' '
      write ( *, '(a,i8)' )  '  Solution index = ', ind
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '   I       W(I)  C(I)'
      write ( *, '(a)' ) ' '
      do i = 1, n
        write ( *, '(2x,i2,2x,i8,5x,i1)' ) i, w(i), c(i)
      end do

      return
      end
      subroutine test02 ( n, w, t, ind_min, ind_max )

c*********************************************************************72
c
cc TEST02 counts solutions to the subset sum problem.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    02 February 2012
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the number of weights.
c
c    Input, integer W(N), a set of weights.  The length of this
c    array must be no more than 31.
c
c    Input, integer T, the target value.
c
c    Input, integer IND_MIN, IND_MAX, the lower and upper
c    limits to be searched.
c
      implicit none

      integer n

      integer i
      integer ind_max
      integer ind_min
      integer solution_num
      integer t
      integer w(n)

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST02:'
      write ( *, '(a)' ) '  Count solutions to the subset sum problem.'
      write ( *, '(a)' ) ' '
      write ( *, '(a,i8)' ) '  Target value T = ', t
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '   I       W(I)'
      write ( *, '(a)' ) ' '
      do i = 1, n
        write ( *, '(2x,i2,2x,i8)' ) i, w(i)
      end do

      call subset_sum_count ( n, w, t, ind_min, ind_max, solution_num )

      write ( *, '(a)' ) ' '
      write ( *, '(a,i8)' ) '  Number of solutions = ', solution_num

      return
      end
