      program main

!*********************************************************************72
!
!! MAIN is the main program for PARTITION_PROBLEM_PRB.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    11 May 2012
!
!  Author:
!
!    John Burkardt
!
      implicit none

      integer n_max
      parameter ( n_max = 10 )

      integer n
      integer test
      integer test_num
      parameter ( test_num = 5 )
      integer w(n_max)

      call timestamp ( )
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'PARTITION_PROBLEM_PRB:'
      write ( *, '(a)' ) '  FORTRAN77 version'
      write ( *, '(a)' ) '  Test the PARTITION_PROBLEM library.'
!
!  Find individual solutions.
!
      do test = 1, test_num

        if ( test .eq. 1 ) then
          n = 5
          w(1) = 19
          w(2) = 17
          w(3) = 13
          w(4) =  9
          w(5) =  6
        else if ( test .eq. 2 ) then
          n = 9
          w(1) = 484
          w(2) = 114
          w(3) = 205
          w(4) = 288
          w(5) = 506
          w(6) = 503
          w(7) = 201
          w(8) = 127
          w(9) = 410
        else if ( test .eq. 3 ) then
          n = 10
          w(1) =  771
          w(2) =  121
          w(3) =  281
          w(4) =  854
          w(5) =  885
          w(6) =  734
          w(7) =  486
          w(8) = 1003
          w(9) =   83
          w(10) =  62
        else if ( test .eq. 4 ) then
          n = 10
          w(1) =  2
          w(2) = 10
          w(3) =  3
          w(4) =  8
          w(5) =  5
          w(6) =  7
          w(7) =  9
          w(8) =  5
          w(9) =  3
          w(10) = 2
        else if ( test .eq. 5 ) then
          n = 9
          w(1) = 3
          w(2) = 4
          w(3) = 3
          w(4) = 1
          w(5) = 3
          w(6) = 2
          w(7) = 3
          w(8) = 2
          w(9) = 1
        end if

        call test01 ( n, w )

      end do
!
!  Count solutions.
!
      do test = 1, test_num

        if ( test .eq. 1 ) then
          n = 5
          w(1) = 19
          w(2) = 17
          w(3) = 13
          w(4) =  9
          w(5) =  6
        else if ( test .eq. 2 ) then
          n = 9
          w(1) = 484
          w(2) = 114
          w(3) = 205
          w(4) = 288
          w(5) = 506
          w(6) = 503
          w(7) = 201
          w(8) = 127
          w(9) = 410
        else if ( test .eq. 3 ) then
          n = 10
          w(1) =  771
          w(2) =  121
          w(3) =  281
          w(4) =  854
          w(5) =  885
          w(6) =  734
          w(7) =  486
          w(8) = 1003
          w(9) =   83
          w(10) =  62
        else if ( test .eq. 4 ) then
          n = 10
          w(1) =  2
          w(2) = 10
          w(3) =  3
          w(4) =  8
          w(5) =  5
          w(6) =  7
          w(7) =  9
          w(8) =  5
          w(9) =  3
          w(10) = 2
        else if ( test .eq. 5 ) then
          n = 9
          w(1) = 3
          w(2) = 4
          w(3) = 3
          w(4) = 1
          w(5) = 3
          w(6) = 2
          w(7) = 3
          w(8) = 2
          w(9) = 1
        end if

        call test02 ( n, w )

      end do
!
!  Terminate.
!
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'PARTITION_PROBLEM_PRB'
      write ( *, '(a)' ) '  Normal end of execution.'
      write ( *, '(a)' ) ' '
      call timestamp ( )

      return
      end
      subroutine test01 ( n, w )

!*********************************************************************72
!
!! TEST01 tests PARTITION_BRUTE.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    11 May 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of weights.
!
!    Input, integer W(N), a set of weights.
!
      implicit none

      integer n

      integer c(n)
      integer discrepancy
      integer i
      integer w(n)
      integer w0_sum
      integer w1_sum

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST01:'
      write ( *, '(a)' ) 
     &  '  Partition a set of N integers W so that the subsets'
      write ( *, '(a)' ) '  have equal sums.'

      call partition_brute ( n, w, c, discrepancy )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '     I        W0        W1'
      write ( *, '(a)' ) ' '
      w0_sum = 0
      w1_sum = 0
      do i = 1, n
        if ( c(i) .eq. 0 ) then
          w0_sum = w0_sum + w(i)
          write ( *, '(2x,i4,2x,i8,2x,8x)' ) i, w(i)
        else
          w1_sum = w1_sum + w(i)
          write ( *, '(2x,i4,2x,8x,2x,i8)' ) i, w(i)
        end if
      end do
      write ( *, '(a)' ) '        --------  --------'
      write ( *, '(2x,4x,2x,i8,2x,i8)' ) w0_sum, w1_sum
      write ( *, '(a)' ) ' '
      write ( *, '(a,i8)' ) '  Discrepancy = ', discrepancy

      return
      end
      subroutine test02 ( n, w )

!*********************************************************************72
!
!! TEST02 tests PARTITION_COUNT.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    11 May 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of weights.
!
!    Input, integer W(N), a set of weights.
!
      implicit none

      integer n

      integer c(n)
      integer count
      integer i
      integer w(n)

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST02:'
      write ( *, '(a)' ) 
     &  '  PARTITION_COUNT counts the number of exact solutions'
      write ( *, '(a)' ) '  of the partition problem.'

      call partition_count ( n, w, count )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '     I        W'
      write ( *, '(a)' ) ' '
      do i = 1, n
        write ( *, '(2x,i4,2x,i8)' ) i, w(i)
      end do
      write ( *, '(a)' ) ' '
      write ( *, '(2x,a,i4)' ) 'Number of solutions = ', count

      return
      end
