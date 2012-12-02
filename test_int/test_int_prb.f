      program main

c*********************************************************************72
c
cc MAIN is the main program for TEST_INT_PRB.
c
c  Discussion:
c
c    TEST_INT_PRB demonstrates the use of the TEST_INT integration
c    test functions.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    01 September 2012
c
c  Author:
c
c    John Burkardt
c
      implicit none

      call timestamp ( )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST_INT_PRB'
      write ( *, '(a)' ) '  FORTRAN77 version'
      write ( *, '(a)' ) '  Test the TEST_INT library.'

      call test01 ( )
      call test02 ( )
      call test03 ( )
      call test04 ( )
      call test05 ( )
      call test06 ( )
      call test07 ( )
c
c  Terminate.
c
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST_INT_PRB'
      write ( *, '(a)' ) '  Normal end of execution.'

      write ( *, '(a)' ) ' '
      call timestamp ( )

      stop
      end
      subroutine test01 ( )

c*********************************************************************72
c
cc TEST01 applies a composite midpoint rule to finite interval 1D problems.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    06 October 2006
c
c  Author:
c
c    John Burkardt
c
      implicit none

      double precision error
      double precision exact
      integer int_log
      integer int_num
      integer prob
      integer prob_num
      double precision result

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST01'
      write ( *, '(a)' ) '  Composite midpoint rule,'
      write ( *, '(a)' ) '  for 1D finite interval problems.'

      call p00_prob_num ( prob_num )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Problem       Exact'
      write ( *, '(a)' ) '         Ints   Approx       Error'
c
c  Pick a problem.
c
      do prob = 1, prob_num

        call p00_exact ( prob, exact )

        write ( *, '(a)' ) ' '
        write ( *, '(i6,2x,4x,2x,g14.6)' ) prob, exact
c
c  Pick a number of subintervals.
c
        do int_log = 0, 7

          int_num = 2**int_log

          call p00_midpoint ( prob, int_num, result )

          error = abs ( exact - result )

          write ( *, '(6x,2x,i4,2x,2g14.6)' ) int_num, result, error

        end do

      end do

      return
      end
      subroutine test02 ( )

c*********************************************************************72
c
cc TEST02 applies a composite Simpson rule to finite interval 1D problems.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    06 October 2006
c
c  Author:
c
c    John Burkardt
c
      implicit none

      double precision error
      double precision exact
      integer int_log
      integer int_num
      integer prob
      integer prob_num
      double precision result

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST02'
      write ( *, '(a)' ) '  Composite Simpson rule,'
      write ( *, '(a)' ) '  for 1D finite interval problems.'

      call p00_prob_num ( prob_num )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Problem       Exact'
      write ( *, '(a)' ) '         Ints   Approx       Error'
c
c  Pick a problem.
c
      do prob = 1, prob_num
c
c  Some problems have singularities that kill the calculation.
c
        call p00_exact ( prob, exact )

        write ( *, '(a)' ) ' '
        write ( *, '(i6,2x,4x,2x,g14.6)' ) prob, exact
c
c  Pick a number of subintervals.
c
        do int_log = 0, 7

          int_num = 2**int_log

          call p00_simpson ( prob, int_num, result )

          error = abs ( exact - result )

          write ( *, '(6x,2x,i4,2x,2g14.6)' ) int_num, result, error

        end do

      end do

      return
      end
      subroutine test03 ( )

c*********************************************************************72
c
cc TEST03 applies a Monte Carlo rule to finite interval 1D problems.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    06 October 2006
c
c  Author:
c
c    John Burkardt
c
      implicit none

      double precision error
      double precision exact
      integer int_log
      integer int_num
      integer prob
      integer prob_num
      double precision result

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST03'
      write ( *, '(a)' ) '  Monte Carlo rule,'
      write ( *, '(a)' ) '  for 1D finite interval problems.'

      call p00_prob_num ( prob_num )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Problem       Exact'
      write ( *, '(a)' ) '         Pts    Approx       Error'
c
c  Pick a problem.
c
      do prob = 1, prob_num

        call p00_exact ( prob, exact )

        write ( *, '(a)' ) ' '
        write ( *, '(i6,2x,4x,2x,g14.6)' ) prob, exact
c
c  Pick a number of points.
c
        do int_log = 0, 7

          int_num = 2**int_log

          call p00_montecarlo ( prob, int_num, result )

          error = abs ( exact - result )

          write ( *, '(6x,2x,i4,2x,2g14.6)' ) int_num, result, error

        end do

      end do

      return
      end
      subroutine test04 ( )

c*********************************************************************72
c
cc TEST04 applies a composite Gauss-Legendre rule to finite interval 1D problems.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    06 October 2006
c
c  Author:
c
c    John Burkardt
c
      implicit none

      double precision error
      double precision exact
      integer int_log
      integer int_num
      integer prob
      integer prob_num
      double precision result

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST04'
      write ( *, '(a)' ) 
     &  '  Use a composite 4 point Gauss-Legendre rule,'
      write ( *, '(a)' ) '  for 1D finite interval problems.'

      call p00_prob_num ( prob_num )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Problem       Exact'
      write ( *, '(a)' ) '         Ints   Approx	   Error'
c
c  Pick a problem.
c
      do prob = 1, prob_num

        call p00_exact ( prob, exact )

        write ( *, '(a)' ) ' '
        write ( *, '(i6,2x,4x,2x,g14.6)' ) prob, exact
c
c  Pick a number of subintervals.
c
        do int_log = 0, 7

          int_num = 2**int_log

          call p00_gauss_legendre ( prob, int_num, result )

          error = abs ( exact - result )

          write ( *, '(6x,2x,i4,2x,2g14.6)' ) int_num, result, error

        end do

      end do

      return
      end
      subroutine test05 ( )

c*********************************************************************72
c
cc TEST05 applies a composite trapezoid rule to finite interval 1D problems.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    14 November 2009
c
c  Author:
c
c    John Burkardt
c
      implicit none

      double precision error
      double precision exact
      integer int_log
      integer int_num
      integer prob
      integer prob_num
      double precision result

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST05'
      write ( *, '(a)' ) '  Composite trapezoid rule,'
      write ( *, '(a)' ) '  for 1D finite interval problems.'

      call p00_prob_num ( prob_num )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Problem       Exact'
      write ( *, '(a)' ) '         Ints   Approx       Error'
c
c  Pick a problem.
c
      do prob = 1, prob_num

        call p00_exact ( prob, exact )

        write ( *, '(a)' ) ' '
        write ( *, '(i6,2x,4x,2x,g14.6)' ) prob, exact
c
c  Pick a number of subintervals.
c
        do int_log = 0, 7

          int_num = 2**int_log

          call p00_trapezoid ( prob, int_num, result )

          error = abs ( exact - result )

          write ( *, '(6x,2x,i4,2x,2g14.6)' ) int_num, result, error

        end do

      end do

      return
      end
      subroutine test06 ( )

c*********************************************************************72
c
cc TEST06 applies a Halton sequence rule to finite interval 1D problems.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    06 October 2006
c
c  Author:
c
c    John Burkardt
c
      implicit none

      double precision error
      double precision exact
      integer int_log
      integer int_num
      integer prob
      integer prob_num
      double precision result

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST06'
      write ( *, '(a)' ) '  Halton sequence rule,'
      write ( *, '(a)' ) '  for 1D finite interval problems.'

      call p00_prob_num ( prob_num )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Problem       Exact'
      write ( *, '(a)' ) '         Pts    Approx       Error'
c
c  Pick a problem.
c
      do prob = 1, prob_num

        call p00_exact ( prob, exact )

        write ( *, '(a)' ) ' '
        write ( *, '(i6,2x,4x,2x,g14.6)' ) prob, exact
c
c  Pick a number of points.
c
        do int_log = 0, 7

          int_num = 2**int_log

          call p00_halton ( prob, int_num, result )

          error = abs ( exact - result )

          write ( *, '(6x,2x,i4,2x,2g14.6)' ) int_num, result, error

        end do

      end do

      return
      end
      subroutine test07 ( )

c*********************************************************************72
c
cc TEST07 applies an evenly spaced point rule to finite interval 1D problems.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    14 November 2009
c
c  Author:
c
c    John Burkardt
c
      implicit none

      double precision error
      double precision exact
      integer int_log
      integer int_num
      integer prob
      integer prob_num
      double precision result

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST07'
      write ( *, '(a)' ) '  Evenly spaced point sequence rule,'
      write ( *, '(a)' ) '  for 1D finite interval problems.'

      call p00_prob_num ( prob_num )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Problem       Exact'
      write ( *, '(a)' ) '         Pts    Approx       Error'
c
c  Pick a problem.
c
      do prob = 1, prob_num

        call p00_exact ( prob, exact )

        write ( *, '(a)' ) ' '
        write ( *, '(i6,2x,4x,2x,g14.6)' ) prob, exact
c
c  Pick a number of points.
c
        do int_log = 0, 7

          int_num = 2**int_log

          call p00_even ( prob, int_num, result )

          error = abs ( exact - result )

          write ( *, '(6x,2x,i4,2x,2g14.6)' ) int_num, result, error

        end do

      end do

      return
      end
