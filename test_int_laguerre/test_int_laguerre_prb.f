      program main

c*********************************************************************72
c
cc MAIN is the main program for TEST_INT_LAGUERRE_PRB.
c
c  Discussion:
c
c    TEST_INT_LAGUERRE_PRB demonstrates the use of the TEST_INT_LAGUERRE
c    integration test functions.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    13 September 2012
c
c  Author:
c
c    John Burkardt
c
      implicit none

      call timestamp ( )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST_INT_LAGUERRE_PRB'
      write ( *, '(a)' ) '  FORTRAN77 version'
      write ( *, '(a)' ) '  Test the TEST_INT_LAGUERRE library.'

      call test01 ( )
      call test02 ( )
      call test03 ( )
      call test04 ( )
      call test05 ( )
c
c  Terminate.
c
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST_INT_LAGUERRE_PRB'
      write ( *, '(a)' ) '  Normal end of execution.'

      write ( *, '(a)' ) ' '
      call timestamp ( )

      stop
      end
      subroutine test01 ( )

c*********************************************************************72
c
cc TEST01 tests P00_PROBLEM_NUM and P00_TITLE.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    27 July 2007
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer problem
      integer problem_num
      character * ( 80 ) title

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST01'
      write ( *, '(a)' ) 
     &  '  P00_PROBLEM_NUM returns the number of problems.'
      write ( *, '(a)' ) '  P00_TITLE returns the title of a problem.'

      call p00_problem_num ( problem_num )

      write ( *, '(a)' ) ' '
      write ( *, '(a,i8)' ) 
     &  '  P00_PROBLEM_NUM: number of problems is ', problem_num
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '   Problem       Title'
      write ( *, '(a)' ) ' '

      do problem = 1, problem_num

        call p00_title ( problem, title )

        write ( *, '(2x,i8,2x,a)' ) 
     &    problem, '"' // trim ( title ) // '"'

      end do

      return
      end
      subroutine test02 ( )

c*********************************************************************72
c
cc TEST02 tests P00_ALPHA and P00_EXACT.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    27 July 2007
c
c  Author:
c
c    John Burkardt
c
      implicit none

      double precision alpha
      double precision exact
      integer problem
      integer problem_num

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST02'
      write ( *, '(a)' ) 
     &  '  P00_ALPHA returns the lower limit of integration.'
      write ( *, '(a)' ) '  P00_EXACT returns the "exact" integral.'

      call p00_problem_num ( problem_num )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '   Problem       ALPHA          EXACT'
      write ( *, '(a)' ) ' '

      do problem = 1, problem_num

        call p00_alpha ( problem, alpha )

        call p00_exact ( problem, exact )

        write ( *, '(2x,i8,2x,g14.6,2x,g24.16)' ) problem, alpha, exact

      end do

      return
      end
      subroutine test03 ( )

c*********************************************************************72
c
cc TEST03 tests P00_GAUSS_LAGUERRE.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    26 December 2011
c
c  Author:
c
c    John Burkardt
c
      implicit none

      double precision error
      double precision estimate
      double precision exact
      integer order
      integer order_log
      integer problem
      integer problem_num

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST03'
      write ( *, '(a)' ) 
     &  '  P00_GAUSS_LAGUERRE applies a Gauss-Laguerre rule'
      write ( *, '(a)' ) '  to estimate an integral on [ALPHA,+oo).'

      call p00_problem_num ( problem_num )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '                          Exact/'
      write ( *, '(a)' ) 
     &  '   Problem     Order      Estimate        Error'

      do problem = 1, problem_num

        call p00_exact ( problem, exact )

        order = 1

        write ( *, '(a)' ) ' '
        write ( *, '(2x,i8,2x,8x,2x,g14.6)' ) problem, exact

        do order_log = 0, 6

          call p00_gauss_laguerre ( problem, order, estimate )

          error = abs ( exact - estimate )

          write ( *, '(2x,8x,2x,i8,2x,g14.6,2x,g14.6)' ) 
     &      order, estimate, error

          order = order * 2

        end do

      end do

      return
      end
      subroutine test04 ( )

c*********************************************************************72
c
cc TEST04 tests P00_EXP_TRANSFORM.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    26 December 2011
c
c  Author:
c
c    John Burkardt
c
      implicit none

      double precision error
      double precision estimate
      double precision exact
      integer order
      integer order_log
      integer problem
      integer problem_num

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST04'
      write ( *, '(a)' ) 
     &  '  P00_EXP_TRANSFORM applies an exponential transform'
      write ( *, '(a)' ) '  to estimate an integral on [ALPHA,+oo)'
      write ( *, '(a)' ) 
     &  '  as a transformed integral on (0,exp(-ALPHA)],'
      write ( *, '(a)' ) '  and applying a Gauss-Legendre rule.'

      call p00_problem_num ( problem_num )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '                          Exact/'
      write ( *, '(a)' ) 
     &  '   Problem     Order      Estimate        Error'

      do problem = 1, problem_num

        call p00_exact ( problem, exact )

        order = 1

        write ( *, '(a)' ) ' '
        write ( *, '(2x,i8,2x,8x,2x,g14.6)' ) problem, exact

        do order_log = 0, 6

          call p00_exp_transform ( problem, order, estimate )

          error = abs ( exact - estimate )

          write ( *, '(2x,8x,2x,i8,2x,g14.6,2x,g14.6)' ) 
     &      order, estimate, error

          order = order * 2

        end do

      end do

      return
      end
      subroutine test05 ( )

c*********************************************************************72
c
cc TEST05 tests P00_RAT_TRANSFORM.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    26 December 2011
c
c  Author:
c
c    John Burkardt
c
      implicit none

      double precision error
      double precision estimate
      double precision exact
      integer order
      integer order_log
      integer problem
      integer problem_num

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST05'
      write ( *, '(a)' ) 
     &  '  P00_RAT_TRANSFORM applies a rational transform'
      write ( *, '(a)' ) '  to estimate an integral on [ALPHA,+oo)'
      write ( *, '(a)' ) 
     &  '  as a transformed integral on (0,1/(1+ALPHA)],'
      write ( *, '(a)' ) '  and applying a Gauss-Legendre rule.'

      call p00_problem_num ( problem_num )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '                          Exact/'
      write ( *, '(a)' ) 
     &  '   Problem     Order      Estimate        Error'

      do problem = 1, problem_num

        call p00_exact ( problem, exact )

        order = 1

        write ( *, '(a)' ) ' '
        write ( *, '(2x,i8,2x,8x,2x,g14.6)' ) problem, exact

        do order_log = 0, 6

          call p00_rat_transform ( problem, order, estimate )

          error = abs ( exact - estimate )

          write ( *, '(2x,8x,2x,i8,2x,g14.6,2x,g14.6)' ) 
     &      order, estimate, error

          order = order * 2

        end do

      end do

      return
      end
