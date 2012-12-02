      program main

c*********************************************************************72
c
cc MAIN is the main program for TEST_CON_PRB.
c
c  Discussion:
c
c    TEST_CON_PRB demonstrates the TEST_CON continuation test problems.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    25 September 2008
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer problem
      integer problem_num

      call timestamp ( )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST_CON_PRB'
      write ( *, '(a)' ) '  FORTRAN77 version'
      write ( *, '(a)' ) '  Test the TEST_CON library.'
c
c  Find out how many problems are available.
c
      call p00_problem_num ( problem_num )

      write ( *, '(a)' ) ' '
      write ( *, '(a,i8,a)' )
     &  '  There are ', problem_num, ' test functions.'
c
c  Print the number of options for each problem.
c
      call p00_option_num_test ( problem_num )
c
c  Print the title of each problem.
c
      call p00_title_test ( problem_num )
c
c  Print the size of each problem.
c
      call p00_nvar_test ( problem_num )
c
c  Get the starting point X0 and the norm of F(X0).
c
c     call p00_start_test ( problem_num )
c
c  Check the jacobian near the starting point.
c
c     call p00_jac_test ( problem_num )
c
c  Solve each problem.
c
c     call p00_newton_test ( problem_num )
c
c  Get the starting stepsize.
c
      call p00_stepsize_test ( problem_num )
c
c  Terminate.
c
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST_CON_PRB'
      write ( *, '(a)' ) '  Normal end of execution.'
      write ( *, '(a)' ) ' '
      call timestamp ( )

      stop
      end
      subroutine p00_nvar_test ( problem_num )

c*********************************************************************72
c
cc P00_NVAR_TEST prints the problem size.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    18 September 2008
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer PROBLEM_NUM, the number of problems.
c
      implicit none

      integer nvar
      integer option
      integer option_num
      integer problem
      integer problem_num

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'P00_NVAR_TEST'
      write ( *, '(a)' ) '  List the problem size.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '   Problem    Option      Size'
      write ( *, '(a)' ) ' '

      do problem = 1, problem_num
        call p00_option_num ( problem, option_num )
        write ( *, '(a)' ) ' '
        do option = 1, option_num
          call p00_nvar ( problem, option, nvar )
          write ( *, '(2x,i8,2x,i8,2x,i8)' ) problem, option, nvar
        end do
      end do

      return
      end
      subroutine p00_option_num_test ( problem_num )

c*********************************************************************72
c
cc P00_OPTION_NUM_TEST lists the number of options.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    17 September 2008
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer PROBLEM_NUM, the number of problems.
c
      implicit none

      integer option_num
      integer problem
      integer problem_num

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'P00_OPTION_NUM_TEST'
      write ( *, '(a)' )
     &  '  List the number of options for each problem.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '   Problem   Options'
      write ( *, '(a)' ) ' '

      do problem = 1, problem_num
        call p00_option_num ( problem, option_num )
        write ( *, '(2x,i8,2x,i8)' ) problem, option_num
      end do

      return
      end
      subroutine p00_stepsize_test ( problem_num )

c*********************************************************************72
c
cc P00_STEPSIZE_TEST prints the stepsizes for each problem.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    25 September 2008
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer PROBLEM_NUM, the number of problems.
c
      implicit none

      double precision h
      double precision hmax
      double precision hmin
      integer option
      integer option_num
      integer problem
      integer problem_num

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'P00_STEPSIZE_TEST'
      write ( *, '(a)' ) '  Print the stepsizes for each problem.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) &
     &'   Problem    Option      H               HMIN             HMAX'
      write ( *, '(a)' ) ' '

      do problem = 1, problem_num

        call p00_option_num ( problem, option_num )
        write ( *, '(a)' ) ' '

        do option = 1, option_num

          call p00_stepsize ( problem, option, h, hmin, hmax )

          write ( *, '(2x,i8,2x,i8,2x,g14.6,2x,g14.6,2x,g14.6)' )
     &      problem, option, h, hmin, hmax

        end do
      end do

      return
      end
      subroutine p00_title_test ( problem_num )

c*********************************************************************72
c
cc P00_TITLE_TEST prints the problem titles.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    18 September 2008
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer PROBLEM_NUM, the number of problems.
c
      implicit none

      integer option
      integer option_num
      integer problem
      integer problem_num
      character*80 title

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'P00_TITLE_TEST'
      write ( *, '(a)' ) '  List the problem title'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '   Problem    Option  Title'
      write ( *, '(a)' ) ' '

      do problem = 1, problem_num
        call p00_option_num ( problem, option_num )
        write ( *, '(a)' ) ' '
        do option = 1, option_num
          call p00_title ( problem, option, title )
          write ( *, '(2x,i8,2x,i8,2x,a)' )
     &      problem, option, trim ( title )
        end do
      end do

      return
      end
