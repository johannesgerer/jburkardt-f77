      program main

c*********************************************************************72
c
cc MAIN is the main program for TEST_OPT_CON_PRB.
c
c  Discussion:
c
c    TEST_OPT_CON_PRB calls the TEST_OPT_CON tests.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    15 February 2012
c
c  Author:
c
c    John Burkardt
c
      implicit none

      call timestamp (  )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST_OPT_CON_PRB'
      write ( *, '(a)' ) '  FORTRAN77 version'
      write ( *, '(a)' ) '  Test the TEST_OPT_CON library.'

      call test01 ( )
      call test02 ( )
c
c  Terminate.
c
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST_OPT_CON_PRB'
      write ( *, '(a)' ) '  Normal end of execution.'

      write ( *, '(a)' ) ' '
      call timestamp (  )

      stop
      end
      subroutine test01 ( )

c*********************************************************************72
c
cc TEST01 simply prints the title of each problem.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    15 February 2012
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer problem_num
      integer problem
      character * ( 50 ) title

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST01'
      write ( *, '(a)' ) '  For each problem, print the title.'
c
c  Get the number of problems.
c
      call p00_problem_num ( problem_num )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Problem   Title'
      write ( *, '(a)' ) ' '

      do problem = 1, problem_num

        call p00_title ( problem, title )

        write ( *, '(2x,i6,2x,a)' ) problem, title

      end do

      return
      end
      subroutine test02 ( )

c*********************************************************************72
c
cc TEST02 evaluates the objective function at each starting point.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    15 February 2012
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer m_max
      parameter ( m_max = 10 )
      integer n
      parameter ( n = 100000 )

      double precision a(m_max)
      double precision b(m_max)
      double precision f(n)
      double precision f_max
      double precision f_min
      double precision fs
      integer i
      integer know
      integer m
      integer problem
      integer problem_num
      integer seed
      character * ( 50 ) title
      double precision x(m_max,n)
      double precision xs(m_max)

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST02'
      write ( *, '(a)' ) 
     &  '  For each problem, evaluate the function at many points.'
      write ( *, '(a,i8)' ) '  Number of sample points = ', n
c
c  Get the number of problems.
c
      call p00_problem_num ( problem_num )

      do problem = 1, problem_num

        write ( *, '(a)' ) ' '
        write ( *, '(a,i8)' ) '  Problem ', problem

        call p00_title ( problem, title )

        write ( *, '(2x,a)' ) trim ( title )

        call p00_m ( problem, m )

        write ( *, '(2x,a,i8)' ) '  M =     ', m
     
        call p00_ab ( problem, m, a, b )

        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) '    I      A(i)      B(i)'
        write ( *, '(a)' ) ' '

        do i = 1, m
          write ( *, '(2x,i4,2x,f10.4,2x,f10.4,2x,f10.4)' ) 
     &      i, a(i), b(i)
        end do

        seed = 123456789
        call r8col_uniform ( m, n, a, b, seed, x )
        call p00_f ( problem, m, n, x, f )
        
        call r8vec_max ( n, f, f_max )
        call r8vec_min ( n, f, f_min )

        write ( *, '(a)' ) ' '
        write ( *, '(a,g14.6)' ) '  Max(F) = ', f_max
        write ( *, '(a,g14.6)' ) '  Min(F) = ', f_min

        know = 0
        call p00_sol ( problem, m, know, xs )
        if ( know .ne. 0 ) then
          call p00_f ( problem, m, 1, xs, fs )
          write ( *, '(a,g14.6)' ) '  F(X*)  = ', fs
        else
          write ( *, '(a)' ) '  X* is not given.'
        end if

      end do

      return
      end
