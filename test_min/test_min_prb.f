      program main

c*********************************************************************72
c
cc MAIN is the main program for TEST_MIN_PRB.
c
c  Discussion:
c
c    TEST_MIN_PRB calls the TEST_MIN tests.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    03 February 2012
c
c  Author:
c
c    John Burkardt
c
      implicit none

      call timestamp (  )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST_MIN_PRB'
      write ( *, '(a)' ) '  FORTRAN77 version'
      write ( *, '(a)' ) '  Test the TEST_MIN library.'

      call test01 ( )
      call test02 ( )
      call test03 ( )
      call test04 ( )
      call test05 ( )
      call test06 ( )
c
c  Terminate.
c
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST_MIN_PRB'
      write ( *, '(a)' ) '  Normal end of execution.'

      write ( *, '(a)' ) ' '
      call timestamp ( )

      stop
      end
      subroutine test01 ( )

c*********************************************************************72
c
cc TEST01 prints the title of each problem.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    03 February 2012
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
      write ( *, '(a)' ) '   Problem Title'
      write ( *, '(a)' ) ' '

      do problem = 1, problem_num

        call p00_title ( problem, title )

        write ( *, '(2x,i8,2x,a)' ) problem, trim ( title )

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
c    03 February 2012
c
c  Author:
c
c    John Burkardt
c
      implicit none

      double precision f_sol
      double precision f_start
      integer know
      integer problem_num
      integer problem
      character * ( 50 ) title
      double precision x

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST02'
      write ( *, '(a)' ) '  For each problem, evaluate the function'
      write ( *, '(a)' ) '  at the starting point and the solution.'
c
c  Get the number of problems.
c
      call p00_problem_num ( problem_num )

      do problem = 1, problem_num

        call p00_title ( problem, title )

        write ( *, '(a)' ) ' '
        write ( *, '(a,i8)' ) '  Problem ', problem
        write ( *, '(2x,a)' ) trim ( title )
        write ( *, '(a)' ) ' '
     
        call p00_start ( problem, x )

        call p00_f ( problem, x, f_start )

        write ( *, '(4x,a,g16.8)' ) 'F(X_START)=', f_start

        call p00_sol ( problem, know, x )

        if ( 0 .lt. know ) then
          call p00_f ( problem, x, f_sol )
          write ( *, '(4x,a,g16.8)' ) 'F(X_SOL)=  ', f_sol
        end if

      end do

      return
      end
      subroutine test03 ( )

c*********************************************************************72
c
cc TEST03 compares the exact and approximate first derivatives.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    03 February 2012
c
c  Author:
c
c    John Burkardt
c
      implicit none

      double precision f1
      double precision f1_dif
      integer problem_num
      integer problem
      character * ( 50 ) title
      double precision x

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST03'
      write ( *, '(a)' ) '  For each problem, compare the exact and'
      write ( *, '(a)' ) 
     &  '  approximate gradients at the starting point.'
c
c  Get the number of problems.
c
      call p00_problem_num ( problem_num )

      do problem = 1, problem_num

        call p00_title ( problem, title )

        write ( *, '(a)' ) ' '
        write ( *, '(a,i8)' ) '  Problem ', problem
        write ( *, '(2x,a)' ) trim ( title )

        call p00_start ( problem, x )

        call p00_f1 ( problem, x, f1 )

        call p00_f1_dif ( problem, x, f1_dif )

        write ( *, '(a)' ) ' '
        write ( *, '(2x,a)' ) 'X'
        write ( *, '(4x,5g16.8)' ) x
        write ( *, '(2x,a)' ) 'F''(X) (exact)'
        write ( *, '(4x,5g16.8)' ) f1
        write ( *, '(2x,a)' ) 'F''(X) (difference)'
        write ( *, '(4x,5g16.8)' ) f1_dif

      end do

      return
      end
      subroutine test04 ( )

c*********************************************************************72
c
cc TEST04 compares the exact and approximate second derivatives.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    03 February 2012
c
c  Author:
c
c    John Burkardt
c
      implicit none

      double precision f2
      double precision f2_dif
      integer problem_num
      integer problem
      character * ( 50 ) title
      double precision x

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST04'
      write ( *, '(a)' ) '  For each problem, compare the exact and'
      write ( *, '(a)' ) 
     &  '  approximate second derivatives at the starting point.'
c
c  Get the number of problems.
c
      call p00_problem_num ( problem_num )

      do problem = 1, problem_num

        call p00_title ( problem, title )

        write ( *, '(a)' ) ' '
        write ( *, '(a,i8)' ) '  Problem ', problem
        write ( *, '(2x,a)' ) trim ( title )

        call p00_start ( problem, x )

        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) '  X:'
        write ( *, '(4x,5g16.8)' ) x

        call p00_f2 ( problem, x, f2 )

        write ( *, '(a)' ) '  F"(X) (exact):'
        write ( *, '(4x,6g13.5)' ) f2

        call p00_f2_dif ( problem, x, f2_dif )

        write ( *, '(a)' ) '  F"(X) (difference):'
        write ( *, '(4x,6g13.5)' ) f2_dif

      end do

      return
      end
      subroutine test05 ( )

c*********************************************************************72
c
cc TEST05 carries out a simple bisection method.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    03 February 2012
c
c  Author:
c
c    John Burkardt
c
      implicit none

      double precision fa
      double precision fb
      double precision fc
      double precision fd
      double precision fe
      integer i
      integer max_step
      parameter ( max_step = 10 )
      integer problem_num
      integer problem
      character * ( 50 ) title
      double precision xa
      double precision xb
      double precision xc
      double precision xd
      double precision xe

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST05'
      write ( *, '(a)' ) '  For each problem, take a few steps of '
      write ( *, '(a)' ) '  the bisection method.'
c
c  Get the number of problems.
c
      call p00_problem_num ( problem_num )

      do problem = 1, problem_num

        call p00_title ( problem, title )

        write ( *, '(a)' ) ' '
        write ( *, '(a,i8)' ) '  Problem ', problem
        write ( *, '(2x,a)' ) trim ( title )

        call p00_interval ( problem, xa, xc )
        xb = 0.5D+00 * ( xa + xc )
        call p00_f ( problem, xa, fa )
        call p00_f ( problem, xc, fc )
        call p00_f ( problem, xb, fb )

        i = 0
        write ( *, '(a)' ) ' '
        write ( *, '(i6)' ) i
        write ( *, '(a,3g16.8)' ) '  X:', xa, xb, xc
        write ( *, '(a,3g16.8)' ) '  F:', fa, fb, fc

        do i = 1, max_step

          xd = 0.5D+00 * ( xa + xb )
          call p00_f ( problem, xd, fd )

          xe = 0.5D+00 * ( xb + xc )
          call p00_f ( problem, xe, fe )

          if ( fd .le. fb ) then
            xc = xb
            fc = fb
            xb = xd
            fb = fd
          else if ( fe .le. fb ) then
            xa = xb
            fa = fb
            xb = xe
            fb = fe
          else
            xa = xd
            fa = fd
            xc = xe
            fc = fe
          end if

          write ( *, '(i6)' ) i
          write ( *, '(a,3g16.8)' ) '  X:', xa, xb, xc
          write ( *, '(a,3g16.8)' ) '  F:', fa, fb, fc

        end do

      end do

      return
      end
      subroutine test06 ( )

c*********************************************************************72
c
cc TEST06 carries out a version of Brent's derivative-free minimizer.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    03 February 2012
c
c  Author:
c
c    John Burkardt
c
      implicit none

      double precision fa
      double precision fb
      double precision fx
      double precision p00_fmin
      integer problem_num
      integer problem
      character * ( 50 ) title
      double precision tol
      parameter ( tol = 0.000001D+00 )
      double precision x
      double precision xa
      double precision xb

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST06'
      write ( *, '(a)' ) '  For each problem, use Brent''s method.'
c
c  Get the number of problems.
c
      call p00_problem_num ( problem_num )

      do problem = 1, problem_num

        call p00_title ( problem, title )

        write ( *, '(a)' ) ' '
        write ( *, '(a,i8)' ) '  Problem ', problem
        write ( *, '(2x,a)' ) trim ( title )

        call p00_interval ( problem, xa, xb )

        call p00_f ( problem, xa, fa )
        call p00_f ( problem, xb, fb )

        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) '  Initial interval [A,B]:'
        write ( *, '(a)' ) ' '
        write ( *, '(a,g16.8,14x,g16.8)' ) 
     &    '   A,       B:', xa,     xb
        write ( *, '(a,g16.8,14x,g16.8)' ) 
     &    '  FA,      FB:', fa,     fb

        x = p00_fmin ( xa, xb, problem, tol )

        call p00_f ( problem, xa, fa )
        call p00_f ( problem, xb, fb )
        call p00_f ( problem, x, fx )

        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) '  Final interval [A,X*,B]:'
        write ( *, '(a)' ) ' '
        write ( *, '(a,g16.8,g16.8,g16.8)' ) 
     &    '   A,  X*,  B:', xa, x,  xb
        write ( *, '(a,g16.8,g16.8,g16.8)' ) 
     &    '  FA, FX*, FB:', fa, fx, fb

      end do

      return
      end
