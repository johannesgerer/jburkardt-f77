      subroutine p00_fun ( problem, option, nvar, x, fx )

c*********************************************************************72
c
cc P00_FUN evaluates the function for any problem.
c
c  Discussion:
c
c    These problems were collected by Professor Werner Rheinboldt, of the
c    University of Pittsburgh, and were used in the development of the
c    PITCON program.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    19 September 2008
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer PROBLEM, the problem index.
c
c    Input, integer OPTION, the option index.
c
c    Input, integer NVAR, the number of variables.
c
c    Input, double precision X(NVAR), the argument of the function.
c
c    Output, double precision FX(NVAR-1), the value of the function at X.
c
      implicit none

      integer nvar

      double precision fx(nvar-1)
      integer option
      integer problem
      double precision x(nvar)

      if ( problem .eq. 1 ) then
        call p01_fun ( option, nvar, x, fx )
      else if ( problem .eq. 2 ) then
        call p02_fun ( option, nvar, x, fx )
      else if ( problem .eq. 3 ) then
        call p03_fun ( option, nvar, x, fx )
      else if ( problem .eq. 4 ) then
        call p04_fun ( option, nvar, x, fx )
      else if ( problem .eq. 5 ) then
        call p05_fun ( option, nvar, x, fx )
      else if ( problem .eq. 6 ) then
        call p06_fun ( option, nvar, x, fx )
      else if ( problem .eq. 7 ) then
        call p07_fun ( option, nvar, x, fx )
      else if ( problem .eq. 8 ) then
        call p08_fun ( option, nvar, x, fx )
      else if ( problem .eq. 9 ) then
        call p09_fun ( option, nvar, x, fx )
      else if ( problem .eq. 10 ) then
        call p10_fun ( option, nvar, x, fx )
      else if ( problem .eq. 11 ) then
        call p11_fun ( option, nvar, x, fx )
      else if ( problem .eq. 12 ) then
        call p12_fun ( option, nvar, x, fx )
      else if ( problem .eq. 13 ) then
        call p13_fun ( option, nvar, x, fx )
      else if ( problem .eq. 14 ) then
        call p14_fun ( option, nvar, x, fx )
      else if ( problem .eq. 15 ) then
        call p15_fun ( option, nvar, x, fx )
      else if ( problem .eq. 16 ) then
        call p16_fun ( option, nvar, x, fx )
      else if ( problem .eq. 17 ) then
        call p17_fun ( option, nvar, x, fx )
      else if ( problem .eq. 18 ) then
        call p18_fun ( option, nvar, x, fx )
      else if ( problem .eq. 19 ) then
        call p19_fun ( option, nvar, x, fx )
      else if ( problem .eq. 20 ) then
        call p20_fun ( option, nvar, x, fx )
      else if ( problem .eq. 21 ) then
        call p21_fun ( option, nvar, x, fx )
      else if ( problem .eq. 22 ) then
        call p22_fun ( option, nvar, x, fx )
      else if ( problem .eq. 23 ) then
        call p23_fun ( option, nvar, x, fx )
      else if ( problem .eq. 24 ) then
        call p24_fun ( option, nvar, x, fx )
      else if ( problem .eq. 25 ) then
        call p25_fun ( option, nvar, x, fx )
      else if ( problem .eq. 26 ) then
        call p26_fun ( option, nvar, x, fx )
      else if ( problem .eq. 27 ) then
        call p27_fun ( option, nvar, x, fx )
      else if ( problem .eq. 28 ) then
        call p28_fun ( option, nvar, x, fx )
      else if ( problem .eq. 29 ) then
        call p29_fun ( option, nvar, x, fx )
      else
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'P00_FUN - Fatal error!'
        write ( *, '(a,i8)' )
     &    '  Unrecognized problem number = ', problem
        stop
      end if

      return
      end
      subroutine p00_jac ( problem, option, nvar, x, jac )

c*********************************************************************72
c
cc P00_JAC evaluates the jacobian for any problem.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    24 September 2008
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer PROBLEM, the problem index.
c
c    Input, integer OPTION, the option index.
c
c    Input, integer NVAR, the number of variables.
c
c    Input, double precision X(NVAR), the argument of the jacobian.
c
c    Output, double precision JAC(NVAR,NVAR), the jacobian matrix evaluated
c    at X.  The NVAR-th row is not set by this routine.
c
      implicit none

      integer nvar

      integer option
      double precision jac(nvar,nvar)
      integer problem
      double precision x(nvar)

      if ( problem .eq. 1 ) then
        call p01_jac ( option, nvar, x, jac )
      else if ( problem .eq. 2 ) then
        call p02_jac ( option, nvar, x, jac )
      else if ( problem .eq. 3 ) then
        call p03_jac ( option, nvar, x, jac )
      else if ( problem .eq. 4 ) then
        call p04_jac ( option, nvar, x, jac )
      else if ( problem .eq. 5 ) then
        call p05_jac ( option, nvar, x, jac )
      else if ( problem .eq. 6 ) then
        call p06_jac ( option, nvar, x, jac )
      else if ( problem .eq. 7 ) then
        call p07_jac ( option, nvar, x, jac )
      else if ( problem .eq. 8 ) then
        call p08_jac ( option, nvar, x, jac )
      else if ( problem .eq. 9 ) then
        call p09_jac ( option, nvar, x, jac )
      else if ( problem .eq. 10 ) then
        call p10_jac ( option, nvar, x, jac )
      else if ( problem .eq. 11 ) then
        call p11_jac ( option, nvar, x, jac )
      else if ( problem .eq. 12 ) then
        call p12_jac ( option, nvar, x, jac )
      else if ( problem .eq. 13 ) then
        call p13_jac ( option, nvar, x, jac )
      else if ( problem .eq. 14 ) then
        call p14_jac ( option, nvar, x, jac )
      else if ( problem .eq. 15 ) then
        call p15_jac ( option, nvar, x, jac )
      else if ( problem .eq. 16 ) then
        call p16_jac ( option, nvar, x, jac )
      else if ( problem .eq. 17 ) then
        call p17_jac ( option, nvar, x, jac )
      else if ( problem .eq. 18 ) then
        call p18_jac ( option, nvar, x, jac )
      else if ( problem .eq. 19 ) then
        call p19_jac ( option, nvar, x, jac )
      else if ( problem .eq. 20 ) then
        call p20_jac ( option, nvar, x, jac )
      else if ( problem .eq. 21 ) then
        call p21_jac ( option, nvar, x, jac )
      else if ( problem .eq. 22 ) then
        call p22_jac ( option, nvar, x, jac )
      else if ( problem .eq. 23 ) then
        call p23_jac ( option, nvar, x, jac )
      else if ( problem .eq. 24 ) then
        call p24_jac ( option, nvar, x, jac )
      else if ( problem .eq. 25 ) then
        call p25_jac ( option, nvar, x, jac )
      else if ( problem .eq. 26 ) then
        call p26_jac ( option, nvar, x, jac )
      else if ( problem .eq. 27 ) then
        call p27_jac ( option, nvar, x, jac )
      else if ( problem .eq. 28 ) then
        call p28_jac ( option, nvar, x, jac )
      else if ( problem .eq. 29 ) then
        call p29_jac ( option, nvar, x, jac )
      else
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'P00_JAC - Fatal errorc'
        write ( *, '(a,i8)' )
     &    '  Unrecognized problem number = ', problem
        stop
      end if

      return
      end
      subroutine p00_nvar ( problem, option, nvar )

c*********************************************************************72
c
cc P00_NVAR returns the number of variables for any problem.
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
c    Input, integer PROBLEM, the problem index.
c
c    Input, integer OPTION, the option index.
c
c    Output, integer NVAR, the number of variables.
c
      implicit none

      integer nvar
      integer option
      integer problem

      if (  problem .eq. 1 ) then
        call p01_nvar ( option, nvar )
      else if ( problem .eq. 2 ) then
        call p02_nvar ( option, nvar )
      else if ( problem .eq. 3 ) then
        call p03_nvar ( option, nvar )
      else if ( problem .eq. 4 ) then
        call p04_nvar ( option, nvar )
      else if ( problem .eq. 5 ) then
        call p05_nvar ( option, nvar )
      else if ( problem .eq. 6 ) then
        call p06_nvar ( option, nvar )
      else if ( problem .eq. 7 ) then
        call p07_nvar ( option, nvar )
      else if ( problem .eq. 8 ) then
        call p08_nvar ( option, nvar )
      else if ( problem .eq. 9 ) then
        call p09_nvar ( option, nvar )
      else if ( problem .eq. 10 ) then
        call p10_nvar ( option, nvar )
      else if ( problem .eq. 11 ) then
        call p11_nvar ( option, nvar )
      else if ( problem .eq. 12 ) then
        call p12_nvar ( option, nvar )
      else if ( problem .eq. 13 ) then
        call p13_nvar ( option, nvar )
      else if ( problem .eq. 14 ) then
        call p14_nvar ( option, nvar )
      else if ( problem .eq. 15 ) then
        call p15_nvar ( option, nvar )
      else if ( problem .eq. 16 ) then
        call p16_nvar ( option, nvar )
      else if ( problem .eq. 17 ) then
        call p17_nvar ( option, nvar )
      else if ( problem .eq. 18 ) then
        call p18_nvar ( option, nvar )
      else if ( problem .eq. 19 ) then
        call p19_nvar ( option, nvar )
      else if ( problem .eq. 20 ) then
        call p20_nvar ( option, nvar )
      else if ( problem .eq. 21 ) then
        call p21_nvar ( option, nvar )
      else if ( problem .eq. 22 ) then
        call p22_nvar ( option, nvar )
      else if ( problem .eq. 23 ) then
        call p23_nvar ( option, nvar )
      else if ( problem .eq. 24 ) then
        call p24_nvar ( option, nvar )
      else if ( problem .eq. 25 ) then
        call p25_nvar ( option, nvar )
      else if ( problem .eq. 26 ) then
        call p26_nvar ( option, nvar )
      else if ( problem .eq. 27 ) then
        call p27_nvar ( option, nvar )
      else if ( problem .eq. 28 ) then
        call p28_nvar ( option, nvar )
      else if ( problem .eq. 29 ) then
        call p29_nvar ( option, nvar )
      else
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'P00_NVAR - Fatal error!'
        write ( *, '(a,i8)' )
     &  '  Unrecognized problem index  = ', problem
        stop
      end if

      return
      end
      subroutine p00_option_num ( problem, option_num )

c*********************************************************************72
c
cc P00_OPTION_NUM returns the number of options available for a problem.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    16 September 2008
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer PROBLEM, the problem index.
c
c    Output, integer OPTION_NUM, the number of options available for
c    this problem.  OPTION_NUM is always at least 1.
c
      implicit none

      integer option_num
      integer problem

      if (  problem .eq. 1 ) then
        call p01_option_num ( option_num )
      else if ( problem .eq. 2 ) then
        call p02_option_num ( option_num )
      else if ( problem .eq. 3 ) then
        call p03_option_num ( option_num )
      else if ( problem .eq. 4 ) then
        call p04_option_num ( option_num )
      else if ( problem .eq. 5 ) then
        call p05_option_num ( option_num )
      else if ( problem .eq. 6 ) then
        call p06_option_num ( option_num )
      else if ( problem .eq. 7 ) then
        call p07_option_num ( option_num )
      else if ( problem .eq. 8 ) then
        call p08_option_num ( option_num )
      else if ( problem .eq. 9 ) then
        call p09_option_num ( option_num )
      else if ( problem .eq. 10 ) then
        call p10_option_num ( option_num )
      else if ( problem .eq. 11 ) then
        call p11_option_num ( option_num )
      else if ( problem .eq. 12 ) then
        call p12_option_num ( option_num )
      else if ( problem .eq. 13 ) then
        call p13_option_num ( option_num )
      else if ( problem .eq. 14 ) then
        call p14_option_num ( option_num )
      else if ( problem .eq. 15 ) then
        call p15_option_num ( option_num )
      else if ( problem .eq. 16 ) then
        call p16_option_num ( option_num )
      else if ( problem .eq. 17 ) then
        call p17_option_num ( option_num )
      else if ( problem .eq. 18 ) then
        call p18_option_num ( option_num )
      else if ( problem .eq. 19 ) then
        call p19_option_num ( option_num )
      else if ( problem .eq. 20 ) then
        call p20_option_num ( option_num )
      else if ( problem .eq. 21 ) then
        call p21_option_num ( option_num )
      else if ( problem .eq. 22 ) then
        call p22_option_num ( option_num )
      else if ( problem .eq. 23 ) then
        call p23_option_num ( option_num )
      else if ( problem .eq. 24 ) then
        call p24_option_num ( option_num )
      else if ( problem .eq. 25 ) then
        call p25_option_num ( option_num )
      else if ( problem .eq. 26 ) then
        call p26_option_num ( option_num )
      else if ( problem .eq. 27 ) then
        call p27_option_num ( option_num )
      else if ( problem .eq. 28 ) then
        call p28_option_num ( option_num )
      else if ( problem .eq. 29 ) then
        call p29_option_num ( option_num )
      else
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'P00_OPTION_NUM - Fatal error!'
        write ( *, '(a,i8)' )
     &  '  Unrecognized problem index  = ', problem
        stop
      end if

      return
      end
      subroutine p00_problem_num ( problem_num )

c*********************************************************************72
c
cc P00_PROBLEM_NUM returns the number of problems available.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    16 September 2008
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, integer PROBLEM_NUM, the number of problems.
c
      implicit none

      integer problem_num

      problem_num = 29

      return
      end
      subroutine p00_start ( problem, option, nvar, x )

c*********************************************************************72
c
cc P00_START returns a starting point for any problem.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c   19 September 2008
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer PROBLEM, the problem index.
c
c    Input, integer OPTION, the option index.
c
c    Input, integer NVAR, the number of variables.
c
c    Output, double precision X(NVAR), the starting point.
c
      implicit none

      integer nvar

      integer option
      integer problem
      double precision x(nvar)

      if ( problem .eq. 1 ) then
        call p01_start ( option, nvar, x )
      else if ( problem .eq. 2 ) then
        call p02_start ( option, nvar, x )
      else if ( problem .eq. 3 ) then
        call p03_start ( option, nvar, x )
      else if ( problem .eq. 4 ) then
        call p04_start ( option, nvar, x )
      else if ( problem .eq. 5 ) then
        call p05_start ( option, nvar, x )
      else if ( problem .eq. 6 ) then
        call p06_start ( option, nvar, x )
      else if ( problem .eq. 7 ) then
        call p07_start ( option, nvar, x )
      else if ( problem .eq. 8 ) then
        call p08_start ( option, nvar, x )
      else if ( problem .eq. 9 ) then
        call p09_start ( option, nvar, x )
      else if ( problem .eq. 10 ) then
        call p10_start ( option, nvar, x )
      else if ( problem .eq. 11 ) then
        call p11_start ( option, nvar, x )
      else if ( problem .eq. 12 ) then
        call p12_start ( option, nvar, x )
      else if ( problem .eq. 13 ) then
        call p13_start ( option, nvar, x )
      else if ( problem .eq. 14 ) then
        call p14_start ( option, nvar, x )
      else if ( problem .eq. 15 ) then
        call p15_start ( option, nvar, x )
      else if ( problem .eq. 16 ) then
        call p16_start ( option, nvar, x )
      else if ( problem .eq. 17 ) then
        call p17_start ( option, nvar, x )
      else if ( problem .eq. 18 ) then
        call p18_start ( option, nvar, x )
      else if ( problem .eq. 19 ) then
        call p19_start ( option, nvar, x )
      else if ( problem .eq. 20 ) then
        call p20_start ( option, nvar, x )
      else if ( problem .eq. 21 ) then
        call p21_start ( option, nvar, x )
      else if ( problem .eq. 22 ) then
        call p22_start ( option, nvar, x )
      else if ( problem .eq. 23 ) then
        call p23_start ( option, nvar, x )
      else if ( problem .eq. 24 ) then
        call p24_start ( option, nvar, x )
      else if ( problem .eq. 25 ) then
        call p25_start ( option, nvar, x )
      else if ( problem .eq. 26 ) then
        call p26_start ( option, nvar, x )
      else if ( problem .eq. 27 ) then
        call p27_start ( option, nvar, x )
      else if ( problem .eq. 28 ) then
        call p28_start ( option, nvar, x )
      else if ( problem .eq. 29 ) then
        call p29_start ( option, nvar, x )
      else
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'P00_START - Fatal error!'
        write ( *, '(a,i8)' )
     &    '  Unrecognized problem index = ', problem
        stop
      end if

      return
      end
      subroutine p00_stepsize ( problem, option, h, hmin, hmax )

c*********************************************************************72
c
cc P00_STEPSIZE returns step sizes for any problem.
c
c  Discussion:
c
c    The routine returns a suggested initial stepsize, and suggestions for
c    the minimum and maximum stepsizes.
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
c    Input, integer PROBLEM, the problem index.
c
c    Input, integer OPTION, the option index.
c
c    Output, double precision H, HMIN, HMAX, suggested values for the
c    initial step, the minimum step, and the maximum step.
c
      implicit none

      double precision h
      double precision hmax
      double precision hmin
      integer option
      integer problem

      if ( problem == 1 ) then
        call p01_stepsize ( option, h, hmin, hmax )
      else if ( problem == 2 ) then
        call p02_stepsize ( option, h, hmin, hmax )
      else if ( problem == 3 ) then
        call p03_stepsize ( option, h, hmin, hmax )
      else if ( problem == 4 ) then
        call p04_stepsize ( option, h, hmin, hmax )
      else if ( problem == 5 ) then
        call p05_stepsize ( option, h, hmin, hmax )
      else if ( problem == 6 ) then
        call p06_stepsize ( option, h, hmin, hmax )
      else if ( problem == 7 ) then
        call p07_stepsize ( option, h, hmin, hmax )
      else if ( problem == 8 ) then
        call p08_stepsize ( option, h, hmin, hmax )
      else if ( problem == 9 ) then
        call p09_stepsize ( option, h, hmin, hmax )
      else if ( problem == 10 ) then
        call p10_stepsize ( option, h, hmin, hmax )
      else if ( problem == 11 ) then
        call p11_stepsize ( option, h, hmin, hmax )
      else if ( problem == 12 ) then
        call p12_stepsize ( option, h, hmin, hmax )
      else if ( problem == 13 ) then
        call p13_stepsize ( option, h, hmin, hmax )
      else if ( problem == 14 ) then
        call p14_stepsize ( option, h, hmin, hmax )
      else if ( problem == 15 ) then
        call p15_stepsize ( option, h, hmin, hmax )
      else if ( problem == 16 ) then
        call p16_stepsize ( option, h, hmin, hmax )
      else if ( problem == 17 ) then
        call p17_stepsize ( option, h, hmin, hmax )
      else if ( problem == 18 ) then
        call p18_stepsize ( option, h, hmin, hmax )
      else if ( problem == 19 ) then
        call p19_stepsize ( option, h, hmin, hmax )
      else if ( problem == 20 ) then
        call p20_stepsize ( option, h, hmin, hmax )
      else if ( problem == 21 ) then
        call p21_stepsize ( option, h, hmin, hmax )
      else if ( problem == 22 ) then
        call p22_stepsize ( option, h, hmin, hmax )
      else if ( problem == 23 ) then
        call p23_stepsize ( option, h, hmin, hmax )
      else if ( problem == 24 ) then
        call p24_stepsize ( option, h, hmin, hmax )
      else if ( problem == 25 ) then
        call p25_stepsize ( option, h, hmin, hmax )
      else if ( problem == 26 ) then
        call p26_stepsize ( option, h, hmin, hmax )
      else if ( problem == 27 ) then
        call p27_stepsize ( option, h, hmin, hmax )
      else if ( problem == 28 ) then
        call p28_stepsize ( option, h, hmin, hmax )
      else if ( problem == 29 ) then
        call p29_stepsize ( option, h, hmin, hmax )
      else
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'P00_STEPSIZE - Fatal error!'
        write ( *, '(a,i8)' )
     &    '  Unrecognized problem number = ', problem
        stop
      end if

      return
      end
      subroutine p00_title ( problem, option, title )

c*********************************************************************72
c
cc P00_TITLE sets the title for any problem.
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
c    Input, integer PROBLEM, the problem index.
c
c    Input, integer OPTION, the option index.
c
c    Output, character*(*) TITLE, the title of the problem.
c    TITLE will never be longer than 80 characters.
c
      implicit none

      integer option
      integer problem
      character*(*) title

      if (  problem .eq. 1 ) then
        call p01_title ( option, title )
      else if ( problem .eq. 2 ) then
        call p02_title ( option, title )
      else if ( problem .eq. 3 ) then
        call p03_title ( option, title )
      else if ( problem .eq. 4 ) then
        call p04_title ( option, title )
      else if ( problem .eq. 5 ) then
        call p05_title ( option, title )
      else if ( problem .eq. 6 ) then
        call p06_title ( option, title )
      else if ( problem .eq. 7 ) then
        call p07_title ( option, title )
      else if ( problem .eq. 8 ) then
        call p08_title ( option, title )
      else if ( problem .eq. 9 ) then
        call p09_title ( option, title )
      else if ( problem .eq. 10 ) then
        call p10_title ( option, title )
      else if ( problem .eq. 11 ) then
        call p11_title ( option, title )
      else if ( problem .eq. 12 ) then
        call p12_title ( option, title )
      else if ( problem .eq. 13 ) then
        call p13_title ( option, title )
      else if ( problem .eq. 14 ) then
        call p14_title ( option, title )
      else if ( problem .eq. 15 ) then
        call p15_title ( option, title )
      else if ( problem .eq. 16 ) then
        call p16_title ( option, title )
      else if ( problem .eq. 17 ) then
        call p17_title ( option, title )
      else if ( problem .eq. 18 ) then
        call p18_title ( option, title )
      else if ( problem .eq. 19 ) then
        call p19_title ( option, title )
      else if ( problem .eq. 20 ) then
        call p20_title ( option, title )
      else if ( problem .eq. 21 ) then
        call p21_title ( option, title )
      else if ( problem .eq. 22 ) then
        call p22_title ( option, title )
      else if ( problem .eq. 23 ) then
        call p23_title ( option, title )
      else if ( problem .eq. 24 ) then
        call p24_title ( option, title )
      else if ( problem .eq. 25 ) then
        call p25_title ( option, title )
      else if ( problem .eq. 26 ) then
        call p26_title ( option, title )
      else if ( problem .eq. 27 ) then
        call p27_title ( option, title )
      else if ( problem .eq. 28 ) then
        call p28_title ( option, title )
      else if ( problem .eq. 29 ) then
        call p29_title ( option, title )
      else
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'P00_TITLE - Fatal error!'
        write ( *, '(a,i8)' ) '  Unrecognized problem index = ', problem
        stop
      end if

      return
      end
      subroutine p01_fun ( option, nvar, x, fx )

c*********************************************************************72
c
cc P01_FUN evaluates the Freudenstein-Roth function.
c
c  Title:
c
c    The Freudenstein-Roth function
c
c  Description:
c
c    One way to use a continuation code as a nonlinear root finder
c    is to start with a set of nonlinear equations G(X), and an
c    approximate root A, and create a "homotopy" function F(X,Y)
c    with the properties that F(A,0.0) = 0 and F(X,1.0) = G(X).
c    Thus, the homotopy function F has a known exact solution
c    from which we can start with no difficulty.  If the continuation
c    code can take us from Y = 0 to Y = 1, then we have found
c    an X so that F(X,1.0) = 0, so we have found a solution to G(X)=0.
c
c    The Freudenstein-Roth function F(X) is derived in this way
c    from a homotopy of G(X):
c
c      F ( X(1), X(2), X(3) ) =
c        G ( X(1), X(2) ) - ( 1 - X(3) ) * G ( Y1, Y2 )
c
c    where Y1 and Y2 are some fixed values, and
c
c      G(1) = X(1) - X(2)*X(2)*X(2) + 5*X(2)*X(2) -  2*X(2) - 13
c      G(2) = X(1) + X(2)*X(2)*X(2) +   X(2)*X(2) - 14*X(2) - 29
c
c  Options 1, 2, 3:
c
c    The starting point is X0 = ( 15, -2, 0 ).
c
c    A great deal of information is available about the homotopy curve
c    generated by this starting point:
c
c    The function F(X) has the form
c
c      F(1) = X(1) - X(2)**3 + 5*X(2)**2 -  2*X(2) - 13 + 34*(X(3)-1)
c      F(2) = X(1) + X(2)**3 +   X(2)**2 - 14*X(2) - 29 + 10*(X(3)-1)
c
c    There is a closed form representation of the curve in terms of the
c    second parameter:
c
c      X(1) = (-11*X(2)**3 + 4*X(2)**2 + 114*X(2) + 214) /  6
c      X(2) = X(2)
c      X(3) = (    X(2)**3 - 2*X(2)**2 -   6*X(2) +   4) / 12
c
c    The first option simply requests the production of solution points
c    along the curve until a point is reached whose third component is
c    exactly 1.
c
c    Options 2 and 3 use the same starting point, and also stop when the
c    third component is 1.  However, these options in addition search
c    for limit points in the first and third components of the solution,
c    respectively.
c
c    The target solution has X(3) = 1, and is ( 5, 4, 1 ).
c
c    Limit points for X1:
c
c      ( 14.28309, -1.741377,  0.2585779 )
c      ( 61.66936,  1.983801, -0.6638797 )
c
c    Limit points for X3:
c
c     (20.48586, -0.8968053, 0.5875873)
c     (61.02031,  2.230139, -0.6863528)
c
c    The curve has several dramatic bends.
c
c
c  Options 4, 5, and 6:
c
c    The starting point is (4, 3, 0).
c
c    The function F(X) has the form
c
c      F(1) = X(1) - X(2)**3 + 5*X(2)**2 -  2*X(2) - 13 +  3*(X(3)-1)
c      F(2) = X(1) + X(2)**3 +   X(2)**2 - 14*X(2) - 29 - 31*(X(3)-1)
c
c    There is a closed form representation of the curve in terms of the
c    second parameter:
c
c      X(1) = (14*X(2)**3 -79*X(2)**2 +52*X(2) + 245) / 17
c      X(2) = X(2)
c      X(3) = (   X(2)**3 - 2*X(2)**2 - 6*X(2) +   9) / 17
c
c    The correct value of the solution at X(3)=1 is:
c
c      (5, 4, 1)
c
c    In option 5, limit points in the first component are sought,
c    and in option 6, limit points in the third component are
c    sought.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    29 September 2008
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Ferdinand Freudenstein, Bernhard Roth,
c    Numerical Solutions of Nonlinear Equations,
c    Journal of the Association for Computing Machinery,
c    Volume 10, 1963, Pages 550-556.
c
c  Parameters:
c
c    Input, integer OPTION, the option index.
c
c    Input, integer NVAR, the number of variables.
c
c    Input, double precision X(NVAR), the argument of the function.
c
c    Output, double precision FX(NVAR-1), the value of the function at X.
c
      integer nvar

      double precision fx(nvar-1)
      double precision gx(2)
      double precision gy(2)
      integer i
      integer option
      double precision x(nvar)
      double precision y(3)
c
c  Get the starting point.
c
      call p01_start ( option, nvar, y )
c
c  G is the function value at the starting point,
c  F the function value at the current point.
c
      call p01_gx ( y, gy )

      call p01_gx ( x, gx )
c
c  The parameter X3 generates the homotopy curve.
c
      do i = 1, nvar - 1
        fx(i) = gx(i) + ( x(3) - 1.0D+00 ) * gy(i)
      end do

      return
      end
      subroutine p01_gx (neqn,x,g)

c*********************************************************************72
c
cc P01_GX is an auxilliary routine for the Freudenstein-Roth function.
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
      implicit none

      integer neqn

      double precision g(2)
      double precision x(2)

      g(1) = x(1) - ( ( x(2) - 5.0D+00 ) * x(2) + 2.0D+00 ) * x(2)
     &  - 13.0D+00

      g(2) = x(1) + ( ( x(2) + 1.0D+00 ) * x(2) - 14.0D+00 ) * x(2)
     &  - 29.0D+00

      return
      end
      subroutine p01_jac ( option, nvar, x, jac )

c*********************************************************************72
c
cc P01_JAC evaluates the Freudenstein-Roth jacobian.
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
c    Input, integer OPTION, the option index.
c
c    Input, integer NVAR, the number of variables.
c
c    Input, double precision X(NVAR), the argument of the jacobian.
c
c    Output, double precision JAC(NVAR,NVAR), the jacobian matrix evaluated
c    at X.  The NVAR-th row is not set by this routine.
c
      implicit none

      integer nvar

      double precision gy(3)
      integer option
      double precision jac(nvar,nvar)
      double precision x(nvar)
      double precision y(3)

      jac(1,1) = 1.0D+00
      jac(2,1) = 1.0D+00

      jac(1,2) = ( - 3.0D+00 * x(2) + 10.0D+00 ) * x(2) - 2.0D+00
      jac(2,2) = ( 3.0D+00 * x(2) + 2.0D+00 ) * x(2) - 14.0D+00
c
c  Get the starting point
c
      call p01_start ( option, nvar, y )
c
c  Get the function value at the starting point
c
      call p01_gx ( y, gy )

      jac(1,3) = gy(1)
      jac(2,3) = gy(2)

      return
      end
      subroutine p01_nvar ( option, nvar )

c*********************************************************************72
c
cc P01_NVAR sets the number of variables for problem 1.
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
c    Input, integer OPTION, the option chosen for this problem.
c    For some problems, several options are available.  At least,
c    OPTION = 1 is always legal.
c
c    Output, integer NVAR, the number of variables.
c
      implicit none

      integer nvar
      integer option

      nvar = 3

      return
      end
      subroutine p01_option_num ( option_num )

c*********************************************************************72
c
cc P01_OPTION_NUM returns the number of options for problem 1.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    16 September 2008
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, integer OPTION_NUM, the number of options.
c
      implicit none

      integer option_num

      option_num = 6

      return
      end
      subroutine p01_start ( option, nvar, x )

c*********************************************************************72
c
cc P01_START returns a starting point for problem 1.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    19 September 2008
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer OPTION, the option index.
c
c    Input, integer NVAR, the number of variables.
c
c    Output, double precision X(NVAR), the starting point.
c
      implicit none

      integer nvar

      integer option
      double precision x(nvar)

      if ( option == 1 .or. option == 2 .or. option == 3 ) then

        x(1) = 15.0D+00
        x(2) = -2.0D+00
        x(3) =  0.0D+00

      else if ( option == 4 .or. option == 5 .or. option == 6 ) then

        x(1) =  4.0D+00
        x(2) =  3.0D+00
        x(3) =  0.0D+00

      else

        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'P01_START - Fatal error!'
        write ( *, '(a,i8)' ) '  Unrecognized option  = ', option
        stop

      end if

      return
      end
      subroutine p01_stepsize ( option, h, hmin, hmax )

c*********************************************************************72
c
cc P01_STEPSIZE returns step sizes for problem 1.
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
c    Input, integer OPTION, the option index.
c
c    Output, double precision H, HMIN, HMAX, suggested values for the
c    initial step, the minimum step, and the maximum step.
c
      implicit none

      double precision h
      double precision hmax
      double precision hmin
      integer option

      h =    0.30000D+00
      hmin = 0.03125D+00
      hmax = 4.00000D+00

      return
      end
      subroutine p01_title ( option, title )

c*********************************************************************72
c
cc P01_TITLE sets the title for problem 1.
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
c    Input, integer OPTION, the option index.
c
c    Output, character*(*) TITLE, the title of the problem.
c    TITLE will never be longer than 80 characters.
c
      implicit none

      integer option
      character*(*) title

      if ( option .eq. 1 ) then
        title = 'Freudenstein-Roth function, (15,-2,0).'
      else if ( option .eq. 2 ) then
        title = 'Freudenstein-Roth function, (15,-2,0), x1 limits.'
      else if ( option .eq. 3 ) then
        title = 'Freudenstein-Roth function, (15,-2,0), x3 limits.'
      else if ( option .eq. 4 ) then
        title = 'Freudenstein-Roth function, (4,3,0).'
      else if ( option .eq. 5 ) then
        title = 'Freudenstein-Roth function, (4,3,0), x1 limits.'
      else if ( option .eq. 6 ) then
        title = 'Freudenstein-Roth function, (4,3,0), x3 limits.'
      else
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'P01_TITLE - Fatal error!'
        write ( *, '(a)' ) '  Unrecognized option number.'
        stop
      end if

      return
      end
      subroutine p02_fun ( option, nvar, x, fx )

c*********************************************************************72
c
cc P02_FUN evaluates the function for problem 2.
c
c  Title:
c
c    The Boggs function
c
c  Description:
c
c    The function F is derived via homotopy from a simpler function:
c
c      F(X(1),X(2),X(3)) = G(X(1),X(2)) + (X(3)-1) * G(Y1,Y2)
c
c    with
c
c      (Y1, Y2) some starting value,
c
c    and
c
c      G(1) = X(1)*X(1) - X(2) + 1
c      G(2) = X(1) - COS(PI*X(2)/2)
c
c  Options:
c
c    OPTION = 1,
c      use starting point (  1,  0, 0 ).
c    OPTION = 2,
c      use starting point (  1, -1, 0 ).
c    OPTION = 3,
c      use starting point ( 10, 10, 0 ).
c
c  Target Points:
c
c    For the target value X(3) = 1.0, the solution is ( 0, 1, 1 ).
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    29 September 2008
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Paul Boggs,
c    The Solution of Nonlinear Systems by A-stable Integration Techniques,
c    SIAM Journal on Numerical Analysis,
c    Volume 8, Number 4, December 1971, pages 767-785.
c
c  Parameters:
c
c    Input, integer OPTION, the option index.
c
c    Input, integer NVAR, the number of variables.
c
c    Input, double precision X(NVAR), the argument of the function.
c
c    Output, double precision FX(NVAR-1), the value of the function at X.
c
      implicit none

      integer nvar

      double precision fx(nvar-1)
      double precision gx(2)
      double precision gy(2)
      integer i
      integer option
      double precision x(nvar)
      double precision y(3)
c
c  Get the starting point
c
      call p02_start ( option, nvar, y )
c
c  Get the function value at the starting point and at the
c  current point.
c
      call p02_gx ( y, gy )
      call p02_gx ( x, gx )
c
c  Use X3 to compute a homotopy.
c
      do i = 1, nvar - 1
        fx(i) = gx(i) + ( x(3) - 1.0D+00 ) * gy(i)
      end do

      return
      end
      subroutine p02_gx ( x, g )

c*********************************************************************72
c
cc P02_GX evaluates the underlying function for problem 2.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    29 September 2008
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, double precision X(2), the point at which the function is to
c    be evaluated.
c
c    Output, double precision G(2), the value of the function at X.
c
      implicit none

      double precision g(2)
      double precision pi
      parameter ( pi = 3.141592653589793D+00 )
      double precision x(2)

      g(1) = x(1) * x(1) - x(2) + 1.0D+00
      g(2) = x(1) - cos ( pi * x(2) / 2.0D+00 )

      return
      end
      subroutine p02_jac ( option, nvar, x, jac )

c*********************************************************************72
c
cc P02_JAC evaluates the jacobian for problem 2.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    29 September 2008
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer OPTION, the option index.
c
c    Input, integer NVAR, the number of variables.
c
c    Input, double precision X(NVAR), the argument of the jacobian.
c
c    Output, double precision JAC(NVAR,NVAR), the jacobian matrix evaluated
c    at X.  The NVAR-th row is not set by this routine.
c
      implicit none

      integer nvar

      double precision gy(2)
      integer option
      double precision jac(nvar,nvar)
      double precision pi
      parameter ( pi = 3.141592653589793D+00 )
      double precision x(nvar)
      double precision y(3)

      jac(1,1) = 2.0D+00 * x(1)
      jac(2,1) = 1.0D+00
      jac(1,2) = - 1.0D+00
      jac(2,2) = 0.5D+00 * pi * sin ( 0.5D+00 * pi * x(2) )
c
c  Get the starting point
c
      call p02_start ( option, nvar, y )
c
c  Get the function value at the starting point
c
      call p02_gx ( y, gy )

      jac(1,3) = gy(1)
      jac(2,3) = gy(2)

      return
      end
      subroutine p02_nvar ( option, nvar )

c*********************************************************************72
c
cc P02_NVAR sets the number of variables for problem 2.
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
c    Input, integer OPTION, the option chosen for this problem.
c    For some problems, several options are available.  At least,
c    OPTION = 1 is always legal.
c
c    Output, integer NVAR, the number of variables.
c
      implicit none

      integer nvar
      integer option

      nvar = 3

      return
      end
      subroutine p02_option_num ( option_num )

c*********************************************************************72
c
cc P02_OPTION_NUM returns the number of options for problem 2.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    16 September 2008
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, integer OPTION_NUM, the number of options.
c
      implicit none

      integer option_num

      option_num = 3

      return
      end
      subroutine p02_start ( option, nvar, x )

c*********************************************************************72
c
cc P02_START returns a starting point for problem 2.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    19 September 2008
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer OPTION, the option index.
c
c    Input, integer NVAR, the number of variables.
c
c    Output, double precision X(NVAR), the starting point.
c
      implicit none

      integer nvar

      integer option
      double precision x(nvar)

      if ( option .eq. 1 ) then

        x(1) = 1.0
        x(2) = 0.0
        x(3) = 0.0

      else if ( option .eq. 2 ) then

        x(1) =  1.0
        x(2) = -1.0
        x(3) =  0.0

      else if ( option .eq. 3 ) then

        x(1) = 10.0
        x(2) = 10.0
        x(3) =  0.0

      else

        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'P02_START - Fatal error!'
        write ( *, '(a,i8)' ) '  Unrecognized option  = ', option
        stop

      end if

      return
      end
      subroutine p02_stepsize ( option, h, hmin, hmax )

c*********************************************************************72
c
cc P02_STEPSIZE returns step sizes for problem 2.
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
c    Input, integer OPTION, the option index.
c
c    Output, double precision H, HMIN, HMAX, suggested values for the
c    initial step, the minimum step, and the maximum step.
c
      implicit none

      double precision h
      double precision hmax
      double precision hmin
      integer option

      h =    0.250D+00
      hmin = 0.001D+00
      hmax = 1.000D+00

      return
      end
      subroutine p02_title ( option, title )

c*********************************************************************72
c
cc P02_TITLE sets the title for problem 2.
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
c    Input, integer OPTION, the option index.
c
c    Output, character*(*) TITLE, the title of the problem.
c    TITLE will never be longer than 80 characters.
c
      implicit none

      integer option
      character*(*) title

      if ( option .eq. 1 ) then
        title = 'Boggs function, (1,0,0).'
      else if ( option .eq. 2 ) then
        title = 'Boggs function, (1,-1,0).'
      else if ( option .eq. 3 ) then
        title = 'Boggs function, (10,10,0).'
      else
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'P02_TITLE - Fatal error!'
        write ( *, '(a)' ) '  Unrecognized option number.'
        stop
      end if

      return
      end
      subroutine p03_data(DIRIPC,HFACT,HMAX,HMIN,HTAN,option,
     & IPC,IPRAM,ISIZE,IT,LIM,MAXCON,
     & MAXLIM,MAXSTP,MAXTAR,NCOL,NROW,NVAR,RELERR,XIT,RWORK)

c*********************************************************************72
c
cc P03_DATA sets parameters for the Powell function.
c
c  3.1 the powell function.
c
c  3.2 reference
c
c  m j d powell,
c  a fortran subroutine for solving systems of nonlinear
c    algebraic equations,
c  numerical methods for nonlinear algebraic equations,
c  edited by p rabinowitz
c  gordon and breach, london 1970, pages 115-161
c
c  3.3 the function
c
c  fx is of the form
c
c  fx(x1,x2,x3) = gx(x1,x2) - (1-x3) * gx(x1(0),x2(0)).
c
c  with gx1 = 10000*x1*x2 - 1.0
c       gx2 = exp(-x1) + exp(-x2) - 1.0001
c
c  3.4 options
c
c  option=1  starting point is (3,6,0), it=3, xit=1.0
c  option=2  starting point is (4,5,0), it=3, xit=1.0
c  option=3  starting point is (6,3,0), it=3, xit=1.0
c  option=4  starting point is (1,1,0), it=3, xit=1.0
c
c  3.5  target point
c
c  for all values of option,
c  the target point is (1.098159e-5, 9.106146, 1.0).
c
c  3.8 comments
c
c  seek limit point in x3.
c
      DIMENSION RWORK(ISIZE)
      COMMON /AUXMEM/ GXZERO(2)

      NVAR=3
      DIRIPC=1.0
      HFACT=2.0
      IPC=2
      IT=3
      LIM=3
      MAXCON=35
      MAXLIM=35
      MAXSTP=35
      MAXTAR=2
      NCOL=NVAR
      NROW=NVAR
      RELERR=1.0E-4
      XIT=1.0
c
c  evaluate function at starting point
c
      NEQN=NVAR-1
      call p03_gx (NEQN,RWORK,GXZERO)
      return
      END
      subroutine p03_fun ( option, nvar, x, fx )

c*********************************************************************72
c
cc P03_FUN evaluates the Powell function.
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
      integer nvar

      double precision fx(nvar-1)
      integer option
      double precision x(nvar)

      COMMON /AUXMEM/ GXZERO(2)
      NEQN=NVAR-1
      call p03_gx (NEQN,X,FX)
      X3=X(3)
      DO I=1,NEQN
        FX(I)=FX(I)-(1.0-X3)*GXZERO(I)
      end do
      return
      END
      subroutine p03_gx (NEQN,X,FX)

c*********************************************************************72
c
cc P03_GX is an auxilliary routine for the Powell function.
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
      DIMENSION X(NEQN),FX(NEQN)
      X1=X(1)
      X2=X(2)
      FX(1)=10000.0*X1*X2-1.0
      FX(2)=EXP(-X1)+EXP(-X2)-1.0001
      return
      END
      subroutine p03_jac ( option, nvar, x, jac )

c*********************************************************************72
c
cc P03_JAC evaluates the Powell jacobian.
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
c    Input, integer OPTION, the option index.
c
c    Input, integer NVAR, the number of variables.
c
c    Input, double precision X(NVAR), the argument of the jacobian.
c
c    Output, double precision JAC(NVAR,NVAR), the jacobian matrix evaluated
c    at X.  The NVAR-th row is not set by this routine.
c
      implicit none

      integer nvar

      double precision gxzero(2)
      integer option
      double precision jac(nvar,nvar)
      double precision x(nvar)

      COMMON /AUXMEM/ GXZERO

      jac(1,1)=10000.0*X(2)
      jac(1,2)=10000.0*X(1)
      jac(1,3)=GXZERO(1)
      jac(2,1)=-EXP(-X(1))
      jac(2,2)=-EXP(-X(2))
      jac(2,3)=GXZERO(2)

      return
      END
      subroutine p03_nvar ( option, nvar )

c*********************************************************************72
c
cc P03_NVAR sets the number of variables for problem 3.
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
c    Input, integer OPTION, the option chosen for this problem.
c    For some problems, several options are available.  At least,
c    OPTION = 1 is always legal.
c
c    Output, integer NVAR, the number of variables.
c
      implicit none

      integer nvar
      integer option

      nvar = 3

      return
      end
      subroutine p03_option_num ( option_num )

c*********************************************************************72
c
cc P03_OPTION_NUM returns the number of options for problem 3.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    16 September 2008
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, integer OPTION_NUM, the number of options.
c
      implicit none

      integer option_num

      option_num = 4

      return
      end
      subroutine p03_start ( option, nvar, x )

c*********************************************************************72
c
cc P03_START returns a starting point for problem 3.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    19 September 2008
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer OPTION, the option index.
c
c    Input, integer NVAR, the number of variables.
c
c    Output, double precision X(NVAR), the starting point.
c
      implicit none

      integer nvar

      integer option
      double precision x(nvar)

      if ( option .eq. 1 ) then

        x(1) = 3.0
        x(2) = 6.0
        x(3) = 0.0

      else if ( option .eq. 2 ) then

        x(1) = 4.0
        x(2) = 5.0
        x(3) = 0.0

      else if ( option .eq. 3 ) then

        x(1) = 6.0
        x(2) = 3.0
        x(3) = 0.0

      else if ( option .eq. 4 ) then

        x(1) = 1.0
        x(2) = 1.0
        x(3) = 0.0

      else

        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'P03_START - Fatal error!'
        write ( *, '(a,i8)' ) '  Unrecognized option  = ', option
        stop

      end if

      return
      end
      subroutine p03_stepsize ( option, h, hmin, hmax )

c*********************************************************************72
c
cc P03_STEPSIZE returns step sizes for problem 3.
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
c    Input, integer OPTION, the option index.
c
c    Output, double precision H, HMIN, HMAX, suggested values for the
c    initial step, the minimum step, and the maximum step.
c
      implicit none

      double precision h
      double precision hmax
      double precision hmin
      integer option

      h =    0.50000D+00
      hmin = 0.00025D+00
      hmax = 3.00000D+00

      return
      end
      subroutine p03_title ( option, title )

c*********************************************************************72
c
cc P03_TITLE sets the title for problem 3.
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
c    Input, integer OPTION, the option index.
c
c    Output, character*(*) TITLE, the title of the problem.
c    TITLE will never be longer than 80 characters.
c
      implicit none

      integer option
      character*(*) title

      if ( option .eq. 1 ) then
        title = 'Powell function, (3,6,0).'
      else if ( option .eq. 2 ) then
        title = 'Powell function, (4,5,0).'
      else if ( option .eq. 3 ) then
        title = 'Powell function, (6,3,0).'
      else if ( option .eq. 4 ) then
        title = 'Powell function, (1,1,0).'
      else
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'P03_TITLE - Fatal error!'
        write ( *, '(a)' ) '  Unrecognized option number.'
        stop
      end if

      return
      end
      subroutine p04_data (DIRIPC,HFACT,HMAX,HMIN,HTAN,option,
     & IPC,IPRAM,ISIZE,IT,LIM,MAXCON,
     & MAXLIM,MAXSTP,MAXTAR,NCOL,NROW,NVAR,RELERR,XIT,RWORK)

c*********************************************************************72
c
cc P04_DATA sets parameters for the Broyden function.
c
c  4.1 the broyden function.
c
c  4.2 reference
c
c  c broyden,
c  a new method of solving nonlinear simultaneous equations,
c  comput j 12, 1969, pages 94-99.
c
c  4.3 the function
c
c  fx is of the form
c
c  fx(x1,x2,x3) = gx(x1,x2) - (1-x3) * gx(x1(0),x2(0)).
c
c  with gx1 = .5*sin(x1*x2) - x2/pi - x1
c       gx2 = (1-1/(4*pi))*(exp(2*x1)-exp(1)) + exp(1)*x2/pi
c           - 2*exp(1)*x1
c
c  4.4 options
c
c  option=1  starting point is (.4,3,0), it=3, xit=1.0.
c
c  4.5 target point
c
c  the target point is (-.2207014, .8207467, 1.0)
c
c  4.8 comments
c
c  program seems ok.
c
      DIMENSION RWORK(ISIZE)
      COMMON /AUXMEM/ GXZERO(2)

      NVAR=3
      DIRIPC=-1.0
      HFACT=3.0
      IPC=2
      IT=3
      LIM=0
      MAXCON=20
      MAXLIM=20
      MAXSTP=20
      MAXTAR=1
      NCOL=NVAR
      NROW=NVAR
      RELERR=1.0E-5
      XIT=1.0
c
c  get function value at starting point
c
      NEQN=NVAR-1
      call p04_gx (NEQN,RWORK,GXZERO)
      return
      END
      subroutine p04_fun ( option, nvar, x, fx )

c*********************************************************************72
c
cc P04_FUN evaluates the Broyden function.
c
      integer nvar

      double precision fx(nvar-1)
      integer option
      double precision x(nvar)

      COMMON /AUXMEM/ GXZERO(2)
      NEQN=NVAR-1
      call p04_gx (NEQN,X,FX)
      X3=X(3)
      DO I=1,NEQN
        FX(I)=FX(I)-(1.0-X3)*GXZERO(I)
      end do
      return
      END
      subroutine p04_gx (NEQN,X,FX)

c*********************************************************************72
c
cc P04_GX is an auxilliary routine for the Broyden function.
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
      real enat
      parameter ( enat = 2.7182818284E+00 )
      real FX(NEQN)
      real pi
      parameter ( pi = 3.1415926535E+00 )
      DIMENSION X(NEQN)

      X1=X(1)
      X2=X(2)
      FX(1)=.5*SIN(X1*X2)-X2/(2.0*PI)-X1
      FX(2)=(1.0-1.0/(4.0*PI))*(EXP(2.0*X1)-ENAT)
     &      +ENAT*X2/PI -2.0*ENAT*X1
      return
      END
      subroutine p04_jac ( option, nvar, x, jac )

c*********************************************************************72
c
cc P04_JAC evaluates the Broyden jacobian.
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
c    Input, integer OPTION, the option index.
c
c    Input, integer NVAR, the number of variables.
c
c    Input, double precision X(NVAR), the argument of the jacobian.
c
c    Output, double precision JAC(NVAR,NVAR), the jacobian matrix evaluated
c    at X.  The NVAR-th row is not set by this routine.
c
      implicit none

      integer nvar

      real enat
      parameter ( enat = 2.7182818284E+00 )
      real gxzero(2)
      double precision jac(nvar,nvar)
      integer option
      real pi
      parameter ( pi = 3.1415926535E+00 )
      double precision x(nvar)

      COMMON /AUXMEM/ GXZERO

      jac(1,1)=0.5*X(1)*COS(X(1)*X(2))-1.0
      jac(1,2)=0.5*X(1)*COS(X(1)*X(2))-1.0/(2.0*PI)
      jac(1,3)=GXZERO(1)
      jac(2,1)=(1.0-1.0/(4.0*PI))*2.0*EXP(2.0*X(1))-2.0*ENAT
      jac(2,2)=ENAT/PI
      jac(2,3)=GXZERO(2)

      return
      END
      subroutine p04_nvar ( option, nvar )

c*********************************************************************72
c
cc P04_NVAR sets the number of variables for problem 4.
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
c    Input, integer OPTION, the option chosen for this problem.
c    For some problems, several options are available.  At least,
c    OPTION = 1 is always legal.
c
c    Output, integer NVAR, the number of variables.
c
      implicit none

      integer nvar
      integer option

      nvar = 3

      return
      end
      subroutine p04_option_num ( option_num )

c*********************************************************************72
c
cc P04_OPTION_NUM returns the number of options for problem 4.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    16 September 2008
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, integer OPTION_NUM, the number of options.
c
      implicit none

      integer option_num

      option_num = 1

      return
      end
      subroutine p04_start ( option, nvar, x )

c*********************************************************************72
c
cc P04_START returns a starting point for problem 4.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    19 September 2008
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer OPTION, the option index.
c
c    Input, integer NVAR, the number of variables.
c
c    Output, double precision X(NVAR), the starting point.
c
      implicit none

      integer nvar

      integer option
      double precision x(nvar)

      if ( option .eq. 1 ) then

        x(1) = 0.4
        x(2) = 3.0
        x(3) = 0.0

      else

        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'P04_START - Fatal error!'
        write ( *, '(a,i8)' ) '  Unrecognized option  = ', option
        stop

      end if

      return
      end
      subroutine p04_stepsize ( option, h, hmin, hmax )

c*********************************************************************72
c
cc P04_STEPSIZE returns step sizes for problem 4.
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
c    Input, integer OPTION, the option index.
c
c    Output, double precision H, HMIN, HMAX, suggested values for the
c    initial step, the minimum step, and the maximum step.
c
      implicit none

      double precision h
      double precision hmax
      double precision hmin
      integer option

      h =     0.300D+00
      hmin =  0.001D+00
      hmax = 25.000D+00

      return
      end
      subroutine p04_title ( option, title )

c*********************************************************************72
c
cc P04_TITLE sets the title for problem 4.
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
c    Input, integer OPTION, the option index.
c
c    Output, character*(*) TITLE, the title of the problem.
c    TITLE will never be longer than 80 characters.
c
      implicit none

      integer option
      character*(*) title

      if ( option .eq. 1 ) then
        title = 'Broyden function.'
      else
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'P04_TITLE - Fatal error!'
        write ( *, '(a)' ) '  Unrecognized option number.'
        stop
      end if

      return
      end
      subroutine p05_data(DIRIPC,HFACT,HMAX,HMIN,HTAN,option,
     & IPC,IPRAM,ISIZE,IT,LIM,MAXCON,
     & MAXLIM,MAXSTP,MAXTAR,NCOL,NROW,NVAR,RELERR,XIT,RWORK)

c*********************************************************************72
c
cc P05_DATA sets parameters for the Wacker function.
c
c  5.1 the wacker function.
c
c  5.2 reference
c
c  h wacker, e zarzer, w zulehner,
c  optimal stepsize control for the globalized newton method,
c  in   continuation methods,
c  edited by h wacker,
c  academic press, new york, 1978,  pages 249-276.
c
c  5.3 the function
c
c  the function is of the form
c
c  fx1= (1-aval*x4)*x1 +x4*exp(-x2)/3.0 + x4*(aval-1-1/(3*e))
c  fx2= (1-aval*x4)*x2 -x4*ln(1+x3*x3)/5 +x4*(aval-1-ln(2)/5)
c  fx3= (1-aval*x3)*x3 +x4*sin(x1) +x4*(aval-1-sin(1))
c
c  with aval a parameter.
c
c  5.4 options
c
c  option=1  aval=0.1, it=3, xit=1.0.
c  option=2  aval=0.5, it=3, xit=1.0
c  option=3  aval=1.0, it=3, xit=1.0
c
c  5.5 target point
c
c  for option=1, 1.147009,  1.431931,  1.000000, 1.084425.
c  for option=2, 0.2412182, 0.4558247, 1.000000, 0.4534797.
c  for option=3, 0.0000000, 0.0000000, 1.000000, 0.000000.
c
c  5.6 limit points
c
c  for lim=4,
c  for option=3, -0.07109918, 0.06921115, 0.5009694, 0.2739685.
c
c  5.8 comments
c
c  program seems ok.
c  perhaps later, add aval as x(5).
c
      DIMENSION RWORK(ISIZE)
      COMMON /AUXMEM/ AVAL

      NVAR=4
      IF(option.EQ.1)AVAL=0.1
      IF(option.EQ.2)AVAL=0.5
      IF(option.EQ.3)AVAL=1.0
      DIRIPC=1.0
      HFACT=3.0
      IPC=3
      IT=3
      LIM=4
      MAXCON=20
      MAXLIM=20
      MAXSTP=20
      MAXTAR=1
      NCOL=NVAR
      NROW=NVAR
      RELERR=1.0E-5
      XIT=1.0
      return
      END
      subroutine p05_fun ( option, nvar, x, fx ))

c*********************************************************************72
c
cc P05_FUN evaluates the Wacker function.
c
      integer nvar

      double precision aval
      real enat
      parameter ( enat = 2.7182818284E+00 )
      double precision fx(nvar-1)
      integer option
      double precision x(nvar)

      COMMON /AUXMEM/ AVAL

      X1=X(1)
      X2=X(2)
      X3=X(3)
      X4=X(4)
      TERM1=(1.0-AVAL*X4)*X1
      TERM2=X4*EXP(-X2)/3.0
      TERM3=X4*(AVAL-1.0-1.0/(3.0*ENAT))
      FX(1)=TERM1+TERM2+TERM3
      TERM1=(1.0-AVAL*X4)*X2
      TERM2=-X4*ALOG(1.0+X3*X3)/5.0
      TERM3=X4*(AVAL-1.0-ALOG(2.0)/5.0)
      FX(2)=TERM1+TERM2+TERM3
      TERM1=(1.0-AVAL*X3)*X3
      TERM2=X4*SIN(X1)
      TERM3=X4*(AVAL-1.0-SIN(1.0))
      FX(3)=TERM1+TERM2+TERM3
      return
      END
      subroutine p05_jac ( option, nvar, x, jac )

c*********************************************************************72
c
cc P05_JAC evaluates the Wacker jacobian.
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
c    Input, integer OPTION, the option index.
c
c    Input, integer NVAR, the number of variables.
c
c    Input, double precision X(NVAR), the argument of the jacobian.
c
c    Output, double precision JAC(NVAR,NVAR), the jacobian matrix evaluated
c    at X.  The NVAR-th row is not set by this routine.
c
      implicit none

      integer nvar

      integer option
      double precision jac(nvar,nvar)
      double precision x(nvar)

      real enat
      parameter ( enat = 2.7182818284E+00 )

      COMMON /AUXMEM/ AVAL

      X1=X(1)
      X2=X(2)
      X3=X(3)
      X4=X(4)
      jac(1,1)=1.0-AVAL*X4
      jac(1,2)=-X4*EXP(-X2)/3.0
      jac(1,3)=0.0
      jac(1,4)=-AVAL*X1+EXP(-X2)/3.0+(AVAL-1.0-1.0/(3.0*ENAT))
      jac(2,1)=0.0
      jac(2,2)=1.0-AVAL*X4
      jac(2,3)=-2.0*X3*X4/(5.0*(1.0+X3*X3))
      jac(2,4)=-AVAL*X2-ALOG(1.0+X3*X3)/5.0+(AVAL-1.0-ALOG(2.0)/5.0)
      jac(3,1)=X4*COS(X1)
      jac(3,2)=0.0
      jac(3,3)=1.0-2.0*AVAL*X3
      jac(3,4)=SIN(X1)+(AVAL-1.0-SIN(1.0))

      return
      END
      subroutine p05_nvar ( option, nvar )

c*********************************************************************72
c
cc P05_NVAR sets the number of variables for problem 5.
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
c    Input, integer OPTION, the option chosen for this problem.
c    For some problems, several options are available.  At least,
c    OPTION = 1 is always legal.
c
c    Output, integer NVAR, the number of variables.
c
      implicit none

      integer nvar
      integer option

      nvar = 4

      return
      end
      subroutine p05_option_num ( option_num )

c*********************************************************************72
c
cc P05_OPTION_NUM returns the number of options for problem 5.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    16 September 2008
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, integer OPTION_NUM, the number of options.
c
      implicit none

      integer option_num

      option_num = 3

      return
      end
      subroutine p05_start ( option, nvar, x )

c*********************************************************************72
c
cc P05_START returns a starting point for problem 5.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    19 September 2008
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer OPTION, the option index.
c
c    Input, integer NVAR, the number of variables.
c
c    Output, double precision X(NVAR), the starting point.
c
      implicit none

      integer nvar

      integer i
      integer option
      double precision x(nvar)

      if ( 1 .le. option .and. option .le. 3 ) then

        do i = 1, nvar
          x(i) = 0.0
        end do

      else

        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'P05_START - Fatal error!'
        write ( *, '(a,i8)' ) '  Unrecognized option  = ', option
        stop

      end if

      return
      end
      subroutine p05_stepsize ( option, h, hmin, hmax )

c*********************************************************************72
c
cc P05_STEPSIZE returns step sizes for problem 5.
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
c    Input, integer OPTION, the option index.
c
c    Output, double precision H, HMIN, HMAX, suggested values for the
c    initial step, the minimum step, and the maximum step.
c
      implicit none

      double precision h
      double precision hmax
      double precision hmin
      integer option

      h =     0.300D+00
      hmin =  0.001D+00
      hmax = 25.000D+00

      return
      end
      subroutine p05_title ( option, title )

c*********************************************************************72
c
cc P05_TITLE sets the title for problem 5.
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
c    Input, integer OPTION, the option index.
c
c    Output, character*(*) TITLE, the title of the problem.
c    TITLE will never be longer than 80 characters.
c
      implicit none

      integer option
      character*(*) title

      if ( option .eq. 1 ) then
        title = 'Wacker function, A = 0.1.'
      else if ( option .eq. 2 ) then
        title = 'Wacker function, A = 0.5.'
      else if ( option .eq. 3 ) then
        title = 'Wacker function, A = 1.0.'
      else
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'P05_TITLE - Fatal error!'
        write ( *, '(a)' ) '  Unrecognized option number.'
        stop
      end if

      return
      end
      subroutine p06_data (DIRIPC,HFACT,HMAX,HMIN,HTAN,option,
     & IPC,IPRAM,ISIZE,IT,LIM,MAXCON,
     & MAXLIM,MAXSTP,MAXTAR,NCOL,NROW,NVAR,RELERR,XIT,RWORK)

c*********************************************************************72
c
cc P06_DATA sets parameters for the Aircraft Stability function.
c
c  The aircraft stability problem.
c
c  Reference:
c
c  a a schy, m e hannah,
c  prediction of jump phenomena in roll-coupled maneuvers
c    of airplanes,
c  journal of aircraft, 14, 1977,  pages 375-382.
c
c  j w young, a a schy, k g johnson
c  prediction of jump phenomena in aircraft maneuvers,
c    including nonlinear aerodynamic effects,
c  journal of guidance and control 1, 1978, pages 26-31.
c
c  The function
c
c  the equations describe the behavior of an aircraft under
c  the control of a pilot. the variables are-
c
c  x(1) = roll
c  x(2) = pitch
c  x(3) = yaw
c  x(4) = angle of attack
c  x(5) = sideslip
c  x(6) = elevator
c  x(7) = aileron
c  x(8) = rudder
c
c  the function is of the following form
c
c  for indices i=1 through 5,
c
c  f(i)=sum on j and k of  b(i,j)*x(j)+phi(i,j,k)*x(j)*x(k)
c
c  f(6)=x(ifix1)-val1
c  f(7)=x(ifix2)-val2
c
c  Options
c
c  option=1  ifix1=6, val1=-.050.
c  option=2  ifix1=6, val1=-.008.
c  option=3  ifix1=6, val1=0.000.
c  option=4  ifix1=6, val1=0.050.
c  option=5  ifix1=6, val1=0.100.
c  option=6  ifix1=6, val1=-0.01250. (check for bifurcation)
c
c  for all choices, lim=7, ifix2=8, val2=0.0.
c
c  Limit points
c
c  melhem lists the following limit points in x7
c  (note that melhem has barray(4,1)=1.0, barray(4,2)=0.0)
c
c     x1       x2        x3        x4        x5      x6     x7      x8
c
c  -2.9691  0.83074  -0.072748  0.41029  -0.26880  -.05  0.50919    0.0
c  -2.8158 -0.17481  -0.089469  0.026319  0.070951 -.008 0.20442    0.0
c  -3.7571 -0.64911  -0.39350   0.091810  0.19685  -.008 -.0038238  0.0
c  -4.1637  0.092284 -0.092610  0.022402 -0.017106 -.008 0.37823    0.0
c  -2.5839 -0.22128  -0.054079  0.013524  0.090871 0.0   0.18608    0.0
c  -3.9007 -1.1421   -0.57863   0.13284   0.32685  0.0  -0.50703    0.0
c  -2.3610 -0.72360   0.032739 -0.039108  0.29347  0.05  0.29272    0.0
c  -2.2982  1.4033    0.063244 -0.079383  0.58336  0.10  0.58336    0.0
c
c  computational results with barray(4,1)=0.0, barray(4,2)=1.0 are
c
c     x1       x2        x3        x4        x5      x6     x7      x8
c
c   2.9648  0.82556   0.07366   0.041309  0.26734 -0.050 -.050481   0.0
c   2.8173 -0.17628   0.08992   0.026429 -0.07147 -0.008 -.204973   0.0
c   3.7579 -0.65541   0.38658   0.092520 -0.19867 -0.008 0.006200   0.0
c   4.1638  0.08913   0.09480   0.022888  0.16232 -0.008 -.377660   0.0
c   2.5873 -0.22354   0.05468   0.013676 -0.09168  0.000 -.186908   0.0
c   3.9005 -1.14815   0.58156   0.133516 -0.32858   0.000  .510158  0.0
c   2.3639 -0.72974  -0.31604  -0.038785 -0.29583   0.050 -.295772  0.0
c   2.2992 -1.41023  -0.06184  -0.079009 -0.58629   0.100 -.689717  0.0
c
c  Bifurcation points
c
c  rheinboldt lists
c
c   4.482   0.1632    0.02373   0.006205  0.03527 -0.0006177 -0.3986 0.0
c   3.319  -0.1869    0.1605    0.04379  -0.06888 -0.01250   -0.2374 0.0
c   4.466   0.1467    0.04045   0.009777  0.03089 -0.006129  -0.3995 0.0
c  -3.325   0.1880   -0.1614    0.04395   0.06911 -0.01247    0.2367 0.0
c
      DIMENSION RWORK(ISIZE)
      COMMON /AUXMEM/ IFIX1,IFIX2,VAL1,VAL2,BARRAY(5,8)

      NVAR=8

      IFIX1=6
      IFIX2=8
      IF(option.EQ.1)VAL1=-.050
      IF(option.EQ.2)VAL1=-.008
      IF(option.EQ.3)VAL1=0.000
      IF(option.EQ.4)VAL1=0.050
      IF(option.EQ.5)VAL1=0.100
      IF(option.EQ.6)VAL1=-0.01250
      VAL2=0.0
      DIRIPC=-1.0
      HFACT=2.0
      IF(option.NE.6)HMAX=0.5
      IF(option.EQ.6)HMAX=0.25
      HMIN=.001
      HTAN=0.25
      IPC=7
      IT=0
      IF(option.NE.6)LIM=7
      IF(option.EQ.6)LIM=0
      MAXCON=25
      IF(option.EQ.1)MAXLIM=1
      IF(option.EQ.2)MAXLIM=3
      IF(option.EQ.3)MAXLIM=2
      IF(option.EQ.4)MAXLIM=1
      IF(option.EQ.5)MAXLIM=1
      IF(option.EQ.6)MAXLIM=25
      MAXSTP=25
      MAXTAR=25
      NCOL=NVAR
      NROW=NVAR
      RELERR=1.0E-5
      XIT=0.0
c
c  set approximate starting point
c
      DO I=1,NVAR
        RWORK(I)=0.0
      enddo

      RWORK(IFIX1)=VAL1
      RWORK(IFIX2)=VAL2
c
c  set b array
c
      BARRAY(1,1)=-3.933
      BARRAY(2,1)=0.0
      BARRAY(3,1)=0.002
      BARRAY(4,1)=0.0
      BARRAY(5,1)=0.0
      BARRAY(1,2)=0.107
      BARRAY(2,2)=-0.987
      BARRAY(3,2)=0.0
      BARRAY(4,2)=1.0
      BARRAY(5,2)=0.0
      BARRAY(1,3)=0.126
      BARRAY(2,3)=0.0
      BARRAY(3,3)=-0.235
      BARRAY(4,3)=0.0
      BARRAY(5,3)=-1.0
      BARRAY(1,4)=0.0
      BARRAY(2,4)=-22.95
      BARRAY(3,4)=0.0
      BARRAY(4,4)=-1.0
      BARRAY(5,4)=0.0
      BARRAY(1,5)=-9.99
      BARRAY(2,5)=0.0
      BARRAY(3,5)=5.67
      BARRAY(4,5)=0.0
      BARRAY(5,5)=-0.196
      BARRAY(1,6)=0.0
      BARRAY(2,6)=-28.37
      BARRAY(3,6)=0.0
      BARRAY(4,6)=-0.168
      BARRAY(5,6)=0.0
      BARRAY(1,7)=-45.83
      BARRAY(2,7)=0.0
      BARRAY(3,7)=-0.921
      BARRAY(4,7)=0.0
      BARRAY(5,7)=-0.0071
      BARRAY(1,8)=-7.64
      BARRAY(2,8)=0.0
      BARRAY(3,8)=-6.51
      BARRAY(4,8)=0.0
      BARRAY(5,8)=0.0
      return
      END
      subroutine p06_fun ( option, nvar, x, fx )

c*********************************************************************72
c
cc P06_FUN evaluates the aircraft stability function.
c
      integer nvar

      double precision fx(nvar-1)
      integer option
      double precision x(nvar)

      COMMON /AUXMEM/ IFIX1,IFIX2,VAL1,VAL2,BARRAY(5,8)
c
c  compute linear terms
c
      DO I=1,5
        FX(I)=0.0
        DO J=1,8
          FX(I)=FX(I)+BARRAY(I,J)*X(J)
        end do
      end do
c
c  compute nonlinear terms
c
      PHI=-.727*X(2)*X(3)+8.39*X(3)*X(4)-684.4*X(4)*X(5)+63.5*X(4)*X(7)
      FX(1)=FX(1)+PHI
      PHI=.949*X(1)*X(3)+.173*X(1)*X(5)
      FX(2)=FX(2)+PHI
      PHI=-.716*X(1)*X(2)-1.578*X(1)*X(4)+1.132*X(4)*X(7)
      FX(3)=FX(3)+PHI
      PHI=-X(1)*X(5)
      FX(4)=FX(4)+PHI
      PHI=X(1)*X(4)
      FX(5)=FX(5)+PHI
c
c  set function values for two fixed variables
c
      FX(6)=X(IFIX1)-VAL1
      FX(7)=X(IFIX2)-VAL2

      return
      END
      subroutine p06_jac ( option, nvar, x, jac )

c*********************************************************************72
c
cc P06_JAC evaluates the aircraft stability jacobian.
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
c    Input, integer OPTION, the option index.
c
c    Input, integer NVAR, the number of variables.
c
c    Input, double precision X(NVAR), the argument of the jacobian.
c
c    Output, double precision JAC(NVAR,NVAR), the jacobian matrix evaluated
c    at X.  The NVAR-th row is not set by this routine.
c
      implicit none

      integer nvar

      integer option
      double precision jac(nvar,nvar)
      double precision x(nvar)

      COMMON /AUXMEM/ IFIX1,IFIX2,VAL1,VAL2,BARRAY(5,8)
c
c  compute terms from linear part of function
c
      DO I=1,5
        DO J=1,8
          jac(I,J)=BARRAY(I,J)
        end do
      end do

      DO I=6,7
        DO J=1,8
          jac(I,J)=0.0
        end do
      end do

      jac(6,IFIX1)=1.0
      jac(7,IFIX2)=1.0
c
c  compute terms from nonlinear part of function
c
      jac(1,2)=jac(1,2)-.727*X(3)
      jac(1,3)=jac(1,3)-.727*X(2)+8.39*X(4)
      jac(1,4)=jac(1,4)+8.39*X(3)-684.4*X(5)+63.5*X(7)
      jac(1,5)=jac(1,5)-684.4*X(4)
      jac(1,7)=jac(1,7)+63.5*X(4)
      jac(2,1)=jac(2,1)+.949*X(3)+.173*X(5)
      jac(2,3)=jac(2,3)+.949*X(1)
      jac(2,5)=jac(2,5)+.173*X(1)
      jac(3,1)=jac(3,1)-.716*X(2)-1.578*X(4)
      jac(3,2)=jac(3,2)-.716*X(1)
      jac(3,4)=jac(3,4)-1.578*X(1)+1.132*X(7)
      jac(3,7)=jac(3,7)+1.132*X(4)
      jac(4,1)=jac(4,1)-X(5)
      jac(4,5)=jac(4,5)-X(1)
      jac(5,1)=jac(5,1)+X(4)
      jac(5,4)=jac(5,4)+X(1)
      return
      END
      subroutine p06_nvar ( option, nvar )

c*********************************************************************72
c
cc P06_NVAR sets the number of variables for problem 6.
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
c    Input, integer OPTION, the option chosen for this problem.
c    For some problems, several options are available.  At least,
c    OPTION = 1 is always legal.
c
c    Output, integer NVAR, the number of variables.
c
      implicit none

      integer nvar
      integer option

      nvar = 8

      return
      end
      subroutine p06_option_num ( option_num )

c*********************************************************************72
c
cc P06_OPTION_NUM returns the number of options for problem 6.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    16 September 2008
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, integer OPTION_NUM, the number of options.
c
      implicit none

      integer option_num

      option_num = 6

      return
      end
      subroutine p06_start ( option, nvar, x )

c*********************************************************************72
c
cc P06_START returns a starting point for problem 6.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    19 September 2008
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer OPTION, the option index.
c
c    Input, integer NVAR, the number of variables.
c
c    Output, double precision X(NVAR), the starting point.
c
      implicit none

      integer nvar

      integer option
      double precision x(nvar)

      if ( option == 1 ) then

        x(1) =  0.000009180674117D+00
        x(3) =  0.051206082777927D+00
        x(3) = -0.000003471121026D+00
        x(4) =  0.059606082627213D+00
        x(5) =  0.000016537664587D+00
        x(6) = -0.050000000000000D+00
        x(7) =  0.000109431378424D+00
        x(8) =  0.000000000000000D+00

      else if ( option == 2 ) then

        x(1) =  0.000001548268247D+00
        x(2) =  0.008192973225663D+00
        x(3) = -0.000000682134573D+00
        x(4) =  0.009536973221178D+00
        x(5) =  0.000002896734870D+00
        x(6) = -0.008000000000000D+00
        x(7) =  0.000018188778989D+00
        x(8) =  0.000000000000000D+00

      else if ( option == 3 ) then

         x(1) = 0.0D+00
         x(2) = 0.0D+00
         x(3) = 0.0D+00
         x(4) = 0.0D+00
         x(5) = 0.0D+00
         x(6) = 0.0D+00
         x(7) = 0.0D+00
         x(8) = 0.0D+00

      else if ( option == 4 ) then

        x(1) = -0.000010655314069D+00
        x(2) = -0.051206082422980D+00
        x(3) =  0.000005600187501D+00
        x(4) = -0.059606082643400D+00
        x(5) = -0.000020891016199D+00
        x(6) =  0.050000000000000D+00
        x(7) = -0.000122595323216D+00
        x(8) =  0.000000000000000D+00

      else if ( option == 5 ) then

        x(1) = -0.000027083319493D+00
        x(2) = -0.102412164106124D+00
        x(3) =  0.000014540858026D+00
        x(4) = -0.119212165322433D+00
        x(5) = -0.000048014067202D+00
        x(6) =  0.100000000000000D+00
        x(7) = -0.000267808407544D+00
        x(8) =  0.000000000000000D+00

      else

        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'P06_START - Fatal error!'
        write ( *, '(a,i8)' ) '  Unrecognized option  = ', option
        stop

      end if

      return
      end
      subroutine p06_stepsize ( option, h, hmin, hmax )

c*********************************************************************72
c
cc P06_STEPSIZE returns step sizes for problem 6.
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
c    Input, integer OPTION, the option index.
c
c    Output, double precision H, HMIN, HMAX, suggested values for the
c    initial step, the minimum step, and the maximum step.
c
      implicit none

      double precision h
      double precision hmax
      double precision hmin
      integer option

      h =    0.250D+00
      hmin = 0.001D+00
      hmax = 0.500D+00

      return
      end
      subroutine p06_title ( option, title )

c*********************************************************************72
c
cc P06_TITLE sets the title for problem 6.
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
c    Input, integer OPTION, the option index.
c
c    Output, character*(*) TITLE, the title of the problem.
c    TITLE will never be longer than 80 characters.
c
      implicit none

      integer option
      character*(*) title

      if ( option .eq. 1 ) then
        title = 'Aircraft function, x(6) = - 0.0500.'
      else if ( option .eq. 2 ) then
        title = 'Aircraft function, x(6) = - 0.0080.'
      else if ( option .eq. 3 ) then
        title = 'Aircraft function, x(6) =   0.0000.'
      else if ( option .eq. 4 ) then
        title = 'Aircraft function, x(6) = + 0.0500.'
      else if ( option .eq. 5 ) then
        title = 'Aircraft function, x(6) = + 0.1000.'
      else if ( option .eq. 6 ) then
        title = 'Aircraft function, x(6) = - 0.0125.'
      else
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'P06_TITLE - Fatal error!'
        write ( *, '(a)' ) '  Unrecognized option number.'
        stop
      end if

      return
      end
      subroutine p07_data (DIRIPC,HFACT,HMAX,HMIN,HTAN,option,
     & IPC,IPRAM,ISIZE,IT,LIM,MAXCON,
     & MAXLIM,MAXSTP,MAXTAR,NCOL,NROW,NVAR,RELERR,XIT,RWORK)

c*********************************************************************72
c
cc P07_DATA sets parameters for the cell kinetics function.
c
c  7.1 cell kinetics problem.
c
c  7.2 reference
c
c  werner c rheinboldt,
c  solution fields of nonlinear equations and continuation methods,
c  siam journal numerical analysis, 17, 1980, pages 221-237.
c
c  7.3 the function
c
c  the function is of the form
c
c  fx(i)=a(i,j)*x(j) + rho(x(i)) - x(nvar)
c
c  with tridiagonal a.
c
c  7.4 options
c
c  option=1 lim=nvar.
c  option=2 lim=0. (check for bifurcation)
c
c  7.6 limit points
c
c  there are two limit points for lim=nvar.
c
c  1.048362, 1.048362, 1.048362, 1.048362, 1.048362, 34.35693.
c  8.822219, 8.822219, 8.822219, 8.822219, 8.822218, 18.88707.
c
c  7.7 bifurcation points.
c
c  there are four bifurcation points.
c
c  7.8 comments
c
c  program seems ok.
c
      integer nvar
      double precision RWORK(ISIZE)

      NVAR=6
      DIRIPC=1.0
      HFACT=2.0
      IPC=NVAR
      IT=0
      IF(option.EQ.1)LIM=NVAR
      IF(option.EQ.2)LIM=0
      MAXCON=30
      IF(option.EQ.1)MAXLIM=2
      IF(option.EQ.2)MAXLIM=30
      MAXSTP=30
      MAXTAR=30
      NCOL=NVAR
      NROW=NVAR
      RELERR=1.0E-5
      XIT=0.0

      return
      END
      subroutine p07_fun ( option, nvar, x, fx )

c*********************************************************************72
c
cc P07_FUN evaluates the cell kinetics function.
c
      integer nvar

      double precision fx(nvar-1)
      integer i
      integer neqn
      integer option
      double precision x(nvar)

      NEQN=NVAR-1
      DO I=1,NEQN
        FX(I)=(100.0*X(I))/(1.0+X(I)+X(I)*X(I)) - X(NVAR)
      end do

      FX(1)=FX(1)+2.0*X(1)-X(2)

      NEQNM1=NVAR-2

      DO I=2, neqn - 1
        FX(I)=FX(I)-X(I-1)+3.0*X(I)-X(I+1)
      end do

      FX(NEQN)=FX(NEQN)-X(NEQNM1)+2.0*X(NEQN)

      return
      END
      subroutine p07_jac ( option, nvar, x, jac )

c*********************************************************************72
c
cc P07_JAC evaluates the cell kinetics jacobian.
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
c    Input, integer OPTION, the option index.
c
c    Input, integer NVAR, the number of variables.
c
c    Input, double precision X(NVAR), the argument of the jacobian.
c
c    Output, double precision JAC(NVAR,NVAR), the jacobian matrix evaluated
c    at X.  The NVAR-th row is not set by this routine.
c
      implicit none

      integer nvar

      integer i
      integer im1
      integer ip1
      integer j
      double precision jac(nvar,nvar)
      integer neqn
      integer neqnm1
      integer option
      double precision x(nvar)
      double precision xi

      NEQN=NVAR-1
      NEQNM1=NVAR-2

      do j = 1, nvar
        do i = 1, nvar
          jac(i,j) = 0.0E+00
        end do
      end do

      DO I=1,neqn
        XI=X(I)
        jac(I,I)=100.0*(1.0-XI*XI)/(1.0+XI+XI*XI)**2
      end do

      jac(1,1)=jac(1,1)+2.0
      jac(1,2)=jac(1,2)-1.0
      jac(1,NVAR)=jac(1,NVAR)-1.0

      DO I=2,NEQNM1
        IM1=I-1
        IP1=I+1
        jac(I,IM1)=jac(I,IM1)-1.0
        jac(I,I)=jac(I,I)+3.0
        jac(I,IP1)=jac(I,IP1)-1.0
        jac(I,NVAR)=jac(I,NVAR)-1.0
      end do

      jac(NEQN,NEQNM1)=jac(NEQN,NEQNM1)-1.0
      jac(NEQN,NEQN)=jac(NEQN,NEQN)+2.0
      jac(NEQN,NVAR)=jac(NEQN,NVAR)-1.0

      return
      END
      subroutine p07_nvar ( option, nvar )

c*********************************************************************72
c
cc P07_NVAR sets the number of variables for problem 7.
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
c    Input, integer OPTION, the option chosen for this problem.
c    For some problems, several options are available.  At least,
c    OPTION = 1 is always legal.
c
c    Output, integer NVAR, the number of variables.
c
      implicit none

      integer nvar
      integer option

      nvar = 6

      return
      end
      subroutine p07_option_num ( option_num )

c*********************************************************************72
c
cc P07_OPTION_NUM returns the number of options for problem 7.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    16 September 2008
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, integer OPTION_NUM, the number of options.
c
      implicit none

      integer option_num

      option_num = 2

      return
      end
      subroutine p07_start ( option, nvar, x )

c*********************************************************************72
c
cc P07_START returns a starting point for problem 7.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    19 September 2008
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer OPTION, the option index.
c
c    Input, integer NVAR, the number of variables.
c
c    Output, double precision X(NVAR), the starting point.
c
      implicit none

      integer nvar

      integer i
      integer option
      double precision x(nvar)

      if ( 1 .le. option .and. option .le. 2 ) then

        do i = 1, nvar
          x(i) = 0.0
        end do

      else

        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'P07_START - Fatal error!'
        write ( *, '(a,i8)' ) '  Unrecognized option  = ', option
        stop

      end if

      return
      end
      subroutine p07_stepsize ( option, h, hmin, hmax )

c*********************************************************************72
c
cc P07_STEPSIZE returns step sizes for problem 7.
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
c    Input, integer OPTION, the option index.
c
c    Output, double precision H, HMIN, HMAX, suggested values for the
c    initial step, the minimum step, and the maximum step.
c
      implicit none

      double precision h
      double precision hmax
      double precision hmin
      integer option

      h =    1.000D+00
      hmin = 0.001D+00
      hmax = 2.000D+00

      return
      end
      subroutine p07_title ( option, title )

c*********************************************************************72
c
cc P07_TITLE sets the title for problem 7.
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
c    Input, integer OPTION, the option index.
c
c    Output, character*(*) TITLE, the title of the problem.
c    TITLE will never be longer than 80 characters.
c
      implicit none

      integer option
      character*(*) title

      if ( option .eq. 1 ) then
        title = 'Cell kinetics problem, seeking limit points.'
      else if ( option .eq. 2 ) then
        title = 'Cell kinetics problem, no search for limit points.'
      else
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'P07_TITLE - Fatal error!'
        write ( *, '(a)' ) '  Unrecognized option number.'
        stop
      end if

      return
      end
      subroutine p08_data (DIRIPC,HFACT,HMAX,HMIN,HTAN,option,
     & IPC,IPRAM,ISIZE,IT,LIM,MAXCON,
     & MAXLIM,MAXSTP,MAXTAR,NCOL,NROW,NVAR,RELERR,XIT,RWORK)

c*********************************************************************72
c
cc P08_DATA sets parameters for Riks's mechanical function.
c
c  8.1 riks mechanical problem.
c
c  8.2 reference
c
c  e riks,
c  the application of newtons method to the problem of elastic
c    stability,
c  trans asme, journal of appl mech, december 1972, pages 1060-1065.
c
c  8.3  the function
c
c  the equations describe the equilibrium state of a structure
c  made of three springs with a common moveable endpoint
c  and the other endpoints fixed.  a load is applied to the common
c  endpoint.
c
c  x1, x2, and x3 are the x, y, and z coordinates of the common point.
c  x4 is the magnitude of the load which is applied in the x direction.
c  if c(i) are spring constants, and a(i) are the three fixed
c  endpoints, then the (vector) equation may be written-
c
c  f(j)=sum(i=1,3) coef(i)*(a(i)-x(j)) + p(j)
c
c  where coef(i)=c(i)*(norm(a(i)-norm(x-a(i)))/norm(x-a(i))
c
c  and p=(x4,x5,x6), the norm is the euclidean norm,
c  and c(1)+c(2)+c(3)=1.0
c
c  two augmenting equations restrict the load vector p.
c
c  fx(4)=x(ifix1)-val1.
c  fx(5)=x(ifix2)-val2.
c
c  8.4 options
c
c  option=1 lim=6, it=0, ifix1=4, val1=0.0, ifix2=5, val2=0.0.
c
c  8.6 limit points
c
c  in riks paper, there seem to be limit points in x6 for
c  x6=4.10 and  x6=-3.84.
c
c  8.8 comments
c
c  current run seems uninteresting.
c  also x7, x8, x9 = origin coordinates.
c
      DIMENSION RWORK(ISIZE)
      COMMON /AUXMEM/ AVAL(3,3),CVAL(3),IFIX1,VAL1,IFIX2,VAL2

      NVAR=6
      IFREE=6
      IFIX1=4
      IFIX2=5
      DIRIPC=1.0
      HFACT=2.0
      IPC=IFREE
      IT=IFREE
      LIM=IFREE
      MAXCON=20
      MAXLIM=2
      MAXSTP=20
      MAXTAR=1
      NCOL=NVAR
      NROW=NVAR
      RELERR=1.0E-5
      XIT=10.0
      VAL1=0.0
      VAL2=0.0
      AVAL(1,1)=2.0
      AVAL(1,2)=0.0
      AVAL(1,3)=0.0
      AVAL(2,1)=-1.0
      AVAL(2,2)=1.0
      AVAL(2,3)=0.0
      AVAL(3,1)=-1.0
      AVAL(3,2)=-2.0
      AVAL(3,3)=1.0
      CVAL(1)=10.0/21.0
      CVAL(2)=6.0/21.0
      CVAL(3)=5.0/21.0
      return
      END
      subroutine p08_fun ( option, nvar, x, fx )

c*********************************************************************72
c
cc P08_FUN evaluates Riks's mechanical function.
c
      integer nvar

      double precision fx(nvar-1)
      integer option
      double precision x(nvar)

      COMMON /AUXMEM/ AVAL(3,3),CVAL(3),IFIX1,VAL1,IFIX2,VAL2

      DO I=1,3

        IP3=I+3
        FX(I)=0.0

        DO J=1,3

          ANRM=0.0
          XMANRM=0.0
          DO K=1,3
            ANRM=ANRM+AVAL(J,K)**2
            XMANRM=XMANRM+(X(K)-AVAL(J,K))**2
          end do
          ANRM=SQRT(ANRM)
          XMANRM=SQRT(XMANRM)
          FX(I)=FX(I)+CVAL(J)*(1.0-ANRM/XMANRM)*(X(I)-AVAL(J,I))
        end do

        FX(I)=FX(I)+X(IP3)

      end do

      FX(4)=X(IFIX1)-VAL1
      FX(5)=X(IFIX2)-VAL2
      return
      END
      subroutine p08_jac ( option, nvar, x, jac )

c*********************************************************************72
c
cc P08_JAC evaluates the Riks's mechanical jacobian.
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
c    Input, integer OPTION, the option index.
c
c    Input, integer NVAR, the number of variables.
c
c    Input, double precision X(NVAR), the argument of the jacobian.
c
c    Output, double precision JAC(NVAR,NVAR), the jacobian matrix evaluated
c    at X.  The NVAR-th row is not set by this routine.
c
      implicit none

      integer nvar

      integer option
      double precision jac(nvar,nvar)
      double precision x(nvar)

      COMMON /AUXMEM/ AVAL(3,3),CVAL(3),IFIX1,VAL1,IFIX2,VAL2

      do j = 1, ncol
        do i = 1, nrow
          jac(i,j) = 0.0E+00
        end do
      end do

      DO I=1,3
        DO J=1,3
          DO K=1,3
            ANRM=0.0
            XMANRM=0.0
            DO L=1,3
              ANRM=ANRM+AVAL(K,L)**2
              XMANRM=XMANRM+(X(L)-AVAL(K,L))**2
            end do
            ANRM=SQRT(ANRM)
            XMANRM=SQRT(XMANRM)
            TERM=CVAL(K)*ANRM*(X(I)-AVAL(K,I))
     &              *(X(J)-AVAL(K,J))/(XMANRM**3)
            jac(I,J)=jac(I,J)+TERM
            IF(I.EQ.J)jac(I,J)=jac(I,J)-CVAL(K)*ANRM/XMANRM
          end do

        end do

        jac(I,I)=jac(I,I)+1.0

      end do

      jac(1,4)=1.0
      jac(2,5)=1.0
      jac(3,6)=1.0
      jac(4,IFIX1)=1.0
      jac(5,IFIX2)=1.0
      jac(2,4)=0.0
      jac(3,4)=0.0
      return
      END
      subroutine p08_nvar ( option, nvar )

c*********************************************************************72
c
cc P08_NVAR sets the number of variables for problem 8.
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
c    Input, integer OPTION, the option chosen for this problem.
c    For some problems, several options are available.  At least,
c    OPTION = 1 is always legal.
c
c    Output, integer NVAR, the number of variables.
c
      implicit none

      integer nvar
      integer option

      nvar = 6

      return
      end
      subroutine p08_option_num ( option_num )

c*********************************************************************72
c
cc P08_OPTION_NUM returns the number of options for problem 8.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    16 September 2008
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, integer OPTION_NUM, the number of options.
c
      implicit none

      integer option_num

      option_num = 1

      return
      end
      subroutine p08_start ( option, nvar, x )

c*********************************************************************72
c
cc P08_START returns a starting point for problem 8.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    19 September 2008
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer OPTION, the option index.
c
c    Input, integer NVAR, the number of variables.
c
c    Output, double precision X(NVAR), the starting point.
c
      implicit none

      integer nvar

      integer option
      double precision x(nvar)

      if ( option .eq. 1 ) then

        do i = 1, nvar
          x(i) = 0.0
        end do

      else

        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'P08_START - Fatal error!'
        write ( *, '(a,i8)' ) '  Unrecognized option  = ', option
        stop

      end if

      return
      end
      subroutine p08_stepsize ( option, h, hmin, hmax )

c*********************************************************************72
c
cc P08_STEPSIZE returns step sizes for problem 8.
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
c    Input, integer OPTION, the option index.
c
c    Output, double precision H, HMIN, HMAX, suggested values for the
c    initial step, the minimum step, and the maximum step.
c
      implicit none

      double precision h
      double precision hmax
      double precision hmin
      integer option

      h =    1.000D+00
      hmin = 0.001D+00
      hmax = 1.000D+00

      return
      end
      subroutine p08_title ( option, title )

c*********************************************************************72
c
cc P08_TITLE sets the title for problem 8.
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
c    Input, integer OPTION, the option index.
c
c    Output, character*(*) TITLE, the title of the problem.
c    TITLE will never be longer than 80 characters.
c
      implicit none

      integer option
      character*(*) title

      if ( option .eq. 1 ) then
        title = 'Riks mechanical function, seeking limit points.'
      else
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'P08_TITLE - Fatal error!'
        write ( *, '(a)' ) '  Unrecognized option number.'
        stop
      end if

      return
      end
      subroutine p09_data (DIRIPC,HFACT,HMAX,HMIN,HTAN,option,
     & IPC,IPRAM,ISIZE,IT,LIM,MAXCON,
     & MAXLIM,MAXSTP,MAXTAR,NCOL,NROW,NVAR,RELERR,XIT,RWORK)

c*********************************************************************72
c
cc P09_DATA sets parameters for Oden's mechanical function.
c
c  9.1 oden mechanical problem.
c
c  9.2 reference
c
c  j t oden,
c  finite elements of nonlinear continua,
c  mcgraw hill, new york, 1972.
c
c  9.3 the function
c
c  the equations describe the equilibrium of a simple two bar
c  framework, with one common endpoint, and the other endpoints
c  fixed.  a load is applied to the common endpoint.
c  the bars are constructed of an isotropic hookean material.
c
c  the function is of the form
c
c  f1= x1**3 - 3*height*x1**2 + 2*height**2*x1
c     +x1*x2**2 - height*x2**2 - x3*cos(x4)
c
c  f2 = x1**2*x2 - 2*height*x1*x2 + x2**3 + 2*x2
c      -x3*sin(x4)
c
c  f3= x(ifix1)-val1
c
c  with height=2.0.
c
c  9.4 options
c
c  option=1  ifix1=4, val1=0.00, it=1, xit=4.0, lim=1
c  option=2  ifix1=4, val1=0.25, it=1, xit=4.0, lim=1
c  option=3  ifix1=4, val1=0.50, it=1, xit=4.0, lim=1
c  option=4  ifix1=4, val1=1.00, it=1, xit=4.0, lim=1
c  option=5  ifix1=4, val1=0.00, it=1, xit=4.0, lim=2
c            (uninteresting, no limit points for x(2)).
c  option=6  ifix1=4, val1=0.25, it=1, xit=4.0, lim=2
c  option=7  ifix1=4, val1=0.50, it=1, xit=4.0, lim=2
c  option=8  ifix1=4, val1=1.00, it=1, xit=4.0, lim=2
c  option=9  ifix1=4, val1=0.00, it=1, xit=4.0, lim=3
c  option=10 ifix1=4, val1=0.25, it=1, xit=4.0, lim=3
c  option=11 ifix1=4, val1=0.50, it=1, xit=4.0, lim=3
c  option=12 ifix1=4, val1=1.00, it=1, xit=4.0, lim=3
c  option=13 ifix1=4, val1=0.00, it=0, lim=0. (check for bifurcation).
c
c  9.5 target point
c
c  for option=1, 5, or 9, the target point is (4,0,0,0).
c
c  9.6 limit points
c
c  for option=9 there are the following limit points in x(3)
c
c  (2+-2/sqrt(3), 0, +-16/sqrt(27), 0)
c
c  for skew loads (x4 nonzero) there are various limit points.
c
c  melhem lists,
c
c  (0.5903206, 0.8391448, 0.9581753, 1.252346)
c  (2.705446,  0.6177675, 0.9581753, 1.252346)
c
c  with x3,x4 corresponding to a load vector of (.30,.91).
c
c  computational results with this program are-
c
c  option=2  limit points in x(1)
c
c  2.816913  0.7396444  -2.348587  0.2500000
c  1.183087 -0.7396445   2.348587  0.2500000
c
c  option=3  limit points in x(1)
c
c  2.520900  0.8598542  -1.774344  0.5000000
c  1.479100 -0.8598521   1.774346  0.5000000
c
c  option=4  limit points in x(1)
c
c  2.210747  0.9241686  -1.209751  1.0000000
c  (limit point finder failed at second limit point)
c
c  option=5  no limit points in x(2)
c
c  option=6  limit points in x(2)
c
c  1.831179  1.424861  0.3392428  0.2500000
c  (apparently did not reach second limit point)
c
c  option=7  limit points in x(2)
c
c  1.697061  1.453503  0.6198216  0.2500000
c  2.302939 -1.453503 -0.6198219  0.2500000
c
c  option=8  limit points in x(2)
c
c  1.534293  1.555364  1.175649  1.0000000
c  2.465706 -1.555364 -1.175648  1.0000000
c
c  option=9  limit points in x(3)
c
c  0.8452995  0.0000000  3.079199  0.0000000
c  3.154701   0.0000000 -3.079197  0.0000000
c
c  option=10  limit points in x(3)
c
c  0.5800046  0.7846684  2.004746  0.2500000
c  2.777765   0.5695726 -2.464886  0.2500000
c
c  option=11  limit points in x(3)
c
c  0.6305253  0.9921379  1.779294  0.5000000
c  2.501894   0.7202593 -1.846869  0.5000000
c
c  option=12  limit points in x(3)
c
c  0.7650624  1.292679   1.837450  1.000000
c  2.204188   0.8010838 -1.253382  1.000000
c
c  9.7 bifurcation points
c
c  for option=13
c
c  (2+-sqrt(2), 0, +-sqrt(2), 0)
c
c  9.8 comments
c
c  program ok.
c  rerun option=4 and option=6 to get second limit point.
c
      DIMENSION RWORK(ISIZE)
      COMMON /AUXMEM/ HEIGHT,IFIX1,VAL1

      NVAR=4
      HEIGHT=2.0
      IF(option.EQ.1.OR.option.EQ.5.OR.option.EQ.9)VAL1=0.00
      IF(option.EQ.13)VAL1=0.00
      IF(option.EQ.2.OR.option.EQ.6.OR.option.EQ.10)VAL1=0.25
      IF(option.EQ.3.OR.option.EQ.7.OR.option.EQ.11)VAL1=0.50
      IF(option.EQ.4.OR.option.EQ.8.OR.option.EQ.12)VAL1=1.00

      DIRIPC=1.0
      HFACT=3.0
      IPC=1
      IF(option.NE.13)IT=1
      IF(option.EQ.13)IT=0
      IF(option.GE.1.AND.option.LE.4)LIM=1
      IF(option.GE.5.AND.option.LE.8)LIM=2
      IF(option.GE.9.AND.option.LE.12)LIM=3
      IF(option.EQ.13)LIM=0
      MAXCON=30
      IF(option.NE.13)MAXLIM=2
      IF(option.EQ.13)MAXLIM=30
      MAXSTP=30
      IF(option.NE.13)MAXTAR=1
      IF(option.EQ.13)MAXTAR=30
      NCOL=NVAR
      NROW=NVAR
      RELERR=1.0E-5
      XIT=4.0

      return
      END
      subroutine p09_fun ( option, nvar, x, fx )

c*********************************************************************72
c
cc P09_FUN evaluates Oden's mechanical function.
c
      integer nvar

      double precision fx(nvar-1)
      integer option
      double precision x(nvar)

      COMMON /AUXMEM/ HEIGHT,IFIX1,VAL1

      X1=X(1)
      X2=X(2)
      X3=X(3)
      X4=X(4)
      FX(1)=X1*X1*X1-3.0*HEIGHT*X1*X1+2.0*HEIGHT*HEIGHT*X1
      FX(1)=FX(1)+X1*X2*X2-HEIGHT*X2*X2-X3*COS(X4)
      FX(2)=X1*X1*X2-2.0*HEIGHT*X1*X2+X2*X2*X2+2.0*X2-X3*SIN(X4)
      FX(3)=X(IFIX1)-VAL1
      return
      END
      subroutine p09_jac ( option, nvar, x, jac )

c*********************************************************************72
c
cc P09_JAC evaluates Oden's mechanical jacobian.
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
c    Input, integer OPTION, the option index.
c
c    Input, integer NVAR, the number of variables.
c
c    Input, double precision X(NVAR), the argument of the jacobian.
c
c    Output, double precision JAC(NVAR,NVAR), the jacobian matrix evaluated
c    at X.  The NVAR-th row is not set by this routine.
c
      implicit none

      integer nvar

      integer option
      double precision jac(nvar,nvar)
      double precision x(nvar)

      COMMON /AUXMEM/ HEIGHT,IFIX1,VAL1

      X1=X(1)
      X2=X(2)
      X3=X(3)
      X4=X(4)
      jac(1,1)=3.0*X1*X1-6.0*HEIGHT*X1+2.0*HEIGHT*HEIGHT+X2*X2
      jac(1,2)=2.0*X1*X2-2.0*HEIGHT*X2
      jac(1,3)=-COS(X4)
      jac(1,4)=X3*SIN(X4)
      jac(2,1)=2.0*X1*X2-2.0*HEIGHT*X2
      jac(2,2)=X1*X1-2.0*HEIGHT*X1+3.0*X2*X2+2.0
      jac(2,3)=-SIN(X4)
      jac(2,4)=-X3*COS(X4)
      jac(3,1)=0.0
      jac(3,2)=0.0
      jac(3,3)=0.0
      jac(3,4)=0.0
      jac(3,IFIX1)=1.0
      return
      END
      subroutine p09_nvar ( option, nvar )

c*********************************************************************72
c
cc P09_NVAR sets the number of variables for problem 9.
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
c    Input, integer OPTION, the option chosen for this problem.
c    For some problems, several options are available.  At least,
c    OPTION = 1 is always legal.
c
c    Output, integer NVAR, the number of variables.
c
      implicit none

      integer nvar
      integer option

      nvar = 4

      return
      end
      subroutine p09_option_num ( option_num )

c*********************************************************************72
c
cc P09_OPTION_NUM returns the number of options for problem 9.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    16 September 2008
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, integer OPTION_NUM, the number of options.
c
      implicit none

      integer option_num

      option_num = 13

      return
      end
      subroutine p09_start ( option, nvar, x )

c*********************************************************************72
c
cc P09_START returns a starting point for problem 9.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    19 September 2008
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer OPTION, the option index.
c
c    Input, integer NVAR, the number of variables.
c
c    Output, double precision X(NVAR), the starting point.
c
      implicit none

      integer nvar

      integer i
      integer option
      double precision x(nvar)

      do i = 1, nvar - 1
        x(i) = 0.0
      end do

      if ( mod ( option, 4 ) .eq. 1 ) then
        x(4) = 0.0
      else if ( mod ( option, 4 ) .eq. 2 ) then
        x(4) = 0.25
      else if ( mod ( option, 4 ) .eq. 3 ) then
        x(4) = 0.50
      else if ( mod ( option, 4 ) .eq. 4 ) then
        x(4) = 1.00
      else

        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'P09_START - Fatal error!'
        write ( *, '(a,i8)' ) '  Unrecognized option  = ', option
        stop

      end if

      return
      end
      subroutine p09_stepsize ( option, h, hmin, hmax )

c*********************************************************************72
c
cc P09_STEPSIZE returns step sizes for problem 9.
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
c    Input, integer OPTION, the option index.
c
c    Output, double precision H, HMIN, HMAX, suggested values for the
c    initial step, the minimum step, and the maximum step.
c
      implicit none

      double precision h
      double precision hmax
      double precision hmin
      integer option

      h =    0.300D+00
      hmin = 0.001D+00
      hmax = 0.600D+00

      return
      end
      subroutine p09_title ( option, title )

c*********************************************************************72
c
cc P09_TITLE sets the title for problem 9.
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
c    Input, integer OPTION, the option index.
c
c    Output, character*(*) TITLE, the title of the problem.
c    TITLE will never be longer than 80 characters.
c
      implicit none

      integer option
      character*(*) title

      if ( option .eq. 1 ) then
        title =
     &    'Oden problem, VAL=0.00, Target X(1)=4.0, Limits in X(1).'
      else if ( option .eq. 2 ) then
        title =
     &    'Oden problem, VAL=0.25, Target X(1)=4.0, Limits in X(1).'
      else if ( option .eq. 3 ) then
        title =
     &    'Oden problem, VAL=0.50, Target X(1)=4.0, Limits in X(1).'
      else if ( option .eq. 4 ) then
        title =
     &    'Oden problem, VAL=1.00, Target X(1)=4.0, Limits in X(1).'
      else if ( option .eq. 5 ) then
        title =
     &    'Oden problem, VAL=0.00, Target X(1)=4.0, Limits in X(2).'
      else if ( option .eq. 6 ) then
        title =
     &    'Oden problem, VAL=0.25, Target X(1)=4.0, Limits in X(2).'
      else if ( option .eq. 7 ) then
        title =
     &    'Oden problem, VAL=0.50, Target X(1)=4.0, Limits in X(2).'
      else if ( option .eq. 8 ) then
        title =
     &    'Oden problem, VAL=1.00, Target X(1)=4.0, Limits in X(2).'
      else if ( option .eq. 9 ) then
        title =
     &    'Oden problem, VAL=0.00, Target X(1)=4.0, Limits in X(3).'
      else if ( option .eq. 10 ) then
        title =
     &    'Oden problem, VAL=0.25, Target X(1)=4.0, Limits in X(3).'
      else if ( option .eq. 11 ) then
        title =
     &    'Oden problem, VAL=0.50, Target X(1)=4.0, Limits in X(3).'
      else if ( option .eq. 12 ) then
        title =
     &    'Oden problem, VAL=1.00, Target X(1)=4.0, Limits in X(3).'
      else if ( option .eq. 13 ) then
        title =
     &    'Oden problem, VAL=0.00, no targets, no limits.'
      else
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'P09_TITLE - Fatal error!'
        write ( *, '(a)' ) '  Unrecognized option number.'
        stop
      end if

      return
      end
      subroutine p10_data (DIRIPC,HFACT,HMAX,HMIN,HTAN,option,
     & IPC,IPRAM,ISIZE,IT,LIM,MAXCON,
     & MAXLIM,MAXSTP,MAXTAR,NCOL,NROW,NVAR,RELERR,XIT,RWORK)

c*********************************************************************72
c
cc P10_DATA sets parameters for the finite difference rod torsion function.
c
c  The rod torsion problem, finite difference version
c
c  Reference:
c
c  werner c rheinboldt,
c  on the solution of some nonlinear equations arising in
c    the application of finite element methods,
c  in mathematics of finite element methods,
c  edited by j whiteman
c  academic press, london, 1976, pages 465-482.
c
c  The function:
c
c  the problem is a boundary value problem on (0,1) x (0,1).
c
c  if del(u)=d2 u/dx2 + d2 u/dy2, then the problem may be written
c
c  -d/dx ( gamma(del u) du/dx) - d/dy ( gamma(del u)* du/dy) = x(nvar)
c
c  a standard finite difference approximation on a uniform
c  mesh is applied to yield the equations.
c
c  Options:
c
c  option=1  it=nvar, xit=5.0, gamma(s)=exp(5*s)
c
c  option=2  it=nvar, xit=15.0,
c
c  gamma(s)=  a                for s.le.0.15
c          =  .5*(a+b)+.25*(b-a)*(3/7*(40*s-13)-(40*s-13)**3/343)
c                              for 0.15.le.s.le.0.5
c          =  b                for 0.5.le.s
c  with a=1.0, b=10.0.
c
c  Comments:
c
c  note that writeup in rheinboldt omits second gamma
c    in equation.
c  program seems ok.
c
      COMMON /AUXMEM/ IROW,ICOL,H,IFUN12
      DIMENSION RWORK(ISIZE)

      IROW=6
      ICOL=6
      NVAR=IROW*ICOL+1
      DIRIPC=1.0
      HFACT=3.0
      IPC=NVAR
      IT=NVAR
      LIM=0
      MAXCON=30
      MAXLIM=30
      MAXSTP=30
      MAXTAR=1
      NCOL=NVAR
      NROW=NVAR
      RELERR=1.0E-5
      IF(option.EQ.1)XIT=5.0
      IF(option.EQ.2)XIT=15.0
c
c  set common parameters
c
      H=1.0/FLOAT(IROW+1)
      IF(option.EQ.1)IFUN12=1
      IF(option.EQ.2)IFUN12=2
      return
      END
      subroutine p10_fun ( option, nvar, x, fx )

c*********************************************************************72
c
cc P10_FUN evaluates the finite difference rod torsion function.
c
      integer nvar

      double precision fx(nvar-1)
      integer option
      double precision x(nvar)

      DIMENSION UC(2),UX(4),CX(2,4)
      COMMON /AUXMEM/ IROW,ICOL,H,IFUN12
c
c  uc contains the two cornerpoints
c  ux contains the four side-points
c  cx contains the elements connected to the side points
c
      NEQN=NVAR-1
      H2=H*H

      DO I=1,IROW

        DO 70 J=1,ICOL
          UC(1)=0.0
          UC(2)=0.0
          DO K=1,4
            UX(K)=0.0
          end do
          IJ=(J-1)*IROW+I
          INDEX=IJ-1
          IF (I.NE.1) UX(1)=X(INDEX)
          INDEX=IJ+1
          IF (I.NE.IROW) UX(2)=X(INDEX)
          IF (J.EQ.1) GO TO 20
          JK=IJ-IROW
          UX(3)=X(JK)
          INDEX=JK-1
          IF (I.NE.1) UC(1)=X(INDEX)
          IF (J.EQ.ICOL) GO TO 30
   20     JK=IJ+IROW
          UX(4)=X(JK)
          INDEX=JK+1
          IF (I.NE.IROW) UC(2)=X(INDEX)
c
c  k=1, 2*qw calculated and stored in ( cx(1,k) + cx(2,k) )
c  k=2, 2*qe calculated and stored in ( cx(1,k) + cx(2,k) )
c  k=3, 2*qs calculated and stored in ( cx(1,k) + cx(2,k) )
c  k=4, 2*qn calculated and stored in ( cx(1,k) + cx(2,k) )
c
   30     SC=0.0

          DO K=1,4
            IF (K.EQ.1.OR.K.EQ.3) K1=1
            IF (K.EQ.2.OR.K.EQ.4) K1=2
            K2=5-K
            DO L=1,2
              RLK=(UX(K)-X(IJ))**2
              IF (L.EQ.1) RLK=(RLK+(UX(K)-UC(K1))**2)/H2
              IF (L.EQ.2) RLK=(RLK+(X(IJ)-UX(K2))**2)/H2
              CX(L,K)= p10_gx (RLK)
              SC=SC+CX(L,K)
            end do

          end do
c
c  sc=qn + qs + qe + qw
c
          FX(IJ)=.5*SC*X(IJ)-X(NVAR)*H2

          DO K=1,4
            FX(IJ)=FX(IJ)-.5*UX(K)*(CX(1,K)+CX(2,K))
          end do

   70   CONTINUE

      end do

      return
      END
      FUNCTION p10_gx (S)

c*********************************************************************72
c
cc P10_GX is an auxilliary function for the rod torsion problem.
c
      COMMON /AUXMEM/ IROW,ICOL,H,IFUN12

      aval = 1.0
      bval = 10.0
      s1val = 0.15
      s2val = 0.50

      IF ( IFUN12 .EQ. 1 ) then
        GX0012=EXP(5.0*S)
      else
        IF (S.GT.S1VAL) GO TO 20
        GX0012=AVAL
        return
   20   IF (S.LT.S2VAL) GO TO 30
        GX0012=BVAL
        return
   30   R=(S-0.5*(S1VAL+S2VAL))/(0.5*(S2VAL-S1VAL))
        GX0012=.5*(AVAL+BVAL)+.25*(BVAL-AVAL)*R*(3.0-R**2)
      end if

      return
      END
      subroutine p10_jac ( option, nvar, x, jac )

c*********************************************************************72
c
cc P10_JAC evaluates the finite difference rod torsion jacobian.
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
c    Input, integer OPTION, the option index.
c
c    Input, integer NVAR, the number of variables.
c
c    Input, double precision X(NVAR), the argument of the jacobian.
c
c    Output, double precision JAC(NVAR,NVAR), the jacobian matrix evaluated
c    at X.  The NVAR-th row is not set by this routine.
c
      implicit none

      integer nvar

      integer option
      double precision jac(nvar,nvar)
      double precision x(nvar)

      COMMON /AUXMEM/ IROW,ICOL,H,IFUN12

      DIMENSION UC(2),UX(4),CX(2,4),DX(2,4)
      NEQN=NVAR-1

      do j = 1, ncol
        do i = 1, nrow
          jac(i,j) = 0.0E+00
        end do
      end do

      H2=H*H
      DO 150 I=1,IROW
        DO 150 J=1,ICOL
          UC(1)=0.0
          UC(2)=0.0
          DO K=1,4
            UX(K)=0.0
          end do
          IJ=(J-1)*IROW+I
          IJM1=IJ-1
          IF (I.NE.1) UX(1)=X(IJM1)
          IJP1=IJ+1
          IF (I.NE.IROW) UX(2)=X(IJP1)
          IF (J.EQ.1) GO TO 30
          JK=IJ-IROW
          UX(3)=X(JK)
          JKM1=JK-1
          IF (I.NE.1) UC(1)=X(JKM1)
          IF (J.EQ.ICOL) GO TO 40
   30     JK=IJ+IROW
          UX(4)=X(JK)
          JKP1=JK+1
          IF (I.NE.IROW) UC(2)=X(JKP1)
   40     SC=0.0

          DO K=1,4

            IF (K.EQ.1.OR.K.EQ.3) K1=1
            IF (K.EQ.2.OR.K.EQ.4) K1=2
            K2=5-K

            DO L=1,2
              RLK=(UX(K)-X(IJ))**2
              IF (L.EQ.1) RLK=(RLK+(UX(K)-UC(K1))**2)/H2
              IF (L.EQ.2) RLK=(RLK+(X(IJ)-UX(K2))**2)/H2
              CX(L,K)=p10_gx ( RLK )
              DX(L,K)=p10_gp ( RLK )
              SC=SC+CX(L,K)
            end do

          end do
c
c  diagonal
c
          XJAC=0.5*SC
          DO K=1,4
            K2=5-K
            XJAC=XJAC+DX(2,K)*(X(IJ)-UX(K))*(2.0*X(IJ)-UX(K)-UX(K2))/H2
            XJAC=XJAC+DX(1,K)*(X(IJ)-UX(K))**2/H2
          end do
          jac(IJ,IJ)=XJAC
c
c  off-diagonals
c
          DO 130 K=1,4
            GO TO (80,90,100,110),K
   80       IF (I.EQ.1) GO TO 130
            JK=IJ-1
            GO TO 120
   90       IF (I.EQ.IROW) GO TO 130
            JK=IJ+1
            GO TO 120
  100       IF (J.EQ.1) GO TO 130
            JK=IJ-IROW
            GO TO 120
  110       IF (J.EQ.ICOL) GO TO 130
            JK=IJ+IROW
  120       IF (K.EQ.1.OR.K.EQ.3) K1=1
            IF (K.EQ.2.OR.K.EQ.4) K1=2
            K2=5-K
            XJAC=(X(IJ)-UX(K))*(DX(1,K)*(2.0*UX(K)-X(IJ)-UC(K1))
     &               +DX(2,K)*(UX(K)-X(IJ))+DX(2,K2)*(UX(K2)-X(IJ)))
            XJAC=XJAC/H2-.5*(CX(1,K)+CX(2,K))
            jac(IJ,JK)=XJAC
  130       CONTINUE
          IF ((I.EQ.1).OR.(J.EQ.1)) GO TO 140
          JK=IJ-IROW-1
          XJAC=(X(IJ)-UX(1))*DX(1,1)*(UC(1)-UX(1))
     &         +(X(IJ)-UX(3))*DX(1,3)*(UC(1)-UX(3))
          XJAC=XJAC/H2
          jac(IJ,JK)=XJAC
  140     IF (I.EQ.IROW.OR.J.EQ.ICOL) GO TO 150
          JK=IJ+IROW+1
          XJAC=(X(IJ)-UX(2))*DX(1,2)*(UC(2)-UX(2))
     &         +(X(IJ)-UX(4))*DX(1,4)*(UC(2)-UX(4))
          XJAC=XJAC/H2
          jac(IJ,JK)=XJAC
  150     CONTINUE

      DO I=1,NEQN
        jac(I,NVAR)=-H2
      end do

      return
      END
      FUNCTION p10_gp (S)

c*********************************************************************72
c
cc P10_GP is an auxilliary function for the rod torsion problem.
c
      COMMON /AUXMEM/ IROW,ICOL,H,IFUN12

      aval = 1.0
      bval = 10.0
      s1val = 0.15
      s2val = 0.50

      IF(IFUN12.EQ.2)GO TO 10
      GP0012=5.0*EXP(5.0*S)
      return
   10 IF(S.GT.S1VAL.AND.S.LT.S2VAL)GO TO 20
      GP0012=0.0
      return
   20 R=(S-0.5*(S1VAL+S2VAL))/(0.5*(S2VAL-S1VAL))
      GP0012=0.75*(BVAL-AVAL)*(1.0-R*R)/(0.5*(S2VAL-S1VAL))
      return
      END
      subroutine p10_nvar ( option, nvar )

c*********************************************************************72
c
cc P10_NVAR sets the number of variables for the finite difference rod torsion function.
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
c    Input, integer OPTION, the option chosen for this problem.
c    For some problems, several options are available.  At least,
c    OPTION = 1 is always legal.
c
c    Output, integer NVAR, the number of variables.
c
      implicit none

      integer nvar
      integer option

      nvar = 37

      return
      end
      subroutine p10_option_num ( option_num )

c*********************************************************************72
c
cc P10_OPTION_NUM returns the number of options for the finite difference rod torsion function.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    16 September 2008
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, integer OPTION_NUM, the number of options.
c
      implicit none

      integer option_num

      option_num = 2

      return
      end
      subroutine p10_start ( option, nvar, x )

c*********************************************************************72
c
cc P10_START returns a starting point for the finite difference rod torsion function.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    19 September 2008
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer OPTION, the option index.
c
c    Input, integer NVAR, the number of variables.
c
c    Output, double precision X(NVAR), the starting point.
c
      implicit none

      integer nvar

      integer i
      integer option
      double precision x(nvar)

      if ( 1 .le. option .and. option .le. 2 ) then

        do i = 1, nvar
          x(i) = 0.0
        end do

      else

        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'P10_START - Fatal error!'
        write ( *, '(a,i8)' ) '  Unrecognized option  = ', option
        stop

      end if

      return
      end
      subroutine p10_stepsize ( option, h, hmin, hmax )

c*********************************************************************72
c
cc P10_STEPSIZE returns step sizes for the finite difference rod torsion function.
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
c    Input, integer OPTION, the option index.
c
c    Output, double precision H, HMIN, HMAX, suggested values for the
c    initial step, the minimum step, and the maximum step.
c
      implicit none

      double precision h
      double precision hmax
      double precision hmin
      integer option

      h =     2.000D+00
      hmin =  0.001D+00
      hmax = 10.000D+00

      return
      end
      subroutine p10_title ( option, title )

c*********************************************************************72
c
cc P10_TITLE sets the title for the finite difference rod torsion function.
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
c    Input, integer OPTION, the option index.
c
c    Output, character*(*) TITLE, the title of the problem.
c    TITLE will never be longer than 80 characters.
c
      implicit none

      integer option
      character*(*) title

      if ( option .eq. 1 ) then
        title = 'Finite difference rod torsion, exponential PHI'
      else if ( option .eq. 2 ) then
        title = 'Finite difference rod torsion, piecewise PHI.'
      else
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'P10_TITLE - Fatal error!'
        write ( *, '(a)' ) '  Unrecognized option number.'
        stop
      end if

      return
      end
      subroutine p11_data (DIRIPC,HFACT,HMAX,HMIN,HTAN,option,
     & IPC,IPRAM,ISIZE,IT,LIM,MAXCON,
     & MAXLIM,MAXSTP,MAXTAR,NCOL,NROW,NVAR,RELERR,XIT,XR)

c*********************************************************************72
c
cc P11_DATA sets parameters for the finite element rod torsion function.
c
c  Torsion of a square rod, finite element
c
c  Reference
c
c  werner c. rheinboldt,
c  on the solution of some nonlinear equations arising in the
c    application of finite element methods,
c  mathematics of finite element methods,
c  edited by j. whiteman,
c  academic press, london, 1976.
c  pages 465-482.
c
c  The problem formulation
c
c  on the 2-dimensional region [0,1] x [0,1], the following
c  square grid is used:
c
c  if k = j * nside + i + 1,
c  and h=1/(nside-1),
c  point k is (i*h, j*h).
c
c  thus, for nside=5,
c
c  k  i  j  x  y
c
c  1  0  0  0  0
c  2  0  1  h  0
c  3  0  2 2h  0
c  4  0  3 3h  0
c  5  0  4 4h  0
c  6  1  0  0  h
c  7  1  1  h  h
c  8  1  2 2h  h
c  9  1  3 3h  h
c 10  1  4 4h  h
c 11  2  0  0 2h
c
c ...
c
c 24  4  3 3h 4h
c 25  4  4 4h 4h.
c
c  21---22---23---24---25
c   !    !    !    !    !
c   ! 13 ! 14 ! 15 ! 16 !
c   !    !    !    !    !
c  16---17---18---19---20
c   !    !    !    !    !
c   ! 09 ! 10 ! 11 ! 12 !
c   !    !    !    !    !
c  11---12---13---14---15
c   !    !    !    !    !
c   ! 05 ! 06 ! 07 ! 08 !
c   !    !    !    !    !
c  06---07---08---09---10
c
c  on a single element, the ordering of nodes and shape functions is
c
c  03---04
c   !    !
c   ! 01 !
c   !    !
c  01---02
c
c  thus, if h is the length of a side, the first shape function
c  is  psi 1 (x,y) = (x-xright)*(y-ytop)/h**2, and so on.
c
c   !    !    !    !    !
c   ! 01 ! 02 ! 03 ! 04 !
c   !    !    !    !    !
c  01---02---03---04---05
c
c  Options:
c
c  option=1  phi1=phi2=exp(5*(dudx**2+dudy**2))
c            f(l,u)=-5*l
c  option=2  phi1=phi2
c           let s=(dudx**2+dudy**2),
c           sbar=(40*s-13)/7
c           if (s.le.0.15) phi1=1.0
c           if (0.15.le.s.le.0.50)
c             phi1=0.5*(11)*0.25*9*(3*sbar-sbar**3)
c           if (sbar.ge.0.50) phi1=10.0
c           f(l,u)=-10*l
c
      DIMENSION XR(26)
      COMMON /CM0029/ HSIDE,IFUN29,NSHAPE,NSIDE
      COMMON /BC0029/ BCVAL
      COMMON /GS0029/ NGAUSS,WGAUSS(4),XGAUSS(4),YGAUSS(4)

      NVAR=26

      DIRIPC=1.0
      HFACT=3.0
      IPC=NVAR
      IT=NVAR
      MAXSTP=10
      MAXCON=MAXSTP
      MAXLIM=MAXSTP
      MAXTAR=1
      NCOL=NVAR
      NROW=NVAR
      RELERR=1.0E-6
      XIT=2.0
c
c  set stuff
c
      NSIDE=5
      IF(option.NE.2)IFUN29=1
      IF(option.EQ.2)IFUN29=2
      BCVAL=0.0
      HSIDE=1.0/FLOAT(NSIDE-1)
      NGAUSS=4
      NSHAPE=4
      return
      END
      subroutine p11_fun ( option, nvar, x, fx )

c*********************************************************************72
c
cc P11_FUN evaluates the finite element rod torsion function.
c
      implicit none

      integer nvar

      double precision fx(nvar-1)
      integer option
      double precision x(nvar)

      DIMENSION NODE(4)
      COMMON /CM0029/ HSIDE,IFUN29,NSHAPE,NSIDE
      COMMON /BC0029/ BCVAL
      COMMON /UV0029/ UVAL,DUDX,DUDY
      COMMON /GS0029/ NGAUSS,WGAUSS(4),XGAUSS(4),YGAUSS(4)
      COMMON /PS0029/ PSI(4),DPSIDX(4),DPSIDY(4)

      DO I=1,NVAR
        FX(I)=0.0
      end do

      NELEM=(NSIDE-1)**2
      DO IELEM=1,NELEM
c
c  from element number, compute 4 node numbers
c  in order sw, se, nw, ne.
c
        IROWM1=(IELEM-1)/(NSIDE-1)
        ICOL=IELEM-IROWM1*(NSIDE-1)
        IROW=IROWM1+1
        XMID=HSIDE*FLOAT(2*ICOL-1)/2.0
        YMID=HSIDE*FLOAT(2*IROW-1)/2.0
        NODE(1)=IROWM1*NSIDE+ICOL
        NODE(2)=NODE(1)+1
        NODE(3)=NODE(1)+NSIDE
        NODE(4)=NODE(3)+1
c
c  get gauss points for this element
c
        call p11_gaus (XMID,YMID)
c
c  for each gauss point in this element, evaluate the
c  integrand
c
        DO JGAUSS=1,NGAUSS

          XVAL=XGAUSS(JGAUSS)
          YVAL=YGAUSS(JGAUSS)
c
c  evaluate shape functions
c
          call p11_sh (XMID,XVAL,YMID,YVAL)
c
c  evaluate u and derivatives
c
          call p11_uval (NODE,NSHAPE,NVAR,X)
c
c  evaluate phi
c
          call p11_phi (PHI1,PHI1P,PHI2,PHI2P,XVAL,YVAL)
c
c  evaluate f(lambda,u)
c
          call p11_flam (FLAM,FLAMPL,FLAMPU,UVAL,X(NVAR),XVAL,YVAL)
c
c  compute inner product of equation with each shape function
c  and add to appropriate function
c
          DO 30 KSHAPE=1,NSHAPE
            NOD=NODE(KSHAPE)
            JROWM1=(NOD-1)/NSIDE
            JCOL=NOD-JROWM1*NSIDE
            JROW=JROWM1+1
            IF(JROW.EQ.1.OR.JROW.EQ.NSIDE)GO TO 20
            IF(JCOL.EQ.1.OR.JCOL.EQ.NSIDE)GO TO 20
c
c  interior node
c
            TERM=PHI1*DUDX*DPSIDX(KSHAPE)
     &          +PHI2*DUDY*DPSIDY(KSHAPE)
     &          +FLAM*PSI(KSHAPE)
            TERM=WGAUSS(JGAUSS)*HSIDE*HSIDE*TERM
            FX(NOD)=FX(NOD)+TERM
            GO TO 30
c
c  boundary node
c
   20       CONTINUE

            FX(NOD)=X(NOD)-BCVAL

   30     CONTINUE

        end do

      end do

      return
      END
      subroutine p11_flam (FLAM,FLAMPL,FLAMPU,UVAL,XNVAR,XVAL,YVAL)

c*********************************************************************72
c
cc P11_FLAM is an auxilliary function for the finite element rod torsion function.
c
      COMMON /CM0029/ HSIDE,IFUN29,NSHAPE,NSIDE
      FLAM=0.0
      FLAMPU=0.0
      FLAMPL=0.0
      IF(IFUN29.EQ.1)FLAM=-5.0*XNVAR
      IF(IFUN29.EQ.1)FLAMPL=-5.0
      IF(IFUN29.EQ.2)FLAM=-10.0*XNVAR
      IF(IFUN29.EQ.2)FLAMPL=-10.0
      return
      END
      subroutine p11_phi (PHI1,PHI1P,PHI2,PHI2P,XVAL,YVAL)

c*********************************************************************72
c
cc P11_PHI is an auxilliary function for the finite element rod torsion function.
c
      COMMON /CM0029/ HSIDE,IFUN29,NSHAPE,NSIDE
      COMMON /GS0029/ NGAUSS,WGAUSS(4),XGAUSS(4),YGAUSS(4)
      COMMON /UV0029/ UVAL,DUDX,DUDY
      COMMON /PS0029/ PSI(4),DPSIDX(4),DPSIDY(4)
      S=DUDX*DUDX+DUDY*DUDY
      PHI1=0.0
      PHI1P=0.0
      PHI2=0.0
      PHI2P=0.0
      IF(IFUN29.EQ.1)GO TO 10
      IF(IFUN29.EQ.2)GO TO 20
      WRITE(6,1000)IFUN29
      STOP
   10 CONTINUE
      PHI1=EXP(5.0*S)
      PHI1P=5.0*PHI1
      PHI2=PHI1
      PHI2P=PHI1P
      return
   20 CONTINUE
      SBAR=(40.0*S-13.0)/7.0
      IF(S.LE.0.15)PHI1=1.0
      IF(S.LE.0.15)PHI1P=0.0
      IF(0.15.LE.S.AND.S.LE.0.5)
     & PHI1=0.5*(11.0)+0.25*9.0*(3.0*SBAR-SBAR**3)
      IF(0.15.LE.S.AND.S.LE.0.5)
     & PHI1P=0.25*9.0*(3.0-3.0*SBAR**2)*40.0/7.0
      IF(0.5.LE.S)PHI1=10.0
      IF(0.5.LE.S)PHI1P=0.0
      PHI2=PHI1
      PHI2P=PHI1P
      return
 1000 FORMAT(' PHI029 - ILLEGAL VALUE OF IFUN29=',I6)
      END
      subroutine p11_uval (NODE,NSHAPE,NVAR,X)

c*********************************************************************72
c
cc P11_UVAL is an auxilliary function for the finite element rod torsion function.
c
c  evaluate u, dudx, dudy at (xval,yval) in the square of
c  side hside with midpoint (xmid,ymid), with x containing the
c  coefficients of the four shape functions.
c
      DIMENSION NODE(NSHAPE)
      DIMENSION X(NVAR)
      COMMON /PS0029/ PSI(4),DPSIDX(4),DPSIDY(4)
      COMMON /UV0029/ UVAL,DUDX,DUDY
c
c  compute u, dudx, dudy
c
      UVAL=0.0
      DUDX=0.0
      DUDY=0.0

      DO I=1,NSHAPE
        UVAL=UVAL+X(NODE(I))*PSI(I)
        DUDX=DUDX+X(NODE(I))*DPSIDX(I)
        DUDY=DUDY+X(NODE(I))*DPSIDY(I)
      end do

      return
      END
      subroutine p11_sh (XMID,XVAL,YMID,YVAL)

c*********************************************************************72
c
cc P11_SH is an auxilliary function for the finite element rod torsion function.
c
c  at the point (xval, yval) in the square of side hside and center
c  (xmid, ymid), evaluate
c
c  psi, the shape functions,
c  dpsidx and dpsidy, the partial derivatives,
c
      COMMON /PS0029/ PSI(4),DPSIDX(4),DPSIDY(4)
      COMMON /CM0029/ HSIDE,IFUN29,NSHAPE,NSIDE
c
c  set coordinates
c
      XLEFT=XMID-0.5*HSIDE
      XRITE=XMID+0.5*HSIDE
      YBOT=YMID-0.5*HSIDE
      YTOP=YMID+0.5*HSIDE
c
c  evaluate shape functions
c  order is sw, se, nw, ne.
c
      PSI(1)= (XVAL-XRITE)*(YVAL-YTOP)/(HSIDE**2)
      PSI(2)=-(XVAL-XLEFT)*(YVAL-YTOP)/(HSIDE**2)
      PSI(3)=-(XVAL-XRITE)*(YVAL-YBOT)/(HSIDE**2)
      PSI(4)= (XVAL-XLEFT)*(YVAL-YBOT)/(HSIDE**2)
c
c  evaluate derivatives
c
      DPSIDX(1)= (YVAL-YTOP)/(HSIDE**2)
      DPSIDX(2)=-(YVAL-YTOP)/(HSIDE**2)
      DPSIDX(3)=-(YVAL-YBOT)/(HSIDE**2)
      DPSIDX(4)= (YVAL-YBOT)/(HSIDE**2)
      DPSIDY(1)= (XVAL-XRITE)/(HSIDE**2)
      DPSIDY(2)=-(XVAL-XLEFT)/(HSIDE**2)
      DPSIDY(3)=-(XVAL-XRITE)/(HSIDE**2)
      DPSIDY(4)= (XVAL-XLEFT)/(HSIDE**2)
      return
      END
      subroutine p11_gaus (XMID,YMID)

c*********************************************************************72
c
cc P11_GAUS returns Gauss points for the finite element rod torsion problem.
c
c  g1 = (xmid-alfa*hside, ymid-alfa*hside)
c  g2 = (xmid-alfa*hside, ymid+alfa*hside)
c  g3 = (xmid+alfa*hside, ymid-alfa*hside)
c  g4 = (xmid+alfa*hside, ymid+alfa*hside)
c
c  where (xmid, ymid) is the center point of the square of
c  side hside, and alfa=1/(2*sqrt(3))
c
      COMMON /CM0029/ HSIDE,IFUN29,NSHAPE,NSIDE
      COMMON /GS0029/ NGAUSS,WGAUSS(4),XGAUSS(4),YGAUSS(4)
      ALFA=1.0/(2.0*SQRT(3.0))
      WGAUSS(1)=0.25
      WGAUSS(2)=0.25
      WGAUSS(3)=0.25
      WGAUSS(4)=0.25
      XGAUSS(1)=XMID-ALFA*HSIDE
      XGAUSS(2)=XMID+ALFA*HSIDE
      XGAUSS(3)=XMID-ALFA*HSIDE
      XGAUSS(4)=XMID+ALFA*HSIDE
      YGAUSS(1)=YMID-ALFA*HSIDE
      YGAUSS(2)=YMID-ALFA*HSIDE
      YGAUSS(3)=YMID+ALFA*HSIDE
      YGAUSS(4)=YMID+ALFA*HSIDE
      return
      END
      subroutine p11_jac ( option, nvar, x, jac )

c*********************************************************************72
c
cc P11_JAC evaluates the finite element rod torsion jacobian.
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
c    Input, integer OPTION, the option index.
c
c    Input, integer NVAR, the number of variables.
c
c    Input, double precision X(NVAR), the argument of the jacobian.
c
c    Output, double precision JAC(NVAR,NVAR), the jacobian matrix evaluated
c    at X.  The NVAR-th row is not set by this routine.
c
      implicit none

      integer nvar

      integer option
      double precision jac(nvar,nvar)
      double precision x(nvar)

      DIMENSION NODE(4)

      COMMON /CM0029/ HSIDE,IFUN29,NSHAPE,NSIDE
      COMMON /BC0029/ BCVAL
      COMMON /UV0029/ UVAL,DUDX,DUDY
      COMMON /GS0029/ NGAUSS,WGAUSS(4),XGAUSS(4),YGAUSS(4)
      COMMON /PS0029/ PSI(4),DPSIDX(4),DPSIDY(4)

      do j = 1, ncol
        do i = 1, nrow
          jac(i,j) = 0.0E+00
        end do
      end do

      NELEM=(NSIDE-1)**2
      DO 70 IELEM=1,NELEM
c
c  from element number, compute 4 node numbers
c  in order sw, se, nw, ne.
c
        IROWM1=(IELEM-1)/(NSIDE-1)
        ICOL=IELEM-IROWM1*(NSIDE-1)
        IROW=IROWM1+1
        XMID=HSIDE*FLOAT(2*ICOL-1)/2.0
        YMID=HSIDE*FLOAT(2*IROW-1)/2.0
        NODE(1)=IROWM1*NSIDE+ICOL
        NODE(2)=NODE(1)+1
        NODE(3)=NODE(1)+NSIDE
        NODE(4)=NODE(3)+1
c
c  get gauss points for this element
c
        call p11_gaus (XMID,YMID)
c
c  for each gauss point in this element, evaluate the
c  integrand
c
        DO 60 JGAUSS=1,NGAUSS

          XVAL=XGAUSS(JGAUSS)
          YVAL=YGAUSS(JGAUSS)
c
c  evaluate shape functions
c
          call p11_sh (XMID,XVAL,YMID,YVAL)
c
c  evaluate u and derivatives
c
          call p11_uval (NODE,NSHAPE,NVAR,X)
c
c  evaluate phi
c
          call p11_phi (PHI1,PHI1P,PHI2,PHI2P,XVAL,YVAL)
c
c  evaluate f(lambda,u)
c
          call p11_flam (FLAM,FLAMPL,FLAMPU,UVAL,X(NVAR),XVAL,YVAL)
c
c  compute inner product of equation with each shape function
c  and add to appropriate function
c
          DO 50 KSHAPE=1,NSHAPE
            NOD=NODE(KSHAPE)
            JROWM1=(NOD-1)/NSIDE
            JCOL=NOD-JROWM1*NSIDE
            JROW=JROWM1+1
            IF(JROW.EQ.1.OR.JROW.EQ.NSIDE)GO TO 40
            IF(JCOL.EQ.1.OR.JCOL.EQ.NSIDE)GO TO 40
c
c  interior node
c
            DO LSHAPE=1,NSHAPE
              NOD2=NODE(LSHAPE)
              TERM1=PHI1*DPSIDX(LSHAPE)
     &             +2.0*PHI1P*DUDX*DUDX*DPSIDX(LSHAPE)
     &             +2.0*PHI1P*DUDX*DUDY*DPSIDY(LSHAPE)
              TERM2=PHI2*DPSIDY(LSHAPE)
     &             +2.0*PHI2P*DUDY*DUDX*DPSIDX(LSHAPE)
     &             +2.0*PHI2P*DUDY*DUDY*DPSIDY(LSHAPE)
              TERM3=FLAMPU*PSI(LSHAPE)
              TERM=TERM1*DPSIDX(KSHAPE)+TERM2*DPSIDY(KSHAPE)
     &            +TERM3*PSI(KSHAPE)
              TERM=WGAUSS(JGAUSS)*HSIDE*HSIDE*TERM
              jac(NOD,NOD2)=jac(NOD,NOD2)+TERM
            end do

            TERM=FLAMPL*PSI(KSHAPE)
            TERM=WGAUSS(JGAUSS)*HSIDE*HSIDE*TERM
            jac(NOD,NVAR)=jac(NOD,NVAR)+TERM
            GO TO 50
c
c  boundary node
c
   40       CONTINUE
            jac(NOD,NOD)=1.0
   50       CONTINUE
   60     CONTINUE
   70   CONTINUE
      return
      END
      subroutine p11_nvar ( option, nvar )

c*********************************************************************72
c
cc P11_NVAR sets the number of variables for the finite element rod torsion function.
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
c    Input, integer OPTION, the option chosen for this problem.
c    For some problems, several options are available.  At least,
c    OPTION = 1 is always legal.
c
c    Output, integer NVAR, the number of variables.
c
      implicit none

      integer nvar
      integer option

      nvar = 26

      return
      end
      subroutine p11_option_num ( option_num )

c*********************************************************************72
c
cc P11_OPTION_NUM returns the number of options for the finite element rod torsion function.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    16 September 2008
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, integer OPTION_NUM, the number of options.
c
      implicit none

      integer option_num

      option_num = 2

      return
      end
      subroutine p11_start ( option, nvar, x )

c*********************************************************************72
c
cc P11_START returns a starting point for the finite element rod torsion function.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    19 September 2008
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer OPTION, the option index.
c
c    Input, integer NVAR, the number of variables.
c
c    Output, double precision X(NVAR), the starting point.
c
      implicit none

      integer nvar

      integer i
      integer option
      double precision x(nvar)

      if ( 1 .le. option .and. option .le. 2 ) then

        do i = 1, nvar
          x(i) = 0.0
        end do

      else

        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'P11_START - Fatal error!'
        write ( *, '(a,i8)' ) '  Unrecognized option  = ', option
        stop

      end if

      return
      end
      subroutine p11_stepsize ( option, h, hmin, hmax )

c*********************************************************************72
c
cc P11_STEPSIZE returns step sizes for the finite element rod torsion function.
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
c    Input, integer OPTION, the option index.
c
c    Output, double precision H, HMIN, HMAX, suggested values for the
c    initial step, the minimum step, and the maximum step.
c
      implicit none

      double precision h
      double precision hmax
      double precision hmin
      integer option

      h =    0.12500D+00
      hmin = 0.03125D+00
      hmax = 4.00000D+00

      return
      end
      subroutine p11_title ( option, title )

c*********************************************************************72
c
cc P11_TITLE sets the title for the finite element rod torsion function.
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
c    Input, integer OPTION, the option index.
c
c    Output, character*(*) TITLE, the title of the problem.
c    TITLE will never be longer than 80 characters.
c
      implicit none

      integer option
      character*(*) title

      if ( option .eq. 1 ) then
        title = 'Finite element rod torsion function'
      else if ( option .eq. 2 ) then
        title = 'Finite element rod torsion function.'
      else
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'P11_TITLE - Fatal error!'
        write ( *, '(a)' ) '  Unrecognized option number.'
        stop
      end if

      return
      end
      subroutine p12_data (DIRIPC,HFACT,HMAX,HMIN,HTAN,option,
     & IPC,IPRAM,ISIZE,IT,LIM,MAXCON,
     & MAXLIM,MAXSTP,MAXTAR,NCOL,NROW,NVAR,RELERR,XIT,RWORK)

c*********************************************************************72
c
cc P12_DATA sets parameters for the materially nonlinear function.
c
c  The materially nonlinear problem.
c
c  Reference:
c
c  ivo babuska, werner c rheinboldt,
c  reliable error estimations and mesh adaptation for
c    the finite element method,
c  in  compu methods in nonlinear mechanics,
c  edited  j t oden,
c  north holland publishing company,
c  amsterdam 1980, pages 67-108.
c
c  The function:
c
c  the general problem is the two point boundary value problem
c
c  -d/dt phi(t,u,du/dt,x(nvar))) + psi(t,u,du/dt,x(nvar)) = 0,
c  u(0)=u0(nvar), u(1)=u1(nvar).
c
c  the particular problem solved here has the form
c
c  u'' + x(nvar)*sin(u+u**2+u**3)=0
c
c  u(0)=u(1)=0.0
c
c  with u approximated by piecewise polynomials whose
c  coefficients are the unknowns x(1),...,x(neqn).
c
c  Options:
c
c  option=1  piecewise linear, 1 continuity condition
c  option=2  piecewise cubic, 1 continuity condition
c  option=3  piecewise cubic, 2 continuity conditions
c  option=4  piecewise quintic, 1 continuity condition
c  option=5  piecewise quintic, 2 continuity conditions
c  option=6  piecewise quintic, 3 continuity conditions
c
c  all options use 8 intervals.
c
c  Comments:
c
c  current program has zero as solution for all x(nvar).
c  must find bifurcation branch and jump on to it.
c  perhaps add x(nvar+1) a perturbation to right hand side.
c
      DIMENSION RWORK(ISIZE)
      COMMON /AUXMEM/ BCZERO(8),DBCZDT(8),BCONE(8),DBCODT(8),
     & NPOLYS,NDRV,NVARY,NVARZ,NBCZ,NBCO,
     & PL(8),PLD(8),GCOEF(8),GPOINT(8),NINT,THETA(10,10)
      DATA GCOEF /
     &  .1012285363,  .2223810345,  .3137066459,  .3626837834,
     &  .3626837834,  .3137066459,  .2223810345,  .1012285363/
      DATA GPOINT /
     &  .9602898565,  .7966664774,  .5255324099,  .1834346425,
     & -.1834346425, -.5255324099, -.7966664774, -.9602898565/
      DATA THETA/
     &        1.0,       0.0,       0.0,       0.0,       0.0,
     &        0.0,       0.0,       0.0,       0.0,       0.0,
     &        1.0,       1.0,       0.0,       0.0,       0.0,
     &        0.0,       0.0,       0.0,       0.0,       0.0,
     &        1.0,       3.0,       3.0,       0.0,       0.0,
     &        0.0,       0.0,       0.0,       0.0,       0.0,
     &        1.0,       6.0,      15.0,      15.0,       0.0,
     &        0.0,       0.0,       0.0,       0.0,       0.0,
     &        1.0,      10.0,      45.0,     105.0,     105.0,
     &        0.0,       0.0,       0.0,       0.0,       0.0,
     &        1.0,      15.0,     105.0,     420.0,     945.0,
     &      945.0,       0.0,       0.0,       0.0,       0.0,
     &        1.0,      21.0,     210.0,    1260.0,    4725.0,
     &    10395.0,   10395.0,       0.0,       0.0,       0.0,
     &        1.0,      28.0,     378.0,    3150.0,   17325.0,
     &    62370.0,  135135.0,  135135.0,       0.0,       0.0,
     &        1.0,      36.0,     630.0,    6930.0,   51975.0,
     &   270270.0,  945945.0, 2027025.0, 2027025.0,       0.0,
     &        1.0,      45.0,     990.0,   13860.0,  135135.0,
     &  945945.0,  4729725.0, 16216200.0,34459425.0,34459425.0/

      NINT=8
      IF(option.EQ.1)NPOLYS=2
      IF(option.EQ.2.OR.option.EQ.3)NPOLYS=4
      IF(option.EQ.4.OR.option.EQ.5.OR.option.EQ.6)NPOLYS=6
      NVARY=NINT*NPOLYS
      NBCZ=1
      IF(option.EQ.1.OR.option.EQ.2.OR.option.EQ.4)NDRV=1
      IF(option.EQ.3.OR.option.EQ.5)NDRV=2
      IF(option.EQ.6)NDRV=3
      NBCO=1
      NVARZ=NBCZ+(NINT-1)*NDRV+NBCO
      NVAR=NVARY+NVARZ+1

      NEQN=NVAR-1
      DIRIPC=1.0
      HFACT=3.0
      IPC=NVAR
      IT=0
      LIM=NVAR
      MAXCON=30
      MAXLIM=30
      MAXSTP=30
      MAXTAR=30
      NCOL=NVAR
      NROW=NVAR
      RELERR=1.0E-5
      XIT=0.0

      call p12_bd ( )
      return
      END
      subroutine p12_bd ( )

c*********************************************************************72
c
cc P12_BD sets boundary values for the materially nonlinear problem.
c
      COMMON /AUXMEM/ BCZERO(8),DBCZDT(8),BCONE(8),DBCODT(8),
     & NPOLYS,NDRV,NVARY,NVARZ,NBCZ,NBCO,
     & PL(8),PLD(8),GCOEF(8),GPOINT(8),NINT,THETA(10,10)
      BCONE(1)=0.0
      BCZERO(1)=0.0
      DBCODT(1)=0.0
      DBCZDT(1)=0.0
      return
      END
      subroutine p12_fun ( option, nvar, x, fx )

c*********************************************************************72
c
cc P12_FUN evaluates the materially nonlinear function.
c
      integer nvar

      double precision fx(nvar-1)
      integer option
      double precision x(nvar)

      COMMON /AUXMEM/ BCZERO(8),DBCZDT(8),BCONE(8),DBCODT(8),
     & NPOLYS,NDRV,NVARY,NVARZ,NBCZ,NBCO,
     & PL(8),PLD(8),GCOEF(8),GPOINT(8),NINT,THETA(10,10)
c
c  zero out fx
c
      call p12_bd ( )
      NEQN=NVAR-1
      DO J=1,NEQN
        FX(J)=0.0
      end do
c
c  1. set up the terms (a*y) involving the bivariate form
c
c  i  for each interval
c
      DO I=1,NINT

        ISKIP=(I-1)*NPOLYS
        XL=FLOAT(I-1)/FLOAT(NINT)
        XR=FLOAT(I)/FLOAT(NINT)
        DTDX=2.0/(XR-XL)
        DXDT=(XR-XL)/2.0
c
c  j  for each gauss point, evaluate the integrand
c
        DO J=1,8

          T=GPOINT(J)
          COEF=GCOEF(J)*DXDT
          XG=XL+(T+1.0)*DXDT
          call p12_vl (T,DTDX,X,NVAR,ISKIP,U,UPRYM,NPOLYS,PL,PLD)
          TCON=X(NVAR)
          call p12_gx (XG,U,UPRYM,PHI,PSI,TCON)
          LSKIP=ISKIP
c
c  l   project onto each test function pl(l) and pld(l)
c
          DO L=1,NPOLYS
            TERM=PSI*PL(L)
            TERM=TERM+PHI*PLD(L)
            IEQN=LSKIP+L
            FX(IEQN)=FX(IEQN)+COEF*TERM
          end do

          LSKIP=LSKIP+NPOLYS

        end do

      end do
c
c  2. add the terms involving the continuity of the test functions
c  which are the terms b*z in  f=a*y + b*z
c
      NCL=NVARY
      NCR=NVARY+NBCZ
c
c  i  for each interval
c
      DO 110 I=1,NINT

        XL=FLOAT(I-1)/FLOAT(NINT)
        XR=FLOAT(I)/FLOAT(NINT)
        DTDX=2.0/(XR-XL)
c
c  j   for the polynomials used in approximating each u,
c  count conditions at left endpoint, lhil, and at right, lhir
c  if we are in the first or last interval, one of
c  these will be boundary conditions
c
        LHIL=NDRV
        LHIR=NDRV
        IF (I.EQ.1) LHIL=NBCZ
        IF (I.EQ.NINT) LHIR=NBCO
c
c  k  for each test function pl(k)
c
        DO K=1,NPOLYS

          TERML=0.0
          TERMR=0.0
          S=(-1.0)**(K+1)
          IEQN=(I-1)*NPOLYS+K
c
c  l  for derivatives 0 thru nbcz specified at 0.0 to be zero
c  or derivatives 0 thru ndrv specified to be continuous
c  at interior nodes
c  or derivatives 0 thru nbco specified at 1.0 to be zero
c
          H2I=1.0
          DO L=1,LHIL
            S=-S
            INDEX=NCL+L
            TERML=TERML+S*X(INDEX)*H2I*THETA(L,K)
            H2I=H2I*DTDX
          end do

          H2I=1.0
          DO L=1,LHIR
            INDEX=NCR+L
            TERMR=TERMR+H2I*THETA(L,K)*X(INDEX)
            H2I=H2I*DTDX
          end do

          FX(IEQN)=FX(IEQN)+TERML+TERMR

        end do

        IF (LHIL.GT.0) NCL=NCL+LHIL
        IF (LHIR.GT.0) NCR=NCR+LHIR

      end do
c
c  3. create the terms for the u functions and their derivatives
c  the matrix terms  ( c*y )
c  one equation is generated for component and condition
c
      NPSUM=0
      DTDXR=0.0
      DTDXL=0.0
c
c  i  for each node
c
      NDSUM=NVARY
      NODES=NINT+1

      DO I=1,NODES

        IF (I.GT.1) XL=FLOAT(I-2)/FLOAT(NINT)
        XC=FLOAT(I-1)/FLOAT(NINT)
        IF (I.LT.NODES) XR=FLOAT(I)/FLOAT(NINT)
        IF (XC.NE.XL) DTDXL=2.0/(XC-XL)
        IF (XR.NE.XC) DTDXR=2.0/(XR-XC)
        H2IL=1.0
        H2IR=1.0
c
c  k   for derivatives 0 thru nbcz specified at 0.0
c  or derivatives 0 thru ndrv specified to be continuous
c  at interior nodes
c  or derivatives 0 thru nbco specified at 1.0
c
        KHI=NDRV
        IF (I.EQ.1) KHI=NBCZ
        IF (I.EQ.NODES) KHI=NBCO

        DO K=1,KHI

          S=(-1.0)**(K+1)
c
c  l   set up the term from the left hand interval
c
          IF ( I .eq. 1 ) then

            TERML=BCZERO(K)

          else

            TERML=0.0
            DO L=1,NPOLYS
              IVAR=NPSUM+L-NPOLYS
              TERML=TERML+X(IVAR)*H2IL*THETA(K,L)
            end do

          end if
c
c  l   set up the term from the right hand interval
c
          IF ( I .eq. NODES ) then

            TERMR=-BCONE(K)

          else

            TERMR=0.0
            DO L=1,NPOLYS
              IVAR=NPSUM+L
              S=-S
              TERMR=TERMR+S*X(IVAR)*H2IR*THETA(K,L)
            end do

          end if

          IEQN=NDSUM+K
          FX(IEQN)=TERML+TERMR
          H2IL=H2IL*DTDXL
          H2IR=H2IR*DTDXR

        end do

        NDSUM=NDSUM+KHI
        NPSUM=NPSUM+NPOLYS

      end do

      return
      end
      subroutine p12_gx ( xg, u, uprym, phi, psi, tcon )

c*********************************************************************72
c
cc P12_GX is an auxilliary routine for the materially nonlinear problem.
c
      implicit none

      double precision phi
      double precision psi
      double precision tcon
      double precision u
      double precision uprym
      double precision xg

      phi = - uprym
      psi = tcon * sin ( u * ( 1.0D+00 + u* ( 1.0D+00 + u )))

      return
      end
      subroutine p12_jac ( option, nvar, x, jac )

c*********************************************************************72
c
cc P12_JAC evaluates the materially nonlinear jacobian.
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
c    Input, integer OPTION, the option index.
c
c    Input, integer NVAR, the number of variables.
c
c    Input, double precision X(NVAR), the argument of the jacobian.
c
c    Output, double precision JAC(NVAR,NVAR), the jacobian matrix evaluated
c    at X.  The NVAR-th row is not set by this routine.
c
      implicit none

      integer nvar

      integer option
      double precision jac(nvar,nvar)
      double precision x(nvar)

      COMMON /AUXMEM/ BCZERO(8),DBCZDT(8),BCONE(8),DBCODT(8),
     & NPOLYS,NDRV,NVARY,NVARZ,NBCZ,NBCO,
     & PL(8),PLD(8),GCOEF(8),GPOINT(8),NINT,THETA(10,10)
c
c  zero out the matrix
c
      do j = 1, ncol
        do i = 1, nrow
          jac(i,j) = 0.0E+00
        end do
      end do

      TCON=X(NVAR)
      call p12_bd ( )
c
c  1.  set up the terms from the bivariate form (a*y )
c
c      i  for each interval (xl,xr)
c
      DO I=1,NINT
        XL=FLOAT(I-1)/FLOAT(NINT)
        XR=FLOAT(I)/FLOAT(NINT)
        ISKIP=(I-1)*NPOLYS
        DTDX=2.0/(XR-XL)
        DXDT=(XR-XL)/2.0
c
c  j  for each gauss point xg in the interval
c  note that the legendre polynomials and derivatives are set up
c  by the call to value
c
        DO J=1,8
          T=GPOINT(J)
          COEF=GCOEF(J)*DXDT
          XG=XL+DXDT*(T+1.0)
          call p12_vl (T,DTDX,X,NVAR,ISKIP,U,UPRYM,NPOLYS,PL,PLD)
            call p12_gp (XG,U,UPRYM,PHIPU,PHIPUP,PHIPT,PSIPU,
     &      PSIPUP,PSIPT,TCON)
c
c  l  for each legendre polynomial coefficient
c
            DO L=1,NPOLYS

              IEQN=ISKIP+L
              TERM=COEF*(PSIPT*PL(L)+PHIPT*PLD(L))
              jac(IEQN,NVAR)=jac(IEQN,NVAR)+TERM
c
c  n   for each y-coefficient of a u
c
              DO N=1,NPOLYS
                IVAR=NPOLYS*(I-1)+N
                TERM1=PSIPU*PL(N)*PL(L)
                TERM2=PSIPUP*PLD(N)*PL(L)
                TERM3=PHIPU*PL(N)*PLD(L)
                TERM4=PHIPUP*PLD(N)*PLD(L)
                SUM=TERM1+TERM2+TERM3+TERM4
                jac(IEQN,IVAR)=jac(IEQN,IVAR)+COEF*SUM
              end do

            end do

        end do
      end do
c
c  2. add the terms involving the continuity of the test functions
c  which are the terms b*z in  f=a*y+b*z
c
c
c     i  for each interval,
c
      DO 150 I=1,NINT

        NCL=NVARY
        IF (I.GT.1) NCL=NVARY+NBCZ+(I-2)*NDRV
        NCR=NVARY+NBCZ+(I-1)*NDRV
        XL=FLOAT(I-1)/FLOAT(NINT)
        XR=FLOAT(I)/FLOAT(NINT)
        DTDX=2.0/(XR-XL)
c
c  j   for the polynomials used in approximating each u,
c      count conditions at left endpoint, lhil, and at right, lhir
c      if we are in the first or last interval, one of
c      these will be boundary conditions
c
          LHIL=NDRV
          LHIR=NDRV
          IF (I.EQ.1) LHIL=NBCZ
          IF (I.EQ.NINT) LHIR=NBCO
          IF (LHIL.LE.0.AND.LHIR.LE.0) GO TO 130
c
c  k  for each test function pl(k)
c
          DO 120 K=1,NPOLYS
            S=(-1.0)**(K+1)
            IEQN=(I-1)*NPOLYS+K
c
c  l  for derivatives 0 thru nbcz specified at 0.0 to be zero
c     or derivatives 0 thru ndrv specified to be continuous
c     at interior nodes
c     or derivatives 0 thru nbco specified at 1.0 to be zero
c
c     evaluate contribution from left endpoint
c
            IF (LHIL.LE.0) GO TO 100
            H2I=1.0
            DO L=1,LHIL
              S=-S
              IVAR=NCL+L
              jac(IEQN,IVAR)=S*H2I*THETA(L,K)
              H2I=H2I*DTDX
            end do
c
c  evaluate contribution from right endpoint
c
  100       IF (LHIR.LE.0) GO TO 120

            H2I=1.0
            DO L=1,LHIR
              IVAR=NCR+L
              jac(IEQN,IVAR)=H2I*THETA(L,K)
              H2I=H2I*DTDX
            end do

  120       CONTINUE
          IF (LHIL.GT.0) NCL=NCL+LHIL
          IF (LHIR.GT.0) NCR=NCR+LHIR
  130     CONTINUE
  140     CONTINUE
  150   CONTINUE
c
c  3. create the terms for the u functions and their derivatives
c  the matrix terms  ( c*y )
c  one equation is generated for component and condition
c
c
c  i  for each interval
c
      DO 230 I=1,NINT

        NCL=NVARY
        IF (I.GT.1) NCL=NVARY+NBCZ+(I-2)*NDRV
        NCR=NVARY+NBCZ+(I-1)*NDRV
        NPSUM=(I-1)*NPOLYS
        XL=FLOAT(I-1)/FLOAT(NINT)
        XR=FLOAT(I)/FLOAT(NINT)
        DTDX=2.0/(XR-XL)
          H2I=1.0
c
c  k   for derivatives 0 thru nbcz specified at 0.0
c      or derivatives 0 thru ndrv specified to be continuous
c      at interior nodes
c      or derivatives 0 thru nbco specified at 1.0
c
          KHIL=NDRV
          IF (I.EQ.1) KHIL=NBCZ
          IF (KHIL.LE.0) GO TO 180
          DO 170 K=1,KHIL
           IEQN=NCL+K
c
c  l   set up the term from the left hand endpoint
c
            IF (I.EQ.1) jac(IEQN,NVAR)=DBCZDT(K)
            S=(-1.0)**(K+1)
            DO L=1,NPOLYS
              IVAR=NPSUM+L
              S=-S
              jac(IEQN,IVAR)=S*H2I*THETA(K,L)
            end do

            H2I=H2I*DTDX
  170       CONTINUE
          NCL=NCL+KHIL
  180     CONTINUE
          H2I=1.0
          KHIR=NDRV
          IF (I.EQ.NINT) KHIR=NBCO
          IF (KHIR.LE.0) GO TO 210
          DO K=1,KHIR
            IEQN=NCR+K
            jac(IEQN,NVAR)=0.0
            IF (I.EQ.NINT) jac(IEQN,NVAR)=-DBCODT(K)
            DO L=1,NPOLYS
              IVAR=NPSUM+L
              jac(IEQN,IVAR)=H2I*THETA(K,L)
            end do
            H2I=H2I*DTDX
          end do
          NCR=NCR+KHIR
  210     NPSUM=NPSUM+NPOLYS
  220     CONTINUE
  230   CONTINUE
      return
      END
      subroutine p12_vl (T,DTDX,X,NVAR,ISKIP,U,UPRYM,NPOLYS,PL,PLD)

c*********************************************************************72
c
cc P12_VL is an auxilliary routine for the materially nonlinear problem.
c
      dimension x(nvar),pl(8),pld(8)
c
c  evaluate legendre polynomials and derivatives
c
      pl(1)=1.0
      pld(1)=0.0
      pl(2)=t
      pld(2)=1.0
      a=0.0
      do i=3,8
        im1=i-1
        im2=i-2
        a=a+1.0
        pl(i)=((a+a+1.0)*t*pl(im1)-a*pl(im2))/(a+1.0)
        pld(i)=((a+a+1.0)*(t*pld(im1)+pl(im1))-a*pld(im2))/(a+1.0)
      end do

      do i=1,8
        pld(i)=dtdx*pld(i)
      end do
c
c  sum on coefficents to evaluate u and uprym
c
      u=0.0
      uprym=0.0
      do j=1,npolys
        index=iskip+j
        u=u+x(index)*pl(j)
        uprym=uprym+x(index)*pld(j)
      end do

      return
      end
      subroutine p12_gp (xg,u,uprym,phipu,phipup,phipt,psipu,
     & psipup,psipt,tcon)

c*********************************************************************72
c
cc P12_GP is an auxilliary function for the materially nonlinear problem.
c
      phipu=0.0
      phipup=-1.0
      phipt=0.0
      psipu=(1.0+u*(2.0+3.0*u))*cos(u*(1.0+u*(1.0+u)))
      psipup=0.0
      psipt=sin(u*(1.0+u*(1.0+u)))
      return
      end
      subroutine p12_nvar ( option, nvar )

c*********************************************************************72
c
cc P12_NVAR sets the number of variables for the materially nonlinear function.
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
c    Input, integer OPTION, the option chosen for this problem.
c    For some problems, several options are available.  At least,
c    OPTION = 1 is always legal.
c
c    Output, integer NVAR, the number of variables.
c
      implicit none

      integer nbco
      integer nbcz
      integer ndrv
      integer nint
      integer npolys
      integer nvar
      integer nvary
      integer nvarz
      integer option

      nint = 8

      if ( option .eq. 1 ) then
        npolys = 2
      else if ( option .eq. 2 .or.
     &          option .eq. 3 ) then
        npolys = 4
      else if ( option .eq. 4 .or.
     &          option .eq. 5 .or.
     &          option .eq. 6 ) then
        npolys = 6
      end if

      nvary = nint * npolys
      nbcz = 1

      if ( option .eq. 1 .or.
     &     option .eq. 2 .or.
     &     option .eq. 4 ) then
        ndrv = 1
      else if ( option .eq. 3 .or.
     &     option .eq. 5 ) then
        ndrv = 2
      else if ( option .eq. 6 ) then
        ndrv = 3
      end if

      nbco = 1
      nvarz = nbcz + ( nint - 1 ) * ndrv + nbco
      nvar = nvary + nvarz + 1

      return
      end
      subroutine p12_option_num ( option_num )

c*********************************************************************72
c
cc P12_OPTION_NUM returns the number of options for the materially nonlinear function.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    16 September 2008
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, integer OPTION_NUM, the number of options.
c
      implicit none

      integer option_num

      option_num = 6

      return
      end
      subroutine p12_start ( option, nvar, x )

c*********************************************************************72
c
cc P12_START returns a starting point for the materially nonlinear function.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    19 September 2008
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer OPTION, the option index.
c
c    Input, integer NVAR, the number of variables.
c
c    Output, double precision X(NVAR), the starting point.
c
      implicit none

      integer nvar

      integer i
      integer option
      double precision x(nvar)

      if ( 1 .le. option .and. option .le. 6 ) then

        do i = 1, nvar
          x(i) = 0.0
        end do

      else

        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'P12_START - Fatal error!'
        write ( *, '(a,i8)' ) '  Unrecognized option  = ', option
        stop

      end if

      return
      end
      subroutine p12_stepsize ( option, h, hmin, hmax )

c*********************************************************************72
c
cc P12_STEPSIZE returns step sizes for the materially nonlinear function.
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
c    Input, integer OPTION, the option index.
c
c    Output, double precision H, HMIN, HMAX, suggested values for the
c    initial step, the minimum step, and the maximum step.
c
      implicit none

      double precision h
      double precision hmax
      double precision hmin
      integer option

      h =     2.000D+00
      hmin =  0.001D+00
      hmax = 10.000D+00

      return
      end
      subroutine p12_title ( option, title )

c*********************************************************************72
c
cc P12_TITLE sets the title for the materially nonlinear function.
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
c    Input, integer OPTION, the option index.
c
c    Output, character*(*) TITLE, the title of the problem.
c    TITLE will never be longer than 80 characters.
c
      implicit none

      integer option
      character*(*) title

      if ( option .eq. 1 ) then
        title = 'Materially nonlinear problem, NPOLYS = 2, NDERIV = 1.'
      else if ( option .eq. 2 ) then
        title = 'Materially nonlinear problem, NPOLYS = 4, NDERIV = 1.'
      else if ( option .eq. 3 ) then
        title = 'Materially nonlinear problem, NPOLYS = 4, NDERIV = 2.'
      else if ( option .eq. 4 ) then
        title = 'Materially nonlinear problem, NPOLYS = 6, NDERIV = 1.'
      else if ( option .eq. 5 ) then
        title = 'Materially nonlinear problem, NPOLYS = 6, NDERIV = 2.'
      else if ( option .eq. 6 ) then
        title = 'Materially nonlinear problem, NPOLYS = 6, NDERIV = 3.'
      else
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'P12_TITLE - Fatal error!'
        write ( *, '(a)' ) '  Unrecognized option number.'
        stop
      end if

      return
      end
      subroutine p13_data (DIRIPC,HFACT,HMAX,HMIN,HTAN,option,
     & IPC,IPRAM,ISIZE,IT,LIM,MAXCON,
     & MAXLIM,MAXSTP,MAXTAR,NCOL,NROW,NVAR,RELERR,XIT,RWORK)

c*********************************************************************72
c
cc P13_DATA sets parameters for Simpson's boundary value function.
c
c  13.1 simpsons mildly nonlinear boundary value problem
c
c  13.2 reference
c
c  r b simpson,
c  a method for the numerical determination of bifurcation
c    states of nonlinear systems of equations,
c  siam journal of numerical analysis, 12, 1975, pages 439-451.
c
c  13.3  the function
c
c  the problem has the general form
c
c  -del(u(x,y)) = x(nvar) * f(u(x,y)).
c  u(x,y)=0 on boundary of unit square.
c
c  a uniform mesh is used.  if del9=nine point discrete laplacian,
c  and del5=five point discrete laplacian, the system has the
c  form
c
c  del9 x + x(nvar) * (f(x) + h*h*del5(f(x))/12) = 0.0
c  where x is the mesh function approximation to u.
c
c  13.4 options
c
c  option=1  8 intervals, right hand side is exp(x)
c  option=2  8 intervals, right hand side is 1 +(x+.5*x*x)/(1+.01*x*x)
c
c  option=3  12 intervals, right hand side is exp(x)
c  option=4  12 intervals, right hand side is 1 +(x+.5*x*x)/(1+.01*x*x)
c
c  note option=3 or option=4 require isize=15,494
c
c  option=5  16 intervals, right hand side is exp(x)
c  option=6  16 intervals, right hand side is 1 +(x+.5*x*x)/(1+.01*x*x)
c
c  note option=5 or option=6 require isize very large.
c
c  13.6 limit points
c
c  melhem lists
c
c  with x(25)=the value of the function at the center point,
c
c  for option=1, xnvar=6.807504, x(25)=1.391598
c  for option=2, xnvar=7.980356, x(25)=2.272364
c  for option=3, xnvar=6.808004, x(25)=1.391657
c  for option=4, xnvar=7.981426, x(25)=2.273045
c  for option=5, xnvar=6.808087, x(25)=1.391656
c  for option=6, xnvar=7.981605, x(25)=2.273159
c
c  computational results
c
c
c  13.8 comments
c
c  on 14 october 1983, the jacobian routine was debugged.  all
c  results before that time are invalid.  a check can be made
c  by comparing runs with ijac=0 and ijac=1.
c
c  options 3 through 6 cannot be run until a band storage mode is used
c  with special treatment of column nvar.
c  or else more storage allocated through isize.
c
      DIMENSION RWORK(ISIZE)
      COMMON /AUXMEM/ IFUN13,M,H

      nvar=m*m+1
      H=1.0/FLOAT(M+1)
      DIRIPC=1.0
      HFACT=3.0
      IPC=NVAR
      IT=NVAR
      LIM=NVAR
      MAXCON=75
      MAXLIM=1
      MAXSTP=75
      MAXTAR=75
      NCOL=NVAR
      NROW=NVAR
      RELERR=1.0E-5
      XIT=1.0
      if ( option .eq. 1 .or.
     &     option .eq. 3 .or.
     &     option .eq. 5 ) then
        ifun13 = 1
      else if ( option .eq. 2 .or.
     &          option .eq. 4 .or.
     &          option .eq. 6 ) then
        ifun13 = 2
      end if

      return
      END
      subroutine p13_fun ( option, nvar, x, fx )

c*********************************************************************72
c
cc P13_FUN evaluates Simpson's mildly nonlinear boundary value function.
c
      integer nvar

      double precision fx(nvar-1)
      integer option
      double precision x(nvar)

      if ( option .eq. 1 .or. option .eq. 2 ) then
        m = 7
      else if ( option .eq. 3 .or. option .eq. 4 ) then
        m = 11
      else if ( option .eq. 5 .or. option .eq. 6 ) then
        m = 15
      end if

      H=1.0/FLOAT(M+1)

      NEQN=NVAR-1
      FAC1=1.0/(6.0*H*H)
      FAC2=X(NVAR)/12.0

      DO I=1,M
        IPM=I+M
        TERM1=FAC1*(4.0*X(IPM)-20.0*X(I))
        TERM2=X(NVAR)*
     &      (p13_gx(option,X(I)) + (p13_gx(option,X(IPM))
     &    - 4.0*p13_gx(option,X(I)))/12.0)
        FX(I)=TERM1+TERM2
      end do

      I=1
      IP1=I+1
      IPMP1=I+M+1
      FX(I)=FX(I)+FAC1*(4.0*X(IP1)+X(IPMP1))
     &  +FAC2*(p13_gx(option,X(IP1))+2.0)
      IHI=M-1
      DO I=2,IHI
        IM1=I-1
        IP1=IP1
        IPMM1=I+M-1
        IPMP1=I+M+1
        TERM1=FAC1*(4.0*X(IM1)+4.0*X(IP1)+X(IPMM1)+X(IPMP1))
        TERM2=FAC2*(p13_gx(option,X(IM1))+p13_gx(option,X(IP1))+1.0)
        FX(I)=FX(I)+TERM1+TERM2
      end do

      I=M
      IM1=I-1
      IPMM1=I+M-1
      FX(I)=FX(I)+FAC1*(4.0*X(IM1)+X(IPMM1))
     &  +FAC2*(p13_gx(option,X(IM1))+2.0)
      JHI=M-1

      DO J=2,JHI

        ILO=M*(J-1)+1
        IHI=M*(J-1)+M
        DO I=ILO,IHI
          IPM=I+M
          IMM=I-M
          TERM1=FAC1*(4.0*X(IMM)-20.0*X(I)+4.0*X(IPM))
          TERM2=(p13_gx(X(IMM))+8.*p13_gx(option,X(I))
     &      +p13_gx(option,X(IPM)))/12.0
          FX(I)=TERM1+X(NVAR)*TERM2
        end do

        I=M*(J-1)+1
        IP1=I+1
        IPMP1=I+M+1
        IMMP1=I-M+1
        TERM1=FAC1*(X(IMMP1)+4.0*X(IP1)+X(IPMP1))
        TERM2=FAC2*(p13_gx(option,X(IP1))+1.0)
        FX(I)=FX(I)+TERM1+TERM2
        ILO=M*(J-1)+2
        IHI=M*(J-1)+M-1

        DO I=ILO,IHI
          IP1=I+1
          IM1=I-1
          IPMM1=I+M-1
          IPMP1=I+M+1
          IMMM1=I-M-1
          IMMP1=I-M+1
          TERM1=X(IMMM1)+X(IMMP1)+4.0*X(IM1)
          TERM2=4.0*X(IP1)+X(IPMM1)+X(IPMP1)
          TERM3=p13_gx(option,X(IM1))+p13_gx(option,X(IP1))
          FX(I)=FX(I)+FAC1*(TERM1+TERM2)+FAC2*TERM3
        end do

        I=M*(J-1)+M
        IM1=I-1
        IPMM1=I+M-1
        IMMM1=I-M-1
        TERM1=FAC1*(X(IMMM1)+4.0*X(IM1)+X(IPMM1))
        TERM2=FAC2*(p13_gx(option,X(IM1))+1.0)
        FX(I)=FX(I)+TERM1+TERM2

      end do

      ILO=NEQN-(M-1)
      IHI=NEQN

      DO I=ILO,IHI
        IMM=I-M
        TERM1=FAC1*(4.0*X(IMM)-20.0*X(I))
        TERM2=X(NVAR)*
     &     (p13_gx(option,X(I))+(p13_gx(option,X(IMM))
     &    -4.0*p13_gx(option,X(I)))/12.0)
        FX(I)=TERM1+TERM2
      end do

      I=NEQN-(M-1)
      IP1=I+1
      IMMP1=I-M+1
      FX(I)=FX(I)+FAC1*(X(IMMP1)+4.0*X(IP1))
     &  +FAC2*(p13_gx(option,X(IP1))+2.0)
      ILO=NEQN-(M-2)
      IHI=NEQN-1

      DO I=ILO,IHI
        IP1=I+1
        IM1=I-1
        IMMM1=I-M-1
        IMMP1=I-M+1
        TERM1=FAC1*(X(IMMM1)+X(IMMP1)+4.0*X(IM1)+4.0*X(IP1))
        TERM2=FAC2*(p13_gx(option,X(IM1))+p13_gx(option,X(IP1))+1.0)
        FX(I)=FX(I)+TERM1+TERM2
      end do

      I=NEQN
      IM1=I-1
      IMMM1=I-M-1
      FX(I)=FX(I)+FAC1*(X(IMMM1)+4.0*X(IM1))
     &  +FAC2*(p13_gx(option,X(IM1))+2.0)
      return
      END
      FUNCTION p13_gp ( option, x )

c*********************************************************************72
c
cc P13_GP is an auxilliary function for simpsons mildly nonlinear boundary function.
c
      integer option

      if ( option .eq. 1 .or.
     &     option .eq. 3 .or.
     &     option .eq. 5 ) then
        ifun13 = 1
      else if ( option .eq. 2 .or.
     &          option .eq. 4 .or.
     &          option .eq. 6 ) then
        ifun13 = 2
      end if

      IF(IFUN13.EQ.1) then
        p13_gp=EXP(X)
      else IF ( IFUN13.EQ.2 ) then
        p13_gp=(1.0+X-.01*X*X)/((1.0+.01*X*X)**2)
      end if

      return
      END
      FUNCTION p13_gx ( option, x )

c*********************************************************************72
c
cc P13_GX is an auxilliary function for simpsons mildly nonlinear boundary function.
c
      integer option

      if ( option .eq. 1 .or.
     &     option .eq. 3 .or.
     &     option .eq. 5 ) then
        ifun13 = 1
      else if ( option .eq. 2 .or.
     &          option .eq. 4 .or.
     &          option .eq. 6 ) then
        ifun13 = 2
      end if


      IF(IFUN13.EQ.1)then
        p13_gx=EXP(X)
      else IF(IFUN13.EQ.2) then
        p13_gx=(100.0+100.0*X+51.0*X*X)/(100.0+X*X)
      end if

      return
      END
      subroutine p13_jac ( option, nvar, x, jac )

c*********************************************************************72
c
cc P13_JAC evaluates the Simpson's boundary value jacobian.
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
c    Input, integer OPTION, the option index.
c
c    Input, integer NVAR, the number of variables.
c
c    Input, double precision X(NVAR), the argument of the jacobian.
c
c    Output, double precision JAC(NVAR,NVAR), the jacobian matrix evaluated
c    at X.  The NVAR-th row is not set by this routine.
c
      implicit none

      integer nvar

      integer option
      double precision jac(nvar,nvar)
      double precision x(nvar)

      COMMON /AUXMEM/ IFUN13,M,H

      if ( option .eq. 1 .or. option .eq. 2 ) then
        m = 7
      else if ( option .eq. 3 .or. option .eq. 4 ) then
        m = 11
      else if ( option .eq. 5 .or. option .eq. 6 ) then
        m = 15
      end if

      H=1.0/FLOAT(M+1)

      do j = 1, ncol
        do i = 1, nrow
          jac(i,j) = 0.0E+00
        end do
      end do

      NEQN=NVAR-1

      DO I=1,M
        IPM=I+M
        jac(I,IPM)=jac(I,IPM)+4.0/(6.0*H*H)
        jac(I,I)=jac(I,I)-20.0/(6.0*H*H)
        jac(I,NVAR)=jac(I,NVAR)+
     &      (p13_gx(option,X(I))+(p13_gx(option,X(IPM))
     &    -4.0*p13_gx(option,X(I)))/12.0)
        jac(I,I)=jac(I,I)+X(NVAR)*p13_gp(option,X(I))
        jac(I,IPM)=jac(I,IPM)+X(NVAR)*p13_gp(option,X(IPM))/12.0
        jac(I,I)=jac(I,I)-X(NVAR)*p13_gp(option,X(I))/3.0
      end do

      I=1
      IP1=I+1
      IPMP1=I+M+1
      jac(I,IP1)=jac(I,IP1)+4.0/(6.0*H*H)
      jac(I,IPMP1)=jac(I,IPMP1)+1.0/(6.0*H*H)
      jac(I,NVAR)=jac(I,NVAR)+(p13_gx(option,X(IP1))+2.0)/12.0
      jac(I,IP1)=jac(I,IP1)+X(NVAR)*p13_gp(option,X(IP1))/12.0
      IHI=M-1

      DO I=2,IHI
        IM1=I-1
        IP1=IP1
        IPMM1=I+M-1
        IPMP1=I+M+1
        jac(I,IM1)=jac(I,IM1)+4.0/(6.0*H*H)
        jac(I,IP1)=jac(I,IP1)+4.0/(6.0*H*H)
        jac(I,IPMM1)=jac(I,IPMM1)+1.0/(6.0*H*H)
        jac(I,IPMP1)=jac(I,IPMP1)+1.0/(6.0*H*H)
        jac(I,NVAR)=jac(I,NVAR)+
     &  (p13_gx(option,X(IM1))+p13_gx(option,X(IP1))+1.0)/12.0
        jac(I,IM1)=jac(I,IM1)+X(NVAR)*p13_gp(option,X(IM1))/12.0
        jac(I,IP1)=jac(I,IP1)+X(NVAR)*p13_gp(option,X(IP1))/12.0
      end do

      I=M
      IM1=I-1
      IPMM1=I+M-1
      jac(I,IM1)=jac(I,IM1)+4.0/(6.0*H*H)
      jac(I,IPMM1)=jac(I,IPMM1)+1.0/(6.0*H*H)
      jac(I,NVAR)=jac(I,NVAR)+(p13_gx(option,X(IM1))+2.0)/12.0
      jac(I,IM1)=jac(I,IM1)+X(NVAR)*p13_gp(option,X(IM1))/12.0
      JHI=M-1

      DO J=2,JHI

        ILO=M*(J-1)+1
        IHI=M*(J-1)+M

        DO I=ILO,IHI
          IPM=I+M
          IMM=I-M
          jac(I,IMM)=jac(I,IMM)+4.0/(6.0*H*H)
          jac(I,I)=jac(I,I)-20.0/(6.0*H*H)
          jac(I,IPM)=jac(I,IPM)+4.0/(6.0*H*H)
          jac(I,NVAR)=jac(I,NVAR)+
     &    (p13_gx(option,X(IMM))+8.0*p13_gx(option,X(I))+p13_gx(option,X(IPM)))/12.0
          jac(I,IMM)=jac(I,IMM)+X(NVAR)*p13_gp(option,X(IMM))/12.0
          jac(I,I)=jac(I,I)+X(NVAR)*8.0*p13_gp(option,X(I))/12.0
          jac(I,IPM)=jac(I,IPM)+X(NVAR)*p13_gp(option,X(IPM))/12.0
        end do

        I=M*(J-1)+1
        IP1=I+1
        IPMP1=I+M+1
        IMMP1=I-M+1
        jac(I,IMMP1)=jac(I,IMMP1)+1.0/(6.0*H*H)
        jac(I,IP1)=jac(I,IP1)+4.0/(6.0*H*H)
        jac(I,IPMP1)=jac(I,IPMP1)+1.0/(6.0*H*H)
        jac(I,NVAR)=jac(I,NVAR)+(p13_gx(option,X(IP1))+1.0)/12.0
        jac(I,IP1)=jac(I,IP1)+X(NVAR)*p13_gp(option,X(IP1))/12.0
        ILO=M*(J-1)+2
        IHI=M*(J-1)+M-1

        DO I=ILO,IHI
          IP1=I+1
          IM1=I-1
          IPMM1=I+M-1
          IPMP1=I+M+1
          IMMM1=I-M-1
          IMMP1=I-M+1
          jac(I,IMMM1)=jac(I,IMMM1)+1.0/(6.0*H*H)
          jac(I,IMMP1)=jac(I,IMMP1)+1.0/(6.0*H*H)
          jac(I,IM1)=jac(I,IM1)+4.0/(6.0*H*H)
          jac(I,IP1)=jac(I,IP1)+4.0/(6.0*H*H)
          jac(I,IPMM1)=jac(I,IPMM1)+1.0/(6.0*H*H)
          jac(I,IPMP1)=jac(I,IPMP1)+1.0/(6.0*H*H)
          jac(I,NVAR)=jac(I,NVAR)+
     &    (p13_gx(option,X(IM1))+p13_gx(option,X(IP1)))/12.0
          jac(I,IM1)=jac(I,IM1)+X(NVAR)*p13_gp(option,X(IM1))/12.0
          jac(I,IP1)=jac(I,IP1)+X(NVAR)*p13_gp(option,X(IP1))/12.0
        end do

        I=M*(J-1)+M
        IM1=I-1
        IPMM1=I+M-1
        IMMM1=I-M-1
        jac(I,IMMM1)=jac(I,IMMM1)+1.0/(6.0*H*H)
        jac(I,IM1)=jac(I,IM1)+4.0/(6.0*H*H)
        jac(I,IPMM1)=jac(I,IPMM1)+1.0/(6.0*H*H)
        jac(I,NVAR)=jac(I,NVAR)+(p13_gx(option,X(IM1))+1.0)/12.0
        jac(I,IM1)=jac(I,IM1)+X(NVAR)*p13_gp(option,X(IM1))/12.0

      end do

      ILO=NEQN-(M-1)
      IHI=NEQN

      DO I=ILO,IHI
        IMM=I-M
        jac(I,IMM)=jac(I,IMM)+4.0/(6.0*H*H)
        jac(I,I)=jac(I,I)-20.0/(6.0*H*H)
        jac(I,NVAR)=jac(I,NVAR)+
     &  (p13_gx(option,X(I))
     &  +(p13_gx(option,X(IMM))-4.0*p13_gx(option,X(I)))/12.0)
        jac(I,I)=jac(I,I)+X(NVAR)*p13_gp(option,X(I))
        jac(I,IMM)=jac(I,IMM)+X(NVAR)*p13_gp(option,X(IMM))/12.0
        jac(I,I)=jac(I,I)-4.0*X(NVAR)*p13_gp(option,X(I))/12.0
      end do

      I=NEQN-(M-1)
      IP1=I+1
      IMMP1=I-M+1
      jac(I,IMMP1)=jac(I,IMMP1)+1.0/(6.0*H*H)
      jac(I,IP1)=jac(I,IP1)+4.0/(6.0*H*H)
      jac(I,NVAR)=jac(I,NVAR)+(p13_gx(option,X(IP1))+2.0)/12.0
      jac(I,IP1)=jac(I,IP1)+X(NVAR)*p13_gp(option,X(IP1))/12.0
      ILO=NEQN-(M-2)
      IHI=NEQN-1

      DO I=ILO,IHI
        IP1=I+1
        IM1=I-1
        IMMM1=I-M-1
        IMMP1=I-M+1
        jac(I,IMMM1)=jac(I,IMMM1)+1.0/(6.0*H*H)
        jac(I,IMMP1)=jac(I,IMMP1)+1.0/(6.0*H*H)
        jac(I,IM1)=jac(I,IM1)+4.0/(6.0*H*H)
        jac(I,IP1)=jac(I,IP1)+4.0/(6.0*H*H)
        jac(I,NVAR)=jac(I,NVAR)+
     &  (p13_gx(option,X(IM1))+p13_gx(option,X(IP1))+1.0)/12.0
        jac(I,IM1)=jac(I,IM1)+X(NVAR)*p13_gp(option,X(IM1))/12.0
        jac(I,IP1)=jac(I,IP1)+X(NVAR)*p13_gp(option,X(IP1))/12.0
      end do

      I=NEQN
      IM1=I-1
      IMMM1=I-M-1
      jac(I,IMMM1)=jac(I,IMMM1)+1.0/(6.0*H*H)
      jac(I,IM1)=jac(I,IM1)+4.0/(6.0*H*H)
      jac(I,NVAR)=jac(I,NVAR)+(p13_gx(option,X(IM1))+2.0)/12.0
      jac(I,IM1)=jac(I,IM1)+X(NVAR)*p13_gp(option,X(IM1))/12.0

      return
      END
      subroutine p13_nvar ( option, nvar )

c*********************************************************************72
c
cc P13_NVAR sets the number of variables for problem 13.
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
c    Input, integer OPTION, the option chosen for this problem.
c    For some problems, several options are available.  At least,
c    OPTION = 1 is always legal.
c
c    Output, integer NVAR, the number of variables.
c
      implicit none

      integer m
      integer nvar
      integer option

      if ( option .eq. 1 .or. option .eq. 2 ) then
        m = 7
      else if ( option .eq. 3 .or. option .eq. 4 ) then
        m = 11
      else if ( option .eq. 5 .or. option .eq. 6 ) then
        m = 15
      end if

      nvar = m * m + 1

      return
      end
      subroutine p13_option_num ( option_num )

c*********************************************************************72
c
cc P13_OPTION_NUM returns the number of options for problem 13.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    16 September 2008
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, integer OPTION_NUM, the number of options.
c
      implicit none

      integer option_num

      option_num = 6

      return
      end
      subroutine p13_start ( option, nvar, x )

c*********************************************************************72
c
cc P13_START returns a starting point for problem 13.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    19 September 2008
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer OPTION, the option index.
c
c    Input, integer NVAR, the number of variables.
c
c    Output, double precision X(NVAR), the starting point.
c
      implicit none

      integer nvar

      integer i
      integer option
      double precision x(nvar)

      if ( 1 .le. option .and. option .le. 6 ) then

        do i = 1, nvar
          x(i) = 0.0
        end do

      else

        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'P13_START - Fatal error!'
        write ( *, '(a,i8)' ) '  Unrecognized option  = ', option
        stop

      end if

      return
      end
      subroutine p13_stepsize ( option, h, hmin, hmax )

c*********************************************************************72
c
cc P13_STEPSIZE returns step sizes for problem 13.
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
c    Input, integer OPTION, the option index.
c
c    Output, double precision H, HMIN, HMAX, suggested values for the
c    initial step, the minimum step, and the maximum step.
c
      implicit none

      double precision h
      double precision hmax
      double precision hmin
      integer option

      h =     2.000D+00
      hmin =  0.001D+00
      hmax = 10.000D+00

      return
      end
      subroutine p13_title ( option, title )

c*********************************************************************72
c
cc P13_TITLE sets the title for problem 13.
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
c    Input, integer OPTION, the option index.
c
c    Output, character*(*) TITLE, the title of the problem.
c    TITLE will never be longer than 80 characters.
c
      implicit none

      integer option
      character*(*) title

      if ( option .eq. 1 ) then
        title = 'Simpson''s BVP, F(U) = EXP(U), M = 8.'
      else if ( option .eq. 2 ) then
        title = 'Simpson''s BVP, F(U) = function 2, M = 8.'
      else if ( option .eq. 3 ) then
        title = 'Simpson''s BVP, F(U) = EXP(U), M = 12.'
      else if ( option .eq. 4 ) then
        title = 'Simpson''s BVP, F(U) = function 2, M = 12.'
      else if ( option .eq. 5 ) then
        title = 'Simpson''s BVP, F(U) = EXP(U), M = 16.'
      else if ( option .eq. 6 ) then
        title = 'Simpson''s BVP, F(U) = function 2, M = 16.'
      else
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'P13_TITLE - Fatal error!'
        write ( *, '(a)' ) '  Unrecognized option number.'
        stop
      end if

      return
      end
      subroutine p14_data (DIRIPC,HFACT,HMAX,HMIN,HTAN,option,
     & IPC,IPRAM,ISIZE,IT,LIM,MAXCON,
     & MAXLIM,MAXSTP,MAXTAR,NCOL,NROW,NVAR,RELERR,XIT,RWORK)

c*********************************************************************72
c
cc P14_DATA sets parameters for Keller's boundary value function.
c
c  14.1 kellers boundary value problem
c
c  14.2 reference
c
c  h b keller,
c  numerical methods for two-point boundary value problems,
c  blaisdell publishing company, 1968
c
c  14.3 the function
c
c  the original problem describes a diffusion-kinetics system,
c  and has the form
c
c  -d/dt (t*t*d(x)*dx/dt) + t*t*g(x)=0
c  dx/dt(0)=0, x(1)=1.
c
c  the difference approximation used is
c
c  neqn points t(1)....t(neqn),
c  with h=1/(neqn-2)
c  and t(1)=-h, t(2)=0, t(3)=h, ..., t(neqn)=1.0.
c
c  first equation-
c
c  x(3)-x(1)=0.0
c
c  neqn-th equation
c
c  x(neqn)-1.0=0.0
c
c  equations j=2 through j=(neqn-1)-
c  centered at point x(j)
c
c  t**2(j-.5) * d(x(j-.5)) * (x(j)-x(j-1)) +
c  t**2(j+.5) * d(x(j+.5)) * (x(j)-x(j+1)) +
c  h**2 * t**2(j) * g(x(j))=0.0
c
c  with t(j-.5)=.5 * (t(j-1)+t(j))
c  and x(j-.5)=.5 * (x(j-1)+x(j)).
c  and the diffusion function d(x)=1 + x(nvar)/(x+alpha)**2
c  and g(x)=x/(beta*(x+gamma))
c
c  for this version alpha=beta=gamma=.1
c
c  14.4 options
c
c  option=1 it=nvar, xit=1.0, lim=nvar.
c
c  14.8 comments
c
c  trying linear starting point.
c
      DIMENSION RWORK(ISIZE)

      NVAR=12
      H=1.0/FLOAT(NVAR-3)
      DIRIPC=1.0
      HFACT=3.0
      IPC=NVAR
      IT=NVAR
      LIM=NVAR
      MAXTAR=1
      NCOL=NVAR
      NROW=NVAR
      RELERR=5.0E-3
      XIT=1.0
      NEQN=NVAR-1

      return
      END
      subroutine p14_fun ( option, nvar, x, fx )

c*********************************************************************72
c
cc P14_FUN evaluates Keller's boundary value function.
c
      integer nvar

      double precision fx(nvar-1)
      integer option
      double precision x(nvar)

      h=1.0/float(nvar-3)
      alpha=0.1
      beta=0.1
      gamma=0.1

      neqn=nvar-1
      fx(1)=x(3)-x(1)
      fx(neqn)=x(neqn)-1.0
      xnvar=x(nvar)
      neqnm1=neqn-1

      do i=2,neqnm1
        ip1=i+1
        im1=i-1
        ti=float(i-2)*h
        tim1=float(i-3)*h
        timh=0.5*(ti+tim1)
        tip1=float(im1)*h
        tiph=0.5*(ti+tip1)
        xi=x(i)
        xim1=x(im1)
        ximh=.5*(xi+xim1)
        xip1=x(ip1)
        xiph=.5*(xi+xip1)
        fact1=timh*timh*(1.0+xnvar/((ximh+alpha)**2))
        fact2=tiph*tiph*(1.0+xnvar/((xiph+alpha)**2))
        term=h*h*ti*ti*gx0014(xi)
        fx(i)=term+fact1*(xi-xim1)+fact2*(xi-xip1)
      end do

      return
      end
      function p14_gx (x)

c*********************************************************************72
c
cc P14_GX is an auxilliary function for Keller's boundary value problem.
c
      H=1.0/FLOAT(NVAR-3)
      ALPHA=0.1
      BETA=0.1
      GAMMA=0.1

      p14_gx =X/(BETA*(X+GAMMA))

      return
      END
      subroutine p14_jac ( option, nvar, x, jac )

c*********************************************************************72
c
cc P14_JAC evaluates the Keller boundary value jacobian.
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
c    Input, integer OPTION, the option index.
c
c    Input, integer NVAR, the number of variables.
c
c    Input, double precision X(NVAR), the argument of the jacobian.
c
c    Output, double precision JAC(NVAR,NVAR), the jacobian matrix evaluated
c    at X.  The NVAR-th row is not set by this routine.
c
      implicit none

      integer nvar

      integer option
      double precision jac(nvar,nvar)
      double precision x(nvar)

      H=1.0/FLOAT(NVAR-3)
      ALPHA=0.1
      BETA=0.1
      GAMMA=0.1

      NEQN=NVAR-1

      do j = 1, ncol
        do i = 1, nrow
          jac(i,j) = 0.0E+00
        end do
      end do

      jac(1,1)=-1.0
      jac(1,3)=1.0
      jac(NEQN,NEQN)=1.0
      XNVAR=X(NVAR)
      NEQNM1=NEQN-1

      DO I=2,NEQNM1

        IM1=I-1
        IP1=I+1
        TI=FLOAT(I-2)*H
        TIM1=FLOAT(I-3)*H
        TIMH=.5*(TI+TIM1)
        TIP1=FLOAT(IM1)*H
        TIPH=.5*(TI+TIP1)
        XI=X(I)
        XIM1=X(IM1)
        XIMH=.5*(XI+XIM1)
        XIP1=X(IP1)
        XIPH=.5*(XI+XIP1)
        TIMHSQ=TIMH*TIMH
        TIPHSQ=TIPH*TIPH
c
c  differentiate -timh*timh*xim1*(1.0+xnvar/((ximh+alpha)**2))
c
        jac(I,IM1)=jac(I,IM1)-TIMHSQ*(1.0+XNVAR/((XIMH+ALPHA)**2))
        jac(I,IM1)=jac(I,IM1)+TIMHSQ*XIM1*XNVAR/((XIMH+ALPHA)**3)
        jac(I,I)=jac(I,I)+TIMHSQ*XNVAR*XIM1/((XIMH+ALPHA)**3)
        jac(I,NVAR)=jac(I,NVAR)-TIMHSQ*XIM1/((XIMH+ALPHA)**2)
c
c  differentiate timh*timh*xi*(1.0+xnvar/((ximh+alpha)**2))
c
        jac(I,I)=jac(I,I)+TIMHSQ*(1.0+XNVAR/((XIMH+ALPHA)**2))
        jac(I,I)=jac(I,I)-TIMHSQ*XI*XNVAR/((XIMH+ALPHA)**3)
        jac(I,IM1)=jac(I,IM1)-TIMHSQ*XI*XNVAR/((XIMH+ALPHA)**3)
        jac(I,NVAR)=jac(I,NVAR)+TIMHSQ*XI/((XIMH+ALPHA)**2)
c
c  differentiate tiph*tiph*xi*(1.0+xnvar/((xiph+alpha)**2))
c
        jac(I,I)=jac(I,I)+TIPHSQ*(1.0+XNVAR/((XIPH+ALPHA)**2))
        jac(I,I)=jac(I,I)-TIPHSQ*XI*XNVAR/((XIPH+ALPHA)**3)
        jac(I,IP1)=jac(I,IP1)-TIPHSQ*XI*XNVAR/((XIPH+ALPHA)**3)
        jac(I,NVAR)=jac(I,NVAR)+TIPHSQ*XI/((XIPH+ALPHA)**2)
c
c  differentiate -tiph*tiph*xip1*(1.0+xnvar/((xiph+alpha)**2))
c
        jac(I,IP1)=jac(I,IP1)-TIPHSQ*(1.0+XNVAR/((XIPH+ALPHA)**2))
        jac(I,I)=jac(I,I)+TIPHSQ*XIP1*XNVAR/((XIPH+ALPHA)**3)
        jac(I,IP1)=jac(I,IP1)+TIPHSQ*XIP1*XNVAR/((XIPH+ALPHA)**3)
        jac(I,NVAR)=jac(I,NVAR)-TIPHSQ*XIP1/((XIPH+ALPHA)**2)
c
c  differentiate h*h*ti*ti*gx0014(xi)
c
        jac(I,I)=jac(I,I)+H*H*TI*TI*GP0014(XI)

      end do

      return
      END
      FUNCTION p14_gp (X)

c*********************************************************************72
c
cc P14_GP is an auxilliary function for Keller's boundary value problem.
c
      H=1.0/FLOAT(NVAR-3)
      ALPHA=0.1
      BETA=0.1
      GAMMA=0.1

      p14_gp = GAMMA/(BETA*(X+GAMMA)**2)

      return
      END
      subroutine p14_nvar ( option, nvar )

c*********************************************************************72
c
cc P14_NVAR sets the number of variables for problem 14.
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
c    Input, integer OPTION, the option chosen for this problem.
c    For some problems, several options are available.  At least,
c    OPTION = 1 is always legal.
c
c    Output, integer NVAR, the number of variables.
c
      implicit none

      integer nvar
      integer option

      nvar = 12

      return
      end
      subroutine p14_option_num ( option_num )

c*********************************************************************72
c
cc P14_OPTION_NUM returns the number of options for problem 14.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    16 September 2008
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, integer OPTION_NUM, the number of options.
c
      implicit none

      integer option_num

      option_num = 1

      return
      end
      subroutine p14_start ( option, nvar, x )

c*********************************************************************72
c
cc P14_START returns a starting point for problem 14.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    29 September 2008
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer OPTION, the option index.
c
c    Input, integer NVAR, the number of variables.
c
c    Output, double precision X(NVAR), the starting point.
c
      implicit none

      integer nvar

      double precision h
      integer i
      integer option
      double precision x(nvar)

      if ( option == 1 ) then

        x(1)  = 0.029742673007439D+00
        x(2)  = 0.029742673007439D+00
        x(3)  = 0.029742673007439D+00
        x(4)  = 0.039933250735582D+00
        x(5)  = 0.061866539016825D+00
        x(6)  = 0.101137641789028D+00
        x(7)  = 0.164623875371221D+00
        x(8)  = 0.258536575943466D+00
        x(9)  = 0.387217701462343D+00
        x(10) = 0.553103336509555D+00
        x(11) = 0.757271228030916D+00
        x(12) = 1.000000000000000D+00
        x(13) = 0.000000000000000D+00

      else

        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'P14_START - Fatal error!'
        write ( *, '(a,i8)' ) '  Unrecognized option  = ', option
        stop

      end if

      return
      end
      subroutine p14_stepsize ( option, h, hmin, hmax )

c*********************************************************************72
c
cc P14_STEPSIZE returns step sizes for problem 14.
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
c    Input, integer OPTION, the option index.
c
c    Output, double precision H, HMIN, HMAX, suggested values for the
c    initial step, the minimum step, and the maximum step.
c
      implicit none

      double precision h
      double precision hmax
      double precision hmin
      integer option

      h =    2.000D+00
      hmin = 0.001D+00
      hmax = 2.000D+00

      return
      end
      subroutine p14_title ( option, title )

c*********************************************************************72
c
cc P14_TITLE sets the title for problem 14.
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
c    Input, integer OPTION, the option index.
c
c    Output, character*(*) TITLE, the title of the problem.
c    TITLE will never be longer than 80 characters.
c
      implicit none

      integer option
      character*(*) title

      if ( option .eq. 1 ) then
        title = 'Keller''s boundary value function.'
      else
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'P14_TITLE - Fatal error!'
        write ( *, '(a)' ) '  Unrecognized option number.'
        stop
      end if

      return
      end
      subroutine p15_data (DIRIPC,HFACT,HMAX,HMIN,HTAN,option,
     & IPC,IPRAM,ISIZE,IT,LIM,MAXCON,
     & MAXLIM,MAXSTP,MAXTAR,NCOL,NROW,NVAR,RELERR,XIT,RWORK)

c*********************************************************************72
c
cc P15_DATA sets parameters for the trigger circuit function.
c
c  The trigger circuit problem.
c
c  Reference
c
c  g poenisch, h schwetlick,
c  computing turning points of curves implicitly defined
c    by nonlinear equations depending on a parameter.
c  computing 26, 1981, pages 107-121.
c
c  The function
c
c  the current flow of a trigger circuit with an operational
c  amplifier is modeled by these equations.  the variables
c  x1 through x7 are voltages, and x6 is the output voltage,
c  x7 the input voltage.
c
c  the function has the form
c
c  f(x)=ax+phi(x)
c
c  Options
c
c  option=1  lim=7.
c
c  Limit points
c
c  melhem lists the following limit points
c
c  (0.04936699 0.5473584 0.04944722 0.04944743 0.1292014
c   1.166020   0.6018530)
c  (0.2357780  0.6629694 0.2375982  0.2376028  0.6208333
c   9.609130   0.3228661)
c
      DIMENSION RWORK(ISIZE)
      COMMON /AUXMEM/ ARRAY(6,7)

      NVAR=7
      DIRIPC=1.0
      HFACT=3.0
      IPC=3
      IT=0
      LIM=7
      MAXCON=30
      MAXLIM=2
      MAXSTP=30
      MAXTAR=30
      NCOL=NVAR
      NROW=NVAR
      RELERR=1.0E-5
      XIT=0.0
c
c  set up array
c
      ARRAY(1,1)=.04534887
      ARRAY(2,1)=-.025641026
      ARRAY(3,1)=-.0001
      ARRAY(4,1)=0.0
      ARRAY(5,1)=0.0
      ARRAY(6,1)=0.0
      ARRAY(1,2)=-.025641026
      ARRAY(2,2)=.125641026
      ARRAY(3,2)=0.0
      ARRAY(4,2)=0.0
      ARRAY(5,2)=0.0
      ARRAY(6,2)=-.1
      ARRAY(1,3)=-.0001
      ARRAY(2,3)=0.0
      ARRAY(3,3)=.039315686
      ARRAY(4,3)=-.039215686
      ARRAY(5,3)=0.0
      ARRAY(6,3)=0.0
      ARRAY(1,4)=0.0
      ARRAY(2,4)=0.0
      ARRAY(3,4)=-.039215686
      ARRAY(4,4)=2.652119
      ARRAY(5,4)=-1.0
      ARRAY(6,4)=0.0
      ARRAY(1,5)=0.0
      ARRAY(2,5)=0.0
      ARRAY(3,5)=0.0
      ARRAY(4,5)=-1.0
      ARRAY(5,5)=1.076923078
      ARRAY(6,5)=-.076923078
      ARRAY(1,6)=0.0
      ARRAY(2,6)=-.1
      ARRAY(3,6)=0.0
      ARRAY(4,6)=0.0
      ARRAY(5,6)=-.076923078
      ARRAY(6,6)=5.1520474
      ARRAY(1,7)=.019607844
      ARRAY(2,7)=0.0
      ARRAY(3,7)=0.0
      ARRAY(4,7)=0.0
      ARRAY(5,7)=0.0
      ARRAY(6,7)=0.0
      return
      END
      subroutine p15_fun ( option, nvar, x, fx )

c*********************************************************************72
c
cc P15_FUN evaluates the trigger circuit function.
c
      integer nvar

      double precision fx(nvar-1)
      integer option
      double precision x(nvar)

      COMMON /AUXMEM/ ARRAY(6,7)
      NEQN=NVAR-1

      DO I=1,NEQN
        FX(I)=0.0
        DO J=1,NVAR
          FX(I)=FX(I)+ARRAY(I,J)*X(J)
        end do
      end do

      X1=X(1)
      X2=X(2)
      X3=X(3)
      X5=X(5)
c
c  add nonlinear terms
c
      FX(2)=FX(2)+(5.6E-8)*(EXP(25.0*X2)-1.0)
      FX(5)=FX(5)+(5.6E-8)*(EXP(25.0*X5)-1.0)
      FX(6)=FX(6)-7.65*ATAN(1962.0*(X3-X1))/.201
      return
      END
      subroutine p15_jac ( option, nvar, x, jac )

c*********************************************************************72
c
cc P15_JAC evaluates the trigger circuit jacobian.
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
c    Input, integer OPTION, the option index.
c
c    Input, integer NVAR, the number of variables.
c
c    Input, double precision X(NVAR), the argument of the jacobian.
c
c    Output, double precision JAC(NVAR,NVAR), the jacobian matrix evaluated
c    at X.  The NVAR-th row is not set by this routine.
c
      implicit none

      integer nvar

      integer option
      double precision jac(nvar,nvar)
      double precision x(nvar)

      COMMON /AUXMEM/ ARRAY(6,7)

      do i = 1, nrow
        do j = 1, ncol
          jac(i,j) = array(i,j)
        end do
      end do

      X1=X(1)
      X2=X(2)
      X3=X(3)
      X5=X(5)
      jac(2,2)=jac(2,2)+(5.6E-8)*25.0*EXP(25.0*X2)
      jac(5,5)=jac(5,5)+(5.6E-8)*25.0*EXP(25.0*X5)
      T=(1962.0*(X3-X1))**2
      jac(6,1)=jac(6,1)+(7.65*1962.0)/(.201*(1.+T))
      jac(6,3)=jac(6,3)-(7.65*1962.0)/(.201*(1.0+T))
      return
      END
      subroutine p15_nvar ( option, nvar )

c*********************************************************************72
c
cc P15_NVAR sets the number of variables for the trigger circuit function.
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
c    Input, integer OPTION, the option chosen for this problem.
c    For some problems, several options are available.  At least,
c    OPTION = 1 is always legal.
c
c    Output, integer NVAR, the number of variables.
c
      implicit none

      integer nvar
      integer option

      nvar = 7

      return
      end
      subroutine p15_option_num ( option_num )

c*********************************************************************72
c
cc P15_OPTION_NUM returns the number of options for the trigger circuit function.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    16 September 2008
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, integer OPTION_NUM, the number of options.
c
      implicit none

      integer option_num

      option_num = 1

      return
      end
      subroutine p15_start ( option, nvar, x )

c*********************************************************************72
c
cc P15_START returns a starting point for the trigger circuit function.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    19 September 2008
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer OPTION, the option index.
c
c    Input, integer NVAR, the number of variables.
c
c    Output, double precision X(NVAR), the starting point.
c
      implicit none

      integer nvar

      integer i
      integer option
      double precision x(nvar)

      if ( option .eq. 1 ) then

        do i = 1, nvar
          x(i) = 0.0
        end do

      else

        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'P15_START - Fatal error!'
        write ( *, '(a,i8)' ) '  Unrecognized option  = ', option
        stop

      end if

      return
      end
      subroutine p15_stepsize ( option, h, hmin, hmax )

c*********************************************************************72
c
cc P15_STEPSIZE returns step sizes for the trigger circuit function.
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
c    Input, integer OPTION, the option index.
c
c    Output, double precision H, HMIN, HMAX, suggested values for the
c    initial step, the minimum step, and the maximum step.
c
      implicit none

      double precision h
      double precision hmax
      double precision hmin
      integer option

      h =    0.300D+00
      hmin = 0.001D+00
      hmax = 0.600D+00

      return
      end
      subroutine p15_title ( option, title )

c*********************************************************************72
c
cc P15_TITLE sets the title for the trigger circuit function.
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
c    Input, integer OPTION, the option index.
c
c    Output, character*(*) TITLE, the title of the problem.
c    TITLE will never be longer than 80 characters.
c
      implicit none

      integer option
      character*(*) title

      if ( option .eq. 1 ) then
        title = 'The trigger circuit function.'
      else
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'P15_TITLE - Fatal error!'
        write ( *, '(a)' ) '  Unrecognized option number.'
        stop
      end if

      return
      end
      subroutine p16_data (DIRIPC,HFACT,HMAX,HMIN,HTAN,option,
     & IPC,IPRAM,ISIZE,IT,LIM,MAXCON,
     & MAXLIM,MAXSTP,MAXTAR,NCOL,NROW,NVAR,RELERR,XIT,RWORK)

c*********************************************************************72
c
cc P16_DATA sets parameters for the Moore-Spence chemical reaction function.
c
c  moore-spence chemical reaction integral
c
c  Reference
c
c  g moore, a spence,
c  the calculation of turning points of nonlinear equations,
c  siam journal of numerical analysis, volume 17, number 4, august 1980.
c  pages 567-576.
c
c  The function
c
c  the equations describe the heat and mass transfer in a plate-
c  shaped porous catalyst.  the integral equation is discretized
c  using the trapezoidal rule.
c
c  the problem is of the form
c
c  y(s) - x(nvar)*integral(0 to 1) (k(s,t)*h(y(t))dt) - 1 = 0
c
c  with k(s,t)=max(s,t)-1
c  h(z)=z*exp(gamma*beta*(1-z)/(1+beta*(1-z)))
c  with gamma=20, beta=0.4.
c
c  Options
c
c  option=1 it=nvar, xit=1.0.
c
c  Limit points
c
c  melhem lists the limit points
c
c  x(nvar)=0.1375390,  x(16)=0.8524311,
c  x(nvar)=0.07791579, x(16)=0.4657826.
c
c  computational results
c
c  x(nvar)=0.1286312,  x(16)=0.8977113,
c  x(nvar)=0.0926850,  x(16)=0.2956740.
c
      DIMENSION RWORK(ISIZE)

      NVAR=33

      DIRIPC=1.0
      HFACT=2.0
      IPC=NVAR
      IT=NVAR
      LIM=NVAR
      MAXCON=45
      MAXLIM=2
      MAXSTP=45
      MAXTAR=45
      NCOL=NVAR
      NROW=NVAR
      RELERR=1.0E-4
      XIT=1.0

      return
      END
      subroutine p16_fun ( option, nvar, x, fx )

c*********************************************************************72
c
cc P16_FUN evaluates the Moore-Spence chemical reaction function.
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
c    Input, integer OPTION, the option index.
c
c    Input, integer NVAR, the number of variables.
c
c    Input, double precision X(NVAR), the argument of the jacobian.
c
c    Output, double precision FX(NVAR), contains F(X) in entries 1 to NVAR-1.
c
      implicit none

      integer nvar

      double precision arg
      double precision beta
      parameter ( beta = 0.4D+00 )
      double precision fx(nvar)
      double precision gamma
      parameter ( gamma = 20.0D+00 )
      double precision hxj
      integer i
      integer j
      integer option
      double precision sum
      double precision x(nvar)

      do i = 1, nvar - 1
        sum = 0.0D+00
        do j = 1, nvar - 1
          arg = beta * gamma * ( 1.0D+00 - x(j) )
     &      / ( 1.0D+00 + beta * ( 1.0D+00 - x(j) ) )
          hxj = x(j) * exp ( arg )
          if ( j .eq. 1 ) then
            hxj = 0.5D+00 * hxj
          end if
          sum = sum + dble ( max ( i, j ) - neqn - 1 ) * hxj
     &      / dble ( nvar - 1 )
        end do
        fx(i) = x(i) - ( x(nvar) * sum / dble ( nvar - 1 ) ) - 1.0D+00
      end do

      return
      end
      subroutine p16_jac ( option, nvar, x, jac )

c*********************************************************************72
c
cc P16_JAC evaluates the Moore-Spence chemical reaction jacobian.
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
c    Input, integer OPTION, the option index.
c
c    Input, integer NVAR, the number of variables.
c
c    Input, double precision X(NVAR), the argument of the jacobian.
c
c    Output, double precision JAC(NVAR,NVAR), the jacobian matrix evaluated
c    at X.  The NVAR-th row is not set by this routine.
c
      implicit none

      integer nvar

      double precision beta
      parameter ( beta = 0.4D+00 )
      double precision gamma
      parameter ( gamma = 20.0D+00 )
      integer i
      integer j
      double precision jac(nvar,nvar)
      integer option
      double precision x(nvar)

      do j = 1, ncol
        do i = 1, nrow
          jac(i,j) = 0.0E+00
        end do
      end do

      NEQN=NVAR-1

      DO I=1,NEQN

        XI=X(I)
        FACTI=EXP(BETA*GAMMA*(1.0-XI/(1.0+BETA*(1.0-XI))))
        HPXI=FACTI*(1.0-(BETA*GAMMA*XI)/((1.0+BETA*(1.0-XI))**2))
        IF(I.EQ.1)HPXI=0.5*HPXI
        T=-X(NVAR)*FLOAT(I-NEQN-1)

        DO J=1,I
          XJ=X(J)
          FACTJ=EXP(BETA*GAMMA*(1.0-XJ/(1.0+BETA*(1.0-XJ))))
          HPXJ=FACTJ*(1.0-(BETA*GAMMA*XJ)/((1.0+BETA*(1.0-XJ))**2))
          IF(J.EQ.1)HPXJ=0.5*HPXJ
          jac(I,J)=T*HPXJ/FLOAT(NEQN*NEQN)
          IF(J.LT.I)jac(J,I)=T*HPXI/FLOAT(NEQN*NEQN)
        end do

      end do

      DO I=1,NEQN
        jac(I,I)=jac(I,I)+1.0
      end do

      DO I=1,NEQN
        SUM=0.0
        DO J=1,NEQN
          XJ=X(J)
          HXJ=XJ*EXP(BETA*GAMMA*(1.0-XJ)/(1.0+BETA*(1.0-XJ)))
          IF(J.EQ.1)HXJ=0.5*HXJ
          INDEX=MAX0(I,J)
          SUM=SUM+FLOAT(INDEX-NEQN-1)*HXJ
        end do
        jac(I,NVAR)=-SUM/FLOAT(NEQN*NEQN)
      end do

      return
      END
      subroutine p16_nvar ( option, nvar )

c*********************************************************************72
c
cc P16_NVAR sets the number of variables for the Moore-Spence chemical reaction function.
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
c    Input, integer OPTION, the option chosen for this problem.
c    For some problems, several options are available.  At least,
c    OPTION = 1 is always legal.
c
c    Output, integer NVAR, the number of variables.
c
      implicit none

      integer nvar
      integer option

      nvar = 33

      return
      end
      subroutine p16_option_num ( option_num )

c*********************************************************************72
c
cc P16_OPTION_NUM returns the number of options for the Moore-Spence chemical reaction function.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    16 September 2008
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, integer OPTION_NUM, the number of options.
c
      implicit none

      integer option_num

      option_num = 1

      return
      end
      subroutine p16_start ( option, nvar, x )

c*********************************************************************72
c
cc P16_START returns a starting point for the Moore-Spence chemical reaction function.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    19 September 2008
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer OPTION, the option index.
c
c    Input, integer NVAR, the number of variables.
c
c    Output, double precision X(NVAR), the starting point.
c
      implicit none

      integer nvar

      integer option
      double precision x(nvar)

      if ( option .eq. 1 ) then

        do i = 1, nvar - 1
          x(i) = 1.0
        end do
        x(nvar) = 0.0D+00

      else

        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'P16_START - Fatal error!'
        write ( *, '(a,i8)' ) '  Unrecognized option  = ', option
        stop

      end if

      return
      end
      subroutine p16_stepsize ( option, h, hmin, hmax )

c*********************************************************************72
c
cc P16_STEPSIZE returns step sizes for the Moore-Spence chemical reaction function.
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
c    Input, integer OPTION, the option index.
c
c    Output, double precision H, HMIN, HMAX, suggested values for the
c    initial step, the minimum step, and the maximum step.
c
      implicit none

      double precision h
      double precision hmax
      double precision hmin
      integer option

      h =    0.200D+00
      hmin = 0.001D+00
      hmax = 2.000D+00

      return
      end
      subroutine p16_title ( option, title )

c*********************************************************************72
c
cc P16_TITLE sets the title for the Moore-Spence chemical reaction function.
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
c    Input, integer OPTION, the option index.
c
c    Output, character*(*) TITLE, the title of the problem.
c    TITLE will never be longer than 80 characters.
c
      implicit none

      integer option
      character*(*) title

      if ( option .eq. 1 ) then
        title = 'The Moore-Spence chemical reaction function.'
      else
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'P16_TITLE - Fatal error!'
        write ( *, '(a)' ) '  Unrecognized option number.'
        stop
      end if

      return
      end
      subroutine p17_data (DIRIPC,HFACT,HMAX,HMIN,HTAN,option,
     & IPC,IPRAM,ISIZE,IT,LIM,MAXCON,
     & MAXLIM,MAXSTP,MAXTAR,NCOL,NROW,NVAR,RELERR,XIT,RWORK)

c*********************************************************************72
c
cc P17_DATA sets parameters for the Bremerman propane combustion system.
c
c  bremerman propane combustion system
c
c  Reference:
c
c  h j bremerman,
c  calculation of equilibrium points for models of ecological and
c    chemical systems,
c  in  proceedings of a conference on the applications of undergraduate
c    mathematics in the engineering, life, managerial and
c    social sciences,
c  georgia institute of technology, june 1973, pages 198-217.
c
c  k l hiebert,
c  a comparison of software which solves systems of nonlinear equations,
c  sandia national laboratory, albuquerque, new mexico,
c  report sand 80-0181, 1980.
c
c  The function
c
c  the equations describe the combustion of propane (c3-h4)
c  in air (o2 and n2) at 2200 degrees fahrenheit.
c  the substances created are
c
c  x1  co2
c  x2  h2o
c  x3  n2
c  x4  co
c  x5  h2
c  x6  h
c  x7  oh
c  x8  o
c  x9  no
c  x10 o2
c
c  with auxilliary variables
c
c  x11 amount of air - x11/2 moles of o2, 2x11 moles of n2.
c  x12 pressure in atmospheres.
c
c  the mass balance and reaction equations become, once square
c  roots are eliminated via squaring,
c
c  fx1  = x1 + x4 - 3.0
c  fx2  = 2*x1 + x2 + x4 + x7 + x8 + x9 + 2*x10 - x12
c  fx3 = 2*x2 + 2*x5 + x6 + x7 - 8.0
c  fx4 = 2*x3 + x9 - 4*x12
c  fx5 = x1*x5 - 0.193*x2*x4
c  fx6 = x11*x1*x6*x6 - (0.002597)**2 *x2*x4*sum(x)
c  fx7 = x11*x4*x7*x7 - (0.003448)**2 *x1*x2*sum(x)
c  fx8 = x11*x4*x8    - (0.0001799)*x1*sum(x)
c  fx9 = x11*x4*x4*x9*x9 - (0.0002155)**2 *x1*x1*x3*sum(x)
c  fx10= x11*x4*x4*x10 - 0.00003846 *x1*x1*sum(x)
c  fx11= x(ifix1)-val1
c
c  where sum(x) = sum(i=1 to 10) x(i)
c  and ifix1=11 or 12 with val1 some fixed value.
c
c  Options
c
c  option=1, ifix1=11, val1=1.0, it=12, xit=15.0 (fixed pressure)
c  option=2, ifix1=12, val1=5.0, it=11, xit=10.0 (fixed concentrations)
c
      DIMENSION RWORK(ISIZE)
      COMMON /AUXMEM/ IFIX1,VAL1
      NVAR=12


      DIRIPC=1.0
      HFACT=2.0
      IPC=IFREE
      IT=IFREE
      LIM=IFREE
      MAXCON=20
      MAXLIM=20
      MAXSTP=20
      MAXTAR=1
      NCOL=NVAR
      NROW=NVAR
      RELERR=1.0E-4
      IF(option.EQ.1)XIT=15.0
      IF(option.EQ.2)XIT=10.0

      return
      END
      subroutine p17_fun ( option, nvar, x, fx )

c*********************************************************************72
c
cc P17_FUN evaluates the Bremerman propane combustion function.
c
      integer nvar

      double precision fx(nvar-1)
      integer option
      double precision x(nvar)

      COMMON /AUXMEM/ IFIX1,VAL1
      X1=X(1)
      X2=X(2)
      X3=X(3)
      X4=X(4)
      X5=X(5)
      X6=X(6)
      X7=X(7)
      X8=X(8)
      X9=X(9)
      X10=X(10)
      X11=X(11)
      X12=X(12)

      XSUM=0.0
      DO I=1,10
        XSUM=XSUM+X(I)
      end do

      FX(1)=X1+X4-3.0
      FX(2)=2.0*X1+X2+X4+X7+X8+X9+2.0*X10-X12
      FX(3)=2.0*X2+2.0*X5+X6+X7-8.0
      FX(4)=2.0*X3+X9-4.0*X12
      FX(5)=X1*X5-0.193*X2*X4
      FX(6)=X11*X1*X6*X6-(0.002597)**2*X2*X4*XSUM
      FX(7)=X11*X4*X7*X7-(0.003448)**2*X1*X2*XSUM
      FX(8)=X11*X4*X8-0.0001799*X1*XSUM
      FX(9)=X11*X4*X4*X9*X9-(0.0002155)**2*X1*X1*X3*XSUM
      FX(10)=X11*X4*X4*X10-0.00003846*X1*X1*XSUM
      FX(11)=X(IFIX1)-VAL1
      return
      END
      subroutine p17_jac ( option, nvar, x, jac )

c*********************************************************************72
c
cc P17_JAC evaluates the Bremerman propane combustion jacobian.
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
c    Input, integer OPTION, the option index.
c
c    Input, integer NVAR, the number of variables.
c
c    Input, double precision X(NVAR), the argument of the jacobian.
c
c    Output, double precision JAC(NVAR,NVAR), the jacobian matrix evaluated
c    at X.  The NVAR-th row is not set by this routine.
c
      implicit none

      integer nvar

      integer option
      double precision jac(nvar,nvar)
      double precision x(nvar)

      COMMON /AUXMEM/ IFIX1,VAL1

      do j = 1, ncol
        do i = 1, nrow
          jac(i,j) = 0.0E+00
        end do
      end do

      X1=X(1)
      X2=X(2)
      X3=X(3)
      X4=X(4)
      X5=X(5)
      X6=X(6)
      X7=X(7)
      X8=X(8)
      X9=X(9)
      X10=X(10)
      X11=X(11)
      X12=X(12)
      XSUM=0.0
      DO I=1,10
        XSUM=XSUM+X(I)
      end do
      jac(1,1)=1.0
      jac(1,4)=1.0
      jac(2,1)=2.0
      jac(2,2)=1.0
      jac(2,4)=1.0
      jac(2,7)=1.0
      jac(2,8)=1.0
      jac(2,9)=1.0
      jac(2,10)=2.0
      jac(2,12)=-1.0
      jac(3,2)=2.0
      jac(3,5)=2.0
      jac(3,6)=1.0
      jac(3,7)=1.0
      jac(4,3)=2.0
      jac(4,9)=1.0
      jac(4,12)=-4.0
      jac(5,1)=X5
      jac(5,2)=-0.193*X4
      jac(5,4)=-0.193*X2
      jac(5,5)=X1
      TERM=(0.002597)**2*X2*X4
      jac(6,1)=X11*X6*X6+TERM
      jac(6,2)=TERM-(0.002597)**2*X4*XSUM
      jac(6,3)=TERM
      jac(6,4)=TERM-(0.002597)**2*X2*XSUM
      jac(6,5)=TERM
      jac(6,6)=TERM+2.0*X11*X1*X6
      jac(6,7)=TERM
      jac(6,8)=TERM
      jac(6,9)=TERM
      jac(6,10)=TERM
      jac(6,11)=X1*X6*X6
      TERM=-(0.003448)**2*X1*X2
      jac(7,1)=TERM-(0.003448)**2*X2*XSUM
      jac(7,2)=TERM-(0.003448)**2*X1*XSUM
      jac(7,3)=TERM
      jac(7,4)=X11*X7*X7+TERM
      jac(7,5)=TERM
      jac(7,6)=TERM
      jac(7,7)=2.0*X11*X4*X7+TERM
      jac(7,8)=TERM
      jac(7,9)=TERM
      jac(7,10)=TERM
      jac(7,11)=X4*X7*X7
      TERM=-0.0001799*X1
      jac(8,1)=TERM-0.0001799*X1*XSUM
      jac(8,2)=TERM
      jac(8,3)=TERM
      jac(8,4)=X11*X8+TERM
      jac(8,5)=TERM
      jac(8,6)=TERM
      jac(8,7)=TERM
      jac(8,8)=X11*X4+TERM
      jac(8,9)=TERM
      jac(8,10)=TERM
      jac(8,11)=X4*X8
      TERM=-(0.0002155)**2*X1*X1*X3
      jac(9,1)=TERM-2.0*(0.0002155)**2*X1*X3*XSUM
      jac(9,2)=TERM
      jac(9,3)=TERM-(0.0002155)**2*X1*X1*XSUM
      jac(9,4)=2.0*X11*X4*X9*X9+TERM
      jac(9,5)=TERM
      jac(9,6)=TERM
      jac(9,7)=TERM
      jac(9,8)=TERM
      jac(9,9)=2.0*X11*X4*X4*X9+TERM
      jac(9,10)=TERM
      jac(9,11)=X4*X4*X9*X9
      TERM=-0.00003846*X1*X1
      jac(10,1)=TERM-2.0*(0.00003846)*X1*XSUM
      jac(10,2)=TERM
      jac(10,3)=TERM
      jac(10,4)=2.0*X11*X4*X10+TERM
      jac(10,5)=TERM
      jac(10,6)=TERM
      jac(10,7)=TERM
      jac(10,8)=TERM
      jac(10,9)=TERM
      jac(10,10)=X11*X4*X4+TERM
      jac(10,11)=X4*X4*X10
      jac(11,IFIX1)=1.0
      return
      END
      subroutine p17_nvar ( option, nvar )

c*********************************************************************72
c
cc P17_NVAR sets the number of variables for the Bremerman propane combustion system.
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
c    Input, integer OPTION, the option chosen for this problem.
c    For some problems, several options are available.  At least,
c    OPTION = 1 is always legal.
c
c    Output, integer NVAR, the number of variables.
c
      implicit none

      integer nvar
      integer option

      nvar = 12

      return
      end
      subroutine p17_option_num ( option_num )

c*********************************************************************72
c
cc P17_OPTION_NUM returns the number of options for the Bremerman propane combustion system.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    16 September 2008
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, integer OPTION_NUM, the number of options.
c
      implicit none

      integer option_num

      option_num = 2

      return
      end
      subroutine p17_start ( option, nvar, x )

c*********************************************************************72
c
cc P17_START returns a starting point for the Bremerman propane combustion system.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    19 September 2008
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer OPTION, the option index.
c
c    Input, integer NVAR, the number of variables.
c
c    Output, double precision X(NVAR), the starting point.
c
      implicit none

      integer nvar

      integer option
      double precision x(nvar)

      if ( 1 .le. option .and. option .le. 2 ) then

        x(1) =  0.3564320
        x(2) =  1.636071
        x(3) =  9.999810
        x(4) =  2.643568
        x(5) =  2.341926
        x(6) =  0.3732447E-01
        x(7) =  0.6681509E-02
        x(8) =  0.4128999E-03
        x(9) =  0.3790901E-03
        x(10) = 0.1190167E-04
        x(11) = 1.0
        x(12) = 5.0

      else

        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'P17_START - Fatal error!'
        write ( *, '(a,i8)' ) '  Unrecognized option  = ', option
        stop

      end if

      return
      end
      subroutine p17_stepsize ( option, h, hmin, hmax )

c*********************************************************************72
c
cc P17_STEPSIZE returns step sizes for the Bremerman propane combustion system.
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
c    Input, integer OPTION, the option index.
c
c    Output, double precision H, HMIN, HMAX, suggested values for the
c    initial step, the minimum step, and the maximum step.
c
      implicit none

      double precision h
      double precision hmax
      double precision hmin
      integer option

      h =    1.000D+00
      hmin = 0.001D+00
      hmax = 2.000D+00

      return
      end
      subroutine p17_title ( option, title )

c*********************************************************************72
c
cc P17_TITLE sets the title for the Bremerman propane combustion system.
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
c    Input, integer OPTION, the option index.
c
c    Output, character*(*) TITLE, the title of the problem.
c    TITLE will never be longer than 80 characters.
c
      implicit none

      integer option
      character*(*) title

      if ( option .eq. 1 ) then
        title = 'Bremermann Propane Combustion System, fixed pressure.'
      else if ( option .eq. 2 ) then
        title =
     &    'Bremermann Propane Combustion System, fixed concentration.'
      else
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'P17_TITLE - Fatal error!'
        write ( *, '(a)' ) '  Unrecognized option number.'
        stop
      end if

      return
      end
      subroutine p18_data (DIRIPC,HFACT,HMAX,HMIN,HTAN,option,
     & IPC,IPRAM,ISIZE,IT,LIM,MAXCON,
     & MAXLIM,MAXSTP,MAXTAR,NCOL,NROW,NVAR,RELERR,XIT,RWORK)

c*********************************************************************72
c
cc P18_DATA sets parameters for the semiconductor function.
c
c  Semiconductor problem.
c
c  Reference:
c
c  s j polak, et al
c  a continuation method for the calculation of electrostatic
c    potentials in semiconductors,
c  n v philips gloeilampen-fabrieken,
c  eindhoven, the netherlands, tech report isa-tis/card, 1979.
c
c  cor den heijer, werner c rheinboldt,
c  on steplength algorithms for a class of continuation methods,
c  siam journal of numerical analysis 18, 1981, pages 925-947.
c
c  The function
c
c  the problem is a two point boundary value problem of the form
c
c  -u'' = f(u,t,x(nvar)) for a.lt.t.lt.b, u(a)=ua, u(b)=ub.
c
c  with f(u,t,x(nvar)) = x(nvar)*ca*exp(x(nvar)*beta*(x(nvar)*ua-u))
c                       -x(nvar)*cb*exp(x(nvar)*beta*(u-x(nvar)*ub))
c                       +x(nvar)*d(t)
c
c  with d(t) = -ca for t.le.0,
c            =  cb for t.gt.0.
c
c  Options
c
c  option=1 it=nvar, xit=1.0.
c
c  Target point
c
c  the target point x(nvar)=1.0 is the solution of the embedded problem.
c
c  24.95389 25.01537 24.99566 25.00044 24.99995
c  25.00001 25.00000 25.00000 25.00000  1.000000
c
      DIMENSION RWORK(ISIZE)
      COMMON /AUXMEM/ DX(10),XX(10),UA,UB,BETA,CA,CB

      NVAR=10
      UA=0.0
      UB=25.0
      BETA=20.0
      CA=1.0E6
      CB=1.0E7
      NEQN=NVAR-1
      NEQNM1=NVAR-2
      XX(1)=2.00E-3
      XX(2)=4.00E-3
      XX(3)=4.50E-3
      XX(4)=4.75E-3
      XX(5)=5.00E-3
      XX(6)=5.25E-3
      XX(7)=5.50E-3
      XX(8)=6.00E-3
      XX(9)=8.00E-3
      XX(10)=10.0E-3
      DX(1)=2.00E-3
      DX(2)=2.00E-3
      DX(3)=0.50E-3
      DX(4)=0.25E-3
      DX(5)=0.25E-3
      DX(6)=0.25E-3
      DX(7)=0.25E-3
      DX(8)=0.50E-3
      DX(9)=2.00E-3
      DX(10)=2.00E-3
      DIRIPC=1.0
      HFACT=2.0
      IPC=NVAR
      IT=NVAR
      LIM=0
      MAXCON=30
      MAXLIM=30
      MAXSTP=30
      MAXTAR=1
      NCOL=NVAR
      NROW=NVAR
      RELERR=1.0E-3
      XIT=1.0
      return
      END
      subroutine p18_fun ( option, nvar, x, fx )

c*********************************************************************72
c
cc P18_FUN evaluates the semiconductor function.
c
      integer nvar

      double precision fx(nvar-1)
      integer option
      double precision x(nvar)

      DIMENSION F(20)
      COMMON /AUXMEM/ DX(10),XX(10),UA,UB,BETA,CA,CB
      NEQN=NVAR-1
      I1HI=2*NVAR-1

      DO I1=1,I1HI
        IF(I1.EQ.1)XST=0.5*(UA*X(NVAR)+X(1))
        I1H=I1/2
        I1P1H=(I1+1)/2
        IF(I1.GT.1.AND.I1.LE.2*NEQN)XST=0.5*(X(I1H)+X(I1P1H))
        IF(I1.EQ.2*NEQN+1)XST=0.5*(X(NEQN)+UB*X(NVAR))
        EM1=BETA*X(NVAR)*(UA*X(NVAR)-XST)
        EM2=BETA*X(NVAR)*(XST-UB*X(NVAR))
        IF(I1.EQ.1)HX=.5*XX(1)
        IF(I1.GT.1)HX=.5*(XX(I1H)+XX(I1P1H))
        DTEMP=-CA
        IF(HX.GT.0.0)DTEMP=CB
        F(I1)=(CA*EXP(EM1)-CB*EXP(EM2)+DTEMP)*X(NVAR)
      end do

      DO I=1,NEQN
        J=2*I
        JM1=J-1
        JP1=J+1
        IP1=I+1
        IM1=I-1
        XI=X(I)
        XIP1=X(IP1)
        IF(I.EQ.NEQN)XIP1=UB*X(NVAR)
        IF(I.GT.1)XIM1=X(IM1)
        IF(I.EQ.1)XIM1=UA*X(NVAR)
        RK=DX(I)*(2.0*F(JM1)+F(J))/6.0+DX(IP1)*(F(J)+2.0*F(JP1))/6.0
        FX(I)=(XI-XIM1)/DX(I)+(XI-XIP1)/DX(IP1)-RK
      end do

      return
      END
      subroutine p18_jac ( option, nvar, x, jac )

c*********************************************************************72
c
cc P18_JAC evaluates the semiconductor jacobian.
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
c    Input, integer OPTION, the option index.
c
c    Input, integer NVAR, the number of variables.
c
c    Input, double precision X(NVAR), the argument of the jacobian.
c
c    Output, double precision JAC(NVAR,NVAR), the jacobian matrix evaluated
c    at X.  The NVAR-th row is not set by this routine.
c
      implicit none

      integer nvar

      integer option
      double precision jac(nvar,nvar)
      double precision x(nvar)

      DIMENSION FV(20),FT(20)
      COMMON /AUXMEM/ DX(10),XX(10),UA,UB,BETA,CA,CB

      do j = 1, ncol
        do i = 1, nrow
          jac(i,j) = 0.0E+00
        end do
      end do

      NEQN=NVAR-1
      NEQNM1=NVAR-2
      I1HI=2*NVAR-1
      DO I1=1,I1HI
        IF(I1.EQ.1)XST=0.5*(UA*X(NVAR)+X(1))
        I1H=I1/2
        I1P1H=(I1+1)/2
        IF(I1.GT.1.AND.I1.LE.2*NEQN)XST=0.5*(X(I1H)+X(I1P1H))
        IF(I1.EQ.2*NEQN+1)XST=0.5*(X(NEQN)+UB*X(NVAR))
        EM1=BETA*X(NVAR)*(UA*X(NVAR)-XST)
        EM2=BETA*X(NVAR)*(XST-UB*X(NVAR))
        IF(I1.EQ.1)HX=.5*XX(1)
        IF(I1.GT.1)HX=.5*(XX(I1H)+XX(I1P1H))
        DTEMP=-CA
        IF(HX.GT.0.0)DTEMP=CB
        FV(I1)=-BETA*X(NVAR)*X(NVAR)*(CA*EXP(EM1)+CB*EXP(EM2))
        FT(I1)=(CA+BETA*X(NVAR)*CA*(2.0*UA*X(NVAR)-XST))*EXP(EM1)
     &        -(CB+BETA*X(NVAR)*CB*(XST-2.0*UB*X(NVAR)))*EXP(EM2)
     &        +DTEMP
      end do

      DO I=1,NEQN
        J=2*I
        JM1=J-1
        JP1=J+1
        IP1=I+1
        IM1=I-1
        RMD=DX(I)*(FV(JM1)+FV(J))/6.0
     &         +DX(IP1)*(FV(J)+FV(JP1))/6.0
        RMH=DX(I)*FV(JM1)/6.0
        jac(I,NVAR)=-DX(I)*(2.0*FT(JM1)+FT(J))/6.0
     &       -DX(IP1)*(FT(J)+2.0*FT(JP1))/6.0
        jac(I,I)=1.0/DX(I)+1.0/DX(IP1)-RMD
        IF(I.GT.1)jac(IM1,I)=-1.0/DX(I)-RMH
        IF(I.GT.1)jac(I,IM1)=-1.0/DX(I)-RMH
      end do

      RMH=DX(1)*FV(1)/6.0
      jac(1,NVAR)=jac(1,NVAR)-UA/DX(1)-RMH*UA
      INDEX=2*NVAR-1
      RMH=DX(NVAR)*FV(INDEX)/6.0
      jac(NEQN,NVAR)=jac(NEQN,NVAR)-UB/DX(NVAR)-RMH*UB
      return
      END
      subroutine p18_nvar ( option, nvar )

c*********************************************************************72
c
cc P18_NVAR sets the number of variables for the semiconductor function.
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
c    Input, integer OPTION, the option chosen for this problem.
c    For some problems, several options are available.  At least,
c    OPTION = 1 is always legal.
c
c    Output, integer NVAR, the number of variables.
c
      implicit none

      integer nvar
      integer option

      nvar = 10

      return
      end
      subroutine p18_option_num ( option_num )

c*********************************************************************72
c
cc P18_OPTION_NUM returns the number of options for the semiconductor function.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    16 September 2008
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, integer OPTION_NUM, the number of options.
c
      implicit none

      integer option_num

      option_num = 1

      return
      end
      subroutine p18_start ( option, nvar, x )

c*********************************************************************72
c
cc P18_START returns a starting point for the semiconductor function.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    19 September 2008
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer OPTION, the option index.
c
c    Input, integer NVAR, the number of variables.
c
c    Output, double precision X(NVAR), the starting point.
c
      implicit none

      integer nvar

      integer i
      integer option
      double precision x(nvar)

      if ( option .eq. 1 ) then

        do i = 1, nvar
          x(i) = 0.0
        end do

      else

        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'P18_START - Fatal error!'
        write ( *, '(a,i8)' ) '  Unrecognized option  = ', option
        stop

      end if

      return
      end
      subroutine p18_stepsize ( option, h, hmin, hmax )

c*********************************************************************72
c
cc P18_STEPSIZE returns step sizes for the semiconductor function.
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
c    Input, integer OPTION, the option index.
c
c    Output, double precision H, HMIN, HMAX, suggested values for the
c    initial step, the minimum step, and the maximum step.
c
      implicit none

      double precision h
      double precision hmax
      double precision hmin
      integer option

      h =    2.500D+00
      hmin = 0.001D+00
      hmax = 5.000D+00

      return
      end
      subroutine p18_title ( option, title )

c*********************************************************************72
c
cc P18_TITLE sets the title for the semiconductor function.
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
c    Input, integer OPTION, the option index.
c
c    Output, character*(*) TITLE, the title of the problem.
c    TITLE will never be longer than 80 characters.
c
      implicit none

      integer option
      character*(*) title

      if ( option .eq. 1 ) then
        title = 'The semiconductor function.'
      else
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'P18_TITLE - Fatal error!'
        write ( *, '(a)' ) '  Unrecognized option number.'
        stop
      end if

      return
      end
      subroutine p19_data (DIRIPC,HFACT,HMAX,HMIN,HTAN,option,
     & IPC,IPRAM,ISIZE,IT,LIM,MAXCON,
     & MAXLIM,MAXSTP,MAXTAR,NCOL,NROW,NVAR,RELERR,XIT,XR)

c*********************************************************************72
c
cc P19_DATA sets parameters for the nitric absorption problem.
c
c  The nitric acid absorption problem
c
c  Reference
c
c  t. w. copeman,
c  air products and chemicals, inc.
c  box 538,
c  allentown, pa 18105.
c
c  glossary
c
c  con    - physical equilibrium constants for the reagents.
c  flow   - flow rates for the five reagants in moles/hour.
c  press  - pressure in atmospheres.
c  temper - temperature in degrees kelvin.
c  xr     - starting value containing an approximate solution
c           of the equations.
c
c  the vector x contains the following quantities-
c
c  x(1)  = relative liquid concentration of no2.
c  x(2)  = relative liquid concentration of n2o4.
c  x(3)  = relative liquid concentration of no.
c  x(4)  = relative liquid concentration of h2o.
c  x(5)  = relative liquid concentration of hno3.
c  x(6)  = relative vapor concentration of no2.
c  x(7)  = relative vapor concentration of n2o4.
c  x(8)  = relative vapor concentration of no.
c  x(9)  = relative vapor concentration of h2o.
c  x(10) = relative vapor concentration of hno3.
c  x(11) = total number of moles of liquid.
c  x(12) = total number of moles of vapor.
c  x(13) = lambda, a multiplier for the flowrates.
c
c
c  The function
c
c  the equations have the following form
c
c  for i=1 to 5,
c
c  fx(i) = x(11)*x(i) + x(12)*x(i+5) - x(13)*flow(i)  (mole balance).
c
c  for i=6 to 10,
c
c  fx(i) = x(i) - con(i-5)*x(i-5)   (liquid-vapor transfer equation)
c
c  fx(11) = x(1) + x(2) + x(3) + x(4) + x(5) - 1.0
c  fx(12) = x(6) + x(7) + x(8) + x(9) + x(10) - 1.0
c
c  Options
c
c  option=1  good starting point
c  option=2  bad starting point
c
c  Target points
c
c  seek x(13)=1.0
c
      DIMENSION FTEMP(13)
      DIMENSION XBAD(13)
      DIMENSION XGOOD(13)
      DIMENSION XR(13)
      COMMON /COMM28/ CON(5),FLOW(5),FXOLD(13),PRESS,SCALE(13),
     & TEMPER
c
c  set flow rates, pressure, temperature,
c  good starting point, rough starting point
c
      FLOW(1)=10.0
      FLOW(2)=10.0
      FLOW(3)=10.0
      FLOW(4)=100.0
      FLOW(5)=100.0
      PRESS=7.0
      TEMPER=323.0
      XGOOD(1)=0.45
      XGOOD(2)=6.56
      XGOOD(3)=0.22
      XGOOD(4)=99.95
      XGOOD(5)=99.90
      XGOOD(6)=9.54
      XGOOD(7)=3.43
      XGOOD(8)=9.97
      XGOOD(9)=0.044
      XGOOD(10)=0.01
      SUMLIQ=0.0
      SUMVAP=0.0

      DO I=1,5
        SUMLIQ=SUMLIQ+XGOOD(I)
        SUMVAP=SUMVAP+XGOOD(I+5)
      end do

      DO I=1,5
        XGOOD(I)=XGOOD(I)/SUMLIQ
        XGOOD(I+5)=XGOOD(I+5)/SUMVAP
      end do

      XGOOD(11)=SUMLIQ
      XGOOD(12)=SUMVAP
      XGOOD(13)=0.0
      XBAD(1)=0.0
      XBAD(2)=10.0
      XBAD(3)=0.0
      XBAD(4)=100.0
      XBAD(5)=100.0
      XBAD(6)=10.0
      XBAD(7)=0.0
      XBAD(8)=10.0
      XBAD(9)=0.0
      XBAD(10)=0.0
      SUMLIQ=0.0
      SUMVAP=0.0
      DO I=1,5
        SUMLIQ=SUMLIQ+XBAD(I)
        SUMVAP=SUMVAP+XBAD(I+5)
      end do

      DO I=1,5
        XBAD(I)=XBAD(I)/SUMLIQ
        XBAD(I+5)=XBAD(I+5)/SUMVAP
      end do

      XBAD(11)=SUMLIQ
      XBAD(12)=SUMVAP
      XBAD(13)=0.0
c
c  write header
c
   60 CONTINUE
c
c  set istart to first call value
c
      ISTART=1
c
c  jacobian is supplied
c
      IJAC=0
c
c  parameterization option 1
c
      IPRAM=1
c
c  continuation process output option
c
      IWRITE=2
c
c  solution point output option
c
      JWRITE=1
c
c  function computation output option
c
      KWRITE=0
c
c  jacobian computation output option
c
      LWRITE=0
c
c  newton method option
c
      MODNEW=0
c
c  limit point index
c
      LIM=0
c
c  set startup quantities
c
   70 CONTINUE
c
c  number of variables
c
      NVAR=13
c
c  initial direction
c
      DIRIPC=1.0
c
c  stepsize increase factor
c
      HFACT=4.0
c
c  starting parameter
c
      IPC=NVAR
c
c  target index
c
      IT=NVAR
c
c  maximum number of steps
c
      MAXSTP=10
c
c  maximum number of continuation points
c
      MAXCON=MAXSTP
c
c  maximum number of limit points
c
      MAXLIM=MAXSTP
c
c  maximum number of target points
c
      MAXTAR=1
c
c  number of columns used for jacobian
c
      NCOL=NVAR
c
c  number of rows used for jacobian
c
      NROW=NVAR
c
c  relative error tolerance
c
      RELERR=1.0E-5
c
c  target value, seeking x(it)=xit
c
      XIT=1.0
c
c  set scaling quantities
c
      DO 99 I=1,NVAR
        AXR=ABS(XR(I))
        SCALE(I)=1.0
98      CONTINUE
        IF(AXR.LE.SCALE(I))GO TO 99
        SCALE(I)=SCALE(I)*2.0
        GO TO 98
99      CONTINUE
c
c  scale xr down
c
      call p19_sc (-1,NVAR,SCALE,XR)
c
c  set fxold
c
      DO I=1,NVAR
        FXOLD(I)=0.0
      end do

      call p19_fx (NVAR,XR,FTEMP,IERROR)

      NEQN=NVAR-1

      DO I=1,NEQN
        FXOLD(I)=FTEMP(I)
      end do

  130 CONTINUE
      return
      END
      subroutine p19_edit (ARCLXC,ARCLXF,ARCLXR,DIRIPC,FPNAME,
     & jac,FXNAME,HFACT,HMAX,HMIN,HTANCF,IERROR,IJAC,IPC,IPIVOT,
     & IPL,IPOINT,IPRAM,ISTART,IT,IWRITE,JPOINT,JWRITE,KWRITE,LIM,
     & LWRITE,MODNEW,NCOL,NROW,NUMCON,NUMLIM,NUMSTR,NUMTAR,
     & NVAR,RELERR,SLNAME,TC,TL,TR,WORK1,WORK2,WORK3,XC,XF,XIT,XR)

c*********************************************************************72
c
cc P19_EDIT is an edit routine for the nitric absorption problem.
c
      EXTERNAL FPNAME
      EXTERNAL FXNAME
      EXTERNAL SLNAME
      DIMENSION jac(NROW,NCOL)
      DIMENSION IPIVOT(NVAR)
      character*4 lab(5)
      DIMENSION SUMMOL(5)
      DIMENSION TC(NVAR)
      DIMENSION TL(NVAR)
      DIMENSION TR(NVAR)
      DIMENSION WORK1(NVAR)
      DIMENSION WORK2(NVAR)
      DIMENSION WORK3(NVAR)
      DIMENSION XC(NVAR)
      DIMENSION XF(NVAR)
      DIMENSION XR(NVAR)
      DIMENSION XEDIT(13)
      COMMON /COMM28/ CON(5),FLOW(5),FXOLD(13),PRESS,SCALE(13),
     & TEMPER

      LAB(1)='NO2 '
      LAB(2)='N2O4'
      LAB(3)='NO  '
      LAB(4)='H20 '
      LAB(5)='HNO3'
c
c  descale xr
c
      call p19_sc (1,NVAR,SCALE,XR)
c
c  compute unscaled quantities
c
      DO I=1,5
        XEDIT(I)=XR(I)*XR(11)
        XEDIT(I+5)=XR(I+5)*XR(12)
        SUMMOL(I)=XEDIT(I)+XEDIT(I+5)
      end do
c
c  get equilibrium constants
c
      call p19_cn (CON,NVAR,PRESS,TEMPER,XR)
c
c  write header
c
      WRITE(*,*)' '
      WRITE(*,*)' '
      WRITE(*,1020)XR(NVAR)
      WRITE(*,1000)PRESS
      WRITE(*,1010)TEMPER
      WRITE(*,*)' '
      WRITE(*,*)' '
      WRITE(*,1030)
      WRITE(*,*)' '
      DO I=1,5
        WRITE(*,1040)LAB(I),XEDIT(I),XEDIT(I+5),SUMMOL(I),
     &  FLOW(I),CON(I)
      end do
      XTOT=XR(11)+XR(12)
      FTOT=FLOW(1)+FLOW(2)+FLOW(3)+FLOW(4)+FLOW(5)
      WRITE(*,*)' '
      WRITE(*,1040)'Total',XR(11),XR(12),XTOT,FTOT
c
c  rescale xr
c
      call p19_sc (-1,NVAR,SCALE,XR)
      return
 1000 FORMAT(' PRESSURE IN ATMOSPHERES=',G14.6)
 1010 FORMAT(' TEMPERATURE IN KELVIN=',G14.6)
 1020 FORMAT(' EDIT OF SOLUTION FOR LAMBDA=',G14.6)
 1030 FORMAT(' REAGENT LIQUID MOLES  VAPOR MOLES  TOTAL MOLES   FLOWRA',
     & 'TE       EQUILIBRIUM')
 1040 FORMAT(1X,A5,5G14.6)
      END
      subroutine p19_fun ( option, nvar, x, fx )

c*********************************************************************72
c
cc P19_FUN evaluates the nitric absorption function.
c
      integer nvar

      double precision fx(nvar-1)
      integer option
      double precision x(nvar)

      COMMON /COMM28/ CON(5),FLOW(5),FXOLD(13),PRESS,SCALE(13),
     & TEMPER
c
c  descale xr
c
      call p19_sc (1,NVAR,SCALE,X)
c
c  get chemical constants
c
      call p19_cn (CON,NVAR,PRESS,TEMPER,X)
c
c  evaluate functions
c
      DO I=1,5
        FX(I)=X(11)*X(I)+X(12)*X(I+5)-FLOW(I)
      end do

      DO I=1,5
        FX(I+5)=X(I+5)-CON(I)*X(I)
      end do
      FX(11)=X(1)+X(2)+X(3)+X(4)+X(5)-1.0
      FX(12)=X(6)+X(7)+X(8)+X(9)+X(10)-1.0
c
c  adjust function for starting point
c
      NEQN=NVAR-1
      DO I=1,NEQN
        FX(I)=FX(I)+(X(NVAR)-1.0)*FXOLD(I)
      end do
c
c  rescale xr
c
      call p19_sc (-1,NVAR,SCALE,X)
      return
      END
      subroutine p19_cn (CON,NVAR,PRESS,TEMPER,X)

c*********************************************************************72
c
cc P19_CN evaluates constitutive relations for the nitric absorption function.
c
      DIMENSION X(NVAR)
      DIMENSION CON(5)
      CON(1)=1333.0/PRESS
      CON(2)=33.0/PRESS
      CON(3)=28780.0/PRESS
      ARG=11.99 - 4004.0/(TEMPER-39.06)
     & -8546.0*X(5)*X(5)/TEMPER +4.0*X(5)*X(5)
     & +6754.0*X(5)*X(5)*X(4)/TEMPER
     & -8.0*X(5)*X(5)*X(4)
     & -ALOG(PRESS)
      CON(4)=EXP(ARG)
      ARG=10.98 - 3362.0/(TEMPER-50.79)
     & -2872.0*X(4)*X(4)/TEMPER -6754.0*X(5)*X(4)*X(4)/TEMPER
     & +8.0*X(5)*X(4)*X(4)
     & -ALOG(PRESS)
      CON(5)=EXP(ARG)
      return
      END
      subroutine p19_jac ( option, nvar, x, jac )

c*********************************************************************72
c
cc P19_JAC evaluates the nitric absorption jacobian.
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
c    Input, integer OPTION, the option index.
c
c    Input, integer NVAR, the number of variables.
c
c    Input, double precision X(NVAR), the argument of the jacobian.
c
c    Output, double precision JAC(NVAR,NVAR), the jacobian matrix evaluated
c    at X.  The NVAR-th row is not set by this routine.
c
      implicit none

      integer nvar

      integer option
      double precision jac(nvar,nvar)
      double precision x(nvar)

      COMMON /COMM28/ CON(5),FLOW(5),FXOLD(13),PRESS,SCALE(13),
     & TEMPER
c
c  descale x
c
      call p19_sc (1,NVAR,SCALE,X)
c
c  zero out matrix
c
      do j = 1, ncol
        do i = 1, nrow
          jac(i,j) = 0.0E+00
        end do
      end do
c
c  get chemical constants
c
      call p19_cn (CON,NVAR,PRESS,TEMPER,X)
c
c  set d con(4)/ d x(4), etc.
c
      DC4DX4=CON(4)*(6754.0*X(5)*X(5)/TEMPER
     & -8.0*X(5)*X(5) )
      DC4DX5=CON(4)*(-8546.0*2.0*X(5)/TEMPER
     & +8.0*X(5) +6754.0*2.0*X(4)*X(5)/TEMPER
     & -16.0*X(5)*X(4) )
      DC5DX4=CON(5)*(-2872.0*2.0*X(4)/TEMPER
     & -6754.0*2.0*X(4)*X(5)/TEMPER + 16.0*X(4)*X(5) )
      DC5DX5=CON(5)*(-6754.0*X(4)*X(4)/TEMPER
     & +8.0*X(4)*X(4) )
c
c  evaluate partials
c
      DO I=1,5
        jac(I,11)=X(I)
        jac(I,I)=X(11)
        jac(I,12)=X(I+5)
        jac(I,I+5)=X(12)
      end do

      DO I=1,5
        jac(I+5,I+5)=1.0
        jac(I+5,I)=-CON(I)
        IF(I.EQ.4)jac(I+5,4)=jac(I+5,4)-DC4DX4*X(4)
        IF(I.EQ.4)jac(I+5,5)=jac(I+5,5)-DC4DX5*X(5)
        IF(I.EQ.5)jac(I+5,4)=jac(I+5,4)-DC5DX4*X(4)
        IF(I.EQ.5)jac(I+5,5)=jac(I+5,5)-DC5DX5*X(5)
      end do

      DO I=1,5
        jac(11,I)=1.0
        jac(12,I+5)=1.0
      end do
c
c  set last column
c
      NEQN=NVAR-1
      DO I=1,NEQN
        jac(I,NVAR)=FXOLD(I)
      end do
c
c  scale jacobian
c
      NEQN=NVAR-1
      DO I=1,NROW
        DO J=1,NEQN
          jac(I,J)=jac(I,J)*SCALE(J)
        end do
      end do
c
c  rescale x
c
      call p19_sc (-1,NVAR,SCALE,X)
      return
      END
      subroutine p19_sc (ISCALE,NVAR,SCALE,X)

c*********************************************************************72
c
cc P19_SC scales and unscales data for the nitric absorption problem.
c
      DIMENSION SCALE(NVAR)
      DIMENSION X(NVAR)
      IF(ISCALE.EQ.1)GO TO 10
      IF(ISCALE.EQ.-1)GO TO 20
      return
c
c  scale up
c
10    CONTINUE
      DO I=1,NVAR
        X(I)=X(I)*SCALE(I)
      end do
      return
c
c  scale down
c
20    CONTINUE
      DO I=1,NVAR
        X(I)=X(I)/SCALE(I)
      end do
      return
      END
      subroutine p19_nvar ( option, nvar )

c*********************************************************************72
c
cc P19_NVAR sets the number of variables for the nitric absorption problem.
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
c    Input, integer OPTION, the option chosen for this problem.
c    For some problems, several options are available.  At least,
c    OPTION = 1 is always legal.
c
c    Output, integer NVAR, the number of variables.
c
      implicit none

      integer nvar
      integer option

      nvar = 13

      return
      end
      subroutine p19_option_num ( option_num )

c*********************************************************************72
c
cc P19_OPTION_NUM returns the number of options for the nitric absorption problem.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    16 September 2008
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, integer OPTION_NUM, the number of options.
c
      implicit none

      integer option_num

      option_num = 2

      return
      end
      subroutine p19_start ( option, nvar, x )

c*********************************************************************72
c
cc P19_START returns a starting point for the nitric absorption problem.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    19 September 2008
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer OPTION, the option index.
c
c    Input, integer NVAR, the number of variables.
c
c    Output, double precision X(NVAR), the starting point.
c
      implicit none

      integer nvar

      integer i
      integer option
      double precision x(nvar)

      if ( option .eq. 1 ) then

        x(1)  =   0.00218216D+00
        x(2)  =   0.03171126D+00
        x(3)  =   0.00010562D+00
        x(4)  =   0.48301846D+00
        x(5)  =   0.48298250D+00
        x(6)  =   0.41554567D+00
        x(7)  =   0.14949595D+00
        x(8)  =   0.43425476D+00
        x(9)  =   0.00018983D+00
        x(10) =   0.00051379D+00
        x(11) = 207.02239583D+00
        x(12) =  22.97760417D+00
        x(13) =   1.00000000D+00

      else if ( option .eq. 2 ) then

        x(1) =   0.0
        x(2) =  10.0
        x(3) =   0.0
        x(4) = 100.0
        x(5) = 100.0
        x(6) =  10.0
        x(7) =   0.0
        x(8) =  10.0
        x(9) =   0.0
        x(10) =  0.0

        x(11) = 0.0
        do i = 1, 5
          x(11) = x(11) + x(i)
        end do

        do i = 1, 5
          x(i) = x(i) / x(11)
        end do

        x(12) = 0.0
        do i = 6, 10
          x(12) = x(12) + x(i)
        end do

        do i = 6, 10
          x(i) = x(i) / x(12)
        end do

        x(13) = 0.0

      else

        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'P19_START - Fatal error!'
        write ( *, '(a,i8)' ) '  Unrecognized option  = ', option
        stop

      end if



      return
      end
      subroutine p19_stepsize ( option, h, hmin, hmax )

c*********************************************************************72
c
cc P19_STEPSIZE returns step sizes for the nitric absorption problem.
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
c    Input, integer OPTION, the option index.
c
c    Output, double precision H, HMIN, HMAX, suggested values for the
c    initial step, the minimum step, and the maximum step.
c
      implicit none

      double precision h
      double precision hmax
      double precision hmin
      integer option

      h =    0.125000D+00
      hmin = 0.015625D+00
      hmax = 4.000000D+00

      return
      end
      subroutine p19_title ( option, title )

c*********************************************************************72
c
cc P19_TITLE sets the title for the nitric absorption problem.
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
c    Input, integer OPTION, the option index.
c
c    Output, character*(*) TITLE, the title of the problem.
c    TITLE will never be longer than 80 characters.
c
      implicit none

      integer option
      character*(*) title

      if ( option .eq. 1 ) then
        title = 'Nitric acid absorption function, good start.'
      else if ( option .eq. 2 ) then
        title = 'Nitric acid absorption function, bad start.'
      else
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'P19_TITLE - Fatal error!'
        write ( *, '(a)' ) '  Unrecognized option number.'
        stop
      end if

      return
      end
      subroutine p20_data (DIRIPC,HFACT,HMAX,HMIN,HTAN,option,
     & IPC,IPRAM,ISIZE,IT,LIM,MAXCON,
     & MAXLIM,MAXSTP,MAXTAR,NCOL,NROW,NVAR,RELERR,XIT,RWORK)

c*********************************************************************72
c
cc P20_DATA sets parameters for Rabinowitz's symmetric integration rule function.
c
c 20.1 rabinowitz symmetric integration formula.
c
c 20.2 reference
c
c  philip rabinowitz, nira richter
c  perfectly symmetric two dimensional integration formulas
c    with minimal numbers of points,
c  mathematics of computation, october 1969,
c  volume 23, number 108, pages 765-780.
c
c 20.3 the function
c
c  on the two-dimensional square between -1 and +1, a number
c  of symmetrically placed quadrature points and their weights
c  are sought for a given order quadrature rule.
c
c  The quadrature points fall into four classes.  each class has a generator
c  and if (x,y) is in a class, so are (+-x,+-y) and (+-y,+-x).
c
c  The following table lists the four classes
c
c  Generator  number of generators  weight  points per generator
c
c  (0,0)      k1 (=0 or 1)          a       1
c  (u,0)      k2                    b       4
c  (v,v)      k3                    c       4
c  (w,z)      k4                    d       8
c
!  By the choice of points, the integral of any odd function will
!  vanish.  Our conditions on the points will therefore involve the
!  approximation of even functions.  The equations fall into three groups.
c  note that i(f) denotes integral over the square of f(x,y)
c  and that for a different (symmetric) region the only change
c  required is the calculation of i(f).
c
c  i)  K1*A + 4*K2*B + 4*K3*C + 8*K4*D = I(1)
c
c  ii) 2*K2*B*u**2i + 4*K3*C*v**2i + 4*K4*d*(w**2i + z**2i) = i(x**2i)
c
c  iii) 4*K3*(c*v**2(i+j)) + 4*K4*d*(w**2j*z**2i + w**2i*z**2j)
c        = i(x**2i*y**2j)
c
c
c  data storage in the vector x is as follows:
c
c  (a,u1,b1,u2,b2,...,v1,c1,v2,c2,...,w1,z1,d1,w2,z2,d2,...)
c   k1   2*k2         2*k3            3*k4
c
c 20.4 options
c
c            k1 k2 k3 k4 nvar norder
c  option=1   0  2  3  1   13     11
c  option=2   1  4  4  3   26     17
c  option=3   0  5  2  4   26     17
c  option=4   0  4  3  4   26     17
c  option=5   0  3  4  4   26     17
c  option=6   0  2  5  4   26     17
c  option=7   1  4  1  5   26     17
c  option=8   1  3  2  5   26     17
c  option=9   1  2  3  5   26     17
c  option=10  1  1  4  5   26     17
c
c 20.5 target points
c
c 20.6 limit points
c
c 20.8 comments
c
c  option=1 is ok.  option=2 thru 10 do not have a starting point
c  and option=2 in particular has not gotten started.
c
      DIMENSION RWORK(ISIZE)
      COMMON /AUXMEM/ NORDER,K1,K2,K3,K4

      DIRIPC=1.0
      HFACT=3.0
      IPC=1
      IT=0
      LIM=0
      MAXCON=20
      MAXLIM=20
      MAXSTP=20
      MAXTAR=1
      NCOL=NVAR
      NROW=NVAR
      RELERR=1.0E-5
      XIT=0.0
      IF(option.EQ.1)GO TO 10
      IF(option.EQ.2)GO TO 20
      IF(option.EQ.3)GO TO 30
      IF(option.EQ.4)GO TO 40
      IF(option.EQ.5)GO TO 50
      IF(option.EQ.6)GO TO 60
      IF(option.EQ.7)GO TO 70
      IF(option.EQ.8)GO TO 80
      IF(option.EQ.9)GO TO 90
      IF(option.EQ.10)GO TO 100
c
c  first starting point
c
   10 NORDER=11
      K1=0
      K2=2
      K3=3
      K4=1

      GO TO 110
c
c  second starting point
c
   20 NORDER=17
      K1=1
      K2=4
      K3=4
      K4=3

      GO TO 110
c
c  third starting point
c
   30 NORDER=17
      K1=0
      K2=5
      K3=2
      K4=4
      GO TO 110
c
c  fourth starting point
c
   40 NORDER=17
      K1=0
      K2=4
      K3=3
      K4=4
      GO TO 110
c
c  fifth starting point
c
   50 NORDER=17
      K1=0
      K2=3
      K3=4
      K4=4
      GO TO 110
c
c  sixth starting point
c
   60 NORDER=17
      K1=0
      K2=2
      K3=5
      K4=4
      GO TO 110
c
c  seventh starting point
c
   70 NORDER=17
      K1=1
      K2=4
      K3=1
      K4=5
      GO TO 110
c
c  eighth starting point
c
   80 NORDER=17
      K1=1
      K2=3
      K3=2
      K4=5
      GO TO 110
c
c  ninth starting point
c
   90 NORDER=17
      K1=1
      K2=2
      K3=3
      K4=5
      GO TO 110
c
c  tenth starting point
c
  100 NORDER=17
      K1=1
      K2=1
      K3=4
      K4=5
      GO TO 110
  110 NVAR=K1+2*K2+2*K3+3*K4
      return
      END
      subroutine p20_fun ( option, nvar, x, fx )

c*********************************************************************72
c
cc P20_FUN evaluates the Rabinowitz symmetric integration formula function.
c
      integer nvar

      double precision fx(nvar-1)
      integer option
      double precision x(nvar)

      COMMON /AUXMEM/ NORDER,K1,K2,K3,K4
      NEQN=NVAR-1
c
c  The first equation, with right hand side Integral(1).
c
      TINTEG=4.0
      IEQN=1
      FX(IEQN)=-TINTEG
      ICOF=2

      IF(K1.eq.1) then
        FX(IEQN)=FX(IEQN)+X(1)
        ICOF=3
      end if

      DO I=1,K2
        FX(IEQN)=FX(IEQN)+4.0*X(ICOF)
        ICOF=ICOF+2
      end do

      DO I=1,K3
        FX(IEQN)=FX(IEQN)+4.0*X(ICOF)
        ICOF=ICOF+2
      end do

      ICOF=ICOF+1
      DO I=1,K4
        FX(IEQN)=FX(IEQN)+8.0*X(ICOF)
        ICOF=ICOF+3
      end do
c
c  Evaluate the second group of equations
c
      ITOP=NORDER/2
      IF(ITOP.LE.0)GO TO 150

      DO INUMBR = 1, ITOP

        IEQN=IEQN+1

        IEX=2*INUMBR
        TINTEG=4.0/(FLOAT(2*INUMBR+1))

        FX(IEQN)=-TINTEG
        ICOF=2
        IPT1=1

        IF(K1.eq.1)then
          ICOF=ICOF+1
          IPT1=IPT1+1
        end if

        DO I=1,K2
          FX(IEQN)=FX(IEQN)+2.0*X(ICOF)*X(IPT1)**IEX
          ICOF=ICOF+2
          IPT1=IPT1+2
        end do

        DO I=1,K3
          FX(IEQN)=FX(IEQN)+4.0*X(ICOF)*X(IPT1)**IEX
          ICOF=ICOF+2
          IPT1=IPT1+2
        end do

  120   IF(K4.LE.0)GO TO 140
        IPT2=IPT1+1
        ICOF=ICOF+1

        DO I=1,K4
          FX(IEQN)=FX(IEQN)+4.0*X(ICOF)*(X(IPT1)**IEX+X(IPT2)**IEX)
          IPT1=IPT1+3
          IPT2=IPT1+1
          ICOF=IPT1+2
        end do

  140   CONTINUE
c
c  Evaluate third group of equations
c
  150 continue

      ITOP=NORDER/4

      DO 210 I = 1, ITOP

        JLOW=I
        JTOP=(NORDER/2)-I

        DO 200 J = I, JTOP

          IEX=2*I
          JEX=2*J
          KEX=IEX+JEX
          IEQN=IEQN+1
          TINTEG=4.0/FLOAT((2*I+1)*(2*J+1))
          FX(IEQN)=-TINTEG
          IPT1=K1+2*K2+1
          ICOF=K1+2*K2+2
  160     continue

          DO K=1,K3
            FX(IEQN)=FX(IEQN)+4.0*X(ICOF)*X(IPT1)**KEX
            ICOF=ICOF+2
            IPT1=IPT1+2
          end do

  180     IF(K4.LE.0)GO TO 210
          IPT2=IPT1+1
          ICOF=ICOF+1

          DO K=1,K4
            TERM1=(X(IPT1)**IEX)*(X(IPT2)**JEX)
            TERM2=(X(IPT1)**JEX)*(X(IPT2)**IEX)
            FX(IEQN)=FX(IEQN)+4.0*X(ICOF)*(TERM1+TERM2)
            IPT1=IPT1+3
            IPT2=IPT1+1
            ICOF=IPT1+2
          end do

  200     CONTINUE

  210   CONTINUE

      end do

      return
      END
      subroutine p20_jac ( option, nvar, x, jac )

c*********************************************************************72
c
cc P20_JAC evaluates the Rabinowitz symmetric integration formula jacobian.
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
c    Input, integer OPTION, the option index.
c
c    Input, integer NVAR, the number of variables.
c
c    Input, double precision X(NVAR), the argument of the jacobian.
c
c    Output, double precision JAC(NVAR,NVAR), the jacobian matrix evaluated
c    at X.  The NVAR-th row is not set by this routine.
c
      implicit none

      integer nvar

      integer option
      double precision jac(nvar,nvar)
      double precision x(nvar)

      COMMON /AUXMEM/ NORDER,K1,K2,K3,K4

      do j = 1, ncol
        do i = 1, nrow
          jac(i,j) = 0.0E+00
        end do
      end do

      NEQN=NVAR-1
c
c  differentiate first equation
c
      IEQN=1
      ICOF=2

      IF ( K1.eq.1)then
        ICOF=1
        jac(IEQN,ICOF)=jac(IEQN,ICOF)+1.0
        ICOF=3
      end if

      DO I=1,K2
        jac(IEQN,ICOF)=jac(IEQN,ICOF)+4.0
        ICOF=ICOF+2
      end do

   50 continue

      DO I=1,K3
        jac(IEQN,ICOF)=jac(IEQN,ICOF)+4.0
        ICOF=ICOF+2
      end do

   70 IF(K4.LE.0)GO TO 90
      ICOF=ICOF+1

      DO I=1,K4
        jac(IEQN,ICOF)=jac(IEQN,ICOF)+8.0
        ICOF=ICOF+3
      end do
c
c  differentiate second group of equations
c
   90 ITOP=NORDER/2
      IF(ITOP.LE.0)GO TO 170
      ILOW=1
      DO 160 INUMBR=ILOW,ITOP
        IEQN=IEQN+1
        IEX=2*INUMBR
        IEXM1=IEX-1
        ICOF=2
        IPT1=1
        IF(K1.NE.1)GO TO 100
        ICOF=ICOF+1
        IPT1=IPT1+1

  100   continue

        DO I=1,K2
          jac(IEQN,ICOF)=jac(IEQN,ICOF)+2.0*X(IPT1)**IEX
          jac(IEQN,IPT1)=jac(IEQN,IPT1)
     &       +2.0*FLOAT(IEX)*X(ICOF)*X(IPT1)**IEXM1
          ICOF=ICOF+2
          IPT1=IPT1+2
        end do

  120   continue

        DO I=1,K3
          jac(IEQN,ICOF)=jac(IEQN,ICOF)+4.0*X(IPT1)**IEX
          jac(IEQN,IPT1)=jac(IEQN,IPT1)
     &      +4.0*FLOAT(IEX)*X(ICOF)*X(IPT1)**IEXM1
          ICOF=ICOF+2
          IPT1=IPT1+2
        end do

  140   IF(K4.LE.0)GO TO 160
        IPT2=IPT1+1
        ICOF=ICOF+1
        DO I=1,K4
          jac(IEQN,ICOF)=jac(IEQN,ICOF)
     &      +4.0*(X(IPT1)**IEX+X(IPT2)**IEX)
          jac(IEQN,IPT1)=jac(IEQN,IPT1)
     &      +4.0*FLOAT(IEX)*X(ICOF)*X(IPT1)**IEXM1
          jac(IEQN,IPT2)=jac(IEQN,IPT2)
     &      +4.0*FLOAT(IEX)*X(ICOF)*X(IPT2)**IEXM1
          IPT1=IPT1+3
          IPT2=IPT1+1
          ICOF=IPT1+2
        end do

  160   CONTINUE
c
c  differentiate third group of equations
c
  170 ILOW=1
      ITOP=NORDER/4

      DO 230 I=ILOW,ITOP

        JLOW=I
        JTOP=(NORDER/2)-I

        DO 220 J=JLOW,JTOP

          IEX=2*I
          IEXM1=IEX-1
          JEX=2*J
          JEXM1=JEX-1
          KEX=IEX+JEX
          KEXM1=KEX-1
          IEQN=IEQN+1
          IPT1=K1+2*K2+1
          ICOF=K1+2*K2+2
  180     IF(K3.LE.0)GO TO 200

          DO K=1,K3
            jac(IEQN,ICOF)=jac(IEQN,ICOF)+4.0*X(IPT1)**KEX
            jac(IEQN,IPT1)=jac(IEQN,IPT1)
     &        +4.0*FLOAT(KEX)*X(ICOF)*X(IPT1)**KEXM1
            ICOF=ICOF+2
            IPT1=IPT1+2
          end do

  200     IF(K4.LE.0)GO TO 230
          IPT2=IPT1+1
          ICOF=ICOF+1

          DO K=1,K4
            TERM1=(X(IPT1)**IEX)*(X(IPT2)**JEX)
            TERM2=(X(IPT1)**JEX)*(X(IPT2)**IEX)
            jac(IEQN,ICOF)=jac(IEQN,ICOF)+4.0*(TERM1+TERM2)
            jac(IEQN,IPT1)=jac(IEQN,IPT1)
     &        +4.0*X(ICOF)*FLOAT(IEX)*X(IPT1)**IEXM1*X(IPT2)**JEX
     &        +4.0*X(ICOF)*FLOAT(JEX)*X(IPT1)**JEXM1*X(IPT2)**IEX
            jac(IEQN,IPT2)=jac(IEQN,IPT2)
     &        +4.0*X(ICOF)*FLOAT(JEX)*X(IPT1)**IEX*X(IPT2)**JEXM1
     &        +4.0*X(ICOF)*FLOAT(IEX)*X(IPT1)**JEX*X(IPT2)**IEXM1
            IPT1=IPT1+3
            IPT2=IPT1+1
            ICOF=IPT1+2
          end do

  220     CONTINUE
  230   CONTINUE
  240 CONTINUE
      return
      END
      subroutine p20_nvar ( option, nvar )

c*********************************************************************72
c
cc P20_NVAR sets the number of variables for problem 20.
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
c    Input, integer OPTION, the option chosen for this problem.
c    For some problems, several options are available.  At least,
c    OPTION = 1 is always legal.
c
c    Output, integer NVAR, the number of variables.
c
      implicit none

      integer nvar
      integer option

      if ( option .eq. 1 ) then
        nvar = 13
      else
        nvar = 26
      end if

      return
      end
      subroutine p20_option_num ( option_num )

c*********************************************************************72
c
cc P20_OPTION_NUM returns the number of options for problem 20.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    16 September 2008
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, integer OPTION_NUM, the number of options.
c
      implicit none

      integer option_num

      option_num = 2

      return
      end
      subroutine p20_start ( option, nvar, x )

c*********************************************************************72
c
cc P20_START returns a starting point for problem 20.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    19 September 2008
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer OPTION, the option index.
c
c    Input, integer NVAR, the number of variables.
c
c    Output, double precision X(NVAR), the starting point.
c
      implicit none

      integer nvar

      double precision coef
      integer npts
      integer option
      double precision x(nvar)

      if ( option .eq. 1 ) then

        x(1) =  0.8989737240828844D+00
        x(2) =  0.0176679598882646D+00
        x(3) =  0.7632367891419969D+00
        x(4) =  0.2322248008989674D+00
        x(5) =  0.8949648832822285D+00
        x(6) =  0.0715516745178401D+00
        x(7) =  0.6322452037101431D+00
        x(8) =  0.2192868905662522D+00
        x(9) =  0.2797353125538562D+00
        x(10) = 0.2965842326220580D+00
        x(11) = 0.9602661668053869D+00
        x(12) = 0.4347413023856830D+00
        x(13) = 0.0813422207533089D+00

      else if ( option .eq. 2 ) then

        npts = 57
        coef = 4.0D+00 / dble ( npts )
        x(1) =  coef
        x(2) =  0.20D+00
        x(3) =  coef
        x(4) =  0.40D+00
        x(5) =  coef
        x(6) =  0.60D+00
        x(7) =  coef
        x(8) =  0.80D+00
        x(9) =  coef
        x(10) = 0.20D+00
        x(11) = coef
        x(12) = 0.40D+00
        x(13) = coef
        x(14) = 0.60D+00
        x(15) = coef
        x(16) = 0.80D+00
        x(17) = coef
        x(18) = 0.20D+00
        x(19) = 0.40D+00
        x(20) = coef
        x(21) = 0.30D+00
        x(22) = 0.50D+00
        x(23) = coef
        x(24) = 0.20D+00
        x(25) = 0.80D+00
        x(26) = coef

      else

        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'P20_START - Fatal error!'
        write ( *, '(a,i8)' ) '  Unrecognized option  = ', option
        stop

      end if

      return
      end
      subroutine p20_stepsize ( option, h, hmin, hmax )

c*********************************************************************72
c
cc P20_STEPSIZE returns step sizes for problem 20.
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
c    Input, integer OPTION, the option index.
c
c    Output, double precision H, HMIN, HMAX, suggested values for the
c    initial step, the minimum step, and the maximum step.
c
      implicit none

      double precision h
      double precision hmax
      double precision hmin
      integer option

      h =    0.250D+00
      hmin = 0.001D+00
      hmax = 1.000D+00

      return
      end
      subroutine p20_title ( option, title )

c*********************************************************************72
c
cc P20_TITLE sets the title for problem 20.
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
c    Input, integer OPTION, the option index.
c
c    Output, character*(*) TITLE, the title of the problem.
c    TITLE will never be longer than 80 characters.
c
      implicit none

      integer option
      character*(*) title

      if ( option .eq. 1 ) then
        title =
     &  'The Rabinowitz integration formula (0,2,3,1), 13 variables.'
      else if ( option .eq. 2 ) then
        title =
     &  'The Rabinowitz integration formula (1,4,4,3), 26 variables.'
      else
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'P20_TITLE - Fatal error!'
        write ( *, '(a)' ) '  Unrecognized option number.'
        stop
      end if

      return
      end
      subroutine p21_data (DIRIPC,HFACT,HMAX,HMIN,HTAN,option,
     & IPC,IPRAM,ISIZE,IT,LIM,MAXCON,
     & MAXLIM,MAXSTP,MAXTAR,NCOL,NROW,NVAR,RELERR,XIT,RWORK)

c*********************************************************************72
c
cc P21_DATA sets parameters for Porsching's tube flow.
c
c 21.1 porsching tube flow
c
c 21.2 reference
c
c  t. a. porsching,
c  work in progress.
c
c 21.3 the function
c
c  a network of links and nodes is given schematically as
c
c  - x - x - x - x
c
c  where - represents a link, and x a node.
c
c  quantities associated with a link are
c
c  alfa   - the void fraction on the link.
c  vliq   - the velocity of the liquid on the link.
c  vgas   - the velocity of the gas on te link.
c
c  quantities associated with a node are
c
c  pres   - the pressure.
c
c  at the first (extreme left) link, alfa, vliq and vgas are
c  specified.  at the last (extreme right) node, pres is specified.
c
c  the equations for the system then have the following form,
c
c  specifications at left link -
c
c  fx1 = alfa(link1) - alfspc
c  fx2 = vliq(link1) - vlqspc
c  fx3 = vgas(link1) - vgsspc
c
c  equations for node with linkl and linkr as left and right link
c  neighbors.
c
c  fx(k)  = alfa(linkr)*vgas(linkr)-alfa(linkl)*vgas(linkl)
c         - dell*gamma(node)*scon/rhogas
c
c  fx(k+1)= (1-alfa(linkr))*vliq(linkr) - (1-alfa(linkl))*vliq(linkl)
c         + dell*gamma(node)*scon/rholiq
c
c  equations for link with nodes nodel and noder as left and right
c  node neighbors
c
c  fx(k+2)= hgas + scon**2*hmix
c           - alfa(link)*(pres(nodel)-pres(noder))/dell
c           - rhogas*alfa(link)*gravty(link)*scon
c  fx(k+3)= hliq - scon**2*hmix
c           - (1-alfa(link))*(pres(nodel)-pres(noder))/dell
c           - rholiq*(1-alfa(link))*gravty(link)*scon
c
c  specification of pressure at last node
c
c  fx(last) = pres(node) - prsspc
c
c  where the following dictionary applies
c
c  dell   - length of a link (uniform value).
c  gamma  - gas generation rate at node.
c  gravty - gravitational acceleration on link.
c  hgas   - a gas friction factor.
c  hliq   - liquid friction factor.
c  hmix   - mixture friction factor.
c  rhogas - density of gas.
c  rholiq - density of liquid.
c  scon   - value of continuation parameter which initially
c           blocks out hmix term.
c
c 21.4 options
c
c  option=1  nodes=2
c  option=2  nodes=4
c  option=3  nodes=8
c  option=4  nodes=16
c
c 21.5 target point
c
c  the target point sought has s=1.0
c  at which point the full system of equations has been
c  solved.
c
c 21.8 comments
c
      real pi
      parameter ( pi = 3.1415926535E+00 )
      DIMENSION RWORK(ISIZE)
      COMMON /AUXMEM/ DELL,RHOLIQ,RHOGAS,RHOMIX,ANGLE,FLIQ,FGAS,FMIX,
     & WCRIT,SIGMA,NODES,GAMMA(16),GRAVTY(16),
     & ALFSPC,VLQSPC,VGSSPC,PRSSPC

      IF(option.EQ.1)NODES=2
      IF(option.EQ.2)NODES=4
      IF(option.EQ.3)NODES=8
      IF(option.EQ.4)NODES=16
      NVAR=4*NODES+1

      DIRIPC=1.0
      HFACT=3.0
      IPC=NVAR
      IT=NVAR
      LIM=0
      MAXCON=15
      MAXLIM=15
      MAXSTP=15
      MAXTAR=1
      NCOL=NVAR
      NROW=NVAR
      RELERR=1.0E-3
      XIT=1.0
      DELL=1.0
      RHOGAS=0.04
      RHOLIQ=60.0
      ANGLE=PI
      FLIQ=.021
      FGAS=0.06
      FMIX=3.43
      RHOMIX=60.0
      WCRIT=1.24
      SIGMA=0.161
      DO I=1,NODES
        GAMMA(I)=0.0
        GRAVTY(I)=0.0
      end do
c
c  boundary specifications of alpha, vliq, vgas and pres
c
      ALFSPC=.25
      VLQSPC=1.0
      VGSSPC=10.0
      PRSSPC=0.0
c
c  computations for starting point
c
      call p21_cn (ALFSPC,VLQSPC,VGSSPC,HLIQ,HGAS,HMIX)
      ALIQ=SQRT((1.0-ALFSPC)*HLIQ)
      AGAS=SQRT(ALFSPC*HGAS)
      ALFVAL=AGAS/(AGAS+ALIQ)
      VLQVAL=(1.0-ALFSPC)*VLQSPC*(AGAS+ALIQ)/ALIQ
      VGSVAL=ALFSPC*VGSSPC*(AGAS+ALIQ)/AGAS
      DELPRS=DELL*(AGAS+ALIQ)**2
c
c  starting point
c
      DO I=1,NODES
        J=4*(I-1)+1
        RWORK(J)=ALFVAL
        RWORK(J+1)=VLQVAL
        RWORK(J+2)=VGSVAL
        RWORK(J+3)=PRSSPC+(NODES-I)*DELPRS
      end do

      RWORK(1)=ALFSPC
      RWORK(2)=VLQSPC
      RWORK(3)=VGSSPC
      RWORK(NVAR-1)=PRSSPC
      RWORK(NVAR)=0.0
      return
      END
      subroutine p21_fun ( option, nvar, x, fx )

c*********************************************************************72
c
cc P21_FUN evaluates the Porsching tube flow function.
c
      integer nvar

      double precision fx(nvar-1)
      integer option
      double precision x(nvar)

      COMMON /AUXMEM/ DELL,RHOLIQ,RHOGAS,RHOMIX,ANGLE,FLIQ,FGAS,FMIX,
     & WCRIT,SIGMA,NODES,GAMMA(16),GRAVTY(16),
     & ALFSPC,VLQSPC,VGSSPC,PRSSPC
      NEQN=NVAR-1
      SCON=X(NVAR)
      LINK=1
      NODE=1
c
c  equations for first link
c
      FX(1)=X(1)-ALFSPC
      FX(2)=X(2)-VLQSPC
      FX(3)=X(3)-VGSSPC
      IEQN=3
c
c  set pointers for left and right variables
c
      IALF=1
      IVLQ=2
      IVGS=3
      IPRS=4
      JALF=5
      JVLQ=6
      JVGS=7
      JPRS=8
c
c  equations for next node
c
   10 IEQN=IEQN+1
      FX(IEQN)=X(JALF)*X(JVGS)-X(IALF)*X(IVGS)
     & -DELL*GAMMA(NODE)*SCON/RHOGAS
      IEQN=IEQN+1
      FX(IEQN)=(1.0-X(JALF))*X(JVLQ)-(1.0-X(IALF))*X(IVLQ)
     & +DELL*GAMMA(NODE)*SCON/RHOLIQ
c
c  equations for next link
c
      LINK=LINK+1
      call p21_cn (X(JALF),X(JVLQ),X(JVGS),HLIQ,HGAS,HMIX)
      IEQN=IEQN+1
      FX(IEQN)=HGAS+SCON**2*HMIX-X(JALF)*(X(IPRS)-X(JPRS))/DELL
     & -RHOGAS*X(JALF)*GRAVTY(LINK)*SCON
      IEQN=IEQN+1
      FX(IEQN)=HLIQ-SCON**2*HMIX-(1.0-X(JALF))*(X(IPRS)-X(JPRS))/DELL
     & -RHOLIQ*(1.0-X(JALF))*GRAVTY(LINK)*SCON
c
c  check whether to finish up or repeat for new pair
c
      NODE=NODE+1
      IF(NODE.GE.NODES)GO TO 20
c
c  advance pointers for next link, node pair
c
      IALF=IALF+4
      IVLQ=IVLQ+4
      IVGS=IVGS+4
      IPRS=IPRS+4
      JALF=JALF+4
      JVLQ=JVLQ+4
      JVGS=JVGS+4
      JPRS=JPRS+4
      GO TO 10
c
c  equation for last node
c
   20 IEQN=IEQN+1
      FX(IEQN)=X(JPRS)-PRSSPC
      return
      END
      subroutine p21_jac ( option, nvar, x, jac )

c*********************************************************************72
c
cc P21_JAC evaluates the Porsching tube flow jacobian.
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
c    Input, integer OPTION, the option index.
c
c    Input, integer NVAR, the number of variables.
c
c    Input, double precision X(NVAR), the argument of the jacobian.
c
c    Output, double precision JAC(NVAR,NVAR), the jacobian matrix evaluated
c    at X.  The NVAR-th row is not set by this routine.
c
      implicit none

      integer nvar

      integer option
      double precision jac(nvar,nvar)
      double precision x(nvar)

      COMMON /AUXMEM/ DELL,RHOLIQ,RHOGAS,RHOMIX,ANGLE,FLIQ,FGAS,FMIX,
     & WCRIT,SIGMA,NODES,GAMMA(16),GRAVTY(16),
     & ALFSPC,VLQSPC,VGSSPC,PRSSPC

      do j = 1, ncol
        do i = 1, nrow
          jac(i,j) = 0.0E+00
        end do
      end do

      NEQN=NVAR-1
      SCON=X(NVAR)
c
c  set pointers for left and right variables
c
      IALF=1
      IVLQ=2
      IVGS=3
      IPRS=4
      JALF=5
      JVLQ=6
      JVGS=7
      JPRS=8
      LINK=1
      NODE=1
c
c  equations for first link
c
      jac(1,IALF)=1.0
      jac(2,IVLQ)=1.0
      jac(3,IVGS)=1.0
      IEQN=3
c
c  equations for next node
c
   30 IEQN=IEQN+1
      jac(IEQN,JALF)=X(JVGS)
      jac(IEQN,IALF)=-X(IVGS)
      jac(IEQN,JVGS)=X(JALF)
      jac(IEQN,IVGS)=-X(IALF)
      jac(IEQN,NVAR)=-DELL*GAMMA(NODE)/RHOGAS
      IEQN=IEQN+1
      jac(IEQN,JALF)=-X(JVLQ)
      jac(IEQN,JVLQ)=1.0-X(JALF)
      jac(IEQN,IALF)=X(IVLQ)
      jac(IEQN,IVLQ)=X(IALF)-1.0
      jac(IEQN,NVAR)=DELL*GAMMA(NODE)/RHOLIQ
c
c  equations for next link
c
      LINK=LINK+1
      call p21_cn (X(JALF),X(JVLQ),X(JVGS),HLIQ,HGAS,HMIX)
      call p21_cp (X(JALF),X(JVLQ),X(JVGS),HLIQAL,HLIQVL,HGASAL,HGASVG,
     & HMIXAL,HMIXVL,HMIXVG)
      IEQN=IEQN+1
      jac(IEQN,JALF)=HGASAL+SCON**2*HMIXAL-(X(IPRS)-X(JPRS))/DELL
     & -RHOGAS*GRAVTY(LINK)*SCON
      jac(IEQN,JVLQ)=SCON**2*HMIXVL
      jac(IEQN,JVGS)=HGASVG+SCON**2*HMIXVG
      jac(IEQN,IPRS)=-X(JALF)/DELL
      jac(IEQN,JPRS)=X(JALF)/DELL
      jac(IEQN,NVAR)=2.0*SCON*HMIX-RHOGAS*X(JALF)*GRAVTY(LINK)
      IEQN=IEQN+1
      jac(IEQN,JALF)=HLIQAL-SCON**2*HMIXAL+(X(IPRS)-X(JPRS))/DELL
     & +RHOLIQ*GRAVTY(LINK)*SCON
      jac(IEQN,JVLQ)=HLIQVL-SCON**2*HMIXVL
      jac(IEQN,JVGS)=-SCON**2*HMIXVG
      jac(IEQN,IPRS)=(X(JALF)-1.0)/DELL
      jac(IEQN,JPRS)=(1.0-X(JALF))/DELL
      jac(IEQN,NVAR)=-2.0*SCON*HMIX
     & -RHOLIQ*(1.0-X(JALF))*GRAVTY(LINK)
c
c  check whether to finish up or repeat for new pair
c
      NODE=NODE+1
      IF(NODE.GE.NODES)GO TO 40
c
c  advance pointers for next link, node pair
c
      IALF=IALF+4
      IVLQ=IVLQ+4
      IVGS=IVGS+4
      IPRS=IPRS+4
      JALF=JALF+4
      JVLQ=JVLQ+4
      JVGS=JVGS+4
      JPRS=JPRS+4
      GO TO 30
c
c  equation for last node
c
   40 IEQN=IEQN+1
      jac(IEQN,JPRS)=1.0
      return
      END
      subroutine p21_cn(ALPHA,VLIQ,VGAS,HLIQ,HGAS,HMIX)

c*********************************************************************72
c
cc P21_CN evaluates constitutive relations for the Porsching tube flow function.
c
      real pi
      parameter ( pi = 3.1415926535E+00 )
      COMMON /AUXMEM/ DELL,RHOLIQ,RHOGAS,RHOMIX,ANGLE,FLIQ,FGAS,FMIX,
     & WCRIT,SIGMA,NODES,GAMMA(16),GRAVTY(16),
     & ALFSPC,VLQSPC,VGSSPC,PRSSPC

      HLIQ=SQRT(PI)*RHOLIQ*FLIQ*(1.0-ALPHA)*ABS(VLIQ)*VLIQ
      HLIQ=HLIQ/(4.0*SQRT(ANGLE))
      HGAS=SQRT(PI)*RHOGAS*FGAS*ALPHA*ABS(VGAS)*VGAS
      HGAS=HGAS/(4.0*SQRT(ANGLE))
      HMIX=3.0*RHOMIX*RHOGAS*FMIX*ALPHA*ABS(VGAS-VLIQ)*(VGAS-VLIQ)**3
      HMIX=HMIX/(4.0*0.06147*WCRIT*SIGMA)

      return
      END
      subroutine p21_cp (ALPHA,VLIQ,VGAS,HLIQAL,HLIQVL,HGASAL,HGASVG,
     & HMIXAL,HMIXVL,HMIXVG)

c*********************************************************************72
c
cc P21_CP differentiates constitutive relations for the Porsching tube flow function.
c
      real pi
      parameter ( pi = 3.1415926535E+00 )
      COMMON /AUXMEM/ DELL,RHOLIQ,RHOGAS,RHOMIX,ANGLE,FLIQ,FGAS,FMIX,
     & WCRIT,SIGMA,NODES,GAMMA(16),GRAVTY(16),
     & ALFSPC,VLQSPC,VGSSPC,PRSSPC

      HLIQAL=-SQRT(PI)*RHOLIQ*FLIQ*ABS(VLIQ)*VLIQ
      HLIQAL=HLIQAL/(4.0*SQRT(ANGLE))
      HLIQVL=2.0*SQRT(PI)*RHOLIQ*FLIQ*(1.0-ALPHA)*SIGN(1.0,VLIQ)*VLIQ
      HLIQVL=HLIQVL/(4.0*SQRT(ANGLE))
      HGASAL=SQRT(PI)*RHOGAS*FGAS*ABS(VGAS)*VGAS
      HGASAL=HGASAL/(4.0*SQRT(ANGLE))
      HGASVG=2.0*SQRT(PI)*RHOGAS*FGAS*ALPHA*SIGN(1.0,VGAS)*VGAS
      HGASVG=HGASVG/(4.0*SQRT(ANGLE))
      HMIXAL=3.0*RHOMIX*RHOGAS*FMIX*ABS(VGAS-VLIQ)*(VGAS-VLIQ)**3
      HMIXAL=HMIXAL/(4.0*0.06147*WCRIT*SIGMA)
      DIF=VGAS-VLIQ
      HMIXVL=-12.0*RHOMIX*RHOGAS*FMIX*ALPHA*SIGN(1.0,DIF)*DIF**3
      HMIXVL=HMIXVL/(4.0*0.06147*WCRIT*SIGMA)
      HMIXVG=12.0*RHOMIX*RHOGAS*FMIX*ALPHA*SIGN(1.0,DIF)*DIF**3
      HMIXVG=HMIXVG/(4.0*0.06147*WCRIT*SIGMA)
      return
      END
      subroutine p21_nvar ( option, nvar )

c*********************************************************************72
c
cc P21_NVAR sets the number of variables for problem 21.
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
c    Input, integer OPTION, the option chosen for this problem.
c    For some problems, several options are available.  At least,
c    OPTION = 1 is always legal.
c
c    Output, integer NVAR, the number of variables.
c
      implicit none

      integer nodes
      integer nvar
      integer option

      if ( option .eq. 1 ) then
        nodes = 2
      else if ( option .eq. 2 ) then
        nodes = 4
      else if ( option .eq. 3 ) then
        nodes = 8
      else if ( option .eq. 4 ) then
        nodes = 16
      end if

      nvar = 4 * nodes + 1

      return
      end
      subroutine p21_option_num ( option_num )

c*********************************************************************72
c
cc P21_OPTION_NUM returns the number of options for problem 21.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    16 September 2008
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, integer OPTION_NUM, the number of options.
c
      implicit none

      integer option_num

      option_num = 4

      return
      end
      subroutine p21_start ( option, nvar, x )

c*********************************************************************72
c
cc P21_START returns a starting point for problem 21.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    19 September 2008
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer OPTION, the option index.
c
c    Input, integer NVAR, the number of variables.
c
c    Output, double precision X(NVAR), the starting point.
c
      implicit none

      integer nvar

      integer option
      double precision x(nvar)

      if ( ? ) then

      else

        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'P21_START - Fatal error!'
        write ( *, '(a,i8)' ) '  Unrecognized option  = ', option
        stop

      end if

      return
      end
      subroutine p21_stepsize ( option, h, hmin, hmax )

c*********************************************************************72
c
cc P21_STEPSIZE returns step sizes for problem 21.
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
c    Input, integer OPTION, the option index.
c
c    Output, double precision H, HMIN, HMAX, suggested values for the
c    initial step, the minimum step, and the maximum step.
c
      implicit none

      double precision h
      double precision hmax
      double precision hmin
      integer option

      h =    0.250D+00
      hmin = 0.001D+00
      hmax = 0.500D+00

      return
      end
      subroutine p21_title ( option, title )

c*********************************************************************72
c
cc P21_TITLE sets the title for problem 21.
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
c    Input, integer OPTION, the option index.
c
c    Output, character*(*) TITLE, the title of the problem.
c    TITLE will never be longer than 80 characters.
c
      implicit none

      integer option
      character*(*) title

      if ( option .eq. 1 ) then
        title = 'Porsching''s tube flow function, 2 nodes.'
      else if ( option .eq. 2 ) then
        title = 'Porsching''s tube flow function, 4 nodes.'
      else if ( option .eq. 3 ) then
        title = 'Porsching''s tube flow function, 8 nodes.'
      else if ( option .eq. 4 ) then
        title = 'Porsching''s tube flow function, 16 nodes.'
      else
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'P21_TITLE - Fatal error!'
        write ( *, '(a)' ) '  Unrecognized option number.'
        stop
      end if

      return
      end
      subroutine p22_data (DIRIPC,HFACT,HMAX,HMIN,HTAN,option,
     & IPC,IPRAM,ISIZE,IT,LIM,MAXCON,
     & MAXLIM,MAXSTP,MAXTAR,NCOL,NROW,NVAR,RELERR,XIT,RWORK)

c*********************************************************************72
c
cc P22_DATA sets parameters for the arch limit run on lambda.
c
c 22.1 arch limit run on lambda
c
c 22.2 reference
c
c  werner c rheinboldt
c  work in progress
c
c 22.3 the function
c
c  an arch consists of a segment of a circle.  the
c  radius of the circle is 10, and the angle subtended
c  is 30 degrees.  a sixth degree polynomial is used
c  to approximate the function u1=radial displacement
c  and u2=tangential displacement of the arch, when
c  subject to a load.  the coefficients of the polynomial
c  are the actual variables.  there are 7 coefficients
c  for u1, and 6 more used for u2.
c  to satisfy the boundary conditions at -1 and +1,
c  which are u1=u2=0, only quadratic and higher terms
c  of the integrated legendre polynomials are used.
c  the load function has three parameters,
c  lambda, mu, and nu, and has the form
c
c  load(x) = lambda*(bload+exp(-mu*(x-nu)**2))/s
c
c  where bload is a base load
c  and s is the integral of bload+exp(-mu*(x-nu)**2))
c  from -1 to 1.
c
c  during continuation, mu and nu are fixed for this problem.
c
c 22.4 options
c
c  option=1  vertical load, mu=0.0, nu=0.0
c  option=2  radial load,   mu=0.0, nu=0.0
c  option=3  vertical load, mu=0.1, nu=.2304583159
c  option=4  radial load,   mu=0.1, nu=.2304583159
c  option=5  vertical load, mu=0.1, nu=.6423493394
c  option=6  radial load,   mu=0.1, nu=.6423493394
c  option=7  vertical load, mu=0.1, nu=.9841830547
c  option=8  radial load,   mu=0.1, nu=.9841830547
c  option=9  vertical load, mu=0.1, nu=0.0
c  option=10 radial load,   mu=0.1, nu=0.0
c
c 22.5 target points
c
c  no target points are sought
c
c 22.6 limit points
c
c  option=1,  lambda=4.34374
c  option=2,  lambda=4.42049
c  option=3,  lambda=1.76144
c  option=5,  lambda=1.63558
c  option=7,  lambda=1.56086
c
c 22.7  bifurcation points
c
c  for mu=0.0, a bifurcation point is passed on
c  the way to the limit point.
c
c 22.8 comments
c
c  for mu=0, the lambda values here represent a uniform load.
c  divided by 2, the value may be compared to the values computed
c  for the walker's arch problem, pitlib number 17.  in particular,
c  the turning point in lambda seems to match.
c
      real pi
      parameter ( pi = 3.1415926535E+00 )
      DIMENSION RWORK(ISIZE)
      COMMON /AUXMEM/ THETA,COEF2,COEF3,COEF4,COEF5,IFUN22,
     & NGAUSS,GAUSCO(13),GAUSAB(13),PL(13),PLD(13),PLDD(13),
     & UX(2),UD(2),UDD(2),NCOF,ILM,IMU,INU,VALLM,VALMU,VALNU,BLOAD,
     & CLOAD,THICK,RADIUS,ALFA0,ALFA1,ALFA2

      NCOF=13
      ILM=14
      IMU=15
      INU=16
      NVAR=NCOF+3

      DIRIPC=1.0
      HFACT=2.0
      IPC=ILM
      IT=0
      LIM=ILM
      MAXCON=30
      MAXLIM=1
      MAXSTP=30
      MAXTAR=30
      NCOL=NVAR
      NROW=NVAR
      RELERR=1.0E-3
      XIT=0.0
      THETA=PI/12.0
      EI=796.0
      EA=2.056*10.0**6
      THICK=.0684
      RADIUS=10.0
      ALFA0=THICK/RADIUS
      ALFA1=EI/(EA*RADIUS*RADIUS)
      ALFA2=EI/(EA*THICK*THICK)
      WRITE(*,*)' '
      WRITE(*,*)' '
      WRITE(*,1050)ALFA0,ALFA1,ALFA2
c
c  ifun22=1 for constant load
c  ifun22=2 for linear load (hat)
c  ifun22=3 for exponential load
c
      IFUN22=3
      WRITE(*,*)' '
      IF(IFUN22.EQ.1)WRITE(*,1060)
      IF(IFUN22.EQ.2)WRITE(*,1070)
      IF(IFUN22.EQ.3)WRITE(*,1080)
      COEF2=2.0/THETA
      COEF3=ALFA0/THETA
      COEF4=2.0*ALFA1/(THETA**3)
      COEF5=2.0*ALFA2*THETA
      WRITE(*,1010)COEF2
      WRITE(*,1020)COEF3
      WRITE(*,1030)COEF4
      WRITE(*,1040)COEF5
      VALLM=0.0
      IF(option.EQ.1.OR.option.EQ.2)VALMU=0.0
      IF(option.GE.3.AND.option.LE.10)VALMU=0.1
      IF(IFUN22.EQ.2)VALMU=1.0
      IF(option.EQ.1.OR.option.EQ.2)VALNU=0.0
      IF(option.EQ.3.OR.option.EQ.4)VALNU=0.2304583159
      IF(option.EQ.5.OR.option.EQ.6)VALNU=0.6423493394
      IF(option.EQ.7.OR.option.EQ.8)VALNU=0.9841830547
      IF(option.EQ.9.OR.option.EQ.10)VALNU=0.0
      BLOAD=0.0
      IF(option.EQ.1.OR.option.EQ.3.OR.option.EQ.5)CLOAD=1.0
      IF(option.EQ.7.OR.option.EQ.9)CLOAD=1.0
      IF(option.EQ.2.OR.option.EQ.4.OR.option.EQ.6)CLOAD=0.0
      IF(option.EQ.8.OR.option.EQ.10)CLOAD=0.0
      DO I=1,NCOF
        RWORK(I)=0.0
      end do
      RWORK(ILM)=VALLM
      RWORK(IMU)=VALMU
      RWORK(INU)=VALNU
c
c  store gauss coefficients and abscissas for 13-th order
c  quadrature on (-1,1)
c
      NGAUSS=13
      GAUSCO(1)=0.4048400476E-01
      GAUSCO(2)=0.9212149983E-01
      GAUSCO(3)=0.1388735102
      GAUSCO(4)=0.1781459807
      GAUSCO(5)=0.2078160475
      GAUSCO(6)=0.2262831802
      GAUSCO(7)=0.2325515532
      GAUSCO(8)=GAUSCO(6)
      GAUSCO(9)=GAUSCO(5)
      GAUSCO(10)=GAUSCO(4)
      GAUSCO(11)=GAUSCO(3)
      GAUSCO(12)=GAUSCO(2)
      GAUSCO(13)=GAUSCO(1)
      GAUSAB(1)=-0.9841830547
      GAUSAB(2)=-0.9175983992
      GAUSAB(3)=-0.8015780907
      GAUSAB(4)=-0.6423493394
      GAUSAB(5)=-0.4484927510
      GAUSAB(6)=-0.2304583159
      GAUSAB(7)= 0.0000000000
      GAUSAB(8)=-GAUSAB(6)
      GAUSAB(9)=-GAUSAB(5)
      GAUSAB(10)=-GAUSAB(4)
      GAUSAB(11)=-GAUSAB(3)
      GAUSAB(12)=-GAUSAB(2)
      GAUSAB(13)=-GAUSAB(1)
      return
 1010 FORMAT(1X,14HCOEF2=2/THETA=,G14.6)
 1020 FORMAT(1X,18HCOEF3=ALFA0/THETA=,G14.6)
 1030 FORMAT(1X,25HCOEF4=2*ALFA1/(THETA**3)=,G14.6)
 1040 FORMAT(1X,20HCOEF5=2*ALFA2*THETA=,G14.6)
 1050 FORMAT(7H ALFA0=,G14.6,7H ALFA1=,G14.6,7H ALFA2=,G14.6)
 1060 FORMAT(28H CONSTANT LOAD FUNCTION USED)
 1070 FORMAT(26H LINEAR LOAD FUNCTION USED)
 1080 FORMAT(31H EXPONENTIAL LOAD FUNCTION USED)
      END
      subroutine p22_fun ( option, nvar, x, fx )

c*********************************************************************72
c
cc P22_FUN evaluates the arch function for limit run on lambda.
c
      integer nvar

      double precision fx(nvar-1)
      integer option
      double precision x(nvar)

      DIMENSION Y(13)
      COMMON /AUXMEM/ THETA,COEF2,COEF3,COEF4,COEF5,IFUN22,
     & NGAUSS,GAUSCO(13),GAUSAB(13),PL(13),PLD(13),PLDD(13),
     & UX(2),UD(2),UDD(2),NCOF,ILM,IMU,INU,VALLM,VALMU,VALNU,BLOAD,
     & CLOAD,THICK,RADIUS,ALFA0,ALFA1,ALFA2
c
c  get terms involving coefficients
c
      call p22_gx (NVAR,X,FX)
c
c  get load terms
c
      call p22_ld (X,1,Y)
c
c  subtract load terms for first ncof functions
c
      DO IEQN=1,NCOF
        FX(IEQN)=FX(IEQN)-Y(IEQN)
      end do
c
c  equation fixing value of mu
c
      IEQN=NCOF+1
      FX(IEQN)=X(IMU)-VALMU
c
c  equation fixing value of nu
c
      IEQN=NCOF+2
      FX(IEQN)=X(INU)-VALNU
      return
      END
      subroutine p22_jac ( option, nvar, x, jac )

c*********************************************************************72
c
cc P22_JAC evaluates the arch limit run on lambda jacobian.
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
c    Input, integer OPTION, the option index.
c
c    Input, integer NVAR, the number of variables.
c
c    Input, double precision X(NVAR), the argument of the jacobian.
c
c    Output, double precision JAC(NVAR,NVAR), the jacobian matrix evaluated
c    at X.  The NVAR-th row is not set by this routine.
c
      implicit none

      integer nvar

      integer option
      double precision jac(nvar,nvar)
      double precision x(nvar)
      double precision y(13)

      common /auxmem/ theta,coef2,coef3,coef4,coef5,ifun22,
     & ngauss,gausco(13),gausab(13),pl(13),pld(13),pldd(13),
     & ux(2),ud(2),udd(2),ncof,ilm,imu,inu,vallm,valmu,valnu,bload,
     & cload,thick,radius,alfa0,alfa1,alfa2
c
c  zero out jacobian
c
      do j = 1, ncol
        do i = 1, nrow
          jac(i,j) = 0.0e+00
        end do
      end do
c
c  compute d fx/d coefficient ivar
c  and store
c
      do i=1,ncof
        ivar=i
        call p22_gp (x,ivar,y)
        do ieqn=1,ncof
          jac(ieqn,ivar)=y(ieqn)
        end do
      end do
c
c  compute d load/d lambda,
c          d load/d mu,
c      and d load/d nu
c  and store
c
      call p22_ld (x,2,y)
      ivar=ilm

      do ieqn=1,ncof
        jac(ieqn,ivar)=-y(ieqn)
      end do

      call p22_ld (x,3,y)
      ivar=imu

      do ieqn=1,ncof
        jac(ieqn,ivar)=-y(ieqn)
      end do

      call p22_ld (x,5,y)
      ivar=inu
      do ieqn=1,ncof
        jac(ieqn,ivar)=-y(ieqn)
      end do
c
c  equation fixing mu
c
      ieqn=14
      ivar=imu
      jac(ieqn,ivar)=1.0
c
c  equation fixing nu
c
      ieqn=15
      ivar=inu
      jac(ieqn,ivar)=1.0
      return
      end
      subroutine p22_edit (abserr,arclxc,arclxf,arclxr,diripc,fpname,
     & jac,fxname,hfact,hmax,hmin,htancf,ierror,ijac,ipc,ipivot,
     & ipl,ipoint,ipram,istart,it,iwrite,jpoint,jwrite,kwrite,lim,
     & lwrite,modnew,ncol,nrow,numcon,numlim,numstr,numtar,
     & nvar,relerr,slname,tc,tl,tr,work1,work2,work3,xc,xf,xit,xr)

c*********************************************************************72
c
cc P22_EDIT is an edit routine for problem 22.
c
c  this routine computes from the coefficients contained in x
c  the corresponding radial and angular displacement functions
c  u1 and u2, computes the corresponding position vector,
c  and calls the terminal plotter for a picture of the current arch.
c
      EXTERNAL FPNAME
      EXTERNAL FXNAME
      EXTERNAL SLNAME
      DIMENSION jac(NROW,NCOL)
      DIMENSION IPIVOT(NVAR)
      DIMENSION TC(NVAR)
      DIMENSION TL(NVAR)
      DIMENSION TR(NVAR)
      DIMENSION WORK1(NVAR)
      DIMENSION WORK2(NVAR)
      DIMENSION WORK3(NVAR)
      DIMENSION XC(NVAR)
      DIMENSION XF(NVAR)
      DIMENSION XR(NVAR)
      DIMENSION XVALUE(31)
      DIMENSION YVALUE(31)
      COMMON /AUXMEM/ THETA,COEF2,COEF3,COEF4,COEF5,IFUN22,
     & NGAUSS,GAUSCO(13),GAUSAB(13),PL(13),PLD(13),PLDD(13),
     & UX(2),UD(2),UDD(2),NCOF,ILM,IMU,INU,VALLM,VALMU,VALNU,BLOAD,
     & CLOAD,THICK,RADIUS,ALFA0,ALFA1,ALFA2
c
c  full edit for limit points and starting points
c  which have ipoint=1 or ipoint=5
c
      XLM=XR(ILM)
      XMU=XR(IMU)
      XNU=XR(INU)
      WRITE(*,*)' '
      WRITE(*,1000)
      WRITE(*,1040)XLM,XMU,XNU
      WRITE(*,*)' '
      WRITE(*,*)' '
      IF(IPOINT.NE.1.AND.IPOINT.NE.5)GO TO 30
      WRITE(*,1050)
      WRITE(*,*)' '
      NPLOT=21
      QLOAD=BLOAD
      DO J=1,NGAUSS
        T=GAUSAB(J)-XNU
        T2=T*T
        QLOAD=QLOAD+GAUSCO(J)*EXP(-XMU*T2)
      end do
      DT=2.0/FLOAT(NPLOT-1)
      RADIUS=10.0
      DTHETA=2.0*THETA/FLOAT(NPLOT-1)

      DO I=1,NPLOT
        T=-1.0+FLOAT(I-1)*DT
        call p22_vl (T,XR,NCOF,PL,PLD,PLDD,UX,UD,UDD)
        IF(IFUN22.EQ.1)XLOAD=BLOAD+XLM
        IF(IFUN22.EQ.2.AND.T.GT.XNU)XLOAD=BLOAD+XLM+XLM*(XNU-T)/XMU
        IF(IFUN22.EQ.2.AND.T.LE.XNU)XLOAD=BLOAD+XLM-XLM*(XNU-T)/XMU
        IF(IFUN22.EQ.3)XLOAD=XLM*(BLOAD+EXP(-XMU*(T-XNU)**2))/QLOAD
        RLOAD=XLOAD*COS(CLOAD*THETA*T)
        TLOAD=XLOAD*SIN(CLOAD*THETA*T)
        WRITE(*,1060)T,UX(1),UX(2),RLOAD,TLOAD
        RAD=RADIUS-UX(1)
        DTH=ATAN(-UX(2)/RAD)
        ANGLE=-THETA+FLOAT(I-1)*DTHETA+DTH
        XVALUE(I)=RAD*COS(ANGLE)
        YVALUE(I)=RAD*SIN(ANGLE)
      end do

      IYMAX=16
      IXMAX=61
      NVAL=NPLOT
      call DOTPLT(IXMAX,IYMAX,NVAL,YVALUE,XVALUE)
      return
c
c  center point value only for ipoint.ne.1.and.ipoint.ne.5
c
   30 T=0.0
      call p22_vl (T,XR,NCOF,PL,PLD,PLDD,UX,UD,UDD)
      WRITE(*,*)' '
      WRITE(*,1010)
      WRITE(*,1020)
      WRITE(*,1030)T,UX(1),UX(2)
      WRITE(*,*)' '
      return
 1000 FORMAT(17H PLOT OF SOLUTION)
 1010 FORMAT(24H CENTER POINT DEFLECTION)
 1020 FORMAT(39H      T           UX(1)           UX(2))
 1030 FORMAT(1X,3G14.6)
 1040 FORMAT(8H LAMBDA=,G14.6,4H MU=,G14.6,4H NU=,G14.6)
 1050 FORMAT(7X,1HT,13X,2HU1,12X,2HU2,7X,11HRADIAL LOAD,4X,9HTANG LOAD)
 1060 FORMAT(1X,5G14.6)
      END
      subroutine p22_bs (T,PL,PLD,PLDD)

c*********************************************************************72
c
cc P22_BS evaluates basis functions for the arch limit run.
c
c  this routine evaluates, for t between -1 and 1, the integrated
c  legendre polynomials used as basis functions for u1 and u2.  note
c  that while 7 terms are used for u1, only 6 are used for u2.
c  both u1 and u2 are approximated by polynomials beginning with
c  the quadratic.
c
      DIMENSION PL(13),PLD(13),PLDD(13)
      QDL=1.0
      QDDL=0.0
      QDC=T
      QDDC=1.0
      PL(1)=0.5*(T*T-1.0)
      PLD(1)=QDC
      PLDD(1)=QDDC
      AK=0.0

      DO K=2,7
        AK=AK+1.0
        PLD(K)=((AK+AK+1.0)*T*QDC-AK*QDL)/(AK+1.0)
        PLDD(K)=((AK+AK+1.0)*(T*QDDC+QDC)-AK*QDDL)/(AK+1.0)
        QDL=QDC
        QDC=PLD(K)
        QDDL=QDDC
        QDDC=PLDD(K)
        PL(K)=(T*QDC-QDL)/(AK+2.0)
      end do

      DO K=8,13
        KM7=K-7
        PL(K)=PL(KM7)
        PLD(K)=PLD(KM7)
        PLDD(K)=PLDD(KM7)
      end do

      return
      END
      subroutine p23_bs (T,PL,PLD,PLDD)

c*********************************************************************72
c
cc P23_BS is a new version of BS0022 with a different formula.
c
c  this routine evaluates, for t between -1 and 1, the integrated
c  legendre polynomials used as basis functions for u1 and u2.  note
c  that while 7 terms are used for u1, only 6 are used for u2.
c  both u1 and u2 are approximated by polynomials beginning with
c  the quadratic.
c
      DIMENSION PL(13),PLD(13),PLDD(13)
      QDL=1.0
      QDDL=0.0
      QDC=T
      QDDC=1.0
      FK=0.0

      DO K=1,7
        FK=FK+1.0
        PLD(K)=QDC
        PLDD(K)=QDDC
        FKFKP1=FK+FK+1.0
        QDN=(FKFKP1*T*QDC-FK*QDL)/(FK+1.0)
        QDDN=QDDL+FKFKP1*QDC
        PL(K)=(QDN-QDL)/(FKFKP1)
        QDL=QDC
        QDC=QDN
        QDDL=QDDC
        QDDC=QDDN
      end do

      DO K=8,13
        KM7=K-7
        PL(K)=PL(KM7)
        PLD(K)=PLD(KM7)
        PLDD(K)=PLDD(KM7)
      end do

      return
      END
      subroutine p22_vl (T,X,NCOF,PL,PLD,PLDD,UX,UD,UDD)

c*********************************************************************72
c
cc P22_VL evaluates finite element functions for the arch limit run.
c
c  this routine evaluates u1 and u2, as well as their first
c  and second derivatives, where u1 is stored as 7 coefficients
c  in x, and u2 as 6 coefficients.  the coefficients multiply
c  the legendre polynomials which are computed by bs0022.
c
      DIMENSION X(NCOF)
      DIMENSION PL(NCOF),PLD(NCOF),PLDD(NCOF),UX(2),UD(2),UDD(2)

      call p22_bs (T,PL,PLD,PLDD)

      DO I=1,2
        UX(I)=0.0
        UD(I)=0.0
        UDD(I)=0.0
      end do

      DO J=1,NCOF
        XJ=X(J)
        I=1
        IF(J.GT.7)I=2
        UX(I)=UX(I)+XJ*PL(J)
        UD(I)=UD(I)+XJ*PLD(J)
        UDD(I)=UDD(I)+XJ*PLDD(J)
      end do

      return
      END
      subroutine p22_ld (X,ID,Y)

c*********************************************************************72
c
cc P22_LD evaluates the load function for the arch limit run.
c
      DIMENSION X(16),Y(13)
      COMMON /AUXMEM/ THETA,COEF2,COEF3,COEF4,COEF5,IFUN22,
     & NGAUSS,GAUSCO(13),GAUSAB(13),PL(13),PLD(13),PLDD(13),
     & UX(2),UD(2),UDD(2),NCOF,ILM,IMU,INU,VALLM,VALMU,VALNU,BLOAD,
     & CLOAD,THICK,RADIUS,ALFA0,ALFA1,ALFA2
      IF(IFUN22.LT.1.OR.IFUN22.GT.3)IFUN22=1
      IF(IFUN22.EQ.1)CALL p22_lcon (X,ID,Y)
      IF(IFUN22.EQ.2)CALL p22_llin (X,ID,Y)
      IF(IFUN22.EQ.3)CALL p22_lexp (X,ID,Y)
      return
      END
      subroutine p22_lexp (X,ID,Y)

c*********************************************************************72
c
cc P22_LEXP evaluates load and load derivatives for the arch limit run.
c
c  id=1  compute load
c  id=2  compute d load/d lambda
c  id=3  compute d load/d mu
c  id=4  compute d2 load/d lambda d mu
c  id=5  compute d load/d nu
c  id=6  compute d2 load/d lambda d nu
c
c  load(x,lambda,mu,nu)=lambda*(bload+exp(-mu*(x-nu)**2))/s
c  where s=integral(-1 to 1) (bload+exp(-mu*(x-nu)**2))
c
      DIMENSION X(16),Y(13)
      COMMON /AUXMEM/ THETA,COEF2,COEF3,COEF4,COEF5,IFUN22,
     & NGAUSS,GAUSCO(13),GAUSAB(13),PL(13),PLD(13),PLDD(13),
     & UX(2),UD(2),UDD(2),NCOF,ILM,IMU,INU,VALLM,VALMU,VALNU,BLOAD,
     & CLOAD,THICK,RADIUS,ALFA0,ALFA1,ALFA2

      DO I=1,NCOF
        Y(I)=0.0
      end do

      XLM=X(ILM)
      XMU=X(IMU)
      XNU=X(INU)
      S1=2.0*BLOAD
      S2=0.0
      S3=0.0
      DO 30 J=1,NGAUSS
        T1=GAUSAB(J)-XNU
        T2=T1*T1
        EXPON=EXP(-XMU*T2)
        S1=S1+GAUSCO(J)*EXPON
        IF(ID.LE.2)GO TO 30
        IF(ID.GT.4)GO TO 20
        S2=S2+T2*GAUSCO(J)*EXPON
        GO TO 30
   20   S3=S3+2.0*XMU*T1*GAUSCO(J)*EXPON
   30   CONTINUE

      DO J=1,NGAUSS

        T=GAUSAB(J)
        GC=GAUSCO(J)
        call p22_bs (T,PL,PLD,PLDD)
        T1=T-XNU
        T2=T1*T1
        EXPON=EXP(-XMU*T2)
        IF(ID.EQ.1)PLOAD=XLM*(BLOAD+EXPON)/S1
        IF(ID.EQ.2)PLOAD=(BLOAD+EXPON)/S1
        IF(ID.EQ.3)PLOAD=XLM*(-T2*EXPON+(BLOAD+EXPON)*(S2/S1))/S1
        IF(ID.EQ.4)PLOAD=(-T2*EXPON+(BLOAD+EXPON)*(S2/S1))/S1
        IF(ID.EQ.5)
     &      PLOAD=XLM*(2.0*XMU*T1*EXPON-(BLOAD+EXPON)*(S3/S1))/S1
        IF(ID.EQ.6)PLOAD=(2.0*XMU*T1*EXPON-(BLOAD+EXPON)*(S3/S1))/S1

        DO K=1,NCOF
          IF(K.LE.7)TERM=PLOAD*PL(K)*COS(CLOAD*THETA*T)
          IF(K.GT.7)TERM=PLOAD*PL(K)*SIN(CLOAD*THETA*T)
          Y(K)=Y(K)+2.0*ALFA2*THETA*TERM*GC
        end do

      end do

      return
      END
      subroutine p22_lcon (X,ID,Y)

c*********************************************************************72
c
cc P22_LCON evaluates the load or load derivatives for the arch limit run.
c
c  id=1  compute load
c  id=2  compute d load/d lambda
c  id=3  compute d load/d mu
c  id=4  compute d2 load/d lambda d mu
c  id=5  compute d load/d nu
c  id=6  compute d2 load/d lambda d nu
c
c  load(x,lambda,mu,nu)=lambda
c
      DIMENSION X(16),Y(13)
      COMMON /AUXMEM/ THETA,COEF2,COEF3,COEF4,COEF5,IFUN22,
     & NGAUSS,GAUSCO(13),GAUSAB(13),PL(13),PLD(13),PLDD(13),
     & UX(2),UD(2),UDD(2),NCOF,ILM,IMU,INU,VALLM,VALMU,VALNU,BLOAD,
     & CLOAD,THICK,RADIUS,ALFA0,ALFA1,ALFA2

      DO I=1,NCOF
        Y(I)=0.0
      end do

      XLM=X(ILM)
      XMU=X(IMU)
      XNU=X(INU)

      DO J=1,NGAUSS

        T=GAUSAB(J)
        GC=GAUSCO(J)
        call p22_bs (T,PL,PLD,PLDD)
        IF(ID.EQ.1)PLOAD=XLM+BLOAD
        IF(ID.EQ.2)PLOAD=1.0
        IF(ID.EQ.3)PLOAD=0.0
        IF(ID.EQ.4)PLOAD=0.0
        IF(ID.EQ.5)PLOAD=0.0
        IF(ID.EQ.6)PLOAD=0.0

        DO K=1,NCOF
          IF(K.LE.7)TERM=PLOAD*PL(K)*COS(CLOAD*THETA*T)
          IF(K.GT.7)TERM=PLOAD*PL(K)*SIN(CLOAD*THETA*T)
          Y(K)=Y(K)+2.0*ALFA2*THETA*TERM*GC
        end do

      end do

      return
      END
      subroutine p22_llin (X,ID,Y)

c*********************************************************************72
c
cc P22_LLIN evaluates the load and load derivatives for the arch limit run.
c
c  simpler version of ld0022
c
c  id=1  compute load
c  id=2  compute d load/d lambda
c  id=3  compute d load/d mu
c  id=4  compute d2 load/d lambda d mu
c  id=5  compute d load/d nu
c  id=6  compute d2 load/d lambda d nu
c
c  load(x,lambda,mu,nu)=bload+xlm+xlm*(xnu-t)/xmu
c  for xnu-xmu.le.t.le.xnu
c  load(x,lambda,mu,nu)=bload+xlm-xlm*(xnu-t)/xmu
c  for xnu.le.t.le.xnu+xmu
c
      DIMENSION X(16),Y(13)
      COMMON /AUXMEM/ THETA,COEF2,COEF3,COEF4,COEF5,IFUN22,
     & NGAUSS,GAUSCO(13),GAUSAB(13),PL(13),PLD(13),PLDD(13),
     & UX(2),UD(2),UDD(2),NCOF,ILM,IMU,INU,VALLM,VALMU,VALNU,BLOAD,
     & CLOAD,THICK,RADIUS,ALFA0,ALFA1,ALFA2

      DO I=1,NCOF
        Y(I)=0.0
      end do

      XLM=X(ILM)
      XMU=X(IMU)
      XNU=X(INU)
      DO 60 J=1,NGAUSS
        T=GAUSAB(J)
        GC=GAUSCO(J)
        call p22_bs (T,PL,PLD,PLDD)
        PLOAD=0.0
        IF(T.GT.XNU.AND.T.LE.(XNU+XMU))GO TO 20
        IF(T.LE.XNU.AND.T.GE.(XNU-XMU))GO TO 30
        GO TO 40
   20   IF(ID.EQ.1)PLOAD=BLOAD+XLM+XLM*(XNU-T)/XMU
        IF(ID.EQ.2)PLOAD=1.0+(XNU-T)/XMU
        IF(ID.EQ.3)PLOAD=-XLM*(XNU-T)/(XMU*XMU)
        IF(ID.EQ.4)PLOAD=-(XNU-T)/(XMU*XMU)
        IF(ID.EQ.5)PLOAD=XLM/XMU
        IF(ID.EQ.6)PLOAD=1.0/XMU
        GO TO 40
   30   IF(ID.EQ.1)PLOAD=BLOAD+XLM-XLM*(XNU-T)/XMU
        IF(ID.EQ.2)PLOAD=1.0-(XNU-T)/XMU
        IF(ID.EQ.3)PLOAD=+XLM*(XNU-T)/(XMU*XMU)
        IF(ID.EQ.4)PLOAD=+(XNU-T)/(XMU*XMU)
        IF(ID.EQ.5)PLOAD=-XLM/XMU
        IF(ID.EQ.6)PLOAD=-1.0/XMU
        GO TO 40

   40   DO K=1,NCOF
          IF(K.LE.7)TERM=PLOAD*PL(K)*COS(CLOAD*THETA*T)
          IF(K.GT.7)TERM=PLOAD*PL(K)*SIN(CLOAD*THETA*T)
          Y(K)=Y(K)+2.0*ALFA2*THETA*TERM*GC
        end do

   60   CONTINUE
      return
      END
      subroutine p22_nvar ( option, nvar )

c*********************************************************************72
c
cc P22_NVAR sets the number of variables for problem 22.
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
c    Input, integer OPTION, the option chosen for this problem.
c    For some problems, several options are available.  At least,
c    OPTION = 1 is always legal.
c
c    Output, integer NVAR, the number of variables.
c
      implicit none

      integer nvar
      integer option

      nvar = 16

      return
      end
      subroutine p22_option_num ( option_num )

c*********************************************************************72
c
cc P22_OPTION_NUM returns the number of options for problem 22.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    16 September 2008
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, integer OPTION_NUM, the number of options.
c
      implicit none

      integer option_num

      option_num = 10

      return
      end
      subroutine p22_start ( option, nvar, x )

c*********************************************************************72
c
cc P22_START returns a starting point for problem 22.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    19 September 2008
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer OPTION, the option index.
c
c    Input, integer NVAR, the number of variables.
c
c    Output, double precision X(NVAR), the starting point.
c
      implicit none

      integer nvar

      integer option
      double precision x(nvar)

      if ( ? ) then

      else

        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'P22_START - Fatal error!'
        write ( *, '(a,i8)' ) '  Unrecognized option  = ', option
        stop

      end if

      return
      end
      subroutine p22_stepsize ( option, h, hmin, hmax )

c*********************************************************************72
c
cc P22_STEPSIZE returns step sizes for problem 22.
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
c    Input, integer OPTION, the option index.
c
c    Output, double precision H, HMIN, HMAX, suggested values for the
c    initial step, the minimum step, and the maximum step.
c
      implicit none

      double precision h
      double precision hmax
      double precision hmin
      integer option

      h =    0.250D+00
      hmin = 0.001D+00
      hmax = 1.000D+00

      return
      end
      subroutine p22_title ( option, title )

c*********************************************************************72
c
cc P22_TITLE sets the title for problem 22.
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
c    Input, integer OPTION, the option index.
c
c    Output, character*(*) TITLE, the title of the problem.
c    TITLE will never be longer than 80 characters.
c
      implicit none

      integer option
      character*(*) title

      if ( option .eq. 1 ) then
        title = 'Arch limit run on lambda'
      else if ( option .eq. 2 ) then
        title = 'Arch limit run on lambda'
      else if ( option .eq. 3 ) then
        title = 'Arch limit run on lambda'
      else if ( option .eq. 4 ) then
        title = 'Arch limit run on lambda'
      else if ( option .eq. 5 ) then
        title = 'Arch limit run on lambda'
      else if ( option .eq. 6 ) then
        title = 'Arch limit run on lambda'
      else if ( option .eq. 7 ) then
        title = 'Arch limit run on lambda'
      else if ( option .eq. 8 ) then
        title = 'Arch limit run on lambda'
      else if ( option .eq. 9 ) then
        title = 'Arch limit run on lambda'
      else if ( option .eq. 10 ) then
        title = 'Arch limit run on lambda'
      else
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'P22_TITLE - Fatal error!'
        write ( *, '(a)' ) '  Unrecognized option number.'
        stop
      end if

      return
      end
      subroutine p23_data (DIRIPC,HFACT,HMAX,HMIN,HTAN,option,
     & IPC,IPRAM,ISIZE,IT,LIM,MAXCON,
     & MAXLIM,MAXSTP,MAXTAR,NCOL,NROW,NVAR,RELERR,XIT,RWORK)

c*********************************************************************72
c
cc P23_DATA sets parameters for the arch limit run on mu after lambda.
c
c 23.1 arch limit run on mu after lambda
c
c 23.2 reference
c
c  werner c rheinboldt
c  work in progress
c
c 23.3 the function
c
c  see 22.3 for the original function and variables.  following
c  a run of problem 22, this program may be run.  it is assumed
c  that on input, rwork contains a solution point (x,lambda,mu,nu)
c  and the tangent at that solution point, tan, such that
c  tan(lambda)=0.0.  the new set of variables is then
c  (x,tan,lambda,mu,nu) and the new equations are
c
c  coefficient equations of problem 22     (ncof equations)
c  jacobian of 22 * tan = 0.0              (ncof equations)
c  tan*tan-1 =0                            (1 equation)
c  nu-valnu=0       (valnu given)          (1 equation)
c
c 23.4 options
c
c  option=1  vertical load, mu=0.0, nu=0.0
c  option=2  radial load,   mu=0.0, nu=0.0
c  option=3  vertical load, mu=0.1, nu=.2304583159
c  option=4  radial load,   mu=0.1, nu=.2304583159
c  option=5  vertical load, mu=0.1, nu=.6423493394
c  option=6  radial load,   mu=0.1, nu=.6423493394
c  option=7  vertical load, mu=0.1, nu=.9841830547
c  option=8  radial load,   mu=0.1, nu=.9841830547
c  option=9  vertical load, mu=0.1, nu=0.0
c  option=10 radial load,   mu=0.1, nu=0.0
c
c 23.5 target points
c
c  no target points are sought
c
c 23.6 limit points
c
c  limit points are sought in mu
c
c 23.7  bifurcation points
c
c 23.8 comments
c
      DIMENSION RWORK(ISIZE)
      COMMON /AUXMEM/ THETA,COEF2,COEF3,COEF4,COEF5,IFUN22,
     & NGAUSS,GAUSCO(13),GAUSAB(13),PL(13),PLD(13),PLDD(13),
     & UX(2),UD(2),UDD(2),NCOF,ILM,IMU,INU,VALLM,VALMU,VALNU,BLOAD,
     & CLOAD,THICK,RADIUS,ALFA0,ALFA1,ALFA2

      NCOF=13
      NCOF2=NCOF+NCOF
      ILM=NCOF2+1
      IMU=NCOF2+2
      INU=NCOF2+3
      NVAR=NCOF2+3

      DIRIPC=1.0
      HFACT=2.0
      IPC=IMU
      IT=0
      LIM=IMU
      MAXCON=20
      MAXLIM=1
      MAXSTP=20
      MAXTAR=20
      NCOL=NVAR
      NROW=NVAR
      RELERR=1.0E-3
      XIT=0.0
c
c  save current values of lambda, mu, nu
c
      VALLM=RWORK(14)
      VALMU=RWORK(15)
      VALNU=RWORK(16)
c
c  copy old tangent as second set of variables
c  make sure this index matches that computed in contrl
c
      NVOLD=NCOF+3
      JTL=3*NVOLD+NVOLD*NVOLD
      DO I=1,NCOF
        RWORK(NCOF+I)=RWORK(JTL+I)
      end do
c
c  lambda, mu, nu stored last
c
      RWORK(ILM)=VALLM
      RWORK(IMU)=VALMU
      RWORK(INU)=VALNU
      return
      END
      subroutine p23_fun ( option, nvar, x, fx )

c*********************************************************************72
c
cc P23_FUN evaluates the arch function for limit run on mu after lambda.
c
      integer nvar

      double precision fx(nvar-1)
      integer option
      double precision x(nvar)

      DIMENSION Y(13)
      COMMON /AUXMEM/ THETA,COEF2,COEF3,COEF4,COEF5,IFUN22,
     & NGAUSS,GAUSCO(13),GAUSAB(13),PL(13),PLD(13),PLDD(13),
     & UX(2),UD(2),UDD(2),NCOF,ILM,IMU,INU,VALLM,VALMU,VALNU,BLOAD,
     & CLOAD,THICK,RADIUS,ALFA0,ALFA1,ALFA2
c
c  compute first ncof entries of function
c  as original function
c
      call p22_gx (NVAR,X,FX)
      call p22_ld (X,1,Y)
      DO IEQN=1,NCOF
        FX(IEQN)=FX(IEQN)-Y(IEQN)
      end do
c
c  next ncof equations are jacobian of original
c  problem times tangent vector
c
      ITAN=NCOF+1
      call p22_gpa (X,X(ITAN),FX(ITAN))
c
c  ncof+ncof+1 equation is tan*tan-1=0
c
      IEQN=NCOF+NCOF+1
      ILO=ITAN
      IHI=NCOF+NCOF
      SUM=-1.0
      DO IVAR=ILO,IHI
        SUM=SUM+X(IVAR)*X(IVAR)
      end do
      FX(IEQN)=SUM
c
c  ncof+ncof+2 equation fixes nu
c
      IEQN=NCOF+NCOF+2
      FX(IEQN)=X(INU)-VALNU
      return
      END
      subroutine p23_jac ( option, nvar, x, jac )

c*********************************************************************72
c
cc P23_JAC evaluates the arch limit run on mu after lambda jacobian.
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
c    Input, integer OPTION, the option index.
c
c    Input, integer NVAR, the number of variables.
c
c    Input, double precision X(NVAR), the argument of the jacobian.
c
c    Output, double precision JAC(NVAR,NVAR), the jacobian matrix evaluated
c    at X.  The NVAR-th row is not set by this routine.
c
      implicit none

      integer nvar

      integer option
      double precision jac(nvar,nvar)
      double precision x(nvar)

      DIMENSION Y(13),Z(13),W(13)
      COMMON /AUXMEM/ THETA,COEF2,COEF3,COEF4,COEF5,IFUN22,
     & NGAUSS,GAUSCO(13),GAUSAB(13),PL(13),PLD(13),PLDD(13),
     & UX(2),UD(2),UDD(2),NCOF,ILM,IMU,INU,VALLM,VALMU,VALNU,BLOAD,
     & CLOAD,THICK,RADIUS,ALFA0,ALFA1,ALFA2
c
c  zero out jacobian
c
      do j = 1, ncol
        do i = 1, nrow
          jac(i,j) = 0.0E+00
        end do
      end do
c
c  store dload/dlambda in y
c        dload/dmu in z
c        dload/dnu in w
c
      call p22_ld (X,2,Y)
      call p22_ld (X,3,Z)
      call p22_ld (X,5,W)
c
c  put load derivatives into jacobian
c
      DO IEQN=1,NCOF
        jac(IEQN,ILM)=-Y(IEQN)
        jac(IEQN,IMU)=-Z(IEQN)
        jac(IEQN,INU)=-W(IEQN)
      end do

      ITAN=NCOF+1

      DO J=1,NCOF

        IVAR=J
        call p22_gp (X,IVAR,Y)
        call p22_g2a (X,IVAR,X(ITAN),Z)
        JVAR=IVAR+13

        DO IEQN=1,NCOF
          jac(IEQN,IVAR)=Y(IEQN)
          JEQN=IEQN+13
          jac(JEQN,IVAR)=Z(IEQN)
          jac(JEQN,JVAR)=Y(IEQN)
        end do

      end do
c
c  tangent normalization equation
c
      ILO=NCOF+1
      IHI=NCOF+NCOF
      IEQN=NCOF+NCOF+1

      DO IVAR=ILO,IHI
        jac(IEQN,IVAR)=2.0*X(IVAR)
      end do
c
c  last equation fixing nu
c
      IEQN=NCOF+NCOF+2
      jac(IEQN,INU)=1.0
      return
      END
      subroutine p23_nvar ( option, nvar )

c*********************************************************************72
c
cc P23_NVAR sets the number of variables for problem 23.
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
c    Input, integer OPTION, the option chosen for this problem.
c    For some problems, several options are available.  At least,
c    OPTION = 1 is always legal.
c
c    Output, integer NVAR, the number of variables.
c
      implicit none

      integer nvar
      integer option

      nvar = 29

      return
      end
      subroutine p23_option_num ( option_num )

c*********************************************************************72
c
cc P23_OPTION_NUM returns the number of options for problem 23.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    16 September 2008
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, integer OPTION_NUM, the number of options.
c
      implicit none

      integer option_num

      option_num = 10

      return
      end
      subroutine p23_start ( option, nvar, x )

c*********************************************************************72
c
cc P23_START returns a starting point for problem 23.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    19 September 2008
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer OPTION, the option index.
c
c    Input, integer NVAR, the number of variables.
c
c    Output, double precision X(NVAR), the starting point.
c
      implicit none

      integer nvar

      integer option
      double precision x(nvar)

      if ( ? ) then

      else

        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'P23_START - Fatal error!'
        write ( *, '(a,i8)' ) '  Unrecognized option  = ', option
        stop

      end if

      return
      end
      subroutine p23_stepsize ( option, h, hmin, hmax )

c*********************************************************************72
c
cc P23_STEPSIZE returns step sizes for problem 23.
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
c    Input, integer OPTION, the option index.
c
c    Output, double precision H, HMIN, HMAX, suggested values for the
c    initial step, the minimum step, and the maximum step.
c
      implicit none

      double precision h
      double precision hmax
      double precision hmin
      integer option

      h =    0.250D+00
      hmin = 0.001D+00
      hmax = 4.000D+00

      return
      end
      subroutine p23_title ( option, title )

c*********************************************************************72
c
cc P23_TITLE sets the title for problem 23.
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
c    Input, integer OPTION, the option index.
c
c    Output, character*(*) TITLE, the title of the problem.
c    TITLE will never be longer than 80 characters.
c
      implicit none

      integer option
      character*(*) title

      if ( option .eq. 1 ) then
        title = 'Arch limit run on mu after lambda'
      else if ( option .eq. 2 ) then
        title = 'Arch limit run on mu after lambda'
      else if ( option .eq. 3 ) then
        title = 'Arch limit run on mu after lambda'
      else if ( option .eq. 4 ) then
        title = 'Arch limit run on mu after lambda'
      else if ( option .eq. 5 ) then
        title = 'Arch limit run on mu after lambda'
      else if ( option .eq. 6 ) then
        title = 'Arch limit run on mu after lambda'
      else if ( option .eq. 7 ) then
        title = 'Arch limit run on mu after lambda'
      else if ( option .eq. 8 ) then
        title = 'Arch limit run on mu after lambda'
      else if ( option .eq. 9 ) then
        title = 'Arch limit run on mu after lambda'
      else if ( option .eq. 10 ) then
        title = 'Arch limit run on mu after lambda'
      else
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'P23_TITLE - Fatal error!'
        write ( *, '(a)' ) '  Unrecognized option number.'
        stop
      end if

      return
      end
      subroutine p24_data (DIRIPC,HFACT,HMAX,HMIN,HTAN,option,
     & IPC,IPRAM,ISIZE,IT,LIM,MAXCON,
     & MAXLIM,MAXSTP,MAXTAR,NCOL,NROW,NVAR,RELERR,XIT,RWORK)

c*********************************************************************72
c
cc P24_DATA sets parameters for the arch limit run on nu after lambda and mu.
c
c 24.1 arch limit run on nu after lambda and mu
c
c 24.2 reference
c
c  werner c rheinboldt
c  work in progress
c
c 24.3 the function
c
c  repeating the process that generated problem 23 from 22,
c  it is assumed on input that rwork contains a solution point
c  (x,tan1,lambda,mu,nu) of problem 23, as well as the tangent
c  tan2 for this problem.  the new set of variables is then
c  (x,tan1,tan2,lambda,mu,nu), and the new equations are
c
c  coefficient equations of 22        ncof equations
c  jacobian of 22 * tan1 = 0.0        ncof equations
c  jacobian of 23 * tan2 = 0.0        2*ncof+1 equations
c  tan1*tan1-1.0=0.0                  1 equation
c  tan2*tan2-1.0=0.0                  1 equation
c
c 24.4 options
c
c  option=1  vertical load, mu=0.0, nu=0.0
c  option=2  radial load,   mu=0.0, nu=0.0
c  option=3  vertical load, mu=0.1, nu=.2304583159
c  option=4  radial load,   mu=0.1, nu=.2304583159
c  option=5  vertical load, mu=0.1, nu=.6423493394
c  option=6  radial load,   mu=0.1, nu=.6423493394
c  option=7  vertical load, mu=0.1, nu=.9841830547
c  option=8  radial load,   mu=0.1, nu=.9841830547
c  option=9  vertical load, mu=0.1, nu=0.0
c  option=10 radial load,   mu=0.1, nu=0.0
c
c 24.5 target points
c
c  no target points are sought
c
c 24.6 limit points
c
c  limit points in nu are sought
c
c 24.7 bifurcation points
c
c 24.8 comments
c
      DIMENSION RWORK(ISIZE)
      COMMON /AUXMEM/ THETA,COEF2,COEF3,COEF4,COEF5,IFUN22,
     & NGAUSS,GAUSCO(13),GAUSAB(13),PL(13),PLD(13),PLDD(13),
     & UX(2),UD(2),UDD(2),NCOF,ILM,IMU,INU,VALLM,VALMU,VALNU,BLOAD,
     & CLOAD,THICK,RADIUS,ALFA0,ALFA1,ALFA2

      NCOF=13
      NCOF2=NCOF+NCOF
      NCOF4=NCOF2+NCOF2
      ILM=NCOF4+2
      IMU=NCOF4+3
      INU=NCOF4+4
      NVAR=NCOF4+4

      DIRIPC=1.0
      HFACT=3.0
      IPC=INU
      IT=0
      LIM=INU
      MAXCON=30
      MAXLIM=1
      MAXSTP=30
      MAXTAR=30
      NCOL=NVAR
      NROW=NVAR
      RELERR=1.0E-3
      XIT=0.0
c
c  save current values of lambda, mu, nu
c
      VALLM=RWORK(27)
      VALMU=RWORK(28)
      VALNU=RWORK(29)
c
c  copy old tangent as second set of variables
c  make sure index jtl matches that computed in contrl
c
      NVOLD=NCOF2+3
      JTL=3*NVOLD+NVOLD*NVOLD
      IHI=NCOF2+1

      DO I=1,IHI
        INDEX1=NCOF2+I
        INDEX2=JTL+I
        RWORK(INDEX1)=RWORK(INDEX2)
      end do

      RWORK(ILM)=VALLM
      RWORK(IMU)=VALMU
      RWORK(INU)=VALNU
      return
      END
      subroutine p24_fun ( option, nvar, x, fx )

c*********************************************************************72
c
cc P24_FUN evaluates the arch function for limit run on nu after mu after lambda.
c
      integer nvar

      double precision fx(nvar-1)
      integer option
      double precision x(nvar)

      DIMENSION W(13),V(13),Y(13)
      COMMON /AUXMEM/ THETA,COEF2,COEF3,COEF4,COEF5,IFUN22,
     & NGAUSS,GAUSCO(13),GAUSAB(13),PL(13),PLD(13),PLDD(13),
     & UX(2),UD(2),UDD(2),NCOF,ILM,IMU,INU,VALLM,VALMU,VALNU,BLOAD,
     & CLOAD,THICK,RADIUS,ALFA0,ALFA1,ALFA2
c
c  set original function
c
      call p22_gx (NVAR,X,FX)
      call p22_ld (X,1,W)

      DO I=1,NCOF
        FX(I)=FX(I)-W(I)
      end do

      ITAN1=NCOF+1
      ITAN2=NCOF+NCOF+1
      ITAN3=NCOF+NCOF+NCOF+1
c
c  set entries 14 through 26
c
      call p22_gpa (X,X(ITAN1),FX(ITAN1))
c
c  set entries 27 through 39
c
      call p22_gpa (X,X(ITAN2),Y)
      call p22_ld (X,2,V)
      DO I=27,39
        FX(I)=Y(I-26)-X(53)*V(I-26)
      end do
c
c  set entries 40 through 52
c
      call p22_gpa (X,X(ITAN3),Y)
      call p22_g2ab (X,X(ITAN1),X(ITAN2),V)
      DO I=40,52
        IM39=I-39
        FX(I)=V(IM39)+Y(IM39)
      end do
c
c  set equation 53
c  tangent1 normalization
c
      SUM=-1.0
      DO I=1,NCOF
        SUM=SUM+X(I+13)**2
      end do
      FX(53)=SUM
c
c  set equation 54
c  normalization equation of 23 * tangent2
c
      SUM=0.0
      DO I=1,13
        SUM=SUM+2.0*X(I+39)*X(I+ITAN1)
      end do
      FX(54)=SUM
c
c  set equation 55
c  tangent2 normalization
c
      SUM=-1.0
      DO I=27,53
        SUM=SUM+X(I)**2
      end do
      FX(55)=SUM
      return
      END
      subroutine p24_jac ( option, nvar, x, jac )

c*********************************************************************72
c
cc P24_JAC evaluates the arch limit run on nu after mu after lambda jacobian.
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
c    Input, integer OPTION, the option index.
c
c    Input, integer NVAR, the number of variables.
c
c    Input, double precision X(NVAR), the argument of the jacobian.
c
c    Output, double precision JAC(NVAR,NVAR), the jacobian matrix evaluated
c    at X.  The NVAR-th row is not set by this routine.
c
      implicit none

      integer nvar

      integer option
      double precision jac(nvar,nvar)
      double precision x(nvar)

      DIMENSION W(13),V(13),Y(13),Z(13)
      COMMON /AUXMEM/ THETA,COEF2,COEF3,COEF4,COEF5,IFUN22,
     & NGAUSS,GAUSCO(13),GAUSAB(13),PL(13),PLD(13),PLDD(13),
     & UX(2),UD(2),UDD(2),NCOF,ILM,IMU,INU,VALLM,VALMU,VALNU,BLOAD,
     & CLOAD,THICK,RADIUS,ALFA0,ALFA1,ALFA2
c
c  zero out jacobian
c
      do j = 1, ncol
        do i = 1, nrow
          jac(i,j) = 0.0E+00
        end do
      end do
c
c  get load derivatives
c
      call p22_ld (X,2,V)
      call p22_ld (X,3,Y)
      call p22_ld (X,5,Z)

      DO I=1,NCOF
        jac(I,54)=-V(I)
        jac(I,55)=-Y(I)
        jac(I,56)=-Z(I)
        IP13=I+13
        W(I)=X(IP13)
        IP26=I+26
        Z(I)=X(IP26)
        jac(IP26,53)=-V(I)
      end do

      DO 50 J=1,13
        IVAR=J
        call p22_gp (X,IVAR,Y)
        call p22_g2a (X,IVAR,W,V)
        JP13=J+13
        JP26=J+26
        JP39=J+39
        T=2.0*W(J)
        jac(53,JP13)=T
        jac(54,JP39)=T
        DO I=1,13
          T=Y(I)
          jac(I,J)=T
          IP13=I+13
          jac(IP13,JP13)=T
          IP26=I+26
          jac(IP26,JP26)=T
          IP39=I+39
          jac(IP39,JP39)=T
          T=V(I)
          jac(IP13,J)=T
          jac(IP26,JP39)=T
        end do

   50   CONTINUE

      DO J=1,13
        JD=J
        call G2A022(X,JD,Z,Y)
        call G3AB22(X,JD,W,Z,V)
        JP13=J+13
        DO I=1,13
          T=Y(I)
          IP26=I+26
          IP39=I+39
          jac(IP26,J)=T
          jac(IP39,JP13)=T
          jac(IP39,J)=V(I)
        end do
        T=2.0*W(J)
        jac(53,JP13)=T
        JP26=J+26
        jac(55,JP26)=2.0*Z(J)
        JP39=J+39
        jac(54,JP39)=T
      end do

      call p22_ld (X,4,W)
      call p22_ld (X,6,V)
      S=X(53)
      jac(55,53)=2.0*S

      DO J=1,13
        JP26=J+26
        jac(JP26,55)=-S*W(J)
        jac(JP26,56)=-S*V(J)
        JP39=J+39
        Z(J)=X(JP39)
      end do

      DO J=1,13
        IVAR=J
        call p22_gp (X,IVAR,Z,Y)
        DO I=1,13
          IP39=I+39
          jac(IP39,J)=jac(IP39,J)+Y(J)
        end do
        JP13=J+13
        JP39=J+39
        T=2.0*Y(J)
        jac(54,JP13)=T
        jac(55,JP39)=T
      end do

      return
      END
      subroutine p22_gx (NVAR,X,Y)

c*********************************************************************72
c
cc P22_GX is an auxilliary function for the arch limit run.
c
c  this routine evaluates the equilibrium integrals involving the
c  coefficients of the polynomials which describe the shape of the arch.
c
      DIMENSION X(NVAR),Y(13)
      COMMON /AUXMEM/ THETA,COEF2,COEF3,COEF4,COEF5,IFUN22,
     & NGAUSS,GAUSCO(13),GAUSAB(13),PL(13),PLD(13),PLDD(13),
     & UX(2),UD(2),UDD(2),NCOF,ILM,IMU,INU,VALLM,VALMU,VALNU,BLOAD,
     & CLOAD,THICK,RADIUS,ALFA0,ALFA1,ALFA2

      DO IEQN=1,NCOF
        Y(IEQN)=0.0
      end do

      DO J=1,NGAUSS
        T=GAUSAB(J)
        call p22_vl (T,X,NCOF,PL,PLD,PLDD,UX,UD,UDD)
        DO IEQN=1,NCOF
          IF(IEQN.LE.7)TERM=
     &    2.0*ALFA1*UDD(1)*PLDD(IEQN)/THETA**3
     &    +2.0*ALFA0*UD(1)*UD(2)*PLD(IEQN)/THETA**2
     &    -2.0*ALFA0*UX(1)*UD(1)*PLD(IEQN)/THETA
     &    +ALFA0**2*UD(1)**3*PLD(IEQN)/THETA**3
     &    +2.0*THETA*UX(1)*PL(IEQN)
     &    -2.0*UD(2)*PL(IEQN)
     &    -ALFA0*UD(1)*UD(1)*PL(IEQN)/THETA
          IF(IEQN.GT.7)TERM=
     &    ALFA0*UD(1)*UD(1)*PLD(IEQN)/THETA**2
     &    +2.0*UD(2)*PLD(IEQN)/THETA
     &    -2.0*UX(1)*PLD(IEQN)
          Y(IEQN)=Y(IEQN)+TERM*GAUSCO(J)
        end do
      end do

      return
      END
      subroutine p22_gp (X,IVAR,Y)

c*********************************************************************72
c
cc P22_GP is an auxilliary function for the arch limit run.
c
      DIMENSION X(16),Y(13)
      COMMON /AUXMEM/ THETA,COEF2,COEF3,COEF4,COEF5,IFUN22,
     & NGAUSS,GAUSCO(13),GAUSAB(13),PL(13),PLD(13),PLDD(13),
     & UX(2),UD(2),UDD(2),NCOF,ILM,IMU,INU,VALLM,VALMU,VALNU,BLOAD,
     & CLOAD,THICK,RADIUS,ALFA0,ALFA1,ALFA2

      DO IEQN=1,NCOF
        Y(IEQN)=0.0
      end do

      DO J=1,NGAUSS
        T=GAUSAB(J)
        GC=GAUSCO(J)
        call p22_vl (T,X,NCOF,PL,PLD,PLDD,UX,UD,UDD)
        DO IEQN=1,NCOF
          IF(IEQN.LE.7.AND.IVAR.LE.7)TERM=
     &    2.0*ALFA1*PLDD(IVAR)*PLDD(IEQN)/THETA**3
     &    +2.0*ALFA0*UD(2)*PLD(IVAR)*PLD(IEQN)/THETA**2
     &    -2.0*ALFA0*UX(1)*PLD(IVAR)*PLD(IEQN)/THETA
     &    -2.0*ALFA0*UD(1)*PL(IVAR)*PLD(IEQN)/THETA
     &    +3.0*ALFA0**2*UD(1)*UD(1)*PLD(IVAR)*PLD(IEQN)/THETA**3
     &    +2.0*THETA*PL(IVAR)*PL(IEQN)
     &    -2.0*ALFA0*UD(1)*PLD(IVAR)*PL(IEQN)/THETA
          IF(IEQN.LE.7.AND.IVAR.GT.7)TERM=
     &    2.0*ALFA0*UD(1)*PLD(IVAR)*PLD(IEQN)/THETA**2
     &    -2.0*PLD(IVAR)*PL(IEQN)
          IF(IEQN.GT.7.AND.IVAR.LE.7)TERM=
     &    2.0*ALFA0*UD(1)*PLD(IVAR)*PLD(IEQN)/THETA**2
     &    -2.0*PL(IVAR)*PLD(IEQN)
          IF(IEQN.GT.7.AND.IVAR.GT.7)TERM=
     &    2.0*PLD(IVAR)*PLD(IEQN)/THETA
          Y(IEQN)=Y(IEQN)+GC*TERM
        end do
      end do

      return
      END
      subroutine p22_gpa (X,W,Y)

c*********************************************************************72
c
cc P22_GPA computes y = J * w for the arch limit run.
c
c  this routine computes the product of the jacobian
c  of the coefficient matrix times the vector w, and
c  stores the result in y
c
      DIMENSION X(16),W(13),Y(13)
      COMMON /AUXMEM/ THETA,COEF2,COEF3,COEF4,COEF5,IFUN22,
     & NGAUSS,GAUSCO(13),GAUSAB(13),PL(13),PLD(13),PLDD(13),
     & UX(2),UD(2),UDD(2),NCOF,ILM,IMU,INU,VALLM,VALMU,VALNU,BLOAD,
     & CLOAD,THICK,RADIUS,ALFA0,ALFA1,ALFA2

      DO IEQN=1,NCOF
        Y(IEQN)=0.0
      end do

      DO J=1,NGAUSS
        T=GAUSAB(J)
        GC=GAUSCO(J)
        call p22_vl (T,X,NCOF,PL,PLD,PLDD,UX,UD,UDD)
        DO IEQN=1,NCOF
          DO IVAR=1,NCOF
            IF(IEQN.LE.7.AND.IVAR.LE.7)TERM=
     &      2.0*ALFA1*PLDD(IVAR)*PLDD(IEQN)/THETA**3
     &      +2.0*ALFA0*UD(2)*PLD(IVAR)*PLD(IEQN)/THETA**2
     &      -2.0*ALFA0*UX(1)*PLD(IVAR)*PLD(IEQN)/THETA
     &      -2.0*ALFA0*UD(1)*PL(IVAR)*PLD(IEQN)/THETA
     &      +3.0*ALFA0**2*UD(1)*UD(1)*PLD(IVAR)*PLD(IEQN)/THETA**3
     &      +2.0*THETA*PL(IVAR)*PL(IEQN)
     &      -2.0*ALFA0*UD(1)*PLD(IVAR)*PL(IEQN)/THETA
            IF(IEQN.LE.7.AND.IVAR.GT.7)TERM=
     &      2.0*ALFA0*UD(1)*PLD(IVAR)*PLD(IEQN)/THETA**2
     &      -2.0*PLD(IVAR)*PL(IEQN)
            IF(IEQN.GT.7.AND.IVAR.LE.7)TERM=
     &      2.0*ALFA0*UD(1)*PLD(IVAR)*PLD(IEQN)/THETA**2
     &      -2.0*PL(IVAR)*PLD(IEQN)
            IF(IEQN.GT.7.AND.IVAR.GT.7)TERM=
     &      2.0*PLD(IVAR)*PLD(IEQN)/THETA
            Y(IEQN)=Y(IEQN)+GC*TERM*W(IVAR)
          end do
        end do
      end do

      return
      END
      subroutine p22_g2a (X,IVAR2,W,Y)

c*********************************************************************72
c
cc P22_G2A is an auxilliary function for the arch limit run.
c
      DIMENSION X(16),W(13),Y(13)
      COMMON /AUXMEM/ THETA,COEF2,COEF3,COEF4,COEF5,IFUN22,
     & NGAUSS,GAUSCO(13),GAUSAB(13),PL(13),PLD(13),PLDD(13),
     & UX(2),UD(2),UDD(2),NCOF,ILM,IMU,INU,VALLM,VALMU,VALNU,BLOAD,
     & CLOAD,THICK,RADIUS,ALFA0,ALFA1,ALFA2

      DO IEQN=1,NCOF
        Y(IEQN)=0.0
      end do

      DO J=1,NGAUSS
        T=GAUSAB(J)
        call p22_vl (T,X,NCOF,PL,PLD,PLDD,UX,UD,UDD)
        DO IEQN=1,NCOF
          DO IVAR1=1,NCOF
            TERM=0.0
            IF(IEQN.LE.7.AND.IVAR1.LE.7.AND.IVAR2.LE.7)TERM=
     &      6.0*ALFA0**2*UD(1)*PLD(IEQN)*PLD(IVAR1)*PLD(IVAR2)/THETA**3
     &      -2.0*ALFA0*PL(IEQN)*PLD(IVAR1)*PLD(IVAR2)/THETA
     &      -2.0*ALFA0*PLD(IEQN)*PL(IVAR1)*PLD(IVAR2)/THETA
     &      -2.0*ALFA0*PLD(IEQN)*PLD(IVAR1)*PL(IVAR2)/THETA
            IF(IEQN.GT.7.AND.IVAR1.LE.7.AND.IVAR2.LE.7)
     &      TERM=2.0*ALFA0*PLD(IEQN)*PLD(IVAR1)*PLD(IVAR2)/THETA**2
            IF(IEQN.LE.7.AND.IVAR1.GT.7.AND.IVAR2.LE.7)
     &      TERM=2.0*ALFA0*PLD(IEQN)*PLD(IVAR1)*PLD(IVAR2)/THETA**2
            IF(IEQN.LE.7.AND.IVAR1.LE.7.AND.IVAR2.GT.7)
     &      TERM=2.0*ALFA0*PLD(IEQN)*PLD(IVAR1)*PLD(IVAR2)/THETA**2
            Y(IEQN)=Y(IEQN)+TERM*GAUSCO(J)*W(IVAR1)
          end do
        end do
      end do

      return
      END
      subroutine p22_g2ab (X,W,V,Y)

c*********************************************************************72
c
cc P22_G2AB is an auxilliary function for the arch limit run.
c
      DIMENSION X(16),W(13),V(13),Y(13)
      COMMON /AUXMEM/ THETA,COEF2,COEF3,COEF4,COEF5,IFUN22,
     & NGAUSS,GAUSCO(13),GAUSAB(13),PL(13),PLD(13),PLDD(13),
     & UX(2),UD(2),UDD(2),NCOF,ILM,IMU,INU,VALLM,VALMU,VALNU,BLOAD,
     & CLOAD,THICK,RADIUS,ALFA0,ALFA1,ALFA2

      DO IEQN=1,NCOF
        Y(IEQN)=0.0
      end do

      DO 50 J=1,NGAUSS
        T=GAUSAB(J)
        GC=GAUSCO(J)
        call p22_vl (T,X,NCOF,PL,PLD,PLDD,UX,UD,UDD)
        DO 40 IEQN=1,NCOF
          DO 30 IVAR1=1,NCOF
            DO 20 IVAR2=1,NCOF
              TERM=0.0
              IF(IEQN.LE.7.AND.IVAR1.LE.7.AND.IVAR2.LE.7)TERM=
     &        6.0*ALFA0**2*
     &        UD(1)*PLD(IEQN)*PLD(IVAR1)*PLD(IVAR2)/THETA**3
     &        -2.0*ALFA0*PL(IEQN)*PLD(IVAR1)*PLD(IVAR2)/THETA
     &        -2.0*ALFA0*PLD(IEQN)*PL(IVAR1)*PLD(IVAR2)/THETA
     &        -2.0*ALFA0*PLD(IEQN)*PLD(IVAR1)*PL(IVAR2)/THETA
              IF(IEQN.GT.7.AND.IVAR1.LE.7.AND.IVAR2.LE.7)TERM=
     &        2.0*ALFA0*PLD(IEQN)*PLD(IVAR1)*PLD(IVAR2)/THETA**2
              IF(IEQN.LE.7.AND.IVAR1.GT.7.AND.IVAR2.LE.7)TERM=
     &        2.0*ALFA0*PLD(IEQN)*PLD(IVAR1)*PLD(IVAR2)/THETA**2
              IF(IEQN.LE.7.AND.IVAR1.LE.7.AND.IVAR2.GT.7)TERM=
     &        2.0*ALFA0*PLD(IEQN)*PLD(IVAR1)*PLD(IVAR2)/THETA**2
   20         Y(IEQN)=Y(IEQN)+TERM*GC*W(IVAR1)*V(IVAR2)
   30       CONTINUE
   40     CONTINUE
   50   CONTINUE
        return
      END
      subroutine p22_g3ab (X,IVAR3,W,Z,Y)

c*********************************************************************72
c
cc P22_G3AB is an auxilliary function for the arch limit run.
c
      DIMENSION X(16),W(13),Z(13),Y(13)
      COMMON /AUXMEM/ THETA,COEF2,COEF3,COEF4,COEF5,IFUN22,
     & NGAUSS,GAUSCO(13),GAUSAB(13),PL(13),PLD(13),PLDD(13),
     & UX(2),UD(2),UDD(2),NCOF,ILM,IMU,INU,VALLM,VALMU,VALNU,BLOAD,
     & CLOAD,THICK,RADIUS,ALFA0,ALFA1,ALFA2

      DO IEQN=1,NCOF
        Y(IEQN)=0.0
      end do

      IF(IVAR3.GT.7)return

      DO J=1,NGAUSS
        T=GAUSAB(J)
        call p22_vl (T,X,NCOF,PL,PLD,PLDD,UX,UD,UDD)
        DO IEQN=1,7
          DO IVAR1=1,7
            DO IVAR2=1,7
              TERM=6.0*ALFA0**2*
     &        PLD(IEQN)*PLD(IVAR1)*PLD(IVAR2)*PLD(IVAR3)/THETA**3
              Y(IEQN)=Y(IEQN)+TERM*GAUSCO(J)*W(IVAR1)*Z(IVAR2)
            end do
          end do
        end do
      end do

      return
      END
      subroutine p24_nvar ( option, nvar )

c*********************************************************************72
c
cc P24_NVAR sets the number of variables for problem 24.
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
c    Input, integer OPTION, the option chosen for this problem.
c    For some problems, several options are available.  At least,
c    OPTION = 1 is always legal.
c
c    Output, integer NVAR, the number of variables.
c
      implicit none

      integer nvar
      integer option

      nvar = 56

      return
      end
      subroutine p24_option_num ( option_num )

c*********************************************************************72
c
cc P24_OPTION_NUM returns the number of options for problem 24.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    16 September 2008
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, integer OPTION_NUM, the number of options.
c
      implicit none

      integer option_num

      option_num = 10

      return
      end
      subroutine p24_start ( option, nvar, x )

c*********************************************************************72
c
cc P24_START returns a starting point for problem 24.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    19 September 2008
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer OPTION, the option index.
c
c    Input, integer NVAR, the number of variables.
c
c    Output, double precision X(NVAR), the starting point.
c
      implicit none

      integer nvar

      integer option
      double precision x(nvar)

      if ( ? ) then

      else

        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'P24_START - Fatal error!'
        write ( *, '(a,i8)' ) '  Unrecognized option  = ', option
        stop

      end if

      return
      end
      subroutine p24_stepsize ( option, h, hmin, hmax )

c*********************************************************************72
c
cc P24_STEPSIZE returns step sizes for problem 24.
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
c    Input, integer OPTION, the option index.
c
c    Output, double precision H, HMIN, HMAX, suggested values for the
c    initial step, the minimum step, and the maximum step.
c
      implicit none

      double precision h
      double precision hmax
      double precision hmin
      integer option

      h =     0.300D+00
      hmin =  0.001D+00
      hmax = 20.000D+00

      return
      end
      subroutine p24_title ( option, title )

c*********************************************************************72
c
cc P24_TITLE sets the title for problem 24.
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
c    Input, integer OPTION, the option index.
c
c    Output, character*(*) TITLE, the title of the problem.
c    TITLE will never be longer than 80 characters.
c
      implicit none

      integer option
      character*(*) title

      if ( option .eq. 1 ) then
        title = 'Arch limit run on nu after mu after lambda'
      else if ( option .eq. 2 ) then
        title = 'Arch limit run on nu after mu after lambda'
      else if ( option .eq. 3 ) then
        title = 'Arch limit run on nu after mu after lambda'
      else if ( option .eq. 4 ) then
        title = 'Arch limit run on nu after mu after lambda'
      else if ( option .eq. 5 ) then
        title = 'Arch limit run on nu after mu after lambda'
      else if ( option .eq. 6 ) then
        title = 'Arch limit run on nu after mu after lambda'
      else if ( option .eq. 7 ) then
        title = 'Arch limit run on nu after mu after lambda'
      else if ( option .eq. 8 ) then
        title = 'Arch limit run on nu after mu after lambda'
      else if ( option .eq. 9 ) then
        title = 'Arch limit run on nu after mu after lambda'
      else if ( option .eq. 10 ) then
        title = 'Arch limit run on nu after mu after lambda'
      else
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'P24_TITLE - Fatal error!'
        write ( *, '(a)' ) '  Unrecognized option number.'
        stop
      end if

      return
      end
      subroutine p25_data (DIRIPC,HFACT,HMAX,HMIN,HTAN,option,
     & IPC,IPRAM,ISIZE,IT,LIM,MAXCON,
     & MAXLIM,MAXSTP,MAXTAR,NCOL,NROW,NVAR,RELERR,XIT,RWORK)

c*********************************************************************72
c
cc P25_DATA sets parameters for the three-element arch.
c
c 25.1 3 element arch
c
c 25.2 reference
c
c  werner c rheinboldt
c  work in progress
c
c 25.3 the function
c
c  an arch consists of a segment of a circle.  the
c  radius of the circle is 10, and the angle subtended
c  is 30 degrees.  a third degree polynomial is used
c  to approximate the function u1=radial displacement
c  and u2=tangential displacement of the arch, when
c  subject to a load.  the coefficients of the polynomial
c  are the actual variables.  there are 2 coefficients
c  for u1, and 1 more used for u2.
c  to satisfy the boundary conditions at -1 and +1,
c  which are u1=u2=0, only quadratic and higher terms
c  of the integrated legendre polynomials are used.
c
c  constant and point loads are used.
c
c
c 25.4 options
c
c   1, point load.
c   2, constant load.
c
c 25.5 target points
c
c  no target points are sought
c
c 25.6 limit points
c
c
c 25.7  bifurcation points
c
c
c 25.8 comments
c
      real pi
      parameter ( pi = 3.1415926535E+00 )
      DIMENSION RWORK(ISIZE)
      COMMON /AUXMEM/ THETA,COEF2,COEF3,COEF4,COEF5,COEF6,IFUN22,VALLM,
     & VALNU

      NCOF=3
      ILM=4
      INU=5
      NVAR=NCOF+2

      DIRIPC=1.0
      HFACT=2.0
      IPC=ILM
      IT=0
      LIM=ILM
      MAXCON=30
      MAXLIM=1
      MAXSTP=30
      MAXTAR=30
      NCOL=NVAR
      NROW=NVAR
      RELERR=1.0E-3
      XIT=0.0
      THETA=PI/12.0
      EI=796.0
      EA=2.056*10.0**6
      H=.0684
      RADIUS=10.0
      ALFA0=H/RADIUS
      ALFA1=EI/(EA*RADIUS*RADIUS)
      ALFA2=EI/(EA*H*H)
      WRITE(*,*)' '
      WRITE(*,*)' '
      WRITE(*,1060)ALFA0,ALFA1,ALFA2
c
c  ifun22=1 for point load
c  ifun22=2 for constant load
c
      IF(option.EQ.1)IFUN22=1
      IF(option.EQ.2)IFUN22=2
      WRITE(*,*)' '
      IF(IFUN22.EQ.1)WRITE(*,1070)
      IF(IFUN22.EQ.2)WRITE(*,1080)
      COEF2=ALFA1/(THETA*THETA*THETA)
      COEF3=ALFA0/THETA
      COEF4=ALFA0/(THETA**2)
      COEF5=ALFA0*ALFA0/(THETA**3)
      COEF6=ALFA2*THETA
      WRITE(*,1010)COEF2
      WRITE(*,1020)COEF3
      WRITE(*,1030)COEF4
      WRITE(*,1040)COEF5
      WRITE(*,1050)COEF6
      VALLM=0.0
      VALNU=0.0
      DO I=1,NCOF
        RWORK(I)=0.0
      end do
      RWORK(ILM)=VALLM
      RWORK(INU)=VALNU
      return
 1010 FORMAT(1X,23HCOEF2=ALFA1/(THETA**3)=,G14.6)
 1020 FORMAT(1X,18HCOEF3=ALFA0/THETA=,G14.6)
 1030 FORMAT(1X,23HCOEF4=ALFA0/(THETA**2)=,G14.6)
 1040 FORMAT(1X,26HCOEF5=ALFA0**2/(THETA**3)=,G14.6)
 1050 FORMAT(1X,18HCOEF6=ALFA2*THETA=,G14.6)
 1060 FORMAT(7H ALFA0=,G14.6,7H ALFA1=,G14.6,7H ALFA2=,G14.6)
 1070 FORMAT(25H POINT LOAD FUNCTION USED)
 1080 FORMAT(28H CONSTANT LOAD FUNCTION USED)
      END
      subroutine p25_fun ( option, nvar, x, fx )

c*********************************************************************72
c
cc P25_FUN evaluates the three-element arch function.
c
      integer nvar

      double precision fx(nvar-1)
      integer option
      double precision x(nvar)

      DIMENSION XLOAD(3)
      COMMON /AUXMEM/ THETA,C2,C3,C4,C5,C6,IFUN22,VALLM,VALNU
      XLM=X(4)
      XNU=X(5)
      IF(IFUN22.EQ.1)XLOAD(1)=-0.5*C6*XLM*(XNU*XNU-1.0)*COS(THETA*XNU)
      IF(IFUN22.EQ.2)XLOAD(1)=-2.0*C6*XLM/3.0
      IF(IFUN22.EQ.1)
     & XLOAD(2)=-0.5*C6*XLM*(XNU*XNU-1.0)*XNU*COS(THETA*XNU)
      IF(IFUN22.EQ.2)XLOAD(2)=-C6*XLM/4.0
      IF(IFUN22.EQ.1)XLOAD(3)=-0.5*C6*XLM*(XNU*XNU-1.0)*SIN(THETA*XNU)
      IF(IFUN22.EQ.2)XLOAD(3)=0.0
      FX(1)=4.0*C2*X(1)
     &      +8.0*THETA*X(1)/15.0
     &      +4.0*C2*X(2)
     &      +0.4*C3*X(1)*X(1)
     &      +0.4*C5*X(1)*X(1)*X(1)
     &      +8.0*C4*X(2)*X(3)/15.0
     &      +14.0*C3*X(2)*X(2)/105.0
     &      +22.0*C5*X(2)*X(2)*X(1)/35.0
     &      +XLOAD(1)
      FX(2)=4.0*C2*X(1)
     &      +36.0*C2*X(2)/5.0
     &      +8.0*THETA*X(2)/105.0
     &      +4.0*X(3)/15.0
     &      +28.0*C3*X(1)*X(2)/105.0
     &      +8.0*C4*X(1)*X(3)/15.0
     &      +22.0*C5*X(1)*X(1)*X(2)/35.0
     &      +4.0*C5*X(2)*X(2)*X(2)/35.0
     &      +XLOAD(2)
      FX(3)=4.0*X(2)/15.0
     &      +8.0*C4*X(1)*X(2)/15.0
     &      +4.0*X(3)/(3.0*THETA)
     &      +XLOAD(3)
      FX(4)=XNU-VALNU
      return
      END
      subroutine p25_jac ( option, nvar, x, jac )

c*********************************************************************72
c
cc P25_JAC evaluates the three-element arch jacobian.
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
c    Input, integer OPTION, the option index.
c
c    Input, integer NVAR, the number of variables.
c
c    Input, double precision X(NVAR), the argument of the jacobian.
c
c    Output, double precision JAC(NVAR,NVAR), the jacobian matrix evaluated
c    at X.  The NVAR-th row is not set by this routine.
c
      implicit none

      integer nvar

      integer option
      double precision jac(nvar,nvar)
      double precision x(nvar)
      COMMON /AUXMEM/ THETA,C2,C3,C4,C5,C6,IFUN22,VALLM,VALNU

      do j = 1, ncol
        do i = 1, nrow
          jac(i,j) = 0.0E+00
        end do
      end do

      XLM=X(4)
      XNU=X(5)
      jac(1,1)=4.0*C2
     &           +8.0*THETA/15.0
     &           +4.0*C3*X(1)/5.0
     &           +6.0*C5*X(1)*X(1)/5.0
     &           +22.0*C5*X(2)*X(2)/35.0
      jac(1,2)=4.0*C2
     &           +28.0*C3*X(2)/105.0
     &           +8.0*C4*X(3)/15.0
     &           +44.0*C5*X(1)*X(2)/35.0
      jac(1,3)=8.0*C4*X(2)/15.0
      IF(IFUN22.EQ.1)
     & jac(1,4)=-0.5*C6*COS(THETA*XNU)*(XNU*XNU-1.0)
      IF(IFUN22.EQ.2)
     & jac(1,4)=-2.0*C6/3.0
      IF(IFUN22.EQ.1)
     & jac(1,5)=-C6*XLM*XNU*COS(THETA*XNU)
     &     +0.5*C6*XLM*THETA*(XNU*XNU-1.0)*SIN(THETA*XNU)
      IF(IFUN22.EQ.2)jac(1,5)=0.0
      jac(2,1)=4.0*C2
     &           +28.0*C3*X(2)/105.0
     &           +8.0*C4*X(3)/15.0
     &           +44.0*C5*X(1)*X(2)/35.0
      jac(2,2)=36.0*C2/5.0
     &           +8.0*THETA/105.0
     &           +28.0*C3*X(1)/105.0
     &           +22.0*C5*X(1)*X(1)/35.0
     &           +12.0*C5*X(2)*X(2)/35.0
      jac(2,3)=4.0/15.0
     &           +8.0*C4*X(1)/15.0
      IF(IFUN22.EQ.1)
     & jac(2,4)=-0.5*C6*COS(THETA*XNU)*(XNU*XNU-1.0)*XNU
      IF(IFUN22.EQ.2)
     & jac(2,4)=-C6/4.0
      IF(IFUN22.EQ.1)
     & jac(2,5)=-0.5*C6*XLM*(3.0*XNU*XNU-1.0)*COS(THETA*XNU)
     &      +0.5*C6*XLM*THETA*(XNU*XNU-1.0)*SIN(THETA*XNU)
      IF(IFUN22.EQ.2)jac(2,5)=0.0
      jac(3,1)=8.0*C4*X(2)/15.0
      jac(3,2)=4.0/15.0
     &           +8.0*C4*X(1)/15.0
      jac(3,3)=4.0/(3.0*THETA)
      IF(IFUN22.EQ.1)
     & jac(3,4)=-0.5*C6*SIN(THETA*XNU)*(XNU*XNU-1.0)
      IF(IFUN22.EQ.2)
     & jac(3,4)=0.0
      IF(IFUN22.EQ.1)
     & jac(3,5)=0.5*C6*XLM*2.0*XNU*SIN(THETA*XNU)
     &           +0.5*C6*XLM*(XNU*XNU-1.0)*COS(THETA*XNU)
      IF(IFUN22.EQ.2)jac(3,5)=0.0
      jac(4,5)=-1.0
      return
      END
      subroutine p25_nvar ( option, nvar )

c*********************************************************************72
c
cc P25_NVAR sets the number of variables for problem 25.
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
c    Input, integer OPTION, the option chosen for this problem.
c    For some problems, several options are available.  At least,
c    OPTION = 1 is always legal.
c
c    Output, integer NVAR, the number of variables.
c
      implicit none

      integer nvar
      integer option

      nvar = 5

      return
      end
      subroutine p25_option_num ( option_num )

c*********************************************************************72
c
cc P25_OPTION_NUM returns the number of options for problem 25.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    16 September 2008
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, integer OPTION_NUM, the number of options.
c
      implicit none

      integer option_num

      option_num = 2

      return
      end
      subroutine p25_start ( option, nvar, x )

c*********************************************************************72
c
cc P25_START returns a starting point for problem 25.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    19 September 2008
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer OPTION, the option index.
c
c    Input, integer NVAR, the number of variables.
c
c    Output, double precision X(NVAR), the starting point.
c
      implicit none

      integer nvar

      integer i
      integer option
      double precision x(nvar)

      if ( 1 .le. option .and. option .le. 2 ) then

        do i = 1, nvar
          x(i) = 0.0
        end do

      else

        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'P25_START - Fatal error!'
        write ( *, '(a,i8)' ) '  Unrecognized option  = ', option
        stop

      end if

      return
      end
      subroutine p25_stepsize ( option, h, hmin, hmax )

c*********************************************************************72
c
cc P25_STEPSIZE returns step sizes for problem 25.
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
c    Input, integer OPTION, the option index.
c
c    Output, double precision H, HMIN, HMAX, suggested values for the
c    initial step, the minimum step, and the maximum step.
c
      implicit none

      double precision h
      double precision hmax
      double precision hmin
      integer option

      h =    0.250D+00
      hmin = 0.001D+00
      hmax = 1.000D+00

      return
      end
      subroutine p25_title ( option, title )

c*********************************************************************72
c
cc P25_TITLE sets the title for problem 25.
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
c    Input, integer OPTION, the option index.
c
c    Output, character*(*) TITLE, the title of the problem.
c    TITLE will never be longer than 80 characters.
c
      implicit none

      integer option
      character*(*) title

      if ( option .eq. 1 ) then
        title = 'Three-element arch with point load.'
      else if ( option .eq. 2 ) then
        title = 'Three-element arch with constant load.'
      else
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'P25_TITLE - Fatal error!'
        write ( *, '(a)' ) '  Unrecognized option number.'
        stop
      end if

      return
      end
      subroutine p26_data (DIRIPC,HFACT,HMAX,HMIN,HTAN,option,
     & IPC,IPRAM,ISIZE,IT,LIM,MAXCON,
     & MAXLIM,MAXSTP,MAXTAR,NCOL,NROW,NVAR,RELERR,XIT,RWORK)

c*********************************************************************72
c
cc P26_DATA sets parameters for the circle-line function.
c
c  26.1 the circle-line problem
c
c  26.2 reference
c
c  26.3 the function
c
c  fx is of the form
c
c  fx(x,y) = y * (x*x + (y-y0)*(y-y0) - 1)
c
c  26.4 options
c
c  option=1  starting point is (1.22,0), y0=0. (check for bifurcation)
c  option=2  starting point is (1.22,0), y0=0.5. (check for bifurcation)
c
c  26.5 target points
c
c  26.6 limit points
c
c  26.7  bifurcation points
c
c  bifurcation points are at (1,0) and (-1,0).
c
c  26.8 comments
c
c  program seems ok.
c
      DIMENSION RWORK(ISIZE)
      COMMON /AUXMEM/ YZERO

      NVAR=2


      DIRIPC=-1.0
      HFACT=2.0
      IPC=1
      IT=0
      LIM=0
      MAXCON=10
      MAXLIM=10
      MAXSTP=15
      MAXTAR=10
      NCOL=2
      NROW=2
      RELERR=1.0E-5
      XIT=0.0

      return
      END
      subroutine p26_fun ( option, nvar, x, fx )

c*********************************************************************72
c
cc P26_FUN evaluates the circle-line function.
c
      integer nvar

      double precision fx(nvar-1)
      integer option
      double precision x(nvar)

      if ( option .eq. 1 ) then
        yzero = 0.0
      else if ( option .eq. 2 ) then
        yzero = 0.5
      end if

      FX(1)=X(2)*(X(1)*X(1)+(X(2)-YZERO)*(X(2)-YZERO)-1.0)

      return
      END
      subroutine p26_jac ( option, nvar, x, jac )

c*********************************************************************72
c
cc P26_JAC evaluates the circle-line jacobian.
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
c    Input, integer OPTION, the option index.
c
c    Input, integer NVAR, the number of variables.
c
c    Input, double precision X(NVAR), the argument of the jacobian.
c
c    Output, double precision JAC(NVAR,NVAR), the jacobian matrix evaluated
c    at X.  The NVAR-th row is not set by this routine.
c
      implicit none

      integer nvar

      integer option
      double precision jac(nvar,nvar)
      double precision x(nvar)

      if ( option .eq. 1 ) then
        yzero = 0.0
      else if ( option .eq. 2 ) then
        yzero = 0.5
      end if

      jac(1,1)=2.0*X(1)*X(2)
      jac(1,2)=X(1)*X(1)+3.0*X(2)*X(2)-4.0*X(2)*YZERO+YZERO*YZERO-1.0

      return
      END
      subroutine p26_nvar ( option, nvar )

c*********************************************************************72
c
cc P26_NVAR sets the number of variables for problem 26.
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
c    Input, integer OPTION, the option chosen for this problem.
c    For some problems, several options are available.  At least,
c    OPTION = 1 is always legal.
c
c    Output, integer NVAR, the number of variables.
c
      implicit none

      integer nvar
      integer option

      nvar = 2

      return
      end
      subroutine p26_option_num ( option_num )

c*********************************************************************72
c
cc P26_OPTION_NUM returns the number of options for problem 26.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    16 September 2008
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, integer OPTION_NUM, the number of options.
c
      implicit none

      integer option_num

      option_num = 2

      return
      end
      subroutine p26_start ( option, nvar, x )

c*********************************************************************72
c
cc P26_START returns a starting point for problem 26.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    19 September 2008
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer OPTION, the option index.
c
c    Input, integer NVAR, the number of variables.
c
c    Output, double precision X(NVAR), the starting point.
c
      implicit none

      integer nvar

      integer option
      double precision x(nvar)

      if ( 1 .le. option .and. option .le. 2 ) then

        x(1) = 1.22D+00
        x(2) = 0.0D+00

      else

        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'P26_START - Fatal error!'
        write ( *, '(a,i8)' ) '  Unrecognized option  = ', option
        stop

      end if

      return
      end
      subroutine p26_stepsize ( option, h, hmin, hmax )

c*********************************************************************72
c
cc P26_STEPSIZE returns step sizes for problem 26.
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
c    Input, integer OPTION, the option index.
c
c    Output, double precision H, HMIN, HMAX, suggested values for the
c    initial step, the minimum step, and the maximum step.
c
      implicit none

      double precision h
      double precision hmax
      double precision hmin
      integer option

      h =    0.100D+00
      hmin = 0.001D+00
      hmax = 0.100D+00

      return
      end
      subroutine p26_title ( option, title )

c*********************************************************************72
c
cc P26_TITLE sets the title for problem 26.
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
c    Input, integer OPTION, the option index.
c
c    Output, character*(*) TITLE, the title of the problem.
c    TITLE will never be longer than 80 characters.
c
      implicit none

      integer option
      character*(*) title

      if ( option .eq. 1 ) then
        title = 'Circle-line function, center (0.0, 0.0).'
      else if ( option .eq. 2 ) then
        title = 'Circle-line function, center (0.0, 0.5).'
      else
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'P26_TITLE - Fatal error!'
        write ( *, '(a)' ) '  Unrecognized option number.'
        stop
      end if

      return
      end
      subroutine p27_data (DIRIPC,HFACT,HMAX,HMIN,HTAN,option,
     & IPC,IPRAM,ISIZE,IT,LIM,MAXCON,
     & MAXLIM,MAXSTP,MAXTAR,NCOL,NROW,NVAR,RELERR,XIT,XR)

c*********************************************************************72
c
cc P27_DATA sets parameters for finite difference 2D Poisson problem.
c
c  27.1  2-d finite difference  del2 u + 10*l*sin(u) + delta*cos(u)=0
c
c  27.2  reference
c
c  professor shi miaogen
c  work in progress
c
c  27.3  the problem
c
c  del2 u + 10*lambda*sin(u) + delta*cos(u) = 0.0
c  u=0 on boundary.
c
c  region is (0,1) x (0,1).
c
c  region is divided into grid lines, u is evaluated at mesh
c  points.  this version has a 7 by 7 grid, 25 interior mesh points.
c
c
c
c  27.4  options
c
c  option=1  delta=0.0  (bifurcations from zero solution)
c  option=2  delta=0.5  (non zero solution)
c
c  27.5  target points
c
c  27.6  limit points
c
c  27.7  bifurcation points
c
c  numerous bifurcation points of various orders split off
c  from the solution u=0.
c
c  using the formula
c
c  10*lambda=4*n*n*(1-0.5*(cos(kh/n)+cos(l*h/n)))
c  n=6, h=1/6, 1.le.k,l.le.n-1.
c  as calculated, they are
c
c  order  lambda
c
c    1     1.9292342
c    2     4.5646171
c    1     7.2
c    2     8.1646171
c    2    10.8
c    2    11.764617
c    5    14.4
c    2    17.037029
c    2    18.0
c    2    20.635383
c    1    21.6
c    2    24.235383
c    1    26.870766
c
c  27.8  comments
c
      COMMON /PERTUR/ DELTA
      DIMENSION XR(26)
      NVAR=26
      DELTA=0.0
      IF(option.EQ.2)DELTA=0.5

      DIRIPC=1.0
      HFACT=3.0
      IF(IPRAM.EQ.3)HFACT=1.5
      IPC=NVAR
      IT=NVAR
      MAXSTP=10
      MAXCON=MAXSTP
      MAXLIM=MAXSTP
      MAXTAR=1
      NCOL=NVAR
      NROW=NVAR
      RELERR=1.0E-6
      XIT=10.0

      return
      END
      subroutine p27_fun ( option, nvar, x, fx )

c*********************************************************************72
c
cc P27_FUN evaluates the finite difference 2D Poisson function.
c
      integer nvar

      double precision fx(nvar-1)
      integer option
      double precision x(nvar)

      COMMON /PERTUR/ DELTA

      NSIDE=IFIX(SQRT(FLOAT(NVAR-1)))
      DX=1.0/(FLOAT(NSIDE+1))
      IEQN=0

      DO I=1,NSIDE
        DO J=1,NSIDE
          IEQN=IEQN+1
          IN=IEQN-NSIDE
          UN=0.0
          IF(I.GT.1.AND.IN.GT.0.AND.IN.LT.NVAR)UN=X(IN)
          IW=IEQN-1
          UW=0.0
          IF(J.GT.1.AND.IW.GT.0.AND.IW.LT.NVAR)UW=X(IW)
          IE=IEQN+1
          UE=0.0
          IF(J.LT.NSIDE.AND.IE.GT.0.AND.IE.LT.NVAR)UE=X(IE)
          IS=IEQN+NSIDE
          US=0.0
          IF(I.LT.NSIDE.AND.IS.GT.0.AND.IS.LT.NVAR)US=X(IS)
          FX(IEQN)=-4.0*X(IEQN)+(UN+US+UE+UW)
          FX(IEQN)=FX(IEQN)/(DX*DX)
          FX(IEQN)=FX(IEQN)+10.0*X(NVAR)*SIN(X(IEQN))+DELTA*COS(X(IEQN))
        end do
      end do

      return
      END
      subroutine p27_jac ( option, nvar, x, jac )

c*********************************************************************72
c
cc P27_JAC evaluates the finite difference 2D Poisson jacobian.
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
c    Input, integer OPTION, the option index.
c
c    Input, integer NVAR, the number of variables.
c
c    Input, double precision X(NVAR), the argument of the jacobian.
c
c    Output, double precision JAC(NVAR,NVAR), the jacobian matrix evaluated
c    at X.  The NVAR-th row is not set by this routine.
c
      implicit none

      integer nvar

      integer option
      double precision jac(nvar,nvar)
      double precision x(nvar)

      COMMON /PERTUR/ DELTA

      do j = 1, ncol
        do i = 1, nrow
          jac(i,j) = 0.0E+00
        end do
      end do

      NSIDE=IFIX(SQRT(FLOAT(NVAR-1)))
      DX=1.0/(FLOAT(NSIDE+1))
      IEQN=0

      DO I=1,NSIDE
        DO J=1,NSIDE
          IEQN=IEQN+1
          IN=IEQN-NSIDE
          IW=IEQN-1
          IE=IEQN+1
          IS=IEQN+NSIDE
          jac(IEQN,IEQN)=-4.0/(DX*DX)+10.0*X(NVAR)*COS(X(IEQN))
     &    -DELTA*SIN(X(IEQN))
          jac(IEQN,NVAR)=10.0*SIN(X(IEQN))
          IF(I.GT.1.AND.IN.GT.0.AND.IN.LT.NVAR)
     &    jac(IEQN,IN)=1.0/(DX*DX)
          IF(J.GT.1.AND.IW.GT.0.AND.IW.LT.NVAR)
     &    jac(IEQN,IW)=1.0/(DX*DX)
          IF(J.LT.NSIDE.AND.IE.GT.0.AND.IE.LT.NVAR)
     &    jac(IEQN,IE)=1.0/(DX*DX)
          IF(I.LT.NSIDE.AND.IS.GT.0.AND.IS.LT.NVAR)
     &    jac(IEQN,IS)=1.0/(DX*DX)
        end do
      end do

      return
      END
      subroutine p27_nvar ( option, nvar )

c*********************************************************************72
c
cc P27_NVAR sets the number of variables for problem 27.
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
c    Input, integer OPTION, the option chosen for this problem.
c    For some problems, several options are available.  At least,
c    OPTION = 1 is always legal.
c
c    Output, integer NVAR, the number of variables.
c
      implicit none

      integer nvar
      integer option

      nvar = 26

      return
      end
      subroutine p27_option_num ( option_num )

c*********************************************************************72
c
cc P27_OPTION_NUM returns the number of options for problem 27.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    16 September 2008
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, integer OPTION_NUM, the number of options.
c
      implicit none

      integer option_num

      option_num = 2

      return
      end
      subroutine p27_start ( option, nvar, x )

c*********************************************************************72
c
cc P27_START returns a starting point for problem 27.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    19 September 2008
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer OPTION, the option index.
c
c    Input, integer NVAR, the number of variables.
c
c    Output, double precision X(NVAR), the starting point.
c
      implicit none

      integer nvar

      integer i
      integer option
      double precision x(nvar)

      if ( 1 .le. option .and. option .le. 2 ) then

        do i = 1, nvar
          x(i) = 0.0
        end do

      else

        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'P27_START - Fatal error!'
        write ( *, '(a,i8)' ) '  Unrecognized option  = ', option
        stop

      end if

      return
      end
      subroutine p27_stepsize ( option, h, hmin, hmax )

c*********************************************************************72
c
cc P27_STEPSIZE returns step sizes for problem 27.
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
c    Input, integer OPTION, the option index.
c
c    Output, double precision H, HMIN, HMAX, suggested values for the
c    initial step, the minimum step, and the maximum step.
c
      implicit none

      double precision h
      double precision hmax
      double precision hmin
      integer option

      h =    1.00000D+00
      hmin = 0.03125D+00
      hmax = 4.00000D+00

      return
      end
      subroutine p27_title ( option, title )

c*********************************************************************72
c
cc P27_TITLE sets the title for problem 27.
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
c    Input, integer OPTION, the option index.
c
c    Output, character*(*) TITLE, the title of the problem.
c    TITLE will never be longer than 80 characters.
c
      implicit none

      integer option
      character*(*) title

      if ( option .eq. 1 ) then
        title = '2D finite difference Poisson function, delta = 0.'
      else if ( option .eq. 2 ) then
        title = '2D finite difference Poisson function, delta = 0.5.'
      else
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'P27_TITLE - Fatal error!'
        write ( *, '(a)' ) '  Unrecognized option number.'
        stop
      end if

      return
      end
      subroutine p28_data (DIRIPC,HFACT,HMAX,HMIN,HTAN,option,
     & IPC,IPRAM,ISIZE,IT,LIM,MAXCON,
     & MAXLIM,MAXSTP,MAXTAR,NCOL,NROW,NVAR,RELERR,XIT,RWORK)

c*********************************************************************72
c
cc P28_DATA sets parameters for Walker's arch.
c
c  walker's arch.
c
c  Reference
c
c  a c walker,
c  a nonlinear finite element analysis of shallow circular arches,
c  international journal of solids and structures 5, 1969, pages 97-107.
c
c  The function
c
c  the function models the deformations of a thin, shallow, circular
c  arch of radius 10, thickness .0684, and subtending 30 degrees.
c  8 intervals are used.  shape functions introduced by
c  a c walker are employed.
c  two control variables, u=x(46) and v=x(47) are added which
c  determine a value of the load function.  an auxilliary
c  equation x(ifix1)=val1 is added, where ifix1=46 or ifix1=47.
c
c  Options
c
c  option=1  ifix1=46, val1=0.1,  load function 1, lim=47.
c  option=2  ifix1=47, val1=0.01, load function 1, lim=46.
c  option=3  ifix1=46, val1=0.0,  load function 2, lim=47.
c  option=4  ifix1=46, val1=1.0,  load function 3, lim=47.
c  option=5  ifix1=47, val1=0.1,  load function 3, lim=46.
c  option=6  ifix1=46, val1=0.0,  load function 4, lim=47.
c  option=7  ifix1=47, val1=0.0, load function 5, lim=46
c  option=8  ifix1=46, val1=input value of rwork(46), load function 5,
c            lim=47, it is assumed that option=7 has just
c            been run, and that this is a continuation of that run.
c  option=9  information set from terminal (ie. timesharing)
c  option=10 bifurcation run.  ifix1=46, val1=0.0,
c            load function 2, lim=0.
c
c  load function 1 - load=u   for intervals 1,2,3,5,6,7,8
c                    load=u+v for interval 4.
c  load function 2 - load=v.
c  load function 3 - load= (1.0-v)*u for intervals 1 thru 4.
c                    load= (1.0+v)*u for intervals 5 thru 8.
c  load function 4 - load= 15*v/8 for interval 2
c                    load= 7*v/8 for intervals 1,3,4,5,6,7,8.
c  load function 5 - load=u-u*(v-x)/wide for v-wide.le.x.le.v
c                        =u+u*(v-x)/wide for v.le.x.le.v+wide
c                        =0 elsewhere
c                    wide=theta/(2*intrvl)=1/4 interval
c
c  Target points
c
c  the target point has x(ifree)=1.7
c
c  option=3, for u=0.0, v=1.7, center deflection=.340011
c
c  Limit points
c
c  melhem lists the following limit points
c
c  x(21) represents the radial displacement of the center node.
c
c  for option=3, x(21)=0.8493987, v=x(ifree)=2.270319, u=x(ifix1)=0.0.
c  for option=5, x(21)=0.7430691, u=x(ifree)=1.637480, v=x(ifix1)=0.1
c  for option=6, x(21)=0.7694139, v=x(ifree)=1.662499, u=x(ifix1)=0.0.
c  for option=7, x(21)=1.81346,   u=x(ifree)=.657347,  v=x(ifix1)=0.0
c
c  computational results with this program
c  limit points in u or v.
c
c  option=3, x(21)=0.8419974, u=0.0, v=2.255799.
c  option=4, x(21)=3.129776,  u=1.0, v=1.675206.
c  option=6, x(21)=0.7467170, u=0.0, v=1.658697.
c
c  Bifurcation points
c
c  for option option=10, there is a bifurcation point.
c
c  Comments
c
c  compare run of option=5 with stated melhem result.
c  note error in melhem/rheinboldt paper for their third
c    load function.
c  note that for option=2, solution will not go beyond u=pi.
c  must verify data.
c  get bifurcation points.
c
      real pi
      parameter ( pi = 3.1415926535E+00 )
      DIMENSION RWORK(ISIZE)

      COMMON /AUXMEM/ SINT(6),COEF(6),DL,ALFA0,THETA,ALFA1,IFUN17,
     & IFREE,IFIX1,VAL1,INTRVL,WIDE,NGAUSS,RADIUS

      NVAR=47
      NEQN=NVAR-1
      NCOF=NVAR-2
c
c  set number of gauss points, gauss points, gauss coefficients
c
      NGAUSS=6
      SINT(1)=.96623476
      SINT(2)=.8306047
      SINT(3)=.6193096
      SINT(4)=.38069040
      SINT(5)=.1693953
      SINT(6)=.033765332
      COEF(1)=.085662246
      COEF(2)=.18038078
      COEF(3)=.23395696
      COEF(4)=.23395696
      COEF(5)=.18038078
      COEF(6)=.085662246
c
c  set problem parameters
c
c  radius=10
c  thick=0.0684
c  ei=796
c  ea=2.056e+06
c  theta=15.*pi/180
c  alfa0=thick/radius
c  alfa1=ei/(ea*radius**2).
c  dl=ei/(ea*thick*radius)
c
      RADIUS=10.0
      THICK=0.0684
      THETA=PI/12.0
      EI=796.0
      EA=2.056E+06
      ALFA0=THICK/RADIUS
      ALFA1=EI/(EA*RADIUS*RADIUS)
      DL=EI/(EA*THICK*RADIUS)
      IF(option.EQ.1)GO TO 10
      IF(option.EQ.2)GO TO 20
      IF(option.EQ.3)GO TO 30
      IF(option.EQ.4)GO TO 40
      IF(option.EQ.5)GO TO 50
      IF(option.EQ.6)GO TO 60
      IF(option.EQ.7)GO TO 70
      IF(option.EQ.8)GO TO 80
      IF(option.EQ.9)GO TO 90
      IF(option.EQ.10)GO TO 100
   10 CONTINUE
      IFUN17=1
      IFIX1=NEQN
      IFREE=NVAR
      VALU=0.1
      VALV=0.0
      VAL1=0.1
      GO TO 110
   20 CONTINUE
      IFUN17=1
      IFIX1=NVAR
      IFREE=NEQN
      VALU=0.0
      VALV=0.01
      VAL1=0.01
      GO TO 110
   30 CONTINUE
      IFUN17=2
      IFIX1=NEQN
      IFREE=NVAR
      VALU=0.0
      VALV=0.0
      VAL1=0.0
      GO TO 110
   40 CONTINUE
      IFUN17=3
      IFIX1=NEQN
      IFREE=NVAR
      VALU=1.0
      VALV=0.0
      VAL1=1.0
      GO TO 110
   50 CONTINUE
      IFUN17=3
      IFIX1=NVAR
      IFREE=NEQN
      VALU=0.0
      VALV=0.1
      VAL1=0.1
      GO TO 110
   60 CONTINUE
      IFUN17=4
      IFIX1=NEQN
      IFREE=NVAR
      VALU=0.0
      VALV=0.0
      VAL1=0.0
      GO TO 110
   70 CONTINUE
      IFUN17=5
      IFIX1=NVAR
      IFREE=NEQN
      VALU=0.0
      VALV=0.0
      VAL1=0.0
      GO TO 110
   80 CONTINUE
      IFUN17=5
      IFIX1=NEQN
      IFREE=NVAR
      VALU=RWORK(NEQN)
      VALV=0.0
      VAL1=VALU
      GO TO 110
c
c  determine load function desired
c
   90 CONTINUE
      WRITE(*,1010)
      call READIN(IFUN17)
      IF(IFUN17.LT.1.OR.IFUN17.GT.6)IFUN17=1
c
c  find whether fixing u or v
c
      WRITE(*,1020)NEQN,NVAR
      call READIN(IFIX1)
      IF(IFIX1.NE.NEQN)IFIX1=NVAR
      IFREE=NVAR
      IF(IFIX1.EQ.NVAR)IFREE=NEQN
c
c  get initial values for u and v
c
      WRITE(6,1030)
      call READRL(VALU)
      WRITE(6,1040)
      call READRL(VALV)
      IF(IFIX1.EQ.NEQN)VAL1=VALU
      IF(IFIX1.EQ.NVAR)VAL1=VALV
      GO TO 110
  100 CONTINUE
      IFUN17=2
      IFIX1=NEQN
      IFREE=NVAR
      VALU=0.0
      VALV=0.0
      VAL1=0.0
      GO TO 110
  110 CONTINUE
c
c  set number of intervals
c
      INTRVL=8
c
c  set width for load function 5
c
      WIDE=2.0*THETA/FLOAT(INTRVL)


      DIRIPC=1.0
      HFACT=3.0
      IPC=IFREE
      IT=IFREE
      LIM=IFREE
      IF(option.EQ.10)LIM=0
      MAXCON=45
      MAXLIM=1
      MAXSTP=45
      MAXTAR=45
      NCOL=NVAR
      NROW=NVAR
      RELERR=1.0E-4
      XIT=1.7
c
c  set starting point
c
      DO I=1,NCOF
        RWORK(I)=0.0
      end do
      RWORK(NEQN)=VALU
      RWORK(NVAR)=VALV
      return
 1010 FORMAT(36H ENTER LOAD FUNCTION BETWEEN 1 AND 6)
 1020 FORMAT(6H TYPE ,I4,11H TO FIX U, ,I4,9H TO FIX V)
 1030 FORMAT(26H ENTER INITIAL VALUE FOR U)
 1040 FORMAT(26H ENTER INITIAL VALUE FOR V)
      END
      subroutine p28_edit (ARCLXC,ARCLXF,ARCLXR,DIRIPC,FPNAME,
     & jac,FXNAME,HFACT,HMAX,HMIN,HTANCF,IERROR,IJAC,IPC,IPIVOT,
     & IPL,IPOINT,IPRAM,ISTART,IT,IWRITE,JPOINT,JWRITE,KWRITE,LIM,
     & LWRITE,MODNEW,NCOL,NROW,NUMCON,NUMLIM,NUMSTR,NUMTAR,
     & NVAR,RELERR,SLNAME,TC,TL,TR,WORK1,WORK2,WORK3,XC,XF,XIT,XR)

c*********************************************************************72
c
cc P28_EDIT is an edit routine for Walker's arch.
c
      EXTERNAL FPNAME
      EXTERNAL FXNAME
      EXTERNAL SLNAME
      DIMENSION jac(NROW,NCOL)
      DIMENSION IPIVOT(NVAR)
      DIMENSION TC(NVAR)
      DIMENSION TL(NVAR)
      DIMENSION TR(NVAR)
      DIMENSION WORK1(NVAR)
      DIMENSION WORK2(NVAR)
      DIMENSION WORK3(NVAR)
      DIMENSION XC(NVAR)
      DIMENSION XF(NVAR)
      DIMENSION XR(NVAR)
      DIMENSION XVALUE(33),YVALUE(33)
      DIMENSION PH(6),PHP(6),PHPP(6),PS(4),PSP(4)
      COMMON /AUXMEM/ SINT(6),COEF(6),DL,ALFA0,THETA,ALFA1,IFUN17,
     & IFREE,IFIX1,VAL1,INTRVL,WIDE,NGAUSS,RADIUS
      DATA ILIST  /0/
      NEQN=NVAR-1
      DEF=XR(21)
      WRITE(*,1010)DEF,IFUN17
      IF(IFUN17.EQ.5)GO TO 10
      VALU=XR(NEQN)
      VALV=XR(NVAR)
      WRITE(*,1020)VALU,VALV
      GO TO 20
   10 HEIGHT=XR(NEQN)
      CENTER=XR(NVAR)
      WRITE(*,1030)HEIGHT,WIDE,CENTER
   20 CONTINUE
c     if(ipoint.ne.5)return
c
c  if (ilist.eq.1) print out solution
c
      IF(ILIST.EQ.1)WRITE(*,*)' '
      IF(ILIST.EQ.1)WRITE(*,1050)
      IF(ILIST.EQ.1)WRITE(*,*)' '
      NPLOT=33
      ANG0=-THETA
      DANG0=2.0*THETA/FLOAT(NPLOT-1)
      IPLOT=0
      XDELT=2.0*THETA/FLOAT(NPLOT-1)
      HSPACE=2.0*THETA/FLOAT(INTRVL)
      XVAL=-THETA
      IF(ILIST.EQ.1)WRITE(*,1000)
      DO I=1,INTRVL
        INT=I
        S=0.0
        DO J=1,4
          call p28_bs (S,PH,PHP,PHPP,PS,PSP)
          call p28_vl (INT,XR,NVAR,PH,PHP,PHPP,PS,PSP,OMGA,OMGAP,OMGAPP,
     &    UX,UP)
          IF(ILIST.EQ.1)WRITE(*,1040)INT,OMGA,UX
          IPLOT=IPLOT+1
          RAD=RADIUS-OMGA
          DTH=ATAN(-UX/RAD)
          ANG=ANG0+DTH
          ANG0=ANG0+DANG0
          XVALUE(IPLOT)=RAD*COS(ANG)
          YVALUE(IPLOT)=RAD*SIN(ANG)
          S=S+0.25
        end do
      end do

      INT=INTRVL
      S=1.0
      call p28_bs (S,PH,PHP,PHPP,PS,PSP)
      call p28_vl (INT,XR,NVAR,PH,PHP,PHPP,PS,PSP,OMGA,OMGAP,OMGAPP,
     & UX,UP)
      IF(ILIST.EQ.1)WRITE(*,1040)INT,OMGA,UX
      IPLOT=IPLOT+1
      RAD=RADIUS-OMGA
      DTH=ATAN(-UX/RAD)
      ANG=ANG0+DTH
      ANG0=ANG0+DANG0
      XVALUE(IPLOT)=RAD*COS(ANG)
      YVALUE(IPLOT)=RAD*SIN(ANG)
c
c  call dotplt
c
      IYMAX=16
      IXMAX=61
      NVAL=IPLOT
      call DOTPLT(IXMAX,IYMAX,NVAL,YVALUE,XVALUE)
      return
 1000 FORMAT(' INTERVAL  RADIAL DISPLACEMENT  TANGENTIAL DISPLACEMENT')
 1010 FORMAT(19H CENTER DEFLECTION=,G14.6,21H LOAD FUNCTION NUMBER,I4)
 1020 FORMAT(3H U=,G14.6,3H V=,G14.6)
 1030 FORMAT(13H LOAD HEIGHT=,G14.6,7H WIDTH=,G14.6,8H CENTER=,G14.6)
 1040 FORMAT(3X,I6,7X,G14.6,11X,G14.6)
 1050 FORMAT(13H PLOT OF ARCH)
      END
      subroutine p28_fun ( option, nvar, x, fx )

c*********************************************************************72
c
cc P28_FUN evaluates Walker's arch function.
c
      integer nvar

      double precision fx(nvar-1)
      integer option
      double precision x(nvar)

      DIMENSION PH(6)
      DIMENSION PHP(6)
      DIMENSION PHPP(6)
      DIMENSION PS(4)
      DIMENSION PSP(4)
      DIMENSION TEMP(5)

      COMMON /AUXMEM/ SINT(6),COEF(6),DL,ALFA0,THETA,ALFA1,IFUN17,
     & IFREE,IFIX1,VAL1,INTRVL,WIDE,NGAUSS,RADIUS
      NCOF=NVAR-2
      NEQN=NVAR-1
      NCOFM1=NCOF-1
      NCOFM2=NCOF-2
      NCOFM3=NCOF-3
      NCOFM4=NCOF-4
      BT=FLOAT(INTRVL)/(2.0*THETA)
      HSPACE=2.0*THETA/FLOAT(INTRVL)

      DO I=1,NCOF
        FX(I)=0.0
      end do
c
c  loop for summation in gauss quadrature.
c
      DO 70 K=1,NGAUSS

        call p28_bs (SINT(K),PH,PHP,PHPP,PS,PSP)

        DO IK=1,5
          TEMP(IK)=0.0
        end do
c
c  loop on all the intervals.
c
        DO L=1,INTRVL
c
c  determine actual position xpos in interval -theta,theta
c  and evaluate load function
c
          XPOS=-THETA+HSPACE*(FLOAT(L-1)+SINT(K))

          call p28_gx (XPOS,X,NVAR,DLOAD)
c
c  evaluate omga, omgap, omgapp and up for this interval
c
          INT=L
          call p28_vl (INT,X,NVAR,PH,PHP,PHPP,PS,PSP,OMGA,OMGAP,OMGAPP,
     &    UX,UP)
c
c  for each interval find 10 equations (integrals), the last five
c  are stored in temp , to be added to the first five of the next
c  interval. a map is needed because in each interval
c  x(1) ...x(5) are  u(1),u(2),u(3),w(1),w(2).
c  x(6) ...x(10) are u(4),u(5),u(6),w(3),w(4).
c
          DO 60 IM=1,10
            IF((IM.LE.3).OR.((IM.GE.6).AND.(IM.LE.8))) GO TO 30
            IF(IM.LE.5)KM=IM-3
            IF(IM.GT.5) KM=IM-6
            FTEMP=(ALFA0*(BT*UP-OMGA)+(BT**2*ALFA0**2*OMGAP**2)/2.0)
     &          *PSP(KM)
            GO TO 40
   30       IF(IM.LE.3)KM=IM
            IF(IM.GT.3) KM=IM-2
            FTEMP=((((BT*UP-OMGA)+(BT**2*ALFA0*OMGAP**2)/2.0)*
     &      (BT**2*OMGAP*ALFA0*PHP(KM)-PH(KM))+
     &      (ALFA1*BT**4*OMGAPP*PHPP(KM)))*ALFA0-(DL*DLOAD*PH(KM)))/BT
   40       IF(IM.LE.5) GO TO 50
            IMM5=IM-5
            TEMP(IMM5)=FTEMP
            GO TO 60
   50       I=(L-1)*5+IM
            FX(I)=FX(I)+COEF(K)*(FTEMP+TEMP(IM))
   60       CONTINUE

        end do
c
c  adjust equations of the last intervals.
c
        FX(NCOFM2)=FX(NCOFM2)+COEF(K)*TEMP(3)
        FX(NCOF)=FX(NCOF)+COEF(K)*TEMP(5)
   70   CONTINUE
c
c  adjust equations corresponding to the initial conditions.
c
      FX(1)=X(1)
      FX(2)=X(2)
      FX(4)=X(4)
      FX(NCOFM1)=X(NCOFM1)
      FX(NCOFM4)=X(NCOFM4)
      FX(NCOFM3)=X(NCOFM3)
      FX(NEQN)=X(IFIX1)-VAL1
      return
      END
      subroutine p28_bs (S,PH,PHP,PHPP,PS,PSP)

c*********************************************************************72
c
cc P28_BS is an auxilliary routine for Walker's arch.
c
      DIMENSION PH(6),PHP(6),PHPP(6),PS(4),PSP(4)
      S2=S*S
      S3=S*S2
      PH(1)=((-6.*S+15.)*S-10.)*S3+1.
      PH(2)=((-3.*S+8.)*S-6.)*S3+S
      PH(3)=(((-0.5*S+1.5)*S-1.5)*S+0.5)*S2
      PH(4)=((6.*S-15.)*S+10.)*S3
      PH(5)=((-3.*S+7.)*S-4.)*S3
      PH(6)=((0.5*S-1.)*S+0.5)*S3
      PHP(1)=((-30.*S+60.)*S-30.)*S2
      PHP(2)=((-15.*S+32.0)*S-18.)*S2+1.
      PHP(3)=(((-2.5*S+6.)*S-4.5)*S+1.)*S
      PHP(4)=((30.*S-60.)*S+30.)*S2
      PHP(5)=((-15.*S+28.)*S-12.)*S2
      PHP(6)=((2.5*S-4.)*S+1.5)*S2
      PHPP(1)=((-120.*S+180.)*S-60.)*S
      PHPP(2)=((-60.*S+96.)*S-36.)*S
      PHPP(3)=((-10.*S+18.)*S-9.)*S+1.
      PHPP(4)=((120.*S-180.)*S+60.)*S
      PHPP(5)=((-60.*S+84.0)*S-24.)*S
      PHPP(6)=((10.*S-12.)*S+3.)*S
      PS(1)=1.0-S*S*(3.0-2.0*S)
      PS(2)=S*(1.0-S)**2
      PS(3)=S*S*(3.0-2.0*S)
      PS(4)=S*S*(S-1.0)
      PSP(1)=6.*(S2-S)
      PSP(2)=(3.0*S-4.0)*S+1.0
      PSP(3)=-PSP(1)
      PSP(4)=(3.*S-2.)*S
      return
      END
      subroutine p28_gx (XPOS,X,NVAR,DLOAD)

c*********************************************************************72
c
cc P28_GX is an auxilliary function for Walker's arch.
c
      DIMENSION X(NVAR)
      COMMON /AUXMEM/ SINT(6),COEF(6),DL,ALFA0,THETA,ALFA1,IFUN17,
     & IFREE,IFIX1,VAL1,INTRVL,WIDE,NGAUSS,RADIUS
      NEQN=NVAR-1
      IF(IFUN17.EQ.1)GO TO 10
      IF(IFUN17.EQ.2)GO TO 20
      IF(IFUN17.EQ.3)GO TO 30
      IF(IFUN17.EQ.4)GO TO 40
      IF(IFUN17.EQ.5)GO TO 50
      IF(IFUN17.EQ.6)GO TO 60
   10 X1=-0.25*THETA
      X2=0.0
      IF(XPOS.LT.X1.OR.XPOS.GT.X2)DLOAD=X(NEQN)
      IF(XPOS.GE.X1.AND.XPOS.LE.X2)DLOAD=X(NEQN)+X(NVAR)
      return
   20 DLOAD=X(NVAR)
      return
   30 IF(XPOS.LE.0.0)DLOAD=(1.0-X(NVAR))*X(NEQN)
      IF(XPOS.GT.0.0)DLOAD=(1.0+X(NVAR))*X(NEQN)
      return
   40 X1=-0.75*THETA
      X2=-0.50*THETA
      IF(XPOS.GT.X1.AND.XPOS.LT.X2)DLOAD=15.0*X(NVAR)/8.0
      IF(XPOS.LE.X1.OR.XPOS.GE.X2)DLOAD=7.0*X(NVAR)/8.0
      return
   50 HEIGHT=X(NEQN)
      CENTER=X(NVAR)
      DLOAD=0.0
      IF(XPOS.GT.CENTER.AND.XPOS.LE.(CENTER+WIDE))
     & DLOAD=HEIGHT*(1.0+(CENTER-XPOS)/WIDE)
      IF(XPOS.LE.CENTER.AND.XPOS.GE.(CENTER-WIDE))
     & DLOAD=HEIGHT*(1.0-(CENTER-XPOS)/WIDE)
      DLOAD=2.0*THETA*DLOAD/WIDE
      return
   60 IF(XPOS.LT.X(NVAR))DLOAD=0.0
      IF(XPOS.GE.X(NVAR))DLOAD=X(NEQN)
      return
      END
      subroutine p28_jac ( option, nvar, x, jac )

c*********************************************************************72
c
cc P28_JAC evaluates the Walker's arch jacobian.
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
c    Input, integer OPTION, the option index.
c
c    Input, integer NVAR, the number of variables.
c
c    Input, double precision X(NVAR), the argument of the jacobian.
c
c    Output, double precision JAC(NVAR,NVAR), the jacobian matrix evaluated
c    at X.  The NVAR-th row is not set by this routine.
c
      implicit none

      integer nvar

      integer option
      double precision jac(nvar,nvar)
      double precision x(nvar)

      DIMENSION PH(6),PHP(6),PHPP(6),
     & PS(4),PSP(4),TEMP(5,5),STOR(3),GTOR(3)
      COMMON /AUXMEM/ SINT(6),COEF(6),DL,ALFA0,THETA,ALFA1,IFUN17,
     & IFREE,IFIX1,VAL1,INTRVL,WIDE,NGAUSS,RADIUS
      NCOF=NVAR-2
      NEQN=NVAR-1
      NCOFM1=NCOF-1
      NCOFM2=NCOF-2
      NCOFM3=NCOF-3
      NCOFM4=NCOF-4
      NCOFM5=NCOF-5
      NCOFM9=NCOF-9
      BT=FLOAT(INTRVL)/(2.0*THETA)
      HSPACE=2.0*THETA/FLOAT(INTRVL)

      do j = 1, ncol
        do i = 1, nrow
          jac(i,j) = 0.0E+00
        end do
      end do
c
c  loop for summing the gauss quadrature.
c
      DO 200 K=1,NGAUSS

        call p28_bs (SINT(K),PH,PHP,PHPP,PS,PSP)

        DO I=1,5
          DO J=1,5
            TEMP(I,J)=0.0
          end do
        end do

        DO I=1,3
          GTOR(I)=0.0
          STOR(I)=0.0
        end do
c
c  loop over all intervals.
c
        DO 190 L=1,INTRVL
c
c  determine actual position xpos in interval -theta,theta
c  and evaluate load function derivatives
c
          XPOS=-THETA+HSPACE*(FLOAT(L-1)+SINT(K))
          call p28_gp (XPOS,X,NVAR,DPDU,DPDV)
c
c  get values of omga, omgap, omgapp, and up
c
          INT=L
          call p28_vl (INT,X,NVAR,PH,PHP,PHPP,PS,PSP,OMGA,OMGAP,OMGAPP,
     &    UX,UP)
c
c  for each interval we have 100=10*10 derivatives, only 55 of
c  which are computed because of the symmetry.
c  the derivatives in the lower right corner of the 10*10 matrix
c  are stored in temp to be added to the derivatives in the
c  upper left corner of the next interval.
c
          IT=(L-1)*5
          DO 70 JM=1,5
            DO 70 KM=1,JM
              IF(KM.GE.4) GO TO 50
              IF(JM.GE.4) GO TO 40
              JMM=JM
              KMM=KM
              FTEMP=(((-PH(KMM)+BT**2*ALFA0*OMGAP*PHP(KMM))*
     &        (-PH(JMM)+BT**2*ALFA0*OMGAP*PHP(JMM))/BT)+ALFA0*
     &        BT*PHP(KMM)*PHP(JMM)*((BT*UP-OMGA)+(BT**2*ALFA0*OMGAP**2)
     &        *.5)+ALFA1*BT**3*PHPP(KMM)*PHPP(JMM))*ALFA0
              GO TO 60
   40         JMM=JM-3
              KMM=KM
              FTEMP=(PSP(JMM)*(-PH(KMM)+BT**2*ALFA0*OMGAP*PHP(KMM)))
              FTEMP=ALFA0*FTEMP
              GO TO 60
   50         KMM=KM-3
              JMM=JM-3
              FTEMP=BT*PSP(KMM)*PSP(JMM)*ALFA0
   60         ITPJM=IT+JM
              ITPKM=IT+KM
              jac(ITPJM,ITPKM)=jac(ITPJM,ITPKM)
     &        +COEF(K)*(FTEMP+TEMP(JM,KM))
   70         CONTINUE

          DO JM=6,10
            DO KM=1,5
              IF(KM.GT.3) GO TO 90
              KMM=KM
              IF(JM.LE.8) GO TO 80
              JMM=JM-6
              FTEMP=(PSP(JMM)*(-PH(KMM)+BT**2*ALFA0*OMGAP*PHP(KMM)))
              FTEMP=ALFA0*FTEMP
              GO TO 110
   80         JMM=JM-2
              FTEMP=(((-PH(KMM)+BT**2*ALFA0*OMGAP*PHP(KMM))*
     &        (-PH(JMM)+BT**2*ALFA0*OMGAP*PHP(JMM))/BT)+ALFA0*
     &        BT*PHP(KMM)*PHP(JMM)*((BT*UP-OMGA)+(BT**2*ALFA0*OMGAP**2)
     &        *.5)+ALFA1*BT**3*PHPP(KMM)*PHPP(JMM))*ALFA0
              GO TO 110
   90         KMM=KM-3
              IF(JM.LE.8) GO TO 100
              JMM=JM-6
              FTEMP=BT*PSP(KMM)*PSP(JMM)*ALFA0
              GO TO 110
  100         JMM=JM-2
              FTEMP=(PSP(KMM)*(-PH(JMM)+BT**2*ALFA0*OMGAP*PHP(JMM)))
              FTEMP=FTEMP*ALFA0
  110         ITPJM=IT+JM
              ITPKM=IT+KM
              jac(ITPJM,ITPKM)=jac(ITPJM,ITPKM)+COEF(K)*FTEMP
            end do
          end do

          DO 160 JM=6,10
            DO 160 KM=6,JM
              IF(KM.GE.9)GO TO 140
              KMM=KM-2
              IF(JM.GE.9) GO TO 130
              JMM=JM-2
              FTEMP=(((-PH(KMM)+BT**2*ALFA0*OMGAP*PHP(KMM))*
     &        (-PH(JMM)+BT**2*ALFA0*OMGAP*PHP(JMM))/BT)+ALFA0*
     &        BT*PHP(KMM)*PHP(JMM)*((BT*UP-OMGA)+(BT**2*ALFA0*OMGAP**2)
     &        *.5)+ALFA1*BT**3*PHPP(KMM)*PHPP(JMM))*ALFA0
              GO TO 150
  130         JMM=JM-6
              FTEMP=(PSP(JMM)*(-PH(KMM)+BT**2*ALFA0*OMGAP*PHP(KMM)))
              FTEMP=ALFA0*FTEMP
              GO TO 150
  140         KMM=KM-6
              JMM=JM-6
              FTEMP=BT*PSP(KMM)*PSP(JMM)*ALFA0
  150         JMM5=JM-5
              KMM5=KM-5
              TEMP(JMM5,KMM5)=FTEMP
  160         CONTINUE

          DO JM=1,3
            FTEMP=-DL*DPDU*PH(JM)/BT
            GTEMP=-DL*DPDV*PH(JM)/BT
            ITPJM=IT+JM
            TERM=COEF(K)*(FTEMP+STOR(JM))
            jac(ITPJM,NEQN)=jac(ITPJM,NEQN)+TERM
            TERM=COEF(K)*(GTEMP+GTOR(JM))
            jac(ITPJM,NVAR)=jac(ITPJM,NVAR)+TERM
          end do

          DO JM=6,8
            JMM2=JM-2
            JMM5=JM-5
            STOR(JMM5)=-DL*DPDU*PH(JMM2)/BT
            GTOR(JMM5)=-DL*DPDV*PH(JMM2)/BT
          end do

  190     CONTINUE
c
c  adjust the lowest right corner of the matrix.
c
        jac(NCOFM2,NCOFM2)=jac(NCOFM2,NCOFM2)+COEF(K)*TEMP(3,3)
        jac(NCOF,NCOF)=jac(NCOF,NCOF)+COEF(K)*TEMP(5,5)
        jac(NCOF,NCOFM2)=jac(NCOF,NCOFM2)+COEF(K)*TEMP(5,3)
        jac(NCOFM2,NEQN)=jac(NCOFM2,NEQN)+COEF(K)*STOR(3)
        jac(NCOFM2,NVAR)=jac(NCOFM2,NVAR)+COEF(K)*GTOR(3)
  200   CONTINUE
c
c  satisfy the initial conditions.
c
      jac(1,1)=1.0
      jac(2,2)=1.0
      jac(4,4)=1.0
      jac(ncofm4,ncofm4)=1.0
      jac(ncofm3,ncofm3)=1.0
      jac(ncofm1,ncofm1)=1.0
      jac(4,3)=0.0
      jac(ncofm2,ncofm4)=0.0
      jac(ncofm2,ncofm3)=0.0
      jac(ncof,ncofm4)=0.0
      jac(ncof,ncofm3)=0.0
      jac(ncof,ncofm1)=0.0

      do i=2,10
        jac(i,1)=0.0
      end do

      do i=3,10
        jac(i,2)=0.0
      end do

      do i=5,10
        jac(i,4)=0.0
      end do

      do i=ncofm9,ncofm5
        jac(ncofm4,i)=0.0
      end do

      do i=ncofm9,ncofm4
        jac(ncofm3,i)=0.0
      end do

      do i=ncofm9,ncofm2
        jac(ncofm1,i)=0.0
      end do

      do i=4,ncofm1,5
        jac(i,neqn)=0.0
        jac(i,nvar)=0.0
        ip1=i+1
        jac(ip1,neqn)=0.0
        jac(ip1,nvar)=0.0
      end do

      jac(1,neqn)=0.0
      jac(1,nvar)=0.0
      jac(2,neqn)=0.0
      jac(2,nvar)=0.0
      jac(ncofm3,neqn)=0.0
      jac(ncofm3,nvar)=0.0
      jac(ncofm4,neqn)=0.0
      jac(ncofm4,nvar)=0.0

      do i=1,ncofm1
        ip1=i+1
        do k=ip1,ncof
          jac(i,k)=jac(k,i)
        end do
      end do
c
c  last equation was fx(neqn)=x(ifix1)-val1=0.0
c
      jac(neqn,ifix1)=1.0
      return
      end
      subroutine p28_gp (xpos,x,nvar,dpdu,dpdv)

c*********************************************************************72
c
cc P28_GP is an auxilliary routine for Walker's arch.
c
      DIMENSION X(NVAR)
      COMMON /AUXMEM/ SINT(6),COEF(6),DL,ALFA0,THETA,ALFA1,IFUN17,
     & IFREE,IFIX1,VAL1,INTRVL,WIDE,NGAUSS,RADIUS
      NEQN=NVAR-1
      IF(IFUN17.EQ.1)GO TO 10
      IF(IFUN17.EQ.2)GO TO 20
      IF(IFUN17.EQ.3)GO TO 30
      IF(IFUN17.EQ.4)GO TO 40
      IF(IFUN17.EQ.5)GO TO 50
      IF(IFUN17.EQ.6)GO TO 60
   10 DPDU=1.0
      X1=-0.25*THETA
      X2=0.0
      IF(XPOS.LT.X1.OR.XPOS.GT.X2)DPDV=0.0
      IF(XPOS.GE.X1.AND.XPOS.LE.X2)DPDV=1.0
      return
   20 DPDU=0.0
      DPDV=1.0
      return
   30 IF(XPOS.LE.0.0)DPDU=1.0-X(NVAR)
      IF(XPOS.GT.0.0)DPDU=1.0+X(NVAR)
      IF(XPOS.LE.0.0)DPDV=-X(NEQN)
      IF(XPOS.GT.0.0)DPDV=X(NEQN)
      return
   40 DPDU=0.0
      X1=-0.75*THETA
      X2=-0.50*THETA
      IF(XPOS.GE.X1.AND.XPOS.LE.X2)DPDV=15.0/8.0
      IF(XPOS.LT.X1.OR.XPOS.GT.X2)DPDV=7.0/8.0
      return
   50 HEIGHT=X(NEQN)
      CENTER=X(NVAR)
      DPDU=0.0
      IF(XPOS.GT.CENTER.AND.XPOS.LE.(CENTER+WIDE))
     & DPDU=1.0+(CENTER-XPOS)/WIDE
      IF(XPOS.LE.CENTER.AND.XPOS.GE.(CENTER-WIDE))
     & DPDU=1.0-(CENTER-XPOS)/WIDE
      DPDU=2.0*THETA*DPDU/WIDE
      DPDV=0.0
      IF(XPOS.GT.CENTER.AND.XPOS.LE.(CENTER+WIDE))
     & DPDV=HEIGHT/WIDE
      IF(XPOS.LE.CENTER.AND.XPOS.GE.(CENTER-WIDE))
     & DPDV=-HEIGHT/WIDE
      DPDV=2.0*THETA*DPDV/WIDE
      return
   60 IF(XPOS.LT.X(NVAR))DPDU=0.0
      IF(XPOS.LT.X(NVAR))DPDV=0.0
      IF(XPOS.GT.X(NVAR))DPDU=1.0
      IF(XPOS.GT.X(NVAR))DPDV=0.0
      return
      END
      subroutine p28_vl (INT,X,NVAR,PH,PHP,PHPP,PS,PSP,OMGA,OMGAP,
     &  OMGAPP,UX,UP)

c*********************************************************************72
c
cc P28_VL is an auxilliary function for Walker's arch.
c
      DIMENSION X(NVAR),PH(6),PHP(6),PHPP(6),PS(4),PSP(4)
c
c  find omga,omgap,omgapp and up in this interval.
c
      omga=0.0
      omgap=0.0
      omgapp=0.0
      ux=0.0
      up=0.0

      do im=1,6
        km=im
        if (im.gt.3) km=im+2
        index=5*(int-1)+km
        z=x(index)
        omga=omga+ph(im)*z
        omgap=omgap+php(im)*z
        omgapp=omgapp+phpp(im)*z
      end do

      do im=1,4
        km=im+3
        if(im.gt.2) km=km+3
        index=5*(int-1)+km
        z=x(index)
        ux=ux+ps(im)*z
        up=up+psp(im)*z
      end do

      return
      end
      subroutine p28_nvar ( option, nvar )

c*********************************************************************72
c
cc P28_NVAR sets the number of variables for Walker's arch.
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
c    Input, integer OPTION, the option chosen for this problem.
c    For some problems, several options are available.  At least,
c    OPTION = 1 is always legal.
c
c    Output, integer NVAR, the number of variables.
c
      implicit none

      integer nvar
      integer option

      nvar = 47

      return
      end
      subroutine p28_option_num ( option_num )

c*********************************************************************72
c
cc P28_OPTION_NUM returns the number of options for Walker's arch.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    16 September 2008
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, integer OPTION_NUM, the number of options.
c
      implicit none

      integer option_num

      option_num = 10

      return
      end
      subroutine p28_start ( option, nvar, x )

c*********************************************************************72
c
cc P28_START returns a starting point for Walker's arch.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    19 September 2008
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer OPTION, the option index.
c
c    Input, integer NVAR, the number of variables.
c
c    Output, double precision X(NVAR), the starting point.
c
      implicit none

      integer nvar

      integer option
      double precision x(nvar)

      if ( ? ) then

      else

        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'P28_START - Fatal error!'
        write ( *, '(a,i8)' ) '  Unrecognized option  = ', option
        stop

      end if

      return
      end
      subroutine p28_stepsize ( option, h, hmin, hmax )

c*********************************************************************72
c
cc P28_STEPSIZE returns step sizes for Walker's arch.
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
c    Input, integer OPTION, the option index.
c
c    Output, double precision H, HMIN, HMAX, suggested values for the
c    initial step, the minimum step, and the maximum step.
c
      implicit none

      double precision h
      double precision hmax
      double precision hmin
      integer option

      h =    0.250D+00
      hmin = 0.001D+00
      hmax = 0.500D+00

      return
      end
      subroutine p28_title ( option, title )

c*********************************************************************72
c
cc P28_TITLE sets the title for Walker's arch.
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
c    Input, integer OPTION, the option index.
c
c    Output, character*(*) TITLE, the title of the problem.
c    TITLE will never be longer than 80 characters.
c
      implicit none

      integer option
      character*(*) title

      if ( option .eq. 1 ) then
        title = 'Walker''s arch, load 1, fix u = 0.1'
      else if ( option .eq. 2 ) then
        title = 'Walker''s arch, load 1, fix v = 0.01'
      else if ( option .eq. 3 ) then
        title = 'Walker''s arch, load 2, fix u = 0.0'
      else if ( option .eq. 4 ) then
        title = 'Walker''s arch, load 3, fix u = 1.0'
      else if ( option .eq. 5 ) then
        title = 'Walker''s arch, load 3, fix v = 0.1'
      else if ( option .eq. 6 ) then
        title = 'Walker''s arch, load 4, fix u = 0.0'
      else if ( option .eq. 7 ) then
        title = 'Walker''s arch, load 5, fix v, 0.0'
      else if ( option .eq. 8 ) then
        title = 'Walker''s arch, load 5, fix u.'
      else if ( option .eq. 9 ) then
        title = 'Walker''s arch, interactive.'
      else if ( option .eq. 10 ) then
        title = 'Walker''s arch, load 2, fix u = 0.0'
      else
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'P28_TITLE - Fatal error!'
        write ( *, '(a)' ) '  Unrecognized option number.'
        stop
      end if

      return
      end
      subroutine p29_data (DIRIPC,HFACT,HMAX,HMIN,HTAN,option,
     & IPC,IPRAM,ISIZE,IT,LIM,MAXCON,
     & MAXLIM,MAXSTP,MAXTAR,NCOL,NROW,NVAR,RELERR,XIT,RWORK)

c*********************************************************************72
c
cc P29_DATA sets parameters for the trussed dome function.
c
c  The trussed dome.
c
c  Reference:
c
c  y hangai, s kawamata,
c  analysis of geometrically nonlinear and stability problems
c    by static perturbation method,
c  department of the institute of industrial science,
c  university of tokyo, volume 22, number 5, january 1973.
c
c  The function:
c
c  the general problem involves a three dimensional structure
c  of truss elements.  small strains and large displacements
c  are assumed.  the cross-sectional areas of the trusses
c  are assumed constant.  the unknowns x describe the
c  displacements of the endpoints from a given starting position.
c  the load is multiplied by x(nvar), so that the null vector
c  is a satisfactory starting point.
c
c  the particular problem under consideration is a 3-dimensional box
c  with eight nodes, of which 7 are free.  there are 12 trusses in the
c  usual positions.  a downward load is applied to the top nodes.
c  the node at zero is fixed.
c
c  the function is of the form
c
c  fx(i)=sum k(x(i)-x(j))*(x(i)-x(j)) - x(nvar)*force(i)
c
c  where the sum is over all j such that there is a truss to node i.
c
c  Options:
c
c  option=1  it=nvar, xit=1.0
c
c  Comments:
c
c  code is not ready, jacobian unfinished, all needs revised.
c  possibly the code may run with ijac=1.
c
      DIMENSION RWORK(ISIZE)
      COMMON /AUXMEM/ EA,NODES,NTRUSS,XYZ(3,20),FORCE(3,20),IEND(2,20)

      NODES=7
      NVAR=3*NODES+1
      NEQN=NVAR-1
      NTRUSS=12
      EA=1.0
      DIRIPC=1.0
      HFACT=3.0
      IPC=NVAR
      IT=NVAR
      LIM=NVAR
      MAXCON=30
      MAXLIM=30
      MAXSTP=30
      MAXTAR=30
      NCOL=NVAR
      NROW=NVAR
      RELERR=1.0E-5
      XIT=1.0
c
c  set node data
c
c  for node 1 (1,1,1)
c
      XYZ(1,1)=1.0
      XYZ(2,1)=1.0
      XYZ(3,1)=1.0
      FORCE(1,1)=-1.0
      FORCE(2,1)=0.0
      FORCE(3,1)=0.0
c
c  for node 2 (1,1,0)
c
      XYZ(1,2)=1.0
      XYZ(2,2)=1.0
      XYZ(3,2)=0.0
      FORCE(1,2)=-1.0
      FORCE(2,2)=0.0
      FORCE(3,2)=0.0
c
c  for node 3 (1,0,1)
c
      XYZ(1,3)=1.0
      XYZ(2,3)=0.0
      XYZ(3,3)=1.0
      FORCE(1,3)=-1.0
      FORCE(2,3)=0.0
      FORCE(3,3)=0.0
c
c  for node 4 (1,0,0)
c
      XYZ(1,4)=1.0
      XYZ(2,4)=0.0
      XYZ(3,4)=0.0
      FORCE(1,4)=-1.0
      FORCE(2,4)=0.0
      FORCE(3,4)=0.0
c
c  for node 5 (0,1,1)
c
      XYZ(1,5)=0.0
      XYZ(2,5)=1.0
      XYZ(3,5)=1.0
      FORCE(1,5)=0.0
      FORCE(2,5)=0.0
      FORCE(3,5)=0.0
c
c  for node 6 (0,1,0)
c
      XYZ(1,6)=0.0
      XYZ(2,6)=1.0
      XYZ(3,6)=0.0
      FORCE(1,6)=0.0
      FORCE(2,6)=0.0
      FORCE(3,6)=0.0
c
c  for node 7 (0,0,1)
c
      XYZ(1,7)=0.0
      XYZ(2,7)=0.0
      XYZ(3,7)=1.0
      FORCE(1,7)=0.0
      FORCE(2,7)=0.0
      FORCE(3,7)=0.0
c
c  for node 8 (0,0,0)
c
      XYZ(1,8)=0.0
      XYZ(2,8)=0.0
      XYZ(3,8)=0.0
c
c  set endpoints of trusses
c
      IEND(1,1)=1
      IEND(2,1)=2
      IEND(1,2)=1
      IEND(2,2)=3
      IEND(1,3)=1
      IEND(2,3)=5
      IEND(1,4)=4
      IEND(2,4)=2
      IEND(1,5)=4
      IEND(2,5)=3
      IEND(1,6)=4
      IEND(2,6)=8
      IEND(1,7)=6
      IEND(2,7)=2
      IEND(1,8)=6
      IEND(2,8)=5
      IEND(1,9)=6
      IEND(2,9)=8
      IEND(1,10)=7
      IEND(2,10)=3
      IEND(1,11)=7
      IEND(2,11)=5
      IEND(1,12)=7
      IEND(2,12)=8
      return
      END
      subroutine p29_fun ( option, nvar, x, fx )

c*********************************************************************72
c
cc P29_FUN evaluates the trussed dome function.
c
      integer nvar

      double precision fx(nvar-1)
      integer option
      double precision x(nvar)

      DIMENSION BARRAY(3,3),Z(3)
      COMMON /AUXMEM/ EA,NODES,NTRUSS,XYZ(3,20),FORCE(3,20),IEND(2,20)
      XNVAR=X(NVAR)
      NEQN=NVAR-1

      DO I=1,NEQN
        FX(I)=0.0
      end do
c
c  for each node i, find node j with truss between.
c
      DO I=1,NODES
        DO K=1,NTRUSS
          J=0
          IF(IEND(1,K).EQ.I)J=IEND(2,K)
          IF(IEND(2,K).EQ.I)J=IEND(1,K)
          IF(J.EQ.0)GO TO 40
c
c  node j is a truss neighbor
c
          INDEXI=3*(I-1)+1
          INDEXJ=3*(J-1)+1
          call p29_gx (INDEXI,X(INDEXI),INDEXJ,X(INDEXJ),Z,BARRAY)
          DO L=1,3
            DO M=1,3
              FX(INDEXI+L-1)=FX(INDEXI+L-1)+BARRAY(L,M)*Z(M)
              FX(INDEXJ+L-1)=FX(INDEXJ+L-1)-BARRAY(L,M)*Z(M)
            end do
          end do

        end do
      end do

      DO I=1,NODES
        IEQN=3*(I-1)+1
        FX(IEQN)=FX(IEQN)+XNVAR*FORCE(1,I)
        IEQN=IEQN+1
        FX(IEQN)=FX(IEQN)+XNVAR*FORCE(2,I)
        IEQN=IEQN+1
        FX(IEQN)=FX(IEQN)+XNVAR*FORCE(3,I)
      end do

      return
      END
      subroutine p29_gx (INDEXI,XI,INDEXJ,XJ,Z,BARRAY)

c*********************************************************************72
c
cc P29_GX is an auxilliary routine for the trussed dome function.
c
      DIMENSION Z(3),ZUNIT(3),DZUNIT(3)
      DIMENSION XI(3),XJ(3)
      DIMENSION BARRAY(3,3),CARRAY(3,3),DARRAY(3,3)
      COMMON /AUXMEM/ EA,NODES,NTRUSS,XYZ(3,20),FORCE(3,20),IEND(2,20)
c
c  some quantities below depend on xyz, others on x.
c
c  set z
c
      DO I=1,3
        Z(I)=XI(I)-XJ(I)
      end do
c
c  set znorm, zunit, delta
c
      ZNORM=SQRT(Z(1)*Z(1)+Z(2)*Z(2)+Z(3)*Z(3))
      IF(ZNORM.EQ.0.0)STOP
      DO I=1,3
        ZUNIT(I)=Z(I)/ZNORM
      end do
      DELTA=SQRT(ZUNIT(1)*ZUNIT(1)+ZUNIT(2)*ZUNIT(2))
c
c  set darray
c
      IF(DELTA.NE.0.0)GO TO 50

      DO I=1,3
        DO J=1,3
          DARRAY(I,J)=0.0
        end do
      end do

      DARRAY(3,1)=-1.0
      DARRAY(2,2)=1.0
      DARRAY(1,3)=1.0
      GO TO 60
   50 DARRAY(1,1)=ZUNIT(1)
      DARRAY(1,2)=ZUNIT(2)
      DARRAY(1,3)=ZUNIT(3)
      DARRAY(2,1)=-ZUNIT(2)/DELTA
      DARRAY(2,2)=ZUNIT(1)/DELTA
      DARRAY(2,3)=0.0
      DARRAY(3,1)=-ZUNIT(1)*ZUNIT(3)/DELTA
      DARRAY(3,2)=-ZUNIT(2)*ZUNIT(3)/DELTA
      DARRAY(3,3)=DELTA
   60 CONTINUE
c
c  accumulate entries of carray=2*m0+m1+m2
c
      DO I=1,3
        DO J=1,3
          CARRAY(I,J)=0.0
        end do
      end do
c
c  contribution from m0
c
      CARRAY(1,1)=2.0
c
c  contribution from m1
c  (first compute dzunit)
c
      DO I=1,3
        DZUNIT(I)=0.0
        DO J=1,3
          DZUNIT(I)=DZUNIT(I)+DARRAY(I,J)*ZUNIT(J)
        end do
      end do

      DO I=1,3
        CARRAY(I,I)=CARRAY(I,I)+DZUNIT(1)
        CARRAY(1,I)=CARRAY(1,I)+DZUNIT(I)
        CARRAY(I,1)=CARRAY(I,1)+DZUNIT(I)
      end do
c
c  contribution from m2
c
      DO I=1,3
        DO J=1,3
          CARRAY(I,J)=CARRAY(I,J)+DZUNIT(I)*DZUNIT(J)
        end do
      end do
c
c  set barray
c
      DO I=1,3
        DO J=1,3
          SUM=0.0
          DO K=1,3
            DO L=1,3
              SUM=SUM+DARRAY(K,I)*CARRAY(K,L)*DARRAY(L,J)
            end do
          end do
          BARRAY(I,J)=.5*EA*SUM
        end do
      end do

      return
      END
      subroutine p29_jac ( option, nvar, x, jac )

c*********************************************************************72
c
cc P29_JAC evaluates the trussed dome jacobian.
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
c    Input, integer OPTION, the option index.
c
c    Input, integer NVAR, the number of variables.
c
c    Input, double precision X(NVAR), the argument of the jacobian.
c
c    Output, double precision JAC(NVAR,NVAR), the jacobian matrix evaluated
c    at X.  The NVAR-th row is not set by this routine.
c
      implicit none

      integer nvar

      integer option
      double precision jac(nvar,nvar)
      double precision x(nvar)

      double precision BARRAY(3,3)
      COMMON /AUXMEM/ EA,NODES,NTRUSS,XYZ(3,20),FORCE(3,20),IEND(2,20)

      do j = 1, ncol
        do i = 1, nrow
          jac(i,j) = 0.0E+00
        end do
      end do

      xnvar=x(nvar)
      neqn=nvar-1
c
c  the derivative has the form k(x)+kp(x)*x
c  can compute the k(x) term now.
c
c  for each node i, find node j with truss between.
c
      do i=1,nodes
        do k=1,ntruss
          j=0
          if(iend(1,k).eq.i)j=iend(2,k)
          if(iend(2,k).eq.i)j=iend(1,k)
          if(j.eq.0)go to 50
c
c  node j is a truss neighbor
c
          indexi=3*(i-1)+1
          indexj=3*(j-1)+1
          call p29_gx (indexi,x(indexi),indexj,x(indexj),z,barray)
          do l=1,3
            do m=1,3
c
c  the indexing is off here, and there should be four
c  contributions, not two
c
              jac(indexi+l-1,m)=jac(indexi+l-1,m)+barray(l,m)
              jac(indexj+l-1,m)=jac(indexj+l-1,m)-barray(l,m)
            end do
          end do
        end do
      end do
c
c  d/d(x(nvar) of x(nvar)*force(i)
c
      do i=1,nodes
        ieqn=3*(i-1)+1
        jac(ieqn,nvar)=jac(ieqn,nvar)+force(1,i)
        ieqn=ieqn+1
        jac(ieqn,nvar)=jac(ieqn,nvar)+force(2,i)
        ieqn=ieqn+1
        jac(ieqn,nvar)=jac(ieqn,nvar)+force(3,i)
      end do

      return
      end
      subroutine p29_nvar ( option, nvar )

c*********************************************************************72
c
cc P29_NVAR sets the number of variables for the trussed dome.
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
c    Input, integer OPTION, the option chosen for this problem.
c    For some problems, several options are available.  At least,
c    OPTION = 1 is always legal.
c
c    Output, integer NVAR, the number of variables.
c
      implicit none

      integer nvar
      integer option

      nvar = 22

      return
      end
      subroutine p29_option_num ( option_num )

c*********************************************************************72
c
cc P29_OPTION_NUM returns the number of options for the trussed dome.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    16 September 2008
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, integer OPTION_NUM, the number of options.
c
      implicit none

      integer option_num

      option_num = 1

      return
      end
      subroutine p29_start ( option, nvar, x )

c*********************************************************************72
c
cc P29_START returns a starting point for the trussed dome.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    19 September 2008
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer OPTION, the option index.
c
c    Input, integer NVAR, the number of variables.
c
c    Output, double precision X(NVAR), the starting point.
c
      implicit none

      integer nvar

      integer i
      integer option
      double precision x(nvar)

      if ( option .eq. 1 ) then

        do i = 1, nvar
          x(i) = 0.0
        end do

      else

        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'P29_START - Fatal error!'
        write ( *, '(a,i8)' ) '  Unrecognized option  = ', option
        stop

      end if

      return
      end
      subroutine p29_stepsize ( option, h, hmin, hmax )

c*********************************************************************72
c
cc P29_STEPSIZE returns step sizes for the trussed dome.
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
c    Input, integer OPTION, the option index.
c
c    Output, double precision H, HMIN, HMAX, suggested values for the
c    initial step, the minimum step, and the maximum step.
c
      implicit none

      double precision h
      double precision hmax
      double precision hmin
      integer option

      h =     2.000D+00
      hmin =  0.001D+00
      hmax = 10.000D+00

      return
      end
      subroutine p29_title ( option, title )

c*********************************************************************72
c
cc P29_TITLE sets the title for the trussed dome.
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
c    Input, integer OPTION, the option index.
c
c    Output, character*(*) TITLE, the title of the problem.
c    TITLE will never be longer than 80 characters.
c
      implicit none

      integer option
      character*(*) title

      if ( option .eq. 1 ) then
        title = 'The trussed dome function.'
      else
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'P29_TITLE - Fatal error!'
        write ( *, '(a)' ) '  Unrecognized option number.'
        stop
      end if

      return
      end
      subroutine timestamp ( )

c*********************************************************************72
c
cc TIMESTAMP prints out the current YMDHMS date as a timestamp.
c
c  Discussion:
c
c    This FORTRAN77 version is made available for cases where the
c    FORTRAN90 version cannot be used.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    12 January 2007
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    None
c
      implicit none

      character * ( 8 ) ampm
      integer d
      character * ( 8 ) date
      integer h
      integer m
      integer mm
      character * ( 9 ) month(12)
      integer n
      integer s
      character * ( 10 ) time
      integer y

      save month

      data month /
     &  'January  ', 'February ', 'March    ', 'April    ',
     &  'May      ', 'June     ', 'July     ', 'August   ',
     &  'September', 'October  ', 'November ', 'December ' /

      call date_and_time ( date, time )

      read ( date, '(i4,i2,i2)' ) y, m, d
      read ( time, '(i2,i2,i2,1x,i3)' ) h, n, s, mm

      if ( h .lt. 12 ) then
        ampm = 'AM'
      else if ( h .eq. 12 ) then
        if ( n .eq. 0 .and. s .eq. 0 ) then
          ampm = 'Noon'
        else
          ampm = 'PM'
        end if
      else
        h = h - 12
        if ( h .lt. 12 ) then
          ampm = 'PM'
        else if ( h .eq. 12 ) then
          if ( n .eq. 0 .and. s .eq. 0 ) then
            ampm = 'Midnight'
          else
            ampm = 'AM'
          end if
        end if
      end if

      write ( *,
     &  '(i2,1x,a,1x,i4,2x,i2,a1,i2.2,a1,i2.2,a1,i3.3,1x,a)' )
     &  d, month(m), y, h, ':', n, ':', s, '.', mm, ampm

      return
      end
