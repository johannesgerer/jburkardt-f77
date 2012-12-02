      subroutine p00_f ( prob, n, x, value )

c*********************************************************************72
c
cc P00_F evaluates the function for any problem.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    28 August 2012
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer PROB, the number of the desired test problem.
c
c    Input, integer N, the number of evaluation points.
c
c    Input, double precision X(N), the evaluation points.
c
c    Output, double precision VALUE(N), the function values.
c
      implicit none

      integer n

      integer prob
      double precision value(n)
      double precision x(n)

      if ( prob .eq. 1 ) then
        call p01_f ( n, x, value )
      else if ( prob .eq. 2 ) then
        call p02_f ( n, x, value )
      else if ( prob .eq. 3 ) then
        call p03_f ( n, x, value )
      else if ( prob .eq. 4 ) then
        call p04_f ( n, x, value )
      else if ( prob .eq. 5 ) then
        call p05_f ( n, x, value )
      else if ( prob .eq. 6 ) then
        call p06_f ( n, x, value )
      else if ( prob .eq. 7 ) then
        call p07_f ( n, x, value )
      else if ( prob .eq. 8 ) then
        call p08_f ( n, x, value )
      else
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'P00_F - Fatal error!'
        write ( *, '(a,i4)' ) '  Illegal problem number = ', prob
        stop
      end if

      return
      end
      subroutine p00_prob_num ( prob_num )

c*********************************************************************72
c
cc P00_PROB_NUM returns the number of problems.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    28 August 2012
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, integer PROB_NUM, the number of problems.
c
      implicit none

      integer prob_num

      prob_num = 8

      return
      end
      subroutine p00_title ( prob, title )

c*********************************************************************72
c
cc P00_TITLE returns the title of any problem.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    28 August 2012
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer PROB, the number of the desired test problem.
c
c    Output, character * ( * ) TITLE, the title of the problem.
c
      implicit none

      integer prob
      character * ( * ) title

      if ( prob .eq. 1 ) then
        call p01_title ( title )
      else if ( prob .eq. 2 ) then
        call p02_title ( title )
      else if ( prob .eq. 3 ) then
        call p03_title ( title )
      else if ( prob .eq. 4 ) then
        call p04_title ( title )
      else if ( prob .eq. 5 ) then
        call p05_title ( title )
      else if ( prob .eq. 6 ) then
        call p06_title ( title )
      else if ( prob .eq. 7 ) then
        call p07_title ( title )
      else if ( prob .eq. 8 ) then
        call p08_title ( title )
      else
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'P00_TITLE - Fatal error!'
        write ( *, '(a,i4)' ) '  Illegal problem number = ', prob
        stop
      end if

      return
      end
      subroutine p01_f ( n, x, value )

c*********************************************************************72
c
cc P01_F evaluates the function for problem p01.
c
c  Discussion:
c
c    Step function.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    28 August 2012
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the number of evaluation points.
c
c    Input, double precision X(N), the evaluation points.
c
c    Output, double precision VALUE(N), the function values.
c
      implicit none

      integer n

      integer i
      double precision value(n)
      double precision x(n)

      do i = 1, n
        if ( x(i) .le. 1.0D+00 / 3.0D+00 ) then
          value(i) = -1.0D+00
        else if ( x(i) .le. 4.0D+00 / 5.0D+00 ) then
          value(i) = 2.0D+00
        else
          value(i) = 1.0D+00
        end if
      end do

      return
      end
      subroutine p01_title ( title )

c*********************************************************************72
c
cc P01_TITLE returns the title of problem p01.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    28 August 2012
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, character * ( * ) TITLE, the title of the problem.
c
      implicit none

      character * ( * ) title

      title = 'f(x) = steps -1/2/1 at [0,1/3], [1/3,4/5], [4/5,1].'

      return
      end
      subroutine p02_f ( n, x, value )

c*********************************************************************72
c
cc P02_F evaluates the function for problem p01.
c
c  Discussion:
c
c    Nondifferentiable function.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    28 August 2012
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the number of evaluation points.
c
c    Input, double precision X(N), the evaluation points.
c
c    Output, double precision VALUE(N), the function values.
c
      implicit none

      integer n

      integer i
      double precision value(n)
      double precision x(n)

      do i = 1, n
        if ( x(i) .le. 1.0D+00 / 3.0D+00 ) then
          value(i) = 1.0D+00 - 3.0D+00 * x(i)
        else
          value(i) = 6.0D+00 * x(i) - 2.0D+00
        end if
      end do

      return
      end
      subroutine p02_title ( title )

c*********************************************************************72
c
cc P02_TITLE returns the title of problem p02.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    28 August 2012
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, character * ( * ) TITLE, the title of the problem.
c
      implicit none

      character * ( * ) title

      title = 'f(x) = (1-3x), x < 1/3 (6x-2) if 1/3 < x'

      return
      end
      subroutine p03_f ( n, x, value )

c*********************************************************************72
c
cc P03_F evaluates the function for problem p03.
c
c  Discussion:
c
c    Step function.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    28 August 2012
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the number of evaluation points.
c
c    Input, double precision X(N), the evaluation points.
c
c    Output, double precision VALUE(N), the function values.
c
      implicit none

      integer n

      integer j
      double precision value(n)
      double precision x(n)

      do j = 1, n
        value(j) = x(j) 
     &    * ( 10.0D+00 * x(j) - 1.0D+00 ) 
     &    * ( 5.0D+00 * x(j) - 2.0D+00 ) 
     &    * ( 5.0D+00 * x(j) - 2.0D+00 ) 
     &    * ( 4.0D+00 * x(j) - 3.4D+00 ) 
     &    * ( x(j) - 1.0D+00 )
      end do

      return
      end
      subroutine p03_title ( title )

c*********************************************************************72
c
cc P03_TITLE returns the title of problem p03.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    28 August 2012
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, character * ( * ) TITLE, the title of the problem.
c
      implicit none

      character * ( * ) title

      title = 'f(x) = x (10*x-1) (5x-2) (5x-2) (4x-3.4) (x-1)'

      return
      end
      subroutine p04_f ( n, x, value )

c*********************************************************************72
c
cc P04_F evaluates the function for problem p04.
c
c  Discussion:
c
c    Step function.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    28 August 2012
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the number of evaluation points.
c
c    Input, double precision X(N), the evaluation points.
c
c    Output, double precision VALUE(N), the function values.
c
      implicit none

      integer n

      integer j 
      double precision value(n)
      double precision x(n)

      do j = 1, n
        value(j) = atan ( 40.0D+00 * x(j) - 15.0D+00 )
      end do

      return
      end
      subroutine p04_title ( title )

c*********************************************************************72
c
cc P04_TITLE returns the title of problem p04.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    28 August 2012
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, character * ( * ) TITLE, the title of the problem.
c
      implicit none

      character * ( * ) title

      title = 'f(x) = atan ( 40 * x - 15 )'

      return
      end
      subroutine p05_f ( n, x, value )

c*********************************************************************72
c
cc P05_F evaluates the function for problem p05.
c
c  Discussion:
c
c    Step function.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    28 August 2012
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the number of evaluation points.
c
c    Input, double precision X(N), the evaluation points.
c
c    Output, double precision VALUE(N), the function values.
c
      implicit none

      integer n

      integer j
      double precision value(n)
      double precision x(n)

      do j = 1, n
         value(j) =          cos (  7.0D+00 * x(j) ) 
     &           + 5.0D+00 * cos ( 11.2D+00 * x(j) ) 
     &           - 2.0D+00 * cos ( 14.0D+00 * x(j) ) 
     &           + 5.0D+00 * cos ( 31.5D+00 * x(j) ) 
     &           + 7.0D+00 * cos ( 63.0D+00 * x(j) )
      end do

      return
      end
      subroutine p05_title ( title )

c*********************************************************************72
c
cc P05_TITLE returns the title of problem p05.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    28 August 2012
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, character * ( * ) TITLE, the title of the problem.
c
      implicit none

      character * ( * ) title

      title = 'f(x) = cos(7*x)+5*cos(11.2*x)-2*cos(14*x)' //
     &  '+5*cos(31.5*x)+7*cos(63*x).'

      return
      end
      subroutine p06_f ( n, x, value )

c*********************************************************************72
c
cc P06_F evaluates the function for problem p06.
c
c  Discussion:
c
c    f(x) = exp ( - (4 * x - 1)^2 )
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    28 August 2012
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Alan Genz,
c    A Package for Testing Multiple Integration Subroutines,
c    in Numerical Integration: Recent Developments, Software
c    and Applications,
c    edited by Patrick Keast and Graeme Fairweather,
c    D Reidel, 1987, pages 337-340,
c    LC: QA299.3.N38.
c
c  Parameters:
c
c    Input, integer N, the number of points.
c
c    Input, double precision X(N), the evaluation points.
c
c    Output, double precision F(N), the function values.
c
      implicit none

      integer n

      integer j 
      double precision value(n)
      double precision x(n)

      do j = 1, n
        value(j) = exp ( - ( 4.0D+00 * x(j) - 1.0D+00 ) ** 2 )
      end do

      return
      end
      subroutine p06_title ( title )

c*********************************************************************72
c
cc P06_TITLE returns the title of problem p06.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    28 August 2012
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, character * ( * ) TITLE, the title of the problem.
c
      implicit none

      character * ( * ) title

      title = 'f(x) = exp ( - ( 4*x-1 )^2 )'

      return
      end
      subroutine p07_f ( n, x, value )

c*********************************************************************72
c
cc P07_F evaluates the function for problem p07.
c
c  Discussion:
c
c    f(x) = exp ( 4 * x ) if x .le. 1/2
c           0                  otherwise
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    28 August 2012
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Alan Genz,
c    A Package for Testing Multiple Integration Subroutines,
c    in Numerical Integration: Recent Developments, Software
c    and Applications,
c    edited by Patrick Keast and Graeme Fairweather,
c    D Reidel, 1987, pages 337-340,
c    LC: QA299.3.N38.
c
c  Parameters:
c
c    Input, integer N, the number of points.
c
c    Input, double precision X(N), the evaluation points.
c
c    Output, double precision F(N), the function values.
c
      implicit none

      integer n

      integer i
      double precision value(n)
      double precision x(n)

      do i = 1, n
        if ( x(i) .lt. 0.5D+00 ) then
          value(i) = exp ( 4.0D+00 * x(i) )
        else
          value(i) = 0.0D+00
        end if
      end do

      return
      end
      subroutine p07_title ( title )

c*********************************************************************72
c
cc P07_TITLE returns the title of problem p07.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    28 August 2012
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, character * ( * ) TITLE, the title of the problem.
c
      implicit none

      character * ( * ) title

      title = 'f(x) = exp ( 2 x ) if x < 0.5, 0 otherwise'

      return
      end
      subroutine p08_f ( n, x, value )

c*********************************************************************72
c
cc P08_F evaluates the function for problem p08.
c
c  Discussion:
c
c    This is a famous example, due to Runge.  If equally spaced
c    abscissas are used, the sequence of interpolating polynomials Pn(X)
c    diverges, in the sense that the max norm of the difference
c    between Pn(X) and F(X) becomes arbitrarily large as N increases.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    28 August 2012
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the number of evaluation points.
c
c    Input, double precision X(N), the evaluation points.
c
c    Output, double precision VALUE(N), the function values.
c
      implicit none

      integer n

      integer j
      double precision value(n)
      double precision x(n)

      do j = 1, n
        value(j) = 1.0D+00 
     &    / ( ( 10.0D+00 * ( x(j) - 0.5D+00 ) ) ** 2 + 1.0D+00 )
      end do

      return
      end
      subroutine p08_title ( title )

c*********************************************************************72
c
cc P08_TITLE returns the title of problem p08.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    28 August 2012
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, character * ( * ) TITLE, the title of the problem.
c
      implicit none

      character * ( * ) title

      title = 'f(x) = 1 / ( 1 + ( 10 * (x-1/2) )^2 )'

      return
      end
