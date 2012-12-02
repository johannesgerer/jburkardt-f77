      subroutine f00_f0 ( fi, n, x, y, f )

c*********************************************************************72
c
cc F00_F0 returns the value of any function.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    03 August 2012
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer FI, the index of the function.
c
c    Input, integer N, the number of evaluation points.
c
c    Input, double precision X(N), Y(N), the evalution points.
c
c    Output, double precision F(N), the function values.
c
      implicit none

      integer n

      double precision f(n)
      integer fi
      double precision x(n)
      double precision y(n)  

      if ( fi .eq. 1 ) then
        call f01_f0 ( n, x, y, f )
      else if ( fi .eq. 2 ) then
        call f02_f0 ( n, x, y, f )
      else if ( fi .eq. 3 ) then
        call f03_f0 ( n, x, y, f )
      else if ( fi .eq. 4 ) then
        call f04_f0 ( n, x, y, f )
      else if ( fi .eq. 5 ) then
        call f05_f0 ( n, x, y, f )
      else if ( fi .eq. 6 ) then
        call f06_f0 ( n, x, y, f )
      else if ( fi .eq. 7 ) then
        call f07_f0 ( n, x, y, f )
      else if ( fi .eq. 8 ) then
        call f08_f0 ( n, x, y, f )
      else if ( fi .eq. 9 ) then
        call f09_f0 ( n, x, y, f )
      else if ( fi .eq. 10 ) then
        call f10_f0 ( n, x, y, f )
      else if ( fi .eq. 11 ) then
        call f11_f0 ( n, x, y, f )
      else if ( fi .eq. 12 ) then
        call f12_f0 ( n, x, y, f )
      else if ( fi .eq. 13 ) then
        call f13_f0 ( n, x, y, f )
      else
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'F00_F0 - Fatal error!'
        write ( *, '(a,i8)' ) '  Illegal function index FI = ', fi
        stop
      end if

      return
      end
      subroutine f00_f1 ( fi, n, x, y, fx, fy )

c*********************************************************************72
c
cc F00_F1 returns first derivatives of any function.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    03 August 2012
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer FI, the index of the function.
c
c    Input, integer N, the number of evaluation points.
c
c    Input, double precision X(N), Y(N), the evalution points.
c
c    Output, double precision FX(N), FY(N), the first derivative values.
c
      implicit none

      integer n

      integer fi
      double precision fx(n)
      double precision fy(n)
      double precision x(n)
      double precision y(n)  

      if ( fi .eq. 1 ) then
        call f01_f1 ( n, x, y, fx, fy )
      else if ( fi .eq. 2 ) then
        call f02_f1 ( n, x, y, fx, fy )
      else if ( fi .eq. 3 ) then
        call f03_f1 ( n, x, y, fx, fy )
      else if ( fi .eq. 4 ) then
        call f04_f1 ( n, x, y, fx, fy )
      else if ( fi .eq. 5 ) then
        call f05_f1 ( n, x, y, fx, fy )
      else if ( fi .eq. 6 ) then
        call f06_f1 ( n, x, y, fx, fy )
      else if ( fi .eq. 7 ) then
        call f07_f1 ( n, x, y, fx, fy )
      else if ( fi .eq. 8 ) then
        call f08_f1 ( n, x, y, fx, fy )
      else if ( fi .eq. 9 ) then
        call f09_f1 ( n, x, y, fx, fy )
      else if ( fi .eq. 10 ) then
        call f10_f1 ( n, x, y, fx, fy )
      else if ( fi .eq. 11 ) then
        call f11_f1 ( n, x, y, fx, fy )
      else if ( fi .eq. 12 ) then
        call f12_f1 ( n, x, y, fx, fy )
      else if ( fi .eq. 13 ) then
        call f13_f1 ( n, x, y, fx, fy )
      else
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'F00_F1 - Fatal error!'
        write ( *, '(a,i8)' ) '  Illegal function index FI = ', fi
        stop
      end if

      return
      end
      subroutine f00_f2 ( fi, n, x, y, fxx, fxy, fyy )

c*********************************************************************72
c
cc F00_F2 returns second derivatives of any function.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    03 August 2012
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer FI, the index of the function.
c
c    Input, integer N, the number of evaluation points.
c
c    Input, double precision X(N), Y(N), the evalution points.
c
c    Output, double precision FXX(N), FXY(N), FYY(N), the second derivatives.
c
      implicit none

      integer n

      integer fi
      double precision fxx(n)
      double precision fxy(n)
      double precision fyy(n)
      double precision x(n)
      double precision y(n)  

      if ( fi .eq. 1 ) then
        call f01_f2 ( n, x, y, fxx, fxy, fyy )
      else if ( fi .eq. 2 ) then
        call f02_f2 ( n, x, y, fxx, fxy, fyy )
      else if ( fi .eq. 3 ) then
        call f03_f2 ( n, x, y, fxx, fxy, fyy )
      else if ( fi .eq. 4 ) then
        call f04_f2 ( n, x, y, fxx, fxy, fyy )
      else if ( fi .eq. 5 ) then
        call f05_f2 ( n, x, y, fxx, fxy, fyy )
      else if ( fi .eq. 6 ) then
        call f06_f2 ( n, x, y, fxx, fxy, fyy )
      else if ( fi .eq. 7 ) then
        call f07_f2 ( n, x, y, fxx, fxy, fyy )
      else if ( fi .eq. 8 ) then
        call f08_f2 ( n, x, y, fxx, fxy, fyy )
      else if ( fi .eq. 9 ) then
        call f09_f2 ( n, x, y, fxx, fxy, fyy )
      else if ( fi .eq. 10 ) then
        call f10_f2 ( n, x, y, fxx, fxy, fyy )
      else if ( fi .eq. 11 ) then
        call f11_f2 ( n, x, y, fxx, fxy, fyy )
      else if ( fi .eq. 12 ) then
        call f12_f2 ( n, x, y, fxx, fxy, fyy )
      else if ( fi .eq. 13 ) then
        call f13_f2 ( n, x, y, fxx, fxy, fyy )
      else
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'F00_F2 - Fatal error!'
        write ( *, '(a,i8)' ) '  Illegal function index FI = ', fi
        stop
      end if

      return
      end
      subroutine f00_num ( fun_num )

c*********************************************************************72
c
cc F00_NUM returns the number of test functions available.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    03 August 2012
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c   Output, integer FUN_NUM, the number of test functions.
c
      implicit none

      integer fun_num

      fun_num = 13

      return
      end
      subroutine f00_title ( fi, ft )

c*********************************************************************72
c
cc F00_TITLE returns the title for any function.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    03 August 2012
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer FI, the index of the function.
c
c    Output, character * ( * ) FT, the function title.
c
      implicit none

      integer fi
      character * ( * ) ft

      if ( fi .eq. 1 ) then
        call f01_title ( ft )
      else if ( fi .eq. 2 ) then
        call f02_title ( ft )
      else if ( fi .eq. 3 ) then
        call f03_title ( ft )
      else if ( fi .eq. 4 ) then
        call f04_title ( ft )
      else if ( fi .eq. 5 ) then
        call f05_title ( ft )
      else if ( fi .eq. 6 ) then
        call f06_title ( ft )
      else if ( fi .eq. 7 ) then
        call f07_title ( ft )
      else if ( fi .eq. 8 ) then
        call f08_title ( ft )
      else if ( fi .eq. 9 ) then
        call f09_title ( ft )
      else if ( fi .eq. 10 ) then
        call f10_title ( ft )
      else if ( fi .eq. 11 ) then
        call f11_title ( ft )
      else if ( fi .eq. 12 ) then
        call f12_title ( ft )
      else if ( fi .eq. 13 ) then
        call f13_title ( ft )
      else
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'F00_TITLE - Fatal error!'
        write ( *, '(a,i8)' ) '  Illegal function index FI = ', fi
        stop
      end if

      return
      end
      subroutine f01_f0 ( n, x, y, f )

c*********************************************************************72
c
cc F01_F0 returns the value of function f01.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    28 January 2012
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the number of evaluation points.
c
c    Input, double precision X(N), Y(N), the evalution points.
c
c    Output, double precision F(N), the function values.
c
      implicit none

      integer n

      double precision f(n)
      integer i
      double precision x(n)
      double precision y(n)  

      do i = 1, n
        f(i) = 
     &      0.75 * exp ( - ( ( 9.0 * x(i) - 2.0 )**2            
     &                     + ( 9.0 * y(i) - 2.0 )**2 ) / 4.0 )  
     &    + 0.75 * exp ( - ( ( 9.0 * x(i) + 1.0 )**2 ) / 49.0   
     &                     - ( 9.0 * y(i) + 1.0 ) / 10.0 )      
     &    + 0.5  * exp ( - ( ( 9.0 * x(i) - 7.0 )**2            
     &                     + ( 9.0 * y(i) - 3.0 )**2 ) / 4.0 )  
     &    - 0.2  * exp ( - (   9.0 * x(i) - 4.0 )**2            
     &                     - ( 9.0 * y(i) - 7.0 )**2 )
      end do

      return
      end
      subroutine f01_f1 ( n, x, y, fx, fy )

c*********************************************************************72
c
cc F01_F1 returns first derivatives of function f01.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    28 January 2012
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the number of evaluation points.
c
c    Input, double precision X(N), Y(N), the evalution points.
c
c    Output, double precision FX(N), FY(N), the derivative values.
c
      implicit none

      integer n

      double precision fx(n)
      double precision fy(n)
      integer i
      double precision t1
      double precision t2
      double precision t3
      double precision t4
      double precision x(n)
      double precision y(n)  

      do i = 1, n

        t1 = exp ( - ( ( 9.0 * x(i) - 2.0 )**2
     &                  + ( 9.0 * y(i) - 2.0 )**2 ) / 4.0 )
        t2 = exp ( - ( ( 9.0 * x(i) + 1.0 )**2 ) / 49.0
     &                  - ( 9.0 * y(i) + 1.0 ) / 10.0 )
        t3 = exp ( - ( ( 9.0 * x(i) - 7.0 )**2
     &                  + ( 9.0 * y(i) - 3.0 )**2 ) / 4.0 )
        t4 = exp ( -   ( 9.0 * x(i) - 4.0 )**2
     &                  - ( 9.0 * y(i) - 7.0 )**2 )

        fx(i) =
     &    - 3.375           * ( 9.0 * x(i) - 2.0 ) * t1
     &    - ( 27.0 / 98.0 ) * ( 9.0 * x(i) + 1.0 ) * t2
     &    - 2.25            * ( 9.0 * x(i) - 7.0 ) * t3
     &    + 3.6             * ( 9.0 * x(i) - 4.0 ) * t4

        fy(i) =
     &    - 3.375 * ( 9.0 * y(i) - 2.0 ) * t1
     &    - 0.675                        * t2
     &    - 2.25  * ( 9.0 * y(i) - 3.0 ) * t3
     &    + 3.6   * ( 9.0 * y(i) - 7.0 ) * t4

      end do

      return
      end
      subroutine f01_f2 ( n, x, y, fxx, fxy, fyy )

c*********************************************************************72
c
cc F01_F2 returns second derivatives of function f01.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    28 January 2012
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the number of evaluation points.
c
c    Input, double precision X(N), Y(N), the evalution points.
c
c    Output, double precision FXX(N), FXY(N), FYY(N), second derivatives.
c
      implicit none

      integer n

      double precision fxx(n)
      double precision fxy(n)
      double precision fyy(n)
      integer i
      double precision t1
      double precision t2
      double precision t3
      double precision t4
      double precision x(n)
      double precision y(n)  

      do i = 1, n

        t1 = exp ( - ( ( 9.0 * x(i) - 2.0 )**2
     &                  + ( 9.0 * y(i) - 2.0 )**2 ) / 4.0 )
        t2 = exp ( - ( ( 9.0 * x(i) + 1.0 )**2 ) / 49.0
     &                  - ( 9.0 * y(i) + 1.0 ) / 10.0 )
        t3 = exp ( - ( ( 9.0 * x(i) - 7.0 )**2
     &                  + ( 9.0 * y(i) - 3.0 )**2 ) / 4.0 )
        t4 = exp ( -   ( 9.0 * x(i) - 4.0 )**2
     &                  - ( 9.0 * y(i) - 7.0 )**2 )

        fxx(i) =
     &      15.1875 * ( ( 9.0 * x(i) - 2.0 )**2 - 2.0 )  * t1
     &    + 60.75   * ( ( 9.0 * x(i) + 1.0 )**2 - 24.5 ) * t2
     &    + 10.125  * ( ( 9.0 * x(i) - 7.0 )**2 - 2.0 )  * t3
     &    - 64.8    * ( ( 9.0 * x(i) - 4.0 )**2 - 0.5 )  * t4

        fxy(i) =
     &      15.1875 * ( 9.0 * x(i) - 2.0 ) * ( 9.0 * y(i) - 2.0 ) * t1
     &    + ( 243.0 / 980.0 ) * ( 9.0 * x(i) + 1.0 ) * t2
     &    + 10.125 * ( 9.0 * x(i) - 7.0 ) * ( 9.0 * y(i) - 3.0 ) * t3
     &    - 64.8 * ( 9.0 * x(i) - 4.0 ) * ( 9.0 * y(i) - 7.0 ) * t4

        fyy(i) =
     &      15.1875 * ( ( 9.0 * y(i) - 2.0 )**2 - 2.0 ) * t1
     &    + 0.6075  *                                     t2
     &    + 10.125  * ( ( 9.0 * y(i) - 3.0 )**2 - 2.0 ) * t3
     &    - 64.8    * ( ( 9.0 * y(i) - 7.0 )**2 - 0.5 ) * t4

      end do

      return
      end
      subroutine f01_title ( ft )

c*********************************************************************72
c
cc F01_TITLE returns the title for function f01.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    28 January 2012
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, character * ( * ) FT, the function title.
c
      implicit none

      character * ( * ) ft

      ft = 'Exponential'

      return
      end
      subroutine f02_f0 ( n, x, y, f )

c*********************************************************************72
c
cc F02_F0 returns the value of function 2.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    28 January 2012
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the number of evaluation points.
c
c    Input, double precision X(N), Y(N), the evalution points.
c
c    Output, double precision F(N), the function values.
c
      implicit none

      integer n

      double precision f(n)
      integer i
      double precision x(n)
      double precision y(n)  

      do i = 1, n 
        f(i) = ( tanh ( 9.0 * ( y(i) - x(i) ) ) + 1.0 ) / 9.0
      end do

      return
      end
      subroutine f02_f1 ( n, x, y, fx, fy )

c*********************************************************************72
c
cc F02_F1 returns first derivatives of function 2.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    28 January 2012
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the number of evaluation points.
c
c    Input, double precision X(N), Y(N), the evalution points.
c
c    Output, double precision FX(N), FY(N), the derivative values.
c
      implicit none

      integer n

      double precision fx(n)
      double precision fy(n)
      integer i
      double precision t1
      double precision x(n)
      double precision y(n)  

      do i = 1, n
        t1 = 18.0 * ( y(i) - x(i) )
        fx(i) = - 4.0 / ( exp ( t1 ) + 2.0 + exp ( - t1 ) )
        fy(i) = - fx(i)
      end do

      return
      end
      subroutine f02_f2 ( n, x, y, fxx, fxy, fyy )

c*********************************************************************72
c
cc F02_F2 returns second derivatives of function 2.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    28 January 2012
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the number of evaluation points.
c
c    Input, double precision X(N), Y(N), the evalution points.
c
c    Output, double precision FXX(N), FXY(N), FYY(N), second derivatives.
c
      implicit none

      integer n

      double precision fxx(n)
      double precision fxy(n)
      double precision fyy(n)
      integer i
      double precision t1
      double precision x(n)
      double precision y(n)  

      do i = 1, n
        t1 = 18.0 * ( y(i) - x(i) )

        fxx(i) = 18.0 * tanh ( 0.5 * t1 )
     &    * ( tanh ( 9.0 * ( y(i) - x(i) ) ) + 1.0 ) / 9.0
        fxy(i) = - fxx(i)
        fyy(i) = fxx(i)
      end do

      return
      end
      subroutine f02_title ( ft )

c*********************************************************************72
c
cc F02_TITLE returns the title for function 2.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    28 January 2012
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, character * ( * ) FT, the function title.
c
      implicit none

      character * ( * ) ft

      ft = 'Cliff'

      return
      end
      subroutine f03_f0 ( n, x, y, f )

c*********************************************************************72
c
cc F03_F0 returns the value of function 3.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    28 January 2012
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the number of evaluation points.
c
c    Input, double precision X(N), Y(N), the evalution points.
c
c    Output, double precision F(N), the function values.
c
      implicit none

      integer n

      double precision f(n)
      integer i
      double precision x(n)
      double precision y(n)  

      do i = 1, n
        f(i) = ( 1.25 + cos ( 5.4 * y(i) ) )
     &    / ( 6.0 + 6.0 * ( 3.0 * x(i) - 1.0 )**2 )
      end do

      return
      end
      subroutine f03_f1 ( n, x, y, fx, fy )

c*********************************************************************72
c
cc F03_F1 returns first derivatives of function 3.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    28 January 2012
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the number of evaluation points.
c
c    Input, double precision X(N), Y(N), the evalution points.
c
c    Output, double precision FX(N), FY(N), the derivative values.
c
      implicit none

      integer n

      double precision fx(n)
      double precision fy(n)
      integer i
      double precision t1
      double precision t2
      double precision x(n)
      double precision y(n)  

      do i = 1, n
        t1 = 5.4 * y(i)
        t2 = 1.0 + ( 3.0 * x(i) - 1.0 )**2
        fx(i) = - ( 3.0 * x(i) - 1.0 )
     &    * ( 1.25 + cos ( t1 ) ) / ( t2 * t2 )
        fy(i) = - 0.9 * sin ( t1 ) / t2
      end do

      return
      end
      subroutine f03_f2 ( n, x, y, fxx, fxy, fyy )

c*********************************************************************72
c
cc F03_F2 returns second derivatives of function 3.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    28 January 2012
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the number of evaluation points.
c
c    Input, double precision X(N), Y(N), the evalution points.
c
c    Output, double precision FXX(N), FXY(N), FYY(N), second derivatives.
c
      implicit none

      integer n

      double precision fxx(n)
      double precision fxy(n)
      double precision fyy(n)
      integer i
      double precision t1
      double precision t2
      double precision x(n)
      double precision y(n)  

      do i = 1, n 
        t1 = 5.4 * y(i)
        t2 = 1.0 + ( 3.0 * x(i) - 1.0 )**2

        fxx(i) = 3.0 * ( 1.25 + cos ( t1 ) ) * ( 3.0 * t2 - 4.0 )
     &    / ( t2**3 )
        fxy(i) = 5.4 * ( 3.0 * x(i) - 1.0 ) * sin ( t1 )
     &    / ( t2 * t2 )
        fyy(i) = - 4.86 * cos ( t1 ) / t2
      end do

      return
      end
      subroutine f03_title ( ft )

c*********************************************************************72
c
cc F03_TITLE returns the title for function 3.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    28 January 2012
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, character * ( * ) FT, the function title.
c
      implicit none

      character * ( * ) ft

      ft = 'Saddle'

      return
      end
      subroutine f04_f0 ( n, x, y, f )

c*********************************************************************72
c
cc F04_F0 returns the value of function 4.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    28 January 2012
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the number of evaluation points.
c
c    Input, double precision X(N), Y(N), the evalution points.
c
c    Output, double precision F(N), the function values.
c
      implicit none

      integer n

      double precision f(n)
      integer i
      double precision x(n)
      double precision y(n)  

      do i = 1, n
        f(i) = exp ( - 5.0625 * ( ( x(i) - 0.5 )**2
     &                          + ( y(i) - 0.5 )**2 ) ) / 3.0
      end do

      return
      end
      subroutine f04_f1 ( n, x, y, fx, fy )

c*********************************************************************72
c
cc F04_F1 returns first derivatives of function 4.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    28 January 2012
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the number of evaluation points.
c
c    Input, double precision X(N), Y(N), the evalution points.
c
c    Output, double precision FX(N), FY(N), the derivative values.
c
      implicit none

      integer n

      double precision fx(n)
      double precision fy(n)
      integer i
      double precision t1
      double precision t2
      double precision t3
      double precision x(n)
      double precision y(n)  

      do i = 1, n
        t1 = x(i) - 0.5
        t2 = y(i) - 0.5
        t3 = - 3.375 * exp ( - 5.0625 * ( t1 * t1 + t2 * t2 ) )
        fx(i) = t1 * t3
        fy(i) = t2 * t3
      end do

      return
      end
      subroutine f04_f2 ( n, x, y, fxx, fxy, fyy )

c*********************************************************************72
c
cc F04_F2 returns second derivatives of function 4.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    28 January 2012
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the number of evaluation points.
c
c    Input, double precision X(N), Y(N), the evalution points.
c
c    Output, double precision FXX(N), FXY(N), FYY(N), second derivatives.
c
      implicit none

      integer n

      double precision fxx(n)
      double precision fxy(n)
      double precision fyy(n)
      integer i
      double precision t1
      double precision t2
      double precision t3
      double precision x(n)
      double precision y(n)  

      do i = 1, n
        t1 = x(i) - 0.5
        t2 = y(i) - 0.5
        t3 = - 3.375
     &    * exp ( - 5.0625 * ( t1 * t1 + t2 * t2 ) )

        fxx(i) = ( 1.0 - 10.125 * t1 * t1 ) * t3
        fxy(i) = - 10.125 * t1 * t2 * t3
        fyy(i) = ( 1.0 - 10.125 * t2 * t2 ) * t3
      end do

      return
      end
      subroutine f04_title ( ft )

c*********************************************************************72
c
cc F04_TITLE returns the title for function 4.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    28 January 2012
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, character * ( * ) FT, the function title.
c
      implicit none

      character * ( * ) ft

      ft = 'Gentle'

      return
      end
      subroutine f05_f0 ( n, x, y, f )

c*********************************************************************72
c
cc F05_F0 returns the value of function 5.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    28 January 2012
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the number of evaluation points.
c
c    Input, double precision X(N), Y(N), the evalution points.
c
c    Output, double precision F(N), the function values.
c
      implicit none

      integer n

      double precision f(n)
      integer i
      double precision x(n)
      double precision y(n)  

      do i = 1, n
        f(i) = exp ( - 20.25 * ( ( x(i) - 0.5 )**2 
     &                         + ( y(i) - 0.5 )**2 ) ) / 3.0
      end do

      return
      end
      subroutine f05_f1 ( n, x, y, fx, fy )

c*********************************************************************72
c
cc F05_F1 returns first derivatives of function 5.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    28 January 2012
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the number of evaluation points.
c
c    Input, double precision X(N), Y(N), the evalution points.
c
c    Output, double precision FX(N), FY(N), the derivative values.
c
      implicit none

      integer n

      double precision fx(n)
      double precision fy(n)
      integer i
      double precision t1
      double precision t2
      double precision t3
      double precision x(n)
      double precision y(n)  

      do i = 1, n
        t1 = x(i) - 0.5
        t2 = y(i) - 0.5
        t3 = - 13.5 * exp ( - 20.25 * ( t1 * t1 + t2 * t2 ) )
        fx(i) = t1 * t3
        fy(i) = t2 * t3
      end do

      return
      end
      subroutine f05_f2 ( n, x, y, fxx, fxy, fyy )

c*********************************************************************72
c
cc F05_F2 returns second derivatives of function 5.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    28 January 2012
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the number of evaluation points.
c
c    Input, double precision X(N), Y(N), the evalution points.
c
c    Output, double precision FXX(N), FXY(N), FYY(N), second derivatives.
c
      implicit none

      integer n

      double precision fxx(n)
      double precision fxy(n)
      double precision fyy(n)
      integer i
      double precision t1
      double precision t2
      double precision t3
      double precision x(n)
      double precision y(n)  

      do i = 1, n
        t1 = x(i) - 0.5
        t2 = y(i) - 0.5
        t3 = - 13.5 * exp ( - 20.25 * ( t1 * t1 + t2 * t2 ) )

        fxx(i) = ( 1.0 - 40.5 * t1 * t1 ) * t3
        fxy(i) = - 40.5 * t1 * t2 * t3
        fyy(i) = ( 1.0 - 40.5 * t2 * t2 ) * t3

      end do

      return
      end
      subroutine f05_title ( ft )

c*********************************************************************72
c
cc F05_TITLE returns the title for function 5.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    28 January 2012
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, character * ( * ) FT, the function title.
c
      implicit none

      character * ( * ) ft

      ft = 'Steep'

      return
      end
      subroutine f06_f0 ( n, x, y, f )

c*********************************************************************72
c
cc F06_F0 returns the value of function 6.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    28 January 2012
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the number of evaluation points.
c
c    Input, double precision X(N), Y(N), the evalution points.
c
c    Output, double precision F(N), the function values.
c
      implicit none

      integer n

      double precision f(n)
      integer i
      double precision t4
      double precision x(n)
      double precision y(n)  

      do i = 1, n

        t4 = 64.0 - 81.0 * ( ( x(i) - 0.5 )**2 + ( y(i) - 0.5 )**2 )

        if ( 0.0 .le. t4 ) then
          f(i) = sqrt ( t4 ) / 9.0 - 0.5
        else
          f(i) = 0.0
        end if

      end do

      return
      end
      subroutine f06_f1 ( n, x, y, fx, fy )

c*********************************************************************72
c
cc F06_F1 returns first derivatives of function 6.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    28 January 2012
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the number of evaluation points.
c
c    Input, double precision X(N), Y(N), the evalution points.
c
c    Output, double precision FX(N), FY(N), the derivative values.
c
      implicit none

      integer n

      double precision fx(n)
      double precision fy(n)
      integer i
      double precision t1
      double precision t2
      double precision t3
      double precision t4
      double precision x(n)
      double precision y(n)  

      do i = 1, n

        t4 = 64.0 - 81.0 * ( ( x(i) - 0.5 )**2 + ( y(i) - 0.5 )**2 )

        if ( 0.0 .lt. t4 ) then
          t1 = x(i) - 0.5
          t2 = y(i) - 0.5
          t3 = - 9.0 / sqrt ( t4 )
          fx(i) = t1 * t3
          fy(i) = t2 * t3
        else
          fx(i) = 0.0
          fy(i) = 0.0
        end if

      end do

      return
      end
      subroutine f06_f2 ( n, x, y, fxx, fxy, fyy )

c*********************************************************************72
c
cc F06_F2 returns second derivatives of function 6.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    28 January 2012
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the number of evaluation points.
c
c    Input, double precision X(N), Y(N), the evalution points.
c
c    Output, double precision FXX(N), FXY(N), FYY(N), second derivatives.
c
      implicit none

      integer n

      double precision fxx(n)
      double precision fxy(n)
      double precision fyy(n)
      integer i
      double precision t1
      double precision t2
      double precision t3
      double precision t4
      double precision x(n)
      double precision y(n)  

      do i = 1, n

        t4 = 64.0 - 81.0 * ( ( x(i) - 0.5 )**2 + ( y(i) - 0.5 )**2 )

        if ( 0.0 .lt. t4 ) then
          t1 = x(i) - 0.5
          t2 = y(i) - 0.5
          t3 = - 9.0 / sqrt ( t4 )
          fxx(i) = ( 1.0 + t1 * t3 * t1 * t3 ) * t3
          fxy(i) = t1 * t3 * t2 * t3
          fyy(i) = ( 1.0 + t2 * t3 * t2 * t3 ) * t3
        else
          fxx(i) = 0.0
          fxy(i) = 0.0
          fyy(i) = 0.0
        end if

      end do

      return
      end
      subroutine f06_title ( ft )

c*********************************************************************72
c
cc F06_TITLE returns the title for function 6.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    28 January 2012
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, character * ( * ) FT, the function title.
c
      implicit none

      character * ( * ) ft

      ft = 'Sphere'

      return
      end
      subroutine f07_f0 ( n, x, y, f )

c*********************************************************************72
c
cc F07_F0 returns the value of function 7.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    28 January 2012
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the number of evaluation points.
c
c    Input, double precision X(N), Y(N), the evalution points.
c
c    Output, double precision F(N), the function values.
c
      implicit none

      integer n

      double precision f(n)
      integer i
      double precision x(n)
      double precision y(n)  

      do i = 1, n
        f(i) = 2.0 * cos ( 10.0 * x(i) ) * sin ( 10.0 * y(i) )
     &    + sin ( 10.0 * x(i) * y(i) )
      end do

      return
      end
      subroutine f07_f1 ( n, x, y, fx, fy )

c*********************************************************************72
c
cc F07_F1 returns first derivatives of function 7.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    28 January 2012
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the number of evaluation points.
c
c    Input, double precision X(N), Y(N), the evalution points.
c
c    Output, double precision FX(N), FY(N), the derivative values.
c
      implicit none

      integer n

      double precision fx(n)
      double precision fy(n)
      integer i
      double precision t1
      double precision t2
      double precision t3
      double precision x(n)
      double precision y(n)  

      do i = 1, n
        t1 = 10.0 * x(i)
        t2 = 10.0 * y(i)
        t3 = 10.0 * cos ( 10.0 * x(i) * y(i) )
        fx(i) = - 20.0 * sin ( t1 ) * sin ( t2 ) + t3 * y(i)
        fy(i) =   20.0 * cos ( t1 ) * cos ( t2 ) + t3 * x(i)
      end do

      return
      end
      subroutine f07_f2 ( n, x, y, fxx, fxy, fyy )

c*********************************************************************72
c
cc F07_F2 returns second derivatives of function 7.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    28 January 2012
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the number of evaluation points.
c
c    Input, double precision X(N), Y(N), the evalution points.
c
c    Output, double precision FXX(N), FXY(N), FYY(N), second derivatives.
c
      implicit none

      integer n

      double precision fxx(n)
      double precision fxy(n)
      double precision fyy(n)
      integer i
      double precision t1
      double precision t2
      double precision t3
      double precision t4
      double precision x(n)
      double precision y(n)  

      do i = 1, n
        t1 = 10.0 * x(i)
        t2 = 10.0 * y(i)
        t3 = 10.0 * cos ( 10.0 * x(i) * y(i) )
        t4 = 100.0 * sin ( 10.0 * x(i) * y(i) )

        fxx(i) = - 200.0 * cos ( t1 ) * sin ( t2 )
     &    - t4 * y(i) * y(i)
        fxy(i) = - 200.0 * sin ( t1 ) * cos ( t2 )
     &    + t3 - t4 * x(i) * y(i)
        fyy(i) = - 200.0 * cos ( t1 ) * sin ( t2 )
     &    - t4 * x(i) * x(i)
      end do

      return
      end
      subroutine f07_title ( ft )

c*********************************************************************72
c
cc F07_TITLE returns the title for function 7.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    28 January 2012
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, character * ( * ) FT, the function title.
c
      implicit none

      character * ( * ) ft

      ft = 'Trig'

      return
      end
      subroutine f08_f0 ( n, x, y, f )

c*********************************************************************72
c
cc F08_F0 returns the value of function 8.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    28 January 2012
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the number of evaluation points.
c
c    Input, double precision X(N), Y(N), the evalution points.
c
c    Output, double precision F(N), the function values.
c
      implicit none

      integer n

      double precision f(n)
      integer i
      double precision t1
      double precision t2
      double precision t3
      double precision t4
      double precision x(n)
      double precision y(n)  

      do i = 1, n
        t1 = 5.0 - 10.0 * x(i)
        t2 = 5.0 - 10.0 * y(i)
        t3 = exp ( - 0.5 * t1 * t1 )
        t4 = exp ( - 0.5 * t2 * t2 )
        f(i) = t3 + 0.75 * t4 * ( 1.0 + t3 )
      end do

      return
      end
      subroutine f08_f1 ( n, x, y, fx, fy )

c*********************************************************************72
c
cc F08_F1 returns first derivatives of function 8.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    28 January 2012
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the number of evaluation points.
c
c    Input, double precision X(N), Y(N), the evalution points.
c
c    Output, double precision FX(N), FY(N), the derivative values.
c
      implicit none

      integer n

      double precision fx(n)
      double precision fy(n)
      integer i
      double precision t1
      double precision t2
      double precision t3
      double precision t4
      double precision x(n)
      double precision y(n)  

      do i = 1, n
        t1 = 5.0 - 10.0 * x(i)
        t2 = 5.0 - 10.0 * y(i)
        t3 = exp ( - 0.5 * t1 * t1 )
        t4 = exp ( - 0.5 * t2 * t2 )

        fx(i) = t1 * t3 * ( 10.0 + 7.5 * t4 )
        fy(i) = t2 * t4 * ( 7.5 + 7.5 * t3 )
      end do

      return
      end
      subroutine f08_f2 ( n, x, y, fxx, fxy, fyy )

c*********************************************************************72
c
cc F08_F2 returns second derivatives of function 8.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    28 January 2012
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the number of evaluation points.
c
c    Input, double precision X(N), Y(N), the evalution points.
c
c    Output, double precision FXX(N), FXY(N), FYY(N), second derivatives.
c
      implicit none

      integer n

      double precision fxx(n)
      double precision fxy(n)
      double precision fyy(n)
      integer i
      double precision t1
      double precision t2
      double precision t3
      double precision t4
      double precision x(n)
      double precision y(n)  

      do i = 1, n
        t1 = 5.0 - 10.0 * x(i)
        t2 = 5.0 - 10.0 * y(i)
        t3 = exp ( - 0.5 * t1 * t1 )
        t4 = exp ( - 0.5 * t2 * t2 )

        fxx(i) = t3 * ( t1 * t1 - 1.0 ) * ( 100.0 + 75.0 * t4 )
        fxy(i) = 75.0 * t1 * t2 * t3 * t4
        fyy(i) = t4 * ( t2 * t2 - 1.0 ) * ( 75.0 + 75.0 * t3 )
      end do

      return
      end
      subroutine f08_title ( ft )

c*********************************************************************72
c
cc F08_TITLE returns the title for function 8.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    28 January 2012
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, character * ( * ) FT, the function title.
c
      implicit none

      character * ( * ) ft

      ft = 'Synergistic Gaussian'

      return
      end
      subroutine f09_f0 ( n, x, y, f )

c*********************************************************************72
c
cc F09_F0 returns the value of function 9.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    28 January 2012
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the number of evaluation points.
c
c    Input, double precision X(N), Y(N), the evalution points.
c
c    Output, double precision F(N), the function values.
c
      implicit none

      integer n

      double precision f(n)
      integer i
      double precision t1
      double precision t2
      double precision t3
      double precision t4
      double precision x(n)
      double precision y(n)  

      do i = 1, n 
        t1 = exp ( ( 10.0 - 20.0 * x(i) ) / 3.0 )
        t2 = exp ( ( 10.0 - 20.0 * y(i) ) / 3.0 )
        t3 = 1.0 / ( 1.0 + t1 )
        t4 = 1.0 / ( 1.0 + t2 )
        f(i) = ( ( 20.0 / 3.0 )**3 * t1 * t2 )**2
     &    * ( t3 * t4 )**5
     &    * ( t1 - 2.0 * t3 ) * ( t2 - 2.0 * t4 )
      end do

      return
      end
      subroutine f09_f1 ( n, x, y, fx, fy )

c*********************************************************************72
c
cc F09_F1 returns first derivatives of function 9.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    28 January 2012
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the number of evaluation points.
c
c    Input, double precision X(N), Y(N), the evalution points.
c
c    Output, double precision FX(N), FY(N), the derivative values.
c
      implicit none

      integer n

      double precision fx(n)
      double precision fy(n)
      integer i
      double precision t1
      double precision t2
      double precision t3
      double precision t4
      double precision x(n)
      double precision y(n)  

      do i = 1, n
        t1 = exp ( ( 10.0 - 20.0 * x(i) ) / 3.0 )
        t2 = exp ( ( 10.0 - 20.0 * y(i) ) / 3.0 )
        t3 = 1.0 / ( 1.0 + t1 )
        t4 = 1.0 / ( 1.0 + t2 )

        fx(i) = ( ( 20.0 / 3.0 ) * t1 )**2 * ( ( 20.0 / 3.0 ) * t3 )**5
     &    * ( 2.0 * t1 - 3.0 * t3 - 5.0 + 12.0 * t3 * t3 )
     &    * t2 * t2 * t4**5
     &    * ( t2 - 2.0 * t4 )

        fy(i) = ( ( 20.0 / 3.0 ) * t1 )**2 * ( ( 20.0 / 3.0 ) * t3 )**5
     &    * ( 2.0 * t2 - 3.0 * t4 - 5.0 + 12.0 * t4 * t4 )
     &    * t2 * t2 * t4**5
     &    * ( t1 - 2.0 * t3 )
      end do

      return
      end
      subroutine f09_f2 ( n, x, y, fxx, fxy, fyy )

c*********************************************************************72
c
cc F09_F2 returns second derivatives of function 9.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    28 January 2012
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the number of evaluation points.
c
c    Input, double precision X(N), Y(N), the evalution points.
c
c    Output, double precision FXX(N), FXY(N), FYY(N), second derivatives.
c
      implicit none

      integer n

      double precision fxx(n)
      double precision fxy(n)
      double precision fyy(n)
      integer i
      double precision t1
      double precision t2
      double precision t3
      double precision t4
      double precision t5
      double precision t6
      double precision x(n)
      double precision y(n)  

      do i = 1, n
        t1 = exp ( ( 10.0 - 20.0 * x(i) ) / 3.0 )
        t2 = exp ( ( 10.0 - 20.0 * y(i) ) / 3.0 )
        t3 = 1.0 / ( 1.0 + t1 )
        t4 = 1.0 / ( 1.0 + t2 )
        t5 = 20.0 / 3.0
        t6 = ( t5 * t1 * t2 )**2 * ( t5 * t3 * t4 )**5

        fxx(i) = t5 * t6 * ( t2 - 2.0 * t4 )
     &    * ( ( ( - 84.0 * t3 + 78.0 ) * t3 + 23.0 ) * t3
     &    + 4.0 * t1 - 25.0 )

        fxy(i) = t5 * t6
     &    * ( ( 12.0 * t4 - 3.0 ) * t4 + 2.0 * t2 - 5.0 )
     &    * ( ( 12.0 * t3 - 3.0 ) * t3 + 2.0 * t1 - 5.0 )

        fyy(i) = t5 * t6 * ( t1 - 2.0 * t3 )
     &    * ( ( ( - 84.0 * t4 + 78.0 ) * t4 + 23.0 ) * t4
     &    + 4.0 * t2 - 25.0 )
      end do

      return
      end
      subroutine f09_title ( ft )

c*********************************************************************72
c
cc F09_TITLE returns the title for function 9.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    28 January 2012
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, character * ( * ) FT, the function title.
c
      implicit none

      character * ( * ) ft

      ft = 'Cloverleaf Asymmetric Peak/Valley'

      return
      end
      subroutine f10_f0 ( n, x, y, f )

c*********************************************************************72
c
cc F10_F0 returns the value of function f10.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    28 January 2012
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the number of evaluation points.
c
c    Input, double precision X(N), Y(N), the evalution points.
c
c    Output, double precision F(N), the function values.
c
      implicit none

      integer n

      double precision f(n)
      integer i
      double precision t1
      double precision t2
      double precision t3
      double precision x(n)
      double precision y(n)  

      do i = 1, n
        t1 = sqrt ( ( 80.0 * x(i) - 40.0 )**2 
     &            + ( 90.0 * y(i) - 45.0 )**2 )
        t2 = exp ( - 0.04 * t1 )
        t3 = cos ( 0.15 * t1 )

        f(i) = t2 * t3
      end do

      return
      end
      subroutine f10_f1 ( n, x, y, fx, fy )

c*********************************************************************72
c
cc F10_F1 returns first derivatives of function f10.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    28 January 2012
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the number of evaluation points.
c
c    Input, double precision X(N), Y(N), the evalution points.
c
c    Output, double precision FX(N), FY(N), the derivative values.
c
      implicit none

      integer n

      double precision fx(n)
      double precision fy(n)
      integer i
      double precision t1
      double precision t2
      double precision t3
      double precision t4
      double precision x(n)
      double precision y(n)  

      do i = 1, n

        t1 = sqrt ( ( 80.0 * x(i) - 40.0 )**2 
     &            + ( 90.0 * y(i) - 45.0 )**2 )

        if ( t1 .eq. 0.0 ) then
          fx(i) = 0.0
          fy(i) = 0.0
        else
          t2 = exp ( - 0.04 * t1 )
          t3 = cos ( 0.15 * t1 )
          t4 = sin ( 0.15 * t1 )
          fx(i) = - t2 * ( 12.0 * t4 + 3.2 * t3 ) 
     &      * ( 80.0 * x(i) - 40.0 ) / t1
          fy(i) = - t2 * ( 13.5 * t4 + 3.6 * t3 ) 
     &      * ( 90.0 * y(i) - 45.0 ) / t1
        end if

      end do

      return
      end
      subroutine f10_f2 ( n, x, y, fxx, fxy, fyy )

c*********************************************************************72
c
cc F10_F2 returns second derivatives of function f10.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    28 January 2012
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the number of evaluation points.
c
c    Input, double precision X(N), Y(N), the evalution points.
c
c    Output, double precision FXX(N), FXY(N), FYY(N), second derivatives.
c
      implicit none

      integer n

      double precision fxx(n)
      double precision fxy(n)
      double precision fyy(n)
      integer i
      double precision t1
      double precision t2
      double precision t3
      double precision t4
      double precision t5
      double precision x(n)
      double precision y(n)  

      do i = 1, n

        t1 = sqrt ( ( 80.0 * x(i) - 40.0 )**2 
     &            + ( 90.0 * y(i) - 45.0 )**2 )

        if ( t1 .eq. 0.0 ) then
          fxx(i) = 0.0
          fxy(i) = 0.0
          fyy(i) = 0.0
        else
          t2 = exp ( - 0.04 * t1 )
          t3 = cos ( 0.15 * t1 )
          t4 = sin ( 0.15 * t1 )
          t5 = t2 / t1**3

          fxx(i) = t5 * ( t1 * ( 76.8 * t4 - 133.76 * t3 )
     &      * ( 80.0 * x(i) - 40.0 )**2
     &      - ( 960.0 * t4 + 256.0 * t3 ) 
     &      * ( 90.0 * y(i) - 45.0 )**2 )

          fxy(i) = t5 * ( t1 * ( 86.4 * t4 - 150.48 * t3 ) 
     &      + 1080.0 * t4 + 288.0 * t3 ) * ( 80.0 * x(i) - 40.0 ) 
     &      * ( 90.0 * y(i) - 45.0 )

          fyy(i) = t5 * ( t1 * ( 97.2 * t4 - 169.29 * t3 )
     &      * ( 90.0 * y(i) - 45.0 )**2 - ( 1215.0 * t4 + 324.0 * t3 ) 
     &      * ( 80.0 * x(i) - 40.0 )**2 )

        end if

      end do

      return
      end
      subroutine f10_title ( ft )

c*********************************************************************72
c
cc F10_TITLE returns the title for function f10.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    28 January 2012
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, character * ( * ) FT, the function title.
c
      implicit none

      character * ( * ) ft

      ft = 'Cosine Peak'

      return
      end
      subroutine f11_f0 ( n, x, y, f )

c*********************************************************************72
c
cc F11_F0 returns the value of function f11.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    03 August 2012
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the number of evaluation points.
c
c    Input, double precision X(N), Y(N), the evalution points.
c
c    Output, double precision F(N), the function values.
c
      implicit none

      integer n

      double precision f(n)
      integer i
      double precision x(n)
      double precision y(n)  

      do i = 1, n
        f(i) = x(i) * ( y(i) + 1.0D+00 )
      end do

      return
      end
      subroutine f11_f1 ( n, x, y, fx, fy )

c*********************************************************************72
c
cc F11_F1 returns first derivatives of function f11.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    03 August 2012
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the number of evaluation points.
c
c    Input, double precision X(N), Y(N), the evalution points.
c
c    Output, double precision FX(N), FY(N), the derivative values.
c
      implicit none

      integer n

      double precision fx(n)
      double precision fy(n)
      integer i
      double precision x(n)
      double precision y(n)  

      do i = 1, n

        fx(i) = y(i) + 1.0D+00
        fy(i) = x(i)

      end do

      return
      end
      subroutine f11_f2 ( n, x, y, fxx, fxy, fyy )

c*********************************************************************72
c
cc F11_F2 returns second derivatives of function f11.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    03 August 2012
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the number of evaluation points.
c
c    Input, double precision X(N), Y(N), the evalution points.
c
c    Output, double precision FXX(N), FXY(N), FYY(N), second derivatives.
c
      implicit none

      integer n

      double precision fxx(n)
      double precision fxy(n)
      double precision fyy(n)
      integer i
      double precision x(n)
      double precision y(n)  

      do i = 1, n

        fxx(i) = 0.0D+00
        fxy(i) = 1.0D+00
        fyy(i) = 0.0D+00

      end do

      return
      end
      subroutine f11_title ( ft )

c*********************************************************************72
c
cc F11_TITLE returns the title for function f11.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    03 August 2012
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, character * ( * ) FT, the function title.
c
      implicit none

      character * ( * ) ft

      ft = 'Bilinear function'

      return
      end
      subroutine f12_f0 ( n, x, y, f )

c*********************************************************************72
c
cc F12_F0 returns the value of function f12.
c
c  Discussion:
c
c    This is an example from Vicente Romero.
c
c      R = sqrt ( X^2 + Y^2 )
c      T = atan ( Y / X )
c      F(X,Y) = ( 0.8 * R + 0.35 * sin ( 2.4 * pi * R / sqrt ( 2 )  ) )
c               * 1.5 * sin ( 1.3 * T )
c
c    The mean and standard deviation of the function over the interval
c    are approximately:
c
c      mu    = 0.581608
c      sigma = 0.343208
c
c    Since the interpolation interval is the unit square, this means the
c    integral of the function over the interval can also be estimated as
c
c      I = 0.581608
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    03 August 2012
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the number of evaluation points.
c
c    Input, double precision X(N), Y(N), the evalution points.
c
c    Output, double precision F(N), the function values.
c
      implicit none

      integer n

      double precision f(n)
      integer i
      double precision pi
      parameter ( pi = 3.141592653589793D+00 )
      double precision r
      double precision t
      double precision x(n)
      double precision y(n)  

      do i = 1, n

        r = sqrt ( x(i)**2 + y(i)**2 )
        t = atan2 ( y(i), x(i) )

        f(i) = 1.5D+00 * ( 0.8D+00 * r 
     &    + 0.35D+00 * sin ( 2.4D+00 * pi * r / sqrt ( 2.0D+00 ) ) ) 
     &    * sin ( 1.3D+00 * t )

      end do

      return
      end
      subroutine f12_f1 ( n, x, y, fx, fy )

c*********************************************************************72
c
cc F12_F1 returns first derivatives of function f12.
c
c  Discussion:
c
c    Currently, the derivative information is of no interest to me.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    03 August 2012
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the number of evaluation points.
c
c    Input, double precision X(N), Y(N), the evalution points.
c
c    Output, double precision FX(N), FY(N), the derivative values.
c
      implicit none

      integer n

      double precision fx(n)
      double precision fy(n)
      integer i
      double precision x(n)
      double precision y(n)  

      do i = 1, n

        fx(i) = 0.0D+00
        fy(i) = 0.0D+00

      end do

      return
      end
      subroutine f12_f2 ( n, x, y, fxx, fxy, fyy )

c*********************************************************************72
c
cc F12_F2 returns second derivatives of function f12.
c
c  Discussion:
c
c    Currently, the derivative information is of no interest to me.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    03 August 2012
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the number of evaluation points.
c
c    Input, double precision X(N), Y(N), the evalution points.
c
c    Output, double precision FXX(N), FXY(N), FYY(N), second derivatives.
c
      implicit none

      integer n

      double precision fxx(n)
      double precision fxy(n)
      double precision fyy(n)
      integer i
      double precision x(n)
      double precision y(n)  

      do i = 1, n

        fxx(i) = 0.0D+00
        fxy(i) = 0.0D+00
        fyy(i) = 0.0D+00

      end do

      return
      end
      subroutine f12_title ( ft )

c*********************************************************************72
c
cc F12_TITLE returns the title for function f12.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    03 August 2012
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, character * ( * ) FT, the function title.
c
      implicit none

      character * ( * ) ft

      ft = 'Vicente Romero function'

      return
      end
      subroutine f13_f0 ( n, x, y, f )

c*********************************************************************72
c
cc F13_F0 returns the value of function f13.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    03 August 2012
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the number of evaluation points.
c
c    Input, double precision X(N), Y(N), the evalution points.
c
c    Output, double precision F(N), the function values.
c
      implicit none

      integer n

      double precision f(n)
      integer i
      double precision x(n)
      double precision y(n)  

      do i = 1, n

        f(i) = 1.0D+00 / ( ( 10.0D+00 * x(i) - 5.0D+00 )**2 
     &                   + ( 10.0D+00 * y(i) - 5.0D+00 )**2 + 1.0D+00 )

      end do

      return
      end
      subroutine f13_f1 ( n, x, y, fx, fy )

c*********************************************************************72
c
cc F13_F1 returns first derivatives of function f13.
c
c  Discussion:
c
c    Currently, the derivative information is of no interest to me.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    03 August 2012
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the number of evaluation points.
c
c    Input, double precision X(N), Y(N), the evalution points.
c
c    Output, double precision FX(N), FY(N), the derivative values.
c
      implicit none

      integer n

      double precision fx(n)
      double precision fy(n)
      integer i
      double precision x(n)
      double precision y(n)  

      do i = 1, n

        fx(i) = 0.0D+00
        fy(i) = 0.0D+00

      end do

      return
      end
      subroutine f13_f2 ( n, x, y, fxx, fxy, fyy )

c*********************************************************************72
c
cc F13_F2 returns second derivatives of function f13.
c
c  Discussion:
c
c    Currently, the derivative information is of no interest to me.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    03 August 2012
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the number of evaluation points.
c
c    Input, double precision X(N), Y(N), the evalution points.
c
c    Output, double precision FXX(N), FXY(N), FYY(N), second derivatives.
c
      implicit none

      integer n

      double precision fxx(n)
      double precision fxy(n)
      double precision fyy(n)
      integer i
      double precision x(n)
      double precision y(n)  

      do i = 1, n

        fxx(i) = 0.0D+00
        fxy(i) = 0.0D+00
        fyy(i) = 0.0D+00

      end do

      return
      end
      subroutine f13_title ( ft )

c*********************************************************************72
c
cc F13_TITLE returns the title for function f13.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    03 August 2012
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, character * ( * ) FT, the function title.
c
      implicit none

      character * ( * ) ft

      ft = 'Rescaled Runge function'

      return
      end
      subroutine g00_num ( grid_num )

c*********************************************************************72
c
cc G00_NUM returns the number of grids available.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    28 January 2012
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c   Output, integer GRID_NUM, the number of grids.
c
      implicit none

      integer grid_num

      grid_num = 5

      return
      end
      subroutine g00_size ( gi, gn )

c*********************************************************************72
c
cc G00_SIZE returns the size for any grid.
c
c  Discussion:
c
c    The "grid size" is simply the number of points in the grid.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    28 January 2012
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer GI, the index of the grid.
c
c    Output, integer GN, the grid size.
c
      implicit none

      integer gi
      integer gn

      if ( gi .eq. 1 ) then
        call g01_size ( gn )
      else if ( gi .eq. 2 ) then
        call g02_size ( gn )
      else if ( gi .eq. 3 ) then
        call g03_size ( gn )
      else if ( gi .eq. 4 ) then
        call g04_size ( gn )
      else if ( gi .eq. 5 ) then
        call g05_size ( gn )
      else
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'G00_SIZE - Fatal error!'
        write ( *, '(a,i8)' ) '  Illegal grid index GI = ', gi
        stop
      end if

      return
      end
      subroutine g00_title ( gi, gt )

c*********************************************************************72
c
cc G00_TITLE returns the title for any grid.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    28 January 2012
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer GI, the index of the grid.
c
c    Output, character * ( * ) GT, the grid title.
c
      implicit none

      integer gi
      character * ( * ) gt

      if ( gi .eq. 1 ) then
        call g01_title ( gt )
      else if ( gi .eq. 2 ) then
        call g02_title ( gt )
      else if ( gi .eq. 3 ) then
        call g03_title ( gt )
      else if ( gi .eq. 4 ) then
        call g04_title ( gt )
      else if ( gi .eq. 5 ) then
        call g05_title ( gt )
      else
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'G00_TITLE - Fatal error!'
        write ( *, '(a,i8)' ) '  Illegal grid index GI = ', gi
        stop
      end if

      return
      end
      subroutine g00_xy ( gi, gn, gx, gy )

c*********************************************************************72
c
cc G00_XY returns the grid points for any grid.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    28 January 2012
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer GI, the index of the grid.
c
c    Input, integer GN, the grid size.
c
c    Output, double precision GX(GN), GY(GN), the grid coordinates.
c
      implicit none

      integer gn

      integer gi
      double precision gx(gn)
      double precision gy(gn)

      if ( gi .eq. 1 ) then
        call g01_xy ( gn, gx, gy )
      else if ( gi .eq. 2 ) then
        call g02_xy ( gn, gx, gy )
      else if ( gi .eq. 3 ) then
        call g03_xy ( gn, gx, gy )
      else if ( gi .eq. 4 ) then
        call g04_xy ( gn, gx, gy )
      else if ( gi .eq. 5 ) then
        call g05_xy ( gn, gx, gy )
      else
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'G00_XY - Fatal error!'
        write ( *, '(a,i8)' ) '  Illegal grid index GI = ', gi
        stop
      end if

      return
      end
      subroutine g01_size ( gn )

c*********************************************************************72
c
cc G01_SIZE returns the size for grid 1.
c
c  Discussion:
c
c    The "grid size" is simply the number of points in the grid.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    28 January 2012
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, integer GN, the grid size.
c
      implicit none

      integer gn

      gn = 100

      return
      end
      subroutine g01_title ( gt )

c*********************************************************************72
c
cc G01_TITLE returns the title for grid 1.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    28 January 2012
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, character * ( * ) GT, the grid title.
c
      implicit none

      character * ( * ) gt

      gt = 'Franke''s 100 node set'

      return
      end
      subroutine g01_xy ( gn, gx, gy )

c*********************************************************************72
c
cc G01_XY returns the grid points for grid 1.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    28 January 2012
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer GN, the grid size.
c
c    Output, double precision GX(GN), GY(GN), the grid coordinates.
c
      implicit none

      integer gn

      double precision gx(gn)
      double precision gy(gn)
      integer i
      double precision x(100)
      double precision y(100)

      save x
      save y

      data x /
     &   0.0227035,  0.0539888,  0.0217008,  0.0175129,  0.0019029,
     &  -0.0509685,  0.0395408, -0.0487061,  0.0315828, -0.0418785,
     &   0.1324189,  0.1090271,  0.1254439,  0.0934540,  0.0767578,
     &   0.1451874,  0.0626494,  0.1452734,  0.0958668,  0.0695559,
     &   0.2645602,  0.2391645,  0.2088990,  0.2767329,  0.1714726,
     &   0.2266781,  0.1909212,  0.1867647,  0.2304634,  0.2426219,
     &   0.3663168,  0.3857662,  0.3832392,  0.3179087,  0.3466321,
     &   0.3776591,  0.3873159,  0.3812917,  0.3795364,  0.2803515,
     &   0.4149771,  0.4277679,  0.4200010,  0.4663631,  0.4855658,
     &   0.4092026,  0.4792578,  0.4812279,  0.3977761,  0.4027321,
     &   0.5848691,  0.5730076,  0.6063893,  0.5013894,  0.5741311,
     &   0.6106955,  0.5990105,  0.5380621,  0.6096967,  0.5026188,
     &   0.6616928,  0.6427836,  0.6396475,  0.6703963,  0.7001181,
     &   0.6333590,  0.6908947,  0.6895638,  0.6718889,  0.6837675,
     &   0.7736939,  0.7635332,  0.7410424,  0.8258981,  0.7306034,
     &   0.8086609,  0.8214531,  0.7290640,  0.8076643,  0.8170951,
     &   0.8424572,  0.8684053,  0.8366923,  0.9418461,  0.8478122,
     &   0.8599583,  0.9175700,  0.8596328,  0.9279871,  0.8512805,
     &   1.0449820,  0.9670631,  0.9857884,  0.9676313,  1.0129299,
     &   0.9657040,  1.0019855,  1.0359297,  1.0414677,  0.9471506 /
      data y /
     &  -0.0310206,   0.1586742,   0.2576924,   0.3414014,   0.4943596,
     &   0.5782854,   0.6993418,   0.7470194,   0.9107649,   0.9962890,
     &   0.0501330,   0.0918555,   0.2592973,   0.3381592,   0.4171125,
     &   0.5615563,   0.6552235,   0.7524066,   0.9146523,   0.9632421,
     &   0.0292939,   0.0602303,   0.2668783,   0.3696044,   0.4801738,
     &   0.5940595,   0.6878797,   0.8185576,   0.9046507,   0.9805412,
     &   0.0396955,   0.0684484,   0.2389548,   0.3124129,   0.4902989,
     &   0.5199303,   0.6445227,   0.8203789,   0.8938079,   0.9711719,
     &  -0.0284618,   0.1560965,   0.2262471,   0.3175094,   0.3891417,
     &   0.5084949,   0.6324247,   0.7511007,   0.8489712,   0.9978728,
     &  -0.0271948,   0.1272430,   0.2709269,   0.3477728,   0.4259422,
     &   0.6084711,   0.6733781,   0.7235242,   0.9242411,   1.0308762,
     &   0.0255959,   0.0707835,   0.2008336,   0.3259843,   0.4890704,
     &   0.5096324,   0.6697880,   0.7759569,   0.9366096,   1.0064516,
     &   0.0285374,   0.1021403,   0.1936581,   0.3235775,   0.4714228,
     &   0.6091595,   0.6685053,   0.8022808,   0.8476790,   1.0512371,
     &   0.0380499,   0.0902048,   0.2083092,   0.3318491,   0.4335632,
     &   0.5910139,   0.6307383,   0.8144841,   0.9042310,   0.9696030,
     &  -0.0120900,   0.1334114,   0.2695844,   0.3795281,   0.4396054,
     &   0.5044425,   0.6941519,   0.7459923,   0.8682081,   0.9801409 /

      do i = 1, gn
        gx(i) = x(i)
        gy(i) = y(i)
      end do

      return
      end
      subroutine g02_size ( gn )

c*********************************************************************72
c
cc G02_SIZE returns the size for grid 2.
c
c  Discussion:
c
c    The "grid size" is simply the number of points in the grid.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    28 January 2012
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, integer GN, the grid size.
c
      implicit none

      integer gn

      gn = 33

      return
      end
      subroutine g02_title ( gt )

c*********************************************************************72
c
cc G02_TITLE returns the title for grid 2.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    28 January 2012
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, character * ( * ) GT, the grid title.
c
      implicit none

      character * ( * ) gt

      gt = 'Franke''s 33 node set'

      return
      end
      subroutine g02_xy ( gn, gx, gy )

c*********************************************************************72
c
cc G02_XY returns the grid points for grid 2.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    28 January 2012
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer GN, the grid size.
c
c    Output, double precision GX(GN), GY(GN), the grid coordinates.
c
      implicit none

      integer gn

      double precision gx(gn)
      double precision gy(gn)
      integer i
      double precision x(33)
      double precision y(33)

      save x
      save y

      data x /
     &  0.05,  0.00,  0.00,  0.00,  0.10,
     &  0.10,  0.15,  0.20,  0.25,  0.30,
     &  0.35,  0.50,  0.50,  0.55,  0.60,
     &  0.60,  0.60,  0.65,  0.70,  0.70,
     &  0.70,  0.75,  0.75,  0.75,  0.80,
     &  0.80,  0.85,  0.90,  0.90,  0.95,
     &  1.00,  1.00,  1.00 /
      data y /
     &  0.45,  0.50,  1.00,  0.00,  0.15,
     &  0.75,  0.30,  0.10,  0.20,  0.35,
     &  0.85,  0.00,  1.00,  0.95,  0.25,
     &  0.65,  0.85,  0.70,  0.20,  0.65,
     &  0.90,  0.10,  0.35,  0.85,  0.40,
     &  0.65,  0.25,  0.35,  0.80,  0.90,
     &  0.00,  0.50,  1.00 /

      do i = 1, gn
        gx(i) = x(i)
        gy(i) = y(i)
      end do

      return
      end
      subroutine g03_size ( gn )

c*********************************************************************72
c
cc G03_SIZE returns the size for grid 3.
c
c  Discussion:
c
c    The "grid size" is simply the number of points in the grid.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    28 January 2012
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, integer GN, the grid size.
c
      implicit none

      integer gn

      gn = 25

      return
      end
      subroutine g03_title ( gt )

c*********************************************************************72
c
cc G03_TITLE returns the title for grid 3.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    28 January 2012
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, character * ( * ) GT, the grid title.
c
      implicit none

      character * ( * ) gt

      gt = 'Lawson''s 25 node set'

      return
      end
      subroutine g03_xy ( gn, gx, gy )

c*********************************************************************72
c
cc G03_XY returns the grid points for grid 3.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    28 January 2012
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer GN, the grid size.
c
c    Output, double precision GX(GN), GY(GN), the grid coordinates.
c
      implicit none

      integer gn

      double precision gx(gn)
      double precision gy(gn)
      integer i
      double precision x(25)
      double precision y(25)

      save x
      save y

      data x /
     &  0.13750,   0.91250,   0.71250,   0.22500,  -0.05000,
     &  0.47500,   0.05000,   0.45000,   1.08750,   0.53750,
     & -0.03750,   0.18750,   0.71250,   0.85000,   0.70000,
     &  0.27500,   0.45000,   0.81250,   0.45000,   1.00000,
     &  0.50000,   0.18750,   0.58750,   1.05000,   0.10000 /
      data y /
     &  0.97500,   0.98750,   0.76250,   0.83750,   0.41250,
     &  0.63750,  -0.05000,   1.03750,   0.55000,   0.80000,
     &  0.75000,   0.57500,   0.55000,   0.43750,   0.31250,
     &  0.42500,   0.28750,   0.18750,  -0.03750,   0.26250,
     &  0.46250,   0.26250,   0.12500,  -0.06125,   0.11250 /

      do i = 1, gn
        gx(i) = x(i)
        gy(i) = y(i)
      end do

      return
      end
      subroutine g04_size ( gn )

c*********************************************************************72
c
cc G04_SIZE returns the size for grid 4.
c
c  Discussion:
c
c    The "grid size" is simply the number of points in the grid.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    28 January 2012
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, integer GN, the grid size.
c
      implicit none

      integer gn

      gn = 100

      return
      end
      subroutine g04_title ( gt )

c*********************************************************************72
c
cc G04_TITLE returns the title for grid 4.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    28 January 2012
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, character * ( * ) GT, the grid title.
c
      implicit none

      character * ( * ) gt

      gt = 'Random 100 node set'

      return
      end
      subroutine g04_xy ( gn, gx, gy )

c*********************************************************************72
c
cc G04_XY returns the grid points for grid 4.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    28 January 2012
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer GN, the grid size.
c
c    Output, double precision GX(GN), GY(GN), the grid coordinates.
c
      implicit none

      integer gn

      double precision gx(gn)
      double precision gy(gn)
      integer i
      double precision x(100)
      double precision y(100)

      save x
      save y

      data x /
     &  0.0096326,  0.0216348,  0.0298360,  0.0417447,  0.0470462,
     &  0.0562965,  0.0646857,  0.0740377,  0.0873907,  0.0934832,
     &  0.1032216,  0.1110176,  0.1181193,  0.1251704,  0.1327330,
     &  0.1439536,  0.1564861,  0.1651043,  0.1786039,  0.1886405,
     &  0.2016706,  0.2099886,  0.2147003,  0.2204141,  0.2343715,
     &  0.2409660,  0.2527740,  0.2570839,  0.2733365,  0.2853833,
     &  0.2901755,  0.2964854,  0.3019725,  0.3125695,  0.3307163,
     &  0.3378504,  0.3439061,  0.3529922,  0.3635507,  0.3766172,
     &  0.3822429,  0.3869838,  0.3973137,  0.4170708,  0.4255588,
     &  0.4299218,  0.4372839,  0.4705033,  0.4736655,  0.4879299,
     &  0.4940260,  0.5055324,  0.5162593,  0.5219219,  0.5348529,
     &  0.5483213,  0.5569571,  0.5638611,  0.5784908,  0.5863950,
     &  0.5929148,  0.5987839,  0.6117561,  0.6252296,  0.6331381,
     &  0.6399048,  0.6488972,  0.6558537,  0.6677405,  0.6814074,
     &  0.6887812,  0.6940896,  0.7061687,  0.7160957,  0.7317445,
     &  0.7370798,  0.7462030,  0.7566957,  0.7699998,  0.7879347,
     &  0.7944014,  0.8164468,  0.8192794,  0.8368405,  0.8500993,
     &  0.8588255,  0.8646496,  0.8792329,  0.8837536,  0.8900077,
     &  0.8969894,  0.9044917,  0.9083947,  0.9203972,  0.9347906,
     &  0.9434519,  0.9490328,  0.9569571,  0.9772067,  0.9983493 /
      data y /
     &  0.3083158,  0.2450434,  0.8613847,  0.0977864,  0.3648355,
     &  0.7156339,  0.5311312,  0.9755672,  0.1781117,  0.5452797,
     &  0.1603881,  0.7837139,  0.9982015,  0.6910589,  0.1049580,
     &  0.8184662,  0.7086405,  0.4456593,  0.1178342,  0.3189021,
     &  0.9668446,  0.7571834,  0.2016598,  0.3232444,  0.4368583,
     &  0.8907869,  0.0647260,  0.5692618,  0.2947027,  0.4332426,
     &  0.3347464,  0.7436284,  0.1066265,  0.8845357,  0.5158730,
     &  0.9425637,  0.4799701,  0.1783069,  0.1146760,  0.8225797,
     &  0.2270688,  0.4073598,  0.8875080,  0.7631616,  0.9972804,
     &  0.4959884,  0.3410421,  0.2498120,  0.6409007,  0.1058690,
     &  0.5411969,  0.0089792,  0.8784268,  0.5515874,  0.4038952,
     &  0.1654023,  0.2965158,  0.3660356,  0.0366554,  0.9502420,
     &  0.2638101,  0.9277386,  0.5377694,  0.7374676,  0.4674627,
     &  0.9186109,  0.0416884,  0.1291029,  0.6763676,  0.8444238,
     &  0.3273328,  0.1893879,  0.0645923,  0.0180147,  0.8904992,
     &  0.4160648,  0.4688995,  0.2174508,  0.5734231,  0.8853319,
     &  0.8018436,  0.6388941,  0.8931002,  0.1000558,  0.2789506,
     &  0.9082948,  0.3259159,  0.8318747,  0.0508513,  0.9708450,
     &  0.5120548,  0.2859716,  0.9581641,  0.6183429,  0.3779934,
     &  0.4010423,  0.9478657,  0.7425486,  0.8883287,  0.5496750 /

      do i = 1, gn
        gx(i) = x(i)
        gy(i) = y(i)
      end do

      return
      end
      subroutine g05_size ( gn )

c*********************************************************************72
c
cc G05_SIZE returns the size for grid 5.
c
c  Discussion:
c
c    The "grid size" is simply the number of points in the grid.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    28 January 2012
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, integer GN, the grid size.
c
      implicit none

      integer gn

      gn = 81

      return
      end
      subroutine g05_title ( gt )

c*********************************************************************72
c
cc G05_TITLE returns the title for grid 5.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    28 January 2012
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, character * ( * ) GT, the grid title.
c
      implicit none

      character * ( * ) gt

      gt = 'Gridded 81 node set'

      return
      end
      subroutine g05_xy ( gn, gx, gy )

c*********************************************************************72
c
cc G05_XY returns the grid points for grid 5.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    01 August 2012
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer GN, the grid size.
c
c    Output, double precision GX(GN), GY(GN), the grid coordinates.
c
      implicit none

      integer gn

      double precision gx(gn)
      double precision gy(gn)
      integer i
      double precision x(81)
      double precision y(81)

      save x
      save y

      data x /
     &  0.000,  0.000,  0.000,  0.000,  0.000,
     &  0.000,  0.000,  0.000,  0.000,  0.125,
     &  0.125,  0.125,  0.125,  0.125,  0.125,
     &  0.125,  0.125,  0.125,  0.250,  0.250,
     &  0.250,  0.250,  0.250,  0.250,  0.250,
     &  0.250,  0.250,  0.375,  0.375,  0.375,
     &  0.375,  0.375,  0.375,  0.375,  0.375,
     &  0.375,  0.500,  0.500,  0.500,  0.500,
     &  0.500,  0.500,  0.500,  0.500,  0.500,
     &  0.625,  0.625,  0.625,  0.625,  0.625,
     &  0.625,  0.625,  0.625,  0.625,  0.750,
     &  0.750,  0.750,  0.750,  0.750,  0.750,
     &  0.750,  0.750,  0.750,  0.875,  0.875,
     &  0.875,  0.875,  0.875,  0.875,  0.875,
     &  0.875,  0.875,  1.000,  1.000,  1.000,
     &  1.000,  1.000,  1.000,  1.000,  1.000,
     &  1.000 /
      data y /
     &  0.000,  0.125,  0.250,  0.375,  0.500,
     &  0.625,  0.750,  0.875,  1.000,  0.000,
     &  0.125,  0.250,  0.375,  0.500,  0.625,
     &  0.750,  0.875,  1.000,  0.000,  0.125,
     &  0.250,  0.375,  0.500,  0.625,  0.750,
     &  0.875,  1.000,  0.000,  0.125,  0.250,
     &  0.375,  0.500,  0.625,  0.750,  0.875,
     &  1.000,  0.000,  0.125,  0.250,  0.375,
     &  0.500,  0.625,  0.750,  0.875,  1.000,
     &  0.000,  0.125,  0.250,  0.375,  0.500,
     &  0.625,  0.750,  0.875,  1.000,  0.000,
     &  0.125,  0.250,  0.375,  0.500,  0.625,
     &  0.750,  0.875,  1.000,  0.000,  0.125,
     &  0.250,  0.375,  0.500,  0.625,  0.750,
     &  0.875,  1.000,  0.000,  0.125,  0.250,
     &  0.375,  0.500,  0.625,  0.750,  0.875,
     &  1.000 /

      do i = 1, gn
        gx(i) = x(i)
        gy(i) = y(i)
      end do

      return
      end
