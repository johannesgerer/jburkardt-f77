      program main

c*********************************************************************72
c
cc MAIN is the main program for ASA047_PRB.
c
c  Discussion:
c
c    ASA047_PRB calls the ASA047 routines.
c
c  Modified:
c
c    02 February 2008
c
c  Author:
c
c    John Burkardt
c
      implicit none

      call timestamp ( )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'ASA047_PRB:'
      write ( *, '(a)' ) '  FORTRAN77 version'
      write ( *, '(a)' ) '  Test the ASA047 library.'

      call test01 ( )
      call test02 ( )
      call test03 ( )
      call test04 ( )
c
c  Terminate.
c
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'ASA047_PRB:'
      write ( *, '(a)' ) '  Normal end of execution.'

      write ( *, '(a)' ) ' '
      call timestamp ( )

      stop
      end
      subroutine test01 ( )

c*********************************************************************72
c
cc TEST01 demonstrates the use of NELMIN on ROSENBROCK.
c
c  Modified:
c
c    01 February 2008
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n
      parameter ( n = 2 )

      integer i
      integer icount
      integer ifault
      integer kcount
      integer konvge
      integer numres
      double precision reqmin
      double precision rosenbrock
      external rosenbrock
      double precision start(n)
      double precision step(n)
      double precision xmin(n)
      double precision ynewlo

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST01'
      write ( *, '(a)' ) '  Apply NELMIN to ROSENBROCK function.'

      start(1) = -1.2D+00
      start(2) =  1.0D+00

      reqmin = 1.0D-08

      step(1) = 1.0D+00
      step(2) = 1.0D+00

      konvge = 10
      kcount = 500

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Starting point X:'
      write ( *, '(a)' ) ' '
      do i = 1, n
        write ( *, '(2x,g14.6)' ) start(i)
      end do

      ynewlo = rosenbrock ( start )

      write ( *, '(a)' ) ' '
      write ( *, '(a,g14.6)' ) '  F(X) = ', ynewlo

      call nelmin ( rosenbrock, n, start, xmin, ynewlo, reqmin, step,
     &  konvge, kcount, icount, numres, ifault )

      write ( *, '(a)' ) ' '
      write ( *, '(a,i8)' ) '  Return code IFAULT = ', ifault
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Estimate of minimizing value X*:'
      write ( *, '(a)' ) ' '
      do i = 1, n
        write ( *, '(2x,g14.6)' ) xmin(i)
      end do

      write ( *, '(a)' ) ' '
      write ( *, '(a,g14.6)' ) '  F(X*) = ', ynewlo

      write ( *, '(a)' ) ' '
      write ( *, '(a,i8)' ) '  Number of iterations = ', icount
      write ( *, '(a,i8)' ) '  Number of restarts =   ', numres

      return
      end
      function rosenbrock ( x )

c*********************************************************************72
c
cc ROSENBROCK evaluates the Rosenbrock parabolic value function.
c
c  Modified:
c
c    01 February 2008
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    R ONeill,
c    Algorithm AS 47:
c    Function Minimization Using a Simplex Procedure,
c    Applied Statistics,
c    Volume 20, Number 3, 1971, pages 338-345.
c
c  Parameters:
c
c    Input, double precision X(2), the argument.
c
c    Output, double precision ROSENBROCK, the value of the function.
c
      implicit none

      double precision fx
      double precision fx1
      double precision fx2
      double precision rosenbrock
      double precision x(3)

      fx1 = x(2) - x(1) * x(1)
      fx2 = 1.0D+00 - x(1)

      fx = 100.0D+00 * fx1 * fx1
     &   +             fx2 * fx2

      rosenbrock = fx

      return
      end
      subroutine test02 ( )

c*********************************************************************72
c
cc TEST02 demonstrates the use of NELMIN on POWELL.
c
c  Modified:
c
c    02 February 2008
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n
      parameter ( n = 4 )

      integer i
      integer icount
      integer ifault
      integer kcount
      integer konvge
      integer numres
      double precision powell
      external powell
      double precision reqmin
      double precision start(n)
      double precision step(n)
      double precision xmin(n)
      double precision ynewlo

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST02'
      write ( *, '(a)' ) '  Apply NELMIN to POWELL quartic function.'

      start(1) =   3.0D+00
      start(2) = - 1.0D+00
      start(3) =   0.0D+00
      start(4) =   1.0D+00

      reqmin = 1.0D-08

      step(1) = 1.0D+00
      step(2) = 1.0D+00
      step(3) = 1.0D+00
      step(4) = 1.0D+00

      konvge = 10
      kcount = 500

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Starting point X:'
      write ( *, '(a)' ) ' '
      do i = 1, n
        write ( *, '(2x,g14.6)' ) start(i)
      end do

      ynewlo = powell ( start )

      write ( *, '(a)' ) ' '
      write ( *, '(a,g14.6)' ) '  F(X) = ', ynewlo

      call nelmin ( powell, n, start, xmin, ynewlo, reqmin, step,
     &  konvge, kcount, icount, numres, ifault )

      write ( *, '(a)' ) ' '
      write ( *, '(a,i8)' ) '  Return code IFAULT = ', ifault
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Estimate of minimizing value X*:'
      write ( *, '(a)' ) ' '
      do i = 1, n
        write ( *, '(2x,g14.6)' ) xmin(i)
      end do

      write ( *, '(a)' ) ' '
      write ( *, '(a,g14.6)' ) '  F(X*) = ', ynewlo

      write ( *, '(a)' ) ' '
      write ( *, '(a,i8)' ) '  Number of iterations = ', icount
      write ( *, '(a,i8)' ) '  Number of restarts =   ', numres

      return
      end
      function powell ( x )

c*********************************************************************72
c
cc POWELL evaluates the Powell quartic function.
c
c  Modified:
c
c    01 February 2008
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    R ONeill,
c    Algorithm AS 47:
c    Function Minimization Using a Simplex Procedure,
c    Applied Statistics,
c    Volume 20, Number 3, 1971, pages 338-345.
c
c  Parameters:
c
c    Input, double precision X(4), the argument.
c
c    Output, double precision POWELL, the value of the function.
c
      implicit none

      double precision fx
      double precision fx1
      double precision fx2
      double precision fx3
      double precision fx4
      double precision powell
      double precision x(4)

      fx1 = x(1) + 10.0D+00 * x(2)
      fx2 = x(3) - x(4)
      fx3 = x(2) - 2.0D+00 * x(3)
      fx4 = x(1) - x(4)

      fx =            fx1 * fx1
     &   +  5.0D+00 * fx2 * fx2
     &   +            fx3 * fx3 * fx3 * fx3
     &   + 10.0D+00 * fx4 * fx4 * fx4 * fx4

      powell = fx

      return
      end
      subroutine test03 ( )

c*********************************************************************72
c
cc TEST03 demonstrates the use of NELMIN on HELICAL.
c
c  Modified:
c
c    02 February 2008
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n
      parameter ( n = 3 )

      double precision helical
      external helical
      integer i
      integer icount
      integer ifault
      integer kcount
      integer konvge
      integer numres
      double precision reqmin
      double precision start(n)
      double precision step(n)
      double precision xmin(n)
      double precision ynewlo

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST03'
      write ( *, '(a)' ) '  Apply NELMIN to the HELICAL function.'

      start(1) = - 1.0D+00
      start(2) =   0.0D+00
      start(3) =   0.0D+00

      reqmin = 1.0D-08

      step(1) = 1.0D+00
      step(2) = 1.0D+00
      step(3) = 1.0D+00

      konvge = 10
      kcount = 500

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Starting point X:'
      write ( *, '(a)' ) ' '
      do i = 1, n
        write ( *, '(2x,g14.6)' ) start(i)
      end do

      ynewlo = helical ( start )

      write ( *, '(a)' ) ' '
      write ( *, '(a,g14.6)' ) '  F(X) = ', ynewlo

      call nelmin ( helical, n, start, xmin, ynewlo, reqmin, step,
     &  konvge, kcount, icount, numres, ifault )

      write ( *, '(a)' ) ' '
      write ( *, '(a,i8)' ) '  Return code IFAULT = ', ifault
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Estimate of minimizing value X*:'
      write ( *, '(a)' ) ' '
      do i = 1, n
        write ( *, '(2x,g14.6)' ) xmin(i)
      end do

      write ( *, '(a)' ) ' '
      write ( *, '(a,g14.6)' ) '  F(X*) = ', ynewlo

      write ( *, '(a)' ) ' '
      write ( *, '(a,i8)' ) '  Number of iterations = ', icount
      write ( *, '(a,i8)' ) '  Number of restarts =   ', numres

      return
      end
      function helical ( x )

c*********************************************************************72
c
cc HELICAL evaluates the Fletcher-Powell helical valley function.
c
c  Modified:
c
c    01 February 2008
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    R ONeill,
c    Algorithm AS 47:
c    Function Minimization Using a Simplex Procedure,
c    Applied Statistics,
c    Volume 20, Number 3, 1971, pages 338-345.
c
c  Parameters:
c
c    Input, double precision X(3), the argument.
c
c    Output, double precision HELICAL, the value of the function.
c
      implicit none

      double precision fx
      double precision fx1
      double precision fx2
      double precision fx3
      double precision helical
      double precision pi
      parameter ( pi = 3.141592653589793D+00 )
      double precision theta
      double precision x(3)

      if ( 0.0D+00 .lt. x(1) ) then
        theta = atan2 ( x(2), x(1) ) / 2.0D+00 / pi
      else if ( x(1) .lt. 0.0D+00 ) then
        theta = 0.5D+00 + atan2 ( x(2), x(1) ) / 2.0D+00 / pi
      else if ( x(1) .eq. 0.0D+00 ) then
        theta = 0.25D+00
      end if

      fx1 = x(3) - 10.0D+00 * theta
      fx2 = sqrt ( x(1) * x(1) + x(2) * x(2) )
      fx3 = x(3)

      fx = 100.0D+00 * fx1 * fx1
     &   +             fx2 * fx2
     &   +             fx3 * fx3

      helical = fx

      return
      end
      subroutine test04 ( )

c*********************************************************************72
c
cc TEST04 demonstrates the use of NELMIN on QUARTIC.
c
c  Modified:
c
c    02 February 2008
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n
      parameter ( n = 10 )

      integer i
      integer icount
      integer ifault
      integer kcount
      integer konvge
      integer numres
      double precision quartic
      external quartic
      double precision reqmin
      double precision start(n)
      double precision step(n)
      double precision xmin(n)
      double precision ynewlo

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST04'
      write ( *, '(a)' ) '  Apply NELMIN to the QUARTIC function.'

      do i = 1, n
        start(i) = 1.0D+00
      end do

      reqmin = 1.0D-08

      do i = 1, n
        step(i) = 1.0D+00
      end do

      konvge = 10
      kcount = 500

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Starting point X:'
      write ( *, '(a)' ) ' '
      do i = 1, n
        write ( *, '(2x,g14.6)' ) start(i)
      end do

      ynewlo = quartic ( start )

      write ( *, '(a)' ) ' '
      write ( *, '(a,g14.6)' ) '  F(X) = ', ynewlo

      call nelmin ( quartic, n, start, xmin, ynewlo, reqmin, step,
     &  konvge, kcount, icount, numres, ifault )

      write ( *, '(a)' ) ' '
      write ( *, '(a,i8)' ) '  Return code IFAULT = ', ifault
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Estimate of minimizing value X*:'
      write ( *, '(a)' ) ' '
      do i = 1, n
        write ( *, '(2x,g14.6)' ) xmin(i)
      end do

      write ( *, '(a)' ) ' '
      write ( *, '(a,g14.6)' ) '  F(X*) = ', ynewlo

      write ( *, '(a)' ) ' '
      write ( *, '(a,i8)' ) '  Number of iterations = ', icount
      write ( *, '(a,i8)' ) '  Number of restarts =   ', numres

      return
      end
      function quartic ( x )

c*********************************************************************72
c
cc QUARTIC evaluates a function defined by a sum of fourth powers.
c
c  Modified:
c
c    01 February 2008
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    R ONeill,
c    Algorithm AS 47:
c    Function Minimization Using a Simplex Procedure,
c    Applied Statistics,
c    Volume 20, Number 3, 1971, pages 338-345.
c
c  Parameters:
c
c    Input, double precision X(10), the argument.
c
c    Output, double precision QUARTIC, the value of the function.
c
      implicit none

      double precision fx
      integer i
      double precision quartic
      double precision x(10)

      fx = 0.0D+00

      do i = 1, 10
        fx = fx + x(i) * x(i) * x(i) * x(i)
      end do

      quartic = fx

      return
      end
