      program main

c*********************************************************************72
c
cc TOMS418_PRB tests FSPL2
c
c  Modified:
c
c    12 January 2006
c
c  Author:
c
c    John Burkardt
c
      implicit none

      call timestamp ( )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TOMS418_PRB:'
      write ( *, '(a)' ) '  FORTRAN77 version'
      write ( *, '(a)' ) '  Test the TOMS418 library.'

      call test01 ( )
      call test02 ( )
c
c  Terminate.
c
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TOMS418_PRB'
      write ( *, '(a)' ) '  Normal end of execution.'

      write ( *, '(a)' ) ' '
      call timestamp ( )

      stop
      end
      subroutine test01 ( )

c*********************************************************************72
c
cc TEST01 tests FSPL2 with integrands of the form F(X)*COS(W*X).
c
      implicit none

      real a
      real b
      real c
      real eps
      real error
      real exact
      external f1
      real f1_p
      real f1_pp
      external f2
      real f2_p
      real f2_pp
      external f3
      real f3_p
      real f3_pp
      real fba
      real fbb
      real fpa
      real fpb
      integer i
      integer j
      integer k
      integer lc
      integer ls
      integer max
      real pi
      real s
      real w

      save pi

      data pi / 3.141592653589793E+00 /

      a = 0.0E+00
      b = 2.0E+00 * pi

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST01'
      write ( *, '(a)' ) '  Use FSPL2 to estimate the integrals of'
      write ( *, '(a)' ) '  the form F(X) * COS ( W * X )'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Use integrand factors:'
      write ( *, '(a)' ) '  F(X) = 1, X, X*X.'
      write ( *, '(a)' ) ' '
      write ( *, '(a,g24.16)' ) '  A = ', a
      write ( *, '(a,g24.16)' ) '  B = ', b
      write ( *, '(a)' ) ' '
      write ( *, '(a,a)' )
     &  '             W                      Approximate',
     &  '             Exact'
      write ( *, '(a)' ) ' '

      max = 12
      lc = 1
      ls = -1

      do k = 1, 3

        if ( k == 1 ) then
          w = 1.0E+00
        else if ( k == 2 ) then
          w = 2.0E+00
        else if ( k == 3 ) then
          w = 10.0E+00
        end if

        do i = 1, 3

          eps = 0.00001E+00

          if ( i == 1 )then

            fpa = f1_p ( a )
            fpb = f1_p ( b )
            fba = f1_pp ( a )
            fbb = f1_pp ( b )

            call fspl2 ( f1, a, b, fpa, fpb, fba, fbb, w, eps, max,
     &        lc, ls, c, s )

            exact = ( sin ( w * b ) - sin ( w * a ) ) / w

            if ( eps < 0.0E+00 ) then
              write ( *, '(a)' ) '  Next result did not converge:'
            end if

            write ( *, '(a,3g24.16)' ) '    1:  ', w, c, exact

          else if ( i == 2 ) then

            fpa = f2_p ( a )
            fpb = f2_p ( b )
            fba = f2_pp ( a )
            fbb = f2_pp ( b )

            call fspl2 ( f2, a, b, fpa, fpb, fba, fbb, w, eps, max,
     &        lc, ls, c, s )

            exact = ( ( cos ( w * b ) + w * b * sin ( w * b ) )
     &               - ( cos ( w * a ) + w * a * sin ( w * a ) ) )
     &               / w**2

            if ( eps < 0.0E+00 ) then
              write ( *, '(a)' ) '  Next result did not converge:'
            end if

            write ( *, '(a,3g24.16)' ) '    X:  ', w, c, exact

          else if ( i == 3 ) then

            fpa = f3_p ( a )
            fpb = f3_p ( b )
            fba = f3_pp ( a )
            fbb = f3_pp ( b )

            call fspl2 ( f3, a, b, fpa, fpb, fba, fbb, w, eps, max,
     &        lc, ls, c, s )

            exact = ( ( 2.0E+00 * w * b * cos ( w * b )
     &            + ( w * w * b**2 - 2.0E+00 ) * sin ( w * b ) )
     &              - ( 2.0E+00 * w * a * cos ( w * a )
     &            + ( w * w * a**2 - 2.0E+00 ) * sin ( w * a ) ) )
     &            / w**3

            if ( eps < 0.0E+00 ) then
              write ( *, '(a)' ) '  Next result did not converge:'
            end if

            write ( *, '(a,3g24.16)' ) '  X*X:  ', w, c, exact

          end if

        end do

        write ( *, '(a)' ) ' '

      end do

      return
      end
      subroutine test02 ( )

c*********************************************************************72
c
cc TEST02 tests FSPL2 with integrands of the form F(X)*SIN(W*X).
c
      implicit none

      real a
      real b
      real c
      real eps
      real error
      real exact
      external f1
      real f1_p
      real f1_pp
      external f2
      real f2_p
      real f2_pp
      external f3
      real f3_p
      real f3_pp
      real fba
      real fbb
      real fpa
      real fpb
      integer i
      integer j
      integer k
      integer lc
      integer ls
      integer max
      real pi
      real s
      real w

      save pi

      data pi / 3.141592653589793E+00 /

      a = 0.0E+00
      b = 2.0E+00 * pi

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST02'
      write ( *, '(a)' ) '  Use FSPL2 to estimate the integrals of'
      write ( *, '(a)' ) '  the form F(X) * SIN ( W * X )'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Use integrand factors:'
      write ( *, '(a)' ) '  F(X) = 1, X, X*X.'
      write ( *, '(a)' ) ' '
      write ( *, '(a,g24.16)' ) '  A = ', a
      write ( *, '(a,g24.16)' ) '  B = ', b
      write ( *, '(a)' ) ' '
      write ( *, '(a,a)' )
     &  '             W                      Approximate',
     &  '             Exact'
      write ( *, '(a)' ) ' '

      max = 12
      lc = -1
      ls = 1

      do k = 1, 3

        if ( k == 1 ) then
          w = 1.0E+00
        else if ( k == 2 ) then
          w = 2.0E+00
        else if ( k == 3 ) then
          w = 10.0E+00
        end if

        do i = 1, 3

          eps = 0.00001E+00

          if ( i == 1 )then

            fpa = f1_p ( a )
            fpb = f1_p ( b )
            fba = f1_pp ( a )
            fbb = f1_pp ( b )

            call fspl2 ( f1, a, b, fpa, fpb, fba, fbb, w, eps, max,
     &        lc, ls, c, s )

            exact = ( - cos ( w * b ) + cos ( w * a ) ) / w

            if ( eps < 0.0E+00 ) then
              write ( *, '(a)' ) '  Next result did not converge:'
            end if

            write ( *, '(a,3g24.16)' ) '    1:  ', w, s, exact

          else if ( i == 2 ) then

            fpa = f2_p ( a )
            fpb = f2_p ( b )
            fba = f2_pp ( a )
            fbb = f2_pp ( b )

            call fspl2 ( f2, a, b, fpa, fpb, fba, fbb, w, eps, max,
     &        lc, ls, c, s )

            exact = ( ( sin ( w * b ) - w * b * cos ( w * b ) )
     &               - ( sin ( w * a ) - w * a * cos ( w * a ) ) )
     &               / w**2

            if ( eps < 0.0E+00 ) then
              write ( *, '(a)' ) '  Next result did not converge:'
            end if

            write ( *, '(a,3g24.16)' ) '    X:  ', w, s, exact

          else if ( i == 3 ) then

            fpa = f3_p ( a )
            fpb = f3_p ( b )
            fba = f3_pp ( a )
            fbb = f3_pp ( b )

            call fspl2 ( f3, a, b, fpa, fpb, fba, fbb, w, eps, max,
     &        lc, ls, c, s )

            exact = ( ( 2.0E+00 * w * b * sin ( w * b )
     &            + ( 2.0E+00 - w**2 * b**2 ) * cos ( w * b ) )
     &            - ( 2.0E+00 * w * a * sin ( w * a )
     &            + ( 2.0E+00 - w**2 * a**2 ) * cos ( w * a ) ) )
     &            / w**3

            if ( eps < 0.0E+00 ) then
              write ( *, '(a)' ) '  Next result did not converge:'
            end if

            write ( *, '(a,3g24.16)' ) '  X*X:  ', w, s, exact

          end if

        end do

        write ( *, '(a)' ) ' '

      end do

      return
      end
      function f1 ( x )

c*********************************************************************72
c
cc F1 evaluates the integrand factor F(X) = 1.
c
      implicit none

      real f1
      real x

      f1 = 1.0E+00

      return
      end
      function f1_p ( x )

c*********************************************************************72
c
cc F1_P evaluates the first derivative of the integrand factor F(X) = 1.
c
      implicit none

      real f1_p
      real x

      f1_p = 0.0E+00

      return
      end
      function f1_pp ( x )

c*********************************************************************72
c
cc F1_PP evaluates the second derivative of the integrand factor F(X) = 1.
c
      implicit none

      real f1_pp
      real x

      f1_pp = 0.0E+00

      return
      end
      function f2 ( x )

c*********************************************************************72
c
cc F2 evaluates the integrand factor F(X) = X.
c
      implicit none

      real f2
      real x

      f2 = x

      return
      end
      function f2_p ( x )

c*********************************************************************72
c
cc F2_P evaluates the first derivative of the integrand factor F(X) = X.
c
      implicit none

      real f2_p
      real x

      f2_p = 1.0E+00

      return
      end
      function f2_pp ( x )

c*********************************************************************72
c
cc F2_PP evaluates the second derivative of the integrand factor F(X) = X.
c
      implicit none

      real f2_pp
      real x

      f2_pp = 0.0E+00

      return
      end
      function f3 ( x )

c*********************************************************************72
c
cc F3 evaluates the integrand factor F(X) = X*X.
c
      implicit none

      real f3
      real x

      f3 = x*x

      return
      end
      function f3_p ( x )

c*********************************************************************72
c
cc F3_P evaluates the first derivative of the integrand factor F(X) = X*X.
c
      implicit none

      real f3_p
      real x

      f3_p = x

      return
      end
      function f3_pp ( x )

c*********************************************************************72
c
cc F3_PP evaluates the second derivative of the integrand factor F(X) = X*X.
c
      implicit none

      real f3_pp
      real x

      f3_pp = 1.0E+00

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
c  Modified:
c
c    16 September 2005
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

      character ( len = 8 ) date
      character ( len = 10 ) time

      call date_and_time ( date, time )

      write ( *, '(a8,2x,a10)' ) date, time

      return
      end
