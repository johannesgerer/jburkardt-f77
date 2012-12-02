      program main

c*********************************************************************72
c
cc TOMS365_PRB tests the CRF routine.
c
c  Modified:
c
c    29 November 2007
c
c  Author:
c
c    John Burkardt
c
      implicit none

      call timestamp ( )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TOMS365_PRB:'
      write ( *, '(a)' ) '  FORTRAN77 version'
      write ( *, '(a)' ) '  Test the TOMS365 library.'

      call test01 ( )
      call test02 ( )
      call test03 ( )
c
c  Terminate.
c
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TOMS365_PRB:'
      write ( *, '(a)' ) '  Normal end of execution.'

      write ( *, '(a)' ) ' '
      call timestamp ( )

      stop
      end
      subroutine test01 ( )

c*********************************************************************72
c
cc TEST01 tests CRF on the function F(Z) = Z + 1.
c
c  Modified:
c
c    08 January 2006
c
c  Author:
c
c    John Burkardt
c
      implicit none

      real de
      real dm
      real ds
      real he
      real hm
      real hs
      external f01
      complex f01
      integer n
      complex ze
      complex zs

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST01:'
      write ( *, '(a)' ) '  CRF uses the downhill method to find'
      write ( *, '(a)' ) '  a root of a complex analytic function.'
      write ( *, '(a)' ) '  Here, we use F(Z) = Z + 1.'

      zs = ( 2.0E+00, 0.5E+00 )
      hs = 0.25E+00
      hm = 0.0001E+00
      dm = 0.00001E+00

      write ( *, '(a)' ) ' '
      write ( *, '(a,2g14.6)' ) '  Initial estimate ZS =     ', zs
      write ( *, '(a,g14.6)' ) '  Initial stepsize HS =     ', hs
      write ( *, '(a,g14.6)' ) '  Minimum stepsize HM =     ', hm
      write ( *, '(a,g14.6)' ) '  Deviation tolerance DM =  ', dm

      call crf ( zs, hs, hm, dm, f01, ds, ze, he, de, n )

      write ( *, '(a)' ) ' '
      write ( *, '(a,2g14.6)' ) '  Final estimate ZE =       ', ze
      write ( *, '(a,g14.6)' ) '  Final stepsize HE =       ', he
      write ( *, '(a,g14.6)' ) '  Initial deviation DS =    ', ds
      write ( *, '(a,g14.6)' ) '  Final deviation DE =      ', de
      write ( *, '(a,i8)' ) '  Number of iterations, N = ', n

      return
      end
      subroutine test02 ( )

c*********************************************************************72
c
cc TEST02 tests CRF on the function F(Z) = Z**5 + 1.
c
c  Modified:
c
c    08 January 2006
c
c  Author:
c
c    John Burkardt
c
      implicit none

      real de
      real dm
      real ds
      real he
      real hm
      real hs
      external f02
      complex f02
      integer n
      complex ze
      complex zs

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST02:'
      write ( *, '(a)' ) '  CRF uses the downhill method to find'
      write ( *, '(a)' ) '  a root of a complex analytic function.'
      write ( *, '(a)' ) '  Here, we use F(Z) = Z**5 + 1.'

      zs = ( 2.0E+00, 0.5E+00 )
      hs = 0.25E+00
      hm = 0.0001E+00
      dm = 0.00001E+00

      write ( *, '(a)' ) ' '
      write ( *, '(a,2g14.6)' ) '  Initial estimate ZS =     ', zs
      write ( *, '(a,g14.6)' ) '  Initial stepsize HS =     ', hs
      write ( *, '(a,g14.6)' ) '  Minimum stepsize HM =     ', hm
      write ( *, '(a,g14.6)' ) '  Deviation tolerance DM =  ', dm

      call crf ( zs, hs, hm, dm, f02, ds, ze, he, de, n )

      write ( *, '(a)' ) ' '
      write ( *, '(a,2g14.6)' ) '  Final estimate ZE =       ', ze
      write ( *, '(a,g14.6)' ) '  Final stepsize HE =       ', he
      write ( *, '(a,g14.6)' ) '  Initial deviation DS =    ', ds
      write ( *, '(a,g14.6)' ) '  Final deviation DE =      ', de
      write ( *, '(a,i8)' ) '  Number of iterations, N = ', n

      return
      end
      subroutine test03 ( )

c*********************************************************************72
c
cc TEST03 tests CRF on W - sqrt ( Z*Z - 1 ) - acosh ( Z )
c
c  Modified:
c
c    29 November 2007
c
c  Author:
c
c    John Burkardt
c
      implicit none

      complex cacosh
      real de
      real dm
      real ds
      real he
      real hm
      real hs
      external f03
      complex f03
      complex fz
      integer n
      complex w
      complex ze
      complex zs

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST03:'
      write ( *, '(a)' ) '  CRF uses the downhill method to find'
      write ( *, '(a)' ) '  a root of a complex analytic function.'
      write ( *, '(a)' ) '  F(Z) = W - sqrt ( Z*Z - 1 ) - ACOSH(Z).'

      zs = ( 2.0E+00, 0.5E+00 )
      hs = 0.25E+00
      hm = 0.0001E+00
      dm = 0.000001E+00

      write ( *, '(a)' ) ' '
      write ( *, '(a,2g14.6)' ) '  Initial estimate ZS =     ', zs
      write ( *, '(a,g14.6)' ) '  Initial stepsize HS =     ', hs
      write ( *, '(a,g14.6)' ) '  Minimum stepsize HM =     ', hm
      write ( *, '(a,g14.6)' ) '  Deviation tolerance DM =  ', dm

      call crf ( zs, hs, hm, dm, f03, ds, ze, he, de, n )

      write ( *, '(a)' ) ' '
      write ( *, '(a,2g14.6)' ) '  Final estimate ZE =       ', ze
      write ( *, '(a,g14.6)' ) '  Final stepsize HE =       ', he
      write ( *, '(a,g14.6)' ) '  Initial deviation DS =    ', ds
      write ( *, '(a,g14.6)' ) '  Final deviation DE =      ', de
      write ( *, '(a,i8)' ) '  Number of iterations, N = ', n
c
c  The value of W here should correspond to the value in F03.
c  I prefer not to use a common block to guarantee this.
c
      w = ( 0.5E+00, 2.0E+00 )
      fz = csqrt ( ze * ze - 1.0E+00 ) + cacosh ( ze )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  We are actually solving W = F(Z), where'
      write ( *, '(a)' ) '  F(Z) = sqrt ( Z * Z - 1 ) + CACOSH ( Z ).'
      write ( *, '(a)' ) ' '
      write ( *, '(a,2g14.6)' ) '  W =    ', w
      write ( *, '(a,2g14.6)' ) '  Z =    ', ze
      write ( *, '(a,2g14.6)' ) '  F(Z) = ', fz

      return
      end
      function f01 ( z )

c*********************************************************************72
c
cc F01 evaluates the function F(Z) = Z + 1.
c
c  Modified:
c
c    06 January 2006
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, complex Z, the aergument.
c
c    Output, complex F01, the function value.
c
      implicit none

      complex f01
      complex z

      f01 = z + 1.0E+00

      return
      end
      function f02 ( z )

c*********************************************************************72
c
cc F02 evaluates the function F(Z) = Z**5 + 1.
c
c  Modified:
c
c    06 January 2006
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, complex Z, the aergument.
c
c    Output, complex F02, the function value.
c
      implicit none

      complex f02
      complex z

      f02 = z**5 + 1.0E+00

      return
      end
      function f03 ( z )

c*********************************************************************72
c
cc F03 evaluates W - sqrt ( z * z - 1 ) - arccosh ( Z ).
c
c  Discussion:
c
c    The complex number W is a parameter, to be set in this routine.
c
c    Fortran77 does not have an inverse hyperbolic cosine function.
c
c    The CACOSH function is used to supply this value.
c
c  Modified:
c
c    29 November 2007
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, complex Z, the aergument.
c
c    Output, complex F01, the function value.
c
      implicit none

      complex cacosh
      complex f03
      complex w
      complex z
c
c  Define W.
c
      w = ( 0.5E+00, 2.0E+00 )
c
c  Evaluate the function.
c
      f03 = w - csqrt ( z * z - 1.0E+00 ) - cacosh ( z )

      return
      end
      function cacos ( z )

c*********************************************************************72
c
cc CACOS evaluates the complex inverse cosine.
c
c  Discussion:
c
c    Fortran77 does not have an intrinsic inverse cosine function.
c
c    Here we use the relationship:
c
c      CACOS ( Z ) = pi/2 - CASIN ( Z ).
c
c  Modified:
c
c    29 November 2007
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, complex Z, the aergument.
c
c    Output, complex CACOS, the function value.
c
      implicit none

      complex cacos
      complex casin
      complex z

      cacos = 1.57079632679489661923E+00 - casin ( z )

      return
      end
      function cacosh ( z )

c*********************************************************************72
c
cc CACOSH evaluates the complex inverse hyperbolic cosine.
c
c  Discussion:
c
c    Fortran77 does not have an intrinsic inverse hyperbolic
c    cosine function.
c
c    Here we use the relationship:
c
c      CACOSH ( Z ) = i * CACOS ( Z ).
c
c  Modified:
c
c    29 November 2007
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, complex Z, the aergument.
c
c    Output, complex CACOSH, the function value.
c
      implicit none

      complex cacos
      complex cacosh
      complex z

      cacosh = ( 0.0E+00, 1.0E+00 ) * cacos ( z )

      return
      end
      function casin ( z )

c*********************************************************************72
c
cc CASIN evaluates the complex inverse sine.
c
c  Discussion:
c
c    Fortran77 does not have an intrinsic inverse sine function.
c
c    Here we use the relationship:
c
c      CASIN ( Z ) = - i * log ( i * z + sqrt ( 1 - z * z ) )
c
c  Modified:
c
c    29 November 2007
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, complex Z, the aergument.
c
c    Output, complex CACOS, the function value.
c
      implicit none

      complex casin
      complex i
      complex z

      i = ( 0.0E+00, 1.0E+00 )

      casin = - i * clog ( i * z + sqrt ( 1.0E+00 - z * z ) )

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
