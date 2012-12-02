      program main

c*********************************************************************72
c
cc TOMS344_PRB tests TTEST
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
      write ( *, '(a)' ) 'TOMS344_PRB:'
      write ( *, '(a)' ) '  FORTRAN77 version'
      write ( *, '(a)' ) '  Test the TOMS344 library.'

      call test01 ( )
c
c  Terminate.
c
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TOMS344_PRB'
      write ( *, '(a)' ) '  Normal end of execution.'

      write ( *, '(a)' ) ' '
      call timestamp ( )

      stop
      end
      subroutine test01

c*********************************************************************72
c
cc TEST01 compares TTEST versus values from STUDENT_CDF_VALUES.
c
c  Modified:
c
c    06 January 2006
c
      implicit none

      real c
      integer df
      real fx
      real fx2
      integer kerr
      integer n_data
      real x

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST01:'
      write ( *, '(a)' ) '  STUDENT_CDF_VALUES returns values of '
      write ( *, '(a)' ) '  the Student T CDF.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' )
     &  '      DF        X         CDF                       CDF'
      write ( *, '(a)' )
     &  '                          exact                     computed'
      write ( *, '(a)' ) ' '

      n_data = 0

 10   continue

      call student_cdf_values ( n_data, c, x, fx )

      if ( n_data == 0 ) then
        go to 20
      end if

      df = int ( c )

      if ( c .ne. real ( df ) ) then
        go to 10
      end if

      call ttest ( x, df, fx2, kerr )
c
c  Transform FX.  TTEST computes the two-tailed distribution, and uses the
c  complementary form of the integral.
c
      fx = 2.0E+00 * ( 1.0E+00 - fx )

      write ( *, '(2x,i6,2x,f10.4,2x,g24.16,2x,g24.16)' ) df, x, fx, fx2

      go to 10

20    continue

      return
      end
      subroutine student_cdf_values ( n_data, c, x, fx )

c*********************************************************************72
c
cc STUDENT_CDF_VALUES returns some values of the Student CDF.
c
c  Discussion:
c
c    In Mathematica, the function can be evaluated by:
c
c      Needs["Statistics`ContinuousDistributions`"]
c      dist = StudentTDistribution [ c ]
c      CDF [ dist, x ]
c
c  Modified:
c
c    06 January 2005
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Milton Abramowitz and Irene Stegun,
c    Handbook of Mathematical Functions,
c    US Department of Commerce, 1964.
c
c    Stephen Wolfram,
c    The Mathematica Book,
c    Fourth Edition,
c    Wolfram Media / Cambridge University Press, 1999.
c
c  Parameters:
c
c    Input/output, integer N_DATA.  The user sets N_DATA to 0 before the
c    first call.  On each call, the routine increments N_DATA by 1, and
c    returns the corresponding data; when there is no more data, the
c    output value of N_DATA will be 0 again.
c
c    Output, real C, is usually called the number of
c    degrees of freedom of the distribution.  C is typically an
c    integer, but that is not essential.  It is required that
c    C be strictly positive.
c
c    Output, real X, the argument of the function.
c
c    Output, real FX, the value of the function.
c
      implicit none

      integer n_max
      parameter ( n_max = 13 )

      real c
      real c_vec(n_max)
      real fx
      real fx_vec(n_max)
      integer n_data
      real x
      real x_vec(n_max)

      save c_vec
      save fx_vec
      save x_vec

      data c_vec /
     &  1.0E+00,
     &  2.0E+00,
     &  3.0E+00,
     &  4.0E+00,
     &  5.0E+00,
     &  2.0E+00,
     &  5.0E+00,
     &  2.0E+00,
     &  5.0E+00,
     &  2.0E+00,
     &  3.0E+00,
     &  4.0E+00,
     &  5.0E+00 /
      data fx_vec /
     &  0.6000231200328521E+00,
     &  0.6001080279134390E+00,
     &  0.6001150934648930E+00,
     &  0.6000995134721354E+00,
     &  0.5999341989834830E+00,
     &  0.7498859393137811E+00,
     &  0.7500879487671045E+00,
     &  0.9500004222186464E+00,
     &  0.9499969138365968E+00,
     &  0.9900012348724744E+00,
     &  0.9900017619355059E+00,
     &  0.9900004567580596E+00,
     &  0.9900007637471291E+00 /
      data x_vec /
     &  0.325E+00,
     &  0.289E+00,
     &  0.277E+00,
     &  0.271E+00,
     &  0.267E+00,
     &  0.816E+00,
     &  0.727E+00,
     &  2.920E+00,
     &  2.015E+00,
     &  6.965E+00,
     &  4.541E+00,
     &  3.747E+00,
     &  3.365E+00 /

      if ( n_data < 0 ) then
        n_data = 0
      end if

      n_data = n_data + 1

      if ( n_max < n_data ) then
        n_data = 0
        c = 0.0E+00
        x = 0.0E+00
        fx = 0.0E+00
      else
        c = c_vec(n_data)
        x = x_vec(n_data)
        fx = fx_vec(n_data)
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
