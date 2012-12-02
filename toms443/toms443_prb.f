      program main

c*********************************************************************72
c
cc TOMS443_PRB tests TOMS443.
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
      write ( *, '(a)' ) 'TOMS443_PRB'
      write ( *, '(a)' ) '  FORTRAN77 version'
      write ( *, '(a)' ) '  Test the TOMS443 library.'

      call test01 ( )
      call test02 ( )
c
c  Terminate.
c
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TOMS435_PRB'
      write ( *, '(a)' ) '  Normal end of execution.'

      write ( *, '(a)' ) ' '
      call timestamp ( )

      stop
      end
      subroutine test01 ( )

c*********************************************************************72
c
cc TEST01 tests WEW_A
c
      implicit none

      real en
      integer n_data
      real w1
      real w2
      real wew_a
      real x

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST01'
      write ( *, '(a)' ) '  Test WEW_A to evaluate'
      write ( *, '(a)' ) '  Lambert''s W function.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '      X         Exact Value      Computed'
      write ( *, '(a)' ) ' '

      n_data = 0

10    continue

      call lambert_w_values ( n_data, x, w1 )

      if ( n_data <= 0 ) then
        go to 20
      end if

      if ( x == 0.0 ) then
        w2 = 0.0
      else
        w2 = wew_a ( x, en )
      end if

      write ( *, '(2x,f12.4,2x,f8.4,2x,g16.8,2x,g16.8)' )
     &  x, w1, w2

      go to 10

20    continue

      return
      end
      subroutine test02 ( )

c*********************************************************************72
c
cc TEST02 tests WEW_B
c
      implicit none

      real en
      integer n_data
      real w1
      real w2
      real wew_b
      real x

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST02'
      write ( *, '(a)' ) '  Test WEW_B to evaluate'
      write ( *, '(a)' ) '  Lambert''s W function.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '      X         Exact Value      Computed'
      write ( *, '(a)' ) ' '

      n_data = 0

10    continue

      call lambert_w_values ( n_data, x, w1 )

      if ( n_data <= 0 ) then
        go to 20
      end if

      if ( x == 0.0 ) then
        w2 = 0.0
      else
        w2 = wew_b ( x, en )
      end if

      write ( *, '(2x,f12.4,2x,f8.4,2x,g16.8,2x,g16.8)' )
     &  x, w1, w2

      go to 10

20    continue

      return
      end
      subroutine lambert_w_values ( n_data, x, fx )

c*********************************************************************72
c
cc LAMBERT_W_VALUES returns some values of the Lambert W function.
c
c  Discussion:
c
c    The function W(X) is defined implicitly by:
c
c      W(X) * e^W(X) = X
c
c    The function is also known as the "Omega" function.
c
c    In Mathematica, the function can be evaluated by:
c
c      W = ProductLog [ X ]
c
c  Modified:
c
c    23 February 2005
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    R M Corless, G H Gonnet, D E Hare, D J Jeffrey, D E Knuth,
c    On the Lambert W Function,
c    Advances in Computational Mathematics,
c    Volume 5, 1996, pages 329-359.
c
c    Brian Hayes,
c    "Why W?",
c    The American Scientist,
c    Volume 93, March-April 2005, pages 104-108.
c
c    Eric Weisstein,
c    "Lambert's W-Function",
c    CRC Concise Encyclopedia of Mathematics,
c    CRC Press, 1998.
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
c    Output, real X, the argument of the function.
c
c    Output, real FX, the value of the function.
c
      implicit none

      integer n_max
      parameter ( n_max = 22 )

      real fx
      real fx_vec(n_max)
      integer n_data
      real x
      real x_vec(n_max)

      save fx_vec
      save x_vec

      data fx_vec /
     &  0.0000000000000000D+00,
     &  0.3517337112491958D+00,
     &  0.5671432904097839D+00,
     &  0.7258613577662263D+00,
     &  0.8526055020137255D+00,
     &  0.9585863567287029D+00,
     &  0.1000000000000000D+01,
     &  0.1049908894964040D+01,
     &  0.1130289326974136D+01,
     &  0.1202167873197043D+01,
     &  0.1267237814307435D+01,
     &  0.1326724665242200D+01,
     &  0.1381545379445041D+01,
     &  0.1432404775898300D+01,
     &  0.1479856830173851D+01,
     &  0.1524345204984144D+01,
     &  0.1566230953782388D+01,
     &  0.1605811996320178D+01,
     &  0.1745528002740699D+01,
     &  0.3385630140290050D+01,
     &  0.5249602852401596D+01,
     &  0.1138335808614005D+02 /
      data x_vec /
     &  0.0000000000000000D+00,
     &  0.5000000000000000D+00,
     &  0.1000000000000000D+01,
     &  0.1500000000000000D+01,
     &  0.2000000000000000D+01,
     &  0.2500000000000000D+01,
     &  0.2718281828459045D+01,
     &  0.3000000000000000D+01,
     &  0.3500000000000000D+01,
     &  0.4000000000000000D+01,
     &  0.4500000000000000D+01,
     &  0.5000000000000000D+01,
     &  0.5500000000000000D+01,
     &  0.6000000000000000D+01,
     &  0.6500000000000000D+01,
     &  0.7000000000000000D+01,
     &  0.7500000000000000D+01,
     &  0.8000000000000000D+01,
     &  0.1000000000000000D+02,
     &  0.1000000000000000D+03,
     &  0.1000000000000000D+04,
     &  0.1000000000000000D+07 /

      if ( n_data < 0 ) then
        n_data = 0
      end if

      n_data = n_data + 1

      if ( n_max < n_data ) then
        n_data = 0
        x = 0.0
        fx = 0.0
      else
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
