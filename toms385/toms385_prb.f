      program main

c*********************************************************************72
c
cc TOMS385_PRB tests TOMS385.
c
c  Modified:
c
c    10 January 2006
c
c  Author:
c
c    John Burkardt
c
      implicit none

      call timestamp ( )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TOMS385_PRB'
      write ( *, '(a)' ) '  FORTRAN77 version'
      write ( *, '(a)' ) '  Test the TOMS385 library.'

      call test01 ( )
c
c  Terminate.
c
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TOMS385_PRB'
      write ( *, '(a)' ) '  Normal end of execution.'

      write ( *, '(a)' ) ' '
      call timestamp ( )

      stop
      end
      subroutine test01 ( )

c*********************************************************************72
c
cc TEST01 tests DEI.
c
c  Modified:
c
c    10 January 2006
c
c  Author:
c
c    John Burkardt
c
      implicit none

      double precision dei
      double precision f1
      double precision f2
      integer n_data
      double precision x

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST01'
      write ( *, '(a)' ) '  Test DEI to evaluate'
      write ( *, '(a)' ) '  the exponential integral function.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '          X         Exact Value      Computed'
      write ( *, '(a)' ) ' '

      n_data = 0

10    continue

      call ei_values ( n_data, x, f1 )

      if ( n_data <= 0 ) then
        go to 20
      end if

      f2 = dei ( x )

      write ( *, '(2x,f12.4,2x,g16.8,2x,g16.8)' ) x, f1, f2

      go to 10

20    continue

      return
      end
      subroutine ei_values ( n_data, x, fx )

c*********************************************************************72
c
cc EI_VALUES returns some values of the exponential integral function EI(X).
c
c  Discussion:
c
c    The exponential integral EI(X) has the formula:
c
c      EI(X) = - integral ( -X <= T <= Infinity ) exp ( -T ) / T dT
c
c    In Mathematica, the function can be evaluated by:
c
c      ExpIntegralEi[x]
c
c  Modified:
c
c    10 January 2006
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
c    Output, double precision X, the argument of the function.
c
c    Output, double precision FX, the value of the function.
c
      implicit none

      integer n_max
      parameter ( n_max = 16 )

      double precision fx
      double precision fx_vec(n_max)
      integer n_data
      double precision x
      double precision x_vec(n_max)

      save fx_vec
      save x_vec

      data fx_vec /
     &  0.4542199048631736D+00,
     &  0.7698812899373594D+00,
     &  0.1064907194624291D+01,
     &  0.1347396548212326D+01,
     &  0.1622811713696867D+01,
     &  0.1895117816355937D+01,
     &  0.2167378279563403D+01,
     &  0.2442092285192652D+01,
     &  0.2721398880232024D+01,
     &  0.3007207464150646D+01,
     &  0.3301285449129798D+01,
     &  0.3605319949019469D+01,
     &  0.3920963201354904D+01,
     &  0.4249867557487934D+01,
     &  0.4593713686953585D+01,
     &  0.4954234356001890D+01 /
      data x_vec /
     &  0.5D+00,
     &  0.6D+00,
     &  0.7D+00,
     &  0.8D+00,
     &  0.9D+00,
     &  1.0D+00,
     &  1.1D+00,
     &  1.2D+00,
     &  1.3D+00,
     &  1.4D+00,
     &  1.5D+00,
     &  1.6D+00,
     &  1.7D+00,
     &  1.8D+00,
     &  1.9D+00,
     &  2.0D+00 /

      if ( n_data < 0 ) then
        n_data = 0
      end if

      n_data = n_data + 1

      if ( n_max < n_data ) then
        n_data = 0
        x = 0.0D+00
        fx = 0.0D+00
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
