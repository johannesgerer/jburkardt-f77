      program main

c*********************************************************************72
c
cc TOMS332_PRB tests the JACOBI routine.
c
c  Modified:
c
c    04 January 2006
c
c  Author:
c
c    John Burkardt
c
      implicit none

      call timestamp ( )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TOMS332_PRB:'
      write ( *, '(a)' ) '  FORTRAN77 version'
      write ( *, '(a)' ) '  Test the TOMS332 library.'

      call test01 ( )
c
c  Terminate.
c
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TOMS332_PRB:'
      write ( *, '(a)' ) '  Normal end of execution.'

      write ( *, '(a)' ) ' '
      call timestamp ( )

      stop
      end
      subroutine test01 ( )

c*********************************************************************72
c
cc TEST01 tests JACOBI against values stored in JACOBI_POLY_VALUES.
c
c  Modified:
c
c    04 January 2006
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer a
      double precision alfa
      integer b
      double precision beta
      double precision e
      double precision ed
      double precision fd
      integer flagf
      integer flagfd
      double precision fx
      double precision fx2
      integer n
      integer n_data
      double precision x

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST01:'
      write ( *, '(a)' ) '  JACOBI_POLY_VALUES returns exact values of '
      write ( *, '(a)' ) '  the Jacobi polynomial.'
      write ( *, '(a)' ) '  JACOBI computes values of '
      write ( *, '(a)' ) '  the Jacobi polynomial.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' )
     &  '   N       A       B      X          Exact           Computed'
      write ( *, '(a)' ) ' '

      n_data = 0

10    continue

        call jacobi_poly_values ( n_data, n, a, b, x, fx )

        if ( n_data == 0 ) then
          go to 20
        end if

        alfa = dble ( a )
        beta = dble ( b )

        call jacobi ( n, alfa, beta, x, fx2, fd, e, ed, flagf, flagfd )

        write ( *, '(2x,i2,2x,i6,2x,i6,2x,f10.6,2x,g16.8,2x,g16.8)' )
     &    n, a, b, x, fx, fx2

      go to 10

20    continue

      return
      end
      subroutine jacobi_poly_values ( n_data, n, a, b, x, fx )

c*********************************************************************72
c
cc JACOBI_POLY_VALUES returns some values of the Jacobi polynomial.
c
c  Discussion:
c
c    In Mathematica, the function can be evaluated by:
c
c      JacobiP[ n, a, b, x ]
c
c  Modified:
c
c    04 January 2006
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
c    Output, integer N, the degree of the polynomial.
c
c    Output, integer A, B, parameters of the function.
c
c    Output, double precision X, the argument of the function.
c
c    Output, double precision FX, the value of the function.
c
      implicit none

      integer n_max
      parameter ( n_max = 26 )

      integer a
      integer a_vec(n_max)
      integer b
      integer b_vec(n_max)
      double precision fx
      double precision fx_vec(n_max)
      integer n
      integer n_data
      integer n_vec(n_max)
      double precision x
      double precision x_vec(n_max)

      save a_vec
      save b_vec
      save fx_vec
      save n_vec
      save x_vec

      data a_vec /
     &    0, 0, 0, 0,
     &    0, 0, 1, 2,
     &    3, 4, 5, 0,
     &    0, 0, 0, 0,
     &    0, 0, 0, 0,
     &    0, 0, 0, 0,
     &    0, 0 /
      data b_vec /
     &   1, 1, 1, 1,
     &   1, 1, 1, 1,
     &   1, 1, 1, 2,
     &   3, 4, 5, 1,
     &   1, 1, 1, 1,
     &   1, 1, 1, 1,
     &   1, 1 /
      data fx_vec /
     &    0.1000000000000000D+01,
     &    0.2500000000000000D+00,
     &   -0.3750000000000000D+00,
     &   -0.4843750000000000D+00,
     &   -0.1328125000000000D+00,
     &    0.2753906250000000D+00,
     &   -0.1640625000000000D+00,
     &   -0.1174804687500000D+01,
     &   -0.2361328125000000D+01,
     &   -0.2616210937500000D+01,
     &    0.1171875000000000D+00,
     &    0.4218750000000000D+00,
     &    0.5048828125000000D+00,
     &    0.5097656250000000D+00,
     &    0.4306640625000000D+00,
     &   -0.6000000000000000D+01,
     &    0.3862000000000000D-01,
     &    0.8118400000000000D+00,
     &    0.3666000000000000D-01,
     &   -0.4851200000000000D+00,
     &   -0.3125000000000000D+00,
     &    0.1891200000000000D+00,
     &    0.4023400000000000D+00,
     &    0.1216000000000000D-01,
     &   -0.4396200000000000D+00,
     &    0.1000000000000000D+01 /
      data n_vec /
     &    0, 1, 2, 3,
     &    4, 5, 5, 5,
     &    5, 5, 5, 5,
     &    5, 5, 5, 5,
     &    5, 5, 5, 5,
     &    5, 5, 5, 5,
     &    5, 5 /
      data x_vec /
     &    0.5D+00,
     &    0.5D+00,
     &    0.5D+00,
     &    0.5D+00,
     &    0.5D+00,
     &    0.5D+00,
     &    0.5D+00,
     &    0.5D+00,
     &    0.5D+00,
     &    0.5D+00,
     &    0.5D+00,
     &    0.5D+00,
     &    0.5D+00,
     &    0.5D+00,
     &    0.5D+00,
     &   -1.0D+00,
     &   -0.8D+00,
     &   -0.6D+00,
     &   -0.4D+00,
     &   -0.2D+00,
     &    0.0D+00,
     &    0.2D+00,
     &    0.4D+00,
     &    0.6D+00,
     &    0.8D+00,
     &    1.0D+00 /

      if ( n_data < 0 ) then
        n_data = 0
      end if

      n_data = n_data + 1

      if ( n_max < n_data ) then
        n_data = 0
        n = 0
        a = 0
        b = 0
        x = 0.0D+00
        fx = 0.0D+00
      else
        n = n_vec(n_data)
        a = a_vec(n_data)
        b = b_vec(n_data)
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
