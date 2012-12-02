      program main

c*********************************************************************72
c
cc TOMS352_PRB tests TOMS352.
c
c  Modified:
c
c    06 February 2006
c
c  Author:
c
c    John Burkardt
c
      implicit none

      call timestamp ( )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TOMS352_PRB'
      write ( *, '(a)' ) '  FORTRAN77 version'
      write ( *, '(a)' ) '  Test the TOMS352 library.'

      call test01 ( )
      call test02 ( )
      call test03 ( )
      call test04 ( )
      call test05 ( )
      call test06 ( )
c
c  Terminate.
c
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TOMS352_PRB'
      write ( *, '(a)' ) '  Normal end of execution.'

      write ( *, '(a)' ) ' '
      call timestamp ( )

      stop
      end
      subroutine test01 ( )

c*********************************************************************72
c
cc TEST01 tests BESSEL.
c
c  Modified:
c
c    03 January 2006
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer order_max
      parameter ( order_max = 1 )

      double precision fx
      double precision fx2
      double precision jy(250)
      integer  n_data
      integer sol
      double precision x

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST01'
      write ( *, '(a)' ) '  Test BESSEL, which can compute'
      write ( *, '(a)' ) '  J0 Bessel functions.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '      X         Exact Value      Computed'
      write ( *, '(a)' ) ' '

      sol = 1
      n_data = 0

10    continue

      call bessel_j0_values ( n_data, x, fx )

      if ( n_data <= 0 ) then
        go to 20
      end if

      call bessel ( sol, x, jy, order_max )

      fx2 = jy(1)

      write ( *, '(2x,f8.4,2x,g16.8,2x,g16.8)' ) x, fx, fx2

      go to 10

20    continue

      return
      end
      subroutine test02 ( )

c*********************************************************************72
c
cc TEST02 tests BESSEL.
c
c  Modified:
c
c    28 January 2006
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer order_max
      parameter ( order_max = 1 )

      double precision fx
      double precision fx2
      double precision jy(250)
      integer  n_data
      integer sol
      double precision x

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST02'
      write ( *, '(a)' ) '  Test BESSEL, which can compute'
      write ( *, '(a)' ) '  J1 Bessel functions.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '      X         Exact Value      Computed'
      write ( *, '(a)' ) ' '

      sol = 1
      n_data = 0

10    continue

      call bessel_j1_values ( n_data, x, fx )

      if ( n_data <= 0 ) then
        go to 20
      end if

      call bessel ( sol, x, jy, order_max )

      fx2 = jy(2)

      write ( *, '(2x,f8.4,2x,g16.8,2x,g16.8)' ) x, fx, fx2

      go to 10

20    continue

      return
      end
      subroutine test03 ( )

c*********************************************************************72
c
cc TEST03 tests BESSEL.
c
c  Modified:
c
c    28 January 2006
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer order_max
      parameter ( order_max = 2 )

      double precision fx
      double precision fx2
      double precision jy(250)
      integer  n_data
      integer sol
      double precision x

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST03'
      write ( *, '(a)' ) '  Test BESSEL, which can compute'
      write ( *, '(a)' ) '  Y0 Bessel functions.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '      X         Exact Value      Computed'
      write ( *, '(a)' ) ' '

      sol = 2
      n_data = 0

10    continue

      call bessel_y0_values ( n_data, x, fx )

      if ( n_data <= 0 ) then
        go to 20
      end if

      call bessel ( sol, x, jy, order_max )

      fx2 = jy(1)

      write ( *, '(2x,f8.4,2x,g16.8,2x,g16.8)' ) x, fx, fx2

      go to 10

20    continue

      return
      end
      subroutine test04 ( )

c*********************************************************************72
c
cc TEST04 tests BESSEL.
c
c  Modified:
c
c    28 January 2006
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer order_max
      parameter ( order_max = 3 )

      double precision fx
      double precision fx2
      double precision jy(250)
      integer  n_data
      integer sol
      double precision x

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST04'
      write ( *, '(a)' ) '  Test BESSEL, which can compute'
      write ( *, '(a)' ) '  Y1 Bessel functions.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '      X         Exact Value      Computed'
      write ( *, '(a)' ) ' '

      sol = 2
      n_data = 0

10    continue

      call bessel_y1_values ( n_data, x, fx )

      if ( n_data <= 0 ) then
        go to 20
      end if

      call bessel ( sol, x, jy, order_max )

      fx2 = jy(2)

      write ( *, '(2x,f8.4,2x,g16.8,2x,g16.8)' ) x, fx, fx2

      go to 10

20    continue

      return
      end
      subroutine test05 ( )

c*********************************************************************72
c
cc TEST05 tests MFCVAL.
c
c  Modified:
c
c    04 February 2006
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer r_max
      parameter ( r_max = 20 )

      double precision a
      double precision a2
      double precision cv(6,r_max)
      integer j
      integer  n_data
      integer q
      double precision q_dble
      integer r

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST05'
      write ( *, '(a)' ) '  Test MFCVAL, which can compute'
      write ( *, '(a)' ) '  the eigenvalues and eigensolutions'
      write ( *, '(a)' ) '  of Mathieu''s differential equation.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Here, we work with EVEN solutions.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' )
     &  '       R       Q      Exact Value      Computed'
      write ( *, '(a)' ) ' '

      n_data = 0

10    continue

      call mathieu_even_values ( n_data, r, q, a )

      if ( n_data <= 0 ) then
        go to 20
      end if

      q_dble = dble ( q )

      call mfcval ( r+1, r, q_dble, cv, j )

      a2 = cv(1,r+1)

      write ( *, '(2x,i6,2x,i6,2x,g16.8,2x,g16.8)' ) r, q, a, a2

      go to 10

20    continue

      return
      end
      subroutine test06 ( )

c*********************************************************************72
c
cc TEST06 tests MFCVAL.
c
c  Modified:
c
c    04 February 2006
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer r_max
      parameter ( r_max = 20 )

      double precision a
      double precision a2
      double precision cv(6,r_max)
      integer j
      integer  n_data
      integer q
      double precision q_dble
      integer r

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST06'
      write ( *, '(a)' ) '  Test MFCVAL, which can compute'
      write ( *, '(a)' ) '  the eigenvalues and eigensolutions'
      write ( *, '(a)' ) '  of Mathieu''s differential equation.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Here, we work with ODD solutions.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' )
     &  '       R       Q      Exact Value      Computed'
      write ( *, '(a)' ) ' '

      n_data = 0

10    continue

      call mathieu_odd_values ( n_data, r, q, a )

      if ( n_data <= 0 ) then
        go to 20
      end if

      q_dble = dble ( q )

      call mfcval ( r+1, r+1, q_dble, cv, j )

      a2 = cv(1,r)

      write ( *, '(2x,i6,2x,i6,2x,g16.8,2x,g16.8)' ) r, q, a, a2

      go to 10

20    continue

      return
      end
      subroutine bessel_j0_values ( n_data, x, fx )

c*********************************************************************72
c
cc BESSEL_J0_VALUES returns some values of the J0 Bessel function.
c
c  Discussion:
c
c    In Mathematica, the function can be evaluated by:
c
c      BesselJ[0,x]
c
c  Modified:
c
c    14 January 2006
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
      parameter ( n_max = 21 )

      double precision fx
      double precision fx_vec(n_max)
      integer n_data
      double precision x
      double precision x_vec(n_max)

      save fx_vec
      save x_vec

      data fx_vec /
     &  -0.1775967713143383D+00,
     &  -0.3971498098638474D+00,
     &  -0.2600519549019334D+00,
     &   0.2238907791412357D+00,
     &   0.7651976865579666D+00,
     &   0.1000000000000000D+01,
     &   0.7651976865579666D+00,
     &   0.2238907791412357D+00,
     &  -0.2600519549019334D+00,
     &  -0.3971498098638474D+00,
     &  -0.1775967713143383D+00,
     &   0.1506452572509969D+00,
     &   0.3000792705195556D+00,
     &   0.1716508071375539D+00,
     &  -0.9033361118287613D-01,
     &  -0.2459357644513483D+00,
     &  -0.1711903004071961D+00,
     &   0.4768931079683354D-01,
     &   0.2069261023770678D+00,
     &   0.1710734761104587D+00,
     &  -0.1422447282678077D-01 /
      data x_vec /
     &  -5.0D+00,
     &  -4.0D+00,
     &  -3.0D+00,
     &  -2.0D+00,
     &  -1.0D+00,
     &   0.0D+00,
     &   1.0D+00,
     &   2.0D+00,
     &   3.0D+00,
     &   4.0D+00,
     &   5.0D+00,
     &   6.0D+00,
     &   7.0D+00,
     &   8.0D+00,
     &   9.0D+00,
     &  10.0D+00,
     &  11.0D+00,
     &  12.0D+00,
     &  13.0D+00,
     &  14.0D+00,
     &  15.0D+00 /

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
      subroutine bessel_j1_values ( n_data, x, fx )

c*********************************************************************72
c
cc BESSEL_J1_VALUES returns some values of the J1 Bessel function.
c
c  Discussion:
c
c    In Mathematica, the function can be evaluated by:
c
c      BesselJ[1,x]
c
c  Modified:
c
c    14 January 2006
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
      parameter ( n_max = 21 )

      double precision fx
      double precision fx_vec(n_max)
      integer n_data
      double precision x
      double precision x_vec(n_max)

      save fx_vec
      save x_vec

      data fx_vec /
     &   0.3275791375914652D+00,
     &   0.6604332802354914D-01,
     &  -0.3390589585259365D+00,
     &  -0.5767248077568734D+00,
     &  -0.4400505857449335D+00,
     &   0.0000000000000000D+00,
     &   0.4400505857449335D+00,
     &   0.5767248077568734D+00,
     &   0.3390589585259365D+00,
     &  -0.6604332802354914D-01,
     &  -0.3275791375914652D+00,
     &  -0.2766838581275656D+00,
     &  -0.4682823482345833D-02,
     &   0.2346363468539146D+00,
     &   0.2453117865733253D+00,
     &   0.4347274616886144D-01,
     &  -0.1767852989567215D+00,
     &  -0.2234471044906276D+00,
     &  -0.7031805212177837D-01,
     &   0.1333751546987933D+00,
     &   0.2051040386135228D+00 /
      data x_vec /
     &  -5.0D+00,
     &  -4.0D+00,
     &  -3.0D+00,
     &  -2.0D+00,
     &  -1.0D+00,
     &   0.0D+00,
     &   1.0D+00,
     &   2.0D+00,
     &   3.0D+00,
     &   4.0D+00,
     &   5.0D+00,
     &   6.0D+00,
     &   7.0D+00,
     &   8.0D+00,
     &   9.0D+00,
     &  10.0D+00,
     &  11.0D+00,
     &  12.0D+00,
     &  13.0D+00,
     &  14.0D+00,
     &   15.0D+00 /

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
      subroutine bessel_y0_values ( n_data, x, fx )

c*********************************************************************72
c
cc BESSEL_Y0_VALUES returns some values of the Y0 Bessel function.
c
c  Discussion:
c
c    In Mathematica, the function can be evaluated by:
c
c      BesselY[0,x]
c
c  Modified:
c
c    14 January 2006
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
     &  -0.1534238651350367D+01,
     &   0.8825696421567696D-01,
     &   0.5103756726497451D+00,
     &   0.3768500100127904D+00,
     &  -0.1694073932506499D-01,
     &  -0.3085176252490338D+00,
     &  -0.2881946839815792D+00,
     &  -0.2594974396720926D-01,
     &   0.2235214893875662D+00,
     &   0.2499366982850247D+00,
     &   0.5567116728359939D-01,
     &  -0.1688473238920795D+00,
     &  -0.2252373126343614D+00,
     &  -0.7820786452787591D-01,
     &   0.1271925685821837D+00,
     &   0.2054642960389183D+00 /
      data x_vec /
     &   0.1D+00,
     &   1.0D+00,
     &   2.0D+00,
     &   3.0D+00,
     &   4.0D+00,
     &   5.0D+00,
     &   6.0D+00,
     &   7.0D+00,
     &   8.0D+00,
     &   9.0D+00,
     &  10.0D+00,
     &  11.0D+00,
     &  12.0D+00,
     &  13.0D+00,
     &  14.0D+00,
     &  15.0D+00 /

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
      subroutine bessel_y1_values ( n_data, x, fx )

c*********************************************************************72
c
cc BESSEL_Y1_VALUES returns some values of the Y1 Bessel function.
c
c  Discussion:
c
c    In Mathematica, the function can be evaluated by:
c
c      BesselY[1,x]
c
c  Modified:
c
c    14 January 2006
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
     &  -0.6458951094702027D+01,
     &  -0.7812128213002887D+00,
     &  -0.1070324315409375D+00,
     &   0.3246744247918000D+00,
     &   0.3979257105571000D+00,
     &   0.1478631433912268D+00,
     &  -0.1750103443003983D+00,
     &  -0.3026672370241849D+00,
     &  -0.1580604617312475D+00,
     &   0.1043145751967159D+00,
     &   0.2490154242069539D+00,
     &   0.1637055374149429D+00,
     &  -0.5709921826089652D-01,
     &  -0.2100814084206935D+00,
     &  -0.1666448418561723D+00,
     &   0.2107362803687351D-01 /
      data x_vec /
     &   0.1D+00,
     &   1.0D+00,
     &   2.0D+00,
     &   3.0D+00,
     &   4.0D+00,
     &   5.0D+00,
     &   6.0D+00,
     &   7.0D+00,
     &   8.0D+00,
     &   9.0D+00,
     &  10.0D+00,
     &  11.0D+00,
     &  12.0D+00,
     &  13.0D+00,
     &  14.0D+00,
     &  15.0D+00 /

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
      subroutine mathieu_even_values ( n_data, r, q, a )

c*********************************************************************72
c
cc MATHIEU_EVEN_VALUES returns eigenvalues of even Mathieu solutions.
c
c  Discussion:
c
c    Mathieu's differential equation can be written
c
c      d2y/dx2 + ( a - 2 * q * cos ( 2 * x ) ) * y = 0
c
c    For each integer Q, there are sets of eigenvalues and
c    associated periodic eigensolutions.  We denote by A(R,Q)
c    the R-th eigenvalue associated with an even periodic
c    solution for Q, and B(R,Q) the R-th eigenvalue associated
c    with an odd periodic solution for Q.
c
c    In Mathematica, the eigenvalues for even functions can
c    be evaluated by
c
c      MathieuCharacteristicA [ r, q ]
c
c  Modified:
c
c    18 January 2006
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
c    Output, integer R, the index of the eigenvalue.
c
c    Output, integer Q, the value of the parameter Q.
c
c    Output, double precision A, the eigenvalue of the even solution
c    of Mathieu's equation, A(Q,R).
c
      implicit none

      integer n_max
      parameter ( n_max = 15 )

      double precision a
      double precision a_vec(n_max)
      integer n_data
      integer q
      integer q_vec(n_max)
      integer r
      integer r_vec(n_max)

      save a_vec
      save q_vec
      save r_vec

      data a_vec /
     &    0.0D+00,
     &    1.0D+00,
     &    4.0D+00,
     &  225.0D+00,
     &   25.0D+00,
     &   25.54997174998161D+00,
     &   27.70376873393928D+00,
     &   31.95782125217288D+00,
     &   36.64498973413284D+00,
     &   40.05019098580771D+00,
     &   -5.800046020851508D+00,
     &  -40.25677954656679D+00,
     &  -14.49130142517482D+00,
     &    5.077983197543472D+00,
     &  100.5067700246813D+00 /
      data q_vec /
     &   0,
     &   0,
     &   0,
     &   0,
     &   0,
     &   5,
     &  10,
     &  15,
     &  20,
     &  25,
     &   5,
     &  25,
     &  20,
     &  15,
     &  10 /
      data r_vec /
     &   0,
     &   1,
     &   2,
     &  15,
     &   5,
     &   5,
     &   5,
     &   5,
     &   5,
     &   5,
     &   0,
     &   0,
     &   1,
     &   2,
     &  10 /

      if ( n_data < 0 ) then
        n_data = 0
      end if

      n_data = n_data + 1

      if ( n_max < n_data ) then
        n_data = 0
        r = 0
        q = 0
        a = 0.0D+00
      else
        r = r_vec(n_data)
        q = q_vec(n_data)
        a = a_vec(n_data)
      end if

      return
      end
      subroutine mathieu_odd_values ( n_data, r, q, b )

c*********************************************************************72
c
cc MATHIEU_ODD_VALUES returns eigenvalues of odd Mathieu solutions.
c
c  Discussion:
c
c    Mathieu's differential equation can be written
c
c      d2y/dx2 + ( a - 2 * q * cos ( 2 * x ) ) * y = 0
c
c    For each integer Q, there are sets of eigenvalues and
c    associated periodic eigensolutions.  We denote by A(R,Q)
c    the R-th eigenvalue associated with an even periodic
c    solution for Q, and B(R,Q) the R-th eigenvalue associated
c    with an odd periodic solution for Q.
c
c    In Mathematica, the eigenvalues for odd functions can
c    be evaluated by
c
c      MathieuCharacteristicB [ r, q ]
c
c  Modified:
c
c    18 January 2006
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Milton Abramowitz, Irene Stegun,
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
c    Output, integer R, the index of the eigenvalue.
c
c    Output, integer Q, the value of the parameter Q.
c
c    Output, double precision B, the eigenvalue of the odd solution
c    of Mathieu's equation, B(Q,R).
c
      implicit none

      integer n_max
      parameter ( n_max = 13 )

      double precision b
      double precision b_vec(n_max)
      integer n_data
      integer q
      integer q_vec(n_max)
      integer r
      integer r_vec(n_max)

      save b_vec
      save q_vec
      save r_vec

      data b_vec /
     &    1.0D+00,
     &  -40.25677898468416D+00,
     &    4.0D+00,
     &    2.099460445486665D+00,
     &   -2.382158235956956D+00,
     &   -8.099346798895896D+00,
     &  -14.49106325598072D+00,
     &  -21.31486062224985D+00,
     &   27.96788059671747D+00,
     &  100.0D+00,
     &  100.5067694628784D+00,
     &  225.0D+00,
     &  225.8951534161768D+00 /
      data q_vec /
     &   0,
     &  25,
     &   0,
     &   5,
     &  10,
     &  15,
     &  20,
     &  25,
     &  15,
     &   0,
     &  10,
     &   0,
     &  20 /
      data r_vec /
     &  1,
     &  1,
     &  2,
     &  2,
     &  2,
     &  2,
     &  2,
     &  2,
     &  5,
     & 10,
     & 10,
     & 15,
     & 15 /

      if ( n_data < 0 ) then
        n_data = 0
      end if

      n_data = n_data + 1

      if ( n_max < n_data ) then
        n_data = 0
        r = 0
        q = 0
        b = 0.0D+00
      else
        r = r_vec(n_data)
        q = q_vec(n_data)
        b = b_vec(n_data)
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
