      program main

c*********************************************************************72
c
cc MAIN is the main program for SUBPAK_PRB.
c
c  Discussion:
c
c    SUBPAK_PRB calls sample problems for the SUBPAK library.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    03 November 2011
c
c  Author:
c
c    John Burkardt
c
      implicit none

      call timestamp ( )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'SUBPAK_PRB'
      write ( *, '(a)' ) '  FORTRAN77 version'
      write ( *, '(a)' ) '  Test the SUBPAK library.'

      call test01 ( )
      call test02 ( )
      call test03 ( )
      call test04 ( )
      call test05 ( )
      call test07 ( )
      call test08 ( )
      call test09 ( )

      call test10 ( )
      call test11 ( )
      call test12 ( )
      call test13 ( )
      call test14 ( )
      call test15 ( )
      call test16 ( )
      call test17 ( )
      call test18 ( )
      call test19 ( )

      call test22 ( )
      call test225 ( )
      call test23 ( )
      call test24 ( )
      call test25 ( )
      call test29 ( )

      call test30 ( )
      call test31 ( )
      call test32 ( )
      call test33 ( )
      call test34 ( )
      call test35 ( )
c
c  Terminate.
c
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'SUBPAK_PRB'
      write ( *, '(a)' ) '  Normal end of execution.'

      write ( *, '(a)' ) ' '
      call timestamp ( )

      stop
      end
      subroutine test01 ( )

c*********************************************************************72
c
cc TEST01 tests ANGLE_SHIFT.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    08 July 2010
c
c  Author:
c
c    John Burkardt
c
      implicit none

      double precision alpha
      double precision angle_hi
      double precision angle_lo
      double precision beta
      double precision gamma
      double precision pi
      parameter ( pi = 3.141592653589793D+00 )
      double precision r8_uniform
      integer seed
      integer test
      integer test_num
      parameter ( test_num = 10 )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST01'
      write ( *, '(a)' ) '  ANGLE_SHIFT shifts an angle by multiples of'
      write ( *, '(a)' ) '  2 Pi til it is between BETA and BETA+2Pi.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '     ALPHA      BETA     GAMMA   BETA+2Pi'
      write ( *, '(a)' ) ' '

      angle_lo = -4.0D+00 * pi
      angle_hi = +4.0D+00 * pi

      seed = 123456789

      do test = 1, test_num

        alpha = r8_uniform ( angle_lo, angle_hi, seed )

        beta = r8_uniform ( angle_lo, angle_hi, seed )

        call angle_shift ( alpha, beta, gamma )

        write ( *, '(2x,f8.1,2x,f8.1,2x,f8.1,2x,f8.1)' ) 
     &    alpha, beta, gamma, beta + 2.0D+00 * pi

      end do

      return
      end
      subroutine test02 ( )

c*********************************************************************72
c
cc TEST02 tests ANGLE_SHIFT_DEG.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    20 July 2010
c
c  Author:
c
c    John Burkardt
c
      implicit none

      double precision alpha
      double precision angle_hi
      double precision angle_lo
      double precision beta
      double precision gamma
      double precision r8_uniform
      integer seed
      integer test
      integer test_num
      parameter ( test_num = 10 )

      angle_hi = +720.0D+00
      angle_lo = -720.0D+00

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST02'
      write ( *, '(a)' ) '  ANGLE_SHIFT_DEG shifts an angle by'
      write ( *, '(a)' ) '  multiples of 360 degrees until it lies'
      write ( *, '(a)' ) '  between BETA and BETA+360.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '     ALPHA      BETA     GAMMA   BETA+360'
      write ( *, '(a)' ) ' '

      seed = 123456789

      do test = 1, test_num

        alpha = r8_uniform ( angle_lo, angle_hi, seed )

        beta = r8_uniform ( angle_lo, angle_hi, seed )

        call angle_shift_deg ( alpha, beta, gamma )

        write ( *, '(2x,f8.1,2x,f8.1,2x,f8.1,2x,f8.1)' ) 
     &    alpha, beta, gamma, beta + 360.0D+00

      end do

      return
      end
      subroutine test03 ( )

c*********************************************************************72
c
cc TEST03 tests ANGLE_TO_RGB.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    20 July 2010
c
c  Author:
c
c    John Burkardt
c
      implicit none

      double precision angle
      double precision angle_hi
      double precision angle_lo
      double precision b
      double precision g
      double precision r
      double precision r8_uniform
      integer seed
      integer test
      integer test_num

      angle_hi = 360.0D+00
      angle_lo = 0.0D+00
      test_num = 10

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST03'
      write ( *, '(a)' ) 
     &  '  ANGLE_TO_RGB converts an angle into an RGB color.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '     ANGLE        R         G         B'
      write ( *, '(a)' ) ' '

      seed = 123456789

      do test = 1, test_num

        angle = r8_uniform ( angle_lo, angle_hi, seed )

        call angle_to_rgb ( angle, r, g, b )

        write ( *, '(2x,f8.1,2x,f8.3,2x,f8.3,2x,f8.3)' ) angle, r, g, b

      end do

      return
      end
      subroutine test04 ( )

c*********************************************************************72
c
cc TEST04 tests AXIS_LIMITS.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    20 July 2010
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer ndivs
      integer nticks
      double precision pxdiv
      double precision pxmax
      double precision pxmin
      double precision xmax
      double precision xmin

      xmin = 67.3D+00
      xmax = 114.7D+00
      ndivs = 6

      call axis_limits ( xmin, xmax, ndivs, pxmin, pxmax, pxdiv, 
     &  nticks )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST04'
      write ( *, '(a)' ) '  AXIS_LIMITS adjusts plot limits'
      write ( *, '(a)' ) '  to "nicer" values.'
      write ( *, '(a)' ) ' '
      write ( *, '(a,g14.6)' ) '  Input XMIN =    ', xmin
      write ( *, '(a,g14.6)' ) '  Input XMAX =    ', xmax
      write ( *, '(a,i8)' ) '  Input NDIVS =   ', ndivs
      write ( *, '(a)' ) ' '
      write ( *, '(a,g14.6)' ) '  Output PXMIN =  ', pxmin
      write ( *, '(a,g14.6)' ) '  Output PXMAX =  ', pxmax
      write ( *, '(a,g14.6)' ) '  Output PXDIV =  ', pxdiv
      write ( *, '(a,i8)' ) '  Output NTICKS = ', nticks

      xmin = - 26.0D+00
      xmax = + 26.0D+00
      ndivs = 10

      call axis_limits ( xmin, xmax, ndivs, pxmin, pxmax, pxdiv, 
     &  nticks )

      write ( *, '(a)' ) ' '
      write ( *, '(a,g14.6)' ) '  Input XMIN =    ', xmin
      write ( *, '(a,g14.6)' ) '  Input XMAX =    ', xmax
      write ( *, '(a,i8)' ) '  Input NDIVS =   ', ndivs
      write ( *, '(a)' ) ' '
      write ( *, '(a,g14.6)' ) '  Output PXMIN =  ', pxmin
      write ( *, '(a,g14.6)' ) '  Output PXMAX =  ', pxmax
      write ( *, '(a,g14.6)' ) '  Output PXDIV =  ', pxdiv
      write ( *, '(a,i8)' ) '  Output NTICKS = ', nticks

      return
      end
      subroutine test05 ( )

c*********************************************************************72
c
cc TEST05 tests AXIS_LIMITS.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    20 July 2010
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer test_num
      parameter ( test_num = 5 )

      integer ndivs
      integer nticks
      double precision pxdiv
      double precision pxmax
      double precision pxmin
      integer test
      double precision xmax
      double precision xmax_test(test_num)
      double precision xmin
      double precision xmin_test(test_num)

      save xmax_test
      save xmin_test

      data xmax_test /
     &  9.0D+00, 4.125D+00, 193.75D+00, 2000.250D+00, 12.0D+00 /
      data xmin_test /
     &  1.0D+00, 1.003D+00, 101.25D+00, 2000.125D+00, -7.0D+00 /

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST05'
      write ( *, '(a)' ) '  AXIS_LIMITS computes "nice" limits for a'
      write ( *, '(a)' ) '  graph that must include a given range.'

      ndivs = 5

      write ( *, '(a)' ) ' '
      write ( *, '(a,i8)' ) '  All tests use NDIVS = ', ndivs
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '      XMIN        XMAX       PXMIN       '
     &  // 'PXMAX       PXDIV    NTICKS'
      write ( *, '(a)') ' '

      do test = 1, test_num

        xmin = xmin_test(test)
        xmax = xmax_test(test)

        call axis_limits ( xmin, xmax, ndivs, pxmin, pxmax, pxdiv, 
     &    nticks )

        write ( *, '(2x,5g12.4,i8)' ) 
     &    xmin, xmax, pxmin, pxmax, pxdiv, nticks

      end do

      return
      end
      subroutine test07 ( )

c*********************************************************************72
c
cc TEST07 tests BMI_ENGLISH.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    21 July 2010
c
c  Author:
c
c    John Burkardt
c
      implicit none

      double precision b
      double precision bmi
      double precision bmi_english
      double precision c
      double precision h
      double precision h_ft
      double precision h_in
      integer i
      double precision r8_uniform
      integer seed
      integer test
      integer test_num
      parameter ( test_num = 10 )
      double precision w

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST07'
      write ( *, '(a)' ) '  BMI_ENGLISH computes the Body Mass Index'
      write ( *, '(a)' ) '  given body measurements in English Units.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '      Weight         Height            BMI'
      write ( *, '(a)' ) '       (LB)     (FT           IN)'
      write ( *, '(a)' ) ' '

      seed = 123456789

      do test = 1, test_num

        b = 100.0D+00
        c = 250.0D+00

        w = r8_uniform ( b, c, seed )

        b = 4.0D+00
        c = 6.75D+00

        h = r8_uniform ( b, c, seed )

        h_ft = int ( h )
        h_in = dble ( nint ( 12.0D+00 * ( h - h_ft ) ) )

        bmi = bmi_english ( w, h_ft, h_in )
        write ( *, '(2x,4f10.2)' ) w, h_ft, h_in, bmi

      end do

      return
      end
      subroutine test08 ( )

c*********************************************************************72
c
cc TEST08 tests FAC_DIV, FAC_GCD, FAC_LCM, FAC_MUL, FAC_TO_I4, and I4_TO_FAC.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    21 July 2010
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer prime_num
      parameter ( prime_num = 5 )

      integer bot
      integer i1
      integer i2
      integer npower1(prime_num)
      integer npower2(prime_num)
      integer npower3(prime_num)
      integer top

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST08'
      write ( *, '(a)' ) '  For products of prime factors:'
      write ( *, '(a)' ) '  FAC_DIV computes a quotient;'
      write ( *, '(a)' ) '  FAC_MUL multiplies;'
      write ( *, '(a)' ) '  FAC_LCM computes the LCM;'
      write ( *, '(a)' ) '  FAC_GCD computes the GCD;'
      write ( *, '(a)' ) '  I4_TO_FAC converts an integer;'
      write ( *, '(a)' ) '  FAC_TO_I4 converts to an integer.'
      write ( *, '(a)' ) '  FAC_TO_RAT converts to a ratio.'

      i1 = 720
      i2 = 42

      call i4_to_fac ( i1, prime_num, npower1 )

      write ( *, '(a)' ) ' '
      write ( *, '(a,i8)' ) '  Representation of I1 = ', i1
      write ( *, '(a)' ) ' '

      call fac_print ( prime_num, npower1 )

      call i4_to_fac ( i2, prime_num, npower2 )

      write ( *, '(a)' ) ' '
      write ( *, '(a,i8)' ) '  Representation of I2 = ', i2
      write ( *, '(a)' ) ' '

      call fac_print ( prime_num, npower2 )

      call fac_lcm ( prime_num, npower1, npower2, npower3 )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  LCM of I1, I2:'
      write ( *, '(a)' ) ' '

      call fac_print ( prime_num, npower3 )

      call fac_gcd ( prime_num, npower1, npower2, npower3 )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  GCD of I1, I2:'
      write ( *, '(a)' ) ' '

      call fac_print ( prime_num, npower3 )

      call fac_mul ( prime_num, npower1, npower2, npower3 )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Product of I1, I2:'
      write ( *, '(a)' ) ' '

      call fac_print ( prime_num, npower3 )

      call fac_div ( prime_num, npower2, npower1, npower3 )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Quotient of I2 / I1:'
      write ( *, '(a)' ) ' '

      call fac_print ( prime_num, npower3 )

      call fac_to_rat ( prime_num, npower3, top, bot )

      write ( *, '(a)' ) ' '
      write ( *, '(a,i8,a,i8)' ) 
     &  '  Quotient as a rational: ', top, ' / ', bot

      return
      end
      subroutine test09 ( )

c*********************************************************************72
c
cc TEST09 tests GAUSS_SUM
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    21 July 2010
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer dim_num
      parameter ( dim_num = 2 )
      integer n
      parameter ( n = 3 )

      double precision amplitude(n)
      double precision center(dim_num,n)
      double precision gauss_sum
      double precision gxy
      integer i
      integer j
      double precision width(n)
      double precision x(dim_num)

      save amplitude
      save center
      save width

      data amplitude /
     &  10.0D+00, 5.0D+00, -3.0D+00 /
      data center /
     &  2.0D+00, 3.0D+00, 
     &  5.0D+00, 8.0D+00, 
     &  7.0D+00, 5.0D+00 /
      data width / 1.0D+00, 2.0D+00, 4.0D+00 /

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST09'
      write ( *, '(a)' ) '  GAUSS_SUM evaluates a function which is'
      write ( *, '(a)' ) '  the sum of Gaussian functions.'
      write ( *, '(a)' ) ' '
      write ( *, '(a,i8)' ) '  Number of component Gaussians = ', n
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '          Center    Amplitude  Width'
      write ( *, '(a)' ) '        X       Y'
      write ( *, '(a)' ) ' '
      do j = 1, n
        write ( *, '(2x,i2,2x,f6.2,2x,f6.2,2x,f6.2,2x,f6.2)' ) 
     &    j, center(1:2,j), amplitude(j), width(j)
      end do

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '      X       Y        Gauss_Sum(X,Y)'
      write ( *, '(a)' ) ' '

      do i = 0, 10
        x(1) = real ( i, kind = 8 )
        do j = 0, 10
         x(2) = real ( j, kind = 8 )
          gxy = gauss_sum ( dim_num, n, amplitude, center, width, x )
          write ( *, '(2x,f6.2,2x,f6.2,2x,g14.6)' ) x(1), x(2), gxy
        end do
      end do

      return
      end
      subroutine test10 ( )

c*********************************************************************72
c
cc TEST10 tests GET_SEED.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    21 July 2010
c
c  Author:
c
c    John Burkardt
c
      implicit none

      double precision r8_uniform_01
      integer seed
      integer seed_0
      integer seed_1
      integer seed_2
      integer seed_3
      integer test
      integer test_num
      parameter ( test_num = 10 )
      double precision x

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST10'
      write ( *, '(a)' ) '  GET_SEED gets a seed for the random number'
      write ( *, '(a)' ) '  generator.  These values are computed '
      write ( *, '(a)' ) '  from the time and date.  Values computed'
      write ( *, '(a)' ) '  nearby in time will be near to each'
      write ( *, '(a)' ) '  other, and should be passed through a '
      write ( *, '(a)' ) '  random number generator a few times '
      write ( *, '(a)' ) '  before use.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '     I      R(I)  R2(I)        R3(I)'
      write ( *, '(a)' ) ' '

      do test = 1, test_num
        call get_seed ( seed )
        seed_0 = seed
        x = r8_uniform_01 ( seed )
        seed_1 = seed
        x = r8_uniform_01 ( seed )
        seed_2 = seed
        x = r8_uniform_01 ( seed )
        seed_3 = seed
        write ( *, '(2x,4i12)' ) seed_0, seed_1, seed_2, seed_3
      end do
 
      return
      end
      subroutine test11 ( )

c*********************************************************************72
c
cc TEST11 tests GRID1.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    22 July 2010
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer dim_num
      parameter ( dim_num = 5 )
      integer nstep 
      parameter ( nstep = 11 )

      double precision x(dim_num,nstep)
      double precision x1(dim_num)
      double precision x2(dim_num)

      save x1
      save x2

      data x1 /
     &  1.0D+00,  0.0D+00, 20.0D+00, -5.0D+00, 1.0D+00 /
      data x2 /
     &  1.0D+00, 10.0D+00,  0.0D+00,  5.0D+00, 2.0D+00 /

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST11'
      write ( *, '(a)' ) '  GRID1 computes a 1D grid between'
      write ( *, '(a)' ) '  two DIM_NUM dimensional points X1 and X2.'
      write ( *, '(a)' ) ' '
      write ( *, '(a,i8,a)' ) '  Here, we will use ', nstep, ' steps'
      write ( *, '(a)' ) '  going from '
      write ( *, '(2x,5g12.4)' ) x1(1:dim_num)
      write ( *, '(a)' ) '  to'
      write ( *, '(2x,5g12.4)' ) x2(1:dim_num)

      call grid1 ( dim_num, nstep, x1, x2, x )

      call r8mat_transpose_print ( dim_num, nstep, x, 
     &  '  The grid matrix:' )

      return
      end
      subroutine test12 ( )

c*********************************************************************72
c
cc TEST12 tests GRID1N.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    22 July 2010
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer dim_num
      parameter ( dim_num = 5 )
      integer nstep
      parameter ( nstep = 11 )

      integer i
      double precision x(dim_num)
      double precision x1(dim_num)
      double precision x2(dim_num)

      save x1
      save x2

      data x1 /
     &  1.0D+00,  0.0D+00, 20.0D+00, -5.0D+00, 1.0D+00 /
      data x2 /
     &  1.0D+00, 10.0D+00,  0.0D+00,  5.0D+00, 2.0D+00 /

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST12'
      write ( *, '(a)' ) '  GRID1N computes a 1D grid between'
      write ( *, '(a)' ) '  two DIM_NUM dimensional points X1 and X2,'
      write ( *, '(a)' ) '  one point at a time.'
      write ( *, '(a)' ) ' '
      write ( *, '(a,i8,a)' ) '  Here, we will use ', nstep, ' steps'
      write ( *, '(a)' ) '  going from '
      write ( *, '(2x,5g12.4)' ) x1(1:dim_num)
      write ( *, '(a)' ) '  to'
      write ( *, '(2x,5g12.4)' ) x2(1:dim_num)
      write ( *, '(a)') ' '

      do i = 1, nstep
        call grid1n ( i, dim_num, nstep, x1, x2, x )
        write ( *, '(2x,i3,5g12.4)' ) i, x(1:dim_num)
      end do

      return
      end
      subroutine test13 ( )

c*********************************************************************72
c
cc TEST13 tests GRID2.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    22 July 2010
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer dim_num
      parameter ( dim_num = 5 )
      integer nstep
      parameter ( nstep = 20 )

      integer j1
      integer j2
      double precision x(dim_num,nstep)
      double precision x1(dim_num)
      double precision x2(dim_num)

      save x1
      save x2

      data x1 /
     &  1.0D+00,  0.0D+00, 20.0D+00, -5.0D+00, 1.0D+00 /
      data x2 /
     &  1.0D+00, 10.0D+00,  0.0D+00,  5.0D+00, 2.0D+00 /
      j1 = 3
      j2 = 13

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST13'
      write ( *, '(a)' ) '  GRID2 computes a 1 D grid between'
      write ( *, '(a)' ) '  two DIM_NUM dimensional points X1 and X2,'
      write ( *, '(a)' ) 
     &  '  computing X1 and X2 at user specified times.'
      write ( *, '(a)' ) ' '
      write ( *, '(a,i8,a)' ) '  Here, we will use ', nstep, ' steps,'
      write ( *, '(a,i8,a)' ) '  and on step ', j1, ' we will compute'
      write ( *,'(2x,5g12.4)' ) x1(1:dim_num)
      write ( *, '(a,i8,a)' ) '  and on step ', j2, ' we will compute'
      write ( *,'(2x,5g12.4)' ) x2(1:dim_num)

      call grid2 ( j1, j2, dim_num, nstep, x1, x2, x )

      call r8mat_print ( dim_num, nstep, x, '  The grid matrix:' )

      return
      end
      subroutine test14 ( )

c*********************************************************************72
c
cc TEST14 tests GRID2N.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    22 July 2010
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer dim_num
      parameter ( dim_num = 5 )

      integer j
      integer j1
      integer j2
      double precision x(dim_num)
      double precision x1(dim_num)
      double precision x2(dim_num)

      save x1
      save x2

      data x1 /
     &  1.0D+00,  0.0D+00, 20.0D+00, -5.0D+00, 1.0D+00 /
      data x2 /
     &  1.0D+00, 10.0D+00,  0.0D+00,  5.0D+00, 2.0D+00 /

      j1 = 3
      j2 = 13

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST14'
      write ( *, '(a)' ) '  GRID2N computes points from a 1D grid'
      write ( *, '(a)' ) '  between two DIM_NUM dimensional points'
      write ( *, '(a)' ) '  X1 and X2, one at a time, with X1 and X2'
      write ( *, '(a)' ) '  having user specified J coordinates.'
      write ( *, '(a)' ) ' '
      write ( *, '(a,i8,a)' ) '  Here, on step ', j1, ' we compute'
      write ( *, '(2x,5g12.4)' ) x1(1:dim_num)
      write ( *, '(a,i8,a)' ) '  and on step ', j2, ' we would compute'
      write ( *, '(2x,5g12.4)' ) x2(1:dim_num)
      write ( *, '(a)' ) ' '

      do j = 1, 20
        call grid2n ( j, j1, j2, dim_num, x1, x2, x )
        write ( *, '(2x,i3,5g12.4)' ) j, x(1:dim_num)
      end do

      return
      end
      subroutine test15 ( )

c*********************************************************************72
c
cc TEST15 tests GRID3.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    22 July 2010
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer dim_num
      parameter ( dim_num = 5 )
      integer nstep1
      parameter ( nstep1 = 3 )
      integer nstep2
      parameter ( nstep2 = 6 )

      integer i
      integer j
      double precision x(dim_num,nstep1,nstep2)
      double precision x1(dim_num)
      double precision x2(dim_num)
      double precision x3(dim_num)

      save x1
      save x2
      save x3

      data x1 /
     &  1.0D+00,  0.0D+00, 20.0D+00, -5.0D+00, 1.0D+00 /
      data x2 /
     &  1.0D+00, 10.0D+00,  0.0D+00,  5.0D+00, 2.0D+00 /
      data x3 /
     &  1.0D+00,  5.0D+00,  0.0D+00,  0.0D+00, 3.0D+00 /

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST15'
      write ( *, '(a)' ) '  GRID3 computes a 2D grid in the plane'
      write ( *, '(a)' ) '  containing the DIM_NUM-dimensional'
      write ( *, '(a)' ) '  points X1, X2 and X3.'
      write ( *, '(a)' ) ' '
      write ( *, '(a,i8,a)' ) '  Here, we will use ', nstep1, ' steps'
      write ( *, '(a)' ) '  going from '
      write ( *, '(2x,5g12.4)' ) x1(1:dim_num)
      write ( *, '(a)' ) '  to'
      write ( *, '(2x,5g12.4)' ) x2(1:dim_num)
      write ( *, '(a,i8,a)' ) '  and ', nstep2,' steps going to '
      write ( *, '(2x,5g12.4)' ) x3(1:dim_num)

      call grid3 ( dim_num, nstep1, nstep2, x1, x2, x3, x )

      do i = 1, nstep1
        write ( *, '(a)' ) ' '
        do j = 1, nstep2

          write ( *, '(2x,i3,i3,5g12.4)' ) i, j, x(1:dim_num,i,j)

        end do
      end do

      return
      end
      subroutine test16 ( )

c*********************************************************************72
c
cc TEST16 tests GRID3N.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    22 July 2010
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer dim_num
      parameter ( dim_num = 5 )
      integer nstep1
      parameter ( nstep1 = 3 )
      integer nstep2
      parameter ( nstep2 = 6 )

      integer j
      integer k
      double precision x(dim_num)
      double precision x1(dim_num)
      double precision x2(dim_num)
      double precision x3(dim_num)

      save x1
      save x2
      save x3

      data x1 /
     &  1.0D+00,  0.0D+00, 20.0D+00, -5.0D+00, 1.0D+00 /
      data x2 /
     &  1.0D+00, 10.0D+00,  0.0D+00,  5.0D+00, 2.0D+00 /
      data x3 /
     &  1.0D+00,  5.0D+00,  0.0D+00,  0.0D+00, 3.0D+00 /

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST16'
      write ( *, '(a)' ) '  GRID3N computes a point from a 2D'
      write ( *, '(a)' ) '  grid in the plane containing the '
      write ( *, '(a)' ) '  DIM_NUM-dimensional points X1, X2 and X3.'
      write ( *, '(a)' ) ' '
      write ( *, '(a,i8,a)' ) '  We use ', nstep1, ' steps from '
      write ( *, '(2x,5g12.4)' ) x1(1:dim_num)
      write ( *, '(a)' ) '  to'
      write ( *, '(2x,5g12.4)' ) x2(1:dim_num)
      write ( *, '(a,i8,a)' ) '  and ', nstep2, ' steps going to '
      write ( *, '(2x,5g12.4)' ) x3(1:dim_num)

      do j = 1, nstep1
        write ( *, '(a)' ) ' '
        do k = 1, nstep2

          call grid3n ( j, k, dim_num, nstep1, nstep2, x1, x2, x3, x )
          write ( *, '(2x,i3,i3,5g12.4)' ) j, k, x(1:dim_num)

        end do
      end do

      return
      end
      subroutine test17 ( )

c*********************************************************************72
c
cc TEST17 tests GRID4.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    22 July 2010
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer dim_num
      parameter ( dim_num = 5 )
      integer nstep1
      parameter ( nstep1 = 6 )
      integer nstep2
      parameter ( nstep2 = 10 )

      integer j
      integer j1
      integer j2
      integer k
      integer k1
      integer k2
      double precision x(dim_num,nstep1,nstep2)
      double precision x1(dim_num)
      double precision x2(dim_num)
      double precision x3(dim_num)

      save x1
      save x2
      save x3

      data x1 /
     &  1.0D+00,  0.0D+00, 20.0D+00, -5.0D+00, 1.0D+00 /
      data x2 /
     &  1.0D+00, 10.0D+00,  0.0D+00,  5.0D+00, 2.0D+00 /
      data x3 /
     &  1.0D+00,  5.0D+00,  0.0D+00,  0.0D+00, 3.0D+00 /
      j1 = 2
      j2 = 5
      k1 = 3
      k2 = 9

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST17'
      write ( *, '(a)' ) '  GRID4 computes a 2D planar grid'
      write ( *, '(a)' ) '  containing the DIM_NUM-dimensional'
      write ( *, '(a)' ) '  points X1, X2 and X3.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 
     &  '  We compute the points on the following steps:'
      write ( *, '(a)' ) ' '
      write ( *, '(a,2i8)' ) '  X1 on step ', j1, k1
      write ( *, '(a,2i8)' ) '  X2 on step ', j2, k1
      write ( *, '(a,2i8)' ) '  X3 on step ', j1, k2
      write ( *, '(a)' ) ' '
      write ( *, '(a,i8,a)' ) 
     &  '  We use ', nstep1, ' steps in the J direction'
      write ( *, '(a,i8,a)' ) 
     &  '  and ', nstep2, ' steps in the K direction.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  The points X1, X2 and X3 are:'
      write ( *, '(a)' ) ' '
      write ( *, '(2x,5g12.4)' ) x1(1:dim_num)
      write ( *, '(a)' ) ' '
      write ( *, '(2x,5g12.4)' ) x2(1:dim_num)
      write ( *, '(a)' ) ' '
      write ( *, '(2x,5g12.4)' ) x3(1:dim_num)

      call grid4 ( j1, j2, k1, k2, dim_num, nstep1, nstep2, x1, x2, 
     &  x3, x )

      do j = 1, nstep1
        write ( *, '(a)' ) ' '
        do k = 1, nstep2

          write ( *, '(2x,i3,i3,5g12.4)' ) j, k, x(1:dim_num,j,k)

        end do
      end do

      return
      end
      subroutine test18 ( )

c*********************************************************************72
c
cc TEST18 tests GRID4N.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    22 July 2010
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer dim_num
      parameter ( dim_num = 5 )
      integer nstep1
      parameter ( nstep1 = 6 )
      integer nstep2
      parameter ( nstep2 = 10 )

      integer j
      integer j1
      integer j2
      integer k
      integer k1
      integer k2
      double precision x(dim_num)
      double precision x1(dim_num)
      double precision x2(dim_num)
      double precision x3(dim_num)

      save x1
      save x2
      save x3

      data x1 /
     &  1.0D+00,  0.0D+00, 20.0D+00, -5.0D+00, 1.0D+00 /
      data x2 /
     &  1.0D+00, 10.0D+00,  0.0D+00,  5.0D+00, 2.0D+00 /
      data x3 /
     &  1.0D+00,  5.0D+00,  0.0D+00,  0.0D+00, 3.0D+00 /

      j1 = 2
      j2 = 5
      k1 = 3
      k2 = 9

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST18'
      write ( *, '(a)' ) '  GRID4N computes, one at a time, points'
      write ( *, '(a)' ) '  on a 2D grid in the plane containing'
      write ( *, '(a)' ) 
     &  '  the DIM_NUM-dimensional points X1, X2 and X3.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 
     &  '  We wish to compute the points on the following'
      write ( *, '(a)' ) '  steps:'
      write ( *, '(a)' ) ' '
      write ( *, '(a,2i8)' ) '  X1 on step ', j1, k1
      write ( *, '(a,2i8)' ) '  X2 on step ', j2, k1
      write ( *, '(a,2i8)' ) '  X3 on step ', j1, k2
      write ( *, '(a)' ) ' '
      write ( *, '(a,i8,a)' ) 
     &  '  We use ', nstep1, ' steps in the J direction'
      write ( *, '(a,i8,a)' ) 
     &  '  and ', nstep2, ' steps in the K direction.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  The points X1, X2 and X3 are:'
      write ( *, '(a)' ) ' '
      write ( *, '(2x,5g12.4)' ) x1(1:dim_num)
      write ( *, '(a)' ) ' '
      write ( *, '(2x,5g12.4)' ) x2(1:dim_num)
      write ( *, '(a)' ) ' '
      write ( *, '(2x,5g12.4)' ) x3(1:dim_num)

      do j = 1, nstep1
        write ( *, '(a)' ) ' '
        do k = 1, nstep2

          call grid4n ( j, j1, j2, k, k1, k2, dim_num, nstep1, nstep2, 
     &      x1, x2, x3, x )

          write ( *, '(2x,i3,i3,5g12.4)' ) j, k, x(1:dim_num)

        end do
      end do

      return
      end
      subroutine test19 ( )

c*********************************************************************72
c
cc TEST19 tests INDEX1_COL, INDEX1_ROW, and related functions.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    03 April 2010
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n_max
      parameter ( n_max = 4 )

      integer i
      integer i_max
      integer i_min
      integer in(n_max)
      integer in_max(n_max)
      integer in_min(n_max)
      integer index1_col
      integer index1_row
      integer index2_col
      integer index2_row
      integer index3_col
      integer index3_row
      integer index4_col
      integer index4_row
      integer indexn_col
      integer indexn_row
      integer index_min
      integer j
      integer j_max
      integer j_min
      integer k
      integer k_max
      integer k_min
      integer l
      integer l_max
      integer l_min
      integer n
      integer value

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST19'
      write ( *, '(a)' ) '  INDEX1_COL column indexes a 1D array,'
      write ( *, '(a)' ) '  INDEX1_ROW row indexes a 1D array,'
      write ( *, '(a)' ) 
     &  '  and there are several more versions of these functions.'

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  By COLS:'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Imin     I  Imax  Xmin Index'
      write ( *, '(a)' ) ' '

      i_min = 1
      i = 3
      i_max = 5
      index_min = 0

      value = index1_col ( i_min, i, i_max, index_min )
      write ( *, '(a)' ) ' '
      write ( *, '(2x,i4,2x,i4,2x,i4)' ) i_min, i, i_max
      write ( *, '(a18,  2x,i4,2x,i4)' ) 'INDEX1_COL', index_min, value

      n = 1
      in_min(1) = 1
      in(1) = 3
      in_max(1) = 5
      index_min = 0
      value = indexn_col ( n, in_min, in, in_max, index_min )
      write ( *, '(a18,  2x,i4,2x,i4)' ) 'INDEXN_COL', index_min, value

      i_min = 1
      i = 3
      i_max = 5
      j_min = 1
      j = 2
      j_max = 4
      index_min = 0
      value = index2_col ( i_min, i, i_max, j_min, j, j_max, index_min )
      write ( *, '(a)' ) ' '
      write ( *, '(2x,i4,2x,i4,2x,i4)' ) i_min, i, i_max
      write ( *, '(2x,i4,2x,i4,2x,i4)' ) j_min, j, j_max
      write ( *, '(a18,  2x,i4,2x,i4)' ) 'INDEX2_COL', index_min, value

      n = 2
      in_min(1) = 1
      in(1) = 3
      in_max(1) = 5
      in_min(2) = 1
      in(2) = 2
      in_max(2) = 4
      index_min = 0
      value = indexn_col ( n, in_min, in, in_max, index_min )
      write ( *, '(a18,  2x,i4,2x,i4)' ) 'INDEXN_COL', index_min, value

      i_min = 1
      i = 3
      i_max = 5
      j_min = 1
      j = 2
      j_max = 4
      k_min = 1
      k = 1
      k_max = 3
      index_min = 0
      value = index3_col ( i_min, i, i_max, j_min, j, j_max, k_min, k, 
     &  k_max, index_min )
      write ( *, '(a)' ) ' '
      write ( *, '(2x,i4,2x,i4,2x,i4)' ) i_min, i, i_max
      write ( *, '(2x,i4,2x,i4,2x,i4)' ) j_min, j, j_max
      write ( *, '(2x,i4,2x,i4,2x,i4)' ) k_min, k, k_max
      write ( *, '(a18,  2x,i4,2x,i4)' ) 'INDEX3_COL', index_min, value

      n = 3
      in_min(1) = 1
      in(1) = 3
      in_max(1) = 5
      in_min(2) = 1
      in(2) = 2
      in_max(2) = 4
      in_min(3) = 1
      in(3) = 1
      in_max(3) = 3
      index_min = 0
      value = indexn_col ( n, in_min, in, in_max, index_min )
      write ( *, '(a18,  2x,i4,2x,i4)' ) 'INDEXN_COL', index_min, value

      i_min = 1
      i = 3
      i_max = 5
      j_min = 1
      j = 2
      j_max = 4
      k_min = 1
      k = 1
      k_max = 3
      l_min = 1
      l = 2
      l_max = 2
      index_min = 0
      value = index4_col ( i_min, i, i_max, j_min, j, j_max, k_min, k, 
     &  k_max, l_min, l, l_max, index_min )
      write ( *, '(a)' ) ' '
      write ( *, '(2x,i4,2x,i4,2x,i4)' ) i_min, i, i_max
      write ( *, '(2x,i4,2x,i4,2x,i4)' ) j_min, j, j_max
      write ( *, '(2x,i4,2x,i4,2x,i4)' ) k_min, k, k_max
      write ( *, '(2x,i4,2x,i4,2x,i4)' ) l_min, l, l_max
      write ( *, '(a18,  2x,i4,2x,i4)' ) 'INDEX4_COL', index_min, value

      n = 4
      in_min(1) = 1
      in(1) = 3
      in_max(1) = 5
      in_min(2) = 1
      in(2) = 2
      in_max(2) = 4
      in_min(3) = 1
      in(3) = 1
      in_max(3) = 3
      in_min(4) = 1
      in(4) = 2
      in_max(4) = 2
      index_min = 0
      value = indexn_col ( n, in_min, in, in_max, index_min )
      write ( *, '(a18,  2x,i4,2x,i4)' ) 'INDEXN_COL', index_min, value

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  By ROWS:'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Imin     I  Imax  Xmin Index'
      write ( *, '(a)' ) ' '

      i_min = 1
      i = 3
      i_max = 5
      index_min = 0
      value = index1_row ( i_min, i, i_max, index_min )
      write ( *, '(a)' ) ' '
      write ( *, '(2x,i4,2x,i4,2x,i4)' ) i_min, i, i_max
      write ( *, '(a18,  2x,i4,2x,i4)' ) 'INDEX1_ROW', index_min, value

      n = 1
      in_min(1) = 1
      in(1) = 3
      in_max(1) = 5
      index_min = 0
      value = indexn_row ( n, in_min, in, in_max, index_min )
      write ( *, '(a18,  2x,i4,2x,i4)' ) 'INDEXN_ROW', index_min, value

      i_min = 1
      i = 3
      i_max = 5
      j_min = 1
      j = 2
      j_max = 4
      index_min = 0
      value = index2_row ( i_min, i, i_max, j_min, j, j_max, index_min )
      write ( *, '(a)' ) ' '
      write ( *, '(2x,i4,2x,i4,2x,i4)' ) i_min, i, i_max
      write ( *, '(2x,i4,2x,i4,2x,i4)' ) j_min, j, j_max
      write ( *, '(a18,  2x,i4,2x,i4)' ) 'INDEX2_ROW', index_min, value

      n = 2
      in_min(1) = 1
      in(1) = 3
      in_max(1) = 5
      in_min(2) = 1
      in(2) = 2
      in_max(2) = 4
      index_min = 0
      value = indexn_row ( n, in_min, in, in_max, index_min )
      write ( *, '(a18,  2x,i4,2x,i4)' ) 'INDEXN_ROW', index_min, value

      i_min = 1
      i = 3
      i_max = 5
      j_min = 1
      j = 2
      j_max = 4
      k_min = 1
      k = 1
      k_max = 3
      index_min = 0
      value = index3_row ( i_min, i, i_max, j_min, j, j_max, k_min, k, 
     &  k_max, index_min )
      write ( *, '(a)' ) ' '
      write ( *, '(2x,i4,2x,i4,2x,i4)' ) i_min, i, i_max
      write ( *, '(2x,i4,2x,i4,2x,i4)' ) j_min, j, j_max
      write ( *, '(2x,i4,2x,i4,2x,i4)' ) k_min, k, k_max
      write ( *, '(a18,  2x,i4,2x,i4)' ) 'INDEX3_ROW', index_min, value

      n = 3
      in_min(1) = 1
      in(1) = 3
      in_max(1) = 5
      in_min(2) = 1
      in(2) = 2
      in_max(2) = 4
      in_min(3) = 1
      in(3) = 1
      in_max(3) = 3
      index_min = 0
      value = indexn_row ( n, in_min, in, in_max, index_min )
      write ( *, '(a18,  2x,i4,2x,i4)' ) 'INDEXN_ROW', index_min, value

      i_min = 1
      i = 3
      i_max = 5
      j_min = 1
      j = 2
      j_max = 4
      k_min = 1
      k = 1
      k_max = 3
      l_min = 1
      l = 2
      l_max = 2
      index_min = 0
      value = index4_row ( i_min, i, i_max, j_min, j, j_max, k_min, k, 
     &  k_max, l_min, l, l_max, index_min )
      write ( *, '(a)' ) ' '
      write ( *, '(2x,i4,2x,i4,2x,i4)' ) i_min, i, i_max
      write ( *, '(2x,i4,2x,i4,2x,i4)' ) j_min, j, j_max
      write ( *, '(2x,i4,2x,i4,2x,i4)' ) k_min, k, k_max
      write ( *, '(2x,i4,2x,i4,2x,i4)' ) l_min, l, l_max
      write ( *, '(a18,  2x,i4,2x,i4)' ) 'INDEX4_ROW', index_min, value

      n = 4
      in_min(1) = 1
      in(1) = 3
      in_max(1) = 5
      in_min(2) = 1
      in(2) = 2
      in_max(2) = 4
      in_min(3) = 1
      in(3) = 1
      in_max(3) = 3
      in_min(4) = 1
      in(4) = 2
      in_max(4) = 2
      index_min = 0
      value = indexn_row ( n, in_min, in, in_max, index_min )
      write ( *, '(a18,  2x,i4,2x,i4)' ) 'INDEXN_ROW', index_min, value

      return
      end
      subroutine test22 ( )

c*********************************************************************72
c
cc TEST22 tests LCM_12N.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    22 July 2010
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer lcm_12n
      integer n

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST22'
      write ( *, '(a)' ) '  LCM_12N computes the least common multiple '
      write ( *, '(a)' ) '  of the integers 1 through N.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  N       LCM_12N ( N )'
      write ( *, '(a)' ) ' '
      do n = 1, 12
        write ( *, '(2x,i3,2x,i8)' ) n, lcm_12n ( n )
      end do

      return
      end
      subroutine test225 ( )

c*********************************************************************72
c
cc TEST225 tests LMAT_PRINT.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    03 November 2011
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer m
      parameter ( m = 20 )
      integer n
      parameter ( n = 50 )

      logical a(m,n)
      integer i
      integer j

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST225'
      write ( *, '(a)' ) '  LMAT_PRINT prints a logical matrix.'

      do i = 1, m
        do j = 1, n
          a(i,j) = ( mod ( i, j ) .eq. 0 )
        end do
      end do

      call lmat_print ( m, n, a, '  A(I,J) = I is divisible by J' )

      return
      end
      subroutine test23 ( )

c*********************************************************************72
c
cc TEST23 tests LUHN_CHECK.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    22 July 2010
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer test_num
      parameter ( test_num = 4 )

      integer check_sum
      integer check_sum_test(test_num)
      integer digit(15)
      integer digit_num
      integer digit_num_test(test_num)
      integer digit_test(32)
      integer i
      integer k
      integer test

      save check_sum_test
      save digit_num_test
      save digit_test

      data check_sum_test /
     &   6, 20, 40, 80 /
      data digit_num_test /
     &   4, 4, 9, 15 /
      data digit_test /
     &  1, 1, 1, 1, 
     &  8, 7, 6, 3, 
     &  4, 4, 6, 6, 6, 7, 6, 5, 1, 
     &  3, 7, 7, 9, 5, 6, 5, 7, 0, 9, 4, 4, 7, 2, 6 /

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST23'
      write ( *, '(a)' ) '  LUHN_CHECK computes the Luhn checksum'
      write ( *, '(a)' ) '  for a string of digits.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 
     &  '  A correct string has a checksum divisible by 10.'

      k = 0

      do test = 1, test_num

        digit_num = digit_num_test(test)

        do i = 1, digit_num
          k = k + 1
          digit(i) = digit_test(k)
        end do

        call luhn_check ( digit_num, digit, check_sum )

        write ( *, '(a)' ) ' '
        write ( *, '(a,i8)' ) '  Test number ', test
        write ( *, '(a,i8)' ) '  Number of digits = ', digit_num
        write ( *, '(a,20i1)' ) '  Digits = ', digit(1:digit_num)
        write ( *, '(a,i8)' ) '  Computed check sum = ', check_sum
        write ( *, '(a,i8)' ) 
     &    '  Correct check sum =  ', check_sum_test(test)

      end do

      return
      end
      subroutine test24 ( )

c*********************************************************************72
c
cc TEST24 tests PERM_INVERSE;
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    22 July 2010
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n
      parameter ( n = 7 )

      integer p(n)

      save p

      data p / 4, 3, 5, 1, 7, 6, 2 /

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST24'
      write ( *, '(a)' ) 
     &  '  PERM_INVERSE inverts a permutation in place;'

      call perm_print ( n, p, '  The original permutation:' )

      call perm_inverse ( n, p )

      call perm_print ( n, p, '  The inverted permutation:' )

      return
      end
      subroutine test25 ( )

c*********************************************************************72
c
cc TEST25 tests PRIME_GE.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    22 July 2010
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n
      integer p
      integer prime_ge

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST25'
      write ( *, '(a)' ) '  PRIME_GE returns the smallest prime number '
      write ( *, '(a)' ) '  greater than or equal to N.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '     N     P'
      write ( *, '(a)' ) ' '

      do n = 1, 20

        p = prime_ge ( n )
        write ( *, '(2i8)' ) n, p

      end do

      return
      end
      subroutine test29 ( )

c*********************************************************************72
c
cc TEST29 tests RAT_FACTOR.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    22 July 2010
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer factor_max
      parameter ( factor_max = 10 )

      integer factor(factor_max)
      integer factor_num
      integer i
      integer m
      integer mleft
      integer n
      integer nleft
      integer power(factor_max)

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST29'
      write ( *, '(a)' ) '  RAT_FACTOR factors a rational value.'

      m = 13 * 7 * 9 * 2
      n = 12

      write ( *, '(a)' ) ' '
      write ( *, '(a,i8,a,i8)' ) '  Rational value is ', m, '/', n

      call rat_factor ( m, n, factor_max, factor_num, factor, power,
     &  mleft, nleft )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Prime representation:'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  I, FACTOR(I), POWER(I)'
      write ( *, '(a)' ) ' '

      if ( mleft /= 1 .or. nleft /= 1 ) then
        write ( *, '(i8,i8,a,i8,a)' ) 0, mleft, ' / ', nleft, 
     &    ' (UNFACTORED PORTION)'
      end if

      do i = 1, factor_num
        write ( *, '(2x,3i8)' ) i, factor(i), power(i)
      end do

      return
      end
      subroutine test30 ( )

c*********************************************************************72
c
cc TEST30 tests ROOTS_TO_R8POLY.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    22 July 2010
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n
      parameter ( n = 4 )

      double precision c(0:n)
      double precision x(n)

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST30'
      write ( *, '(a)' ) 
     &  '  ROOTS_TO_R8POLY computes the coefficients of'
      write ( *, '(a)' ) '  a polynomial from its roots.'
      write ( *, '(a)' ) '  R8POLY_PRINT prints a polynomial.'

      call r8vec_indicator ( n, x )

      call r8vec_print ( n, x, '  Roots:' )

      call roots_to_r8poly ( n, x, c )

      call r8poly_print ( n, c, '  The polynomial' )

      return
      end
      subroutine test31 ( )

c*********************************************************************72
c
cc TEST31 tests SORT_HEAP_EXTERNAL.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    22 July 2010
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n
      parameter ( n = 20 )

      integer a(n)
      integer b
      integer c
      integer i
      integer indx
      integer isgn
      integer j
      integer seed
      integer t

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST31'
      write ( *, '(a)' ) 
     &  '  SORT_HEAP_EXTERNAL sorts objects externally.'

      indx = 0
      i = 0
      j = 0
      isgn = 0

      b = 1
      c = n
      seed = 123456789

      call i4vec_uniform ( n, b, c, seed, a )

      call i4vec_print ( n, a, '  Unsorted array:' )

10    continue

        call sort_heap_external ( n, indx, i, j, isgn )

        if ( indx .lt. 0 ) then

          isgn = 1
          if ( a(i) .le. a(j) ) then
            isgn = -1
          end if

        else if ( 0 .lt. indx ) then

          t = a(i)
          a(i) = a(j)
          a(j) = t

        else

          go to 20

        end if

      go to 10

20    continue

      call i4vec_print ( n, a, '  Sorted array:' )

      return
      end
      subroutine test32 ( )

c*********************************************************************72
c
cc TEST32 tests TVEC_EVEN, TVEC_EVEN2 and TVEC_EVEN3.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    22 July 2010
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer maxt
      parameter ( maxt = 5 )

      integer nt
      double precision t(maxt)

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST32'
      write ( *, '(a)' ) 
     &  '  For evenly spaced angles between 0 and 2*PI:'
      write ( *, '(a)' ) '  TVEC_EVEN'
      write ( *, '(a)' ) '  TVEC_EVEN2'
      write ( *, '(a)' ) '  TVEC_EVEN3'

      nt = 4

      call tvec_even ( nt, t )

      call r8vec_print ( nt, t, '  TVEC_EVEN' )

      nt = 4

      call tvec_even2 ( nt, t )

      call r8vec_print ( nt, t, '  TVEC_EVEN2' )

      nt = 4

      call tvec_even3 ( nt, t )

      call r8vec_print ( nt, t, '  TVEC_EVEN3' )

      return
      end
      subroutine test33 ( )

c*********************************************************************72
c
cc TEST33 tests TVEC_EVEN_BRACKET, TVEC_EVEN_BRACKET2, TVEC_EVEN_BRACKET3.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    22 July 2010
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer maxt
      parameter ( maxt = 5 )

      integer nt
      double precision t(maxt)
      double precision theta1
      double precision theta2

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST33'
      write ( *, '(a)' ) 
     &  '  For evenly spaced angles between THETA1 and THETA2:'
      write ( *, '(a)' ) '  TVEC_EVEN_BRACKET'
      write ( *, '(a)' ) '  TVEC_EVEN_BRACKET2.'
      write ( *, '(a)' ) '  TVEC_EVEN_BRACKET3.'

      nt = 4
      theta1 = 30.0D+00
      theta2 = 90.0D+00

      call tvec_even_bracket ( nt, theta1, theta2, t )

      write ( *, '(a)' ) ' '
      write ( *, '(a,i8)' ) '    NT = ', nt
      write ( *, '(a,g14.6)' ) '    THETA1 = ', theta1
      write ( *, '(a,g14.6)' ) '    THETA2 = ', theta2

      call r8vec_print ( nt, t, '  TVEC_EVEN_BRACKET' )

      nt = 5

      call tvec_even_bracket2 ( nt, theta1, theta2, t )

      write ( *, '(a)' ) ' '
      write ( *, '(a,i8)' ) '    NT = ', nt
      write ( *, '(a,g14.6)' ) '    THETA1 = ', theta1
      write ( *, '(a,g14.6)' ) '    THETA2 = ', theta2

      call r8vec_print ( nt, t, '  TVEC_EVEN_BRACKET2' )

      nt = 3

      call tvec_even_bracket3 ( nt, theta1, theta2, t )

      write ( *, '(a)' ) ' '
      write ( *, '(a,i8)' ) '    NT = ', nt
      write ( *, '(a,g14.6)' ) '    THETA1 = ', theta1
      write ( *, '(a,g14.6)' ) '    THETA2 = ', theta2

      call r8vec_print ( nt, t, '  TVEC_EVEN_BRACKET3' )

      return
      end
      subroutine test34 ( )

c*********************************************************************72
c
cc TEST34 tests UPC_CHECK_DIGIT.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    08 October 2005
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer test_num
      parameter ( test_num = 2 )

      integer c
      integer l
      integer l_test(test_num)
      integer p
      integer p_test(test_num)
      integer r
      integer r_test(test_num)
      integer test

      save l_test
      save p_test
      save r_test

      data l_test / 72890, 12345 /
      data p_test /     0,     0 /
      data r_test / 00011, 67890 /

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST34'
      write ( *, '(a)' ) 
     &  '  UPC_CHECK_DIGIT determines the check digit for a UPC.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  P-LLLLL-RRRRR-C'
      write ( *, '(a)' ) ' '

      do test = 1, test_num

        p = p_test(test)
        l = l_test(test)
        r = r_test(test)

        call upc_check_digit ( p, l, r, c )

        write ( *, '(2x,i1,''-'',i5.5,''-'',i5.5,''-'',i1 )' ) 
     &    p, l, r, c

      end do

      return
      end
      subroutine test35 ( )

c*********************************************************************72
c
cc TEST35 calls VERSINE_PULSE.
c
c  Modified:
c
c    21 July 2010
c
c  Author:
c
c    John Burkardt
c
      implicit none

      double precision amp
      integer i
      double precision t
      double precision ta
      double precision tb
      double precision v
      double precision versine_pulse
      double precision v1

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST35'
      write ( *, '(a)' ) '  VERSINE_PULSE adds a versine pulse '
      write ( *, '(a)' ) '  to a constant signal.'
      write ( *, '(a)' ) ' '

      ta = 2.0D+00
      tb = 4.0D+00
      v1 = 1.0D+00
      amp = 3.0D+00

      do i = 0, 100
        t = dble ( i ) / 10.0D+00
        v = versine_pulse ( t, ta, tb, v1, amp )
        write ( *, '(2x,i4,2x,f10.4,2x,f10.4)' ) i, t, v
      end do

      return
      end
