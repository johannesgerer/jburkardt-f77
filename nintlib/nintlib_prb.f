      program main

c*********************************************************************72
c
cc NINTLIB_PRB runs the NINTLIB tests.
c
c  Modified:
c
c    25 February 2007
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer test_num
      parameter ( test_num = 3 )

      double precision a
      double precision b
      integer dim_num
      integer dim_num_test(test_num)
      external f1dn
      double precision f1dn
      external fbdn
      double precision fbdn
      external fedn
      double precision fedn
      external fxdn
      double precision fxdn
      external fx2dn
      double precision fx2dn
      external fx3dn
      double precision fx3dn
      integer test

      save dim_num_test

      data dim_num_test / 2, 3, 4 /

      call timestamp ( )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'NINTLIB_PRB'
      write ( *, '(a)' ) '  FORTRAN77 version'
      write ( *, '(a)' ) '  Test the NINTLIB library.'

      a = 0.0D+00
      b = 1.0D+00

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TESTND'
      write ( *, '(a)' ) '  Test N-D quadrature codes'
      write ( *, '(a)' ) '  for integral of F(X) on [A,B]**NDIM.'
      write ( *, '(a)' ) ' '

      do test = 1, test_num

        dim_num = dim_num_test(test)

        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) ' '
        write ( *, '(a,i8)' ) '  NDIM = ', dim_num
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) ' '
        write ( *, '(a,g20.12)' ) '  A(1:NDIM) = ', a
        write ( *, '(a,g20.12)' ) '  B(1:NDIM) = ', b

        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) '  F(X)=1'
        write ( *, '(a)' ) ' '

        call testnd ( dim_num, f1dn )

        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) '  F(X)=X'
        write ( *, '(a)' ) ' '

        call testnd ( dim_num, fxdn )

        write ( *, '(a)') ' '
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) '  F(X)=X*2'
        write ( *, '(a)' ) ' '

        call testnd ( dim_num, fx2dn )

        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) '  F(X)=X*3'
        write ( *, '(a)' ) ' '

        call testnd ( dim_num, fx3dn )

        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) '  F(X)=EXP(X)'
        write ( *, '(a)' ) ' '

        call testnd ( dim_num, fedn )

        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) '  F(X)=1/(1+X*X)'
        write ( *, '(a)' ) ' '

        call testnd ( dim_num, fbdn )

      end do
c
c  Terminate.
c
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'NINTLIB_PRB'
      write ( *, '(a)' ) '  Normal end of execution.'

      write ( *, '(a)' ) ' '
      call timestamp ( )

      stop
      end
      subroutine testnd ( dim_num, func )

c****&****************************************************************72
c
cc TESTND tests the integrators on a particular function.
c
c  Modified:
c
c    25 February 2007
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer dim_num
      external func
      double precision func

      call test01 ( dim_num, func )
      call test02 ( dim_num, func )
      call test03 ( dim_num, func )
      call test04 ( dim_num, func )
      if ( dim_num .eq. 2 ) then
        call test05 ( dim_num, func )
      end if
      call test06 ( dim_num, func )

      return
      end
      subroutine test01 ( dim_num, func )

c****&****************************************************************72
c
cc TEST01 tests BOX_ND.
c
c  Modified:
c
c    25 February 2007
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer DIM_NUM, the spatial dimension.
c
c    Input, double precision, external FUNC, the name of the function
c    to be integrated.
c
      implicit none

      integer dim_num
      integer order
      parameter ( order = 5 )

      integer eval_num
      double precision func
      external func
      integer i
      double precision result
      double precision wtab(order)
      double precision wtab2(order)
      double precision xtab(order)
      double precision xtab2(order)

      save wtab
      save xtab

      data wtab /
     &  0.236926885056189087514264040720D+00,
     &  0.478628670499366468041291514836D+00,
     &  0.568888888888888888888888888889D+00,
     &  0.478628670499366468041291514836D+00,
     &  0.236926885056189087514264040720D+00 /
      data xtab /
     &  -0.906179845938663992797626878299D+00,
     &  -0.538469310105683091036314420700D+00,
     &   0.0D+00,
     &   0.538469310105683091036314420700D+00,
     &   0.906179845938663992797626878299D+00 /
c
c  Adjust the quadrature rule from [-1,1] to [0,1]:
c
      do i = 1, order
        xtab2(i) = ( xtab(i) + 1.0D+00 ) / 2.0D+00
      end do

      do i = 1, order
        wtab2(i) = 0.5D+00 * wtab(i)
      end do

      call box_nd ( func, dim_num, order, xtab2, wtab2, result,
     &  eval_num )

      write ( *, '(a,g20.12,2x,i8)' )
     &  '  BOX_ND:         ', result, eval_num

      return
      end
      subroutine test02 ( dim_num, func )

c****&****************************************************************72
c
cc TEST02 tests P5_ND.
c
c  Modified:
c
c    25 February 2007
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer DIM_NUM, the spatial dimension.
c
c    Input, double precision, external FUNC, the name of the function
c    to be integrated.
c
      implicit none

      integer dim_num

      double precision a(dim_num)
      double precision b(dim_num)
      integer dim
      integer eval_num
      double precision func
      external func
      double precision result
c
c  Set the integration limits.
c
      do dim = 1, dim_num
        a(dim) = 0.0D+00
        b(dim) = 1.0D+00
      end do

      call p5_nd ( func, dim_num, a, b, result, eval_num )

      write ( *, '(a,g20.12,2x,i8)' )
     &  '  P5_ND:          ', result, eval_num

      return
      end
      subroutine test03 ( dim_num, f )

c****&****************************************************************72
c
cc TEST03 tests ROMBERG_ND.
c
c  Modified:
c
c    11 September 2006
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer dim_num

      double precision a(dim_num)
      double precision al
      double precision b(dim_num)
      double precision e
      integer eval_num
      external f
      double precision f
      double precision h(dim_num)
      integer i
      integer ind
      integer key
      integer m
      double precision result

      al = 1.5D+00
      m = 3
c
c  Set the integration limits.
c
      do i = 1, dim_num
        a(i) = 0.0D+00
        b(i) = 1.0D+00
        h(i) = 1.0D+00 / 10.0D+00
      end do

      e = 0.001D+00

      call romberg_nd ( dim_num, a, b, h, al, m, e, f, result, key,
     &  eval_num )

      write ( *, '(a,g20.12,2x,i8)' )
     &  '  ROMBERG_ND:     ', result, eval_num

      return
      end
      subroutine test04 ( dim_num, func )

c****&****************************************************************72
c
cc TEST04 tests SAMPLE_ND.
c
c  Modified:
c
c    25 February 2007
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer DIM_NUM, the spatial dimension.
c
c    Input, double precision, external FUNC, the name of the function
c    to be integrated.
c
      implicit none

      integer dim_num

      integer k2
      parameter ( k2 = 4 )

      double precision dev1(k2)
      double precision dev2(k2)
      double precision err1(k2)
      double precision est1(k2)
      double precision est2(k2)
      double precision err2(k2)
      integer eval_num
      double precision func
      external func
      integer k1

      k1 = 1

      call sample_nd ( func, k1, k2, dim_num, est1, err1, dev1, est2,
     &  err2, dev2, eval_num )

      write ( *, '(a,g20.12,2x,i8)' )
     &  '  SAMPLE_ND:      ', est2(k2), eval_num

      return
      end
      subroutine test05 ( dim_num, func )

c****&****************************************************************72
c
cc TEST05 demonstrates how to refine N-dimensional integration results.
c
c  Discussion:
c
c    This routine is only set up for DIM_NUM = 2 for now.
c
c    We are given a routine, NDP5, which will integrate over a
c    DIM_NUM dimensional hypercube using a fixed method.  In order to
c    improve the approximation to an integral, we can subdivide
c    the hypercube and call NDP5 to integrate again over each of
c    these regions.
c
c    The information that we gather can be used to tell us when
c    to expect that we have achieved a certain degree of accuracy.
c
c    With a little more work, we could make this code adaptive.
c    That is, it would only refine SOME of the subregions, where
c    the approximation to the integral was still not good enough.
c
c  Modified:
c
c    25 February 2007
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer DIM_NUM, the spatial dimension.
c
c    Input, double precision, external FUNC, the name of the function
c    to be integrated.
c
      implicit none

      integer dim_num

      double precision a(dim_num)
      double precision b(dim_num)
      integer dim
      integer eval_num
      integer eval_total
      double precision func
      external func
      integer i
      integer igrid
      integer j
      integer ngrid
      double precision result
      double precision result_total
      double precision xlo(dim_num)
      double precision xhi(dim_num)

      do dim = 1, dim_num
        xlo(dim) = 0.0D+00
        xhi(dim) = 1.0D+00
      end do

      do igrid = 1, 6

        ngrid = 2**( igrid - 1 )

        result_total = 0.0D+00
        eval_total = 0

        do i = 1, ngrid

          a(1) = ( dble ( ngrid - i + 1 ) * xlo(1)
     &           + dble (         i - 1 ) * xhi(1) )
     &           / dble ( ngrid         )

          b(1) = ( dble ( ngrid - i ) * xlo(1)
     &           + dble (         i ) * xhi(1) )
     &           / dble ( ngrid     )

          do j = 1, ngrid

            a(2) = ( dble ( ngrid - j + 1 ) * xlo(2)
     &             + dble (         j - 1 ) * xhi(2) )
     &             / dble ( ngrid         )

            b(2) = ( dble ( ngrid - j ) * xlo(2)
     &             + dble (         j ) * xhi(2) )
     &             / dble ( ngrid     )

            call p5_nd ( func, dim_num, a, b, result, eval_num )

            result_total = result_total + result
            eval_total = eval_total + eval_num

          end do

        end do

        write ( *, '(a,g20.12,2x,i8)' )
     &    '  P5_ND+:         ', result_total, eval_total

      end do

      return
      end
      subroutine test06 ( dim_num, func )

c****&****************************************************************72
c
cc TEST06 tests MONTE_CARLO_ND.
c
c  Modified:
c
c    25 February 2007
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer DIM_NUM, the spatial dimension.
c
c    Input, double precision, external FUNC, the name of the function
c    to be integrated.
c
      implicit none

      integer dim_num
      integer test_num
      parameter ( test_num = 3 )

      double precision a(dim_num)
      double precision b(dim_num)
      integer dim
      integer eval_num
      double precision func
      external func
      double precision result
      integer seed
      integer test

      seed = 123456789
c
c  Set the integration limits.
c
      do dim = 1, dim_num
        a(dim) = 0.0D+00
        b(dim) = 1.0D+00
      end do

      do test = 1, test_num

        eval_num = 8**test * 10000

        call monte_carlo_nd ( func, dim_num, a, b, eval_num, seed,
     &    result )

        write ( *, '(a,g20.12,2x,i8)' )
     &    '  MONTE_CARLO_ND: ', result, eval_num

      end do

      return
      end
      function fbdn ( dim_num, x )

c****&****************************************************************72
c
cc FBDN(X(1:DIM_NUM)) = 1 / ( 1 + SUM ( X(1:DIM_NUM)**2 ) )
c
c  Modified:
c
c    25 February 2007
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer DIM_NUM, the spatial dimension.
c
c    Input, double precision X(DIM_NUM), the argument.
c
c    Output, double precision FBDN, the value of the function at X.
c
      implicit none

      integer dim_num

      double precision fbdn
      integer i
      double precision temp
      double precision x(dim_num)

      temp = 0.0D+00
      do i = 1, dim_num
        temp = temp + x(i) * x(i)
      end do

      fbdn = 1.0D+00 / ( 1.0D+00 + temp )

      return
      end
      function f1dn ( dim_num, x )

c****&****************************************************************72
c
cc F1DN(X(1:DIM_NUM)) = 1.
c
c  Modified:
c
c    25 February 2007
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer DIM_NUM, the spatial dimension.
c
c    Input, double precision X(DIM_NUM), the argument.
c
c    Output, double precision F1DN, the value of the function at X.
c
      implicit none

      integer dim_num

      double precision f1dn
      double precision x(dim_num)

      f1dn = 1.0D+00

      return
      end
      function fxdn ( dim_num, x )

c****&****************************************************************72
c
cc FXDN(X(1:DIM_NUM)) = SUM ( X(1:DIM_NUM) )
c
c  Modified:
c
c    25 February 2007
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer DIM_NUM, the spatial dimension.
c
c    Input, double precision X(DIM_NUM), the argument.
c
c    Output, double precision FXDN, the value of the function at X.
c
      implicit none

      integer dim_num

      double precision fxdn
      integer i
      double precision temp
      double precision x(dim_num)

      temp = 0.0D+00
      do i = 1, dim_num
        temp = temp + x(i)
      end do

      fxdn = temp

      return
      end
      function fx2dn ( dim_num, x )

c****&****************************************************************72
c
cc FX2DN(X(1:DIM_NUM)) = SUM ( X(1:DIM_NUM)**2 )
c
c  Modified:
c
c    25 February 2007
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer DIM_NUM, the spatial dimension.
c
c    Input, double precision X(DIM_NUM), the argument.
c
c    Output, double precision FX2DN, the value of the function at X.
c
      implicit none

      integer dim_num

      double precision fx2dn
      integer i
      double precision temp
      double precision x(dim_num)

      temp = 0.0D+00
      do i = 1, dim_num
        temp = temp + x(i) * x(i)
      end do

      fx2dn = temp

      return
      end
      function fx3dn ( dim_num, x )

c****&****************************************************************72
c
cc FX3DN(X(1:DIM_NUM)) = SUM ( X(1:DIM_NUM)**3 )
c
c  Modified:
c
c    25 February 2007
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer DIM_NUM, the spatial dimension.
c
c    Input, double precision X(DIM_NUM), the argument.
c
c    Output, double precision FX3DN, the value of the function at X.
c
      implicit none

      integer dim_num

      double precision fx3dn
      integer i
      double precision temp
      double precision x(dim_num)

      temp = 0.0D+00
      do i = 1, dim_num
        temp = temp + x(i) * x(i) * x(i)
      end do

      fx3dn = temp

      return
      end
      function fedn ( dim_num, x )

c****&****************************************************************72
c
cc FEDN(X(1:DIM_NUM)) = EXP ( SUM ( X(1:DIM_NUM) ) )
c
c  Modified:
c
c    25 February 2007
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer DIM_NUM, the spatial dimension.
c
c    Input, double precision X(DIM_NUM), the argument.
c
c    Output, double precision FEDN, the value of the function at X.
c
      implicit none

      integer dim_num

      double precision fedn
      integer i
      double precision temp
      double precision x(dim_num)

      temp = 0.0D+00
      do i = 1, dim_num
        temp = temp + x(i) * x(i)
      end do

      fedn = exp ( temp )

      return
      end
