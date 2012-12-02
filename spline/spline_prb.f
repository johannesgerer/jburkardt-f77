      program main

c*********************************************************************72
c
cc MAIN is the main program for SPLINE_PRB.
c
c  Discussion:
c
c    SPLINE_PRB calls the SPLINE tests.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    06 February 2012
c
c  Author:
c
c    John Burkardt
c
      implicit none

      call timestamp ( )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'SPLINE_PRB'
      write ( *, '(a)' ) '  FORTRAN77 version:'
      write ( *, '(a)' ) '  Test the SPLINE library.'

      call test001 ( )
      call test002 ( )
      call test003 ( )
      call test004 ( )
      call test005 ( )
      call test006 ( )

      call test01 ( )
      call test02 ( )
      call test03 ( )
      call test04 ( )
      call test05 ( )
      call test06 ( )
      call test07 ( )
      call test08 ( )
      call test09 ( )
      call test10 ( )

      call test11 ( )
      call test115 ( )
      call test116 ( )
      call test12 ( )
      call test125 ( )
      call test126 ( )
      call test127 ( )
      call test13 ( )
      call test14 ( )
      call test143 ( )
      call test144 ( )
      call test145 ( )
      call test15 ( )
      call test16 ( )
      call test17 ( )
      call test18 ( )
      call test19 ( )
      call test20 ( )
      call test205 ( )

      call test21 ( )
      call test215 ( )
      call test22 ( )
      call test225 ( )
      call test23 ( )
      call test235 ( )
      call test24 ( )
c
c  Terminate.
c
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'SPLINE_PRB'
      write ( *, '(a)' ) '  Normal end of execution.'

      write ( *, '(a)' ) ' '
      call timestamp ( )

      stop
      end
      subroutine test001 ( )

c*********************************************************************72
c
cc TEST001 tests PARABOLA_VAL2.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    17 June 2006
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer dim_num
      parameter ( dim_num = 1 )
      integer ndata
      parameter ( ndata = 5 )

      integer i
      integer left
      double precision xdata(ndata)
      double precision xval
      double precision ydata(dim_num,ndata)
      double precision yval(dim_num)
      double precision zdata(ndata)
      double precision zval(dim_num)

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST001'
      write ( *, '(a)' ) '  PARABOLA_VAL2 evaluates parabolas through'
      write ( *, '(a)' ) '    3 points in a table'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Our data tables will actually be parabolas:'
      write ( *, '(a)' ) '    Y: 2*x**2 + 3 * x + 1.'
      write ( *, '(a)' ) '    Z: 4*x**2 - 2 * x + 5.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '         I        X         Y           Z'
      write ( *, '(a)' ) ' '

      do i = 1, ndata

        xval = 2.0D+00 * dble ( i )
        xdata(i) = xval
        ydata(1,i) = 2.0D+00 * xval * xval + 3.0 * xval + 1.0D+00
        zdata(i) = 4.0D+00 * xval * xval - 2.0D+00 * xval + 5.0D+00
        write ( *, '(2x,i8,2x,f10.4,2x,f10.4,2x,f10.4)' ) 
     &    i, xdata(i), ydata(1,i), zdata(i)

      end do

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Interpolated data:'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '      LEFT        X         Y           Z'
      write ( *, '(a)' ) ' '

      do i = 1, 5

        xval = dble ( 2 * i - 1 )
        left = i

        if ( ndata - 2 .lt. left ) then
          left = ndata - 2
        end if

        if ( left .lt. 1 ) then
          left = 1
        end if

        call parabola_val2 ( dim_num, ndata, xdata, ydata, left, 
     &    xval, yval )

        call parabola_val2 ( dim_num, ndata, xdata, zdata, left, 
     &    xval, zval )

        write ( *, '(2x,i8,2x,f10.4,2x,f10.4,2x,f10.4)' ) 
     &    left, xval, yval(1), zval(1)

      end do

      return
      end
      subroutine test002 ( )

c*********************************************************************72
c
cc TEST002 tests R8VEC_BRACKET.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    17 June 2006
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n
      parameter ( n = 10 )
      integer test_num
      parameter ( test_num = 6 )

      integer left
      integer right
      integer test
      double precision x(n)
      double precision xtest(test_num)
      double precision xval

      save xtest

      data xtest /
     & -10.0D+00, 1.0D+00, 4.5D+00, 5.0D+00, 10.0D+00, 12.0D+00 /

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST002'
      write ( *, '(a)' ) '  R8VEC_BRACKET finds a pair of entries in a'
      write ( *, '(a)' ) '    sorted real array which bracket a value.'

      call r8vec_indicator ( n, x )
      x(6) = x(5)

      call r8vec_print ( n, x, '  Sorted array:' )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '    LEFT             RIGHT'
      write ( *, '(a)' ) '  X(LEFT)   XVAL   X(RIGHT)'
      write ( *, '(a)' ) ' '

      do test = 1, test_num

        xval = xtest(test)

        call r8vec_bracket ( n, x, xval, left, right )

        write ( *, '(2x,i14,14x,i14)' ) left, right

        if ( 1 .le. left .and. 1 .le. right ) then
          write ( *, '(2x,3g14.6)' ) x(left), xval, x(right)
        else if ( left .lt. 1 .and. 1 .le. right ) then
          write ( *, '(2x,14x,2g14.6)' )          xval, x(right)
        else if ( 1 .le. left .and. right .lt. 1 ) then
          write ( *, '(2x,2g14.6)' ) x(left), xval
        else if ( left .lt. 1 .and. right .lt. 1 ) then
          write ( *, '(2x,14x,g14.6)' )          xval
        end if

      end do

      return
      end
      subroutine test003 ( )

c*********************************************************************72
c
cc TEST003 tests R8VEC_BRACKET3.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    17 June 2006
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n
      parameter ( n = 10 )
      integer test_num
      parameter ( test_num = 6 )

      integer left
      integer test
      double precision x(n)
      double precision xtest(test_num)
      double precision xval

      save xtest

      data xtest /
     & -10.0D+00, 1.0D+00, 4.5D+00, 5.0D+00, 10.0D+00, 12.0D+00 /

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST003'
      write ( *, '(a)' ) '  R8VEC_BRACKET3 finds a pair of entries in a'
      write ( *, '(a)' ) '    sorted real array which bracket a value.'

      call r8vec_indicator ( n, x )
      x(6) = x(5)

      call r8vec_print ( n, x, '  Sorted array:' )

      left = ( n + 1 ) / 2

      do test = 1, test_num

        xval = xtest(test)

        write ( *, '(a)' ) ' '
        write ( *, '(a,g14.6)' ) '  Search for XVAL = ', xval

        write ( *, '(a,i8)' ) 
     &    '  Starting guess for interval is = ', left

        call r8vec_bracket3 ( n, x, xval, left )

        write ( *, '(a)' ) '  Nearest interval:'
        write ( *, '(2x,a,i8,a,g14.6)' ) 
     &    '    X[', left,' ]= ', x(left)
        write ( *, '(2x,a,i8,a,g14.6)' ) 
     &    '    X[', left+1, ' ]= ', x(left+1)

      end do

      return
      end
      subroutine test004 ( )

c*********************************************************************72
c
cc TEST004 tests R8VEC_ORDER_TYPE.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    17 June 2006
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n
      parameter ( n = 4 )
      integer test_num
      parameter ( test_num = 6 )

      integer j
      integer order
      integer test
      double precision x(n,test_num)

      x(1,1) = 1.0D+00
      x(2,1) = 3.0D+00
      x(3,1) = 2.0D+00
      x(4,1) = 4.0D+00

      x(1,2) = 2.0D+00
      x(2,2) = 2.0D+00
      x(3,2) = 2.0D+00
      x(4,2) = 2.0D+00

      x(1,3) = 1.0D+00
      x(2,3) = 2.0D+00
      x(3,3) = 2.0D+00
      x(4,3) = 4.0D+00

      x(1,4) = 1.0D+00
      x(2,4) = 2.0D+00
      x(3,4) = 3.0D+00
      x(4,4) = 4.0D+00

      x(1,5) = 4.0D+00
      x(2,5) = 4.0D+00
      x(3,5) = 3.0D+00
      x(4,5) = 1.0D+00

      x(1,6) = 9.0D+00
      x(2,6) = 7.0D+00
      x(3,6) = 3.0D+00
      x(4,6) = 0.0D+00

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST004'
      write ( *, '(a)' ) 
     &  '  R8VEC_ORDER_TYPE classifies a real vector as'
      write ( *, '(a)' ) '  -1: no order'
      write ( *, '(a)' ) '   0: all equal;'
      write ( *, '(a)' ) '   1: ascending;'
      write ( *, '(a)' ) '   2: strictly ascending;'
      write ( *, '(a)' ) '   3: descending;'
      write ( *, '(a)' ) '   4: strictly descending.'
      write ( *, '(a)' ) ' '

      do test = 1, test_num

        call r8vec_order_type ( n, x(1,test), order )

        write ( *, '(a)' ) ' '
        write ( *, '(a,i8)' ) 
     &    '  The following vector has order type ', order
        write ( *, '(a)' ) ' '
        do j = 1, n
          write ( *, '(2x,i8,g14.6)' ) j, x(j,test)
        end do

      end do

      return
      end
      subroutine test005 ( )

c*********************************************************************72
c
cc TEST005 tests R83_NP_FS.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    17 June 2006
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n
      parameter ( n = 10 )

      double precision a(3,n)
      double precision b(n)
      integer seed
      double precision x(n)

      seed = 123456789

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST005'
      write ( *, '(a)' ) '  R83_NP_FS factors and solves a tridiagonal'
      write ( *, '(a)' ) '  linear system.'
      write ( *, '(a)' ) ' '
      write ( *, '(a,i8)' ) '  Matrix order N = ', n
c
c  Set the matrix elements.
c
      call r83_uniform ( n, seed, a )
c
c  Set the desired solution.
c
      call r8vec_indicator ( n, x )
c
c  Compute b = A * x.
c
      call r83_mxv ( n, a, x, b )
c
c  Wipe out the solution.
c
      x(1:n) = 0.0D+00
c
c  Solve the system.
c
      call r83_np_fs ( n, a, b, x )

      call r8vec_print ( n, x, '  Solution:' )

      return
      end
      subroutine test006 ( )

c*********************************************************************72
c
cc TEST006 tests DATA_TO_DIF and DIF_VAL.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    17 June 2006
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer maxtab
      parameter ( maxtab = 8 )

      double precision diftab(maxtab)
      double precision err
      double precision exact
      integer j
      integer ntab
      double precision xtab(maxtab)
      double precision xval
      double precision ytab(maxtab)
      double precision yval

      xval = 2.5D+00
      exact = exp ( xval )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST006'
      write ( *, '(a,i8)' ) 
     &  '  Approximate Y = EXP(X) using orders 1 to ', maxtab
      write ( *, '(a,g14.6)' ) '  Evaluate at X = ', xval
      write ( *, '(a,g14.6)' ) '  where EXP(X)=   ', exact

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '    Order  Approximate Y     Error'
      write ( *, '(a)' ) ' '

      do ntab = 1, maxtab

        do j = 1, ntab
          xtab(j) = dble ( j - 1 )
          ytab(j) = exp ( xtab(j) )
        end do

        call data_to_dif ( ntab, xtab, ytab, diftab )

        call dif_val ( ntab, xtab, diftab, xval, yval )

        err = yval - exact
        write ( *, ' ( 2x, i8, 2g14.6 )' ) ntab, yval, err

      end do

      return
      end
      subroutine test01 ( )

c*********************************************************************72
c
cc TEST01 tests BASIS_FUNCTION_B_VAL.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    17 June 2006
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer ndata
      parameter ( ndata = 5 )

      integer i
      integer j
      integer jhi
      character mark
      integer nsample
      parameter ( nsample = 4 )
      double precision tdata(ndata)
      double precision thi
      double precision tlo
      double precision tval
      double precision yval

      save tdata

      data tdata /
     &  0.0D+00, 1.0D+00, 4.0D+00, 6.0D+00, 10.0D+00 /

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST01'
      write ( *, '(a)' ) '  BASIS_FUNCTION_B_VAL evaluates the '
      write ( *, '(a)' ) '    B spline basis function.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '    T            B(T)'
      write ( *, '(a)' ) ' '

      do i = 0, ndata

        if ( i .eq. 0 ) then
          tlo = tdata(1) - 0.5D+00 * ( tdata(2) - tdata(1) )
          thi = tdata(1)
        else if ( i .lt. ndata ) then
          tlo = tdata(i)
          thi = tdata(i+1)
        else if ( ndata .le. i ) then
          tlo = tdata(ndata)
          thi = tdata(ndata) 
     &      + 0.5D+00 * ( tdata(ndata) - tdata(ndata-1) )
        end if

        if ( i .lt. ndata ) then
          jhi = nsample - 1
        else
          jhi = nsample
        end if

        do j = 0, jhi

          tval = ( dble ( nsample - j ) * tlo   
     &           + dble (           j ) * thi ) 
     &           / dble ( nsample     )

          call basis_function_b_val ( tdata, tval, yval )

          if ( 0 .lt. i .and. j .eq. 0 ) then
            mark = '*'
          else
            mark = ' '
          end if

          write ( *, '(2x,a1,2g14.6)' ) mark, tval, yval

        end do

      end do

      return
      end
      subroutine test02 ( )

c*********************************************************************72
c
cc TEST02 tests BASIS_FUNCTION_BETA_VAL.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    17 June 2006
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer ndata
      parameter ( ndata = 5 )

      double precision beta1
      double precision beta2
      integer i
      integer j
      integer jhi
      character mark
      integer nsample
      parameter ( nsample = 4 )
      double precision tdata(ndata)
      double precision thi
      double precision tlo
      double precision tval
      double precision yval

      save tdata

      data tdata /
     & 0.0D+00, 1.0D+00, 4.0D+00, 6.0D+00, 10.0D+00 /

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST02'
      write ( *, '(a)' ) '  BASIS_FUNCTION_BETA_VAL evaluates the '
      write ( *, '(a)' ) '    Beta spline basis function.'

      beta1 = 1.0D+00
      beta2 = 0.0D+00

      write ( *, '(a)' ) ' '
      write ( *, '(a,g14.6)' ) '  BETA1 = ', beta1
      write ( *, '(a,g14.6)' ) '  BETA2 = ', beta2
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '       T            B(T)'
      write ( *, '(a)' ) ' '

      do i = 0, ndata

        if ( i .eq. 0 ) then
          tlo = tdata(1) - 0.5D+00 * ( tdata(2) - tdata(1) )
          thi = tdata(1)
        else if ( i .lt. ndata ) then
          tlo = tdata(i)
          thi = tdata(i+1)
        else if ( ndata .le. i ) then
          tlo = tdata(ndata)
          thi = tdata(ndata) 
     &      + 0.5D+00 * ( tdata(ndata) - tdata(ndata-1) )
        end if

        if ( i .lt. ndata ) then
          jhi = nsample - 1
        else
          jhi = nsample
        end if

        do j = 0, jhi

          tval = ( dble ( nsample - j ) * tlo   
     &           + dble (           j ) * thi ) 
     &           / dble ( nsample     )

          call basis_function_beta_val ( beta1, beta2, tdata, tval,
     &      yval )

          if ( 0 .lt. i .and. j .eq. 0 ) then
            mark = '*'
          else
            mark = ' '
          end if

          write ( *, '(2x,a1,2g14.6)' ) mark, tval, yval

        end do

      end do

      beta1 = 1.0D+00
      beta2 = 100.0D+00

      write ( *, '(a)' ) ' '
      write ( *, '(a,g14.6)' ) '  BETA1 = ', beta1
      write ( *, '(a,g14.6)' ) '  BETA2 = ', beta2
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '       T            B(T)'
      write ( *, '(a)' ) ' '

      do i = 0, ndata

        if ( i .eq. 0 ) then
          tlo = tdata(1) - 0.5D+00 * ( tdata(2) - tdata(1) )
          thi = tdata(1)
        else if ( i .lt. ndata ) then
          tlo = tdata(i)
          thi = tdata(i+1)
        else if ( ndata .le. i ) then
          tlo = tdata(ndata)
          thi = tdata(ndata) 
     &      + 0.5D+00 * ( tdata(ndata) - tdata(ndata-1) )
        end if

        if ( i .lt. ndata ) then
          jhi = nsample - 1
        else
          jhi = nsample
        end if

        do j = 0, jhi

          tval = ( dble ( nsample - j ) * tlo   
     &           + dble (           j ) * thi ) 
     &           / dble ( nsample     )

          call basis_function_beta_val ( beta1, beta2, tdata, tval, 
     &      yval )

          if ( 0 .lt. i .and. j .eq. 0 ) then
            mark = '*'
          else
            mark = ' '
          end if

          write ( *, '(2x,a1,2g14.6)' ) mark, tval, yval

        end do

      end do

      beta1 = 100.0D+00
      beta2 = 0.0D+00

      write ( *, '(a)' ) ' '
      write ( *, '(a,g14.6)' ) '  BETA1 = ', beta1
      write ( *, '(a,g14.6)' ) '  BETA2 = ', beta2
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '       T            B(T)'
      write ( *, '(a)' ) ' '

      do i = 0, ndata

        if ( i .eq. 0 ) then
          tlo = tdata(1) - 0.5D+00 * ( tdata(2) - tdata(1) )
          thi = tdata(1)
        else if ( i .lt. ndata ) then
          tlo = tdata(i)
          thi = tdata(i+1)
        else if ( ndata .le. i ) then
          tlo = tdata(ndata)
          thi = tdata(ndata) 
     &      + 0.5D+00 * ( tdata(ndata) - tdata(ndata-1) )
        end if

        if ( i .lt. ndata ) then
          jhi = nsample - 1
        else
          jhi = nsample
        end if

        do j = 0, jhi

          tval = ( dble ( nsample - j ) * tlo   
     &           + dble (           j ) * thi ) 
     &           / dble ( nsample     )

          call basis_function_beta_val ( beta1, beta2, tdata, tval, 
     &      yval )

          if ( 0 .lt. i .and. j .eq. 0 ) then
            mark = '*'
          else
            mark = ' '
          end if

          write ( *, '(2x,a1,2g14.6)' ) mark, tval, yval

        end do

      end do

      return
      end
      subroutine test03 ( )

c*********************************************************************72
c
cc TEST03 tests BASIS_MATRIX_B_UNI and BASIS_MATRIX_TMP.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    17 June 2006
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n
      parameter ( n = 4 )
      integer ndata
      parameter ( ndata = 4 )

      integer i
      integer j
      integer jhi
      integer left
      character mark
      double precision mbasis(n,n)
      integer nsample
      parameter ( nsample = 4 )
      double precision tdata(ndata)
      double precision thi
      double precision tlo
      double precision tval
      double precision ydata(ndata)
      double precision yval

      save tdata
      save ydata

      data tdata /
     &  -1.0D+00, 0.0D+00, 1.0D+00, 2.0D+00 /
      data ydata /
     &  4.0D+00, 7.0D+00, 12.0D+00, 19.0D+00 /

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST03'
      write ( *, '(a)' ) '  BASIS_MATRIX_B_UNI sets up the basis matrix'
      write ( *, '(a)' ) '    for the uniform B spline.'

      call basis_matrix_b_uni ( mbasis )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '       TDATA         YDATA'
      write ( *, '(a)' ) ' '
      do i = 1, ndata
        write ( *, '(2x,2g14.6)' ) tdata(i), ydata(i)
      end do

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '        T            Spline(T)'
      write ( *, '(a)' ) ' '

      left = 2

      do i = 0, ndata

        if ( i .eq. 0 ) then
          tlo = tdata(1) - 0.5D+00 * ( tdata(2) - tdata(1) )
          thi = tdata(1)
        else if ( i .lt. ndata ) then
          tlo = tdata(i)
          thi = tdata(i+1)
        else if ( ndata .le. i ) then
          tlo = tdata(ndata)
          thi = tdata(ndata) 
     &      + 0.5D+00 * ( tdata(ndata) - tdata(ndata-1) )
        end if

        if ( i .lt. ndata ) then
          jhi = nsample - 1
        else
          jhi = nsample
        end if

        do j = 0, jhi

          tval = ( dble ( nsample - j ) * tlo   
     &           + dble (           j ) * thi ) 
     &           / dble ( nsample     )

          call basis_matrix_tmp ( left, n, mbasis, ndata, tdata, 
     &      ydata, tval, yval )

          if ( 0 .lt. i .and. j .eq. 0 ) then
            mark = '*'
          else
            mark = ' '
          end if

          write ( *, '(2x,a1,2g14.6)' ) mark, tval, yval

        end do

      end do

      return
      end
      subroutine test04 ( )

c*********************************************************************72
c
cc TEST04 tests BASIS_MATRIX_BETA_UNI and BASIS_MATRIX_TMP.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    17 June 2006
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n
      parameter ( n = 4 )
      integer ndata
      parameter ( ndata = 4 )

      double precision beta1
      double precision beta2
      integer i
      integer j
      integer jhi
      integer left
      character mark
      double precision mbasis(n,n)
      integer nsample
      parameter ( nsample = 4 )
      double precision tdata(ndata)
      double precision thi
      double precision tlo
      double precision tval
      double precision ydata(ndata)
      double precision yval

      save tdata
      save ydata

      data tdata /
     &  -1.0D+00, 0.0D+00, 1.0D+00, 2.0D+00 /
      data ydata /
     &  4.0D+00, 7.0D+00, 12.0D+00, 19.0D+00 /

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST04'
      write ( *, '(a)' ) 
     &  '  BASIS_MATRIX_BETA_UNI sets up the basis matrix'
      write ( *, '(a)' ) '  for the uniform beta spline.'
c
c  First test
c
      beta1 = 1.0D+00
      beta2 = 0.0D+00

      write ( *, '(a)' ) ' '
      write ( *, '(a,g14.6)' ) '  BETA1 = ', beta1
      write ( *, '(a,g14.6)' ) '  BETA2 = ', beta2

      call basis_matrix_beta_uni ( beta1, beta2, mbasis )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '    TDATA, YDATA'
      write ( *, '(a)' ) ' '
      do i = 1, ndata
        write ( *, '(2g14.6)' ) tdata(i), ydata(i)
      end do

      left = 2

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '    T, Spline(T)'
      write ( *, '(a)' ) ' '

      do i = 0, ndata

        if ( i .eq. 0 ) then
          tlo = tdata(1) - 0.5D+00 * ( tdata(2) - tdata(1) )
          thi = tdata(1)
        else if ( i .lt. ndata ) then
          tlo = tdata(i)
          thi = tdata(i+1)
        else if ( ndata .le. i ) then
          tlo = tdata(ndata)
          thi = tdata(ndata) 
     &      + 0.5D+00 * ( tdata(ndata) - tdata(ndata-1) )
        end if

        if ( i .lt. ndata ) then
          jhi = nsample - 1
        else
          jhi = nsample
        end if

        do j = 0, jhi

          tval = ( dble ( nsample - j ) * tlo   
     &           + dble (           j ) * thi ) 
     &           / dble ( nsample     )

          call basis_matrix_tmp ( left, n, mbasis, ndata, tdata, ydata, 
     &      tval, yval )

          if ( 0 .lt. i .and. j .eq. 0 ) then
            mark = '*'
          else
            mark = ' '
          end if

          write ( *, '(2x,a1,2g14.6)' ) mark, tval, yval

        end do

      end do
c
c  Second test
c
      beta1 = 1.0D+00
      beta2 = 100.0D+00

      write ( *, '(a)' ) ' '
      write ( *, '(a,g14.6)' ) '  BETA1 = ', beta1
      write ( *, '(a,g14.6)' ) '  BETA2 = ', beta2

      call basis_matrix_beta_uni ( beta1, beta2, mbasis )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '    TDATA, YDATA'
      write ( *, '(a)' ) ' '
      do i = 1, ndata
        write ( *, '(2x,2g14.6)' ) tdata(i), ydata(i)
      end do

      left = 2

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '    T, Spline(T)'
      write ( *, '(a)' ) ' '

      do i = 0, ndata

        if ( i .eq. 0 ) then
          tlo = tdata(1) - 0.5D+00 * ( tdata(2) - tdata(1) )
          thi = tdata(1)
        else if ( i .lt. ndata ) then
          tlo = tdata(i)
          thi = tdata(i+1)
        else if ( ndata .le. i ) then
          tlo = tdata(ndata)
          thi = tdata(ndata) 
     &      + 0.5D+00 * ( tdata(ndata) - tdata(ndata-1) )
        end if

        if ( i .lt. ndata ) then
          jhi = nsample - 1
        else
          jhi = nsample
        end if

        do j = 0, jhi

          tval = ( dble ( nsample - j ) * tlo   
     &           + dble (           j ) * thi ) 
     &           / dble ( nsample     )

          call basis_matrix_tmp ( left, n, mbasis, ndata, tdata, ydata, 
     &      tval, yval )

          if ( 0 .lt. i .and. j .eq. 0 ) then
            mark = '*'
          else
            mark = ' '
          end if

          write ( *, '(2x,a1,2g14.6)' ) mark, tval, yval

        end do

      end do
c
c  Third test
c
      beta1 = 100.0D+00
      beta2 = 0.0D+00

      write ( *, '(a)' ) ' '
      write ( *, '(a,g14.6)' ) '  BETA1 = ', beta1
      write ( *, '(a,g14.6)' ) '  BETA2 = ', beta2

      call basis_matrix_beta_uni ( beta1, beta2, mbasis )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '     TDATA        YDATA'
      write ( *, '(a)' ) ' '
      do i = 1, ndata
        write ( *, '(2x,2g14.6)' ) tdata(i), ydata(i)
      end do

      left = 2

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '        T           Spline(T)'
      write ( *, '(a)' ) ' '

      do i = 0, ndata

        if ( i .eq. 0 ) then
          tlo = tdata(1) - 0.5D+00 * ( tdata(2) - tdata(1) )
          thi = tdata(1)
        else if ( i .lt. ndata ) then
          tlo = tdata(i)
          thi = tdata(i+1)
        else if ( ndata .le. i ) then
          tlo = tdata(ndata)
          thi = tdata(ndata) 
     &      + 0.5D+00 * ( tdata(ndata) - tdata(ndata-1) )
        end if

        if ( i .lt. ndata ) then
          jhi = nsample - 1
        else
          jhi = nsample
        end if

        do j = 0, jhi

          tval = ( dble ( nsample - j ) * tlo   
     &           + dble (           j ) * thi ) 
     &           / dble ( nsample     )

          call basis_matrix_tmp ( left, n, mbasis, ndata, tdata, ydata, 
     &      tval, yval )

          if ( 0 .lt. i .and. j .eq. 0 ) then
            mark = '*'
          else
            mark = ' '
          end if

          write ( *, '(2x,a1,2g14.6)' ) mark, tval, yval

        end do

      end do

      return
      end
      subroutine test05 ( )

c*********************************************************************72
c
cc TEST05 tests BASIS_MATRIX_BEZIER and BASIS_MATRIX_TMP.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    17 June 2006
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n
      parameter ( n = 4 )
      integer ndata
      parameter ( ndata = 4 )

      integer i
      integer j
      integer jhi
      integer left
      character mark
      double precision mbasis(n,n)
      integer nsample
      parameter ( nsample = 4 )
      double precision tdata(ndata)
      double precision thi
      double precision tlo
      double precision tval
      double precision ydata(ndata)
      double precision yval

      save tdata
      save ydata

      data tdata /
     &  0.0D+00, 0.0D+00, 1.0D+00, 1.0D+00 /
      data ydata /
     &  7.0D+00,  8.3333333D+00,   10.0D+00, 12.0D+00 /

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST05'
      write ( *, '(a)' ) '  BASIS_MATRIX_BEZIER sets up the basis'
      write ( *, '(a)' ) '    matrix for the uniform Bezier spline.'

      call basis_matrix_bezier ( mbasis )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '       TDATA          YDATA'
      write ( *, '(a)' ) ' '
      do i = 1, ndata
        write ( *, '(2x,2g14.6)' ) tdata(i), ydata(i)
      end do

      left = 2

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '        T            Spline(T)'
      write ( *, '(a)' ) ' '

      do i = 0, ndata

        if ( i .eq. 0 ) then
          tlo = tdata(1) - 0.5D+00 * ( tdata(2) - tdata(1) )
          thi = tdata(1)
        else if ( i .lt. ndata ) then
          tlo = tdata(i)
          thi = tdata(i+1)
        else if ( ndata .le. i ) then
          tlo = tdata(ndata)
          thi = tdata(ndata) 
     &      + 0.5D+00 * ( tdata(ndata) - tdata(ndata-1) )
        end if

        if ( i .lt. ndata ) then
          jhi = nsample - 1
        else
          jhi = nsample
        end if

        do j = 0, jhi

          tval = ( dble ( nsample - j ) * tlo   
     &           + dble (           j ) * thi ) 
     &           / dble ( nsample     )

          call basis_matrix_tmp ( left, n, mbasis, ndata, tdata, ydata, 
     &      tval, yval )

          if ( 0 .lt. i .and. j .eq. 0 ) then
            mark = '*'
          else
            mark = ' '
          end if

          write ( *, '(2x,a1,2g14.6)' ) mark, tval, yval

        end do

      end do

      return
      end
      subroutine test06 ( )

c*********************************************************************72
c
cc TEST06 tests BASIS_MATRIX_HERMITE and BASIS_MATRIX_TMP.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    17 June 2006
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n
      parameter ( n = 4 )
      integer ndata
      parameter ( ndata = 4 )

      integer i
      integer j
      integer jhi
      integer left
      character mark
      double precision mbasis(n,n)
      integer nsample
      parameter ( nsample = 4 )
      double precision tdata(ndata)
      double precision thi
      double precision tlo
      double precision tval
      double precision ydata(ndata)
      double precision yval

      save tdata
      save ydata

      data tdata /
     &  0.0D+00, 0.0D+00, 1.0D+00, 1.0D+00 /
      data ydata /
     &  7.0D+00, 12.0D+00, 4.0D+00, 6.0D+00 /

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST06'
      write ( *, '(a)' ) 
     &  '  BASIS_MATRIX_HERMITE sets up the basis matrix'
      write ( *, '(a)' ) '    for the Hermite spline.'

      call basis_matrix_hermite ( mbasis )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '       TDATA        YDATA'
      write ( *, '(a)' ) ' '
      do i = 1, ndata
        write ( *, '(2x,2g14.6)' ) tdata(i), ydata(i)
      end do

      left = 2

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '        T           Spline(T)'
      write ( *, '(a)' ) ' '

      do i = 0, ndata

        if ( i .eq. 0 ) then
          tlo = tdata(1) - 0.5D+00 * ( tdata(2) - tdata(1) )
          thi = tdata(1)
        else if ( i .lt. ndata ) then
          tlo = tdata(i)
          thi = tdata(i+1)
        else if ( ndata .le. i ) then
          tlo = tdata(ndata)
          thi = tdata(ndata) 
     &      + 0.5D+00 * ( tdata(ndata) - tdata(ndata-1) )
        end if

        if ( i .lt. ndata ) then
          jhi = nsample - 1
        else
          jhi = nsample
        end if

        do j = 0, jhi

          tval = ( dble ( nsample - j ) * tlo   
     &           + dble (           j ) * thi ) 
     &           / dble ( nsample     )

          call basis_matrix_tmp ( left, n, mbasis, ndata, tdata, ydata, 
     &      tval, yval )

          if ( 0 .lt. i .and. j .eq. 0 ) then
            mark = '*'
          else
            mark = ' '
          end if

          write ( *, '(2x,a1,2g14.6)' ) mark, tval, yval

        end do

      end do

      return
      end
      subroutine test07 ( )

c*********************************************************************72
c
cc TEST07 tests BASIS_MATRIX_OVERHAUSER_UNI and BASIS_MATRIX_TMP.
c
c  Discussion:
c
c   YDATA(1:NDATA) = ( TDATA(1:NDATA) + 2 )**2 + 3
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    17 June 2006
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n
      parameter ( n = 4 )
      integer ndata
      parameter ( ndata = 4 )

      integer i
      integer j
      integer jhi
      integer left
      character mark
      double precision mbasis(n,n)
      integer nsample
      parameter ( nsample = 4 )
      double precision tdata(ndata)
      double precision thi
      double precision tlo
      double precision tval
      double precision ydata(ndata)
      double precision yval

      save tdata
      save ydata

      data tdata /
     &  -1.0D+00, 0.0D+00, 1.0D+00, 2.0D+00 /
      data ydata /
     &  4.0D+00, 7.0D+00, 12.0D+00, 19.0D+00 /

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST07'
      write ( *, '(a)' ) 
     &  '  BASIS_MATRIX_OVERHAUSER_UNI sets up the basis'
      write ( *, '(a)' ) '    matrix for the uniform Overhauser spline.'

      call basis_matrix_overhauser_uni ( mbasis )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '       TDATA         YDATA'
      write ( *, '(a)' ) ' '
      do i = 1, ndata
        write ( *, '(2x,2g14.6)' ) tdata(i), ydata(i)
      end do

      left = 2

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '        T            Spline(T)'
      write ( *, '(a)' ) ' '

      do i = 0, ndata

        if ( i .eq. 0 ) then
          tlo = tdata(1) - 0.5D+00 * ( tdata(2) - tdata(1) )
          thi = tdata(1)
        else if ( i .lt. ndata ) then
          tlo = tdata(i)
          thi = tdata(i+1)
        else if ( ndata .le. i ) then
          tlo = tdata(ndata)
          thi = tdata(ndata) 
     &      + 0.5D+00 * ( tdata(ndata) - tdata(ndata-1) )
        end if

        if ( i .lt. ndata ) then
          jhi = nsample - 1
        else
          jhi = nsample
        end if

        do j = 0, jhi

          tval = ( dble ( nsample - j ) * tlo   
     &           + dble (           j ) * thi ) 
     &           / dble ( nsample     )

          call basis_matrix_tmp ( left, n, mbasis, ndata, tdata, ydata, 
     &      tval, yval )

          if ( 0 .lt. i .and. j .eq. 0 ) then
            mark = '*'
          else
            mark = ' '
          end if

          write ( *, '(2x,a1,2g14.6)' ) mark, tval, yval

        end do

      end do

      return
      end
      subroutine test08 ( )

c*********************************************************************72
c
cc TEST08 tests BASIS_MATRIX_OVERHAUSER_NONUNI and BASIS_MATRIX_TMP.
c
c  Discussion:
c
c    YDATA(1:NDATA) = ( TDATA(1:NDATA) - 2 )**2 + 3
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    17 June 2006
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n
      parameter ( n = 4 )
      integer ndata
      parameter ( ndata = 4 )

      double precision alpha
      double precision beta
      integer i
      integer j
      integer jhi
      integer left
      character mark
      double precision mbasis(n,n)
      integer nsample
      parameter ( nsample = 4 )
      double precision tdata(ndata)
      double precision thi
      double precision tlo
      double precision tval
      double precision ydata(ndata)
      double precision yval

      save tdata

      data tdata /
     & 0.0D+00, 1.0D+00, 2.0D+00, 3.0D+00 /

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST08'
      write ( *, '(a)' ) 
     &  '  BASIS_MATRIX_OVERHAUSER_NONUNI sets up the'
      write ( *, '(a)' ) 
     &  '    basis matrix for the nonuniform Overhauser'
      write ( *, '(a)' ) '    spline.'

      alpha = ( tdata(3) - tdata(2) ) / ( tdata(3) - tdata(1) )
      beta = ( tdata(3) - tdata(2) ) / ( tdata(4) - tdata(2) )

      write ( *, '(a)' ) ' '
      write ( *, '(a,g14.6)' ) '  ALPHA = ', alpha
      write ( *, '(a,g14.6)' ) '  BETA = ', beta

      call basis_matrix_overhauser_nonuni ( alpha, beta, mbasis )

      do i = 1, ndata
        ydata(i) = ( tdata(i) - 2.0D+00 )**2 + 3.0D+00
      end do

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '       TDATA         YDATA'
      write ( *, '(a)' ) ' '
      do i = 1, ndata
        write ( *, '(2x,2g14.6)' ) tdata(i), ydata(i)
      end do

      left = 2

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '       T           Spline(T)'
      write ( *, '(a)' ) ' '

      do i = 0, ndata

        if ( i .eq. 0 ) then
          tlo = tdata(1) - 0.5D+00 * ( tdata(2) - tdata(1) )
          thi = tdata(1)
        else if ( i .lt. ndata ) then
          tlo = tdata(i)
          thi = tdata(i+1)
        else if ( ndata .le. i ) then
          tlo = tdata(ndata)
          thi = tdata(ndata) 
     &      + 0.5D+00 * ( tdata(ndata) - tdata(ndata-1) )
        end if

        if ( i .lt. ndata ) then
          jhi = nsample - 1
        else
          jhi = nsample
        end if

        do j = 0, jhi

          tval = ( dble ( nsample - j ) * tlo   
     &           + dble (           j ) * thi ) 
     &           / dble ( nsample     )

          call basis_matrix_tmp ( left, n, mbasis, ndata, tdata, ydata, 
     &      tval, yval )

          if ( 0 .lt. i .and. j .eq. 0 ) then
            mark = '*'
          else
            mark = ' '
          end if

          write ( *, '(2x,a1,2g14.6)' ) mark, tval, yval

        end do

      end do

      tdata(1:4) = (/ 0.0D+00, 1.0D+00, 2.0D+00, 5.0D+00 /)

      alpha = ( tdata(3) - tdata(2) ) / ( tdata(3) - tdata(1) )
      beta = ( tdata(3) - tdata(2) ) / ( tdata(4) - tdata(2) )

      write ( *, '(a)' ) ' '
      write ( *, '(a,g14.6)' ) '  ALPHA = ', alpha
      write ( *, '(a,g14.6)' ) '  BETA = ', beta

      call basis_matrix_overhauser_nonuni ( alpha, beta, mbasis )

      do i = 1, ndata
        ydata(i) = ( tdata(i) - 2.0D+00 )**2 + 3.0D+00
      end do

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '    TDATA, YDATA'
      write ( *, '(a)' ) ' '
      do i = 1, ndata
        write ( *, '(2x,2g14.6)' ) tdata(i), ydata(i)
      end do

      left = 2

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '    T, Spline(T)'
      write ( *, '(a)' ) ' '

      do i = 0, ndata

        if ( i .eq. 0 ) then
          tlo = tdata(1) - 0.5D+00 * ( tdata(2) - tdata(1) )
          thi = tdata(1)
        else if ( i .lt. ndata ) then
          tlo = tdata(i)
          thi = tdata(i+1)
        else if ( ndata .le. i ) then
          tlo = tdata(ndata)
          thi = tdata(ndata) 
     &      + 0.5D+00 * ( tdata(ndata) - tdata(ndata-1) )
        end if

        if ( i .lt. ndata ) then
          jhi = nsample - 1
        else
          jhi = nsample
        end if

        do j = 0, jhi

          tval = ( dble ( nsample - j ) * tlo   
     &           + dble (           j ) * thi ) 
     &           / dble ( nsample     )

          call basis_matrix_tmp ( left, n, mbasis, ndata, tdata, ydata, 
     &      tval, yval )

          if ( 0 .lt. i .and. j .eq. 0 ) then
            mark = '*'
          else
            mark = ' '
          end if

          write ( *, '(2x,a1,2g14.6)' ) mark, tval, yval

        end do

      end do

      tdata(1:4) = (/ 0.0D+00, 3.0D+00, 4.0D+00, 5.0D+00 /)

      alpha = ( tdata(3) - tdata(2) ) / ( tdata(3) - tdata(1) )
      beta = ( tdata(3) - tdata(2) ) / ( tdata(4) - tdata(2) )

      write ( *, '(a)' ) ' '
      write ( *, '(a,g14.6)' ) '  ALPHA = ', alpha
      write ( *, '(a,g14.6)' ) '  BETA = ', beta

      call basis_matrix_overhauser_nonuni ( alpha, beta, mbasis )

      do i = 1, ndata
        ydata(i) = ( tdata(i) - 2.0D+00 )**2 + 3.0D+00
      end do

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '    TDATA, YDATA'
      write ( *, '(a)' ) ' '
      do i = 1, ndata
        write ( *, '(2x,2g14.6)' ) tdata(i), ydata(i)
      end do

      left = 2

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '    T, Spline(T)'
      write ( *, '(a)' ) ' '

      do i = 0, ndata

        if ( i .eq. 0 ) then
          tlo = tdata(1) - 0.5D+00 * ( tdata(2) - tdata(1) )
          thi = tdata(1)
        else if ( i .lt. ndata ) then
          tlo = tdata(i)
          thi = tdata(i+1)
        else if ( ndata .le. i ) then
          tlo = tdata(ndata)
          thi = tdata(ndata) 
     &      + 0.5D+00 * ( tdata(ndata) - tdata(ndata-1) )
        end if

        if ( i .lt. ndata ) then
          jhi = nsample - 1
        else
          jhi = nsample
        end if

        do j = 0, jhi

          tval = ( dble ( nsample - j ) * tlo   
     &           + dble (           j ) * thi ) 
     &           / dble ( nsample     )

          call basis_matrix_tmp ( left, n, mbasis, ndata, tdata, ydata, 
     &      tval, yval )

          if ( 0 .lt. i .and. j .eq. 0 ) then
            mark = '*'
          else
            mark = ' '
          end if

          write ( *, '(2x,a1,2g14.6)' ) mark, tval, yval

        end do

      end do

      return
      end
      subroutine test09 ( )

c*********************************************************************72
c
cc TEST09 tests BASIS_MATRIX_OVERHAUSER_NONUNI and BASIS_MATRIX_TMP.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    17 June 2006
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n
      parameter ( n = 4 )
      integer ndata
      parameter ( ndata = 4 )

      double precision alpha
      double precision beta
      integer i
      integer j
      integer jhi
      integer left
      character mark
      double precision mbasis(n,n)
      integer nsample
      parameter ( nsample = 4 )
      double precision tdata(ndata)
      double precision thi
      double precision tlo
      double precision tval
      double precision ydata(ndata)
      double precision yval

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST09'
      write ( *, '(a)' ) 
     &  '  BASIS_MATRIX_OVERHAUSER_NONUNI sets up the'
      write ( *, '(a)' ) 
     &  '    basis matrix for the nonuniform Overhauser '
      write ( *, '(a)' ) '    spline.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 
     &  '  First test that the nonuniform code can match'
      write ( *, '(a)' ) 
     &  '  the uniform code.  Compare these results with'
      write ( *, '(a)' ) '  the uniform output.'
      write ( *, '(a)' ) ' '

      tdata(1:4) = (/ -1.0D+00, 0.0D+00, 1.0D+00, 2.0D+00 /)

      alpha = ( tdata(3) - tdata(2) ) / ( tdata(3) - tdata(1) )
      beta =  ( tdata(3) - tdata(2) ) / ( tdata(4) - tdata(2) )

      write ( *, '(a)' ) ' '
      write ( *, '(a,g14.6)' ) '  ALPHA = ', alpha
      write ( *, '(a,g14.6)' ) '  BETA = ', beta

      call basis_matrix_overhauser_nonuni ( alpha, beta, mbasis )

      do i = 1, ndata
        ydata(i) = ( tdata(i) + 2.0D+00 )**2 + 3.0D+00
      end do

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '    TDATA, YDATA'
      write ( *, '(a)' ) ' '
      do i = 1, ndata
        write ( *, '(2x,2g14.6)' ) tdata(i), ydata(i)
      end do

      left = 2

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '    T, Spline(T)'
      write ( *, '(a)' ) ' '

      do i = 0, ndata

        if ( i .eq. 0 ) then
          tlo = tdata(1) - 0.5D+00 * ( tdata(2) - tdata(1) )
          thi = tdata(1)
        else if ( i .lt. ndata ) then
          tlo = tdata(i)
          thi = tdata(i+1)
        else if ( ndata .le. i ) then
          tlo = tdata(ndata)
          thi = tdata(ndata) 
     &      + 0.5D+00 * ( tdata(ndata) - tdata(ndata-1) )
        end if

        if ( i .lt. ndata ) then
          jhi = nsample - 1
        else
          jhi = nsample
        end if

        do j = 0, jhi

          tval = ( dble ( nsample - j ) * tlo   
     &           + dble (           j ) * thi ) 
     &           / dble ( nsample     )

          call basis_matrix_tmp ( left, n, mbasis, ndata, tdata, ydata, 
     &      tval, yval )

          if ( 0 .lt. i .and. j .eq. 0 ) then
            mark = '*'
          else
            mark = ' '
          end if

          write ( *, '(2x,a1,2g14.6)' ) mark, tval, yval

        end do

      end do

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Now test that the nonuniform code on a'
      write ( *, '(a)' ) '  nonuniform grid.'
      write ( *, '(a)' ) ' '

      tdata(1:4) = (/ -4.0D+00, -3.0D+00, -1.0D+00, 2.0D+00 /)

      alpha = ( tdata(3) - tdata(2) ) / ( tdata(3) - tdata(1) )
      beta =  ( tdata(3) - tdata(2) ) / ( tdata(4) - tdata(2) )

      write ( *, '(a)' ) ' '
      write ( *, '(a,g14.6)' ) '  ALPHA = ', alpha
      write ( *, '(a,g14.6)' ) '  BETA =  ', beta

      call basis_matrix_overhauser_nonuni ( alpha, beta, mbasis )

      ydata(1:ndata) = ( tdata(1:ndata) + 2.0D+00 )**2 + 3.0D+00

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '       TDATA         YDATA'
      write ( *, '(a)' ) ' '
      do i = 1, ndata
        write ( *, '(2x,2g14.6)' ) tdata(i), ydata(i)
      end do

      left = 2

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '        T            Spline(T)'
      write ( *, '(a)' ) ' '

      do i = 0, ndata

        if ( i .eq. 0 ) then
          tlo = tdata(1) - 0.5D+00 * ( tdata(2) - tdata(1) )
          thi = tdata(1)
        else if ( i .lt. ndata ) then
          tlo = tdata(i)
          thi = tdata(i+1)
        else if ( ndata .le. i ) then
          tlo = tdata(ndata)
          thi = tdata(ndata) 
     &      + 0.5D+00 * ( tdata(ndata) - tdata(ndata-1) )
        end if

        if ( i .lt. ndata ) then
          jhi = nsample - 1
        else
          jhi = nsample
        end if

        do j = 0, jhi

          tval = ( dble ( nsample - j ) * tlo   
     &           + dble (           j ) * thi ) 
     &           / dble ( nsample     )

          call basis_matrix_tmp ( left, n, mbasis, ndata, tdata, ydata, 
     &      tval, yval )

          if ( 0 .lt. i .and. j .eq. 0 ) then
            mark = '*'
          else
            mark = ' '
          end if

          write ( *, '(2x,a1,2g14.6)' ) mark, tval, yval

        end do

      end do

      return
      end
      subroutine test10 ( )

c*********************************************************************72
c
cc TEST10 tests BC_VAL.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    17 June 2006
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n
      parameter ( n = 2 )

      integer i
      integer nsample
      parameter ( nsample = 101 )
      double precision t
      double precision xcon(0:n)
      double precision xval
      double precision ycon(0:n)
      double precision yval

      save xcon
      save ycon

      data xcon /
     & 0.0D+00, 0.75D+00, 1.0D+00 /
      data ycon /
     & 1.0D+00, 0.0D+00,  1.0D+00 /

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST10'
      write ( *, '(a)' ) '  BC_VAL evaluates a general Bezier function.'
c
c  One point on the curve should be about (0.75, 0.536).
c
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '       T            X(T)          Y(T)'
      write ( *, '(a)' ) ' '

      do i = 1, nsample
        t = dble ( i - 1 ) / dble ( nsample - 1 )
        call bc_val ( n, t, xcon, ycon, xval, yval )
        write ( *, '(2x,3g14.6)' ) t, xval, yval
      end do

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 
     &  '  The point ( 0.75, 0.536 ) should be on the curve.'

      return
      end
      subroutine test11 ( )

c*********************************************************************72
c
cc TEST11 tests BEZ_VAL.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    17 June 2006
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n
      parameter ( n = 2 )

      double precision a
      parameter ( a = 0.0D+00 )
      double precision b
      parameter ( b = 1.0D+00 )
      double precision bez_val
      integer i
      integer nsample
      parameter ( nsample = 21 )
      double precision x
      double precision y(0:n)

      save y

      data y /
     & 1.0D+00, 0.0D+00, 1.0D+00 /

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST11'
      write ( *, '(a)' ) '  BEZ_VAL evaluates a Bezier function.'
c
c  One point on the curve should be (0.75, 20/32).
c
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '         T    X(T)          Y(T)'
      write ( *, '(a)' ) ' '

      do i = 1, nsample

        x = ( dble ( nsample - i     ) * a   
     &      + dble (           i - 1 ) * b ) 
     &      / dble ( nsample     - 1 )

        write ( *, '(2x,i8,2g14.6)' ) i, x, bez_val ( n, x, a, b, y )

      end do

      write ( *, '(a)' ) ' '
      write ( *, '(a,g14.6)' ) '  When X = ', 0.75D+00
      write ( *, '(a,g14.6)' ) '  BEZ_VAL(X) should be ', 0.625D+00

      return
      end
      subroutine test115 ( )

c*********************************************************************72
c
cc TEST115 tests BP01.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    17 June 2006
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n_max 
      parameter ( n_max = 3 )

      double precision a
      parameter ( a = 0.0D+00 )
      double precision b
      parameter ( b = 1.0D+00 )
      double precision bern(0:n_max)
      integer i
      integer n
      integer nsample
      parameter ( nsample = 11 )
      double precision x

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST115'
      write ( *, '(a)' ) 
     &  '  BP01 evaluates the Bernstein basis polynomials'
      write ( *, '(a)' ) '  for the interval [0,1].'

      do n = 0, n_max

        write ( *, '(a)' ) ' '
        write ( *, '(a,i8)' ) '  Degree N = ', n
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 
     &    '   X         BERN(N,0,X)  BERN(N,1,X)  ' // 
     &    'BERN(N,2,X)  BERN(N,3,X)'
        write ( *, '(a)' ) ' '

        do i = 1, nsample

          x = ( dble ( nsample - i     ) * a   
     &        + dble (           i - 1 ) * b ) 
     &        / dble ( nsample     - 1 )

          call bp01 ( n, x, bern )

          write ( *, '(2x,f8.4,4x,5g14.6)' ) x, bern(0:n)

        end do

      end do

      return
      end
      subroutine test116 ( )

c*********************************************************************72
c
cc TEST116 tests BPAB.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    17 June 2006
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n_max
      parameter ( n_max = 3 )

      double precision a
      parameter ( a = 1.0D+00 )
      double precision b
      parameter ( b = 3.0D+00 )
      double precision bern(0:n_max)
      integer i
      integer n
      double precision x
      integer nsample
      parameter ( nsample = 11 )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST116'
      write ( *, '(a)' ) 
     &  '  BPAB evaluates the Bernstein basis polynomials'
      write ( *, '(a)' ) '  for the interval [A,B].'
      write ( *, '(a)' ) ' '
      write ( *, '(a,g14.6)' ) '  A = ', a
      write ( *, '(a,g14.6)' ) '  B = ', b

      do n = 0, n_max

        write ( *, '(a)' ) ' '
        write ( *, '(a,i8)' ) '  Degree N = ', n
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 
     &    '   X         BERN(N,0,X)  BERN(N,1,X)  ' // 
     &    'BERN(N,2,X)  BERN(N,3,X)'
        write ( *, '(a)' ) ' '

        do i = 1, nsample

          x = ( dble ( nsample - i     ) * a   
     &        + dble (           i - 1 ) * b ) 
     &        / dble ( nsample     - 1 )

          call bpab ( n, a, b, x, bern )

          write ( *, '(2x,f8.4,4x,5g14.6)' ) x, bern(0:n)

        end do

      end do

      return
      end
      subroutine test12 ( )

c*********************************************************************72
c
cc TEST12 tests BPAB_APPROX.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    17 June 2006
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer maxdata
      parameter ( maxdata = 10 )

      double precision a
      double precision b
      integer i
      integer ndata
      integer nsample
      double precision xdata(0:maxdata)
      double precision xval
      double precision ydata(0:maxdata)
      double precision yval

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST12'
      write ( *, '(a)' ) 
     &  '  BPAB_APPROX evaluates the Bernstein polynomial'
      write ( *, '(a)' ) '  approximant to a function F(X).'

      a = 1.0D+00
      b = 3.0D+00

      do ndata = 0, 9, 3

        do i = 0, ndata

          if ( ndata .eq. 0 ) then
            xdata(i) = 0.5D+00 * ( a + b )
          else
            xdata(i) = ( dble ( ndata - i ) * a   
     &                 + dble (         i ) * b ) 
     &                 / dble ( ndata     )
          end if

          ydata(i) = sin ( xdata(i) )

        end do

        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) '       XDATA        YDATA'
        write ( *, '(a)' ) ' '
        do i = 0, ndata
          write ( *, '(2x,2g14.6)' ) xdata(i), ydata(i)
        end do

        write ( *, '(a)' ) ' '
        write ( *, '(a,i8)' ) 
     &    '  Bernstein approximant of degree N = ', ndata
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 
     &    '       X            F(X)          BERN(X)        ERROR'
        write ( *, '(a)' ) ' '

        nsample = 2 * ndata + 1

        do i = 1, nsample

          if ( nsample .eq. 1 ) then
            xval = 0.5D+00 * ( a + b )
          else
            xval = ( dble ( nsample - i     ) * a   
     &             + dble (           i - 1 ) * b ) 
     &             / dble ( nsample     - 1 )
          end if

          call bpab_approx ( ndata, a, b, ydata, xval, yval )

          write ( *, '(2x,4g14.6)' ) 
     &      xval, sin(xval), yval, yval - sin(xval)

        end do

      end do

      return
      end
      subroutine test125 ( )

c*********************************************************************72
c
cc TEST125 tests LEAST_SET_OLD and LEAST_VAL_OLD.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    17 June 2006
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer maxdeg
      parameter ( maxdeg = 6 )
      integer ntab
      parameter ( ntab = 21 )

      double precision b(1:maxdeg)
      double precision c(0:maxdeg)
      double precision d(2:maxdeg)
      double precision eps
      double precision error
      integer i
      integer ierror
      integer j
      integer jhi
      integer ndeg
      double precision ptab(ntab)
      double precision xtab(ntab)
      double precision xval
      double precision ytab(ntab)
      double precision ytrue
      double precision yval

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST125'
      write ( *, '(a)' ) 
     &  '  LEAST_SET_OLD sets a least squares polynomial,'
      write ( *, '(a)' ) '  LEAST_VAL_OLD evaluates it.'

      do i = 1, ntab
        xtab(i) = ( dble ( ntab - i     ) * ( -1.0D+00 )   
     &            + dble (        i - 1 ) * ( +1.0D+00 ) ) 
     &            / dble ( ntab     - 1 )
        ytab(i) = dble ( int ( exp ( xtab(i) ) * 100.0D+00 + 0.5D+00 ) ) 
     &    / 100.0D+00
      end do

      write ( *, '(a)'    ) ' '
      write ( *, '(a)'    ) '  The data to be interpolated:'
      write ( *, '(a)'    ) ' '
      write ( *, '(a,i8)' ) '  Number of data values = ', ntab
      write ( *, '(a)'    ) ' '
      write ( *, '(a)'    ) '       X             Y'
      write ( *, '(a)'    ) ' '

      do i = 1, ntab
        write ( *, '(2x,2g14.6)' ) xtab(i), ytab(i)
      end do

      do ndeg = 1, maxdeg

        write ( *, '(a)' ) ' '
        write ( *, '(a,i8)' ) '  Using a polynomial of degree: ', ndeg
        write ( *, '(a)' ) ' '

        call least_set_old ( ntab, xtab, ytab, ndeg, ptab, b, c, d, 
     &    eps, ierror )

        write ( *, '(a)' ) ' '
        write ( *, '(a,g14.6)' ) '  Total approximation error = ', eps
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 
     &    '      X            F(X)          P(X)         Error'
        write ( *, '(a)' ) ' '

        do i = 1, ntab

          if ( i .lt. ntab ) then
            jhi = 2
          else
            jhi = 0
          end if

          do j = 0, jhi

            if ( i .lt. ntab ) then

              xval = ( dble ( 3 - j ) * xtab(i)     
     &               + dble (     j ) * xtab(i+1) ) 
     &               / dble ( 3     )

            else

              xval = xtab(i)

            end if

            call least_val_old ( xval, ndeg, b, c, d, yval )

            ytrue = dble ( int ( exp ( xval ) * 100.0D+00 + 0.5D+00 ) ) 
     &        / 100.0D+00

            error = yval - ytrue
            write ( *, '(2x,5g14.6)' ) xval, yval, ytrue, error
          end do

        end do

      end do

      return
      end
      subroutine test126 ( )

c*********************************************************************72
c
cc TEST126 tests LEAST_SET and LEAST_VAL.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    09 December 2006
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer point_num
      parameter ( point_num = 21 )
      integer nterms
      parameter ( nterms = 4 )

      double precision b(nterms)
      double precision c(nterms)
      double precision d(nterms)
      double precision f(point_num)
      integer i
      integer nterms2
      double precision px
      double precision w(point_num)
      double precision x(point_num)

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST126'
      write ( *, '(a)' ) '  LEAST_SET sets a least squares polynomial,'
      write ( *, '(a)' ) '  LEAST_VAL evaluates it.'

      w(1:point_num) = 1.0D+00

      do i = 1, point_num
        x(i) = - 1.0D+00 + dble ( i - 1 ) / 10.0D+00
        f(i) = dble ( int ( exp ( x(i) ) * 100.0D+00 + 0.5D+00 ) ) 
     &    / 100.0D+00
      end do

      call least_set ( point_num, x, f, w, nterms, b, c, d )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  X, F(X), P(X), Error'
      write ( *, '(a)' ) ' '

      do nterms2 = 1, nterms
        write ( *, '(a)' ) ' '
        write ( *, '(a,i8)' ) '  Using polynomial order = ', nterms2
        write ( *, '(a)' ) ' '
        do i = 1, point_num
          call least_val ( nterms2, b, c, d, x(i), px )
          write ( *, '(5g14.6)' ) x(i), f(i), px, px - f(i)
        end do
      end do

      return
      end
      subroutine test127 ( )

c*********************************************************************72
c
cc TEST127 tests LEAST_SET and LEAST_VAL2.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    09 December 2006
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer point_num
      parameter ( point_num = 21 )
      integer nterms
      parameter ( nterms = 4 )

      double precision b(nterms)
      double precision c(nterms)
      double precision d(nterms)
      double precision f(point_num)
      double precision fp(point_num)
      integer i
      integer nterms2
      double precision px
      double precision pxp
      double precision w(point_num)
      double precision x(point_num)

      w(1:point_num) = 1.0D+00

      do i = 1, point_num
        x(i) = -1.0D+00 + dble ( i - 1 ) / 10.0D+00
        f(i) = x(i)**2 - x(i) - 6.0D+00
        fp(i) = 2.0D+00 * x(i) - 1.0D+00
      end do

      call least_set ( point_num, x, f, w, nterms, b, c, d )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST127'
      write ( *, '(a)' ) '  LEAST_SET sets a least squares polynomial,'
      write ( *, '(a)' ) '  LEAST_VAL2 evaluates it.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  X, F(X), P(X), FP(X), PP(X)'
      write ( *, '(a)' ) ' '

      do nterms2 = 1, nterms
        write ( *, '(a)' ) ' '
        write ( *, '(a,i8)' ) '  Using polynomial order = ', nterms2
        write ( *, '(a)' ) ' '
        do i = 1, point_num
          call least_val2 ( nterms2, b, c, d, x(i), px, pxp )
          write ( *, '(5g14.6)' ) x(i), f(i), px, fp(i), pxp
        end do
      end do

      return
      end
      subroutine test13 ( )

c*********************************************************************72
c
cc TEST13 tests SPLINE_B_VAL.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    17 June 2006
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer ndata 
      parameter ( ndata = 11 )

      integer i
      integer j
      integer jhi
      character mark
      integer nsample
      parameter ( nsample = 4 )
      double precision pi
      parameter ( pi = 3.141592653589793D+00 )
      double precision tdata(ndata)
      double precision thi
      double precision tlo
      double precision tval
      double precision ydata(ndata)
      double precision yval

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST13'
      write ( *, '(a)' ) '  SPLINE_B_VAL evaluates the B spline.'

      do i = 1, ndata
        tdata(i) = dble ( i - 1 )
        ydata(i) = sin ( 2.0D+00 * pi * tdata(i) / dble ( ndata - 1 ) )
      end do

      write ( *, '(a)'    ) ' '
      write ( *, '(a)'    ) '  The data to be interpolated:'
      write ( *, '(a)'    ) ' '
      write ( *, '(a,i8)' ) '  Number of data values = ', ndata
      write ( *, '(a)'    ) ' '
      write ( *, '(a)'    ) '       T             Y'
      write ( *, '(a)'    ) ' '

      do i = 1, ndata
        write ( *, '(2x,2g14.6)' ) tdata(i), ydata(i)
      end do

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '       T           Spline(T)'
      write ( *, '(a)' ) ' '

      do i = 0, ndata

        if ( i .eq. 0 ) then
          tlo = tdata(1) - 0.5D+00 * ( tdata(2) - tdata(1) )
          thi = tdata(1)
        else if ( i .lt. ndata ) then
          tlo = tdata(i)
          thi = tdata(i+1)
        else if ( ndata .le. i ) then
          tlo = tdata(ndata)
          thi = tdata(ndata) 
     &      + 0.5D+00 * ( tdata(ndata) - tdata(ndata-1) )
        end if

        if ( i .lt. ndata ) then
          jhi = nsample - 1
        else
          jhi = nsample
        end if

        do j = 0, jhi

          tval = ( dble ( nsample - j ) * tlo   
     &           + dble (           j ) * thi ) 
     &           / dble ( nsample     )

          call spline_b_val ( ndata, tdata, ydata, tval, yval )

          if ( 0 .lt. i .and. j .eq. 0 ) then
            mark = '*'
          else
            mark = ' '
          end if

          write ( *, '(2x,a1,2g14.6)' ) mark, tval, yval

        end do

      end do

      return
      end
      subroutine test14 ( )

c*********************************************************************72
c
cc TEST14 tests SPLINE_BETA_VAL.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    17 June 2006
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer ndata
      parameter ( ndata = 11 )

      double precision beta1
      double precision beta2
      integer i
      integer j
      integer jhi
      character mark
      integer nsample
      parameter ( nsample = 4 )
      double precision pi
      parameter ( pi = 3.141592653589793D+00 )
      double precision tdata(ndata)
      double precision thi
      double precision tlo
      double precision tval
      double precision ydata(ndata)
      double precision yval

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST14'
      write ( *, '(a)' ) '  SPLINE_BETA_VAL evaluates the BETA spline.'

      do i = 1, ndata
        tdata(i) = dble ( i - 1 )
        ydata(i) = sin ( 2.0D+00 * pi * tdata(i) / dble ( ndata - 1 ) )
      end do

      write ( *, '(a)'    ) ' '
      write ( *, '(a)'    ) '  The data to be interpolated:'
      write ( *, '(a)'    ) ' '
      write ( *, '(a,i8)' ) '  Number of data values = ', ndata
      write ( *, '(a)'    ) ' '
      write ( *, '(a)'    ) '       T             Y'
      write ( *, '(a)'    ) ' '

      do i = 1, ndata
        write ( *, '(2x,2g14.6)' ) tdata(i), ydata(i)
      end do

      beta1 = 1.0D+00
      beta2 = 0.0D+00

      write ( *, '(a)' ) ' '
      write ( *, '(a,g14.6)' ) '  BETA1 = ', beta1
      write ( *, '(a,g14.6)' ) '  BETA2 = ', beta2
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '    T, Spline(T)'
      write ( *, '(a)' ) ' '

      do i = 0, ndata

        if ( i .eq. 0 ) then
          tlo = tdata(1) - 0.5D+00 * ( tdata(2) - tdata(1) )
          thi = tdata(1)
        else if ( i .lt. ndata ) then
          tlo = tdata(i)
          thi = tdata(i+1)
        else if ( ndata .le. i ) then
          tlo = tdata(ndata)
          thi = tdata(ndata) 
     &      + 0.5D+00 * ( tdata(ndata) - tdata(ndata-1) )
        end if

        if ( i .lt. ndata ) then
          jhi = nsample - 1
        else
          jhi = nsample
        end if

        do j = 0, jhi

          tval = ( dble ( nsample - j ) * tlo   
     &           + dble (           j ) * thi ) 
     &           / dble ( nsample     )

          call spline_beta_val ( beta1, beta2, ndata, tdata, ydata, 
     &      tval, yval )

          if ( 0 .lt. i .and. j .eq. 0 ) then
            mark = '*'
          else
            mark = ' '
          end if

          write ( *, '(2x,a1,2g14.6)' ) mark, tval, yval

        end do

      end do

      beta1 = 1.0D+00
      beta2 = 100.0D+00

      write ( *, '(a)' ) ' '
      write ( *, '(a,g14.6)' ) '  BETA1 = ', beta1
      write ( *, '(a,g14.6)' ) '  BETA2 = ', beta2
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '    T, Spline(T)'
      write ( *, '(a)' ) ' '

      do i = 0, ndata

        if ( i .eq. 0 ) then
          tlo = tdata(1) - 0.5D+00 * ( tdata(2) - tdata(1) )
          thi = tdata(1)
        else if ( i .lt. ndata ) then
          tlo = tdata(i)
          thi = tdata(i+1)
        else if ( ndata .le. i ) then
          tlo = tdata(ndata)
          thi = tdata(ndata) 
     &      + 0.5D+00 * ( tdata(ndata) - tdata(ndata-1) )
        end if

        if ( i .lt. ndata ) then
          jhi = nsample - 1
        else
          jhi = nsample
        end if

        do j = 0, jhi

          tval = ( dble ( nsample - j ) * tlo   
     &           + dble (           j ) * thi ) 
     &           / dble ( nsample     )

          call spline_beta_val ( beta1, beta2, ndata, tdata, ydata, 
     &      tval, yval )

          if ( 0 .lt. i .and. j .eq. 0 ) then
            mark = '*'
          else
            mark = ' '
          end if

          write ( *, '(2x,a1,2g14.6)' ) mark, tval, yval

        end do

      end do

      beta1 = 100.0D+00
      beta2 = 0.0D+00

      write ( *, '(a)' ) ' '
      write ( *, '(a,g14.6)' ) '  BETA1 = ', beta1
      write ( *, '(a,g14.6)' ) '  BETA2 = ', beta2
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '    T, Spline(T)'
      write ( *, '(a)' ) ' '

      do i = 0, ndata

        if ( i .eq. 0 ) then
          tlo = tdata(1) - 0.5D+00 * ( tdata(2) - tdata(1) )
          thi = tdata(1)
        else if ( i .lt. ndata ) then
          tlo = tdata(i)
          thi = tdata(i+1)
        else if ( ndata .le. i ) then
          tlo = tdata(ndata)
          thi = tdata(ndata) 
     &      + 0.5D+00 * ( tdata(ndata) - tdata(ndata-1) )
        end if

        if ( i .lt. ndata ) then
          jhi = nsample - 1
        else
          jhi = nsample
        end if

        do j = 0, jhi

          tval = ( dble ( nsample - j ) * tlo   
     &           + dble (           j ) * thi ) 
     &           / dble ( nsample     )

          call spline_beta_val ( beta1, beta2, ndata, tdata, ydata, 
     &      tval, yval )

          if ( 0 .lt. i .and. j .eq. 0 ) then
            mark = '*'
          else
            mark = ' '
          end if

          write ( *, '(2x,a1,2g14.6)' ) mark, tval, yval

        end do

      end do

      return
      end
      subroutine test143 ( )

c*********************************************************************72
c
cc TEST143 tests SPLINE_BEZIER_VAL.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    17 June 2006
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer dim_num
      parameter ( dim_num = 1 )
      integer interval_num
      parameter ( interval_num = 3 )
      integer point_num
      parameter ( point_num = 6 * interval_num + 1 )

      double precision data_val(dim_num,3*interval_num+1)
      double precision dxdt
      integer interval
      integer j
      character mark
      integer point
      double precision point_t(point_num)
      double precision point_val(dim_num,point_num)
      double precision t
      double precision t_max
      double precision t_min
      double precision x
      double precision x_max
      parameter ( x_max = 2.0D+00 * 3.141592653589793D+00 )
      double precision x_min 
      parameter ( x_min = 0.0D+00 )
      double precision xdata(0:interval_num)
      double precision ydata(0:interval_num)
      double precision ypdata(0:interval_num)

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST143'
      write ( *, '(a)' ) 
     &  '  SPLINE_BEZIER_VAL evaluates a cubic Bezier spline.'
c
c  Construct the data.
c
      do interval = 0, interval_num

        x = ( dble ( interval_num - interval ) * x_min   
     &      + dble (                interval ) * x_max ) 
     &      / dble ( interval_num            )

        xdata(interval) = x
        ydata(interval) =  sin ( x )
        ypdata(interval) = cos ( x )

      end do

      write ( *, '(a)'    ) ' '
      write ( *, '(a)'    ) '  The data to be interpolated:'
      write ( *, '(a)'    ) ' '
      write ( *, '(a,i8)' ) '  Number of intervals = ', interval_num
      write ( *, '(a)'    ) ' '
      write ( *, '(a)'    ) '       X             Y            dYdX'
      write ( *, '(a)'    ) ' '

      do interval = 0, interval_num
        write ( *, '(2x,3g14.6)' ) 
     &    xdata(interval), ydata(interval), ypdata(interval)
      end do
c
c  Construct the Bezier control data.
c
      dxdt = ( x_max - x_min ) / dble ( interval_num )
      j = 0

      do interval = 1, interval_num

        if ( interval .eq. 1 ) then
          j = j + 1
          data_val(1,j) = ydata(interval-1)
        end if

        j = j + 1
        data_val(1,j) = ydata(interval-1) 
     &    + ypdata(interval-1) * dxdt / 3.0D+00

        j = j + 1
        data_val(1,j) = ydata(interval)   
     &    - ypdata(interval) * dxdt  / 3.0D+00

        j = j + 1
        data_val(1,j) = ydata(interval)

      end do

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  The control points'
      write ( *, '(a)' ) '  Interpolation points are marked with a "*".'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '                     T            P(T)'
      write ( *, '(a)' ) ' '

      do j = 1, 3 * interval_num + 1

        t = dble ( j - 1 ) / 3.0

        if ( abs ( nint ( t ) - t ) .lt. 0.00001D+00 ) then
          mark = '*'
        else
          mark = ' '
        end if

        write ( *, '(2x,a,2x,i8,2g14.6)' ) mark, j, t, data_val(1,j)
      end do

      t_min = 0.0D+00
      t_max = dble ( interval_num )

      do point = 1, point_num

        t = ( dble ( point_num - point     ) * t_min   
     &      + dble (             point - 1 ) * t_max ) 
     &      / dble ( point_num         - 1 )

        point_t(point) = t

      end do

      call spline_bezier_val ( dim_num, interval_num, data_val, 
     &  point_num, point_t, point_val )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 
     &  '  The Bezier spline, sampled at various points.'
      write ( *, '(a)' ) 
     &  '  Interpolation points are marked with a "*".'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 
     &  '                   T        Spline(T)      F(X(T))'
      write ( *, '(a)' ) ' '

      do point = 1, point_num

        if ( abs ( nint ( point_t(point) ) - point_t(point) ) 
     &    .lt. 0.00001D+00 ) then
          mark = '*'
        else
          mark = ' '
        end if

        x = ( dble ( point_num - point     ) * x_min   
     &      + dble (             point - 1 ) * x_max ) 
     &      / dble ( point_num         - 1 )

        write ( *, '(2x,a,2x,i8,2x,f10.6,2g14.6)' ) 
     &    mark, point, point_t(point), point_val(1,point), sin ( x )

      end do

      return
      end
      subroutine test144 ( )

c*********************************************************************72
c
cc TEST144 tests SPLINE_BEZIER_VAL.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    18 June 2006
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer dim_num
      parameter ( dim_num = 1 )
      integer interval_num
      parameter ( interval_num = 3 )
      integer point_num
      parameter ( point_num = 6 * interval_num + 1 )

      double precision data_val(dim_num,3*interval_num+1)
      integer interval
      integer j
      character mark
      integer point
      double precision point_t(point_num)
      double precision point_val(dim_num,point_num)
      double precision t
      double precision t_max
      double precision t_min
      double precision x
      double precision x_max
      parameter ( x_max = 2.0D+00 * 3.141592653589793D+00 )
      double precision x_min
      parameter ( x_min = 0.0D+00 )
      double precision xdata(0:3*interval_num)
      double precision ydata(0:3*interval_num)
      double precision y0
      double precision y1
      double precision y2
      double precision y3
      double precision ypdata(0:interval_num)

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST144'
      write ( *, '(a)' ) 
     &  '  SPLINE_BEZIER_VAL evaluates a cubic Bezier spline.'
      write ( *, '(a)' ) 
     &  '  Normally, the "interior" points of a Bezier spline'
      write ( *, '(a)' ) '  are not interpolating points.  Instead, the'
      write ( *, '(a)' ) 
     &  '  derivatives at the interval endpoints are used.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 
     &  '  This example shows, however, that it is possible'
      write ( *, '(a)' ) 
     &  '  to start with only data values, and to "massage"'
      write ( *, '(a)' ) 
     &  '  the data so that we can define a cubic Bezier spline'
      write ( *, '(a)' ) '  which interpolates ALL the data.'
c
c  Construct the data.
c
      do interval = 0, 3 * interval_num

        x = ( dble ( 3 * interval_num - interval ) * x_min   
     &      + dble (                    interval ) * x_max ) 
     &      / dble ( 3 * interval_num            )

        xdata(interval) = x
        ydata(interval) = sin ( x )

      end do

      write ( *, '(a)'    ) ' '
      write ( *, '(a)'    ) '  The data to be interpolated:'
      write ( *, '(a)'    ) ' '
      write ( *, '(a,i8)' ) '  Number of intervals = ', interval_num
      write ( *, '(a)'    ) ' '
      write ( *, '(a)'    ) '       X             Y'
      write ( *, '(a)'    ) ' '

      do interval = 0, 3 * interval_num
        write ( *, '(2x,3g14.6)' ) xdata(interval), ydata(interval)
      end do
c
c  Construct the Bezier control data.
c  The middle points must be "massaged".
c
      j = 0

      do interval = 1, interval_num

        y0 = ydata(3*(interval-1))
        y1 = ydata(3*(interval-1)+1)
        y2 = ydata(3*(interval-1)+2)
        y3 = ydata(3*(interval-1)+3)

        if ( interval .eq. 1 ) then
          j = j + 1
          data_val(1,j) = y0
        end if

        j = j + 1
        data_val(1,j) = 
     &   ( -5.0D+00 * y0 + 18.0D+00 * y1 - 9.0D+00 * y2 + 2.0D+00 * y3 ) 
     &   / 6.0D+00

        j = j + 1
        data_val(1,j) = 
     &    ( 2.0D+00 * y0 - 9.0D+00 * y1 + 18.0D+00 * y2 - 5.0D+00 * y3 ) 
     &    / 6.0D+00

        j = j + 1
        data_val(1,j) = y3

      end do

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  The control points'
      write ( *, '(a)' ) 
     &  '  ALL control points will be interpolation points!'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '                     T            P(T)'
      write ( *, '(a)' ) ' '

      do j = 1, 3 * interval_num + 1

        t = dble ( j - 1 ) / 3.0

        mark = '*'

        write ( *, '(2x,a,2x,i8,2g14.6)' ) mark, j, t, data_val(1,j)
      end do

      t_min = 0.0D+00
      t_max = dble ( interval_num )

      do point = 1, point_num

        t = ( dble ( point_num - point     ) * t_min   
     &      + dble (             point - 1 ) * t_max ) 
     &      / dble ( point_num         - 1 )

        point_t(point) = t

      end do

      call spline_bezier_val ( dim_num, interval_num, data_val, 
     &  point_num, point_t, point_val )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 
     &  '  The Bezier spline, sampled at various points.'
      write ( *, '(a)' ) '  Interpolation points are marked with a "*".'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 
     &  '                   T        Spline(T)      F(X(T))'
      write ( *, '(a)' ) ' '

      do point = 1, point_num

        if ( abs ( nint ( 3 * point_t(point) ) - 3 * point_t(point) ) 
     &    .lt. 0.00001D+00 ) then
          mark = '*'
        else
          mark = ' '
        end if

        x = ( dble ( point_num - point     ) * x_min   
     &      + dble (             point - 1 ) * x_max ) 
     &      / dble ( point_num         - 1 )

        write ( *, '(2x,a,2x,i8,2x,f10.6,2g14.6)' ) 
     &    mark, point, point_t(point), point_val(1,point), sin ( x )

      end do

      return
      end
      subroutine test145 ( )

c*********************************************************************72
c
cc TEST145 tests SPLINE_CONSTANT_VAL.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    17 June 2006
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer ndata
      parameter ( ndata = 12 )
      integer n_test
      parameter ( n_test = 20 )

      double precision ahi
      double precision alo
      double precision frunge
      double precision fval
      integer i
      integer j
      integer seed
      double precision tdata(ndata-1)
      double precision thi
      double precision tlo
      double precision t_test(n_test)
      double precision tval
      double precision ydata(ndata)
      double precision yval

      seed = 123456789

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST145'
      write ( *, '(a)' ) '  SPLINE_CONSTANT_VAL evaluates a piecewise '
      write ( *, '(a)' ) '  constant spline.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Runge''s function, evenly spaced knots.'
c
c  Set the data.
c
      tlo = -1.0D+00
      thi = +1.0D+00
      call r8vec_even ( ndata-1, tlo, thi, tdata )

      do i = 1, ndata

        if ( i .eq. 1 ) then
          tval = tdata(1)
        else if ( i .lt. ndata ) then
          tval = 0.5D+00 * ( tdata(i-1) + tdata(i) )
        else if ( i .eq. ndata ) then
          tval = tdata(i-1)
        end if

        ydata(i) = frunge ( tval )

      end do

      write ( *, '(a)'    ) ' '
      write ( *, '(a)'    ) '  The data to be interpolated:'
      write ( *, '(a)'    ) ' '
      write ( *, '(a,i8)' ) '  Number of data values = ', ndata
      write ( *, '(a)'    ) ' '
      write ( *, '(a)'    ) '       T             Y'
      write ( *, '(a)'    ) ' '

      do i = 1, ndata
        write ( *, '(2x,a1,14x,g14.6)' ) '*', ydata(i)
        if ( i .lt. ndata ) then
          write ( *, '(2x,a1, g14.6)' ) '*', tdata(i)
        end if
      end do
c
c  Sample the spline.
c
      write ( *, * ) 'DEBUG: TLO = ', tlo
      write ( *, * ) 'DEBUG: THI = ', thi

      alo = tlo - 1.0D+00
      ahi = thi + 1.0D+00

      call r8vec_uniform_01 ( n_test, seed, t_test )

      t_test(1:n_test) = ( 1.0D+00 - t_test(1:n_test) ) * alo 
     &                 +             t_test(1:n_test)   * ahi

      call r8vec_sort_bubble_a ( n_test, t_test )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '     T     Y(interp)    Y(exact)'
      write ( *, '(a)' ) ' '

      j = 0
      write ( *, '(2x,a1,14x,g14.6)' ) '*', ydata(j+1)
      j = j + 1

      do i = 1, n_test

        tval = t_test(i)

        call spline_constant_val ( ndata, tdata, ydata, tval, yval )

        if ( j .le. ndata - 1 ) then
          do while ( tdata(j) .le. tval )
            fval = frunge ( tdata(j) )
            write ( *, '(2x,a1,g14.6,14x,g14.6)' ) '*', tdata(j), fval
            write ( *, '(2x,a1,14x,g14.6)' ) '*', ydata(j+1)
            j = j + 1
            if ( ndata .le. j ) then
              exit
            end if
          end do
        end if

        fval = frunge ( tval )

        write ( *, '(2x,a1,3g14.6)' ) ' ', tval, yval, fval

      end do

      return
      end
      subroutine test15 ( )

c*********************************************************************72
c
cc TEST15 tests SPLINE_CUBIC_SET and SPLINE_CUBIC_VAL.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    17 June 2006
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n
      parameter ( n = 11 )

      double precision frunge
      double precision fprunge
      double precision fpprunge
      integer i
      integer ibcbeg
      integer ibcend
      integer j
      integer jhi
      integer k
      double precision t(n)
      double precision tval
      double precision y(n)
      double precision ybcbeg
      double precision ybcend
      double precision ypp(n)
      double precision yppval
      double precision ypval
      double precision yval

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST15'
      write ( *, '(a)' ) '  SPLINE_CUBIC_SET sets up a cubic spline;'
      write ( *, '(a)' ) '  SPLINE_CUBIC_VAL evaluates it.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Runge''s function, evenly spaced knots.'

      do i = 1, n

        t(i) = ( dble ( n - i     ) * (-1.0D+00)   
     &         + dble (     i - 1 ) * (+1.0D+00) ) 
     &         / dble ( n     - 1 )

        y(i) =  frunge ( t(i) )

      end do

      write ( *, '(a)'    ) ' '
      write ( *, '(a)'    ) '  The data to be interpolated:'
      write ( *, '(a)'    ) ' '
      write ( *, '(a,i8)' ) '  Number of data values = ', n
      write ( *, '(a)'    ) ' '
      write ( *, '(a)'    ) '       T             Y'
      write ( *, '(a)'    ) ' '

      do i = 1, n
        write ( *, '(2x,g14.6,2x,g14.6)' ) t(i), y(i)
      end do
c
c  Try boundary condition types 0, 1 and 2.
c
      do k = 0, 3

        if ( k .eq. 0 ) then

          ibcbeg = 0
          ybcbeg = 0.0D+00

          ibcend = 0
          ybcend = 0.0D+00

          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) '  Boundary condition 0 at both ends:'
          write ( *, '(a)' ) 
     &      '  Spline is quadratic in boundary intervals.'

        else if ( k .eq. 1 ) then

          ibcbeg = 1
          ybcbeg = fprunge ( t(1) )

          ibcend = 1
          ybcend = fprunge ( t(n) )

          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) '  Boundary condition 1 at both ends:'
          write ( *, '(a,g14.6)' ) '  Y''(left) =  ', ybcbeg
          write ( *, '(a,g14.6)' ) '  Y''(right) = ', ybcend

        else if ( k .eq. 2 ) then

          ibcbeg = 2
          ybcbeg = fpprunge ( t(1) )

          ibcend = 2
          ybcend = fpprunge ( t(n) )

          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) '  Boundary condition 2 at both ends:'
          write ( *, '(a,g14.6)' ) '  YP"(left) =  ', ybcbeg
          write ( *, '(a,g14.6)' ) '  YP"(right) = ', ybcend

        else if ( k .eq. 3 ) then

          ibcbeg = 2
          ybcbeg = 0.0D+00

          ibcend = 2
          ybcend = 0.0D+00

          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) '  "Natural" spline:'
          write ( *, '(a)' ) '  Boundary condition 2 at both ends:'
          write ( *, '(a,g14.6)' ) '  YP"(left) =  ', ybcbeg
          write ( *, '(a,g14.6)' ) '  YP"(right) = ', ybcend

        end if

        call spline_cubic_set ( n, t, y, ibcbeg, ybcbeg, ibcend, 
     &    ybcend, ypp )

        if ( k .eq. 3 ) then
          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) '       I      Y(I)          YPP(I)'
          write ( *, '(a)' ) ' '
          do i = 1, n
            write ( *, '(2x,i8,2g14.6)' ) i, y(i), ypp(i)
          end do
        end if

        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) '    SPLINE"(T)      F"(T):'
        write ( *, '(a)' ) ' '
        do i = 1, n
          write ( *, '(2x,2g14.6)' ) ypp(i), fpprunge(t(i))
        end do

        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) '      T            SPLINE(T)      F(T)'
        write ( *, '(a)' ) ' '

        do i = 0, n

          if ( i .eq. 0 ) then
            jhi = 1
          else if ( i .lt. n ) then
            jhi = 2
          else
            jhi = 2
          end if

          do j = 1, jhi

            if ( i .eq. 0 ) then
              tval = t(1) - 1.0D+00
            else if ( i .lt. n ) then
              tval = ( dble ( jhi - j + 1 ) * t(i)     
     &               + dble (       j - 1 ) * t(i+1) ) 
     &               / dble ( jhi         )
            else
              if ( j .eq. 1 ) then
                tval = t(n)
              else
                tval = t(n) + 1.0D+00
              end if
            end if

            call spline_cubic_val ( n, t, y, ypp, tval, yval, ypval, 
     &        yppval )

            write ( *, '(2x,3g14.6)' ) tval, yval, frunge ( tval )

          end do
        end do

      end do

      return
      end
      subroutine test16 ( )

c*********************************************************************72
c
cc TEST16 tests SPLINE_CUBIC_SET and SPLINE_CUBIC_VAL2.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    17 June 2006
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n
      parameter ( n = 11 )

      double precision frunge
      double precision fprunge
      double precision fpprunge
      integer i
      integer ibcbeg
      integer ibcend
      integer j
      integer jhi
      integer k
      integer left
      integer left_in
      double precision t(n)
      double precision tval
      double precision y(n)
      double precision ybcbeg
      double precision ybcend
      double precision ypp(n)
      double precision yppval
      double precision ypval
      double precision yval

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST16'
      write ( *, '(a)' ) '  SPLINE_CUBIC_SET sets up a cubic spline;'
      write ( *, '(a)' ) '  SPLINE_CUBIC_VAL2 evaluates it.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Runge''s function, evenly spaced knots.'
      do i = 1, n

        t(i) = ( dble ( n - i     ) * (-1.0D+00)   
     &         + dble (     i - 1 ) * (+1.0D+00) ) 
     &         / dble ( n     - 1 )

        y(i) =  frunge ( t(i) )

      end do

      write ( *, '(a)'    ) ' '
      write ( *, '(a)'    ) '  The data to be interpolated:'
      write ( *, '(a)'    ) ' '
      write ( *, '(a,i8)' ) '  Number of data values = ', n
      write ( *, '(a)'    ) ' '
      write ( *, '(a)'    ) '       T             Y'
      write ( *, '(a)'    ) ' '

      do i = 1, n
        write ( *, '(2x,2g14.6)' ) t(i), y(i)
      end do
c
c  Try boundary condition types 0, 1 and 2.
c
      do k = 0, 2

        if ( k .eq. 0 ) then

          ibcbeg = 0
          ybcbeg = 0.0D+00

          ibcend = 0
          ybcend = 0.0D+00

          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) '  Boundary condition 0 at both ends:'
          write ( *, '(a)' ) 
     &      '  Spline is quadratic in boundary intervals.'

        else if ( k .eq. 1 ) then

          ibcbeg = 1
          ybcbeg = fprunge ( t(1) )

          ibcend = 1
          ybcend = fprunge ( t(n) )

          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) '  Boundary condition 1 at both ends:'
          write ( *, '(a,g14.6)' ) '  Y''(left) =  ', ybcbeg
          write ( *, '(a,g14.6)' ) '  Y''(right) = ', ybcend

        else if ( k .eq. 2 ) then

          ibcbeg = 2
          ybcbeg = fpprunge ( t(1) )

          ibcend = 2
          ybcend = fpprunge ( t(n) )

          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) '  Boundary condition 2 at both ends:'
          write ( *, '(a,g14.6)' ) '  YP"(left) =  ', ybcbeg
          write ( *, '(a,g14.6)' ) '  YP"(right) = ', ybcend

        end if

        call spline_cubic_set ( n, t, y, ibcbeg, ybcbeg, ibcend, 
     &    ybcend, ypp )

        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) '  SPLINE"(T)      F"(T)'
        write ( *, '(a)' ) ' '
        do i = 1, n
          write ( *, '(2x,2g14.6)' ) ypp(i), fpprunge(t(i))
        end do

        left = 0

        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 
     &    '      T             SPLINE(T)       F(T)       ' // 
     &    'LEFT_IN  LEFT_OUT'
        write ( *, '(a)' ) ' '

        do i = 0, n

          if ( i .eq. 0 ) then
            jhi = 1
          else if ( i .lt. n ) then
            jhi = 2
          else
            jhi = 2
          end if

          do j = 1, jhi

            if ( i .eq. 0 ) then
              tval = t(1) - 1.0D+00
            else if ( i .lt. n ) then
              tval = ( dble ( jhi - j + 1 ) * t(i)     
     &               + dble (       j - 1 ) * t(i+1) ) 
     &               / dble ( jhi         )
            else
              if ( j .eq. 1 ) then
                tval = t(n)
              else
                tval = t(n) + 1.0D+00
              end if
            end if

            left_in = left

            call spline_cubic_val2 ( n, t, y, ypp, left, tval, 
     &        yval, ypval, yppval )

            write ( *, '(2x,3g14.6,2i8)' ) tval, yval, frunge ( tval ), 
     &        left_in, left

          end do
        end do

      end do

      return
      end
      subroutine test17 ( )

c*********************************************************************72
c
cc TEST17 tests SPLINE_CUBIC_SET and SPLINE_CUBIC_VAL.
c
c  Discussion:
c
c    For boundary condition 0, the spline should come very close within
c    the interpolation interval.
c
c    For conditions 1 and 2, the spline should be essentially exactly equal
c    to the data, inside and outside the interpolation interval.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    17 June 2006
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n
      parameter ( n = 11 )

      double precision fcube
      double precision fpcube
      double precision fppcube
      integer i
      integer ibcbeg
      integer ibcend
      integer j
      integer jhi
      integer k
      double precision t(n)
      double precision tval
      double precision y(n)
      double precision ybcbeg
      double precision ybcend
      double precision ypp(n)
      double precision yppval
      double precision ypval
      double precision yval

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST17'
      write ( *, '(a)' ) '  SPLINE_CUBIC_SET sets up a cubic spline;'
      write ( *, '(a)' ) '  SPLINE_CUBIC_VAL evaluates it.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Cubic data, unevenly spaced knots.'

      do i = 1, n
        t(i) = ( dble ( i - 1 ) / dble ( n - 1 ) )**2
      end do

      do i = 1, n
        y(i) = fcube ( t(i) )
      end do

      write ( *, '(a)'    ) ' '
      write ( *, '(a)'    ) '  The data to be interpolated:'
      write ( *, '(a)'    ) ' '
      write ( *, '(a,i8)' ) '  Number of data values = ', n
      write ( *, '(a)'    ) ' '
      write ( *, '(a)'    ) '       T             Y'
      write ( *, '(a)'    ) ' '

      do i = 1, n
        write ( *, '(2x,g14.6,2x,g14.6)' ) t(i), y(i)
      end do
c
c  Try boundary condition types 0, 1 and 2.
c
      do k = 0, 2

        if ( k .eq. 0 ) then

          ibcbeg = 0
          ybcbeg = 0.0D+00

          ibcend = 0
          ybcend = 0.0D+00

          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) '  Boundary condition 0 at both ends:'
          write ( *, '(a)' ) 
     &      '  Spline is quadratic in boundary intervals.'

        else if ( k .eq. 1 ) then

          ibcbeg = 1
          ybcbeg = fpcube ( t(1) )

          ibcend = 1
          ybcend = fpcube ( t(n) )

          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) '  Boundary condition 1 at both ends:'
          write ( *, '(a,g14.6)' ) '  Y''(left) =  ', ybcbeg
          write ( *, '(a,g14.6)' ) '  Y''(right) = ', ybcend

        else if ( k .eq. 2 ) then

          ibcbeg = 2
          ybcbeg = fppcube ( t(1) )

          ibcend = 2
          ybcend = fppcube ( t(n) )

          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) '  Boundary condition 2 at both ends:'
          write ( *, '(a,g14.6)' ) '  YP"(left) =  ', ybcbeg
          write ( *, '(a,g14.6)' ) '  YP"(right) = ', ybcend

        end if

        call spline_cubic_set ( n, t, y, ibcbeg, ybcbeg, ibcend, 
     &    ybcend, ypp )

        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) '       SPLINE"(T)    F"(T):'
        write ( *, '(a)' ) ' '
        do i = 1, n
          write ( *, '(2x,2g14.6)' ) ypp(i), fppcube(t(i))
        end do

        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) '       T           SPLINE(T)      F(T)'
        write ( *, '(a)' ) ' '

        do i = 0, n

          if ( i .eq. 0 ) then
            jhi = 1
          else if ( i .lt. n ) then
            jhi = 2
          else
            jhi = 2
          end if

          do j = 1, jhi

            if ( i .eq. 0 ) then
              tval = t(1) - 1.0D+00
            else if ( i .lt. n ) then
              tval = ( dble ( jhi - j + 1 ) * t(i)     
     &               + dble (       j - 1 ) * t(i+1) ) 
     &               / dble ( jhi         )
            else
              if ( j .eq. 1 ) then
                tval = t(n)
              else
                tval = t(n) + 1.0D+00
              end if
            end if

            call spline_cubic_val ( n, t, y, ypp, tval, yval, ypval, 
     &        yppval )

            write ( *, '(2x,3g14.6)' ) tval, yval, fcube ( tval )

          end do
        end do

      end do

      return
      end
      subroutine test18 ( )

c*********************************************************************72
c
cc TEST18 tests SPLINE_CUBIC_SET and SPLINE_CUBIC_VAL.
c
c  Discussion:
c
c    For boundary condition 0, the spline should come very close within
c    the interpolation interval.
c
c    For conditions 1 and 2, the spline should be essentially exactly equal
c    to the data, inside and outside the interpolation interval.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    17 June 2006
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n
      parameter ( n = 11 )

      double precision fcube
      double precision fpcube
      double precision fppcube
      integer i
      integer ibcbeg
      integer ibcend
      integer j
      integer jhi
      integer k
      double precision t(n)
      double precision tval
      double precision y(n)
      double precision ybcbeg
      double precision ybcend
      double precision ypp(n)
      double precision yppval
      double precision ypval
      double precision yval

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST18'
      write ( *, '(a)' ) '  SPLINE_CUBIC_SET sets up a cubic spline;'
      write ( *, '(a)' ) '  SPLINE_CUBIC_VAL evaluates it.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Cubic data, evenly spaced knots.'
      do i = 1, n
        t(i) = dble ( i - 1 ) / dble ( n - 1 )
        y(i) =  fcube ( t(i) )
      end do

      write ( *, '(a)'    ) ' '
      write ( *, '(a)'    ) '  The data to be interpolated:'
      write ( *, '(a)'    ) ' '
      write ( *, '(a,i8)' ) '  Number of data values = ', n
      write ( *, '(a)'    ) ' '
      write ( *, '(a)'    ) '       T             Y'
      write ( *, '(a)'    ) ' '

      do i = 1, n
        write ( *, '(2x,g14.6,2x,g14.6)' ) t(i), y(i)
      end do
c
c  Try boundary condition types 0, 1 and 2.
c
      do k = 0, 2

        if ( k .eq. 0 ) then

          ibcbeg = 0
          ybcbeg = 0.0D+00

          ibcend = 0
          ybcend = 0.0D+00

          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) '  Boundary condition 0 at both ends:'
          write ( *, '(a)' ) 
     &      '  Spline is quadratic in boundary intervals.'

        else if ( k .eq. 1 ) then

          ibcbeg = 1
          ybcbeg = fpcube ( t(1) )

          ibcend = 1
          ybcend = fpcube ( t(n) )

          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) '  Boundary condition 1 at both ends:'
          write ( *, '(a,g14.6)' ) '  Y''(left) =  ', ybcbeg
          write ( *, '(a,g14.6)' ) '  Y''(right) = ', ybcend

        else if ( k .eq. 2 ) then

          ibcbeg = 2
          ybcbeg = fppcube ( t(1) )

          ibcend = 2
          ybcend = fppcube ( t(n) )

          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) '  Boundary condition 2 at both ends:'
          write ( *, '(a,g14.6)' ) '  YP"(left) =  ', ybcbeg
          write ( *, '(a,g14.6)' ) '  YP"(right) = ', ybcend

        end if

        call spline_cubic_set ( n, t, y, ibcbeg, ybcbeg, ibcend, 
     &    ybcend, ypp )

        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) '   SPLINE"(T)      F"(T):'
        write ( *, '(a)' ) ' '
        do i = 1, n
          write ( *, '(2x,2g14.6)' ) ypp(i), fppcube(t(i))
        end do

        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) '        T      SPLINE(T)    F(T)'
        write ( *, '(a)' ) ' '

        do i = 0, n

          if ( i .eq. 0 ) then
            jhi = 1
          else if ( i .lt. n ) then
            jhi = 2
          else
            jhi = 2
          end if

          do j = 1, jhi

            if ( i .eq. 0 ) then
              tval = t(1) - 1.0D+00
            else if ( i .lt. n ) then
              tval = ( dble ( jhi - j + 1 ) * t(i)     
     &               + dble (       j - 1 ) * t(i+1) ) 
     &               / dble ( jhi         )
            else
              if ( j .eq. 1 ) then
                tval = t(n)
              else
                tval = t(n) + 1.0D+00
              end if
            end if

            call spline_cubic_val ( n, t, y, ypp, tval, yval, 
     &        ypval, yppval )

            write ( *, '(2x,f10.4)' ) tval
            write ( *, '(2x,10x,2f10.4)' )   yval, fcube ( tval )
            write ( *, '(2x,10x,2f10.4)' )   ypval, fpcube ( tval )
            write ( *, '(2x,10x,2f10.4)' )   yppval, fppcube ( tval )

          end do
        end do

      end do

      return
      end
      subroutine test19 ( )

c*********************************************************************72
c
cc TEST19 tests SPLINE_CUBIC_SET and SPLINE_CUBIC_VAL.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    17 June 2006
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n
      parameter ( n = 2 )

      double precision fcube
      double precision fpcube
      double precision fppcube
      integer i
      integer ibcbeg
      integer ibcend
      integer j
      integer jhi
      integer k1
      integer k2
      double precision t(n)
      double precision tval
      double precision y(n)
      double precision ybcbeg
      double precision ybcend
      double precision ypp(n)
      double precision yppval
      double precision ypval
      double precision yval

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST19'
      write ( *, '(a)' ) '  SPLINE_CUBIC_SET sets up a cubic spline;'
      write ( *, '(a)' ) '  SPLINE_CUBIC_VAL evaluates it.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Cubic data, evenly spaced knots.'
      write ( *, '(a)' ) '  ONLY TWO KNOTS!'

      do i = 1, n
        t(i) = dble ( i - 1 ) / dble ( n - 1 )
        y(i) = fcube ( t(i) )
      end do

      write ( *, '(a)'    ) ' '
      write ( *, '(a)'    ) '  The data to be interpolated:'
      write ( *, '(a)'    ) ' '
      write ( *, '(a,i8)' ) '  Number of data values = ', n
      write ( *, '(a)'    ) ' '
      write ( *, '(a)'    ) '       T             Y'
      write ( *, '(a)'    ) ' '

      do i = 1, n
        write ( *, '(2x,g14.6,2x,g14.6)' ) t(i), y(i)
      end do
c
c  Try all 9 pairs of boundary condition types 0, 1 and 2.
c
      do k1 = 0, 2

        if ( k1 .eq. 0 ) then

          ibcbeg = 0
          ybcbeg = 0.0D+00

          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) '  Boundary condition 0 at left end.'

        else if ( k1 .eq. 1 ) then

          ibcbeg = 1
          ybcbeg = fpcube ( t(1) )

          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) '  Boundary condition 1 at left end.'
          write ( *, '(a,g14.6)' ) '  Y''(left) =  ', ybcbeg

        else if ( k1 .eq. 2 ) then

          ibcbeg = 2
          ybcbeg = fppcube ( t(1) )

          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) '  Boundary condition 2 at left ends:'
          write ( *, '(a,g14.6)' ) '  YP"(left) =  ', ybcbeg

        end if

        do k2 = 0, 2

          if ( k2 .eq. 0 ) then

            ibcend = 0
            ybcend = 0.0D+00

            write ( *, '(a)' ) ' '
            write ( *, '(a)' ) '  Boundary condition 0 at right end.'

          else if ( k2 .eq. 1 ) then

            ibcend = 1
            ybcend = fpcube ( t(n) )

            write ( *, '(a)' ) ' '
            write ( *, '(a)' ) '  Boundary condition 1 at right end.'
            write ( *, '(a,g14.6)' ) '  Y''(right) = ', ybcend

          else if ( k2 .eq. 2 ) then

            ibcend = 2
            ybcend = fppcube ( t(n) )

            write ( *, '(a)' ) ' '
            write ( *, '(a)' ) '  Boundary condition 2 at right end.'
            write ( *, '(a,g14.6)' ) '  YP"(right) = ', ybcend

          end if

          call spline_cubic_set ( n, t, y, ibcbeg, ybcbeg, ibcend, 
     &      ybcend, ypp )

          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) '  SPLINE"(T), F"(T):'
          write ( *, '(a)' ) ' '
          do i = 1, n
            write ( *, '(2x,2g14.6)' ) ypp(i), fppcube(t(i))
          end do

          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) '  T, SPLINE(T), F(T)'
          write ( *, '(a)' ) ' '

          do i = 0, n

            if ( i .eq. 0 ) then
              jhi = 1
            else if ( i .lt. n ) then
              jhi = 2
            else
              jhi = 2
            end if

            do j = 1, jhi

              if ( i .eq. 0 ) then
                tval = t(1) - 1.0D+00
              else if ( i .lt. n ) then
                tval = ( dble ( jhi - j + 1 ) * t(i)     
     &                 + dble (       j - 1 ) * t(i+1) ) 
     &                 / dble ( jhi         )
              else
                if ( j .eq. 1 ) then
                  tval = t(n)
                else
                  tval = t(n) + 1.0D+00
                end if
              end if

              call spline_cubic_val ( n, t, y, ypp, tval, yval, ypval, 
     &          yppval )

              write ( *, '(2x,f10.4)' ) tval
              write ( *, '(2x,10x,2f10.4)' )   yval, fcube ( tval )
              write ( *, '(2x,10x,2f10.4)' )   ypval, fpcube ( tval )
              write ( *, '(2x,10x,2f10.4)' )   yppval, fppcube ( tval )

            end do
          end do

        end do

      end do

      return
      end
      subroutine test20 ( )

c*********************************************************************72
c
cc TEST20 tests SPLINE_HERMITE_SET and SPLINE_HERMITE_VAL.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    17 June 2006
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer ndata
      parameter ( ndata = 4 )

      double precision c(4,ndata)
      double precision fpval
      double precision fval
      integer i
      integer j
      integer jhi
      character mark
      double precision pi
      parameter ( pi = 3.141592653589793D+00 )
      double precision tdata(ndata)
      double precision tval
      double precision ydata(ndata)
      double precision ypdata(ndata)
      double precision ypval
      double precision yval

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST20'
      write ( *, '(a)' ) 
     &  '  SPLINE_HERMITE_SET sets up a Hermite spline;'
      write ( *, '(a)' ) '  SPLINE_HERMITE_VAL evaluates it.'
c
c  Set the data.
c
      do i = 1, ndata
        tdata(i) = ( dble ( ndata - i     ) *   0.0D+00          
     &             + dble (         i - 1 ) * ( 0.5D+00 * pi ) ) 
     &             / dble ( ndata     - 1 )
        ydata(i) = sin ( tdata(i) )
        ypdata(i) = cos ( tdata(i) )
      end do

      write ( *, '(a)'    ) ' '
      write ( *, '(a)'    ) '  The data to be interpolated:'
      write ( *, '(a)'    ) ' '
      write ( *, '(a,i8)' ) '  Number of data values = ', ndata
      write ( *, '(a)'    ) ' '
      write ( *, '(a)'    ) '       T             Y             Y'''
      write ( *, '(a)'    ) ' '
      do i = 1, ndata
        write ( *, '(2x,3g14.6)' ) tdata(i), ydata(i), ypdata(i)
      end do
c
c  Set up the spline.
c
      call spline_hermite_set ( ndata, tdata, ydata, ypdata, c )
c
c  Now evaluate the spline all over the place.
c
      write ( *, '(a)'    ) ' '
      write ( *, '(a)'    ) '       T        Y(hermite)     ' // 
     &  'Y(exact)      Y''(hermite)   Y''(exact)'
      write ( *, '(a)'    ) ' '

      do i = 1, ndata

        if ( i .eq. ndata ) then
          jhi = 0
        else
          jhi = 2
        end if

        do j = 0, jhi

          tval = dble ( 3 * ( i - 1 ) + j ) * ( 0.5D+00 * pi ) 
     &         / dble ( 3 * ( ndata - 1 ) )

          fval = sin ( tval )
          fpval = cos ( tval )

          call spline_hermite_val ( ndata, tdata, c, tval, yval, ypval )

          if (j .eq. 0 ) then
            mark = '*'
          else
            mark = ' '
          end if

          write ( *, '(2x,a1,5g14.6)' ) 
     &      mark, tval, yval, fval, ypval, fpval

        end do

      end do

      return
      end
      subroutine test205 ( )

c*********************************************************************72
c
cc TEST205 tests SPLINE_LINEAR_INT and SPLINE_LINEAR_INTSET.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    17 June 2006
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n
      parameter ( n = 4 )

      double precision a
      double precision b
      double precision data_x(n)
      double precision data_y(n)
      integer i
      double precision int_x(n+1)
      double precision int_v(n)
      double precision value

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST205'
      write ( *, '(a)' ) 
     &  '  SPLINE_LINEAR_INTSET is given some interval endpoints,'
      write ( *, '(a)' ) '  and a value associated with each interval.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 
     &  '  It determines a linear spline, with breakpoints'
      write ( *, '(a)' ) 
     &  '  at the centers of each interval, whose integral'
      write ( *, '(a)' ) 
     &  '  over each interval is equal to the given value.'

      int_x(1:n+1) = (/ 0.0D+00, 1.0D+00, 4.0D+00, 5.0D+00, 10.0D+00 /)
      int_v(1:n) = (/ 10.0D+00, 2.0D+00, 8.0D+00, 27.5D+00 /)

      call r8vec_print ( n+1, int_x, '  The interval end points:' )

      call r8vec_print ( n, int_v, 
     &  '  The desired interval integral values:' )

      call spline_linear_intset ( n, int_x, int_v, data_x, data_y )

      call r8vec_print ( n, data_x, '  The spline break points:' )
      call r8vec_print ( n, data_y, '  The spline data values: ' )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 
     &  '  As a check, call SPLINE_LINEAR_INT to compute'
      write ( *, '(a)' ) 
     &  '  the integral of the spline over each interval,'
      write ( *, '(a)' ) '  and compare to the desired value.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '     A       B    Desired      Computed'
      write ( *, '(a)' ) ' '

      do i = 1, n
        a = int_x(i)
        b = int_x(i+1)
        call spline_linear_int ( n, data_x, data_y, a, b, value )
        write ( *, '(2x,2f8.2,2g14.6)' ) a, b, int_v(i), value
      end do

      return
      end
      subroutine test21 ( )

c*********************************************************************72
c
cc TEST21 tests SPLINE_LINEAR_VAL.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    17 June 2006
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n
      parameter ( n = 11 )

      double precision frunge
      double precision fval
      integer i
      integer j
      integer jhi
      double precision t(n)
      double precision tval
      double precision y(n)
      double precision ypval
      double precision yval

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST21'
      write ( *, '(a)' ) 
     &  '  SPLINE_LINEAR_VAL evaluates a linear spline.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Runge''s function, evenly spaced knots.'

      do i = 1, n

        t(i) = ( dble ( n - i     ) * (-1.0D+00)   
     &         + dble (     i - 1 ) * (+1.0D+00) ) 
     &         / dble ( n     - 1 )

        y(i) =  frunge ( t(i) )

      end do

      write ( *, '(a)'    ) ' '
      write ( *, '(a)'    ) '  The data to be interpolated:'
      write ( *, '(a)'    ) ' '
      write ( *, '(a,i8)' ) '  Number of data values = ', n
      write ( *, '(a)'    ) ' '
      write ( *, '(a)'    ) '       T             Y'
      write ( *, '(a)'    ) ' '
      do i = 1, n
        write ( *, '(2x,2g14.6)' ) t(i), y(i)
      end do

      write ( *, '(a)'    ) ' '
      write ( *, '(a)'    ) '  Interpolation:'
      write ( *, '(a)'    ) ' '
      write ( *, '(a)'    ) '       T             Y            Yexact'
      write ( *, '(a)'    ) ' '

      do i = 0, n

        if ( i .eq. 0 ) then
          jhi = 1
        else if ( i .lt. n ) then
          jhi = 2
        else
          jhi = 2
        end if

        do j = 1, jhi

          if ( i .eq. 0 ) then
            tval = t(1) - 1.0D+00
          else if ( i .lt. n ) then
            tval = ( dble ( jhi - j + 1 ) * t(i)     
     &             + dble (       j - 1 ) * t(i+1) ) 
     &             / dble ( jhi         )
          else
            if ( j .eq. 1 ) then
              tval = t(n)
            else
              tval = t(n) + 1.0D+00
            end if
          end if

          call spline_linear_val ( n, t, y, tval, yval, ypval )

          fval = frunge ( tval )

          write ( *, '(2x,3g14.6)' ) tval, yval, fval

        end do

      end do

      return
      end
      subroutine test215 ( )

c*********************************************************************72
c
cc TEST215 tests SPLINE_LINEAR_INT.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    17 June 2006
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n
      parameter ( n = 3 )

      double precision a
      double precision b
      integer i
      double precision int_val
      double precision t(n)
      double precision y(n)

      save t
      save y

      data t /
     & 2.0D+00, 4.5D+00, 7.5D+00 /
      data y /
     & 3.0D+00, 3.75D+00, 5.5D+00 /

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST215'
      write ( *, '(a)' ) '  SPLINE_LINEAR_INT computes the integral '
      write ( *, '(a)' ) '  of a linear spline.'
      write ( *, '(a)' ) ' '

      write ( *, '(a)'    ) ' '
      write ( *, '(a)'    ) '  The data to be interpolated:'
      write ( *, '(a)'    ) ' '
      write ( *, '(a,i8)' ) '  Number of data values = ', n
      write ( *, '(a)'    ) ' '
      write ( *, '(a)'    ) '       T             Y'
      write ( *, '(a)'    ) ' '
      do i = 1, n
        write ( *, '(2x,2g14.6)' ) t(i), y(i)
      end do

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '    A             B           Integral'
      write ( *, '(a)' ) ' '

      do i = 1, 5

        if ( i .eq. 1 ) then
          a = 0.0D+00
          b = 4.0D+00
        else if ( i .eq. 2 ) then
          a = 4.0D+00
          b = 5.0D+00
        else if ( i .eq. 3 ) then
          a = 5.0D+00
          b = 10.0D+00
        else if ( i .eq. 4 ) then
          a = 0.0D+00
          b = 10.0D+00
        else
          a = 10.0D+00
          b = 0.0D+00
        end if

        call spline_linear_int ( n, t, y, a, b, int_val )

        write ( *, '(2x,3g14.6)' ) a, b, int_val

      end do

      return
      end
      subroutine test22 ( )

c*********************************************************************72
c
cc TEST22 tests SPLINE_OVERHAUSER_UNI_VAL.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    17 June 2006
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer ndata
      parameter ( ndata = 11 )

      integer i
      integer j
      integer jhi
      character mark
      integer nsample
      parameter ( nsample = 4 )
      double precision pi
      parameter ( pi = 3.141592653589793D+00 )
      double precision tdata(ndata)
      double precision thi
      double precision tlo
      double precision tval
      double precision ydata(ndata)
      double precision yval

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST22'
      write ( *, '(a)' ) '  SPLINE_OVERHAUSER_UNI_VAL evaluates the'
      write ( *, '(a)' ) '    uniform Overhauser spline.'

      do i = 1, ndata
        tdata(i) = dble ( i - 1 )
        ydata(i) = sin ( 2.0D+00 * pi * tdata(i) / dble ( ndata - 1 ) )
      end do

      write ( *, '(a)'    ) ' '
      write ( *, '(a)'    ) '  The data to be interpolated:'
      write ( *, '(a)'    ) ' '
      write ( *, '(a,i8)' ) '  Number of data values = ', ndata
      write ( *, '(a)'    ) ' '
      write ( *, '(a)'    ) '         T             Y'
      write ( *, '(a)'    ) ' '
      do i = 1, ndata
        write ( *, '(2x,2g14.6)' ) tdata(i), ydata(i)
      end do

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '      T, Spline(T)'
      write ( *, '(a)' ) ' '

      do i = 0, ndata

        if ( i .eq. 0 ) then
          tlo = tdata(1) - 0.5D+00 * ( tdata(2) - tdata(1) )
          thi = tdata(1)
        else if ( i .lt. ndata ) then
          tlo = tdata(i)
          thi = tdata(i+1)
        else if ( ndata .le. i ) then
          tlo = tdata(ndata)
          thi = tdata(ndata) 
     &      + 0.5D+00 * ( tdata(ndata) - tdata(ndata-1) )
        end if

        if ( i .lt. ndata ) then
          jhi = nsample - 1
        else
          jhi = nsample
        end if

        do j = 0, jhi

          tval = ( dble ( nsample - j ) * tlo   
     &           + dble (           j ) * thi ) 
     &           / dble ( nsample     )

          call spline_overhauser_uni_val ( ndata, tdata, ydata, tval, 
     &      yval )

          if ( 0 .lt. i .and. j .eq. 0 ) then
            mark = '*'
          else
            mark = ' '
          end if

          write ( *, '(2x,a1,2g14.6)' ) mark, tval, yval

        end do

      end do

      return
      end
      subroutine test225 ( )

c*********************************************************************72
c
cc TEST225 tests SPLINE_OVERHAUSER_NONUNI_VAL.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    17 June 2006
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer ndata
      parameter ( ndata = 11 )

      integer i
      integer j
      integer jhi
      character mark
      integer nsample
      parameter ( nsample = 4 )
      double precision pi
      parameter ( pi = 3.141592653589793D+00 )
      double precision tdata(ndata)
      double precision thi
      double precision tlo
      double precision tval
      double precision ydata(ndata)
      double precision yval

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST225'
      write ( *, '(a)' ) '  SPLINE_OVERHAUSER_NONUNI_VAL evaluates the'
      write ( *, '(a)' ) '    nonuniform Overhauser spline.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  In this initial draft of a test, we simply'
      write ( *, '(a)' ) '  use uniform nodes.'

      do i = 1, ndata
        tdata(i) = dble ( i - 1 )
        ydata(i) = sin ( 2.0D+00 * pi * tdata(i) / dble ( ndata - 1 ) )
      end do

      write ( *, '(a)'    ) ' '
      write ( *, '(a)'    ) '  The data to be interpolated:'
      write ( *, '(a)'    ) ' '
      write ( *, '(a,i8)' ) '  Number of data values = ', ndata
      write ( *, '(a)'    ) ' '
      write ( *, '(a)'    ) '         T             Y'
      write ( *, '(a)'    ) ' '
      do i = 1, ndata
        write ( *, '(2x,2g14.6)' ) tdata(i), ydata(i)
      end do

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '      T, Spline(T)'
      write ( *, '(a)' ) ' '

      do i = 0, ndata

        if ( i .eq. 0 ) then
          tlo = tdata(1) - 0.5D+00 * ( tdata(2) - tdata(1) )
          thi = tdata(1)
        else if ( i .lt. ndata ) then
          tlo = tdata(i)
          thi = tdata(i+1)
        else if ( ndata .le. i ) then
          tlo = tdata(ndata)
          thi = tdata(ndata) 
     &      + 0.5D+00 * ( tdata(ndata) - tdata(ndata-1) )
        end if

        if ( i .lt. ndata ) then
          jhi = nsample - 1
        else
          jhi = nsample
        end if

        do j = 0, jhi

          tval = ( dble ( nsample - j ) * tlo   
     &           + dble (           j ) * thi ) 
     &           / dble ( nsample     )

          call spline_overhauser_nonuni_val ( ndata, tdata, ydata, 
     &      tval, yval )

          if ( 0 .lt. i .and. j .eq. 0 ) then
            mark = '*'
          else
            mark = ' '
          end if

          write ( *, '(2x,a1,2g14.6)' ) mark, tval, yval

        end do

      end do

      return
      end
      subroutine test23 ( )

c*********************************************************************72
c
cc TEST23 tests SPLINE_OVERHAUSER_VAL.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    17 June 2006
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer ndata
      parameter ( ndata = 4 )
      integer dim_num
      parameter ( dim_num = 2 )

      integer i
      double precision tdata(ndata)
      double precision tval
      double precision ydata(dim_num,ndata)
      double precision yval(dim_num)

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST23'
      write ( *, '(a)' ) '  SPLINE_OVERHAUSER_VAL evaluates the'
      write ( *, '(a)' ) '    Overhauser spline.'
c
c  Set the data.
c
      tdata(1) = 1.0D+00
      ydata(1,1) =   0.0D+00
      ydata(2,1) =   0.0D+00

      tdata(2) = 2.0D+00
      ydata(1,2) =   1.0D+00
      ydata(2,2) =   1.0D+00

      tdata(3) = 3.0D+00
      ydata(1,3) =   2.0D+00
      ydata(2,3) = - 1.0D+00

      tdata(4) = 4.0D+00
      ydata(1,4) =   3.0D+00
      ydata(2,4) =   0.0D+00

      write ( *, '(a)'    ) ' '
      write ( *, '(a)'    ) '  The data to be interpolated:'
      write ( *, '(a)'    ) ' '
      write ( *, '(a,i8)' ) '  Number of data values = ', ndata
      write ( *, '(a)'    ) ' '
      write ( *, '(a)'    ) '       T             Y'
      write ( *, '(a)'    ) ' '

      do i = 1, ndata
        write ( *, '(2x,3g14.6)' ) tdata(i), ydata(1:dim_num,i)
      end do
c
c  Now evaluate the spline all over the place.
c
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  T, Spline value'
      write ( *, '(a)' ) ' '

      do i = 0, 6 * ndata + 3

        tval = dble ( i ) / 6.0D+00
        call spline_overhauser_val ( dim_num, ndata, tdata, ydata, 
     &    tval, yval )
        write ( *, '(2x,3g14.6)' ) tval, yval(1:dim_num)

      end do

      return
      end
      subroutine test235 ( )

c*********************************************************************72
c
cc TEST235 tests SPLINE_PCHIP_SET and SPLINE_PCHIP_VAL.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    14 August 2005
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n
      parameter ( n = 21 )
      integer ne
      parameter ( ne = 101 )

      double precision d(n)
      double precision diff
      double precision f(n)
      double precision fd(ne)
      double precision fe(ne)
      integer i
      double precision frunge
      external frunge
      double precision x(n)
      double precision xe(ne)

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST235'
      write ( *, '(a)' ) 
     &  '  SPLINE_PCHIP_SET carries out piecewise cubic '
      write ( *, '(a)' ) '    Hermite interpolation.'
      write ( *, '(a)' ) '  SPLINE_PCHIP_VAL evaluates the interpolant.'
      write ( *, '(a)' ) ' '
c
c  Compute Runge's function at N points in [-1,1].
c
      do i = 1, n
        x(i) = -1.0D+00 + dble ( i - 1 ) / 10.0D+00
        f(i) = frunge ( x(i) )
      end do
c
c  SPLINE_PCHIP_SET takes the data in X and F, and constructs a table in D
c  that defines the interpolant.
c
      call spline_pchip_set ( n, x, f, d )
c
c  Evaluate the interpolant and derivative at NE points from -1 to 0.
c
      do i = 1, ne
        xe(i) = -1.0D+00 + dble ( i - 1 ) / dble ( ne - 1 )
      end do

      call spline_pchip_val ( n, x, f, d, ne, xe, fe )
c
c  Print the table of X, F(exact) and F(interpolated)
c
      do i = 1, ne
        diff = fe(i) - frunge ( xe(i) )
        write ( *, '(2x,f8.4,2x,f10.6,2x,f10.6,2x,g14.6)' ) 
     &    xe(i), frunge ( xe(i) ), fe(i), diff
      end do

      return
      end
      subroutine test24 ( )

c*********************************************************************72
c
cc TEST24 tests SPLINE_QUADRATIC_VAL.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    17 June 2006
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n
      parameter ( n = 11 )

      double precision frunge
      double precision fval
      integer i
      integer j
      integer jhi
      double precision t(n)
      double precision tval
      double precision y(n)
      double precision ypval
      double precision yval

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST24'
      write ( *, '(a)' ) '  SPLINE_QUADRATIC_VAL evaluates a '
      write ( *, '(a)' ) '    quadratic spline.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Runge''s function, evenly spaced knots.'

      do i = 1, n

        t(i) = ( dble ( n - i     ) * (-1.0D+00)   
     &         + dble (     i - 1 ) * (+1.0D+00) ) 
     &         / dble ( n     - 1 )

        y(i) =  frunge ( t(i) )

      end do

      write ( *, '(a)'    ) ' '
      write ( *, '(a)'    ) '  The data to be interpolated:'
      write ( *, '(a)'    ) ' '
      write ( *, '(a,i8)' ) '  Number of data values = ', n
      write ( *, '(a)'    ) ' '
      write ( *, '(a)'    ) '         T             Y'
      write ( *, '(a)'    ) ' '
      do i = 1, n
        write ( *, '(2x,2g14.6)' ) t(i), y(i)
      end do

      write ( *, '(a)'    ) ' '
      write ( *, '(a)'    ) '  Interpolated values'
      write ( *, '(a)'    ) ' '
      write ( *, '(a)'    ) '       T             Y           Y(exact)'
      write ( *, '(a)'    ) ' '

      do i = 0, n

        if ( i .eq. 0 ) then
          jhi = 1
        else if ( i .lt. n ) then
          jhi = 2
        else
          jhi = 2
        end if

        do j = 1, jhi

          if ( i .eq. 0 ) then
            tval = t(1) - 1.0D+00
          else if ( i .lt. n ) then
            tval = ( dble ( jhi - j + 1 ) * t(i)     
     &             + dble (       j - 1 ) * t(i+1) ) 
     &             / dble ( jhi         )
          else
            if ( j .eq. 1 ) then
              tval = t(n)
            else
              tval = t(n) + 1.0D+00
            end if
          end if

          call spline_quadratic_val ( n, t, y, tval, yval, ypval )

          fval = frunge ( tval )

          write ( *, '(2x,3g14.6)' ) tval, yval, fval

        end do

      end do

      return
      end
      function fcube ( x )

c*********************************************************************72
c
cc FCUBE evaluates a cubic function.
c
c  Discussion:
c
c    Y(X) = ( ( X + 2 ) * X + 3 ) * X + 4
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    10 February 2004
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, real X, the argument.
c
c    Output, real FCUBE, the value of the function.
c
      implicit none

      double precision fcube
      double precision x

      fcube = ( ( (       1.0D+00 ) 
     &              * x + 2.0D+00 ) 
     &              * x + 3.0D+00 ) 
     &              * x + 4.0D+00

      return
      end
      function fpcube ( x )

c*********************************************************************72
c
cc FPCUBE evaluates the derivative of a cubic function.
c
c  Discussion:
c
c    Y(X) = ( ( X + 2 ) * X + 3 ) * X + 4
c
c  Modified:
c
c    10 February 2004
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, real X, the argument.
c
c    Output, real FPCUBE, the value of the derivative of the cubic function.
c
      implicit none

      double precision fpcube
      double precision x

      fpcube = ( 3.0D+00 * x + 4.0D+00 ) * x + 3.0D+00

      return
      end
      function fppcube ( x )

c*********************************************************************72
c
cc FPPCUBE evaluates the second derivative of a cubic function.
c
c  Discussion:
c
c    Y(X) = ( ( X + 2 ) * X + 3 ) * X + 4
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    10 February 2004
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, real X, the argument.
c
c    Output, real FPPCUBE, the second derivative of the cubic function.
c
      implicit none

      double precision fppcube
      double precision x

      fppcube = 6.0D+00 * x + 4.0D+00

      return
      end
      function frunge ( x )

c*********************************************************************72
c
cc FRUNGE evaluates the Runge function.
c
c  Discussion:
c
c    Interpolation of the Runge function at evenly spaced points in [-1,1]
c    is a common test.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    24 January 2004
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, real X, the argument.
c
c    Output, real FRUNGE, the value of the function.
c
      implicit none

      double precision frunge
      double precision x

      frunge = 1.0D+00 / ( 1.0D+00 + 25.0D+00 * x * x )

      return
      end
      function fprunge ( x )

c*********************************************************************72
c
cc FPRUNGE evaluates the derivative of the Runge function.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    24 January 2004
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, real X, the argument.
c
c    Output, real FPRUNGE, the value of the derivative of the Runge function.
c
      implicit none

      double precision fprunge
      double precision x

      fprunge = - 50.0D+00 * x / ( 1.0D+00 + 25.0D+00 * x * x )**2

      return
      end
      function fpprunge ( x )

c*********************************************************************72
c
cc FPPRUNGE evaluates the second derivative of the Runge function.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    24 January 2004
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, real X, the argument.
c
c    Output, real FPPRUNGE, the value of the second derivative of
c    the Runge function.
c
      implicit none

      double precision fpprunge
      double precision x

      fpprunge = ( - 50.0D+00 + 3750.0D+00 * x * x ) 
     &  / ( 1.0D+00 + 25.0D+00 * x * x )**3

      return
      end
