      program main

c*********************************************************************72
c
cc MAIN is the main program for BERNSTEIN_PRB.
c
c  Discussion:
c
c    BERNSTEIN_PRB calls the BERNSTEIN test routines.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    10 February 2012
c
c  Author:
c
c    John Burkardt
c
      implicit none

      call timestamp ( )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'BERNSTEIN_PRB'
      write ( *, '(a)' ) '  FORTRAN77 version'
      write ( *, '(a)' ) '  Test the BERNSTEIN library.'

      call test01 ( )
      call test02 ( )
      call test03 ( )
      call test04 ( )
      call test05 ( )
c
c  Terminate.
c
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'BERNSTEIN_PRB'
      write ( *, '(a)' ) '  Normal end of execution.'

      write ( *, '(a)' ) ' '
      call timestamp ( )

      stop
      end
      subroutine test01 ( )

c*********************************************************************72
c
cc TEST01 tests BERNSTEIN_POLY and BERNSTEIN_POLY_VALUES.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    10 February 2012
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n_max
      parameter ( n_max = 20 )

      double precision b
      double precision bvec(0:n_max)
      integer k
      integer n
      integer n_data
      double precision x

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST01:'
      write ( *, '(a)' ) 
     &  '  BERNSTEIN_POLY evaluates the Bernstein polynomials'
      write ( *, '(a)' ) '  based on the interval [0,1].'
      write ( *, '(a)' ) 
     &  '  BERNSTEIN_POLY_VALUES returns some exact values.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 
     &  '     N     K     X       Exact         BP01(N,K)(X)'
      write ( *, '(a)' ) ' '

      n_data = 0

10    continue

        call bernstein_poly_values ( n_data, n, k, x, b )

        if ( n_data .eq. 0 ) then
          go to 20
        end if

        call bernstein_poly ( n, x, bvec )

        write ( *, '(2x,i4,2x,i4,2x,f7.4,g14.6,2x,g14.6)' ) 
     &    n, k, x, b, bvec(k)

      go to 10

20    continue

      return
      end
      subroutine test02 ( )

c*********************************************************************72
c
cc TEST02 tests BPAB.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    10 February 2012
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n
      parameter ( n = 10 )

      double precision a
      double precision b
      double precision bern(0:n)
      integer k
      double precision x

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST02'
      write ( *, '(a)' ) 
     &  '  BPAB evaluates Bernstein polynomials over an'
      write ( *, '(a)' ) '  arbitrary interval [A,B].'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Here, we demonstrate that '
      write ( *, '(a)' ) '    BPAB(N,K,A1,B1)(X1) = BPAB(N,K,A2,B2)(X2)'
      write ( *, '(a)' ) '  provided only that'
      write ( *, '(a)' ) '    (X1-A1)/(B1-A1) = (X2-A2)/(B2-A2).'

      x = 0.3D+00
      a = 0.0D+00
      b = 1.0D+00
      call bpab ( n, a, b, x, bern )
     
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '     N     K     A        B        ' // 
     &  'X       BPAB(N,K,A,B)(X)'
      write ( *, '(a)' ) ' ' 
      do k = 0, n
        write ( *, '(2x,i4,2x,i4,2x,f7.4,2x,f7.4,2x,f7.4,2x,g14.6)' ) 
     &    n, k, a, b, x, bern(k)
      end do
     
      x = 1.3D+00
      a = 1.0D+00
      b = 2.0D+00
      call bpab ( n, a, b, x, bern )
     
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '     N     K     A        B        X' // 
     &  '       BPAB(N,K,A,B)(X)'
      write ( *, '(a)' ) ' ' 
      do k = 0, n
        write ( *, '(2x,i4,2x,i4,2x,f7.4,2x,f7.4,2x,f7.4,2x,g14.6)' ) 
     &    n, k, a, b, x, bern(k)
      end do

      x = 2.6D+00
      a = 2.0D+00
      b = 4.0D+00
      call bpab ( n, a, b, x, bern )
     
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '     N     K     A        B        X' // 
     &  '       BPAB(N,K,A,B)(X)'
      write ( *, '(a)' ) ' '
     
      do k = 0, n
        write ( *, '(2x,i4,2x,i4,2x,f7.4,2x,f7.4,2x,f7.4,2x,g14.6)' ) 
     &    n, k, a, b, x, bern(k)
      end do

      return
      end
      subroutine test03 ( )

c*********************************************************************72
c
cc TEST03 tests the Partition-of-Unity property.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    10 February 2012
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n_max
      parameter ( n_max = 10 )

      double precision bvec(0:n_max)
      integer n
      integer n_data
      double precision r8_uniform_01
      integer seed
      double precision x

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST03:'
      write ( *, '(a)' ) 
     &  '  BERNSTEIN_POLY evaluates the Bernstein polynomials'
      write ( *, '(a)' ) '  based on the interval [0,1].'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 
     &  '  Here we test the partition of unity property.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 
     &  '     N     X          Sum ( 0 <= K <= N ) BP01(N,K)(X)'
      write ( *, '(a)' ) ' '

      seed = 123456789

      do n = 0, 10

        x = r8_uniform_01 ( seed )

        call bernstein_poly ( n, x, bvec )

        write ( *, '(2x,i4,2x,f7.4,2x,g14.6)' ) n, x, sum ( bvec(0:n) )

      end do

      return
      end
      subroutine test04 ( )

c*********************************************************************72
c
cc TEST04 tests BPAB_APPROX.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    10 February 2012
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer maxdata
      parameter ( maxdata = 20 )
      integer nval
      parameter ( nval = 501 )

      double precision a
      double precision b
      double precision error_max
      integer i
      integer ndata
      integer nsample
      double precision xdata(0:maxdata)
      double precision xval(nval)
      double precision ydata(0:maxdata)
      double precision yval(nval)

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST04'
      write ( *, '(a)' ) 
     &  '  BPAB_APPROX evaluates the Bernstein polynomial'
      write ( *, '(a)' ) '  approximant to a function F(X).'

      a = 1.0D+00
      b = 3.0D+00

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '     N      Max Error'
      write ( *, '(a)' ) ' '

      do ndata = 0, maxdata
c
c  Generate data values.
c
        do i = 0, ndata

          if ( ndata == 0 ) then
            xdata(i) = 0.5D+00 * ( a + b )
          else
            xdata(i) = ( dble ( ndata - i ) * a   
     &                 + dble (         i ) * b ) 
     &                 / dble ( ndata )
          end if

          ydata(i) = sin ( xdata(i) )

        end do
c
c  Compare the true function and the approximant.
c
        call r8vec_linspace ( nval, a, b, xval )

        error_max = 0.0D+00

        call bpab_approx ( ndata, a, b, ydata, nval, xval, yval )

        error_max = 0.0D+00
        do i = 1, nval
          error_max = 
     &      max ( error_max, abs ( yval(i) - sin ( xval(i) ) ) )
        end do

        write ( *, '(2x,i4,2x,g14.6)' ) ndata, error_max

      end do

      return
      end
      subroutine test05 ( )

c*********************************************************************72
c
cc TEST05 tests BERNSTEIN_MATRIX and BERNSTEIN_MATRIX_INVERSE.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    10 February 2012
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n_max
      parameter ( n_max = 15 )

      double precision a(n_max,n_max)
      double precision a_norm_frobenius
      double precision b(n_max,n_max)
      double precision b_norm_frobenius
      double precision c(n_max,n_max)
      double precision error_norm_frobenius
      integer n
      double precision r8mat_norm_fro

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST05'
      write ( *, '(a)' ) 
     &  '  BERNSTEIN_MATRIX returns a matrix A which transforms a'
      write ( *, '(a)' ) 
     &  '  polynomial coefficient vector from the power basis to'
      write ( *, '(a)' ) '  the Bernstein basis.'
      write ( *, '(a)' ) 
     &  '  BERNSTEIN_MATRIX_INVERSE computes the inverse B.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 
     &  '     N     ||A||            ||B||      ||I-A*B||'
      write ( *, '(a)' ) ' '

      do n = 5, 15

        call bernstein_matrix ( n, a )
        a_norm_frobenius = r8mat_norm_fro ( n, n, a )

        call bernstein_matrix_inverse ( n, b )
        b_norm_frobenius = r8mat_norm_fro ( n, n, b )

        call r8mat_mm ( n, n, n, a, b, c )

        call r8mat_is_identity ( n, c, error_norm_frobenius )

        write ( *, '(2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) 
     &    n, a_norm_frobenius, b_norm_frobenius, error_norm_frobenius

      end do

      return
      end

