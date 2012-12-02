      program main

c*********************************************************************72
c
cc MAIN is the main program for SPHERE_LEBEDEV_RULE_PRB.
c
c  Discussion:
c
c    SPHERE_LEBEDEV_RULE_PRB tests routines from the SPHERE_LEBEDEV_RULE library.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    13 September 2010
c
c  Author:
c
c    John Burkardt
c
      implicit none

      call timestamp ( )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'SPHERE_LEBEDEV_RULE_PRB'
      write ( *, '(a)' ) '  FORTRAN77 version'
      write ( *, '(a)' ) '  Test the SPHERE_LEBEDEV_RULE library.'

      call test01 ( )
      call test02 ( )
c
c  Terminate.
c
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'SPHERE_LEBEDEV_RULE_PRB'
      write ( *, '(a)' ) '  Normal end of execution.'

      write ( *, '(a)' ) ' '
      call timestamp ( )

      stop
      end
      subroutine test01 ( )

c*********************************************************************72
c
cc TEST01 tests AVAILABLE_TABLE, ORDER_TABLE, PRECISION_TABLE.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c   12 September 2010
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer available
      integer available_table
      integer order
      integer order_table
      integer precision
      integer precision_table
      integer rule
      integer rule_max
      parameter ( rule_max = 65 )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST01'
      write ( *, '(a)' ) '  List Lebedev rule properties.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Rule Avail Order  Prec'
      write ( *, '(a)' ) ' '
      do rule = 1, rule_max
        available = available_table ( rule )
        order = order_table ( rule )
        precision = precision_table ( rule )
        write ( *, '(2x,i4,2x,i4,2x,i4,2x,i4)' ) 
     &    rule, available, order, precision
      end do

      return
      end
      subroutine test02 ( )

c*********************************************************************72
c
cc TEST02 tests the SPHERE_LEBEDEV_RULE functions.
c
c  Modified:
c
c    13 September 2010
c
c  Author:
c
c    Dmitri Laikov
c
c  Reference:
c
c    Vyacheslav Lebedev, Dmitri Laikov,
c    A quadrature formula for the sphere of the 131st
c    algebraic order of accuracy,
c    Russian Academy of Sciences Doklady Mathematics,
c    Volume 59, Number 3, 1999, pages 477-481.
c
      implicit none

      integer nmax
      parameter ( nmax = 65 )
      integer mmax
      parameter ( mmax = ( ( nmax * 2 + 3 ) * ( nmax * 2 + 3 ) / 3 ) )

      double precision alpha
      integer available
      integer available_table
      double precision beta
      double precision err
      double precision err_max
      integer i
      double precision integral_exact
      double precision integral_approx
      integer j
      integer k
      integer m
      integer n
      integer order
      integer order_table
      integer precision_table
      double precision s(0:nmax+1)
      double precision w(mmax)
      double precision x(mmax)
      double precision xn(mmax,0:nmax)
      double precision y(mmax)
      double precision yn(mmax,0:nmax)
      double precision z(mmax)
      double precision zn(mmax,0:nmax)

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST02'
      write ( *, '(a)' ) 
     &  '  Generate each available rule and test for accuracy.'

      do n = 1, nmax

        available = available_table ( n )

        if ( available .eq. 1 ) then

          order = order_table ( n )

          call ld_by_order ( order, x, y, z, w )

          s(0) = 1.0D+00
          do k = 1, n + 1
            s(k) = ( 2 * k - 1 ) * s(k-1)
          end do
c
c  For each abscissa X(M), compute the values 1, X(M)^2, X(M)^4, ..., X(M)^2*N.
c
          do m = 1, order
            xn(m,0) = 1.0D+00
            yn(m,0) = 1.0D+00
            zn(m,0) = 1.0D+00
            do k = 1, n
              xn(m,k) = xn(m,k-1) * x(m) * x(m)
              yn(m,k) = yn(m,k-1) * y(m) * y(m)
              zn(m,k) = zn(m,k-1) * z(m) * z(m)
            end do
          end do

          err_max = 0.0D+00
          do i = 0, n
            do j = 0, n - i
              k = n - i - j
c
c  Apply Lebedev rule to x^2i y^2j z^2k.
c
              integral_approx = 0.0D+00
              do m = 1, order
                integral_approx = integral_approx 
     &            + w(m) * xn(m,i) * yn(m,j) * zn(m,k)
              end do
c
c  Compute exact value of integral (aside from factor of 4 pi!).
c
              integral_exact = s(i) * s(j) * s(k) / s(1+i+j+k)
c
c  Record the maximum error for this rule.
c
              err = abs ( ( integral_approx - integral_exact ) 
     &          / integral_exact )

              if ( err_max .lt. err ) then
                err_max = err
              end if

            end do
          end do

          write ( *, '(a)' ) ' '
          write ( *, '(a,i4,a,i4,a,g14.6)' ) 
     &      '  Order = ', order,
     &      '  LMAXW = ', precision_table ( n ),
     &      '  max error = ', err_max
c
c  Convert (X,Y,Z) to (Theta,Phi) and print the data.
c
          if ( order .le. 50 ) then
            do m = 1, order
              call xyz_to_tp ( x(m), y(m), z(m), alpha, beta )
              write ( *, '(g24.15,g24.15,g24.15)' ) alpha, beta, w(m)
            end do
          end if

        end if

      end do

      return
      end
