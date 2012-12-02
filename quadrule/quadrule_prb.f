      program main

c*********************************************************************72
c
cc MAIN is the main program for QUADRULE_PRB.
c
c  Discussion:
c
c    QUADRULE_PRB calls a set of tests for the QUADRULE library.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    14 September 2012
c
c  Author:
c
c    John Burkardt
c
      implicit none

      double precision alpha
      integer n
      integer order

      call timestamp ( )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'QUADRULE_PRB'
      write ( *, '(a)' ) '  FORTRAN77 version'
      write ( *, '(a)' ) '  Test the QUADRULE library.'

      call test0725 ( )
      order = 31
      call test087 ( order )
      order = 63
      call test087 ( order )
      call test095 ( )
c
c  Terminate.
c
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'QUADRULE_PRB'
      write ( *, '(a)' ) '  Normal end of execution.'

      write ( *, '(a)' ) ' '
      call timestamp ( )
 
      stop
      end
      subroutine test0725 ( )

c*********************************************************************72
c
cc TEST0725 tests CLENSHAW_CURTIS_COMPUTE
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    23 April 2011
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n_max
      parameter ( n_max = 16 )

      integer i
      integer n
      double precision w(n_max)
      double precision x(n_max)

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST0725'
      write ( *, '(a)' ) '  CLENSHAW_CURTIS_COMPUTE computes'
      write ( *, '(a)' ) 
     &  '  a Clenshaw-Curtis quadrature rule over [-1,1]'
      write ( *, '(a)' ) '  of given order.'

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 
     &  '     Order       W                       X'
      write ( *, '(a)' ) ' '

      do n = 1, 10

        call clenshaw_curtis_compute ( n, x, w )

        write ( *, '(a)' ) ' '
        write ( *, '(2x,i8)' ) n

        do i = 1, n
          write ( *, '(10x,2x,g25.16,2x,f24.16)' ) w(i), x(i)
        end do

      end do

      return
      end
      subroutine test087 ( order )

c*********************************************************************72
c
cc TEST087 tests HERMITE_EK_COMPUTE.
c
c  Discussion:
c
c    I used this test to generate tabular values of weights and
c    abscissas for Gauss-Hermite quadrature.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    14 September 2012
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer ORDER, the order of the rule.
c
      implicit none

      integer order

      integer i
      double precision w(order)
      double precision x(order)

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST087'
      write ( *, '(a)' ) 
     &  '  HERMITE_EK_COMPUTE computes a Gauss-Hermite rule;'
      write ( *, '(a)' ) ' '
      write ( *, '(a,i8)' ) '  Compute the data for ORDER = ', order

      call hermite_ek_compute ( order, x, w )
 
      write ( *, '(a)' ) ' '
      do i = 1, order
        write ( *, '(a,i3,a,g40.32)' ) '    x(', i, ') = ', x(i)
      end do
      write ( *, '(a)' ) ' '
      do i = 1, order
        write ( *, '(a,i3,a,g40.32)' ) '    w(', i, ') = ', w(i)
      end do

      return
      end
      subroutine test095 ( )

c*********************************************************************72
c
cc TEST095 tests HERMITE_GENZ_KEISTER_SET and SUMMER.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    02 June 2010
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer l_max
      parameter ( l_max = 8 )
      integer n_max
      parameter ( n_max = 35 )
 
      double precision error
      double precision estimate
      double precision exact
      double precision f(n_max)
      integer i
      integer l
      integer m
      integer n
      integer n_list(0:l_max)
      integer p
      integer p_list(0:l_max)
      double precision r8vec_dot_product
      double precision w(n_max)
      double precision x(n_max)

      save n_list
      save p_list

      data n_list / 1, 3, 7, 9, 17, 19, 31, 33, 35 /
      data p_list / 1, 5, 7, 15, 17, 29, 31, 33, 51 /

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST095'
      write ( *, '(a)' ) 
     &  '  HERMITE_GENZ_KEISTER_SET sets up a nested rule'
      write ( *, '(a)' ) '  for the Hermite integration problem.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  The integration interval is ( -oo, +oo ).'
      write ( *, '(a)' ) '  The weight function is exp ( - x * x )'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  HERMITE_INTEGRAL determines the exact value'
      write ( *, '(a)' ) '  of the integal when f(x) = x^m.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '         M         N       Estimate' //
     &  '       Exact            Error'

      do l = 0, l_max

        write ( *, '(a)' ) ' '
        n = n_list(l)
        p = p_list(l)

        call hermite_genz_keister_set ( n, x, w )

        do m = 0, min ( p + 2, 20 ), 2

          call hermite_integral ( m, exact )

          if ( m .eq. 0 ) then
            do i = 1, n
              f(i) = 1.0D+00
            end do
          else
            do i = 1, n
              f(i) = x(i)**m
            end do
          end if

          estimate = r8vec_dot_product ( n, w, f )

          error = abs ( exact - estimate )
  
          write ( *, '(2x,i8,2x,i8,2x,g14.6,2x,g14.6,2x,g14.6)' ) 
     &      m, n, estimate, exact, error

        end do

      end do

      return
      end
