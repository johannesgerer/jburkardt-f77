      program main

c*********************************************************************72
c
cc TEST_INTERP_1D_TEST tests the TEST_INTERP_1D library.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    07 October 2012
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer nd

      call timestamp ( )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST_INTERP_1D_TEST'
      write ( *, '(a)' ) '  FORTRAN77 version'
      write ( *, '(a)' ) '  Test the TEST_INTERP_1D library.'
      write ( *, '(a)' ) '  The R8LIB library is needed.'

      call test01 ( )

      nd = 11
      call test02 ( nd )
c
c  Terminate.
c
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST_INTERP_1D_TEST'
      write ( *, '(a)' ) '  Normal end of execution.'

      write ( *, '(a)' ) ' '
      call timestamp ( )

      return;
      end
      subroutine test01 ( )

c*********************************************************************72
c
cc TEST01 simply prints the title of each function.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    28 August 2012
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer prob
      integer prob_num
      character * ( 80 ) title

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST01'
      write ( *, '(a)' ) '  Print the title of each function.'

      call p00_prob_num ( prob_num )
      
      write ( *, '(a)' ) ' '
      write ( *, '(a,i2,a)' ) 
     &  '  There are ', prob_num, ' functions available:'
      write ( *, '(a)' ) '  Index  Title'
      write ( *, '(a)' ) ' '

      do prob = 1, prob_num

        call p00_title ( prob, title )

        write ( *, '(2x,i2,2x,a)' ) prob, trim ( title )

      end do

      return
      end
      subroutine test02 ( nd )

c*********************************************************************72
c
cc TEST_INTERP_1D_TEST02 evaluates each test function at ND sample points.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    28 August 2012
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer ND, the number of sample points.
c
      implicit none

      integer nd

      double precision a
      double precision b
      integer prob
      integer prob_num
      character * ( 80 ) title
      double precision x(nd)
      double precision y(nd)

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST_INTERP_1D_TEST02'
      write ( *, '(a)' ) '  Use P00_F to sample each function.'

      call p00_prob_num ( prob_num )

      a = 0.0D+00
      b = 1.0D+00
      call r8vec_linspace ( nd, a, b, x )

      write ( *, '(a)' ) ' '

      do prob = 1, prob_num

        call p00_f ( prob, nd, x, y )
        write ( title, '(a,i2)' ) 'X, Y for problem ', prob
        call r8vec2_print ( nd, x, y, trim ( title ) )

      end do

      return
      end
