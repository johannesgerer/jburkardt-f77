      program main

c*********************************************************************72
c
cc MAIN is the main program for SPHERE_QUAD_PRB.
c
c  Discussion:
c
c    SPHERE_QUAD_PRB tests SPHERE_QUAD.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    23 September 2010
c
c  Author:
c
c    John Burkardt
c
      implicit none

      call timestamp ( )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'SPHERE_QUAD_PRB'
      write ( *, '(a)' ) '  FORTRAN77 version'
      write ( *, '(a)' ) '  Test the SPHERE_QUAD library.'

      call test01 ( )
      call test02 ( )
      call test03 ( )
      call test04 ( )
      call test05 ( )
c
c  Terminate.
c
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'SPHERE_QUAD_PRB'
      write ( *, '(a)' ) '  Normal end of execution.'

      write ( *, '(a)' ) ' '
      call timestamp ( )

      stop
      end

      subroutine test01 ( )

c*********************************************************************72
c
cc TEST01 tests SPHERE01_QUAD_LL*.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    23 September 2010
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer e(3)
      double precision exact
      double precision polyterm_value_3d
      external polyterm_value_3d
      double precision h
      integer h_test
      integer i
      integer n_llc
      integer n_llm
      integer n_llv
      integer n_mc
      double precision result_llc
      double precision result_llm
      double precision result_llv
      double precision result_mc
      integer seed

      seed = 123456789

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST01'
      write ( *, '(a)' ) 
     &  '  Approximate the integral of a function on the unit sphere.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  SPHERE01_QUAD_MC uses a Monte Carlo method.'
      write ( *, '(a)' ) 
     &  '  SPHERE01_QUAD_LLC uses centroids of spherical triangles.'
      write ( *, '(a)' ) 
     &  '  SPHERE01_QUAD_LLM uses midsides of spherical triangles.'
      write ( *, '(a)' ) 
     &  '  SPHERE01_QUAD_LLV uses vertices of spherical triangles.'
      write ( *, '(a)' ) ' '
      write ( *, '(a,a)' ) 
     &  '  H              QUAD_MC       QUAD_LLC      QUAD_LLM      ', 
     &  'QUAD_LLV         EXACT'

      do i = 0, 17

        if ( i == 0 ) then
          e(1) = 0
          e(2) = 0
          e(3) = 0
        else if ( i == 1 ) then
          e(1) = 0
          e(2) = 0
          e(3) = 0
        else if ( i == 2 ) then
          e(1) = 1
          e(2) = 0
          e(3) = 0
        else if ( i == 3 ) then
          e(1) = 0
          e(2) = 1
          e(3) = 0
        else if ( i == 4 ) then
          e(1) = 0
          e(2) = 0
          e(3) = 1
        else if ( i == 5 ) then
          e(1) = 2
          e(2) = 0
          e(3) = 0
        else if ( i == 6 ) then
          e(1) = 0
          e(2) = 2
          e(3) = 2
        else if ( i == 7 ) then
          e(1) = 2
          e(2) = 2
          e(3) = 2
        else if ( i == 8 ) then
          e(1) = 0
          e(2) = 2
          e(3) = 4
        else if ( i == 9 ) then
          e(1) = 0
          e(2) = 0
          e(3) = 6
        else if ( i == 10 ) then
          e(1) = 1
          e(2) = 2
          e(3) = 4
        else if ( i == 11 ) then
          e(1) = 2
          e(2) = 4
          e(3) = 2
        else if ( i == 12 ) then
          e(1) = 6
          e(2) = 2
          e(3) = 0
        else if ( i == 13 ) then
          e(1) = 0
          e(2) = 0
          e(3) = 8
        else if ( i == 14 ) then
          e(1) = 6
          e(2) = 0
          e(3) = 4
        else if ( i == 15 ) then
          e(1) = 4
          e(2) = 6
          e(3) = 2
        else if ( i == 16 ) then
          e(1) = 2
          e(2) = 4
          e(3) = 8
        else if ( i == 17 ) then
          e(1) = 16
          e(2) = 0
          e(3) = 0
        end if

        call polyterm_exponent ( 'SET', e )

        if ( i .eq. 0 ) then
          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) 'Point counts per method:'
        else
          call polyterm_exponent ( 'PRINT', e )
        end if

        do h_test = 1, 3

          if ( h_test .eq. 1 ) then
            h = 1.0D+00
          else if (  h_test .eq. 2 ) then
            h = 0.1D+00
          else if (  h_test .eq. 3 ) then
            h = 0.01D+00
          end if

          call sphere01_quad_mc_size ( h, n_mc )

          call sphere01_quad_mc ( polyterm_value_3d, h, seed, n_mc, 
     &      result_mc )

          call sphere01_quad_llc ( polyterm_value_3d, h, n_llc, 
     &      result_llc )

          call sphere01_quad_llm ( polyterm_value_3d, h, n_llm, 
     &      result_llm )

          call sphere01_quad_llv ( polyterm_value_3d, h, n_llv, 
     &      result_llv )

          call sphere01_monomial_int ( e, exact )

          if ( i .eq. 0 ) then
            write ( *, '(g14.6,5i14)' ) h, n_mc, n_llc, n_llm, n_llv
          else
            write ( *, '(6g14.6)' ) 
     &        h, result_mc, result_llc, result_llm, result_llv, exact

          end if

        end do

      end do

      return
      end
      subroutine test02 ( )

c*********************************************************************72
c
cc TEST02 tests SPHERE01_QUAD_ICOS1C.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    23 September 2010
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer e(3)
      double precision error
      double precision exact
      double precision polyterm_value_3d
      external polyterm_value_3d
      integer factor
      integer factor_log
      integer i
      integer n
      double precision result

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST02'
      write ( *, '(a)' ) 
     &  '  Approximate the integral of a function on the unit sphere.'
      write ( *, '(a)' ) 
     &  '  SPHERE01_QUAD_ICOS1C uses centroids of spherical triangles.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 
     &  'FACTOR         N        QUAD          EXACT         ERROR'

      do i = 1, 17
     
        if ( i == 1 ) then
          e(1) = 0
          e(2) = 0
          e(3) = 0
        else if ( i .eq. 2 ) then
          e(1) = 1
          e(2) = 0
          e(3) = 0
        else if ( i .eq. 3 ) then
          e(1) = 0
          e(2) = 1
          e(3) = 0
        else if ( i .eq. 4 ) then
          e(1) = 0
          e(2) = 0
          e(3) = 1
        else if ( i .eq. 5 ) then
          e(1) = 2
          e(2) = 0
          e(3) = 0
        else if ( i .eq. 6 ) then
          e(1) = 0
          e(2) = 2
          e(3) = 2
        else if ( i .eq. 7 ) then
          e(1) = 2
          e(2) = 2
          e(3) = 2
        else if ( i .eq. 8 ) then
          e(1) = 0
          e(2) = 2
          e(3) = 4
        else if ( i .eq. 9 ) then
          e(1) = 0
          e(2) = 0
          e(3) = 6
        else if ( i .eq. 10 ) then
          e(1) = 1
          e(2) = 2
          e(3) = 4
        else if ( i .eq. 11 ) then
          e(1) = 2
          e(2) = 4
          e(3) = 2
        else if ( i .eq. 12 ) then
          e(1) = 6
          e(2) = 2
          e(3) = 0
        else if ( i .eq. 13 ) then
          e(1) = 0
          e(2) = 0
          e(3) = 8
        else if ( i .eq. 14 ) then
          e(1) = 6
          e(2) = 0
          e(3) = 4
        else if ( i .eq. 15 ) then
          e(1) = 4
          e(2) = 6
          e(3) = 2
        else if ( i .eq. 16 ) then
          e(1) = 2
          e(2) = 4
          e(3) = 8
        else if ( i .eq. 17 ) then
          e(1) = 16
          e(2) = 0
          e(3) = 0
        end if

        call polyterm_exponent ( 'SET', e )

        call polyterm_exponent ( 'PRINT', e )

        factor = 1
        do factor_log = 0, 5

          call sphere01_quad_icos1c ( factor, polyterm_value_3d, n, 
     &      result )

          call sphere01_monomial_int ( e, exact )

          error = abs ( exact - result )

          write ( *, '(2x,i4,2x,i8,2x,3g14.6)' ) 
     &        factor, n, result, exact, error

          factor = factor * 2

        end do

      end do

      return
      end
      subroutine test03 ( )

c*********************************************************************72
c
cc TEST03 tests SPHERE01_QUAD_ICOS1M.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    23 September 2010
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer e(3)
      double precision error
      double precision exact
      double precision polyterm_value_3d
      external polyterm_value_3d
      integer factor
      integer factor_log
      integer i
      integer n
      double precision result

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST03'
      write ( *, '(a)' ) 
     &  '  Approximate the integral of a function on the unit sphere.'
      write ( *, '(a)' ) 
     &  '  SPHERE01_QUAD_ICOS1M uses midpoints of spherical triangles.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 
     &  'FACTOR         N        QUAD          EXACT         ERROR'

      do i = 1, 17
     
        if ( i == 1 ) then
          e(1) = 0
          e(2) = 0
          e(3) = 0
        else if ( i .eq. 2 ) then
          e(1) = 1
          e(2) = 0
          e(3) = 0
        else if ( i .eq. 3 ) then
          e(1) = 0
          e(2) = 1
          e(3) = 0
        else if ( i .eq. 4 ) then
          e(1) = 0
          e(2) = 0
          e(3) = 1
        else if ( i .eq. 5 ) then
          e(1) = 2
          e(2) = 0
          e(3) = 0
        else if ( i .eq. 6 ) then
          e(1) = 0
          e(2) = 2
          e(3) = 2
        else if ( i .eq. 7 ) then
          e(1) = 2
          e(2) = 2
          e(3) = 2
        else if ( i .eq. 8 ) then
          e(1) = 0
          e(2) = 2
          e(3) = 4
        else if ( i .eq. 9 ) then
          e(1) = 0
          e(2) = 0
          e(3) = 6
        else if ( i .eq. 10 ) then
          e(1) = 1
          e(2) = 2
          e(3) = 4
        else if ( i .eq. 11 ) then
          e(1) = 2
          e(2) = 4
          e(3) = 2
        else if ( i .eq. 12 ) then
          e(1) = 6
          e(2) = 2
          e(3) = 0
        else if ( i .eq. 13 ) then
          e(1) = 0
          e(2) = 0
          e(3) = 8
        else if ( i .eq. 14 ) then
          e(1) = 6
          e(2) = 0
          e(3) = 4
        else if ( i .eq. 15 ) then
          e(1) = 4
          e(2) = 6
          e(3) = 2
        else if ( i .eq. 16 ) then
          e(1) = 2
          e(2) = 4
          e(3) = 8
        else if ( i .eq. 17 ) then
          e(1) = 16
          e(2) = 0
          e(3) = 0
        end if

        call polyterm_exponent ( 'SET', e )

        call polyterm_exponent ( 'PRINT', e )

        factor = 1
        do factor_log = 0, 5

          call sphere01_quad_icos1m ( factor, polyterm_value_3d, n, 
     &      result )

          call sphere01_monomial_int ( e, exact )

          error = abs ( exact - result )

          write ( *, '(2x,i4,2x,i8,2x,3g14.6)' ) 
     &        factor, n, result, exact, error

          factor = factor * 2

        end do

      end do

      return
      end
      subroutine test04 ( )

c*********************************************************************72
c
cc TEST04 tests SPHERE01_QUAD_ICOS1V.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    23 September 2010
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer e(3)
      double precision error
      double precision exact
      double precision polyterm_value_3d 
      external polyterm_value_3d
      integer factor
      integer factor_log
      integer i
      integer n
      double precision result

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST04'
      write ( *, '(a)' ) 
     &  '  Approximate the integral of a function on the unit sphere.'
      write ( *, '(a)' ) 
     &  '  SPHERE01_QUAD_ICOS1V uses vertices of spherical triangles.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 
     &  'FACTOR         N        QUAD          EXACT         ERROR'

      do i = 1, 17
     
        if ( i == 1 ) then
          e(1) = 0
          e(2) = 0
          e(3) = 0
        else if ( i .eq. 2 ) then
          e(1) = 1
          e(2) = 0
          e(3) = 0
        else if ( i .eq. 3 ) then
          e(1) = 0
          e(2) = 1
          e(3) = 0
        else if ( i .eq. 4 ) then
          e(1) = 0
          e(2) = 0
          e(3) = 1
        else if ( i .eq. 5 ) then
          e(1) = 2
          e(2) = 0
          e(3) = 0
        else if ( i .eq. 6 ) then
          e(1) = 0
          e(2) = 2
          e(3) = 2
        else if ( i .eq. 7 ) then
          e(1) = 2
          e(2) = 2
          e(3) = 2
        else if ( i .eq. 8 ) then
          e(1) = 0
          e(2) = 2
          e(3) = 4
        else if ( i .eq. 9 ) then
          e(1) = 0
          e(2) = 0
          e(3) = 6
        else if ( i .eq. 10 ) then
          e(1) = 1
          e(2) = 2
          e(3) = 4
        else if ( i .eq. 11 ) then
          e(1) = 2
          e(2) = 4
          e(3) = 2
        else if ( i .eq. 12 ) then
          e(1) = 6
          e(2) = 2
          e(3) = 0
        else if ( i .eq. 13 ) then
          e(1) = 0
          e(2) = 0
          e(3) = 8
        else if ( i .eq. 14 ) then
          e(1) = 6
          e(2) = 0
          e(3) = 4
        else if ( i .eq. 15 ) then
          e(1) = 4
          e(2) = 6
          e(3) = 2
        else if ( i .eq. 16 ) then
          e(1) = 2
          e(2) = 4
          e(3) = 8
        else if ( i .eq. 17 ) then
          e(1) = 16
          e(2) = 0
          e(3) = 0
        end if

        call polyterm_exponent ( 'SET', e )

        call polyterm_exponent ( 'PRINT', e )

        factor = 1
        do factor_log = 0, 5

          call sphere01_quad_icos1v ( factor, polyterm_value_3d, n, 
     &      result )

          call sphere01_monomial_int ( e, exact )

          error = abs ( exact - result )

          write ( *, '(2x,i4,2x,i8,2x,3g14.6)' ) 
     &        factor, n, result, exact, error

          factor = factor * 2

        end do

      end do

      return
      end
      subroutine test05 ( )

c*********************************************************************72
c
cc TEST05 tests SPHERE01_QUAD_ICOS2V.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    23 September 2010
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer e(3)
      double precision error
      double precision exact
      double precision polyterm_value_3d
      external polyterm_value_3d
      integer factor
      integer factor_log
      integer i
      integer n
      double precision result

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST05'
      write ( *, '(a)' ) 
     &  '  Approximate the integral of a function on the unit sphere.'
      write ( *, '(a)' ) 
     &  '  SPHERE01_QUAD_ICOS2V uses vertices of spherical triangles.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 
     &  'FACTOR         N        QUAD          EXACT         ERROR'

      do i = 1, 17
     
        if ( i == 1 ) then
          e(1) = 0
          e(2) = 0
          e(3) = 0
        else if ( i .eq. 2 ) then
          e(1) = 1
          e(2) = 0
          e(3) = 0
        else if ( i .eq. 3 ) then
          e(1) = 0
          e(2) = 1
          e(3) = 0
        else if ( i .eq. 4 ) then
          e(1) = 0
          e(2) = 0
          e(3) = 1
        else if ( i .eq. 5 ) then
          e(1) = 2
          e(2) = 0
          e(3) = 0
        else if ( i .eq. 6 ) then
          e(1) = 0
          e(2) = 2
          e(3) = 2
        else if ( i .eq. 7 ) then
          e(1) = 2
          e(2) = 2
          e(3) = 2
        else if ( i .eq. 8 ) then
          e(1) = 0
          e(2) = 2
          e(3) = 4
        else if ( i .eq. 9 ) then
          e(1) = 0
          e(2) = 0
          e(3) = 6
        else if ( i .eq. 10 ) then
          e(1) = 1
          e(2) = 2
          e(3) = 4
        else if ( i .eq. 11 ) then
          e(1) = 2
          e(2) = 4
          e(3) = 2
        else if ( i .eq. 12 ) then
          e(1) = 6
          e(2) = 2
          e(3) = 0
        else if ( i .eq. 13 ) then
          e(1) = 0
          e(2) = 0
          e(3) = 8
        else if ( i .eq. 14 ) then
          e(1) = 6
          e(2) = 0
          e(3) = 4
        else if ( i .eq. 15 ) then
          e(1) = 4
          e(2) = 6
          e(3) = 2
        else if ( i .eq. 16 ) then
          e(1) = 2
          e(2) = 4
          e(3) = 8
        else if ( i .eq. 17 ) then
          e(1) = 16
          e(2) = 0
          e(3) = 0
        end if

        call polyterm_exponent ( 'SET', e )

        call polyterm_exponent ( 'PRINT', e )

        factor = 1
        do factor_log = 0, 5

          call sphere01_quad_icos2v ( factor, polyterm_value_3d, n, 
     &      result )

          call sphere01_monomial_int ( e, exact )

          error = abs ( exact - result )

          write ( *, '(2x,i4,2x,i8,2x,3g14.6)' ) 
     &        factor, n, result, exact, error

          factor = factor * 2

        end do

      end do

      return
      end
      subroutine polyterm_exponent ( action, e )

c*********************************************************************72
c
cc POLYTERM_EXPONENT gets or sets the exponents for the polynomial term.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    08 June 2002
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, character * ( * ) ACTION.
c    'GET' asks the routine to return the current values in E.
c    'SET' asks the routine to set the current values to E.
c
c    Input/output, integer E(3), storage used to set or get values.
c
      implicit none

      character * ( * ) action
      integer e(3)
      integer e_save(3)
      integer i

      save e_save

      data e_save / 0, 0, 0 /

      if ( action(1:1) .eq. 'G' ) then

        do i = 1, 3
          e(i) = e_save(i)
        end do

      else if ( action(1:1) .eq. 'P' ) then

        write ( *, '(a)' ) ' '
        write ( *, '(a,i2,a,i2,a,i2)' ) 'P(X,Y,Z) = X^', e_save(1), 
     &                                            ' Y^', e_save(2), 
     &                                            ' Z^', e_save(3)
        
      else if ( action(1:1) .eq. 'S' ) then

        do i = 1, 3
          e_save(i) = e(i)
        end do

      end if

      return
      end
      subroutine polyterm_value_3d ( n, x, f )

c*********************************************************************72
c
cc POLYTERM_VALUE_3D evaluates a single polynomial term in 3D.
c
c  Discussion:
c
c    The polynomial term has the form:
c
c      F(X) = X(1)^E(1) * X(2)^E(2) * X(3)^E(3)
c
c    The exponents E(1:3) are set by calling POLYTERM_EXPONENT.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    23 September 2010
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the number of points.
c
c    Input, double precision X(3,N), the points where the polynomial term 
c    is to be evaluated.
c
c    Output, double precision F(N), the value of the polynomial term.
c
      implicit none

      integer n

      integer e(3)
      double precision f(n)
      integer i
      integer j
      double precision x(3,n)

      call polyterm_exponent ( 'GET', e )

      do j = 1, n
        f(j) = 1.0D+00
      end do

      do i = 1, 3

        if ( e(i) .ne. 0 ) then
          do j = 1, n
            f(j) = f(j) * x(i,j)**e(i)
          end do
        end if

      end do
  
      return
      end
