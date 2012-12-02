      program main

c*********************************************************************72
c
cc MAIN is the main program for BURGERS_SOLUTION_PRB.
c
c  Discussion:
c
c    BURGERS_SOLUTION_PRB tests the BURGERS_SOLUTION library.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    17 November 2011
c
c  Author:
c
c    John Burkardt
c
      implicit none

      call timestamp ( )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'BURGERS_SOLUTION_PRB'
      write ( *, '(a)' ) '  FORTRAN77 version'
      write ( *, '(a)' ) '  Test the BURGERS_SOLUTION library.'

      call test01 ( )
      call test02 ( )
c
c  Terminate.
c
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'BURGERS_SOLUTION_PRB'
      write ( *, '(a)' ) '  Normal end of execution.'

      write ( *, '(a)' ) ' '
      call timestamp ( )

      stop
      end
      subroutine test01 ( )

c*********************************************************************72
c
cc TEST01 tests sets up a small test case.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    17 November 2011
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer vtn
      parameter ( vtn = 11 )
      integer vxn 
      parameter ( vxn = 11 )

      character * ( 80 ) filename
      double precision nu
      double precision pi
      parameter ( pi = 3.141592653589793D+00 )
      double precision thi
      double precision tlo
      double precision vu(vxn,vtn)
      double precision vt(vtn)
      double precision vx(vxn)
      double precision xhi
      double precision xlo

      nu = 0.01D+00 / pi

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST01'
      write ( *, '(a)' ) 
     &  '  Compute an analytic solution to the Burgers equation.'
      write ( *, '(a)' ) ' '
      write ( *, '(a,g14.6)' ) '  Viscosity NU = ', nu
      write ( *, '(a,i4)' ) '  NX = ', vxn
      write ( *, '(a,i4)' ) '  NT = ', vtn

      xlo = -1.0D+00
      xhi = +1.0D+00
      call r8vec_even ( vxn, xlo, xhi, vx )
      call r8vec_print ( vxn, vx, '  X grid points:' )

      tlo = 0.0D+00
      thi = 3.0D+00 / pi
      call r8vec_even ( vtn, tlo, thi, vt )
      call r8vec_print ( vtn, vt, '  T grid points:' )

      call burgers_solution ( nu, vxn, vx, vtn, vt, vu )

      call r8mat_print ( vxn, vtn, vu, '  U(X,T) at grid points:' )

      filename = 'burgers_test01.txt'

      call r8mat_write ( filename, vxn, vtn, vu )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 
     &  '  Data written to file "' // trim ( filename ) // '".'

      return
      end
      subroutine test02 ( )

c*********************************************************************72
c
cc TEST02 tests sets up a finer test case.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    17 November 2011
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer vtn
      parameter ( vtn = 41 )
      integer vxn
      parameter ( vxn = 41 )

      character * ( 80 ) filename
      double precision nu
      double precision pi
      parameter ( pi = 3.141592653589793D+00 )
      double precision thi
      double precision tlo
      double precision vu(vxn,vtn)
      double precision vt(vtn)
      double precision vx(vxn)
      double precision xhi
      double precision xlo

      nu = 0.01D+00 / pi

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST02'
      write ( *, '(a)' ) 
     &  '  Compute an analytic solution to the Burgers equation.'
      write ( *, '(a)' ) ' '
      write ( *, '(a,g14.6)' ) '  Viscosity NU = ', nu
      write ( *, '(a,i4)' ) '  NX = ', vxn
      write ( *, '(a,i4)' ) '  NT = ', vtn

      xlo = -1.0D+00
      xhi = +1.0D+00
      call r8vec_even ( vxn, xlo, xhi, vx )
      call r8vec_print ( vxn, vx, '  X grid points:' )

      tlo = 0.0D+00
      thi = 3.0D+00 / pi
      call r8vec_even ( vtn, tlo, thi, vt )
      call r8vec_print ( vtn, vt, '  T grid points:' )

      call burgers_solution ( nu, vxn, vx, vtn, vt, vu )

      filename = 'burgers_test02.txt'

      call r8mat_write ( filename, vxn, vtn, vu )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 
     &  '  Data written to file "' // trim ( filename ) // '".'

      return
      end
