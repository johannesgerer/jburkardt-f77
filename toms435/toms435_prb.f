      program main

c*********************************************************************72
c
cc TOMS435_PRB tests TOMS435.
c
c  Modified:
c
c    12 January 2006
c
c  Author:
c
c    John Burkardt
c
      implicit none

      call timestamp ( )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TOMS435_PRB'
      write ( *, '(a)' ) '  FORTRAN77 version'
      write ( *, '(a)' ) '  Test the TOMS435 library.'

      call test01 ( )
c
c  Terminate.
c
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TOMS435_PRB'
      write ( *, '(a)' ) '  Normal end of execution.'

      write ( *, '(a)' ) ' '
      call timestamp ( )

      stop
      end
      subroutine test01 ( )

c*********************************************************************72
c
cc TEST01 tests GAMINC.
c
c  Modified:
c
c    19 January 2008
c
c  Author:
c
c    John Burkardt
c
      implicit none

      real a
      real fx
      real fx2
      real gam
      real gaminc
      real gamma
      integer n_data
      real x
      real x1
      real x2

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST01'
      write ( *, '(a)' ) '  GAMAINC evaluates'
      write ( *, '(a)' ) '  the modified incomplete Gamma function.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Note that, for comparison, we must'
      write ( *, '(a)' ) '  divide GAMINC by Gamma(A).'
      write ( *, '(a)' ) ' '
      write ( *, '(a,a)' )
     &  '      A         X         Exact Value       ',
     &  'Computed        Diff'
      write ( *, '(a)' ) ' '

      n_data = 0

10    continue

      call gamma_inc_values ( n_data, a, x, fx )

      if ( n_data <= 0 ) then
        go to 20
      end if

      gam = gamma ( a )

      x1 = 0.0
      x2 = x

      fx2 = gaminc ( a, x1, x2, gam ) / gam

      write ( *, '(2x,f8.4,2x,f8.4,2x,g16.8,2x,g16.8,2x,g10.4)' )
     &  a, x, fx, fx2, abs ( fx - fx2 )

      go to 10

20    continue

      return
      end
