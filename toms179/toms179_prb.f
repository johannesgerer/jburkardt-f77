      program main

c*********************************************************************72
c
cc TOMS179_PRB tests TOMS179.
c
c  Modified:
c
c    03 January 2006
c
c  Author:
c
c    John Burkardt
c
      implicit none

      call timestamp ( )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TOMS179_PRB'
      write ( *, '(a)' ) '  FORTRAN77 version'
      write ( *, '(a)' ) '  Test the TOMS179 library.'

      call test01 ( )
      call test02 ( )
c
c  Terminate.
c
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TOMS179_PRB'
      write ( *, '(a)' ) '  Normal end of execution.'

      write ( *, '(a)' ) ' '
      call timestamp ( )

      stop
      end
      subroutine test01 ( )

c*********************************************************************72
c
cc TEST01 tests ALOGAM.
c
c  Modified:
c
c    30 January 2006
c
c  Author:
c
c    John Burkardt
c
      implicit none

      double precision alogam
      double precision fx
      double precision fx2
      integer ifault
      integer n_data
      double precision x

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST01'
      write ( *, '(a)' ) '  Test ALOGAM, which estimates the logarithm'
      write ( *, '(a)' ) '  of the Gamma function.'
      write ( *, '(a)' ) ' '
      write ( *, '(a,a)' ) '      X         Exact Value               ',
     &  'Computed                Diff'
      write ( *, '(a)' ) ' '

      n_data = 0

10    continue

      call gamma_log_values ( n_data, x, fx )

      if ( n_data <= 0 ) then
        go to 20
      end if

      fx2 = alogam ( x, ifault )

      write ( *, '(2x,f8.4,2x,g24.16,2x,g24.16,2x,g10.4)' )
     &  x, fx, fx2, dabs ( fx - fx2 )

      go to 10

20    continue

      return
      end
      subroutine test02 ( )

c*********************************************************************72
c
cc TEST02 tests MDBETA.
c
c  Modified:
c
c    03 January 2006
c
c  Author:
c
c    John Burkardt
c
      implicit none

      double precision fx
      double precision fx2
      integer ier
      integer n_data
      double precision p
      double precision q
      double precision x

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST02'
      write ( *, '(a)' ) '  Test MDBETA, which estimates the value of'
      write ( *, '(a)' ) '  the modified Beta function.'
      write ( *, '(a)' ) ' '
      write ( *, '(a,a)' ) '      X         P         Q         ',
     &  'Exact Value               Computed                Diff'
      write ( *, '(a)' ) ' '

      n_data = 0

10    continue

      call beta_cdf_values ( n_data, p, q, x, fx )

      if ( n_data <= 0 ) then
        go to 20
      end if

      call mdbeta ( x, p, q, fx2, ier )

      write ( *,
     &  '(2x,f8.4,2x,f8.4,2x,f8.4,2x,g24.16,2x,g24.16,2x,g10.4)' )
     &  x, p, q, fx, fx2, dabs ( fx - fx2 )

      go to 10

20    continue

      return
      end

