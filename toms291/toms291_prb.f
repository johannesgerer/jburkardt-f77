      program main

c*********************************************************************72
c
cc MAIN is the main program for TOMS291_PRB.
c
c  Modified:
c
c    28 December 2007
c
c  Author:
c
c    John Burkardt
c
      implicit none

      write ( *, '(a)' ) ' '
      call timestamp ( )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TOMS291_PRB'
      write ( *, '(a)' ) '  FORTRAN77 version'
      write ( *, '(a)' ) '  Test the TOMS291 library.'

      call test01 ( )
c
c  Terminate.
c
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TOMS291_PRB:'
      write ( *, '(a)' ) '  Normal end of execution.'
      write ( *, '(a)' ) ' '
      call timestamp ( )

      stop
      end
      subroutine test01 ( )

c*********************************************************************72
c
cc TEST01 compares results from ALOGAM with tabulated data.
c
c  Modified:
c
c    07 January 2008
c
c  Author:
c
c    John Burkardt
c
      implicit none

      double precision alogam
      double precision fx1
      double precision fx2
      integer ifault
      integer n_data
      double precision x

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST01:'
      write ( *, '(a)' )
     &'  ALOGAM computes the logarithm of the Gamma function.'
      write ( *, '(a)' ) '  Compare against tabulated data.'
      write ( *, '(a)' ) ' '
      write ( *, '(a,a)' ) '          X        ',
     &  '   FX                        FX'
      write ( *, '(a,a)' ) '                   ',
     &  'Tabulated                  ALOGAM'
      write ( *, '(a)' ) ' '

      n_data = 0

10    continue

        call gamma_log_values ( n_data, x, fx1 )

        if ( n_data .eq. 0 ) then
          go to 20
        end if

        fx2 = alogam ( x, ifault )

        write ( *, '(2x,f14.6,2x,g24.16,2x,g24.16)' ) x, fx1, fx2

      go to 10

20    continue

      return
      end
