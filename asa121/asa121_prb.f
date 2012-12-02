      program main

c*********************************************************************72
c
cc MAIN is the main program for ASA121_PRB.
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

      write ( *, '(a)' ) ' '
      call timestamp ( )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'ASA121_PRB:'
      write ( *, '(a)' ) '  FORTRAN77 version'
      write ( *, '(a)' ) '  Test the ASA121 library.'

      call test01 ( )
c
c  Terminate.
c
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'ASA121_PRB:'
      write ( *, '(a)' ) '  Normal end of execution.'

      write ( *, '(a)' ) ' '
      call timestamp ( )

      stop
      end
      subroutine test01 ( )

c*********************************************************************72
c
cc TEST01 compare TRIGAM against tabulated values.
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

      double precision fx
      double precision fx2
      integer ifault
      integer n_data
      double precision trigam
      double precision x

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST01:'
      write ( *, '(a)' ) '  Compare tabulated values of the trigamma'
      write ( *, '(a)' ) '  function against values computed'
      write ( *, '(a)' ) '  computed by TRIGAM.'
      write ( *, '(a)' ) ' '
      write ( *, '(a,a)' )
     &  '        X          FX                        FX     ',
     &  '                DIFF'
      write ( *, '(a)' )
     &  '                  (Tabulated)               (TRIGAM)'
      write ( *, '(a)' ) ' '

      n_data = 0

10    continue

        call trigamma_values ( n_data, x, fx )

        if ( n_data .eq. 0 ) then
          go to 20
        end if

        fx2 = trigam ( x, ifault )

        write ( *, '(2x,f10.4,2x,g24.16,2x,g24.16,2x,g10.4)' )
     &  x, fx, fx2, dabs ( fx - fx2 )

      go to 10

20    continue

      return
      end
