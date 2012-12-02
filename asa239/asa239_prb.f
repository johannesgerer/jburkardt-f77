      program main

c*********************************************************************72
c
cc MAIN is the main program for ASA239_PRB.
c
c  Discussion:
c
c    ASA239_PRB calls the ASA239 routines.
c
c  Modified:
c
c    04 January 2008
c
c  Author:
c
c    John Burkardt
c
      implicit none

      call timestamp ( )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'ASA239_PRB:'
      write ( *, '(a)' ) '  FORTRAN77 version'
      write ( *, '(a)' ) '  Test the ASA239 library.'

      call test01 ( )
c
c  Terminate.
c
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'ASA239_PRB:'
      write ( *, '(a)' ) '  Normal end of execution.'

      write ( *, '(a)' ) ' '
      call timestamp ( )

      stop
      end
      subroutine test01 ( )

c*********************************************************************72
c
cc TEST01 demonstrates the use of GAMMAD.
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

      double precision a
      double precision fx
      double precision fx2
      double precision gammad
      integer ifault
      integer n_data
      double precision x

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST01:'
      write ( *, '(a)' ) '  GAMMAD computes the incomplete Gamma '
      write ( *, '(a)' ) '  function.  Compare to tabulated values.'
      write ( *, '(a)' ) ' '
      write ( *, '(a,a)' )
     &  '      A             X          ',
     &  '   FX                        FX                      DIFF'
      write ( *, '(a,a)' )
     &  '                               ',
     &  ' (tabulated)                (GAMMAD)'
      write ( *, '(a)' ) ' '

      n_data = 0

10    continue

        call gamma_inc_values ( n_data, a, x, fx )

        if ( n_data .eq. 0 ) then
          go to 20
        end if

        fx2 = gammad ( x, a, ifault )

        write ( *, '(2x,f12.8,2x,f12.8,2x,g24.16,2x,g24.16,2x,g10.4)' )
     &    a, x, fx, fx2, dabs ( fx - fx2 )

      go to 10

20    continue

      return
      end
