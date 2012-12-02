      program main

c*********************************************************************72
c
cc MAIN is the main program for ASA245_PRB.
c
c  Discussion:
c
c    ASA245_PRB calls the ASA245 routines.
c
c  Modified:
c
c    13 January 2008
c
c  Author:
c
c    John Burkardt
c
      implicit none

      call timestamp ( )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'ASA245_PRB:'
      write ( *, '(a)' ) '  FORTRAN77 version'
      write ( *, '(a)' ) '  Test the ASA245 library.'

      call test01 ( )
      call test02 ( )
c
c  Terminate.
c
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'ASA245_PRB:'
      write ( *, '(a)' ) '  Normal end of execution.'

      write ( *, '(a)' ) ' '
      call timestamp ( )

      stop
      end
      subroutine test01 ( )

c*********************************************************************72
c
cc TEST01 demonstrates the use of ALNGAM.
c
c  Modified:
c
c    04 January 2006
c
c  Author:
c
c    John Burkardt
c
      implicit none

      double precision alngam
      double precision fx
      double precision fx2
      integer ifault
      integer n_data
      double precision x

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST01:'
      write ( *, '(a)' ) '  ALNGAM computes the logarithm of the '
      write ( *, '(a)' ) '  Gamma function.  We compare the result'
      write ( *, '(a)' ) '  to tabulated values.'
      write ( *, '(a)' ) ' '
      write ( *, '(a,a)' ) '          X                     ',
     &  'FX                        FX2'
      write ( *, '(a,a)' ) '                                ',
     &  '(Tabulated)               (ALNGAM)               DIFF'
      write ( *, '(a)' ) ' '

      n_data = 0

10    continue

        call gamma_log_values ( n_data, x, fx )

        if ( n_data .eq. 0 ) then
          go to 20
        end if

        fx2 = alngam ( x, ifault )

        write ( *, '(2x,f24.16,2x,g24.16,2x,g24.16,2x,g10.4)' )
     &  x, fx, fx2, dabs ( fx - fx2 )

      go to 10

20    continue

      return
      end
      subroutine test02 ( )

c*********************************************************************72
c
cc TEST02 demonstrates the use of LNGAMMA.
c
c  Modified:
c
c    04 January 2006
c
c  Author:
c
c    John Burkardt
c
      implicit none

      double precision fx
      double precision fx2
      integer ier
      double precision lngamma
      integer n_data
      double precision x

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST02:'
      write ( *, '(a)' ) '  LNGAMMA computes the logarithm of the '
      write ( *, '(a)' ) '  Gamma function.  We compare the result'
      write ( *, '(a)' ) '  to tabulated values.'
      write ( *, '(a)' ) ' '
      write ( *, '(a,a)' ) '          X                     ',
     &  'FX                        FX2'
      write ( *, '(a,a)' ) '                                ',
     &  '(Tabulated)               (LNGAMMA)               DIFF'
      write ( *, '(a)' ) ' '

      n_data = 0

10    continue

        call gamma_log_values ( n_data, x, fx )

        if ( n_data .eq. 0 ) then
          go to 20
        end if

        fx2 = lngamma ( x, ier )

        write ( *, '(2x,f24.16,2x,g24.16,2x,g24.16,2x,g10.4)' )
     &  x, fx, fx2, dabs ( fx - fx2 )

      go to 10

20    continue

      return
      end
