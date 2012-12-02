      program main

c*********************************************************************72
c
cc MAIN is the main program for ASA063_PRB.
c
c  Discussion:
c
c    ASA063_PRB calls the ASA063 routines.
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
      write ( *, '(a)' ) 'ASA063_PRB:'
      write ( *, '(a)' ) '  FORTRAN77 version'
      write ( *, '(a)' ) '  Test the ASA063 library.'

      call test01 ( )
c
c  Terminate.
c
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'ASA063_PRB:'
      write ( *, '(a)' ) '  Normal end of execution.'

      write ( *, '(a)' ) ' '
      call timestamp ( )

      stop
      end
      subroutine test01 ( )

c*********************************************************************72
c
cc TEST01 demonstrates the use of BETAIN.
c
c  Modified:
c
c    10 January 2006
c
c  Author:
c
c    John Burkardt
c
      implicit none

      double precision a
      double precision alogam
      double precision b
      double precision beta_log
      double precision betain
      double precision fx
      double precision fx2
      integer ifault
      integer n_data
      double precision x

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST01:'
      write ( *, '(a)' )
     &  '  BETAIN evaluates the incomplete Beta function.'
      write ( *, '(a)' ) '  Compare with tabulated values.'
      write ( *, '(a)' ) ' '
      write ( *, '(a,a)' ) '      A       B       X      ',
     &  ' BetaInc                   BetaInc                 DIFF'
      write ( *, '(a,a)' ) '                             ',
     &  '(tabulated)               (BETAIN)'
      write ( *, '(a)' ) ' '

      n_data = 0

10    continue

        call beta_inc_values ( n_data, a, b, x, fx )

        if ( n_data .eq. 0 ) then
          go to 20
        end if

        beta_log = alogam ( a, ifault )
     &           + alogam ( b, ifault )
     &           - alogam ( a + b, ifault )

        fx2 = betain ( x, a, b, beta_log, ifault )

        write ( *,
     &  '(2x,f6.2,2x,f6.2,2x,f6.2,2x,g24.16,2x,g24.16,2x,g10.4)' )
     &  a, b, x, fx, fx2, dabs ( fx - fx2 )

      go to 10

20    continue

      return
      end
