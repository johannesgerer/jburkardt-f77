      program main

c*********************************************************************72
c
cc MAIN is the main program for BETA_NC_PRB.
c
c  Discussion:
c
c    BETA_NC_PRB calls the BETA_NC routines.
c
c  Modified:
c
c    29 January 2008
c
c  Author:
c
c    John Burkardt
c
      implicit none

      call timestamp ( )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'BETA_NC_PRB:'
      write ( *, '(a)' ) '  FORTRAN77 version'
      write ( *, '(a)' ) '  Test the BETA_NC library.'

      call test01 ( )
c
c  Terminate.
c
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'BETA_NC_PRB:'
      write ( *, '(a)' ) '  Normal end of execution.'

      write ( *, '(a)' ) ' '
      call timestamp ( )

      stop
      end
      subroutine test01 ( )

c*********************************************************************72
c
cc TEST01 tests BETA_NONCENTRAL_CDF against tabulated values.
c
c  Modified:
c
c    29 January 2008
c
c  Author:
c
c    John Burkardt
c
      implicit none

      double precision a
      double precision b
      double precision error_max
      double precision fx
      double precision fx2
      double precision lambda
      integer n_data
      double precision x

      error_max = 1.0D-10

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST01:'
      write ( *, '(a)' ) '  Compare tabulated values of the noncentral'
      write ( *, '(a)' ) '  incomplete Beta Function against values'
      write ( *, '(a)' ) '  computed by BETA_NONCENTRAL_CDF.'
      write ( *, '(a)' ) ' '
      write ( *, '(a,a)' )
     &  '      A        B     LAMBDA        X       ',
     &  ' CDF                         CDF                    DIFF'
      write ( *, '(a,a)' )
     &  '                                           ',
     &  '(tabulated)                 (BETA_NC)'
      write ( *, '(a)' ) ' '

      n_data = 0

10    continue

        call beta_noncentral_cdf_values ( n_data, a, b, lambda,
     &  x, fx )

        if ( n_data .eq. 0 ) then
          go to 20
        end if

        call beta_noncentral_cdf ( a, b, lambda, x, error_max, fx2 )

        write ( *,
     &  '(2x,f7.1,2x,f7.1,2x,f7.1,2x,f10.4,
     &  2x,g24.6,2x,g24.6,2x,g10.4)' )
     &  a, b, lambda, x, fx, fx2, dabs ( fx - fx2 )

      go to 10

20    continue

      return
      end
