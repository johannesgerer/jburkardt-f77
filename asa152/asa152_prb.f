      program main

c*********************************************************************72
c
cc MAIN is the main program for ASA152_PRB.
c
c  Discussion:
c
c    ASA152_PRB calls the ASA152 routines.
c
c  Modified:
c
c    08 January 2008
c
c  Author:
c
c    John Burkardt
c
      implicit none

      call timestamp ( )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'ASA152_PRB:'
      write ( *, '(a)' ) '  FORTRAN77 version'
      write ( *, '(a)' ) '  Test the ASA152 library.'

      call test01 ( )
      call test02 ( )
c
c  Terminate.
c
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'ASA152_PRB:'
      write ( *, '(a)' ) '  Normal end of execution.'

      write ( *, '(a)' ) ' '
      call timestamp ( )

      stop
      end
      subroutine test01 ( )

c*********************************************************************72
c
cc TEST01 demonstrates CHYPER for cumulative probabilities.
c
c  Modified:
c
c    27 January 2008
c
c  Author:
c
c    John Burkardt
c
      implicit none

      double precision chyper
      double precision fx
      double precision fx2
      integer ifault
      integer n_data
      logical point
      integer pop
      integer sam
      integer suc
      integer x

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST01:'
      write ( *, '(a)' ) '  CHYPER computes cumulative probabilities'
      write ( *, '(a)' ) '  of the hypergeometric PDF.'
      write ( *, '(a)' ) '  Compare with tabulated values.'
      write ( *, '(a)' ) ' '
      write ( *, '(a,a)' ) '   SAM   SUC   POP     X    ',
     &  '  CDF                       CDF                     DIFF'
      write ( *, '(a,a)' ) '                            ',
     &  ' (tabulated)               (CHYPER)'
      write ( *, '(a)' ) ' '

      n_data = 0

10    continue

        call hypergeometric_cdf_values ( n_data, sam, suc, pop, x, fx )

        if ( n_data .eq. 0 ) then
          go to 20
        end if

        point = .false.
        fx2 = chyper ( point, sam, x, pop, suc, ifault )

        write ( *,
     &  '(2x,i4,2x,i4,2x,i4,2x,i4,2x,g24.16,2x,g24.16,2x,g10.4)' )
     &  sam, suc, pop, x, fx, fx2, dabs ( fx - fx2 )

      go to 10

20    continue

      return
      end
      subroutine test02 ( )

c*********************************************************************72
c
cc TEST02 demonstrates CHYPER for point probabilities.
c
c  Modified:
c
c    08 January 2008
c
c  Author:
c
c    John Burkardt
c
      implicit none

      double precision chyper
      double precision fx
      double precision fx2
      integer ifault
      integer n_data
      logical point
      integer pop
      integer sam
      integer suc
      integer x

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST02:'
      write ( *, '(a)' ) '  CHYPER computes point probabilities'
      write ( *, '(a)' ) '  of the hypergeometric PDF.'
      write ( *, '(a)' ) '  Compare with tabulated values.'
      write ( *, '(a)' ) ' '
      write ( *, '(a,a)' ) '   SAM   SUC   POP     X',
     &  '      PDF                       PDF                    DIFF'
      write ( *, '(a,a)' ) '                        ',
     &  '     (tabulated)                (CHYPER)'
      write ( *, '(a)' ) ' '

      n_data = 0

10    continue

        call hypergeometric_pdf_values ( n_data, sam, suc, pop, x, fx )

        if ( n_data .eq. 0 ) then
          go to 20
        end if

        point = .true.
        fx2 = chyper ( point, sam, x, pop, suc, ifault )

        write ( *,
     &  '(2x,i4,2x,i4,2x,i4,2x,i4,2x,g24.16,2x,g24.16,2x,g10.4)' )
     &  sam, suc, pop, x, fx, fx2, dabs ( fx - fx2 )

      go to 10

20    continue

      return
      end
