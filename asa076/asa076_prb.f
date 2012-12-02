      program main

c*********************************************************************72
c
cc MAIN is the main program for ASA076_PRB.
c
c  Discussion:
c
c    ASA076_PRB calls the ASA076 routines.
c
c  Modified:
c
c    16 January 2008
c
c  Author:
c
c    John Burkardt
c
      implicit none

      call timestamp ( )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'ASA076_PRB:'
      write ( *, '(a)' ) '  FORTRAN77 version'
      write ( *, '(a)' ) '  Test the ASA076 library.'

      call test01 ( )
      call test02 ( )
c
c  Terminate.
c
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'ASA076_PRB:'
      write ( *, '(a)' ) '  Normal end of execution.'

      write ( *, '(a)' ) ' '
      call timestamp ( )

      stop
      end
      subroutine test01 ( )

c*********************************************************************72
c
cc TEST01 demonstrates the use of TFN.
c
c  Modified:
c
c    16 January 2008
c
c  Author:
c
c    John Burkardt
c
      implicit none

      double precision a
      double precision h
      integer n_data
      double precision t1
      double precision t2
      double precision tfn

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST01:'
      write ( *, '(a)' ) '  TFN computes the Owen T function.'
      write ( *, '(a)' ) '  Compare with tabulated values.'
      write ( *, '(a)' ) ' '
      write ( *, '(a,a)' ) '          H               A           ',
     &'T                         T                       DIFF'
      write ( *, '(a,a)' ) '                                     ',
     &'(Tabulated)               (TFN)'
      write ( *, '(a)' ) ' '

      n_data = 0

10    continue

        call owen_values ( n_data, h, a, t1 )

        if ( n_data .eq. 0 ) then
          go to 20
        end if

        t2 = tfn ( h, a )

        write ( *, '(2x,f14.6,2x,f14.6,2x,g24.16,2x,g24.16,2x,g10.4)' )
     &  h, a, t1, t2, dabs ( t1 - t2 )

      go to 10

20    continue

      return
      end
      subroutine test02 ( )

c*********************************************************************72
c
cc TEST02 demonstrates the use of THA.
c
c  Modified:
c
c    16 January 2008
c
c  Author:
c
c    John Burkardt
c
      implicit none

      double precision a
      double precision h
      integer n_data
      double precision t1
      double precision t2
      double precision tha

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST02:'
      write ( *, '(a)' ) '  THA computes the Owen T function.'
      write ( *, '(a)' ) '  Compare with tabulated values.'
      write ( *, '(a)' ) ' '
      write ( *, '(a,a)' ) '          H               A           ',
     &'T                         T                       DIFF'
      write ( *, '(a,a)' ) '                                     ',
     &'(Tabulated)               (THA)'
      write ( *, '(a)' ) ' '

      n_data = 0

10    continue

        call owen_values ( n_data, h, a, t1 )

        if ( n_data .eq. 0 ) then
          go to 20
        end if

        t2 = tha ( h, 1.0D+00, a, 1.0D+00 )

        write ( *, '(2x,f14.6,2x,f14.6,2x,g24.16,2x,g24.16,2x,g10.4)' )
     &  h, a, t1, t2, dabs ( t1 - t2 )

      go to 10

20    continue

      return
      end
