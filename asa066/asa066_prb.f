      program main

c*********************************************************************72
c
cc MAIN is the main program for ASA066_PRB.
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

      write ( *, '(a)' ) ' '
      call timestamp ( )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'ASA066_PRB:'
      write ( *, '(a)' ) '  FORTRAN77 version'
      write ( *, '(a)' ) '  Test the ASA066 library.'

      call test01 ( )
      call test02 ( )
      call test03 ( )
c
c  Terminate.
c
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'ASA066_PRB:'
      write ( *, '(a)' ) '  Normal end of execution.'

      write ( *, '(a)' ) ' '
      call timestamp ( )

      stop
      end
      subroutine test01 ( )

c*********************************************************************72
c
cc TEST01 compares ALNORM against tabulated values.
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

      double precision alnorm
      double precision fx
      double precision fx2
      integer n_data
      logical upper
      double precision x

      upper = .false.

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST01:'
      write ( *, '(a)' ) '  Compare tabulated values of the normal '
      write ( *, '(a)' ) '  Cumulative Density Function against values'
      write ( *, '(a)' ) '  computed by ALNORM.'
      write ( *, '(a)' ) ' '
      write ( *, '(a,a)' )
     &  '         X        CDF                       CDF',
     &  '                    DIFF'
      write ( *, '(a)' )
     &  '               (tabulated)                 (ALNORM)'
      write ( *, '(a)' ) ' '

      n_data = 0

10    continue

        call normal_01_cdf_values ( n_data, x, fx )

        if ( n_data .eq. 0 ) then
          go to 20
        end if

        fx2 = alnorm ( x, upper )

        write ( *, '(2x,f10.4,2x,g24.16,2x,g24.16,2x,g10.4)' )
     &  x, fx, fx2, dabs ( fx - fx2 )

      go to 10

20    continue

      return
      end
      subroutine test02 ( )

c*********************************************************************72
c
cc TEST02 compares NORMP against tabulated values.
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

      double precision fx
      double precision fx2
      double precision fx3
      integer n_data
      double precision pdf
      logical upper
      double precision x

      upper = .false.

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST02:'
      write ( *, '(a)' ) '  Compare tabulated values of the normal '
      write ( *, '(a)' ) '  Cumulative Density Function against values'
      write ( *, '(a)' ) '  computed by NORMP.'
      write ( *, '(a)' ) ' '
      write ( *, '(a,a)' )
     &  '         X        CDF                       CDF',
     &  '                    DIFF'
      write ( *, '(a)' )
     &  '               (tabulated)                 (NORMP)'
      write ( *, '(a)' ) ' '

      n_data = 0

10    continue

        call normal_01_cdf_values ( n_data, x, fx )

        if ( n_data .eq. 0 ) then
          go to 20
        end if

        call normp ( x, fx2, fx3, pdf )

        write ( *, '(2x,f10.4,2x,g24.16,2x,g24.16,2x,g10.4)' )
     &  x, fx, fx2, dabs ( fx - fx2 )


      go to 10

20    continue

      return
      end
      subroutine test03 ( )

c*********************************************************************72
c
cc TEST03 compares NPROB against tabulated values.
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

      double precision fx
      double precision fx2
      double precision fx3
      integer n_data
      double precision pdf
      logical upper
      double precision x

      upper = .false.

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST03:'
      write ( *, '(a)' ) '  Compare tabulated values of the normal '
      write ( *, '(a)' ) '  Cumulative Density Function against values'
      write ( *, '(a)' ) '  computed by NPROB.'
      write ( *, '(a)' ) ' '
      write ( *, '(a,a)' )
     &  '         X        CDF                       CDF',
     &  '                    DIFF'
      write ( *, '(a)' )
     &  '               (tabulated)                 (NPROB)'
      write ( *, '(a)' ) ' '

      n_data = 0

10    continue

        call normal_01_cdf_values ( n_data, x, fx )

        if ( n_data .eq. 0 ) then
          go to 20
        end if

        call nprob ( x, fx2, fx3, pdf )

        write ( *, '(2x,f10.4,2x,g24.16,2x,g24.16,2x,g10.4)' )
     &  x, fx, fx2, dabs ( fx - fx2 )


      go to 10

20    continue

      return
      end
