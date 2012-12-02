      program main

c*********************************************************************72
c
cc MAIN is the main program for OWENS_PRB.
c
c  Discussion:
c
c    OWENS_PRB calls the OWENS routines.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    23 November 2010
c
c  Author:
c
c    John Burkardt
c
      implicit none

      call timestamp ( )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'OWENS_PRB:'
      write ( *, '(a)' ) '  FORTRAN77 version'
      write ( *, '(a)' ) '  Test the OWENS library.'

      call test01 ( )
      call test02 ( )
      call test03 ( )
      call test04 ( )
c
c  Terminate.
c
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'OWENS_PRB:'
      write ( *, '(a)' ) '  Normal end of execution.'

      write ( *, '(a)' ) ' '
      call timestamp ( )

      stop
      end
      subroutine test01 ( )

c*********************************************************************72
c
cc TEST01 demonstrates the use of T.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
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
      double precision t
      double precision t1
      double precision t2

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST01:'
      write ( *, '(a)' ) '  T computes the Owen T function.'
      write ( *, '(a)' ) '  Compare with tabulated values.'
      write ( *, '(a)' ) ' '
      write ( *, '(a,a)' ) '          H               A           ',
     &  'T                         T                       DIFF'
      write ( *, '(a,a)' ) '                                     ',
     &  '(Tabulated)               (TFN)'
      write ( *, '(a)' ) ' '

      n_data = 0

10    continue

        call owen_values ( n_data, h, a, t1 )

        if ( n_data .eq. 0 ) then
          go to 20
        end if

        t2 = t ( h, a )

        write ( *, '(2x,f14.6,2x,f14.6,2x,g24.16,2x,g24.16,2x,g10.4)' )
     &    h, a, t1, t2, dabs ( t1 - t2 )

      go to 10

20    continue

      return
      end
      subroutine test02 ( )

c*********************************************************************72
c
cc TEST02 demonstrates the use of BIVNOR.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    23 May 2009
c
c  Author:
c
c    John Burkardt
c
      implicit none

      double precision bivnor
      double precision fxy1
      double precision fxy2
      integer n_data
      double precision r
      double precision x
      double precision y

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST02:'
      write ( *, '(a)' )
     &  '  BIVNOR computes the bivariate normal probability.'
      write ( *, '(a)' ) '  Compare with tabulated values.'
      write ( *, '(a)' ) ' '
      write ( *, '(a,a,a)' )
     &  '          X               Y               ',
     &  'R           P                         P',
     &  '                       DIFF'
      write ( *, '(a,a)' )
     &  '                                          ',
     &  '           (Tabulated)               (BIVNOR)'
      write ( *, '(a)' ) ' '

      n_data = 0

10    continue

        call bivariate_normal_cdf_values ( n_data, x, y, r, fxy1 )

        if ( n_data .eq. 0 ) then
          go to 20
        end if
c
c  BIVNOR computes the "tail" of the probability, and we want the
c  initial part!  To get that value, negate X and Y.
c
        fxy2 = bivnor ( - x, - y, r )

        write ( *,
     &    '(2x,f14.6,2x,f14.6,2x,f14.6,2x,g24.16,2x,g24.16,2x,g10.4)' )
     &    x, y, r, fxy1, fxy2, dabs ( fxy1 - fxy2 )

      go to 10

20    continue

      return
      end
      subroutine test03 ( )

c*********************************************************************72
c
cc TEST03 demonstrates the use of ZNORM1.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    24 May 2009
c
c  Author:
c
c    John Burkardt
c
      implicit none

      double precision fx1
      double precision fx2
      integer n_data
      double precision x
      double precision znorm1

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST03:'
      write ( *, '(a)' )
     &  '  ZNORM1 computes the normal CDF starting at 0.'
      write ( *, '(a)' ) '  Compare with tabulated values.'
      write ( *, '(a)' ) ' '
      write ( *, '(a,a)' )
     &  '          X           P                         ',
     &  'P                       DIFF',
     &  '                       DIFF'
      write ( *, '(a)' )
     &  '                     (Tabulated)               (ZNORM1)'
      write ( *, '(a)' ) ' '

      n_data = 0

10    continue

        call normal_01_cdf_values ( n_data, x, fx1 )

        if ( n_data .eq. 0 ) then
          go to 20
        end if

        fx1 = fx1 - 0.5D+00

        fx2 = znorm1 ( x )

        write ( *,
     &    '(2x,f14.6,2x,g24.16,2x,g24.16,2x,g10.4)' )
     &    x, fx1, fx2, dabs ( fx1 - fx2 )

      go to 10

20    continue

      return
      end
      subroutine test04 ( )

c*********************************************************************72
c
cc TEST04 demonstrates the use of ZNORM2.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    23 November 2010
c
c  Author:
c
c    John Burkardt
c
      implicit none

      double precision fx1
      double precision fx2
      integer n_data
      double precision x
      double precision znorm2

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST04:'
      write ( *, '(a)' )
     &  '  ZNORM2 computes the complementary normal CDF.'
      write ( *, '(a)' ) '  Compare with tabulated values.'
      write ( *, '(a)' ) ' '
      write ( *, '(a,a)' )
     &  '          X           P                         ',
     &  'P                       DIFF',
     &  '                       DIFF'
      write ( *, '(a)' )
     &  '                     (Tabulated)               (ZNORM2)'
      write ( *, '(a)' ) ' '

      n_data = 0

10    continue

        call normal_01_cdf_values ( n_data, x, fx1 )

        if ( n_data .eq. 0 ) then
          go to 20
        end if

        fx1 = 1.0D+00 - fx1

        fx2 = znorm2 ( x )

        write ( *,
     &    '(2x,f14.6,2x,g24.16,2x,g24.16,2x,g10.4)' )
     &    x, fx1, fx2, dabs ( fx1 - fx2 )

      go to 10

20    continue

      return
      end
