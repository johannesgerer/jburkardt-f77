      program main

c*********************************************************************72
c
cc MAIN is the main program for REVNEW_PRB.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    06 December 2008
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n
      parameter ( n = 4 )

      double precision fx(n)
      double precision fxtol
      integer i
      integer ido
      double precision x(n)

      call timestamp ( )
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'REVNEW_PRB'
      write ( *, '(a)' ) '  FORTRAN77 version'
      write ( *, '(a)' ) '  Test the REVNEW library.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  REVNEW solves a nonlinear equation'
      write ( *, '(a)' ) '  using reverse communication.'
c
c  Initialization.
c
      ido = 0
      fxtol = 0.000001D+00
      do i = 1, n
        x(i) = 0.0D+00
      end do

      do i = 1, n
        fx(i) = ( x(i) - dble ( i ) ) * ( x(i) - dble ( i ) )
      end do

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Initial Values:'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '     I      X               FX'
      write ( *, '(a)' ) ' '
      do i = 1, n
        write ( *, '(2x,i4,2x,g14.6,2x,g14.6)' ) i, x(i), fx(i)
      end do

      write ( *, '(a)' ) ' '

      fx(1) = fxtol
c
c  The solution loop.
c
10    continue

      call revnew ( fx, ido, n, x )

      if ( ido .eq. 0 ) then

        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) '  Convergence:'
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) '     I      X               FX'
        write ( *, '(a)' ) ' '
        do i = 1, n
          write ( *, '(2x,i4,2x,g14.6,2x,g14.6)' ) i, x(i), fx(i)
        end do

      else if ( ido .eq. 1 .or. ido .eq. 2 ) then

        do i = 1, n
          fx(i) = ( x(i) - dble ( i ) ) * ( x(i) - dble ( i ) )
        end do

        go to 10

      end if

20    continue
c
c  Terminate.
c
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'REVNEW_PRB:'
      write ( *, '(a)' ) '  Normal end of execution.'
      write ( *, '(a)' ) ' '
      call timestamp ( )

      stop
      end
