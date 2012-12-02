      program main

c*********************************************************************72
c
cc MAIN is the main program for TEST_ZERO_PRB.
c
c  Discussion:
c
c    TEST_ZERO_PRB demonstrates the use of the TEST_ZERO scalar test functions.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    16 May 2011
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer max_root
      parameter ( max_root = 4 )
      integer max_start
      parameter ( max_start = 4 )

      double precision fatol
      parameter ( fatol = 1.0D-06 )
      double precision fx
      double precision fxa
      double precision fxb
      double precision fxc
      integer i
      integer max_step
      parameter ( max_step = 25 )
      integer prob
      integer prob_num
      double precision r8_sign
      double precision range(2)
      integer root_num
      integer start_num
      character*80 title
      double precision x
      double precision xa
      double precision xatol
      parameter ( xatol = 1.0D-06 )
      double precision xb
      double precision xc
      double precision xmax
      double precision xmin
      double precision xrtol
      parameter ( xrtol = 1.0D-06 )
      double precision xstart(max_start)

      call timestamp ( )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST_ZERO_PRB'
      write ( *, '(a)' ) '  FORTRAN77 version'
      write ( *, '(a)' ) '  Test the TEST_ZERO library.'
      write ( *, '(a)' ) ' '
      write ( *, '(a,g14.6)' ) '  Function value tolerance = ', fatol
      write ( *, '(a,g14.6)' ) '  Root absolute tolerance =  ', xatol
      write ( *, '(a,g14.6)' ) '  Root relative tolerance =  ', xrtol
      write ( *, '(a,i4)' ) '  Maximum number of steps =  ', max_step
c
c  Find out how many problems there are
c
      call p00_prob_num ( prob_num )
      write ( *, '(a)' ) ' '
      write ( *, '(a,i4)' ) '  Number of problems is ', prob_num

      do prob = 1, prob_num
c
c  Get the problem title.
c
        call p00_title ( prob, title )
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) ' '
        write ( *, '(a,i4)' ) '  Problem number ', prob
        write ( *, '(2x,a)' ) trim ( title )
c
c  Get the problem interval.
c
        call p00_range ( prob, range )
        xmin = range(1)
        xmax = range(2)
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) '  We seek roots between'
        write ( *, '(2x,g14.6)' ) xmin
        write ( *, '(a)' ) '  and'
        write ( *, '(2x,g14.6)' ) xmax
c
c  Get the number of roots.
c
        call p00_root_num ( prob, root_num )

        write ( *, '(a)' ) ' '
        write ( *, '(a,i4)' ) '  Number of known roots = ', root_num
c
c  Print the exact solution.
c
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) '  Tabulated solutions:'
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) '        X            F(X)'
        write ( *, '(a)' ) ' '
        do i = 1, root_num
          call p00_root ( prob, i, x )
          call p00_fx ( prob, x, fx )
          write ( *, '(2x,2g16.8)' ) x, fx
        end do
c
c  Get the number of starting points.
c
        call p00_start_num ( prob, start_num )

        write ( *, '(a)' ) ' '
        write ( *, '(a,i4)' ) 
     &    '  Number of starting points = ', start_num
c
c  Get the starting points.
c
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) '  I  XSTART(I), F(XSTART(I))'
        write ( *, '(a)' ) ' '
        do i = 1, start_num
          call p00_start ( prob, i, xstart(i) )
          call p00_fx ( prob, xstart(i), fx )
          write ( *, '(2x,i2,2x,2g16.8)' ) i, xstart(i), fx
        end do
c
c  Bisection.
c
        call p00_start ( prob, 1, xa )
        call p00_fx ( prob, xa, fxa )

        do i = 2, start_num
          call p00_start ( prob, i, xb )
          call p00_fx ( prob, xb, fxb )
          if ( r8_sign ( fxa ) .ne. r8_sign ( fxb ) ) then
            call bisection ( fatol, max_step, prob, xatol, xa, xb, 
     &        fxa, fxb )
            go to 10
          end if

        end do

10      continue
c
c  Brent's method.
c
        call p00_start ( prob, 1, xa )
        call p00_fx ( prob, xa, fxa )

        do i = 2, start_num
          call p00_start ( prob, i, xb )
          call p00_fx ( prob, xb, fxb )
          if ( r8_sign ( fxa ) .ne. r8_sign ( fxb ) ) then
            call brent ( fatol, max_step, prob, xatol, xrtol, xa, xb, 
     &        fxa, fxb )
            go to 20
          end if

        end do

20      continue
c
c  Muller.
c
        if ( 3 <= start_num ) then

          call p00_start ( prob, 1, xa )
          call p00_fx ( prob, xa, fxa )
          call p00_start ( prob, 2, xb )
          call p00_fx ( prob, xb, fxb )
          call p00_start ( prob, 3, xc )
          call p00_fx ( prob, xc, fxc )

          call muller ( fatol, max_step, prob, xatol, xrtol, xa, xb, xc,
     &      fxa, fxb, fxc )

        end if
c
c  Newton's method.
c
        do i = 1, start_num

          call p00_start ( prob, i, xa )
          call p00_fx ( prob, xa, fxa )

          call newton ( fatol, max_step, prob, xatol, xmin, xmax, 
     &      xa, fxa )

        end do
c
c  Regula Falsi.
c
        call p00_start ( prob, 1, xa )
        call p00_fx ( prob, xa, fxa )

        do i = 2, start_num

          call p00_start ( prob, i, xb )
          call p00_fx ( prob, xb, fxb )
          if ( r8_sign ( fxa ) .ne. r8_sign ( fxb ) ) then
            call regula_falsi ( fatol, max_step, prob, xatol, xa, xb, 
     &        fxa, fxb )
            go to 30
          end if

        end do

30      continue
c
c  Secant.
c
        do i = 1, start_num - 1

          call p00_start ( prob, i, xa )
          call p00_fx ( prob, xa, fxa )

          call p00_start ( prob, i + 1, xb )
          call p00_fx ( prob, xb, fxb )

          call secant ( fatol, max_step, prob, xatol, xmin, xmax, xa, 
     &      xb, fxa, fxb )

        end do

      end do
c
c  Terminate.
c
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST_ZERO_PRB'
      write ( *, '(a)' ) '  Normal end of execution.'

      write ( *, '(a)' ) ' '
      call timestamp ( )

      stop
      end
