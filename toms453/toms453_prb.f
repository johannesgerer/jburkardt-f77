      program main

c*********************************************************************72
c
cc TOMS453_PRB tests BROMIN, ACM TOMS algorithm 453.
c
c  Modified:
c
c    12 January 2006
c
c  Author:
c
c    John Burkardt
c
      implicit none

      call timestamp ( )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TOMS453_PRB'
      write ( *, '(a)' ) '  FORTRAN77 version'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Test the TOMS453 library.'

      call test01 ( )
c
c  Terminate.
c
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TOMS453_PRB'
      write ( *, '(a)' ) '  Normal end of execution.'

      write ( *, '(a)' ) ' '
      call timestamp ( )

      stop
      end
      subroutine test01 ( )

c*********************************************************************72
c
cc TEST01 tests BROMIN, ACM TOMS algorithm 453.
c
c  Modified:
c
c    11 July 2008
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer nhalf_max
      integer n_num
      integer s_num

      parameter ( nhalf_max = 6 )
      parameter ( n_num = 4 )
      parameter ( s_num = 4 )

      double precision eps
      integer i
      integer ier
      integer j
      integer k
      integer kk
      integer n
      integer n_half
      integer n_vec(n_num)
      double precision s
      double precision s_vec(s_num)
      double precision tol
      double precision wi(nhalf_max)
      double precision wr(nhalf_max)
      double precision xi(nhalf_max)
      double precision xr(nhalf_max)
      double precision total

      save n_vec
      save s_vec

      data n_vec / 3, 6, 9, 12 /
      data s_vec / 0.0D+00, 0.1D+00, 1.0D+00, 4.0D+00 /

      tol = 0.1D-8

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST01'
      write ( *, '(a)' ) '  Determine abscissas and weights for'
      write ( *, '(a)' ) '  a variety of values of S and N.'

      do i = 1, n_num

        n = n_vec(i)
        n_half = ( n + 1 ) / 2

        do j = 1, s_num

          s = s_vec(j)

          write ( *, '(a)' ) ' '
          write ( *, '(a,i6)' ) '  N = ', n
          write ( *, '(a,g14.6)' ) '  S = ', s

          call bromin ( n, s, tol, xr, xi, wr, wi, eps, ier )

          if ( 0 .lt. ier ) then

            write ( *, '(a)' ) ' '
            write ( *, '(a,i6)' ) 'BROMIN returned IER = ', ier

          else

            if ( ier .eq. -1 ) then
              write ( *, '(a)' ) ' '
              write ( *, '(a)' ) '  Note that the requested accuracy'
              write ( *, '(a)' ) '  was not achieved.'
            end if

            write ( *, '(a)' ) ' '
            write ( *, '(a,a)' ) '                           ',
     &        'XR              XI              WR              WI'
            write ( *, '(a)' ) ' '
            total = 0.0D+00
            do kk = 1, n
              if ( kk .le. ( n - n_half ) ) then
                k = n_half + 1 - kk
                write ( *, '(2x,i8,2x,i8,4(2x,g14.6))' )
     &            kk, k, xr(k), - xi(k), wr(k), - wi(k)
                total = total + wr(k)
              else
                k = kk - ( n - n_half )
                write ( *, '(2x,i8,2x,i8,4(2x,g14.6))' )
     &            kk, k, xr(k),   xi(k), wr(k),   wi(k)
                total = total + wr(k)
              end if
            end do
            write ( *, '(2x,a8,2x,8x,16x,16x,2x,g14.6)' )
     &        'WR total', total

          end if

        end do

      end do

      return
      end
