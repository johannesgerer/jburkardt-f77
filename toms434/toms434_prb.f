      program main

c*********************************************************************72
c
cc TOMS434_PRB tests CONP.
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

      integer c
      integer r
      parameter ( c = 3 )
c     parameter ( r = 3 )
      parameter ( r = 2 )

      integer nc
      integer nr
      parameter ( nc = c + 1 )
      parameter ( nr = r + 1 )

      integer i
      integer j
      integer matrix(nr,nc)
      real pc
      real ps
      real pt

      save matrix
c
c  Matrix data is by COLUMNS.
c
      data matrix /
     &   1,  1,  2,
     &   2,  0,  2,
     &   3,  2,  5,
     &   6,  3,  9 /

c    &  14, 10,  5, 29,
c    &   2,  4, 10, 16,
c    &   4,  6,  5, 15,
c    &  20, 20, 20, 60 /

      call timestamp ( )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TOMS434_PRB'
      write ( *, '(a)' ) '  FORTRAN77 version'
      write ( *, '(a)' ) '  Test the TOMS434 library.'

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  The contingency table:'
      write ( *, '(a)' ) ' '
      do i = 1, nr
        write ( *, '(6(2x,i4))' ) ( matrix(i,j), j = 1, nc )
      end do

      call conp ( matrix, nr, nc, pt, ps, pc )

      write ( *, '(a)' ) ' '
      write ( *, '(a,g14.6)' ) '  PT, prob of this table ', pt
      write ( *, '(a,g14.6)' )
     &  '  PS, prob of this or less likely table ', ps
      write ( *, '(a,g14.6)' ) '  PC, prob of some table ', pc
c
c  Terminate.
c
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TOMS434_PRB'
      write ( *, '(a)' ) '  Normal end of execution.'

      write ( *, '(a)' ) ' '
      call timestamp ( )

      stop
      end
