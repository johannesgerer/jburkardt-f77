      program main

!*********************************************************************72
!
!! MAIN is the main program for SPARSEKIT_PRB06.
!
!  Discussion:
!
!    SPARSEKIT_PRB06 demonstrates how to read a Harwell-Boeing 
!    sparse matrix file.
!
      implicit none

      integer nmax
      parameter ( nmax = 500 )
      integer nzmax
      parameter ( nzmax = 7000 )

      double precision a(nzmax)
      double precision a1(nzmax)
      character * ( 80 ) filnam
      character * ( 2 ) guesol
      integer ia(nmax+1)
      integer ia1(nmax+1)
      integer ierr
      integer iin
      integer ios
      integer iout
      integer ja(nzmax)
      integer ja1(nzmax)
      integer job
      character * ( 8 ) key
      integer ncol
      integer nnz
      integer nrhs
      integer nrow
      double precision rhs(1)
      character * ( 72 ) title
      character * ( 3 ) type
      logical valued

      iout = 6
      job = 2
      nrhs = 0

      call timestamp ( )

      write ( *, * ) ' '
      write ( *, * ) 'SPARSEKIT_PRB06'
      write ( *, '(a)' ) '  FORTRAN77 version'
      write ( *, '(a)' ) ' '
      write ( *, * ) '  This program demonstrates the use of the'
      write ( *, * ) '  routines READMT and DINFO1 to read and report'
      write ( *, * ) '  on a sparse matrix stored in a file in the'
      write ( *, * ) '  format used by the Harwell-Boeing Sparse Matrix'
      write ( *, * ) '  Collection or "HBSMC".'
      write ( *, * ) ' '

      filnam = 'saylor_hb.txt'
      iin = 20

      open ( unit = iin, file = filnam, status = 'old', err = 99 )

      call readmt ( nmax, nzmax, job, iin, a, ja, ia, rhs, nrhs, 
     &  guesol, nrow, ncol, nnz, title, key, type, ierr )

      close ( unit = iin )
!
!  If not readable, return.
!
      if ( ierr .ne. 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'Error!'
        write ( *, '(a)' ) '  Unable to read matrix.'
        write ( *, '(a,i6)' ) '  READMT returned IERR = ', ierr
        stop
      end if

      valued = ( 2 <= job )

      call dinfo1 ( ncol, iout, a, ja, ia, valued, title, key, type, 
     &  a1, ja1, ia1 )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'SPARSEKIT_PRB06'
      write ( *, '(a)' ) '  Normal end of execution.'

      write ( *, '(a)' ) ' '
      call timestamp ( )

      stop

99    continue

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'SPARSEKIT_PRB06 - Error!'
      write ( *, '(a)' ) '  Unable to open file.'
      write ( *, '(a)' ) '  Abnormal end of execution.'
      write ( *, '(a)' ) ' '
      call timestamp ( )
      stop
      end
