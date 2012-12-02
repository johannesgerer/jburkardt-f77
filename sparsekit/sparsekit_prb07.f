      program main

c*********************************************************************72
c
cc MAIN generates a Markov chain matrix to test eigenvalue routines.
c
c  The matrix A produced by this routine can also be used to generate
c  a related singular matrix, I-A.
c  The matrix models  simple random walk on a triangular grid.
c  see additional comments in subroutine.
c
c   will create a matrix in the HARWELL/BOEING format and put it in
c   the file markov.mat
c
      implicit none

      integer nmax
      parameter ( nmax = 5000 )
      integer nzmax
      parameter ( nzmax= 4 * nmax )

      double precision a(nzmax)
      integer ia(nmax+1)
      integer ifmt
      integer ios
      integer iout
      integer ja(nzmax)
      integer job
      character * ( 8 ) key
      integer m
      integer n
      character * ( 72 ) title
      character * ( 3 ) type
      double precision rhs(1)

      call timestamp ( )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'SPARSEKIT_PRB07'
      write ( *, '(a)' ) '  FORTRAN77 version'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Generate a Markov chain matrix to test the'
      write ( *, '(a)' ) '  eigenvalue routine.'

      open ( unit = 11, file = 'markov.mat', status = 'replace', 
     &  err = 99 )
c
c  Set grid size - will not accept too large grids.
c
      m = 5

      if ( 2 * nmax < m * ( m + 1 ) ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'SPARSEKIT_PRB07 - Fatal error!'
        write ( *, '(a)' ) 
     &    '  M is too large - unable to produce matrix.'
        stop
      end if
c
c  Call the matrix generator.
c
      call markgen ( m, n, a, ja, ia )
c
c  Store result in file.
c
      title = ' Test matrix from SPARSKIT - Markov chain model'
      key = 'randwk01'
      type = 'rua'
      iout = 11
      job = 2
      ifmt = 10

      call prtmt ( n, n, a, ja, ia, rhs, 'NN', title, key, 
     &  type, ifmt, job, iout )

      close ( unit = iout )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 
     &  '  The matrix has been stored in Harwell/Boeing format'
      write ( *, '(a)' ) '  in the file "markov.mat".'

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'SPARSEKIT_PRB07'
      write ( *, '(a)' ) '  Normal end of execution.'

      write ( *, '(a)' ) ' '
      call timestamp ( )

      stop

99    continue

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'SPARSEKIT_PRB07 - Fatal error!'
      write ( *, '(a)' ) '  Could not open the file.'
      write ( *, '(a)' ) '  Abnormal end of execution.'

      write ( *, '(a)' ) ' '
      call timestamp ( )
      stop

      end

