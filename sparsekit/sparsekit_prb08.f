      program main

c*********************************************************************72
c
cc MAIN is the main program for SPARSEKIT_PRB08.
c
c  Discussion:
c
c    SPARSEKIT_PRB08 runs the Zlatev test suite.
c
c    Three matrices are generated, written to separate files
c    in the Harwell-Boeing format.
c
      implicit none

      call timestamp ( )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'SPARSEKIT_PRB08'
      write ( *, '(a)' ) '  FORTRAN77 version'

      call test01
      call test02
      call test03

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'SPARSEKIT_PRB08'
      write ( *, '(a)' ) '  Normal end of execution.'

      write ( *, '(a)' ) ' '
      call timestamp ( )

      stop
      end
      subroutine test01

c*********************************************************************72
c
cc TEST01 runs the first test.
c
      implicit none

      integer nmax
      parameter ( nmax = 1000 )

      integer nzmax
      parameter ( nzmax = 20 * nmax )

      double precision a(nzmax)
      double precision alpha
      character * ( 2 ) guesol
      integer ia(nzmax)
      integer ic
      integer ierr
      integer ifmt
      integer indx
      integer ios
      integer iout
      integer iwk(nmax)
      integer ja(nzmax)
      integer job
      character * ( 3 ) key
      integer m
      integer n
      integer nn
      integer nz
      double precision rhs(1)
      character * ( 72 ) title
      character * ( 8 ) type

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST01'
      write ( *, '(a)' ) '  Write Zlatev example to file.'
c
c  Call MATRF2 to set up the matrix in COO (coordinate) format.
c
      m = 100
      n = m
      ic = n / 2
      indx = 10
      alpha = 5.0D+00
      nn = nzmax

      call matrf2 ( m, n, ic, indx, alpha, nn, nz, a, ia, ja, ierr )
c
c  Convert the matrix from COO to CSR format.
c
      job = 1
      call coicsr ( n, nz, job, a, ia, ja, iwk )
c
c  Write the matrix to a file using Harwell-Boeing format.
c
      title = 'First matrix from Zlatev examples'
      type = 'RUA'
      key = ' ZLATEV1'
      guesol = 'NN'
      ifmt = 3
      job = 2
      iout = 7

      open ( unit = iout, file = 'zlatev1_hb.txt', status = 'replace',
     &  err = 99 )

      call prtmt ( n, n, a, ia, ja, rhs, guesol, title, type, key, 
     &  ifmt, job, iout )

      close ( unit = iout )

      return

99    continue

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST01 - Fatal error!'
      write ( *, '(a)' ) '  Could not open the output file.'
      return
      end
      subroutine test02

c*********************************************************************72
c
cc TEST02 runs the second test, with DCN matrices.
c
      implicit none

      integer nmax
      parameter ( nmax = 1000 )
      integer nzmax
      parameter ( nzmax = 20 * nmax )

      double precision a(nzmax)
      character * ( 2 ) guesol
      integer ia(nzmax)
      integer ic
      integer ierr
      integer ifmt
      integer ios
      integer iout
      integer iwk(nmax)
      integer ja(nzmax)
      integer job
      character * ( 3 ) key
      integer n
      integer ne
      integer nn
      double precision rhs(1)
      character * ( 72 ) title
      character * ( 8 ) type 

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST02'
      write ( *, '(a)' ) '  Write DCN example to file.'
c
c  Call DCN to set up the matrix in COO format.
c
      n = 200
      nn = nzmax
      ic = 20
      call dcn ( a, ia, ja, n, ne, ic, nn, ierr )
c
c  Convert from COO to CSR format.
c
      job = 1

      call coicsr ( n, ne, job, a, ia, ja, iwk )
c
c  Write the matrix to file using Harwell Boeing format.
c
      title = 'second matrix from Zlatev examples'
      guesol = 'NN'
      type = 'RUA'
      key = ' ZLATEV2'
      ifmt = 3
      job = 2
      iout = 7

      open ( unit = iout, file = 'zlatev2_hb.txt', status = 'replace', 
     &  err = 99 )

      call prtmt ( n, n, a, ia, ja, rhs, guesol, title, type, key, 
     &  ifmt, job, iout )

      close ( unit = iout )

      return

99    continue

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST02 - Fatal error!'
      write ( *, '(a)' ) '  Could not open the output file.'

      return
      end
      subroutine test03

c*********************************************************************72
c
cc TEST03 runs the test with ECN matrices.
c
      implicit none

      integer nmax
      parameter ( nmax = 1000 )
      integer nzmax
      parameter ( nzmax = 20 * nmax )

      double precision a(nzmax)
      character * ( 2 ) guesol
      integer ia(nzmax)
      integer ic
      integer ierr
      integer ifmt
      integer ios
      integer iout
      integer iwk(nmax)
      integer ja(nzmax)
      integer job
      character * ( 3 ) key
      integer n
      integer ne
      integer nn
      double precision rhs(1)
      character * ( 72 ) title
      character * ( 8 ) type 

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST03'
      write ( *, '(a)' ) '  Write ECN example to file.'
c
c  Call ECN to set up the matrix in COO format.
c
      n = 200
      ic = 20
      nn = nzmax
      call ecn ( n, ic, ne, ia, ja, a, nn, ierr )
c
c  Convert the matrix from COO to CSR format.
c
      job = 1
      call coicsr ( n, ne, job, a, ia, ja, iwk )
c
c  Store matrix to a file using Harwell Boeing format.
c
      title = 'Third matrix from Zlatev examples'
      guesol = 'NN'
      type = 'RUA'
      key = ' ZLATEV3'
      ifmt = 3
      job = 2
      iout = 7

      open ( unit = iout, file = 'zlatev3_hb.txt', status = 'replace', 
     &  err = 99 )

      call prtmt ( n, n, a, ia, ja, rhs, guesol, title, type, key, 
     &  ifmt, job, iout )

      close ( unit = iout )

      return

99    continue

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST03 - Fatal errorc'
      write ( *, '(a)' ) '  Could not open the output file.'
      return
      end
