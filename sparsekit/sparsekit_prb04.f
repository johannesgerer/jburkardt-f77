      program main

c*********************************************************************72
c
cc MAIN is the main program for SPARSEKIT_PRB04.
c
c  test suite for the unary routines.
c  tests some of the routines in the module unary. Still needs to tests
c  many other routines.
c
      integer nxmax
      parameter ( nxmax = 10 )
      integer nmx
      parameter ( nmx = nxmax * nxmax )
      integer nnzmax
      parameter ( nnzmax = 10 * nmx )
      integer ndns
      parameter ( ndns = 20 )

      double precision a(nnzmax)
      double precision a1(nnzmax)
      double precision dns(ndns,ndns)
      integer i
      integer ia(nmx+1)
      integer ia1(nnzmax)
      integer iwk(nmx*2+1)
      integer iwork(nnzmax*2)
      integer ja(nnzmax)
      integer ja1(nnzmax)
      integer perm(16)
      integer qperm(16)
      double precision stencil(100)

      data qperm / 1, 3, 6, 8, 9, 11, 14, 16, 2, 4, 5, 7, 
     &  10, 12, 13, 15 /

      call timestamp ( )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'SPARSEKIT_PRB04'
      write ( *, '(a)' ) '  FORTRAN77 version'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  A set of tests for SPARSEKIT.'
c
c  define correct permutation
c
      do i=1, 16
        perm(qperm(i)) = i
      end do
c
c  dimension of grid
c
          nx = 4
          ny = 4
          nz = 1
c
c  generate grid problem.
c
          call gen57pt(nx,ny,nz,a,ja,ia,iwk,stencil)
          n = nx*ny*nz
c
c  write out the matrix
c
          nnz = ia(n+1)-1
          WRITE(*,*) '-----------------------------------------'
          WRITE(*,*) '  +++  initial matrix in CSR format +++ '
          WRITE(*,*) '-----------------------------------------'
          call dump ( n, a, ja, ia, 6 )
c
c  call csrdns
c
          call csrdns(n,n,a,ja,ia,dns,ndns,ierr)
c
c  write it out as a dense matrix.
c
          WRITE(*,*) '-----------------------------------------'
          WRITE(*,*) '  +++ initial matrix in DENSE format+++ '
          WRITE(*,*) '-----------------------------------------'
          call dmpdns ( n, n, ndns, dns )
c
c  red black ordering
c
          job = 1
          call dperm (n,a,ja,ia,a1,ja1,ia1,perm,perm,job)
          nnz = ia(n+1)-1
          WRITE(*,*) '-----------------------------------------'
          WRITE(*,*) '  +++ red-black matrix in CSR format +++ '
          WRITE(*,*) '-----------------------------------------'
          call dump (n,a1,ja1,ia1,6)
c
c  sort matrix
c
          call csort (n,a1,ja1,ia1,iwork,.true.)
          nnz = ia(n+1)-1
          WRITE(*,*) '-----------------------------------------'
          WRITE(*,*) '  +++     matrix after sorting    +++ '
          WRITE(*,*) '-----------------------------------------'
          call dump (n,a1,ja1,ia1,6)
c
c  convert into dense format
c
           call csrdns(n, n, a1,ja1,ia1,dns,ndns,ierr)
           WRITE(*,*) '-----------------------------------------'
           WRITE(*,*) '  +++ red-black matrix in DENSE format+++ '
           WRITE(*,*) '-----------------------------------------'
           call dmpdns(n,n, ndns, dns)

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'SPARSEKIT_PRB04'
      write ( *, '(a)' ) '  Normal end of execution.'

      write ( *, '(a)' ) ' '
      call timestamp ( )

      stop
      end
      subroutine dmpdns ( nrow, ncol, ndns, dns )

c*********************************************************************72
c
c this subroutine prints out a dense matrix in a simple format.
c the zero elements of the matrix are omitted. The format for the
c nonzero elements is f4.1, i.e., very little precision is provided.
c
c on entry
c
c nrow = row dimension of matrix
c ncol = column dimension of matrix
c ndns = first dimension of array dns.
c dns  = double dimensional array of size n x n containing the matrix
c
c on return
c ---------
c matrix will be printed out.
c
      integer ndns

      double precision dns(ndns,*)
      character ( len = 80 ) fmt
      integer i
      integer j
      integer j1
      integer j2
      integer last
      integer ncol
      integer nrow
c
c prints out a dense matrix -- without the zeros.
c
      write ( *, '(4x,16i4)' ) ( j, j = 1, ncol )

          fmt(1:5) = '    |'
          j1 = 6
          do j=1, ncol
            j2 = j1+4
            fmt(j1:j2) = '----'
            j1 = j2
          end do
          last = j1
          fmt(last:last) = '|'
          write ( *, '(a)' ) fmt(1:last)
c
c undo loop 1 ---
c
          j1 = 6
          do j=1,ncol
            j2 = j1+4
            fmt(j1:j2) = '   '
            j1 = j2
          end do

      do i=1, nrow
            j1 = 6
            write (fmt,101) i
101         format(1x,i2,' |')
            do j=1, ncol
              j2= j1+4
              if ( dns(i,j) /= 0.0 ) then
                write ( fmt(j1:j2), '(f4.1)' ) dns(i,j)
102             format(f4.1)
                endif
              j1 = j2
              end do
            fmt(last:last) = '|'
            write ( *, '(a)' ) fmt(1:last)
      end do

          fmt(1:5) = '    |'
          j1 = 6
          do j=1, ncol
            j2 = j1+4
            fmt(j1:j2) = '----'
            j1 = j2
          end do
          fmt(last:last) = '|'
          write ( *, '(a)' ) fmt(1:last)

      return
      end
      function afun ( x, y, z )

c*********************************************************************72
c
      double precision afun, x,y, z

      afun = -1.0D+00

      return
      end
      function bfun (x,y,z)

c*********************************************************************72
c
      double precision bfun, x,y, z
      bfun = -1.0D+00
      return
      end
      function cfun (x,y,z)

c*********************************************************************72
c
      double precision cfun, x,y, z
      cfun = -1.0D+00
      return
      end
      function dfun (x,y,z)

c*********************************************************************72
c
      double precision dfun, x,y, z
      data gamma /100.0D+00/
      dfun = 10.0D+00
      return
      end
      function efun (x,y,z)

c*********************************************************************72
c
      double precision efun, x,y, z
      data gamma /100.0D+00/
      efun = 0.0D+00
      return
      end
      function ffun (x,y,z)

c*********************************************************************72
c
      double precision ffun
      double precision x
      double precision y
      double precision z

      ffun = 0.0D+00

      return
      end
      function gfun (x,y,z)

c*********************************************************************72
c
      double precision gfun, x,y, z
      gfun = 0.0D+00
      return
      end
