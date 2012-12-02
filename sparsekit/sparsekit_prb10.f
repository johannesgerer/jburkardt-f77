      program main

c*********************************************************************72
c
cc MAIN is the main program for SPARSEKIT_PRB10
c
c main program for generating BLOCK 5 point and 7-point matrices in the
c Harwell-Boeing format.  Creates a file containing a
c harwell-boeing matrix.
c
      implicit none

      integer nxmax
      parameter ( nxmax = 30 )
      integer nmx
      parameter ( nmx = nxmax * nxmax )

      double precision a(9,5*nmx)
      double precision ao(45*nmx)
      character * ( 2 ) guesol
      integer ia(nmx)
      integer iao(5*nmx)
      integer iau(nmx)
      integer ifmt
      integer iout
      integer ja(7*nmx)
      integer jao(15*nmx)
      integer job
      character * ( 8 ) key
      character * ( 50 ) matfil
      integer n
      integer na
      integer nfree
      integer nfree2
      integer nr
      integer nx
      integer ny
      integer nz
      double precision rhs(1)
      double precision stencil(7,100)
      character * ( 72 ) title
      character * ( 3 ) type

      call timestamp ( )

      WRITE(*,*)' '
      write ( *, '(a)' ) 'SPARSEKIT_PRB10'
      write ( *, '(a)' ) '  FORTRAN77 version'
      write ( *, '(a)' ) ' '
      WRITE(*,*)'This program demonstrates the use of GEN57BL'
      WRITE(*,*)'to generate a block sparse matrix derived from a 5 or'
      WRITE(*,*)'7 point stencil on an NX by NY by NZ grid.'
      WRITE(*,*)' '
      WRITE(*,*)'In this example, the block size is 3.'
      WRITE(*,*)' '
      WRITE(*,*)'The matrix is then stored in Harwell-Boeing format'
      WRITE(*,*)'in a file, using routine PRTMT.'
      WRITE(*,*)' '
c
c  NFREE is the block size.
c
      nx = 10
      ny = 10
      nz = 1
      nfree = 3
      na = 9
      call gen57bl ( nx, ny, nz, nfree, na, n, a, ja, ia, iau, stencil )
c 
c  Convert from BSR (block sparse row) to CSR (column sparse row) format.
c
      nfree2 = nfree * nfree
      nr = n / nfree
      call bsrcsr ( n, nfree, na, a, ja, ia, ao, jao, iao )
c
c  Set other data needed for Harwell Boeing format.
c
      guesol = 'NN'
      title = ' BLOCK 5-POINT TEST MATRIX FROM SPARSKIT               '
      type = 'RUA'
      key = 'BLOCK5PT'
      ifmt = 104
      job = 2
c
c  Store matrix in file, using Harwell Boeing format.
c
      matfil = 'test.mat'
      iout = 7

      open ( unit = iout, file = matfil, status = 'replace' )

      call prtmt ( n, n, ao, jao, iao, rhs, guesol, title, key, type, 
     &  ifmt, job, iout )

      close ( unit = iout )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'SPARSEKIT_PRB10'
      write ( *, '(a)' ) '  Normal end of execution.'

      write ( *, '(a)' ) ' '
      call timestamp ( )

      stop
      end
      subroutine afunbl (nfree,x,y,z,coeff)

c*********************************************************************72
c
      implicit none

      integer i
      integer j
      integer nfree
      double precision x, y, z, coeff(100)

      do j=1, nfree
        do i=1, nfree
          coeff((j-1)*nfree+i) = 0.0
        end do
        coeff((j-1)*nfree+j) = -1.0
      end do

      return
      end
      subroutine bfunbl (nfree,x,y,z,coeff)

c*********************************************************************72
c
      implicit none

      integer i
      integer j
      integer nfree
      double precision x, y, z, coeff(100)

      do j=1, nfree
        do i=1, nfree
          coeff((j-1)*nfree+i) = 0.0
        end do
        coeff((j-1)*nfree+j) = -1.0
      end do

      return
      end
      subroutine cfunbl (nfree,x,y,z,coeff)

c*********************************************************************72
c
      implicit none

      integer i
      integer j
      integer nfree
      double precisionx, y, z, coeff(100)

      do j=1, nfree
        do i=1, nfree
          coeff((j-1)*nfree+i) = 0.0
        end do
        coeff((j-1)*nfree+j) = -1.0
      end do

      return
      end
      subroutine dfunbl (nfree,x,y,z,coeff)

c*********************************************************************72
c
      implicit none

      integer i
      integer j
      integer nfree
      double precisionx, y, z, coeff(100)

      do j=1, nfree
        do i=1, nfree
          coeff((j-1)*nfree+i) = 0.0
        end do
      end do

      return
      end
      subroutine efunbl (nfree,x,y,z,coeff)

c*********************************************************************72
c
      implicit none

      integer i
      integer j
      integer nfree
      double precisionx, y, z, coeff(100)

      do j=1, nfree
        do i=1, nfree
          coeff((j-1)*nfree+i) = 0.0
        end do
      end do

      return
      end
      subroutine ffunbl (nfree,x,y,z,coeff)

c*********************************************************************72
c
      implicit none

      integer i
      integer j
      integer nfree
      double precision x, y, z, coeff(100)

      do j=1, nfree
        do i=1, nfree
          coeff((j-1)*nfree+i) = 0.0
        end do
      end do

      return
      end
      subroutine gfunbl (nfree,x,y,z,coeff)

c*********************************************************************72
c
      implicit none

      integer i
      integer j
      integer nfree
      double precision x, y, z, coeff(100)

      do j=1, nfree
        do i=1, nfree
          coeff((j-1)*nfree+i) = 0.0
        end do
      end do

      return
      end
