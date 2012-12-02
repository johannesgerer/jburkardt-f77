      program main

c*********************************************************************72
c
cc MAIN is a test suite for Part I : I/O routines.
c
c tests the following : gen5pt.f, prtmt, readmt, amd pltmt.
c 1) generates a 100 x 100 5pt matrix,
c 2) prints it with a given format in file 'first.mat'
c 3) reads the matrix from 'first.mat' using readmat
c 4) prints it again in file 'second.mat' in  a different format
c 5) makes 4 pic files to show the different options of pltmt.
c    these are in job0.pic, job01.pic, job10.pic, job11.pic 
c                          coded by Y. Saad, RIACS, 08/31/1989.
c
      implicit none

      integer nxmax
      parameter ( nxmax = 20 )
      integer nmx
      parameter ( nmx = nxmax * nxmax )

      double precision a(7*nmx)
      character * ( 2 ) guesol
      integer i
      integer ia(nmx)
      integer iau(nmx)
      integer ierr
      integer ifmt
      integer iout
      integer j
      integer ja(7*nmx)
      integer job
      integer k
      character * ( 8 ) key
      integer mode
      integer n
      integer ncol
      integer nmax
      integer nnz
      integer nrhs
      integer nrow
      integer nx
      integer ny
      integer nz
      integer nzmax
      double precision rhs(3*nmx)
      double precision stencil(7)
      character * ( 72 ) title
      character * ( 3 ) type

      call timestamp ( )

      WRITE(*,*)' '
      WRITE(*,*)'SPARSEKIT_PRB05'
      write ( *, '(a)' ) '  FORTRAN77 version'
      write ( *, '(a)' ) ' '
      WRITE(*,*)'A set of tests for SPARSKIT'
      WRITE(*,*)' '

      open (unit=7,file='first.mat',STATUS='replace')
      open (unit=8,file='second.mat',STATUS='replace')
      open (unit=20,file='job00.pic',STATUS='replace')
      open (unit=21,file='job01.pic',STATUS='replace')
      open (unit=22,file='job10.pic',STATUS='replace')
      open (unit=23,file='job11.pic',STATUS='replace')
c
c  dimension of grid
c
      nx = 10
      ny = 10
      nz = 1
c
c  generate grid problem.
c
      call gen57pt (nx,ny,nz,a,ja,ia,iau,stencil)
c
c  create the Harwell-Boeing matrix. Start by defining title,
c  and type. them define format and print it.
c
      write (title,9) nx, ny
 9    format('Five-point matrix on a square region', 
     &  ' using a ',I2,' by ',I2,' grid *SPARSKIT*')
      key = 'Fivept10'
      type= 'RSA'
      n = nx*ny*nz
      ifmt = 5
      job = 3
      guesol = 'GX'
c
c  define a right hand side of ones, an initial guess of two's
c  and an exact solution of three's.
c
      do k=1, 3*n
        rhs(k) = real( 1 +  (k-1)/n )
      end do

      call prtmt (n,n,a,ja,ia,rhs,guesol,title,key,type,ifmt,job,7)
c
c  read it again in same matrix a, ja, ia
c
      nmax = nmx
      nzmax = 7*nmx
      do k=1, 3*n
        rhs(k) = 0.0
      end do
      job = 3

      rewind 7
      nrhs = 3*n

      call readmt (nmax,nzmax,job,7,a,ja,ia,rhs,nrhs,guesol, 
     &             nrow,ncol,nnz,title,key,type,ierr)
      WRITE(*,*)'ierr = ',ierr,' nrhs ', nrhs
c
c  matrix read.  print it again in a different format
c
      ifmt = 102
      ncol = nrow
      job = 3

      call prtmt (nrow,ncol,a,ja,ia,rhs,guesol,title,key,type, 
     &                        ifmt,job,8)
c
c  print four pic files
c
      mode = 0
      do i=1, 2
        do j=1, 2
          job = (i-1)*10 +j-1
          iout = 20+(i-1)*2+j-1
          call pltmt(nrow,ncol,mode,ja,ia,title,key,type,job,iout)
        end do
      end do

      CLOSE(UNIT=7)
      CLOSE(UNIT=8)
      CLOSE(UNIT=20)
      CLOSE(UNIT=21)
      CLOSE(UNIT=22)
      CLOSE(UNIT=23)

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'SPARSEKIT_PRB05'
      write ( *, '(a)' ) '  Normal end of execution.'

      write ( *, '(a)' ) ' '
      call timestamp ( )

      stop
      end
      function afun (x,y,z)

c*********************************************************************72
c
cc AFUN ???
c
      implicit none

      double precision afun, x,y, z

      afun = -1.0

      return
      end
      function bfun (x,y,z)

c*********************************************************************72
c
cc BFUN ???
c
      implicit none

      double precision bfun, x,y, z

      bfun = -1.0

      return
      end
      function cfun (x,y,z)

c*********************************************************************72
c
cc CFUN ???
c
      implicit none

      double precision cfun, x,y, z

      cfun = -1.0

      return
      end
      function dfun (x,y,z)

c*********************************************************************72
c
cc DFUN ???
c
      implicit none

      double precision dfun, gamma, x,y, z

      data gamma /100.0/

      dfun = 10.0

      return
      end
      function efun (x,y,z)

c*********************************************************************72
c
cc EFUN ???
c
      implicit none

      double precision efun, gamma, x,y, z

      data gamma /100.0/

      efun = 0.0

      return
      end
      function ffun (x,y,z)

c*********************************************************************72
c
cc FFUN ???
c
      implicit none

      double precision ffun, x,y, z

      ffun = 0.0

      return
      end
      function gfun (x,y,z)

c*********************************************************************72
c
cc GFUN ???
c
      implicit none

      double precision gfun, x,y, z

      gfun = 0.0

      return
      end
