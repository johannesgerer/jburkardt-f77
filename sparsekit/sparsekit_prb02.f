      program main

c*********************************************************************72
c
cc MAIN tests all the subroutines in matvec.
c
c  Discussion:
c
c    This test generates matrices and transforms them in appropriate formats
c    and then calls the appropriate routines.
c
c  NMAX is the maximum number of rows in the matrix, plus 1.
c
      implicit none

      integer nmax
      parameter ( nmax=10*10*10+1 )
c
c  NZMAX is the maximum number of nonzero entries in the matrix.
c  There are at most 7 nonzeroes in the 3-D Laplacian operator.
c
      integer nzmax
      parameter ( nzmax = 7 * nmax )

      double precision a1(nzmax)
      double precision a2(nzmax)
      integer ia1(nmax)
      integer ia2(nmax)
      integer idiag
      integer idim(2)
      integer ierr
      integer ii
      integer ioff(10)
      integer iwk1(nmax)
      integer iwk2(nmax)
      integer j
      integer ja1(nzmax)
      integer ja2(nzmax)
      integer jad(nzmax)
      integer jdiag
      integer jj
      integer k
      integer n
      integer na
      integer ndiag
      integer nfree
      integer nlev
      integer nx
      integer ny
      integer nz
      double precision scal
      double precision stencil(100)
      double precision x(nmax)
      double precision y(nmax)
      double precision y0(nmax)
      double precision y1(nmax)

      idim(1) = 4
      idim(2) = 10

      call timestamp ( )

      WRITE(*,*)' '
      WRITE(*,*)'SPARSEKIT_PRB02'
      write ( *, '(a)' ) '  FORTRAN77 version'
      WRITE(*,*)' '
      WRITE(*,*)'A set of tests for SPARSKIT.'
      WRITE(*,*)'These test the conversion routines, and the'
      WRITE(*,*)'lower and upper triangular solution routines.'
      WRITE(*,*)' '
c
c  ii loop controls the size of the grid.
c
      do ii = 1, 2

        WRITE(*,*) '---------------- ii ',ii,'--------------------'
        nfree = 1
        NA=NFREE*NFREE
        nx = idim(ii)
        ny = nx
c
c  jj loop corresponds to 2-D and 3-D problems.
c
        do jj=1, 2

           WRITE(*,*) '     ----------- jj ',jj,' -------------'
           nz = 1
           if (jj .eq. 2) nz = 10
c
c  call matrix generation routine to create matrix in CSR format.
c
            WRITE(*,*)'nx=',nx
            WRITE(*,*)'ny=',ny
            WRITE(*,*)'nz=',nz
            call gen57bl(nx,ny,nz,nfree,na,n,a1,ja1,ia1,ia2,stencil)
            WRITE(*,*)'The solution vector will have length N=',N
c
c  initialize the solution vector x
c
            do j=1, n
               x(j) = real ( j, kind = 8 )
            end do
c
c  Get exact answer in y0 for later checks.
c
            call amux(n,x,y0,a1,ja1,ia1)
c
c  Convert CSR format to ELL format, and test AMUXE.
c
            call csrell(n,a1,ja1,ia1,7,a2,jad,n,ndiag,ierr)
            call amuxe(n, x, y, n, ndiag, a2,jad)
            call errpr(n, y, y0,'amuxe ')
c
c  Convert CSR format to DIA format, and test AMUXD.
c
            idiag = 7
            call csrdia (n, idiag,10,a1, ja1, ia1, nmax, a2, 
     &           ioff, a2, ja2, ia2, jad)
            call amuxd (n,x,y,a2,nmax,idiag,ioff)
            call errpr (n, y, y0,'amuxd ')
c
c  Convert CSR format to CSC format, and test ATMUX.
c
            call csrcsc(n,1,1,a1,ja1,ia1,a2,ja2,ia2)
            call atmux(n,x,y,a2,ja2,ia2)
            call errpr(n, y, y0,'atmux ')
c
c  Convert CSR format to JAD format, and test AMUXJ.
c
            call csrjad (n,a1,ja1,ia1, jdiag, jad, a2, ja2, ia2)
            call amuxj (n, x, y, jdiag, a2, ja2, ia2)
            call dvperm (n, y, jad)
            call errpr (n, y, y0,'amuxj ')
c
c  convert JAD format back to CSR, and test JADCSR and AMUX.
c
            call jadcsr (n, jdiag, a2, ja2, ia2, jad, a1, ja1, ia1)
            call amux (n, x, y, a1, ja1, ia1)
            call errpr (n, y, y0,'jadcsr')
c
c  triangular systems solutions
c
c TESTING LDSOL
c
            call getl (n, a1, ja1, ia1, a2, ja2, ia2)
            call amux (n,x,y0, a2, ja2, ia2)
            call atmux(n,x,y1, a2, ja2, ia2)
            call csrmsr (n, a2, ja2, ia2, a2, ja2, y, iwk2)
            do k=1,n
               a2(k) = 1.0/ a2(k)
            end do
            call ldsol (n, y, y0, a2, ja2)
            call errpr (n, x, y, 'ldsol ')
c
c TESTING LDSOLL
c
            call levels (n, ja2, ja2, nlev, jad, iwk1, iwk2)
            call ldsoll (n, y, y0, a2, ja2, nlev, jad, iwk1)
            call errpr (n, x, y,'ldsoll')
c
c TESTING UDSOLC
c
c here we take advantage of the fact that the MSR format for U
c is the MSC format for L
c
            call udsolc (n, y, y1, a2, ja2)
            call errpr (n, x, y,'udsolc')
c
c TESTING LSOL
c
c here we exploit the fact that with MSR format a, ja, ja is actually
c the correct data structure for the strict lower triangular part of
c the CSR format. First rescale matrix.
c
            scal = 0.1
            do k=ja2(1), ja2(n+1)-1
               a2(k)=a2(k)*scal
            end do
            call amux(n, x, y0, a2, ja2, ja2)
            do j=1,n
               y0(j) = x(j) + y0(j)
            end do
            call lsol (n, y, y0, a2, ja2, ja2)
            call errpr (n, x, y,'lsol  ')
c
c TESTING UDSOL
c
            call getu (n, a1, ja1, ia1, a2, ja2, ia2)
            call amux (n,x,y0, a2, ja2, ia2)
            call atmux(n,x,y1, a2, ja2, ia2)
            call csrmsr (n, a2, ja2, ia2, a2, ja2, y, jad)
            do k=1,n
               a2(k) = 1.0/ a2(k)
            end do
            call udsol (n, y, y0, a2, ja2)
            call errpr (n, x, y,'udsol ')
c
c TESTING LDSOLC
c
c here we take advantage of the fact that the MSR format for L
c is the MSC format for U
c
            call ldsolc (n, y, y1, a2, ja2)
            call errpr (n, x, y,'ldsolc')
c
c TESTING USOL
c
c here we exploit the fact that with MSR format a, ja, ja is actually
c the correct data structure for the strict lower triangular part of
c the CSR format. First rescale matrix.
c
            scal = 0.1
            do k=ja2(1), ja2(n+1)-1
               a2(k)=a2(k)*scal
            end do
            call amux(n, x, y1, a2, ja2, ja2)
            do j=1,n
               y1(j) = x(j) + y1(j)
            end do
            call usol (n, y, y1, a2, ja2, ja2)
            call errpr (n, x, y,'usol  ')
        end do
      end do

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'SPARSEKIT_PRB02'
      write ( *, '(a)' ) '  Normal end of execution.'

      write ( *, '(a)' ) ' '
      call timestamp ( )

      stop
      end
      subroutine afunbl ( nfree, x, y, z, coeff )

c*********************************************************************72
c
cc AFUNBL ???
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
      double precision x, y, z, coeff(100)

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
      double precision x, y, z, coeff(100)

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
      double precision x, y, z, coeff(100)

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
      function afun ()

c*********************************************************************72
c
      implicit none

      double precision afun

      afun = 0.0

      return
      end
      function bfun ()

c*********************************************************************72
c  
      implicit none

      double precision bfun

      bfun = 0.0

      return
      end
      function cfun ()

c*********************************************************************72
c
      implicit none

      double precision cfun

      cfun = 0.0

      return
      end
      function dfun ( )

c*********************************************************************72
c
      implicit none

      double precision dfun

      dfun = 0.0

      return
      end
      function efun ( )

c*********************************************************************72
c
      implicit none

      double precision efun

      efun = 0.0

      return
      end
      function ffun ( )

c*********************************************************************72
c
      implicit none

      double precision ffun

      ffun = 0.0

      return
      end
      function gfun ( )

c*********************************************************************72
c 
      implicit none

      double precision gfun

      gfun = 0.0

      return
      end
      subroutine errpr ( n, y, y1, msg )

c*********************************************************************72
c
      implicit none

      integer n

      integer k
      double precision y(n), y1(n), t
      character*6 msg

      t = 0.0
      do k=1,n
        t = t+(y(k)-y1(k))**2
      end do

      t = sqrt(t)
      WRITE(*,*)'RMS error in ',msg,' =', t

      return
      end
