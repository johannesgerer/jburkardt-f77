      program main

c*********************************************************************72
c
cc SPARSEKIT_PRB01
c
c  Modified:
c
c    17 August 2008
c
      implicit none

      integer nxmax
      parameter ( nxmax = 30 )
      integer nmx
      parameter ( nmx = nxmax * nxmax )
      integer nzmax
      parameter ( nzmax = 7 * nmx )

      double precision a(nzmax)
      double precision alpha
      double precision b(nzmax)
      double precision beta
      double precision c(nzmax)
      integer i
      integer ia(nmx+1)
      integer ib(nmx+1)
      integer ic(nzmax)
      integer ierr
      integer ifmt
      integer iw(nmx)
      integer j
      integer ja(nzmax)
      integer jb(nzmax)
      integer jc(nzmax)
      integer jj
      integer job
      integer k
      character*8 key
      integer n
      integer na
      integer nx
      integer ny
      integer nz
      double precision s
      double precision stencil(100)
      character*7 title
      character*3 type
      double precision x(nmx)
      double precision y(nmx)
      double precision y1(nmx)

      call timestamp ( )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'SPARSEKIT_PRB01'
      write ( *, '(a)' ) '  FORTRAN77 version'
      write ( *, '(a)' ) '  A test problem for SPARSEKIT.'

      nx=10
      ny=10
      nz = 1
      na = 5*nmx
      call gen57pt ( nx, ny, nz, a, ja, ia, iw, stencil )
      beta  = alpha
      alpha = 0.0D+00

      call gen57pt(ny,nx,nz,b,jb,ib,iw,stencil)
      n = nx*ny*nz

      s = 3.812D+00

      call aplsb1(n,n,a,ja,ia,s,b,jb,ib,c,jc,ic,nzmax,ierr)
      if (ierr .ne. 0)WRITE(*,*)' ierr = ',ierr

      call dump(n,c,jc,ic,6)

      do k=1,n
        x(k) = dble ( k ) / dble ( n )
      end do

      call ope (n,x,y1,a,ja,ia)
      call ope (n,x,y,b,jb,ib)
      do j=1, n
        y1(j) = s*y(j) + y1(j)
      end do

      call ope (n,x,y,c,jc,ic)

      write (*,*) ' ------------ checking APLSB --------------'
      call ydfnorm(n,y1,y)

      type = '--------'
      title=' test matrix for blassm c = a+b '
      key = 'rua'

      ifmt = 103

      job = -1

      do jj=1,2
        call apmbt(n,n,job,a,ja,ia,b,jb,ib,c,jc,ic,nzmax,iw,ierr)
        if (ierr .ne. 0) WRITE(*,*)' ierr = ',ierr
        call ope (n,x,y1,a,ja,ia)
        call opet (n,x,y,b,jb,ib)
        s = real ( job, kind = 8 )

        do i = 1, n
          y1(i) = y1(i) + s * y(i)
        end do

        call ope (n,x,y,c,jc,ic)

        WRITE(*,*) '  '
        WRITE(*,*) ' ------------ checking APMBT---------------'
        WRITE(*,*) ' ------------ with JOB = ',job,' -------------'
        call ydfnorm(n,y1,y)

        job = job + 2
      end do

      type = '--------'
      title=' test matrix for blassm c = a+b^T '

      s = 0.1232445D+00
      call aplsbt(n,n,a,ja,ia,s,b,jb,ib,c,jc,ic,nzmax,iw,ierr)

      if (ierr .ne. 0) WRITE(*,*)' ierr = ',ierr
      call ope (n,x,y1,a,ja,ia)
      call opet (n,x,y,b,jb,ib)

      do i = 1, n
        y1(i) = y1(i) + s * y(i)
      end do

      call ope (n,x,y,c,jc,ic)

      WRITE(*,*) '  '
      WRITE(*,*) ' ------------ checking APLSBT---------------'
      call ydfnorm(n,y1,y)
c
c  Testing matrix products
c
      job = 1
      call amub (n,n,job,a,ja,ia,b,jb,ib,c,jc,ic,nzmax,iw,ierr)

      if (ierr .ne. 0) WRITE(*,*)' ierr = ',ierr
      call ope(n,x,y,b,jb,ib)
      call ope(n,y,y1,a,ja,ia)

      call ope(n,x,y,c,jc,ic)

      WRITE(*,*) '  '
      WRITE(*,*) ' ------------ checking AMUB  ---------------'
      call ydfnorm(n,y1,y)

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'SPARSEKIT_PRB01'
      write ( *, '(a)' ) '  Normal end of execution.'

      write ( *, '(a)' ) ' '
      call timestamp ( )

      stop
      end
      function afun ( x, y, z )

c*********************************************************************72
c
cc AFUN
c
      implicit none

      double precision afun
      double precision x
      double precision y
      double precision z

      afun = -1.0D+00

      return
      end
      function bfun ( x, y, z )

c*********************************************************************72
c
cc BFUN
c
      implicit none

      double precision bfun
      double precision x
      double precision y
      double precision z

      bfun = -1.0D+00

      return
      end
      function cfun ( x, y, z )

c*********************************************************************72
c
cc CFUN
c
      implicit none

      double precision cfun
      double precision x
      double precision y
      double precision z

      cfun = -1.0D+00

      return
      end
      function dfun ( x, y, z )

c*********************************************************************72
c
cc DFUN
c
      implicit none

      double precision dfun
      double precision x
      double precision y
      double precision z

      dfun = 0.0D+00

      return
      end
      function efun ( x, y, z )

c*********************************************************************72
c
cc EFUN
c
      implicit none

      double precision efun
      double precision x
      double precision y
      double precision z

      efun = 0.0D+00

      return
      end
      function ffun ( x, y, z )

c*********************************************************************72
c
cc FFUN
c
      implicit none

      double precision ffun
      double precision x
      double precision y
      double precision z

      ffun = 0.0D+00

      return
      end
      function gfun ( x, y, z )

c*********************************************************************72
c
cc GFUN
c
      implicit none

      double precision gfun
      double precision x
      double precision y
      double precision z

      gfun = 0.0D+00

      return
      end
      subroutine ope ( n, x, y, a, ja, ia )
     
c*********************************************************************72
c
cc OPE computes A * x for a sparse matrix A.
c
      implicit none

      integer n

      double precision a(*)
      integer i
      integer ia(n+1)
      integer ja(*)
      integer k
      integer k1
      integer k2
      double precision x(*)
      double precision y(*)
c
c sparse matrix * vector multiplication
c
      do i=1,n
        k1 = ia(i)
        k2 = ia(i+1) -1
        y(i) = 0.0D+00
        do k = k1, k2
          y(i) = y(i) + a(k) * x(ja(k))
        end do
      end do

      return
      end
      subroutine opet ( n, x, y, a, ja, ia )

c*********************************************************************72
c
cc OPET computes A' * x for a sparse matrix A.
c
      implicit none

      double precision a(*)
      integer i
      integer ia(*)
      integer ja(*)
      integer k
      integer n
      double precision x(*)
      double precision y(*)
c
c sparse matrix * vector multiplication
c
      y(1:n) = 0.0D+00

      do  i=1,n
        do k=ia(i), ia(i+1)-1
          y(ja(k)) = y(ja(k)) + x(i)*a(k)
        end do
      end do

      return
      end
      subroutine ydfnorm ( n, y1, y )

c*********************************************************************72
c
cc YDFNORM prints the L2 norm of the difference of two vectors.
c
      implicit none

      integer n

      integer i
      double precision t
      double precision y(n)
      double precision y1(n)

      t = 0.0D+00
      do i = 1, n
        t = t + ( y(i) - y1(i) )**2
      end do
      t = sqrt ( t )

      write(*,*) '2-norm of error (exact answer-tested answer)=',t

      return
      end
      subroutine dump0 ( n, a, ja, ia )

c*********************************************************************72
c
cc DUMP0
c
      implicit none

      double precision a(*)
      integer i
      integer ia(*)
      integer ja(*)
      integer n
      integer k
      integer k1
      integer k2

      do i = 1, n
        write(*,100) i
        k1=ia(i)
        k2 = ia(i+1)-1
        write(*,101) (ja(k),k=k1,k2)
        write(*,102) (a(k),k=k1,k2)
      end do

 100  format ('row :',i2,20(2h -))
 101  format('     column indices:',10i5)
 102  format('             values:',10f5.1)

      return
      end
      subroutine afunbl ( nfree, x, y, z, coeff )

c*********************************************************************72
c
cc AFUNBL
c
      implicit none

      double precision coeff(100)
      integer i
      integer j
      integer nfree
      double precision x
      double precision y
      double precision z

      do j=1, nfree
        do i=1, nfree
          coeff((j-1)*nfree+i) = 0.0D+00
        end do
        coeff((j-1)*nfree+j) = -1.0D+00
      end do

      return
      end
      subroutine bfunbl ( nfree, x, y, z, coeff )

c*********************************************************************72
c
cc BFUNBL
c
      implicit none

      double precision coeff(100)
      integer i
      integer j
      integer nfree
      double precision x
      double precision y
      double precision z

      do j=1, nfree
        do i=1, nfree
          coeff((j-1)*nfree+i) = 0.0D+00
        end do
        coeff((j-1)*nfree+j) = -1.0D+00
      end do

      return
      end
      subroutine cfunbl ( nfree, x, y, z, coeff )

c*********************************************************************72
c
cc CFUNBL
c
      implicit none

      double precision coeff(100)
      integer i
      integer j
      integer nfree
      double precision x
      double precision y
      double precision z

      do j=1, nfree
        do i=1, nfree
          coeff((j-1)*nfree+i) = 0.0D+00
        end do
        coeff((j-1)*nfree+j) = -1.0D+00
      end do

      return
      end
      subroutine dfunbl ( nfree, x, y, z, coeff )

c*********************************************************************72
c
cc DFUNBL 
c
      implicit none

      double precision coeff(100)
      integer i
      integer j
      integer nfree
      double precision x
      double precision y
      double precision z

      do j=1, nfree
        do i=1, nfree
          coeff((j-1)*nfree+i) = 0.0D+00
        end do
      end do

      return
      end
      subroutine efunbl ( nfree, x, y, z, coeff )

c*********************************************************************72
c
cc EFUNBL
c
      implicit none

      double precision coeff(100)
      integer i
      integer j
      integer nfree
      double precision x
      double precision y
      double precision z

      do j=1, nfree
        do i=1, nfree
          coeff((j-1)*nfree+i) = 0.0D+00
        end do
      end do

      return
      end
      subroutine ffunbl ( nfree, x, y, z, coeff )

c*********************************************************************72
c
cc FFUNBL
c
      implicit none

      double precision coeff(100)
      integer i
      integer j
      integer nfree
      double precision x
      double precision y
      double precision z

      do j=1, nfree
        do i=1, nfree
          coeff((j-1)*nfree+i) = 0.0D+00
        end do
      end do

      return
      end
      subroutine gfunbl ( nfree, x, y, z, coeff )

c*********************************************************************72
c
cc GFUNBL
c
      implicit none

      double precision coeff(100)
      integer i
      integer j
      integer nfree
      double precision x
      double precision y
      double precision z

      do j=1, nfree
        do i=1, nfree
          coeff((j-1)*nfree+i) = 0.0D+00
        end do
      end do

      return
      end
      subroutine xyk ( nel, xyke, x, y, ijk, node )

c*********************************************************************72
c
cc  XYK evaluates the material property function xyk
c
c  Discussion:
c
c    In this version of the routine, the matrix returned is the identity matrix.
c
      implicit none

      integer node

      integer ijk(node,*)
      integer nel
      double precision x(*)
      double precision xyke(2,2)
      double precision y(*)

      xyke(1,1) = 1.0D+00
      xyke(2,2) = 1.0D+00
      xyke(1,2) = 0.0D+00
      xyke(2,1) = 0.0D+00

      return
      end
