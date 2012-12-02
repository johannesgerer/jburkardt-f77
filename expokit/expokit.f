      function clock ( )

c*********************************************************************72
c
cc CLOCK returns a reading from the real time clock.
c
      implicit none

      double precision clock
      real*4 etime, tm(2)
      clock = etime( tm )

*---  IBM & CRAY --
*      double precision function clock()
*      real*8 timef
*      clock = 1000.0d0 * timef()
*      end

*---  others ??


*---  if in trouble, use this to get out of the hook!
*      double precision function clock()
*      clock = 0.0d0
*
      end
      function dasum ( n, dx, incx )

c*********************************************************************72
c
cc DASUM sums the absolute value of the entries of a vector.
c
c     takes the sum of the absolute values.
c     jack dongarra, linpack, 3/11/78.
c     modified 3/93 to return if incx .le. 0.
c
      implicit none

      double precision dasum
      double precision dx(1),dtemp
      integer i,incx,m,mp1,n,nincx
c
      dasum = 0.0d0
      dtemp = 0.0d0
      if( n.le.0 .or. incx.le.0 )return
      if(incx.eq.1)go to 20
c
c        code for increment not equal to 1
c
      nincx = n*incx
      do 10 i = 1,nincx,incx
        dtemp = dtemp + dabs(dx(i))
   10 continue
      dasum = dtemp
      return
c
c        code for increment equal to 1
c
c
c        clean-up loop
c
   20 m = mod(n,6)
      if( m .eq. 0 ) go to 40
      do 30 i = 1,m
        dtemp = dtemp + dabs(dx(i))
   30 continue
      if( n .lt. 6 ) go to 60
   40 mp1 = m + 1
      do 50 i = mp1,n,6
        dtemp = dtemp + dabs(dx(i)) + dabs(dx(i + 1)) + dabs(dx(i + 2))
     *  + dabs(dx(i + 3)) + dabs(dx(i + 4)) + dabs(dx(i + 5))
   50 continue
   60 dasum = dtemp
      return
      end
      subroutine daxpy ( n, da, dx, incx, dy, incy )

c*********************************************************************72
c
cc DAXPY adds a multiple of one vector to another.
c
c     constant times a vector plus a vector.
c     uses unrolled loops for increments equal to one.
c     jack dongarra, linpack, 3/11/78.
c
      implicit none

      double precision dx(1),dy(1),da
      integer i,incx,incy,ix,iy,m,mp1,n
c
      if(n.le.0)return
      if (da .eq. 0.0d0) return
      if(incx.eq.1.and.incy.eq.1)go to 20
c
c        code for unequal increments or equal increments
c          not equal to 1
c
      ix = 1
      iy = 1
      if(incx.lt.0)ix = (-n+1)*incx + 1
      if(incy.lt.0)iy = (-n+1)*incy + 1
      do 10 i = 1,n
        dy(iy) = dy(iy) + da*dx(ix)
        ix = ix + incx
        iy = iy + incy
   10 continue
      return
c
c        code for both increments equal to 1
c
c
c        clean-up loop
c
   20 m = mod(n,4)
      if( m .eq. 0 ) go to 40
      do 30 i = 1,m
        dy(i) = dy(i) + da*dx(i)
   30 continue
      if( n .lt. 4 ) return
   40 mp1 = m + 1
      do 50 i = mp1,n,4
        dy(i) = dy(i) + da*dx(i)
        dy(i + 1) = dy(i + 1) + da*dx(i + 1)
        dy(i + 2) = dy(i + 2) + da*dx(i + 2)
        dy(i + 3) = dy(i + 3) + da*dx(i + 3)
   50 continue
      return
      end
      function dcabs1 ( z )

c*********************************************************************72
c
cc DCABS1 computes the L1 norm of a complex number.
c
      implicit none

      double precision dcabs1
      complex*16 z,zz
      double precision t(2)
      equivalence (zz,t(1))
      zz = z
      dcabs1 = dabs(t(1)) + dabs(t(2))
      return
      end
      subroutine dcmpac ( n, nx, ix, ixx, xx, iwsp )

c*********************************************************************72
c
cc DCMPAC compacts the array IX and sorts IXX and XX.
c
c  This is a gateway routine for DGCNVR) ...
c
      implicit none

      integer          n, nx, ix(nx), ixx(nx), iwsp(n)
      double precision xx(nx)
      integer          k
*
*---  sort ix and carry ixx and xx along ...
*
      call idsrt2( nx, ix, ixx, xx )
*
*---  adjust pointers ...
*
      do k = 1,n
         iwsp(k) = 0
      end do
      do k = 1,nx
         iwsp(ix(k)) = iwsp(ix(k)) + 1
      end do
      ix(n+1) = nx + 1
      do k = n,1,-1
         ix(k) = ix(k+1)-iwsp(k)
      end do
*
*---  sort ixx in increasing order and carry xx along ...
*
      do k = 1,n
         call idsrt1( iwsp(k), ixx(ix(k)), xx(ix(k)) )
      end do
      END
      subroutine dcopy ( n, dx, incx, dy, incy )

c*********************************************************************72
c
cc DCOPY copies a vector.
c
c     copies a vector, x, to a vector, y.
c     uses unrolled loops for increments equal to one.
c     jack dongarra, linpack, 3/11/78.
c
      implicit none

      double precision dx(1),dy(1)
      integer i,incx,incy,ix,iy,m,mp1,n
c
      if(n.le.0)return
      if(incx.eq.1.and.incy.eq.1)go to 20
c
c        code for unequal increments or equal increments
c          not equal to 1
c
      ix = 1
      iy = 1
      if(incx.lt.0)ix = (-n+1)*incx + 1
      if(incy.lt.0)iy = (-n+1)*incy + 1
      do 10 i = 1,n
        dy(iy) = dx(ix)
        ix = ix + incx
        iy = iy + incy
   10 continue
      return
c
c        code for both increments equal to 1
c
c
c        clean-up loop
c
   20 m = mod(n,7)
      if( m .eq. 0 ) go to 40
      do 30 i = 1,m
        dy(i) = dx(i)
   30 continue
      if( n .lt. 7 ) return
   40 mp1 = m + 1
      do 50 i = mp1,n,7
        dy(i) = dx(i)
        dy(i + 1) = dx(i + 1)
        dy(i + 2) = dx(i + 2)
        dy(i + 3) = dx(i + 3)
        dy(i + 4) = dx(i + 4)
        dy(i + 5) = dx(i + 5)
        dy(i + 6) = dx(i + 6)
   50 continue
      return
      end
      function ddot(n,dx,incx,dy,incy)

c*********************************************************************72
c
cc DDOT computes the dot product of two vectors.
c
c     forms the dot product of two vectors.
c     uses unrolled loops for increments equal to one.
c     jack dongarra, linpack, 3/11/78.
c
      implicit none

      double precision ddot
      double precision dx(1),dy(1),dtemp
      integer i,incx,incy,ix,iy,m,mp1,n
c
      ddot = 0.0d0
      dtemp = 0.0d0
      if(n.le.0)return
      if(incx.eq.1.and.incy.eq.1)go to 20
c
c        code for unequal increments or equal increments
c          not equal to 1
c
      ix = 1
      iy = 1
      if(incx.lt.0)ix = (-n+1)*incx + 1
      if(incy.lt.0)iy = (-n+1)*incy + 1
      do 10 i = 1,n
        dtemp = dtemp + dx(ix)*dy(iy)
        ix = ix + incx
        iy = iy + incy
   10 continue
      ddot = dtemp
      return
c
c        code for both increments equal to 1
c
c
c        clean-up loop
c
   20 m = mod(n,5)
      if( m .eq. 0 ) go to 40
      do 30 i = 1,m
        dtemp = dtemp + dx(i)*dy(i)
   30 continue
      if( n .lt. 5 ) go to 60
   40 mp1 = m + 1
      do 50 i = mp1,n,5
        dtemp = dtemp + dx(i)*dy(i) + dx(i + 1)*dy(i + 1) +
     *   dx(i + 2)*dy(i + 2) + dx(i + 3)*dy(i + 3) + dx(i + 4)*dy(i + 4)
   50 continue
   60 ddot = dtemp
      return
      end
      subroutine dgccsv ( x, y )

c*********************************************************************72
c
cc DGCCSV computes y = A * x for A in compressed column storage (CCS).
c
*
*---  Computes y = A*x. A is passed via a fortran `common statement'.
*---  A is assumed to be under the Compress Column Storage (CCS) format.
*
      implicit none

      double precision x(*), y(*)

      integer n, nz, nzmax
      parameter( nzmax = 600000 )
      integer ia(nzmax), ja(nzmax)
      double precision a(nzmax)
      common /RMAT/ a, ia, ja, nz, n
      integer i, j

      do i = 1,n
         y(i) = 0.0d0
      end do
      do j = 1,n
         do i = ja(j),ja(j+1)-1
            y(ia(i)) = y(ia(i)) + a(i)*x(j)
         end do
      end do
      end
      subroutine dgchbv ( m, t, h,ldh, y, wsp, iwsp, iflag )

c*********************************************************************72
c
cc DGCHBV computes y=exp(t*H)*y, general matrix H, using rational Chebyshev.
c
*-----Purpose----------------------------------------------------------|
*
*---  DGCHBV computes y = exp(t*H)*y using the partial fraction
*     expansion of the uniform rational Chebyshev approximation
*     to exp(-x) of type (14,14). H is a General matrix.
*     About 14-digit accuracy is expected if the matrix H is negative
*     definite. The algorithm may behave poorly otherwise.
*
*-----Arguments--------------------------------------------------------|
*
*     m       : (input) order of the matrix H
*
*     t       : (input) time-scaling factor (can be < 0).
*
*     H(ldh,m): (input) argument matrix.
*
*     y(m)    : (input/output) on input the operand vector,
*               on output the resulting vector exp(t*H)*y.
*
*     iwsp(m) : (workspace)
*
*     wsp     : (workspace). Observe that a double precision vector of
*               length 2*m*(m+2) can be used as well when calling this
*               routine (thus avoiding an idle complex array elsewhere)
*
*----------------------------------------------------------------------|
*     Roger B. Sidje (rbs@maths.uq.edu.au)
*     EXPOKIT: Software Package for Computing Matrix Exponentials.
*     ACM - Transactions On Mathematical Software, 24(1):130-156, 1998
*----------------------------------------------------------------------|
*
      implicit none

      integer          m, ldh, iflag, iwsp(m)
      double precision t, H(ldh,m), y(m)
      complex*16       wsp(m*(m+2))

      integer ndeg, i, j, ip, ih, iy, iz
      parameter ( ndeg=7 )
      double precision alpha0
      complex*16 alpha(ndeg), theta(ndeg)

      intrinsic DBLE

*---  Pointers ...

      ih = 1
      iy = ih + m*m
      iz = iy + m

*---  Coefficients and poles of the partial fraction expansion ...

      alpha0  =  0.183216998528140087D-11
      alpha(1)=( 0.557503973136501826D+02,-0.204295038779771857D+03)
      alpha(2)=(-0.938666838877006739D+02, 0.912874896775456363D+02)
      alpha(3)=( 0.469965415550370835D+02,-0.116167609985818103D+02)
      alpha(4)=(-0.961424200626061065D+01,-0.264195613880262669D+01)
      alpha(5)=( 0.752722063978321642D+00, 0.670367365566377770D+00)
      alpha(6)=(-0.188781253158648576D-01,-0.343696176445802414D-01)
      alpha(7)=( 0.143086431411801849D-03, 0.287221133228814096D-03)

      theta(1)=(-0.562314417475317895D+01, 0.119406921611247440D+01)
      theta(2)=(-0.508934679728216110D+01, 0.358882439228376881D+01)
      theta(3)=(-0.399337136365302569D+01, 0.600483209099604664D+01)
      theta(4)=(-0.226978543095856366D+01, 0.846173881758693369D+01)
      theta(5)=( 0.208756929753827868D+00, 0.109912615662209418D+02)
      theta(6)=( 0.370327340957595652D+01, 0.136563731924991884D+02)
      theta(7)=( 0.889777151877331107D+01, 0.166309842834712071D+02)
*
*---  Accumulation of the contribution of each pole ...
*
      do j = 1,m
         wsp(iz+j-1) = y(j)
         y(j) = y(j)*alpha0
      end do
      do ip = 1,ndeg
*---     Solve each fraction using Gaussian elimination with pivoting...
         do j = 1,m
            do i = 1,m
               wsp(ih+(j-1)*m+i-1) = -t*H(i,j)
            end do
            wsp(ih+(j-1)*m+j-1) = wsp(ih+(j-1)*m+j-1)-theta(ip)
            wsp(iy+j-1) = wsp(iz+j-1)
         end do
         call ZGESV( M, 1, WSP(iH),M, IWSP, WSP(iY),M, IFLAG )
         if ( IFLAG.ne.0 ) stop 'Error in DGCHBV'
c
*---     Accumulate the partial result in y.
c
         do j = 1,m
            y(j) = y(j) + DBLE( alpha(ip)*wsp(iy+j-1) )
         end do
      end do
      END
      subroutine dgcnvr ( from, to, diag, nrow, ncol, nz, ia, ja, a,
     &  iwsp )

c*********************************************************************72
c
cc DGCNVR converts between sparse matrix storage formats.
c
*-----Purpose----------------------------------------------------------|
*
*---  DGCNVR converts a sparse storage format into another sparse
*     storage format. The matrix can be rectangular.
*
*-----Arguments--------------------------------------------------------|
*
*     from   : (input, character*3) the storage format holding the
*              matrix on input. Accepted values are:
*              'COO' : COOrdinates
*              'CRS' : Compressed Row Storage
*              'CCS' : Compressed Column Storage.
*
*     to     : (input, character*3) the storage format holding the
*              matrix on output. Same accepted values as above.
*
*     diag   : (input, character*1) specifies whether or not the
*              entire diagonal should be lodged, i.e., including
*              null diagonal entries. (This may be required by
*              some applications that make use of the locations
*              of diagonal elements.)
*              diag = 'N', no specific attention is given to diagonal
*                          elements.
*              diag = 'D', the entire diagonal is lodged, including
*                          null elements on the diagonal.
*              if from=to & diag='D', null diagonal entries are
*              explicitly inserted in the actual matrix.
*
*     nrow   : (input) number of rows in the matrix.
*
*     ncol   : (input) number of columns in the matrix.
*
*     nz     : (input/output) number of non-zeros elements.
*              If diag='D' and null diagonal entries are inserted, then
*              nz is updated on exit and contains the effective number
*              of entries stored. In what follows, nz' (read nz prime)
*              denotes the updated value of nz, nz <= nz' <= nz+n.
*              If diag='N' then nz'=nz.
*
*     ia(*) : (input/output) of declared length .ge. nz'.
*             On input,
*                if from='CRS', ia(1:nrow+1) contains pointers for the
*                beginning of each row.
*                if from='COO', or 'CCS', ia(1:nz) contains row indices.
*             On output,
*                if to='CRS', ia(1:nrow+1) contains pointers for the
*                beginning of each row.
*                if to='COO' or 'CCS', ia(1:nz') contains row indices
*                in increasing order in each column.
*
*     ja(*) : (input/output) of declared length .ge. nz'.
*             On input,
*                if from='CRS', ja(1:ncol+1) contains pointers for the
*                beginning of each column.
*                if from='COO' or 'CCS', ja(1:nz) contains col. indices.
*             On output,
*                if to='CRS', ja(1:ncol+1) contains pointers for the
*                beginning of each column.
*                if to='COO' or 'CCS', ja(1:nz') contains column indices
*                in increasing order in each row.
*
*     a(*)  : (input/output) On input, a(1:nz) fits the input format
*             and on output, a(1:nz') fits the output format.
*
*   iwsp(*) : (workspace) of declared length .ge. max(nrow,ncol).
*
*----------------------------------------------------------------------|
*     Roger B. Sidje (rbs@maths.uq.edu.au)
*     Department of Mathematics, University of Queensland.
*     Brisbane QLD 4072, Australia. 1996.
*----------------------------------------------------------------------|

      implicit none

      character        from*3, to*3, diag*1
      integer          nrow, ncol, nz, ia(*), ja(*), iwsp(*)
      double precision a(*)


      integer i, j, k, nn, iflag
      character cpy_from*3, cpy_to*3, cpy_diag*1, c*1, S1*3, S2*3
      logical LSAME

      intrinsic CHAR,ICHAR,LLE,LGE,MIN
      LSAME(S1,S2) = LLE(S1,S2).and.LGE(S1,S2)
*
*---  upper case strings ...
*
      do k = 1,3
         c = from(k:k)
         if ( c.ge.'a' .and. c.le.'z' ) c = CHAR( ICHAR(c)-32 )
         cpy_from(k:k) = c
         c = to(k:k)
         if ( c.ge.'a' .and. c.le.'z' ) c = CHAR( ICHAR(c)-32 )
         cpy_to(k:k) = c
      end do
      c = diag(1:1)
      if ( c.ge.'a' .and. c.le.'z' ) c = CHAR( ICHAR(c)-32 )
      cpy_diag = c

*---  quick return if possible ...
      if ( LSAME(cpy_from,cpy_to) .and. cpy_diag.eq.'N' ) return
      iflag = 1
      if ( LSAME(cpy_from,'COO') ) iflag = 0
      if ( LSAME(cpy_from,'CRS') ) iflag = 0
      if ( LSAME(cpy_from,'CCS') ) iflag = 0
      if ( iflag.ne.0 ) stop 'unexpected i/o formats (in DGCNVR)'
      iflag = 1
      if ( LSAME(cpy_to,'COO') )   iflag = 0
      if ( LSAME(cpy_to,'CRS') )   iflag = 0
      if ( LSAME(cpy_to,'CCS') )   iflag = 0
      if ( iflag.ne.0 ) stop 'unexpected i/o formats (in DGCNVR)'

*
*---  transit via COOrdinates format if input is not in COO ...
*
      if ( LSAME(cpy_from,'CRS') ) then
*---     expand ia indices ...
         nn = nz
         do i = nrow,1,-1
            do k = 1,ia(i+1)-ia(i)
               ia(nn) = i
               nn = nn-1
            end do
         end do
      endif
      if ( LSAME(cpy_from,'CCS') ) then
*---     expand ja indices ...
         nn = nz
         do j = ncol,1,-1
            do k = 1,ja(j+1)-ja(j)
               ja(nn) = j
               nn = nn-1
            end do
         end do
      endif
*
*--   if requested, insert diagonal elements even if they are zero...
*
      if ( cpy_diag.eq.'D' ) then
         nn = MIN( nrow,ncol )
         do i = 1,nn
            iwsp(i) = 0
         end do
         do k = 1,nz
            if ( ia(k).eq.ja(k) ) iwsp(ia(k)) = 1
         end do
         do i = 1,nn
            if ( iwsp(i).eq.0 ) then
               nz = nz + 1
               ia(nz) = i
               ja(nz) = i
               a(nz) = 0.0d0
            endif
         end do
      endif

*---  COO convertion ...
      if ( LSAME(cpy_to,'COO') ) return

*---  CRS convertion ...
      if ( LSAME(cpy_to,'CRS') ) call dcmpac( nrow,nz,ia,ja,a, iwsp )

*---  CCS convertion ...
      if ( LSAME(cpy_to,'CCS') ) call dcmpac( ncol,nz,ja,ia,a, iwsp )

      END
      subroutine dgcoov ( x, y )

c*********************************************************************72
c
cc DGCOOV computes y = A * x, with A stored in "COO" format.
c
*
*---  Computes y = A*x. A is passed via a fortran `common statement'.
*---  A is assumed here to be under the COOrdinates storage format.
*
*----------------------------------------------------------------------|
*-------------------------------NOTE-----------------------------------*
*     This is an accessory to Expokit and it is not intended to be     *
*     complete. It is supplied primarily to ensure an unconstrained    *
*     distribution and portability of the package. The matrix-vector   *
*     multiplication routines supplied here fit the non symmetric      *
*     storage and for a symmetric matrix, the entire (not half) matrix *
*     is required.  If the sparsity pattern is known a priori, it is   *
*     recommended to use the most advantageous format and to devise    *
*     the most advantageous matrix-vector multiplication routine.      *
*----------------------------------------------------------------------*
*----------------------------------------------------------------------*
      implicit none

      double precision x(*), y(*)

      integer n, nz, nzmax
      parameter( nzmax = 600000 )
      integer ia(nzmax), ja(nzmax)
      double precision a(nzmax)
      common /RMAT/ a, ia, ja, nz, n
      integer i, j

      do j = 1,n
         y(j) = 0.0d0
      end do
      do i = 1,nz
         y(ia(i)) = y(ia(i)) + a(i)*x(ja(i))
      end do
      END
      subroutine dgcrsv ( x, y )

c*********************************************************************72
c
cc DGCRSV computes y = A * x, for A stored in Compressed Row Storage (CRS) format.
c
      implicit none

      double precision x(*), y(*)
*
*---  Computes y = A*x. A is passed via a fortran `common statement'.
*---  A is assumed to be under the Compress Row Storage (CRS) format.
*
      integer n, nz, nzmax
      parameter( nzmax = 600000 )
      integer ia(nzmax), ja(nzmax)
      double precision a(nzmax)
      common /RMAT/ a, ia, ja, nz, n
      integer i, j

      do i = 1,n
         y(i) = 0.0d0
         do j = ia(i),ia(i+1)-1
            y(i) = y(i) + a(j)*x(ja(j))
         end do
      end do
      END
      subroutine dgefa ( a, lda, n, ipvt, info )

c*********************************************************************72
c
cc DGEFA factors a double precision matrix by gaussian elimination.
c
      implicit none

      integer lda,n,ipvt(n),info
      double precision a(lda,n)
c
c     dgefa is usually called by dgeco, but it can be called
c     directly with a saving in time if  rcond  is not needed.
c     (time for dgeco) = (1 + 9/n)*(time for dgefa) .
c
c     on entry
c
c        a       double precision(lda, n)
c                the matrix to be factored.
c
c        lda     integer
c                the leading dimension of the array  a .
c
c        n       integer
c                the order of the matrix  a .
c
c     on return
c
c        a       an upper triangular matrix and the multipliers
c                which were used to obtain it.
c                the factorization can be written  a = l*u  where
c                l  is a product of permutation and unit lower
c                triangular matrices and  u  is upper triangular.
c
c        ipvt    integer(n)
c                an integer vector of pivot indices.
c
c        info    integer
c                = 0  normal value.
c                = k  if  u(k,k) .eq. 0.0 .  this is not an error
c                     condition for this routine, but it does
c                     indicate that dgesl or dgedi will divide by zero
c                     if called.  use  rcond  in dgeco for a reliable
c                     indication of singularity.
c
c     linpack. this version dated 08/14/78 .
c     cleve moler, university of new mexico, argonne national lab.
c
c     subroutines and functions
c
c     blas daxpy,dscal,idamax
c
c     internal variables
c
      double precision t
      integer idamax,j,k,kp1,l,nm1
c
c
c     gaussian elimination with partial pivoting
c
      info = 0
      nm1 = n - 1
      if (nm1 .lt. 1) go to 70
      do 60 k = 1, nm1
         kp1 = k + 1
c
c        find l = pivot index
c
         l = idamax(n-k+1,a(k,k),1) + k - 1
         ipvt(k) = l
c
c        zero pivot implies this column already triangularized
c
         if (a(l,k) .eq. 0.0d0) go to 40
c
c           interchange if necessary
c
            if (l .eq. k) go to 10
               t = a(l,k)
               a(l,k) = a(k,k)
               a(k,k) = t
   10       continue
c
c           compute multipliers
c
            t = -1.0d0/a(k,k)
            call dscal(n-k,t,a(k+1,k),1)
c
c           row elimination with column indexing
c
            do 30 j = kp1, n
               t = a(l,j)
               if (l .eq. k) go to 20
                  a(l,j) = a(k,j)
                  a(k,j) = t
   20          continue
               call daxpy(n-k,t,a(k+1,k),1,a(k+1,j),1)
   30       continue
         go to 50
   40    continue
            info = k
   50    continue
   60 continue
   70 continue
      ipvt(n) = n
      if (a(n,n) .eq. 0.0d0) info = n
      return
      end
      subroutine dgemm ( transa, transb, m, n, k, alpha, a, lda, b, ldb,
     &  beta, c, ldc )

c*********************************************************************72
c
cc DGEMM performs matrix-matrix operations of the form C=alpha*A*B+beta*C.
c
      implicit none

      CHARACTER*1        TRANSA, TRANSB
      INTEGER            M, N, K, LDA, LDB, LDC
      DOUBLE PRECISION   ALPHA, BETA
*     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * ), B( LDB, * ), C( LDC, * )
*     ..
*
*  Purpose
*  =======
*
*  DGEMM  performs one of the matrix-matrix operations
*
*     C := alpha*op( A )*op( B ) + beta*C,
*
*  where  op( X ) is one of
*
*     op( X ) = X   or   op( X ) = X',
*
*  alpha and beta are scalars, and A, B and C are matrices, with op( A )
*  an m by k matrix,  op( B )  a  k by n matrix and  C an m by n matrix.
*
*  Parameters
*  ==========
*
*  TRANSA - CHARACTER*1.
*           On entry, TRANSA specifies the form of op( A ) to be used in
*           the matrix multiplication as follows:
*
*              TRANSA = 'N' or 'n',  op( A ) = A.
*
*              TRANSA = 'T' or 't',  op( A ) = A'.
*
*              TRANSA = 'C' or 'c',  op( A ) = A'.
*
*           Unchanged on exit.
*
*  TRANSB - CHARACTER*1.
*           On entry, TRANSB specifies the form of op( B ) to be used in
*           the matrix multiplication as follows:
*
*              TRANSB = 'N' or 'n',  op( B ) = B.
*
*              TRANSB = 'T' or 't',  op( B ) = B'.
*
*              TRANSB = 'C' or 'c',  op( B ) = B'.
*
*           Unchanged on exit.
*
*  M      - INTEGER.
*           On entry,  M  specifies  the number  of rows  of the  matrix
*           op( A )  and of the  matrix  C.  M  must  be at least  zero.
*           Unchanged on exit.
*
*  N      - INTEGER.
*           On entry,  N  specifies the number  of columns of the matrix
*           op( B ) and the number of columns of the matrix C. N must be
*           at least zero.
*           Unchanged on exit.
*
*  K      - INTEGER.
*           On entry,  K  specifies  the number of columns of the matrix
*           op( A ) and the number of rows of the matrix op( B ). K must
*           be at least  zero.
*           Unchanged on exit.
*
*  ALPHA  - DOUBLE PRECISION.
*           On entry, ALPHA specifies the scalar alpha.
*           Unchanged on exit.
*
*  A      - DOUBLE PRECISION array of DIMENSION ( LDA, ka ), where ka is
*           k  when  TRANSA = 'N' or 'n',  and is  m  otherwise.
*           Before entry with  TRANSA = 'N' or 'n',  the leading  m by k
*           part of the array  A  must contain the matrix  A,  otherwise
*           the leading  k by m  part of the array  A  must contain  the
*           matrix A.
*           Unchanged on exit.
*
*  LDA    - INTEGER.
*           On entry, LDA specifies the first dimension of A as declared
*           in the calling (sub) program. When  TRANSA = 'N' or 'n' then
*           LDA must be at least  max( 1, m ), otherwise  LDA must be at
*           least  max( 1, k ).
*           Unchanged on exit.
*
*  B      - DOUBLE PRECISION array of DIMENSION ( LDB, kb ), where kb is
*           n  when  TRANSB = 'N' or 'n',  and is  k  otherwise.
*           Before entry with  TRANSB = 'N' or 'n',  the leading  k by n
*           part of the array  B  must contain the matrix  B,  otherwise
*           the leading  n by k  part of the array  B  must contain  the
*           matrix B.
*           Unchanged on exit.
*
*  LDB    - INTEGER.
*           On entry, LDB specifies the first dimension of B as declared
*           in the calling (sub) program. When  TRANSB = 'N' or 'n' then
*           LDB must be at least  max( 1, k ), otherwise  LDB must be at
*           least  max( 1, n ).
*           Unchanged on exit.
*
*  BETA   - DOUBLE PRECISION.
*           On entry,  BETA  specifies the scalar  beta.  When  BETA  is
*           supplied as zero then C need not be set on input.
*           Unchanged on exit.
*
*  C      - DOUBLE PRECISION array of DIMENSION ( LDC, n ).
*           Before entry, the leading  m by n  part of the array  C must
*           contain the matrix  C,  except when  beta  is zero, in which
*           case C need not be set on entry.
*           On exit, the array  C  is overwritten by the  m by n  matrix
*           ( alpha*op( A )*op( B ) + beta*C ).
*
*  LDC    - INTEGER.
*           On entry, LDC specifies the first dimension of C as declared
*           in  the  calling  (sub)  program.   LDC  must  be  at  least
*           max( 1, m ).
*           Unchanged on exit.
*
*
*  Level 3 Blas routine.
*
*  -- Written on 8-February-1989.
*     Jack Dongarra, Argonne National Laboratory.
*     Iain Duff, AERE Harwell.
*     Jeremy Du Croz, Numerical Algorithms Group Ltd.
*     Sven Hammarling, Numerical Algorithms Group Ltd.
*
*
*     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
*     .. External Subroutines ..
      EXTERNAL           XERBLA
*     .. Intrinsic Functions ..
      INTRINSIC          MAX
*     .. Local Scalars ..
      LOGICAL            NOTA, NOTB
      INTEGER            I, INFO, J, L, NCOLA, NROWA, NROWB
      DOUBLE PRECISION   TEMP
*     .. Parameters ..
      DOUBLE PRECISION   ONE         , ZERO
      PARAMETER        ( ONE = 1.0D+0, ZERO = 0.0D+0 )
*     ..
*     .. Executable Statements ..
*
*     Set  NOTA  and  NOTB  as  true if  A  and  B  respectively are not
*     transposed and set  NROWA, NCOLA and  NROWB  as the number of rows
*     and  columns of  A  and the  number of  rows  of  B  respectively.
*
      NOTA  = LSAME( TRANSA, 'N' )
      NOTB  = LSAME( TRANSB, 'N' )
      IF( NOTA )THEN
         NROWA = M
         NCOLA = K
      ELSE
         NROWA = K
         NCOLA = M
      END IF
      IF( NOTB )THEN
         NROWB = K
      ELSE
         NROWB = N
      END IF
*
*     Test the input parameters.
*
      INFO = 0
      IF(      ( .NOT.NOTA                 ).AND.
     $         ( .NOT.LSAME( TRANSA, 'C' ) ).AND.
     $         ( .NOT.LSAME( TRANSA, 'T' ) )      )THEN
         INFO = 1
      ELSE IF( ( .NOT.NOTB                 ).AND.
     $         ( .NOT.LSAME( TRANSB, 'C' ) ).AND.
     $         ( .NOT.LSAME( TRANSB, 'T' ) )      )THEN
         INFO = 2
      ELSE IF( M  .LT.0               )THEN
         INFO = 3
      ELSE IF( N  .LT.0               )THEN
         INFO = 4
      ELSE IF( K  .LT.0               )THEN
         INFO = 5
      ELSE IF( LDA.LT.MAX( 1, NROWA ) )THEN
         INFO = 8
      ELSE IF( LDB.LT.MAX( 1, NROWB ) )THEN
         INFO = 10
      ELSE IF( LDC.LT.MAX( 1, M     ) )THEN
         INFO = 13
      END IF
      IF( INFO.NE.0 )THEN
         CALL XERBLA( 'DGEMM ', INFO )
         RETURN
      END IF
*
*     Quick return if possible.
*
      IF( ( M.EQ.0 ).OR.( N.EQ.0 ).OR.
     $    ( ( ( ALPHA.EQ.ZERO ).OR.( K.EQ.0 ) ).AND.( BETA.EQ.ONE ) ) )
     $   RETURN
*
*     And if  alpha.eq.zero.
*
      IF( ALPHA.EQ.ZERO )THEN
         IF( BETA.EQ.ZERO )THEN
            DO 20, J = 1, N
               DO 10, I = 1, M
                  C( I, J ) = ZERO
   10          CONTINUE
   20       CONTINUE
         ELSE
            DO 40, J = 1, N
               DO 30, I = 1, M
                  C( I, J ) = BETA*C( I, J )
   30          CONTINUE
   40       CONTINUE
         END IF
         RETURN
      END IF
*
*     Start the operations.
*
      IF( NOTB )THEN
         IF( NOTA )THEN
*
*           Form  C := alpha*A*B + beta*C.
*
            DO 90, J = 1, N
               IF( BETA.EQ.ZERO )THEN
                  DO 50, I = 1, M
                     C( I, J ) = ZERO
   50             CONTINUE
               ELSE IF( BETA.NE.ONE )THEN
                  DO 60, I = 1, M
                     C( I, J ) = BETA*C( I, J )
   60             CONTINUE
               END IF
               DO 80, L = 1, K
                  IF( B( L, J ).NE.ZERO )THEN
                     TEMP = ALPHA*B( L, J )
                     DO 70, I = 1, M
                        C( I, J ) = C( I, J ) + TEMP*A( I, L )
   70                CONTINUE
                  END IF
   80          CONTINUE
   90       CONTINUE
         ELSE
*
*           Form  C := alpha*A'*B + beta*C
*
            DO 120, J = 1, N
               DO 110, I = 1, M
                  TEMP = ZERO
                  DO 100, L = 1, K
                     TEMP = TEMP + A( L, I )*B( L, J )
  100             CONTINUE
                  IF( BETA.EQ.ZERO )THEN
                     C( I, J ) = ALPHA*TEMP
                  ELSE
                     C( I, J ) = ALPHA*TEMP + BETA*C( I, J )
                  END IF
  110          CONTINUE
  120       CONTINUE
         END IF
      ELSE
         IF( NOTA )THEN
*
*           Form  C := alpha*A*B' + beta*C
*
            DO 170, J = 1, N
               IF( BETA.EQ.ZERO )THEN
                  DO 130, I = 1, M
                     C( I, J ) = ZERO
  130             CONTINUE
               ELSE IF( BETA.NE.ONE )THEN
                  DO 140, I = 1, M
                     C( I, J ) = BETA*C( I, J )
  140             CONTINUE
               END IF
               DO 160, L = 1, K
                  IF( B( J, L ).NE.ZERO )THEN
                     TEMP = ALPHA*B( J, L )
                     DO 150, I = 1, M
                        C( I, J ) = C( I, J ) + TEMP*A( I, L )
  150                CONTINUE
                  END IF
  160          CONTINUE
  170       CONTINUE
         ELSE
*
*           Form  C := alpha*A'*B' + beta*C
*
            DO 200, J = 1, N
               DO 190, I = 1, M
                  TEMP = ZERO
                  DO 180, L = 1, K
                     TEMP = TEMP + A( L, I )*B( J, L )
  180             CONTINUE
                  IF( BETA.EQ.ZERO )THEN
                     C( I, J ) = ALPHA*TEMP
                  ELSE
                     C( I, J ) = ALPHA*TEMP + BETA*C( I, J )
                  END IF
  190          CONTINUE
  200       CONTINUE
         END IF
      END IF
*
      RETURN
      END
      subroutine dgemv ( trans, m, n, alpha, a, lda, x, incx,
     &  beta, y, incy )

c*********************************************************************72
c
cc DGEMV performs a matrix-vector operation of the form y=alpha*A*x+beta*y.
c
      implicit none

      DOUBLE PRECISION   ALPHA, BETA
      INTEGER            INCX, INCY, LDA, M, N
      CHARACTER*1        TRANS
*     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * ), X( * ), Y( * )
*     ..
*
*  Purpose
*  =======
*
*  DGEMV  performs one of the matrix-vector operations
*
*     y := alpha*A*x + beta*y,   or   y := alpha*A'*x + beta*y,
*
*  where alpha and beta are scalars, x and y are vectors and A is an
*  m by n matrix.
*
*  Parameters
*  ==========
*
*  TRANS  - CHARACTER*1.
*           On entry, TRANS specifies the operation to be performed as
*           follows:
*
*              TRANS = 'N' or 'n'   y := alpha*A*x + beta*y.
*
*              TRANS = 'T' or 't'   y := alpha*A'*x + beta*y.
*
*              TRANS = 'C' or 'c'   y := alpha*A'*x + beta*y.
*
*           Unchanged on exit.
*
*  M      - INTEGER.
*           On entry, M specifies the number of rows of the matrix A.
*           M must be at least zero.
*           Unchanged on exit.
*
*  N      - INTEGER.
*           On entry, N specifies the number of columns of the matrix A.
*           N must be at least zero.
*           Unchanged on exit.
*
*  ALPHA  - DOUBLE PRECISION.
*           On entry, ALPHA specifies the scalar alpha.
*           Unchanged on exit.
*
*  A      - DOUBLE PRECISION array of DIMENSION ( LDA, n ).
*           Before entry, the leading m by n part of the array A must
*           contain the matrix of coefficients.
*           Unchanged on exit.
*
*  LDA    - INTEGER.
*           On entry, LDA specifies the first dimension of A as declared
*           in the calling (sub) program. LDA must be at least
*           max( 1, m ).
*           Unchanged on exit.
*
*  X      - DOUBLE PRECISION array of DIMENSION at least
*           ( 1 + ( n - 1 )*abs( INCX ) ) when TRANS = 'N' or 'n'
*           and at least
*           ( 1 + ( m - 1 )*abs( INCX ) ) otherwise.
*           Before entry, the incremented array X must contain the
*           vector x.
*           Unchanged on exit.
*
*  INCX   - INTEGER.
*           On entry, INCX specifies the increment for the elements of
*           X. INCX must not be zero.
*           Unchanged on exit.
*
*  BETA   - DOUBLE PRECISION.
*           On entry, BETA specifies the scalar beta. When BETA is
*           supplied as zero then Y need not be set on input.
*           Unchanged on exit.
*
*  Y      - DOUBLE PRECISION array of DIMENSION at least
*           ( 1 + ( m - 1 )*abs( INCY ) ) when TRANS = 'N' or 'n'
*           and at least
*           ( 1 + ( n - 1 )*abs( INCY ) ) otherwise.
*           Before entry with BETA non-zero, the incremented array Y
*           must contain the vector y. On exit, Y is overwritten by the
*           updated vector y.
*
*  INCY   - INTEGER.
*           On entry, INCY specifies the increment for the elements of
*           Y. INCY must not be zero.
*           Unchanged on exit.
*
*
*  Level 2 Blas routine.
*
*  -- Written on 22-October-1986.
*     Jack Dongarra, Argonne National Lab.
*     Jeremy Du Croz, Nag Central Office.
*     Sven Hammarling, Nag Central Office.
*     Richard Hanson, Sandia National Labs.
*
*
*     .. Parameters ..
      DOUBLE PRECISION   ONE         , ZERO
      PARAMETER        ( ONE = 1.0D+0, ZERO = 0.0D+0 )
*     .. Local Scalars ..
      DOUBLE PRECISION   TEMP
      INTEGER            I, INFO, IX, IY, J, JX, JY, KX, KY, LENX, LENY
*     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
*     .. External Subroutines ..
      EXTERNAL           XERBLA
*     .. Intrinsic Functions ..
      INTRINSIC          MAX
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters.
*
      INFO = 0
      IF     ( .NOT.LSAME( TRANS, 'N' ).AND.
     $         .NOT.LSAME( TRANS, 'T' ).AND.
     $         .NOT.LSAME( TRANS, 'C' )      )THEN
         INFO = 1
      ELSE IF( M.LT.0 )THEN
         INFO = 2
      ELSE IF( N.LT.0 )THEN
         INFO = 3
      ELSE IF( LDA.LT.MAX( 1, M ) )THEN
         INFO = 6
      ELSE IF( INCX.EQ.0 )THEN
         INFO = 8
      ELSE IF( INCY.EQ.0 )THEN
         INFO = 11
      END IF
      IF( INFO.NE.0 )THEN
         CALL XERBLA( 'DGEMV ', INFO )
         RETURN
      END IF
*
*     Quick return if possible.
*
      IF( ( M.EQ.0 ).OR.( N.EQ.0 ).OR.
     $    ( ( ALPHA.EQ.ZERO ).AND.( BETA.EQ.ONE ) ) )
     $   RETURN
*
*     Set  LENX  and  LENY, the lengths of the vectors x and y, and set
*     up the start points in  X  and  Y.
*
      IF( LSAME( TRANS, 'N' ) )THEN
         LENX = N
         LENY = M
      ELSE
         LENX = M
         LENY = N
      END IF
      IF( INCX.GT.0 )THEN
         KX = 1
      ELSE
         KX = 1 - ( LENX - 1 )*INCX
      END IF
      IF( INCY.GT.0 )THEN
         KY = 1
      ELSE
         KY = 1 - ( LENY - 1 )*INCY
      END IF
*
*     Start the operations. In this version the elements of A are
*     accessed sequentially with one pass through A.
*
*     First form  y := beta*y.
*
      IF( BETA.NE.ONE )THEN
         IF( INCY.EQ.1 )THEN
            IF( BETA.EQ.ZERO )THEN
               DO 10, I = 1, LENY
                  Y( I ) = ZERO
   10          CONTINUE
            ELSE
               DO 20, I = 1, LENY
                  Y( I ) = BETA*Y( I )
   20          CONTINUE
            END IF
         ELSE
            IY = KY
            IF( BETA.EQ.ZERO )THEN
               DO 30, I = 1, LENY
                  Y( IY ) = ZERO
                  IY      = IY   + INCY
   30          CONTINUE
            ELSE
               DO 40, I = 1, LENY
                  Y( IY ) = BETA*Y( IY )
                  IY      = IY           + INCY
   40          CONTINUE
            END IF
         END IF
      END IF
      IF( ALPHA.EQ.ZERO )
     $   RETURN
      IF( LSAME( TRANS, 'N' ) )THEN
*
*        Form  y := alpha*A*x + y.
*
         JX = KX
         IF( INCY.EQ.1 )THEN
            DO 60, J = 1, N
               IF( X( JX ).NE.ZERO )THEN
                  TEMP = ALPHA*X( JX )
                  DO 50, I = 1, M
                     Y( I ) = Y( I ) + TEMP*A( I, J )
   50             CONTINUE
               END IF
               JX = JX + INCX
   60       CONTINUE
         ELSE
            DO 80, J = 1, N
               IF( X( JX ).NE.ZERO )THEN
                  TEMP = ALPHA*X( JX )
                  IY   = KY
                  DO 70, I = 1, M
                     Y( IY ) = Y( IY ) + TEMP*A( I, J )
                     IY      = IY      + INCY
   70             CONTINUE
               END IF
               JX = JX + INCX
   80       CONTINUE
         END IF
      ELSE
*
*        Form  y := alpha*A'*x + y.
*
         JY = KY
         IF( INCX.EQ.1 )THEN
            DO 100, J = 1, N
               TEMP = ZERO
               DO 90, I = 1, M
                  TEMP = TEMP + A( I, J )*X( I )
   90          CONTINUE
               Y( JY ) = Y( JY ) + ALPHA*TEMP
               JY      = JY      + INCY
  100       CONTINUE
         ELSE
            DO 120, J = 1, N
               TEMP = ZERO
               IX   = KX
               DO 110, I = 1, M
                  TEMP = TEMP + A( I, J )*X( IX )
                  IX   = IX   + INCX
  110          CONTINUE
               Y( JY ) = Y( JY ) + ALPHA*TEMP
               JY      = JY      + INCY
  120       CONTINUE
         END IF
      END IF
*
      RETURN
*
*     End of DGEMV .
*
      END
      subroutine dgesl ( a, lda, n, ipvt, b, job )

c*********************************************************************72
c
cc DGESL solves a factored linear system of the form A*x=b.
c
      implicit none

      integer lda,n,ipvt(n),job
      double precision a(lda,n),b(n)
c
c     dgesl solves the double precision system
c     a * x = b  or  trans(a) * x = b
c     using the factors computed by dgeco or dgefa.
c
c     on entry
c
c        a       double precision(lda, n)
c                the output from dgeco or dgefa.
c
c        lda     integer
c                the leading dimension of the array  a .
c
c        n       integer
c                the order of the matrix  a .
c
c        ipvt    integer(n)
c                the pivot vector from dgeco or dgefa.
c
c        b       double precision(n)
c                the right hand side vector.
c
c        job     integer
c                = 0         to solve  a*x = b ,
c                = nonzero   to solve  trans(a)*x = b  where
c                            trans(a)  is the transpose.
c
c     on return
c
c        b       the solution vector  x .
c
c     error condition
c
c        a division by zero will occur if the input factor contains a
c        zero on the diagonal.  technically this indicates singularity
c        but it is often caused by improper arguments or improper
c        setting of lda .  it will not occur if the subroutines are
c        called correctly and if dgeco has set rcond .gt. 0.0
c        or dgefa has set info .eq. 0 .
c
c     to compute  inverse(a) * c  where  c  is a matrix
c     with  p  columns
c           call dgeco(a,lda,n,ipvt,rcond,z)
c           if (rcond is too small) go to ...
c           do 10 j = 1, p
c              call dgesl(a,lda,n,ipvt,c(1,j),0)
c        10 continue
c
c     linpack. this version dated 08/14/78 .
c     cleve moler, university of new mexico, argonne national lab.
c
c     subroutines and functions
c
c     blas daxpy,ddot
c
c     internal variables
c
      double precision ddot,t
      integer k,kb,l,nm1
c
      nm1 = n - 1
      if (job .ne. 0) go to 50
c
c        job = 0 , solve  a * x = b
c        first solve  l*y = b
c
         if (nm1 .lt. 1) go to 30
         do 20 k = 1, nm1
            l = ipvt(k)
            t = b(l)
            if (l .eq. k) go to 10
               b(l) = b(k)
               b(k) = t
   10       continue
            call daxpy(n-k,t,a(k+1,k),1,b(k+1),1)
   20    continue
   30    continue
c
c        now solve  u*x = y
c
         do 40 kb = 1, n
            k = n + 1 - kb
            b(k) = b(k)/a(k,k)
            t = -b(k)
            call daxpy(k-1,t,a(1,k),1,b(1),1)
   40    continue
      go to 100
   50 continue
c
c        job = nonzero, solve  trans(a) * x = b
c        first solve  trans(u)*y = b
c
         do 60 k = 1, n
            t = ddot(k-1,a(1,k),1,b(1),1)
            b(k) = (b(k) - t)/a(k,k)
   60    continue
c
c        now solve trans(l)*x = y
c
         if (nm1 .lt. 1) go to 90
         do 80 kb = 1, nm1
            k = n - kb
            b(k) = b(k) + ddot(n-k,a(k+1,k),1,b(k+1),1)
            l = ipvt(k)
            if (l .eq. k) go to 70
               t = b(l)
               b(l) = b(k)
               b(k) = t
   70       continue
   80    continue
   90    continue
  100 continue
      return
      end
      subroutine dgesv ( n, m, a,lda, ipiv, b,ldb, iflag )

c*********************************************************************72
c
cc DGESV solves multiple general linear systems of the form A*X=B.
c
*     This is a lightweight substitute to the external LAPACK routines
*     used by EXPOKIT. It is supplied to ensure that EXPOKIT is
*     self-contained and can still run if LAPACK is not yet installed
*     in your environement.
c
      implicit none

      integer N, M, LDA, LDB, IPIV(N), IFLAG
      double precision A(LDA,N), B(LDB,M)
      integer j

      call DGEFA( A,LDA, N, IPIV, IFLAG )
      if ( IFLAG.ne.0 ) stop "Error in DGESV (LU factorisation)"

      do j = 1,M
         call DGESL( A,LDA, N, IPIV,B(1,j), 0 )
      end do

      end
      subroutine dgexpv ( n, m, t, v, w, tol, anorm,
     .  wsp, lwsp, iwsp, liwsp, matvec, itrace, iflag )

c*********************************************************************72
c
cc DGEXPV computes w=exp(t*A)*v for a general matrix A.
c
      implicit none

      integer n, m, lwsp, liwsp, itrace, iflag, iwsp(liwsp)
      double precision t, tol, anorm, v(n), w(n), wsp(lwsp)
      external matvec

*-----Purpose----------------------------------------------------------|
*
*---  DGEXPV computes w = exp(t*A)*v - for a General matrix A.
*
*     It does not compute the matrix exponential in isolation but
*     instead, it computes directly the action of the exponential
*     operator on the operand vector. This way of doing so allows
*     for addressing large sparse problems.
*
*     The method used is based on Krylov subspace projection
*     techniques and the matrix under consideration interacts only
*     via the external routine `matvec' performing the matrix-vector
*     product (matrix-free method).
*
*-----Arguments--------------------------------------------------------|
*
*     n      : (input) order of the principal matrix A.
*
*     m      : (input) maximum size for the Krylov basis.
*
*     t      : (input) time at wich the solution is needed (can be < 0).
*
*     v(n)   : (input) given operand vector.
*
*     w(n)   : (output) computed approximation of exp(t*A)*v.
*
*     tol    : (input/output) the requested accuracy tolerance on w.
*              If on input tol=0.0d0 or tol is too small (tol.le.eps)
*              the internal value sqrt(eps) is used, and tol is set to
*              sqrt(eps) on output (`eps' denotes the machine epsilon).
*              (`Happy breakdown' is assumed if h(j+1,j) .le. anorm*tol)
*
*     anorm  : (input) an approximation of some norm of A.
*
*   wsp(lwsp): (workspace) lwsp .ge. n*(m+1)+n+(m+2)^2+4*(m+2)^2+ideg+1
*                                   +---------+-------+---------------+
*              (actually, ideg=6)        V        H      wsp for PADE
*
* iwsp(liwsp): (workspace) liwsp .ge. m+2
*
*     matvec : external routine for matrix-vector multiplication.
*              synopsis: matvec( x, y )
*                        double precision x(*), y(*)
*              computes: y(1:n) <- A*x(1:n)
*                        where A is the principal matrix.
*
*     itrace : (input) running mode. 0=silent, 1=print step-by-step info
*
*     iflag  : (output) exit flag.
*              <0 - bad input arguments
*               0 - no problem
*               1 - maximum number of steps reached without convergence
*               2 - requested tolerance was too high
*
*-----Accounts on the computation--------------------------------------|
*     Upon exit, an interested user may retrieve accounts on the
*     computations. They are located in wsp and iwsp as indicated below:
*
*     location  mnemonic                 description
*     -----------------------------------------------------------------|
*     iwsp(1) = nmult, number of matrix-vector multiplications used
*     iwsp(2) = nexph, number of Hessenberg matrix exponential evaluated
*     iwsp(3) = nscale, number of repeated squaring involved in Pade
*     iwsp(4) = nstep, number of integration steps used up to completion
*     iwsp(5) = nreject, number of rejected step-sizes
*     iwsp(6) = ibrkflag, set to 1 if `happy breakdown' and 0 otherwise
*     iwsp(7) = mbrkdwn, if `happy brkdown', basis-size when it occured
*     -----------------------------------------------------------------|
*     wsp(1)  = step_min, minimum step-size used during integration
*     wsp(2)  = step_max, maximum step-size used during integration
*     wsp(3)  = dummy
*     wsp(4)  = dummy
*     wsp(5)  = x_error, maximum among all local truncation errors
*     wsp(6)  = s_error, global sum of local truncation errors
*     wsp(7)  = tbrkdwn, if `happy breakdown', time when it occured
*     wsp(8)  = t_now, integration domain successfully covered
*     wsp(9)  = hump, i.e., max||exp(sA)||, s in [0,t] (or [t,0] if t<0)
*     wsp(10) = ||w||/||v||, scaled norm of the solution w.
*     -----------------------------------------------------------------|
*     The `hump' is a measure of the conditioning of the problem. The
*     matrix exponential is well-conditioned if hump = 1, whereas it is
*     poorly-conditioned if hump >> 1. However the solution can still be
*     relatively fairly accurate even when the hump is large (the hump
*     is an upper bound), especially when the hump and the scaled norm
*     of w [this is also computed and returned in wsp(10)] are of the
*     same order of magnitude (further details in reference below).
*
*----------------------------------------------------------------------|
*-----The following parameters may also be adjusted herein-------------|
*
      integer mxstep, mxreject, ideg
      double precision delta, gamma
      parameter( mxstep   = 1000,
     .           mxreject = 0,
     .           ideg     = 6,
     .           delta    = 1.2d0,
     .           gamma    = 0.9d0 )

*     mxstep  : maximum allowable number of integration steps.
*               The value 0 means an infinite number of steps.
*
*     mxreject: maximum allowable number of rejections at each step.
*               The value 0 means an infinite number of rejections.
*
*     ideg    : the Pade approximation of type (ideg,ideg) is used as
*               an approximation to exp(H). The value 0 switches to the
*               uniform rational Chebyshev approximation of type (14,14)
*
*     delta   : local truncation error `safety factor'
*
*     gamma   : stepsize `shrinking factor'
*
*----------------------------------------------------------------------|
*     Roger B. Sidje (rbs@maths.uq.edu.au)
*     EXPOKIT: Software Package for Computing Matrix Exponentials.
*     ACM - Transactions On Mathematical Software, 24(1):130-156, 1998
*----------------------------------------------------------------------|
*
      integer i, j, k1, mh, mx, iv, ih, j1v, ns, ifree, lfree, iexph,
     .        ireject,ibrkflag,mbrkdwn, nmult, nreject, nexph, nscale,
     .        nstep
      double precision sgn, t_out, tbrkdwn, step_min,step_max, err_loc,
     .                 s_error, x_error, t_now, t_new, t_step, t_old,
     .                 xm, beta, break_tol, p1, p2, p3, eps, rndoff,
     .                 vnorm, avnorm, hj1j, hij, hump, SQR1

      intrinsic AINT,ABS,DBLE,LOG10,MAX,MIN,NINT,SIGN,SQRT
      double precision DDOT, DNRM2

*---  check restrictions on input parameters ...
      iflag = 0
      if ( lwsp.lt.n*(m+2)+5*(m+2)**2+ideg+1 ) iflag = -1
      if ( liwsp.lt.m+2 ) iflag = -2
      if ( m.ge.n .or. m.le.0 ) iflag = -3
      if ( iflag.ne.0 ) stop 'bad sizes (in input of DGEXPV)'
*
*---  initialisations ...
*
      k1 = 2
      mh = m + 2
      iv = 1
      ih = iv + n*(m+1) + n
      ifree = ih + mh*mh
      lfree = lwsp - ifree + 1

      ibrkflag = 0
      mbrkdwn  = m
      nmult    = 0
      nreject  = 0
      nexph    = 0
      nscale   = 0

      t_out    = ABS( t )
      tbrkdwn  = 0.0d0
      step_min = t_out
      step_max = 0.0d0
      nstep    = 0
      s_error  = 0.0d0
      x_error  = 0.0d0
      t_now    = 0.0d0
      t_new    = 0.0d0

      p1 = 4.0d0/3.0d0
 1    p2 = p1 - 1.0d0
      p3 = p2 + p2 + p2
      eps = ABS( p3-1.0d0 )
      if ( eps.eq.0.0d0 ) go to 1
      if ( tol.le.eps ) tol = SQRT( eps )
      rndoff = eps*anorm

      break_tol = 1.0d-7
*>>>  break_tol = tol
*>>>  break_tol = anorm*tol

      sgn = SIGN( 1.0d0,t )
      call DCOPY( n, v,1, w,1 )
      beta = DNRM2( n, w,1 )
      vnorm = beta
      hump = beta
*
*---  obtain the very first stepsize ...
*
      SQR1 = SQRT( 0.1d0 )
      xm = 1.0d0/DBLE( m )
      p1 = tol*(((m+1)/2.72D0)**(m+1))*SQRT(2.0D0*3.14D0*(m+1))
      t_new = (1.0d0/anorm)*(p1/(4.0d0*beta*anorm))**xm
      p1 = 10.0d0**(NINT( LOG10( t_new )-SQR1 )-1)
      t_new = AINT( t_new/p1 + 0.55d0 ) * p1
*
*---  step-by-step integration ...
*
 100  if ( t_now.ge.t_out ) goto 500

      nstep = nstep + 1
      t_step = MIN( t_out-t_now, t_new )

      p1 = 1.0d0/beta
      do i = 1,n
         wsp(iv + i-1) = p1*w(i)
      end do
      do i = 1,mh*mh
         wsp(ih+i-1) = 0.0d0
      end do
*
*---  Arnoldi loop ...
*
      j1v = iv + n
      do 200 j = 1,m
         nmult = nmult + 1
         call matvec( wsp(j1v-n), wsp(j1v) )
         do i = 1,j
            hij = DDOT( n, wsp(iv+(i-1)*n),1, wsp(j1v),1 )
            call DAXPY( n, -hij, wsp(iv+(i-1)*n),1, wsp(j1v),1 )
            wsp(ih+(j-1)*mh+i-1) = hij
         end do
         hj1j = DNRM2( n, wsp(j1v),1 )
*---     if `happy breakdown' go straightforward at the end ...
         if ( hj1j.le.break_tol ) then
            print*,'happy breakdown: mbrkdwn =',j,' h =',hj1j
            k1 = 0
            ibrkflag = 1
            mbrkdwn = j
            tbrkdwn = t_now
            t_step = t_out-t_now
            goto 300
         endif
         wsp(ih+(j-1)*mh+j) = hj1j
         call DSCAL( n, 1.0d0/hj1j, wsp(j1v),1 )
         j1v = j1v + n
 200  continue
      nmult = nmult + 1
      call matvec( wsp(j1v-n), wsp(j1v) )
      avnorm = DNRM2( n, wsp(j1v),1 )
*
*---  set 1 for the 2-corrected scheme ...
*
 300  continue
      wsp(ih+m*mh+m+1) = 1.0d0
*
*---  loop while ireject<mxreject until the tolerance is reached ...
*
      ireject = 0
 401  continue

*
*---  compute w = beta*V*exp(t_step*H)*e1 ...
*
      nexph = nexph + 1
      mx = mbrkdwn + k1
      if ( ideg.ne.0 ) then
*---     irreducible rational Pade approximation ...
         call DGPADM( ideg, mx, sgn*t_step, wsp(ih),mh,
     .                wsp(ifree),lfree, iwsp, iexph, ns, iflag )
         iexph = ifree + iexph - 1
         nscale = nscale + ns
      else
*---     uniform rational Chebyshev approximation ...
         iexph = ifree
         do i = 1,mx
            wsp(iexph+i-1) = 0.0d0
         end do
         wsp(iexph) = 1.0d0
         call DNCHBV(mx,sgn*t_step,wsp(ih),mh,wsp(iexph),wsp(ifree+mx))
      endif

 402  continue
*
*---  error estimate ...
*
      if ( k1.eq.0 ) then
         err_loc = tol
      else
         p1 = ABS( wsp(iexph+m) )   * beta
         p2 = ABS( wsp(iexph+m+1) ) * beta * avnorm
         if ( p1.gt.10.0d0*p2 ) then
            err_loc = p2
            xm = 1.0d0/DBLE( m )
         elseif ( p1.gt.p2 ) then
            err_loc = (p1*p2)/(p1-p2)
            xm = 1.0d0/DBLE( m )
         else
            err_loc = p1
            xm = 1.0d0/DBLE( m-1 )
         endif
      endif
*
*---  reject the step-size if the error is not acceptable ...
*
      if ( (k1.ne.0) .and. (err_loc.gt.delta*t_step*tol) .and.
     .     (mxreject.eq.0 .or. ireject.lt.mxreject) ) then
         t_old = t_step
         t_step = gamma * t_step * (t_step*tol/err_loc)**xm
         p1 = 10.0d0**(NINT( LOG10( t_step )-SQR1 )-1)
         t_step = AINT( t_step/p1 + 0.55d0 ) * p1
         if ( itrace.ne.0 ) then
            print*,'t_step =',t_old
            print*,'err_loc =',err_loc
            print*,'err_required =',delta*t_old*tol
            print*,'stepsize rejected, stepping down to:',t_step
         endif
         ireject = ireject + 1
         nreject = nreject + 1
         if ( mxreject.ne.0 .and. ireject.gt.mxreject ) then
            print*,"Failure in DGEXPV: ---"
            print*,"The requested tolerance is too high."
            Print*,"Rerun with a smaller value."
            iflag = 2
            return
         endif
         goto 401
      endif
*
*---  now update w = beta*V*exp(t_step*H)*e1 and the hump ...
*
      mx = mbrkdwn + MAX( 0,k1-1 )
      call DGEMV( 'n', n,mx,beta,wsp(iv),n,wsp(iexph),1,0.0d0,w,1 )
      beta = DNRM2( n, w,1 )
      hump = MAX( hump, beta )
*
*---  suggested value for the next stepsize ...
*
      t_new = gamma * t_step * (t_step*tol/err_loc)**xm
      p1 = 10.0d0**(NINT( LOG10( t_new )-SQR1 )-1)
      t_new = AINT( t_new/p1 + 0.55d0 ) * p1

      err_loc = MAX( err_loc,rndoff )
*
*---  update the time covered ...
*
      t_now = t_now + t_step
*
*---  display and keep some information ...
*
      if ( itrace.ne.0 ) then
         print*,'integration',nstep,'---------------------------------'
         print*,'scale-square =',ns
         print*,'step_size =',t_step
         print*,'err_loc   =',err_loc
         print*,'next_step =',t_new
      endif

      step_min = MIN( step_min, t_step )
      step_max = MAX( step_max, t_step )
      s_error = s_error + err_loc
      x_error = MAX( x_error, err_loc )

      if ( mxstep.eq.0 .or. nstep.lt.mxstep ) goto 100
      iflag = 1

 500  continue

      iwsp(1) = nmult
      iwsp(2) = nexph
      iwsp(3) = nscale
      iwsp(4) = nstep
      iwsp(5) = nreject
      iwsp(6) = ibrkflag
      iwsp(7) = mbrkdwn

      wsp(1)  = step_min
      wsp(2)  = step_max
      wsp(3)  = 0.0d0
      wsp(4)  = 0.0d0
      wsp(5)  = x_error
      wsp(6)  = s_error
      wsp(7)  = tbrkdwn
      wsp(8)  = sgn*t_now
      wsp(9)  = hump/vnorm
      wsp(10) = beta/vnorm
      END
      subroutine dgpadm ( ideg, m, t, h, ldh, wsp, lwsp, ipiv, iexph,
     &  ns, iflag )

c*********************************************************************72
c
cc DGPADM computes exp(t*H) using an irreducible rational Pade approximation.
c
*-----Purpose----------------------------------------------------------|
*
*     Computes exp(t*H), the matrix exponential of a general matrix in
*     full, using the irreducible rational Pade approximation to the
*     exponential function exp(x) = r(x) = (+/-)( I + 2*(q(x)/p(x)) ),
*     combined with scaling-and-squaring.
*
*-----Arguments--------------------------------------------------------|
*
*     ideg      : (input) the degre of the diagonal Pade to be used.
*                 a value of 6 is generally satisfactory.
*
*     m         : (input) order of H.
*
*     H(ldh,m)  : (input) argument matrix.
*
*     t         : (input) time-scale (can be < 0).
*
*     wsp(lwsp) : (workspace/output) lwsp .ge. 4*m*m+ideg+1.
*
*     ipiv(m)   : (workspace)
*
*>>>> iexph     : (output) number such that wsp(iexph) points to exp(tH)
*                 i.e., exp(tH) is located at wsp(iexph ... iexph+m*m-1)
*                       ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
*                 NOTE: if the routine was called with wsp(iptr),
*                       then exp(tH) will start at wsp(iptr+iexph-1).
*
*     ns        : (output) number of scaling-squaring used.
*
*     iflag     : (output) exit flag.
*                      0 - no problem
*                     <0 - problem
*
*----------------------------------------------------------------------|
*     Roger B. Sidje (rbs@maths.uq.edu.au)
*     EXPOKIT: Software Package for Computing Matrix Exponentials.
*     ACM - Transactions On Mathematical Software, 24(1):130-156, 1998
*----------------------------------------------------------------------|
*
      implicit none

      integer ideg, m, ldh, lwsp, iexph, ns, iflag, ipiv(m)
      double precision t, H(ldh,m), wsp(lwsp)

      integer mm,i,j,k,ih2,ip,iq,iused,ifree,iodd,icoef,iput,iget
      double precision hnorm,scale,scale2,cp,cq

      intrinsic INT,ABS,DBLE,LOG,MAX

*---  check restrictions on input parameters ...
      mm = m*m
      iflag = 0
      if ( ldh.lt.m ) iflag = -1
      if ( lwsp.lt.4*mm+ideg+1 ) iflag = -2
      if ( iflag.ne.0 ) stop 'bad sizes (in input of DGPADM)'
*
*---  initialise pointers ...
*
      icoef = 1
      ih2 = icoef + (ideg+1)
      ip  = ih2 + mm
      iq  = ip + mm
      ifree = iq + mm
*
*---  scaling: seek ns such that ||t*H/2^ns|| < 1/2;
*     and set scale = t/2^ns ...
*
      do i = 1,m
         wsp(i) = 0.0d0
      end do
      do j = 1,m
         do i = 1,m
            wsp(i) = wsp(i) + ABS( H(i,j) )
         end do
      end do
      hnorm = 0.0d0
      do i = 1,m
         hnorm = MAX( hnorm,wsp(i) )
      end do
      hnorm = ABS( t*hnorm )
      if ( hnorm.eq.0.0d0 ) stop 'Error - null H in input of DGPADM.'
      ns = MAX( 0,INT(LOG(hnorm)/LOG(2.0d0))+2 )
      scale = t / DBLE(2**ns)
      scale2 = scale*scale
*
*---  compute Pade coefficients ...
*
      i = ideg+1
      j = 2*ideg+1
      wsp(icoef) = 1.0d0
      do k = 1,ideg
         wsp(icoef+k) = (wsp(icoef+k-1)*DBLE( i-k ))/DBLE( k*(j-k) )
      end do
*
*---  H2 = scale2*H*H ...
*
      call DGEMM( 'n','n',m,m,m,scale2,H,ldh,H,ldh,0.0d0,wsp(ih2),m )
*
*---  initialize p (numerator) and q (denominator) ...
*
      cp = wsp(icoef+ideg-1)
      cq = wsp(icoef+ideg)
      do j = 1,m
         do i = 1,m
            wsp(ip + (j-1)*m + i-1) = 0.0d0
            wsp(iq + (j-1)*m + i-1) = 0.0d0
         end do
         wsp(ip + (j-1)*(m+1)) = cp
         wsp(iq + (j-1)*(m+1)) = cq
      end do
*
*---  Apply Horner rule ...
*
      iodd = 1
      k = ideg - 1
 100  continue
      iused = iodd*iq + (1-iodd)*ip
      call DGEMM( 'n','n',m,m,m, 1.0d0,wsp(iused),m,
     .             wsp(ih2),m, 0.0d0,wsp(ifree),m )
      do j = 1,m
         wsp(ifree+(j-1)*(m+1)) = wsp(ifree+(j-1)*(m+1))+wsp(icoef+k-1)
      end do
      ip = (1-iodd)*ifree + iodd*ip
      iq = iodd*ifree + (1-iodd)*iq
      ifree = iused
      iodd = 1-iodd
      k = k-1
      if ( k.gt.0 )  goto 100
*
*---  Obtain (+/-)(I + 2*(p\q)) ...
*
      if ( iodd .eq. 1 ) then
         call DGEMM( 'n','n',m,m,m, scale,wsp(iq),m,
     .                H,ldh, 0.0d0,wsp(ifree),m )
         iq = ifree
      else
         call DGEMM( 'n','n',m,m,m, scale,wsp(ip),m,
     .                H,ldh, 0.0d0,wsp(ifree),m )
         ip = ifree
      endif
      call DAXPY( mm, -1.0d0,wsp(ip),1, wsp(iq),1 )
      call DGESV( m,m, wsp(iq),m, ipiv, wsp(ip),m, iflag )
      if ( iflag.ne.0 ) stop 'Problem in DGESV (within DGPADM)'
      call DSCAL( mm, 2.0d0, wsp(ip), 1 )
      do j = 1,m
         wsp(ip+(j-1)*(m+1)) = wsp(ip+(j-1)*(m+1)) + 1.0d0
      end do
      iput = ip
      if ( ns.eq.0 .and. iodd.eq.1 ) then
         call DSCAL( mm, -1.0d0, wsp(ip), 1 )
         goto 200
      endif
*
*--   squaring : exp(t*H) = (exp(t*H))^(2^ns) ...
*
      iodd = 1
      do k = 1,ns
         iget = iodd*ip + (1-iodd)*iq
         iput = (1-iodd)*ip + iodd*iq
         call DGEMM( 'n','n',m,m,m, 1.0d0,wsp(iget),m, wsp(iget),m,
     .                0.0d0,wsp(iput),m )
         iodd = 1-iodd
      end do
 200  continue
      iexph = iput
      END
      subroutine dgphiv( n, m, t, u, v, w, tol, anorm,
     .                   wsp,lwsp, iwsp,liwsp, matvec, itrace,iflag )

c*********************************************************************72
c
cc DGPHIV computes the solution of the ODE w'=A*w+u, w(0)=v.
c
*-----Purpose----------------------------------------------------------|
*
*---  DGPHIV computes w = exp(t*A)v + t*phi(tA)u which is the solution
*     of the nonhomogeneous linear ODE problem w' = Aw + u, w(0) = v.
*     phi(z) = (exp(z)-1)/z and A is a General matrix.
*
*     The method used is based on Krylov subspace projection
*     techniques and the matrix under consideration interacts only
*     via the external routine `matvec' performing the matrix-vector
*     product (matrix-free method).
*
*-----Arguments--------------------------------------------------------|
*
*     n      : (input) order of the principal matrix A.
*
*     m      : (input) maximum size for the Krylov basis.
*
*     t      : (input) time at wich the solution is needed (can be < 0).
*
*     u(n)   : (input) operand vector with respect to the phi function
*              (forcing term of the ODE problem).
*
*     v(n)   : (input) operand vector with respect to the exp function
*              (initial condition of the ODE problem).
*
*     w(n)   : (output) computed approximation of exp(t*A)v + t*phi(tA)u
*
*     tol    : (input/output) the requested accuracy tolerance on w.
*              If on input tol=0.0d0 or tol is too small (tol.le.eps)
*              the internal value sqrt(eps) is used, and tol is set to
*              sqrt(eps) on output (`eps' denotes the machine epsilon).
*              (`Happy breakdown' is assumed if h(j+1,j) .le. anorm*tol)
*
*     anorm  : (input) an approximation of some norm of A.
*
*   wsp(lwsp): (workspace) lwsp .ge. n*(m+1)+n+(m+3)^2+4*(m+3)^2+ideg+1
*                                   +---------+-------+---------------+
*              (actually, ideg=6)        V        H      wsp for PADE
*
* iwsp(liwsp): (workspace) liwsp .ge. m+3
*
*     matvec : external routine for matrix-vector multiplication.
*              synopsis: matvec( x, y )
*                        double precision x(*), y(*)
*              computes: y(1:n) <- A*x(1:n)
*                        where A is the principal matrix.
*
*     itrace : (input) running mode. 0=silent, 1=print step-by-step info
*
*     iflag  : (output) exit flag.
*              <0 - bad input arguments
*               0 - no problem
*               1 - maximum number of steps reached without convergence
*               2 - requested tolerance was too high
*
*-----Accounts on the computation--------------------------------------|
*     Upon exit, an interested user may retrieve accounts on the
*     computations. They are located in the workspace arrays wsp and
*     iwsp as indicated below:
*
*     location  mnemonic                 description
*     -----------------------------------------------------------------|
*     iwsp(1) = nmult, number of matrix-vector multiplications used
*     iwsp(2) = nexph, number of Hessenberg matrix exponential evaluated
*     iwsp(3) = nscale, number of repeated squaring involved in Pade
*     iwsp(4) = nstep, number of integration steps used up to completion
*     iwsp(5) = nreject, number of rejected step-sizes
*     iwsp(6) = ibrkflag, set to 1 if `happy breakdown' and 0 otherwise
*     iwsp(7) = mbrkdwn, if `happy brkdown', basis-size when it occured
*     -----------------------------------------------------------------|
*     wsp(1)  = step_min, minimum step-size used during integration
*     wsp(2)  = step_max, maximum step-size used during integration
*     wsp(3)  = dummy
*     wsp(4)  = dummy
*     wsp(5)  = x_error, maximum among all local truncation errors
*     wsp(6)  = s_error, global sum of local truncation errors
*     wsp(7)  = tbrkdwn, if `happy breakdown', time when it occured
*     wsp(8)  = t_now, integration domain successfully covered
*
*----------------------------------------------------------------------|
*-----The following parameters may also be adjusted herein-------------|
*
      implicit none

      integer n, m, lwsp, liwsp, itrace, iflag, iwsp(liwsp)
      double precision t, tol, anorm, u(n), v(n), w(n), wsp(lwsp)
      external matvec

      integer mxstep, mxreject, ideg
      double precision delta, gamma
      parameter( mxstep   = 1000,
     .           mxreject = 0,
     .           ideg     = 6,
     .           delta    = 1.2d0,
     .           gamma    = 0.9d0 )

*     mxstep  : maximum allowable number of integration steps.
*               The value 0 means an infinite number of steps.
*
*     mxreject: maximum allowable number of rejections at each step.
*               The value 0 means an infinite number of rejections.
*
*     ideg    : the Pade approximation of type (ideg,ideg) is used as
*               an approximation to exp(H).
*
*     delta   : local truncation error `safety factor'
*
*     gamma   : stepsize `shrinking factor'
*
*----------------------------------------------------------------------|
*     Roger B. Sidje (rbs@maths.uq.edu.au)
*     EXPOKIT: Software Package for Computing Matrix Exponentials.
*     ACM - Transactions On Mathematical Software, 24(1):130-156, 1998
*----------------------------------------------------------------------|
*
      integer i, j, k1, mh, mx, iv, ih, j1v, ns, ifree, lfree, iexph,
     .        ireject,ibrkflag,mbrkdwn, nmult, nreject, nexph, nscale,
     .        nstep, iphih
      double precision sgn, t_out, tbrkdwn, step_min,step_max, err_loc,
     .                 s_error, x_error, t_now, t_new, t_step, t_old,
     .                 xm, beta, break_tol, p1, p2, p3, eps, rndoff,
     .                 avnorm, hj1j, hij, SQR1

      intrinsic AINT,ABS,DBLE,LOG10,MAX,MIN,NINT,SIGN,SQRT
      double precision DDOT, DNRM2

*---  check restrictions on input parameters ...
      iflag = 0
      if ( lwsp.lt.n*(m+3)+5*(m+3)**2+ideg+1 ) iflag = -1
      if ( liwsp.lt.m+3 ) iflag = -2
      if ( m.ge.n .or. m.le.0 ) iflag = -3
      if ( iflag.ne.0 ) stop 'bad sizes (in input of DGPHIV)'
*
*---  initialisations ...
*
      k1 = 3
      mh = m + 3
      iv = 1
      ih = iv + n*(m+1) + n
      ifree = ih + mh*mh
      lfree = lwsp - ifree + 1

      ibrkflag = 0
      mbrkdwn  = m
      nmult    = 0
      nreject  = 0
      nexph    = 0
      nscale   = 0

      t_out    = ABS( t )
      tbrkdwn  = 0.0d0
      step_min = t_out
      step_max = 0.0d0
      nstep    = 0
      s_error  = 0.0d0
      x_error  = 0.0d0
      t_now    = 0.0d0
      t_new    = 0.0d0

      p1 = 4.0d0/3.0d0
 1    p2 = p1 - 1.0d0
      p3 = p2 + p2 + p2
      eps = ABS( p3-1.0d0 )
      if ( eps.eq.0.0d0 ) go to 1
      if ( tol.le.eps ) tol = SQRT( eps )
      rndoff = eps*anorm

      break_tol = 1.0d-7
*>>>  break_tol = tol
*>>>  break_tol = anorm*tol

*
*---  step-by-step integration ...
*
      sgn = SIGN( 1.0d0,t )
      SQR1 = SQRT( 0.1d0 )
      call DCOPY( n, v,1, w,1 )

 100  if ( t_now.ge.t_out ) goto 500

      nmult =  nmult + 1
      call matvec( w, wsp(iv) )
      call DAXPY( n, 1.0d0, u,1, wsp(iv),1 )
      beta = DNRM2( n, wsp(iv),1 )
      if ( beta.eq.0.0d0 ) goto 500
      call DSCAL( n, 1.0d0/beta, wsp(iv),1 )
      do i = 1,mh*mh
         wsp(ih+i-1) = 0.0d0
      end do

      if ( nstep.eq.0 ) then
*---     obtain the very first stepsize ...
         xm = 1.0d0/DBLE( m )
         p1 = tol*(((m+1)/2.72D0)**(m+1))*SQRT(2.0D0*3.14D0*(m+1))
         t_new = (1.0d0/anorm)*(p1/(4.0d0*beta*anorm))**xm
         p1 = 10.0d0**(NINT( LOG10( t_new )-SQR1 )-1)
         t_new = AINT( t_new/p1 + 0.55d0 ) * p1
      endif
      nstep = nstep + 1
      t_step = MIN( t_out-t_now, t_new )
*
*---  Arnoldi loop ...
*
      j1v = iv + n
      do 200 j = 1,m
         nmult = nmult + 1
         call matvec( wsp(j1v-n), wsp(j1v) )
         do i = 1,j
            hij = DDOT( n, wsp(iv+(i-1)*n),1, wsp(j1v),1 )
            call DAXPY( n, -hij, wsp(iv+(i-1)*n),1, wsp(j1v),1 )
            wsp(ih+(j-1)*mh+i-1) = hij
         end do
         hj1j = DNRM2( n, wsp(j1v),1 )
*---     if `happy breakdown' go straightforward at the end ...
         if ( hj1j.le.break_tol ) then
            print*,'happy breakdown: mbrkdwn =',j,' h =',hj1j
            k1 = 0
            ibrkflag = 1
            mbrkdwn = j
            tbrkdwn = t_now
            t_step = t_out-t_now
            goto 300
         endif
         wsp(ih+(j-1)*mh+j) = hj1j
         call DSCAL( n, 1.0d0/hj1j, wsp(j1v),1 )
         j1v = j1v + n
 200  continue
      nmult = nmult + 1
      call matvec( wsp(j1v-n), wsp(j1v) )
      avnorm = DNRM2( n, wsp(j1v),1 )
*
*---  set 1's for the 3-extended scheme ...
*
 300  continue
      wsp(ih+mh*mbrkdwn) = 1.0d0
      wsp(ih+(m-1)*mh+m) = 0.0d0
      do i = 1,k1-1
         wsp(ih+(m+i)*mh+m+i-1) = 1.0d0
      end do
*
*---  loop while ireject<mxreject until the tolerance is reached ...
*
      ireject = 0
 401  continue
*
*---  compute w = beta*t_step*V*phi(t_step*H)*e1 + w
*
      nexph = nexph + 1
*---  irreducible rational Pade approximation ...
      mx = mbrkdwn + MAX( 1,k1 )
      call DGPADM( ideg, mx, sgn*t_step, wsp(ih),mh,
     .             wsp(ifree),lfree, iwsp, iexph, ns, iflag )
      iexph = ifree + iexph - 1
      iphih = iexph + mbrkdwn*mx
      nscale = nscale + ns
      wsp(iphih+mbrkdwn)   = hj1j*wsp(iphih+mx+mbrkdwn-1)
      wsp(iphih+mbrkdwn+1) = hj1j*wsp(iphih+2*mx+mbrkdwn-1)

 402  continue
*---  error estimate ...
      if ( k1.eq.0 ) then
         err_loc = tol
      else
         p1 = ABS( wsp(iphih+m) )   * beta
         p2 = ABS( wsp(iphih+m+1) ) * beta * avnorm
         if ( p1.gt.10.0d0*p2 ) then
            err_loc = p2
            xm = 1.0d0/DBLE( m+1 )
         elseif ( p1.gt.p2 ) then
            err_loc = (p1*p2)/(p1-p2)
            xm = 1.0d0/DBLE( m+1 )
         else
            err_loc = p1
            xm = 1.0d0/DBLE( m )
         endif
      endif

*---  reject the step-size if the error is not acceptable ...
      if ( (k1.ne.0) .and. (err_loc.gt.delta*t_step*tol) .and.
     .     (mxreject.eq.0 .or. ireject.lt.mxreject) ) then
         t_old = t_step
         t_step = gamma * t_step * (t_step*tol/err_loc)**xm
         p1 = 10.0d0**(NINT( LOG10( t_step )-SQR1 )-1)
         t_step = AINT( t_step/p1 + 0.55d0 ) * p1
         if ( itrace.ne.0 ) then
            print*,'t_step =',t_old
            print*,'err_loc =',err_loc
            print*,'err_required =',delta*t_old*tol
            print*,'stepsize rejected, stepping down to:',t_step
         endif
         ireject = ireject + 1
         nreject = nreject + 1
         if ( mxreject.ne.0 .and. ireject.gt.mxreject ) then
            print*,"Failure in DGPHIV: ---"
            print*,"The requested tolerance is too high."
            Print*,"Rerun with a smaller value."
            iflag = 2
            return
         endif
         goto 401
      endif
*
      mx = mbrkdwn + MAX( 0,k1-2 )
      call DGEMV( 'n', n,mx,beta,wsp(iv),n,wsp(iphih),1,1.0d0,w,1 )
*
*---  suggested value for the next stepsize ...
*
      t_new = gamma * t_step * (t_step*tol/err_loc)**xm
      p1 = 10.0d0**(NINT( LOG10( t_new )-SQR1 )-1)
      t_new = AINT( t_new/p1 + 0.55d0 ) * p1

      err_loc = MAX( err_loc,rndoff )
*
*---  update the time covered ...
*
      t_now = t_now + t_step
*
*---  display and keep some information ...
*
      if ( itrace.ne.0 ) then
         print*,'integration',nstep,'---------------------------------'
         print*,'scale-square =',ns
         print*,'step_size =',t_step
         print*,'err_loc   =',err_loc
         print*,'next_step =',t_new
      endif

      step_min = MIN( step_min, t_step )
      step_max = MAX( step_max, t_step )
      s_error = s_error + err_loc
      x_error = MAX( x_error, err_loc )

      if ( mxstep.eq.0 .or. nstep.lt.mxstep ) goto 100
      iflag = 1

 500  continue

      iwsp(1) = nmult
      iwsp(2) = nexph
      iwsp(3) = nscale
      iwsp(4) = nstep
      iwsp(5) = nreject
      iwsp(6) = ibrkflag
      iwsp(7) = mbrkdwn

      wsp(1)  = step_min
      wsp(2)  = step_max
      wsp(3)  = 0.0d0
      wsp(4)  = 0.0d0
      wsp(5)  = x_error
      wsp(6)  = s_error
      wsp(7)  = tbrkdwn
      wsp(8)  = sgn*t_now
      END
      function dlaran ( iseed )

c*********************************************************************72
c
cc DLARAN returns a random real number uniformly selected from (0,1).
*
*  -- LAPACK auxiliary routine (version 2.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     February 29, 1992
*
*  Arguments
*  =========
*
*  ISEED   (input/output) INTEGER array, dimension (4)
*          On entry, the seed of the random number generator; the array
*          elements must be between 0 and 4095, and ISEED(4) must be
*          odd.
*          On exit, the seed is updated.
*
*  Further Details
*  ===============
*
*  This routine uses a multiplicative congruential method with modulus
*  2**48 and multiplier 33952834046453 (see G.S.Fishman,
*  'Multiplicative congruential random number generators with modulus
*  2**b: an exhaustive analysis for b = 32 and a partial analysis for
*  b = 48', Math. Comp. 189, pp 331-344, 1990).
*
*  48-bit integers are stored in 4 integer array elements with 12 bits
*  per element. Hence the routine is portable across machines with
*  integers of 32 bits or more.
*
      implicit none

      double precision dlaran
      integer            iseed( 4 )

      integer            m1, m2, m3, m4
      parameter          ( m1 = 494, m2 = 322, m3 = 2508, m4 = 2549 )
      double precision   one
      parameter          ( one = 1.0d+0 )
      integer            ipw2
      double precision   r
      parameter          ( ipw2 = 4096, r = one / ipw2 )
*     ..
*     .. Local Scalars ..
      integer            it1, it2, it3, it4
*     ..
*     .. Intrinsic Functions ..
      intrinsic          dble, mod
*     ..
*     .. Executable Statements ..
*
*     multiply the seed by the multiplier modulo 2**48
*
      it4 = iseed( 4 )*m4
      it3 = it4 / ipw2
      it4 = it4 - ipw2*it3
      it3 = it3 + iseed( 3 )*m4 + iseed( 4 )*m3
      it2 = it3 / ipw2
      it3 = it3 - ipw2*it2
      it2 = it2 + iseed( 2 )*m4 + iseed( 3 )*m3 + iseed( 4 )*m2
      it1 = it2 / ipw2
      it2 = it2 - ipw2*it1
      it1 = it1 + iseed( 1 )*m4 + iseed( 2 )*m3 + iseed( 3 )*m2 +
     $      iseed( 4 )*m1
      it1 = mod( it1, ipw2 )
*
*     return updated seed
*
      iseed( 1 ) = it1
      iseed( 2 ) = it2
      iseed( 3 ) = it3
      iseed( 4 ) = it4
*
*     convert 48-bit integer to a real number in the interval (0,1)
*
      dlaran = r*( dble( it1 )+r*( dble( it2 )+r*( dble( it3 )+r*
     $         ( dble( it4 ) ) ) ) )
      return
*
*     end of dlaran
*
      end
      subroutine dmexpv( n, m, t, v, w, tol, anorm,
     .                   wsp,lwsp, iwsp,liwsp, matvec, itrace,iflag )

c*********************************************************************72
c
cc DMEXPV computes w = exp(t*A)*v, customized for Markov chains.
c
      implicit none

      integer n, m, lwsp, liwsp, itrace, iflag, iwsp(liwsp)
      double precision t, tol, anorm, v(n), w(n), wsp(lwsp)
      external matvec

*-----Purpose----------------------------------------------------------|
*
*---  DMEXPV computes w = exp(t*A)*v - Customised for MARKOV CHAINS.
*
*     It does not compute the matrix exponential in isolation but
*     instead, it computes directly the action of the exponential
*     operator on the operand vector. This way of doing so allows
*     for addressing large sparse problems.
*
*     The method used is based on Krylov subspace projection
*     techniques and the matrix under consideration interacts only
*     via the external routine `matvec' performing the matrix-vector
*     product (matrix-free method).
*
*     This is a customised version for Markov Chains. This means that a
*     check is done within this code to ensure that the resulting vector
*     w is a probability vector, i.e., w must have all its components
*     in [0,1], with sum equal to 1. This check is done at some expense
*     and the user may try DGEXPV which is cheaper since it ignores
*     probability constraints.
*
*     IMPORTANT: The check assumes that the transition rate matrix Q
*                satisfies Qe = 0, where e=(1,...,1)'. Don't use DMEXPV
*                if this condition does not hold. Use DGEXPV instead.
*                DMEXPV/DGEXPV require the matrix-vector product
*                y = A*x = Q'*x, i.e, the TRANSPOSE of Q times a vector.
*                Failure to remember this leads to wrong results.
*
*-----Arguments--------------------------------------------------------|
*
*     n      : (input) order of the principal matrix A.
*
*     m      : (input) maximum size for the Krylov basis.
*
*     t      : (input) time at wich the solution is needed (can be < 0).
*
*     v(n)   : (input) given operand vector.
*
*     w(n)   : (output) computed approximation of exp(t*A)*v.
*
*     tol    : (input/output) the requested acurracy tolerance on w.
*              If on input tol=0.0d0 or tol is too small (tol.le.eps)
*              the internal value sqrt(eps) is used, and tol is set to
*              sqrt(eps) on output (`eps' denotes the machine epsilon).
*              (`Happy breakdown' is assumed if h(j+1,j) .le. anorm*tol)
*
*     anorm  : (input) an approximation of some norm of A.
*
*   wsp(lwsp): (workspace) lwsp .ge. n*(m+1)+n+(m+2)^2+4*(m+2)^2+ideg+1
*                                   +---------+-------+---------------+
*              (actually, ideg=6)        V        H     wsp for PADE
*
* iwsp(liwsp): (workspace) liwsp .ge. m+2
*
*     matvec : external routine for matrix-vector multiplication.
*              synopsis: matvec( x, y )
*                        double precision x(*), y(*)
*              computes: y(1:n) <- A*x(1:n)
*                        where A is the principal matrix.
*
*              IMPORTANT: DMEXPV requires the product y = Ax = Q'x, i.e.
*              the TRANSPOSE of the transition rate matrix.
*
*     itrace : (input) running mode. 0=silent, 1=print step-by-step info
*
*     iflag  : (output) exit flag.
*              <0 - bad input arguments
*               0 - no problem
*               1 - maximum number of steps reached without convergence
*               2 - requested tolerance was too high
*
*-----Accounts on the computation--------------------------------------|
*     Upon exit, an interested user may retrieve accounts on the
*     computations. They are located in the workspace arrays wsp and
*     iwsp as indicated below:
*
*     location  mnemonic                 description
*     -----------------------------------------------------------------|
*     iwsp(1) = nmult, number of matrix-vector multiplications used
*     iwsp(2) = nexph, number of Hessenberg matrix exponential evaluated
*     iwsp(3) = nscale, number of repeated squaring involved in Pade
*     iwsp(4) = nstep, number of integration steps used up to completion
*     iwsp(5) = nreject, number of rejected step-sizes
*     iwsp(6) = ibrkflag, set to 1 if `happy breakdown' and 0 otherwise
*     iwsp(7) = mbrkdwn, if `happy brkdown', basis-size when it occured
*     -----------------------------------------------------------------|
*     wsp(1)  = step_min, minimum step-size used during integration
*     wsp(2)  = step_max, maximum step-size used during integration
*     wsp(3)  = x_round, maximum among all roundoff errors (lower bound)
*     wsp(4)  = s_round, sum of roundoff errors (lower bound)
*     wsp(5)  = x_error, maximum among all local truncation errors
*     wsp(6)  = s_error, global sum of local truncation errors
*     wsp(7)  = tbrkdwn, if `happy breakdown', time when it occured
*     wsp(8)  = t_now, integration domain successfully covered
*     wsp(9)  = hump, i.e., max||exp(sA)||, s in [0,t] (or [t,0] if t<0)
*     wsp(10) = ||w||/||v||, scaled norm of the solution w.
*     -----------------------------------------------------------------|
*     The `hump' is a measure of the conditioning of the problem. The
*     matrix exponential is well-conditioned if hump = 1, whereas it is
*     poorly-conditioned if hump >> 1. However the solution can still be
*     relatively fairly accurate even when the hump is large (the hump
*     is an upper bound), especially when the hump and the scaled norm
*     of w [this is also computed and returned in wsp(10)] are of the
*     same order of magnitude (further details in reference below).
*     Markov chains are usually well-conditioned problems.
*
*----------------------------------------------------------------------|
*-----The following parameters may also be adjusted herein-------------|
*
      integer mxstep, mxreject, ideg
      double precision delta, gamma
      parameter( mxstep   = 500,
     .           mxreject = 0,
     .           ideg     = 6,
     .           delta    = 1.2d0,
     .           gamma    = 0.9d0 )

*     mxstep  : maximum allowable number of integration steps.
*               The value 0 means an infinite number of steps.
*
*     mxreject: maximum allowable number of rejections at each step.
*               The value 0 means an infinite number of rejections.
*
*     ideg    : the Pade approximation of type (ideg,ideg) is used as
*               an approximation to exp(H). The value 0 switches to the
*               uniform rational Chebyshev approximation of type (14,14)
*
*     delta   : local truncation error `safety factor'
*
*     gamma   : stepsize `shrinking factor'
*
*----------------------------------------------------------------------|
*     Roger B. Sidje (rbs@maths.uq.edu.au)
*     EXPOKIT: Software Package for Computing Matrix Exponentials.
*     ACM - Transactions On Mathematical Software, 24(1):130-156, 1998
*----------------------------------------------------------------------|
*
      integer i, j, k1, mh, mx, iv, ih, j1v, ns, ifree, lfree, iexph,
     .        ireject,ibrkflag,mbrkdwn, nmult, nreject, nexph, nscale,
     .        nstep
      double precision sgn, t_out, tbrkdwn, step_min,step_max, err_loc,
     .                 s_error, x_error, t_now, t_new, t_step, t_old,
     .                 xm, beta, break_tol, p1, p2, p3, eps, rndoff,
     .                 vnorm, avnorm, hj1j, hij, hump, SQR1,
     .                 roundoff, s_round, x_round

      intrinsic AINT,ABS,DBLE,LOG10,MAX,MIN,NINT,SIGN,SQRT
      double precision DDOT, DNRM2, DASUM

*---  check restrictions on input parameters ...
      iflag = 0
      if ( lwsp.lt.n*(m+2)+5*(m+2)**2+ideg+1 ) iflag = -1
      if ( liwsp.lt.m+2 ) iflag = -2
      if ( m.ge.n .or. m.le.0 ) iflag = -3
      if ( iflag.ne.0 ) stop 'bad sizes (in input of DMEXPV)'
*
*---  initialisations ...
*
      k1 = 2
      mh = m + 2
      iv = 1
      ih = iv + n*(m+1) + n
      ifree = ih + mh*mh
      lfree = lwsp - ifree + 1

      ibrkflag = 0
      mbrkdwn  = m
      nmult    = 0
      nreject  = 0
      nexph    = 0
      nscale   = 0

      sgn      = SIGN( 1.0d0,t )
      t_out    = ABS( t )
      tbrkdwn  = 0.0d0
      step_min = t_out
      step_max = 0.0d0
      nstep    = 0
      s_error  = 0.0d0
      s_round  = 0.0d0
      x_error  = 0.0d0
      x_round  = 0.0d0
      t_now    = 0.0d0
      t_new    = 0.0d0

      p1 = 4.0d0/3.0d0
 1    p2 = p1 - 1.0d0
      p3 = p2 + p2 + p2
      eps = ABS( p3-1.0d0 )
      if ( eps.eq.0.0d0 ) go to 1
      if ( tol.le.eps ) tol = SQRT( eps )
      rndoff = eps*anorm

      break_tol = 1.0d-7
*>>>  break_tol = tol
*>>>  break_tol = anorm*tol

      call DCOPY( n, v,1, w,1 )
      beta = DNRM2( n, w,1 )
      vnorm = beta
      hump = beta
*
*---  obtain the very first stepsize ...
*
      SQR1 = SQRT( 0.1d0 )
      xm = 1.0d0/DBLE( m )
      p1 = tol*(((m+1)/2.72D0)**(m+1))*SQRT(2.0D0*3.14D0*(m+1))
      t_new = (1.0d0/anorm)*(p1/(4.0d0*beta*anorm))**xm
      p1 = 10.0d0**(NINT( LOG10( t_new )-SQR1 )-1)
      t_new = AINT( t_new/p1 + 0.55d0 ) * p1
*
*---  step-by-step integration ...
*
 100  if ( t_now.ge.t_out ) goto 500

      nstep = nstep + 1
      t_step = MIN( t_out-t_now, t_new )

      p1 = 1.0d0/beta
      do i = 1,n
         wsp(iv + i-1) = p1*w(i)
      end do
      do i = 1,mh*mh
         wsp(ih+i-1) = 0.0d0
      end do
*
*---  Arnoldi loop ...
*
      j1v = iv + n
      do 200 j = 1,m
         nmult = nmult + 1
         call matvec( wsp(j1v-n), wsp(j1v) )
         do i = 1,j
            hij = DDOT( n, wsp(iv+(i-1)*n),1, wsp(j1v),1 )
            call DAXPY( n, -hij, wsp(iv+(i-1)*n),1, wsp(j1v),1 )
            wsp(ih+(j-1)*mh+i-1) = hij
         end do
         hj1j = DNRM2( n, wsp(j1v),1 )
*---     if `happy breakdown' go straightforward at the end ...
         if ( hj1j.le.break_tol ) then
            print*,'happy breakdown: mbrkdwn =',j,' h =',hj1j
            k1 = 0
            ibrkflag = 1
            mbrkdwn = j
            tbrkdwn = t_now
            t_step = t_out-t_now
            goto 300
         endif
         wsp(ih+(j-1)*mh+j) = hj1j
         call DSCAL( n, 1.0d0/hj1j, wsp(j1v),1 )
         j1v = j1v + n
 200  continue
      nmult = nmult + 1
      call matvec( wsp(j1v-n), wsp(j1v) )
      avnorm = DNRM2( n, wsp(j1v),1 )
*
*---  set 1 for the 2-corrected scheme ...
*
 300  continue
      wsp(ih+m*mh+m+1) = 1.0d0
*
*---  loop while ireject<mxreject until the tolerance is reached ...
*
      ireject = 0
 401  continue

*
*---  compute w = beta*V*exp(t_step*H)*e1 ..
*
      nexph = nexph + 1
      mx = mbrkdwn + k1
      if ( ideg.ne.0 ) then
*---     irreducible rational Pade approximation ...
         call DGPADM( ideg, mx, sgn*t_step, wsp(ih),mh,
     .                wsp(ifree),lfree, iwsp, iexph, ns, iflag )
         iexph = ifree + iexph - 1
         nscale = nscale + ns
      else
*---     uniform rational Chebyshev approximation ...
         iexph = ifree
         do i = 1,mx
            wsp(iexph+i-1) = 0.0d0
         end do
         wsp(iexph) = 1.0d0
         call DNCHBV(mx,sgn*t_step,wsp(ih),mh,wsp(iexph),wsp(ifree+mx))
      endif

 402  continue
*
*---  error estimate ...
*
      if ( k1.eq.0 ) then
         err_loc = tol
      else
         p1 = ABS( wsp(iexph+m) )   * beta
         p2 = ABS( wsp(iexph+m+1) ) * beta * avnorm
         if ( p1.gt.10.0d0*p2 ) then
            err_loc = p2
            xm = 1.0d0/DBLE( m )
         elseif ( p1.gt.p2 ) then
            err_loc = (p1*p2)/(p1-p2)
            xm = 1.0d0/DBLE( m )
         else
            err_loc = p1
            xm = 1.0d0/DBLE( m-1 )
         endif
      endif
*
*---  reject the step-size if the error is not acceptable ...
*
      if ( (k1.ne.0) .and. (err_loc.gt.delta*t_step*tol) .and.
     .     (mxreject.eq.0 .or. ireject.lt.mxreject) ) then
         t_old = t_step
         t_step = gamma * t_step * (t_step*tol/err_loc)**xm
         p1 = 10.0d0**(NINT( LOG10( t_step )-SQR1 )-1)
         t_step = AINT( t_step/p1 + 0.55d0 ) * p1
         if ( itrace.ne.0 ) then
            print*,'t_step =',t_old
            print*,'err_loc =',err_loc
            print*,'err_required =',delta*t_old*tol
            print*,'stepsize rejected, stepping down to:',t_step
         endif
         ireject = ireject + 1
         nreject = nreject + 1
         if ( mxreject.ne.0 .and. ireject.gt.mxreject ) then
            print*,"Failure in DMEXPV: ---"
            print*,"The requested tolerance is too high."
            Print*,"Rerun with a smaller value."
            iflag = 2
            return
         endif
         goto 401
      endif
*
*---  now update w = beta*V*exp(t_step*H)*e1 and the hump ...
*
      mx = mbrkdwn + MAX( 0,k1-1 )
      call DGEMV( 'n', n,mx,beta,wsp(iv),n,wsp(iexph),1,0.0d0,w,1 )
      beta = DNRM2( n, w,1 )
      hump = MAX( hump, beta )
*
*---  Markov model constraints ...
*
      j = 0
      do i = 1,n
         if ( w(i).lt.0.0d0 ) then
            w(i) = 0.0d0
            j = j + 1
         endif
      end do
      p1 = DASUM( n, w,1 )
      if ( j.gt.0 ) call DSCAL( n, 1.0d0/p1, w,1 )
      roundoff = DABS( 1.0d0-p1 ) / DBLE( n )
*
*---  suggested value for the next stepsize ...
*
      t_new = gamma * t_step * (t_step*tol/err_loc)**xm
      p1 = 10.0d0**(NINT( LOG10( t_new )-SQR1 )-1)
      t_new = AINT( t_new/p1 + 0.55d0 ) * p1

      err_loc = MAX( err_loc, roundoff )
      err_loc = MAX( err_loc, rndoff )
*
*---  update the time covered ...
*
      t_now = t_now + t_step
*
*---  display and keep some information ...
*
      if ( itrace.ne.0 ) then
         print*,'integration',nstep,'---------------------------------'
         print*,'scale-square =',ns
         print*,'wnorm     =',beta
         print*,'step_size =',t_step
         print*,'err_loc   =',err_loc
         print*,'roundoff  =',roundoff
         print*,'next_step =',t_new
      endif

      step_min = MIN( step_min, t_step )
      step_max = MAX( step_max, t_step )
      s_error = s_error + err_loc
      s_round = s_round + roundoff
      x_error = MAX( x_error, err_loc )
      x_round = MAX( x_round, roundoff )

      if ( mxstep.eq.0 .or. nstep.lt.mxstep ) goto 100
      iflag = 1

 500  continue

      iwsp(1) = nmult
      iwsp(2) = nexph
      iwsp(3) = nscale
      iwsp(4) = nstep
      iwsp(5) = nreject
      iwsp(6) = ibrkflag
      iwsp(7) = mbrkdwn

      wsp(1)  = step_min
      wsp(2)  = step_max
      wsp(3)  = x_round
      wsp(4)  = s_round
      wsp(5)  = x_error
      wsp(6)  = s_error
      wsp(7)  = tbrkdwn
      wsp(8)  = sgn*t_now
      wsp(9)  = hump/vnorm
      wsp(10) = beta/vnorm
      END
      subroutine dnchbv( m, t, h,ldh, y, wsp )

c*********************************************************************72
c
cc DNCHBV computes y = exp(t*H)*y using a partial fraction expansion.
c
      implicit none

      integer          m, ldh
      double precision t, H(ldh,m), y(m)
      complex*16       wsp(m*(m+2))

*-----Purpose----------------------------------------------------------|
*
*---  DNCHBV computes y = exp(t*H)*y using the partial fraction
*     expansion of the uniform rational Chebyshev approximation
*     to exp(-x) of type (14,14). H is assumed to be upper-Hessenberg.
*     About 14-digit accuracy is expected if the matrix H is negative
*     definite. The algorithm may behave poorly otherwise.
*
*-----Arguments--------------------------------------------------------|
*
*     m       : (input) order of the Hessenberg matrix H
*
*     t       : (input) time-scaling factor (can be < 0).
*
*     H(ldh,m): (input) upper Hessenberg matrix.
*
*     y(m)    : (input/output) on input the operand vector,
*               on output the resulting vector exp(t*H)*y.
*
*     wsp     : (workspace). Observe that a double precision vector of
*               length 2*m*(m+2) can be used as well when calling this
*               routine (thus avoiding an idle complex array elsewhere)
*
*----------------------------------------------------------------------|
*     Roger B. Sidje (rbs@maths.uq.edu.au)
*     EXPOKIT: Software Package for Computing Matrix Exponentials.
*     ACM - Transactions On Mathematical Software, 24(1):130-156, 1998
*----------------------------------------------------------------------|
*
      complex*16 ZERO
      integer ndeg, i, j, k, ip, ih, iy, iz
      parameter ( ndeg=7, ZERO=(0.0d0,0.0d0) )
      double precision alpha0
      complex*16 alpha(ndeg), theta(ndeg), tmpc

      intrinsic ABS,DBLE,MIN

*---  Pointers ...

      ih = 1
      iy = ih + m*m
      iz = iy + m

*---  Coefficients and poles of the partial fraction expansion...

      alpha0  =  0.183216998528140087D-11
      alpha(1)=( 0.557503973136501826D+02,-0.204295038779771857D+03)
      alpha(2)=(-0.938666838877006739D+02, 0.912874896775456363D+02)
      alpha(3)=( 0.469965415550370835D+02,-0.116167609985818103D+02)
      alpha(4)=(-0.961424200626061065D+01,-0.264195613880262669D+01)
      alpha(5)=( 0.752722063978321642D+00, 0.670367365566377770D+00)
      alpha(6)=(-0.188781253158648576D-01,-0.343696176445802414D-01)
      alpha(7)=( 0.143086431411801849D-03, 0.287221133228814096D-03)

      theta(1)=(-0.562314417475317895D+01, 0.119406921611247440D+01)
      theta(2)=(-0.508934679728216110D+01, 0.358882439228376881D+01)
      theta(3)=(-0.399337136365302569D+01, 0.600483209099604664D+01)
      theta(4)=(-0.226978543095856366D+01, 0.846173881758693369D+01)
      theta(5)=( 0.208756929753827868D+00, 0.109912615662209418D+02)
      theta(6)=( 0.370327340957595652D+01, 0.136563731924991884D+02)
      theta(7)=( 0.889777151877331107D+01, 0.166309842834712071D+02)
*
*---  Accumulation of the contribution of each pole ...
*
      do j = 1,m
         wsp(iz+j-1) = y(j)
         y(j) = y(j)*alpha0
      end do
      do ip = 1,ndeg
*---     Solve each fraction using Gaussian elimination with pivoting...
         do j = 1,m
            wsp(iy+j-1) = wsp(iz+j-1)
            do i = 1,MIN( j+1,m )
               wsp(ih+(j-1)*m+i-1) = -t*H(i,j)
            end do
            wsp(ih+(j-1)*m+j-1) = wsp(ih+(j-1)*m+j-1)-theta(ip)
            do k = i,m
               wsp(ih+(j-1)*m+k-1) = ZERO
            end do
         end do
         do i = 1,m-1
*---        Get pivot and exchange rows ...
            if (ABS(wsp(ih+(i-1)*m+i-1)).lt.ABS(wsp(ih+(i-1)*m+i))) then
               call ZSWAP( m-i+1, wsp(ih+(i-1)*m+i-1),m,
     .                     wsp(ih+(i-1)*m+i),m )
               call ZSWAP( 1, wsp(iy+i-1),1, wsp(iy+i),1 )
            endif
*---        Forward eliminiation ...
            tmpc = wsp(ih+(i-1)*m+i) / wsp(ih+(i-1)*m+i-1)
            call ZAXPY( m-i, -tmpc, wsp(ih+i*m+i-1),m, wsp(ih+i*m+i),m )
            wsp(iy+i) = wsp(iy+i) - tmpc*wsp(iy+i-1)
         end do
*---     Backward substitution ...
         do i = m,1,-1
            tmpc = wsp(iy+i-1)
            do j = i+1,m
               tmpc = tmpc - wsp(ih+(j-1)*m+i-1)*wsp(iy+j-1)
            end do
            wsp(iy+i-1) = tmpc / wsp(ih+(i-1)*m+i-1)
         end do
*---     Accumulate the partial result in y ...
         do j = 1,m
            y(j) = y(j) + DBLE( alpha(ip)*wsp(iy+j-1) )
         end do
      end do
      END
      function dnrm2 ( n, dx, incx)

c*********************************************************************72
c
cc DNRM2 computes the Euclidean norm of a vector.
c
      implicit none

      double precision dnrm2
      integer i, incx, ix, j, n, next
      double precision   dx(1), cutlo, cuthi, hitest, sum, xmax,zero,one
      data   zero, one /0.0d0, 1.0d0/
c
c     euclidean norm of the n-vector stored in dx() with storage
c     increment incx .
c     if    n .le. 0 return with result = 0.
c     if n .ge. 1 then incx must be .ge. 1
c
c           c.l.lawson, 1978 jan 08
c     modified to correct failure to update ix, 1/25/92.
c     modified 3/93 to return if incx .le. 0.
c
c     four phase method     using two built-in constants that are
c     hopefully applicable to all machines.
c         cutlo = maximum of  dsqrt(u/eps)  over all known machines.
c         cuthi = minimum of  dsqrt(v)      over all known machines.
c     where
c         eps = smallest no. such that eps + 1. .gt. 1.
c         u   = smallest positive no.   (underflow limit)
c         v   = largest  no.            (overflow  limit)
c
c     brief outline of algorithm..
c
c     phase 1    scans zero components.
c     move to phase 2 when a component is nonzero and .le. cutlo
c     move to phase 3 when a component is .gt. cutlo
c     move to phase 4 when a component is .ge. cuthi/m
c     where m = n for x() real and m = 2*n for complex.
c
c     values for cutlo and cuthi..
c     from the environmental parameters listed in the imsl converter
c     document the limiting values are as follows..
c     cutlo, s.p.   u/eps = 2**(-102) for  honeywell.  close seconds are
c                   univac and dec at 2**(-103)
c                   thus cutlo = 2**(-51) = 4.44089e-16
c     cuthi, s.p.   v = 2**127 for univac, honeywell, and dec.
c                   thus cuthi = 2**(63.5) = 1.30438e19
c     cutlo, d.p.   u/eps = 2**(-67) for honeywell and dec.
c                   thus cutlo = 2**(-33.5) = 8.23181d-11
c     cuthi, d.p.   same as s.p.  cuthi = 1.30438d19
c     data cutlo, cuthi / 8.232d-11,  1.304d19 /
c     data cutlo, cuthi / 4.441e-16,  1.304e19 /
      data cutlo, cuthi / 8.232d-11,  1.304d19 /
c
      if(n .gt. 0 .and. incx.gt.0) go to 10
         dnrm2  = zero
         go to 300
c
   10 assign 30 to next
      sum = zero
      i = 1
      ix = 1
c                                                 begin main loop
   20    go to next,(30, 50, 70, 110)
   30 if( dabs(dx(i)) .gt. cutlo) go to 85
      assign 50 to next
      xmax = zero
c
c                        phase 1.  sum is zero
c
   50 if( dx(i) .eq. zero) go to 200
      if( dabs(dx(i)) .gt. cutlo) go to 85
c
c                                prepare for phase 2.
      assign 70 to next
      go to 105
c
c                                prepare for phase 4.
c
  100 continue
      ix = j
      assign 110 to next
      sum = (sum / dx(i)) / dx(i)
  105 xmax = dabs(dx(i))
      go to 115
c
c                   phase 2.  sum is small.
c                             scale to avoid destructive underflow.
c
   70 if( dabs(dx(i)) .gt. cutlo ) go to 75
c
c                     common code for phases 2 and 4.
c                     in phase 4 sum is large.  scale to avoid overflow.
c
  110 if( dabs(dx(i)) .le. xmax ) go to 115
         sum = one + sum * (xmax / dx(i))**2
         xmax = dabs(dx(i))
         go to 200
c
  115 sum = sum + (dx(i)/xmax)**2
      go to 200
c
c
c                  prepare for phase 3.
c
   75 sum = (sum * xmax) * xmax
c
c
c     for real or d.p. set hitest = cuthi/n
c     for complex      set hitest = cuthi/(2*n)
c
   85 hitest = cuthi/float( n )
c
c                   phase 3.  sum is mid-range.  no scaling.
c
      do 95 j = ix,n
      if(dabs(dx(i)) .ge. hitest) go to 100
         sum = sum + dx(i)**2
         i = i + incx
   95 continue
      dnrm2 = dsqrt( sum )
      go to 300
c
  200 continue
      ix = ix + 1
      i = i + incx
      if( ix .le. n ) go to 20
c
c              end of main loop.
c
c              compute square root and adjust for scaling.
c
      dnrm2 = xmax * dsqrt(sum)
  300 continue
      return
      end
      subroutine dscal ( n, da, dx, incx )

c*********************************************************************72
c
cc DSCAL scales a vector by a constant.
c
c     scales a vector by a constant.
c     uses unrolled loops for increment equal to one.
c     jack dongarra, linpack, 3/11/78.
c     modified 3/93 to return if incx .le. 0.
c
      implicit none

      double precision da,dx(1)
      integer i,incx,m,mp1,n,nincx

      if( n.le.0 .or. incx.le.0 )return
      if(incx.eq.1)go to 20
c
c        code for increment not equal to 1
c
      nincx = n*incx
      do 10 i = 1,nincx,incx
        dx(i) = da*dx(i)
   10 continue
      return
c
c        code for increment equal to 1
c
c
c        clean-up loop
c
   20 m = mod(n,5)
      if( m .eq. 0 ) go to 40
      do 30 i = 1,m
        dx(i) = da*dx(i)
   30 continue
      if( n .lt. 5 ) return
   40 mp1 = m + 1
      do 50 i = mp1,n,5
        dx(i) = da*dx(i)
        dx(i + 1) = da*dx(i + 1)
        dx(i + 2) = da*dx(i + 2)
        dx(i + 3) = da*dx(i + 3)
        dx(i + 4) = da*dx(i + 4)
   50 continue
      return
      end
      subroutine dschbv ( m, t, h,ldh, y, wsp, iwsp, iflag )

c*********************************************************************72
c
cc DSCHBV computes y = exp(t*H)*y using a partial fraction expansion.
c
      implicit none

      integer          m, ldh, iflag, iwsp(m)
      double precision t, H(ldh,m), y(m)
      complex*16       wsp(m*(m+2))

*-----Purpose----------------------------------------------------------|
*
*---  DSCHBV computes y = exp(t*H)*y using the partial fraction
*     expansion of the uniform rational Chebyshev approximation
*     to exp(-x) of type (14,14). H is assumed to be symmetric.
*     About 14-digit accuracy is expected if the matrix H is negative
*     definite. The algorithm may behave poorly otherwise.
*
*-----Arguments--------------------------------------------------------|
*
*     m       : (input) order of matrix H
*
*     t       : (input) time-scaling factor (can be < 0).
*
*     H(ldh,m): (input) symmetric matrix.
*
*     y(m)    : (input/output) on input the operand vector,
*               on output the resulting vector exp(t*H)*y.
*
*     iwsp(m) : (workspace)
*
*     wsp     : (workspace). Observe that a double precision vector of
*               length 2*m*(m+2) can be used as well when calling this
*               routine (thus avoiding an idle complex array elsewhere)
*
*----------------------------------------------------------------------|
*     Roger B. Sidje (rbs@maths.uq.edu.au)
*     EXPOKIT: Software Package for Computing Matrix Exponentials.
*     ACM - Transactions On Mathematical Software, 24(1):130-156, 1998
*----------------------------------------------------------------------|
*
      integer ndeg, i, j, ip, ih, iy, iz
      parameter ( ndeg=7 )
      double precision alpha0
      complex*16 alpha(ndeg), theta(ndeg), w

      intrinsic ABS,CMPLX,DBLE,MIN

*---  Pointers ...

      ih = 1
      iy = ih + m*m
      iz = iy + m

*---  Coefficients and poles of the partial fraction expansion ...

      alpha0  =  0.183216998528140087D-11
      alpha(1)=( 0.557503973136501826D+02,-0.204295038779771857D+03)
      alpha(2)=(-0.938666838877006739D+02, 0.912874896775456363D+02)
      alpha(3)=( 0.469965415550370835D+02,-0.116167609985818103D+02)
      alpha(4)=(-0.961424200626061065D+01,-0.264195613880262669D+01)
      alpha(5)=( 0.752722063978321642D+00, 0.670367365566377770D+00)
      alpha(6)=(-0.188781253158648576D-01,-0.343696176445802414D-01)
      alpha(7)=( 0.143086431411801849D-03, 0.287221133228814096D-03)

      theta(1)=(-0.562314417475317895D+01, 0.119406921611247440D+01)
      theta(2)=(-0.508934679728216110D+01, 0.358882439228376881D+01)
      theta(3)=(-0.399337136365302569D+01, 0.600483209099604664D+01)
      theta(4)=(-0.226978543095856366D+01, 0.846173881758693369D+01)
      theta(5)=( 0.208756929753827868D+00, 0.109912615662209418D+02)
      theta(6)=( 0.370327340957595652D+01, 0.136563731924991884D+02)
      theta(7)=( 0.889777151877331107D+01, 0.166309842834712071D+02)
*
*---  Accumulation of the contribution of each pole ...
*
      do j = 1,m
         wsp(iz+j-1) = y(j)
         y(j) = y(j)*alpha0
      end do
      do ip = 1,ndeg
*---     Solve each fraction using Gaussian elimination with pivoting...
         do j = 1,m
            do i = 1,m
               wsp(ih+(j-1)*m+i-1) = -t*H(i,j)
            end do
            wsp(ih+(j-1)*m+j-1) = wsp(ih+(j-1)*m+j-1)-theta(ip)
            wsp(iy+j-1) = wsp(iz+j-1)
         end do
         call ZSYSV('U', M, 1, WSP(iH),M, IWSP, WSP(iY),M, W,1, IFLAG )
         if ( IFLAG.ne.0 ) stop 'Error in DSCHBV'
*---     Accumulate the partial result in y ...
         do i = 1,m
            y(i) = y(i) + DBLE( alpha(ip)*wsp(iy+i-1) )
         end do
      end do
      END
      subroutine dsexpv ( n, m, t, v, w, tol, anorm,
     .                   wsp,lwsp, iwsp,liwsp, matvec, itrace,iflag )

c*********************************************************************72
c
cc DSEXPV computes w = exp(t*A)*v, for a symmetric matrix A.
c
      implicit none

      integer n, m, lwsp, liwsp, itrace, iflag, iwsp(liwsp)
      double precision t, tol, anorm, v(n), w(n), wsp(lwsp)
      external matvec

*-----Purpose----------------------------------------------------------|
*
*---  DSEXPV computes w = exp(t*A)*v - for a Symmetric matrix A.
*
*     It does not compute the matrix exponential in isolation but
*     instead, it computes directly the action of the exponential
*     operator on the operand vector. This way of doing so allows
*     for addressing large sparse problems.
*
*     The method used is based on Krylov subspace projection
*     techniques and the matrix under consideration interacts only
*     via the external routine `matvec' performing the matrix-vector
*     product (matrix-free method).
*
*-----Arguments--------------------------------------------------------|
*
*     n      : (input) order of the principal matrix A.
*
*     m      : (input) maximum size for the Krylov basis.
*
*     t      : (input) time at wich the solution is needed (can be < 0).
*
*     v(n)   : (input) given operand vector.
*
*     w(n)   : (output) computed approximation of exp(t*A)*v.
*
*     tol    : (input/output) the requested accuracy tolerance on w.
*              If on input tol=0.0d0 or tol is too small (tol.le.eps)
*              the internal value sqrt(eps) is used, and tol is set to
*              sqrt(eps) on output (`eps' denotes the machine epsilon).
*              (`Happy breakdown' is assumed if h(j+1,j) .le. anorm*tol)
*
*     anorm  : (input) an approximation of some norm of A.
*
*   wsp(lwsp): (workspace) lwsp .ge. n*(m+1)+n+(m+2)^2+4*(m+2)^2+ideg+1
*                                   +---------+-------+---------------+
*              (actually, ideg=6)        V        H      wsp for PADE
*
* iwsp(liwsp): (workspace) liwsp .ge. m+2
*
*     matvec : external routine for matrix-vector multiplication.
*              synopsis: matvec( x, y )
*                        double precision x(*), y(*)
*              computes: y(1:n) <- A*x(1:n)
*                        where A is the principal matrix.
*
*     itrace : (input) running mode. 0=silent, 1=print step-by-step info
*
*     iflag  : (output) exit flag.
*              <0 - bad input arguments
*               0 - no problem
*               1 - maximum number of steps reached without convergence
*               2 - requested tolerance was too high
*
*-----Accounts on the computation--------------------------------------|
*     Upon exit, an interested user may retrieve accounts on the
*     computations. They are located in the workspace arrays wsp and
*     iwsp as indicated below:
*
*     location  mnemonic                 description
*     -----------------------------------------------------------------|
*     iwsp(1) = nmult, number of matrix-vector multiplications used
*     iwsp(2) = nexph, nbr of Tridiagonal matrix exponential evaluated
*     iwsp(3) = nscale, number of repeated squaring involved in Pade
*     iwsp(4) = nstep, nbr of integration steps used up to completion
*     iwsp(5) = nreject, number of rejected step-sizes
*     iwsp(6) = ibrkflag, set to 1 if `happy breakdown' and 0 otherwise
*     iwsp(7) = mbrkdwn, if `happy brkdown', basis-size when it occured
*     -----------------------------------------------------------------|
*     wsp(1)  = step_min, minimum step-size used during integration
*     wsp(2)  = step_max, maximum step-size used during integration
*     wsp(3)  = dummy
*     wsp(4)  = dummy
*     wsp(5)  = x_error, maximum among all local truncation errors
*     wsp(6)  = s_error, global sum of local truncation errors
*     wsp(7)  = tbrkdwn, if `happy breakdown', time when it occured
*     wsp(8)  = t_now, integration domain successfully covered
*     wsp(9)  = hump, i.e., max||exp(sA)||, s in [0,t] (or [t,0] if t<0)
*     wsp(10) = ||w||/||v||, scaled norm of the solution w.
*     -----------------------------------------------------------------|
*     The `hump' is a measure of the conditioning of the problem. The
*     matrix exponential is well-conditioned if hump = 1, whereas it is
*     poorly-conditioned if hump >> 1. However the solution can still be
*     relatively fairly accurate even when the hump is large (the hump
*     is an upper bound), especially when the hump and the scaled norm
*     of w [this is also computed and returned in wsp(10)] are of the
*     same order of magnitude (further details in reference below).
*
*----------------------------------------------------------------------|
*-----The following parameters may also be adjusted herein-------------|
*
      integer mxstep, mxreject, ideg
      double precision delta, gamma
      parameter( mxstep   = 500,
     .           mxreject = 0,
     .           ideg     = 6,
     .           delta    = 1.2d0,
     .           gamma    = 0.9d0 )

*     mxstep  : maximum allowable number of integration steps.
*               The value 0 means an infinite number of steps.
*
*     mxreject: maximum allowable number of rejections at each step.
*               The value 0 means an infinite number of rejections.
*
*     ideg    : the Pade approximation of type (ideg,ideg) is used as
*               an approximation to exp(H). The value 0 switches to the
*               uniform rational Chebyshev approximation of type (14,14)
*
*     delta   : local truncation error `safety factor'
*
*     gamma   : stepsize `shrinking factor'
*
*----------------------------------------------------------------------|
*     Roger B. Sidje (rbs@maths.uq.edu.au)
*     EXPOKIT: Software Package for Computing Matrix Exponentials.
*     ACM - Transactions On Mathematical Software, 24(1):130-156, 1998
*----------------------------------------------------------------------|
*
      integer i, j, k1, mh, mx, iv, ih, j1v, ns, ifree, lfree, iexph,
     .        ireject,ibrkflag,mbrkdwn, nmult, nreject, nexph, nscale,
     .        nstep
      double precision sgn, t_out, tbrkdwn, step_min,step_max, err_loc,
     .                 s_error, x_error, t_now, t_new, t_step, t_old,
     .                 xm, beta, break_tol, p1, p2, p3, eps, rndoff,
     .                 vnorm, avnorm, hj1j, hjj, hump, SQR1

      intrinsic AINT,ABS,DBLE,LOG10,MAX,MIN,NINT,SIGN,SQRT
      double precision DDOT, DNRM2

*---  check restrictions on input parameters ...
      iflag = 0
      if ( lwsp.lt.n*(m+2)+5*(m+2)**2+ideg+1 ) iflag = -1
      if ( liwsp.lt.m+2 ) iflag = -2
      if ( m.ge.n .or. m.le.0 ) iflag = -3
      if ( iflag.ne.0 ) stop 'bad sizes (in input of DSEXPV)'
*
*---  initialisations ...
*
      k1 = 2
      mh = m + 2
      iv = 1
      ih = iv + n*(m+1) + n
      ifree = ih + mh*mh
      lfree = lwsp - ifree + 1

      ibrkflag = 0
      mbrkdwn  = m
      nmult    = 0
      nreject  = 0
      nexph    = 0
      nscale   = 0

      t_out    = ABS( t )
      tbrkdwn  = 0.0d0
      step_min = t_out
      step_max = 0.0d0
      nstep    = 0
      s_error  = 0.0d0
      x_error  = 0.0d0
      t_now    = 0.0d0
      t_new    = 0.0d0

      p1 = 4.0d0/3.0d0
 1    p2 = p1 - 1.0d0
      p3 = p2 + p2 + p2
      eps = ABS( p3-1.0d0 )
      if ( eps.eq.0.0d0 ) go to 1
      if ( tol.le.eps ) tol = SQRT( eps )
      rndoff = eps*anorm

      break_tol = 1.0d-7
*>>>  break_tol = tol
*>>>  break_tol = anorm*tol

      sgn = SIGN( 1.0d0,t )
      call DCOPY( n, v,1, w,1 )
      beta = DNRM2( n, w,1 )
      vnorm = beta
      hump = beta
*
*---  obtain the very first stepsize ...
*
      SQR1 = SQRT( 0.1d0 )
      xm = 1.0d0/DBLE( m )
      p1 = tol*(((m+1)/2.72D0)**(m+1))*SQRT(2.0D0*3.14D0*(m+1))
      t_new = (1.0d0/anorm)*(p1/(4.0d0*beta*anorm))**xm
      p1 = 10.0d0**(NINT( LOG10( t_new )-SQR1 )-1)
      t_new = AINT( t_new/p1 + 0.55d0 ) * p1
*
*---  step-by-step integration ...
*
 100  if ( t_now.ge.t_out ) goto 500

      nstep = nstep + 1
      t_step = MIN( t_out-t_now, t_new )

      p1 = 1.0d0/beta
      do i = 1,n
         wsp(iv + i-1) = p1*w(i)
      end do
      do i = 1,mh*mh
         wsp(ih+i-1) = 0.0d0
      end do
*
*---  Lanczos loop ...
*
      j1v = iv + n
      do 200 j = 1,m
         nmult = nmult + 1
         call matvec( wsp(j1v-n), wsp(j1v) )
         if ( j.gt.1 )
     .     call DAXPY(n,-wsp(ih+(j-1)*mh+j-2),wsp(j1v-2*n),1,wsp(j1v),1)
         hjj = DDOT( n, wsp(j1v-n),1, wsp(j1v),1 )
         call DAXPY( n, -hjj, wsp(j1v-n),1, wsp(j1v),1 )
         hj1j = DNRM2( n, wsp(j1v),1 )
         wsp(ih+(j-1)*(mh+1)) = hjj
*---     if `happy breakdown' go straightforward at the end ...
         if ( hj1j.le.break_tol ) then
            print*,'happy breakdown: mbrkdwn =',j,' h =',hj1j
            k1 = 0
            ibrkflag = 1
            mbrkdwn = j
            tbrkdwn = t_now
            t_step = t_out-t_now
            goto 300
         endif
         wsp(ih+(j-1)*mh+j) = hj1j
         wsp(ih+j*mh+j-1) = hj1j
         call DSCAL( n, 1.0d0/hj1j, wsp(j1v),1 )
         j1v = j1v + n
 200  continue
      nmult = nmult + 1
      call matvec( wsp(j1v-n), wsp(j1v) )
      avnorm = DNRM2( n, wsp(j1v),1 )
*
*---  set 1 for the 2-corrected scheme ...
*
 300  continue
      wsp(ih+m*mh+m-1) = 0.0d0
      wsp(ih+m*mh+m+1) = 1.0d0
*
*---  loop while ireject<mxreject until the tolerance is reached ...
*
      ireject = 0
 401  continue
*
*---  compute w = beta*V*exp(t_step*H)*e1 ...
*
      nexph = nexph + 1
      mx = mbrkdwn + k1
      if ( ideg.ne.0 ) then
*---     irreducible rational Pade approximation ...
         call DGPADM( ideg, mx, sgn*t_step, wsp(ih),mh,
     .                wsp(ifree),lfree, iwsp, iexph, ns, iflag )
         iexph = ifree + iexph - 1
         nscale = nscale + ns
      else
*---     uniform rational Chebyshev approximation ...
         iexph = ifree
         do i = 1,mx
            wsp(iexph+i-1) = 0.0d0
         end do
         wsp(iexph) = 1.0d0
         call DNCHBV(mx,sgn*t_step,wsp(ih),mh,wsp(iexph),wsp(ifree+mx))
      endif
 402  continue
*
*---  error estimate ...
*
      if ( k1.eq.0 ) then
         err_loc = tol
      else
         p1 = ABS( wsp(iexph+m) )   * beta
         p2 = ABS( wsp(iexph+m+1) ) * beta * avnorm
         if ( p1.gt.10.0d0*p2 ) then
            err_loc = p2
            xm = 1.0d0/DBLE( m )
         elseif ( p1.gt.p2 ) then
            err_loc = (p1*p2)/(p1-p2)
            xm = 1.0d0/DBLE( m )
         else
            err_loc = p1
            xm = 1.0d0/DBLE( m-1 )
         endif
      endif
*
*---  reject the step-size if the error is not acceptable ...
*
      if ( (k1.ne.0) .and. (err_loc.gt.delta*t_step*tol) .and.
     .     (mxreject.eq.0 .or. ireject.lt.mxreject) ) then
         t_old = t_step
         t_step = gamma * t_step * (t_step*tol/err_loc)**xm
         p1 = 10.0d0**(NINT( LOG10( t_step )-SQR1 )-1)
         t_step = AINT( t_step/p1 + 0.55d0 ) * p1
         if ( itrace.ne.0 ) then
            print*,'t_step =',t_old
            print*,'err_loc =',err_loc
            print*,'err_required =',delta*t_old*tol
            print*,'stepsize rejected, stepping down to:',t_step
         endif
         ireject = ireject + 1
         nreject = nreject + 1
         if ( mxreject.ne.0 .and. ireject.gt.mxreject ) then
            print*,"Failure in DSEXPV: ---"
            print*,"The requested tolerance is too high."
            Print*,"Rerun with a smaller value."
            iflag = 2
            return
         endif
         goto 401
      endif
*
*---  now update w = beta*V*exp(t_step*H)*e1 and the hump ...
*
      mx = mbrkdwn + MAX( 0,k1-1 )
      call DGEMV( 'n', n,mx,beta,wsp(iv),n,wsp(iexph),1,0.0d0,w,1 )
      beta = DNRM2( n, w,1 )
      hump = MAX( hump, beta )
*
*---  suggested value for the next stepsize ...
*
      t_new = gamma * t_step * (t_step*tol/err_loc)**xm
      p1 = 10.0d0**(NINT( LOG10( t_new )-SQR1 )-1)
      t_new = AINT( t_new/p1 + 0.55d0 ) * p1

      err_loc = MAX( err_loc,rndoff )
*
*---  update the time covered ...
*
      t_now = t_now + t_step
*
*---  display and keep some information ...
*
      if ( itrace.ne.0 ) then
         print*,'integration',nstep,'---------------------------------'
         print*,'scale-square =',ns
         print*,'step_size =',t_step
         print*,'err_loc   =',err_loc
         print*,'next_step =',t_new
      endif

      step_min = MIN( step_min, t_step )
      step_max = MAX( step_max, t_step )
      s_error = s_error + err_loc
      x_error = MAX( x_error, err_loc )

      if ( mxstep.eq.0 .or. nstep.lt.mxstep ) goto 100
      iflag = 1

 500  continue

      iwsp(1) = nmult
      iwsp(2) = nexph
      iwsp(3) = nscale
      iwsp(4) = nstep
      iwsp(5) = nreject
      iwsp(6) = ibrkflag
      iwsp(7) = mbrkdwn

      wsp(1)  = step_min
      wsp(2)  = step_max
      wsp(3)  = 0.0d0
      wsp(4)  = 0.0d0
      wsp(5)  = x_error
      wsp(6)  = s_error
      wsp(7)  = tbrkdwn
      wsp(8)  = sgn*t_now
      wsp(9)  = hump/vnorm
      wsp(10) = beta/vnorm
      END
      subroutine dsifa(a,lda,n,kpvt,info)

c*********************************************************************72
c
cc DSIFA factors a symmetric matrix.
c
      implicit none

      integer lda,n,kpvt(1),info
      double precision a(lda,1)
c
c     dsifa factors a double precision symmetric matrix by elimination
c     with symmetric pivoting.
c
c     to solve  a*x = b , follow dsifa by dsisl.
c     to compute  inverse(a)*c , follow dsifa by dsisl.
c     to compute  determinant(a) , follow dsifa by dsidi.
c     to compute  inertia(a) , follow dsifa by dsidi.
c     to compute  inverse(a) , follow dsifa by dsidi.
c
c     on entry
c
c        a       double precision(lda,n)
c                the symmetric matrix to be factored.
c                only the diagonal and upper triangle are used.
c
c        lda     integer
c                the leading dimension of the array  a .
c
c        n       integer
c                the order of the matrix  a .
c
c     on return
c
c        a       a block diagonal matrix and the multipliers which
c                were used to obtain it.
c                the factorization can be written  a = u*d*trans(u)
c                where  u  is a product of permutation and unit
c                upper triangular matrices , trans(u) is the
c                transpose of  u , and  d  is block diagonal
c                with 1 by 1 and 2 by 2 blocks.
c
c        kpvt    integer(n)
c                an integer vector of pivot indices.
c
c        info    integer
c                = 0  normal value.
c                = k  if the k-th pivot block is singular. this is
c                     not an error condition for this subroutine,
c                     but it does indicate that dsisl or dsidi may
c                     divide by zero if called.
c
c     linpack. this version dated 08/14/78 .
c     james bunch, univ. calif. san diego, argonne nat. lab.
c
c     subroutines and functions
c
c     blas daxpy,dswap,idamax
c     fortran dabs,dmax1,dsqrt
c
c     internal variables
c
      double precision ak,akm1,bk,bkm1,denom,mulk,mulkm1,t
      double precision absakk,alpha,colmax,rowmax
      integer imax,imaxp1,j,jj,jmax,k,km1,km2,kstep,idamax
      logical swap
c
c
c     initialize
c
c     alpha is used in choosing pivot block size.
      alpha = (1.0d0 + dsqrt(17.0d0))/8.0d0
c
      info = 0
c
c     main loop on k, which goes from n to 1.
c
      k = n
   10 continue
c
c        leave the loop if k=0 or k=1.
c
c     ...exit
         if (k .eq. 0) go to 200
         if (k .gt. 1) go to 20
            kpvt(1) = 1
            if (a(1,1) .eq. 0.0d0) info = 1
c     ......exit
            go to 200
   20    continue
c
c        this section of code determines the kind of
c        elimination to be performed.  when it is completed,
c        kstep will be set to the size of the pivot block, and
c        swap will be set to .true. if an interchange is
c        required.
c
         km1 = k - 1
         absakk = dabs(a(k,k))
c
c        determine the largest off-diagonal element in
c        column k.
c
         imax = idamax(k-1,a(1,k),1)
         colmax = dabs(a(imax,k))
         if (absakk .lt. alpha*colmax) go to 30
            kstep = 1
            swap = .false.
         go to 90
   30    continue
c
c           determine the largest off-diagonal element in
c           row imax.
c
            rowmax = 0.0d0
            imaxp1 = imax + 1
            do 40 j = imaxp1, k
               rowmax = dmax1(rowmax,dabs(a(imax,j)))
   40       continue
            if (imax .eq. 1) go to 50
               jmax = idamax(imax-1,a(1,imax),1)
               rowmax = dmax1(rowmax,dabs(a(jmax,imax)))
   50       continue
            if (dabs(a(imax,imax)) .lt. alpha*rowmax) go to 60
               kstep = 1
               swap = .true.
            go to 80
   60       continue
            if (absakk .lt. alpha*colmax*(colmax/rowmax)) go to 70
               kstep = 1
               swap = .false.
            go to 80
   70       continue
               kstep = 2
               swap = imax .ne. km1
   80       continue
   90    continue
         if (dmax1(absakk,colmax) .ne. 0.0d0) go to 100
c
c           column k is zero.  set info and iterate the loop.
c
            kpvt(k) = k
            info = k
         go to 190
  100    continue
         if (kstep .eq. 2) go to 140
c
c           1 x 1 pivot block.
c
            if (.not.swap) go to 120
c
c              perform an interchange.
c
               call dswap(imax,a(1,imax),1,a(1,k),1)
               do 110 jj = imax, k
                  j = k + imax - jj
                  t = a(j,k)
                  a(j,k) = a(imax,j)
                  a(imax,j) = t
  110          continue
  120       continue
c
c           perform the elimination.
c
            do 130 jj = 1, km1
               j = k - jj
               mulk = -a(j,k)/a(k,k)
               t = mulk
               call daxpy(j,t,a(1,k),1,a(1,j),1)
               a(j,k) = mulk
  130       continue
c
c           set the pivot array.
c
            kpvt(k) = k
            if (swap) kpvt(k) = imax
         go to 190
  140    continue
c
c           2 x 2 pivot block.
c
            if (.not.swap) go to 160
c
c              perform an interchange.
c
               call dswap(imax,a(1,imax),1,a(1,k-1),1)
               do 150 jj = imax, km1
                  j = km1 + imax - jj
                  t = a(j,k-1)
                  a(j,k-1) = a(imax,j)
                  a(imax,j) = t
  150          continue
               t = a(k-1,k)
               a(k-1,k) = a(imax,k)
               a(imax,k) = t
  160       continue
c
c           perform the elimination.
c
            km2 = k - 2
            if (km2 .eq. 0) go to 180
               ak = a(k,k)/a(k-1,k)
               akm1 = a(k-1,k-1)/a(k-1,k)
               denom = 1.0d0 - ak*akm1
               do 170 jj = 1, km2
                  j = km1 - jj
                  bk = a(j,k)/a(k-1,k)
                  bkm1 = a(j,k-1)/a(k-1,k)
                  mulk = (akm1*bk - bkm1)/denom
                  mulkm1 = (ak*bkm1 - bk)/denom
                  t = mulk
                  call daxpy(j,t,a(1,k),1,a(1,j),1)
                  t = mulkm1
                  call daxpy(j,t,a(1,k-1),1,a(1,j),1)
                  a(j,k) = mulk
                  a(j,k-1) = mulkm1
  170          continue
  180       continue
c
c           set the pivot array.
c
            kpvt(k) = 1 - k
            if (swap) kpvt(k) = -imax
            kpvt(k-1) = kpvt(k)
  190    continue
         k = k - kstep
      go to 10
  200 continue
      return
      end
      subroutine dsisl(a,lda,n,kpvt,b)

c*********************************************************************72
c
cc DSISL solves a factored symmetric linear system.
c
      implicit none

      integer lda,n,kpvt(1)
      double precision a(lda,1),b(1)
c
c     dsisl solves the double precision symmetric system
c     a * x = b
c     using the factors computed by dsifa.
c
c     on entry
c
c        a       double precision(lda,n)
c                the output from dsifa.
c
c        lda     integer
c                the leading dimension of the array  a .
c
c        n       integer
c                the order of the matrix  a .
c
c        kpvt    integer(n)
c                the pivot vector from dsifa.
c
c        b       double precision(n)
c                the right hand side vector.
c
c     on return
c
c        b       the solution vector  x .
c
c     error condition
c
c        a division by zero may occur if  dsico  has set rcond .eq. 0.0
c        or  dsifa  has set info .ne. 0  .
c
c     to compute  inverse(a) * c  where  c  is a matrix
c     with  p  columns
c           call dsifa(a,lda,n,kpvt,info)
c           if (info .ne. 0) go to ...
c           do 10 j = 1, p
c              call dsisl(a,lda,n,kpvt,c(1,j))
c        10 continue
c
c     linpack. this version dated 08/14/78 .
c     james bunch, univ. calif. san diego, argonne nat. lab.
c
c     subroutines and functions
c
c     blas daxpy,ddot
c     fortran iabs
c
c     internal variables.
c
      double precision ak,akm1,bk,bkm1,ddot,denom,temp
      integer k,kp
c
c     loop backward applying the transformations and
c     d inverse to b.
c
      k = n
   10 if (k .eq. 0) go to 80
         if (kpvt(k) .lt. 0) go to 40
c
c           1 x 1 pivot block.
c
            if (k .eq. 1) go to 30
               kp = kpvt(k)
               if (kp .eq. k) go to 20
c
c                 interchange.
c
                  temp = b(k)
                  b(k) = b(kp)
                  b(kp) = temp
   20          continue
c
c              apply the transformation.
c
               call daxpy(k-1,b(k),a(1,k),1,b(1),1)
   30       continue
c
c           apply d inverse.
c
            b(k) = b(k)/a(k,k)
            k = k - 1
         go to 70
   40    continue
c
c           2 x 2 pivot block.
c
            if (k .eq. 2) go to 60
               kp = iabs(kpvt(k))
               if (kp .eq. k - 1) go to 50
c
c                 interchange.
c
                  temp = b(k-1)
                  b(k-1) = b(kp)
                  b(kp) = temp
   50          continue
c
c              apply the transformation.
c
               call daxpy(k-2,b(k),a(1,k),1,b(1),1)
               call daxpy(k-2,b(k-1),a(1,k-1),1,b(1),1)
   60       continue
c
c           apply d inverse.
c
            ak = a(k,k)/a(k-1,k)
            akm1 = a(k-1,k-1)/a(k-1,k)
            bk = b(k)/a(k-1,k)
            bkm1 = b(k-1)/a(k-1,k)
            denom = ak*akm1 - 1.0d0
            b(k) = (akm1*bk - bkm1)/denom
            b(k-1) = (ak*bkm1 - bk)/denom
            k = k - 2
   70    continue
      go to 10
   80 continue
c
c     loop forward applying the transformations.
c
      k = 1
   90 if (k .gt. n) go to 160
         if (kpvt(k) .lt. 0) go to 120
c
c           1 x 1 pivot block.
c
            if (k .eq. 1) go to 110
c
c              apply the transformation.
c
               b(k) = b(k) + ddot(k-1,a(1,k),1,b(1),1)
               kp = kpvt(k)
               if (kp .eq. k) go to 100
c
c                 interchange.
c
                  temp = b(k)
                  b(k) = b(kp)
                  b(kp) = temp
  100          continue
  110       continue
            k = k + 1
         go to 150
  120    continue
c
c           2 x 2 pivot block.
c
            if (k .eq. 1) go to 140
c
c              apply the transformation.
c
               b(k) = b(k) + ddot(k-1,a(1,k),1,b(1),1)
               b(k+1) = b(k+1) + ddot(k-1,a(1,k+1),1,b(1),1)
               kp = iabs(kpvt(k))
               if (kp .eq. k) go to 130
c
c                 interchange.
c
                  temp = b(k)
                  b(k) = b(kp)
                  b(kp) = temp
  130          continue
  140       continue
            k = k + 2
  150    continue
      go to 90
  160 continue
      return
      end
      subroutine dspadm ( ideg,m,t,h,ldh,wsp,lwsp,ipiv,iexph,ns,iflag )

c*********************************************************************72
c
cc DSPADM computes exp(t*H), the matrix exponential of a symmetric matrix.
c
      implicit none

      integer ideg, m, ldh, lwsp, iexph, ns, iflag, ipiv(m)
      double precision t, H(ldh,m), wsp(lwsp)

*-----Purpose----------------------------------------------------------|
*
*     Computes exp(t*H), the matrix exponential of a symmetric matrix
*     in full, using the irreducible rational Pade approximation to the
*     exponential function exp(x) = r(x) = (+/-)( I + 2*(q(x)/p(x)) ),
*     combined with scaling-and-squaring.
*
*-----Arguments--------------------------------------------------------|
*
*     ideg      : (input) the degre of the diagonal Pade to be used.
*                 a value of 6 is generally satisfactory.
*
*     m         : (input) order of H.
*
*     H(ldh,m)  : (input) argument matrix (both lower and upper parts).
*
*     t         : (input) time-scale (can be < 0).
*
*     wsp(lwsp) : (workspace/output) lwsp .ge. 4*m*m+ideg+1.
*
*     ipiv(m)   : (workspace)
*
*>>>> iexph     : (output) number such that wsp(iexph) points to exp(tH)
*                 i.e., exp(tH) is located at wsp(iexph ... iexph+m*m-1)
*                       ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
*                 NOTE: if the routine was called with wsp(iptr),
*                       then exp(tH) will start at wsp(iptr+iexph-1).
*
*     ns        : (output) number of scaling-squaring used.
*
*     iflag     : (output) exit flag.
*                      0 - no problem
*                     <0 - problem
*
*----------------------------------------------------------------------|
*     Roger B. Sidje (rbs@maths.uq.edu.au)
*     EXPOKIT: Software Package for Computing Matrix Exponentials.
*     ACM - Transactions On Mathematical Software, 24(1):130-156, 1998
*----------------------------------------------------------------------|
*
      integer mm,i,j,k,ih2,ip,iq,iused,ifree,iodd,icoef,iput,iget
      double precision hnorm,scale,scale2,cp,cq

      intrinsic INT,ABS,DBLE,LOG,MAX

*---  check restrictions on input parameters ...
      mm = m*m
      iflag = 0
      if ( ldh.lt.m ) iflag = -1
      if ( lwsp.lt.4*mm+ideg+1 ) iflag = -2
      if ( iflag.ne.0 ) stop 'bad sizes (in input of DSPADM)'
*
*---  initialise pointers ...
*
      icoef = 1
      ih2 = icoef + (ideg+1)
      ip  = ih2 + mm
      iq  = ip + mm
      ifree = iq + mm
*
*---  scaling: seek ns such that ||t*H/2^ns|| < 1/2;
*     and set scale = t/2^ns ...
*
      do i = 1,m
         wsp(i) = 0.0d0
      end do
      do j = 1,m
         do i = 1,m
            wsp(i) = wsp(i) + ABS( H(i,j) )
         end do
      end do
      hnorm = 0.0d0
      do i = 1,m
         hnorm = MAX( hnorm,wsp(i) )
      end do
      hnorm = ABS( t*hnorm )
      if ( hnorm.eq.0.0d0 ) stop 'Error - null H in input of DSPADM.'
      ns = MAX( 0,INT(LOG(hnorm)/LOG(2.0d0))+2 )
      scale = t / DBLE(2**ns)
      scale2 = scale*scale
*
*---  compute Pade coefficients ...
*
      i = ideg+1
      j = 2*ideg+1
      wsp(icoef) = 1.0d0
      do k = 1,ideg
         wsp(icoef+k) = (wsp(icoef+k-1)*DBLE( i-k ))/DBLE( k*(j-k) )
      end do
*
*---  H2 = scale2*H*H ...
*
      call DGEMM( 'n','n',m,m,m,scale2,H,ldh,H,ldh,0.0d0,wsp(ih2),m )
*
*---  initialize p (numerator) and q (denominator) ...
*
      cp = wsp(icoef+ideg-1)
      cq = wsp(icoef+ideg)
      do j = 1,m
         do i = 1,m
            wsp(ip + (j-1)*m + i-1) = 0.0d0
            wsp(iq + (j-1)*m + i-1) = 0.0d0
         end do
         wsp(ip + (j-1)*(m+1)) = cp
         wsp(iq + (j-1)*(m+1)) = cq
      end do
*
*---  Apply Horner rule ...
*
      iodd = 1
      k = ideg - 1
 100  continue
      iused = iodd*iq + (1-iodd)*ip
      call DGEMM( 'n','n',m,m,m, 1.0d0,wsp(iused),m,
     .             wsp(ih2),m, 0.0d0,wsp(ifree),m )
      do j = 1,m
         wsp(ifree+(j-1)*(m+1)) = wsp(ifree+(j-1)*(m+1))+wsp(icoef+k-1)
      end do
      ip = (1-iodd)*ifree + iodd*ip
      iq = iodd*ifree + (1-iodd)*iq
      ifree = iused
      iodd = 1-iodd
      k = k-1
      if ( k.gt.0 )  goto 100
*
*---  Obtain (+/-)(I + 2*(p\q)) ...
*
      if ( iodd .eq. 1 ) then
         call DGEMM( 'n','n',m,m,m, scale,wsp(iq),m,
     .                H,ldh, 0.0d0,wsp(ifree),m )
         iq = ifree
      else
         call DGEMM( 'n','n',m,m,m, scale,wsp(ip),m,
     .                H,ldh, 0.0d0,wsp(ifree),m )
         ip = ifree
      endif
      call DAXPY( mm, -1.0d0,wsp(ip),1, wsp(iq),1 )
      call DSYSV( 'U',m,m,wsp(iq),m,ipiv,wsp(ip),m,wsp(ih2),mm,iflag )
      if ( iflag.ne.0 ) stop 'Problem in DSYSV (within DSPADM)'
      call DSCAL( mm, 2.0d0, wsp(ip), 1 )
      do j = 1,m
         wsp(ip+(j-1)*(m+1)) = wsp(ip+(j-1)*(m+1)) + 1.0d0
      end do
      iput = ip
      if ( ns.eq.0 .and. iodd.eq.1 ) then
         call DSCAL( mm, -1.0d0, wsp(ip), 1 )
         goto 200
      endif
*
*--   squaring : exp(t*H) = (exp(t*H))^(2^ns) ...
*
      iodd = 1
      do k = 1,ns
         iget = iodd*ip + (1-iodd)*iq
         iput = (1-iodd)*ip + iodd*iq
         call DGEMM( 'n','n',m,m,m, 1.0d0,wsp(iget),m, wsp(iget),m,
     .                0.0d0,wsp(iput),m )
         iodd = 1-iodd
      end do
 200  continue
      iexph = iput
      END
      subroutine dsphiv( n, m, t, u, v, w, tol, anorm,
     .                   wsp,lwsp, iwsp,liwsp, matvec, itrace,iflag )

c*********************************************************************72
c
cc DSPHIV solves the ODE w' = Aw + u, w(0) = v, for A symmetric.
c
      implicit none

      integer n, m, lwsp, liwsp, itrace, iflag, iwsp(liwsp)
      double precision t, tol, anorm, u(n), v(n), w(n), wsp(lwsp)
      external matvec

*-----Purpose----------------------------------------------------------|
*
*---  DSPHIV computes w = exp(t*A)v + t*phi(tA)u which is the solution
*     of the nonhomogeneous linear ODE problem w' = Aw + u, w(0) = v.
*     phi(z) = (exp(z)-1)/z and A is a Symmetric matrix.
*
*     The method used is based on Krylov subspace projection
*     techniques and the matrix under consideration interacts only
*     via the external routine `matvec' performing the matrix-vector
*     product (matrix-free method).
*
*-----Arguments--------------------------------------------------------|
*
*     n      : (input) order of the principal matrix A.
*
*     m      : (input) maximum size for the Krylov basis.
*
*     t      : (input) time at wich the solution is needed (can be < 0).
*
*     u(n)   : (input) operand vector with respect to the phi function
*              (forcing term of the ODE problem).
*
*     v(n)   : (input) operand vector with respect to the exp function
*              (initial condition of the ODE problem).
*
*     w(n)   : (output) computed approximation of exp(t*A)v + t*phi(tA)u
*
*     tol    : (input/output) the requested accuracy tolerance on w.
*              If on input tol=0.0d0 or tol is too small (tol.le.eps)
*              the internal value sqrt(eps) is used, and tol is set to
*              sqrt(eps) on output (`eps' denotes the machine epsilon).
*              (`Happy breakdown' is assumed if h(j+1,j) .le. anorm*tol)
*
*     anorm  : (input) an approximation of some norm of A.
*
*   wsp(lwsp): (workspace) lwsp .ge. n*(m+1)+n+(m+3)^2+4*(m+3)^2+ideg+1
*                                   +---------+-------+---------------+
*              (actually, ideg=6)        V        H      wsp for PADE
*
* iwsp(liwsp): (workspace) liwsp .ge. m+3
*
*     matvec : external routine for matrix-vector multiplication.
*              synopsis: matvec( x, y )
*                        double precision x(*), y(*)
*              computes: y(1:n) <- A*x(1:n)
*                        where A is the principal matrix.
*
*     itrace : (input) running mode. 0=silent, 1=print step-by-step info
*
*     iflag  : (output) exit flag.
*              <0 - bad input arguments
*               0 - no problem
*               1 - maximum number of steps reached without convergence
*               2 - requested tolerance was too high
*
*-----Accounts on the computation--------------------------------------|
*     Upon exit, an interested user may retrieve accounts on the
*     computations. They are located in the workspace arrays wsp and
*     iwsp as indicated below:
*
*     location  mnemonic                 description
*     -----------------------------------------------------------------|
*     iwsp(1) = nmult, number of matrix-vector multiplications used
*     iwsp(2) = nexph, number of Hessenberg matrix exponential evaluated
*     iwsp(3) = nscale, number of repeated squaring involved in Pade
*     iwsp(4) = nstep, number of integration steps used up to completion
*     iwsp(5) = nreject, number of rejected step-sizes
*     iwsp(6) = ibrkflag, set to 1 if `happy breakdown' and 0 otherwise
*     iwsp(7) = mbrkdwn, if `happy brkdown', basis-size when it occured
*     -----------------------------------------------------------------|
*     wsp(1)  = step_min, minimum step-size used during integration
*     wsp(2)  = step_max, maximum step-size used during integration
*     wsp(3)  = dummy
*     wsp(4)  = dummy
*     wsp(5)  = x_error, maximum among all local truncation errors
*     wsp(6)  = s_error, global sum of local truncation errors
*     wsp(7)  = tbrkdwn, if `happy breakdown', time when it occured
*     wsp(8)  = t_now, integration domain successfully covered
*
*----------------------------------------------------------------------|
*-----The following parameters may also be adjusted herein-------------|
*
      integer mxstep, mxreject, ideg
      double precision delta, gamma
      parameter( mxstep   = 500,
     .           mxreject = 0,
     .           ideg     = 6,
     .           delta    = 1.2d0,
     .           gamma    = 0.9d0 )

*     mxstep  : maximum allowable number of integration steps.
*               The value 0 means an infinite number of steps.
*
*     mxreject: maximum allowable number of rejections at each step.
*               The value 0 means an infinite number of rejections.
*
*     ideg    : the Pade approximation of type (ideg,ideg) is used as
*               an approximation to exp(H).
*
*     delta   : local truncation error `safety factor'
*
*     gamma   : stepsize `shrinking factor'
*
*----------------------------------------------------------------------|
*     Roger B. Sidje (rbs@maths.uq.edu.au)
*     EXPOKIT: Software Package for Computing Matrix Exponentials.
*     ACM - Transactions On Mathematical Software, 24(1):130-156, 1998
*----------------------------------------------------------------------|
*
      integer i, j, k1, mh, mx, iv, ih, j1v, ns, ifree, lfree, iexph,
     .        ireject,ibrkflag,mbrkdwn, nmult, nreject, nexph, nscale,
     .        nstep, iphih
      double precision sgn, t_out, tbrkdwn, step_min,step_max, err_loc,
     .                 s_error, x_error, t_now, t_new, t_step, t_old,
     .                 xm, beta, break_tol, p1, p2, p3, eps, rndoff,
     .                 avnorm, hj1j, hjj, SQR1

      intrinsic AINT,ABS,DBLE,LOG10,MAX,MIN,NINT,SIGN,SQRT
      double precision DDOT, DNRM2

*---  check restrictions on input parameters ...
      iflag = 0
      if ( lwsp.lt.n*(m+3)+5*(m+3)**2+ideg+1 ) iflag = -1
      if ( liwsp.lt.m+3 ) iflag = -2
      if ( m.ge.n .or. m.le.0 ) iflag = -3
      if ( iflag.ne.0 ) stop 'bad sizes (in input of DSPHIV)'
*
*---  initialisations ...
*
      k1 = 3
      mh = m + 3
      iv = 1
      ih = iv + n*(m+1) + n
      ifree = ih + mh*mh
      lfree = lwsp - ifree + 1

      ibrkflag = 0
      mbrkdwn  = m
      nmult    = 0
      nreject  = 0
      nexph    = 0
      nscale   = 0

      t_out    = ABS( t )
      tbrkdwn  = 0.0d0
      step_min = t_out
      step_max = 0.0d0
      nstep    = 0
      s_error  = 0.0d0
      x_error  = 0.0d0
      t_now    = 0.0d0
      t_new    = 0.0d0

      p1 = 4.0d0/3.0d0
 1    p2 = p1 - 1.0d0
      p3 = p2 + p2 + p2
      eps = ABS( p3-1.0d0 )
      if ( eps.eq.0.0d0 ) go to 1
      if ( tol.le.eps ) tol = SQRT( eps )
      rndoff = eps*anorm

      break_tol = 1.0d-7
*>>>  break_tol = tol
*>>>  break_tol = anorm*tol

*
*---  step-by-step integration ...
*
      sgn = SIGN( 1.0d0,t )
      SQR1 = SQRT( 0.1d0 )
      call DCOPY( n, v,1, w,1 )

 100  if ( t_now.ge.t_out ) goto 500

      nmult =  nmult + 1
      call matvec( w, wsp(iv) )
      call DAXPY( n, 1.0d0, u,1, wsp(iv),1 )
      beta = DNRM2( n, wsp(iv),1 )
      if ( beta.eq.0.0d0 ) goto 500
      call DSCAL( n, 1.0d0/beta, wsp(iv),1 )
      do i = 1,mh*mh
         wsp(ih+i-1) = 0.0d0
      end do

      if ( nstep.eq.0 ) then
*---     obtain the very first stepsize ...
         xm = 1.0d0/DBLE( m )
         p1 = tol*(((m+1)/2.72D0)**(m+1))*SQRT(2.0D0*3.14D0*(m+1))
         t_new = (1.0d0/anorm)*(p1/(4.0d0*beta*anorm))**xm
         p1 = 10.0d0**(NINT( LOG10( t_new )-SQR1 )-1)
         t_new = AINT( t_new/p1 + 0.55d0 ) * p1
      endif
      nstep = nstep + 1
      t_step = MIN( t_out-t_now, t_new )
*
*---  Lanczos loop ...
*
      j1v = iv + n
      do 200 j = 1,m
         nmult = nmult + 1
         call matvec( wsp(j1v-n), wsp(j1v) )
         if ( j.gt.1 )
     .     call DAXPY(n,-wsp(ih+(j-1)*mh+j-2),wsp(j1v-2*n),1,wsp(j1v),1)
         hjj = DDOT( n, wsp(j1v-n),1, wsp(j1v),1 )
         call DAXPY( n, -hjj, wsp(j1v-n),1, wsp(j1v),1 )
         hj1j = DNRM2( n, wsp(j1v),1 )
         wsp(ih+(j-1)*(mh+1)) = hjj
*---     if `happy breakdown' go straightforward at the end ...
         if ( hj1j.le.break_tol ) then
            print*,'happy breakdown: mbrkdwn =',j,' h =',hj1j
            k1 = 0
            ibrkflag = 1
            mbrkdwn = j
            tbrkdwn = t_now
            t_step = t_out-t_now
            goto 300
         endif
         wsp(ih+(j-1)*mh+j) = hj1j
         wsp(ih+j*mh+j-1) = hj1j
         call DSCAL( n, 1.0d0/hj1j, wsp(j1v),1 )
         j1v = j1v + n
 200  continue
      nmult = nmult + 1
      call matvec( wsp(j1v-n), wsp(j1v) )
      avnorm = DNRM2( n, wsp(j1v),1 )
*
*---  set 1's for the 3-extended scheme ...
*
 300  continue
      wsp(ih+mh*mbrkdwn) = 1.0d0
      wsp(ih+m*mh+m-1)   = 0.0d0
      wsp(ih+(m-1)*mh+m) = 0.0d0
      do i = 1,k1-1
         wsp(ih+(m+i)*mh+m+i-1) = 1.0d0
      end do
*
*---  loop while ireject<mxreject until the tolerance is reached ...
*
      ireject = 0
 401  continue
*
*---  compute w = beta*t_step*V*phi(t_step*H)*e1 + w
*
      nexph = nexph + 1
      mx = mbrkdwn + MAX( 1,k1 )

*---  irreducible rational Pade approximation ...
      call DGPADM( ideg, mx, sgn*t_step, wsp(ih),mh,
     .              wsp(ifree),lfree, iwsp, iexph, ns, iflag )
      iexph = ifree + iexph - 1
      iphih = iexph + mbrkdwn*mx
      nscale = nscale + ns
      wsp(iphih+mbrkdwn)   = hj1j*wsp(iphih+mx+mbrkdwn-1)
      wsp(iphih+mbrkdwn+1) = hj1j*wsp(iphih+2*mx+mbrkdwn-1)

 402  continue
*
*---  error estimate ...
*
      if ( k1.eq.0 ) then
         err_loc = tol
      else
         p1 = ABS( wsp(iphih+m) )   * beta
         p2 = ABS( wsp(iphih+m+1) ) * beta * avnorm
         if ( p1.gt.10.0d0*p2 ) then
            err_loc = p2
            xm = 1.0d0/DBLE( m+1 )
         elseif ( p1.gt.p2 ) then
            err_loc = (p1*p2)/(p1-p2)
            xm = 1.0d0/DBLE( m+1 )
         else
            err_loc = p1
            xm = 1.0d0/DBLE( m )
         endif
      endif
*
*---  reject the step-size if the error is not acceptable ...
*
      if ( (k1.ne.0) .and. (err_loc.gt.delta*t_step*tol) .and.
     .     (mxreject.eq.0 .or. ireject.lt.mxreject) ) then
         t_old = t_step
         t_step = gamma * t_step * (t_step*tol/err_loc)**xm
         p1 = 10.0d0**(NINT( LOG10( t_step )-SQR1 )-1)
         t_step = AINT( t_step/p1 + 0.55d0 ) * p1
         if ( itrace.ne.0 ) then
            print*,'t_step =',t_old
            print*,'err_loc =',err_loc
            print*,'err_required =',delta*t_old*tol
            print*,'stepsize rejected, stepping down to:',t_step
         endif
         ireject = ireject + 1
         nreject = nreject + 1
         if ( mxreject.ne.0 .and. ireject.gt.mxreject ) then
            print*,"Failure in DSPHIV: ---"
            print*,"The requested tolerance is too high."
            Print*,"Rerun with a smaller value."
            iflag = 2
            return
         endif
         goto 401
      endif
*
      mx = mbrkdwn + MAX( 0,k1-2 )
      call DGEMV( 'n', n,mx,beta,wsp(iv),n,wsp(iphih),1,1.0d0,w,1 )
*
*---  suggested value for the next stepsize ...
*
      t_new = gamma * t_step * (t_step*tol/err_loc)**xm
      p1 = 10.0d0**(NINT( LOG10( t_new )-SQR1 )-1)
      t_new = AINT( t_new/p1 + 0.55d0 ) * p1

      err_loc = MAX( err_loc,rndoff )
*
*---  update the time covered ...
*
      t_now = t_now + t_step
*
*---  display and keep some information ...
*
      if ( itrace.ne.0 ) then
         print*,'integration',nstep,'---------------------------------'
         print*,'scale-square =',ns
         print*,'step_size =',t_step
         print*,'err_loc   =',err_loc
         print*,'next_step =',t_new
      endif

      step_min = MIN( step_min, t_step )
      step_max = MAX( step_max, t_step )
      s_error = s_error + err_loc
      x_error = MAX( x_error, err_loc )

      if ( mxstep.eq.0 .or. nstep.lt.mxstep ) goto 100
      iflag = 1

 500  continue

      iwsp(1) = nmult
      iwsp(2) = nexph
      iwsp(3) = nscale
      iwsp(4) = nstep
      iwsp(5) = nreject
      iwsp(6) = ibrkflag
      iwsp(7) = mbrkdwn

      wsp(1)  = step_min
      wsp(2)  = step_max
      wsp(3)  = 0.0d0
      wsp(4)  = 0.0d0
      wsp(5)  = x_error
      wsp(6)  = s_error
      wsp(7)  = tbrkdwn
      wsp(8)  = sgn*t_now
      END
      subroutine dswap (n,dx,incx,dy,incy)

c*********************************************************************72
c
cc DSWAP interchanges two vectors.
c
c     interchanges two vectors.
c     uses unrolled loops for increments equal one.
c     jack dongarra, linpack, 3/11/78.
c
      implicit none

      double precision dx(1),dy(1),dtemp
      integer i,incx,incy,ix,iy,m,mp1,n
c
      if(n.le.0)return
      if(incx.eq.1.and.incy.eq.1)go to 20
c
c       code for unequal increments or equal increments not equal
c         to 1
c
      ix = 1
      iy = 1
      if(incx.lt.0)ix = (-n+1)*incx + 1
      if(incy.lt.0)iy = (-n+1)*incy + 1
      do 10 i = 1,n
        dtemp = dx(ix)
        dx(ix) = dy(iy)
        dy(iy) = dtemp
        ix = ix + incx
        iy = iy + incy
   10 continue
      return
c
c       code for both increments equal to 1
c
c
c       clean-up loop
c
   20 m = mod(n,3)
      if( m .eq. 0 ) go to 40
      do 30 i = 1,m
        dtemp = dx(i)
        dx(i) = dy(i)
        dy(i) = dtemp
   30 continue
      if( n .lt. 3 ) return
   40 mp1 = m + 1
      do 50 i = mp1,n,3
        dtemp = dx(i)
        dx(i) = dy(i)
        dy(i) = dtemp
        dtemp = dx(i + 1)
        dx(i + 1) = dy(i + 1)
        dy(i + 1) = dtemp
        dtemp = dx(i + 2)
        dx(i + 2) = dy(i + 2)
        dy(i + 2) = dtemp
   50 continue
      return
      end
      subroutine dsysv(uplo, n,m, a,lda, ipiv, b,ldb, wrk,lwrk, iflag )

c*********************************************************************72
c
cc DSYSV solves multiple symmetric linear systems of the form A*X=B.
c
      implicit none

      character UPLO*1
      integer N, M, LDA, LDB, LWRK, IFLAG, IPIV(N)
      double precision A(LDA,N), B(LDB,M), WRK(LWRK)
      integer j

      call DSIFA( A,LDA, N, IPIV, IFLAG )
      if ( IFLAG.ne.0 ) stop "Error in DSYSV (LDL' factorisation)"
      do j = 1,M
         call DSISL( A,LDA, N, IPIV, B(1,j) )
      end do
      end
      function dznrm2( n, zx, incx)

c*********************************************************************72
c
cc DZNRM2 computes the unitary norm of a complex vector.
c
      implicit none

      double precision dznrm2
      logical imag, scale
      integer i, incx, ix, n, next
      double precision cutlo, cuthi, hitest, sum, xmax, absx, zero, one
      double complex      zx(1)
      double precision dreal,dimag
      double complex zdumr,zdumi
      dreal(zdumr) = zdumr
      dimag(zdumi) = (0.0d0,-1.0d0)*zdumi
      data         zero, one /0.0d0, 1.0d0/
c
c     unitary norm of the complex n-vector stored in zx() with storage
c     increment incx .
c     if    n .le. 0 return with result = 0.
c     if n .ge. 1 then incx must be .ge. 1
c
c           c.l.lawson , 1978 jan 08
c     modified 3/93 to return if incx .le. 0.
c
c     four phase method     using two built-in constants that are
c     hopefully applicable to all machines.
c         cutlo = maximum of  sqrt(u/eps)  over all known machines.
c         cuthi = minimum of  sqrt(v)      over all known machines.
c     where
c         eps = smallest no. such that eps + 1. .gt. 1.
c         u   = smallest positive no.   (underflow limit)
c         v   = largest  no.            (overflow  limit)
c
c     brief outline of algorithm..
c
c     phase 1    scans zero components.
c     move to phase 2 when a component is nonzero and .le. cutlo
c     move to phase 3 when a component is .gt. cutlo
c     move to phase 4 when a component is .ge. cuthi/m
c     where m = n for x() real and m = 2*n for complex.
c
c     values for cutlo and cuthi..
c     from the environmental parameters listed in the imsl converter
c     document the limiting values are as follows..
c     cutlo, s.p.   u/eps = 2**(-102) for  honeywell.  close seconds are
c                   univac and dec at 2**(-103)
c                   thus cutlo = 2**(-51) = 4.44089e-16
c     cuthi, s.p.   v = 2**127 for univac, honeywell, and dec.
c                   thus cuthi = 2**(63.5) = 1.30438e19
c     cutlo, d.p.   u/eps = 2**(-67) for honeywell and dec.
c                   thus cutlo = 2**(-33.5) = 8.23181d-11
c     cuthi, d.p.   same as s.p.  cuthi = 1.30438d19
c     data cutlo, cuthi / 8.232d-11,  1.304d19 /
c     data cutlo, cuthi / 4.441e-16,  1.304e19 /
      data cutlo, cuthi / 8.232d-11,  1.304d19 /
c
      if(n .gt. 0 .and. incx.gt.0) go to 10
         dznrm2  = zero
         go to 300
c
   10 assign 30 to next
      sum = zero
      i = 1
c                                                 begin main loop
      do 220 ix = 1,n
         absx = dabs(dreal(zx(i)))
         imag = .false.
         go to next,(30, 50, 70, 90, 110)
   30 if( absx .gt. cutlo) go to 85
      assign 50 to next
      scale = .false.
c
c                        phase 1.  sum is zero
c
   50 if( absx .eq. zero) go to 200
      if( absx .gt. cutlo) go to 85
c
c                                prepare for phase 2.
      assign 70 to next
      go to 105
c
c                                prepare for phase 4.
c
  100 assign 110 to next
      sum = (sum / absx) / absx
  105 scale = .true.
      xmax = absx
      go to 115
c
c                   phase 2.  sum is small.
c                             scale to avoid destructive underflow.
c
   70 if( absx .gt. cutlo ) go to 75
c
c                     common code for phases 2 and 4.
c                     in phase 4 sum is large.  scale to avoid overflow.
c
  110 if( absx .le. xmax ) go to 115
         sum = one + sum * (xmax / absx)**2
         xmax = absx
         go to 200
c
  115 sum = sum + (absx/xmax)**2
      go to 200
c
c
c                  prepare for phase 3.
c
   75 sum = (sum * xmax) * xmax
c
   85 assign 90 to next
      scale = .false.
c
c     for real or d.p. set hitest = cuthi/n
c     for complex      set hitest = cuthi/(2*n)
c
      hitest = cuthi/dble( 2*n )
c
c                   phase 3.  sum is mid-range.  no scaling.
c
   90 if(absx .ge. hitest) go to 100
         sum = sum + absx**2
  200 continue
c                  control selection of real and imaginary parts.
c
      if(imag) go to 210
         absx = dabs(dimag(zx(i)))
         imag = .true.
      go to next,(  50, 70, 90, 110 )
c
  210 continue
      i = i + incx
  220 continue
c
c              end of main loop.
c              compute square root and adjust for scaling.
c
      dznrm2 = dsqrt(sum)
      if(scale) dznrm2 = dznrm2 * xmax
  300 continue
      return
      end
      subroutine getpat ( filename, n, nz, ia, ja )

c*********************************************************************72
c
cc GETPAT loads a Harwell-Boeing pattern.
c
      implicit none

      character*(*) filename
      integer          n, nz
      integer          ia(*), ja(*)
c
      character        title*72, key*8, type*3, ptrfmt*16,
     .                 indfmt*16, valfmt*20, rhsfmt*20, rhstyp*1
      integer          totcrd, ptrcrd, indcrd, valcrd, rhscrd, nrow,
     .                 nrhs, nrhsix, i, io, nn, nnz
c
      open( UNIT=7, STATUS='old', IOstat=io, FILE=filename )
      if ( io.ne.0 ) stop 'Could not access Harwell-Boeing matrix'
      read( UNIT=7,FMT=10 ) title, key,
     .     totcrd, ptrcrd, indcrd, valcrd, rhscrd,
     .     type, nrow, nn, nnz, nrhs,
     .     ptrfmt, indfmt, valfmt, rhsfmt
      print*, title,'type :',type,' size :',nrow,nn
      print*,'order :',nn,' number of nonzero :',nnz

      if ( nn.gt.n ) stop 'Please increase nmax.'
      if ( nnz.gt.nz ) stop 'Please increase nzmax.'
      n = nn
      nz = nnz
      if ( rhscrd.gt.0 ) then
         read( UNIT=7,FMT=11 ) rhstyp, nrhs, nrhsix
         print*,'There is a second hand'
      endif

c  read data...
      read( UNIT=7,FMT=ptrfmt ) (ja(i), i = 1,n+1)
      read( UNIT=7,FMT=indfmt ) (ia(i), i = 1,nz)

      close( UNIT=7 )
      print*,'Harwell-Boeing pattern loaded'
10    format(A72, A8/ 5i14 / A3, 11x, 4i14 / 2a16, 2a20)
11    format(A1, 13x, 2i14)
      end
      function idamax(n,dx,incx)

c*********************************************************************72
c
cc IDAMAX indexes the vector entry of maximum absolute value.
c
c     jack dongarra, linpack, 3/11/78.
c     modified 3/93 to return if incx .le. 0.
c
      implicit none

      integer idamax
      double precision dx(1),dmax
      integer i,incx,ix,n
c
      idamax = 0
      if( n.lt.1 .or. incx.le.0 ) return
      idamax = 1
      if(n.eq.1)return
      if(incx.eq.1)go to 20
c
c        code for increment not equal to 1
c
      ix = 1
      dmax = dabs(dx(1))
      ix = ix + incx
      do 10 i = 2,n
         if(dabs(dx(ix)).le.dmax) go to 5
         idamax = i
         dmax = dabs(dx(ix))
    5    ix = ix + incx
   10 continue
      return
c
c        code for increment equal to 1
c
   20 dmax = dabs(dx(1))
      do 30 i = 2,n
         if(dabs(dx(i)).le.dmax) go to 30
         idamax = i
         dmax = dabs(dx(i))
   30 continue
      return
      end
      subroutine idsrt1( nx, ix, xx )

c*********************************************************************72
c
cc IDSRT1 is an indirect sort; sorts IX and carries XX along.
c
*---  adapted from a SLAP (Sparse Linear Algebra Package) code.
c

      implicit none
      integer          nx, ix(nx)
      double precision xx(nx)

      integer          M,I,J,K,IL(21),IU(21), IT,IIT,IJ,L
      double precision TX, TTX, R

      if ( nx.le.1 ) return

*---  And now...Just a little black magic...
      M = 1
      I = 1
      J = NX
      R = .375
 210  IF( R.LE.0.5898437 ) THEN
         R = R + 3.90625E-2
      ELSE
         R = R-.21875
      ENDIF
 225  K = I
*
*---  Select a central element of the array and save it in location
*---  IT, TX.
*
      IJ = I + IDINT( DBLE(J-I)*R )
      IT = IX(IJ)
      TX = XX(IJ)
*
*---  If first element of array is greater than IT, interchange with IT.
*
      IF( IX(I).GT.IT ) THEN
         IX(IJ) = IX(I)
         IX(I)  = IT
         IT     = IX(IJ)
         XX(IJ)  = XX(I)
         XX(I)   = TX
         TX     = XX(IJ)
      ENDIF
      L=J
*
*---  If last element of array is less than IT, swap with IT.
*
      IF( IX(J).LT.IT ) THEN
         IX(IJ) = IX(J)
         IX(J)  = IT
         IT     = IX(IJ)
         XX(IJ)  = XX(J)
         XX(J)   = TX
         TX     = XX(IJ)
*
*---  If first element of array is greater than IT, swap with IT.
*
         IF ( IX(I).GT.IT ) THEN
            IX(IJ) = IX(I)
            IX(I)  = IT
            IT     = IX(IJ)
            XX(IJ)  = XX(I)
            XX(I)   = TX
            TX     = XX(IJ)
         ENDIF
      ENDIF
*
*---  Find an element in the second half of the array which is
*---  smaller than IT.
*
 240  L=L-1
      IF( IX(L).GT.IT ) GO TO 240
*
*---  Find an element in the first half of the array which is
*---  greater than IT.
*
 245  K=K+1
      IF( IX(K).LT.IT ) GO TO 245
*
*---  Interchange these elements.
*
      IF( K.LE.L ) THEN
         IIT   = IX(L)
         IX(L) = IX(K)
         IX(K) = IIT
         TTX   = XX(L)
         XX(L)  = XX(K)
         XX(K)  = TTX
         GOTO 240
      ENDIF
*
*---  Save upper and lower subscripts of the array yet to be sorted.
*
      IF( L-I.GT.J-K ) THEN
         IL(M) = I
         IU(M) = L
         I = K
         M = M+1
      ELSE
         IL(M) = K
         IU(M) = J
         J = L
         M = M+1
      ENDIF
      GO TO 260
*
*---  Begin again on another portion of the unsorted array.
*
 255  M = M-1
      IF( M.EQ.0 ) GO TO 300
      I = IL(M)
      J = IU(M)
 260  IF( J-I.GE.1 ) GO TO 225
      IF( I.EQ.J ) GO TO 255
      IF( I.EQ.1 ) GO TO 210
      I = I-1
 265  I = I+1
      IF( I.EQ.J ) GO TO 255
      IT = IX(I+1)
      TX =  XX(I+1)
      IF( IX(I).LE.IT ) GO TO 265
      K=I
 270  IX(K+1) = IX(K)
      XX(K+1)  =  XX(K)
      K = K-1
      IF( IT.LT.IX(K) ) GO TO 270
      IX(K+1) = IT
      XX(K+1)  = TX
      GO TO 265

 300  CONTINUE
      RETURN
      END
      subroutine idsrt2( nx, ix, ixx, xx )

c*********************************************************************72
c
cc IDSRT2 is an indirect sort: sorts IX and carries IXX and XX along.
c
*---  adapted from a SLAP (Sparse Linear Algebra Package) code.
c
      implicit none

      integer          nx, ix(nx), ixx(nx)
      double precision xx(nx)

      integer          M,I,J,K,IL(21),IU(21), IT,IIT,IJ,JT,JJT,L
      double precision TX, TTX, R

      if ( nx.le.1 ) return

*---  And now...Just a little black magic...
      M = 1
      I = 1
      J = NX
      R = .375
 210  IF( R.LE.0.5898437 ) THEN
         R = R + 3.90625E-2
      ELSE
         R = R-.21875
      ENDIF
 225  K = I
*
*---  Select a central element of the array and save it in location
*---  IT, JT, TX.
*
      IJ = I + IDINT( DBLE(J-I)*R )
      IT = IX(IJ)
      JT = IXX(IJ)
      TX = XX(IJ)
*
*---  If first element of array is greater than IT, interchange with IT.
*
      IF( IX(I).GT.IT ) THEN
         IX(IJ) = IX(I)
         IX(I)  = IT
         IT     = IX(IJ)
         IXX(IJ)= IXX(I)
         IXX(I) = JT
         JT     = IXX(IJ)
         XX(IJ) = XX(I)
         XX(I)  = TX
         TX     = XX(IJ)
      ENDIF
      L=J
*
*---  If last element of array is less than IT, swap with IT.
*
      IF( IX(J).LT.IT ) THEN
         IX(IJ) = IX(J)
         IX(J)  = IT
         IT     = IX(IJ)
         IXX(IJ)= IXX(J)
         IXX(J) = JT
         JT     = IXX(IJ)
         XX(IJ) = XX(J)
         XX(J)  = TX
         TX     = XX(IJ)
*
*---  If first element of array is greater than IT, swap with IT.
*
         IF ( IX(I).GT.IT ) THEN
            IX(IJ) = IX(I)
            IX(I)  = IT
            IT     = IX(IJ)
            IXX(IJ)= IXX(I)
            IXX(I) = JT
            JT     = IXX(IJ)
            XX(IJ) = XX(I)
            XX(I)  = TX
            TX     = XX(IJ)
         ENDIF
      ENDIF
*
*---  Find an element in the second half of the array which is
*---  smaller than IT.
*
 240  L=L-1
      IF( IX(L).GT.IT ) GO TO 240
*
*---  Find an element in the first half of the array which is
*---  greater than IT.
*
 245  K=K+1
      IF( IX(K).LT.IT ) GO TO 245
*
*---  Interchange these elements.
*
      IF( K.LE.L ) THEN
         IIT   = IX(L)
         IX(L) = IX(K)
         IX(K) = IIT
         JJT   = IXX(L)
         IXX(L)= IXX(K)
         IXX(K)= JJT
         TTX   = XX(L)
         XX(L) = XX(K)
         XX(K) = TTX
         GOTO 240
      ENDIF
*
*---  Save upper and lower subscripts of the array yet to be sorted.
*
      IF( L-I.GT.J-K ) THEN
         IL(M) = I
         IU(M) = L
         I = K
         M = M+1
      ELSE
         IL(M) = K
         IU(M) = J
         J = L
         M = M+1
      ENDIF
      GO TO 260
*
*---  Begin again on another portion of the unsorted array.
*
 255  M = M-1
      IF( M.EQ.0 ) GO TO 300
      I = IL(M)
      J = IU(M)
 260  IF( J-I.GE.1 ) GO TO 225
      IF( I.EQ.J ) GO TO 255
      IF( I.EQ.1 ) GO TO 210
      I = I-1
 265  I = I+1
      IF( I.EQ.J ) GO TO 255
      IT = IX(I+1)
      JT = IXX(I+1)
      TX =  XX(I+1)
      IF( IX(I).LE.IT ) GO TO 265
      K=I
 270  IX(K+1) = IX(K)
      IXX(K+1) = IXX(K)
      XX(K+1)  =  XX(K)
      K = K-1
      IF( IT.LT.IX(K) ) GO TO 270
      IX(K+1) = IT
      IXX(K+1) = JT
      XX(K+1)  = TX
      GO TO 265

 300  CONTINUE
      RETURN
      END
      function izamax(n,zx,incx)

c*********************************************************************72
c
cc IZAMAX indexes the complex vector entry of maximum magnitude.
c
c     finds the index of element having max. absolute value.
c     jack dongarra, 1/15/85.
c     modified 3/93 to return if incx .le. 0.
c
      implicit none

      integer izamax
      double complex zx(1)
      double precision smax
      integer i,incx,ix,n
      double precision dcabs1
c
      izamax = 0
      if( n.lt.1 .or. incx.le.0 )return
      izamax = 1
      if(n.eq.1)return
      if(incx.eq.1)go to 20
c
c        code for increment not equal to 1
c
      ix = 1
      smax = dcabs1(zx(1))
      ix = ix + incx
      do 10 i = 2,n
         if(dcabs1(zx(ix)).le.smax) go to 5
         izamax = i
         smax = dcabs1(zx(ix))
    5    ix = ix + incx
   10 continue
      return
c
c        code for increment equal to 1
c
   20 smax = dcabs1(zx(1))
      do 30 i = 2,n
         if(dcabs1(zx(i)).le.smax) go to 30
         izamax = i
         smax = dcabs1(zx(i))
   30 continue
      return
      end
      subroutine izsrt1( nx, ix, xx )

c*********************************************************************72
c
cc IZSRT1 is an indirect sort; sorts IX and carries XX along.
c
*---  adapted from a SLAP (Sparse Linear Algebra Package) code.
c
      implicit none

      integer          nx, ix(nx)
      complex*16       xx(nx)

      integer          M,I,J,K,IL(21),IU(21), IT,IIT,IJ,L
      complex*16       TX, TTX
      double precision R

      if ( nx.le.1 ) return

*---  And now...Just a little black magic...
      M = 1
      I = 1
      J = NX
      R = .375
 210  IF( R.LE.0.5898437 ) THEN
         R = R + 3.90625E-2
      ELSE
         R = R-.21875
      ENDIF
 225  K = I
*
*---  Select a central element of the array and save it in location
*---  IT, TX.
*
      IJ = I + IDINT( DBLE(J-I)*R )
      IT = IX(IJ)
      TX = XX(IJ)
*
*---  If first element of array is greater than IT, interchange with IT.
*
      IF( IX(I).GT.IT ) THEN
         IX(IJ) = IX(I)
         IX(I)  = IT
         IT     = IX(IJ)
         XX(IJ)  = XX(I)
         XX(I)   = TX
         TX     = XX(IJ)
      ENDIF
      L=J
*
*---  If last element of array is less than IT, swap with IT.
*
      IF( IX(J).LT.IT ) THEN
         IX(IJ) = IX(J)
         IX(J)  = IT
         IT     = IX(IJ)
         XX(IJ)  = XX(J)
         XX(J)   = TX
         TX     = XX(IJ)
*
*---  If first element of array is greater than IT, swap with IT.
*
         IF ( IX(I).GT.IT ) THEN
            IX(IJ) = IX(I)
            IX(I)  = IT
            IT     = IX(IJ)
            XX(IJ)  = XX(I)
            XX(I)   = TX
            TX     = XX(IJ)
         ENDIF
      ENDIF
*
*---  Find an element in the second half of the array which is
*---  smaller than IT.
*
 240  L=L-1
      IF( IX(L).GT.IT ) GO TO 240
*
*---  Find an element in the first half of the array which is
*---  greater than IT.
*
 245  K=K+1
      IF( IX(K).LT.IT ) GO TO 245
*
*---  Interchange these elements.
*
      IF( K.LE.L ) THEN
         IIT   = IX(L)
         IX(L) = IX(K)
         IX(K) = IIT
         TTX   = XX(L)
         XX(L)  = XX(K)
         XX(K)  = TTX
         GOTO 240
      ENDIF
*
*---  Save upper and lower subscripts of the array yet to be sorted.
*
      IF( L-I.GT.J-K ) THEN
         IL(M) = I
         IU(M) = L
         I = K
         M = M+1
      ELSE
         IL(M) = K
         IU(M) = J
         J = L
         M = M+1
      ENDIF
      GO TO 260
*
*---  Begin again on another portion of the unsorted array.
*
 255  M = M-1
      IF( M.EQ.0 ) GO TO 300
      I = IL(M)
      J = IU(M)
 260  IF( J-I.GE.1 ) GO TO 225
      IF( I.EQ.J ) GO TO 255
      IF( I.EQ.1 ) GO TO 210
      I = I-1
 265  I = I+1
      IF( I.EQ.J ) GO TO 255
      IT = IX(I+1)
      TX =  XX(I+1)
      IF( IX(I).LE.IT ) GO TO 265
      K=I
 270  IX(K+1) = IX(K)
      XX(K+1)  =  XX(K)
      K = K-1
      IF( IT.LT.IX(K) ) GO TO 270
      IX(K+1) = IT
      XX(K+1)  = TX
      GO TO 265

 300  CONTINUE
      RETURN
      END
      subroutine izsrt2( nx, ix, ixx, xx )

c*********************************************************************72
c
cc IZSRT2 is an indirect sort: sorts IX and carries IXX and XX along.
c
*---  adapted from a SLAP (Sparse Linear Algebra Package) code.
c
      implicit none

      integer          nx, ix(nx), ixx(nx)
      complex*16       xx(nx)

      integer          M,I,J,K,IL(21),IU(21), IT,IIT,IJ,JT,JJT,L
      complex*16       TX, TTX
      double precision R

      if ( nx.le.1 ) return

*---  And now...Just a little black magic...
      M = 1
      I = 1
      J = NX
      R = .375
 210  IF( R.LE.0.5898437 ) THEN
         R = R + 3.90625E-2
      ELSE
         R = R-.21875
      ENDIF
 225  K = I
*
*---  Select a central element of the array and save it in location
*---  IT, JT, TX.
*
      IJ = I + IDINT( DBLE(J-I)*R )
      IT = IX(IJ)
      JT = IXX(IJ)
      TX = XX(IJ)
*
*---  If first element of array is greater than IT, interchange with IT.
*
      IF( IX(I).GT.IT ) THEN
         IX(IJ) = IX(I)
         IX(I)  = IT
         IT     = IX(IJ)
         IXX(IJ)= IXX(I)
         IXX(I) = JT
         JT     = IXX(IJ)
         XX(IJ) = XX(I)
         XX(I)  = TX
         TX     = XX(IJ)
      ENDIF
      L=J
*
*---  If last element of array is less than IT, swap with IT.
*
      IF( IX(J).LT.IT ) THEN
         IX(IJ) = IX(J)
         IX(J)  = IT
         IT     = IX(IJ)
         IXX(IJ)= IXX(J)
         IXX(J) = JT
         JT     = IXX(IJ)
         XX(IJ) = XX(J)
         XX(J)  = TX
         TX     = XX(IJ)
*
*---  If first element of array is greater than IT, swap with IT.
*
         IF ( IX(I).GT.IT ) THEN
            IX(IJ) = IX(I)
            IX(I)  = IT
            IT     = IX(IJ)
            IXX(IJ)= IXX(I)
            IXX(I) = JT
            JT     = IXX(IJ)
            XX(IJ) = XX(I)
            XX(I)  = TX
            TX     = XX(IJ)
         ENDIF
      ENDIF
*
*---  Find an element in the second half of the array which is
*---  smaller than IT.
*
 240  L=L-1
      IF( IX(L).GT.IT ) GO TO 240
*
*---  Find an element in the first half of the array which is
*---  greater than IT.
*
 245  K=K+1
      IF( IX(K).LT.IT ) GO TO 245
*
*---  Interchange these elements.
*
      IF( K.LE.L ) THEN
         IIT   = IX(L)
         IX(L) = IX(K)
         IX(K) = IIT
         JJT   = IXX(L)
         IXX(L)= IXX(K)
         IXX(K)= JJT
         TTX   = XX(L)
         XX(L) = XX(K)
         XX(K) = TTX
         GOTO 240
      ENDIF
*
*---  Save upper and lower subscripts of the array yet to be sorted.
*
      IF( L-I.GT.J-K ) THEN
         IL(M) = I
         IU(M) = L
         I = K
         M = M+1
      ELSE
         IL(M) = K
         IU(M) = J
         J = L
         M = M+1
      ENDIF
      GO TO 260
*
*---  Begin again on another portion of the unsorted array.
*
 255  M = M-1
      IF( M.EQ.0 ) GO TO 300
      I = IL(M)
      J = IU(M)
 260  IF( J-I.GE.1 ) GO TO 225
      IF( I.EQ.J ) GO TO 255
      IF( I.EQ.1 ) GO TO 210
      I = I-1
 265  I = I+1
      IF( I.EQ.J ) GO TO 255
      IT = IX(I+1)
      JT = IXX(I+1)
      TX =  XX(I+1)
      IF( IX(I).LE.IT ) GO TO 265
      K=I
 270  IX(K+1) = IX(K)
      IXX(K+1) = IXX(K)
      XX(K+1)  =  XX(K)
      K = K-1
      IF( IT.LT.IX(K) ) GO TO 270
      IX(K+1) = IT
      IXX(K+1) = JT
      XX(K+1)  = TX
      GO TO 265

 300  CONTINUE
      RETURN
      END
      subroutine loadhb ( filename, spformat, n, nz, ia, ja, a, iwsp )

c*********************************************************************72
c
cc LOADHB loads a Harwell-Boeing sparse matrix from a file.
c
      implicit none

      character*(*) filename
      character *3 spformat*3
      integer          n, nz, ia(*), ja(*), iwsp(*)
      double precision a(*)

*---  Purpose ---------------------------------------------------------|
*
*---  LOADHB loads a matrix stored under the Harwell-Boeing format
*     and renders it into the sparse format specified by spformat.
*
*---  Arguments -------------------------------------------------------|
*
*     filename (input)
*           name of the file containing the matrix.
*
*     spformat (input)
*           sparse format in which the matrix is forced to be on output
*           accepted values are:
*              'COO' : COOrdinates
*              'CRS' : Compressed Row Storage
*              'CCS' : Compressed Column Storage (default H-B format)
*
*     n (input/output)
*           On input,  the maximum allowable order
*           On output, the actual order of the matrix loaded
*
*     nz (input/output)
*           On input,  the maximum allowable number of non zero entries
*           On output, the actual number of non zero entries
*
*     ia,ja,a (output)
*           sparse matrix data stored in the format given in spformat
*           sufficient room is needed to achieve this: each component
*           must be of length >= nz. If the matrix is symmetric, both
*           lower and upper parts are included explicitly
*
*     iwsp (workspace) of length >= n
*
*----------------------------------------------------------------------|
*
      character        title*72, key*8, type*3, ptrfmt*16,
     .                 indfmt*16, valfmt*20, rhsfmt*20, rhstyp*1
      integer          totcrd, ptrcrd, indcrd, valcrd, rhscrd, nrow,
     .                 nrhs, nrhsix, i, j, k, io, nn, nnz
      intrinsic        INDEX


      open( UNIT=7, STATUS='old', IOstat=io, FILE=filename )

      if ( io.ne.0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'LOADHB - Fatal error!'
        write ( *, '(a)' ) '  Could not access Harwell-Boeing matrix'
        write ( *, '(a)' ) '  "' // trim ( filename ) // '".'
        stop
      end if

      read( UNIT=7,FMT=10 ) title, key,
     .     totcrd, ptrcrd, indcrd, valcrd, rhscrd,
     .     type, nrow, nn, nnz, nrhs,
     .     ptrfmt, indfmt, valfmt, rhsfmt
      print*, title,'type :',type,' size :',nrow,nn
      print*,'order :',nn,' number of nonzero :',nnz

      if ( nn.gt.n ) stop 'in LOADHB. Please increase n'
      if ( nnz.gt.nz ) stop 'in LOADHB. Please increase nz'

*---  leave if there is no values ...
      if ( valcrd.le.0 ) stop 'Empty Harwell-Boeing matrix'
      if ( rhscrd.gt.0 ) then
         read( UNIT=7,FMT=11 ) rhstyp, nrhs, nrhsix
         print*,'There is a second hand'
      endif

*---  read data...
      read( UNIT=7,FMT=ptrfmt ) (ja(i), i = 1,nn+1)
      read( UNIT=7,FMT=indfmt ) (ia(i), i = 1,nnz)
      read( UNIT=7,FMT=valfmt ) (a(i),  i = 1,nnz)
      close( UNIT=7 )

*---  for the sake of experiments, store both parts if symmetric matrix
      if ( type.eq.'RSA' ) then
*---     expand ja indices ...
         k = nnz
         do j = nn,1,-1
            do i = 1,ja(j+1)-ja(j)
               ja(k) = j
               k = k-1
            end do
         end do
*---     insert the other half ...
         k = nnz
         do i = 1,k
            if ( ia(i).ne.ja(i) ) then
               nnz = nnz + 1
               if ( nnz.gt.nz ) stop 'in LOADHB. Please increase nz'
               ia(nnz) = ja(i)
               ja(nnz) = ia(i)
               a(nnz) = a(i)
            endif
         end do
         type = 'COO'
      else
         type = 'CCS'
      endif
      call dgcnvr( type,spformat,'n', nn,nn, nnz, ia, ja, a, iwsp )
      n = nn
      nz = nnz

      print*,'Harwell-Boeing matrix loaded'
10    format(A72, A8/ 5i14 / A3, 11x, 4i14 / 2a16, 2a20)
11    format(A1, 13x, 2i14)
      end
      function lsame( ca, cb )

c*********************************************************************72
c
cc LSAME compares two characters for identity.
*
*  -- LAPACK auxiliary routine (version 1.1) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     February 29, 1992
*
      implicit none

      logical lsame
      CHARACTER          CA, CB
*     ..
*
*  Purpose
*  =======
*
*  LSAME returns .TRUE. if CA is the same letter as CB regardless of
*  case.
*
*  Arguments
*  =========
*
*  CA      (input) CHARACTER*1
*  CB      (input) CHARACTER*1
*          CA and CB specify the single characters to be compared.
*
*     .. Intrinsic Functions ..
      INTRINSIC          ICHAR
*     ..
*     .. Local Scalars ..
      INTEGER            INTA, INTB, ZCODE
*     ..
*     .. Executable Statements ..
*
*     Test if the characters are equal
*
      LSAME = CA.EQ.CB
      IF( LSAME )
     $   RETURN
*
*     Now test for equivalence if both characters are alphabetic.
*
      ZCODE = ICHAR( 'Z' )
*
*     Use 'Z' rather than 'A' so that ASCII can be detected on Prime
*     machines, on which ICHAR returns a value with bit 8 set.
*     ICHAR('A') on Prime machines returns 193 which is the same as
*     ICHAR('A') on an EBCDIC machine.
*
      INTA = ICHAR( CA )
      INTB = ICHAR( CB )
*
      IF( ZCODE.EQ.90 .OR. ZCODE.EQ.122 ) THEN
*
*        ASCII is assumed - ZCODE is the ASCII code of either lower or
*        upper case 'Z'.
*
         IF( INTA.GE.97 .AND. INTA.LE.122 ) INTA = INTA - 32
         IF( INTB.GE.97 .AND. INTB.LE.122 ) INTB = INTB - 32
*
      ELSE IF( ZCODE.EQ.233 .OR. ZCODE.EQ.169 ) THEN
*
*        EBCDIC is assumed - ZCODE is the EBCDIC code of either lower or
*        upper case 'Z'.
*
         IF( INTA.GE.129 .AND. INTA.LE.137 .OR.
     $       INTA.GE.145 .AND. INTA.LE.153 .OR.
     $       INTA.GE.162 .AND. INTA.LE.169 ) INTA = INTA + 64
         IF( INTB.GE.129 .AND. INTB.LE.137 .OR.
     $       INTB.GE.145 .AND. INTB.LE.153 .OR.
     $       INTB.GE.162 .AND. INTB.LE.169 ) INTB = INTB + 64
*
      ELSE IF( ZCODE.EQ.218 .OR. ZCODE.EQ.250 ) THEN
*
*        ASCII is assumed, on Prime machines - ZCODE is the ASCII code
*        plus 128 of either lower or upper case 'Z'.
*
         IF( INTA.GE.225 .AND. INTA.LE.250 ) INTA = INTA - 32
         IF( INTB.GE.225 .AND. INTB.LE.250 ) INTB = INTB - 32
      END IF
      LSAME = INTA.EQ.INTB
*
*     RETURN
*
*     End of LSAME
*
      END
      subroutine timestamp ( )

c*********************************************************************72
c
cc TIMESTAMP prints out the current YMDHMS date as a timestamp.
c
c  Discussion:
c
c    This FORTRAN77 version is made available for cases where the
c    FORTRAN90 version cannot be used.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    12 January 2007
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    None
c
      implicit none

      character * ( 8 ) ampm
      integer d
      character * ( 8 ) date
      integer h
      integer m
      integer mm
      character * ( 9 ) month(12)
      integer n
      integer s
      character * ( 10 ) time
      integer y

      save month

      data month /
     &  'January  ', 'February ', 'March    ', 'April    ',
     &  'May      ', 'June     ', 'July     ', 'August   ',
     &  'September', 'October  ', 'November ', 'December ' /

      call date_and_time ( date, time )

      read ( date, '(i4,i2,i2)' ) y, m, d
      read ( time, '(i2,i2,i2,1x,i3)' ) h, n, s, mm

      if ( h .lt. 12 ) then
        ampm = 'AM'
      else if ( h .eq. 12 ) then
        if ( n .eq. 0 .and. s .eq. 0 ) then
          ampm = 'Noon'
        else
          ampm = 'PM'
        end if
      else
        h = h - 12
        if ( h .lt. 12 ) then
          ampm = 'PM'
        else if ( h .eq. 12 ) then
          if ( n .eq. 0 .and. s .eq. 0 ) then
            ampm = 'Midnight'
          else
            ampm = 'AM'
          end if
        end if
      end if

      write ( *,
     &  '(i2,1x,a,1x,i4,2x,i2,a1,i2.2,a1,i2.2,a1,i3.3,1x,a)' )
     &  d, month(m), y, h, ':', n, ':', s, '.', mm, ampm

      return
      end
      subroutine tnspos ( n, nz, ia, ja, a, iwsp )

c*********************************************************************72
c
cc TNSPOS transposes a CRS matrix.
c
      implicit none

      integer          n, nz, ia(nz), ja(nz), iwsp(n)
      double precision a(nz)

      integer i, j

      call dgcnvr( 'crs','ccs','n', n,n, nz, ia, ja, a, iwsp )
      do i = 1,nz
         j = ja(i)
         ja(i) = ia(i)
         ia(i) = j
      end do
      END
      subroutine xerbla ( srname, info )

c*********************************************************************72
c
cc XERBLA is an error handler for the level 2 BLAS routines.
c
      implicit none

      INTEGER            INFO
      CHARACTER*6        SRNAME
*     ..
*
*  Purpose
*  =======
*
*  XERBLA  is an error handler for the Level 2 BLAS routines.
*
*  It is called by the Level 2 BLAS routines if an input parameter is
*  invalid.
*
*  Installers should consider modifying the STOP statement in order to
*  call system-specific exception-handling facilities.
*
*  Parameters
*  ==========
*
*  SRNAME - CHARACTER*6.
*           On entry, SRNAME specifies the name of the routine which
*           called XERBLA.
*
*  INFO   - INTEGER.
*           On entry, INFO specifies the position of the invalid
*           parameter in the parameter-list of the calling routine.
*
*
*  Auxiliary routine for Level 2 Blas.
*
*  Written on 20-July-1986.
*
*     .. Executable Statements ..
*
      WRITE (*,99999) SRNAME, INFO
*
      STOP
*
99999 FORMAT ( ' ** On entry to ', A6, ' parameter number ', I2,
     $         ' had an illegal value' )
*
*     End of XERBLA.
*
      END
      subroutine zaxpy ( n, za, zx, incx, zy, incy )

c*********************************************************************72
c
cc ZAXPY adds a multiple of one vector to another.
c
c     jack dongarra, 3/11/78.
c
      implicit none

      double precision dcabs1
      integer i
      integer incx
      integer incy
      integer ix
      integer iy
      integer n
      complex*16 za
      complex*16 zx(1)
      complex*16 zy(1)

      if(n.le.0)return
      if (dcabs1(za) .eq. 0.0d0) return
      if (incx.eq.1.and.incy.eq.1)go to 20
c
c        code for unequal increments or equal increments
c          not equal to 1
c
      ix = 1
      iy = 1
      if(incx.lt.0)ix = (-n+1)*incx + 1
      if(incy.lt.0)iy = (-n+1)*incy + 1
      do 10 i = 1,n
        zy(iy) = zy(iy) + za*zx(ix)
        ix = ix + incx
        iy = iy + incy
   10 continue
      return
c
c        code for both increments equal to 1
c
   20 do 30 i = 1,n
        zy(i) = zy(i) + za*zx(i)
   30 continue
      return
      end
      subroutine zcmpac ( n, nx, ix, ixx, xx, iwsp )

c*********************************************************************72
c
cc ZCMPAC compacts the array IX and sorts IXX and XX.
c
*--   (This is a gateway routine for ZGCNVR) ...
c
      implicit none

      integer          n, nx, ix(nx), ixx(nx), iwsp(n)
      complex*16       xx(nx)
      integer          k
*
*---  sort ix and carry ixx and xx along ...
*
      call izsrt2( nx, ix, ixx, xx )
*
*---  adjust pointers ...
*
      do k = 1,n
         iwsp(k) = 0
      end do
      do k = 1,nx
         iwsp(ix(k)) = iwsp(ix(k)) + 1
      end do
      ix(n+1) = nx + 1
      do k = n,1,-1
         ix(k) = ix(k+1)-iwsp(k)
      end do
*
*---  sort ixx in increasing order and carry xx along ...
*
      do k = 1,n
         call izsrt1( iwsp(k), ixx(ix(k)), xx(ix(k)) )
      end do
      END
      subroutine zcopy ( n, zx, incx, zy, incy )

c*********************************************************************72
c
cc ZCOPY copies a vector.
c
c     jack dongarra, linpack, 4/11/78.
c
      implicit none

      double complex zx(1),zy(1)
      integer i,incx,incy,ix,iy,n
c
      if(n.le.0)return
      if(incx.eq.1.and.incy.eq.1)go to 20
c
c        code for unequal increments or equal increments
c          not equal to 1
c
      ix = 1
      iy = 1
      if(incx.lt.0)ix = (-n+1)*incx + 1
      if(incy.lt.0)iy = (-n+1)*incy + 1
      do 10 i = 1,n
        zy(iy) = zx(ix)
        ix = ix + incx
        iy = iy + incy
   10 continue
      return
c
c        code for both increments equal to 1
c
   20 do 30 i = 1,n
        zy(i) = zx(i)
   30 continue
      return
      end
      function zdotc(n,zx,incx,zy,incy)

c*********************************************************************72
c
cc ZDOTC computes the conjugated dot product of complex vectors.
c
c     jack dongarra, 3/11/78.
c
      implicit none

      integer i
      integer incx
      integer incy
      integer ix
      integer iy
      integer n
      double complex zdotc
      double complex zx(1),zy(1),ztemp

      ztemp = (0.0d0,0.0d0)
      zdotc = (0.0d0,0.0d0)
      if(n.le.0)return
      if(incx.eq.1.and.incy.eq.1)go to 20
c
c        code for unequal increments or equal increments
c          not equal to 1
c
      ix = 1
      iy = 1
      if(incx.lt.0)ix = (-n+1)*incx + 1
      if(incy.lt.0)iy = (-n+1)*incy + 1
      do 10 i = 1,n
        ztemp = ztemp + dconjg(zx(ix))*zy(iy)
        ix = ix + incx
        iy = iy + incy
   10 continue
      zdotc = ztemp
      return
c
c        code for both increments equal to 1
c
   20 do 30 i = 1,n
        ztemp = ztemp + dconjg(zx(i))*zy(i)
   30 continue
      zdotc = ztemp
      return
      end
      function zdotu(n,zx,incx,zy,incy)

c*********************************************************************72
c
cc ZDOTU computes the unconjugated dot product of complex vectors.
c
c     jack dongarra, 3/11/78.
c
      implicit none

      integer i
      integer incx
      integer incy
      integer ix
      integer iy
      integer n
      double complex zdotu
      double complex zx(1),zy(1),ztemp
      ztemp = (0.0d0,0.0d0)
      zdotu = (0.0d0,0.0d0)
      if(n.le.0)return
      if(incx.eq.1.and.incy.eq.1)go to 20
c
c        code for unequal increments or equal increments
c          not equal to 1
c
      ix = 1
      iy = 1
      if(incx.lt.0)ix = (-n+1)*incx + 1
      if(incy.lt.0)iy = (-n+1)*incy + 1
      do 10 i = 1,n
        ztemp = ztemp + zx(ix)*zy(iy)
        ix = ix + incx
        iy = iy + incy
   10 continue
      zdotu = ztemp
      return
c
c        code for both increments equal to 1
c
   20 do 30 i = 1,n
        ztemp = ztemp + zx(i)*zy(i)
   30 continue
      zdotu = ztemp
      return
      end
      subroutine zdscal ( n, da, zx, incx )

c*********************************************************************72
c
cc ZDSCAL scales a complex vector by a real constant.
c
c     jack dongarra, 3/11/78.
c     modified 3/93 to return if incx .le. 0.
c
      implicit none

      double complex zx(1)
      double precision da
      integer i,incx,ix,n
c
      if( n.le.0 .or. incx.le.0 )return
      if(incx.eq.1)go to 20
c
c        code for increment not equal to 1
c
      ix = 1
      do 10 i = 1,n
        zx(ix) = dcmplx(da,0.0d0)*zx(ix)
        ix = ix + incx
   10 continue
      return
c
c        code for increment equal to 1
c
   20 do 30 i = 1,n
        zx(i) = dcmplx(da,0.0d0)*zx(i)
   30 continue
      return
      end
      subroutine zgccsv ( x, y )

c*********************************************************************72
c
cc ZGCCSV computes y = A * x, for A in compressed column storage (CCS).
c
      implicit none

      complex*16 x(*), y(*)
*
*---  Computes y = A*x. A is passed via a fortran `common statement'.
*---  A is assumed to be under the Compress Column Storage (CCS) format.
*
      integer n, nz, nzmax
      parameter( nzmax = 50000 )
      integer ia(nzmax), ja(nzmax)
      complex*16 a(nzmax)
      common /CMAT/ a, ia, ja, nz, n

      integer i, j
      complex*16 ZERO
      parameter( ZERO=(0.0d0,0.0d0) )

      do i = 1,n
         y(i) = ZERO
      end do
      do j = 1,n
         do i = ja(j),ja(j+1)-1
            y(ia(i)) = y(ia(i)) + a(i)*x(j)
         end do
      end do
      end
      subroutine zgchbv ( m, t, h,ldh, y, wsp, iwsp, iflag )

c*********************************************************************72
c
cc ZGCHBV computes y = exp(t*H)*y using a partial fraction expansion.
c
      implicit none

      integer          m, ldh, iflag, iwsp(m)
      double precision t
      complex*16       H(ldh,m), y(m), wsp(m*(m+2))

*-----Purpose----------------------------------------------------------|
*
*---  ZGCHBV computes y = exp(t*H)*y using the partial fraction
*     expansion of the uniform rational Chebyshev approximation
*     to exp(-x) of type (14,14). H is a General matrix.
*     About 14-digit accuracy is expected if the matrix H is negative
*     definite. The algorithm may behave poorly otherwise.
*
*-----Arguments--------------------------------------------------------|
*
*     m       : (input) order of the matrix H.
*
*     t       : (input) time-scaling factor (can be < 0).
*
*     H(ldh,m): (input) argument matrix.
*
*     y(m)    : (input/output) on input the operand vector,
*               on output the resulting vector exp(t*H)*y.
*
*     iwsp(m) : (workspace)
*
*     wsp     : (workspace). Observe that a double precision vector of
*               length 2*m*(m+2) can be used as well when calling this
*               routine.
*
*----------------------------------------------------------------------|
*     Roger B. Sidje (rbs@maths.uq.edu.au)
*     EXPOKIT: Software Package for Computing Matrix Exponentials.
*     ACM - Transactions On Mathematical Software, 24(1):130-156, 1998
*----------------------------------------------------------------------|
*
      integer     ndeg, i, j, ip, ih, iy, iz
      parameter ( ndeg=7 )
      double      precision alpha0
      complex*16  alpha(2*ndeg), theta(2*ndeg)

*---  Pointers ...

      ih = 1
      iy = ih + m*m
      iz = iy + m

*---  Coefficients and poles of the partial fraction expansion ...

      alpha0  =  0.183216998528140087D-11
      alpha(1)=( 0.557503973136501826D+02,-0.204295038779771857D+03)
      alpha(2)=(-0.938666838877006739D+02, 0.912874896775456363D+02)
      alpha(3)=( 0.469965415550370835D+02,-0.116167609985818103D+02)
      alpha(4)=(-0.961424200626061065D+01,-0.264195613880262669D+01)
      alpha(5)=( 0.752722063978321642D+00, 0.670367365566377770D+00)
      alpha(6)=(-0.188781253158648576D-01,-0.343696176445802414D-01)
      alpha(7)=( 0.143086431411801849D-03, 0.287221133228814096D-03)

      theta(1)=(-0.562314417475317895D+01, 0.119406921611247440D+01)
      theta(2)=(-0.508934679728216110D+01, 0.358882439228376881D+01)
      theta(3)=(-0.399337136365302569D+01, 0.600483209099604664D+01)
      theta(4)=(-0.226978543095856366D+01, 0.846173881758693369D+01)
      theta(5)=( 0.208756929753827868D+00, 0.109912615662209418D+02)
      theta(6)=( 0.370327340957595652D+01, 0.136563731924991884D+02)
      theta(7)=( 0.889777151877331107D+01, 0.166309842834712071D+02)
*
      do ip = 1,ndeg
         theta(ndeg+ip) = CONJG( theta(ip) )
         alpha(ndeg+ip) = CONJG( alpha(ip) )
      end do
*
*---  Accumulation of the contribution of each pole ...
*
      do j = 1,m
         wsp(iz+j-1) = y(j)
         y(j) = y(j)*alpha0
      end do
      do ip = 1,2*ndeg
         alpha(ip) = 0.5d0*alpha(ip)
*---     Solve each fraction using Gaussian elimination with pivoting...
         do j = 1,m
            do i = 1,m
               wsp(ih+(j-1)*m+i-1) = -t*H(i,j)
            end do
            wsp(ih+(j-1)*m+j-1) = wsp(ih+(j-1)*m+j-1)-theta(ip)
            wsp(iy+j-1) = wsp(iz+j-1)
         end do
         call ZGESV( M, 1, WSP(iH),M, IWSP, WSP(iY),M, IFLAG )
         if ( IFLAG.ne.0 ) stop 'Error in ZGCHBV'
*---     Accumulate the partial result in y ...
         do i = 1,m
            y(i) = y(i) + alpha(ip)*wsp(iy+i-1)
         end do
      end do
      END
      subroutine zgcnvr ( from,to,diag, nrow,ncol, nz, ia, ja, a, iwsp )

c*********************************************************************72
c
cc ZGCNVR transforms a sparse storage format into another sparse storage format.
c
      implicit none

      character        from*3, to*3, diag*1
      integer          nrow, ncol, nz, ia(*), ja(*), iwsp(*)
      complex*16       a(*)

*-----Purpose----------------------------------------------------------|
*
*---  ZGCNVR transforms a sparse storage format into another sparse
*     storage format. The matrix can be rectangular.
*     This is the complex counterpart of DGCNVR.
*
*-----Arguments--------------------------------------------------------|
*
*     from   : (input, character*3) the storage format holding the
*              matrix on input. Accepted values are:
*              'COO' : COOrdinates
*              'CRS' : Compressed Row Storage
*              'CCS' : Compressed Column Storage.
*
*     to     : (input, character*3) the storage format holding the
*              matrix on output. Same accepted values as above.
*
*     diag   : (input, character*1) specifies whether or not the
*              entire diagonal should be lodged, i.e., including
*              null diagonal entries. (This may be required by
*              some applications that make use of the locations
*              of diagonal elements.)
*              diag = 'N', no specific attention is given to diagonal
*                          elements.
*              diag = 'D', the entire diagonal is lodged, including
*                          null elements on the diagonal.
*              if from=to & diag='D', null diagonal entries are
*              explicitly inserted in the actual matrix.
*
*     nrow   : (input) number of rows in the matrix.
*
*     ncol   : (input) number of columns in the matrix.
*
*     nz     : (input/output) number of non-zeros elements.
*              If diag='D' and null diagonal entries are inserted, then
*              nz is updated on exit and contains the effective number
*              of entries stored. In what follows, nz' (read nz prime)
*              denotes the updated value of nz, nz <= nz' <= nz+n.
*              If diag='N' then nz'=nz.
*
*     ia(*) : (input/output) of declared length .ge. nz'.
*             On input,
*                if from='CRS', ia(1:nrow+1) contains pointers for the
*                beginning of each row.
*                if from='COO', or 'CCS', ia(1:nz) contains row indices.
*             On output,
*                if to='CRS', ia(1:nrow+1) contains pointers for the
*                beginning of each row.
*                if to='COO' or 'CCS', ia(1:nz') contains row indices
*                in increasing order in each column.
*
*     ja(*) : (input/output) of declared length .ge. nz'.
*             On input,
*                if from='CRS', ja(1:ncol+1) contains pointers for the
*                beginning of each column.
*                if from='COO' or 'CCS', ja(1:nz) contains col. indices.
*             On output,
*                if to='CRS', ja(1:ncol+1) contains pointers for the
*                beginning of each column.
*                if to='COO' or 'CCS', ja(1:nz') contains column indices
*                in increasing order in each row.
*
*     a(*)  : (input/output) On input, a(1:nz) fits the input format
*             and on output, a(1:nz') fits the output format.
*
*   iwsp(*) : (workspace) of declared length .ge. max(nrow,ncol).
*
*----------------------------------------------------------------------|
*     Roger B. Sidje (rbs@maths.uq.edu.au)
*     Department of Mathematics, University of Queensland.
*     Brisbane QLD 4072, Australia. 1996.
*----------------------------------------------------------------------|

      integer i, j, k, nn, iflag
      character cpy_from*3, cpy_to*3, cpy_diag*1, c*1, S1*3, S2*3
      logical LSAME

      intrinsic CHAR,ICHAR,LLE,LGE,MIN
      LSAME(S1,S2) = LLE(S1,S2).and.LGE(S1,S2)
*
*---  upper case strings ...
*
      do k = 1,3
         c = from(k:k)
         if ( c.ge.'a' .and. c.le.'z' ) c = CHAR( ICHAR(c)-32 )
         cpy_from(k:k) = c
         c = to(k:k)
         if ( c.ge.'a' .and. c.le.'z' ) c = CHAR( ICHAR(c)-32 )
         cpy_to(k:k) = c
      end do
      c = diag(1:1)
      if ( c.ge.'a' .and. c.le.'z' ) c = CHAR( ICHAR(c)-32 )
      cpy_diag = c

*---  quick return if possible ...
      if ( LSAME(cpy_from,cpy_to) .and. cpy_diag.eq.'N' ) return
      iflag = 1
      if ( LSAME(cpy_from,'COO') ) iflag = 0
      if ( LSAME(cpy_from,'CRS') ) iflag = 0
      if ( LSAME(cpy_from,'CCS') ) iflag = 0
      if ( iflag.ne.0 ) stop 'unexpected i/o formats (in ZGCNVR)'
      iflag = 1
      if ( LSAME(cpy_to,'COO') )   iflag = 0
      if ( LSAME(cpy_to,'CRS') )   iflag = 0
      if ( LSAME(cpy_to,'CCS') )   iflag = 0
      if ( iflag.ne.0 ) stop 'unexpected i/o formats (in ZGCNVR)'

*
*---  transit via COOrdinates format if input is not in COO ...
*
      if ( LSAME(cpy_from,'CRS') ) then
*---     expand ia indices ...
         nn = nz
         do i = nrow,1,-1
            do k = 1,ia(i+1)-ia(i)
               ia(nn) = i
               nn = nn-1
            end do
         end do
      endif
      if ( LSAME(cpy_from,'CCS') ) then
*---     expand ja indices ...
         nn = nz
         do j = ncol,1,-1
            do k = 1,ja(j+1)-ja(j)
               ja(nn) = j
               nn = nn-1
            end do
         end do
      endif
*
*--   if requested, insert diagonal elements even if they are zero...
*
      if ( cpy_diag.eq.'D' ) then
         nn = MIN( nrow,ncol )
         do i = 1,nn
            iwsp(i) = 0
         end do
         do k = 1,nz
            if ( ia(k).eq.ja(k) ) iwsp(ia(k)) = 1
         end do
         do i = 1,nn
            if ( iwsp(i).eq.0 ) then
               nz = nz + 1
               ia(nz) = i
               ja(nz) = i
               a(nz) = 0.0d0
            endif
         end do
      endif

*---  COO convertion ...
      if ( LSAME(cpy_to,'COO') ) return

*---  CRS convertion ...
      if ( LSAME(cpy_to,'CRS') ) call zcmpac( nrow,nz,ia,ja,a, iwsp )

*---  CCS convertion ...
      if ( LSAME(cpy_to,'CCS') ) call zcmpac( ncol,nz,ja,ia,a, iwsp )

      END
      subroutine zgcoov ( x, y )

c*********************************************************************72
c
cc ZGCOOV computes y = A * x, with A stored in COO storage format.
c
      implicit none

      complex*16 x(*), y(*)
*
*---  Computes y = A*x. A is passed via a fortran `common statement'.
*---  A is assumed here to be under the COOrdinates storage format.
*
*----------------------------------------------------------------------|
*-------------------------------NOTE-----------------------------------*
*     This is an accessory to Expokit and it is not intended to be     *
*     complete. It is supplied primarily to ensure an unconstrained    *
*     distribution and portability of the package. The matrix-vector   *
*     multiplication routines supplied here fit the non symmetric      *
*     storage and for a symmetric matrix, the entire (not half) matrix *
*     is required.  If the sparsity pattern is known a priori, it is   *
*     recommended to use the most advantageous format and to devise    *
*     the most advantageous matrix-vector multiplication routine.      *
*----------------------------------------------------------------------*
*----------------------------------------------------------------------*
      integer n, nz, nzmax
      parameter( nzmax = 50000 )
      integer ia(nzmax), ja(nzmax)
      complex*16 a(nzmax)
      common /CMAT/ a, ia, ja, nz, n

      integer i, j
      complex*16 ZERO
      parameter( ZERO=(0.0d0,0.0d0) )

      do j = 1,n
         y(j) = ZERO
      end do
      do i = 1,nz
         y(ia(i)) = y(ia(i)) + a(i)*x(ja(i))
      end do
      END
      subroutine zgcrsv ( x, y )

c*********************************************************************72
c
cc ZGCRSV computes y = A * x, with A stored in compressed row storage (CRS).
c
      implicit none

      complex*16 x(*), y(*)
*
*---  Computes y = A*x. A is passed via a fortran `common statement'.
*---  A is assumed to be under the Compress Row Storage (CRS) format.
*
      integer n, nz, nzmax
      parameter( nzmax = 50000 )
      integer ia(nzmax), ja(nzmax)
      complex*16 a(nzmax)
      common /CMAT/ a, ia, ja, nz, n

      integer i, j
      complex*16 ZERO
      parameter( ZERO=(0.0d0,0.0d0) )

      do i = 1,n
         y(i) = ZERO
         do j = ia(i),ia(i+1)-1
            y(i) = y(i) + a(j)*x(ja(j))
         end do
      end do
      END
      subroutine zgefa ( a, lda, n, ipvt, info )

c*********************************************************************72
c
cc ZGEFA factors a complex general linear system.c
c
      implicit none

      integer lda,n,ipvt(1),info
      complex*16 a(lda,1)
c
c     zgefa factors a complex*16 matrix by gaussian elimination.
c
c     zgefa is usually called by zgeco, but it can be called
c     directly with a saving in time if  rcond  is not needed.
c     (time for zgeco) = (1 + 9/n)*(time for zgefa) .
c
c     on entry
c
c        a       complex*16(lda, n)
c                the matrix to be factored.
c
c        lda     integer
c                the leading dimension of the array  a .
c
c        n       integer
c                the order of the matrix  a .
c
c     on return
c
c        a       an upper triangular matrix and the multipliers
c                which were used to obtain it.
c                the factorization can be written  a = l*u  where
c                l  is a product of permutation and unit lower
c                triangular matrices and  u  is upper triangular.
c
c        ipvt    integer(n)
c                an integer vector of pivot indices.
c
c        info    integer
c                = 0  normal value.
c                = k  if  u(k,k) .eq. 0.0 .  this is not an error
c                     condition for this subroutine, but it does
c                     indicate that zgesl or zgedi will divide by zero
c                     if called.  use  rcond  in zgeco for a reliable
c                     indication of singularity.
c
c     linpack. this version dated 08/14/78 .
c     cleve moler, university of new mexico, argonne national lab.
c
c     subroutines and functions
c
c     blas zaxpy,zscal,izamax
c     fortran dabs
c
c     internal variables
c
      complex*16 t
      integer izamax,j,k,kp1,l,nm1
c
      complex*16 zdum
      double precision cabs1
      double precision dreal,dimag
      complex*16 zdumr,zdumi
      dreal(zdumr) = zdumr
      dimag(zdumi) = (0.0d0,-1.0d0)*zdumi
      cabs1(zdum) = dabs(dreal(zdum)) + dabs(dimag(zdum))
c
c     gaussian elimination with partial pivoting
c
      info = 0
      nm1 = n - 1
      if (nm1 .lt. 1) go to 70
      do 60 k = 1, nm1
         kp1 = k + 1
c
c        find l = pivot index
c
         l = izamax(n-k+1,a(k,k),1) + k - 1
         ipvt(k) = l
c
c        zero pivot implies this column already triangularized
c
         if (cabs1(a(l,k)) .eq. 0.0d0) go to 40
c
c           interchange if necessary
c
            if (l .eq. k) go to 10
               t = a(l,k)
               a(l,k) = a(k,k)
               a(k,k) = t
   10       continue
c
c           compute multipliers
c
            t = -(1.0d0,0.0d0)/a(k,k)
            call zscal(n-k,t,a(k+1,k),1)
c
c           row elimination with column indexing
c
            do 30 j = kp1, n
               t = a(l,j)
               if (l .eq. k) go to 20
                  a(l,j) = a(k,j)
                  a(k,j) = t
   20          continue
               call zaxpy(n-k,t,a(k+1,k),1,a(k+1,j),1)
   30       continue
         go to 50
   40    continue
            info = k
   50    continue
   60 continue
   70 continue
      ipvt(n) = n
      if (cabs1(a(n,n)) .eq. 0.0d0) info = n
      return
      end
      subroutine zgemm ( transa, transb, m, n, k, alpha, a, lda, b, ldb,
     $                   beta, c, ldc )

c*********************************************************************72
c
cc ZGEMM performs a matrix-matrix operation of the form C=alpha*A*B+beta*C.
c
      implicit none

      CHARACTER*1        TRANSA, TRANSB
      INTEGER            M, N, K, LDA, LDB, LDC
      COMPLEX*16         ALPHA, BETA
*     .. Array Arguments ..
      COMPLEX*16         A( LDA, * ), B( LDB, * ), C( LDC, * )
*     ..
*
*  Purpose
*  =======
*
*  ZGEMM  performs one of the matrix-matrix operations
*
*     C := alpha*op( A )*op( B ) + beta*C,
*
*  where  op( X ) is one of
*
*     op( X ) = X   or   op( X ) = X'   or   op( X ) = conjg( X' ),
*
*  alpha and beta are scalars, and A, B and C are matrices, with op( A )
*  an m by k matrix,  op( B )  a  k by n matrix and  C an m by n matrix.
*
*  Parameters
*  ==========
*
*  TRANSA - CHARACTER*1.
*           On entry, TRANSA specifies the form of op( A ) to be used in
*           the matrix multiplication as follows:
*
*              TRANSA = 'N' or 'n',  op( A ) = A.
*
*              TRANSA = 'T' or 't',  op( A ) = A'.
*
*              TRANSA = 'C' or 'c',  op( A ) = conjg( A' ).
*
*           Unchanged on exit.
*
*  TRANSB - CHARACTER*1.
*           On entry, TRANSB specifies the form of op( B ) to be used in
*           the matrix multiplication as follows:
*
*              TRANSB = 'N' or 'n',  op( B ) = B.
*
*              TRANSB = 'T' or 't',  op( B ) = B'.
*
*              TRANSB = 'C' or 'c',  op( B ) = conjg( B' ).
*
*           Unchanged on exit.
*
*  M      - INTEGER.
*           On entry,  M  specifies  the number  of rows  of the  matrix
*           op( A )  and of the  matrix  C.  M  must  be at least  zero.
*           Unchanged on exit.
*
*  N      - INTEGER.
*           On entry,  N  specifies the number  of columns of the matrix
*           op( B ) and the number of columns of the matrix C. N must be
*           at least zero.
*           Unchanged on exit.
*
*  K      - INTEGER.
*           On entry,  K  specifies  the number of columns of the matrix
*           op( A ) and the number of rows of the matrix op( B ). K must
*           be at least  zero.
*           Unchanged on exit.
*
*  ALPHA  - COMPLEX*16      .
*           On entry, ALPHA specifies the scalar alpha.
*           Unchanged on exit.
*
*  A      - COMPLEX*16       array of DIMENSION ( LDA, ka ), where ka is
*           k  when  TRANSA = 'N' or 'n',  and is  m  otherwise.
*           Before entry with  TRANSA = 'N' or 'n',  the leading  m by k
*           part of the array  A  must contain the matrix  A,  otherwise
*           the leading  k by m  part of the array  A  must contain  the
*           matrix A.
*           Unchanged on exit.
*
*  LDA    - INTEGER.
*           On entry, LDA specifies the first dimension of A as declared
*           in the calling (sub) program. When  TRANSA = 'N' or 'n' then
*           LDA must be at least  max( 1, m ), otherwise  LDA must be at
*           least  max( 1, k ).
*           Unchanged on exit.
*
*  B      - COMPLEX*16       array of DIMENSION ( LDB, kb ), where kb is
*           n  when  TRANSB = 'N' or 'n',  and is  k  otherwise.
*           Before entry with  TRANSB = 'N' or 'n',  the leading  k by n
*           part of the array  B  must contain the matrix  B,  otherwise
*           the leading  n by k  part of the array  B  must contain  the
*           matrix B.
*           Unchanged on exit.
*
*  LDB    - INTEGER.
*           On entry, LDB specifies the first dimension of B as declared
*           in the calling (sub) program. When  TRANSB = 'N' or 'n' then
*           LDB must be at least  max( 1, k ), otherwise  LDB must be at
*           least  max( 1, n ).
*           Unchanged on exit.
*
*  BETA   - COMPLEX*16      .
*           On entry,  BETA  specifies the scalar  beta.  When  BETA  is
*           supplied as zero then C need not be set on input.
*           Unchanged on exit.
*
*  C      - COMPLEX*16       array of DIMENSION ( LDC, n ).
*           Before entry, the leading  m by n  part of the array  C must
*           contain the matrix  C,  except when  beta  is zero, in which
*           case C need not be set on entry.
*           On exit, the array  C  is overwritten by the  m by n  matrix
*           ( alpha*op( A )*op( B ) + beta*C ).
*
*  LDC    - INTEGER.
*           On entry, LDC specifies the first dimension of C as declared
*           in  the  calling  (sub)  program.   LDC  must  be  at  least
*           max( 1, m ).
*           Unchanged on exit.
*
*
*  Level 3 Blas routine.
*
*  -- Written on 8-February-1989.
*     Jack Dongarra, Argonne National Laboratory.
*     Iain Duff, AERE Harwell.
*     Jeremy Du Croz, Numerical Algorithms Group Ltd.
*     Sven Hammarling, Numerical Algorithms Group Ltd.
*
*
*     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
*     .. External Subroutines ..
      EXTERNAL           XERBLA
*     .. Intrinsic Functions ..
      INTRINSIC          DCONJG, MAX
*     .. Local Scalars ..
      LOGICAL            CONJA, CONJB, NOTA, NOTB
      INTEGER            I, INFO, J, L, NCOLA, NROWA, NROWB
      COMPLEX*16         TEMP
*     .. Parameters ..
      COMPLEX*16         ONE
      PARAMETER        ( ONE  = ( 1.0D+0, 0.0D+0 ) )
      COMPLEX*16         ZERO
      PARAMETER        ( ZERO = ( 0.0D+0, 0.0D+0 ) )
*     ..
*     .. Executable Statements ..
*
*     Set  NOTA  and  NOTB  as  true if  A  and  B  respectively are not
*     conjugated or transposed, set  CONJA and CONJB  as true if  A  and
*     B  respectively are to be  transposed but  not conjugated  and set
*     NROWA, NCOLA and  NROWB  as the number of rows and  columns  of  A
*     and the number of rows of  B  respectively.
*
      NOTA  = LSAME( TRANSA, 'N' )
      NOTB  = LSAME( TRANSB, 'N' )
      CONJA = LSAME( TRANSA, 'C' )
      CONJB = LSAME( TRANSB, 'C' )
      IF( NOTA )THEN
         NROWA = M
         NCOLA = K
      ELSE
         NROWA = K
         NCOLA = M
      END IF
      IF( NOTB )THEN
         NROWB = K
      ELSE
         NROWB = N
      END IF
*
*     Test the input parameters.
*
      INFO = 0
      IF(      ( .NOT.NOTA                 ).AND.
     $         ( .NOT.CONJA                ).AND.
     $         ( .NOT.LSAME( TRANSA, 'T' ) )      )THEN
         INFO = 1
      ELSE IF( ( .NOT.NOTB                 ).AND.
     $         ( .NOT.CONJB                ).AND.
     $         ( .NOT.LSAME( TRANSB, 'T' ) )      )THEN
         INFO = 2
      ELSE IF( M  .LT.0               )THEN
         INFO = 3
      ELSE IF( N  .LT.0               )THEN
         INFO = 4
      ELSE IF( K  .LT.0               )THEN
         INFO = 5
      ELSE IF( LDA.LT.MAX( 1, NROWA ) )THEN
         INFO = 8
      ELSE IF( LDB.LT.MAX( 1, NROWB ) )THEN
         INFO = 10
      ELSE IF( LDC.LT.MAX( 1, M     ) )THEN
         INFO = 13
      END IF
      IF( INFO.NE.0 )THEN
         CALL XERBLA( 'ZGEMM ', INFO )
         RETURN
      END IF
*
*     Quick return if possible.
*
      IF( ( M.EQ.0 ).OR.( N.EQ.0 ).OR.
     $    ( ( ( ALPHA.EQ.ZERO ).OR.( K.EQ.0 ) ).AND.( BETA.EQ.ONE ) ) )
     $   RETURN
*
*     And when  alpha.eq.zero.
*
      IF( ALPHA.EQ.ZERO )THEN
         IF( BETA.EQ.ZERO )THEN
            DO 20, J = 1, N
               DO 10, I = 1, M
                  C( I, J ) = ZERO
   10          CONTINUE
   20       CONTINUE
         ELSE
            DO 40, J = 1, N
               DO 30, I = 1, M
                  C( I, J ) = BETA*C( I, J )
   30          CONTINUE
   40       CONTINUE
         END IF
         RETURN
      END IF
*
*     Start the operations.
*
      IF( NOTB )THEN
         IF( NOTA )THEN
*
*           Form  C := alpha*A*B + beta*C.
*
            DO 90, J = 1, N
               IF( BETA.EQ.ZERO )THEN
                  DO 50, I = 1, M
                     C( I, J ) = ZERO
   50             CONTINUE
               ELSE IF( BETA.NE.ONE )THEN
                  DO 60, I = 1, M
                     C( I, J ) = BETA*C( I, J )
   60             CONTINUE
               END IF
               DO 80, L = 1, K
                  IF( B( L, J ).NE.ZERO )THEN
                     TEMP = ALPHA*B( L, J )
                     DO 70, I = 1, M
                        C( I, J ) = C( I, J ) + TEMP*A( I, L )
   70                CONTINUE
                  END IF
   80          CONTINUE
   90       CONTINUE
         ELSE IF( CONJA )THEN
*
*           Form  C := alpha*conjg( A' )*B + beta*C.
*
            DO 120, J = 1, N
               DO 110, I = 1, M
                  TEMP = ZERO
                  DO 100, L = 1, K
                     TEMP = TEMP + DCONJG( A( L, I ) )*B( L, J )
  100             CONTINUE
                  IF( BETA.EQ.ZERO )THEN
                     C( I, J ) = ALPHA*TEMP
                  ELSE
                     C( I, J ) = ALPHA*TEMP + BETA*C( I, J )
                  END IF
  110          CONTINUE
  120       CONTINUE
         ELSE
*
*           Form  C := alpha*A'*B + beta*C
*
            DO 150, J = 1, N
               DO 140, I = 1, M
                  TEMP = ZERO
                  DO 130, L = 1, K
                     TEMP = TEMP + A( L, I )*B( L, J )
  130             CONTINUE
                  IF( BETA.EQ.ZERO )THEN
                     C( I, J ) = ALPHA*TEMP
                  ELSE
                     C( I, J ) = ALPHA*TEMP + BETA*C( I, J )
                  END IF
  140          CONTINUE
  150       CONTINUE
         END IF
      ELSE IF( NOTA )THEN
         IF( CONJB )THEN
*
*           Form  C := alpha*A*conjg( B' ) + beta*C.
*
            DO 200, J = 1, N
               IF( BETA.EQ.ZERO )THEN
                  DO 160, I = 1, M
                     C( I, J ) = ZERO
  160             CONTINUE
               ELSE IF( BETA.NE.ONE )THEN
                  DO 170, I = 1, M
                     C( I, J ) = BETA*C( I, J )
  170             CONTINUE
               END IF
               DO 190, L = 1, K
                  IF( B( J, L ).NE.ZERO )THEN
                     TEMP = ALPHA*DCONJG( B( J, L ) )
                     DO 180, I = 1, M
                        C( I, J ) = C( I, J ) + TEMP*A( I, L )
  180                CONTINUE
                  END IF
  190          CONTINUE
  200       CONTINUE
         ELSE
*
*           Form  C := alpha*A*B'          + beta*C
*
            DO 250, J = 1, N
               IF( BETA.EQ.ZERO )THEN
                  DO 210, I = 1, M
                     C( I, J ) = ZERO
  210             CONTINUE
               ELSE IF( BETA.NE.ONE )THEN
                  DO 220, I = 1, M
                     C( I, J ) = BETA*C( I, J )
  220             CONTINUE
               END IF
               DO 240, L = 1, K
                  IF( B( J, L ).NE.ZERO )THEN
                     TEMP = ALPHA*B( J, L )
                     DO 230, I = 1, M
                        C( I, J ) = C( I, J ) + TEMP*A( I, L )
  230                CONTINUE
                  END IF
  240          CONTINUE
  250       CONTINUE
         END IF
      ELSE IF( CONJA )THEN
         IF( CONJB )THEN
*
*           Form  C := alpha*conjg( A' )*conjg( B' ) + beta*C.
*
            DO 280, J = 1, N
               DO 270, I = 1, M
                  TEMP = ZERO
                  DO 260, L = 1, K
                     TEMP = TEMP +
     $                      DCONJG( A( L, I ) )*DCONJG( B( J, L ) )
  260             CONTINUE
                  IF( BETA.EQ.ZERO )THEN
                     C( I, J ) = ALPHA*TEMP
                  ELSE
                     C( I, J ) = ALPHA*TEMP + BETA*C( I, J )
                  END IF
  270          CONTINUE
  280       CONTINUE
         ELSE
*
*           Form  C := alpha*conjg( A' )*B' + beta*C
*
            DO 310, J = 1, N
               DO 300, I = 1, M
                  TEMP = ZERO
                  DO 290, L = 1, K
                     TEMP = TEMP + DCONJG( A( L, I ) )*B( J, L )
  290             CONTINUE
                  IF( BETA.EQ.ZERO )THEN
                     C( I, J ) = ALPHA*TEMP
                  ELSE
                     C( I, J ) = ALPHA*TEMP + BETA*C( I, J )
                  END IF
  300          CONTINUE
  310       CONTINUE
         END IF
      ELSE
         IF( CONJB )THEN
*
*           Form  C := alpha*A'*conjg( B' ) + beta*C
*
            DO 340, J = 1, N
               DO 330, I = 1, M
                  TEMP = ZERO
                  DO 320, L = 1, K
                     TEMP = TEMP + A( L, I )*DCONJG( B( J, L ) )
  320             CONTINUE
                  IF( BETA.EQ.ZERO )THEN
                     C( I, J ) = ALPHA*TEMP
                  ELSE
                     C( I, J ) = ALPHA*TEMP + BETA*C( I, J )
                  END IF
  330          CONTINUE
  340       CONTINUE
         ELSE
*
*           Form  C := alpha*A'*B' + beta*C
*
            DO 370, J = 1, N
               DO 360, I = 1, M
                  TEMP = ZERO
                  DO 350, L = 1, K
                     TEMP = TEMP + A( L, I )*B( J, L )
  350             CONTINUE
                  IF( BETA.EQ.ZERO )THEN
                     C( I, J ) = ALPHA*TEMP
                  ELSE
                     C( I, J ) = ALPHA*TEMP + BETA*C( I, J )
                  END IF
  360          CONTINUE
  370       CONTINUE
         END IF
      END IF
*
      RETURN
*
*     End of ZGEMM .
*
      END
      subroutine zgemv ( trans, m, n, alpha, a, lda, x, incx,
     $                   beta, y, incy )

c*********************************************************************72
c
cc ZGEMV performs a matrix-vector operation of the form y=alpha*A*x+beta*y.
c
      implicit none

      COMPLEX*16         ALPHA, BETA
      INTEGER            INCX, INCY, LDA, M, N
      CHARACTER*1        TRANS
*     .. Array Arguments ..
      COMPLEX*16         A( LDA, * ), X( * ), Y( * )
*     ..
*
*  Purpose
*  =======
*
*  ZGEMV  performs one of the matrix-vector operations
*
*     y := alpha*A*x + beta*y,   or   y := alpha*A'*x + beta*y,   or
*
*     y := alpha*conjg( A' )*x + beta*y,
*
*  where alpha and beta are scalars, x and y are vectors and A is an
*  m by n matrix.
*
*  Parameters
*  ==========
*
*  TRANS  - CHARACTER*1.
*           On entry, TRANS specifies the operation to be performed as
*           follows:
*
*              TRANS = 'N' or 'n'   y := alpha*A*x + beta*y.
*
*              TRANS = 'T' or 't'   y := alpha*A'*x + beta*y.
*
*              TRANS = 'C' or 'c'   y := alpha*conjg( A' )*x + beta*y.
*
*           Unchanged on exit.
*
*  M      - INTEGER.
*           On entry, M specifies the number of rows of the matrix A.
*           M must be at least zero.
*           Unchanged on exit.
*
*  N      - INTEGER.
*           On entry, N specifies the number of columns of the matrix A.
*           N must be at least zero.
*           Unchanged on exit.
*
*  ALPHA  - COMPLEX*16      .
*           On entry, ALPHA specifies the scalar alpha.
*           Unchanged on exit.
*
*  A      - COMPLEX*16       array of DIMENSION ( LDA, n ).
*           Before entry, the leading m by n part of the array A must
*           contain the matrix of coefficients.
*           Unchanged on exit.
*
*  LDA    - INTEGER.
*           On entry, LDA specifies the first dimension of A as declared
*           in the calling (sub) program. LDA must be at least
*           max( 1, m ).
*           Unchanged on exit.
*
*  X      - COMPLEX*16       array of DIMENSION at least
*           ( 1 + ( n - 1 )*abs( INCX ) ) when TRANS = 'N' or 'n'
*           and at least
*           ( 1 + ( m - 1 )*abs( INCX ) ) otherwise.
*           Before entry, the incremented array X must contain the
*           vector x.
*           Unchanged on exit.
*
*  INCX   - INTEGER.
*           On entry, INCX specifies the increment for the elements of
*           X. INCX must not be zero.
*           Unchanged on exit.
*
*  BETA   - COMPLEX*16      .
*           On entry, BETA specifies the scalar beta. When BETA is
*           supplied as zero then Y need not be set on input.
*           Unchanged on exit.
*
*  Y      - COMPLEX*16       array of DIMENSION at least
*           ( 1 + ( m - 1 )*abs( INCY ) ) when TRANS = 'N' or 'n'
*           and at least
*           ( 1 + ( n - 1 )*abs( INCY ) ) otherwise.
*           Before entry with BETA non-zero, the incremented array Y
*           must contain the vector y. On exit, Y is overwritten by the
*           updated vector y.
*
*  INCY   - INTEGER.
*           On entry, INCY specifies the increment for the elements of
*           Y. INCY must not be zero.
*           Unchanged on exit.
*
*
*  Level 2 Blas routine.
*
*  -- Written on 22-October-1986.
*     Jack Dongarra, Argonne National Lab.
*     Jeremy Du Croz, Nag Central Office.
*     Sven Hammarling, Nag Central Office.
*     Richard Hanson, Sandia National Labs.
*
*
*     .. Parameters ..
      COMPLEX*16         ONE
      PARAMETER        ( ONE  = ( 1.0D+0, 0.0D+0 ) )
      COMPLEX*16         ZERO
      PARAMETER        ( ZERO = ( 0.0D+0, 0.0D+0 ) )
*     .. Local Scalars ..
      COMPLEX*16         TEMP
      INTEGER            I, INFO, IX, IY, J, JX, JY, KX, KY, LENX, LENY
      LOGICAL            NOCONJ
*     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
*     .. External Subroutines ..
      EXTERNAL           XERBLA
*     .. Intrinsic Functions ..
      INTRINSIC          DCONJG, MAX
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters.
*
      INFO = 0
      IF     ( .NOT.LSAME( TRANS, 'N' ).AND.
     $         .NOT.LSAME( TRANS, 'T' ).AND.
     $         .NOT.LSAME( TRANS, 'C' )      )THEN
         INFO = 1
      ELSE IF( M.LT.0 )THEN
         INFO = 2
      ELSE IF( N.LT.0 )THEN
         INFO = 3
      ELSE IF( LDA.LT.MAX( 1, M ) )THEN
         INFO = 6
      ELSE IF( INCX.EQ.0 )THEN
         INFO = 8
      ELSE IF( INCY.EQ.0 )THEN
         INFO = 11
      END IF
      IF( INFO.NE.0 )THEN
         CALL XERBLA( 'ZGEMV ', INFO )
         RETURN
      END IF
*
*     Quick return if possible.
*
      IF( ( M.EQ.0 ).OR.( N.EQ.0 ).OR.
     $    ( ( ALPHA.EQ.ZERO ).AND.( BETA.EQ.ONE ) ) )
     $   RETURN
*
      NOCONJ = LSAME( TRANS, 'T' )
*
*     Set  LENX  and  LENY, the lengths of the vectors x and y, and set
*     up the start points in  X  and  Y.
*
      IF( LSAME( TRANS, 'N' ) )THEN
         LENX = N
         LENY = M
      ELSE
         LENX = M
         LENY = N
      END IF
      IF( INCX.GT.0 )THEN
         KX = 1
      ELSE
         KX = 1 - ( LENX - 1 )*INCX
      END IF
      IF( INCY.GT.0 )THEN
         KY = 1
      ELSE
         KY = 1 - ( LENY - 1 )*INCY
      END IF
*
*     Start the operations. In this version the elements of A are
*     accessed sequentially with one pass through A.
*
*     First form  y := beta*y.
*
      IF( BETA.NE.ONE )THEN
         IF( INCY.EQ.1 )THEN
            IF( BETA.EQ.ZERO )THEN
               DO 10, I = 1, LENY
                  Y( I ) = ZERO
   10          CONTINUE
            ELSE
               DO 20, I = 1, LENY
                  Y( I ) = BETA*Y( I )
   20          CONTINUE
            END IF
         ELSE
            IY = KY
            IF( BETA.EQ.ZERO )THEN
               DO 30, I = 1, LENY
                  Y( IY ) = ZERO
                  IY      = IY   + INCY
   30          CONTINUE
            ELSE
               DO 40, I = 1, LENY
                  Y( IY ) = BETA*Y( IY )
                  IY      = IY           + INCY
   40          CONTINUE
            END IF
         END IF
      END IF
      IF( ALPHA.EQ.ZERO )
     $   RETURN
      IF( LSAME( TRANS, 'N' ) )THEN
*
*        Form  y := alpha*A*x + y.
*
         JX = KX
         IF( INCY.EQ.1 )THEN
            DO 60, J = 1, N
               IF( X( JX ).NE.ZERO )THEN
                  TEMP = ALPHA*X( JX )
                  DO 50, I = 1, M
                     Y( I ) = Y( I ) + TEMP*A( I, J )
   50             CONTINUE
               END IF
               JX = JX + INCX
   60       CONTINUE
         ELSE
            DO 80, J = 1, N
               IF( X( JX ).NE.ZERO )THEN
                  TEMP = ALPHA*X( JX )
                  IY   = KY
                  DO 70, I = 1, M
                     Y( IY ) = Y( IY ) + TEMP*A( I, J )
                     IY      = IY      + INCY
   70             CONTINUE
               END IF
               JX = JX + INCX
   80       CONTINUE
         END IF
      ELSE
*
*        Form  y := alpha*A'*x + y  or  y := alpha*conjg( A' )*x + y.
*
         JY = KY
         IF( INCX.EQ.1 )THEN
            DO 110, J = 1, N
               TEMP = ZERO
               IF( NOCONJ )THEN
                  DO 90, I = 1, M
                     TEMP = TEMP + A( I, J )*X( I )
   90             CONTINUE
               ELSE
                  DO 100, I = 1, M
                     TEMP = TEMP + DCONJG( A( I, J ) )*X( I )
  100             CONTINUE
               END IF
               Y( JY ) = Y( JY ) + ALPHA*TEMP
               JY      = JY      + INCY
  110       CONTINUE
         ELSE
            DO 140, J = 1, N
               TEMP = ZERO
               IX   = KX
               IF( NOCONJ )THEN
                  DO 120, I = 1, M
                     TEMP = TEMP + A( I, J )*X( IX )
                     IX   = IX   + INCX
  120             CONTINUE
               ELSE
                  DO 130, I = 1, M
                     TEMP = TEMP + DCONJG( A( I, J ) )*X( IX )
                     IX   = IX   + INCX
  130             CONTINUE
               END IF
               Y( JY ) = Y( JY ) + ALPHA*TEMP
               JY      = JY      + INCY
  140       CONTINUE
         END IF
      END IF
*
      RETURN
      END
      subroutine zgesl ( a, lda, n, ipvt, b, job )

c*********************************************************************72
c
cc ZGESL solves a factored complex general linear system.
c
      implicit none

      integer lda,n,ipvt(1),job
      complex*16 a(lda,1),b(1)
c
c     zgesl solves the complex*16 system
c     a * x = b  or  ctrans(a) * x = b
c     using the factors computed by zgeco or zgefa.
c
c     on entry
c
c        a       complex*16(lda, n)
c                the output from zgeco or zgefa.
c
c        lda     integer
c                the leading dimension of the array  a .
c
c        n       integer
c                the order of the matrix  a .
c
c        ipvt    integer(n)
c                the pivot vector from zgeco or zgefa.
c
c        b       complex*16(n)
c                the right hand side vector.
c
c        job     integer
c                = 0         to solve  a*x = b ,
c                = nonzero   to solve  ctrans(a)*x = b  where
c                            ctrans(a)  is the conjugate transpose.
c
c     on return
c
c        b       the solution vector  x .
c
c     error condition
c
c        a division by zero will occur if the input factor contains a
c        zero on the diagonal.  technically this indicates singularity
c        but it is often caused by improper arguments or improper
c        setting of lda .  it will not occur if the subroutines are
c        called correctly and if zgeco has set rcond .gt. 0.0
c        or zgefa has set info .eq. 0 .
c
c     to compute  inverse(a) * c  where  c  is a matrix
c     with  p  columns
c           call zgeco(a,lda,n,ipvt,rcond,z)
c           if (rcond is too small) go to ...
c           do 10 j = 1, p
c              call zgesl(a,lda,n,ipvt,c(1,j),0)
c        10 continue
c
c     linpack. this version dated 08/14/78 .
c     cleve moler, university of new mexico, argonne national lab.
c
c     subroutines and functions
c
c     blas zaxpy,zdotc
c     fortran dconjg
c
c     internal variables
c
      complex*16 zdotc,t
      integer k,kb,l,nm1
      double precision dreal,dimag
      complex*16 zdumr,zdumi
      dreal(zdumr) = zdumr
      dimag(zdumi) = (0.0d0,-1.0d0)*zdumi
c
      nm1 = n - 1
      if (job .ne. 0) go to 50
c
c        job = 0 , solve  a * x = b
c        first solve  l*y = b
c
         if (nm1 .lt. 1) go to 30
         do 20 k = 1, nm1
            l = ipvt(k)
            t = b(l)
            if (l .eq. k) go to 10
               b(l) = b(k)
               b(k) = t
   10       continue
            call zaxpy(n-k,t,a(k+1,k),1,b(k+1),1)
   20    continue
   30    continue
c
c        now solve  u*x = y
c
         do 40 kb = 1, n
            k = n + 1 - kb
            b(k) = b(k)/a(k,k)
            t = -b(k)
            call zaxpy(k-1,t,a(1,k),1,b(1),1)
   40    continue
      go to 100
   50 continue
c
c        job = nonzero, solve  ctrans(a) * x = b
c        first solve  ctrans(u)*y = b
c
         do 60 k = 1, n
            t = zdotc(k-1,a(1,k),1,b(1),1)
            b(k) = (b(k) - t)/dconjg(a(k,k))
   60    continue
c
c        now solve ctrans(l)*x = y
c
         if (nm1 .lt. 1) go to 90
         do 80 kb = 1, nm1
            k = n - kb
            b(k) = b(k) + zdotc(n-k,a(k+1,k),1,b(k+1),1)
            l = ipvt(k)
            if (l .eq. k) go to 70
               t = b(l)
               b(l) = b(k)
               b(k) = t
   70       continue
   80    continue
   90    continue
  100 continue
      return
      end
      subroutine zgesv ( n, m, a,lda, ipiv, b,ldb, iflag )

c*********************************************************************72
c
cc ZGESV solves multiple general linear systems of the form A*X=B.
c
      implicit none

      integer N, M, LDA, LDB, IPIV(N), IFLAG
      complex*16 A(LDA,N), B(LDB,M)
      integer j

      call ZGEFA( A,LDA, N, IPIV, IFLAG )
      if ( IFLAG.ne.0 ) stop "Error in ZGESV (LU factorisation)"
      do j = 1,M
         call ZGESL( A,LDA, N, IPIV,B(1,j), 0 )
      end do
      end
      subroutine zgexpv ( n, m, t, v, w, tol, anorm,
     .                   wsp,lwsp, iwsp,liwsp, matvec, itrace,iflag )

c*********************************************************************72
c
cc ZGEXPV computes w = exp(t*A)*v.
c
      implicit none

      integer          n, m, lwsp, liwsp, itrace, iflag, iwsp(liwsp)
      double precision t, tol, anorm
      complex*16       v(n), w(n), wsp(lwsp)
      external         matvec

*-----Purpose----------------------------------------------------------|
*
*---  ZGEXPV computes w = exp(t*A)*v
*     for a Zomplex (i.e., complex double precision) matrix A
*
*     It does not compute the matrix exponential in isolation but
*     instead, it computes directly the action of the exponential
*     operator on the operand vector. This way of doing so allows
*     for addressing large sparse problems.
*
*     The method used is based on Krylov subspace projection
*     techniques and the matrix under consideration interacts only
*     via the external routine `matvec' performing the matrix-vector
*     product (matrix-free method).
*
*-----Arguments--------------------------------------------------------|
*
*     n      : (input) order of the principal matrix A.
*
*     m      : (input) maximum size for the Krylov basis.
*
*     t      : (input) time at wich the solution is needed (can be < 0).
*
*     v(n)   : (input) given operand vector.
*
*     w(n)   : (output) computed approximation of exp(t*A)*v.
*
*     tol    : (input/output) the requested accuracy tolerance on w.
*              If on input tol=0.0d0 or tol is too small (tol.le.eps)
*              the internal value sqrt(eps) is used, and tol is set to
*              sqrt(eps) on output (`eps' denotes the machine epsilon).
*              (`Happy breakdown' is assumed if h(j+1,j) .le. anorm*tol)
*
*     anorm  : (input) an approximation of some norm of A.
*
*   wsp(lwsp): (workspace) lwsp .ge. n*(m+1)+n+(m+2)^2+4*(m+2)^2+ideg+1
*                                   +---------+-------+---------------+
*              (actually, ideg=6)        V        H      wsp for PADE
*
* iwsp(liwsp): (workspace) liwsp .ge. m+2
*
*     matvec : external routine for matrix-vector multiplication.
*              synopsis: matvec( x, y )
*                        complex*16 x(*), y(*)
*              computes: y(1:n) <- A*x(1:n)
*                        where A is the principal matrix.
*
*     itrace : (input) running mode. 0=silent, 1=print step-by-step info
*
*     iflag  : (output) exit flag.
*              <0 - bad input arguments
*               0 - no problem
*               1 - maximum number of steps reached without convergence
*               2 - requested tolerance was too high
*
*-----Accounts on the computation--------------------------------------|
*     Upon exit, an interested user may retrieve accounts on the
*     computations. They are located in the workspace arrays wsp and
*     iwsp as indicated below:
*
*     location  mnemonic                 description
*     -----------------------------------------------------------------|
*     iwsp(1) = nmult, number of matrix-vector multiplications used
*     iwsp(2) = nexph, number of Hessenberg matrix exponential evaluated
*     iwsp(3) = nscale, number of repeated squaring involved in Pade
*     iwsp(4) = nstep, number of integration steps used up to completion
*     iwsp(5) = nreject, number of rejected step-sizes
*     iwsp(6) = ibrkflag, set to 1 if `happy breakdown' and 0 otherwise
*     iwsp(7) = mbrkdwn, if `happy brkdown', basis-size when it occured
*     -----------------------------------------------------------------|
*     wsp(1)  = step_min, minimum step-size used during integration
*     wsp(2)  = step_max, maximum step-size used during integration
*     wsp(3)  = x_round, maximum among all roundoff errors (lower bound)
*     wsp(4)  = s_round, sum of roundoff errors (lower bound)
*     wsp(5)  = x_error, maximum among all local truncation errors
*     wsp(6)  = s_error, global sum of local truncation errors
*     wsp(7)  = tbrkdwn, if `happy breakdown', time when it occured
*     wsp(8)  = t_now, integration domain successfully covered
*     wsp(9)  = hump, i.e., max||exp(sA)||, s in [0,t] (or [t,0] if t<0)
*     wsp(10) = ||w||/||v||, scaled norm of the solution w.
*     -----------------------------------------------------------------|
*     The `hump' is a measure of the conditioning of the problem. The
*     matrix exponential is well-conditioned if hump = 1, whereas it is
*     poorly-conditioned if hump >> 1. However the solution can still be
*     relatively fairly accurate even when the hump is large (the hump
*     is an upper bound), especially when the hump and the scaled norm
*     of w [this is also computed and returned in wsp(10)] are of the
*     same order of magnitude (further details in reference below).
*
*----------------------------------------------------------------------|
*-----The following parameters may also be adjusted herein-------------|
*
      integer mxstep, mxreject, ideg
      double precision delta, gamma
      parameter( mxstep   = 500,
     .           mxreject = 0,
     .           ideg     = 6,
     .           delta    = 1.2d0,
     .           gamma    = 0.9d0 )

*     mxstep  : maximum allowable number of integration steps.
*               The value 0 means an infinite number of steps.
*
*     mxreject: maximum allowable number of rejections at each step.
*               The value 0 means an infinite number of rejections.
*
*     ideg    : the Pade approximation of type (ideg,ideg) is used as
*               an approximation to exp(H). The value 0 switches to the
*               uniform rational Chebyshev approximation of type (14,14)
*
*     delta   : local truncation error `safety factor'
*
*     gamma   : stepsize `shrinking factor'
*
*----------------------------------------------------------------------|
*     Roger B. Sidje (rbs@maths.uq.edu.au)
*     EXPOKIT: Software Package for Computing Matrix Exponentials.
*     ACM - Transactions On Mathematical Software, 24(1):130-156, 1998
*----------------------------------------------------------------------|

      complex*16 ZERO, ONE
      parameter( ZERO=(0.0d0,0.0d0), ONE=(1.0d0,0.0d0) )
      integer i, j, k1, mh, mx, iv, ih, j1v, ns, ifree, lfree, iexph,
     .        ireject,ibrkflag,mbrkdwn, nmult, nreject, nexph, nscale,
     .        nstep
      double precision sgn, t_out, tbrkdwn, step_min,step_max, err_loc,
     .                 s_error, x_error, t_now, t_new, t_step, t_old,
     .                 xm, beta, break_tol, p1, p2, p3, eps, rndoff,
     .                 vnorm, avnorm, hj1j, hump, SQR1
      complex*16 hij

      intrinsic AINT,ABS,CMPLX,DBLE,INT,LOG10,MAX,MIN,NINT,SIGN,SQRT
      complex*16 ZDOTC
      double precision DZNRM2
*
*---  check restrictions on input parameters ...
*
      iflag = 0
      if ( lwsp.lt.n*(m+2)+5*(m+2)**2+ideg+1 ) iflag = -1
      if ( liwsp.lt.m+2 ) iflag = -2
      if ( m.ge.n .or. m.le.0 ) iflag = -3
      if ( iflag.ne.0 ) stop 'bad sizes (in input of ZGEXPV)'
*
*---  initialisations ...
*
      k1 = 2
      mh = m + 2
      iv = 1
      ih = iv + n*(m+1) + n
      ifree = ih + mh*mh
      lfree = lwsp - ifree + 1

      ibrkflag = 0
      mbrkdwn  = m
      nmult    = 0
      nreject  = 0
      nexph    = 0
      nscale   = 0

      t_out    = ABS( t )
      tbrkdwn  = 0.0d0
      step_min = t_out
      step_max = 0.0d0
      nstep    = 0
      s_error  = 0.0d0
      x_error  = 0.0d0
      t_now    = 0.0d0
      t_new    = 0.0d0

      p1 = 4.0d0/3.0d0
 1    p2 = p1 - 1.0d0
      p3 = p2 + p2 + p2
      eps = ABS( p3-1.0d0 )
      if ( eps.eq.0.0d0 ) go to 1
      if ( tol.le.eps ) tol = SQRT( eps )
      rndoff = eps*anorm

      break_tol = 1.0d-7
*>>>  break_tol = tol
*>>>  break_tol = anorm*tol

      sgn = SIGN( 1.0d0,t )
      call ZCOPY( n, v,1, w,1 )
      beta = DZNRM2( n, w,1 )
      vnorm = beta
      hump = beta
*
*---  obtain the very first stepsize ...
*
      SQR1 = SQRT( 0.1d0 )
      xm = 1.0d0/DBLE( m )
      p2 = tol*(((m+1)/2.72D0)**(m+1))*SQRT(2.0D0*3.14D0*(m+1))
      t_new = (1.0d0/anorm)*(p2/(4.0d0*beta*anorm))**xm
      p1 = 10.0d0**(NINT( LOG10( t_new )-SQR1 )-1)
      t_new = AINT( t_new/p1 + 0.55d0 ) * p1
*
*---  step-by-step integration ...
*
 100  if ( t_now.ge.t_out ) goto 500

      nstep = nstep + 1
      t_step = MIN( t_out-t_now, t_new )
      p1 = 1.0d0/beta
      do i = 1,n
         wsp(iv + i-1) = p1*w(i)
      end do
      do i = 1,mh*mh
         wsp(ih+i-1) = ZERO
      end do
*
*---  Arnoldi loop ...
*
      j1v = iv + n
      do 200 j = 1,m
         nmult = nmult + 1
         call matvec( wsp(j1v-n), wsp(j1v) )
         do i = 1,j
            hij = ZDOTC( n, wsp(iv+(i-1)*n),1, wsp(j1v),1 )
            call ZAXPY( n, -hij, wsp(iv+(i-1)*n),1, wsp(j1v),1 )
            wsp(ih+(j-1)*mh+i-1) = hij
         end do
         hj1j = DZNRM2( n, wsp(j1v),1 )
*---     if `happy breakdown' go straightforward at the end ...
         if ( hj1j.le.break_tol ) then
            print*,'happy breakdown: mbrkdwn =',j,' h =',hj1j
            k1 = 0
            ibrkflag = 1
            mbrkdwn = j
            tbrkdwn = t_now
            t_step = t_out-t_now
            goto 300
         endif
         wsp(ih+(j-1)*mh+j) = CMPLX( hj1j )
         call ZDSCAL( n, 1.0d0/hj1j, wsp(j1v),1 )
         j1v = j1v + n
 200  continue
      nmult = nmult + 1
      call matvec( wsp(j1v-n), wsp(j1v) )
      avnorm = DZNRM2( n, wsp(j1v),1 )
*
*---  set 1 for the 2-corrected scheme ...
*
 300  continue
      wsp(ih+m*mh+m+1) = ONE
*
*---  loop while ireject<mxreject until the tolerance is reached ...
*
      ireject = 0
 401  continue
*
*---  compute w = beta*V*exp(t_step*H)*e1 ...
*
      nexph = nexph + 1
      mx = mbrkdwn + k1
      if ( ideg.ne.0 ) then
*---     irreducible rational Pade approximation ...
         call ZGPADM( ideg, mx, sgn*t_step, wsp(ih),mh,
     .                wsp(ifree),lfree, iwsp, iexph, ns, iflag )
         iexph = ifree + iexph - 1
         nscale = nscale + ns
      else
*---     uniform rational Chebyshev approximation ...
         iexph = ifree
         do i = 1,mx
            wsp(iexph+i-1) = ZERO
         end do
         wsp(iexph) = ONE
         call ZNCHBV(mx,sgn*t_step,wsp(ih),mh,wsp(iexph),wsp(ifree+mx))
      endif
 402  continue
*
*---  error estimate ...
*
      if ( k1.eq.0 ) then
         err_loc = tol
      else
         p1 = ABS( wsp(iexph+m) )   * beta
         p2 = ABS( wsp(iexph+m+1) ) * beta * avnorm
         if ( p1.gt.10.0d0*p2 ) then
            err_loc = p2
            xm = 1.0d0/DBLE( m )
         elseif ( p1.gt.p2 ) then
            err_loc = (p1*p2)/(p1-p2)
            xm = 1.0d0/DBLE( m )
         else
            err_loc = p1
            xm = 1.0d0/DBLE( m-1 )
         endif
      endif
*
*---  reject the step-size if the error is not acceptable ...
*
      if ( (k1.ne.0) .and. (err_loc.gt.delta*t_step*tol) .and.
     .     (mxreject.eq.0 .or. ireject.lt.mxreject) ) then
         t_old = t_step
         t_step = gamma * t_step * (t_step*tol/err_loc)**xm
         p1 = 10.0d0**(NINT( LOG10( t_step )-SQR1 )-1)
         t_step = AINT( t_step/p1 + 0.55d0 ) * p1
         if ( itrace.ne.0 ) then
            print*,'t_step =',t_old
            print*,'err_loc =',err_loc
            print*,'err_required =',delta*t_old*tol
            print*,'stepsize rejected, stepping down to:',t_step
         endif
         ireject = ireject + 1
         nreject = nreject + 1
         if ( mxreject.ne.0 .and. ireject.gt.mxreject ) then
            print*,"Failure in ZGEXPV: ---"
            print*,"The requested tolerance is too high."
            Print*,"Rerun with a smaller value."
            iflag = 2
            return
         endif
         goto 401
      endif
*
*---  now update w = beta*V*exp(t_step*H)*e1 and the hump ...
*
      mx = mbrkdwn + MAX( 0,k1-1 )
      hij = CMPLX( beta )
      call ZGEMV( 'n', n,mx,hij,wsp(iv),n,wsp(iexph),1,ZERO,w,1 )
      beta = DZNRM2( n, w,1 )
      hump = MAX( hump, beta )
*
*---  suggested value for the next stepsize ...
*
      t_new = gamma * t_step * (t_step*tol/err_loc)**xm
      p1 = 10.0d0**(NINT( LOG10( t_new )-SQR1 )-1)
      t_new = AINT( t_new/p1 + 0.55d0 ) * p1

      err_loc = MAX( err_loc,rndoff )
*
*---  update the time covered ...
*
      t_now = t_now + t_step
*
*---  display and keep some information ...
*
      if ( itrace.ne.0 ) then
         print*,'integration',nstep,'---------------------------------'
         print*,'scale-square =',ns
         print*,'step_size =',t_step
         print*,'err_loc   =',err_loc
         print*,'next_step =',t_new
      endif

      step_min = MIN( step_min, t_step )
      step_max = MAX( step_max, t_step )
      s_error = s_error + err_loc
      x_error = MAX( x_error, err_loc )

      if ( mxstep.eq.0 .or. nstep.lt.mxstep ) goto 100
      iflag = 1

 500  continue

      iwsp(1) = nmult
      iwsp(2) = nexph
      iwsp(3) = nscale
      iwsp(4) = nstep
      iwsp(5) = nreject
      iwsp(6) = ibrkflag
      iwsp(7) = mbrkdwn

      wsp(1)  = CMPLX( step_min )
      wsp(2)  = CMPLX( step_max )
      wsp(3)  = CMPLX( 0.0d0 )
      wsp(4)  = CMPLX( 0.0d0 )
      wsp(5)  = CMPLX( x_error )
      wsp(6)  = CMPLX( s_error )
      wsp(7)  = CMPLX( tbrkdwn )
      wsp(8)  = CMPLX( sgn*t_now )
      wsp(9)  = CMPLX( hump/vnorm )
      wsp(10) = CMPLX( beta/vnorm )
      end
      subroutine zgpadm ( ideg, m, t, h, ldh, wsp, lwsp, ipiv, iexph,
     &  ns, iflag )

c*********************************************************************72
c
cc ZGPADM computes exp(t*H), the matrix exponential of a general complex matrix.
c
      implicit none

      double precision t
      integer          ideg, m, ldh, lwsp, iexph, ns, iflag, ipiv(m)
      complex*16       H(ldh,m), wsp(lwsp)

*-----Purpose----------------------------------------------------------|
*
*     Computes exp(t*H), the matrix exponential of a general complex
*     matrix in full, using the irreducible rational Pade approximation
*     to the exponential exp(z) = r(z) = (+/-)( I + 2*(q(z)/p(z)) ),
*     combined with scaling-and-squaring.
*
*-----Arguments--------------------------------------------------------|
*
*     ideg      : (input) the degre of the diagonal Pade to be used.
*                 a value of 6 is generally satisfactory.
*
*     m         : (input) order of H.
*
*     H(ldh,m)  : (input) argument matrix.
*
*     t         : (input) time-scale (can be < 0).
*
*     wsp(lwsp) : (workspace/output) lwsp .ge. 4*m*m+ideg+1.
*
*     ipiv(m)   : (workspace)
*
*>>>> iexph     : (output) number such that wsp(iexph) points to exp(tH)
*                 i.e., exp(tH) is located at wsp(iexph ... iexph+m*m-1)
*                       ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
*                 NOTE: if the routine was called with wsp(iptr),
*                       then exp(tH) will start at wsp(iptr+iexph-1).
*
*     ns        : (output) number of scaling-squaring used.
*
*     iflag     : (output) exit flag.
*                       0 - no problem
*                      <0 - problem
*
*----------------------------------------------------------------------|
*     Roger B. Sidje (rbs@maths.uq.edu.au)
*     EXPOKIT: Software Package for Computing Matrix Exponentials.
*     ACM - Transactions On Mathematical Software, 24(1):130-156, 1998
*----------------------------------------------------------------------|
*
      integer i,j,k,icoef,mm,ih2,iodd,iused,ifree,iq,ip,iput,iget
      double precision hnorm
      complex*16 cp, cq, scale, scale2, ZERO, ONE

      parameter( ZERO=(0.0d0,0.0d0), ONE=(1.0d0,0.0d0) )
      intrinsic ABS, CMPLX, DBLE, INT, LOG, MAX

*---  check restrictions on input parameters ...
      mm = m*m
      iflag = 0
      if ( ldh.lt.m ) iflag = -1
      if ( lwsp.lt.4*mm+ideg+1 ) iflag = -2
      if ( iflag.ne.0 ) stop 'bad sizes (in input of ZGPADM)'
*
*---  initialise pointers ...
*
      icoef = 1
      ih2 = icoef + (ideg+1)
      ip  = ih2 + mm
      iq  = ip + mm
      ifree = iq + mm
*
*---  scaling: seek ns such that ||t*H/2^ns|| < 1/2;
*     and set scale = t/2^ns ...
*
      do i = 1,m
         wsp(i) = ZERO
      end do
      do j = 1,m
         do i = 1,m
            wsp(i) = wsp(i) + ABS( H(i,j) )
         end do
      end do
      hnorm = 0.0d0
      do i = 1,m
         hnorm = MAX( hnorm,DBLE(wsp(i)) )
      end do
      hnorm = ABS( t*hnorm )
      if ( hnorm.eq.0.0d0 ) stop 'Error - null H in input of ZGPADM.'
      ns = MAX( 0,INT(LOG(hnorm)/LOG(2.0d0))+2 )
      scale =  CMPLX( t/DBLE(2**ns),0.0d0 )
      scale2 = scale*scale
*
*---  compute Pade coefficients ...
*
      i = ideg+1
      j = 2*ideg+1
      wsp(icoef) = ONE
      do k = 1,ideg
         wsp(icoef+k) = (wsp(icoef+k-1)*DBLE( i-k ))/DBLE( k*(j-k) )
      end do
*
*---  H2 = scale2*H*H ...
*
      call ZGEMM( 'n','n',m,m,m,scale2,H,ldh,H,ldh,ZERO,wsp(ih2),m )
*
*---  initialise p (numerator) and q (denominator) ...
*
      cp = wsp(icoef+ideg-1)
      cq = wsp(icoef+ideg)
      do j = 1,m
         do i = 1,m
            wsp(ip + (j-1)*m + i-1) = ZERO
            wsp(iq + (j-1)*m + i-1) = ZERO
         end do
         wsp(ip + (j-1)*(m+1)) = cp
         wsp(iq + (j-1)*(m+1)) = cq
      end do
*
*---  Apply Horner rule ...
*
      iodd = 1
      k = ideg - 1
 100  continue
      iused = iodd*iq + (1-iodd)*ip
      call ZGEMM( 'n','n',m,m,m, ONE,wsp(iused),m,
     .             wsp(ih2),m, ZERO,wsp(ifree),m )
      do j = 1,m
         wsp(ifree+(j-1)*(m+1)) = wsp(ifree+(j-1)*(m+1))+wsp(icoef+k-1)
      end do
      ip = (1-iodd)*ifree + iodd*ip
      iq = iodd*ifree + (1-iodd)*iq
      ifree = iused
      iodd = 1-iodd
      k = k-1
      if ( k.gt.0 )  goto 100
*
*---  Obtain (+/-)(I + 2*(p\q)) ...
*
      if ( iodd.ne.0 ) then
         call ZGEMM( 'n','n',m,m,m, scale,wsp(iq),m,
     .                H,ldh, ZERO,wsp(ifree),m )
         iq = ifree
      else
         call ZGEMM( 'n','n',m,m,m, scale,wsp(ip),m,
     .                H,ldh, ZERO,wsp(ifree),m )
         ip = ifree
      endif
      call ZAXPY( mm, -ONE,wsp(ip),1, wsp(iq),1 )
      call ZGESV( m,m, wsp(iq),m, ipiv, wsp(ip),m, iflag )
      if ( iflag.ne.0 ) stop 'Problem in ZGESV (within ZGPADM)'
      call ZDSCAL( mm, 2.0d0, wsp(ip), 1 )
      do j = 1,m
         wsp(ip+(j-1)*(m+1)) = wsp(ip+(j-1)*(m+1)) + ONE
      end do
      iput = ip
      if ( ns.eq.0 .and. iodd.ne.0 ) then
         call ZDSCAL( mm, -1.0d0, wsp(ip), 1 )
         goto 200
      endif
*
*--   squaring : exp(t*H) = (exp(t*H))^(2^ns) ...
*
      iodd = 1
      do k = 1,ns
         iget = iodd*ip + (1-iodd)*iq
         iput = (1-iodd)*ip + iodd*iq
         call ZGEMM( 'n','n',m,m,m, ONE,wsp(iget),m, wsp(iget),m,
     .                ZERO,wsp(iput),m )
         iodd = 1-iodd
      end do
 200  continue
      iexph = iput
      end
      subroutine zgphiv ( n, m, t, u, v, w, tol, anorm,
     .                   wsp,lwsp, iwsp,liwsp, matvec, itrace,iflag )

c*********************************************************************72
c
cc ZGPHIV solves the ODE w' = Aw + u, w(0) = v.
c
      implicit none

      integer          n, m, lwsp, liwsp, itrace, iflag, iwsp(liwsp)
      double precision t, tol, anorm
      complex*16       u(n), v(n), w(n), wsp(lwsp)
      external         matvec

*-----Purpose----------------------------------------------------------|
*
*---  ZGPHIV computes w = exp(t*A)v + t*phi(tA)u which is the solution
*     of the nonhomogeneous linear ODE problem w' = Aw + u, w(0) = v.
*     phi(z) = (exp(z)-1)/z and A is a Zomplex (i.e., complex double
*     precision matrix).
*
*     The method used is based on Krylov subspace projection
*     techniques and the matrix under consideration interacts only
*     via the external routine `matvec' performing the matrix-vector
*     product (matrix-free method).
*
*-----Arguments--------------------------------------------------------|
*
*     n      : (input) order of the principal matrix A.
*
*     m      : (input) maximum size for the Krylov basis.
*
*     t      : (input) time at wich the solution is needed (can be < 0).
*
*     u(n)   : (input) operand vector with respect to the phi function
*              (forcing term of the ODE problem).
*
*     v(n)   : (input) operand vector with respect to the exp function
*              (initial condition of the ODE problem).
*
*     w(n)   : (output) computed approximation of exp(t*A)v + t*phi(tA)u
*
*     tol    : (input/output) the requested accuracy tolerance on w.
*              If on input tol=0.0d0 or tol is too small (tol.le.eps)
*              the internal value sqrt(eps) is used, and tol is set to
*              sqrt(eps) on output (`eps' denotes the machine epsilon).
*              (`Happy breakdown' is assumed if h(j+1,j) .le. anorm*tol)
*
*     anorm  : (input) an approximation of some norm of A.
*
*   wsp(lwsp): (workspace) lwsp .ge. n*(m+1)+n+(m+3)^2+4*(m+3)^2+ideg+1
*                                   +---------+-------+---------------+
*              (actually, ideg=6)        V        H      wsp for PADE
*
* iwsp(liwsp): (workspace) liwsp .ge. m+3
*
*     matvec : external routine for matrix-vector multiplication.
*              synopsis: matvec( x, y )
*                        double precision x(*), y(*)
*              computes: y(1:n) <- A*x(1:n)
*                        where A is the principal matrix.
*
*     itrace : (input) running mode. 0=silent, 1=print step-by-step info
*
*     iflag  : (output) exit flag.
*              <0 - bad input arguments
*               0 - no problem
*               1 - maximum number of steps reached without convergence
*               2 - requested tolerance was too high
*
*-----Accounts on the computation--------------------------------------|
*     Upon exit, an interested user may retrieve accounts on the
*     computations. They are located in the workspace arrays wsp and
*     iwsp as indicated below:
*
*     location  mnemonic                 description
*     -----------------------------------------------------------------|
*     iwsp(1) = nmult, number of matrix-vector multiplications used
*     iwsp(2) = nexph, number of Hessenberg matrix exponential evaluated
*     iwsp(3) = nscale, number of repeated squaring involved in Pade
*     iwsp(4) = nstep, number of integration steps used up to completion
*     iwsp(5) = nreject, number of rejected step-sizes
*     iwsp(6) = ibrkflag, set to 1 if `happy breakdown' and 0 otherwise
*     iwsp(7) = mbrkdwn, if `happy brkdown', basis-size when it occured
*     -----------------------------------------------------------------|
*     wsp(1)  = step_min, minimum step-size used during integration
*     wsp(2)  = step_max, maximum step-size used during integration
*     wsp(3)  = dummy
*     wsp(4)  = dummy
*     wsp(5)  = x_error, maximum among all local truncation errors
*     wsp(6)  = s_error, global sum of local truncation errors
*     wsp(7)  = tbrkdwn, if `happy breakdown', time when it occured
*     wsp(8)  = t_now, integration domain successfully covered
*
*----------------------------------------------------------------------|
*-----The following parameters may also be adjusted herein-------------|
*
      integer mxstep, mxreject, ideg
      double precision delta, gamma
      parameter( mxstep   = 500,
     .           mxreject = 0,
     .           ideg     = 6,
     .           delta    = 1.2d0,
     .           gamma    = 0.9d0 )

*     mxstep  : maximum allowable number of integration steps.
*               The value 0 means an infinite number of steps.
*
*     mxreject: maximum allowable number of rejections at each step.
*               The value 0 means an infinite number of rejections.
*
*     ideg    : the Pade approximation of type (ideg,ideg) is used as
*               an approximation to exp(H).
*
*     delta   : local truncation error `safety factor'
*
*     gamma   : stepsize `shrinking factor'
*
*----------------------------------------------------------------------|
*     Roger B. Sidje (rbs@maths.uq.edu.au)
*     EXPOKIT: Software Package for Computing Matrix Exponentials.
*     ACM - Transactions On Mathematical Software, 24(1):130-156, 1998
*----------------------------------------------------------------------|
*
      complex*16 ZERO, ONE
      parameter( ZERO=(0.0d0,0.0d0), ONE=(1.0d0,0.0d0) )
      integer i, j, k1, mh, mx, iv, ih, j1v, ns, ifree, lfree, iexph,
     .        ireject,ibrkflag,mbrkdwn, nmult, nreject, nexph, nscale,
     .        nstep, iphih
      double precision sgn, t_out, tbrkdwn, step_min,step_max, err_loc,
     .                 s_error, x_error, t_now, t_new, t_step, t_old,
     .                 xm, beta, break_tol, p1, p2, p3, eps, rndoff,
     .                 avnorm, hj1j, SQR1
      complex*16 hij

      intrinsic AINT,ABS,CMPLX,DBLE,INT,LOG10,MAX,MIN,NINT,SIGN,SQRT
      complex*16 ZDOTC
      double precision DZNRM2
*
*---  check restrictions on input parameters ...
*
      iflag = 0
      if ( lwsp.lt.n*(m+3)+5*(m+3)**2+ideg+1 ) iflag = -1
      if ( liwsp.lt.m+3 ) iflag = -2
      if ( m.ge.n .or. m.le.0 ) iflag = -3
      if ( iflag.ne.0 ) stop 'bad sizes (in input of ZGPHIV)'
*
*---  initialisations ...
*
      k1 = 3
      mh = m + 3
      iv = 1
      ih = iv + n*(m+1) + n
      ifree = ih + mh*mh
      lfree = lwsp - ifree + 1

      ibrkflag = 0
      mbrkdwn  = m
      nmult    = 0
      nreject  = 0
      nexph    = 0
      nscale   = 0

      t_out    = ABS( t )
      tbrkdwn  = 0.0d0
      step_min = t_out
      step_max = 0.0d0
      nstep    = 0
      s_error  = 0.0d0
      x_error  = 0.0d0
      t_now    = 0.0d0
      t_new    = 0.0d0

      p1 = 4.0d0/3.0d0
 1    p2 = p1 - 1.0d0
      p3 = p2 + p2 + p2
      eps = ABS( p3-1.0d0 )
      if ( eps.eq.0.0d0 ) go to 1
      if ( tol.le.eps ) tol = SQRT( eps )
      rndoff = eps*anorm

      break_tol = 1.0d-7
*>>>  break_tol = tol
*>>>  break_tol = anorm*tol

      sgn = SIGN( 1.0d0,t )
      SQR1 = SQRT( 0.1d0 )
      call ZCOPY( n, v,1, w,1 )
*
*---  step-by-step integration ...
*
 100  if ( t_now.ge.t_out ) goto 500

      nmult =  nmult + 1
      call matvec( w, wsp(iv) )
      call ZAXPY( n, ONE, u,1, wsp(iv),1 )
      beta = DZNRM2( n, wsp(iv),1 )
      if ( beta.eq.0.0d0 ) goto 500
      call ZDSCAL( n, 1.0d0/beta, wsp(iv),1 )
      do i = 1,mh*mh
         wsp(ih+i-1) = ZERO
      end do

      if ( nstep.eq.0 ) then
*---     obtain the very first stepsize ...
         xm = 1.0d0/DBLE( m )
         p1 = tol*(((m+1)/2.72D0)**(m+1))*SQRT(2.0D0*3.14D0*(m+1))
         t_new = (1.0d0/anorm)*(p1/(4.0d0*beta*anorm))**xm
         p1 = 10.0d0**(NINT( LOG10( t_new )-SQR1 )-1)
         t_new = AINT( t_new/p1 + 0.55d0 ) * p1
      endif
      nstep = nstep + 1
      t_step = MIN( t_out-t_now, t_new )
*
*---  Arnoldi loop ...
*
      j1v = iv + n
      do 200 j = 1,m
         nmult = nmult + 1
         call matvec( wsp(j1v-n), wsp(j1v) )
         do i = 1,j
            hij = ZDOTC( n, wsp(iv+(i-1)*n),1, wsp(j1v),1 )
            call ZAXPY( n, -hij, wsp(iv+(i-1)*n),1, wsp(j1v),1 )
            wsp(ih+(j-1)*mh+i-1) = hij
         end do
         hj1j = DZNRM2( n, wsp(j1v),1 )
*---     if `happy breakdown' go straightforward at the end ...
         if ( hj1j.le.break_tol ) then
            print*,'happy breakdown: mbrkdwn =',j,' h =',hj1j
            k1 = 0
            ibrkflag = 1
            mbrkdwn = j
            tbrkdwn = t_now
            t_step = t_out-t_now
            goto 300
         endif
         wsp(ih+(j-1)*mh+j) = CMPLX( hj1j )
         call ZDSCAL( n, 1.0d0/hj1j, wsp(j1v),1 )
         j1v = j1v + n
 200  continue
      nmult = nmult + 1
      call matvec( wsp(j1v-n), wsp(j1v) )
      avnorm = DZNRM2( n, wsp(j1v),1 )
*
*---  set 1's for the 3-extended scheme ...
*
 300  continue
      wsp(ih+mh*mbrkdwn) = ONE
      wsp(ih+(m-1)*mh+m) = ZERO
      do i = 1,k1-1
         wsp(ih+(m+i)*mh+m+i-1) = ONE
      end do
*
*---  loop while ireject<mxreject until the tolerance is reached ...
*
      ireject = 0
 401  continue
*
*---  compute w = beta*t_step*V*phi(t_step*H)*e1 + w
*
      nexph = nexph + 1
      mx = mbrkdwn + MAX( 1,k1 )

*---  irreducible rational Pade approximation ...
      call ZGPADM( ideg, mx, sgn*t_step, wsp(ih),mh,
     .             wsp(ifree),lfree, iwsp, iexph, ns, iflag )
      iexph = ifree + iexph - 1
      iphih = iexph + mbrkdwn*mx
      nscale = nscale + ns
      wsp(iphih+mbrkdwn)   = hj1j*wsp(iphih+mx+mbrkdwn-1)
      wsp(iphih+mbrkdwn+1) = hj1j*wsp(iphih+2*mx+mbrkdwn-1)

 402  continue
*
*---  error estimate ...
*
      if ( k1.eq.0 ) then
         err_loc = tol
      else
         p1 = ABS( wsp(iphih+m) )   * beta
         p2 = ABS( wsp(iphih+m+1) ) * beta * avnorm
         if ( p1.gt.10.0d0*p2 ) then
            err_loc = p2
            xm = 1.0d0/DBLE( m+1 )
         elseif ( p1.gt.p2 ) then
            err_loc = (p1*p2)/(p1-p2)
            xm = 1.0d0/DBLE( m+1 )
         else
            err_loc = p1
            xm = 1.0d0/DBLE( m )
         endif
      endif
*
*---  reject the step-size if the error is not acceptable ...
*
      if ( (k1.ne.0) .and. (err_loc.gt.delta*t_step*tol) .and.
     .     (mxreject.eq.0 .or. ireject.lt.mxreject) ) then
         t_old = t_step
         t_step = gamma * t_step * (t_step*tol/err_loc)**xm
         p1 = 10.0d0**(NINT( LOG10( t_step )-SQR1 )-1)
         t_step = AINT( t_step/p1 + 0.55d0 ) * p1
         if ( itrace.ne.0 ) then
            print*,'t_step =',t_old
            print*,'err_loc =',err_loc
            print*,'err_required =',delta*t_old*tol
            print*,'stepsize rejected, stepping down to:',t_step
         endif
         ireject = ireject + 1
         nreject = nreject + 1
         if ( mxreject.ne.0 .and. ireject.gt.mxreject ) then
            print*,"Failure in ZGPHIV: ---"
            print*,"The requested tolerance is too high."
            Print*,"Rerun with a smaller value."
            iflag = 2
            return
         endif
         goto 401
      endif
*
      mx = mbrkdwn + MAX( 0,k1-2 )
      hij = CMPLX( beta )
      call ZGEMV( 'n', n,mx,hij,wsp(iv),n,wsp(iphih),1,ONE,w,1 )
*
*---  suggested value for the next stepsize ...
*
      t_new = gamma * t_step * (t_step*tol/err_loc)**xm
      p1 = 10.0d0**(NINT( LOG10( t_new )-SQR1 )-1)
      t_new = AINT( t_new/p1 + 0.55d0 ) * p1

      err_loc = MAX( err_loc,rndoff )
*
*---  update the time covered ...
*
      t_now = t_now + t_step
*
*---  display and keep some information ...
*
      if ( itrace.ne.0 ) then
         print*,'integration',nstep,'---------------------------------'
         print*,'scale-square =',ns
         print*,'step_size =',t_step
         print*,'err_loc   =',err_loc
         print*,'next_step =',t_new
      endif

      step_min = MIN( step_min, t_step )
      step_max = MAX( step_max, t_step )
      s_error = s_error + err_loc
      x_error = MAX( x_error, err_loc )

      if ( mxstep.eq.0 .or. nstep.lt.mxstep ) goto 100
      iflag = 1

 500  continue

      iwsp(1) = nmult
      iwsp(2) = nexph
      iwsp(3) = nscale
      iwsp(4) = nstep
      iwsp(5) = nreject
      iwsp(6) = ibrkflag
      iwsp(7) = mbrkdwn

      wsp(1)  = CMPLX( step_min )
      wsp(2)  = CMPLX( step_max )
      wsp(3)  = CMPLX( 0.0d0 )
      wsp(4)  = CMPLX( 0.0d0 )
      wsp(5)  = CMPLX( x_error )
      wsp(6)  = CMPLX( s_error )
      wsp(7)  = CMPLX( tbrkdwn )
      wsp(8)  = CMPLX( sgn*t_now )
      end
      subroutine zhesv ( uplo, n, m, a, lda, ipiv, b, ldb, wrk, lwrk,
     &  iflag )

c*********************************************************************72
c
cc ZHESV solves multiple hermitian linear systems of the form A*X=B.
c
      implicit none

      character UPLO*1
      integer N, M, LDA, LDB, LWRK, IFLAG, IPIV(N)
      complex*16 A(LDA,N), B(LDB,M), WRK(LWRK)
      integer j

      call ZHIFA( A,LDA, N, IPIV, IFLAG )
      if ( IFLAG.ne.0 ) stop "Error in ZHESV (LDL' factorisation)"
      do j = 1,M
         call ZHISL( A,LDA, N, IPIV,B(1,j) )
      end do
      end
      subroutine zhexpv ( n, m, t, v, w, tol, anorm,
     .                   wsp,lwsp, iwsp,liwsp, matvec, itrace,iflag )

c*********************************************************************72
c
cc ZHEXPV computes w = exp(t*A)*v for a Hermitian matrix A.
c
      implicit none

      integer          n, m, lwsp, liwsp, itrace, iflag, iwsp(liwsp)
      double precision t, tol, anorm
      complex*16       v(n), w(n), wsp(lwsp)
      external         matvec

*-----Purpose----------------------------------------------------------|
*
*---  ZHEXPV computes w = exp(t*A)*v for a Hermitian matrix A.
*
*     It does not compute the matrix exponential in isolation but
*     instead, it computes directly the action of the exponential
*     operator on the operand vector. This way of doing so allows
*     for addressing large sparse problems.
*
*     The method used is based on Krylov subspace projection
*     techniques and the matrix under consideration interacts only
*     via the external routine `matvec' performing the matrix-vector
*     product (matrix-free method).
*
*-----Arguments--------------------------------------------------------|
*
*     n      : (input) order of the principal matrix A.
*
*     m      : (input) maximum size for the Krylov basis.
*
*     t      : (input) time at wich the solution is needed (can be < 0).
*
*     v(n)   : (input) given operand vector.
*
*     w(n)   : (output) computed approximation of exp(t*A)*v.
*
*     tol    : (input/output) the requested accuracy tolerance on w.
*              If on input tol=0.0d0 or tol is too small (tol.le.eps)
*              the internal value sqrt(eps) is used, and tol is set to
*              sqrt(eps) on output (`eps' denotes the machine epsilon).
*              (`Happy breakdown' is assumed if h(j+1,j) .le. anorm*tol)
*
*     anorm  : (input) an approximation of some norm of A.
*
*   wsp(lwsp): (workspace) lwsp .ge. n*(m+1)+n+(m+2)^2+4*(m+2)^2+ideg+1
*                                   +---------+-------+---------------+
*              (actually, ideg=6)        V        H      wsp for PADE
*
* iwsp(liwsp): (workspace) liwsp .ge. m+2
*
*     matvec : external routine for matrix-vector multiplication.
*              synopsis: matvec( x, y )
*                        complex*16 x(*), y(*)
*              computes: y(1:n) <- A*x(1:n)
*                        where A is the principal matrix.
*
*     itrace : (input) running mode. 0=silent, 1=print step-by-step info
*
*     iflag  : (output) exit flag.
*              <0 - bad input arguments
*               0 - no problem
*               1 - maximum number of steps reached without convergence
*               2 - requested tolerance was too high
*
*-----Accounts on the computation--------------------------------------|
*     Upon exit, an interested user may retrieve accounts on the
*     computations. They are located in the workspace arrays wsp and
*     iwsp as indicated below:
*
*     location  mnemonic                 description
*     -----------------------------------------------------------------|
*     iwsp(1) = nmult, number of matrix-vector multiplications used
*     iwsp(2) = nexph, number of Hessenberg matrix exponential evaluated
*     iwsp(3) = nscale, number of repeated squaring involved in Pade
*     iwsp(4) = nstep, number of integration steps used up to completion
*     iwsp(5) = nreject, number of rejected step-sizes
*     iwsp(6) = ibrkflag, set to 1 if `happy breakdown' and 0 otherwise
*     iwsp(7) = mbrkdwn, if `happy brkdown', basis-size when it occured
*     -----------------------------------------------------------------|
*     wsp(1)  = step_min, minimum step-size used during integration
*     wsp(2)  = step_max, maximum step-size used during integration
*     wsp(3)  = dummy
*     wsp(4)  = dummy
*     wsp(5)  = x_error, maximum among all local truncation errors
*     wsp(6)  = s_error, global sum of local truncation errors
*     wsp(7)  = tbrkdwn, if `happy breakdown', time when it occured
*     wsp(8)  = t_now, integration domain successfully covered
*     wsp(9)  = hump, i.e., max||exp(sA)||, s in [0,t] (or [t,0] if t<0)
*     wsp(10) = ||w||/||v||, scaled norm of the solution w.
*     -----------------------------------------------------------------|
*     The `hump' is a measure of the conditioning of the problem. The
*     matrix exponential is well-conditioned if hump = 1, whereas it is
*     poorly-conditioned if hump >> 1. However the solution can still be
*     relatively fairly accurate even when the hump is large (the hump
*     is an upper bound), especially when the hump and the scaled norm
*     of w [this is also computed and returned in wsp(10)] are of the
*     same order of magnitude (further details in reference below).
*
*----------------------------------------------------------------------|
*-----The following parameters may also be adjusted herein-------------|
*
      integer mxstep, mxreject, ideg
      double precision delta, gamma
      parameter( mxstep   = 500,
     .           mxreject = 0,
     .           ideg     = 6,
     .           delta    = 1.2d0,
     .           gamma    = 0.9d0 )

*     mxstep  : maximum allowable number of integration steps.
*               The value 0 means an infinite number of steps.
*
*     mxreject: maximum allowable number of rejections at each step.
*               The value 0 means an infinite number of rejections.
*
*     ideg    : the Pade approximation of type (ideg,ideg) is used as
*               an approximation to exp(H). The value 0 switches to the
*               uniform rational Chebyshev approximation of type (14,14)
*
*     delta   : local truncation error `safety factor'
*
*     gamma   : stepsize `shrinking factor'
*
*----------------------------------------------------------------------|
*     Roger B. Sidje (rbs@maths.uq.edu.au)
*     EXPOKIT: Software Package for Computing Matrix Exponentials.
*     ACM - Transactions On Mathematical Software, 24(1):130-156, 1998
*----------------------------------------------------------------------|
*
      complex*16 ZERO, ONE
      parameter( ZERO=(0.0d0,0.0d0), ONE=(1.0d0,0.0d0) )
      integer i, j, k1, mh, mx, iv, ih, j1v, ns, ifree, lfree, iexph,
     .        ireject,ibrkflag,mbrkdwn, nmult, nreject, nexph, nscale,
     .        nstep
      double precision sgn, t_out, tbrkdwn, step_min,step_max, err_loc,
     .                 s_error, x_error, t_now, t_new, t_step, t_old,
     .                 xm, beta, break_tol, p1, p2, p3, eps, rndoff,
     .                 vnorm, avnorm, hj1j, hump, SQR1
      complex*16 hjj

      intrinsic AINT,ABS,CMPLX,DBLE,INT,LOG10,MAX,MIN,NINT,SIGN,SQRT
      complex*16 ZDOTC
      double precision DZNRM2
*
*---  check restrictions on input parameters ...
*
      iflag = 0
      if ( lwsp.lt.n*(m+2)+5*(m+2)**2+ideg+1 ) iflag = -1
      if ( liwsp.lt.m+2 ) iflag = -2
      if ( m.ge.n .or. m.le.0 ) iflag = -3
      if ( iflag.ne.0 ) stop 'bad sizes (in input of DHEXPV)'
*
*---  initialisations ...
*
      k1 = 2
      mh = m + 2
      iv = 1
      ih = iv + n*(m+1) + n
      ifree = ih + mh*mh
      lfree = lwsp - ifree + 1

      ibrkflag = 0
      mbrkdwn  = m
      nmult    = 0
      nreject  = 0
      nexph    = 0
      nscale   = 0

      t_out    = ABS( t )
      tbrkdwn  = 0.0d0
      step_min = t_out
      step_max = 0.0d0
      nstep    = 0
      s_error  = 0.0d0
      x_error  = 0.0d0
      t_now    = 0.0d0
      t_new    = 0.0d0

      p1 = 4.0d0/3.0d0
 1    p2 = p1 - 1.0d0
      p3 = p2 + p2 + p2
      eps = ABS( p3-1.0d0 )
      if ( eps.eq.0.0d0 ) go to 1
      if ( tol.le.eps ) tol = SQRT( eps )
      rndoff = eps*anorm

      break_tol = 1.0d-7
*>>>  break_tol = tol
*>>>  break_tol = anorm*tol

      sgn = SIGN( 1.0d0,t )
      call ZCOPY( n, v,1, w,1 )
      beta = DZNRM2( n, w,1 )
      vnorm = beta
      hump = beta
*
*---  obtain the very first stepsize ...
*
      SQR1 = SQRT( 0.1d0 )
      xm = 1.0d0/DBLE( m )
      p2 = tol*(((m+1)/2.72D0)**(m+1))*SQRT(2.0D0*3.14D0*(m+1))
      t_new = (1.0d0/anorm)*(p2/(4.0d0*beta*anorm))**xm
      p1 = 10.0d0**(NINT( LOG10( t_new )-SQR1 )-1)
      t_new = AINT( t_new/p1 + 0.55d0 ) * p1
*
*---  step-by-step integration ...
*
 100  if ( t_now.ge.t_out ) goto 500

      nstep = nstep + 1
      t_step = MIN( t_out-t_now, t_new )
      beta = DZNRM2( n, w,1 )
      p1 = 1.0d0/beta
      do i = 1,n
         wsp(iv + i-1) = p1*w(i)
      end do
      do i = 1,mh*mh
         wsp(ih+i-1) = ZERO
      end do
*
*---  Lanczos loop ...
*
      j1v = iv + n
      do 200 j = 1,m
         nmult = nmult + 1
         call matvec( wsp(j1v-n), wsp(j1v) )
         if ( j.gt.1 )
     .     call ZAXPY(n,-wsp(ih+(j-1)*mh+j-2),wsp(j1v-2*n),1,wsp(j1v),1)
         hjj = ZDOTC( n, wsp(j1v-n),1, wsp(j1v),1 )
         call ZAXPY( n, -hjj, wsp(j1v-n),1, wsp(j1v),1 )
         hj1j = DZNRM2( n, wsp(j1v),1 )
         wsp(ih+(j-1)*(mh+1)) = hjj
*---     if `happy breakdown' go straightforward at the end ...
         if ( hj1j.le.break_tol ) then
            print*,'happy breakdown: mbrkdwn =',j,' h =',hj1j
            k1 = 0
            ibrkflag = 1
            mbrkdwn = j
            tbrkdwn = t_now
            t_step = t_out-t_now
            goto 300
         endif
         wsp(ih+(j-1)*mh+j) = CMPLX( hj1j )
         wsp(ih+j*mh+j-1) = CMPLX( hj1j )
         call ZDSCAL( n, 1.0d0/hj1j, wsp(j1v),1 )
         j1v = j1v + n
 200  continue
      nmult = nmult + 1
      call matvec( wsp(j1v-n), wsp(j1v) )
      avnorm = DZNRM2( n, wsp(j1v),1 )
*
*---  set 1 for the 2-corrected scheme ...
*
 300  continue
      wsp(ih+m*mh+m-1) = ZERO
      wsp(ih+m*mh+m+1) = ONE
*
*---  loop while ireject<mxreject until the tolerance is reached ...
*
      ireject = 0
 401  continue
*
*---  compute w = beta*V*exp(t_step*H)*e1 ...
*
      nexph = nexph + 1
      mx = mbrkdwn + k1
      if ( ideg.ne.0 ) then
*---     irreducible rational Pade approximation ...
         call ZGPADM( ideg, mx, sgn*t_step, wsp(ih),mh,
     .                wsp(ifree),lfree, iwsp, iexph, ns, iflag )
         iexph = ifree + iexph - 1
         nscale = nscale + ns
      else
*---     uniform rational Chebyshev approximation ...
         iexph = ifree
         do i = 1,mx
            wsp(iexph+i-1) = ZERO
         end do
         wsp(iexph) = ONE
         call ZNCHBV(mx,sgn*t_step,wsp(ih),mh,wsp(iexph),wsp(ifree+mx))
      endif

 402  continue
*
*---  error estimate ...
*
      if ( k1.eq.0 ) then
         err_loc = tol
      else
         p1 = ABS( wsp(iexph+m) )   * beta
         p2 = ABS( wsp(iexph+m+1) ) * beta * avnorm
         if ( p1.gt.10.0d0*p2 ) then
            err_loc = p2
            xm = 1.0d0/DBLE( m )
         elseif ( p1.gt.p2 ) then
            err_loc = (p1*p2)/(p1-p2)
            xm = 1.0d0/DBLE( m )
         else
            err_loc = p1
            xm = 1.0d0/DBLE( m-1 )
         endif
      endif
*
*---  reject the step-size if the error is not acceptable ...
*
      if ( (k1.ne.0) .and. (err_loc.gt.delta*t_step*tol) .and.
     .     (mxreject.eq.0 .or. ireject.lt.mxreject) ) then
         t_old = t_step
         t_step = gamma * t_step * (t_step*tol/err_loc)**xm
         p1 = 10.0d0**(NINT( LOG10( t_step )-SQR1 )-1)
         t_step = AINT( t_step/p1 + 0.55d0 ) * p1
         if ( itrace.ne.0 ) then
            print*,'t_step =',t_old
            print*,'err_loc =',err_loc
            print*,'err_required =',delta*t_old*tol
            print*,'stepsize rejected, stepping down to:',t_step
         endif
         ireject = ireject + 1
         nreject = nreject + 1
         if ( mxreject.ne.0 .and. ireject.gt.mxreject ) then
            print*,"Failure in ZHEXPV: ---"
            print*,"The requested tolerance is too high."
            Print*,"Rerun with a smaller value."
            iflag = 2
            return
         endif
         goto 401
      endif
*
*---  now update w = beta*V*exp(t_step*H)*e1 and the hump ...
*
      mx = mbrkdwn + MAX( 0,k1-1 )
      hjj = CMPLX( beta )
      call ZGEMV( 'n', n,mx,hjj,wsp(iv),n,wsp(iexph),1,ZERO,w,1 )
      beta = DZNRM2( n, w,1 )
      hump = MAX( hump, beta )
*
*---  suggested value for the next stepsize ...
*
      t_new = gamma * t_step * (t_step*tol/err_loc)**xm
      p1 = 10.0d0**(NINT( LOG10( t_new )-SQR1 )-1)
      t_new = AINT( t_new/p1 + 0.55d0 ) * p1

      err_loc = MAX( err_loc,rndoff )
*
*---  update the time covered ...
*
      t_now = t_now + t_step
*
*---  display and keep some information ...
*
      if ( itrace.ne.0 ) then
         print*,'integration',nstep,'---------------------------------'
         print*,'scale-square =',ns
         print*,'step_size =',t_step
         print*,'err_loc   =',err_loc
         print*,'next_step =',t_new
      endif

      step_min = MIN( step_min, t_step )
      step_max = MAX( step_max, t_step )
      s_error = s_error + err_loc
      x_error = MAX( x_error, err_loc )

      if ( mxstep.eq.0 .or. nstep.lt.mxstep ) goto 100
      iflag = 1

 500  continue

      iwsp(1) = nmult
      iwsp(2) = nexph
      iwsp(3) = nscale
      iwsp(4) = nstep
      iwsp(5) = nreject
      iwsp(6) = ibrkflag
      iwsp(7) = mbrkdwn

      wsp(1)  = CMPLX( step_min )
      wsp(2)  = CMPLX( step_max )
      wsp(3)  = CMPLX( 0.0d0 )
      wsp(4)  = CMPLX( 0.0d0 )
      wsp(5)  = CMPLX( x_error )
      wsp(6)  = CMPLX( s_error )
      wsp(7)  = CMPLX( tbrkdwn )
      wsp(8)  = CMPLX( sgn*t_now )
      wsp(9)  = CMPLX( hump/vnorm )
      wsp(10) = CMPLX( beta/vnorm )
      END
      subroutine zhifa ( a, lda, n, kpvt, info )

c*********************************************************************72
c
cc ZHIFA factors a complex Hermitian linear system.
c
      implicit none

      integer lda,n,kpvt(1),info
      complex*16 a(lda,1)
c
c     zhifa factors a complex*16 hermitian matrix by elimination
c     with symmetric pivoting.
c
c     to solve  a*x = b , follow zhifa by zhisl.
c     to compute  inverse(a)*c , follow zhifa by zhisl.
c     to compute  determinant(a) , follow zhifa by zhidi.
c     to compute  inertia(a) , follow zhifa by zhidi.
c     to compute  inverse(a) , follow zhifa by zhidi.
c
c     on entry
c
c        a       complex*16(lda,n)
c                the hermitian matrix to be factored.
c                only the diagonal and upper triangle are used.
c
c        lda     integer
c                the leading dimension of the array  a .
c
c        n       integer
c                the order of the matrix  a .
c
c     on return
c
c        a       a block diagonal matrix and the multipliers which
c                were used to obtain it.
c                the factorization can be written  a = u*d*ctrans(u)
c                where  u  is a product of permutation and unit
c                upper triangular matrices , ctrans(u) is the
c                conjugate transpose of  u , and  d  is block diagonal
c                with 1 by 1 and 2 by 2 blocks.
c
c        kpvt    integer(n)
c                an integer vector of pivot indices.
c
c        info    integer
c                = 0  normal value.
c                = k  if the k-th pivot block is singular. this is
c                     not an error condition for this subroutine,
c                     but it does indicate that zhisl or zhidi may
c                     divide by zero if called.
c
c     linpack. this version dated 08/14/78 .
c     james bunch, univ. calif. san diego, argonne nat. lab.
c
c     subroutines and functions
c
c     blas zaxpy,zswap,izamax
c     fortran dabs,dmax1,dcmplx,dconjg,dsqrt
c
c     internal variables
c
      complex*16 ak,akm1,bk,bkm1,denom,mulk,mulkm1,t
      double precision absakk,alpha,colmax,rowmax
      integer imax,imaxp1,j,jj,jmax,k,km1,km2,kstep,izamax
      logical swap
c
      complex*16 zdum
      double precision cabs1
      double precision dreal,dimag
      complex*16 zdumr,zdumi
      dreal(zdumr) = zdumr
      dimag(zdumi) = (0.0d0,-1.0d0)*zdumi
      cabs1(zdum) = dabs(dreal(zdum)) + dabs(dimag(zdum))
c
c     initialize
c
c     alpha is used in choosing pivot block size.
      alpha = (1.0d0 + dsqrt(17.0d0))/8.0d0
c
      info = 0
c
c     main loop on k, which goes from n to 1.
c
      k = n
   10 continue
c
c        leave the loop if k=0 or k=1.
c
c     ...exit
         if (k .eq. 0) go to 200
         if (k .gt. 1) go to 20
            kpvt(1) = 1
            if (cabs1(a(1,1)) .eq. 0.0d0) info = 1
c     ......exit
            go to 200
   20    continue
c
c        this section of code determines the kind of
c        elimination to be performed.  when it is completed,
c        kstep will be set to the size of the pivot block, and
c        swap will be set to .true. if an interchange is
c        required.
c
         km1 = k - 1
         absakk = cabs1(a(k,k))
c
c        determine the largest off-diagonal element in
c        column k.
c
         imax = izamax(k-1,a(1,k),1)
         colmax = cabs1(a(imax,k))
         if (absakk .lt. alpha*colmax) go to 30
            kstep = 1
            swap = .false.
         go to 90
   30    continue
c
c           determine the largest off-diagonal element in
c           row imax.
c
            rowmax = 0.0d0
            imaxp1 = imax + 1
            do 40 j = imaxp1, k
               rowmax = dmax1(rowmax,cabs1(a(imax,j)))
   40       continue
            if (imax .eq. 1) go to 50
               jmax = izamax(imax-1,a(1,imax),1)
               rowmax = dmax1(rowmax,cabs1(a(jmax,imax)))
   50       continue
            if (cabs1(a(imax,imax)) .lt. alpha*rowmax) go to 60
               kstep = 1
               swap = .true.
            go to 80
   60       continue
            if (absakk .lt. alpha*colmax*(colmax/rowmax)) go to 70
               kstep = 1
               swap = .false.
            go to 80
   70       continue
               kstep = 2
               swap = imax .ne. km1
   80       continue
   90    continue
         if (dmax1(absakk,colmax) .ne. 0.0d0) go to 100
c
c           column k is zero.  set info and iterate the loop.
c
            kpvt(k) = k
            info = k
         go to 190
  100    continue
         if (kstep .eq. 2) go to 140
c
c           1 x 1 pivot block.
c
            if (.not.swap) go to 120
c
c              perform an interchange.
c
               call zswap(imax,a(1,imax),1,a(1,k),1)
               do 110 jj = imax, k
                  j = k + imax - jj
                  t = dconjg(a(j,k))
                  a(j,k) = dconjg(a(imax,j))
                  a(imax,j) = t
  110          continue
  120       continue
c
c           perform the elimination.
c
            do 130 jj = 1, km1
               j = k - jj
               mulk = -a(j,k)/a(k,k)
               t = dconjg(mulk)
               call zaxpy(j,t,a(1,k),1,a(1,j),1)
               a(j,j) = dcmplx(dreal(a(j,j)),0.0d0)
               a(j,k) = mulk
  130       continue
c
c           set the pivot array.
c
            kpvt(k) = k
            if (swap) kpvt(k) = imax
         go to 190
  140    continue
c
c           2 x 2 pivot block.
c
            if (.not.swap) go to 160
c
c              perform an interchange.
c
               call zswap(imax,a(1,imax),1,a(1,k-1),1)
               do 150 jj = imax, km1
                  j = km1 + imax - jj
                  t = dconjg(a(j,k-1))
                  a(j,k-1) = dconjg(a(imax,j))
                  a(imax,j) = t
  150          continue
               t = a(k-1,k)
               a(k-1,k) = a(imax,k)
               a(imax,k) = t
  160       continue
c
c           perform the elimination.
c
            km2 = k - 2
            if (km2 .eq. 0) go to 180
               ak = a(k,k)/a(k-1,k)
               akm1 = a(k-1,k-1)/dconjg(a(k-1,k))
               denom = 1.0d0 - ak*akm1
               do 170 jj = 1, km2
                  j = km1 - jj
                  bk = a(j,k)/a(k-1,k)
                  bkm1 = a(j,k-1)/dconjg(a(k-1,k))
                  mulk = (akm1*bk - bkm1)/denom
                  mulkm1 = (ak*bkm1 - bk)/denom
                  t = dconjg(mulk)
                  call zaxpy(j,t,a(1,k),1,a(1,j),1)
                  t = dconjg(mulkm1)
                  call zaxpy(j,t,a(1,k-1),1,a(1,j),1)
                  a(j,k) = mulk
                  a(j,k-1) = mulkm1
                  a(j,j) = dcmplx(dreal(a(j,j)),0.0d0)
  170          continue
  180       continue
c
c           set the pivot array.
c
            kpvt(k) = 1 - k
            if (swap) kpvt(k) = -imax
            kpvt(k-1) = kpvt(k)
  190    continue
         k = k - kstep
      go to 10
  200 continue
      return
      end
      subroutine zhisl ( a, lda, n, kpvt, b )

c*********************************************************************72
c
cc ZHISL solves a factored complex Hermitian linear system.
c
      implicit none

      integer lda,n,kpvt(1)
      complex*16 a(lda,1),b(1)
c
c     zhisl solves the complex*16 hermitian system
c     a * x = b
c     using the factors computed by zhifa.
c
c     on entry
c
c        a       complex*16(lda,n)
c                the output from zhifa.
c
c        lda     integer
c                the leading dimension of the array  a .
c
c        n       integer
c                the order of the matrix  a .
c
c        kpvt    integer(n)
c                the pivot vector from zhifa.
c
c        b       complex*16(n)
c                the right hand side vector.
c
c     on return
c
c        b       the solution vector  x .
c
c     error condition
c
c        a division by zero may occur if  zhico  has set rcond .eq. 0.0
c        or  zhifa  has set info .ne. 0  .
c
c     to compute  inverse(a) * c  where  c  is a matrix
c     with  p  columns
c           call zhifa(a,lda,n,kpvt,info)
c           if (info .ne. 0) go to ...
c           do 10 j = 1, p
c              call zhisl(a,lda,n,kpvt,c(1,j))
c        10 continue
c
c     linpack. this version dated 08/14/78 .
c     james bunch, univ. calif. san diego, argonne nat. lab.
c
c     subroutines and functions
c
c     blas zaxpy,zdotc
c     fortran dconjg,iabs
c
c     internal variables.
c
      complex*16 ak,akm1,bk,bkm1,zdotc,denom,temp
      integer k,kp
c
c     loop backward applying the transformations and
c     d inverse to b.
c
      k = n
   10 if (k .eq. 0) go to 80
         if (kpvt(k) .lt. 0) go to 40
c
c           1 x 1 pivot block.
c
            if (k .eq. 1) go to 30
               kp = kpvt(k)
               if (kp .eq. k) go to 20
c
c                 interchange.
c
                  temp = b(k)
                  b(k) = b(kp)
                  b(kp) = temp
   20          continue
c
c              apply the transformation.
c
               call zaxpy(k-1,b(k),a(1,k),1,b(1),1)
   30       continue
c
c           apply d inverse.
c
            b(k) = b(k)/a(k,k)
            k = k - 1
         go to 70
   40    continue
c
c           2 x 2 pivot block.
c
            if (k .eq. 2) go to 60
               kp = iabs(kpvt(k))
               if (kp .eq. k - 1) go to 50
c
c                 interchange.
c
                  temp = b(k-1)
                  b(k-1) = b(kp)
                  b(kp) = temp
   50          continue
c
c              apply the transformation.
c
               call zaxpy(k-2,b(k),a(1,k),1,b(1),1)
               call zaxpy(k-2,b(k-1),a(1,k-1),1,b(1),1)
   60       continue
c
c           apply d inverse.
c
            ak = a(k,k)/dconjg(a(k-1,k))
            akm1 = a(k-1,k-1)/a(k-1,k)
            bk = b(k)/dconjg(a(k-1,k))
            bkm1 = b(k-1)/a(k-1,k)
            denom = ak*akm1 - 1.0d0
            b(k) = (akm1*bk - bkm1)/denom
            b(k-1) = (ak*bkm1 - bk)/denom
            k = k - 2
   70    continue
      go to 10
   80 continue
c
c     loop forward applying the transformations.
c
      k = 1
   90 if (k .gt. n) go to 160
         if (kpvt(k) .lt. 0) go to 120
c
c           1 x 1 pivot block.
c
            if (k .eq. 1) go to 110
c
c              apply the transformation.
c
               b(k) = b(k) + zdotc(k-1,a(1,k),1,b(1),1)
               kp = kpvt(k)
               if (kp .eq. k) go to 100
c
c                 interchange.
c
                  temp = b(k)
                  b(k) = b(kp)
                  b(kp) = temp
  100          continue
  110       continue
            k = k + 1
         go to 150
  120    continue
c
c           2 x 2 pivot block.
c
            if (k .eq. 1) go to 140
c
c              apply the transformation.
c
               b(k) = b(k) + zdotc(k-1,a(1,k),1,b(1),1)
               b(k+1) = b(k+1) + zdotc(k-1,a(1,k+1),1,b(1),1)
               kp = iabs(kpvt(k))
               if (kp .eq. k) go to 130
c
c                 interchange.
c
                  temp = b(k)
                  b(k) = b(kp)
                  b(kp) = temp
  130          continue
  140       continue
            k = k + 2
  150    continue
      go to 90
  160 continue
      return
      end
      subroutine zhpadm(ideg,m,t,h,ldh,wsp,lwsp,ipiv,iexph,ns,iflag)

c*********************************************************************72
c
cc ZHPADM computes exp(t*H), the matrix exponential of a Hermitian matrix.
c
      implicit none

      double precision t
      integer          ideg, m, ldh, lwsp, iexph, ns, iflag, ipiv(m)
      complex*16       H(ldh,m), wsp(lwsp)

*-----Purpose----------------------------------------------------------|
*
*     Computes exp(t*H), the matrix exponential of a Hermitian matrix
*     in full, using the irreducible rational Pade approximation to the
*     exponential function exp(z) = r(z) = (+/-)( I + 2*(q(z)/p(z)) ),
*     combined with scaling-and-squaring.
*
*-----Arguments--------------------------------------------------------|
*
*     ideg      : (input) the degre of the diagonal Pade to be used.
*                 a value of 6 is generally satisfactory.
*
*     m         : (input) order of H.
*
*     H(ldh,m)  : (input) argument matrix (both lower and upper parts).
*
*     t         : (input) time-scale (can be < 0).
*
*     wsp(lwsp) : (workspace/output) lwsp .ge. 4*m*m+ideg+1.
*
*     ipiv(m)   : (workspace)
*
*>>>> iexph     : (output) number such that wsp(iexph) points to exp(tH)
*                 i.e., exp(tH) is located at wsp(iexph ... iexph+m*m-1)
*                       ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
*                 NOTE: if the routine was called with wsp(iptr),
*                       then exp(tH) will start at wsp(iptr+iexph-1).
*
*     ns        : (output) number of scaling-squaring used.
*
*     iflag     : (output) exit flag.
*                       0 - no problem
*                      <0 - problem
*
*----------------------------------------------------------------------|
*     Roger B. Sidje (rbs@maths.uq.edu.au)
*     EXPOKIT: Software Package for Computing Matrix Exponentials.
*     ACM - Transactions On Mathematical Software, 24(1):130-156, 1998
*----------------------------------------------------------------------|
*
      integer i,j,k,icoef,mm,ih2,iodd,iused,ifree,iq,ip,iput,iget
      double precision hnorm
      complex*16 cp, cq, scale, scale2, ZERO, ONE

      parameter( ZERO=(0.0d0,0.0d0), ONE=(1.0d0,0.0d0) )
      intrinsic ABS, CMPLX, DBLE, INT, LOG, MAX

*---  check restrictions on input parameters ...
      mm = m*m
      iflag = 0
      if ( ldh.lt.m ) iflag = -1
      if ( lwsp.lt.4*mm+ideg+1 ) iflag = -2
      if ( iflag.ne.0 ) stop 'bad sizes (in input of ZHPADM)'
*
*---  initialise pointers ...
*
      icoef = 1
      ih2 = icoef + (ideg+1)
      ip  = ih2 + mm
      iq  = ip + mm
      ifree = iq + mm
*
*---  scaling: seek ns such that ||t*H/2^ns|| < 1/2;
*     and set scale = t/2^ns ...
*
      do i = 1,m
         wsp(i) = ZERO
      end do
      do j = 1,m
         do i = 1,m
            wsp(i) = wsp(i) + ABS( H(i,j) )
         end do
      end do
      hnorm = 0.0d0
      do i = 1,m
         hnorm = MAX( hnorm,DBLE(wsp(i)) )
      end do
      hnorm = ABS( t*hnorm )
      if ( hnorm.eq.0.0d0 ) stop 'Error - null H in input of ZHPADM.'
      ns = MAX( 0,INT(LOG(hnorm)/LOG(2.0d0))+2 )
      scale =  CMPLX( t/DBLE(2**ns),0.0d0 )
      scale2 = scale*scale
*
*---  compute Pade coefficients ...
*
      i = ideg+1
      j = 2*ideg+1
      wsp(icoef) = ONE
      do k = 1,ideg
         wsp(icoef+k) = (wsp(icoef+k-1)*DBLE( i-k ))/DBLE( k*(j-k) )
      end do
*
*---  H2 = scale2*H*H ...
*
      call ZGEMM( 'n','n',m,m,m,scale2,H,ldh,H,ldh,ZERO,wsp(ih2),m )
*
*---  initialise p (numerator) and q (denominator) ...
*
      cp = wsp(icoef+ideg-1)
      cq = wsp(icoef+ideg)
      do j = 1,m
         do i = 1,m
            wsp(ip + (j-1)*m + i-1) = ZERO
            wsp(iq + (j-1)*m + i-1) = ZERO
         end do
         wsp(ip + (j-1)*(m+1)) = cp
         wsp(iq + (j-1)*(m+1)) = cq
      end do
*
*---  Apply Horner rule ...
*
      iodd = 1
      k = ideg - 1
 100  continue
      iused = iodd*iq + (1-iodd)*ip
      call ZGEMM( 'n','n',m,m,m, ONE,wsp(iused),m,
     .             wsp(ih2),m, ZERO,wsp(ifree),m )
      do j = 1,m
         wsp(ifree+(j-1)*(m+1)) = wsp(ifree+(j-1)*(m+1))+wsp(icoef+k-1)
      end do
      ip = (1-iodd)*ifree + iodd*ip
      iq = iodd*ifree + (1-iodd)*iq
      ifree = iused
      iodd = 1-iodd
      k = k-1
      if ( k.gt.0 )  goto 100
*
*---  Obtain (+/-)(I + 2*(p\q)) ...
*
      if ( iodd.ne.0 ) then
         call ZGEMM( 'n','n',m,m,m, scale,wsp(iq),m,
     .                H,ldh, ZERO,wsp(ifree),m )
         iq = ifree
      else
         call ZGEMM( 'n','n',m,m,m, scale,wsp(ip),m,
     .                H,ldh, ZERO,wsp(ifree),m )
         ip = ifree
      endif
      call ZAXPY( mm, -ONE,wsp(ip),1, wsp(iq),1 )
      call ZHESV( 'U',m,m,wsp(iq),m,ipiv,wsp(ip),m,wsp(ih2),mm,iflag )
      if ( iflag.ne.0 ) stop 'Problem in ZHESV (within ZHPADM)'
      call ZDSCAL( mm, 2.0d0, wsp(ip), 1 )
      do j = 1,m
         wsp(ip+(j-1)*(m+1)) = wsp(ip+(j-1)*(m+1)) + ONE
      end do
      iput = ip
      if ( ns.eq.0 .and. iodd.ne.0 ) then
         call ZDSCAL( mm, -1.0d0, wsp(ip), 1 )
         goto 200
      endif
*
*--   squaring : exp(t*H) = (exp(t*H))^(2^ns) ...
*
      iodd = 1
      do k = 1,ns
         iget = iodd*ip + (1-iodd)*iq
         iput = (1-iodd)*ip + iodd*iq
         call ZGEMM( 'n','n',m,m,m, ONE,wsp(iget),m, wsp(iget),m,
     .                ZERO,wsp(iput),m )
         iodd = 1-iodd
      end do
 200  continue
      iexph = iput
      end
      subroutine zhphiv ( n, m, t, u, v, w, tol, anorm,
     .                   wsp,lwsp, iwsp,liwsp, matvec, itrace,iflag )

c*********************************************************************72
c
cc ZHPHIV solves the linear ODE w' = Aw + u, w(0) = v for A hermitian.
c
      implicit none

      integer          n, m, lwsp, liwsp, itrace, iflag, iwsp(liwsp)
      double precision t, tol, anorm
      complex*16       u(n), v(n), w(n), wsp(lwsp)
      external         matvec

*-----Purpose----------------------------------------------------------|
*
*---  ZHPHIV computes w = exp(t*A)v + t*phi(tA)u which is the solution
*     of the nonhomogeneous linear ODE problem w' = Aw + u, w(0) = v.
*     phi(z) = (exp(z)-1)/z and A is an Hermitian matrix.
*
*     The method used is based on Krylov subspace projection
*     techniques and the matrix under consideration interacts only
*     via the external routine `matvec' performing the matrix-vector
*     product (matrix-free method).
*
*-----Arguments--------------------------------------------------------|
*
*     n      : (input) order of the principal matrix A.
*
*     m      : (input) maximum size for the Krylov basis.
*
*     t      : (input) time at wich the solution is needed (can be < 0).
*
*     u(n)   : (input) operand vector with respect to the phi function
*              (forcing term of the ODE problem).
*
*     v(n)   : (input) operand vector with respect to the exp function
*              (initial condition of the ODE problem).
*
*     w(n)   : (output) computed approximation of exp(t*A)v + t*phi(tA)u
*
*     tol    : (input/output) the requested accuracy tolerance on w.
*              If on input tol=0.0d0 or tol is too small (tol.le.eps)
*              the internal value sqrt(eps) is used, and tol is set to
*              sqrt(eps) on output (`eps' denotes the machine epsilon).
*              (`Happy breakdown' is assumed if h(j+1,j) .le. anorm*tol)
*
*     anorm  : (input) an approximation of some norm of A.
*
*   wsp(lwsp): (workspace) lwsp .ge. n*(m+1)+n+(m+3)^2+4*(m+3)^2+ideg+1
*                                   +---------+-------+---------------+
*              (actually, ideg=6)        V        H      wsp for PADE
*
* iwsp(liwsp): (workspace) liwsp .ge. m+3
*
*     matvec : external routine for matrix-vector multiplication.
*              synopsis: matvec( x, y )
*                        double precision x(*), y(*)
*              computes: y(1:n) <- A*x(1:n)
*                        where A is the principal matrix.
*
*     itrace : (input) running mode. 0=silent, 1=print step-by-step info
*
*     iflag  : (output) exit flag.
*              <0 - bad input arguments
*               0 - no problem
*               1 - maximum number of steps reached without convergence
*               2 - requested tolerance was too high
*
*-----Accounts on the computation--------------------------------------|
*     Upon exit, an interested user may retrieve accounts on the
*     computations. They are located in the workspace arrays wsp and
*     iwsp as indicated below:
*
*     location  mnemonic                 description
*     -----------------------------------------------------------------|
*     iwsp(1) = nmult, number of matrix-vector multiplications used
*     iwsp(2) = nexph, number of Hessenberg matrix exponential evaluated
*     iwsp(3) = nscale, number of repeated squaring involved in Pade
*     iwsp(4) = nstep, number of integration steps used up to completion
*     iwsp(5) = nreject, number of rejected step-sizes
*     iwsp(6) = ibrkflag, set to 1 if `happy breakdown' and 0 otherwise
*     iwsp(7) = mbrkdwn, if `happy brkdown', basis-size when it occured
*     -----------------------------------------------------------------|
*     wsp(1)  = step_min, minimum step-size used during integration
*     wsp(2)  = step_max, maximum step-size used during integration
*     wsp(3)  = dummy
*     wsp(4)  = dummy
*     wsp(5)  = x_error, maximum among all local truncation errors
*     wsp(6)  = s_error, global sum of local truncation errors
*     wsp(7)  = tbrkdwn, if `happy breakdown', time when it occured
*     wsp(8)  = t_now, integration domain successfully covered
*
*----------------------------------------------------------------------|
*-----The following parameters may also be adjusted herein-------------|
*
      integer mxstep, mxreject, ideg
      double precision delta, gamma
      parameter( mxstep   = 500,
     .           mxreject = 0,
     .           ideg     = 6,
     .           delta    = 1.2d0,
     .           gamma    = 0.9d0 )

*     mxstep  : maximum allowable number of integration steps.
*               The value 0 means an infinite number of steps.
*
*     mxreject: maximum allowable number of rejections at each step.
*               The value 0 means an infinite number of rejections.
*
*     ideg    : the Pade approximation of type (ideg,ideg) is used as
*               an approximation to exp(H).
*
*     delta   : local truncation error `safety factor'
*
*     gamma   : stepsize `shrinking factor'
*
*----------------------------------------------------------------------|
*     Roger B. Sidje (rbs@maths.uq.edu.au)
*     EXPOKIT: Software Package for Computing Matrix Exponentials.
*     ACM - Transactions On Mathematical Software, 24(1):130-156, 1998
*----------------------------------------------------------------------|
*
      complex*16 ZERO, ONE
      parameter( ZERO=(0.0d0,0.0d0), ONE=(1.0d0,0.0d0) )
      integer i, j, k1, mh, mx, iv, ih, j1v, ns, ifree, lfree, iexph,
     .        ireject,ibrkflag,mbrkdwn, nmult, nreject, nexph, nscale,
     .        nstep, iphih
      double precision sgn, t_out, tbrkdwn, step_min,step_max, err_loc,
     .                 s_error, x_error, t_now, t_new, t_step, t_old,
     .                 xm, beta, break_tol, p1, p2, p3, eps, rndoff,
     .                 avnorm, hj1j, SQR1
      complex*16 hjj

      intrinsic AINT,ABS,CMPLX,DBLE,INT,LOG10,MAX,MIN,NINT,SIGN,SQRT
      complex*16 ZDOTC
      double precision DZNRM2
*
*---  check restrictions on input parameters ...
*
      iflag = 0
      if ( lwsp.lt.n*(m+3)+5*(m+3)**2+ideg+1 ) iflag = -1
      if ( liwsp.lt.m+3 ) iflag = -2
      if ( m.ge.n .or. m.le.0 ) iflag = -3
      if ( iflag.ne.0 ) stop 'bad sizes (in input of ZHPHIV)'
*
*---  initialisations ...
*
      k1 = 3
      mh = m + 3
      iv = 1
      ih = iv + n*(m+1) + n
      ifree = ih + mh*mh
      lfree = lwsp - ifree + 1

      ibrkflag = 0
      mbrkdwn  = m
      nmult    = 0
      nreject  = 0
      nexph    = 0
      nscale   = 0

      t_out    = ABS( t )
      tbrkdwn  = 0.0d0
      step_min = t_out
      step_max = 0.0d0
      nstep    = 0
      s_error  = 0.0d0
      x_error  = 0.0d0
      t_now    = 0.0d0
      t_new    = 0.0d0

      p1 = 4.0d0/3.0d0
 1    p2 = p1 - 1.0d0
      p3 = p2 + p2 + p2
      eps = ABS( p3-1.0d0 )
      if ( eps.eq.0.0d0 ) go to 1
      if ( tol.le.eps ) tol = SQRT( eps )
      rndoff = eps*anorm

      break_tol = 1.0d-7
*>>>  break_tol = tol
*>>>  break_tol = anorm*tol

      sgn = SIGN( 1.0d0,t )
      SQR1 = SQRT( 0.1d0 )
      call ZCOPY( n, v,1, w,1 )
*
*---  step-by-step integration ...
*
 100  if ( t_now.ge.t_out ) goto 500

      nmult =  nmult + 1
      call matvec( w, wsp(iv) )
      call ZAXPY( n, ONE, u,1, wsp(iv),1 )
      beta = DZNRM2( n, wsp(iv),1 )
      if ( beta.eq.0.0d0 ) goto 500
      call ZDSCAL( n, 1.0d0/beta, wsp(iv),1 )
      do i = 1,mh*mh
         wsp(ih+i-1) = ZERO
      end do

      if ( nstep.eq.0 ) then
*---     obtain the very first stepsize ...
         xm = 1.0d0/DBLE( m )
         p1 = tol*(((m+1)/2.72D0)**(m+1))*SQRT(2.0D0*3.14D0*(m+1))
         t_new = (1.0d0/anorm)*(p1/(4.0d0*beta*anorm))**xm
         p1 = 10.0d0**(NINT( LOG10( t_new )-SQR1 )-1)
         t_new = AINT( t_new/p1 + 0.55d0 ) * p1
      endif
      nstep = nstep + 1
      t_step = MIN( t_out-t_now, t_new )
*
*---  Lanczos loop ...
*
      j1v = iv + n
      do 200 j = 1,m
         nmult = nmult + 1
         call matvec( wsp(j1v-n), wsp(j1v) )
         if ( j.gt.1 )
     .     call ZAXPY(n,-wsp(ih+(j-1)*mh+j-2),wsp(j1v-2*n),1,wsp(j1v),1)
         hjj = ZDOTC( n, wsp(j1v-n),1, wsp(j1v),1 )
         call ZAXPY( n, -hjj, wsp(j1v-n),1, wsp(j1v),1 )
         hj1j = DZNRM2( n, wsp(j1v),1 )
         wsp(ih+(j-1)*(mh+1)) = hjj
         if ( hj1j.le.break_tol ) then
            print*,'happy breakdown: mbrkdwn =',j,' h =',hj1j
            k1 = 0
            ibrkflag = 1
            mbrkdwn = j
            tbrkdwn = t_now
            t_step = t_out-t_now
            goto 300
         endif
         wsp(ih+(j-1)*mh+j) = CMPLX( hj1j )
         wsp(ih+j*mh+j-1) = CMPLX( hj1j )
         call ZDSCAL( n, 1.0d0/hj1j, wsp(j1v),1 )
         j1v = j1v + n
 200  continue
      nmult = nmult + 1
      call matvec( wsp(j1v-n), wsp(j1v) )
      avnorm = DZNRM2( n, wsp(j1v),1 )
*
*---  set 1's for the 3-extended scheme ...
*
 300  continue
      wsp(ih+mh*mbrkdwn) = ONE
      wsp(ih+m*mh+m-1)   = ZERO
      wsp(ih+(m-1)*mh+m) = ZERO
      do i = 1,k1-1
         wsp(ih+(m+i)*mh+m+i-1) = ONE
      end do
*
*---  loop while ireject<mxreject until the tolerance is reached ...
*
      ireject = 0
 401  continue
*
*---  compute w = beta*t_step*V*phi(t_step*H)*e1 + w
*
      nexph = nexph + 1
      mx = mbrkdwn + MAX( 1,k1 )

*---  irreducible rational Pade approximation ...
      call ZGPADM( ideg, mx, sgn*t_step, wsp(ih),mh,
     .              wsp(ifree),lfree, iwsp, iexph, ns, iflag )
      iexph = ifree + iexph - 1
      iphih = iexph + mbrkdwn*mx
      nscale = nscale + ns
      wsp(iphih+mbrkdwn)   = hj1j*wsp(iphih+mx+mbrkdwn-1)
      wsp(iphih+mbrkdwn+1) = hj1j*wsp(iphih+2*mx+mbrkdwn-1)

 402  continue
*
*---  error estimate ...
*
      if ( k1.eq.0 ) then
         err_loc = tol
      else
         p1 = ABS( wsp(iphih+m) )   * beta
         p2 = ABS( wsp(iphih+m+1) ) * beta * avnorm
         if ( p1.gt.10.0d0*p2 ) then
            err_loc = p2
            xm = 1.0d0/DBLE( m+1 )
         elseif ( p1.gt.p2 ) then
            err_loc = (p1*p2)/(p1-p2)
            xm = 1.0d0/DBLE( m+1 )
         else
            err_loc = p1
            xm = 1.0d0/DBLE( m )
         endif
      endif
*
*---  reject the step-size if the error is not acceptable ...
*
      if ( (k1.ne.0) .and. (err_loc.gt.delta*t_step*tol) .and.
     .     (mxreject.eq.0 .or. ireject.lt.mxreject) ) then
         t_old = t_step
         t_step = gamma * t_step * (t_step*tol/err_loc)**xm
         p1 = 10.0d0**(NINT( LOG10( t_step )-SQR1 )-1)
         t_step = AINT( t_step/p1 + 0.55d0 ) * p1
         if ( itrace.ne.0 ) then
            print*,'t_step =',t_old
            print*,'err_loc =',err_loc
            print*,'err_required =',delta*t_old*tol
            print*,'stepsize rejected, stepping down to:',t_step
         endif
         ireject = ireject + 1
         nreject = nreject + 1
         if ( mxreject.ne.0 .and. ireject.gt.mxreject ) then
            print*,"Failure in ZHPHIV: ---"
            print*,"The requested tolerance is too high."
            Print*,"Rerun with a smaller value."
            iflag = 2
            return
         endif
         goto 401
      endif
*
      mx = mbrkdwn + MAX( 0,k1-2 )
      hjj = CMPLX( beta )
      call ZGEMV( 'n', n,mx,hjj,wsp(iv),n,wsp(iphih),1,ONE,w,1 )
*
*---  suggested value for the next stepsize ...
*
      t_new = gamma * t_step * (t_step*tol/err_loc)**xm
      p1 = 10.0d0**(NINT( LOG10( t_new )-SQR1 )-1)
      t_new = AINT( t_new/p1 + 0.55d0 ) * p1

      err_loc = MAX( err_loc,rndoff )
*
*---  update the time covered ...
*
      t_now = t_now + t_step
*
*---  display and keep some information ...
*
      if ( itrace.ne.0 ) then
         print*,'integration',nstep,'---------------------------------'
         print*,'scale-square =',ns
         print*,'step_size =',t_step
         print*,'err_loc   =',err_loc
         print*,'next_step =',t_new
      endif

      step_min = MIN( step_min, t_step )
      step_max = MAX( step_max, t_step )
      s_error = s_error + err_loc
      x_error = MAX( x_error, err_loc )

      if ( mxstep.eq.0 .or. nstep.lt.mxstep ) goto 100
      iflag = 1

 500  continue

      iwsp(1) = nmult
      iwsp(2) = nexph
      iwsp(3) = nscale
      iwsp(4) = nstep
      iwsp(5) = nreject
      iwsp(6) = ibrkflag
      iwsp(7) = mbrkdwn

      wsp(1)  = CMPLX( step_min )
      wsp(2)  = CMPLX( step_max )
      wsp(3)  = CMPLX( 0.0d0 )
      wsp(4)  = CMPLX( 0.0d0 )
      wsp(5)  = CMPLX( x_error )
      wsp(6)  = CMPLX( s_error )
      wsp(7)  = CMPLX( tbrkdwn )
      wsp(8)  = CMPLX( sgn*t_now )
      end
      subroutine znchbv ( m, t, h,ldh, y, wsp )

c*********************************************************************72
c
cc ZNCHBV computes y = exp(t*H)*y using a partial fraction expansion.
c
      implicit none

      integer          m, ldh
      double precision t
      complex*16       H(ldh,m), y(m), wsp(m*(m+2))

*-----Purpose----------------------------------------------------------|
*
*---  ZNCHBV computes y = exp(t*H)*y using the partial fraction
*     expansion of the uniform rational Chebyshev approximation
*     to exp(-x) of type (14,14). H is assumed to be upper-Hessenberg.
*     About 14-digit accuracy is expected if the matrix H is negative
*     definite. The algorithm may behave poorly otherwise.
*
*-----Arguments--------------------------------------------------------|
*
*     m       : (input) order of the Hessenberg matrix H
*
*     t       : (input) time-scaling factor (can be < 0).
*
*     H(ldh,m): (input) upper Hessenberg matrix.
*
*     y(m)    : (input/output) on input the operand vector,
*               on output the resulting vector exp(t*H)*y.
*
*     wsp     : (workspace). Observe that a double precision vector of
*               length 2*m*(m+2) can be used as well when calling this
*               routine.
*
*----------------------------------------------------------------------|
*     Roger B. Sidje (rbs@maths.uq.edu.au)
*     EXPOKIT: Software Package for Computing Matrix Exponentials.
*     ACM - Transactions On Mathematical Software, 24(1):130-156, 1998
*----------------------------------------------------------------------|
*
      complex*16 ZERO
      integer ndeg, i, j, k, ip, ih, iy, iz
      parameter ( ndeg=7, ZERO=(0.0d0,0.0d0) )
      double precision alpha0
      complex*16 alpha(ndeg), theta(ndeg), tmpc

      intrinsic ABS,DBLE,CONJG,MIN

*---  Pointers ...

      ih = 1
      iy = ih + m*m
      iz = iy + m

*---  Coefficients and poles of the partial fraction expansion...

      alpha0  =  0.183216998528140087D-11
      alpha(1)=( 0.557503973136501826D+02,-0.204295038779771857D+03)
      alpha(2)=(-0.938666838877006739D+02, 0.912874896775456363D+02)
      alpha(3)=( 0.469965415550370835D+02,-0.116167609985818103D+02)
      alpha(4)=(-0.961424200626061065D+01,-0.264195613880262669D+01)
      alpha(5)=( 0.752722063978321642D+00, 0.670367365566377770D+00)
      alpha(6)=(-0.188781253158648576D-01,-0.343696176445802414D-01)
      alpha(7)=( 0.143086431411801849D-03, 0.287221133228814096D-03)

      theta(1)=(-0.562314417475317895D+01, 0.119406921611247440D+01)
      theta(2)=(-0.508934679728216110D+01, 0.358882439228376881D+01)
      theta(3)=(-0.399337136365302569D+01, 0.600483209099604664D+01)
      theta(4)=(-0.226978543095856366D+01, 0.846173881758693369D+01)
      theta(5)=( 0.208756929753827868D+00, 0.109912615662209418D+02)
      theta(6)=( 0.370327340957595652D+01, 0.136563731924991884D+02)
      theta(7)=( 0.889777151877331107D+01, 0.166309842834712071D+02)
*
      do ip = 1,ndeg
         theta(ndeg+ip) = CONJG( theta(ip) )
         alpha(ndeg+ip) = CONJG( alpha(ip) )
      end do
*
*---  Accumulation of the contribution of each pole ...
*
      do j = 1,m
         wsp(iz+j-1) = y(j)
         y(j) = y(j)*alpha0
      end do
      do ip = 1,2*ndeg
         alpha(ip) = 0.5d0*alpha(ip)
*---     Solve each fraction using Gaussian elimination with pivoting...
         do j = 1,m
            wsp(iy+j-1) = wsp(iz+j-1)
            do i = 1,MIN( j+1,m )
               wsp(ih+(j-1)*m+i-1) = -t*H(i,j)
            end do
            wsp(ih+(j-1)*m+j-1) = wsp(ih+(j-1)*m+j-1)-theta(ip)
            do k = i,m
               wsp(ih+(j-1)*m+k-1) = ZERO
            end do
         end do
         do i = 1,m-1
*---        Get pivot and exchange rows ...
            if (ABS(wsp(ih+(i-1)*m+i-1)).lt.ABS(wsp(ih+(i-1)*m+i))) then
               call ZSWAP( m-i+1, wsp(ih+(i-1)*m+i-1),m,
     .                     wsp(ih+(i-1)*m+i),m )
               call ZSWAP( 1, wsp(iy+i-1),1, wsp(iy+i),1 )
            endif
*---        Forward eliminiation ...
            tmpc = wsp(ih+(i-1)*m+i) / wsp(ih+(i-1)*m+i-1)
            call ZAXPY( m-i, -tmpc, wsp(ih+i*m+i-1),m, wsp(ih+i*m+i),m )
            wsp(iy+i) = wsp(iy+i) - tmpc*wsp(iy+i-1)
         end do
*---     Backward substitution ...
         do i = m,1,-1
            tmpc = wsp(iy+i-1)
            do j = i+1,m
               tmpc = tmpc - wsp(ih+(j-1)*m+i-1)*wsp(iy+j-1)
            end do
            wsp(iy+i-1) = tmpc / wsp(ih+(i-1)*m+i-1)
         end do
*---     Accumulate the partial result in y ...
         do j = 1,m
            y(j) = y(j) + alpha(ip)*wsp(iy+j-1)
         end do
      end do
      END
      subroutine zscal ( n, za, zx, incx )

c*********************************************************************72
c
cc ZSCAL scales a vector by a constant.
c
c     jack dongarra, 3/11/78.
c     modified 3/93 to return if incx .le. 0.
c
      implicit none

      double complex za,zx(1)
      integer i,incx,ix,n
c
      if( n.le.0 .or. incx.le.0 )return
      if(incx.eq.1)go to 20
c
c        code for increment not equal to 1
c
      ix = 1
      do 10 i = 1,n
        zx(ix) = za*zx(ix)
        ix = ix + incx
   10 continue
      return
c
c        code for increment equal to 1
c
   20 do 30 i = 1,n
        zx(i) = za*zx(i)
   30 continue
      return
      end
      subroutine zsifa ( a, lda, n, kpvt, info )

c*********************************************************************72
c
cc ZSIFA factors a complex symmetric matrix.
c
      implicit none

      integer lda,n,kpvt(1),info
      complex*16 a(lda,1)
c
c     zsifa factors a complex*16 symmetric matrix by elimination
c     with symmetric pivoting.
c
c     to solve  a*x = b , follow zsifa by zsisl.
c     to compute  inverse(a)*c , follow zsifa by zsisl.
c     to compute  determinant(a) , follow zsifa by zsidi.
c     to compute  inverse(a) , follow zsifa by zsidi.
c
c     on entry
c
c        a       complex*16(lda,n)
c                the symmetric matrix to be factored.
c                only the diagonal and upper triangle are used.
c
c        lda     integer
c                the leading dimension of the array  a .
c
c        n       integer
c                the order of the matrix  a .
c
c     on return
c
c        a       a block diagonal matrix and the multipliers which
c                were used to obtain it.
c                the factorization can be written  a = u*d*trans(u)
c                where  u  is a product of permutation and unit
c                upper triangular matrices , trans(u) is the
c                transpose of  u , and  d  is block diagonal
c                with 1 by 1 and 2 by 2 blocks.
c
c        kpvt    integer(n)
c                an integer vector of pivot indices.
c
c        info    integer
c                = 0  normal value.
c                = k  if the k-th pivot block is singular. this is
c                     not an error condition for this subroutine,
c                     but it does indicate that zsisl or zsidi may
c                     divide by zero if called.
c
c     linpack. this version dated 08/14/78 .
c     james bunch, univ. calif. san diego, argonne nat. lab.
c
c     subroutines and functions
c
c     blas zaxpy,zswap,izamax
c     fortran dabs,dmax1,dsqrt
c
c     internal variables
c
      complex*16 ak,akm1,bk,bkm1,denom,mulk,mulkm1,t
      double precision absakk,alpha,colmax,rowmax
      integer imax,imaxp1,j,jj,jmax,k,km1,km2,kstep,izamax
      logical swap
c
      complex*16 zdum
      double precision cabs1
      double precision dreal,dimag
      complex*16 zdumr,zdumi
      dreal(zdumr) = zdumr
      dimag(zdumi) = (0.0d0,-1.0d0)*zdumi
      cabs1(zdum) = dabs(dreal(zdum)) + dabs(dimag(zdum))
c
c     initialize
c
c     alpha is used in choosing pivot block size.
      alpha = (1.0d0 + dsqrt(17.0d0))/8.0d0
c
      info = 0
c
c     main loop on k, which goes from n to 1.
c
      k = n
   10 continue
c
c        leave the loop if k=0 or k=1.
c
c     ...exit
         if (k .eq. 0) go to 200
         if (k .gt. 1) go to 20
            kpvt(1) = 1
            if (cabs1(a(1,1)) .eq. 0.0d0) info = 1
c     ......exit
            go to 200
   20    continue
c
c        this section of code determines the kind of
c        elimination to be performed.  when it is completed,
c        kstep will be set to the size of the pivot block, and
c        swap will be set to .true. if an interchange is
c        required.
c
         km1 = k - 1
         absakk = cabs1(a(k,k))
c
c        determine the largest off-diagonal element in
c        column k.
c
         imax = izamax(k-1,a(1,k),1)
         colmax = cabs1(a(imax,k))
         if (absakk .lt. alpha*colmax) go to 30
            kstep = 1
            swap = .false.
         go to 90
   30    continue
c
c           determine the largest off-diagonal element in
c           row imax.
c
            rowmax = 0.0d0
            imaxp1 = imax + 1
            do 40 j = imaxp1, k
               rowmax = dmax1(rowmax,cabs1(a(imax,j)))
   40       continue
            if (imax .eq. 1) go to 50
               jmax = izamax(imax-1,a(1,imax),1)
               rowmax = dmax1(rowmax,cabs1(a(jmax,imax)))
   50       continue
            if (cabs1(a(imax,imax)) .lt. alpha*rowmax) go to 60
               kstep = 1
               swap = .true.
            go to 80
   60       continue
            if (absakk .lt. alpha*colmax*(colmax/rowmax)) go to 70
               kstep = 1
               swap = .false.
            go to 80
   70       continue
               kstep = 2
               swap = imax .ne. km1
   80       continue
   90    continue
         if (dmax1(absakk,colmax) .ne. 0.0d0) go to 100
c
c           column k is zero.  set info and iterate the loop.
c
            kpvt(k) = k
            info = k
         go to 190
  100    continue
         if (kstep .eq. 2) go to 140
c
c           1 x 1 pivot block.
c
            if (.not.swap) go to 120
c
c              perform an interchange.
c
               call zswap(imax,a(1,imax),1,a(1,k),1)
               do 110 jj = imax, k
                  j = k + imax - jj
                  t = a(j,k)
                  a(j,k) = a(imax,j)
                  a(imax,j) = t
  110          continue
  120       continue
c
c           perform the elimination.
c
            do 130 jj = 1, km1
               j = k - jj
               mulk = -a(j,k)/a(k,k)
               t = mulk
               call zaxpy(j,t,a(1,k),1,a(1,j),1)
               a(j,k) = mulk
  130       continue
c
c           set the pivot array.
c
            kpvt(k) = k
            if (swap) kpvt(k) = imax
         go to 190
  140    continue
c
c           2 x 2 pivot block.
c
            if (.not.swap) go to 160
c
c              perform an interchange.
c
               call zswap(imax,a(1,imax),1,a(1,k-1),1)
               do 150 jj = imax, km1
                  j = km1 + imax - jj
                  t = a(j,k-1)
                  a(j,k-1) = a(imax,j)
                  a(imax,j) = t
  150          continue
               t = a(k-1,k)
               a(k-1,k) = a(imax,k)
               a(imax,k) = t
  160       continue
c
c           perform the elimination.
c
            km2 = k - 2
            if (km2 .eq. 0) go to 180
               ak = a(k,k)/a(k-1,k)
               akm1 = a(k-1,k-1)/a(k-1,k)
               denom = 1.0d0 - ak*akm1
               do 170 jj = 1, km2
                  j = km1 - jj
                  bk = a(j,k)/a(k-1,k)
                  bkm1 = a(j,k-1)/a(k-1,k)
                  mulk = (akm1*bk - bkm1)/denom
                  mulkm1 = (ak*bkm1 - bk)/denom
                  t = mulk
                  call zaxpy(j,t,a(1,k),1,a(1,j),1)
                  t = mulkm1
                  call zaxpy(j,t,a(1,k-1),1,a(1,j),1)
                  a(j,k) = mulk
                  a(j,k-1) = mulkm1
  170          continue
  180       continue
c
c           set the pivot array.
c
            kpvt(k) = 1 - k
            if (swap) kpvt(k) = -imax
            kpvt(k-1) = kpvt(k)
  190    continue
         k = k - kstep
      go to 10
  200 continue
      return
      end
      subroutine zsisl ( a, lda, n, kpvt, b )

c*********************************************************************72
c
cc ZSISL solves a factored complex symmetric linear system.
c
      implicit none

      integer lda,n,kpvt(1)
      complex*16 a(lda,1),b(1)
c
c     zsisl solves the complex*16 symmetric system
c     a * x = b
c     using the factors computed by zsifa.
c
c     on entry
c
c        a       complex*16(lda,n)
c                the output from zsifa.
c
c        lda     integer
c                the leading dimension of the array  a .
c
c        n       integer
c                the order of the matrix  a .
c
c        kpvt    integer(n)
c                the pivot vector from zsifa.
c
c        b       complex*16(n)
c                the right hand side vector.
c
c     on return
c
c        b       the solution vector  x .
c
c     error condition
c
c        a division by zero may occur if  zsico  has set rcond .eq. 0.0
c        or  zsifa  has set info .ne. 0  .
c
c     to compute  inverse(a) * c  where  c  is a matrix
c     with  p  columns
c           call zsifa(a,lda,n,kpvt,info)
c           if (info .ne. 0) go to ...
c           do 10 j = 1, p
c              call zsisl(a,lda,n,kpvt,c(1,j))
c        10 continue
c
c     linpack. this version dated 08/14/78 .
c     james bunch, univ. calif. san diego, argonne nat. lab.
c
c     subroutines and functions
c
c     blas zaxpy,zdotu
c     fortran iabs
c
c     internal variables.
c
      complex*16 ak,akm1,bk,bkm1,zdotu,denom,temp
      integer k,kp
c
c     loop backward applying the transformations and
c     d inverse to b.
c
      k = n
   10 if (k .eq. 0) go to 80
         if (kpvt(k) .lt. 0) go to 40
c
c           1 x 1 pivot block.
c
            if (k .eq. 1) go to 30
               kp = kpvt(k)
               if (kp .eq. k) go to 20
c
c                 interchange.
c
                  temp = b(k)
                  b(k) = b(kp)
                  b(kp) = temp
   20          continue
c
c              apply the transformation.
c
               call zaxpy(k-1,b(k),a(1,k),1,b(1),1)
   30       continue
c
c           apply d inverse.
c
            b(k) = b(k)/a(k,k)
            k = k - 1
         go to 70
   40    continue
c
c           2 x 2 pivot block.
c
            if (k .eq. 2) go to 60
               kp = iabs(kpvt(k))
               if (kp .eq. k - 1) go to 50
c
c                 interchange.
c
                  temp = b(k-1)
                  b(k-1) = b(kp)
                  b(kp) = temp
   50          continue
c
c              apply the transformation.
c
               call zaxpy(k-2,b(k),a(1,k),1,b(1),1)
               call zaxpy(k-2,b(k-1),a(1,k-1),1,b(1),1)
   60       continue
c
c           apply d inverse.
c
            ak = a(k,k)/a(k-1,k)
            akm1 = a(k-1,k-1)/a(k-1,k)
            bk = b(k)/a(k-1,k)
            bkm1 = b(k-1)/a(k-1,k)
            denom = ak*akm1 - 1.0d0
            b(k) = (akm1*bk - bkm1)/denom
            b(k-1) = (ak*bkm1 - bk)/denom
            k = k - 2
   70    continue
      go to 10
   80 continue
c
c     loop forward applying the transformations.
c
      k = 1
   90 if (k .gt. n) go to 160
         if (kpvt(k) .lt. 0) go to 120
c
c           1 x 1 pivot block.
c
            if (k .eq. 1) go to 110
c
c              apply the transformation.
c
               b(k) = b(k) + zdotu(k-1,a(1,k),1,b(1),1)
               kp = kpvt(k)
               if (kp .eq. k) go to 100
c
c                 interchange.
c
                  temp = b(k)
                  b(k) = b(kp)
                  b(kp) = temp
  100          continue
  110       continue
            k = k + 1
         go to 150
  120    continue
c
c           2 x 2 pivot block.
c
            if (k .eq. 1) go to 140
c
c              apply the transformation.
c
               b(k) = b(k) + zdotu(k-1,a(1,k),1,b(1),1)
               b(k+1) = b(k+1) + zdotu(k-1,a(1,k+1),1,b(1),1)
               kp = iabs(kpvt(k))
               if (kp .eq. k) go to 130
c
c                 interchange.
c
                  temp = b(k)
                  b(k) = b(kp)
                  b(kp) = temp
  130          continue
  140       continue
            k = k + 2
  150    continue
      go to 90
  160 continue
      return
      end
      subroutine zswap ( n, zx, incx, zy, incy )

c*********************************************************************72
c
cc ZSWAP interchanges two vectors.
c
c     jack dongarra, 3/11/78.
c
      implicit none

      integer i
      integer incx
      integer incy
      integer ix
      integer iy
      integer n
      complex*16 zx(1),zy(1),ztemp

      if(n.le.0)return
      if(incx.eq.1.and.incy.eq.1)go to 20
c
c  code for unequal increments or equal increments not equal to 1
c
      ix = 1
      iy = 1
      if(incx.lt.0)ix = (-n+1)*incx + 1
      if(incy.lt.0)iy = (-n+1)*incy + 1
      do 10 i = 1,n
        ztemp = zx(ix)
        zx(ix) = zy(iy)
        zy(iy) = ztemp
        ix = ix + incx
        iy = iy + incy
   10 continue
      return
c
c  code for both increments equal to 1
c
   20 do 30 i = 1,n
        ztemp = zx(i)
        zx(i) = zy(i)
        zy(i) = ztemp
   30 continue
      return
      end
      subroutine zsysv ( uplo, n, m, a, lda, ipiv, b, ldb, wrk, lwrk,
     &  iflag )

c*********************************************************************72
c
cc ZSYSV solves multiple symmetric linear systems of the form A*X=B.
c
      implicit none

      character UPLO*1
      integer N, M, LDA, LDB, LWRK, IFLAG, IPIV(N)
      complex*16 A(LDA,N), B(LDB,M), WRK(LWRK)
      integer j

      call ZSIFA( A,LDA, N, IPIV, IFLAG )
      if ( IFLAG.ne.0 ) stop "Error in ZSYSV (LDL' factorisation)"
      do j = 1,M
         call ZSISL( A,LDA, N, IPIV, B(1,j) )
      end do
      end
