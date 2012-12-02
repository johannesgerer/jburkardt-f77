c  SPARSKIT.F  03 December 1991
c
      subroutine amask (nrow,ncol,a,ja,ia,jmask,imask,
     *                  c,jc,ic,iw,nzmax,ierr)

c*********************************************************************72
c
cc AMASK extracts a sparse matrix from a masked input matrix.
c
c  Discussion:
c
c    AMASK extracts a sparse matrix from the input matrix
c    by looking at positions defined by the mask jmask, imask
c
c    the  algorithm is in place: c, jc, ic can be the same as
c    a, ja, ia in which cas the code will overwrite the matrix c
c    on a, ja, ia
c
c On entry:
c
c nrow  = integer. row dimension of input matrix
c ncol      = integer. Column dimension of input matrix.
c
c a,
c ja,
c ia      = matrix in Compressed Sparse Row format
c
c jmask,
c imask = matrix defining mask (pattern only) stored in compressed
c         sparse row format.
c
c nzmax = length of arrays c and jc. see ierr.
c
c On return:
c
c
c a, ja, ia and jmask, imask are unchanged.
c
c c
c jc,
c ic      = the output matrix in Compressed Sparse Row format.
c
c ierr  = integer. serving as error message.c
c         ierr = 1  means normal return
c         ierr .gt. 1 means that amask stopped when processing
c         row number ierr, because there was not enough space in
c         c, jc according to the value of nzmax.
c
c work arrays:
c
c iw      = logical work array of length ncol.
c
      double precision a(*),c(*)
      integer ia(nrow+1),ja(*),jc(*),ic(nrow+1),jmask(*),imask(nrow+1)
      logical iw(ncol)

      ierr = 0
      len = 0

      do 1 j=1, ncol
         iw(j) = .false.
 1    continue
c
c     unpack the mask for row ii in iw
c
      do ii=1, nrow
c
c     save pointer in order to be able to do things in place
c
         do k=imask(ii), imask(ii+1)-1
            iw(jmask(k)) = .true.
         end do
c
c     add umasked elemnts of row ii
c
         k1 = ia(ii)
         k2 = ia(ii+1)-1
         ic(ii) = len+1
         do k=k1,k2
            j = ja(k)
            if (iw(j)) then
               len = len+1
               if (len .gt. nzmax) then
                  ierr = ii
                  return
               end if
               jc(len) = j
               c(len) = a(k)
            end if
         end do

         do k=imask(ii), imask(ii+1)-1
            iw(jmask(k)) = .false.
         end do

      end do

      ic(nrow+1)=len+1

      return
      end
      subroutine amub (nrow,ncol,job,a,ja,ia,b,jb,ib,
     *                 c,jc,ic,nzmax,iw,ierr)

c*********************************************************************72
c
cc AMUB performs the matrix product C = A * B.
c
c on entry:
c
c nrow  = integer. The row dimension of A
c ncol  = integer. The column dimension of A
c job   = integer. Job indicator. When job = 0, only the structure
c                  (i.e. the arrays jc, ic) is computed and the
c                  real values are ignored.
c
c a,
c ja,
c ia   = Matrix A in compressed sparse row format.
c
c b,
c jb,
c ib    =  Matrix B in compressed sparse row format.
c
c nzmax = integer. The  length of the arrays c and jc.
c         amub will stop if the result matrix C  has a number
c         of elements that exceeds exceeds nzmax. See ierr.
c
c on return:
c
c c,
c jc,
c ic    = resulting matrix C in compressed sparse row sparse format.
c
c ierr      = integer. serving as error message.
c         ierr = 0 means normal return,
c         ierr .gt. 0 means that amub stopped while computing the
c         i-th row  of C with i=ierr, because the number
c         of elements in C exceeds nzmax.
c
c work arrays:
c
c iw      = integer work array of length equal to the number of
c         columns in A.
c Notes:
c
c         The column dimension of B is not needed.
c
      double precision a(*), b(*), c(*)
      integer ja(*),jb(*),jc(*),ia(nrow+1),ib(ncol+1),
     *     ic(ncol+1),iw(ncol)

      double precision scal
      logical values
      values = (job .ne. 0)
      len = 0
      ic(1) = 1
      ierr = 0
c
c  Initialize array iw.
c
      do j=1, ncol
         iw(j) = 0
      end do

      do 500 ii=1, nrow
c     row i
         do 200 ka=ia(ii), ia(ii+1)-1
          if (values) scal = a(ka)
          jj   = ja(ka)
          do 100 kb=ib(jj),ib(jj+1)-1
               jcol = jb(kb)
               jpos = iw(jcol)
               if (jpos .eq. 0) then
                  len = len+1
                  if (len .gt. nzmax) then
                     ierr = ii
                     return
                  end if
                  jc(len) = jcol
                  iw(jcol)= len
                  if (values) c(len)  = scal*b(kb)
               else
                  if (values) c(jpos) = c(jpos) + scal*b(kb)
               end if
 100          continue
 200     continue
         do 201 k=ic(ii), len
          iw(jc(k)) = 0
 201     continue
         ic(ii+1) = len+1
 500  continue
      return
      end
      subroutine amubdg(nrow,ncol,ncolb,ja,ia,jb,ib,ndegr,nnz,iw)

c*********************************************************************72
c
cc AMUBDG gets the number of nonzero elements in each row of A * B.
c
c on entry:
c
c nrow  = integer.  row dimension of matrix A
c ncol  = integer.  column dimension of matrix A = row dimension of
c                   matrix B.
c ncolb = integer. the colum dimension of the matrix B.
c
c ja, ia= row structure of input matrix A: ja = column indices of
c         the nonzero elements of A stored by rows.
c         ia = pointer to beginning of each row  in ja.
c
c jb, ib= row structure of input matrix B: jb = column indices of
c         the nonzero elements of A stored by rows.
c         ib = pointer to beginning of each row  in jb.
c
c on return:
c
c ndegr      = integer array of length nrow containing the degrees (i.e.,
c         the number of nonzeros in  each row of the matrix A * B
c
c nnz   = total number of nonzero elements found in A * B
c
c work arrays:
c
c iw      = integer work array of length ncolb.
c
      integer ja(*),jb(*),ia(nrow+1),ib(ncol+1),ndegr(nrow),iw(ncolb)

      do 1 k=1, ncolb
         iw(k) = 0
 1    continue

      do 2 k=1, nrow
         ndegr(k) = 0
 2    continue
c
c     method used: Transp(A) * A = sum [over i=1, nrow]  a(i)^T a(i)
c     where a(i) = i-th row of  A. We must be careful not to add  the
c     elements already accounted for.
c
c
      do 7 ii=1,nrow
c
c     for each row of A
c
         ldg = 0
c
c    end-of-linked list
c
         last = -1
         do 6 j = ia(ii),ia(ii+1)-1
c
c     row number to be added:
c
            jr = ja(j)
            do 5 k=ib(jr),ib(jr+1)-1
               jc = jb(k)
               if (iw(jc) .eq. 0) then
c
c     add one element to the linked list
c
                  ldg = ldg + 1
                  iw(jc) = last
                  last = jc
               end if
 5          continue
 6       continue
         ndegr(ii) = ldg
c
c     reset iw to zero
c
         do 61 k=1,ldg
            j = iw(last)
            iw(last) = 0
            last = j
 61      continue

 7    continue

      nnz = 0
      do 8 ii=1, nrow
         nnz = nnz+ndegr(ii)
 8    continue

      return
      end
      subroutine amudia (nrow,job, a, ja, ia, diag, b, jb, ib)

c*********************************************************************72
c
cc AMUDIA performs the matrix by matrix product B = A * Diag  (in place)
c
c on entry:
c
c nrow      = integer. The row dimension of A
c
c job   = integer. job indicator. Job=0 means get array b only
c         job = 1 means get b, and the integer arrays ib, jb.
c
c a,
c ja,
c ia   = Matrix A in compressed sparse row format.
c
c diag = diagonal matrix stored as a vector dig(1:n)
c
c on return:
c
c
c b,
c jb,
c ib      = resulting matrix B in compressed sparse row sparse format.
c
c Notes:
c
c 1)        The column dimension of A is not needed.
c 2)        algorithm in place (B can take the place of A).
c
      double precision a(*), b(*), diag(nrow)
      integer ja(*),jb(*), ia(nrow+1),ib(nrow+1)

      do 1 ii=1,nrow
c
c     scale each element
c
         k1 = ia(ii)
         k2 = ia(ii+1)-1
         do 2 k=k1, k2
            b(k) = a(k)*diag(ja(k))
 2       continue
 1    continue

      if (job .eq. 0) return

      ib(1) = ia(1)
      do 3 ii=1, nrow
         ib(ii) = ia(ii)
         do 31 k=ia(ii),ia(ii+1)-1
            jb(k) = ja(k)
 31      continue
 3    continue

      return
      end
      subroutine amux(n, x, y, a,ja,ia)

c*********************************************************************72
c
cc AMUX multiplies a CSR matrix A times a vector.
c
c multiplies a matrix by a vector using the dot product form
c Matrix A is stored in compressed sparse row storage.
c
c on entry:
c
c n     = row dimension of A
c x     = double precision array of length equal to the column dimension of
c         the A matrix.
c a, ja,
c    ia = input matrix in compressed sparse row format.
c
c on return:
c
c y     = double precision array of length n, containing the product y=Ax
c
      double precision  x(*), y(*), a(*)
      integer n, ja(*), ia(*)
      double precision t
      integer i, k

      do 100 i = 1,n
c
c     compute the inner product of row i with vector x
c
         t = 0.0
         do 99 k=ia(i), ia(i+1)-1
            t = t + a(k)*x(ja(k))
 99      continue
c
c     store result in y(i)
c
         y(i) = t
 100  continue

      return
      end
      subroutine amuxd (n,x,y,diag,ndiag,idiag,ioff)

c*********************************************************************72
c
cc AMUXD multiplies a DIA matrix times a vector.
c
c multiplies a matrix by a vector when the original matrix is stored
c in the diagonal storage format.
c
c on entry:
c
c n     = row dimension of A
c x     = double precision array of length equal to the column dimension of
c         the A matrix.
c ndiag  = integer. The first dimension of array adiag as declared in
c         the calling program.
c idiag  = integer. The number of diagonals in the matrix.
c diag   = double precision array containing the diagonals stored of A.
c idiag  = number of diagonals in matrix.
c diag   = double precision array of size (ndiag x idiag) containing the diagonals
c
c ioff   = integer array of length idiag, containing the offsets of the
c            diagonals of the matrix:
c          diag(i,k) contains the element a(i,i+ioff(k)) of the matrix.
c
c on return:
c
c y     = double precision array of length n, containing the product y=A*x
c
      integer n, ndiag, idiag, ioff(idiag)
      double precision x(n), y(n), diag(ndiag,idiag)
      integer j, k, io, i1, i2

      do 1 j=1, n
         y(j) = 0.0
 1    continue
      do 10 j=1, idiag
         io = ioff(j)
         i1 = max0(1,1-io)
         i2 = min0(n,n-io)
         do 9 k=i1, i2
            y(k) = y(k)+diag(k,j)*x(k+io)
 9       continue
 10   continue

      return
      end
      subroutine amuxe(n,x,y,na,ncol,a,ja)

c*********************************************************************72
c
cc AMUXE multiplies an ELL matrix times a vector.
c
c multiplies a matrix by a vector when the original matrix is stored
c in the ellpack-itpack sparse format.
c
c on entry:
c
c n     = row dimension of A
c x     = double precision array of length equal to the column dimension of
c         the A matrix.
c na    = integer. The first dimension of arrays a and ja
c         as declared by the calling program.
c ncol  = integer. The number of active columns in array a.
c         (i.e., the number of generalized diagonals in matrix.)
c a, ja = the real and integer arrays of the itpack format
c         (a(i,k),k=1,ncol contains the elements of row i in matrix
c          ja(i,k),k=1,ncol contains their column numbers)
c
c on return:
c
c y     = double precision array of length n, containing the product y=y=A*x
c
      INTEGER N,NA
      double precision x(n), y(n), a(na,*)
      integer ncol, ja(na,*)
      integer i, j

      do 1 i=1, n
         y(i) = 0.0
 1    continue
      do 10 j=1,ncol
         do 25 i = 1,n
            y(i) = y(i)+a(i,j)*x(ja(i,j))
 25      continue
 10   continue
      return
      end
      subroutine amuxj (n, x, y, jdiag, a, ja, ia)

c*********************************************************************72
c
cc AMUXJ multiplies a JAD matrix times a vector.
c
c multiplies a matrix by a vector when the original matrix is stored
c in the jagged diagonal storage format.
c
c on entry:
c
c n      = row dimension of A
c x      = double precision array of length equal to the column dimension of
c         the A matrix.
c jdiag  = integer. The number of jadded-diagonals in the data-structure.
c a      = double precision array containing the jadded diagonals of A stored
c          in succession (in decreasing lengths)
c j      = integer array containing the colum indices of the
c          corresponding elements in a.
c ia     = integer array containing the lengths of the  jagged diagonals
c
c on return:
c
c y      = double precision array of length n, containing the product y=A*x
c
c Note:
c
c Permutation related to the JAD format is not performed.
c this can be done by:
c     call permvec (n,y,y,iperm)
c after the call to amuxj, where iperm is the permutation produced
c by csrjad.
c
      integer n, jdiag, ja(*), ia(*)
      double precision x(n), y(n), a(*)
      integer i, ii, k1, len, j

      do 1 i=1, n
         y(i) = 0.0
 1    continue
      do 70 ii=1, jdiag
         k1 = ia(ii)-1
         len = ia(ii+1)-k1-1
         do 60 j=1,len
            y(j)= y(j)+a(k1+j)*x(ja(k1+j))
 60      continue
 70   continue
      return
      end
      subroutine aplb (nrow,ncol,job,a,ja,ia,b,jb,ib,
     *     c,jc,ic,nzmax,iw,ierr)

c*********************************************************************72
c
cc APLB performs the CSR matrix sum C = A + B.
c
c on entry:
c
c nrow      = integer. The row dimension of A and B
c ncol  = integer. The column dimension of A and B.
c job   = integer. Job indicator. When job = 0, only the structure
c                  (i.e. the arrays jc, ic) is computed and the
c                  real values are ignored.
c
c a,
c ja,
c ia   = Matrix A in compressed sparse row format.
c
c b,
c jb,
c ib      =  Matrix B in compressed sparse row format.
c
c nzmax      = integer. The  length of the arrays c and jc.
c         amub will stop if the result matrix C  has a number
c         of elements that exceeds exceeds nzmax. See ierr.
c
c on return:
c
c c,
c jc,
c ic      = resulting matrix C in compressed sparse row sparse format.
c
c ierr      = integer. serving as error message.
c         ierr = 0 means normal return,
c         ierr .gt. 0 means that amub stopped while computing the
c         i-th row  of C with i=ierr, because the number
c         of elements in C exceeds nzmax.
c
c work arrays:
c
c iw      = integer work array of length equal to the number of
c         columns in A.
c
      double precision a(*), b(*), c(*)
      integer ja(*),jb(*),jc(*),ia(nrow+1),ib(nrow+1),ic(nrow+1),
     *     iw(ncol)
      logical values

      values = (job .ne. 0)
      ierr = 0
      len = 0
      ic(1) = 1
      do 1 j=1, ncol
         iw(j) = 0
 1    continue

      do 500 ii=1, nrow
c     row i
         do 200 ka=ia(ii), ia(ii+1)-1
            len = len+1
            jcol    = ja(ka)
            if (len .gt. nzmax) goto 999
            jc(len) = jcol
            if (values) c(len)  = a(ka)
            iw(jcol)= len
 200     continue

         do 300 kb=ib(ii),ib(ii+1)-1
            jcol = jb(kb)
            jpos = iw(jcol)
            if (jpos .eq. 0) then
               len = len+1
               if (len .gt. nzmax) goto 999
               jc(len) = jcol
               if (values) c(len)  = b(kb)
               iw(jcol)= len
            else
               if (values) c(jpos) = c(jpos) + b(kb)
            end if
 300     continue
         do 301 k=ic(ii), len
          iw(jc(k)) = 0
 301     continue
         ic(ii+1) = len+1
 500  continue
      return
 999  ierr = ii
      return
      end
      subroutine aplb1(nrow,ncol,job,a,ja,ia,b,jb,ib,c,jc,ic,nzmax,ierr)

c*********************************************************************72
c
cc APLB1 performs the sum C = A + B for sorted CSR matrices.
c
c the difference with aplb  is that the resulting matrix is such that
c the elements of each row are sorted with increasing column indices in
c each row, provided the original matrices are sorted in the same way.
c
c on entry:
c
c nrow      = integer. The row dimension of A and B
c ncol  = integer. The column dimension of A and B.
c job   = integer. Job indicator. When job = 0, only the structure
c                  (i.e. the arrays jc, ic) is computed and the
c                  real values are ignored.
c
c a,
c ja,
c ia   = Matrix A in compressed sparse row format with entries sorted
c
c b,
c jb,
c ib      =  Matrix B in compressed sparse row format with entries sorted
c        ascendly in each row
c
c nzmax      = integer. The  length of the arrays c and jc.
c         amub will stop if the result matrix C  has a number
c         of elements that exceeds exceeds nzmax. See ierr.
c
c on return:
c
c c,
c jc,
c ic      = resulting matrix C in compressed sparse row sparse format
c         with entries sorted ascendly in each row.
c
c ierr      = integer. serving as error message.
c         ierr = 0 means normal return,
c         ierr .gt. 0 means that amub stopped while computing the
c         i-th row  of C with i=ierr, because the number
c         of elements in C exceeds nzmax.
c
c Notes:
c
c     this will not work if any of the two input matrices is not sorted
c
      double precision a(*), b(*), c(*)
      integer ja(*),jb(*),jc(*),ia(nrow+1),ib(nrow+1),ic(nrow+1)
      logical values

      values = (job .ne. 0)
      ierr = 0
      kc = 1
      ic(1) = kc

      do 6 i=1, nrow
         ka = ia(i)
         kb = ib(i)
         kamax = ia(i+1)-1
         kbmax = ib(i+1)-1
 5       continue
         if (ka .le. kamax) then
            j1 = ja(ka)
         else
            j1 = ncol+1
         end if
         if (kb .le. kbmax) then
            j2 = jb(kb)
         else
            j2 = ncol+1
         end if
c
c     three cases
c
         if (j1 .eq. j2) then
            if (values) c(kc) = a(ka)+b(kb)
            jc(kc) = j1
            ka = ka+1
            kb = kb+1
            kc = kc+1
         else if (j1 .lt. j2) then
            jc(kc) = j1
            if (values) c(kc) = a(ka)
            ka = ka+1
            kc = kc+1
         else if (j1 .gt. j2) then
            jc(kc) = j2
            if (values) c(kc) = b(kb)
            kb = kb+1
            kc = kc+1
         end if
         if (kc .gt. nzmax) goto 999
         if (ka .le. kamax .or. kb .le. kbmax) goto 5
         ic(i+1) = kc
 6    continue
      return
 999  ierr = i
      return
      end
      subroutine aplbdg (nrow,ncol,ja,ia,jb,ib,ndegr,nnz,iw)
      integer ja(*),jb(*),ia(nrow+1),ib(nrow+1),iw(ncol),ndegr(nrow)

c*********************************************************************72
c
cc APLBDG gets the number of nonzero elements in each row of A + B.
c
c on entry:
c
c nrow      = integer. The row dimension of A and B
c ncol  = integer. The column dimension of A and B.
c
c a,
c ja,
c ia   = Matrix A in compressed sparse row format.
c
c b,
c jb,
c ib      =  Matrix B in compressed sparse row format.
c
c on return:
c
c ndegr      = integer array of length nrow containing the degrees (i.e.,
c         the number of nonzeros in  each row of the matrix A + B.
c
c nnz   = total number of nonzero elements found in A * B
c
c work arrays:
c
c iw      = integer work array of length equal to ncol.
c
      do 1 k=1, ncol
         iw(k) = 0
 1    continue

      do 2 k=1, nrow
         ndegr(k) = 0
 2    continue

      do 7 ii=1,nrow
         ldg = 0
c
c    end-of-linked list
c
         last = -1
c
c     row of A
c
         do 5 j = ia(ii),ia(ii+1)-1
            jr = ja(j)
c
c     add element to the linked list
c
            ldg = ldg + 1
            iw(jr) = last
            last = jr
 5       continue
c
c     row of B
c
         do 6 j=ib(ii),ib(ii+1)-1
            jc = jb(j)
            if (iw(jc) .eq. 0) then
c
c     add one element to the linked list
c
               ldg = ldg + 1
               iw(jc) = last
               last = jc
            end if
 6       continue
c     done with row ii.
         ndegr(ii) = ldg
c
c     reset iw to zero
c
         do 61 k=1,ldg
            j = iw(last)
            iw(last) = 0
            last = j
 61      continue

 7    continue
c
      nnz = 0
      do 8 ii=1, nrow
         nnz = nnz+ndegr(ii)
 8    continue
      return
      end
      subroutine apldia (nrow, job, a, ja, ia, diag, b, jb, ib, iw)

c*********************************************************************72
c
cc APLDIA adds a diagonal matrix to a general sparse matrix:  B = A + Diag.
c
c on entry:
c
c nrow      = integer. The row dimension of A
c
c job   = integer. job indicator. Job=0 means get array b only
c         (i.e. assume that a has already been copied into array b,
c         or that algorithm is used in place. ) For all practical
c         puposes enter job=0 for an in-place call and job=1 otherwise.
c
c         Note: in case there are missing diagonal elements in A,
c         then the option job =0 will be ignored, since the algorithm
c         must modify the data structure (i.e. jb, ib) in this
c         situation.
c
c a,
c ja,
c ia   = Matrix A in compressed sparse row format.
c
c diag = diagonal matrix stored as a vector dig(1:n)
c
c on return:
c
c
c b,
c jb,
c ib      = resulting matrix B in compressed sparse row sparse format.
c
c
c iw    = integer work array of length n. On return iw will
c         contain  the positions of the diagonal entries in the
c         output matrix. (i.e., a(iw(k)), ja(iw(k)), k=1,...n,
c         are the values/column indices of the diagonal elements
c         of the output matrix. ).
c
c Notes:
c
c 1)        The column dimension of A is not needed.
c 2)        algorithm in place (b, jb, ib, can be the same as
c           a, ja, ia, on entry). See comments for parameter job.
c
c coded by Y. Saad. Latest version July, 19, 1990
c
      double precision a(*), b(*), diag(nrow)
      integer ja(*),jb(*), ia(nrow+1),ib(nrow+1), iw(*)

      logical test
c
c     copy integer arrays into b's data structure if required
c
      if (job .ne. 0) then
         nnz = ia(nrow+1)-1
         do 2  k=1, nnz
            jb(k) = ja(k)
 2       continue
         do 3 k=1, nrow+1
            ib(k) = ia(k)
 3       continue
      end if
c
c     get positions of diagonal elements in data structure.
c
      call diapos (nrow,ja,ia,iw)
c
c     count number of holes in diagonal and add diag(*) elements to
c     valid diagonal entries.
c
      icount = 0
      do 1 j=1, nrow
         if (iw(j) .eq. 0) then
            icount = icount+1
         else
            b(iw(j)) = a(iw(j)) + diag(j)
         end if
 1    continue
c
c     if no diagonal elements to insert return
c
      if (icount .eq. 0) return
c
c     shift the nonzero elements if needed, to allow for created
c     diagonal elements.
c
      ko = ib(nrow+1)+icount
c
c     copy rows backward
c
      do 5 ii=nrow, 1, -1
c
c     go through  row ii
c
         k1 = ib(ii)
         k2 = ib(ii+1)-1
         ib(ii+1) = ko
         test = (iw(ii) .eq. 0)
         do 4 k = k2,k1,-1
            j = jb(k)
            if (test .and. (j .lt. ii)) then
               test = .false.
               ko = ko - 1
               b(ko) = diag(ii)
               jb(ko) = ii
               iw(ii) = ko
            end if
            ko = ko-1
            b(ko) = a(k)
            jb(ko) = j
 4       continue
c     diagonal element has not been added yet.
         if (test) then
            ko = ko-1
            b(ko) =  diag(ii)
            jb(ko) = ii
            iw(ii) = ko
         end if
 5    continue
      ib(1) = ko
      return
      end
      subroutine aplsb (nrow,ncol,a,ja,ia,s,b,jb,ib,
     *     c,jc,ic,nzmax,iw,ierr)

c*********************************************************************72
c
cc APLSB performs the matrix linear combination C = A + s * B.
c
c on entry:
c
c nrow      = integer. The row dimension of A
c ncol  = integer. The column dimension of A
c
c a,
c ja,
c ia   = Matrix A in compressed sparse row format.
c
c s      = real. scalar factor for B.
c
c b,
c jb,
c ib      =  Matrix B in compressed sparse row format.
c
c nzmax      = integer. The  length of the arrays c and jc.
c         amub will stop if the result matrix C  has a number
c         of elements that exceeds exceeds nzmax. See ierr.
c
c on return:
c
c c,
c jc,
c ic      = resulting matrix C in compressed sparse row sparse format.
c
c ierr      = integer. serving as error message.
c         ierr = 0 means normal return,
c         ierr .gt. 0 means that amub stopped while computing the
c         i-th row  of C with i=ierr, because the number
c         of elements in C exceeds nzmax.
c
c work arrays:
c
c iw      = integer work array of length equal to the number of
c         columns in A.
c
      double precision a(*),b(*),c(*),s
      integer ja(*),jb(*),jc(*),ia(nrow+1),ib(nrow+1),
     *     ic(nrow+1),iw(ncol)

      ierr = 0
      len = 0
      ic(1) = 1

      do 1 j=1, ncol
         iw(j) = 0
 1    continue

      do 500 ii=1, nrow
c     row i
         do 200 ka=ia(ii), ia(ii+1)-1
            len = len+1
            jcol    = ja(ka)
            if (len .gt. nzmax) goto 999
            jc(len) = jcol
            c(len)  = a(ka)
            iw(jcol)= len
 200     continue
         do 300 kb=ib(ii),ib(ii+1)-1
            jcol = jb(kb)
            jpos = iw(jcol)
            if (jpos .eq. 0) then
               len = len+1
               if (len .gt. nzmax) goto 999
               jc(len) = jcol
               c(len)  = s*b(kb)
               iw(jcol)= len
            else
               c(jpos) = c(jpos) + s*b(kb)
            end if
 300     continue
         do 301 k=ic(ii), len
            iw(jc(k)) = 0
 301     continue
         ic(ii+1) = len+1
 500  continue
      return
 999  ierr = ii
      return
      end
      subroutine aplsb1 (nrow,ncol,a,ja,ia,s,b,jb,ib,c,jc,ic,
     *     nzmax,ierr)

c*********************************************************************72
c
cc APLSB1 performs the operation C = A + s * B for sorted CSR matrices.
c
c the difference with aplsb is that the resulting matrix is such that
c the elements of each row are sorted with increasing column indices in
c each row, provided the original matrices are sorted in the same way.
c
c on entry:
c
c nrow      = integer. The row dimension of A and B
c ncol  = integer. The column dimension of A and B.
c
c a,
c ja,
c ia   = Matrix A in compressed sparse row format with entries sorted
c
c s      = real. scalar factor for B.
c
c b,
c jb,
c ib      =  Matrix B in compressed sparse row format with entries sorted
c        ascendly in each row
c
c nzmax      = integer. The  length of the arrays c and jc.
c         amub will stop if the result matrix C  has a number
c         of elements that exceeds exceeds nzmax. See ierr.
c
c on return:
c
c c,
c jc,
c ic      = resulting matrix C in compressed sparse row sparse format
c         with entries sorted ascendly in each row.
c
c ierr      = integer. serving as error message.
c         ierr = 0 means normal return,
c         ierr .gt. 0 means that amub stopped while computing the
c         i-th row  of C with i=ierr, because the number
c         of elements in C exceeds nzmax.
c
c Notes:
c
c     this will not work if any of the two input matrices is not sorted
c
      double precision a(*), b(*), c(*), s
      integer ja(*),jb(*),jc(*),ia(nrow+1),ib(nrow+1),ic(nrow+1)

      ierr = 0
      kc = 1
      ic(1) = kc

      do 6 i=1, nrow
         ka = ia(i)
         kb = ib(i)
         kamax = ia(i+1)-1
         kbmax = ib(i+1)-1
 5       continue
         if (ka .le. kamax) then
            j1 = ja(ka)
         else
            j1 = ncol+1
         end if
         if (kb .le. kbmax) then
            j2 = jb(kb)
         else
            j2 = ncol+1
         end if
c
c     three cases
c
         if (j1 .eq. j2) then
            c(kc) = a(ka)+s*b(kb)
            jc(kc) = j1
            ka = ka+1
            kb = kb+1
            kc = kc+1
         else if (j1 .lt. j2) then
            jc(kc) = j1
            c(kc) = a(ka)
            ka = ka+1
            kc = kc+1
         else if (j1 .gt. j2) then
            jc(kc) = j2
            c(kc) = s*b(kb)
            kb = kb+1
            kc = kc+1
         end if
         if (kc .gt. nzmax) goto 999
         if (ka .le. kamax .or. kb .le. kbmax) goto 5
         ic(i+1) = kc
 6    continue
      return
 999  ierr = i
      return
      end
      subroutine aplsbt(nrow,ncol,a,ja,ia,s,b,jb,ib,
     *     c,jc,ic,nzmax,iw,ierr)

c*********************************************************************72
c
cc APLSBT performs the matrix sum C = A + B'.
c
c on entry:
c
c nrow      = integer. The row dimension of A and transp(B)
c ncol  = integer. The column dimension of A. Also the row
c                  dimension of B.
c
c a,
c ja,
c ia    = Matrix A in compressed sparse row format.
c
c s      = real. scalar factor for B.
c
c
c b,
c jb,
c ib      =  Matrix B in compressed sparse row format.
c
c nzmax      = integer. The  length of the arrays c, jc, and ic.
c         amub will stop if the result matrix C  has a number
c         of elements that exceeds exceeds nzmax. See ierr.
c
c on return:
c
c c,
c jc,
c ic      = resulting matrix C in compressed sparse row format.
c
c ierr      = integer. serving as error message.
c         ierr = 0 means normal return.
c         ierr = -1 means that nzmax was .lt. either the number of
c         nonzero elements of A or the number of nonzero elements in B.
c         ierr .gt. 0 means that amub stopped while computing the
c         i-th row  of C with i=ierr, because the number
c         of elements in C exceeds nzmax.
c
c work arrays:
c
c iw      = integer work array of length equal to the number of
c         columns in A.
c
c Notes:
c        It is important to note that here all of three arrays c, ic,
c        and jc are assumed to be of length nnz(c). This is because
c        the matrix is internally converted in coordinate format.
c
      double precision a(*), b(*), c(*), s
      integer ja(*),jb(*),jc(*),ia(nrow+1),ib(ncol+1),ic(*),
     *     iw(ncol)

      ierr = 0
      do 1 j=1, ncol
         iw(j) = 0
 1    continue

      nnza = ia(nrow+1)-1
      nnzb = ib(ncol+1)-1
      len = nnzb
      if (nzmax .lt. nnzb .or. nzmax .lt. nnza) then
         ierr = -1
         return
      end if
c
c     trasnpose matrix b into c
c
      ljob = 1
      ipos = 1
      call csrcsc (ncol,ljob,ipos,b,jb,ib,c,jc,ic)
      do 2 k=1,len
 2       c(k) = c(k)*s
c
c     main loop. add rows from ii = 1 to nrow.
c
         do 500 ii=1, nrow
c     iw is used as a system to recognize whether there
c     was a nonzero element in c.
            do 200 k = ic(ii),ic(ii+1)-1
               iw(jc(k)) = k
 200        continue
c
            do 300 ka = ia(ii), ia(ii+1)-1
               jcol = ja(ka)
               jpos = iw(jcol)
           if (jpos .eq. 0) then
c
c     if fill-in append in coordinate format to matrix.
c
              len = len+1
              if (len .gt. nzmax) goto 999
              jc(len) = jcol
              ic(len) = ii
              c(len)  = a(ka)
           else
c     else do addition.
              c(jpos) = c(jpos) + a(ka)
           end if
 300    continue
        do 301 k=ic(ii), ic(ii+1)-1
           iw(jc(k)) = 0
 301    continue
 500  continue
c
c     convert matrix without fill-ins into coordinate format
c
      ljob = 3
      call csrcoo (nrow,ljob,nnzb,c,jc,ic,nnzb,c,ic,jc,ierr)
      if (ierr .ne. 0) ierr = -ierr
c     convert the whole thing back to csr format.
      ljob = 1
      call coicsr (nrow,len,1,c,jc,ic,iw)
      return
 999  ierr = ii
      return
      end
      subroutine aplsca (nrow, a, ja, ia, scal,iw)

c*********************************************************************72
c
cc APLSCA adds a scalar to the diagonal entries of a sparse matrix A :=A + s I
c
c on entry:
c
c nrow      = integer. The row dimension of A
c
c a,
c ja,
c ia    = Matrix A in compressed sparse row format.
c
c scal  = real. scalar to add to the diagonal entries.
c
c on return:
c
c
c a,
c ja,
c ia      = matrix A with diagonal elements shifted (or created).
c
c iw    = integer work array of length n. On return iw will
c         contain  the positions of the diagonal entries in the
c         output matrix. (i.e., a(iw(k)), ja(iw(k)), k=1,...n,
c         are the values/column indices of the diagonal elements
c         of the output matrix. ).
c
c Notes:
c
c     The column dimension of A is not needed.
c     important: the matrix a may be expanded slightly to allow for
c     additions of nonzero elements to previously nonexisting diagonals.
c     The is no checking as to whether there is enough space appended
c     to the arrays a and ja. if not sure allow for n additional
c     elemnts.
c coded by Y. Saad. Latest version July, 19, 1990
c
      double precision a(*), scal
      integer ja(*), ia(nrow+1),iw(*)
      logical test

      call diapos (nrow,ja,ia,iw)
      icount = 0
      do 1 j=1, nrow
         if (iw(j) .eq. 0) then
            icount = icount+1
         else
            a(iw(j)) = a(iw(j)) + scal
         end if
 1    continue
c
c     if no diagonal elements to insert in data structure return.
c
      if (icount .eq. 0) return
c
c shift the nonzero elements if needed, to allow for created
c diagonal elements.
c
      ko = ia(nrow+1)+icount
c
c     copy rows backward
c
      do 5 ii=nrow, 1, -1
c
c     go through  row ii
c
         k1 = ia(ii)
         k2 = ia(ii+1)-1
         ia(ii+1) = ko
         test = (iw(ii) .eq. 0)
         do 4 k = k2,k1,-1
            j = ja(k)
            if (test .and. (j .lt. ii)) then
               test = .false.
               ko = ko - 1
               a(ko) = scal
               ja(ko) = ii
               iw(ii) = ko
            end if
            ko = ko-1
            a(ko) = a(k)
            ja(ko) = j
 4       continue
c     diagonal element has not been added yet.
         if (test) then
            ko = ko-1
            a(ko) = scal
            ja(ko) = ii
            iw(ii) = ko
         end if
 5    continue
      ia(1) = ko
      return
      end
      subroutine apmbt (nrow,ncol,job,a,ja,ia,b,jb,ib,
     *     c,jc,ic,nzmax,iw,ierr)

c*********************************************************************72
c
cc APMBT performs the matrix sum  C = A + B' or C = A - B'.
c
c on entry:
c
c nrow      = integer. The row dimension of A and transp(B)
c ncol  = integer. The column dimension of A. Also the row
c                  dimension of B.
c
c job      = integer. if job = -1, apmbt will compute C= A - transp(B)
c         (structure + values)
c         if (job .eq. 1)  it will compute C=A+transp(A)
c         (structure+ values)
c         if (job .eq. 0) it will compute the structure of
c         C= A+/-transp(B) only (ignoring all real values).
c         any other value of job will be treated as  job=1
c a,
c ja,
c ia    = Matrix A in compressed sparse row format.
c
c b,
c jb,
c ib      =  Matrix B in compressed sparse row format.
c
c nzmax      = integer. The  length of the arrays c, jc, and ic.
c         amub will stop if the result matrix C  has a number
c         of elements that exceeds exceeds nzmax. See ierr.
c
c on return:
c
c c,
c jc,
c ic      = resulting matrix C in compressed sparse row format.
c
c ierr      = integer. serving as error message.
c         ierr = 0 means normal return.
c         ierr = -1 means that nzmax was .lt. either the number of
c         nonzero elements of A or the number of nonzero elements in B.
c         ierr .gt. 0 means that amub stopped while computing the
c         i-th row  of C with i=ierr, because the number
c         of elements in C exceeds nzmax.
c
c work arrays:
c
c iw      = integer work array of length equal to the number of
c         columns in A.
c
c Notes:
c  It is important to note that here all of three arrays c, ic,
c        and jc are assumed to be of length nnz(c). This is because
c        the matrix is internally converted in coordinate format.
c
      double precision a(*), b(*), c(*)
      integer ja(*),jb(*),jc(*),ia(nrow+1),ib(ncol+1),ic(*),
     *     iw(ncol)

      logical values
      values = (job .ne. 0)

      ierr = 0
      do 1 j=1, ncol
         iw(j) = 0
 1    continue

      nnza = ia(nrow+1)-1
      nnzb = ib(ncol+1)-1
      len = nnzb
      if (nzmax .lt. nnzb .or. nzmax .lt. nnza) then
         ierr = -1
         return
      end if
c
c trasnpose matrix b into c
c
      ljob = 0
      if (values) ljob = 1
      ipos = 1
      call csrcsc (ncol,ljob,ipos,b,jb,ib,c,jc,ic)

      if (job .eq. -1) then
         do 2 k=1,len
          c(k) = -c(k)
 2       continue
      end if
c
c  main loop
c
      do 500 ii=1, nrow
         do 200 k = ic(ii),ic(ii+1)-1
            iw(jc(k)) = k
 200     continue

         do 300 ka = ia(ii), ia(ii+1)-1
            jcol = ja(ka)
            jpos = iw(jcol)
            if (jpos .eq. 0) then
c
c     if fill-in append in coordinate format to matrix.
c
               len = len+1
               if (len .gt. nzmax) goto 999
               jc(len) = jcol

               ic(len) = ii
               if (values) c(len)  = a(ka)
            else
c     else do addition.
               if (values) c(jpos) = c(jpos) + a(ka)
            end if
 300     continue
         do 301 k=ic(ii), ic(ii+1)-1
          iw(jc(k)) = 0
 301     continue
 500  continue
c
c     convert first part of matrix (without fill-ins) into coo format
c
      ljob = 2
      if (values) ljob = 3
      call csrcoo (nrow,ljob,nnzb,c,jc,ic,nnzb,c,ic,jc,ierr)
      if (ierr .ne. 0) ierr = -ierr
c     convert the whole thing back to csr format.
      ljob = 0
      if (values) ljob = 1
      call coicsr (nrow,len,ljob,c,jc,ic,iw)
      return
 999  ierr = ii
      return
      end
      subroutine assmb1 (u,nu,a,ja,ia,fu,f,nx,nelx,ijk,nodcode,
     *     node,iwk,jwk)

c*********************************************************************72
c
cc ASSMB1 assembles a finite element matrix in the CSR format.
c
c u       = unassembled matrix u(na,node,node)
c nu       = 1-st dimension of u
c a,ja,ia= assembled matrix on output
c fu       = unassembled right hand side
c f      = right hand side (global load vector) assembled
c nx     = number of nodes at input
c nelx       = number of elements at input
c ijk       = connectivity matrix: for node k, ijk(*,k) point to the
c          nodes of element k.
c node       = total number of nodal points in each element
c
c nodcode= boundary information list for each node with the
c         following meaning:
c      nodcode(i) = 0 -->  node i is internal
c      nodcode(i) = 1 -->  node i is a boundary but not a corner point
c      nodcode(i) = 2 -->  node i is a corner point (corner points
c
c x,y   = double precision arrays containing the $x$ and $y$ coordinates
c        resp. of the nodes.
c         K11, K22, and K12 at that element.
c iwk,jwk = two integer work arrays.
c ierr      = error message integer .
c        ierr = 0 --> normal return
c        ierr = 1 --> negative area encountered (due to bad
c                 numbering of nodes of an element- see
c               message printed in unit iout). not used..
c iout      = output unit (not used here).
c
      implicit double precision (a-h,o-z)
      double precision u(nu,node,node),a(*),fu(node,*),f(*)
      integer ja(*),ia(*),ijk(node,*),iwk(*),jwk(*),nodcode(*)
c     max number of nonzeros per row allowed  = 200
c
c     initialize
c
      do 100 i=1,nx
         f(i) = 0.0
 100  continue
c
c     initialize  pointer arrays.
c
      do 5 k=1,nx+1
         ia(k) = 1
         jwk(k) = 0
 5    continue
      do 6 k=1,nelx
         do 59 j=1,node
            knod = ijk(j,k)
            ia(knod) = ia(knod) + 1
 59      continue
 6    continue

      do 7 k=1, nx
         if (nodcode(k) .ge.1 ) ia(k)=ia(k)+1
 7    continue

      ksav = ia(1)
      ia(1) = 1
      do 101 j=2, nx+1
         ksavn = ia(j)
         ia(j) = ia(j-1) +  ksav
         iwk(j-1) = ia(j-1)-1
         ksav = ksavn
 101  continue
c
c     main loop
c
      do 102 nel=1, nelx
c
c     get nodal points
c
         do 120 ka=1, node
            ii = ijk(ka,nel)
            f(ii) = f(ii) + fu(ka,nel)
c
c     unpack row into jwk1
c
            irowst = ia(ii)
            ilast  = iwk(ii)
            do 109 k=irowst,ilast
               jwk(ja(k)) = k
 109        continue

            do 108 kb = 1,node
c
c     column number = jj
c
               jj = ijk(kb,nel)
               k = jwk(jj)
               if (k .eq. 0) then
                  ilast = ilast+1
                  jwk(jj) = ilast
                  ja(ilast) = jj
                  a(ilast) = u(nel,ka,kb)
               else
                  a(k) = a(k) + u(nel,ka,kb)
               end if
 108        continue
c     refresh jwk
            do 119 k=irowst,ilast
               jwk(ja(k)) = 0
 119        continue
            iwk(ii) = ilast
 120     continue

 102  continue
      return
      end
      subroutine assmbo (nx,nelx,node,ijk,nodcode,x,y,
     *             a,ja,ia,f,iwk,jwk,ierr,xyk)

c*********************************************************************72
c
cc ASSMBO assembles a finite element matrix.
c
c
c a,ja,ia= assembled matrix on output
c f      = right hand side (global load vector)
c nx     = number of nodes at input
c nelx       = number of elements at input
c ijk       = connectivity matrix: for node k, ijk(*,k) point to the
c          nodes of element k.
c node       = total number of nodal points in each element
c
c nodcode= boundary information list for each node with the
c         following meaning:
c      nodcode(i) = 0 -->  node i is internal
c      nodcode(i) = 1 -->  node i is a boundary but not a corner point
c      nodcode(i) = 2 -->  node i is a corner point (corner points
c
c mp      = group material number for each element.
c
c x,y   = double precision arrays containing the $x$ and $y$ coordinates
c        resp. of the nodes.
c iwk,jwk = two integer work arrays.
c ierr      = error message integer .
c        ierr = 0 --> normal return
c        ierr = 1 --> negative area encountered (due to bad
c                 numbering of nodes of an element- see
c               message printed in unit iout ).
c iout      = output unit (not used here).
c
c xyk      = routine defining the material properties at each
c         element. Form:
c       call xyk(nel,xyke,x,y,ijk,node) with on return
c         xyke =  material constant matrices.
c         for each element nel, xyke(1,nel),xyke(2,nel)
c         and xyke(3,nel) represent the constants
c         K11, K22, and K12 at that element.
c
      implicit double precision  (a-h,o-z)
      dimension a(1),ijk(node,1),x(1),y(1),f(1),ske(3,3),fe(3),
     *      xe(3),ye(3),xyke(2,2),iwk(1),jwk(1)
      integer ia(1), ja(1), nodcode(1)
      external xyk
c max number of nonzeros allowed  = 200
c
c   initialize
c
      do 100 i=1,nx
 100          f(i) = 0.0
c initialize  pointer arrays.
      do 5 k=1,nx+1
      ia(k) = 1
      jwk(k) = 0
 5      continue
      do 6 k=1,nelx
      do 59 j=1,node
      knod = ijk(j,k)
 59      ia(knod) = ia(knod) + 1
 6      continue

      do 7 k=1, nx
      if (nodcode(k) .ge.1 ) ia(k)=ia(k)+1
 7      continue

      ksav = ia(1)
      ia(1) = 1
      do 101 j=2, nx+1
            ksavn = ia(j)
            ia(j) = ia(j-1) +  ksav
            iwk(j-1) = ia(j-1)-1
 101            ksav = ksavn
c
c main loop
c
      do 102 nel=1, nelx
c
c get coordinates of nodal points
c
      do 104 i=1, node
      j = ijk(i,nel)
      xe(i) = x(j)
      ye(i) = y(j)
 104      continue
c
c compute determinant
c
       det=xe(2)*(ye(3)-ye(1))+xe(3)*(ye(1)-ye(2))+xe(1)*(ye(2)-ye(3))
c
c set material properties
c
      call xyk(nel,xyke,x,y,ijk,node)
c
c construct element stiffness matrix
c
      ierr = 0

      call estif3(nel,ske,fe,det,xe,ye,xyke,ierr)
      if (ierr .ne. 0) return
c
c assemble: add element stiffness matrix to global matrix
c
      do 120 ka=1, node
            ii = ijk(ka,nel)
          f(ii) = f(ii) + fe(ka)
c
c unpack row into jwk1
c
         irowst = ia(ii)
         ilast  = iwk(ii)
           do 109 k=irowst,ilast
         jwk(ja(k)) = k
 109         continue
c
          do 108 kb = 1,node
c
c column number = jj
c
            jj = ijk(kb,nel)
           k = jwk(jj)
          if (k .eq. 0) then
             ilast = ilast+1
             jwk(jj) = ilast
             ja(ilast) = jj
             a(ilast) = ske(ka,kb)
          else
              a(k) = a(k) + ske(ka,kb)
          end if
 108         continue
c refresh jwk
           do 119 k=irowst,ilast
         jwk(ja(k)) = 0
 119         continue
           iwk(ii) = ilast
 120    continue
c
 102      continue
        return
      end
      subroutine atmux (n, x, y, a, ja, ia)

c*********************************************************************72
c
cc ATMUX computes A' * x for a CSR matrix A.
c
c multiplies the transpose of a matrix by a vector when the original
c matrix is stored in compressed sparse row storage. Can also be
c viewed as the product of a matrix by a vector when the original
c matrix is stored in the compressed sparse column format.
c
c on entry:
c
c n     = row dimension of A
c x     = double precision array of length equal to the column dimension of
c         the A matrix.
c a, ja,
c    ia = input matrix in compressed sparse row format.
c
c on return:
c
c y     = double precision array of length n, containing the product y=transp(A)*x
c
c
      double precision x(*), y(*), a(*)
      integer n, ia(*), ja(*)
      integer i, k
c
c     zero out output vector
c
      do 1 i=1,n
         y(i) = 0.0
 1    continue
c
c loop over the rows
c
      do 100 i = 1,n
         do 99 k=ia(i), ia(i+1)-1
            y(ja(k)) = y(ja(k)) + x(i)*a(k)
 99      continue
 100  continue

      return
      end
      subroutine blkchk (nrow,ja,ia,nblk,imsg)

c*********************************************************************72
c
cc BLKCHK checks whether the input matrix is a block matrix.
c
c BLKCHK checks whether the input matrix is a block
c matrix with block size of nblk. A block matrix is one which is
c comprised of small square dense blocks. If there are zero
c elements within the square blocks and the data structure
c takes them into account then blkchk may fail to find the
c correct block size.
c
c on entry
c
c nrow      = integer equal to the row dimension of the matrix.
c ja    = integer array containing the column indices of the entries
c         nonzero entries of the matrix stored by row.
c ia    = integer array of length nrow + 1 containing the pointers
c         beginning of each row in array ja.
c
c nblk  = integer containing the value of nblk to be checked.
c
c on return
c
c
c imsg  = integer containing a message  with the following meaning.
c          imsg = 0 means that the output value of nblk is a correct
c                   block size. nblk .lt. 0 means nblk not correct
c                   block size.
c          imsg = -1 : nblk does not divide nrow
c          imsg = -2 : a starting element in a row is at wrong position
c             (j .ne. mult*nblk +1 )
c          imsg = -3 : nblk does divide a row length -
c          imsg = -4 : an element is isolated outside a block or
c             two rows in same group have different lengths
c
c           Y. Saad, Sep. 21 1989                                      c
c
      integer ia(nrow+1),ja(*)
c
c first part of code will find candidate block sizes.
c this is not guaranteed to work . so a check is done at the end
c the criterion used here is a simple one:
c scan rows and determine groups of rows that have the same length
c and such that the first column number and the last column number
c are identical.
c
      imsg = 0
      if (nblk .le. 1) return
      nr = nrow/nblk
      if (nr*nblk .ne. nrow) goto 101
c   main loop
      irow = 1
      do 20 ii=1, nr
c     i1= starting position for group of nblk rows in original matrix
         i1 = ia(irow)
         j2 = i1
c     lena = length of each row in that group  in the original matrix
         lena = ia(irow+1)-i1
c     len = length of each block-row in that group in the output matrix
         len = lena/nblk
         if (len* nblk .ne. lena) goto 103
c
c     for each row
c
         do 6 i = 1, nblk
            irow = irow + 1
            if (ia(irow)-ia(irow-1) .ne. lena ) goto 104
c
c     for each block
c
            do 7 k=0, len-1
               jstart = ja(i1+nblk*k)-1
               if ( (jstart/nblk)*nblk .ne. jstart) goto 102
c
c     for each column
c
               do 5 j=1, nblk
                  if (jstart+j .ne. ja(j2) )  goto 104
                  j2 = j2+1
 5             continue
 7          continue
 6       continue
 20   continue
c     went through all loops successfully:
      return
 101  imsg = -1
      return
 102  imsg = -2
      return
 103  imsg = -3
      return
 104  imsg = -4
      return
      end
      subroutine blkfnd (nrow,ja,ia,nblk)

c*********************************************************************72
c
cc BLKFND determines the block structure of a matrix.
c
c This routine attemptps to determine whether or not  the input
c matrix has a block structure and finds the blocks size
c if it does. A block matrix is one which is
c comprised of small square dense blocks. If there are zero
c elements within the square blocks and the original data structure
c takes these zeros into account then blkchk may fail to find the
c correct block size.
c
c on entry
c
c nrow      = integer equal to the row dimension of the matrix.
c ja    = integer array containing the column indices of the entries
c         nonzero entries of the matrix stored by row.
c ia    = integer array of length nrow + 1 containing the pointers
c         beginning of each row in array ja.
c
c nblk  = integer containing the assumed value of nblk if job = 0
c
c on return
c
c nblk  = integer containing the value found for nblk when job = 1.
c         if imsg .ne. 0 this value is meaningless however.
c
c
c           Y. Saad, Sep. 21 1989                                      c
c
      integer ia(nrow+1),ja(*)
c
c first part of code will find candidate block sizes.
c criterion used here is a simple one: scan rows and  determine groups
c of rows that have the same length and such that the first column
c number and the last column number are identical.
c
      minlen = ia(2)-ia(1)
      irow   = 1
      do 1 i=2,nrow
         len = ia(i+1)-ia(i)
         if (len .lt. minlen) then
            minlen = len
            irow = i
         end if
 1    continue
c
c  candidates are all dividers of minlen
c
      nblk = 1
      if (minlen .le. 1) return

      do 99 iblk = minlen, 1, -1
         if (mod(minlen,iblk) .ne. 0) goto 99
         len = ia(2) - ia(1)
         len0 = len
         jfirst = ja(1)
         jlast = ja(ia(2)-1)
         do 10 jrow = irow+1,irow+nblk-1
            i1 = ia(jrow)
            i2 = ia(jrow+1)-1
            len = i2+1-i1
            jf = ja(i1)
            jl = ja(i2)
            if (len .ne. len0 .or. jf .ne. jfirst .or.
     *           jl .ne. jlast) goto 99
 10      continue
c
c     check for this candidate
c
         call blkchk (nrow,ja,ia,iblk,imsg)
         if (imsg .eq. 0) then
c
c     block size found
c
            nblk = iblk
            return
         end if
 99   continue
      end
      subroutine bndcsr (n,abd,nabd,lowd,ml,mu,a,ja,ia,len,ierr)

c*********************************************************************72
c
cc BNDCSR converts Banded Linpack format to Compressed Sparse Row format.
c
c Banded (Linpack ) format   to    Compressed Sparse Row  format.
c
c on entry:
c
c n      = integer,the actual row dimension of the matrix.
c
c nabd  = first dimension of array abd.
c
c abd   = double precision array containing the values of the matrix stored in
c         banded form. The j-th column of abd contains the elements
c         of the j-th column of  the original matrix,comprised in the
c         band ( i in (j-ml,j+mu) ) with the lowest diagonal located
c         in row lowd (see below).
c
c lowd  = integer. this should be set to the row number in abd where
c         the lowest diagonal (leftmost) of A is located.
c         lowd should be s.t.  ( 1  .le.  lowd  .le. nabd).
c         The routines dgbco, ... of linpack use lowd=2*ml+mu+1.
c
c ml      = integer. equal to the bandwidth of the strict lower part of A
c mu      = integer. equal to the bandwidth of the strict upper part of A
c         thus the total bandwidth of A is ml+mu+1.
c         if ml+mu+1 is found to be larger than nabd then an error
c         message is set. see ierr.
c
c len   = integer. length of arrays a and ja. bndcsr will stop if the
c         length of the arrays a and ja is insufficient to store the
c         matrix. see ierr.
c
c on return:
c
c a,
c ja,
c ia    = input matrix stored in compressed sparse row format.
c
c lowd  = if on entry lowd was zero then lowd is reset to the default
c         value ml+mu+l.
c
c ierr  = integer. used for error message output.
c         ierr .eq. 0 :means normal return
c         ierr .eq. -1 : means invalid value for lowd.
c        ierr .gt. 0 : means that there was not enough storage in a and ja
c         for storing the ourput matrix. The process ran out of space
c         (as indicated by len) while trying to fill row number ierr.
c         This should give an idea of much more storage might be required.
c         Moreover, the first irow-1 rows are correctly filled.
c
c notes:  the values in abd found to be equal to zero
c         (actual test: if (abd(...) .eq. 0.0) are removed.
c         The resulting may not be identical to a csr matrix
c         originally transformed to a bnd format.
c
      double precision a(*),abd(nabd,*), t
      integer ia(n+1),ja(*)

      ierr = 0

      if (lowd .gt. nabd .or. lowd .le. 0) then
         ierr = -1
         return
      end if

      ko = 1
      ia(1) = 1
      do 30 irow=1,n

         i = lowd
          do  20 j=irow-ml,irow+mu
             if (j .le. 0 ) goto 19
             if (j .gt. n) goto 21
             t = abd(i,j)
             if (t .eq. 0.0) goto 19
             if (ko .gt. len) then
               ierr = irow
               return
            end if
            a(ko) = t
            ja(ko) = j
            ko = ko+1
 19         i = i-1
 20      continue
c     end for row irow
 21      ia(irow+1) = ko
 30   continue
      return
      end
      subroutine bound (nx,nelx,ijk,nodcode,node,nint,iperm,
     *             x,y,wk,iwk)

c*********************************************************************72
c
cc BOUND counts the number of boundary points.
c
c  this routine counts the number of boundary points and
c reorders the points in such a way that the boundary nodes
c are last.
c
c nx, nelx, ijk, nodcode, node: see other routines
c nint = on return the number of points on the boundary
c iperm = permutation array from old orderin to new ordering,
c iwk   = reverse permutation array or return.
c wk      = double precision work array
c On return
c x, y, nodecode, are permuted
c ijk  is updated according to new oerdering.
c nint = number of interior points.
c
      implicit double precision  (a-h,o-z)
      dimension ijk(node,1),x(1),y(1),wk(1),iwk(1),iperm(1),
     *              nodcode(1)
c max number of nonzeros allowed  = 200
c put all boundary points at the end, backwards
      nint = 1
      nbound = nx
      do 1 j=1, nx
      if (nodcode(j) .eq. 0) then
        iperm(nint) = j
        nint = nint+1
      else
      iperm(nbound) = j
      nbound = nbound-1
       end if
 1      continue

      nint = nint-1
c
c permute x's
c
      do 2 k=1, nx
      wk(k) = x(k)
 2      continue
      do 3 k=1,nx
        x(k) = wk(iperm(k))
 3      continue
c
c permute the y's
c
      do 4 k=1, nx
      wk(k) = y(k)
 4      continue
      do 5 k=1, nx
      y(k) = wk(iperm(k))
 5      continue
c
c permute the boundary information
c
      do 6 k=1, nx
      iwk(k) = nodcode(k)
 6      continue
      do 7 k=1,nx
      nodcode(k) = iwk(iperm(k))
 7      continue
c
c get reverse permutation
c
      do 8 k=1, nx
        iwk(iperm(k)) = k
 8      continue
c
c update the elements connectivity matrix
c
      do 10 nel = 1, nelx
         do 9 j=1, node
            knod = ijk(j,nel)
              ijk(j,nel) = iwk(knod)
 9         continue
 10      continue
      return
      end
      subroutine bsort2 (w, ind, n, ncut)

c*********************************************************************72
c
cc BSORT2 returns the NCUT largest elements of an array, using bubble sort.
c
c simple bubble sort for getting the ncut largest
c elements in modulus, in array w. ind is sorted accordingly.
c (Ought to be replaced by a more efficient sort especially
c if ncut is not that small).
c
      integer n, ncut, ind(*)
      double precision w(*)
      logical test
      integer i, j, iswp
      double precision wswp

      i = 1
 1    test = .false.
      do 2 j = n-1,i,-1
         if (abs(w(j+1))  .gt. abs(w(j)) ) then
c  swap.
            wswp = w(j)
            w(j) = w(j+1)
            w(j+1) = wswp
c
c  reorder original ind array accordingly.
c
            iswp = ind(j)
            ind(j) = ind(j+1)
            ind(j+1) = iswp
c  set indicator that sequence is still unsorted--
            test = .true.
         end if
 2    continue
      i = i+ 1
      if (test .and. i .le. ncut) goto 1
      return
      end
      subroutine bsrcsr (n, nblk, na, a, ja, ia, ao, jao, iao)

c*********************************************************************72
c
cc BSRCSR converts Block Sparse Row to Compressed Sparse Row (CSR) format.
c
c this routine converts a matrix stored in block-reduced
c a, ja, ia format to the general sparse row a, ja, ia format.
c A matrix that has a block structure is a matrix whose entries
c are blocks of the same size nblk (e.g. 3 x 3). Then it is often
c preferred to work with the reduced graph of the matrix, i.e.,
c Instead of storing one element at a time one can store the whole
c block. In this storage scheme a row of the array a will
c hold the nblk**2 entries of a block.
c
c on entry:
c
c n      = integer, the actual row dimension of the matrix.
c nblk  = integer equal to the dimension of each block.
c         nblk must divide n.
c na      = first dimension of array a as declared in calling program
c a      = double precision array containing the values of the matrix. For details
c         on the format see below. Each row of a contains the nblk x nblk
c         block matrix unpacked column-wise (this allows the user to
c         declare the array a as a(na,nblk,nblk) on entry if desired).
c         the block rows are stored in sequence just as for the compressed
c         sparse row format.
c ja      = integer array of length n/nblk. ja(k) contains the column index
c         of the leading element, i.e., the element (1,1) of the block
c         that is held in the row a(k,*) of the value array.
c ia    = integer array of length n/nblk+1. ia(i) points to the beginning
c        of block row number i in the arrays a and ja.
c
c on return:
c
c ao, jao,
c     iao = matrix stored in compressed sparse row format.
c
c Notes: this code is not in place.
c
c
c general picture: (nblk = 2)
c     --- A ---                                --- JA --  -- IA --
c A=  x x x x   1st block in block row 1           x         x
c     x x x x  2-nd block in block row 1           x
c     . . . .                                      .
c     x x x x  last block in block row 1           x
c     -------                                     ---
c     x x x x   1st block in block row 2           x          x
c     x x x x  2-nd block in block row 2           x
c     . . . .                                      x
c     x x x x   last block in block row 2          x
c     -------                                     ---
c     .......                                     ...         .
c     -------                                     ---
c     x x x x   1st block in block row n/nblk      x          x
c     x x x x  2-nd block in block row n/nblk      x
c     . . . .                                      x
c     x x x x  last block in block row n/nblk      x
c     -------                                     ---
c                                               end + 1       x
c
c
c example  with nblk = 2:
c
c
c             1   2   0   0   3   4
c             5   6   0   0   7   8
c             0   0   9  10  11  12
c             0   0  13  14  15  16
c             17 18   0   0   0   0
c             22 23   0   0   0   0
c THEN:
c
c  ---- A ----                                     -- JA --   -- IA --
c
c  1   5   2  6  Block row 1 (2 block matrices)      | 1  <--- | 1
c  3   7   4  8                                      | 5       |
c  ------------                                      | --      |
c  9  13  10 14  block row 2 (2 block matrices)      | 3  <--- | 3
c 11  15  12 16                                      | 5       |
c  ------------                                      | --      |
c 17  22  18 23  Block row 3 (1 block matrix)        | 1  <--- | 5
c  ------------                                      | --      |
c                                                   end+1 <--- | 6
c
c JA  =  1  5 | 3  5 | 1       column numbers of (1,1) entries of blocks
c IA  =  1      3      5  6    pointers to beginnings of BLOCK-rows
c
c
c get ia, ja data structure for output matrix
c
      double precision a(na,*), ao(*)
      integer ia(*), ja(*), jao(*), iao(n+1)

      nr = n/nblk
      do 1 k=1,n+1
         iao(k) = 0
 1    continue

      irow = 0
      krow = 1
      do 2 ii=1, nr
c     nr is the dimension of the reduced matrix.
         i1 = ia(ii)
         i2 = ia(ii+1)-1
c     create nblk rows for each k
         do 23 i=1,nblk
            do 21 k=i1, i2
               jst = ja(k)-1
               do 22  j=1,nblk
                  ij = (j-1)*nblk + i
                  ao(krow) = a(k,ij)
                  jao(krow) = jst+j
                  krow = krow+1
 22            continue
 21          continue
          iao(irow+i) = krow
 23      continue
         irow = irow + nblk
 2    continue
      do 3 jj=1,n
         j = n-jj+1
         iao(j+1)=iao(j)
 3    continue
      iao(1) = 1
      return
      end
      subroutine bsten (nx,ny,nz,kx,ky,kz,nfree,stencil,h)

c*********************************************************************72
c
cc BSTEN calculates block stencil values.
c
c  This routine calculates the correct block-stencil values for
c     centered difference discretization of the elliptic operator
c     (block version of stencil)
c
c L u = delx( a delx u ) + dely ( b dely u) + delz ( c delz u ) +
c       d delx ( u ) + e dely (u) + f delz( u ) + g u
c
c   For 2-D problems the discretization formula that is used is:
c
c h**2 * Lu == a(i+1/2,j)*{u(i+1,j) - u(i,j)} +
c             a(i-1/2,j)*{u(i-1,j) - u(i,j)} +
c              b(i,j+1/2)*{u(i,j+1) - u(i,j)} +
c              b(i,j-1/2)*{u(i,j-1) - u(i,j)} +
c              (h/2)*d(i,j)*{u(i+1,j) - u(i-1,j)} +
c              (h/2)*e(i,j)*{u(i,j+1) - u(i,j-1)} +
c              (h/2)*e(i,j)*{u(i,j+1) - u(i,j-1)} +
c              (h**2)*g(i,j)*u(i,j)
c
c      implicit double precision (a-h,o-z)
      double precision stencil(7,*)
      double precision cntr(225), coeff(225),h, hhalf, x, y, z

      if (nfree .gt. 15) then
        WRITE(*,*)'BSTEN  - FATAL ERROR'
        WRITE(*,*)'         Input value of NFREE is greater than 15.'
        STOP
      end if

      nfree2 = nfree*nfree
      do 200 k=1, nfree2
         cntr(k) = 0.0
         do 199 i=1,7
            stencil(i,k) = 0.0
 199     continue
 200  continue

      hhalf = h*0.5
      h2 = h*h
      x = h*dble(kx)
      y = h*dble(ky)
      z = h*dble(kz)
c differentiation wrt x:
      call afunbl(nfree,x+hhalf,y,z,coeff)
      do 1 k=1, nfree2
      stencil(3,k) = stencil(3,k) + coeff(k)
      cntr(k) = cntr(k) + coeff(k)
 1    continue

      call afunbl(nfree,x-hhalf,y,z,coeff)
      do 2 k=1, nfree2
         stencil(2,k) = stencil(2,k) + coeff(k)
         cntr(k) = cntr(k) + coeff(k)
 2    continue

      call dfunbl(nfree,x,y,z,coeff)
      do 3 k=1, nfree2
         stencil(3,k) = stencil(3,k) + coeff(k)*hhalf
         stencil(2,k) = stencil(2,k) - coeff(k)*hhalf
 3    continue
      if (ny .le. 1) goto 99
c
c differentiation wrt y:
c
      call bfunbl(nfree,x,y+hhalf,z,coeff)
      do 4 k=1,nfree2
         stencil(5,k) = stencil(5,k) + coeff(k)
         cntr(k) = cntr(k) + coeff(k)
 4    continue
c
      call bfunbl(nfree,x,y-hhalf,z,coeff)
      do 5 k=1, nfree2
         stencil(4,k) = stencil(4,k) + coeff(k)
         cntr(k) = cntr(k) + coeff(k)
 5    continue
c
      call efunbl(nfree,x,y,z,coeff)
      do 6 k=1, nfree2
         stencil(5,k) = stencil(5,k) + coeff(k)*hhalf
         stencil(4,k) = stencil(4,k) - coeff(k)*hhalf
 6    continue
      if (nz .le. 1) goto 99
c
c differentiation wrt z:
c
      call cfunbl(nfree,x,y,z+hhalf,coeff)
      do 7 k=1, nfree2
         stencil(7,k) = stencil(7,k) + coeff(k)
         cntr(k) = cntr(k) + coeff(k)
 7    continue
c
      call cfunbl(nfree,x,y,z-hhalf,coeff)
      do 8 k=1, nfree2
         stencil(6,k) = stencil(6,k) + coeff(k)
         cntr(k) = cntr(k) + coeff(k)
 8    continue
c
      call ffunbl(nfree,x,y,z,coeff)
      do 9 k=1, nfree2
         stencil(7,k) = stencil(7,k) + coeff(k)*hhalf
         stencil(6,k) = stencil(6,k) - coeff(k)*hhalf
 9    continue
c
c discretization of  product by g:
c
 99   call gfunbl(nfree,x,y,z,coeff)
      do 10 k=1, nfree2
         stencil(1,k) = h2*coeff(k) - cntr(k)
 10   continue
c
      return
      end
      subroutine checkref(nx,nelx,ijk,node,nodcode,
     *               nbound,  nxnew,nelxnew)

c*********************************************************************72
c
cc CHECKREF returns the expected number of new nodes and elements.
c
c returns the expected the new number of nodes and
c elemnts if refall is applied to current grid once.
c
c nx      = number of nodes at input
c nelx      = number of elements at input
c ijk      = connectivity matrix: for node k, ijk(*,k) point to the
c         nodes of element k.
c nbound  = number of boundary points on entry - enter zero if
c           unknown
c
c nodcode= boundary information list for each node with the
c         following meaning:
c      nodcode(i) = 0 -->  node i is internal
c      nodcode(i) = 1 -->  node i is a boundary but not a corner point
c      nodcode(i) = 2 -->  node i is a corner point.
c
c nxnew  = new number of nodes if refall were to be applied
c nelxnew = same for nelx.
c
       integer ijk(node,1),nodcode(nx)

      nelxnew = nelx*4
c
c count the number of boundary nodes
c
      if (nbound .ne. 0) goto 2
      do 1 j=1, nx
      if (nodcode(j) .ge. 1) nbound = nbound+1
 1      continue
c number of edges=[3*(number of elmts) + number of bound nodes ]/ 2
 2      continue
      nxnew = nx + (3*nelx+nbound)/2
      nbound = 2*nbound
      return
      end
      subroutine chkelmt (nx, x, y, nelx, ijk, node)

c*********************************************************************72
c
cc CHKELMT checks the labeling within each element and reorders if necessary.
c
c  CHKELMT checks the labeling within each element and reorders
c the nodes if they are not correctly ordered.
c
      implicit double precision (a-h,o-z)
      dimension ijk(node,1),x(1),y(1)

      do 1 nel =1, nelx
       det = x(ijk(2,nel))*(y(ijk(3,nel))-y(ijk(1,nel)))+
     *        x(ijk(3,nel))*(y(ijk(1,nel))-y(ijk(2,nel)))+
     *        x(ijk(1,nel))*(y(ijk(2,nel))-y(ijk(3,nel)))
c
c if determinant negative exchange last two nodes of elements.
c
      if (det .lt. 0.0) then
          j = ijk(2,nel)
          ijk(2,nel) = ijk(3,nel)
          ijk(3,nel) = j
      end if
 1      continue

      return
      end
      subroutine cnrms   (nrow, nrm, a, ja, ia, diag)

c*********************************************************************72
c
cc CNRMS gets the norms of each column of A.
c
c gets the norms of each column of A. (choice of three norms)
c
c on entry:
c
c nrow      = integer. The row dimension of A
c
c nrm   = integer. norm indicator. nrm = 1, means 1-norm, nrm =2
c                  means the 2-nrm, nrm = 0 means max norm
c
c a,
c ja,
c ia   = Matrix A in compressed sparse row format.
c
c on return:
c
c
c diag = double precision vector of length nrow containing the norms
c
      double precision a(*), diag(nrow)
      integer ja(*), ia(nrow+1)

      do 10 k=1, nrow
         diag(k) = 0.0
 10   continue
      do 1 ii=1,nrow
         k1 = ia(ii)
         k2 = ia(ii+1)-1
         do 2 k=k1, k2
            j = ja(k)
c     update the norm of each column
            if (nrm .eq. 0) then
               diag(j) = max(diag(j),abs(a(k) ) )
            elseif (nrm .eq. 1) then
               diag(j) = diag(j) + abs(a(k) )
            else
               diag(j) = diag(j)+a(k)**2
            end if
 2       continue
 1    continue
      if (nrm .ne. 2) return
      do 3 k=1, nrow
         diag(k) = sqrt(diag(k))
 3    continue
      return
      end
      subroutine coicsr (n,nnz,job,a,ja,ia,iwk)

c*********************************************************************72
c
cc COICSR converts COO to CSR in place.
c
c IN-PLACE coo-csr conversion routine.
c
c COICSR converts a matrix stored in coordinate format into
c the csr format. The conversion is done in place in that the arrays
c a,ja,ia of the result are overwritten onto the original arrays.
c
c on entry:
c
c n      = integer. row dimension of A.
c nnz      = integer. number of nonzero elements in A.
c job   = integer. Job indicator. when job=1, the real values in a are
c         filled. Otherwise a is not touched and the structure of the
c         array only (i.e. ja, ia)  is obtained.
c a      = double precision array of size nnz (number of nonzero elements in A)
c         containing the nonzero elements
c ja      = integer array of length nnz containing the column positions
c         of the corresponding elements in a.
c ia      = integer array of length nnz containing the row positions
c         of the corresponding elements in a.
c iwk      = integer work array of length n.
c on return:
c
c a
c ja
c ia      = contains the compressed sparse row data structure for the
c         resulting matrix.
c Note:
c
c         the entries of the output matrix are not sorted (the column
c         indices in each are not in increasing order) use coocsr
c         if you want them sorted.
c
c  Coded by Y. Saad, Sep. 26 1989
c
      integer ia(nnz),ja(nnz),iwk(n)
      double precision a(*)
      double precision t,tnext
      logical values

      values = (job .eq. 1)
c
c find pointer array for resulting matrix.
c
      do 35 i=1,n+1
         iwk(i) = 0
 35   continue
      do 4 k=1,nnz
         i = ia(k)
         iwk(i+1) = iwk(i+1)+1
 4    continue

      iwk(1) = 1
      do 44 i=2,n
         iwk(i) = iwk(i-1) + iwk(i)
 44   continue
c
c  loop for a cycle in chasing process.
c
      init = 1
      k = 0
 5    if (values) t = a(init)
      i = ia(init)
      j = ja(init)
      ia(init) = -1

 6    k = k+1
c
c  current row number is i.  determine  where to go.
c
      ipos = iwk(i)
c
c  save the chased element.
c
      if (values) tnext = a(ipos)
      inext = ia(ipos)
      jnext = ja(ipos)
c
c  then occupy its location.
c
      if (values) a(ipos)  = t
      ja(ipos) = j
c
c  update pointer information for next element to come in row i.
c
      iwk(i) = ipos+1
c
c  determine  next element to be chased,
c
      if (ia(ipos) .lt. 0) goto 65
      t = tnext
      i = inext
      j = jnext
      ia(ipos) = -1
      if (k .lt. nnz) goto 6
      goto 70
 65   init = init+1
      if (init .gt. nnz) goto 70
      if (ia(init) .lt. 0) goto 65
c
c  restart chasing.
c
      goto 5
 70   do 80 i=1,n
         ia(i+1) = iwk(i)
 80   continue
      ia(1) = 1
      return
      end
      subroutine coocsr(nrow,nnz,a,ir,jc,ao,jao,iao)

c*********************************************************************72
c
cc COOCSR converts COO to CSR.
c
c converts a matrix that is stored in coordinate format
c  a, ir, jc into a row general sparse ao, jao, iao format.
c
c on entry:
c
c nrow      = dimension of the matrix
c nnz      = number of nonzero elements in matrix
c a,
c ir,
c jc    = matrix in coordinate format. a(k), ir(k), jc(k) store the nnz
c         nonzero elements of the matrix with a(k) = actual real value of
c         the elements, ir(k) = its row number and jc(k) = its column
c        number. The order of the elements is arbitrary.
c
c on return:
c
c ir       is destroyed
c
c ao, jao, iao = matrix in general sparse matrix format with ao
c       continung the real values, jao containing the column indices,
c      and iao being the pointer to the beginning of the row,
c      in arrays ao, jao.
c
c Notes:
c This routine is NOT in place.  See coicsr
c
      double precision a(*),ao(*),x
      integer ir(*),jc(*),jao(*),iao(*)

      do 1 k=1,nrow+1
         iao(k) = 0
 1    continue
c determine row-lengths.
      do 2 k=1, nnz
         iao(ir(k)) = iao(ir(k))+1
 2    continue
c starting position of each row.
      k = 1
      do 3 j=1,nrow+1
         k0 = iao(j)
         iao(j) = k
         k = k+k0
 3    continue
c go through the structure  once more. Fill in output matrix.
      do 4 k=1, nnz
         i = ir(k)
         j = jc(k)
         x = a(k)
         iad = iao(i)
         ao(iad) =  x
         jao(iad) = j
         iao(i) = iad+1
 4    continue
c shift back iao
      do 5 j=nrow,1,-1
         iao(j+1) = iao(j)
 5    continue
      iao(1) = 1
      return
      end
      subroutine cooell(n,nnz,a,ja,ia,ac,jac,nac,ner,ncmax,ierr)

c*********************************************************************72
c
cc COOELL converts coordinate format to Ellpack/Itpack format.
c
c   DATE WRITTEN: June 4, 1989.
c
c   PURPOSE
c
c  COOELL takes a sparse matrix in coordinate format and
c  converts it into the Ellpack-Itpack storage.
c
c  Example:
c
c       (   11   0   13    0     0     0  )
c       |   21  22    0   24     0     0  |
c       |    0  32   33    0    35     0  |
c   A = |    0   0   43   44     0    46  |
c       |   51   0    0   54    55     0  |
c       (   61  62    0    0    65    66  )
c
c   Coordinate storage scheme:
c
c    A  = (11,22,33,44,55,66,13,21,24,32,35,43,46,51,54,61,62,65)
c    IA = (1, 2, 3, 4, 5, 6, 1, 2, 2, 3, 3, 4, 4, 5, 5, 6, 6, 6 )
c    JA = ( 1, 2, 3, 4, 5, 6, 3, 1, 4, 2, 5, 3, 6, 1, 4, 1, 2, 5)
c
c   Ellpack-Itpack storage scheme:
c
c       (   11  13    0    0   )          (   1   3   *    *  )
c       |   22  21   24    0   |          |   2   1   4    *  |
c  AC = |   33  32   35    0   |    JAC = |   3   2   5    *  |
c       |   44  43   46    0   |          |   4   3   6    *  |
c       |   55  51   54    0   |          |   5   1   4    *  |
c       (   66  61   62   65   )          (   6   1   2    5  )
c
c   Note: * means that you can store values from 1 to 6 (1 to n, where
c         n is the order of the matrix) in that position in the array.
c
c   Contributed by:
c
c   Ernest E. Rothman
c   Cornell Thoery Center/Cornell National Supercomputer Facility
c   e-mail address: BITNET:   EER@CORNELLF.BITNET
c                   INTERNET: eer@cornellf.tn.cornell.edu
c
c   checked and modified  04/13/90 Y.Saad.
c
c   REFERENCES
c
c   Kincaid, D. R.; Oppe, T. C.; Respess, J. R.; Young, D. M. 1984.
c   ITPACKV 2C User's Guide, CNA-191. Center for Numerical Analysis,
c   University of Texas at Austin.
c
c   "Engineering and Scientific Subroutine Library; Guide and
c   Reference; Release 3 (SC23-0184-3). Pp. 79-86.
c
c   INPUT PARAMETERS
c
c  N       - Integer. The size of the square matrix.
c
c  NNZ     - Integer. Must be greater than or equal to the number of
c            nonzero elements in the sparse matrix. Dimension of A, IA
c            and JA.
c
c  NCA     - Integer. First dimension of output arrays ca and jac.
c
c  A(NNZ)  - Real array.
c            Stored entries of the sparse matrix A.
c            NNZ is the number of nonzeros.
c
c  IA(NNZ) - Integer array.
c            Pointers to specify rows for the stored nonzero entries
c            in A.
c
c  JA(NNZ) - Integer array.
c            Pointers to specify columns for the stored nonzero
c            entries in A.
c
c  NER     - Integer. Must be set greater than or equal to the maximum
c            number of nonzeros in any row of the sparse matrix.
c
c  OUTPUT PARAMETERS
c
c  AC(NAC,*)  - Real array.
c               Stored entries of the sparse matrix A in compressed
c               storage mode.
c
c  JAC(NAC,*) - Integer array.
c               Contains the column numbers of the sparse matrix
c               elements stored in the corresponding positions in
c               array AC.
c
c  NCMAX   -  Integer. Equals the maximum number of nonzeros in any
c             row of the sparse matrix.
c
c  IERR    - Error parameter is returned as zero on successful
c             execution of the subroutin<e.
c             Error diagnostics are given by means of positive values
c             of this parameter as follows:
c
c             IERR = -1   -  NER is too small and should be set equal
c                            to NCMAX. The array AC may not be large
c                            enough to accomodate all the non-zeros of
c                            of the sparse matrix.
c             IERR =  1   -  The array AC has a zero column. (Warning)
c             IERR =  2   -  The array AC has a zero row.    (Warning)
c
c
      double precision a(nnz), ac(nac,ner)
      integer ja(nnz), ia(nnz), jac(nac,ner), ierr, ncmax, icount
c
c   Initial error parameter to zero:
c
      ierr = 0
c
c   Initial output arrays to zero:
c
      do 4 in = 1,ner
         do 4 innz =1,n
            jac(innz,in) = n
            ac(innz,in) = 0.0
 4    continue
c
c   Assign nonzero elements of the sparse matrix (stored in the one
c   dimensional array A to the two dimensional array AC.
c   Also, assign the correct values with information about their
c   column indices to the two dimensional array KA. And at the same
c   time count the number of nonzeros in each row so that the
c   parameter NCMAX equals the maximum number of nonzeros in any row
c   of the sparse matrix.
c
      ncmax = 1
      do 10 is = 1,n
         k = 0
         do 30 ii = 1,nnz
            if(ia(ii).eq.is)then
               k = k + 1
               if (k .le. ner) then
                  ac(is,k) = a(ii)
                  jac(is,k) = ja(ii)
               end if
            end if
 30      continue
         if (k.ge.ncmax) ncmax = k
 10   continue
c
c     Perform some simple error checks:
c
c  check maximum number of nonzeros in each row:
c
      if (ncmax.eq.ner) ierr = 0
      if (ncmax.gt.ner) then
         ierr = -1
         return
      end if
c
c  check if there are any zero columns in AC:
c
      do 45 in = 1,ncmax
         icount = 0
         do 44 inn =1,n
            if (ac(inn,in).ne.0.0) icount = 1
 44      continue
         if (icount.eq.0) then
            ierr = 1
            return
         end if
 45   continue
c
c  check if there are any zero rows in AC:
c
      do 55 inn = 1,n
         icount = 0
         do 54 in =1,ncmax
            if (ac(inn,in).ne.0.0) icount = 1
 54      continue
         if (icount.eq.0) then
            ierr = 2
            return
         end if
 55   continue
      return
      end
      subroutine copmat (nrow,a,ja,ia,ao,jao,iao,ipos)

c*********************************************************************72
c
cc COPMAT copies the matrix A, JA, IA, into the matrix AO, JAO, IAO.
c
c on entry:
c
c nrow      = row dimension of the matrix
c a,
c ja,
c ia    = input matrix in compressed sparse row format.
c ipos  = integer. indicates the position in the array ao, jao
c         where the first element should be copied. Thus
c         iao(1) = ipos on return.
c
c on return:
c
c ao,
c jao,
c iao   = output matrix containing the same data as a, ja, ia.
c
c           Y. Saad, March 1990.
c
      double precision a(*),ao(*)
      integer nrow, ia(*),ja(*),jao(*),iao(*), ipos
      integer kst, i, k

      kst    = ipos -ia(1)
      do 100 i = 1, nrow+1
         iao(i) = ia(i) + kst
 100  continue

      do 200 k=ia(1), ia(nrow+1)-1
         ao(kst+k) = a(k)
         jao(kst+k)= ja(k)
 200  continue

      return
      end
      subroutine cperm (nrow,a,ja,ia,ao,jao,iao,perm,job)

c*********************************************************************72
c
cc CPERM permutes the columns of a matrix.
c
c  This routine permutes the columns of a matrix a, ja, ia.
c the result is written in the output matrix  ao, jao, iao.
c cperm computes B = A P, where  P is a permutation matrix
c that maps column j into column perm(j), i.e., on return
c      a(i,j) becomes a(i,perm(j)) in new matrix
c Y. Saad, May 2, 1990 / modified Jan. 28, 1991.
c
c on entry:
c
c nrow       = row dimension of the matrix
c
c a, ja, ia = input matrix in csr format.
c
c perm      = integer array of length ncol (number of columns of A
c         containing the permutation array  the columns:
c         a(i,j) in the original matrix becomes a(i,perm(j))
c         in the output matrix.
c
c job      = integer indicating the work to be done:
c             job = 1      permute a, ja, ia into ao, jao, iao
c                       (including the copying of real values ao and
c                       the array iao).
c             job .ne. 1 :  ignore real values ao and ignore iao.
c
c
c on return:
c
c ao, jao, iao = input matrix in a, ja, ia format (array ao not needed)
c
c Notes:
c
c 1. if job=1 then ao, iao are not used.
c 2. This routine is in place: ja, jao can be the same.
c 3. If the matrix is initially sorted (by increasing column number)
c    then ao,jao,iao  may not be on return.
c
c
      integer nrow,ja(*),ia(nrow+1),jao(*),iao(nrow+1),perm(*), job
      double precision a(*), ao(*)
      integer k, i, nnz

      nnz = ia(nrow+1)-1
      do 100 k=1,nnz
         jao(k) = perm(ja(k))
 100  continue
c
c     done with ja array. return if no need to touch values.
c
      if (job .ne. 1) return
c
c else get new pointers -- and copy values too.
c
      do 1 i=1, nrow+1
         iao(i) = ia(i)
 1    continue

      do 2 k=1, nnz
         ao(k) = a(k)
 2    continue

      return
      end
      subroutine cscal(nrow, job, nrm, a, ja, ia, diag, b, jb, ib)

c*********************************************************************72
c
cc CSCAL scales the columns of A such that their norms are one.
c
c result matrix written on b, or overwritten on A.
c 3 choices of norms: 1-norm, 2-norm, max-norm. in place.
c
c on entry:
c
c nrow      = integer. The row dimension of A
c
c job   = integer. job indicator. Job=0 means get array b only
c         job = 1 means get b, and the integer arrays ib, jb.
c
c nrm   = integer. norm indicator. nrm = 1, means 1-norm, nrm =2
c                  means the 2-nrm, nrm = 0 means max norm
c
c a,
c ja,
c ia   = Matrix A in compressed sparse row format.
c
c on return:
c
c
c diag = diagonal matrix stored as a vector containing the matrix
c        by which the columns have been scaled, i.e., on return
c        we have B = A * Diag
c
c b,
c jb,
c ib      = resulting matrix B in compressed sparse row sparse format.
c
c Notes:
c
c 1)        The column dimension of A is not needed.
c 2)       algorithm in place (B can take the place of A).
c
      double precision a(*), b(*), diag(nrow)
      integer ja(*),jb(*), ia(nrow+1),ib(nrow+1)
      call cnrms (nrow,nrm,a,ja,ia,diag)
      do 1 j=1, nrow
         diag(j) = 1.0/diag(j)
 1    continue
      call amudia (nrow,job,a,ja,ia,diag,b,jb,ib)
      return
      end
      subroutine csort (n,a,ja,ia,iwork,values)

c*********************************************************************72
c
cc CSORT sorts the elements of a CSR matrix.
c
c This routine sorts the elements of  a matrix (stored in Compressed
c Sparse Row Format) in increasing order of their column indices within
c each row. It uses a form of bucket sort with a cost of O(nnz) where
c nnz = number of nonzero elements.
c requires an integer work array of size length 2*nnz.
c
c on entry:
c
c n     = the row dimension of the matrix
c a     = the matrix A in compressed sparse row format.
c ja    = the array of column indices of the elements in array a.
c ia    = the array of pointers to the rows.
c iwork = integer work array of length max ( n+1, 2*nnz )
c         where nnz = 2* (ia(n+1)-ia(1))  ) .
c values= logical indicating whether or not the real values a(*) must
c         also be permuted. if (.not. values) then the array a is not
c         touched by csort and can be a dummy array.
c
c on return:
c
c the matrix stored in the structure a, ja, ia is permuted in such a
c way that the column indices are in increasing order within each row.
c iwork(1:nnz) contains the permutation used  to rearrange the elements.
c
c Y. Saad - Feb. 1, 1991.
c
      logical values
      integer n, ja(*), ia(n+1), iwork(*)
      double precision a(*)
      integer i, k, j, ifirst, nnz, next
c
c count the number of elements in each column
c
      do 1 i=1,n+1
         iwork(i) = 0
 1    continue
      do 3 i=1, n
         do 2 k=ia(i), ia(i+1)-1
            j = ja(k)+1
            iwork(j) = iwork(j)+1
 2       continue
 3    continue
c
c compute pointers from lengths.
c
      iwork(1) = 1
      do 4 i=1,n
         iwork(i+1) = iwork(i) + iwork(i+1)
 4    continue
c
c get the positions of the nonzero elements in order of columns.
c
      ifirst = ia(1)
      nnz = ia(n+1)-ifirst
      do 5 i=1,n
         do 51 k=ia(i),ia(i+1)-1
            j = ja(k)
            next = iwork(j)
            iwork(nnz+next) = k
            iwork(j) = next+1
 51      continue
 5    continue
c
c convert to coordinate format
c
      do 6 i=1, n
         do 61 k=ia(i), ia(i+1)-1
            iwork(k) = i
 61      continue
 6    continue
c
c loop to find permutation: for each element find the correct
c position in (sorted) arrays a, ja. Record this in iwork.
c
      do 7 k=1, nnz
         ko = iwork(nnz+k)
         irow = iwork(ko)
         next = ia(irow)
c
c the current element should go in next position in row. iwork
c records this position.
c
         iwork(ko) = next
         ia(irow)  = next+1
 7       continue
c
c perform an in-place permutation of the  arrays.
c
         call ivperm (nnz, ja(ifirst), iwork)
         if (values) call dvperm (nnz, a(ifirst), iwork)
c
c reshift the pointers of the original matrix back.
c
      do 8 i=n,1,-1
         ia(i+1) = ia(i)
 8    continue
      ia(1) = ifirst
c
      return
      end
      subroutine csrbnd (n,a,ja,ia,job,abd,nabd,lowd,ml,mu,ierr)

c*********************************************************************72
c
cc CSRBND converts Compressed Sparse Row to Banded Linpack format.
c
c CSRBND converts a general sparse matrix stored in
c compressed sparse row format into the banded format. for the
c banded format,the Linpack conventions are assumed (see below).
c
c on entry:
c
c n      = integer,the actual row dimension of the matrix.
c
c a,
c ja,
c ia    = input matrix stored in compressed sparse row format.
c
c job      = integer. if job=1 then the values of the lower bandwith ml
c         and the upper bandwidth mu are determined internally.
c         otherwise it is assumed that the values of ml and mu
c         are the correct bandwidths on input. See ml and mu below.
c
c nabd  = integer. first dimension of array abd.
c
c lowd  = integer. this should be set to the row number in abd where
c         the lowest diagonal (leftmost) of A is located.
c         lowd should be  ( 1  .le.  lowd  .le. nabd).
c         if it is not known in advance what lowd should be
c         enter lowd = 0 and the default value lowd = ml+mu+1
c         will be chosen. Alternative: call routine getbwd from unary
c         first to detrermione ml and mu then define lowd accordingly.
c         (Note: the banded solvers in linpack use lowd=2*ml+mu+1. )
c
c ml      = integer. equal to the bandwidth of the strict lower part of A
c mu      = integer. equal to the bandwidth of the strict upper part of A
c         thus the total bandwidth of A is ml+mu+1.
c         if ml+mu+1 is found to be larger than lowd then an error
c         flag is raised (unless lowd = 0). see ierr.
c
c note:   ml and mu are assumed to have       the correct bandwidth values
c         as defined above if job is set to zero on entry.
c
c on return:
c
c
c abd   = double precision array of dimension abd(nabd,n).
c         on return contains the values of the matrix stored in
c         banded form. The j-th column of abd contains the elements
c         of the j-th column of  the original matrix comprised in the
c         band ( i in (j-ml,j+mu) ) with the lowest diagonal at
c         the bottom row (row lowd). See details below for this format.
c
c ml      = integer. equal to the bandwidth of the strict lower part of A
c mu      = integer. equal to the bandwidth of the strict upper part of A
c         if job=1 on entry then these two values are internally computed.
c
c lowd  = integer. row number in abd where the lowest diagonal
c         (leftmost) of A is located on return. In case lowd = 0
c         on return, then it is defined to ml+mu+1 on return and the
c         lowd will contain this value on return. `
c
c ierr  = integer. used for error messages. On return:
c         ierr .eq. 0  :means normal return
c         ierr .eq. -1 : means invalid value for lowd. (either .lt. 0
c         or larger than nabd).
c         ierr .eq. -2 : means that lowd is not large enough and as
c         result the matrix cannot be stored in array abd.
c         lowd should be at least ml+mu+1, where ml and mu are as
c         provided on output.
c
c Additional details on banded format.  (this closely follows the
c format used in linpack. may be useful for converting a matrix into
c this storage format in order to use the linpack  banded solvers).
c
c band storage format  for matrix abd
c uses ml+mu+1 rows of abd(nabd,*) to store the diagonals of
c a in rows of abd starting from the lowest (sub)-diagonal  which  is
c stored in row number lowd of abd. the minimum number of rows needed
c in abd is ml+mu+1, i.e., the minimum value for lowd is ml+mu+1. the
c j-th  column  of  abd contains the elements of the j-th column of a,
c from bottom to top: the element a(j+ml,j) is stored in  position
c abd(lowd,j), then a(j+ml-1,j) in position abd(lowd-1,j) and so on.
c Generally, the element a(j+k,j) of original matrix a is stored in
c position abd(lowd+k-ml,j), for k=ml,ml-1,..,0,-1, -mu.
c The first dimension nabd of abd must be .ge. lowd
c
c     example [from linpack ]:   if the original matrix is
c
c              11 12 13  0  0  0
c              21 22 23 24  0  0
c               0 32 33 34 35  0     original banded matrix
c               0  0 43 44 45 46
c               0  0  0 54 55 56
c               0  0  0  0 65 66
c
c then  n = 6, ml = 1, mu = 2. lowd should be .ge. 4 (=ml+mu+1)  and
c if lowd = 5 for example, abd  should be:
c
c untouched --> x  x  x  x  x  x
c               *  * 13 24 35 46
c               * 12 23 34 45 56    resulting abd matrix in banded
c              11 22 33 44 55 66    format
c  row lowd--> 21 32 43 54 65  *
c
c * = not used
c
      double precision a(*),abd(nabd,n)
      integer ia(n+1),ja(*)
c
c first determine ml and mu.
c
      ierr = 0

      if (job .eq. 1) call getbwd(n,a,ja,ia,ml,mu)
      m = ml+mu+1
      if (lowd .eq. 0) lowd = m
      if (m .gt. lowd)  ierr = -2
      if (lowd .gt. nabd .or. lowd .lt. 0) ierr = -1
      if (ierr .lt. 0) return

      do 15  i=1,m
         ii = lowd -i+1
         do 10 j=1,n
          abd(ii,j) = 0.0
 10      continue
 15   continue

      mdiag = lowd-ml
      do 30 i=1,n
         do 20 k=ia(i),ia(i+1)-1
            j = ja(k)
            abd(i-j+mdiag,j) = a(k)
 20      continue
 30   continue
      return
      end
      subroutine csrbsr (n,nblk,na,a,ja,ia,ao,jao,iao)

c*********************************************************************72
c
cc CSRBSR converts Compressed Sparse Row to Block Sparse Row.
c
c CSRBSR does the reverse of bsrcsr. It converts
c a matrix stored in a general compressed a, ja, ia format into a
c a block reduced matrix a(*,*),ja(*),ia(*) format. The code
c assumes that the original matrix is indeed a block format
c and that the elements are ordered in such a way that their
c column numbers are increasing. (This can be achieved
c by transposing a, ja, ia twice, putting the resulting matrix
c into a, ja, ia).
c
c See routine bsrcsr for more details on data structure for blocked
c matrices. The input matrix is a, ja, ia (in compressed format) and
c the output matrix is the matrix ao, jao, iao in block-reduced
c format.
c
c on entry:
c
c n      = integer, the actual row dimension of the matrix.
c nblk  = integer equal to the dimension of each block.
c         nblk must divide n.
c na      = first dimension of array a as declared in calling program
c
c a, ja,
c    ia = input matrix stored in compressed sparse row format.
c
c on return:
c
c ao    = double precision array containing the values of the matrix. For details
c         on the format see below. Each row of a contains the nblk x nblk
c         block matrix unpacked column-wise (this allows the user to
c         declare the array a as a(na,nblk,nblk) on entry if desired).
c         the block rows are stored in sequence just as for the compressed
c         sparse row format.
c jao   = integer array of length n/nblk. ja(k) contains the column index
c         of the leading element, i.e., the element (1,1) of the block
c         that is held in the row a(k,*) of the value array.
c iao   = integer array of length n/nblk+1. ia(i) points to the beginning
c        of block row number i in the arrays a and ja.
c
c Notes:
c  1) this code is not in place.
c        2) see routine bsrcsr for details on data sctructure for
c           block sparse row format.
c        3) The routine assumes that  the input matrix has been
c           sorted in such a way that the column indices are always
c           in increasing order for the same row.
c           for all k "in the SAME ROW."
c        4) THERE IS NO CHECKING AS TO WHETHER the input is correct.
c           it is recommended to use the routine blchk to check
c           if the matrix is a block-matrix before calling csrbsr.
c
      double precision a(*),ao(na,*)
      integer ia(n+1),ja(*),jao(*),iao(*)
c
c     nr is the dimension of the reduced matrix.
c
      nr = n/nblk
      iao(1) = 1
      ibrow = 1
      irow  = 1
c  main loop
      do 2 ii=1, nr
c     i1= starting position for group of nblk rows in original matrix
         i1 = ia(irow)
c     lena = length of each row in that group  in the original matrix
         lena = ia(irow+1)-i1
c     len = length of each block-row in that group in the output matrix
         len = lena/nblk
         k1 = iao(ibrow)
c  copy the real values of A
c     for each block
         do 7 k=0, len-1
c     store column positions of the (1,1) elements of each block
            jao(k1+k) = ja(i1+nblk*k)
c     for each column
            do 5 j=1, nblk
               j1 = (j-1)*nblk
               j2 = i1+k*nblk+j-1
c     for each row
               do 6 i = 1, nblk
                  ao(k1+k,j1+i) = a(j2+(i-1)*lena)
 6             continue
 5          continue
 7       continue
c     done with a whole block row. now update iao(*), ibrow and irow
         iao(ibrow+1) = iao(ibrow)+len
         ibrow = ibrow + 1
         irow  = irow + nblk
 2    continue
      return
      end
      subroutine csrcoo (nrow,job,nzmax,a,ja,ia,nnz,ao,ir,jc,ierr)

c*********************************************************************72
c
cc CSRCOO converts Compressed Sparse Row to Coordinate format.
c
c converts a matrix that is stored in coordinate format
c  a, ir, jc into a row general sparse ao, jao, iao format.
c
c on entry:
c
c nrow      = dimension of the matrix.
c job   = integer serving as a job indicator.
c         if job = 1 fill in only the array ir, ignore jc, and ao.
c         if job = 2 fill in ir, and jc but not ao
c         if job = 3 fill in everything.
c         The reason why these options are provided is that on return
c         ao and jc are the same as a, ja. So when job = 3, a and ja are
c         simply copied into ao, jc.  When job=2, only jc and ir are
c         returned. With job=1 only the array ir is returned. Moreover,
c         the algorithm is in place:
c           call csrcoo (nrow,1,nzmax,a,ja,ia,nnz,a,ia,ja,ierr)
c         will write the output matrix in coordinate format on a, ja,ia.
c         (Important: note the order in the output arrays a, ja, ia. )
c         i.e., ao can be the same as a, ir can be the same as ia
c         and jc can be the same as ja.
c
c a,
c ja,
c ia    = matrix in compressed sparse row format.
c nzmax = length of space available in ao, ir, jc.
c         the code will stop immediatly if the number of
c         nonzero elements found in input matrix exceeds nzmax.
c
c on return:
c
c ao, ir, jc = matrix in coordinate format.
c
c nnz        = number of nonzero elements in matrix.
c ierr       = integer error indicator.
c         ierr .eq. 0 means normal retur
c         ierr .eq. 1 means that the the code stopped
c         because there was no space in ao, ir, jc
c         (according to the value of  nzmax).
c
      double precision a(*),ao(*)
      integer ir(*),jc(*),ja(*),ia(*)

      ierr = 0
      nnz = ia(nrow+1)-1
      if (nnz .gt. nzmax) then
         ierr = 1
         return
      end if

      goto (3,2,1) job
 1    do 10 k=1,nnz
         ao(k) = a(k)
 10   continue
 2    do 11 k=1,nnz
         jc(k) = ja(k)
 11   continue
c copy backward to allow
 3    do 13 i=nrow,1,-1
         k1 = ia(i+1)-1
         k2 = ia(i)
         do 12 k=k1,k2,-1
            ir(k) = i
 12      continue
 13   continue
      return
      end
      subroutine csrcsc(n,job,ipos,a,ja,ia,ao,jao,iao)

c*********************************************************************72
c
cc CSRCSC converts Compressed Sparse Row to Compressed Sparse Column.
c
c (transposition operation)   Not in place.
c
c on entry:
c
c n      = dimension of A.
c job      = integer to indicate whether or not to fill the values of the
c         matrix ao or only the pattern (ia, and ja). Enter 1 for yes.
c ipos  = starting position in ao, jao of the transposed matrix.
c         the iao array takes this into account (thus iao(1) is set to ipos.)
c         Note: this may be useful if one needs to append the data structure
c         of the transpose to that of A. In this case use
c                call csrcsc (n,1,n+2,a,ja,ia,a,ja,ia(n+2))
c        for any other normal usage, enter ipos=1.
c a      = double precision array of length nnz (nnz=number of nonzero elements in input
c         matrix) containing the nonzero elements.
c ja      = integer array of length nnz containing the column positions
c         of the corresponding elements in a.
c ia      = integer of size n+1. ia(k) contains the position in a, ja of
c        the beginning of the k-th row.
c
c on return:
c
c output arguments:
c ao      = double precision array of size nzz containing the "a" part of the transpose
c jao      = integer array of size nnz containing the column indices.
c iao      = integer array of size n+1 containing the "ia" index array of
c        the transpose.
c
      integer ia(n+1),iao(n+1),ja(*),jao(*)
      double precision  a(*),ao(*)
c
c  compute lengths of rows of transp(A)
      do 1 i=1,n+1
         iao(i) = 0
 1    continue
      do 3 i=1, n
         do 2 k=ia(i), ia(i+1)-1
            j = ja(k)+1
            iao(j) = iao(j)+1
 2       continue
 3    continue
c compute pointers from lengths
      iao(1) = ipos
      do 4 i=1,n
         iao(i+1) = iao(i) + iao(i+1)
 4    continue
c  now do the actual copying
      do 6 i=1,n
         do 62 k=ia(i),ia(i+1)-1
            j = ja(k)
            next = iao(j)
            if (job .eq. 1)  ao(next) = a(k)
            jao(next) = i
            iao(j) = next+1
 62      continue
 6    continue
c  reshift iao and leave
      do 7 i=n,1,-1
         iao(i+1) = iao(i)
 7    continue
      iao(1) = ipos
      end
      subroutine csrdia (n,idiag,job,a,ja,ia,ndiag,
     *                   diag,ioff,ao,jao,iao,ind)

c*********************************************************************72
c
cc CSRDIA converts Compressed Sparse Row to diagonal format.
c
c CSRDIA extracts  idiag diagonals  from the  input matrix a,
c a, ia, and puts the rest of  the matrix  in the  output matrix ao,
c jao, iao.  The diagonals to be extracted depend  on the  value of job
c (see below for details.)  In  the first  case, the  diagonals to be
c extracted are simply identified by  their offsets  provided in ioff
c by the caller.  In the second case, the  code internally determines
c the idiag most significant diagonals, i.e., those  diagonals of the
c matrix which  have  the  largest  number  of  nonzero elements, and
c extracts them.
c
c on entry:
c
c n      = dimension of the matrix a.
c idiag = integer equal to the number of diagonals to be extracted.
c         Note: on return idiag may be modified.
c a, ja,
c    ia = matrix stored in a, ja, ia, format
c job      = integer. serves as a job indicator.  Job is better thought
c         of as a two-digit number job=xy. If the first (x) digit
c         is one on entry then the diagonals to be extracted are
c         internally determined. In this case csrdia exctracts the
c         idiag most important diagonals, i.e. those having the largest
c         number on nonzero elements. If the first digit is zero
c         then csrdia assumes that ioff(*) contains the offsets
c         of the diagonals to be extracted. there is no verification
c         that ioff(*) contains valid entries.
c         The second (y) digit of job determines whether or not
c         the remainder of the matrix is to be written on ao,jao,iao.
c         If it is zero  then ao, jao, iao is not filled, i.e.,
c         the diagonals are found  and put in array diag and the rest is
c         is discarded. if it is one, ao, jao, iao contains matrix
c         of the remaining elements.
c         Thus:
c         job= 0 means do not select diagonals internally (pick those
c                defined by ioff) and do not fill ao,jao,iao
c         job= 1 means do not select diagonals internally
c                      and fill ao,jao,iao
c         job=10 means  select diagonals internally
c                      and do not fill ao,jao,iao
c         job=11 means select diagonals internally
c                      and fill ao,jao,iao
c
c ndiag = integer equal to the first dimension of array diag.
c
c on return:
c
c
c idiag = number of diagonals found. This may be smaller than its value
c         on entry.
c diag  = double precision array of size (ndiag x idiag) containing the diagonals
c         of A on return
c
c ioff  = integer array of length idiag, containing the offsets of the
c           diagonals to be extracted.
c ao, jao
c  iao  = remainder of the matrix in a, ja, ia format.
c work arrays:
c
c ind   = integer array of length 2*n-1 used as integer work space.
c         needed only when job.ge.10 i.e., in case the diagonals are to
c         be selected internally.
c
c Notes:
c
c    1) The algorithm is in place: ao, jao, iao can be overwritten on
c       a, ja, ia if desired
c    2) When the code is required to select the diagonals (job .ge. 10)
c       the selection of the diagonals is done from left to right
c       as a result if several diagonals have the same weight (number
c       of nonzero elemnts) the leftmost one is selected first.
c
      double precision diag(ndiag,idiag), a(*), ao(*)
      integer ia(*), ind(*), ja(*), jao(*), iao(*), ioff(*)

      job1 = job/10
      job2 = job-job1*10
      if (job1 .eq. 0) goto 50
      n2 = n+n-1
      call infdia(n,ja,ia,ind,idum)
c  determine diagonals to  accept.
c
      ii = 0
 4    ii=ii+1
      jmax = 0
      do 41 k=1, n2
         j = ind(k)
         if (j .le. jmax) goto 41
         i = k
         jmax = j
 41   continue
      if (jmax .le. 0) then
         ii = ii-1
         goto 42
      end if
      ioff(ii) = i-n
      ind(i) = - jmax
      if (ii .lt.  idiag) goto 4
 42   idiag = ii
c  initialize diago to zero
 50   continue
      do 55 j=1,idiag
         do 54 i=1,n
            diag(i,j) = 0.0
 54      continue
 55   continue

      ko = 1
c
c extract diagonals and accumulate remaining matrix.
c
      do 6 i=1, n
         do 51 k=ia(i),ia(i+1)-1
            j = ja(k)
            do 52 l=1,idiag
               if (j-i .ne. ioff(l)) goto 52
               diag(i,l) = a(k)
               goto 51
 52         continue
c  append element not in any diagonal to ao,jao,iao
            if (job2 .eq. 0) goto 51
            ao(ko) = a(k)
            jao(ko) = j
            ko = ko+1
 51      continue
         if (job2 .ne. 0 ) ind(i+1) = ko
 6    continue
      if (job2 .eq. 0) return
c     finish with iao
      iao(1) = 1
      do 7 i=2,n+1
         iao(i) = ind(i)
 7    continue
      return
      end
      subroutine csrdns(nrow,ncol,a,ja,ia,dns,ndns,ierr)

c*********************************************************************72
c
cc CSRDNS converts Compressed Sparse Row to Dense format.
c
c converts a row-stored sparse matrix into a densely stored one
c
c On entry:
c
c
c nrow      = row-dimension of a
c ncol      = column dimension of a
c a,
c ja,
c ia    = input matrix in compressed sparse row format.
c         (a=value array, ja=column array, ia=pointer array)
c dns   = array where to store dense matrix
c ndns      = first dimension of array dns
c
c on return:
c
c dns   = the sparse matrix a, ja, ia has been stored in dns(ndns,*)
c
c ierr  = integer error indicator.
c         ierr .eq. 0  means normal return
c         ierr .eq. i  means that the code has stopped when processing
c         row number i, because it found a column number .gt. ncol.
c
      double precision dns(ndns,*),a(*)
      integer ja(*),ia(*)

      ierr = 0
      do 1 i=1, nrow
         do 2 j=1,ncol
          dns(i,j) = 0.0
 2       continue
 1    continue

      do 4 i=1,nrow
         do 3 k=ia(i),ia(i+1)-1
            j = ja(k)
          if (j .gt. ncol) then
               ierr = i
               return
          end if
          dns(i,j) = a(k)
 3       continue
 4    continue
      return
      end
      subroutine csrell (nrow,a,ja,ia,maxcol,coef,jcoef,ncoef,
     *                   ndiag,ierr)

c*********************************************************************72
c
cc CSRELL converts Compressed Sparse Row to Ellpack/Itpack format.
c
c on entry:
c
c nrow         = row dimension of the matrix A.
c
c a,
c ia,
c ja      = input matrix in compressed sparse row format.
c
c ncoef  = first dimension of arrays coef, and jcoef.
c
c maxcol = integer equal to the number of columns available in coef.
c
c on return:
c
c coef      = double precision array containing the values of the matrix A in
c         itpack-ellpack format.
c jcoef = integer array containing the column indices of coef(i,j)
c         in A.
c ndiag = number of active 'diagonals' found.
c
c ierr       = error message. 0 = correct return. If ierr .ne. 0 on
c        return this means that the number of diagonals found
c         (ndiag) exceeds maxcol.
c
      integer ia(nrow+1), ja(*), jcoef(ncoef,1)
      double precision a(*), coef(ncoef,1)
c
c first determine the length of each row of lower-part-of(A)
      ierr = 0
      ndiag = 0
      do 3 i=1, nrow
         k = ia(i+1)-ia(i)
         ndiag = max0(ndiag,k)
 3    continue
c  check whether sufficient columns are available.
      if (ndiag .gt. maxcol) then
         ierr = 1
         return
      end if
c
c fill coef with zero elements and jcoef with row numbers.
c
      do 4 j=1,ndiag
         do 41 i=1,nrow
            coef(i,j) = 0.0
            jcoef(i,j) = i
 41      continue
 4    continue
c
c  copy elements row by row.
c
      do 6 i=1, nrow
         k1 = ia(i)
         k2 = ia(i+1)-1
         do 5 k=k1,k2
            coef(i,k-k1+1) = a(k)
            jcoef(i,k-k1+1) = ja(k)
 5       continue
 6    continue
      return
      end
      subroutine csrjad (nrow, a, ja, ia, idiag, iperm, ao, jao, iao)

c*********************************************************************72
c
cc CSRJAD converts Compressed Sparse Row to Jagged Diagonal storage.
c
c CSRJAG converts  matrix stored in the compressed sparse
c row format to the jagged diagonal format. The data structure
c for the JAD (Jagged Diagonal storage) is as follows. The rows of
c the matrix are (implicitly) permuted so that their lengths are in
c decreasing order. The real entries ao(*) and their column indices
c jao(*) are stored in succession. The number of such diagonals is idiag.
c the lengths of each of these diagonals is stored in iao(*).
c For more details see [E. Anderson and Y. Saad,
c ``Solving sparse triangular systems on parallel computers'' in
c Inter. J. of High Speed Computing, Vol 1, pp. 73-96 (1989).]
c or  [Y. Saad, ``Krylov Subspace Methods on Supercomputers''
c SIAM J. on  Stat. Scient. Comput., volume 10, pp. 1200-1232 (1989).]
c
c on entry:
c
c nrow         = row dimension of the matrix A.
c
c a,
c ia,
c ja      = input matrix in compressed sparse row format.
c
c on return:
c
c
c idiag = integer. The number of jagged diagonals in the matrix.
c
c iperm = integer array of length nrow containing the permutation
c         of the rows that leads to a decreasing order of the
c         number of nonzero elements.
c
c ao    = double precision array containing the values of the matrix A in
c         jagged diagonal storage. The j-diagonals are stored
c         in ao in sequence.
c
c jao   = integer array containing the column indices of the
c         entries in ao.
c
c iao   = integer array containing pointers to the beginning
c         of each j-diagonal in ao, jao. iao is also used as
c         a work array and it should be of length n at least.
c
      integer ja(*), jao(*), ia(nrow+1), iperm(nrow), iao(nrow)
      double precision a(*), ao(*)
c
c  define initial iperm and get lengths of each row
c  jao is used a work vector to store tehse lengths
c
      idiag = 0
      ilo = nrow
      do 10 j=1, nrow
         iperm(j) = j
         len = ia(j+1) - ia(j)
         ilo = min(ilo,len)
         idiag = max(idiag,len)
         jao(j) = len
 10   continue
c
c     call sorter to get permutation. use iao as work array.
c
      call dcsort (jao, nrow, iao, iperm, ilo, idiag)
c
c     define output data structure. first lengths of j-diagonals
c
      do 20 j=1, nrow
         iao(j) = 0
 20   continue
      do 40 k=1, nrow
         len = jao(iperm(k))
         do 30 i=1,len
            iao(i) = iao(i)+1
 30      continue
 40   continue
c
c     get the output matrix itself
c
      k1 = 1
      k0 = k1
      do 60 jj=1, idiag
         len = iao(jj)
         do 50 k=1,len
            i = ia(iperm(k))+jj-1
            ao(k1) = a(i)
            jao(k1) = ja(i)
            k1 = k1+1
 50      continue
         iao(jj) = k0
         k0 = k1
 60   continue
      iao(idiag+1) = k1
      return
      end
      subroutine csrlnk (n,a,ja,ia,link)

c*********************************************************************72
c
cc CSRLNK converts Compressed Sparse Row to Linked storage format.
c
c CSRLNK translates a matrix stored in compressed sparse
c row into one with a linked list storage format. Only the link
c array needs to be obtained since the arrays a, ja, and ia may
c be unchanged and have carry the same meaning for the output matrix.
c in  other words a, ja, ia, link   ia the output linked list data
c structure with a, ja, ia being the same.
c
c Coded by Y. Saad, Feb 21, 1991.
c
c
c on entry:
c
c n      = integer equal to the dimension of A.
c
c a      = double precision array of size nna containing the nonzero elements
c ja      = integer array of size      nnz containing the column positions
c         of the corresponding elements in a.
c ia      = integer of size n+1 containing the pointers to the beginning
c         of each row. ia(k) contains the position in a, ja of the
c         beginning of the k-th row.
c
c on return:
c
c a, ja, ia are not changed.
c
c a     = nonzero elements.
c ja    = column positions.
c ia    = points to the first row of matrix in structure.
c link      = integer array of size containing the linked list information.
c         link(k) points to the next element of the row after element
c         ao(k), jcol(k). if link(k) = 0, then there is no next element,
c         i.e., ao(k), jcol(k) is the last element of the current row.
c
      double precision a(*)
      integer n, ja(*), ia(n+1), link(*)
      integer i, k
c
c loop through all rows
c
      do 100 i =1, n
         do 99  k=ia(i), ia(i+1)-2
            link(k) = k+1
 99      continue
         link(ia(i+1)-1) = 0
 100  continue

      return
      end
      subroutine csrmsr (n,a,ja,ia,ao,jao,wk,iwk)

c*********************************************************************72
c
cc CSRMSR converts Compressed Sparse Row to Modified Sparse Row.
c
c converts a general sparse matrix a, ja, ia into
c a compressed matrix using a separated diagonal (referred to as
c the bell-labs format as it is used by bell labs semi conductor
c group. We refer to it here as the modified sparse row format.
c Note: this has been coded in such a way that one can overwrite
c the output matrix onto the input matrix if desired by a call of
c the form
c
c     call csrmsr (n, a, ja, ia, a, ja, wk,iwk)
c
c In case ao, jao, are different from a, ja, then one can
c use ao, jao as the work arrays in the calling sequence:
c
c     call csrmsr (n, a, ja, ia, ao, jao, ao,jao)
c
c
c
c on entry :
c
c a, ja, ia = matrix in csr format. note that the
c           algorithm is in place: ao, jao can be the same
c            as a, ja, in which case it will be overwritten on it
c            upon return.
c
c on return :
c
c ao, jao  = sparse matrix in modified sparse row storage format:
c         +  ao(1:n) contains the diagonal of the matrix.
c         +  ao(n+2:nnz) contains the nondiagonal elements of the
c             matrix, stored rowwise.
c         +  jao(n+2:nnz) : their column indices
c         +  jao(1:n+1) contains the pointer array for the nondiagonal
c             elements in ao(n+1:nnz) and jao(n+2:nnz).
c             i.e., for i .le. n+1 jao(i) points to beginning of row i
c            in arrays ao, jao.
c             here nnz = number of nonzero elements+1
c work arrays:
c
c wk      = real work array of length n
c iwk   = integer work array of length n+1
c
c notes:
c
c        Algorithm is in place.  i.e. both:
c
c          call csrmsr (n, a, ja, ia, ao, jao, ao,jao)
c          (in which  ao, jao, are different from a, ja)
c           and
c          call csrmsr (n, a, ja, ia, a, ja, wk,iwk)
c          (in which  wk, jwk, are different from a, ja)
c        are OK.
c
c coded by Y. Saad Sep. 1989. Rechecked Feb 27, 1990.
c
      double precision a(*),ao(*),wk(n)
      integer ia(n+1),ja(*),jao(*),iwk(n+1)

      icount = 0
c
c store away diagonal elements and count nonzero diagonal elements.
c
      do 1 i=1,n
         wk(i) = 0.0
         iwk(i+1) = ia(i+1)-ia(i)
         do 2 k=ia(i),ia(i+1)-1
            if (ja(k) .eq. i) then
               wk(i) = a(k)
               icount = icount + 1
               iwk(i+1) = iwk(i+1)-1
            end if
 2       continue
 1    continue
c
c compute total length
c
      iptr = n + ia(n+1) - icount
c
c     copy backwards (to avoid collisions)
c
      do 500 ii=n,1,-1
         do 100 k=ia(ii+1)-1,ia(ii),-1
            j = ja(k)
            if (j .ne. ii) then
               ao(iptr) = a(k)
               jao(iptr) = j
               iptr = iptr-1
            end if
 100     continue
 500  continue
c
c compute pointer values and copy wk(*)
c
      jao(1) = n+2
      do 600 i=1,n
         ao(i) = wk(i)
         jao(i+1) = jao(i)+iwk(i+1)
 600  continue
      return
      end
      subroutine csrncf ( nrow, a, ja, ia, maxnz, nonz, coef, jcoef,
     &  ierr )

c*********************************************************************72
c
cc CSRNCF converts CSR to NSPCG NCF format.
c
c  Discussion:
c
c    This routine converts a matrix stored in the general A, JA, IA
c    compressed sparse row format into the Nonsymmetric Coordinate Format
c    used as storage format 5 by NSPCG.
c
c  Modified:
c
c    25 October 2008
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer NROW, the row dimension of the matrix A.
c
c    Input, real A(*), integer JA(*), IA(NROW+1), the matrix, stored in
c    compressed sparse row format.
c
c    Input, integer MAXNZ, the maximum number of nonzeros allowed for
c    in the storage of COEF and JCOEF.
c
c    Output, integer NONZ, the actual number of nonzeros encountered.
c
c    Output, real COEF(MAXNZ), the values of the matrix A in NCF format.
c
c    Output, integer JCOEF(MAXNZ,2), the row and column indices of each
c    entry in COEF.
c
c    Output, integer IERR, an error flag.
c    0 = correct return.
c    nonzero means that MAXNZ < NONZ.
c
      implicit none

      integer maxnz
      integer nrow

      double precision a(*)
      double precision coef(maxnz)
      integer i
      integer ia(nrow+1)
      integer ierr
      integer j
      integer ja(*)
      integer jcoef(maxnz,2)
      integer k
      integer k1
      integer k2
      integer nonz

      ierr = 0
c
c  Initialize COEF and JCOEF.
c
      do i = 1, maxnz
        coef(i) = 0.0D+00
      end do

      do j = 1, 2
        do i = 1, maxnz
          jcoef(i,j) = 0
        end do
      end do
c
c  The first N entries are reserved for the diagonals.
c
      do i = 1, nrow
        jcoef(i,1:2) = i
      end do

      nonz = nrow

      if ( maxnz .lt. nonz ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'CSRNCF - Fatal error!'
        write ( *, '(a)' ) '  MAXNZ < NONZ.'
        ierr = 1
        return
      end if

      do i = 1, nrow

        k1 = ia(i)
        k2 = ia(i+1) - 1

        do k = k1, k2

          if ( ja(k) .eq. i ) then

            coef(i) = coef(i) + a(k)

          else if ( 0.0D+00 .ne. a(k) ) then

            nonz = nonz + 1

            if ( maxnz .lt. nonz ) then
              write ( *, '(a)' ) ' '
              write ( *, '(a)' ) 'CSRNCF - Fatal error!'
              write ( *, '(a)' ) '  MAXNZ < NONZ.'
              ierr = 1
              return
            end if

            coef(nonz) = a(k)
            jcoef(nonz,1) = i
            jcoef(nonz,2) = ja(k)

          end if

        end do

      end do

      return
      end
      subroutine csrssk(n,imod,a,ja,ia,asky,isky,nzmax,ierr)

c*********************************************************************72
c
cc CSRSSK converts Compressed Sparse Row to Symmetric Skyline Format.
c
c CSRSSK translates a compressed sparse row or a symmetric
c sparse row format into a symmetric skyline format.
c the input matrix can be in either compressed sparse row or the
c symmetric sparse row format. The output matrix is in a symmetric
c skyline format: a double precision array containing the (active portions) of the
c rows in  sequence and a pointer to the beginning of each row.
c
c This module is NOT  in place.
c
c Coded by Y. Saad, Oct 5, 1989. Revised Feb. 18, 1991.
c
c on entry:
c
c n      = integer equal to the dimension of A.
c imod  = integer indicating the variant of skyline format wanted:
c         imod = 0 means the pointer isky points to the `zeroth'
c         element of the row, i.e., to the position of the diagonal
c         element of previous row (for i=1, isky(1)= 0)
c         imod = 1 means that itpr points to the beginning of the row.
c         imod = 2 means that isky points to the end of the row (diagonal
c                  element)
c
c a      = double precision array of size nna containing the nonzero elements
c ja      = integer array of size      nnz containing the column positions
c         of the corresponding elements in a.
c ia      = integer of size n+1. ia(k) contains the position in a, ja of
c        the beginning of the k-th row.
c nzmax = integer. must be set to the number of available locations
c         in the output array asky.
c
c on return:
c
c
c asky    = double precision array containing the values of the matrix stored in skyline
c         format. asky contains the sequence of active rows from
c         i=1, to n, an active row being the row of elemnts of
c         the matrix contained between the leftmost nonzero element
c         and the diagonal element.
c isky      = integer array of size n+1 containing the pointer array to
c         each row. The meaning of isky depends on the input value of
c         imod (see above).
c ierr  =  integer.  Error message. If the length of the
c         output array asky exceeds nzmax. ierr returns the minimum value
c         needed for nzmax. otherwise ierr=0 (normal return).
c
c Notes:
c         1) This module is NOT  in place.
c         2) even when imod = 2, length of  isky is  n+1, not n.
c
c
c first determine individial bandwidths and pointers.
c
      INTEGER NZMAX
      double precision a(*),asky(nzmax)
      integer n, imod, ierr, ia(n+1), isky(n+1), ja(*)
c
      ierr = 0
      isky(1) = 0
      do 3 i=1,n
         ml = 0
         do 31 k=ia(i),ia(i+1)-1
            ml = max(ml,i-ja(k)+1)
 31      continue
         isky(i+1) = isky(i)+ml
 3    continue
c
c     test if there is enough space  asky to do the copying.
c
      nnz = isky(n+1)
      if (nnz .gt. nzmax) then
         ierr = nnz
         return
      end if
c
c   fill asky with zeros.
c
      do 1 k=1, nnz
         asky(k) = 0.0
 1    continue
c
c     copy nonzero elements.
c
      do 4 i=1,n
         kend = isky(i+1)
         do 41 k=ia(i),ia(i+1)-1
            j = ja(k)
            if (j .le. i) asky(kend+j-i) = a(k)
 41      continue
 4    continue
c
c modify pointer according to imod if necessary.
c
      if (imod .eq. 0) return
      if (imod .eq. 1) then
         do 50 k=1, n+1
            isky(k) = isky(k)+1
 50      continue
      end if
      if (imod .eq. 2) then
         do 60 k=1, n
            isky(k) = isky(k+1)
 60      continue
      end if

      return
      end
      subroutine csrssr (nrow,a,ja,ia,nzmax,ao,jao,iao,ierr)

c*********************************************************************72
c
cc CSRSSR converts Compressed Sparse Row to Symmetric Sparse Row.
c
c CSRSSR extracts the lower triangular part of a matrix.
c It can used as a means for converting a symmetric matrix for
c which all the entries are stored in sparse format into one
c in which only the lower part is stored. The routine is in place in
c that the output matrix ao, jao, iao can be overwritten on
c the  input matrix  a, ja, ia if desired. Csrssr has been coded to
c put the diagonal elements of the matrix in the last position in
c each row (i.e. in position  ao(ia(i+1)-1   of ao and jao)
c
c On entry
c
c nrow  = dimension of the matrix a.
c a, ja,
c    ia = matrix stored in compressed row sparse format
c
c nzmax = length of arrays ao,  and jao.
c
c On return:
c
c ao, jao,
c     iao = lower part of input matrix (a,ja,ia) stored in compressed sparse
c          row format format.
c
c ierr   = integer error indicator.
c          ierr .eq. 0  means normal return
c          ierr .eq. i  means that the code has stopped when processing
c          row number i, because there is not enough space in ao, jao
c          (according to the value of nzmax)
c
      double precision a(*), ao(*), t
      integer ia(*), ja(*), iao(*), jao(*)

      ierr = 0
      ko = 0

      do  7 i=1, nrow
         kold = ko
         kdiag = 0
         do 71 k = ia(i), ia(i+1) -1
            if (ja(k)  .gt. i) goto 71
            ko = ko+1
            if (ko .gt. nzmax) then
               ierr = i
               return
            end if
            ao(ko) = a(k)
            jao(ko) = ja(k)
            if (ja(k)  .eq. i) kdiag = ko
 71      continue
         if (kdiag .eq. 0 .or. kdiag .eq. ko) goto 72
c
c     exchange
c
         t = ao(kdiag)
         ao(kdiag) = ao(ko)
         ao(ko) = t

         k = jao(kdiag)
         jao(kdiag) = jao(ko)
         jao(ko) = k
 72      iao(i) = kold+1
 7    continue
c     redefine iao(n+1)
      iao(nrow+1) = ko+1
      return
      end
      subroutine daxpy ( n, da, dx, incx, dy, incy )

c*********************************************************************72
c
cc DAXPY computes constant times a vector plus a vector.
c
c  Discussion:
c
c    This routine uses double precision real arithmetic.
c
c    This routine uses unrolled loops for increments equal to one.
c
c  Modified:
c
c    08 July 2007
c
c  Author:
c
c    Jack Dongarra
c
c  Reference:
c
c    Jack Dongarra, Jim Bunch, Cleve Moler, Pete Stewart,
c    LINPACK User's Guide,
c    SIAM, 1979,
c    ISBN13: 978-0-898711-72-1,
c    LC: QA214.L56.
c
c    Charles Lawson, Richard Hanson, David Kincaid, Fred Krogh,
c    Basic Linear Algebra Subprograms for FORTRAN usage,
c    ACM Transactions on Mathematical Software,
c    Volume 5, Number 3, pages 308-323, 1979.
c
c  Parameters:
c
c    Input, integer N, the number of elements in DX and DY.
c
c    Input, double precision DA, the multiplier of DX.
c
c    Input, double precision DX(*), the first vector.
c
c    Input, integer INCX, the increment between successive entries of DX.
c
c    Input/output, double precision DY(*), the second vector.
c    On output, DY(*) has been replaced by DY(*) + DA * DX(*).
c
c    Input, integer INCY, the increment between successive entries of DY.
c
      implicit none

      double precision da
      double precision dx(*)
      double precision dy(*)
      integer i
      integer incx
      integer incy
      integer ix
      integer iy
      integer m
      integer n

      if ( n .le. 0 ) then
        return
      end if

      if ( da .eq. 0.0d0 ) then
        return
      end if

      if(incx.eq.1.and.incy.eq.1)go to 20
c
c  code for unequal increments or equal increments
c  not equal to 1
c
      ix = 1
      iy = 1
      if(incx.lt.0)ix = (-n+1)*incx + 1
      if(incy.lt.0)iy = (-n+1)*incy + 1

      do i = 1,n
        dy(iy) = dy(iy) + da*dx(ix)
        ix = ix + incx
        iy = iy + incy
      end do

      return
c
c  code for both increments equal to 1
c
c
c  clean-up loop
c
   20 m = mod(n,4)
      if( m .eq. 0 ) go to 40
      do i = 1,m
        dy(i) = dy(i) + da*dx(i)
      end do

      if( n .lt. 4 ) return

   40 continue

      do i = m + 1, n, 4
        dy(i) = dy(i) + da*dx(i)
        dy(i + 1) = dy(i + 1) + da*dx(i + 1)
        dy(i + 2) = dy(i + 2) + da*dx(i + 2)
        dy(i + 3) = dy(i + 3) + da*dx(i + 3)
      end do

      return
      end
      subroutine dcn(AR,IA,JA,N,NE,IC,NN,IERR)

c*********************************************************************72
c
cc DCN generates sparse square matrices of type D(N,C).
c
c   PURPOSE
c
c   DCN generates sparse (square) matrices of the type
c   D(N,C).  This type of matrix has the following characteristics:
c   1's in the diagonal, three bands at the distance C above the
c   diagonal (and reappearing cyclicly under it), and a 10 x 10
c   triangle of elements in the upper right hand corner.
c   Different software libraries require different storage schemes.
c   This routine generates the matrix in the storage  by
c   indices mode.
c
c
c   Note: If A is the sparse matrix of type D(N,C), then
c
c       min|A(i,j)| = 1,     max|A(i,j)| = max(1000,N + 1)
c
c
c
c   CONTRIBUTOR: Ernest E. Rothman
c                Cornell Theory Center/Cornell National Supercomputer
c                Facility.
c                e-mail address: BITNET:   eer@cornellf
c                                INTERNET: eer@cornellf.tn.cornell.edu
c
c
c   REFERENCE
c
c   1) Zlatev, Zahari; Schaumburg, Kjeld; Wasniewski, Jerzy;
c      "A Testing Scheme for Subroutines Solving Large Linear Problems",
c       Computers and Chemistry, Vol. 5, No. 2-3, pp. 91-100, 1981.
c   2) Osterby, Ole and Zletev, Zahari;
c      "Direct Methods for Sparse Matrices";
c       Springer-Verlag 1983.
c
c
c
c   INPUT PARAMETERS
c
c   N    - Integer. The size of the square matrix.
c          N > 13 must be specified.
c
c   NN   - Integer. The dimension of integer arrays IA and JA and
c          double precision array AR. Must be at least NE.
c
c   IC   - Integer. The sparsity pattern can be changed by means of this
c          parameter.  0 < IC < N-12  must be specified.
c
c
c   OUTPUT PARAMETERS
c
c   NE   - Integer. The number of nonzero elements in the sparse matrix
c          of the type D(N,C). NE = 4*N + 55.
c
c   AR(NN) - Real array. (Double precision)
c            Stored entries of a sparse matrix to be generated by this
c            routine.
c            NN is greater then or equal to, NE, the number of
c            nonzeros including a mandatory diagonal entry for
c            each row. Entries are stored by indices.
c
c   IA(NN) - Integer array.
c            Pointers to specify rows for the stored nonzero entries
c            in AR.
c
c   JA(NN) - Integer array.
c            Pointers to specify columns for the stored nonzero entries
c            in AR.
c
c   IERR   - Error parameter is returned as zero on successful
c             execution of the routine.
c             Error diagnostics are given by means of positive values
c             of this parameter as follows:
c             IERR = 1    -  N       is out of range.
c             IERR = 2    -  IC      is out of range.
c             IERR = 3    -  NN      is out of range.
c
c
c
      double precision ar(nn)
      integer ia(nn), ja(nn), ierr
      ierr = 0
c
c
c  check the input parameters:
c
      if(n.le.13)then
         ierr = 1
         return
      end if
      if(ic .le. 0 .or. ic .ge. n-12)then
         ierr = 2
         return
      end if
      ne = 4*n+55
      if(nn.lt.ne)then
         ierr = 3
         return
      end if
c
c Begin to generate the nonzero elements as well as the row and column
c pointers:
c
      do 20 i=1,n
        ar(i) = 1.0
        ia(i) = i
        ja(i) = i
20    continue
      ilast = n
      do 30 i=1,n-ic
        it = ilast + i
        ar(it) = 1.0 + dble(i)
        ia(it) = i
        ja(it) = i+ic
30    continue
      ilast = ilast + n-ic
      do 40 i=1,n-ic-1
        it = ilast + i
        ar(it) = -dble(i)
        ia(it) = i
        ja(it) = i+ic+1
40    continue
      ilast = ilast + n-ic-1
      do 50 i=1,n-ic-2
        it = ilast + i
        ar(it) = 16.0
        ia(it) = i
        ja(it) = i+ic+2
50    continue
      ilast = ilast + n-ic-2
      icount = 0
      do 70 j=1,10
        do 60 i=1,11-j
         icount = icount + 1
         it = ilast + icount
         ar(it) = 100.0 * dble(j)
         ia(it) = i
         ja(it) = n-11+i+j
60    continue
70    continue
      icount = 0
      ilast = 55 + ilast
      do 80 i=n-ic+1,n
        icount = icount + 1
        it = ilast + icount
        ar(it) = 1.0 + dble(i)
        ia(it) = i
        ja(it) = i-n+ic
80    continue
      ilast = ilast + ic
      icount = 0
      do 90 i=n-ic,n
        icount = icount + 1
        it = ilast + icount
        ar(it) = -dble(i)
        ia(it) = i
        ja(it) = i-n+ic+1
90    continue
      ilast = ilast + ic + 1
      icount = 0
      do 100 i=n-ic-1,n
        icount = icount + 1
        it = ilast + icount
        ar(it) = 16.0
        ia(it) = i
        ja(it) = i-n+ic+2
100   continue
c     ilast = ilast + ic + 2
c     if(ilast.ne.4*n+55) then
c     write(*,*)' ilast equal to ', ilast
c     write(*,*)' ILAST, the number of nonzeros, should = ', 4*n + 55
c     stop
c     end if
c
      return
      end
      subroutine dcsort(ival, n, icnt, index, ilo, ihi)

c*********************************************************************72
c
cc DCSORT computes a sorting permutation for a vector.
c
c     Specifications for arguments:
c
c    This routine computes a permutation which, when applied to the
c    input vector ival, sorts the integers in ival in descending
c    order.  The permutation is represented by the vector index.  The
c    permuted ival can be interpreted as follows:
c      ival(index(i-1)) .ge. ival(index(i)) .ge. ival(index(i+1))
c
c    A specialized sort, the distribution counting sort, is used
c    which takes advantage of the knowledge that
c        1)  The values are in the (small) range [ ilo, ihi ]
c        2)  Values are likely to be repeated often
c
c    contributed to SPARSKIT by Mike Heroux. (Cray Research)
c
c Usage:
c
c     call dcsort( ival, n, icnt, index, ilo, ihi )
c
c Arguments:
c
c    ival  integer array (input)
c          On entry, ia is an n dimensional array that contains
c          the values to be sorted.  ival is unchanged on exit.
c
c    n     integer (input)
c          On entry, n is the number of elements in ival and index.
c
c    icnt  integer (work)
c          On entry, is an integer work vector of length
c          (ihi - ilo + 1).
c
c    index integer array (output)
c          On exit, index is an n-length integer vector containing
c          the permutation which sorts the vector ival.
c
c    ilo   integer (input)
c          On entry, ilo is .le. to the minimum value in ival.
c
c    ihi   integer (input)
c          On entry, ihi is .ge. to the maximum value in ival.
c
c Remarks:
c
c         The permutation is NOT applied to the vector ival.
c
c Author:
c
c    Michael Heroux
c    Sandra Carney
c       Mathematical Software Research Group
c       Cray Research, Inc.
c
c References:
c    Knuth, Donald E., "The Art of Computer Programming, Volume 3:
c    Sorting and Searching," Addison-Wesley, Reading, Massachusetts,
c    1973, pp. 78-79.
c
c Revision history:
c    05/09/90: Original implementation.  A variation of the
c              Distribution Counting Sort recommended by
c              Sandra Carney. (Mike Heroux)
c
c
      integer n, ilo, ihi, ival(n), icnt(ilo:ihi), index(n)
      integer i, j, ivalj

      do i = ilo, ihi
        icnt(i) = 0
      end do

      do i = 1, n
        icnt(ival(i)) = icnt(ival(i)) + 1
      end do

      do i = ihi-1,ilo,-1
        icnt(i) = icnt(i) + icnt(i+1)
      end do

      do j = n, 1, -1
        ivalj = ival(j)
        index(icnt(ivalj)) = j
        icnt(ivalj) = icnt(ivalj) - 1
      end do

      return
      end
      function ddot ( n, dx, incx, dy, incy )

c*********************************************************************72
c
cc DDOT forms the dot product of two vectors.
c
c  Discussion:
c
c    This routine uses double precision real arithmetic.
c
c    This routine uses unrolled loops for increments equal to one.
c
c  Modified:
c
c    07 July 2007
c
c  Author:
c
c    Jack Dongarra
c
c  Reference:
c
c    Jack Dongarra, Jim Bunch, Cleve Moler, Pete Stewart,
c    LINPACK User's Guide,
c    SIAM, 1979,
c    ISBN13: 978-0-898711-72-1,
c    LC: QA214.L56.
c
c    Charles Lawson, Richard Hanson, David Kincaid, Fred Krogh,
c    Basic Linear Algebra Subprograms for FORTRAN usage,
c    ACM Transactions on Mathematical Software,
c    Volume 5, Number 3, pages 308-323, 1979.
c
c  Parameters:
c
c    Input, integer N, the number of entries in the vectors.
c
c    Input, double precision DX(*), the first vector.
c
c    Input, integer INCX, the increment between successive entries in DX.
c
c    Input, double precision DY(*), the second vector.
c
c    Input, integer INCY, the increment between successive entries in DY.
c
c    Output, double precision DDOT, the sum of the product of the
c    corresponding entries of DX and DY.
c
      implicit none

      double precision ddot
      double precision dx(*),dy(*),dtemp
      integer i,incx,incy,ix,iy,m,n

      ddot = 0.0d0
      dtemp = 0.0d0
      if(n.le.0)return
      if(incx.eq.1.and.incy.eq.1)go to 20
c
c  code for unequal increments or equal increments not equal to 1
c
      ix = 1
      iy = 1
      if(incx.lt.0)ix = (-n+1)*incx + 1
      if(incy.lt.0)iy = (-n+1)*incy + 1
      do i = 1,n
        dtemp = dtemp + dx(ix)*dy(iy)
        ix = ix + incx
        iy = iy + incy
      end do
      ddot = dtemp
      return
c
c  code for both increments equal to 1
c
c
c  clean-up loop
c
   20 m = mod(n,5)
      if( m .eq. 0 ) go to 40
      do i = 1,m
        dtemp = dtemp + dx(i)*dy(i)
      end do
      if( n .lt. 5 ) go to 60
   40 continue
      do i = m+1, n, 5
        dtemp = dtemp + dx(i)*dy(i) + dx(i + 1)*dy(i + 1) +
     &   dx(i + 2)*dy(i + 2) + dx(i + 3)*dy(i + 3) + dx(i + 4)*dy(i + 4)
      end do

   60 ddot = dtemp

      return
      end
      subroutine diacsr (n,job,idiag,diag,ndiag,ioff,a,ja,ia)

c*********************************************************************72
c
cc DIACSR converts diagonal format to compressed sparse row.
c
c DIACSR extracts the idiag most important diagonals from the
c input matrix a, ja, ia, i.e, those diagonals of the matrix which have
c the largest number of nonzero elements. If requested (see job),
c the rest of the matrix is put in a the output matrix ao, jao, iao
c
c on entry:
c
c n      = integer. dimension of the matrix a.
c job      = integer. job indicator with the following meaning.
c         if (job .eq. 0) then check for each entry in diag
c         whether this entry is zero. If it is then do not include
c         in the output matrix. Note that the test is a test for
c         an exact arithmetic zero. Be sure that the zeros are
c         actual zeros otherwise this would not
c         work.
c
c idiag = integer equal to the number of diagonals to be extracted.
c         Note: on return idiag may be modified.
c
c diag  = double precision array of size (ndiag x idiag) containing the diagonals
c         of A on return.
c
c ndiag = integer equal to the first dimension of array diag.
c
c ioff  = integer array of length idiag, containing the offsets of the
c           diagonals to be extracted.
c
c on return:
c
c a,
c ja,
c ia    = matrix stored in a, ja, ia, format
c
c Note:
c the arrays a and ja should be of length n*idiag.
c
      double precision diag(ndiag,idiag), a(*), t
      integer ia(*), ja(*), ioff(*)

      ia(1) = 1
      ko = 1
      do 80 i=1, n
         do 70 jj = 1, idiag
            j = i+ioff(jj)
            if (j .lt. 1 .or. j .gt. n) goto 70
            t = diag(i,jj)
            if (job .eq. 0 .and. t .eq. 0.0) goto 70
            a(ko) = t
            ja(ko) = j
            ko = ko+1
 70      continue
         ia(i+1) = ko
 80   continue
      return
      end
      subroutine diamua (nrow,job, a, ja, ia, diag, b, jb, ib)

c*********************************************************************72
c
cc DIAMUA performs the matrix by matrix product B = Diag * A.
c
c on entry:
c
c nrow      = integer. The row dimension of A
c
c job   = integer. job indicator. Job=0 means get array b only
c         job = 1 means get b, and the integer arrays ib, jb.
c
c a,
c ja,
c ia   = Matrix A in compressed sparse row format.
c
c diag = diagonal matrix stored as a vector dig(1:n)
c
c on return:
c
c
c b,
c jb,
c ib      = resulting matrix B in compressed sparse row sparse format.
c
c Notes:
c
c 1)        The column dimension of A is not needed.
c 2)        algorithm in place (B can take the place of A).
c           in this case use job=0.
c
      double precision a(*), b(*), diag(nrow), scal
      integer ja(*),jb(*), ia(nrow+1),ib(nrow+1)

      do 1 ii=1,nrow
c
c     normalize each row
c
         k1 = ia(ii)
         k2 = ia(ii+1)-1
         scal = diag(ii)
         do 2 k=k1, k2
            b(k) = a(k)*scal
 2       continue
 1    continue

      if (job .eq. 0) return

      ib(1) = ia(1)
      do 3 ii=1, nrow
         ib(ii) = ia(ii)
         do 31 k=ia(ii),ia(ii+1)-1
            jb(k) = ja(k)
 31      continue
 3    continue
      return
      end
      subroutine diapos  (n,ja,ia,idiag)

c*********************************************************************72
c
cc DIAPOS returns the positions of the diagonal elements of a sparse matrix.
c
c on entry:
c
c
c n      = integer. row dimension of the matrix a.
c a,ja,
c    ia = matrix stored compressed sparse row format. a array skipped.
c
c on return:
c
c idiag  = integer array of length n. The i-th entry of idiag
c          points to the diagonal element a(i,i) in the arrays
c          a, ja. (i.e., a(idiag(i)) = element A(i,i) of matrix A)
c          if no diagonal element is found the entry is set to 0.
c
c           Y. Saad, March, 1990
c
      integer ia(n+1), ja(*), idiag(n)

      do 1 i=1, n
         idiag(i) = 0
 1    continue
c
c     sweep through data structure.
c
      do  6 i=1,n
         do 51 k= ia(i),ia(i+1) -1
            if (ja(k) .eq. i) idiag(i) = k
 51      continue
 6    continue
      return
      end
      subroutine dinfo1(n,iout,a,ja,ia,valued,
     *                      title,key,type,ao,jao,iao)

c*********************************************************************72
c
cc DINFO1 computes and prints matrix statistics.
c
c  SPARSKIT:  ELEMENTARY INFORMATION ROUTINE.
c
c info1 obtains a number of statistics on a sparse matrix and writes
c it into the output unit iout. The matrix is assumed
c to be stored in the compressed sparse COLUMN format sparse a, ja, ia
c
c Modified Nov 1, 1989. 1) Assumes A is stored in column
c format. 2) Takes symmetry into account, i.e., handles Harwell-Boeing
c            matrices correctly.
c          ***  (Because of the recent modification the words row and
c            column may be mixed-up at occasions... to be checked...
c
c bug-fix July 25: 'upper' 'lower' mixed up in formats 108-107.
c
c On entry :
c
c n      = integer. column dimension of matrix
c iout  = integer. unit number where the information it to be output.
c a      = double precision array containing the nonzero elements of the matrix
c        the elements are stored by columns in order
c        (i.e. column i comes before column i+1, but the elements
c         within each column can be disordered).
c ja      = integer array containing the row indices of elements in a
c ia      = integer array containing of length n+1 containing the
c         pointers to the beginning of the columns in arrays a and ja.
c        It is assumed that ia(*) = 1 and ia(n+1) = nzz+1.
c
c valued= logical equal to .true. if values are provided and .false.
c         if only the pattern of the matrix is provided. (in that
c         case a(*) and ao(*) are dummy arrays.
c
c title = a 72-character title describing the matrix
c         NOTE: The first character in title is ignored (it is often
c         a one).
c
c key   = an 8-character key for the matrix
c type  = a 3-character string to describe the type of the matrix.
c         see harwell/Boeing documentation for more details on the
c         above three parameters.
c
c on return
c
c 1) elementary statistics on the matrix is written on output unit
c    iout. See below for detailed explanation of typical output.
c 2) the entries of a, ja, ia are sorted.
c
c
c
c ao      = double precision array of length nnz used as work array.
c jao      = integer work array of length max(2*n+1,nnz)
c iao   = integer work array of length n+1
c
c Note  : title, key, type are the same paramaters as those
c         used for Harwell-Bowing matrices.
c
c Output description:
c
c *** The following info needs to be updated.
c
c + A header containing the Title, key, type of the matrix and, if values
c   are not provided a message to that effect.
c
c    * SYMMETRIC STRUCTURE MEDIEVAL RUSSIAN TOWNS
c    *                    Key = RUSSIANT , Type = SSA
c    * No values provided - Information of pattern only
c
c
c +  dimension n, number of nonzero elements nnz, average number of
c    nonzero elements per column, standard deviation for this average.
c +  if the matrix is upper or lower triangular a message to that effect
c    is printed. Also the number of nonzeros in the strict upper
c    (lower) parts and the main diagonal are printed.
c +  weight of longest column. This is the largest number of nonzero
c    elements in a column encountered. Similarly for weight of
c    largest/smallest row.
c +  lower dandwidth as defined by
c          ml = max ( i-j, / all  a(i,j).ne. 0 )
c +  upper bandwidth as defined by
c          mu = max ( j-i, / all  a(i,j).ne. 0 )
c    NOTE that ml or mu can be negative. ml .lt. 0 would mean
c    that A is confined to the strict upper part above the diagonal
c    number -ml. Similarly for mu.
c
c +  maximun bandwidth as defined by
c    Max (  Max [ j ; a(i,j) .ne. 0 ] - Min [ j ; a(i,j) .ne. 0 ] )
c     i
c +  average bandwidth = average over all columns of the widths each column.
c
c +  If there are zero columns /or rows a message is printed
c    giving the number of such columns/rows.
c
c +  matching elements in A and transp(A) :this counts the number of
c    positions (i,j) such that if a(i,j) .ne. 0 then a(j,i) .ne. 0.
c    if this number is equal to nnz then the matrix is symmetric.
c +  Relative symmetry match : this is the ratio of the previous integer
c    over nnz. If this ratio is equal to one then the matrix has a
c    symmetric structure.
c
c +  average distance of a given element from the diagonal, standard dev.
c    the distance of a(i,j) is defined as iabs(j-i).
c
c +  Frobenious norm of A
c    Frobenious norm of 0.5*(A + transp(A))
c    Frobenious norm of 0.5*(A - transp(A))
c    these numbers provide information on the degree of symmetry
c    of the matrix. If the norm of the nonsymmetric part is
c    zero then the matrix is symmetric.
c
c + 90% of matrix is in the band of width k, means that
c   by moving away and in a symmetric manner from the main
c   diagonal you would have to include exactly k diagonals
c   (k is always odd), in order to include 90% of the nonzero
c   elements of A.  The same thing is then for 80%.
c
c + The total number of nonvoid diagonals, i.e., among the
c   2n-1 diagonals of the matrix which have at least one nonxero
c   element.
c
c +  Most important diagonals. The code selects a number of k
c    (k .le. 10) diagonals that are the most important ones, i.e.
c    that have the largest number of nonzero elements. Any diagonal
c    that has fewer than 1% of the nonzero elements of A is dropped.
c    the numbers printed are the offsets with respect to the
c    main diagonal, going from left tp right.
c    Thus 0 means the main diagonal -1 means the subdiagonal, and
c    +10 means the 10th upper diagonal.
c +  The accumulated percentages in the next line represent the
c    percentage of the nonzero elements represented by the diagonals
c    up the current one put together.
c    Thus:
c    *  The 10 most important diagonals are (offsets)    :
c    *     0     1     2    24    21     4    23    22    20    19
c    *  The accumulated percentages they represent are   :
c    *  40.4  68.1  77.7  80.9  84.0  86.2  87.2  88.3  89.4  90.4
c
c    shows the offsets of the most important  diagonals and
c    40.4 represent ratio of the number of nonzero elements in the
c    diagonal zero (main diagonal) over the total number of nonzero
c    elements. the second number indicates that the diagonal 0 and the
c    diagonal 1 together hold 68.1% of the matrix, etc..
c
c +  Block structure:
c    if the matrix has a block structure then the block size is found
c    and printed. Otherwise the info1 will say that the matrix
c    does not have a block structure. Note that block structure has
c    a very specific meaning here. the matrix has a block structure
c    if it consists of square blocks that are dense. even if there
c    are zero elements in the blocks  they should be represented
c    otherwise it would be possible to determine the block size.
c
      implicit double precision (a-h,o-z)
        double precision a(*),ao(*)
        integer ja(*),ia(n+1),jao(*),iao(n+1)
        character title*72,key*8,type*3
        logical valued

      PARAMETER (IPAR1=1)

      double precision dcount(20),amx
      integer ioff(20)
      character*61 tmpst
      logical sym

      write (iout,99)
        write (iout,97) title(2:72), key, type
 97     format(2x,' * ',a71,' *'/,
     *         2x,' *',20x,'Key = ',a8,' , Type = ',a3,25x,' *')
      if (.not. valued) write (iout,98)
 98    format(2x,' * No values provided - Information on pattern only',
     *   23x,' *')

      nnz = ia(n+1)-ia(1)
      sym = (type(2:2) .eq. 'S')

        write (iout, 99)
        write(iout, 100) n, nnz
c first average and standard deviation
        av = dble(nnz)/dble(n)
c av will be coorrected later.
      if (sym) av = 2.0*av-1.0
      job = 0
      if (valued) job = 1
      ipos = 1
        call csrcsc(n, job, ipos, a, ja, ia, ao, jao, iao)
        call csrcsc(n, job, ipos, ao, jao, iao, a, ja, ia)
c bandwidth
        iband = 0
c number of nonzero elements in lower part
      nupper = 0
c number of nonzero elements in diagonal
      ndiag  = 0
c distance of an element from diagonal.
      dist   = 0.0
c number of diagonally dominant columns
      nddomc = 0
c number of diagonally dominant rows
      nddomr = 0
c max length of rows
      nzmaxc = 0
c min length of column
      nzminc = n
c max length of rows
      nzmaxr = 0
c min length of rows
      nzminr = n
c number of zero columns
      nzcol  = 0
c number of zero column
      nzrow  = 0
c lower and upper band
      ml = -n
      mu = -n
c average bandwidth
        bndav = 0.0
c standard deviation for average nonzero elements --
      st = 0.0
c dianrm = Frobenius norm of the diagonal (used only in symmetric case)
        dianrm = 0.0
c nskyu = skyline storage for upper part
      nskyu = 0
c nskyl = skyline storage for lower part
      nskyl = 0
c
c computing max bandwith, max number of nonzero elements per column
c min nonzero elements per column, column and col. diagonal dominance
c occurences, average distance of an element from diagonal, number of
c elemnts in lower and upper parts, ...
c
        do 3 i=1,n
                j0 = ia(i)
                j1 = ia(i+1) - 1
            j0r = iao(i)
            j1r = iao(i+1)-1
c
c bandwidth info:
c
              jminc = ja(j0)
              jmaxc = ja(j1)
            jminr = jao(j0r)
            if (sym) jminc = jminr
c
            nskyl = nskyl + i-jminr + 1
            nskyu = nskyu + i-jminc + 1
c
            ml = max0(ml,i-jminc)
              mu = max0(mu,jmaxc-i)
                iband = max0(iband,jmaxc-jminc+1)
              bndav = bndav+dble( jmaxc-jminc+1)
c max min number of nonzero elements per column .
            lenr = j1r+1-j0r
                if (lenr .le. 0) nzrow = nzrow+1
            nzmaxr = max0(nzmaxr,lenr)
            nzminr = min0(nzminr,lenr)
c indiag = nonzero diagonal element indicator
            indiag = 0
          do 31 k=j0, j1
                j=ja(k)
                if (j .lt. i) nupper = nupper+1
               if (j .eq. i) indiag = 1
                dist = dist + dble(iabs(j-i) )
 31        continue
              ndiag = ndiag + indiag
c max/ min number of nonzero elements per column --
            lenc = j1+1-j0
            if (sym) lenc = lenc + lenr - indiag
                if (lenc .le. 0) nzcol = nzcol +1
            nzmaxc = max0(nzmaxc,lenc)
            nzminc = min0(nzminc,lenc)
              st = st + (dble(lenc) - av)**2
c  diagonal dominance
c
          if (valued) then
        dsumc = 0.0
        aii = 0.0
          do 32 k=j0,j1
        j = ja(k)
        if (j .eq. i) aii = abs(a(k))
        dsumc = dsumc + abs(a(k))
 32       continue

        dianrm = dianrm + aii*aii

        dsumr = 0.0
          do 33 k=iao(i), iao(i+1)-1
 33        dsumr = dsumr+abs(ao(k))
        if (sym) then
           if (dsumr+dsumc .le. 3.0*aii) nddomc = nddomc + 1
           else
             if (dsumc .le. 2.0*aii) nddomc = nddomc + 1
             if (dsumr .le. 2.0*aii) nddomr = nddomr+1
        end if
          end if
 3        continue
        if (sym) nddomr = nddomc

        nlower = nnz-nupper-ndiag
c write bandwidth info
        dist = dist/dble(nnz)
c
c if ndiag.ne.n then we should correct av and std in symmetric case
c
        if ((sym) .and. ndiag .ne. n) then
           eps = dble(ndiag-n)/dble(n)
           av = av + eps
           st = st-eps*eps
          end if

          st = sqrt( st / dble(n) )
        bndav = bndav/dble(n)
c  write out info
        if (sym)  nupper = nlower
          write(iout, 101) av, st
          if (nlower .eq. 0 ) write(iout, 105)
 1        if (nupper .eq. 0) write(iout, 106)
          write(iout, 107) nlower
          write(iout, 108) nupper
          write(iout, 109) ndiag
c
          write(iout, 1020) nzmaxc, nzminc
        if (.not. sym) write(iout, 1021) nzmaxc, nzminc
c
        if (nzcol .ne. 0) write(iout,116) nzcol
        if (nzrow .ne. 0) write(iout,115) nzrow
c
c normalize various results of above loop
c
        ddomr = dble(nddomc)/ dble(n)
        ddomc = dble(nddomr)/ dble(n)
c
c symmetry and near symmetry - Frobenius  norms
c
          st  = 0.0
          tan = 0.0
          tas = 0.0
          std = 0.0
          imatch = 0
c
c main loop for symmetry detection and frobenius norms.
c
          do 6 i=1,n
             k1 = ia(i)
             k2 = iao(i)
             k1max = ia(i+1) - 1
             k2max = iao(i+1) - 1
             do 61 k=k1, k1max
                std=std+(dist-dble(iabs(ja(k)-i)))**2
 61          continue

             if (sym) goto 6
 5           if (k1 .gt. k1max .or. k2 .gt. k2max) goto 6

             j1 = ja(k1)
             j2 = jao(k2)
             if (j1 .ne. j2 ) goto 51
             imatch = imatch + 1
             if (valued) then
                tas = tas + (a(k1)+ao(k2))**2
                tan = tan + (a(k1)-ao(k2))**2
                st  = st+a(k1)**2
             end if

 51          k1 = k1+1
             k2 = k2+1
             if (j1 .lt. j2)  k2 = k2 - 1
             if (j1. gt. j2)  k1 = k1 - 1
             goto 5
 6        continue

          if (sym) imatch = nnz
          av = dble(imatch)/dble(nnz)
          std = sqrt(std/ dble(nnz))
c  max abs value in a
          if (valued) then
             amx = 0.0
             ta = 0.0
             do 7 k=1, nnz
                ta = ta+a(k)**2
                amx = amax1(amx, abs(a(k)) )
 7           continue
             if (sym) then
                ta = sqrt( 2.0*ta - dianrm )
                tas = ta
                tan = 0.0
             else
                st = ta-st
                tas = 0.5*sqrt(tas + st)
                tan = 0.5*sqrt(tan + st)
                ta = sqrt(ta)
             end if
          end if

          write (iout,103) imatch, av, dist, std
          write(iout,96)
          if (valued) then
             write(iout,104) ta, tas, tan, amx, ddomr, ddomc
             write (iout,96)
          end if
c
c  bandedness- main diagonals
c
      n2 = n+n-1
      do 8 i=1, n2
             jao(i) = 0
 8        continue
          do 9 i=1, n
             k1 = ia(i)
             k2 = ia(i+1) -1
             do 91 k=k1, k2
                j = ja(k)
                jao(n+i-j) = jao(n+i-j) +1
 91          continue
 9        continue
          iacc = jao(n)
          jb1 = 0
          jb2 = 0
          j = 0
 92       j = j+1
          iacc = iacc + jao(n+j)+ jao(n-j)
          if (iacc*100 .le. nnz*80)  jb1 = jb1+1
 93       if (iacc*100 .le. nnz*90) then
             jb2 = jb2+1
             goto 92
          end if
c     write bandwidth information .
c
          write(iout,117)  ml, mu, iband, bndav

          nsky = nskyl+nskyu-n
          if (sym) nsky = nskyl

          write(iout,1175) nsky

          write (iout,112) 2*jb2+1, 2*jb1+1
c
c     count the number of nonzero diagonals.
          nzdiag = 0
          do 42 i=1, n2
             if (jao(i) .ne. 0) nzdiag=nzdiag+1
 42       continue

          ndiag = 10
          ndiag = min0(n2,ndiag)
          itot  = 0
          ii = 0
          idiag = 0
c     sort diagonals by decreasing order of weights.
 40       jmax = 0
          i    = 1
          do 41 k=1, n2
             j = jao(k)
             if (j .lt. jmax) goto 41
             i = k
             jmax = j
 41       continue
c     permute
c     save offsets and accumulated count if diagonal is acceptable
c     (if it has at least ipar1*nnz/100 nonzero elements)
c     quite if no more acceptable diagonals --
c
          if (jmax*100 .lt. ipar1*nnz) goto 4
          ii = ii+1
          ioff(ii) = i-n
          jao(i)   = - jmax
          itot = itot + jmax
          dcount(ii) = dble(100*itot)/dble(nnz)
          if (ii .lt. ndiag) goto 40
 4        continue
          ndiag = ii
c
c         t = dble (icount) / dble (nnz)
c     write in internat file tmpst first.
          write (iout,118) nzdiag
          write (tmpst,'(10i6)') (ioff(j),j=1,ndiag)
          write (iout,110) ndiag,tmpst
          write (tmpst,'(10f6.1)')(dcount(j), j=1,ndiag)
          write (iout,111) tmpst
          write (iout, 96)
c     jump to next page -- optional //
c     write (iout,'(1h1)')
c
c     determine block size if matrix is a block matrix..
c
          call blkfnd (n, ja, ia, nblk)
          if (nblk .le. 1) then
             write(iout,113)
          else
             write(iout,114) nblk
          end if
          write (iout,96)
c
c done. Next define all the formats
c
 99    format (2x,38(2h *))
 96    format (6x,' *',65(1h-),'*')

 100   format(
     * 6x,' *  Dimension N                                      = ',
     * i10,'  *'/
     * 6x,' *  Number of nonzero elements                       = ',
     * i10,'  *')
 101   format(
     * 6x,' *  Average number of nonzero elements/Column        = ',
     * f10.4,'  *'/
     * 6x,' *  Standard deviation for above average             = ',
     * f10.4,'  *')

 1020       format(
     * 6x,' *  Weight of longest column                         = ',
     * i10,'  *'/
     * 6x,' *  Weight of shortest column                        = ',
     * i10,'  *')
 1021       format(
     * 6x,' *  Weight of longest row                            = ',
     * i10,'  *'/
     * 6x,' *  Weight of shortest row                           = ',
     * i10,'  *')
 117        format(
     * 6x,' *  Lower bandwidth  (max: i-j, a(i,j) .ne. 0)       = ',
     * i10,'  *'/
     * 6x,' *  Upper bandwidth  (max: j-i, a(i,j) .ne. 0)       = ',
     * i10,'  *'/
     * 6x,' *  Maximum Bandwidth                                = ',
     * i10,'  *'/
     * 6x,' *  Average Bandwidth                                = ',
     * e10.3,'  *')
 1175       format(
     * 6x,' *  Number of nonzeros in skyline storage            = ',
     * i10,'  *')
 103   format(
     * 6x,' *  Matching elements in symmetry                    = ',
     * i10,'  *'/
     * 6x,' *  Relative Symmetry Match (symmetry=1)             = ',
     * f10.4,'  *'/
     * 6x,' *  Average distance of a(i,j)  from diag.           = ',
     * e10.3,'  *'/
     * 6x,' *  Standard deviation for above average             = ',
     * e10.3,'  *')
 104   format(
     * 6x,' *  Frobenius norm of A                              = ',
     * e10.3,'  *'/
     * 6x,' *  Frobenius norm of symmetric part                 = ',
     * e10.3,'  *'/
     * 6x,' *  Frobenius norm of nonsymmetric part              = ',
     * e10.3,'  *'/
     * 6x,' *  Maximum element in A                             = ',
     * e10.3,'  *'/
     * 6x,' *  Percentage of weakly diagonally dominant rows    = ',
     * e10.3,'  *'/
     * 6x,' *  Percentage of weakly diagonally dominant columns = ',
     * e10.3,'  *')
 105        format(
     * 6x,' *  The matrix is lower triangular ...       ',21x,' *')
 106        format(
     * 6x,' *  The matrix is upper triangular ...       ',21x,' *')
 107        format(
     * 6x,' *  Nonzero elements in strict lower part            = ',
     * i10,'  *')
 108       format(
     * 6x,' *  Nonzero elements in strict upper part            = ',
     * i10,'  *')
 109       format(
     * 6x,' *  Nonzero elements in main diagonal                = ',
     * i10,'  *')
 110   format(6x,' *  The ', i2, ' most important',
     *         ' diagonals are (offsets)    : ',10x,'  *',/,
     * 6x,' *',a61,3x,' *')

 111   format(6x,' *  The accumulated percentages they represent are ',
     * '  : ', 10x,'  *',/,
     * 6x,' *',a61,3x,' *')
c 111      format(
c     * 6x,' *  They constitute the following % of A             = ',
c     * f8.1,' %  *')
 112      format(
     * 6x,' *  90% of matrix is in the band of width            = ',
     * i10,'  *',/,
     * 6x,' *  80% of matrix is in the band of width            = ',
     * i10,'  *')
 113       format(
     * 6x,' *  The matrix does not have a block structure ',19x,
     *    ' *')
 114       format(
     * 6x,' *  Block structure found with block size            = ',
     * i10,'  *')
 115       format(
     * 6x ' *  There are zero rows. Number of such rows         = ',
     * i10,'  *')
 116       format(
     * 6x ' *  There are zero columns. Number of such columns   = ',
     * i10,'  *')
 118       format(
     * 6x ' *  The total number of nonvoid diagonals is         = ',
     * i10,'  *')
      RETURN
      end
      subroutine diric (nx,nint,a,ja,ia, f)

c*********************************************************************72
c
cc DIRIC accounts for Dirichlet boundary conditions.
c
c
      implicit double precision  (a-h,o-z)
      dimension a(*),ia(*),ja(*),f(*)
c call extract from UNARY
      call submat (nx,1,1,nint,1,nint,a,ja,ia,nr,nc,a,ja,ia)
      return
      END
      subroutine dlauny(x,y,nodes,elmnts,nemax,nelmnt)

c*********************************************************************72
c
cc DLAUNY is a simple, nonoptimal Delaunay triangulation code.
c
c code written by P.K. Sweby
c
c     Performs a Delaunay triangularisation of a region given a set
c     of mesh points.
c       X,Y    :- 1D arrays holding coordinates of mesh points.
c                 dimensioned AT LEAST NODES+3.
c       NODES  :- number of mesh points.
c       ELMNTS :- INTEGER array, dimensioned NEMAX x 3, which on exit
c                 contains the index of global nodes associated with
c                 each element.
c       NELMNT :- on exit contains the number of elements in the
c                 triangularisation.
c
c                                       P.K.Sweby
c
      IMPLICIT double precision (A-H,O-Z)

      INTEGER ELMNTS
      DIMENSION X(NODES),Y(NODES),ELMNTS(NEMAX,3)

      PI=4.0*ATAN(1.0)
c
c     Calculate artificial nodes NODES+i i=1,2,3,4 and construct first
c     two (artificial) elements.
c
      XMIN=X(1)
      XMAX=X(1)
      YMIN=Y(1)
      YMAX=Y(1)
      DO 10 I=2,NODES
      XMIN=MIN(XMIN,X(I))
      XMAX=MAX(XMAX,X(I))
      YMIN=MIN(YMIN,Y(I))
      YMAX=MAX(YMAX,Y(I))
 10   CONTINUE
      DX=XMAX-XMIN
      DY=YMAX-YMIN
      XL=XMIN-4.0*DX
      XR=XMAX+4.0*DX
      YL=YMIN-4.0*DY
      YR=YMAX+4.0*DY
      X(NODES+1)=XL
      Y(NODES+1)=YL
      X(NODES+2)=XL
      Y(NODES+2)=YR
      X(NODES+3)=XR
      Y(NODES+3)=YR
      X(NODES+4)=XR
      Y(NODES+4)=YL
      ELMNTS(1,1)=NODES+1
      ELMNTS(1,2)=NODES+2
      ELMNTS(1,3)=NODES+3
      ELMNTS(2,1)=NODES+3
      ELMNTS(2,2)=NODES+4
      ELMNTS(2,3)=NODES+1
      NELMNT=2
      DO 90 IN=1,NODES
c
c     Add one mesh point at a time and remesh locally if necessary
c
      NDEL=0
      NEWEL=0
      DO 40 IE=1,NELMNT
c
c     Is point IN insided circumcircle of element IE ?
c
      I1=ELMNTS(IE,1)
      I2=ELMNTS(IE,2)
      I3=ELMNTS(IE,3)
      X2=X(I2)-X(I1)
      X3=X(I3)-X(I1)
      Y2=Y(I2)-Y(I1)
      Y3=Y(I3)-Y(I1)
      Z=(X2*(X2-X3)+Y2*(Y2-Y3))/(Y2*X3-Y3*X2)
      CX=0.5*(X3-Z*Y3)
      CY=0.5*(Y3+Z*X3)
      R2=CX**2+CY**2
      RN2=((X(IN)-X(I1)-CX)**2+(Y(IN)-Y(I1)-CY)**2)
      IF(RN2.GT.R2)GOTO 40
c
c     Yes it is inside,create new elements and mark old for deletion.
c
      DO 30 J=1,3
      DO 20 K=1,3
      ELMNTS(NELMNT+NEWEL+J,K)=ELMNTS(IE,K)
 20   CONTINUE
      ELMNTS(NELMNT+NEWEL+J,J)=IN
 30   CONTINUE
      NEWEL=NEWEL+3
      ELMNTS(IE,1)=0
      NDEL=NDEL+1

 40   CONTINUE
c
c     If IN was inside circumcircle of more than 1 element then will
c     have created 2 identical new elements: delete them both.
c
      IF(NDEL.GT.1)THEN
          DO 60 IE=NELMNT+1,NELMNT+NEWEL-1
          DO 60 JE=IE+1,NELMNT+NEWEL
          MATCH=0
          DO 50 K=1,3
          DO 50 L=1,3
          IF(ELMNTS(IE,K).EQ.ELMNTS(JE,L))MATCH=MATCH+1
 50       CONTINUE
          IF(MATCH.EQ.3)THEN
              ELMNTS(IE,1)=0
              ELMNTS(JE,1)=0
              NDEL=NDEL+2
          end if
 60       CONTINUE
      end if
c
c     Delete any elements
c
      NN=NELMNT+NEWEL
      IE=1
 70   CONTINUE
      IF(ELMNTS(IE,1).EQ.0)THEN
          DO 80 J=IE,NN-1
          DO 80 K=1,3
          ELMNTS(J,K)=ELMNTS(J+1,K)
 80       CONTINUE
          NN=NN-1
          IE=IE-1
      end if
      IE=IE+1
      IF(IE.LE.NN)GOTO 70
      NELMNT=NN
 90   CONTINUE
c
c     Finally remove elements containing artificial nodes
c
      IE=1
 100  CONTINUE
      NART=0
      DO 110 L=1,3
      IF(ELMNTS(IE,L).GT.NODES)NART=NART+1
 110  CONTINUE
      IF(NART.GT.0)THEN
          DO 120 J=IE,NN-1
          DO 120 K=1,3
          ELMNTS(J,K)=ELMNTS(J+1,K)
 120      CONTINUE
          NELMNT=NELMNT-1
          IE=IE-1
      end if
      IE=IE+1
      IF(IE.LE.NELMNT)GOTO 100
      RETURN
      END
      subroutine dnscsr(nrow,ncol,nzmax,dns,ndns,a,ja,ia,ierr)

c*********************************************************************72
c
cc DNSCSR converts Dense to Compressed Row Sparse format.
c
c converts a densely stored matrix into a row orientied
c compactly sparse matrix. ( reverse of csrdns )
c Note: this routine does not check whether an element
c is small. It considers that a(i,j) is zero if it is exactly
c equal to zero: see test below.
c
c on entry:
c
c
c nrow      = row-dimension of a
c ncol      = column dimension of a
c nzmax = maximum number of nonzero elements allowed. This
c         should be set to be the lengths of the arrays a and ja.
c dns   = input nrow x ncol (dense) matrix.
c ndns      = first dimension of dns.
c
c on return:
c
c
c a, ja, ia = value, column, pointer  arrays for output matrix
c
c ierr      = integer error indicator:
c         ierr .eq. 0 means normal retur
c         ierr .eq. i means that the the code stopped while
c         processing row number i, because there was no space left in
c         a, and ja (as defined by parameter nzmax).
c
      double precision dns(ndns,*),a(*)
      integer ia(*),ja(*)

      ierr = 0
      next = 1
      ia(1) = 1
      do 4 i=1,nrow
         do 3 j=1, ncol
            if (dns(i,j) .eq. 0.0) goto 3
            if (next .gt. nzmax) then
               ierr = i
               return
            end if
            ja(next) = j
            a(next) = dns(i,j)
            next = next+1
 3       continue
         ia(i+1) = next
 4    continue
      return
      end
      subroutine dperm (nrow,a,ja,ia,ao,jao,iao,perm,qperm,job)

c*********************************************************************72
c
cc DPERM permutes the rows and columns of a matrix stored in CSR format.
c
c This routine permutes the rows and columns of a matrix stored in CSR
c format. i.e., it computes P A Q, where P, Q are permutation matrices.
c P maps row i into row perm(i) and Q maps column j into column qperm(j):
c      a(i,j)    becomes   a(perm(i),qperm(j)) in new matrix
c In the particular case where Q is the transpose of P (symmetric
c permutation of A) then qperm is not needed.
c note that qperm should be of length ncol (number of columns) but this
c is not checked.
c
c Y. Saad, Sep. 21 1989 / recoded Jan. 28 1991.
c
c on entry:
c
c n       = dimension of the matrix
c a, ja,
c    ia = input matrix in a, ja, ia format
c perm       = integer array of length n containing the permutation arrays
c        for the rows: perm(i) is the destination of row i in the
c         permuted matrix -- also the destination of column i in case
c         permutation is symmetric (job .le. 2)
c
c qperm      = same thing for the columns. This should be provided only
c         if job=3 or job=4, i.e., only in the case of a nonsymmetric
c        permutation of rows and columns. Otherwise qperm is a dummy
c
c job      = integer indicating the work to be done:
c * job = 1,2 permutation is symmetric  Ao :== P * A * transp(P)
c             job = 1      permute a, ja, ia into ao, jao, iao
c             job = 2 permute matrix ignoring real values.
c * job = 3,4 permutation is non-symmetric  Ao :== P * A * Q
c             job = 3      permute a, ja, ia into ao, jao, iao
c             job = 4 permute matrix ignoring real values.
c
c on return:
c
c ao, jao, iao = input matrix in a, ja, ia format
c
c in case job .eq. 2 or job .eq. 4, a and ao are never referred to
c and can be dummy arguments.
c Notes:
c
c  1) algorithm is in place
c  2) column indices may not be sorted on return even  though they may be
c     on entry.
c
      integer nrow,ja(*),ia(nrow+1),jao(*),iao(nrow+1),perm(nrow),
     +        qperm(*),job
      double precision a(*),ao(*)
      integer locjob, mod
c
c     locjob indicates whether or not real values must be copied.
c
      locjob = mod(job,2)
c
c permute rows first
c
      call rperm (nrow,a,ja,ia,ao,jao,iao,perm,locjob)
c
c then permute columns
c
      locjob = 0

      if (job .le. 2) then
         call cperm (nrow,ao,jao,iao,ao,jao,iao,perm,locjob)
      else
         call cperm (nrow,ao,jao,iao,ao,jao,iao,qperm,locjob)
      end if

      return
      end
      subroutine dscaldg (n,a,ja,ia,diag,job)

c*********************************************************************72
c
cc DSCALDG scales rows by a diagonal factor.
c
c scales rows by diag where diag is either given (job=0)
c or to be computed:
c  job = 1 ,scale row i by  by  +/- max |a(i,j) | and put inverse of
c       scaling factor in diag(i),where +/- is the sign of a(i,i).
c  job = 2 scale by 2-norm of each row..
c if diag(i) = 0,then diag(i) is replaced by one
c (no scaling)..
c
c           Y. Saad, Sep. 21 1989
c
      double precision a(*), diag(*),t
      integer ia(*),ja(*)

      goto (12,11,10) job+1
 10   do 110 j=1,n
         k1= ia(j)
         k2 = ia(j+1)-1
         t = 0.0
         do 111 k = k1,k2
 111        t = t+a(k)*a(k)
 110        diag(j) = sqrt(t)
            goto 12
 11   continue
      call retmx (n,a,ja,ia,diag)

 12   do 1 j=1,n
         if (diag(j) .ne. 0.0) then
            diag(j) = 1.0/diag(j)
         else
            diag(j) = 1.0
         end if
 1    continue
      do 2 i=1,n
         t = diag(i)
         do 21 k=ia(i),ia(i+1) -1
            a(k) = a(k)*t
 21      continue
 2    continue
      return
      end
      subroutine dump(n,a,ja,ia,iout)

c*********************************************************************72
c
cc DUMP writes the matrix to a file.
c
c writes the matrix in a file, one row at a time in a nice readable
c format. This is a simple routine which is useful for debugging.
c
c on entry:
c
c n     = integer = size of matrix
c a,
c ja,
c ia    =  matrix in CSR format
c iout  = output unit number.
c
c on return:
c
c the output file iout will have written in it the matrix in
c one of two possible formats (depending on the max number of
c elements per row. the values are output with only two digits
c of accuracy (D9.2).
c
      integer ia(*),ja(*)
      double precision a(*)
c
c select mode horizontal or vertical
c
        maxr = 0
        do 1 i=1, n
           maxr = max0(maxr,ia(i+1)-ia(i))
 1      continue
        if (maxr .le. 8) then
c
c able to one row across line
c
      do 2 i=1, n
           write(iout,100) i
         k1=ia(i)
         k2 = ia(i+1)-1
         write (iout,101) (ja(k),k=k1,k2)
         write (iout,102) (a(k),k=k1,k2)
 2      continue
      else
c
c unable to one row across line. do three items at a time across line.
c
         do 3 i=1, n
            write(iout,200) i
            k1=ia(i)
            k2 = ia(i+1)-1
            write (iout,201) (ja(k),a(k),k=k1,k2)
 3       continue
      end if

 100  format(1X,35(1h-),' row',i3,1x,35(1h-) )
 101  format(' col:',8(i5,6h     :))
 102  format(' val:',8(E9.2,2h :) )
 200  format(1h ,31(1h-),' row',i3,1x,31(1h-),/
     *       3('  columns :   values   *') )
 201  format(3(1h ,i5,6h    : ,D9.2,3h  *) )
      return
      end
      subroutine dvperm (n, x, perm)

c*********************************************************************72
c
cc DVPERM performs an in-place permutation of a double precision vector.
c
c this subroutine performs an in-place permutation of a double precision vector x
c according to the permutation array perm(*), i.e., on return,
c the vector x satisfies,
c
c      x(perm(j)) :== x(j), j=1,2,.., n
c
c
c on entry:
c
c n       = length of vector x.
c perm       = integer array of length n containing the permutation  array.
c x      = input vector
c
c on return:
c
c x      = vector x permuted according to x(perm(*)) :=  x(*)
c
c
c           Y. Saad, Sep. 21 1989
c
      integer n, perm(n)
      double precision x(n)
      double precision tmp, tmp1

      init      = 1
      tmp      = x(init)
      ii        = perm(init)
      perm(init)= -perm(init)
      k         = 0
c
c loop
c
 6    k = k+1
c
c save the chased element --
c
      tmp1        = x(ii)
      x(ii)     = tmp
      next        = perm(ii)
      if (next .lt. 0 ) goto 65
c
c test for end
c
      if (k .gt. n) goto 101
      tmp       = tmp1
      perm(ii)  = - perm(ii)
      ii        = next
c
c end loop
c
      goto 6
c
c reinitilaize cycle --
c
 65   init      = init+1
      if (init .gt. n) goto 101
      if (perm(init) .lt. 0) goto 65
      tmp      = x(init)
      ii      = perm(init)
      perm(init)=-perm(init)
      goto 6

 101  continue
      do 200 j=1, n
         perm(j) = -perm(j)
 200  continue

      return
      end
      subroutine ecn(N,IC,NE,IA,JA,AR,NN,IERR)

c*********************************************************************72
c
cc ECN generates sparse (square) matrices of the type E(N,C).
c
c   PURPOSE
c
c   The subroutine generates sparse (square) matrices of the type
c   E(N,C).  This type of matrix has the following characteristics:
c   Symmetric, positive-definite, N x N matrices with 4 in the diagonal
c   and -1 in the two sidediagonal and in the two bands at the distance
c   C from the diagonal. These matrices are similar to matrices obtained
c   from using the five-point formula in the discretization of the
c   elliptic PDE.
c
c
c   Note: If A is the sparse matrix of type E(N,C), then
c
c       min|A(i,j)| = 1,     max|A(i,j)| = 4
c
c
c   CONTRIBUTOR: Ernest E. Rothman
c                Cornell Theory Center/Cornell National Supercomputer
c                Facility.
c                e-mail address: BITNET:   eer@cornellf
c                                INTERNET: eer@cornellf.tn.cornell.edu
c
c   REFERENCE
c
c   1) Zlatev, Zahari; Schaumburg, Kjeld; Wasniewski, Jerzy;
c      "A Testing Scheme for Subroutines Solving Large Linear Problems",
c       Computers and Chemistry, Vol. 5, No. 2-3, pp. 91-100, 1981.
c   2) Osterby, Ole and Zletev, Zahari;
c      "Direct Methods for Sparse Matrices";
c       Springer-Verlag 1983.
c
c
c   INPUT PARAMETERS
c
c   N    - Integer. The size of the square matrix.
c          N > 2 must be specified.
c
c   NN   - Integer. The dimension of integer arrays IA and JA and
c          double precision array AR. Must be at least NE.
c
c   NN  - Integer. The dimension of integer array JA. Must be at least
c          NE.
c
c   IC   - Integer. The sparsity pattern can be changed by means of this
c          parameter.  1 < IC < N   must be specified.
c
c
c
c   OUTPUT PARAMETERS
c
c   NE   - Integer. The number of nonzero elements in the sparse matrix
c          of the type E(N,C). NE = 5*N - 2*IC - 2 .
c
c   AR(NN)  - Real array.
c             Stored entries of the sparse matrix A.
c             NE is the number of nonzeros including a mandatory
c             diagonal entry for each row.
c
c   IA(NN)  - Integer array.(Double precision)
c             Pointers to specify rows for the stored nonzero entries
c             in AR.
c
c   JA(NN) - Integer array.
c             Pointers to specify columns for the stored nonzero entries
c             in AR.
c
c   IERR    - Error parameter is returned as zero on successful
c             execution of the subroutine.
c             Error diagnostics are given by means of positive values
c             of this parameter as follows:
c             IERR = 1    -  N       is out of range.
c             IERR = 2    -  IC      is out of range.
c             IERR = 3    -  NN      is out of range.
c
      double precision ar(nn)
      integer ia(nn), ja(nn), n, ne, ierr
      ierr = 0
c
c  check the input parameters:
c
      if(n.le.2)then
         ierr = 1
         return
      end if
      if(ic.le.1.or.ic.ge.n)then
         ierr = 2
         return
      end if

      ne = 5*n-2*ic-2
      if(nn.lt.ne)then
         ierr = 3
         return
      end if
c
c Begin to generate the nonzero elements as well as the row and column
c pointers:
c
      do 20 i=1,n
      ar(i) = 4.0
      ia(i) = i
      ja(i) = i
20    continue
      ilast = n
      do 30 i=1,n-1
      it = ilast + i
      ar(it) = -1.0
      ia(it) = i+1
      ja(it) = i
30    continue
      ilast = ilast + n - 1
      do 40 i=1,n-1
      it = ilast + i
      ar(it) = -1.0
      ia(it) = i
      ja(it) = i+1
40    continue
      ilast = ilast + n-1
      do 50 i=1,n-ic
      it = ilast + i
      ar(it) = -1.0
      ia(it) = i+ic
      ja(it) = i
50    continue
      ilast = ilast + n-ic
      do 60 I=1,n-ic
      it = ilast + i
      ar(it) = -1.0
      ia(it) = i
      ja(it) = i+ic
60    continue
c      ilast = ilast + n-ic
c      if(ilast.ne.5*n-2*ic-2) then
c      write(*,*)' ilast equal to ', ilast
c      write(*,*)' ILAST, the no. of nonzeros, should = ', 5*n-2*ic-2
c      stop
c      end if
c
      return
      end
      subroutine ellcsr(nrow,coef,jcoef,ncoef,ndiag,a,ja,ia,nzmax,ierr)

c*********************************************************************72
c
cc ELLCSR converts Ellpack/Itpack to Compressed Sparse Row.
c
c this subroutine converts a matrix stored in ellpack-itpack format
c coef-jcoef into the compressed sparse row format. It actually checks
c whether an entry in the input matrix is a nonzero element before
c putting it in the output matrix. The test does not account for small
c values but only for exact zeros.
c
c on entry:
c
c
c nrow       = row dimension of the matrix A.
c coef      = array containing the values of the matrix A in ellpack format.
c jcoef = integer arraycontains the column indices of coef(i,j) in A.
c ncoef = first dimension of arrays coef, and jcoef.
c ndiag = number of active columns in coef, jcoef.
c
c ndiag = on entry the number of columns made available in coef.
c
c on return:
c
c a, ia,
c    ja = matrix in a, ia, ja format where.
c
c nzmax      = size of arrays a and ja. ellcsr will abort if the storage
c         provided in a, ja is not sufficient to store A. See ierr.
c
c ierr       = integer. serves are output error message.
c         ierr = 0 means normal return.
c         ierr = 1 means that there is not enough space in
c         a and ja to store output matrix.
c
      integer ia(nrow+1), ja(*), jcoef(ncoef,1)
      double precision a(*), coef(ncoef,1)
c
c first determine the length of each row of lower-part-of(A)
      ierr = 0
c  check whether sufficient columns are available.
c
c copy elements row by row.
      kpos = 1
      do 6 i=1, nrow
         do 5 k=1,ndiag
            if (coef(i,k) .ne. 0.0) then
               if (kpos .gt. nzmax) then
                  ierr = kpos
                  return
               end if
               a(kpos) = coef(i,k)
               ja(kpos) = jcoef(i,k)
               kpos = kpos+1
          end if
 5       continue
         ia(i+1) = kpos
 6    continue
      return
      end
      subroutine estif3(nel,ske,fe,det,xe,ye,xyke,ierr)

c*********************************************************************72
c
cc ESTIF3 constructs an element stiffness matrix using 3 node triangles.
c
c constructs the element stiffness matrix using 3-node triangular elements
c arguments:
c nel      = element number
c ske      = element stiffness matrix
c fe      = element load vector
c det      = 2*area of the triangle
c xy, ye= coordinates of the three nodal points in an element.
c xyke  = material constants (kxx, kxy, kyx, kyy)
c
      implicit double precision (a-h,o-z)
      dimension ske(3,3), fe(3), xe(3), ye(3), dn(3,2),xyke(2,2)
c
c initialize
c
      area = 0.5*det

      do 200 i=1,3
        fe(i) = 0.0
      do 200 j=1,3
      ske(i,j) = 0.0
 200      continue
c
c get first gradient of shape function
c
      call gradi3(nel,xe,ye,dn,det,ierr)
      if (ierr .ne. 0) return

      do 100 i=1,3
      do 100 j=1,3
      t = 0.0
      do 102 k=1,2
      do 102 l=1,2
 102      t = t+xyke(k,l)*dn(i,k)*dn(j,l)
 100      ske(i,j) = t*area

      return
      end
      subroutine exphes (n,m,dt,eps,u,w,job,z,wkc,beta,errst,hh,ih,
     *                   x, y, indic,ierr)

c*********************************************************************72
c
cc EXPHES computes the Arnoldi basis.
c
c this subroutine computes the Arnoldi basis and the corresponding
c coeffcient vector in the approximation
c
c              w  ::= beta  Vm  ym
c               where ym = exp(- Hm *dt) * e1
c
c to the vector exp(-A dt) w where A is an arbitary matrix and
c w is a given input vector. In case job = 0 the arnoldi basis
c is recomputed. Otherwise the
c code assumes assumes that  u(*) contains an already computed
c arnoldi basis and computes only the y-vector (which is stored in v(*))
c
c en entry:
c
c n      = dimension of matrix
c
c m      = dimension of Krylov subspace (= degree of polynomial
c         approximation to the exponential used. )
c
c dt      = scalar by which to multiply matrix. Can be viewed
c         as a time step. dt must be positive [to be fixed].
c
c eps   = scalar indicating the relative error tolerated for the result.
c         the code will try to compute an answer such that
c         norm2(exactanswer-approximation) / norm2(w) .le. eps
c
c u      = work array of size n*(m+1) to contain the Arnoldi basis
c
c w      = double precision array of length n = input vector to  which exp(-A) is
c         to be applied.
c
c y     = duble precision work array of  size (m+1)
c wkc   = double complex work array of size (m+1)
c
c job      = integer. job indicator. If job .lt.  0 then the Arnoldi
c         basis is recomputed. If job .gt. 0 then it is assumed
c         that the user wants to use a previously computed Krylov
c         subspace but a different dt. Thus the Arnoldi basis and
c         the Hessenberg matrix Hm are not recomputed.
c        In that case the user should not modify the values of beta
c         and the matrices hh and u(n,*) when recalling phipro.
c         job = -1 : recompute basis and get an initial estimate for
c                    time step dt to be used.
c         job = 0  : recompute basis and do not alter dt.
c         job = 1  : do not recompute arnoldi basis.
c
c hh    = work array of size size at least (m+1)*m
c
c ih      = first dimension of hh as declared in the calling program.
c         ih must be .ge. m.
c
c  entries specific to the matrix
c
c diagonal storage is used :
c         a(n,ndiag) is a rectangular array with a(*,k) containing the
c         the diagonal offset by ioff(k) (negative or positive or zero)
c         i.e.,
c        a(i,jdiag) contains the element A(i,i+ioff(jdiag)) in the
c         usual dense storage scheme.
c
c a      = matrix in diagonal storage form
c ioff      = offsets  of diagonals
c ndiag = number of diagonals.
c
c on return:
c
c w2      = resulting vector w2 = exp(-A *dt) * w
c beta  = real equal to the 2-norm of w. Needed if exppro will
c         be recalled with the same Krylov subspace and a different
c         dt.
c errst = rough estimates of the 2-norm of the error.
c hh      = work array of dimension at least (m+1) x m
c
      parameter (ndmx=20)
      implicit double precision (a-h,o-z)
      double precision hh(ih,ih), u(n,*), w(*),z(m+1), errst,x(*), y(*)
      double precision alp0
      double complex alp(ndmx+1), rd(ndmx+1),wkc(ih)
      save
c  use degree 14 chebyshev all the time
      if (indic .ge. 3) goto 60
c
c  input fraction expansion of rational function
c
      ldg= 7
      alp0 =  0.183216998528140087E-11
      alp(1)=( 0.557503973136501826E+02,-0.204295038779771857E+03)
      rd(1)=(-0.562314417475317895E+01, 0.119406921611247440E+01)
      alp(2)=(-0.938666838877006739E+02, 0.912874896775456363E+02)
      rd(2)=(-0.508934679728216110E+01, 0.358882439228376881E+01)
      alp(3)=( 0.469965415550370835E+02,-0.116167609985818103E+02)
      rd(3)=(-0.399337136365302569E+01, 0.600483209099604664E+01)
      alp(4)=(-0.961424200626061065E+01,-0.264195613880262669E+01)
      rd(4)=(-0.226978543095856366E+01, 0.846173881758693369E+01)
      alp(5)=( 0.752722063978321642E+00, 0.670367365566377770E+00)
      rd(5)=( 0.208756929753827868E+00, 0.109912615662209418E+02)
      alp(6)=(-0.188781253158648576E-01,-0.343696176445802414E-01)
      rd(6)=( 0.370327340957595652E+01, 0.136563731924991884E+02)
      alp(7)=( 0.143086431411801849E-03, 0.287221133228814096E-03)
      rd(7)=( 0.889777151877331107E+01, 0.166309842834712071E+02)
c
c     if job .gt. 0 skip arnoldi process:
c
      if (job .gt. 0) goto 2
c
c  normalize vector w and put in first column of u --
c
      beta = sqrt( ddot (n,w,1,w,1))
      if (beta .eq. 0.0) then
         ierr = -1
         indic = 1
         return
      end if

      t = 1.0/beta
      do 25 j=1, n
         u(j,1) = w(j)*t
 25   continue
c  Arnoldi loop
      i1 = 1
 58   i = i1
      i1 = i + 1
      do 59 k=1, n
         x(k) = u(k,i)
 59   continue
      indic = 3
      return
 60   continue
      do 61 k=1, n
         u(k,i1) = y(k)
 61   continue
      i0 =1
c
c switch  for Lanczos version
c
c     i0 = max0(1, i-1)
      call mgsr (n, i0, i1, u, hh(1,i))
      fnorm = fnorm + ddot(i1, hh(1,i),1, hh(1,i),1)
      if (hh(i1,i) .eq. 0.0) m = i
      if  (i .lt. m) goto 58
c  done with arnoldi loop
      rm = dble(m)
      fnorm = sqrt( fnorm / rm )
c  get  : beta*e1 into z
      m1 = m+1
      do 4 i=1,m1
         hh(i,m1) = 0.0
 4    continue
c     compute initial dt when  job .lt.
      if (job .ge. 0) goto 2
c
c     t = eps / beta
c
      t = eps
      do 41 k=1, m-1
         t = t*(1.0 - dble(m-k)/rm )
 41   continue
      t = 2.0*rm* (t**(1.0/rm) )  / fnorm
c      dt = min(dt,t)
      t = min(abs(dt),t)
      dt = sign(t, dt)

 2    continue
      z(1) = beta
      do 3 k=2, m1
         z(k) = 0.0
 3    continue
c
c  get  : exp(H) * beta*e1
c
      call hes(ldg,m1,hh,ih,dt,z,rd,alp,alp0,wkc)
c  error estimate
      errst = abs(z(m1))

      indic = 2
      return
      end
      subroutine exppro(n, m, eps, tn, u, w, x, y, indic, ierr)

c*********************************************************************72
c
cc EXPPRO computes an approximation to the vector
c
c              w :=  exp( - A * tn ) * w
c
c where A is an arbitary matrix and w is a given input vector
c uses a dynamic estimation of internal time advancement (dt)
c
c THIS IS A REVERSE COMMUNICATION IMPLEMENTATION.
c
c USAGE: (see also comments on indic below).
c
c
c      indic = 0
c 1    continue
c      call exppro (n, m, eps, tn, u, w, x, y, indic)
c      if (indic .eq. 1) goto 2 <-- indic .eq.1 means job is finished
c      call matvec(n, x, y)     <--- user's matrix-vec. product
c                                    with x = input vector, and
c                                     y = result = A * x.
c      goto 1
c 2    continue
c      .....
c
c en entry:
c
c n      = dimension of matrix
c
c m      = dimension of Krylov subspace (= degree of polynomial
c         approximation to the exponential used. )
c
c eps   = scalar indicating the relative error tolerated for the result.
c         the code will try to compute an answer such that
c         norm2(exactanswer-approximation) / norm2(w) .le. eps
c
c tn      = scalar by which to multiply matrix. (may be .lt. 0)
c         the code will compute an approximation to exp(- tn * A) w
c         and overwrite the result onto w.
c
c u      = work array of size n*(m+1) (used to hold the Arnoldi basis )
c
c w      = double precision array of length n = input vector to  which exp(-A) is
c         to be applied. this is also an output argument
c
c x, y  = two double precision work vectors of length at least  n each.
c         see indic for usage.
c
c indic = integer used as indicator for the reverse communication.
c         in the first call enter indic = 0. See below for more.
c
c on return:
c
c w     = contains the resulting vector exp(-A * tn ) * w when
c         exppro has finished (see indic)
c
c indic = indicator for the reverse communication protocole.
c       * INDIC .eq. 1  means that exppro has finished and w contains the
c         result.
c       * INDIC .gt. 1 ,  means that exppro has not finished and that
c         it is requesting another matrix vector product before
c         continuing. The user must compute Ax where A is the matrix
c         and x is the vector provided by exppro, and return the
c         result in y. Then exppro must be called again without
c         changing any other argument. typically this must be
c         implemented in a loop with exppro being called as long
c         indic is returned with a value .ne. 1.
c
c ierr  = error indicator.
c         ierr = 1 means phipro was called with indic=1 (not allowed)
c         ierr = -1 means that the input is zero the solution has been
c         unchanged.
c
c NOTES:  im should not exceed 60 in this version  (see ih0 below)
c
c written by Y. Saad -- version feb, 1991.
c
c For reference see fololowing papers :
c (1) E. Gallopoulos and Y. Saad: Efficient solution of parabolic
c     equations by Krylov approximation methods. RIACS technical
c     report 90-14.
c (2) Y.Saad: Analysis of some Krylov subspace approximations to the
c     matrix exponential operator. RIACS Tech report. 90-14
c
      integer n, m, indic, ierr
      double precision eps, tn, u(*), w(n), x(*), y(*)
      parameter (ih0=60)
      double precision hh(ih0,ih0)
      double precision z(ih0),errst, tcur, told, dtl, beta, red
      double complex wkc(ih0)
      save
c
c indic = 3  means  passing through only with result of y= Ax to exphes
c indic = 2  means exphes has finished its job
c indic = 1  means exppro has finished its job (real end)/
c
      ierr = 0
      if (indic .eq. 3) goto 101
      if (indic .eq. 1) then
         ierr = 1
         return
      end if
      ih = ih0
      m  = min(m,ih0)
      tcur = 0.0
      dtl = tn-tcur
      job = -1
c outer loop
 100  continue
c
c  call exponential propagator
c
      told = tcur
 101  continue
c     if (told + dtl .gt. tn) dtl = tn-told
      call  exphes (n,m,dtl,eps,u,w,job,z,wkc,beta,errst,hh,ih,
     *              x,y,indic,ierr)

      if (ierr .ne. 0) return
      if (indic .ge. 3) return
      tcur = told + dtl
c
c     relative error
c
      errst = errst / beta

      if ((errst .le. eps) .and. ( (errst .gt. eps/100.0) .or.
     *     (tcur .eq. tn))) goto 102
c
c     use approximation :  [ new err ] = fact**m  * [cur. error]
c
      red =  (0.5*eps / errst)**(1.0 /dble(m) )
      dtl = dtl*red
      if (abs(told+dtl) .gt. abs(tn) )  dtl = tn-told
      job = 1
      goto 101

 102  continue
      call project(n,m,u,z,w)
      job = 0
      dtl = min(dtl, tn-tcur)
      if (abs(tcur+dtl) .gt. abs(tn)) dtl = tn-tcur
      if (abs(tcur) .lt. abs(tn)) goto 100
      indic = 1

      return
      end
      subroutine expprod(n, m, eps, tn, u, w, x, y, a, ioff, ndiag)

c*********************************************************************72
c
cc EXPPROD computes an approximation to the vector
c
c              w :=  exp( - A * tn ) * w
c
c for matrices stored in diagonal (DIA) format.
c
c this routine constitutes an interface for the routine exppro for
c matrices stored in diagonal (DIA) format.
c
c ARGUMENTS
c
c see exppro for meaning of parameters n, m, eps, tn, u, w, x, y.
c
c a, ioff, and ndiag are the arguments of the matrix:
c
c a(n,ndiag) = a rectangular array with a(*,k) containing the diagonal
c              offset by ioff(k) (negative or positive or zero), i.e.,
c              a(i,jdiag) contains the element A(i,i+ioff(jdiag)) in
c              the usual dense storage scheme.
c
c ioff           = integer array containing the offsets  of the ndiag diagonals
c ndiag      = integer. the number of diagonals.
c
c
      double precision a(*), u(*), w(n), x(n), y(n), tn
      integer ioff(ndiag)

      indic = 0
 101  continue
      call exppro (n, m, eps, tn, u, w, x, y, indic,ierr)
      if (indic .eq. 1) goto 102
c
c     matrix vector-product for diagonal storage --
c
      call oped(n, x, y, a, ioff, ndiag)
      goto 101
 102  continue
      return
      end
      subroutine extbdg (n,a,ja,ia,bdiag,nblk,ao,jao,iao)

c*********************************************************************72
c
cc EXTBDG extracts the main diagonal blocks of a matrix.
c
c this subroutine extracts the main diagonal blocks of a
c matrix stored in compressed sparse row format and puts the result
c into the array bdiag and the remainder in ao,jao,iao.
c
c on entry:
c
c n      = integer. The row dimension of the matrix a.
c a,
c ja,
c ia    = matrix stored in csr format
c nblk  = dimension of each diagonal block. The diagonal blocks are
c         stored in compressed format rowwise,i.e.,we store in
c        succession the i nonzeros of the i-th row after those of
c        row number i-1..
c
c on return:
c
c bdiag = double precision array of size (n x nblk) containing the diagonal
c        blocks of A on return
c ao,
c jao,
c iao   = remainder of the matrix stored in csr format.
c
c           Y. Saad, Sep. 21 1989                                      c
c
      implicit double precision (a-h,o-z)
      double precision bdiag(*),a(*),ao(*)
      integer ia(*),ja(*),jao(*),iao(*)

      m = 1 + (n-1)/nblk
c this version is sequential -- there is a more parallel version
c that goes through the structure twice ....
      ltr =  ((nblk-1)*nblk)/2
      l = m * ltr
      do 1 i=1,l
         bdiag(i) = 0.0
 1    continue
      ko = 0
      kb = 1
      iao(1) = 1

      do 11 jj = 1,m
         j1 = (jj-1)*nblk+1
         j2 =  min0 (n,j1+nblk-1)
         do 12 j=j1,j2
            do 13 i=ia(j),ia(j+1) -1
               k = ja(i)
               if (k .lt. j1) then
                  ko = ko+1
                  ao(ko) = a(i)
                  jao(ko) = k
               else if (k .lt. j) then
c     kb = (jj-1)*ltr+((j-j1)*(j-j1-1))/2+k-j1+1
c     bdiag(kb) = a(i)
                  bdiag(kb+k-j1) = a(i)
               end if
 13         continue
            kb = kb + j-j1
            iao(j+1) = ko+1
 12      continue
 11   continue
      return
      end
      subroutine filter(n,job,drptol,a,ja,ia,b,jb,ib,len,ierr)

c*********************************************************************72
c
cc FILTER copies a matrix, dropping small elements.
c
c     This module removes any elements whose absolute value
c     is small from an input matrix A and puts the resulting
c     matrix in B.  The input parameter job selects a definition
c     of small.
c
c on entry:
c
c  n       = integer. row dimension of matrix
c  job   = integer. used to determine strategy chosen by caller to
c         drop elements from matrix A.
c          job = 1
c              Elements whose absolute value is less than the
c              drop tolerance are removed.
c          job = 2
c              Elements whose absolute value is less than the
c              product of the drop tolerance and the Euclidean
c              norm of the row are removed.
c          job = 3
c              Elements whose absolute value is less that the
c              product of the drop tolerance and the largest
c              element in the row are removed.
c
c drptol = real. drop tolerance used for dropping strategy.
c a
c ja
c ia     = input matrix in compressed sparse format
c len       = integer. the amount of space in arrays a and ja.
c
c on return:
c
c b
c jb
c ib    = resulting matrix in compressed sparse format.
c
c ierr      = integer. containing error message.
c         ierr .eq. 0 indicates normal return
c         ierr .gt. 0 indicates that there is'nt enough
c         space is a and ja to store the resulting matrix.
c         ierr then contains the row number where filter stopped.
c note:
c This module is in place. (b,jb,ib can ne the same as
c       a, ja, ia in which case the result will be overwritten).
c
c           contributed by David Day,  Sep 19, 1989.                   c
c
      double precision a(*),b(*),drptol
      integer ja(*),jb(*),ia(*),ib(*),n,job,len,ierr
      double precision norm,loctol
      integer index,row,k,k1,k2

      index = 1
      do 10 row= 1,n
         k1 = ia(row)
         k2 = ia(row+1) - 1
         ib(row) = index
       goto (100,200,300) job
 100     norm = 1.0
         goto 400
 200     norm = 0.0
         do 22 k = k1,k2
            norm = norm + a(k) * a(k)
 22      continue
         norm = sqrt(norm)
         goto 400
 300     norm = 0.0
         do 23 k = k1,k2
            if( abs(a(k))  .gt. norm) then
               norm = abs(a(k))
            end if
 23      continue
 400     loctol = drptol * norm
       do 30 k = k1,k2
          if( abs(a(k)) .gt. loctol)then
               if (index .gt. len) then
               ierr = row
               return
            end if
            b(index) =  a(k)
            jb(index) = ja(k)
            index = index + 1
         end if
 30   continue
 10   continue
      ib(n+1) = index
      return
      end
      subroutine gen57bl(nx,ny,nz,nfree,na,n,a,ja,ia,iau,stencil)

c*********************************************************************72
c
cc GEN57BL computes the sparse matrix for an elliptic operator.
c
c This subroutine computes the sparse matrix in compressed
c format for the elliptic operator
c
c L u = delx( a . delx u ) + dely ( b . dely u) + delz ( c . delz u ) +
c      delx ( d . u ) + dely (e . u) + delz( f . u ) + g . u
c
c Here u is a vector of nfree components and each of the functions
c a, b, c, d, e, f, g   is an (nfree x nfree) matrix depending of
c the coordinate (x,y,z).
c with Dirichlet Boundary conditions, on a rectangular 1-D,
c 2-D or 3-D grid using centered difference schemes.
c
c The functions a, b, ..., g are known through the
c subroutines  afunbl, bfunbl, ..., gfunbl. (user supplied) .
c
c uses natural ordering, first x direction, then y, then z
c mesh size h is uniform and determined by grid points
c in the x-direction.
c
c parameters:
c
c nx      = number of points in x direction
c ny        = number of points in y direction
c nz        = number of points in z direction
c nfree   = number of degrees of freedom per point
c n        = dimension of matrix (output)
c na        = first dimension of array a as declared in calling
c           program. Must be .ge. nfree**2
c
c a, ja, ia = resulting matrix in  row-sparse block-reduced format
c           a(1:nfree**2, j ) contains a nonzero block.
c           ja(j) contains the column number of (1,1) entry of the block.
c
c iau     = integer*n containing the position of the diagonal element
c           in the a, ja, ia structure
c
c stencil =  work array of size (7,nfree**2), used to store
c            local stencils.
c
c
c     stencil (1:7,*) has the following meaning:
c
c     center point = stencil(1)
c     west point   = stencil(2)
c     east point   = stencil(3)
c     south point  = stencil(4)
c     north point  = stencil(5)
c     front point  = stencil(6)
c     back point   = stencil(7)
c
c
c                           st(5)
c                            |
c                            |
c                            |
c                            |          .st(7)
c                            |     .
c                            | .
c         st(2) ----------- st(1) ---------- st(3)
c                       .    |
c                   .        |
c               .            |
c            st(6)           |
c                            |
c                            |
c                           st(4)
c
c     implicit double precision (a-h,o-z)
      integer ja(*),ia(*),iau(*)
      double precision a(na,1), stencil(7,1)

      h = 1.0/dble(nx+1)
      kx = 1
      ky = nx
      kz = nx*ny
      nfree2 = nfree*nfree
      iedge = 1
      node = 1
      do 100 iz = 1,nz
         do 90 iy = 1,ny
            do 80 ix = 1,nx
               ia(node) = iedge
               call bsten(nx,ny,nz,ix,iy,iz,nfree,stencil,h)
c     west
               if (ix.gt.1) then
                  ja(iedge)=node-kx
                do 4 k=1,nfree2
                 a(iedge,k) = stencil(2,k)
 4              continue
                  iedge=iedge + 1
               end if
c     south
               if (iy.gt.1) then
                  ja(iedge)=node-ky
                do 5 k=1,nfree2
                 a(iedge,k) = stencil(4,k)
 5              continue
                  iedge=iedge + 1
               end if
c     front plane
               if (iz.gt.1) then
                  ja(iedge)=node-kz
                do 6 k=1,nfree2
                 a(iedge,k) = stencil(6,k)
 6              continue
                  iedge=iedge + 1
               end if
c     center node
               ja(iedge) = node
               iau(node) = iedge
               do 7 k=1,nfree2
                  a(iedge,k) = stencil(1,k)
 7             continue
               iedge = iedge + 1
c     -- upper part
c     east
               if (ix.lt.nx) then
                  ja(iedge)=node+kx
                do 8 k=1,nfree2
                 a(iedge,k) = stencil(3,k)
 8              continue
                  iedge=iedge + 1
               end if
c     north
               if (iy.lt.ny) then
                  ja(iedge)=node+ky
                do 9 k=1,nfree2
                 a(iedge,k) = stencil(5,k)
 9              continue
                  iedge=iedge + 1
               end if
c     back plane
               if (iz.lt.nz) then
                  ja(iedge)=node+kz
                do 10 k=1,nfree2
                     a(iedge,k) = stencil(7,k)
 10              continue
                  iedge=iedge + 1
               end if
c  next node
               node=node+1
 80         continue
 90      continue
 100  continue
c     change numbering of nodes so that each ja(k) will contain the
c     actual column number in the original matrix of entry (1,1) of each
c     block (k).
      do 101 k=1,iedge-1
         ja(k) = (ja(k)-1)*nfree+1
 101  continue

      n = (node-1)*nfree
      ia(node)=iedge
      return
      end
      subroutine gen57pt(nx,ny,nz,a,ja,ia,iau,stencil)

c*********************************************************************72
c
cc GEN57PT computes the compressed sparse matrix for an elliptic operator.
c
c L u = delx( a delx u ) + dely ( b dely u) + delz ( c delz u ) +
c       d delx ( u ) + e dely (u) + f delz( u ) + g u
c
c with Dirichlet Boundary conditions, on a rectangular 1-D,
c 2-D or 3-D grid using centered difference schemes.
c
c The functions a, b, ..., g are known through the
c subroutines  afun, bfun, ..., gfun.
c note that to obtain the correct matrix, any function that is not
c needed should be set to zero. For example for two-dimensional
c problems, nz should be set to 1 and the functions cfun and ffun
c should be zero functions.
c
c uses natural ordering, first x direction, then y, then z
c mesh size h is uniform and determined by grid points
c in the x-direction.
c
c parameters:
c
c nx      = number of points in x direction
c ny        = number of points in y direction
c nz        = number of points in z direction
c
c a, ja, ia =  resulting matrix in row-sparse format
c
c iau     = integer*n containing the poisition of the diagonal element
c           in the a, ja, ia structure
c
c stencil =  work array of size 7, used to store local stencils.
c
c
c     stencil [1:7] has the following meaning:
c
c     center point = stencil(1)
c     west point = stencil(2)
c     east point = stencil(3)
c     south point = stencil(4)
c     north point = stencil(5)
c     front point = stencil(6)
c     back point = stencil(7)
c
c
c                           st(5)
c                            |
c                            |
c                            |
c                            |          .st(7)
c                            |     .
c                            | .
c         st(2) ----------- st(1) ---------- st(3)
c                       .    |
c                   .        |
c               .            |
c            st(6)           |
c                            |
c                            |
c                           st(4)
c
c
      integer ja(*),ia(*),iau(*)
      double precision a(*), stencil(*), h

      h = 1.0/dble(nx+1)
      kx = 1
      ky = nx
      kz = nx*ny
      iedge = 1
      node = 1
      do 100 iz = 1,nz
         do 90 iy = 1,ny
            do 80 ix = 1,nx
               ia(node) = iedge
               call getsten(nx,ny,nz,ix,iy,iz,stencil,h)
c     west
               if (ix.gt.1) then
                  ja(iedge)=node-kx
              a(iedge) = stencil(2)
                  iedge=iedge + 1
               end if
c     south
               if (iy.gt.1) then
                  ja(iedge)=node-ky
              a(iedge) = stencil(4)
                  iedge=iedge + 1
               end if
c     front plane
               if (iz.gt.1) then
                  ja(iedge)=node-kz
              a(iedge) = stencil(6)
                  iedge=iedge + 1
               end if
c     center node
               ja(iedge) = node
               iau(node) = iedge
               a(iedge) = stencil(1)
               iedge = iedge + 1
c     -- upper part
c     east
               if (ix.lt.nx) then
                  ja(iedge)=node+kx
              a(iedge) = stencil(3)
                  iedge=iedge + 1
               end if
c     north
               if (iy.lt.ny) then
                  ja(iedge)=node+ky
              a(iedge) = stencil(5)
                  iedge=iedge + 1
               end if
c     back plane
               if (iz.lt.nz) then
                  ja(iedge)=node+kz
                  a(iedge) = stencil(7)
                  iedge=iedge + 1
               end if
c  next node
               node=node+1
 80         continue
 90      continue
 100  continue
      ia(node)=iedge
      return
      end
      subroutine genfea (nx,nelx,node,job,x,y,ijk,nodcode,fs,nint,
     *     a,ja,ia,f,iwk,jwk,ierr,xyk)

c*********************************************************************72
c
cc GENFEA generates finite element matrices for heat conduction problems.
c
c                  - Div ( K(x,y) Grad u ) = f
c                    u = 0 on boundary
c
c (with Dirichlet boundary conditions). The matrix is returned
c assembled in compressed sparse row format. See genfeu for
c matrices in unassembled form. The user must provide the grid,
c (coordinates x, y and connectivity matrix ijk) as well as some
c information on the nodes (nodcode) and the material properties
c (the function K(x,y) above) in the form of a subroutine xyk.
c
c
c
c on entry:
c
c nx          = integer . the number of nodes in the grid .
c nelx          = integer . the number of elements in the grid.
c node      = integer = the number of nodes per element (should be
c             set to three in this version). also the first dimension
c             of ijk
c job          = integer. If job=0, it is assumed that there is no heat
c             source (i.e. fs = 0) and the right hand side
c             produced will therefore be a zero vector.
c             If job = 1 on entry then the contributions from the
c             heat source in each element are taken into account.
c
c x, y      = two double precision arrays containing the coordinates of the nodes.
c
c ijk       =  an integer array containing the connectivity matrix.
c              ijk(i,nel), i=1,2,..node, is the list of the nodes
c              constituting the element nel, ans listed in
c              counter clockwise order.
c
c nodcode   = an integer array containing the boundary information for
c             each node with the following meaning.
c      nodcode(i) = 0 -->  node i is internal
c      nodcode(i) = 1 -->  node i is a boundary but not a corner point
c      nodcode(i) = 2 -->  node i is a corner node. [This node and the
c             corresponmding element are discarded.]
c
c fs          = double precision array of length nelx on entry containing the heat
c             source for each element (job = 1 only)
c
c xyk          = subroutine defining the material properties at each
c            element. Form:
c             call xyk(nel,xyke,x,y,ijk,node) with on return
c             xyke =  material constant matrices.
c            for each element nel, xyke(1,nel),xyke(2,nel)
c             and xyke(3,nel) represent the constants
c             K11, K22, and K12 at that element.
c
c on return
c
c nint          = integer. The number of active (nonboundary) nodes. Also
c             equal to the dimension of the assembled matrix.
c
c a, ja, ia = assembled matrix in compressed sparse row format.
c
c f          = double precision array containing the right hand for the linears
c             system to solve.
c
c ierr          = integer. Error message. If (ierr .ne. 0) on return
c             it means that one of the elements has a negative or zero
c             area probably because of a bad ordering of the nodes
c             (see ijk above). Use the subroutine chkelmt to reorder
c             the nodes properly if necessary.
c iwk, jwk  = two integer work arrays of length nx each.
c
      double precision a(*),x(*),y(*),f(*),fs(*)
      integer ijk(node,*), nodcode(*),ia(*),ja(*),iwk(*),jwk(*)
      external xyk

      ierr = 0
c
c     take into boundary conditions to remove boundary nodes.
c
      call bound (nx,nelx,ijk,nodcode,node,nint,jwk,
     *     x,y,f,iwk)
c
c     assemble the matrix
c
      call assmbo (nx,nelx,node,ijk,nodcode,x,y,
     *     a,ja,ia,f,iwk,jwk,ierr,xyk)
c
c     if applicable (job .eq. 1) get heat source function
c
      indic = 1
      if (job .eq. 1)
     *     call hsourc (indic,nx,nelx,node,x,y,ijk,fs,f)
c
c     call diric for Dirichlet conditions
c
      call diric(nx,nint,a,ja,ia,f)
c     done
      return
      end
      subroutine genfeu (nx,nelx,node,job,x,y,ijk,nodcode,fs,nint,
     *     a,na,f,iwk,jwk,ierr,xyk)

c*********************************************************************72
c
cc GENFEU generates finite element matrices for heat conduction problems.
c
c                  - Div ( K(x,y) Grad u ) = f
c                    u = 0 on boundary
c
c (with Dirichlet boundary conditions). The matrix is returned
c in unassembled form. The user must provide the grid,
c (coordinates x, y and connectivity matrix ijk) as well as some
c information on the nodes (nodcode) and the material properties
c (the function K(x,y) above) in the form of a subroutine xyk.
c
c
c
c on entry:
c
c
c nx          = integer . the number of nodes in the grid .
c nelx          = integer . the number of elements in the grid.
c node      = integer = the number of nodes per element (should be
c             set to three in this version). also the first dimension
c             of ijk
c job          = integer. If job=0, it is assumed that there is no heat
c             source (i.e. fs = 0) and the right hand side
c             produced will therefore be a zero vector.
c             If job = 1 on entry then the contributions from the
c             heat source in each element are taken into account.
c
c na          = integer. The first dimension of the array a.
c             a is declared as an array of dimension a(na,node,node).
c
c x, y      = two double precision arrays containing the coordinates of the nodes.
c
c ijk       =  an integer array containing the connectivity matrix.
c              ijk(i,nel), i=1,2,..node, is the list of the nodes
c              constituting the element nel, ans listed in
c              counter clockwise order.
c
c xyk          = subroutine defining the material properties at each
c            element. Form:
c             call xyk(nel,xyke,x,y,ijk,node) with on return
c             xyke =  material constant matrices.
c            for each element nel, xyke(1,nel),xyke(2,nel)
c             and xyke(3,nel) represent the constants
c             K11, K22, and K12 at that element.
c
c nodcode   = an integer array containing the boundary information for
c             each node with the following meaning.
c      nodcode(i) = 0 -->  node i is internal
c      nodcode(i) = 1 -->  node i is a boundary but not a corner point
c      nodcode(i) = 2 -->  node i is a corner node. [This node and the
c             corresponmding element are discarded.]
c
c fs          = double precision array of length nelx on entry containing the heat
c             source for each element (job = 1 only)
c
c on return
c
c nint          = integer. The number of active (nonboundary) nodes. Also
c             equal to the dimension of the assembled matrix.
c
c a         = matrix in unassembled form. a(nel,*,*) contains the
c             element matrix for element nel.
c
c f          = double precision array containing the right hand for the linears
c             system to solve, in assembled form.
c
c ierr          = integer. Error message. If (ierr .ne. 0) on return
c             it means that one of the elements has a negative or zero
c             area probably because of a bad ordering of the nodes
c             (see ijk above). Use the subroutine chkelmt to reorder
c             the nodes properly if necessary.
c iwk, jwk  = two integer work arrays of length nx each.
c
      double precision a(na,node,node),x(*),y(*),f(*), fs(*)
      integer ijk(node,*), nodcode(*),iwk(*),jwk(*)
      external xyk

      ierr = 0
c
c     take boundary conditions into account to move boundary nodes to
c     the end.
c
      call bound (nx,nelx,ijk,nodcode,node,nint,jwk,
     *     x,y,f,iwk)
c
c     assemble the matrix
c
      call unassbl (a,na,f,nx,nelx,ijk,nodcode,
     *     node,x,y,ierr,xyk)
c
c     if applicable (job .eq. 1) get heat source function
c
      indic = 0

      if (job .eq. 1) then
          call hsourc (indic,nx,nelx,node,x,y,ijk,fs,f)
      end if

      return
      end
      subroutine getbwd(n,a,ja,ia,ml,mu)

c*********************************************************************72
c
cc GETBWD gets the bandwidth of lower part and upper part of A.
c
c  Discussion:
c
c    This routine does not assume that A is sorted.
c
c on entry:
c
c n      = integer = the row dimension of the matrix
c a, ja,
c    ia = matrix in compressed sparse row format.
c
c on return:
c
c ml      = integer. The bandwidth of the strict lower part of A
c mu      = integer. The bandwidth of the strict upper part of A
c
c Notes:
c ===== ml and mu are allowed to be negative or return. This may be
c       useful since it will tell us whether a band is confined
c       in the strict  upper/lower triangular part.
c       indeed the definitions of ml and mu are
c
c       ml = max ( (i-j)  s.t. a(i,j) .ne. 0  )
c       mu = max ( (j-i)  s.t. a(i,j) .ne. 0  )
c
c Y. Saad, Sep. 21 1989                                                c
c
      double precision a(*)
      integer ja(*),ia(n+1),ml,mu,ldist,i,k
      ml = - n
      mu = - n
      do 3 i=1,n
         do 31 k=ia(i),ia(i+1)-1
            ldist = i-ja(k)
            ml = max(ml,ldist)
            mu = max(mu,-ldist)
 31      continue
 3    continue
      return
      end
      subroutine getdia (nrow,ncol,job,a,ja,ia,len,diag,idiag,ioff)

c*********************************************************************72
c
cc GETDIA extracts a given diagonal from a matrix stored in CSR format.
c
c this subroutine extracts a given diagonal from a matrix stored in CSR
c format. The output matrix may be transformed with the diagonal removed
c from it if desired (as indicated by job.)
c
c Our definition of a diagonal of matrix is a vector of length nrow
c (always) which contains the elements in rows 1 to nrow of
c the matrix that are contained in the diagonal offset by ioff
c with respect to the main diagonal. If the diagonal element
c falls outside the matrix then it is defined as a zero entry.
c Thus the proper definition of diag(*) with offset ioff is
c
c     diag(k) = a(k,ioff+k) k=1,2,...,nrow
c     with elements falling outside the matrix being defined as zero.
c
c on entry:
c
c
c nrow      = integer. The row dimension of the matrix A.
c ncol      = integer. The column dimension of the matrix A.
c job   = integer. Job indicator.  If job = 0 then
c         the matrix a, ja, ia, is not altered on return.
c         if job.ne.1  then getdia will remove the entries
c         collected in diag from the original matrix.
c         This is done in place.
c
c a,ja,
c    ia = matrix stored in compressed sparse row a,ja,ia,format
c ioff  = integer,containing the offset of the wanted diagonal
c        the diagonal extracted is the one corresponding to the
c        entries a(i,j) with j-i = ioff.
c        thus ioff = 0 means the main diagonal
c
c on return:
c
c len   = number of nonzero elements found in diag.
c         (len .le. min(nrow,ncol-ioff)-max(1,1-ioff) + 1 )
c
c diag  = double precision array of length nrow containing the wanted diagonal.
c        diag contains the diagonal (a(i,j),j-i = ioff ) as defined
c         above.
c
c idiag = integer array of  length len, containing the poisitions
c         in the original arrays a and ja of the diagonal elements
c         collected in diag. A zero entry in idiag(i) means that
c         there was no entry found in row i belonging to the diagonal.
c
c a, ja,
c    ia = if job .ne. 0 the matrix is unchanged. otherwise the nonzero
c         diagonal entries collected in diag are removed from the
c         matrix. the structure is modified since the diagonal elements
c        are removed from a,ja,ia. Thus, the  returned matrix will
c         have len fewer elements if the diagonal is full.
c
c           Y. Saad, Sep. 21 1989 - modified and tested May 9, 1990.
c
      implicit double precision (a-h,o-z)
      double precision diag(*),a(*)
      integer nrow, ncol, job, len, ia(*), ja(*), idiag(*)

      integer istart, max, iend, i, kold, k, kdiag, ko

      istart = max(0,-ioff)
      iend = min0(nrow,ncol-ioff)
      len = 0
      do 1 i=1,nrow
         idiag(i) = 0
       diag(i) = 0.0
 1    continue
c
c     extract  diagonal elements
c
      do 6 i=istart+1, iend
         do 51 k= ia(i),ia(i+1) -1
            if (ja(k)-i .eq. ioff) then
               diag(i)= a(k)
               idiag(i) = k
               len = len+1
               goto 6
            end if
 51      continue
 6    continue
      if (job .eq. 0 .or. len .eq.0) return
c
c  rewind the structure
c
      ko = 0
      do  7 i=istart+1,iend
         kold = ko
         kdiag = idiag(i)
         if (kdiag .eq. 0) goto 7
         do 71 k= ia(i), ia(i+1)-1
            if (ja(k) .eq. kdiag) goto 71
            ko = ko+1
            a(ko) = a(k)
            ja(ko) = ja(k)
 71      continue
         ia(i) = kold+1
 7    continue
c
c redefine ia(nrow+1)
c
      ia(nrow+1) = ko+1
      return
      end
      function getelm (i,j,a,ja,ia,iadd,sorted)

c*********************************************************************72
c
cc GETELM returns the element A(I,J) of a CSR matrix A.
c
c     Purpose:
c
c     This function returns the element a(i,j) of a matrix A,
c     for any pair (i,j).  The matrix is assumed to be stored
c     in Compressed Sparse Row (CSR) format. GETELM Performs a
c     binary search in the case where it is known that the elements
c     are sorted so that the column indices are in increasing order.
c     Also returns (in iadd) the address of the element a(i,j) in
c     arrays A and JA when the search is successsful (zero if not).
c
c     First contributed by Noel Nachtigal (MIT).
c     Recoded Jan. 20, 1991, by Y. Saad [In particular
c     added handling of the non-sorted case + the IADD output]
c
c     Parameters:
c
c On entry:
c
c     I      = the row index of the element sought (input).
c     J      = the column index of the element sought (input).
c     A      = the matrix A in compressed sparse row format (input).
c     JA     = the array of column indices (input).
c     IA     = the array of pointers to the rows' data (input).
c     SORTED = logical indicating whether the matrix is knonw to
c              have its column indices sorted in increasing order
c              (sorted=.true.) or not (sorted=.false.).
c              (input).
c On return:
c
c     GETELM = value of a(i,j).
c     IADD   = address of element a(i,j) in arrays a, ja if found,
c              zero if not found. (output)
c
c     Note: the inputs I and J are not checked for validity.
c
c     Noel M. Nachtigal October 28, 1990 -- Youcef Saad Jan 20, 1991.
c
      INTEGER I, IA(*), IADD, J, JA(*)
      double precision A(*)
      double precision getelm
      LOGICAL SORTED
      INTEGER IBEG, IEND, IMID, K
c
c     Initialization
c
      IADD = 0
      GETELM = 0.0
      IBEG = IA(I)
      IEND = IA(I+1)-1
c
c  case where matrix is not necessarily sorted
c
      IF (.NOT. SORTED) THEN
c
c scan the row - exit as soon as a(i,j) is found
c
         DO 5  K=IBEG, IEND
            IF (JA(K) .EQ.  J) THEN
               IADD = K
               GOTO 20
            end if
 5       CONTINUE
c
c     end unsorted case. begin sorted case
c
      ELSE
c
c     begin binary search.   Compute the middle index.
c
 10      IMID = ( IBEG + IEND ) / 2
c
c     test if  found
c
         IF (JA(IMID).EQ.J) THEN
            IADD = IMID
            GOTO 20
         end if
         IF (IBEG .GE. IEND) GOTO 20
c
c     else     Update the interval bounds.
c
         IF (JA(IMID).GT.J) THEN
            IEND = IMID -1
         ELSE
            IBEG = IMID +1
         end if
         GOTO 10
c
c     end both cases
c
      end if

 20   IF (IADD .NE. 0) GETELM = A(IADD)

      RETURN
      END
      subroutine getl (n,a,ja,ia,ao,jao,iao)

c*********************************************************************72
c
cc GETL extracts the lower triangular part of a matrix.
c
c This routine extracts the lower triangular part of a matrix
c and writes the result ao, jao, iao. The routine is in place in
c that ao, jao, iao can be the same as a, ja, ia if desired.
c
c on input:
c
c n     = dimension of the matrix a.
c a, ja,
c    ia = matrix stored in compressed sparse row format.
c On return:
c ao, jao,
c    iao = lower triangular matrix (lower part of a)
c      stored in a, ja, ia, format
c note: the diagonal element is the last element in each row.
c i.e. in  a(ia(i+1)-1 )
c ao, jao, iao may be the same as a, ja, ia on entry -- in which case
c getl will overwrite the result on a, ja, ia.
c
      integer n, ia(*), ja(*), iao(*), jao(*)
      double precision a(*), ao(*)
      double precision t
      integer ko, kold, kdiag, k, i
c
c  inititialize ko (pointer for output matrix)
c
      ko = 0
      do  7 i=1, n
         kold = ko
         kdiag = 0
         do 71 k = ia(i), ia(i+1) -1
            if (ja(k)  .gt. i) goto 71
            ko = ko+1
            ao(ko) = a(k)
            jao(ko) = ja(k)
            if (ja(k)  .eq. i) kdiag = ko
 71      continue
         if (kdiag .eq. 0 .or. kdiag .eq. ko) goto 72
c
c  exchange
c
         t = ao(kdiag)
         ao(kdiag) = ao(ko)
         ao(ko) = t

         k = jao(kdiag)
         jao(kdiag) = jao(ko)
         jao(ko) = k
 72      iao(i) = kold+1
 7    continue
c     redefine iao(n+1)
      iao(n+1) = ko+1
      return
      end
      subroutine getsten(nx,ny,nz,kx,ky,kz,stencil,h)

c*********************************************************************72
c
cc GETSTEN calculates the stencil for centered elliptic discretization.
c
c     This subroutine calculates the correct stencil values for
c     centered difference discretization of the elliptic operator
c
c L u = delx( a delx u ) + dely ( b dely u) + delz ( c delz u ) +
c      delx ( d u ) + dely (e u) + delz( f u ) + g u
c
c   For 2-D problems the discretization formula that is used is:
c
c h**2 * Lu == a(i+1/2,j)*{u(i+1,j) - u(i,j)} +
c             a(i-1/2,j)*{u(i-1,j) - u(i,j)} +
c              b(i,j+1/2)*{u(i,j+1) - u(i,j)} +
c              b(i,j-1/2)*{u(i,j-1) - u(i,j)} +
c              (h/2)*d(i,j)*{u(i+1,j) - u(i-1,j)} +
c              (h/2)*e(i,j)*{u(i,j+1) - u(i,j-1)} +
c              (h/2)*e(i,j)*{u(i,j+1) - u(i,j-1)} +
c               (h**2)*g(i,j)*u(i,j)
c
      double precision stencil(*), h, hhalf,cntr, afun, bfun, cfun, dfun,
     *      efun, ffun, gfun, x, y, z, coeff

      do k=1,7
         stencil(k) = 0.0
      end do

      hhalf = h*0.5
      x = h*dble(kx)
      y = h*dble(ky)
      z = h*dble(kz)
      cntr = 0.0
c     differentiation wrt x:
      coeff = afun(x+hhalf,y,z)
      stencil(3) = stencil(3) + coeff
      cntr = cntr + coeff

      coeff = afun(x-hhalf,y,z)
      stencil(2) = stencil(2) + coeff
      cntr = cntr + coeff

      coeff = dfun(x,y,z)*hhalf
      stencil(3) = stencil(3) + coeff
      stencil(2) = stencil(2) - coeff
      if (ny .le. 1) goto 99
c
c     differentiation wrt y:
c
      coeff = bfun(x,y+hhalf,z)
      stencil(5) = stencil(5) + coeff
      cntr = cntr + coeff

      coeff = bfun(x,y-hhalf,z)
      stencil(4) = stencil(4) + coeff
      cntr = cntr + coeff

      coeff = efun(x,y,z)*hhalf
      stencil(5) = stencil(5) + coeff
      stencil(4) = stencil(4) - coeff
      if (nz .le. 1) goto 99
c
c differentiation wrt z:
c
      coeff = cfun(x,y,z+hhalf)
      stencil(7) = stencil(7) + coeff
      cntr = cntr + coeff

      coeff = cfun(x,y,z-hhalf)
      stencil(6) = stencil(6) + coeff
      cntr = cntr + coeff

      coeff = ffun(x,y,z)*hhalf
      stencil(7) = stencil(7) + coeff
      stencil(6) = stencil(6) - coeff
c
c discretization of  product by g:
c
 99   coeff = gfun(x,y,z)
      stencil(1) = h*h*coeff - cntr
      return
      end
      subroutine getu (n,a,ja,ia,ao,jao,iao)

c*********************************************************************72
c
cc GETU extracts the upper triangular part of a matrix.
c
c This routine extracts the upper triangular part of a matrix
c and writes the result ao, jao, iao. The routine is in place in
c that ao, jao, iao can be the same as a, ja, ia if desired.
c
c on input:
c
c n     = dimension of the matrix a.
c a, ja,
c    ia = matrix stored in a, ja, ia, format
c On return:
c ao, jao,
c    iao = upper triangular matrix (upper part of a)
c      stored in compressed sparse row format
c note: the diagonal element is the last element in each row.
c i.e. in  a(ia(i+1)-1 )
c ao, jao, iao may be the same as a, ja, ia on entry -- in which case
c getu will overwrite the result on a, ja, ia.
c
      integer n, ia(*), ja(*), iao(*), jao(*)
      double precision a(*), ao(*)
      double precision t
      integer ko, k, i, kdiag, kfirst

      ko = 0
      do  7 i=1, n
         kfirst = ko+1
         kdiag = 0
         do 71 k = ia(i), ia(i+1) -1
            if (ja(k)  .lt. i) goto 71
            ko = ko+1
            ao(ko) = a(k)
            jao(ko) = ja(k)
            if (ja(k)  .eq. i) kdiag = ko
 71      continue
         if (kdiag .eq. 0 .or. kdiag .eq. kfirst) goto 72
c     exchange
         t = ao(kdiag)
         ao(kdiag) = ao(kfirst)
         ao(kfirst) = t

         k = jao(kdiag)
         jao(kdiag) = jao(kfirst)
         jao(kfirst) = k
 72      iao(i) = kfirst
 7    continue
c     redefine iao(n+1)
      iao(n+1) = ko+1
      return
      end
      subroutine gradi3(nel, xe, ye, dn, det,ierr)

c*********************************************************************72
c
cc GRADI3 constructs the first derivative of the shape functions.
c
c arguments:
c nel      = element nuumber
c xy, ye= coordinates of the three nodal points in an element.
c dn      = gradients (1-st derivatives) of the shape functions.
c area      = area of the triangle
c
c
      PARAMETER (TOL=1.0E-17)
c
      dimension xe(3), ye(3), dn(3,2)
c
c compute area
c
      ierr = 0
      if (det .le. TOL) goto 100
c
      dn(1,1) = (ye(2)-ye(3))/det
      dn(2,1) = (ye(3)-ye(1))/det
      dn(3,1) = (ye(1)-ye(2))/det
      dn(1,2) = (xe(3)-xe(2))/det
      dn(2,2) = (xe(1)-xe(3))/det
      dn(3,2) = (xe(2)-xe(1))/det
c
      return
c
 100      continue
      ierr = 3
c       write(iout,*) ' ** Error-negative area encountered at elmt: '
c      write(iout,*) nel,(xe(i),ye(i),i=1,3)
      return
      end
      subroutine hes (ndg,m,hh,ih,dt,y,root,coef,coef0,w2)

c*********************************************************************72
c
cc HES computes exp ( H dt) * y  where H = Hessenberg matrix (hh)
c
c ndg      = number of poles as determined by getrat
c m     = dimension of hessenberg matrix
c hh      = hessenberg matrix (real)
c ih      = first dimenbsion of hh
c dt      = scaling factor used for hh (see (1))
c y      = double precision vector. on return exp(H dt ) y is computed
c         and overwritten on y.
c root  = poles of the rational approximation to exp as
c         computed by getrat
c coef,
c coef0 = coefficients of partial fraction expansion
c
c  exp(t) ~ coef0 +  sum     Real [   coef(i) / (t - root(i)  ]
c                  i=1,ndg
c
c valid for double precision real t.
c coef0 is double precision, coef(*) is a double complex array.
c
      parameter (mmax=70)
      implicit double precision (a-h,o-z)
      double precision  hh(ih,*), y(*)
      double complex hloc(mmax+1,mmax), coef(*),root(*),w2(*),t,zpiv
      double precision yloc(mmax)
c
c  loop associated with the poles.
c
      do 10 j=1,m
         yloc(j) = y(j)
         y(j) = y(j)*coef0
 10   continue

      do 8 ii = 1, ndg
c
c  copy Hessenberg matrix into temporary
c
         do 2 j=1, m
            do 1 i=1, j+1
               hloc(i,j) = CMPLX( dt*hh(i,j) )
 1          continue
            hloc(j,j) = hloc(j,j) - root(ii)
            w2(j)     = CMPLX(yloc(j))
 2       continue
c
c  forward solve
c
         do 4 i=2,m
            zpiv  = hloc(i,i-1) / hloc(i-1,i-1)
            do 3 j=i,m
               hloc(i,j) = hloc(i,j) - zpiv*hloc(i-1,j)
 3          continue
            w2(i)     = w2(i) - zpiv*w2(i-1)
 4       continue
c
c     backward solve
c
         do 6 i=m,1,-1
            t=w2(i)
            do 5 j=i+1,m
               t = t-hloc(i,j)*w2(j)
 5          continue
            w2(i) = t/hloc(i,i)
 6       continue
c
c     accumulate result in y.
c
         do 7 i=1,m
            y(i) = y(i) + coef(ii) * w2(i)
 7       continue
 8    continue
      return
      end
      subroutine hsourc (indic,nx,nelx,node,x,y,ijk,fs,f)

c*********************************************************************72
c
cc HSOURC assembles the load vector F from element contributions in FS.
c
c generates the load vector f in assembled form from the
c the element contributions fs.
c indic = indicates if f is to be assembled (1) or not (zero)
c note: f(*) not initilazed. because might use values from boundary
c conditions.
c
      implicit double precision (a-h,o-z)
        double precision x(*),y(*),fs(*),f(*),xe(3),ye(3),det,areao3
      integer ijk(node,*)

      jnod = 0
      do 130 nel = 1,nelx
c
c get coordinates of nodal points
c
      do 104 i=1, node
      j = ijk(i,nel)
      xe(i) = x(j)
      ye(i) = y(j)
 104      continue
c
c compute determinant
c
      det=xe(2)*(ye(3)-ye(1))+xe(3)*(ye(1)-ye(2))+xe(1)*(ye(2)-ye(3))
c area3 = area/3
      areao3 = det/6.0
c
c contributions to nodes in the element
c
      if (indic .eq. 0) then
         do 115 ka=1,node
         jnod = jnod+1
         f(jnod) = fs(nel)*areao3
 115      continue
      else
      do 120 ka=1, node
            ii = ijk(ka,nel)
          f(ii) = f(ii) + fs(nel)*areao3
 120      continue
      end if

 130      continue
      return
      end
      subroutine ilu0(n, a, ja, ia, alu, jlu, ju, iw, ierr)

c*********************************************************************72
c
cc ILU0 is an ILU(0) preconditioner.
c
c Note that this has been coded in such a way that it can be used
c with pgmres. Normally, since the data structure of a, ja, ia is
c the same as that of a, ja, ia, savings can be made. In fact with
c some definitions (not correct for general sparse matrices) all we
c need in addition to a, ja, ia is an additional diagonal.
c Ilu0 is not recommended for serious problems. It is only provided
c here for comparison purposes.
c
c on entry:
c
c n       = dimension of matrix
c a, ja,
c ia      = original matrix in compressed sparse row storage.
c
c on return:
c
c alu,jlu = matrix stored in Modified Sparse Row (MSR) format containing
c           the L and U factors together. The diagonal (stored in
c           alu(1:n) ) is inverted. Each i-th row of the alu,jlu matrix
c           contains the i-th row of L (excluding the diagonal entry=1)
c           followed by the i-th row of U.
c
c ju        = pointer to the diagonal elements in alu, jlu.
c
c ierr        = integer indicating error code on return
c           ierr = 0 --> normal return
c           ierr = k --> code encountered a zero pivot at step k.
c work arrays:
c
c iw          = integer work array of length n.
c
c IMPORTANT
c
c it is assumed that the the elements in the input matrix are stored
c    in such a way that in each row the lower part comes first and
c    then the upper part. To get the correct ILU factorization, it is
c    also necessary to have the elements of L sorted by increasing
c    column number. It may therefore be necessary to sort the
c    elements of a, ja, ia prior to calling ilu0. This can be
c    achieved by transposing the matrix twice using csrcsc.
c
      implicit double precision (a-h,o-z)
      double precision a(*), alu(*)
        integer ja(*), ia(*), ju(*), jlu(*), iw(*)

        ju0 = n+2
        jlu(1) = ju0
c
c initialize work vector to zero's
c
      do 31 i=1, n
           iw(i) = 0
 31     continue
c
c main loop
c
      do 500 ii = 1, n
           js = ju0
c
c generating row number ii of L and U.
c
           do 100 j=ia(ii),ia(ii+1)-1
c
c     copy row ii of a, ja, ia into row ii of alu, jlu (L/U) matrix.
c
              jcol = ja(j)
              if (jcol .eq. ii) then
                 alu(ii) = a(j)
                 iw(jcol) = ii
                 ju(ii)  = ju0
              else
                 alu(ju0) = a(j)
                 jlu(ju0) = ja(j)
                 iw(jcol) = ju0
                 ju0 = ju0+1
              end if
 100       continue
           jlu(ii+1) = ju0
           jf = ju0-1
           jm = ju(ii)-1
c
c exit if diagonal element is reached.
c
           do 150 j=js, jm
              jrow = jlu(j)
              tl = alu(j)*alu(jrow)
              alu(j) = tl
c
c perform  linear combination
c
              do 140 jj = ju(jrow), jlu(jrow+1)-1
                 jw = iw(jlu(jj))
                 if (jw .ne. 0) alu(jw) = alu(jw) - tl*alu(jj)
 140          continue
 150       continue
c
c     invert  and store diagonal element.
c
           if (alu(ii) .eq. 0.0) goto 600
           alu(ii) = 1.0/alu(ii)
c
c     reset pointer iw to zero
c
           iw(ii) = 0
           do 201 i = js, jf
 201          iw(jlu(i)) = 0
 500       continue
           ierr = 0
           return
c
c     zero pivot :
c
 600       ierr = ii

           return
           end
       subroutine ilut (n,a,ja,ia,lfil,tol,alu,jlu,ju,iwk,
     *                  wu,wl,jr,jwl,jwu,ierr)

c*********************************************************************72
c
cc ILUT is an ILUT preconditioner.
c
c      incomplete LU factorization with dual truncation mechanism
c      VERSION 2 : sorting  done for both L and U.
c
c
c  coded by Youcef Saad May, 5, 1990.
c  Dual drop-off strategy works as follows.
c
c     1) Theresholding in L and U as set by tol. Any element whose size
c        is less than some tolerance (relative to the norm of current
c        row in u) is dropped.
c
c     2) Keeping only the largest lenl0+lfil elements in L and the
c        largest lenu0+lfil elements in U, where lenl0=initial number
c        of nonzero elements in a given row of lower part of A
c        and lenlu0 is similarly defined.
c
c Flexibility: one can use tol=0 to get a strategy based on keeping the
c largest elements in each row of L and U. Taking tol .ne. 0 but lfil=n
c will give the usual threshold strategy (however, fill-in is then
c impredictible).
c
c*
c PARAMETERS
c
c on entry:
c
c n       = integer. The dimension of the matrix A.
c
c a,ja,ia = matrix stored in Compressed Sparse Row format.
c
c lfil    = integer. The fill-in parameter. Each row of L and
c           each row of U will have a maximum of lfil elements
c           in addition to the original number of nonzero elements.
c           Thus storage can be determined beforehand.
c           lfil must be .ge. 0.
c
c iwk     = integer. The minimum length of arrays alu and jlu
c
c On return:
c
c
c alu,jlu = matrix stored in Modified Sparse Row (MSR) format containing
c           the L and U factors together. The diagonal (stored in
c           alu(1:n) ) is inverted. Each i-th row of the alu,jlu matrix
c           contains the i-th row of L (excluding the diagonal entry=1)
c           followed by the i-th row of U.
c
c ju      = integer array of length n containing the pointers to
c           the beginning of each row of U in the matrix alu,jlu.
c
c ierr    = integer. Error message with the following meaning.
c           ierr  = 0    --> successful return.
c           ierr .gt. 0  --> zero pivot encountered at step number ierr.
c           ierr  = -1   --> Error. input matrix may be wrong.
c                            (The elimination process has generated a
c                            row in L or U whose length is .gt.  n.)
c           ierr  = -2   --> The matrix L overflows the array al.
c           ierr  = -3   --> The matrix U overflows the array alu.
c           ierr  = -4   --> Illegal value for lfil.
c           ierr  = -5   --> zero pivot encountered.
c
c work arrays:
c
c jr,jwu,jwl         = integer work arrays of length n.
c wu, wl          = double precision work arrays of length n+1, and n resp.
c
c Notes:
c
c A must have all nonzero diagonal elements.
c
       implicit double precision (a-h,o-z)
      INTEGER N
       double precision a(*), alu(*), wu(n), wl(n), tol
       integer ja(*),ia(n+1),jlu(*),ju(n),jr(n), jwu(n),
     *      jwl(n), lfil, iwk, ierr

        if (lfil .lt. 0) goto 998
c
c initialize ju0 (points to next element to be added to alu,jlu)
c and pointer.
c
      ju0 = n+2
      jlu(1) = ju0
c
c  integer double pointer array.
c
      do 1 j=1, n
            jr(j)  = 0
 1           continue
c
c  beginning of main loop.
c
      do 500 ii = 1, n
           j1 = ia(ii)
           j2 = ia(ii+1) - 1
           lenu = 0
           lenl = 0
           tnorm = 0.0
           do 501 k=j1,j2
              tnorm = tnorm+abs(a(k))
 501          continue
              tnorm = tnorm/dble(j2-j1+1)
c
c  unpack L-part and U-part of row of A in arrays wl, wu --
c
      do 170  j = j1, j2
           k = ja(j)
           t = a(j)
           if (abs(t) .lt. tol*tnorm) goto 170
           if (k .lt. ii) then
              lenl = lenl+1
              jwl(lenl) = k
              wl(lenl) = t
              jr(k) = lenl
           else
              lenu = lenu+1
              jwu(lenu) = k
              wu(lenu) = t
              jr(k) = lenu
           end if
 170      continue
c      tnorm = 0.0
c      do 171 k=j1,j2
c           tnorm = tnorm + abs(a(k))
c 171      continue
c
c        tnorm = tnorm/dble(j2-j1+1)
        lenl0 = lenl
        lenu0 = lenu
        jj = 0
        nl = 0
c
c  eliminate previous rows
c
 150    jj = jj+1
        if (jj .gt. lenl) goto 160
c
c in order to do the elimination in the correct order we need to
c exchange the current row number with the one that has
c smallest column number, among jj,jj+1,...,lenl.
c
        jrow = jwl(jj)
        k = jj
c
c determine smallest column index
c
        do 151 j=jj+1,lenl
           if (jwl(j) .lt. jrow) then
              jrow = jwl(j)
              k = j
           end if
 151    continue
c
c exchange in jwl
c
        j = jwl(jj)
        jwl(jj) = jrow
        jwl(k) = j
c
c exchange in jr
c
        jr(jrow) = jj
        jr(j) = k
c
c exchange in wl
c
        s = wl(k)
        wl(k) = wl(jj)
        wl(jj) = s
c
        if (jrow .ge. ii) goto 160
c  get the multiplier for row to be eliminated: jrow
        fact = wl(jj)*alu(jrow)
        jr(jrow) = 0
        if (abs(fact)*wu(n+2-jrow) .le. tol*tnorm) goto 150
c
c  combine current row and row jrow
c
        do 203 k = ju(jrow), jlu(jrow+1)-1
           s = fact*alu(k)
           j = jlu(k)
           jpos = jr(j)
c
c if fill-in element and small disregard:
c
           if (abs(s) .lt. tol*tnorm .and. jpos .eq. 0) goto 203
           if (j .ge. ii) then
c
c     dealing with upper part.
c
              if (jpos .eq. 0) then
c     this is a fill-in element
                 lenu = lenu+1
                 if (lenu .gt. n) goto 995
                 jwu(lenu) = j
                 jr(j) = lenu
                 wu(lenu) = - s
              else
c     no fill-in element --
                 wu(jpos) = wu(jpos) - s
              end if
           else
c
c     dealing with lower part.
c
              if (jpos .eq. 0) then
c     this is a fill-in element
                 lenl = lenl+1
                 if (lenl .gt. n) goto 995
                 jwl(lenl) = j
                 jr(j) = lenl
                 wl(lenl) = - s
              else
c     no fill-in element --
                 wl(jpos) = wl(jpos) - s
              end if
           end if
 203      continue
        nl = nl+1
        wl(nl) = fact
        jwl(nl)  = jrow
      goto 150
c
c  update l-matrix
c
 160    len = min0(nl,lenl0+lfil)
        call bsort2 (wl,jwl,nl,len)
c
        do 204 k=1, len
           if (ju0 .gt. iwk) goto 996
           alu(ju0) =  wl(k)
           jlu(ju0) =  jwl(k)
           ju0 = ju0+1
 204    continue
c
c  save pointer to beginning of row ii of U
c
        ju(ii) = ju0
c
c  reset double-pointer jr to zero (L-part - except first
c  jj-1 elements which have already been reset)
c
      do 306 k= jj, lenl
              jr(jwl(k)) = 0
 306      continue
c
c
c be sure that the diagonal element is first in w, jw
c
        idiag = 0
        idiag = jr(ii)
        if (idiag .eq. 0) goto 900

        if (idiag .ne. 1) then
           s = wu(1)
           wu(j) = wu(idiag)
           wu(idiag) = s

           j = jwu(1)
           jwu(1) = jwu(idiag)
           jwu(idiag) = j

        end if

        len = min0(lenu,lenu0+lfil)
      call bsort2 (wu(2), jwu(2), lenu-1,len)
c
c  update u-matrix
c
        t = 0.0
        do 302 k=2, len
           if (ju0 .gt. iwk) goto 997
           jlu(ju0) = jwu(k)
           alu(ju0)  = wu(k)
           t = t+ abs(wu(k) )
           ju0 = ju0+1
 302      continue
c
c     save norm in wu (backwards). Norm is in fact average abs value
c
        wu(n+2-ii) = t / dble(len+1)
c
c     store inverse of diagonal element of u
c
        if (wu(1) .eq. 0.0) goto 999

        alu(ii) = 1.0/ wu(1)
c
c     update pointer to beginning of next row of U.
c
      jlu(ii+1) = ju0
c
c     reset double-pointer jr to zero (U-part)
c
      do 308 k=1, lenu
           jr(jwu(k)) = 0
 308      continue
c
c     end main loop
c
 500  continue
        ierr = 0
        return
c
c     zero pivot :
c
 900    ierr = ii
        return
c
c     incomprehensible error. Matrix must be wrong.
c
 995    ierr = -1
        return
c
c     insufficient storage in L.
c
 996    ierr = -2
        return
c
c     insufficient storage in U.
c
 997    ierr = -3
        return
c
c     illegal lfil entered.
c
 998    ierr = -4
        return
c
c     zero pivot encountered
c
 999    ierr = -5
        return
        end
      subroutine infdia (n,ja,ia,ind,idiag)

c*********************************************************************72
c
cc INFDIA obtains information on the diagonals of A.
c
c this subroutine finds the lengths of each of the 2*n-1 diagonals of A
c it also outputs the number of nonzero diagonals found.
c
c on entry:
c
c n      = dimension of the matrix a.
c
c a,    not needed here.
c ja,
c ia    = matrix stored in csr format
c
c on return:
c
c
c idiag = integer. number of nonzero diagonals found.
c
c ind   = integer array of length at least 2*n-1. The k-th entry in
c         ind contains the number of nonzero elements in the diagonal
c         number k, the numbering beeing from the lowermost diagonal
c         (bottom-left). In other words ind(k) = length of diagonal
c         whose offset wrt the main diagonal is = - n + k.
c
c           Y. Saad, Sep. 21 1989
c
      integer ia(*), ind(*), ja(*)

      n2= n+n-1
      do 1 i=1,n2
         ind(i) = 0
 1    continue
      do 3 i=1, n
         do 2 k=ia(i),ia(i+1)-1
            j = ja(k)
            ind(n+j-i) = ind(n+j-i) +1
 2       continue
 3    continue
c     count the nonzero ones.
      idiag = 0
      do 41 k=1, n2
         if (ind(k) .ne. 0) idiag = idiag+1
 41   continue
      return
      end
      subroutine ivperm (n, ix, perm)

c*********************************************************************72
c
cc IVPERM performs an in-place permutation of an integer vector.
c
c this subroutine performs an in-place permutation of an integer vector
c ix according to the permutation array perm(*), i.e., on return,
c the vector x satisfies,
c
c      ix(perm(j)) :== ix(j), j=1,2,.., n
c
c on entry:
c
c n       = length of vector x.
c perm       = integer array of length n containing the permutation  array.
c ix      = input vector
c
c on return:
c
c ix      = vector x permuted according to ix(perm(*)) :=  ix(*)
c
c
c           Y. Saad, Sep. 21 1989
c
      integer n, perm(n), ix(n)
      integer tmp, tmp1

      init      = 1
      tmp      = ix(init)
      ii        = perm(init)
      perm(init)= -perm(init)
      k         = 0
c
c loop
c
 6    k = k+1
c
c save the chased element --
c
      tmp1        = ix(ii)
      ix(ii)     = tmp
      next        = perm(ii)
      if (next .lt. 0 ) goto 65
c
c test for end
c
      if (k .gt. n) goto 101
      tmp       = tmp1
      perm(ii)  = - perm(ii)
      ii        = next
c
c end loop
c
      goto 6
c
c reinitilaize cycle --
c
 65   init      = init+1
      if (init .gt. n) goto 101
      if (perm(init) .lt. 0) goto 65
      tmp      = ix(init)
      ii      = perm(init)
      perm(init)=-perm(init)
      goto 6

 101  continue
      do 200 j=1, n
         perm(j) = -perm(j)
 200  continue

      return
      end
      subroutine jadcsr (nrow, idiag, a, ja, ia, iperm, ao, jao, iao)

c*********************************************************************72
c
cc JADSCR converts Jagged Diagonal Storage to Compressed Sparse Row.
c
c this subroutine converts a matrix stored in the jagged diagonal format
c to the compressed sparse row format.
c
c on entry:
c
c nrow         = integer. the row dimension of the matrix A.
c
c idiag   = integer. The  number of jagged diagonals in the data
c           structure a, ja, ia.
c
c a,
c ja,
c ia      = input matrix in jagged diagonal format.
c
c iperm   = permutation of the rows used to obtain the JAD ordering.
c
c on return:
c
c
c ao, jao,
c iao     = matrix in CSR format.
c
c determine first the pointers for output matrix. Go through the
c structure once:
c
      integer ja(*), jao(*), ia(idiag+1), iperm(nrow), iao(nrow+1)
      double precision a(*), ao(*)

      do 137 j=1,nrow
         jao(j) = 0
 137  continue
c
c     compute the lengths of each row of output matrix
c
      do 140 i=1, idiag
         len = ia(i+1)-ia(i)
         do 138 k=1,len
            jao(iperm(k)) = jao(iperm(k))+1
 138     continue
 140  continue
c
c     remember to permute
c
      kpos = 1
      iao(1) = 1
      do 141 i=1, nrow
         kpos = kpos+jao(i)
         iao(i+1) = kpos
 141  continue
c
c     copy elemnts one at a time.
c
      do 200 jj = 1, idiag
         k1 = ia(jj)-1
         len = ia(jj+1)-k1-1
         do 160 k=1,len
            kpos = iao(iperm(k))
            ao(kpos) = a(k1+k)
            jao(kpos) = ja(k1+k)
            iao(iperm(k)) = kpos+1
 160     continue
 200  continue
c
c     rewind pointers
c
      do 5 j=nrow,1,-1
         iao(j+1) = iao(j)
 5    continue
      iao(1) = 1
      return
      end
      subroutine ldsol (n,x,y,al,jal)

c*********************************************************************72
c
cc LDSOL solves L * x = y, for L a triangular matrix in MSR format.
c
c solves a (non-unit) lower triangular system by standard (sequential)
c forward elimination - matrix stored in MSR format
c with diagonal elements already inverted (otherwise do inversion,
c al(1:n) = 1.0/al(1:n),  before calling ldsol).
c
c On entry:
c
c n      = integer. dimension of problem.
c y      = double precision array containg the right hand side.
c
c al,
c jal,   = Lower triangular matrix stored in Modified Sparse Row
c          format.
c
c On return:
c
c      x = The solution of  L x = y .
c
      integer n, jal(*)
      double precision x(n), y(n), al(*)
      integer k, j
      double precision t

      x(1) = y(1)*al(1)
      do 150 k = 2, n
         t = y(k)
         do 100 j = jal(k), jal(k+1)-1
            t = t - al(j)*x(jal(j))
 100     continue
         x(k) = al(k)*t
 150  continue
      return
      end
      subroutine ldsolc (n,x,y,al,jal)

c*********************************************************************72
c
cc LDSOLC solves L*x = y;    L = nonunit Low. Triang. MSC format
c
c solves a (non-unit) lower triangular system by standard (sequential)
c forward elimination - matrix stored in Modified Sparse Column format
c with diagonal elements already inverted (otherwise do inversion,
c al(1:n) = 1.0/al(1:n),  before calling ldsol).
c
c
c On entry:
c
c n      = integer. dimension of problem.
c y      = double precision array containg the right hand side.
c
c al,
c jal,
c ial,    = Lower triangular matrix stored in Modified Sparse Column
c           format.
c
c On return:
c
c      x = The solution of  L x = y .
c
      integer n, jal(*)
      double precision x(n), y(n), al(*)
      integer k, j
      double precision t

      do 140 k=1,n
         x(k) = y(k)
 140  continue
      do 150 k = 1, n
         x(k) = x(k)*al(k)
         t = x(k)
         do 100 j = jal(k), jal(k+1)-1
            x(jal(j)) = x(jal(j)) - t*al(j)
 100     continue
 150  continue
c
      return
      end
      subroutine ldsoll (n,x,y,al,jal,nlev,lev,ilev)

c*********************************************************************72
c
cc LDSOLL solves L*x = y; L = triangular.
c
c On entry:
c
c n      = integer. dimension of problem.
c y      = double precision array containg the right hand side.
c
c al,
c jal,   = Lower triangular matrix stored in Modified Sparse Row
c          format.
c nlev   = number of levels in matrix
c lev    = integer array of length n, containing the permutation
c          that defines the levels in the level scheduling ordering.
c ilev   = pointer to beginning of levels in lev.
c          the numbers lev(i) to lev(i+1)-1 contain the row numbers
c          that belong to level number i, in the level shcheduling
c          ordering.
c
c On return:
c
c      x = The solution of  L x = y .
      integer n, nlev, jal(*), ilev(nlev+1), lev(n)
      double precision x(n), y(n), al(*)
      integer ii, jrow, i
      double precision t
c
c     outer loop goes through the levels. (SEQUENTIAL loop)
c
      do 150 ii=1, nlev
c
c     next loop executes within the same level. PARALLEL loop
c
         do 100 i=ilev(ii), ilev(ii+1)-1
            jrow = lev(i)
c
c compute inner product of row jrow with x
c
            t = y(jrow)
            do 130 k=jal(jrow), jal(jrow+1)-1
               t = t - al(k)*x(jal(k))
 130        continue
            x(jrow) = t*al(jrow)
 100     continue
 150  continue
      return
      end
      subroutine levels (n, jal, ial, nlev, lev, ilev, levnum)

c*********************************************************************72
c
cc LEVELS gets the level structure of a lower triangular matrix.
c
c levels gets the level structure of a lower triangular matrix
c for level scheduling in the parallel solution of triangular systems
c strict lower matrices (e.g. unit) as well matrices with their main
c diagonal are accepted.
c
c on entry:
c
c n        = integer. The row dimension of the matrix
c jal, ial =
c
c on return:
c
c nlev     = integer. number of levels found
c lev      = integer array of length n containing the level
c            scheduling permutation.
c ilev     = integer array. pointer to beginning of levels in lev.
c            the numbers lev(i) to lev(i+1)-1 contain the row numbers
c            that belong to level number i, in the level scheduling
c            ordering. The equations of the same level can be solved
c            in parallel, once those of all the previous levels have
c            been solved.
c work arrays:
c
c levnum   = integer array of length n (containing the level numbers
c            of each unknown on return)
      integer jal(*),ial(*), levnum(*), ilev(*), lev(*)

      do 10 i = 1, n
         levnum(i) = 0
 10   continue
c
c     compute level of each node --
c
      nlev = 0
      do 20 i = 1, n
         levi = 0
         do 15 j = ial(i), ial(i+1) - 1
            levi = max (levi, levnum(jal(j)))
 15      continue
         levi = levi+1
         levnum(i) = levi
         nlev = max(nlev,levi)
 20   continue
c  set data structure.
      do 21 j=1, nlev+1
         ilev(j) = 0
 21   continue
c  count  number   of elements in each level.
      do 22 j=1, n
         i = levnum(j)+1
         ilev(i) = ilev(i)+1
 22   continue
c  set up pointer for  each  level.
      ilev(1) = 1
      do 23 j=1, nlev
         ilev(j+1) = ilev(j)+ilev(j+1)
 23   continue

c  determine elements of each level.
      do 30 j=1,n
         i = levnum(j)
         lev(ilev(i)) = j
         ilev(i) = ilev(i)+1
 30   continue
c     reset pointers backwards
      do 35 j=nlev, 1, -1
         ilev(j+1) = ilev(j)
 35   continue
      return
      end
      subroutine lnkcsr (n, a, jcol, istart, link, ao, jao, iao)

c*********************************************************************72
c
cc LNKCSR converts linked list storage to Compressed Sparse Row format.
c
c this subroutine translates a matrix stored in linked list storage
c format into the compressed sparse row format.
c
c Coded by Y. Saad, Feb 21, 1991.
c
c
c on entry:
c
c n      = integer equal to the dimension of A.
c
c a      = double precision array of size nna containing the nonzero elements
c jcol      = integer array of size      nnz containing the column positions
c         of the corresponding elements in a.
c istart= integer array of size n poiting to the beginning of the rows.
c         istart(i) contains the position of the first element of
c         row i in data structure. (a, jcol, link).
c         if a row is empty istart(i) must be zero.
c link      = integer array of size nnz containing the links in the linked
c         list data structure. link(k) points to the next element
c         of the row after element ao(k), jcol(k). if link(k) = 0,
c         then there is no next element, i.e., ao(k), jcol(k) is
c         the last element of the current row.
c
c on return:
c
c ao, jao, iao = matrix stored in csr format:
c
c ao    = double precision array containing the values of the nonzero elements of
c         the matrix stored row-wise.
c jao      = integer array of size nnz containing the column indices.
c iao      = integer array of size n+1 containing the pointers array to the
c         beginning of each row. iao(i) is the address in ao,jao of
c         first element of row i.
c
c  NZMAX is not provided in the calling sequence, so the following line
c  had to be commented out.
c
c     double precision a(*), ao(nzmax)
      double precision A(*), AO(*)
      integer n, jcol(*), istart(n), link(*), jao(*), iao(*)
      integer irow, ipos, next
c
c first determine individual bandwidths and pointers.
c
      ipos = 1
      iao(1) = ipos
c
c     loop through all rows
c
      do 100 irow =1, n
c
c     unroll i-th row.
c
         next = istart(irow)
 10      if (next .eq. 0) goto 99
         jao(ipos) = jcol(next)
         ao(ipos)  = a(next)
         ipos = ipos+1
         next = link(next)
         goto 10
 99      iao(irow+1) = ipos
 100  continue
c
      return
      end
      subroutine lsol (n,x,y,al,jal,ial)

c*********************************************************************72
c
cc LSOL solves L*x = y ; L = lower unit triang. /  CSR format
c
c solves a unit lower triangular system by standard (sequential )
c forward elimination - matrix stored in CSR format.
c
c On entry:
c
c n      = integer. dimension of problem.
c y      = double precision array containg the right side.
c
c al,
c jal,
c ial,    = Lower triangular matrix stored in compressed sparse row
c          format.
c
c On return:
c
c      x  = The solution of  L x  = y.
c
      integer n, jal(*),ial(n+1)
      double precision  x(n), y(n), al(*)
      integer k, j
      double precision  t

      x(1) = y(1)
      do 150 k = 2, n
         t = y(k)
         do 100 j = ial(k), ial(k+1)-1
            t = t-al(j)*x(jal(j))
 100     continue
         x(k) = t
 150  continue

      return
      end
      subroutine lsolc (n,x,y,al,jal,ial)

c*********************************************************************72
c
cc LSOLC solves L*x = y where L = unit lower triang. CSC format
c
c solves a unit lower triangular system by standard (sequential )
c forward elimination - matrix stored in CSC format.
c
c On entry:
c
c n      = integer. dimension of problem.
c y      = double precision array containg the right side.
c
c al,
c jal,
c ial,    = Lower triangular matrix stored in compressed sparse column
c          format.
c
c On return:
c
c      x  = The solution of  L x  = y.
c
      integer n, jal(*),ial(*)
      double precision  x(n), y(n), al(*)
      integer k, j
      double precision t

      do 140 k=1,n
         x(k) = y(k)
 140  continue
      do 150 k = 1, n-1
         t = x(k)
         do 100 j = ial(k), ial(k+1)-1
            x(jal(j)) = x(jal(j)) - t*al(j)
 100     continue
 150  continue
c
      return
      end
      subroutine lusol0 (n, y, x, alu, jlu, ju)

c*********************************************************************72
c
cc LUSOL0 performs forward and backward solves for LU matrix produced by ILUT.
c
c performs a forward followed by a backward solve
c for LU matrix as produced by  ILUT
c
      INTEGER N
        double precision x(n), y(n), alu(*)
      integer jlu(*), ju(*)
        integer i,k
c
c forward solve
c
        do 40 i = 1, n
           x(i) = y(i)
           do 41 k=jlu(i),ju(i)-1
              x(i) = x(i) - alu(k)* x(jlu(k))
 41        continue
 40     continue
c
c     backward solve.
c
      do 90 i = n, 1, -1
         do 91 k=ju(i),jlu(i+1)-1
              x(i) = x(i) - alu(k)*x(jlu(k))
 91         continue
           x(i) = alu(i)*x(i)
 90     continue

        return
      end
      subroutine markgen (m, n, a, ja, ia)

c*********************************************************************72
c
cc MARKGEN is a matrix generator for a Markov random walk on a triang. grid
c
c this subroutine generates a test matrix that models a random
c walk on a triangular grid. This test example was used by
c G. W. Stewart ["{SRRIT} - a FORTRAN subroutine to calculate the
c dominant invariant subspaces of a real matrix",
c Tech. report. TR-514, University of Maryland (1978).] and in a few
c papers on eigenvalue problems by Y. Saad [see e.g. LAA, vol. 34,
c pp. 269-295 (1980) ]. These matrices provide reasonably easy
c test problems for eigenvalue algorithms. The transpose of the
c matrix  is stochastic and so it is known that one is an exact
c eigenvalue. One seeks the eigenvector of the transpose associated
c with the eigenvalue unity. The problem is to calculate the
c steady state probability distribution of the system, which is
c the eigevector associated with the eigenvalue one and scaled in
c such a way that the sum all the components is equal to one.
c
c parameters
c
c on entry :
c
c m     = integer. number of points in each direction.
c
c on return:
c
c n     = integer. The dimension of the matrix. (In fact n is known
c         to be equal to (m(m+1))/2      )
c a,
c ja,
c ia    = the matrix stored in CSR format.
c
c
c Notes: 1) the code will actually compute the transpose of the
c stochastic matrix that contains the transition probibilities.
c        2) It should also be possible to have a matrix generator
c with an additional parameter (basically redefining `half' below
c to be another parameter and changing the rest accordingly, but
c this is not as simple as it sounds). This is not likely to provide
c any more interesting matrices.
c
      double precision a(*), cst, pd, pu, half
      integer ja(*), ia(*)
c
      cst = 0.5/dble(m-1)
c
c  ix counts the grid point (natural ordering used), i.e.,
c  the row number of the matrix.
c
      ix = 0
      jax = 1
      ia(1) = jax
c
c     sweep y coordinates
c
      do 20 i=1,m
         jmax = m-i+1
c
c     sweep x coordinates
c
         do 10 j=1,jmax
            ix = ix + 1
            if (j .eq. jmax) goto 2
            pd = cst*dble(i+j-1)
c
c     north
c
            a(jax) = pd
            if (i.eq. 1) a(jax) = a(jax)+pd
            ja(jax) =  ix + 1
            jax = jax+1
c     east
            a(jax) = pd
            if (j .eq. 1) a(jax) = a(jax)+pd
            ja(jax) = ix + jmax
            jax = jax+1
c     south
 2          pu = 0.5 - cst*dble(i+j-3)
            if ( j .gt. 1) then
               a(jax) = pu
               ja(jax) = ix-1
               jax = jax+1
            end if
c     west
            if ( i .gt. 1) then
               a(jax) = pu
               ja(jax) = ix - jmax - 1
               jax = jax+1
            end if
            ia(ix+1) = jax
 10      continue
 20   continue
      n = ix
      return
      end
      subroutine matrf2(M,N,C,INDEX,ALPHA,NN,NZ,A,SNR,RNR,FEJLM)

c*********************************************************************72
c
cc MATRF2 generates sparse (rectangular or square) matrices.
c
c   PURPOSE
c
c   The subroutine generates sparse (rectangular or square) matrices.
c   The dimensions of the matrix and the average number of nonzero
c   elements per row can be specified by the user. Moreover, the user
c   can also change the sparsity pattern and the condition number of the
c   matrix. The non-zero elements of the desired matrix will be
c   accumulated (in an arbitrary order) in the first NZ positions of
c   array A. The column and the row numbers of the non-zero element
c   stored in A(I), I=1,...,NZ, will be found in SNR(I) and RNR(I),
c   respectively. The matrix generated by this subroutine is of the
c   class F(M,N,C,R,ALPHA) (see reference).
c
c   Note: If A is the sparse matrix of type F(M,N,C,R,ALPHA), then
c
c           min|A(i,j)| = 1/ALPHA,
c
c           max|A(i,j)| = max(INDEX*N - N,10*ALPHA).
c
c
c   CONTRIBUTOR: Ernest E. Rothman
c                Cornell Theory Center/Cornell National Supercomputer
c                Facility.
c                e-mail address: BITNET:   eer@cornellf
c                                INTERNET: eer@cornellf.tn.cornell.edu
c
c   minor modifications by Y. Saad. April 26, 1990.
c
c   Note: This subroutine has been copied from the following reference.
c         The allowable array sizes have been changed.
c
c   REFERENCE: Zlatev, Zahari; Schaumburg, Kjeld; Wasniewski, Jerzy;
c      "A testing Scheme for Subroutines Solving Large Linear Problems",
c      Computers and Chemistry, Vol. 5, No. 2-3, pp. 91-100, 1981.
c
c
c   INPUT PARAMETERS
c
c   M    - Integer. The number of rows in the desired matrix.
c          N < M+1 < 9000001 must be specified.
c
c   N    - Integer. The number of columns in the desired matrix.
c          21 < N < 9000001 must be specified.
c
c   C    - Integer. The sparsity pattern can be changed by means of this
c          parameter.  10 < C < N-10  must be specified.
c
c   INDEX - Integer.  The average number of non-zero elements per row in
c           the matrix will be equal to INDEX.
c           1 < INDEX < N-C-8 must be specified.
c
c   ALPHA - Real. The condition number of the matrix can be changed
c           BY THIS PARAMETER. ALPHA > 0.0 MUST BE SPECIFIED.
c           If ALPHA is approximately equal to 1.0 then the generated
c           matrix is well-conditioned. Large values of ALPHA will
c           usually produce ill-conditioned matrices. Note that no
c           round-off errors during the computations in this subroutine
c           are made if ALPHA = 2**I (where I is an arbitrary integer
c           which produces numbers in the machine range).
c
c   NN    - Integer. The length of arrays A, RNR, and SNR (see below).
c           INDEX*M+109 < NN < 9000001 must be specified.
c
c
c   OUTPUT PARAMETERS
c
c   NZ    - Integer. The number of non-zero elements in the matrix.
c
c   A(NN) - Real array. The non-zero elements of the matrix generated
c           are accumulated in the first NZ locations of array A.
c
c   SNR(NN) - INTEGER array. The column number of the non-zero element
c           kept in A(I), I=1,...NZ, is stored in SNR(I).
c
c   RNR(NN) - Integer array. The row number of the non-zero element
c           kept in A(I), I=1,...NZ, is stored in RNR(I).
c
c   FEJLM - Integer. FEJLM=0 indicates that the call is successful.
c           Error diagnostics are given by means of positive values of
c           this parameter as follows:
c             FEJLM = 1    -  N       is out of range.
c             FEJLM = 2    -  M       is out of range.
c             FEJLM = 3    -  C       is out of range.
c             FEJLM = 4    -  INDEX   is out of range.
c             FEJLM = 5    -  NN      is out of range.
c             FEJLM = 7    -  ALPHA   is out of range.
c
c
c
c
      double precision A, ALPHA, ALPHA1
      INTEGER M, N, NZ, C, NN, FEJLM, M1, NZ1, RR1, RR2, RR3, K
      INTEGER M2, N2
      INTEGER SNR, RNR
      DIMENSION A(NN), SNR(NN), RNR(NN)
      M1 = M
      FEJLM = 0
      NZ1 = INDEX*M + 110
      K = 1
      ALPHA1 = ALPHA
      INDEX1 = INDEX - 1
c
c  Check the parameters.
c
      IF(N.GE.22) GO TO 1
2     FEJLM = 1
      RETURN
1     IF(N.GT.9000000) GO TO 2
      IF(M.GE.N) GO TO 3
4     FEJLM = 2
      RETURN
3     IF(M.GT.9000000) GO TO 4
      IF(C.LT.11)GO TO 6
      IF(N-C.GE.11)GO TO 5
6     FEJLM = 3
      RETURN
5     IF(INDEX.LT.1) GO TO 12
      IF(N-C-INDEX.GE.9)GO TO 13
12    FEJLM = 4
13    IF(NN.GE.NZ1)GO TO 7
8     FEJLM = 5
      RETURN
7     IF(NN.GT.9000000)GO TO 8
      IF(ALPHA.GT.0.0)GO TO 9
      FEJLM = 6
      RETURN
9     CONTINUE
c
c  End of the error check. Begin to generate the non-zero elements of
c  the required matrix.
c
      DO 20 I=1,N
      A(I) = 1.0
      SNR(I) = I
20    RNR(I) = I
      NZ = N
      J1 = 1
      IF(INDEX1.EQ.0) GO TO 81
      DO 21 J = 1,INDEX1
      J1 = -J1
      DO 22 I=1,N
      A(NZ+I) = dble ( J1 * J * I )
      IF(I+C+J-1.LE.N)SNR(NZ+I) = I + C + J - 1
      IF(I+C+J-1.GT.N)SNR(NZ+I) = C + I + J - 1 - N
22    RNR(NZ + I) = I
21    NZ = NZ + N
81    RR1 = 10
      RR2 = NZ
      RR3 = 1
25    CONTINUE
      DO 26 I=1,RR1
      A(RR2 + I) = ALPHA* dble ( I )
      SNR(RR2+I) = N - RR1 + I
      RNR(RR2+I) = RR3
26    CONTINUE
      IF(RR1.EQ.1) GO TO 27
      RR2 = RR2 + RR1
      RR1 = RR1 - 1
      RR3 = RR3 + 1
      GO TO 25
27    NZ = NZ + 55
29    M1 = M1 - N
      ALPHA = 1.0/ALPHA
      IF(M1.LE.0) GO TO 28
      N2 = K*N
      IF(M1.GE.N)M2 = N
      IF(M1.LT.N)M2 = M1
      DO 30 I=1,M2
      A(NZ+I) = ALPHA* dble ( K + 1 )
      SNR(NZ + I) = I
30    RNR(NZ + I) = N2 + I
      NZ = NZ + M2
      IF(INDEX1.EQ.0) GO TO 82
      J1 = 1
      DO 41 J = 1,INDEX1
      J1 = -J1
      DO 42 I = 1,M2
      A(NZ+I) = ALPHA* dble (J*J1)*( dble ((K+1)*I)+1.0)
      IF(I+C+J-1.LE.N)SNR(NZ+I) = I + C + J - 1
      IF(I+C+J-1.GT.N)SNR(NZ+I) = C + I + J - 1 - N
42    RNR(NZ + I) = N2 + I
41    NZ = NZ +M2
82    K = K + 1
      GO TO 29
28    CONTINUE
      ALPHA = 1.0/ALPHA1
      RR1 = 1
      RR2 = NZ
35    CONTINUE
      DO 36 I = 1,RR1
      A(RR2+I) = ALPHA* dble (RR1+1-I)
      SNR(RR2+I) = I
      RNR(RR2+I) = N - 10 + RR1
36    CONTINUE
      IF(RR1.EQ.10) GO TO 34
      RR2 = RR2 + RR1
      RR1 = RR1 + 1
      GO TO 35
34    NZ = NZ + 55
      ALPHA = ALPHA1
      RETURN
      END
      subroutine mgsr (n, i0, i1, ss, r)

c*********************************************************************72
c
cc MGSR is a modified Gram - Schmidt with partial reorthogonalization.
c
c modified gram - schmidt  with  partial  reortho. the vector ss(*,i1) is
c orthogonalized against the first i vectors  of ss  (which  are  already
c orthogonal).  the coefficients of the orthogonalization are returned in
c the array r
c
      implicit double precision (a-h,o-z)
      double precision ss(n,1), r(1), hinorm, tet, ddot, t, sqrt
      do 53 j=1, i1
         r(j) = 0.0
 53   continue
      i = i1-1
      it = 0
 54   hinorm = 0.0
      it = it +1
      if (i .eq. 0) goto 56
c
      do 55 j=i0, i
         t = ddot(n, ss(1,j),1,ss(1,i1),1)
         hinorm = hinorm + t**2
         r(j) = r(j) + t
         call daxpy(n,-t,ss(1,j),1,ss(1,i1),1)
 55   continue
      t = ddot(n, ss(1,i1), 1, ss(1,i1), 1)
 56   continue
c
c     test for reorthogonalization see daniel et. al.
c     two reorthogonalization allowed.
c
      if (t*10.0 .le. hinorm .and. it .lt. 2) goto 54
      t =sqrt(t)
      r(i1)= t
      if (t .eq. 0.0) return
      t = 1.0/t
      do 57  k=1,n
         ss(k,i1) = ss(k,i1)*t
 57   continue
      return
      end
      subroutine milu0 (n, a, ja, ia, alu, jlu, ju, iw, ierr)

c*********************************************************************72
c
cc MILU0 is a simple milu(0) preconditioner.
c
c Note that this has been coded in such a way that it can be used
c with pgmres. Normally, since the data structure of a, ja, ia is
c the same as that of a, ja, ia, savings can be made. In fact with
c some definitions (not correct for general sparse matrices) all we
c need in addition to a, ja, ia is an additional diagonal.
c Ilu0 is not recommended for serious problems. It is only provided
c here for comparison purposes.
c
c on entry:
c
c n       = dimension of matrix
c a, ja,
c ia      = original matrix in compressed sparse row storage.
c
c on return:
c
c alu,jlu = matrix stored in Modified Sparse Row (MSR) format containing
c           the L and U factors together. The diagonal (stored in
c           alu(1:n) ) is inverted. Each i-th row of the alu,jlu matrix
c           contains the i-th row of L (excluding the diagonal entry=1)
c           followed by the i-th row of U.
c
c ju        = pointer to the diagonal elements in alu, jlu.
c
c ierr        = integer indicating error code on return
c           ierr = 0 --> normal return
c           ierr = k --> code encountered a zero pivot at step k.
c work arrays:
c
c iw          = integer work array of length n.
c
c Note (IMPORTANT):
c
c it is assumed that the the elements in the input matrix are ordered
c    in such a way that in each row the lower part comes first and
c    then the upper part. To get the correct ILU factorization, it is
c    also necessary to have the elements of L ordered by increasing
c    column number. It may therefore be necessary to sort the
c    elements of a, ja, ia prior to calling milu0. This can be
c    achieved by transposing the matrix twice using csrcsc.
c
      implicit double precision (a-h,o-z)
      double precision a(*), alu(*)
      integer ja(*), ia(*), ju(*), jlu(*), iw(*)

          ju0 = n+2
          jlu(1) = ju0
c initialize work vector to zero's
      do 31 i=1, n
 31           iw(i) = 0
c
c  MAIN LOOP
c
      do 500 ii = 1, n
           js = ju0
c
c generating row number ii or L and U.
c
           do 100 j=ia(ii),ia(ii+1)-1
c
c copy row ii of a, ja, ia into row ii of alu, jlu (L/U) matrix.
c
              jcol = ja(j)
              if (jcol .eq. ii) then
                 alu(ii) = a(j)
                 iw(jcol) = ii
                 ju(ii)  = ju0
              else
                 alu(ju0) = a(j)
                 jlu(ju0) = ja(j)
                 iw(jcol) = ju0
                 ju0 = ju0+1
              end if
 100       continue
           jlu(ii+1) = ju0
           jf = ju0-1
           jm = ju(ii)-1
c s accumulates fill-in values
           s = 0.0
           do 150 j=js, jm
              jrow = jlu(j)
              tl = alu(j)*alu(jrow)
              alu(j) = tl
c  perform linear combination
              do 140 jj = ju(jrow), jlu(jrow+1)-1
                 jw = iw(jlu(jj))
                 if (jw .ne. 0) then
                       alu(jw) = alu(jw) - tl*alu(jj)
                    else
                       s = s + tl*alu(jj)
                    end if
 140          continue
 150       continue
c  invert and store diagonal element.
           alu(ii) = alu(ii)-s
           if (alu(ii) .eq. 0.0) goto 600
           alu(ii) = 1.0/alu(ii)
c  reset pointer iw to zero
           iw(ii) = 0
           do 201 i = js, jf
 201          iw(jlu(i)) = 0
 500       continue
           ierr = 0
           return
c     zero pivot :
 600       ierr = ii
           return
           end
      subroutine msrcsr (n,a,ja,ao,jao,iao,wk)

c*********************************************************************72
c
cc MSRCSR converts Modified Sparse Row to Compressed Sparse Row.
c
c converts a compressed matrix using a separated diagonal
c (modified sparse row format) in the Compressed Sparse Row
c format.
c does not check for zero elements in the diagonal.
c
c
c on entry :
c
c n        = row dimension of matrix
c ao, jao  = sparse matrix in msr sparse storage format
c           see routine csrmsr for details
c
c on return :
c
c a, ja, ia = matrix in csr format. note that the
c           algorithm is in place: ao, jao can be the same
c            as a, ja, in which case it will be overwritten on it
c            upon return.
c
c             here nnz = number of nonzero elements+1
c work arrays:
c
c wk      = ouble precision work array of length n
c
c notes:
c  In place algorithm (see a, ja, ia).
c

      double precision a(*),ao(*),wk(n)
      integer ja(*),jao(*),iao(n+1)

      logical added
      do 1 i=1,n
         wk(i) = a(i)
 1    continue
      iao(1) = 1
      iptr = 1

      do 500 ii=1,n
         added = .false.
         idiag = iptr + (ja(ii+1)-ja(ii))
         do 100 k=ja(ii),ja(ii+1)-1
            j = ja(k)
            if (j .lt. ii) then
               ao(iptr) = a(k)
               jao(iptr) = j
               iptr = iptr+1
            elseif (added) then
               ao(iptr) = a(k)
               jao(iptr) = j
               iptr = iptr+1
            else
c add diag element - only reserve a position for it.
               idiag = iptr
               iptr = iptr+1
               added = .true.
c     then other element
               ao(iptr) = a(k)
               jao(iptr) = j
               iptr = iptr+1
            end if
 100     continue
         ao(idiag) = wk(ii)
         jao(idiag) = ii
         if (.not. added) iptr = iptr+1
         iao(ii+1) = iptr
 500  continue
      return
      end
      subroutine ope (n, x, y, a, ja, ia)

c*********************************************************************72
c
cc OPE sparse matrix * vector multiplication
c

      INTEGER N
      double precision  x(n), y(n), a(*)
      integer k1, k2, ja(*), ia(n+1)

      do 100 i=1,n
           k1 = ia(i)
           k2 = ia(i+1) -1
           y(i) = 0.0
           do 99 k=k1, k2
 99           y(i) = y(i) + a(k) * x(ja(k))
 100      continue
      return
      end
      subroutine pgmres(n, im, rhs, sol, vv, eps, maxits, iout,
     *                    aa, ja, ia, alu, jlu, ju, ierr)

c*********************************************************************72
c
cc PGMRES is an ILUT - Preconditioned GMRES solver.
c
c                 *** ILUT - Preconditioned GMRES ***
c
c*
c This is a simple version of the ILUT preconditioned GMRES algorithm.
c The ILUT preconditioner uses a dual strategy for dropping elements
c instead  of the usual level of-fill-in approach. See details in ILUT
c subroutine documentation. PGMRES uses the L and U matrices generated
c from the subroutine ILUT to precondition the GMRES algorithm.
c The preconditioning is applied to the right. The stopping criterion
c utilized is based simply on reducing the residual norm by epsilon.
c This preconditioning is more reliable than ilu0 but requires more
c storage. It seems to be much less prone to difficulties related to
c strong nonsymmetries in the matrix. We recommend using a nonzero tol
c (tol=.005 or .001 usually give good results) in ILUT. Use a large
c lfil whenever possible (e.g. lfil = 5 to 10). The higher lfil the
c more reliable the code is. Efficiency may also be much improved.
c Note that lfil=n and tol=0.0 in ILUT  will yield the same factors as
c Gaussian elimination without pivoting.
c
c ILU(0) and MILU(0) are also provided for comparison purposes
c USAGE: first call ILUT or ILU0 or MILU0 to set up preconditioner and
c then call pgmres.
c
c Coded by Y. Saad - This version dated May, 7, 1990.
c*
c parameters
c
c on entry:
c
c
c n     == integer. The dimension of the matrix.
c im    == size of krylov subspace:  should not exceed 50 in this
c          version (can be reset by changing parameter command for
c          kmax below)
c rhs   == double precision vector of length n containing the right hand side.
c          Destroyed on return.
c sol   == double precision vector of length n containing an initial guess to the
c          solution on input. approximate solution on output
c eps   == tolerance for stopping criterion. process is stopped
c          as soon as ( ||.|| is the euclidean norm):
c          || current residual||/||initial residual|| <= eps
c maxits== maximum number of iterations allowed
c iout  == output unit number number for printing intermediate results
c          if (iout .le. 0) nothing is printed out.
c
c aa, ja,
c ia    == the input matrix in compressed sparse row format:
c          aa(1:nnz)  = nonzero elements of A stored row-wise in order
c          ja(1:nnz) = corresponding column indices.
c          ia(1:n+1) = pointer to beginning of each row in aa and ja.
c          here nnz = number of nonzero elements in A = ia(n+1)-ia(1)
c
c alu,jlu== A matrix stored in Modified Sparse Row format containing
c           the L and U factors, as computed by subroutine ilut.
c
c ju     == integer array of length n containing the pointers to
c           the beginning of each row of U in alu, jlu as computed
c           by subroutine ILUT.
c
c on return:
c
c sol   == contains an approximate solution (upon successful return).
c ierr  == integer. Error message with the following meaning.
c          ierr = 0 --> successful return.
c          ierr = 1 --> convergence not achieved in itmax iterations.
c          ierr =-1 --> the initial guess seems to be the exact
c                       solution (initial residual computed was zero)
c
c
c
c work arrays:
c
c vv    == work array of length  n x (im+1) (used to store the Arnoli
c          basis)
c
c subroutines called :
c ope    : matrix by vector multiplication delivers y=ax, given x
c lusol0 : combined forward and backward solves (Preconditioning ope.)
c BLAS2  routines.
c
       implicit double precision (a-h,o-z)
       integer n, im, maxits, iout, ierr, ja(*), ia(n+1), jlu(*), ju(n)
       double precision vv(n,*), rhs(n), sol(n), aa(*), alu(*), eps

       parameter (kmax=50)
      PARAMETER (EPSMAC=1.0E-16)
c
       double precision hh(kmax+1,kmax), c(kmax), s(kmax), rs(kmax+1),t
c
c arnoldi size should not exceed kmax=50 in this version..
c to reset modify paramter kmax accordingly.
c
       n1 = n + 1
       its = 0
c
c outer loop starts here..
c  compute initial residual vector
       call ope (n, sol, vv, aa, ja, ia)
       do 21 j=1,n
          vv(j,1) = rhs(j) - vv(j,1)
 21    continue

 20    ro = sqrt( ddot(n, vv, 1, vv, 1) )
       if (iout .gt. 0 .and. its .eq. 0)
     *      write(iout, 199) its, ro
       if (ro .eq. 0.0) goto 999
       t = 1.0/ ro
       do 210 j=1, n
          vv(j,1) = vv(j,1)*t
 210   continue
       if (its .eq. 0) eps1=eps*ro
c     ** initialize 1-st term  of rhs of hessenberg system..
       rs(1) = ro
       i = 0
 4     i=i+1
       its = its + 1
       i1 = i + 1
       call lusol0 (n, vv(1,i), rhs, alu, jlu, ju)
       call ope (n, rhs, vv(1,i1), aa, ja, ia)
c
c     modified gram - schmidt...
c
       do 55 j=1, i
          t = ddot(n, vv(1,j),1,vv(1,i1),1)
          hh(j,i) = t
          call daxpy(n, -t, vv(1,j), 1, vv(1,i1), 1)
 55    continue
       t = sqrt(ddot(n, vv(1,i1), 1, vv(1,i1), 1))
       hh(i1,i) = t
       if ( t .eq. 0.0) goto 58
       t = 1.0/t
       do 57  k=1,n
          vv(k,i1) = vv(k,i1)*t
 57    continue
c
c     done with modified gram schimd and arnoldi step..
c     now  update factorization of hh
c
 58    if (i .eq. 1) goto 121
c
c  perform previous transformations  on i-th column of h
c
       do 66 k=2,i
          k1 = k-1
          t = hh(k1,i)
          hh(k1,i) = c(k1)*t + s(k1)*hh(k,i)
          hh(k,i) = -s(k1)*t + c(k1)*hh(k,i)
 66    continue
 121   gam = sqrt(hh(i,i)**2 + hh(i1,i)**2)
c
c     if gamma is zero then any small value will do...
c     will affect only residual estimate
c
       if (gam .eq. 0.0) gam = epsmac
c
c     get  next plane rotation
c
       c(i) = hh(i,i)/gam
       s(i) = hh(i1,i)/gam
       rs(i1) = -s(i)*rs(i)
       rs(i) =  c(i)*rs(i)
c
c     detrermine residual norm and test for convergence-
c
       hh(i,i) = c(i)*hh(i,i) + s(i)*hh(i1,i)
       ro = abs(rs(i1))
 131   format(1h ,2e14.4)
       if (iout .gt. 0)
     *      write(iout, 199) its, ro
       if (i .lt. im .and. (ro .gt. eps1))  goto 4
c
c     now compute solution. first solve upper triangular system.
c
       rs(i) = rs(i)/hh(i,i)
       do 30 ii=2,i
          k=i-ii+1
          k1 = k+1
          t=rs(k)
          do 40 j=k1,i
             t = t-hh(k,j)*rs(j)
 40       continue
          rs(k) = t/hh(k,k)
 30    continue
c
c     form linear combination of v(*,i)'s to get solution
c
       t = rs(1)
       do 15 k=1, n
          rhs(k) = vv(k,1)*t
 15    continue
       do 16 j=2, i
          t = rs(j)
          do 161 k=1, n
             rhs(k) = rhs(k)+t*vv(k,j)
 161      continue
 16    continue
c
c     call preconditioner.
c
       call lusol0 (n, rhs, rhs, alu, jlu, ju)
       do 17 k=1, n
          sol(k) = sol(k) + rhs(k)
 17    continue
c
c     restart outer loop  when necessary
c
       if (ro .le. eps1) goto 990
       if (its .gt. maxits) goto 991
c
c     else compute residual vector and continue..
c
       do 24 j=1,i
          jj = i1-j+1
          rs(jj-1) = -s(jj-1)*rs(jj)
          rs(jj) = c(jj-1)*rs(jj)
 24    continue
       do 25  j=1,i1
          t = rs(j)
          if (j .eq. 1)  t = t-1.0
          call daxpy (n, t, vv(1,j), 1,  vv, 1)
 25    continue
 199   format(' its =', i4, ' res. norm =', G14.6)
c     restart outer loop.
       goto 20
 990   ierr = 0
       return
 991   ierr = 1
       return
 999   continue
       ierr = -1
       return
       end
      subroutine pltmt (nrow,ncol,mode,ja,ia,title,key,type,
     1     job, iounit)

c*********************************************************************72
c
cc PLTMT creates a 'pic' plot of a matrix.
c
c this subroutine creates a 'pic' file for plotting the pattern of
c a sparse matrix stored in general sparse format. It is not intended
c to be a means of plotting large matrices (It is very inefficient).
c It is however useful for small matrices and can be used for example
c for inserting matrix plots in a text. The size of the plot can be
c 7in x 7in or 5 in x 5in .. There is also an option for writing a
c 3-line header in troff (see description of parameter job).
c Author: Youcef Saad - Date: Sept., 1989
c See SPARSKIT/UNSUPP/ for a version of this to produce a post-script
c file.
c
c nrow   = number of rows in matrix
c
c ncol       = number of columns in matrix
c
c mode   = integer indicating whether the matrix is stored
c          row-wise (mode = 0) or column-wise (mode=1)
c
c ja     = column indices of nonzero elements when matrix is
c         stored rowise. Row indices if stores column-wise.
c ia     = integer array of containing the pointers to the
c         beginning of the columns in arrays a, ja.
c
c title  = character*71 = title of matrix test ( character a*71 ).
c key    = character*8  = key of matrix
c type   = character*3  = type of matrix.
c
c job    = this integer parameter allows to set a few minor
c          options. First it tells pltmt whether or not to
c          reduce the plot. The standard size of 7in is then
c          replaced by a 5in plot. It also tells pltmt whether or
c          not to append to the pic file a few 'troff' lines that
c          produce a centered caption includingg the title, key and
c          types as well as the size and number of nonzero elements.
c          job =  0 : do not reduce and do not make caption.
c          job =  1 : reduce and do not make caption.
c          job = 10 : do not reduce and make caption
c          job = 11 : reduce and make caption.
c          (i.e. trailing digit for reduction, leading digit for caption)
c
c iounit = logical unit number where to write the matrix into.
c
c example of usage .
c
c In the fortran code:
c  a) read a Harwell/Boeing matrix
c          call readmt (.....)
c         iout = 13
c  b) generate pic file:
c          call  pltmt (nrow,ncol,mode,ja,ia,title,key,type,iout)
c         stop
c
c Then in a unix environment plot the matrix by the command
c
c      pic FOR013.DAT | troff -me | lpr -Ppsx
c
c notes: 1) Plots square as well as rectangular matrices.
c            (however not as much tested with rectangular matrices.)
c        2) the dot-size is adapted according to the size of the
c            matrix.
c        3) This is not meant at all as a way of plotting large
c            matrices. The pic file generaled will have one line for
c            each nonzero element. It is  only meant for use in
c           such things as document poreparations etc..
c         4) The caption written will print the 71 character long
c            title. This may not be centered correctly if the
c            title has trailing blanks (a problem with Troff).
c            if you want the title centered then you can center
c            the string in title before calling pltmt.
c
      integer ja(*), ia(*)
      character key*8,title*72,type*3
      double precision x, y

      n = ncol
      if (mode .eq. 0) n = nrow
      nnz = ia(n+1) - ia(1)
      maxdim = max0 (nrow, ncol)
      xnrow = dble(nrow)
      xncol = dble(ncol)
      ptsize = 0.08
      hscale = (7.0 -2.0*ptsize)/dble(maxdim-1)
      vscale = hscale
      xwid  = ptsize + dble(ncol-1)*hscale + ptsize
      xht   = ptsize + dble(nrow-1)*vscale + ptsize
      xshift = (7.0-xwid)/2.0
      yshift = (7.0-xht)/2.0

      if (mod(job,10) .eq. 1) then
         write (iounit,88)
      else
         write (iounit,89)
      end if
 88   format('.PS 5in',/,'.po 1.8i')
 89   format('.PS',/,'.po 0.7i')
      write(iounit,90)
 90   format('box invisible wid 7.0 ht 7.0 with .sw at (0.0,0.0) ')
      write(iounit,91) xwid, xht, xshift, yshift
 91   format('box wid ',f5.2,' ht ',f5.2,
     *     ' with .sw at (',f5.2,',',f5.2,')' )
c
c     shift points slightly to account for size of dot , etc..
c
      tiny = 0.03
      if (mod(job,10) .eq. 1) tiny = 0.05
      xshift = xshift + ptsize - tiny
      yshift = yshift + ptsize + tiny

      ips = 8
      if (maxdim .le. 500) ips = 10
      if (maxdim .le. 300) ips = 12
      if (maxdim .le. 100) ips = 16
      if (maxdim .lt. 50) ips = 24
      write(iounit,92) ips
 92   format('.ps ',i2)
c
c  plottingloop
c
      do 1 ii=1, n
         istart = ia(ii)
         ilast  = ia(ii+1)-1
         if (mode .ne. 0) then
            x = dble(ii-1)
            do 2 k=istart, ilast
               y = xnrow-dble(ja(k))
               write(iounit,128) xshift+x*hscale, yshift+y*vscale
 2          continue
         else
            y = xnrow - dble(ii)
            do 3 k=istart, ilast
               x = dble(ja(k)-1)
               write(iounit,128) xshift+x*hscale, yshift+y*vscale
 3          continue
         end if
 1    continue

 128  format(7h"." at ,f6.3,',',f6.3,8h ljust  )
      write (iounit, 129)
 129  format('.PE')
c     quit if caption not desired.
      if ( (job/10) .ne. 1)  return

      write(iounit,127) key, type, title
      write(iounit,130) nrow,ncol,nnz
 127  format('.sp 4'/'.ll 7i'/'.ps 12'/'.po 0.7i'/'.ce 3'/,
     *     'Matrix:  ',a8,',  Type:  ',a3,/,a72)
 130  format('Dimension: ',i4,' x ',i4,',  Nonzero elements: ',i5)
      return
      end
      subroutine pltmtps (nrow,ncol,mode,ja,ia,title,key,type,
     1                        job, iounit)

c*********************************************************************72
c
cc PLTMTPS creates a PostScript plot of a sparse matrix.
c
c this subroutine creates a 'PS' file for plotting the pattern of
c a sparse matrix stored in general sparse format. It can be used
c for inserting matrix plots in a text. The size of the plot can be
c 7in x 7in or 5 in x 5in ..
c
c Adapted from pltmt in module INOUT by Paul Frederickson. March, 1990
c         + slight modifications by Y. Saad.
c
c nrow   = number of rows in matrix
c
c ncol       = number of columns in matrix
c
c mode   = integer indicating whether the matrix is stored
c          row-wise (mode = 0) or column-wise (mode=1)
c
c ja     = column indices of nonzero elements when matrix is
c         stored rowise. Row indices if stores column-wise.
c ia     = integer array of containing the pointers to the
c         beginning of the columns in arrays a, ja.
c
c title  = character*72 = title of matrix test ( character a*72 ).
c key    = character*8  = key of matrix
c type   = character*3  = type of matrix.
c
c job    =  integer. tells pltmt whether or not to reduce the plot.
c           if enabled then the standard size of 7in will be
c           replaced by a 5in plot.
c          job =  0 : do not reduce
c          job =  1 : reduce plot to 5 inches.
c
c iounit = logical unit number where to write the matrix into.
c
c notes: 1) Plots square as well as rectangular matrices.
c        2) Does not writer a caption yet.
c        3) No bounding box put in yet
c
      integer ja(*), ia(*)
      character key*8,title*72,type*3
      double precision x, y, delta

      n = ncol
      if (mode .eq. 0) n = nrow
      nnz = ia(n+1) - ia(1)
      maxdim = max0 (nrow, ncol)
      m = 1 + maxdim
c keep this test as in old pltmt (for future changes).
       if (mod(job,10) .eq. 1) then
        delta = 72*5.0/(2.0+maxdim)
       else
        delta = 72*7.0/(2.0+maxdim)
       end if

      write(iounit,*)'%!PS'
       write(iounit,*)' gsave 50 50 translate'
      write(iounit,*) delta, delta, ' scale'
      write(iounit,*) ' 0.25 setlinewidth'

       if (mod(job,10) .eq. 1) then
          write (iounit,*) ' 23 55 translate'
       else
          write (iounit,*) ' 2 35 translate'
       end if

      write(iounit,*) ' newpath'
      write(iounit,*) 0,0,' moveto'
      write(iounit,*) m,0,' lineto'
      write(iounit,*) m,m,' lineto'
      write(iounit,*) 0,m,' lineto'
      write(iounit,*) ' closepath stroke'
      write(iounit,*) ' 1 1 translate'
      write(iounit,*) ' 0.5 setlinewidth'
      write(iounit,*) ' /p {moveto 0 -.25 rmoveto '
      write(iounit,*) '            0  .50 rlineto stroke} def'
c
c  plotting loop
c
           do 1 ii=1, n
           istart = ia(ii)
           ilast  = ia(ii+1)-1
             if (mode .ne. 0) then
                do 2 k=istart, ilast
            write(iounit,*) ii-1, nrow-ja(k), ' p'
 2              continue
          else
c             y = xnrow - dble(ii)
              do 3 k=istart, ilast
c               x = dble(ja(k)-1)
            write(iounit,*) ja(k)-1, nrow-ii, ' p'
 3              continue
            end if
 1      continue

      write(iounit,*)' showpage grestore'
c
c quit if caption not desired.
c      if ( (job/10) .ne. 1)  return
c
c      write(iounit,127) key, type, title
c      write(iounit,130) nrow,ncol,nnz
c 127         format('.sp 4'/'.ll 7i'/'.ps 12'/'.po 0.7i'/'.ce 3'/,
c     *        'Matrix:  ',a8,',  Type:  ',a3,/,a71)
 130    format('Dimension: ',i4,' x ',i4',  Nonzero elements: ',i5)
      return
      end
      subroutine project(n,m,u,v,w)

c*********************************************************************72
c
cc PROJECT computes the matrix-vector product w = U * v.
c
      implicit double precision (a-h,o-z)
      double precision u(n,*),v(*),w(*)

      do k=1,n
         w(k) = 0.0
      end do

      do j=1,m
         do k=1,n
            w(k) = w(k) + v(j) * u(k,j)
         end do
      end do

      return
      end
      subroutine prtmt (nrow,ncol,a,ja,ia,rhs,guesol,title,key,type,
     1     ifmt,job,iounit)

c*********************************************************************72
c
cc PRTMT writes a matrix in Harwell-Boeing format into a file.
c
c writes a matrix in Harwell-Boeing format into a file.
c assumes that the matrix is stored in COMPRESSED SPARSE COLUMN FORMAT.
c some limited functionality for right hand sides.
c Author: Youcef Saad - Date: Sept., 1989 - updated Oct. 31, 1989 to
c cope with new format.
c
c on entry:
c
c nrow   = number of rows in matrix
c ncol       = number of columns in matrix
c a       = double precision array containing the values of the matrix stored
c          columnwise
c ja        = integer array of the same length as a containing the row indices
c          of the corresponding matrix elements of array a.
c ia     = integer array of containing the pointers to the beginning of
c         the columns in arrays a, ja.
c rhs    = double precision array  containing the right hand side (s) and optionally
c          the associated initial guesses and/or exact solutions
c          in this order. See also guesol for details. the vector rhs will
c          be used only if job .gt. 2 (see below). Only full storage for
c          the right hand sides is supported.
c
c guesol = a 2-character string indicating whether an initial guess
c          (1-st character) and / or the exact solution (2-nd)
c          character) is provided with the right hand side.
c         if the first character of guesol is 'G' it means that an
c          an intial guess is provided for each right hand sides.
c          These are assumed to be appended to the right hand sides in
c          the array rhs.
c         if the second character of guesol is 'X' it means that an
c          exact solution is provided for each right hand side.
c          These are assumed to be appended to the right hand sides
c          and the initial guesses (if any) in the array rhs.
c
c title  = character*71 = title of matrix test ( character a*71 ).
c key    = character*8  = key of matrix
c type   = charatcer*3  = type of matrix.
c
c ifmt       = integer specifying the format chosen for the real values
c         to be output (i.e., for a, and for rhs-guess-sol if
c          applicable). the meaning of ifmt is as follows.
c        * if (ifmt .lt. 100) then the E descriptor is used,
c           format Ed.m, in which the length (m) of the mantissa is
c           precisely the integer ifmt (and d = ifmt+6)
c        * if (ifmt .gt. 100) then prtmt will use the
c           F- descriptor (format Fd.m) in which the length of the
c           mantissa (m) is the integer mod(ifmt,100) and the length
c           of the integer part is k=ifmt/100 (and d = k+m+2)
c          Thus  ifmt= 4   means  E10.4  +.xxxxD+ee    while
c                ifmt=104  means  F7.4   +x.xxxx
c                ifmt=205  means  F9.5   +xx.xxxxx
c          Note: formats for ja, and ia are internally computed.
c
c job       = integer to indicate whether matrix values and
c         a right hand side is available to be written
c          job = 1   write srtucture only, i.e., the arrays ja and ia.
c          job = 2   write matrix including values, i.e., a, ja, ia
c          job = 3   write matrix and one right hand side: a,ja,ia,rhs.
c         job = nrhs+2 write matrix and nrhs successive right hand sides
c         Note that there cannot be any right hand side if the matrix
c         has no values. Also the initial guess and exact solutions when
c          provided are for each right hand side. For example if nrhs=2
c          and guesol='GX' there are 6 vectors to write.
c
c
c iounit = logical unit number where to write the matrix into.
c
c on return:
c
c the matrix a, ja, ia will be written in output unit iounit
c in the Harwell-Boeing format. Noe of the inputs is modofied.
c
c Notes: 1) This code attempts to pack as many elements as possible per
c        80-character line.
c        2) this code attempts to avoid as much as possible to put
c        blanks in the formats that are written in the 4-line header
c       (This is done for purely esthetical reasons since blanks
c        are ignored in format descriptors.)
c        3) sparse formats for righr hand sides and guesses not
c        suported.
c
      character title*72,key*8,type*3,ptrfmt*16,indfmt*16,valfmt*20,
     *              guesol*2, rhstyp*3
      integer totcrd, ptrcrd, indcrd, valcrd, rhscrd, nrow, ncol,
     1     nnz, nrhs, len, nperli
      integer ja(*), ia(*)
      double precision a(*),rhs(*)
c
c     compute pointer format
c
      nnz    = ia(ncol+1) -1
      len    = int ( dlog10(0.1+dble(nnz+1))) + 1
      nperli = 80/len
      ptrcrd = ncol/nperli + 1
      if (len .gt. 9) then
         assign 101 to ix
      else
         assign 100 to ix
      end if
      write (ptrfmt,ix) nperli,len
 100  format(1h(,i2,1HI,i1,1h) )
 101  format(1h(,i2,1HI,i2,1h) )
c
c compute ROW index format
c
      len    = int ( dlog10(0.1+dble(nrow) )) + 1
      nperli = min0(80/len,nnz)
      indcrd = (nnz-1)/nperli+1
      write (indfmt,100) nperli,len
c
c compute values and rhs format (using the same for both)
c
      valcrd      = 0
      rhscrd  = 0
c quit this part if no values provided.
      if (job .le. 1) goto 20

      if (ifmt .ge. 100) then
         ihead = ifmt/100
         ifmt = ifmt-100*ihead
         len = ihead+ifmt+2
         nperli = 80/len

         if (len .le. 9 ) then
            assign 102 to ix
         elseif (ifmt .le. 9) then
            assign 103 to ix
         else
            assign 104 to ix
         end if

         write(valfmt,ix) nperli,len,ifmt
 102     format(1h(,i2,1hF,i1,1h.,i1,1h) )
 103     format(1h(,i2,1hF,i2,1h.,i1,1h) )
 104     format(1h(,i2,1hF,i2,1h.,i2,1h) )

      else
         len = ifmt + 6
         nperli = 80/len
c     try to minimize the blanks in the format strings.
         if (nperli .le. 9) then
          if (len .le. 9 ) then
             assign 105 to ix
          elseif (ifmt .le. 9) then
             assign 106 to ix
          else
             assign 107 to ix
          end if
       else
          if (len .le. 9 ) then
             assign 108 to ix
          elseif (ifmt .le. 9) then
             assign 109 to ix
          else
               assign 110 to ix
            end if
         end if

         write(valfmt,ix) nperli,len,ifmt
 105     format(1h(,i1,1hE,i1,1h.,i1,1h) )
 106     format(1h(,i1,1hE,i2,1h.,i1,1h) )
 107     format(1h(,i1,1hE,i2,1h.,i2,1h) )
 108     format(1h(,i2,1hE,i1,1h.,i1,1h) )
 109     format(1h(,i2,1hE,i2,1h.,i1,1h) )
 110     format(1h(,i2,1hE,i2,1h.,i2,1h) )

      end if
      valcrd = (nnz-1)/nperli+1
      nrhs   = job -2
      if (nrhs .ge. 1) then
         i = (nrhs*nrow-1)/nperli+1
         rhscrd = i
         if (guesol(1:1) .eq. 'G') rhscrd = rhscrd+i
         if (guesol(2:2) .eq. 'X') rhscrd = rhscrd+i
         rhstyp = 'F'//guesol
      end if
 20   continue
c
      totcrd = ptrcrd+indcrd+valcrd+rhscrd
c     write 4-line or five line header
      write(iounit,10) title,key,totcrd,ptrcrd,indcrd,valcrd,
     1     rhscrd,type,nrow,ncol,nnz,nrhs,ptrfmt,indfmt,valfmt,valfmt
c
      if (nrhs .ge. 1) write (iounit,11) rhstyp, nrhs
 10   format (a72, a8 / 5i14 / a3, 11x, 4i14 / 2a16, 2a20)
 11   format(A3,11x,i4)
c
      write(iounit,ptrfmt) (ia (i), i = 1, ncol+1)
      write(iounit,indfmt) (ja (i), i = 1, nnz)
      if (job .le. 1) return
      write(iounit,valfmt) (a(i), i = 1, nnz)
      if (job .le. 2) return
      len = nrow*nrhs
      next = 1
      iend = len
      write(iounit,valfmt) (rhs(i), i = next, iend)
c
c     write initial guesses if available
c
      if (guesol(1:1) .eq. 'G') then
         next = next+len
         iend = iend+ len
         write(iounit,valfmt) (rhs(i), i = next, iend)
      end if
c
c     write exact solutions if available
c
      if (guesol(2:2) .eq. 'X')then
         next = next+len
         iend = iend+ len
         write(iounit,valfmt) (rhs(i), i = next, iend)
      end if

      return
      end
      subroutine readmt (nmax,nzmax,job,iounit,a,ja,ia,rhs,nrhs,
     *                     guesol,nrow,ncol,nnz,title,key,type,ierr)

c*********************************************************************72
c
cc READMT reads a Harwell/Boeing sparse matrix file.
c
c this subroutine reads  a boeing/harwell matrix. handles right hand
c sides in full format only (no sparse right hand sides).
c Author: Youcef Saad - Date: Sept., 1989
c         updated Oct 31, 1989.
c
c on entry:
c
c nmax        =  max column dimension  allowed for matrix. The array ia should
c          be of length at least ncol+1 (see below) if job.gt.0
c nzmax       = max number of nonzeros elements allowed. the arrays a,
c          and ja should be of length equal to nnz (see below) if these
c          arrays are to be read (see job).
c
c job       = integer to indicate what is to be read. (note: job is an
c          input and output parameter, it can be modified on return)
c          job = 0    read the values of ncol, nrow, nnz title, key,
c                     type and return. matrix is not read and arrays
c                     a, ja, ia, rhs are not touched.
c          job = 1    read srtucture only, i.e., the arrays ja and ia.
c          job = 2    read matrix including values, i.e., a, ja, ia
c          job = 3    read matrix and right hand sides: a,ja,ia,rhs.
c                  rhs may contain initial guesses and exact
c                     solutions appended to the actual right hand sides.
c                  this will be indicated by the output parameter
c                     guesol [see below].
c
c nrhs   = integer. nrhs is an input as well as ouput parameter.
c          at input nrhs contains the total length of the array rhs.
c          See also ierr and nrhs in output parameters.
c
c iounit = logical unit number where to read the matrix from.
c
c on return:
c
c job    = on return job may be modified to the highest job it could
c          do: if job=2 on entry but no matrix values are available it
c          is reset to job=1 on return. Similarly of job=3 but no rhs
c          is provided then it is rest to job=2 or job=1 depending on
c          whether or not matrix values are provided.
c          Note that no error message is triggered (i.e. ierr = 0
c          on return in these cases. It is therefore important to
c          compare the values of job on entry and return ).
c
c a       = the a matrix in the a, ia, ja (column) storage format
c ja        = row number of element a(i,j) in array a.
c ia     = pointer  array. ia(i) points to the beginning of column i.
c
c rhs    = double precision array of size nrow + 1 if available (see job)
c
c nrhs   = integer containing the number of right hand sides found
c          each right hand side may be accompanied with an intial guess
c          and also the exact solution.
c
c guesol = a 2-character string indicating whether an initial guess
c          (1-st character) and / or the exact solution (2-nd
c          character) is provided with the right hand side.
c         if the first character of guesol is 'G' it means that an
c          an intial guess is provided for each right hand side.
c          These are appended to the right hand sides in the array rhs.
c         if the second character of guesol is 'X' it means that an
c          exact solution is provided for each right hand side.
c          These are  appended to the right hand sides
c          and the initial guesses (if any) in the array rhs.
c
c nrow   = number of rows in matrix
c ncol       = number of columns in matrix
c nnz       = number of nonzero elements in A. This info is returned
c          even if there is not enough space in a, ja, ia, in order
c          to determine the minimum storage needed.
c
c title  = character*72 = title of matrix test ( character a*72).
c key    = character*8  = key of matrix
c type   = charatcer*3  = type of matrix.
c          for meaning of title, key and type refer to documentation
c          Harwell/Boeing matrices.
c
c ierr   = integer used for error messages
c         * ierr  =  0 means that  the matrix has been read normally.
c         * ierr  =  1 means that  the array matrix could no be read
c         because ncol+1 .gt. nmax
c         * ierr  =  2 means that  the array matrix could no be read
c         because nnz .gt. nzmax
c         * ierr  =  3 means that  the array matrix could no be read
c         because both (ncol+1 .gt. nmax) and  (nnz .gt. nzmax )
c         * ierr  =  4 means that  the right hand side (s) initial
c         guesse (s) and exact solution (s)   could  not be
c         read because they are stored in sparse format (not handled
c         by this routine ...)
c         * ierr  =  5 means that the right hand sides, initial guesses
c         and exact solutions could not be read because the length of
c         rhs as specified by the input value of nrhs is not
c         insufficient to store them. The rest of the matrix may have
c         been read normally.
c
c Notes:
c
c 1) The file inout must be open (and possibly rewound if necessary)
c    prior to calling readmt.
c 2) Refer to the documentation on the Harwell-Boeing formats
c    for details on the format assumed by readmt.
c    We summarize the format here for convenience.
c
c    a) all lines in inout are assumed to be 80 character long.
c    b) the file consists of a header followed by the block of the
c       column start pointers followed by the block of the
c       row indices, followed by the block of the real values and
c       finally the numerical values of the right hand side if a
c       right hand side is supplied.
c    c) the file starts by a header which contains four lines if no
c       right hand side is supplied and five lines otherwise.
c       * first line contains the title (72 characters long) followed by
c         the 8-character identifier (name of the matrix, called key)
c        [ A72,A8 ]
c       * second line contains the number of lines for each
c         of the following data blocks (4 of them) and the total number
c         of lines excluding the header.
c        [5i4]
c       * the third line contains a three character string identifying
c         the type of matrices as they are referenced in the Harwell
c         Boeing documentation [e.g., rua, rsa,..] and the number of
c         rows, columns, nonzero entries.
c         [A3,11X,4I14]
c       * The fourth line contains the variable fortran format
c         for the following data blocks.
c         [2A16,2A20]
c       * The fifth line is only present if right hand sides are
c         supplied. It consists of three one character-strings containing
c         the storage format for the right hand sides
c         ('F'= full,'M'=sparse=same as matrix), an initial guess
c         indicator ('G' for yes), an exact solution indicator
c         ('X' for yes), followed by the number of right hand sides
c         and then the number of row indices.
c         [A3,11X,2I14]
c     d) The three following blocks follow the header as described
c        above.
c     e) In case the right hand side are in sparse formats then
c        the fourth block uses the same storage format as for the matrix
c        to describe the NRHS right hand sides provided, with a column
c        being replaced by a right hand side.
c
      character title*72, key*8, type*3, ptrfmt*16, indfmt*16,
     1       valfmt*20, rhsfmt*20, rhstyp*3, guesol*2
      integer totcrd, ptrcrd, indcrd, valcrd, rhscrd, nrow, ncol,
     1     nnz, neltvl, nrhs, nmax, nzmax
      integer ia (nmax+1), ja (nzmax)
      double precision a(nzmax), rhs(*)

      lenrhs = nrhs
c
      read(iounit,2010,end=10)title,key
      read(iounit,2011,end=10)totcrd,ptrcrd,indcrd,valcrd,rhscrd
      read(iounit,2012,end=10)type,nrow,ncol,nnz,neltvl
      read(iounit,2013,end=10)ptrfmt,indfmt,valfmt,rhsfmt
2010  format(a72,a8)
2011  format(5i14)
2012  format(a3,11x,4i14)
2013  format(2a16, 2a20)
c
      if (rhscrd .gt. 0) read (iounit,2014,end=10) rhstyp, nrhs
2014  format (a3,11x,i4)
c
c anything else to read ?
c
      if (job .le. 0) return
      ierr = 0
c  check whether matrix is readable.
      n = ncol
      if (ncol .gt. nmax) ierr = 1
      if (nnz .gt. nzmax) ierr = ierr + 2
      if (ierr .ne. 0) return
c  read pointer and row numbers.
      read (iounit,ptrfmt,end=10) (ia (i), i = 1, n+1)
      read (iounit,indfmt,end=10) (ja (i), i = 1, nnz)
c  reading values of matrix if required....
      if (job .le. 1)  return
c  and if available.
      if (valcrd .le. 0) then
       job = 1
       return
      end if
      read (iounit,valfmt,end=10) (a(i), i = 1, nnz)
c  reading rhs if required
      if (job .le. 2)  return
c  and if available
      if ( rhscrd .le. 0) then
       job = 2
       return
      end if
c
c  read right hand side.
c
      if (rhstyp(1:1) .eq. 'M') then
         ierr = 4
         return
      end if

      guesol = rhstyp(2:3)

      nvec = 1
      if (guesol(1:1) .eq. 'G') nvec=nvec+1
      if (guesol(2:2) .eq. 'X') nvec=nvec+1

      len = nrhs*nrow

      if (len*nvec .gt. lenrhs) then
         ierr = 5
         return
      end if
c
c read right hand sides
c
      next = 1
      iend = len
      read(iounit,rhsfmt,end=10) (rhs(i), i = next, iend)
c
c read initial guesses if available
c
      if (guesol(1:1) .eq. 'G') then
         next = next+len
         iend = iend+ len
         read(iounit,valfmt,end=10) (rhs(i), i = next, iend)
      end if
c
c read exact solutions if available
c
      if (guesol(2:2) .eq. 'X')then
         next = next+len
         iend = iend+ len
         read(iounit,valfmt,end=10) (rhs(i), i = next, iend)
      end if
c
      return
10    continue
      WRITE(*,*)' '
      WRITE(*,*)'READMT - Fatal error.'
      WRITE(*,*)'         End of file while reading information!'
      WRITE(*,*)'         Results are unreliable!'
      return
      end
      subroutine refall(nx, nelx,ijk,node,ndeg,x,y,
     *              ichild,iparnts,nodcode,nxmax,nelmax,ierr)

c*********************************************************************72
c
cc REFALL refines a finite element grid using triangular elements.
c
c REFALL refines a finite element grid using triangular elements.
c uses mid points to refine all the elements of the grid.
c
c nx      = number of nodes at input
c nelx      = number of elements at input
c ijk      = connectivity matrix: for node k, ijk(*,k) point to the
c         nodes of element k.
c ndeg      = first dimension of array ichild which is at least as large
c         as the max degree of each node
c x,y   = double precision arrays containing the x(*) and y(*) coordinates
c        resp. of the nodes.
c ichild= list of the children of a node: ichild(1,k) stores
c         the position in ichild(*,k)  of the last child so far.
c         (local use)
c iparnts= list of the 2 parents of each node.
c         (local use)
c nodcode= boundary information list for each node with the
c         following meaning:
c      nodcode(i) = 0 -->  node i is internal
c      nodcode(i) = 1 -->  node i is a boundary but not a corner point
c      nodcode(i) = 2 -->  node i is a corner point.
c corner elements are used only to generate the grid by refinement
c since they do not  correspond to real elements.
c nxmax  = maximum number of nodes allowed. If during the algorithm
c          the number of nodes being created exceeds nxmax then
c         refall  quits without modifying the (x,y) xoordinates
c         and nx, nelx. ijk is modified. Also ierr is set to 1.
c nelmax = same as above for number of elements allowed. See ierr..
c ierr       = error message:
c         0 --> normal return
c         1 --> refall quit because nxmax  was exceeded.
c         2 --> refall quit because nelmax was exceeded.
c
       implicit double precision  (a-h,o-z)
       integer ichild(ndeg,1),iparnts(2,nx),ijk(node,1),nodcode(nx)
       integer midnode(10),inod(10)
       double precision  x(1),y(1)
c
c inilitialize lists of children and parents --
c data structure is as follows
c ichild(1,k) stores the position of last child of node k so far in list
c ichild(j,k) , j .ge. 2 = list of children of node k.
c iparnts(1,k) and iparnts(2,k) are the two parents of node k.
c
c  do a first check :
      if (nx .ge. nxmax) goto 800
      if (nelx .ge. nelmax) goto 900
c  initialize
        do 1 k=1,nx
      do 2 j=2,ndeg
       ichild(j,k) = 0
 2      continue
      ichild(1,k) = 1
      iparnts(1,k)= 0
      iparnts(2,k)= 0
 1      continue
c  initialize nelxnew and nxnew
      nelxnew = nelx
      nxnew   = nx
      ierr    = 0
c
c main loop: scan all elements
c
c      do 100 nel = nelx,1,-1
      do 100 nel = 1, nelx
c note : interesting question which order is best for parallelism?
c alternative order: do 100 nel = nelx, 1, -1
c
c  unpack nodes of element
      do 101 i=1,node
      inod(i) = ijk(i,nel)
c convention: node after last node = first node.
      inod(node+i) = inod(i)
      midnode(i) = 0
 101      continue
c
c for each new potential node determine if it has already been
c numbered. a potential node is the middle of any two nodes ..
c
      do 80 ii=1,node
              k1 = inod(ii)
              k2 = inod(ii+1)
c  test for current pair :
      last = ichild(1,k1)
      do 21 k=2,last
            jchild = ichild(k,k1)
            ipar1 = iparnts(1,jchild)
            ipar2 = iparnts(2,jchild)
              if( (ipar1 .eq. k1 .and. ipar2 .eq. k2) .or.
     *                (ipar2 .eq. k1 .and. ipar1 .eq. k2)) then
c node has already been created and numbered ....
          midnode(ii) = jchild
c therefore it must be an internal node
          nodcode(jchild) = 0
c  and no new node to create.
          goto 80
              end if
 21      continue
c
c else  create a new node
c
      nxnew = nxnew + 1
      if (nxnew .gt. nxmax) goto 800

      x(nxnew) = (x(k1) + x(k2))*0.5
      y(nxnew) = (y(k1) + y(k2))*0.5
      midnode(ii) = nxnew
c
c update nodcode information -- normally min0(nodcode(k1),nodcode(k2))
c
       nodcode(nxnew) = min0(1,nodcode(k1),nodcode(k2))
c
c update parents and children's lists
c
      iparnts(1,nxnew) = k1
      iparnts(2,nxnew) = k2
c
      last = last+1
      ichild(last,k1) = nxnew
      ichild(1,k1) = last

      last = ichild(1,k2)+1
      ichild(last,k2) = nxnew
      ichild(1,k2) = last
c
 80     continue
c
c  replace current element by new one
c
      do 81 i=1,node
      jnod = midnode(i)
        ijk(i,nel) = jnod
 81     continue
c  create new elements
      do 82 ii=1, node
      nelxnew = nelxnew+1
      if (nelxnew .gt. nelmax) goto 900
      ijk(1,nelxnew) = inod(ii)
      k = ii
      do 82 jj=2,node
      ijk(jj,nelxnew) = midnode(k)
       k = k+2
      if (k .gt. node) k =  k-node
 82      continue
c  done !
 100      continue
      nx = nxnew
      nelx = nelxnew
        return
 800       ierr = 1
      return
 900    ierr = 2
       return
      end
      subroutine retmx(n,a,ja,ia,dd)

c*********************************************************************72
c
cc RETMX returns in dd(*) the max absolute value of elements in row *.
c
c RETMX returns in dd(*) the max absolute value of elements in row *.
c used for scaling purposes. superseded by rnrms  .
c
c on entry:
c n      = dimension of A
c a,ja,ia
c      = matrix stored in compressed sparse row format
c dd      = double precision array of length n. On output,entry dd(i) contains
c        the element of row i that has the largest absolute value.
c        Moreover the sign of dd is modified such that it is the
c        same as that of the diagonal element in row i.
c
c           Y. Saad, Sep. 21 1989                                      c
c
      double precision a(*),dd(*)
      integer n,ia(*),ja(*)
      integer k2, i, k1, k
      double precision t, t1, t2
c
c initialize
c
      k2 = 1
      do 11 i=1,n
         k1 = k2
         k2 = ia(i+1) - 1
         t = 0.0
         do 101  k=k1,k2
            t1 = abs(a(k))
            if (t1 .gt. t) t = t1
            if (ja(k) .eq. i)then
              if(a(k).lt.0.0)then
                t2=-1.0
              elseif(a(k).eq.0.0)then
                t2=0.0
              else
                t2=1.0
                end if
              end if
 101     continue
         dd(i) =  t2*t
c     we do not invert diag here..
 11   continue
      return
      end
      subroutine rnrms(nrow, nrm, a, ja, ia, diag)

c*********************************************************************72
c
cc RNRMS gets the norms of each row of A. (choice of three norms)
c
c on entry:
c
c nrow      = integer. The row dimension of A
c
c nrm   = integer. norm indicator. nrm = 1, means 1-norm, nrm =2
c                  means the 2-nrm, nrm = 0 means max norm
c
c a,
c ja,
c ia   = Matrix A in compressed sparse row format.
c
c on return:
c
c
c diag = double precision vector of length nrow containing the norms
c
c
      double precision a(*), diag(nrow), scal
      integer ja(*), ia(nrow+1)
c
      do 1 ii=1,nrow
c
c     compute the norm if each element.
c
         scal = 0.0
         k1 = ia(ii)
         k2 = ia(ii+1)-1
         if (nrm .eq. 0) then
            do 2 k=k1, k2
               scal = max(scal,abs(a(k) ) )
 2          continue
         elseif (nrm .eq. 1) then
            do 3 k=k1, k2
               scal = scal + abs(a(k) )
 3          continue
         else
            do 4 k=k1, k2
               scal = scal+a(k)**2
 4          continue
         end if
         if (nrm .eq. 2) scal = sqrt(scal)
         diag(ii) = scal
 1    continue
      return
      end
      subroutine rperm (nrow,a,ja,ia,ao,jao,iao,perm,job)

c*********************************************************************72
c
cc RPERM permutes the rows of a matrix in CSR format.
c
c rperm  computes B = P A  where P is a permutation matrix.
c the permutation P is defined through the array perm: for each j,
c perm(j) represents the destination row number of row number j.
c Youcef Saad -- recoded Jan 28, 1991.
c
c on entry:
c
c n       = dimension of the matrix
c a, ja, ia = input matrix in csr format
c perm       = integer array of length nrow containing the permutation arrays
c        for the rows: perm(i) is the destination of row i in the
c         permuted matrix.
c         ---> a(i,j) in the original matrix becomes a(perm(i),j)
c         in the output  matrix.
c
c job      = integer indicating the work to be done:
c             job = 1      permute a, ja, ia into ao, jao, iao
c                       (including the copying of real values ao and
c                       the array iao).
c             job .ne. 1 :  ignore real values.
c                     (in which case arrays a and ao are not needed nor
c                      used).
c
c
c on return:
c
c ao, jao, iao = input matrix in a, ja, ia format
c note :
c        if (job.ne.1)  then the arrays a and ao are not used.
c
c           Y. Saad, May  2, 1990                                      c
c
      integer nrow,ja(*),ia(nrow+1),jao(*),iao(nrow+1),perm(nrow),job
      double precision a(*),ao(*)
c
      logical values
      values = (job .eq. 1)
c
c     determine pointers for output matix.
c
      do j=1,nrow
         i = perm(j)
         iao(i+1) = ia(j+1) - ia(j)
      end do
c
c get pointers from lengths
c
      iao(1) = 1
      do 51 j=1,nrow
         iao(j+1)=iao(j+1)+iao(j)
 51   continue
c
c copying
c
      do 100 ii=1,nrow
c
c old row = ii  -- new row = iperm(ii) -- ko = new pointer
c
         ko = iao(perm(ii))
         do 60 k=ia(ii), ia(ii+1)-1
            jao(ko) = ja(k)
            if (values) ao(ko) = a(k)
            ko = ko+1
 60      continue
 100  continue
c
      return
      end
      subroutine rscal(nrow, job, nrm, a, ja, ia, diag, b, jb, ib)

c*********************************************************************72
c
cc RSCAL normalizes the rows of A.
c
c scales the rows of A such that their norms are one on return
c 3 choices of norms: 1-norm, 2-norm, max-norm.
c
c on entry:
c
c nrow      = integer. The row dimension of A
c
c job   = integer. job indicator. Job=0 means get array b only
c         job = 1 means get b, and the integer arrays ib, jb.
c
c nrm   = integer. norm indicator. nrm = 1, means 1-norm, nrm =2
c                  means the 2-nrm, nrm = 0 means max norm
c
c a,
c ja,
c ia   = Matrix A in compressed sparse row format.
c
c on return:
c
c
c diag = diagonal matrix stored as a vector containing the matrix
c        by which the rows have been scaled, i.e., on return
c        we have B = Diag*A.
c
c b,
c jb,
c ib      = resulting matrix B in compressed sparse row sparse format.
c
c Notes:
c
c 1)        The column dimension of A is not needed.
c 2)        algorithm in place (B can take the place of A).
c
      double precision a(*), b(*), diag(nrow)
      integer ja(*),jb(*), ia(nrow+1),ib(nrow+1)

      call rnrms (nrow,nrm,a,ja,ia,diag)

      do j=1, nrow
         diag(j) = 1.0/diag(j)
      end do

      call diamua(nrow,job,a,ja,ia,diag,b,jb,ib)
      return
      end
      subroutine sskssr (n,imod,asky,isky,ao,jao,iao,nzmax,ierr)

c*********************************************************************72
c
cc SSKSSR converts Symmetric Skyline Format to Symmetric Sparse Row format.
c
c  tests for exact zeros in skyline matrix (and ignores them in
c  output matrix).  In place routine (a, isky :: ao, iao)
c
c this subroutine translates a  symmetric skyline format into a
c symmetric sparse row format. Each element is tested to see if it is
c a zero element. Only the actual nonzero elements are retained. Note
c that the test used is simple and does take into account the smallness
c of a value. the subroutine filter (see unary module) can be used
c for this purpose.
c
c Coded by Y. Saad, Oct 5, 1989. Revised Feb 18, 1991./
c
c on entry:
c
c n      = integer equal to the dimension of A.
c imod  = integer indicating the variant of skyline format used:
c         imod = 0 means the pointer iao points to the `zeroth'
c         element of the row, i.e., to the position of the diagonal
c         element of previous row (for i=1, iao(1)= 0)
c         imod = 1 means that itpr points to the beginning of the row.
c         imod = 2 means that iao points to the end of the row
c                  (diagonal element)
c asky  = double precision array containing the values of the matrix. asky contains
c         the sequence of active rows from i=1, to n, an active row
c         being the row of elemnts of the matrix contained between the
c         leftmost nonzero element and the diagonal element.
c isky       = integer array of size n+1 containing the pointer array to
c         each row. isky (k) contains the address of the beginning of the
c         k-th active row in the array asky.
c nzmax = integer. equal to the number of available locations in the
c         output array ao.
c
c on return:
c
c ao      = double precision array of size nna containing the nonzero elements
c jao      = integer array of size      nnz containing the column positions
c         of the corresponding elements in a.
c iao      = integer of size n+1. iao(k) contains the position in a, ja of
c        the beginning of the k-th row.
c ierr  = integer. Serving as error message. If the length of the
c         output arrays ao, jao exceeds nzmax then ierr returns
c         the row number where the algorithm stopped: rows
c         i, to ierr-1 have been processed succesfully.
c         ierr = 0 means normal return.
c         ierr = -1  : illegal value for imod
c Notes:
c
c This module is in place: ao and iao can be the same as asky, and isky.
c
      INTEGER NZMAX
      double precision asky(*),ao(nzmax)
      integer n, imod,ierr, isky(n+1),iao(n+1),jao(nzmax)
      integer next, kend, kstart, i, j
      ierr = 0
c
c check for validity of imod
c
      if (imod.ne.0 .and. imod.ne.1 .and. imod .ne. 2) then
         ierr =-1
         return
      end if
c
c next  = pointer to next available position in output matrix
c kend  = pointer to end of current row in skyline matrix.
c
      next = 1
c
c set kend = start position -1 in  skyline matrix.
c
      kend = 0
      if (imod .eq. 1) kend = isky(1)-1
      if (imod .eq. 0) kend = isky(1)
c
c loop through all rows
c
      do 50 i=1,n
c
c save value of pointer to ith row in output matrix
c
         iao(i) = next
c
c get beginnning and end of skyline  row
c
         kstart = kend+1
         if (imod .eq. 0) kend = isky(i+1)
         if (imod .eq. 1) kend = isky(i+1)-1
         if (imod .eq. 2) kend = isky(i)
c
c copy element into output matrix unless it is a zero element.
c
         do 40 k=kstart,kend
            if (asky(k) .eq. 0.0) goto 40
            j = i-(kend-k)
            jao(next) = j
            ao(next)  = asky(k)
            next=next+1
            if (next .gt. nzmax+1) then
               ierr = i
               return
            end if
 40      continue
 50    continue
      iao(n+1) = next
      return
      end
      subroutine ssrcsr (nrow,a,ja,ia,nzmax,ao,jao,iao,indu,ierr)

c*********************************************************************72
c
cc SSRCSR converts Symmetric Sparse Row to (regular) Compressed Sparse Row.
c
c this subroutine converts  a symmetric  matrix in which only the lower
c part is  stored in compressed sparse row format, i.e.,
c a matrix stored in symmetric sparse format, into a fully stored matrix
c i.e., a matrix where both the lower and upper parts are stored in
c compressed sparse row format. the algorithm is in place (i.e. result
c may be overwritten onto the input matrix a, ja, ia ).
c the output matrix delivered by ssrcsr is such that each row starts with
c the elements of the lower part followed by those of the upper part.
c
c on entry:
c
c
c nrow  = row dimension of inout matrix
c a,
c ia,
c ja    = matrix in compressed sparse row format. This is assumed to be
c         a lower triangular matrix.
c
c nzmax      = size of arrays ao and jao. ssrcsr will abort if the storage
c         provided in a, ja is not sufficient to store A. See ierr.
c
c on return:
c
c ao, iao,
c   jao = output matrix in compressed sparse row format. The resulting
c         matrix is symmetric and is equal to A+A**T - D, if
c         A is the original lower triangular matrix. ao, jao, iao,
c         can be the same as a, ja, ia in the calling sequence.
c
c indu  = integer array of length nrow+1. If the input matrix is such
c         that the last element in each row is its diagonal element then
c         on return, indu will contain the pointers to the diagonal
c         element in each row of the output matrix. Otherwise used as
c         work array.
c ierr  = integer. Serving as error message. If the length of the arrays
c         ao, jao exceeds nzmax, ierr returns the minimum value
c         needed for nzmax. otherwise ierr=0 (normal return).
c
      integer ia(nrow+1),iao(nrow+1),ja(*),jao(nzmax),indu(nrow+1)
      double precision a(*),ao(nzmax)

      ierr = 0
      do 1 i=1,nrow+1
         indu(i) = 0
 1    continue
c
c     compute  number of elements in each row of strict upper part.
c     put result in indu(i+1)  for row i.
c
      do 3 i=1, nrow
         do 2 k=ia(i),ia(i+1)-1
            j = ja(k)
            if (j .lt. i) indu(j+1) = indu(j+1)+1
 2       continue
 3    continue
c
c     find addresses of first elements of ouput matrix. result in indu
c
      indu(1) = 1
      do 4 i=1,nrow
         lenrow = ia(i+1)-ia(i)
         indu(i+1) = indu(i) + indu(i+1) + lenrow
 4    continue
c  enough storage in a, ja ?
c
      nnz = indu(nrow+1)-1
      if (nnz .gt. nzmax) then
         ierr = nnz
         return
      end if
c
c     now copy lower part (backwards).
c
      kosav = indu(nrow+1)
      do 6 i=nrow,1,-1
         klast = ia(i+1)-1
         kfirst = ia(i)
         iao(i+1) = kosav
         ko = indu(i)
         kosav = ko
         do 5 k = kfirst, klast
            ao(ko) = a(k)
            jao(ko) = ja(k)
          ko = ko+1
 5       continue
         indu(i) = ko
 6    continue
      iao(1) = 1
c
c     now copy upper part. Go through the structure of ao, jao, iao
c     that has already been copied (lower part). indu(i) is the address
c     of the next free location in row i for ao, jao.
c
      do 8 i=1,nrow
c     i-th row is now in ao, jao, iao structure: lower half part
         do 9 k=iao(i), iao(i+1)-1
            j = jao(k)
            if (j .ge. i)  goto 8
            ipos = indu(j)
            ao(ipos) = ao(k)
            jao(ipos) = i
            indu(j) = indu(j) + 1
 9       continue
 8    continue
      return
      end
      subroutine submat (n,job,i1,i2,j1,j2,a,ja,ia,nr,nc,ao,jao,iao)

c*********************************************************************72
c
cc SUBMAT extracts the submatrix A(i1:i2,j1:j2).
c
c extracts the submatrix A(i1:i2,j1:j2) and puts the result in
c matrix ao,iao,jao
c  In place: ao,jao,iao may be the same as a,ja,ia.
c
c on input
c
c n      = row dimension of the matrix
c i1,i2 = two integers with i2 .ge. i1 indicating the range of rows to be
c          extracted.
c j1,j2 = two integers with j2 .ge. j1 indicating the range of columns
c         to be extracted.
c         * There is no checking whether the input values for i1, i2, j1,
c           j2 are between 1 and n.
c a,
c ja,
c ia    = matrix in compressed sparse row format.
c
c job      = job indicator: if job .ne. 1 then the real values in a are NOT
c         extracted, only the column indices (i.e. data structure) are.
c         otherwise values as well as column indices are extracted...
c
c on output
c
c nr      = number of rows of submatrix
c nc      = number of columns of submatrix
c        * if either of nr or nc is nonpositive the code will quit.
c
c ao,
c jao,iao = extracted matrix in general sparse format with jao containing
c      the column indices,and iao being the pointer to the beginning
c      of the row,in arrays a,ja.
c
c           Y. Saad, Sep. 21 1989                                      c
c
      integer n,job,i1,i2,j1,j2,nr,nc,ia(*),ja(*),jao(*),iao(*)
      double precision a(*),ao(*)

      nr = i2-i1+1
      nc = j2-j1+1
c
      if ( nr .le. 0 .or. nc .le. 0) return
c
      klen = 0
c
c     simple procedure that proceeds row-wise...
c
      do 100 i = 1,nr
         ii = i1+i-1
         k1 = ia(ii)
         k2 = ia(ii+1)-1
         iao(i) = klen+1

         do k=k1,k2
            j = ja(k)
            if (j .ge. j1 .and. j .le. j2) then
               klen = klen+1
               if (job .eq. 1) ao(klen) = a(k)
               jao(klen) =j
            end if
         end do

 100  continue
      iao(nr+1) = klen+1
      return
      end
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
      subroutine transp(nrow,ncol,a,ja,ia,iwk,ierr)

c*********************************************************************72
c
cc TRANSP carries out in-place transposition routine.
c
c this subroutine transposes a matrix stored in compressed sparse row
c format. the transposition is done in place in that the arrays a,ja,ia
c of the transpose are overwritten onto the original arrays.
c
c on entry:
c
c nrow      = integer. The row dimension of A.
c ncol      = integer. The column dimension of A.
c a      = double precision array of size nnz (number of nonzero elements in A).
c         containing the nonzero elements
c ja      = integer array of length nnz containing the column positions
c         of the corresponding elements in a.
c ia      = integer of size n+1, where n = max(nrow,ncol). On entry
c         ia(k) contains the position in a,ja of  the beginning of
c         the k-th row.
c
c iwk      = integer work array of same length as ja.
c
c on return:
c
c
c ncol      = actual row dimension of the transpose of the input matrix.
c         Note that this may be .le. the input value for ncol, in
c         case some of the last columns of the input matrix are zero
c         columns. In the case where the actual number of rows found
c         in transp(A) exceeds the input value of ncol, transp will
c         return without completing the transposition. see ierr.
c a,
c ja,
c ia      = contains the transposed matrix in compressed sparse
c         row format. The row dimension of a, ja, ia is now ncol.
c
c ierr      = integer. error message. If the number of rows for the
c         transposed matrix exceeds the input value of ncol,
c         then ierr is  set to that number and transp quits.
c         Otherwise ierr is set to 0 (normal return).
c
c Note:
c      1) If you do not need the transposition to be done in place
c         it is preferrable to use the conversion routine csrcsc
c         (see conversion routines in formats).
c      2) the entries of the output matrix are not sorted (the column
c         indices in each are not in increasing order) use csrcsc
c         if you want them sorted.
c
c           Y. Saad, Sep. 21 1989                                      c
c  modified Oct. 11, 1989.                                             c
c
      integer nrow, ncol, ia(*), ja(*), iwk(*), ierr
      double precision a(*)
      double precision t, t1
      ierr = 0
      nnz = ia(nrow+1)-1
c
c     determine column dimension
c
      jcol = 0

      do k=1, nnz
         jcol = max(jcol,ja(k))
      end do

      if (jcol .gt. ncol) then
         ierr = jcol
         return
      end if
c
c     convert to coordinate format. use iwk for row indices.
c
      ncol = jcol
c
      do 3 i=1,nrow
         do 2 k=ia(i),ia(i+1)-1
            iwk(k) = i
 2       continue
 3    continue
c     find pointer array for transpose.
      do 35 i=1,ncol+1
         ia(i) = 0
 35   continue
      do 4 k=1,nnz
         i = ja(k)
         ia(i+1) = ia(i+1)+1
 4    continue
      ia(1) = 1

      do i=1,ncol
         ia(i+1) = ia(i) + ia(i+1)
      end do
c
c  loop for a cycle in chasing process.
c
      init = 1
      k = 0
 5    t = a(init)
      i = ja(init)
      j = iwk(init)
      iwk(init) = -1

 6    k = k+1
c
c     current row number is i.  determine  where to go.
c
      l = ia(i)
c
c     save the chased element.
c
      t1 = a(l)
      inext = ja(l)
c
c     then occupy its location.
c
      a(l)  = t
      ja(l) = j
c
c     update pointer information for next element to be put in row i.
c
      ia(i) = l+1
c
c     determine  next element to be chased
c
      if (iwk(l) .lt. 0) goto 65
      t = t1
      i = inext
      j = iwk(l)
      iwk(l) = -1
      if (k .lt. nnz) goto 6
      goto 70
 65   init = init+1
      if (init .gt. nnz) goto 70
      if (iwk(init) .lt. 0) goto 65
c
c     restart chasing
c
      goto 5
 70   continue
      do 80 i=ncol,1,-1
         ia(i+1) = ia(i)
 80   continue
      ia(1) = 1

      return
      end
      subroutine udsol (n,x,y,au,jau)

c*********************************************************************72
c
cc UDSOL solves U*x = y;   U = upper triangular in MSR format
c
c solves a non-unit upper triangular matrix by standard (sequential )
c backward elimination - matrix stored in MSR format.
c with diagonal elements already inverted (otherwise do inversion,
c au(1:n) = 1.0/au(1:n),  before calling).
c
c On entry:
c
c n      = integer. dimension of problem.
c y      = double precision array containg the right side.
c
c au,
c jau,    = Lower triangular matrix stored in modified sparse row
c          format.
c
c On return:
c
c      x = The solution of  U x = y .
c
      integer n, jau(*)
      double precision  x(n), y(n),au(*)
      integer k, j
      double precision t

      x(n) = y(n)*au(n)
      do k = n-1,1,-1
         t = y(k)
         do j = jau(k), jau(k+1)-1
            t = t - au(j)*x(jau(j))
         end do
         x(k) = au(k)*t
      end do

      return
      end
      subroutine udsolc (n,x,y,au,jau)

c*********************************************************************72
c
cc UDSOLC solves U * x = y, for upper triangular U in MSC format.
c
c solves a (non-unit) upper triangular system by standard (sequential)
c forward elimination - matrix stored in Modified Sparse Column format
c with diagonal elements already inverted (otherwise do inversion,
c auuuul(1:n) = 1.0/au(1:n),  before calling ldsol).
c
c
c On entry:
c
c n      = integer. dimension of problem.
c y      = double precision array containg the right hand side.
c
c au,
c jau,   = Upper triangular matrix stored in Modified Sparse Column
c          format.
c
c On return:
c
c      x = The solution of  U x = y .
c
      integer n, jau(*)
      double precision x(n), y(n), au(*)
      integer k, j
      double precision t

      do k=1,n
         x(k) = y(k)
      end do

      do 150 k = n,1,-1
         x(k) = x(k)*au(k)
         t = x(k)
         do 100 j = jau(k), jau(k+1)-1
            x(jau(j)) = x(jau(j)) - t*au(j)
 100     continue
 150  continue
c
      return
      end
      subroutine unassbl (a,na,f,nx,nelx,ijk,nodcode,
     *                     node,x,y,ierr,xyk)

c*********************************************************************72
c
cc UNASSBL ?
c
c a      = un-assembled matrix on output
c na       = 1-st dimension of a.  a(na,node,node)
c
c f      = right hand side (global load vector) in un-assembled form
c nx     = number of nodes at input
c nelx       = number of elements at input
c ijk       = connectivity matrix: for node k, ijk(*,k) point to the
c          nodes of element k.
c node       = total number of nodal points in each element
c         also second dimension of a.
c
c nodcode= boundary information list for each node with the
c         following meaning:
c      nodcode(i) = 0 -->  node i is internal
c      nodcode(i) = 1 -->  node i is a boundary but not a corner point
c      nodcode(i) = 2 -->  node i is a corner point (corner points
c
c x,y   = double precision arrays containing the $x$ and $y$ coordinates
c        resp. of the nodes.
c         K11, K22, and K12 at that element.
c ierr      = error message integer .
c        ierr = 0 --> normal return
c        ierr = 1 --> negative area encountered (due to bad
c                 numbering of nodes of an element-
c               message printed in unit iout).
c iout      = output unit (not used here).
c
c xyk      = subroutine defining the material properties at each
c         element. Form:
c       call xyk(nel,xyke,x,y,ijk,node)
c
      implicit double precision (a-h,o-z)
        dimension a(na,node,node),ijk(node,1),x(1),y(1),f(node,1),
     *          ske(3,3),fe(3),xe(3),ye(3),xyke(2,2)
            integer nodcode(1)
       external xyk
c max number of nonzeros allowed  = 200
c
c   initialize
c
      do i=1, node
        do j=1, nx
          f(i,j) = 0.0
        end do
      end do
c
c main loop
c
      do 102 nel=1, nelx
c
c get coordinetes of nodal points
c
      do 104 i=1, node
      j = ijk(i,nel)
      xe(i) = x(j)
      ye(i) = y(j)
 104      continue
c
c compute determinant
c
       det=xe(2)*(ye(3)-ye(1))+xe(3)*(ye(1)-ye(2))+xe(1)*(ye(2)-ye(3))
c        print *, 'nel', nel, ' det = ' , del
c
c set material properties
c
      call xyk(nel,xyke,x,y,ijk,node)
c
c construct element stiffness matrix
c
      ierr = 0
      call estif3(nel,ske,fe,det,xe,ye,xyke,ierr)
      if (ierr .ne. 0) return
c      write (8,'(9f8.4)') ((ske(i,j),j=1,3),i=1,3)
c assemble: add element stiffness matrix to global matrix
c
      do 120 ka=1, node
            f(ka,nel) = fe(ka)
        do 108 kb = 1,node
            a(nel,ka,kb) = ske(ka,kb)
 108      continue
 120      continue
 102      continue
        return
      end
      subroutine usol (n,x,y,au,jau,iau)

c*********************************************************************72
c
cc USOL solves   U x = y    U = unit upper triangular.
c
c solves a unit upper triangular system by standard (sequential )
c backward elimination - matrix stored in CSR format.
c
c
c On entry:
c
c n      = integer. dimension of problem.
c y      = double precision array containg the right side.
c
c au,
c jau,
c iau,    = Lower triangular matrix stored in compressed sparse row
c          format.
c
c On return:
c
c      x = The solution of  U x = y .
c
      integer n, jau(*),iau(n+1)
      double precision  x(n), y(n), au(*)
      integer k, j
      double precision  t

      x(n) = y(n)
      do 150 k = n-1,1,-1
         t = y(k)
         do 100 j = iau(k), iau(k+1)-1
            t = t - au(j)*x(jau(j))
 100     continue
         x(k) = t
 150  continue

      return
      end
      subroutine usolc (n,x,y,au,jau,iau)

c*********************************************************************72
c
cc USOLC solves U * x = y for unit upper triangular U in CSC format.
c
c solves a unit upper triangular system by standard (sequential )
c forward elimination - matrix stored in CSC format.
c
c
c On entry:
c
c n      = integer. dimension of problem.
c y      = double precision array containg the right side.
c
c au,
c jau,
c iau,    = Uower triangular matrix stored in compressed sparse column
c          format.
c
c On return:
c
c      x  = The solution of  U x  = y.
c
      double precision  x(*), y(*), au(*)
      integer n, jau(*),iau(*)
      integer k, j
      double precision t
c
      do 140 k=1,n
         x(k) = y(k)
 140  continue
      do 150 k = n,1,-1
         t = x(k)
         do 100 j = iau(k), iau(k+1)-1
            x(jau(j)) = x(jau(j)) - t*au(j)
 100     continue
 150  continue

      return
      end
