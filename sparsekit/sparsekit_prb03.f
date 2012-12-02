      program main

c*********************************************************************72
c
cc MAIN is a test suite for the Formats routines.
c
c tests all of the routines in the module formats.
c
c Note: the comments may not have been updated.
c
c Here is the sequence of what is done by this program.
c 1) generates a block matrix associated with a simple 5-point
c    matrix on a 4 x 2 grid (2-D) with 2 degrees of freedom per
c    grid point. Thus N = 16. This is produced in block format
c    using the generation routine genpbl.
c 2) the block format is translated into a compressed sparse row
c    format by bsrcsr. The result is dumped in file csr.mat
c 3) the matrix is translated in dense format by csrdns.
c    result a 16 x 16 matrix is written in unit dns.mat.
c    This is a good file to look at to see what the matrix is
c    and to compare results of other formats with.
c 4) the dense matrix obtained in 3) is reconverted back to
c    csr format using dnscsr. Result appended to file csr.mat
c 5) The matrix obtained in 4) is converted in coordinate format
c    and the resulting matrix is written in file coo.mat
c 6) the result is converted back to csr format. matrix
c    appended to csr.mat.
c 7) result of 6) is converted to symmetric sparse row storage
c    (ssr) and the result is appended to csr.mat
c 8) result of 7) converted back to csr format and result is
c    appended to csr.mat
c 9) matrix resulting from 8) is converted to modified sparse
c    row format using csrmsr and result is written in msr.mat.
c10) the resulting matrix is converted back to csrformat and
c    result is appended to csr.mat
c11) result is converted to ellpack-itpack format with
c    csrell and result is printed in ell.mat
c12) result is converted back to csr format and appended to csr.mat
c12) result converted to csc format (transposition) using csrcsc
c    which should produce the same matrix here. result appended
c    to csr.mat. A second call to csrcsc is made on resulting
c    matrix.
c13) the subroutine csrdia is used to extract two diagonals
c    (offsets -1 and 0) and then all the diagonals of matrix.
c    results in dia.mat
c14) diacsr is then called to convert the diagonally stored matrix
c    back to csr format. result appended to csr.mat
c15) result is converted to band format (bnd) by calling
c    csrbnd. result dumped to bnd.mat
c16) result is converted back to csr format and appended to csr.mat
c17) result sorted by a call to csrcsc and then converted to
c    block format (csrbsr) and then back to csr format again.
c    result appedned to csr.mat.
c18) matrix converted to symmetric skyline format. result appended
c    to file band.mat
c19) matrix converted back to csr format and result appended to
c    csr.mat.
c20) result converted to jad format. result output in jad.mat
c21) result concverted back to csr fromat. appended to csr.mat
c
      implicit none

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
      integer idiag
      integer idiag0
      integer ierr
      integer imod
      integer ioff(20)
      integer iout
      integer iwk(nmx*2+1)
      integer j
      integer ja(nnzmax)
      integer ja1(nnzmax)
      integer job
      integer k
      integer k1
      integer k2
      integer kdiag
      integer kend
      integer kstart
      integer len
      integer lowd
      integer maxcol
      integer ml
      integer mu
      integer n
      integer na
      integer ndiag
      integer nel
      integer nfree
      integer nnz
      integer nx
      integer ny
      integer nz
      double precision stencil(7,100)
      double precision wk(nmx)

      call timestamp ( )

      WRITE(*,*)' '
      WRITE(*,*)'SPARSEKIT_PRB03'
      write ( *, '(a)' ) '  FORTRAN77 version'
      write ( *, '(a)' ) ' '
      WRITE(*,*)'A set of tests for SPARSKIT'
      WRITE(*,*)' '

      open (unit=7,file='csr.mat',STATUS='replace')
      open (unit=8,file='dns.mat',STATUS='replace')
      open (unit=9,file='coo.mat',STATUS='replace')
      open (unit=10,file='msr.mat',STATUS='replace')
      open (unit=11,file='ell.mat',STATUS='replace')
      open (unit=12,file='dia.mat',STATUS='replace')
      open (unit=13,file='bnd.mat',STATUS='replace')
      open (unit=14,file='jad.mat',STATUS='replace')
c
c  dimension of grid
c
      nx = 4
      ny = 2
      nz = 1
      nfree = 2
c
c  generate grid problem.
c
      na = nx*ny*nz*5
      call gen57bl (nx,ny,nz,nfree,na,n,a1,ja1,ia1,iwk,stencil)
c
c  write out the matrix
c
      call bsrcsr (n,nfree,na,a1,ja1,ia1,a,ja,ia)
      iout = 7
      nnz = ia(n+1)-1
      write (iout,*) '-----------------------------------------'
      write (iout,*) '  +++  initial matrix in CSR format +++ '
      write (iout,*) '-----------------------------------------'
        call dump (n,a,ja,ia,6)
c
c  call csrdns
c
       call csrdns(n,n,a,ja,ia,dns,ndns,ierr)
       iout = iout+1
       write (iout,*) '-----------------------------------------'
         write (iout,*) '  +++  initial matrix in DENSE format+++ '
       write (iout,*) '-----------------------------------------'
       write (iout,'(4x,16i4)') (j,j=1,n)
       write (iout,'(3x,65(1h-))')
       do i=1, n
          write (8,102) i,(dns(i,j), j=1,n)
 102          format(1h ,i2,1h|,16f4.1)
       end do
c
c  convert it back to sparse format.
c
        call dnscsr(n,n,nnzmax,dns,ndns,a1,ja1,ia1,ierr)
          write(*,*) '-----------------------------------------'
        write(*,*) '  +++ matrix after conversion from dnscsr +++ '
         write(*,*) '-----------------------------------------'
        if (ierr .ne. 0) write(*,*) ' ***** ERROR FROM DNSCSR'
        if (ierr .ne. 0) write(*,*)  '     IERR = ', ierr
        call dump (n,a1,ja1,ia1,6)
c
c  convert it to coordinate format.
c
          call csrcoo(n,3,nnzmax,a,ja,ia,nnz,a1,ia1,ja1,ierr)
        iout = iout+ 1
        if (ierr .ne. 0) write (iout,*) ' ***** ERROR IN CSRCOO'
        if (ierr .ne. 0) write (iout,*)   '     IERR = ', ierr
          write (iout,*) '-----------------------------------------'
        write(iout,*) ' +++ Matrix in coordinate format +++ '
         write (iout,*) '-----------------------------------------'
          write(iout,103) (ia1(j),ja1(j),a1(j),j=1,nnz)
 103      format (' i =', i3,'    j = ',i3,'     a(i,j) = ',f4.1)
c
c  convert it back again to csr format
c
          call coocsr(n,nnz,a1,ia1,ja1,a,ja,ia)

         write(*,*) '-----------------------------------------'
        write(*,*) '  +++ matrix after conversion from coocsr +++ '
         write(*,*) '-----------------------------------------'
          call dump(n,a,ja,ia,6)
c
c  going to srs format
c
          call csrssr(n,a,ja,ia,nnzmax,a1,ja1,ia1,ierr)
         write(*,*) '-----------------------------------------'
        write(*,*) '  +++ matrix after conversion to ssr format +++ '
        write(*,*) '      (lower part only stored in csr format)    '
         write(*,*) '-----------------------------------------'
         call dump(n,a1,ja1,ia1,6)
c
c  back to csr
c
          call ssrcsr (n,a1,ja1,ia1,nnzmax,a,ja,ia,iwk,ierr)
          if (ierr .ne. 0) write(7,*) ' error in ssrcsr-IERR=',ierr
         write(*,*) '-----------------------------------------'
        write(*,*) '  +++ matrix after conversion from ssrcsr +++ '
         write(*,*) '-----------------------------------------'
         call dump(n,a,ja,ia,6)
c
c  msr format
c
         iout = iout+1
         call csrmsr (n,a,ja,ia,a1,ja1,a1,ja1)
        write (iout,*) '-----------------------------------------'
         write (iout,*) '  +++ matrix in modified sparse row format +++'
        write (iout,*) '-----------------------------------------'
         write (iout,*) ' ** MAIN DIAGONAL '
         write (iout,'(16f4.1)') (a1(k),k=1,n)
            write (iout,*) ' ** POINTERS: '
         write (iout,'(17i4)') (ja1(k),k=1,n+1)
             write (iout,*) ' ** REMAINDER :'
         call dump(n,a1,ja1,ja1,iout)

         call msrcsr (n,a1,ja1,a,ja,ia,wk)
        write(*,*) '-----------------------------------------'
        write(*,*) ' +++ matrix after conversion from msrcsr +++'
        write(*,*) '-----------------------------------------'

           call dump(n,a,ja,ia,6)

        maxcol = 13
c
        call csrell (n,a,ja,ia,maxcol,a1,ja1,n,ndiag,ierr)
        iout = iout+1
        if (ierr .ne. 0) write (iout,*) ' ***** ERROR IN CSRELL'
        if (ierr .ne. 0) write (iout,*)   '     IERR = ', ierr
        write (iout,*) '-----------------------------------------'
        write (iout,*) '  +++ matrix in ELLPACK-ITPACK format +++ '
        write (iout,*) '-----------------------------------------'
         do i=1,ndiag
         write (iout,*) ' Column number: ', i
         write (iout,104) (a1(n*(i-1)+k),k=1,n)
 104         format(9h COEF  = ,16f4.0)
         write (iout,105) (ja1(n*(i-1)+k),k=1,n)
 105         format (9h JCOEF = ,16i4)
        end do
        call ellcsr (n,a1,ja1,n,ndiag,a,ja,ia,nnzmax,ierr)
        if (ierr .ne. 0) write(*,*) ' ***** ERROR IN ELLCSR'
        if (ierr .ne. 0) write(*,*)   '     IERR = ', ierr
        write(*,*) '-----------------------------------------'
         write(*,*) '  +++ matrix after conversion from ellcsr +++'
        write(*,*) '-----------------------------------------'
         call dump(n,a,ja,ia,6)

       call csrcsc(n,1,1,a,ja,ia, a1,ja1,ia1)
        write(*,*) '-----------------------------------------'
         write(*,*) '  +++ matrix after conversion from csrcsc  +++ '
        write(*,*) '-----------------------------------------'
       call dump(n,a1,ja1,ia1,6)
       call csrcsc(n,1,1,a1,ja1,ia1, a,ja,ia)
c
c  test 1:
c  get main diagonal and subdiagonal
c  get some info on diagonals
c
       call infdia(n,ja,ia,iwk,idiag0)
       job = 0
       ioff(1) = 0
       ioff(2) = -1
       idiag = 2
       call csrdia (n,idiag,job,a,ja,ia,ndns,dns,ioff,a1,ja1,ia1,iwk)
       iout = iout+1
              write (iout,*) '-----------------------------------------'
         write (iout,*) '  +++  diagonal format +++ '
        write (iout,*) '-----------------------------------------'
         write (iout,*) '  diagonals ioff = 0 and ioff = -1 '
       write (iout,*) ' number of diag.s returned from csrdia=',idiag
       do kdiag = 1, idiag
         write (iout,*) ' diagonal offset = ', ioff(kdiag)
         write (iout,'(16f4.1)') (dns(k,kdiag),k=1,n)
       end do
c
c reverse conversion
c
       ndiag = ndns
       idiag = idiag0
       job = 10
       call csrdia (n,idiag,job,a,ja,ia,ndns,dns,ioff,a1,ja1,ia1,iwk)
       write (iout,*) '-----------------------------------------'
         write (iout,*) '  +++  second test diagonal format +++ '
         write (iout,*) '         ** all diagonals of A  ** '
       write (iout,*) '-----------------------------------------'
       write (iout,*) ' number of diagonals on return from csrdia=', 
     &                 idiag
       do kdiag = 1, idiag
         write (iout,*) ' diagonal offset = ', ioff(kdiag)
         write (iout,'(16f4.1)') (dns(k,kdiag),k=1,n)
       end do
c
c  reverse conversion
c
       job = 0
       call diacsr (n,job,idiag,dns,ndns,ioff,a,ja,ia)

        write(*,*) '-----------------------------------------'
         write(*,*) '  +++ matrix after conversion from diacsr  +++ '
        write(*,*) '-----------------------------------------'
                call dump(n,a,ja,ia,6)
c
c  checking the banded format
c
      lowd = 0
      job = 1
      call csrbnd(n,a,ja,ia,job,dns,ndns,lowd,ml,mu,ierr)
        iout = iout+1
        if (ierr .ne. 0) write (iout,*) ' ***** ERROR IN CSRBND'
        if (ierr .ne. 0) write (iout,*)   '     IERR = ', ierr
        write (iout,*) '-----------------------------------------'
         write (iout,*) '       +++  banded  format +++ '
         write (iout,*) ' bandwidth values found ml=',ml,'  mu=',mu
       write (iout,*) '-----------------------------------------'
       write (iout,'(4x,16i4)') (j,j=1,n)
       write (iout,'(3x,65(1h-))')
       do i=1, lowd
         write (iout,102) i, (dns(i,j), j=1,n)
       end do
c
c  convert back to a, ja, ia format.
c
      len = nnzmax
      call bndcsr(n,dns,ndns,lowd,ml,mu,a,ja,ia,len,ierr)
        write(*,*) ' IERR IN BNDCSR = ', ierr
        write(*,*) '-----------------------------------------'
        write(*,*) '  +++ matrix after conversion from bndcsr +++'
        write(*,*) '-----------------------------------------'
      call dump(n,a,ja,ia,6)
c
c  make sure it is sorted
c
        call csrcsc(n,1,1,a,ja,ia, a1,ja1,ia1)
c
c  checking skyline format.
c
         imod = 1
       call csrssk (n,imod,a1,ja1,ia1,a,ia,nnzmax,ierr)
       if (ierr .ne. 0) write (iout,*)   '     IERR = ', ierr
        write (iout,*) '-----------------------------------------'
         write (iout,*) '    +++ Sym. Skyline format +++ '
       write (iout,*) '-----------------------------------------'
       write (iout,'(3x,65(1h-))')
c
c  create column values.
c
       do i=1, n
         kend = ia(i+1)-1
         kstart = ia(i)
         do k=kstart,kend
           ja(k) =  i-(kend-k)
         end do
      end do

         call dump(n,a,ja,ia,iout)
c
c  back to ssr format.
c
      call sskssr (n,imod,a,ia,a1,ja1,ia1,nnzmax,ierr)
        write(*,*) '-----------------------------------------'
        write(*,*) '  +++ matrix after conversion from sskcsr +++'
       write(*,*) '-----------------------------------------'
      call dump(n,a1,ja1,ia1,6)
c
c  checking jad format.
c
c  first go back to the csr format ----
c
       call ssrcsr (n,a1,ja1,ia1,nnzmax,a,ja,ia,iwk,ierr)

       call csrjad (n, a, ja, ia, ndiag, iwk, a1, ja1, ia1)

       iout = iout+1
        write (iout,*) '-----------------------------------------'
        write (iout,*) '   +++   matrix in JAD format +++ '
        write (iout,*) '-----------------------------------------'
c
c  permutation array
c
                write (iout,*) ' ** PERMUTATION ARRAY '
         write (iout,'(17i4)') (iwk(k),k=1,n)
c
c  diagonals.
c
         do i=1,ndiag
           write (iout,*) ' J-diagonal number: ', i
           k1 = ia1(i)
           k2 = ia1(i+1)-1
           write (iout,104) (a1(k),k=k1,k2)
           write (iout,105) (ja1(k),k=k1,k2)
         end do
c
c  back to csr format.
c
        call jadcsr (n, ndiag, a1, ja1, ia1, iwk, a, ja, ia)

        write(*,*) '-----------------------------------------'
        write(*,*) '  +++ matrix after conversion from jadcsr +++'
       write(*,*) '-----------------------------------------'
      call dump(n,a,ja,ia,6)
c
c checking the linked list format
c
       nnz = ia(n+1) - ia(1)
       call csrlnk (n, a, ja, ia, iwk)
c
c  print links in file 7 (no need for another file)
c
       iout = 7
       write (iout,*) '-----------------------------------------'
       write (iout,*) '   +++   matrix in LNK format +++ '
       write (iout,*) '-----------------------------------------'
c
c  permutation array
c
                write (iout,*) ' LINK ARRAY '
           write (iout,*) ' ---------- '
         write (iout,'(17i4)') (iwk(k),k=1,nnz)
c
c  back to csr format..
c
        call lnkcsr (n, a, ja, ia, iwk, a1, ja1, ia1)

        write(*,*) '-----------------------------------------'
        write(*,*) '  +++ matrix after conversion from lnkcsr +++'
       write(*,*) '-----------------------------------------'
      call dump(n,a,ja,ia,6)

      CLOSE(UNIT=7)
      CLOSE(UNIT=8)
      CLOSE(UNIT=9)
      CLOSE(UNIT=10)
      CLOSE(UNIT=11)
      CLOSE(UNIT=12)
      CLOSE(UNIT=13)
      CLOSE(UNIT=14)

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'SPARSEKIT_PRB03'
      write ( *, '(a)' ) '  Normal end of execution.'

      write ( *, '(a)' ) ' '
      call timestamp ( )

      stop
      end
      function afun (x,y,z)

c*********************************************************************72
c
      implicit none

      double precision afun, x,y, z

      afun = -1.0

      return
      end
      function bfun (x,y,z)

c*********************************************************************72
c
      implicit none

      double precision bfun, x,y, z

      bfun = -1.0

      return
      end
      function cfun (x,y,z)

c*********************************************************************72
c
      implicit none

      double precision cfun, x,y, z

      cfun = -1.0

      return
      end
      function dfun (x,y,z)

c*********************************************************************72
c
      implicit none

      double precision dfun
      double precision gamma
      double precision x,y, z
      data gamma /100.0/
c        dfun = gamma*dexp(x*y)

      dfun = 10.0

      return
      end
      function efun (x,y,z)

c*********************************************************************72
c
      implicit none

      double precision efun
      double precision gamma
      double precision x,y, z

      data gamma /100.0/
c        efun = gamma*dexp(-x*y)

      efun = 0.0

      return
      end
      function ffun (x,y,z)

c*********************************************************************72
c
      implicit none

      double precision ffun, x,y, z

      ffun = 0.0

      return
      end
      function gfun (x,y,z)

c*********************************************************************72
c
      implicit none

      double precision gfun, x,y, z

      gfun = 0.0

      return
      end
      subroutine afunbl (nfree,x,y,z,coeff)

c*********************************************************************72
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

      double precision coeff(100)
      integer i
      integer j
      integer nfree
      double precision x
      double precision y
      double precision z

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

      double precision coeff(100)
      integer i
      integer j
      integer nfree
      double precision x
      double precision y
      double precision z

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

      double precision coeff(100)
      integer i
      integer j
      integer nfree
      double precision x
      double precision y
      double precision z

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

      double precision coeff(100)
      integer i
      integer j
      integer nfree
      double precision x
      double precision y
      double precision z

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

      double precision coeff(100)
      integer i
      integer j
      integer nfree
      double precision x
      double precision y
      double precision z

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

      double precision coeff(100)
      integer i
      integer j
      integer nfree
      double precision x
      double precision y
      double precision z

      do j=1, nfree
        do i=1, nfree
          coeff((j-1)*nfree+i) = 0.0
        end do
      end do

      return
      end
      subroutine xyk(nel,xyke,x,y,ijk,node)

c*********************************************************************72
c
c The material property function xyk for the finite element problem
c
      implicit none

      integer node

      integer ijk(node,*)
      integer nel
      double precision x(*)
      double precision xyke(2,2)
      double precision y(*)
c
c  this is the identity matrix.
c
      xyke(1,1) = 1.0
      xyke(2,2) = 1.0
      xyke(1,2) = 0.0
      xyke(2,1) = 0.0

      return
      end
