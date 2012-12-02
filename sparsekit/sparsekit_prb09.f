      program main

c*********************************************************************72
c
cc MAIN generates 5 point and 7-point matrices in Harwell-Boeing format. 
c
c  Creates a file containing a harwell-boeing matrix.
c
c  nz = 1 will create a 2-D problem
c
      implicit none

      integer nxmax
      parameter ( nxmax = 50 )
      integer nmx
      parameter ( nmx = nxmax * nxmax )

      double precision a(7*nmx)
      character * ( 2 ) guesol
      integer ia(nmx)
      integer iau(nmx)
      integer ifmt
      integer iout
      integer ja(7*nmx)
      integer job
      character * ( 8 ) key
      character * ( 50 ) matfil
      integer n
      integer nx
      integer ny
      integer nz
      double precision rhs(1)
      double precision stencil(7)
      character * ( 72 ) title
      character * ( 3 ) type

      call timestamp ( )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'SPARSEKIT_PRB09:'
      write ( *, '(a)' ) '  FORTRAN77 version'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 
     &  '  This program demonstrates the use of GEN57PT'
      write ( *, '(a)' ) 
     &  '  to generate a sparse matrix derived from a 5 or'
      write ( *, '(a)' ) 
     &  '  7 point stencil on an NX by NY by NZ grid.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 
     &  '  The matrix is then stored in Harwell-Boeing format'
      write ( *, '(a)' ) '  in a file, using routine PRTMT.'
c
c  Set data defining the matrix.
c
      nx = 10
      ny = 10
      nz = 1
c
c  Call GEN57PT to generate the matrix.
c
      call gen57pt(nx,ny,nz,a,ja,ia,iau,stencil)
c
c  Set parameters required for the Harwell-Boeing format.
c
      n = nx * ny * nz
      guesol = 'NN'
      title = ' 5-POINT TEST MATRIX FROM SPARSKIT                    '
      type = 'RUA'
      key = 'SC5POINT'
      ifmt = 104
      job = 2
c
c  Write matrix to file.
c
      iout = 7
      matfil = 'test.mat'

      open ( unit = IOUT, file = MATFIL, STATUS = 'replace' )

      call prtmt ( n, n, a, ja, ia, rhs, guesol,title,key,type,
     &  ifmt,job,iout)

      close ( unit = iout )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'SPARSEKIT_PRB09'
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

      double precision dfun, gamma, x,y, z
      data gamma /100.0/

      dfun = 10.0

      return
      end
      function efun (x,y,z)

c*********************************************************************72
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
