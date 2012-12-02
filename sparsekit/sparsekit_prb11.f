      program main

c*********************************************************************72
c
cc MAIN is a finite element matrix generator.
c
c  Discussion:
c
c    This driver will generate a finite element matrix for the
c    conduction problem
c
c      -Div ( K(x,y) Grad u ) = f
c      u = 0 on boundary
c
c    (Dirichlet boundary conditions). The matrix is returned
c    assembled in compressed sparse row format. Unassembled matrices
c    can be generated (using genfeu) but this is not supported yet.
c
c    This driver will provide a few grids if wanted, with an
c    arbitrary number of levels of refinement (as permitted by the
c    sizes of the arrays as declared below).
c
c  Modified:
c
c    02 July 2005
c
c  Reference:
c
c    Noborou Kikuchi
c    Finite element methods in mechanics,
c    Cambridge University Press, 1986.
c
      implicit none

      integer maxnx
      parameter ( maxnx = 2000 )
      integer maxnel
      parameter ( maxnel = 4000 )

      double precision a(7*maxnx)
      double precision f(3*maxnx)
      double precision fs(3*maxnx)
      character * ( 50 ) gridfile 
      character * ( 2 ) guesol
      integer ia(maxnx)
      integer ichild(8,maxnel)
      integer ierr
      integer ifmt
      integer ii
      integer iin
      parameter ( iin = 7 )
      integer ijk(3,maxnel)
      integer iout
      parameter ( iout = 8 )
      integer iparnts(2,maxnx)
      integer iu
      integer iwk(maxnx)
      integer ja(7*maxnx)
      integer job
      integer jwk(maxnx)
      character * ( 8 ) key
      character * ( 50 ) matfile
      integer n
      integer n2
      integer na
      parameter ( na = 3000 )
      integer nb
      integer ndeg 
       parameter ( ndeg = 8 )
      integer nelmax
      integer nelx
      integer nelxnew
      integer ngrid
      integer nodcode(maxnx)
      integer node
      parameter ( node = 3 )
      integer nref
      integer nx
      integer nxmax
      integer nxnew
      character * ( 72 ) title
      character * ( 3 ) type
      double precision x(maxnx)
      external xyk
      double precision y(maxnx)      
c
c  Choose starting grid.
c
      iu = 10

      call timestamp ( )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'SPARSEKIT_PRB11'
      write ( *, '(a)' ) '  FORTRAN77 version'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 
     &  '  Demonstrate the SPARSKIT routines that can generate'
      write ( *, '(a)' ) 
     &  '  test matrices based on finite element grids.'
      write ( *, '(a)' ) ' '
c
c  Force NGRID to be 1.  Normally, this would be chosen by the user.
c
      ngrid = 1
c
c  Generate the grid.
c
      if ( ngrid .eq. 7 ) then
        write(*,*)'Grid type 7 : user provided initial grid '
        write(*,*)'Filename for initial grid :'
        read(*,'(A)') gridfile
        open ( unit = iin, file = gridfile, STATUS = 'old' )
      end if

      nelx = 0
      nx = 0

      call ingrid(ngrid,iin,nx,nelx,node,x,y,nodcode,ijk)

      if ( ngrid .eq. 7 ) then
        close ( unit = iin )
      end if
c
c  Refine the grid.  Here we choose NREF=2, although this would
c  normally be interactively chosen.
c
      nref = 2
      nxmax = maxnx
      nelmax= maxnel
      nb = 0

      do ii = 1, nref
c
c estimate the number nx and nelx at next refinement level.
c
        call checkref(nx,nelx,ijk,node,nodcode,nb,nxnew,nelxnew)

        if ( nxmax .lt. nxnew .or. nelmax .lt. nelxnew ) then
          WRITE ( *, * ) 'Was able to do only ', ii-1 ,'  refinements'
          exit
        end if

        call refall(nx,nelx,ijk,node,ndeg,x,y,ichild,iparnts,nodcode, 
     &    nxmax, nelmax, ierr)

        if (ierr .ne. 0) then 
          WRITE ( *, * ) '** ERROR IN REFALL : ierr =',ierr
        end if

      end do

      job = 0

      call genfea ( nx, nelx, node, job, x, y, ijk, nodcode, fs, n2, 
     &  a, ja, ia, f, iwk, jwk, ierr, xyk )
c
c  Store matrix as a Harwell Boeing Matrix, by a call to prtmt.
c
      guesol = 'NN'
      title = ' FINITE ELEMENT TEST MATRIX FROM SPARSKIT            '
      type = 'RSA'
      key = 'SPARSKIT'
      ifmt = 104
      job = 2
      n = n2
c
c  Set the filename for the matrix data.
c
      matfile = 'test.mat'

      open ( unit = iout, file = matfile, STATUS = 'replace' )

      call prtmt(n,n,a,ja,ia,f,guesol,title,key,type, 
     &  ifmt,job,iout)

      close ( UNIT = IOUT )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 
     &  '  The information about the matrix generated has been'
      write ( *, '(a)' ) 
     &  '  stored in a file, using Harwell-Boeing format. '
      write ( *, '(a)' ) '  The file is "' // trim ( matfile ) // '".'

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'SPARSEKIT_PRB11'
      write ( *, '(a)' ) '  Normal end of execution.'

      write ( *, '(a)' ) ' '
      call timestamp ( )

      stop
      end
      subroutine ingrid ( ngrid, iin, nx, nelx, node, x, y, nodcode, 
     &  ijk )

c*********************************************************************72
c
cc INGRID initializes the grid according to the choice NGRID.
c
c  Discussion:
c
c    There are 6 initial grids provided and the user can
c   also enter his own grid as a seventh option.
c
c on entry:
c
c ngrid          = integer indicating the grid chosen. ngrid=1...6
c             corresponds to one of the 6 examples supplied by
c             SPARSKIT. ngrid = 7 is a user supplied initial grid.
c             see below for additional information for the format.
c iin       = integer containing the I/O unit number where to read
c             the data from in case ngrid = 7. A dummy integer
c             otherwise.
c node      = integer = the number of nodes per element (should be
c             set to three in this version). also the first dimension
c             of array ijk.
c
c on return
c 
c nx          = integer . the number of nodes
c nelx          = integer . the number of elements
c x, y      = two real arrays containing the coordinates of the nodes.
c nodcode   = an integer array containing the boundary information for
c             each node with the following meaning.
c      nodcode(i) = 0 -->  node i is internal
c      nodcode(i) = 1 -->  node i is a boundary but not a corner point
c      nodcode(i) = 2 -->  node i is a corner node.
c
c ijk(node,*)= an integer array containing the connectivity matrix.
c
c
c format for user supplied grid (when ngrid = 7)
c
c option 7 is a user defined initial grid.
c
c format is as follows:
c line 1: two integers, the first containing the number of nodes
c         the second the number of elements.
c line 2 to line nx+1:  node information
c        enter the following one line per node:
c        * the number of the node in the numbering chosen (integer from
c        taking the values 1 to nx), followed by
c        * the coordinates of the nodes (2 reals)  followed by
c       the boundary information, an integer taking one of the
c        values 0, 1, or 2,  with the meaning explained above.
c
c line nx+2 to nx+nelx+1: connectivity matrix
c       enter the following one line per element:
c       * number of the element in the numbering chosen, followed by
c       * The three numbers of the nodes (according to the above numbering
c       of the nodes) that constitute the element, in a counter clock-wise
c       order (this is in fact not important since it is checked by the
c       subroutine chkelemt).
c
c AN EXAMPLE: consisting on one single element (a triangle)
c------------
c    3    1
c    1    0.0000    0.0000    2
c    2    4.0000    0.0000    2
c    3    0.0000    4.0000    2
c    1    1    2    3
c
      implicit none

      integer node

      integer i
      integer ii
      integer iin
      integer ijk(node,*)
      integer j
      integer nelx
      integer ngrid
      integer nodcode(*)
      integer nx
      double precision x(*)
      double precision y(*)

      if ( ngrid .eq. 1 ) then

        call fgrid1 ( nx, nelx, node, x, y, nodcode, ijk )

      else if ( ngrid .eq. 2 ) then

        call fgrid2 ( nx, nelx, node, x, y, nodcode, ijk )

      else if ( ngrid .eq. 3 ) then

        call fgrid3 ( nx, nelx, node, x, y, nodcode, ijk )

      else if ( ngrid .eq. 4 ) then

        call fgrid4 ( nx, nelx, node, x, y, nodcode, ijk )

      else if ( ngrid .eq. 5 ) then

        call fgrid5 ( nx, nelx, node, x, y, nodcode, ijk )

      else if ( ngrid .eq. 6 ) then

        call fgrid6 ( nx, nelx, node, x, y, nodcode, ijk )

      else if ( ngrid .eq. 7 ) then

        read (iin,*) nx, nelx

        do i = 1, nx
          read(iin,*) ii, x(ii), y(ii), nodcode(ii)
        end do

        do i = 1, nelx

          read(iin,*) ii, (ijk(j,ii),j=1,node)
          nelx = max ( nelx, ii )

        end do

      end if

      call chkelmt ( nx, x, y, nelx, ijk, node )

      return
      end
      subroutine fgrid1 (nx,nelx,node,x,y,nodcode,ijk)

c*********************************************************************72
c
cc FGRID1: initial grid for a simple square with two elements.
c
c      3             4
c       --------------
c       |          . |
c       |   2    .   |
c       |      .     |
c       |   .    1   |
c       | .          |
c       --------------
c      1              2
c
c input parameters: node = first dimensoin of ijk (must be .ge. 3)
c output parameters:
c    nx    = number of nodes
c    nelx = number of elemnts
c    (x(1:nx), y(1:nx)) = coordinates of nodes
c    nodcode(1:nx) = integer code for each node with the
c          following meening:
c      nodcode(i) = 0 -->  node i is internal
c      nodcode(i) = 1 -->  node i is a boundary but not a corner point
c      nodcode(i) = 2 -->  node i is a corner point.
c   ijk(1:3,1:nelx) = connectivity matrix. for a given element
c          number nel, ijk(k,nel), k=1,2,3 represent the nodes
c          composing the element nel.
c
      implicit none

      integer node

      integer ijk(node,*)
      integer ijk1(2)
      integer ijk2(2)
      integer ijk3(2)
      integer k
      integer nelx
      integer nodcode(*)
      integer nx
      double precision x(*)
      double precision x1(4)
      double precision y(*)
      double precision y1(4)

      data ijk1 / 1, 1 /
      data ijk2 / 2, 4 /
      data ijk3 / 4, 3 /
      data x1 / 0.0, 1.0, 0.0, 1.0 /
      data y1 / 0.0, 0.0, 1.0, 1.0 /

      nx = 4

      do k = 1, nx
        x(k) = x1(k)
        y(k) = y1(k)
        nodcode(k) = 1
      end do

      nodcode(2) = 2
      nodcode(3) = 2

      nelx = 2

      do k = 1, nelx
        ijk(1,k) = ijk1(k)
        ijk(2,k) = ijk2(k)
        ijk(3,k) = ijk3(k)
      end do

      return
      end
      subroutine fgrid2 (nx,nelx,node,x,y,nodcode,ijk)

c*********************************************************************72
c
cc FGRID2: initial grid for a simple D-shaped region with 4 elemnts
c       6
c       | .
c       |    .
c       |      .
c       |   4     .
c       |           .
c     4 -------------- 5
c       |          . |
c       |   3    .   |
c       |      .     |
c       |   .    2   |
c       | .          |
c       --------------
c       | 2         . 3
c       |         .
c       |   1   .
c       |     .
c       |  .
c       |.
c       1
c
c input parameters: node = first dimensoin of ijk (must be .ge. 3)
c output parameters:
c    nx    = number of nodes
c    nelx = number of elemnts
c    (x(1:nx), y(1:nx)) = coordinates of nodes
c    nodcode(1:nx) = integer code for each node with the
c          following meening:
c      nodcode(i) = 0 -->  node i is internal
c      nodcode(i) = 1 -->  node i is a boundary but not a corner point
c      nodcode(i) = 2 -->  node i is a corner point.
c   ijk(1:3,1:nelx) = connectivity matrix. for a given element
c          number nel, ijk(k,nel), k=1,2,3 represent the nodes
c          composing the element nel.
c
      implicit none

      integer node

      integer ijk(node,*)
      integer ijk1(4)
      integer ijk2(4)
      integer ijk3(4)
      integer k
      integer nelx
      integer nodcode(*)
      integer nx
      double precision x(*)
      double precision x1(6)
      double precision y(*)
      double precision y1(6)
c
c coordinates of nodal points
c
      data x1/0.0, 0.0, 1.0, 0.0, 1.0, 0.0/
      data y1/0.0, 1.0, 1.0, 2.0, 2.0, 3.0/
c
c------------------|--|--|--|
c elements         1  2  3  4
c------------------|--|--|--|
      data ijk1   /1, 2, 2, 4/
      data ijk2   /3, 3, 5, 5/
      data ijk3   /2, 5, 4, 6/

      nx = 6

      do k = 1, nx
        x(k) = x1(k)
        y(k) = y1(k)
        nodcode(k) = 1
      end do

      nelx = 4

      do k = 1, nelx
        ijk(1,k) = ijk1(k)
        ijk(2,k) = ijk2(k)
        ijk(3,k) = ijk3(k)
      end do

      return
      end
      subroutine fgrid3 (nx,nelx,node,x,y,nodcode,ijk)

c*********************************************************************72
c
cc FGRID3: initial grid for a C-shaped region composed of 10 elements --
c
c
c      10           11            12
c       ---------------------------
c       |          . |          . |
c       |  7     .   |   9    .   |
c       |      .     |      .     |
c       |   .    8   |   .   10   |
c       | .          | .          |
c     7 ---------------------------
c       |          . |8           9
c       |   5    .   |
c       |      .     |
c       |   .    6   |
c     4 | .          |5           6
c       ---------------------------
c       |          . |          . |
c       |   1    .   |  3     .   |
c       |      .     |      .     |
c       |   .    2   |   .   4    |
c       | .          | .          |
c       ---------------------------
c      1             2            3
c
c
c input parameters: node = first dimensoin of ijk (must be .ge. 3)
c    nx    = number of nodes
c    nelx = number of elemnts
c    (x(1:nx), y(1:nx)) = coordinates of nodes
c    nodcode(1:nx) = integer code for each node with the
c          following meening:
c      nodcode(i) = 0 -->  node i is internal
c      nodcode(i) = 1 -->  node i is a boundary but not a corner point
c      nodcode(i) = 2 -->  node i is a corner point.
c   ijk(1:3,1:nelx) = connectivity matrix. for a given element
c          number nel, ijk(k,nel), k=1,2,3 represent the nodes
c          composing the element nel.
c
      implicit none

      integer node

      integer ijk(node,*)
      integer ijk1(10)
      integer ijk2(10)
      integer ijk3(10)
      integer k
      integer nelx
      integer nodcode(*)
      integer nx
      double precision x(*)
      double precision x1(12)
      double precision y(*)
      double precision y1(12)
c
c coordinates of nodal points
c
      data x1/0.0,1.0,2.0,0.0,1.0,2.0,0.0,1.0,2.0,0.0,1.0,2.0/
      data y1/0.0,0.0,0.0,1.0,1.0,1.0,2.0,2.0,2.0,3.0,3.0,3.0/
c
c------------------|--|--|--|--|--|--|---|---|---|
c elements         1  2  3  4  5  6  7   8   9  10
c------------------|--|--|--|--|--|--|---|---|---|
      data ijk1   /1, 1, 2, 2, 4, 4, 7,  7,  8, 8/
      data ijk2   /5, 2, 6, 3, 8, 5, 11, 8, 12, 9/
      data ijk3   /4, 5, 5, 6, 7, 8, 10, 11,11, 12/

      nx = 12

      do k = 1, nx
        x(k) = x1(k)
        y(k) = y1(k)
        nodcode(k) = 1
      end do

      nodcode(3) = 2
      nodcode(10) = 2
      nodcode(9) = 2

      nelx = 10

      do k = 1, nelx
        ijk(1,k) = ijk1(k)
        ijk(2,k) = ijk2(k)
        ijk(3,k) = ijk3(k)
      end do

      return
      end
      subroutine fgrid4 (nx,nelx,node,x,y,nodcode,ijk)

c*********************************************************************72
c
cc FGRID4: initial grid for a C-shaped region composed of 10 elements --
c      10                   11
c       +------------------+ .
c       | .                |    .
c       |    .       8     |       . 12
c       |        .         |  9   . |
c       |     7      .     |   .    |
c     7 |                . | .   10 |
c       -------------------+--------+ 9
c       |                 .| 8
c       |     5       .    |
c       |         .        |
c       |    .       6     |
c       |.                 | 5      6
c    4  +------------------+--------+
c       |               .  | .   4  |
c       |    1       .     |    .   |
c       |        .         |  3    .| 3
c       |    .        2    |    .
c       | .                | .
c       --------------------
c       1                  2
c
c input parameters: node = first dimensoin of ijk (must be .ge. 3)
c    nx    = number of nodes
c    nelx = number of elemnts
c    (x(1:nx), y(1:nx)) = coordinates of nodes
c    nodcode(1:nx) = integer code for each node with the
c          following meening:
c      nodcode(i) = 0 -->  node i is internal
c      nodcode(i) = 1 -->  node i is a boundary but not a corner point
c      nodcode(i) = 2 -->  node i is a corner point.
c   ijk(1:3,1:nelx) = connectivity matrix. for a given element
c          number nel, ijk(k,nel), k=1,2,3 represent the nodes
c          composing the element nel.
c
      implicit none

      integer node

      integer ijk(node,*)
      integer ijk1(10)
      integer ijk2(10)
      integer ijk3(10)
      integer k
      integer nelx
      integer nodcode(*)
      integer nx
      double precision x(*)
      double precision x1(12)
      double precision y(*)
      double precision y1(12)
c
c coordinates of nodal points
c
      data x1/0.0,1.0,1.5,0.0,1.0,1.5,0.0,1.0,1.5,0.0,1.0,1.5/
      data y1/0.0,0.0,0.5,1.0,1.0,1.0,2.0,2.0,2.0,3.0,3.0,2.5/
c
c------------------|--|--|--|--|--|--|---|---|---|
c elements         1  2  3  4  5  6  7   8   9  10
c------------------|--|--|--|--|--|--|---|---|---|
      data ijk1   /1, 1, 2, 5, 4, 4, 7, 10,  8, 8/
      data ijk2   /5, 2, 3, 3, 8, 5, 8,  8, 12, 9/
      data ijk3   /4, 5, 5, 6, 7, 8, 10, 11,11, 12/

      nx = 12

      do k = 1, nx
        x(k) = x1(k)
        y(k) = y1(k)
        nodcode(k) = 1
      end do

      nodcode(6) = 2
      nodcode(9) = 2

      nelx = 10

      do k = 1, nelx
        ijk(1,k) = ijk1(k)
        ijk(2,k) = ijk2(k)
        ijk(3,k) = ijk3(k)
      end do

      return
      end
      subroutine fgrid5 (nx,nelx,node,x,y,nodcode,ijk)

c*********************************************************************72
c
cc FGRID5: initial grid for a whrench shaped region composed of 14 elements --
c
c                                      13            15
c                                        . ----------.           |-3
c                                      .   .   13  .   .         |
c                                   .   12   .   .  14    .      |
c 9        10        11       12  .            . 14        . 16  |
c ----------------------------------------------------------     |-2
c |       . |       . |       . |            . |                 |
c | 1   .   |  3  .   |  5  .   |    7   .     |                 |
c |   .  2  |   .  4  |   .  6  |     .    8   |                 |
c |.        |.        |.        | .            |                 |
c -----------------------------------------------------------    |-1
c 1         2         3       4  .           6 .           . 8   |
c                                   .   9    .   .   11   .      |
c                                      .   .  10    .   .        |
c                                        .___________.           |-0
c                                       5             7
c
c 0---------1--------2----------3--------------4-------------5
c
c input parameters: node = first dimensoin of ijk (must be .ge. 3)
c    nx    = number of nodes
c    nelx = number of elemnts
c    (x(1:nx), y(1:nx)) = coordinates of nodes
c    nodcode(1:nx) = integer code for each node with the
c          following meening:
c      nodcode(i) = 0 -->  node i is internal
c      nodcode(i) = 1 -->  node i is a boundary but not a corner point
c      nodcode(i) = 2 -->  node i is a corner point.
c   ijk(1:3,1:nelx) = connectivity matrix. for a given element
c          number nel, ijk(k,nel), k=1,2,3 represent the nodes
c          composing the element nel.
c
      implicit none

      integer node

      integer ijk(node,*)
      integer ijk1(14)
      integer ijk2(14)
      integer ijk3(14)
      integer k
      integer nelx
      integer nodcode(*)
      integer nx
      double precision x(*)
      double precision x1(16)
      double precision y(*)
      double precision y1(16)
c
c coordinates of nodal points
c
      data x1/0.,1.,2.,3.,3.5,4.,4.5,5.,0.,1.,2.,3.,3.5,4.,4.5,5./
      data y1/1.,1.,1.,1.,0.,1.,0.,1.,2.,2.,2.,2.,3.,2.,3.,2./
c
c------------------|--|--|--|--|--|--|---|---|---|--|---|---|---|
c elements         1  2  3  4  5  6  7   8   9  10  11  12  13  14
c------------------|--|--|--|--|--|--|---|---|---|--|---|---|---|
      data ijk1   /1, 1, 2, 2, 3, 3, 4,  4,  4,  5, 6, 12, 14, 14/
      data ijk2   /10,2,11, 3,12, 4,14,  6,  5,  7, 7, 14, 15, 16/
      data ijk3   /9,10,10,11,11,12,12, 14,  6,  6, 8, 13, 13, 15/

      nx = 16

      do k=1, nx
        x(k) = x1(k)
        y(k) = y1(k)
        nodcode(k) = 1
      end do

      nodcode(9) = 2
      nodcode(8) = 2
      nodcode(16) = 2

      nelx = 14

      do k=1,nelx
        ijk(1,k) = ijk1(k)
        ijk(2,k) = ijk2(k)
        ijk(3,k) = ijk3(k)
      end do

      return
      end
      subroutine fgrid6 (nx,nelx,node,x,y,nodcode,ijk)

c*********************************************************************72
c
cc FGRID6 generates a random finite element grid. 
c
c  Discussion:
c
c    Random coordinates are generated by using the library random number
c    generator and then a Delauney triangulation is used to generate the grid.
c
c    The algorithm used for the triangulation is coded by Sweby.
c
      implicit none

      integer node

      integer adjlist(200,12)
      integer i
      integer i1
      integer i2
      integer ijk(node,*)
      integer ijktr(200,3)
      integer il(6)
      integer j
      integer jj
      integer k
      integer nadj(200)
      integer nbr
      integer nel(200)
      integer nelx
      integer nemax
      integer nod
      integer nodcode(*)
      integer nx
      double precision random
      double precision x(*)
      double precision y(*)

      nx = 20

      do j = 1, nx
        x(j) = random()
      end do

      do j = 1, nx
        y(j) = random()
      end do

      nemax = 200
      call dlauny ( x, y, nx, ijktr, nemax, nelx )

      print *, ' delauny -- nx, nelx ', nx, nelx

      do j = 1, nx
        nel(j) = 0
        nadj(j) = 0
      end do
c
c  transpose ijktr into ijk and count the number of
c  elemnts to which each node belongs.
c
      do j = 1, nelx
        do k = 1, node
          i = ijktr(j,k)
          ijk(k,j) = i
          nel(i) = nel(i)+1
        end do
      end do
c
c  Take care of ordering within each element.
c
      call chkelmt (nx, x, y, nelx, ijk, node)
c
c  The next blocks are to determine  the nature of each point.
c  (interior point, boundary point, corner point.
c
c  List and count the neighbors of each node.
c
      do j = 1, nelx

        do k=1, node
          il(k) = ijk(k,j)
          il(k+node) = il(k)
        end do
c
c  neighbors of node il(k) are il(k+1), il(k+2)
c
        do k = 1, node

          nod = il(k)
          i1 = il(k+1)
          i2 = il(k+2)
c
c  see if already there
c
          nbr = nadj(nod)

          do jj = 1, nbr

            if ( adjlist(nod,nbr) .eq. i1 ) then 
              i1 = 0
            end if

            if ( adjlist(nod,nbr) .eq. i2 ) then
              i2 = 0
            end if

          end do

          if ( i1 .ne. 0 ) then
            nbr = nbr + 1
            adjlist(nod,nbr) = i1
          end if

          if ( i2 .ne. 0 ) then
            nbr = nbr + 1
            adjlist(nod,nbr) = i2
          end if

          nadj(nod) = nbr
     
        end do

      end do
c
c  Boundary info:
c  if number of neighbors = number of elemnts to which it belongs then it is
c  an internal point
c  if not but number of neighnors >= 2 then boundary point
c  if nadj(k) = 2 then corner point.
c
      do j=1, nx

        nodcode(j) = 0
        nbr = nadj(j)

        if (nel(j) .lt. nbr) then
          nodcode(j) = 1
        end if

        if ( nbr .eq. 2 ) then
          nodcode(j) = 2
        end if

      end do

      return
      end
      function random ( )

c*********************************************************************72
c
cc RANDOM returns a pseudorandom value.
c
c  This routine was extracted from ELEFUNT.
c
      implicit none

      integer iy
      double precision random

      save iy

      data iy / 100001 /

      iy = iy * 125
      iy = iy - ( iy / 2796203 ) * 2796203
      random = dble ( iy ) / 2796203.0

      return
      end
      subroutine xyk ( nel, xyke, x, y, ijk, node )

c*********************************************************************72
c
cc XYK evaluates the material property matrix function.
c
c  Discussion:
c
c    In this example, the function is just the identity matrix.
c
      implicit none

      integer node

      integer ijk(node,*)
      integer nel
      double precision x(*)
      double precision xyke(2,2)
      double precision y(*)

      xyke(1,1) = 1.0
      xyke(2,2) = 1.0
      xyke(1,2) = 0.0
      xyke(2,1) = 0.0

      return
      end
