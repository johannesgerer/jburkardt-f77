      subroutine addcom ( isub, jsub, value, invprm, diag, xlnz,
     &  ixlnz, nzsub, xnzsub, neqn )

c*********************************************************************72
c
cc ADDCOM adds values to a matrix stored in compressed storage scheme.
c
c  Discussion:
c
c    ADDCOM is called once the equations and variables have
c    been reordered, to add numbers to the matrix A, which is
c    assumed to be stored according to the compressed storage
c    scheme used by the ND and QMD methods.  
c
c    The arrays DIAG and XLNZ, which are used as storage for array entries,
c    should be zeroed out before this routine is first called for a
c    particular problem.
c
c  Modified:
c
c    01 January 2009
c
c  Author:
c
c    Alan George, Joseph Liu
c
c  Reference:
c
c    Alan George, Joseph Liu,
c    Computer Solution of Large Sparse Positive Definite Systems,
c    Prentice Hall, 1981,
c    ISBN: 0131652745,
c    LC: QA188.G46.
c
c  Parameters:
c
c    Input, integer ISUB, JSUB, the row and column of the matrix to
c    which the value is to be added.  Note that the matrix is assumed
c    to be symmetric, and only the lower triangle is stored.  Hence,
c    if ISUB < JSUB, the routine ignores the request.
c
c    Input, double precision VALUE, the number to be added to A(ISUB,JSUB).
c
c    Input, integer INVPRM(NEQN), the inverse ordering, which should have
c    been created by calling PERM_INVERSE once the reordering is set.
c
c    Input/output, double precision DIAG(NEQN), the diagonal elements of the
c    matrix.
c
c    Input/output, double precision XLNZ(*), the nonzero subdiagonal elements
c    of the matrix.
c
c    Input, integer IXLNZ(NEQN+1), NZSUB(*), XNZSUB(NEQN), data structures which
c    define the compressed storage scheme.
c
c    Input, integer NEQN, the order of the matrix.
c
      integer neqn

      double precision diag(neqn)
      integer i
      integer invprm(neqn)
      integer isub
      integer ixlnz(neqn+1)
      integer j
      integer jsub
      integer k
      integer kstop
      integer kstrt
      integer ksub
      integer nzsub(*)
      double precision value
      double precision xlnz(*)
      integer xnzsub(neqn)
c
c  Figure out the current locations of the given row and column.
c
      i = invprm(isub)
      j = invprm(jsub)
c
c  if the entry is on the diagonal, update DIAG and return.
c
      if ( i .eq. j ) then
        diag(i) = diag(i) + value
        return
      endif
c
c  if the entry is above the diagonal, don't store it at all.
c
      if ( i .lt. j ) then
        return
       end if

      kstrt = ixlnz(j)
      kstop = ixlnz(j+1) - 1
 
      if ( kstop .lt. kstrt ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'ADDCOM - Fatal error!'
        write ( *, '(a)' ) '  The IXLNZ array is incorrect.'
        write ( *, '(a)' ) ' '
        write ( *, '(a,i8)' ) '  ISUB = ', isub
        write ( *, '(a,i8)' ) '  JSUB = ', jsub
        write ( *, '(a)' ) ' '
        write ( *, '(a,i8)' ) '  I = ', i
        write ( *, '(a,i8)' ) '  J = ', j
        write ( *, '(a)' ) ' '
        write ( *, '(a,i8)' ) '  IXLNZ(J) =     ', ixlnz(j)
        write ( *, '(a,i8)' ) '  IXLNZ(J+1)-1 = ', ixlnz(j+1) - 1
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) '  However, we must have that'
        write ( *, '(a)' ) '  IXLNZ(J) <= IXLNZ(J+1)-1.'
        stop
      end if
 
      ksub = xnzsub(j)
 
      do k = kstrt, kstop
 
        if ( nzsub(ksub) .eq. i ) then
          xlnz(k) = xlnz(k) + value
          return
        end if
 
        ksub = ksub + 1
 
      end do
 
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'ADDCOM - Fatal error!'
      write ( *, '(a)' ) '  No storage for entry ISUB, JSUB'
      write ( *, '(a,i8)' ) '  ISUB = ', isub
      write ( *, '(a,i8)' ) '  JSUB = ', jsub

      stop
      end
      subroutine addrcm ( isub, jsub, value, invprm, diag, xenv, env,
     &  ierror, neqn )

c*********************************************************************72
c
cc ADDRCM adds values to a matrix stored in the RCM scheme.
c
c  Discussion:
c
c    ADDRCM can be called, once the equations and variables have
c    been reordered for the RCM method, in order to enter values
c    into the matrix A.  Note that values are added, so that
c    the matrix storage should be zeroed out before calling
c    ADDRCM.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    31 December 2008
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Alan George, Joseph Liu,
c    Computer Solution of Large Sparse Positive Definite Systems,
c    Prentice Hall, 1981,
c    ISBN: 0131652745,
c    LC: QA188.G46.
c
c  Parameters:
c
c    Input, integer ISUB, the original row to which the entry is to be added.
c
c    Input, integer JSUB, the original column to which the entry is to be added.
c
c    Input, double precision VALUE, the value to be added to entry A(ISUB,JSUB).
c
c    Input, integer INVPRM(NEQN), the reverse ordering produced by
c    subroutine INVERSE.
c
c    Input/output, double precision DIAG(NEQN), the diagonal elements of 
c    the matrix.
c
c    Input/output, integer XENV(NEQN+1), double precision ENV(*), the envelope 
c    structure of the matrix.
c
c    Output, integer IERROR, error flag.
c    0 for no errors detected.
c    5 for error (isub,jsub) lies outside envelope.
c
c    Input, integer NEQN, the number of equations.
c
      integer neqn

      double precision diag(neqn)
      double precision env(*)
      integer i
      integer ierror
      integer invprm(neqn)
      integer isub
      integer j
      integer jsub
      integer k
      double precision value
      integer xenv(neqn+1)

      ierror = 0
      i = invprm(isub)
      j = invprm(jsub)

      if ( i .lt. j ) then
        return
       end if

      if ( i .eq. j ) then
        diag(i) = diag(i) + value
        return
      endif
 
      k = xenv(i+1) - i + j
 
      if ( k .lt. xenv(i) ) then
        ierror = 5
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'ADDRCM - Fatal error!'
        write ( *, '(a)' ) '  Indices outside of envelope.'
        write ( *, '(a)' ) '  ISUB = ', isub, ' JSUB = ', jsub
        stop
      end if

      env(k) = env(k) + value
 
      return
      end
      subroutine addrhs ( invprm, isub, neqn, rhs, value )

c*********************************************************************72
c
cc ADDRHS adds a quantity to a specific entry of the right hand side.
c
c  Discussion:
c
c    ADDRHS is called after the variables and equations have been
c    reordered.  addrhs allows the user to add the number value
c    to the i-th entry of the right hand side, which has been
c    permuted to a new location.  please zero out rhs once before
c    calling this routine to build up a right hand side.
c
c  Modified:
c
c    01 January 2009
c
c  Parameters:
c
c    Input, integer INVPRM(NEQN), the inverse ordering vector which must 
c    have been set up before calling ADDRHS.  INVPRM is constructed
c    by subroutine INVERSE.
c
c  isub   - the original index, before reordering, of the location
c           in rhs to which the number value is to be added.
c
c  neqn   - the number of entries in rhs and invprm.
c
c  rhs    - the right hand side vector.
c
c  value  - the number to be added to rhs.
c
      integer neqn

      integer i
      integer invprm(neqn)
      integer isub
      double precision rhs(neqn)
      double precision value

      i = invprm(isub)
 
      if ( i .lt. 0 ) then
        return
      end if

      if ( neqn .lt. i ) then
        return
      end if

      rhs(i) = rhs(i) + value
 
      return
      end
      subroutine addrqt ( isub, jsub, value, invprm, 
     &  diag, xenv, env, xnonz, nonz, nzsubs, neqn )

c*********************************************************************72
c
cc ADDRQT adds values to a matrix stored in the implicit block storage scheme.
c
c  Discussion:
c
c    addrqt adds a number to the (i,j)-th position of a matrix
c    stored using the implicit block storage scheme.  since addrqt
c    only adds new values to those currently in storage, the
c    space used to store the matrix must be initialized to 0
c    before numerical values are supplied.  addrqt can be used with the
c    rqt and 1wd methods.
c
c  Modified:
c
c    01 January 2009
c
c  Reference:
c
c    Alan George, Joseph Liu,
c    Computer Solution of Large Sparse Positive Definite Systems,
c    Prentice Hall, 1981,
c    ISBN: 0131652745,
c    LC: QA188.G46.
c
c  input
c
c    Input, integer ISUB, the row index of the number to be added.
c
c    Input, integer JSUB, the column index of the number to be added.
c
c    diag   - an array containing the diagonal elements of the matrix.
c
c    value  - the number to be added.
c
c    invprm - invprm(i) is the new position of the variable whose
c    original number is i.  invp can be set up by
c    calling subroutine inverse after the reordering
c    routines are called.
c
c    xenv, env - a pair of arrays which contain the envelope
c    structure of the diagonal blocks.
c
c    xnonz, nonz, nzsubs - levels structure containing the off-
c    block diagonal parts of the rows of the lower
c    triangle of the original matrix.
c
c  output
c
c    also, diag, env or nonz has been modified by having the number
c    stored in value added in the (i,j)-th position.
c
      integer neqn

      double precision diag(neqn)
      double precision env(*)
      integer i
      integer invprm(neqn)
      integer isub
      integer j
      integer jsub
      integer k
      integer kstop
      integer kstrt
      double precision nonz(*)
      integer nzsubs(*)
      double precision value
      integer xenv(neqn+1)
      integer xnonz(neqn)

      i = invprm(isub)
      j = invprm(jsub)
 
      if ( i .eq. j ) then
        diag(i) = diag(i) + value
        return
      end if
 
      if ( i .lt. j ) then
        return
      end if
c
c  The value goes within the diagonal envelope.
c
      k = xenv(i+1) - i + j
 
      if ( k .ge. xenv(i) ) then
        env(k) = env(k) + value
        return
      endif
c
c  The value goes outside the diagonal blocks.
c
      kstrt = xnonz(i)
      kstop = xnonz(i+1)-1
 
      do k = kstrt, kstop
 
        if ( nzsubs(k) .eq. j ) then
          nonz(k) = nonz(k) + value
          return
        end if
 
      end do
 
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'ADDRQT - Fatal error!'
      write ( *, '(a)' ) '  lack of storage!'
      write ( *, * ) '  i=',i
      write ( *, * ) '  isub=',isub
      write ( *, * ) '  j=',j
      write ( *, * ) '  jsub=',jsub

      stop
      end
      subroutine adj_env_size ( n, adj_row, adj_num, adj, perm, 
     &  perm_inv, env_size )

c*********************************************************************72
c
cc ADJ_ENV_SIZE computes the envelope size for an adjacency structure.
c
c  Discussion:
c
c    The matrix is assumed to be symmetric.
c
c    The variables are assumed to have been permuted.
c
c  Modified:
c
c    24 February 2007
c
c  Author:
c
c    Alan George, Joseph Liu
c
c  Reference:
c
c    Alan George, Joseph Liu,
c    Computer Solution of Large Sparse Positive Definite Systems,
c    Prentice Hall, 1981,
c    ISBN: 0131652745,
c    LC: QA188.G46.
c
c  Parameters:
c
c    Input, integer N, the order of the matrix.
c
c    Input, integer ADJ_ROW(N+1).  Information about row I is stored
c    in entries ADJ_ROW(I) through ADJ_ROW(I+1)-1 of ADJ.
c
c    Input, integer ADJ_NUM, the number of adjacency entries.
c
c    Input, integer ADJ(ADJ_NUM), the adjacency structure. 
c    For each row, the column indices of the nonzero entries.
c
c    Input, integer PERM(N), integer PERM_INV(N), the permutation
c    and inverse permutation.
c
c    Output, integer ENV_SIZE, the number of cells in the envelope.
c    You could think of this number as the sum of the
c    bandwidths of each row.
c
      implicit none

      integer adj_num
      integer n

      integer add
      integer adj(adj_num)
      integer adj_row(n+1)
      integer col
      integer env_size
      integer i
      integer j
      integer perm(n)
      integer perm_inv(n)
      integer row

      env_size = 0
 
      do i = 1, n

        row = perm(i)

        add = 0

        do j = adj_row(row), adj_row(row+1) - 1
          col = perm_inv(adj(j))
          if ( col .lt. i ) then
            add = max ( add, i - col )
          end if
        end do

        env_size = env_size + add

      end do

      return
      end
      subroutine adj_print ( n, nadj, xadj, adj )

c*********************************************************************72
c
cc ADJ_PRINT prints the adjacency information stored in XADJ and ADJNCY.
c
c  Discussion:
c
c    The list has the form:
c
c    row   nonzeros
c
c    1       2   5   9
c    2       7   8   9   15   78   79   81  86  91  99
c          100 103
c    3      48  49  53
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    18 December 2002
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    input, integer N, the number of equations.
c
c    input, integer NADJ, the dimension of adjncy.
c
c    input, integer XADJ(N+1), organizes the entries of adjncy
c    into rows.  The entries for row i are in entries xadj(i)
c    through XADJ(I+1)-1.
c
c    input, integer ADJ(NADJ), the adjacency structure, which contains,
c    for each row, the column indices of the nonzero entries.
c
       implicit none

       integer nadj
       integer n

       integer adj(nadj)
       integer i
       integer j
       integer jhi
       integer jlo
       integer jmax
       integer jmin
       integer xadj(n+1)

       write ( *, '(a)' ) ' '
       write ( *, '(a)' ) 'ADJ_PRINT'
       write ( *, '(a)' ) '  show adjacency structure of sparse matrix.'
       write ( *, '(a,i6)' ) '  the matrix order is ', n
       write ( *, '(a,i6)' ) '  the number of entries is ', nadj
       write ( *, '(a)' ) ' '
       write ( *, '(a)' ) '  row       nonzeros '
       write ( *, '(a)' ) ' '

       do i = 1, n

         jmin = xadj(i)
         jmax = xadj(i+1) - 1
         do jlo = jmin, jmax, 10
           jhi = min ( jlo + 9, jmax )
           if ( jlo .eq. jmin ) then
             write ( *, '(i6,6x,10i6)' ) i, (adj(j),j=jlo,jhi)
           else
             write ( *, '(6x,6x,10i6)' ) ( adj(j),j=jlo,jhi)
           end if
         end do

       end do

       return
      end
      subroutine adj_show ( adj, iband, invprm, nadj, n, perm, xadj )

c*********************************************************************72
c
cc ADJ_SHOW displays a symbolic picture of a matrix.
c
c  Discussion:
c
c    The matrix is defined by the adjacency information in xadj and adj,
c    with a possible permutation through perm and invprm.  the routine
c    also computes the bandwidth and the size of the envelope.
c
c    if no permutation has been done, you must set invprm(i) = perm(i) = i
c    before calling this routine.  otherwise, you must call perm_inverse to
c    get the inverse permutation invprm before calling showmat.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    01 January 2009
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Alan George, Joseph Liu,
c    Computer Solution of Large Sparse Positive Definite Systems,
c    Prentice Hall, 1981,
c    ISBN: 0131652745,
c    LC: QA188.G46.
c
c  Parameters:
c
c    Input, integer ADJ(NADJ), the adjacency structure.
c    for each row, it contains the column indices of the nonzero entries.
c
c    Output, integer IBAND, the bandwidth of the matrix.
c
c    Input, integer PERM(N), INVPRM(N), the permutation and inverse permutation.
c
c    Input, integer NADJ, the number of adjacency entries in ADJ.
c
c    Input, integer N, the order of the matrix.
c
c    Input, integer XADJ(N+1).  information about row I is stored
c    in entries XADJ(I) through XADJ(I+1)-1 of ADJ.
c
      implicit none

      integer n_max
      parameter ( n_max = 100 )

      integer nadj
      integer n

      integer add
      integer adj(nadj)
      character band(n_max)
      integer col
      integer i
      integer iband
      integer invprm(n)
      integer j
      integer jhi
      integer jlo
      integer k
      integer nonz
      integer perm(n)
      integer xadj(n+1)

      iband = 0
      nonz = 0

      if ( n_max .lt. n ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'ADJ_SHOW - Fatal error!'
        write ( *, '(a)' ) '  N is too large!'
        write ( *, '(a,i6)' ) '  Maximum legal value is ', n_max
        write ( *, '(a,i6)' ) '  Your input value was ', n
        stop
      end if

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'ADJ_SHOW:'
      write ( *, '(a)' ) '  display nonzero structure of matrix.'
      write ( *, '(a)' ) ' '

      do i = 1, n

        do k = 1, n
          band(k) = ' '
        end do

        band(i) = 'x'

        do j = xadj(perm(i)), xadj(perm(i)+1)-1
          col = invprm(adj(j))
          if ( col .lt. i ) then
            nonz = nonz + 1
          end if
          iband = max ( iband, i-col )
          band(col) = 'x'
        end do

        write ( *, '(i6,1x,100a1)' ) i, ( band(j), j = 1, n )

      end do

      write ( *, '(a)' ) ' '
      write ( *, '(a,i6)' ) 'lower bandwidth = ', iband
      write ( *, '(a,i6,a)' ) 
     &  'lower envelope contains ', nonz, ' nonzeros.'

      return
      end
      subroutine bshufl ( xadj, adjncy, perm, nblks, xblk,
     &  bnum, mask, subg, xls )

c*********************************************************************72
c
cc BSHUFL renumbers the nodes of each block to reduce its envelope.  
c
c  Discussion:
c
c    This routine renumbers the nodes of each block
c    so as to reduce its envelope.
c    nodes in a block with no neighbors in previous
c    blocks are renumbered by subrcm before the others.
c
c  Modified:
c
c    01 January 2009
c
c  Reference:
c
c    Alan George, Joseph Liu,
c    Computer Solution of Large Sparse Positive Definite Systems,
c    Prentice Hall, 1981,
c    ISBN: 0131652745,
c    LC: QA188.G46.
c
c  input parameters 
c
c  (xadj, adjncy) - the graph adjacency structure.
c
c  (nblks, xblk ) - the tree partitioning.
c
c  updated parameter 
c  perm - the permutation vector. on return, it contains
c  the new permutation.
c
c  working vectors -
c  bnum - stores the block number of each variable.
c  mask - mask vector used to prescribe a subgraph.
c  subg - vector used to contain a subgraph.
c  xls - index vector to a level structure.
c
      integer adjncy(*)
      integer bnum(*), mask(*), perm(*)
      integer subg(*), xblk(*), xls(*)
      integer xadj(*), i, ip, istop, istrt, j
      integer jstrt, jstop, k, nabor, nblks, nbrblk
      integer node, nsubg

      if ( nblks .le. 0 ) then
        return
      end if
c
c  initialization.  find the block number for each
c  variable and initialize the vector mask.
c
      do k = 1, nblks
        istrt = xblk(k)
        istop = xblk(k+1) - 1
        do i = istrt, istop
          node = perm(i)
          bnum(node) = k
          mask(node) = 0
        end do
      end do
c
c  for each block, find those nodes with no neighbors
c  in previous blocks and accumulate them in subg.
c  they will be renumbered before others in the block.
c
      do k = 1, nblks

        istrt = xblk(k)
        istop = xblk(k+1) - 1
        nsubg = 0

        do i = istrt, istop

          node = perm(i)
          jstrt = xadj(node)
          jstop = xadj(node+1) - 1

          if ( jstrt .le. jstop ) then

            do j = jstrt, jstop
              nabor = adjncy(j)
              nbrblk = bnum(nabor)
              if (nbrblk .lt. k) go to 400
            end do
            nsubg = nsubg + 1
            subg(nsubg) = node
            ip = istrt + nsubg - 1
            perm(i) = perm(ip)

          end if

400       continue

        end do
c
c  call subrcm to renumber the subgraph stored
c  in (nsubg, subg).
c
        if ( 0 .lt. nsubg ) then
          call subrcm ( xadj, adjncy, mask, nsubg,
     &      subg, perm(istrt),xls)
        endif
        
      end do

      return
      end
      subroutine copysi ( n, a, b )

c*********************************************************************72
c
cc COPYSI copies the N integer elements from the vector A to B.
c
c  Modified:
c
c    01 January 2009
c
c  Parameters:
c
c    Input, integer N, the number of entries to copy.
c
c    Input, integer A(N), the vector to copy.
c
c    Output, integer B(N), the copied vector.
c
      integer n

      integer a(n)
      integer b(n)
      integer i

      do i = 1, n
        b(i) = a(i)
      end do

      return
      end
      subroutine degree ( root, xadj, adjncy, mask, deg, ccsize, ls )

c*********************************************************************72
c
cc DEGREE computes node degrees in a connected component, for the RCM method.
c
c  Discussion:
c
c    This routine computes the degrees of the nodes in the connected component 
c    specified by MASK and ROOT.
c
c    Nodes for which MASK is zero are ignored.
c
c  Modified:
c
c    01 January 2009
c
c  Reference:
c
c    Alan George, Joseph Liu,
c    Computer Solution of Large Sparse Positive Definite Systems,
c    Prentice Hall, 1981,
c    ISBN: 0131652745,
c    LC: QA188.G46.
c
c  Parameters:
c
c    Input, integer ROOT, the node that defines the component.
c
c    (xadj, adjncy) - adjacency structure pair.
c
c    Input, integer MASK(N), specifies a section subgraph.
c
c    Output, integer DEGREE(N), the degrees of the nodes in the component.
c
c    Output, integer CCSIZE, the size of the component.
c
c    Workspace, integer LS(N), used to store the nodes of the
c    component level by level.
c
      integer adjncy(*)
      integer ccsize
      integer deg(*)
      integer i
      integer ideg
      integer j
      integer jstop
      integer jstrt
      integer lbegin
      integer ls(*)
      integer lvlend
      integer lvsize
      integer mask(*)
      integer nbr
      integer node
      integer root
      integer xadj(*)
c
c  Initialization.
c
c  the array xadj is used as a temporary marker to
c  indicate which nodes have been considered so far
c
      ls(1) = root
      xadj(root) = -xadj(root)
      lvlend = 0
      ccsize = 1
c
c  lbegin is the pointer to the beginning of the curi
c  level, and lvlend points to the end of this level.
c
100   continue

      lbegin = lvlend + 1
      lvlend = ccsize
c
c  find the degrees of nodes in the current level,
c  and at the same time, generate the next level.
c
      do i = lbegin, lvlend
        node = ls(i)
        jstrt = -xadj(node)
        jstop = iabs(xadj(node+ 1) ) - 1
        ideg = 0
        if ( jstop .lt. jstrt ) go to 300
        do 200 j = jstrt, jstop
          nbr = adjncy(j)
          if ( mask(nbr) .eq. 0 ) go to 200
          ideg = ideg + 1
          if ( xadj(nbr) .lt. 0 ) go to 200
          xadj(nbr) = -xadj(nbr)
          ccsize = ccsize + 1
          ls(ccsize) = nbr
200     continue
300     deg(node) = ideg
      end do
c
c  Compute the current level width.
c  if it is nonzero , generate another level.
c
      lvsize = ccsize - lvlend
      if ( lvsize .gt. 0 ) go to 100
c
c  Reset xadj to its correct sign and return.
c
      do i = 1, ccsize
        node = ls(i)
        xadj(node) = - xadj(node)
      end do

      return
      end
      subroutine elslv ( neqns, xenv, env, diag, rhs )

c*********************************************************************72
c
cc ELSLV solves a lower triangular system stored in the envelope format.
c
c  Discussion:
c
c    This routine solves a lower triangular system l x - rhs.
c
c    The factor l is stored in the envelope format.
c
c  Modified:
c
c    01 January 2009
c
c  input parameters 
c
c    Input, integer NEQNS, the number of equations.
c
c  (xenv, env) - array pair for the envelope of l.
c  diag - array for the diagonal of l.
c
c  updated parameter 
c  rhs - on input, it contains the right hand vector.
c  on return, it contains the solution vector.
c
      double precision diag(*)
      double precision env(*)
      double precision rhs(*)
      double precision s
      integer xenv(*), i, iband, ifirst, k, kstop
      integer kstrt, l, last, neqns
c
c  find the position of the first nonzero in rhs and put it in ifirst.
c
      ifirst = 0
100   ifirst = ifirst + 1
      if ( rhs(ifirst) .ne. 0.0D+00 ) go to 200
      if ( ifirst .lt. neqns ) go to 100
      return
200   last = 0
c
c  last contains the position of the most recently
c  computed nonzero component of the solution.
c
      do i = ifirst, neqns

        iband = xenv(i+1) - xenv(i)
        if ( iband .ge. i ) iband = i - 1
        s = rhs(i)
        l = i - iband
        rhs(i) = 0.0D+00
c
c  row of the envelope is empty, or corresponding
c  components of the solution are all zeros.
c
        if ( iband .eq. 0 .or. last .lt. l ) go to 400
        kstrt = xenv(i+1) - iband
        kstop = xenv(i+1) - 1
        do k = kstrt, kstop
          s = s - env(k)*rhs(l)
          l = l + 1
        end do

400     continue

        if ( s .ne. 0.0D+00 ) then
          rhs(i) = s / diag(i)
          last = i
        end if

      end do

      return
      end
      subroutine esfct ( neqns, xenv, env, diag, flag )

c*********************************************************************72
c
cc ESFCT factors a positive definite envelope matrix into L * L'.
c
c  Discussion:
c
c    The algorithm used is the standard bordering method.
c
c  Modified:
c
c    31 December 2008
c
c  Author:
c
c    Alan George, Joseph Liu
c
c  Reference:
c
c    Alan George, Joseph Liu,
c    Computer Solution of Large Sparse Positive Definite Systems,
c    Prentice Hall, 1981,
c    ISBN: 0131652745,
c    LC: QA188.G46.
c
c  Parameters:
c
c    Input, integer NEQNS, the number of equations.
c
c    Input, integer XENV(N+1), the envelope index vector.
c
c    Input/output, double precision ENV(*), on input, the envelope of A, 
c    and on output the envelope of the factor L.
c
c    Input/output, double precision DIAG(NEQNS), on input the diagonal of A, 
c    and on output the diagonal of L.
c
c    Output, integer FLAG, error flag.
c    0, no error, the factorization was carried out.
c    1, the matrix is not positive definite.
c
      implicit none

      integer neqns

      double precision diag(neqns)
      double precision env(*)
      integer flag
      integer i
      integer iband
      integer ifirst
      integer ixenv
      integer j
      integer jstop
      double precision temp
      integer xenv(neqns+1)

      flag = 0

      if ( diag(1) .le. 0.0D+00 ) then
        write ( *, '(a)' ) ' '
        write ( *, * ) 'ESFCT: Fail for DIAG(1)'
        flag = 1
        return
      end if

      diag(1) = sqrt ( diag(1) )
c
c  Loop over rows 2:NEQNS of the matrix
c
      do i = 2, neqns

        ixenv = xenv(i)
        iband = xenv(i+1) - ixenv
c
c  Compute row I of the triangular factor.
c
        temp = diag(i)

        if ( iband .ne. 0 ) then

          ifirst = i - iband

          call elslv ( iband, xenv(ifirst), env, 
     &      diag(ifirst), env(ixenv) )

          jstop = xenv(i+1) - 1

          do j = ixenv, jstop
            temp = temp - env(j) * env(j)
          end do

        end if

        if ( temp .le. 0.0D+00 ) then
          write ( *, '(a)' ) ' '
          write ( *, * ) 'ESFCT: Fail for DIAG(', i ,')'
          flag = 1
          return
        end if

        diag(i) = sqrt ( temp )

      end do

      return
      end
      subroutine euslv ( neqns, xenv, env, diag, rhs )

c*********************************************************************72
c
cc EUSLV solves an upper triangular system u x - rhs. 
c
c  Discussion:
c
c    the factor u is stored in the envelope format.
c
c  Modified:
c
c    01 January 2009
c
c  input parameters -
c
c    Input, integer NEQNS, the number of equations.
c
c  (xenv, env) - array pair for the envelope of u.
c  diag - array for the diagonal of u.
c
c  updated parameter -
c
c  rhs - on input, it contains the right hand side.
c  on output, it contains the solution vector.
c
      double precision diag(*)
      double precision env(*), rhs(*), s
      integer xenv(*), i, iband, k, kstop, kstrt, l
      integer neqns

      i = neqns + 1

100   continue

      i = i - 1 

      if ( i .eq. 0 ) then
        return
      end if

      if ( rhs(i) .eq. 0.0D+00 ) then
        go to 100
      end if

      s = rhs(i) / diag(i)
      rhs(i) = s

      iband = xenv(i+1) - xenv(i)

      if ( iband .ge. i ) then
        iband = i - 1
      end if

      if ( iband .eq. 0 ) then
        go to 100
      end if

      kstrt = i - iband
      kstop = i - 1
      l = xenv(i+1) - iband

      do k = kstrt, kstop
        rhs(k) = rhs(k) - s * env(l)
        l = l + 1
      end do

      go to 100
      end
      subroutine fn1wd ( root, xadj, adjncy, mask, nsep, sep, nlvl, 
     &  xls, ls )

c*********************************************************************72
c
cc FN1WD finds one-way dissectors of
c  a connected component specified by mask and root.
c
c  Modified:
c
c    01 January 2009
c
c  input parameters
c
c  root - a node that defines (along with mask) the
c  component to be processed.
c
c  (xadj, adjncy) - the adjacency structure.
c
c  output parameters 
c  nsep - number of nodes in the one-way dissectors.
c  sep - vector containing the dissector nodes.
c
c  updated parameter 
c  mask - nodes in the dissector have their mask values
c  set to zero.
c
c  working parameters
c  (xls, ls) - level structure used by the routine fnroot.
c
      integer adjncy(*)
      integer ls(*)
      integer mask(*), sep(*), xls(*)
      integer xadj(*), i, j, k, kstop, kstrt, lp1beg, lp1end
      integer lvl, lvlbeg, lvlend, nbr, nlvl, node
      integer nsep, root
      double precision deltp1, fnlvl, width

      call fnroot ( root, xadj, adjncy, mask,
     &  nlvl, xls, ls )

      fnlvl = dble ( nlvl )
      nsep = xls(nlvl + 1) - 1
      width = dble ( nsep ) / fnlvl
      deltp1 = 1.0D+00 
     &  + sqrt ( ( 3.0D+00 * width + 13.0D+00 ) / 2.0D+00 )

      if ( nsep .ge. 50 .and. deltp1 .le. 0.5D+00 * fnlvl ) then
        go to 300
      end if
c
c  the component is too small, or the level structure
c  is very long and narrow. return the whole component.
c
      do i = 1, nsep
        node = ls(i)
        sep(i) = node
        mask(node) = 0
      end do

      return
c
c  find the parallel dissectors.
c
300   nsep = 0
      i = 0
400   i = i + 1
      lvl = int ( dble ( i ) * deltp1 + 0.5D+00 )

      if ( lvl .ge. nlvl ) then
        return
      end if

      lvlbeg = xls(lvl)
      lp1beg = xls(lvl + 1)
      lvlend = lp1beg - 1
      lp1end = xls(lvl + 2) - 1

      do j = lp1beg, lp1end
        node = ls(j)
        xadj(node) = - xadj(node)
      end do
c
c  nodes in level lvl are chosen to form dissector.
c  include only those with neighbors in lvl+1 level.
c  xadj is used temporarily to mark nodes in lvl+1
c
      do j = lvlbeg, lvlend
        node = ls(j)
        kstrt = xadj(node)
        kstop = iabs(xadj(node+1)) - 1
        do 600 k = kstrt, kstop
          nbr = adjncy(k)
          if ( xadj(nbr) .gt. 0 ) go to 600
          nsep = nsep + 1
          sep(nsep) = node
          mask(node) = 0
          go to 700
600     continue

700     continue

      end do

      do j = lp1beg, lp1end
        node = ls(j)
        xadj(node) = - xadj(node)
      end do

      go to 400
      end
      subroutine fnbenv ( xadj, adjncy, perm, invp, nblks, xblk,
     &  xenv, envsze, smask, marker, rchset )

c*********************************************************************72
c
cc FNBENV finds the exact envelope structure of the diagonal blocks of 
c  the cholesky factor of a permuted partitioned matrix.
c
c  Modified:
c
c    01 January 2009
c
c  input parameters 
c
c  (xadj, adjncy) - adjacency structure of the graph.
c  (perm, invp) - the permutation vector and its inverse.
c  (nblks, xblk) - the partitioning.
c
c  output parameters _
c  xenv - the envelope index vector.
c  envsze - the size of the envelope.
c
c  working parameters 
c  smask - marks nodes that have been considered.
c  marker - is used by routine reach.
c  rchset - is used by the subroutine reach.
c  stores both reachable and neighborhood sets.
c
      integer adjncy(*)
      integer invp(*)
      integer marker(*)
      integer perm(*)
      integer rchset(*)
      integer smask(*), xblk(*)
      integer xadj(*), xenv(*), blkbeg, blkend, i
      integer ifirst, inhd, k, envsze, nblks, neqns
      integer newnhd, nhdsze, node, rchsze
c
c  Initialization.
c
      neqns = xblk(nblks+1) - 1
      envsze = 1

      do i = 1, neqns
        smask(i) = 0
        marker(i) = 1
      end do
c
c  loop over the blocks 
c
      do k = 1, nblks

        nhdsze = 0
        blkbeg = xblk(k)
        blkend = xblk(k+1) - 1

        do i = blkbeg, blkend
          node = perm(i)
          marker(node) = 0
        end do
c
c  loop through the nodes in current block ...
c
        do i = blkbeg, blkend
          node = perm(i)
          call reach ( node, xadj, adjncy, smask,
     &      marker, rchsze, rchset(blkbeg),
     &      newnhd, rchset(nhdsze+1) )
          nhdsze = nhdsze + newnhd
          ifirst = marker(node)
          ifirst = invp(ifirst)
          xenv(i) = envsze
          envsze = envsze + i - ifirst
         end do
c
c  reset marker values of nodes in nbrhd set.
c
        do inhd = 1, nhdsze
          node = rchset(inhd)
          marker(node) = 0
        end do
c
c  reset marker and smask values of nodes in
c  the current block.
c
        do i = blkbeg, blkend
          node = perm(i)
          marker(node) = 0
          smask(node) = 1
        end do

      end do

      xenv(neqns+1) = envsze
      envsze = envsze - 1

      return
      end
      subroutine fndsep ( root, xadj, adjncy, mask, nsep, sep,
     &  xls, ls )

c*********************************************************************72
c
cc FNDSEP finds a small separator for a connected component.
c
c  Discussion:
c
c    The connected component is specified by MASK.
c
c  Modified:
c
c    24 July 2012
c
c  Reference:
c
c    Alan George, Joseph Liu,
c    Computer Solution of Large Sparse Positive Definite Systems,
c    Prentice Hall, 1981,
c    ISBN: 0131652745,
c    LC: QA188.G46.
c
c  Parameters:
c
c    Input, integer ROOT, the node that determines the masked component.
c
c    Input, integer XADJ(*), integer ADJNCY(*), the adjacency structure pair.
c
c    Input/output, integer MASK(*), nodes in the separator will have their mask
c    values set to zero.
c
c    Output, integer NSEP, the number of variables in the separator.
c
c    Output, integer SEP(*), the separator nodes.
c
c    Workspace, integer XLS(*), LS(*), the level structure pair for level 
c    structure found by FNROOT.
c
      implicit none

      integer adjncy(*)
      integer i
      integer j
      integer jstop
      integer jstrt
      integer ls(*)
      integer mask(*)
      integer midbeg
      integer midend
      integer midlvl
      integer mp1beg
      integer mp1end
      integer nbr
      integer nlvl
      integer node
      integer nsep
      integer root
      integer sep(*)
      integer xls(*)
      integer xadj(*)
c
c  Determine the level structure associated with ROOT.
c
      call rootls ( root, xadj, adjncy, mask, nlvl, xls, ls )
c
c  If the number of levels is less than 3, return the
c  whole component as the separator.
c
      if ( nlvl .lt. 3 ) then

        nsep = xls(nlvl+1) - 1
        do i = 1, nsep
          node = ls(i)
          sep(i) = node
          mask(node) = 0
        end do

        return

      end if
c
c  Find the middle level of the rooted level structure.
c
      midlvl = ( nlvl + 2 ) / 2
      midbeg = xls(midlvl)
      mp1beg = xls(midlvl + 1)
      midend = mp1beg - 1
      mp1end = xls(midlvl+2)
c
c  The separator is obtained by including only those
c  MIDDLE-level nodes with neighbors in the MIDDLE+1
c  level.  XADJ is used temporarily to mark those
c  nodes in the MIDDLE+1 level.
c
      do i = mp1beg, mp1end
        node = ls(i)
        xadj(node) = - xadj(node)
      end do

      nsep = 0
      
      do i = midbeg, midend
        node = ls(i)
        jstrt = xadj(node)
        jstop = iabs ( xadj(node+1) ) - 1

        do j = jstrt, jstop

          nbr = adjncy(j)

          if ( xadj(nbr) .le. 0 ) then
            nsep = nsep + 1
            sep(nsep) = node
            mask(node) = 0
            go to 10
          end if

        end do

10      continue

      end do
c
c  Reset XADJ to its correct sign.
c
      do i = mp1beg, mp1end
        node = ls(i)
        xadj(node) = - xadj(node)
      end do

      return
      end
      subroutine fnenv ( neqns, xadj, adjncy, perm, invp, xenv, envsze, 
     &  bandw )

c*********************************************************************72
c
cc FNENV finds the envelope structure of a permuted matrix.
c
c  Modified:
c
c    01 January 2009
c
c  Parameters:
c
c    Input, integer NEQNS, the number of equations.
c
c    (xadj, adjncy) - array pair containing the adjacency
c    structure of the graph of the matrix.
c
c    perm,invp - arrays containing permutation data about
c    the reordered matrix.
c
c  output parameters 
c
c    xenv - index vector for the level structure
c    to be used to store the lower (or upper)
c    envelope of the reordered matrix.
c
c    envsze - is equal to xenv(neqns+1) - 1.
c
c    bandw - bandwidth of the reordered matrix.
c
      integer adjncy(*)
      integer invp(*)
      integer perm(*)
      integer xadj(*), xenv(*), bandw, i, iband
      integer ifirst, iperm, j, jstop, jstrt, envsze
      integer nabor, neqns

      bandw = 0
      envsze = 1

      do i = 1, neqns

        xenv(i) = envsze
        iperm = perm(i)
        jstrt = xadj(iperm)
        jstop = xadj(iperm + 1) - 1
c
c  find the first nonzero in row i.
c
        ifirst = i

        do j = jstrt, jstop

          nabor = adjncy(j)
          nabor = invp(nabor)
          ifirst = min ( ifirst, nabor )

        end do

        iband = i - ifirst
        envsze = envsze + iband
        bandw = max ( bandw, iband )

      end do

      xenv(neqns+1) = envsze
      envsze = envsze - 1

      return
      end
      subroutine fnlvls ( root, xadj, adjncy, nodlvl,
     &  nlvl, xls, ls )

c*********************************************************************72
c
cc FNLVLS generates a rooted level structure for
c  a masked connected subgraph, rooted at a pseudo
c  peripheral node. the level numbers are recorded.
c
c  Modified:
c
c    01 January 2009
c
c  input parameters 
c
c  (xadj, adjncy) - the adjacency structure.
c
c  output parameters 
c  nlvl - number of levels in the level structure found
c  (xls, ls) - the level structure returned.
c
c  updated parameters 
c  root - on input, with the array nodlvl, specifies
c  the component whose pseudo-peripheral node is
c  to be found. on output, it contains that node
c  nodlvl - an inpht it cpfrlftes a sfnttan subgrapu
c  on return, it contains the node level numbers.
c
      integer adjncy(*)
      integer ls(*)
      integer nodlvl(*)
      integer xls(*)
      integer xadj(*), j, lbegin, lvl, lvlend, nlvl
      integer node, root

      call fnroot ( root, xadj, adjncy, nodlvl,
     &  nlvl, xls, ls )
     
      do lvl = 1, nlvl
        lbegin = xls(lvl)
        lvlend = xls(lvl + 1) - 1
        do j = lbegin, lvlend
          node = ls(j)
          nodlvl(node) = lvl
        end do
      end do

      return
      end
      subroutine fnofnz ( xadj, adjncy, perm, invp,
     &  nblks, xblk, xnonz, nzsubs, maxnz )

c*********************************************************************72
c
cc FNOFNZ finds the column subscripts of
c  the off-block-diagonal nonzeros in the lower triangle
c  of a partitioned matrix.
c
c  Modified:
c
c    01 January 2009
c
c  input parameters 
c
c  (xadj, adjncy) - adjacency structure pair for the graph.
c  (perm, invp) - the permutation vectors.
c  (nblks, xblk) - the block partitioning.
c
c  output parameters 
c  (xnonz, nzsubs) - the column subscripts of the nonzeros
c  of a to the left of the diagonal blocks are
c  stored row by row in continguous locations in the
c  array nzsubs. xnonz is the ~
c
c  index vector to it.
c
c  updated parameter 
c  maxnz - on input, it contains the size of the vector
c  nzsubs; and on output, the number of nonzeros
c
      integer adjncy(*)
      integer invp(*)
      integer nzsubs(*)
      integer perm(*), xblk(*)
      integer xadj(*), xnonz(*), blkbeg, blkend, i, j
      integer jperm, jxnonz, k, kstop, kstrt, maxnz
      integer nabor, nblks, nzcnt

      nzcnt = 1
c
c  loop over the blocks.
c
      do i = 1, nblks

        blkbeg = xblk(i)
        blkend = xblk(i+1) - 1
c
c  loop over the rows of the i-th block.
c
        do j = blkbeg, blkend

          xnonz(j) = nzcnt
          jperm = perm(j)
          kstrt = xadj(jperm)
          kstop = xadj(jperm+1) - 1
          if ( kstrt .gt. kstop ) go to 200
c
c  loop over the nonzeros of row j.
c
          do k = kstrt, kstop
            nabor = adjncy(k)
            nabor = invp(nabor)
c
c  check to see if it is to the left of the i-th diagonal block.
c
            if ( nabor .lt. blkbeg ) then
              if ( nzcnt .le. maxnz ) nzsubs(nzcnt) = nabor
              nzcnt = nzcnt + 1
            end if

          end do
c
c  sort the subscripts of row j
c
          jxnonz = xnonz(j)
          if ( nzcnt - 1 .le. maxnz ) then
            call sorts1 (nzcnt - jxnonz, nzsubs(jxnonz))
          endif

200       continue

        end do

      end do

      if ( 0 < nblks ) then
        xnonz(blkend+1) = nzcnt
      end if

      maxnz = nzcnt - 1

      return
      end
      subroutine fnroot ( root, xadj, adjncy, mask, nlvl, xls, ls)

c*********************************************************************72
c
cc FNROOT implements a modified version of the
c  scheme by gibbs, poole, and stockmeyer to find pseudo-
c  peripheral nodes it determines such a node for the
c  section subgraph specified by mask and root.
c
c  Modified:
c
c    01 January 2009
c
c  input parameters -
c
c  (xadj, adjncy) - adjacency structure pair for the graph.
c  mask - specifies a section subgraph. nodes for which
c  mask is zero are ignored by fnroot.
c
c  updated parameter 
c  root - on input, it (along with mask) defines the
c
c
c  component for which a pseudo-peripheral node is
c  to be found. on output, it is the node obtained
c
c  output parameters 
c  nlvl - is the number of levels in the level structure
c  rooted at the node root.
c  (xls,ls) - the level structure array pair containing
c  the level structure found.
c
      integer adjncy(*)
      integer ls(*)
      integer mask(*), xls(*)
      integer xadj(*), ccsize, j, jstrt, k, kstop, kstrt
      integer mindeg, nabor, ndeg, nlvl, node, nunlvl, root
c
c  determine the level structure rooted at root.
c
      call rootls ( root, xadj, adjncy, mask, nlvl, xls, ls)
      ccsize = xls(nlvl+1) - 1

      if ( nlvl .eq. 1 .or. nlvl .eq. ccsize ) then
        return
      end if
c
c  pick a node with minimum degree from the last level.
c
100   continue

      jstrt = xls(nlvl)
      mindeg = ccsize
      root = ls(jstrt)
      if ( ccsize .eq. jstrt ) go to 400

      do j = jstrt, ccsize

        node = ls(j)
        ndeg = 0
        kstrt = xadj(node)
        kstop = xadj(node+1) - 1

        do k = kstrt, kstop
          nabor = adjncy(k)
          if ( mask(nabor) .gt. 0 ) then
            ndeg = ndeg + 1
          end if
        end do

        if ( ndeg .lt. mindeg ) then
          root = node
          mindeg = ndeg
        end if

      end do
c
c  and generate its rooted level structure.
c
400   continue

      call rootls ( root, xadj, adjncy, mask, nunlvl, xls, ls)
      if (nunlvl .le. nlvl) return
      nlvl = nunlvl
      if (nlvl .lt. ccsize ) go to 100
      return
      end
      subroutine fnspan ( xadj, adjncy, nodlvl, nspan, set,
     &  level, nadjs, adjs, leaf )

c*********************************************************************72
c
cc FNSPAN finds the span of a given subset
c  in a given level subgraph in a level structure.
c  the adjacent set of the span in the lower level is
c  also determined. if the span has an unnumbered node
c  in the higher level, an unnumbered leaf node (i.e. one
c  with no neighbor in next level) will be returned.
c
c  Modified:
c
c    01 January 2009
c
c  input parameters 
c
c  (xadj, adjncy) - the adjacent structure.
c  level - level number of the current set.
c
c  updated parameters 
c  (nspan, set) - the input set. on return, it contains
c  the resulting span set.
c  nodlvl - the level number vector. nodes considered
c  will have their nodlvl changed to zero.
c
c  output parameters 
c  (nadjs, adjs) - the adjacent set of the span in the
c  lower level.
c  leaf - if the span has an unnumbered higher level nodi
c  leaf returns an unnumbered leaf node in the le'
c  structure, otherwise, leaf is zero.
c
      integer adjncy(*)
      integer adjs(*)
      integer nodlvl(*), set(*)
      integer xadj(*), i, j, jstop, jstrt, leaf, level
      integer lvl, lvlm1, nadjs, nbr, nbrlvl, node
      integer nspan, setptr
c
c  initialization ...
c
      leaf = 0
      nadjs = 0
      setptr = 0

100   continue

      setptr = setptr + 1
      if ( setptr .gt. nspan ) then
        return
      end if
c
c  for each node in the partially spanned set ...
c
      node = set(setptr)
      jstrt = xadj(node)
      jstop = xadj(node + 1) - 1
      if ( jstop .lt. jstrt ) go to 100
c
c  for each neighbor of node, test its nodlvl value
c
      do 500 j = jstrt, jstop
        nbr = adjncy(j)
        nbrlvl = nodlvl(nbr)
        if (nbrlvl .le. 0) go to 500
        if (nbrlvl - level) 200, 300, 600
c
c  nbr is in level-1, add it to adjs.
c
200     nadjs = nadjs + 1
        adjs(nadjs) = nbr
        go to 400
c
c  nbr is in level, add it to the span set.
c
300     nspan = nspan + 1
        set(nspan) = nbr
400     nodlvl(nbr) = 0
500   continue
      go to 100
c
c  nbr is in level+1. find an unnumbered leaf node by
c  tracing a path up the level structure. then
c  reset the nodlvl values of nodes in adjs.
c
600   leaf = nbr
      lvl = level + 1
700   jstrt = xadj(leaf)
      jstop = xadj(leaf+1) - 1
      do 800 j = jstrt, jstop
        nbr = adjncy(j)
        if ( nodlvl(nbr) .le. lvl ) go to 800
          leaf = nbr
          lvl = lvl + 1
          go to 700
800   continue

      do i = 1, nadjs
        node = adjs(i)
        nodlvl(node) = level - 1
      end do

      return
      end
      subroutine fntadj ( xadj, adjncy, perm, invp,
     &  nblks, xblk, father, bnum )

c*********************************************************************72
c
cc FNTADJ determines the quotient tree adjacency structure of a graph. 
c
c  Discussion:
c
c    The structure is represented by the father vector.
c
c  Modified:
c
c    01 January 2009
c
c  input parameters 
c
c  (xadj, adjncy) - adjacency structure pair for the
c  (perm, invp) - the permutation vectors.
c  (nblks, xblk) - the tree partitioning. 
c
c  output parameters 
c  father - the father vector of the quotient tree.
c
c  working parameters 
c  bnum - temporary vector to store the block number
c  of each variable.
c
      integer adjncy(*)
      integer bnum(*)
      integer father(*)
      integer invp(*)
      integer perm(*), xblk(*)
      integer xadj(*), i, istop, istrt, j, jstop, jstrt
      integer k, nabor, nblks, nbm1, nbrblk, node
c
c  initialize the block number vector.
c
      do k = 1, nblks
        istrt = xblk(k)
        istop = xblk(k+1) - 1
        do i = istrt, istop
          node = perm(i)
          bnum(node) = k
        end do
      end do
c
c  for each block
c
      father(nblks) = 0
      nbm1 = nblks - 1

      do k = 1, nbm1

        istrt = xblk(k)
        istop = xblk(k+1) - 1
c
c  find its father block in the tree structure.
c
        do i = istrt, istop

          node = perm(i)
          jstrt = xadj(node)
          jstop = xadj(node+1) -1

	      do j = jstrt, jstop
            nabor = adjncy(j)
            nbrblk = bnum(nabor)
            if ( nbrblk .gt. k ) go to 500
          end do

          father(k) = 0

        end do

        go to 600
500     father(k) = nbrblk
600     continue

      end do

      return
      end
      subroutine fntenv ( xadj, adjncy, perm, invp,
     &  nblks, xblk, xenv, envsze )

c*********************************************************************72
c
cc FNTENV determines the envelope index
c  vector for the envelope of the diagonal blocks of a
c  tree partitioned system.
c
c  Modified:
c
c    01 January 2009
c
c  input parameters -
c
c  (xadj, adjncy) - adjacency structure pair for the graph.
c  (perm, invp) - the permutation vectors.
c  (nblks, xblk) - the tree partitioning.
c
c  output parameters -
c  xenv - the envelope index vector.
c  envsze - the size of the envelope found.
c
      integer adjncy(*)
      integer invp(*)
      integer perm(*), xblk(*)
      integer xadj(*), xenv(*), blkbeg, blkend
      integer i, ifirst, j, jstop, jstrt, k, kfirst
      integer envsze, nblks, nbr, node

      envsze = 1
c
c  loop through each block in the partitioning .
c
      do 400 k = 1, nblks

        blkbeg = xblk(k)
        blkend = xblk(k+1) - 1
c
c  kfirst stores the first node in the k-th block
c  that has a neighbour in the previous blocks.
c
        kfirst = blkend

        do 300 i = blkbeg, blkend

          xenv(i) = envsze
          node = perm(i)
          jstrt = xadj(node)
          jstop = xadj(node+1) - 1
          if ( jstop .lt. jstrt ) go to 300
c
c  ifirst stores the first nonzero in the
c  i-th row within the k-th block.
c
          ifirst = i
          do j = jstrt, jstop
            nbr = adjncy(j)
            nbr = invp(nbr)
            if ( nbr .lt. blkbeg ) go to 100
            if ( nbr .lt. ifirst ) ifirst = nbr 
            go to 200
100         if ( kfirst .lt. ifirst ) ifirst = kfirst
            if (i .lt. kfirst ) kfirst = i
200         continue
          end do

          envsze = envsze + i - ifirst
300     continue
400   continue
      xenv(blkend+1) = envsze
      envsze = envsze - 1
      return
      end
      subroutine gen1wd ( neqns, xadj, adjncy, mask,
     &  nblks, xblk, perm, xls, ls )

c*********************************************************************72
c
cc GEN1WD finds a one-way dissection partitioning
c  for a general graph. fnlwd is used for each connected
c  component.
c
c  Modified:
c
c    01 January 2009
c
c  input parameters 
c
c    Input, integer NEQNS, the number of equations.
c
c  (xadj, adjncy) - the adjacency structure pair.
c
c  output parameters 
c  (nblks, xblk) - the partitioning found.
c  perm - the one-way dissection ordering.
c
c  working vectors 
c  mask - is used to mark variables that have
c  been numbered during the ordering process.
c  (xls, ls) - level structure used by rootls
c
      integer adjncy(*)
      integer ls(*)
      integer mask(*), perm(*)
      integer xblk(*), xls(*)
      integer xadj(*), ccsize, i, j, k, lnum
      integer nblks, neqns, nlvl, node, nsep
      integer num, root

      do i = 1, neqns
        mask(i) = 1
      end do

      nblks = 0
      num = 0

      do i = 1, neqns

        if ( mask(i) .eq. 0 ) go to 400
c
c  find a one-way dissector for each component.
c
        root = i
        call fn1wd ( root, xadj, adjncy, mask,
     &    nsep, perm(num+1), nlvl, xls, ls )
        num = num + nsep
        nblks = nblks + 1
        xblk(nblks) = neqns - num + 1
        ccsize = xls(nlvl+1) - 1
c
c  number the remaining nodes in the component.
c  each component in the remaining subgraph forms
c  a new block in the partitioning.
c
        do 300 j = 1, ccsize
          node = ls(j)
          if ( mask(node) .eq. 0 ) go to 300
          call rootls ( node, xadj, adjncy, mask,
     &      nlvl, xls, perm(num+1) )
          lnum = num + 1
          num = num + xls(nlvl+1) - 1
          nblks = nblks + 1
          xblk(nblks) = neqns - num + 1
          do k = lnum, num
            node = perm(k)
            mask(node) = 0
          end do
          if ( num .gt. neqns ) go to 500
300     continue

400     continue          

      end do
c
c  since dissectors found first should be ordered last,
c  routine revrse is called to adjust the ordering
c  vector, and the block index vector.
c
500   continue

      call i4vec_reverse ( neqns, perm )
      call i4vec_reverse ( nblks, xblk )
      xblk(nblks+1) = neqns + 1

      return
      end
      subroutine gennd ( neqns, xadj, adjncy, mask,
     &  perm, xls, ls )

c*********************************************************************72
c
cc GENND finds a nested dissection ordering for a general graph.
c
c  Modified:
c
c    01 January 2009
c
c  input parameters 
c
c    Input, integer NEQNS, the number of equations.
c
c  (xadj, adjncy) - adjacency structure pair.
c
c  output parameters 
c  perm - the nested dissection ordering.
c  
c  working parameters 
c  mask - is used to mask off variables that have
c  been numbered during the orderng process,
c  (xls, ls) - this level structure pair is used as
c  temporary storage by fnroot.
c
      integer adjncy(*)
      integer ls(*)
      integer mask(*), perm(*)
      integer xls(*)
      integer xadj(*), i, neqns, nsep, num, root

      do i = 1, neqns
        mask(i) = 1
      end do

      num = 0

      do 300 i = 1, neqns
c
c  for each masked component ...
c
200     if ( mask(i) .eq. 0 ) go to 300
        root = i
c
c  find a separator and number the nodes next,
c
        call fndsep ( root, xadj, adjncy, mask,
     &    nsep, perm(num+1), xls, ls )
        num = num + nsep
        if ( num .ge. neqns ) go to 400
        go to 200
300   continue
c
c  since separators found first should be ordered
c  last, routine revrse is called to adjust the
c  ordering vector,
c
400   continue

      call i4vec_reverse ( neqns, perm )

      return
      end      
      subroutine genqmd ( neqns, xadj, adjncy, perm, invp, deg,
     &  marker, rchset, nbrhd, qsize, qlink, nofsub )

c*********************************************************************72
c
cc GENQMD implements the minimum degre]
c  algorithm. it makes use of the implicit represei
c  ation of the elimination graphs by quotient grapi
c  and the notion of indistinguishable nodes.
c  caution - the adjacency vector adjncy will be
c  destroyed.
c
c  Modified:
c
c    01 January 2009
c
c  input parameters 
c
c    Input, integer NEQNS, the number of equations.
c
c  (xadj, adjncy) - the adjacency structure.
c
c  output parameters 
c  perm - the minimum degree ordering
c  invp - the inverse of perm.
c
c  jrking parameters 
c  deg - the degree vector. deg(i) is negative mean
c  node i has been numbered.
c  marker - a marker vector, where marker(i) is
c  negative means node i has been merged wit
c  another node and thus can be ignored.
c  rchset - vector used for the reachable set.
c  nbrhd - vector used for the neighborhood set.
c  qsize - vector used to store the size of
c  indistinguishable supernodes.
c  qlink - vector to store indistinguishable nodes,
c  i, qlink(i), qlink(qlink(i)) ... are the
c  members of the supernode represented by i.
c
      integer adjncy(*)
      integer invp(*)
      integer perm(*), deg(*), marker(*)
      integer rchset(*), nbrhd(*), qsize(*), qlink(*)
      integer xadj(*), inode, ip, irch, j, mindeg, ndeg
      integer neqns, nhdsze, node, nofsub, np, num, nump1
      integer nxnode, rchsze, search, thresh
c
c  initialize degree vector and other working variab
c
      mindeg = neqns
      nofsub = 0
      
      do node = 1, neqns
        perm(node) = node
        invp(node) = node
        marker(node) = 0
        qsize(node) = 1
        qlink(node) = 0
        ndeg = xadj(node+1) - xadj(node)
        deg(node) = ndeg
        if ( ndeg .lt. mindeg ) mindeg = ndeg
      end do

      num = 0
c
c  perform threshold search to get a node of min degree.
c  variable search points to where search should start.
c
200   continue

      search = 1
      thresh = mindeg
      mindeg = neqns

300   continue

      nump1 = num + 1
      if ( nump1 .gt. search ) search = nump1
      do 400 j = search, neqns
        node = perm(j)
        if ( marker(node) .lt. 0 ) go to 400
        ndeg = deg(node)
        if ( ndeg .le. thresh ) go to 500
        if ( ndeg .lt. mindeg ) mindeg = ndeg
400   continue

      go to 200
c
c  node has minimum degree. find its reachable sets by
c  calling qmdrch.
c
500   continue

      search = j
      nofsub = nofsub + deg(node)
      marker ( node ) = 1
      call qmdrch (node, xadj, adjncy, deg, marker,
     &  rchsze, rchset, nhdsze, nbrhd )
c
c  eliminate all nodes indistinguishable from node.
c  they are given by node, qlink(node), ....
c
      nxnode = node
600   continue
      num = num + 1
      np = invp(nxnode)
      ip = perm(num)
      perm(np) = ip
      invp(ip) = np
      perm(num) = nxnode
      invp(nxnode) = num
      deg(nxnode) = - 1
      nxnode = qlink(nxnode)
      if (nxnode .gt. 0) go to 600

      if ( rchsze .le. 0 ) go to 800
c
c  update the degrees of the nodes in the reachable
c  set and identify indistinguishable nodes.
c
      call qmdupd ( xadj, adjncy, rchsze, rchset, deg,
     &  qsize, qlink, marker,
     &  rchset(rchsze+1), nbrhd(nhdsze+1) )
c
c  reset marker value of nodes in reach set.
c  update threshold value for cyclic search.
c  also call qmdqt to form new quotient graph.
c
      marker(node) = 0

      do irch = 1, rchsze
        inode = rchset(irch)
        if ( marker(inode) .lt. 0 ) go to 700
        marker(inode) = 0
        ndeg = deg(inode)
        if ( ndeg .lt. mindeg ) mindeg = ndeg
        if ( ndeg .gt. thresh ) goto 700
        mindeg = thresh
        thresh = ndeg
        search = invp(inode)
700     continue
      end do

      if ( nhdsze .gt. 0 ) then
        call qmdqt ( node, xadj,
     &    adjncy, marker, rchsze, rchset, nbrhd )
      endif

800   if ( num .lt. neqns ) go to 300

      return
      end
      subroutine genrcm ( neqns, xadj, adjncy, perm, mask, xls )

c*********************************************************************72
c
cc GENRCM finds the reverse Cuthill-Mckee ordering for a general graph.
c
c  Discussion:
c
c    For each connected component in the graph, GENRCM obtains the order
c    by calling the subroutine RCM.
c
c  Modified:
c
c    31 December 2008
c
c  Parameters:
c
c    Input, integer NEQNS, the number of equations.
c
c    Input, (xadj, adjncy) - array pair containing the adjacency
c    structure of the graph of the matrix.
c
c    Output, integer PERM(NEQNS), the RCM ordering.
c
c    Workspace, integer MASK(NEQNS), used to mark variables that have been
c    numbered during the ordering process.  It is initialized to 1, and set 
c    to zero as each item is numbered.
c
c    Workspace, xls - the index vector for a level structure. this
c    level structure is stored in the currently
c    unused spaces in the permutation vector perm.
c
      implicit none

      integer neqns

      integer adjncy(*)
      integer ccsize
      integer i
      integer mask(neqns)
      integer nlvl
      integer num
      integer perm(neqns)
      integer root
      integer xls(*)
      integer xadj(*)

      do i = 1, neqns
        mask(i) = 1
      end do

      num = 1

      do i = 1, neqns
c
c  for each masked connected component ...
c
        if ( mask(i) .ne. 0 ) then

          root = i
c
c  first find a pseudo-peripheral node root.
c
c  note that the level structure found by
c  fnroot is stored starting at perm(num).
c
          call fnroot ( root, xadj, adjncy, mask,
     &      nlvl, xls, perm(num) )
c
c  RCM is called to order the component
c  using root as the starting node.
c
          call rcm ( root, xadj, adjncy, mask,
     &      perm(num), ccsize, xls )

          num = num + ccsize

          if ( neqns .lt. num ) then
            return
          end if

        end if

      end do

      return
      end
      subroutine genrqt ( neqns, xadj, adjncy, nblks,
     &  xblk, perm, xls, ls, nodlvl )

c*********************************************************************72
c
cc GENRQT determines a
c  partitioned ordering for a possibly disconnected
c  graph using the refined quotient tree algorithm.
c
c  Modified:
c
c    01 January 2009
c
c  input parameters
c
c    Input, integer NEQNS, the number of equations.
c
c  (xadj, adjncy) - the adjacency structure.
c
c	output parameters 
c	(nblks, xblk) - the quotient tree partitioning.
c	perm - the permutation vector.
c
c	working parameters 
c	(xls, ls) - this level structure pair is used by
c	fnroot to find a pseudo-peripheral node.
c	nodlvl - a temporary vector to store the level
c	number of each node in a level structure.
c
      integer adjncy(*)
      integer ls(*), nodlvl(*), perm(*)
      integer xblk(*), xls(*)
      integer xadj(*), i, ixls, leaf, nblks, neqns, nlvl
      integer root
c
c  initialization.
c
      do i = 1, neqns
        nodlvl(i) = 1
      end do

      nblks = 0
      xblk(1) = 1
c
c  for each connected component, find a rooted level
c  structure, and then call rqtree for its block order.
c
      do i = 1, neqns

        if ( 0 .lt. nodlvl(i) ) then
          root = i
          call fnlvls ( root, xadj, adjncy, nodlvl, nlvl, xls, ls )
          ixls = xls(nlvl)
          leaf = ls(ixls)
          call rqtree ( leaf, xadj, adjncy, perm,
     &      nblks, xblk, nodlvl, xls, ls )
         end if

      end do

      return
      end
      subroutine gsfct ( neqns, xlnz, lnz, xnzsub, nzsub, diag,
     &  link, first, temp, flag )

c*********************************************************************72
c
cc GSFCT performs the symmetric
c  factorization for a general sparse system, stor
c  the compressed subscript data format.
c
c  Modified:
c
c    01 January 2009
c
c  input parameters 
c
c    Input, integer NEQNS, the number of equations.
c
c  xlnz - index vector for lnz. xlnz(i) points to the
c  start of nonzeros in column i of factor L.
c
c  (xnzsub, nzsub) - the compressed subscript data
c  structure for factor l.
c
c  updated parameters 
c
c  Input/output, LNZ(*). on input, contains nonzeros of a, and on
c  return, the nonzeros of l.
c
c  Input/output, double precision DIAG(*).  the diagonal of l overwrites 
c  that of a.
c
c  Output, integer FLAG, the error flag. it is set to 1 if a zero
c  or negative square root occurs during the factorization.
c
c  parameters 
c  LINK, at step j. the list in
c  link(j), link(link(j)), ...........
c  consists of those columns that will modify
c  the column l(*,j).
c
c  FIRST - temporary vector to point to the first
c  nonzero in each column that will be used
c  next for modification.
c
c  TEMP - a temporary vector to accumulate modifica'
c
      double precision diag(*)
      double precision lnz(*), temp(*), diagj, ljk
      integer link(*), nzsub(*)
      integer first(*), xlnz(*), xnzsub(*)
      integer i, flag, ii, istop, istrt, isub, j
      integer k, kfirst, neqns, newk

      flag = 0
c
c  Initialize working vectors.
c
      do i = 1, neqns
        link(i) = 0
      end do

      do i = 1, neqns
        temp(i) = 0.0D+00
      end do
c
c  compute column l(*,j) for j - 1,...., neqns.
c
      do 600 j = 1, neqns
c
c  for each column l(*,k) that affects l(s,j).
c
        diagj = 0.0D+00
        newk = link(j)
200     k = newk
        if (k.eq.0) go to 400
        newk = link(k)
c
c  outer product modification of l($ j) by
c  l(e,k) starting at first(k) of l(;,k).
c
        kfirst = first(k)
        ljk = lnz(kfirst)
        diagj = diagj + ljk*ljk
        istrt = kfirst + 1
        istop = xlnz(k+1) - 1
        if ( istop .lt. istrt ) go to 200
c
c  before modification, update vectors first,
c  and link for future modification steps.
c
        first(k) = istrt
        i = xnzsub(k) + (kfirst-xlnz(k)) + 1
        isub = nzsub(i)
        link(k) = link(isub)
        link(isub) = k
c
c  the actual mod is saved in vector temp.
c
        do ii = istrt, istop
          isub = nzsub(i)
          temp(isub) = temp(isub) + lnz(ii)*ljk
          i = i + 1
        end do

        go to 200
c
c  apply the modifications accumulated in temp to column l(*,j).
c
400     continue

        diagj = diag(j) - diagj

        if ( diagj .le. 0.0D+00 ) then
          flag = 1
          return
        end if

        diagj = sqrt(diagj)
        diag(j) = diagj
        istrt = xlnz(j)
        istop = xlnz(j+1) - 1
        if ( istop .lt. istrt ) go to 600
        first(j) = istrt
        i = xnzsub(j)
        isub = nzsub(i)
        link(j) = link(isub)
        link(isub) = j

        do ii = istrt, istop
          isub = nzsub(i)
          lnz(ii) = ( lnz(ii) - temp(isub) ) / diagj
          temp(isub) = 0.0D+00
          i = i + 1
        end do

600   continue

      return
      end
      subroutine gsslv ( neqns, xlnz, lnz, xnzsub, nzsub,
     &  diag, rhs )

c*********************************************************************72
c
cc GSSLV solves a factored system in compressed subscript sparse format.
c
c  Modified:
c
c    01 January 2009
c
c  input parameters 
c
c    Input, integer NEQNS, the number of equations.
c
c  {xlnz, lnz) - structure of nonzeros in l.
c  (xnzsub, nzsub) - compressed subscript structure.
c  diag - diagonal components of l.
c
c  updated parameter. 
c
c  rhs - on input, it contains the rhs vector, and on
c  output, the solution vector.
c
      double precision diag(*)
      integer i
      integer ii
      integer istop
      integer istrt
      integer isub
      integer j
      integer jj
      double precision lnz(*)
      integer neqns
      integer nzsub(*)
      double precision rhs(*)
      double precision rhsj
      double precision s
      integer xlnz(*)
      integer xnzsub(*)
c
c  Forward substitution.
c
      do j = 1, neqns
      
        rhsj = rhs(j) / diag(j)
        rhs(j) = rhsj
        istrt = xlnz(j)
        istop = xlnz(j+1) - 1
        i = xnzsub(j)

        do ii = istrt, istop
          isub = nzsub(i)
          rhs(isub) = rhs(isub) - lnz(ii) * rhsj
          i = i + 1
        end do

      end do
c
c  Backward substitution.
c
      j = neqns

      do jj = 1, neqns

        s = rhs(j)
        istrt = xlnz(j)
        istop = xlnz(j+1) - 1

        i = xnzsub(j)

        do ii = istrt, istop
          isub = nzsub(i)
          s = s - lnz(ii) * rhs(isub)
          i = i + 1
        end do

        rhs(j) = s / diag(j)
        j = j - 1

      end do

      return
      end
      subroutine i4vec_reverse ( n, a )

c*********************************************************************72
c
cc I4VEC_REVERSE reverses the elements of an I4VEC.
c
c  Discussion:
c
c    An I4VEC is a vector of I4's.
c
c    In FORTRAN90, call I4VEC_REVERSE is equivalent to:
c
c      A(1:N) = A(N:1:-1)
c
c  Example:
c
c    Input:
c
c      N = 5,
c      A = ( 11, 12, 13, 14, 15 ).
c
c    Output:
c
c      A = ( 15, 14, 13, 12, 11 ).
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    25 June 2009
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the number of entries in the array.
c
c    Input/output, integer A(N), the array to be reversed.
c
      implicit none

      integer n

      integer a(n)
      integer i
      integer t

      do i = 1, n / 2
        t        = a(i)
        a(i)     = a(n+1-i)
        a(n+1-i) = t
      end do

      return
      end
      subroutine perm_inverse ( neqn, perm, invprm ) 

c*********************************************************************72
c
cc PERM_INVERSE produces the inverse permutation.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    01 January 2009
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer NEQN, the number of equations.
c
c    Input, integer PERM(NEQN), the reordering of the variables and equations.
c
c    Output, integer INVPRM(NEQN), the inverse ordering, with the property 
c    that INVPRM(PERM(K))=K.
c
      integer neqn

      integer i
      integer invprm(neqn)
      integer k
      integer perm(neqn)

      do i = 1, neqn
        k = perm(i)
        invprm(k) = i
      end do

      return
      end
      subroutine permrv ( neqn, rhs, perm )

c*********************************************************************72
c
cc PERMRV should be called once the linear system has been solved and the
c  solution returned in rhs.  permrv then undoes the permutation of rhs, 
c  restoring the original ordering.  to do this, it needs the perm 
c  vector which defined the reordering used by the solver.
c
c  Modified:
c
c    01 January 2009
c
c  Parameters:
c
c    Input, integer NEQN, the number of equations.
c 
c  rhs    input/output, double precision rhs(neqn).  
c         on input, rhs contains the solution of the permuted linear 
c         system.
c         on output, rhs contains the solution of the original linear 
c         system.
c 
c  perm   input, integer perm(neqn), the permutation information.    
c         perm(i)=k means that the k-th equation and variable in the 
c         original ordering became the i-th equation and variable in the 
c         reordering.
c
      integer neqn

      integer i
      integer iput
      integer istart
      integer perm(neqn)
      double precision pull
      double precision put
      double precision rhs(neqn)
c
c  mark perm with negative signs which will be removed
c  as each permuted element is restored to its rightful place
c
      do i = 1, neqn
        perm(i) = - perm(i)
      enddo
c
c  search for the next element of perm which is the first
c  element of a permutation cycle
c
      istart = 0

   20 continue

      istart = istart + 1
      if ( istart .gt. neqn ) then
        return
      end if
      if(perm(istart).gt.0)go to 20
      if(iabs(perm(istart)).ne.istart)go to 30
      perm(istart)=iabs(perm(istart))
      go to 20
c
c  begin a cycle
c
   30 continue

      perm(istart)=iabs(perm(istart))
      iput=istart
      pull=rhs(iput)

   40 continue

      iput=iabs(perm(iput))
      put=rhs(iput)
      rhs(iput)=pull
      pull=put
      if(perm(iput).gt.0)go to 20
      perm(iput)=iabs(perm(iput))
      go to 40

      end
      subroutine qmdmrg ( xadj, adjncy, deg, qsize, qlink,
     &  marker, dego, nhdsze, nbrhd, rchset, ovrlp)

c*********************************************************************72
c
cc QMDMRG merges indistinguishable nodes in
c  the minimum degree ordering algorithm.
c  it also computes the new degrees of these
c  new supernodes.
c
c  Modified:
c
c    01 January 2009
c
c  input parameters 
c
c  (xadj, adjncy) - the adjacency structure.
c  deg0 - the number of nodes in the given set.
c  (nhdsze, nbrhd) - the set of eliminated supernodes
c  adjacent to some nodes in the set.
c
c  updated parameters 
c  deg - the degree vector.
c  qsize - size of indistinguishable nodes.
c  qlink - linked list for indistinguishable nodes.
c  marker - the given set is given by those nodes with
c  marker value set to 1. those nodes with degree
c  updated will have marker value set to 2.
c
c  working parameters 
c  rchset - the reachable set.
c  ovrlp - temp vector to store the intersection of two
c  reachable sets.
c
      integer adjncy(*), deg(*), qsize(*), qlink(*)
      integer marker(*), rchset(*), nbrhd(*), ovrlp(*)
      integer xadj(*), dego, deg1, head, inhd, iov, irch
      integer j, jstrt, jstop, link, lnode, mark, mrgsze
      integer nabor, nhdsze, node, novrlp, rchsze, root
c
c  initialization ...
c
      if ( nhdsze .le. 0 ) return

      do inhd = 1, nhdsze
        root = nbrhd(inhd)
        marker(root) = 0
      end do
c
c  loop through each eliminated supernode in the
c  (nhdsze, nbrhd).
c
      do 1400 inhd = 1, nhdsze
        root = nbrhd(inhd)
        marker (root) = - 1
        rchsze = 0
        novrlp = 0
        deg1 = 0
200     jstrt = xadj(root)
        jstop = xadj(root+1) - 1
c
c  determine the reachable set and its inters
c  ion with the input reachable set.
c
        do 600 j = jstrt, jstop
          nabor = adjncy(j)
          root = - nabor
          if (nabor) 200, 700, 300
300       mark = marker(nabor)
          if ( mark ) 600, 400, 500
400       rchsze = rchsze + 1
          rchset(rchsze) = nabor
          deg1 = deg1 + qsize(nabor)
          marker(nabor) = 1
          go to 600
500       if ( mark .gt. 1 ) goto 600
          novrlp = novrlp + 1
          ovrlp(novrlp) = nabor
          marker(nabor) = 2
600     continue
c
c  from the overlapped set, determine the node
c  that can be merged together.
c
700     head = 0
        mrgsze = 0
        do 1100 iov = 1, novrlp
          node = ovrlp(iov)
          jstrt = xadj(node)
          jstop = xadj(node+1) - 1
          do j = jstrt, jstop
            nabor = adjncy(j)
            if ( marker(nabor) .eq. 0 ) then
              marker(node) = 1
              go to 1100
            end if
          end do
c
c  node belongs to the new merged supernode.
c  update the vectors qlink and qsize.
c
          mrgsze = mrgsze + qsize(node)
          marker(node) = - 1
          lnode = node
900       link = qlink(lnode)
          if ( link .le. 0 ) goto 1000
          lnode = link
          go to 900
1000      qlink(lnode) = head
          head = node
1100    continue

        if ( head .le. 0 ) goto 1200
        qsize(head) = mrgsze
        deg(head) = dego + deg1 - 1
        marker (head) = 2
c
c  reset marker values.
c
1200    root = nbrhd(inhd)
        marker(root) = 0
        
        do irch = 1, rchsze
          node = rchset(irch)
          marker(node) = 0
        end do

1400  continue

      return
      end
      subroutine qmdqt ( root, xadj, adjncy, marker,
     &  rchsze, rchset, nbrhd )

c*********************************************************************72
c
cc QMDQT performs the quotient
c  transformation after a node has been elimin}
c
c  Modified:
c
c    01 January 2009
c
c  input parameters 
c  root - the node just eliminated. it becomes
c  representative of the new supernode.
c  (xadj, adjncy) - the adjacency structure.
c  (rchsze, rchset) - the reachable set of root
c  old quotient graph.
c  nbrhd - the neighborhood set which will be b
c  with root to form the new supernode.
c  marker - the marker vector.
c
c  updated parameter 
c  adjncy - becomes the adjncy of the quotient
c
      integer adjncy(*)
      integer marker(*), rchset(*), nbrhd(*)
      integer xadj(*), inhd, irch, j, jstrt, jstop, link
      integer nabor, node, rchsze, root

      irch = 0
      inhd = 0
      node = root
100   jstrt = xadj(node)
      jstop = xadj(node+1) - 2
      if ( jstop .lt. jstrt ) go to 300
c
c  place reach nodes into the adjacent list 4
c
      do j = jstrt, jstop
        irch = irch + 1
        adjncy(j) = rchset(irch)
        if ( irch .ge. rchsze ) go to 400
      end do
c
c  link to other space provided by the neighborhood set.
c
300   link = adjncy(jstop+1)
      node = - link
      if ( link .lt. 0 ) go to 100
      inhd = inhd + 1
      node = nbrhd(inhd)
      adjncy(jstop+1) = - node
      go to 100
c
c  all reachable nodes have been saved. end the adj list.
c  add root to the nbr list of each node in the reach set.
c
400   adjncy(j+1) = 0
      do 600 irch = 1, rchsze
        node = rchset(irch)
        if ( marker(node) .lt. 0 ) go to 600
        jstrt = xadj(node)
        jstop = xadj(node+1) - 1
        do 500 j = jstrt, jstop
          nabor = adjncy(j)
          if ( marker(nabor) .ge. 0 ) go to 500
          adjncy(j) = root
          go to 600
500     continue
600   continue

      return
      end
      subroutine qmdrch ( root, xadj, adjncy, deg, marker,
     &  rchsze, rchset, nhdsze, nbrhd )

c*********************************************************************72
c
cc QMDRCH determines the reachable set
c  a node through a given subset. the adjacency stre
c  is assumed to be stored in a quotient graph forma1
c
c  Modified:
c
c    01 January 2009
c
c input parameters 
c
c  root - the given node not in the subset.
c  (xadj, adjncy) - the adjacency structure pair.
c  deg - the degree vector. deg( i ) lt 0 means the n(
c  belongs to the given subset.
c
c  output parameters 
c  (rchsze, rchset) - the reachable set.
c  (nhdsze, nbrhd) - the neighborhood set.
c
c  updated parameters 
c  marker - the marker vector for reach and nbrhd se1
c  gt 0 means the node is in reach set.
c  lt 0 means the node has been merged with
c  others in the quotient or it is in nbrhd sl
c
      integer adjncy(*)
      integer deg(*), marker(*)
      integer rchset(*), nbrhd(*)
      integer xadj(*), i, istrt, istop, j, jstrt, jstop
      integer nabor, nhdsze, node, rchsze, root
c
c  loop through the neighbors of root in the quotient graph.
c
      nhdsze = 0
      rchsze = 0
      istrt = xadj(root)
      istop = xadj(root+1) - 1
      if ( istop .lt. istrt ) return
      
      do 600 i = istrt, istop
        nabor = adjncy(i)
        if ( nabor .eq. 0 ) return
        if ( marker(nabor) .ne. 0 ) go to 600
        if ( deg(nabor) .lt. 0 ) go to 200
c
c  include nabor into the reachable set.
c
        rchsze = rchsze + 1
        rchset(rchsze) = nabor
        marker(nabor) = 1
        go to 600
c
c  nabor has been eliminated. find nodes
c  reachable from it.
c
200     marker(nabor) = -1
        nhdsze = nhdsze + 1
        nbrhd(nhdsze) = nabor
300     jstrt = xadj(nabor)
        jstop = xadj(nabor+1) - 1
        
        do 500 j = jstrt, jstop
          node = adjncy(j)
          nabor = - node
          if (node) 300, 600, 400
400       if ( marker(node) .ne. 0 ) go to 500
          rchsze = rchsze + 1
          rchset(rchsze) = node
          marker(node) = 1
500     continue
       
600   continue

      return
      end
      subroutine qmdupd ( xadj, adjncy, nlist, list, deg,
     &  qsize, qlink, marker, rchset, nbrhd)

c*********************************************************************72
c
cc QMDUPD performs degree update for a set
c  of nodes in the minimum degree algorithm.
c
c  Modified:
c
c    01 January 2009
c
c  input parameters 
c
c  (xadj, adjncy) - the adjacency structure.
c  (nlist, list) - the list of nodes whose degree has to
c  be updated.
c
c  updated parameters 
c  deg - the degree vector.
c  qsize - size of indistinguishable supernodes.
c  qlink - linked list for indistinguishable nodes.
c  marker - used to mark those nodes in reach/nbrhd sets.
c
c  working parameters 
c  rchset - the reachable set.
c  nbrhd - the neighborhood set.
c
      integer adjncy(*)
      integer list(*), deg(*), marker(*)
      integer rchset(*), nbrhd(*), qsize(*), qlink(*)
      integer xadj(*), deg0, deg1, il, inhd, inode, irch
      integer j, jstrt, jstop, mark, nabor, nhdsze, nlist
      integer node, rchsze 
c
c  find all eliminated supernodes that are adjm
c  to some nodes in the given list. put them in'
c  (nhdsze, nbrhd). deg0 contains the number of
c  nodes in the list.
c
      if ( nlist .le. 0 ) then
        return
      end if

      deg0 = 0
      nhdsze = 0
      do 200 il = 1, nlist
        node = list(il)
        deg0 = deg0 + qsize(node)
        jstrt = xadj(node)
        jstop = xadj(node+1) - 1
        do 100 j = jstrt, jstop
          nabor = adjncy(j)
          if ( marker(nabor) . ne . 0 . or .
     &      deg(nabor) .ge. 0 ) go to 100
          marker(nabor) = - 1
          nhdsze = nhdsze + 1
          nbrhd(nhdsze) = nabor
100     continue
200   continue
c
c  merge indistinguishable nodes in the list by
c  calling the subroutine qmdmrg.
c
      if ( nhdsze .gt. 0 ) then
        call qmdmrg ( xadj, adjncy, deg, qsize, qlink,
     &    marker, deg0, nhdsze, nbrhd, rchset,
     &    nbrhd(nhdsze+1) )
      endif
c
c  find the new degrees of the nodes that have not been
c  merged.
c
      do 600 il = 1, nlist

        node = list(il)
        mark = marker(node)

        if ( mark .gt. 1 .or. mark .lt. 0 ) go to 600

        marker(node) = 2

        call qmdrch ( node, xadj, adjncy, deg, marker,
     &    rchsze, rchset, nhdsze, nbrhd)

        deg1 = deg0

        do irch = 1, rchsze
          inode = rchset(irch)
          deg1 = deg1 + qsize(inode)
          marker(inode) = 0
        end do

        deg(node) = deg1 - 1

        do inhd = 1, nhdsze
          inode = nbrhd(inhd)
          marker(inode) = 0
        end do

600   continue

      return
      end
      subroutine rcm ( root, xadj, adjncy, mask, perm, ccsize, deg )

c*********************************************************************72
c
cc RCM renumbers a connected component by the reverse Cuthill McKee algorithm.
c
c  Discussion:
c
c    The connected component is specified by mask and root.
c    The numbering is done starting at the node ROOT, and proceeding
c    using the Reverse Cuthill-McKee algorithm.
c
c    An outline of the algorithm is as follows:
c
c    X(1) = ROOT.
c
c    for ( I = 1 to N-1)
c      Find all unlabeled neighbors of X(I), 
c      assign them the next available labels, in order of increasing degree.
c
c    When done, reverse the ordering.
c
c  Modified:
c
c    31 December 2008
c
c  Author:
c
c    Alan George, Joseph Liu
c
c  Reference:
c
c    Alan George, Joseph Liu,
c    Computer Solution of Large Sparse Positive Definite Systems,
c    Prentice Hall, 1981,
c    ISBN: 0131652745,
c    LC: QA188.G46.
c
c  input parameters 
c
c    Input, integer ROOT, the node that defines the connected component.
c    It is used as the starting point for the RCM ordering.
c
c    Input, integer XADJ(N+1).  Information about row I is stored
c    in entries XADJ(I) through XADJ(I+1)-1 of ADJNCY.
c
c    Input, integer ADJNCY(*), the adjacency structure. 
c    For each row, the column indices of the nonzero entries.
c
c    Input/output, integer MASK(N), a mask for the nodes.  Only those nodes 
c    with nonzero input mask values are considered by the routine.  The
c    nodes numbered by RCM will have their mask values set to zero.
c
c    Output, integer PERM(N), the RCM ordering.
c
c    Output, integer CCSZE, the size of the connected component
c    that has been numbered.
c
c    Workspace, integer DEG(N), a temporary vector used to hold the degree
c    of the nodes in the section graph specified by MASK and ROOT.
c
      integer adjncy(*)
      integer ccsize
      integer deg(*)
      integer fnbr
      integer i
      integer j
      integer jstop
      integer jstrt
      integer k
      integer l
      integer lbegin
      integer lnbr
      integer lperm
      integer lvlend
      integer mask(*)
      integer nbr
      integer node
      integer perm(*)
      integer root
      integer xadj(*)
c
c  Find the degrees of the nodes in the component specified by MASK and ROOT.
c
      call degree ( root, xadj, adjncy, mask, deg, ccsize, perm )

      mask(root) = 0

      if ( ccsize .le. 1 ) then
        return
      end if

      lvlend = 0
      lnbr = 1
c
c  lbegin and lvlend point to the beginning and
c  the end of the current level respectively.
c
100   continue

      lbegin = lvlend + 1
      lvlend = lnbr

      do i = lbegin, lvlend
c
c  For each node in current level.
c
        node = perm(i)
        jstrt = xadj(node)
        jstop = xadj(node+i) - 1
c
c  Find the unnumbered neighbors of node.
c
c  fnbr and lnbr point to the first and last
c  nundneumibnred neighbors respectively of the curr
c
        fnbr = lnbr + 1

        do j = jstrt, jstop
          nbr = adjncy(j)
          if ( mask(nbr) .ne. 0 ) then
            lnbr = lnbr + 1
            mask(nbr) = 0
            perm(lnbr) = nbr
          end if
        end do
c
c  sort the neighbors of node in increasing
c  order by degree. linear insertion is used
c
        if ( fnbr .lt. lnbr ) then

          k = fnbr

300       continue

          l = k
          k = k + 1
          nbr = perm(k)

400       if ( l .lt. fnbr ) go to 500
          lperm = perm(l)
          if(deg(lperm) .le. deg(nbr) ) go to 500
          perm(l+1)=lperm
          l=l-1
          go to 400

500       continue

          perm(l+1)=nbr
 
          if ( k .lt. lnbr ) then
            go to 300
          end if

        end if

      end do

      if ( lnbr .gt. lvlend ) then
        go to 100
      end if
c
c  we now have the cuthill mckee ordering.
c  reverse it.
c
      k = ccsize / 2
      l = ccsize

      do i=1, k
        lperm = perm(l)
        perm(l) = perm(i)
        perm(i) = lperm
        l = l - 1
      end do

      return
      end
      subroutine reach ( root, xadj, adjncy, smask, marker,
     &  rchsze, rchset, nhdsze, nbrhd )

c*********************************************************************72
c
cc REACH determines the reachable set of a node.
c
c  Discussion:
c
c    The routine determines the reachable set of a node y through a subset s
c    (i.e. reach(y,s) ) in a given subgraph. moreover
c    it returns the neighborhood set of y in s. i.e
c    nbrhd(y,s), the set of nodes in s that can be
c    reached from y through a subset of s.
c
c  Modified:
c
c    01 January 2009
c
c  input parameters -
c
c    Input, integer ROOT, the given node not in the subset s.
c
c  (xadj, adjncy) - the adjacency structure pair.
c  smask - the mask vector for the set s.
c  - 0, if the node is not in s.
c  > 0, if the node is in s.
c
c  output parameters -
c  (nhdsze, nbrhd) - the neighborhood set.
c  (rchsze, rchset) - the reachable set.
c
c  updated parameters -
c  marker - the marker vector used to define the subgraph
c  nodes in the subgraph have marker value 0.
c  on return, the reachable and neighborhood node
c  sets have their marker values reset to root.
c
      integer adjncy(*)
      integer marker(*)
      integer nbrhd(*)
      integer rchset(*)
      integer smask(*)
      integer xadj(*), i, istop, istrt, j, jstop, jstrt
      integer nabor, nbr, nhdbeg, nhdptr, nhdsze, node
      integer rchsze, root
c
c  initialization ...
c
      nhdsze = 0
      rchsze = 0

      if ( marker(root) .le. 0 ) then
        rchsze = 1
        rchset(1) = root
        marker(root) = root
      end if

      istrt = xadj(root)
      istop = xadj(root+1) - 1
c
c  loop through the neighbors of root ...
c
      do 600 i = istrt, istop

        nabor = adjncy(i)
        if ( marker(nabor) .ne. 0 ) go to 600
        if ( smask(nabor) .gt. 0 ) go to 200
c
c  nabor is not in s. include it in the reach set.
c
        rchsze = rchsze + 1
        rchset(rchsze) = nabor
        marker(nabor) = root
        go to 600
c
c  nabor is in subset s. and has not been considerei
c  include it into the nbrhd set and find the nodes
c  reachable from root through this nabor.
c
200     nhdsze = nhdsze + 1
        nbrhd(nhdsze) = nabor
        marker(nabor) = root
        nhdbeg = nhdsze
        nhdptr = nhdsze
300     node = nbrhd(nhdptr)
        jstrt = xadj(node)
        jstop = xadj(node+1) - 1
        do 500 j = jstrt, jstop
          nbr = adjncy(j)
          if ( marker(nbr) .ne. 0 ) go to 500
          if ( smask(nbr) .eq. 0 ) go to 400
          nhdsze = nhdsze + 1
          nbrhd(nhdsze) = nbr
          marker(nbr) = root
          go to 500
400       rchsze = rchsze + 1
          rchset(rchsze) = nbr
          marker(nbr) = root
500     continue
        nhdptr = nhdptr + 1
        if ( nhdptr .le. nhdsze ) go to 300

600   continue

      return
      end
      subroutine rootls ( root, xadj, adjncy, mask, nlvl, xls, ls )

c*********************************************************************72
c
cc ROOTLS generates the level structure rooted at a given node.
c
c  Discussion:
c
c    Only those nodes for which MASK is nonzero will be considered.
c
c  Modified:
c
c    01 January 2009
c
c  input parameters 
c
c  root, the node at which the level structure is to be rooted.
c
c  (xadj, adjncy) - adjacency structure pair for the
c  given graph.
c
c  mask - is used to specify a section subgraph. nodes
c  with mask(i)=0 are ignored.
c
c  output parameters -
c  nlvl - is the number of levels in the level structure.
c  (xls, ls) - array pair for the rooted level structure.
c
      integer adjncy(*)
      integer ls(*), mask(*), xls(*)
      integer xadj(*), i, j, jstop, jstrt, lbegin
      integer ccsize, lvlend, lvsize, nbr, nlvl
      integer node, root
c
c  Initialization.
c
      mask(root) = 0
      ls(1) = root
      nlvl = 0
      lvlend = 0
      ccsize = 1
c
c  LBEGIN is the pointer to the beginning of the current
c  level, and lvlend points to the end of this level.
c
200   lbegin = lvlend + 1
      lvlend = ccsize
      nlvl = nlvl + 1
      xls(nlvl) = lbegin
c
c  generate the next level by finding all the masked
c  neighbors of nodes in the current level.
c
      do 400 i = lbegin, lvlend
        node = ls(i)
        jstrt = xadj(node)
        jstop = xadj(node + 1) - 1
        if ( jstop .lt. jstrt ) go to 400
        do 300 j = jstrt, jstop
          nbr = adjncy(j)
          if (mask(nbr) .eq. 0) go to 300
          ccsize = ccsize + 1
          ls(ccsize) = nbr
          mask(nbr) = 0
300     continue
400   continue
c
c  compute the current level width.
c  if it is nonzero, generate the next level.
c
      lvsize = ccsize - lvlend
      if (lvsize .gt. 0 ) go to 200
c
c  reset mask to one for the nodes in the level structure.
c
      xls(nlvl+1) = lvlend + 1
      do i = 1, ccsize
        node = ls(i)
        mask(node) = 1
      end do

      return
      end
      subroutine rqtree ( leaf, xadj, adjncy, perm,
     &  nblks, xblk, nodlvl, adjs, stack )

c*********************************************************************72
c
cc RQTREE finds a quotient tree ordering
c  for the component specified by leaf and nodlvl.
c
c  Modified:
c
c    01 January 2009
c
c  input parameters 
c
c  (xadj, adjncy) - the adjacency structure.
c
c  leaf - the input node that defines the connected
c  component. it is also a leaf node in the
c  rooted level structure passed to rqtree.
c  i.e. it has no neighbor in the next level.
c
c  output parameters 
c  perm - the permutation vector containing the ordering.
c  (nblks, xblk) - the quotient tree partitioning.
c
c  updated parameters 
c  nodlvl - the node level number vector. nodes in the
c  component have their nodlvl set to zero as
c  as they are nwbered.
c
c  working parameters 
c  adjs - temporary vector to store the adjacent set
c  of nodes in a particular level.
c  stack - temporary vector used to maintain the stack
c  of node subsets. it is organised as 
c  ( subset nodes, subset size, subset level )
c
      integer adjncy(*)
      integer adjs(*)
      integer nodlvl(*)
      integer perm(*)
      integer stack(*), xblk(*)
      integer xadj(*), blksze, ip, j, jp, leaf, level
      integer nadjs, nblks, node, npop, nuleaf
      integer num, toplvl, topstk
c
c  initialize the stack vector and its pointers.
c
      stack(1) = 0
      stack(2) = 0
      topstk = 2
      toplvl = 0
      num = xblk(nblks+1) - 1
c
c  form a leaf block, that is, one with no neighbors
c  in its next higher level.
c
100   level = nodlvl(leaf)
      nodlvl(leaf) = 0
      perm(num+1) = leaf
      blksze = 1
      call fnspan ( xadj, adjncy, nodlvl, blksze, perm(num+1),
     &  level, nadjs, adjs, nuleaf )
      if ( nuleaf .le. 0 ) go to 300
      jp = num

      do j = 1, blksze
        jp = jp + 1
        node = perm(jp)
        nodlvl(node) = level
      end do

      leaf = nuleaf
      go to 100
c
c  a new block has been found ...
c
300   nblks = nblks + 1
      xblk(nblks) = num + 1
      num = num + blksze
c
c  find the next possible block by using the adjacent
c  set in the lower level and the top node subset (if
c  appropriate) in the stack.
c
      level = level - 1
      if ( level .le. 0 ) go to 500
      call copysi ( nadjs, adjs, perm(num+1) )
      blksze = nadjs
      if ( level .ne. toplvl ) go to 400
c
c  the level of the node subset at the top of the
c  stack is the same as that of the adjacent set.
c  pop the node subset from the stack.
c
      npop = stack(topstk-1)
      topstk = topstk - npop - 2
      ip = num + blksze + 1
      call copysi ( npop, stack(topstk+1), perm(ip) )
      blksze = blksze + npop
      toplvl = stack(topstk)

400   continue

      call fnspan ( xadj, adjncy, nodlvl, blksze,
     &  perm(num+1), level, nadjs, adjs, nuleaf)

      if ( nuleaf .le. 0 ) go to 300
c
c  push the current node set into the stack.
c
      call copysi ( blksze, perm(num+1), stack(topstk+1) )
      topstk = topstk + blksze + 2
      stack(topstk-1) = blksze
      stack(topstk) = level
      toplvl = level
      leaf = nuleaf
      go to 100
c
c  before exit 
c
500   xblk(nblks+1) = num + 1
      return
      end
      subroutine setadj ( adj, irow, jcol, maxadj, nadj, neqn, xadj )

c*********************************************************************72
c
cc SETADJ sets up adjacency information in XADJ and ADJNCY.
c
c  the first call to setadj for a given problem should be with
c  irow=jcol=0, which is a flag telling setadj to clear out
c  the old information, and set xadj(i)=1 for i=1,neqn+1.
c
c  thereafter, just call setadj with values of irow and jcol
c  for each entry a(irow,jcol) of the matrix that is nonzero.
c  you only have to give the subdiagonal or superdiagonal
c  entries (irow.lt.jcol) or (jcol.lt.irow).  setadj automatically
c  'reflects' the information for the other half of the matrix.
c  repeated calls with the same values of irow and jcol do not
c  hurt.  no extra storage will be allocated.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    01 January 2009
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input/output, integer ADJ(MAXADJ), stores adjacency information.
c
c  irow   - the row index of a nonzero entry.
c           irow should be between 1 and neqn, except for the
c           initialization option, when irow=0 should be used.
c
c  jcol   - the column index of a nonzero entry.
c           jcol should be between 1 and neqn, except for the
c           initialization option, when jcol=0 should be used.
c
c  maxadj - the dimension of adjncy, the maximum amount of
c           storage you have allocated.
c
c    input/output, integer nadj, the number of adjaceny entries in adjncy.
c
c  neqn   - the number of equations represented by the matrix.
c           (the number of rows and columns).
c
c  xadj   - an integer vector, indexed by a row or column of the
c           matrix, and itself used to index adjncy.
c
      integer maxadj
      integer neqn

      integer adj(maxadj)
      integer i
      integer icall
      integer irow
      integer j
      integer jcol
      integer k
      integer kback
      integer khi
      integer klo
      integer nadj
      integer xadj(neqn+1)

      save icall

      data icall /0/

      if(irow.eq.0.or.jcol.eq.0.or.icall.eq.0)then
 
        icall = 1
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'SETADJ:'
        write ( *, '(a)' ) '  Initialize XADJ adjacency pointer array.'

        nadj = 0
        do i = 1, neqn + 1
          xadj(i) = 1
        end do
 
        return
 
      end if
 
      if ( irow .eq. jcol ) then
        return
      end if

      if(xadj(neqn+1).gt.maxadj+1)then
        write(*,*)' '
        write(*,*)'setadj - fatal error!'
        write(*,*)'  all available storage in adjncy has been used.'
        write(*,*)'  no more information can be stored!'
        stop
      endif
 
      if(irow.gt.neqn.or.irow.le.0.or.jcol.gt.neqn.or.jcol.le.0)then
        write(*,*)' '
        write(*,*)'setadj - fatal error!'
        write(*,*)'  illegal indices irow=',irow,' jcol=',jcol
        stop
      endif
 
      i=irow
      j=jcol

   20 continue

      klo=xadj(i)
      khi=xadj(i+1)-1
 
      do k=klo,khi

        if(adj(k).eq.j)then

          if(i.eq.irow)then
            i=jcol
            j=irow
            go to 20
          endif

          return

        endif

      enddo
 
      do k=xadj(i+1),xadj(neqn+1)
        kback=xadj(neqn+1)+xadj(i+1)-k
        adj(kback+1)=adj(kback)
      enddo
 
      adj(xadj(i+1))=j
 
      do k=i+1,neqn+1
        xadj(k)=xadj(k)+1
      end do
 
      nadj = xadj(neqn+1) - 1

      if(i.eq.irow)then
        i=jcol
        j=irow
        go to 20
      end if
 
      return
      end
      subroutine shomat ( adjncy, iband, ienv, invprm, maxadj, neqn,
     &  perm, xadj )

c*********************************************************************72
c
cc SHOMAT can display a symbolic picture of the matrix defined by
c  the adjacency information in xadj and adjncy, with a possible
c  permutation through perm and invprm.  shomat also computes
c  the bandwidth and the size of the envelope.
c
c  if no permutation has been done, you must set invprm(i)=perm(i)=i
c  before calling shomat.  otherwise, you must call invrse to
c  get the inverse permutation invprm before calling showmat.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    01 January 2009
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c  adjncy - an integer vector containing adjacency information.
c
c  iband  - compute by shomat, the bandwidth of the matrix.
c
c  ienv   - computed by shomat, the number of cells in the envelope.
c           you could think of this number as the sum of the
c           bandwidths of each row.
c
c  invprm - the inverse permutation vector.
c
c  maxadj - the dimension of adjncy.
c
c  neqn   - the number of equations.
c
c  xadj   - integer vector of dimension neqn+1, used to index
c           the vector adjncy.
c
      integer maxadj
      integer neqn

      integer adjncy(maxadj)
      character*1 band(100)
      integer i
      integer iadd
      integer iband
      integer ienv
      integer invprm(neqn)
      integer itemp
      integer j
      integer jhi
      integer jlo
      integer k
      integer perm(neqn)
      integer xadj(neqn+1)

      iband=0
      ienv=0
 
      if(neqn.gt.100)then
        write(*,*)' '
        write(*,*)'shomat - fatal error!'
        write(*,*)'  neqn is too large!'
        write(*,*)'  maximum legal value is 100.'
        write(*,*)'  your input value was ',neqn
        stop
      end if
 
      write(*,*)' '
      write(*,*)'shomat - display nonzero structure of matrix.'
      write(*,*)' '

        do i = 1, neqn

          do k = 1, neqn
            band(k) = ' '
          end do

          band(i) = 'x'

          iadd = 0
          jlo = xadj(perm(i))
          jhi = xadj(perm(i)+1)-1

          do j = jlo, jhi
            itemp = invprm(adjncy(j))
            if(i-itemp.gt.iband) iband=i-itemp
            band(itemp) = 'x'
            if(i-itemp.gt.iadd) iadd=i-itemp
          end do

          write(*,'(1x,i6,1x,100a1)')i,(band(j),j=1,neqn)
          ienv = ienv + iadd

        end do

      write(*,*)' '
      write(*,*)'lower bandwidth=',iband
      write(*,*)'lower envelope contains ',ienv,' entries.'

      return
      end
      subroutine smbfct ( neqns, xadj, adjncy, perm, invp,
     &  xlnz, maxlnz, xnzsub, nzsub, maxsub,
     &  rchlnk, mrglnk, marker, flag )

c*********************************************************************72
c
cc SMBFCT performs symbolic factorization on a permuted linear system.
c
c  Discussion:
c
c    The routine also sets up the compressed data structure for the system.
c
c  Modified:
c
c    01 January 2009
c
c  input parameters -
c
c    Input, integer NEQNS, the number of equations.
c
c  (xadj, adjncy) - the adjacency structure.
c
c  (perm, invp) - the permutation vector and its inverse.
c
c  updated parameters 
c  maxsub - size of the subscript array nzsub. on return,
c  it contains the number of subscripts used
c
c  output parameters 
c  xlnz - index into the nonzero storage vector lnz.
c  (xnzsub, nzsub) - the compressed subscript vectors.
c  maxlnz - the number of nonzeros found.
c  flag - error flag. positive value indicates that.
c  nzsub array is too small.
c
c  working parameters 
c  mrglnk - a vector of size neqns. at the kth step,
c  mrglnk(k), mrglnk(mrglnk(k)) , .........
c  is a list containing all those columns l($,j)
c  with j less than k, such that its first off-
c  diagonal nonzero is l(k,j). thus, the
c  nonzero structure of column l(*,k) can be found
c  by merging that of such columns l(t,j) with
c  the structure of a(*,k).
c  k - a vector of size neqns. it is used
c  the structure of each column l(-,k). at the
c  end of the kth step,
c  rchlnk(k), rchlnk(rchlnk(k)), ........
c  is the list of positions of nonzeros in column k
c  of the factoi2 t
c
c  to accumulate
c
c  marker - an integer vector of length neqns. it is used
c  to test if mass symbolic elimination can be
c  performed. that is, it is used to check whether
c  the structure of the current column k being
c  processed is completely determined by the single
c  column mrglnk(k).
c
      integer adjncy(*), invp(*), mrglnk(*), nzsub(*)
      integer perm(*), rchlnk(*), marker(*)
      integer xadj(*), xlnz(*), xnzsub(*)
      integer flag, i, inz, j, jstop, jstrt, k, knz
      integer kxsub, mrgk, lmax, m, maxlnz, maxsub
      integer nabor, neqns, node, np1, nzbeg, nzend
      integer rchm, mrkflg

      flag = 0
c
c  initialization ...
c
      nzbeg = 1
      nzend = 0
      xlnz(1) = 1

      do k = 1, neqns
        mrglnk(k) = 0
        marker(k) = 0
      end do
c
c  for each column knz counts the number
c  of nonzeros in column k accumulated in rchlnk.
c
      np1 = neqns + 1

      do 1500 k = 1, neqns

        knz = 0
        mrgk = mrglnk(k)
        mrkflg = 0
        marker(k) = k
        if (mrgk .ne. 0 ) marker(k) = marker(mrgk)
        xnzsub(k) = nzend
        node = perm(k)
        jstrt = xadj(node)
        jstop = xadj(node+1) - 1
        if (jstrt.gt.jstop) go to 1500
c
c  use rchlnk to link through the structure of
c  a(*,k) below diagonal
c
        rchlnk(k) = np1

        do 300 j = jstrt, jstop
          nabor = adjncy(j)
          nabor = invp(nabor)
          if ( nabor .le. k ) go to 300
          rchm = k
200       continue
          m = rchm
          rchm = rchlnk(m)
          if ( rchm .le. nabor ) go to 200
          knz = knz+1
          rchlnk(m) = nabor
          rchlnk(nabor) = rchm
          if ( marker(nabor) .ne. marker(k) ) mrkflg = 1
300     continue
c
c  test for mass symbolic elimination ...
c
        lmax = 0
        if ( mrkflg .ne. 0 .or. mrgk .eq. 0 ) go to 350
        if ( mrglnk(mrgk) .ne. 0 ) go to 350
        xnzsub(k) = xnzsub(mrgk) + 1
        knz = xlnz(mrgk+1) - (xlnz(mrgk) + 1)
        go to 1400
c
c  link through each column i that affects l(*,k
c
350     i = k
400     i = mrglnk(i)
        if (i.eq.0) go to 800
        inz = xlnz(i+1) - (xlnz(i)+1)
        jstrt = xnzsub(i) + 1
        jstop = xnzsub(i) + inz
        if (inz.le.lmax) go to 500
        lmax = inz
        xnzsub(k) = jstrt
c
c  merge structure of l(*,i) in nzsub into rci
c
500     rchm = k

        do 700 j = jstrt, jstop
          nabor = nzsub(j)
600       m = rchm
          rchm = rchlnk(m)
          if (rchm.lt.nabor) go to 600
          if (rchm.eq.nabor) go to 700
          knz = knz+1
          rchlnk(m) = nabor
          rchlnk(nabor) = rchm
          rchm = nabor
700     continue

        go to 400
c
c  check if subscripts duplicate those of another column.
c
800     if (knz.eq.lmax) go to 1400
c
c  or if tail of k-1st column matches head of kth.
c
        if (nzbeg.gt.nzend) go to 1200
        i = rchlnk(k)
        
        do 900 jstrt=nzbeg,nzend
          if (nzsub(jstrt)-i) 900, 1000, 1200
900     continue

        go to 1200
1000    xnzsub(k) = jstrt

        do j = jstrt, nzend
          if (nzsub(j).ne.i) go to 1200
          i = rchlnk(i)
          if (i.gt.neqns) go to 1400
        end do

        nzend = jstrt - 1
c
c  copy the structure of l($,k) from rchlnk
c  to the data structure (xnzsub. nzsiir~
c
1200    nzbeg = nzend + 1
        nzend = nzend + knz
        if (nzend.gt.maxsub) go to 1600
        i = k
        do j=nzbeg,nzend
          i = rchlnk(i)
          nzsub(j) = i
          marker(i) = k
        end do

        xnzsub(k) = nzbeg
        marker(k) = k
c
c  update the vector mrglnk. note column l(e,k) just found
c  is required to determine column l(*,j), where
c  l(j,k) is the first nonzero in l(e,k) below diagonal.
c
1400    if (knz.le.1) go to 1500
        kxsub = xnzsub(k)
        i = nzsub(kxsub)
        mrglnk(k) = mrglnk(i)
        mrglnk(i) = k
1500    xlnz(k+1) = xlnz(k) + knz

      maxlnz = xlnz(neqns) - 1
      maxsub = xnzsub(neqns)
      xnzsub(neqns+1) = xnzsub(neqns)
      flag = 0
      return
c
c  error - insufficient storage for nonzero subscripts.
c
1600  flag=1
      return
      end
      subroutine sorts1 ( na, array )

c*********************************************************************72
c
cc SORTS1 uses linear insertion to sort integers into increasing order.
c
c  Modified:
c
c    01 January 2009
c
c  Parameters:
c
c    Input, integer NA, the array size.
c
c    Input/output, integer ARRAY(NA).  On output, the entries have been
c    sorted in increasing order.
c
      integer na

      integer array(na)
      integer k
      integer l
      integer node

      do k = 2, na
        node = array(k)
        l = k - 1
100     if(l.lt. 1) go to 200
        if ( array(l) .le. node) go to 200
        array(l+1) = array(l)
        l = l - 1
        go to 100
200	    array(l+1) = node
      end do

      return
      end
      subroutine subrcm ( xadj, adjncy, mask, nsubg, subg, perm, xls )

c*********************************************************************72
c
cc SUBRCM finds the rcm ordering for a given subgraph (possibly disconnected).
c
c  Modified:
c
c    01 January 2009
c
c  input parameters 
c
c  (xadj, adjncy) - adjacency structure pair for the graph.
c
c  (nsubg, subg) - the given subgraph. nsubg is the
c  the size of the subgraph, and subg contains
c  the nodes in it.
c
c  output parameter 
c
c  perm - the permutation vector. it is also used
c  temporarily to store a level structure.
c
c  working parameters 
c
c  mask - mask vector with all zeros. it is used to
c  specify nodes in the subgraph.
c
c  xls - index to a level structure. note that the level
c  structure is stored in part of perm
c
      integer adjncy(*)
      integer ccsize
      integer i
      integer mask(*)
      integer nlvl
      integer node
      integer nsubg
      integer num
      integer perm(*)
      integer subg(*)
      integer xls(*)
      integer xadj(*)

      do i = 1, nsubg
        node = subg(i)
        mask(node) = 1
      end do

      num = 0

      do i = 1, nsubg

        node = subg(i)
c
c  For each connected component in the subgraph,
c  call FNROOT and RCM for the ordering.
c
        if ( 0 .lt. mask(node) ) then

          call fnroot ( node, xadj, adjncy, mask,
     &      nlvl, xls, perm(num+1) )

          call rcm ( node, xadj, adjncy, mask,
     &      perm(num+1), ccsize, xls )

          num = num + ccsize

          if ( nsubg .le. num ) then
            return
          end if

        end if

      end do

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
      subroutine tsfct ( nblks, xblk, father, diag, xenv, env,
     &  xnonz, nonz, nzsubs, temp, first, flag )

c*********************************************************************72
c
cc TSFCT performs the symmetric factorization of a tree-partitioned system.
c
c  Modified:
c
c    01 January 2009
c
c  input parameters 
c
c  (nblks, xblk, father) - the tree partitioning.
c
c  xenv - the envelope index vector.
c  (xnonz, nonz, nzsubs) - the off-diagonal nonzeros in
c  the original matrix.
c
c  updated parameters 
c  (diag, env) - storage arrays for the envelope of
c  the diagonal blocks of the matrix. on output,
c  contains the diagonal blocks of the factor.
c  flag - the error flag. it is set to 1 if a zero or
c  negative square root is detected during the
c  factorization.
c
c  working parameter 
c  temp - temporary array required to implement the
c  asymmetric version of the factorization.
c  first - temporary vector used to facilitate the
c  indexing to the vector nonz (or nzsubs)
c  for non-null subcolumns in off-diagonal
c  blocks.
c
      double precision diag(*)
      double precision env(*), nonz(*), temp(*), s
      integer father(*), nzsubs(*), xblk(*)
      integer first(*), xenv(*), xnonz(*)
      integer blksze, col, col1, colbeg, colend
      integer colsze, fnz, fnz1, i, flag, istrt, istop
      integer isub, j, jstop, jstrt, k, kenv, kenv0, kfathr
      integer nblks, neqns, row, rowbeg, rowend

      flag = 0
c
c  initialization.
c
      neqns = xblk(nblks+1) - 1

      do i = 1,neqns
        temp(i) = 0.0D+00
      end do

      do i = 1, neqns
        first(i) = xnonz(i)
      end do
c
c  loop through the blocks ...
c
      do 1600 k = 1, nblks
        rowbeg = xblk(k)
        rowend = xblk(k+1) - 1
        blksze = rowend - rowbeg + 1
        call esfct ( blksze, xenv(rowbeg), env,
     &    diag(rowbeg), flag )
        if ( flag .gt. 0 ) return
c
c  perform modification of the father diagonal block
c  a(father(k),father(k)) from the off-diagonal block
c  a(k,father(k)).
c
        kfathr = father(k)
        if ( kfathr .le. 0 ) go to 1600
        colbeg = xblk(kfathr)
        colend = xblk(kfathr+1) - 1
c
c  find the first and last non-null column in
c  the off-diagonal block. reset colbeg,colend.
c
        do col = colbeg, colend
          jstrt = first(col)
          jstop = xnonz(col+1) - 1
          if ( jstop .ge. jstrt .and.
     &      nzsubs(jstrt) .le. rowend ) go to 300
        end do

300     colbeg = col
        col = colend

        do col1 = colbeg, colend
          jstrt = first(col)
          jstop = xnonz(col+1) - 1
          if ( jstop .ge. jstrt .and.
     &      nzsubs(jstrt) .le. rowend ) go to 500
          col = col - 1
        end do

        colend = col
500     do 1300 col = colbeg, colend
          jstrt = first(col)
          jstop = xnonz(col+1) - 1
c
c  test for null subcolumn. fnz stores the
c  first nonzero subscript in the block column.
c
          if ( jstop .lt. jstrt ) go to 1300
          fnz = nzsubs(jstrt)
          if ( fnz .gt. rowend ) go to 1300
c
c  unpack a column in the off-diagonal block
c  and perform upper and lower solves on the
c  unpacked column.
c
          do j = jstrt, jstop
            row = nzsubs(j)
            if ( row .gt. rowend ) go to 700
            temp(row) = nonz(j)
          end do

700       colsze = rowend - fnz + 1
          call elslv ( colsze, xenv(fnz), env,
     &      diag(fnz), temp(fnz) )
          call euslv ( colsze, xenv(fnz), env,
     &      diag(fnz), temp(fnz) )
c
c  do the modification by looping through
c  the columns and forming inner products.
c
          kenv0 = xenv(col+1) - col
          do 1100 col1 = colbeg, colend
            istrt = first(col1)
            istop = xnonz(col1+1) - 1
c
c  check to see if subcolumn is null.
c
            fnz1 = nzsubs(istrt)
            if ( istop .lt. istrt .or.
     &        fnz1 .gt. rowend ) go to 1100
c
c  check if inner product should be done
c
            if ( fnz1 .lt. fnz ) go to 1100
            if ( fnz1 .eq. fnz .and.
     &        col1 .lt. col ) go to 1100

            s = 0.0D+00
            do i = istrt, istop
              isub = nzsubs(i)
              if ( isub .gt. rowend ) go to 900
              s = s + temp(isub) * nonz(i)
            end do
c
c  modify the env or the diag entry.
c
900         if ( col1 .eq. col ) go to 1000
            kenv = kenv0 + col1
            if ( col1 .gt. col )then
              kenv = xenv(col1+1) - col1 + col
            endif
            env(kenv) = env(kenv) - s
            go to 1100
1000        diag(col1) = diag(col1) - s
1100      continue
c
c  reset part of the temp vector to zero.
c
          do row = fnz, rowend
            temp(row) = 0.0D+00
          end do

1300    continue
c
c  update the first vector for columns in
c  father(k) block, so that it will index to
c  the beginning of the next off-diagonal
c  block to be considered.
c
        do 1500 col = colbeg, colend

          jstrt = first(col)
          jstop = xnonz(col+1) - 1
          if ( jstop .lt. jstrt ) go to 1500

          do 1400 j = jstrt, jstop
            row = nzsubs(j)
            if ( row .le. rowend ) go to 1400
            first(col) = j
            go to 1500
1400      continue

        first(col) = jstop + 1
1500    continue

1600  continue

      return
      end
      subroutine tsslv ( nblks, xblk, diag, xenv, env,
     &  xnonz, nonz, nzsubs, rhs, temp )

c*********************************************************************72
c
cc TSSLV solves a tree-partitioned factored system by implicit back substitution.
c
c  Modified:
c
c    01 January 2009
c
c  input parameters 
c
c  (nblks, xblk) - the partitioning.
c  (xenv, env) - envelope of the diagonal blocks.
c  (xnonz, nonz, nzsubs) - data structure for the off
c  block diagonal nonzeros.
c
c  updated parameters 
c  rhs - on input it contains the right hand vector
c  on output, the solution vector.
c
c  working vector 
c  temp - temporary vector used in back substitution.
c
      double precision diag(*)
      double precision env(*)
      double precision nonz(*)
      double precision rhs(*)
      double precision s
      double precision temp(*)

      integer nzsubs(*), xblk(*)
      integer xenv(*), xnonz(*), col, col1, col2, i, j
      integer jstop, jstrt, last, nblks, ncol, nrow, row
      integer row1, row2
c
c  forward substitution 
c
      do 400 i = 1, nblks

        row1 = xblk(i)
        row2 = xblk(i+1) - 1
        last = xnonz(row2+1)
        if ( i .eq. 1 .or. last .eq. xnonz(row1) ) go to 300
c
c  modify rhs vector by the product of the off-
c  diagonal block with the corresponding part of rhs.
c
        do 200 row = row1, row2

          jstrt = xnonz(row)
          if ( jstrt .eq. last ) go to 300
          jstop = xnonz(row+1) - 1
          if ( jstop .lt. jstrt ) go to 200
          s = 0.0D+00

          do j = jstrt, jstop
            col = nzsubs(j)
            s = s + rhs(col)*nonz(j)
          end do

          rhs(row) = rhs(row) - s

200     continue

300     nrow = row2 - row1 + 1
        call elslv ( nrow, xenv(row1), env, diag(row1),
     &    rhs(row1) )
        call euslv ( nrow, xenv(row1), env, diag(row1),
     &    rhs(row1) )

400   continue
c
c  backward solution 
c
      if ( nblks .eq. 1 ) return
      last = xblk(nblks) - 1
      do i = 1, last
        temp(i) = 0.0D+00
      end do

      i = nblks
      col1 = xblk(i)
      col2 = xblk(i+1)-1
600   if (i.eq.1)return
      last = xnonz(col2+1)
      if ( last .eq. xnonz(col1) ) go to 900
c
c  multiply off-diagonal block by the corresponding
c  part of the solution vector and store in temp.
c
      do 800 col = col1, col2

        s = rhs(col) 
        if (s.eq.0.0D+00 ) go to 800
        jstrt = xnonz(col)
        if ( jstrt .eq. last) go to 900
        jstop = xnonz(col+1) - 1
        if ( jstop .lt. jstrt ) go to 800

        do j = jstrt, jstop
          row = nzsubs(j)
          temp(row) = temp(row) + s*nonz(j)
        end do

800   continue

900   i = i - 1
      col1 = xblk(i)
      col2 = xblk(i+1) - 1
      ncol = col2 - col1 + 1
      call elslv ( ncol, xenv(col1), env,
     &  diag(col1), temp(col1) )

      call euslv ( ncol, xenv(col1), env, diag(col1),
     &  temp(col1) )

      do j = col1, col2
        rhs(j) = rhs(j) - temp(j)
      end do

      go to 600
      end
