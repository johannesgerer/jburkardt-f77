      program main

c*********************************************************************72
c
cc MAIN is the main program for SPARSEPAK_PRB.
c
c  Discussion:
c
c    SPARSEPAK_PRB runs the SPARSEPAK tests.
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
      implicit none

      call timestamp ( )
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'SPARSEPAK_PRB'
      write ( *, '(a)' ) '  FORTRAN77 version.'
      write ( *, '(a)' ) '  Test the SPARSEPAK library.'

      call test05 ( )
c
c  Terminate.
c
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'SPARSEPAK_PRB'
      write ( *, '(a)' ) '  Normal end of execution.'
      write ( *, '(a)' ) ' '
      call timestamp ( )

      stop
      end
      subroutine test05 ( )

c*********************************************************************72
c
cc TEST05 tests the RQT method.
c
c  Discussion:
c
c    This example involves the following equations:
c
c    2*x1  - x10       = 0
c    2*x2  - x9  - x10 = 0
c    2*x3  - x8  - x9  = 0
c    2*x4  - x7  - x8  = 0
c    2*x5  - x6  - x7  = 0
c    2*x6  - x5        = 11
c    2*x7  - x4  - x5  = 0
c    2*x8  - x3  - x4  = 0
c    2*x9  - x2  - x3  = 0
c    2*x10 - x1  - x2  = 0
c
c    with solution (1,3,5,7,9,10,8,6,4,2).
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
      implicit none

      integer maxadj
      integer maxblk
      integer maxenv
      integer maxnon
      integer n

      parameter ( maxadj = 300 )
      parameter ( maxblk = 10 )
      parameter ( maxenv = 300 )
      parameter ( maxnon = 300 )
      parameter ( n = 10 )

      integer adj(maxadj)
      integer bnum(n)
      double precision diag(n)
      double precision env(maxenv)
      integer father(n)
      integer i
      integer iband
      integer env_size
      integer first(n)
      integer ierror
      integer invprm(n)
      integer ls(n)
      integer mask(n)
      integer nadj
      integer nblks
      integer nodlvl(n)
      integer nofnz
      double precision nonz(maxnon)
      integer nzsub(maxnon)
      integer perm(n)
      double precision rhs(n)
      integer subg(n)
      double precision temp(n)
      integer xadj(n+1)
      integer xblk(maxblk+1)
      integer xenv(n+1)
      integer xls(n)
      integer xnonz(n+1)

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST05'
      write ( *, '(a)' ) '  Use the RQT method.'
      write ( *, '(a)' ) ' '
      write ( *, '(a,i6)' ) '  The matrix order is N = ', n
c
c  Initialize the permutation vectors.
c
      do i = 1, n
        perm(i) = i
      end do

      do i = 1, n
        invprm(i) = i
      end do
c
c  Store the adjacency information.
c
      call adj_set_5 ( adj, maxadj, nadj, n, xadj )

      write ( *, '(a)' ) ' '
      write ( *, '(a,i12)' )
     &  '  Number of adjacency entries NADJ = ', nadj
c
c  Display the adjacency information.
c
      call adj_print ( n, nadj, xadj, adj )
c
c  Get a picture of the matrix.
c
      call adj_show ( adj, iband, invprm, nadj, n, perm, xadj )
c
c  Generate the RQT ordering.
c
      call genrqt ( n, xadj, adj, nblks, xblk, perm, xls, ls, nodlvl )

      write ( *, '(a)' ) ' '
      write ( *, '(a,i12)' )
     &  '  After GENRQT, the  number of blocks is ', nblks
c
c  Get a picture of the matrix.
c
      call adj_show ( adj, iband, invprm, nadj, n, perm, xadj )

      write ( *, '(a)' ) ' '
      write ( *, '(a,i12)' ) '  The envelope size is ', env_size

      call bshufl ( xadj, adj, perm, nblks, xblk, bnum, mask,
     &  subg, xls )
c
c  Compute the inverse ordering.
c
      call perm_inverse ( n, perm, invprm )
c
c  Print orderings
c
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '    I    Perm(I)   InvPerm(I)'
      write ( *, '(a)' ) ' '

      do i = 1, n
        write ( *, '(3i6)' ) i, perm(i), invprm(i)
      end do
c
c  Get a picture of the matrix.
c
      call adj_show ( adj, iband, invprm, nadj, n, perm, xadj )
c
c  Determine the quotient tree adjacency structure.
c
      call fntadj ( xadj, adj, perm, invprm, nblks, xblk, father, bnum )
c
c  Determine the envelope index vector.
c
      call fntenv ( xadj, adj, perm, invprm, nblks, xblk, xenv,
     &    env_size )

      write ( *, '(a)' ) ' '
      write ( *, '(a,i12)' )
     &  '  The reordered envelope size is ', env_size

      nofnz = maxnon

      call fnofnz ( xadj, adj, perm, invprm, nblks, xblk, xnonz,
     &  nzsub, nofnz )
c
c  Set RHS, DIAG, ENV.
c
      call setsy5 ( diag, env, adj, invprm, xadj, maxadj, maxenv,
     &  n, nonz, nzsub, rhs, xenv, xnonz )
c
c  Factor the system.
c
      write ( *, * ) 'NBLKS = ', nblks

      call tsfct ( nblks, xblk, father, diag, xenv, env, xnonz,
     &   nonz, nzsub, temp, first, ierror )

      if ( ierror .ne. 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'TEST05 - Fatal error!'
        write ( *, '(a)' ) '  The matrix is not positive definite.'
        write ( *, '(a,i8)' ) '  TSFCT returns IERROR = ', ierror
        return
      end if
c
c  Solve the system.
c
      call tsslv ( nblks, xblk, diag, xenv, env, xnonz, nonz,
     &   nzsub, rhs, temp )
c
c  Unpermute the solution.
c
      call permrv ( n, rhs, perm )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Solution:'
      write ( *, '(a)' ) ' '
      do i = 1, n
        write ( *, '(5g14.6)' ) rhs(i)
      end do

      return
      end
      subroutine setsy5 ( diag, env, adj, invprm, xadj, maxadj,
     &  maxenv, n, nonz, nzsub, rhs, xenv, xnonz )

c*********************************************************************72
c
cc SETSY5 stores the numerical values defining problem 5.
c
c  Discussion:
c
c    There is only one nonzero right hand side entry.
c    The matrix diagonal entries are all 2.
c    The nonzero offdiagonal entries are all -1.
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
      implicit none

      integer maxadj
      integer maxenv
      integer n

      integer adj(maxadj)
      double precision diag(n)
      double precision env(maxenv)
      integer i
      integer invprm(n)
      integer isub
      integer j
      integer xadj(n+1)
      integer jsub
      double precision nonz(*)
      integer nzsub(*)
      double precision rhs(n)
      double precision value
      integer xenv(n+1)
      integer xnonz(n+1)
c
c  Zero out storage.
c
      do i = 1, n
        rhs(i) = 0.0D+00
      end do

      do i = 1, n
        diag(i) = 0.0D+00
      end do

      do i = 1, maxenv
        env(i) = 0.0D+00
      end do
c
c  Set the nonzero elements of the right hand side vector.
c
      isub = 6
      value = 11.0D+00

      call addrhs ( invprm, isub, n, rhs, value )
c
c  Set the diagonal entries of the matrix.
c
      do i = 1, n
        diag(i) = 2.0D+00
      end do
c
c  Set the off diagonal terms of the matrix.
c
      do i = 1, n

        isub = i

        do j = xadj(i), xadj(i+1) - 1

          jsub = adj(j)
          value = - 1.0D+00

          call addrqt ( isub, jsub, value, invprm, diag,
     &        xenv, env, xnonz, nonz, nzsub, n )

        end do

      end do

      return
      end
      subroutine adj_set_5 ( adj, maxadj, nadj, n, xadj )

c*********************************************************************72
c
cc ADJ_SET_5 sets up the adjacency structure for problem 5.
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
      implicit none

      integer nonz
      parameter ( nonz = 19 )

      integer maxadj
      integer n

      integer adj(maxadj)
      integer i
      integer ilist(nonz)
      integer jlist(nonz)
      integer nadj
      integer xadj(n+1)

      data ilist /0,1,2,2,3,3,4,4,5,5,6,7,7,8,8,9,9,10,10/
      data jlist /0,10,10,9,9,8,7,8,7,6,5,4,5,3,4,2,3,1,2/

      do i = 1, nonz
        call setadj ( adj, ilist(i), jlist(i), maxadj, nadj, n, xadj )
      end do

      return
      end
