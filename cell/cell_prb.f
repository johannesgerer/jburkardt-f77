      program main

c*********************************************************************72
c
cc MAIN tests CELL.
c
c  Discussion:
c
c    An R8CVV is a "cell vector of vectors" of R8's.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    16 November 2012
c
c  Author:
c
c    John Burkardt
c
      implicit none

      call timestamp ( )
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'CELL_PRB:'
      write ( *, '(a)' ) '  FORTRAN77 version'
      write ( *, '(a)' ) '  Test the CELL library.'

      call test01 ( )
c
c  Terminate.
c
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'CELL_PRB:'
      write ( *, '(a)' ) '  Normal end of execution.'
      write ( *, '(a)' ) ' '
      call timestamp ( )

      stop
      end
      subroutine test01 ( )

c*********************************************************************72
c
cc TEST01 stores some of Pascal's triangle in an R8CVV.
c
c  Discussion:
c
c    An R8CVV is a "cell array vector of vectors" of R8's.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    16 November 2012
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer m
      parameter ( m = 5 )
      integer nn
      parameter ( nn = 4 )

      double precision a(30)
      double precision ai(8)
      double precision aij
      integer col
      integer i
      integer in(nn)
      integer j
      integer jn(nn)
      integer mn
      integer nr(m)
      integer nv
      integer roff(m+1)
      integer row
      double precision vn(nn)

      save in
      save jn
      save nr

      data in / 1, 2, 5, 5 /
      data jn / 2, 3, 4, 8 /
      data nr / 4, 5, 6, 7, 8 /

      write ( *, '(a)' ) ''
      write ( *, '(a)' ) 'TEST01'
      write ( *, '(a)' ) 
     &  '  Use a cell array (vector of vectors) to store rows 3:7'
      write ( *, '(a)' ) '  of Pascal''s triangle.'

      call i4vec_print ( m, nr, '  The row sizes:' )
c
c  From the NR information:
c  * determine the total size, MN
c
      call r8cvv_size ( m, nr, mn )
      write ( *, '(a)' ) ''
      write ( *, '(a,i4)' ) '  The storage for the cell array is ', mn
c
c  Zero out the cell array.
c
      do i = 1, mn
        a(i) = 0.0D+00
      end do
c
c  From the NR information:
c  * determine the offsets.
c
      call r8cvv_offset ( m, nr, roff )
      call i4vec_print ( m + 1, roff, '  The row offsets:' )
c
c  Rows 1 through 5 of A will contain rows 3 through 7 of Pascal's triangle.
c  Set these values one row at a time.
c
      ai(1) = 1.0D+00

      do row = 1, 7

        col = row + 1
        ai(col) = ai(col-1)
        do col = row, 2, -1
          ai(col) = ai(col) + ai(col-1)
        end do

        if ( 3 .le. row ) then
          i = row - 2
          call r8cvv_rset ( mn, a, m, roff, i, ai )
        end if

      end do
c
c  Print the cell array.
c
      call r8cvv_print ( mn, a, m, roff, 
     &  '  Rows 3:7 of Pascal''s Triangle:' )
c
c  Retrieve the entry from cell array row 3, column 4:
c
      i = 3
      j = 4
      call r8cvv_iget ( mn, a, m, roff, i, j, aij )
      write ( *, '(a)' ) ' '
      write ( *, '(a,i1,a,i1,a,g14.6)' ) '  A(', i, ',', j, ') = ', aij
c
c  Retrieve row 4:
c
      i = 4
      call r8cvv_rget ( mn, a, m, roff, i, ai )
      nv = roff(i+1) - roff(i)
      call r8vec_transpose_print ( nv, ai, '  A(4,*):' )
c
c  Retrieve a list of entries.
c
      call r8cvv_nget ( mn, a, m, roff, nn, in, jn, vn )
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Retrieve an arbitrary list of items:'
      write ( *, '(a)' ) ' '
      do i = 1, nn
        write ( *, '(a,i1,a,i1,a,g14.6)' ) 
     &    '  A(', in(i), ',', jn(i), ') = ', vn(i)
      end do

      return
      end

