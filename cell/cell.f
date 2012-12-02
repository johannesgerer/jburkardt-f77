      subroutine i4vec_print ( n, a, title )

c*********************************************************************72
c
cc I4VEC_PRINT prints an I4VEC.
c
c  Discussion:
c
c    An I4VEC is a vector of integer values.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    27 November 2006
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the number of components of the vector.
c
c    Input, integer A(N), the vector to be printed.
c
c    Input, character*(*) TITLE, a title.
c
      implicit none

      integer n

      integer a(n)
      integer i
      character*(*) title

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) trim ( title )
      write ( *, '(a)' ) ' '
      do i = 1, n
        write ( *, '(2x,i8,a,1x,i12)' ) i, ':', a(i)
      end do

      return
      end
      subroutine r8cvv_iget ( mn, a, m, roff, i, j, aij )

c*********************************************************************72
c
cc R8CVV_IGET gets item J from row I in an R8CVV.
c
c  Discussion:
c
c    An R8CVV is a "vector of vectors" of R8's.
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
c  Parameters:
c
c    Input, integer MN, the size of the cell array.
c
c    Input, double precision A(MN), the cell array.
c
c    Input, integer M, the number of rows in the array.
c
c    Input, integer ROFF(M+1), the row offsets.
c
c    Input, integer I, the row of the item.
c    1 <= I <= M.
c
c    Input, integer J, the column of the item.
c    1 <= J <= NR(I).
c
c    Output, double precision AIJ, the value of item A(I,J).
c
      implicit none

      integer m
      integer mn

      double precision a(mn)
      double precision aij
      integer i
      integer j
      integer k
      integer roff(m+1)

      k = roff(i) + j
      aij = a(k)

      return
      end
      subroutine r8cvv_iinc ( mn, a, m, roff, i, j, daij )

c*********************************************************************72
c
cc R8CVV_IINC increments item J from row I in an R8CVV.
c
c  Discussion:
c
c    An R8CVV is a "vector of vectors" of R8's.
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
c  Parameters:
c
c    Input, integer MN, the size of the cell array.
c
c    Input/output, double precision A(MN), the cell array.
c
c    Input, integer M, the number of rows in the array.
c
c    Input, integer ROFF(M+1), the row offsets.
c
c    Input, integer I, the row of the item.
c    1 <= I <= M.
c
c    Input, integer J, the column of the item.
c    1 <= J <= NR(I).
c
c    Input, double precision DAIJ, the increment to the value of item A(I,J).
c
      implicit none

      integer m
      integer mn

      double precision a(mn)
      double precision daij
      integer i
      integer j
      integer k
      integer roff(m+1)

      k = roff(i) + j
      a(k) = a(k) + daij

      return
      end
      subroutine r8cvv_iset ( mn, a, m, roff, i, j, aij )

c*********************************************************************72
c
cc R8CVV_ISET sets item J from row I in an R8CVV.
c
c  Discussion:
c
c    An R8CVV is a "vector of vectors" of R8's.
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
c  Parameters:
c
c    Input, integer MN, the size of the cell array.
c
c    Input/output, double precision A(MN), the cell array.
c
c    Input, integer M, the number of rows in the array.
c
c    Input, integer ROFF(M+1), the row offsets.
c
c    Input, integer I, the row of the item.
c    1 <= I <= M.
c
c    Input, integer J, the column of the item.
c    1 <= J <= NR(I).
c
c    Input, double precision AIJ, the new value of item A(I,J).
c
      implicit none

      integer m
      integer mn

      double precision a(mn)
      double precision aij
      integer i
      integer j
      integer k
      integer roff(m+1)

      k = roff(i) + j
      a(k) = aij

      return
      end
      subroutine r8cvv_nget ( mn, a, m, roff, nn, in, jn, vn )

c*********************************************************************72
c
cc R8CVV_NGET gets N items JN(*) from row IN(*) in an R8CVV.
c
c  Discussion:
c
c    An R8CVV is a "vector of vectors" of R8's.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    15 November 2012
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer MN, the size of the cell array.
c
c    Input, double precision A(MN), the cell array.
c
c    Input, integer M, the number of rows in the array.
c
c    Input, integer ROFF(M+1), the row offsets.
c
c    Input, integer NN, the number of items.
c
c    Input, integer IN(NN), the rows of the items.
c    1 <= IN(*) <= M.
c
c    Input, integer JN(NN), the columns of the items.
c    1 <= JN(*) <= NR(IN(*)).
c
c    Output, double precision VN(NN), the value of items A(IN(*),JN(*)).
c
      implicit none

      integer m
      integer mn
      integer nn

      double precision a(mn)
      integer i
      integer in(nn)
      integer jn(nn)
      integer k
      integer roff(m+1)
      double precision vn(nn)

      do i = 1, nn
        k = roff(in(i)) + jn(i)
        vn(i) = a(k)
      end do

      return
      end
      subroutine r8cvv_ninc ( mn, a, m, roff, nn, in, jn, dvn )

c*********************************************************************72
c
cc R8CVV_NINC increments items JN(*) from row IN(*) in an R8CVV.
c
c  Discussion:
c
c    An R8CVV is a "vector of vectors" of R8's.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    15 November 2012
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer MN, the size of the cell array.
c
c    Input/output, double precision A(MN), the cell array.
c
c    Input, integer M, the number of rows in the array.
c
c    Input, integer ROFF(M+1), the row offsets.
c
c    Input, integer NN, the number of items.
c
c    Input, integer IN(NN), the rows of the items.
c    1 <= IN(*) <= M.
c
c    Input, integer JN(NN), the columns of the items.
c    1 <= JN(*) <= NR(IN(*)).
c
c    Input, double precision DVN(NN), the increments of items A(IN(*),JN(*)).
c
      implicit none

      integer m
      integer mn
      integer nn

      double precision a(mn)
      double precision dvn(nn)
      integer i
      integer in(nn)
      integer jn(nn)
      integer k
      integer roff(m+1)

      do i = 1, nn
        k = roff(in(i)) + jn(i)
        a(k) = a(k) + dvn(i)
      end do

      return
      end
      subroutine r8cvv_nset ( mn, a, m, roff, nn, in, jn, vn )

c*********************************************************************72
c
cc R8CVV_NSET sets items JN(*) from row IN(*) in an R8CVV.
c
c  Discussion:
c
c    An R8CVV is a "vector of vectors" of R8's.
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
c  Parameters:
c
c    Input, integer MN, the size of the cell array.
c
c    Input/output, double precision A(MN), the cell array.
c
c    Input, integer M, the number of rows in the array.
c
c    Input, integer ROFF(M+1), the row offsets.
c
c    Input, integer NN, the number of items.
c
c    Input, integer IN(NN), the rows of the items.
c    1 <= IN(*) <= M.
c
c    Input, integer JN(NN), the columns of the items.
c    1 <= JN(*) <= NR(IN(*)).
c
c    Input, double precision VN(NN), the new value of items A(IN(*),JN(*)).
c
      implicit none

      integer m
      integer mn
      integer nn

      double precision a(mn)
      integer i
      integer in(nn)
      integer jn(nn)
      integer k
      integer roff(m+1)
      double precision vn(nn)

      do i = 1, nn
        k = roff(in(i)) + jn(i)
        a(k) = vn(i)
      end do

      return
      end
      subroutine r8cvv_offset ( m, nr, roff )

c*********************************************************************72
c
cc R8CVV_OFFSET determines the row offsets of an R8CVV.
c
c  Discussion:
c
c    An R8CVV is a "vector of vectors" of R8's.
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
c  Parameters:
c
c    Input, integer M, the number of rows in the array.
c
c    Input, integer NR(M), the row sizes.
c
c    Output, integer ROFF(M+1), the row offsets.
c
      implicit none

      integer m

      integer i
      integer roff(m+1)
      integer nr(m)

      roff(1) = 0
      do i = 1, m
        roff(i+1) = roff(i) + nr(i)
      end do

      return
      end
      subroutine r8cvv_print ( mn, a, m, roff, title )

c*********************************************************************72
c
cc R8CVV_PRINT prints an R8CVV.
c
c  Discussion:
c
c    An R8CVV is a "vector of vectors" of R8's.
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
c  Parameters:
c
c    Input, integer MN, the size of the cell array.
c
c    Input, double precision A(MN), the cell array.
c
c    Input, integer M, the number of rows in the array.
c
c    Input, integer ROFF(M+1), the row offsets.
c
c    Input, character ( len = * ) TITLE, a title.
c
      implicit none

      integer m
      integer mn

      double precision a(mn)
      integer i
      integer k1
      integer k2
      integer khi
      integer klo
      integer roff(m+1)
      character ( len = * ) title

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) trim ( title )
      write ( *, '(a)' ) ' '

      do i = 1, m

        k1 = roff(i) + 1
        k2 = roff(i+1)

        do klo = k1, k2, 5
          khi = min ( klo + 5 - 1, k2 )
          if ( klo .eq. k1 ) then
            write ( *, '(i5,2x, 5g14.6)' ) i, a(klo:khi)
          else
            write ( *, '(5x,2x, 5g14.6)' )    a(klo:khi)
          end if
        end do

      end do

      return
      end
      subroutine r8cvv_rget ( mn, a, m, roff, i, ai )

c*********************************************************************72
c
cc R8CVV_RGET gets row I from an R8CVV.
c
c  Discussion:
c
c    An R8CVV is a "vector of vectors" of R8's.
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
c  Parameters:
c
c    Input, integer MN, the size of the cell array.
c
c    Input, double precision A(MN), the cell array.
c
c    Input, integer M, the number of rows in the array.
c
c    Input, integer ROFF(M+1), the row offsets.
c
c    Input, integer I, the row.
c    1 <= I <= M.
c
c    Output, double precision AI(NR(I)), the value of A(I,*).
c
      implicit none

      integer m
      integer mn

      double precision a(mn)
      double precision ai(*)
      integer i
      integer j
      integer k
      integer k1
      integer k2
      integer nv
      integer roff(m+1)

      k1 = roff(i) + 1
      k2 = roff(i+1)
      nv = k2 + 1 - k1
      do j = 1, nv
        k = roff(i) + j
        ai(j) = a(k)
      end do

      return
      end
      subroutine r8cvv_rinc ( mn, a, m, roff, i, dai )

c*********************************************************************72
c
cc R8CVV_RINC increments row I in an R8CVV.
c
c  Discussion:
c
c    An R8CVV is a "vector of vectors" of R8's.
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
c  Parameters:
c
c    Input, integer MN, the size of the cell array.
c
c    Input/output, double precision A(MN), the cell array.
c
c    Input, integer M, the number of rows in the array.
c
c    Input, integer ROFF(M+1), the row offsets.
c
c    Input, integer I, the row.
c    1 <= I <= M.
c
c    Input, double precision DAI(NR(I)), the increment for A(I,*).
c
      implicit none

      integer m
      integer mn

      double precision a(mn)
      double precision dai(*)
      integer i
      integer j
      integer k
      integer k1
      integer k2
      integer nv
      integer roff(m+1)

      k1 = roff(i) + 1
      k2 = roff(i+1)
      nv = k2 + 1 - k1
      do j = 1, nv
        k = roff(i) + j
        a(k) = a(k) + dai(j)
      end do

      return
      end
      subroutine r8cvv_rset ( mn, a, m, roff, i, ai )

c*********************************************************************72
c
cc R8CVV_RSET sets row I from an R8CVV.
c
c  Discussion:
c
c    An R8CVV is a "vector of vectors" of R8's.
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
c  Parameters:
c
c    Input, integer MN, the size of the cell array.
c
c    Input/output, double precision A(MN), the cell array.
c
c    Input, integer M, the number of rows in the array.
c
c    Input, integer ROFF(M+1), the row offsets.
c
c    Input, integer I, the row.
c    1 <= I <= M.
c
c    Input, double precision AI(NR(I)), the new value of A(I,*).
c
      implicit none

      integer m
      integer mn

      double precision a(mn)
      double precision ai(*)
      integer i
      integer j
      integer k
      integer k1
      integer k2
      integer nv
      integer roff(m+1)

      k1 = roff(i) + 1
      k2 = roff(i+1)
      nv = k2 + 1 - k1
      do j = 1, nv
        k = roff(i) + j
        a(k) = ai(j)
      end do

      return
      end
      subroutine r8cvv_size ( m, nr, mn )

c*********************************************************************72
c
cc R8CVV_SIZE determines the size of an R8CVV.
c
c  Discussion:
c
c    An R8CVV is a "vector of vectors" of R8's.
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
c  Parameters:
c
c    Input, integer M, the number of rows in the array.
c
c    Input, integer NR(M), the size of each row.
c
c    Output, integer MN, the size of the cell array.
c
      implicit none

      integer m

      integer i
      integer mn
      integer nr(m)

      mn = 0
      do i = 1, m
        mn = mn + nr(i)
      end do

      return
      end
      subroutine r8vec_print ( n, a, title )

c*********************************************************************72
c
cc R8VEC_PRINT prints an R8VEC.
c
c  Discussion:
c
c    An R8VEC is a vector of R8's.
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
c    Input, integer N, the number of components of the vector.
c
c    Input, double precision A(N), the vector to be printed.
c
c    Input, character * ( * ) TITLE, a title.
c
      implicit none

      integer n

      double precision a(n)
      integer i
      character ( len = * ) title

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) trim ( title )
      write ( *, '(a)' ) ' '
      do i = 1, n
        write ( *, '(2x,i8,a,1x,g16.8)' ) i, ':', a(i)
      end do

      return
      end
      subroutine r8vec_transpose_print ( n, a, title )

c*********************************************************************72
c
cc R8VEC_TRANSPOSE_PRINT prints an R8VEC "transposed".
c
c  Discussion:
c
c    An R8VEC is a vector of R8's.
c
c  Example:
c
c    A = (/ 1.0, 2.1, 3.2, 4.3, 5.4, 6.5, 7.6, 8.7, 9.8, 10.9, 11.0 /)
c    TITLE = 'My vector:  '
c
c    My vector:
c
c        1.0    2.1    3.2    4.3    5.4
c        6.5    7.6    8.7    9.8   10.9
c       11.0
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    12 November 2010
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the number of components of the vector.
c
c    Input, double precision A(N), the vector to be printed.
c
c    Input, character * ( * ) TITLE, a title.
c
      implicit none

      integer n

      double precision a(n)
      integer ihi
      integer ilo
      character * ( * )  title

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) trim ( title )
      write ( *, '(a)' ) ' '
      do ilo = 1, n, 5
        ihi = min ( ilo + 5 - 1, n )
        write ( *, '(5g14.6)' ) a(ilo:ihi)
      end do

      return
      end
      subroutine timestamp ( )

c*********************************************************************72
c
cc TIMESTAMP prints out the current YMDHMS date as a timestamp.
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
