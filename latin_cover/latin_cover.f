      function i4_modp ( i, j )

c*********************************************************************72
c
cc I4_MODP returns the nonnegative remainder of integer division.
c
c  Discussion:
c
c    If
c      NREM = I4_MODP ( I, J )
c      NMULT = ( I - NREM ) / J
c    then
c      I = J * NMULT + NREM
c    where NREM is always nonnegative.
c
c    The MOD function computes a result with the same sign as the
c    quantity being divided.  Thus, suppose you had an angle A,
c    and you wanted to ensure that it was between 0 and 360.
c    Then mod(A,360) would do, if A was positive, but if A
c    was negative, your result would be between -360 and 0.
c
c    On the other hand, I4_MODP(A,360) is between 0 and 360, always.
c
c  Example:
c
c        I     J     MOD I4_MODP    Factorization
c
c      107    50       7       7    107 =  2 *  50 + 7
c      107   -50       7       7    107 = -2 * -50 + 7
c     -107    50      -7      43   -107 = -3 *  50 + 43
c     -107   -50      -7      43   -107 =  3 * -50 + 43
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    30 December 2006
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer I, the number to be divided.
c
c    Input, integer J, the number that divides I.
c
c    Output, integer I4_MODP, the nonnegative remainder when I is
c    divided by J.
c
      implicit none

      integer i
      integer i4_modp
      integer j
      integer value

      if ( j .eq. 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'I4_MODP - Fatal error!'
        write ( *, '(a,i8)' ) '  Illegal divisor J = ', j
        stop
      end if

      value = mod ( i, j )

      if ( value .lt. 0 ) then
        value = value + abs ( j )
      end if

      i4_modp = value

      return
      end
      function i4_uniform ( a, b, seed )

c*********************************************************************72
c
cc I4_UNIFORM returns a scaled pseudorandom I4.
c
c  Discussion:
c
c    An I4 is an integer value.
c
c    The pseudorandom number should be uniformly distributed
c    between A and B.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    12 November 2006
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Paul Bratley, Bennett Fox, Linus Schrage,
c    A Guide to Simulation,
c    Springer Verlag, pages 201-202, 1983.
c
c    Pierre L'Ecuyer,
c    Random Number Generation,
c    in Handbook of Simulation,
c    edited by Jerry Banks,
c    Wiley Interscience, page 95, 1998.
c
c    Bennett Fox,
c    Algorithm 647:
c    Implementation and Relative Efficiency of Quasirandom
c    Sequence Generators,
c    ACM Transactions on Mathematical Software,
c    Volume 12, Number 4, pages 362-376, 1986.
c
c    Peter Lewis, Allen Goodman, James Miller
c    A Pseudo-Random Number Generator for the System/360,
c    IBM Systems Journal,
c    Volume 8, pages 136-143, 1969.
c
c  Parameters:
c
c    Input, integer A, B, the limits of the interval.
c
c    Input/output, integer SEED, the "seed" value, which should NOT be 0.
c    On output, SEED has been updated.
c
c    Output, integer I4_UNIFORM, a number between A and B.
c
      implicit none

      integer a
      integer b
      integer i4_uniform
      integer k
      real r
      integer seed
      integer value

      if ( seed .eq. 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'I4_UNIFORM - Fatal error!'
        write ( *, '(a)' ) '  Input value of SEED = 0.'
        stop
      end if

      k = seed / 127773

      seed = 16807 * ( seed - k * 127773 ) - k * 2836

      if ( seed .lt. 0 ) then
        seed = seed + 2147483647
      end if

      r = real ( seed ) * 4.656612875E-10
c
c  Scale R to lie between A-0.5 and B+0.5.
c
      r = ( 1.0E+00 - r ) * ( real ( min ( a, b ) ) - 0.5E+00 )
     &  +             r   * ( real ( max ( a, b ) ) + 0.5E+00 )
c
c  Use rounding to convert R to an integer between A and B.
c
      value = nint ( r )

      value = max ( value, min ( a, b ) )
      value = min ( value, max ( a, b ) )

      i4_uniform = value

      return
      end
      function i4_wrap ( ival, ilo, ihi )

c*********************************************************************72
c
cc I4_WRAP forces an I4 to lie between given limits by wrapping.
c
c  Example:
c
c    ILO = 4, IHI = 8
c
c    I  Value
c
c    -2     8
c    -1     4
c     0     5
c     1     6
c     2     7
c     3     8
c     4     4
c     5     5
c     6     6
c     7     7
c     8     8
c     9     4
c    10     5
c    11     6
c    12     7
c    13     8
c    14     4
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    31 December 2006
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer IVAL, an integer value.
c
c    Input, integer ILO, IHI, the desired bounds for the integer value.
c
c    Output, integer I4_WRAP, a "wrapped" version of IVAL.
c
      implicit none

      integer i4_modp
      integer i4_wrap
      integer ihi
      integer ilo
      integer ival
      integer jhi
      integer jlo
      integer value
      integer wide

      jlo = min ( ilo, ihi )
      jhi = max ( ilo, ihi )

      wide = jhi - jlo + 1

      if ( wide .eq. 1 ) then
        value = jlo
      else
        value = jlo + i4_modp ( ival - jlo, wide )
      end if

      i4_wrap = value

      return
      end
      subroutine i4block_print ( l, m, n, a, title )

c*********************************************************************72
c
cc I4BLOCK_PRINT prints an I4BLOCK.
c
c  Discussion:
c
c    An I4BLOCK is a 3D array of I4 values.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    25 June 2012
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer L, M, N, the dimensions of the block.
c
c    Input, integer A(L,M,N), the matrix to be printed.
c
c    Input, character * ( * ) TITLE, a title.
c
      implicit none

      integer l
      integer m
      integer n

      integer a(l,m,n)
      integer i
      integer j
      integer jhi
      integer jlo
      integer k
      character * ( * ) title

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) trim ( title )

      do k = 1, n

        write ( *, '(a)' ) ' '
        write ( *, '(a,i8)' ) '  K = ', k

        do jlo = 1, m, 10
          jhi = min ( jlo + 10 - 1, m )
          write ( *, '(a)' ) ' '
          write ( *, '(8x,a2,10(2x,i6))' ) 'J:', ( j, j = jlo, jhi )
          write ( *, '(7x,a2)' ) 'I:'
          do i = 1, l
            write ( *, '(2x,i6,a1,1x,10(2x,i6))' ) 
     &        i, ':', a(i,jlo:jhi,k)
          end do
        end do

      end do

      return
      end
      subroutine i4mat_print ( m, n, a, title )

c*********************************************************************72
c
cc I4MAT_PRINT prints an I4MAT.
c
c  Discussion:
c
c    An I4MAT is an array of I4's.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    30 June 2003
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer M, the number of rows in A.
c
c    Input, integer N, the number of columns in A.
c
c    Input, integer A(M,N), the matrix to be printed.
c
c    Input, character*(*) TITLE, a title.
c
      implicit none

      integer m
      integer n

      integer a(m,n)
      integer ihi
      integer ilo
      integer jhi
      integer jlo
      character*(*) title

      ilo = 1
      ihi = m
      jlo = 1
      jhi = n

      call i4mat_print_some ( m, n, a, ilo, jlo, ihi, jhi, title )

      return
      end
      subroutine i4mat_print_some ( m, n, a, ilo, jlo, ihi, jhi, title )

c*********************************************************************72
c
cc I4MAT_PRINT_SOME prints some of an I4MAT.
c
c  Discussion:
c
c    An I4MAT is an array of I4's.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    04 November 2003
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer M, N, the number of rows and columns.
c
c    Input, integer A(M,N), an M by N matrix to be printed.
c
c    Input, integer ILO, JLO, the first row and column to print.
c
c    Input, integer IHI, JHI, the last row and column to print.
c
c    Input, character*(*) TITLE, a title.
c
      implicit none

      integer incx
      parameter ( incx = 10 )
      integer m
      integer n

      integer a(m,n)
      character*(8) ctemp(incx)
      integer i
      integer i2hi
      integer i2lo
      integer ihi
      integer ilo
      integer inc
      integer j
      integer j2
      integer j2hi
      integer j2lo
      integer jhi
      integer jlo
      character*(*) title

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) trim ( title )

      if ( m .le. 0 .or. n .le. 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) '  (None)'
        return
      end if

      do j2lo = max ( jlo, 1 ), min ( jhi, n ), incx

        j2hi = j2lo + incx - 1
        j2hi = min ( j2hi, n )
        j2hi = min ( j2hi, jhi )

        inc = j2hi + 1 - j2lo

        write ( *, '(a)' ) ' '

        do j = j2lo, j2hi
          j2 = j + 1 - j2lo
          write ( ctemp(j2), '(i8)' ) j
        end do

        write ( *, '(''  Col '',10a8)' ) ( ctemp(j), j = 1, inc )
        write ( *, '(a)' ) '  Row'
        write ( *, '(a)' ) ' '

        i2lo = max ( ilo, 1 )
        i2hi = min ( ihi, m )

        do i = i2lo, i2hi

          do j2 = 1, inc

            j = j2lo - 1 + j2

            write ( ctemp(j2), '(i8)' ) a(i,j)

          end do

          write ( *, '(i5,a,10a8)' ) i, ':', ( ctemp(j), j = 1, inc )

        end do

      end do

      return
      end
      subroutine latin_cover ( n, p, a )

c*********************************************************************72
c
cc LATIN_COVER returns a 2D Latin Square Covering.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    24 June 2012
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the number of points.
c
c    Input, integer P(N), a permutation which describes the
c    first Latin square.
c
c    Output, integer A(N,N), the Latin cover.  A(I,J) = K
c    means that (I,J) is one element of the K-th Latin square.
c
      implicit none

      integer n

      integer a(n,n)
      integer i
      integer i4_wrap
      integer ik
      integer k
      integer p(n)

      call perm_check ( n, p )

      do i = 1, n
        do k = 1, n
          ik = i4_wrap ( i + k - 1, 1, n )
          a(i,p(ik)) = k
        end do
      end do

      return
      end
      subroutine latin_cover_2d ( n, p, a )

c*********************************************************************72
c
cc LATIN_COVER_2D returns a 2D Latin Square Covering.
c
c  Discussion:
c
c    This procedure has a chance of being extended to M dimensions.
c
c    A basic solution is computed, and the user is permitted to permute
c    both the I and J coordinates.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    24 June 2012
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the number of points.
c
c    Input, integer P(2,N), permutations to be applied
c    to the spatial dimensions.
c
c    Output, integer A(N,N), the Latin cover.  A(I,J) = K
c    means that (I,J) is one element of the K-th Latin square.
c
      implicit none

      integer n

      integer a(n,n)
      integer b(n,n)
      integer base
      integer i
      integer i4_wrap
      integer j
      integer p(2,n)

      base = 1

      call perm_check ( n, p(1,1:n) )
      call perm_check ( n, p(2,1:n) )
c
c  Set up the basic solution.
c
      do i = 1, n
        do j = 1, n
          a(i,j) = i4_wrap ( i - j + base, 0 + base, n - 1 + base )
        end do
      end do
c
c  Apply permutation to dimension I.
c
      do i = 1, n
        b(p(1,i),1:n) = a(i,1:n) 
      end do
c
c  Apply permutation to dimension J.
c
      do j = 1, n
        a(1:n,p(2,j)) = b(1:n,j) 
      end do

      return
      end
      subroutine latin_cover_3d ( n, p, a )

c*********************************************************************72
c
cc LATIN_COVER_3D returns a 3D Latin Square Covering.
c
c  Discussion:
c
c    A basic solution is computed, and the user is permitted to permute
c    both I, J and K coordinates.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    26 June 2012
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the number of points.
c
c    Input, integer P(3,N), permutations to be applied
c    to the spatial dimensions.
c
c    Output, integer A(N,N,N), the Latin cover.  A(I,J,K) = L
c    means that (I,J,K) is one element of the L-th Latin square.
c
      implicit none

      integer n

      integer a(n,n,n)
      integer b(n,n,n)
      integer base
      integer i
      integer i4_wrap
      integer ik
      integer j
      integer jk
      integer k
      integer p(3,n)

      base = 1

      call perm_check ( n, p(1,1:n) )
      call perm_check ( n, p(2,1:n) )
      call perm_check ( n, p(3,1:n) )
c
c  Set up the basic solution.
c
      do i = 1, n
        do j = 1, n
          do k = 1, n
            ik = i4_wrap ( i + 1 - k, 1, n )
            jk = i4_wrap ( j + 1 - k, 1, n )
            b(i,j,k) = ik + ( jk - 1 ) * n
          end do
        end do
      end do
c
c  Apply permutation to dimension I.
c
      do i = 1, n
        a(p(1,i),1:n,1:n) = b(i,1:n,1:n) 
      end do
c
c  Apply permutation to dimension J.
c
      do j = 1, n
        b(1:n,p(2,j),1:n) = a(1:n,j,1:n) 
      end do
c
c  Apply permutation to dimension K.
c
      do k = 1, n
        a(1:n,1:n,p(3,k)) = b(1:n,1:n,k) 
      end do

      return
      end
      subroutine perm_check ( n, p )

c*********************************************************************72
c
cc PERM_CHECK checks that a vector represents a permutation.
c
c  Discussion:
c
c    The routine verifies that each of the integers from 1
c    to N occurs among the N entries of the permutation.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    25 January 2007
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the number of entries.
c
c    Input, integer P(N), the permutation, in standard index form.
c
      implicit none

      integer n

      integer ierror
      integer ifind
      integer iseek
      integer p(n)

      ierror = 0

      do iseek = 1, n

        ierror = iseek

        do ifind = 1, n
          if ( p(ifind) .eq. iseek ) then
            ierror = 0
            exit
          end if
        end do

        if ( ierror .ne. 0 ) then
          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) 'PERM_CHECK - Fatal error!'
          write ( *, '(a)' ) '  The input array does not represent'
          write ( *, '(a)' ) '  a proper permutation.  In particular,'
          write ( *, '(a,i8)' ) '  it is missing the value ', ierror
          stop
        end if

      end do

      return
      end
      subroutine perm_print ( n, p, title )

c*********************************************************************72
c
cc PERM_PRINT prints a permutation.
c
c  Example:
c
c    Input:
c
c      P = 7 2 4 1 5 3 6
c
c    Printed output:
c
c      "This is the permutation:"
c
c      1 2 3 4 5 6 7
c      7 2 4 1 5 3 6
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    25 January 2007
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the number of objects permuted.
c
c    Input, integer P(N), the permutation, in standard index form.
c
c    Input, character * ( * ) TITLE, an optional title.
c    If no title is supplied, then only the permutation is printed.
c
      implicit none

      integer n

      integer i
      integer ihi
      integer ilo
      integer inc
      parameter ( inc = 20 )
      integer p(n)
      character * ( * ) title
      integer title_length

      title_length = len_trim ( title )

      if ( 0 .lt. title_length ) then

        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) title(1:title_length)

        do ilo = 1, n, inc
          ihi = min ( n, ilo + inc - 1 )
          write ( *, '(a)' ) ' '
          write ( *, '(2x,20i4)' ) ( i, i = ilo, ihi )
          write ( *, '(2x,20i4)' ) ( p(i), i = ilo, ihi )
        end do

      else

        do ilo = 1, n, inc
          ihi = min ( n, ilo + inc - 1 )
          write ( *, '(2x,20i4)' ) ( p(i), i = ilo, ihi )
        end do

      end if

      return
      end
      subroutine perm_uniform ( n, base, seed, p )

c**********************************************************************72
c
cc PERM_UNIFORM selects a random permutation of N objects.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    26 June 2012
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Albert Nijenhuis, Herbert Wilf,
c    Combinatorial Algorithms,
c    Academic Press, 1978, second edition,
c    ISBN 0-12-519260-6.
c
c  Parameters:
c
c    Input, integer N, the number of objects to be permuted.
c
c    Input, integer BASE, is 0 for a 0-based permutation and 1 for 
c    a 1-based permutation.
c
c    Input/output, integer SEED, a seed for the random
c    number generator.
c
c    Output, integer P(N), the permutation.  P(I) is the "new"
c    location of the object originally at I.
c
      implicit none

      integer n

      integer base
      integer i
      integer i4_uniform
      integer j
      integer k
      integer p(n)
      integer seed

      do i = 1, n
        p(i) = ( i - 1 ) + base
      end do

      do i = 1, n
        j = i4_uniform ( i, n, seed )
        k    = p(i)
        p(i) = p(j)
        p(j) = k
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
