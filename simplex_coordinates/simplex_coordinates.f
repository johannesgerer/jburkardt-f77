      function r8_factorial ( n )

c*********************************************************************72
c
cc R8_FACTORIAL computes the factorial of N.
c
c  Discussion:
c
c    factorial ( N ) = product ( 1 <= I <= N ) I
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    07 June 2008
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the argument of the factorial function.
c    If N is less than 1, the function value is returned as 1.
c
c    Output, double precision R8_FACTORIAL, the factorial of N.
c
      implicit none

      integer i
      integer n
      double precision r8_factorial

      r8_factorial = 1.0D+00

      do i = 1, n
        r8_factorial = r8_factorial * dble ( i )
      end do

      return
      end
      subroutine r8mat_det ( n, a, det )

c*********************************************************************72
c
cc R8MAT_DET computes the determinant of an R8MAT.
c
c  Discussion:
c
c    An R8MAT is an array of R8 values.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    19 August 2010
c
c  Author:
c
c    Original FORTRAN77 version by Helmut Spaeth.
c    This FORTRAN77 version by John Burkardt.
c
c  Reference:
c
c    Helmut Spaeth,
c    Cluster Analysis Algorithms
c    for Data Reduction and Classification of Objects,
c    Ellis Horwood, 1980, page 125-127.
c
c  Parameters:
c
c    Input, integer N, the order of the matrix.
c
c    Input, double precision A(N,N), the matrix whose determinant is desired.
c
c    Output, double precision DET, the determinant of the matrix.
c
      implicit none

      integer n

      double precision a(n,n)
      double precision b(n,n)
      double precision det
      integer i
      integer j
      integer k
      integer m
      integer piv
      double precision t

      do j = 1, n
        do i = 1, n
          b(i,j) = a(i,j)
        end do
      end do

      det = 1.0D+00

      do k = 1, n

        piv = k
        do i = k + 1, n
          if ( abs ( b(piv,k) ) .lt. abs ( b(i,k) ) ) then
            piv = i
          end if
        end do

        m = piv

        if ( m .ne. k ) then
          det = - det
          t      = b(m,k)
          b(m,k) = b(k,k)
          b(k,k) = t
        end if

        det = det * b(k,k)

        if ( b(k,k) .ne. 0.0D+00 ) then

          do i = k + 1, n
            b(i,k) = - b(i,k) / b(k,k)
          end do

          do j = k + 1, n

            if ( m .ne. k ) then
              t      = b(m,j)
              b(m,j) = b(k,j)
              b(k,j) = t
            end if

            do i = k + 1, n
              b(i,j) = b(i,j) + b(i,k) * b(k,j)
            end do

          end do

        end if

      end do

      return
      end
      subroutine r8mat_transpose_print ( m, n, a, title )

c*********************************************************************72
c
cc R8MAT_TRANSPOSE_PRINT prints an R8MAT, transposed.
c
c  Discussion:
c
c    An R8MAT is an array of R8's.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    28 April 2008
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer M, N, the number of rows and columns.
c
c    Input, double precision A(M,N), an M by N matrix to be printed.
c
c    Input, character*(*) TITLE, a title.
c
      implicit none

      integer m
      integer n

      double precision a(m,n)
      character*(*) title

      call r8mat_transpose_print_some ( m, n, a, 1, 1, m, n, title )

      return
      end
      subroutine r8mat_transpose_print_some ( m, n, a, ilo, jlo, ihi, 
     &  jhi, title )

c*********************************************************************72
c
cc R8MAT_TRANSPOSE_PRINT_SOME prints some of an R8MAT transposed.
c
c  Discussion:
c
c    An R8MAT is an array of R8's.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    28 April 2008
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer M, N, the number of rows and columns.
c
c    Input, double precision A(M,N), an M by N matrix to be printed.
c
c    Input, integer ILO, JLO, the first row and column to print.
c
c    Input, integer IHI, JHI, the last row and column to print.
c
c    Input, character * ( * ) TITLE, a title.
c
      implicit none

      integer incx
      parameter ( incx = 5 )
      integer m
      integer n

      double precision a(m,n)
      character * ( 14 ) ctemp(incx)
      integer i
      integer i2
      integer i2hi
      integer i2lo
      integer ihi
      integer ilo
      integer inc
      integer j
      integer j2hi
      integer j2lo
      integer jhi
      integer jlo
      character * ( * ) title

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) trim ( title )

      do i2lo = max ( ilo, 1 ), min ( ihi, m ), incx

        i2hi = i2lo + incx - 1
        i2hi = min ( i2hi, m )
        i2hi = min ( i2hi, ihi )

        inc = i2hi + 1 - i2lo

        write ( *, '(a)' ) ' '

        do i = i2lo, i2hi
          i2 = i + 1 - i2lo
          write ( ctemp(i2), '(i8,6x)') i
        end do

        write ( *, '(''       Row'',5a14)' ) ctemp(1:inc)
        write ( *, '(a)' ) '       Col'

        j2lo = max ( jlo, 1 )
        j2hi = min ( jhi, n )

        do j = j2lo, j2hi

          do i2 = 1, inc
            i = i2lo - 1 + i2
            write ( ctemp(i2), '(g14.6)' ) a(i,j)
          end do

          write ( *, '(2x,i8,a,5a14)' ) j, ':', ( ctemp(i), i = 1, inc )

        end do

      end do

      return
      end
      function r8vec_dot_product ( n, v1, v2 )

c*********************************************************************72
c
cc R8VEC_DOT_PRODUCT finds the dot product of a pair of R8VEC's.
c
c  Discussion:
c
c    An R8VEC is a vector of R8 values.
c
c    In FORTRAN90, the system routine DOT_PRODUCT should be called
c    directly.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    27 May 2008
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the dimension of the vectors.
c
c    Input, double precision V1(N), V2(N), the vectors.
c
c    Output, double precision R8VEC_DOT_PRODUCT, the dot product.
c
      implicit none

      integer n

      integer i
      double precision r8vec_dot_product
      double precision v1(n)
      double precision v2(n)
      double precision value

      value = 0.0D+00
      do i = 1, n
        value = value + v1(i) * v2(i)
      end do

      r8vec_dot_product = value

      return
      end
      function r8vec_norm ( n, a )

c*********************************************************************72
c
cc R8VEC_NORM returns the L2 norm of an R8VEC.
c
c  Discussion:
c
c    An R8VEC is a vector of R8 values.
c
c    The vector L2 norm is defined as:
c
c      R8VEC_NORM = sqrt ( sum ( 1 <= I <= N ) A(I)**2 ).
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    27 May 2008
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the number of entries in A.
c
c    Input, double precision A(N), the vector whose L2 norm is desired.
c
c    Output, double precision R8VEC_NORM, the L2 norm of A.
c
      implicit none

      integer n

      double precision a(n)
      integer i
      double precision r8vec_norm
      double precision value

      value = 0.0D+00
      do i = 1, n
        value = value + a(i) * a(i)
      end do
      value = sqrt ( value )

      r8vec_norm = value

      return
      end
      function r8vec_sum ( n, v1 )

c*********************************************************************72
c
cc R8VEC_SUM sums the entries of an R8VEC.
c
c  Discussion:
c
c    An R8VEC is a vector of R8's.
c
c    In FORTRAN90, the system routine SUM should be called
c    directly.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    22 July 2008
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the dimension of the vectors.
c
c    Input, double precision V1(N), the vector.
c
c    Output, double precision R8VEC_SUM, the sum of the entries.
c
      implicit none

      integer n

      integer i
      double precision r8vec_sum
      double precision v1(n)
      double precision value

      value = 0.0D+00
      do i = 1, n
        value = value + v1(i)
      end do

      r8vec_sum = value

      return
      end
      subroutine simplex_coordinates1 ( n, x )

c*********************************************************************72
c
cc SIMPLEX_COORDINATES1 computes the Cartesian coordinates of simplex vertices.
c
c  Discussion:
c
c    The simplex will have its centroid at 0;
c
c    The sum of the vertices will be zero.
c
c    The distance of each vertex from the origin will be 1.
c
c    The length of each edge will be constant.
c
c    The dot product of the vectors defining any two vertices will be - 1 / N.
c    This also means the angle subtended by the vectors from the origin
c    to any two distinct vertices will be arccos ( - 1 / N ).
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    19 September 2010
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the spatial dimension.
c
c    Output, double precision X(N,N+1), the coordinates of the vertices
c    of a simplex in N dimensions.  
c
      implicit none

      integer n

      integer i
      integer ii
      integer j
      double precision r8vec_dot_product
      double precision s
      double precision x(n,n+1)

      do j = 1, n + 1
        do i = 1, n
          x(i,j) = 0.0D+00
        end do
      end do

      do i = 1, n
c
c  Set X(I,I) so that sum ( X(1:I,I)**2 ) = 1.
c
        s = 1.0D+00
        do ii = 1, i - 1
          s = s - x(ii,i)**2
        end do

        x(i,i) = sqrt ( s )
c
c  Set X(I,J) for J = I+1 to N+1 by using the fact that XI dot XJ = - 1 / N 
c
        do j = i + 1, n + 1
          x(i,j) = ( - 1.0D+00 / dble ( n ) 
     &      - r8vec_dot_product ( i - 1, x(1,i), x(1,j) ) ) / x(i,i)
        end do

      end do

      return
      end
      subroutine simplex_coordinates2 ( n, x )

c*********************************************************************72
c
cc SIMPLEX_COORDINATES2 computes the Cartesian coordinates of simplex vertices.
c
c  Discussion:
c
c    This routine uses a simple approach to determining the coordinates of
c    the vertices of a regular simplex in n dimensions.
c
c    We want the vertices of the simplex to satisfy the following conditions:
c
c    1) The centroid, or average of the vertices, is 0.
c    2) The distance of each vertex from the centroid is 1.
c       By 1), this is equivalent to requiring that the sum of the squares
c       of the coordinates of any vertex be 1.
c    3) The distance between any pair of vertices is equal (and is not zero.)
c    4) The dot product of any two coordinate vectors for distinct vertices
c       is -1/N; equivalently, the angle subtended by two distinct vertices
c       from the centroid is arccos ( -1/N).
c
c    Note that if we choose the first N vertices to be the columns of the
c    NxN identity matrix, we are almost there.  By symmetry, the last column
c    must have all entries equal to some value A.  Because the square of the
c    distance between the last column and any other column must be 2 (because
c    that's the distance between any pair of columns), we deduce that
c    (A-1)^2 + (N-1)*A^2 = 2, hence A = (1-sqrt(1+N))/N.  Now compute the 
c    centroid C of the vertices, and subtract that, to center the simplex 
c    around the origin.  Finally, compute the norm of one column, and rescale 
c    the matrix of coordinates so each vertex has unit distance from the origin.
c
c    This approach devised by John Burkardt, 19 September 2010.  What,
c    I'm not the first?
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    19 September 2010
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the spatial dimension.
c
c    Output, double precision X(N,N+1), the coordinates of the vertices
c    of a simplex in N dimensions.  
c
      implicit none

      integer n

      double precision a
      double precision c(n)
      integer i
      integer j
      double precision r8vec_norm
      double precision s
      double precision x(n,n+1)

      do j = 1, n + 1
        do i = 1, n
          x(i,j) = 0.0D+00
        end do
      end do

      do j = 1, n
        x(j,j) = 1.0D+00
      end do

      a = ( 1.0D+00 - sqrt ( 1.0D+00 + dble ( n ) ) ) / dble ( n )

      do i = 1, n
        x(i,n+1) = a
      end do
c
c  Now adjust coordinates so the centroid is at zero.
c
      do i = 1, n
        c(i) = 0.0D+00
        do j = 1, n + 1
          c(i) = c(i) + x(i,j)
        end do
        c(i) = c(i) / dble ( n + 1 )
      end do

      do j = 1, n + 1
        do i = 1, n
          x(i,j) = x(i,j) - c(i)
        end do
      end do
c
c  Now scale so each column has norm 1.
c
      s = r8vec_norm ( n, x(1,1) )

      do j = 1, n + 1
        do i = 1, n
          x(i,j) = x(i,j) / s
        end do
      end do

      return
      end
      subroutine simplex_volume ( n, x, volume )

c*********************************************************************72
c
cc SIMPLEX_VOLUME computes the volume of a simplex.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    19 September 2010
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the spatial dimension.
c
c    Input, double precision X(N,N+1), the coordinates of the vertices
c    of a simplex in N dimensions.  
c
c    Output, double precision VOLUME, the volume of the simplex.
c
      implicit none

      integer n

      double precision a(n,n)
      double precision det
      integer i
      integer j
      double precision volume
      double precision x(n,n+1)

      do j = 1, n
        do i = 1, n
          a(i,j) = x(i,j) - x(i,n+1)
        end do
      end do

      call r8mat_det ( n, a, det )

      volume = abs ( det )
      do i = 1, n
        volume = volume / dble ( i )
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
