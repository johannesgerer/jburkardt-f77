      subroutine c4mat_print_some ( m, n, a, ilo, jlo, ihi, jhi, 
     &  title )

c*********************************************************************72
c
cc C4MAT_PRINT_SOME prints some of a C4MAT.
c
c  Discussion:
c
c    A C4MAT is a matrix of C4's.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    13 December 2008
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer M, N, the number of rows and columns 
c    in the matrix.
c
c    Input, complex A(M,N), the matrix.
c
c    Input, integer ILO, JLO, IHI, JHI, the first row and
c    column, and the last row and column to be printed.
c
c    Input, character * ( * ) TITLE, a title.
c
      implicit none

      integer incx
      parameter ( incx = 4 )
      integer m
      integer n

      complex a(m,n)
      character * ( 20 ) ctemp(incx)
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
      character * ( * ) title
      complex zero

      zero = cmplx ( 0.0E+00, 0.0E+00 )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) trim ( title )
c
c  Print the columns of the matrix, in strips of INCX.
c
      do j2lo = jlo, min ( jhi, n ), incx

        j2hi = j2lo + incx - 1
        j2hi = min ( j2hi, n )
        j2hi = min ( j2hi, jhi )

        inc = j2hi + 1 - j2lo

        write ( *, '(a)' ) ' '

        do j = j2lo, j2hi
          j2 = j + 1 - j2lo
          write ( ctemp(j2), '(i10,10x)' ) j
        end do

        write ( *, '(a,4a20)' ) '  Col: ', ( ctemp(j2), j2 = 1, inc )
        write ( *, '(a)' ) '  Row'
        write ( *, '(a)' ) '  ---'
c
c  Determine the range of the rows in this strip.
c
        i2lo = max ( ilo, 1 )
        i2hi = min ( ihi, m )

        do i = i2lo, i2hi
c
c  Print out (up to) INCX entries in row I, that lie in the current strip.
c
          do j2 = 1, inc

            j = j2lo - 1 + j2

            if ( a(i,j) .eq. zero ) then
              ctemp(j2) = '       0.0          '
            else if ( aimag ( a(i,j) ) .eq. 0.0E+00 ) then
              write ( ctemp(j2), '(g10.3,10x)' ) real ( a(i,j) )
            else
              write ( ctemp(j2), '(2g10.3)' ) a(i,j)
            end if

          end do

          write ( *, '(i5,a1,4a20)' ) i, ':', ( ctemp(j2), j2 = 1, inc )

        end do

      end do

      return
      end
      subroutine c4mat_sftb ( n1, n2, y, x )

c*********************************************************************72
c
cc C4MAT_SFTB computes a "slow" backward Fourier transform of a C4MAT.
c
c  Discussion:
c
c    SFTF and SFTB are inverses of each other.  If we begin with data
c    X and apply SFTF to get Y, and then apply SFTB to Y,
c    we should get back the original X.
c
c    This routine is provided for illustration and testing.  It is inefficient
c    relative to optimized routines that use fast Fourier techniques.
c
c    For 0 <= I1 <= N1 - 1, 
c        0 <= I2 <= N2 - 1,
c
c      X(I1,I2) = Sum ( 0 <= K2 <= N2 - 1 ) Sum ( 0 <= K1 <= N1 - 1 ) 
c        Y(K1,K2) * exp ( 2 pi i I1 K1 / N1 ) * exp ( 2 pi i I2 K2 / N2 )
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    21 June 2010
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N1, N2, the number of rows and columns of data.
c
c    Input, complex Y(0:N1-1,0:N2-1), the Fourier coefficients.
c
c    Output, complex X(0:N1-1,0:N2-1), the data.
c
      implicit none

      integer n1
      integer n2

      complex cs1
      complex cs2
      integer i1
      integer i2
      integer j1
      integer j2
      real pi
      parameter ( pi = 3.141592653589793E+00 )
      real theta1
      real theta2
      complex x(0:n1-1,0:n2-1)
      complex y(0:n1-1,0:n2-1)

      do i2 = 0, n2 - 1
        do i1 = 0, n1 - 1
          x(i1,i2) = cmplx ( 0.0E+00, 0.0E+00 )
        end do
      end do

      do i2 = 0, n2 - 1
        do j2 = 0, n2 - 1
          theta2 = 2.0E+00 * pi * real ( i2 * j2 ) / real ( n2 )
          cs2 = cmplx ( cos ( theta2 ), - sin ( theta2 ) )
          do i1 = 0, n1 - 1
            do j1 = 0, n1 - 1
              theta1 = 2.0E+00 * pi * real ( i1 * j1 ) / real ( n1 )
              cs1 = cmplx ( cos ( theta1 ), - sin ( theta1 ) )
              x(i1,i2) = x(i1,i2) + y(j1,j2) * cs1 * cs2
            end do
          end do
        end do
      end do

      do i2 = 0, n2 - 1
        do i1 = 0, n1 - 1
          x(i1,i2) = x(i1,i2) / real ( n1 * n2 )
        end do
      end do

      return
      end
      subroutine c4mat_sftf ( n1, n2, x, y )

c*********************************************************************72
c
cc C4MAT_SFTF computes a "slow" forward Fourier transform of a C4MAT.
c
c  Discussion:
c
c    SFTF and SFTB are inverses of each other.  If we begin with data
c    X and apply SFTF to get Y, and then apply SFTB to Y, 
c    we should get back the original X.
c
c    This routine is provided for illustration and testing.  It is inefficient
c    relative to optimized routines that use fast Fourier techniques.
c
c    For 0 <= I1 <= N1 - 1, 
c        0 <= I2 <= N2 - 1,
c
c      Y(I1,I2) = Sum ( 0 <= K2 <= N2 - 1 ) Sum ( 0 <= K1 <= N1 - 1 ) 
c        X(K1,K2) * exp ( - 2 pi i I1 K1 / N1 ) * exp ( - 2 pi i I2 K2 / N2 )
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    21 June 2010
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N1, N2, the number of rows and columns of data.
c
c    Input, complex X(0:N1-1,0:N2-1), the data to be transformed.
c
c    Output, complex Y(0:N1-1,0:N2-1), the Fourier coefficients.
c
      implicit none

      integer n1
      integer n2

      complex cs1
      complex cs2
      integer i1
      integer i2
      integer j1
      integer j2
      real pi
      parameter ( pi = 3.141592653589793E+00 )
      real theta1
      real theta2
      complex x(0:n1-1,0:n2-1)
      complex y(0:n1-1,0:n2-1)

      do i2 = 0, n2 - 1
        do i1 = 0, n1 - 1
          y(i1,i2) = cmplx ( 0.0E+00, 0.0E+00, kind = 4 )
        end do
      end do

      do i2 = 0, n2 - 1
        do j2 = 0, n2 - 1
          theta2 = - 2.0E+00 * pi * real ( i2 * j2 ) / real ( n2 )
          cs2 = cmplx ( cos ( theta2 ), - sin ( theta2 ) )
          do i1 = 0, n1 - 1
            do j1 = 0, n1 - 1
              theta1 = - 2.0E+00 * pi * real ( i1 * j1 ) / real ( n1 )
              cs1 = cmplx ( cos ( theta1 ), - sin ( theta1 ) )
              y(i1,i2) = y(i1,i2) + x(j1,j2) * cs1 * cs2
            end do
          end do
        end do
      end do

      return
      end
      subroutine c4mat_uniform_01 ( m, n, seed, c )

c*********************************************************************72
c
cc C4MAT_UNIFORM_01 returns a unit pseudorandom C4MAT.
c
c  Discussion:
c
c    A C4MAT is a matrix of C4's.
c
c    The angles should be uniformly distributed between 0 and 2 * PI,
c    the square roots of the radius uniformly distributed between 0 and 1.
c
c    This results in a uniform distribution of values in the unit circle.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    13 December 2008
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer M, N, the number of rows and columns in the matrix.
c
c    Input/output, integer SEED, the "seed" value, which should NOT be 0.
c    On output, SEED has been updated.
c
c    Output, complex C(M,N), the pseudorandom complex matrix.
c
      implicit none

      integer m
      integer n

      complex c(m,n)
      integer i
      integer j
      integer k
      real r
      real r4_pi
      parameter ( r4_pi = 3.141592653589793E+00 )
      integer seed
      real theta

      if ( seed .eq. 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'C4MAT_UNIFORM_01 - Fatal errorc'
        write ( *, '(a)' ) '  Input value of SEED = 0.'
        stop
      end if

      do j = 1, n
        do i = 1, m

          k = seed / 127773

          seed = 16807 * ( seed - k * 127773 ) - k * 2836

          if ( seed .lt. 0 ) then
            seed = seed + 2147483647
          end if

          r = sqrt ( real ( seed ) * 4.656612875E-10 )

          k = seed / 127773

          seed = 16807 * ( seed - k * 127773 ) - k * 2836

          if ( seed .lt. 0 ) then
            seed = seed + 2147483647
          end if

          theta = 2.0E+00 * r4_pi * ( real ( seed ) * 4.656612875E-10 )

          c(i,j) = r * cmplx ( cos ( theta ), sin ( theta ) )

        end do

      end do

      return
      end
      subroutine c4vec_print_part ( n, a, max_print, title )

c*********************************************************************72
c
cc C4VEC_PRINT_PART prints "part" of a C4VEC.
c
c  Discussion:
c
c    The user specifies MAX_PRINT, the maximum number of lines to print.
c
c    If N, the size of the vector, is no more than MAX_PRINT, then
c    the entire vector is printed, one entry per line.
c
c    Otherwise, if possible, the first MAX_PRINT-2 entries are printed,
c    followed by a line of periods suggesting an omission,
c    and the last entry.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    22 June 2010
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the number of entries of the vector.
c
c    Input, complex A(N), the vector to be printed.
c
c    Input, integer MAX_PRINT, the maximum number of lines
c    to print.
c
c    Input, character * ( * ) TITLE, a title.
c
      implicit none

      integer n

      complex a(n)
      integer i
      integer max_print
      character * ( * )  title

      if ( max_print .le. 0 ) then
        return
      end if

      if ( n .le. 0 ) then
        return
      end if

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) trim ( title )
      write ( *, '(a)' ) ' '

      if ( n .le. max_print ) then

        do i = 1, n
          write ( *, '(2x,i8,a,1x,g14.6,2x,g14.6)' ) i, ':', a(i)
        end do

      else if ( 3 .le. max_print ) then

        do i = 1, max_print - 2
          write ( *, '(2x,i8,a,1x,g14.6,2x,g14.6)' ) i, ':', a(i)
        end do
        write ( *, '(a)' ) '  ........  ..............  ..............'
        i = n
        write ( *, '(2x,i8,a,1x,g14.6,2x,g14.6)' ) i, ':', a(i)

      else

        do i = 1, max_print - 1
          write ( *, '(2x,i8,a,1x,g14.6,2x,g14.6)' ) i, ':', a(i)
        end do
        i = max_print
        write ( *, '(2x,i8,a,1x,g14.6,2x,g14.6,2x,a)' ) i, ':', a(i), 
     &    '...more entries...'

      end if

      return
      end
      subroutine c4vec_sftb ( n, y, x )

c*********************************************************************72
c
cc C4VEC_SFTB computes a "slow" backward Fourier transform of a C4VEC.
c
c  Discussion:
c
c    SFTF and SFTB are inverses of each other.  If we begin with data
c    X and apply SFTF to get Y, and then apply SFTB to Y,
c    we should get back the original X.
c
c    This routine is provided for illustration and testing.  It is inefficient
c    relative to optimized routines that use fast Fourier techniques.
c
c    For 0 <= I <= N - 1
c
c      X(I) = 1/N * Sum ( 0 <= J <= N - 1 ) Y(J) * exp ( 2 pi i I J / N )
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    10 June 2010
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the number of data values.
c
c    Input, complex Y(0:N-1), the Fourier coefficients.
c
c    Output, complex X(0:N-1), the data.
c
      implicit none

      integer n

      integer i
      integer j
      real pi
      parameter ( pi = 3.141592653589793E+00 )
      real theta
      complex x(0:n-1)
      complex y(0:n-1)

      do i = 0, n - 1
        x(i) = cmplx ( 0.0E+00, 0.0E+00 )
        do j = 0, n - 1
          theta = - 2.0E+00 * pi * real ( i * j ) / real ( n )
          x(i) = x(i) + y(j) * cmplx ( cos ( theta ), sin ( theta ) )
        end do
        x(i) = x(i) / real ( n )
      end do

      return
      end
      subroutine c4vec_sftf ( n, x, y )

c*********************************************************************72
c
cc C4VEC_SFTF computes a "slow" forward Fourier transform of a C4VEC.
c
c  Discussion:
c
c    SFTF and SFTB are inverses of each other.  If we begin with data
c    X and apply SFTF to get Y, and then apply SFTB to Y, 
c    we should get back the original X.
c
c    This routine is provided for illustration and testing.  It is inefficient
c    relative to optimized routines that use fast Fourier techniques.
c
c    For 0 <= I <= N - 1
c
c      Y(I) = Sum ( 0 <= J <= N - 1 ) X(J) * exp ( - 2 pi i I J / N )
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    10 June 2010
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the number of data values.
c
c    Input, complex X(0:N-1), the data to be transformed.
c
c    Output, complex Y(0:N-1), the Fourier coefficients.
c
      implicit none

      integer n

      integer i
      integer j
      real pi
      parameter ( pi = 3.141592653589793E+00 )
      real theta
      complex x(0:n-1)
      complex y(0:n-1)

      do i = 0, n - 1
        y(i) = cmplx ( 0.0E+00, 0.0E+00 )
        do j = 0, n - 1
          theta = - 2.0E+00 * pi * real ( i * j ) / real ( n )
          y(i) = y(i) + x(j) * cmplx ( cos ( theta ), - sin ( theta ) )
        end do
      end do

      return
      end
      subroutine c4vec_uniform_01 ( n, seed, c )

c*********************************************************************72
c
cc C4VEC_UNIFORM_01 returns a unit pseudorandom C4VEC.
c
c  Discussion:
c
c    A C4VEC is a vector of complex single precision values.
c
c    The angles should be uniformly distributed between 0 and 2 * PI,
c    the square roots of the radius uniformly distributed between 0 and 1.
c
c    This results in a uniform distribution of values in the unit circle.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    06 January 2006
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Paul Bratley, Bennett Fox, Linus Schrage,
c    A Guide to Simulation,
c    Second Edition,
c    Springer, 1987,
c    ISBN: 0387964673,
c    LC: QA76.9.C65.B73.
c
c    Bennett Fox,
c    Algorithm 647:
c    Implementation and Relative Efficiency of Quasirandom
c    Sequence Generators,
c    ACM Transactions on Mathematical Software,
c    Volume 12, Number 4, December 1986, pages 362-376.
c
c    Pierre L'Ecuyer,
c    Random Number Generation,
c    in Handbook of Simulation,
c    edited by Jerry Banks,
c    Wiley, 1998,
c    ISBN: 0471134031,
c    LC: T57.62.H37.
c
c    Peter Lewis, Allen Goodman, James Miller,
c    A Pseudo-Random Number Generator for the System/360,
c    IBM Systems Journal,
c    Volume 8, Number 2, 1969, pages 136-143.
c
c  Parameters:
c
c    Input, integer N, the number of values to compute.
c
c    Input/output, integer SEED, the "seed" value, which should NOT be 0.
c    On output, SEED has been updated.
c
c    Output, complex C(N), the pseudorandom complex vector.
c
      implicit none

      integer n

      complex c(n)
      integer i
      integer i4_huge
      parameter ( i4_huge = 2147483647 )
      real r
      integer k
      real pi
      parameter ( pi = 3.141526E+00 )
      integer seed
      real theta

      if ( seed .eq. 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'C4VEC_UNIFORM_01 - Fatal error!'
        write ( *, '(a)' ) '  Input value of SEED = 0.'
        stop
      end if

      do i = 1, n

        k = seed / 127773

        seed = 16807 * ( seed - k * 127773 ) - k * 2836

        if ( seed .lt. 0 ) then
          seed = seed + i4_huge
        end if

        r = sqrt ( real ( dble ( seed ) * 4.656612875D-10 ) )

        k = seed / 127773

        seed = 16807 * ( seed - k * 127773 ) - k * 2836

        if ( seed .lt. 0 ) then
          seed = seed + i4_huge
        end if

        theta = 2.0E+00 * pi * real ( dble ( seed ) * 4.656612875D-10 )

        c(i) = r * cmplx ( cos ( theta ), sin ( theta ) )

      end do

      return
      end
      subroutine c8mat_print_some ( m, n, a, ilo, jlo, ihi, jhi, 
     &  title )

c*********************************************************************72
c
cc C8MAT_PRINT_SOME prints some of a C8MAT.
c
c  Discussion:
c
c    A C8MAT is a matrix of C8's.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    21 June 2010
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer M, N, the number of rows and columns 
c    in the matrix.
c
c    Input, double complex A(M,N), the matrix.
c
c    Input, integer ILO, JLO, IHI, JHI, the first row and
c    column, and the last row and column to be printed.
c
c    Input, character * ( * ) TITLE, a title.
c
      implicit none

      integer incx
      parameter ( incx = 4 )
      integer m
      integer n

      double complex a(m,n)
      character * ( 20 ) ctemp(incx)
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
      character * ( * ) title
      double complex zero

      zero = dcmplx ( 0.0D+00, 0.0D+00 )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) trim ( title )
c
c  Print the columns of the matrix, in strips of INCX.
c
      do j2lo = jlo, min ( jhi, n ), incx

        j2hi = j2lo + incx - 1
        j2hi = min ( j2hi, n )
        j2hi = min ( j2hi, jhi )

        inc = j2hi + 1 - j2lo

        write ( *, '(a)' ) ' '

        do j = j2lo, j2hi
          j2 = j + 1 - j2lo
          write ( ctemp(j2), '(i10,10x)' ) j
        end do

        write ( *, '(a,4a20)' ) '  Col: ', ( ctemp(j2), j2 = 1, inc )
        write ( *, '(a)' ) '  Row'
        write ( *, '(a)' ) '  ---'
c
c  Determine the range of the rows in this strip.
c
        i2lo = max ( ilo, 1 )
        i2hi = min ( ihi, m )

        do i = i2lo, i2hi
c
c  Print out (up to) INCX entries in row I, that lie in the current strip.
c
          do j2 = 1, inc

            j = j2lo - 1 + j2

            if ( a(i,j) .eq. zero ) then
              ctemp(j2) = '       0.0          '
            else if ( dimag ( a(i,j) ) .eq. 0.0D+00 ) then
              write ( ctemp(j2), '(g10.3,10x)' ) dreal ( a(i,j) )
            else
              write ( ctemp(j2), '(2g10.3)' ) a(i,j)
            end if

          end do

          write ( *, '(i5,a,4a20)' ) i, ':', ( ctemp(j2), j2 = 1, inc )

        end do

      end do

      return
      end
      subroutine c8mat_sftb ( n1, n2, y, x )

c*********************************************************************72
c
cc C8MAT_SFTB computes a "slow" backward Fourier transform of a C8MAT.
c
c  Discussion:
c
c    SFTF and SFTB are inverses of each other.  If we begin with data
c    X and apply SFTF to get Y, and then apply SFTB to Y,
c    we should get back the original X.
c
c    This routine is provided for illustration and testing.  It is inefficient
c    relative to optimized routines that use fast Fourier techniques.
c
c    For 0 <= I1 <= N1 - 1, 
c        0 <= I2 <= N2 - 1,
c
c      X(I1,I2) = Sum ( 0 <= K2 <= N2 - 1 ) Sum ( 0 <= K1 <= N1 - 1 ) 
c        Y(K1,K2) * exp ( 2 pi i I1 K1 / N1 ) * exp ( 2 pi i I2 K2 / N2 )
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    21 June 2010
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N1, N2, the number of rows and columns of data.
c
c    Input, double complex Y(0:N1-1,0:N2-1), the Fourier coefficients.
c
c    Output, double complex X(0:N1-1,0:N2-1), the data.
c
      implicit none

      integer n1
      integer n2

      double complex cs1
      double complex cs2
      integer i1
      integer i2
      integer j1
      integer j2
      double precision pi
      parameter ( pi = 3.141592653589793D+00 )
      double precision theta1
      double precision theta2
      double complex x(0:n1-1,0:n2-1)
      double complex y(0:n1-1,0:n2-1)

      do i2 = 0, n2 - 1
        do i1 = 0, n1 - 1
          x(i1,i2) = cmplx ( 0.0D+00, 0.0D+00 )
        end do
      end do

      do i2 = 0, n2 - 1
        do j2 = 0, n2 - 1
          theta2 = 2.0D+00 * pi * dble ( i2 * j2 ) / dble ( n2 )
          cs2 = dcmplx ( dcos ( theta2 ), - dsin ( theta2 ) )
          do i1 = 0, n1 - 1
            do j1 = 0, n1 - 1
              theta1 = 2.0D+00 * pi * dble ( i1 * j1 ) / dble ( n1 )
              cs1 = dcmplx ( dcos ( theta1 ), - dsin ( theta1 ) )
              x(i1,i2) = x(i1,i2) + y(j1,j2) * cs1 * cs2
            end do
          end do
        end do
      end do

      do i2 = 0, n2 - 1
        do i1 = 0, n1 - 1
          x(i1,i2) = x(i1,i2) / dble ( n1 * n2 )
        end do
      end do

      return
      end
      subroutine c8mat_sftf ( n1, n2, x, y )

c*********************************************************************72
c
cc C8MAT_SFTF computes a "slow" forward Fourier transform of a C8MAT.
c
c  Discussion:
c
c    SFTF and SFTB are inverses of each other.  If we begin with data
c    X and apply SFTF to get Y, and then apply SFTB to Y, 
c    we should get back the original X.
c
c    This routine is provided for illustration and testing.  It is inefficient
c    relative to optimized routines that use fast Fourier techniques.
c
c    For 0 <= I1 <= N1 - 1, 
c        0 <= I2 <= N2 - 1,
c
c      Y(I1,I2) = Sum ( 0 <= K2 <= N2 - 1 ) Sum ( 0 <= K1 <= N1 - 1 ) 
c        X(K1,K2) * exp ( - 2 pi i I1 K1 / N1 ) * exp ( - 2 pi i I2 K2 / N2 )
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    21 June 2010
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N1, N2, the number of rows and columns of data.
c
c    Input, double complex X(0:N1-1,0:N2-1), the data to be transformed.
c
c    Output, double complex Y(0:N1-1,0:N2-1), the Fourier coefficients.
c
      implicit none

      integer n1
      integer n2

      double complex cs1
      double complex cs2
      integer i1
      integer i2
      integer j1
      integer j2
      double precision pi
      parameter ( pi = 3.141592653589793D+00 )
      double precision theta1
      double precision theta2
      double complex x(0:n1-1,0:n2-1)
      double complex y(0:n1-1,0:n2-1)

      do i2 = 0, n2 - 1
        do i1 = 0, n1 - 1
          y(i1,i2) = cmplx ( 0.0D+00, 0.0D+00, kind = 4 )
        end do
      end do

      do i2 = 0, n2 - 1
        do j2 = 0, n2 - 1
          theta2 = - 2.0D+00 * pi * dble ( i2 * j2 ) / dble ( n2 )
          cs2 = dcmplx ( dcos ( theta2 ), - dsin ( theta2 ) )
          do i1 = 0, n1 - 1
            do j1 = 0, n1 - 1
              theta1 = - 2.0D+00 * pi * dble ( i1 * j1 ) / dble ( n1 )
              cs1 = dcmplx ( dcos ( theta1 ), - dsin ( theta1 ) )
              y(i1,i2) = y(i1,i2) + x(j1,j2) * cs1 * cs2
            end do
          end do
        end do
      end do

      return
      end
      subroutine c8mat_uniform_01 ( m, n, seed, c )

c*********************************************************************72
c
cc C8MAT_UNIFORM_01 returns a unit pseudorandom C8MAT.
c
c  Discussion:
c
c    A C8MAT is a matrix of complex double precision values.
c
c    The angles should be uniformly distributed between 0 and 2 * PI,
c    the square roots of the radius uniformly distributed between 0 and 1.
c
c    This results in a uniform distribution of values in the unit circle.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    05 March 2006
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Paul Bratley, Bennett Fox, Linus Schrage,
c    A Guide to Simulation,
c    Second Edition,
c    Springer, 1987,
c    ISBN: 0387964673,
c    LC: QA76.9.C65.B73.
c
c    Bennett Fox,
c    Algorithm 647:
c    Implementation and Relative Efficiency of Quasirandom
c    Sequence Generators,
c    ACM Transactions on Mathematical Software,
c    Volume 12, Number 4, December 1986, pages 362-376.
c
c    Pierre L'Ecuyer,
c    Random Number Generation,
c    in Handbook of Simulation,
c    edited by Jerry Banks,
c    Wiley, 1998,
c    ISBN: 0471134031,
c    LC: T57.62.H37.
c
c    Peter Lewis, Allen Goodman, James Miller,
c    A Pseudo-Random Number Generator for the System/360,
c    IBM Systems Journal,
c    Volume 8, Number 2, 1969, pages 136-143.
c
c  Parameters:
c
c    Input, integer M, N, the number of rows and columns in the matrix.
c
c    Input/output, integer SEED, the "seed" value, which should NOT be 0.
c    On output, SEED has been updated.
c
c    Output, double complex C(M,N), the pseudorandom complex matrix.
c
      implicit none

      integer m
      integer n

      double complex c(m,n)
      integer i
      integer i4_huge
      parameter ( i4_huge = 2147483647 )
      integer j
      double precision r
      integer k
      double precision pi
      parameter ( pi = 3.141592653589793D+00 )
      integer seed
      double precision theta

      if ( seed .eq. 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'C8MAT_UNIFORM_01 - Fatal error!'
        write ( *, '(a)' ) '  Input value of SEED = 0.'
        stop
      end if

      do j = 1, n
        do i = 1, m

          k = seed / 127773

          seed = 16807 * ( seed - k * 127773 ) - k * 2836

          if ( seed .lt. 0 ) then
            seed = seed + i4_huge
          end if

          r = sqrt ( dble ( seed ) * 4.656612875D-10 )

          k = seed / 127773

          seed = 16807 * ( seed - k * 127773 ) - k * 2836

          if ( seed .lt. 0 ) then
            seed = seed + i4_huge
          end if

          theta = 2.0D+00 * pi * ( dble ( seed ) * 4.656612875D-10 )

          c(i,j) = r * dcmplx ( dcos ( theta ), dsin ( theta ) )

        end do

      end do

      return
      end
      subroutine c8vec_print_part ( n, a, max_print, title )

c*********************************************************************72
c
cc C8VEC_PRINT_PART prints "part" of a C8VEC.
c
c  Discussion:
c
c    The user specifies MAX_PRINT, the maximum number of lines to print.
c
c    If N, the size of the vector, is no more than MAX_PRINT, then
c    the entire vector is printed, one entry per line.
c
c    Otherwise, if possible, the first MAX_PRINT-2 entries are printed,
c    followed by a line of periods suggesting an omission,
c    and the last entry.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    22 June 2010
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the number of entries of the vector.
c
c    Input, double complex A(N), the vector to be printed.
c
c    Input, integer MAX_PRINT, the maximum number of lines
c    to print.
c
c    Input, character * ( * ) TITLE, a title.
c
      implicit none

      integer n

      double complex a(n)
      integer i
      integer max_print
      character * ( * )  title

      if ( max_print .le. 0 ) then
        return
      end if

      if ( n .le. 0 ) then
        return
      end if

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) trim ( title )
      write ( *, '(a)' ) ' '

      if ( n .le. max_print ) then

        do i = 1, n
          write ( *, '(2x,i8,a,1x,g14.6,2x,g14.6)' ) i, ':', a(i)
        end do

      else if ( 3 .le. max_print ) then

        do i = 1, max_print - 2
          write ( *, '(2x,i8,a,1x,g14.6,2x,g14.6)' ) i, ':', a(i)
        end do
        write ( *, '(a)' ) '  ........  ..............  ..............'
        i = n
        write ( *, '(2x,i8,a,1x,g14.6,2x,g14.6)' ) i, ':', a(i)

      else

        do i = 1, max_print - 1
          write ( *, '(2x,i8,a,1x,g14.6,2x,g14.6)' ) i, ':', a(i)
        end do
        i = max_print
        write ( *, '(2x,i8,a,1x,g14.6,2x,g14.6,2x,a)' ) i, ':', a(i), 
     &    '...more entries...'

      end if

      return
      end
      subroutine c8vec_sftb ( n, y, x )

c*********************************************************************72
c
cc C8VEC_SFTB computes a "slow" backward Fourier transform of a C8VEC.
c
c  Discussion:
c
c    SFTF and SFTB are inverses of each other.  If we begin with data
c    X and apply SFTF to get Y, and then apply SFTB to Y,
c    we should get back the original X.
c
c    This routine is provided for illustration and testing.  It is inefficient
c    relative to optimized routines that use fast Fourier techniques.
c
c    For 0 <= I <= N - 1
c
c      X(I) = 1/N * Sum ( 0 <= J <= N - 1 ) Y(J) * exp ( 2 pi i I J / N )
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    10 June 2010
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the number of data values.
c
c    Input, double complex Y(0:N-1), the Fourier coefficients.
c
c    Output, double complex X(0:N-1), the data.
c
      implicit none

      integer n

      integer i
      integer j
      double precision pi
      parameter ( pi = 3.141592653589793D+00 )
      double precision theta
      double complex x(0:n-1)
      double complex y(0:n-1)

      do i = 0, n - 1
        x(i) = dcmplx ( 0.0D+00, 0.0D+00 )
        do j = 0, n - 1
          theta = - 2.0D+00 * pi * dble ( i * j ) / dble ( n )
          x(i) = x(i) + y(j) * dcmplx ( dcos ( theta ), dsin ( theta ) )
        end do
        x(i) = x(i) / dble ( n )
      end do

      return
      end
      subroutine c8vec_sftf ( n, x, y )

c*********************************************************************72
c
cc C8VEC_SFTF computes a "slow" forward Fourier transform of a C8VEC.
c
c  Discussion:
c
c    SFTF and SFTB are inverses of each other.  If we begin with data
c    X and apply SFTF to get Y, and then apply SFTB to Y, 
c    we should get back the original X.
c
c    This routine is provided for illustration and testing.  It is inefficient
c    relative to optimized routines that use fast Fourier techniques.
c
c    For 0 <= I <= N - 1
c
c      Y(I) = Sum ( 0 <= J <= N - 1 ) X(J) * exp ( - 2 pi i I J / N )
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    10 June 2010
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the number of data values.
c
c    Input, complex X(0:N-1), the data to be transformed.
c
c    Output, complex Y(0:N-1), the Fourier coefficients.
c
      implicit none

      integer n

      integer i
      integer j
      double precision pi
      parameter ( pi = 3.141592653589793D+00 )
      double precision theta
      double complex x(0:n-1)
      double complex y(0:n-1)

      do i = 0, n - 1
        y(i) = dcmplx ( 0.0D+00, 0.0D+00 )
        do j = 0, n - 1
          theta = - 2.0D+00 * pi * dble ( i * j ) / dble ( n )
          y(i) = y(i) + x(j) 
     &      * dcmplx ( dcos ( theta ), - dsin ( theta ) )
        end do
      end do

      return
      end
      subroutine c8vec_uniform_01 ( n, seed, c )

c*********************************************************************72
c
cc C8VEC_UNIFORM_01 returns a unit pseudorandom C8VEC.
c
c  Discussion:
c
c    A C8VEC is a vector of complex double precision values.
c
c    The angles should be uniformly distributed between 0 and 2 * PI,
c    the square roots of the radius uniformly distributed between 0 and 1.
c
c    This results in a uniform distribution of values in the unit circle.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    05 March 2006
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Paul Bratley, Bennett Fox, Linus Schrage,
c    A Guide to Simulation,
c    Second Edition,
c    Springer, 1987,
c    ISBN: 0387964673,
c    LC: QA76.9.C65.B73.
c
c    Bennett Fox,
c    Algorithm 647:
c    Implementation and Relative Efficiency of Quasirandom
c    Sequence Generators,
c    ACM Transactions on Mathematical Software,
c    Volume 12, Number 4, December 1986, pages 362-376.
c
c    Pierre L'Ecuyer,
c    Random Number Generation,
c    in Handbook of Simulation,
c    edited by Jerry Banks,
c    Wiley, 1998,
c    ISBN: 0471134031,
c    LC: T57.62.H37.
c
c    Peter Lewis, Allen Goodman, James Miller,
c    A Pseudo-Random Number Generator for the System/360,
c    IBM Systems Journal,
c    Volume 8, Number 2, 1969, pages 136-143.
c
c  Parameters:
c
c    Input, integer N, the number of values to compute.
c
c    Input/output, integer SEED, the "seed" value, which should NOT be 0.
c    On output, SEED has been updated.
c
c    Output, double complex C(N), the pseudorandom complex vector.
c
      implicit none

      integer n

      double complex c(n)
      integer i
      integer i4_huge
      parameter ( i4_huge = 2147483647 )
      integer k
      double precision pi
      parameter ( pi = 3.141592653589793D+00 )
      double precision r
      integer seed
      double precision theta

      if ( seed .eq. 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'C8VEC_UNIFORM_01 - Fatal error!'
        write ( *, '(a)' ) '  Input value of SEED = 0.'
        stop
      end if

      do i = 1, n

        k = seed / 127773

        seed = 16807 * ( seed - k * 127773 ) - k * 2836

        if ( seed .lt. 0 ) then
          seed = seed + i4_huge
        end if

        r = sqrt ( dble ( seed ) * 4.656612875D-10 )

        k = seed / 127773

        seed = 16807 * ( seed - k * 127773 ) - k * 2836

        if ( seed .lt. 0 ) then
          seed = seed + i4_huge
        end if

        theta = 2.0D+00 * pi * ( dble ( seed ) * 4.656612875D-10 )

        c(i) = r * dcmplx ( dcos ( theta ), dsin ( theta ) )

      end do

      return
      end
      subroutine r4vec_print_part ( n, a, max_print, title )

c*********************************************************************72
c
cc R4VEC_PRINT_PART prints "part" of an R4VEC.
c
c  Discussion:
c
c    The user specifies MAX_PRINT, the maximum number of lines to print.
c
c    If N, the size of the vector, is no more than MAX_PRINT, then
c    the entire vector is printed, one entry per line.
c
c    Otherwise, if possible, the first MAX_PRINT-2 entries are printed,
c    followed by a line of periods suggesting an omission,
c    and the last entry.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    22 June 2010
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the number of entries of the vector.
c
c    Input, real A(N), the vector to be printed.
c
c    Input, integer MAX_PRINT, the maximum number of lines to print.
c
c    Input, character*(*) TITLE, a title.
c
      implicit none

      integer n

      real a(n)
      integer i
      integer max_print
      character*(*) title

      if ( max_print .le. 0 ) then
        return
      end if

      if ( n .le. 0 ) then
        return
      end if

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) trim ( title )
      write ( *, '(a)' ) ' '

      if ( n .le. max_print ) then

        do i = 1, n
          write ( *, '(2x,i8,a,1x,g14.6)' ) i, ':', a(i)
        end do

      else if ( 3 .le. max_print ) then

        do i = 1, max_print - 2
          write ( *, '(2x,i8,a,1x,g14.6)' ) i, ':', a(i)
        end do

        write ( *, '(a)' ) '  ........  ..............'
        i = n

        write ( *, '(2x,i8,a,1x,g14.6)' ) i, ':', a(i)

      else

        do i = 1, max_print - 1
          write ( *, '(2x,i8,a,1x,g14.6)' ) i, ':', a(i)
        end do

        i = max_print

        write ( *, '(2x,i8,a,1x,g14.6,a)' ) 
     &    i, ':', a(i), '...more entries...'

      end if

      return
      end
      subroutine r4vec_sct ( n, x, y )

c*********************************************************************72
c
cc R4VEC_SCT computes a "slow" cosine transform of an R4VEC.
c
c  Discussion:
c
c    This routine is provided for illustration and testing.  It is inefficient
c    relative to optimized routines that use fast Fourier techniques.
c
c      Y(1) = Sum ( 1 <= J <= N ) X(J) 
c
c      For 2 <= I <= N-1:
c
c        Y(I) = 2 * Sum ( 1 <= J <= N ) X(J) 
c          * cos ( PI * ( I - 1 ) * ( J - 1 ) / ( N - 1 ) )
c
c      Y(N) = Sum ( X(1:N:2) ) - Sum ( X(2:N:2) )
c
c    Applying the routine twice in succession should yield the original data,
c    multiplied by 2 * ( N + 1 ).  This is a good check for correctness 
c    and accuracy.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    10 June 2010
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the number of data values.
c
c    Input, real X(N), the data sequence.
c
c    Output, real Y(N), the transformed data.
c
      implicit none

      integer n

      real angle
      integer i
      integer j
      real pi
      parameter ( pi = 3.141592653589793E+00 )
      real x(n)
      real y(n)

      do i = 1, n

        y(i) = x(1) / 2.0E+00

        do j = 2, n - 1
          angle = pi 
     &      * real ( mod ( ( i - 1 ) * ( j - 1 ), 2 * ( n - 1 ) ) ) 
     &      / real ( n - 1 )
          y(i) = y(i) + x(j) * cos ( angle )
        end do

        j = n

        angle = pi 
     &    * real ( mod ( ( i - 1 ) * ( j - 1 ), 2 * ( n - 1 ) ) ) 
     &    / real ( n - 1 )

        y(i) = y(i) + x(n) * cos ( angle ) / 2.0E+00

      end do

      do i = 1, n
        y(i) = 2.0E+00 * y(i) * sqrt ( real ( n ) / real ( n - 1 ) )
      end do

      return
      end
      subroutine r4vec_sftb ( n, azero, a, b, r )

c*********************************************************************72
c
cc R4VEC_SFTB computes a "slow" backward Fourier transform of an R4VEC.
c
c  Discussion:
c
c    SFTB and SFTF are inverses of each other.  If we begin with data
c    AZERO, A, and B, and apply SFTB to it, and then apply SFTF to the
c    resulting R vector, we should get back the original AZERO, A and B.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    10 June 2010
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the number of data values. 
c
c    Input, real AZERO, the constant Fourier coefficient.
c
c    Input, real A(N/2), B(N/2), the Fourier coefficients.
c
c    Output, real R(N), the reconstructed data sequence.
c
      implicit none

      integer n

      real a(n/2)
      real azero
      real b(n/2)
      integer i
      integer k
      real r(n)
      real pi
      parameter ( pi = 3.141592653589793E+00 )
      real theta

      do i = 1, n
        r(i) = azero
      end do

      do i = 1, n
        do k = 1, n / 2
          theta = real ( k * ( i - 1 ) * 2 ) * pi / real ( n )
          r(i) = r(i) + a(k) * cos ( theta ) + b(k) * sin ( theta )
        end do
      end do

      return
      end
      subroutine r4vec_sftf ( n, r, azero, a, b )

c*********************************************************************72
c
cc R4VEC_SFTF computes a "slow" forward Fourier transform of an R4VEC.
c
c  Discussion:
c
c    SFTF and SFTB are inverses of each other.  If we begin with data
c    R and apply SFTB to it, and then apply SFTB to the resulting AZERO, 
c    A, and B, we should get back the original R.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    10 June 2010
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the number of data values.
c
c    Input, real R(N), the data to be transformed.
c
c    Output, real AZERO, = sum ( 1 <= I <= N ) R(I) / N.
c
c    Output, real A(N/2), B(N/2), the Fourier coefficients.
c
      implicit none

      integer n

      real a(1:n/2)
      real azero
      real b(1:n/2)
      integer i
      integer j
      real pi
      parameter ( pi = 3.141592653589793E+00 )
      real r(n)
      real theta

      azero = 0.0E+00
      do i = 1, n
        azero = azero + r(i)
      end do
      azero = azero / real ( n )

      do i = 1, n / 2

        a(i) = 0.0E+00
        b(i) = 0.0E+00

        do j = 1, n
          theta = real ( 2 * i * ( j - 1 ) ) * pi 
     &      / real ( n )
          a(i) = a(i) + r(j) * cos ( theta )
          b(i) = b(i) + r(j) * sin ( theta )
        end do

        a(i) = a(i) / real ( n )
        b(i) = b(i) / real ( n )

        if ( i .ne. ( n / 2 ) ) then
          a(i) = 2.0E+00 * a(i)
          b(i) = 2.0E+00 * b(i)
        end if

      end do

      return
      end
      subroutine r4vec_sht ( n, a, b  )

c*********************************************************************72
c
cc R4VEC_SHT computes a "slow" Hartley transform of an R4VEC.
c
c  Discussion:
c
c    The discrete Hartley transform B of a set of data A is 
c
c      B(I) = 1/sqrt(N) * Sum ( 0 <= J <= N-1 ) A(J) * CAS(2*PI*I*J/N)
c
c    Here, the data and coefficients are indexed from 0 to N-1.
c
c    With the above normalization factor of 1/sqrt(N), the Hartley
c    transform is its own inverse.
c
c    This routine is provided for illustration and testing.  It is inefficient
c    relative to optimized routines.  
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    10 June 2010
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Ralph Hartley,
c    A More Symmetrical Fourier Analysis Applied to Transmission Problems,
c    Proceedings of the Institute of Radio Engineers,
c    Volume 30, pages 144-150, 1942.
c
c  Parameters:
c
c    Input, integer N, the number of data values.
c
c    Input, real A(0:N-1), the data to be transformed.
c
c    Output, real B(0:N-1), the transformed data.
c
      implicit none

      integer n

      real a(0:n-1)
      real b(0:n-1)
      integer i
      integer j
      real pi
      parameter ( pi = 3.141592653589793E+00 )
      real theta

      do i = 0, n - 1
        b(i) = 0.0E+00
      end do

      do i = 0, n - 1
        do j = 0, n - 1
          theta = 2.0E+00 * pi * real ( mod ( i * j, n ) ) 
     &      / real ( n )
          b(i) = b(i) + a(j) * ( cos ( theta ) + sin ( theta ) )
        end do
      end do

      do i = 0, n - 1
        b(i) = b(i) / sqrt ( real ( n ) )
      end do

      return
      end
      subroutine r4vec_sqctb ( n, x, y )

c*********************************************************************72
c
cc R4VEC_SQCTB computes a "slow" quarter cosine transform backward of an R4VEC.
c
c  Discussion:
c
c    This routine is provided for illustration and testing.  It is inefficient
c    relative to optimized routines that use fast Fourier techniques.
c
c    For 0 <= I <= N-1,
c
c      Y(I) = X(0) + 2 Sum ( 1 <= J <= N-1 ) X(J) * cos ( PI * J * (I+1/2) / N )
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    10 June 2010
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    William Briggs, Van Emden Henson,
c    The Discrete Fourier Transform, 
c    SIAM, 1995,
c    LC: QA403.5 B75
c
c  Parameters:
c
c    Input, integer N, the number of data values.
c
c    Input, real X(0:N-1), the data sequence.
c
c    Output, real Y(0:N-1), the transformed data.
c
      implicit none

      integer n

      integer i
      integer j
      real pi
      parameter ( pi = 3.141592653589793E+00 )
      real theta
      real x(0:n-1)
      real y(0:n-1)

      do i = 0, n - 1
        y(i) = x(0)
      end do

      do i = 0, n - 1
        do j = 1, n - 1

          theta = 0.5E+00 * pi * real ( j * ( 2 * i + 1 ) )
     &      / real ( n )
          y(i) = y(i) + 2.0E+00 * x(j) * cos ( theta  )

        end do

      end do

      return
      end
      subroutine r4vec_sqctf ( n, x, y )

c*********************************************************************72
c
cc R4VEC_SQCTF computes a "slow" quarter cosine transform forward of an R4VEC.
c
c  Discussion:
c
c    This routine is provided for illustration and testing.  It is inefficient
c    relative to optimized routines that use fast Fourier techniques.
c
c    For 0 <= I <= N-1,
c
c      Y(I) = (1/N) Sum ( 0 <= J <= N-1 ) X(J) * cos ( PI * I * (J+1/2) / N )
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    10 June 2010
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    William Briggs, Van Emden Henson,
c    The Discrete Fourier Transform, 
c    SIAM, 1995,
c    QA403.5 B75
c
c  Parameters:
c
c    Input, integer N, the number of data values.
c
c    Input, real X(0:N-1), the data sequence.
c
c    Output, real Y(0:N-1), the transformed data.
c
      implicit none

      integer n

      integer i
      integer j
      real pi
      parameter ( pi = 3.141592653589793E+00 )
      real theta
      real x(0:n-1)
      real y(0:n-1)

      do i = 0, n - 1
        y(i) = 0.0E+00
      end do

      do i = 0, n - 1
        do j = 0, n - 1
          theta = 0.5E+00 * pi * real ( i * ( 2 * j + 1 ) ) 
     &      / real ( n )
          y(i) = y(i) + x(j) * cos ( theta  )
        end do
      end do

      do i = 0, n - 1
        y(i) = y(i) / real ( n )
      end do

      return
      end
      subroutine r4vec_sqstb ( n, x, y )

c*********************************************************************72
c
cc R4VEC_SQSTB computes a "slow" quarter sine transform backward of an R4VEC.
c
c  Discussion:
c    
c    This routine is provided for illustration and testing.  It is inefficient
c    relative to optimized routines that use fast Fourier techniques.  
c
c    For 0 <= I <= N-1,
c
c      Y(I) = -2 Sum ( 1 <= J <= N-1 ) X(J) * sin ( PI * J * (I+1/2) / N )
c             - X(N) * cos ( pi * I )
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    10 June 2010
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    William Briggs, Van Emden Henson,
c    The Discrete Fourier Transform, 
c    SIAM, 1995,
c    QA403.5 B75
c
c  Parameters:
c
c    Input, integer N, the number of data values. 
c
c    Input, real X(N), the data sequence.
c
c    Output, real Y(0:N-1), the transformed data.
c
      implicit none

      integer n

      integer i
      integer j
      real pi
      parameter ( pi = 3.141592653589793E+00 )
      real theta
      real x(1:n)
      real y(0:n-1)

      do i = 0, n - 1
        y(i) = 0.0E+00
      end do

      do i = 0, n - 1
        do j = 1, n - 1

          theta = 0.5E+00 * pi * real ( j * ( 2 * i + 1 ) ) 
     &      / real ( n )
          y(i) = y(i) - 2.0E+00 * x(j) * sin ( theta  )

        end do

        theta = pi * real ( i )
        y(i) = y(i) - x(n) * cos ( theta )

      end do

      return
      end
      subroutine r4vec_sqstf ( n, x, y )

c*********************************************************************72
c
cc R4VEC_SQSTF computes a "slow" quarter sine transform forward of an R4VEC.
c
c  Discussion:
c    
c    This routine is provided for illustration and testing.  It is inefficient
c    relative to optimized routines that use fast Fourier techniques.  
c
c    For 1 <= I <= N,
c
c      Y(I) = -(1/N) Sum ( 0 <= J <= N-1 ) X(J) * sin ( PI * I * (J+1/2) / N )
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    10 June 2010
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    William Briggs, Van Emden Henson,
c    The Discrete Fourier Transform, 
c    SIAM, 1995,
c    QA403.5 B75
c
c  Parameters:
c
c    Input, integer N, the number of data values. 
c
c    Input, real X(0:N-1), the data sequence.
c
c    Output, real Y(N), the transformed data.
c
      implicit none

      integer n

      integer i
      integer j
      real pi
      parameter ( pi = 3.141592653589793E+00 )
      real theta
      real x(0:n-1)
      real y(n)

      do i = 1, n
        y(i) = 0.0E+00
      end do

      do i = 1, n
        do j = 0, n - 1
          theta = 0.5E+00 * pi * real ( i * ( 2 * j + 1 ) ) 
     &      / real ( n )
          y(i) = y(i) + x(j) * sin ( theta  )
        end do
      end do

      do i = 1, n
        y(i) = - y(i) / real ( n )
      end do

      return
      end
      subroutine r4vec_sst ( n, x, y )

c*********************************************************************72
c
cc R4VEC_SST computes a "slow" sine transform of an R4VEC.
c
c  Discussion:
c
c    This routine is provided for illustration and testing.  It is inefficient
c    relative to optimized routines that use fast Fourier techniques.  
c
c    For 1 <= I <= N,
c
c      Y(I) = Sum ( 1 <= J <= N ) X(J) * sin ( PI * I * J / ( N + 1 ) )
c
c    Applying the routine twice in succession should yield the original data,
c    multiplied by N / 2.  This is a good check for correctness and accuracy.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    10 June 2010
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the number of data values. 
c
c    Input, real X(N), the data sequence.
c
c    Output, real Y(N), the transformed data.
c
      implicit none

      integer n

      integer i
      integer j
      real pi
      parameter ( pi = 3.141592653589793E+00 )
      real theta(n)
      real x(n)
      real y(n)

      do i = 1, n
        theta(i) = pi * real ( i ) / real ( n + 1 )
      end do

      do i = 1, n
        y(i) = 0.0E+00
      end do

      do i = 1, n
        do j = 1, n
          y(j) = y(j) + 2.0E+00 * x(i) * sin ( real ( i ) * theta(j) )
        end do
      end do

      return
      end
      subroutine r4vec_uniform ( n, a, b, seed, r )

c*********************************************************************72
c
cc R4VEC_UNIFORM returns a scaled pseudorandom R4VEC.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    22 March 2005
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Paul Bratley, Bennett Fox, Linus Schrage,
c    A Guide to Simulation,
c    Second Edition,
c    Springer, 1987,
c    ISBN: 0387964673,
c    LC: QA76.9.C65.B73.
c
c    Bennett Fox,
c    Algorithm 647:
c    Implementation and Relative Efficiency of Quasirandom
c    Sequence Generators,
c    ACM Transactions on Mathematical Software,
c    Volume 12, Number 4, December 1986, pages 362-376.
c
c    Pierre L'Ecuyer,
c    Random Number Generation,
c    in Handbook of Simulation,
c    edited by Jerry Banks,
c    Wiley, 1998,
c    ISBN: 0471134031,
c    LC: T57.62.H37.
c
c    Peter Lewis, Allen Goodman, James Miller,
c    A Pseudo-Random Number Generator for the System/360,
c    IBM Systems Journal,
c    Volume 8, Number 2, 1969, pages 136-143.
c
c  Parameters:
c
c    Input, integer M, the number of entries in the vector.
c
c    Input, real A, B, the lower and upper limits.
c
c    Input/output, integer SEED, the "seed" value, which should NOT be 0.
c    On output, SEED has been updated.
c
c    Output, real R(N), the vector of pseudorandom values.
c
      implicit none

      integer n

      real a
      real b
      integer i
      integer i4_huge
      parameter ( i4_huge = 2147483647 )
      integer k
      integer seed
      real r(n)

      if ( seed .eq. 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'R4VEC_UNIFORM - Fatal error!'
        write ( *, '(a)' ) '  Input value of SEED = 0.'
        stop
      end if

      do i = 1, n

        k = seed / 127773

        seed = 16807 * ( seed - k * 127773 ) - k * 2836

        if ( seed .lt. 0 ) then
          seed = seed + i4_huge
        end if
    
        r(i) = a + ( b - a ) * real ( seed ) * 4.656612875E-10

      end do

      return
      end
      subroutine r8vec_print_part ( n, a, max_print, title )

c*********************************************************************72
c
cc R8VEC_PRINT_PART prints "part" of an R8VEC.
c
c  Discussion:
c
c    The user specifies MAX_PRINT, the maximum number of lines to print.
c
c    If N, the size of the vector, is no more than MAX_PRINT, then
c    the entire vector is printed, one entry per line.
c
c    Otherwise, if possible, the first MAX_PRINT-2 entries are printed,
c    followed by a line of periods suggesting an omission,
c    and the last entry.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    22 June 2010
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the number of entries of the vector.
c
c    Input, double precision A(N), the vector to be printed.
c
c    Input, integer MAX_PRINT, the maximum number of lines to print.
c
c    Input, character*(*) TITLE, a title.
c
      implicit none

      integer n

      double precision a(n)
      integer i
      integer max_print
      integer s_len_trim
      character*(*) title

      if ( max_print .le. 0 ) then
        return
      end if

      if ( n .le. 0 ) then
        return
      end if

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) title
      write ( *, '(a)' ) ' '

      if ( n .le. max_print ) then

        do i = 1, n
          write ( *, '(2x,i8,a1,1x,g14.6)' ) i, ':', a(i)
        end do

      else if ( 3 .le. max_print ) then

        do i = 1, max_print - 2
          write ( *, '(2x,i8,a1,1x,g14.6)' ) i, ':', a(i)
        end do

        write ( *, '(a)' ) '  ........  ..............'
        i = n

        write ( *, '(2x,i8,a1,1x,g14.6)' ) i, ':', a(i)

      else

        do i = 1, max_print - 1
          write ( *, '(2x,i8,a1,1x,g14.6)' ) i, ':', a(i)
        end do

        i = max_print

        write ( *, '(2x,i8,a1,1x,g14.6,a)' ) 
     &    i, ':', a(i), '...more entries...'

      end if

      return
      end
      subroutine r8vec_sct ( n, x, y )

c*********************************************************************72
c
cc R8VEC_SCT computes a "slow" cosine transform of an R8VEC.
c
c  Discussion:
c
c    This routine is provided for illustration and testing.  It is inefficient
c    relative to optimized routines that use fast Fourier techniques.
c
c      Y(1) = Sum ( 1 <= J <= N ) X(J) 
c
c      For 2 <= I <= N-1:
c
c        Y(I) = 2 * Sum ( 1 <= J <= N ) X(J) 
c          * cos ( PI * ( I - 1 ) * ( J - 1 ) / ( N - 1 ) )
c
c      Y(N) = Sum ( X(1:N:2) ) - Sum ( X(2:N:2) )
c
c    Applying the routine twice in succession should yield the original data,
c    multiplied by 2 * ( N + 1 ).  This is a good check for correctness 
c    and accuracy.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    12 July 2008
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the number of data values.
c
c    Input, double precision X(N), the data sequence.
c
c    Output, double precision Y(N), the transformed data.
c
      implicit none

      integer n

      double precision angle
      integer i
      integer j
      double precision pi
      parameter ( pi = 3.141592653589793D+00 )
      double precision x(n)
      double precision y(n)

      do i = 1, n

        y(i) = x(1) / 2.0D+00

        do j = 2, n - 1
          angle = pi 
     &      * dble ( mod ( ( i - 1 ) * ( j - 1 ), 2 * ( n - 1 ) ) ) 
     &      / dble ( n - 1 )
          y(i) = y(i) + x(j) * cos ( angle )
        end do

        j = n

        angle = pi 
     &    * dble ( mod ( ( i - 1 ) * ( j - 1 ), 2 * ( n - 1 ) ) ) 
     &    / dble ( n - 1 )

        y(i) = y(i) + x(n) * cos ( angle ) / 2.0D+00

      end do

      do i = 1, n
        y(i) = 2.0D+00 * y(i) * sqrt ( dble ( n ) / dble ( n - 1 ) )
      end do

      return
      end
      subroutine r8vec_sftb ( n, azero, a, b, r )

c*********************************************************************72
c
cc R8VEC_SFTB computes a "slow" backward Fourier transform of an R8VEC.
c
c  Discussion:
c
c    SFTB and SFTF are inverses of each other.  If we begin with data
c    AZERO, A, and B, and apply SFTB to it, and then apply SFTF to the
c    resulting R vector, we should get back the original AZERO, A and B.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    15 July 2008
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the number of data values. 
c
c    Input, double precision AZERO, the constant Fourier coefficient.
c
c    Input, double precision A(N/2), B(N/2), the Fourier coefficients.
c
c    Output, double precision R(N), the reconstructed data sequence.
c
      implicit none

      integer n

      double precision a(n/2)
      double precision azero
      double precision b(n/2)
      integer i
      integer k
      double precision r(n)
      double precision pi
      parameter ( pi = 3.141592653589793D+00 )
      double precision theta

      do i = 1, n
        r(i) = azero
      end do

      do i = 1, n
        do k = 1, n / 2
          theta = dble ( k * ( i - 1 ) * 2 ) * pi / dble ( n )
          r(i) = r(i) + a(k) * cos ( theta ) + b(k) * sin ( theta )
        end do
      end do

      return
      end
      subroutine r8vec_sftf ( n, r, azero, a, b )

c*********************************************************************72
c
cc R8VEC_SFTF computes a "slow" forward Fourier transform of an R8VEC.
c
c  Discussion:
c
c    SFTF and SFTB are inverses of each other.  If we begin with data
c    R and apply SFTB to it, and then apply SFTB to the resulting AZERO, 
c    A, and B, we should get back the original R.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    15 July 2008
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the number of data values.
c
c    Input, double precision R(N), the data to be transformed.
c
c    Output, double precision AZERO, = sum ( 1 <= I <= N ) R(I) / N.
c
c    Output, double precision A(N/2), B(N/2), the Fourier coefficients.
c
      implicit none

      integer n

      double precision a(1:n/2)
      double precision azero
      double precision b(1:n/2)
      integer i
      integer j
      double precision pi
      parameter ( pi = 3.141592653589793D+00 )
      double precision r(n)
      double precision theta

      azero = 0.0D+00
      do i = 1, n
        azero = azero + r(i)
      end do
      azero = azero / dble ( n )

      do i = 1, n / 2

        a(i) = 0.0D+00
        b(i) = 0.0D+00

        do j = 1, n
          theta = dble ( 2 * i * ( j - 1 ) ) * pi 
     &      / dble ( n )
          a(i) = a(i) + r(j) * cos ( theta )
          b(i) = b(i) + r(j) * sin ( theta )
        end do

        a(i) = a(i) / dble ( n )
        b(i) = b(i) / dble ( n )

        if ( i .ne. ( n / 2 ) ) then
          a(i) = 2.0D+00 * a(i)
          b(i) = 2.0D+00 * b(i)
        end if

      end do

      return
      end
      subroutine r8vec_sht ( n, a, b  )

c*********************************************************************72
c
cc R8VEC_SHT computes a "slow" Hartley transform of an R8VEC.
c
c  Discussion:
c
c    The discrete Hartley transform B of a set of data A is 
c
c      B(I) = 1/sqrt(N) * Sum ( 0 <= J <= N-1 ) A(J) * CAS(2*PI*I*J/N)
c
c    Here, the data and coefficients are indexed from 0 to N-1.
c
c    With the above normalization factor of 1/sqrt(N), the Hartley
c    transform is its own inverse.
c
c    This routine is provided for illustration and testing.  It is inefficient
c    relative to optimized routines.  
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    15 July 2008
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Ralph Hartley,
c    A More Symmetrical Fourier Analysis Applied to Transmission Problems,
c    Proceedings of the Institute of Radio Engineers,
c    Volume 30, pages 144-150, 1942.
c
c  Parameters:
c
c    Input, integer N, the number of data values.
c
c    Input, double precision A(0:N-1), the data to be transformed.
c
c    Output, double precision B(0:N-1), the transformed data.
c
      implicit none

      integer n

      double precision a(0:n-1)
      double precision b(0:n-1)
      integer i
      integer j
      double precision pi
      parameter ( pi = 3.141592653589793D+00 )
      double precision theta

      do i = 0, n - 1
        b(i) = 0.0D+00
      end do

      do i = 0, n - 1
        do j = 0, n - 1
          theta = 2.0D+00 * pi * dble ( mod ( i * j, n ) ) 
     &      / dble ( n )
          b(i) = b(i) + a(j) * ( cos ( theta ) + sin ( theta ) )
        end do
      end do

      do i = 0, n - 1
        b(i) = b(i) / sqrt ( dble ( n ) )
      end do

      return
      end
      subroutine r8vec_sqctb ( n, x, y )

c*********************************************************************72
c
cc R8VEC_SQCTB computes a "slow" quarter cosine transform backward of an R8VEC.
c
c  Discussion:
c
c    This routine is provided for illustration and testing.  It is inefficient
c    relative to optimized routines that use fast Fourier techniques.
c
c    For 0 <= I <= N-1,
c
c      Y(I) = X(0) + 2 Sum ( 1 <= J <= N-1 ) X(J) * cos ( PI * J * (I+1/2) / N )
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    16 July 2008
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    William Briggs, Van Emden Henson,
c    The Discrete Fourier Transform, 
c    SIAM, 1995,
c    LC: QA403.5 B75
c
c  Parameters:
c
c    Input, integer N, the number of data values.
c
c    Input, double precision X(0:N-1), the data sequence.
c
c    Output, double precision Y(0:N-1), the transformed data.
c
      implicit none

      integer n

      integer i
      integer j
      double precision pi
      parameter ( pi = 3.141592653589793D+00 )
      double precision theta
      double precision x(0:n-1)
      double precision y(0:n-1)

      do i = 0, n - 1
        y(i) = x(0)
      end do

      do i = 0, n - 1
        do j = 1, n - 1

          theta = 0.5D+00 * pi * dble ( j * ( 2 * i + 1 ) )
     &      / dble ( n )
          y(i) = y(i) + 2.0D+00 * x(j) * cos ( theta  )

        end do

      end do

      return
      end
      subroutine r8vec_sqctf ( n, x, y )

c*********************************************************************72
c
cc R8VEC_SQCTF computes a "slow" quarter cosine transform forward of an R8VEC.
c
c  Discussion:
c
c    This routine is provided for illustration and testing.  It is inefficient
c    relative to optimized routines that use fast Fourier techniques.
c
c    For 0 <= I <= N-1,
c
c      Y(I) = (1/N) Sum ( 0 <= J <= N-1 ) X(J) * cos ( PI * I * (J+1/2) / N )
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    16 July 2008
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    William Briggs, Van Emden Henson,
c    The Discrete Fourier Transform, 
c    SIAM, 1995,
c    QA403.5 B75
c
c  Parameters:
c
c    Input, integer N, the number of data values.
c
c    Input, double precision X(0:N-1), the data sequence.
c
c    Output, double precision Y(0:N-1), the transformed data.
c
      implicit none

      integer n

      integer i
      integer j
      double precision pi
      parameter ( pi = 3.141592653589793D+00 )
      double precision theta
      double precision x(0:n-1)
      double precision y(0:n-1)

      do i = 0, n - 1
        y(i) = 0.0D+00
      end do

      do i = 0, n - 1
        do j = 0, n - 1
          theta = 0.5D+00 * pi * dble ( i * ( 2 * j + 1 ) ) 
     &      / dble ( n )
          y(i) = y(i) + x(j) * cos ( theta  )
        end do
      end do

      do i = 0, n - 1
        y(i) = y(i) / dble ( n )
      end do

      return
      end
      subroutine r8vec_sqstb ( n, x, y )

c*********************************************************************72
c
cc R8VEC_SQSTB computes a "slow" quarter sine transform backward of an R8VEC.
c
c  Discussion:
c    
c    This routine is provided for illustration and testing.  It is inefficient
c    relative to optimized routines that use fast Fourier techniques.  
c
c    For 0 <= I <= N-1,
c
c      Y(I) = -2 Sum ( 1 <= J <= N-1 ) X(J) * sin ( PI * J * (I+1/2) / N )
c             - X(N) * cos ( pi * I )
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    16 July 2008
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    William Briggs, Van Emden Henson,
c    The Discrete Fourier Transform, 
c    SIAM, 1995,
c    QA403.5 B75
c
c  Parameters:
c
c    Input, integer N, the number of data values. 
c
c    Input, double precision X(N), the data sequence.
c
c    Output, double precision Y(0:N-1), the transformed data.
c
      implicit none

      integer n

      integer i
      integer j
      double precision pi
      parameter ( pi = 3.141592653589793D+00 )
      double precision theta
      double precision x(1:n)
      double precision y(0:n-1)

      do i = 0, n - 1
        y(i) = 0.0D+00
      end do

      do i = 0, n - 1
        do j = 1, n - 1

          theta = 0.5D+00 * pi * dble ( j * ( 2 * i + 1 ) ) 
     &      / dble ( n )
          y(i) = y(i) - 2.0D+00 * x(j) * sin ( theta  )

        end do

        theta = pi * dble ( i )
        y(i) = y(i) - x(n) * cos ( theta )

      end do

      return
      end
      subroutine r8vec_sqstf ( n, x, y )

c*********************************************************************72
c
cc R8VEC_SQSTF computes a "slow" quarter sine transform forward of an R8VEC.
c
c  Discussion:
c    
c    This routine is provided for illustration and testing.  It is inefficient
c    relative to optimized routines that use fast Fourier techniques.  
c
c    For 1 <= I <= N,
c
c      Y(I) = -(1/N) Sum ( 0 <= J <= N-1 ) X(J) * sin ( PI * I * (J+1/2) / N )
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    16 July 2008
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    William Briggs, Van Emden Henson,
c    The Discrete Fourier Transform, 
c    SIAM, 1995,
c    QA403.5 B75
c
c  Parameters:
c
c    Input, integer N, the number of data values. 
c
c    Input, double precision X(0:N-1), the data sequence.
c
c    Output, double precision Y(N), the transformed data.
c
      implicit none

      integer n

      integer i
      integer j
      double precision pi
      parameter ( pi = 3.141592653589793D+00 )
      double precision theta
      double precision x(0:n-1)
      double precision y(n)

      do i = 1, n
        y(i) = 0.0D+00
      end do

      do i = 1, n
        do j = 0, n - 1
          theta = 0.5D+00 * pi * dble ( i * ( 2 * j + 1 ) ) 
     &      / dble ( n )
          y(i) = y(i) + x(j) * sin ( theta  )
        end do
      end do

      do i = 1, n
        y(i) = - y(i) / dble ( n )
      end do

      return
      end
      subroutine r8vec_sst ( n, x, y )

c*********************************************************************72
c
cc R8VEC_SST computes a "slow" sine transform of an R8VEC.
c
c  Discussion:
c
c    This routine is provided for illustration and testing.  It is inefficient
c    relative to optimized routines that use fast Fourier techniques.  
c
c    For 1 <= I <= N,
c
c      Y(I) = Sum ( 1 <= J <= N ) X(J) * sin ( PI * I * J / ( N + 1 ) )
c
c    Applying the routine twice in succession should yield the original data,
c    multiplied by N / 2.  This is a good check for correctness and accuracy.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    16 July 2008
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the number of data values. 
c
c    Input, double precision X(N), the data sequence.
c
c    Output, double precision Y(N), the transformed data.
c
      implicit none

      integer n

      integer i
      integer j
      double precision pi
      parameter ( pi = 3.141592653589793D+00 )
      double precision theta(n)
      double precision x(n)
      double precision y(n)

      do i = 1, n
        theta(i) = pi * dble ( i ) / dble ( n + 1 )
      end do

      do i = 1, n
        y(i) = 0.0D+00
      end do

      do i = 1, n
        do j = 1, n
          y(j) = y(j) + 2.0D+00 * x(i) * sin ( dble ( i ) * theta(j) )
        end do
      end do

      return
      end
      subroutine r8vec_uniform ( n, a, b, seed, r )

c*********************************************************************72
c
cc R8VEC_UNIFORM returns a scaled pseudorandom R8VEC.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    29 January 2005
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
c    Bennett Fox,
c    Algorithm 647:
c    Implementation and Relative Efficiency of Quasirandom
c    Sequence Generators,
c    ACM Transactions on Mathematical Software,
c    Volume 12, Number 4, pages 362-376, 1986.
c
c    Peter Lewis, Allen Goodman, James Miller,
c    A Pseudo-Random Number Generator for the System/360,
c    IBM Systems Journal,
c    Volume 8, pages 136-143, 1969.
c
c  Parameters:
c
c    Input, integer M, the number of entries in the vector.
c
c    Input, double precision A, B, the lower and upper limits.
c
c    Input/output, integer SEED, the "seed" value, which should NOT be 0.
c    On output, SEED has been updated.
c
c    Output, double precision R(N), the vector of pseudorandom values.
c
      implicit none

      integer n

      double precision a
      double precision b
      integer i
      integer k
      integer seed
      double precision r(n)

      if ( seed .eq. 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'R8VEC_UNIFORM - Fatal error!'
        write ( *, '(a)' ) '  Input value of SEED = 0.'
        stop
      end if

      do i = 1, n

        k = seed / 127773

        seed = 16807 * ( seed - k * 127773 ) - k * 2836

        if ( seed .lt. 0 ) then
          seed = seed + 2147483647
        end if

        r(i) = a + ( b - a ) * dble ( seed ) * 4.656612875D-10

      end do

      return
      end
      function s_len_trim ( s )

c*********************************************************************72
c
cc S_LEN_TRIM returns the length of a string to the last nonblank.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    05 March 2004
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, character*(*) S, a string.
c
c    Output, integer S_LEN_TRIM, the length of the string to the last nonblank.
c
      implicit none

      integer i
      character*(*) s
      integer s_len_trim

      do i = len ( s ), 1, -1

        if ( s(i:i) .ne. ' ' ) then
          s_len_trim = i
          return
        end if

      end do

      s_len_trim = 0

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
