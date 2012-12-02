      subroutine asm_enum ( n, asm_num )

c*********************************************************************72
c
cc ASM_ENUM returns the number of alternating sign matrices of a given order.
c
c  Discussion:
c
c    N     ASM_NUM
c
c    0       1
c    1       1
c    2       2
c    3       7
c    4      42
c    5     429
c    6    7436
c    7  218348
c
c    A direct formula is
c
c      ASM_NUM ( N ) = product ( 0 <= I <= N-1 ) ( 3 * I + 1 )! / ( N + I )!
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    03 January 2007
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the order of the matrices.
c
c    Output, integer ASM_NUM, the number of alternating sign matrices of
c    order N.
c
      implicit none

      integer n

      integer a(n+1)
      integer asm_num
      integer b(n+1)
      integer c(n+1)
      integer i
      integer nn

      asm_num = 0

      nn = n + 1

      if ( n + 1 .le. 0 ) then
        return
      end if
c
c  Row 1
c
      a(1) = 1

      if ( n + 1 .eq. 1 ) then
        asm_num = a(1)
        return
      end if
c
c  Row 2
c
      a(1) = 1

      if ( n + 1 .eq. 2 ) then
        asm_num = a(1)
        return
      end if

      b(1) = 2
      c(1) = 2
      a(2) = 1
c
c  Row 3 and on.
c
      do nn = 3, n

        b(nn-1) = nn
        do i = nn-2, 2, -1
          b(i) = b(i) + b(i-1)
        end do
        b(1) = 2

        c(nn-1) = 2
        do i = nn-2, 2, -1
          c(i) = c(i) + c(i-1)
        end do
        c(1) = nn

        do i = 2, nn-1
          a(1) = a(1) + a(i)
        end do

        do i = 2, nn
          a(i) = a(i-1) * c(i-1) / b(i-1)
        end do

      end do

      asm_num = 0
      do i = 1, n
        asm_num = asm_num + a(i)
      end do

      return
      end
      subroutine asm_triangle ( n, a )

c*********************************************************************72
c
cc ASM_TRIANGLE returns a row of the alternating sign matrix triangle.
c
c  Discussion:
c
c    The first seven rows of the triangle are as follows:
c
c          1      2      3      4      5      6     7
c
c    0     1
c    1     1      1
c    2     2      3      2
c    3     7     14     14      7
c    4    42    105    135    105     42
c    5   429   1287   2002   2002   1287    429
c    6  7436  26026  47320  56784  47320  26026  7436
c
c    For a given N, the value of A(J) represents entry A(I,J) of
c    the triangular matrix, and gives the number of alternating sign matrices
c    of order N in which the (unique) 1 in row 1 occurs in column J.
c
c    Thus, of alternating sign matrices of order 3, there are
c    2 with a leading 1 in column 1:
c
c      1 0 0  1 0 0
c      0 1 0  0 0 1
c      0 0 1  0 1 0
c
c    3 with a leading 1 in column 2, and
c
c      0 1 0  0 1 0  0 1 0
c      1 0 0  0 0 1  1-1 1
c      0 0 1  1 0 0  0 1 0
c
c    2 with a leading 1 in column 3:
c
c      0 0 1  0 0 1
c      1 0 0  0 1 0
c      0 1 0  1 0 0
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    03 January 2007
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the desired row.
c
c    Output, integer A(N+1), the entries of the row.
c
      implicit none

      integer n

      integer a(n+1)
      integer b(n+1)
      integer c(n+1)
      integer i
      integer i4vec_sum
      integer nn

      if ( n + 1 .le. 0 ) then
        return
      end if
c
c  Row 1
c
      a(1) = 1

      if ( n + 1 .eq. 1 ) then
        return
      end if
c
c  Row 2
c
      nn = 2
      b(1) = 2
      c(1) = nn

      a(1) = i4vec_sum ( nn-1, a )
      do i = 2, nn
        a(i) = a(i-1) * c(i-1) / b(i-1)
      end do

      if ( n + 1 .eq. 2 ) then
        return
      end if
c
c  Row 3 and on.
c
      do nn = 3, n + 1

        b(nn-1) = nn
        do i = nn-2, 2, -1
          b(i) = b(i) + b(i-1)
        end do
        b(1) = 2

        c(nn-1) = 2
        do i = nn-2, 2, -1
          c(i) = c(i) + c(i-1)
        end do
        c(1) = nn

        a(1) = i4vec_sum ( nn-1, a )
        do i = 2, nn
          a(i) = a(i-1) * c(i-1) / b(i-1)
        end do

      end do

      return
      end
      subroutine bell ( n, b )

c*********************************************************************72
c
cc BELL returns the Bell numbers from 0 to N.
c
c  Discussion:
c
c    The Bell number B(N) is defined as the number of partitions (of
c    any size) of a set of N distinguishable objects.
c
c    A partition of a set is a division of the objects of the set into
c    subsets.
c
c    The Bell number B(N) is the number of restricted growth functions
c    on N.
c
c    Note that the Stirling numbers of the second kind, S^m_n, count the
c    number of partitions of N objects into M classes, and so it is
c    true that
c
c      B(N) = S^1_N + S^2_N + ... + S^N_N.
c
c  Example:
c
c    There are 15 partitions of a set of 4 objects:
c
c      (1234), (123)(4), (124)(3), (12)(34), (12)(3)(4),
c      (134)(2), (13)(24), (13)(2)(4), (14)(23), (1)(234),
c      (1)(23)(4), (14)(2)(3), (1)(24)(3), (1)(2)(34), (1)(2)(3)(4)
c
c    and so B(4) = 15.
c
c  First values:
c
c     N         B(N)
c     0           1
c     1           1
c     2           2
c     3           5
c     4          15
c     5          52
c     6         203
c     7         877
c     8        4140
c     9       21147
c    10      115975
c
c  Recursion:
c
c    B(I) = sum ( 1 <= J <= I ) Binomial ( I-1, J-1 ) * B(I-J)
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    03 January 2007
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the number of Bell numbers desired.
c
c    Output, integer B(0:N), the Bell numbers from 0 to N.
c
      implicit none

      integer n

      integer b(0:n)
      integer combo
      integer i
      integer j

      b(0) = 1

      do i = 1, n
        b(i) = 0
        do j = 1, i
          call combin2 ( i-1, j-1, combo )
          b(i) = b(i) + combo * b(i-j)
        end do
      end do

      return
      end
      subroutine bell_values ( n_data, n, c )

c*********************************************************************72
c
cc BELL_VALUES returns some values of the Bell numbers for testing.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    03 January 2007
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Milton Abramowitz, Irene Stegun,
c    Handbook of Mathematical Functions,
c    US Department of Commerce, 1964,
c    ISBN: 0-486-61272-4,
c    LC: QA47.A34.
c
c  Parameters:
c
c    Input/output, integer N_DATA.
c    On input, if N_DATA is 0, the first test data is returned, and N_DATA
c    is set to 1.  On each subsequent call, the input value of N_DATA is
c    incremented and that test data item is returned, if available.  When
c    there is no more test data, N_DATA is set to 0.
c
c    Output, integer N, the order of the Bell number.
c
c    Output, integer C, the value of the Bell number.
c
      implicit none

      integer nmax
      parameter ( nmax = 11 )

      integer c
      integer c_vec(nmax)
      integer n
      integer n_data
      integer n_vec(nmax)

      save c_vec
      save n_vec 

      data c_vec /
     &  1, 1, 2, 5, 15, 52, 203, 877, 4140, 21147, 115975 /
      data n_vec /
     &   0,  1,  2,  3,  4, 5,  6,  7,  8,  9,  10 /

      if ( n_data .lt. 0 ) then
        n_data = 0
      end if

      n_data = n_data + 1

      if ( nmax .lt. n_data ) then
        n_data = 0
        n = 0
        c = 0
      else
        n = n_vec(n_data)
        c = c_vec(n_data)
      end if

      return
      end
      subroutine binary_vector_next ( n, bvec )

c*********************************************************************72
c
cc BINARY_VECTOR_NEXT generates the next binary vector.
c
c  Discussion:
c
c    A binary vector is a vector whose entries are 0 or 1.
c
c    The user inputs an initial zero vector to start.  The program returns
c    the "next" vector.
c
c    The vectors are produced in the order:
c
c    ( 0, 0, 0, ..., 0 )
c    ( 1, 0, 0, ..., 0 ) 
c    ( 0, 1, 0, ..., 0 )
c    ( 1, 1, 0, ..., 0 )
c    ( 0, 0, 1, ..., 0 )
c    ( 1, 0, 1, ..., 0 )
c               ...
c    ( 1, 1, 1, ..., 1)
c
c    and the "next" vector after (1,1,...,1) is (0,0,...,0).  That is,
c    we allow wrap around.
c
c  Example:
c
c    N = 3
c
c    Input      Output
c    -----      ------
c    0 0 0  =>  1 0 0
c    1 0 0  =>  0 1 0
c    0 1 0  =>  1 1 0
c    1 1 0  =>  0 0 1
c    0 0 1  =>  1 0 1
c    1 0 1  =>  0 1 1
c    0 1 1  =>  1 1 1
c    1 1 1  =>  0 0 0
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    15 November 2009
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the dimension of the vectors.
c
c    Input/output, integer BVEC(N), on output, the successor 
c    to the input vector.
c
      implicit none

      integer n
  
      integer bvec(n)
      integer i
  
      do i = 1, n

        if ( bvec(i) .eq. 1 ) then
          bvec(i) = 0
        else 
          bvec(i) = 1
          return
        end if

      end do

      return
      end
      subroutine bvec_add ( n, bvec1, bvec2, bvec3 )

c*********************************************************************72
c
cc BVEC_ADD adds two (signed) binary vectors.
c
c  Discussion:
c
c    A BVEC is a vector of binary digits representing an integer.  
c
c    BVEC(1) is 0 for positive values and 1 for negative values, which
c    are stored in 2's complement form.
c
c    For positive values, BVEC(N) contains the units digit, BVEC(N-1)
c    the coefficient of 2, BVEC(N-2) the coefficient of 4 and so on,
c    so that printing the digits in order gives the binary form of the number.
c
c  Example:
c
c    N = 5
c
c      BVEC1       +   BVEC2       =   BVEC3
c
c    ( 0 0 0 0 1 ) + ( 0 0 0 1 1 ) = ( 0 0 1 0 0 )
c
c              1   +           3   =           4
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    20 March 2009
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the length of the vectors.
c
c    Input, integer BVEC1(N), BVEC2(N), the vectors to be added.
c
c    Output, integer BVEC3(N), the sum of the two input vectors.
c
      implicit none

      integer n

      integer base
      parameter ( base = 2 )
      integer bvec1(n)
      integer bvec2(n)
      integer bvec3(n)
      integer i
      logical overflow

      overflow = .false.
    
      do i = 1, n
        bvec3(i) = bvec1(i) + bvec2(i)
      end do

      do i = 1, n

10      continue

        if ( base .le. bvec3(i) ) then

          bvec3(i) = bvec3(i) - base

          if ( 1 .lt. i ) then
            bvec3(i-1) = bvec3(i-1) + 1
          else
            overflow = .true.
          end if

        end if

      end do

      return
      end
      subroutine bvec_and ( n, bvec1, bvec2, bvec3 )

c*********************************************************************72
c
cc BVEC_AND computes the AND of two binary vectors.
c
c  Discussion:
c
c    A BVEC is a vector of binary digits representing an integer.  
c
c    BVEC(1) is 0 for positive values and 1 for negative values, which
c    are stored in 2's complement form.
c
c    For positive values, BVEC(N) contains the units digit, BVEC(N-1)
c    the coefficient of 2, BVEC(N-2) the coefficient of 4 and so on,
c    so that printing the digits in order gives the binary form of the number.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    03 January 2007
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the length of the vectors.
c
c    Input, integer BVEC1(N), BVEC2(N), the binary vectors.
c
c    Input, integer BVEC3(N), the AND of the two vectors.
c
      implicit none

      integer n

      integer bvec1(n)
      integer bvec2(n)
      integer bvec3(n)
      integer i

      do i = 1, n
        bvec3(i) = min ( bvec1(i), bvec2(i) )
      end do

      return
      end
      subroutine bvec_check ( n, bvec, ierror )

c*********************************************************************72
c
cc BVEC_CHECK checks a binary vector.
c
c  Discussion:
c
c    A BVEC is a vector of binary digits representing an integer.  
c
c    BVEC(1) is 0 for positive values and 1 for negative values, which
c    are stored in 2's complement form.
c
c    For positive values, BVEC(N) contains the units digit, BVEC(N-1)
c    the coefficient of 2, BVEC(N-2) the coefficient of 4 and so on,
c    so that printing the digits in order gives the binary form of the number.
c
c    The only check made is that the entries are all 0 or 1.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    03 January 2007
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the length of the vectors.
c
c    Input, integer BVEC(N), the vector to be checked.
c
c    Output, integer IERROR, is nonzero if an error occurred.
c
      implicit none

      integer n

      integer base
      parameter ( base = 2 )
      integer bvec(n)
      integer i
      integer ierror

      ierror = 0

      do i = 1, n
        if ( bvec(i) .lt. 0 .or. base .le. bvec(i) ) then
          ierror = i
          return
        end if
      end do

      return
      end
      subroutine bvec_complement2 ( n, bvec1, bvec2 )

c*********************************************************************72
c
cc BVEC_COMPLEMENT2 computes the two's complement of a binary vector.
c
c  Discussion:
c
c    A BVEC is a vector of binary digits representing an integer.  
c
c    BVEC(1) is 0 for positive values and 1 for negative values, which
c    are stored in 2's complement form.
c
c    For positive values, BVEC(N) contains the units digit, BVEC(N-1)
c    the coefficient of 2, BVEC(N-2) the coefficient of 4 and so on,
c    so that printing the digits in order gives the binary form of the number.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    20 March 2009
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the length of the vectors.
c
c    Input, integer BVEC1(N), the vector to be complemented.
c
c    Output, integer BVEC2(N), the two's complemented vector.
c
      implicit none

      integer n
      integer n_max
      parameter ( n_max = 100 )

      integer base
      parameter ( base = 2 )
      integer bvec1(n)
      integer bvec2(n)
      integer bvec3(n_max)
      integer bvec4(n_max)
      integer i

      if ( n_max .lt. n ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'BVEC_COMPLEMENT2 - Fatal error!'
        write ( *, '(a)' ) '  Internal size limit N_MAX exceeded.'
        stop
      end if

      do i = 1, n
        bvec3(i) = ( base - 1 ) - bvec1(i)
      end do

      do i = 1, n - 1
        bvec4(i) = 0
      end do
      bvec4(n) = 1

      call bvec_add ( n, bvec3, bvec4, bvec2 )

      return
      end
      subroutine bvec_mul ( n, bvec1, bvec2, bvec3 )

c*********************************************************************72
c
cc BVEC_MUL computes the product of two binary vectors.
c
c  Discussion:
c
c    A BVEC is a vector of binary digits representing an integer.  
c
c    BVEC(1) is 0 for positive values and 1 for negative values, which
c    are stored in 2's complement form.
c
c    For positive values, BVEC(N) contains the units digit, BVEC(N-1)
c    the coefficient of 2, BVEC(N-2) the coefficient of 4 and so on,
c    so that printing the digits in order gives the binary form of the number.
c
c    Since the user may want to make calls like
c
c      call bvec_mul ( n, bvec1, bvec1, bvec3 )
c    or even
c      call bvec_mul ( n, bvec1, bvec1, bvec1 )
c
c    we need to copy the arguments, work on them, and then copy out the result.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    20 March 2009
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the length of the vectors.
c
c    Input, integer BVEC1(N), BVEC2(N), the vectors to be multiplied.
c
c    Output, integer BVEC3(N), the product of the two input vectors.
c
      implicit none

      integer n
      integer n_max
      parameter ( n_max = 100 )

      integer base
      parameter ( base = 2 )
      integer carry
      integer bvec1(n)
      integer bvec2(n)
      integer bvec3(n)
      integer bveca(n_max)
      integer bvecb(n_max)
      integer bvecc(n_max)
      integer i
      integer j
      integer product_sign

      if ( n_max .lt. n ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'BVEC_MUL - Fatal error!'
        write ( *, '(a)' ) '  Internal size limit N_MAX exceeded.'
        stop
      end if
c
c  Copy the input.
c
      do i = 1, n
        bveca(i) = bvec1(i)
      end do

      do i = 1, n
        bvecb(i) = bvec2(i)
      end do
c
c  Record the sign of the product.
c  Make the factors positive.
c
      product_sign = 1

      if ( bveca(n) .ne. 0 ) then
        product_sign = - product_sign
        call bvec_complement2 ( n, bveca, bveca )
      end if

      if ( bvecb(n) .ne. 0 ) then
        product_sign = - product_sign
        call bvec_complement2 ( n, bvecb, bvecb )
      end if

      do i = 1, n
        bvecc(i) = 0
      end do
c
c  Multiply.
c
      do i = 2, n
        do j = 2, n + 2 - i
          bvecc(j) = bvecc(j) + bveca(n+2-i) * bvecb(j+i-2)
        end do
      end do
c
c  Take care of carries.
c
      do i = n, 2, -1

        carry = bvecc(i) / base
        bvecc(i) = bvecc(i) - carry * base
c
c  Unlike the case of BVEC_ADD, we do NOT allow carries into
c  the sign position when multiplying.
c
        if ( 2 .lt. i ) then
          bvecc(i-1) = bvecc(i-1) + carry
        end if

      end do
c
c  Take care of the sign of the product.
c
      if ( product_sign .lt. 0 ) then
        call bvec_complement2 ( n, bvecc, bvecc )
      end if
c
c  Copy the output.
c
      do i = 1, n
        bvec3(i) = bvecc(i)
      end do

      return
      end
      subroutine bvec_next ( n, bvec )

c*********************************************************************72
c
cc BVEC_NEXT generates the next binary vector.
c
c  Discussion:
c
c    A BVEC is a vector of binary digits representing an integer.  
c
c    BVEC(1) is 0 for positive values and 1 for negative values, which
c    are stored in 2's complement form.
c
c    For positive values, BVEC(N) contains the units digit, BVEC(N-1)
c    the coefficient of 2, BVEC(N-2) the coefficient of 4 and so on,
c    so that printing the digits in order gives the binary form of the number.
c
c    The vectors have the order
c
c      (0,0,...,0),
c      (0,0,...,1),
c      ...
c      (1,1,...,1)
c
c    and the "next" vector after (1,1,...,1) is (0,0,...,0).  That is,
c    we allow wrap around.
c
c  Example:
c
c    N = 3
c
c    Input      Output
c    -----      ------
c    0 0 0  =>  0 0 1
c    0 0 1  =>  0 1 0
c    0 1 0  =>  0 1 1
c    0 1 1  =>  1 0 0
c    1 0 0  =>  1 0 1
c    1 0 1  =>  1 1 0
c    1 1 0  =>  1 1 1
c    1 1 1  =>  0 0 0
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    02 May 2008
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the dimension of the vectors.
c
c    Input/output, integer BVEC(N), on output, the successor to the
c    input vector.
c
      implicit none

      integer n

      integer bvec(n)
      integer i

      do i = n, 1, -1

        if ( bvec(i) == 0 ) then
          bvec(i) = 1
          return
        end if

        bvec(i) = 0

      end do

      return
      end
      subroutine bvec_not ( n, bvec1, bvec2 )

c*********************************************************************72
c
cc BVEC_NOT "negates" or takes the 1's complement of a binary vector.
c
c  Discussion:
c
c    A BVEC is a vector of binary digits representing an integer.  
c
c    BVEC(1) is 0 for positive values and 1 for negative values, which
c    are stored in 2's complement form.
c
c    For positive values, BVEC(N) contains the units digit, BVEC(N-1)
c    the coefficient of 2, BVEC(N-2) the coefficient of 4 and so on,
c    so that printing the digits in order gives the binary form of the number.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    03 January 2007
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the length of the vectors.
c
c    Input, integer BVEC1(N), the vector to be negated.
c
c    Output, integer BVEC2(N), the negated vector.
c
      implicit none

      integer n

      integer base
      parameter ( base = 2 )
      integer bvec1(n)
      integer bvec2(n)
      integer i
 
      do i = 1, n
        bvec2(i) = ( base - 1 ) - bvec1(i)
      end do

      return
      end
      subroutine bvec_or ( n, bvec1, bvec2, bvec3 )

c*********************************************************************72
c
cc BVEC_OR computes the inclusive OR of two binary vectors.
c
c  Discussion:
c
c    A BVEC is a vector of binary digits representing an integer.  
c
c    BVEC(1) is 0 for positive values and 1 for negative values, which
c    are stored in 2's complement form.
c
c    For positive values, BVEC(N) contains the units digit, BVEC(N-1)
c    the coefficient of 2, BVEC(N-2) the coefficient of 4 and so on,
c    so that printing the digits in order gives the binary form of the number.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    03 January 2007
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the length of the vectors.
c
c    Input, integer BVEC1(N), BVEC2(N), the binary vectors.
c
c    Input, integer BVEC3(N), the inclusive OR of the two vectors.
c
      implicit none

      integer n

      integer bvec1(n)
      integer bvec2(n)
      integer bvec3(n)
      integer i

      do i = 1, n
        bvec3(i) = max ( bvec1(i), bvec2(i) )
      end do

      return
      end
      subroutine bvec_print ( n, bvec, title )

c*********************************************************************72
c
cc BVEC_PRINT prints a binary integer vector, with an optional title.
c
c  Discussion:
c
c    A BVEC is a vector of binary digits representing an integer.  
c
c    BVEC(1) is 0 for positive values and 1 for negative values, which
c    are stored in 2's complement form.
c
c    For positive values, BVEC(N) contains the units digit, BVEC(N-1)
c    the coefficient of 2, BVEC(N-2) the coefficient of 4 and so on,
c    so that printing the digits in order gives the binary form of the number.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    20 March 2009
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the number of components of the vector.
c
c    Input, integer BVEC(N), the vector to be printed.
c
c    Input, character*(*) TITLE, a title to be printed first.
c    TITLE may be blank.
c
      implicit none

      integer n

      integer bvec(n)
      integer ihi
      integer ilo
      integer i
      integer s_len
      integer s_len_trim
      character * ( * ) title

      s_len = s_len_trim ( title )

      if ( 0 .lt. s_len ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) title(1:s_len)
        write ( *, '(a)' ) ' '
      end if

      do ilo = 1, n, 70
        ihi = min ( ilo + 70 - 1, n )
        write ( *, '(2x,80i1)' ) ( bvec(i), i = ilo, ihi )
      end do

      return
      end
      subroutine bvec_reverse ( n, bvec1, bvec2 )

c*********************************************************************72
c
cc BVEC_REVERSE reverses a binary vector.
c
c  Discussion:
c
c    A BVEC is a vector of binary digits representing an integer.  
c
c    BVEC(1) is 0 for positive values and 1 for negative values, which
c    are stored in 2's complement form.
c
c    For positive values, BVEC(N) contains the units digit, BVEC(N-1)
c    the coefficient of 2, BVEC(N-2) the coefficient of 4 and so on,
c    so that printing the digits in order gives the binary form of the number.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    03 January 2007
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the length of the vectors.
c
c    Input, integer BVEC1(N), the vector to be reversed.
c
c    Output, integer BVEC2(N), the reversed vector.
c
      implicit none

      integer n

      integer bvec1(n)
      integer bvec2(n)
      integer i

      do i = 1, n
        bvec2(i) = bvec1(n+1-i)
      end do

      return
      end
      subroutine bvec_sub ( n, bvec1, bvec2, bvec3 )

c*********************************************************************72
c
cc BVEC_SUB subtracts two binary vectors.
c
c  Discussion:
c
c    A BVEC is a vector of binary digits representing an integer.  
c
c    BVEC(1) is 0 for positive values and 1 for negative values, which
c    are stored in 2's complement form.
c
c    For positive values, BVEC(N) contains the units digit, BVEC(N-1)
c    the coefficient of 2, BVEC(N-2) the coefficient of 4 and so on,
c    so that printing the digits in order gives the binary form of the number.
c
c  Example:
c
c    N = 4
c
c    BVEC1         BVEC2         BVEC3
c    -------       -------       -------
c    0 1 0 0   -   0 0 0 1   =   0 0 1 1
c          4             1             3
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    20 March 2009
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the length of the vectors.
c
c    Input, integer BVEC1(N), BVEC2(N), the vectors to be subtracted.
c
c    Output, integer BVEC3(N), the value of BVEC1 - BVEC2.
c
      implicit none

      integer n
      integer n_max
      parameter ( n_max = 100 )

      integer bvec1(n)
      integer bvec2(n)
      integer bvec3(n)
      integer bvec4(n_max)

      if ( n_max .lt. n ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'BVEC_SUB - Fatal error!'
        write ( *, '(a)' ) '  Internal size limit N_MAX exceeded.'
        stop
      end if

      call bvec_complement2 ( n, bvec2, bvec4 )

      call bvec_add ( n, bvec1, bvec4, bvec3 )

      return
      end
      subroutine bvec_to_i4 ( n, bvec, i4 )

c*********************************************************************72
c
cc BVEC_TO_I4 makes an integer from a (signed) binary vector.
c
c  Discussion:
c
c    A BVEC is a vector of binary digits representing an integer.  
c
c    BVEC(1) is 0 for positive values and 1 for negative values, which
c    are stored in 2's complement form.
c
c    For positive values, BVEC(N) contains the units digit, BVEC(N-1)
c    the coefficient of 2, BVEC(N-2) the coefficient of 4 and so on,
c    so that printing the digits in order gives the binary form of the number.
c
c  Example:
c
c         BVEC   binary  I
c    ----------  -----  --
c    1  2  3  4
c    ----------
c    0  0  0  1       1  1
c    0  0  1  0      10  2
c    1  1  0  0    -100 -4
c    0  1  0  0     100  4
c    1  0  0  1    -111 -9
c    1  1  1  1      -0  0
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    20 March 2009
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the dimension of the vector.
c
c    Input, integer BVEC(N), the binary representation.
c
c    Output, integer I4, the integer.
c
      implicit none

      integer n
      integer n_max
      parameter ( n_max = 100 )

      integer base
      parameter ( base = 2 )
      integer bvec(n)
      integer bvec2(n_max)
      integer i
      integer i_sign
      integer i4

      if ( n_max .lt. n ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'BVEC_TO_I4 - Fatal error!'
        write ( *, '(a)' ) '  Internal size limit N_MAX exceeded.'
        stop
      end if

      do i = 1, n
        bvec2(i) = bvec(i)
      end do

      if ( bvec2(1) .eq. base - 1 ) then
        i_sign = -1
        bvec2(1) = 0
        call bvec_complement2 ( n - 1, bvec2(2:n), bvec2(2:n) )
      else 
        i_sign = 1
      end if

      i4 = 0
      do i = 2, n
        i4 = base * i4 + bvec2(i)
      end do

      i4 = i_sign * i4

      return
      end
      subroutine bvec_xor ( n, bvec1, bvec2, bvec3 )

c*********************************************************************72
c
cc BVEC_XOR computes the exclusive OR of two binary vectors.
c
c  Discussion:
c
c    A BVEC is a vector of binary digits representing an integer.  
c
c    BVEC(1) is 0 for positive values and 1 for negative values, which
c    are stored in 2's complement form.
c
c    For positive values, BVEC(N) contains the units digit, BVEC(N-1)
c    the coefficient of 2, BVEC(N-2) the coefficient of 4 and so on,
c    so that printing the digits in order gives the binary form of the number.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    03 January 2007
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the length of the vectors.
c
c    Input, integer BVEC1(N), BVEC2(N), the binary vectors to be XOR'ed.
c
c    Input, integer BVEC3(N), the exclusive OR of the two vectors.
c
      implicit none

      integer n

      integer bvec1(n)
      integer bvec2(n)
      integer bvec3(n)
      integer i

      do i = 1, n
        bvec3(i) = mod ( bvec1(i) + bvec2(i), 2 )
      end do

      return
      end
      subroutine catalan ( n, c )

c*********************************************************************72
c
cc CATALAN computes the Catalan numbers, from C(0) to C(N).
c
c  First values:
c
c     C(0)     1
c     C(1)     1
c     C(2)     2
c     C(3)     5
c     C(4)    14
c     C(5)    42
c     C(6)   132
c     C(7)   429
c     C(8)  1430
c     C(9)  4862
c    C(10) 16796
c
c  Formula:
c
c    C(N) = (2*N)c / ( (N+1) * (Nc) * (Nc) )
c         = 1 / (N+1) * COMB ( 2N, N )
c         = 1 / (2N+1) * COMB ( 2N+1, N+1).
c
c  Recursion:
c
c    C(N) = 2 * (2*N-1) * C(N-1) / (N+1)
c    C(N) = sum ( 1 <= I <= N-1 ) C(I) * C(N-I)
c
c  Discussion:
c
c    The Catalan number C(N) counts:
c
c    1) the number of binary trees on N vertices;
c    2) the number of ordered trees on N+1 vertices;
c    3) the number of full binary trees on 2N+1 vertices;
c    4) the number of well formed sequences of 2N parentheses;
c    5) number of ways 2N ballots can be counted, in order,
c       with N positive and N negative, so that the running sum
c       is never negative;
c    6) the number of standard tableaus in a 2 by N rectangular Ferrers diagram;
c    7) the number of monotone functions from [1..N} to [1..N} which
c       satisfy f(i) <= i for all i;
c    8) the number of ways to triangulate a polygon with N+2 vertices.
c
c  Example:
c
c    N = 3
c
c    ()()()
c    ()(())
c    (()())
c    (())()
c    ((()))
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    02 January 2006
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Dennis Stanton, Dennis White,
c    Constructive Combinatorics,
c    Springer, 1986,
c    ISBN: 0387963472.
c
c  Parameters:
c
c    Input, integer N, the number of Catalan numbers desired.
c
c    Output, integer C(N+1), the Catalan numbers from C(0) to C(N).
c
      implicit none

      integer n

      integer c(n+1)
      integer i

      if ( n .lt. 0 ) then
        return
      end if

      c(1) = 1
c
c  The extra parentheses ensure that the integer division is
c  done AFTER the integer multiplication.
c
      do i = 1, n
        c(i+1) = ( c(i) * 2 * ( 2 * i - 1 ) ) / ( i + 1 )
      end do

      return
      end
      subroutine catalan_row_next ( ido, n, irow )

c*********************************************************************72
c
cc CATALAN_ROW_NEXT computes row N of Catalan's triangle.
c
c  Example:
c
c    I\J 0   1   2   3   4   5   6
c
c    0   1
c    1   1   1
c    2   1   2   2
c    3   1   3   5   5
c    4   1   4   9  14  14
c    5   1   5  14  28  42  42
c    6   1   6  20  48  90 132 132
c
c  Recursion:
c
c    C(0,0) = 1
c    C(I,0) = 1
c    C(I,J) = 0 for I < J
c    C(I,J) = C(I,J-1) + C(I-1,J)
c    C(I,I) is the I-th Catalan number.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    03 January 2007
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer IDO, indicates whether this is a call for
c    the 'next' row of the triangle.
c    IDO = 0, this is a startup call.  Row N is desired, but
c    presumably this is a first call, or row N-1 was not computed
c    on the previous call.
c    IDO = 1, this is not the first call, and row N-1 was computed
c    on the previous call.  In this case, much work can be saved
c    by using the information from the previous values of IROW
c    to build the next values.
c
c    Input, integer N, the index of the row of the triangle desired.  
c
c    Input/output, integer IROW(0:N), the row of coefficients.
c    If IDO = 0, then IROW is not required to be set on input.
c    If IDO = 1, then IROW must be set on input to the value of
c    row N-1.
c
      implicit none

      integer n

      integer i
      integer ido
      integer irow(0:n)
      integer j

      if ( n .lt. 0 ) then
        return
      end if

      if ( ido .eq. 0 ) then
 
        irow(0) = 1

        do i = 1, n
          irow(i) = 0
        end do

        do i = 1, n

          irow(0) = 1

          do j = 1, i-1
            irow(j) = irow(j) + irow(j-1)
          end do

          irow(i) = irow(i-1)

        end do
 
      else
 
        irow(0) = 1

        do j = 1, n-1
          irow(j) = irow(j) + irow(j-1)
        end do

        if ( 1 .le. n ) then
          irow(n) = irow(n-1)
        end if
 
      end if
 
      return
      end
      subroutine catalan_values ( n_data, n, c )

c*********************************************************************72
c
cc CATALAN_VALUES returns some values of the Catalan numbers for testing.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    02 January 2007
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Milton Abramowitz, Irene Stegun,
c    Handbook of Mathematical Functions,
c    US Department of Commerce, 1964,
c    ISBN: 0-486-61272-4,
c    LC: QA47.A34.
c
c  Parameters:
c
c    Input/output, integer N_DATA.
c    On input, if N_DATA is 0, the first test data is returned, and N_DATA
c    is set to 1.  On each subsequent call, the input value of N_DATA is
c    incremented and that test data item is returned, if available.  When 
c    there is no more test data, N_DATA is set to 0.
c
c    Output, integer N, the order of the Catalan number.
c
c    Output, integer C, the value of the Catalan number.
c
      implicit none

      integer nmax
      parameter ( nmax = 11 )

      integer c
      integer c_vec(nmax)
      integer n
      integer n_data
      integer n_vec(nmax)

      save c_vec
      save n_vec

      data c_vec /
     &  1, 1, 2, 5, 14, 42, 132, 429, 1430, 4862, 16796 /

      data n_vec /
     &   0,  1,  2,  3,  4, 5,  6,  7,  8,  9,  10 /

      if ( n_data .lt. 0 ) then
        n_data = 0
      end if

      n_data = n_data + 1

      if ( nmax .lt. n_data ) then
        n_data = 0
        n = 0
        c = 0
      else
        n = n_vec(n_data)
        c = c_vec(n_data)
      end if

      return
      end
      subroutine cfrac_to_rat ( n, a, p, q )

c*********************************************************************72
c
cc CFRAC_TO_RAT converts a monic continued fraction to an ordinary fraction.
c
c  Discussion:
c
c    The routine is given the monic or "simple" continued fraction with
c    integer coefficients:
c
c      A(1) + 1 / ( A(2) + 1 / ( A(3) ... + 1 / A(N) ) )
c
c    and returns the N successive approximants P(I)/Q(I)
c    to the value of the rational number represented by the continued
c    fraction, with the value exactly equal to the final ratio P(N)/Q(N).
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    03 January 2007
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    John Hart, Ward Cheney, Charles Lawson, Hans Maehly, Charles Mesztenyi, 
c    John Rice, Henry Thatcher, Christoph Witzgall,
c    Computer Approximations,
c    Wiley, 1968.
c
c  Parameters:
c
c    Input, integer N, the number of continued fraction coefficients.
c
c    Input, integer A(N), the continued fraction coefficients.
c
c    Output, integer P(N), Q(N), the N successive approximations
c    to the value of the continued fraction.
c
      implicit none

      integer n
   
      integer a(n)
      integer i
      integer p(n)
      integer q(n)

      do i = 1, n

        if ( i .eq. 1 ) then
          p(i) = a(i) * 1 + 0
          q(i) = a(i) * 0 + 1
        else if ( i .eq. 2 ) then
          p(i) = a(i) * p(i-1) + 1
          q(i) = a(i) * q(i-1) + 0
        else
          p(i) = a(i) * p(i-1) + p(i-2)
          q(i) = a(i) * q(i-1) + q(i-2)
        end if

      end do

      return
      end
      subroutine cfrac_to_rfrac ( m, g, h, p, q )

c*********************************************************************72
c
cc CFRAC_TO_RFRAC converts polynomial fractions from continued to rational form.
c
c  Discussion:
c
c    The routine accepts a continued polynomial fraction:
c
c      G(1)     / ( H(1) +
c      G(2) * X / ( H(2) +
c      G(3) * X / ( H(3) + ...
c      G(M) * X / ( H(M) )...) ) )
c
c    and returns the equivalent rational polynomial fraction:
c
c      P(1) + P(2) * X + ... + P(L1) * X**(L1)
c      -------------------------------------------------------
c      Q(1) + Q(2) * X + ... + Q(L2) * X**(L2-1)
c
c    where
c
c      L1 = (M+1)/2
c      L2 = (M+2)/2.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    04 January 2007
c
c  Author:
c
c    Original FORTRAN77 version by John Hart, Ward Cheney, Charles Lawson, 
c    Hans Maehly, Charles Mesztenyi, John Rice, Henry Thatcher, Christoph Witzgall.
c    This FORTRAN77 version by John Burkardt.
c
c  Reference:
c
c    John Hart, Ward Cheney, Charles Lawson, Hans Maehly, Charles Mesztenyi,
c    John Rice, Henry Thatcher, Christoph Witzgall,
c    Computer Approximations,
c    Wiley, 1968.
c
c  Parameters:
c
c    Input, integer M, the number of continued fraction polynomial coefficients.
c
c    Input, double precision G(M), H(M), the continued polynomial 
c    fraction coefficients.
c
c    Output, double precision P((M+1)/2), Q((M+2)/2), the rational 
c    polynomial fraction coefficients.
c
      implicit none

      integer m

      double precision a(m,(m+2)/2)
      double precision g(m)
      double precision h(m)
      integer i
      integer j
      double precision p((m+1)/2)
      double precision q((m+2)/2)

      if ( m .eq. 1 ) then
        p(1) = g(1)
        q(1) = h(1)
        return
      end if

      do i = 1, m
        do j = 1, (m+2)/2
          a(i,j) = 0.0D+00
        end do
      end do
c
c  Solve for P's.
c
      a(1,1) = g(1)
      a(2,1) = g(1) * h(2)

      do i = 3, m
        a(i,1) = h(i) * a(i-1,1)
        do j = 2, (i+1)/2
          a(i,j) = h(i) * a(i-1,j) + g(i) * a(i-2,j-1)
        end do
      end do

      do j = 1, (m+1)/2
        p(j) = a(m,j)
      end do
c
c  Solve for Q's.
c
      a(1,1) = h(1)
      a(2,1) = h(1) * h(2)
      a(2,2) = g(2)

      do i = 3, m
        a(i,1) = h(i) * a(i-1,1)
        do j = 2, (i+2) / 2
          a(i,j) = h(i) * a(i-1,j) + g(i) * a(i-2,j-1)
        end do
      end do

      do j = 1, (m+2)/2
        q(j) = a(m,j)
      end do

      return
      end
      subroutine ch_cap ( ch )

c*********************************************************************72
c
cc CH_CAP capitalizes a single character.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    03 January 2007
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input/output, character CH, the character to capitalize.
c
      implicit none

      character ch
      integer itemp

      itemp = ichar ( ch )
 
      if ( 97 .le. itemp .and. itemp .le. 122 ) then
        ch = char ( itemp - 32 )
      end if
 
      return
      end
      subroutine change_greedy ( total, coin_num, coin_value, 
     &  change_num, change )

c*********************************************************************72
c
cc CHANGE_GREEDY makes change for a given total using the biggest coins first.
c
c  Discussion:
c
c    The algorithm is simply to use as many of the largest coin first,
c    then the next largest, and so on.
c
c    It is assumed that there is always a coin of value 1.  The
c    algorithm will otherwise fail!
c
c  Example:
c
c    Total = 17
c    COIN_NUM = 3
c    COIN_VALUE = (/ 1, 5, 10 /)
c
c
c    #  CHANGE              COIN_VALUE(CHANGE)
c
c    4  3 2 1 1             10 5 1 1
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    04 January 2007
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer TOTAL, the total for which change is to be made.
c
c    Input, integer COIN_NUM, the number of types of coins.
c
c    Input, integer COIN_VALUE(COIN_NUM), the value of each coin.
c    The values should be in ascending order, and if they are not,
c    they will be sorted.
c
c    Output, integer CHANGE_NUM, the number of coins given in change.
c
c    Output, integer CHANGE(TOTAL), the indices of the coins will be
c    in entries 1 through CHANGE_NUM.
c
      implicit none

      integer coin_num
      integer total

      integer change(total)
      integer change_num
      integer coin_value(coin_num)
      integer j
      integer total_copy

      change_num = 0
c
c  Find the largest coin smaller than the total.
c
      j = coin_num

10    continue

      if ( 0 .lt. j ) then
        if ( coin_value(j) .le. total ) then
          go to 20
        end if
        j = j - 1
        go to 10
      end if

20    continue

      if ( j .le. 0 ) then
        return
      end if
c
c  Subtract the current coin from the total.
c  Once that coin is too big, use the next coin.
c
      total_copy = total

30    continue

      if ( 0 .lt. total_copy ) then

        if ( coin_value(j) .le. total_copy ) then

          total_copy = total_copy - coin_value(j)
          change_num = change_num + 1
          change(change_num) = j

        else

          j = j - 1
          if ( j .le. 0 ) then
            go to 40
          end if

        end if

        go to 30

      end if

40    continue

      return
      end
      subroutine change_next ( total, coin_num, coin_value, change_num, 
     & change, done  )

c*********************************************************************72
c
cc CHANGE_NEXT computes the next set of change for a given sum.
c
c  Example:
c
c    Total = 17
c    COIN_NUM = 3
c    COIN_VALUE = (/ 1, 5, 10 /)
c
c
c        #  CHANGE              COIN_VALUE(CHANGE)
c
c    1   4  3 2 1 1             10 5 1 1
c    2   8  3 1 1 1 1 1 1 1     10 1 1 1 1 1 1 1
c    3   5  2 2 2 1 1            5 5 5 1 1
c    4   9  2 2 1 1 1 1 1 1 1    5 5 1 1 1 1 1 1 1
c    5  13  2 1 1 1 1 1 1 1 1 1  5 1 1 1 1 1 1 1 1 1
c           1 1 1                1 1 1
c    6  17  1 1 1 1 1 1 1 1 1 1  1 1 1 1 1 1 1 1 1 1 1
c           1 1 1 1 1 1 1        1 1 1 1 1 1
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    06 January 2007
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer TOTAL, the total for which change is to be made.
c
c    Input, integer COIN_NUM, the number of types of coins.
c
c    Input, integer COIN_VALUE(COIN_NUM), the value of each coin.
c    The values must be in ascending order.
c
c    Input/output, integer CHANGE_NUM, the number of coins given in change
c    for this form of the change.
c
c    Input/output, integer CHANGE(CHANGE_NUM), the indices of the coins.
c    The user must dimension this array to have dimension TOTAL!
c
c    Input/output, logical DONE.  The user sets DONE = TRUE on
c    first call to tell the routine this is the beginning of a computation.
c    The program resets DONE to FALSE and it stays that way until
c    the last possible change combination is made, at which point the
c    program sets DONE to TRUE again.
c
      implicit none
 
      integer coin_num
      integer total

      integer change(total)
      integer change_num
      integer change_num2
      integer coin_num2
      integer coin_value(coin_num)
      logical done
      integer i
      logical i4vec_ascends
      integer last
      integer total2

      if ( done ) then
c
c  Make sure the coin values are sorted.
c
        if ( .not. i4vec_ascends ( coin_num, coin_value ) ) then
          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) 'CHANGE_NEXT - Fatal error!'
          write ( *, '(a)' ) 
     &      '  The COIN_VALUE array is not in ascending order.'
          stop
        end if
c
c  Start with the greedy change.
c
        call change_greedy ( total, coin_num, coin_value, change_num, 
     &    change )
c
c  In a few cases, like change for 4 cents, we're done after the first call.
c
        if ( change_num .eq. total ) then
          done = .true.
        else
          done = .false.
        end if

      else
c
c  Find the last location in the input change which is NOT a penny.
c
        last = 0

        do i = change_num, 1, -1

          if ( change(i) .ne. 1 ) then
            last = i
            go to 10
          end if

        end do

10      continue
c
c  If that location is still 0, an error was made.
c
        if ( last .eq. 0 ) then
          done = .true.
          return
        end if
c
c  Sum the entries from that point to the end.
c
        total2 = 0

        do i = last, change_num
          total2 = total2 + coin_value ( change(i) )
        end do
c
c  Make greedy change for the partial sum using coins smaller than that one.
c
        coin_num2 = change(last) - 1

        call change_greedy ( total2, coin_num2, coin_value, 
     &    change_num2, change(last) )

        change_num = ( last - 1 ) + change_num2

      end if

      return
      end
      subroutine chinese_check ( n, m, ierror )

c*********************************************************************72
c
cc CHINESE_CHECK checks the Chinese remainder moduluses.
c
c  Discussion:
c
c    For a Chinese remainder representation, the moduluses M(I) must
c    be positive and pairwise prime.  Also, in case this is not obvious,
c    no more than one of the moduluses may be 1.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    06 January 2007
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the number of moduluses.
c
c    Input, integer M(N), the moduluses.  These should be positive
c    and pairwise prime.
c
c    Output, integer IERROR, an error flag.
c    0, no error was detected.
c    nonzero, an error was detected.
c
      implicit none

      integer n

      integer i
      integer ierror
      logical i4vec_pairwise_prime
      integer j
      integer m(n)

      ierror = 0
c
c  Do not allow nonpositive entries.
c
      do i = 1, n
        if ( m(i) .le. 0 ) then
          ierror = 1
          return
        end if
      end do
c
c  Allow one entry to be 1, but not two entries.
c
      do i = 1, n
        do j = i+1, n
          if ( m(i) .eq. 1 .and. m(j) .eq. 1 ) then
            ierror = 2
            return
          end if
        end do
      end do
c
c  Now check pairwise primeness.
c
      if ( .not. i4vec_pairwise_prime ( n, m ) ) then
        ierror = 3
        return
      end if

      return
      end
      subroutine chinese_to_i4 ( n, m, r, j )

c*********************************************************************72
c
cc CHINESE_TO_I4 converts a set of Chinese remainders to an equivalent integer.
c
c  Discussion:
c
c    Given a set of N pairwise prime, positive moduluses M(I), and
c    a corresponding set of remainders R(I), this routine finds an
c    integer J such that, for all I,
c
c      J = R(I) mod M(I)
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    07 January 2007
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the number of moduluses.
c
c    Input, integer M(N), the moduluses.  These should be positive
c    and pairwise prime.
c
c    Input, integer R(N), the Chinese remainder representation of the integer.
c
c    Output, integer J, the corresponding integer.
c
      implicit none

      integer n

      integer a
      integer b(n)
      integer big_m
      integer c
      integer i
      integer ierror
      integer j
      integer m(n)
      integer r(n)

      call chinese_check ( n, m, ierror )

      if ( ierror .ne. 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'CHINESE_TO_I4 - Fatal error!'
        write ( *, '(a)' ) '  The moduluses are not legal.'
        stop
      end if
c
c  Set BIG_M.
c
      big_m = 1
      do i = 1, n
        big_m = big_m * m(i)
      end do
c
c  Solve BIG_M / M(I) * B(I) = 1, mod M(I)
c
      do i = 1, n
        a = big_m / m(i)
        c = 1
        call congruence ( a, m(i), c, ierror, b(i) )
      end do
c
c  Set J = sum ( 1 <= I <= N ) ( R(I) * B(I) * BIG_M / M(I) ) mod M
c
      j = 0
      do i = 1, n
        j = mod ( j + r(i) * b(i) * ( big_m / m(i) ), big_m )
      end do

      return
      end
      subroutine comb_next ( n, k, a, done )

c*********************************************************************72
c
cc COMB_NEXT computes combinations of K things out of N.
c
c  Discussion:
c
c    The combinations are computed one at a time, in lexicographical order.
c
c    10 April 1009: Thanks to "edA-qa mort-ora-y" for supplying a 
c    correction to this code!
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    10 April 2009
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Charles Mifsud,
c    Algorithm 154:
c    Combination in Lexicographic Order,
c    Communications of the ACM, 
c    March 1963.
c
c  Parameters:
c
c    Input, integer N, the total number of things.
c
c    Input, integer K, the number of things in each combination.
c
c    Input/output, integer A(K), contains the list of elements in
c    the current combination.
c
c    Input/output, logical DONE.  On first call, set DONE to TRUE,
c    and therafter, its input value should be the output value from
c    the previous call.  The output value will normally be FALSE,
c    indicating that there are further combinations that can be
c    returned.  When DONE is returned TRUE, the sequence is exhausted.
c
      implicit none

      integer k

      integer a(k)
      logical done
      integer i
      integer j
      integer n

      if ( done ) then

        if ( k <= 0 ) then
          return
        end if 

        call i4vec_indicator ( k, a )

        done = .false.

      else

        if ( a(k) .lt. n ) then
          a(k) = a(k) + 1
          return
        end if

        do i = k, 2, -1

          if ( a(i-1) .lt. n-k+i-1 ) then

            a(i-1) = a(i-1) + 1

            do j = i, k
              a(j) = a(i-1) + j - ( i-1 )
            end do

            return

          end if

        end do

        done = .true.

      end if

      return
      end
      subroutine comb_row ( ido, n, irow )

c*********************************************************************72
c
cc COMB_ROW computes row N of Pascal's triangle.
c
c  Discussion:
c
c    Row N contains the combinatorial coefficients
c
c      C(N,0), C(N,1), C(N,2), ... C(N,N)
c
c    The sum of the elements of row N is equal to 2**N.
c
c  Formula:
c
c    C(N,K) = N! / ( K! * (N-K)! )
c
c  First terms:
c
c     N K:0  1   2   3   4   5   6   7  8  9 10
c
c     0   1
c     1   1  1
c     2   1  2   1
c     3   1  3   3   1
c     4   1  4   6   4   1
c     5   1  5  10  10   5   1
c     6   1  6  15  20  15   6   1
c     7   1  7  21  35  35  21   7   1
c     8   1  8  28  56  70  56  28   8  1
c     9   1  9  36  84 126 126  84  36  9  1
c    10   1 10  45 120 210 252 210 120 45 10  1
c
c  Recursion:
c
c    C(N,K) = C(N-1,K-1)+C(N-1,K)
c
c  Special values:
c
c    C(N,0) = C(N,N) = 1
c    C(N,1) = C(N,N-1) = N
c    C(N,N-2) = sum ( 1 <= I <= N ) N
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    07 January 2007
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer IDO, indicates whether this is a call for
c    the 'next' row of the triangle.
c
c    0 means this is a startup call.  Row N is desired, but
c    presumably this is a first call, or row N-1 was not computed
c    on the previous call.
c
c    1 means this is not the first call, and row N-1 was computed
c    on the previous call.  In this case, much work can be saved
c    by using the information from the previous values of IROW
c    to build the next values.
c
c    Input, integer N, the row of the triangle desired.  The triangle
c    begins with row 0.
c
c    Output, integer IROW(N+1), the row of coefficients.
c    IROW(I) = C(N,I-1).
c
      implicit none

      integer n

      integer i
      integer ido
      integer irow(n+1)
      integer j

      if ( n .lt. 0 ) then
        return
      end if

      if ( ido .eq. 1 ) then

        do i = n, 2, -1
          irow(i) = irow(i) + irow(i-1)
        end do

        irow(n+1) = 1

      else

        irow(1) = 1
        do i = 2, n+1
          irow(i) = 0
        end do

        do j = 1, n
          do i = j+1, 2, -1
            irow(i) = irow(i) + irow(i-1)
          end do
        end do

      end if

      return
      end
      subroutine comb_unrank ( m, n, rank, a )

c*********************************************************************72
c
cc COMB_UNRANK returns the RANK-th combination of N things out of M.
c
c  Discussion:
c
c    Going from a rank to a thing is called "unranking".
c
c    The combinations are ordered lexically.
c
c    Lexical order can be illustrated for the general case of N and M as
c    follows:
c
c    1:       1,     2,     3,     ..., N-2, N-1, N
c    2:       1,     2,     3,     ..., N-2, N-1, N+1
c    3:       1,     2,     3,     ..., N-2, N-1, N+2
c    ...
c    M-N+1:   1,     2,     3,     ..., N-2, N-1, M
c    M-N+2:   1,     2,     3,     ..., N-2, N,   N+1
c    M-N+3:   1,     2,     3,     ..., N-2, N,   N+2
c    ...
c    LAST-2:  M-N,   M-N+1, M-N+3, ..., M-2, M-1, M
c    LAST-1:  M-N,   M-N+2, M-N+3, ..., M-2, M-1, M
c    LAST:    M-N+1, M-N+2, M-N+3, ..., M-2, M-1, M
c
c    There are a total of M!/(N!*(M-N)!) combinations of M
c    things taken N at a time.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    07 January 2007
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Bill Buckles, Matthew Lybanon,
c    Algorithm 515,
c    Generation of a Vector from the Lexicographical Index,
c    ACM Transactions on Mathematical Software,
c    Volume 3, Number 2, pages 180-182, June 1977.
c
c  Parameters:
c
c    Input, integer M, the size of the set.
c
c    Input, integer N, the number of things in the combination.
c    N must be greater than 0, and no greater than M.
c
c    Input, integer RANK, the lexicographical rank of the combination
c    sought.  RANK must be at least 1, and no greater than M!/(N!*(M-N)!).
c
c    Output, integer A(N), the combination.
c
      implicit none

      integer n

      integer a(n)
      integer i
      integer j
      integer k
      integer m
      integer rank
c
c  Initialize the lower bound index.
c
      k = 0
c
c  Select elements in ascending order.
c
      do i = 1, n - 1
c
c  Set the lower bound element number for next element value.
c
        a(i) = 0

        if ( 1 .lt. i ) then
          a(i) = a(i-1)
        end if
c
c  Check each element value.
c
10      continue

          a(i) = a(i) + 1
          call combin2 ( m-a(i), n-i, j )
          k = k + j

          if ( rank .le. k ) then
            go to 20
          end if

        go to 10

20      continue

        k = k - j

      end do

      a(n) = a(n-1) + rank - k

      return
      end
      subroutine combin ( n, k, cnk )

c*********************************************************************72
c
cc COMBIN computes the combinatorial coefficient C(N,K).
c
c  Discussion:
c
c    C(N,K) is the number of distinct combinations of K objects
c    chosen from a set of N distinct objects.  A combination is
c    like a set, in that order does not matter.
c
c    Real arithmetic is used in this routine, and C(N,K) is computed
c    directly, via Gamma functions, rather than recursively.
c
c  Example:
c
c    The number of combinations of 2 things chosen from 5 is 10.
c
c    C(5,2) = ( 5 * 4 * 3 * 2 * 1 ) / ( ( 3 * 2 * 1 ) * ( 2 * 1 ) ) = 10.
c
c    The actual combinations may be represented as:
c
c      (1,2), (1,3), (1,4), (1,5), (2,3),
c      (2,4), (2,5), (3,4), (3,5), (4,5).
c
c  Formula:
c
c    C(N,K) = N! / ( (N-K)! * K! )
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    07 January 2007
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the value of N.
c
c    Input, integer K, the value of K.
c
c    Output, double precision CNK, the value of C(N,K)
c
      implicit none

      double precision arg
      double precision cnk
      double precision fack
      double precision facn
      double precision facnmk
      double precision gamma_log
      integer k
      integer n

      if ( n .lt. 0 ) then

        cnk = 0.0D+00

      else if ( k .eq. 0 ) then

        cnk = 1.0D+00

      else if ( k .eq. 1 ) then

        cnk = dble ( n )

      else if ( 1 .lt. k .and. k .lt. n-1 ) then

        arg = dble ( n + 1 )
        facn = gamma_log ( arg )

        arg = dble ( k + 1 )
        fack = gamma_log ( arg )

        arg = dble ( n - k + 1 )
        facnmk = gamma_log ( arg )

        cnk = anint ( exp ( facn - fack - facnmk ) )

      else if ( k .eq. n - 1 ) then

        cnk = dble ( n )

      else if ( k .eq. n ) then

        cnk = 1.0D+00

      else

        cnk = 0.0D+00

      end if

      return
      end
      subroutine combin2 ( n, k, icnk )

c*********************************************************************72
c
cc COMBIN2 computes the binomial coefficient C(N,K).
c
c  Discussion:
c
c    The value is calculated in such a way as to avoid overflow and
c    roundoff.  The calculation is done in integer arithmetic.
c
c    The formula used is:
c
c      C(N,K) = N! / ( K! * (N-K)! )
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    03 January 2007
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    ML Wolfson, HV Wright,
c    Algorithm 160:
c    Combinatorial of M Things Taken N at a Time,
c    Communications of the ACM,
c    April, 1963.
c
c  Parameters:
c
c    Input, integer N, K, are the values of N and K.
c
c    Output, integer ICNK, the number of combinations of N
c    things taken K at a time.
c
      implicit none

      integer i
      integer icnk
      integer k
      integer mn
      integer mx
      integer n

      mn = min ( k, n - k )

      if ( mn .lt. 0 ) then

        icnk = 0

      else if ( mn .eq. 0 ) then

        icnk = 1

      else

        mx = max ( k, n - k )
        icnk = mx + 1

        do i = 2, mn
          icnk = ( icnk * ( mx + i ) ) / i
        end do

      end if

      return
      end
      subroutine comp_enum ( n, k, number )

c*********************************************************************72
c
cc COMP_ENUM returns the number of compositions of the integer N into K parts.
c
c  Discussion:
c
c    A composition of the integer N into K parts is an ordered sequence
c    of K nonnegative integers which sum to N.  The compositions (1,2,1)
c    and (1,1,2) are considered to be distinct.
c
c    The 28 compositions of 6 into three parts are:
c
c      6 0 0,  5 1 0,  5 0 1,  4 2 0,  4 1 1,  4 0 2,
c      3 3 0,  3 2 1,  3 1 2,  3 0 3,  2 4 0,  2 3 1,
c      2 2 2,  2 1 3,  2 0 4,  1 5 0,  1 4 1,  1 3 2,
c      1 2 3,  1 1 4,  1 0 5,  0 6 0,  0 5 1,  0 4 2,
c      0 3 3,  0 2 4,  0 1 5,  0 0 6.
c
c    The formula for the number of compositions of N into K parts is
c
c      Number = ( N + K - 1 )! / ( N! * ( K - 1 )! )
c
c    (Describe the composition using N '1's and K-1 dividing lines '|'.
c    The number of distinct permutations of these symbols is the number
c    of compositions.  This is equal to the number of permutations of 
c    N+K-1 things, with N identical of one kind and K-1 identical of another.)
c
c    Thus, for the above example, we have:
c
c      Number = ( 6 + 3 - 1 )! / ( 6! * (3-1)! ) = 28
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    16 January 2007
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Albert Nijenhuis, Herbert Wilf,
c    Combinatorial Algorithms for Computers and Calculators,
c    Second Edition,
c    Academic Press, 1978,
c    ISBN: 0-12-519260-6,
c    LC: QA164.N54.
c
c  Parameters:
c
c    Input, integer N, the integer whose compositions are desired.
c
c    Input, integer K, the number of parts in the composition.
c
c    Output, integer NUMBER, the number of compositions of N into K parts.
c
      implicit none

      integer k
      integer n
      integer number

      call combin2 ( n + k - 1, n, number )

      return
      end
      subroutine comp_next ( n, k, a, more, h, t )

c*********************************************************************72
c
cc COMP_NEXT computes the compositions of the integer N into K parts.
c
c  Discussion:
c
c    A composition of the integer N into K parts is an ordered sequence
c    of K nonnegative integers which sum to N.  The compositions (1,2,1)
c    and (1,1,2) are considered to be distinct.
c
c    The routine computes one composition on each call until there are no more.
c    For instance, one composition of 6 into 3 parts is
c    3+2+1, another would be 6+0+0.
c
c    On the first call to this routine, set MORE = FALSE.  The routine
c    will compute the first element in the sequence of compositions, and
c    return it, as well as setting MORE = TRUE.  If more compositions
c    are desired, call again, and again.  Each time, the routine will
c    return with a new composition.
c
c    However, when the LAST composition in the sequence is computed 
c    and returned, the routine will reset MORE to FALSE, signaling that
c    the end of the sequence has been reached.
c
c    This routine originally used a SAVE statement to maintain the
c    variables H and T.  I have decided (based on an wasting an
c    entire morning trying to track down a problem) that it is safer
c    to pass these variables as arguments, even though the user should
c    never alter them.  This allows this routine to safely shuffle
c    between several ongoing calculations.
c
c
c    There are 28 compositions of 6 into three parts.  This routine will
c    produce those compositions in the following order:
c
c     I         A
c     -     ---------
c     1     6   0   0
c     2     5   1   0
c     3     4   2   0
c     4     3   3   0
c     5     2   4   0
c     6     1   5   0
c     7     0   6   0
c     8     5   0   1
c     9     4   1   1
c    10     3   2   1
c    11     2   3   1
c    12     1   4   1
c    13     0   5   1
c    14     4   0   2
c    15     3   1   2
c    16     2   2   2
c    17     1   3   2
c    18     0   4   2
c    19     3   0   3
c    20     2   1   3
c    21     1   2   3
c    22     0   3   3
c    23     2   0   4
c    24     1   1   4
c    25     0   2   4
c    26     1   0   5
c    27     0   1   5
c    28     0   0   6
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    02 July 2008
c
c  Author:
c
c    Original FORTRAN77 version by Albert Nijenhuis, Herbert Wilf.
c    This FORTRAN77 version by John Burkardt.
c
c  Reference:
c
c    Albert Nijenhuis, Herbert Wilf,
c    Combinatorial Algorithms for Computers and Calculators,
c    Second Edition,
c    Academic Press, 1978,
c    ISBN: 0-12-519260-6,
c    LC: QA164.N54.
c
c  Parameters:
c
c    Input, integer N, the integer whose compositions are desired.
c
c    Input, integer K, the number of parts in the composition.
c
c    Input/output, integer A(K), the parts of the composition.
c
c    Input/output, logical MORE, set by the user to start the computation,
c    and by the routine to terminate it.
c
c    Input/output, integer H, T, two internal parameters needed for the
c    computation.  The user should allocate space for these in the calling
c    program, include them in the calling sequence, but never alter them!
c
      implicit none

      integer k

      integer a(k)
      integer h
      integer i
      logical more
      integer n
      integer t
c
c  The first computation.
c
      if ( .not. more ) then

        t = n
        h = 0
        a(1) = n
        do i = 2, k
          a(i) = 0
        end do
c
c  The next computation.
c
      else

        if ( 1 .lt. t ) then
          h = 0
        end if

        h = h + 1
        t = a(h)
        a(h) = 0
        a(1) = t - 1
        a(h+1) = a(h+1) + 1

      end if
c
c  This is the last element of the sequence if all the
c  items are in the last slot.
c
      more = ( a(k) .ne. n )

      return
      end
      subroutine comp_random ( n, k, seed, a )

c*********************************************************************72
c
cc COMP_RANDOM selects a random composition of the integer N into K parts.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    07 January 2007
c
c  Author:
c
c    Original FORTRAN77 version by Albert Nijenhuis, Herbert Wilf.
c    This FORTRAN77 version by John Burkardt.
c
c  Reference:
c
c    Albert Nijenhuis, Herbert Wilf,
c    Combinatorial Algorithms for Computers and Calculators,
c    Second Edition,
c    Academic Press, 1978,
c    ISBN: 0-12-519260-6,
c    LC: QA164.N54.
c
c  Parameters:
c
c    Input, integer N, the integer to be decomposed.
c
c    Input, integer K, the number of parts in the composition.
c
c    Input/output, integer SEED, a seed for the random number generator.
c
c    Output, integer A(K), the parts of the composition.
c
      implicit none

      integer k

      integer a(k)
      integer i
      integer l
      integer m
      integer n
      integer seed

      call ksub_random ( n+k-1, k-1, seed, a )

      a(k) = n + k
      l = 0

      do i = 1, k
        m = a(i)
        a(i) = a(i) - l - 1
        l = m
      end do

      return
      end
      subroutine compnz_enum ( n, k, number )

c*********************************************************************72
c
cc COMPNZ_ENUM returns the number of nonzero compositions of the N into K parts.
c
c  Discussion:
c
c    A composition of the integer N into K nonzero parts is an ordered sequence
c    of K positive integers which sum to N.  The compositions (1,2,1)
c    and (1,1,2) are considered to be distinct.
c
c    The routine computes one composition on each call until there are no more.
c    For instance, one composition of 6 into 3 parts is 3+2+1, another would
c    be 4+1+1 but 5+1+0 is not allowed since it includes a zero part.
c
c    On the first call to this routine, set MORE = FALSE.  The routine
c    will compute the first element in the sequence of compositions, and
c    return it, as well as setting MORE = TRUE.  If more compositions
c    are desired, call again, and again.  Each time, the routine will
c    return with a new composition.
c
c    However, when the LAST composition in the sequence is computed 
c    and returned, the routine will reset MORE to FALSE, signaling that
c    the end of the sequence has been reached.
c
c    The 10 compositions of 6 into three nonzero parts are:
c
c      4 1 1,  3 2 1,  3 1 2,  2 3 1,  2 2 2,  2 1 3,  
c      1 4 1,  1 3 2,  1 2 3,  1 1 4.
c
c    The formula for the number of compositions of N into K nonzero
c    parts is
c
c      Number = ( N - 1 )! / ( ( N - K )! * ( K - 1 )! )
c
c    (Describe the composition using N-K '1's and K-1 dividing lines '|'.
c    The number of distinct permutations of these symbols is the number
c    of compositions into nonzero parts.  This is equal to the number of 
c    permutations of  N-1 things, with N-K identical of one kind 
c    and K-1 identical of another.)
c
c    Thus, for the above example, we have:
c
c      Number = ( 6 - 1 )! / ( ( 6 - 3 )! * ( 3 - 1 )! ) = 10
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    16 January 2007
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Albert Nijenhuis, Herbert Wilf,
c    Combinatorial Algorithms for Computers and Calculators,
c    Second Edition,
c    Academic Press, 1978,
c    ISBN: 0-12-519260-6,
c    LC: QA164.N54.
c
c  Parameters:
c
c    Input, integer N, the integer whose compositions are desired.
c
c    Input, integer K, the number of parts in the composition.
c
c    Output, integer NUMBER, the number of compositions of N into 
c    K nonzero parts.
c
      implicit none

      integer k
      integer n
      integer number

      call combin2 ( n - 1, n - k, number )

      return
      end
      subroutine compnz_next ( n, k, a, more )

c*********************************************************************72
c
cc COMPNZ_NEXT computes the compositions of the integer N into K nonzero parts.
c
c  Discussion:
c
c    A composition of the integer N into K nonzero parts is an ordered sequence
c    of K positive integers which sum to N.  The compositions (1,2,1)
c    and (1,1,2) are considered to be distinct.
c
c    The routine computes one composition on each call until there are no more.
c    For instance, one composition of 6 into 3 parts is 3+2+1, another would
c    be 4+1+1 but 5+1+0 is not allowed since it includes a zero part.
c
c    On the first call to this routine, set MORE = FALSE.  The routine
c    will compute the first element in the sequence of compositions, and
c    return it, as well as setting MORE = TRUE.  If more compositions
c    are desired, call again, and again.  Each time, the routine will
c    return with a new composition.
c
c    However, when the LAST composition in the sequence is computed 
c    and returned, the routine will reset MORE to FALSE, signaling that
c    the end of the sequence has been reached.
c
c  Example:
c
c    The 10 compositions of 6 into three nonzero parts are:
c
c      4 1 1,  3 2 1,  3 1 2,  2 3 1,  2 2 2,  2 1 3,  
c      1 4 1,  1 3 2,  1 2 3,  1 1 4.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    07 January 2007
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Albert Nijenhuis, Herbert Wilf,
c    Combinatorial Algorithms for Computers and Calculators,
c    Second Edition,
c    Academic Press, 1978,
c    ISBN: 0-12-519260-6,
c    LC: QA164.N54.
c
c  Parameters:
c
c    Input, integer N, the integer whose compositions are desired.
c
c    Input, integer K, the number of parts in the composition.  K must
c    be less than or equal to N.
c
c    Input/output, integer A(K), the parts of the composition.
c
c    Input/output, logical MORE, set by the user to start the computation,
c    and by the routine to terminate it.
c
      implicit none

      integer k

      integer a(k)
      integer h
      integer i
      logical more
      integer n
      integer t

      save h
      save t

      data h / 0 /
      data t / 0 /
c
c  We use the trick of computing ordinary compositions of (N-K)
c  into K parts, and adding 1 to each part.
c
      if ( n .lt. k ) then
        more = .false.
        do i = 1, k
          a(i) = -1
        end do
        return
      end if
c
c  The first computation.
c
      if ( .not. more ) then

        t = n - k
        h = 0
        a(1) = n - k
        do i = 2, k
          a(i) = 0
        end do
c
c  The next computation.
c
      else

        do i = 1, k
          a(i) = a(i) - 1
        end do

        if ( 1 .lt. t ) then
          h = 0
        end if

        h = h + 1
        t = a(h)
        a(h) = 0
        a(1) = t - 1
        a(h+1) = a(h+1) + 1

      end if
c
c  This is the last element of the sequence if all the
c  items are in the last slot.
c
      more = ( a(k) .ne. ( n - k ) )

      do i = 1, k
        a(i) = a(i) + 1
      end do

      return
      end
      subroutine compnz_random ( n, k, seed, a )

c*********************************************************************72
c
cc COMPNZ_RANDOM selects a random composition of N into K nonzero parts.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    07 January 2007
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Albert Nijenhuis, Herbert Wilf,
c    Combinatorial Algorithms for Computers and Calculators,
c    Second Edition,
c    Academic Press, 1978,
c    ISBN: 0-12-519260-6,
c    LC: QA164.N54.
c
c  Parameters:
c
c    Input, integer N, the integer to be decomposed.
c
c    Input, integer K, the number of parts in the composition.  K must
c    be no greater than N.
c
c    Input/output, integer SEED, a seed for the random number generator.
c
c    Output, integer A(K), the parts of the composition.
c
      implicit none

      integer k

      integer a(k)
      integer i
      integer l
      integer m
      integer n
      integer seed

      if ( n .lt. k ) then
        do i = 1, k
          a(i) = -1
        end do
        return
      end if

      call ksub_random ( n-1, k-1, seed, a )

      a(k) = n
      l = 0

      do i = 1, k
        m = a(i)
        a(i) = a(i) - l - 1
        l = m
      end do

      do i = 1, k
        a(i) = a(i) + 1
      end do

      return
      end
      subroutine congruence ( a, b, c, ierror, x )

c*********************************************************************72
c
cc CONGRUENCE solves a congruence of the form A * X = C ( mod B ).
c
c  Discussion:
c
c    A, B and C are given integers.  The equation is solvable if and only
c    if the greatest common divisor of A and B also divides C.
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
c  Reference:
c
c    Eric Weisstein, editor,
c    CRC Concise Encylopedia of Mathematics,
c    CRC Press, 1998, page 446.
c
c  Parameters:
c
c    Input, integer A, B, C, the coefficients of the Diophantine equation.
c
c    Output, integer IERROR, error flag.
c    0, no error, X was computed.
c    1, A = B = 0, C is nonzero.
c    2, A = 0, B and C nonzero, but C is not a multiple of B.
c    3, A nonzero, B zero, C nonzero, but C is not a multiple of A.
c    4, A, B, C nonzero, but GCD of A and B does not divide C.
c    5, algorithm ran out of internal space.
c
c    Output, integer X, the solution of the Diophantine equation.
c    X will be between 0 and B-1.
c
      implicit none

      integer nmax
      parameter ( nmax = 100 )

      integer a
      integer a_copy
      integer a_mag
      integer a_sign
      integer b
      integer b_copy
      integer b_mag
      integer b_sign
      integer c
      integer c_copy
      integer g
      integer i4_gcd
      integer ierror
      integer k
      integer n
      integer q(nmax)
      logical swap
      integer x
      integer y
      integer z
c
c  Defaults for output parameters.
c
      ierror = 0
      x = 0
      y = 0
c
c  Special cases.
c
      if ( a .eq. 0 .and. b .eq. 0 .and. c .eq. 0 ) then
        x = 0
        return
      else if ( a .eq. 0 .and. b .eq. 0 .and. c .ne. 0 ) then
        ierror = 1
        x = 0
        return
      else if ( a .eq. 0 .and. b .ne. 0 .and. c .eq. 0 ) then
        x = 0
        return
      else if ( a .eq. 0 .and. b .ne. 0 .and. c .ne. 0 ) then
        x = 0
        if ( mod ( c, b ) .ne. 0 ) then
          ierror = 2
        end if
        return
      else if ( a .ne. 0 .and. b .eq. 0 .and. c .eq. 0 ) then
        x = 0
        return
      else if ( a .ne. 0 .and. b .eq. 0 .and. c .ne. 0 ) then
        x = c / a
        if ( mod ( c, a ) .ne. 0 ) then
          ierror = 3
        end if
        return
      else if ( a .ne. 0 .and. b .ne. 0 .and. c .eq. 0 ) then
c       g = i4_gcd ( a, b )
c       x = b / g
        x = 0
        return
      end if
c
c  Handle the "general" case: A, B and C are nonzero.
c
c  Step 1: Compute the GCD of A and B, which must also divide C.
c
      g = i4_gcd ( a, b )

      if ( mod ( c, g ) .ne. 0 ) then
        ierror = 4
        return
      end if

      a_copy = a / g
      b_copy = b / g
      c_copy = c / g
c
c  Step 2: Split A and B into sign and magnitude.
c
      a_mag = abs ( a_copy )
      a_sign = sign ( 1, a_copy )
      b_mag = abs ( b_copy )
      b_sign = sign ( 1, b_copy )
c
c  Another special case, A_MAG = 1 or B_MAG = 1.
c
      if ( a_mag .eq. 1 ) then
        x = a_sign * c_copy
        return
      else if ( b_mag .eq. 1 ) then
        x = 0
        return
      end if
c
c  Step 3: Produce the Euclidean remainder sequence.
c
      if ( b_mag .le. a_mag ) then

        swap = .false.
        q(1) = a_mag
        q(2) = b_mag

      else

        swap = .true.
        q(1) = b_mag
        q(2) = a_mag

      end if

      n = 3

10    continue

        q(n) = mod ( q(n-2), q(n-1) )

        if ( q(n) .eq. 1 ) then
          go to 20
        end if

        n = n + 1

        if ( nmax .lt. n ) then
          ierror = 5
          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) 'CONGRUENCE - Fatal error!'
          write ( *, '(a)' ) '  Exceeded number of iterations.'
          stop
        end if

      go to 10

20    continue
c
c  Step 4: Go backwards to solve X * A_MAG + Y * B_MAG = 1.
c
      y = 0
      do k = n, 2, -1
        x = y
        y = ( 1 - x * q(k-1) ) / q(k)
      end do
c
c  Step 5: Undo the swapping.
c
      if ( swap ) then
        z = x
        x = y
        y = z
      end if
c
c  Step 6: Apply signs to X and Y so that X * A + Y * B = 1.
c
      x = x * a_sign
c
c  Step 7: Multiply by C, so that X * A + Y * B = C.
c
      x = x * c_copy
c
c  Step 8: Force 0 <= X < B.
c
      x = mod ( x, b )
c
c  Step 9: Force positivity.
c
      if ( x .lt. 0 ) then
        x = x + b
      end if

      return
      end
      subroutine count_pose_random ( seed, blocks, goal )

c*********************************************************************72
c
cc COUNT_POSE_RANDOM poses a problem for the game "The Count is Good"
c
c  Discussion:
c
c    The French television show "The Count is Good" has a game that goes
c    as follows:
c
c      A number is chosen at random between 100 and 999.  This is the GOAL.
c
c      Six numbers are randomly chosen from the set 1, 2, 3, 4, 5, 6, 7, 8,
c      9, 10, 25, 50, 75, 100.  These numbers are the BLOCKS.
c
c      The player must construct a formula, using some or all of the blocks,
c      (but not more than once), and the operations of addition, subtraction,
c      multiplication and division.  Parentheses should be used to remove
c      all ambiguity.  However, it is forbidden to use subtraction in a
c      way that produces a negative result, and all division must come out
c      exactly, with no remainder.
c
c    This routine poses a sample problem from the show.  The point is,
c    to determine how to write a program that can solve such a problem.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    08 January 2007
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Raymond Seroul,
c    Programming for Mathematicians,
c    Springer Verlag, 2000, pages 355-357.
c
c  Parameters:
c
c    Input/output, integer SEED, a seed for the random number generator.
c
c    Output, integer BLOCKS(6), the six numbers available for the formula.
c
c    Output, integer GOAL, the goal number.
c
      implicit none

      integer blocks(6)
      integer goal
      integer i
      integer i4_uniform
      integer ind(6)
      integer seed
      integer stuff(14)

      data stuff /
     &  1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 25, 50, 75, 100 /

      goal = i4_uniform ( 100, 999, seed )

      call ksub_random ( 14, 6, seed, ind )

      do i = 1, 6
        blocks(i) = stuff(ind(i))
      end do

      return
      end
      subroutine debruijn ( m, n, string )

c*********************************************************************72
c
cc DEBRUIJN constructs a de Bruijn sequence.
c
c  Discussion:
c
c    Suppose we have an alphabet of M letters, and we are interested in
c    all possible strings of length N.  If M = 2 and N = 3, then we are
c    interested in the M**N strings:
c
c      000
c      001
c      010
c      011
c      100
c      101
c      110
c      111
c
c    Now, instead of making a list like this, we prefer, if possible, to
c    write a string of letters, such that every consecutive sequence of
c    N letters is one of the strings, and every string occurs once, if
c    we allow wraparound.
c
c    For the above example, a suitable sequence would be the 8 characters:
c
c      00011101(00...
c
c    where we have suggested the wraparound feature by repeating the first
c    two characters at the end.
c
c    Such a sequence is called a de Bruijn sequence.  It can easily be
c    constructed by considering a directed graph, whose nodes are all
c    M**(N-1) strings of length N-1.  A node I has a directed edge to
c    node J (labeled with character K) if the string at node J can
c    be constructed by beheading the string at node I and adding character K.
c
c    In this setting, a de Bruijn sequence is simply an Eulerian circuit
c    of the directed graph, with the edge labels being the entries of the
c    sequence.  In general, there are many distinct de Bruijn sequences
c    for the same parameter M and N.  This program will only find one
c    of them.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    08 January 2007
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer M, the number of letters in the alphabet.
c
c    Input, integer N, the number of letters in a codeword.
c
c    Output, integer STRING(M**N), a deBruijn string.
c
      implicit none

      integer m
      integer n

      integer i
      integer iedge
      integer inode(m**n)
      integer ivec(n-1)
      integer j
      integer jnode(m**n)
      integer jvec(n-1)
      integer k
      integer knode(m**n)
      integer nedge
      integer nnode
      logical success
      integer string(m**n)
      integer trail(m**n)
c
c  Construct the adjacency information.
c
      nnode = m**(n-1)
      nedge = m**n

      iedge = 0

      do i = 1, nnode

        call index_unrank0 ( n-1, m, i, ivec )

        do k = 1, m
          do j = 1, n-2
            jvec(j) = ivec(j+1)
          end do
          jvec(n-1) = k
          call index_rank0 ( n-1, m, jvec, j )
          iedge = iedge + 1
          inode(iedge) = i
          jnode(iedge) = j
          knode(iedge) = k
        end do

      end do
c
c  Determine a circuit.
c
      call digraph_arc_euler ( nnode, nedge, inode, jnode, success, 
     &  trail )
c
c  The string is constructed from the labels of the edges in the circuit.
c
      do i = 1, nedge
        string(i) = knode(trail(i))
      end do

      return
      end
      subroutine dec_add ( mantissa1, exponent1, mantissa2, exponent2, 
     &  dec_digit, mantissa, exponent )

c*********************************************************************72
c
cc DEC_ADD adds two decimal quantities.
c
c  Discussion:
c
c    A decimal value is represented by MANTISSA * 10**EXPONENT.
c
c    The routine computes
c
c      MANTISSA * 10**EXPONENT = 
c        MANTISSA1 * 10**EXPONENT1 
c      + MANTISSA2 * 10**EXPONENT2
c
c    while trying to avoid integer overflow.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    08 January 2006
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer MANTISSA1, EXPONENT1, the first number to be added.
c
c    Input, integer MANTISSA2, EXPONENT2, the second number to be added.
c
c    Input, integer DEC_DIGIT, the number of decimal digits.
c
c    Output, integer MANTISSA, EXPONENT, the sum.
c
      implicit none

      integer dec_digit
      integer exponent
      integer exponent1
      integer exponent2
      integer mantissa
      integer mantissa1
      integer mantissa2
      integer mantissa3
      integer mantissa4

      if ( mantissa1 .eq. 0 ) then
        mantissa = mantissa2
        exponent = exponent2
        return
      else if ( mantissa2 .eq. 0 ) then
        mantissa = mantissa1
        exponent = exponent1
        return
      else if ( exponent1 .eq. exponent2 ) then
        mantissa = mantissa1 + mantissa2
        exponent = exponent1
        call dec_round ( mantissa, exponent, dec_digit, mantissa, 
     &    exponent )
        return
      end if
c
c  Line up the exponents.
c
      mantissa3 = mantissa1
      mantissa4 = mantissa2

      if ( exponent1 .lt. exponent2 ) then
        mantissa4 = mantissa4 * 10**( exponent2 - exponent1 )
      else
        mantissa3 = mantissa3 * 10**( exponent1 - exponent2 )
      end if
c
c  Add the coefficients.
c
      mantissa = mantissa3 + mantissa4
      exponent = min ( exponent1, exponent2 )
c
c  Clean up the result.
c
      call dec_round ( mantissa, exponent, dec_digit, mantissa, 
     &  exponent )

      return
      end
      subroutine dec_div ( mantissa1, exponent1, mantissa2, exponent2, 
     &  dec_digit, mantissa, exponent, ierror )

c*********************************************************************72
c
cc DEC_DIV divides two decimal values.
c
c  Discussion:
c
c    A decimal value is represented by MANTISSA * 10**EXPONENT.
c
c    The routine computes
c
c      MANTISSA * 10**EXPONENT
c      = ( MANTISSA1 * 10**EXPONENT1 ) / ( MANTISSA2 * 10**EXPONENT2 )
c      = ( MANTISSA1 / MANTISSA2 ) * 10**( EXPONENT1 - EXPONENT2 )
c
c    while avoiding integer overflow.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    09 January 2007
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer MANTISSA1, EXPONENT1, the numerator.
c
c    Input, integer MANTISSA2, EXPONENT2, the denominator.
c
c    Input, integer DEC_DIGIT, the number of decimal digits.
c
c    Output, integer MANTISSA, EXPONENT, the result.
c
c    Output, integer IERROR.
c    0, no error occurred.
c    1, an error occurred.
c
      implicit none

      integer dec_digit
      double precision dval
      integer exponent
      integer exponent1
      integer exponent2
      integer exponent3
      integer ierror
      integer mantissa
      integer mantissa1
      integer mantissa2
      integer mantissa3
c
c  First special case, top fraction is 0.
c
      if ( mantissa1 .eq. 0 ) then
        mantissa = 0
        exponent = 0
        return
      end if
c
c  First error, bottom of fraction is 0.
c
      if ( mantissa2 .eq. 0 ) then
        ierror = 1
        mantissa = 0
        exponent = 0
        return
      end if
c
c  Second special case, result is 1.
c
      if ( mantissa1 .eq. mantissa2 .and. 
     &     exponent1 .eq. exponent2 ) then
        mantissa = 1
        exponent = 0
        return
      end if
c
c  Third special case, result is power of 10.
c
      if ( mantissa1 .eq. mantissa2 ) then
        mantissa = 1
        exponent = exponent1 - exponent2
        return
      end if
c
c  Fourth special case: MANTISSA1/MANTISSA2 is exact.
c
      if ( ( mantissa1 / mantissa2 ) * mantissa2 .eq. mantissa1 ) then
        mantissa = mantissa1 / mantissa2
        exponent = exponent1 - exponent2
        return
      end if
c
c  General case.
c
      dval = dble ( mantissa1 ) / dble ( mantissa2 )

      call r8_to_dec ( dval, dec_digit, mantissa3, exponent3 )

      mantissa = mantissa3
      exponent = exponent3 + exponent1 - exponent2

      return
      end
      subroutine dec_mul ( mantissa1, exponent1, mantissa2, exponent2, 
     &  dec_digit, mantissa, exponent )

c*********************************************************************72
c
cc DEC_MUL multiplies two decimals.
c
c  Discussion:
c
c    A decimal value is represented by MANTISSA * 10**EXPONENT.
c
c    The routine computes
c
c      MANTISSA * 10**EXPONENT 
c      = ( MANTISSA1 * 10**EXPONENT1) * (MANTISSA2 * 10**EXPONENT2)
c      = ( MANTISSA1 * MANTISSA2 ) * 10**( EXPONENT1 + EXPONENT2 )
c
c    while avoiding integer overflow.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    09 January 2007
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer MANTISSA1, EXPONENT1, the first multiplier.
c
c    Input, integer MANTISSA2, EXPONENT2, the second multiplier.
c
c    Input, integer DEC_DIGIT, the number of decimal digits.
c
c    Output, integer MANTISSA, EXPONENT, the product.
c
      implicit none

      integer dec_digit
      double precision dval
      integer exponent
      integer exponent1
      integer exponent2
      integer exponent3
      integer i_max
      integer i4_huge
      integer mantissa
      integer mantissa1
      integer mantissa2
      integer mantissa3
      double precision temp

      i_max = i4_huge ( )
c
c  The result is zero if either MANTISSA1 or MANTISSA2 is zero.
c
      if ( mantissa1 .eq. 0 .or. mantissa2 .eq. 0 ) then
        mantissa = 0
        exponent = 0
        return
      end if
c
c  The result is simple if either MANTISSA1 or MANTISSA2 is one.
c
      if ( abs ( mantissa1 ) .eq. 1 .or. 
     &     abs ( mantissa2 ) .eq. 1 ) then
        mantissa = mantissa1 * mantissa2
        exponent = exponent1 + exponent2
        return
      end if

      temp = log ( dble ( abs ( mantissa1 ) ) ) 
     &     + log ( dble ( abs ( mantissa2 ) ) )

      if ( temp .lt. log ( dble ( i_max ) ) ) then

        mantissa = mantissa1 * mantissa2
        exponent = exponent1 + exponent2

      else

        dval = dble ( mantissa1 ) * dble ( mantissa2 )

        call r8_to_dec ( dval, dec_digit, mantissa3, exponent3 )

        mantissa = mantissa3
        exponent = exponent3 + ( exponent1 + exponent2 )

      end if

      call dec_round ( mantissa, exponent, dec_digit, mantissa, 
     &  exponent )

      return
      end
      subroutine dec_round ( mantissa1, exponent1, dec_digit, 
     &  mantissa2, exponent2 )

c*********************************************************************72
c
cc DEC_ROUND rounds a decimal fraction to a given number of digits.
c
c  Discussion:
c
c    A decimal value is represented by MANTISSA * 10**EXPONENT.
c
c    The routine takes an arbitrary decimal fraction makes sure that 
c    MANTISSA has no more than DEC_DIGIT digits.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    09 January 2007
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer MANTISSA1, EXPONENT1, the coefficient and exponent
c    of a decimal fraction to be rounded.
c
c    Input, integer DEC_DIGIT, the number of decimal digits.
c
c    Output, integer MANTISSA2, EXPONENT2, the rounded coefficient and 
c    exponent of a decimal fraction.  MANTISSA2 has no more than
c    DEC_DIGIT decimal digits.
c
      implicit none

      integer dec_digit
      integer exponent1
      integer exponent2
      integer mantissa1
      integer mantissa2

      mantissa2 = mantissa1
      exponent2 = exponent1

      if ( mantissa2 .eq. 0 ) then
        exponent2 = 0
        return
      end if

10    continue

      if ( 10**dec_digit .le. abs ( mantissa2 ) ) then
        mantissa2 = nint ( dble ( mantissa2 ) / 10.0D+00 )
        exponent2 = exponent2 + 1
        go to 10
      end if
c
c  Absorb trailing 0's into the exponent.
c
20    continue

      if ( ( mantissa2 / 10 ) * 10 .eq. mantissa2 ) then
        mantissa2 = mantissa2 / 10
        exponent2 = exponent2 + 1
        go to 20
      end if

      return
      end
      subroutine dec_to_r8 ( mantissa, exponent, r )

c*********************************************************************72
c
cc DEC_TO_R8 converts a decimal to an R8.
c
c  Discussion:
c
c    A decimal value is represented by MANTISSA * 10**EXPONENT.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    09 January 2007
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer MANTISSA, EXPONENT, the coefficient and exponent
c    of the decimal value.
c
c    Output, double precision R, the equivalent real value.
c
      implicit none

      integer exponent
      integer mantissa
      double precision r

      r = mantissa * 10.0D+00**exponent

      return
      end
      subroutine dec_to_rat ( mantissa, exponent, rat_top, rat_bot )

c*********************************************************************72
c
cc DEC_TO_RAT converts a decimal to a rational representation.
c
c  Discussion:
c
c    A decimal value is represented by MANTISSA * 10**EXPONENT.
c
c    A rational value is represented by RAT_TOP / RAT_BOT.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    09 January 2007
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer MANTISSA, EXPONENT, the decimal number.
c
c    Output, integer RAT_TOP, RAT_BOT, the rational value.
c
      implicit none

      integer gcd
      integer exponent
      integer i4_gcd
      integer mantissa
      integer rat_bot
      integer rat_top

      if ( exponent .eq. 0 ) then
        rat_top = mantissa
        rat_bot = 1
      else if ( 0 .lt. exponent ) then
        rat_top = mantissa * 10**exponent
        rat_bot = 1
      else
        rat_top = mantissa
        rat_bot = 10**( - exponent )
        gcd = i4_gcd ( rat_top, rat_bot )
        rat_top = rat_top / gcd
        rat_bot = rat_bot / gcd
      end if

      return
      end
      subroutine dec_to_s ( mantissa, exponent, s )

c*********************************************************************72
c
cc DEC_TO_S returns a string representation of a decimal.
c
c  Discussion:
c
c    A decimal value is represented by MANTISSA * 10**EXPONENT.
c
c  Example:
c
c    MANTISSA EXPONENT   S
c    ----     ----       ------
c       0        0       0
c      21        3       21000
c      -3        0       -3
c     147       -2       14.7
c      16       -5       0.00016
c      34       30       Inf
c     123      -21       0.0000000000000000012
c      34      -30       0.0
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    09 January 2007
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer MANTISSA, EXPONENT, integers which represent the decimal.
c
c    Output, character*(*) S, the representation of the value.
c    The string is 'Inf' or '0.0' if the value was too large
c    or small to represent with a fixed point format.
c
      implicit none

      character*(22) chrrep
      integer exponent
      integer i
      integer iget1
      integer iget2
      integer iput1
      integer iput2
      integer mantissa
      integer maxdigit
      integer ndigit
      integer nleft
      character*(*) s
      integer s_len_trim

      s = ' '

      if ( mantissa .eq. 0 ) then
        s = '0'
        return
      end if

      maxdigit = len ( s )
c
c  Store a representation of MANTISSA in CHRREP.
c
      write ( chrrep, '(i22)' ) mantissa
      call s_blank_delete ( chrrep )
      ndigit = s_len_trim ( chrrep )
c
c  Overflow if EXPONENT is positive, and MAXDIGIT .lt. NDIGIT + EXPONENT.
c
      if ( 0 .lt. exponent ) then
        if ( maxdigit .lt. ndigit + exponent ) then
          s = 'Inf'
          return
        end if
      end if
c
c  Underflow if EXPONENT is negative, and MAXDIGIT .lt. 3 + NDIGIT - EXPONENT.
c
      if ( exponent .lt. 0 ) then
        if ( 0 .lt. mantissa ) then
          if ( maxdigit .lt. 3 - ndigit - exponent ) then
            s = '0.0'
            return
          end if
        else
          if ( maxdigit .lt. 5 - ndigit - exponent ) then
            s = '0.0'
            return
          end if
        end if
      end if
c
c  If EXPONENT is nonnegative, insert trailing zeros.
c
      if ( 0 .le. exponent ) then

        s(1:ndigit) = chrrep(1:ndigit)

        do i = ndigit + 1, ndigit + exponent
          s(i:i) = '0'
        end do

      else if ( exponent .lt. 0 ) then

        iput2 = 0
        iget2 = 0
c
c  Sign.
c
        if ( mantissa .lt. 0 ) then
          iput1 = 1
          iput2 = 1
          iget2 = 1
          s(iput1:iput2) = '-'
          ndigit = ndigit - 1
        end if
c
c  Digits of the integral part.
c
        if ( 0 .lt. ndigit + exponent ) then
          iput1 = iput2 + 1
          iput2 = iput1 + ndigit + exponent -1
          iget1 = iget2 + 1
          iget2 = iget1 + ndigit + exponent - 1
          s(iput1:iput2) = chrrep(iget1:iget2)
        else
          iput1 = iput2 + 1
          iput2 = iput1
          s(iput1:iput2) = '0'
        end if
c
c  Decimal point.
c
        iput1 = iput2 + 1
        iput2 = iput1
        s(iput1:iput2) = '.'
c
c  Leading zeroes.
c
        do i = 1, - exponent - ndigit
          iput1 = iput2 + 1
          iput2 = iput1
          s(iput1:iput2) = '0'
        end do

        nleft = min ( -exponent, ndigit )
        nleft = min ( nleft, maxdigit - iput2 )
        iput1 = iput2 + 1
        iput2 = iput1 + nleft - 1
        iget1 = iget2 + 1
        iget2 = iget1 + nleft - 1
        s(iput1:iput2) = chrrep(iget1:iget2)

      end if

      return
      end
      function dec_width ( mantissa, exponent )

c*********************************************************************72
c
cc DEC_WIDTH returns the "width" of a decimal number.
c
c  Discussion:
c
c    A decimal value is represented by MANTISSA * 10**EXPONENT.
c
c    The "width" of a decimal number is the number of characters
c    required to print it.
c
c  Example:
c
c    Mantissa  Exponent Width  Representation:
c
c         523      -1       4           5.23
c         134       2       5       13400
c           0      10       1           0
c      123456     -10      12           0.0000123456
c      123456      -3       7         123.456
c      123456       0       6      123456
c      123456       3       9   123456000
c 
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    09 January 2007
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer MANTISSA, EXPONENT, the decimal number.
c
c    Output, integer DEC_WIDTH, the "width" of the decimal number.
c
      implicit none

      integer dec_width
      integer exponent
      integer mantissa
      integer mantissa_abs
      integer ten_pow
      integer value

      value = 1
      ten_pow = 10

      if ( mantissa .eq. 0 ) then
        dec_width = value
        return
      end if
      
      mantissa_abs = abs ( mantissa )

10    continue

      if ( ten_pow .le. mantissa_abs ) then
        value = value + 1
        ten_pow = ten_pow * 10
        go to 10
      end if

      if ( 0 .lt. exponent ) then
        value = value + exponent
      else if ( exponent .lt. 0 ) then
        value = max ( value, 1-exponent )
c
c  An internal decimal point adds one position.
c
        if ( 0 .lt. value ) then
          value = value + 1
c
c  A leading "0." adds two positions.
c
        else
          value = 2 - value
        end if
      end if

      if ( mantissa .lt. 0 ) then
        value = value + 1
      end if

      dec_width = value

      return
      end
      subroutine decmat_det ( n, atop, abot, dec_digit, dtop, dbot, 
     &  ierror )

c*********************************************************************72
c
cc DECMAT_DET finds the determinant of an N by N matrix of decimal entries.
c
c  Discussion:
c
c    The brute force method is used.  The routine should only be used for
c    small matrices, since this calculation requires the summation of N!
c    products of N numbers.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    09 January 2007
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the number of rows and columns of A.
c
c    Input, integer ATOP(N,N), ABOT(N,N), the decimal
c    representation of the matrix.
c
c    Output, integer DTOP, DBOT, the decimal determinant of the matrix.
c
c    Output, integer IERROR.
c    0, no error occurred.
c    1, an error occurred (probably overflow).
c
      implicit none

      integer n

      integer dec_digit
      logical even
      integer i
      integer abot(n,n)
      integer atop(n,n)
      integer iarray(n)
      integer ibot
      integer ibot1
      integer ibot2
      integer dbot
      integer dtop
      integer ierror
      integer itop
      integer itop1
      integer itop2
      logical more

      ierror = 0
      more = .false.
      dtop = 0
      dbot = 1
c
c  Compute the next permutation.
c
10    continue

        call perm_next ( n, iarray, more, even )
c
c  The sign of this term depends on the sign of the permutation.
c
        if ( even ) then
          itop = 1
        else
          itop = -1
        end if
c
c  Choose one item from each row, as specified by the permutation,
c  and multiply them
c
        ibot = 0

        do i = 1, n

          itop1 = itop
          ibot1 = ibot
          itop2 = atop(i,iarray(i))
          ibot2 = abot(i,iarray(i))

          call dec_mul ( itop1, ibot1, itop2, ibot2, dec_digit, 
     &      itop, ibot )

        end do
c
c  Add this term to the total.
c
        itop1 = itop
        ibot1 = ibot

        call dec_add ( itop1, ibot1, dtop, dbot, dec_digit, 
     &    itop, ibot )

        dtop = itop
        dbot = ibot

        if ( .not. more ) then
          go to 20
        end if

      go to 10

20    continue

      return
      end
      subroutine decmat_print ( m, n, a, b, title )

c*********************************************************************72
c
cc DECMAT_PRINT prints out decimal vectors and matrices.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    09 January 2007
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer M, N, the number of rows and columns in the matrix.
c
c    Input, integer A(M,N), B(M,N), the decimal matrix.
c
c    Input, character*(*) TITLE, a label for the object being printed.
c
      implicit none

      integer m
      integer n

      integer a(m,n)
      integer b(m,n)
      character * ( 22 ) chrtmp
      character * ( 10 ) chrtmp2
      character * ( 10 ) chrtmp3
      character * ( 40 ) format2
      integer i
      integer imax
      integer j
      integer jmax
      integer jmin
      integer khi
      integer klo
      integer kmax
      integer lenc
      integer ncolum
      parameter ( ncolum = 80 )
      integer npline
      character * ( 100 ) output
      integer s_len_trim
      character * ( * ) title
      integer title_length
c
c  Figure out how wide we must make each column.
c
      imax = 0
      jmax = 0

      do i = 1, m
        do j = 1, n

          call dec_to_s ( a(i,j), b(i,j), chrtmp )
          lenc = s_len_trim ( chrtmp )
          jmax = max ( jmax, lenc )

        end do
      end do

      kmax = 2 + imax + 1 + jmax
      npline = ncolum / kmax
c
c  Set up the format for the heading.
c
      call i4_to_s_left ( npline, chrtmp2 )
      call i4_to_s_left ( kmax, chrtmp3 )
      format2 = '(' // chrtmp2 // 'i' // chrtmp3 // ')'

      call s_blank_delete ( format2 )

      do jmin = 1, n, npline

        jmax = min ( jmin + npline - 1, n )

        title_length = s_len_trim ( title )
        if ( 0 .lt. title_length ) then
          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) title(1:title_length)
        end if

        write ( *, '(a)' ) ' '

        if ( 1 .lt. jmin .or. jmax .lt. n ) then
          write ( output, * ) 'Columns ', jmin, ' to ', jmax
          call s_blanks_delete ( output )
          write ( *, '(a)' ) output
          write ( *, '(a)' ) ' '
        end if

        do i = 1, m

          output = ' '

          do j = jmin, jmax
            klo = 4 + ( j - jmin ) * kmax + 1
            khi = 4 + ( j - jmin ) * kmax + kmax
            chrtmp = ' '
            call dec_to_s ( a(i,j), b(i,j), chrtmp )
            call s_adjustr ( chrtmp(1:kmax) )
            output(klo:khi) = chrtmp(1:kmax)
          end do

          write ( *, '(a)' ) output

        end do

      end do

      return
      end
      subroutine derange_back_candidate ( n, maxstack, a, k, nstack, 
     &  stack, ncan )

c*********************************************************************72
c
cc DERANGE_BACK_CANDIDATE finds possible values for the K-th entry of a derangement.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    10 January 2007
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the order of the derangement.
c
c    Input, integer MAXSTACK, the maximum stack length.
c
c    Input, integer A(N).  The first K-1 entries of A
c    record the currently set values of the derangement.
c
c    Input, integer K, the entry of the derangement for which candidates
c    are to be found.
c
c    Input/output, integer NSTACK, the length of the stack.
c
c    Input/output, integer STACK(MAXSTACK).  On output, we have added
c    the candidates for entry K to the end of the stack.
c
c    Input/output, integer NCAN(N), the number of candidates for each level.
c
      implicit none

      integer maxstack
      integer n

      integer a(n)
      integer ican
      integer ifree(n)
      integer k
      integer ncan(n)
      integer nfree
      integer nstack
      integer stack(maxstack)
c
c  Consider all the integers from 1 through N that have not been used yet.
c
      nfree = n - k + 1

      call perm_free ( k-1, a, nfree, ifree )
c
c  Everything but K is a legitimate candidate for the K-th entry.
c
      ncan(k) = 0

      do ican = 1, nfree

        if ( ifree(ican) .ne. k ) then

          ncan(k) = ncan(k) + 1
          nstack = nstack + 1

          if ( maxstack .lt. nstack ) then
            write ( *, '(a)' ) ' '
            write ( *, '(a)' ) 'DERANGE_BACK_CANDIDATE - Fatal error!'
            write ( *, '(a,i8)' ) 
     &      '  Exceeding stacksize limit of ', maxstack
            stop
          end if

          stack(nstack) = ifree(ican)

        end if

      end do

      return
      end
      subroutine derange_back_next ( n, a, more )

c*********************************************************************72
c
cc DERANGE_BACK_NEXT returns the next derangement of N items.
c
c  Discussion:
c
c    This routine uses backtracking.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    27 January 2007
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the number of items to be deranged.  N should
c    be at least 2.
c
c    Input/output, integer A(N).
c    On the first call, the input value of A is not important.
c    On return with MORE = TRUE, A contains the next derangement.
c    On subsequent input, A should not be changed.
c
c    Input/output, logical MORE.
c    On first call, set MORE to FALSE, and do not alter it after.
c    On return, MORE is TRUE if another derangement is being treturned in A,
c    and FALSE if no more can be found.
c
      implicit none

      integer n
      integer n_max
      parameter ( n_max = 100 )
      integer n2_max
      parameter ( n2_max = 50 * 101 )

      integer a(n)
      integer i
      integer indx
      integer k
      integer maxstack
      logical more
      integer ncan(n_max)
      integer nstack
      integer stack(n2_max)

      save indx
      save k
      save maxstack
      save ncan
      save nstack
      save stack

      data indx / 0 /
      data k / 0 /
      data maxstack / 0 /
      data nstack / 0 /

      if ( .not. more ) then

        if ( n_max .lt. n ) then
          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) 'DERANGE_BACK_NEXT - Fatal error!'
          write ( *, '(a)' ) '  Input N exceeds internal limit.'
          stop
        end if

        if ( n .lt. 2 ) then
          more = .false.
          return
        end if

        indx = 0
        k = 0
        maxstack = ( n * ( n + 1 ) ) / 2
        nstack = 0

        do i = 1, maxstack
          stack(i) = 0
        end do

        do i = 1, n
          ncan(i) = 0
        end do

        more = .true.

      end if

10    continue

        call i4vec_backtrack ( n, maxstack, stack, a, indx, k, 
     &    nstack, ncan )

        if ( indx .eq. 1 ) then

          go to 20

        else if ( indx .eq. 2 ) then

          call derange_back_candidate ( n, maxstack, a, k, nstack, 
     &      stack, ncan )

        else

          more = .false.
          go to 20

        end if

      go to 10

20    continue

      return
      end
      subroutine derange_check ( n, a, deranged )

c*********************************************************************72
c
cc DERANGE_CHECK determines whether a permutation is a derangement.
c
c  Discussion:
c
c    A derangement of the integers 1 through N is a permutation of the
c    integers such that the first value is not 1, the second is not 2,
c    and so on.
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
c    Input, integer N, the number of objects permuted.
c
c    Input, integer A(N), a permutation of the integers 1 through N.
c
c    Output, logical DERANGED, is TRUE if A is a derangement, and
c    FALSE otherwise.
c
      implicit none

      integer n

      integer a(n)
      logical deranged
      integer i

      do i = 1, n
        if ( a(i) .eq. i ) then
          deranged = .false.
          return
        end if
      end do

      deranged = .true.

      return
      end
      function derange_enum ( n )

c*********************************************************************72
c
cc DERANGE_ENUM returns the number of derangements of N objects.
c
c  Discussion:
c
c    A derangement of N objects is a permutation with no fixed
c    points.  If we symbolize the permutation operation by "P",
c    then for a derangment, P(I) is never equal to I.
c
c    D(N) is the number of ways of placing N non-attacking rooks on
c    an N by N chessboard with one diagonal deleted.
c
c    Limit ( N -> Infinity ) D(N)/Nc = 1 / e.
c
c    The number of permutations with exactly K items in the right
c    place is COMB(N,K) * D(N-K).
c
c  Recursion:
c
c      D(0) = 1
c      D(1) = 0
c      D(2) = 1
c      D(N) = (N-1) * ( D(N-1) + D(N-2) )
c
c    or
c
c      D(0) = 1
c      D(1) = 0
c      D(N) = N * D(N-1) + (-1)**N
c
c  Formula:
c
c    D(N) = Nc * ( 1 - 1/1c + 1/2c - 1/3c ... 1/Nc )
c
c    Based on the inclusion/exclusion law.
c
c    D(N) = nint ( Nc / E )
c
c  First values:
c
c     N         D(N)
c     0           1
c     1           0
c     2           1
c     3           2
c     4           9
c     5          44
c     6         265
c     7        1854
c     8       14833
c     9      133496
c    10     1334961
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
c    Input, integer N, the number of objects to be permuted.
c
c    Output, integer DERANGE_ENUM, the number of derangements of N objects.
c
      implicit none

      integer derange_enum
      integer dn
      integer dnm1
      integer dnm2
      integer i
      integer n

      if ( n .lt. 0 ) then

        dn = 0

      else if ( n .eq. 0 ) then

        dn = 1

      else if ( n .eq. 1 ) then

        dn = 0

      else if ( n .eq. 2 ) then

        dn = 1

      else

        dnm1 = 0
        dn = 1

        do i = 3, n
          dnm2 = dnm1
          dnm1 = dn
          dn = ( i - 1 ) * ( dnm1 + dnm2 )
        end do

      end if

      derange_enum = dn

      return
      end
      subroutine derange_enum2 ( n, d )

c*********************************************************************72
c
cc DERANGE_ENUM2 returns the number of derangements of 0 through N objects.
c
c  Discussion:
c
c    A derangement of N objects is a permutation with no fixed
c    points.  If we symbolize the permutation operation by "P",
c    then for a derangment, P(I) is never equal to I.
c
c    D(N) is the number of ways of placing N non-attacking rooks on
c    an N by N chessboard with one diagonal deleted.
c
c    Limit ( N -> Infinity ) D(N)/Nc = 1 / e.
c
c    The number of permutations with exactly K items in the right
c    place is COMB(N,K) * D(N-K).
c
c  Recursion:
c
c      D(0) = 1
c      D(1) = 0
c      D(2) = 1
c      D(N) = (N-1) * ( D(N-1) + D(N-2) )
c
c    or
c
c      D(0) = 1
c      D(1) = 0
c      D(N) = N * D(N-1) + (-1)**N
c
c  Formula:
c
c    D(N) = Nc * ( 1 - 1/1c + 1/2c - 1/3c ... 1/Nc )
c
c    Based on the inclusion/exclusion law.
c
c    D(N) = nint ( Nc / E )
c
c  First values:
c
c     N         D(N)
c     0           1
c     1           0
c     2           1
c     3           2
c     4           9
c     5          44
c     6         265
c     7        1854
c     8       14833
c     9      133496
c    10     1334961
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
c    Input, integer N, the maximum number of objects to be permuted.
c
c    Output, integer D(0:N); D(I) is the number of derangements of
c    I objects.
c
      implicit none

      integer n

      integer d(0:n)
      integer i

      d(0) = 1
      d(1) = 0

      do i = 2, n
        d(i) = ( i - 1 ) * ( d(i-1) + d(i-2) )
      end do

      return
      end
      function derange_enum3 ( n )

c*********************************************************************72
c
cc DERANGE_ENUM3 returns the number of derangements of 0 through N objects.
c
c  Discussion:
c
c    A derangement of N objects is a permutation with no fixed
c    points.  If we symbolize the permutation operation by "P",
c    then for a derangment, P(I) is never equal to I.
c
c    D(N) is the number of ways of placing N non-attacking rooks on
c    an N by N chessboard with one diagonal deleted.
c
c    Limit ( N -> Infinity ) D(N)/Nc = 1 / e.
c
c    The number of permutations with exactly K items in the right
c    place is COMB(N,K) * D(N-K).
c
c  Recursion:
c
c      D(0) = 1
c      D(1) = 0
c      D(2) = 1
c      D(N) = (N-1) * ( D(N-1) + D(N-2) )
c
c    or
c
c      D(0) = 1
c      D(1) = 0
c      D(N) = N * D(N-1) + (-1)**N
c
c  Formula:
c
c    D(N) = Nc * ( 1 - 1/1c + 1/2c - 1/3c ... 1/Nc )
c
c    Based on the inclusion/exclusion law.
c
c    D(N) = nint ( Nc / E )
c
c  First values:
c
c     N         D(N)
c     0           1
c     1           0
c     2           1
c     3           2
c     4           9
c     5          44
c     6         265
c     7        1854
c     8       14833
c     9      133496
c    10     1334961
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
c    Input, integer N, the maximum number of objects to be permuted.
c
c    Output, integer DERANGE_ENUM3, the number of derangements of N objects.
c
      implicit none

      integer derange_enum3
      double precision e
      parameter ( e = 2.718281828459045D+00 )
      integer n
      double precision r8_factorial

      if ( n .lt. 0 ) then
        derange_enum3 = -1
      else if ( n .eq. 0 ) then
        derange_enum3 = 1
      else if ( n .eq. 1 ) then
        derange_enum3 = 0
      else
        derange_enum3 = nint ( r8_factorial ( n ) / e )
      end if

      return
      end
      subroutine derange_weed_next ( n, a, more )

c*********************************************************************72
c
cc DERANGE_WEED_NEXT computes all of the derangements of N objects, one at a time.
c
c  Discussion:
c
c    A derangement of the integers 1 through N is a permutation of the
c    integers such that the first value is not 1, the second is not 2,
c    and so on.
c
c    This routine simply generates all permutations, one at a time,
c    and weeds out those that are not derangements.
c
c  Example:
c
c    Here are the derangements when N = 4:
c
c    2143  3142  4123
c    2341  3412  4312
c    2413  3421  4321
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    16 January 2007
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the number of objects being permuted.
c
c    Input/output, integer A(N).
c    On first call, the input contents of A are unimportant.  But
c    on the second and later calls, the input value of A should be
c    the output value returned on the previous call.
c    On output, A contains the next derangement.
c
c    Input/output, logical MORE.
c    Set MORE = FALSE before the first call.
c    MORE will be reset to TRUE and a derangement will be returned.
c    Each new call produces a new derangement until MORE is returned FALSE.
c
      implicit none

      integer n

      integer a(n)
      logical deranged
      integer derange_enum
      integer maxder
      logical more
      integer numder

      save maxder
      save numder

      data maxder / 0 /
      data numder / 0 /
c
c  Initialization on call with MORE = FALSE.
c
      if ( .not. more ) then

        maxder = derange_enum ( n )
        numder = 0

      end if
c
c  Watch out for cases where there are no derangements.
c
      if ( maxder .eq. 0 ) then
        more = .false.
        return
      end if
c
c  Get the next permutation.
c
10    continue

        call perm_lex_next ( n, a, more )
c
c  See if it is a derangment.
c
        call derange_check ( n, a, deranged )

        if ( deranged ) then
          go to 20
        end if

      go to 10

20    continue

      numder = numder + 1

      if ( maxder .le. numder ) then
        more = .false.
      end if

      return
      end
      subroutine digit_to_ch ( digit, c )

c*********************************************************************72
c
cc DIGIT_TO_CH returns the character representation of a decimal digit.
c
c  Example:
c
c    DIGIT   C
c    -----  ---
c      0    '0'
c      1    '1'
c    ...    ...
c      9    '9'
c     17    '*'
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    04 August 1999
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer DIGIT, the digit value between 0 and 9.
c
c    Output, character C, the corresponding character, or '*' if DIGIT
c    was illegal.
c
      implicit none

      character c
      integer digit

      if ( 0 .le. digit .and. digit .le. 9 ) then

        c = char ( digit + 48 )

      else

        c = '*'

      end if

      return
      end
      subroutine digraph_arc_euler ( nnode, nedge, inode, jnode, 
     &  success, trail )

c*********************************************************************72
c
cc DIGRAPH_ARC_EULER returns an Euler circuit in a digraph.
c
c  Discussion:
c
c    An Euler circuit of a digraph is a path which starts and ends at
c    the same node and uses each directed edge exactly once.  A digraph is
c    eulerian if it has an Euler circuit.  The problem is to decide whether
c    a given digraph is eulerian and to find an Euler circuit if the
c    answer is affirmative.
c
c
c    A digraph has an Euler circuit if and only if the number of incoming
c    edges is equal to the number of outgoing edges at each node.
c
c    This characterization gives a straightforward procedure to decide whether
c    a digraph is eulerian.  Furthermore, an Euler circuit in an eulerian
c    digraph G of NEDGE edges can be determined by the following method:
c
c      STEP 1: Choose any node U as the starting node, and traverse any edge
c      ( U, V ) incident to node U, and than traverse any unused edge incident
c      to node U.  Repeat this process of traversing unused edges until the
c      starting node U is reached.  Let P be the resulting walk consisting of
c      all used edges.  If all edges of G are in P, than stop.
c
c      STEP 2: Choose any unused edge ( X,  Y) in G such that X is
c      in P and Y is not in P.  Use node X as the starting node and
c      find another walk Q using all unused edges as in step 1.
c
c      STEP 3: Walk P and walk Q share a common node X, they can be merged
c      to form a walk R by starting at any node S of P and to traverse P
c      until node X is reached; than, detour and traverse all edges of Q
c      until node X is reached and continue to traverse the edges of P until
c      the starting node S is reached.  Set P = R.
c
c      STEP 4: Repeat steps 2 and 3 until all edges are used.
c
c    The running time of the algorithm is O ( NEDGE ).
c
c    The digraph is assumed to be connected.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    16 January 2007
c
c  Author:
c
c    Original FORTRAN77 version by Hang Tong Lau.
c    This FORTRAN77 version by John Burkardt.
c
c  Reference:
c
c    Hang Tong Lau,
c    Algorithms on Graphs,
c    Tab Books, 1989.
c
c  Parameters:
c
c    Input, integer NNODE, the number of nodes.
c
c    Input, integer NEDGE, the number of edges.
c
c    Input, integer INODE(NEDGE), JNODE(NEDGE); the I-th edge starts at node
c    INODE(I) and ends at node JNODE(I).
c
c    Output, logical SUCCESS, is TRUE if an Euler circuit was found,
c    and FALSE otherwise.
c
c    Output, integer TRAIL(NEDGE).  TRAIL(I) is the edge number of the I-th
c    edge in the Euler circuit.
c
      implicit none

      integer nedge
    
      logical candid(nedge)
      integer edge
      integer endnod(nedge)
      integer i
      integer inode(nedge)
      integer istak
      integer j
      integer jnode(nedge)
      integer k
      integer l
      integer len
      integer lensol
      integer lenstk
      integer nnode
      integer stack(2*nedge)
      logical success
      integer trail(nedge)
c
c  Check if the digraph is eulerian.
c
      do edge = 1, nedge
        trail(edge) = 0
      end do

      do edge = 1, nedge
        endnod(edge) = 0
      end do

      do i = 1, nedge
        j = inode(i)
        trail(j) = trail(j) + 1
        j = jnode(i)
        endnod(j) = endnod(j) + 1
      end do

      do i = 1, nnode
        if ( trail(i) .ne. endnod(i) ) then
          success = .false.
          return
        end if
      end do
c
c  The digraph is eulerian; find an Euler circuit.
c
      success = .true.
      lensol = 1
      lenstk = 0
c
c  Find the next edge.
c
10    continue

        if ( lensol .eq. 1 ) then

          endnod(1) = inode(1)
          stack(1) = 1
          stack(2) = 1
          lenstk = 2

        else

          l = lensol - 1

          if ( lensol .ne. 2 ) then
            endnod(l) = inode(trail(l)) 
     &        + jnode(trail(l)) - endnod(l-1)
          end if

          k = endnod(l)

          do i = 1, nedge
            candid(i) = ( k .eq. jnode(i) )
          end do

          do i = 1, l
            candid(trail(i)) = .false.
          end do

          len = lenstk

          do i = 1, nedge

            if ( candid(i) ) then
              len = len + 1
              stack(len) = i
            end if

          end do

          stack(len+1) = len - lenstk
          lenstk = len + 1

        end if

20      continue

          istak = stack(lenstk)
          lenstk = lenstk - 1

          if ( istak .ne. 0 ) then
            go to 30
          end if

          lensol = lensol - 1

          if ( lensol .eq. 0 ) then
            call i4vec_reverse ( nedge, trail )
            return
          end if

        go to 20

        trail(lensol) = stack(lenstk)
        stack(lenstk) = istak - 1

        if ( lensol .eq. nedge ) then
          go to 30
        end if

        lensol = lensol + 1

      go to 10

30    continue

      call i4vec_reverse ( nedge, trail )

      return
      end
      subroutine digraph_arc_print ( nedge, inode, jnode, title )

c*********************************************************************72
c
cc DIGRAPH_ARC_PRINT prints out a digraph from an edge list.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    17 January 2007
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer NEDGE, the number of edges.
c
c    Input, integer INODE(NEDGE), JNODE(NEDGE), the beginning and end
c    nodes of the edges.
c
c    Input, character * ( * ) TITLE, a title.
c
      implicit none

      integer nedge

      integer i
      integer inode(nedge)
      integer jnode(nedge)
      integer s_len_trim
      integer s_length
      character * ( * ) title

      s_length = s_len_trim ( title )

      if ( 0 .lt. s_length ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) title(1:s_length)
      end if

      write ( *, '(a)' ) ' '

      do i = 1, nedge
        write ( *, '(i8,4x,2i8)' ) i, inode(i), jnode(i)
      end do

      return
      end
      subroutine diophantine ( a, b, c, ierror, x, y )

c*********************************************************************72
c
cc DIOPHANTINE solves a Diophantine equation A * X + B * Y = C.
c
c  Discussion:
c
c    Given integers A, B and C, produce X and Y so that
c
c      A * X + B * Y = C.
c
c    In general, the equation is solvable if and only if the
c    greatest common divisor of A and B also divides C.
c
c    A solution (X,Y) of the Diophantine equation also gives the solution
c    X to the congruence equation:
c
c      A * X = C mod ( B ).
c
c    Generally, if there is one nontrivial solution, there are an infinite
c    number of solutions to a Diophantine problem.
c    If (X0,Y0) is a solution, then so is ( X0+T*B/D, Y0-T*A/D ) where
c    T is any integer, and D is the greatest common divisor of A and B.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    17 January 2007
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Eric Weisstein, editor,
c    CRC Concise Encylopedia of Mathematics,
c    CRC Press, 1998, page 446.
c
c  Parameters:
c
c    Input, integer A, B, C, the coefficients of the Diophantine equation.
c
c    Output, integer IERROR, error flag.
c    0, no error, X and Y were computed.
c    1, A = B = 0, C is nonzero.
c    2, A = 0, B and C nonzero, but C is not a multiple of B.
c    3, A nonzero, B zero, C nonzero, but C is not a multiple of A.
c    4, A, B, C nonzero, but GCD of A and B does not divide C.
c    5, the algorithm ran out of internal space.
c
c    Output, integer X, Y, the solution of the Diophantine equation.
c    Note that the algorithm will attempt to return a solution with
c    smallest Euclidean norm.
c
      implicit none

      integer nmax
      parameter ( nmax = 100 )

      integer a
      integer a_copy
      integer a_mag
      integer a_sign
      integer b
      integer b_copy
      integer b_mag
      integer b_sign
      integer c
      integer c_copy
      integer g
      integer i4_gcd
      integer ierror
      integer k
      integer n
      integer q(nmax)
      logical swap
      integer x
      integer y
c
c  Defaults for output parameters.
c
      ierror = 0
      x = 0
      y = 0
c
c  Special cases.
c
      if ( a .eq. 0 .and. b .eq. 0 .and. c .eq. 0 ) then
        x = 0
        y = 0
        return
      else if ( a .eq. 0 .and. b .eq. 0 .and. c .ne. 0 ) then
        ierror = 1
        x = 0
        y = 0
        return
      else if ( a .eq. 0 .and. b .ne. 0 .and. c .eq. 0 ) then
        x = 0
        y = 0
        return
      else if ( a .eq. 0 .and. b .ne. 0 .and. c .ne. 0 ) then
        x = 0
        y = c / b
        if ( mod ( c, b ) .ne. 0 ) then
          ierror = 2
        end if
        return
      else if ( a .ne. 0 .and. b .eq. 0 .and. c .eq. 0 ) then
        x = 0
        y = 0
        return
      else if ( a .ne. 0 .and. b .eq. 0 .and. c .ne. 0 ) then
        x = c / a
        y = 0
        if ( mod ( c, a ) .ne. 0 ) then
          ierror = 3
        end if
        return
      else if ( a .ne. 0 .and. b .ne. 0 .and. c .eq. 0 ) then
        g = i4_gcd ( a, b )
        x = b / g
        y = -a / g
        return
      end if
c
c  Now handle the "general" case: A, B and C are nonzero.
c
c  Step 1: Compute the GCD of A and B, which must also divide C.
c
      g = i4_gcd ( a, b )

      if ( mod ( c, g ) .ne. 0 ) then
        ierror = 4
        return
      end if

      a_copy = a / g
      b_copy = b / g
      c_copy = c / g
c
c  Step 2: Split A and B into sign and magnitude.
c
      a_mag = abs ( a_copy )
      a_sign = sign ( 1, a_copy )
      b_mag = abs ( b_copy )
      b_sign = sign ( 1, b_copy )
c
c  Another special case, A_MAG = 1 or B_MAG = 1.
c
      if ( a_mag .eq. 1 ) then
        x = a_sign * c_copy
        y = 0
        return
      else if ( b_mag .eq. 1 ) then
        x = 0
        y = b_sign * c_copy
        return
      end if
c
c  Step 3: Produce the Euclidean remainder sequence.
c
      if ( b_mag .le. a_mag ) then

        swap = .false.
        q(1) = a_mag
        q(2) = b_mag

      else

        swap = .true.
        q(1) = b_mag
        q(2) = a_mag

      end if

      n = 3

10    continue

        q(n) = mod ( q(n-2), q(n-1) )

        if ( q(n) .eq. 1 ) then
          go to 20
        end if

        n = n + 1

        if ( nmax .lt. n ) then
          ierror = 5
          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) 'DIOPHANTINE - Fatal error!'
          write ( *, '(a)' ) '  Exceeded number of iterations.'
          stop
        end if

      go to 10

20    continue
c
c  Step 4: Now go backwards to solve X * A_MAG + Y * B_MAG = 1.
c
      y = 0
      do k = n, 2, -1
        x = y
        y = ( 1 - x * q(k-1) ) / q(k)
      end do
c
c  Step 5: Undo the swapping.
c
      if ( swap ) then
        call i4_swap ( x, y )
      end if
c
c  Step 6: Now apply signs to X and Y so that X * A + Y * B = 1.
c
      x = x * a_sign
      y = y * b_sign
c
c  Step 7: Multiply by C, so that X * A + Y * B = C.
c
      x = x * c_copy
      y = y * c_copy
c
c  Step 8: Given a solution (X,Y), try to find the solution of
c  minimal magnitude.
c
      call diophantine_solution_minimize ( a_copy, b_copy, x, y )

      return
      end
      subroutine diophantine_solution_minimize ( a, b, x, y )

c*********************************************************************72
c
cc DIOPHANTINE_SOLUTION_MINIMIZE seeks a minimal solution of a Diophantine equation.
c
c  Discussion:
c
c    Given a solution (X,Y) of a Diophantine equation:
c
c      A * X + B * Y = C.
c
c    then there are an infinite family of solutions of the form
c
c      ( X(i), Y(i) ) = ( X + i * B, Y - i * A )
c
c    An integral solution of minimal Euclidean norm can be found by
c    tentatively moving along the vectors (B,-A) and (-B,A) one step
c    at a time.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    17 January 2007
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Eric Weisstein, editor,
c    CRC Concise Encylopedia of Mathematics,
c    CRC Press, 1998, page 446.
c
c  Parameters:
c
c    Input, integer A, B, the coefficients of the Diophantine equation.
c    A and B are assumed to be relatively prime.
c
c    Input/output, integer X, Y, on input, a solution of the Diophantine
c    equation.  On output, a solution of minimal Euclidean norm.
c
      implicit none

      integer a
      integer b
      double precision norm
      double precision norm_new
      double precision t
      integer x
      integer xnew
      integer y
      integer ynew
c
c  Compute the minimum for T real, and then look nearby.
c
      t = ( - dble ( b ) * dble ( x ) 
     &  + dble ( a ) * dble ( y ) ) 
     &  / ( dble ( a )**2 + dble ( b )**2 )

      x = x + nint ( t ) * b
      y = y - nint ( t ) * a
c
c  Now look nearby.
c
      norm = ( dble ( x ) )**2 + ( dble ( y ) )**2

10    continue

        xnew = x + b
        ynew = y - a

        norm_new = ( dble ( xnew ) )**2 + ( dble ( ynew ) )**2

        if ( norm .le. norm_new ) then
          go to 20
        end if

        x = xnew
        y = ynew
        norm = norm_new

      go to 10

20    continue

        xnew = x - b
        ynew = y + a

        norm_new = ( dble ( xnew ) )**2 + ( dble ( ynew ) )**2

        if ( norm .le. norm_new ) then
          go to 30
        end if

        x = xnew
        y = ynew
        norm = norm_new

      go to 20

30    continue

      return
      end
      subroutine dvec_add ( n, dvec1, dvec2, dvec3 )

c*********************************************************************72
c
cc DVEC_ADD adds two (signed) decimal vectors.
c
c  Discussion:
c
c    A DVEC is an integer vector of decimal digits, intended to
c    represent an integer.  DVEC(1) is the units digit, DVEC(N-1)
c    is the coefficient of 10**(N-2), and DVEC(N) contains sign
c    information.  It is 0 if the number is positive, and 9 if
c    the number is negative.
c
c  Example:
c
c    N = 4
c
c      DVEC1     +   DVEC2     =   DVEC3
c
c    ( 0 0 1 7 ) + ( 0 1 0 4 ) = ( 0 0 1 2 1 )
c
c          17    +       104   =         121
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    18 January 2007
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the length of the vectors.
c
c    Input, integer DVEC1(N), DVEC2(N), the vectors to be added.
c
c    Output, integer DVEC3(N), the sum of the two input vectors.
c
      implicit none

      integer n

      integer base
      parameter ( base = 10 )
      integer dvec1(n)
      integer dvec2(n)
      integer dvec3(n)
      integer i
      logical overflow

      overflow = .false.

      do i = 1, n
        dvec3(i) = dvec1(i) + dvec2(i)
      end do

      do i = 1, n

10      continue

        if ( base .le. dvec3(i) ) then

          dvec3(i) = dvec3(i) - base

          if ( i .lt. n ) then
            dvec3(i+1) = dvec3(i+1) + 1
          else
            overflow = .true.
          end if

          go to 10

        end if

      end do

      return
      end
      subroutine dvec_complementx ( n, dvec1, dvec2 )

c*********************************************************************72
c
cc DVEC_COMPLEMENTX computes the ten's complement of a decimal vector.
c
c  Discussion:
c
c    A DVEC is an integer vector of decimal digits, intended to
c    represent an integer.  DVEC(1) is the units digit, DVEC(N-1)
c    is the coefficient of 10**(N-2), and DVEC(N) contains sign
c    information.  It is 0 if the number is positive, and 9 if
c    the number is negative.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    18 January 2007
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the length of the vectors.
c
c    Input, integer DVEC1(N), the vector to be complemented.
c
c    Output, integer DVEC2(N), the complemented vector.
c
      implicit none

      integer n

      integer base
      parameter ( base = 10 )
      integer dvec1(n)
      integer dvec2(n)
      integer dvec3(n)
      integer dvec4(n)
      integer i

      do i = 1, n
        dvec3(i) = ( base - 1 ) - dvec1(i)
      end do

      dvec4(1) = 1
      do i = 2, n
        dvec4(i) = 0
      end do

      call dvec_add ( n, dvec3, dvec4, dvec2 )

      return
      end
      subroutine dvec_mul ( n, dvec1, dvec2, dvec3 )

c*********************************************************************72
c
cc DVEC_MUL computes the product of two decimal vectors.
c
c  Discussion:
c
c    A DVEC is an integer vector of decimal digits, intended to
c    represent an integer.  DVEC(1) is the units digit, DVEC(N-1)
c    is the coefficient of 10**(N-2), and DVEC(N) contains sign
c    information.  It is 0 if the number is positive, and 9 if
c    the number is negative.
c
c    Since the user may want to make calls like
c
c      call dvec_mul ( n, dvec1, dvec1, dvec3 )
c    or even
c      call dvec_mul ( n, dvec1, dvec1, dvec1 )
c
c    we need to copy the arguments, work on them, and then copy out the result.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    18 January 2007
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the length of the vectors.
c
c    Input, integer DVEC1(N), DVEC2(N), the vectors to be multiplied.
c
c    Output, integer DVEC3(N), the product of the two input vectors.
c
      implicit none

      integer n

      integer base
      parameter ( base = 10 )
      integer carry
      integer dvec1(n)
      integer dvec2(n)
      integer dvec3(n)
      integer dveca(n)
      integer dvecb(n)
      integer dvecc(n)
      integer i
      integer j
      integer product_sign
c
c  Copy the input.
c
      call i4vec_copy ( n, dvec1, dveca )
      call i4vec_copy ( n, dvec2, dvecb )
c
c  Record the sign of the product.
c  Make the factors positive.
c
      product_sign = 1

      if ( dveca(n) .ne. 0 ) then
        product_sign = - product_sign
        call dvec_complementx ( n, dveca, dveca )
      end if

      if ( dvecb(n) .ne. 0 ) then
        product_sign = - product_sign
        call dvec_complementx ( n, dvecb, dvecb )
      end if

      call i4vec_zero ( n, dvecc )
c
c  Multiply.
c
      do i = 1, n-1
        do j = 1, n-i
          dvecc(j+i-1) = dvecc(j+i-1) + dveca(i) * dvecb(j)
        end do
      end do
c
c  Take care of carries.
c  Unlike the DVEC_ADD routine, we do NOT allow carries into the
c  N-th position.
c
      do i = 1, n-1

        carry = dvecc(i) / base
        dvecc(i) = dvecc(i) - carry * base

        if ( i .lt. n - 1 ) then
          dvecc(i+1) = dvecc(i+1) + carry
        end if

      end do
c
c  Take care of the sign of the product.
c
      if ( product_sign .lt. 0 ) then
        call dvec_complementx ( n, dvecc, dvecc )
      end if
c
c  Copy the output.
c
      call i4vec_copy ( n, dvecc, dvec3 )

      return
      end
      subroutine dvec_print ( n, dvec, title )

c*********************************************************************72
c
cc DVEC_PRINT prints a decimal integer vector, with an optional title.
c
c  Discussion:
c
c    A DVEC is an integer vector of decimal digits, intended to
c    represent an integer.  DVEC(1) is the units digit, DVEC(N-1)
c    is the coefficient of 10**(N-2), and DVEC(N) contains sign
c    information.  It is 0 if the number is positive, and 9 if
c    the number is negative.
c
c    The vector is printed "backwards", that is, the first entry
c    printed is DVEC(N).
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    18 January 2007
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the number of components of the vector.
c
c    Input, integer DVEC(N), the vector to be printed.
c
c    Input, character * ( * ) TITLE, a title to be printed first.
c    TITLE may be blank.
c
      implicit none

      integer n

      integer dvec(n)
      integer i
      integer ihi
      integer ilo
      integer s_len_trim
      character * ( * ) title
      integer title_length

      title_length = s_len_trim ( title )

      if ( 0 .lt. title_length ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) title(1:title_length)
        write ( *, '(a)' ) ' '
      end if

      do ihi = n, 1, -80
        ilo = max ( ihi - 79, 1 )
        write ( *, '(2x,80i1)' ) ( dvec(i), i = ihi, ilo, -1 )
      end do

      return
      end
      subroutine dvec_sub ( n, dvec1, dvec2, dvec3 )

c*********************************************************************72
c
cc DVEC_SUB subtracts two decimal vectors.
c
c  Discussion:
c
c    A DVEC is an integer vector of decimal digits, intended to
c    represent an integer.  DVEC(1) is the units digit, DVEC(N-1)
c    is the coefficient of 10**(N-2), and DVEC(N) contains sign
c    information.  It is 0 if the number is positive, and 9 if
c    the number is negative.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    18 January 2007
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the length of the vectors.
c
c    Input, integer DVEC1(N), DVEC2(N), the vectors to be subtracted.
c
c    Output, integer DVEC3(N), the value of DVEC1 - DVEC2.
c
      implicit none

      integer n

      integer dvec1(n)
      integer dvec2(n)
      integer dvec3(n)

      call i4vec_copy ( n, dvec2, dvec3 )

      call dvec_complementx ( n, dvec3, dvec3 )

      call dvec_add ( n, dvec1, dvec3, dvec3 )

      return
      end
      subroutine dvec_to_i4 ( n, dvec, i4 )

c*********************************************************************72
c
cc DVEC_TO_I4 makes an integer from a (signed) decimal vector.
c
c  Discussion:
c
c    A DVEC is an integer vector of decimal digits, intended to
c    represent an integer.  DVEC(1) is the units digit, DVEC(N-1)
c    is the coefficient of 10**(N-2), and DVEC(N) contains sign
c    information.  It is 0 if the number is positive, and 9 if
c    the number is negative.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    18 January 2007
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the dimension of the vector.
c
c    Input, integer DVEC(N), the decimal vector.
c
c    Output, integer I4, the integer.
c
      implicit none

      integer n

      integer base
      parameter ( base = 10 )
      integer dvec(n)
      integer dvec2(n)
      integer i
      integer i_sign
      integer i4

      call i4vec_copy ( n, dvec, dvec2 )

      i_sign = 1

      if ( dvec2(n) .eq. base - 1 ) then
        i_sign = -1
        call dvec_complementx ( n-1, dvec2, dvec2 )
      end if

      i4 = 0
      do i = n-1, 1, -1
        i4 = base * i4 + dvec2(i)
      end do

      i4 = i_sign * i4

      return
      end
      subroutine equiv_next ( n, npart, jarray, iarray, more )

c*********************************************************************72
c
cc EQUIV_NEXT computes the partitions of a set one at a time.
c
c  Discussion:
c
c    A partition of a set assigns each element to exactly one subset.
c
c    The number of partitions of a set of size N is the Bell number B(N).
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    19 January 2007
c
c  Author:
c
c    Original FORTRAN77 version by Albert Nijenhuis, Herbert Wilf.
c    This FORTRAN77 version by John Burkardt.
c
c  Reference:
c
c    Albert Nijenhuis, Herbert Wilf,
c    Combinatorial Algorithms for Computers and Calculators,
c    Second Edition,
c    Academic Press, 1978,
c    ISBN: 0-12-519260-6,
c    LC: QA164.N54.
c
c  Parameters:
c
c    Input, integer N, the number of elements in the set to be partitioned.
c
c    Input/output, integer NPART, the number of subsets in the partition.
c
c    Input/output, integer JARRAY(N), the number of elements
c    in each subset of the partition.
c
c    Input/output, integer IARRAY(N), the subset to which each element belongs.
c
c    Input/output, logical MORE.  Set MORE = FALSE before first call.
c    It is reset and held at TRUE as long as
c    the partition returned is not the last one.
c    When MORE is returned FALSE, all the partitions
c    have been computed and returned.
c
      implicit none

      integer n

      integer i
      integer iarray(n)
      integer jarray(n)
      integer l
      integer m
      logical more
      integer npart

      if ( .not. more ) then

        npart = 1
        do i = 1, n
          iarray(i) = 1
        end do
        jarray(1) = n

      else

        m = n

10      continue

        if ( jarray(iarray(m)) .eq. 1 ) then
          iarray(m) = 1
          m = m - 1
          go to 10
        end if

        l = iarray(m)
        npart = npart + m - n
        jarray(1) = jarray(1) + n - m

        if ( l .eq. npart ) then
          npart = npart + 1
          jarray(npart) = 0
        end if

        iarray(m) = l + 1
        jarray(l) = jarray(l) - 1
        jarray(l+1) = jarray(l+1) + 1

      end if

      more = ( npart .ne. n )

      return
      end
      subroutine equiv_next2 ( n, a, done )

c*********************************************************************72
c
cc EQUIV_NEXT2 computes, one at a time, the partitions of a set.
c
c  Discussion:
c
c    A partition of a set assigns each element to exactly one subset.
c
c    The number of partitions of a set of size N is the Bell number B(N).
c
c    The entries of A are the partition subset to which each
c    element of the original set belongs.  If there are NPART distinct
c    parts of the partition, then each entry of A will be a
c    number between 1 and NPART.  Every number from 1 to NPART will
c    occur somewhere in the list.  If the entries of A are
c    examined in order, then each time a new partition subset occurs,
c    it will be the next unused integer.
c
c    For instance, for N = 4, the program will describe the set
c    where each element is in a separate subset as 1, 2, 3, 4,
c    even though such a partition might also be described as
c    4, 3, 2, 1 or even 1, 5, 8, 19.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    19 January 2007
c
c  Parameters:
c
c    Input, integer N, the number of elements in the set.
c
c    Input/output, integer A(N), contains the information
c    defining the current partition.  The user should not alter
c    A between calls.  Except for the very first
c    call, the routine uses the previous output value of A to compute
c    the next value.
c
c    Input/output, logical DONE.  Before the very first call, the
c    user should set DONE to TRUE, which prompts the program
c    to initialize its data, and return the first partition.
c    Thereafter, the user should call again, for the next
c    partition, and so on, until the routine returns with DONE
c    equal to TRUE, at which point there are no more partitions
c    to compute.
c
      implicit none

      integer n

      integer a(n)
      logical done
      integer i
      integer imax
      integer j
      integer jmax

      if ( done ) then

        done = .false.

        do i = 1, n
          a(i) = 1
        end do

      else
c
c  Find the last element J that can be increased by 1.
c  This is the element that is not equal to its maximum possible value,
c  which is the maximum value of all preceding elements +1.
c
        jmax = a(1)
        imax = 1

        do j = 2, n

          if ( jmax .lt. a(j) ) then
            jmax = a(j)
          else
            imax = j
          end if

        end do
c
c  If no element can be increased by 1, we are done.
c
        if ( imax .eq. 1 ) then
          done = .true.
          return
        end if
c
c  Increase the value of the IMAX-th element by 1, set its successors to 1.
c
        done = .false.
        a(imax) = a(imax) + 1
        do i = imax+1, n
          a(i) = 1
        end do

      end if

      return
      end
      subroutine equiv_print ( n, a, title )

c*********************************************************************72
c
cc EQUIV_PRINT prints a partition of a set.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    19 January 2007
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, number of elements in set to be partitioned.
c
c    Input, integer A(N), defines the partition or set of equivalence
c    classes.  Element I belongs to subset A(I).
c
c    Input, character * ( * ) TITLE, a title to be printed first.
c    TITLE may be blank.
c
      implicit none

      integer n

      integer a(n)
      integer i
      integer j
      integer k
      integer karray(n)
      integer s
      integer s_max
      integer s_min
      character * ( * ) title

      if ( title .ne. ' ' ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) title
      end if

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '   Set  Size'

      call i4vec_min ( n, a, s_min )
      call i4vec_max ( n, a, s_max )

      do s = s_min, s_max

        k = 0

        do j = 1, n

          if ( a(j) .eq. s ) then
            k = k + 1
            karray(k) = j
          end if

        end do

        if ( 0 .lt. k ) then
          write ( *, '(i8,i8,a,(10i4))' ) 
     &    s, k, ' :: ', ( karray(i), i = 1, k )
        end if

      end do

      return
      end
      subroutine equiv_random ( n, seed, npart, a, b )

c*********************************************************************72
c
cc EQUIV_RANDOM selects a random partition of a set.
c
c  Discussion:
c
c    The user does not control the number of parts in the partition.
c
c    The equivalence classes are numbered in no particular order.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    19 January 2007
c
c  Author:
c
c    Original FORTRAN77 version by Albert Nijenhuis, Herbert Wilf.
c    This FORTRAN77 version by John Burkardt.
c
c  Reference:
c
c    Albert Nijenhuis, Herbert Wilf,
c    Combinatorial Algorithms for Computers and Calculators,
c    Second Edition,
c    Academic Press, 1978,
c    ISBN: 0-12-519260-6,
c    LC: QA164.N54.
c
c  Parameters:
c
c    Input, integer N, the number of elements in the set to be partitioned.
c
c    Input/output, integer SEED, a seed for the random number generator.
c
c    Output, integer NPART, the number of classes or parts in the 
c    partition.  NPART will be between 1 and N.
c
c    Output, integer A(N), indicates the class to which each element
c    is assigned.
c
c    Output, double precision B(N).  B(K) = C(K)/(Kc), where
c    C(K) = number of partitions of a set of K objects.
c
      implicit none

      integer n

      integer a(n)
      double precision b(n)
      integer k
      integer l
      integer m
      integer npart
      double precision r8_uniform_01
      integer seed
      double precision sum1
      double precision z

      b(1) = 1.0D+00

      do l = 1, n-1

        sum1 = 1.0D+00 / dble ( l )
        do k = 1, l-1
          sum1 = ( sum1 + b(k) ) / dble ( l - k )
        end do

        b(l+1) = ( sum1 + b(l) ) / dble ( l + 1 )

      end do

      m = n
      npart = 0

10    continue

        z = r8_uniform_01 ( seed )
        z = dble ( m ) * b(m) * z
        k = 0
        npart = npart + 1

20      continue

        if ( 0.0D+00 .le. z ) then

          a(m) = npart
          m = m - 1

          if ( m .eq. 0 ) then
            go to 30
          end if

          z = z - b(m)
          k = k + 1
          z = z * k

          go to 20

        end if

30      continue

        if ( m .eq. 0 ) then
          go to 40
        end if

      go to 10
c
c  Randomly permute the assignments.
c
40    continue

      call perm_random2 ( n, seed, a )

      return
      end
      subroutine euler ( n, ieuler )

c*********************************************************************72
c
cc EULER returns the N-th row of Euler's triangle.
c
c  Discussion:
c
c    E(N,K) counts the number of permutations of the N digits that have
c    exactly K "ascents", that is, K places where the Ith digit is
c    less than the (I+1)th digit.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    19 January 2007
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the row of Euler's triangle desired.
c
c    Output, integer IEULER(0:N), the N-th row of Euler's
c    triangle, IEULER(K) contains the value of E(N,K).  Note
c    that IEULER(0) should be 1 and IEULER(N) should be 0.
c
      implicit none

      integer n

      integer ieuler(0:n)
      integer irow
      integer k

      ieuler(0) = 1

      if ( 0 .lt. n ) then

        ieuler(1) = 0

        do irow = 2, n

          ieuler(irow) = 0

          do k = irow-1, 1, -1

            ieuler(k) = ( k + 1 ) * ieuler(k) 
     &      + ( irow - k ) * ieuler(k-1)

          end do

          ieuler(0) = 1

        end do

      end if

      return
      end
      function frobenius_number_order2 ( c1, c2 )

c*********************************************************************72
c
cc FROBENIUS_NUMBER_ORDER2 returns the Frobenius number for order 2.
c
c  Discussion:
c
c    The Frobenius number of order N is the solution of the Frobenius
c    coin sum problem for N coin denominations.
c
c    The Frobenius coin sum problem assumes the existence of
c    N coin denominations, and asks for the largest value that cannot
c    be formed by any combination of coins of these denominations.
c
c    The coin denominations are assumed to be distinct positive integers.
c
c    For general N, this problem is fairly difficult to handle.
c
c    For N = 2, it is known that:
c
c    * if C1 and C2 are not relatively prime, then
c      there are infinitely large values that cannot be formed.
c
c    * otherwise, the largest value that cannot be formed is
c      C1 * C2 - C1 - C2, and that exactly half the values between
c      1 and C1 * C2 - C1 - C2 + 1 cannot be represented.
c
c    As a simple example, if C1 = 2 and C2 = 7, then the largest
c    unrepresentable value is 5, and there are (5+1)/2 = 3
c    unrepresentable values, namely 1, 3, and 5.
c
c    For a general N, and a set of coin denominations C1, C2, ..., CN,
c    the Frobenius number F(N, C(1:N) ) is defined as the largest value
c    B for which the equation
c
c      C1*X1 + C2*X2 + ... + CN*XN = B
c
c    has no nonnegative integer solution X(1:N).
c
c    In the Mathematica Package "NumberTheory", the Frobenius number
c    can be determined by
c
c    <<NumberTheory`Frobenius`
c    FrobeniusF[ {C1,...,CN} ]
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    19 January 2007
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    James Sylvester,
c    Question 7382,
c    Mathematical Questions with their Solutions,
c    Educational Times,
c    Volume 41, page 21, 1884.
c
c    Stephen Wolfram,
c    The Mathematica Book,
c    Fourth Edition,
c    Cambridge University Press, 1999,
c    ISBN: 0-521-64314-7,
c    LC: QA76.95.W65.
c
c  Parameters:
c
c    Input, integer C1, C2, the coin denominations. C1 and C2
c    should be positive and relatively prime.
c
c    Output, integer FROBENIUS_NUMBER_ORDER2, the Frobenius number of (C1,C2).
c
      implicit none

      integer c1
      integer c2
      integer frobenius_number_order2
      integer i4_gcd
      integer i4_huge
      integer value

      if ( c1 .le. 0 ) then
        value = i4_huge ( )
      else if ( c2 .le. 0 ) then
        value = i4_huge ( )
      else if ( i4_gcd ( c1, c2 ) .ne. 1 ) then
        value = i4_huge ( )
      else
        value = c1 * c2 - c1 - c2
      end if

      frobenius_number_order2 = value

      return
      end
      subroutine frobenius_number_order2_values ( n_data, c1, c2, f )

c*********************************************************************72
c
cc FROBENIUS_NUMBER_ORDER2_VALUES returns values of the order 2 Frobenius number.
c
c  Discussion:
c
c    The Frobenius number of order N is the solution of the Frobenius
c    coin sum problem for N coin denominations.
c
c    The Frobenius coin sum problem assumes the existence of
c    N coin denominations, and asks for the largest value that cannot
c    be formed by any combination of coins of these denominations.
c
c    The coin denominations are assumed to be distinct positive integers.
c
c    For general N, this problem is fairly difficult to handle.
c
c    For N = 2, it is known that:
c
c    * if C1 and C2 are not relatively prime, then
c      there are infinitely large values that cannot be formed.
c
c    * otherwise, the largest value that cannot be formed is
c      C1 * C2 - C1 - C2, and that exactly half the values between
c      1 and C1 * C2 - C1 - C2 + 1 cannot be represented.
c
c    As a simple example, if C1 = 2 and C2 = 7, then the largest
c    unrepresentable value is 5, and there are (5+1)/2 = 3
c    unrepresentable values, namely 1, 3, and 5.
c
c    For a general N, and a set of coin denominations C1, C2, ..., CN,
c    the Frobenius number F(N, C(1:N) ) is defined as the largest value
c    B for which the equation
c
c      C1*X1 + C2*X2 + ... + CN*XN = B
c
c    has no nonnegative integer solution X(1:N).
c
c    In the Mathematica Package "NumberTheory", the Frobenius number
c    can be determined by
c
c    <<NumberTheory`Frobenius`
c    FrobeniusF[ {C1,...,CN} ]
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    20 January 2007
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    James Sylvester,
c    Question 7382,
c    Mathematical Questions with their Solutions,
c    Educational Times,
c    Volume 41, page 21, 1884.
c
c    Stephen Wolfram,
c    The Mathematica Book,
c    Fourth Edition,
c    Cambridge University Press, 1999,
c    ISBN: 0-521-64314-7,
c    LC: QA76.95.W65.
c
c  Parameters:
c
c    Input/output, integer N_DATA.  The user sets N_DATA to 0 before the
c    first call.  On each call, the routine increments N_DATA by 1, and
c    returns the corresponding data; when there is no more data, the
c    output value of N_DATA will be 0 again.
c
c    Output, integer C1, C2, the parameters of the function.
c
c    Output, integer F, the value of the function.
c
      implicit none

      integer n_max
      parameter ( n_max = 6 )

      integer c1
      integer c1_vec(n_max)
      integer c2
      integer c2_vec(n_max)
      integer f
      integer f_vec(n_max)
      integer n_data

      save c1_vec
      save c2_vec
      save f_vec

      data c1_vec /
     &   2,
     &   3,
     &   4,
     &   5,
     &  12,
     &  99 /
      data c2_vec /
     &    5,
     &   17,
     &   19,
     &   13,
     &   11,
     &  100 /
      data f_vec /
     &     3,
     &    31,
     &    23,
     &    47,
     &   109,
     &  9701 /

      if ( n_data .lt. 0 ) then
        n_data = 0
      end if

      n_data = n_data + 1

      if ( n_max .lt. n_data ) then
        n_data = 0
        c1 = 0
        c2 = 0
        f = 0
      else
        c1 = c1_vec(n_data)
        c2 = c2_vec(n_data)
        f = f_vec(n_data)
      end if

      return
      end
      function gamma_log ( x )

c*********************************************************************72
c
cc GAMMA_LOG calculates the natural logarithm of GAMMA ( X ) for positive X.
c
c  Discussion:
c
c    The program uses rational functions that theoretically approximate
c    LOG(GAMMA(X)) to at least 18 significant decimal digits.  The
c    approximation for 12 < X is from Hart et al, while approximations
c    for X < 12.0 are similar to those in Cody and Hillstrom, but are
c    unpublished.  The accuracy achieved depends on the arithmetic system,
c    the compiler, intrinsic functions, and proper selection of the
c    machine-dependent constants.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    16 June 1999
c
c  Author:
c
c    Original FORTRAN77 version by William Cody, Laura Stoltz.
c    This FORTRAN77 version by John Burkardt.
c
c  Reference:
c
c    William Cody, Kenneth Hillstrom,
c    Chebyshev Approximations for the Natural Logarithm of the Gamma Function,
c    Mathematics of Computation,
c    Volume 21, 1967, pages 198-203.
c
c    Kenneth Hillstrom,
c    ANL/AMD Program ANLC366S, DGAMMA/DLGAMA,
c    May 1969.
c
c    John Hart, Ward Cheney, Charles Lawson, Hans Maehly, Charles Mesztenyi,
c    John Rice, Henry Thatcher, Christoph Witzgall,
c    Computer Approximations,
c    Wiley, 1968.
c
c  Parameters:
c
c    Input, double precision X, the argument of the Gamma function.
c    X must be positive.
c
c    Output, double precision GAMMA_LOG, the logarithm of the Gamma
c    function of X.  If X <= 0.0, or if overflow would occur, the program
c    returns the value XINF, the largest representable floating point number.
c
c  Local Parameters:
c
c    BETA   - radix for the floating-point representation.
c
c    MAXEXP - the smallest positive power of BETA that overflows.
c
c    XBIG   - largest argument for which LN(GAMMA(X)) is representable
c           in the machine, i.e., the solution to the equation
c             LN(GAMMA(XBIG)) = BETA**MAXEXP.
c
c    FRTBIG - Rough estimate of the fourth root of XBIG
c
c
c  Approximate values for some important machines are:
c
c                            BETA      MAXEXP         XBIG
c
c  CRAY-1        (S.P.)        2        8191       9.62D+2461
c  Cyber 180/855
c    under NOS   (S.P.)        2        1070       1.72D+319
c  IEEE (IBM/XT,
c    SUN, etc.)  (S.P.)        2         128       4.08D+36
c  IEEE (IBM/XT,
c    SUN, etc.)  (D.P.)        2        1024       2.55D+305
c  IBM 3033      (D.P.)       16          63       4.29D+73
c  VAX D-Format  (D.P.)        2         127       2.05D+36
c  VAX G-Format  (D.P.)        2        1023       1.28D+305
c
c
c                           FRTBIG
c
c  CRAY-1        (S.P.)   3.13D+615
c  Cyber 180/855
c    under NOS   (S.P.)   6.44D+79
c  IEEE (IBM/XT,
c    SUN, etc.)  (S.P.)   1.42D+9
c  IEEE (IBM/XT,
c    SUN, etc.)  (D.P.)   2.25D+76
c  IBM 3033      (D.P.)   2.56D+18
c  VAX D-Format  (D.P.)   1.20D+9
c  VAX G-Format  (D.P.)   1.89D+76
c
      implicit none

      double precision c(7)
      double precision corr
      double precision d1
      double precision d2
      double precision d4
      double precision frtbig
      integer i
      double precision gamma_log
      double precision p1(8)
      double precision p2(8)
      double precision p4(8)
      double precision pnt68
      double precision q1(8)
      double precision q2(8)
      double precision q4(8)
      double precision r8_epsilon
      double precision r8_huge
      double precision res
      double precision sqrtpi
      double precision x
      double precision xbig
      double precision xden
      double precision xm1
      double precision xm2
      double precision xm4
      double precision xnum
      double precision xsq

      data c /
     &  -1.910444077728D-03, 
     &   8.4171387781295D-04, 
     &  -5.952379913043012D-04, 
     &   7.93650793500350248D-04, 
     &  -2.777777777777681622553D-03, 
     &   8.333333333333333331554247D-02, 
     &   5.7083835261D-03 /
      data d1 /  -5.772156649015328605195174D-01 /
      data d2 /  4.227843350984671393993777D-01 /
      data d4 /  1.791759469228055000094023D+00 /
      data frtbig / 1.42D+09 /
      data p1 /
     &   4.945235359296727046734888D+00, 
     &   2.018112620856775083915565D+02, 
     &   2.290838373831346393026739D+03, 
     &   1.131967205903380828685045D+04, 
     &   2.855724635671635335736389D+04, 
     &   3.848496228443793359990269D+04, 
     &   2.637748787624195437963534D+04, 
     &   7.225813979700288197698961D+03 /
      data p2 /
     &   4.974607845568932035012064D+00, 
     &   5.424138599891070494101986D+02, 
     &   1.550693864978364947665077D+04, 
     &   1.847932904445632425417223D+05, 
     &   1.088204769468828767498470D+06, 
     &   3.338152967987029735917223D+06, 
     &   5.106661678927352456275255D+06, 
     &   3.074109054850539556250927D+06 /
      data p4 /
     &   1.474502166059939948905062D+04, 
     &   2.426813369486704502836312D+06, 
     &   1.214755574045093227939592D+08, 
     &   2.663432449630976949898078D+09, 
     &   2.940378956634553899906876D+010, 
     &   1.702665737765398868392998D+011, 
     &   4.926125793377430887588120D+011, 
     &   5.606251856223951465078242D+011 /
      data pnt68 / 0.6796875D+00 /
      data q1 / 
     &   6.748212550303777196073036D+01, 
     &   1.113332393857199323513008D+03, 
     &   7.738757056935398733233834D+03, 
     &   2.763987074403340708898585D+04, 
     &   5.499310206226157329794414D+04, 
     &   6.161122180066002127833352D+04, 
     &   3.635127591501940507276287D+04, 
     &   8.785536302431013170870835D+03 /
      data q2 /
     &   1.830328399370592604055942D+02, 
     &   7.765049321445005871323047D+03, 
     &   1.331903827966074194402448D+05, 
     &   1.136705821321969608938755D+06, 
     &   5.267964117437946917577538D+06, 
     &   1.346701454311101692290052D+07, 
     &   1.782736530353274213975932D+07, 
     &   9.533095591844353613395747D+06 /
      data q4 /
     &   2.690530175870899333379843D+03, 
     &   6.393885654300092398984238D+05, 
     &   4.135599930241388052042842D+07, 
     &   1.120872109616147941376570D+09, 
     &   1.488613728678813811542398D+010, 
     &   1.016803586272438228077304D+011, 
     &   3.417476345507377132798597D+011, 
     &   4.463158187419713286462081D+011 /
      data sqrtpi / 0.9189385332046727417803297D+00 /
      data xbig / 4.08D+36 /
c
c  Return immediately if the argument is out of range.
c
      if ( x .le. 0.0D+00 .or. xbig .lt. x ) then
        gamma_log = r8_huge ( )
        return
      end if

      if ( x .le. r8_epsilon ( ) ) then

        res = -log ( x )

      else if ( x .le. 1.5D+00 ) then

        if ( x .lt. pnt68 ) then
          corr = -log ( x )
          xm1 = x
        else
          corr = 0.0D+00
          xm1 = ( x - 0.5D+00 ) - 0.5D+00
        end if

        if ( x .le. 0.5D+00 .or. pnt68 .le. x ) then

          xden = 1.0D+00
          xnum = 0.0D+00

          do i = 1, 8
            xnum = xnum * xm1 + p1(i)
            xden = xden * xm1 + q1(i)
          end do

          res = corr + ( xm1 * ( d1 + xm1 * ( xnum / xden ) ) )

        else

          xm2 = ( x - 0.5D+00 ) - 0.5D+00
          xden = 1.0D+00
          xnum = 0.0D+00
          do i = 1, 8
            xnum = xnum * xm2 + p2(i)
            xden = xden * xm2 + q2(i)
          end do

          res = corr + xm2 * ( d2 + xm2 * ( xnum / xden ) )

        end if

      else if ( x .le. 4.0D+00 ) then

        xm2 = x - 2.0D+00
        xden = 1.0D+00
        xnum = 0.0D+00
        do i = 1, 8
          xnum = xnum * xm2 + p2(i)
          xden = xden * xm2 + q2(i)
        end do

        res = xm2 * ( d2 + xm2 * ( xnum / xden ) )

      else if ( x .le. 12.0D+00 ) then

        xm4 = x - 4.0D+00
        xden = -1.0D+00
        xnum = 0.0D+00
        do i = 1, 8
          xnum = xnum * xm4 + p4(i)
          xden = xden * xm4 + q4(i)
        end do

        res = d4 + xm4 * ( xnum / xden )

      else

        res = 0.0D+00

        if ( x .le. frtbig ) then

          res = c(7)
          xsq = x * x

          do i = 1, 6
            res = res / xsq + c(i)
          end do

        end if

        res = res / x
        corr = log ( x )
        res = res + sqrtpi - 0.5D+00 * corr
        res = res + x * ( corr - 1.0D+00 )

      end if

      gamma_log = res

      return
      end
      subroutine gamma_log_values ( n_data, x, fx )

c*********************************************************************72
c
cc GAMMA_LOG_VALUES returns some values of the Log Gamma function.
c
c  Discussion:
c
c    In Mathematica, the function can be evaluated by:
c
c      Log[Gamma[x]]
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    19 January 2007
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Milton Abramowitz, Irene Stegun,
c    Handbook of Mathematical Functions,
c    National Bureau of Standards, 1964,
c    ISBN: 0-486-61272-4,
c    LC: QA47.A34.
c
c    Stephen Wolfram,
c    The Mathematica Book,
c    Fourth Edition,
c    Wolfram Media / Cambridge University Press, 1999.
c
c  Parameters:
c
c    Input/output, integer N_DATA.  The user sets N_DATA to 0 before the
c    first call.  On each call, the routine increments N_DATA by 1, and
c    returns the corresponding data; when there is no more data, the
c    output value of N_DATA will be 0 again.
c
c    Output, double precision X, the argument of the function.
c
c    Output, double precision FX, the value of the function.
c
      implicit none

      integer n_max
      parameter ( n_max = 20 )

      double precision fx
      double precision fx_vec(n_max)
      integer n_data
      double precision x
      double precision x_vec(n_max)

      save fx_vec
      save x_vec

      data fx_vec /
     &   0.1524063822430784D+01, 
     &   0.7966778177017837D+00, 
     &   0.3982338580692348D+00, 
     &   0.1520596783998375D+00, 
     &   0.0000000000000000D+00, 
     &  -0.4987244125983972D-01, 
     &  -0.8537409000331584D-01, 
     &  -0.1081748095078604D+00, 
     &  -0.1196129141723712D+00, 
     &  -0.1207822376352452D+00, 
     &  -0.1125917656967557D+00, 
     &  -0.9580769740706586D-01, 
     &  -0.7108387291437216D-01, 
     &  -0.3898427592308333D-01, 
     &  0.00000000000000000D+00, 
     &  0.69314718055994530D+00, 
     &  0.17917594692280550D+01, 
     &  0.12801827480081469D+02, 
     &  0.39339884187199494D+02, 
     &  0.71257038967168009D+02 /
      data x_vec /
     &  0.20D+00, 
     &  0.40D+00, 
     &  0.60D+00, 
     &  0.80D+00, 
     &  1.00D+00, 
     &  1.10D+00, 
     &  1.20D+00, 
     &  1.30D+00, 
     &  1.40D+00, 
     &  1.50D+00, 
     &  1.60D+00, 
     &  1.70D+00, 
     &  1.80D+00, 
     &  1.90D+00, 
     &  2.00D+00, 
     &  3.00D+00, 
     &  4.00D+00, 
     & 10.00D+00, 
     & 20.00D+00, 
     & 30.00D+00 /

      if ( n_data .lt. 0 ) then
        n_data = 0
      end if

      n_data = n_data + 1

      if ( n_max .lt. n_data ) then
        n_data = 0
        x = 0.0D+00
        fx = 0.0D+00
      else
        x = x_vec(n_data)
        fx = fx_vec(n_data)
      end if

      return
      end
      subroutine get_seed ( seed )

c*********************************************************************72
c
cc GET_SEED returns a seed for the random number generator.
c
c  Discussion:
c
c    The seed depends on the current time, and ought to be (slightly)
c    different every millisecond.  Thus, calling this routine several
c    times in succession will probably return the SAME seed, but
c    calling it a few minutes or days apart will turn a suitably
c    "random" seed.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    11 January 2007
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, integer SEED, a pseudorandom seed value.
c
      implicit none

      integer day
      integer hour
      integer i4_huge
      integer milli
      integer minute
      integer month
      integer second
      integer seed
      double precision temp
      character * ( 10 ) time
      character * ( 8 ) date
      integer year

      call date_and_time ( date, time )

      read ( date, '(i4,i2,i2)' ) year, month, day
      read ( time, '(i2,i2,i2,1x,i3)' ) hour, minute, second, milli

      temp = 0.0D+00
      temp = temp + dble ( month - 1 ) / 11.0D+00
      temp = temp + dble ( day   - 1 ) / 30.0D+00
      temp = temp + dble ( hour      ) / 23.0D+00
      temp = temp + dble ( minute    ) / 59.0D+00
      temp = temp + dble ( second    ) / 59.0D+00
      temp = temp + dble ( milli     ) / 999.0D+00

      temp = temp / 6.0D+00
c
c  Force 0 < TEMP <= 1.
c
10    continue

      if ( temp .le. 0.0D+00 ) then
        temp = temp + 1.0D+00
        go to 10
      end if

20    continue

      if ( 1.0D+00 .lt. temp ) then
        temp = temp - 1.0D+00
        go to 20
      end if

      seed = int ( dble ( i4_huge ( ) ) * temp )
c
c  Never use a seed of 0 or maximum integer.
c
      if ( seed .eq. 0 ) then
        seed = 1
      end if

      if ( seed .eq. i4_huge ( ) ) then
        seed = seed - 1
      end if

      return
      end
      subroutine gray_next ( n, change, a )

c*********************************************************************72
c
cc GRAY_NEXT generates the next Gray code by switching one item at a time.
c
c  Discussion:
c
c    On the first call only, the user must set CHANGE = -N.
c    This initializes the routine to the Gray code for N zeroes.
c
c    Each time it is called thereafter, it returns in CHANGE the index
c    of the item to be switched in the Gray code.  The sign of CHANGE
c    indicates whether the item is to be added or subtracted (or
c    whether the corresponding bit should become 1 or 0).  When
c    CHANGE is equal to N+1 on output, all the Gray codes have been
c    generated.
c
c    The routine has internal memory that is set up on call with
c    CHANGE = -N, and released on final return.
c
c  Example:
c
c    N  CHANGE         Subset in/out   Binary Number
c                      Interpretation  Interpretation
c                       1 2 4 8
c   --  ---------      --------------  --------------
c
c    4   -4 / 0         0 0 0 0         0
c
c        +1             1 0 0 0         1
c          +2           1 1 0 0         3
c        -1             0 1 0 0         2
c            +3         0 1 1 0         6
c        +1             1 1 1 0         7
c          -2           1 0 1 0         5
c        -1             0 0 1 0         4
c              +4       0 0 1 1        12
c        +1             1 0 1 1        13
c          +2           1 1 1 1        15
c        -1             0 1 1 1        14
c            -3         0 1 0 1        10
c        +1             1 1 0 1        11
c          -2           1 0 0 1         9
c        -1             0 0 0 1         8
c              -4       0 0 0 0         0
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    19 January 2007
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Albert Nijenhuis, Herbert Wilf,
c    Combinatorial Algorithms for Computers and Calculators,
c    Second Edition,
c    Academic Press, 1978,
c    ISBN: 0-12-519260-6,
c    LC: QA164.N54.
c
c  Parameters:
c
c    Input, integer N, the order of the total set from which
c    subsets will be drawn.
c
c    Input/output, integer CHANGE.  This is an input item only
c    on the first call for a particular sequence of Gray codes,
c    at which time it must be set to -N.  This corresponds to
c    all items being excluded, or all bits being 0, in the Gray code.
c    On output, CHANGE indicates which of the N items must be "changed", 
c    and the sign indicates whether the item is to be added or removed
c    (or the bit is to become 1 or 0).  On return from the 
c    first call, CHANGE is set to 0, indicating that the first set
c    is the empty set.
c
c    Workspace, integer A(N).
c
      implicit none

      integer n

      integer a(n)
      integer change
      integer i
      integer k
      integer n_save

      save k
      save n_save

      data k / 0 /
      data n_save / 0 /

      if ( n .le. 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'GRAY_NEXT - Fatal error!'
        write ( *, '(a)' ) '  Input value of N <= 0.'
        stop
      end if

      if ( change .eq. -n ) then

        do i = 1, n
          a(i) = 0
        end do

        k = 1
        change = 0

        return
      end if
c
c  First determine WHICH item is to be changed.
c
      if ( mod ( k, 2 ) .eq. 1 ) then

        change = 1

      else

        do i = 1, n
          if ( a(i) .eq. 1 ) then
            change = i + 1
            go to 10
          end if
        end do

10      continue

      end if
c
c  Take care of the terminal case CHANGE = N + 1.
c
      if ( change .eq. n + 1 ) then
        change = n
      end if
c
c  Now determine HOW the item is to be changed.
c
      if ( a(change) .eq. 0 ) then
        a(change) = 1
      else if ( a(change) .eq. 1 ) then
        a(change) = 0
        change = -change
      end if
c
c  Update the counter.
c
      k = k + 1
c
c  If the output CHANGE = -N_SAVE, then we're done.
c
      if ( change .eq. -n ) then

        k = 0

      end if

      return
      end
      subroutine gray_rank ( gray, rank )

c*****************************************************************************80
c
cc GRAY_RANK ranks a Gray code.
c
c  Discussion:
c
c    Given the number GRAY, its ranking is the order in which it would be
c    visited in the Gray code ordering.  The Gray code ordering begins
c
c    Rank  Gray  Gray
c          (Dec) (Bin)
c
c       0     0  0000
c       1     1  0001
c       2     3  0011
c       3     2  0010
c       4     6  0110
c       5     7  0111
c       6     5  0101
c       7     4  0100
c       8    12  0110
c       etc
c
c   This routine is given a Gray code, and has to return the rank.
c
c  Example:
c
c    Gray  Gray  Rank
c    (Dec) (Bin)
c
c     0       0     0
c     1       1     1
c     2      10     3
c     3      11     2
c     4     100     7
c     5     101     6
c     6     110     4
c     7     111     5
c     8    1000    15
c     9    1001    14
c    10    1010    12
c    11    1011    13
c    12    1100     8
c    13    1101     9
c    14    1110    11
c    15    1111    10
c    16   10000    31
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    20 January 2007
c
c  Parameters:
c
c    Input, integer GRAY, the Gray code to be ranked.
c
c    Output, integer RANK, the rank of GRAY, and the integer whose Gray
c    code is GRAY.
c
      implicit none

      integer gray
      integer i
      integer i4_bset
      logical i4_btest
      integer j
      integer nbits
      parameter ( nbits = 32 )
      integer rank

      j = 0

      if ( i4_btest ( gray, nbits-1 ) ) then
        j = i4_bset ( j, nbits-1 )
      end if

      do i = nbits-2, 0, -1

        if ( i4_btest ( j, i+1 ) .and. 
     &    ( .not. i4_btest ( gray, i ) ) ) then
          j = i4_bset ( j, i )
        end if

        if ( ( .not. i4_btest ( j, i+1 ) ) .and. 
     &    i4_btest ( gray, i ) ) then
          j = i4_bset ( j, i )
        end if

      end do

      rank = j

      return
      end
      subroutine gray_rank2 ( gray, rank )

c*********************************************************************72
c
cc GRAY_RANK2 ranks a Gray code.
c
c  Discussion:
c
c    In contrast to GRAY_RANK, this routine is entirely arithmetical,
c    and does not require access to bit testing and setting routines.
c
c    Given the number GRAY, its ranking is the order in which it would be
c    visited in the Gray code ordering.  The Gray code ordering begins
c
c    Rank  Gray  Gray
c          (Dec) (Bin)
c
c       0     0  0000
c       1     1  0001
c       2     3  0011
c       3     2  0010
c       4     6  0110
c       5     7  0111
c       6     5  0101
c       7     4  0100
c       8    12  0110
c       etc
c
c   This routine is given a Gray code, and has to return the rank.
c
c  Example:
c
c    Gray  Gray  Rank
c    (Dec) (Bin) 
c
c     0       0     0
c     1       1     1
c     2      10     3
c     3      11     2
c     4     100     7
c     5     101     6
c     6     110     4
c     7     111     5
c     8    1000    15
c     9    1001    14
c    10    1010    12
c    11    1011    13
c    12    1100     8
c    13    1101     9
c    14    1110    11
c    15    1111    10
c    16   10000    31
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    20 January 2007
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer GRAY, the Gray code to be ranked.
c
c    Output, integer RANK, the rank of GRAY, and the integer whose Gray
c    code is GRAY.
c
      implicit none

      integer gray
      integer gray_copy
      integer k
      logical last
      logical next
      integer rank
      integer two_k

      gray_copy = gray

      if ( gray_copy .lt. 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'GRAY_RANK2 - Fatal error!'
        write ( *, '(a)' ) '  Input value of GRAY < 0.'
        stop
      end if

      if ( gray_copy .eq. 0 ) then
        rank = 0
        return
      end if
c
c  Find TWO_K, the largest power of 2 less than or equal to GRAY.
c
      k = 0
      two_k = 1

10    continue

      if ( 2 * two_k .le. gray_copy ) then
        two_k = two_k * 2
        k = k + 1
        go to 10
      end if

      rank = two_k
      last = .true.
      gray_copy = gray_copy - two_k

20    continue

      if ( 0 .lt. k ) then

        two_k = two_k / 2
        k = k - 1

        next = ( two_k .le. gray_copy .and. 
     &           gray_copy .lt. two_k * 2 )

        if ( next ) then
          gray_copy = gray_copy - two_k
        end if

        if ( next .neqv. last ) then
          rank = rank + two_k
          last = .true.
        else
          last = .false.
        end if

        go to 20

      end if

      return
      end
      subroutine gray_unrank ( rank, gray )

c*********************************************************************72
c
cc GRAY_UNRANK unranks a Gray code.
c
c  Discussion:
c
c    The binary values of the Gray codes of successive integers differ in
c    just one bit.
c
c    The sequence of Gray codes for 0 to (2**N)-1 can be interpreted as a
c    Hamiltonian cycle on a graph of the cube in N dimensions.
c
c  Example:
c
c    Rank  Gray  Gray
c          (Dec) (Bin)
c
c     0     0       0
c     1     1       1
c     2     3      11
c     3     2      10
c     4     6     110
c     5     7     111
c     6     5     101
c     7     4     100
c     8    12    1100
c     9    14    1001
c    10    12    1010
c    11    13    1011
c    12     8    1100
c    13     9    1101
c    14    11    1110
c    15    10    1111
c    16    31   10000
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    20 January 2007
c
c  Parameters:
c
c    Input, integer RANK, the integer whose Gray code is desired.
c
c    Output, integer GRAY, the Gray code of the given rank.
c
      implicit none

      integer gray
      integer i
      integer i4_bset
      logical i4_btest
      integer j
      integer nbits
      parameter ( nbits = 32 )
      integer rank

      j = 0

      if ( i4_btest ( rank, nbits-1 ) ) then
        j = i4_bset ( j, nbits-1 )
      end if

      do i = nbits-2, 0, -1

        if ( i4_btest ( rank, i+1 ) .neqv. i4_btest ( rank, i ) ) then
          j = i4_bset ( j, i )
        end if

      end do

      gray = j

      return
      end
      subroutine gray_unrank2 ( rank, gray )

c*********************************************************************72
c
cc GRAY_UNRANK2 unranks a Gray code.
c
c  Discussion:
c
c    In contrast to GRAY_UNRANK, this routine is entirely arithmetical,
c    and does not require access to bit testing and setting routines.
c
c    The binary values of the Gray codes of successive integers differ in
c    just one bit.
c
c    The sequence of Gray codes for 0 to (2**N)-1 can be interpreted as a
c    Hamiltonian cycle on a graph of the cube in N dimensions.
c
c  Example:
c
c    Rank  Gray  Gray
c          (Dec) (Bin)
c
c     0     0       0
c     1     1       1
c     2     3      11
c     3     2      10
c     4     6     110
c     5     7     111
c     6     5     101
c     7     4     100
c     8    12    1100
c     9    14    1001
c    10    12    1010
c    11    13    1011
c    12     8    1100
c    13     9    1101
c    14    11    1110
c    15    10    1111
c    16    31   10000
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    20 January 2007
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer RANK, the integer whose Gray code is desired.
c
c    Output, integer GRAY, the Gray code of the given rank.
c
      implicit none

      integer gray
      integer k
      logical last
      logical next
      integer rank
      integer rank_copy
      integer two_k

      if ( rank .le. 0 ) then
        gray = 0
        return
      end if

      rank_copy = rank
      k = 0
      two_k = 1

10    continue

      if ( 2 * two_k .le. rank_copy ) then
        two_k = two_k * 2
        k = k + 1
      end if

      gray = two_k
      rank_copy = rank_copy - two_k
      next = .true.

20    continue

      if ( 0 .lt. k ) then

        two_k = two_k / 2
        k = k - 1

        last = next
        next = ( two_k .le. rank_copy .and. rank_copy .le. two_k * 2 )

        if ( next .neqv. last ) then
          gray = gray + two_k
        end if

        if ( next ) then
          rank_copy = rank_copy - two_k
        end if

      end if

      return
      end
      function i4_bclr ( i4, pos )

c*********************************************************************72
c
cc I4_BCLR returns a copy of an I4 in which the POS-th bit is set to 0.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    23 January 2007
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Military Standard 1753,
c    FORTRAN, DoD Supplement To American National Standard X3.9-1978,
c    9 November 1978.
c
c  Parameters:
c
c    Input, integer I4, the integer.
c
c    Input, integer POS, the bit position, between 0 and 31.
c
c    Output, integer I4_BCLR, a copy of I4, but with the POS-th bit
c    set to 0.
c
      implicit none

      integer i4
      integer i4_bclr
      integer i4_huge
      integer j
      integer k
      integer pos
      integer sub
      integer value

      value = i4

      if ( pos .lt. 0 ) then

      else if ( pos .lt. 31 ) then

        sub = 1

        if ( 0 .le. i4 ) then
          j = i4
        else
          j = ( i4_huge ( ) + i4 ) + 1
        end if

        do k = 1, pos
          j = j / 2
          sub = sub * 2
        end do

        if ( mod ( j, 2 ) .eq. 1 ) then
          value = i4 - sub
        end if

      else if ( pos .eq. 31 ) then

        if ( i4 .lt. 0 ) then
          value = ( i4_huge ( ) + i4 ) + 1
        end if

      else if ( 31 .lt. pos ) then

        value = i4

      end if

      i4_bclr = value

      return
      end
      function i4_bset ( i4, pos )

c*********************************************************************72
c
cc I4_BSET returns a copy of an I4 in which the POS-th bit is set to 1.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    23 January 2007
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Military Standard 1753,
c    FORTRAN, DoD Supplement To American National Standard X3.9-1978,
c    9 November 1978.
c
c  Parameters:
c
c    Input, integer I4, the integer to be tested.
c
c    Input, integer POS, the bit position, between 0 and 31.
c
c    Output, integer I4_BSET, a copy of I4, but with the POS-th bit
c    set to 1.
c
      implicit none

      integer add
      integer i4
      integer i4_bset
      integer i4_huge
      integer j
      integer k
      integer pos
      integer value

      value = i4

      if ( pos .lt. 0 ) then

      else if ( pos .lt. 31 ) then

        add = 1

        if ( 0 .le. i4 ) then
          j = i4
        else
          j = ( i4_huge ( ) + i4 ) + 1
        end if

        do k = 1, pos
          j = j / 2
          add = add * 2
        end do

        if ( mod ( j, 2 ) .eq. 0 ) then
          value = i4 + add
        end if

      else if ( pos .eq. 31 ) then

        if ( 0 .lt. i4 ) then
          value = - ( i4_huge ( ) - i4 ) - 1
        end if

      else if ( 31 .lt. pos ) then

        value = i4

      end if

      i4_bset = value

      return
      end
      function i4_btest ( i4, pos )

c*********************************************************************72
c
cc I4_BTEST returns TRUE if the POS-th bit of an I4 is 1.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    23 January 2007
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Military Standard 1753,
c    FORTRAN, DoD Supplement To American National Standard X3.9-1978,
c    9 November 1978.
c
c  Parameters:
c
c    Input, integer I4, the integer to be tested.
c
c    Input, integer POS, the bit position, between 0 and 31.
c
c    Output, logical I4_BTEST, is TRUE if the POS-th bit is 1.
c
      implicit none

      integer i4
      logical i4_btest
      integer i4_huge
      integer j
      integer k
      integer pos

      if ( pos .lt. 0 ) then

        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'I4_BTEST - Fatal error!'
        write ( *, '(a)' ) '  POS < 0.'
        stop

      else if ( pos .lt. 31 ) then

        if ( 0 .le. i4 ) then
          j = i4
        else
          j = ( i4_huge ( ) + i4 ) + 1
        end if

        do k = 1, pos
          j = j / 2
        end do

        if ( mod ( j, 2 ) .eq. 0 ) then
          i4_btest = .false.
        else
          i4_btest = .true.
        end if

      else if ( pos .eq. 31 ) then

        if ( i4 .lt. 0 ) then
          i4_btest = .true.
        else
          i4_btest = .false.
        end if

      else if ( 31 .lt. pos ) then

        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'I4_BTEST - Fatal error!'
        write ( *, '(a)' ) '  31 < POS.'
        stop

      end if

      return
      end
      subroutine i4_factor ( n, factor_max, factor_num, factor, power, 
     &  nleft )

c*********************************************************************72
c
cc I4_FACTOR factors an I4 into prime factors.
c
c  Discussion:
c
c    The formula used is:
c
c      N = NLEFT * product ( 1 <= I <= FACTOR_NUM ) FACTOR(I)**POWER(I).
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    23 January 2007
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the integer to be factored.  N may be positive,
c    negative, or 0.
c
c    Input, integer FACTOR_MAX, the maximum number of prime factors for
c    which storage has been allocated.
c
c    Output, integer FACTOR_NUM, the number of prime factors of N discovered
c    by the routine.
c
c    Output, integer FACTOR(FACTOR_MAX), the prime factors of N.
c
c    Output, integer POWER(FACTOR_MAX).  POWER(I) is the power of
c    the FACTOR(I) in the representation of N.
c
c    Output, integer NLEFT, the factor of N that the routine could not
c    divide out.  If NLEFT is 1, then N has been completely factored.
c    Otherwise, NLEFT represents factors of N involving large primes.
c
      implicit none

      integer factor_max

      integer factor(factor_max)
      integer factor_num
      integer i
      integer n
      integer nleft
      integer p
      integer power(factor_max)
      integer prime
      integer prime_max

      factor_num = 0

      do i = 1, factor_max
        factor(i) = 0
      end do

      do i = 1, factor_max
        power(i) = 0
      end do

      nleft = n

      if ( n .eq. 0 ) then
        return
      end if

      if ( abs ( n ) .eq. 1 ) then
        factor_num = 1
        factor(1) = 1
        power(1) = 1
        return
      end if
c
c  Find out how many primes we stored.
c
      prime_max = prime ( -1 )
c
c  Try dividing the remainder by each prime.
c
      do i = 1, prime_max

        p = prime ( i )

        if ( mod ( abs ( nleft ), p ) .eq. 0 ) then
    
          if ( factor_num .lt. factor_max ) then

            factor_num = factor_num + 1
            factor(factor_num) = p
            power(factor_num) = 0

10          continue

              power(factor_num) = power(factor_num) + 1
              nleft = nleft / p

              if ( mod ( abs ( nleft ), p ) .ne. 0 ) then
                go to 20
              end if

            go to 10

20          continue

            if ( abs ( nleft ) .eq. 1 ) then
              go to 30
            end if

          end if

        end if

      end do

30    continue

      return
      end
      function i4_factorial ( n )

c*********************************************************************72
c
cc I4_FACTORIAL computes the factorial of N.
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
c    26 June 2008
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the argument of the factorial function.
c    If N is less than 1, the function value is returned as 1.
c    0 <= N <= 13 is required.
c
c    Output, integer I4_FACTORIAL, the factorial of N.
c
      implicit none

      integer i
      integer i4_factorial
      integer n

      i4_factorial = 1

      if ( 13 .lt. n ) then
        i4_factorial = - 1
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'I4_FACTORIAL - Fatal error!'
        write ( *, '(a)' ) 
     &  '  I4_FACTORIAL(N) cannot be computed as an integer'
        write ( *, '(a)' ) '  for 13 < N.'
        write ( *, '(a,i8)' ) '  Input value N = ', n
        stop
      end if

      do i = 1, n
        i4_factorial = i4_factorial * i
      end do

      return
      end
      function i4_gcd ( i, j )

c*********************************************************************72
c
cc I4_GCD finds the greatest common divisor of I and J.
c
c  Discussion:
c
c    Only the absolute values of I and J are
c    considered, so that the result is always nonnegative.
c
c    If I or J is 0, I4_GCD is returned as max ( 1, abs ( I ), abs ( J ) ).
c
c    If I and J have no common factor, I4_GCD is returned as 1.
c
c    Otherwise, using the Euclidean algorithm, I4_GCD is the
c    largest common factor of I and J.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    03 March 1999
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer I, J, two numbers whose greatest common divisor
c    is desired.
c
c    Output, integer I4_GCD, the greatest common divisor of I and J.
c
      implicit none

      integer i
      integer i4_gcd
      integer ip
      integer iq
      integer ir
      integer j

      i4_gcd = 1
c
c  Return immediately if either I or J is zero.
c
      if ( i .eq. 0 ) then
        i4_gcd = max ( 1, abs ( j ) )
        return
      else if ( j .eq. 0 ) then
        i4_gcd = max ( 1, abs ( i ) )
        return
      end if
c
c  Set IP to the larger of I and J, IQ to the smaller.
c  This way, we can alter IP and IQ as we go.
c
      ip = max ( abs ( i ), abs ( j ) )
      iq = min ( abs ( i ), abs ( j ) )
c
c  Carry out the Euclidean algorithm.
c
10    continue

        ir = mod ( ip, iq )

        if ( ir .eq. 0 ) then
          go to 20
        end if

        ip = iq
        iq = ir

      go to 10

20    continue

      i4_gcd = iq

      return
      end
      function i4_huge ( )

c*********************************************************************72
c
cc I4_HUGE returns a "huge" I4.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    13 November 2006
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, integer I4_HUGE, a huge number.
c
      implicit none

      integer i4_huge

      i4_huge = 2147483647

      return
      end
      function i4_log_10 ( i )

c*********************************************************************72
c
cc I4_LOG_10 returns the integer part of the logarithm base 10 of an I4.
c
c  Example:
c
c        I  I4_LOG_10
c    -----  --------
c        0    0
c        1    0
c        2    0
c        9    0
c       10    1
c       11    1
c       99    1
c      100    2
c      101    2
c      999    2
c     1000    3
c     1001    3
c     9999    3
c    10000    4
c
c  Discussion:
c
c    I4_LOG_10 ( I ) + 1 is the number of decimal digits in I.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    04 January 2007
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer I, the number whose logarithm base 10 is desired.
c
c    Output, integer I4_LOG_10, the integer part of the logarithm base 10 of
c    the absolute value of X.
c
      implicit none

      integer i
      integer i_abs
      integer i4_log_10
      integer ten_pow

      if ( i .eq. 0 ) then

        i4_log_10 = 0

      else

        i4_log_10 = 0
        ten_pow = 10

        i_abs = abs ( i )

10      continue

        if ( ten_pow .le. i_abs ) then
          i4_log_10 = i4_log_10 + 1
          ten_pow = ten_pow * 10
          go to 10
        end if

      end if

      return
      end
      function i4_modp ( i, j )

c*********************************************************************72
c
cc I4_MODP returns the nonnegative remainder of I4 division.
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
      subroutine i4_moebius ( n, mu )

c*********************************************************************72
c
cc I4_MOEBIUS returns the value of MU(N), the Moebius function of N.
c
c  Discussion:
c
c    MU(N) is defined as follows:
c
c      MU(N) = 1 if N = 1;
c              0 if N is divisible by the square of a prime;
c              (-1)**K, if N is the product of K distinct primes.
c
c    As special cases, MU(N) is -1 if N is a prime, and MU(N) is 0
c    if N is a square, cube, etc.
c
c    The Moebius function MU(D) is related to Euler's totient 
c    function PHI(N):
c
c      PHI(N) = sum ( D divides N ) MU(D) * ( N / D ).
c
c  First values:
c
c     N  MU(N)
c
c     1    1
c     2   -1
c     3   -1
c     4    0
c     5   -1
c     6    1
c     7   -1
c     8    0
c     9    0
c    10    1
c    11   -1
c    12    0
c    13   -1
c    14    1
c    15    1
c    16    0
c    17   -1
c    18    0
c    19   -1
c    20    0
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    20 January 2007
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the value to be analyzed.
c
c    Output, integer MU, the value of MU(N).
c    If N is less than or equal to 0, MU will be returned as -2.
c    If there was not enough internal space for factoring, MU
c    is returned as -3.
c
      implicit none

      integer maxfactor
      parameter ( maxfactor = 20 )

      integer exponent(maxfactor)
      integer factor(maxfactor)
      integer i
      integer mu
      integer n
      integer nfactor
      integer nleft

      if ( n .le. 0 ) then
        mu = -2
        return
      end if

      if ( n .eq. 1 ) then
        mu = 1
        return
      end if
c
c  Factor N.
c
      call i4_factor ( n, maxfactor, nfactor, factor, exponent, nleft )

      if ( nleft .ne. 1 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'I4_MOEBIUS - Fatal error!'
        write ( *, '(a)' ) '  Not enough factorization space.'
        mu = -3
        stop
      end if

      mu = 1

      do i = 1, nfactor

        mu = -mu

        if ( 1 .lt. exponent(i) ) then
          mu = 0
          return
        end if

      end do

      return
      end
      subroutine i4_partition_conj ( n, iarray1, mult1, npart1, 
     &  iarray2, mult2, npart2 )

c*********************************************************************72
c
cc I4_PARTITION_CONJ computes the conjugate of a partition.
c
c  Discussion:
c
c    A partition of an integer N is a representation of the integer
c    as the sum of nonzero positive integers.  The order of the summands
c    does not matter.  Thus, the number 5 has the following partitions
c    and no more:
c
c    5 = 5
c      = 4 + 1 
c      = 3 + 2 
c      = 3 + 1 + 1 
c      = 2 + 2 + 1 
c      = 2 + 1 + 1 + 1 
c      = 1 + 1 + 1 + 1 + 1
c
c    so the number of partitions of 5 is 7.
c
c    The conjugate of a partition P1 of N is another partition P2 of N 
c    obtained in the following way:
c
c      The first element of P2 is the number of parts of P1 greater than
c      or equal to 1.
c
c      The K-th element of P2 is the number of parts of P1 greater than
c      or equal to K.
c
c    Clearly, P2 will have no more than N elements; it may be surprising
c    to find that P2 is guaranteed to be a partition of N.  However, if
c    we symbolize the initial partition P1 by rows of X's, then we can
c    see that P2 is simply produced by grouping by columns:
c
c        6 3 2 2 1
c      5 X X X X X
c      4 X X X X
c      2 X X
c      1 X
c      1 X
c      1 X
c
c  Example:
c
c    14 = 5 + 4 + 2 + 1 + 1 + 1
c
c    The conjugate partition is:
c
c    14 = 6 + 3 + 2 + 2 + 1
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    28 January 2007
c
c  Author:
c
c    John Burkardt
c
c  Parameters
c
c    Input, integer N, the integer to be partitioned.
c
c    Input, integer IARRAY1(NPART1).  IARRAY1 contains the parts of
c    the partition.  The value of N is represented by
c      sum ( 1 <= I <= NPART1 ) MULT1(I) * IARRAY1(I).
c    Input, integer MULT1(NPART1).  MULT1 counts the multiplicity of
c    the parts of the partition.  MULT1(I) is the multiplicity
c    of the part IARRAY1(I), for I = 1 to NPART1.
c
c    Input, integer N, the integer to be partitioned.
c
c    Input, integer NPART1, the number of "parts" in the partition.
c
c    Output, integer IARRAY2(N).  IARRAY contains the parts of
c    the conjugate partition in entries 1 through NPART2.
c
c    Output, integer MULT2(N).  MULT2 counts the multiplicity of
c    the parts of the conjugate partition in entries 1 through NPART2.
c
c    Output, integer NPART2, the number of "parts" in the conjugate partition.
c
      implicit none

      integer n
      integer npart1

      integer i
      integer iarray1(npart1)
      integer iarray2(n)
      integer itemp
      integer itest
      integer mult1(npart1)
      integer mult2(n)
      integer npart2

      do i = 1, n
        iarray2(i) = 0
      end do

      do i = 1, n
        mult2(i) = 0
      end do

      npart2 = 0

      itest = 0

10    continue

        itest = itest + 1

        itemp = 0

        do i = 1, npart1
          if ( itest .le. iarray1(i) ) then
            itemp = itemp + mult1(i)
          end if
        end do

        if ( itemp .le. 0 ) then
          go to 20
        end if

        if ( 0 .lt. npart2 ) then

          if ( itemp .eq. iarray2(npart2) ) then
            mult2(npart2) = mult2(npart2) + 1
          else
            npart2 = npart2 + 1
            iarray2(npart2) = itemp
            mult2(npart2) = 1
          end if

        else

          npart2 = npart2 + 1
          iarray2(npart2) = itemp
          mult2(npart2) = 1

        end if

      go to 10

20    continue

      return
      end
      subroutine i4_partition_count ( n, p )

c*********************************************************************72
c
cc I4_PARTITION_COUNT computes the number of partitions of an I4.
c
c  Discussion:
c
c    A partition of an integer N is a representation of the integer
c    as the sum of nonzero positive integers.  The order of the summands
c    does not matter.  Thus, the number 5 has the following partitions
c    and no more:
c
c    5 = 5
c      = 4 + 1 
c      = 3 + 2 
c      = 3 + 1 + 1 
c      = 2 + 2 + 1 
c      = 2 + 1 + 1 + 1 
c      = 1 + 1 + 1 + 1 + 1
c
c    so the number of partitions of 5 is 7.
c
c    Partition numbers are difficult to compute.  This routine uses
c    Euler's method, which observes that:
c
c      P(0) = 1
c      P(N) =   P(N-1)  + P(N-2)
c             - P(N-5)  - P(N-7)
c             + P(N-12) + P(N-15)
c             - ...
c
c      where the numbers 1, 2, 5, 7, ... to be subtracted from N in the
c      indices are the successive pentagonal numbers, (with both positive 
c      and negative indices) with the summation stopping when a negative 
c      index is reached.
c
c  First values:
c
c    N   P
c
c    0   1
c    1   1
c    2   2
c    3   3
c    4   5
c    5   7
c    6  11
c    7  15
c    8  22
c    9  30
c   10  42
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    20 January 2007
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    John Conway, Richard Guy,
c    The Book of Numbers,
c    Springer Verlag, 1996, page 95.
c
c  Parameters:
c
c    Input, integer N, the index of the highest partition number desired.
c
c    Output, integer P(0:N), the partition numbers.
c
      implicit none

      integer n

      integer i
      integer j
      integer p(0:n)
      integer pj
      integer sgn

      p(0) = 1

      do i = 1, n

        p(i) = 0

        j = 0
        sgn = 1

10      continue

          j = j + 1
          call pent_enum ( j, pj )

          if ( i .lt. pj ) then
            go to 20
          end if

          p(i) = p(i) + sgn * p(i-pj)
          sgn = -sgn

        go to 10

20      continue

        j = 0
        sgn = 1

30      continue

          j = j - 1
          call pent_enum ( j, pj )

          if ( i .lt. pj ) then
            go to 40
          end if

          p(i) = p(i) + sgn * p(i-pj)
          sgn = -sgn

        go to 30

40      continue

      end do

      return
      end
      subroutine i4_partition_count2 ( n, p )

c*********************************************************************72
c
cc I4_PARTITION_COUNT2 computes the number of partitions of an I4.
c
c  Discussion:
c
c    A partition of an integer N is a representation of the integer
c    as the sum of nonzero positive integers.  The order of the summands
c    does not matter.  Thus, the number 5 has the following partitions
c    and no more:
c
c    5 = 5
c      = 4 + 1 
c      = 3 + 2 
c      = 3 + 1 + 1 
c      = 2 + 2 + 1 
c      = 2 + 1 + 1 + 1 
c      = 1 + 1 + 1 + 1 + 1
c
c    so the number of partitions of 5 is 7.
c
c  First values:
c
c    N   P
c
c    0   1
c    1   1
c    2   2
c    3   3
c    4   5
c    5   7
c    6  11
c    7  15
c    8  22
c    9  30
c   10  42
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    20 January 2007
c
c  Author:
c
c    Original FORTRAN77 version by Albert Nijenhuis, Herbert Wilf.
c    This FORTRAN77 version by John Burkardt.
c
c  Reference:
c
c    Albert Nijenhuis, Herbert Wilf,
c    Combinatorial Algorithms for Computers and Calculators,
c    Second Edition,
c    Academic Press, 1978,
c    ISBN: 0-12-519260-6,
c    LC: QA164.N54.
c
c  Parameters:
c
c    Input, integer N, the largest integer to be considered.
c
c    Output, integer P(0:N), the partition numbers.
c
      implicit none

      integer n

      integer i
      integer j
      integer p(0:n)
      integer s
      integer t
      integer total

      if ( n .lt. 0 ) then
        return
      end if

      p(0) = 1

      if ( n .lt. 1 ) then
        return
      end if

      p(1) = 1

      do i = 2, n

        total = 0

        do t = 1, i

          s = 0
          j = i

10        continue

            j = j - t

            if ( j .lt. 0 ) then
              go to 20
            end if

            s = s + p(j)

          go to 10

20        continue

          total = total + s * t

        end do

        p(i) = total / i

      end do

      return
      end
      subroutine i4_partition_count_values ( n_data, n, c )

c*********************************************************************72
c
cc I4_PARTITION_COUNT_VALUES returns some values of the integer partition count.
c
c  Discussion:
c
c    A partition of an integer N is a representation of the integer
c    as the sum of nonzero positive integers.  The order of the summands
c    does not matter.  The number of partitions of N is symbolized
c    by P(N).  Thus, the number 5 has P(N) = 7, because it has the
c    following partitions:
c
c    5 = 5
c      = 4 + 1
c      = 3 + 2
c      = 3 + 1 + 1
c      = 2 + 2 + 1
c      = 2 + 1 + 1 + 1
c      = 1 + 1 + 1 + 1 + 1
c
c    In Mathematica, the function can be evaluated by
c
c      PartitionsP[n]
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    20 January 2007
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Milton Abramowitz, Irene Stegun,
c    Handbook of Mathematical Functions,
c    US Department of Commerce, 1964,
c    ISBN: 0-486-61272-4,
c    LC: QA47.A34.
c
c    Stephen Wolfram,
c    The Mathematica Book,
c    Fourth Edition,
c    Wolfram Media / Cambridge University Press, 1999.
c
c  Parameters:
c
c    Input/output, integer N_DATA.  The user sets N_DATA to 0 before the
c    first call.  On each call, the routine increments N_DATA by 1, and
c    returns the corresponding data; when there is no more data, the
c    output value of N_DATA will be 0 again.
c
c    Output, integer N, the integer.
c
c    Output, integer C, the number of partitions of the integer.
c
      implicit none

      integer n_max
      parameter ( n_max = 21 )

      integer c
      integer c_vec(n_max)
      integer n
      integer n_data
      integer n_vec(n_max)

      save c_vec
      save n_vec

      data c_vec / 
     &    1, 
     &    1,   2,   3,   5,   7,  11,  15,  22,  30,  42, 
     &   56,  77, 101, 135, 176, 231, 297, 385, 490, 627 /
      data n_vec /
     &   0,  
     &   1,  2,  3,  4,  5,  6,  7,  8,  9, 10, 
     &  11, 12, 13, 14, 15, 16, 17, 18, 19, 20 /


      if ( n_data .lt. 0 ) then
        n_data = 0
      end if

      n_data = n_data + 1

      if ( n_max .lt. n_data ) then
        n_data = 0
        n = 0
        c = 0
      else
        n = n_vec(n_data)
        c = c_vec(n_data)
      end if

      return
      end
      subroutine i4_partition_next ( n, npart, a, mult, done )

c*********************************************************************72
c
cc I4_PARTITION_NEXT generates the partitions of an I4, one at a time.
c
c  Discussion:
c
c    A partition of an integer N is a representation of the integer
c    as the sum of nonzero positive integers.  The order of the summands
c    does not matter.  Thus, the number 5 has the following partitions
c    and no more:
c
c    5 = 5
c      = 4 + 1 
c      = 3 + 2 
c      = 3 + 1 + 1 
c      = 2 + 2 + 1 
c      = 2 + 1 + 1 + 1 
c      = 1 + 1 + 1 + 1 + 1
c
c    so the number of partitions of 5 is 7.
c
c    The number of partitions of N is:
c
c      1     1
c      2     2
c      3     3
c      4     5
c      5     7
c      6    11
c      7    15
c      8    22
c      9    30
c     10    42
c     11    56
c     12    77
c     13   101
c     14   135
c     15   176
c     16   231
c     17   297
c     18   385
c     19   490
c     20   627
c     21   792
c     22  1002
c     23  1255
c     24  1575
c     25  1958
c     26  2436
c     27  3010
c     28  3718
c     29  4565
c     30  5604
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    21 January 2007
c
c  Author:
c
c    Original FORTRAN77 version by Albert Nijenhuis, Herbert Wilf.
c    This FORTRAN77 version by John Burkardt.
c
c  Reference:
c
c    Albert Nijenhuis, Herbert Wilf,
c    Combinatorial Algorithms for Computers and Calculators,
c    Second Edition,
c    Academic Press, 1978,
c    ISBN: 0-12-519260-6,
c    LC: QA164.N54.
c
c  Parameters:
c
c    Input, integer N, the integer to be partitioned.
c
c    Input/output, integer NPART, the number of "parts" in the partition.
c
c    Input/output, integer A(N).  A contains the parts of
c    the partition.  The value of N is represented by
c      N = sum ( 1 <= I <= NPART ) MULT(I) * A(I).
c
c    Input/output, integer MULT(N).  MULT counts the multiplicity of
c    the parts of the partition.  MULT(I) is the multiplicity
c    of the part A(I), for I = 1 to NPART.
c
c    Input/output, logical DONE.
c    On first call, the user should set DONE to TRUE to signal
c    that the program should initialize data.
c    On each return, the programs sets DONE to FALSE if it
c    has another partition to return.  If the program returns
c    with DONE TRUE, then there are no more partitions.
c
      implicit none

      integer n

      integer a(n)
      logical done
      integer i
      integer is
      integer iu
      integer iv
      integer iw
      integer k
      integer k1
      integer mult(n)
      integer npart

      if ( n .le. 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'I4_PARTITION_NEXT - Fatal error!'
        write ( *, '(a)' ) '  N must be positive.'
        write ( *, '(a,i8)' ) '  The input value of N was ', n
        stop
      end if

      if ( done ) then

        a(1) = n
        do i = 2, n
          a(i) = 0
        end do

        mult(1) = 1
        do i = 2, n
          mult(i) = 0
        end do

        npart = 1
        done = .false.

      else

        if ( 1 .lt. a(npart) .or. 1 .lt. npart ) then

          done = .false.

          if ( a(npart) .eq. 1 ) then
            is = a(npart-1) + mult(npart)
            k = npart - 1
          else
            is = a(npart)
            k = npart
          end if

          iw = a(k) - 1
          iu = is / iw
          iv = mod ( is, iw )
          mult(k) = mult(k) - 1

          if ( mult(k) .eq. 0 ) then
            k1 = k
          else
            k1 = k + 1
          end if

          mult(k1) = iu
          a(k1) = iw

          if ( iv .eq. 0 ) then
            npart = k1
          else
            mult(k1+1) = 1
            a(k1+1) = iv
            npart = k1 + 1
          end if

        else
          done = .true.
        end if

      end if

      return
      end
      subroutine i4_partition_next2 ( n, npart, a, mult, more )

c*********************************************************************72
c
cc I4_PARTITION_NEXT2 computes the partitions of the integer N one at a time.
c
c  Discussion:
c
c    A partition of an integer N is a representation of the integer
c    as the sum of nonzero positive integers.  The order of the summands
c    does not matter.  Thus, the number 5 has the following partitions
c    and no more:
c
c    5 = 5
c      = 4 + 1 
c      = 3 + 2 
c      = 3 + 1 + 1 
c      = 2 + 2 + 1 
c      = 2 + 1 + 1 + 1 
c      = 1 + 1 + 1 + 1 + 1
c
c    so the number of partitions of 5 is 7.
c
c    Unlike compositions, order is not important in a partition.  Thus
c    the sequences 3+2+1 and 1+2+3 represent distinct compositions, but
c    not distinct partitions.  Also 0 is never returned as one of the
c    elements of the partition.
c
c  Example:
c
c    Sample partitions of 6 include:
c
c      6 = 4+1+1 = 3+2+1 = 2+2+2
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    21 January 2007
c
c  Parameters:
c
c    Input, integer N, the integer whose partitions are desired.
c
c    Input/output, integer NPART, the number of distinct, nonzero parts in the
c    output partition.
c
c    Input/output, integer A(N).  A(I) is the I-th distinct part
c    of the partition, for I = 1, NPART.  Note that if a certain number
c    shows up several times in the partition, it is listed only
c    once in A, and its multiplicity is counted in MULT.
c
c    Input/output, integer MULT(N).  MULT(I) is the multiplicity of A(I)
c    in the partition, for I = 1, NPART; that is, the number of repeated
c    times that A(I) is used in the partition.
c
c    Input/output, logical MORE.  Set MORE = FALSE on first call.  It
c    will be reset TRUE on return with the first partition.
c    Keep calling for more partitions until MORE
c    is returned FALSE.
c
      implicit none

      integer n

      integer a(n)
      integer iff
      integer is
      integer isum
      logical more
      integer mult(n)
      integer nlast
      integer npart

      save nlast

      data nlast / 0 /
c
c  On the first call, set NLAST to 0.
c
      if ( .not. more ) then
        nlast = 0
      end if

      if ( n .ne. nlast .or. ( .not. more ) ) then
        nlast = n
        npart = 1
        a(npart) = n
        mult(npart) = 1
        more = ( mult(npart) .ne. n )
        return
      end if

      isum = 1

      if ( a(npart) .le. 1 ) then
        isum = mult(npart) + 1
        npart = npart - 1
      end if

      iff = a(npart) - 1

      if ( mult(npart) .ne. 1 ) then
        mult(npart) = mult(npart) - 1
        npart = npart + 1
      end if

      a(npart) = iff
      mult(npart) = 1 + ( isum / iff )
      is = mod ( isum, iff )

      if ( 0 .lt. is ) then
        npart = npart + 1
        a(npart) = is
        mult(npart) = 1
      end if
c
c  There are more partitions, as long as we haven't just computed
c  the last one, which is N copies of 1.
c
      more = ( mult(npart) .ne. n )

      return
      end
      subroutine i4_partition_print ( n, npart, a, mult )

c*********************************************************************72
c
cc I4_PARTITION_PRINT prints a partition of an I4.
c
c  Discussion:
c
c    A partition of an integer N is a representation of the integer as
c    the sum of nonzero integers:
c
c      N = A1 + A2 + A3 + ...
c
c    It is standard practice to gather together all the values that 
c    are equal, and replace them in the sum by a single term, multiplied
c    by its "multiplicity":
c
c      N = M1 * A1 + M2 * A2 + ... + M(NPART) * A(NPART)
c    
c    In this representation, every A is a unique positive number, and 
c    no M is zero (or negative).
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    21 January 2007
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the integer to be partitioned.
c
c    Input, integer NPART, the number of "parts" in the partition.
c
c    Input, integer A(NPART), the parts of the partition.  
c
c    Input, integer MULT(NPART), the multiplicies of the parts.
c
      implicit none

      integer npart

      integer a(npart)
      integer i
      integer j
      integer mult(npart)
      integer n
      character * ( 80 ) string

      string(1:2) = '  '
      call i4_to_s_left ( n, string(3:) )
      j = len_trim ( string )
      string (j+1:j+3) = ' = '
      j = j + 3

      do i = 1, npart

        if ( 1 .lt. i ) then
          string(j+1:j+3) = ' + '
          j = j + 3
        end if

        call i4_to_s_left ( mult(i), string(j+1:) )
        j = len_trim ( string )

        string(j+1:j+3) = ' * '
        j = j + 3

        call i4_to_s_left ( a(i), string(j+1:) )
        j = len_trim ( string )

      end do

      j = len_trim ( string )
      write ( *, '(a)' ) string(1:j)

      return
      end
      subroutine i4_partition_random ( n, table, seed, a, mult, npart )

c*********************************************************************72
c
cc I4_PARTITION_RANDOM selects a random partition of an integer.
c
c  Discussion:
c
c    A partition of an integer N is a representation of the integer
c    as the sum of nonzero positive integers.  The order of the summands
c    does not matter.  Thus, the number 5 has the following partitions
c    and no more:
c
c    5 = 5
c      = 4 + 1 
c      = 3 + 2 
c      = 3 + 1 + 1 
c      = 2 + 2 + 1 
c      = 2 + 1 + 1 + 1 
c      = 1 + 1 + 1 + 1 + 1
c
c    so the number of partitions of 5 is 7.
c
c    Note that some elements of the partition may be 0.  The partition is
c    returned as (MULT(I),I), with NPART nonzero entries in MULT, and
c
c      N = sum ( 1 <= I <= N ) MULT(I) * I.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    21 January 2007
c
c  Author:
c
c    Original FORTRAN77 version by Albert Nijenhuis, Herbert Wilf.
c    This FORTRAN77 version by John Burkardt.
c
c  Reference:
c
c    Albert Nijenhuis, Herbert Wilf,
c    Combinatorial Algorithms for Computers and Calculators,
c    Second Edition,
c    Academic Press, 1978,
c    ISBN: 0-12-519260-6,
c    LC: QA164.N54.
c
c  Parameters:
c
c    Input, integer N, the integer to be partitioned.
c
c    Input, integer TABLE(N), the number of partitions of each integer 
c    from 1 to N.  This table may be computed by I4_PARTITION_COUNT2.
c
c    Input/output, integer SEED, a seed for the random number generator.
c
c    Output, integer A(N), contains in A(1:NPART) the parts of the partition.
c
c    Output, integer MULT(N), contains in MULT(1:NPART) the multiplicity
c    of the parts.
c
c    Output, integer NPART, the number of parts in the partition chosen,
c    that is, the number of integers I with nonzero multiplicity MULT(I).
c
      implicit none

      integer n

      integer a(n)
      integer i
      integer i1
      integer id
      integer j
      integer m
      integer mult(n)
      integer npart
      double precision r8_uniform_01
      integer seed
      integer table(n)
      double precision z

      m = n
      npart = 0
      do i = 1, n
        mult(i) = 0
      end do

10    continue

      if ( 0 .lt. m ) then

        z = r8_uniform_01 ( seed )
        z = m * table(m) * z
        id = 1
        i1 = m
        j = 0

20      continue

          j = j + 1
          i1 = i1 - id

          if ( i1 .lt. 0 ) then
            id = id + 1
            i1 = m
            j = 0
            go to 20
          end if

          if ( i1 .eq. 0 ) then
            z = z - id
            if ( 0.0D+00 .lt. z ) then
              id = id + 1
              i1 = m
              j = 0
              go to 20
            else
              go to 30
            end if
          end if

          if ( 0 .lt. i1 ) then
            z = z - id * table(i1)
            if ( z .le. 0.0D+00 ) then
              go to 30
            end if
          end if

        go to 20

30      continue

        mult(id) = mult(id) + j
        npart = npart + j
        m = i1

      end if
c
c  Reformulate the partition in the standard form.
c  NPART is the number of distinct parts.
c
      npart = 0

      do i = 1, n
        if ( mult(i) .ne. 0 ) then
          npart = npart + 1
          a(npart) = i
          mult(npart) = mult(i)
        end if
      end do

      do i = npart + 1, n
        mult(i) = 0
      end do

      return
      end
      subroutine i4_partitions_next ( s, m )

c*********************************************************************72
c
cc I4_PARTITIONS_NEXT: next partition into S parts.
c
c  Discussion:
c
c    This function generates, one at a time, entries from the list of
c    nondecreasing partitions of the integers into S or fewer parts.
c
c    The list is ordered first by the integer that is partitioned
c    (the sum of the entries), and second by decreasing lexical order
c    in the partition vectors.
c
c    The first value returned is the only such partition of 0.
c
c    Next comes the only partition of 1.
c
c    There follow two partitions of 2, and so on.
c
c    Typical use of this function begins with an initialization call,
c    and then repeated calls in which the output from the previous call
c    is used as input to the next call:
c
c    m = [ 0, 0, 0 ];
c
c    while ( condition )
c      m = i4_partitions_next ( s, m );
c    end
c
c  Example:
c
c    S = 3
c
c    P  D    M
c    _  _  _____
c    1  0  0 0 0
c    2  1  1 0 0
c    3  2  2 0 0
c    4  2  1 1 0
c    5  3  3 0 0
c    6  3  2 1 0
c    7  3  1 1 1
c    8  4  4 0 0
c    9  4  3 1 0
c   10  4  2 2 0
c   11  4  2 1 1
c   12  5  5 0 0
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    03 August 2010
c
c  Author:
c
c    Original MATLAB version by Alan Genz.
c    FORTRAN77 version by John Burkardt.
c
c  Parameters:
c
c    Input, integer S, the number of items in the partition.
c
c    Input/output, integer M(S).  On input, the current partition.  
c    On first call, this should be a nondecreasing partition.  Thereafter, it 
c    should be the output partition from the previous call.  On output, the
c    next partition.
c
      implicit none

      integer s

      integer i
      integer j
      integer m(s)
      integer msum

      msum = m(1)

      do i = 2, s

        msum = msum + m(i)

        if ( m(1) .le. m(i) + 1 ) then
          m(i) = 0
        else
          m(1) = msum - ( i - 1 ) * ( m(i) + 1 )
          do j = 2, i
            m(j) = m(i) + 1
          end do
          return
        end if

      end do
c
c  If we failed to find a suitable index I, put
c  the entire sum into M(1), increment by 1, and
c  prepare to partition the next integer.
c
      m(1) = msum + 1

      return
      end
      function i4_rise ( x, n )

c*********************************************************************72
c
cc I4_RISE computes the rising factorial function [X]^N.
c
c  Discussion:
c
c    [X]^N = X * ( X + 1 ) * ( X + 2 ) * ... * ( X + N - 1 ).
c
c    Note that the number of ways of arranging N objects in M ordered
c    boxes is [M]^N.  (Here, the ordering of the objects in each box matters).  
c    Thus, 2 objects in 2 boxes have the following 6 possible arrangements:
c
c      -|12, 1|2, 12|-, -|21, 2|1, 21|-.
c
c    Moreover, the number of non-decreasing maps from a set of
c    N to a set of M ordered elements is [M]^N / N!.  Thus the set of
c    nondecreasing maps from (1,2,3) to (a,b,c,d) is the 20 elements:
c
c      aaa, abb, acc, add, aab, abc, acd, aac, abd, aad
c      bbb, bcc, bdd, bbc, bcd, bbd, ccc, cdd, ccd, ddd.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    20 January 2007
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer X, the argument of the rising factorial function.
c
c    Input, integer N, the order of the rising factorial function.
c    If N = 0, RISE = 1, if N = 1, RISE = X.  Note that if N is
c    negative, a "falling" factorial will be computed.
c
c    Output, integer I4_RISE, the value of the rising factorial 
c    function.
c
      implicit none

      integer arg
      integer i
      integer i4_rise
      integer n
      integer value
      integer x

      value = 1

      arg = x

      if ( 0 .lt. n ) then

        do i = 1, n
          value = value * arg
          arg = arg + 1
        end do

      else if ( n .lt. 0 ) then

        do i = -1, n, -1
          value = value * arg
          arg = arg - 1
        end do

      end if

      i4_rise = value

      return
      end
      subroutine i4_sqrt ( n, q, r )

c*********************************************************************72
c
cc I4_SQRT finds the integer square root of N by solving N = Q**2 + R.
c
c  Discussion:
c
c    The integer square root of N is an integer Q such that
c    Q**2 <= N but N < (Q+1)**2.
c
c    A simpler calculation would be something like
c
c      Q = int ( sqrt ( real ( N ) ) )
c
c    but this calculation has the virtue of using only integer arithmetic.
c
c    To avoid the tedium of worrying about negative arguments, the routine
c    automatically considers the absolute value of the argument.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    20 January 2007
c
c  Author:
c
c   John Burkardt
c
c  Reference:
c
c    Mark Herkommer,
c    Number Theory, A Programmer's Guide,
c    McGraw Hill, 1999, pages 294-307.
c
c  Parameters:
c
c    Input, integer N, the number whose integer square root is desired.
c    Actually, only the absolute value of N is considered.
c
c    Output, integer Q, R, the integer square root, and positive remainder,
c    of N.
c
      implicit none

      integer n
      integer n_abs
      integer q
      integer r

      n_abs = abs ( n )

      q = n_abs

      if ( 0 .lt. n_abs ) then

10      continue

        if ( ( n_abs / q ) .lt. q ) then
          q = ( q + ( n_abs / q ) ) / 2
          go to 10
        end if

      end if

      r = n_abs - q**2

      return
      end
      subroutine i4_sqrt_cf ( n, max_term, n_term, b )

c*********************************************************************72
c
cc I4_SQRT_CF finds the continued fraction representation of a square root of an integer.
c
c  Discussion:
c
c    The continued fraction representation of the square root of an integer
c    has the form
c
c      [ B0, (B1, B2, B3, ..., BM), ... ]
c
c    where
c
c      B0 = int ( sqrt ( real ( N ) ) )
c      BM = 2 * B0
c      the sequence ( B1, B2, B3, ..., BM ) repeats in the representation.
c      the value M is termed the period of the representation.
c
c  Example:
c
c     N  Period  Continued Fraction
c
c     2       1  [ 1, 2, 2, 2, ... ]
c     3       2  [ 1, 1, 2, 1, 2, 1, 2... ]
c     4       0  [ 2 ]
c     5       1  [ 2, 4, 4, 4, ... ]
c     6       2  [ 2, 2, 4, 2, 4, 2, 4, ... ]
c     7       4  [ 2, 1, 1, 1, 4, 1, 1, 4, 1, 1, 4... ]
c     8       2  [ 2, 1, 4, 1, 4, 1, 4, 1, 4, ... ]
c     9       0  [ 3 ]
c    10       1  [ 3, 6, 6, 6, ... ]
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    20 January 2007
c
c  Author:
c
c   John Burkardt
c
c  Reference:
c
c    Mark Herkommer,
c    Number Theory, A Programmer's Guide,
c    McGraw Hill, 1999, pages 294-307.
c
c  Parameters:
c
c    Input, integer N, the number whose continued fraction square root
c    is desired.
c
c    Input, integer MAX_TERM, the maximum number of terms that may
c    be computed.
c
c    Output, integer N_TERM, the number of terms computed beyond the
c    0 term.  The routine should stop if it detects that the period
c    has been reached.
c
c    Output, integer B(0:MAX_TERM), contains the continued fraction
c    coefficients in entries B(0:N_TERM).
c
      implicit none

      integer max_term

      integer b(0:max_term)
      integer n
      integer n_term
      integer p
      integer q
      integer r
      integer s

      n_term = 0

      call i4_sqrt ( n, s, r )
      b(0) = s

      if ( 0 .lt. r ) then

        p = 0
        q = 1

10      continue

          p = b(n_term) * q - p
          q = ( n - p**2 ) / q

          if ( max_term .le. n_term ) then
            return
          end if

          n_term = n_term + 1
          b(n_term) = ( p + s ) / q

          if ( q .eq. 1 ) then
            go to 20
          end if

        go to 10

20      continue

      end if

      return
      end
      subroutine i4_swap ( i, j )

c*********************************************************************72
c
cc I4_SWAP switches two I4's.
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
c  Parameters:
c
c    Input/output, integer I, J.  On output, the values of I and
c    J have been interchanged.
c
      implicit none

      integer i
      integer j
      integer k

      k = i
      i = j
      j = k

      return
      end
      subroutine i4_to_bvec ( i4, n, bvec )

c*********************************************************************72
c
cc I4_TO_BVEC makes a signed binary vector from an I4.
c
c  Discussion:
c
c    A BVEC is a vector of binary digits representing an integer.  
c
c    BVEC(1) is 0 for positive values and 1 for negative values, which
c    are stored in 2's complement form.
c
c    For positive values, BVEC(N) contains the units digit, BVEC(N-1)
c    the coefficient of 2, BVEC(N-2) the coefficient of 4 and so on,
c    so that printing the digits in order gives the binary form of the number.
c
c    Negative values have a two's complement operation applied.
c
c    To guarantee that there will be enough space for any
c    value of I, it would be necessary to set N = 32.
c
c  Example:
c
c    I4       IVEC         binary
c    --  ----------------  ------
c     1  1  0  0  0  0  1      1
c     2  0  0  0  0  1  0     10
c     3  0  0  0  0  1  1     11
c     4  0  0  0  1  0  0    100
c     9  0  0  1  0  0  1   1001
c    -9  1  1  0  1  1  1  -1001 = 110111 (2's complement)
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    20 March 2009
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer I4, an integer to be represented.
c
c    Input, integer N, the dimension of the vector.
c
c    Output, integer BVEC(N), the signed binary representation.
c
      implicit none

      integer n

      integer base
      parameter ( base = 2 )
      integer bvec(n)
      integer i4
      integer i4_copy
      integer j

      i4_copy = abs ( i4 )

      do j = n, 2, -1

        bvec(j) = mod ( i4_copy, base )

        i4_copy = i4_copy / base

      end do

      bvec(1) = 0

      if ( i4 .lt. 0 ) then
        call bvec_complement2 ( n, bvec, bvec )
      end if

      return
      end
      subroutine i4_to_chinese ( j, n, m, r )

c*********************************************************************72
c
cc I4_TO_CHINESE converts an I4 to its Chinese remainder form.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    20 January 2007
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer J, the integer to be converted.
c
c    Input, integer N, the number of moduluses.
c
c    Input, integer M(N), the moduluses.  These should be positive
c    and pairwise prime.
c
c    Output, integer R(N), the Chinese remainder representation of the integer.
c
      implicit none

      integer n

      integer i
      integer ierror
      integer i4_modp
      integer j
      integer m(n)
      integer r(n)

      call chinese_check ( n, m, ierror )

      if ( ierror .ne. 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'I4_TO_CHINESE - Fatal error!'
        write ( *, '(a)' ) '  The moduluses are not legal.'
        stop
      end if

      do i = 1, n
        r(i) = i4_modp ( j, m(i) )
      end do

      return
      end
      subroutine i4_to_dvec ( i4, n, dvec )

c*********************************************************************72
c
cc I4_TO_DVEC makes a signed decimal vector from an I4.
c
c  Discussion:
c
c    A DVEC is an integer vector of decimal digits, intended to
c    represent an integer.  DVEC(1) is the units digit, DVEC(N-1)
c    is the coefficient of 10**(N-2), and DVEC(N) contains sign
c    information.  It is 0 if the number is positive, and 9 if
c    the number is negative.
c
c    Negative values have a ten's complement operation applied.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    20 January 2007
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer I4, an integer to be represented.
c
c    Input, integer N, the dimension of the vector.
c
c    Output, integer DVEC(N), the signed decimal representation.
c
      implicit none

      integer n

      integer base
      parameter ( base = 10 )
      integer dvec(n)
      integer i
      integer i4
      integer i4_copy

      i4_copy = abs ( i4 )

      do i = 1, n-1

        dvec(i) = mod ( i4_copy, base )

        i4_copy = i4_copy / base

      end do

      dvec(n) = 0

      if ( i4 .lt. 0 ) then
        call dvec_complementx ( n, dvec, dvec )
      end if

      return
      end
      subroutine i4_to_i4poly ( intval, base, degree_max, degree, a )

c*********************************************************************72
c
cc I4_TO_I4POLY converts an I4 to an I4POLY in a given base.
c
c  Example:
c
c    INTVAL  BASE  Degree     A (in reverse order!)
c
c         1     2       0     1
c         6     2       2     1  1  0
c        23     2       5     1  0  1  1  1
c        23     3       3     2  1  2
c        23     4       3     1  1  3
c        23     5       2     4  3
c        23     6       2     3  5
c        23    23       1     1  0
c        23    24       0    23
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    20 January 2007
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer INTVAL, an integer to be converted.
c
c    Input, integer BASE, the base, which should be greater than 1.
c
c    Input, integer DEGREE_MAX, the maximum degree.
c
c    Output, integer DEGREE, the degree of the polynomial.
c
c    Output, integer A(0:DEGREE_MAX), contains the coefficients
c    of the polynomial expansion of INTVAL in base BASE.
c
      implicit none

      integer degree_max

      integer a(0:degree_max)
      integer base
      integer degree
      integer i
      integer intval
      integer j

      do i = 0, degree_max
        a(i) = 0
      end do

      j = abs ( intval )

      degree = 0

      a(degree) = mod ( j, base )

      j = j - a(degree) 
      j = j / base

10    continue

      if ( 0 .lt. j ) then

        degree = degree + 1

        if ( degree .lt. degree_max ) then
          a(degree) = mod ( j, base )
        end if

        j = j - a(degree)
        j = j / base

        go to 10

      end if

      if ( intval .lt. 0 ) then
        do i = 0, degree_max
          a(i) = -a(i)
        end do
      end if

      return
      end
      subroutine i4_to_s_left ( intval, s )

c*********************************************************************72
c
cc I4_TO_S_LEFT converts an I4 to a left-justified string.
c
c  Example:
c
c    Assume that S is 6 characters long:
c
c    INTVAL  S
c
c         1  1
c        -1  -1
c         0  0
c      1952  1952
c    123456  123456
c   1234567  ******  <-- Not enough room!
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    20 January 2007
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer INTVAL, an integer to be converted.
c
c    Output, character * ( * ) S, the representation of the integer.
c    The integer will be left-justified.  If there is not enough space,
c    the string will be filled with stars.
c 
      implicit none

      character c
      integer i
      integer idig
      integer ihi
      integer ilo
      integer intval
      integer ipos
      integer ival
      character * ( * ) s

      s = ' '

      ilo = 1
      ihi = len ( s )

      if ( ihi .le. 0 ) then
        return
      end if
c
c  Make a copy of the integer.
c
      ival = intval
c
c  Handle the negative sign.
c
      if ( ival .lt. 0 ) then

        if ( ihi .le. 1 ) then
          s(1:1) = '*'
          return
        end if

        ival = -ival
        s(1:1) = '-'
        ilo = 2

      end if
c
c  The absolute value of the integer goes into S(ILO:IHI).
c
      ipos = ihi
c
c  Find the last digit of IVAL, strip it off, and stick it into the string.
c
10    continue

        idig = mod ( ival, 10 )
        ival = ival / 10

        if ( ipos .lt. ilo ) then
          do i = 1, ihi
            s(i:i) = '*'
          end do
          return
        end if

        call digit_to_ch ( idig, c )

        s(ipos:ipos) = c
        ipos = ipos - 1

        if ( ival .eq. 0 ) then
          go to 20
        end if

      go to 10

20    continue
c
c  Shift the string to the left.
c
      s(ilo:ilo+ihi-ipos-1) = s(ipos+1:ihi)
      s(ilo+ihi-ipos:ihi) = ' '

      return
      end
      subroutine i4_to_van_der_corput ( seed, base, r )

c*********************************************************************72
c
cc I4_TO_VAN_DER_CORPUT computes an element of a van der Corput sequence.
c
c  Discussion:
c
c    The van der Corput sequence is often used to generate a "subrandom"
c    sequence of points which have a better covering property
c    than pseudorandom points.
c
c    The van der Corput sequence generates a sequence of points in [0,1]
c    which (theoretically) never repeats.  Except for SEED = 0, the
c    elements of the van der Corput sequence are strictly between 0 and 1.
c
c    The van der Corput sequence writes an integer in a given base B,
c    and then its digits are "reflected" about the decimal point.
c    This maps the numbers from 1 to N into a set of numbers in [0,1],
c    which are especially nicely distributed if N is one less
c    than a power of the base.
c
c    Hammersley suggested generating a set of N nicely distributed
c    points in two dimensions by setting the first component of the
c    Ith point to I/N, and the second to the van der Corput 
c    value of I in base 2.  
c
c    Halton suggested that in many cases, you might not know the number 
c    of points you were generating, so Hammersley's formulation was
c    not ideal.  Instead, he suggested that to generate a nicely
c    distributed sequence of points in M dimensions, you simply
c    choose the first M primes, P(1:M), and then for the J-th component of
c    the I-th point in the sequence, you compute the van der Corput
c    value of I in base P(J).
c
c    Thus, to generate a Halton sequence in a 2 dimensional space,
c    it is typical practice to generate a pair of van der Corput sequences,
c    the first with prime base 2, the second with prime base 3.
c    Similarly, by using the first K primes, a suitable sequence
c    in K-dimensional space can be generated.
c
c    The generation is quite simple.  Given an integer SEED, the expansion
c    of SEED in base BASE is generated.  Then, essentially, the result R
c    is generated by writing a decimal point followed by the digits of
c    the expansion of SEED, in reverse order.  This decimal value is actually
c    still in base BASE, so it must be properly interpreted to generate
c    a usable value.
c
c  Example:
c
c    BASE = 2
c
c    SEED     SEED      van der Corput
c    decimal  binary    binary   decimal
c    -------  ------    ------   -------
c        0  =     0  =>  .0     = 0.0
c        1  =     1  =>  .1     = 0.5
c        2  =    10  =>  .01    = 0.25
c        3  =    11  =>  .11    = 0.75
c        4  =   100  =>  .001   = 0.125
c        5  =   101  =>  .101   = 0.625
c        6  =   110  =>  .011   = 0.375
c        7  =   111  =>  .111   = 0.875
c        8  =  1000  =>  .0001  = 0.0625
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    03 January 2007
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    John Halton,
c    On the efficiency of certain quasi-random sequences of points
c    in evaluating multi-dimensional integrals,
c    Numerische Mathematik,
c    Volume 2, pages 84-90, 1960.
c 
c    John Hammersley,
c    Monte Carlo methods for solving multivariable problems,
c    Proceedings of the New York Academy of Science,
c    Volume 86, pages 844-874, 1960.
c
c    Johannes van der Corput,
c    Verteilungsfunktionen I & II,
c    Nederl. Akad. Wetensch. Proc.,
c    Volume 38, 1935, pages 813-820, pages 1058-1066.
c
c  Parameters:
c
c    Input, integer SEED, the seed or index of the desired element.
c    SEED should be nonnegative.
c    SEED = 0 is allowed, and returns R = 0.
c
c    Input, integer BASE, the van der Corput base, which is typically
c    a prime number.  BASE must be greater than 1.
c
c    Output, double precision R, the SEED-th element of the van der 
c    Corput sequence for base BASE.
c
      implicit none

      integer base
      double precision base_inv
      integer digit
      double precision r
      integer seed
      integer seed2

      if ( base .le. 1 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'I4_TO_VAN_DER_CORPUT - Fatal error!'
        write ( *, '(a)' ) '  The input base BASE is <= 1!'
        write ( *, '(a,i8)' ) '  BASE = ', base
        stop
      end if

      if ( seed .lt. 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'I4_TO_VAN_DER_CORPUT - Fatal error!'
        write ( *, '(a)' ) '  The input base SEED is < 0!'
        write ( *, '(a,i8)' ) '  SEED = ', seed
        stop
      end if

      seed2 = seed

      r = 0.0D+00

      base_inv = 1.0D+00 / dble ( base )

10    continue

      if ( seed2 .ne. 0 ) then
        digit = mod ( seed2, base )
        r = r + dble ( digit ) * base_inv
        base_inv = base_inv / dble ( base )
        seed2 = seed2 / base
        go to 10
      end if

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
      integer i4_huge
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
        seed = seed + i4_huge ( )
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
      subroutine i4mat_01_rowcolsum ( m, n, r, c, a, ierror )

c*********************************************************************72
c
cc I4MAT_01_ROWCOLSUM creates a 0/1 I4MAT with given row and column sums.
c
c  Discussion:
c
c    Given an M vector R and N vector C, there may exist one or more
c    M by N matrices with entries that are 0 or 1, whose row sums are R
c    and column sums are C.
c
c    For convenience, this routine requires that the entries of R and C
c    be given in nonincreasing order.
c
c    There are several requirements on R and C.  The simple requirements
c    are that the entries of R and C must be nonnegative, that the entries
c    of R must each be no greater than N, and those of C no greater than M,
c    and that the sum of the entries of R must equal the sum of the entries 
c    of C.
c
c    The final technical requirement is that if we form R*, the conjugate
c    partition of R, then C is majorized by R*, that is, that every partial
c    sum from 1 to K of the entries of C is no bigger than the sum of the same
c    entries of R*, for every K from 1 to N.
c
c    Given these conditions on R and C, there is at least one 0/1 matrix
c    with the given row and column sums.
c
c    The conjugate partition of R is constructed as follows:
c      R*(1) is the number of entries of R that are 1 or greater.
c      R*(2) is the number of entries of R that are 2 or greater.
c      ...
c      R*(N) is the number of entries of R that are N (can't be greater).
c
c  Example:
c
c    M = N = 5
c    R = ( 3, 2, 2, 1, 1 )
c    C = ( 2, 2, 2, 2, 1 )
c
c    A =
c      1 0 1 0 1
c      1 0 0 1 0
c      0 1 0 1 0
c      0 1 0 0 0
c      0 0 1 0 0
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    21 January 2007
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Jack van Lint, Richard Wilson,
c    A Course in Combinatorics,
c    Oxford, 1992, pages 148-156.
c
c    James Sandeson,
c    Testing Ecological Patterns,
c    American Scientist,
c    Volume 88, July-August 2000, pages 332-339.
c
c    Ian Saunders,
c    Algorithm AS 205,
c    Enumeration of R x C Tables with Repeated Row Totals,
c    Applied Statistics,
c    Volume 33, Number 3, pages 340-352, 1984.
c
c  Parameters:
c
c    Input, integer M, N, the number of rows and columns in the array.
c
c    Input, integer R(M), C(N), the row and column sums desired for the array.
c    Both vectors must be arranged in descending order.
c    The elements of R must be between 0 and N.
c    The elements of C must be between 0 and M.
c
c    Output, integer A(M,N), the M by N matrix with the given row and
c    column sums.
c    Each entry of A is 0 or 1.
c
c    Output, integer IERROR, an error flag.
c    0, no error was encountered, and the array was computed.
c    1, R and C do not have the same total.
c    2, R is not monotone decreasing, or has illegal entries.
c    3, C is not monotone decreasing, or has illegal entries.
c    4, R and C are not a possible set of row and column sums.
c
      implicit none

      integer m
      integer n

      integer a(m,n)
      integer c(n)
      integer c_sum
      integer i
      integer ierror
      logical i4vec_descends
      integer i4vec_maxloc_last
      integer i4vec_sum
      integer j
      integer k
      integer r(m)
      integer r_conj(n)
      integer r_sum
      integer r2(m)

      do j = 1, n
        do i = 1, m
          a(i,j) = 0
        end do
      end do
c
c  Check conditions.
c
      ierror = 0

      if ( i4vec_sum ( m, r ) .ne. i4vec_sum ( n, c ) ) then
        ierror = 1
        return
      end if

      if ( .not. i4vec_descends ( m, r ) ) then
        ierror = 2
        return
      end if

      if ( n .lt. r(1) .or. r(m) .lt. 0 ) then
        ierror = 2
        return
      end if

      if ( .not. i4vec_descends ( n, c ) ) then
        ierror = 3
        return
      end if

      if ( m .lt. c(1) .or. c(n) .lt. 0 ) then
        ierror = 3
        return
      end if
c
c  Compute the conjugate of R.
c
      do i = 1, n
        r_conj(i) = 0
      end do
    
      do i = 1, m
        do j = 1, r(i)
          r_conj(j) = r_conj(j) + 1
        end do
      end do
c
c  C must be majorized by R_CONJ.
c
      r_sum = 0
      c_sum = 0
      do i = 1, n
        r_sum = r_sum + r_conj(i)
        c_sum = c_sum + c(i)
        if ( r_sum .lt. c_sum ) then
          ierror = 4
          return
        end if
      end do
c
c  We need a temporary copy of R that we can decrement.
c
      do i = 1, m
        r2(i) = r(i)
      end do

      do j = n, 1, -1

        i = i4vec_maxloc_last ( m, r2 )

        do k = 1, c(j)
c
c  By adding 1 rather than setting A(I,J) to 1, we were able to spot
c  an error where the index was "sticking".
c
          a(i,j) = a(i,j) + 1

          r2(i) = r2(i) - 1
c
c  There's a special case you have to watch out for.
c  If I was 1, and when you decrement R2(1), I is going to be 1 again,
c  and you're staying in the same column, that's not good.
c
          if ( 1 .lt. i ) then
            i = i - 1
          else
            i = i4vec_maxloc_last ( m, r2 )
            if ( i .eq. 1 .and. k .lt. c(j) ) then
              i = 1 + i4vec_maxloc_last ( m-1, r2(2) )
            end if
          end if

        end do

      end do

      return
      end
      subroutine i4mat_01_rowcolsum2 ( m, n, r, c, a, ierror )

c*********************************************************************72
c
cc I4MAT_01_ROWCOLSUM2 creates a 0/1 I4MAT with given row and column sums.
c
c  Discussion:
c
c    This routine uses network flow optimization to compute the results.
c
c    Given an M vector R and N vector C, there may exist one or more
c    M by N matrices with entries that are 0 or 1, whose row sums are R
c    and column sums are C.
c
c    For convenience, this routine requires that the entries of R and C
c    be given in nonincreasing order.
c
c    There are several requirements on R and C.  The simple requirements
c    are that the entries of R and C must be nonnegative, that the entries
c    of R must each no greater than N, and those of C no greater than M,
c    and that the sum of the entries of R must equal the sum of the 
c    entries of C.
c
c    The final technical requirement is that if we form R*, the conjugate
c    partition of R, then C is majorized by R*, that is, that every partial
c    sum from 1 to K of the entries of C is no bigger than the sum of the same
c    entries of R*, for every K from 1 to N.
c
c    Given these conditions on R and C, there is at least one 0/1 matrix
c    with the given row and column sums.
c
c    The conjugate partition of R is constructed as follows:
c      R*(1) is the number of entries of R that are 1 or greater.
c      R*(2) is the number of entries of R that are 2 or greater.
c      ...
c      R*(N) is the number of entries of R that are N (can't be greater).
c
c  Example:
c
c    M = N = 5
c    R = ( 3, 2, 2, 1, 1 )
c    C = ( 2, 2, 2, 2, 1 )
c
c    A =
c      1 0 1 0 1
c      1 0 0 1 0
c      0 1 0 1 0
c      0 1 0 0 0
c      0 0 1 0 0
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    21 January 2007
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Albert Nijenhuis, Herbert Wilf,
c    Combinatorial Algorithms for Computers and Calculators,
c    Second Edition,
c    Academic Press, 1978,
c    ISBN: 0-12-519260-6,
c    LC: QA164.N54.
c
c    Jack van Lint, Richard Wilson,
c    A Course in Combinatorics,
c    Oxford, 1992, pages 148-156.
c
c    James Sandeson,
c    Testing Ecological Patterns,
c    American Scientist,
c    Volume 88, July-August 2000, pages 332-339.
c
c  Parameters:
c
c    Input, integer M, N, the number of rows and columns in the array.
c    These values do not have to be equal.
c
c    Input, integer R(M), C(N), the row and column sums desired for the array.
c    Both vectors must be arranged in descending order.
c    The elements of R must be between 0 and N.
c    The elements of C must be between 0 and M.
c    One of the conditions for a solution to exist is that the sum of the
c    elements in R equal the sum of the elements in C.
c
c    Output, integer A(M,N), the matrix with the given row and column sums.
c    Each entry of A is 0 or 1.
c
c    Output, integer IERROR, an error flag.
c    0, no error was encountered, and the array was computed.
c    1, R and C are not consistent.  A partial solution may be constructed.
c
      implicit none

      integer m
      integer n

      integer a(m,n)
      integer c(n)
      integer capflo(2,2*(m+m*n+n))
      integer cut(m+n+2)
      integer i
      integer iendpt(2,2*(m+m*n+n))
      integer ierror
      integer j
      integer k
      integer nedge
      integer nnode
      integer node_flow(m+n+2)
      integer r(m)
      integer sink
      integer source

      ierror = 0
c
c  There are M + N + 2 nodes.  The last two are the special source and sink.
c
      source = m + n + 1
      sink = m + n + 2
      nnode = m + n + 2
c
c  The source is connected to each of the R nodes.
c
      k = 0

      do i = 1, m

        k = k + 1
        iendpt(1,k) = source
        iendpt(2,k) = i
        capflo(1,k) = r(i)
        capflo(2,k) = 0

        k = k + 1
        iendpt(1,k) = i
        iendpt(2,k) = source
        capflo(1,k) = r(i)
        capflo(2,k) = 0

      end do
c
c  Every R node is connected to every C node, with capacity 1.
c
      do i = 1, m
        do j = 1, n

          k = k + 1
          iendpt(1,k) = i
          iendpt(2,k) = j+m
          capflo(1,k) = 1
          capflo(2,k) = 0

          k = k + 1
          iendpt(1,k) = j+m
          iendpt(2,k) = i
          capflo(1,k) = 1
          capflo(2,k) = 0

        end do
      end do
c
c  Every C node is connected to the sink.
c
      do j = 1, n

        k = k + 1
        iendpt(1,k) = j+m
        iendpt(2,k) = sink
        capflo(1,k) = c(j)
        capflo(2,k) = 0

        k = k + 1
        iendpt(1,k) = sink
        iendpt(2,k) = j+m
        capflo(1,k) = c(j)
        capflo(2,k) = 0

      end do
c
c  Determine the maximum flow on the network.
c
      nedge = k

      call network_flow_max ( nnode, nedge, iendpt, capflo, source, 
     &  sink, cut, node_flow )
c
c  We have a perfect solution if, and only if, the edges leading from the
c  source, and the edges leading to the sink, are all saturated.
c
      do k = 1, nedge

        i = iendpt(1,k)
        j = iendpt(2,k) - m

        if ( i .le. m .and. 1 .le. j .and. j .le. n ) then
          if ( capflo(2,k) .ne. 0 .and. capflo(2,k) .ne. 1 ) then
            ierror = 1
          end if
        end if

        if ( iendpt(1,k) .eq. source ) then
          if ( capflo(1,k) .ne. capflo(2,k) ) then
            ierror = 1
          end if
        end if

        if ( iendpt(2,k) .eq. sink ) then
          if ( capflo(1,k) .ne. capflo(2,k) ) then
            ierror = 1
          end if
        end if

      end do
c
c  If we have a solution, then A(I,J) = the flow on the edge from
c  R node I to C node J.
c
      do j = 1, n
        do i = 1, m
          a(i,j) = 0
        end do
      end do

      do k = 1, nedge

        i = iendpt(1,k)
        j = iendpt(2,k) - m

        if ( i .le. m .and. 1 .le. j .and. j .le. n ) then
          a(i,j) = capflo(2,k)
        end if

      end do

      return
      end
      subroutine i4mat_perm ( n, a, p )

c*********************************************************************72
c
cc I4MAT_PERM permutes the rows and columns of a square I4MAT.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    21 January 2007
c
c  Author:
c
c    Original FORTRAN77 version by Albert Nijenhuis, Herbert Wilf.
c    This FORTRAN77 version by John Burkardt.
c
c  Reference:
c
c    Albert Nijenhuis, Herbert Wilf,
c    Combinatorial Algorithms for Computers and Calculators,
c    Second Edition,
c    Academic Press, 1978,
c    ISBN: 0-12-519260-6,
c    LC: QA164.N54.
c
c  Parameters:
c
c    Input, integer N, the order of the matrix.
c
c    Input/output, integer A(N,N).
c    On input, the matrix to be permuted.
c    On output, the permuted matrix.
c
c    Input, integer P(N), the permutation.  P(I) is the new number of row
c    and column I.
c
      implicit none

      integer n

      integer a(n,n)
      integer i
      integer i1
      integer is
      integer it
      integer j
      integer j1
      integer j2
      integer k
      integer lc
      integer nc
      integer p(n)

      call perm_cycle ( n, 1, p, is, nc )

      do i = 1, n

        i1 = -p(i)

        if ( 0 .lt. i1 ) then

          lc = 0

10        continue

            i1 = p(i1)
            lc = lc + 1

            if ( i1 .le. 0 ) then
              go to 20
            end if

          go to 10

20        continue

          i1 = i

          do j = 1, n

            if ( p(j) .le. 0 ) then

              j2 = j
              k = lc

30            continue

                j1 = j2
                it = a(i1,j1)

40              continue

                  i1 = abs ( p(i1) )
                  j1 = abs ( p(j1) )

                  call i4_swap ( a(i1,j1), it )

                  if ( j1 .ne. j2 ) then
                    go to 40
                  end if

                  k = k - 1

                  if ( i1 .eq. i ) then
                    go to 50
                  end if

                go to 40

50              continue

                j2 = abs ( p(j2) )

                if ( k .eq. 0 ) then
                  go to 60
                end if

              go to 30

60            continue

            end if

          end do

        end if

      end do
c
c  Restore the positive signs of the data.
c
      do i = 1, n
        p(i) = abs ( p(i) )
      end do

      return
      end
      subroutine i4mat_perm2 ( m, n, a, p, q )

c*********************************************************************72
c
cc I4MAT_PERM2 permutes the rows and columns of a rectangular I4MAT.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    21 January 2007
c
c  Author:
c
c    Original FORTRAN77 version by Albert Nijenhuis, Herbert Wilf.
c    This FORTRAN77 version by John Burkardt.
c
c  Reference:
c
c    Albert Nijenhuis, Herbert Wilf,
c    Combinatorial Algorithms for Computers and Calculators,
c    Second Edition,
c    Academic Press, 1978,
c    ISBN: 0-12-519260-6,
c    LC: QA164.N54.
c
c  Parameters:
c
c    Input, integer M, the number of rows in the matrix.
c
c    Input, integer N, the number of columns in the matrix.
c
c    Input/output, integer A(M,N).
c    On input, the matrix to be permuted.
c    On output, the permuted matrix.
c
c    Input, integer P(M), the row permutation.  P(I) is the new number of row I.
c
c    Input, integer Q(N), the column permutation.  Q(I) is the new number
c    of column I.  Note that this routine allows you to pass a single array
c    as both P and Q.
c
      implicit none

      integer m
      integer n

      integer a(m,n)
      integer i
      integer i1
      integer is
      integer it
      integer j
      integer j1
      integer j2
      integer k
      integer lc
      integer nc
      integer p(m)
      integer q(n)

      call perm_cycle ( m, 1, p, is, nc )

      if ( 0 .lt. q(1) ) then
        call perm_cycle ( n, 1, q, is, nc )
      end if

      do i = 1, m

        i1 = -p(i)

        if ( 0 .lt. i1 ) then

          lc = 0

10        continue

            i1 = p(i1)
            lc = lc + 1

            if ( i1 .le. 0 ) then
              go to 20
            end if

          go to 10

20        continue

          i1 = i

          do j = 1, n

            if ( q(j) .le. 0 ) then

              j2 = j
              k = lc

30            continue

                j1 = j2
                it = a(i1,j1)

40              continue

                  i1 = abs ( p(i1) )
                  j1 = abs ( q(j1) )

                  call i4_swap ( a(i1,j1), it )

                  if ( j1 .ne. j2 ) then
                    go to 40
                  end if

                  k = k - 1

                  if ( i1 .eq. i ) then
                    go to 50
                  end if

                go to 40

50              continue

                j2 = abs ( q(j2) )
    
                if ( k .eq. 0 ) then
                  go to 60
                end if

              go to 30

60            continue

            end if

          end do

        end if

      end do
c
c  Restore the positive signs of the data.
c
      do i = 1, m
        p(i) = abs ( p(i) )
      end do

      if ( q(1) .le. 0 ) then
        do i = 1, n
          q(i) = abs ( q(i) )
        end do
      end if

      return
      end
      subroutine i4mat_print ( m, n, a, title )

c*********************************************************************72
c
cc I4MAT_PRINT prints an I4MAT.
c
c  Discussion:
c
c    An I4MAT is a rectangular array of integer values.
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
c    Input, character*(*) TITLE, a title to be printed first.
c    TITLE may be blank.
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
c    An I4MAT is a rectangular array of integer values.
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
      integer s_len_trim
      character*(*) title

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) title

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

          write ( *, '(i5,1x,10a8)' ) i, ( ctemp(j), j = 1, inc )

        end do

      end do

      write ( *, '(a)' ) ' '

      return
      end
      subroutine i4mat_u1_inverse ( n, a, b )

c*********************************************************************72
c
cc I4MAT_U1_INVERSE inverts a unit upper triangular matrix.
c
c  Discussion:
c
c    A unit upper triangular matrix is a matrix with only 1's on the main
c    diagonal, and only 0's below the main diagonal.  Above the main
c    diagonal, the entries may be assigned any value.
c
c    It may be surprising to note that the inverse of an integer unit upper
c    triangular matrix is also an integer unit upper triangular matrix.
c
c    Note that this routine can invert a matrix in place, that is, with no
c    extra storage.  If the matrix is stored in A, then the call
c
c      call i4mat_u1_inverse ( n, a, a )
c
c    will result in A being overwritten by its inverse, which can
c    save storage if the original value of A is not needed later.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    03 January 2007
c
c  Author:
c
c    Original FORTRAN77 version by Albert Nijenhuis, Herbert Wilf.
c    This FORTRAN77 version by John Burkardt.
c
c  Reference:
c
c    Albert Nijenhuis, Herbert Wilf,
c    Combinatorial Algorithms for Computers and Calculators,
c    Second Edition,
c    Academic Press, 1978,
c    ISBN: 0-12-519260-6,
c    LC: QA164.N54.
c
c  Parameters:
c
c    Input, integer N, the number of rows and columns in the matrix.
c
c    Input, integer A(N,N), the unit upper triangular matrix
c    to be inverted.
c
c    Output, integer B(N,N), the inverse matrix.
c
      implicit none

      integer n

      integer a(n,n)
      integer b(n,n)
      integer i
      integer isum
      integer j
      integer k

      do j = n, 1, -1

        do i = n, 1, -1

          if ( i .eq. j ) then
           isum = 1
          else
            isum = 0
          end if

          do k = i + 1, j
            isum = isum - a(i,k) * b(k,j)
          end do

          b(i,j) = isum

        end do
      end do

      return
      end
      subroutine i4poly ( n, a, x0, iopt, value )

c*********************************************************************72
c
cc I4POLY performs operations on I4POLY's in power or factorial form.
c
c  Discussion:
c
c    The power sum form of a polynomial is
c
c      P(X) = A1 + A2*X + A3*X**2 + ... + (AN+1)*X**N
c
c    The Taylor expansion at C has the form
c
c      P(X) = A1 + A2*(X-C) + A3*(X-C)**2 + ... + (AN+1)*(X-C)**N
c
c    The factorial form of a polynomial is
c
c      P(X) = A1 + A2*X + A3*(X)*(X-1) + A4*(X)*(X-1)*(X-2)+...
c        + (AN+1)*(X)*(X-1)*...*(X-N+1)
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    21 January 2007
c
c  Author:
c
c    Original FORTRAN77 version by Albert Nijenhuis, Herbert Wilf.
c    This FORTRAN77 version by John Burkardt.
c
c  Reference:
c
c    Albert Nijenhuis, Herbert Wilf,
c    Combinatorial Algorithms for Computers and Calculators,
c    Second Edition,
c    Academic Press, 1978,
c    ISBN: 0-12-519260-6,
c    LC: QA164.N54.
c
c  Parameters:
c
c    Input, integer N, the number of coefficients in the polynomial
c    (in other words, the polynomial degree + 1)
c
c    Input/output, integer A(N), the coefficients of the polynomial.  Depending
c    on the option chosen, these coefficients may be overwritten by those
c    of a different form of the polynomial.
c
c    Input, integer X0, for IOPT = -1, 0, or positive, the value of the
c    argument at which the polynomial is to be evaluated, or the
c    Taylor expansion is to be carried out.
c
c    Input, integer IOPT, a flag describing which algorithm is to
c    be carried out:
c    -3: Reverse Stirling.  Input the coefficients of the polynomial in
c    factorial form, output them in power sum form.
c    -2: Stirling.  Input the coefficients in power sum form, output them
c    in factorial form.
c    -1: Evaluate a polynomial which has been input in factorial form.
c    0:  Evaluate a polynomial input in power sum form.
c    1 or more:  Given the coefficients of a polynomial in
c    power sum form, compute the first IOPT coefficients of
c    the polynomial in Taylor expansion form.
c
c    Output, integer VALUE, for IOPT = -1 or 0, the value of the
c    polynomial at the point X0.
c
      implicit none

      integer n

      integer a(n)
      integer eps
      integer i
      integer iopt
      integer m
      integer n1
      integer value
      integer w
      integer x0
      integer z

      n1 = min ( n, iopt )
      n1 = max ( 1, n1 )

      if ( iopt .lt. -1 ) then
        n1 = n
      end if

      eps = mod ( max ( -iopt, 0 ), 2 )

      w = -n * eps

      if ( -2 .lt. iopt ) then
        w = w + x0
      end if

      do m = 1, n1

        value = 0
        z = w

        do i = m, n
          z = z + eps
          value = a(n+m-i) + z * value
          if ( iopt .ne. 0 .and. iopt .ne. -1 ) then
            a(n+m-i) = value
          end if
        end do

        if ( iopt .lt. 0 ) then
          w = w + 1
        end if

      end do

      return
      end
      subroutine i4poly_cyclo ( n, phi )

c*********************************************************************72
c
cc I4POLY_CYCLO computes a cyclotomic polynomial.
c
c  Discussion:
c
c    For 1 <= N, let
c
c      I = SQRT ( - 1 )
c      L = EXP ( 2 * PI * I / N )
c
c    Then the N-th cyclotomic polynomial is defined by
c
c      PHI(N;X) = Product ( 1 <= K <= N and GCD(K,N) = 1 ) ( X - L**K )
c
c    We can use the Moebius MU function to write
c
c      PHI(N;X) = Product ( mod ( D, N ) = 0 ) ( X**D - 1 )**MU(N/D)
c
c    There is a sort of inversion formula:
c
c      X**N - 1 = Product ( mod ( D, N ) = 0 ) PHI(D;X)
c
c  Example:
c
c     N  PHI
c
c     0  1
c     1  X - 1
c     2  X + 1
c     3  X**2 + X + 1
c     4  X**2 + 1
c     5  X**4 + X**3 + X**2 + X + 1
c     6  X**2 - X + 1
c     7  X**6 + X**5 + X**4 + X**3 + X**2 + X + 1
c     8  X**4 + 1
c     9  X**6 + X**3 + 1
c    10  X**4 - X**3 + X**2 - X + 1
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    21 January 2007
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Raymond Seroul,
c    Programming for Mathematicians,
c    Springer Verlag, 2000, page 269.
c
c  Parameters:
c
c    Input, integer N, the index of the cyclotomic polynomial desired.
c
c    Output, integer PHI(0:N), the N-th cyclotomic polynomial.
c
      implicit none

      integer n
      integer max_poly
      parameter ( max_poly = 100 )

      integer d
      integer den(0:max_poly)
      integer den_n
      integer factor(0:n)
      integer i
      integer mu
      integer nq
      integer nr
      integer num(0:max_poly)
      integer num_n
      integer phi(0:n)
      integer rem(0:n)

      num(0) = 1
      do i = 1, max_poly
        num(i) = 0
      end do
      num_n = 0

      den(0) = 1
      do i = 1, max_poly
        den(i) = 0
      end do
      den_n = 0

      do i = 0, n
        phi(i) = 0
      end do

      do d = 1, n
c
c  For each divisor D of N, ...
c
        if ( mod ( n, d ) .eq. 0 ) then

          call i4_moebius ( n / d, mu )
c
c  ...multiply the numerator or denominator by (X^D-1).
c
          factor(0) = -1
          do i = 1, d-1
            factor(i) = 0
          end do
          factor(d) = 1

          if ( mu .eq. + 1 ) then

            if ( max_poly .lt. num_n + d ) then
              write ( *, '(a)' ) ' '
              write ( *, '(a)' ) 'I4POLY_CYCLO - Fatal error!'
              write ( *, '(a)' ) 
     &          '  Numerator polynomial degree too high.'
              stop
            end if

            call i4poly_mul ( num_n, num, d, factor, num )

            num_n = num_n + d

          else if ( mu .eq. -1 ) then

            if ( max_poly .lt. den_n + d ) then
              write ( *, '(a)' ) ' '
              write ( *, '(a)' ) 'I4POLY_CYCLO - Fatal error!'
              write ( *, '(a)' ) 
     &          '  Denominator polynomial degree too high.'
              stop
            end if

            call i4poly_mul ( den_n, den, d, factor, den )

            den_n = den_n + d

          end if

        end if

      end do
c
c  PHI = NUM / DEN
c
      call i4poly_div ( num_n, num, den_n, den, nq, phi, nr, rem )

      return
      end
      subroutine i4poly_degree ( na, a, degree )

c*********************************************************************72
c
cc I4POLY_DEGREE returns the degree of an I4POLY.
c
c  Discussion:
c
c    The degree of a polynomial is the index of the highest power
c    of X with a nonzero coefficient.
c
c    The degree of a constant polynomial is 0.  The degree of the
c    zero polynomial is debatable, but this routine returns the
c    degree as 0.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    20 January 2007
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer NA, the dimension of A.
c
c    Input, integer A(0:NA), the coefficients of the polynomials.
c
c    Output, integer DEGREE, the degree of A.
c
      implicit none

      integer na

      integer a(0:na)
      integer degree

      degree = na

10    continue

      if ( 0 .lt. degree ) then

        if ( a(degree) .ne. 0 ) then
          return
        end if

        degree = degree - 1

        go to 10

      end if

      return
      end
      subroutine i4poly_div ( na, a, nb, b, nq, q, nr, r )

c*********************************************************************72
c
cc I4POLY_DIV computes the quotient and remainder of two I4POLY's.
c
c  Discussion:
c
c    Normally, the quotient and remainder would have rational coefficients.
c    This routine assumes that the special case applies that the quotient
c    and remainder are known beforehand to be integral.
c
c    The polynomials are assumed to be stored in power sum form.
c
c    The power sum form is:
c
c      p(x) = a(0) + a(1)*x + ... + a(n-1)*x**(n-1) + a(n)*x**(n)
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    21 January 2007
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer NA, the degree of polynomial A.
c
c    Input, integer A(0:NA), the coefficients of the polynomial to be divided.
c
c    Input, integer NB, the degree of polynomial B.
c
c    Input, integer B(0:NB), the coefficients of the divisor polynomial.
c
c    Output, integer NQ, the degree of polynomial Q.
c    If the divisor polynomial is zero, NQ is returned as -1.
c
c    Output, integer Q(0:NA-NB), contains the quotient of A/B.
c    If A and B have full degree, Q should be dimensioned Q(0:NA-NB).
c    In any case, Q(0:NA) should be enough.
c
c    Output, integer NR, the degree of polynomial R.
c    If the divisor polynomial is zero, NR is returned as -1.
c
c    Output, integer R(0:NB-1), contains the remainder of A/B.
c    If B has full degree, R should be dimensioned R(0:NB-1).
c    Otherwise, R will actually require less space.
c
      implicit none

      integer na
      integer nb

      integer a(0:na)
      integer a2(0:na)
      integer b(0:nb)
      integer i
      integer j
      integer na2
      integer nb2
      integer nq
      integer nr
      integer q(0:*)
      integer r(0:*)

      call i4poly_degree ( na, a, na2 )

      call i4poly_degree ( nb, b, nb2 )

      if ( b(nb2) .eq. 0 ) then
        nq = -1
        nr = -1
        return
      end if

      do i = 0, na2
        a2(i) = a(i)
      end do

      nq = na2 - nb2
      nr = nb2 - 1

      do i = nq, 0, -1
        q(i) = a2(i+nb2) / b(nb2)
        a2(i+nb2) = 0
        do j = 0, nb2 - 1
          a2(i+j) = a2(i+j) - q(i) * b(j)
        end do
      end do

      do i = 0, nr
        r(i) = a2(i)
      end do

      return
      end
      subroutine i4poly_mul ( na, a, nb, b, c )

c*********************************************************************72
c
cc I4POLY_MUL computes the product of two I4POLY's.
c
c  Discussion:
c
c    The polynomials are in power sum form.
c
c    The power sum form is:
c
c      p(x) = a(0) + a(1)*x + ... + a(n-1)*x**(n-1) + a(n)*x**(n)
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    21 January 2007
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer NA, the degree of polynomial A.
c
c    Input, integer A(0:NA), the coefficients of the first polynomial factor.
c
c    Input, integer NB, the degree of polynomial B.
c
c    Input, integer B(0:NB), the coefficients of the second polynomial factor.
c
c    Output, integer C(0:NA+NB), the coefficients of A * B.
c
      implicit none

      integer na
      integer nb

      integer a(0:na)
      integer b(0:nb)
      integer c(0:na+nb)
      integer d(0:na+nb)
      integer i
      integer j

      do i = 0, na+nb
        d(i) = 0
      end do

      do i = 0, na
        do j = 0, nb
          d(i+j) = d(i+j) + a(i) * b(j)
        end do
      end do

      do i = 0, na+nb
        c(i) = d(i)
      end do

      return
      end
      subroutine i4poly_print ( n, a, title )

c*********************************************************************72
c
cc I4POLY_PRINT prints an I4POLY.
c
c  Discussion:
c
c    The power sum form is:
c
c      p(x) = a(0) + a(1)*x + ... + a(n-1)*x**(n-1) + a(n)*x**(n)
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    21 January 2007
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the degree of the polynomial of A.
c
c    Input, integer A(0:N), the polynomial coefficients.
c    A(0) is the constant term and
c    A(N) is the coefficient of X**N.
c
c    Input, character * ( * ) TITLE, an optional title.
c
      implicit none

      integer n

      integer a(0:n)
      integer i
      integer mag
      integer n2
      character plus_minus
      character * ( * ) title
      integer title_length

      title_length = len_trim ( title )

      if ( 0 .lt. title_length ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) title(1:title_length)
      end if

      write ( *, '(a)' ) ' '

      call i4poly_degree ( n, a, n2 )

      if ( a(n2) .lt. 0 ) then
        plus_minus = '-'
      else
        plus_minus = ' '
      end if

      mag = abs ( a(n2) )

      if ( 2 .le. n2 ) then
        write ( *, '( ''  p(x) = '', a1, i8, '' * x ^ '', i3 )' ) 
     &    plus_minus, mag, n2
      else if ( n2 .eq. 1 ) then
        write ( *, '( ''  p(x) = '', a1, i8, '' * x'' )' ) 
     &    plus_minus, mag
      else if ( n2 .eq. 0 ) then
        write ( *, '( ''  p(x) = '', a1, i8 )' ) plus_minus, mag
      end if

      do i = n2-1, 0, -1

        if ( a(i) .lt. 0.0D+00 ) then
          plus_minus = '-'
        else
          plus_minus = '+'
        end if

        mag = abs ( a(i) )

        if ( mag .ne. 0 ) then

          if ( 2 .le. i ) then
            write ( *, ' ( ''         '', a1, i8, '' * x ^ '', i3 )' ) 
     &        plus_minus, mag, i
          else if ( i .eq. 1 ) then
            write ( *, ' ( ''         '', a1, i8, '' * x'' )' ) 
     &        plus_minus, mag
          else if ( i .eq. 0 ) then
            write ( *, ' ( ''         '', a1, i8 )' ) plus_minus, mag
          end if
        end if

      end do

      return
      end
      subroutine i4poly_to_i4 ( n, a, x, value )

c*********************************************************************72
c
cc I4POLY_TO_I4 evaluates an I4POLY.
c
c  Discussion:
c
c    The power sum form is:
c
c      p(x) = a(0) + a(1)*x + ... + a(n-1)*x**(n-1) + a(n)*x**(n)
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    21 January 2007
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the degree of the polynomial.
c
c    Input, integer A(0:N), the polynomial coefficients.
c    A(0) is the constant term and
c    A(N) is the coefficient of X**N.
c
c    Input, integer X, the point at which the polynomial is to be evaluated.
c
c    Output, integer VALUE, the value of the polynomial.
c
      implicit none

      integer n

      integer a(0:n)
      integer i
      integer value
      integer x

      value = 0

      do i = n, 0, -1
        value = value * x + a(i)
      end do

      return
      end
      function i4vec_ascends ( n, x )

c*********************************************************************72
c
cc I4VEC_ASCENDS determines if an I4VEC is (weakly) ascending.
c
c  Discussion:
c
c    An I4VEC is a vector of integer values.
c
c  Example:
c
c    X = ( -8, 1, 2, 3, 7, 7, 9 )
c
c    I4VEC_ASCENDS = TRUE
c
c    The sequence is not required to be strictly ascending.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    07 January 2007
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the size of the array.
c
c    Input, integer X(N), the array to be examined.
c
c    Output, logical I4VEC_ASCENDS, is TRUE if the entries of X ascend.
c
      implicit none

      integer n

      integer i
      logical i4vec_ascends
      integer x(n)

      i4vec_ascends = .false.

      do i = 1, n - 1
        if ( x(i+1) .lt. x(i) ) then
          return
        end if
      end do

      i4vec_ascends = .true.

      return
      end
      subroutine i4vec_backtrack ( n, maxstack, stack, x, indx, k, 
     &  nstack, ncan )

c*********************************************************************72
c
cc I4VEC_BACKTRACK supervises a backtrack search for an I4VEC.
c
c  Discussion:
c
c    The routine tries to construct an integer vector one index at a time,
c    using possible candidates as supplied by the user.
c
c    At any time, the partially constructed vector may be discovered to be
c    unsatisfactory, but the routine records information about where the
c    last arbitrary choice was made, so that the search can be
c    carried out efficiently, rather than starting out all over again.
c
c    First, call the routine with INDX = 0 so it can initialize itself.
c
c    Now, on each return from the routine, if INDX is:
c      1, you've just been handed a complete candidate vector;
c         Admire it, analyze it, do what you like.
c      2, please determine suitable candidates for position X(K).
c         Return the number of candidates in NCAN(K), adding each
c         candidate to the end of STACK, and increasing NSTACK.
c      3, you're done.  Stop calling the routine;
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    23 January 2007
c
c  Author:
c
c    Original FORTRAN77 version by Albert Nijenhuis, Herbert Wilf.
c    This FORTRAN77 version by John Burkardt.
c
c  Reference:
c
c    Albert Nijenhuis, Herbert Wilf,
c    Combinatorial Algorithms for Computers and Calculators,
c    Second Edition,
c    Academic Press, 1978,
c    ISBN: 0-12-519260-6,
c    LC: QA164.N54.
c
c  Parameters:
c
c    Input, integer N, the number of positions to be filled in the vector.
c
c    Input, integer MAXSTACK, the maximum length of the stack.
c
c    Input, integer STACK(MAXSTACK), a list of all current candidates for
c    all positions 1 through K.
c
c    Input/output, integer X(N), the partially filled in candidate vector.
c
c    Input/output, integer INDX, a communication flag.
c    On input,
c      0, to begin a backtracking search.
c      2, the requested candidates for position K have been added to 
c      STACK, and NCAN(K) was updated.
c    On output:
c      1, a complete output vector has been determined and returned in X(1:N);
c      2, candidates are needed for position X(K);
c      3, no more possible vectors exist.
c
c    Input/output, integer K, the index in X that we are trying to fill.
c
c    Input/output, integer NSTACK, the current length of the stack.
c
c    Input/output, integer NCAN(N), lists the current number of candidates for
c    all positions 1 through K.
c
      implicit none

      integer n
      integer maxstack

      integer indx
      integer k
      integer ncan(n)
      integer nstack
      integer stack(maxstack)
      integer x(n)
c
c  If this is the first call, request a candidate for position 1.
c
      if ( indx .eq. 0 ) then
        k = 1
        nstack = 0
        indx = 2
        return
      end if
c
c  Examine the stack.
c
10    continue
c
c  If there are candidates for position K, take the first available
c  one off the stack, and increment K.
c
c  This may cause K to reach the desired value of N, in which case
c  we need to signal the user that a complete set of candidates
c  is being returned.
c
        if ( 0 .lt. ncan(k) ) then

          x(k) = stack(nstack)
          nstack = nstack - 1

          ncan(k) = ncan(k) - 1

          if ( k .ne. n ) then
            k = k + 1
            indx = 2
          else
            indx = 1
          end if

          go to 20
c
c  If there are no candidates for position K, then decrement K.
c  If K is still positive, repeat the examination of the stack.
c
        else

          k = k - 1

          if ( k .le. 0 ) then
            indx = 3
            go to 20
          end if

        end if

      go to 10

20    continue

      return
      end
      subroutine i4vec_copy ( n, a1, a2 )

c*********************************************************************72
c
cc I4VEC_COPY copies an I4VEC.
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
c    02 January 2007
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the length of the vectors.
c
c    Input, integer A1(N), the vector to be copied.
c
c    Output, integer A2(N), a copy of A1.
c
      implicit none

      integer n

      integer a1(n)
      integer a2(n)
      integer i

      do i = 1, n
        a2(i) = a1(i)
      end do

      return
      end
      function i4vec_descends ( n, x )

c*********************************************************************72
c
cc I4VEC_DESCENDS determines if an I4VEC is (weakly) descending.
c
c  Discussion:
c
c    An I4VEC is a vector of integer values.
c
c  Example:
c
c    X = ( 9, 7, 7, 3, 2, 1, -8 )
c
c    I4VEC_DESCENDS = TRUE
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    19 January 2007
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the size of the array.
c
c    Input, integer X(N), the array to be examined.
c
c    Output, logical I4VEC_DESCENDS, is TRUE if the entries of X descend.
c
      implicit none

      integer n

      integer i
      logical i4vec_descends
      integer x(n)

      i4vec_descends = .false.

      do i = 1, n - 1
        if ( x(i) .lt. x(i+1) ) then
          return
        end if
      end do

      i4vec_descends = .true.

      return
      end
      subroutine i4vec_frac ( n, a, k, afrac )

c*********************************************************************72
c
cc I4VEC_FRAC searches for the K-th smallest element in an I4VEC.
c
c  Discussion:
c
c    Hoare's algorithm is used.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    23 January 2007
c
c  Parameters:
c
c    Input, integer N, the number of elements of A.
c
c    Input/output, integer A(N), array to search.  On output,
c    the elements of A have been somewhat rearranged.
c
c    Input, integer K, the fractile to be sought.  If K = 1, the
c    minimum entry is sought.  If K = N, the maximum is sought.
c    Other values of K search for the entry which is K-th in size.
c    K must be at least 1, and no greater than N.
c
c    Output, integer AFRAC, the value of the K-th fractile of A.
c
      implicit none

      integer n

      integer i
      integer a(n)
      integer afrac
      integer iryt
      integer ix
      integer j
      integer k
      integer left

      if ( n .le. 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'I4VEC_FRAC  - Fatal error!'
        write ( *, '(a,i8)' ) 
     &    '  Illegal nonpositive value of N = ', n
        stop
      end if

      if ( k .le. 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'I4VEC_FRAC  - Fatal error!'
        write ( *, '(a,i8)' ) 
     &    '  Illegal nonpositive value of K = ', k
        stop
      end if

      if ( n .lt. k ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'I4VEC_FRAC  - Fatal error!'
        write ( *, '(a,i8)' ) '  Illegal N < K, K = ', k
        stop
      end if
    
      left = 1
      iryt = n

10    continue

        if ( iryt .le. left ) then
          afrac = a(k)
          go to 60
        end if

        ix = a(k)
        i = left
        j = iryt

20      continue

          if ( j .lt. i ) then

            if ( j .lt. k ) then
              left = i
            end if

            if ( k .lt. i ) then
              iryt = j
            end if

            go to 50

          end if
c
c  Find I so that IX <= A(I).
c
30        continue

          if ( a(i) .lt. ix ) then
            i = i + 1
            go to 30
          end if
c
c  Find J so that A(J) <= IX.
c
40        continue

          if ( ix .lt. a(j) ) then
            j = j - 1
            go to 40
          end if

          if ( i .le. j ) then
            call i4_swap ( a(i), a(j) )
            i = i + 1
            j = j - 1
          end if

        go to 20

50      continue
    
      go to 10

60    continue

      return
      end
      subroutine i4vec_heap_d ( n, a )

c*********************************************************************72
c
cc I4VEC_HEAP_D reorders an I4VEC into an descending heap.
c
c  Discussion:
c
c    An I4VEC is a vector of I4's.
c
c    A descending heap is an array A with the property that, for every index J,
c    A(J) >= A(2*J) and A(J) >= A(2*J+1), (as long as the indices
c    2*J and 2*J+1 are legal).
c
c                  A(1)
c                /      \
c            A(2)         A(3)
c          /     \        /  \
c      A(4)       A(5)  A(6) A(7)
c      /  \       /   \
c    A(8) A(9) A(10) A(11)
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    23 July 2008
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Albert Nijenhuis, Herbert Wilf,
c    Combinatorial Algorithms for Computers and Calculators,
c    Academic Press, 1978,
c    ISBN: 0-12-519260-6,
c    LC: QA164.N54.
c
c  Parameters:
c
c    Input, integer N, the size of the input array.
c
c    Input/output, integer A(N).
c    On input, an unsorted array.
c    On output, the array has been reordered into a heap.
c
      implicit none

      integer n

      integer a(n)
      integer i
      integer ifree
      integer key
      integer m
c
c  Only nodes N/2 down to 1 can be "parent" nodes.
c
      do i = n/2, 1, -1
c
c  Copy the value out of the parent node.
c  Position IFREE is now "open".
c
        key = a(i)
        ifree = i

10      continue
c
c  Positions 2*IFREE and 2*IFREE + 1 are the descendants of position
c  IFREE.  (One or both may not exist because they exceed N.)
c
          m = 2 * ifree
c
c  Does the first position exist?
c
          if ( n .lt. m ) then
            go to 20
          end if
c
c  Does the second position exist?
c
          if ( m + 1 .le. n ) then
c
c  If both positions exist, take the larger of the two values,
c  and update M if necessary.
c
            if ( a(m) .lt. a(m+1) ) then
              m = m + 1
            end if

          end if
c
c  If the large descendant is larger than KEY, move it up,
c  and update IFREE, the location of the free position, and
c  consider the descendants of THIS position.
c
          if ( a(m) .le. key ) then
            go to 20
          end if

          a(ifree) = a(m)
          ifree = m

        go to 10
c
c  Once there is no more shifting to do, KEY moves into the free spot IFREE.
c
20      continue

        a(ifree) = key

      end do

      return
      end
      subroutine i4vec_indicator ( n, a )

c*********************************************************************72
c
cc I4VEC_INDICATOR sets an I4VEC to the indicator vector.
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
c    02 January 2007
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the number of elements of A.
c
c    Output, integer A(N), the array to be initialized.
c
      implicit none

      integer n

      integer a(n)
      integer i

      do i = 1, n
        a(i) = i
      end do

      return
      end
      function i4vec_index ( n, a, aval )

c*********************************************************************72
c
cc I4VEC_INDEX returns the location of the first occurrence of a given value.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    19 January 2007
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the number of entries in the vector.
c
c    Input, integer A(N), the vector to be searched.
c
c    Input, integer AVAL, the value to be indexed.
c
c    Output, integer I4VEC_INDEX, the first location in A which has the
c    value AVAL, or -1 if no such index exists.
c
      implicit none

      integer n

      integer a(n)
      integer aval
      integer i
      integer i4vec_index

      do i = 1, n
        if ( a(i) .eq. aval ) then
          i4vec_index = i
          return
        end if
      end do

      i4vec_index = -1

      return
      end
      subroutine i4vec_max ( n, a, amax )

c*********************************************************************72
c
cc I4VEC_MAX computes the maximum element of an I4VEC.
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
c    19 January 2007
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the number of entries in the array.
c
c    Input, integer A(N), the array.
c
c    Output, integer AMAX, the value of the largest entry.
c
      implicit none

      integer n

      integer a(n)
      integer amax
      integer i

      amax = a(1)

      do i = 2, n
        amax = max ( amax, a(i) )
      end do

      return
      end
      function i4vec_maxloc_last ( n, x )

c*********************************************************************72
c
cc I4VEC_MAXLOC_LAST returns the index of the last maximal I4VEC entry.
c
c  Example:
c
c    X = ( 5, 1, 2, 5, 0, 5, 3 )
c
c    I4VEC_MAXLOC_LAST = 6
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    23 January 2007
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the size of the array.
c
c    Input, integer X(N), the array to be examined.
c
c    Output, integer I4VEC_MAXLOC_LAST, the index of the last element of
c    X of maximal value.
c
      implicit none

      integer n

      integer i
      integer i4vec_maxloc_last
      integer i4vec_maxval_last
      integer x(n)

      i4vec_maxloc_last = 0

      do i = 1, n
        if ( i .eq. 1 ) then
          i4vec_maxloc_last = 1
          i4vec_maxval_last = x(1)
        else if ( i4vec_maxval_last .le. x(i) ) then
          i4vec_maxloc_last = i
          i4vec_maxval_last = x(i)
        end if
      end do

      return
      end
      subroutine i4vec_min ( n, a, amin )

c*********************************************************************72
c
cc I4VEC_MIN computes the minimum element of an I4VEC.
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
c    19 January 2007
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the number of entries in the array.
c
c    Input, integer A(N), the array.
c
c    Output, integer AMIN, the value of the smallest entry.
c
      implicit none

      integer n

      integer a(n)
      integer amin
      integer i

      amin = a(1)

      do i = 2, n
        amin = min ( amin, a(i) )
      end do

      return
      end
      function i4vec_pairwise_prime ( n, a )

c*******************************************************************************
c
cc I4VEC_PAIRWISE_PRIME checks whether an I4VEC's entries are pairwise prime.
c
c  Discussion:
c
c    An I4VEC is a vector of integer values.
c
c    Two positive integers I and J are pairwise prime if they have no common
c    factor greater than 1.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    07 January 2007
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the number of values to check.
c
c    Input, integer A(N), the vector of integers.
c
c    Output, logical I4VEC_PAIRWISE_PRIME, is TRUE if the vector of integers
c    is pairwise prime.
c
      implicit none

      integer n

      integer a(n)
      integer i
      integer i4_gcd
      logical i4vec_pairwise_prime
      integer j

      i4vec_pairwise_prime = .false.

      do i = 1, n
        do j = i + 1, n
          if ( i4_gcd ( a(i), a(j) ) .ne. 1 ) then
            return
          end if
        end do
      end do

      i4vec_pairwise_prime = .true.

      return
      end
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
      integer s_len_trim
      character*(*) title
      integer title_length

      title_length = s_len_trim ( title )
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) title(1:title_length)
      write ( *, '(a)' ) ' '
      do i = 1, n
        write ( *, '(2x,i8,2x,i12)' ) i, a(i)
      end do

      return
      end
      subroutine i4vec_print_some ( n, a, max_print, title )

c*********************************************************************72
c
cc I4VEC_PRINT_SOME prints "some" of an I4VEC.
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
c    17 December 2001
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the number of entries of the vector.
c
c    Input, integer A(N), the vector to be printed.
c
c    Input, integer MAX_PRINT, the maximum number of lines to print.
c
c    Input, character*(*) TITLE, a title.
c
      implicit none

      integer n

      integer a(n)
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
          write ( *, '(i8,2x,i10)' ) i, a(i)
        end do

      else if ( 3 .le. max_print ) then

        do i = 1, max_print-2
          write ( *, '(i8,2x,i10)' ) i, a(i)
        end do
        write ( *, '(a)' ) '........  ..............'
        i = n
        write ( *, '(i8,2x,i10)' ) i, a(i)

      else

        do i = 1, max_print - 1
          write ( *, '(i8,2x,i10)' ) i, a(i)
        end do
        i = max_print
        write ( *, '(i8,2x,i10,2x,a)' ) i, a(i), '...more entries...'

      end if

      return
      end
      function i4vec_product ( n, a )

c*********************************************************************72
c
cc I4VEC_PRODUCT returns the product of the entries of an I4VEC.
c
c  Discussion:
c
c    An I4VEC is a vector of integer values.
c
c    In FORTRAN90, this facility is offered by the built in
c    PRODUCT function:
c
c      I4VEC_PRODUCT ( N, A ) = PRODUCT ( A(1:N) )
c
c    In MATLAB, this facility is offered by the built in
c    PROD function:
c
c      I4VEC_PRODUCT ( N, A ) = PROD ( A(1:N) )
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    29 January 2007
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the number of entries in the array.
c
c    Input, integer A(N), the array.
c
c    Output, integer I4VEC_PRODUCT, the product of the entries.
c
      implicit none

      integer n

      integer a(n)
      integer i
      integer i4vec_product

      i4vec_product = 1
      do i = 1, n
        i4vec_product = i4vec_product * a(i)
      end do

      return
      end
      subroutine i4vec_reverse ( n, a )

c*********************************************************************72
c
cc I4VEC_REVERSE reverses the elements of an I4VEC.
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
c    19 January 2007
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

      do i = 1, n/2
        call i4_swap ( a(i), a(n+1-i) )
      end do

      return
      end
      subroutine i4vec_sort_bubble_a ( n, a )

c*********************************************************************72
c
cc I4VEC_SORT_BUBBLE_A ascending sorts an I4VEC using bubble sort.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    19 January 2007
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the number of entries in the array.
c
c    Input/output, integer A(N).
c    On input, the array to be sorted;
c    On output, the array has been sorted.
c
      implicit none

      integer n

      integer a(n)
      integer i
      integer j

      do i = 1, n-1
        do j = i+1, n
          if ( a(j) .lt. a(i) ) then
            call i4_swap ( a(i), a(j) )
          end if
        end do
      end do

      return
      end
      subroutine i4vec_sort_heap_a ( n, a )

c*********************************************************************72
c
cc I4VEC_SORT_HEAP_A ascending sorts an I4VEC using heap sort.
c
c  Discussion:
c
c    An I4VEC is a vector of I4's.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    23 July 2008
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Albert Nijenhuis, Herbert Wilf,
c    Combinatorial Algorithms for Computers and Calculators,
c    Academic Press, 1978,
c    ISBN: 0-12-519260-6,
c    LC: QA164.N54.
c
c  Parameters:
c
c    Input, integer N, the number of entries in the array.
c
c    Input/output, integer A(N).
c    On input, the array to be sorted;
c    On output, the array has been sorted.
c
      implicit none

      integer n

      integer a(n)
      integer n1

      if ( n .le. 1 ) then
        return
      end if
c
c  1: Put A into descending heap form.
c
      call i4vec_heap_d ( n, a )
c
c  2: Sort A.
c
c  The largest object in the heap is in A(1).
c  Move it to position A(N).
c
      call i4_swap ( a(1), a(n) )
c
c  Consider the diminished heap of size N1.
c
      do n1 = n - 1, 2, -1
c
c  Restore the heap structure of A(1) through A(N1).
c
        call i4vec_heap_d ( n1, a )
c
c  Take the largest object from A(1) and move it to A(N1).
c
        call i4_swap ( a(1), a(n1) )

      end do

      return
      end
      subroutine i4vec_sort_heap_index_d ( n, a, indx )

c*********************************************************************72
c
cc I4VEC_SORT_HEAP_INDEX_D does an indexed heap descending sort of an I4VEC.
c
c  Discussion:
c
c    The sorting is not actually carried out.  Rather an index array is
c    created which defines the sorting.  This array may be used to sort
c    or index the array, or to sort or index related arrays keyed on the
c    original array.
c
c    Once the index array is computed, the sorting can be carried out
c    "implicitly:
c
c      A(INDX(I)), I = 1 to N is sorted,
c
c    after which A(I), I = 1 to N is sorted.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    23 January 2007
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the number of entries in the array.
c
c    Input, integer A(N), an array to be index-sorted.
c
c    Output, integer INDX(N), contains the sort index.  The
c    I-th element of the sorted array is A(INDX(I)).
c
      implicit none

      integer n

      integer a(n)
      integer aval
      integer i
      integer indx(n)
      integer indxt
      integer ir
      integer j
      integer l

      call i4vec_indicator ( n, indx )

      l = n / 2 + 1
      ir = n

10    continue

        if ( 1 .lt. l ) then

          l = l - 1
          indxt = indx(l)
          aval = a(indxt)

        else

          indxt = indx(ir)
          aval = a(indxt)
          indx(ir) = indx(1)
          ir = ir - 1

          if ( ir .eq. 1 ) then
            indx(1) = indxt
            go to 30
          end if

        end if

        i = l
        j = l + l

20      continue

        if ( j .le. ir ) then

          if ( j .lt. ir ) then
            if ( a(indx(j+1)) .lt. a(indx(j)) ) then
              j = j + 1
            end if
          end if
    
          if ( a(indx(j)) .lt. aval ) then
            indx(i) = indx(j)
            i = j
            j = j + j
          else
            j = ir + 1
          end if

          go to 20

        end if

        indx(i) = indxt

      go to 10

30    continue

      return
      end
      function i4vec_sum ( n, a )

c*********************************************************************72
c
cc I4VEC_SUM returns the sum of the entries of an I4VEC.
c
c  Discussion:
c
c    An I4VEC is a vector of integer values.
c
c    In FORTRAN90, this facility is offered by the built in
c    SUM function:
c
c      I4VEC_SUM ( N, A ) = SUM ( A(1:N) )
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    03 January 2007
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the number of entries in the array.
c
c    Input, integer A(N), the array.
c
c    Output, integer I4VEC_SUM, the sum of the entries.
c
      implicit none

      integer n

      integer a(n)
      integer i
      integer i4vec_sum

      i4vec_sum = 0

      do i = 1, n
        i4vec_sum = i4vec_sum + a(i)
      end do

      return
      end
      subroutine i4vec_transpose_print ( n, a, title )

c*********************************************************************72
c
cc I4VEC_TRANSPOSE_PRINT prints an I4VEC "transposed".
c
c  Example:
c
c    A = (/ 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11 /)
c    TITLE = 'My vector:  '
c
c    My vector:      1    2    3    4    5
c                    6    7    8    9   10
c                   11
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    23 January 2007
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
c    Input, character * ( * ) TITLE, a title to be printed first.
c    TITLE may be blank.
c
      implicit none

      integer n

      integer a(n)
      integer i
      integer ihi
      integer ilo
      character * ( 11 ) string
      character * ( * ) title
      integer title_len

      if ( 0 .lt. len_trim ( title ) ) then

        title_len = len ( title )

        write ( string, '(a,i3,a)' ) '(', title_len, 'x,5i12)'

        do ilo = 1, n, 5
          ihi = min ( ilo + 5 - 1, n )
          if ( ilo .eq. 1 ) then
            write ( *, '(a, 5i12)' ) title, ( a(i), i = ilo, ihi )
          else
            write ( *, string      )        ( a(i), i = ilo, ihi )
          end if
        end do

      else

        do ilo = 1, n, 5
          ihi = min ( ilo + 5 - 1, n )
          write ( *, '(5i12)' ) ( a(i), i = ilo, ihi )
        end do

      end if

      return
      end
      subroutine i4vec_uniform ( n, a, b, seed, x )

c*********************************************************************72
c
cc I4VEC_UNIFORM returns a scaled pseudorandom I4VEC.
c
c  Discussion:
c
c    An I4VEC is a vector of integer values.
c
c    The pseudorandom numbers should be uniformly distributed
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
c  Parameters:
c
c    Input, integer N, the dimension of the vector.
c
c    Input, integer A, B, the limits of the interval.
c
c    Input/output, integer SEED, the "seed" value, which should NOT be 0.
c    On output, SEED has been updated.
c
c    Output, integer X(N), a vector of numbers between A and B.
c
      implicit none

      integer n

      integer a
      integer b
      integer i
      integer i4_huge
      integer k
      real r
      integer seed
      integer value
      integer x(n)

      if ( seed .eq. 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'I4VEC_UNIFORM - Fatal error!'
        write ( *, '(a)' ) '  Input value of SEED = 0.'
        stop
      end if

      do i = 1, n

        k = seed / 127773

        seed = 16807 * ( seed - k * 127773 ) - k * 2836

        if ( seed .lt. 0 ) then
          seed = seed + i4_huge ( )
        end if

        r = real ( seed ) * 4.656612875E-10
c
c  Scale R to lie between A-0.5 and B+0.5.
c
        r = ( 1.0E+00 - r ) * ( real ( min ( a, b ) ) - 0.5E+00 )
     &    +             r   * ( real ( max ( a, b ) ) + 0.5E+00 )
c
c  Use rounding to convert R to an integer between A and B.
c
        value = nint ( r )

        value = max ( value, min ( a, b ) )
        value = min ( value, max ( a, b ) )

        x(i) = value

      end do

      return
      end
      subroutine i4vec_zero ( n, a )

c*********************************************************************72
c
cc I4VEC_ZERO sets the entries of an I4VEC to 0.
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
c    18 January 2007
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the number of entries in the vector.
c
c    Output, integer A(N), the vector, which has been set to zero.
c
      implicit none

      integer n

      integer a(n)
      integer i

      do i = 1, n
        a(i) = 0
      end do

      return
      end
      subroutine index_box2_next_2d ( n1, n2, ic, jc, i, j, more )

c*********************************************************************72
c
cc INDEX_BOX2_NEXT_2D produces index vectors on the surface of a box in 2D.
c
c  Discussion:
c
c    The box is has center at (IC,JC), and has half-widths N1 and N2.
c    The index vectors are exactly those which are between (IC-N1,JC-N1) and
c    (IC+N1,JC+N2) with the property that at least one of I and J
c    is an "extreme" value.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    23 January 2007
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N1, N2, the half-widths of the box, that is, the
c    maximum distance allowed between (IC,JC) and (I,J).
c
c    Input, integer IC, JC, the central cell of the box.
c
c    Input/output, integer I, J.  On input, the previous index set.
c    On output, the next index set.  On the first call, MORE should
c    be set to FALSE, and the input values of I and J are ignored.
c
c    Input/output, logical MORE.
c    On the first call for a given box, the user should set MORE to FALSE.
c    On return, the routine sets MORE to TRUE.
c    When there are no more indices, the routine sets MORE to FALSE.
c
      implicit none

      integer i
      integer ic
      integer j
      integer jc
      logical more
      integer n1
      integer n2

      if ( .not. more ) then
        more = .true.
        i = ic - n1
        j = jc - n2
        return
      end if

      if ( i .eq. ic + n1 .and. j .eq. jc + n2 ) then
        more = .false.
        return
      end if
c
c  Increment J.
c
      j = j + 1
c
c  Check J.
c
      if ( jc + n2 .lt. j ) then
        j = jc - n2
        i = i + 1
      else if ( j .lt. jc + n2 .and. 
     &  ( i .eq. ic - n1 .or. i .eq. ic + n1 ) ) then
        return
      else
        j = jc + n2
        return
      end if

      return
      end
      subroutine index_box2_next_3d ( n1, n2, n3, ic, jc, kc, i, 
     &  j, k, more )

c*********************************************************************72
c
cc INDEX_BOX2_NEXT_3D produces index vectors on the surface of a box in 3D.
c
c  Discussion:
c
c    The box has a central cell of (IC,JC,KC), with a half widths of
c    (N1,N2,N3).  The index vectors are exactly those between
c    (IC-N1,JC-N2,KC-N3) and (IC+N1,JC+N2,KC+N3) with the property that 
c    at least one of I, J, and K is an "extreme" value.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    23 January 2007
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N1, N2, N3, the "half widths" of the box, that is, the
c    maximum distances from the central cell allowed for I, J and K.
c
c    Input, integer IC, JC, KC, the central cell of the box.
c
c    Input/output, integer I, J, K.  On input, the previous index set.
c    On output, the next index set.  On the first call, MORE should
c    be set to FALSE, and the input values of I, J, and K are ignored.
c
c    Input/output, logical MORE.
c    On the first call for a given box, the user should set MORE to FALSE.
c    On return, the routine sets MORE to TRUE.
c    When there are no more indices, the routine sets MORE to FALSE.
c
      implicit none

      integer i
      integer ic
      integer j
      integer jc
      integer k
      integer kc
      logical more
      integer n1
      integer n2
      integer n3

      if ( .not. more ) then
        more = .true.
        i = ic - n1
        j = jc - n2
        k = kc - n3
        return
      end if

      if ( i .eq. ic + n1 .and. 
     &     j .eq. jc + n2 .and. 
     &     k .eq. kc + n3 ) then
        more = .false.
        return
      end if
c
c  Increment K.
c
      k = k + 1
c
c  Check K.
c
      if ( kc + n3 .lt. k ) then
        k = kc - n3
        j = j + 1
      else if ( k .lt. kc + n3 .and. 
     &  ( i .eq. ic - n1 .or. i .eq. ic + n1 .or. 
     &    j .eq. jc - n2 .or. j .eq. jc + n2 ) ) then
        return
      else
        k = kc + n3
        return
      end if
c
c  Check J.
c
      if ( jc + n2 .lt. j ) then
        j = jc - n2
        i = i + 1
      else if ( j .lt. jc + n2 .and. 
     &  ( i .eq. ic - n1 .or. i .eq. ic + n1 .or. 
     &    k .eq. kc - n3 .or. k .eq. kc + n3 ) ) then
        return
      else
        j = jc + n2
        return
      end if

      return
      end
      subroutine index_box_next_2d ( n1, n2, i, j, more )

c*********************************************************************72
c
cc INDEX_BOX_NEXT_2D produces index vectors on the surface of a box in 2D.
c
c  Discussion:
c
c    The index vectors are exactly those which are between (1,1) and
c    (N1,N2) with the property that at least one of I and J
c    is an "extreme" value.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    23 January 2007
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N1, N2, the "dimensions" of the box, that is, the
c    maximum values allowed for I and J.  The minimum values are
c    assumed to be 1.
c
c    Input/output, integer I, J.  On input, the previous index set.
c    On output, the next index set.  On the first call, MORE should
c    be set to FALSE, and the input values of I and J are ignored.
c
c    Input/output, logical MORE.
c    On the first call for a given box, the user should set MORE to FALSE.
c    On return, the routine sets MORE to TRUE.
c    When there are no more indices, the routine sets MORE to FALSE.
c
      implicit none

      integer i
      integer j
      logical more
      integer n1
      integer n2

      if ( .not. more ) then
        more = .true.
        i = 1
        j = 1
        return
      end if

      if ( i .eq. n1 .and. j .eq. n2 ) then
        more = .false.
        return
      end if
c
c  Increment J.
c
      j = j + 1
c
c  Check J.
c
      if ( n2 .lt. j ) then
        j = 1
        i = i + 1
      else if ( j .lt. n2 .and. ( i .eq. 1 .or. i .eq. n1 ) ) then
        return
      else
        j = n2
        return
      end if

      return
      end
      subroutine index_box_next_3d ( n1, n2, n3, i, j, k, more )

c*********************************************************************72
c
cc INDEX_BOX_NEXT_3D produces index vectors on the surface of a box in 3D.
c
c  Discussion:
c
c    The index vectors are exactly those which are between (1,1,1) and
c    (N1,N2,N3) with the property that at least one of I, J, and K
c    is an "extreme" value.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    23 January 2007
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N1, N2, N3, the "dimensions" of the box, that is, the
c    maximum values allowed for I, J and K.  The minimum values are
c    assumed to be 1.
c
c    Input/output, integer I, J, K.  On input, the previous index set.
c    On output, the next index set.  On the first call, MORE should
c    be set to FALSE, and the input values of I, J, and K are ignored.
c
c    Input/output, logical MORE.
c    On the first call for a given box, the user should set MORE to FALSE.
c    On return, the routine sets MORE to TRUE.
c    When there are no more indices, the routine sets MORE to FALSE.
c
      implicit none

      integer i
      integer j
      integer k
      logical more
      integer n1
      integer n2
      integer n3

      if ( .not. more ) then
        more = .true.
        i = 1
        j = 1
        k = 1
        return
      end if

      if ( i .eq. n1 .and. j .eq. n2 .and. k .eq. n3 ) then
        more = .false.
        return
      end if
c
c  Increment K.
c
      k = k + 1
c
c  Check K.
c
      if ( n3 .lt. k ) then
        k = 1
        j = j + 1
      else if ( k .lt. n3 .and. 
     &  ( i .eq. 1 .or. i .eq. n1 .or. j .eq. 1 .or. j .eq. n2 ) ) then
        return
      else
        k = n3
        return
      end if
c
c  Check J.
c
      if ( n2 .lt. j ) then
        j = 1
        i = i + 1
      else if ( j .lt. n2 .and. 
     &  ( i .eq. 1 .or. i .eq. n1 .or. k .eq. 1 .or. k .eq. n3 ) ) then
        return
      else
        j = n2
        return
      end if

      return
      end
      subroutine index_next0 ( n, hi, a, more )

c*********************************************************************72
c
cc INDEX_NEXT0 generates all index vectors within given upper limits.
c
c  Discussion:
c
c    The index vectors are generated in such a way that the reversed
c    sequences are produced in lexicographic order.
c
c  Example:
c
c    N = 3,
c    HI = 3
c
c    1   2   3
c    ---------
c    1   1   1
c    2   1   1
c    3   1   1
c    1   2   1
c    2   2   1
c    3   2   1
c    1   3   1
c    2   3   1
c    3   3   1
c    1   1   2
c    2   1   2
c    3   1   2
c    1   2   2
c    2   2   2
c    3   2   2
c    1   3   2
c    2   3   2
c    3   3   2
c    1   1   3
c    2   1   3
c    3   1   3
c    1   2   3
c    2   2   3
c    3   2   3
c    1   3   3
c    2   3   3
c    3   3   3
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    23 January 2007
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the number of entries in A.
c
c    Input, integer HI, the upper limit for the array indices.
c    The lower limit is implicitly 1 and HI must be at least 1.
c
c    Input/output, integer A(N).
c    On startup calls, with MORE = FALSE, the input value of A
c    doesn't matter, because the routine initializes it.
c    On calls with MORE = TRUE, the input value of A must be
c    the output value of A from the previous call.  (In other words,
c    just leave it alonec).
c    On output, A contains the successor set of indices to the input
c    value.
c
c    Input/output, logical MORE.  Set this variable FALSE before
c    the first call.  Normally, MORE will be returned TRUE but
c    once all the vectors have been generated, MORE will be
c    reset to FALSE and you should stop calling the program.
c
      implicit none

      integer n

      integer a(n)
      integer hi
      integer i
      integer inc
      logical more

      if ( .not. more ) then

        do i = 1, n
          a(i) = 1
        end do

        if ( hi .lt. 1 ) then
          more = .false.
          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) 'INDEX_NEXT0 - Fatal error!'
          write ( *, '(a,i8)' ) '  HI is ', hi
          write ( *, '(a)' ) '  but HI must be at least 1.'
          stop
        end if

      else

        inc = 1

10      continue

        if ( hi .le. a(inc) ) then
          a(inc) = 1
          inc = inc + 1
          go to 10
        end if

        a(inc) = a(inc) + 1

      end if
c
c  See if there are more entries to compute.
c
      more = .false.

      do i = 1, n
        if ( a(i) .lt. hi ) then
          more = .true.
        end if
      end do

      return
      end
      subroutine index_next1 ( n, hi, a, more )

c*********************************************************************72
c
cc INDEX_NEXT1 generates all index vectors within given upper limits.
c
c  Discussion:
c
c    The index vectors are generated in such a way that the reversed
c    sequences are produced in lexicographic order.
c
c  Example:
c
c    N = 3,
c    HI(1) = 4, HI(2) = 2, HI(3) = 3
c
c    1   2   3
c    ---------
c    1   1   1
c    2   1   1
c    3   1   1
c    4   1   1
c    1   2   1
c    2   2   1
c    3   2   1
c    4   2   1
c    1   1   2
c    2   1   2
c    3   1   2
c    4   1   2
c    1   2   2
c    2   2   2
c    3   2   2
c    4   2   2
c    1   1   3
c    2   1   3
c    3   1   3
c    4   1   3
c    1   2   3
c    2   2   3
c    3   2   3
c    4   2   3
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    23 January 2007
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the number of entries in A.
c
c    Input, integer HI(N), the upper limits for the array indices.
c    The lower limit is implicitly 1, and each HI(I) should be at least 1.
c
c    Input/output, integer A(N).
c    On startup calls, with MORE = FALSE, the input value of A
c    doesn't matter, because the routine initializes it.
c    On calls with MORE = TRUE, the input value of A must be
c    the output value of A from the previous call.  (In other words,
c    just leave it alonec).
c    On output, A contains the successor set of indices to the input
c    value.
c
c    Input/output, logical MORE.  Set this variable FALSE before
c    the first call.  Normally, MORE will be returned TRUE but
c    once all the vectors have been generated, MORE will be
c    reset FALSE and you should stop calling the program.
c
      implicit none

      integer n

      integer a(n)
      integer hi(n)
      integer i
      integer inc
      logical more

      if ( .not. more ) then

        do i = 1, n
          a(i) = 1
        end do

        do i = 1, n
          if ( hi(i) .lt. 1 ) then
            more = .false.
            write ( *, '(a)' ) ' '
            write ( *, '(a)' ) 'INDEX_NEXT1 - Fatal error!'
            write ( *, '(a,i8,a,i8)' ) 
     &        '  Entry ', i, ' of HI is ', hi(i)
            write ( *, '(a)' ) '  but all entries must be at least 1.'
            stop
          end if
        end do

      else

        inc = 1

10      continue

        if ( hi(inc) .le. a(inc) ) then
          a(inc) = 1
          inc = inc + 1
          go to 10
        end if

        a(inc) = a(inc) + 1

      end if
c
c  See if there are more entries to compute.
c
      more = .false.

      do i = 1, n
        if ( a(i) .lt. hi(i) ) then
          more = .true.
        end if
      end do

      return
      end
      subroutine index_next2 ( n, lo, hi, a, more )

c*********************************************************************72
c
cc INDEX_NEXT2 generates all index vectors within given lower and upper limits.
c
c  Example:
c
c    N = 3,
c    LO(1) = 1, LO(2) = 10, LO(3) = 4
c    HI(1) = 2, HI(2) = 11, HI(3) = 6
c
c    1   2   3
c    ---------
c    1  10   4
c    2  10   4
c    1  11   4
c    2  11   4
c    1  10   5
c    2  10   5
c    1  11   5
c    2  11   5
c    1  10   6
c    2  10   6
c    1  11   6
c    2  11   6
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    23 January 2007
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the number of entries in A.  The rank of
c    the object being indexed.
c
c    Input, integer LO(N), HI(N), the lower and upper limits for the array
c    indices.  LO(I) should be less than or equal to HI(I), for each I.
c
c    Input/output, integer A(N).
c    On startup calls, with MORE = FALSE, the input value of A
c    doesn't matter, because the routine initializes it.
c    On calls with MORE = TRUE, the input value of A must be
c    the output value of A from the previous call.  (In other words,
c    just leave it alonec).
c    On output, A contains the successor set of indices to the input
c    value.
c
c    Input/output, logical MORE.  Set this variable FALSE before
c    the first call.  Normally, MORE will be returned TRUE but
c    once all the vectors have been generated, MORE will be
c    reset FALSE and you should stop calling the program.
c
      implicit none

      integer n

      integer a(n)
      integer hi(n)
      integer i
      integer inc
      integer lo(n)
      logical more

      if ( .not. more ) then

        do i = 1, n
          a(i) = lo(i)
        end do

        do i = 1, n
          if ( hi(i) .lt. lo(i) ) then
            more = .false.
            write ( *, '(a)' ) ' '
            write ( *, '(a)' ) 'INDEX_NEXT2 - Fatal error!'
            write ( *, '(a,i8,a,i8)' ) 
     &        '  Entry ', i, ' of HI is ', hi(i)
            write ( *, '(a,i8,a,i8)' ) 
     &        '  Entry ', i, ' of LO is ', lo(i)
            write ( *, '(a)' ) '  but LO(I) <= HI(I) is required.'
            stop
          end if
        end do

      else

        inc = 1

10      continue

        if ( hi(inc) .le. a(inc) ) then
          a(inc) = lo(inc)
          inc = inc + 1
          go to 10
        end if

        a(inc) = a(inc) + 1

      end if
c
c  See if there are more entries to compute.
c
      more = .false.

      do i = 1, n
        if ( a(i) .lt. hi(i) ) then
          more = .true.
        end if
      end do

      return
      end
      subroutine index_rank0 ( n, hi, a, rank )

c*********************************************************************72
c
cc INDEX_RANK0 ranks an index vector within given upper limits.
c
c  Example:
c
c    N = 3,
c    HI = 3
c    A = ( 3, 1, 2 )
c
c    RANK = 12
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    23 January 2007
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the number of entries in A.
c
c    Input, integer HI, the upper limit for the array indices.
c    The lower limit is implicitly 1, and HI should be at least 1.
c
c    Input, integer A(N), the index vector to be ranked.
c
c    Output, integer RANK, the rank of the index vector, or -1 if A
c    is not a legal index.
c
      implicit none

      integer n

      integer a(n)
      integer hi
      integer i
      integer range
      integer rank

      rank = -1
      do i = 1, n
        if ( a(i) .lt. 1 .or. hi .lt. a(i) ) then
          return
        end if
      end do

      rank = 0
      do i = n, 1, -1
        rank = hi * rank + a(i)
      end do

      rank = 1
      range = 1
      do i = 1, n
        rank = rank + ( a(i) - 1 ) * range
        range = range * hi
      end do

      return
      end
      subroutine index_rank1 ( n, hi, a, rank )

c*********************************************************************72
c
cc INDEX_RANK1 ranks an index vector within given upper limits.
c
c  Example:
c
c    N = 3,
c    HI(1) = 4, HI(2) = 2, HI(3) = 3
c    A = ( 4, 1, 2 )
c
c    RANK = 12
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    24 January 2007
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the number of entries in A.
c
c    Input, integer HI(N), the upper limits for the array indices.
c    The lower limit is implicitly 1, and each HI(I) should be at least 1.
c
c    Input, integer A(N), the index to be ranked.
c
c    Output, integer RANK, the rank of the index vector, or -1 if A
c    is not a legal index.
c
      implicit none

      integer n

      integer a(n)
      integer hi(n)
      integer i
      integer range
      integer rank

      rank = -1
      do i = 1, n
        if ( a(i) .lt. 1 .or. hi(i) .lt. a(i) ) then
          return
        end if
      end do

      rank = 0
      do i = n, 1, -1
        rank = hi(i) * rank + a(i)
      end do

      rank = 1
      range = 1
      do i = 1, n
        rank = rank + ( a(i) - 1 ) * range
        range = range * hi(i)
      end do

      return
      end
      subroutine index_rank2 ( n, lo, hi, a, rank )

c*********************************************************************72
c
cc INDEX_RANK2 ranks an index vector within given lower and upper limits.
c
c  Example:
c
c    N = 3,
c    LO(1) = 1, LO(2) = 10, LO(3) = 4
c    HI(1) = 2, HI(2) = 11, HI(3) = 6
c    A = ( 1, 11, 5 )
c
c    RANK = 7
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    24 January 2007
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the number of entries in A.
c
c    Input, integer LO(N), HI(N), the lower and upper limits for the array
c    indices.  LO(I) should be less than or equal to HI(I), for each I.
c
c    Input, integer A(N), the index vector to be ranked.
c
c    Output, integer RANK, the rank of the index vector, or -1 if A
c    is not a legal index vector.
c
      implicit none

      integer n

      integer a(n)
      integer hi(n)
      integer i
      integer lo(n)
      integer range
      integer rank

      do i = 1, n
        if ( a(i) .lt. lo(i) .or. hi(i) .lt. a(i) ) then
          rank = -1
          return
        end if
      end do

      rank = 1
      range = 1
      do i = 1, n
        rank = rank + ( a(i) - lo(i) ) * range
        range = range * ( hi(i) + 1 - lo(i) )
      end do

      return
      end
      subroutine index_unrank0 ( n, hi, rank, a )

c*********************************************************************72
c
cc INDEX_UNRANK0 unranks an index vector within given upper limits.
c
c  Example:
c
c    N = 3,
c    HI = 3
c    RANK = 12
c
c    A = ( 3, 1, 2 )
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    23 January 2007
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the number of entries in A.
c
c    Input, integer HI, the upper limit for the array indices.
c    The lower limit is implicitly 1, and HI should be at least 1.
c
c    Input, integer RANK, the rank of the desired index vector.
c
c    Output, integer A(N), the index vector of the given rank.
c
      implicit none

      integer n

      integer a(n)
      integer hi
      integer i
      integer j
      integer k
      integer range
      integer rank

      do i = 1, n
        a(i) = 0
      end do
c
c  The rank might be too small.
c
      if ( rank .lt. 1 ) then
        return
      end if

      range = hi**n
c
c  The rank might be too large.
c
      if ( range .lt. rank ) then
        return
      end if

      k = rank - 1
      do i = n, 1, -1
        range = range / hi
        j = k / range
        a(i) = j + 1
        k = k - j * range
      end do

      return
      end
      subroutine index_unrank1 ( n, hi, rank, a )

c*********************************************************************72
c
cc INDEX_UNRANK1 unranks an index vector within given upper limits.
c
c  Example:
c
c    N = 3,
c    HI(1) = 4, HI(2) = 2, HI(3) = 3
c    RANK = 11
c
c    A = ( 3, 1, 2 )
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    24 January 2007
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the number of entries in A.
c
c    Input, integer HI(N), the upper limits for the array indices.
c    The lower limit is implicitly 1, and each HI(I) should be at least 1.
c
c    Input, integer RANK, the rank of the desired index vector.
c
c    Output, integer A(N), the index vector of the given rank.
c
      implicit none

      integer n

      integer a(n)
      integer hi(n)
      integer i
      integer j
      integer k
      integer range
      integer rank

      do i = 1, n
        a(i) = 0
      end do
c
c  The rank might be too small.
c
      if ( rank .lt. 1 ) then
        return
      end if

      range = 1
      do i = 1, n
        range = range * hi(i)
      end do
c
c  The rank might be too large.
c
      if ( range .lt. rank ) then
        return
      end if

      k = rank - 1
      do i = n, 1, -1
        range = range / hi(i)
        j = k / range
        a(i) = j + 1
        k = k - j * range
      end do

      return
      end
      subroutine index_unrank2 ( n, lo, hi, rank, a )

c*********************************************************************72
c
cc INDEX_UNRANK2 unranks an index vector within given lower and upper limits.
c
c  Example:
c
c    N = 3,
c    LO(1) = 1, LO(2) = 10, LO(3) = 4
c    HI(1) = 2, HI(2) = 11, HI(3) = 6
c    RANK = 7
c
c    A = ( 1, 11, 5 )
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    24 January 2007
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the number of entries in A.
c
c    Input, integer LO(N), HI(N), the lower and upper limits for the array
c    indices.  It should be the case that LO(I) <= HI(I) for each I.
c
c    Input, integer RANK, the rank of the desired index.
c
c    Output, integer A(N), the index vector of the given rank.
c
      implicit none

      integer n

      integer a(n)
      integer hi(n)
      integer i
      integer j
      integer k
      integer lo(n)
      integer range
      integer rank

      do i = 1, n
        a(i) = 0
      end do
c
c  The rank might be too small.
c
      if ( rank .lt. 1 ) then
        return
      end if

      range = 1
      do i = 1, n
        range = range * ( hi(i) + 1 - lo(i) )
      end do
c
c  The rank might be too large.
c
      if ( range .lt. rank ) then
        return
      end if

      k = rank - 1
      do i = n, 1, -1
        range = range / ( hi(i) + 1 - lo(i) )
        j = k / range
        a(i) = j + lo(i)
        k = k - j * range
      end do

      return
      end
      subroutine ins_perm ( n, ins, p )

c*********************************************************************72
c
cc INS_PERM computes a permutation from its inversion sequence.
c
c  Discussion:
c
c    For a given permutation P acting on objects 1 through N, the
c    inversion sequence INS is defined as:
c
c      INS(1) = 0
c      INS(I) = number of values J < I for which P(I) < P(J).
c
c  Example:
c
c    Input:
c
c      ( 0, 0, 2, 1, 3 )
c
c    Output:
c
c      ( 3, 5, 1, 4, 2 )
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    19 January 2007
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Dennis Stanton, Dennis White,
c    Constructive Combinatorics,
c    Springer, 1986,
c    ISBN: 0387963472.
c
c  Parameters:
c
c    Input, integer N, the number of objects being permuted.
c
c    Input, integer INS(N), the inversion sequence of a permutation.
c    It must be the case that 0 <= INS(I) < I for I = 1 to N.
c
c    Output, integer P(N), the permutation.
c
      implicit none

      integer n

      integer i
      integer ins(n)
      integer itemp
      integer j
      integer p(n)

      call i4vec_indicator ( n, p )

      do i = n, 2, -1

        itemp = p(i-ins(i))

        do j = i-ins(i), i-1
          p(j) = p(j+1)
        end do

        p(i) = itemp

      end do

      return
      end
      subroutine inverse_mod_n ( b, n, y )

c*********************************************************************72
c
cc INVERSE_MOD_N computes the inverse of B mod N.
c
c  Discussion:
c
c    If 
c
c      Y = inverse_mod_n ( B, N )
c
c    then
c
c      mod ( B * Y, N ) = 1
c
c    The value Y will exist if and only if B and N are relatively prime.
c
c  Examples:
c
c    B  N  Y
c
c    1  2  1
c
c    1  3  1
c    2  3  2
c
c    1  4  1
c    2  4  0
c    3  4  3
c
c    1  5  1
c    2  5  3
c    3  5  2
c    4  5  4
c
c    1  6  1
c    2  6  0
c    3  6  0
c    4  6  0
c    5  6  5
c
c    1  7  1
c    2  7  4
c    3  7  5
c    4  7  2
c    5  7  3
c    6  7  6
c
c    1  8  1
c    2  8  0
c    3  8  3
c    4  8  0
c    5  8  5
c    6  8  0
c    7  8  7
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    03 November 2009
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer B, the number whose inverse mod N is desired.
c    B should be positive.  Normally, B < N, but this is not required.
c
c    Input, integer N, the number with respect to which the
c    modulus is computed.  N should be positive.
c
c    Output, integer Y, the inverse of B mod N, or 0 if there
c    is not inverse for B mode N.  1 <= Y < N if the inverse exists.
c
      implicit none

      integer b
      integer b0
      integer n
      integer n0
      integer q
      integer r
      integer t
      integer t0
      integer temp
      integer y

      n0 = n
      b0 = b
      t0 = 0
      t = 1

      q = n / b
      r = n - q * b

10    continue

      if ( 0 .lt. r ) then

        temp = t0 - q * t

        if ( 0 .le. temp ) then
          temp = mod ( temp, n )
        end if

        if ( temp .lt. 0 ) then
          temp = n - mod ( - temp, n )
        end if

        t0 = t
        t = temp
        n0 = b0
        b0 = r
        q = n0 / b0
        r = n0 - q * b0

        go to 10

      end if

      if ( b0 .ne. 1 ) then
        y = 0
        return
      end if

      y = mod ( t, n )

      return
      end
      subroutine involute_enum ( n, s )

c*******************************************************************************
c
cc INVOLUTE_ENUM enumerates the involutions of N objects.
c
c  Discussion:
c
c    An involution is a permutation consisting only of fixed points and
c    pairwise transpositions.
c
c    An involution is its own inverse permutation.
c
c  Recursion:
c
c    S(0) = 1
c    S(1) = 1
c    S(N) = S(N-1) + (N-1) * S(N-2)
c
c  First values:
c
c     N         S(N)
c     0           1
c     1           1
c     2           2
c     3           4
c     4          10
c     5          26
c     6          76
c     7         232
c     8         764
c     9        2620
c    10        9496
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    04 January 2007
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the number of objects to be permuted.
c
c    Output, integer S(0:N), the number of involutions of 0, 1, 2, ... N
c    objects.
c
      implicit none

      integer n

      integer i
      integer s(0:n)

      if ( n .lt. 0 ) then
        return
      end if

      s(0) = 1

      if ( n .le. 0 ) then
        return
      end if

      s(1) = 1

      do i = 2, n
        s(i) = s(i-1) + ( i - 1 ) * s(i-2)
      end do

      return
      end
      subroutine jfrac_to_rfrac ( m, r, s, p, q )

c*********************************************************************72
c
cc JFRAC_TO_RFRAC converts a J-fraction into a rational polynomial fraction.
c
c  Discussion:
c
c    The routine accepts a J-fraction:
c
c        R(1) / ( X + S(1)
c      + R(2) / ( X + S(2)
c      + R(3) / ...
c      + R(M) / ( X + S(M) )... ))
c
c    and returns the equivalent rational polynomial fraction:
c
c      P(1) + P(2) * X + ... + P(M) * X**(M-1)
c      -------------------------------------------------------
c      Q(1) + Q(2) * X + ... + Q(M) * X**(M-1) + Q(M+1) * X**M
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    24 January 2007
c
c  Author:
c
c    Original FORTRAN77 version by John Hart, Ward Cheney, Charles Lawson, 
c    Hans Maehly, Charles Mesztenyi, John Rice, Henry Thatcher, 
c    Christoph Witzgall.
c    This FORTRAN77 version by John Burkardt.
c
c  Reference:
c
c    John Hart, Ward Cheney, Charles Lawson, Hans Maehly, Charles Mesztenyi,
c    John Rice, Henry Thatcher, Christoph Witzgall,
c    Computer Approximations,
c    Wiley, 1968.
c
c  Parameters:
c
c    Input, integer M, defines the number of P, R, and S
c    coefficients, and is one less than the number of Q
c    coefficients.
c
c    Input, double precision R(M), S(M), the coefficients defining 
c    the J-fraction.
c
c    Output, double precision P(M), Q(M+1), the coefficients defining
c    the rational polynomial fraction.  The algorithm used normalizes 
c    the coefficients so that Q(M+1) = 1.0.
c
      implicit none

      integer m
      integer m_max
      parameter ( m_max = 15 )

      double precision a(m_max,m_max)
      double precision b(m_max,m_max)
      integer i
      integer k
      double precision p(m)
      double precision q(m+1)
      double precision r(m)
      double precision s(m)

      if ( m_max .lt. m ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'JFRAC_TO_RFRAC - Fatal error!'
        write ( *, '(a)' ) '  Internal value of M_MAX exceeded by M.'
        stop
      end if

      a(1,1) = r(1)
      b(1,1) = s(1)

      if ( 1 .lt. m ) then

        do k = 2, m
          a(k,k) = r(1)
          b(k,k) = b(k-1,k-1) + s(k)
        end do

        a(1,2) = r(1) * s(2)
        b(1,2) = r(2) + s(1) * s(2)

        do k = 3, m
          a(1,k) = s(k) * a(1,k-1) + r(k) * a(1,k-2)
          a(k-1,k) = a(k-2,k-1) + s(k) * r(1)
          b(1,k) = s(k) * b(1,k-1) + r(k) * b(1,k-2)
          b(k-1,k) = b(k-2,k-1) + s(k) * b(k-1,k-1) + r(k)
        end do

        do k = 4, m
          do i = 2, k-2
            a(i,k) = a(i-1,k-1) + s(k) * a(i,k-1) + r(k) * a(i,k-2)
            b(i,k) = b(i-1,k-1) + s(k) * b(i,k-1) + r(k) * b(i,k-2)
          end do
        end do

      end if

      do i = 1, m
        p(i) = a(i,m)
      end do

      do i = 1, m
        q(i) = b(i,m)
      end do

      q(m+1) = 1.0D+00

      return
      end
      subroutine josephus ( n, m, k, x )

c*******************************************************************************
c
cc JOSEPHUS returns the position X of the K-th man to be executed.
c
c  Discussion:
c
c    The classic Josephus problem concerns a circle of 41 men.
c    Every third man is killed and removed from the circle.  Counting
c    and executing continues until all are dead.  Where was the last
c    survivor sitting?
c
c    Note that the first person killed was sitting in the third position.
c    Moreover, when we get down to 2 people, and we need to count the
c    "third" one, we just do the obvious thing, which is to keep counting
c    around the circle until our count is completed.
c
c    The process may be regarded as generating a permutation of
c    the integers from 1 to N.  The permutation would be the execution
c    list, that is, the list of the executed men, by position number.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    04 January 2007
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    WW Rouse Ball,
c    Mathematical Recreations and Essays,
c    Macmillan, 1962, pages 32-36.
c
c    Donald Knuth,
c    The Art of Computer Programming,
c    Volume 1, Fundamental Algorithms,
c    Addison Wesley, 1968, pages 158-159.
c
c    Donald Knuth,
c    The Art of Computer Programming,
c    Volume 3, Sorting and Searching,
c    Addison Wesley, 1968, pages 18-19.
c
c  Parameters:
c
c    Input, integer N, the number of men.
c    N must be positive.
c
c    Input, integer M, the counting index.
c    M must not be zero.  Ordinarily, M is positive, and no greater than N.
c
c    Input, integer K, the index of the executed man of interest.
c    K must be between 1 and N.
c
c    Output, integer X, the position of the K-th man.
c    X will be between 1 and N.
c
      implicit none

      integer i4_modp
      integer k
      integer m
      integer m2
      integer n
      integer x

      if ( n .le. 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'JOSEPHUS - Fatal error!'
        write ( *, '(a)' ) '  N <= 0.'
        stop
      end if

      if ( m .eq. 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'JOSEPHUS - Fatal error!'
        write ( *, '(a)' ) '  M = 0.'
        stop
      end if

      if ( k .le. 0 .or. n .lt. k ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'JOSEPHUS - Fatal error!'
        write ( *, '(a)' ) '  J <= 0 or N < K.'
        stop
      end if
c
c  In case M is bigger than N, or negative, get the
c  equivalent positive value between 1 and N.
c  You can skip this operation if 1 <= M <= N.
c
      m2 = i4_modp ( m, n )

      x = k * m2

10    continue

      if ( n .lt. x ) then
        x = ( m2 * ( x - n ) - 1 ) / ( m2 - 1 )
        go to 10
      end if

      return
      end
      subroutine ksub_next ( n, k, a, more )

c*******************************************************************************
c
cc KSUB_NEXT generates the subsets of size K from a set of size N.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    04 January 2007
c
c  Author:
c
c    Original FORTRAN77 version by Albert Nijenhuis, Herbert Wilf.
c    This FORTRAN77 version by John Burkardt.
c
c  Reference:
c
c    Albert Nijenhuis, Herbert Wilf,
c    Combinatorial Algorithms for Computers and Calculators,
c    Second Edition,
c    Academic Press, 1978,
c    ISBN: 0-12-519260-6,
c    LC: QA164.N54.
c
c  Parameters:
c
c    Input, integer N, the size of the set from which subsets are drawn.
c
c    Input, integer K, the desired size of the subsets.  K must
c    be between 0 and N.
c
c    Input/output, integer A(K).  A(I) is the I-th element of the
c    subset.  Thus A(I) will be an integer between 1 and N.
c    Note that the routine will return the values in A
c    in sorted order: 1 <= A(1) < A(2) < ... < A(K) <= N
c
c    Input/output, logical MORE.  Set MORE = FALSE before first call
c    for a new sequence of subsets.  It then is set and remains
c    TRUE as long as the subset computed on this call is not the
c    final one.  When the final subset is computed, MORE is set to
c    FALSE as a signal that the computation is done.
c
      implicit none

      integer k

      integer a(k)
      integer j
      integer m
      integer m2
      logical more
      integer n

      save m
      save m2

      data m / 0 /
      data m2 / 0 /

      if ( k .lt. 0 .or. n .lt. k ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'KSUB_NEXT - Fatal error!'
        write ( *, '(a,i8)' ) 'N = ', n
        write ( *, '(a,i8)' ) 'K = ', k
        write ( *, '(a)' ) 'but 0 <= K <= N is required!'
        stop
      end if

      if ( .not. more ) then
        m2 = 0
        m = k
      else
        if ( m2 .lt. n - m ) then
          m = 0
        end if
        m = m + 1
        m2 = a(k+1-m)
      end if

      do j = 1, m
        a(k+j-m) = m2 + j
      end do

      more = ( a(1) .ne. ( n - k + 1 ) )

      return
      end
      subroutine ksub_next2 ( n, k, a, in, iout )

c*********************************************************************72
c
cc KSUB_NEXT2 generates the subsets of size K from a set of size N.
c
c  Discussion:
c
c    This routine uses the revolving door method.  It has no "memory".
c    It simply calculates the successor of the input set,
c    and will start from the beginning after the last set.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    24 January 2007
c
c  Author:
c
c    Original FORTRAN77 version by Albert Nijenhuis, Herbert Wilf.
c    This FORTRAN77 version by John Burkardt.
c
c  Reference:
c
c    Albert Nijenhuis, Herbert Wilf,
c    Combinatorial Algorithms for Computers and Calculators,
c    Second Edition,
c    Academic Press, 1978,
c    ISBN: 0-12-519260-6,
c    LC: QA164.N54.
c
c  Parameters:
c
c    Input, integer N, the size of the set from which subsets are drawn.
c    N must be positive.
c
c    Input, integer K, the size of the desired subset.  K must be
c    between 0 and N.
c
c    Input/output, integer A(K).  On input, the user must
c    supply a subset of size K in A.  That is, A must
c    contain K unique numbers, in order, between 1 and N.  On
c    output, A(I) is the I-th element of the output subset.
c    The output array is also in sorted order.
c
c    Output, integer IN, the element of the output subset which
c    was not in the input set.  Each new subset differs from the
c    last one by adding one element and deleting another.
c
c    Output, integer IOUT, the element of the input subset which
c    is not in the output subset.
c
      implicit none

      integer k

      integer a(k)
      integer in
      integer iout
      integer j
      integer m
      integer n

      if ( n .le. 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'KSUB_NEXT2 - Fatal error!'
        write ( *, '(a,i8)' ) '  N = ', n
        write ( *, '(a)' ) '  but 0 .lt. N is required!'
        stop
      end if

      if ( k .lt. 0 .or. n .lt. k ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'KSUB_NEXT2 - Fatal error!'
        write ( *, '(a,i8)' ) '  N = ', n
        write ( *, '(a,i8)' ) '  K = ', k
        write ( *, '(a)' ) '  but 0 <= K <= N is required!'
        stop
      end if

      j = 0

10    continue

        if ( 0 .lt. j .or. mod ( k, 2 ) .eq. 0 ) then

          j = j + 1

          if ( k .lt. j ) then
            a(k) = k
            in = k
            iout = n
            return
          end if

          if ( a(j) .ne. j ) then

            iout = a(j)
            in = iout - 1
            a(j) = in

            if ( j .ne. 1 ) then
              in = j - 1
              a(j-1) = in
            end if

            return

          end if

        end if

        j = j + 1
        m = n

        if ( j .lt. k ) then
          m = a(j+1) - 1
        end if

        if ( m .ne. a(j) ) then
          go to 20
        end if

      go to 10

20    continue

      in = a(j) + 1
      a(j) = in
      iout = in - 1

      if ( j .ne. 1 ) then
        a(j-1) = iout
        iout = j - 1
      end if

      return
      end
      subroutine ksub_next3 ( n, k, a, more, in, iout )

c*********************************************************************72
c
cc KSUB_NEXT3 generates the subsets of size K from a set of size N.
c
c  Discussion:
c
c    The routine uses the revolving door method.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    24 January 2007
c
c  Author:
c
c    Original FORTRAN77 version by Albert Nijenhuis, Herbert Wilf.
c    This FORTRAN77 version by John Burkardt.
c
c  Reference:
c
c    Albert Nijenhuis, Herbert Wilf,
c    Combinatorial Algorithms for Computers and Calculators,
c    Second Edition,
c    Academic Press, 1978,
c    ISBN: 0-12-519260-6,
c    LC: QA164.N54.
c
c  Parameters:
c
c    Input, integer N, the size of the set from which subsets are drawn.
c    N must be positive.
c
c    Input, integer K, the size of the desired subsets.  K must be
c    between 0 and N.
c
c    Input/output, integer A(K).  A(I) is the I-th element of the
c    output subset.  The elements of A are sorted.
c
c    Input/output, logical MORE.  On first call, set MORE = FALSE
c    to signal the beginning.  MORE will be set to TRUE, and on
c    each call, the routine will return another K-subset.
c    Finally, when the last subset has been returned,
c    MORE will be set FALSE and you may stop calling.
c
c    Output, integer IN, the element of the output subset which
c    was not in the input set.  Each new subset differs from the
c    last one by adding one element and deleting another.  IN is not
c    defined the first time that the routine returns, and is
c    set to zero.
c
c    Output, integer IOUT, the element of the input subset which is
c    not in the output subset.  IOUT is not defined the first time
c    the routine returns, and is set to zero.
c
      implicit none

      integer k

      integer a(k)
      integer in
      integer iout
      integer j
      integer m
      logical more
      integer n

      if ( n .le. 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'KSUB_NEXT3 - Fatal error!'
        write ( *, '(a,i8)' ) '  N = ', n
        write ( *, '(a)' ) '  but 0 .lt. N is required!'
        stop
      end if

      if ( k .lt. 0 .or. n .lt. k ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'KSUB_NEXT3 - Fatal error!'
        write ( *, '(a,i8)' ) '  N = ', n
        write ( *, '(a,i8)' ) '  K = ', k
        write ( *, '(a)' ) '  but 0 <= K <= N is required!'
        stop
      end if

      if ( .not. more ) then
        in = 0
        iout = 0
        call i4vec_indicator ( k, a )
        more = ( k .ne. n )
        return
      end if

      j = 0

10    continue

        if ( 0 .lt. j .or. mod ( k, 2 ) .eq. 0 ) then

          j = j + 1

          if ( a(j) .ne. j ) then

            iout = a(j)
            in = iout - 1
            a(j) = in

            if ( j .ne. 1 ) then
              in = j - 1
              a(j-1) = in
            end if

            if ( k .ne. 1 ) then
              more = ( a(k-1) .eq. k - 1 )
            end if

            more = ( .not. more ) .or. ( a(k) .ne. n )

            return

          end if

        end if

        j = j + 1
        m = n

        if ( j .lt. k ) then
          m = a(j+1) - 1
        end if

        if ( m .ne. a(j) ) then
          go to 20
        end if

      go to 10

20    continue

      in = a(j) + 1
      a(j) = in
      iout = in - 1

      if ( j .ne. 1 ) then
        a(j-1) = iout
        iout = j - 1
      end if

      if ( k .ne. 1 ) then
        more = ( a(k-1) .eq. k-1 )
      end if

      more = ( .not. more ) .or. ( a(k) .ne. n )

      return
      end
      subroutine ksub_next4 ( n, k, a, done )

c*********************************************************************72
c
cc KSUB_NEXT4 generates the subsets of size K from a set of size N.
c
c  Discussion:
c
c    The subsets are generated one at a time.
c
c    The routine should be used by setting DONE to TRUE, and then calling
c    repeatedly.  Each call returns with DONE equal to FALSE, the array
c    A contains information defining a new subset.  When DONE returns
c    equal to TRUE, there are no more subsets.
c
c    There are ( N*(N-1)*...*(N+K-1)) / ( K*(K-1)*...*2*1) such subsets.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    24 January 2007
c
c  Parameters:
c
c    Input, integer N, the size of the entire set.
c
c    Input, integer K, the size of the desired subset.  K must be
c    between 0 and N.
c
c    Input/output, integer A(K), contains information about
c    the subsets.  On the first call with DONE = TRUE, the input contents
c    of A don't matter.  Thereafter, the input value of A
c    should be the same as the output value of the previous call.
c    In other words, leave the array alonec
c    On output, as long as DONE is returned FALSE, A contains
c    information defining a subset of K elements of a set of N elements.
c    In other words, A will contain K distinct numbers (in order)
c    between 1 and N.
c
c    Input/output, logical DONE.
c    On the first call, DONE is an input quantity with a value
c    of TRUE which tells the program to initialize data and
c    return the first subset.
c    On return, DONE is an output quantity that is TRUE as long as
c    the routine is returning another subset, and FALSE when
c    there are no more.
c
      implicit none

      integer k

      integer a(k)
      logical done
      integer j
      integer jsave
      integer n

      if ( k .lt. 0 .or. n .lt. k ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'KSUB_NEXT4 - Fatal error!'
        write ( *, '(a,i8)' ) '  N = ', n
        write ( *, '(a,i8)' ) '  K = ', k
        write ( *, '(a)' ) '  but 0 <= K <= N is required!'
        stop
      end if
c
c  First call:
c
      if ( done ) then

        call i4vec_indicator ( k, a )

        if ( 0 .lt. n ) then
          done = .false.
        else
          done = .true.
        end if
c
c  Next call.
c
      else

        if ( a(1) .lt. n-k+1 ) then

          done = .false.

          jsave = k

          do j = 1, k-1

            if ( a(j) + 1 .lt. a(j+1) ) then
              jsave = j
              go to 10
            end if

          end do

10        continue

          call i4vec_indicator ( jsave-1, a )
          a(jsave) = a(jsave) + 1

        else

          done = .true.

        end if

      end if

      return
      end
      subroutine ksub_random ( n, k, seed, a )

c*********************************************************************72
c
cc KSUB_RANDOM selects a random subset of size K from a set of size N.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    18 January 2007
c
c  Author:
c
c    Original FORTRAN77 version by Albert Nijenhuis, Herbert Wilf.
c    This FORTRAN77 version by John Burkardt.
c
c  Reference:
c
c    Albert Nijenhuis, Herbert Wilf,
c    Combinatorial Algorithms for Computers and Calculators,
c    Second Edition,
c    Academic Press, 1978,
c    ISBN: 0-12-519260-6,
c    LC: QA164.N54.
c
c  Parameters:
c
c    Input, integer N, the size of the set from which subsets are drawn.
c
c    Input, integer K, number of elements in desired subsets.  K must
c    be between 0 and N.
c
c    Input/output, integer SEED, a seed for the random number generator.
c
c    Output, integer A(K).  A(I) is the I-th element of the
c    output set.  The elements of A are in order.
c
      implicit none

      integer k

      integer a(k)
      integer i
      integer i4_uniform
      integer ids
      integer ihi
      integer ip
      integer ir
      integer is
      integer ix
      integer l
      integer ll
      integer m
      integer m0
      integer n
      integer seed

      if ( k .lt. 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'KSUB_RANDOM - Fatal error!'
        write ( *, '(a,i8)' ) '  K = ', k
        write ( *, '(a)' ) '  but 0 <= K is required!'
        stop
      else if ( n .lt. k ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'KSUB_RANDOM - Fatal error!'
        write ( *, '(a,i8)' ) '  N = ', n
        write ( *, '(a,i8)' ) '  K = ', k
        write ( *, '(a)' ) '  K <= N is required!'
        stop
      end if

      if ( k .eq. 0 ) then
        return
      end if

      do i = 1, k
        a(i) = ( ( i - 1 ) * n ) / k
      end do

      do i = 1, k

10      continue

          ix = i4_uniform ( 1, n, seed )

          l = 1 + ( ix * k - 1 ) / n

          if ( a(l) .lt. ix ) then
            go to 20
          end if

        go to 10

20      continue

        a(l) = a(l) + 1

      end do

      ip = 0
      is = k

      do i = 1, k

        m = a(i)
        a(i) = 0

        if ( m .ne. ( ( i - 1 ) * n ) / k ) then
          ip = ip + 1
          a(ip) = m
        end if

      end do

      ihi = ip

      do i = 1, ihi
        ip = ihi + 1 - i
        l = 1 + ( a(ip) * k - 1 ) / n
        ids = a(ip) - ( ( l - 1 ) * n ) / k
        a(ip) = 0
        a(is) = l
        is = is - ids
      end do

      do ll = 1, k

        l = k + 1 - ll

        if ( a(l) .ne. 0 ) then
          ir = l
          m0 = 1 + ( ( a(l) - 1 ) * n ) / k
          m = ( a(l) * n ) / k - m0 + 1
        end if

        ix = i4_uniform ( m0, m0 + m - 1, seed )

        i = l + 1

30      continue

        if ( i .le. ir ) then

          if ( ix .lt. a(i) ) then
            go to 40
          end if

          ix = ix + 1
          a(i-1) = a(i)
          i = i + 1

          go to 30

        end if

40      continue

        a(i-1) = ix
        m = m - 1

      end do

      return
      end
      subroutine ksub_random2 ( n, k, seed, a )

c*********************************************************************72
c
cc KSUB_RANDOM2 selects a random subset of size K from a set of size N.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    24 January 2007
c
c  Author:
c
c    Original FORTRAN77 version by Albert Nijenhuis, Herbert Wilf.
c    This FORTRAN77 version by John Burkardt.
c
c  Reference:
c
c    Albert Nijenhuis, Herbert Wilf,
c    Combinatorial Algorithms for Computers and Calculators,
c    Second Edition,
c    Academic Press, 1978,
c    ISBN: 0-12-519260-6,
c    LC: QA164.N54.
c
c  Parameters:
c
c    Input, integer N, the size of the set.
c
c    Input, integer K, the size of the subset, between 0 and N.
c
c    Input/output, integer SEED, a seed for the random number generator.
c
c    Output, integer A(K), the indices of the selected elements.
c
      implicit none

      integer k

      integer a(k)
      integer available
      integer candidate
      integer have
      integer n
      integer need
      double precision r
      double precision r8_uniform_01
      integer seed

      if ( k .lt. 0 .or. n .lt. k ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'KSUB_RANDOM2 - Fatal error!'
        write ( *, '(a,i8)' ) '  N = ', n
        write ( *, '(a,i8)' ) '  K = ', k
        write ( *, '(a)' ) '  but 0 <= K <= N is required!'
        stop
      end if

      if ( k .eq. 0 ) then
        return
      end if

      need = k
      have = 0

      available = n
      candidate = 0

10    continue

        candidate = candidate + 1

        r = r8_uniform_01 ( seed )

        if ( dble ( available ) * r .le. dble ( need ) ) then

          need = need - 1
          have = have + 1
          a(have) = candidate

          if ( need .le. 0 ) then
            go to 20
          end if

        end if

        available = available - 1

      go to 10

20    continue

      return
      end
      subroutine ksub_random3 ( n, k, seed, a )

c*********************************************************************72
c
cc KSUB_RANDOM3 selects a random subset of size K from a set of size N.
c
c  Discussion:
c
c    This routine uses Floyd's algorithm.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    24 January 2007
c
c  Parameters:
c
c    Input, integer N, the size of the set from which subsets are drawn.
c
c    Input, integer K, number of elements in desired subsets.  K must
c    be between 0 and N.
c
c    Input/output, integer SEED, a seed for the random number generator.
c
c    Output, integer A(N).  I is an element of the subset
c    if A(I) = 1, and I is not an element if A(I)=0.
c
      implicit none

      integer n

      integer a(n)
      integer i
      integer i4_uniform
      integer j
      integer k
      integer seed

      if ( k .lt. 0 .or. n .lt. k ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'KSUB_RANDOM3 - Fatal error!'
        write ( *, '(a,i8)' ) '  N = ', n
        write ( *, '(a,i8)' ) '  K = ', k
        write ( *, '(a)' ) '  but 0 <= K <= N is required!'
        stop
      end if

      do i = 1, n
        a(i) = 0
      end do

      if ( k .eq. 0 ) then
        return
      end if

      do i = n-k+1, n

        j = i4_uniform ( 1, i, seed )

        if ( a(j) .eq. 0 ) then
          a(j) = 1
        else
          a(i) = 1
        end if

      end do

      return
      end
      subroutine ksub_random4 ( n, k, seed, a )

c*********************************************************************72
c
cc KSUB_RANDOM4 selects a random subset of size K from a set of size N.
c
c  Discussion:
c
c    This routine is somewhat impractical for the given problem, but
c    it is included for comparison, because it is an interesting
c    approach that is superior for certain applications.
c
c    The approach is mainly interesting because it is "incremental";
c    it proceeds by considering every element of the set, and does not
c    need to know how many elements there are.
c
c    This makes this approach ideal for certain cases, such as the
c    need to pick 5 lines at random from a text file of unknown length,
c    or to choose 6 people who call a certain telephone number on a
c    given day.  Using this technique, it is possible to make the
c    selection so that, whenever the input stops, a valid uniformly
c    random subset has been chosen.
c
c    Obviously, if the number of items is known in advance, and 
c    it is easy to extract K items directly, there is no need for
c    this approach, and it is less efficient since, among other costs,
c    it has to generate a random number for each item, and make an
c    acceptance/rejection test.
c
c    This routine is based on "8.6: Picking a Random Line from a File",
c    in the Perl Cookbook.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    24 January 2007
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Tom Christiansen, Nathan Torkington,
c    Perl Cookbook,
c    OReilly, 1999.
c
c  Parameters:
c
c    Input, integer N, the size of the set from which subsets are drawn.
c
c    Input, integer K, number of elements in desired subsets.  K must
c    be between 0 and N.
c
c    Input/output, integer SEED, a seed for the random number generator.
c
c    Output, integer A(K), contains the indices of the selected items.
c
      implicit none

      integer n

      integer a(n)
      integer i
      integer i4_uniform
      integer k
      integer next
      double precision r
      double precision r8_uniform_01
      integer seed

      next = 0
c
c  Here, we use a DO WHILE to suggest that the algorithm
c  proceeds to the next item, without knowing how many items
c  there are in total.  
c
c  Note that this is really the only place where N occurs,
c  so other termination criteria could be used, and we really
c  don't need to know the value of Nc
c
10    continue

      if ( next .lt. n ) then

        next = next + 1 

        if ( next .le. k ) then

          i = next
          a(i) = next

        else

          r = r8_uniform_01 ( seed )

          if ( r * dble ( next ) .le. dble ( k ) ) then
            i = i4_uniform ( 1, k, seed )
            a(i) = next
          end if

        end if

        go to 10

      end if

      return
      end
      subroutine ksub_random5 ( n, k, seed, a )

c*********************************************************************72
c
cc KSUB_RANDOM5 selects a random subset of size K from a set of size N.
c
c  Discussion:
c
c    Consider the set A(1:N) = 1, 2, 3, ... N.  
c    Choose a random index I1 between 1 and N, and swap items A(1) and A(I1).
c    Choose a random index I2 between 2 and N, and swap items A(2) and A(I2).
c    repeat K times.
c    A(1:K) is your random K-subset.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    09 June 2011
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the size of the set from which subsets
c    are drawn.
c
c    Input, integer K, number of elements in desired subsets.
c    1 <= K <= N.
c
c    Input/output, integer SEED, a seed for the random
c    number generator.
c
c    Output, integer A(K), the indices of the randomly
c    chosen elements.
c
      implicit none

      integer k
      integer n

      integer a(k)
      integer b(n)
      integer i
      integer i4_uniform
      integer j
      integer seed
      integer t
c
c  Let B index the set.
c
      do i = 1, n
        b(i) = i
      end do
c
c  Choose item 1 from N things,
c  choose item 2 from N-1 things,
c  choose item K from N-K+1 things.
c
      do i = 1, k

        j = i4_uniform ( i, n, seed )

        t    = b(i)
        b(i) = b(j)
        b(j) = t

      end do
c
c  Copy the first K elements.
c
      do i = 1, k
        a(i) = b(i)
      end do
c
c  Put the elements in ascending order.
c
      call i4vec_sort_heap_a ( k, a )

      return
      end
      subroutine ksub_rank ( k, a, rank )

c*********************************************************************72
c
cc KSUB_RANK computes the rank of a K subset of an N set.
c
c  Discussion:
c
c    The routine accepts an array representing a subset of size K from a set
c    of size N, and returns the rank (or order) of that subset. 
c
c    It uses the same ranking that KSUB_NEXT2 uses to generate all the subsets 
c    one at a time.  
c
c    Note the value of N is not input, and is not, in fact,
c    needed.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    24 January 2007
c
c  Author:
c
c    Original FORTRAN77 version by Albert Nijenhuis, Herbert Wilf.
c    This FORTRAN77 version by John Burkardt.
c
c  Reference:
c
c    Albert Nijenhuis, Herbert Wilf,
c    Combinatorial Algorithms for Computers and Calculators,
c    Second Edition,
c    Academic Press, 1978,
c    ISBN: 0-12-519260-6,
c    LC: QA164.N54.
c
c  Parameters:
c
c    Input, integer K, the number of elements in the subset.
c
c    Input, integer A(K), contains K distinct numbers between
c    1 and N, in order.
c
c    Output, integer RANK, the rank of this subset.
c
      implicit none

      integer k

      integer a(k)
      integer i
      integer iprod
      integer j
      integer rank

      rank = 0

      do i = 1, k

        iprod = 1

        do j = i+1, a(i)-1
          iprod = iprod * j
        end do

        do j = 1, a(i)-i-1
          iprod = iprod / j
        end do

        if ( a(i) .eq. 1 ) then
          iprod = 0
        end if

        rank = rank + iprod

      end do

      rank = rank + 1

      return
      end
      subroutine ksub_unrank ( k, rank, a )

c*********************************************************************72
c
cc KSUB_UNRANK returns the subset of a given rank.
c
c  Discussion:
c
c    The routine is given a rank and returns the corresponding subset of K
c    elements of a set of N elements.  
c
c    It uses the same ranking that KSUB_NEXT2 uses to generate all the subsets 
c    one at a time.  
c
c    Note that the value of N itself is not input, nor is it needed.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    24 January 2007
c
c  Author:
c
c    Original FORTRAN77 version by Albert Nijenhuis, Herbert Wilf.
c    This FORTRAN77 version by John Burkardt.
c
c  Reference:
c
c    Albert Nijenhuis, Herbert Wilf,
c    Combinatorial Algorithms for Computers and Calculators,
c    Second Edition,
c    Academic Press, 1978,
c    ISBN: 0-12-519260-6,
c    LC: QA164.N54.
c
c  Parameters:
c
c    Input, integer K, the number of elements in the subset.
c
c    Input, integer RANK, the rank of the desired subset.
c    There are ( N*(N-1)*...*(N+K-1)) / ( K*(K-1)*...*2*1) such
c    subsets, so RANK must be between 1 and that value.
c
c    Output, integer A(K), K distinct integers in order between
c    1 and N, which define the subset.
c
      implicit none

      integer k

      integer a(k)
      integer i
      integer ip
      integer iprod
      integer jrank
      integer rank

      jrank = rank - 1

      do i = k, 1, -1

        ip = i - 1
        iprod = 1

10      continue

          ip = ip + 1

          if ( ip .ne. i ) then
            iprod = ( ip * iprod ) / ( ip - i )
          end if

          if ( jrank .lt. iprod ) then
            go to 20
          end if

        go to 10

20      continue

        if ( ip .ne. i ) then
          iprod = ( ( ip - i ) * iprod ) / ip
        end if

        jrank = jrank - iprod
        a(i) = ip

      end do

      return
      end
      subroutine lvec_next ( n, lvec )

c*********************************************************************72
c
cc  Purpose:
c
c    LVEC_NEXT generates the next logical vector.
c
c  Discussion:
c
c    In the following discussion, we will let '0' stand for FALSE and
c    '1' for TRUE.
c
c    The vectors have the order
c
c      (0,0,...,0),
c      (0,0,...,1),
c      ...
c      (1,1,...,1)
c
c    and the "next" vector after (1,1,...,1) is (0,0,...,0).  That is,
c    we allow wrap around.
c
c  Example:
c
c    N = 3
c
c    Input      Output
c    -----      ------
c    0 0 0  =>  0 0 1
c    0 0 1  =>  0 1 0
c    0 1 0  =>  0 1 1
c    0 1 1  =>  1 0 0
c    1 0 0  =>  1 0 1
c    1 0 1  =>  1 1 0
c    1 1 0  =>  1 1 1
c    1 1 1  =>  0 0 0
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    31 May 2008
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the dimension of the vectors.
c
c    Input/output, logical LVEC(N), on output, the successor to the
c    input vector.
c
      implicit none

      integer n

      integer i
      logical lvec(n)

      do i = n, 1, -1

        if ( .not. lvec(i) ) then
          lvec(i) = .true.
          return
        end if

        lvec(i) = .false.

      end do

      return
      end
      subroutine matrix_product_opt ( n, rank, cost, order )

c*********************************************************************72
c
cc MATRIX_PRODUCT_OPT determines the optimal cost of a matrix product.
c
c  Discussion:
c
c    The cost of multiplying an LxM matrix by an M by N matrix is
c    assessed as L*M*N.
c
c    Any particular order of multiplying a set of N matrices is equivalent
c    to parenthesizing an expression of N objects.
c
c    The actual number of ways of parenthesizing an expression
c    of N objects is C(N), the N-th Catalan number.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    24 January 2007
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Robert Sedgewick,
c    Algorithms in C,
c    Addison-Wesley, 1990,
c    ISBN: 0-201-51425-7,
c    LC: QA76.73.C15S43.
c
c  Parameters:
c
c    Input, integer N, the number of matrices to be multiplied.
c
c    Input, integer RANK(N+1), the rank information for the matrices.
c    Matrix I has RANK(I) rows and RANK(I+1) columns.
c
c    Output, integer COST, the cost of the multiplication if the optimal
c    order is used.
c
c    Output, integer ORDER(N-1), indicates the order in which the N-1
c    multiplications are to be carried out.  ORDER(1) is the first
c    multiplication to do, and so on.
c
      implicit none

      integer stack_max
      parameter ( stack_max = 100 )
      integer n

      integer best(n,n)
      integer cost
      integer cost2(n,n)
      integer cost3
      integer i
      integer i1
      integer i2
      integer i3
      integer i4_huge
      integer j
      integer k
      integer order(n-1)
      integer rank(n+1)
      integer stack(stack_max)
      integer stack_num
      integer step
c
c  Initialize the cost matrix.
c
      do i = 1, n
        do j = 1, i
          cost2(i,j) = 0
        end do
        do j = i+1, n
          cost2(i,j) = i4_huge ( )
        end do
      end do
c
c  Initialize the BEST matrix.
c
      do j = 1, n
        do i = 1, n
          best(i,j) = 0
        end do
      end do
c
c  Compute the cost and best matrices.
c
      do j = 1, n-1
        do i = 1, n-j
          do k = i+1, i+j
            cost3 = cost2(i,k-1) + cost2(k,i+j) 
     &        + rank(i) * rank(k) * rank(i+j+1)
            if ( cost3 .lt. cost2(i,i+j) ) then
              cost2(i,i+j) = cost3
              best(i,i+j) = k
            end if
          end do
        end do
      end do
c
c  Pick off the optimal cost.
c
      cost = cost2(1,n)
c
c  Backtrack to determine the optimal order.
c
      stack_num = 0

      i1 = 1
      i2 = n

      if ( i1 + 1 .lt. i2 ) then
        stack_num = stack_num + 1
        stack(stack_num) = i1
        stack_num = stack_num + 1
        stack(stack_num) = i2
      end if

      step = n - 1
c
c  Take an item off the stack.
c
10    continue

      if ( 0 .lt. stack_num ) then

        i3 = stack(stack_num)
        stack_num = stack_num - 1
        i1 = stack(stack_num)
        stack_num = stack_num - 1

        i2 = best(i1,i3)

        order(step) = i2 - 1
        step = step - 1
c
c  The left chunk is matrices (I1...I2-1)
c
        if ( i1 .eq. i2 - 1 ) then

        else if ( i1 + 1 .eq. i2 - 1 ) then
          order(step) = i2 - 2
          step = step - 1
        else
          stack_num = stack_num + 1
          stack(stack_num) = i1
          stack_num = stack_num + 1
          stack(stack_num) = i2 - 1
        end if
c
c  The right chunk is matrices (I2...I3)
c
        if ( i2 .eq. i3 ) then

        else if ( i2 + 1 .eq. i3 ) then
          order(step) = i2
          step = step - 1
        else
          stack_num = stack_num + 1
          stack(stack_num) = i2
          stack_num = stack_num + 1
          stack(stack_num) = i3
        end if

        go to 10

      end if

      return
      end
      subroutine moebius_matrix ( n, a, mu )

c*********************************************************************72
c
cc MOEBIUS_MATRIX finds the Moebius matrix from a covering relation.
c
c  Discussion:
c
c    This routine can be called with A and MU being the same matrix.
c    The routine will correctly compute the Moebius matrix, which
c    will, in this case, overwrite the input matrix.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    01 July 2004
c
c  Author:
c
c    Original FORTRAN77 version by Albert Nijenhuis, Herbert Wilf.
c    This FORTRAN77 version by John Burkardt.
c
c  Reference:
c
c    Albert Nijenhuis, Herbert Wilf,
c    Combinatorial Algorithms for Computers and Calculators,
c    Second Edition,
c    Academic Press, 1978,
c    ISBN: 0-12-519260-6,
c    LC: QA164.N54.
c
c  Parameters:
c
c    Input, integer N, number of elements in the partially ordered set.
c
c    Input, integer A(N,N).  A(I,J) = 1 if I is covered by J,
c    0 otherwise.
c
c    Output, integer MU(N,N), the Moebius matrix as computed by the routine.
c
      implicit none

      integer n

      integer a(n,n)
      integer i
      integer j
      integer mu(n,n)
      integer p(n)
      integer q(n)
c
c  Compute a reordering P of the elements of the partially ordered matrix.
c
      call triang ( n, a, p )
c
c  Copy the matrix.
c
        do j = 1, n
          do i = 1, n
            mu(i,j) = a(i,j)
          end do
        end do
c
c  Apply the reordering to MU.
c
      call i4mat_perm2 ( n, n, mu, p, p )
c
c  Negate the (strict) upper triangular elements of MU.
c
      do i = 1, n-1
        do j = i+1, n
          mu(i,j) = -mu(i,j)
        end do
      end do
c
c  Compute the inverse of MU.
c
      call i4mat_u1_inverse ( n, mu, mu )
c
c  All nonzero elements are reset to 1.
c
      do i = 1, n
        do j = i, n
          if ( mu(i,j) .ne. 0 ) then
            mu(i,j) = 1
          end if
        end do
      end do
c
c  Invert the matrix again.
c
      call i4mat_u1_inverse ( n, mu, mu )
c
c  Compute the inverse permutation.
c
      do i = 1, n
        q(p(i)) = i
      end do
c
c  Unpermute the rows and columns of MU.
c
      call i4mat_perm2 ( n, n, mu, q, q )

      return
      end
      subroutine morse_thue ( i, s )

c*********************************************************************72
c
cc MORSE_THUE generates a Morse_Thue number.
c
c  Discussion:
c
c    The Morse_Thue sequence can be defined in a number of ways.
c
c    A) Start with the string containing the single letter '0'; then
c       repeatedly apply the replacement rules '0' -> '01' and
c       '1' -> '10' to the letters of the string.  The Morse_Thue sequence
c       is the resulting letter sequence.
c
c    B) Starting with the string containing the single letter '0',
c       repeatedly append the binary complement of the string to itself.
c       Thus, '0' becomes '0' + '1' or '01', then '01' becomes
c       '01' + '10', which becomes '0110' + '1001', and so on.
c
c    C) Starting with I = 0, the I-th Morse-Thue number is determined
c       by taking the binary representation of I, adding the digits,
c       and computing the remainder modulo 2.
c
c  Example:
c
c     I  binary   S
c    --  ------  --
c     0       0   0
c     1       1   1
c     2      10   1
c     3      11   0
c     4     100   1
c     5     101   0
c     6     110   0
c     7     111   1
c     8    1000   1
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
c    Input, integer I, the index of the Morse-Thue number.
c    Normally, I is 0 or greater, but any value is allowed.
c
c    Output, integer S, the Morse-Thue number of index I.
c
      implicit none

      integer nbits
      parameter ( nbits = 32 )

      integer b(nbits)
      integer i
      integer i_copy
      integer i4vec_sum
      integer s

      i_copy = abs ( i )
c
c  Expand I into binary form.
c
      call ui4_to_ubvec ( i_copy, nbits, b )
c
c  Sum the 1's in the binary representation.
c
      s = i4vec_sum ( nbits, b )
c
c  Take the value modulo 2.
c
      s = mod ( s, 2 )

      return
      end
      subroutine multinomial_coef1 ( nfactor, factor, ncomb )

c*********************************************************************72
c
cc MULTINOMIAL_COEF1 computes a multinomial coefficient.
c
c  Discussion:
c
c    The multinomial coefficient is a generalization of the binomial
c    coefficient.  It may be interpreted as the number of combinations of
c    N objects, where FACTOR(1) objects are indistinguishable of type 1,
c    ... and FACTOR(NFACTOR) are indistinguishable of type NFACTOR,
c    and N is the sum of FACTOR(1) through FACTOR(NFACTOR).
c
c    NCOMB = N! / ( FACTOR(1)! FACTOR(2)! ... FACTOR(NFACTOR)! )
c
c    The logarithm of the Gamma function is used, to avoid overflow.
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
c    Input, integer NFACTOR, the number of factors.
c
c    Input, integer FACTOR(NFACTOR), contains the factors.
c    0 <= FACTOR(I)
c
c    Output, integer NCOMB, the value of the multinomial coefficient.
c
      implicit none

      integer nfactor

      double precision arg
      double precision fack
      double precision facn
      integer factor(nfactor)
      double precision gamma_log
      integer i
      integer i4vec_sum
      integer n
      integer ncomb
c
c  Each factor must be nonnegative.
c
      do i = 1, nfactor

        if ( factor(i) .lt. 0 ) then
          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) 'MULTINOMIAL_COEF1 - Fatal error!'
          write ( *, '(a,i8,a,i8)' ) '  Factor ', i, ' = ', factor(i)
          write ( *, '(a)' ) '  But this value must be nonnegative.'
          stop
        end if

      end do
c
c  The factors sum to N.
c
      n = i4vec_sum ( nfactor, factor )

      arg = dble ( n + 1 )
      facn = gamma_log ( arg )

      do i = 1, nfactor

        arg = dble ( factor(i) + 1 )
        fack = gamma_log ( arg )
        facn = facn - fack

      end do

      ncomb = nint ( exp ( facn ) )

      return
      end
      subroutine multinomial_coef2 ( nfactor, factor, ncomb )

c*********************************************************************72
c
cc MULTINOMIAL_COEF2 computes a multinomial coefficient.
c
c  Discussion:
c
c    The multinomial coefficient is a generalization of the binomial
c    coefficient.  It may be interpreted as the number of combinations of
c    N objects, where FACTOR(1) objects are indistinguishable of type 1,
c    ... and FACTOR(NFACTOR) are indistinguishable of type NFACTOR,
c    and N is the sum of FACTOR(1) through FACTOR(NFACTOR).
c
c    NCOMB = N! / ( FACTOR(1)! FACTOR(2)! ... FACTOR(NFACTOR)! )
c
c    A direct method is used, which should be exact.  However, there
c    is a possibility of intermediate overflow of the result.
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
c    Input, integer NFACTOR, the number of factors.
c
c    Input, integer FACTOR(NFACTOR), contains the factors.
c    0 <= FACTOR(I)
c
c    Output, integer NCOMB, the value of the multinomial coefficient.
c
      implicit none

      integer nfactor

      integer factor(nfactor)
      integer i
      integer j
      integer k
      integer ncomb
c
c  Each factor must be nonnegative.
c
      do i = 1, nfactor

        if ( factor(i) .lt. 0 ) then
          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) 'MULTINOMIAL_COEF2 - Fatal error!'
          write ( *, '(a,i8,a,i8)' ) '  Factor ', i, ' = ', factor(i)
          write ( *, '(a)' ) '  But this value must be nonnegative.'
          stop
        end if

      end do

      ncomb = 1
      k = 0

      do i = 1, nfactor

        do j = 1, factor(i)
          k = k + 1
          ncomb = ( ncomb * k ) / j
        end do

      end do

      return
      end
      subroutine network_flow_max ( nnode, nedge, iendpt, icpflo, 
     &  source, sink, cut, node_flow )

c*********************************************************************72
c
cc NETWORK_FLOW_MAX finds the maximal flow and a minimal cut in a network.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    26 January 2007
c
c  Author:
c
c    Original FORTRAN77 version by Albert Nijenhuis, Herbert Wilf.
c    This FORTRAN77 version by John Burkardt.
c
c  Reference:
c
c    Albert Nijenhuis, Herbert Wilf,
c    Combinatorial Algorithms for Computers and Calculators,
c    Second Edition,
c    Academic Press, 1978,
c    ISBN: 0-12-519260-6,
c    LC: QA164.N54.
c
c  Parameters:
c
c    Input, integer NNODE, the number of nodes.
c
c    Input, integer NEDGE, the number of edges.
c
c    Input/output, integer IENDPT(2,NEDGE), the edges of the network,
c    defined as pairs of nodes.  Each edge should be listed TWICE,
c    the second time in reverse order.  On output, the edges have
c    been reordered, and so the columns of IENDPT have been rearranged.
c
c    Input/output, integer ICPFLO(2,NEDGE).
c    On input, ICPFLO(1,I) is the capacity of edge I.  On output,
c    ICPFLO(2,I) is the flow on edge I and ICPFLO(1,I) has
c    been rearranged to match the reordering of IENDPT.
c
c    Input, integer SOURCE, the designated source node.
c
c    Input, integer SINK, the designated sink node.
c
c    Output, integer CUT(NNODE).  CUT(I) = 1 if node I is in the
c    minimal cut set, otherwise 0.
c
c    Output, integer NODE_FLOW(NNODE), the flow through each node.
c
      implicit none

      integer nedge
      integer nnode

      integer cut(nnode)
      integer del
      integer i
      integer iarray(nnode)
      integer icpflo(2,nedge)
      integer ien1
      integer ien2
      integer iendpt(2,nedge)
      integer indx
      integer ip
      integer iparm
      integer iq
      integer ir
      integer iread
      integer irite
      integer is
      integer isort
      integer it
      integer j
      integer kz
      integer lst
      integer m
      integer node_flow(nnode)
      integer sink
      integer source
      integer work1(nnode)
      integer work2(nnode)
c
c  Initialization.
c
      do i = 1, nnode
        iarray(i) = 0
      end do

      ien1 = 0
      ien2 = 0
      del = 0

      do j = 1, nedge
        icpflo(2,j) = 0
      end do

      do i = 1, nedge

        ip = iendpt(1,i)

        if ( ip .eq. source ) then
          del = del + icpflo(1,i)
        end if

        iarray(ip) = iarray(ip) + 1

      end do

      node_flow(source) = del
      is = 1

      do i = 1, nnode
        it = iarray(i)
        iarray(i) = is
        work1(i) = is
        is = is + it
      end do

      isort = 0
c
c  Sorting.
c
10    continue

      indx = 0

50    continue

        call sort_heap_external ( nedge, indx, ien1, ien2, is )

        if ( indx .lt. 0 ) then

          is = iendpt(1,ien1) - iendpt(1,ien2)

          if ( is .eq. 0 ) then
            is = iendpt(2,ien1) - iendpt(2,ien2)
          end if

        else if ( 0 .lt. indx ) then

          do ir = 1, 2
            call i4_swap ( iendpt(ir,ien1), iendpt(ir,ien2) )
            call i4_swap ( icpflo(ir,ien1), icpflo(ir,ien2) )
          end do

        else

          if ( 0 .lt. isort ) then
            return
          end if

          do i = 1, nedge
            iq = iendpt(2,i)
            iendpt(1,i) = work1(iq)
            work1(iq) = work1(iq) + 1
          end do

          go to 100

        end if

      go to 50

80    continue

      iendpt(1,iendpt(1,ien1)) = ien2
      iendpt(1,iendpt(1,ien2)) = ien1

      do ir = 1, 2
        call i4_swap ( iendpt(ir,ien1), iendpt(ir,ien2) )
        call i4_swap ( icpflo(ir,ien1), icpflo(ir,ien2) )
      end do

      if ( indx .lt. 0 ) then
        work2(iq) = ien2
        go to 280
      end if

      if ( indx .eq. 0 ) then
        go to 170
      end if

      go to 50

100   continue

      indx = 0

      do i = 1, nnode

        if ( i .ne. source ) then
          node_flow(i) = 0
        end if

        work2(i) = nedge + 1

        if ( i .lt. nnode ) then
          work2(i) = iarray(i+1)
        end if

        cut(i) = 0

      end do

      iread = 0
      irite = 1
      work1(1) = source
      cut(source) = -1

120   continue

      iread = iread + 1

      if ( iread .le. irite ) then

        ip = work1(iread)
        lst = work2(ip) - 1
        i = iarray(ip) - 1

130     continue

          i = i + 1

          if ( lst .lt. i ) then
            go to 120
          end if

          iq = iendpt(2,i)
          del = icpflo(1,i) - icpflo(2,i)

          if ( cut(iq) .eq. 0 .and. del .ne. 0 ) then

            if ( iq .ne. sink ) then
              irite = irite + 1
              work1(irite) = iq
            end if

            cut(iq) = -1

          end if

        go to 130

      end if

      if ( cut(sink) .eq. 0 ) then

        do i = 1, nnode
          cut(i) = -cut(i)
        end do

        do i = 1, nedge
          ip = iendpt(2,iendpt(1,i))
          if ( icpflo(2,i) .lt. 0 ) then
            node_flow(ip) = node_flow(ip) - icpflo(2,i)
          end if
          iendpt(1,i) = ip
        end do

        node_flow(source) = node_flow(sink)
        isort = 1
        go to 10

      end if

      cut(sink) = 1

160   continue

      iread = iread - 1

      if ( iread .eq. 0 ) then
        go to 180
      end if

      ip = work1(iread)
      ien1 = iarray(ip) - 1
      ien2 = work2(ip) - 1

170   continue

      if ( ien1 .ne. ien2 ) then

        iq = iendpt(2,ien2)

        if ( cut(iq) .le. 0 .or. 
     &    icpflo(1,ien2) .eq. icpflo(2,ien2) ) then
          ien2 = ien2 - 1
          go to 170
        end if

        iendpt(2,ien2) = -iq
        icpflo(1,ien2) = icpflo(1,ien2) - icpflo(2,ien2)
        icpflo(2,ien2) = 0
        ien1 = ien1 + 1

        if ( ien1 .lt. ien2 ) then
          go to 80
        end if

      end if

      if ( iarray(ip) .le. ien1 ) then
        cut(ip) = ien1
      end if

      go to 160

180   continue

      kz = 0

      do ir = 1, irite
        if ( 0 .lt. cut(work1(ir)) ) then
          kz = kz + 1
          work1(kz) = work1(ir)
        end if
      end do

      indx = -1
      m = 1

200   continue

      ip = work1(m)

      if ( 0 .lt. node_flow(ip) ) then
        go to 250
      end if

210   continue

      m = m + 1

      if ( m .le. kz ) then
        go to 200
      end if

      iparm = 0

220   continue

      m = m - 1

      if ( m .eq. 1 ) then

        do i = 1, nedge

          iq = -iendpt(2,i)

          if ( 0 .le. iq ) then

            iendpt(2,i) = iq
            j = iendpt(1,i)
            icpflo(1,i) = icpflo(1,i) - icpflo(2,j)

            del = icpflo(2,i) - icpflo(2,j)
            icpflo(2,i) = del
            icpflo(2,j) = -del

          end if

        end do

        go to 100

      end if

      ip = work1(m)

      if ( node_flow(ip) .lt. 0 ) then
        go to 220
      end if

      if ( node_flow(ip) .eq. 0 ) then

        lst = nedge + 1

        if ( ip .lt. nnode ) then
          lst = iarray(ip+1)
        end if

        i = work2(ip)
        work2(ip) = lst

240     continue

          if ( i .eq. lst ) then
            go to 220
          end if

          j = iendpt(1,i)
          del = icpflo(2,j)
          icpflo(2,j) = 0
          icpflo(1,j) = icpflo(1,j) - del
          icpflo(2,i) = icpflo(2,i) - del
          i = i + 1

        go to 240

      end if

      if ( cut(ip) .lt. iarray(ip) ) then
        go to 300
      end if

250   continue

      i = cut(ip) + 1

260   continue

        i = i - 1

        if ( i .lt. iarray(ip) ) then
          go to 290
        end if

        iq = -iendpt(2,i)

        if ( 0 .le. node_flow(iq) ) then
          go to 270
        end if

      go to 260

270   continue

      del = min ( icpflo(1,i) - icpflo(2,i), node_flow(ip) )
      icpflo(2,i) = icpflo(2,i) + del
      node_flow(ip) = node_flow(ip) - del
      node_flow(iq) = node_flow(iq) + del
      iparm = 1
      ien1 = iendpt(1,i)
      ien2 = work2(iq) - 1

      if ( ien1 .lt. ien2 ) then
        go to 80
      end if

      if ( ien1 .eq. ien2 ) then
        work2(iq) = ien2
      end if

280   continue

      if ( 0 .lt. node_flow(ip) ) then
        go to 260
      end if

      if ( icpflo(1,i) .eq. icpflo(2,i) ) then
        i = i - 1
      end if

290   continue

      cut(ip) = i

      if ( iparm .ne. 0 ) then
        go to 210
      end if

300   continue

      i = work2(ip)

310   continue

        j = iendpt(1,i)
        del = min ( icpflo(2,j), node_flow(ip) )
        icpflo(2,j) = icpflo(2,j) - del
        node_flow(ip) = node_flow(ip) - del
        iq = iendpt(2,i)
        node_flow(iq) = node_flow(iq) + del
        i = i + 1

        if ( node_flow(ip) .le. 0 ) then
          go to 320
        end if

      go to 310

320   continue

      node_flow(ip) = -1
      go to 220

      end
      subroutine nim_sum ( i, j, k )

c*********************************************************************72
c
cc NIM_SUM computes the Nim sum of two integers.
c
c  Discussion:
c
c    If K is the Nim sum of I and J, then each bit of K is the exclusive
c    OR of the corresponding bits of I and J.
c
c  Example:
c
c     I     J     K     I base 2    J base 2    K base 2
c   ----  ----  ----  ----------  ----------  ----------
c      0     0     0           0           0           0
c      1     0     1           1           0           1
c      1     1     0           1           1           0
c      2     7     5          10         111         101
c     11    28    23        1011       11100       10111
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    20 January 2007
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer I, J, the integers to be Nim-summed.
c
c    Output, integer K, the Nim sum of I and J.
c
      implicit none

      integer nbits
      parameter ( nbits = 32 )

      integer i
      integer ivec(nbits)
      integer j
      integer jvec(nbits)
      integer k
      integer kvec(nbits)

      call ui4_to_ubvec ( i, nbits, ivec )
      call ui4_to_ubvec ( j, nbits, jvec )
      call bvec_xor ( nbits, ivec, jvec, kvec )
      call ubvec_to_ui4 ( nbits, kvec, k )

      return
      end
      subroutine padovan ( n, p )

c*********************************************************************72
c
cc PADOVAN returns the first N values of the Padovan sequence.
c
c  Discussion:
c
c    The Padovan sequence has the initial values:
c
c      P(0) = 1
c      P(1) = 1
c      P(2) = 1
c
c    and subsequent entries are generated by the recurrence
c
c      P(I+1) = P(I-1) + P(I-2)
c
c  Example:
c
c    0   1
c    1   1
c    2   1
c    3   2
c    4   2
c    5   3
c    6   4
c    7   5
c    8   7
c    9   9
c   10  12
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    04 January 2007
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Ian Stewart,
c    "A Neglected Number",
c    Scientific American, 
c    Volume 274, pages 102-102, June 1996.
c
c    Ian Stewart,
c    Math Hysteria,
c    Oxford, 2004.
c
c  Parameters:
c
c    Input, integer N, the number of terms.
c
c    Output, integer P(N), terms 0 though N-1 of the sequence.
c
      implicit none

      integer n

      integer i
      integer p(n)

      if ( n .lt. 1 ) then
        return
      end if

      p(1) = 1

      if ( n .lt. 2 ) then
        return
      end if

      p(2) = 1

      if ( n .lt. 3 ) then
        return
      end if
 
      p(3) = 1

      do i = 4, n
        p(i) = p(i-2) + p(i-3)
      end do

      return
      end
      subroutine pell_basic ( d, x0, y0 )

c*********************************************************************72
c
cc PELL_BASIC returns the fundamental solution for Pell's basic equation.
c
c  Discussion:
c
c    Pell's equation has the form:
c
c      X**2 - D * Y**2 = 1
c
c    where D is a given non-square integer, and X and Y may be assumed
c    to be positive integers.
c
c  Example:
c
c     D   X0   Y0
c
c     2    3    2
c     3    2    1
c     5    9    4
c     6    5    2
c     7    8    3
c     8    3    1
c    10   19    6
c    11   10    3
c    12    7    2
c    13  649  180
c    14   15    4
c    15    4    1
c    17   33    8
c    18   17    4
c    19  170   39
c    20    9    2
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
c   John Burkardt
c
c  Reference:
c
c    Mark Herkommer,
c    Number Theory, A Programmer's Guide,
c    McGraw Hill, 1999,
c    ISBN: 0-07-913074-7.
c
c  Parameters:
c
c    Input, integer D, the coefficient in Pell's equation.  D should be
c    positive, and not a perfect square.
c
c    Output, integer X0, Y0, the fundamental or 0'th solution.
c    If X0 = Y0 = 0, then the calculation was canceled because of an error.
c    Both X0 and Y0 will be nonnegative.
c
      implicit none

      integer max_term
      parameter ( max_term = 100 )

      integer b(0:max_term)
      integer d
      integer i
      integer n_term
      integer p
      integer pm1
      integer pm2
      integer q
      integer qm1
      integer qm2
      integer r
      integer x0
      integer y0
c
c  If these values are returned, an error has occurred.
c
      x0 = 0
      y0 = 0
c
c  Check D.
c
      if ( d .le. 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'PELL_BASIC - Fatal error!'
        write ( *, '(a)' ) '  Pell coefficient D <= 0.'
        stop
      end if

      call i4_sqrt ( d, q, r )

      if ( r .eq. 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'PELL_BASIC - Fatal error!'
        write ( *, '(a)' ) '  Pell coefficient is a perfect square.'
        stop
      end if
c
c  Find the continued fraction representation of sqrt ( D ).
c
      call i4_sqrt_cf ( d, max_term, n_term, b )
c
c  If necessary, go for two periods.
c
      if ( mod ( n_term, 2 ) .eq. 1 ) then

        do i = n_term + 1, 2*n_term
          b(i) = b(i-n_term)
        end do

        n_term = 2 * n_term

      end if
c
c  Evaluate the continued fraction using the forward recursion algorithm.
c
      pm2 = 0
      pm1 = 1
      qm2 = 1
      qm1 = 0

      do i = 0, n_term-1
        p = b(i) * pm1 + pm2
        q = b(i) * qm1 + qm2
        pm2 = pm1
        pm1 = p
        qm2 = qm1
        qm1 = q
      end do
c
c  Get the fundamental solution.
c
      x0 = p
      y0 = q

      return
      end
      subroutine pell_next ( d, x0, y0, xn, yn, xnp1, ynp1 )

c*****************************************************************************80
c
cc PELL_NEXT returns the next solution of Pell's equation.
c
c  Discussion:
c
c    Pell's equation has the form:
c
c      X**2 - D * Y**2 = 1
c
c    where D is a given non-square integer, and X and Y may be assumed
c    to be positive integers.
c
c    To compute X0, Y0, call PELL_BASIC.
c    To compute X1, Y1, call this routine, with XN and YN set to X0 and Y0.
c    To compute further solutions, call again with X0, Y0 and the previous
c    solution.
c
c  Example:
c
c    ------INPUT--------  --OUTPUT--
c
c    D  X0  Y0   XN   YN  XNP1  YNP1
c
c    2   3   2    3    2    17    12
c    2   3   2   17   12    99    70
c    2   3   2   99   70   577   408
c    2   3   2  577  408  3363  2378
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
c   John Burkardt
c
c  Reference:
c
c    Mark Herkommer,
c    Number Theory, A Programmer's Guide,
c    McGraw Hill, 1999,
c    ISBN: 0-07-913074-7.
c
c  Parameters:
c
c    Input, integer D, the coefficient in Pell's equation.
c
c    Input, integer X0, Y0, the fundamental or 0'th solution.
c
c    Input, integer XN, YN, the N-th solution.
c
c    Output, integer XNP1, YNP1, the N+1-th solution.
c
      implicit none

      integer d
      integer x0
      integer xn
      integer xnp1
      integer y0
      integer yn
      integer ynp1

      xnp1 = x0 * xn + d * y0 * yn
      ynp1 = x0 * yn +     y0 * xn

      return
      end
      subroutine pent_enum ( n, p )

c*********************************************************************72
c
cc PENT_ENUM computes the N-th pentagonal number.
c
c  Discussion:
c
c    The pentagonal number P(N) counts the number of dots in a figure of
c    N nested pentagons.  The pentagonal numbers are defined for both
c    positive and negative N.
c
c    The pentagonal numbers are also useful in determining the
c    number of partitions of an integer.
c
c    P(N) = ( N * ( 3 * N - 1 ) ) / 2
c
c  First values:
c
c     N    P
c
c    -5   40
c    -4   26
c    -3   15
c    -2    7
c    -1    2
c     0    0
c     1    1
c     2    5
c     3   12
c     4   22
c     5   35
c     6   51
c     7   70
c     8   92
c     9  117
c    10  145
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    04 January 2007
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the index of the pentagonal number desired.
c
c    Output, integer P, the value of the N-th pentagonal number.
c
      implicit none

      integer n
      integer p

      p = ( n * ( 3 * n - 1 ) ) / 2

      return
      end
      subroutine perm_ascend ( n, a, length, sub )

c*********************************************************************72
c
cc PERM_ASCEND computes the longest ascending subsequence of a permutation.
c
c  Discussion:
c
c    Although this routine is intended to be applied to a permutation,
c    it will work just as well for an arbitrary vector.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    20 January 2007
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the order of the permutation.
c
c    Input, integer A(N), the permutation to be examined.
c
c    Output, integer LENGTH, the length of the longest increasing subsequence.
c
c    Output, integer SUB(N), contains in entries 1 through LENGTH
c    a longest increasing subsequence of A.
c
      implicit none

      integer n

      integer a(n)
      integer i
      integer j
      integer k
      integer length
      integer sub(n)
      integer top(n)
      integer top_prev(n)

      do i = 1, n
        top(i) = 0
      end do

      do i = 1, n
        top_prev(i) = 0
      end do

      do i = 1, n
        sub(i) = 0
      end do

      if ( n .le. 0 ) then
        length = 0
        return
      end if

      length = 0

      do i = 1, n

        k = 0

        do j = 1, length
          if ( a(i) .le. a(top(j)) ) then
            k = j
            go to 10
          end if
        end do

10      continue

        if ( k .eq. 0 ) then
          length = length + 1
          k = length
        end if

        top(k) = i

        if ( 1 .lt. k ) then
          top_prev(i) = top(k-1)
        else
          top_prev(i) = 0
        end if

      end do

      j = top(length)
      sub(length) = a(j)

      do i = length-1, 1, -1
        j = top_prev(j)
        sub(i) = a(j)
      end do

      return
      end
      subroutine perm_break_count ( n, p, break_count )

c*********************************************************************72
c
cc PERM_BREAK_COUNT counts the number of "breaks" in a permutation.
c
c  Discussion:
c
c    We begin with a permutation of order N.  We prepend an element
c    labeled "0" and append an element labeled "N+1".  There are now
c    N+1 pairs of neighbors.  A "break" is a pair of neighbors whose
c    value differs by more than 1.  
c
c    The identity permutation has a break count of 0.  The maximum
c    break count is N+1.
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
c    Input, integer N, the order of the permutation.
c
c    Input, integer P(N), a permutation, in standard index form.
c
c    Output, integer BREAK_COUNT, the number of breaks in the permutation.
c
      implicit none

      integer n

      integer break_count
      integer i
      integer ierror
      integer p(n)

      break_count = 0
c
c  Make sure the permutation is a legal one.
c  (This is not an efficient way to do soc)
c
      call perm_check ( n, p, ierror )

      if ( p(1) .ne. 1 ) then
        break_count = break_count + 1
      end if

      do i = 1, n-1
        if ( abs ( p(i+1) - p(i) ) .ne. 1 ) then
          break_count = break_count + 1
        end if
      end do

      if ( p(n) .ne. n ) then
        break_count = break_count + 1
      end if

      return
      end
      subroutine perm_canon_to_cycle ( n, p1, p2 )

c*********************************************************************72
c
cc PERM_CANON_TO_CYCLE converts a permutation from canonical to cycle form.
c
c  Example:
c
c    Input:
c
c      4 5 2 1 6 3
c
c    Output:
c
c      -4 5 -2 -1 6 3,
c      indicating the cycle structure
c      ( 4, 5 ) ( 2 ) ( 1, 6, 3 )
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
c  Reference:
c
c    Donald Knuth,
c    The Art of Computer Programming,
c    Volume 1, Fundamental Algorithms,
c    Addison Wesley, 1968, page 176.
c
c  Parameters:
c
c    Input, integer N, the number of objects permuted.
c
c    Input, integer P1(N), the permutation, in canonical form.
c
c    Output, integer P2(N), the permutation, in cycle form.
c
      implicit none

      integer n

      integer i
      integer p1(n)
      integer p2(n)
      integer pmin

      do i = 1, n
        p2(i) = p1(i)
      end do

      pmin = p2(1) + 1

      do i = 1, n

        if ( p2(i) .lt. pmin ) then
          pmin = p2(i)
          p2(i) = -p2(i)
        end if

      end do

      return
      end
      subroutine perm_check ( n, p, ierror )

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
c    Output, integer IERROR, error flag.
c    0, the array does represent a permutation.
c    nonzero, the array does not represent a permutation.  The smallest
c    missing value is equal to IERROR.
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
      subroutine perm_cycle ( n, iopt, p, isgn, ncycle )

c*********************************************************************72
c
cc PERM_CYCLE analyzes a permutation.
c
c  Discussion:
c
c    The routine will count cycles, find the sign of a permutation,
c    and tag a permutation.
c
c  Example:
c
c    Input:
c
c      N = 9
c      IOPT = 1
c      P = 2, 3, 9, 6, 7, 8, 5, 4, 1
c
c    Output:
c
c      NCYCLE = 3
c      ISGN = +1
c      P = -2, 3, 9, -6, -7, 8, 5, 4, 1
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
c    Original FORTRAN77 version by Albert Nijenhuis, Herbert Wilf.
c    This FORTRAN77 version by John Burkardt.
c
c  Reference:
c
c    Albert Nijenhuis, Herbert Wilf,
c    Combinatorial Algorithms for Computers and Calculators,
c    Second Edition,
c    Academic Press, 1978,
c    ISBN: 0-12-519260-6,
c    LC: QA164.N54.
c
c  Parameters:
c
c    Input, integer N, the number of objects being permuted.
c
c    Input, integer IOPT, requests tagging.
c    0, the permutation will not be tagged.
c    1, the permutation will be tagged.
c
c    Input/output, integer P(N).  On input, P describes a
c    permutation, in the sense that entry I is to be moved to P(I).
c    If IOPT = 0, then P will not be changed by this routine.
c    If IOPT = 1, then on output, P will be "tagged".  That is,
c    one element of every cycle in P will be negated.  In this way,
c    a user can traverse a cycle by starting at any entry I1 of P
c    which is negative, moving to I2 = ABS(P(I1)), then to
c    P(I2), and so on, until returning to I1.
c
c    Output, integer ISGN, the "sign" of the permutation, which is
c    +1 if the permutation is even, -1 if odd.  Every permutation
c    may be produced by a certain number of pairwise switches.
c    If the number of switches is even, the permutation itself is
c    called even.
c
c    Output, integer NCYCLE, the number of cycles in the permutation.
c
      implicit none

      integer n

      integer i
      integer i1
      integer i2
      integer ierror
      integer iopt
      integer is
      integer isgn
      integer ncycle
      integer p(n)

      call perm_check ( n, p, ierror )

      is = 1
      ncycle = n

      do i = 1, n

        i1 = p(i)

10      continue

        if ( i .lt. i1 ) then
          ncycle = ncycle - 1
          i2 = p(i1)
          p(i1) = -i2
          i1 = i2
          go to 10
        end if

        if ( iopt .ne. 0 ) then
          is = -sign ( 1, p(i) )
        end if

        p(i) = sign ( p(i), is )

      end do

      isgn = 1 - 2 * mod ( n - ncycle, 2 )

      return
      end
      subroutine perm_cycle_to_canon ( n, p1, p2 )

c*********************************************************************72
c
cc PERM_CYCLE_TO_CANON converts a permutation from cycle to canonical form.
c
c  Discussion:
c
c    The procedure is to "rotate" the elements of each cycle so that
c    the smallest element is first:
c
c      ( 1, 6, 3 ) ( 4, 5 ) ( 2 )
c
c    and then to sort the cycles in decreasing order of their first
c    (and lowest) element:
c
c      ( 4, 5 ) ( 2 ) ( 1, 6, 3 )
c
c    and then to drop the parentheses:
c
c      4 5 2 1 6 3
c
c  Example:
c
c    Input:
c
c      -6 3 1 -5, 4 -2,
c      indicating the cycle structure
c      ( 6, 3, 1 ) ( 5, 4 ) ( 2 )
c
c    Output:
c
c      4 5 2 1 6 3
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
c  Reference:
c
c    Donald Knuth,
c    The Art of Computer Programming,
c    Volume 1, Fundamental Algorithms,
c    Addison Wesley, 1968, pages 176.
c
c  Parameters:
c
c    Input, integer N, the number of objects permuted.
c
c    Input, integer P1(N), the permutation, in cycle form.
c
c    Output, integer P2(N), the permutation, in canonical form.
c
      implicit none

      integer n

      integer hi(n)
      integer i
      integer indx(n)
      integer j
      integer k
      integer lo(n)
      integer ncycle
      integer next
      integer nhi
      integer nlo
      integer nmin
      integer p1(n)
      integer p2(n)
      integer pmin(n)
      integer ptemp(n)

      do i = 1, n
        p2(i) = p1(i)
      end do
c
c  Work on the next cycle.
c
      nlo = 1
      ncycle = 0

10    continue

      if ( nlo .le. n ) then
c
c  Identify NHI, the last index in this cycle.
c
        ncycle = ncycle + 1

        nhi = nlo

20      continue

        if ( nhi .lt. n ) then
          if ( p2(nhi+1) .lt. 0 ) then
            go to 30
          end if
          nhi = nhi + 1
          go to 20
        end if

30      continue
c
c  Identify the smallest value in this cycle.
c
        p2(nlo) = -p2(nlo)
        pmin(ncycle) = p2(nlo)
        nmin = nlo

        do i = nlo+1, nhi
          if ( p2(i) .lt. pmin(ncycle) ) then
            pmin(ncycle) = p2(i)
            nmin = i
          end if
        end do
c
c  Rotate the cycle so A_MIN occurs first.
c
        do i = nlo, nmin-1
          ptemp(nhi+1+i-nmin) = p2(i)
        end do

        do i = nmin, nhi
          ptemp(nlo+i-nmin) = p2(i)
        end do

        lo(ncycle) = nlo
        hi(ncycle) = nhi
c
c  Prepare to operate on the next cycle.
c
        nlo = nhi + 1

        go to 10

      end if
c
c  Compute a sorting index for the cycle minima.
c
      call i4vec_sort_heap_index_d ( ncycle, pmin, indx )
c
c  Copy the cycles out of the temporary array in sorted order.
c
      j = 0
      do i = 1, ncycle
        next = indx(i)
        nlo = lo(next)
        nhi = hi(next)
        do k = nlo, nhi
          j = j + 1
          p2(j) = ptemp(k)
        end do
      end do

      return
      end
      subroutine perm_cycle_to_index ( n, p1, p2 )

c*********************************************************************72
c
cc PERM_CYCLE_TO_INDEX converts a permutation from cycle to standard index form.
c
c  Example:
c
c    Input:
c
c      N = 9
c      P1 = -1, 2, 3, 9, -4, 6, 8, -5, 7
c
c    Output:
c
c      P2 = 2, 3, 9, 6, 7, 8, 5, 4, 1
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
c  Reference:
c
c    Albert Nijenhuis, Herbert Wilf,
c    Combinatorial Algorithms for Computers and Calculators,
c    Second Edition,
c    Academic Press, 1978,
c    ISBN: 0-12-519260-6,
c    LC: QA164.N54.
c
c  Parameters:
c
c    Input, integer N, the number of objects being permuted.
c
c    Input, integer P1(N), the permutation, in cycle form.
c
c    Output, integer P2(N), the permutation, in standard index form.
c
      implicit none

      integer n

      integer j
      integer k1
      integer k2
      integer k3
      integer p1(n)
      integer p2(n)

      do j = 1, n

        k1 = p1(j)

        if ( k1 .lt. 0 ) then
          k1 = -k1
          k3 = k1
        end if

        if ( j + 1 .le. n ) then
          k2 = p1(j+1)
          if ( k2 .lt. 0 ) then
            k2 = k3
          end if
        else
          k2 = k3
        end if

        p2(k1) = k2

      end do

      return
      end
      subroutine perm_distance ( n, a, b, k )

c*********************************************************************72
c
cc PERM_DISTANCE computes the Ulam metric distance of two permutations.
c
c  Discussion:
c
c    If we let N be the order of the permutations A and B, and L(P) be
c    the length of the longest ascending subsequence of a permutation P,
c    then the Ulam metric distance between A and B is
c
c      N - L ( A * inverse ( B ) ).
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
c    Input, integer N, the order of the permutation.
c
c    Input, integer A(N), B(N), the permutations to be examined.
c
c    Output, integer K, the Ulam metric distance between A and B.
c
      implicit none

      integer n

      integer a(n)
      integer b(n)
      integer binv(n)
      integer c(n)
      integer i
      integer k
      integer length
      integer sub(n)

      do i = 1, n
        binv(i) = b(i)
      end do

      call perm_inverse ( n, binv )

      call perm_mul ( n, a, binv, c )

      call perm_ascend ( n, c, length, sub )

      k = n - length

      return
      end
      subroutine perm_fixed_enum ( n, m, fnm )

c*********************************************************************72
c
cc PERM_FIXED_ENUM enumerates the permutations of N objects with M fixed.
c
c  Discussion:
c
c    A permutation of N objects with M fixed is a permutation in which
c    exactly M of the objects retain their original positions.  If
c    M = 0, the permutation is a "derangement".  If M = N, the
c    permutation is the identity.
c
c    F(N,M) = ( Nc / Mc ) * ( 1 - 1/1c + 1/2c - 1/3c ... 1/(N-M)c )
c           = COMB(N,M) * D(N-M)
c
c    where D(N-M) is the number of derangements of N-M objects.
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
c    Input, integer N, the number of objects to be permuted.
c    N should be at least 1.
c
c    Input, integer M, the number of objects that retain their
c    position.  M should be between 0 and N.
c
c    Output, integer FNM, the number of derangements of N objects
c    in which M objects retain their positions.
c
      implicit none

      integer derange_enum
      integer fnm
      integer icnm
      integer m
      integer n

      if ( n .le. 0 ) then

        fnm = 1

      else if ( m .lt. 0 ) then

        fnm = 0

      else if ( n .lt. m ) then

        fnm = 0

      else if ( m .eq. n ) then

        fnm = 1

      else if ( n .eq. 1 ) then

        if ( m .eq. 1 ) then
          fnm = 1
        else
          fnm = 0
        end if

      else

        call combin2 ( n, m, icnm )

        fnm = icnm * derange_enum ( n - m )

      end if

      return
      end
      subroutine perm_free ( npart, ipart, nfree, ifree )

c*********************************************************************72
c
cc PERM_FREE reports the unused items in a partial permutation.
c
c  Discussion:
c
c    It is assumed that the N objects being permuted are the integers
c    from 1 to N, and that IPART contains a "partial" permutation, that
c    is, the NPART entries of IPART represent the beginning of a
c    permutation of all N items.
c
c    The routine returns in IFREE the items that have not been used yet.
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
c    Input, integer NPART, the number of entries in IPART.  NPART may be 0.
c
c    Input, integer IPART(NPART), the partial permutation, which should
c    contain, at most once, some of the integers between 1 and
c    NPART+NFREE.
c
c    Input, integer NFREE, the number of integers that have not been
c    used in IPART.  This is simply N - NPART.  NFREE may be zero.
c
c    Output, integer IFREE(NFREE), the integers between 1 and NPART+NFREE
c    that were not used in IPART.
c
      implicit none

      integer nfree
      integer npart

      integer i
      integer ifree(nfree)
      integer ipart(npart)
      integer j
      integer k
      integer match
      integer n

      n = npart + nfree

      if ( npart .lt. 0 ) then

        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'PERM_FREE - Fatal error!'
        write ( *, '(a)' ) '  NPART .lt. 0.'
        write ( *, '(a,i8)' ) '  NPART = ', npart
        stop

      else if ( npart .eq. 0 ) then

        call i4vec_indicator ( n, ifree )

      else if ( nfree .lt. 0 ) then

        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'PERM_FREE - Fatal error!'
        write ( *, '(a)' ) '  NFREE .lt. 0.'
        write ( *, '(a,i8)' ) '  NFREE = ', nfree
        stop

      else if ( nfree .eq. 0 ) then

        return

      else

        k = 0

        do i = 1, n

          match = 0

          do j = 1, npart
            if ( ipart(j) .eq. i ) then
              match = j
              exit
            end if
          end do

          if ( match .eq. 0 ) then

            k = k + 1

            if ( nfree .lt. k ) then
              write ( *, '(a)' ) ' '
              write ( *, '(a)' ) 'PERM_FREE - Fatal error!'
              write ( *, '(a)' ) '  The partial permutation is illegal.'
              write ( *, '(a)' ) '  It should contain, at most once,'
              write ( *, '(a,i8)' ) 
     &          '  some of the integers between 1 and ', n
              write ( *, '(a)' ) ' '
              write ( *, '(a)' ) '  Our error is that NFREE .lt. K,'
              write ( *, '(a)' ) '  We have TOO MANY missing values.'
              write ( *, '(a,i8)' ) '  Value of NFREE = ', nfree
              write ( *, '(a,i8)' ) '  Value of K =     ', k
              call i4vec_print ( npart, ipart, 
     &          '  Partial permutation:' )
              stop
            end if

            ifree(k) = i

          end if

        end do

      end if

      return
      end
      subroutine perm_index_to_cycle ( n, p1, p2 )

c*********************************************************************72
c
cc PERM_INDEX_TO_CYCLE converts a permutation from standard index to cycle form.
c
c  Example:
c
c    Input:
c
c      N = 9
c      P1 = 2, 3, 9, 6, 7, 8, 5, 4, 1
c
c    Output:
c
c      P2 = -1, 2, 3, 9, -4, 6, 8, -5, 7
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
c  Reference:
c
c    Albert Nijenhuis, Herbert Wilf,
c    Combinatorial Algorithms for Computers and Calculators,
c    Second Edition,
c    Academic Press, 1978,
c    ISBN: 0-12-519260-6,
c    LC: QA164.N54.
c
c  Parameters:
c
c    Input, integer N, the number of objects being permuted.
c
c    Input, integer P1(N), the permutation, in standard index form.
c
c    Output, integer P2(N), the permutation, in cycle form.
c
      implicit none

      integer n

      integer i
      integer j
      integer k
      integer p1(n)
      integer p2(n)

      i = 0
      j = 1

10    continue

      if ( j .le. n ) then

        if ( p1(j) .lt. 0 ) then

          j = j + 1

        else

          k = j

          i = i + 1
          p2(i) = -k

20        continue

          if ( p1(k) .ne. j ) then
            i = i + 1
            p2(i) = p1(k)
            p1(k) = -p1(k)
            k = abs ( p1(k) )
            go to 20
          end if

          p1(k) = -p1(k)

        end if

        go to 10

      end if

      do i = 1, n
        p1(i) = abs ( p1(i) )
      end do

      return
      end
      subroutine perm_ins ( n, p, ins )

c*********************************************************************72
c
cc PERM_INS computes the inversion sequence of a permutation.
c
c  Discussion:
c
c    For a given permutation P acting on objects 1 through N, the inversion
c    sequence INS is defined as:
c
c      INS(1) = 0
c      INS(I) = number of values J .lt. I for which P(I) .lt. P(J).
c
c    The original permutation can be recovered from the inversion sequence.
c
c  Example:
c
c    Input:
c
c      ( 3, 5, 1, 4, 2 )
c
c    Output:
c
c      ( 0, 0, 2, 1, 3 )
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
c  Reference:
c
c    Dennis Stanton, Dennis White,
c    Constructive Combinatorics,
c    Springer, 1986,
c    ISBN: 0387963472.
c
c  Parameters:
c
c    Input, integer N, the number of objects being permuted.
c
c    Input, integer P(N), the permutation, in standard index form.
c    The I-th item has been mapped to P(I).
c
c    Output, integer INS(N), the inversion sequence of the permutation.
c
      implicit none

      integer n

      integer i
      integer ierror
      integer ins(n)
      integer j
      integer p(n)

      call perm_check ( n, p, ierror )

      do i = 1, n
        ins(i) = 0
      end do

      do i = 1, n
        do j = 1, i-1
          if ( p(i) .lt. p(j) ) then
            ins(i) = ins(i) + 1
          end if
        end do
      end do

      return
      end
      subroutine perm_inverse ( n, p )

c*********************************************************************72
c
cc PERM_INVERSE inverts a permutation "in place".
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    25 January 2007
c
c  Parameters:
c
c    Input, integer N, the number of objects being permuted.
c
c    Input/output, integer P(N), the permutation, in standard index form.
c    On output, P describes the inverse permutation
c
      implicit none

      integer n

      integer i
      integer i0
      integer i1
      integer i2
      integer ierror
      integer is
      integer p(n)

      if ( n .le. 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'PERM_INVERSE - Fatal error!'
        write ( *, '(a,i8)' ) '  Input value of N = ', n
        stop
      end if

      call perm_check ( n, p, ierror )

      is = 1

      do i = 1, n

        i1 = p(i)

10      continue

        if ( i .lt. i1 ) then
          i2 = p(i1)
          p(i1) = -i2
          i1 = i2
          go to 10
        end if

        is = -sign ( 1, p(i) )
        p(i) = sign ( p(i), is )

      end do

      do i = 1, n

        i1 = -p(i)

        if ( 0 .le. i1 ) then

          i0 = i

20        continue

            i2 = p(i1)
            p(i1) = i0

            if ( i2 .lt. 0 ) then
              go to 30
            end if

            i0 = i1
            i1 = i2

          go to 20

30        continue

        end if

      end do

      return
      end
      subroutine perm_inverse2 ( n, p )

c*****************************************************************************
c
cc PERM_INVERSE2 inverts a permutation "in place".
c
c  Discussion:
c
c    The routine needs no extra vector storage in order to compute the
c    inverse of a permutation.
c
c    This feature might be useful if the permutation is large.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    25 January 2007
c
c  Parameters:
c
c    Input, integer N, the number of objects in the permutation.
c
c    Input/output, integer P(N), the permutation, in standard index form.
c    On output, the inverse permutation.
c
      implicit none

      integer n

      integer i
      integer ierror
      integer ii
      integer j
      integer k
      integer m
      integer p(n)

      call perm_check ( n, p, ierror )

      do ii = 1, n

        m = n + 1 - ii
        i = p(m)

        if ( i .lt. 0 ) then

          p(m) = -i

        else if ( i .ne. m ) then

          k = m

10        continue

            j = p(i)
            p(i) = -k

            if ( j .eq. m ) then
              p(m) = i
              go to 20
            end if

            k = i
            i = j

          go to 10

20        continue

        end if

      end do

      return
      end
      subroutine perm_inverse3 ( n, perm, perm_inv )

c*********************************************************************72
c
cc PERM_INVERSE3 produces the inverse of a given permutation.
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
c    Input, integer N, the number of items permuted.
c
c    Input, integer PERM(N), a permutation.
c
c    Output, integer PERM_INV(N), the inverse permutation.
c
      implicit none

      integer n

      integer i
      integer perm(n)
      integer perm_inv(n)

      do i = 1, n
        perm_inv(perm(i)) = i
      end do

      return
      end
      subroutine perm_lex_next ( n, p, more )

c*********************************************************************72
c
cc PERM_LEX_NEXT generates permutations in lexical order, one at a time.
c
c  Example:
c
c    N = 3
c
c    1   1 2 3
c    2   1 3 2
c    3   2 1 3
c    4   2 3 1
c    5   3 1 2
c    6   3 2 1
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    25 January 2007
c
c  Reference:
c
c    Mok-Kong Shen,
c    Algorithm 202: Generation of Permutations in Lexicographical Order,
c    Communications of the ACM,
c    Volume 6, September 1963, page 517.
c
c  Parameters:
c
c    Input, integer N, the number of elements being permuted.
c
c    Input/output, integer P(N); on first call with MORE = FALSE,
c    this value is not used.  Otherwise, the input value is the previous
c    permutation.  The output value is the next permutation.
c
c    Input/output, logical MORE.
c    On the first call, set MORE = FALSE, to request initialization.
c    On return, if MORE is TRUE, another permutation has been
c    computed and returned, while if MORE is FALSE, there are no more
c    permutations.
c
      implicit none

      integer n

      integer j
      integer k
      logical more
      integer p(n)
      integer u
      integer w
c
c  Initialization.
c
      if ( .not. more ) then

        call i4vec_indicator ( n, p )
        more = .true.

      else

        if ( n .le. 1 ) then
          more = .false.
          return
        end if

        w = n

10      continue

        if ( p(w) .lt. p(w-1) ) then

          if ( w .eq. 2 ) then
            more = .false.
            return
          end if

          w = w - 1

          go to 10

        end if

        u = p(w-1)

        do j = n, w, -1

          if ( u .lt. p(j) ) then

            p(w-1) = p(j)
            p(j) = u

            do k = 0, ( n - w - 1 ) / 2
              call i4_swap ( p(n-k), p(w+k) )
            end do

            return

          end if

        end do

      end if

      return
      end
      subroutine perm_mul ( n, p1, p2, p3 )

c*********************************************************************72
c
cc PERM_MUL "multiplies" two permutations.
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
c    Input, integer N, the order of the permutations.
c
c    Input, integer P1(N), P2(N), the permutations, in standard index form.
c
c    Output, integer P3(N), the product permutation.
c
      implicit none

      integer n

      integer i
      integer ierror
      integer p1(n)
      integer p2(n)
      integer p3(n)

      call perm_check ( n, p1, ierror )

      call perm_check ( n, p2, ierror )

      do i = 1, n
        p3(i) = p2(p1(i))
      end do

      return
      end
      subroutine perm_next ( n, p, more, even )

c*********************************************************************72
c
cc PERM_NEXT computes all of the permutations of N objects, one at a time.
c
c  Discussion:
c
c    The routine is initialized by calling with MORE = TRUE, in which case
c    it returns the identity permutation.
c
c    If the routine is called with MORE = FALSE, then the successor of the
c    input permutation is computed.
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
c    Original FORTRAN77 version by Albert Nijenhuis, Herbert Wilf.
c    This FORTRAN77 version by John Burkardt.
c
c  Reference:
c
c    Albert Nijenhuis, Herbert Wilf,
c    Combinatorial Algorithms for Computers and Calculators,
c    Second Edition,
c    Academic Press, 1978,
c    ISBN: 0-12-519260-6,
c    LC: QA164.N54.
c
c  Parameters:
c
c    Input, integer N, the number of objects being permuted.
c
c    Input/output, integer P(N), the permutation, in standard index form.
c    On the first call, the input value is unimportant.
c    On subsequent calls, the input value should be the same
c    as the output value from the previous call.  In other words, the
c    user should just leave P alone.
c    On output, contains the "next" permutation.
c
c    Input/output, logical MORE.
c    Set MORE = FALSE before the first call.
c    MORE will be reset to TRUE and a permutation will be returned.
c    Each new call produces a new permutation until
c    MORE is returned FALSE.
c
c    Input/output, logical EVEN.
c    The input value of EVEN should simply be its output value from the
c    previous call; (the input value on the first call doesn't matter.)
c    On output, EVEN is TRUE if the output permutation is even, that is,
c    involves an even number of transpositions.
c
      implicit none

      integer n

      logical even
      integer i
      integer i1
      integer ia
      integer id
      integer is
      integer j
      integer l
      integer m
      logical more
      integer p(n)

      if ( .not. more ) then

        call i4vec_indicator ( n, p )
        more = .true.
        even = .true.

        if ( n .eq. 1 ) then
          more = .false.
          return
        end if

        if ( p(n) .ne. 1 .or. p(1) .ne. 2 + mod ( n, 2 ) ) then
          return
        end if

        do i = 1, n-3
          if ( p(i+1) .ne. p(i)+1 ) then
            return
          end if
        end do

        more = .false.

      else

        if ( n .eq. 1 ) then
          p(1) = 0
          more = .false.
          return
        end if

        if ( even ) then

          ia = p(1)
          p(1) = p(2)
          p(2) = ia
          even = .false.

          if ( p(n) .ne. 1 .or. p(1) .ne. 2 + mod ( n, 2 ) ) then
            return
          end if

          do i = 1, n-3
            if ( p(i+1) .ne. p(i)+1 ) then
              return
            end if
          end do

          more = .false.
          return

        else

          more = .false.

          is = 0

          do i1 = 2, n

            ia = p(i1)
            i = i1 - 1
            id = 0

            do j = 1, i
              if ( ia .lt. p(j) ) then
                id = id + 1
              end if
            end do

            is = id + is
            if ( id .ne. i * mod ( is, 2 ) ) then
              more = .true.
              exit
            end if

          end do

          if ( .not. more ) then
            p(1) = 0
            return
          end if

        end if

        m = mod ( is+1, 2 ) * ( n + 1 )

        do j = 1, i

          if ( sign ( 1, p(j)-ia ) .ne. sign ( 1, p(j)-m ) ) then
            m = p(j)
            l = j
          end if

        end do

        p(l) = ia
        p(i1) = m
        even = .true.

      end if

      return
      end
      subroutine perm_next2 ( n, p, done, invers )

c*********************************************************************72
c
cc PERM_NEXT2 generates all the permutations of N objects.
c
c  Discussion:
c
c    The routine generates the permutations one at a time.  It uses a
c    particular ordering of permutations, generating them from the first
c    (which is the identity permutation) to the Nc-th.  The same ordering
c    is used by the routines PERM_RANK and PERM_UNRANK.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    26 January 2007
c
c  Reference:
c
c    Dennis Stanton, Dennis White,
c    Constructive Combinatorics,
c    Springer, 1986,
c    ISBN: 0387963472.
c
c  Parameters:
c
c    Input, integer N, the number of elements in the set to be permuted.
c
c    Input/output, integer P(N), the permutation, in standard index form.
c
c    Input/output, logical DONE.  The user should set the input value of
c    DONE only once, before the first call to compute the permutations.
c    The user should set DONE to TRUE, which signals the routine
c    that it is to initialize itself.
c    Thereafter, the routine will set DONE to FALSE and will
c    compute a new permutation on each call.
c    However, when there are no more permutations to compute, the
c    routine will not return a new permutation, but instead will
c    return DONE with the value TRUE.  At this point, all the
c    permutations have been computed.
c
c    Output, integer INVERS(N), the inverse permutation of P.
c
c  Local Parameters:
c
c    Local, integer ACTIVE(N), DIR(N).
c
      implicit none

      integer n
      integer n_max
      parameter ( n_max = 100 )

      integer active(n_max)
      integer dir(n_max)
      logical done
      integer i
      integer invers(n)
      integer j
      integer nactiv
      integer p(n)

      save active
      save dir

      if ( n_max < n ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'PERM_NEXT2 - Fatal error!'
        write ( *, '(a)' ) '  Input N exceeds internal limit.'
        stop
      end if
c
c  An input value of TRUE for DONE is assumed to mean a new
c  computation is beginning.
c
      if ( done ) then

        call i4vec_indicator ( n, p )

        do i = 1, n
          invers(i) = p(i)
        end do

        do i = 1, n
          dir(i) = -1
        end do

        active(1) = 0
        do i = 2, n
          active(i) = 1
        end do
c
c  Set the DONE flag to FALSE, signifying there are more permutations
c  to come.  Except, of course, that we must take care of the trivial casec
c
        if ( 1 .lt. n ) then
          done = .false.
        else
          done = .true.
        end if
c
c  Otherwise, assume we are in a continuing computation
c
      else

        nactiv = 0

        do i = 1, n
          if ( active(i) .ne. 0 ) then
            nactiv = i
          end if
        end do

        if ( nactiv .le. 0 ) then

          done = .true.

        else

          j = invers(nactiv)

          p(j) = p(j+dir(nactiv))
          p(j+dir(nactiv)) = nactiv

          invers(nactiv) = invers(nactiv) + dir(nactiv)
          invers(p(j)) = j

          if ( j + 2 * dir(nactiv) .lt. 1 .or. 
     &         n .lt. j + 2 * dir(nactiv) ) then
            dir(nactiv) = -dir(nactiv)
            active(nactiv) = 0
          else if ( nactiv .lt. p(j+2*dir(nactiv)) ) then
            dir(nactiv) = -dir(nactiv)
            active(nactiv) = 0
          end if

          do i = nactiv+1, n
            active(i) = 1
          end do

        end if

      end if

      return
      end
      subroutine perm_next3 ( n, p, more )

c*********************************************************************72
c
cc PERM_NEXT3 computes all of the permutations of N objects, one at a time.
c
c  Discussion:
c
c    The routine is initialized by calling with MORE = TRUE, in which case
c    it returns the identity permutation.
c
c    If the routine is called with MORE = FALSE, then the successor of the
c    input permutation is computed.
c
c    Trotter's algorithm is used.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    26 January 2007
c
c  Reference:
c
c    Hale Trotter,
c    Algorithm 115:
c    PERM,
c    Communications of the Association for Computing Machinery,
c    Volume 5, 1962, pages 434-435.
c
c  Parameters:
c
c    Input, integer N, the number of objects being permuted.
c
c    Input/output, integer P(N), the permutation, in standard index form.
c    If MORE is TRUE, then P is assumed to contain the
c    "previous" permutation, and on P(I) is the value
c    of the I-th object under the next permutation.
c    Otherwise, P will be set to the "first" permutation.
c
c    Input/output, logical MORE.
c    Set MORE = FALSE before first calling this routine.
c    MORE will be reset to TRUE and a permutation will be returned.
c    Each new call produces a new permutation until MORE is returned FALSE.
c
      implicit none

      integer n

      integer i4_factorial
      integer m2
      logical more
      integer n2
      integer nfact
      integer p(n)
      integer q
      integer rank
      integer s
      integer t

      save nfact
      save rank

      data nfact / 0 /
      data rank / 0 /

      if ( .not. more ) then

        call i4vec_indicator ( n, p )
        more = .true.
        rank = 1

        nfact = i4_factorial ( n )

      else

        n2 = n
        m2 = rank
        s = n

10      continue

          q = mod ( m2, n2 )
          t = mod ( m2, 2 * n2 )

          if ( q .ne. 0 ) then
            go to 20
          end if

          if ( t .eq. 0 ) then
            s = s - 1
          end if

          m2 = m2 / n2
          n2 = n2 - 1

        go to 10

20      continue

        if ( q .eq. t ) then
          s = s - q
        else
          s = s + q - n2
        end if

        call i4_swap ( p(s), p(s+1) )

        rank = rank + 1

        if ( rank .eq. nfact ) then
          more = .false.
        end if

      end if

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
      subroutine perm_random ( n, seed, p )

c*********************************************************************72
c
cc PERM_RANDOM selects a random permutation of N objects.
c
c  Discussion:
c
c    The routine assumes the objects are labeled 1, 2, ... N.
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
c    Original FORTRAN77 version by Albert Nijenhuis, Herbert Wilf.
c    This FORTRAN77 version by John Burkardt.
c
c  Reference:
c
c    Albert Nijenhuis, Herbert Wilf,
c    Combinatorial Algorithms for Computers and Calculators,
c    Second Edition,
c    Academic Press, 1978,
c    ISBN: 0-12-519260-6,
c    LC: QA164.N54.
c
c  Parameters:
c
c    Input, integer N, the number of objects to be permuted.
c
c    Input/output, integer SEED, a seed for the random number generator.
c
c    Output, integer P(N), a permutation of ( 1, 2, ..., N ), in standard 
c    index form.
c
      implicit none

      integer n

      integer i
      integer i4_uniform
      integer j
      integer p(n)
      integer seed

      call i4vec_indicator ( n, p )

      do i = 1, n
        j = i4_uniform ( i, n, seed )
        call i4_swap ( p(i), p(j) )
      end do

      return
      end
      subroutine perm_random2 ( n, seed, p )

c*********************************************************************72
c
cc PERM_RANDOM2 selects a random permutation of N objects.
c
c  Discussion:
c
c    The input values of P are used as labels; that is, the I-th object 
c    is labeled P(I).
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
c    Original FORTRAN77 version by Albert Nijenhuis, Herbert Wilf.
c    This FORTRAN77 version by John Burkardt.
c
c  Reference:
c
c    Albert Nijenhuis, Herbert Wilf,
c    Combinatorial Algorithms for Computers and Calculators,
c    Second Edition,
c    Academic Press, 1978,
c    ISBN: 0-12-519260-6,
c    LC: QA164.N54.
c
c  Parameters:
c
c    Input, integer N, the number of objects to be permuted.
c
c    Input/output, integer SEED, a seed for the random number generator.
c
c    Input/output, integer P(N), on input, a list of labels.  On output,
c    the list has been permuted randomly.
c
      implicit none

      integer n

      integer i
      integer i4_uniform
      integer j
      integer p(n)
      integer seed

      do i = 1, n
        j = i4_uniform ( i, n, seed )
        call i4_swap ( p(i), p(j) )
      end do

      return
      end
      subroutine perm_random3 ( n, seed, p )

c*********************************************************************72
c
cc PERM_RANDOM3 selects a random permutation of N elements.
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
c    Original FORTRAN77 version by James Filliben.
c    This FORTRAN77 version by John Burkardt.
c
c  Reference:
c
c    Karla Hoffman, Douglas Shier,
c    Algorithm 564:
c    A Test Problem Generator for Discrete Linear L1 Approximation Problems,
c    ACM Transactions on Mathematical Software,
c    Volume 6, Number 4, December 1980, pages 615-617.
c
c  Parameters:
c
c    Input, integer N, the number of elements of the array.
c
c    Input/output, integer SEED, a seed for the random number generator.
c
c    Output, integer P(N), a permutation, in standard index form.
c
      implicit none

      integer n

      integer i
      integer i4_uniform
      integer iadd
      integer j
      integer p(n)
      integer seed

      if ( n .lt. 1 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'PERM_RANDOM3 - Fatal error!'
        write ( *, '(a,i8)' ) '  Illegal input value of N  = ', n
        write ( *, '(a)' ) '  N must be at least 1!'
        stop
      end if

      if ( n .eq. 1 ) then
        p(1) = 1
        return
      end if

      call i4vec_indicator ( n, p )

      do i = 1, n

        iadd = i4_uniform ( 1, n, seed )

        j = i + iadd

        if ( n .lt. j ) then
          j = j - n
        end if

        if ( i .ne. j ) then
          call i4_swap ( p(j), p(i) )
        end if

      end do

      return
      end
      subroutine perm_rank ( n, p, rank )

c*********************************************************************72
c
cc PERM_RANK computes the rank of a given permutation.
c
c  Discussion:
c
c    This is the same as asking for the step at which PERM_NEXT2
c    would compute the permutation.  The value of the rank will be
c    between 1 and Nc.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    25 January 2007
c
c  Reference:
c
c    Dennis Stanton, Dennis White,
c    Constructive Combinatorics,
c    Springer, 1986,
c    ISBN: 0387963472.
c
c  Parameters:
c
c    Input, integer N, the number of elements in the set that
c    is permuted by P.
c
c    Input, integer P(N), a permutation, in standard index form.
c
c    Output, integer RANK, the rank of the permutation.  This
c    gives the order of the given permutation in the set of all
c    the permutations on N elements.
c
      implicit none

      integer n

      integer count
      integer i
      integer ierror
      integer invers(n)
      integer j
      integer p(n)
      integer rank
      integer rem
c
c  Make sure the permutation is a legal one.
c
      call perm_check ( n, p, ierror )
c
c  Compute the inverse permutation.
c
      do i = 1, n
        invers(i) = p(i)
      end do

      call perm_inverse2 ( n, invers )

      rank = 0

      do i = 1, n

        count = 0

        do j = 1, invers(i)
          if ( p(j) .lt. i ) then
            count = count + 1
          end if
        end do

        if ( mod ( rank, 2 ) .eq. 1 ) then
          rem = count
        else
          rem = i - 1 - count
        end if

        rank = i * rank + rem

      end do

      rank = rank + 1

      return
      end
      subroutine perm_sign ( n, p, p_sign )

c*********************************************************************72
c
cc PERM_SIGN returns the sign of a permutation.
c
c  Discussion:
c
c    A permutation can always be replaced by a sequence of pairwise
c    transpositions.  A given permutation can be represented by
c    many different such transposition sequences, but the number of
c    such transpositions will always be odd or always be even.
c    If the number of transpositions is even or odd, the permutation is
c    said to be even or odd.
c
c  Example:
c
c    Input:
c
c      N = 9
c      P = 2, 3, 9, 6, 7, 8, 5, 4, 1
c
c    Output:
c
c      P_SIGN = +1
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
c    Original FORTRAN77 version by Albert Nijenhuis, Herbert Wilf.
c    This FORTRAN77 version by John Burkardt.
c
c  Reference:
c
c    Albert Nijenhuis, Herbert Wilf,
c    Combinatorial Algorithms for Computers and Calculators,
c    Second Edition,
c    Academic Press, 1978,
c    ISBN: 0-12-519260-6,
c    LC: QA164.N54.
c
c  Parameters:
c
c    Input, integer N, the number of objects permuted.
c
c    Input, integer P(N), a permutation, in standard index form.
c
c    Output, integer P_SIGN, the "sign" of the permutation.
c    +1, the permutation is even,
c    -1, the permutation is odd.
c
      implicit none

      integer n

      integer i
      integer ierror
      integer i4vec_index
      integer j
      integer p(n)
      integer p_sign
      integer q(n)

      call perm_check ( n, p, ierror )
c
c  Make a temporary copy of the permutation.
c
      do i = 1, n
        q(i) = p(i)
      end do
c
c  Start with P_SIGN indicating an even permutation.
c  Restore each element of the permutation to its correct position,
c  updating P_SIGN as you go.
c
      p_sign = 1

      do i = 1, n-1

        j = i4vec_index ( n, q, i )

        if ( j .ne. i ) then
          call i4_swap ( q(i), q(j) )
          p_sign = -p_sign
        end if

      end do

      return
      end
      subroutine perm_to_equiv ( n, p, npart, jarray, iarray )

c*********************************************************************72
c
cc PERM_TO_EQUIV computes the partition induced by a permutation.
c
c  Example:
c
c    Input:
c
c      N = 9
c      P = 2, 3, 9, 6, 7, 8, 5, 4, 1
c
c    Output:
c
c      NPART = 3
c      JARRAY = 4, 3, 2
c      IARRAY = 1, 1, 1, 2  3  2  3  2, 1
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    26 January 2007
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the number of objects being permuted.
c
c    Input, integer P(N), a permutation, in standard index form.
c
c    Output, integer NPART, number of subsets in the partition.
c
c    Output, integer JARRAY(N).  JARRAY(I) is the number of elements
c    in the I-th subset of the partition.
c
c    Output, integer IARRAY(N).  IARRAY(I) is the class to which
c    element I belongs.
c
      implicit none

      integer n

      integer i
      integer iarray(n)
      integer ierror
      integer j
      integer jarray(n)
      integer k
      integer npart
      integer p(n)

      call perm_check ( n, p, ierror )
c
c  Initialize.
c
      do i = 1, n
        iarray(i) = 0
      end do

      do i = 1, n
        jarray(i) = 0
      end do

      npart = 0
c
c  Search for the next item J which has not been assigned a subset/orbit.
c
      do j = 1, n

        if ( iarray(j) .ne. 0 ) then
          go to 20
        end if
c
c  Begin a new subset/orbit.
c
        npart = npart + 1
        k = j
c
c  Add the item to the subset/orbit.
c
10      continue

          jarray(npart) = jarray(npart) + 1
          iarray(k) = npart
c
c  Apply the permutation.  If the permuted object isn't already in the
c  subset/orbit, add it.
c
          k = p(k)

          if ( iarray(k) .ne. 0 ) then
            go to 20
          end if

        go to 10

20      continue

      end do

      return
      end
      subroutine perm_to_ytb ( n, p, lambda, a )

c*********************************************************************72
c
cc PERM_TO_YTB converts a permutation to a Young table.
c
c  Discussion:
c
c    The mapping is not invertible.  In most cases, several permutations
c    correspond to the same table.
c
c  Example:
c
c    N = 7
c    P = 7 2 4 1 5 3 6
c
c    YTB =
c      1 2 3 6
c      4 5
c      7
c
c    LAMBDA = 4 2 1 0 0 0 0
c
c    A = 1 1 1 2 2 1 3
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    11 April 2001
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Albert Nijenhuis, Herbert Wilf,
c    Combinatorial Algorithms for Computers and Calculators,
c    Second Edition,
c    Academic Press, 1978,
c    ISBN: 0-12-519260-6,
c    LC: QA164.N54.
c
c  Parameters:
c
c    Input, integer N, the integer to be partitioned.
c
c    Input, integer P(N), a permutation, in standard index form.
c
c    Output, integer LAMBDA(N).  LAMBDA(I) is the length of the I-th row.
c
c    Output, integer A(N).  A(I) is the row containing I.
c
      implicit none

      integer n

      integer a(n)
      logical another
      integer compare
      integer i
      integer lambda(n)
      integer p(n)
      integer put_index
      integer put_row
      integer put_value
c
c  Initialize.
c
      do i = 1, n
        lambda(i) = 0
      end do

      do i = 1, n
        a(i) = 0
      end do
c
c  Now insert each item of the permutation.
c
      do put_index = 1, n

        put_value = p(put_index)
        put_row = 1

10      continue

          another = .false.

          do compare = put_value + 1, n

            if ( a(compare) .eq. put_row ) then
              another = .true.
              a(put_value) = put_row
              a(compare) = 0
              put_value = compare
              put_row = put_row + 1
              go to 20
            end if

          end do

20        continue

          if ( .not. another ) then
            go to 30
          end if

        go to 10

30      continue

        a(put_value) = put_row
        lambda(put_row) = lambda(put_row) + 1

      end do

      return
      end
      subroutine perm_unrank ( n, rank, p )

c*********************************************************************72
c
cc PERM_UNRANK "unranks" a permutation.
c
c  Discussion:
c
c    That is, given a rank, it computes the corresponding permutation.
c    This is the same as asking for the permutation which PERM_NEXT2
c    would compute at the RANK-th step.
c
c    The value of the rank should be between 1 and Nc.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    22 July 2004
c
c  Reference:
c
c    Dennis Stanton, Dennis White,
c    Constructive Combinatorics,
c    Springer, 1986,
c    ISBN: 0387963472.
c
c  Parameters:
c
c    Input, integer N, the number of elements in the set.
c
c    Input, integer RANK, the desired rank of the permutation.  This
c    gives the order of the given permutation in the set of all
c    the permutations on N elements, using the ordering of PERM_NEXT2.
c
c    Output, integer P(N), the permutation, in standard index form.
c
      implicit none

      integer n

      integer i
      integer icount
      integer iprev
      integer irem
      integer j
      integer jdir
      integer jrank
      integer nfact
      integer p(n)
      integer rank

      do i = 1, n
        p(i) = 0
      end do

      nfact = 1

      do i = 1, n
        nfact = nfact * i
      end do

      if ( rank .lt. 1 .or. nfact .lt. rank ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'PERM_UNRANK - Fatal error!'
        write ( *, '(a)' ) '  Illegal input value for RANK.'
        write ( *, '(a,i8)' ) '  RANK must be between 1 and ', nfact
        write ( *, '(a,i8)' ) '  but the input value is ', rank
        stop
      end if

      jrank = rank - 1

      do i = 1, n

        iprev = n + 1 - i
        irem = mod ( jrank, iprev )
        jrank = jrank / iprev

        if ( mod ( jrank, 2 ) .eq. 1 ) then
          j = 0
          jdir = 1
        else
          j = n + 1
          jdir = -1
        end if

        icount = 0

10      continue

          j = j + jdir

          if ( p(j) .eq. 0 ) then
            icount = icount + 1
          end if

          if ( irem .lt. icount ) then
            go to 20
          end if

        go to 10

20      continue

        p(j) = iprev

      end do

      return
      end
      subroutine perrin ( n, p )

c*********************************************************************72
c
cc PERRIN returns the first N values of the Perrin sequence.
c
c  Discussion:
c
c    The Perrin sequence has the initial values:
c
c      P(0) = 3
c      P(1) = 0
c      P(2) = 2
c
c    and subsequent entries are generated by the recurrence
c
c      P(I+1) = P(I-1) + P(I-2)
c
c    Note that if N is a prime, then N must evenly divide P(N).
c
c  Example:
c
c    0   3
c    1   0
c    2   2
c    3   3
c    4   2
c    5   5
c    6   5
c    7   7
c    8  10
c    9  12
c   10  17
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
c  Reference:
c
c    Ian Stewart,
c    "A Neglected Number",
c    Scientific American, 
c    Volume 274, pages 102-102, June 1996.
c
c    Ian Stewart,
c    Math Hysteria,
c    Oxford, 2004.
c
c    Eric Weisstein,
c    CRC Concise Encyclopedia of Mathematics,
c    CRC Press, 2002,
c    Second edition,
c    ISBN: 1584883472,
c    LC: QA5.W45.
c
c  Parameters:
c
c    Input, integer N, the number of terms.
c
c    Output, integer P(N), the terms 0 through N-1 of the sequence.
c
      implicit none

      integer n

      integer i
      integer p(n)

      if ( n .lt. 1 ) then
        return
      end if

      p(1) = 3

      if ( n .lt. 2 ) then
        return
      end if

      p(2) = 0

      if ( n .lt. 3 ) then
        return
      end if
     
      p(3) = 2

      do i = 4, n
        p(i) = p(i-2) + p(i-3)
      end do

      return
      end
      subroutine pord_check ( n, a, ierror )

c*********************************************************************72
c
cc PORD_CHECK checks a matrix representing a partial ordering.
c
c  Discussion:
c
c    The array A is supposed to represent a partial ordering of
c    the elements of a set of N objects.
c
c    For distinct indices I and J, the value of A(I,J) is:
c
c      1, if I << J
c      0, otherwise ( I and J may be unrelated, or perhaps J << I).
c
c    Diagonal elements of A are ignored.
c
c    This routine checks that the values of A do represent
c    a partial ordering.
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
c    Input, integer N, the number of elements in the set.
c
c    Input, integer A(N,N), the partial ordering.
c    1 if I is less than J in the partial ordering, 
c    0 otherwise.
c
c    Output, integer IERROR, error flag.
c    0, no errors detected.  A is a partial ordering.
c    1, N .le. 0.
c    2, 0 .lt. A(I,J) and 0 .lt. A(J,I) for some I and J.
c
      implicit none

      integer n

      integer a(n,n)
      integer i
      integer ierror
      integer j

      ierror = 0

      if ( n .le. 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'PORD_CHECK - Fatal error!'
        write ( *, '(a,i8)' ) '  N must be positive, but N = ', n
        ierror = 1
        stop
      end if

      do i = 1, n
        do j = i+1, n

          if ( 0 .lt. a(i,j) ) then
            if ( 0 .lt. a(j,i) ) then
              write ( *, '(a)' ) ' '
              write ( *, '(a)' ) 'PORD_CHECK - Fatal error!'
              write ( *, '(a,i8)' ) '  For indices I = ', i
              write ( *, '(a,i8)' ) '  and J = ', j
              write ( *, '(a,i8)' ) '  A(I,J) = ', a(i,j)
              write ( *, '(a,i8)' ) '  A(J,I) = ', a(j,i)
              ierror = 2
              stop
            end if
          end if

        end do
      end do

      return
      end
      subroutine power_mod ( a, n, m, x )

c*********************************************************************72
c
cc POWER_MOD computes mod ( A**N, M ).
c
c  Discussion:
c
c    Some programming tricks are used to speed up the computation, and to
c    allow computations in which the value A**N is much too large to 
c    store in an integer word.
c
c    First, for efficiency, the power A**N is computed by determining
c    the binary expansion of N, then computing A, A**2, A**4, and so on
c    by repeated squaring, and multiplying only those factors that
c    contribute to A**N.
c
c    Secondly, the intermediate products are immediately "mod'ed", which
c    keeps them small.
c
c    For instance, to compute mod ( A**13, 11 ), we essentially compute
c
c       13 = 1 + 4 + 8
c
c       A**13 = A * A**4 * A**8
c
c       mod ( A**13, 11 ) = mod ( A, 11 ) * mod ( A**4, 11 ) * mod ( A**8, 11 ).
c
c    Fermat's little theorem says that if P is prime, and A is not divisible
c    by P, then ( A**(P-1) - 1 ) is divisible by P.
c
c  Example:
c
c     A  N  M  X
c
c    13  0 31  1
c    13  1 31 13
c    13  2 31 14
c    13  3 31 27
c    13  4 31 10
c    13  5 31  6
c    13  6 31 16
c    13  7 31 22
c    13  8 31  7 
c    13  9 31 29
c    13 10 31  5
c    13 11 31  3
c    13 12 31  8
c    13 13 31 11
c    13 14 31 19
c    13 15 31 30
c    13 16 31 18
c    13 17 31 17
c    13 18 31  4
c    13 19 31 21
c    13 20 31 25
c    13 21 31 15
c    13 22 31  9
c    13 23 31 24
c    13 24 31  2
c    13 25 31 26
c    13 26 31 28
c    13 27 31 23
c    13 28 31 20
c    13 29 31 12
c    13 30 31  1
c    13 31 31 13
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    26 January 2007
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer A, the base of the expression to be tested.
c    0 .le. A is required.
c
c    Input, integer N, the power to which the base is raised.
c    0 .le. N is required.
c
c    Input, integer M, the divisor against which the expression is tested.
c    0 .lt. M is required.
c
c    Output, integer X, the remainder when A**N is divided by M.
c    If any input quantity is unacceptable, then the nonsensical value
c    X = -1 is returned.
c
      implicit none

      integer a
      integer*8 a_square2
      integer d
      integer m
      integer*8 m2
      integer n
      integer ncopy
      integer x
      integer*8 x2

      if ( a .lt. 0 ) then
        x = -1
        return
      end if

      if ( m .le. 0 ) then
        x = -1
        return
      end if

      if ( n .lt. 0 ) then
        x = -1
        return
      end if
c
c  A_SQUARE2 contains the successive squares of A.
c
      a_square2 = a
      m2 = m
      x2 = 1

      ncopy = n

10    continue

      if ( 0 .lt. ncopy ) then

        d = mod ( ncopy, 2 )

        if ( d .eq. 1 ) then
          x2 = mod ( x2 * a_square2, m2 )
        end if

        a_square2 = mod ( a_square2 * a_square2, m2 )
        ncopy = ( ncopy - d ) / 2

        go to 10

      end if
c
c  Ensure that X is nonnegative.
c
20    continue

      if ( x2 .lt. 0 ) then
        x2 = x2 + m
        go to 20
      end if

      x = x2

      return
      end
      subroutine power_series1 ( n, alpha, a, b )

c*********************************************************************72
c
cc POWER_SERIES1 computes the power series for G(Z) = (1+F(Z))**ALPHA.
c
c  Discussion:
c
c    The power series for F(Z) is given.
c
c    The form of the power series are:
c
c      F(Z) = A1*Z + A2*Z**2 + A3*Z**3 + ... + AN*Z**N
c
c      G(Z) = B1*Z + B2*Z**2 + B3*Z**3 + ... + BN*Z**N
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
c    Original FORTRAN77 version by Albert Nijenhuis, Herbert Wilf.
c    This FORTRAN77 version by John Burkardt.
c
c  Reference:
c
c    Albert Nijenhuis, Herbert Wilf,
c    Combinatorial Algorithms for Computers and Calculators,
c    Second Edition,
c    Academic Press, 1978,
c    ISBN: 0-12-519260-6,
c    LC: QA164.N54.
c
c  Parameters:
c
c    Input, integer N, the number of terms in the power series.
c
c    Input, double precision ALPHA, the exponent of 1+F(Z) in the 
c    definition of G(Z).
c
c    Input, double precision A(N), the power series coefficients for F(Z).
c
c    Output, double precision B(N), the power series coefficients for G(Z).
c
      implicit none

      integer n

      double precision a(n)
      double precision alpha
      double precision b(n)
      integer i
      integer j
      double precision v

      do j = 1, n

        v = 0.0D+00

        do i = 1, j-1
          v = v + b(i) * a(j-i) * ( alpha * ( j - i ) - i )
        end do

        b(j) = ( alpha * a(j) + v / dble ( j ) )

      end do

      return
      end
      subroutine power_series2 ( n, a, b )

c*********************************************************************72
c
cc POWER_SERIES2 computes the power series for G(Z) = exp(F(Z)) - 1.
c
c  Discussion:
c
c    The power series for F(Z) is given.
c
c    The power series have the form:
c
c      F(Z) = A1*Z + A2*Z**2 + A3*Z**3 + ... + AN*Z**N
c
c      G(Z) = B1*Z + B2*Z**2 + B3*Z**3 + ... + BN*Z**N
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
c    Original FORTRAN77 version by Albert Nijenhuis, Herbert Wilf.
c    This FORTRAN77 version by John Burkardt.
c
c  Reference:
c
c    Albert Nijenhuis, Herbert Wilf,
c    Combinatorial Algorithms for Computers and Calculators,
c    Second Edition,
c    Academic Press, 1978,
c    ISBN: 0-12-519260-6,
c    LC: QA164.N54.
c
c  Parameters:
c
c    Input, integer N, the number of terms in the power series.
c
c    Input, double precision A(N), the power series coefficients for F(Z).
c
c    Output, double precision B(N), the power series coefficients for G(Z).
c
      implicit none

      integer n

      double precision a(n)
      double precision b(n)
      integer i
      integer j
      double precision v

      do j = 1, n

        v = 0.0D+00

        do i = 1, j-1
          v = v + b(i) * a(j-i) * dble ( j - i )
        end do

        b(j) = a(j) + v / dble ( j )

      end do

      return
      end
      subroutine power_series3 ( n, a, b, c )

c*********************************************************************72
c
cc POWER_SERIES3 computes the power series for H(Z) = G(F(Z)).
c
c  Discussion:
c
c    The power series for G and H are given.
c
c    We assume that
c
c      F(Z) = A1*Z + A2*Z**2 + A3*Z**3 + ... + AN*Z**N
c      G(Z) = B1*Z + B2*Z**2 + B3*Z**3 + ... + BN*Z**N
c      H(Z) = C1*Z + C2*Z**2 + C3*Z**3 + ... + CN*Z**N
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    26 January 2007
c
c  Author:
c
c    Original FORTRAN77 version by Albert Nijenhuis, Herbert Wilf.
c    This FORTRAN77 version by John Burkardt.
c
c  Reference:
c
c    Albert Nijenhuis, Herbert Wilf,
c    Combinatorial Algorithms for Computers and Calculators,
c    Second Edition,
c    Academic Press, 1978,
c    ISBN: 0-12-519260-6,
c    LC: QA164.N54.
c
c  Parameters:
c
c    Input, integer N, the number of terms in the power series.
c
c    Input, double precision A(N), the power series for F.
c
c    Input, double precision B(N), the power series for G.
c
c    Output, double precision C(N), the power series for H.
c
      implicit none

      integer n
      integer n_max
      parameter ( n_max = 100 )

      double precision a(n)
      double precision b(n)
      double precision c(n)
      integer i
      integer iq
      integer j
      integer m
      double precision r
      double precision v
      double precision work(n_max)

      if ( n_max < n ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'POWER_SERIES3 - Fatal error!'
        write ( *, '(a)' ) '  Input N exceeds internal limit.'
        stop
      end if

      do i = 1, n
        work(i) = b(1) * a(i)
      end do
c
c  Search for IQ, the index of the first nonzero entry in A.
c
      iq = 0

      do i = 1, n

        if ( a(i) .ne. 0.0D+00 ) then
          iq = i
          go to 10
        end if

      end do

10    continue

      if ( iq .ne. 0 ) then

        m = 1

20      continue

          m = m + 1

          if ( n .lt. m * iq ) then
            go to 30
          end if

          if ( b(m) .eq. 0.0D+00 ) then
            go to 20
          end if

          r = b(m) * a(iq)**m
          work(m*iq) = work(m*iq) + r

          do j = 1, n-m*iq

            v = 0.0D+00
            do i = 1, j-1
              v = v + c(i) * a(j-i+iq) * dble ( m * ( j - i ) - i )
            end do

            c(j) = ( dble ( m ) * a(j) + v / dble ( j ) ) / a(iq)

          end do

          do i = 1, n-m*iq
            work(i+m*iq) = work(i+m*iq) + c(i) * r
          end do

        go to 20

30      continue

      end if

      do i = 1, n
        c(i) = work(i)
      end do

      return
      end
      subroutine power_series4 ( n, a, b, c )

c*********************************************************************72
c
cc POWER_SERIES4 computes the power series for H(Z) = G ( 1/F(Z) ).
c
c  Discussion:
c
c    The routine is given the power series for the functions F and G.
c
c    We assume that
c
c      F(Z) = A1*Z + A2*Z**2 + A3*Z**3 + ... + AN*Z**N
c      G(Z) = B1*Z + B2*Z**2 + B3*Z**3 + ... + BN*Z**N
c      H(Z) = C1*Z + C2*Z**2 + C3*Z**3 + ... + CN*Z**N
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    26 January 2007
c
c  Author:
c
c    Original FORTRAN77 version by Albert Nijenhuis, Herbert Wilf.
c    This FORTRAN77 version by John Burkardt.
c
c  Reference:
c
c    Albert Nijenhuis, Herbert Wilf,
c    Combinatorial Algorithms for Computers and Calculators,
c    Second Edition,
c    Academic Press, 1978,
c    ISBN: 0-12-519260-6,
c    LC: QA164.N54.
c
c  Parameters:
c
c    Input, integer N, the number of terms in the power series.
c
c    Input, double precision A(N), the power series for F. 
c    A(1) may not be 0.0.
c
c    Input, double precision B(N), the power series for G.
c
c    Output, double precision C(N), the power series for H.
c
      implicit none

      integer n
      integer n_max
      parameter ( n_max = 100 )

      double precision a(n)
      double precision b(n)
      double precision c(n)
      integer i
      integer j
      integer k
      double precision s
      double precision t
      double precision work(n_max)

      if ( n_max < n ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'POWER_SERIES4 - Fatal error!'
        write ( *, '(a)' ) '  Input N exceeds internal limit.'
        stop
      end if

      if ( a(1) .eq. 0.0D+00 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'POWER_SERIES4 - Fatal error!'
        write ( *, '(a)' ) '  A(1) is zero.'
        stop
      end if

      t = 1.0D+00

      do i = 1, n
        t = t / a(1)
        c(i) = b(i) * t
        work(i) = a(i) * t
      end do

      do k = 2, n
        s = -work(k)
        do i = k, n
          do j = i, n
            c(j) = c(j) + s * c(j+1-k)
            work(j) = work(j) + s * work(j+1-k)
          end do
        end do
      end do

      return
      end
      function prime ( n )

c*********************************************************************72
c
cc PRIME returns any of the first PRIME_MAX prime numbers.
c
c  Discussion:
c
c    PRIME_MAX is 1600, and the largest prime stored is 13499.
c
c    Thanks to Bart Vandewoestyne for pointing out a typo, 18 February 2005.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    21 January 2007
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Milton Abramowitz, Irene Stegun,
c    Handbook of Mathematical Functions,
c    National Bureau of Standards, 1964,
c    ISBN: 0-486-61272-4,
c    LC: QA47.A34.
c
c    Daniel Zwillinger,
c    CRC Standard Mathematical Tables and Formulae,
c    30th Edition,
c    CRC Press, 1996, pages 95-98.
c
c  Parameters:
c
c    Input, integer N, the index of the desired prime number.
c    In general, is should be true that 0 <= N <= PRIME_MAX.
c    N = -1 returns PRIME_MAX, the index of the largest prime available.
c    N = 0 is legal, returning PRIME = 1.
c
c    Output, integer PRIME, the N-th prime.  If N is out of range,
c    PRIME is returned as -1.
c
      implicit none

      integer prime_max
      parameter ( prime_max = 1600 )

      integer i
      integer n
      integer npvec(prime_max)
      integer prime

      save npvec

      data ( npvec(i), i = 1, 100 ) / 
     &      2,    3,    5,    7,   11,   13,   17,   19,   23,   29, 
     &     31,   37,   41,   43,   47,   53,   59,   61,   67,   71, 
     &     73,   79,   83,   89,   97,  101,  103,  107,  109,  113, 
     &    127,  131,  137,  139,  149,  151,  157,  163,  167,  173, 
     &    179,  181,  191,  193,  197,  199,  211,  223,  227,  229, 
     &    233,  239,  241,  251,  257,  263,  269,  271,  277,  281, 
     &    283,  293,  307,  311,  313,  317,  331,  337,  347,  349, 
     &    353,  359,  367,  373,  379,  383,  389,  397,  401,  409, 
     &    419,  421,  431,  433,  439,  443,  449,  457,  461,  463, 
     &    467,  479,  487,  491,  499,  503,  509,  521,  523,  541 /

      data ( npvec(i), i = 101, 200 ) / 
     &    547,  557,  563,  569,  571,  577,  587,  593,  599,  601, 
     &    607,  613,  617,  619,  631,  641,  643,  647,  653,  659, 
     &    661,  673,  677,  683,  691,  701,  709,  719,  727,  733, 
     &    739,  743,  751,  757,  761,  769,  773,  787,  797,  809, 
     &    811,  821,  823,  827,  829,  839,  853,  857,  859,  863, 
     &    877,  881,  883,  887,  907,  911,  919,  929,  937,  941, 
     &    947,  953,  967,  971,  977,  983,  991,  997, 1009, 1013, 
     &   1019, 1021, 1031, 1033, 1039, 1049, 1051, 1061, 1063, 1069, 
     &   1087, 1091, 1093, 1097, 1103, 1109, 1117, 1123, 1129, 1151, 
     &   1153, 1163, 1171, 1181, 1187, 1193, 1201, 1213, 1217, 1223 /

      data ( npvec(i), i = 201, 300 ) / 
     &   1229, 1231, 1237, 1249, 1259, 1277, 1279, 1283, 1289, 1291, 
     &   1297, 1301, 1303, 1307, 1319, 1321, 1327, 1361, 1367, 1373, 
     &   1381, 1399, 1409, 1423, 1427, 1429, 1433, 1439, 1447, 1451, 
     &   1453, 1459, 1471, 1481, 1483, 1487, 1489, 1493, 1499, 1511, 
     &   1523, 1531, 1543, 1549, 1553, 1559, 1567, 1571, 1579, 1583, 
     &   1597, 1601, 1607, 1609, 1613, 1619, 1621, 1627, 1637, 1657, 
     &   1663, 1667, 1669, 1693, 1697, 1699, 1709, 1721, 1723, 1733, 
     &   1741, 1747, 1753, 1759, 1777, 1783, 1787, 1789, 1801, 1811, 
     &   1823, 1831, 1847, 1861, 1867, 1871, 1873, 1877, 1879, 1889, 
     &   1901, 1907, 1913, 1931, 1933, 1949, 1951, 1973, 1979, 1987 /

      data ( npvec(i), i = 301, 400 ) / 
     &   1993, 1997, 1999, 2003, 2011, 2017, 2027, 2029, 2039, 2053, 
     &   2063, 2069, 2081, 2083, 2087, 2089, 2099, 2111, 2113, 2129, 
     &   2131, 2137, 2141, 2143, 2153, 2161, 2179, 2203, 2207, 2213, 
     &   2221, 2237, 2239, 2243, 2251, 2267, 2269, 2273, 2281, 2287, 
     &   2293, 2297, 2309, 2311, 2333, 2339, 2341, 2347, 2351, 2357, 
     &   2371, 2377, 2381, 2383, 2389, 2393, 2399, 2411, 2417, 2423, 
     &   2437, 2441, 2447, 2459, 2467, 2473, 2477, 2503, 2521, 2531, 
     &   2539, 2543, 2549, 2551, 2557, 2579, 2591, 2593, 2609, 2617, 
     &   2621, 2633, 2647, 2657, 2659, 2663, 2671, 2677, 2683, 2687, 
     &   2689, 2693, 2699, 2707, 2711, 2713, 2719, 2729, 2731, 2741 /

      data ( npvec(i), i = 401, 500 ) / 
     &   2749, 2753, 2767, 2777, 2789, 2791, 2797, 2801, 2803, 2819, 
     &   2833, 2837, 2843, 2851, 2857, 2861, 2879, 2887, 2897, 2903, 
     &   2909, 2917, 2927, 2939, 2953, 2957, 2963, 2969, 2971, 2999, 
     &   3001, 3011, 3019, 3023, 3037, 3041, 3049, 3061, 3067, 3079, 
     &   3083, 3089, 3109, 3119, 3121, 3137, 3163, 3167, 3169, 3181, 
     &   3187, 3191, 3203, 3209, 3217, 3221, 3229, 3251, 3253, 3257, 
     &   3259, 3271, 3299, 3301, 3307, 3313, 3319, 3323, 3329, 3331, 
     &   3343, 3347, 3359, 3361, 3371, 3373, 3389, 3391, 3407, 3413, 
     &   3433, 3449, 3457, 3461, 3463, 3467, 3469, 3491, 3499, 3511, 
     &   3517, 3527, 3529, 3533, 3539, 3541, 3547, 3557, 3559, 3571 /

      data ( npvec(i), i = 501, 600 ) / 
     &   3581, 3583, 3593, 3607, 3613, 3617, 3623, 3631, 3637, 3643, 
     &   3659, 3671, 3673, 3677, 3691, 3697, 3701, 3709, 3719, 3727, 
     &   3733, 3739, 3761, 3767, 3769, 3779, 3793, 3797, 3803, 3821, 
     &   3823, 3833, 3847, 3851, 3853, 3863, 3877, 3881, 3889, 3907, 
     &   3911, 3917, 3919, 3923, 3929, 3931, 3943, 3947, 3967, 3989, 
     &   4001, 4003, 4007, 4013, 4019, 4021, 4027, 4049, 4051, 4057, 
     &   4073, 4079, 4091, 4093, 4099, 4111, 4127, 4129, 4133, 4139, 
     &   4153, 4157, 4159, 4177, 4201, 4211, 4217, 4219, 4229, 4231, 
     &   4241, 4243, 4253, 4259, 4261, 4271, 4273, 4283, 4289, 4297, 
     &   4327, 4337, 4339, 4349, 4357, 4363, 4373, 4391, 4397, 4409 /

      data ( npvec(i), i = 601, 700 ) / 
     &   4421, 4423, 4441, 4447, 4451, 4457, 4463, 4481, 4483, 4493, 
     &   4507, 4513, 4517, 4519, 4523, 4547, 4549, 4561, 4567, 4583, 
     &   4591, 4597, 4603, 4621, 4637, 4639, 4643, 4649, 4651, 4657, 
     &   4663, 4673, 4679, 4691, 4703, 4721, 4723, 4729, 4733, 4751, 
     &   4759, 4783, 4787, 4789, 4793, 4799, 4801, 4813, 4817, 4831, 
     &   4861, 4871, 4877, 4889, 4903, 4909, 4919, 4931, 4933, 4937, 
     &   4943, 4951, 4957, 4967, 4969, 4973, 4987, 4993, 4999, 5003, 
     &   5009, 5011, 5021, 5023, 5039, 5051, 5059, 5077, 5081, 5087, 
     &   5099, 5101, 5107, 5113, 5119, 5147, 5153, 5167, 5171, 5179, 
     &   5189, 5197, 5209, 5227, 5231, 5233, 5237, 5261, 5273, 5279 /

      data ( npvec(i), i = 701, 800 ) / 
     &   5281, 5297, 5303, 5309, 5323, 5333, 5347, 5351, 5381, 5387, 
     &   5393, 5399, 5407, 5413, 5417, 5419, 5431, 5437, 5441, 5443, 
     &   5449, 5471, 5477, 5479, 5483, 5501, 5503, 5507, 5519, 5521, 
     &   5527, 5531, 5557, 5563, 5569, 5573, 5581, 5591, 5623, 5639, 
     &   5641, 5647, 5651, 5653, 5657, 5659, 5669, 5683, 5689, 5693, 
     &   5701, 5711, 5717, 5737, 5741, 5743, 5749, 5779, 5783, 5791, 
     &   5801, 5807, 5813, 5821, 5827, 5839, 5843, 5849, 5851, 5857, 
     &   5861, 5867, 5869, 5879, 5881, 5897, 5903, 5923, 5927, 5939, 
     &   5953, 5981, 5987, 6007, 6011, 6029, 6037, 6043, 6047, 6053, 
     &   6067, 6073, 6079, 6089, 6091, 6101, 6113, 6121, 6131, 6133 /

      data ( npvec(i), i = 801, 900 ) / 
     &   6143, 6151, 6163, 6173, 6197, 6199, 6203, 6211, 6217, 6221, 
     &   6229, 6247, 6257, 6263, 6269, 6271, 6277, 6287, 6299, 6301, 
     &   6311, 6317, 6323, 6329, 6337, 6343, 6353, 6359, 6361, 6367, 
     &   6373, 6379, 6389, 6397, 6421, 6427, 6449, 6451, 6469, 6473, 
     &   6481, 6491, 6521, 6529, 6547, 6551, 6553, 6563, 6569, 6571, 
     &   6577, 6581, 6599, 6607, 6619, 6637, 6653, 6659, 6661, 6673, 
     &   6679, 6689, 6691, 6701, 6703, 6709, 6719, 6733, 6737, 6761, 
     &   6763, 6779, 6781, 6791, 6793, 6803, 6823, 6827, 6829, 6833, 
     &   6841, 6857, 6863, 6869, 6871, 6883, 6899, 6907, 6911, 6917, 
     &   6947, 6949, 6959, 6961, 6967, 6971, 6977, 6983, 6991, 6997 /

      data ( npvec(i), i = 901, 1000 ) / 
     &   7001, 7013, 7019, 7027, 7039, 7043, 7057, 7069, 7079, 7103, 
     &   7109, 7121, 7127, 7129, 7151, 7159, 7177, 7187, 7193, 7207, 
     &   7211, 7213, 7219, 7229, 7237, 7243, 7247, 7253, 7283, 7297, 
     &   7307, 7309, 7321, 7331, 7333, 7349, 7351, 7369, 7393, 7411, 
     &   7417, 7433, 7451, 7457, 7459, 7477, 7481, 7487, 7489, 7499, 
     &   7507, 7517, 7523, 7529, 7537, 7541, 7547, 7549, 7559, 7561, 
     &   7573, 7577, 7583, 7589, 7591, 7603, 7607, 7621, 7639, 7643, 
     &   7649, 7669, 7673, 7681, 7687, 7691, 7699, 7703, 7717, 7723, 
     &   7727, 7741, 7753, 7757, 7759, 7789, 7793, 7817, 7823, 7829, 
     &   7841, 7853, 7867, 7873, 7877, 7879, 7883, 7901, 7907, 7919 /

      data ( npvec(i), i = 1001, 1100 ) / 
     &   7927, 7933, 7937, 7949, 7951, 7963, 7993, 8009, 8011, 8017, 
     &   8039, 8053, 8059, 8069, 8081, 8087, 8089, 8093, 8101, 8111, 
     &   8117, 8123, 8147, 8161, 8167, 8171, 8179, 8191, 8209, 8219, 
     &   8221, 8231, 8233, 8237, 8243, 8263, 8269, 8273, 8287, 8291, 
     &   8293, 8297, 8311, 8317, 8329, 8353, 8363, 8369, 8377, 8387, 
     &   8389, 8419, 8423, 8429, 8431, 8443, 8447, 8461, 8467, 8501, 
     &   8513, 8521, 8527, 8537, 8539, 8543, 8563, 8573, 8581, 8597, 
     &   8599, 8609, 8623, 8627, 8629, 8641, 8647, 8663, 8669, 8677, 
     &   8681, 8689, 8693, 8699, 8707, 8713, 8719, 8731, 8737, 8741, 
     &   8747, 8753, 8761, 8779, 8783, 8803, 8807, 8819, 8821, 8831 /

      data ( npvec(i), i = 1101, 1200 ) / 
     &   8837, 8839, 8849, 8861, 8863, 8867, 8887, 8893, 8923, 8929, 
     &   8933, 8941, 8951, 8963, 8969, 8971, 8999, 9001, 9007, 9011, 
     &   9013, 9029, 9041, 9043, 9049, 9059, 9067, 9091, 9103, 9109, 
     &   9127, 9133, 9137, 9151, 9157, 9161, 9173, 9181, 9187, 9199, 
     &   9203, 9209, 9221, 9227, 9239, 9241, 9257, 9277, 9281, 9283, 
     &   9293, 9311, 9319, 9323, 9337, 9341, 9343, 9349, 9371, 9377, 
     &   9391, 9397, 9403, 9413, 9419, 9421, 9431, 9433, 9437, 9439, 
     &   9461, 9463, 9467, 9473, 9479, 9491, 9497, 9511, 9521, 9533, 
     &   9539, 9547, 9551, 9587, 9601, 9613, 9619, 9623, 9629, 9631, 
     &   9643, 9649, 9661, 9677, 9679, 9689, 9697, 9719, 9721, 9733 /

      data ( npvec(i), i = 1201, 1300 ) / 
     &   9739, 9743, 9749, 9767, 9769, 9781, 9787, 9791, 9803, 9811, 
     &   9817, 9829, 9833, 9839, 9851, 9857, 9859, 9871, 9883, 9887, 
     &   9901, 9907, 9923, 9929, 9931, 9941, 9949, 9967, 9973,10007, 
     &  10009,10037,10039,10061,10067,10069,10079,10091,10093,10099, 
     &  10103,10111,10133,10139,10141,10151,10159,10163,10169,10177, 
     &  10181,10193,10211,10223,10243,10247,10253,10259,10267,10271, 
     &  10273,10289,10301,10303,10313,10321,10331,10333,10337,10343, 
     &  10357,10369,10391,10399,10427,10429,10433,10453,10457,10459, 
     &  10463,10477,10487,10499,10501,10513,10529,10531,10559,10567, 
     &  10589,10597,10601,10607,10613,10627,10631,10639,10651,10657 /

      data ( npvec(i), i = 1301, 1400 ) / 
     &  10663,10667,10687,10691,10709,10711,10723,10729,10733,10739, 
     &  10753,10771,10781,10789,10799,10831,10837,10847,10853,10859, 
     &  10861,10867,10883,10889,10891,10903,10909,10937,10939,10949, 
     &  10957,10973,10979,10987,10993,11003,11027,11047,11057,11059, 
     &  11069,11071,11083,11087,11093,11113,11117,11119,11131,11149, 
     &  11159,11161,11171,11173,11177,11197,11213,11239,11243,11251, 
     &  11257,11261,11273,11279,11287,11299,11311,11317,11321,11329, 
     &  11351,11353,11369,11383,11393,11399,11411,11423,11437,11443, 
     &  11447,11467,11471,11483,11489,11491,11497,11503,11519,11527, 
     &  11549,11551,11579,11587,11593,11597,11617,11621,11633,11657 /

      data ( npvec(i), i = 1401, 1500 ) / 
     &  11677,11681,11689,11699,11701,11717,11719,11731,11743,11777, 
     &  11779,11783,11789,11801,11807,11813,11821,11827,11831,11833, 
     &  11839,11863,11867,11887,11897,11903,11909,11923,11927,11933, 
     &  11939,11941,11953,11959,11969,11971,11981,11987,12007,12011, 
     &  12037,12041,12043,12049,12071,12073,12097,12101,12107,12109, 
     &  12113,12119,12143,12149,12157,12161,12163,12197,12203,12211, 
     &  12227,12239,12241,12251,12253,12263,12269,12277,12281,12289, 
     &  12301,12323,12329,12343,12347,12373,12377,12379,12391,12401, 
     &  12409,12413,12421,12433,12437,12451,12457,12473,12479,12487, 
     &  12491,12497,12503,12511,12517,12527,12539,12541,12547,12553 /

      data ( npvec(i), i = 1501, 1600 ) / 
     &  12569,12577,12583,12589,12601,12611,12613,12619,12637,12641, 
     &  12647,12653,12659,12671,12689,12697,12703,12713,12721,12739, 
     &  12743,12757,12763,12781,12791,12799,12809,12821,12823,12829, 
     &  12841,12853,12889,12893,12899,12907,12911,12917,12919,12923, 
     &  12941,12953,12959,12967,12973,12979,12983,13001,13003,13007, 
     &  13009,13033,13037,13043,13049,13063,13093,13099,13103,13109, 
     &  13121,13127,13147,13151,13159,13163,13171,13177,13183,13187, 
     &  13217,13219,13229,13241,13249,13259,13267,13291,13297,13309, 
     &  13313,13327,13331,13337,13339,13367,13381,13397,13399,13411, 
     &  13417,13421,13441,13451,13457,13463,13469,13477,13487,13499 /

      if ( n .eq. -1 ) then
        prime = prime_max
      else if ( n .eq. 0 ) then
        prime = 1
      else if ( n .le. prime_max ) then
        prime = npvec(n)
      else
        prime = -1
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'PRIME - Fatal error!'
        write ( *, '(a,i8)' ) '  Illegal prime index N = ', n
        write ( *, '(a,i8)' ) 
     &    '  N should be between 1 and PRIME_MAX =', prime_max
        stop
      end if

      return
      end
      subroutine pythag_triple_next ( i, j, a, b, c )

c*********************************************************************72
c
cc PYTHAG_TRIPLE_NEXT computes the next Pythagorean triple.
c
c  Example:
c
c     I       J       A       B       C    A^2+B^2 = C^2
c
c     2       1       3       4       5      25
c     3       2       5      12      13     169
c     4       1      15       8      17     289
c     4       3       7      24      25     625
c     5       2      21      20      29     841
c     5       4       9      40      41    1681
c     6       1      35      12      37    1369
c     6       3      27      36      45    2025
c     6       5      11      60      61    3721
c     7       2      45      28      53    2809
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    24 January 2007
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input/output, integer I, J, the generators.
c    On first call, set I = J = 0.  On repeated calls, leave I and J
c    at their output values from the previous call.
c
c    Output, integer A, B, C, the next Pythagorean triple.
c    A, B, and C are positive integers which have no common factors,
c    and A**2 + B**2 = C**2.
c
      implicit none

      integer a
      integer b
      integer c
      integer i
      integer j
c
c  I starts at 2, and when it increases, increases by 1 and resets J;
c
c  When I is reset, J starts out at 2 if I is odd, or 1 if I is even;
c  Then I is held fixed and J increases by 2, as long as it remains less than I.
c
      if ( i .eq. 0 .and. j .eq. 0 ) then
        i = 2
        j = 1
      else if ( j + 2 .lt. i ) then
        j = j + 2
      else
        i = i + 1
        j = mod ( i, 2 ) + 1
      end if

      a = i**2 - j**2
      b = 2 * i * j
      c = i**2 + j**2

      return
      end
      function r8_agm ( a, b )

c*********************************************************************72
c
cc R8_AGM finds the arithmetic-geometric mean of two numbers.
c
c  Discussion:
c
c    The AGM of (A,B) is produced by the following iteration:
c
c      (A,B) -> ( (A+B)/2, SQRT(A*B) ).
c
c    The sequence of successive values of (A,B) quickly converge to the AGM.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    24 January 2007
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, double precision A, B, the numbers whose AGM is desired.  
c    A and B should both be non-negative.
c
c    Output, double precision R8_AGM, the AGM of the two numbers.
c
      implicit none

      double precision a
      double precision a1
      double precision a2
      double precision b
      double precision b1
      double precision b2
      double precision r8_agm
      double precision r8_epsilon
      double precision tol

      if ( a .lt. 0.0D+00 ) then
        r8_agm = -1.0D+00
        return
      end if

      if ( b .lt. 0.0D+00 ) then
        r8_agm = -1.0D+00
        return
      end if

      if ( a .eq. 0.0D+00 .or. b .eq. 0.0D+00 ) then
        r8_agm = 0.0D+00
        return
      end if

      if ( a .eq. b ) then
        r8_agm = a
        return
      end if

      tol = r8_epsilon ( tol ) * ( a + b + 1.0D+00 )

      a1 = a
      b1 = b

10    continue

        a2 = ( a1 + b1 ) / 2.0D+00
        b2 = sqrt ( a1 * b1 )

        if ( abs ( a2 - b2 ) .lt. tol ) then
          r8_agm = ( a2 + b2 ) / 2.0D+00
          go to 20
        end if

        a1 = a2
        b1 = b2

      go to 10

20    continue

      return
      end
      function r8_epsilon ( )

c*********************************************************************72
c
cc R8_EPSILON returns the R8 roundoff unit.
c
c  Discussion:
c
c    The roundoff unit is a number R which is a power of 2 with the
c    property that, to the precision of the computer's arithmetic,
c      1 .lt. 1 + R
c    but
c      1 = ( 1 + R / 2 )
c
c    FORTRAN90 provides the superior library routine
c
c      EPSILON ( X )
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    06 March 2006
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, double precision R8_EPSILON, the R8 roundoff unit.
c
      implicit none

      double precision r8
      double precision r8_epsilon
      double precision r8_test

      r8 = 1.0D+00
      r8_test = 1.0D+00 + ( r8 / 2.0D+00 )

10    continue

      if ( 1.0D+00 .lt. r8_test ) then
        r8 = r8 / 2.0D+00
        r8_test = 1.0D+00 + ( r8 / 2.0D+00 )
        go to 10
      end if

      r8_epsilon = r8

      return
      end
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
      function r8_fall ( x, n )

c*********************************************************************72
c
cc R8_FALL computes the falling factorial function [X]_N.
c
c  Discussion:
c
c    Note that the number of "injections" or 1-to-1 mappings from
c    a set of N elements to a set of M elements is [M]_N.
c
c    The number of permutations of N objects out of M is [M]_N.
c
c    Moreover, the Stirling numbers of the first kind can be used
c    to convert a falling factorial into a polynomial, as follows:
c
c      [X]_N = S^0_N + S^1_N * X + S^2_N * X^2 + ... + S^N_N X^N.
c
c    The formula used is:
c
c      [X]_N = X * ( X - 1 ) * ( X - 2 ) * ... * ( X - N + 1 ).
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    03 January 2007
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, double precision X, the argument of the falling factorial function.
c
c    Input, integer N, the order of the falling factorial function.
c    If N = 0, FALL = 1, if N = 1, FALL = X.  Note that if N is
c    negative, a "rising" factorial will be computed.
c
c    Output, double precision R8_FALL, the value of the falling 
c    factorial function.
c
      implicit none

      double precision arg
      integer i
      integer n
      double precision r8_fall
      double precision value
      double precision x

      value = 1.0D+00

      arg = x

      if ( 0 .lt. n ) then

        do i = 1, n
          value = value * arg
          arg = arg - 1.0D+00
        end do

      else if ( n .lt. 0 ) then

        do i = -1, n, -1
          value = value * arg
          arg = arg + 1.0D+00
        end do

      end if

      r8_fall = value

      return
      end
      function r8_huge ( )

c*********************************************************************72
c
cc R8_HUGE returns a "huge" R8.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    13 April 2004
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, double precision R8_HUGE, a huge number.
c
      implicit none

      double precision r8_huge

      r8_huge = 1.0D+30

      return
      end
      function r8_is_int ( r )

c*********************************************************************72
c
cc R8_IS_INT determines if an R8 represents an integer value.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    21 March 2001
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, double precision R, the number to be checked.
c
c    Output, logical R8_IS_INT, is TRUE if R is an integer value.
c
      implicit none

      integer i
      integer i4_huge
      double precision r
      logical r8_is_int

      if ( dble ( i4_huge ( ) ) .lt. r ) then
        r8_is_int = .false.
      else if ( r .lt. - dble ( i4_huge ( ) ) ) then
        r8_is_int = .false.
      else if ( r .eq. dble ( int ( r ) ) ) then
        r8_is_int = .true.
      else
        r8_is_int = .false.
      end if

      return
      end
      function r8_pi ( )

c*********************************************************************72
c
cc R8_PI returns the value of pi as an R8.
c
c  Modified:
c
c    01 January 2007
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, double precision R8_PI, the value of pi.
c
      implicit none

      double precision r8_pi
  
      r8_pi = 3.141592653589793D+00

      return
      end
      function r8_rise ( x, n )

c*********************************************************************72
c
cc R8_RISE computes the rising factorial function [X]^N.
c
c  Discussion:
c
c    [X]^N = X * ( X + 1 ) * ( X + 2 ) * ... * ( X + N - 1 ).
c
c    Note that the number of ways of arranging N objects in M ordered
c    boxes is [M]^N.  (Here, the ordering of the objects in each box matters).  
c    Thus, 2 objects in 2 boxes have the following 6 possible arrangements:
c
c      -|12, 1|2, 12|-, -|21, 2|1, 21|-.
c
c    Moreover, the number of non-decreasing maps from a set of
c    N to a set of M ordered elements is [M]^N / Nc.  Thus the set of
c    nondecreasing maps from (1,2,3) to (a,b,c,d) is the 20 elements:
c
c      aaa, abb, acc, add, aab, abc, acd, aac, abd, aad
c      bbb, bcc, bdd, bbc, bcd, bbd, ccc, cdd, ccd, ddd.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    24 January 2007
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, double precision X, the argument of the rising factorial function.
c
c    Input, integer N, the order of the rising factorial function.
c    If N = 0, RISE = 1, if N = 1, RISE = X.  Note that if N is
c    negative, a "falling" factorial will be computed.
c
c    Output, double precision R8_RISE, the value of the rising factorial 
c    function.
c
      implicit none

      double precision arg
      integer i
      integer n
      double precision r8_rise
      double precision value
      double precision x

      value = 1.0D+00

      arg = x

      if ( 0 .lt. n ) then

        do i = 1, n
          value = value * arg
          arg = arg + 1.0D+00
        end do

      else if ( n .lt. 0 ) then

        do i = -1, n, -1
          value = value * arg
          arg = arg - 1.0D+00
        end do

      end if

      r8_rise = value

      return
      end
      subroutine r8_swap ( x, y )

c*********************************************************************72
c
cc R8_SWAP switches two R8's.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    30 November 1998
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input/output, double precision X, Y.  On output, the values of X and
c    Y have been interchanged.
c
      implicit none

      double precision x
      double precision y
      double precision z

      z = x
      x = y
      y = z

      return
      end
      subroutine r8_to_cfrac ( r, n, a, p, q )

c*********************************************************************72
c
cc R8_TO_CFRAC converts a real value to a continued fraction.
c
c  Discussion:
c
c    The routine is given a real number R.  It computes a sequence of
c    continued fraction approximations to R, returning the results as
c    simple fractions of the form P(I) / Q(I).
c
c  Example:
c
c    X = 2 * PI
c    N = 7
c
c    A = [ *, 6,  3,  1,  1,   7,   2,    146,      3 ]
c    P = [ 1, 6, 19, 25, 44, 333, 710, 103993, 312689 ]
c    Q = [ 0, 1,  3,  4,  7,  53, 113,  16551,  49766 ]
c
c    (This ignores roundoff error, which will cause later terms to differ).
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    20 January 2007
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Norman Richert,
c    Strang's Strange Figures,
c    American Mathematical Monthly,
c    Volume 99, Number 2, February 1992, pages 101-107.
c
c  Parameters:
c
c    Input, double precision R, the real value.
c
c    Input, integer N, the number of convergents to compute.
c
c    Output, integer A(0:N), the partial quotients.
c
c    Output, integer P(-1:N), Q(-1:N), the numerators and denominators
c    of the continued fraction approximations.
c
      implicit none

      integer n

      integer a(0:n)
      integer i
      integer p(-1:n)
      integer q(-1:n)
      double precision r
      double precision r_copy
      double precision x(0:n)

      if ( r .eq. 0.0D+00 ) then
        do i = 0, n
          a(i) = 0
        end do
        do i = -1, n
          p(i) = 0
        end do
        do i = -1, n
          q(i) = 1
        end do
        return
      end if

      r_copy = abs ( r )

      p(-1) = 1
      q(-1) = 0

      p(0) = int ( r_copy )
      q(0) = 1
      x(0) = r_copy
      a(0) = int ( x(0) )

      do i = 1, n
        x(i) = 1.0D+00 / ( x(i-1) - dble ( a(i-1) ) )
        a(i) = int ( x(i) )
        p(i) = a(i) * p(i-1) + p(i-2)
        q(i) = a(i) * q(i-1) + q(i-2)
      end do

      if ( r .lt. 0.0D+00 ) then
        do i = -1, n
          p(i) = -p(i)
        end do
      end if

      return
      end
      subroutine r8_to_dec ( dval, dec_digit, mantissa, exponent )

c*********************************************************************72
c
cc R8_TO_DEC converts a real quantity to a decimal representation.
c
c  Discussion:
c
c    Given the real ( kind = 8 ) value DVAL, the routine computes integers
c    MANTISSA and EXPONENT so that it is approximately true that:
c
c      DVAL = MANTISSA * 10 ** EXPONENT
c
c    In particular, only DEC_DIGIT digits of DVAL are used in constructing the
c    representation.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    24 January 2007
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, double precision DVAL, the value whose decimal representation
c    is desired.
c
c    Input, integer DEC_DIGIT, the number of decimal digits to use.
c
c    Output, integer MANTISSA, EXPONENT, the approximate decimal 
c    representation of DVAL.
c
      implicit none

      integer dec_digit
      double precision dval
      integer exponent
      integer mantissa
      double precision mantissa_double
      double precision ten1
      double precision ten2
c
c  Special cases.
c
      if ( dval .eq. 0.0D+00 ) then
        mantissa = 0
        exponent = 0
        return
      end if
c
c  Factor DVAL = MANTISSA_DOUBLE * 10**EXPONENT
c
      mantissa_double = dval
      exponent = 0
c
c  Now normalize so that 
c  10**(DEC_DIGIT-1) <= ABS(MANTISSA_DOUBLE) < 10**(DEC_DIGIT)
c
      ten1 = 10.0D+00**( dec_digit - 1 )
      ten2 = 10.0D+00**dec_digit

10    continue

      if ( abs ( mantissa_double ) .lt. ten1 ) then
        mantissa_double = mantissa_double * 10.0D+00
        exponent = exponent - 1
        go to 10
      end if

20    continue

      if ( ten2 .le. abs ( mantissa_double ) ) then
        mantissa_double = mantissa_double / 10.0D+00
        exponent = exponent + 1
        go to 20
      end if
c
c  MANTISSA is the integer part of MANTISSA_DOUBLE, rounded.
c
      mantissa = nint ( mantissa_double )
c
c  Now divide out any factors of ten from MANTISSA.
c
      if ( mantissa .ne. 0 ) then

30      continue

        if ( 10 * ( mantissa / 10 ) .eq. mantissa ) then
          mantissa = mantissa / 10
          exponent = exponent + 1
          go to 30
        end if

      end if

      return
      end
      subroutine r8_to_rat ( r, ndig, iatop, iabot )

c*********************************************************************72
c
cc R8_TO_RAT converts a real value to a rational value.
c
c  Discussion:
c
c    The rational value (IATOP/IABOT) is essentially computed by truncating
c    the decimal representation of the real value after a given number of
c    decimal digits.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    24 January 2007
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, double precision R, the real value to be converted.
c
c    Input, integer NDIG, the number of decimal digits used.
c
c    Output, integer IATOP, IABOT, the numerator and denominator
c    of the rational value that approximates the real number.
c
      implicit none

      double precision factor
      integer i4_gcd
      integer iabot
      integer iatop
      integer ibot
      integer ifac
      integer itemp
      integer itop
      integer jfac
      integer ndig
      double precision r

      factor = 10.0D+00**ndig

      if ( 0 .lt. ndig ) then
        ifac = 10**ndig
        jfac = 1
      else
        ifac = 1
        jfac = 10**(-ndig)
      end if

      itop = nint ( r * factor ) * jfac
      ibot = ifac
c
c  Factor out the greatest common factor.
c
      itemp = i4_gcd ( itop, ibot )

      iatop = itop / itemp
      iabot = ibot / itemp

      return
      end
      subroutine r8_to_s_left ( rval, s )

c*********************************************************************72
c
cc R8_TO_S_LEFT represents a real using 14 left_justified characters.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    24 January 2007
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, real ( kind = 8 ) RVAL, a real number.
c
c    Output, character * ( * ) S, a left-justified character variable
c    containing the representation of RVAL.
c
      implicit none

      character * ( 14 ) chrtmp
      integer i
      double precision rval
      character * ( * ) s
c
c  We can't seem to write directly into the string because of compiler
c  quibbles.
c
      if ( dble ( int ( rval ) ) .eq. rval .and. 
     &     abs ( rval ) .lt. 1.0D+13 ) then

        write ( chrtmp, '(i14)' ) int ( rval )

      else

        write ( chrtmp, '(g14.6)' ) rval

      end if

      do i = 1, len ( chrtmp )
        if ( chrtmp(i:i) .ne. ' ' ) then
          s = chrtmp(i:)
          return
        end if
      end do

      s = ' '

      return
      end
      function r8_uniform ( a, b, seed )

c*********************************************************************72
c
cc R8_UNIFORM returns a scaled pseudorandom R8.
c
c  Discussion:
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
c    06 January 2006
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, double precision A, B, the limits of the interval.
c
c    Input/output, integer SEED, the "seed" value, which should NOT be 0.
c    On output, SEED has been updated.
c
c    Output, double precision R8_UNIFORM, a number strictly between A and B.
c
      implicit none

      double precision a
      double precision b
      integer i4_huge
      integer k
      double precision r8_uniform
      integer seed

      if ( seed .eq. 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'R8_UNIFORM - Fatal error!'
        write ( *, '(a)' ) '  Input value of SEED = 0.'
        stop
      end if

      k = seed / 127773

      seed = 16807 * ( seed - k * 127773 ) - k * 2836

      if ( seed .lt. 0 ) then
        seed = seed + i4_huge ( )
      end if
c
c  Although SEED can be represented exactly as a 32 bit integer,
c  it generally cannot be represented exactly as a 32 bit real number!
c
      r8_uniform = a + ( b - a ) * dble ( seed ) * 4.656612875D-10

      return
      end
      function r8_uniform_01 ( seed )

c*********************************************************************72
c
cc R8_UNIFORM_01 returns a unit pseudorandom R8.
c
c  Discussion:
c
c    This routine implements the recursion
c
c      seed = 16807 * seed mod ( 2**31 - 1 )
c      r8_uniform_01 = seed / ( 2**31 - 1 )
c
c    The integer arithmetic never requires more than 32 bits,
c    including a sign bit.
c
c    If the initial seed is 12345, then the first three computations are
c
c      Input     Output      R8_UNIFORM_01
c      SEED      SEED
c
c         12345   207482415  0.096616
c     207482415  1790989824  0.833995
c    1790989824  2035175616  0.947702
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    11 August 2004
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
c    Peter Lewis, Allen Goodman, James Miller,
c    A Pseudo-Random Number Generator for the System/360,
c    IBM Systems Journal,
c    Volume 8, pages 136-143, 1969.
c
c  Parameters:
c
c    Input/output, integer SEED, the "seed" value, which should NOT be 0.
c    On output, SEED has been updated.
c
c    Output, double precision R8_UNIFORM_01, a new pseudorandom variate,
c    strictly between 0 and 1.
c
      implicit none

      integer i4_huge
      integer k
      double precision r8_uniform_01
      integer seed

      if ( seed .eq. 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'R8_UNIFORM_01 - Fatal error!'
        write ( *, '(a)' ) '  Input value of SEED = 0.'
        stop
      end if

      k = seed / 127773

      seed = 16807 * ( seed - k * 127773 ) - k * 2836

      if ( seed .lt. 0 ) then
        seed = seed + i4_huge ( )
      end if
c
c  Although SEED can be represented exactly as a 32 bit integer,
c  it generally cannot be represented exactly as a 32 bit real number!
c
      r8_uniform_01 = dble ( seed ) * 4.656612875D-10

      return
      end
      subroutine r8mat_det ( n, a, det )

c*********************************************************************72
c
cc R8MAT_DET finds the determinant of an N by N R8MAT.
c
c  Discussion:
c
c    A brute force calculation is made.
c
c    This routine should only be used for small matrices, since this
c    calculation requires the summation of Nc products of N numbers.
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
c    Input, integer N, the number of rows and columns of A.
c
c    Input, double precision A(N,N), the matrix whose determinant is desired.
c
c    Output, double precision DET, the determinant of the matrix.
c
      implicit none

      integer n
      integer n_max
      parameter ( n_max = 100 )

      double precision a(n,n)
      double precision det
      logical even
      integer i
      integer iarray(n_max)
      logical more
      double precision term

      if ( n_max .lt. n ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'R8MAT_DET - Fatal error!'
        write ( *, '(a)' ) '  Input N exceeds internal limit.'
        stop
      end if

      more = .false.
      det = 0.0D+00

10    continue

        call perm_next ( n, iarray, more, even )

        if ( even ) then
          term = 1.0D+00
        else
          term = -1.0D+00
        end if

        do i = 1, n
          term = term * a(i,iarray(i))
        end do

        det = det + term

        if ( .not. more ) then
          go to 20
        end if

      go to 10

20    continue

      return
      end
      subroutine r8mat_perm ( n, a, p )

c*********************************************************************72
c
cc R8MAT_PERM permutes the rows and columns of a square R8MAT.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    26 January 2007
c
c  Author:
c
c    Original FORTRAN77 version by Albert Nijenhuis, Herbert Wilf.
c    This FORTRAN77 version by John Burkardt.
c
c  Reference:
c
c    Albert Nijenhuis, Herbert Wilf,
c    Combinatorial Algorithms for Computers and Calculators,
c    Second Edition,
c    Academic Press, 1978,
c    ISBN: 0-12-519260-6,
c    LC: QA164.N54.
c
c  Parameters:
c
c    Input, integer N, the order of the matrix.
c
c    Input/output, double precision A(N,N).
c    On input, the matrix to be permuted.
c    On output, the permuted matrix.
c
c    Input, integer P(N), a permutation to be applied to the rows
c    and columns.  P(I) is the new number of row and column I.
c
      implicit none

      integer n

      double precision a(n,n)
      integer i
      integer i1
      integer is
      double precision it
      integer j
      integer j1
      integer j2
      integer k
      integer lc
      integer nc
      integer p(n)

      call perm_cycle ( n, 1, p, is, nc )

      do i = 1, n

        i1 = -p(i)

        if ( 0 .lt. i1 ) then

          lc = 0

10        continue

            i1 = p(i1)
            lc = lc + 1

            if ( i1 .le. 0 ) then
              go to 20
            end if

          go to 10

20        continue

          i1 = i

          do j = 1, n

            if ( p(j) .le. 0 ) then

              j2 = j
              k = lc

30            continue

                j1 = j2
                it = a(i1,j1)

40              continue

                  i1 = abs ( p(i1) )
                  j1 = abs ( p(j1) )

                  call r8_swap ( a(i1,j1), it )

                  if ( j1 .ne. j2 ) then
                    go to 40
                  end if

                  k = k - 1

                  if ( i1 .eq. i ) then
                    go to 50
                  end if

                go to 40

50              continue

                j2 = abs ( p(j2) )

                if ( k .eq. 0 ) then
                  go to 60
                end if

              go to 30

60            continue

            end if

          end do

        end if

      end do
c
c  Restore the positive signs of the data.
c
      do i = 1, n
        p(i) = abs ( p(i) )
      end do

      return
      end
      subroutine r8mat_perm2 ( m, n, a, p, q )

c*********************************************************************72
c
cc R8MAT_PERM2 permutes rows and columns of a rectangular R8MAT, in place.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    26 January 2007
c
c  Author:
c
c    Original FORTRAN77 version by Albert Nijenhuis, Herbert Wilf.
c    This FORTRAN77 version by John Burkardt.
c
c  Reference:
c
c    Albert Nijenhuis, Herbert Wilf,
c    Combinatorial Algorithms for Computers and Calculators,
c    Second Edition,
c    Academic Press, 1978,
c    ISBN: 0-12-519260-6,
c    LC: QA164.N54.
c
c  Parameters:
c
c    Input, integer M, number of rows in the matrix.
c
c    Input, integer N, number of columns in the matrix.
c
c    Input/output, double precision A(M,N).
c    On input, the matrix to be permuted.
c    On output, the permuted matrix.
c
c    Input, integer P(M), the row permutation.  P(I) is the new number of row I.
c
c    Input, integer Q(N), the column permutation.  Q(I) is the new number of
c    column I.  Note that the routine allows you to pass a single array as both
c    P and Q.
c
      implicit none

      integer m
      integer n

      double precision a(m,n)
      integer i
      integer i1
      integer is
      integer j
      integer j1
      integer j2
      integer k
      integer lc
      integer nc
      integer p(m)
      integer q(n)
      double precision t

      call perm_cycle ( m, 1, p, is, nc )

      if ( 0 .lt. q(1) ) then
        call perm_cycle ( n, 1, q, is, nc )
      end if

      do i = 1, m

        i1 = -p(i)

        if ( 0 .lt. i1 ) then

          lc = 0

10        continue

            i1 = p(i1)
            lc = lc + 1

            if ( i1 .le. 0 ) then
              go to 20
            end if

          go to 10

20        continue

          i1 = i

          do j = 1, n

            if ( q(j) .le. 0 ) then

              j2 = j
              k = lc

30            continue

                j1 = j2
                t = a(i1,j1)

40              continue

                  i1 = abs ( p(i1) )
                  j1 = abs ( q(j1) )

                  call r8_swap ( a(i1,j1), t )

                  if ( j1 .ne. j2 ) then
                    go to 40
                  end if

                  k = k - 1

                  if ( i1 .eq. i ) then
                    go to 50
                  end if

                go to 40

50              continue

                j2 = abs ( q(j2) )

                if ( k .eq. 0 ) then
                  go to 60
                end if

              go to 30

60            continue

            end if

          end do

        end if

      end do

      do i = 1, m
        p(i) = abs ( p(i) )
      end do

      if ( q(1) .le. 0 ) then

        do i = 1, n
          q(i) = abs ( q(i) )
        end do

      end if

      return
      end
      subroutine r8mat_permanent ( n, a, perm )

c*********************************************************************72
c
cc R8MAT_PERMANENT computes the permanent of an R8MAT.
c
c  Discussion:
c
c    The permanent function is similar to the determinant.  Recall that
c    the determinant of a matrix may be defined as the sum of all the
c    products:
c
c      S * A(1,I(1)) * A(2,I(2)) * ... * A(N,I(N))
c
c    where I is any permutation of the columns of the matrix, and S is the
c    sign of the permutation.  By contrast, the permanent function is
c    the (unsigned) sum of all the products
c
c      A(1,I(1)) * A(2,I(2)) * ... * A(N,I(N))
c
c    where I is any permutation of the columns of the matrix.  The only
c    difference is that there is no permutation sign multiplying each summand.
c
c    Symbolically, then, the determinant of a 2 by 2 matrix
c
c      a b
c      c d
c
c    is a*d-b*c, whereas the permanent of the same matrix is a*d+b*c.
c
c
c    The permanent is invariant under row and column permutations.
c    If a row or column of the matrix is multiplied by S, then the
c      permanent is likewise multiplied by S.
c    If the matrix is square, then the permanent is unaffected by
c      transposing the matrix.
c    Unlike the determinant, however, the permanent does change if
c      one row is added to another, and it is not true that the
c      permanent of the product is the product of the permanents.
c
c
c    Note that if A is a matrix of all 1's and 0's, then the permanent
c    of A counts exactly which permutations hit exactly 1's in the matrix.
c    This fact can be exploited for various combinatorial purposes.
c
c    For instance, setting the diagonal of A to 0, and the offdiagonals
c    to 1, the permanent of A counts the number of derangements of N
c    objects.
c
c    Setting the diagonal of A to 0, and ensuring that the offdiagonal
c    entries are symmetric, then A is the adjacency matrix of a graph,
c    and its permanent counts the number of perfect matchings.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    26 January 2007
c
c  Author:
c
c    Original FORTRAN77 version by Albert Nijenhuis, Herbert Wilf.
c    This FORTRAN77 version by John Burkardt.
c
c  Reference:
c
c    Albert Nijenhuis, Herbert Wilf,
c    Combinatorial Algorithms for Computers and Calculators,
c    Second Edition,
c    Academic Press, 1978,
c    ISBN: 0-12-519260-6,
c    LC: QA164.N54.
c
c  Parameters:
c
c    Input, integer N, the order of the matrix.
c
c    Input, double precision A(N,N), the value of the matrix.
c
c    Output, double precision PERM, the value of the permanent of the matrix.
c
      implicit none

      integer n
      integer n_max
      parameter ( n_max = 100 )

      double precision a(n,n)
      integer i
      integer iadd
      integer iwork(n_max)
      integer j
      logical more
      integer ncard
      double precision p
      double precision perm
      double precision prod
      double precision sgn
      double precision total
      double precision work(n_max)
      double precision z

      if ( n_max < n ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'R8MAT_PERMANENT - Fatal error!'
        write ( *, '(a)' ) '  Input N exceeds internal limit.'
        stop
      end if

      more = .false.

      do i = 1, n
        total = 0.0D+00
        do j = 1, n
          total = total + a(i,j)
        end do
        work(i) = a(i,n) - 0.5D+00 * total
      end do

      p = 0.0D+00
      sgn = -1.0D+00

10    continue

        sgn = -sgn
        call sub_gray_next ( n-1, iwork, more, ncard, iadd )

        if ( ncard .ne. 0 ) then
          z = dble ( 2 * iwork(iadd) - 1 )
          do i = 1, n
            work(i) = work(i) + z * a(i,iadd)
          end do
        end if

        prod = 1.0D+00
        do i = 1, n
          prod = prod * work(i)
        end do

        p = p + sgn * prod

        if ( .not. more ) then
          go to 20
        end if

      go to 10

20    continue

      perm = p * dble ( 4 * mod ( n, 2 ) - 2 )

      return
      end
      subroutine r8mat_print ( m, n, a, title )

c*********************************************************************72
c
cc R8MAT_PRINT prints an R8MAT.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    20 May 2004
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
c    Input, double precision A(M,N), the matrix.
c
c    Input, character * ( * ) TITLE, a title to be printed.
c
      implicit none

      integer m
      integer n

      double precision a(m,n)
      character * ( * ) title

      call r8mat_print_some ( m, n, a, 1, 1, m, n, title )

      return
      end
      subroutine r8mat_print_some ( m, n, a, ilo, jlo, ihi, jhi, 
     &  title )

c*********************************************************************72
c
cc R8MAT_PRINT_SOME prints some of an R8MAT.
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
c    Input, integer M, N, the number of rows and columns.
c
c    Input, double precision A(M,N), an M by N matrix to be printed.
c
c    Input, integer ILO, JLO, the first row and column to print.
c
c    Input, integer IHI, JHI, the last row and column to print.
c
c    Input, character * ( * ) TITLE, an optional title.
c
      implicit none

      integer incx
      parameter ( incx = 5 )
      integer m
      integer n

      double precision a(m,n)
      character * ( 14 ) ctemp(incx)
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
      integer title_length

      title_length = len_trim ( title )

      if ( 0 .lt. title_length ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) title(1:title_length)
      end if

      do j2lo = max ( jlo, 1 ), min ( jhi, n ), incx

        j2hi = j2lo + incx - 1
        j2hi = min ( j2hi, n )
        j2hi = min ( j2hi, jhi )

        inc = j2hi + 1 - j2lo

        write ( *, '(a)' ) ' '

        do j = j2lo, j2hi
          j2 = j + 1 - j2lo
          write ( ctemp(j2), '(i7,7x)') j
        end do

        write ( *, '(''  Col   '',5a14)' ) ( ctemp(j), j = 1, inc )
        write ( *, '(a)' ) '  Row'
        write ( *, '(a)' ) ' '

        i2lo = max ( ilo, 1 )
        i2hi = min ( ihi, m )

        do i = i2lo, i2hi

          do j2 = 1, inc

            j = j2lo - 1 + j2

            write ( ctemp(j2), '(g14.6)' ) a(i,j)

          end do

          write ( *, '(i5,1x,5a14)' ) i, ( ctemp(j), j = 1, inc )

        end do

      end do

      write ( *, '(a)' ) ' '

      return
      end
      subroutine r8poly ( n, a, x0, iopt, val )

c*********************************************************************72
c
cc R8POLY performs operations on real polynomials in power or factorial form.
c
c  Discussion:
c
c    The power sum form of a polynomial is
c
c      P(X) = A1 + A2 * X + A3 * X**2 + ... + (AN+1) * X**N
c
c    The Taylor expansion at C has the form
c
c      P(X) = A1 + A2 * (X-C) + A3 * (X-C)**2+... + (AN+1) * (X-C)**N
c
c    The factorial form of a polynomial is
c
c      P(X) = A1 + A2 * X + A3 * (X) * (X-1) + A4 * (X) * (X-1) * (X-2) + ...
c        + (AN+1) * (X) * (X-1) *...* (X-N+1)
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    26 January 2007
c
c  Author:
c
c    Original FORTRAN77 version by Albert Nijenhuis, Herbert Wilf.
c    This FORTRAN77 version by John Burkardt.
c
c  Reference:
c
c    Albert Nijenhuis, Herbert Wilf,
c    Combinatorial Algorithms for Computers and Calculators,
c    Second Edition,
c    Academic Press, 1978,
c    ISBN: 0-12-519260-6,
c    LC: QA164.N54.
c
c  Parameters:
c
c    Input, integer N, the number of coefficients in the polynomial
c    (in other words, the polynomial degree + 1)
c
c    Input/output, double precision A(N), the coefficients of the polynomial.
c    Depending on the option chosen, these coefficients may be overwritten
c    by those of a different form of the polynomial.
c
c    Input, double precision X0, for IOPT = -1, 0, or positive, the value of
c    the argument at which the polynomial is to be evaluated, or the
c    Taylor expansion is to be carried out.
c
c    Input, integer IOPT, a flag describing which algorithm is to
c    be carried out:
c
c    -3: Reverse Stirling.  Input the coefficients of the polynomial in
c    factorial form, output them in power sum form.
c
c    -2: Stirling.  Input the coefficients in power sum
c    form, output them in factorial form.
c
c    -1: Evaluate a polynomial which has been input
c    in factorial form.
c
c    0:  Evaluate a polynomial input in power sum form.
c
c    1 or more:  Given the coefficients of a polynomial in
c    power sum form, compute the first IOPT coefficients of
c    the polynomial in Taylor expansion form.
c
c    Output, double precision VAL, for IOPT = -1 or 0, the value of the
c    polynomial at the point X0.
c
      implicit none

      integer n

      double precision a(n)
      double precision eps
      integer i
      integer iopt
      integer m
      integer n1
      double precision val
      double precision w
      double precision x0
      double precision z

      n1 = min ( n, iopt )
      n1 = max ( 1, n1 )

      if ( iopt .lt. -1 ) then
        n1 = n
      end if

      eps = dble ( mod ( max ( -iopt, 0 ), 2 ) )

      w = - dble ( n ) * eps

      if ( -2 .lt. iopt ) then
        w = w + x0
      end if

      do m = 1, n1

        val = 0.0D+00
        z = w

        do i = m, n
          z = z + eps
          val = a(n+m-i) + z * val
          if ( iopt .ne. 0 .and. iopt .ne. -1 ) then
            a(n+m-i) = val
          end if
        end do

        if ( iopt .lt. 0 ) then
          w = w + 1.0D+00
        end if

      end do

      return
      end
      subroutine r8poly_degree ( na, a, degree )

c*********************************************************************72
c
cc R8POLY_DEGREE returns the degree of a polynomial in power sum form.
c
c  Discussion:
c
c    The power sum form of a polynomial is:
c
c      p(x) = a(0) + a(1) * x + ... + a(n-1) * x**(n-1) + a(n) * x**(n)
c
c    The degree of a polynomial is the index of the highest power
c    of X with a nonzero coefficient.
c
c    The degree of a constant polynomial is 0.  The degree of the
c    zero polynomial is debatable, but this routine returns the
c    degree as 0.
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
c    Input, integer NA, the dimension of A.
c
c    Input, double precision A(0:NA), the coefficients of the polynomials.
c
c    Output, integer DEGREE, the degree of A.
c
      implicit none

      integer na

      double precision a(0:na)
      integer degree

      degree = na

10    continue

      if ( 0 .lt. degree ) then

        if ( a(degree) .ne. 0.0D+00 ) then
          go to 20
        end if

        degree = degree - 1

        go to 10

      end if

20    continue

      return
      end
      subroutine r8poly_div ( na, a, nb, b, nq, q, nr, r )

c*********************************************************************72
c
cc R8POLY_DIV computes the quotient and remainder of two polynomials.
c
c  Discussion:
c
c    The polynomials are assumed to be stored in power sum form.
c
c    The power sum form of a polynomial is:
c
c      p(x) = a(0) + a(1) * x + ... + a(n-1) * x**(n-1) + a(n) * x**(n)
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    26 January 2007
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer NA, the dimension of A.
c
c    Input, double precision A(0:NA), the coefficients of the polynomial
c    to be divided.
c
c    Input, integer NB, the dimension of B.
c
c    Input, double precision B(0:NB), the coefficients of the divisor
c    polynomial.
c
c    Output, integer NQ, the degree of Q.
c    If the divisor polynomial is zero, NQ is returned as -1.
c
c    Output, double precision Q(0:NA-NB), contains the quotient of A/B.
c    If A and B have full degree, Q should be dimensioned Q(0:NA-NB).
c    In any case, Q(0:NA) should be enough.
c
c    Output, integer NR, the degree of R.
c    If the divisor polynomial is zero, NR is returned as -1.
c
c    Output, double precision R(0:NB-1), contains the remainder of A/B.
c    If B has full degree, R should be dimensioned R(0:NB-1).
c    Otherwise, R will actually require less space.
c
      implicit none

      integer na
      integer nb

      double precision a(0:na)
      double precision a2(0:na)
      double precision b(0:nb)
      integer i
      integer j
      integer na2
      integer nb2
      integer nq
      integer nr
      double precision q(0:*)
      double precision r(0:*)

      call r8poly_degree ( na, a, na2 )
      call r8poly_degree ( nb, b, nb2 )

      if ( b(nb2) .eq. 0.0D+00 ) then
        nq = -1
        nr = -1
        return
      end if

      do i = 0, na
        a2(i) = a(i)
      end do

      nq = na2 - nb2
      nr = nb2 - 1

      do i = nq, 0, -1
        q(i) = a2(i+nb2) / b(nb2)
        a2(i+nb2) = 0.0D+00
        do j = 0, nb2-1
          a2(i+j) = a2(i+j) - q(i) * b(j)
        end do
      end do

      do i = 0, nr
        r(i) = a2(i)
      end do

      return
      end
      subroutine r8poly_f2p ( n, a )

c*********************************************************************72
c
cc R8POLY_F2P converts a real polynomial from factorial form to power sum form.
c
c  Discussion:
c
c    The (falling) factorial form is
c
c      p(x) =   a(1)
c             + a(2) * x
c             + a(3) * x*(x-1)
c             ...
c             + a(n) * x*(x-1)*...*(x-(n-2))
c
c    The power sum form is
c
c      p(x) = a(1) + a(2)*x + a(3)*x**2 + ... + a(n)*x**(n-1)
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    26 January 2007
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the dimension of A.
c
c    Input/output, double precision A(N), on input, the polynomial
c    coefficients in factorial form.  On output, the polynomial
c    coefficients in power sum form.
c
      implicit none

      integer n

      double precision a(n)
      integer i
      integer m
      double precision val
      double precision w
      double precision z

      w = - dble ( n )

      do m = 1, n

        val = 0.0D+00
        z = w

        do i = m, n
          z = z + 1.0D+00
          val = a(n+m-i) + z * val
          a(n+m-i) = val
        end do

        w = w + 1.0D+00

      end do

      return
      end
      subroutine r8poly_fval ( n, a, x, val )

c*********************************************************************72
c
cc R8POLY_FVAL evaluates a real polynomial in factorial form.
c
c  Discussion:
c
c    The (falling) factorial form of a polynomial is:
c
c      p(x) = a(1)
c           + a(2)  *x
c           + a(3)  *x*(x-1)
c           +...
c           + a(n-1)*x*(x-1)*(x-2)...*(x-(n-3))
c           + a(n)  *x*(x-1)*(x-2)...*(x-(n-3))*(x-(n-2))
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    26 January 2007
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the dimension of A.
c
c    Input, double precision A(N), the coefficients of the polynomial.
c    A(1) is the constant term.
c
c    Input, double precision X, the point at which the polynomial is
c    to be evaluated.
c
c    Output, double precision VAL, the value of the polynomial at X.
c
      implicit none

      integer n

      double precision a(n)
      integer i
      double precision val
      double precision x

      val = 0.0D+00
      do i = 1, n
        val = a(n+1-i) + ( x - n + i ) * val
      end do

      return
      end
      subroutine r8poly_mul ( na, a, nb, b, c )

c*********************************************************************72
c
cc R8POLY_MUL computes the product of two real polynomials A and B.
c
c  Discussion:
c
c    The polynomials are in power sum form.
c
c    The power sum form is:
c
c      p(x) = a(0) + a(1) * x + ... + a(n-1) * x**(n-1) + a(n) * x**(n)
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    26 January 2007
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer NA, the dimension of A.
c
c    Input, double precision A(0:NA), the coefficients of the first
c    polynomial factor.
c
c    Input, integer NB, the dimension of B.
c
c    Input, double precision B(0:NB), the coefficients of the second
c    polynomial factor.
c
c    Output, double precision C(0:NA+NB), the coefficients of A * B.
c
      implicit none

      integer na
      integer nb
      integer nd
      parameter ( nd = 100 )

      double precision a(0:na)
      double precision b(0:nb)
      double precision c(0:na+nb)
      double precision d(0:nd)
      integer i
      integer j

      if ( nd < na + nb ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'R8POLY_MUL - Fatal error!'
        write ( *, '(a)' ) '  Input NA+NB exceeds internal limit.'
        stop
      end if

      do i = 0, na+nb
        d(i) = 0.0D+00
      end do

      do i = 0, na
        do j = 0, nb
          d(i+j) = d(i+j) + a(i) * b(j)
        end do
      end do

      do i = 0, na+nb
        c(i) = d(i)
      end do

      return
      end
      subroutine r8poly_n2p ( n, a, xarray )

c*********************************************************************72
c
cc R8POLY_N2P converts a real polynomial from Newton form to power sum form.
c
c  Discussion:
c
c    This is done by shifting all the Newton abscissas to zero.
c
c    Actually, what happens is that the abscissas of the Newton form
c    are all shifted to zero, which means that A is the power sum
c    polynomial description and A, XARRAY is the Newton polynomial
c    description.  It is only because all the abscissas are shifted to
c    zero that A can be used as both a power sum and Newton polynomial
c    coefficient array.
c
c    The Newton form of a polynomial is described by an array of N coefficients
c    A and N abscissas X:
c
c      p(x) =   a(1)
c             + a(2) * (x-x(1))
c             + a(3) * (x-x(1)) * (x-x(2))
c             ...
c             + a(n) * (x-x(1)) * (x-x(2)) * ... * (x-x(n-1))
c
c    X(N) does not occur explicitly in the formula for the evaluation of p(x),
c    although it is used in deriving the coefficients A.
c
c    The power sum form of a polynomial is:
c
c      p(x) = a(1) + a(2)*x + ... + a(n-1)*x**(n-2) + a(n)*x**(n-1)
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    26 January 2007
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the dimension of A.
c
c    Input/output, double precision A(N).  On input, the coefficients
c    of the polynomial in Newton form, and on output, the coefficients
c    in power sum form.
c
c    Input/output, double precision XARRAY(N).  On input, the abscissas of
c    the Newton form of the polynomial.  On output, these values
c    have all been set to zero.
c
      implicit none

      integer n

      double precision a(n)
      integer i
      double precision xarray(n)
      double precision zero
      parameter ( zero = 0.0D+00 )

      do i = 1, n
        call r8poly_nx ( n, a, xarray, zero )
      end do

      return
      end
      subroutine r8poly_nval ( n, a, xarray, x, val )

c*********************************************************************72
c
cc R8POLY_NVAL evaluates a real polynomial in Newton form.
c
c  Discussion:
c
c    The Newton form of a polynomial is;
c
c      p(x) = a(1)
c           + a(2)  *(x-x1)
c           + a(3)  *(x-x1)*(x-x2)
c           +...
c           + a(n-1)*(x-x1)*(x-x2)*(x-x3)...*(x-x(n-2))
c           + a(n)  *(x-x1)*(x-x2)*(x-x3)...*(x-x(n-2))*(x-x(n-1))
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    26 January 2007
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the dimension of A.
c
c    Input, double precision A(N), the coefficients of the polynomial.
c    A(1) is the constant term.
c
c    Input, double precision XARRAY(N-1), the N-1 points X which are part
c    of the definition of the polynomial.
c
c    Input, double precision X, the point at which the polynomial 
c    is to be evaluated.
c
c    Output, double precision VAL, the value of the polynomial at X.
c
      implicit none

      integer n

      double precision a(n)
      integer i
      double precision val
      double precision x
      double precision xarray(n-1)

      val = a(n)
      do i = n-1, 1, -1
        val = a(i) + ( x - xarray(i) ) * val
      end do

      return
      end
      subroutine r8poly_nx ( n, a, xarray, x )

c*********************************************************************72
c
cc R8POLY_NX replaces one of the base points in a polynomial in Newton form.
c
c  Discussion:
c
c    The Newton form of a polynomial is described by an array of N coefficients
c    A and N abscissas X:
c
c      p(x) =   a(1)
c             + a(2) * (x-x(1))
c             + a(3) * (x-x(1)) * (x-x(2))
c             ...
c             + a(n) * (x-x(1)) * (x-x(2)) * ... * (x-x(n-1))
c
c    X(N) does not occur explicitly in the formula for the evaluation of p(x),
c    although it is used in deriving the coefficients A.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    26 January 2007
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the dimension of A.
c
c    Input/output, double precision A(N), the polynomial coefficients
c    of the Newton form.
c
c    Input/output, double precision XARRAY(N), the set of abscissas that
c    are part of the Newton form of the polynomial.  On output,
c    the abscissas have been shifted up one index, so that
c    the first location now holds X, and the original contents
c    of XARRAY(N) have been discarded.
c
c    Input, double precision X, the new point to be shifted into XARRAY.
c
      implicit none

      integer n

      double precision a(n)
      integer i
      double precision x
      double precision xarray(n)

      do i = n-1, 1, -1
        a(i) = a(i) + ( x - xarray(i) ) * a(i+1)
      end do

      do i = n, 2, -1
        xarray(i) = xarray(i-1)
      end do

      xarray(1) = x

      return
      end
      subroutine r8poly_p2f ( n, a )

c*********************************************************************72
c
cc R8POLY_P2F converts a real polynomial from power sum form to factorial form.
c
c  Discussion:
c
c    The power sum form is
c
c      p(x) = a(1) + a(2) * x + a(3) * x**2 + ... + a(n) * x**(n-1)
c
c    The (falling) factorial form is
c
c      p(x) =   a(1)
c             + a(2) * x
c             + a(3) * x * (x-1)
c             ...
c             + a(n) * x * (x-1) *...* (x-(n-2))
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    26 January 2007
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the dimension of A.
c
c    Input/output, double precision A(N), on input, the polynomial
c    coefficients in the power sum form, on output, the polynomial
c    coefficients in factorial form.
c
      implicit none

      integer n

      double precision a(n)
      integer i
      integer m
      double precision val

      do m = 1, n
        val = 0.0D+00
        do i = m, n
          val = a(n+m-i) + dble ( m - 1 ) * val
          a(n+m-i) = val
        end do
      end do

      return
      end
      subroutine r8poly_p2n ( n, a, xarray )

c*********************************************************************72
c
cc R8POLY_P2N converts a real polynomial from power sum form to Newton form.
c
c  Discussion:
c
c    This is done by shifting all the Newton abscissas from zero.
c
c    The power sum form of a polynomial is:
c
c      p(x) = a(1) + a(2) * x + ... + a(n-1) * x**(n-2) + a(n) * x**(n-1)
c
c    The Newton form of a polynomial is described by an array of N coefficients
c    A and N abscissas X:
c
c      p(x) =   a(1)
c             + a(2) * (x-x(1))
c             + a(3) * (x-x(1)) * (x-x(2))
c             ...
c             + a(n) * (x-x(1)) * (x-x(2)) * ... * (x-x(n-1))
c
c    X(N) does not occur explicitly in the formula for the evaluation of p(x),
c    although it is used in deriving the coefficients A.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    26 January 2007
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the dimension of A.
c
c    Input/output, double precision A(N).  On input, the coefficients
c    of the polynomial in power sum form, and on output, the
c    coefficients in Newton form.
c
c    Input/output, double precision XARRAY(N).  On input, the desired
c    abscissas of the Newton form of the polynomial.
c
      implicit none

      integer n
      integer n_max
      parameter ( n_max = 100 )

      double precision a(n)
      integer i
      double precision xarray(n)
      double precision work(n_max)

      if ( n_max < n ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'R8POLY_P2N - Fatal error!'
        write ( *, '(a)' ) '  Input N exceeds internal limit.'
        stop
      end if

      do i = 1, n
        work(i) = 0.0D+00
      end do

      do i = n, 1, -1
        call r8poly_nx ( n, a, work, xarray(i) )
      end do

      return
      end
      subroutine r8poly_p2t ( n, a, x )

c*********************************************************************72
c
cc R8POLY_P2T converts a real polynomial from power sum form to Taylor form.
c
c  Discussion:
c
c    The power sum form is
c
c      p(x) = a(1) + a(2)*x + a(3)*x**2 + ... + a(n)*x**(n-1)
c
c    The Taylor form of a polynomial based at X0 is
c
c      p(x) =   a(1)
c             + a(2) * (x-x0)
c             + a(3) * (x-x0)**2
c             ...
c             + a(n) * (x-x0)**(n-1)
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    26 January 2007
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the dimension of A.
c
c    Input/output, double precision A(N), on input, the coefficients in
c    power sum form, and on output, the coefficients in Taylor form.
c
c    Input, double precision X, the point at which the Taylor form of the
c    polynomial is to be based.
c
      implicit none

      integer n

      double precision a(n)
      integer i
      integer m
      double precision val
      double precision x

      do m = 1, n
        val = 0.0D+00
        do i = m, n
          val = a(n+m-i) + x * val
          a(n+m-i) = val
        end do
      end do

      return
      end
      subroutine r8poly_power ( na, a, p, b )

c*********************************************************************72
c
cc R8POLY_POWER computes a positive integer power of a polynomial.
c
c  Discussion:
c
c    The power sum form is:
c
c      p(x) = a(0) + a(1)*x + ... + a(n-1)*x**(n-1) + a(n)*x**(n)
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    20 March 2001
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer NA, the dimension of A.
c
c    Input, double precision A(0:NA), the polynomial to be raised to the power.
c
c    Input, integer P, the nonnegative power to which A is raised.
c
c    Output, double precision B(0:P*NA), the power of the polynomial.
c
      implicit none

      integer na
      integer p

      double precision a(0:na)
      double precision b(0:p*na)
      integer i
      integer j
      integer nonzer
c
c  Zero out B.
c
      do i = 0, p*na
        b(i) = 0.0D+00
      end do
c
c  Search for the first nonzero element in A.
c
      nonzer = 0

      do i = 0, na
        if ( a(i) .ne. 0.0D+00 ) then
          nonzer = i
          go to 10
        end if
      end do

10    continue

      if ( nonzer .eq. 0 ) then
        return
      end if

      b(0) = a(nonzer)**p

      do i = 1, p*(na-nonzer)

        if ( i + nonzer .le. na ) then
          b(i) = dble ( i * p ) * b(0) * a(i+nonzer)
        else
          b(i) = 0.0D+00
        end if

        do j = 1, i-1

          if ( j+nonzer .le. na ) then
            b(i) = b(i) - dble ( i - j ) * a(j+nonzer) * b(i-j)
          end if

          if ( i-j+nonzer .le. na ) then
            b(i) = b(i) + dble ( i - j ) * dble ( p ) 
     &        * b(j) * a(i-j+nonzer)
          end if

        end do

        b(i) = b(i) / ( dble ( i ) * a(nonzer) )

      end do
c
c  Shift B up.
c
      do i = p*nonzer, p*na
        b(i) = b(i-p*nonzer)
      end do

      do i = 0, p * nonzer-1
        b(i) = 0.0D+00
      end do

      return
      end
      subroutine r8poly_print ( n, a, title )

c*********************************************************************72
c
cc R8POLY_PRINT prints out a polynomial.
c
c  Discussion:
c
c    The power sum form is:
c
c      p(x) = a(0) + a(1) * x + ... + a(n-1) * x^(n-1) + a(n) * x^(n)
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    26 January 2007
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the dimension of A.
c
c    Input, double precision A(0:N), the polynomial coefficients.
c    A(0) is the constant term and
c    A(N) is the coefficient of X^N.
c
c    Input, character * ( * ) TITLE, a title.
c
      implicit none

      integer n

      double precision a(0:n)
      integer i
      double precision mag
      integer n2
      character plus_minus
      character * ( * ) title
      integer title_length

      title_length = len_trim ( title )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) title(1:title_length)
      write ( *, '(a)' ) ' '

      call r8poly_degree ( n, a, n2 )

      if ( a(n2) .lt. 0.0D+00 ) then
        plus_minus = '-'
      else
        plus_minus = ' '
      end if

      mag = abs ( a(n2) )

      if ( 2 .le. n2 ) then
        write ( *, '(a,a1,g14.6,a,i3)' ) 
     &    '  p(x) = ', plus_minus, mag, ' * x ^ ', n2
      else if ( n2 .eq. 1 ) then
        write ( *, '(a,a1,g14.6,a)' ) 
     &    '  p(x) = ', plus_minus, mag, ' * x'
      else if ( n2 .eq. 0 ) then
        write ( *, '(a,a1,g14.6)' ) '  p(x) = ', plus_minus, mag
      end if

      do i = n2 - 1, 0, -1

        if ( a(i) .lt. 0.0D+00 ) then
          plus_minus = '-'
        else
          plus_minus = '+'
        end if

        mag = abs ( a(i) )

        if ( mag .ne. 0.0D+00 ) then

          if ( 2 .le. i ) then
            write ( *, ' (9x,a1,g14.6,a,i3)' ) 
     &        plus_minus, mag, ' * x ^ ', i
          else if ( i .eq. 1 ) then
            write ( *, ' (9x,a1,g14.6,a)' ) plus_minus, mag, ' * x'
          else if ( i .eq. 0 ) then
            write ( *, ' (9x,a1,g14.6)' ) plus_minus, mag
          end if
        end if

      end do

      return
      end
      subroutine r8poly_pval ( n, a, x, val )

c*********************************************************************72
c
cc R8POLY_PVAL evaluates a real polynomial in power sum form.
c
c  Discussion:
c
c    The power sum form is:
c
c      p(x) = a(0) + a(1) * x + ... + a(n-1) * x**(n-1) + a(n) * x**(n)
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    26 January 2007
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the dimension of A.
c
c    Input, double precision A(0:N), the coefficients of the polynomial.
c    A(0) is the constant term.
c
c    Input, double precision X, the point at which the polynomial 
c    is to be evaluated.
c
c    Output, double precision VAL, the value of the polynomial at X.
c
      implicit none

      integer n

      double precision a(0:n)
      integer i
      double precision val
      double precision x

      val = 0.0D+00
      do i = n, 0, -1
        val = val * x + a(i)
      end do

      return
      end
      subroutine r8poly_t2p ( n, a, x )

c*********************************************************************72
c
cc R8POLY_T2P converts a real polynomial from Taylor form to power sum form.
c
c  Discussion:
c
c    The Taylor form of a polynomial based at X0 is
c
c      p(x) =   a(1)
c             + a(2) * (x-x0)
c             + a(3) * (x-x0)**2
c             ...
c             + a(n) * (x-x0)**(n-1)
c
c    The power sum form is
c
c      p(x) = a(1) + a(2)*x + a(3)*x**2 + ... + a(n)*x**(n-1)
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    26 January 2007
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the dimension of A.
c
c    Input/output, double precision A(N).  On input, the coefficients 
c    in Taylor form, and on output, the coefficients in power sum form.
c
c    Input, double precision X, the point at which the Taylor form 
c    polynomial is based.
c
      implicit none

      integer n

      double precision a(n)
      integer i
      integer j
      double precision x

      do i = n, 1, -1
        do j = i, n-1
          a(j) = a(j) - a(j+1) * x
        end do
      end do

      return
      end
      subroutine r8vec_backtrack ( n, maxstack, stack, x, indx, k, 
     &  nstack, ncan )

c*********************************************************************72
c
cc R8VEC_BACKTRACK supervises a backtrack search for an R8VEC.
c
c  Discussion:
c
c    The routine tries to construct a real vector one index at a time,
c    using possible candidates as supplied by the user.
c
c    At any time, the partially constructed vector may be discovered to be
c    unsatisfactory, but the routine records information about where the
c    last arbitrary choice was made, so that the search can be
c    carried out efficiently, rather than starting out all over again.
c
c    First, call the routine with INDX = 0 so it can initialize itself.
c
c    Now, on each return from the routine, if INDX is:
c      1, you've just been handed a complete candidate vector;
c         Admire it, analyze it, do what you like.
c      2, please determine suitable candidates for position X(K).
c         Return the number of candidates in NCAN(K), adding each
c         candidate to the end of STACK, and increasing NSTACK.
c      3, you're done.  Stop calling the routine;
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    26 January 2007
c
c  Author:
c
c    Original FORTRAN77 version by Albert Nijenhuis, Herbert Wilf.
c    This FORTRAN77 version by John Burkardt.
c
c  Reference:
c
c    Albert Nijenhuis, Herbert Wilf,
c    Combinatorial Algorithms for Computers and Calculators,
c    Second Edition,
c    Academic Press, 1978,
c    ISBN: 0-12-519260-6,
c    LC: QA164.N54.
c
c  Parameters:
c
c    Input, integer N, the number of positions to be filled in the vector.
c
c    Input, integer MAXSTACK, the maximum length of the stack.
c
c    Input, double precision STACK(MAXSTACK), a list of all current 
c    candidates for all positions 1 through K.
c
c    Input/output, double precision X(N), the partially filled in 
c    candidate vector.
c
c    Input/output, integer INDX, a communication flag.
c    On input,
c      0, to begin a backtracking search.
c      2, the requested candidates for position K have been added to 
c      STACK, and NCAN(K) was updated.
c    On output:
c      1, a complete output vector has been determined and returned in X(1:N);
c      2, candidates are needed for position X(K);
c      3, no more possible vectors exist.
c
c    Input/output, integer K, the index in X that we are trying to fill.
c
c    Input/output, integer NSTACK, the current length of the stack.
c
c    Input/output, integer NCAN(N), lists the current number of candidates for
c    all positions 1 through K.
c
      implicit none

      integer n
      integer maxstack

      integer indx
      integer k
      integer ncan(n)
      integer nstack
      double precision stack(maxstack)
      double precision x(n)
c
c  If this is the first call, request a candidate for position 1.
c
      if ( indx .eq. 0 ) then
        k = 1
        nstack = 0
        indx = 2
        return
      end if
c
c  Examine the stack.
c
10    continue
c
c  If there are candidates for position K, take the first available
c  one off the stack, and increment K.
c
c  This may cause K to reach the desired value of N, in which case
c  we need to signal the user that a complete set of candidates
c  is being returned.
c
        if ( 0 .lt. ncan(k) ) then

          x(k) = stack(nstack)
          nstack = nstack - 1

          ncan(k) = ncan(k) - 1

          if ( k .ne. n ) then
            k = k + 1
            indx = 2
          else
            indx = 1
          end if

          go to 20
c
c  If there are no candidates for position K, then decrement K.
c  If K is still positive, repeat the examination of the stack.
c
        else

          k = k - 1

          if ( k .le. 0 ) then
            indx = 3
            go to 20
          end if

        end if

      go to 10

20    continue

      return
      end
      subroutine r8vec_frac ( n, a, k, afrac )

c*********************************************************************72
c
cc R8VEC_FRAC searches for the K-th smallest entry in an R8VEC.
c
c  Discussion:
c
c    Hoare's algorithm is used.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    26 January 2007
c
c  Parameters:
c
c    Input, integer N, the number of elements of A.
c
c    Input/output, double precision A(N).
c    On input, A is the array to search.
c    On output, the elements of A have been somewhat rearranged.
c
c    Input, integer K, the fractile to be sought.  If K = 1, the minimum
c    entry is sought.  If K = N, the maximum is sought.  Other values
c    of K search for the entry which is K-th in size.  K must be at
c    least 1, and no greater than N.
c
c    Output, double precision AFRAC, the value of the K-th fractile of A.
c
      implicit none

      integer n

      double precision a(n)
      double precision afrac
      integer i
      integer iryt
      integer j
      integer k
      integer left
      double precision x

      if ( n .le. 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'R8VEC_FRAC - Fatal error!'
        write ( *, '(a,i8)' ) '  Illegal nonpositive value of N = ', n
        stop
      end if

      if ( k .le. 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'R8VEC_FRAC - Fatal error!'
        write ( *, '(a,i8)' ) '  Illegal nonpositive value of K = ', k
        stop
      end if

      if ( n .lt. k ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'R8VEC_FRAC - Fatal error!'
        write ( *, '(a,i8)' ) '  Illegal N .lt. K, K = ', k
        stop
      end if

      left = 1
      iryt = n

10    continue

        if ( iryt .le. left ) then
          afrac = a(k)
          go to 60
        end if

        x = a(k)
        i = left
        j = iryt

20      continue

          if ( j .lt. i ) then
            if ( j .lt. k ) then
              left = i
            end if
            if ( k .lt. i ) then
              iryt = j
            end if
            go to 50
          end if
c
c  Find I so that X .le. A(I)
c
30        continue

          if ( a(i) .lt. x ) then
            i = i + 1
            go to 30
          end if
c
c  Find J so that A(J) .le. X
c
40        continue

          if ( x .lt. a(j) ) then
            j = j - 1
            go to 40
          end if

          if ( i .le. j ) then
            call r8_swap ( a(i), a(j) )
            i = i + 1
            j = j - 1
          end if

        go to 20

50      continue

      go to 10

60    continue

      return
      end
      subroutine r8vec_indicator ( n, a )

c*********************************************************************72
c
cc R8VEC_INDICATOR sets an R8VEC to the indicator vector.
c
c  Discussion:
c
c    An R8VEC is an array of double precision real values.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    22 January 2007
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the number of elements of A.
c
c    Output, double precision A(N), the array to be initialized.
c
      implicit none

      integer n

      double precision a(n)
      integer i

      do i = 1, n
        a(i) = dble ( i )
      end do

      return
      end
      subroutine r8vec_mirror_next ( n, a, done )

c*********************************************************************72
c
cc R8VEC_MIRROR_NEXT steps through all sign variations of an R8VEC.
c
c  Discussion:
c
c    In normal use, the user would set every element of A to be positive.
c    The routine will take the input value of A, and output a copy in
c    which the signs of one or more entries have been changed.  Repeatedly
c    calling the routine with the output from the previous call will generate
c    every distinct "variation" of A; that is, all possible sign variations.
c
c    When the output variable DONE is TRUE (or equal to 1), then the
c    output value of A_NEW is the last in the series.
c
c    Note that A may have some zero values.  The routine will essentially
c    ignore such entries; more exactly, it will not stupidly assume that -0
c    is a proper "variation" of 0c
c
c    Also, it is possible to call this routine with the signs of A set
c    in any way you like.  The routine will operate properly, but it
c    will nonethess terminate when it reaches the value of A in which
c    every nonzero entry has negative sign.
c
c
c    More efficient algorithms using the Gray code seem to require internal
c    memory in the routine, which is not one of MATLAB's strong points,
c    or the passing back and forth of a "memory array", or the use of
c    global variables, or unnatural demands on the user.  This form of
c    the routine is about as clean as I can make it.
c
c  Example:
c
c      Input         Output
c    ---------    --------------
c    A            A_NEW     DONE
c    ---------    --------  ----
c     1  2  3     -1  2  3  false
c    -1  2  3      1 -2  3  false
c     1 -2  3     -1 -2  3  false
c    -1 -2  3      1  2 -3  false
c     1  2 -3     -1  2 -3  false
c    -1  2 -3      1 -2 -3  false
c     1 -2 -3     -1 -2 -3  false
c    -1 -2 -3      1  2  3  true
c
c     1  0  3     -1  0  3  false
c    -1  0  3      1  0 -3  false
c     1  0 -3     -1  0 -3  false
c    -1  0 -3      1  0  3  true
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
c  Reference:
c
c    Albert Nijenhuis, Herbert Wilf,
c    Combinatorial Algorithms for Computers and Calculators,
c    Second Edition,
c    Academic Press, 1978,
c    ISBN: 0-12-519260-6,
c    LC: QA164.N54.
c
c  Parameters:
c
c    Input, integer N, the number of entries in the vector.
c
c    Input/output, double precision A(N), a vector of real numbers.  
c    On output, some signs have been changed.
c
c    Output, logical DONE, is TRUE if the input vector A was the last element
c    in the series (every entry was nonpositive); the output vector is reset 
c    so that all entries are nonnegative, but presumably the ride is overc
c
      implicit none

      integer n

      double precision a(n)
      logical done
      integer i
      integer positive
c
c  Seek the first strictly positive entry of A.
c
      positive = 0
      do i = 1, n
        if ( 0.0D+00 .lt. a(i) ) then
          positive = i
          go to 10
        end if
      end do

10    continue
c
c  If there is no strictly positive entry of A, there is no successor.
c
      if ( positive .eq. 0 ) then
        do i = 1, n
          a(i) = -a(i)
        end do
        done = .true.
        return
      end if
c
c  Otherwise, negate A up to the positive entry.
c
      do i = 1, positive
        a(i) = -a(i)
      end do
      done = .false.

      return
      end
      subroutine r8vec_print ( n, a, title )

c*********************************************************************72
c
cc R8VEC_PRINT prints an R8VEC.
c
c  Discussion:
c
c    An R8VEC is an array of double precision real values.
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
c    Input, character * ( * ) TITLE, an optional title.
c
      implicit none

      integer n

      double precision a(n)
      integer i
      integer s_len_trim
      character * ( * ) title
      integer title_length

      title_length = s_len_trim ( title )
      if ( 0 .lt. title_length ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) title(1:title_length)
      end if

      write ( *, '(a)' ) ' '
      do i = 1, n
        write ( *, '(2x,i8,2x,g16.8)' ) i, a(i)
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
      integer i4_huge
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
          seed = seed + i4_huge ( )
        end if

        r(i) = a + ( b - a ) * dble ( seed ) * 4.656612875D-10

      end do

      return
      end
      subroutine r8vec_uniform_01 ( n, seed, r )

c*********************************************************************72
c
cc R8VEC_UNIFORM_01 returns a unit pseudorandom R8VEC.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    19 August 2004
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
c    Input, integer N, the number of entries in the vector.
c
c    Input/output, integer SEED, the "seed" value, which should NOT be 0.
c    On output, SEED has been updated.
c
c    Output, double precision R(N), the vector of pseudorandom values.
c
      implicit none

      integer n

      integer i
      integer i4_huge
      integer k
      integer seed
      double precision r(n)

      if ( seed .eq. 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'R8VEC_UNIFORM_01 - Fatal error!'
        write ( *, '(a)' ) '  Input value of SEED = 0.'
        stop
      end if

      do i = 1, n

        k = seed / 127773

        seed = 16807 * ( seed - k * 127773 ) - k * 2836

        if ( seed .lt. 0 ) then
          seed = seed + i4_huge ( )
        end if

        r(i) = dble ( seed ) * 4.656612875D-10

      end do

      return
      end
      subroutine rat_add ( itop1, ibot1, itop2, ibot2, itop, ibot, 
     &  ierror )

c*********************************************************************72
c
cc RAT_ADD adds two rational values.
c
c  Discussion:
c
c    The routine computes
c
c      ITOP/IBOT = ITOP1/IBOT1 + ITOP2/IBOT2
c
c    while trying to avoid integer overflow.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    26 January 2007
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer ITOP1, IBOT1, the first value to add.
c
c    Input, integer ITOP2, IBOT2, the second value to add.
c
c    Output, integer ITOP, IBOT, the sum.
c
c    Output, integer IERROR.
c    0, no error occurred.
c    1, an error occurred.  The addition of the two values
c    requires a numerator or denominator larger than the
c    maximum legal integer.
c
      implicit none

      integer i4_gcd
      integer ibot
      integer ibot1
      integer ibot2
      integer ierror
      integer i_max
      integer i4_huge
      integer itemp
      integer itop
      integer itop1
      integer itop2
      integer jbot1
      integer jbot2
      integer jbot3
      integer jtop1
      integer jtop2

      i_max = i4_huge ( )

      ierror = 0

      if ( itop1 .eq. 0 ) then
        itop = itop2
        ibot = ibot2
        return
      else if ( itop2 .eq. 0 ) then
        itop = itop1
        ibot = ibot1
        return
      end if
c
c  Make copies of the input arguments, since we will change them.
c
      jbot1 = ibot1
      jbot2 = ibot2
      jtop1 = itop1
      jtop2 = itop2
c
c  Compute the greatest common factor of the two denominators,
c  and factor it out.
c
      jbot3 = i4_gcd ( jbot1, jbot2 )
      jbot1 = jbot1 / jbot3
      jbot2 = jbot2 / jbot3
c
c  The fraction may now be formally written as:
c
c    (jtop1*jbot2 + jtop2*jbot1) / (jbot1*jbot2*jbot3)
c
c  Check the tops for overflow.
c
      if ( dble ( i_max ) 
     &  .lt. abs ( dble ( jtop1 ) * dble ( jbot2 ) ) ) then
        ierror = 1
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'RAT_ADD - Fatal error!'
        write ( *, '(a)' ) '  Overflow of top of rational sum.'
        itop = 0
        stop
      end if

      jtop1 = jtop1 * jbot2

      if ( dble ( i_max ) 
     &  .lt. abs ( dble ( jtop2 ) * dble ( jbot1 ) ) ) then
        ierror = 1
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'RAT_ADD - Fatal error!'
        write ( *, '(a)' ) '  Overflow of top of rational sum.'
        itop = 0
        stop
      end if

      jtop2 = jtop2 * jbot1

      if ( dble ( i_max ) 
     &  .lt. abs ( dble ( jtop1 ) + dble ( jtop2 ) ) ) then
        ierror = 1
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'RAT_ADD - Fatal error!'
        write ( *, '(a)' ) '  Overflow of top of rational sum.'
        itop = 0
        stop
      end if

      itop = jtop1 + jtop2
c
c  Check the bottom for overflow.
c
      if ( dble ( i_max ) .lt. abs ( dble ( jbot1 ) 
     &  * dble ( jbot2 ) * dble ( jbot3 ) ) ) then
        ierror = 1
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'RAT_ADD - Fatal error!'
        write ( *, '(a)' ) '  Overflow of bottom of rational sum.'
        ibot = 1
        stop
      end if

      ibot = jbot1 * jbot2 * jbot3
c
c  Put the fraction in lowest terms.
c
      itemp = i4_gcd ( itop, ibot )
      itop = itop / itemp
      ibot = ibot / itemp
c
c  The bottom should be positive.
c
      if ( ibot .lt. 0 ) then
        ibot = -ibot
        itop = -itop
      end if

      return
      end
      subroutine rat_div ( itop1, ibot1, itop2, ibot2, itop, ibot, 
     &  ierror )

c*********************************************************************72
c
cc RAT_DIV divides one rational value by another.
c
c  Discussion:
c
c    The routine computes
c
c      ITOP / IBOT = ( ITOP1 / IBOT1 ) / ( ITOP2 / IBOT2 ).
c
c    while avoiding integer overflow.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    26 January 2007
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer ITOP1, IBOT1, the numerator.
c
c    Input, integer ITOP2, IBOT2, the denominator.
c
c    Output, integer ITOP, IBOT, the result.
c
c    Output, integer IERROR.
c    0, no error occurred.
c    1, an error occurred.  One of the quantities IBOT1, IBOT2,
c    or ITOP2 is zero, or the result of the division
c    requires a numerator or denominator larger than the
c    maximum legal integer.
c
      implicit none

      integer i_max
      integer i4_gcd
      integer i4_huge
      integer ibot
      integer ibot1
      integer ibot2
      integer ierror
      integer itemp
      integer itop
      integer itop1
      integer itop2
      integer jbot1
      integer jbot2
      integer jtop1
      integer jtop2

      ierror = 0

      i_max = i4_huge ( )

      if ( ibot1 .eq. 0 .or. itop2 .eq. 0 .or. ibot2 .eq. 0 ) then
        ierror = 1
        return
      end if

      if ( itop1 .eq. 0 ) then
        itop = 0
        ibot = 1
        return
      end if
c
c  Make copies of the input arguments, since we will change them.
c  Implicitly invert the divisor fraction here.  The rest of
c  the code will be a multiply operation.
c
      jbot1 = ibot1
      jbot2 = itop2
      jtop1 = itop1
      jtop2 = ibot2
c
c  Get rid of all common factors in top and bottom.
c
      itemp = i4_gcd ( jtop1, jbot1 )
      jtop1 = jtop1 / itemp
      jbot1 = jbot1 / itemp
      itemp = i4_gcd ( jtop1, jbot2 )
      jtop1 = jtop1 / itemp
      jbot2 = jbot2 / itemp
      itemp = i4_gcd ( jtop2, jbot1 )
      jtop2 = jtop2 / itemp
      jbot1 = jbot1 / itemp
      itemp = i4_gcd ( jtop2, jbot2 )
      jtop2 = jtop2 / itemp
      jbot2 = jbot2 / itemp
c
c  The fraction (ITOP1*IBOT2)/(IBOT1*ITOP2) is in lowest terms.
c
c  Check the top for overflow.
c
      if ( dble ( i_max ) 
     &  .lt. abs ( dble ( jtop1 ) * dble ( jtop2 ) ) ) then
        ierror = 1
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'RAT_DIV - Fatal error!'
        write ( *, '(a)' ) '  Overflow of top of rational fraction.'
        itop = 0
        stop
      end if

      itop = jtop1 * jtop2
c
c  Check the bottom IBOT1*ITOP2 for overflow.
c
      if ( dble ( i_max ) 
     &  .lt. abs ( dble ( jbot1 ) * dble ( jbot2 ) ) ) then
        ierror = 1
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'RAT_DIV - Fatal error!'
        write ( *, '(a)' ) '  Overflow of bottom of rational fraction.'
        ibot = 1
        stop
      end if

      ibot = jbot1 * jbot2
c
c  The bottom should be positive.
c
      if ( ibot .lt. 0 ) then
        ibot = -ibot
        itop = -itop
      end if
c
c  The fraction is ITOP/IBOT with no loss of accuracy.
c
      return
      end
      subroutine rat_farey ( n, max_frac, num_frac, a, b )

c*********************************************************************72
c
cc RAT_FAREY computes the N-th row of the Farey fraction table.
c
c  Example:
c
c    N = 5
c
c    NUM_FRAC = 11
c    A =  0  1  1  1  2  1  3  2  3  4  1
c    B =  1  5  4  3  5  2  5  3  4  5  1
c
c  Discussion:
c
c    In this form of the Farey fraction table, fractions in row N lie between
c    0 and 1, are in lowest terms, and have a denominator that is no greater
c    than N.  Row N is computed directly, and does not require the computation
c    of previous rows.
c
c    The data satisfy the relationship:
c
c      A(K+1) * B(K) - A(K) * B(K+1) = 1
c
c    The number of items in the N-th row is roughly N**2 / PI**2.
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
c  Reference:
c
c    Donald Knuth,
c    The Art of Computer Programming,
c    Volume 1, Fundamental Algorithms,
c    Addison Wesley, 1968, page 157.
c
c  Parameters:
c
c    Input, integer N, the desired row number.  N must be positive.
c
c    Input, integer MAX_FRAC, the maximum number of fractions to compute.
c
c    Output, integer NUM_FRAC, the number of fractions computed.
c
c    Output, integer A(MAX_FRAC), B(MAX_FRAC), contains the NUM_FRAC
c    numerators and denominators of the N-th row of the Farey fraction table.
c
      implicit none

      integer max_frac

      integer a(max_frac)
      integer b(max_frac)
      integer c
      integer k
      integer n
      integer num_frac

      if ( n .le. 0 ) then
        num_frac = 0
        return
      end if

      if ( max_frac .le. 0 ) then
        num_frac = 0
        return
      end if

      k = 1
      a(k) = 0
      b(k) = 1

      if ( max_frac .le. 1 ) then
        num_frac = k
        return
      end if

      k = 2
      a(k) = 1
      b(k) = n

10    continue

      if ( k .lt. max_frac ) then

        if ( a(k) .eq. 1 .and. b(k) .eq. 1 ) then
          go to 20
        end if

        k = k + 1
        c = ( b(k-2) + n ) / b(k-1)
        a(k) = c * a(k-1) - a(k-2)
        b(k) = c * b(k-1) - b(k-2)

        go to 10

      end if

20    continue

      num_frac = k

      return
      end
      subroutine rat_farey2 ( n, a, b )

c*********************************************************************72
c
cc RAT_FAREY2 computes the next row of the Farey fraction table.
c
c  Example:
c
c    Input:
c
c      N = 3
c      A =  0  1  1  2  1
c      B =  1  3  2  3  1
c
c    Output:
c
c      A =  0  1  1  2  1  3  2  3  1
c      B =  1  4  3  5  2  5  3  4  1
c
c  Discussion:
c
c    In this form of the Farey fraction table, fractions in row N lie between
c    0 and 1, and are in lowest terms.  For every adjacent pair of input
c    fractions, A1/B1 and A2/B2, the mediant (A1+A2)/(B1+B2) is computed
c    and inserted between them.
c
c    The number of items in the N-th row is 1+2**(N-1).
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    26 January 2007
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the input row number.  N must be nonnegative.
c    If N is zero, then the input is ignored, and the entries of
c    row 1 are computed directly.
c
c    Input/output, integer A(1+2**N), B(1+2**N).
c    On input, entries 1 through 1+2**(N-1) contain the entries of row N.
c    On output, entries 1 through 1+2**N contain the entries of row N+1.
c
      implicit none

      integer n

      integer a(1+2**n)
      integer b(1+2**n)
      integer i

      if ( n .eq. 0 ) then
        a(1) = 0
        b(1) = 1
        a(2) = 1
        b(2) = 1
        return
      end if
c
c  Shift the current data.
c
      do i = 1+2**(n-1), 1, -1
        a(2*i-1) = a(i)
        b(2*i-1) = b(i)
      end do
c
c  Compute the mediants.
c
      do i = 2, 2**n, 2
        a(i) = a(i-1) + a(i+1)
        b(i) = b(i-1) + b(i+1)
      end do

      return
      end
      subroutine rat_mul ( itop1, ibot1, itop2, ibot2, itop, ibot, 
     &  ierror )

c*********************************************************************72
c
cc RAT_MUL multiplies two fractions.
c
c  Discussion:
c
c    The routine computes
c
c      ITOP / IBOT = ( ITOP1 / IBOT1 ) * ( ITOP2 / IBOT2 ).
c
c    while avoiding integer overflow.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    26 January 2007
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer ITOP1, IBOT1, the first factor.
c
c    Input, integer ITOP2, IBOT2, the second factor.
c
c    Output, integer ITOP, IBOT, the product.
c
c    Output, integer IERROR.
c    0, no error occurred.
c    1, an error occurred.  The multiplication of the two values
c    requires a numerator or denominator larger than the
c    maximum legal integer.
c
      implicit none

      integer i_max
      integer i4_gcd
      integer i4_huge
      integer ibot
      integer ibot1
      integer ibot2
      integer ierror
      integer itemp
      integer itop
      integer itop1
      integer itop2
      integer jbot1
      integer jbot2
      integer jtop1
      integer jtop2

      ierror = 0

      i_max = i4_huge ( )

      if ( itop1 .eq. 0 .or. itop2 .eq. 0 ) then
        itop = 0
        ibot = 1
        return
      end if
c
c  Make copies of the input arguments, since we will change them.
c
      jbot1 = ibot1
      jbot2 = ibot2
      jtop1 = itop1
      jtop2 = itop2
c
c  Get rid of all common factors in top and bottom.
c
      itemp = i4_gcd ( jtop1, jbot1 )
      jtop1 = jtop1 / itemp
      jbot1 = jbot1 / itemp
      itemp = i4_gcd ( jtop1, jbot2 )
      jtop1 = jtop1 / itemp
      jbot2 = jbot2 / itemp
      itemp = i4_gcd ( jtop2, jbot1 )
      jtop2 = jtop2 / itemp
      jbot1 = jbot1 / itemp
      itemp = i4_gcd ( jtop2, jbot2 )
      jtop2 = jtop2 / itemp
      jbot2 = jbot2 / itemp
c
c  The fraction (ITOP1*ITOP2)/(IBOT1*IBOT2) is in lowest terms.
c
c  Check the top ITOP1*ITOP2 for overflow.
c
      if ( dble ( i_max ) 
     &  .lt. abs ( dble ( jtop1 ) * dble ( jtop2 ) ) ) then
        ierror = 1
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'RAT_MUL - Fatal error!'
        write ( *, '(a)' ) '  Overflow of top of rational product.'
        itop = 0
        stop
      end if

      itop = jtop1 * jtop2
c
c  Check the bottom IBOT1*IBOT2 for overflow.
c
      if ( dble ( i_max ) 
     &  .lt. abs ( dble ( jbot1 ) * dble ( jbot2 ) ) ) then
        ierror = 1
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'RAT_MUL - Fatal error!'
        write ( *, '(a)' ) '  Overflow of bottom of rational product.'
        ibot = 1
        stop
      end if

      ibot = jbot1 * jbot2
c
c  The bottom should be positive.
c
      if ( ibot .lt. 0 ) then
        ibot = -ibot
        itop = -itop
      end if
c
c  The fraction is ITOP/IBOT with no loss of accuracy.
c
      return
      end
      subroutine rat_normalize ( a, b )

c*********************************************************************72
c
cc RAT_NORMALIZE normalizes a rational number.
c
c  Discussion:
c
c    If A = B = 0, return.
c
c    If A = 0 (and B nonzero) set B => 1 and return.
c
c    If A nonzero, and B = 0, then A => 1 and return.
c
c    If A = B, then set A => 1, B => 1 and return.
c
c    If B .lt. 0, then A => -A, B => -B.
c
c    If 1 .lt. C = GCD(|A|,|B|), A => A/C, B => B/C.
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
c    Input/output, integer A, B, the rational number.
c
      implicit none

      integer a
      integer b
      integer c
      integer i4_gcd
c
c  Cases where one or both is 0.
c
      if ( a .eq. 0 .and. b .eq. 0 ) then
        return
      else if ( a .eq. 0 .and. b .ne. 0 ) then
        b = 1
        return
      else if ( a .ne. 0 .and. b .eq. 0 ) then
        a = 1
        return
      end if

      if ( a .eq. b ) then
        a = 1
        b = 1
        return
      end if
        
      if ( b .lt. 0 ) then
        a = -a
        b = -b
      end if

      c = i4_gcd ( abs ( a ), abs ( b ) )

      if ( 1 .lt. c ) then
        a = a / c
        b = b / c
      end if

      return
      end
      subroutine rat_sum_formula ( n, a, b )

c*********************************************************************72
c
cc RAT_SUM_FORMULA computes the formulas for sums of powers of integers.
c
c  Example:
c
c    N = 6
c
c        1    2    3    4    5    6    7
c    -----------------------------------
c    0 | 1    0    0    0    0    0    0
c      |
c    1 | 1    1    0    0    0    0    0
c      | 2    2
c      |
c    2 | 1    1    1    0    0    0    0
c      | 3    2    6
c      |
c    3 | 1    1    1    0    0    0    0
c      | 4    2    4
c      | 
c    4 | 1    1    1    0   -1    0    0
c      | 5    2    3        30
c      |
c    5 | 1    1    5    0   -1    0    0
c      | 6    2   12        12
c      |
c    6 | 1    1    1    0   -1    0    1
c      | 7    2    2         6        42
c
c    The interpretation of row 2, for instance, is:
c
c      sum ( 1 .le. I .le. N ) I**2 = 1/3 N**3 + 1/2 N**2 + 1/6 N
c
c    This suggests that a more sensible way to display the table
c    is to reverse the order of the entries in the row, so that
c    the entry in column J is the coeficient of N**J, which is
c    not the case now.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    26 January 2007
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Robert Owens,
c    Sums of Powers of Integers,
c    Mathematics Magazine,
c    Volume 65, Number 1, February 1992, pages 38-40.
c
c  Parameters:
c
c    Input, integer N, the number of rows of coefficients to compute.
c
c    Output, integer A(0:N,N+1), B(0:N,N+1), the numerator and denominator
c    of the coefficients.
c
      implicit none

      integer n

      integer a(0:n,1:n+1)
      integer asum
      integer b(0:n,1:n+1)
      integer bsum
      integer i
      integer ierror
      integer j

      a(0,1) = 1
      do j = 2, n+1
        a(0,j) = 0
      end do

      b(0,1) = 1
      do j = 2, n+1
        b(0,j) = 0
      end do

      do i = 1, n

        asum = 0
        bsum = 0
c
c  Subdiagonal entries are multiples of entries above them.
c
        do j = 1, i

          call rat_mul ( a(i-1,j), b(i-1,j), i, i+2-j, a(i,j), b(i,j), 
     &      ierror )

          call rat_add ( asum, bsum, a(i,j), b(i,j), asum, bsum, 
     &      ierror )

        end do
c
c  Diagonal entry is 1 - sum of previous entries in row.
c
        asum = -asum
        call rat_add ( 1, 1, asum, bsum, a(i,i+1), b(i,i+1), ierror )
c
c  Superdiagonal entries are zero.
c
        do j = i+2, n+1
          a(i,j) = 0
        end do

        do j = i+2, n+1
          b(i,j) = 1
        end do

      end do

      return
      end
      subroutine rat_to_cfrac ( ip, iq, m, n, a, ierror )

c*********************************************************************72
c
cc RAT_TO_CFRAC converts a rational value to a continued fraction.
c
c  Discussion:
c
c    The routine is given a rational number represented by IP/IQ, and
c    computes the monic or "simple" continued fraction representation
c    with integer coefficients of the number:
c
c      A(1) + 1/ (A(2) + 1/ (A(3) + ... + 1/A(N) ...))
c
c    The user must dimension A to a value M which is "large enough".
c    The actual number of terms needed in the continued fraction
c    representation cannot be known beforehand.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    09 January 2007
c
c  Author:
c
c    Original FORTRAN77 version by John Hart, Ward Cheney, Charles Lawson, 
c    Hans Maehly, Charles Mesztenyi, John Rice, Henry Thatcher, 
c    Christoph Witzgall.
c    This FORTRAN77 version by John Burkardt.
c
c  Reference:
c
c    John Hart, Ward Cheney, Charles Lawson, Hans Maehly, Charles Mesztenyi,
c    John Rice, Henry Thatcher, Christoph Witzgall,
c    Computer Approximations,
c    Wiley, 1968.
c
c  Parameters:
c
c    Input, integer IP, IQ, the numerator and denominator of the
c    rational value whose continued fraction representation is
c    desired.
c
c    Input, integer M, the dimension of A.  If M is not great
c    enough, the algorithm may run out of space.
c
c    Output, integer N, the actual number of entries used in A.
c
c    Output, integer A(M), contains the continued fraction
c    representation of the number.
c
c    Output, integer IERROR, error indicator.  0 if no error,
c    1 if there was an error, namely, M is not large enough.
c
      implicit none

      integer m

      integer a(m)
      integer ierror
      integer ip
      integer iq
      integer jp
      integer jq
      integer n

      jp = ip
      jq = iq

      n = 0

10    continue

        n = n + 1

        if ( m .lt. n ) then
          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) 'RAT_TO_CFRAC - Fatal error!'
          write ( *, '(a)' ) '  M < N.'
          write ( *, '(a)' ) '  M = ', m
          write ( *, '(a)' ) '  N = ', n
          ierror = 1
          stop
        end if

        a(n) = jp / jq
        jp = mod ( jp, jq )

        if ( jp .eq. 0 ) then
          return
        end if

        n = n + 1

        if ( m .lt. n ) then
          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) 'RAT_TO_CFRAC - Fatal error!'
          write ( *, '(a)' ) '  M < N.'
          write ( *, '(a)' ) '  M = ', m
          write ( *, '(a)' ) '  N = ', n
          ierror = 1
          stop
        end if

        a(n) = jq / jp
        jq = mod ( jq, jp )

        if ( jq .eq. 0 ) then
          go to 20
        end if

      go to 10

20    continue

      return
      end
      subroutine rat_to_dec ( rat_top, rat_bot, mantissa, exponent )

c*********************************************************************72
c
cc RAT_TO_DEC converts a rational to a decimal representation.
c
c  Discussion:
c
c    A rational value is represented by RAT_TOP / RAT_BOT.
c
c    A decimal value is represented as MANTISSA * 10**EXPONENT.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    22 November 2004
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer RAT_TOP, RAT_BOT, the rational value.
c
c    Output, integer MANTISSA, EXPONENT, the decimal number.
c
      implicit none

      integer exponent
      integer gcd
      integer i4_gcd
      integer i4_huge
      integer mantissa
      double precision r
      double precision r_max
      integer rat_bot
      integer rat_bot2
      integer rat_top
      integer rat_top2
      integer s

      if ( rat_top .eq. 0 ) then
        mantissa = 0
        exponent = 0
        return
      end if

      gcd = i4_gcd ( rat_top, rat_bot )
      rat_top2 = rat_top / gcd
      rat_bot2 = rat_bot / gcd

      if ( rat_bot2 .lt. 0 ) then
        rat_top2 = -rat_top2
        rat_bot2 = -rat_bot2
      end if

      if ( rat_bot2 .eq. 1 ) then
        mantissa = rat_top2
        exponent = 0
        return
      end if

      exponent = 0

10    continue

      if ( mod ( rat_bot2, 10 ) .eq. 0 ) then
        exponent = exponent - 1
        rat_bot2 = rat_bot2 / 10
        go to 10
      end if

20    continue

      if ( mod ( rat_top2, 10 ) .eq. 0 ) then
        exponent = exponent + 1
        rat_top2 = rat_top2 / 10
        go to 20
      end if

      r = dble ( rat_top2 ) / dble ( rat_bot2 )

      if ( r .lt. 0.0D+00 ) then
        s = -1
        r = -r
      else
        s = 1
      end if

      r_max = dble ( i4_huge ( ) ) / 10.0D+00

30    continue

      if ( r .ne. dble ( int ( r ) ) .and. r .lt. r_max ) then
        r = r * 10.0D+00
        exponent = exponent - 1
        go to 30
      end if

      mantissa = s * int ( r )

      return
      end
      subroutine rat_to_r8 ( a, b, r )

c*********************************************************************72
c
cc RAT_TO_R8 converts rational values to real values.
c
c  Example:
c
c    A    B    R
c   --   --    ---
c    1    2    0.5
c    7    5    1.4
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
c    Input, integer A, B, the rational quantity to be converted.
c
c    Output, double precision R, the value of the rational quantity 
c    as a real number.
c
      implicit none

      integer a
      integer b
      double precision r

      if ( b .eq. 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'RAT_TO_R8 - Warning!'
        write ( *, '(a)' ) '  The input fraction to be converted had a'
        write ( *, '(a)' ) '  zero denominator.'
        r = 0.0D+00
      else
        r = dble ( a ) / dble ( b )
      end if

      return
      end
      subroutine rat_to_s_left ( a, b, s )

c*********************************************************************72
c
cc RAT_TO_S_LEFT returns a left-justified representation of A/B.
c
c  Discussion:
c
c    If the ratio is negative, a minus sign precedes A.
c    A slash separates A and B.
c
c    Note that if A is nonzero and B is 0, S will
c    be returned as "Inf" or "-Inf" (Infinity), and if both
c    A and B are zero, S will be returned as "NaN"
c    (Not-a-Number).
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    03 June 2003
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer A, B, the numerator and denominator.
c
c    Output, character * ( * ) S, a left-justified string
c    containing the representation of A/B.
c
      implicit none

      integer a
      integer b
      character * ( * ) s
      character * ( 25 ) s2
c
c  Take care of simple cases right away.
c
      if ( a .eq. 0 ) then

        if ( b .ne. 0 ) then
          s2 = '0'
        else
          s2= 'NaN'
        end if

      else if ( b .eq. 0 ) then

        if ( 0 .lt. a ) then
          s2 = 'Inf'
        else
          s2 = '-Inf'
        end if
c
c  Make copies of A and B.
c
      else

        if ( b .eq. 1 ) then
          write ( s2, '(i12)' ) a
        else
          write ( s2, '(i12, ''/'', i12)' ) a, b
        end if

        call s_blank_delete ( s2 )

      end if

      s = s2

      return
      end
      function rat_width ( a, b )

c*********************************************************************72
c
cc RAT_WIDTH returns the "width" of a rational number.
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
c    Input, integer A, B, the rational number.
c
c    Output, integer RAT_WIDTH, the "width" of the rational number.
c
      implicit none

      integer a
      integer abs_a
      integer abs_b
      integer b
      integer rat_width
      integer ten_pow
      integer value

      value = 1
      ten_pow = 10

      if ( a .eq. 0 ) then
        rat_width = 1
        return
      end if
      
      abs_a = abs ( a )

10    continue

      if ( ten_pow .le. abs_a ) then
        value = value + 1
        ten_pow = ten_pow * 10
        go to 10
      end if
c
c  If the fraction is negative, a minus sign will be prepended to the
c  numerator.
c
      if ( a * b .lt. 0 ) then
        value = value + 1
        ten_pow = ten_pow * 10
      end if

      abs_b = abs ( b )

20    continue

      if ( ten_pow .le. abs_b ) then
        value = value + 1
        ten_pow = ten_pow * 10
        go to 20
      end if

      rat_width = value

      return
      end
      subroutine ratmat_det ( n, iatop, iabot, idtop, idbot, ierror )

c*********************************************************************72
c
cc RATMAT_DET finds the determinant of an N by N matrix of rational entries.
c
c  Discussion:
c
c    The brute force method is used.
c
c    This routine should only be used for small matrices, since this
c    calculation requires the summation of Nc products of N numbers.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    26 January 2007
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the number of rows and columns of A.
c
c    Input, integer IATOP(N,N), IABOT(N,N), the numerators
c    and denominators of the entries of the matrix.
c
c    Output, integer IDTOP, IDBOT, the determinant of the matrix,
c    expressed as IDTOP/IDBOT.
c
c    Output, integer IERROR.
c    0, the determinant was computed.
c    1, an overflow error occurred, and the determinant was not
c    computed.
c
      implicit none

      integer n

      logical even
      integer i
      integer iabot(n,n)
      integer iatop(n,n)
      integer iarray(n)
      integer ibot
      integer ibot1
      integer ibot2
      integer idbot
      integer idtop
      integer ierror
      integer itop
      integer itop1
      integer itop2
      logical more

      ierror = 0

      more = .false.
      idtop = 0
      idbot = 1

10    continue

        call perm_next ( n, iarray, more, even )

        if ( even ) then
          itop = 1
        else
          itop = -1
        end if

        ibot = 1

        do i = 1, n

          itop1 = itop
          ibot1 = ibot
          itop2 = iatop(i,iarray(i))
          ibot2 = iabot(i,iarray(i))

          call rat_mul ( itop1, ibot1, itop2, ibot2, itop, ibot, 
     &      ierror )

          if ( ierror .ne. 0 ) then
            write ( *, '(a)' ) ' '
            write ( *, '(a)' ) 'RATMAT_DET - Fatal error!'
            write ( *, '(a)' ) '  An overflow occurred.'
            write ( *, '(a)' ) 
     &        '  The determinant calculation cannot be done'
            write ( *, '(a)' ) '  for this matrix.'
            idtop = 0
            idbot = 1
            stop
          end if

        end do

        itop1 = itop
        ibot1 = ibot

        itop2 = idtop
        ibot2 = idbot

        call rat_add ( itop1, ibot1, itop2, ibot2, itop, ibot, ierror )

        if ( ierror .eq. 0 ) then
          idtop = itop
          idbot = ibot
        else
          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) 'RATMAT_DET - Fatal error!'
          write ( *, '(a)' ) '  An overflow occurred.'
          write ( *, '(a)' ) 
     &      '  The determinant calculation cannot be done'
          write ( *, '(a)' ) '  for this matrix.'
          idtop = 0
          idbot = 1
          stop
        end if

        if ( .not. more ) then
          go to 20
        end if

      go to 10

20    continue
c
c  The bottom should be positive.
c
      if ( idbot .lt. 0 ) then
        idbot = -idbot
        idtop = -idtop
      end if

      return
      end
      subroutine ratmat_print ( m, n, a, b, title )

c*********************************************************************72
c
cc RATMAT_PRINT prints out rational vectors or matrices.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    26 January 2007
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer M, N, the number of rows and columns in the matrix.
c
c    Input, integer A(M,N), B(M,N), the current rational or decimal
c    matrix.
c
c    Input, character * ( * ) TITLE, a label for the object being printed.
c
      implicit none

      integer m
      integer n

      integer a(m,n)
      integer b(m,n)
      character * ( 10 ) chrtmp2
      character * ( 10 ) chrtmp3
      character * ( 40 ) format1
      character * ( 40 ) format2
      integer i
      integer ione
      integer itemp
      integer j
      integer jmax
      integer jmin
      integer kmax
      integer ncolum
      parameter ( ncolum = 80 )
      integer none
      integer npline
      character * ( 100 ) output
      integer output_length
      character * ( * ) title
      integer title_length
c
c  Figure out how many rationals we can get in NCOLUM columns.
c
      kmax = 3

      do i = 1, m
        do j = 1, n

          itemp = abs ( a(i,j) )

10        continue

          if ( 10**(kmax-2) .le. itemp ) then
            kmax = kmax + 1
            go to 10
          end if

          itemp = abs ( b(i,j) )

20        continue

          if ( 10**(kmax-2) .lt. itemp ) then
            kmax = kmax + 1
            go to 20
          end if

        end do
      end do

      kmax = kmax + 1
      npline = ncolum / kmax
c
c  Create the formats.
c
      call i4_to_s_left ( npline, chrtmp2 )
      call i4_to_s_left ( kmax, chrtmp3 )

      format1(1:1) = '('
      format1(2:11) = chrtmp2
      format1(12:12) = 'i'
      format1(13:22) = chrtmp3
      format1(23:23) = ')'

      call s_blank_delete ( format1 )

      format2(1:1) = '('
      format2(2:11) = chrtmp2
      format2(12:12) = 'i'
      format2(13:22) = chrtmp3
      format2(23:23) = ')'

      call s_blank_delete ( format2 )

      do jmin = 1, n, npline

        jmax = min ( jmin + npline - 1, n )

        write ( *, '(a)' ) ' '

        if ( jmin .eq. 1 ) then
          title_length = len_trim ( title )
          write ( *, '(a)' ) title(1:title_length)
          write ( *, '(a)' ) ' '
        end if

        if ( 1 .lt. jmin .or. jmax .lt. n ) then
          write ( output, * ) 'Columns ', jmin, ' to ', jmax
          call s_blanks_delete ( output )
          output_length = len_trim ( output )
          write ( *, '(a)' ) output(1:output_length)
          write ( *, '(a)' ) ' '
        end if

        do i = 1, m

          write ( *, format1 ) ( a(i,j), j = jmin, jmax )
          write ( output, format1 ) ( b(i,j), j = jmin, jmax )
c
c  Delete each denominator that is 1.  If all are 1, don't
c  even print out the line.
c
          none = 0

          do j = jmin, jmax

            if ( b(i,j) .eq. 1 ) then
              ione = ( j - jmin + 1 ) * kmax
              output(ione:ione) = ' '
            else
              none = 1
            end if

          end do

          output_length = len_trim ( output )
          write ( *, '(a)' ) output(1:output_length)

          if ( jmax .eq. n .and. i .eq. m ) then
          else
            write ( *, '(a)' ) ' '
          end if

        end do

      end do

      return
      end
      subroutine regro_next ( n, v, vmax, done )

c*********************************************************************72
c
cc REGRO_NEXT computes restricted growth functions one at a time.
c
c  Discussion:
c
c    A restricted growth function on N is a vector (V(1), ..., V(N) )
c    of values V(I) between 1 and N, satisfying the requirements:
c      V(1) = 1;
c      V(I) .le. 1 + max ( V(1), V(2), ..., V(I-1) ).
c
c    The number of restricted growth functions on N is equal to
c    the Bell number B(N).
c
c    There is a bijection between restricted growth functions on N
c    and set partitions of N.
c
c  Example:
c
c    The 15 restricted growth functions for N = 4 are:
c
c    (1111), (1112), (1121), (1122), (1123),
c    (1211), (1212), (1213), (1221), (1222),
c    (1223), (1231), (1232), (1233), (1234).
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    26 January 2007
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Dennis Stanton, Dennis White,
c    Constructive Combinatorics,
c    Springer, 1986,
c    ISBN: 0387963472.
c
c  Parameters:
c
c    Input, integer N, the number of components in the restricted 
c    growth function.
c
c    Input/output, integer V(N).  The user need not set this quantity
c    before the initial call, and should not alter it between successive
c    calls.  On each return from the routine, with DONE = FALSE,
c    V will contain the componentwise values of the next restricted
c    growth function.
c
c    Input/output, integer VMAX(N).  The user need not set this quantity
c    before the initial call, and should not alter it between calls.
c    VMAX(I) records the largest value that component V(I) could take,
c    given the values of components 1 through I-1.
c
c    Input/output, logical DONE.
c    On first call, set DONE to TRUE, and then do not alter it.
c    On output, DONE will be FALSE if the routine has computed another
c    restricted growth function, or TRUE if all the restricted
c    growth functions have been returned.
c
      implicit none

      integer n

      logical done
      integer i
      integer j
      integer v(n)
      integer vmax(n)
c
c  First call:
c
      if ( done ) then

        do i = 1, n
          v(i) = 1
        end do

        vmax(1) = 1
        do i = 2, n
          vmax(i) = 2
        end do

        done = .false.
c
c  Later calls.
c
      else

        j = n

10      continue

          if ( j .eq. 1 ) then
            done = .true.
            return
          end if

          if ( v(j) .ne. vmax(j) ) then
            go to 20
          end if

          j = j - 1

        go to 10

20      continue

        v(j) = v(j) + 1

        do i = j+1, n

          v(i) = 1

          if ( v(j) .eq. vmax(j) ) then
            vmax(i) = vmax(j) + 1
          else
            vmax(i) = vmax(j)
          end if

        end do

      end if

      return
      end
      subroutine rfrac_to_cfrac ( m, p, q, t, ierror )

c*********************************************************************72
c
cc RFRAC_TO_CFRAC converts a rational polynomial fraction to a continued fraction.
c
c  Discussion:
c
c    That is, it accepts
c
c      P(1) + P(2) * X + ... + P(M) * X**(M-1)
c      -------------------------------------------------------
c      Q(1) + Q(2) * X + ... + Q(M) * X**(M-1) + Q(M+1) * X**M
c
c    and returns the equivalent continued fraction:
c
c      1 / ( T(1) + X / ( T(2) + X / (...T(2*M-1) + X / ( T(2*M) ... )))
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
c    Original FORTRAN77 version by John Hart, Ward Cheney, Charles Lawson, 
c    Hans Maehly, Charles Mesztenyi, John Rice, Henry Thatcher, 
c    Christoph Witzgall.
c    This FORTRAN77 version by John Burkardt.
c
c  Reference:
c
c    John Hart, Ward Cheney, Charles Lawson, Hans Maehly, Charles Mesztenyi, 
c    John Rice, Henry Thatcher, Christoph Witzgall,
c    Computer Approximations,
c    Wiley, 1968.
c
c  Parameters:
c
c    Input, integer M, defines the number of P coefficients,
c    and is one less than the number of Q coefficients, and one
c    half the number of T coefficients.
c
c    Input, double precision P(M), Q(M+1), the coefficients defining 
c    the rational polynomial fraction.
c
c    Output, double precision T(2*M), the coefficients defining the 
c    continued fraction.
c
c    Output, integer IERROR, error flag.
c    0, no error;
c    nonzero, the algorithm broke down at some point with a zero divisor.
c
      implicit none

      integer m

      double precision a(m+1,2*m+1)
      integer i
      integer ierror
      integer k
      double precision p(m)
      double precision q(m+1)
      double precision t(2*m)
      double precision ta

      ierror = 0

      do i = 1, m + 1
        a(i,1) = q(i)
      end do

      do i = 1, m
        a(i,2) = p(i)
      end do

      t(1) = a(1,1) / a(1,2)
      ta = a(m+1,1)

      do i = 1, m
        a(m-i+1,2*i+1) = ta
      end do

      do k = 1, 2*m-2

        do i = 1, (2*m-k)/2
          a(i,k+2) = a(i+1,k) - t(k) * a(i+1,k+1)
        end do

        if ( a(1,k+2) .eq. 0.0D+00 ) then
          ierror = 1
          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) 'RFRAC_TO_CFRAC - Fatal error!'
          write ( *, '(a,i8)' ) '  A(1,K+2) is zero for K = ', k
          stop
        end if

        t(k+1) = a(1,k+1) / a(1,k+2)

      end do

      t(2*m) = a(1,2*m) / a(1,2*m+1)

      return
      end
      subroutine rfrac_to_jfrac ( m, p, q, r, s )

c*********************************************************************72
c
cc RFRAC_TO_JFRAC converts a rational polynomial fraction to a J fraction.
c
c  Discussion:
c
c    The routine accepts
c
c    P(1) + P(2) * X + ... + P(M) * X**(M-1)
c    -------------------------------------------------------
c    Q(1) + Q(2) * X + ... + Q(M) * X**(M-1) + Q(M+1) * X**M
c
c    and returns the equivalent J-fraction:
c
c    R(1)/ ( X + S(1) + R(2) / ( X + S(2) + R(3) / ... + R(M) / ( X + S(M) )... ))
c
c    Thanks to Henry Amuasi for noticing and correcting an error in a
c    previous formulation of this routine, 02 October 2010.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    02 October 2010
c
c  Author:
c
c    Original FORTRAN77 version by John Hart, Ward Cheney, Charles Lawson, 
c    Hans Maehly, Charles Mesztenyi, John Rice, Henry Thatcher, 
c    Christoph Witzgall.
c    This FORTRAN77 version by John Burkardt.
c
c  Reference:
c
c    John Hart, Ward Cheney, Charles Lawson, Hans Maehly, Charles Mesztenyi, 
c    John Rice, Henry Thatcher, Christoph Witzgall,
c    Computer Approximations,
c    Wiley, 1968.
c
c  Parameters:
c
c    Input, integer M, defines the number of P, R, and S coefficients,
c    and is one less than the number of Q coefficients.
c    1 <= M.
c
c    Input, double precision P(M), Q(M+1), the coefficients defining 
c    the rational polynomial fraction.
c
c    Output, double precision R(M), S(M), the coefficients defining the
c    J-fraction.
c
      implicit none

      integer m

      double precision a(m+1,m+1)
      integer i
      integer k
      double precision p(m)
      double precision q(m+1)
      double precision r(m)
      double precision s(m)

      if ( m .lt. 1 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'RFRAC_TO_JFRAC - Fatal error!'
        write ( *, '(a)' ) '  Input M < 1.'
        stop
      end if

      do i = 1, m + 1
        a(i,1) = q(i)
      end do

      do i = 1, m
        a(i,2) = p(i)
      end do

      if ( 1 .lt. m ) then

        r(1) = a(m,2) / a(m+1,1)
        s(1) = ( r(1) * a(m,1) - a(m-1,2) ) / a(m,2)

        do k = 1, m - 2

          a(1,k+2) = r(k) * a(1,k) - s(k) * a(1,k+1)

          do i = 2, m - k
            a(i,k+2) = r(k) * a(i,k) - a(i-1,k+1) - s(k) * a(i,k+1)
          end do

          if ( a(m-k,k+2) .eq. 0.0D+00 ) then
            write ( *, '(a)' ) ' '
            write ( *, '(a)' ) 'RFRAC_TO_JFRAC - Fatal error!'
            write ( *, '(a,i8)' ) '  A(M-K,K+2) = 0 for K=', k
            stop
          end if

          r(k+1) = a(m-k,k+2) / a(m-k+1,k+1)
          s(k+1) = ( r(k+1) * a(m-k,k+1) - a(m-k-1,k+2) ) / a(m-k,k+2)

        end do

        a(1,m+1) = r(m-1) * a(1,m-1) - s(m-1) * a(1,m)

      end if

      r(m) = a(1,m+1) / a(2,m)
      s(m) = a(1,m) / a(2,m)

      return
      end
      subroutine s_adjustr ( s )

c*********************************************************************72
c
cc S_ADJUSTR flushes a string right.
c
c  Example:
c
c    Input             Output
c    'Hello     '      '     Hello'
c    ' Hi there!  '    '   Hi there!'
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    14 April 1999
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input/output, character * ( * ) S, on output, trailing blank
c    characters have been cut, and pasted back onto the front.
c
      implicit none

      integer i
      integer nonb
      character * ( * ) s
      integer s_length
c
c  Check the full length of the string.
c
      s_length = len ( s )
c
c  Find the occurrence of the last nonblank.
c
      nonb = len_trim ( s )
c
c  Shift the string right.
c
      do i = s_length, s_length + 1 - nonb, -1
        s(i:i) = s(i-s_length+nonb:i-s_length+nonb)
      end do
c
c  Blank out the beginning of the string.
c
      s(1:s_length-nonb) = ' '

      return
      end
      subroutine s_blank_delete ( s )

c*********************************************************************72
c
cc S_BLANK_DELETE removes blanks from a string, left justifying the remainder.
c
c  Discussion:
c
c    All TAB characters are also removed.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    26 July 1998
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input/output, character*(*) S, the string to be transformed.
c
      implicit none

      character ch
      integer get
      integer put
      character*(*) s
      integer s_len_trim
      integer s_length
      character tab

      tab = char ( 9 )

      put = 0
      s_length = s_len_trim ( s )

      do get = 1, s_length

        ch = s(get:get)

        if ( ch .ne. ' ' .and. ch .ne. tab ) then
          put = put + 1
          s(put:put) = ch
        end if

      end do

      s(put+1:s_length) = ' '

      return
      end
      subroutine s_blanks_delete ( s )

c*********************************************************************72
c
cc S_BLANKS_DELETE replaces consecutive blanks by one blank.
c
c  Discussion:
c
c    Thanks to Bill Richmond for pointing out a programming flaw which
c    meant that, as characters were slid to the left through multiple
c    blanks, their original images were not blanked out.  This problem
c    is easiest resolved by using a copy of the string.
c
c    The remaining characters are left justified and right padded with blanks.
c    TAB characters are converted to spaces.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    03 January 2007
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input/output, character*(*) S, the string to be transformed.
c
      implicit none

      integer i
      integer j
      character newchr
      character oldchr
      character*(*) s
      character*255 s_copy
      integer s_length
      character tab

      tab = char ( 9 )

      s_length = len ( s )

      j = 0
      s_copy(1:s_length) = s(1:s_length)
      s(1:s_length) = ' '

      newchr = ' '

      do i = 1, s_length

        oldchr = newchr
        newchr = s_copy(i:i)

        if ( newchr .eq. tab ) then
          newchr = ' '
        end if
  
        if ( oldchr .ne. ' ' .or. newchr .ne. ' ' ) then
          j = j + 1
          s(j:j) = newchr
        end if

      end do

      return
      end
      function s_eqi ( s1, s2 )

c*********************************************************************72
c
cc S_EQI is a case insensitive comparison of two strings for equality.
c
c  Example:
c
c    S_EQI ( 'Anjana', 'ANJANA' ) is TRUE.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    03 January 2007
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, character*(*) S1, S2, the strings to compare.
c
c    Output, logical S_EQI, the result of the comparison.
c
      implicit none

      character c1
      character c2
      integer i
      integer lenc
      logical s_eqi
      character*(*) s1
      integer s1_length
      character*(*) s2
      integer s2_length

      s1_length = len ( s1 )
      s2_length = len ( s2 )
      lenc = min ( s1_length, s2_length )

      s_eqi = .false.

      do i = 1, lenc

        c1 = s1(i:i)
        c2 = s2(i:i)
        call ch_cap ( c1 )
        call ch_cap ( c2 )

        if ( c1 .ne. c2 ) then
          return
        end if

      end do

      do i = lenc + 1, s1_length
        if ( s1(i:i) .ne. ' ' ) then
          return
        end if
      end do

      do i = lenc + 1, s2_length
        if ( s2(i:i) .ne. ' ' ) then
          return
        end if
      end do

      s_eqi = .true.

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
      subroutine schroeder ( n, s )

c*********************************************************************72
c
cc SCHROEDER generates the Schroeder numbers.
c
c  Discussion:
c
c    The Schroeder number S(N) counts the number of ways to insert
c    parentheses into an expression of N items, with two or more items within
c    a parenthesis.
c
c    Note that the Catalan number C(N) counts the number of ways
c    to legally arrange a set of N left and N right parentheses.
c
c  Example:
c
c    N = 4
c
c    1234
c    12(34)
c    1(234)
c    1(2(34))
c    1(23)4
c    1((23)4)
c    (123)4
c    (12)34
c    (12)(34)
c    (1(23))4
c    ((12)3)4
c
c  First Values:
c
c           1
c           1
c           3
c          11
c          45
c         197
c         903
c        4279
c       20793
c      103049
c      518859
c     2646723
c    13648869
c    71039373
c
c  Formula:
c
c    S(N) = ( P(N)(3.0) - 3 P(N-1)(3.0) ) / ( 4 * ( N - 1 ) )
c    where P(N)(X) is the N-th Legendre polynomial.
c
c  Recursion:
c
c    S(1) = 1
c    S(2) = 1
c    S(N) = ( ( 6 * N - 9 ) * S(N-1) - ( N - 3 ) * S(N-2) ) / N
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
c  Reference:
c
c    RP Stanley,
c    Hipparchus, Plutarch, Schroeder, and Hough,
c    American Mathematical Monthly,
c    Volume 104, Number 4, 1997, pages 344-350.
c
c    Laurent Habsieger, Maxim Kazarian, Sergei Lando,
c    On the Second Number of Plutarch,
c    American Mathematical Monthly,
c    May 1998, page 446.
c
c  Parameters:
c
c    Input, integer N, the number of Schroeder numbers desired.
c
c    Output, integer S(N), the Schroeder numbers.
c
      implicit none

      integer n

      integer i
      integer s(n)

      if ( n .le. 0 ) then
        return
      end if

      s(1) = 1

      if ( n .le. 1 ) then
        return
      end if

      s(2) = 1

      if ( n .le. 2 ) then
        return
      end if

      do i = 3, n
        s(i) = ( ( 6 * i - 9 ) * s(i-1) - ( i - 3 ) * s(i-2) ) / i
      end do

      return
      end
      subroutine sort_heap_external ( n, indx, i, j, isgn )

c*********************************************************************72
c
cc SORT_HEAP_EXTERNAL externally sorts a list of items into ascending order.
c
c  Discussion:
c
c    The actual list of data is not passed to the routine.  Hence this
c    routine may be used to sort integers, reals, numbers, names,
c    dates, shoe sizes, and so on.  After each call, the routine asks
c    the user to compare or interchange two items, until a special
c    return value signals that the sorting is completed.
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
c    Original FORTRAN77 version by Albert Nijenhuis, Herbert Wilf.
c    This FORTRAN77 version by John Burkardt.
c
c  Reference:
c
c    Albert Nijenhuis, Herbert Wilf,
c    Combinatorial Algorithms for Computers and Calculators,
c    Second Edition,
c    Academic Press, 1978,
c    ISBN: 0-12-519260-6,
c    LC: QA164.N54.
c
c  Parameters:
c
c    Input, integer N, the number of items to be sorted.
c
c    Input/output, integer INDX, the main communication signal.
c
c    The user must set INDX to 0 before the first call.
c    Thereafter, the user should not change the value of INDX until
c    the sorting is done.
c
c    On return, if INDX is
c
c      greater than 0,
c      * interchange items I and J;
c      * call again.
c
c      less than 0,
c      * compare items I and J;
c      * set ISGN = -1 if I .lt. J, ISGN = +1 if J .lt. I;
c      * call again.
c
c      equal to 0, the sorting is done.
c
c    Output, integer I, J, the indices of two items.
c    On return with INDX positive, elements I and J should be interchanged.
c    On return with INDX negative, elements I and J should be compared, and
c    the result reported in ISGN on the next call.
c
c    Input, integer ISGN, results of comparison of elements I and J.
c    (Used only when the previous call returned INDX less than 0).
c    ISGN .le. 0 means I is less than or equal to J;
c    0 .le. ISGN means I is greater than or equal to J.
c
      implicit none

      integer i
      integer i_save
      integer indx
      integer isgn
      integer j
      integer j_save
      integer k
      integer k1
      integer n
      integer n1

      save i_save
      save j_save
      save k
      save k1
      save n1

      data i_save / 0 /
      data j_save / 0 /
      data k / 0 /
      data k1 / 0 /
      data n1 / 0 /
c
c  INDX = 0: This is the first call.
c
      if ( indx .eq. 0 ) then

        i_save = 0
        j_save = 0
        k = n / 2
        k1 = k
        n1 = n
c
c  INDX .lt. 0: The user is returning the results of a comparison.
c
      else if ( indx .lt. 0 ) then

        if ( indx .eq. -2 ) then

          if ( isgn .lt. 0 ) then
            i_save = i_save + 1
          end if

          j_save = k1
          k1 = i_save
          indx = -1
          i = i_save
          j = j_save
          return

        end if

        if ( 0 .lt. isgn ) then
          indx = 2
          i = i_save
          j = j_save
          return
        end if

        if ( k .le. 1 ) then

          if ( n1 .eq. 1 ) then
            i_save = 0
            j_save = 0
            indx = 0
          else
            i_save = n1
            n1 = n1 - 1
            j_save = 1
            indx = 1
          end if

          i = i_save
          j = j_save
          return

        end if

        k = k - 1
        k1 = k
c
c  0 .lt. INDX, the user was asked to make an interchange.
c
      else if ( indx .eq. 1 ) then

        k1 = k

      end if

10    continue

        i_save = 2 * k1

        if ( i_save .eq. n1 ) then
          j_save = k1
          k1 = i_save
          indx = -1
          i = i_save
          j = j_save
          return
        else if ( i_save .le. n1 ) then
          j_save = i_save + 1
          indx = -2
          i = i_save
          j = j_save
          return
        end if

        if ( k .le. 1 ) then
          go to 20
        end if

        k = k - 1
        k1 = k

      go to 10

20    continue

      if ( n1 .eq. 1 ) then
        i_save = 0
        j_save = 0
        indx = 0
        i = i_save
        j = j_save
      else
        i_save = n1
        n1 = n1 - 1
        j_save = 1
        indx = 1
        i = i_save
        j = j_save
      end if

      return
      end
      subroutine sub_by_size_next ( n, a, size, more )

c*********************************************************************72
c
cc SUB_BY_SIZE_NEXT returns all subsets of an N set, in order of size.
c
c  Example:
c
c    N = 4:
c
c    1 2 3 4
c    1 2 3
c    1 2 4
c    1 3 4
c    1 3
c    1 4
c    2 3
c    1
c    2
c    3
c    (the empty set)
c
c  Discussion:
c
c    The subsets are returned in decreasing order of size, with the
c    empty set last.
c
c    For a given size K, the K subsets are returned in lexicographic order.
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
c    Input, integer N, the size of the set.
c
c    Input/output, integer A(N).  The entries A(1:SIZE) contain
c    the elements of the subset.  The elements are given in ascending
c    order.
c
c    Input/output, integer SIZE, the number of elements in the subset.
c
c    Input/output, logical MORE.  Set MORE = FALSE before first call
c    for a new sequence of subsets.  It then is set and remains
c    TRUE as long as the subset computed on this call is not the
c    final one.  When the final subset is computed, MORE is set to
c    FALSE as a signal that the computation is done.
c
      implicit none

      integer n

      integer a(n)
      logical more
      logical more2
      integer size

      save more2

      data more2 / .false. /

      if ( .not. more ) then
        more = .true.
        more2 = .false.
        size = n
      else if ( .not. more2 ) then
        size = size - 1
      end if
c
c  Compute the next subset of size SIZE.
c
      if ( 0 .lt. size ) then
        call ksub_next ( n, size, a, more2 )
      else if ( size .eq. 0 ) then
        more = .false.
      end if

      return
      end
      subroutine sub_gray_next ( n, a, more, ncard, iadd )

c*********************************************************************72
c
cc SUB_GRAY_NEXT generates all subsets of a set of order N, one at a time.
c
c  Discussion:
c
c    It generates the subsets one at a time, by adding or subtracting
c    exactly one element on each step.
c
c    This uses a Gray code ordering of the subsets.
c
c    The user should set MORE = FALSE and the value of N before
c    the first call.  On return, the user may examine A which contains
c    the definition of the new subset, and must check MORE, because
c    as soon as it is FALSE on return, all the subsets have been
c    generated and the user probably should cease calling.
c
c    The first set returned is the empty set.
c
c  Example:
c
c    N = 4
c
c    0 0 0 0
c    1 0 0 0
c    1 1 0 0
c    0 1 0 0
c    0 1 1 0
c    1 1 1 0
c    1 0 1 0
c    0 0 1 0
c    0 0 1 1
c    1 0 1 1
c    1 1 1 1
c    0 1 1 1
c    0 1 0 1
c    1 1 0 1
c    1 0 0 1
c    0 0 0 1
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    26 January 2007
c
c  Author:
c
c    Original FORTRAN77 version by Albert Nijenhuis, Herbert Wilf.
c    This FORTRAN77 version by John Burkardt.
c
c  Reference:
c
c    Albert Nijenhuis, Herbert Wilf,
c    Combinatorial Algorithms for Computers and Calculators,
c    Second Edition,
c    Academic Press, 1978,
c    ISBN: 0-12-519260-6,
c    LC: QA164.N54.
c
c  Parameters:
c
c    Input, integer N, the order of the total set from which
c    subsets will be drawn.
c
c    Input/output, integer A(N).  On each return, the Gray code for the newly
c    generated subset.  A(I) = 0 if element I is in the subset, 1 otherwise.
c
c    Input/output, logical MORE.  Set this variable FALSE before
c    the first call.  Normally, MORE will be returned TRUE but once
c    all the subsets have been generated, MORE will be
c    reset FALSE on return and you should stop calling the program.
c
c    Input/output, integer NCARD, the cardinality of the set returned,
c    which may be any value between 0 (the empty set) and N (the
c    whole set).
c
c    Output, integer IADD, the element which was added or removed to the
c    previous subset to generate the current one.  Exception:
c    the empty set is returned on the first call, and IADD is set to 0.
c
      implicit none

      integer n

      integer a(n)
      integer i
      integer iadd
      logical more
      integer ncard
c
c  The first set returned is the empty set.
c
      if ( .not. more ) then

        do i = 1, n
          a(i) = 0
        end do

        iadd = 0
        ncard = 0
        more = .true.

      else

        iadd = 1

        if ( mod ( ncard, 2 ) .ne. 0 ) then

10        continue

            iadd = iadd + 1
            if ( a(iadd-1) .ne. 0 ) then
              go to 20
            end if

          go to 10

20        continue

        end if

        a(iadd) = 1 - a(iadd)
        ncard = ncard + 2 * a(iadd) - 1
c
c  The last set returned is the singleton A(N).
c
        if ( ncard .eq. a(n) ) then
          more = .false.
        end if

      end if

      return
      end
      subroutine sub_gray_rank ( n, a, rank )

c*********************************************************************72
c
cc SUB_GRAY_RANK ranks a subset of an N set, using the Gray code ordering.
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
c    Input, integer N, the order of the total set from which
c    subsets will be drawn.
c
c    Input, integer A(N); A(I) is 1 if element I is in the set,
c    and 0 otherwise.
c
c    Output, integer RANK, the rank of the subset in the Gray code ordering.
c
      implicit none

      integer n

      integer a(n)
      integer gray
      integer rank

      call ubvec_to_ui4 ( n, a, gray )

      call gray_rank ( gray, rank )

      rank = rank + 1

      return
      end
      subroutine sub_gray_unrank ( rank, n, a )

c*********************************************************************72
c
cc SUB_GRAY_UNRANK produces a subset of an N set of the given Gray code rank.
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
c    Input, integer RANK, the rank of the subset in the Gray code ordering.
c
c    Input, integer N, the order of the total set from which
c    subsets will be drawn.
c
c    Output, integer A(N); A(I) is 1 if element I is in the set,
c    and 0 otherwise.
c
      implicit none

      integer n

      integer a(n)
      integer gray
      integer rank

      call gray_unrank ( rank-1, gray )

      call ui4_to_ubvec ( gray, n, a )

      return
      end
      subroutine sub_lex_next ( n, jmp, ndim, k, a )

c*********************************************************************72
c
cc SUB_LEX_NEXT generates the subsets of a set of N elements, one at a time.
c
c  Discussion:
c
c    The subsets are generated in lexicographical order.  
c
c    The routine can also be forced to generate only those subsets whose 
c    size is no greater than some user-specified maximum.
c
c  Example:
c
c    N = 5, JMP = ( K .eq. 3 )
c
c    1
c    1 2
c    1 2 3
c    1 2 4
c    1 2 5
c    1 3
c    1 3 4
c    1 3 5
c    1 4
c    1 4 5
c    1 5
c    2
c    2 3
c    2 3 4
c    2 3 5
c    2 4
c    2 4 5
c    2 5
c    3
c    3 4
c    3 4 5
c    3 5
c    4
c    4 5
c    5
c    empty set.
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
c    Original FORTRAN77 version by Albert Nijenhuis, Herbert Wilf.
c    This FORTRAN77 version by John Burkardt.
c
c  Reference:
c
c    Albert Nijenhuis, Herbert Wilf,
c    Combinatorial Algorithms for Computers and Calculators,
c    Second Edition,
c    Academic Press, 1978,
c    ISBN: 0-12-519260-6,
c    LC: QA164.N54.
c
c  Parameters:
c
c    Input, integer N, the order of the main set from which subsets
c    are chosen.
c
c    Input, logical JMP.  In the simplest case, set JMP = FALSE for
c    a normal computation.  But to jump over supersets of the input set,
c    set JMP = TRUE.  Setting JMP = ( K .eq. 3 ) before every new call
c    will, for example, force all the subsets returned
c    to have cardinality 3 or less.
c
c    Input, integer NDIM, the allowed storage for A.  If NDIM .lt. N,
c    JMP must be used to avoid creation of a subset too large to store in A.
c
c    Input/output, integer K.  On first call, the user must set K = 0 as
c    a startup signal to the program.  Thereafter, the routine returns
c    the size of the computed subset in K.  On the last return,
c    the empty set is returned and K is 0, which is a signal to
c    the user that the computation is complete.
c
c    Input/output, integer A(NDIM).  A(I) is the I-th element of the
c    subset, listed in increasing order, with 0's in entries
c    beyond entry K.
c
      implicit none

      integer ndim

      integer a(ndim)
      integer is
      logical jmp
      integer k
      integer n

      if ( k .le. 0 ) then

        if ( jmp ) then
          return
        end if

        is = 0
        k = 1
        a(1) = 1

      else if ( a(k) .ne. n ) then

        is = a(k)

        if ( .not. jmp ) then
          k = k + 1
        end if

        a(k) = is + 1

      else

        k = k - 1

        if ( k .ne. 0 ) then
          a(k) = a(k) + 1
        end if

      end if

      return
      end
      subroutine sub_random ( n, seed, a )

c*********************************************************************72
c
cc SUB_RANDOM selects a random subset of an N-set.
c
c  Example:
c
c    N = 4
c
c    0 0 1 1
c    0 1 0 1
c    1 1 0 1
c    0 0 1 0
c    0 0 0 1
c    1 1 0 0
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
c    Original FORTRAN77 version by Albert Nijenhuis, Herbert Wilf.
c    This FORTRAN77 version by John Burkardt.
c
c  Reference:
c
c    Albert Nijenhuis, Herbert Wilf,
c    Combinatorial Algorithms for Computers and Calculators,
c    Second Edition,
c    Academic Press, 1978,
c    ISBN: 0-12-519260-6,
c    LC: QA164.N54.
c
c  Parameters:
c
c    Input, integer N, the size of the full set.
c
c    Input/output, integer SEED, a seed for the random number generator.
c
c    Output, integer A(N).  A vector to hold the information about
c    the set chosen.  On return, if A(I) = 1, then
c    I is in the random subset, otherwise, A(I) = 0
c    and I is not in the random subset.
c
      implicit none

      integer n

      integer a(n)
      integer i
      integer i4_uniform
      integer seed

      do i = 1, n
        a(i) = i4_uniform ( 0, 1, seed )
      end do

      return
      end
      subroutine subcomp_next ( n, k, a, more, h, t )

c*********************************************************************72
c
cc SUBCOMP_NEXT computes the next subcomposition of N into K parts.
c
c  Discussion:
c
c    A composition of the integer N into K parts is an ordered sequence
c    of K nonnegative integers which sum to a value of N.
c
c    A subcomposition of the integer N into K parts is a composition
c    of M into K parts, where 0 .le. M .le. N.
c
c    A subcomposition of the integer N into K parts is also a lattice
c    point in the simplex whose vertices are the origin, and the K direction
c    vectors N*E(I) for I = 1 to K.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    02 July 2008
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the integer whose subcompositions are desired.
c
c    Input, integer K, the number of parts in the subcomposition.
c
c    Input/output, integer A(K), the parts of the subcomposition.
c
c    Input/output, logical MORE, set by the user to start the computation,
c    and by the routine to terminate it.
c
c    Input/output, integer H, T, two internal parameters needed for the
c    computation.  The user should allocate space for these in the calling
c    program, include them in the calling sequence, but never alter them!
c
      implicit none

      integer k

      integer a(k)
      integer h
      integer i
      logical more
      logical more2
      integer n
      integer n2
      integer t

      save more2
      save n2

      data more2 / .false. /
      data n2 / 0 /
c
c  The first computation.
c
      if ( .not. more ) then

        more = .true.

        do i = 1, k
          a(i) = 0
        end do

        n2 = 0
        more2 = .false.
c
c  Do the next element at the current value of N.
c    
      else if ( more2 ) then

        call comp_next ( n2, k, a, more2, h, t )

      else

        more2 = .false.
        n2 = n2 + 1

        call comp_next ( n2, k, a, more2, h, t )
        
      end if
c
c  Termination occurs if MORE2 = FALSE and N2 = N.
c
      if ( .not. more2 .and. n2 .eq. n ) then
        more = .false.
      end if

      return
      end
      subroutine subcompnz_next ( n, k, a, more )

c*********************************************************************72
c
cc SUBCOMPNZ_NEXT computes the next subcomposition of N into K nonzero parts.
c
c  Discussion:
c
c    A composition of the integer N into K nonzero parts is an ordered sequence
c    of K positive integers which sum to a value of N.
c
c    A subcomposition of the integer N into K nonzero parts is a composition
c    of M into K nonzero parts, where 0 .lt. M .le. N.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    26 January 2007
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the integer whose subcompositions are desired.
c
c    Input, integer K, the number of parts in the subcomposition.
c    K must be no greater than N.
c
c    Input/output, integer A(K), the parts of the subcomposition.
c
c    Input/output, logical MORE, set by the user to start the computation,
c    and by the routine to terminate it.
c
      implicit none

      integer k

      integer a(k)
      integer i
      logical more
      logical more2
      integer n
      integer n2

      save more2
      save n2

      data more2 / .false. /
      data n2 / 0 /

      if ( n .lt. k ) then
        more = .false.
        do i = 1, k
          a(i) = -1
        end do
        return
      end if
c
c  The first computation.
c
      if ( .not. more ) then

        more = .true.

        do i = 1, k
          a(i) = 1
        end do
        n2 = k
        more2 = .false.
c
c  Do the next element at the current value of N.
c    
      else if ( more2 ) then

        call compnz_next ( n2, k, a, more2 )

      else

        more2 = .false.
        n2 = n2 + 1

        call compnz_next ( n2, k, a, more2 )
        
      end if
c
c  Termination occurs if MORE2 = FALSE and N2 = N.
c
      if ( .not. more2 .and. n2 .eq. n ) then
        more = .false.
      end if

      return
      end
      subroutine subcompnz2_next ( n_lo, n_hi, k, a, more )

c*********************************************************************72
c
cc SUBCOMPNZ2_NEXT computes the next subcomposition of N into K nonzero parts.
c
c  Discussion:
c
c    A composition of the integer N into K nonzero parts is an ordered sequence
c    of K positive integers which sum to a value of N.
c
c    A subcomposition of the integer N into K nonzero parts is a composition
c    of M into K nonzero parts, where 0 .lt. M .le. N.
c
c    This routine computes all compositions of K into nonzero parts which sum
c    to values between N_LO and N_HI.
c
c    The routine SUBCOMPNZ_NEXT can be regarded as a special case where 
c    N_LO = K.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    26 January 2007
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N_LO, N_HI, the range of values N for which
c    compositions are desired.
c
c    Input, integer K, the number of parts in the subcomposition.
c    K must be no greater than N_HI.
c
c    Input/output, integer A(K), the parts of the subcomposition.
c
c    Input/output, logical MORE, set by the user to start the computation,
c    and by the routine to terminate it.
c
      implicit none

      integer k

      integer a(k)
      integer i
      logical more
      logical more2
      integer n_hi
      integer n_lo
      integer n2

      save more2
      save n2

      data more2 / .false. /
      data n2 / 0 /

      if ( n_hi .lt. k ) then
        more = .false.
        do i = 1, k
          a(i) = -1
        end do
        return
      end if

      if ( n_hi .lt. n_lo ) then
        more = .false.
        do i = 1, k
          a(i) = -1
        end do
        return
      end if
c
c  The first computation.
c
      if ( .not. more ) then

        more = .true.

        n2 = max ( k, n_lo )
        more2 = .false.

        call compnz_next ( n2, k, a, more2 )
c
c  Do the next element at the current value of N.
c    
      else if ( more2 ) then

        call compnz_next ( n2, k, a, more2 )

      else

        n2 = n2 + 1

        call compnz_next ( n2, k, a, more2 )
        
      end if
c
c  Termination occurs if MORE2 = FALSE and N2 = N_HI.
c
      if ( .not. more2 .and. n2 .eq. n_hi ) then
        more = .false.
      end if

      return
      end
      subroutine thue_binary_next ( n, thue )

c*********************************************************************72
c
cc THUE_BINARY_NEXT returns the next element in a binary Thue sequence.
c
c  Discussion:
c
c    Thue demonstrated that arbitrarily long sequences of 0's and
c    1's could be generated which had the "cubefree" property.  In
c    other words, for a given string S, there was no substring W
c    such that S contained "WWW".  In fact, a stronger result holds:
c    if "a" is the first letter of W, it is never the case that S
c    contains the substring "WWa". 
c
c    In this example, the digits allowed are binary, that is, just
c    "0" and "1".  The replacement rules are:
c
c    "0" -> "01"
c    "1" -> "10"
c
c    This routine produces the next binary Thue sequence in a given series.
c    However, the input sequence must be a Thue sequence in order for
c    us to guarantee that the output sequence will also have the
c    cubic nonrepetition property.
c
c    Also, enough space must be set aside in THUE to hold the
c    output sequence.  This will always be twice the input
c    value of N.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    26 January 2007
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input/output, integer N.  On input, the length of the input sequence.
c    On output, the length of the output sequence.
c
c    Input/output, integer THUE(N).  On input, the initial Thue sequence, and on
c    output, the result of applying the substitution rules once.
c
      implicit none

      integer n
      integer n_max
      parameter ( n_max = 100 )

      integer i
      integer n_out
      integer thue(n)
      integer thue_out(2*n_max)

      if ( n_max < n ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'THUE_TERNARY_NEXT - Fatal error!'
        write ( *, '(a)' ) '  Input N exceeded internal limit N_MAX.'
        stop
      end if

      n_out = 0

      do i = 1, n

        if ( thue(i) .eq. 0 ) then
          n_out = n_out + 1
          thue_out(n_out) = 0
          n_out = n_out + 1
          thue_out(n_out) = 1
        else if ( thue(i) .eq. 1 ) then
          n_out = n_out + 1
          thue_out(n_out) = 1
          n_out = n_out + 1
          thue_out(n_out) = 0
        else
          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) 'THUE_BINARY_NEXT - Fatal error!'
          write ( *, '(a)' ) 
     &      '  The input sequence contains a non-binary digit'
          write ( *, '(a,i8,a,i8)' ) '  THUE(', i, ') = ', thue(i)
          stop
        end if

      end do

      n = n_out

      do i = 1, n
        thue(i) = thue_out(i)
      end do

      return
      end
      subroutine thue_ternary_next ( n, thue )

c*********************************************************************72
c
cc THUE_TERNARY_NEXT returns the next element in a ternary Thue sequence.
c
c  Discussion:
c
c    Thue was interested in showing that there were arbitrarily long
c    sequences of digits which never displayed a pair of contiguous
c    repetitions of any length.  That is, there was no occurrence of 
c    "00" or "1010" or "121121", anywhere in the string.  This makes
c    the string "squarefree".
c
c    To do this, he demonstrated a way to start with a single digit,
c    and to repeatedly apply a series of transformation rules to each 
c    digit of the sequence, deriving nonrepeating sequences of ever 
c    greater length.
c
c    In this example, the digits allowed are ternary, that is, just
c    "0", "1" and "2".  The replacement rules are:
c
c    "0" -> "12"
c    "1" -> "102"
c    "2" -> "0"
c
c    This routine produces the next Thue sequence in a given series.
c    However, the input sequence must be a Thue sequence in order for
c    us to guarantee that the output sequence will also have the 
c    nonrepetition property.
c
c    Also, enough space must be set aside in THUE to hold the
c    output sequence.  This will never be more than 3 times the input
c    value of N.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    26 January 2007
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Brian Hayes,
c    Third Base,
c    American Scientist, 
c    Volume 89, Number 6, pages 490-494, November-December 2001.
c
c  Parameters:
c
c    Input/output, integer N.  On input, the length of the input sequence.
c    On output, the length of the output sequence.
c
c    Input/output, integer THUE(N).  On input, the initial Thue sequence, and on
c    output, the result of applying the substitution rules once.
c
      implicit none

      integer n
      integer n_max
      parameter ( n_max = 100 )

      integer i
      integer n_out
      integer thue(n)
      integer thue_out(3*n_max)

      if ( n_max < n ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'THUE_TERNARY_NEXT - Fatal error!'
        write ( *, '(a)' ) '  Input N exceeded internal limit N_MAX.'
        stop
      end if

      n_out = 0

      do i = 1, n

        if ( thue(i) .eq. 0 ) then
          n_out = n_out + 1
          thue_out(n_out) = 1
          n_out = n_out + 1
          thue_out(n_out) = 2
        else if ( thue(i) .eq. 1 ) then
          n_out = n_out + 1
          thue_out(n_out) = 1
          n_out = n_out + 1
          thue_out(n_out) = 0
          n_out = n_out + 1
          thue_out(n_out) = 2
        else if ( thue(i) .eq. 2 ) then
          n_out = n_out + 1
          thue_out(n_out) = 0
        else
          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) 'THUE_TERNARY_NEXT - Fatal error!'
          write ( *, '(a)' ) 
     &      '  The input sequence contains a non-ternary digit'
          write ( *, '(a,i8,a,i8)' ) '  THUE(', i, ') = ', thue(i)
          stop
        end if

      end do

      n = n_out
      do i = 1, n
        thue(i) = thue_out(i)
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
      subroutine timestring ( string )

c*********************************************************************72
c
cc TIMESTRING writes the current YMDHMS date into a string.
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
c    Output, character * ( * ) STRING, contains the date information.
c    A character length of 40 should always be sufficient.
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
      character * ( * ) string
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

      write ( string, 
     &  '(i2,1x,a,1x,i4,2x,i2,a1,i2.2,a1,i2.2,a1,i3.3,1x,a)' ) 
     &  d, month(m), y, h, ':', n, ':', s, '.', mm, ampm


      return
      end
      subroutine triang ( n, zeta, p )

c*********************************************************************72
c
cc TRIANG renumbers elements in accordance with a partial ordering.
c
c  Discussion:
c
c    TRIANG is given a partially ordered set.  The partial ordering
c    is defined by a matrix ZETA, where element I is partially less than
c    or equal to element J if and only if ZETA(I,J) = 1.
c
c    TRIANG renumbers the elements with a permutation P so that if
c    element I is partially less than element J in the partial ordering,
c    then P(I) .lt. P(J) in the usual, numerical ordering.
c
c    In other words, the elements are relabeled so that their labels
c    reflect their ordering.  This is equivalent to relabeling the
c    matrix so that, on unscrambling it, the matrix would be upper
c    triangular.
c
c    Calling I4MAT_PERM or R8MAT_PERM with P used for both the row
c    and column permutations applied to matrix ZETA will result in
c    an upper triangular matrix.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    26 January 2007
c
c  Author:
c
c    Original FORTRAN77 version by Albert Nijenhuis, Herbert Wilf.
c    This FORTRAN77 version by John Burkardt.
c
c  Reference:
c
c    Albert Nijenhuis, Herbert Wilf,
c    Combinatorial Algorithms for Computers and Calculators,
c    Second Edition,
c    Academic Press, 1978,
c    ISBN: 0-12-519260-6,
c    LC: QA164.N54.
c
c  Parameters:
c
c    Input, integer N, the number of elements in the set.  
c
c    Input, integer ZETA(N,N), describes the partial ordering.  
c    ZETA(I,J) =:
c      0, for diagonal elements (I = J), or 
c         for unrelated elements, or
c         if J << I.
c      1, if I << J.
c
c    Output, integer P(N), a permutation of the elements that reflects
c    their partial ordering.  P(I) is the new label of element I, with
c    the property that if ZETA(I,J) = 1, that is, I << J,
c    then P(I) .lt. P(J) (in the usual ordering).
c
      implicit none

      integer n

      integer i
      integer ierror
      integer iq
      integer ir
      integer it
      integer l
      integer m
      integer p(n)
      integer zeta(n,n)
c
c  Make sure ZETA represents a partially ordered set.  In other words,
c  if ZETA(I,J) = 1, then ZETA(J,I) must NOT be 1.
c
      call pord_check ( n, zeta, ierror )

      if ( ierror .ne. 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'TRIANG - Fatal error!'
        write ( *, '(a)' ) '  The matrix ZETA does not represent a'
        write ( *, '(a)' ) '  partial ordering.'
        stop
      end if

      m = 0
      l = 0
      do i = 1, n
        p(i) = 0
      end do
c
c  Find the next value of M for which P(M) is 0.
c
10    continue

        m = m + 1

        if ( p(m) .eq. 0 ) then
          go to 20
        end if

        if ( m .eq. n ) then
          return
        end if

      go to 10

20    continue

      it = m + 1
      ir = m + 1

30    continue

        if ( ir .le. n ) then

          if ( p(ir) .eq. 0 .and. zeta(ir,m) .ne. 0 ) then
            p(ir) = m
            m = ir
            ir = it
          else
            ir = ir + 1
          end if

        else

          l = l + 1
          iq = p(m)
          p(m) = l

          if ( iq .ne. 0 ) then

            ir = m + 1
            m = iq

          else if ( m .eq. n ) then

            go to 60

          else

40          continue

              m = m + 1

              if ( p(m) .eq. 0 ) then
                go to 50
              end if

              if ( m .eq. n ) then
                return
              end if

            go to 40

50          continue

            it = m + 1
            ir = m + 1

          end if

        end if

      go to 30

60    continue

      return
      end
      subroutine tuple_next ( m1, m2, n, rank, x )

c*********************************************************************72
c
cc TUPLE_NEXT computes the next element of a tuple space.
c
c  Discussion:
c
c    The elements are N vectors.  Each entry is constrained to lie
c    between M1 and M2.  The elements are produced one at a time.
c    The first element is
c      (M1,M1,...,M1),
c    the second element is
c      (M1,M1,...,M1+1),
c    and the last element is
c      (M2,M2,...,M2)
c    Intermediate elements are produced in lexicographic order.
c
c  Example:
c
c    N = 2, M1 = 1, M2 = 3
c
c    INPUT        OUTPUT
c    -------      -------
c    Rank  X      Rank   X
c    ----  ---    -----  --- 
c    0     * *    1      1 1
c    1     1 1    2      1 2
c    2     1 2    3      1 3
c    3     1 3    4      2 1
c    4     2 1    5      2 2
c    5     2 2    6      2 3
c    6     2 3    7      3 1
c    7     3 1    8      3 2
c    8     3 2    9      3 3
c    9     3 3    0      0 0
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
c    Input, integer M1, M2, the minimum and maximum entries.
c
c    Input, integer N, the number of components.
c
c    Input/output, integer RANK, counts the elements.
c    On first call, set RANK to 0.  Thereafter, the output value of RANK
c    will indicate the order of the element returned.  When there are no 
c    more elements, RANK will be returned as 0.
c
c    Input/output, integer X(N), on input the previous tuple.
c    On output, the next tuple.
c
      implicit none

      integer n

      integer i
      integer m1
      integer m2
      integer rank
      integer x(n)

      if ( m2 .lt. m1 ) then
        rank = 0
        return
      end if

      if ( rank .le. 0 ) then

        do i = 1, n
          x(i) = m1
        end do
        rank = 1

      else

        rank = rank + 1
        i = n

10      continue

          if ( x(i) .lt. m2 ) then
            x(i) = x(i) + 1
            go to 20
          end if

          x(i) = m1

          if ( i .eq. 1 ) then
            rank = 0
            do i = 1, n
              x(i) = m1
            end do
            go to 20
          end if

          i = i - 1

        go to 10

20      continue

      end if

      return
      end
      subroutine tuple_next_fast ( m, n, rank, x )

c*********************************************************************72
c
cc TUPLE_NEXT_FAST computes the next element of a tuple space, "fast".
c
c  Discussion:
c
c    The elements are N vectors.  Each entry is constrained to lie
c    between 1 and M.  The elements are produced one at a time.
c    The first element is
c      (1,1,...,1)
c    and the last element is
c      (M,M,...,M)
c    Intermediate elements are produced in lexicographic order.
c
c    This code was written as a possibly faster version of TUPLE_NEXT.
c
c  Example:
c
c    N = 2,
c    M = 3
c
c    INPUT        OUTPUT
c    -------      -------
c    Rank          X
c    ----          ----
c   -1            -1 -1
c
c    0             1  1
c    1             1  2
c    2             1  3
c    3             2  1
c    4             2  2
c    5             2  3
c    6             3  1
c    7             3  2
c    8             3  3
c    9             1  1
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
c    Input, integer M, the maximum entry in any component.
c    M must be greater than 0.
c
c    Input, integer N, the number of components.
c    N must be greater than 0.
c
c    Input, integer RANK, indicates the rank of the tuple.
c    Typically, 0 <= RANK .lt. N**M.  Values of RANK greater than
c    N**M are legal and meaningful; they are equivalent to the
c    corresponding value mod (N**M).  If RANK .lt. 0, this indicates 
c    that this is the first call for the given values of (M,N).  
c    Initialization is done, and X is set to a dummy value.
c
c    Output, integer X(N), the next tuple, or a dummy value if
c    initialization has just been done.
c
      implicit none

      integer n
      integer n_max
      parameter ( n_max = 100 )

      integer base(n_max)
      integer i
      integer m
      integer rank
      integer x(n)

      save base

      if ( rank .lt. 0 ) then

        if ( m .le. 0 ) then
          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) 'TUPLE_NEXT_FAST - Fatal error!'
          write ( *, '(a)' ) '  The value M <= 0 is not allowed.'
          write ( *, '(a,i8)' ) '  M = ', m
          stop
        end if

        if ( n .le. 0 ) then
          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) 'TUPLE_NEXT_FAST - Fatal error!'
          write ( *, '(a)' ) '  The value N <= 0 is not allowed.'
          write ( *, '(a,i8)' ) '  N = ', n
          stop
        end if

        if ( n_max < n ) then
          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) 'TUPLE_NEXT_FAST - Fatal error!'
          write ( *, '(a)' ) '  N_MAX < N.'
          stop
        end if

        base(n) = 1
        do i = n-1, 1, -1
          base(i) = base(i+1) * m
        end do

        do i = 1, n
          x(i) = -1
        end do

      else

        do i = 1, n
          x(i) = mod ( rank / base(i), m ) + 1
        end do

      end if

      return
      end
      subroutine tuple_next_ge ( m, n, rank, x )

c*********************************************************************72
c
cc TUPLE_NEXT_GE computes the next "nondecreasing" element of a tuple space.
c
c  Discussion:
c
c    The elements are N vectors.  Each element is constrained to lie
c    between 1 and M, and to have components that are nondecreasing.
c    That is, for an element X, and any positive RANK,
c      X(I) <= X(I+RANK)
c
c    The elements are produced one at a time.
c    The first element is
c      (1,1,...,1)
c    and the last element is
c      (M,M,...,M)
c    Intermediate elements are produced in lexicographic order.
c
c  Example:
c
c    N = 3, M = 3
c
c    RANK  X
c    ----  -----
c       1  1 1 1
c       2  1 1 2
c       3  1 1 3
c       4  1 2 2
c       5  1 2 3
c       6  1 3 3
c       7  2 2 2
c       8  2 2 3
c       9  2 3 3
c      10  3 3 3
c       0  0 0 0
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
c    Input, integer M, the maximum entry.
c
c    Input, integer N, the number of components.
c
c    Input/output, integer RANK, counts the elements.
c    On first call, set RANK to 0.  Thereafter, RANK will indicate the
c    order of the element returned.  When there are no more elements,
c    RANK will be returned as 0.
c
c    Input/output, integer X(N), on input the previous tuple (except
c    on the first call, when the input value of X is not needed.)
c    On output, the next tuple.
c
      implicit none

      integer n

      integer i
      integer j
      integer m
      integer rank
      integer x(n)

      if ( m .lt. 1 ) then
        return
      end if

      if ( rank .le. 0 ) then
        do i = 1, n
          x(i) = 1
        end do
        rank = 1
        return
      end if

      do i = n, 1, -1

        if ( x(i) .lt. m ) then
          x(i) = x(i) + 1
          do j = i+1, n
            x(j) = x(i)
          end do
          rank = rank + 1
          return
        end if

      end do

      rank = 0
      do i = 1, n
        x(i) = 0
      end do

      return
      end
      subroutine tuple_next2 ( n, xmin, xmax, rank, x )

c*********************************************************************72
c
cc TUPLE_NEXT2 computes the next element of an integer tuple space.
c
c  Discussion:
c
c    The elements X are N vectors.
c
c    Each entry X(I) is constrained to lie between XMIN(I) and XMAX(I).
c
c    The elements are produced one at a time.
c
c    The first element is
c      (XMIN(1), XMIN(2), ..., XMIN(N)),
c    the second is (probably)
c      (XMIN(1), XMIN(2), ..., XMIN(N)+1),
c    and the last element is
c      (XMAX(1), XMAX(2), ..., XMAX(N))
c
c    Intermediate elements are produced in a lexicographic order, with
c    the first index more important than the last, and the ordering of
c    values at a fixed index implicitly defined by the sign of
c    XMAX(I) - XMIN(I).
c
c  Example:
c
c    N = 2,
c    XMIN = (/ 1, 10 /)
c    XMAX = (/ 3,  8 /)
c
c    RANK    X
c    ----  -----
c      1   1 10
c      2   1  9
c      3   1  8
c      4   2 10
c      5   2  9
c      6   2  8
c      7   3 10
c      8   3  9
c      9   3  8
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
c    Input, integer N, the number of components.
c
c    Input, integer XMIN(N), XMAX(N), the "minimum" and "maximum" entry values.
c    These values are minimum and maximum only in the sense of the lexicographic
c    ordering.  In fact, XMIN(I) may be less than, equal to, or greater
c    than XMAX(I).
c
c    Input/output, integer RANK, the rank of the item.  On first call,
c    set RANK to 0 to start up the sequence.  On return, if RANK is zero,
c    there are no more items in the sequence.
c
c    Input/output, integer X(N), on input the previous tuple.
c    On output, the next tuple.
c
      implicit none

      integer n

      integer i
      integer prod
      integer rank
      integer x(n)
      integer xmin(n)
      integer xmax(n)

      if ( rank .lt. 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'TUPLE_NEXT2 - Fatal error!'
        write ( *, '(a,i8)' ) '  Illegal value of RANK = ', rank
        stop
      end if

      prod = 1
      do i = 1, n
        prod = prod * ( 1 + abs ( xmax(i) - xmin(i) ) )
      end do

      if ( prod .lt. rank ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'TUPLE_NEXT2 - Fatal error!'
        write ( *, '(a,i8)' ) '  Illegal value of RANK = ', rank
        stop
      end if

      if ( rank .eq. 0 ) then
        do i = 1, n
          x(i) = xmin(i)
        end do
        rank = 1
        return
      end if

      rank = rank + 1
      i = n

10    continue

        if ( x(i) .ne. xmax(i) ) then
          x(i) = x(i) + sign ( 1, xmax(i) - xmin(i) )
          go to 20
        end if

        x(i) = xmin(i)

        if ( i .eq. 1 ) then
          rank = 0
          go to 20
        end if

        i = i - 1

      go to 10

20    continue

      return
      end
      subroutine ubvec_add ( n, bvec1, bvec2, bvec3 )

c*********************************************************************72
c
cc UBVEC_ADD adds two unsigned binary vectors.
c
c  Discussion:
c
c    A UBVEC is a vector of binary digits representing an unsigned integer.  
c
c    UBVEC(N) contains the units digit, UBVEC(N-1)
c    the coefficient of 2, UBVEC(N-2) the coefficient of 4 and so on,
c    so that printing the digits in order gives the binary form of the number.
c
c  Example:
c
c    N = 4
c
c     UBVEC1       +  UBVEC2       =  UBVEC3
c
c    ( 0 0 0 1 )   + ( 0 0 1 1 )   = ( 0 1 0 0 )
c
c      1           +   3           =   4
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    20 March 2009
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the length of the vectors.
c
c    Input, integer BVEC1(N), BVEC2(N), the vectors to be added.
c
c    Output, integer BVEC3(N), the sum of the two input vectors.
c
      implicit none

      integer n

      integer base
      parameter ( base = 2 )
      integer bvec1(n)
      integer bvec2(n)
      integer bvec3(n)
      integer i
      logical overflow

      overflow = .false.

      do i = 1, n
        bvec3(i) = bvec1(i) + bvec2(i)
      end do

      do i = n, 1, -1

10      continue

        if ( base .le. bvec3(i) ) then
          bvec3(i) = bvec3(i) - base
          if ( 1 .lt. i ) then
            bvec3(i-1) = bvec3(i-1) + 1
          else
            overflow = .true.
          end if
          go to 10
        end if

      end do

      return
      end
      subroutine ubvec_to_ui4 ( n, bvec, ui4 )

c*********************************************************************72
c
cc UBVEC_TO_UI4 makes an unsigned integer from an unsigned binary vector.
c
c  Discussion:
c
c    A UBVEC is a vector of binary digits representing an unsigned integer.  
c
c    UBVEC(N) contains the units digit, UBVEC(N-1)
c    the coefficient of 2, UBVEC(N-2) the coefficient of 4 and so on,
c    so that printing the digits in order gives the binary form of the number.
c
c  Example:
c
c    N = 4
c
c         BVEC   binary UI4
c    ----------  -----  --
c    1  2  3  4
c    ----------
c    0  0  0  1       1  1
c    0  0  1  0      10  2
c    0  0  1  1      11  3
c    0  1  0  0     100  4
c    1  0  0  1    1001  9
c    1  1  1  1    1111 15
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    19 January 2007
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the dimension of the vector.
c
c    Input, integer BVEC(N), the binary representation.
c
c    Output, integer UI4, the integer.
c
      implicit none

      integer n

      integer base
      parameter ( base = 2 )
      integer bvec(n)
      integer i
      integer ui4

      ui4 = 0
      do i = 1, n
        ui4 = base * ui4 + bvec(i)
      end do

      return
      end
      subroutine ui4_to_ubvec ( ui4, n, bvec )

c*********************************************************************72
c
cc UI4_TO_UBVEC makes an unsigned binary vector from an unsigned integer.
c
c  Discussion:
c
c    A UBVEC is a vector of binary digits representing an unsigned integer.  
c
c    UBVEC(N) contains the units digit, UBVEC(N-1)
c    the coefficient of 2, UBVEC(N-2) the coefficient of 4 and so on,
c    so that printing the digits in order gives the binary form of the number.
c
c    To guarantee that there will be enough space for any
c    value of I, it would be necessary to set N = 32.
c
c  Example:
c
c     I       BVEC         binary
c    --  ----------------  ------
c     1  1  0  0  0  0  1       1
c     2  0  1  0  0  1  0      10
c     3  1  1  0  0  1  1      11
c     4  0  0  0  1  0  0     100
c     9  0  0  1  0  0  1    1001
c    57  1  1  1  0  1  1  110111
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    20 March 2009
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer UI4, an integer to be represented.
c
c    Input, integer N, the dimension of the vector.
c
c    Output, integer BVEC(N), the binary representation.
c
      implicit none

      integer n

      integer base
      parameter ( base = 2 )
      integer bvec(n)
      integer i
      integer ui4
      integer ui4_copy

      ui4_copy = ui4

      do i = n, 1, -1

        bvec(i) = mod ( ui4_copy, base )

        ui4_copy = ui4_copy / base

      end do

      return
      end
      subroutine vec_gray_next ( n, base, a, done, change )

c*********************************************************************72
c
cc VEC_GRAY_NEXT computes the elements of a product space.
c
c  Discussion:
c
c    The elements are produced one at a time.
c
c    This routine handles the case where the number of degrees of freedom may
c    differ from one component to the next.
c
c    A method similar to the Gray code is used, so that successive
c    elements returned by this routine differ by only a single element.
c
c    The routine uses internal static memory.
c
c  Example:
c
c    N = 2, BASE = ( 2, 3 ), DONE = TRUE
c
c     A    DONE  CHANGE
c    ---  -----  ------
c    0 0  FALSE    1
c    0 1  FALSE    2
c    0 2  FALSE    2
c    1 2  FALSE    1
c    1 1  FALSE    2
c    1 0  FALSE    2
c    1 0   TRUE   -1  
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    26 January 2007
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Dennis Stanton, Dennis White,
c    Constructive Combinatorics,
c    Springer, 1986,
c    ISBN: 0387963472.
c
c  Parameters:
c
c    Input, integer N, the number of components.
c
c    Input, integer BASE(N), contains the number of degrees of
c    freedom of each component.  The output values of A will
c    satisfy 0 .le. A(I) .lt. BASE(I).
c
c    Input/output, integer A(N).  On the first call, the input value
c    of A doesn't matter.  Thereafter, it should be the same as
c    its output value from the previous call.  On output, if DONE
c    is FALSE, then A contains the next element of the space.
c
c    Input/output, logical DONE.  On the first call, the user must
c    set DONE to TRUE.  This signals the program to initialize data.
c    On every return, if DONE is FALSE, the program has computed
c    another entry, which is contained in A.  If DONE is TRUE,
c    then there are no more entries, and the program should not be
c    called for any more.
c
c    Output, integer CHANGE, is set to the index of the element whose
c    value was changed.  On return from the first call, CHANGE
c    is 1, even though all the elements have been "changed".  On
c    return with DONE equal to TRUE, CHANGE is -1.
c
      implicit none

      integer n
      integer n_max
      parameter ( n_max = 100 )

      integer a(n)
      integer active(n_max)
      integer base(n)
      integer change
      integer dir(n_max)
      logical done
      integer i

      save active
      save dir

      if ( n_max < n ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'VEC_GRAY_NEXT - Fatal error!'
        write ( *, '(a)' ) '  Input N exceeds internal limit.'
        stop
      end if
c
c  The user is calling for the first time.
c
      if ( done ) then

        done = .false.

        do i = 1, n
          a(i) = 0
        end do

        do i = 1, n
          dir(i) = 1
        end do

        do i = 1, n
          active(i) = 1
        end do

        do i = 1, n

          if ( base(i) .lt. 1 ) then
            write ( *, '(a)' ) ' '
            write ( *, '(a)' ) 'VEC_GRAY_NEXT - Warning!'
            write ( *, '(a,i8)' ) '  For index I = ',i
            write ( *, '(a,i8)' ) 
     &        '  the nonpositive value of BASE(I) = ', base(i)
            write ( *, '(a)' ) '  which was reset to 1c'
            base(i) = 1
            active(i) = 0
          else if ( base(i) .eq. 1 ) then
            active(i) = 0
          end if

        end do

        change = 1

        return

      end if
c
c  Seek the maximum active index.
c
      change = -1

      do i = n, 1, -1
        if ( active(i) .eq. 1 ) then
          change = i
          go to 10
        end if
      end do

10    continue
c
c  If there are NO active indices, we have generated all vectors.
c
      if ( change .eq. -1 ) then
        done = .true.
        return
      end if
c
c  Increment the element with maximum active index.
c
      a(change) = a(change) + dir(change)
c
c  If we attained a minimum or maximum value, reverse the direction
c  vector, and deactivate the index.
c
      if ( a(change) .eq. 0 .or. 
     &     a(change) .eq. base(change) - 1 ) then
        dir(change) = -dir(change)
        active(change) = 0
      end if
c
c  Activate all subsequent indices.
c
      do i = change + 1, n
        if ( 1 .lt. base(i) ) then
          active(i) = 1
        end if
      end do

      return
      end
      subroutine vec_lex_next ( dim_num, base, a, more )

c*********************************************************************72
c
cc VEC_LEX_NEXT generates vectors in lex order.
c
c  Discussion:
c
c    The vectors are produced in lexical order, starting with
c    (0,0,...,0),
c    (0,0,...,1),
c    ...
c    (BASE-1,BASE-1,...,BASE-1).
c
c  Example:
c
c    DIM_NUM = 2,
c    BASE = 3
c
c    0   0
c    0   1
c    0   2
c    1   0
c    1   1
c    1   2
c    2   0
c    2   1
c    2   2
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    01 June 2012
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer DIM_NUM, the size of the vectors to be used.
c
c    Input, integer BASE, the base to be used.  BASE = 2 will
c    give vectors of 0's and 1's, for instance.
c
c    Input/output, integer A(DIM_NUM), the next vector.
c
c    Input/output, logical MORE.  Set this variable FALSE before
c    the first call.  On return, MORE is TRUE if another vector has
c    been computed.  If MORE is returned FALSE, ignore the output
c    vector and stop calling the routine.
c
      implicit none

      integer dim_num

      integer a(dim_num)
      integer base
      integer i
      logical more

      if ( .not. more ) then

        do i = 1, dim_num
          a(i) = 0
        end do
        more = .true.

      else

        do i = dim_num, 1, -1

          a(i) = a(i) + 1

          if ( a(i) < base ) then
            return
          end if

          a(i) = 0

        end do

        more = .false.

      end if

      return
      end
      subroutine vec_next ( n, base, a, more )

c*********************************************************************72
c
cc VEC_NEXT generates all N-vectors of integers modulo a given base.
c
c  Discussion:
c
c    The vectors are produced in lexical order, starting with
c    (0,0,...,0), (0,0,...,1), ... through (BASE-1,BASE-1,...,BASE-1).
c
c  Example:
c
c    N = 2, 
c    BASE = 3
c
c    0   0
c    0   1
c    0   2
c    1   0
c    1   1
c    1   2
c    2   0
c    2   1
c    2   2
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    19 January 2007
c
c  Parameters:
c
c    Input, integer N, the size of the vectors to be used.
c
c    Input, integer BASE, the base to be used.  BASE = 2 will
c    give vectors of 0's and 1's, for instance.
c
c    Input/output, integer A(N).  On each return, A
c    will contain entries in the range 0 to BASE-1.
c
c    Input/output, logical MORE.  Set this variable FALSE before
c    the first call.  Normally, MORE will be returned TRUE but
c    once all the vectors have been generated, MORE will be
c    reset FALSE and you should stop calling the program.
c
      implicit none

      integer n

      integer a(n)
      integer base
      integer i
      integer kount
      integer last
      logical more
      integer nn

      save kount
      save last

      data kount / 0 /
      data last / 0 /

      if ( .not. more ) then

        kount = 1
        last = base**n
        more = .true.

        do i = 1, n
          a(i) = 0
        end do

      else

        kount = kount + 1

        if ( kount .eq. last ) then
          more = .false.
        end if

        a(n) = a(n) + 1

        do i = 1, n

          nn = n - i

          if ( a(nn+1) .lt. base ) then
            return
          end if

          a(nn+1) = 0

          if ( nn .ne. 0 ) then
            a(nn) = a(nn) + 1
          end if

        end do

      end if

      return
      end
      subroutine vec_random ( n, base, seed, a )

c*********************************************************************72
c
cc VEC_RANDOM selects a random N-vector of integers modulo a given base.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    19 January 2007
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the size of the vector to be generated.
c
c    Input, integer BASE, the base to be used.
c
c    Input/output, integer SEED, a random number seed.
c
c    Output, integer A(N), a list of N random values between
c    0 and BASE-1.
c
      implicit none

      integer n

      integer a(n)
      integer base
      integer i
      integer i4_uniform
      integer seed

      do i = 1, n
        a(i) = i4_uniform ( 0, base-1, seed )
      end do

      return
      end
      subroutine vec_rank ( n, base, a, rank )

c*********************************************************************72
c
cc VEC_RANK computes the rank of a product space element.
c
c  Discussion:
c
c    The rank applies only to the elements as produced by the routine
c    VEC_NEXT2.
c
c  Example:
c
c    N = 2, BASE = (/ 2, 3 /), A = ( 1, 2 ),
c
c    RANK = 4.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    18 January 2007
c
c  Reference:
c
c    Dennis Stanton, Dennis White,
c    Constructive Combinatorics,
c    Springer, 1986,
c    ISBN: 0387963472.
c
c  Parameters:
c
c    Input, integer N, the number of components.
c
c    Input, integer BASE(N), contains the number of degrees of
c    freedom of each component.  The output values of A will
c    satisfy 0 <= A(I) < BASE(I).
c
c    Input, integer A(N), the product space element, with the
c    property that 0 <= A(I) < BASE(I) for each entry I.
c
c    Output, integer RANK, the rank, or order, of the element in
c    the list of all elements.  The rank count begins at 1.
c
      implicit none

      integer n

      integer a(n)
      integer base(n)
      integer c
      integer i
      integer rank

      rank = 0

      do i = 1, n

        if ( mod ( rank, 2 ) .eq. 1 ) then
          c = base(i) - a(i) - 1
        else
          c = a(i)
        end if

        rank = base(i) * rank + c

      end do

      rank = rank + 1

      return
      end
      subroutine vec_unrank ( n, base, rank, a )

c*********************************************************************72
c
cc VEC_UNRANK computes the product space element of a given rank.
c
c  Discussion:
c
c    The rank applies only to the elements as produced by the routine
c    VEC_NEXT2.
c
c  Example:
c
c    N = 2, BASE = ( 2, 3 ), RANK = 4.
c
c    A = ( 1, 2 ).
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    19 January 2007
c
c  Reference:
c
c    Dennis Stanton, Dennis White,
c    Constructive Combinatorics,
c    Springer, 1986,
c    ISBN: 0387963472.
c
c  Parameters:
c
c    Input, integer N, the number of components.
c
c    Input, integer BASE(N), contains the number of degrees of
c    freedom of each component.  The output values of A will
c    satisfy 0 <= A(I) < BASE(I).
c
c    Input, integer RANK, the desired rank, or order, of the element in
c    the list of all elements.  The rank count begins at 1 and extends
c    to RANK_MAX = Product ( 1 <= I <= N ) BASE(I).
c
c    Output, integer A(N), the product space element of the given rank.
c
      implicit none

      integer n

      integer a(n)
      integer base(n)
      integer i
      integer rank
      integer s

      s = rank - 1

      do i = n, 1, -1

        a(i) = mod ( s, base(i) )
        s = s / base(i)

        if ( mod ( s, 2 ) .eq. 1 ) then
          a(i) = base(i) - a(i) - 1
        end if

      end do

      return
      end
      subroutine vector_constrained_next ( n, x_min, x_max, x, 
     &  constraint, more )

c*********************************************************************72
c
cc VECTOR_CONSTRAINED_NEXT returns the "next" constrained vector.
c
c  Discussion:
c
c    We consider all vectors of dimension N whose components
c    satisfy X_MIN(1:N) .le. X(1:N) .le. X_MAX(1:N).
c
c    We are only interested in the subset of these vectors which
c    satisfy the following constraint:
c
c      sum ( 1 .le. I .le. N ) ( ( X(I) - 1 ) / X_MAX(I) ) .le. 1
c
c    We can carry out this check using integer arithmetic if we
c    multiply through by P = product ( X_MAX(1:N) ):
c
c      sum ( 1 .le. I .le. N ) ( ( X(I) - 1 ) * ( P / X_MAX(I) ) ) .le. P.
c
c    This routine returns, one at a time, and in right-to-left
c    lexicographic order, exactly those vectors which satisfy
c    the constraint.
c
c  Example:
c
c    N = 3
c    X_MIN:   2   2   1
c    X_MAX:   4   5   3
c 
c    P = 60
c
c    #  X(1)  X(2)  X(3)  CONSTRAINT
c
c    1    2     2     1       27
c    2    3     2     1       42
c    3    4     2     1       57
c    4    2     3     1       39
c    5    3     3     1       54
c    6    2     4     1       51
c    7    2     2     2       47
c    8    2     3     2       59
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    26 January 2007
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the number of components in the vector.
c
c    Input, integer X_MIN(N), X_MAX(N), the minimum and maximum
c    values allowed in each component.
c
c    Input/output, integer X(N).  On first call (with MORE = FALSE), 
c    the input value of X is not important.  On subsequent calls, the
c    input value of X should be the output value from the previous call.
c    On output, (with MORE = TRUE), the value of X will be the "next"
c    vector in the reverse lexicographical list of vectors that satisfy
c    the condition.  However, on output with MORE = FALSE, the vector
c    X is meaningless, because there are no more vectors in the list.
c
c    Output, integer CONSTRAINT, the constraint value for X.  Valid vectors X
c    will have a CONSTRAINT value between product(X_MIN(1:N)) (automatically)
c    and product(X_MAX(1:N)) (because we skip over vectors with a
c    constraint larger than this value).
c
c    Input/output, logical MORE.  On input, if the user has set MORE
c    FALSE, the user is requesting the initiation of a new sequence
c    of values.  If MORE is TRUE, then the user is requesting "more"
c    values in the current sequence.  On output, if MORE is TRUE,
c    then another value was found and returned in X, but if MORE is
c    FALSE, then there are no more values in the sequence, and X is
c    NOT the next value.
c
      implicit none

      integer n

      integer constraint
      integer i
      integer j
      logical more
      integer x(n)
      integer x_max(n)
      integer x_min(n)
      integer x_prod

      save x_prod

      if ( .not. more ) then

        do j = 1, n
          x(j) = x_min(j)
        end do

        call i4vec_product ( n, x_max, x_prod )

        constraint = 0
        do j = 1, n
          constraint = constraint 
     &      + ( x(j) - 1 ) * ( x_prod / x_max(j) )
        end do

        if ( x_prod .lt. constraint ) then
          more = .false.
        else
          more = .true.
        end if

        return

      else

        i = 1

10      continue

          if ( x(i) .lt. x_max(i) ) then

            x(i) = x(i) + 1

            constraint = 0
            do j = 1, n
              constraint = constraint 
     &          + ( x(j) - 1 ) * ( x_prod / x_max(j) )
            end do

            if ( constraint .le. x_prod ) then
              go to 20
            end if

          end if

          x(i) = x_min(i)

          i = i + 1

          if ( n .lt. i ) then
            more = .false.
            go to 20
          end if

        go to 10

20      continue

      end if

      return
      end
      subroutine vector_constrained_next2 ( n, x_min, x_max, x, 
     &  constraint, more )

c*********************************************************************72
c
cc VECTOR_CONSTRAINED_NEXT2 returns the "next" constrained vector.
c
c  Discussion:
c
c    We consider all vectors of dimension N whose components
c    satisfy X_MIN(1:N) .le. X(1:N) .le. X_MAX(1:N).
c
c    We are only interested in the subset of these vectors which
c    satisfy the following constraint:
c
c      sum ( 1 .le. I .le. N ) ( X(I) / X_MAX(I) ) .le. 1
c
c    We can carry out this check using integer arithmetic if we
c    multiply through by P = product ( X_MAX(1:N) ):
c
c      sum ( 1 .le. I .le. N ) ( X(I) * ( P / X_MAX(I) ) ) .le. P.
c
c    This routine returns, one at a time, and in right-to-left
c    lexicographic order, exactly those vectors which satisfy
c    the constraint.
c
c  Example:
c
c    N = 3
c    X_MIN:   1   1   1
c    X_MAX:   5   6   4
c
c    P = 120
c
c    #  X(1)  X(2)  X(3)  CONSTRAINT
c
c    1    1     1     1       74
c    2    2     1     1       98
c    3    1     2     1       94
c    4    2     2     1      119
c    5    1     3     1      114
c    6    1     1     2      104
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    26 January 2007
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the number of components in the vector.
c
c    Input, integer X_MIN(N), X_MAX(N), the minimum and maximum
c    values allowed in each component.
c
c    Input/output, integer X(N).  On first call (with MORE = FALSE), 
c    the input value of X is not important.  On subsequent calls, the
c    input value of X should be the output value from the previous call.
c    On output, (with MORE = TRUE), the value of X will be the "next"
c    vector in the reverse lexicographical list of vectors that satisfy
c    the condition.  However, on output with MORE = FALSE, the vector
c    X is meaningless, because there are no more vectors in the list.
c
c    Output, integer CONSTRAINT, the constraint value for X.  Valid vectors X
c    will have a CONSTRAINT value between product(X_MIN(1:N)) (automatically)
c    and product(X_MAX(1:N)) (because we skip over vectors with a
c    constraint larger than this value).
c
c    Input/output, logical MORE.  On input, if the user has set MORE
c    FALSE, the user is requesting the initiation of a new sequence
c    of values.  If MORE is TRUE, then the user is requesting "more"
c    values in the current sequence.  On output, if MORE is TRUE,
c    then another value was found and returned in X, but if MORE is
c    FALSE, then there are no more values in the sequence, and X is
c    NOT the next value.
c
      implicit none

      integer n

      integer constraint
      integer i
      integer j
      logical more
      integer x(n)
      integer x_max(n)
      integer x_min(n)
      integer x_prod

      save x_prod

      if ( .not. more ) then

        do j = 1, n
          x(j) = x_min(j)
        end do

        call i4vec_product ( n, x_max, x_prod )

        constraint = 0.0D+00
        do j = 1, n
          constraint = constraint + x(j) * ( x_prod / x_max(j) )
        end do

        if ( x_prod .lt. constraint ) then
          more = .false.
        else
          more = .true.
        end if

        return

      else

        i = 1

10      continue

          if ( x(i) .lt. x_max(i) ) then

            x(i) = x(i) + 1

            constraint = 0.0D+00
            do j = 1, n
              constraint = constraint + x(j) * ( x_prod / x_max(j) )
            end do

            if ( constraint .le. x_prod ) then
              go to 20
            end if

          end if

          x(i) = x_min(i)

          i = i + 1

          if ( n .lt. i ) then
            more = .false.
            go to 20
          end if

        go to 10

20      continue

      end if

      return
      end
      subroutine vector_constrained_next3 ( n, x_min, x_max, x, 
     &  constraint, more )

c*********************************************************************72
c
cc VECTOR_CONSTRAINED_NEXT3 returns the "next" constrained vector.
c
c  Discussion:
c
c    This routine addresses the same problem as VECTOR_CONSTRAINED_NEXT2,
c    and differs only in that real arithmetic is used, rather than
c    integer arithmetic.  Integer arithmetic allows us to do an exact
c    calculation, but we run into overflow problems in simple cases
c    where N is 10 and the X_MAX entries are of order 10, for instance.
c
c    We consider all vectors of dimension N whose components
c    satisfy X_MIN(1:N) .le. X(1:N) .le. X_MAX(1:N).
c
c    We are only interested in the subset of these vectors which
c    satisfy the following constraint:
c
c      sum ( 1 .le. I .le. N ) ( X(I) / X_MAX(I) ) .le. 1
c
c    This routine returns, one at a time, and in right-to-left
c    lexicographic order, exactly those vectors which satisfy
c    the constraint.
c
c  Example:
c
c    N = 3
c    X_MIN:   1   1   1
c    X_MAX:   5   6   4
c
c    P = 120
c
c    #  X(1)  X(2)  X(3)  CONSTRAINT
c
c    1    1     1     1       0.62
c    2    2     1     1       0.82
c    3    1     2     1       0.78
c    4    2     2     1       0.98
c    5    1     3     1       0.95
c    6    1     1     2       0.87
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    26 January 2007
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the number of components in the vector.
c
c    Input, integer X_MIN(N), X_MAX(N), the minimum and maximum
c    values allowed in each component.
c
c    Input/output, integer X(N).  On first call (with MORE = FALSE), 
c    the input value of X is not important.  On subsequent calls, the
c    input value of X should be the output value from the previous call.
c    On output, (with MORE = TRUE), the value of X will be the "next"
c    vector in the reverse lexicographical list of vectors that satisfy
c    the condition.  However, on output with MORE = FALSE, the vector
c    X is meaningless, because there are no more vectors in the list.
c
c    Output, double precision CONSTRAINT, the constraint value for X.  
c    Valid vectors X will have a CONSTRAINT value between 
c      product(X_MIN(1:N)) / product(X_MAX(1:N))
c    and 1.0.
c
c    Input/output, logical MORE.  On input, if the user has set MORE
c    FALSE, the user is requesting the initiation of a new sequence
c    of values.  If MORE is TRUE, then the user is requesting "more"
c    values in the current sequence.  On output, if MORE is TRUE,
c    then another value was found and returned in X, but if MORE is
c    FALSE, then there are no more values in the sequence, and X is
c    NOT the next value.
c
      implicit none

      integer n

      double precision constraint
      integer i
      integer j
      logical more
      integer x(n)
      integer x_max(n)
      integer x_min(n)

      if ( .not. more ) then

        do j = 1, n
          x(j) = x_min(j)
        end do

        constraint = 0.0D+00
        do j = 1, n
          constraint = constraint + dble ( x(j) ) / dble ( x_max(j) )
        end do

        if ( 1.0D+00 .lt. constraint ) then
          more = .false.
        else
          more = .true.
        end if

        return

      else

        i = 1

10      continue

          if ( x(i) .lt. x_max(i) ) then

            x(i) = x(i) + 1

            constraint = 0.0D+00
            do j = 1, n
              constraint = constraint 
     &          + dble ( x(j) ) / dble ( x_max(j) )
            end do

            if ( constraint .le. 1.0D+00 ) then
              go to 20
            end if

          end if

          x(i) = x_min(i)

          i = i + 1

          if ( n .lt. i ) then
            more = .false.
            go to 20
          end if

        go to 10

20      continue

      end if

      return
      end
      subroutine vector_constrained_next4 ( n, alpha, x_min, x_max, 
     &  x, q, more )

c*********************************************************************72
c
cc VECTOR_CONSTRAINED_NEXT4 returns the "next" constrained vector.
c
c  Discussion:
c
c    This routine is similar to VECTOR_CONSTRAINED_NEXT2 and 
c    VECTOR_CONSTRAINED_NEXT3.
c
c    We consider all vectors X of dimension N whose components
c    satisfy X_MIN(1:N) .le. X(1:N) .le. X_MAX(1:N).
c
c    We are only interested in the subset of these vectors which
c    satisfy the following constraint:
c
c      sum ( 1 .le. I .le. N ) ALPHA(I) * X(I) .le. Q
c
c    This routine returns, one at a time, and in right-to-left
c    lexicographic order, exactly those vectors which satisfy
c    the constraint.
c
c  Example:
c
c    N = 3
c    ALPHA    4.0  3.0  5.0
c    Q       20.0
c    X_MIN:   1   1   1
c    X_MAX:   5   6   4
c
c    #  X(1)  X(2)  X(3)     Total
c
c    1    1     1     1       12.0
c    2    2     1     1       20.0
c    3    1     2     1       15.0
c    4    2     2     1       19.0
c    5    1     3     1       18.0
c    6    1     1     2       17.0
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    26 January 2007
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the number of components in the vector.
c
c    Input, double precision ALPHA(N), the coefficient vector.
c
c    Input, integer X_MIN(N), X_MAX(N), the minimum and maximum
c    values allowed in each component.
c
c    Input/output, integer X(N).  On first call (with MORE = FALSE),
c    the input value of X is not important.  On subsequent calls, the
c    input value of X should be the output value from the previous call.
c    On output, (with MORE = TRUE), the value of X will be the "next"
c    vector in the reverse lexicographical list of vectors that satisfy
c    the condition.  However, on output with MORE = FALSE, the vector
c    X is meaningless, because there are no more vectors in the list.
c
c    Input, double precision Q, the limit on the sum.
c
c    Input/output, logical MORE.  On input, if the user has set MORE
c    FALSE, the user is requesting the initiation of a new sequence
c    of values.  If MORE is TRUE, then the user is requesting "more"
c    values in the current sequence.  On output, if MORE is TRUE,
c    then another value was found and returned in X, but if MORE is
c    FALSE, then there are no more values in the sequence, and X is
c    NOT the next value.
c
      implicit none

      integer n

      double precision alpha(n)
      integer i
      integer j
      logical more
      double precision q
      double precision total
      integer x(n)
      integer x_max(n)
      integer x_min(n)

      if ( .not. more ) then

        do j = 1, n
          x(j) = x_min(j)
        end do

        total = 0.0D+00
        do j = 1, n
          total = total + alpha(j) * x(j)
        end do

        if ( q .lt. total ) then
          more = .false.
        else
          more = .true.
        end if

        return

      else

        i = 1

10      continue

          if ( x(i) .lt. x_max(i) ) then

            x(i) = x(i) + 1

            total = 0.0D+00
            do j = 1, n
              total = total + alpha(j) * x(j)
            end do

            if ( total .le. q ) then
              go to 20
            end if

          end if

          x(i) = x_min(i)

          i = i + 1

          if ( n .lt. i ) then
            more = .false.
            go to 20
          end if

        go to 10

20      continue

      end if

      return
      end
      subroutine vector_constrained_next5 ( n, x, sum_min, sum_max, 
     &  more )

c*********************************************************************72
c
cc VECTOR_CONSTRAINED_NEXT5 returns the "next" constrained vector.
c
c  Discussion:
c
c    We consider all positive integer vectors of dimension N whose 
c    components satisfy SUM_MIN .le. X(1:N) .le. SUM_MAX.
c
c    This routine returns, one at a time, and in right-to-left
c    lexicographic order, exactly those vectors which satisfy
c    the constraint.
c
c  Example:
c
c    N = 3
c    SUM_MIN = 5
c    SUM_MAX = 6
c
c    #  X(1)  X(2)  X(3)     SUM
c
c    1    3     1     1        5
c    2    2     2     1        5
c    3    2     1     2        5
c    4    1     3     1        5
c    5    1     2     2        5
c    6    1     1     3        5
c
c    7    4     1     1        6
c    8    3     2     1        6
c    9    3     1     2        6
c   10    2     3     1        6
c   11    2     2     2        6
c   12    2     1     3        6
c   13    1     4     1        6
c   14    1     3     2        6
c   15    1     2     3        6
c   16    1     1     4        6
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    26 January 2007
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the number of components in the vector.
c
c    Input, integer SUM_MIN, SUM_MAX, the minimum and maximum sums..
c
c    Input/output, integer X(N).  On first call (with MORE = FALSE), 
c    the input value of X is not important.  On subsequent calls, the
c    input value of X should be the output value from the previous call.
c    On output, (with MORE = TRUE), the value of X will be the "next"
c    vector in the reverse lexicographical list of vectors that satisfy
c    the condition.  However, on output with MORE = FALSE, the vector
c    X is meaningless, because there are no more vectors in the list.
c
c    Input/output, logical MORE.  On input, if the user has set MORE
c    FALSE, the user is requesting the initiation of a new sequence
c    of values.  If MORE is TRUE, then the user is requesting "more"
c    values in the current sequence.  On output, if MORE is TRUE,
c    then another value was found and returned in X, but if MORE is
c    FALSE, then there are no more values in the sequence, and X is
c    NOT the next value.
c
      implicit none

      integer n

      integer base
      integer i
      integer j
      logical more
      integer sum_max
      integer sum_min
      integer x(n)

      save base

      data base / 0 /
c
c  Initialization.
c
      if ( .not. more ) then

        if ( sum_max .lt. n ) then
          more = .false.
          return
        end if

        if ( sum_max .lt. sum_min ) then
          more = .false.
          return
        end if

        more = .true.

        base = max ( sum_min, n )

        x(1) = base - n + 1
        do i = 2, n
          x(i) = 1
        end do

        return
c
c  Next element.
c
      else
c
c  Search from the right, seeking an index I .lt. N for which 1 .lt. X(I).
c
        do i = n-1, 1, -1
c
c  If you find such an I, decrease X(I) by 1, and add that to X(I+1).
c
          if ( 1 .lt. x(i) ) then

            x(i)   = x(i)   - 1
            x(i+1) = x(i+1) + 1
c
c  Now grab all the "excess" 1's from the entries to the right of X(I+1).
c
            do j = i+2, n
              if ( 1 .lt. x(j) ) then
                x(i+1) = x(i+1) + x(j) - 1
                x(j) = 1
              end if
            end do

            return

          end if

        end do
c
c  The current vector is (1,1,1,...,BASE-N+1).
c  If BASE .lt. SUM_MAX, then increase BASE by 1, and start the new series.
c
        if ( base .lt. sum_max ) then
          base = base + 1
          x(1) = base - n + 1
          do i = 2, n
            x(i) = 1
          end do
          return
        end if
c
c  We returned the last legal vector on the previous call.
c  The calculation is done.
c
        more = .false.

      end if

      return
      end
      subroutine vector_constrained_next6 ( n, alpha, x_min, x_max, 
     &  x, q_min, q_max, more )

c*********************************************************************72
c
cc VECTOR_CONSTRAINED_NEXT6 returns the "next" constrained vector.
c
c  Discussion:
c
c    This routine is similar to VECTOR_CONSTRAINED_NEXT2,
c    VECTOR_CONSTRAINED_NEXT3, and VECTOR_CONSTRAINED_NEXT4.
c
c    We consider all vectors X of dimension N whose components
c    satisfy X_MIN(1:N) <= X(1:N) <= X_MAX(1:N).
c
c    We are only interested in the subset of these vectors which
c    satisfy the following constraint:
c
c      Q_MIN <= sum ( 1 <= I <= N ) ALPHA(I) * X(I) <= Q_MAX
c
c    This routine returns, one at a time, and in right-to-left
c    lexicographic order, exactly those vectors which satisfy
c    the constraint.
c
c  Example:
c
c    N = 3
c    ALPHA    4.0  3.0  5.0
c    Q_MIN   16.0
c    Q_MAX   20.0
c    X_MIN:   1   1   1
c    X_MAX:   5   6   4
c
c    #  X(1)  X(2)  X(3)     Total
c
c    1    2     1     1       20.0
c    2    2     2     1       19.0
c    3    1     3     1       18.0
c    4    1     1     2       17.0
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    08 February 2007
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the number of components in the vector.
c
c    Input, double precision ALPHA(N), the coefficient vector.
c
c    Input, integer X_MIN(N), X_MAX(N), the minimum and maximum
c    values allowed in each component.
c
c    Input/output, integer X(N).  On first call (with MORE = FALSE),
c    the input value of X is not important.  On subsequent calls, the
c    input value of X should be the output value from the previous call.
c    On output, (with MORE = TRUE), the value of X will be the "next"
c    vector in the reverse lexicographical list of vectors that satisfy
c    the condition.  However, on output with MORE = FALSE, the vector
c    X is meaningless, because there are no more vectors in the list.
c
c    Input, double precision Q_MIN, Q_MAX, the lower and upper
c    limits on the sum.
c
c    Input/output, logical MORE.  On input, if the user has set MORE
c    FALSE, the user is requesting the initiation of a new sequence
c    of values.  If MORE is TRUE, then the user is requesting "more"
c    values in the current sequence.  On output, if MORE is TRUE,
c    then another value was found and returned in X, but if MORE is
c    FALSE, then there are no more values in the sequence, and X is
c    NOT the next value.
c
      implicit none

      integer n

      double precision alpha(n)
      integer i
      integer j
      logical more
      double precision q_max
      double precision q_min
      double precision total
      integer x(n)
      integer x_max(n)
      integer x_min(n)

      if ( .not. more ) then

        more = .true.

        do i = 1, n
          x(i) = x_min(i)
        end do

        total = 0.0D+00
        do i = 1, n
          total = total + alpha(i) * dble ( x(i) )
        end do

        if ( q_min .le. total .and. total .le. q_max ) then
          return
        end if

      end if

10    continue

        j = n

20      continue

          if ( x(j) .lt. x_max(j) ) then
            go to 30
          end if

          if ( j <= 1 ) then
            more = .false.
            return
          end if

          j = j - 1

        go to 20

30      continue

        x(j) = x(j) + 1
        do i = j + 1, n
          x(i) = x_min(i)
        end do

        total = 0.0D+00
        do i = 1, n
          total = total + alpha(i) * dble ( x(i) )
        end do

        if ( q_min .le. total .and. total .le. q_max ) then
          go to 40
        end if

      go to 10

40    continue

      return
      end
      subroutine vector_constrained_next7 ( n, level_weight, x_max, x, 
     & q_min, q_max, more )

c*********************************************************************72
c
cc VECTOR_CONSTRAINED_NEXT7 returns the "next" constrained vector.
c
c  Discussion:
c
c    We consider vectors X of dimension N satisfying:
c
c      0 <= X(1:N) <= X_MAX(1:N).
c
c    and the following constraint:
c
c      Q_MIN < sum ( 1 <= I <= N ) LEVEL_WEIGHT(I) * X(I) <= Q_MAX
c
c    This routine returns, one at a time, and in right-to-left
c    lexicographic order, exactly those vectors which satisfy
c    the constraint.
c
c  Example:
c
c    N = 3
c    LEVEL_WEIGHT    4.0  3.0  5.0
c    Q_MIN   16.0
c    Q_MAX   20.0
c    X_MAX:   5   6   4
c
c    #  X(1)  X(2)  X(3)     Total
c
c    1    3     1     1       20.0
c    2    2     2     1       19.0
c    3    1     3     1       18.0
c    4    1     1     2       17.0
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    01 July 2010
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the number of components in the vector.
c
c    Input, double precision LEVEL_WEIGHT(N), the coefficient vector.
c
c    Input, integer X_MAX(N), the maximum
c    values allowed in each component.
c
c    Input/output, integer X(N).  On first call, with 
c    MORE = FALSE, the input value of X is not important.  On subsequent calls, 
c    the input value of X should be the output value from the previous call.
c    On output, (with MORE = TRUE), the value of X will be the "next"
c    vector in the reverse lexicographical list of vectors that satisfy
c    the condition.  However, on output with MORE = FALSE, the vector
c    X is meaningless, because there are no more vectors in the list.
c
c    Input, double precision Q_MIN, Q_MAX, the lower and upper
c    limits on the sum.
c
c    Input/output, logical MORE.  On input, if the user has set MORE
c    FALSE, the user is requesting the initiation of a new sequence
c    of values.  If MORE is TRUE, then the user is requesting "more"
c    values in the current sequence.  On output, if MORE is TRUE,
c    then another value was found and returned in X, but if MORE is
c    FALSE, then there are no more values in the sequence, and X is
c    NOT the next value.
c
      implicit none

      integer n

      integer i
      integer j
      double precision level_weight(n)
      logical more
      double precision q_max
      double precision q_min
      double precision total
      integer x(n)
      integer x_max(n)

      if ( .not. more ) then

        more = .true.

        do j = 1, n
          x(j) = 0
        end do

        total = 0.0D+00
        do j = 1, n
          total = total + level_weight(j) * dble ( x(j) )
        end do

        if ( q_min .lt. total .and. total .le. q_max ) then
          return
        end if

      end if

10    continue

        i = n

20      continue

          if ( x(i) .lt. x_max(i) ) then
            go to 30
          end if

          if ( i .le. 1 ) then
            more = .false.
            return
          end if

          i = i - 1

        go to 20

30      continue

        x(i) = x(i) + 1
        do j = i + 1, n
          x(j) = 0
        end do

        total = 0.0D+00
        do j = 1, n
          total = total + level_weight(j) * dble ( x(j) )
        end do

        if ( q_min .lt. total .and. total .le. q_max ) then
          go to 40
        end if

      go to 10

40    continue

      return
      end
      subroutine vector_next ( n, x_min, x_max, x, more )

c*********************************************************************72
c
cc VECTOR_NEXT returns the "next" vector between two ranges.
c
c  Discussion:
c
c    We consider all integer vectors of dimension N satisfying:
c
c      X_MIN(1:N) <= X(1:N) <= X_MAX(1:N).
c
c    This routine returns, one at a time, and in right-to-left
c    lexicographic order, all these vectors.
c
c  Example:
c
c    N = 3
c    X_MIN:   2   2   0
c    X_MAX:   4   3   1
c 
c    #  X(1)  X(2)  X(3)
c
c    1    2     2     0
c    2    3     2     0
c    3    4     2     0
c    4    2     3     0
c    5    3     3     0
c    6    4     3     0
c    7    2     2     1
c    8    3     2     1
c    9    4     2     1
c   10    2     3     1
c   11    3     3     1
c   12    4     3     1
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    01 July 2010
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the number of components in the vector.
c
c    Input, integer X_MIN(N), X_MAX(N), the minimum and maximum
c    values allowed in each component.
c
c    Input/output, integer X(N).  On first call, with 
c    MORE = FALSE, the input value of X is not important.  On subsequent calls,
c    the input value of X should be the output value from the previous call.
c    On output, with MORE = TRUE, the value of X will be the "next"
c    vector in the reverse lexicographical list of vectors.  However, on 
c    output with MORE = FALSE, the vector X is meaningless, because there 
c    are no more vectors in the list.
c
c    Input/output, logical MORE.  On input, if the user has set MORE
c    FALSE, the user is requesting the initiation of a new sequence
c    of values.  If MORE is TRUE, then the user is requesting "more"
c    values in the current sequence.  On output, if MORE is TRUE,
c    then another value was found and returned in X, but if MORE is
c    FALSE, then there are no more values in the sequence, and X is
c    NOT the next value.
c
      implicit none

      integer n

      integer i
      logical more
      integer x(n)
      integer x_max(n)
      integer x_min(n)

      if ( .not. more ) then

        do i = 1, n
          x(i) = x_min(i)
        end do
        more = .true.

      else

        i = 1

10      continue

          if ( x(i) .lt. x_max(i) ) then
            x(i) = x(i) + 1
            go to 20
          end if

          x(i) = x_min(i)

          i = i + 1

          if ( n .lt. i ) then
            more = .false.
            go to 20
          end if

        go to 10

20      continue

      end if

      return
      end
      subroutine ytb_enum ( n, ytb_num )

c*********************************************************************72
c
cc YTB_ENUM enumerates the Young tables of size N.
c
c  Discussion:
c
c    If A(N) is the number of Young table of size N, then A(1) = 1,
c    A(2) = 2, and
c
c    A(N) = A(N-1) + (N-1) * A(N-2).
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    19 January 2007
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the integer which is to be partitioned.
c
c    Output, integer YTB_NUM, the number of Young table of N.
c
      implicit none

      integer a1
      integer a2
      integer a3
      integer i
      integer n
      integer ytb_num

      if ( n .le. 0 ) then
        ytb_num = 0
      else if ( n .eq. 1 ) then
        ytb_num = 1
      else if ( n .eq. 2 ) then
        ytb_num = 2
      else
        a2 = 1
        a3 = 2
        do i = 3, n
          a1 = a2
          a2 = a3
          a3 = a2 + ( i - 1 ) * a1
        end do
        ytb_num = a3
      end if

      return
      end
      subroutine ytb_next ( n, lambda, a, more )

c*********************************************************************72
c
cc YTB_NEXT computes the next Young table for a given shape.
c
c  Discussion:
c
c    When the routine is called with MORE = FALSE (the first time), and
c    if LAMBDA on this call has M parts, with M < N, then the user
c    must also make sure that LAMBDA(M+1) = 0.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    23 January 2007
c
c  Author:
c
c    Original FORTRAN77 version by Albert Nijenhuis, Herbert Wilf.
c    This FORTRAN77 version by John Burkardt.
c
c  Reference:
c
c    Albert Nijenhuis, Herbert Wilf,
c    Combinatorial Algorithms for Computers and Calculators,
c    Second Edition,
c    Academic Press, 1978,
c    ISBN: 0-12-519260-6,
c    LC: QA164.N54.
c
c  Parameters:
c
c    Input, integer N, the integer which is to be partitioned.
c
c    Input/output, integer LAMBDA(N), contains a partition of N. 
c    The elements of LAMBDA are nonnegative integers that sum to N. 
c    On the first call, with MORE = FALSE, the user sets LAMBDA.
c    After the first call, the input value of LAMBDA is not important.
c    On output, the value of LAMBDA is the partition corresponding
c    to the Young table.
c
c    Input/output, integer A(N).  On the first call, with MORE = FALSE,
c    no value of A needs to be set.  After the first call, the input
c    value of A should be its output value from the previous call.
c    The output value of A is the next Young table.  A(I) is the 
c    row containing I in the output table.
c
c    Input/output, logical MORE.  Set MORE to FALSE before the first call.
c    It is reset to TRUE as the program returns a new table
c    on each call, until the last table is computed, when
c    the program also sets MORE = FALSE.
c
      implicit none

      integer n

      integer a(n)
      integer i
      integer ir
      integer it
      integer j
      integer k
      logical more
      integer lambda(n)
      integer isave

      it = n

      if ( more ) then

        lambda(1) = 1
        do i = 2, n
          lambda(i) = 0
        end do

        isave = 0

        do i = 2, n

          lambda(a(i)) = lambda(a(i)) + 1

          if ( a(i) .lt. a(i-1) ) then
            isave = i
            go to 10
          end if

        end do

10      continue

        if ( isave == 0 ) then
          more = .false.
          return
        end if

        it = lambda(1+a(isave))

        do i = n, 1, -1

          if ( lambda(i) .eq. it ) then
            a(isave) = i
            lambda(i) = lambda(i) - 1
            it = isave - 1
            go to 20
          end if

        end do

20      continue

      end if

      k = 1
      ir = 1

30    continue

        if ( n .lt. ir ) then
          go to 40
        end if

        if ( lambda(ir) .ne. 0 ) then
          a(k) = ir
          lambda(ir) = lambda(ir) - 1
          k = k + 1
          ir = ir + 1
          go to 30
        end if

        if ( it .lt. k ) then
          go to 40
        end if

        ir = 1

      go to 30

40    continue

      if ( n .eq. 1 ) then
        more = .false.
        return
      end if

      do j = 2, n
        if ( a(j) .lt. a(j-1) ) then
          more = .true.
          return
        end if
      end do

      more = .false.

      return
      end
      subroutine ytb_print ( n, a, title )

c*********************************************************************72
c
cc YTB_PRINT prints a Young table.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    23 January 2007
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the integer that is partitioned.
c
c    Input, integer A(N), describes the Young table.
c    A(I) is the row of the table on which I occurs.
c
c    Input, character * ( * ) TITLE, an optional title.
c
      implicit none

      integer n

      integer a(n)
      integer j
      integer jarray(n)
      integer row_i
      integer row_length
      character * ( * ) title
      integer title_length

      title_length = len_trim ( title )

      if ( 0 .lt. title_length ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) title(1:title_length)
      end if

      write ( *, '(a)' ) ' '

      row_i = 0

10    continue

        row_i = row_i + 1

        row_length = 0

        do j = 1, n

          if ( a(j) == row_i ) then
            row_length = row_length + 1
            jarray(row_length) = j
          end if

        end do

        if ( row_length <= 0 ) then
          go to 20
        end if

        write ( *, '(20i4)' ) ( jarray(j), j = 1, row_length )

      go to 10

20    continue

      return
      end
      subroutine ytb_random ( n, lambda, seed, a )

c*********************************************************************72
c
cc YTB_RANDOM selects a random Young table of a given shape.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    23 January 2007
c
c  Author:
c
c    Original FORTRAN77 version by Albert Nijenhuis, Herbert Wilf.
c    This FORTRAN77 version by John Burkardt.
c
c  Reference:
c
c    Albert Nijenhuis, Herbert Wilf,
c    Combinatorial Algorithms for Computers and Calculators,
c    Second Edition,
c    Academic Press, 1978,
c    ISBN: 0-12-519260-6,
c    LC: QA164.N54.
c
c  Parameters:
c
c    Input, integer N, the integer which has been partitioned.
c
c    Input, integer LAMBDA(N), the partition of N.
c    N = sum ( 1 <= I <= N ) LAMBDA(I).
c
c    Input/output, integer SEED, a seed for the random number generator.
c
c    Output, integer A(N), the vector describing the Young table.
c
      implicit none

      integer n

      integer a(n)
      integer i
      integer i4_uniform
      integer ih
      integer j
      integer k
      integer lambda(n)
      integer m
      integer seed

      do j = 1, n
        a(j) = 0
      end do

      i = 0
      k = 0

10    continue

        i = i + 1
        do j = 1, lambda(i)
          a(j) = a(j) + 1
          k = k + 1
        end do

        if ( n .le. k ) then
          go to 20
        end if

      go to 10

20    continue

      do m = 1, n

30      continue

          i = i4_uniform ( 1, a(1), seed )
          j = i4_uniform ( 1, lambda(1), seed )

          if ( i .le. a(j) .and. j .le. lambda(i) ) then
            go to 40
          end if

        go to 30

40      continue

          ih = a(j) + lambda(i) - i - j

          if ( ih .eq. 0 ) then
            go to 50
          end if

          k = i4_uniform ( 1, ih, seed )

          if ( k .le. lambda(i)-j ) then
            j = j + k
          else
            i = k - lambda(i) + i + j
          end if

        go to 40

50      continue

        lambda(i) = lambda(i) - 1
        a(j) = a(j) - 1
        a(n+1-m) = i

      end do

      do i = 1, n
        lambda(a(i)) = lambda(a(i)) + 1
      end do

      return
      end
