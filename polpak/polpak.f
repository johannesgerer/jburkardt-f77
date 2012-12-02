      function agm ( a, b )

c*********************************************************************72
c
cc AGM computes the arithmetic-geometric mean of A and B.
c
c  Discussion:
c
c    The AGM is defined for nonnegative A and B.
c
c    The AGM of numbers A and B is defined by setting
c
c      A(0) = A,
c      B(0) = B
c
c      A(N+1) = ( A(N) + B(N) ) / 2
c      B(N+1) = sqrt ( A(N) * B(N) )
c
c    The two sequences both converge to AGM(A,B).
c
c    In Mathematica, the AGM can be evaluated by
c
c      ArithmeticGeometricMean [ a, b ]
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    09 February 2008
c
c  Author:
c
c    John Burkardt
c
c  Reference:
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
c    Input, double precision A, B, the arguments whose AGM is to be computed.
c
c    Output, double precision AGM, the arithmetic-geometric mean of A and B.
c
      implicit none

      double precision a
      double precision agm
      double precision a2
      double precision b
      double precision b2
      double precision c
      double precision d
      integer it
      integer it_max
      parameter ( it_max = 1000 )
      double precision r8_epsilon
      double precision tol

      if ( a .lt. 0.0D+00 ) then
        stop
      end if

      if ( b .lt. 0.0D+00 ) then
        stop
      end if

      if ( a .eq. 0.0D+00 .or. b .eq. 0.0D+00 ) then
        agm = 0.0D+00
        return
      end if

      it = 0
      tol = 100.0D+00 * r8_epsilon ( )

      a2 = a
      b2 = b

10    continue

        it = it + 1

        c = ( a2 + b2 ) / 2.0D+00
        d = sqrt ( a2 * b2 )

        if ( abs ( c - d ) .le. tol * ( c + d ) ) then
          go to 20
        end if

        if ( it_max .lt. it ) then
          go to 20
        end if

        a2 = c
        b2 = d

      go to 10

20    continue

      agm = c

      return
      end
      subroutine agm_values ( n_data, a, b, fx )

c*********************************************************************72
c
cc AGM_VALUES returns some values of the arithmetic geometric mean.
c
c  Discussion:
c
c    The AGM is defined for nonnegative A and B.
c
c    The AGM of numbers A and B is defined by setting
c
c      A(0) = A,
c      B(0) = B
c
c      A(N+1) = ( A(N) + B(N) ) / 2
c      B(N+1) = sqrt ( A(N) * B(N) )
c
c    The two sequences both converge to AGM(A,B).
c
c    In Mathematica, the AGM can be evaluated by
c
c      ArithmeticGeometricMean [ a, b ]
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    09 February 2008
c
c  Author:
c
c    John Burkardt
c
c  Reference:
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
c    Output, double precision A, B, the numbers whose AGM is desired.
c
c    Output, double precision FX, the value of the function.
c
      implicit none

      integer n_max
      parameter ( n_max = 15 )

      double precision a
      double precision a_vec(n_max) 
      double precision b
      double precision b_vec(n_max) 
      double precision fx
      double precision fx_vec(n_max) 
      integer n_data

      save a_vec
      save b_vec
      save fx_vec

      data a_vec /
     &   22.0D+00, 
     &   83.0D+00, 
     &   42.0D+00, 
     &   26.0D+00, 
     &    4.0D+00, 
     &    6.0D+00, 
     &   40.0D+00,
     &   80.0D+00,
     &   90.0D+00,
     &    9.0D+00,
     &   53.0D+00,
     &    1.0D+00,
     &    1.0D+00,
     &    1.0D+00,
     &    1.5D+00 /
      data b_vec /
     &   96.0D+00,
     &   56.0D+00,
     &    7.0D+00,
     &   11.0D+00,
     &   63.0D+00,
     &   45.0D+00,
     &   75.0D+00,
     &    0.0D+00,
     &   35.0D+00,
     &    1.0D+00,
     &   53.0D+00,
     &    2.0D+00,
     &    4.0D+00,
     &    8.0D+00,
     &    8.0D+00 /
      data fx_vec / 
     &   52.274641198704240049D+00,
     &   68.836530059858524345D+00,
     &   20.659301196734009322D+00,
     &   17.696854873743648823D+00,
     &   23.867049721753300163D+00,
     &   20.717015982805991662D+00,
     &   56.127842255616681863D+00,
     &    0.000000000000000000D+00,
     &   59.269565081229636528D+00,
     &   3.9362355036495554780D+00,
     &   53.000000000000000000D+00,
     &   1.4567910310469068692D+00,
     &   2.2430285802876025701D+00,
     &   3.6157561775973627487D+00,
     &   4.0816924080221632670D+00 /

      if ( n_data .lt. 0 ) then
        n_data = 0
      end if

      n_data = n_data + 1

      if ( n_max .lt. n_data ) then
        n_data = 0
        a = 0.0D+00
        b = 0.0D+00
        fx = 0.0D+00
      else
        a = a_vec(n_data)
        b = b_vec(n_data)
        fx = fx_vec(n_data)
      end if

      return
      end
      function agud ( gamma )

c*********************************************************************72
c
cc AGUD evaluates the inverse Gudermannian function.
c
c  Discussion:
c
c    The Gudermannian function relates the hyperbolic and trigonometric
c    functions.  For any argument X, there is a corresponding value
c    GAMMA so that
c
c      SINH(X) = TAN(GAMMA).
c
c    This value GAMMA(X) is called the Gudermannian of X.  The inverse
c    Gudermannian function is given as input a value GAMMA and computes
c    the corresponding value X.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    01 December 2007
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, double precision GAMMA, the value of the Gudermannian.
c
c    Output, double precision AGUD, the argument of the Gudermannian.
c
      implicit none

      double precision agud
      double precision gamma
      double precision pi
      parameter ( pi = 3.141592653589793D+00 )

      agud = log ( tan ( 0.25D+00 * pi + 0.5D+00 * gamma ) )

      return
      end
      function align_enum ( m, n )

c*********************************************************************72
c
cc ALIGN_ENUM counts the alignments of two sequences of M and N elements.
c
c  Discussion:
c
c    We assume that we have sequences A and B of M and N characters each.
c    An alignment of the two sequences is a rule matching corresponding
c    elements of one sequence to another.  Some elements of either sequence
c    can be matched to a null element.  If A(I1) and A(I2) are matched
c    to B(J1) and B(J2), and I1 < I2, then it must be the case that J1 < J2.
c
c    The 5 alignments of a sequence of 1 to a sequence of 2 are:
c
c          _1_   _2_   __3__   __4__   __5__
c
c      A:  1 -   - 1   - 1 -   - - 1   1 - -
c      B:  1 2   1 2   1 - 2   1 2 -   - 1 2
c
c    The formula is:
c
c      F(0,0) = 1
c      F(1,0) = 1
c      F(0,1) = 1
c      F(M,N) = F(M-1,N) + F(M-1,N-1) + F(M,N-1)
c
c    To compute F(M,N), it is not necessary to keep an M+1 by N+1
c    array in memory.  A vector of length N will do.
c
c    F(N,N) is approximately ( 1 + sqrt(2) )^(2*N+1) / sqrt ( N )
c
c  Example:
c
c    The initial portion of the table is:
c
c
c  M/N   0    1    2    3    4       5       6       7       8       9      10
c
c   0    1    1    1    1    1       1       1       1       1       1       1
c   1    1    3    5    7    9      11      13      15      17      19      21
c   2    1    5   13   25   41      61      85     113     145     181     221
c   3    1    7   25   63  129     231     377     575     833    1159    1561
c   4    1    9   41  129  321     681    1289    2241    3649    5641    8361
c   5    1   11   61  231  681    1683    3653    7183   13073   22363   36365
c   6    1   13   85  377 1289    3653    8989   19825   40081   75517  134245
c   7    1   15  113  575 2241    7183   19825   48639  108545  224143  433905
c   8    1   17  145  833 3649   13073   40081  108545  265729  598417 1256465
c   9    1   19  181 1159 5641   22363   75517  224143  598417 1462563 3317445
c  10    1   21  221 1561 8361   36365  134245  433905 1256465 3317445 8097453
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    01 December 2007
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Michael Waterman,
c    Introduction to Computational Biology,
c    Chapman and Hall, 1995, pages 186-190.
c
c  Parameters:
c
c    Input, integer M, N, the number of elements of the
c    two sequences.
c
c    Output, integer ALIGN_ENUM, the number of possible
c    alignments of the sequences.
c
      implicit none

      integer n

      integer align_enum
      integer fi(0:n)
      integer fim1j
      integer fim1jm1
      integer i
      integer j
      integer m

      if ( m .lt. 0 ) then
        align_enum = 0
        return
      else if ( n .lt. 0 ) then
        align_enum = 0
        return
      else if ( m .eq. 0 ) then
        align_enum = 1
        return
      else if ( n .eq. 0 ) then
        align_enum = 1
        return
      end if

      fi(0:n) = 1

      do i = 1, m

        fim1jm1 = 1

        do j = 1, n

          fim1j = fi(j)

          fi(j) = fi(j) + fi(j-1) + fim1jm1

          fim1jm1 = fim1j

        end do
      end do

      align_enum = fi(n)

      return
      end
      function arc_cosine ( c )

c*********************************************************************72
c
cc ARC_COSINE computes the arc cosine function, with argument truncation.
c
c  Discussion:
c
c    If you call your system ACOS routine with an input argument that is
c    outside the range [-1.0, 1.0 ], you may get an unpleasant surprise.
c
c    In particular, you may get the value NaN returned.
c
c    This routine truncates arguments outside the range, avoiding the problem.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    01 December 2007
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, double precision C, the argument.
c
c    Output, double precision ARC_COSINE, an angle whose cosine is C.
c
      implicit none

      double precision arc_cosine
      double precision c
      double precision c2

      c2 = c
      c2 = max ( c2, -1.0D+00 )
      c2 = min ( c2, +1.0D+00 )

      arc_cosine = acos ( c2 )

      return
      end
      function arc_sine ( s )

c*********************************************************************72
c
cc ARC_SINE computes the arc sine function, with argument truncation.
c
c  Discussion:
c
c    If you call your system ASIN routine with an input argument that is
c    outside the range [-1.0, 1.0 ], you may get an unpleasant surprise.
c
c    In particular, you may get the value NaN returned.
c
c    This routine truncates arguments outside the range, avoiding the problem.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    03 July 2008
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, double precision S, the argument.
c
c    Output, double precision ARC_SINE, an angle whose sine is S.
c
      implicit none

      double precision arc_sine
      double precision s
      double precision s2

      s2 = s
      s2 = max ( s2, -1.0D+00 )
      s2 = min ( s2, +1.0D+00 )

      arc_sine = asin ( s2 )

      return
      end
      function atan4 ( y, x )

c*********************************************************************72
c
cc ATAN4 computes the inverse tangent of the ratio Y / X.
c
c  Discussion:
c
c    ATAN4 returns an angle whose tangent is ( Y / X ), a job which
c    the built in functions ATAN and ATAN2 already do.
c
c    However:
c
c    * ATAN4 always returns a positive angle, between 0 and 2 PI,
c      while ATAN and ATAN2 return angles in the interval [-PI/2,+PI/2]
c      and [-PI,+PI] respectively;
c
c    * ATAN4 accounts for the signs of X and Y, (as does ATAN2).  The ATAN
c     function by contrast always returns an angle in the first or fourth
c     quadrants.
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
c    Input, double precision Y, X, two quantities which represent the
c    tangent of an angle.  If Y is not zero, then the tangent is (Y/X).
c
c    Output, double precision ATAN4, an angle between 0 and 2 * PI, whose
c    tangent is (Y/X), and which lies in the appropriate quadrant so that
c    the signs of its cosine and sine match those of X and Y.
c
      implicit none

      double precision abs_x
      double precision abs_y
      double precision atan4
      double precision pi
      parameter ( pi = 3.141592653589793D+00 )
      double precision theta
      double precision theta_0
      double precision x
      double precision y
c
c  Special cases:
c
      if ( x .eq. 0.0D+00 ) then

        if ( 0.0D+00 .lt. y ) then
          theta = pi / 2.0D+00
        else if ( y .lt. 0.0D+00 ) then
          theta = 3.0D+00 * pi / 2.0D+00
        else if ( y .eq. 0.0D+00 ) then
          theta = 0.0D+00
        end if

      else if ( y .eq. 0.0D+00 ) then

        if ( 0.0D+00 .lt. x ) then
          theta = 0.0D+00
        else if ( x .lt. 0.0D+00 ) then
          theta = pi
        end if
c
c  We assume that ATAN2 is correct when both arguments are positive.
c
      else

        abs_y = dabs ( y )
        abs_x = dabs ( x )

        theta_0 = atan2 ( abs_y, abs_x )

        if ( 0.0D+00 .lt. x .and. 0.0D+00 .lt. y ) then
          theta = theta_0
        else if ( x .lt. 0.0D+00 .and. 0.0D+00 .lt. y ) then
          theta = pi - theta_0
        else if ( x .lt. 0.0D+00 .and. y .lt. 0.0D+00 ) then
          theta = pi + theta_0
        else if ( 0.0D+00 .lt. x .and. y .lt. 0.0D+00 ) then
          theta = 2.0D+00 * pi - theta_0
        end if

      end if

      atan4 = theta

      return
      end
      subroutine bell ( n, b )

c*********************************************************************72
c
cc BELL returns the Bell numbers from 0 to N.
c
c  Discussion:
c
c    The Bell number B(N) is the number of restricted growth functions on N.
c
c    Note that the Stirling numbers of the second kind, S^m_n, count the
c    number of partitions of N objects into M classes, and so it is
c    true that
c
c      B(N) = S^1_N + S^2_N + ... + S^N_N.
c
c    The Bell numbers were named for Eric Temple Bell.
c
c  Definition:
c
c    The Bell number B(N) is defined as the number of partitions (of
c    any size) of a set of N distinguishable objects.
c
c    A partition of a set is a division of the objects of the set into
c    subsets.
c
c  Examples:
c
c    There are 15 partitions of a set of 4 objects:
c
c      (1234),
c      (123) (4),
c      (124) (3),
c      (12) (34),
c      (12) (3) (4),
c      (134) (2),
c      (13) (24),
c      (13) (2) (4),
c      (14) (23),
c      (1) (234),
c      (1) (23) (4),
c      (14) (2) (3),
c      (1) (24) (3),
c      (1) (2) (34),
c      (1) (2) (3) (4).
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
c    B(I) = sum ( 1 <= J <=I ) Binomial ( I-1, J-1 ) * B(I-J)
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    01 December 2007
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
      integer i4_choose
      integer j

      if ( n .lt. 0 ) then
        return
      end if

      b(0) = 1

      do i = 1, n
        b(i) = 0
        do j = 1, i
          combo = i4_choose ( i-1, j-1 )
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
c    National Bureau of Standards, 1964,
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

      integer n_max
      parameter ( n_max = 11 )

      integer c
      integer c_vec(n_max)
      integer n
      integer n_data
      integer n_vec(n_max)

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
      function benford ( ival )

c*********************************************************************72
c
cc BENFORD returns the Benford probability of one or more significant digits.
c
c  Discussion:
c
c    Benford's law is an empirical formula explaining the observed
c    distribution of initial digits in lists culled from newspapers,
c    tax forms, stock market prices, and so on.  It predicts the observed
c    high frequency of the initial digit 1, for instance.
c
c    Note that the probabilities of digits 1 through 9 are guaranteed
c    to add up to 1, since
c      LOG10 ( 2/1 ) + LOG10 ( 3/2) + LOG10 ( 4/3 ) + ... + LOG10 ( 10/9 )
c      = LOG10 ( 2/1 * 3/2 * 4/3 * ... * 10/9 ) = LOG10 ( 10 ) = 1.
c
c    The formula is:
c
c      Prob ( First significant digits are IVAL ) =
c        LOG10 ( ( IVAL + 1 ) / IVAL ).
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    02 December 2007
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Frank Benford,
c    The Law of Anomalous Numbers,
c    Proceedings of the American Philosophical Society,
c    Volume 78, pages 551-572, 1938.
c
c    Ted Hill,
c    The First Digit Phenomenon,
c    American Scientist,
c    Volume 86, July/August 1998, pages 358 - 363.
c
c    Ralph Raimi,
c    The Peculiar Distribution of First Digits,
c    Scientific American,
c    December 1969, pages 109-119.
c
c  Parameters:
c
c    Input, integer IVAL, the string of significant digits to
c    be checked.  If IVAL is 1, then we are asking for the Benford probability
c    that a value will have first digit 1.  If IVAL is 123, we are asking for
c    the probability that the first three digits will be 123, and so on.
c    Note that IVAL must not be 0 or negative.
c
c    Output, double precision BENFORD, the Benford probability that an
c    item taken from a real world distribution will have the initial
c    digits IVAL.
c
      implicit none

      double precision benford
      integer ival

      if ( ival <= 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'BENFORD - Fatal errorc'
        write ( *, '(a)' ) '  The input argument must be positive.'
        write ( *, '(a,i8)' ) '  Your value was ', ival
        stop
      end if

      benford = log10 ( dble ( ival + 1 ) / dble ( ival ) )

      return
      end
      subroutine bernoulli_number ( n, b )

c*********************************************************************72
c
cc BERNOULLI_NUMBER computes the value of the Bernoulli numbers B(0) through B(N).
c
c  Discussion:
c
c    The Bernoulli numbers are rational.
c
c    If we define the sum of the M-th powers of the first N integers as:
c
c      SIGMA(M,N) = sum ( 0 <= I <= N ) I**M
c
c    and let C(I,J) be the combinatorial coefficient:
c
c      C(I,J) = I! / ( ( I - J )! * J! )
c
c    then the Bernoulli numbers B(J) satisfy:
c
c      SIGMA(M,N) = 1/(M+1) * sum ( 0 <= J <= M ) C(M+1,J) B(J) * (N+1)**(M+1-J)
c
c  First values:
c
c   B0  1                   =         1.00000000000
c   B1 -1/2                 =        -0.50000000000
c   B2  1/6                 =         1.66666666666
c   B3  0                   =         0
c   B4 -1/30                =        -0.03333333333
c   B5  0                   =         0
c   B6  1/42                =         0.02380952380
c   B7  0                   =         0
c   B8 -1/30                =        -0.03333333333
c   B9  0                   =         0
c  B10  5/66                =         0.07575757575
c  B11  0                   =         0
c  B12 -691/2730            =        -0.25311355311
c  B13  0                   =         0
c  B14  7/6                 =         1.16666666666
c  B15  0                   =         0
c  B16 -3617/510            =        -7.09215686274
c  B17  0                   =         0
c  B18  43867/798           =        54.97117794486
c  B19  0                   =         0
c  B20 -174611/330          =      -529.12424242424
c  B21  0                   =         0
c  B22  854,513/138         =      6192.123
c  B23  0                   =         0
c  B24 -236364091/2730      =    -86580.257
c  B25  0                   =         0
c  B26  8553103/6           =   1425517.16666
c  B27  0                   =         0
c  B28 -23749461029/870     = -27298231.0678
c  B29  0                   =         0
c  B30  8615841276005/14322 = 601580873.901
c
c  Recursion:
c
c    With C(N+1,K) denoting the standard binomial coefficient,
c
c    B(0) = 1.0
c    B(N) = - ( sum ( 0 <= K < N ) C(N+1,K) * B(K) ) / C(N+1,N)
c
c  Warning:
c
c    This recursion, which is used in this routine, rapidly results
c    in significant errors.
c
c  Special Values:
c
c    Except for B(1), all Bernoulli numbers of odd index are 0.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    03 December 2007
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the order of the highest Bernoulli
c    number to compute.
c
c    Output, double precision B(0:N), B(I) contains the I-th Bernoulli number.
c
      implicit none

      integer n

      double precision b(0:n)
      double precision b_sum
      integer i
      integer ido
      integer iwork(0:n+1)
      integer j

      if ( n .lt. 0 ) then
        return
      end if

      b(0) = 1.0D+00

      if ( n .lt. 1 ) then
        return
      end if

      b(1) = -0.5D+00

      ido = 0

      do i = 2, n

        call comb_row ( ido, i+1, iwork )
        ido = 1

        if ( mod ( i, 2 ) .eq. 1 ) then

          b(i) = 0.0D+00

        else

          b_sum = 0.0D+00
          do j = 0, i - 1
            b_sum = b_sum + b(j) * dble ( iwork(j) )
          end do

          b(i) = -b_sum / dble ( iwork(i) )

        end if

      end do

      return
      end
      subroutine bernoulli_number2 ( n, b )

c*********************************************************************72
c
cc BERNOULLI_NUMBER2 evaluates the Bernoulli numbers.
c
c  Discussion:
c
c    The Bernoulli numbers are rational.
c
c    If we define the sum of the M-th powers of the first N integers as:
c
c      SIGMA(M,N) = sum ( 0 <= I <= N ) I**M
c
c    and let C(I,J) be the combinatorial coefficient:
c
c      C(I,J) = Ic / ( ( I - J )c * Jc )
c
c    then the Bernoulli numbers B(J) satisfy:
c
c      SIGMA(M,N) = 1/(M+1) * sum ( 0 <= J <= M ) C(M+1,J) B(J) * (N+1)**(M+1-J)
c
c    Note that the Bernoulli numbers grow rapidly.  Bernoulli number
c    62 is probably the last that can be computed on the VAX without
c    overflow.
c
c    A different method than that used in BERN is employed.
c
c  First values:
c
c   B0  1                   =         1.00000000000
c   B1 -1/2                 =        -0.50000000000
c   B2  1/6                 =         1.66666666666
c   B3  0                   =         0
c   B4 -1/30                =        -0.03333333333
c   B5  0                   =         0
c   B6  1/42                =         0.02380952380
c   B7  0                   =         0
c   B8 -1/30                =        -0.03333333333
c   B9  0                   =         0
c  B10  5/66                =         0.07575757575
c  B11  0                   =         0
c  B12 -691/2730            =        -0.25311355311
c  B13  0                   =         0
c  B14  7/6                 =         1.16666666666
c  B15  0                   =         0
c  B16 -3617/510            =        -7.09215686274
c  B17  0                   =         0
c  B18  43867/798           =        54.97117794486
c  B19  0                   =         0
c  B20 -174611/330          =      -529.12424242424
c  B21  0                   =         0
c  B22  854,513/138         =      6192.123
c  B23  0                   =         0
c  B24 -236364091/2730      =    -86580.257
c  B25  0                   =         0
c  B26  8553103/6           =   1425517.16666
c  B27  0                   =         0
c  B28 -23749461029/870     = -27298231.0678
c  B29  0                   =         0
c  B30  8615841276005/14322 = 601580873.901
c
c  Recursion:
c
c    With C(N+1,K) denoting the standard binomial coefficient,
c
c    B(0) = 1.0
c    B(N) = - ( sum ( 0 <= K < N ) C(N+1,K) * B(K) ) / C(N+1,N)
c
c  Special Values:
c
c    Except for B(1), all Bernoulli numbers of odd index are 0.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    04 December 2007
c
c  Parameters:
c
c    Input, integer N, the highest order Bernoulli number
c    to compute.
c
c    Output, real B(0:N), the requested Bernoulli numbers.
c
      implicit none

      integer n

      real altpi
      real b(0:n)
      integer i
      integer k
      integer kmax
      parameter ( kmax = 400 )
      real pi
      parameter ( pi = 3.141592653589793D+00 )
      real sgn
      real sum2
      real t
      real term
      real tol
      parameter ( tol = 1.0D-06 )

      if ( n .lt. 0 ) then
        return
      end if

      b(0) = 1.0D+00

      if ( n .lt. 1 ) then
        return
      end if

      b(1) = -0.5D+00

      if ( n .lt. 2 ) then
        return
      end if

      altpi = log ( 2.0D+00 * pi )
c
c  Initial estimates for B(I), I = 2 to N
c
      b(2) = log ( 2.0D+00 )
      do i = 3, n
        if ( mod ( i, 2 ) .eq. 1 ) then
          b(i) = 0.0D+00
        else
          b(i) = log ( dble ( i * ( i - 1 ) ) ) + b(i-2)
        end if
      end do

      b(2) = 1.0D+00 / 6.0D+00

      if ( n .le. 3 ) then
        return
      end if

      b(4) = -1.0D+00 / 30.0D+00

      sgn = -1.0D+00

      do i = 6, n, 2

        sgn = -sgn
        t = 2.0D+00 * sgn * exp ( b(i) - dble ( i ) * altpi )

        sum2 = 1.0D+00

        do k = 2, kmax

          term = dble ( k )**(-i)
          sum2 = sum2 + term

          if ( term .le. tol * sum2 ) then
            exit
          end if

        end do

        b(i) = t * sum2

      end do

      return
      end
      subroutine bernoulli_number3 ( n, b )

c*********************************************************************72
c
cc BERNOULLI_NUMBER3 computes the value of the Bernoulli number B(N).
c
c  Discussion:
c
c    The Bernoulli numbers are rational.
c
c    If we define the sum of the M-th powers of the first N integers as:
c
c      SIGMA(M,N) = sum ( 0 <= I <= N ) I**M
c
c    and let C(I,J) be the combinatorial coefficient:
c
c      C(I,J) = Ic / ( ( I - J )c * Jc )
c
c    then the Bernoulli numbers B(J) satisfy:
c
c      SIGMA(M,N) = 1/(M+1) * sum ( 0 <= J <= M ) C(M+1,J) B(J) * (N+1)**(M+1-J)
c
c  First values:
c
c     B0  1                   =         1.00000000000
c     B1 -1/2                 =        -0.50000000000
c     B2  1/6                 =         1.66666666666
c     B3  0                   =         0
c     B4 -1/30                =        -0.03333333333
c     B5  0                   =         0
c     B6  1/42                =         0.02380952380
c     B7  0                   =         0
c     B8 -1/30                =        -0.03333333333
c     B9  0                   =         0
c    B10  5/66                =         0.07575757575
c    B11  0                   =         0
c    B12 -691/2730            =        -0.25311355311
c    B13  0                   =         0
c    B14  7/6                 =         1.16666666666
c    B15  0                   =         0
c    B16 -3617/510            =        -7.09215686274
c    B17  0                   =         0
c    B18  43867/798           =        54.97117794486
c    B19  0                   =         0
c    B20 -174611/330          =      -529.12424242424
c    B21  0                   =         0
c    B22  854513/138          =      6192.123
c    B23  0                   =         0
c    B24 -236364091/2730      =    -86580.257
c    B25  0                   =         0
c    B26  8553103/6           =   1425517.16666
c    B27  0                   =         0
c    B28 -23749461029/870     = -27298231.0678
c    B29  0                   =         0
c    B30  8615841276005/14322 = 601580873.901
c
c  Recursion:
c
c    With C(N+1,K) denoting the standard binomial coefficient,
c
c    B(0) = 1.0
c    B(N) = - ( sum ( 0 <= K < N ) C(N+1,K) * B(K) ) / C(N+1,N)
c
c  Special Values:
c
c    Except for B(1), all Bernoulli numbers of odd index are 0.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    05 February 2003
c
c  Parameters:
c
c    Input, integer N, the order of the Bernoulli number
c    to compute.
c
c    Output, double precision B, the desired Bernoulli number.
c
      implicit none

      double precision b
      integer it
      integer it_max
      parameter ( it_max = 1000 )
      integer n
      double precision pi
      parameter ( pi = 3.141592653589793D+00 )
      double precision r8_factorial
      double precision sum2
      double precision term
      double precision tol
      parameter ( tol = 5.0D-07 )

      if ( n .lt. 0 ) then

        b = 0.0D+00

      else if ( n .eq. 0 ) then

        b = 1.0D+00

      else if ( n .eq. 1 ) then

        b = -0.5D+00

      else if ( n .eq. 2 ) then

        b = 1.0D+00 / 6.0D+00

      else if ( mod ( n, 2 ) .eq. 1 ) then

        b = 0.0D+00

      else

        sum2 = 0.0D+00

        do it = 1, it_max

          term = 1.0D+00 / dble ( it**n )
          sum2 = sum2 + term

          if ( abs ( term ) .lt. tol .or. 
     &         abs ( term ) .lt. tol * abs ( sum2 ) ) then
            go to 10
          end if

        end do

10      continue

        b = 2.0D+00 * sum2 * r8_factorial ( n ) / ( 2.0D+00 * pi )**n

        if ( mod ( n, 4 ) .eq. 0 ) then
          b = - b
        end if

      end if

      return
      end
      subroutine bernoulli_number_values ( n_data, n, c )

c*********************************************************************72
c
cc BERNOULLI_NUMBER_VALUES returns some values of the Bernoulli numbers.
c
c  Discussion:
c
c    The Bernoulli numbers are rational.
c
c    If we define the sum of the M-th powers of the first N integers as:
c
c      SIGMA(M,N) = sum ( 0 <= I <= N ) I**M
c
c    and let C(I,J) be the combinatorial coefficient:
c
c      C(I,J) = Ic / ( ( I - J )c * Jc )
c
c    then the Bernoulli numbers B(J) satisfy:
c
c      SIGMA(M,N) = 1/(M+1) * sum ( 0 <= J <= M ) C(M+1,J) B(J) * (N+1)**(M+1-J)
c
c    In Mathematica, the function can be evaluated by:
c
c      BernoulliB[n]
c
c    With C(N+1,K) denoting the standard binomial coefficient,
c
c      B(0) = 1.0
c      B(N) = - ( sum ( 0 <= K .lt. N ) C(N+1,K) * B(K) ) / C(N+1,N)
c
c    Except for B(1), all Bernoulli numbers of odd index are 0.
c
c  First values:
c
c   B0  1                   =         1.00000000000
c   B1 -1/2                 =        -0.50000000000
c   B2  1/6                 =         1.66666666666
c   B3  0                   =         0
c   B4 -1/30                =        -0.03333333333
c   B5  0                   =         0
c   B6  1/42                =         0.02380952380
c   B7  0                   =         0
c   B8 -1/30                =        -0.03333333333
c   B9  0                   =         0
c  B10  5/66                =         0.07575757575
c  B11  0                   =         0
c  B12 -691/2730            =        -0.25311355311
c  B13  0                   =         0
c  B14  7/6                 =         1.16666666666
c  B15  0                   =         0
c  B16 -3617/510            =        -7.09215686274
c  B17  0                   =         0
c  B18  43867/798           =        54.97117794486
c  B19  0                   =         0
c  B20 -174611/330          =      -529.12424242424
c  B21  0                   =         0
c  B22  854,513/138         =      6192.123
c  B23  0                   =         0
c  B24 -236364091/2730      =    -86580.257
c  B25  0                   =         0
c  B26  8553103/6           =   1425517.16666
c  B27  0                   =         0
c  B28 -23749461029/870     = -27298231.0678
c  B29  0                   =         0
c  B30  8615841276005/14322 = 601580873.901
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    19 March 2007
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
c    Output, integer N, the order of the Bernoulli number.
c
c    Output, double precision C, the value of the Bernoulli number.
c
      implicit none

      integer n_max
      parameter ( n_max = 10 )

      double precision c
      double precision c_vec(n_max)
      integer n
      integer n_data
      integer n_vec(n_max)

      save c_vec
      save n_vec

      data c_vec /
     &   0.1000000000000000D+01,
     &  -0.5000000000000000D+00,
     &   0.1666666666666667D+00,
     &   0.0000000000000000D+00,
     &  -0.3333333333333333D-01,
     &  -0.2380952380952380D-01,
     &  -0.3333333333333333D-01,
     &   0.7575757575757575D-01,
     &  -0.5291242424242424D+03,
     &   0.6015808739006424D+09 /
      data n_vec /
     &   0,  1,  2,  3,  4, 6,  8, 10, 20, 30 /

      if ( n_data .lt. 0 ) then
        n_data = 0
      end if

      n_data = n_data + 1

      if ( n_max .lt. n_data ) then
        n_data = 0
        n = 0
        c = 0.0D+00
      else
        n = n_vec(n_data)
        c = c_vec(n_data)
      end if

      return
      end
      subroutine bernoulli_poly ( n, x, bx )

c*********************************************************************72
c
cc BERNOULLI_POLY evaluates the Bernoulli polynomial of order N at X.
c
c  Discussion:
c
c    B(N,0) = B(N,1) = B(N), the N-th Bernoulli number.
c
c    B'(N,X) = N * B(N-1,X)
c
c    B(N,X+1) - B(N,X) = N * X^(N-1)
c    B(N,X) = (-1)^N * B(N,1-X)
c
c    The formula is:
c
c      B(N,X) = sum ( 1 <= K <= N ) B(K) * C(N,K) * X^(N-K)
c
c  First values:
c
c    B(0,X)  1
c    B(1,X)  X    - 1/2
c    B(2,X)  X^2 -   X      +  1/6
c    B(3,X)  X^3 - 3/2*X^2 +  1/2*X
c    B(4,X)  X^4 - 2*X^3   +      X^2 - 1/30
c    B(5,X)  X^5 - 5/2*X^4 +  5/3*X^3 - 1/6*X
c    B(6,X)  X^6 - 3*X^5   +  5/2*X^4 - 1/2*X^2 + 1/42
c    B(7,X)  X^7 - 7/2*X^6 +  7/2*X^5 - 7/6*X^3 + 1/6*X
c    B(8,X)  X^8 - 4*X^7   + 14/3*X^6 - 7/3*X^4 + 2/3*X^2 - 1/30
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    05 December 2007
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the order of the Bernoulli polynomial to
c    be evaluated.  N must be 0 or greater.
c
c    Input, double precision X, the value of X at which the polynomial is to
c    be evaluated.
c
c    Output, double precision BX, the value of B(N,X).
c
      implicit none

      integer n

      double precision bx
      integer i
      integer ido
      integer iwork(0:n)
      double precision work(0:n)
      double precision x

      call bernoulli_number ( n, work )

      ido = 0
      call comb_row ( ido, n, iwork )

      bx = 1.0D+00
      do i = 1, n
        bx = bx * x + work(i) * dble ( iwork(i) )
      end do

      return
      end
      subroutine bernoulli_poly2 ( n, x, bx )

c*********************************************************************72
c
cc BERNOULLI_POLY2 evaluates the N-th Bernoulli polynomial at X.
c
c  Discussion:
c
c    BERN(N,0) = BERN(N,1) = B(N), the N-th Bernoulli number.
c
c    B'(N,X) = N*B(N-1,X).
c
c    B(N,X+1) - B(N,X) = N*X^(N-1)
c    B(N,X) = (-1)**N * B(N,1-X)
c
c    The formula is:
c
c      B(N,X) = sum ( 1 <= K <= N ) B(K)*C(N,K)*X^(N-K)
c
c  First values:
c
c    B(0,X)  1
c    B(1,X)  X    - 1/2
c    B(2,X)  X^2 -   X      +  1/6
c    B(3,X)  X^3 - 3*X^2/2 +    X/2
c    B(4,X)  X^4 - 2*X^3   +    X^2   - 1/30
c    B(5,X)  X^5 - 5*X^4/2 +  5*X^3/3 -   X/6
c    B(6,X)  X^6 - 3*X^5   +  5*X^4/2 -   X^2/2 + 1/42
c    B(7,X)  X^7 - 7*X^6/2 +  7*X^5/2 - 7*X^3/6 +   X/6
c    B(8,X)  X^8 - 4*X^7   + 14*X^6/3 - 7*X^4/3 + 2*X^2/3 - 1/30
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    06 December 2007
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the order of the Bernoulli polynomial to
c    be evaluated.  N must be 0 or greater.
c
c    Input, double precision X, the value at which the polynomial is to
c    be evaluated.
c
c    Output, double precision BX, the value of B(N,X).
c
      implicit none

      double precision b
      double precision bx
      double precision fact
      integer i
      integer n
      double precision x

      fact = 1.0D+00

      call bernoulli_number3 ( 0, b )

      bx = b

      do i = 1, n
        fact = fact * dble ( n + 1 - i ) / dble ( i )
        call bernoulli_number3 ( i, b )
        bx = bx * x + fact * b
      end do

      return
      end
      subroutine bernstein_poly ( n, x, bern )

c*********************************************************************72
c
cc BERNSTEIN_POLY evaluates the Bernstein polynomials at a point X.
c
c  Discussion:
c
c    The Bernstein polynomials are assumed to be based on [0,1].
c
c    The formula is:
c
c      B(N,I,X) = [N!/(I!*(N-I)!)] * (1-X)**(N-I) * X^I
c
c    B(N,I,X) has a unique maximum value at X = I/N.
c
c    B(N,I,X) has an I-fold zero at 0 and and N-I fold zero at 1.
c
c    B(N,I,1/2) = C(N,K) / 2**N
c
c    For a fixed X and N, the polynomials add up to 1:
c
c      Sum ( 0 <= I <= N ) B(N,I,X) = 1
c
c  First values:
c
c    B(0,0,X) = 1
c
c    B(1,0,X) =      1-X
c    B(1,1,X) =                X
c
c    B(2,0,X) =     (1-X)^2
c    B(2,1,X) = 2 * (1-X)    * X
c    B(2,2,X) =                X^2
c
c    B(3,0,X) =     (1-X)**3
c    B(3,1,X) = 3 * (1-X)^2 * X
c    B(3,2,X) = 3 * (1-X)    * X^2
c    B(3,3,X) =                X^3
c
c    B(4,0,X) =     (1-X)**4
c    B(4,1,X) = 4 * (1-X)**3 * X
c    B(4,2,X) = 6 * (1-X)^2 * X^2
c    B(4,3,X) = 4 * (1-X)    * X^3
c    B(4,4,X) =                X^4
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    03 July 2008
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the degree of the Bernstein polynomials 
c    to be used.  For any N, there is a set of N+1 Bernstein polynomials,
c    each of degree N, which form a basis for polynomials on [0,1].
c
c    Input, double precision X, the evaluation point.
c
c    Output, double precision BERN(0:N), the values of the N+1 
c    Bernstein polynomials at X.
c
      implicit none

      integer  n

      double precision bern(0:n)
      integer i
      integer j
      double precision x

      if ( n .eq. 0 ) then
 
        bern(0) = 1.0D+00
 
      else if ( 0 .lt. n ) then
 
        bern(0) = 1.0D+00 - x
        bern(1) = x
 
        do i = 2, n
          bern(i) = x * bern(i-1)
          do j = i - 1, 1, -1
            bern(j) =             x   * bern(j-1) 
     &              + ( 1.0D+00 - x ) * bern(j)
          end do
          bern(0) = ( 1.0D+00 - x ) * bern(0)
        end do
 
      end if
 
      return
      end
      subroutine bernstein_poly_values ( n_data, n, k, x, b )

c*********************************************************************72
c
cc BERNSTEIN_POLY_VALUES returns some values of the Bernstein polynomials.
c
c  Discussion:
c
c    The Bernstein polynomials are assumed to be based on [0,1].
c
c    The formula for the Bernstein polynomials is
c
c      B(N,I,X) = [Nc/(Ic*(N-I)c)] * (1-X)**(N-I) * X^I
c
c    In Mathematica, the function can be evaluated by:
c
c      Binomial[n,i] * (1-x)^(n-i) * x^i
c
c    B(N,I,X) has a unique maximum value at X = I/N.
c
c    B(N,I,X) has an I-fold zero at 0 and and N-I fold zero at 1.
c
c    B(N,I,1/2) = C(N,K) / 2**N
c
c    For a fixed X and N, the polynomials add up to 1:
c
c      Sum ( 0 <= I <= N ) B(N,I,X) = 1
c
c  First values:
c
c    B(0,0,X) = 1
c
c    B(1,0,X) =      1-X
c    B(1,1,X) =                X
c
c    B(2,0,X) =     (1-X)^2
c    B(2,1,X) = 2 * (1-X)    * X
c    B(2,2,X) =                X^2
c
c    B(3,0,X) =     (1-X)**3
c    B(3,1,X) = 3 * (1-X)^2 * X
c    B(3,2,X) = 3 * (1-X)    * X^2
c    B(3,3,X) =                X^3
c
c    B(4,0,X) =     (1-X)**4
c    B(4,1,X) = 4 * (1-X)**3 * X
c    B(4,2,X) = 6 * (1-X)^2 * X^2
c    B(4,3,X) = 4 * (1-X)    * X^3
c    B(4,4,X) =                X^4
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    20 March 2007
c
c  Author:
c
c    John Burkardt
c
c  Reference:
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
c    Output, integer N, the degree of the polynomial.
c
c    Output, integer K, the index of the polynomial.
c
c    Output, double precision X, the argument of the polynomial.
c
c    Output, double precision B, the value of the polynomial B(N,K,X).
c
      implicit none

      integer n_max
      parameter ( n_max = 15 )

      double precision b
      double precision b_vec(n_max)
      integer k
      integer k_vec(n_max)
      integer n
      integer n_data
      integer n_vec(n_max)
      double precision x
      double precision x_vec(n_max)

      save b_vec
      save k_vec
      save n_vec
      save x_vec

      data b_vec /
     &  0.1000000000000000D+01,
     &  0.7500000000000000D+00,
     &  0.2500000000000000D+00,
     &  0.5625000000000000D+00,
     &  0.3750000000000000D+00,
     &  0.6250000000000000D-01,
     &  0.4218750000000000D+00,
     &  0.4218750000000000D+00,
     &  0.1406250000000000D+00,
     &  0.1562500000000000D-01,
     &  0.3164062500000000D+00,
     &  0.4218750000000000D+00,
     &  0.2109375000000000D+00,
     &  0.4687500000000000D-01,
     &  0.3906250000000000D-02 /
      data k_vec /
     &  0,
     &  0, 1,
     &  0, 1, 2,
     &  0, 1, 2, 3,
     &  0, 1, 2, 3, 4 /
      data n_vec /
     &  0,
     &  1, 1,
     &  2, 2, 2,
     &  3, 3, 3, 3,
     &  4, 4, 4, 4, 4 /
      data x_vec /
     &  0.25D+00,
     &  0.25D+00,
     &  0.25D+00,
     &  0.25D+00,
     &  0.25D+00,
     &  0.25D+00,
     &  0.25D+00,
     &  0.25D+00,
     &  0.25D+00,
     &  0.25D+00,
     &  0.25D+00,
     &  0.25D+00,
     &  0.25D+00,
     &  0.25D+00,
     &  0.25D+00 /

      if ( n_data .lt. 0 ) then
        n_data = 0
      end if

      n_data = n_data + 1

      if ( n_max .lt. n_data ) then
        n_data = 0
        n = 0
        k = 0
        x = 0.0D+00
        b = 0.0D+00
      else
        n = n_vec(n_data)
        k = k_vec(n_data)
        x = x_vec(n_data)
        b = b_vec(n_data)
      end if

      return
      end
      function beta ( x, y )

c*********************************************************************72
c
cc BETA returns the value of the Beta function.
c
c  Discussion:
c
c    The Beta function can be defined in terms of the Gamma function:
c
c      BETA(X,Y) = ( GAMMA(X) * GAMMA(Y) ) / GAMMA(X+Y)
c
c      Both X and Y must be greater than 0.
c
c    The function has the following properties:
c
c      BETA(X,Y) = BETA(Y,X).
c      BETA(X,Y) = Integral ( 0 <= T <= 1 ) T**(X-1) (1-T)**(Y-1) dT.
c      BETA(X,Y) = GAMMA(X) * GAMMA(Y) / GAMMA(X+Y)
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
c    John Burkardt
c
c  Parameters:
c
c    Input, double precision X, Y, the two parameters that define 
c    the Beta function.  X and Y must be greater than 0.
c
c    Output, double precision BETA, the value of the Beta function.
c
      implicit none

      double precision beta
      double precision r8_gamma_log
      double precision x
      double precision y

      if ( x .le. 0.0D+00 .or. y .le. 0.0D+00 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'BETA - Fatal error!'
        write ( *, '(a)' ) '  Both X and Y must be greater than 0.'
        stop
      end if

      beta = exp ( r8_gamma_log ( x ) + r8_gamma_log ( y ) 
     &  - r8_gamma_log ( x + y ) )

      return
      end
      subroutine beta_values ( n_data, x, y, fxy )

c*********************************************************************72
c
cc BETA_VALUES returns some values of the Beta function.
c
c  Discussion:
c
c    Beta(X,Y) = ( Gamma(X) * Gamma(Y) ) / Gamma(X+Y)
c
c    Both X and Y must be greater than 0.
c
c    In Mathematica, the function can be evaluated by:
c
c      Beta[X,Y]
c
c    Beta(X,Y) = Beta(Y,X).
c    Beta(X,Y) = Integral ( 0 .lt.= T .lt.= 1 ) T**(X-1) (1-T)**(Y-1) dT.
c    Beta(X,Y) = Gamma(X) * Gamma(Y) / Gamma(X+Y)
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    20 March 2007
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
c    Output, double precision X, Y, the arguments of the function.
c
c    Output, double precision FXY, the value of the function.
c
      implicit none

      integer n_max
      parameter ( n_max = 17 )

      double precision b_vec(n_max)
      double precision fxy
      integer n_data
      double precision x
      double precision x_vec(n_max)
      double precision y
      double precision y_vec(n_max)

      save b_vec
      save x_vec
      save y_vec

      data b_vec /
     &  0.5000000000000000D+01,
     7  0.2500000000000000D+01,
     &  0.1666666666666667D+01,
     &  0.1250000000000000D+01,
     &  0.5000000000000000D+01,
     &  0.2500000000000000D+01,
     &  0.1000000000000000D+01,
     &  0.1666666666666667D+00,
     &  0.3333333333333333D-01,
     &  0.7142857142857143D-02,
     &  0.1587301587301587D-02,
     &  0.2380952380952381D-01,
     &  0.5952380952380952D-02,
     &  0.1984126984126984D-02,
     &  0.7936507936507937D-03,
     &  0.3607503607503608D-03,
     &  0.8325008325008325D-04 /
      data x_vec /
     &  0.2D+00,
     &  0.4D+00,
     &  0.6D+00,
     &  0.8D+00,
     &  1.0D+00,
     &  1.0D+00,
     &  1.0D+00,
     &  2.0D+00,
     &  3.0D+00,
     &  4.0D+00,
     &  5.0D+00,
     &  6.0D+00,
     &  6.0D+00,
     &  6.0D+00,
     &  6.0D+00,
     &  6.0D+00,
     &  7.0D+00 /
      data y_vec /
     &  1.0D+00,
     &  1.0D+00,
     &  1.0D+00,
     &  1.0D+00,
     &  0.2D+00,
     &  0.4D+00,
     &  1.0D+00,
     &  2.0D+00,
     &  3.0D+00,
     &  4.0D+00,
     &  5.0D+00,
     &  2.0D+00,
     &  3.0D+00,
     &  4.0D+00,
     &  5.0D+00,
     &  6.0D+00,
     &  7.0D+00 /

      if ( n_data .lt. 0 ) then
        n_data = 0
      end if

      n_data = n_data + 1

      if ( n_max .lt. n_data ) then
        n_data = 0
        x = 0.0D+00
        y = 0.0D+00
        fxy = 0.0D+00
      else
        x = x_vec(n_data)
        y = y_vec(n_data)
        fxy = b_vec(n_data)
      end if

      return
      end
      subroutine bpab ( n, x, a, b, bern )

c*********************************************************************72
c
cc BPAB evaluates at X the Bernstein polynomials based in [A,B].
c
c  Discussion:
c
c    The formula is:
c
c      BERN(N,I,X) = [N!/(I!*(N-I)!)] * (B-X)**(N-I) * (X-A)**I / (B-A)**N
c
c  First values:
c
c    B(0,0,X) =   1
c
c    B(1,0,X) = (      B-X                ) / (B-A)
c    B(1,1,X) = (                 X-A     ) / (B-A)
c
c    B(2,0,X) = (     (B-X)^2            ) / (B-A)^2
c    B(2,1,X) = ( 2 * (B-X)    * (X-A)    ) / (B-A)^2
c    B(2,2,X) = (                (X-A)^2 ) / (B-A)^2
c
c    B(3,0,X) = (     (B-X)**3            ) / (B-A)**3
c    B(3,1,X) = ( 3 * (B-X)^2 * (X-A)    ) / (B-A)**3
c    B(3,2,X) = ( 3 * (B-X)    * (X-A)^2 ) / (B-A)**3
c    B(3,3,X) = (                (X-A)**3 ) / (B-A)**3
c
c    B(4,0,X) = (     (B-X)**4            ) / (B-A)**4
c    B(4,1,X) = ( 4 * (B-X)**3 * (X-A)    ) / (B-A)**4
c    B(4,2,X) = ( 6 * (B-X)^2 * (X-A)^2 ) / (B-A)**4
c    B(4,3,X) = ( 4 * (B-X)    * (X-A)**3 ) / (B-A)**4
c    B(4,4,X) = (                (X-A)**4 ) / (B-A)**4
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    03 July 2008
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the degree of the Bernstein polynomials 
c    to be used.  For any N, there is a set of N+1 Bernstein polynomials, 
c    each of degree N, which form a basis for polynomials on [A,B].
c
c    Input, double precision X, the point at which the polynomials 
c    are to be evaluated.
c
c    Input, double precision A, B, the endpoints of the interval on which the
c    polynomials are to be based.  A and B should not be equal.
c
c    Output, double precision BERN(0:N), the values of the N+1
c    Bernstein polynomials at X.
c
      implicit none

      integer n

      double precision a
      double precision b
      double precision bern(0:n)
      integer i
      integer j
      double precision x

      if ( b .eq. a ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'BPAB - Fatal error!'
        write ( *, '(a,g14.6)' ) '  A = B = ', a
        stop
      end if

      if ( n .eq. 0 ) then
 
        bern(0) = 1.0D+00
 
      else if ( 0 .lt. n ) then
 
        bern(0) = ( b - x ) / ( b - a )
        bern(1) = ( x - a ) / ( b - a )
 
        do i = 2, n
          bern(i) = ( x - a ) * bern(i-1) / ( b - a )
          do j = i - 1, 1, -1
            bern(j) = ( ( b - x     ) * bern(j)     
     &                + (     x - a ) * bern(j-1) ) 
     &                / ( b     - a )
          end do
          bern(0) = ( b - x ) * bern(0) / ( b - a )
        end do
 
      end if
 
      return
      end
      function c4_acos ( z )

c*********************************************************************72
c
cc C4_ACOS evaluates the inverse cosine of a C4 argument.
c
c  Discussion:
c
c    FORTRAN77 does not have an intrinsic inverse cosine function for C4 arguments.
c
c    Here we use the relationship:
c
c      C4_ACOS ( Z ) = pi/2 - C4_ASIN ( Z ).
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    29 November 2007
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, complex Z, the argument.
c
c    Output, complex C4_ACOS, the function value.
c
      implicit none

      complex c4_acos
      complex c4_asin
      complex z

      c4_acos = 1.57079632679489661923E+00 - c4_asin ( z )

      return
      end
      function c4_acosh ( z )

c*********************************************************************72
c
cc C4_ACOSH evaluates the inverse hyperbolic cosine of a C4 argument.
c
c  Discussion:
c
c    FORTRAN77 does not have an intrinsic inverse hyperbolic 
c    cosine function for C4 arguments.
c
c    Here we use the relationship:
c
c      C4_ACOSH ( Z ) = i * C4_ACOS ( Z ).
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    29 November 2007
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, complex Z, the argument.
c
c    Output, complex C4_ACOSH, the function value.
c
      implicit none

      complex c4_acos
      complex c4_acosh
      complex i
      complex z

      i = cmplx ( 0.0E+00, 1.0E+00 )

      c4_acosh = i * c4_acos ( z )

      return
      end
      function c4_argument ( x )

c*********************************************************************72
c
cc C4_ARGUMENT returns the argument of a C4.
c
c  Discussion:
c
c    A C4 is a complex value.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    30 November 2008
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, complex X, the value whose argument is desired.
c
c    Output, real C4_ARGUMENT, the argument of X.
c
      implicit none

      real    c4_argument
      complex x
      real    xi
      real    xr

      xr = real ( x )
      xi = aimag ( x )

      if ( xi .eq. 0.0E+00 .and. xr .eq. 0.0E+00 ) then

        c4_argument = 0.0E+00

      else

        c4_argument = atan2 ( xi, xr )

      end if

      return
      end
      function c4_asin ( z )

c*********************************************************************72
c
cc C4_ASIN evaluates the inverse sine of a C4 argument.
c
c  Discussion:
c
c    FORTRAN77 does not have an intrinsic inverse sine function for C4 arguments.
c
c    Here we use the relationship:
c
c      C4_ASIN ( Z ) = - i * log ( i * z + sqrt ( 1 - z * z ) )
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    29 November 2007
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, complex Z, the argument.
c
c    Output, complex C4_ASIN, the function value.
c
      implicit none

      complex c4_asin
      complex i
      complex z

      i = cmplx ( 0.0E+00, 1.0E+00 )

      c4_asin = - i * clog ( i * z + csqrt ( 1.0E+00 - z * z ) )

      return
      end
      function c4_asinh ( z )

c*********************************************************************72
c
cc C4_ASINH evaluates the inverse hyperbolic sine of a C4 argument.
c
c  Discussion:
c
c    FORTRAN77 does not have an intrinsic inverse hyperbolic 
c    sine function for C4 arguments.
c
c    Here we use the relationship:
c
c      C4_ASINH ( Z ) = - i * C4_ASIN ( i * Z ).
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    29 November 2007
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, complex Z, the argument.
c
c    Output, complex C4_ASINH, the function value.
c
      implicit none

      complex c4_asin
      complex c4_asinh
      complex i
      complex z

      i = cmplx ( 0.0E+00, 1.0E+00 )

      c4_asinh = - i * c4_asin ( i * z )

      return
      end
      function c4_atan ( z )

c*********************************************************************72
c
cc C4_ATAN evaluates the inverse tangent of a C4 argument.
c
c  Discussion:
c
c    FORTRAN77 does not have an intrinsic inverse tangent function
c    for C4 arguments.
c
c    Here we use the relationship:
c
c      C4_ATAN ( Z ) = ( i / 2 ) * log ( ( 1 - i * z ) / ( 1 + i * z ) )
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    29 November 2007
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, complex Z, the argument.
c
c    Output, complex C4_ATAN, the function value.
c
      implicit none

      complex arg
      complex c4_atan
      complex i
      complex z

      i = cmplx ( 0.0E+00, 1.0E+00 )

      arg = ( 1.0E+00 - i * z ) / ( 1.0E+00 + i * z )

      c4_atan = 0.5E+00 * i * clog ( arg )

      return
      end
      function c4_atanh ( z )

c*********************************************************************72
c
cc C4_ATANH evaluates the inverse hyperbolic tangent of a C4 argument.
c
c  Discussion:
c
c    FORTRAN77 does not have an intrinsic inverse hyperbolic 
c    tangent function for C4 arguments.
c
c    Here we use the relationship:
c
c      C4_ATANH ( Z ) = - i * C4_ATAN ( i * Z ).
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    29 November 2007
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, complex Z, the argument.
c
c    Output, complex C4_ATANH, the function value.
c
      implicit none

      complex c4_atan
      complex c4_atanh
      complex i
      complex z

      i = cmplx ( 0.0E+00, 1.0E+00 )

      c4_atanh = - i * c4_atan ( i * z )

      return
      end
      function c4_cosh ( z )

c*********************************************************************72
c
cc C4_COSH evaluates the hyperbolic cosine of a C4 argument.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    30 November 2007
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, complex Z, the argument.
c
c    Output, complex C4_COSH, the function value.
c
      implicit none

      complex c4_cosh
      complex z

      c4_cosh =  ( cexp ( z ) + cexp ( - z ) ) / 2.0E+00

      return
      end
      function c4_log ( z )

c*********************************************************************72
c
cc C4_LOG evaluates the logarithm of a C4.
c
c  Discussion:
c
c    FORTRAN77 has a logarithm function for C4 arguments, "clog ( z )".
c
c    Here we use the relationship:
c
c      C4_LOG ( Z ) = LOG ( R ) + i * ARG ( R )
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    30 November 2008
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, complex Z, the argument.
c
c    Output, real C4_LOG, the function value.
c
      implicit none

      real    arg
      real    c4_argument
      complex c4_log
      real    c4_magnitude
      complex i
      real    r
      complex z

      i = cmplx ( 0.0E+00, 1.0E+00 )

      arg = c4_argument ( arg )
      r = c4_magnitude ( arg )

      c4_log = alog ( r ) + i * arg
 
      return
      end
      function c4_magnitude ( x )

c*****************************************************************************80
c
cc C4_MAGNITUDE returns the magnitude of a C4.
c
c  Discussion:
c
c    A C4 is a complex value.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    30 November 2008
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, complex X, the value whose magnitude is desired.
c
c    Output, real C4_MAGNITUDE, the magnitude of X.
c
      implicit none

      real    c4_magnitude
      complex x
      real    xi
      real    xr

      xr = real ( x )
      xi = aimag ( x )

      c4_magnitude = sqrt ( xr * xr + xi * xi )

      return
      end
      function c4_sinh ( z )

c*********************************************************************72
c
cc C4_SINH evaluates the hyperbolic sine of a C4 argument.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    30 November 2007
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, complex Z, the argument.
c
c    Output, complex C4_SINH, the function value.
c
      implicit none

      complex c4_sinh
      complex z

      c4_sinh =  ( cexp ( z ) - cexp ( - z ) ) / 2.0E+00

      return
      end
      function c4_sqrt ( x )

c*********************************************************************72
c
cc C4_SQRT returns the principal square root of a C4.
c
c  Discussion:
c
c    A C4 is a complex value.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    30 November 2008
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, complex X, the number whose square root is desired.
c
c    Output, complex C4_SQRT, the square root of X.
c
      implicit none

      real    argument
      real    c4_argument
      real    c4_magnitude
      complex c4_sqrt
      real    magnitude
      complex x

      argument = c4_argument ( x )
      magnitude = c4_magnitude ( x )

      if ( magnitude .eq. 0.0E+00 ) then

        c4_sqrt = cmplx ( 0.0E+00, 0.0E+00 )

      else

        c4_sqrt = sqrt ( magnitude ) 
     &    * cmplx ( cos ( argument / 2.0E+00 ), 
     &              sin ( argument / 2.0E+00 ) )

      end if

      return
      end
      function c4_tanh ( z )

c*********************************************************************72
c
cc C4_TANH evaluates the complex hyperbolic tangent.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    30 November 2007
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, complex Z, the argument.
c
c    Output, complex C4_TANH, the function value.
c
      implicit none

      complex c4_tanh
      complex z

      c4_tanh =  ( cexp ( z ) - cexp ( - z ) ) 
     &         / ( cexp ( z ) + cexp ( - z ) )

      return
      end
      function c4_uniform_01 ( seed )

c*********************************************************************72
c
cc C4_UNIFORM_01 returns a unit pseudorandom C4.
c
c  Discussion:
c
c    A C4 is a complex single precision value.
c
c    The angle should be uniformly distributed between 0 and 2 * PI,
c    the square root of the radius uniformly distributed between 0 and 1.
c
c    This results in a uniform distribution of values in the unit circle.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    15 March 2005
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
c    Input/output, integer SEED, the "seed" value, which should NOT be 0.
c    On output, SEED has been updated.
c
c    Output, complex C4_UNIFORM_01, a pseudorandom complex value.
c
      implicit none

      complex c4_uniform_01
      integer i4_huge
      integer k
      real pi
      parameter ( pi = 3.1415926E+00 )
      real r
      integer seed
      real theta

      if ( seed .eq. 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'C4_UNIFORM_01 - Fatal error!'
        write ( *, '(a)' ) '  Input value of SEED = 0.'
        stop
      end if

      k = seed / 127773

      seed = 16807 * ( seed - k * 127773 ) - k * 2836

      if ( seed .lt. 0 ) then
        seed = seed + i4_huge ( )
      end if

      r = sqrt ( real ( dble ( seed ) * 4.656612875D-10 ) )

      k = seed / 127773

      seed = 16807 * ( seed - k * 127773 ) - k * 2836

      if ( seed .lt. 0 ) then
        seed = seed + i4_huge ( )
      end if

      theta = 2.0E+00 * pi
     &  * real ( dble ( seed ) * 4.656612875D-10 )

      c4_uniform_01 = r * cmplx ( cos ( theta ), sin ( theta ) )

      return
      end
      function c8_argument ( x )

c*********************************************************************72
c
cc C8_ARGUMENT returns the argument of a C8.
c
c  Discussion:
c
c    A C8 is a double complex value.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    30 November 2008
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, double complex X, the value whose argument is desired.
c
c    Output, double precision C8_ARGUMENT, the argument of X.
c
      implicit none

      double precision c8_argument
      double complex   x
      double precision xi
      double precision xr

      xr = dreal ( x )
      xi = dimag ( x )

      if ( xi .eq. 0.0D+00 .and. xr .eq. 0.0D+00 ) then

        c8_argument = 0.0D+00

      else

        c8_argument = datan2 ( xi, xr )

      end if

      return
      end
      function c8_atan ( z )

c*********************************************************************72
c
cc C8_ATAN evaluates the inverse tangent of a C8.
c
c  Discussion:
c
c    FORTRAN77 does not have an intrinsic inverse tangent function
c    for C8 arguments.
c
c    FORTRAN77 does not have a logarithm function for C8 arguments!
c
c    Here we use the relationship:
c
c      C8_ATAN ( Z ) = ( i / 2 ) * log ( ( 1 - i * z ) / ( 1 + i * z ) )
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    30 November 2008
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, double complex Z, the argument.
c
c    Output, double complex C8_ATAN, the function value.
c
      implicit none

      double complex arg
      double complex c8_atan
      double complex c8_log
      double complex i
      double complex z

      i = dcmplx ( 0.0D+00, 1.0D+00 )

      arg = ( 1.0D+00 - i * z ) / ( 1.0D+00 + i * z )

      c8_atan = 0.5D+00 * i * c8_log ( arg ) 

      return
      end
      function c8_log ( z )

c*********************************************************************72
c
cc C8_LOG evaluates the logarithm of a C8.
c
c  Discussion:
c
c    FORTRAN77 does not have a logarithm function for C8 arguments!
c
c    Here we use the relationship:
c
c      C8_LOG ( Z ) = LOG ( R ) + i * ARG ( R )
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    30 November 2008
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, double complex Z, the argument.
c
c    Output, double complex C8_LOG, the function value.
c
      implicit none

      double precision arg
      double precision c8_argument
      double complex   c8_log
      double precision c8_magnitude
      double complex   i
      double precision r
      double complex   z

      i = dcmplx ( 0.0D+00, 1.0D+00 )

      arg = c8_argument ( arg )
      r = c8_magnitude ( arg )

      c8_log = dlog ( r ) + i * arg
 
      return
      end
      function c8_magnitude ( x )

c*****************************************************************************80
c
cc C8_MAGNITUDE returns the magnitude of a C8.
c
c  Discussion:
c
c    A C8 is a double complex value.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    30 November 2008
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, double complex X, the value whose magnitude is desired.
c
c    Output, double precision C8_MAGNITUDE, the magnitude of X.
c
      implicit none

      double precision c8_magnitude
      double complex   x
      double precision xi
      double precision xr

      xr = dreal ( x )
      xi = dimag ( x )

      c8_magnitude = sqrt ( xr * xr + xi * xi )

      return
      end
      function c8_sqrt ( x )

c*********************************************************************72
c
cc C8_SQRT returns the principal square root of a C8.
c
c  Discussion:
c
c    A C8 is a double complex value.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    30 November 2008
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, double complex X, the number whose square root is desired.
c
c    Output, double complex C8_SQRT, the square root of X.
c
      implicit none

      double precision argument
      double precision c8_argument
      double precision c8_magnitude
      double complex   c8_sqrt
      double precision magnitude
      double complex   x

      argument = c8_argument ( x )
      magnitude = c8_magnitude ( x )

      if ( magnitude .eq. 0.0D+00 ) then

        c8_sqrt = dcmplx ( 0.0D+00, 0.0D+00 )

      else

        c8_sqrt = dsqrt ( magnitude ) 
     &    * dcmplx ( dcos ( argument / 2.0D+00 ), 
     &               dsin ( argument / 2.0D+00 ) )

      end if

      return
      end
      subroutine cardan ( n, x, s, cx )

c*********************************************************************72
c
cc CARDAN evaluates the Cardan polynomials.
c
c  Discussion:
c
c    Writing the N-th polynomial in terms of its coefficients:
c
c      C(N,S,X) = sum ( 0 <= I <= N ) D(N,I) * S**(N-I)/2 * X^I
c
c    then
c
c    D(0,0) = 1
c
c    D(1,1) = 1
c    D(1,0) = 0
c
c    D(N,N) = 1
c    D(N,K) = D(N-1,K-1) - D(N-2,K)
c
c  First terms:
c
c     N  C(N,S,X)
c
c     0  2
c     1  X
c     2  X^2  -  2 S
c     3  X^3  -  3 S X
c     4  X^4  -  4 S X^2 +  2 S^2
c     5  X^5  -  5 S X^3 +  5 S^2 X
c     6  X^6  -  6 S X^4 +  9 S^2 X^2 -  2 S**3
c     7  X^7  -  7 S X^5 + 14 S^2 X^3 -  7 S**3 X
c     8  X^8  -  8 S X^6 + 20 S^2 X^4 - 16 S**3 X^2 +  2 S**4
c     9  X^9  -  9 S X^7 + 27 S^2 X^5 - 30 S**3 X^3 +  9 S**4 X
c    10  X^10 - 10 S X^8 + 35 S^2 X^6 - 50 S**3 X^4 + 25 S**4 X^2 -  2 S**5
c    11  X^11 - 11 S X^9 + 44 S^2 X^7 - 77 S**3 X^5 + 55 S**4 X^3 - 11 S**5 X
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    03 July 2008
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Thomas Osler,
c    Cardan Polynomials and the Reduction of Radicals,
c    Mathematics Magazine, 
c    Volume 74, Number 1, February 2001, pages 26-32.
c
c  Parameters:
c
c    Input, integer N, the highest polynomial to compute.
c
c    Input, double precision X, the point at which the polynomials 
c    are to be computed.
c
c    Input, double precision S, the value of the parameter, which 
c    must be positive.
c
c    Output, double precision CX(0:N), the values of the Cardan 
c    polynomials at X.
c
      implicit none

      integer n

      double precision cx(0:n)
      double precision fact
      integer i
      double precision s
      double precision s2
      double precision x
      double precision x2

      s2 = sqrt ( s )
      x2 = 0.5D+00 * x / s2

      call cheby_t_poly ( 1, n, x2, cx )

      fact = 1.0D+00

      do i = 0, n
        cx(i) = 2.0D+00 * fact * cx(i)
        fact = fact * s2
      end do
     
      return
      end
      subroutine cardan_poly_coef ( n, s, c )

c*********************************************************************72
c
cc CARDAN_POLY_COEF computes the coefficients of the N-th Cardan polynomial.
c
c  First terms:
c
c    2
c    0       1
c   -2 S     0       1
c    0      -3 S     0       1
c    2 S^2  0      -4 S     0       1
c    0       5 S^2  0      -5 S     0       1
c   -2 S**3  0       9 S^2  0      -6 S     0       1
c    0       7 S**3  0      14 S^2  0      -7 S     0       1
c    2 S**4  0     -16 S**3  0      20 S^2  0      -8 S     0        1
c    0       9 S**4  0     -30 S**3  0      27 S^2  0      -9 S      0     1
c   -2 S**5  0      25 S**4  0     -50 S**3  0      35 S^2  0      -10 S   0   1
c    0     -11 S**5  0      55 S**4  0     -77 S**3  0     +44 S^2   0   -11 S 0 1
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    03 July 2008
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Thomas Osler,
c    Cardan Polynomials and the Reduction of Radicals,
c    Mathematics Magazine, 
c    Volume 74, Number 1, February 2001, pages 26-32.
c
c  Parameters:
c
c    Input, integer N, the order of the polynomial
c
c    Input, double precision S, the value of the parameter, which 
c    must be positive.
c
c    Output, double precision C(0:N), the coefficients.  C(0) is the 
c    constant term, and C(N) is the coefficient of X^N.
c
      implicit none

      integer n

      double precision c(0:n)
      double precision cm1(0:n)
      double precision cm2(0:n)
      integer i
      integer j
      double precision s

      if ( n .lt. 0 ) then
        return
      end if

      c(0) = 2.0D+00
      do i = 1, n
        c(i) = 0.0D+00
      end do

      if ( n .eq. 0 ) then
        return
      end if

      do i = 0, n
        cm1(i) = c(i)
      end do

      c(0) = 0.0D+00
      c(1) = 1.0D+00
      do i = 2, n
        c(i) = 0.0D+00
      end do

      do i = 2, n

        do j = 0, i - 2
          cm2(j) = cm1(j)
        end do

        do j = 0, i - 1
          cm1(j) = c(j)
        end do

        c(0) = 0.0D+00
        do j = 1, i
          c(j) = cm1(j-1)
        end do

        do j = 0, i - 2
          c(j) = c(j) - s * cm2(j)
        end do

      end do

      return
      end
      subroutine catalan ( n, c )

c*********************************************************************72
c
cc CATALAN computes the Catalan numbers, from C(0) to C(N).
c
c  Discussion:
c
c    The Catalan number C(N) counts:
c
c    1) the number of binary trees on N vertices;
c    2) the number of ordered trees on N+1 vertices;
c    3) the number of full binary trees on 2N+1 vertices;
c    4) the number of well formed sequences of 2N parentheses;
c    5) the number of ways 2N ballots can be counted, in order,
c       with N positive and N negative, so that the running sum
c       is never negative;
c    6) the number of standard tableaus in a 2 by N rectangular Ferrers diagram;
c    7) the number of monotone functions from [1..N} to [1..N} which 
c       satisfy f(i) <= i for all i;
c    8) the number of ways to triangulate a polygon with N+2 vertices.
c
c    The formula is:
c
c      C(N) = (2*N)c / ( (N+1) * (Nc) * (Nc) ) 
c           = 1 / (N+1) * COMB ( 2N, N )
c           = 1 / (2N+1) * COMB ( 2N+1, N+1).
c
c    C(N) = 2 * (2*N-1) * C(N-1) / (N+1)
c    C(N) = sum ( 1 <= I <= N-1 ) C(I) * C(N-I)
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
c    03 July 2008
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
c    Output, integer C(0:N), the Catalan numbers from C(0) to C(N).
c
      implicit none

      integer n

      integer c(0:n)
      integer i

      if ( n .lt. 0 ) then
        return
      end if

      c(0) = 1
c
c  The extra parentheses ensure that the integer division is
c  done AFTER the integer multiplication.
c
      do i = 1, n
        c(i) = ( c(i-1) * 2 * ( 2 * i - 1 ) ) / ( i + 1 )
      end do
     
      return
      end
      function catalan_constant ( )

c*********************************************************************72
c
cc CATALAN_CONSTANT returns the value of Catalan's constant.
c
c  Discussion:
c
c    Catalan's constant, which may be denoted by G, is defined as
c
c      G = sum ( 0 <= K ) ( -1 )**K / ( 2 * K + 1 )^2
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    03 July 2008
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Eric Weisstein,
c    CRC Concise Encyclopedia of Mathematics,
c    CRC Press, 2002,
c    Second edition,
c    ISBN: 1584883472,
c    LC: QA5.W45
c
c  Parameters:
c
c    Output, double precision CATALAN_CONSTANT, the value of Catalan's
c    constant.
c
      implicit none

      double precision catalan_constant

      catalan_constant = 0.915965594177D+00

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
c    C(I,J) = 0 for I .lt. J
c    C(I,J) = C(I,J-1) + C(I-1,J)
c    C(I,I) is the I-th Catalan number.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    03 July 2008
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
c    Input, integer N, the index of the row of the triangle 
c    desired.  
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

          do j = 1, i - 1
            irow(j) = irow(j) + irow(j-1)
          end do

          irow(i) = irow(i-1)

        end do
     
      else
     
        irow(0) = 1

        do j = 1, n - 1
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
c    National Bureau of Standards, 1964,
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

      integer n_max
      parameter ( n_max = 11 )

      integer c
      integer c_vec(n_max)
      integer n
      integer n_data
      integer n_vec(n_max)

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
      subroutine charlier ( n, a, x, value )

c*********************************************************************72
c
cc CHARLIER evaluates Charlier polynomials at a point.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    17 March 2009
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c   J Simoes Pereira,
c    Algorithm 234: Poisson-Charliers Polynomials,
c    Communications of the ACM,
c    Volume 7, Number 7, page 420, July 1964.
c
c    Walter Gautschi,
c    Orthogonal Polynomials: Computation and Approximation,
c    Oxford, 2004,
c    ISBN: 0-19-850672-4,
c    LC: QA404.5 G3555.
c
c    Gabor Szego,
c    Orthogonal Polynomials,
c    American Mathematical Society, 1975,
c    ISBN: 0821810235,
c    LC: QA3.A5.v23.
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
c    Input, integer N, the maximum order of the polynomial.  
c    N must be at least 0.
c
c    Input, double precision A, the parameter.  A must not be 0.
c
c    Input, double precision X, the evaluation point.
c
c    Output, double precision VALUE(0:N), the value of the polynomials at X.
c
      implicit none

      integer n

      double precision a
      integer i
      double precision value(0:n)
      double precision x

      if ( a .eq. 0.0D+00 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'CHARLIER - Fatal error!'
        write ( *, '(a)' ) '  Parameter A cannot be zero.'
        stop
      end if

      if ( n .lt. 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'CHARLIER - Fatal error!'
        write ( *, '(a)' ) '  Parameter N must be nonnegative.'
        stop
      end if

      value(0) = 1.0D+00

      if ( n == 0 ) then
        return
      end if

      value(1) = - x / a

      if ( n == 1 ) then 
        return
      end if
  
      do i = 1, n - 1
        value(i+1) = ( ( dble ( i ) + a - x ) * value(i) 
     &                 - dble ( i ) * value(i-1) ) / a
      end do

      return
      end
      subroutine cheby_t_poly ( m, n, x, cx )

c*********************************************************************72
c
cc CHEBY_T_POLY evaluates Chebyshev polynomials T(n,x).
c
c  Discussion:
c
c    Chebyshev polynomials are useful as a basis for representing the
c    approximation of functions since they are well conditioned, in the sense
c    that in the interval [-1,1] they each have maximum absolute value 1.
c    Hence an error in the value of a coefficient of the approximation, of
c    size epsilon, is exactly reflected in an error of size epsilon between
c    the computed approximation and the theoretical approximation.
c
c    Typical usage is as follows, where we assume for the moment
c    that the interval of approximation is [-1,1].  The value
c    of N is chosen, the highest polynomial to be used in the
c    approximation.  Then the function to be approximated is
c    evaluated at the N+1 points XJ which are the zeroes of the N+1-th
c    Chebyshev polynomial.  Let these values be denoted by F(XJ).
c
c    The coefficients of the approximation are now defined by
c
c      C(I) = 2/(N+1) * sum ( 1 <= J <= N+1 ) F(XJ) T(I),XJ)
c
c    except that C(0) is given a value which is half that assigned
c    to it by the above formula,
c
c    and the representation is
c
c    F(X) approximated by sum ( 0 <= J <= N ) C(J) T(J,X)
c
c    Now note that, again because of the fact that the Chebyshev polynomials
c    have maximum absolute value 1, if the higher order terms of the
c    coefficients C are small, then we have the option of truncating
c    the approximation by dropping these terms, and we will have an
c    exact value for maximum perturbation to the approximation that
c    this will cause.
c
c    It should be noted that typically the error in approximation
c    is dominated by the first neglected basis function (some multiple of
c    T(N+1,X) in the example above).  If this term were the exact error,
c    then we would have found the minimax polynomial, the approximating
c    polynomial of smallest maximum deviation from the original function.
c    The minimax polynomial is hard to compute, and another important
c    feature of the Chebyshev approximation is that it tends to behave
c    like the minimax polynomial while being easy to compute.
c
c    To evaluate a sum like 
c
c      sum ( 0 <= J <= N ) C(J) T(J,X), 
c
c    Clenshaw's recurrence formula is recommended instead of computing the
c    polynomial values, forming the products and summing.
c
c    Assuming that the coefficients C(J) have been computed
c    for J = 0 to N, then the coefficients of the representation of the
c    indefinite integral of the function may be computed by
c
c      B(I) = ( C(I-1) - C(I+1))/2*(I-1) for I=1 to N+1, 
c
c    with
c 
c      C(N+1)=0
c      B(0) arbitrary.  
c
c    Also, the coefficients of the representation of the derivative of the 
c    function may be computed by:
c
c      D(I) = D(I+2)+2*I*C(I) for I=N-1, N-2, ..., 0, 
c
c    with
c
c      D(N+1) = D(N)=0.
c
c    Some of the above may have to adjusted because of the irregularity of C(0).
c
c    The formula is:
c
c      T(N,X) = COS(N*ARCCOS(X))
c
c  Differential equation:
c
c    (1-X*X) Y'' - X Y' + N N Y = 0
c
c  First terms:
c
c    T(0,X) =  1
c    T(1,X) =  1 X
c    T(2,X) =  2 X^2 -   1
c    T(3,X) =  4 X^3 -   3 X
c    T(4,X) =  8 X^4 -   8 X^2 +  1
c    T(5,X) = 16 X^5 -  20 X^3 +  5 X
c    T(6,X) = 32 X^6 -  48 X^4 + 18 X^2 - 1
c    T(7,X) = 64 X^7 - 112 X^5 + 56 X^3 - 7 X
c
c  Inequality:
c
c    abs ( T(N,X) ) <= 1 for -1 <= X <= 1
c
c  Orthogonality:
c
c    For integration over [-1,1] with weight
c
c      W(X) = 1 / sqrt(1-X*X), 
c
c    if we write the inner product of T(I,X) and T(J,X) as
c
c      < T(I,X), T(J,X) > = integral ( -1 <= X <= 1 ) W(X) T(I,X) T(J,X) dX
c
c    then the result is:
c
c      < T(I,X), T(J,X) > = 0    if I /= J
c      < T(I,X), T(J,X) > = PI/2 if I == J /= 0
c      < T(I,X), T(J,X) > = PI   if I == J == 0
c
c    A discrete orthogonality relation is also satisfied at each of
c    the N zeroes of T(N,X):  sum ( 1 <= K <= N ) T(I,X) * T(J,X)
c                              = 0 if I /= J
c                              = N/2 if I == J /= 0
c                              = N if I == J == 0
c
c  Recursion:
c
c    T(0,X) = 1,
c    T(1,X) = X,
c    T(N,X) = 2 * X * T(N-1,X) - T(N-2,X)
c
c    T'(N,X) = N * ( -X * T(N,X) + T(N-1,X) ) / ( 1 - X^2 )
c
c  Special values:
c
c    T(N,1) = 1
c    T(N,-1) = (-1)^N
c    T(2N,0) = (-1)^N
c    T(2N+1,0) = 0
c    T(N,X) = (-1)^N * T(N,-X)
c
c  Zeroes:
c
c    M-th zero of T(N,X) is X = cos((2*M-1)*PI/(2*N)), M = 1 to N.
c
c  Extrema:
c
c    M-th extremum of T(N,X) is X = cos(PI*M/N), M = 0 to N.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    28 March 2012
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
c  Parameters:
c
c    Input, integer M, the number of evaluation points.
c
c    Input, integer N, the highest polynomial to compute.
c
c    Input, double precision X(M), the evaluation points.
c
c    Output, double precision CX(M,0:N), the values of the N+1 
c    Chebyshev polynomials.
c
      implicit none

      integer m
      integer n

      double precision cx(m,0:n)
      integer i
      integer j
      double precision x(m)

      if ( n .lt. 0 ) then
        return
      end if

      do i = 1, m
        cx(i,0) = 1.0D+00
      end do

      if ( n .lt. 1 ) then
        return
      end if

      do i = 1, m
        cx(i,1) = x(i)
      end do

      do j = 2, n
        do i = 1, m
          cx(i,j) = 2.0D+00 * x(i) * cx(i,j-1) - cx(i,j-2)
        end do
      end do
     
      return
      end
      subroutine cheby_t_poly_coef ( n, c )

c*********************************************************************72
c
cc CHEBY_T_POLY_COEF evaluates coefficients of Chebyshev polynomials T(n,x).
c
c  First terms:
c
c    N/K     0     1      2      3       4     5      6    7      8    9   10
c
c     0      1
c     1      0     1
c     2     -1     0      2
c     3      0    -3      0      4
c     4      1     0     -8      0       8
c     5      0     5      0    -20       0    16
c     6     -1     0     18      0     -48     0     32
c     7      0    -7      0     56       0  -112      0    64
c
c  Recursion:
c
c    T(0,X) = 1,
c    T(1,X) = X,
c    T(N,X) = 2 * X * T(N-1,X) - T(N-2,X)
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    04 July 2008
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
c  Parameters:
c
c    Input, integer N, the highest order polynomial to compute.
c    Note that polynomials 0 through N will be computed.
c
c    Output, double precision C(0:N,0:N), the coefficients of the Chebyshev T
c    polynomials.
c
      implicit none

      integer n

      double precision c(0:n,0:n)
      integer i
      integer j

      if ( n .lt. 0 ) then
        return
      end if

      do j = 0, n
        do i = 0, n
          c(i,j) = 0.0D+00
        end do
      end do

      c(0,0) = 1.0D+00

      if ( n == 0 ) then
        return
      end if

      c(1,1) = 1.0D+00
 
      do i = 2, n
        c(i,0)     =                        - c(i-2,0)
        do j = 1, i - 2
          c(i,j) = 2.0D+00 * c(i-1,j-1) - c(i-2,j-1)
        end do
        c(i,  i-1) = 2.0D+00 * c(i-1,  i-2)
        c(i,  i  ) = 2.0D+00 * c(i-1,  i-1)
      end do
 
      return
      end
      subroutine cheby_t_poly_values ( n_data, n, x, fx )

c*********************************************************************72
c
cc CHEBY_T_POLY_VALUES returns values of Chebyshev polynomials T(n,x).
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    21 March 2007
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
c    Output, integer N, the order of the function.
c
c    Output, double precision X, the point where the function is evaluated.
c
c    Output, double precision FX, the value of the function.
c
      implicit none

      integer n_max
      parameter ( n_max = 13 )

      double precision fx
      double precision fx_vec(n_max)
      integer n
      integer n_data
      integer n_vec(n_max)
      double precision x
      double precision x_vec(n_max)

      save fx_vec
      save n_vec
      save x_vec

      data fx_vec /
     &   0.1000000000000000D+01,
     &   0.8000000000000000D+00,
     &   0.2800000000000000D+00,
     &  -0.3520000000000000D+00,
     &  -0.8432000000000000D+00,
     &  -0.9971200000000000D+00,
     &  -0.7521920000000000D+00,
     &  -0.2063872000000000D+00,
     &   0.4219724800000000D+00,
     &   0.8815431680000000D+00,
     &   0.9884965888000000D+00,
     &   0.7000513740800000D+00,
     &   0.1315856097280000D+00 /
      data n_vec /
     &   0,  1,  2,
     &   3,  4,  5,
     &   6,  7,  8,
     &   9, 10, 11,
     &  12 /
      data x_vec /
     &  0.8D+00,
     &  0.8D+00,
     &  0.8D+00,
     &  0.8D+00,
     &  0.8D+00,
     &  0.8D+00,
     &  0.8D+00,
     &  0.8D+00,
     &  0.8D+00,
     &  0.8D+00,
     &  0.8D+00,
     &  0.8D+00,
     &  0.8D+00 /

      if ( n_data .lt. 0 ) then
        n_data = 0
      end if

      n_data = n_data + 1

      if ( n_max .lt. n_data ) then
        n_data = 0
        n = 0
        x = 0.0D+00
        fx = 0.0D+00
      else
        n = n_vec(n_data)
        x = x_vec(n_data)
        fx = fx_vec(n_data)
      end if

      return
      end
      subroutine cheby_t_poly_zero ( n, z )

c*********************************************************************72
c
cc CHEBY_T_POLY_ZERO returns zeroes of Chebyshev polynomials T(n,x).
c
c  Discussion:
c
c    The I-th zero of T(N,X) is cos((2*I-1)*PI/(2*N)), I = 1 to N
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    04 July 2008
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the order of the polynomial.
c
c    Output, double precision Z(N), the zeroes of T(N,X).
c
      implicit none

      integer n

      double precision angle
      integer i
      double precision pi
      parameter ( pi = 3.141592653589793D+00 )
      double precision z(n)

      do i = 1, n
        angle = dble ( 2 * i - 1 ) * pi / dble ( 2 * n )
        z(i) = cos ( angle )
      end do

      return
      end
      subroutine cheby_u_poly ( n, x, cx )

c*********************************************************************72
c
cc CHEBY_U_POLY evaluates Chebyshev polynomials U(n,x).
c
c  Differential equation:
c
c    (1-X*X) Y'' - 3 X Y' + N (N+2) Y = 0
c
c    The formula is:
c
c      If |X| <= 1, then
c
c        U(N,X) = sin ( (N+1) * arccos(X) ) / sqrt ( 1 - X^2 )
c                = sin ( (N+1) * arccos(X) ) / sin ( arccos(X) )
c
c      else
c
c        U(N,X) = sinh ( (N+1) * arccosh(X) ) / sinh ( arccosh(X) )
c
c  First terms:
c
c    U(0,X) =   1
c    U(1,X) =   2 X
c    U(2,X) =   4 X^2 -   1
c    U(3,X) =   8 X^3 -   4 X
c    U(4,X) =  16 X^4 -  12 X^2 +  1
c    U(5,X) =  32 X^5 -  32 X^3 +  6 X
c    U(6,X) =  64 X^6 -  80 X^4 + 24 X^2 - 1
c    U(7,X) = 128 X^7 - 192 X^5 + 80 X^3 - 8X
c
c  Orthogonality:
c
c    For integration over [-1,1] with weight
c
c      W(X) = sqrt(1-X*X), 
c
c    we have
c
c      < U(I,X), U(J,X) > = integral ( -1 <= X <= 1 ) W(X) U(I,X) U(J,X) dX 
c
c    then the result is:
c
c      < U(I,X), U(J,X) >  =  0    if I /= J
c      < U(I,X), U(J,X) >  =  PI/2 if I == J
c
c  Recursion:
c
c    U(0,X) = 1,
c    U(1,X) = 2 * X,
c    U(N,X) = 2 * X * U(N-1,X) - U(N-2,X)
c
c  Special values:
c
c    U(N,1) = N + 1
c    U(2N,0) = (-1)^N
c    U(2N+1,0) = 0
c    U(N,X) = (-1)^N * U(N,-X)
c
c  Zeroes:
c
c    M-th zero of U(N,X) is X = cos( M*PI/(N+1)), M = 1 to N
c
c  Extrema:
c
c    M-th extremum of U(N,X) is X = cos( M*PI/N), M = 0 to N
c
c  Norm:
c
c    Integral ( -1 <= X <= 1 ) ( 1 - X^2 ) * U(N,X)^2 dX = PI/2
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    01 October 2002
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
c  Parameters:
c
c    Input, integer N, the highest polynomial to compute.
c
c    Input, double precision X, the point at which the polynomials 
c    are to be computed.
c
c    Output, double precision CX(0:N), the values of the N+1 Chebyshev
c    polynomials.
c
      implicit none

      integer n

      double precision cx(0:n)
      integer i
      double precision x

      if ( n .lt. 0 ) then
        return
      end if

      cx(0) = 1.0D+00

      if ( n .lt. 1 ) then
       return
      end if

      cx(1) = 2.0D+00 * x

      do i = 2, n
        cx(i) = 2.0D+00 * x * cx(i-1) - cx(i-2)
      end do
 
      return
      end
      subroutine cheby_u_poly_coef ( n, c )

c*********************************************************************72
c
cc CHEBY_U_POLY_COEF evaluates coefficients of Chebyshev polynomials U(n,x).
c
c  First terms:
c
c    N/K     0     1      2      3       4     5      6    7      8    9   10
c
c     0      1
c     1      0     2
c     2     -1     0      4
c     3      0    -4      0      8
c     4      1     0    -12      0      16
c     5      0     6      0    -32       0    32
c     6     -1     0     24      0     -80     0     64
c     7      0    -8      0     80       0  -192      0   128
c
c  Recursion:
c
c    U(0,X) = 1,
c    U(1,X) = 2*X,
c    U(N,X) = 2 * X * U(N-1,X) - U(N-2,X)
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    07 July 2008
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
c  Parameters:
c
c    Input, integer N, the highest order polynomial to compute.
c    Note that polynomials 0 through N will be computed.
c
c    Output, double precision C(0:N,0:N), the coefficients of the Chebyshev U
c    polynomials.
c
      implicit none

      integer n

      double precision c(0:n,0:n)
      integer i
      integer j

      if ( n .lt. 0 ) then
        return
      end if

      do j = 0, n
        do i = 0, n
          c(i,j) = 0.0D+00
        end do
      end do

      c(0,0) = 1.0D+00

      if ( n .eq. 0 ) then
        return
      end if

      c(1,1) = 2.0D+00
 
      do i = 2, n
        c(i,0)     =                        - c(i-2,0)
        do j = 1, i - 2
          c(i,j) = 2.0D+00 * c(i-1,j-1) - c(i-2,j)
        end do
        c(i,  i-1) = 2.0D+00 * c(i-1,  i-2)
        c(i,  i  ) = 2.0D+00 * c(i-1,  i-1)
      end do
 
      return
      end
      subroutine cheby_u_poly_values ( n_data, n, x, fx )

c*********************************************************************72
c
cc CHEBY_U_POLY_VALUES returns values of Chebyshev polynomials U(n,x).
c
c  Discussion:
c
c    In Mathematica, the function can be evaluated by:
c
c      ChebyshevU[n,x]
c
c    The Chebyshev U polynomial is a solution to the differential equation:
c
c    (1-X*X) Y'' - 3 X Y' + N (N+2) Y = 0
c
c  First terms:
c
c    U(0,X) =   1
c    U(1,X) =   2 X
c    U(2,X) =   4 X^2 -   1
c    U(3,X) =   8 X^3 -   4 X
c    U(4,X) =  16 X^4 -  12 X^2 +  1
c    U(5,X) =  32 X^5 -  32 X^3 +  6 X
c    U(6,X) =  64 X^6 -  80 X^4 + 24 X^2 - 1
c    U(7,X) = 128 X^7 - 192 X^5 + 80 X^3 - 8X
c
c  Recursion:
c
c    U(0,X) = 1,
c    U(1,X) = 2 * X,
c    U(N,X) = 2 * X * U(N-1,X) - U(N-2,X)
c
c  Norm:
c
c    Integral ( -1 <= X <= 1 ) ( 1 - X^2 ) * U(N,X)^2 dX = PI/2
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    25 April 2012
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
c    Output, integer N, the order of the function.
c
c    Output, double precision X, the point where the function is evaluated.
c
c    Output, double precision FX, the value of the function.
c
      implicit none

      integer n_max
      parameter ( n_max = 13 )

      double precision fx
      double precision fx_vec(n_max)
      integer n
      integer n_data
      integer n_vec(n_max)
      double precision x
      double precision x_vec(n_max)

      save fx_vec
      save n_vec
      save x_vec

      data fx_vec /
     &   0.1000000000000000D+01,
     &   0.1600000000000000D+01,
     &   0.1560000000000000D+01,
     &   0.8960000000000000D+00,
     &  -0.1264000000000000D+00,
     &  -0.1098240000000000D+01,
     &  -0.1630784000000000D+01,
     &  -0.1511014400000000D+01,
     &  -0.7868390400000000D+00,
     &   0.2520719360000000D+00,
     &   0.1190154137600000D+01,
     &   0.1652174684160000D+01,
     &   0.1453325357056000D+01 /
      data n_vec /
     &   0,  1,  2,
     &   3,  4,  5,
     &   6,  7,  8,
     &   9, 10, 11,
     &  12 /
      data x_vec /
     &  0.8D+00,
     &  0.8D+00,
     &  0.8D+00,
     &  0.8D+00,
     &  0.8D+00,
     &  0.8D+00,
     &  0.8D+00,
     &  0.8D+00,
     &  0.8D+00,
     &  0.8D+00,
     &  0.8D+00,
     &  0.8D+00,
     &  0.8D+00 /

      if ( n_data .lt. 0 ) then
        n_data = 0
      end if

      n_data = n_data + 1

      if ( n_max .lt. n_data ) then
        n_data = 0
        n = 0
        x = 0.0D+00
        fx = 0.0D+00
      else
        n = n_vec(n_data)
        x = x_vec(n_data)
        fx = fx_vec(n_data)
      end if

      return
      end
      subroutine cheby_u_poly_zero ( n, z )

c*********************************************************************72
c
cc CHEBY_U_POLY_ZERO returns zeroes of Chebyshev polynomials U(n,x).
c
c  Discussion:
c
c    The I-th zero of U(N,X) is cos((I-1)*PI/(N-1)), I = 1 to N
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    07 July 2008
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the order of the polynomial.
c
c    Output, double precision Z(N), the zeroes of U(N,X).
c
      implicit none

      integer n

      double precision angle
      integer i
      double precision pi
      parameter ( pi = 3.141592653589793D+00 )
      double precision z(n)

      do i = 1, n
        angle = dble ( i ) * pi / dble ( n + 1 )
        z(i) = cos ( angle )
      end do

      return
      end
      subroutine chebyshev_discrete ( n, m, x, v )

c*********************************************************************72
c
cc CHEBYSHEV_DISCRETE evaluates discrete Chebyshev polynomials at a point.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    18 March 2009
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Walter Gautschi,
c    Orthogonal Polynomials: Computation and Approximation,
c    Oxford, 2004,
c    ISBN: 0-19-850672-4,
c    LC: QA404.5 G3555.
c
c  Parameters:
c
c    Input, integer N, the highest order of the polynomials to 
c    be evaluated.  0 <= N <= M.
c
c    Input, integer M, the maximum order of the polynomials.
c    0 <= M.
c
c    Input, double precision X, the evaluation point.
c
c    Output, double precision V(0:N), the value of the polynomials at X.
c
      implicit none

      integer n

      integer i
      integer m
      double precision x
      double precision v(0:n)

      if ( m .lt. 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'CHEBYSHEV_DISCRETE - Fatal error!'
        write ( *, '(a)' ) '  Parameter M must be nonnegative.'
        stop
      end if

      if ( n .lt. 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'CHEBYSHEV_DISCRETE - Fatal error!'
        write ( *, '(a)' ) '  Parameter N must be nonnegative.'
        stop
      end if

      if ( m .lt. n ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'CHEBYSHEV_DISCRETE - Fatal error!'
        write ( *, '(a)' ) '  Parameter N must be no greater than M.'
        stop
      end if

      v(0) = 1.0D+00

      if ( n .eq. 0 ) then
        return
      end if

      v(1) = 2.0D+00 * x + dble ( 1 - m )

      if ( n .eq. 1 ) then
        return
      end if

      do i = 1, n - 1
        v(i+1) = ( 
     &    dble ( 2 * i + 1 ) 
     &    * ( 2.0D+00 * x + dble ( 1 - m ) ) * v(i) 
     &    - dble ( i * ( m + i ) * ( m - i ) ) * v(i-1) 
     &  ) / dble ( i + 1 )
      end do

      return
      end
      function collatz_count ( n )

c*****************************************************************************80
c
cc COLLATZ_COUNT counts the number of terms in a Collatz sequence.
c
c  Discussion:
c
c    The rules for generation of the Collatz sequence are recursive.
c    If T is the current entry of the sequence, (T is
c    assumed to be a positive integer), then the next
c    entry, U is determined as follows:
c    
c      if T is 1 (or less)
c        terminate the sequence;
c      else if T is even
c        U = T/2.
c      else (if T is odd and not 1)
c        U = 3*T+1;
c
c     N  Sequence                                                Length
c
c     1                                                               1
c     2   1                                                           2
c     3  10,  5, 16,  8,  4,  2,  1                                   8
c     4   2   1                                                       3
c     5  16,  8,  4,  2,  1                                           6
c     6   3, 10,  5, 16,  8,  4,  2,  1                               9
c     7  22, 11, 34, 17, 52, 26, 13, 40, 20, 10, 5, 16, 8, 4, 2, 1   17
c     8   4,  2,  1                                                   4
c     9  28, 14,  7, ...                                             20
c    10   5, 16,  8,  4,  2,  1                                       7
c    11  34, 17, 52, 26, 13, 40, 20, 10,  5, 16, 8, 4, 2, 1          15
c    12   6,  3, 10,  5, 16,  8,  4,  2,  1                          10
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    03 July 2008
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Eric Weisstein,
c    CRC Concise Encyclopedia of Mathematics,
c    CRC Press, 2002,
c    Second edition,
c    ISBN: 1584883472,
c    LC: QA5.W45
c
c  Parameters:
c
c    Input, integer N, the first element of the sequence.
c
c    Output, integer COLLATZ_COUNT, the number of elements in
c    the Collatz sequence that begins with N.
c
      implicit none

      integer collatz_count
      integer count
      integer n
      integer n_local

      count = 1
      n_local = n

10    continue

        if ( n_local .le. 1 ) then
          go to 20
        else if ( mod ( n_local, 2 ) == 0 ) then
          n_local = n_local / 2
        else
          n_local = 3 * n_local + 1
        end if

        count = count + 1

      go to 10

20    continue

      collatz_count = count
 
      return
      end
      subroutine collatz_count_max ( n, i_max, j_max )

c*********************************************************************72
c
cc COLLATZ_COUNT_MAX seeks the maximum Collatz count for 1 through N.
c
c  Discussion:
c
c    For each integer I, we compute a sequence of values that 
c    terminate when we reach 1.  The number of steps required to
c    reach 1 is the "rank" of I, and we are searching the numbers
c    from 1 to N for the number with maximum rank.
c
c    For a given I, the sequence is produced by:
c
c    1) J = 1, X(J) = I;
c    2) If X(J) = 1, stop.
c    3) J = J + 1; 
c       if X(J-1) was even, X(J) = X(J-1)/2;
c       else                X(J) = 3 * X(J-1) + 1;
c    4) Go to 3
c
c  Example:
c
c            N     I_MAX J_MAX
c
c           10         9    20
c          100        97   119
c        1,000       871   179
c       10,000     6,171   262
c      100,000    77,031   351
c    1,000,000   837,799   525
c   10,000,000 8,400,511   686
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    12 April 2009
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the maximum integer to check.
c
c    Output, integer I_MAX, J_MAX, an integer I with the maximum 
c    rank, and the value of the maximum rank.
c
      implicit none

      integer i
      integer i_max
      integer j
      integer j_max
      integer n
      integer x

      i_max = -1
      j_max = -1

      do i = 1, n

        j = 1
        x = i

10      continue

        if ( x .ne. 1 ) then

          j = j + 1

          if ( mod ( x, 2 ) .eq. 0 ) then
            x = x / 2
          else
            x = 3 * x + 1
          end if

          go to 10

        end if

        if ( j_max .lt. j ) then
          i_max = i
          j_max = j
        end if

      end do

      return
      end
      subroutine collatz_count_values ( n_data, n, count )

c*********************************************************************72
c
cc COLLATZ_COUNT_VALUES returns some values of the Collatz count function.
c
c  Discussion:
c
c    The rules for generation of the Collatz sequence are recursive.
c    If T is the current entry of the sequence, (T is
c    assumed to be a positive integer), then the next
c    entry, U is determined as follows:
c
c      if T is 1 (or less)
c        terminate the sequence;
c      else if T is even
c        U = T/2.
c      else (if T is odd and not 1)
c        U = 3*T+1;
c
c    The Collatz count is the length of the Collatz sequence for a given
c    starting value.  By convention, we include the initial value in the
c    count, so the minimum value of the count is 1.
c
c     N  Sequence                                                 Count
c
c     1                                                               1
c     2   1                                                           2
c     3  10,  5, 16,  8,  4,  2,  1                                   8
c     4   2   1                                                       3
c     5  16,  8,  4,  2,  1                                           6
c     6   3, 10,  5, 16,  8,  4,  2,  1                               9
c     7  22, 11, 34, 17, 52, 26, 13, 40, 20, 10, 5, 16, 8, 4, 2, 1   17
c     8   4,  2,  1                                                   4
c     9  28, 14,  7, ...                                             20
c    10   5, 16,  8,  4,  2,  1                                       7
c    11  34, 17, 52, 26, 13, 40, 20, 10,  5, 16, 8, 4, 2, 1          15
c    12   6,  3, 10,  5, 16,  8,  4,  2,  1                          10
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    07 March 2006
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Eric Weisstein,
c    "The Collatz Problem",
c    CRC Concise Encyclopedia of Mathematics,
c    CRC 1998.
c
c  Parameters:
c
c    Input/output, integer N_DATA.  The user sets N_DATA to 0 before the
c    first call.  On each call, the routine increments N_DATA by 1, and
c    returns the corresponding data; when there is no more data, the
c    output value of N_DATA will be 0 again.
c
c    Output, integer N, the initial value of a Collatz sequence.
c
c    Output, integer COUNT, the length of the Collatz sequence starting
c    with N.
c
      implicit none

      integer n_max
      parameter ( n_max = 20 )

      integer count
      integer count_vec(n_max)
      integer n
      integer n_data
      integer n_vec(n_max)

      save count_vec
      save n_vec

      data count_vec /
     &     1,   2,   8,   3,   6,   9,   17,   4,  20,   7,
     &  112,  25,  26,  27,  17,  28,  111,  18,  83,  29 /
      data n_vec /
     &    1,   2,   3,   4,   5,   6,   7,   8,   9,  10,
     &   27,  50, 100, 200, 300, 400, 500, 600, 700, 800 /

      if ( n_data .lt. 0 ) then
        n_data = 0
      end if

      n_data = n_data + 1

      if ( n_max .lt. n_data ) then
        n_data = 0
        n = 0
        count = 0
      else
        n = n_vec(n_data)
        count = count_vec(n_data)
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
c    The formula is:
c
c      C(N,K) = Nc / ( Kc * (N-K)c )
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
c    03 December 2007
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer IDO, indicates whether this is a call for
c    the 'next' row of the triangle.
c    * 0 means this is a startup call.  Row N is desired, but
c      presumably this is a first call, or row N-1 was not computed
c      on the previous call.
c    * 1 means this is not the first call, and row N-1 was computed
c      on the previous call.  In this case, much work can be saved
c      by using the information from the previous values of IROW
c      to build the next values.
c
c    Input, integer N, the row of the triangle desired.  The 
c    triangle begins with row 0.
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

      if ( n < 0 ) then
        return
      end if
 
      if ( ido .eq. 1 ) then
 
        do i = n, 2, -1
          irow(i) = irow(i) + irow(i-1)
        end do
     
        irow(n+1) = 1
 
      else
 
        irow(1) = 1
        do i = 2, n+ 1
          irow(i) = 0
        end do

        do j = 1, n
          do i = j + 1, 2, -1
            irow(i) = irow(i) + irow(i-1)
          end do
        end do
 
      end if
 
      return
      end
      subroutine commul ( n, nfactor, factor, ncomb )

c*********************************************************************72
c
cc COMMUL computes a multinomial combinatorial coefficient.
c
c  Discussion:
c
c    The multinomial coefficient is a generalization of the binomial
c    coefficient.  It may be interpreted as the number of combinations of
c    N objects, where FACTOR(1) objects are indistinguishable of type 1,
c    ... and FACTOR(K) are indistinguishable of type NFACTOR.
c
c    The formula is:
c
c      NCOMB = N! / ( FACTOR(1)! FACTOR(2)! ... FACTOR(NFACTOR)! )
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    07 July 2008
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, determines the numerator.
c
c    Input, integer NFACTOR, the number of factors in the 
c    numerator.
c
c    Input, integer FACTOR(NFACTOR).
c    FACTOR contains the NFACTOR values used in the denominator.
c    Note that the sum of these entries should be N,
c    and that all entries should be nonnegative.
c
c    Output, integer NCOMB, the value of the multinomial 
c    coefficient.
c
      implicit none

      integer nfactor

      double precision arg
      double precision fack
      double precision facn
      integer factor(nfactor)
      integer i
      integer isum
      integer j
      integer n
      integer ncomb
      double precision r8_gamma_log

      if ( nfactor .lt. 1 ) then
        ncomb = 1
        return
      end if

      do i = 1, nfactor

        if ( factor(i) .lt. 0 ) then
          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) 'COMMUL - Fatal error!'
          write ( *, '(a,i8,a,i8)' ) 
     &      '  Entry ', I, ' of FACTOR = ', factor(i)
          write ( *, '(a)' ) '  But this value must be nonnegative.'
          stop
        end if

      end do
 
      isum = 0
      do j = 1, nfactor
        isum = isum + factor(j)
      end do

      if ( isum .ne. n ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'COMMUL - Fatal error!'
        write ( *, '(a,i8)' ) 
     &    '  The sum of the FACTOR entries is ', isum
        write ( *, '(a,i8)' ) '  But it must equal N = ', n
        stop
      end if
 
      arg = dble ( n + 1 )
      facn = r8_gamma_log ( arg )
 
      do i = 1, nfactor
 
        arg = dble ( factor(i) + 1 )
        fack = r8_gamma_log ( arg )
        facn = facn - fack
 
      end do
 
      ncomb = nint ( exp ( facn ) )
 
      return
      end
      function cos_deg ( angle )

c*********************************************************************72
c
cc COS_DEG returns the cosine of an angle given in degrees.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    04 July 2008
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, double precision ANGLE, the angle, in degrees.
c
c    Output, double precision COS_DEG, the cosine of the angle.
c
      implicit none

      double precision angle
      double precision cos_deg
      double precision degrees_to_radians
      parameter 
     &  ( degrees_to_radians = 3.141592653589793D+00 / 180.0D+00 )

      cos_deg = cos ( degrees_to_radians * angle )

      return
      end
      function cos_power_int ( a, b, n )

c*********************************************************************72
c
cc COS_POWER_INT evaluates the cosine power integral.
c
c  Discussion:
c
c    The function is defined by
c
c      COS_POWER_INT(A,B,N) = Integral ( A <= T <= B ) ( cos ( t ))^n dt
c
c    The algorithm uses the following fact:
c
c      Integral cos^n ( t ) = -(1/n) * (
c        cos^(n-1)(t) * sin(t) + ( n-1 ) * Integral cos^(n-2) ( t ) dt )
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    31 March 2012
c
c  Author:
c
c    John Burkardt
c
c  Parameters
c
c    Input, double precision A, B, the limits of integration.
c
c    Input, integer N, the power of the sine function.
c
c    Output, double precision COS_POWER_INT, the value of the integral.
c
      implicit none

      double precision a
      double precision b
      double precision ca
      double precision cb
      double precision cos_power_int
      integer m
      integer mlo
      integer n
      double precision sa
      double precision sb
      double precision value

      if ( n .lt. 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'COS_POWER_INT - Fatal error!'
        write ( *, '(a)' ) '  Power N < 0.'
        value = 0.0
        stop
      end if

      sa = sin ( a )
      sb = sin ( b )
      ca = cos ( a )
      cb = cos ( b )

      if ( mod ( n, 2 ) .eq. 0 ) then
        value = b - a
        mlo = 2
      else
        value = sb - sa
        mlo = 3
      end if

      do m = mlo, n, 2
        value = ( dble ( m - 1 ) * value 
     &            - ca**(m-1) * sa + cb**(m-1) * sb ) 
     &    / dble ( m )
      end do

      cos_power_int = value

      return
      end
      subroutine cos_power_int_values ( n_data, a, b, n, fx )

c*********************************************************************72
c
cc COS_POWER_INT_VALUES returns some values of the cosine power integral.
c
c  Discussion:
c
c    The function has the form
c
c      COS_POWER_INT(A,B,N) = Integral ( A <= T <= B ) ( cos(T) )^N dt
c
c    In Mathematica, the function can be evaluated by:
c
c      Integrate [ ( Cos[x] )^n, { x, a, b } ]
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    30 March 2012
c
c  Author:
c
c    John Burkardt
c
c  Reference:
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
c    Input/output, integer N_DATA.  The user sets N_DATA to 0
c    before the first call.  On each call, the routine increments N_DATA by 1,
c    and returns the corresponding data; when there is no more data, the
c    output value of N_DATA will be 0 again.
c
c    Output, double precision A, B, the limits of integration.
c
c    Output, integer ( kind = 4 ) N, the power.
c
c    Output, double precision FX, the function value.
c
      implicit none

      integer n_max
      parameter ( n_max = 11 )

      double precision a
      double precision a_vec(n_max)
      double precision b
      double precision b_vec(n_max)
      double precision fx
      double precision fx_vec(n_max)
      integer n
      integer n_data
      integer n_vec(n_max)

      save a_vec
      save b_vec
      save fx_vec
      save n_vec

      data a_vec /
     &   0.00D+00, 
     &   0.00D+00, 
     &   0.00D+00, 
     &   0.00D+00, 
     &   0.00D+00, 
     &   0.00D+00, 
     &   0.00D+00, 
     &   0.00D+00, 
     &   0.00D+00, 
     &   0.00D+00, 
     &   0.00D+00 /
      data b_vec /
     &   3.141592653589793D+00, 
     &   3.141592653589793D+00, 
     &   3.141592653589793D+00, 
     &   3.141592653589793D+00, 
     &   3.141592653589793D+00, 
     &   3.141592653589793D+00, 
     &   3.141592653589793D+00, 
     &   3.141592653589793D+00, 
     &   3.141592653589793D+00, 
     &   3.141592653589793D+00, 
     &   3.141592653589793D+00 /
      data fx_vec / 
     &   3.141592653589793D+00, 
     &   0.0D+00, 
     &   1.570796326794897D+00, 
     &   0.0D+00, 
     &   1.178097245096172D+00, 
     &   0.0D+00, 
     &   0.9817477042468104D+00, 
     &   0.0D+00, 
     &   0.8590292412159591D+00, 
     &   0.0D+00, 
     &   0.7731263170943632D+00 /
      data n_vec / 
     &   0, 
     &   1, 
     &   2, 
     &   3, 
     &   4, 
     &   5, 
     &   6, 
     &   7, 
     &   8, 
     &   9, 
     &  10 /

      if ( n_data .lt. 0 ) then
        n_data = 0
      end if

      n_data = n_data + 1

      if ( n_max .lt. n_data ) then
        n_data = 0
        a = 0.0D+00
        b = 0.0D+00
        n = 0
        fx = 0.0D+00
      else
        a = a_vec(n_data)
        b = b_vec(n_data)
        n = n_vec(n_data)
        fx = fx_vec(n_data)
      end if

      return
      end
      subroutine delannoy ( m, n, a )

c*********************************************************************72
c
cc DELANNOY returns the Delannoy numbers up to orders (M,N).
c
c  Discussion:
c
c    The Delannoy number A(M,N) counts the number of distinct paths
c    from (0,0) to (M,N) in which the only steps used are
c    (1,1), (1,0) and (0,1).
c
c  First values:
c
c      \N 0  1   2   3    4     5     6      7      8
c     M-+--------------------------------------------
c     0 | 1  1   1   1    1     1     1      1      1
c     1 | 1  3   5   7    9    11    13     15     17
c     2 | 1  5  13  25   41    61    85    113    145
c     3 | 1  7  25  63  129   231   377    575    833
c     4 | 1  9  41 129  321   681  1289   2241   3649
c     5 | 1 11  61 231  681  1683  3653   7183  13073
c     6 | 1 13  85 377 1289  3653  8989  19825  40081
c     7 | 1 15 113 575 2241  7183 19825  48639 108545
c     8 | 1 17 145 833 3649 13073 40081 108545 265729
c
c  Recursion:
c
c    A(0,0) = 1
c    A(M,0) = 1
c    A(0,N) = 1
c    A(M,N) = A(M-1,N) + A(M,N-1) + A(M-1,N-1)
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    03 July 2008
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Eric Weisstein,
c    CRC Concise Encyclopedia of Mathematics,
c    CRC Press, 2002,
c    Second edition,
c    ISBN: 1584883472,
c    LC: QA5.W45
c
c  Parameters:
c
c    Input, integer M, N, define the highest order number to 
c    compute.
c
c    Output, integer A(0:M,0:N), the Delannoy numbers.
c
      implicit none

      integer m
      integer n

      integer a(0:m,0:n)
      integer i
      integer j

      if ( m .lt. 0 ) then
        return
      end if

      if ( n .lt. 0 ) then
        return
      end if
   
      a(0,0) = 1

      do i = 1, m
        a(i,0) = 1
      end do

      do j = 1, n
        a(0,j) = 1
      end do

      do i = 1, m
        do j = 1, n
          a(i,j) = a(i-1,j) + a(i,j-1) + a(i-1,j-1)
        end do
      end do

      return
      end
      function e_constant ( )

c*********************************************************************72
c
cc E_CONSTANT returns the value of the base of the natural logarithm system.
c
c  Discussion:
c
c    E = Limit ( N -> +oo ) ( 1 + 1 / N )**N
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    03 July 2008
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, double precision E_CONSTANT, the base of the natural 
c    logarithm system.
c
      implicit none

      double precision e_constant

      e_constant = 2.718281828459045235360287D+00
 
      return
      end
      subroutine erf_values ( n_data, x, fx )

c*********************************************************************72
c
cc ERF_VALUES returns some values of the ERF or "error" function for testing.
c
c  Discussion:
c
c    The error function is defined by:
c
c      ERF(X) = ( 2 / sqrt ( PI ) * integral ( 0 <= T <= X ) exp ( - T^2 ) dT
c
c    In Mathematica, the function can be evaluated by:
c
c      Erf[x]
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    29 March 2007
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
c    Cambridge University Press, 1999,
c    ISBN: 0-521-64314-7,
c    LC: QA76.95.W65.
c
c  Parameters:
c
c    Input/output, integer N_DATA.
c    On input, if N_DATA is 0, the first test data is returned, and
c    N_DATA is set to the index of the test data.  On each subsequent
c    call, N_DATA is incremented and that test data is returned.  When
c    there is no more test data, N_DATA is set to 0.
c
c    Output, double precision X, the argument of the function.
c
c    Output, double precision FX, the value of the function.
c
      implicit none

      integer n_max
      parameter ( n_max = 21 )

      double precision bvec ( n_max )
      double precision fx
      integer n_data
      double precision x
      double precision xvec ( n_max )

      data bvec /
     & 0.0000000000000000D+00,
     & 0.1124629160182849D+00,
     & 0.2227025892104785D+00,
     & 0.3286267594591274D+00,
     & 0.4283923550466685D+00,
     & 0.5204998778130465D+00,
     & 0.6038560908479259D+00,
     & 0.6778011938374185D+00,
     & 0.7421009647076605D+00,
     & 0.7969082124228321D+00,
     & 0.8427007929497149D+00,
     & 0.8802050695740817D+00,
     & 0.9103139782296354D+00,
     & 0.9340079449406524D+00,
     & 0.9522851197626488D+00,
     & 0.9661051464753107D+00,
     & 0.9763483833446440D+00,
     & 0.9837904585907746D+00,
     & 0.9890905016357307D+00,
     & 0.9927904292352575D+00,
     & 0.9953222650189527D+00 /
      data xvec /
     &  0.0D+00,
     &  0.1D+00,
     &  0.2D+00,
     &  0.3D+00,
     &  0.4D+00,
     &  0.5D+00,
     &  0.6D+00,
     &  0.7D+00,
     &  0.8D+00,
     &  0.9D+00,
     &  1.0D+00,
     &  1.1D+00,
     &  1.2D+00,
     &  1.3D+00,
     &  1.4D+00,
     &  1.5D+00,
     &  1.6D+00,
     &  1.7D+00,
     &  1.8D+00,
     &  1.9D+00,
     &  2.0D+00 /

      if ( n_data .lt. 0 ) then
        n_data = 0
      end if

      n_data = n_data + 1

      if ( n_max .lt. n_data ) then
        n_data = 0
        x = 0.0D+00
        fx = 0.0D+00
      else
        x = xvec(n_data)
        fx = bvec(n_data)
      end if

      return
      end
      function error_f ( x )

c*********************************************************************72
c
cc ERROR_F evaluates the error function ERF(X).
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    04 July 2008
c
c  Author:
c
c    Original FORTRAN77 version by William Cody.
c    Modifications by John Burkardt.
c
c  Reference:
c
c    William Cody,
c    Rational Chebyshev approximations for the error function,
c    Mathematics of Computation, 
c    1969, pages 631-638.
c
c  Parameters:
c
c    Input, double precision X, the argument of ERF.
c
c    Output, double precision ERROR_F, the value of ERF(X).
c
      implicit none

      double precision a(5)
      double precision b(4)
      double precision c(9)
      double precision d(8)
      double precision del
      double precision error_f
      integer i
      double precision p(6)
      double precision q(5)
      double precision r8_epsilon
      double precision sqrpi
      parameter ( sqrpi = 0.56418958354775628695D+00 )
      double precision thresh
      parameter ( thresh = 0.46875D+00 )
      double precision x
      double precision xabs
      double precision xbig
      parameter ( xbig = 26.543D+00 )
      double precision xden
      double precision xnum
      double precision xsq

      save a
      save b
      save c
      save d
      save p
      save q

      data a /
     &  3.16112374387056560D+00, 
     &  1.13864154151050156D+02, 
     &  3.77485237685302021D+02, 
     &  3.20937758913846947D+03, 
     &  1.85777706184603153D-01 /
      data b /
     &  2.36012909523441209D+01, 
     &  2.44024637934444173D+02, 
     &  1.28261652607737228D+03, 
     &  2.84423683343917062D+03 /
      data c /
     &  5.64188496988670089D-01, 
     &  8.88314979438837594D+00, 
     &  6.61191906371416295D+01, 
     &  2.98635138197400131D+02, 
     &  8.81952221241769090D+02, 
     &  1.71204761263407058D+03, 
     &  2.05107837782607147D+03, 
     &  1.23033935479799725D+03, 
     &  2.15311535474403846D-08 /
      data d /
     &  1.57449261107098347D+01, 
     &  1.17693950891312499D+02, 
     &  5.37181101862009858D+02, 
     &  1.62138957456669019D+03, 
     &  3.29079923573345963D+03, 
     &  4.36261909014324716D+03, 
     &  3.43936767414372164D+03, 
     &  1.23033935480374942D+03 /
      data p /
     &  3.05326634961232344D-01, 
     &  3.60344899949804439D-01, 
     &  1.25781726111229246D-01, 
     &  1.60837851487422766D-02, 
     &  6.58749161529837803D-04, 
     &  1.63153871373020978D-02 /
      data q /
     &  2.56852019228982242D+00, 
     &  1.87295284992346047D+00, 
     &  5.27905102951428412D-01, 
     &  6.05183413124413191D-02, 
     &  2.33520497626869185D-03 /

      xabs = abs ( x )
c
c  Evaluate ERF(X) for |X| <= 0.46875.
c
      if ( xabs .le. thresh ) then

        if ( r8_epsilon ( ) .lt. xabs ) then
          xsq = xabs * xabs
        else
          xsq = 0.0D+00
        end if

        xnum = a(5) * xsq
        xden = xsq
        do i = 1, 3
          xnum = ( xnum + a(i) ) * xsq
          xden = ( xden + b(i) ) * xsq
        end do

        error_f = x * ( xnum + a(4) ) / ( xden + b(4) )
c
c  Evaluate ERFC(X) for 0.46875 <= |X| <= 4.0.
c
      else if ( xabs .le. 4.0D+00 ) then

        xnum = c(9) * xabs
        xden = xabs
        do i = 1, 7
          xnum = ( xnum + c(i) ) * xabs
          xden = ( xden + d(i) ) * xabs
        end do

        error_f = ( xnum + c(8) ) / ( xden + d(8) )
        xsq = aint ( xabs * 16.0D+00 ) / 16.0D+00
        del = ( xabs - xsq ) * ( xabs + xsq )
        error_f = exp ( - xsq * xsq ) * exp ( - del ) * error_f

        error_f = ( 0.5D+00 - error_f ) + 0.5D+00

        if ( x .lt. 0.0D+00 ) then
          error_f = - error_f
        end if
c
c  Evaluate ERFC(X) for 4.0 < |X|.
c
      else

        if ( xbig .le. xabs ) then

          if ( 0.0D+00 .lt. x ) then
            error_f = 1.0D+00
          else
            error_f = -1.0D+00
          end if

        else

          xsq = 1.0D+00 / ( xabs * xabs )

          xnum = p(6) * xsq
          xden = xsq
          do i = 1, 4
            xnum = ( xnum + p(i) ) * xsq
            xden = ( xden + q(i) ) * xsq
          end do

          error_f = xsq * ( xnum + p(5) ) / ( xden + q(5) )
          error_f = ( sqrpi - error_f ) / xabs
          xsq = aint ( xabs * 16.0D+00 ) / 16.0D+00
          del = ( xabs - xsq ) * ( xabs + xsq )
          error_f = exp ( - xsq * xsq ) * exp ( - del ) * error_f

          error_f = ( 0.5D+00 - error_f ) + 0.5D+00
          if ( x .lt. 0.0D+00 ) then
            error_f = - error_f
          end if

        end if

      end if

      return
      end
      function error_f_inverse ( y )

c*********************************************************************72
c
cc ERROR_F_INVERSE inverts the error function ERF.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    05 August 2010
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, double precision Y, the value of the error function.
c
c    Output, double precision ERROR_F_INVERSE, the value X such that
c    ERROR_F(X) = Y.
c
      implicit none

      double precision error_f_inverse
      double precision normal_01_cdf_inv
      double precision x
      double precision y
      double precision z

      z = ( y + 1.0D+00 ) / 2.0D+00

      x = normal_01_cdf_inv ( z )

      error_f_inverse = x / sqrt ( 2.0D+00 )

      return
      end
      function euler_constant ( )

c*********************************************************************72
c
cc EULER_CONSTANT returns the value of the Euler-Mascheroni constant.
c
c  Discussion:
c
c    The Euler-Mascheroni constant is often denoted by a lower-case
c    Gamma.  Gamma is defined as
c
c      Gamma = limit ( M -> +oo )
c        ( sum ( 1 <= N <= M ) 1 / N ) - log ( M )
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    03 July 2008
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, double precision EULER_CONSTANT, the value of the 
c    Euler-Mascheroni constant.
c
      implicit none

      double precision euler_constant

      euler_constant = 0.577215664901532860606512090082402431042D+00

      return
      end
      subroutine euler_number ( n, e )

c*********************************************************************72
c
cc EULER_NUMBER computes the Euler numbers.
c
c  Discussion:
c
c    The Euler numbers can be evaluated in Mathematica by:
c
c      EulerE[n]
c
c    These numbers rapidly get too big to store in an ordinary integer!
c
c    The terms of odd index are 0.
c
c    E(N) = -C(N,N-2) * E(N-2) - C(N,N-4) * E(N-4) - ... - C(N,0) * E(0).
c
c  First terms:
c
c    E0  = 1
c    E1  = 0
c    E2  = -1
c    E3  = 0
c    E4  = 5
c    E5  = 0
c    E6  = -61
c    E7  = 0
c    E8  = 1385
c    E9  = 0
c    E10 = -50521
c    E11 = 0
c    E12 = 2702765
c    E13 = 0
c    E14 = -199360981
c    E15 = 0
c    E16 = 19391512145
c    E17 = 0
c    E18 = -2404879675441
c    E19 = 0
c    E20 = 370371188237525
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    05 July 2008
c
c  Author:
c
c    John Burkardt
c
c  Reference:
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
c    Input, integer N, the index of the last Euler number
c    to compute.
c
c    Output, integer E(0:N), the Euler numbers.
c
      implicit none

      integer n

      integer e(0:n)
      integer i
      integer i4_choose
      integer j

      if ( n .lt. 0 ) then
        return
      end if

      e(0) = 1

      if ( n .eq. 0 ) then
        return
      end if

      e(1) = 0
 
      if ( n .eq. 1 ) then
        return
      end if

      e(2) = -1
    
      do i = 3, n

        e(i) = 0

        if ( mod ( i, 2 ) .eq. 0 ) then

          do j = 2, i, 2
            e(i) = e(i) - i4_choose ( i, j ) * e(i-j)
          end do

        end if

      end do
 
      return
      end
      function euler_number2 ( n )

c*********************************************************************72
c
cc EULER_NUMBER2 computes the Euler numbers.
c
c  Discussion:
c
c    The Euler numbers can be evaluated in Mathematica by:
c
c      EulerE[n]
c
c  First terms:
c
c    E0  = 1
c    E1  = 0
c    E2  = -1
c    E3  = 0
c    E4  = 5
c    E5  = 0
c    E6  = -61
c    E7  = 0
c    E8  = 1385
c    E9  = 0
c    E10 = -50521
c    E11 = 0
c    E12 = 2702765
c    E13 = 0
c    E14 = -199360981
c    E15 = 0
c    E16 = 19391512145
c    E17 = 0
c    E18 = -2404879675441
c    E19 = 0
c    E20 = 370371188237525
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    08 July 2008
c
c  Author:
c
c    John Burkardt
c
c  Reference:
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
c    Input, integer N, the index of the Euler number to compute.
c
c    Output, double precision EULER_NUMBER2, the value of E(N).
c
      implicit none

      double precision euler_number2
      double precision e(0:6)
      integer i
      integer itmax
      parameter ( itmax = 1000 )
      integer n
      double precision pi
      parameter ( pi = 3.141592653589793D+00 )
      double precision r8_factorial
      double precision sum1
      double precision term

      save e

      data e /
     &  1.0D+00, -1.0D+00, 5.0D+00, -61.0D+00, 1385.0D+00, 
     &  -50521.0D+00, 2702765.0D+00 /

      if ( n .lt. 0 ) then
        euler_number2 = 0.0D+00
        return
      end if

      if ( n .eq. 0 ) then
        euler_number2 = e(0)
        return
      end if

      if ( mod ( n, 2 ) .eq. 1 ) then
        euler_number2 = 0.0D+00
        return
      end if

      if ( n .le. 12 ) then
        euler_number2 = e(n/2)
        return
      end if

      sum1 = 0.0D+00
      do i = 1, itmax

        term = 1.0D+00 / dble ( ( 2 * i - 1 )**( n + 1 ) )

        if ( mod ( i, 2 ) .eq. 1 ) then
          sum1 = sum1 + term
        else
          sum1 = sum1 - term
        end if

        if ( abs ( term ) .lt. 1.0D-10 ) then
          go to 10
        else if ( abs ( term ) .lt. 1.0D-08 * abs ( sum1 ) ) then
          go to 10
        end if

      end do

10    continue

      euler_number2 = 2.0D+00**( n + 2 ) * sum1 * r8_factorial ( n ) 
     &  / pi**( n + 1 )

      if ( mod ( n, 4 ) .ne. 0 ) then
        euler_number2 = - euler_number2
      end if

      return
      end
      subroutine euler_number_values ( n_data, n, c )

c*********************************************************************72
c
cc EULER_NUMBER_VALUES returns some values of the Euler numbers.
c
c  Discussion:
c
c    In Mathematica, the function can be evaluated by:
c
c      EulerE[n]
c
c    These numbers rapidly get too big to store in an ordinary integer.
c
c    The terms of odd index are 0.
c
c    E(N) = -C(N,N-2) * E(N-2) - C(N,N-4) * E(N-4) - ... - C(N,0) * E(0).
c
c  First terms:
c
c    E0  = 1
c    E1  = 0
c    E2  = -1
c    E3  = 0
c    E4  = 5
c    E5  = 0
c    E6  = -61
c    E7  = 0
c    E8  = 1385
c    E9  = 0
c    E10 = -50521
c    E11 = 0
c    E12 = 2702765
c    E13 = 0
c    E14 = -199360981
c    E15 = 0
c    E16 = 19391512145
c    E17 = 0
c    E18 = -2404879675441
c    E19 = 0
c    E20 = 370371188237525
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    22 March 2007
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
c    Output, integer N, the order of the Euler number.
c
c    Output, integer C, the value of the Euler number.
c
      implicit none

      integer n_max
      parameter ( n_max = 8 )

      integer c
      integer c_vec(n_max)
      integer n
      integer n_data
      integer n_vec(n_max)

      save c_vec
      save n_vec

      data c_vec /
     &  1, 0, -1, 5, 61, 1385, -50521, 2702765 /
      data n_vec /
     &   0, 1, 2, 4, 6, 8, 10, 12 /

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
      function euler_poly ( n, x )

c*********************************************************************72
c
cc EULER_POLY evaluates the N-th Euler polynomial at X.
c
c  First values:
c
c    E(0,X) = 1
c    E(1,X) = X - 1/2
c    E(2,X) = X^2 - X 
c    E(3,X) = X^3 - 3/2 X^2 + 1/4
c    E(4,X) = X^4 - 2*X^3 + X
c    E(5,X) = X^5 - 5/2 X^4 + 5/2 X^2 - 1/2
c    E(6,X) = X^6 - 3 X^5 + 5 X^3 - 3 X
c    E(7,X) = X^7 - 7/2 X^6 + 35/4 X^4 - 21/2 X^2 + 17/8
c    E(8,X) = X^8 - 4 X^7 + 14 X^5 - 28 X^3 + 17 X
c
c  Special values:
c
c    E'(N,X) = N * E(N-1,X)
c
c    E(N,1/2) = E(N) / 2**N, where E(N) is the N-th Euler number.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    08 July 2008
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the order of the Euler polynomial to
c    be evaluated.  N must be 0 or greater.
c
c    Input, double precision X, the value at which the polynomial is to
c    be evaluated.
c
c    Output, double precision EULER_POLY, the value of E(N,X).
c
      implicit none

      double precision bx1
      double precision bx2
      double precision euler_poly
      integer n
      double precision x

      call bernoulli_poly2 ( n+1, x, bx1 )
      call bernoulli_poly2 ( n+1, 0.5D+00 * x, bx2 )

      euler_poly = 2.0D+00 * ( bx1 - bx2 * 2.0D+00**( n + 1 ) ) 
     &  / dble ( n + 1 )

      return
      end
      subroutine eulerian ( n, e )

c*********************************************************************72
c
cc EULERIAN computes the Eulerian number E(N,K).
c
c  Discussion:
c
c    A run in a permutation is a sequence of consecutive ascending values.
c
c    E(N,K) is the number of permutations of N objects which contain
c    exactly K runs.
c
c  Examples:
c
c     N = 7
c
c     1     0     0     0     0     0     0
c     1     1     0     0     0     0     0
c     1     4     1     0     0     0     0
c     1    11    11     1     0     0     0
c     1    26    66    26     1     0     0
c     1    57   302   302    57     1     0
c     1   120  1191  2416  1191   120     1
c
c  Recursion:
c
c    E(N,K) = K * E(N-1,K) + (N-K+1) * E(N-1,K-1).
c
c  Properties:
c
c    E(N,1) = E(N,N) = 1.
c    E(N,K) = 0 if K <= 0 or N < K.
c    sum ( 1 <= K <= N ) E(N,K) = N!.
c    X^N = sum ( 0 <= K <= N ) COMB(X+K-1, N ) E(N,K)
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    08 July 2008
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Dennis Stanton, Dennis White,
c    Constructive Combinatorics,
c    Springer Verlag, 1986
c
c  Parameters:
c
c    Input, integer N, the number of rows desired.
c
c    Output, integer E(N,N), the first N rows of Eulerian numbers.
c
      implicit none

      integer n

      integer e(n,n)
      integer i
      integer j

      if ( n .lt. 1 ) then
        return
      end if
!
!  Construct rows 1, 2, ..., N of the Eulerian triangle.
!
      e(1,1) = 1
      do j = 2, n
        e(1,j) = 0
      end do

      do i = 2, n
        e(i,1) = 1
        do j = 2, n
          e(i,j) = j * e(i-1,j) + ( i - j + 1 ) * e(i-1,j-1)
        end do
      end do

      return
      end
      subroutine fibonacci_direct ( n, f )

c*********************************************************************72
c
cc FIBONACCI_DIRECT computes the N-th Fibonacci number directly.
c
c  Discussion:
c
c    A direct formula for the N-th Fibonacci number is:
c
c      F(N) = ( PHIP**N - PHIM**N ) / sqrt(5)
c
c    where 
c
c      PHIP = ( 1 + sqrt(5) ) / 2, 
c      PHIM = ( 1 - sqrt(5) ) / 2.
c
c  Example:
c
c     N   F
c    --  --
c     0   0
c     1   1
c     2   1
c     3   2
c     4   3
c     5   5
c     6   8
c     7  13
c     8  21
c     9  34
c    10  55
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    13 July 2008
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the index of the Fibonacci number 
c    to compute.  N should be nonnegative.
c
c    Output, integer F, the value of the N-th Fibonacci number.
c
      implicit none

      integer f
      integer n
      double precision sqrt5
      parameter ( sqrt5 = 2.236068D+00 )
      double precision phim
      parameter ( phim = ( 1.0D+00 - sqrt5 ) / 2.0D+00 )
      double precision phip
      parameter ( phip = ( 1.0D+00 + sqrt5 ) / 2.0D+00 )

      if ( n .lt. 0 ) then
       f = 0
      else
        f = nint ( ( phip**n - phim**n ) / sqrt ( 5.0D+00 ) )
      end if
 
      return
      end
      subroutine fibonacci_floor ( n, f, i )

c*********************************************************************72
c
cc FIBONACCI_FLOOR returns the largest Fibonacci number less than or equal to N.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    13 July 2008
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the positive integer whose Fibonacci 
c    "floor" is desired.
c
c    Output, integer F, the largest Fibonacci number less 
c    than or equal to N.
c
c    Output, integer ( kind = 4 ) I, the index of the F.
c
      implicit none

      integer f
      integer i
      integer n

      if ( n .le. 0 ) then

        i = 0
        f = 0

      else

        i = int ( 
     &      log ( 0.5D+00 * dble ( 2 * n + 1 ) * sqrt ( 5.0D+00 ) ) 
     &    / log ( 0.5D+00 * ( 1.0D+00 + sqrt ( 5.0D+00 ) ) ) )

        call fibonacci_direct ( i, f )

        if ( n .lt. f ) then
          i = i - 1
          call fibonacci_direct ( i, f )
        end if

      end if

      return
      end
      subroutine fibonacci_recursive ( n, f )

c*********************************************************************72
c
cc FIBONACCI_RECURSIVE computes the first N Fibonacci numbers.
c
c  Discussion:
c
c    The 'golden ratio' 
c
c      PHI = (1+sqrt(5))/2 
c
c    satisfies the algebraic equation:
c
c      X*X-X-1=0
c
c    which is often written as:
c
c       X        1
c      --- =  ------
c       1      X - 1
c
c    expressing the fact that a rectangle, whose sides are in proportion X:1,
c    is similar to the rotated rectangle after a square of side 1 is removed.
c
c      <----X---->
c 
c      +-----*---*
c      |     |   |  1
c      |     |   | 
c      +-----*---+
c      <--1-><X-1>
c
c    A direct formula for the N-th Fibonacci number can be found.
c
c    Let
c
c      PHIP = ( 1 + sqrt(5) ) / 2
c      PHIM = ( 1 - sqrt(5) ) / 2
c
c    Then
c
c      F(N) = ( PHIP**N + PHIM**N ) / sqrt(5)
c
c    Moreover, F(N) can be computed by computing PHIP**N / sqrt(5) and rounding
c    to the nearest whole number.
c
c    The function 
c
c      F(X) = X / ( 1 - X - X^2 )
c
c    has a power series whose coefficients are the Fibonacci numbers:
c
c      F(X) = 0 + 1*X + 1*X^2 + 2*X^3 + 3*X^4 + 5*X^5+...
c
c  First terms:
c
c      0
c      1
c      1
c      2
c      3
c      5
c      8
c     13
c     21
c     34
c     55
c     89
c    144
c
c    The 40th number is                  102,334,155.
c    The 50th number is               12,586,269,025.
c    The 100th number is 354,224,848,179,261,915,075.
c
c  Recursion:
c
c    F(0) = 0
c    F(1) = 1
c
c    F(N) = F(N-1) + F(N-2)
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    13 July 2008
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the highest Fibonacci number to compute.
c
c    Output, integer F(N), the first N Fibonacci numbers.
c
      implicit none

      integer n

      integer f(n)
      integer i

      if ( n .le. 0 ) then
        return
      end if

      f(1) = 1

      if ( n .le. 1 ) then
        return
      end if

      f(2) = 1

      do i = 3, n
        f(i) = f(i-1) + f(i-2)
      end do
 
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
c    03 January 2006
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
     &  0.1524063822430784D+01,
     &  0.7966778177017837D+00,
     &  0.3982338580692348D+00,
     &  0.1520596783998375D+00,
     &  0.0000000000000000D+00,
     & -0.4987244125983972D-01,
     & -0.8537409000331584D-01,
     & -0.1081748095078604D+00,
     & -0.1196129141723712D+00,
     & -0.1207822376352452D+00,
     & -0.1125917656967557D+00,
     & -0.9580769740706586D-01,
     & -0.7108387291437216D-01,
     & -0.3898427592308333D-01,
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
      subroutine gamma_values ( n_data, x, fx )

c*********************************************************************72
c
cc GAMMA_VALUES returns some values of the Gamma function.
c
c  Discussion:
c
c    The Gamma function is defined as:
c
c      Gamma(Z) = Integral ( 0 <= T .lt. +oo) T**(Z-1) exp(-T) dT
c
c    It satisfies the recursion:
c
c      Gamma(X+1) = X * Gamma(X)
c
c    Gamma is undefined for nonpositive integral X.
c    Gamma(0.5) = sqrt(PI)
c    For N a positive integer, Gamma(N+1) = the standard factorial.
c
c    In Mathematica, the function can be evaluated by:
c
c      Gamma[x]
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    18 January 2008
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
c    Output, double precision X, the argument of the function.
c
c    Output, double precision FX, the value of the function.
c
      implicit none

      integer n_max
      parameter ( n_max = 25 )

      double precision fx
      double precision fx_vec(n_max)
      integer n_data
      double precision x
      double precision x_vec(n_max)

      save fx_vec
      save x_vec

      data fx_vec /
     &  -0.3544907701811032D+01,
     &  -0.1005871979644108D+03,
     &   0.9943258511915060D+02,
     &   0.9513507698668732D+01,
     &   0.4590843711998803D+01,
     &   0.2218159543757688D+01,
     &   0.1772453850905516D+01,
     &   0.1489192248812817D+01,
     &   0.1164229713725303D+01,
     &   0.1000000000000000D+01,
     &   0.9513507698668732D+00,
     &   0.9181687423997606D+00,
     &   0.8974706963062772D+00,
     &   0.8872638175030753D+00,
     &   0.8862269254527580D+00,
     &   0.8935153492876903D+00,
     &   0.9086387328532904D+00,
     &   0.9313837709802427D+00,
     &   0.9617658319073874D+00,
     &   0.1000000000000000D+01,
     &   0.2000000000000000D+01,
     &   0.6000000000000000D+01,
     &   0.3628800000000000D+06,
     &   0.1216451004088320D+18,
     &   0.8841761993739702D+31 /
      data x_vec /
     &  -0.50D+00,
     &  -0.01D+00,
     &   0.01D+00,
     &   0.10D+00,
     &   0.20D+00,
     &   0.40D+00,
     &   0.50D+00,
     &   0.60D+00,
     &   0.80D+00,
     &   1.00D+00,
     &   1.10D+00,
     &   1.20D+00,
     &   1.30D+00,
     &   1.40D+00,
     &   1.50D+00,
     &   1.60D+00,
     &   1.70D+00,
     &   1.80D+00,
     &   1.90D+00,
     &   2.00D+00,
     &   3.00D+00,
     &   4.00D+00,
     &  10.00D+00,
     &  20.00D+00,
     &  30.00D+00 /

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
      subroutine gegenbauer_poly ( n, alpha, x, cx )

c*********************************************************************72
c
cc GEGENBAUER_POLY computes the Gegenbauer polynomials C(I,ALPHA,X).
c
c  Discussion:
c
c    The Gegenbauer polynomial can be evaluated in Mathematica with
c    the command 
c
c      GegenbauerC[n,m,x]
c
c  Differential equation:
c
c    (1-X*X) Y'' - (2 ALPHA + 1) X Y' + N (N + 2 ALPHA) Y = 0
c
c  Recursion:
c
c    C(0,ALPHA,X) = 1,
c    C(1,ALPHA,X) = 2*ALPHA*X
c    C(N,ALPHA,X) = ( (2*N-2+2*ALPHA) * X * C(N-1,ALPHA,X) 
c                   + ( -N+2-2*ALPHA)     * C(N-2,ALPHA,X) ) / N
c
c  Restrictions:
c
c    ALPHA must be greater than -0.5.
c
c  Special values:
c
c    If ALPHA = 1, the Gegenbauer polynomials reduce to the Chebyshev
c    polynomials of the second kind.
c
c  Norm:
c
c    Integral ( -1 <= X <= 1 ) 
c      ( 1 - X^2 )**( ALPHA - 0.5 ) * C(N,ALPHA,X)^2 dX
c
c    = PI * 2**( 1 - 2 * ALPHA ) * Gamma ( N + 2 * ALPHA ) 
c      / ( N! * ( N + ALPHA ) * ( Gamma ( ALPHA ) )^2 )
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    09 July 2008
c
c  Author:
c
c    John Burkardt
c
c  Reference:
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
c    Input, integer N, the highest order polynomial to compute.
c    Note that polynomials 0 through N will be computed.
c
c    Input, double precision ALPHA, a parameter which is part of the 
c    definition of the Gegenbauer polynomials.  It must be greater than -0.5.
c
c    Input, double precision X, the point at which the polynomials 
c    are to be evaluated.
c
c    Output, double precision CX(0:N), the values of the first N+1 Gegenbauer
c    polynomials at the point X.  
c
      implicit none

      integer n

      double precision alpha
      double precision cx(0:n)
      integer i
      double precision x

      if ( alpha .le. -0.5D+00 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'GEGENBAUER_POLY - Fatal error!'
        write ( *, '(a,g14.6)' ) '  Illegal value of ALPHA = ', alpha
        write ( *, '(a)' ) '  but ALPHA must be greater than -0.5.'
        return
      end if

      if ( n .lt. 0 ) then
        return
      end if

      cx(0) = 1.0D+00

      if ( n .eq. 0 ) then
        return
      end if

      cx(1) = 2.0D+00 * alpha * x

      do i = 2, n
        cx(i) = 
     &    ( ( dble ( 2 * i - 2 ) + 2.0D+00 * alpha ) * x * cx(i-1)   
     &    + ( dble (   - i + 2 ) - 2.0D+00 * alpha )     * cx(i-2) ) 
     &    /   dble (     i     )
      end do
 
      return
      end
      subroutine gegenbauer_poly_values ( n_data, n, a, x, fx )

c*********************************************************************72
c
cc GEGENBAUER_POLY_VALUES returns some values of the Gegenbauer polynomials.
c
c  Discussion:
c
c    The Gegenbauer polynomials are also known as the "spherical
c    polynomials" or "ultraspherical polynomials".
c
c    In Mathematica, the function can be evaluated by:
c
c      GegenbauerC[n,m,x]
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    24 March 2007
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
c    Output, integer N, the order parameter of the function.
c
c    Output, double precision A, the real parameter of the function.
c
c    Output, double precision X, the argument of the function.
c
c    Output, double precision FX, the value of the function.
c
      implicit none

      integer n_max
      parameter ( n_max = 38 )

      double precision a
      double precision a_vec(n_max)
      double precision fx
      double precision fx_vec(n_max)
      integer n
      integer n_data
      integer n_vec(n_max)
      double precision x
      double precision x_vec(n_max)

      save a_vec
      save fx_vec
      save n_vec
      save x_vec

      data a_vec /
     &   0.5D+00,
     &   0.5D+00,
     &   0.5D+00,
     &   0.5D+00,
     &   0.5D+00,
     &   0.5D+00,
     &   0.5D+00,
     &   0.5D+00,
     &   0.5D+00,
     &   0.5D+00,
     &   0.5D+00,
     &   0.0D+00,
     &   1.0D+00,
     &   2.0D+00,
     &   3.0D+00,
     &   4.0D+00,
     &   5.0D+00,
     &   6.0D+00,
     &   7.0D+00,
     &   8.0D+00,
     &   9.0D+00,
     &  10.0D+00,
     &   3.0D+00,
     &   3.0D+00,
     &   3.0D+00,
     &   3.0D+00,
     &   3.0D+00,
     &   3.0D+00,
     &   3.0D+00,
     &   3.0D+00,
     &   3.0D+00,
     &   3.0D+00,
     &   3.0D+00,
     &   3.0D+00,
     &   3.0D+00,
     &   3.0D+00,
     &   3.0D+00,
     &   3.0D+00 /
      data fx_vec /
     &    1.0000000000D+00,
     &    0.2000000000D+00,
     &   -0.4400000000D+00,
     &   -0.2800000000D+00,
     &    0.2320000000D+00,
     &    0.3075200000D+00,
     &   -0.0805760000D+00,
     &   -0.2935168000D+00,
     &   -0.0395648000D+00,
     &    0.2459712000D+00,
     &    0.1290720256D+00,
     &    0.0000000000D+00,
     &   -0.3600000000D+00,
     &   -0.0800000000D+00,
     &    0.8400000000D+00,
     &    2.4000000000D+00,
     &    4.6000000000D+00,
     &    7.4400000000D+00,
     &   10.9200000000D+00,
     &   15.0400000000D+00,
     &   19.8000000000D+00,
     &   25.2000000000D+00,
     &   -9.0000000000D+00,
     &   -0.1612800000D+00,
     &   -6.6729600000D+00,
     &   -8.3750400000D+00,
     &   -5.5267200000D+00,
     &    0.0000000000D+00,
     &    5.5267200000D+00,
     &    8.3750400000D+00,
     &    6.6729600000D+00,
     &    0.1612800000D+00,
     &   -9.0000000000D+00,
     &  -15.4252800000D+00,
     &   -9.6969600000D+00,
     &   22.4409600000D+00,
     &  100.8892800000D+00,
     &  252.0000000000D+00 /
      data n_vec /
     &   0,  1,  2,
     &   3,  4,  5,
     &   6,  7,  8,
     &   9, 10,  2,
     &   2,  2,  2,
     &   2,  2,  2,
     &   2,  2,  2,
     &   2,  5,  5,
     &   5,  5,  5,
     &   5,  5,  5,
     &   5,  5,  5,
     &   5,  5,  5,
     &   5,  5 /
      data x_vec /
     &   0.20D+00,
     &   0.20D+00,
     &   0.20D+00,
     &   0.20D+00,
     &   0.20D+00,
     &   0.20D+00,
     &   0.20D+00,
     &   0.20D+00,
     &   0.20D+00,
     &   0.20D+00,
     &   0.20D+00,
     &   0.40D+00,
     &   0.40D+00,
     &   0.40D+00,
     &   0.40D+00,
     &   0.40D+00,
     &   0.40D+00,
     &   0.40D+00,
     &   0.40D+00,
     &   0.40D+00,
     &   0.40D+00,
     &   0.40D+00,
     &  -0.50D+00,
     &  -0.40D+00,
     &  -0.30D+00,
     &  -0.20D+00,
     &  -0.10D+00,
     &   0.00D+00,
     &   0.10D+00,
     &   0.20D+00,
     &   0.30D+00,
     &   0.40D+00,
     &   0.50D+00,
     &   0.60D+00,
     &   0.70D+00,
     &   0.80D+00,
     &   0.90D+00,
     &   1.00D+00 /

      if ( n_data .lt. 0 ) then
        n_data = 0
      end if

      n_data = n_data + 1

      if ( n_max .lt. n_data ) then
        n_data = 0
        n = 0
        a = 0.0D+00
        x = 0.0D+00
        fx = 0.0D+00
      else
        n = n_vec(n_data)
        a = a_vec(n_data)
        x = x_vec(n_data)
        fx = fx_vec(n_data)
      end if

      return
      end
      subroutine gen_hermite_poly ( n, x, mu, p )

c*********************************************************************72
c
cc GEN_HERMITE_POLY evaluates the generalized Hermite polynomials at X.
c
c  Discussion:
c
c    The generalized Hermite polynomials are orthogonal under the weight
c    function:
c
c      w(x) = |x|^(2*MU) * exp ( - x^2 )
c
c    over the interval (-oo,+oo).
c
c    When MU = 0, the generalized Hermite polynomial reduces to the standard
c    Hermite polynomial.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    26 February 2010
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Theodore Chihara,
c    An Introduction to Orthogonal Polynomials,
c    Gordon and Breach, 1978,
c    ISBN: 0677041500,
c    LC: QA404.5 C44.
c
c  Parameters:
c
c    Input, integer N, the highest order polynomial to compute.
c
c    Input, double precision X, the point at which the polynomials are 
c    to be evaluated.
c
c    Input, double precision MU, the parameter.
c    - 1 / 2 < MU.
c
c    Output, double precision P(0:N), the values of the first N+1
c    polynomials at the point X.
c
      implicit none

      integer n

      integer i
      double precision mu
      double precision p(0:n)
      double precision theta
      double precision x

      if ( n .lt. 0 ) then
        return
      end if

      p(0) = 1.0D+00

      if ( n .eq. 0 ) then
        return
      end if

      p(1) = 2.0D+00 * x
     
      do i = 1, n - 1

        if ( mod ( i, 2 ) .eq. 0 ) then
          theta = 0.0D+00
        else
          theta = 2.0D+00 * mu
        end if

        p(i+1) = 2.0D+00 * x * p(i) 
     &    - 2.0D+00 * ( dble ( i ) + theta ) * p(i-1)

      end do
     
      return
      end
      subroutine gen_laguerre_poly ( n, alpha, x, cx )

c*********************************************************************72
c
cc GEN_LAGUERRE_POLY evaluates generalized Laguerre polynomials.
c
c  Differential equation:
c
c    X * Y'' + (ALPHA+1-X) * Y' + N * Y = 0
c
c  Recursion:
c
c    L(0,ALPHA,X) = 1
c    L(1,ALPHA,X) = 1+ALPHA-X
c
c    L(N,ALPHA,X) = ( (2*N-1+ALPHA-X) * L(N-1,ALPHA,X) 
c                   - (N-1+ALPHA) * L(N-2,ALPHA,X) ) / N
c
c  Restrictions:
c
c    -1 < ALPHA
c
c  Special values:
c
c    For ALPHA = 0, the generalized Laguerre polynomial L(N,ALPHA,X)
c    is equal to the Laguerre polynomial L(N,X).
c
c    For ALPHA integral, the generalized Laguerre polynomial
c    L(N,ALPHA,X) equals the associated Laguerre polynomial L(N,ALPHA,X).
c
c  Norm:
c
c    Integral ( 0 <= X < +oo ) exp ( - X ) * L(N,ALPHA,X)^2 dX
c    = Gamma ( N + ALPHA + 1 ) / Nc
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    28 February 2010
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
c  Parameters:
c
c    Input, integer N, the highest order function to compute.
c
c    Input, double precision ALPHA, the parameter.  -1 < ALPHA is required.
c
c    Input, double precision X, the point at which the functions are to be
c    evaluated.
c
c    Output, double precision CX(0:N), the polynomials of 
c    degrees 0 through N evaluated at the point X.
c
      implicit none

      integer n

      double precision alpha
      double precision cx(0:n)
      integer i
      double precision x

      if ( alpha .le. -1.0D+00 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'GEN_LAGUERRE_POLY - Fatal error!'
        write ( *, '(a,g14.6)' ) 
     &    '  The input value of ALPHA is ', alpha
        write ( *, '(a)' ) '  but ALPHA must be greater than -1.'
        stop
      end if
     
      if ( n .lt. 0 ) then
        return
      end if

      cx(0) = 1.0D+00

      if ( n .eq. 0 ) then
        return
      end if

      cx(1) = 1.0D+00 + alpha - x

      do i = 2, n
        cx(i) = ( ( dble ( 2 * i - 1 ) + alpha - x ) * cx(i-1)   
     &          + ( dble (   - i + 1 ) - alpha     ) * cx(i-2) ) 
     &            / dble (     i     )
      end do

      return
      end
      function gud ( x )

c*********************************************************************72
c
cc GUD evaluates the Gudermannian function.
c
c  Discussion:
c
c    The Gudermannian function relates the hyperbolic and trigonometric
c    functions.  For any argument X, there is a corresponding value
c    GAMMA so that
c
c      sinh(x) = tan(gamma).
c
c    The value GAMMA is called the Gudermannian of X.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    02 March 1999
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, double precision X, the argument of the Gudermannian.
c
c    Output, double precision GUD, the value of the Gudermannian.
c
      implicit none

      double precision gud
      double precision x

      gud = 2.0D+00 * atan ( tanh ( 0.5D+00 * x ) )

      return
      end
      subroutine gud_values ( n_data, x, fx )

c*********************************************************************72
c
cc GUD_VALUES returns some values of the Gudermannian function.
c
c  Discussion:
c
c    The Gudermannian function relates the hyperbolic and trigonomentric
c    functions.  For any argument X, there is a corresponding value
c    GD so that
c
c      SINH(X) = TAN(GD).
c
c    This value GD is called the Gudermannian of X and symbolized
c    GD(X).  The inverse Gudermannian function is given as input a value
c    GD and computes the corresponding value X.
c
c    GD(X) = 2 * arctan ( exp ( X ) ) - PI / 2
c
c    In Mathematica, the function can be evaluated by:
c
c      2 * Atan[Exp[x]] - Pi/2
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    24 March 2007
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Stephen Wolfram,
c    The Mathematica Book,
c    Fourth Edition,
c    Cambridge University Press, 1999,
c    ISBN: 0-521-64314-7,
c    LC: QA76.95.W65.
c
c    Daniel Zwillinger, editor,
c    CRC Standard Mathematical Tables and Formulae,
c    30th Edition,
c    CRC Press, 1996,
c    ISBN: 0-8493-2479-3,
c    LC: QA47.M315.
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
      parameter ( n_max = 13 )

      double precision fx
      double precision fx_vec(n_max)
      integer n_data
      double precision x
      double precision x_vec(n_max)

      save fx_vec
      save x_vec

      data fx_vec /
     &  -0.1301760336046015D+01,
     &  -0.8657694832396586D+00,
     &   0.0000000000000000D+00,
     &   0.9983374879348662D-01,
     &   0.1986798470079397D+00,
     &   0.4803810791337294D+00,
     &   0.8657694832396586D+00,
     &   0.1131728345250509D+01,
     &   0.1301760336046015D+01,
     &   0.1406993568936154D+01,
     &   0.1471304341117193D+01,
     &   0.1510419907545700D+01,
     &   0.1534169144334733D+01 /
      data x_vec /
     &  -2.0D+00,
     &  -1.0D+00,
     &   0.0D+00,
     &   0.1D+00,
     &   0.2D+00,
     &   0.5D+00,
     &   1.0D+00,
     &   1.5D+00,
     &   2.0D+00,
     &   2.5D+00,
     &   3.0D+00,
     &   3.5D+00,
     &   4.0D+00 /

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
      subroutine hermite_poly ( n, x, cx )

c*********************************************************************72
c
cc HERMITE_POLY evaluates the Hermite polynomials at X.
c
c  Differential equation:
c
c    Y'' - 2 X Y' + 2 N Y = 0
c
c  First terms:
c
c      1
c      2 X
c      4 X^2     -  2
c      8 X^3     - 12 X
c     16 X^4     - 48 X^2     + 12
c     32 X^5    - 160 X^3    + 120 X
c     64 X^6    - 480 X^4    + 720 X^2    - 120
c    128 X^7   - 1344 X^5   + 3360 X^3   - 1680 X
c    256 X^8   - 3584 X^6  + 13440 X^4  - 13440 X^2   + 1680
c    512 X^9   - 9216 X^7  + 48384 X^5  - 80640 X^3  + 30240 X
c   1024 X^10 - 23040 X^8 + 161280 X^6 - 403200 X^4 + 302400 X^2 - 30240
c
c  Recursion:
c
c    H(0,X) = 1,
c    H(1,X) = 2*X,
c    H(N,X) = 2*X * H(N-1,X) - 2*(N-1) * H(N-2,X)
c
c  Norm:
c
c    Integral ( -oo < X < oo ) exp ( - X^2 ) * H(N,X)^2 dX
c    = sqrt ( PI ) * 2^N * N!
c
c    H(N,X) = (-1)^N * exp ( X^2 ) * dn/dXn ( exp(-X^2 ) )
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    10 July 2008
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
c    Larry Andrews,
c    Special Functions of Mathematics for Engineers,
c    Second Edition, 
c    Oxford University Press, 1998.
c
c  Parameters:
c
c    Input, integer N, the highest order polynomial to compute.
c    Note that polynomials 0 through N will be computed.
c
c    Input, double precision X, the point at which the polynomials are 
c    to be evaluated.
c
c    Output, double precision CX(0:N), the values of the first N+1 Hermite
c    polynomials at the point X.
c
      implicit none

      integer n

      double precision cx(0:n)
      integer i
      double precision x

      if ( n .lt. 0 ) then
        return
      end if

      cx(0) = 1.0D+00

      if ( n .eq. 0 ) then
        return
      end if

      cx(1) = 2.0D+00 * x
 
      do i = 2, n
        cx(i) = 2.0D+00 * x * cx(i-1) 
     &    - 2.0D+00 * dble ( i - 1 ) * cx(i-2)
      end do
 
      return
      end
      subroutine hermite_poly_coef ( n, c )

c*********************************************************************72
c
cc HERMITE_POLY_COEF evaluates the Hermite polynomial coefficients.
c
c  First terms:
c
c    N/K     0     1      2      3       4     5      6    7      8    9   10
c
c     0      1
c     1      0     2
c     2     -2     0      4
c     3      0   -12      0      8
c     4     12     0    -48      0      16
c     5      0   120      0   -160       0    32
c     6   -120     0    720      0    -480     0     64
c     7      0 -1680      0   3360       0 -1344      0   128
c     8   1680     0 -13440      0   13440     0  -3584     0    256
c     9      0 30240      0 -80640       0 48384      0 -9216      0 512
c    10 -30240     0 302400      0 -403200     0 161280     0 -23040   0 1024 
c
c  Recursion:
c
c    H(0,X) = 1,
c    H(1,X) = 2*X,
c    H(N,X) = 2*X * H(N-1,X) - 2*(N-1) * H(N-2,X)
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    10 July 2008
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
c  Parameters:
c
c    Input, integer N, the highest order polynomial to compute.
c    Note that polynomials 0 through N will be computed.
c
c    Output, double precision C(0:N,0:N), the coefficients of the Hermite
c    polynomials.
c
      implicit none

      integer n

      double precision c(0:n,0:n)
      integer i
      integer j

      if ( n .lt. 0 ) then
        return
      end if

      do j = 1, n
        do i = 1, n
          c(i,j) = 0.0D+00
        end do
      end do

      c(0,0) = 1.0D+00

      if ( n == 0 ) then
        return
      end if

      c(1,1) = 2.0D+00
 
      do i = 2, n
        c(i,0)     =  -2.0D+00 * dble ( i - 1 ) * c(i-2,0)
        do j = 1, i - 2
          c(i,j) =   2.0D+00                  * c(i-1,j-1)  
     &              -2.0D+00 * dble ( i - 1 ) * c(i-2,j)
        end do
        c(i,  i-1) =   2.0D+00                  * c(i-1,  i-2)
        c(i,  i  ) =   2.0D+00                  * c(i-1,  i-1)
      end do
 
      return
      end
      subroutine hermite_poly_values ( n_data, n, x, fx )

c*********************************************************************72
c
cc HERMITE_POLY_VALUES returns some values of the Hermite polynomial.
c
c  Discussion:
c
c    In Mathematica, the function can be evaluated by:
c
c      HermiteH[n,x]
c
c  Differential equation:
c
c    Y'' - 2 X Y' + 2 N Y = 0
c
c  First terms:
c
c      1
c      2 X
c      4 X^2     -  2
c      8 X^3     - 12 X
c     16 X^4     - 48 X^2     + 12
c     32 X^5    - 160 X^3    + 120 X
c     64 X^6    - 480 X^4    + 720 X^2    - 120
c    128 X^7   - 1344 X^5   + 3360 X^3   - 1680 X
c    256 X^8   - 3584 X^6  + 13440 X^4  - 13440 X^2   + 1680
c    512 X^9   - 9216 X^7  + 48384 X^5  - 80640 X^3  + 30240 X
c   1024 X^10 - 23040 X^8 + 161280 X^6 - 403200 X^4 + 302400 X^2 - 30240
c
c  Recursion:
c
c    H(0,X) = 1,
c    H(1,X) = 2*X,
c    H(N,X) = 2*X * H(N-1,X) - 2*(N-1) * H(N-2,X)
c
c  Norm:
c
c    Integral ( -oo .lt. X .lt. +oo ) exp ( - X^2 ) * H(N,X)^2 dX
c    = sqrt ( PI ) * 2**N * Nc
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    24 March 2007
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
c    Output, integer N, the order of the polynomial.
c
c    Output, double precision X, the point where the polynomial is evaluated.
c
c    Output, double precision FX, the value of the function.
c
      implicit none

      integer n_max
      parameter ( n_max = 17 )

      double precision fx
      double precision fx_vec(n_max)
      integer n
      integer n_data
      integer n_vec(n_max)
      double precision x
      double precision x_vec(n_max)

      save fx_vec
      save n_vec
      save x_vec

      data fx_vec /
     &   0.1000000000000000D+01,
     &   0.1000000000000000D+02,
     &   0.9800000000000000D+02,
     &   0.9400000000000000D+03,
     &   0.8812000000000000D+04,
     &   0.8060000000000000D+05,
     &   0.7178800000000000D+06,
     &   0.6211600000000000D+07,
     &   0.5206568000000000D+08,
     &   0.4212712000000000D+09,
     &   0.3275529760000000D+10,
     &   0.2432987360000000D+11,
     &   0.1712370812800000D+12,
     &   0.4100000000000000D+02,
     &  -0.8000000000000000D+01,
     &   0.3816000000000000D+04,
     &   0.3041200000000000D+07 /
      data n_vec /
     &   0,  1,  2,
     &   3,  4,  5,
     &   6,  7,  8,
     &   9, 10, 11,
     &  12,  5,  5,
     &   5,  5 /
      data x_vec /
     &  5.0D+00,
     &  5.0D+00,
     &  5.0D+00,
     &  5.0D+00,
     &  5.0D+00,
     &  5.0D+00,
     &  5.0D+00,
     &  5.0D+00,
     &  5.0D+00,
     &  5.0D+00,
     &  5.0D+00,
     &  5.0D+00,
     &  5.0D+00,
     &  0.5D+00,
     &  1.0D+00,
     &  3.0D+00,
     &  1.0D+01 /

      if ( n_data .lt. 0 ) then
        n_data = 0
      end if

      n_data = n_data + 1

      if ( n_max .lt. n_data ) then
        n_data = 0
        n = 0
        x = 0.0D+00
        fx = 0.0D+00
      else
        n = n_vec(n_data)
        x = x_vec(n_data)
        fx = fx_vec(n_data)
      end if

      return
      end
      subroutine hyper_2f1_values ( n_data, a, b, c, x, fx )

c*********************************************************************72
c
cc HYPER_2F1_VALUES returns some values of the hypergeometric function 2F1.
c
c  Discussion:
c
c    In Mathematica, the function can be evaluated by:
c
c      fx = Hypergeometric2F1 [ a, b, c, x ]
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    09 September 2007
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
c    Cambridge University Press, 1999,
c    ISBN: 0-521-64314-7,
c    LC: QA76.95.W65.
c
c    Shanjie Zhang, Jianming Jin,
c    Computation of Special Functions,
c    Wiley, 1996,
c    ISBN: 0-471-11963-6,
c    LC: QA351.C45
c
c    Daniel Zwillinger, editor,
c    CRC Standard Mathematical Tables and Formulae,
c    30th Edition,
c    CRC Press, 1996,
c    ISBN: 0-8493-2479-3,
c    LC: QA47.M315.
c
c  Parameters:
c
c    Input/output, integer N_DATA.  The user sets N_DATA to 0 
c    before the first call.  On each call, the routine increments N_DATA by 1,
c    and returns the corresponding data; when there is no more data, the
c    output value of N_DATA will be 0 again.
c
c    Output, double precision A, B, C, X, the parameters.
c
c    Output, double precision FX, the value of the function.
c
      implicit none

      integer n_max
      parameter ( n_max = 24 )

      double precision a
      double precision a_vec(n_max)
      double precision b
      double precision b_vec(n_max)
      double precision c
      double precision c_vec(n_max)
      double precision fx
      double precision fx_vec(n_max)
      integer n_data
      double precision x
      double precision x_vec(n_max)

      save a_vec
      save b_vec
      save c_vec
      save fx_vec
      save x_vec

      data a_vec /
     &   -2.5D+00, 
     &   -0.5D+00, 
     &    0.5D+00, 
     &    2.5D+00, 
     &   -2.5D+00, 
     &   -0.5D+00, 
     &    0.5D+00, 
     &    2.5D+00, 
     &   -2.5D+00, 
     &   -0.5D+00, 
     &    0.5D+00, 
     &    2.5D+00, 
     &    3.3D+00, 
     &    1.1D+00, 
     &    1.1D+00, 
     &    3.3D+00, 
     &    3.3D+00, 
     &    1.1D+00, 
     &    1.1D+00, 
     &    3.3D+00, 
     &    3.3D+00, 
     &    1.1D+00, 
     &    1.1D+00, 
     &    3.3D+00 /
      data b_vec /
     &    3.3D+00, 
     &    1.1D+00, 
     &    1.1D+00, 
     &    3.3D+00, 
     &    3.3D+00, 
     &    1.1D+00, 
     &    1.1D+00, 
     &    3.3D+00, 
     &    3.3D+00, 
     &    1.1D+00, 
     &    1.1D+00, 
     &    3.3D+00, 
     &    6.7D+00, 
     &    6.7D+00, 
     &    6.7D+00, 
     &    6.7D+00, 
     &    6.7D+00, 
     &    6.7D+00, 
     &    6.7D+00, 
     &    6.7D+00, 
     &    6.7D+00, 
     &    6.7D+00, 
     &    6.7D+00, 
     &    6.7D+00 /
      data c_vec /
     &    6.7D+00, 
     &    6.7D+00, 
     &    6.7D+00, 
     &    6.7D+00, 
     &    6.7D+00, 
     &    6.7D+00, 
     &    6.7D+00, 
     &    6.7D+00, 
     &    6.7D+00, 
     &    6.7D+00, 
     &    6.7D+00, 
     &    6.7D+00, 
     &   -5.5D+00, 
     &   -0.5D+00, 
     &    0.5D+00, 
     &    4.5D+00, 
     &   -5.5D+00, 
     &   -0.5D+00, 
     &    0.5D+00, 
     &    4.5D+00, 
     &   -5.5D+00, 
     &   -0.5D+00, 
     &    0.5D+00, 
     &    4.5D+00 /
      data fx_vec /
     &    0.72356129348997784913D+00, 
     &    0.97911109345277961340D+00, 
     &    1.0216578140088564160D+00, 
     &    1.4051563200112126405D+00, 
     &    0.46961431639821611095D+00, 
     &    0.95296194977446325454D+00, 
     &    1.0512814213947987916D+00, 
     &    2.3999062904777858999D+00, 
     &    0.29106095928414718320D+00, 
     &    0.92536967910373175753D+00, 
     &    1.0865504094806997287D+00, 
     &    5.7381565526189046578D+00, 
     &    15090.669748704606754D+00, 
     &   -104.31170067364349677D+00, 
     &    21.175050707768812938D+00, 
     &    4.1946915819031922850D+00, 
     &    1.0170777974048815592D+10, 
     &   -24708.635322489155868D+00, 
     &    1372.2304548384989560D+00, 
     &    58.092728706394652211D+00, 
     &    5.8682087615124176162D+18, 
     &   -4.4635010147295996680D+08, 
     &    5.3835057561295731310D+06, 
     &    20396.913776019659426D+00 /
      data x_vec /
     &    0.25D+00, 
     &    0.25D+00, 
     &    0.25D+00, 
     &    0.25D+00, 
     &    0.55D+00, 
     &    0.55D+00, 
     &    0.55D+00, 
     &    0.55D+00, 
     &    0.85D+00, 
     &    0.85D+00, 
     &    0.85D+00, 
     &    0.85D+00, 
     &    0.25D+00, 
     &    0.25D+00, 
     &    0.25D+00, 
     &    0.25D+00, 
     &    0.55D+00, 
     &    0.55D+00, 
     &    0.55D+00, 
     &    0.55D+00, 
     &    0.85D+00, 
     &    0.85D+00, 
     &    0.85D+00, 
     &    0.85D+00 /

      if ( n_data .lt. 0 ) then
        n_data = 0
      end if
 
      n_data = n_data + 1

      if ( n_max .lt. n_data ) then
        n_data = 0
        a = 0.0D+00
        b = 0.0D+00
        c = 0.0D+00
        x = 0.0D+00
        fx = 0.0D+00
      else
        a = a_vec(n_data)
        b = b_vec(n_data)
        c = c_vec(n_data)
        x = x_vec(n_data)
        fx = fx_vec(n_data)
      end if

      return
      end
      function i4_choose ( n, k )

c*********************************************************************72
c
cc I4_CHOOSE computes the binomial coefficient C(N,K).
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
c    02 June 2007
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
c    Volume 6, Number 4, April 1963, page 161.
c
c  Parameters:
c
c    Input, integer N, K, are the values of N and K.
c
c    Output, integer I4_CHOOSE, the number of combinations of N
c    things taken K at a time.
c
      implicit none

      integer i
      integer i4_choose
      integer k
      integer mn
      integer mx
      integer n
      integer value

      mn = min ( k, n - k )

      if ( mn .lt. 0 ) then

        value = 0

      else if ( mn .eq. 0 ) then

        value = 1

      else

        mx = max ( k, n - k )
        value = mx + 1

        do i = 2, mn
          value = ( value * ( mx + i ) ) / i
        end do

      end if

      i4_choose = value

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
      subroutine i4_factorial_values ( n_data, n, fn )

c*********************************************************************72
c
cc I4_FACTORIAL_VALUES returns values of the factorial function.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    24 March 2007
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
c    Output, integer N, the argument of the function.
c
c    Output, integer FN, the value of the function.
c
      implicit none

      integer n_max
      parameter ( n_max = 13 )

      integer fn_vec(n_max)
      integer fn
      integer n
      integer n_data
      integer n_vec(n_max)

      save fn_vec
      save n_vec

      data fn_vec /
     &          1,
     &          1,
     &          2,
     &          6,
     &         24,
     &        120,
     &        720,
     &       5040,
     &      40320,
     &     362880,
     &    3628800,
     &   39916800,
     &  479001600 /
      data n_vec /
     &   0,  1,  2,  3,
     &   4,  5,  6,  7,
     &   8,  9, 10, 11,
     &  12 /

      if ( n_data .lt. 0 ) then
        n_data = 0
      end if

      n_data = n_data + 1

      if ( n_max .lt. n_data ) then
        n_data = 0
        n = 0
        fn = 0
      else
        n = n_vec(n_data)
        fn = fn_vec(n_data)
      end if

      return
      end
      function i4_factorial2 ( n )

c*********************************************************************72
c
cc I4_FACTORIAL2 computes the double factorial function.
c
c  Discussion:
c
c    The formula is:
c
c      FACTORIAL2( N ) = Product ( N * (N-2) * (N-4) * ... * 2 )  (N even)
c                      = Product ( N * (N-2) * (N-4) * ... * 1 )  (N odd)
c
c  Example:
c
c     N    Factorial2(N)
c
c     0     1
c     1     1
c     2     2
c     3     3
c     4     8
c     5    15
c     6    48
c     7   105
c     8   384
c     9   945
c    10  3840
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    03 July 2008
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the argument of the double factorial 
c    function.  If N is less than 1, I4_FACTORIAL2 is returned as 1.
c
c    Output, integer I4_FACTORIAL2, the value of the function.
c
      implicit none

      integer i4_factorial2
      integer n
      integer n_copy

      if ( n .lt. 1 ) then
        i4_factorial2 = 1
        return
      end if

      n_copy = n
      i4_factorial2 = 1

10    continue

      if ( 1 .lt. n_copy ) then
        i4_factorial2 = i4_factorial2 * n_copy
        n_copy = n_copy - 2
        go to 10
      end if

      return
      end
      subroutine i4_factorial2_values ( n_data, n, fn )

c*********************************************************************72
c
cc I4_FACTORIAL2_VALUES returns values of the double factorial function.
c
c  Discussion:
c
c    FACTORIAL2( N ) = Product ( N * (N-2) * (N-4) * ... * 2 )  (N even)
c                    = Product ( N * (N-2) * (N-4) * ... * 1 )  (N odd)
c
c  Example:
c
c     N    Fctorial2(N)
c
c     0     1
c     1     1
c     2     2
c     3     3
c     4     8
c     5    15
c     6    48
c     7   105
c     8   384
c     9   945
c    10  3840
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    24 March 2007
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
c    Cambridge University Press, 1999,
c    ISBN: 0-521-64314-7,
c    LC: QA76.95.W65.
c
c    Daniel Zwillinger, editor,
c    CRC Standard Mathematical Tables and Formulae,
c    30th Edition,
c    CRC Press, 1996,
c    ISBN: 0-8493-2479-3,
c    LC: QA47.M315.
c
c  Parameters:
c
c    Input/output, integer N_DATA.  The user sets N_DATA to 0 before the
c    first call.  On each call, the routine increments N_DATA by 1, and
c    returns the corresponding data; when there is no more data, the
c    output value of N_DATA will be 0 again.
c
c    Output, integer N, the argument of the function.
c
c    Output, integer FN, the value of the function.
c
      implicit none

      integer n_max
      parameter ( n_max = 16 )

      integer fn_vec(n_max)
      integer fn
      integer n_data
      integer n
      integer n_vec(n_max)

      save fn_vec
      save n_vec

      data fn_vec /
     &        1,
     &        1,
     &        2,
     &        3,
     &        8,
     &       15,
     &       48,
     &      105,
     &      384,
     &      945,
     &     3840,
     &    10395,
     &    46080,
     &   135135,
     &   645120,
     &  2027025 /
      data n_vec /
     &   0,
     &   1,  2,  3,  4,  5,
     &   6,  7,  8,  9, 10,
     &  11, 12, 13, 14, 15 /

      if ( n_data .lt. 0 ) then
        n_data = 0
      end if

      n_data = n_data + 1

      if ( n_max .lt. n_data ) then
        n_data = 0
        n = 0
        fn = 0
      else
        n = n_vec(n_data)
        fn = fn_vec(n_data)
      end if

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
      function i4_is_prime ( n )

c*********************************************************************72
c
cc I4_IS_PRIME reports whether an I4 is prime.
c
c  Discussion:
c
c    A simple, unoptimized sieve of Erasthosthenes is used to
c    check whether N can be divided by any integer between 2
c    and SQRT(N).
c
c    Note that negative numbers, 0 and 1 are not considered prime.
c
c    An I4 is an integer value.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    26 October 2007
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the integer to be tested.
c
c    Output, logical I4_IS_PRIME, is TRUE if N is prime, and FALSE
c    otherwise.
c
      implicit none

      integer i
      logical i4_is_prime
      integer n
      integer nhi

      if ( n .le. 0 ) then
        i4_is_prime = .false.
        return
      end if

      if ( n .eq. 1 ) then
        i4_is_prime = .false.
        return
      end if

      if ( n .le. 3 ) then
        i4_is_prime = .true.
        return
      end if

      nhi = int ( sqrt ( real ( n ) ) )

      do i = 2, nhi
        if ( mod ( n, i ) .eq. 0 ) then
          i4_is_prime = .false.
          return
        end if
      end do

      i4_is_prime = .true.

      return
      end
      function i4_is_triangular ( i )

c*********************************************************************72
c
cc I4_IS_TRIANGULAR determines whether an integer is triangular.
c
c  Discussion:
c
c    The N-th triangular number is equal to the sum of the first
c    N integers.
c
c  First Values:
c
c    Index  Value
c     0      0
c     1      1
c     2      3
c     3      6
c     4     10
c     5     15
c     6     21
c     7     28
c     8     36
c     9     45
c    10     55
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    19 February 2003
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer I, the integer to be checked.
c
c    Output, logical I4_IS_TRIANGULAR, is TRUE if I is triangular.
c
      implicit none

      integer i
      logical i4_is_triangular
      integer j
      integer k

      if ( i .lt. 0 ) then

        i4_is_triangular = .false.

      else if ( i .eq. 0 ) then

        i4_is_triangular = .true.

      else

        call i4_to_triangle ( i, j, k )

        if ( j .eq. k ) then
          i4_is_triangular = .true.
        else
          i4_is_triangular = .false.
        end if

      end if

      return
      end
      subroutine i4_partition_distinct_count ( n, q )

c*********************************************************************72
c
cc I4_PARTITION_DISTINCT_COUNT returns any value of Q(N).
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
c    However, if we require that each member of the partition
c    be distinct, we are computing something symbolized by Q(N).
c    The number 5 has Q(N) = 3, because it has the following partitions 
c    into distinct parts:
c
c    5 = 5
c      = 4 + 1 
c      = 3 + 2 
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    06 July 2008
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
c  Parameters:
c
c    Input, integer N, the integer to be partitioned.
c
c    Output, integer Q, the number of partitions of the integer
c    into distinct parts.
c
      implicit none

      integer n

      integer c(0:n)
      integer i
      logical i4_is_triangular
      integer k
      integer k2
      integer k_sign
      integer q

      c(0) = 1

      do i = 1, n

        if ( i4_is_triangular ( i ) ) then
          c(i) = 1
        else
          c(i) = 0
        end if

        k = 0
        k_sign = -1

10      continue

          k = k + 1
          k_sign = - k_sign
          k2 = k * ( 3 * k + 1 )

          if ( i .lt. k2 ) then
            go to 20
          end if

          c(i) = c(i) + k_sign * c(i-k2)

        go to 10

20      continue

        k = 0
        k_sign = -1

30      continue

          k = k + 1
          k_sign = - k_sign
          k2 = k * ( 3 * k - 1 )

          if ( i .lt. k2 ) then
            go to 40
          end if

          c(i) = c(i) + k_sign * c(i-k2)

        go to 30

40      continue

      end do

      q = c(n)

      return
      end
      function i4_pochhammer ( i, j )

c*********************************************************************72
c
cc I4_POCHHAMMER returns the value of ( I * (I+1) * ... * (J-1) * J ).
c
c  Discussion:
c
c    Pochhammer's symbol (A)_N is the value
c
c      (A)_N = Gamma ( A + N ) / Gamma ( A )
c
c    or, for integer arguments,
c
c      (I)_N = I * ( I + 1 ) * ... * ( I + N - 1 )
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    03 July 2008
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer I, J, values that define the product.
c
c    Output, integer I4_POCHHAMMER, the value of the product.
c
      implicit none

      integer i
      integer i4_pochhammer
      integer j
      integer k

      i4_pochhammer = 1
      do k = i, j
        i4_pochhammer = i4_pochhammer * k
      end do

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
      subroutine i4_to_triangle ( k, i, j )

c*********************************************************************72
c
cc I4_TO_TRIANGLE converts an integer to triangular coordinates.
c
c  Discussion:
c
c    Triangular coordinates are handy when storing a naturally triangular
c    array (such as the lower half of a matrix) in a linear array.
c
c    Thus, for example, we might consider storing 
c
c    (1,1)
c    (2,1) (2,2)
c    (3,1) (3,2) (3,3)
c    (4,1) (4,2) (4,3) (4,4)
c
c    as the linear array
c
c    (1,1) (2,1) (2,2) (3,1) (3,2) (3,3) (4,1) (4,2) (4,3) (4,4)    
c
c    Here, the quantities in parenthesis represent the natural row and
c    column indices of a single number when stored in a rectangular array.
c
c    In this routine, we are given the location K of an item in the 
c    linear array, and wish to determine the row I and column J
c    of the item when stored in the triangular array.
c 
c  First Values:
c
c     K  I  J
c
c     0  0  0
c     1  1  1
c     2  2  1
c     3  2  2
c     4  3  1
c     5  3  2
c     6  3  3
c     7  4  1
c     8  4  2
c     9  4  3
c    10  4  4
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    03 July 2008
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer K, the linear index of the (I,J) element, 
c    which must be nonnegative.
c
c    Output, integer I, J, the row and column indices.
c
      implicit none

      integer i
      integer j
      integer k

      if ( k .lt. 0 ) then

        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'I4_TO_TRIANGLE - Fatal error!'
        write ( *, '(a)' ) '  K < 0.'
        write ( *, '(a,i8)' ) '  K = ', k
        stop
  
      else if ( k .eq. 0 ) then

        i = 0
        j = 0
        return

      end if

      i = int ( sqrt ( real ( 2 * k ) ) )

      if ( i * i + i .lt. 2 * k ) then
        i = i + 1
      end if

      j = k - ( i * ( i - 1 ) ) / 2

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
      subroutine jacobi_poly ( n, alpha, beta, x, cx )

c*********************************************************************72
c
cc JACOBI_POLY evaluates the Jacobi polynomials at X.
c
c  Differential equation:
c
c    (1-X*X) Y'' + (BETA-ALPHA-(ALPHA+BETA+2) X) Y' + N (N+ALPHA+BETA+1) Y = 0
c
c  Recursion:
c
c    P(0,ALPHA,BETA,X) = 1,
c
c    P(1,ALPHA,BETA,X) = ( (2+ALPHA+BETA)*X + (ALPHA-BETA) ) / 2
c
c    P(N,ALPHA,BETA,X)  = 
c      ( 
c        (2*N+ALPHA+BETA-1) 
c        * ((ALPHA^2-BETA^2)+(2*N+ALPHA+BETA)*(2*N+ALPHA+BETA-2)*X) 
c        * P(N-1,ALPHA,BETA,X)
c        -2*(N-1+ALPHA)*(N-1+BETA)*(2*N+ALPHA+BETA) * P(N-2,ALPHA,BETA,X)
c      ) / 2*N*(N+ALPHA+BETA)*(2*N-2+ALPHA+BETA)
c
c  Restrictions:
c
c    -1 < ALPHA
c    -1 < BETA
c
c  Norm:
c
c    Integral ( -1 <= X <= 1 ) ( 1 - X )^ALPHA * ( 1 + X )^BETA 
c      * P(N,ALPHA,BETA,X)^2 dX 
c    = 2^(ALPHA+BETA+1) * Gamma ( N + ALPHA + 1 ) * Gamma ( N + BETA + 1 ) /
c      ( 2 * N + ALPHA + BETA ) * N! * Gamma ( N + ALPHA + BETA + 1 )
c
c  Special values:
c
c    P(N,ALPHA,BETA,1) = (N+ALPHA)!/(N!*ALPHA!) for integer ALPHA.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    06 July 2008
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
c  Parameters:
c
c    Input, integer N, the highest order polynomial to compute.  
c    Note that polynomials 0 through N will be computed.
c
c    Input, double precision ALPHA, one of the parameters defining the Jacobi
c    polynomials, ALPHA must be greater than -1.
c
c    Input, double precision BETA, the second parameter defining the Jacobi
c    polynomials, BETA must be greater than -1.
c
c    Input, double precision X, the point at which the polynomials are 
c    to be evaluated.
c
c    Output, double precision CX(0:N), the values of the first N+1 Jacobi
c    polynomials at the point X.
c
      implicit none

      integer n

      double precision alpha
      double precision beta
      double precision cx(0:n)
      double precision c1
      double precision c2
      double precision c3
      double precision c4
      integer i
      double precision r_i
      double precision x

      if ( alpha .le. -1.0D+00 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'JACOBI_POLY - Fatal error!'
        write ( *, '(a,g14.6)' ) 
     &    '  Illegal input value of ALPHA = ', alpha
        write ( *, '(a)' ) '  But ALPHA must be greater than -1.'
        stop
      end if
 
      if ( beta .le. -1.0D+00 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'JACOBI_POLY - Fatal error!'
        write ( *, '(a,g14.6)' ) 
     &    '  Illegal input value of BETA = ', beta
        write ( *, '(a)' ) '  But BETA must be greater than -1.'
        stop
      end if
  
      if ( n .lt. 0 ) then
        return
      end if

      cx(0) = 1.0D+00

      if ( n .eq. 0 ) then
        return
      end if

      cx(1) = ( 1.0D+00 + 0.5D+00 * ( alpha + beta ) ) * x 
     &  + 0.5D+00 * ( alpha - beta )
 
      do i = 2, n

        r_i = dble ( i ) 

        c1 = 2.0D+00 * r_i * ( r_i + alpha + beta ) 
     &    * ( 2.0D+00 * r_i - 2.0D+00 + alpha + beta )

        c2 = ( 2.0D+00 * r_i - 1.0D+00 + alpha + beta ) 
     &    * ( 2.0D+00 * r_i  + alpha + beta ) 
     &    * ( 2.0D+00 * r_i - 2.0D+00 + alpha + beta )

        c3 = ( 2.0D+00 * r_i - 1.0D+00 + alpha + beta ) 
     &    * ( alpha + beta ) * ( alpha - beta )

        c4 = - 2.0D+00 * ( r_i - 1.0D+00 + alpha ) 
     &    * ( r_i - 1.0D+00 + beta ) 
     &    * ( 2.0D+00 * r_i + alpha + beta )

        cx(i) = ( ( c3 + c2 * x ) * cx(i-1) + c4 * cx(i-2) ) / c1

      end do

      return
      end
      subroutine jacobi_poly_values ( n_data, n, a, b, x, fx )

c*********************************************************************72
c
cc JACOBI_POLY_VALUES returns some values of the Jacobi polynomial.
c
c  Discussion:
c
c    In Mathematica, the function can be evaluated by:
c
c      JacobiP[ n, a, b, x ]
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    19 April 2012
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
c    Output, integer N, the degree of the polynomial.
c
c    Output, integer A, B, parameters of the function.
c
c    Output, double precision X, the argument of the function.
c
c    Output, double precision FX, the value of the function.
c
      implicit none

      integer n_max
      parameter ( n_max = 26 )

      double precision a
      double precision a_vec(n_max)
      double precision b
      double precision b_vec(n_max)
      double precision fx
      double precision fx_vec(n_max)
      integer n
      integer n_data
      integer n_vec(n_max)
      double precision x
      double precision x_vec(n_max)

      save a_vec
      save b_vec
      save fx_vec
      save n_vec
      save x_vec

      data a_vec /
     &   0.0D+00, 0.0D+00, 0.0D+00, 0.0D+00,
     &   0.0D+00, 0.0D+00, 1.0D+00, 2.0D+00,
     &   3.0D+00, 4.0D+00, 5.0D+00, 0.0D+00,
     &   0.0D+00, 0.0D+00, 0.0D+00, 0.0D+00,
     &   0.0D+00, 0.0D+00, 0.0D+00, 0.0D+00,
     &   0.0D+00, 0.0D+00, 0.0D+00, 0.0D+00,
     &   0.0D+00, 0.0D+00 /
      data b_vec /
     &   1.0D+00, 1.0D+00, 1.0D+00, 1.0D+00,
     &   1.0D+00, 1.0D+00, 1.0D+00, 1.0D+00,
     &   1.0D+00, 1.0D+00, 1.0D+00, 2.0D+00,
     &   3.0D+00, 4.0D+00, 5.0D+00, 1.0D+00,
     &   1.0D+00, 1.0D+00, 1.0D+00, 1.0D+00,
     &   1.0D+00, 1.0D+00, 1.0D+00, 1.0D+00,
     &   1.0D+00, 1.0D+00 /
      data fx_vec /
     &    0.1000000000000000D+01,
     &    0.2500000000000000D+00,
     &   -0.3750000000000000D+00,
     &   -0.4843750000000000D+00,
     &   -0.1328125000000000D+00,
     &    0.2753906250000000D+00,
     &   -0.1640625000000000D+00,
     &   -0.1174804687500000D+01,
     &   -0.2361328125000000D+01,
     &   -0.2616210937500000D+01,
     &    0.1171875000000000D+00,
     &    0.4218750000000000D+00,
     &    0.5048828125000000D+00,
     &    0.5097656250000000D+00,
     &    0.4306640625000000D+00,
     &   -0.6000000000000000D+01,
     &    0.3862000000000000D-01,
     &    0.8118400000000000D+00,
     &    0.3666000000000000D-01,
     &   -0.4851200000000000D+00,
     &   -0.3125000000000000D+00,
     &    0.1891200000000000D+00,
     &    0.4023400000000000D+00,
     &    0.1216000000000000D-01,
     &   -0.4396200000000000D+00,
     &    0.1000000000000000D+01 /
      data n_vec /
     &    0, 1, 2, 3,
     &    4, 5, 5, 5,
     &    5, 5, 5, 5,
     &    5, 5, 5, 5,
     &    5, 5, 5, 5,
     &    5, 5, 5, 5,
     &    5, 5 /
      data x_vec /
     &    0.5D+00,
     &    0.5D+00,
     &    0.5D+00,
     &    0.5D+00,
     &    0.5D+00,
     &    0.5D+00,
     &    0.5D+00,
     &    0.5D+00,
     &    0.5D+00,
     &    0.5D+00,
     &    0.5D+00,
     &    0.5D+00,
     &    0.5D+00,
     &    0.5D+00,
     &    0.5D+00,
     &   -1.0D+00,
     &   -0.8D+00,
     &   -0.6D+00,
     &   -0.4D+00,
     &   -0.2D+00,
     &    0.0D+00,
     &    0.2D+00,
     &    0.4D+00,
     &    0.6D+00,
     &    0.8D+00,
     &    1.0D+00 /

      if ( n_data .lt. 0 ) then
        n_data = 0
      end if

      n_data = n_data + 1

      if ( n_max .lt. n_data ) then
        n_data = 0
        n = 0
        a = 0.0D+00
        b = 0.0D+00
        x = 0.0D+00
        fx = 0.0D+00
      else
        n = n_vec(n_data)
        a = a_vec(n_data)
        b = b_vec(n_data)
        x = x_vec(n_data)
        fx = fx_vec(n_data)
      end if

      return
      end
      subroutine jacobi_symbol ( q, p, j )

c*********************************************************************72
c
cc JACOBI_SYMBOL evaluates the Jacobi symbol (Q/P).
c
c  Discussion:
c
c    If P is prime, then
c
c      Jacobi Symbol (Q/P) = Legendre Symbol (Q/P)
c
c    Else 
c
c      let P have the prime factorization
c
c        P = Product ( 1 <= I <= N ) P(I)^E(I)
c
c      Jacobi Symbol (Q/P) =
c
c        Product ( 1 <= I <= N ) Legendre Symbol (Q/P(I))^E(I)
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    11 July 2008
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Daniel Zwillinger,
c    CRC Standard Mathematical Tables and Formulae,
c    30th Edition,
c    CRC Press, 1996, pages 86-87.
c
c  Parameters:
c
c    Input, integer Q, an integer whose Jacobi symbol with
c    respect to P is desired.
c
c    Input, integer P, the number with respect to which the Jacobi
c    symbol of Q is desired.  P should be 2 or greater.
c
c    Output, integer J, the Jacobi symbol (Q/P).
c    Ordinarily, J will be -1, 0 or 1.
c    -2, not enough factorization space.
c    -3, an error during Legendre symbol calculation.
c    
      implicit none

      integer maxfactor
      parameter ( maxfactor = 20 )

      integer factor(maxfactor)
      integer i
      integer j
      integer l
      integer nfactor
      integer nleft
      integer p
      integer power(maxfactor)
      integer pp
      integer q
      integer qq
c
c  P must be greater than 1.
c
      if ( p .le. 1 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'JACOBI_SYMBOL - Fatal error!'
        write ( *, '(a)' ) '  P must be greater than 1.'
        l = -2
        return
      end if
c
c  Decompose P into factors of prime powers.
c
      call i4_factor ( p, maxfactor, nfactor, factor, power, nleft )

      if ( nleft .ne. 1 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'JACOBI_SYMBOL - Fatal error!'
        write ( *, '(a)' ) '  Not enough factorization space.'
        j = -2
        return
      end if
c
c  Force Q to be nonnegative.
c
      qq = q

10    continue

      if ( qq .lt. 0 ) then
        qq = qq + p
        go to 10
      end if
c
c  For each prime factor, compute the Legendre symbol, and
c  multiply the Jacobi symbol by the appropriate factor.
c
      j = 1
      do i = 1, nfactor
        pp = factor(i)
        call legendre_symbol ( qq, pp, l )
        if ( l .lt. -1 ) then
          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) 'JACOBI_SYMBOL - Fatal error!'
          write ( *, '(a)' ) 
     &      '  Error during Legendre symbol calculation.'
          j = -3
          return
        end if
        j = j * l**power(i)
      end do

      return
      end
      subroutine krawtchouk ( n, p, x, m, v )

c*********************************************************************72
c
cc KRAWTCHOUK evaluates the Krawtchouk polynomials at X.
c
c  Discussion:
c
c    The polynomial has a parameter P, which must be striclty between
c    0 and 1, and a parameter M which must be a nonnegative integer.
c
c    The Krawtchouk polynomial of order N, with parameters P and M,
c    evaluated at X, may be written K(N,P,X,M).
c
c    The first two terms are:
c
c      K(0,P,X,M) = 1
c      K(1,P,X,M) = X - P * M
c
c    and the recursion, for fixed P and M is
c
c                             ( N + 1 ) * K(N+1,P,X,M) =
c        ( X - ( N + P * ( M - 2 * N))) * K(N,  P,X,M)
c       - ( M - N + 1 ) * P * ( 1 - P ) * K(N-1,P,X,M)
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    18 March 2009
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Walter Gautschi,
c    Orthogonal Polynomials: Computation and Approximation,
c    Oxford, 2004,
c    ISBN: 0-19-850672-4,
c    LC: QA404.5 G3555.
c
c  Parameters:
c
c    Input, integer N, the highest order polynomial to evaluate.
c    0 <= N.
c
c    Input, double precision P, the parameter.  0 < P < 1.
c
c    Input, double precision X, the evaluation parameter.
c
c    Input, integer M, the parameter.  0 <= M.
c
c    Output, double precision V(0:N), the values of the Krawtchouk polynomials
c    of orders 0 through N at X.
c
      implicit none

      integer  n

      integer i
      integer m
      double precision p
      double precision x
      double precision v(0:n)

      if ( n .lt. 0 ) then
        write ( * , '(a)' ) ' '
        write ( * , '(a)' ) 'KRAWTCHOUK - Fatal error!'
        write ( * , '(a)' ) '  0 <= N is required.'
        stop
      end if

      if ( p .le. 0.0 .or. 1.0 .le. p ) then
        write ( * , '(a)' ) ' '
        write ( * , '(a)' ) 'KRAWTCHOUK - Fatal error!'
        write ( * , '(a)' ) '  0 < P < 1 is required.'
        stop
      end if

      if ( m .lt. 0 ) then
        write ( * , '(a)' ) ' '
        write ( * , '(a)' ) 'KRAWTCHOUK - Fatal error!'
        write ( * , '(a)' ) '  0 <= M is required.'
        stop
      end if

      v(0) = 1.0D+00

      if ( 1 <= n ) then
        v(1) = x - p * real ( m, kind = 8 )
      end if

      do i = 1, n - 1
        v(i+1) = ( 
     &    ( x - ( dble ( i ) + p * dble ( m - 2 * i ) ) )
     &      * v(i)
     &    - dble ( m - i + 1 ) * p * ( 1.0D+00 - p ) * v(i-1)
     &    ) / dble ( i + 1 )
      end do

      return
      end
      subroutine laguerre_associated ( n, m, x, cx )

c*********************************************************************72
c
cc LAGUERRE_ASSOCIATED evaluates associated Laguerre polynomials L(N,M,X).
c
c  Differential equation:
c
c    X Y'' + (M+1-X) Y' + (N-M) Y = 0
c
c  First terms:
c
c    M = 0
c
c    L(0,0,X) =   1
c    L(1,0,X) =  -X    +  1
c    L(2,0,X) =   X^2  -  4 X     +  2
c    L(3,0,X) =  -X^3 +  9 X^2  -  18 X    +    6
c    L(4,0,X) =   X^4 - 16 X^3 +  72 X^2  -   96 X +      24
c    L(5,0,X) =  -X^5 + 25 X^4 - 200 X^3 +  600 X^2  -  600 x    +  120
c    L(6,0,X) =   X^6 - 36 X^5 + 450 X^4 - 2400 X^3 + 5400 X^2  - 4320 X + 720
c
c    M = 1
c
c    L(0,1,X) =    0
c    L(1,1,X) =   -1,
c    L(2,1,X) =    2 X - 4,
c    L(3,1,X) =   -3 X^2  + 18 X - 18,
c    L(4,1,X) =    4 X^3 - 48 X^2 + 144 X - 96
c
c    M = 2
c
c    L(0,2,X) =    0
c    L(1,2,X) =    0,
c    L(2,2,X) =    2,
c    L(3,2,X) =   -6 X + 18,
c    L(4,2,X) =   12 X^2 - 96 X + 144
c
c    M = 3
c
c    L(0,3,X) =    0
c    L(1,3,X) =    0,
c    L(2,3,X) =    0,
c    L(3,3,X) =   -6,
c    L(4,3,X) =   24 X - 96
c
c    M = 4
c
c    L(0,4,X) =    0
c    L(1,4,X) =    0
c    L(2,4,X) =    0
c    L(3,4,X) =    0
c    L(4,4,X) =   24
c
c  Recursion:
c
c    if N = 0:
c
c      L(N,M,X)   = 0 
c
c    if N = 1:
c
c      L(N,M,X)   = (M+1-X)
c
c    if 2 <= N:
c
c      L(N,M,X)   = ( (M+2*N-1-X) * L(N-1,M,X) 
c                  +   (1-M-N)     * L(N-2,M,X) ) / N
c
c  Special values:
c
c    For M = 0, the associated Laguerre polynomials L(N,M,X) are equal 
c    to the Laguerre polynomials L(N,X).
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    13 July 2008
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
c  Parameters:
c
c    Input, integer N, the highest order polynomial to compute.
c    Note that polynomials 0 through N will be computed.
c
c    Input, integer M, the parameter.  M must be nonnegative.
c
c    Input, double precision X, the point at which the polynomials are 
c    to be evaluated.
c
c    Output, double precision CX(0:N), the associated Laguerre polynomials of 
c    degrees 0 through N evaluated at the point X.
c
      implicit none

      integer n

      double precision cx(0:n)
      integer i
      integer m
      double precision x

      if ( m .lt. 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'LAGUERRE_ASSOCIATED - Fatal error!'
        write ( *, '(a,i8)' ) '  Input value of M = ', m
        write ( *, '(a)' ) '  but M must be nonnegative.'
        stop
      end if

      if ( n .lt. 0 ) then
        return
      end if

      cx(0) = 1.0D+00

      if ( n .eq. 0 ) then
        return
      end if

      cx(1) = dble ( m + 1 ) - x

      do i = 2, n
        cx(i) = ( ( dble (   m + 2 * i - 1 ) - x ) * cx(i-1)   
     &            + dble ( - m     - i + 1 )       * cx(i-2) ) 
     &            / dble (           i     )
      end do

      return
      end
      subroutine laguerre_poly ( n, x, cx )

c*********************************************************************72
c
cc LAGUERRE_POLY evaluates the Laguerre polynomials at X.
c
c  Differential equation:
c
c    X * Y'' + (1-X) * Y' + N * Y = 0
c
c  First terms:
c
c      1
c     -X    +  1
c   (  X^2  -  4 X     +  2 ) / 2
c   ( -X^3 +  9 X^2  -  18 X    +    6 ) / 6
c   (  X^4 - 16 X^3 +  72 X^2  -   96 X +      24 ) / 24
c   ( -X^5 + 25 X^4 - 200 X^3 +  600 X^2  -  600 X    +  120 ) / 120
c   (  X^6 - 36 X^5 + 450 X^4 - 2400 X^3 + 5400 X^2 - 4320 X + 720 ) / 720
c   ( -X^7 + 49 X^6 - 882 X^5 + 7350 X^4 - 29400 X^3 
c      + 52920 X^2 - 35280 X + 5040 ) / 5040
c
c  Recursion:
c
c    L(0,X) = 1,
c    L(1,X) = 1-X,
c    N * L(N,X) = (2*N-1-X) * L(N-1,X) - (N-1) * L(N-2,X)
c
c  Orthogonality:
c
c    Integral ( 0 <= X < +oo ) exp ( - X ) * L(N,X) * L(M,X) dX
c    = 0 if N /= M
c    = 1 if N == M
c
c  Special values:
c
c    L(N,0) = 1.
c
c  Relations:
c
c    L(N,X) = (-1)**N / Nc * exp ( x ) * (d/dx)**n ( exp ( - x ) * x**n )  
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    14 July 2008
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
c  Parameters:
c
c    Input, integer N, the highest order polynomial to compute.
c    Note that polynomials 0 through N will be computed.
c
c    Input, double precision X, the point at which the polynomials are 
c    to be evaluated.
c
c    Output, double precision CX(0:N), the Laguerre polynomials of 
c    degree 0 through N evaluated at the point X.
c
      implicit none

      integer n

      double precision cx(0:n)
      integer i
      double precision x

      if ( n .lt. 0 ) then
        return
      end if

      cx(0) = 1.0D+00

      if ( n .eq. 0 ) then
        return
      end if

      cx(1) = 1.0D+00 - x
     
      do i = 2, n

        cx(i) = ( ( dble ( 2 * i - 1 ) - x ) * cx(i-1)   
     &            - dble (     i - 1 )       * cx(i-2) ) 
     &            / dble (     i     )

      end do

      return
      end
      subroutine laguerre_poly_coef ( n, c )

c*****************************************************************************80
c
cc LAGUERRE_POLY_COEF evaluates the Laguerre polynomial coefficients.
c
c  First terms:
c
c    0: 1
c    1: 1  -1
c    2: 1  -2  1/2
c    3: 1  -3  3/2  1/6
c    4: 1  -4  4   -2/3  1/24
c    5: 1  -5  5   -5/3  5/24  -1/120
c
c  Recursion:
c
c    L(0) = ( 1,  0, 0, ..., 0 )
c    L(1) = ( 1, -1, 0, ..., 0 )
c    L(N) = (2*N-1-X) * L(N-1) - (N-1) * L(N-2) / N
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    14 July 2008
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
c  Parameters:
c
c    Input, integer N, the highest order polynomial to compute.
c    Note that polynomials 0 through N will be computed.
c
c    Output, double precision C(0:N,0:N), the coefficients of the
c    Laguerre polynomials of degree 0 through N.  Each polynomial 
c   is stored as a row.
c
      implicit none

      integer n

      double precision c(0:n,0:n)
      integer i
      integer j

      if ( n .lt. 0 ) then
        return
      end if

      do i = 0, n
        c(i,0) = 1.0D+00
        do j = 1, n
          c(i,j) = 0.0D+00
        end do
      end do

      if ( n .eq. 0 ) then
        return
      end if

      c(1,1) = -1.0D+00
     
      do i = 2, n

        do j = 1, n
          c(i,j) = ( 
     &        dble ( 2 * i - 1 ) * c(i-1,j)     
     &      + dble (   - i + 1 ) * c(i-2,j)     
     &      -                      c(i-1,j-1) ) 
     &      / dble (     i     )
        end do
      end do

      return
      end
      subroutine laguerre_polynomial_values ( n_data, n, x, fx )

c*********************************************************************72
c
cc LAGUERRE_POLYNOMIAL_VALUES returns some values of the Laguerre polynomial.
c
c  Discussion:
c
c    In Mathematica, the function can be evaluated by:
c
c      LaguerreL[n,x]
c
c  Differential equation:
c
c    X * Y'' + (1-X) * Y' + N * Y = 0
c
c  First terms:
c
c      1
c     -X    +  1
c   (  X^2 -  4 X     +  2 ) / 2
c   ( -X^3 +  9 X^2 -  18 X    +    6 ) / 6
c   (  X^4 - 16 X^3 +  72 X^2 -   96 X +      24 ) / 24
c   ( -X^5 + 25 X^4 - 200 X^3 +  600 X^2 -  600 x    +  120 ) / 120
c   (  X^6 - 36 X^5 + 450 X^4 - 2400 X^3 + 5400 X^2 - 4320 X + 720 ) / 720
c   ( -X^7 + 49 X^6 - 882 X^5 + 7350 X^4 - 29400 X^3
c      + 52920 X^2 - 35280 X + 5040 ) / 5040
c
c  Recursion:
c
c    L(0,X) = 1,
c    L(1,X) = 1-X,
c    N * L(N,X) = (2*N-1-X) * L(N-1,X) - (N-1) * L(N-2,X)
c
c  Orthogonality:
c
c    Integral ( 0 <= X .lt. +oo ) exp ( - X ) * L(N,X) * L(M,X) dX
c    = 0 if N /= M
c    = 1 if N == M
c
c  Special values:
c
c    L(N,0) = 1.
c
c  Relations:
c
c    L(N,X) = (-1)**N / Nc * exp ( x ) * (d/dx)**n ( exp ( - x ) * x**n )
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    24 March 2007
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
c    Output, integer N, the order of the polynomial.
c
c    Output, double precision X, the point where the polynomial is evaluated.
c
c    Output, double precision FX, the value of the function.
c
      implicit none

      integer n_max
      parameter ( n_max = 17 )

      double precision fx
      double precision fx_vec(n_max)
      integer n
      integer n_data
      integer n_vec(n_max)
      double precision x
      double precision x_vec(n_max)

      save fx_vec
      save n_vec
      save x_vec

      data fx_vec /
     &   0.1000000000000000D+01,
     &   0.0000000000000000D+00,
     &  -0.5000000000000000D+00,
     &  -0.6666666666666667D+00,
     &  -0.6250000000000000D+00,
     &  -0.4666666666666667D+00,
     &  -0.2569444444444444D+00,
     &  -0.4047619047619048D-01,
     &   0.1539930555555556D+00,
     &   0.3097442680776014D+00,
     &   0.4189459325396825D+00,
     &   0.4801341790925124D+00,
     &   0.4962122235082305D+00,
     &  -0.4455729166666667D+00,
     &   0.8500000000000000D+00,
     &  -0.3166666666666667D+01,
     &   0.3433333333333333D+02 /
      data n_vec /
     &   0,  1,  2,
     &   3,  4,  5,
     &   6,  7,  8,
     &   9, 10, 11,
     &  12,  5,  5,
     &   5,  5 /
      data x_vec /
     &  1.0D+00,
     &  1.0D+00,
     &  1.0D+00,
     &  1.0D+00,
     &  1.0D+00,
     &  1.0D+00,
     &  1.0D+00,
     &  1.0D+00,
     &  1.0D+00,
     &  1.0D+00,
     &  1.0D+00,
     &  1.0D+00,
     &  1.0D+00,
     &  0.5D+00,
     &  3.0D+00,
     &  5.0D+00,
     &  1.0D+01 /

      if ( n_data .lt. 0 ) then
        n_data = 0
      end if

      n_data = n_data + 1

      if ( n_max .lt. n_data ) then
        n_data = 0
        n = 0
        x = 0.0D+00
        fx = 0.0D+00
      else
        n = n_vec(n_data)
        x = x_vec(n_data)
        fx = fx_vec(n_data)
      end if

      return
      end
      function lambert_w ( x )

c*********************************************************************72
c
cc LAMBERT_W estimates the Lambert W function.
c
c  Discussion:
c
c    The function W(X) is defined implicitly by:
c
c      W(X) * e^W(X) = X
c
c    The function is also known as the "Omega" function.
c
c    In Mathematica, the function can be evaluated by:
c
c      W = ProductLog [ X ]
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    11 July 2008
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Robert Corless, Gaston Gonnet, David Hare, David Jeffrey, Donald Knuth,
c    On the Lambert W Function,
c    Advances in Computational Mathematics,
c    Volume 5, 1996, pages 329-359.
c
c    Brian Hayes,
c    "Why W?",
c    The American Scientist,
c    Volume 93, March-April 2005, pages 104-108.
c
c    Eric Weisstein,
c    CRC Concise Encyclopedia of Mathematics,
c    CRC Press, 2002,
c    Second edition,
c    ISBN: 1584883472,
c    LC: QA5.W45
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
c    Input, double precision X, the argument of the function.
c
c    Output, double precision LAMBERT_W, an approximation to the 
c    Lambert W function.
c
      implicit none

      double precision lambert_w
      double precision lambert_w_crude
      integer it
      integer it_max 
      parameter ( it_max = 100 )
      double precision tol
      parameter ( tol = 1.0D-10 )
      double precision w
      double precision x

      w = lambert_w_crude ( x )
      it = 0

10    continue

        if ( it_max .lt. it ) then
          go to 20
        end if

        if ( abs ( ( x - w * exp ( w ) ) ) .lt.
     &    tol * abs ( ( w + 1.0D+00 ) * exp ( w ) ) ) then
          go to 20
        end if

        w = w - ( w * exp ( w ) - x ) 
     &    / ( ( w + 1.0D+00 ) * exp ( w ) 
     &    - ( w + 2.0D+00 ) * ( w * exp ( w ) - x ) 
     &    / ( 2.0D+00 * w + 2.0D+00 ) )

        it = it + 1

      go to 10

20    continue

      lambert_w = w

      return
      end
      function lambert_w_crude ( x )

c*********************************************************************72
c
cc LAMBERT_W_CRUDE is a crude estimate of the Lambert W function.
c
c  Discussion:
c
c    This crude approximation can be used as a good starting point
c    for an iterative process.
c
c    The function W(X) is defined implicitly by:
c
c      W(X) * e^W(X) = X
c
c    The function is also known as the "Omega" function.
c
c    In Mathematica, the function can be evaluated by:
c
c      W = ProductLog [ X ]
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    13 July 2008
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Robert Corless, Gaston Gonnet, David Hare, David Jeffrey, Donald Knuth,
c    On the Lambert W Function,
c    Advances in Computational Mathematics,
c    Volume 5, 1996, pages 329-359.
c
c    Brian Hayes,
c    "Why W?",
c    The American Scientist,
c    Volume 93, March-April 2005, pages 104-108.
c
c    Eric Weisstein,
c    CRC Concise Encyclopedia of Mathematics,
c    CRC Press, 2002,
c    Second edition,
c    ISBN: 1584883472,
c    LC: QA5.W45
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
c    Input, double precision X, the argument of the function.
c
c    Output, double precision LAMBERT_W_CRUDE, a crude approximation 
c    to the Lambert W function.
c
      implicit none

      double precision lambert_w_crude
      double precision value
      double precision x

      if ( x .le. 500.0D+00 ) then

        value = 0.04D+00 + 0.665D+00 
     &    * ( 1.0D+00 + 0.0195D+00 * log ( x + 1.0D+00 ) ) 
     &    * log ( x + 1.0D+00 )

      else

        value = log ( x - 4.0D+00 ) 
     &    - ( 1.0D+00 - 1.0D+00 / log ( x ) ) * log ( log ( x ) )

      end if

      lambert_w_crude = value

      return
      end
      subroutine lambert_w_values ( n_data, x, fx )

c*********************************************************************72
c
cc LAMBERT_W_VALUES returns some values of the Lambert W function.
c
c  Discussion:
c
c    The function W(X) is defined implicitly by:
c
c      W(X) * e^W(X) = X
c
c    The function is also known as the "Omega" function.
c
c    In Mathematica, the function can be evaluated by:
c
c      W = ProductLog [ X ]
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    23 February 2005
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    R M Corless, G H Gonnet, D E Hare, D J Jeffrey, D E Knuth,
c    On the Lambert W Function,
c    Advances in Computational Mathematics,
c    Volume 5, 1996, pages 329-359.
c
c    Brian Hayes,
c    "Why W?",
c    The American Scientist,
c    Volume 93, March-April 2005, pages 104-108.
c
c    Eric Weisstein,
c    "Lambert's W-Function",
c    CRC Concise Encyclopedia of Mathematics,
c    CRC Press, 1998.
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
c    Output, real X, the argument of the function.
c
c    Output, real FX, the value of the function.
c
      implicit none

      integer n_max
      parameter ( n_max = 22 )

      real fx
      real fx_vec(n_max)
      integer n_data
      real x
      real x_vec(n_max)

      save fx_vec
      save x_vec

      data fx_vec /
     &  0.0000000000000000D+00,
     &  0.3517337112491958D+00,
     &  0.5671432904097839D+00,
     &  0.7258613577662263D+00,
     &  0.8526055020137255D+00,
     &  0.9585863567287029D+00,
     &  0.1000000000000000D+01,
     &  0.1049908894964040D+01,
     &  0.1130289326974136D+01,
     &  0.1202167873197043D+01,
     &  0.1267237814307435D+01,
     &  0.1326724665242200D+01,
     &  0.1381545379445041D+01,
     &  0.1432404775898300D+01,
     &  0.1479856830173851D+01,
     &  0.1524345204984144D+01,
     &  0.1566230953782388D+01,
     &  0.1605811996320178D+01,
     &  0.1745528002740699D+01,
     &  0.3385630140290050D+01,
     &  0.5249602852401596D+01,
     &  0.1138335808614005D+02 /
      data x_vec /
     &  0.0000000000000000D+00,
     &  0.5000000000000000D+00,
     &  0.1000000000000000D+01,
     &  0.1500000000000000D+01,
     &  0.2000000000000000D+01,
     &  0.2500000000000000D+01,
     &  0.2718281828459045D+01,
     &  0.3000000000000000D+01,
     &  0.3500000000000000D+01,
     &  0.4000000000000000D+01,
     &  0.4500000000000000D+01,
     &  0.5000000000000000D+01,
     &  0.5500000000000000D+01,
     &  0.6000000000000000D+01,
     &  0.6500000000000000D+01,
     &  0.7000000000000000D+01,
     &  0.7500000000000000D+01,
     &  0.8000000000000000D+01,
     &  0.1000000000000000D+02,
     &  0.1000000000000000D+03,
     &  0.1000000000000000D+04,
     &  0.1000000000000000D+07 /

      if ( n_data .lt. 0 ) then
        n_data = 0
      end if

      n_data = n_data + 1

      if ( n_max .lt. n_data ) then
        n_data = 0
        x = 0.0
        fx = 0.0
      else
        x = x_vec(n_data)
        fx = fx_vec(n_data)
      end if

      return
      end
      subroutine legendre_associated ( n, m, x, cx )

c*********************************************************************72
c
cc LEGENDRE_ASSOCIATED evaluates the associated Legendre functions.
c
c  Differential equation:
c
c    (1-X*X) * Y'' - 2 * X * Y + ( N (N+1) - (M*M/(1-X*X)) * Y = 0
c
c  First terms:
c
c    M = 0  ( = Legendre polynomials of first kind P(N,X) )
c
c    P00 =    1
c    P10 =    1 X
c    P20 = (  3 X^2 -   1)/2
c    P30 = (  5 X^3 -   3 X)/2
c    P40 = ( 35 X^4 -  30 X^2 +   3)/8
c    P50 = ( 63 X^5 -  70 X^3 +  15 X)/8
c    P60 = (231 X^6 - 315 X^4 + 105 X^2 -  5)/16
c    P70 = (429 X^7 - 693 X^5 + 315 X^3 - 35 X)/16
c
c    M = 1
c
c    P01 =   0
c    P11 =   1 * SQRT(1-X*X)
c    P21 =   3 * SQRT(1-X*X) * X
c    P31 = 1.5 * SQRT(1-X*X) * (5*X*X-1)
c    P41 = 2.5 * SQRT(1-X*X) * (7*X*X*X-3*X)
c
c    M = 2
c
c    P02 =   0
c    P12 =   0
c    P22 =   3 * (1-X*X)
c    P32 =  15 * (1-X*X) * X
c    P42 = 7.5 * (1-X*X) * (7*X*X-1)
c
c    M = 3
c
c    P03 =   0
c    P13 =   0
c    P23 =   0
c    P33 =  15 * (1-X*X)**1.5
c    P43 = 105 * (1-X*X)**1.5 * X
c
c    M = 4
c
c    P04 =   0
c    P14 =   0
c    P24 =   0
c    P34 =   0
c    P44 = 105 * (1-X*X)^2
c
c  Recursion:
c
c    if N < M:
c      P(N,M) = 0
c    if N = M:
c      P(N,M) = (2*M-1)!! * (1-X*X)**(M/2) where N!! means the product of
c      all the odd integers less than or equal to N.
c    if N = M+1:
c      P(N,M) = X*(2*M+1)*P(M,M)
c    if M+1 < N:
c      P(N,M) = ( X*(2*N-1)*P(N-1,M) - (N+M-1)*P(N-2,M) )/(N-M)
c
c  Special values:
c
c    P(N,0,X) = P(N,X), that is, for M=0, the associated Legendre
c    function of the first kind equals the Legendre polynomial of the
c    first kind.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    17 July 2008
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
c  Parameters:
c
c    Input, integer N, the maximum first index of the Legendre
c    function, which must be at least 0.
c
c    Input, integer M, the second index of the Legendre function,
c    which must be at least 0, and no greater than N.
c
c    Input, double precision X, the point at which the function is to be
c    evaluated.  X must satisfy -1 <= X <= 1.
c
c    Output, double precision CX(0:N), the values of the first N+1 functions.
c
      implicit none

      integer n

      double precision cx(0:n)
      double precision fact
      integer i
      integer m
      double precision somx2
      double precision x

      if ( m .lt. 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'LEGENDRE_ASSOCIATED - Fatal error!'
        write ( *, '(a,i8)' ) '  Input value of M is ', m
        write ( *, '(a)' ) '  but M must be nonnegative.'
        stop
      end if
 
      if ( n .lt. m ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'LEGENDRE_ASSOCIATED - Fatal error!'
        write ( *, '(a,i8)' ) '  Input value of M = ', m
        write ( *, '(a,i8)' ) '  Input value of N = ', n
        write ( *, '(a)' ) '  but M must be less than or equal to N.'
        stop
      end if
 
      if ( x .lt. -1.0D+00 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'LEGENDRE_ASSOCIATED - Fatal error!'
        write ( *, '(a,g14.6)' ) '  Input value of X = ', x
        write ( *, '(a)' ) '  but X must be no less than -1.'
        stop
      end if
 
      if ( 1.0D+00 .lt. x ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'LEGENDRE_ASSOCIATED - Fatal error!'
        write ( *, '(a,g14.6)' ) '  Input value of X = ', x
        write ( *, '(a)' ) '  but X must be no more than 1.'
        stop
      end if
  
      do i = 0, m - 1
        cx(i) = 0.0D+00
      end do

      cx(m) = 1.0D+00
      somx2 = sqrt ( 1.0D+00 - x * x )
 
      fact = 1.0D+00
      do i = 1, m
        cx(m) = -cx(m) * fact * somx2
        fact = fact + 2.0D+00
      end do
 
      if ( m + 1 .le. n ) then
        cx(m+1) = x * dble ( 2 * m + 1 ) * cx(m)
      end if

      do i = m+2, n
        cx(i) = ( dble ( 2 * i     - 1 ) * x * cx(i-1) 
     &          + dble (   - i - m + 1 ) *     cx(i-2) ) 
     &          / dble (     i - m     )
      end do

      return
      end
      subroutine legendre_associated_normalized ( n, m, x, cx )

c*********************************************************************72
c
cc LEGENDRE_ASSOCIATED_NORMALIZED: normalized associated Legendre functions.
c
c  Discussion:
c
c    The unnormalized associated Legendre functions P_N^M(X) have
c    the property that
c
c      Integral ( -1 <= X <= 1 ) ( P_N^M(X) )^2 dX 
c      = 2 * ( N + M )c / ( ( 2 * N + 1 ) * ( N - M )c )
c
c    By dividing the function by the square root of this term,
c    the normalized associated Legendre functions have norm 1.
c
c    However, we plan to use these functions to build spherical
c    harmonics, so we use a slightly different normalization factor of
c
c      sqrt ( ( ( 2 * N + 1 ) * ( N - M )! ) / ( 4 * pi * ( N + M )! ) ) 
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    18 July 2008
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
c  Parameters:
c
c    Input, integer N, the maximum first index of the Legendre
c    function, which must be at least 0.
c
c    Input, integer M, the second index of the Legendre function,
c    which must be at least 0, and no greater than N.
c
c    Input, double precision X, the point at which the function is to be
c    evaluated.  X must satisfy -1 <= X <= 1.
c
c    Output, double precision CX(0:N), the values of the first N+1 functions.
c
      implicit none

      integer n

      double precision cx(0:n)
      double precision factor
      integer i
      integer m
      integer mm
      double precision pi
      parameter ( pi = 3.141592653589793D+00 )
      double precision r8_factorial
      double precision somx2
      double precision x

      if ( m .lt. 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 
     &    'LEGENDRE_ASSOCIATED_NORMALIZED - Fatal error!'
        write ( *, '(a,i8)' ) '  Input value of M is ', m
        write ( *, '(a)' ) '  but M must be nonnegative.'
        stop
      end if
     
      if ( n .lt. m ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 
     &    'LEGENDRE_ASSOCIATED_NORMALIZED - Fatal error!'
        write ( *, '(a,i8)' ) '  Input value of M = ', m
        write ( *, '(a,i8)' ) '  Input value of N = ', n
        write ( *, '(a)' ) '  but M must be less than or equal to N.'
        stop
      end if
     
      if ( x .lt. -1.0D+00 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 
     &    'LEGENDRE_ASSOCIATED_NORMALIZED - Fatal error!'
        write ( *, '(a,g14.6)' ) '  Input value of X = ', x
        write ( *, '(a)' ) '  but X must be no less than -1.'
        stop
      end if
     
      if ( 1.0D+00 .lt. x ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 
     &    'LEGENDRE_ASSOCIATED_NORMALIZED - Fatal error!'
        write ( *, '(a,g14.6)' ) '  Input value of X = ', x
        write ( *, '(a)' ) '  but X must be no more than 1.'
        stop
      end if
c
c  Entries 0 through M-1 are zero.
c
      do i = 0, m - 1
        cx(i) = 0.0D+00
      end do
      cx(m) = 1.0D+00
      somx2 = sqrt ( 1.0D+00 - x * x )
     
      factor = 1.0D+00
      do i = 1, m
        cx(m) = - cx(m) * factor * somx2
        factor = factor + 2.0D+00
      end do
     
      if ( m + 1 .le. n ) then
        cx(m+1) = x * dble ( 2 * m + 1 ) * cx(m)
      end if

      do i = m+2, n
        cx(i) = ( dble ( 2 * i     - 1 ) * x * cx(i-1) 
     &          + dble (   - i - m + 1 ) *     cx(i-2) ) 
     &          / dble (     i - m     )
      end do
c
c  Normalization.
c
      do mm = m, n
        factor = sqrt ( ( dble ( 2 * mm + 1 ) 
     &    * r8_factorial ( mm - m ) ) 
     &    / ( 4.0D+00 * pi * r8_factorial ( mm + m ) ) )
        cx(mm) = cx(mm) * factor
      end do

      return
      end
      subroutine legendre_associated_normalized_values ( n_data, n, m, 
     &  x, fx )

c*********************************************************************72
c
cc LEGENDRE_ASSOCIATED_NORMALIZED_VALUES: normalized associated Legendre.
c
c  Discussion:
c
c    The function considered is the associated Legendre polynomial P^M_N(X).
c
c    In Mathematica, the function can be evaluated by:
c
c      LegendreP [ n, m, x ]
c
c    The function is normalized by dividing by 
c
c      sqrt ( 4 * pi * ( n + m )! / ( 2 * n + 1 ) / ( n - m )! )
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    01 September 2010
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
c    Cambridge University Press, 1999,
c    ISBN: 0-521-64314-7,
c    LC: QA76.95.W65.
c
c  Parameters:
c
c    Input/output, integer N_DATA.  The user sets N_DATA to 0 
c    before the first call.  On each call, the routine increments N_DATA by 1,
c    and returns the corresponding data; when there is no more data, the
c    output value of N_DATA will be 0 again.
c
c    Output, integer N, integer M, double precision X, 
c    the arguments of the function.
c
c    Output, double precision FX, the value of the function.
c
      implicit none

      integer n_max 
      parameter ( n_max = 21 )

      double precision fx
      double precision fx_vec(n_max)
      integer m
      integer m_vec(n_max)
      integer n
      integer n_data
      integer n_vec(n_max)
      double precision x
      double precision x_vec(n_max)

      save fx_vec
      save m_vec
      save n_vec
      save x_vec

      data fx_vec /
     &   0.2820947917738781D+00, 
     &   0.2443012559514600D+00, 
     &  -0.2992067103010745D+00, 
     &  -0.07884789131313000D+00, 
     &  -0.3345232717786446D+00, 
     &   0.2897056515173922D+00, 
     &  -0.3265292910163510D+00, 
     &  -0.06997056236064664D+00, 
     &   0.3832445536624809D+00, 
     &  -0.2709948227475519D+00, 
     &  -0.2446290772414100D+00, 
     &   0.2560660384200185D+00, 
     &   0.1881693403754876D+00, 
     &  -0.4064922341213279D+00, 
     &   0.2489246395003027D+00, 
     &   0.08405804426339821D+00, 
     &   0.3293793022891428D+00, 
     &  -0.1588847984307093D+00, 
     &  -0.2808712959945307D+00, 
     &   0.4127948151484925D+00, 
     &  -0.2260970318780046D+00 /
      data m_vec /
     &  0, 0, 1, 0, 
     &  1, 2, 0, 1, 
     &  2, 3, 0, 1, 
     &  2, 3, 4, 0, 
     &  1, 2, 3, 4, 
     &  5 /
      data n_vec /
     &  0,  1,  1,  2, 
     &  2,  2,  3,  3, 
     &  3,  3,  4,  4, 
     &  4,  4,  4,  5, 
     &  5,  5,  5,  5, 
     &  5 /
      data x_vec /
     &  0.50D+00, 
     &  0.50D+00,
     &  0.50D+00,
     &  0.50D+00, 
     &  0.50D+00, 
     &  0.50D+00, 
     &  0.50D+00, 
     &  0.50D+00, 
     &  0.50D+00, 
     &  0.50D+00, 
     &  0.50D+00, 
     &  0.50D+00, 
     &  0.50D+00, 
     &  0.50D+00, 
     &  0.50D+00, 
     &  0.50D+00, 
     &  0.50D+00, 
     &  0.50D+00, 
     &  0.50D+00, 
     &  0.50D+00, 
     &  0.50D+00 /

      if ( n_data .lt. 0 ) then
        n_data = 0
      end if

      n_data = n_data + 1

      if ( n_max .lt. n_data ) then
        n_data = 0
        n = 0
        m = 0
        x = 0.0D+00
        fx = 0.0D+00
      else
        n = n_vec(n_data)
        m = m_vec(n_data)
        x = x_vec(n_data)
        fx = fx_vec(n_data)
      end if

      return
      end
      subroutine legendre_associated_values ( n_data, n, m, x, fx )

c*********************************************************************72
c
cc LEGENDRE_ASSOCIATED_VALUES returns values of associated Legendre functions.
c
c  Discussion:
c
c    The function considered is the associated Legendre polynomial P^M_N(X).
c
c    In Mathematica, the function can be evaluated by:
c
c      LegendreP [ n, m, x ]
c
c  Differential equation:
c
c    (1-X*X) * Y'' - 2 * X * Y + ( N (N+1) - (M*M/(1-X*X)) * Y = 0
c
c  First terms:
c
c    M = 0  ( = Legendre polynomials of first kind P(N,X) )
c
c    P00 =    1
c    P10 =    1 X
c    P20 = (  3 X^2 -   1)/2
c    P30 = (  5 X^3 -   3 X)/2
c    P40 = ( 35 X^4 -  30 X^2 +   3)/8
c    P50 = ( 63 X^5 -  70 X^3 +  15 X)/8
c    P60 = (231 X^6 - 315 X^4 + 105 X^2 -  5)/16
c    P70 = (429 X^7 - 693 X^5 + 315 X^3 - 35 X)/16
c
c    M = 1
c
c    P01 =   0
c    P11 =   1 * SQRT(1-X*X)
c    P21 =   3 * SQRT(1-X*X) * X
c    P31 = 1.5 * SQRT(1-X*X) * (5*X*X-1)
c    P41 = 2.5 * SQRT(1-X*X) * (7*X*X*X-3*X)
c
c    M = 2
c
c    P02 =   0
c    P12 =   0
c    P22 =   3 * (1-X*X)
c    P32 =  15 * (1-X*X) * X
c    P42 = 7.5 * (1-X*X) * (7*X*X-1)
c
c    M = 3
c
c    P03 =   0
c    P13 =   0
c    P23 =   0
c    P33 =  15 * (1-X*X)**1.5
c    P43 = 105 * (1-X*X)**1.5 * X
c
c    M = 4
c
c    P04 =   0
c    P14 =   0
c    P24 =   0
c    P34 =   0
c    P44 = 105 * (1-X*X)^2
c
c  Recursion:
c
c    if N .lt. M:
c      P(N,M) = 0
c    if N = M:
c      P(N,M) = (2*M-1)cc * (1-X*X)**(M/2) where Ncc means the product of
c      all the odd integers less than or equal to N.
c    if N = M+1:
c      P(N,M) = X*(2*M+1)*P(M,M)
c    if M+1 .lt. N:
c      P(N,M) = ( X*(2*N-1)*P(N-1,M) - (N+M-1)*P(N-2,M) )/(N-M)
c
c  Restrictions:
c
c    -1 <= X <= 1
c     0 <= M <= N
c
c  Special values:
c
c    P(N,0,X) = P(N,X), that is, for M=0, the associated Legendre
c    polynomial of the first kind equals the Legendre polynomial of the
c    first kind.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    24 March 2007
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
c    Output, integer N, integer M, double precision X,
c    the arguments of the function.
c
c    Output, double precision FX, the value of the function.
c
      implicit none

      integer n_max
      parameter ( n_max = 20 )

      double precision fx
      double precision fx_vec(n_max)
      integer m
      integer m_vec(n_max)
      integer n
      integer n_data
      integer n_vec(n_max)
      double precision x
      double precision x_vec(n_max)

      save fx_vec
      save m_vec
      save n_vec
      save x_vec

      data fx_vec /
     &   0.0000000000000000D+00,
     &  -0.5000000000000000D+00,
     &   0.0000000000000000D+00,
     &   0.3750000000000000D+00,
     &   0.0000000000000000D+00,
     &  -0.8660254037844386D+00,
     &  -0.1299038105676658D+01,
     &  -0.3247595264191645D+00,
     &   0.1353164693413185D+01,
     &  -0.2800000000000000D+00,
     &   0.1175755076535925D+01,
     &   0.2880000000000000D+01,
     &  -0.1410906091843111D+02,
     &  -0.3955078125000000D+01,
     &  -0.9997558593750000D+01,
     &   0.8265311444100484D+02,
     &   0.2024442836815152D+02,
     &  -0.4237997531890869D+03,
     &   0.1638320624828339D+04,
     &  -0.2025687389227225D+05 /
      data m_vec /
     &  0, 0, 0, 0,
     &  0, 1, 1, 1,
     &  1, 0, 1, 2,
     &  3, 2, 2, 3,
     &  3, 4, 4, 5 /
      data n_vec /
     &  1,  2,  3,  4,
     &  5,  1,  2,  3,
     &  4,  3,  3,  3,
     &  3,  4,  5,  6,
     &  7,  8,  9, 10 /
      data x_vec /
     &  0.00D+00,
     &  0.00D+00,
     &  0.00D+00,
     &  0.00D+00,
     &  0.00D+00,
     &  0.50D+00,
     &  0.50D+00,
     &  0.50D+00,
     &  0.50D+00,
     &  0.20D+00,
     &  0.20D+00,
     &  0.20D+00,
     &  0.20D+00,
     &  0.25D+00,
     &  0.25D+00,
     &  0.25D+00,
     &  0.25D+00,
     &  0.25D+00,
     &  0.25D+00,
     &  0.25D+00 /

      if ( n_data .lt. 0 ) then
        n_data = 0
      end if

      n_data = n_data + 1

      if ( n_max .lt. n_data ) then
        n_data = 0
        n = 0
        m = 0
        x = 0.0D+00
        fx = 0.0D+00
      else
        n = n_vec(n_data)
        m = m_vec(n_data)
        x = x_vec(n_data)
        fx = fx_vec(n_data)
      end if

      return
      end
      subroutine legendre_function_q ( n, x, cx )

c*********************************************************************72
c
cc LEGENDRE_FUNCTION_Q evaluates the Legendre Q functions.
c
c  Differential equation:
c
c    (1-X*X) Y'' - 2 X Y' + N (N+1) = 0
c
c  First terms:
c
c    Q(0,X) = 0.5 * log((1+X)/(1-X))
c    Q(1,X) = Q(0,X)*X - 1 
c    Q(2,X) = Q(0,X)*(3*X*X-1)/4 - 1.5*X
c    Q(3,X) = Q(0,X)*(5*X*X*X-3*X)/4 - 2.5*X^2 + 2/3
c    Q(4,X) = Q(0,X)*(35*X^4-30*X^2+3)/16 - 35/8 * X^3 + 55/24 * X
c    Q(5,X) = Q(0,X)*(63*X^5-70*X^3+15*X)/16 - 63/8*X^4 + 49/8*X^2 - 8/15
c
c  Recursion:
c
c    Q(0) = 0.5 * log ( (1+X) / (1-X) )
c    Q(1) = 0.5 * X * log ( (1+X) / (1-X) ) - 1.0
c
c    Q(N) = ( (2*N-1) * X * Q(N-1) - (N-1) * Q(N-2) ) / N
c
c  Restrictions:
c
c    -1 < X < 1
c
c  Special values:
c
c    Note that the Legendre function Q(N,X) is equal to the
c    associated Legendre function of the second kind,
c    Q(N,M,X) with M = 0.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    18 July 2008
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
c  Parameters:
c
c    Input, integer N, the highest order function to evaluate.
c
c    Input, double precision X, the point at which the functions are to be
c    evaluated.  X must satisfy -1 < X < 1.
c
c    Output, double precision CX(0:N), the values of the first N+1 Legendre
c    functions at the point X.
c
      implicit none

      integer n

      double precision cx(0:n)
      integer i
      double precision x
c
c  Check the value of X.
c
      if ( x .le. -1.0D+00 .or. 1.0D+00 .le. x ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'LEGENDRE_FUNCTION_Q - Fatal error!'
        write ( *, '(a,g14.6)' ) '  Illegal input value of X = ', x
        write ( *, '(a)' ) '  But X must be between -1 and 1.'
        stop
      end if
     
      if ( n .lt. 0 ) then
        return
      end if

      cx(0) = 0.5D+00 * log ( ( 1.0D+00 + x ) / ( 1.0D+00 - x ) )

      if ( n .eq. 0 ) then
        return
      end if

      cx(1) = x * cx(0) - 1.0D+00

      do i = 2, n
        cx(i) = ( dble ( 2 * i - 1 ) * x * cx(i-1) 
     &          + dble (   - i + 1 )     * cx(i-2) ) 
     &          / dble (     i     )
      end do 
     
      return
      end
      subroutine legendre_function_q_values ( n_data, n, x, fx )

c*********************************************************************72
c
cc LEGENDRE_FUNCTION_Q_VALUES returns values of the Legendre Q function.
c
c  Discussion:
c
c    In Mathematica, the function can be evaluated by:
c
c      LegendreQ[n,x]
c
c  Differential equation:
c
c    (1-X*X) Y'' - 2 X Y' + N (N+1) = 0
c
c  First terms:
c
c    Q(0,X) = 0.5 * log((1+X)/(1-X))
c    Q(1,X) = Q(0,X)*X - 1
c    Q(2,X) = Q(0,X)*(3*X*X-1)/4 - 1.5*X
c    Q(3,X) = Q(0,X)*(5*X*X*X-3*X)/4 - 2.5*X^2 + 2/3
c    Q(4,X) = Q(0,X)*(35*X^4-30*X^2+3)/16 - 35/8 * X^3 + 55/24 * X
c    Q(5,X) = Q(0,X)*(63*X^5-70*X^3+15*X)/16 - 63/8*X^4 + 49/8*X^2 - 8/15
c
c  Recursion:
c
c    Q(0) = 0.5 * log ( (1+X) / (1-X) )
c    Q(1) = 0.5 * X * log ( (1+X) / (1-X) ) - 1.0
c
c    Q(N) = ( (2*N-1) * X * Q(N-1) - (N-1) * Q(N-2) ) / N
c
c  Restrictions:
c
c    -1 .lt. X .lt. 1
c
c  Special values:
c
c    Note that the Legendre function Q(N,X) is equal to the
c    associated Legendre function of the second kind,
c    Q(N,M,X) with M = 0.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    24 March 2007
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
c    Output, integer N, the order of the function.
c
c    Output, double precision X, the point where the function is evaluated.
c
c    Output, double precision FX, the value of the function.
c
      implicit none

      integer n_max
      parameter ( n_max = 21 )

      double precision fx
      double precision fx_vec(n_max)
      integer n
      integer n_data
      integer n_vec(n_max)
      double precision x
      double precision x_vec(n_max)

      save fx_vec
      save n_vec
      save x_vec

      data fx_vec /
     &   0.2554128118829953D+00,
     &  -0.9361467970292512D+00,
     &  -0.4787614548274669D+00,
     &   0.4246139251747229D+00,
     &   0.5448396833845414D+00,
     &  -0.9451328261673470D-01,
     &  -0.4973516573531213D+00,
     &  -0.1499018843853194D+00,
     &   0.3649161918783626D+00,
     &   0.3055676545072885D+00,
     &  -0.1832799367995643D+00,
     &   0.6666666666666667D+00,
     &   0.6268672028763330D+00,
     &   0.5099015515315237D+00,
     &   0.3232754180589764D+00,
     &   0.8026113738148187D-01,
     &  -0.1986547714794823D+00,
     &  -0.4828663183349136D+00,
     &  -0.7252886849144386D+00,
     &  -0.8454443502398846D+00,
     &  -0.6627096245052618D+00 /
      data n_vec /
     &   0,  1,  2,
     &   3,  4,  5,
     &   6,  7,  8,
     &   9, 10,  3,
     &   3,  3,  3,
     &   3,  3,  3,
     &   3,  3,  3 /
      data x_vec /
     &  0.25D+00,
     &  0.25D+00,
     &  0.25D+00,
     &  0.25D+00,
     &  0.25D+00,
     &  0.25D+00,
     &  0.25D+00,
     &  0.25D+00,
     &  0.25D+00,
     &  0.25D+00,
     &  0.25D+00,
     &  0.00D+00,
     &  0.10D+00,
     &  0.20D+00,
     &  0.30D+00,
     &  0.40D+00,
     &  0.50D+00,
     &  0.60D+00,
     &  0.70D+00,
     &  0.80D+00,
     &  0.90D+00 /

      if ( n_data .lt. 0 ) then
        n_data = 0
      end if

      n_data = n_data + 1

      if ( n_max .lt. n_data ) then
        n_data = 0
        n = 0
        x = 0.0D+00
        fx = 0.0D+00
      else
        n = n_vec(n_data)
        x = x_vec(n_data)
        fx = fx_vec(n_data)
      end if

      return
      end
      subroutine legendre_poly ( n, x, cx, cpx )

c*********************************************************************72
c
cc LEGENDRE_POLY evaluates the Legendre polynomials P(N,X) at X.
c
c  Discussion:
c
c    P(N,1) = 1.
c    P(N,-1) = (-1)**N.
c    | P(N,X) | <= 1 in [-1,1].
c
c    P(N,0,X) = P(N,X), that is, for M=0, the associated Legendre
c    function of the first kind and order N equals the Legendre polynomial
c    of the first kind and order N.
c
c    The N zeroes of P(N,X) are the abscissas used for Gauss-Legendre
c    quadrature of the integral of a function F(X) with weight function 1
c    over the interval [-1,1].
c
c    The Legendre polynomials are orthonormal under the inner product defined
c    as integration from -1 to 1:
c
c      Integral ( -1 <= X <= 1 ) P(I,X) * P(J,X) dX 
c        = 0 if I =/= J
c        = 2 / ( 2*I+1 ) if I = J.
c
c    Except for P(0,X), the integral of P(I,X) from -1 to 1 is 0.
c
c    A function F(X) defined on [-1,1] may be approximated by the series
c      C0*P(0,X) + C1*P(1,X) + ... + CN*P(N,X)
c    where
c      C(I) = (2*I+1)/(2) * Integral ( -1 <= X <= 1 ) F(X) P(I,X) dx.
c
c    The formula is:
c
c      P(N,X) = (1/2**N) * sum ( 0 <= M <= N/2 ) C(N,M) C(2N-2M,N) X^(N-2*M)
c
c  Differential equation:
c
c    (1-X*X) * P(N,X)'' - 2 * X * P(N,X)' + N * (N+1) = 0
c
c  First terms:
c
c    P( 0,X) =       1
c    P( 1,X) =       1 X
c    P( 2,X) =  (    3 X^2 -       1)/2
c    P( 3,X) =  (    5 X^3 -     3 X)/2
c    P( 4,X) =  (   35 X^4 -    30 X^2 +     3)/8
c    P( 5,X) =  (   63 X^5 -    70 X^3 +    15 X)/8
c    P( 6,X) =  (  231 X^6 -   315 X^4 +   105 X^2 -     5)/16
c    P( 7,X) =  (  429 X^7 -   693 X^5 +   315 X^3 -    35 X)/16
c    P( 8,X) =  ( 6435 X^8 - 12012 X^6 +  6930 X^4 -  1260 X^2 +   35)/128
c    P( 9,X) =  (12155 X^9 - 25740 X^7 + 18018 X^5 -  4620 X^3 +  315 X)/128
c    P(10,X) =  (46189 X^10-109395 X^8 + 90090 X^6 - 30030 X^4 + 3465 X^2
c                 -63 ) /256
c
c  Recursion:
c
c    P(0,X) = 1
c    P(1,X) = X
c    P(N,X) = ( (2*N-1)*X*P(N-1,X)-(N-1)*P(N-2,X) ) / N
c
c    P'(0,X) = 0
c    P'(1,X) = 1
c    P'(N,X) = ( (2*N-1)*(P(N-1,X)+X*P'(N-1,X)-(N-1)*P'(N-2,X) ) / N
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    18 July 2008
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
c    Daniel Zwillinger, editor,
c    CRC Standard Mathematical Tables and Formulae,
c    30th Edition,
c    CRC Press, 1996.
c
c  Parameters:
c
c    Input, integer N, the highest order polynomial to evaluate.
c    Note that polynomials 0 through N will be evaluated.
c
c    Input, double precision X, the point at which the polynomials
c    are to be evaluated.
c
c    Output, double precision CX(0:N), the values of the Legendre polynomials 
c    of order 0 through N at the point X.
c
c    Output, double precision CPX(0:N), the values of the derivatives of the
c    Legendre polynomials of order 0 through N at the point X.
c
      implicit none

      integer n

      double precision cx(0:n)
      double precision cpx(0:n)
      integer i
      double precision x

      if ( n .lt. 0 ) then
        return
      end if

      cx(0) = 1.0D+00
      cpx(0) = 0.0D+00

      if ( n .lt. 1 ) then
        return
      end if

      cx(1) = x
      cpx(1) = 1.0D+00
     
      do i = 2, n
     
        cx(i) = ( dble ( 2 * i - 1 ) * x * cx(i-1)   
     &          - dble (     i - 1 ) *     cx(i-2) ) 
     &          / dble (     i     )
     
        cpx(i) = ( dble ( 2 * i - 1 ) * ( cx(i-1) + x * cpx(i-1) ) 
     &           - dble (     i - 1 ) *   cpx(i-2)               ) 
     &           / dble (     i     )
     
      end do
     
      return
      end
      subroutine legendre_poly_coef ( n, c )

c*********************************************************************72
c
cc LEGENDRE_POLY_COEF evaluates the Legendre polynomial coefficients.
c
c  First terms:
c
c     1
c     0     1
c    -1/2   0      3/2
c     0    -3/2    0     5/2
c     3/8   0    -30/8   0     35/8
c     0    15/8    0   -70/8    0     63/8
c    -5/16  0    105/16  0   -315/16   0    231/16
c     0   -35/16   0   315/16   0   -693/16   0    429/16
c
c     1.00000
c     0.00000  1.00000
c    -0.50000  0.00000  1.50000
c     0.00000 -1.50000  0.00000  2.5000
c     0.37500  0.00000 -3.75000  0.00000  4.37500
c     0.00000  1.87500  0.00000 -8.75000  0.00000  7.87500
c    -0.31250  0.00000  6.56250  0.00000 -19.6875  0.00000  14.4375
c     0.00000 -2.1875   0.00000  19.6875  0.00000 -43.3215  0.00000  26.8125
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    18 July 2008
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
c    Daniel Zwillinger, editor,
c    CRC Standard Mathematical Tables and Formulae,
c    30th Edition,
c    CRC Press, 1996.
c
c  Parameters:
c
c    Input, integer N, the highest order polynomial to evaluate.
c    Note that polynomials 0 through N will be evaluated.
c
c    Output, double precision C(0:N,0:N), the coefficients of the 
c    Legendre polynomials of degree 0 through N.  Each polynomial is 
c    stored as a row.
c
      implicit none

      integer n

      double precision c(0:n,0:n)
      integer i
      integer j

      if ( n .lt. 0 ) then
        return
      end if

      do j = 0, n
        do i = 0, n
          c(i,j) = 0.0D+00
        end do
      end do

      c(0,0) = 1.0D+00

      if ( n .le. 0 ) then
        return
      end if

      c(1,1) = 1.0D+00
     
      do i = 2, n
        do j = 0, i - 2
          c(i,j) =          dble (   - i + 1 ) * c(i-2,j) 
     &                    / dble (     i     )
        end do
        do j = 1, i
          c(i,j) = c(i,j) + dble ( i + i - 1 ) * c(i-1,j-1) 
     &                    / dble (     i     )
        end do
      end do
     
      return
      end
      subroutine legendre_poly_values ( n_data, n, x, fx )

c*********************************************************************72
c
cc LEGENDRE_POLY_VALUES returns values of the Legendre polynomials.
c
c  Discussion:
c
c    In Mathematica, the function can be evaluated by:
c
c      LegendreP [ n, x ]
c
c    The formula is:
c
c      P(N,X) = (1/2**N) * sum ( 0 <= M <= N/2 ) C(N,M) C(2N-2M,N) X^(N-2*M)
c
c  Differential equation:
c
c    (1-X*X) * P(N,X)'' - 2 * X * P(N,X)' + N * (N+1) = 0
c
c  First terms:
c
c    P( 0,X) =       1
c    P( 1,X) =       1 X
c    P( 2,X) =  (    3 X^2 -       1)/2
c    P( 3,X) =  (    5 X^3 -     3 X)/2
c    P( 4,X) =  (   35 X^4 -    30 X^2 +     3)/8
c    P( 5,X) =  (   63 X^5 -    70 X^3 +    15 X)/8
c    P( 6,X) =  (  231 X^6 -   315 X^4 +   105 X^2 -     5)/16
c    P( 7,X) =  (  429 X^7 -   693 X^5 +   315 X^3 -    35 X)/16
c    P( 8,X) =  ( 6435 X^8 - 12012 X^6 +  6930 X^4 -  1260 X^2 +   35)/128
c    P( 9,X) =  (12155 X^9 - 25740 X^7 + 18018 X^5 -  4620 X^3 +  315 X)/128
c    P(10,X) =  (46189 X^10-109395 X^8 + 90090 X^6 - 30030 X^4 + 3465 X^2
c                 -63 ) /256
c
c  Recursion:
c
c    P(0,X) = 1
c    P(1,X) = X
c    P(N,X) = ( (2*N-1)*X*P(N-1,X)-(N-1)*P(N-2,X) ) / N
c
c    P'(0,X) = 0
c    P'(1,X) = 1
c    P'(N,X) = ( (2*N-1)*(P(N-1,X)+X*P'(N-1,X)-(N-1)*P'(N-2,X) ) / N
c
c  Orthogonality:
c
c    Integral ( -1 <= X <= 1 ) P(I,X) * P(J,X) dX
c      = 0 if I =/= J
c      = 2 / ( 2*I+1 ) if I = J.
c
c  Approximation:
c
c    A function F(X) defined on [-1,1] may be approximated by the series
c
c      C0*P(0,X) + C1*P(1,X) + ... + CN*P(N,X)
c
c    where
c
c      C(I) = (2*I+1)/(2) * Integral ( -1 <= X <= 1 ) F(X) P(I,X) dx.
c
c  Special values:
c
c    P(N,1) = 1.
c    P(N,-1) = (-1)**N.
c    | P(N,X) | <= 1 in [-1,1].
c
c    P(N,0,X) = P(N,X), that is, for M=0, the associated Legendre
c    function of the first kind and order N equals the Legendre polynomial
c    of the first kind and order N.
c
c    The N zeroes of P(N,X) are the abscissas used for Gauss-Legendre
c    quadrature of the integral of a function F(X) with weight function 1
c    over the interval [-1,1].
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    24 March 2007
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
c    Output, integer N, the order of the function.
c
c    Output, double precision X, the point where the function is evaluated.
c
c    Output, double precision FX, the value of the function.
c
      implicit none

      integer n_max
      parameter ( n_max = 22 )

      double precision fx
      double precision fx_vec(n_max)
      integer n
      integer n_data
      integer n_vec(n_max)
      double precision x
      double precision x_vec(n_max)

      save fx_vec
      save n_vec
      save x_vec

      data fx_vec /
     &   0.1000000000000000D+01,
     &   0.2500000000000000D+00,
     &  -0.4062500000000000D+00,
     &  -0.3359375000000000D+00,
     &   0.1577148437500000D+00,
     &   0.3397216796875000D+00,
     &   0.2427673339843750D-01,
     &  -0.2799186706542969D+00,
     &  -0.1524540185928345D+00,
     &   0.1768244206905365D+00,
     &   0.2212002165615559D+00,
     &   0.0000000000000000D+00,
     &  -0.1475000000000000D+00,
     &  -0.2800000000000000D+00,
     &  -0.3825000000000000D+00,
     &  -0.4400000000000000D+00,
     &  -0.4375000000000000D+00,
     &  -0.3600000000000000D+00,
     &  -0.1925000000000000D+00,
     &   0.8000000000000000D-01,
     &   0.4725000000000000D+00,
     &   0.1000000000000000D+01 /
      data n_vec /
     &   0,  1,  2,
     &   3,  4,  5,
     &   6,  7,  8,
     &   9, 10,  3,
     &   3,  3,  3,
     &   3,  3,  3,
     &   3,  3,  3,
     &   3 /
      data x_vec /
     &  0.25D+00,
     &  0.25D+00,
     &  0.25D+00,
     &  0.25D+00,
     &  0.25D+00,
     &  0.25D+00,
     &  0.25D+00,
     &  0.25D+00,
     &  0.25D+00,
     &  0.25D+00,
     &  0.25D+00,
     &  0.00D+00,
     &  0.10D+00,
     &  0.20D+00,
     &  0.30D+00,
     &  0.40D+00,
     &  0.50D+00,
     &  0.60D+00,
     &  0.70D+00,
     &  0.80D+00,
     &  0.90D+00,
     &  1.00D+00 /

      if ( n_data .lt. 0 ) then
        n_data = 0
      end if

      n_data = n_data + 1

      if ( n_max .lt. n_data ) then
        n_data = 0
        n = 0
        x = 0.0D+00
        fx = 0.0D+00
      else
        n = n_vec(n_data)
        x = x_vec(n_data)
        fx = fx_vec(n_data)
      end if

      return
      end
      subroutine legendre_symbol ( q, p, l )

c*********************************************************************72
c
cc LEGENDRE_SYMBOL evaluates the Legendre symbol (Q/P).
c
c  Discussion:
c
c    Let P be an odd prime.  Q is a QUADRATIC RESIDUE modulo P
c    if there is an integer R such that R*R = Q ( mod P ).
c    The Legendre symbol ( Q / P ) is defined to be:
c
c      + 1 if Q ( mod P ) /= 0 and Q is a quadratic residue modulo P,
c      - 1 if Q ( mod P ) /= 0 and Q is not a quadratic residue modulo P,
c        0 if Q ( mod P ) .eq. 0.
c
c    We can also define ( Q / P ) for P = 2 by:
c
c      + 1 if Q ( mod P ) /= 0
c        0 if Q ( mod P ) .eq. 0
c
c  Example:
c
c    (0/7) =   0
c    (1/7) = + 1  ( 1*1 = 1 mod 7 )
c    (2/7) = + 1  ( 3*3 = 2 mod 7 )
c    (3/7) = - 1
c    (4/7) = + 1  ( 2*2 = 4 mod 7 )
c    (5/7) = - 1
c    (6/7) = - 1
c
c    Note that for any prime P, exactly half of the integers from 1 to P-1
c    are quadratic residues.
c
c    ( 0 / P ) = 0.
c
c    ( Q / P ) = ( mod ( Q, P ) / P ).
c
c    ( Q / P ) = ( Q1 / P ) * ( Q2 / P ) if Q = Q1 * Q2.
c
c    If Q is prime, and P is prime and greater than 2, then:
c
c      if ( Q .eq. 1 ) then
c
c        ( Q / P ) = 1
c
c      else if ( Q .eq. 2 ) then
c
c        ( Q / P ) = + 1 if mod ( P, 8 ) = 1 or mod ( P, 8 ) = 7,
c        ( Q / P ) = - 1 if mod ( P, 8 ) = 3 or mod ( P, 8 ) = 5.
c
c      else
c
c        ( Q / P ) = - ( P / Q ) if Q = 3 ( mod 4 ) and P = 3 ( mod 4 ),
c                  =   ( P / Q ) otherwise.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    18 July 2008
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Charles Pinter,
c    A Book of Abstract Algebra,
c    McGraw Hill, 1982, pages 236-237.
c
c    Daniel Zwillinger,
c    CRC Standard Mathematical Tables and Formulae,
c    30th Edition,
c    CRC Press, 1996, pages 86-87.
c
c  Parameters:
c
c    Input, integer Q, an integer whose Legendre symbol with
c    respect to P is desired.
c
c    Input, integer P, a prime number, greater than 1, with respect
c    to which the Legendre symbol of Q is desired.
c
c    Output, integer L, the Legendre symbol (Q/P).
c    Ordinarily, L will be -1, 0 or 1.
c    L = -2, P is less than or equal to 1.
c    L = -3, P is not prime.
c    L = -4, the internal stack of factors overflowed.
c    L = -5, not enough factorization space.
c
      implicit none

      integer maxfactor
      parameter ( maxfactor = 20 )
      integer maxstack
      parameter ( maxstack = 50 )

      integer factor(maxfactor)
      integer i
      logical i4_is_prime
      integer l
      integer nfactor
      integer nleft
      integer nmore
      integer nstack
      integer p
      integer power(maxfactor)
      integer pp
      integer pstack(maxstack)
      integer q
      integer qq
      integer qstack(maxstack)
      integer t
c
c  P must be greater than 1.
c
      if ( p .le. 1 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'LEGENDRE_SYMBOL - Fatal error!'
        write ( *, '(a)' ) '  P must be greater than 1.'
        l = -2
        return
      end if
c
c  P must be prime.
c
      if ( .not. i4_is_prime ( p ) ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'LEGENDRE_SYMBOL - Fatal error!'
        write ( *, '(a)' ) '  P is not prime.'
        l = -3
        return
      end if
c
c  ( k*P / P ) = 0.
c
      if ( mod ( q, p ) .eq. 0 ) then
        l = 0
        return
      end if
c
c  For the special case P = 2, (Q/P) = 1 for all odd numbers.
c
      if ( p .eq. 2 ) then
        l = 1
        return
      end if
c
c  Make a copy of Q, and force it to be nonnegative.
c
      qq = q

10    continue

      if ( qq .lt. 0 ) then
        qq = qq + p
        go to 10
      end if

      nstack = 0
      pp = p
      l = 1

20    continue

        qq = mod ( qq, pp )
c
c  Decompose QQ into factors of prime powers.
c
        call i4_factor ( qq, maxfactor, nfactor, factor, power, nleft )

        if ( nleft .ne. 1 ) then
          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) 'LEGENDRE_SYMBOL - Fatal error!'
          write ( *, '(a)' ) '  Not enough factorization space.'
          l = - 5
          return
        end if
c
c  Each factor which is an odd power is added to the stack.
c
        nmore = 0

        do i = 1, nfactor

          if ( mod ( power(i), 2 ) .eq. 1 ) then

            nmore = nmore + 1
            nstack = nstack + 1

            if ( maxstack .lt. nstack ) then
              write ( *, '(a)' ) ' '
              write ( *, '(a)' ) 'LEGENDRE_SYMBOL - Fatal error!'
              write ( *, '(a)' ) '  Stack overflowc'
              l = - 4
              return
            end if

            pstack(nstack) = pp
            qstack(nstack) = factor(i)

          end if

        end do

        if ( nmore .ne. 0 ) then

          qq = qstack(nstack)
          nstack = nstack - 1
c
c  Check for a QQ of 1 or 2.
c
          if ( qq .eq. 1 ) then

            l = + 1 * l

          else if ( qq .eq. 2 .and. 
     &            ( mod ( pp, 8 ) .eq. 1 .or. 
     &              mod ( pp, 8 ) .eq. 7 ) ) then

            l = + 1 * l

          else if ( qq .eq. 2 .and. 
     &            ( mod ( pp, 8 ) .eq. 3 .or. 
     &              mod ( pp, 8 ) .eq. 5 ) ) then

            l = - 1 * l

          else

            if ( mod ( pp, 4 ) .eq. 3 .and. 
     &           mod ( qq, 4 ) .eq. 3 ) then
              l = - 1 * l
            end if

            t  = pp
            pp = qq
            qq = t

            go to 20

          end if

        end if
c
c  If the stack is empty, we're done.
c
        if ( nstack .eq. 0 ) then
          go to 30
        end if
c
c  Otherwise, get the last P and Q from the stack, and process them.
c
        pp = pstack(nstack)
        qq = qstack(nstack)
        nstack = nstack - 1

      go to 20

30    continue

      return
      end
      function lerch ( z, s, a )

c*********************************************************************72
c
cc LERCH estimates the Lerch transcendent function.
c
c  Discussion:
c
c    The Lerch transcendent function is defined as:
c
c      LERCH ( Z, S, A ) = Sum ( 0 <= K < +oo ) Z**K / ( A + K )**S
c
c    excluding any term with ( A + K ) = 0.
c
c    In Mathematica, the function can be evaluated by:
c
c      LerchPhi[z,s,a]
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    11 July 2008
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Eric Weisstein,
c    CRC Concise Encyclopedia of Mathematics,
c    CRC Press, 2002,
c    Second edition,
c    ISBN: 1584883472,
c    LC: QA5.W45
c
c  Thanks:
c
c    Oscar van Vlijmen
c
c  Parameters:
c
c    Input, double precision Z, integer S, double precision A, 
c    the parameters of the function.
c
c    Output, double precision LERCH, an approximation to the Lerch
c    transcendent function.
c
      implicit none

      double precision a
      double precision eps
      integer k
      double precision lerch
      integer s
      double precision term
      double precision total
      double precision z
      double precision z_k

      if ( z .le. 0.0D+00 ) then
        lerch = 0.0D+00
        return
      end if

      eps = 1.0D-10
      total = 0.0D+00
      k = 0
      z_k = 1.0D+00

10    continue

        if ( a + dble ( k ) .ne. 0.0D+00 ) then

          term = z_k / ( a + dble ( k ) )**s
          total = total + term

          if ( abs ( term ) <= eps * ( 1.0D+00 + abs ( total ) ) ) then
            go to 20
          end if

        end if

        k = k + 1
        z_k = z_k * z

        go to 10

20    continue

      lerch = total

      return
      end
      subroutine lerch_values ( n_data, z, s, a, fx )

c*********************************************************************72
c
cc LERCH_VALUES returns some values of the Lerch transcendent function.
c
c  Discussion:
c
c    The Lerch function is defined as
c
c      Phi(z,s,a) = Sum ( 0 <= k .lt. +oo ) z^k / ( a + k )^s
c
c    omitting any terms with ( a + k ) = 0.
c
c    In Mathematica, the function can be evaluated by:
c
c      LerchPhi[z,s,a]
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    24 March 2007
c
c  Author:
c
c    John Burkardt
c
c  Reference:
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
c    Output, double precision Z, the parameters of the function.
c
c    Output, integer S, the parameters of the function.
c
c    Output, double precision A, the parameters of the function.
c
c    Output, double precision FX, the value of the function.
c
      implicit none

      integer n_max
      parameter ( n_max = 12 )

      double precision a
      double precision a_vec(n_max)
      double precision fx
      double precision fx_vec(n_max)
      integer n_data
      integer s
      integer s_vec(n_max)
      double precision z
      double precision z_vec(n_max)

      save a_vec
      save fx_vec
      save s_vec
      save z_vec

      data a_vec /
     &  0.0D+00,
     &  0.0D+00,
     &  0.0D+00,
     &  1.0D+00,
     &  1.0D+00,
     &  1.0D+00,
     &  2.0D+00,
     &  2.0D+00,
     &  2.0D+00,
     &  3.0D+00,
     &  3.0D+00,
     &  3.0D+00 /
      data fx_vec /
     &  0.1644934066848226D+01,
     &  0.1202056903159594D+01,
     &  0.1000994575127818D+01,
     &  0.1164481052930025D+01,
     &  0.1074426387216080D+01,
     &  0.1000492641212014D+01,
     &  0.2959190697935714D+00,
     &  0.1394507503935608D+00,
     &  0.9823175058446061D-03,
     &  0.1177910993911311D+00,
     &  0.3868447922298962D-01,
     &  0.1703149614186634D-04 /
      data s_vec /
     &   2, 3, 10,
     &   2, 3, 10,
     &   2, 3, 10,
     &   2, 3, 10 /
      data z_vec /
     &  0.1000000000000000D+01,
     &  0.1000000000000000D+01,
     &  0.1000000000000000D+01,
     &  0.5000000000000000D+00,
     &  0.5000000000000000D+00,
     &  0.5000000000000000D+00,
     &  0.3333333333333333D+00,
     &  0.3333333333333333D+00,
     &  0.3333333333333333D+00,
     &  0.1000000000000000D+00,
     &  0.1000000000000000D+00,
     &  0.1000000000000000D+00 /

      if ( n_data .lt. 0 ) then
        n_data = 0
      end if

      n_data = n_data + 1

      if ( n_max .lt. n_data ) then
        n_data = 0
        z = 0.0D+00
        s = 0
        a = 0.0D+00
        fx = 0.0D+00
      else
        z = z_vec(n_data)
        s = s_vec(n_data)
        a = a_vec(n_data)
        fx = fx_vec(n_data)
      end if

      return
      end
      subroutine lock ( n, a )

c*********************************************************************72
c
cc LOCK returns the number of codes for a lock with N buttons.
c
c  Discussion:
c
c    A button lock has N numbered buttons.  To open the lock, groups
c    of buttons must be pressed in the correct order.  Each button
c    may be pushed no more than once.  Thus, a code for the lock is 
c    an ordered list of the groups of buttons to be pushed.
c
c    For this discussion, we will assume that EVERY button is pushed
c    at some time, as part of the code.  To count the total number
c    of codes, including those which don't use all the buttons, then
c    the number is 2 * A(N), or 2 * A(N) - 1 if we don't consider the
c    empty code to be valid.
c
c  Examples:
c
c    If there are 3 buttons, then there are 13 possible "full button" codes:
c
c      (123)
c      (12) (3)
c      (13) (2)
c      (23) (1)
c      (1) (23)
c      (2) (13)
c      (3) (12)
c      (1) (2) (3)
c      (1) (3) (2)
c      (2) (1) (3)
c      (2) (3) (1)
c      (3) (1) (2)
c      (3) (2) (1)
c
c    and, if we don't need to push all the buttons, every "full button" code above
c    yields a distinct "partial button" code by dropping the last set of buttons:
c
c      ()
c      (12)
c      (13)
c      (23)
c      (1)
c      (2)
c      (3)
c      (1) (2)
c      (1) (3)
c      (2) (1)
c      (2) (3)
c      (3) (1)
c      (3) (2)
c
c  First values:
c
c     N         A(N)
c     0           1
c     1           1
c     2           3
c     3          13
c     4          75
c     5         541
c     6        4683
c     7       47293
c     8      545835
c     9     7087261
c    10   102247563
c
c  Recursion:
c
c    A(I) = sum ( 0 <= J < I ) Binomial ( I, N-J ) * A(J)
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    11 July 2008
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Daniel Velleman, Gregory Call,
c    Permutations and Combination Locks,
c    Mathematics Magazine,
c    Volume 68, Number 4, October 1995, pages 243-253.
c
c  Parameters:
c
c    Input, integer N, the maximum number of lock buttons.
c
c    Output, integer A(0:N), the number of lock codes.
c
      implicit none

      integer n

      integer a(0:n)
      integer i
      integer i4_choose
      integer j

      if ( n .lt. 0 ) then
        return
      end if

      a(0) = 1

      do i = 1, n
        a(i) = 0
        do j = 0, i - 1
          a(i) = a(i) + i4_choose ( i, i - j ) * a(j)
        end do
      end do

      return
      end
      subroutine meixner ( n, beta, c, x, v )

c*********************************************************************72
c
cc MEIXNER evaluates Meixner polynomials at a point.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    18 March 2009
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Walter Gautschi,
c    Orthogonal Polynomials: Computation and Approximation,
c    Oxford, 2004,
c    ISBN: 0-19-850672-4,
c    LC: QA404.5 G3555.
c
c  Parameters:
c
c    Input, integer N, the maximum order of the polynomial.  
c    N must be at least 0.
c
c    Input, double precision BETA, the Beta parameter.  0 < BETA.
c
c    Input, double precision C, the C parameter.  0 < C < 1.
c
c    Input, double precision X, the evaluation point.
c
c    Output, double precision V(0:N), the value of the polynomials at X.
c
      implicit none

      integer n

      double precision beta
      double precision c
      integer i
      double precision v(0:n)
      double precision x

      if ( beta .le. 0.0D+00 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'MEIXNER - Fatal error!'
        write ( *, '(a)' ) '  Parameter BETA must be positive.'
        stop
      end if

      if ( c .le. 0.0D+00 .or. 1.0D+00 .le. c ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'MEIXNER - Fatal error!'
        write ( *, '(a)' ) 
     &  '  Parameter C must be strictly between 0 and 1.'
        stop
      end if

      if ( n .lt. 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'MEIXNER - Fatal error!'
        write ( *, '(a)' ) '  Parameter N must be nonnegative.'
        stop
      end if

      v(0) = 1.0D+00

      if ( n .eq. 0 ) then
        return
      end if

      v(1) = ( c - 1.0D+00 ) * x / beta / c + 1.0D+00

      if ( n == 1 ) then
        return
      end if

      do i = 1, n - 1
        v(i+1) = ( 
     &    ( ( c - 1.0D+00 ) * x + ( 1.0D+00 + c ) 
     &    * dble ( i ) + beta * c ) * v(i) 
     &    - dble ( i ) * v(i-1) 
     &    ) / ( dble ( i ) + beta )
      end do

      return
      end
      function mertens ( n )

c*********************************************************************72
c
cc MERTENS evaluates the Mertens function.
c
c  Discussion:
c
c    The Mertens function M(N) is the sum from 1 to N of the Moebius
c    function MU.  That is,
c
c    M(N) = sum ( 1 <= I <= N ) MU(I)
c
c        N   M(N)
c        --  ----
c         1     1
c         2     0
c         3    -1
c         4    -1
c         5    -2
c         6    -1
c         7    -2
c         8    -2
c         9    -2
c        10    -1
c        11    -2
c        12    -2
c       100     1
c      1000     2
c     10000   -23
c    100000   -48
c
c    The determinant of the Redheffer matrix of order N is equal
c    to the Mertens function M(N).
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    05 July 2008
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    M Deleglise, J Rivat,
c    Computing the Summation of the Moebius Function,
c    Experimental Mathematics,
c    Volume 5, 1996, pages 291-295.
c
c    Eric Weisstein,
c    CRC Concise Encyclopedia of Mathematics,
c    CRC Press, 2002,
c    Second edition,
c    ISBN: 1584883472,
c    LC: QA5.W45
c
c  Parameters:
c
c    Input, integer N, the argument.
c
c    Output, integer MERTENS, the value.
c
      implicit none

      integer i
      integer mertens
      integer mu_i
      integer n
      integer value

      value = 0

      do i = 1, n
        call moebius ( i, mu_i )
        value = value + mu_i
      end do

      mertens = value

      return
      end
      subroutine mertens_values ( n_data, n, c )

c*********************************************************************72
c
cc MERTENS_VALUES returns some values of the Mertens function.
c
c  Discussion:
c
c    The Mertens function M(N) is the sum from 1 to N of the Moebius
c    function MU.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    01 Decemberr 2007
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Marc Deleglise, Joel Rivat,
c    Computing the Summation of the Moebius Function,
c    Experimental Mathematics,
c    Volume 5, 1996, pages 291-295.
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
c    Input/output, integer N_DATA.
c    On input, if N_DATA is 0, the first test data is returned, and N_DATA
c    is set to 1.  On each subsequent call, the input value of N_DATA is
c    incremented and that test data item is returned, if available.  When
c    there is no more test data, N_DATA is set to 0.
c
c    Output, integer N, the argument of the Mertens function.
c
c    Output, integer C, the value of the Mertens function.
c
      implicit none

      integer nmax
      parameter ( nmax = 15 )

      integer c
      integer c_vec(nmax)
      integer n
      integer n_data
      integer n_vec(nmax)


      save c_vec
      save n_vec

      data c_vec /
     &    1,   0,  -1,   -1,  -2,  -1,  -2,  -2,   -2,  -1, 
     &   -2,  -2,   1,    2, -23 /
      data n_vec /
     &    1,   2,   3,   4,   5,   6,   7,   8,   9,  10, 
     &   11,  12,  100, 1000, 10000 /
      
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
      subroutine moebius ( n, mu )

c*********************************************************************72
c
cc MOEBIUS returns the value of MU(N), the Moebius function of N.
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
        write ( *, '(a)' ) 'MOEBIUS - Fatal error!'
        write ( *, '(a)' ) '  Not enough factorization space.'
        mu = -3
        return
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
      subroutine moebius_values ( n_data, n, c )

c*********************************************************************72
c
cc MOEBIUS_VALUES returns some values of the Moebius function.
c
c  Discussion:
c
c    MU(N) is defined as follows:
c
c      MU(N) = 1 if N = 1;
c              0 if N is divisible by the square of a prime;
c              (-1)**K, if N is the product of K distinct primes.
c
c    In Mathematica, the function can be evaluated by:
c
c      MoebiusMu[n]
c
c    The Moebius function is related to Euler's totient function:
c
c      PHI(N) = Sum ( D divides N ) MU(D) * ( N / D ).
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
c    Note that, as special cases, MU(N) is -1 if N is a prime, and MU(N) is 0
c    if N is a square, cube, etc.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    24 March 2007
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
c    Output, integer N, the argument of the Moebius function.
c
c    Output, integer C, the value of the Moebius function.
c
      implicit none

      integer n_max
      parameter ( n_max = 20 )

      integer c
      integer c_vec(n_max)
      integer n
      integer n_data
      integer n_vec(n_max)

      save c_vec
      save n_vec

      data c_vec /
     &    1,  -1,  -1,   0,  -1,   1,  -1,   0,   0,   1,
     &   -1,   0,  -1,   1,   1,   0,  -1,   0,  -1,   0 /
      data n_vec /
     &    1,   2,   3,   4,   5,   6,   7,   8,   9,  10,
     &   11,  12,  13,  14,  15,  16,  17,  18,  19,  20 /

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
      subroutine motzkin ( n, a )

c*********************************************************************72
c
cc MOTZKIN returns the Motzkin numbers up to order N.
c
c  Discussion:
c
c    The Motzkin number A(N) counts the number of distinct paths
c    from (0,0) to (0,N) in which the only steps used are
c    (1,1), (1,-1) and (1,0), and the path is never allowed to
c    go below the X axis.
c
c  First values:
c
c     N  A(N)
c
c     0    1
c     1    1
c     2    2
c     3    4
c     4    9
c     5   21
c     6   51
c     7  127
c     8  323
c     9  835
c    10 2188
c
c  Recursion:
c
c    A(N) = A(N-1) + sum ( 0 <= K <= N-2 ) A(K) * A(N-2-K)
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    05 July 2008
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Eric Weisstein,
c    CRC Concise Encyclopedia of Mathematics,
c    CRC Press, 2002,
c    Second edition,
c    ISBN: 1584883472,
c    LC: QA5.W45
c
c  Parameters:
c
c    Input, integer N, the highest order Motzkin number to compute.
c
c    Output, integer A(0:N), the Motzkin numbers.
c
      implicit none

      integer n

      integer a(0:n)
      integer i
      integer j

      if ( n .lt. 0 ) then
        return
      end if

      a(0) = 1

      do i = 1, n
        a(i) = a(i-1)
        do j = 0, i - 2
          a(i) = a(i) + a(j) * a(i-2-j)
        end do
      end do

      return
      end
      function normal_01_cdf_inv ( p )

c*********************************************************************72
c
cc NORMAL_01_CDF_INV inverts the standard normal CDF.
c
c  Discussion:
c
c    The result is accurate to about 1 part in 10**16.
c
c  Modified:
c
c    13 January 2008
c
c  Author:
c
c    Michael Wichura
c
c  Reference:
c
c    Michael Wichura,
c    Algorithm AS 241:
c    The Percentage Points of the Normal Distribution,
c    Applied Statistics,
c    Volume 37, Number 3, 1988, pages 477-484.
c
c  Parameters:
c
c    Input, double precision P, the value of the cumulative probability 
c    densitity function.  0 < P < 1.
c
c    Output, integer IFAULT, error flag.
c    0, no error.
c    1, P <= 0 or P >= 1.
c
c    Output, double precision PPND16, the normal deviate value with the 
c    property that the probability of a standard normal deviate being 
c    less than or equal to PPND16 is P.
c
      implicit none

      double precision a0
      double precision a1
      double precision a2
      double precision a3
      double precision a4
      double precision a5
      double precision a6
      double precision a7
      double precision b1
      double precision b2
      double precision b3
      double precision b4
      double precision b5
      double precision b6
      double precision b7
      double precision c0
      double precision c1
      double precision c2
      double precision c3
      double precision c4
      double precision c5
      double precision c6
      double precision c7
      double precision const1
      double precision const2
      double precision d1
      double precision d2
      double precision d3
      double precision d4
      double precision d5
      double precision d6
      double precision d7
      double precision e0
      double precision e1
      double precision e2
      double precision e3
      double precision e4
      double precision e5
      double precision e6
      double precision e7
      double precision f1
      double precision f2
      double precision f3
      double precision f4
      double precision f5
      double precision f6
      double precision f7
      double precision normal_01_cdf_inv
      double precision p
      double precision q
      double precision r
      double precision split1
      double precision split2

      parameter ( a0 = 3.3871328727963666080D+00 )
      parameter ( a1 = 1.3314166789178437745D+02 )
      parameter ( a2 = 1.9715909503065514427D+03 )
      parameter ( a3 = 1.3731693765509461125D+04 )
      parameter ( a4 = 4.5921953931549871457D+04 )
      parameter ( a5 = 6.7265770927008700853D+04 )
      parameter ( a6 = 3.3430575583588128105D+04 )
      parameter ( a7 = 2.5090809287301226727D+03 )
      parameter ( b1 = 4.2313330701600911252D+01 )
      parameter ( b2 = 6.8718700749205790830D+02 )
      parameter ( b3 = 5.3941960214247511077D+03 )
      parameter ( b4 = 2.1213794301586595867D+04 )
      parameter ( b5 = 3.9307895800092710610D+04 )
      parameter ( b6 = 2.8729085735721942674D+04 )
      parameter ( b7 = 5.2264952788528545610D+03 )
      parameter ( c0 = 1.42343711074968357734D+00 )
      parameter ( c1 = 4.63033784615654529590D+00 )
      parameter ( c2 = 5.76949722146069140550D+00 )
      parameter ( c3 = 3.64784832476320460504D+00 )
      parameter ( c4 = 1.27045825245236838258D+00 )
      parameter ( c5 = 2.41780725177450611770D-01 )
      parameter ( c6 = 2.27238449892691845833D-02 )
      parameter ( c7 = 7.74545014278341407640D-04 )
      parameter ( const1 = 0.180625D+00 )
      parameter ( const2 = 1.6D+00 )
      parameter ( d1 = 2.05319162663775882187D+00 )
      parameter ( d2 = 1.67638483018380384940D+00 )
      parameter ( d3 = 6.89767334985100004550D-01 )
      parameter ( d4 = 1.48103976427480074590D-01 )
      parameter ( d5 = 1.51986665636164571966D-02 )
      parameter ( d6 = 5.47593808499534494600D-04 )
      parameter ( d7 = 1.05075007164441684324D-09 )
      parameter ( e0 = 6.65790464350110377720D+00 )
      parameter ( e1 = 5.46378491116411436990D+00 )
      parameter ( e2 = 1.78482653991729133580D+00 )
      parameter ( e3 = 2.96560571828504891230D-01 )
      parameter ( e4 = 2.65321895265761230930D-02 )
      parameter ( e5 = 1.24266094738807843860D-03 )
      parameter ( e6 = 2.71155556874348757815D-05 )
      parameter ( e7 = 2.01033439929228813265D-07 )
      parameter ( f1 = 5.99832206555887937690D-01 )
      parameter ( f2 = 1.36929880922735805310D-01 )
      parameter ( f3 = 1.48753612908506148525D-02 )
      parameter ( f4 = 7.86869131145613259100D-04 )
      parameter ( f5 = 1.84631831751005468180D-05 )
      parameter ( f6 = 1.42151175831644588870D-07 )
      parameter ( f7 = 2.04426310338993978564D-15 )
      parameter ( split1 = 0.425D+00 )
      parameter ( split2 = 5.D+00 )

      q = p - 0.5D+00

      if ( dabs ( q ) .le. split1 ) then

        r = const1 - q * q

        normal_01_cdf_inv = q * (((((((
     &      a7   * r 
     &    + a6 ) * r 
     &    + a5 ) * r 
     &    + a4 ) * r 
     &    + a3 ) * r 
     &    + a2 ) * r 
     &    + a1 ) * r 
     &    + a0 ) / (((((((
     &      b7   * r 
     &    + b6 ) * r 
     &    + b5 ) * r 
     &    + b4 ) * r 
     &    + b3 ) * r 
     &    + b2 ) * r 
     &    + b1 ) * r 
     &    + 1.0D+00 )

      else

        if ( q .lt. 0.0D+00 ) then
          r = p
        else
          r = 1.0D+00 - p
        end if

        if ( r .le. 0.0D+00 ) then
          normal_01_cdf_inv = 0.0D+00
          return
        end if

        r = dsqrt ( - dlog ( r ) )

        if ( r .le. split2 ) then

          r = r - const2

          normal_01_cdf_inv = (((((((
     &        c7   * r 
     &      + c6 ) * r 
     &      + c5 ) * r 
     &      + c4 ) * r 
     &      + c3 ) * r 
     &      + c2 ) * r 
     &      + c1 ) * r 
     &      + c0 ) / (((((((
     &        d7   * r 
     &      + d6 ) * r 
     &      + d5 ) * r 
     &      + d4 ) * r 
     &      + d3 ) * r 
     &      + d2 ) * r 
     &      + d1 ) * r 
     &      + 1.0D+00 )

        else

          r = r - split2

          normal_01_cdf_inv = (((((((
     &        e7   * r 
     &      + e6 ) * r 
     &      + e5 ) * r 
     &      + e4 ) * r 
     &      + e3 ) * r 
     &      + e2 ) * r 
     &      + e1 ) * r 
     &      + e0 ) / (((((((
     &        f7   * r 
     &      + f6 ) * r 
     &      + f5 ) * r 
     &      + f4 ) * r 
     &      + f3 ) * r 
     &      + f2 ) * r 
     &      + f1 ) * r 
     &      + 1.0D+00 )

        end if

        if ( q .lt. 0.0D+00 ) then
          normal_01_cdf_inv = - normal_01_cdf_inv
        end if

      end if

      return
      end
      subroutine omega ( n, ndiv )

c*********************************************************************72
c
cc OMEGA returns OMEGA(N), the number of distinct prime divisors of N.
c
c  Discussion:
c
c    If N = 1, then
c
c      OMEGA(N) = 1
c
c    else if the prime factorization of N is
c
c      N = P1**E1 * P2**E2 * ... * PM**EM,
c
c    then
c
c      OMEGA(N) = M
c
c  Example:
c
c     N   OMEGA(N)
c
c     1    1
c     2    1
c     3    1
c     4    1
c     5    1
c     6    2
c     7    1
c     8    1
c     9    1
c    10    2
c    11    1
c    12    2
c    13    1
c    14    2
c    15    2
c    16    1
c    17    1
c    18    2
c    19    1
c    20    2
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    05 July 2008
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the value to be analyzed.  N must be 1 or
c    greater.
c
c    Output, integer NDIV, the value of OMEGA(N).  But if N is 0 or
c    less, NDIV is returned as 0, a nonsense value.  If there is
c    not enough room for factoring, NDIV is returned as -1.
c
      implicit none

      integer maxfactor
      parameter ( maxfactor = 20 )

      integer factor(maxfactor)
      integer n
      integer ndiv
      integer nfactor
      integer nleft
      integer power(maxfactor)

      if ( n .le. 0 ) then
        ndiv = 0
        return
      end if

      if ( n .eq. 1 ) then
        ndiv = 1
        return
      end if
c
c  Factor N.
c
      call i4_factor ( n, maxfactor, nfactor, factor, power, nleft )

      if ( nleft .ne. 1 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'OMEGA - Fatal error!'
        write ( *, '(a)' ) '  Not enough factorization space.'
        ndiv = -1
        return
      end if

      ndiv = nfactor

      return
      end
      subroutine omega_values ( n_data, n, c )

c*********************************************************************72
c
cc OMEGA_VALUES returns some values of the OMEGA function.
c
c  Discussion:
c
c    In Mathematica, the function can be evaluated by
c
c      Length [ FactorInteger [ n ] ]
c
c    If N = 1, then
c
c      OMEGA(N) = 1
c
c    else if the prime factorization of N is
c
c      N = P1**E1 * P2**E2 * ... * PM**EM,
c
c    then
c
c      OMEGA(N) = M
c
c  Example:
c
c     N   OMEGA(N)
c
c     1    1
c     2    1
c     3    1
c     4    1
c     5    1
c     6    2
c     7    1
c     8    1
c     9    1
c    10    2
c    11    1
c    12    2
c    13    1
c    14    2
c    15    2
c    16    1
c    17    1
c    18    2
c    19    1
c    20    2
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    25 March 2007
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
c    Output, integer N, the argument of the OMEGA function.
c
c    Output, integer C, the value of the OMEGA function.
c
      implicit none

      integer n_max
      parameter ( n_max = 23 )

      integer c
      integer c_vec(n_max)
      integer n
      integer n_data
      integer n_vec(n_max)

      save c_vec
      save n_vec

      data c_vec /
     &    1,   1,   1,   1,   1,
     &    2,   1,   1,   1,   2,
     &    3,   1,   4,   4,   3,
     &    1,   5,   2,   2,   1,
     &    6,   7,   8 /
      data n_vec /
     &         1,
     &         2,
     &         3,
     &         4,
     &         5,
     &         6,
     &         7,
     &         8,
     &         9,
     &        10,
     &        30,
     &       101,
     &       210,
     &      1320,
     &      1764,
     &      2003,
     &      2310,
     &      2827,
     &      8717,
     &     12553,
     &     30030,
     &    510510,
     &   9699690 /

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
      subroutine partition_count_values ( n_data, n, c )

c*********************************************************************72
c
cc PARTITION_COUNT_VALUES returns some values of the integer partition count.
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
c    National Bureau of Standards, 1964,
c    ISBN: 0-486-61272-4,
c    LC: QA47.A34.
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
      subroutine partition_distinct_count_values ( n_data, n, c )

c*********************************************************************72
c
cc PARTITION_DISTINCT_COUNT_VALUES returns some values of Q(N).
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
c    However, if we require that each member of the partition
c    be distinct, so that no nonzero summand occurs more than once,
c    we are computing something symbolized by Q(N).
c    The number 5 has Q(N) = 3, because it has the following partitions
c    into distinct parts:
c
c    5 = 5
c      = 4 + 1
c      = 3 + 2
c
c    In Mathematica, the function can be evaluated by
c
c      PartitionsQ[n]
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    25 March 2007
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
c    Output, integer N, the integer.
c
c    Output, integer C, the number of partitions of the integer
c    into distinct parts.
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
     &    1,   1,   2,   2,   3,   4,   5,   6,   8,  10,
     &   12,  15,  18,  22,  27,  32,  38,  46,  54,  64 /
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
      subroutine pentagon_num ( n, p )

c*********************************************************************72
c
cc PENTAGON_NUM computes the N-th pentagonal number.
c
c  Discussion:
c
c    The pentagonal number P(N) counts the number of dots in a figure of
c    N nested pentagons.  The pentagonal numbers are defined for both
c    positive and negative N.
c
c    The formula is:
c
c      P(N) = ( N * ( 3 * N - 1 ) ) / 2
c
c  Example:
c
c    N   P
c
c   -5  40
c   -4  26
c   -3  15
c   -2   7
c   -1   2
c    0   0
c    1   1
c    2   5
c    3  12
c    4  22
c    5  35
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    03 July 2008
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
      subroutine phi ( n, phin )

c*********************************************************************72
c
cc PHI computes the number of relatively prime predecessors of an integer.
c
c  Discussion:
c
c    PHI(N) is the number of integers between 1 and N which are
c    relatively prime to N.  I and J are relatively prime if they
c    have no common factors.  The function PHI(N) is known as
c    "Euler's totient function".
c
c    By convention, 1 and N are relatively prime.
c
c    The formula is:
c
c      PHI(U*V) = PHI(U) * PHI(V) if U and V are relatively prime.
c
c      PHI(P**K) = P**(K-1) * ( P - 1 ) if P is prime.
c
c      PHI(N) = N * Product ( P divides N ) ( 1 - 1 / P )
c
c      N = Sum ( D divides N ) PHI(D).
c
c  Example:
c
c     N  PHI(N)
c
c     1    1
c     2    1
c     3    2
c     4    2
c     5    4
c     6    2
c     7    6
c     8    4
c     9    6
c    10    4
c    11   10
c    12    4
c    13   12
c    14    6
c    15    8
c    16    8
c    17   16
c    18    6
c    19   18
c    20    8
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    11 July 2008
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the value to be analyzed.
c
c    Output, integer PHIN, the value of PHI(N).  If N is less than
c    or equal to 0, PHI will be returned as 0.  If there is not enough
c    room for full factoring of N, PHI will be returned as -1.
c
      implicit none

      integer maxfactor
      parameter ( maxfactor = 20 )

      integer factor(maxfactor)
      integer i
      integer n
      integer nfactor
      integer nleft
      integer phin
      integer power(maxfactor)

      if ( n .le. 0 ) then
        phin = 0
        return
      end if

      if ( n .eq. 1 ) then
        phin = 1
        return
      end if
c
c  Factor N.
c
      call i4_factor ( n, maxfactor, nfactor, factor, power, nleft )

      if ( nleft .ne. 1 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'PHI - Fatal error!'
        write ( *, '(a)' ) '  Not enough factorization space!'
        phin = -1
        return
      end if

      phin = 1
      do i = 1, nfactor
        phin = phin * factor(i)**( power(i) - 1 ) * ( factor(i) - 1 )
      end do

      return
      end
      subroutine phi_values ( n_data, n, c )

c*********************************************************************72
c
cc PHI_VALUES returns some values of the PHI function.
c
c  Discussion:
c
c    PHI(N) is the number of integers between 1 and N which are
c    relatively prime to N.  I and J are relatively prime if they
c    have no common factors.  The function PHI(N) is known as
c    "Euler's totient function".
c
c    By convention, 1 and N are relatively prime.
c
c    In Mathematica, the function can be evaluated by:
c
c      EulerPhi[n]
c
c    The formula is:
c
c      PHI(U*V) = PHI(U) * PHI(V) if U and V are relatively prime.
c
c      PHI(P**K) = P**(K-1) * ( P - 1 ) if P is prime.
c
c      PHI(N) = N * Product ( P divides N ) ( 1 - 1 / P )
c
c      N = Sum ( D divides N ) PHI(D).
c
c  Example:
c
c     N  PHI(N)
c
c     1    1
c     2    1
c     3    2
c     4    2
c     5    4
c     6    2
c     7    6
c     8    4
c     9    6
c    10    4
c    11   10
c    12    4
c    13   12
c    14    6
c    15    8
c    16    8
c    17   16
c    18    6
c    19   18
c    20    8
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    25 March 2007
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
c    Output, integer N, the argument of the PHI function.
c
c    Output, integer C, the value of the PHI function.
c
      implicit none

      integer n_max
      parameter ( n_max = 20 )

      integer c
      integer c_vec(n_max)
      integer n
      integer n_data
      integer n_vec(n_max)

      save c_vec
      save n_vec

      data c_vec /
     &    1,   1,   2,   2,   4,   2,   6,   4,   6,   4,
     &    8,   8,  16,  20,  16,  40, 148, 200, 200, 648 /
      data n_vec /
     &    1,   2,   3,   4,   5,   6,   7,   8,   9,  10,
     &   20,  30,  40,  50,  60, 100, 149, 500, 750, 999 /

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
      subroutine poly_bernoulli ( n, k, b )

c*********************************************************************72
c
cc POLY_BERNOULLI evaluates the poly-Bernolli numbers with negative index.
c
c  Discussion:
c
c    The poly-Bernoulli numbers B_n^k were defined by M Kaneko
c    formally as the coefficients of X^n/nc in a particular power 
c    series.  He also showed that, when the super-index is negative,
c    we have
c
c      B_n^(-k) = Sum ( 0 <= j <= min ( n, k ) ) 
c        (jc)^2 * S(n+1,j+1) * S(k+1,j+1)
c
c    where S(n,k) is the Stirling number of the second kind, the number of
c    ways to partition a set of size n into k nonempty subset.
c
c    B_n^(-k) is also the number of "lonesum matrices", that is, 0-1
c    matrices of n rows and k columns which are uniquely reconstructable
c    from their row and column sums.
c
c    The poly-Bernoulli numbers get large very quickly.
c
c  Table:
c
c    \ K 0  1    2     3      4       5        6
c    N
c    0   1  1    1     1      1       1        1
c    1   1  2    4     8     16      32       64
c    2   1  4   14    46    146     454     1394
c    3   1  8   46   230   1066    4718    20266
c    4   1 16  146  1066   6902   41506   237686
c    5   1 32  454  4718  41506  329462  2441314
c    6   1 64 1394 20266 237686 2441314 22934774
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    18 July 2008
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Chad Brewbaker,
c    Lonesum (0,1) Matrices and Poly-Bernoulli Numbers of Negative Index,
c    MS Thesis,
c    Iowa State University, 2005.
c
c    M Kaneko,
c    Poly-Bernoulli Numbers,
c    Journal Theorie des Nombres Bordeaux,
c    Volume 9, 1997, pages 221-228.
c
c  Parameters:
c
c    Input, integer N, K, the indices.  N and K should be 
c    nonnegative.
c
c    Output, integer B, the value of B_N^(-K).
c
      implicit none

      integer m_max
      parameter ( m_max = 20 )

      integer b
      integer j
      integer jfact
      integer jhi
      integer k
      integer m
      integer n
      integer s(m_max*m_max)

      if ( n .lt. 0 ) then
        b = 0
        return
      else if ( n .eq. 0 ) then
        b = 1
        return
      end if

      if ( k .lt. 0 ) then
        b = 0
        return
      else if ( k .eq. 0 ) then
        b = 1
        return
      end if

      jhi = min ( n, k )
      m = max ( n, k ) + 1

      if ( m_max < m ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'POLY_BERNOULLI - Fatal error!'
        write ( *, '(a)' ) '  Internal storage M_MAX = ', m_max
        write ( *, '(a)' ) '  exceeded by value M = ', m
        stop
      end if

      call stirling2 ( m, m, s )

      jfact = 1
      b = 0

      do j = 0, jhi

        b = b + jfact * jfact * s(n+1+j*m_max) * s(k+1+j*m_max)

        jfact = jfact * ( j + 1 )

      end do

      return
      end
      function poly_coef_count ( dim, degree )

c*********************************************************************72
c
cc POLY_COEF_COUNT: polynomial coefficient count given dimension and degree.
c
c  Discussion:
c
c    To count all monomials of degree 5 or less in dimension 3,
c    we can count all monomials of degree 5 in dimension 4.
c
c    To count all monomials of degree 5 in dimension 4, we imagine
c    that each of the variables X, Y, Z and W is a "box" and that
c    we need to drop 5 pebbles into these boxes.  Every distinct
c    way of doing this represents a degree 5 monomial in dimension 4.
c    Ignoring W gives us monomials up to degree five in dimension 3.
c
c    To count them, we draw 3 lines as separators to indicate the
c    4 boxes, and then imagine all distinct sequences involving
c    the three lines and the 5 pebbles.  Indicate the lines by 1's
c    and the pebbles by 0's and we're asking for the number of
c    permutations of 3 1's and 5 0's, which is 8! / (3! 5!)
c
c    In other words, 56 = 8! / (3! 5!) is:
c    * the number of monomials of degree exactly 5 in dimension 4, 
c    * the number of monomials of degree 5 or less in dimension 3, 
c    * the number of polynomial coefficients of a polynomial of 
c      degree 5 in (X,Y,Z).
c
c    In general, the formula for the number of monomials of degree DEG
c    or less in dimension DIM is
c
c      (DEG+DIM)! / (DEG! * DIM!)
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
c    Input, integer DIM, the dimension of the polynomial.
c    0 <= DIM.
c
c    Input, integer DEGREE, the degree of the polynomnial
c    0 <= DEGREE
c
c    Output, integer POLY_COEF_COUNT, the number of coefficients 
c    in the general polynomial of dimension DIM and degree DEGREE.
c
      implicit none

      integer degree
      integer dim
      integer i4_choose
      integer poly_coef_count

      if ( dim .lt. 0 ) then
        poly_coef_count = -1
      else if ( degree .lt. 0 ) then
        poly_coef_count = -1
      else
        poly_coef_count = i4_choose ( degree + dim, degree )
      end if

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
      subroutine psi_values ( n_data, x, fx )

c*********************************************************************72
c
cc PSI_VALUES returns some values of the Psi or Digamma function for testing.
c
c  Discussion:
c
c    PSI(X) = d LN ( GAMMA ( X ) ) / d X = GAMMA'(X) / GAMMA(X)
c
c    PSI(1) = - Euler's constant.
c
c    PSI(X+1) = PSI(X) + 1 / X.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    31 March 2007
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
      parameter ( n_max = 11 )

      double precision fx
      double precision fxvec ( n_max )
      integer n_data
      double precision x
      double precision xvec ( n_max )

      data fxvec /
     &  -0.5772156649015329D+00,
     &  -0.4237549404110768D+00,
     &  -0.2890398965921883D+00,
     &  -0.1691908888667997D+00,
     &  -0.6138454458511615D-01,
     &   0.3648997397857652D-01,
     &   0.1260474527734763D+00,
     &   0.2085478748734940D+00,
     &   0.2849914332938615D+00,
     &   0.3561841611640597D+00,
     &   0.4227843350984671D+00 /

      data xvec /
     &  1.0D+00,
     &  1.1D+00,
     &  1.2D+00,
     &  1.3D+00,
     &  1.4D+00,
     &  1.5D+00,
     &  1.6D+00,
     &  1.7D+00,
     &  1.8D+00,
     &  1.9D+00,
     &  2.0D+00 /

      if ( n_data .lt. 0 ) then
        n_data = 0
      end if

      n_data = n_data + 1

      if ( n_max .lt. n_data ) then
        n_data = 0
        x = 0.0D+00
        fx = 0.0D+00
      else
        x = xvec(n_data)
        fx = fxvec(n_data)
      end if

      return
      end
      function pyramid_num ( n )

c*********************************************************************72
c
cc PYRAMID_NUM returns the N-th pyramidal number.
c
c  Discussion:
c
c    The N-th pyramidal number P(N) is formed by the sum of the first
c    N triangular numbers T(J):
c
c      T(J) = sum ( 1 <= J <= N ) J
c
c      P(N) = sum ( 1 <= I <= N ) T(I)
c
c    By convention, T(0) = 0.
c
c    The formula is:
c
c      P(N) = ( (N+1)**3 - (N+1) ) / 6
c
c  Example:
c
c    0   0
c    1   1
c    2   4
c    3  10
c    4  20
c    5  35
c    6  56
c    7  84
c    8 120
c    9 165
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    03 July 2008
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the index of the desired number, which 
c    must be at least 0.
c
c    Output, integer PYRAMID_NUM, the N-th pyramidal number.
c
      implicit none

      integer n
      integer pyramid_num

      pyramid_num = ( ( n + 1 )**3 - ( n + 1 ) ) / 6

      return
      end
      function r8_acosh ( x )

c*********************************************************************72
c
cc R8_ACOSH returns the inverse hyperbolic cosine of a number.
c
c  Discussion:
c
c    One formula is:
c
c      R8_ACOSH = LOG ( X + SQRT ( X * X - 1.0 ) )
c
c    but this formula suffers from roundoff and overflow problems.
c    The formula used here was recommended by W Kahan, as discussed
c    by Moler.
c
c    Applying the inverse function
c
c      Y = R8_ACOSH ( X )
c
c    implies that
c
c      X = COSH(Y) = 0.5 * ( EXP(Y) + EXP(-Y) ).
c
c    For every X greater than or equal to 1, there are two possible
c    choices Y such that X = COSH(Y), differing only in sign.  It
c    is usual to resolve this choice by taking the value of 
c    R8_ACOSH ( X ) to be nonnegative.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    29 November 2007
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Cleve Moler,
c    Trigonometry is a Complex Subject,
c    MATLAB News and Notes,
c    Summer 1998.
c
c  Parameters:
c
c    Input, double precision X, the number whose inverse hyperbolic 
c    cosine is desired.  X should be greater than or equal to 1.
c
c    Output, double precision R8_ACOSH, the inverse hyperbolic cosine of 
c    X.  The principal value (that is, the positive value of the two ) 
c    is returned.
c
      implicit none

      double precision r8_acosh
      double precision x

      if ( x .lt. 1.0D+00 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'R8_ACOSH - Fatal error!'
        write ( *, '(a)' ) '  Argument X must satisfy 1 <= X.'
        write ( *, '(a,g14.6)' ) '  The input X = ', x
        stop
      end if

      r8_acosh = 2.0D+00 * log ( 
     &    sqrt ( 0.5D+00 * ( x + 1.0D+00 ) ) 
     &  + sqrt ( 0.5D+00 * ( x - 1.0D+00 ) ) )

      return
      end
      function r8_asinh ( x )

c*********************************************************************72
c
cc R8_ASINH returns the inverse hyperbolic sine of a number.
c
c  Definition:
c
c    The assertion that:
c
c      Y = R8_ASINH ( X )
c
c    implies that
c
c      X = SINH(Y) = 0.5 * ( EXP(Y) - EXP(-Y) ).
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    29 November 2007
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, double precision X, the number whose inverse hyperbolic 
c    sine is desired.
c
c    Output, double precision R8_ASINH, the inverse hyperbolic sine of X.
c
      implicit none

      double precision r8_asinh
      double precision x

      r8_asinh = log ( x + sqrt ( x * x + 1.0D+00 ) )

      return
      end
      function r8_atanh ( x )

c*********************************************************************72
c
cc R8_ATANH returns the inverse hyperbolic tangent of a number.
c
c  Definition:
c
c    Y = R8_ATANH ( X )
c
c    implies that
c
c    X = TANH(Y) = ( EXP(Y) - EXP(-Y) ) / ( EXP(Y) + EXP(-Y) )
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    29 November 2007
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, double precision X, the number whose inverse hyperbolic 
c    tangent is desired.  The absolute value of X should be less than 
c    or equal to 1.
c
c    Output, double precision R8_ATANH, the inverse hyperbolic tangent of X.
c
      implicit none

      double precision r8_atanh
      double precision x

      if ( 1.0D+00 .le. abs ( x ) ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'R8_ATANH - Fatal error!'
        write ( *, '(a)' ) '  ABS(X) must be < 1.'
        write ( *, '(a,g14.6)' ) '  Your input is X = ', x
        stop
      end if

      r8_atanh = 0.5D+00 * log ( ( 1.0D+00 + x ) / ( 1.0D+00 - x ) )

      return
      end
      function r8_cas ( x )

c*********************************************************************72
c
cc R8_CAS returns the "casine" of an R8 value.
c
c  Discussion:
c
c    The "casine", used in the discrete Hartley transform, is abbreviated
c    CAS(X), and defined by:
c
c      CAS(X) = cos ( X ) + sin( X )
c             = sqrt ( 2 ) * sin ( X + pi/4 )
c             = sqrt ( 2 ) * cos ( X - pi/4 )
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    01 January 2007
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
c    Input, double precision X, the number whose casine is desired.
c
c    Output, double precision R8_CAS, the casine of X, which will be between
c    plus or minus the square root of 2.
c
      implicit none

      double precision r8_cas
      double precision x

      r8_cas = cos ( x ) + sin ( x )

      return
      end
      function r8_choose ( n, k )

c*********************************************************************72
c
cc R8_CHOOSE computes the binomial coefficient C(N,K) as an R8.
c
c  Discussion:
c
c    The value is calculated in such a way as to avoid overflow and
c    roundoff.  The calculation is done in R8 arithmetic.
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
c    07 June 2008
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
c    Volume 6, Number 4, April 1963, page 161.
c
c  Parameters:
c
c    Input, integer N, K, are the values of N and K.
c
c    Output, double precision R8_CHOOSE, the number of combinations of N
c    things taken K at a time.
c
      implicit none

      integer i
      integer k
      integer mn
      integer mx
      integer n
      double precision r8_choose
      double precision value

      mn = min ( k, n - k )

      if ( mn .lt. 0 ) then

        value = 0.0D+00

      else if ( mn .eq. 0 ) then

        value = 1.0D+00

      else

        mx = max ( k, n - k )
        value = dble ( mx + 1 )

        do i = 2, mn
          value = ( value * dble ( mx + i ) ) / dble ( i )
        end do

      end if

      r8_choose = value

      return
      end
      function r8_cot ( angle )

c*********************************************************************72
c
cc R8_COT returns the cotangent of an angle.
c
c  Discussion:
c
c    R8_COT ( THETA ) = COS ( THETA ) / SIN ( THETA )
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    04 July 2008
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, double precision ANGLE, the angle, in radians.
c
c    Output, double precision R8_COT, the cotangent of the angle.
c
      implicit none

      double precision angle
      double precision r8_cot

      r8_cot = cos ( angle ) / sin ( angle )

      return
      end
      function r8_cot_deg ( angle )

c*********************************************************************72
c
cc R8_COT_DEG returns the cotangent of an angle given in degrees.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    04 July 2008
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, double precision ANGLE, the angle, in degrees.
c
c    Output, double precision R8_COT_DEG, the cotangent of the angle.
c
      implicit none

      double precision angle
      double precision r8_cot_deg
      double precision degrees_to_radians
      parameter 
     & ( degrees_to_radians = 3.141592653589793D+00 / 180.0D+00 )

      r8_cot_deg = cos ( degrees_to_radians * angle ) 
     &           / sin ( degrees_to_radians * angle )

      return
      end
      function r8_csc ( theta )

c*********************************************************************72
c
cc R8_CSC returns the cosecant of X.
c
c  Discussion:
c
c    R8_CSC ( THETA ) = 1.0 / SIN ( THETA )
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    04 July 2008
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, double precision THETA, the angle, in radians, whose 
c    cosecant is desired.  It must be the case that SIN ( THETA ) is not zero.
c
c    Output, double precision R8_CSC, the cosecant of THETA.
c
      implicit none

      double precision r8_csc
      double precision theta

      r8_csc = sin ( theta )

      if ( r8_csc .eq. 0.0D+00 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'R8_CSC - Fatal error!'
        write ( *, '(a,g14.6)' ) 
     &  '  Cosecant undefined for THETA = ', theta
        stop
      end if

      r8_csc = 1.0D+00 / r8_csc

      return
      end
      function r8_csc_deg ( angle )

c*********************************************************************72
c
cc R8_CSC_DEG returns the cosecant of an angle given in degrees.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    04 July 2008
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, double precision ANGLE, the angle, in degrees.
c
c    Output, double precision R8_CSC_DEG, the cosecant of the angle.
c
      implicit none

      double precision angle
      double precision r8_csc_deg
      double precision degrees_to_radians
      parameter 
     &  ( degrees_to_radians = 3.141592653589793D+00 / 180.0D+00 )

      r8_csc_deg = 1.0D+00 / cos ( degrees_to_radians * angle )

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
c    05 December 2007
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
      function r8_factorial_log ( n )

c*********************************************************************72
c
cc R8_FACTORIAL_LOG computes log(factorial(N)).
c
c  Discussion:
c
c    The formula is:
c
c      LOG ( FACTORIAL ( N ) ) 
c        = LOG ( product ( 1 <= I <= N ) I )
c        = sum ( ( 1 <= I <= N ) LOG ( I ) )
c  
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    03 July 2008
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the argument of the factorial function.
c    If N is less than 1, the value is returned as 0.
c
c    Output, double precision R8_FACTORIAL_LOG, the logarithm of
c    the factorial of N.
c
      implicit none

      integer i
      integer n
      double precision r8_factorial_log

      r8_factorial_log = 0.0D+00

      do i = 1, n
        r8_factorial_log = r8_factorial_log + log ( dble ( i ) )
      end do

      return
      end
      subroutine r8_factorial_log_values ( n_data, n, fn )

c*********************************************************************72
c
cc R8_FACTORIAL_LOG_VALUES returns values of log(factorial(n)).
c
c  Discussion:
c
c    The function log(factorial(n)) can be written as
c
c     log(factorial(n)) = sum ( 1 <= i <= n ) log ( i )
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    25 March 2007
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Stephen Wolfram,
c    The Mathematica Book,
c    Fourth Edition,
c    Cambridge University Press, 1999,
c    ISBN: 0-521-64314-7,
c    LC: QA76.95.W65.
c
c    Daniel Zwillinger, editor,
c    CRC Standard Mathematical Tables and Formulae,
c    30th Edition,
c    CRC Press, 1996,
c    ISBN: 0-8493-2479-3,
c    LC: QA47.M315.
c
c  Parameters:
c
c    Input/output, integer N_DATA.  The user sets N_DATA to 0 before the
c    first call.  On each call, the routine increments N_DATA by 1, and
c    returns the corresponding data; when there is no more data, the
c    output value of N_DATA will be 0 again.
c
c    Output, integer N, the argument of the function.
c
c    Output, double precision FN, the value of the function.
c
      implicit none

      integer n_max
      parameter ( n_max = 27 )

      double precision fn
      double precision fn_vec(n_max)
      integer n
      integer n_data
      integer n_vec(n_max)

      save fn_vec
      save n_vec

      data fn_vec /
     &  0.0000000000000000D+00,
     &  0.0000000000000000D+00,
     &  0.6931471805599453D+00,
     &  0.1791759469228055D+01,
     &  0.3178053830347946D+01,
     &  0.4787491742782046D+01,
     &  0.6579251212010101D+01,
     &  0.8525161361065414D+01,
     &  0.1060460290274525D+02,
     &  0.1280182748008147D+02,
     &  0.1510441257307552D+02,
     &  0.1750230784587389D+02,
     &  0.1998721449566189D+02,
     &  0.2255216385312342D+02,
     &  0.2519122118273868D+02,
     &  0.2789927138384089D+02,
     &  0.3067186010608067D+02,
     &  0.3350507345013689D+02,
     &  0.3639544520803305D+02,
     &  0.3933988418719949D+02,
     &  0.4233561646075349D+02,
     &  0.5800360522298052D+02,
     &  0.1484777669517730D+03,
     &  0.3637393755555635D+03,
     &  0.6050201058494237D+03,
     &  0.2611330458460156D+04,
     &  0.5912128178488163D+04 /
      data n_vec /
     &     0,
     &     1,
     &     2,
     &     3,
     &     4,
     &     5,
     &     6,
     &     7,
     &     8,
     &     9,
     &    10,
     &    11,
     &    12,
     &    13,
     &    14,
     &    15,
     &    16,
     &    17,
     &    18,
     &    19,
     &    20,
     &    25,
     &    50,
     &   100,
     &   150,
     &   500,
     &  1000 /

      if ( n_data .lt. 0 ) then
        n_data = 0
      end if

      n_data = n_data + 1

      if ( n_max .lt. n_data ) then
        n_data = 0
        n = 0
        fn = 0.0D+00
      else
        n = n_vec(n_data)
        fn = fn_vec(n_data)
      end if

      return
      end
      subroutine r8_factorial_values ( n_data, n, fn )

c*********************************************************************72
c
cc R8_FACTORIAL_VALUES returns values of the real factorial function.
c
c  Discussion:
c
c    Factorial(N) = Product ( 1 <= I <= N ) I
c
c    Although the factorial is an integer valued function, it quickly
c    becomes too large for an integer to hold.  This routine still accepts
c    an integer as the input argument, but returns the function value
c    as a real number.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    25 March 2007
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
c    Output, integer N, the argument of the function.
c
c    Output, double precision FN, the value of the function.
c
      implicit none

      integer n_max
      parameter ( n_max = 25 )

      double precision fn
      double precision fn_vec(n_max)
      integer n
      integer n_data
      integer n_vec(n_max)

      save fn_vec
      save n_vec

      data fn_vec /
     &  0.1000000000000000D+01,
     &  0.1000000000000000D+01,
     &  0.2000000000000000D+01,
     &  0.6000000000000000D+01,
     &  0.2400000000000000D+02,
     &  0.1200000000000000D+03,
     &  0.7200000000000000D+03,
     &  0.5040000000000000D+04,
     &  0.4032000000000000D+05,
     &  0.3628800000000000D+06,
     &  0.3628800000000000D+07,
     &  0.3991680000000000D+08,
     &  0.4790016000000000D+09,
     &  0.6227020800000000D+10,
     &  0.8717829120000000D+11,
     &  0.1307674368000000D+13,
     &  0.2092278988800000D+14,
     &  0.3556874280960000D+15,
     &  0.6402373705728000D+16,
     &  0.1216451004088320D+18,
     &  0.2432902008176640D+19,
     &  0.1551121004333099D+26,
     &  0.3041409320171338D+65,
     &  0.9332621544394415D+158,
     &  0.5713383956445855D+263 /
      data n_vec /
     &     0,
     &     1,
     &     2,
     &     3,
     &     4,
     &     5,
     &     6,
     &     7,
     &     8,
     &     9,
     &    10,
     &    11,
     &    12,
     &    13,
     &    14,
     &    15,
     &    16,
     &    17,
     &    18,
     &    19,
     &    20,
     &    25,
     &    50,
     &   100,
     &   150 /

      if ( n_data .lt. 0 ) then
        n_data = 0
      end if

      n_data = n_data + 1

      if ( n_max .lt. n_data ) then
        n_data = 0
        n = 0
        fn = 0.0D+00
      else
        n = n_vec(n_data)
        fn = fn_vec(n_data)
      end if

      return
      end
      function r8_factorial2 ( n )

c*********************************************************************72
c
cc R8_FACTORIAL2 computes the double factorial function.
c
c  Discussion:
c
c    The formula is:
c
c      FACTORIAL2( N ) = Product ( N * (N-2) * (N-4) * ... * 2 )  (N even)
c                      = Product ( N * (N-2) * (N-4) * ... * 1 )  (N odd)
c
c  Example:
c
c     N    Factorial2(N)
c
c     0     1
c     1     1
c     2     2
c     3     3
c     4     8
c     5    15
c     6    48
c     7   105
c     8   384
c     9   945
c    10  3840
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    03 July 2008
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the argument of the double factorial 
c    function.  If N is less than 1, R8_FACTORIAL2 is returned as 1.0.
c
c    Output, double precision R8_FACTORIAL2, the value of the function.
c
      implicit none

      integer n
      double precision r8_factorial2
      double precision r8_n

      if ( n .lt. 1 ) then
        r8_factorial2 = 1.0D+00
        return
      end if

      r8_n = dble ( n )
      r8_factorial2 = 1.0D+00

10    continue

      if ( 1.0D+00 .lt. r8_n ) then
        r8_factorial2 = r8_factorial2 * r8_n
        r8_n = r8_n - 2.0D+00
        go to 10
      end if

      return
      end
      function r8_gamma ( x )

c*********************************************************************72
c
cc R8_GAMMA evaluates Gamma(X) for a real argument.
c
c  Discussion:
c
c    This routine calculates the gamma function for a real argument X.
c    Computation is based on an algorithm outlined in reference 1.
c    The program uses rational functions that approximate the gamma
c    function to at least 20 significant decimal digits.  Coefficients
c    for the approximation over the interval (1,2) are unpublished.
c    Those for the approximation for 12 <= X are from reference 2.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    18 January 2008
c
c  Author:
c
c    Original FORTRAN77 version by William Cody, Laura Stoltz.
c    This FORTRAN77 version by John Burkardt.
c
c  Reference:
c
c    William Cody,
c    An Overview of Software Development for Special Functions,
c    in Numerical Analysis Dundee, 1975,
c    edited by GA Watson,
c    Lecture Notes in Mathematics 506,
c    Springer, 1976.
c
c    John Hart, Ward Cheney, Charles Lawson, Hans Maehly, 
c    Charles Mesztenyi, John Rice, Henry Thatcher, 
c    Christoph Witzgall,
c    Computer Approximations,
c    Wiley, 1968,
c    LC: QA297.C64.
c
c  Parameters:
c
c    Input, double precision X, the argument of the function.
c
c    Output, double precision R8_GAMMA, the value of the function.
c
      implicit none

      double precision c(7)
      double precision eps
      double precision fact
      integer i
      integer n
      double precision p(8)
      logical parity
      double precision pi
      double precision q(8)
      double precision r8_gamma
      double precision res
      double precision sqrtpi
      double precision sum
      double precision x
      double precision xbig
      double precision xden
      double precision xinf
      double precision xminin
      double precision xnum
      double precision y
      double precision y1
      double precision ysq
      double precision z
c
c  Mathematical constants
c
      data sqrtpi /0.9189385332046727417803297D+00/
      data pi /3.1415926535897932384626434D+00/
c
c  Machine dependent parameters
c
      data xbig / 171.624D+00 /
      data xminin / 2.23D-308 /
      data eps /2.22D-16/
      data xinf /1.79D+308/
c
c  Numerator and denominator coefficients for rational minimax
c  approximation over (1,2).
c
      data p/
     & -1.71618513886549492533811d+00,
     &  2.47656508055759199108314d+01,
     & -3.79804256470945635097577d+02,
     &  6.29331155312818442661052d+02,
     &  8.66966202790413211295064d+02,
     & -3.14512729688483675254357d+04,
     & -3.61444134186911729807069d+04,
     &  6.64561438202405440627855d+04/

      data q/
     & -3.08402300119738975254353d+01,
     &  3.15350626979604161529144d+02,
     & -1.01515636749021914166146d+03,
     & -3.10777167157231109440444d+03,
     &  2.25381184209801510330112d+04,
     &  4.75584627752788110767815d+03,
     & -1.34659959864969306392456d+05,
     & -1.15132259675553483497211d+05/
c
c  Coefficients for minimax approximation over (12, INF).
c
      data c/
     & -1.910444077728D-03,
     &  8.4171387781295D-04,
     & -5.952379913043012D-04,
     &  7.93650793500350248D-04,
     & -2.777777777777681622553D-03,
     &  8.333333333333333331554247D-02,
     &  5.7083835261D-03/

      parity = .false.
      fact = 1.0D+00
      n = 0
      y = x
c
c  Argument is negative.
c
      if ( y .le. 0.0D+00 ) then

        y = - x
        y1 = aint ( y )
        res = y - y1

        if ( res .ne. 0.0D+00 ) then

          if ( y1 .ne. aint ( y1 * 0.5D+00 ) * 2.0D+00 ) then
            parity = .true.
          end if

          fact = - pi / sin ( pi * res )
          y = y + 1.0D+00

        else

          res = xinf
          r8_gamma = res
          return

        end if

      end if
c
c  Argument is positive.
c
      if ( y .lt. eps ) then
c
c  Argument < EPS.
c
        if ( xminin .le. y ) then
          res = 1.0D+00 / y
        else
          res = xinf
          r8_gamma = res
          return
        end if

      else if ( y .lt. 12.0D+00 ) then

        y1 = y
c
c  0.0 < argument < 1.0.
c
        if ( y .lt. 1.0D+00 ) then

          z = y
          y = y + 1.0D+00
c
c  1.0 < argument < 12.0.
c  Reduce argument if necessary.
c
        else

          n = int ( y ) - 1
          y = y - dble ( n )
          z = y - 1.0D+00

        end if
c
c  Evaluate approximation for 1.0 < argument < 2.0.
c
        xnum = 0.0D+00
        xden = 1.0D+00
        do i = 1, 8
          xnum = ( xnum + p(i) ) * z
          xden = xden * z + q(i)
        end do

        res = xnum / xden + 1.0D+00
c
c  Adjust result for case  0.0 < argument < 1.0.
c
        if ( y1 .lt. y ) then

          res = res / y1
c
c  Adjust result for case 2.0 < argument < 12.0.
c
        else if ( y .lt. y1 ) then

          do i = 1, n
            res = res * y
            y = y + 1.0D+00
          end do

        end if

      else
c
c  Evaluate for 12.0 <= argument.
c
        if ( y .le. xbig ) then

          ysq = y * y
          sum = c(7)
          do i = 1, 6
            sum = sum / ysq + c(i)
          end do
          sum = sum / y - y + sqrtpi
          sum = sum + ( y - 0.5D+00 ) * log ( y )
          res = exp ( sum )

        else

          res = xinf
          r8_gamma = res
          return

        end if

      end if
c
c  Final adjustments and return.
c
      if ( parity ) then
        res = - res
      end if

      if ( fact .ne. 1.0D+00 ) then
        res = fact / res
      end if

      r8_gamma = res

      return
      end
      function r8_gamma_log ( x )

c*********************************************************************72
c
cc R8_GAMMA_LOG evaluates log ( Gamma ( X ) ) for a real argument.
c
c  Discussion:
c
c    This routine calculates the LOG(GAMMA) function for a positive real
c    argument X.  Computation is based on an algorithm outlined in
c    references 1 and 2.  The program uses rational functions that
c    theoretically approximate LOG(GAMMA) to at least 18 significant
c    decimal digits.  The approximation for X > 12 is from reference
c    3, while approximations for X < 12.0 are similar to those in
c    reference 1, but are unpublished.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    08 July 2008
c
c  Author:
c
c    Original FORTRAN77 version by William Cody, Laura Stoltz.
c    This FORTRAN77 version by John Burkardt.
c
c  Reference:
c
c    William Cody, Kenneth Hillstrom,
c    Chebyshev Approximations for the Natural Logarithm of the 
c    Gamma Function,
c    Mathematics of Computation,
c    Volume 21, Number 98, April 1967, pages 198-203.
c
c    Kenneth Hillstrom,
c    ANL/AMD Program ANLC366S, DGAMMA/DLGAMA,
c    May 1969.
c
c    John Hart, Ward Cheney, Charles Lawson, Hans Maehly, 
c    Charles Mesztenyi, John Rice, Henry Thatcher, 
c    Christoph Witzgall,
c    Computer Approximations,
c    Wiley, 1968,
c    LC: QA297.C64.
c
c  Parameters:
c
c    Input, double precision X, the argument of the function.
c
c    Output, double precision R8_GAMMA_LOG, the value of the function.
c
      implicit none

      double precision c(7)
      double precision corr
      double precision d1
      double precision d2
      double precision d4
      double precision eps
      double precision frtbig
      integer i
      double precision pnt68
      double precision p1(8)
      double precision p2(8)
      double precision p4(8)
      double precision q1(8)
      double precision q2(8)
      double precision q4(8)
      double precision r8_gamma_log
      double precision res
      double precision sqrtpi
      double precision x
      double precision xbig
      double precision xden
      double precision xinf
      double precision xm1
      double precision xm2
      double precision xm4
      double precision xnum
      double precision y
      double precision ysq
c
c  Mathematical constants
c
      data pnt68 /0.6796875D+00/
      data sqrtpi /0.9189385332046727417803297D+00/
c
c  Machine dependent parameters
c
      data xbig /2.55D+305/
      data xinf /1.79D+308/
      data eps /2.22D-16/
      data frtbig /2.25D+76/
c
c  Numerator and denominator coefficients for rational minimax
c  approximation over (0.5,1.5).
c
      data d1/-5.772156649015328605195174D-01/
      data p1/
     &   4.945235359296727046734888D+00,
     &   2.018112620856775083915565D+02,
     &   2.290838373831346393026739D+03,
     &   1.131967205903380828685045D+04,
     &   2.855724635671635335736389D+04,
     &   3.848496228443793359990269D+04,
     &   2.637748787624195437963534D+04,
     &   7.225813979700288197698961D+03/
      data q1/
     &   6.748212550303777196073036D+01,
     &   1.113332393857199323513008D+03,
     &   7.738757056935398733233834D+03,
     &   2.763987074403340708898585D+04,
     &   5.499310206226157329794414D+04,
     &   6.161122180066002127833352D+04,
     &   3.635127591501940507276287D+04,
     &   8.785536302431013170870835D+03/
c
c  Numerator and denominator coefficients for rational minimax
c  Approximation over (1.5,4.0).
c
      data d2/4.227843350984671393993777D-01/
      data p2/
     &   4.974607845568932035012064D+00,
     &   5.424138599891070494101986D+02,
     &   1.550693864978364947665077D+04,
     &   1.847932904445632425417223D+05,
     &   1.088204769468828767498470D+06,
     &   3.338152967987029735917223D+06,
     &   5.106661678927352456275255D+06,
     &   3.074109054850539556250927D+06/
      data q2/
     &   1.830328399370592604055942D+02,
     &   7.765049321445005871323047D+03,
     &   1.331903827966074194402448D+05,
     &   1.136705821321969608938755D+06,
     &   5.267964117437946917577538D+06,
     &   1.346701454311101692290052D+07,
     &   1.782736530353274213975932D+07,
     &   9.533095591844353613395747D+06/
c
c  Numerator and denominator coefficients for rational minimax
c  Approximation over (4.0,12.0).
c
      data d4/1.791759469228055000094023D+00/
      data p4/
     &   1.474502166059939948905062D+04,
     &   2.426813369486704502836312D+06,
     &   1.214755574045093227939592D+08,
     &   2.663432449630976949898078D+09,
     &   2.940378956634553899906876D+10,
     &   1.702665737765398868392998D+11,
     &   4.926125793377430887588120D+11,
     &   5.606251856223951465078242D+11/
      data q4/
     &   2.690530175870899333379843D+03,
     &   6.393885654300092398984238D+05,
     &   4.135599930241388052042842D+07,
     &   1.120872109616147941376570D+09,
     &   1.488613728678813811542398D+10,
     &   1.016803586272438228077304D+11,
     &   3.417476345507377132798597D+11,
     &   4.463158187419713286462081D+11/
c
c  Coefficients for minimax approximation over (12, INF).
c
      data c/
     &  -1.910444077728D-03,
     &   8.4171387781295D-04,
     &  -5.952379913043012D-04,
     &   7.93650793500350248D-04,
     &  -2.777777777777681622553D-03,
     &   8.333333333333333331554247D-02,
     &   5.7083835261D-03/

      y = x

      if ( 0.0D+00 .lt. y .and. y .le. xbig ) then

        if ( y .le. eps ) then

          res = - dlog ( y )
c
c  EPS < X <= 1.5.
c
        else if ( y .le. 1.5D+00 ) then

          if ( y .lt. pnt68 ) then
            corr = - dlog ( y )
            xm1 = y
          else
            corr = 0.0D+00
            xm1 = ( y - 0.5D+00 ) - 0.5D+00
          end if

          if ( y .le. 0.5D+00 .or. pnt68 .le. y ) then

            xden = 1.0D+00
            xnum = 0.0D+00
            do i = 1, 8
              xnum = xnum * xm1 + p1(i)
              xden = xden * xm1 + q1(i)
            end do

            res = corr + ( xm1 * ( d1 + xm1 * ( xnum / xden ) ) )

          else

            xm2 = ( y - 0.5D+00 ) - 0.5D+00
            xden = 1.0D+00
            xnum = 0.0D+00
            do i = 1, 8
              xnum = xnum * xm2 + p2(i)
              xden = xden * xm2 + q2(i)
            end do

            res = corr + xm2 * ( d2 + xm2 * ( xnum / xden ) )

          end if
c
c  1.5 < X <= 4.0.
c
        else if ( y .le. 4.0D+00 ) then

          xm2 = y - 2.0D+00
          xden = 1.0D+00
          xnum = 0.0D+00
          do i = 1, 8
            xnum = xnum * xm2 + p2(i)
            xden = xden * xm2 + q2(i)
          end do

          res = xm2 * ( d2 + xm2 * ( xnum / xden ) )
c
c  4.0 < X <= 12.0.
c
        else if ( y .le. 12.0D+00 ) then

          xm4 = y - 4.0D+00
          xden = - 1.0D+00
          xnum = 0.0D+00
          do i = 1, 8
            xnum = xnum * xm4 + p4(i)
            xden = xden * xm4 + q4(i)
          end do

          res = d4 + xm4 * ( xnum / xden )
c
c  Evaluate for 12 <= argument.
c
        else

          res = 0.0D+00

          if ( y .le. frtbig ) then

            res = c(7)
            ysq = y * y

            do i = 1, 6
              res = res / ysq + c(i)
            end do

          end if

          res = res / y
          corr = dlog ( y )
          res = res + sqrtpi - 0.5D+00 * corr
          res = res + y * ( corr - 1.0D+00 )

        end if
c
c  Return for bad arguments.
c
      else

        res = xinf

      end if
c
c  Final adjustments and return.
c
      r8_gamma_log = res

      return
      end
      function r8_huge ( )

c*********************************************************************72
c
cc R8_HUGE returns a "huge" R8.
c
c  Discussion:
c
c    The value returned by this function is NOT required to be the
c    maximum representable R8.  This value varies from machine to machine,
c    from compiler to compiler, and may cause problems when being printed.
c    We simply want a "very large" but non-infinite number.
c
c    FORTRAN90 provides a built-in routine HUGE ( X ) that
c    can return the maximum representable number of the same datatype
c    as X, if that is what is really desired.
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
      subroutine r8_hyper_2f1 ( a_input, b_input, c_input, x_input, hf )

c*********************************************************************72
c
cc R8_HYPER_2F1 evaluates the hypergeometric function F(A,B,C,X).
c
c  Discussion:
c
c    A minor bug was corrected.  The HW variable, used in several places as
c    the "old" value of a quantity being iteratively improved, was not
c    being initialized.  JVB, 11 February 2008.
c
c    The original version of this program allowed the input arguments to
c    be modified, although they were restored to their input values before exit.
c    This is unacceptable if the input arguments are allowed to be constants.
c    The code has been modified so that the input arguments are never modified.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    21 March 2009
c
c  Author:
c
c    Original FORTRAN77 version by Shanjie Zhang, Jianming Jin.
c    This FORTRAN77 version by John Burkardt.
c
c    The original FORTRAN77 version of this routine is copyrighted by
c    Shanjie Zhang and Jianming Jin.  However, they give permission to
c    incorporate this routine into a user program provided that the copyright
c    is acknowledged.
c
c  Reference:
c
c    Shanjie Zhang, Jianming Jin,
c    Computation of Special Functions,
c    Wiley, 1996,
c    ISBN: 0-471-11963-6,
c    LC: QA351.C45
c
c  Parameters:
c
c    Input, double precision A_INPUT, B_INPUT, C_INPUT, X_INPUT, 
c    the arguments of the function.  The user is allowed to pass these
c    values as constants or variables.
c    C_INPUT must not be equal to a nonpositive integer.
c    X_INPUT .lt. 1.
c
c    Output, double precision HF, the value of the function.
c
      implicit none

      double precision a
      double precision a_input
      double precision a0
      double precision aa
      double precision b
      double precision b_input
      double precision bb
      double precision c
      double precision c_input
      double precision c0
      double precision c1
      double precision el
      parameter ( el = 0.5772156649015329D+00 )
      double precision eps
      double precision f0
      double precision f1
      double precision g0
      double precision g1
      double precision g2
      double precision g3
      double precision ga
      double precision gabc
      double precision gam
      double precision gb
      double precision gbm
      double precision gc
      double precision gca
      double precision gcab
      double precision gcb
      double precision gm
      double precision hf
      double precision hw
      integer j
      integer k
      logical l0
      logical l1
      logical l2
      logical l3
      logical l4
      logical l5
      integer m
      integer nm
      double precision pa
      double precision pb
      double precision pi
      parameter ( pi = 3.141592653589793D+00 )
      double precision r
      double precision r0
      double precision r1
      double precision r8_gamma
      double precision r8_psi
      double precision rm
      double precision rp
      double precision sm
      double precision sp
      double precision sp0
      double precision x
      double precision x_input
      double precision x1
c
c  Immediately copy the input argumentsc
c
      a = a_input
      b = b_input
      c = c_input
      x = x_input

      l0 = ( c .eq. aint ( c ) ) .and. ( c .lt. 0.0D+00 )
      l1 = ( 1.0D+00 - x .lt. 1.0D-15 ) .and. ( c - a - b .le. 0.0D+00 )
      l2 = ( a .eq. aint ( a ) ) .and. ( a .lt. 0.0D+00 )
      l3 = ( b .eq. aint ( b ) ) .and. ( b .lt. 0.0D+00 )
      l4 = ( c - a .eq. aint ( c - a ) ) .and. ( c - a .le. 0.0D+00 )
      l5 = ( c - b .eq. aint ( c - b ) ) .and. ( c - b .le. 0.0D+00 )

      if ( l0 .or. l1 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'R8_HYPER_2F1 - Fatal error!'
        write ( *, '(a)' ) '  The hypergeometric series is divergent.'
        return
      end if

      if ( 0.95D+00 .lt. x ) then
        eps = 1.0D-08
      else
        eps = 1.0D-15
      end if

      if ( x .eq. 0.0D+00 .or. a .eq. 0.0D+00 .or. b .eq. 0.0D+00 ) then

        hf = 1.0D+00
        return

      else if ( 1.0D+00 - x .eq. eps .and. 0.0D+00 .lt. c - a - b ) then

        gc = r8_gamma ( c )
        gcab = r8_gamma ( c - a - b )
        gca = r8_gamma ( c - a )
        gcb = r8_gamma ( c - b )
        hf = gc * gcab / ( gca * gcb )
        return

      else if ( 1.0D+00 + x .le. eps .and. 
     &  abs ( c - a + b - 1.0D+00 ) .le. eps ) then

        g0 = sqrt ( pi ) * 2.0D+00**( - a )
        g1 = r8_gamma ( c )
        g2 = r8_gamma ( 1.0D+00 + a / 2.0D+00 - b )
        g3 = r8_gamma ( 0.5D+00 + 0.5D+00 * a )
        hf = g0 * g1 / ( g2 * g3 )
        return

      else if ( l2 .or. l3 ) then

        if ( l2 ) then
          nm = int ( abs ( a ) )
        end if

        if ( l3 ) then
          nm = int ( abs ( b ) )
        end if

        hf = 1.0D+00
        r = 1.0D+00

        do k = 1, nm
          r = r * ( a + k - 1.0D+00 ) * ( b + k - 1.0D+00 ) 
     &      / ( k * ( c + k - 1.0D+00 ) ) * x
          hf = hf + r
        end do

        return

      else if ( l4 .or. l5 ) then

        if ( l4 ) then
          nm = int ( abs ( c - a ) )
        end if

        if ( l5 ) then
          nm = int ( abs ( c - b ) )
        end if

        hf = 1.0D+00
        r  = 1.0D+00
        do k = 1, nm
          r = r * ( c - a + k - 1.0D+00 ) * ( c - b + k - 1.0D+00 ) 
     &      / ( k * ( c + k - 1.0D+00 ) ) * x
          hf = hf + r
        end do
        hf = ( 1.0D+00 - x )**( c - a - b ) * hf
        return

      end if

      aa = a
      bb = b
      x1 = x

      if ( x .lt. 0.0D+00 ) then
        x = x / ( x - 1.0D+00 )
        if ( a .lt. c .and. b .lt. a .and. 0.0D+00 .lt. b ) then
          a = bb
          b = aa
        end if
        b = c - b
      end if

      if ( 0.75D+00 .le. x ) then

        gm = 0.0D+00

        if ( abs ( c - a - b - aint ( c - a - b ) ) .lt. 1.0D-15 ) then

          m = int ( c - a - b )
          ga = r8_gamma ( a )
          gb = r8_gamma ( b )
          gc = r8_gamma ( c )
          gam = r8_gamma ( a + m )
          gbm = r8_gamma ( b + m )

          pa = r8_psi ( a )
          pb = r8_psi ( b )

          if ( m /= 0 ) then
            gm = 1.0D+00
          end if

          do j = 1, abs ( m ) - 1
            gm = gm * j
          end do

          rm = 1.0D+00
          do j = 1, abs ( m )
            rm = rm * j
          end do

          f0 = 1.0D+00
          r0 = 1.0D+00
          r1 = 1.0D+00
          sp0 = 0.0D+00
          sp = 0.0D+00

          if ( 0 .le. m ) then

            c0 = gm * gc / ( gam * gbm )
            c1 = - gc * ( x - 1.0D+00 )**m / ( ga * gb * rm )

            do k = 1, m - 1
              r0 = r0 * ( a + k - 1.0D+00 ) * ( b + k - 1.0D+00 ) 
     &          / ( k * ( k - m ) ) * ( 1.0D+00 - x )
              f0 = f0 + r0
            end do

            do k = 1, m
              sp0 = sp0 + 1.0D+00 / ( a + k - 1.0D+00 ) 
     &          + 1.0D+00 / ( b + k - 1.0D+00 ) - 1.0D+00 / dble ( k )
            end do

            f1 = pa + pb + sp0 + 2.0D+00 * el + log ( 1.0D+00 - x )
            hw = f1

            do k = 1, 250

              sp = sp + ( 1.0D+00 - a ) / ( k * ( a + k - 1.0D+00 ) ) 
     &          + ( 1.0D+00 - b ) / ( k * ( b + k - 1.0D+00 ) )

              sm = 0.0D+00
              do j = 1, m
                sm = sm + ( 1.0D+00 - a ) 
     &            / ( ( j + k ) * ( a + j + k - 1.0D+00 ) ) 
     &            + 1.0D+00 / ( b + j + k - 1.0D+00 )
              end do

              rp = pa + pb + 2.0D+00 * el + sp + sm 
     &          + log ( 1.0D+00 - x )

              r1 = r1 * ( a + m + k - 1.0D+00 ) 
     &          * ( b + m + k - 1.0D+00 ) 
     &          / ( k * ( m + k ) ) * ( 1.0D+00 - x )

              f1 = f1 + r1 * rp

              if ( abs ( f1 - hw ) .lt. abs ( f1 ) * eps ) then
                exit
              end if

              hw = f1

            end do

            hf = f0 * c0 + f1 * c1

          else if ( m .lt. 0 ) then

            m = - m
            c0 = gm * gc / ( ga * gb * ( 1.0D+00 - x )**m )
            c1 = - ( - 1 )**m * gc / ( gam * gbm * rm )

            do k = 1, m - 1
              r0 = r0 * ( a - m + k - 1.0D+00 ) 
     &          * ( b - m + k - 1.0D+00 ) 
     &          / ( k * ( k - m ) ) * ( 1.0D+00 - x )
              f0 = f0 + r0
            end do

            do k = 1, m
              sp0 = sp0 + 1.0D+00 / dble ( k )
            end do

            f1 = pa + pb - sp0 + 2.0D+00 * el + log ( 1.0D+00 - x )
            hw = f1

            do k = 1, 250

              sp = sp + ( 1.0D+00 - a ) 
     &          / ( k * ( a + k - 1.0D+00 ) ) 
     &          + ( 1.0D+00 - b ) / ( k * ( b + k - 1.0D+00 ) )

              sm = 0.0D+00
              do j = 1, m
                sm = sm + 1.0D+00 / dble ( j + k )
              end do

              rp = pa + pb + 2.0D+00 * el + sp - sm 
     &          + log ( 1.0D+00 - x )

              r1 = r1 * ( a + k - 1.0D+00 ) * ( b + k - 1.0D+00 ) 
     &          / ( k * ( m + k ) ) * ( 1.0D+00 - x )

              f1 = f1 + r1 * rp

              if ( abs ( f1 - hw ) .lt. abs ( f1 ) * eps ) then
                exit
              end if

              hw = f1

            end do

            hf = f0 * c0 + f1 * c1

          end if

        else

          ga = r8_gamma ( a )
          gb = r8_gamma ( b )
          gc = r8_gamma ( c )
          gca = r8_gamma ( c - a )
          gcb = r8_gamma ( c - b )
          gcab = r8_gamma ( c - a - b )
          gabc = r8_gamma ( a + b - c )
          c0 = gc * gcab / ( gca * gcb )
          c1 = gc * gabc / ( ga * gb ) * ( 1.0D+00 - x )**( c - a - b )
          hf = 0.0D+00
          hw = hf
          r0 = c0
          r1 = c1

          do k = 1, 250

            r0 = r0 * ( a + k - 1.0D+00 ) * ( b + k - 1.0D+00 ) 
     &        / ( k * ( a + b - c + k ) ) * ( 1.0D+00 - x )

            r1 = r1 * ( c - a + k - 1.0D+00 ) * ( c - b + k - 1.0D+00 )
     &        / ( k * ( c - a - b + k ) ) * ( 1.0D+00 - x )

            hf = hf + r0 + r1

            if ( abs ( hf - hw ) .lt. abs ( hf ) * eps ) then
              exit
            end if

            hw = hf

          end do

          hf = hf + c0 + c1

        end if

      else

        a0 = 1.0D+00

        if ( a .lt. c .and. c .lt. 2.0D+00 * a .and. 
     &       b .lt. c .and. c .lt. 2.0D+00 * b ) then

          a0 = ( 1.0D+00 - x )**( c - a - b )
          a = c - a
          b = c - b

        end if

        hf = 1.0D+00
        hw = hf
        r = 1.0D+00

        do k = 1, 250

          r = r * ( a + k - 1.0D+00 ) * ( b + k - 1.0D+00 ) 
     &      / ( k * ( c + k - 1.0D+00 ) ) * x

          hf = hf + r

          if ( abs ( hf - hw ) .le. abs ( hf ) * eps ) then
            exit
          end if

          hw = hf

        end do

        hf = a0 * hf

      end if

      if ( x1 .lt. 0.0D+00 ) then
        x = x1
        c0 = 1.0D+00 / ( 1.0D+00 - x )**aa
        hf = c0 * hf
      end if

      a = aa
      b = bb

      if ( 120 .lt. k ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'R8_HYPER_2F1 - Warning!'
        write ( *, '(a)' ) '  A large number of iterations were needed.'
        write ( *, '(a)' ) 
     &    '  The accuracy of the results should be checked.'
      end if

      return
      end
      function r8_nint ( x )

c*****************************************************************************80
c
cc R8_NINT returns the nearest integer to an R8.
c
c  Example:
c
c        X        R8_NINT
c
c      1.3         1
c      1.4         1
c      1.5         1 or 2
c      1.6         2
c      0.0         0
c     -0.7        -1
c     -1.1        -1
c     -1.6        -2
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    08 September 2005
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, double precision X, the value.
c
c    Output, integer R8_NINT, the nearest integer to X.
c
      implicit none

      integer r8_nint
      integer s
      double precision x

      if ( x .lt. 0.0D+00 ) then
        s = -1
      else
        s = 1
      end if

      r8_nint = s * int ( abs ( x ) + 0.5D+00 )

      return
      end
      function r8_pi ( )

c*********************************************************************72
c
cc R8_PI returns the value of pi as an R8.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
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
      function r8_psi ( xx )

c*********************************************************************72
c
cc R8_PSI evaluates the function Psi(X).
c
c  Discussion:
c
c    This routine evaluates the logarithmic derivative of the
c    GAMMA function,
c
c      PSI(X) = d/dX (GAMMA(X)) / GAMMA(X) 
c             = d/dX LN ( GAMMA(X) )
c
c    for real X, where either
c
c      -XMAX1 < X < -XMIN  and X is not a negative integer), 
c
c    or
c
c      XMIN < X.
c
c  Modified:
c
c    23 January 2008
c
c  Author:
c
c    William Cody
c
c  Reference:
c
c    William Cody, Anthony Strecok, Henry Thacher,
c    Chebyshev Approximations for the Psi Function,
c    Mathematics of Computation,
c    Volume 27, Number 121, January 1973, pages 123-127.
c
c  Parameters:
c
c    Input, double precision XX, the argument of the function.
c
c    Output, double precision R8_PSI, the value of the function.
c
      implicit none

      double precision aug
      double precision den
      integer i
      integer n
      integer nq
      double precision p1(9)
      double precision p2(7)
      double precision piov4
      double precision q1(8)
      double precision q2(6)
      double precision r8_psi
      double precision sgn
      double precision xlarge
      double precision upper
      double precision w
      double precision x
      double precision xinf
      double precision xmax1
      double precision xmin1
      double precision xsmall
      double precision x01
      double precision x01d
      double precision x02
      double precision xx
      double precision z
c
c  Mathematical constants.  PIOV4 = pi / 4
c
      data piov4 /7.8539816339744830962d-01/
c
c  Machine-dependent constants
c
      data xinf /1.70d+38/
      data xmin1 /5.89d-39/
      data xmax1 /3.60d+16/
      data xsmall /2.05d-09/
      data xlarge /2.04d+15/
c
c  Zero of psi(x)
c
      data x01 /187.0d0/
      data x01d /128.0d0/
      data x02 /6.9464496836234126266d-04/
c
c  Coefficients for approximation to  psi(x)/(x-x0)  over [0.5, 3.0]
c
      data p1/4.5104681245762934160d-03,5.4932855833000385356d+00,
     &        3.7646693175929276856d+02,7.9525490849151998065d+03,
     &        7.1451595818951933210d+04,3.0655976301987365674d+05,
     &        6.3606997788964458797d+05,5.8041312783537569993d+05,
     &        1.6585695029761022321d+05/
      data q1/9.6141654774222358525d+01,2.6287715790581193330d+03,
     &        2.9862497022250277920d+04,1.6206566091533671639d+05,
     &        4.3487880712768329037d+05,5.4256384537269993733d+05,
     &        2.4242185002017985252d+05,6.4155223783576225996d-08/
c
c  Coefficients for approximation to  psi(x) - ln(x) + 1/(2x)
c  for 3.0 < x.
c
      data p2/-2.7103228277757834192d+00,-1.5166271776896121383d+01,
     &        -1.9784554148719218667d+01,-8.8100958828312219821d+00,
     &        -1.4479614616899842986d+00,-7.3689600332394549911d-02,
     &        -6.5135387732718171306d-21/
      data q2/ 4.4992760373789365846d+01, 2.0240955312679931159d+02,
     &         2.4736979003315290057d+02, 1.0742543875702278326d+02,
     &         1.7463965060678569906d+01, 8.8427520398873480342d-01/

      x = xx
      w = abs ( x )
      aug = 0.0D+00
c
c  Check for valid arguments, then branch to appropriate algorithm.
c
      if ( - x .ge. xmax1 .or. w .lt. xmin1 ) then
        r8_psi = xinf
        if ( 0.0D+00 .lt. x ) then
          r8_psi = -xinf
        end if
        return
      end if

      if ( x .ge. 0.5D+00 ) then
        go to 200
c
c  X < 0.5, use reflection formula: psi(1-x) = psi(x) + pi * cot(pi*x)
c  Use 1/X for PI*COTAN(PI*X)  when  XMIN1 < |X| <= XSMALL.
c
      else if ( w .le. xsmall ) then
        aug = - 1.0D+00 / x
        go to 150
      end if
c
c  Argument reduction for cotangent.
c
  100 continue

      if ( x .lt. 0.0D+00 ) then
        sgn = piov4
      else
        sgn = - piov4
      end if

      w = w - aint ( w )
      nq = int ( w * 4.0D+00 )
      w = 4.0D+00 * ( w - dble ( nq ) * 0.25D+00 )
c
c  W is now related to the fractional part of 4.0 * X.
c  Adjust argument to correspond to values in the first
c  quadrant and determine the sign.
c
      n = nq / 2

      if ( n + n .ne. nq ) then
        w = 1.0D+00 - w
      end if

      z = piov4 * w

      if ( mod ( n, 2 ) .ne. 0 ) then
        sgn = - sgn
      end if
c
c  Determine the final value for  -pi * cotan(pi*x).
c
      n = ( nq + 1 ) / 2
      if ( mod ( n, 2 ) .eq. 0 ) then
c
c  Check for singularity.
c
        if ( z .eq. 0.0D+00 ) then
          r8_psi = xinf
          if ( 0.0D+00 .lt. x ) then
            r8_psi = -xinf
          end if
          return
        end if

        aug = sgn * ( 4.0D+00 / tan ( z ) )

      else
        aug = sgn * ( 4.0D+00 * tan ( z ) )
      end if

  150 continue

      x = 1.0D+00 - x

  200 continue
c
c  0.5 <= X <= 3.0.
c
      if ( x .le. 3.0D+00 ) then

        den = x
        upper = p1(1) * x
        do i = 1, 7
          den = ( den + q1(i) ) * x
          upper = ( upper + p1(i+1) ) * x
        end do
        den = ( upper + p1(9) ) / ( den + q1(8) )
        x = ( x - x01 / x01d ) - x02
        r8_psi = den * x + aug
        return

      end if
c
c  3.0 < X.
c
      if ( x .lt. xlarge ) then
        w = 1.0D+00 / ( x * x )
        den = w
        upper = p2(1) * w
        do i = 1, 5
          den = ( den + q2(i) ) * w
          upper = ( upper + p2(i+1) ) * w
        end do
        aug = ( upper + p2(7) ) / ( den + q2(6) ) - 0.5D+00 / x + aug
      end if

      r8_psi = aug + log ( x )

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
      subroutine r8poly_print ( n, a, title )

c*********************************************************************72
c
cc R8POLY_PRINT prints out a polynomial.
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
c    Input, double precision A(0:N), the polynomial coefficients.
c    A(0) is the constant term and
c    A(N) is the coefficient of X^N.
c
c    Input, character * ( * ) TITLE, an optional title.
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

      if ( 0 .lt. title_length ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) title(1:title_length)
      end if

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

      do i = n2-1, 0, -1

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
      subroutine r8poly_val_horner ( n, c, x, cx )

c*********************************************************************72
c
cc R8POLY_VAL_HORNER evaluates a polynomial using Horner's method.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    06 July 2008
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the degree of the polynomial.
c
c    Input, double precision C(0:N), the polynomial coefficients.
c    C(I) is the coefficient of X^I.
c
c    Input, double precision X, the point at which the polynomial 
c    is to be evaluated.
c
c    Output, double precision CX, the value of the polynomial at X.
c
      implicit none

      integer n

      double precision c(0:n)
      double precision cx
      integer i
      double precision x

      cx = c(n)
      do i = n - 1, 0, -1
        cx = cx * x + c(i)
      end do

      return
      end
      function r8poly_value ( n, a, x )

c*********************************************************************72
c
cc R8POLY_VALUE evaluates an R8POLY
c
c  Discussion:
c
c    For sanity's sake, the value of N indicates the NUMBER of 
c    coefficients, or more precisely, the ORDER of the polynomial,
c    rather than the DEGREE of the polynomial.  The two quantities
c    differ by 1, but cause a great deal of confusion.
c
c    Given N and A, the form of the polynomial is:
c
c      p(x) = a(1) + a(2) * x + ... + a(n-1) * x^(n-2) + a(n) * x^(n-1)
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    05 August 2010
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the order of the polynomial.
c
c    Input, double precision A(N), the coefficients of the polynomial.
c    A(1) is the constant term.
c
c    Input, double precision X, the point at which the polynomial is 
c    to be evaluated.
c
c    Output, double precision R8POLY_VALUE, the value of the polynomial at X.
c
      implicit none

      integer n

      double precision a(n)
      integer i
      double precision r8poly_value
      double precision x

      r8poly_value = a(n)
      do i = n - 1, 1, -1
        r8poly_value = r8poly_value * x + a(i)
      end do

      return
      end
      subroutine r8vec_print_some ( n, a, max_print, title )

c*********************************************************************72
c
cc R8VEC_PRINT_SOME prints "some" of an R8VEC.
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
c    16 September 2003
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
c    Input, character*(*) TITLE, an optional title.
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

      if ( 0 .lt. s_len_trim ( title ) ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) title
        write ( *, '(a)' ) ' '
      end if

      if ( n .le. max_print ) then

        do i = 1, n
          write ( *, '(i6,2x,g14.6)' ) i, a(i)
        end do

      else if ( 3 .le. max_print ) then

        do i = 1, max_print-2
          write ( *, '(i6,2x,g14.6)' ) i, a(i)
        end do

        write ( *, '(a)' ) '......  ..............'
        i = n

        write ( *, '(i6,2x,g14.6)' ) i, a(i)

      else

        do i = 1, max_print-1
          write ( *, '(i6,2x,g14.6)' ) i, a(i)
        end do

        i = max_print

        write ( *, '(i6,2x,g14.6,a)' ) i, a(i), '...more entries...'

      end if

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
      function sec_deg ( angle )

c*********************************************************************72
c
cc SEC_DEG returns the secant of an angle given in degrees.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    06 July 2008
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, double precision ANGLE, the angle, in degrees.
c
c    Output, double precision SEC_DEG, the secant of the angle.
c
      implicit none

      double precision angle
      double precision degrees_to_radians
      parameter 
     & ( degrees_to_radians = 3.141592653589793D+00 / 180.0D+00 )
      double precision sec_deg

      sec_deg = 1.0D+00 / sin ( degrees_to_radians * angle )

      return
      end
      subroutine sigma ( n, sigma_n )

c*********************************************************************72
c
cc SIGMA returns the value of SIGMA(N), the divisor sum.
c
c  Discussion:
c
c    SIGMA(N) is the sum of the distinct divisors of N, including 1 and N.
c
c    The formula is:
c
c      SIGMA(U*V) = SIGMA(U) * SIGMA(V) if U and V are relatively prime.
c
c      SIGMA(P**K) = ( P**(K+1) - 1 ) / ( P - 1 ) if P is prime.
c
c  Example:
c
c     N  SIGMA(N)
c
c     1    1
c     2    3
c     3    4
c     4    7
c     5    6
c     6   12
c     7    8
c     8   15
c     9   13
c    10   18
c    11   12
c    12   28
c    13   14
c    14   24
c    15   24
c    16   31
c    17   18
c    18   39
c    19   20
c    20   42
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
c    Input, integer N, the value to be analyzed.
c
c    Output, integer SIGMA_N, the value of SIGMA(N).  If N is
c    less than or equal to 0, SIGMA_N will be returned as 0.  If there is not
c    enough room for factoring N, SIGMA_N is returned as -1.
c
      implicit none

      integer maxfactor
      parameter ( maxfactor = 20 )

      integer factor(maxfactor)
      integer i
      integer n
      integer nfactor
      integer nleft
      integer power(maxfactor)
      integer sigma_n

      if ( n .le. 0 ) then
        sigma_n = 0
        return
      end if

      if ( n .eq. 1 ) then
        sigma_n = 1
        return
      end if
!
!  Factor N.
!
      call i4_factor ( n, maxfactor, nfactor, factor, power, nleft )

      if ( nleft .ne. 1 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'SIGMA - Fatal error!'
        write ( *, '(a)' ) '  Not enough factorization space.'
        sigma_n = -1
        return
      end if

      sigma_n = 1
      do i = 1, nfactor
        sigma_n = ( sigma_n * ( factor(i)**( power(i) + 1 ) - 1 ) ) 
     &    / ( factor(i) - 1 )
      end do

      return
      end
      subroutine sigma_values ( n_data, n, c )

c*********************************************************************72
c
cc SIGMA_VALUES returns some values of the Sigma function.
c
c  Discussion:
c
c    SIGMA(N) is the sum of the distinct divisors of N, including 1 and N.
c
c    In Mathematica, the function can be evaluated by:
c
c      DivisorSigma[1,n]
c
c    The formula is:
c
c    SIGMA(U*V) = SIGMA(U) * SIGMA(V) if U and V are relatively prime.
c
c    SIGMA(P**K) = ( P**(K+1) - 1 ) / ( P - 1 ) if P is prime.
c
c  First values:
c
c     N  SIGMA(N)
c
c     1    1
c     2    3
c     3    4
c     4    7
c     5    6
c     6   12
c     7    8
c     8   15
c     9   13
c    10   18
c    11   12
c    12   28
c    13   14
c    14   24
c    15   24
c    16   31
c    17   18
c    18   39
c    19   20
c    20   42
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    25 March 2007
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
c    Output, integer N, the argument of the Sigma function.
c
c    Output, integer C, the value of the Sigma function.
c
      implicit none

      integer n_max
      parameter ( n_max = 20 )

      integer c
      integer c_vec(n_max)
      integer n
      integer n_data
      integer n_vec(n_max)

      save c_vec
      save n_vec

      data c_vec /
     &   1,    3,    4,    7,    6,   12,    8,   15,   13,   18,
     &  72,  128,  255,  176,  576, 1170,  618,  984, 2232, 2340 /
      data n_vec /
     &    1,   2,   3,   4,   5,   6,   7,   8,   9,   10,
     &   30, 127, 128, 129, 210, 360, 617, 815, 816, 1000 /

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
      function sin_deg ( angle )

c*********************************************************************72
c
cc SIN_DEG returns the sine of an angle given in degrees.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    06 July 2008
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, double precision ANGLE, the angle, in degrees.
c
c    Output, double precision SIN_DEG, the sine of the angle.
c
      implicit none

      double precision angle
      double precision degrees_to_radians
      parameter 
     & ( degrees_to_radians = 3.141592653589793D+00 / 180.0D+00 )
      double precision sin_deg

      sin_deg = sin ( degrees_to_radians * angle ) 

      return
      end
      function sin_power_int ( a, b, n )

c*********************************************************************72
c
cc SIN_POWER_INT evaluates the sine power integral.
c
c  Discussion:
c
c    The function is defined by
c
c      SIN_POWER_INT(A,B,N) = Integral ( A <= T <= B ) ( sin ( t ))^n dt
c
c    The algorithm uses the following fact:
c
c      Integral sin^n ( t ) = (1/n) * (
c        sin^(n-1)(t) * cos(t) + ( n-1 ) * Integral sin^(n-2) ( t ) dt )
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    06 July 2008
c
c  Author:
c
c    John Burkardt
c
c  Parameters
c
c    Input, double precision A, B, the limits of integration.
c
c    Input, integer N, the power of the sine function.
c
c    Output, double precision SIN_POWER_INT, the value of the integral.
c
      implicit none

      double precision a
      double precision b
      double precision ca
      double precision cb
      integer m
      integer mlo
      integer n
      double precision sa
      double precision sb
      double precision sin_power_int
      double precision value

      if ( n .lt. 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'SIN_POWER_INT - Fatal error!'
        write ( *, '(a)' ) '  Power N < 0.'
        value = 0.0
        stop
      end if

      sa = sin ( a )
      sb = sin ( b )
      ca = cos ( a )
      cb = cos ( b )

      if ( mod ( n, 2 ) .eq. 0 ) then
        value = b - a
        mlo = 2
      else
        value = ca - cb
        mlo = 3
      end if

      do m = mlo, n, 2
        value = ( dble ( m - 1 ) * value 
     &            + sa**(m-1) * ca - sb**(m-1) * cb ) 
     &    / dble ( m )
      end do

      sin_power_int = value

      return
      end
      subroutine sin_power_int_values ( n_data, a, b, n, fx )

c*********************************************************************72
c
cc SIN_POWER_INT_VALUES returns some values of the sine power integral.
c
c  Discussion:
c
c    The function has the form
c
c      SIN_POWER_INT(A,B,N) = Integral ( A <= T <= B ) ( sin(T) )^N dt
c
c    In Mathematica, the function can be evaluated by:
c
c      Integrate [ ( Sin[x] )^n, { x, a, b } ]
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    25 March 2007
c
c  Author:
c
c    John Burkardt
c
c  Reference:
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
c    Output, double precision A, B, the limits of integration.
c
c    Output, integer N, the power.
c
c    Output, double precision FX, the value of the function.
c
      implicit none

      integer n_max
      parameter ( n_max = 10 )

      double precision a
      double precision a_vec(n_max)
      double precision b
      double precision b_vec(n_max)
      double precision fx
      double precision fx_vec(n_max)
      integer n
      integer n_data
      integer n_vec(n_max)

      save a_vec
      save b_vec
      save fx_vec
      save n_vec

      data a_vec /
     &   0.10D+02,
     &   0.00D+00,
     &   0.00D+00,
     &   0.00D+00,
     &   0.00D+00,
     &   0.00D+00,
     &   0.00D+00,
     &   0.10D+01,
     &   0.00D+00,
     &   0.00D+00 /
      data b_vec /
     &   0.20D+02,
     &   0.10D+01,
     &   0.10D+01,
     &   0.10D+01,
     &   0.10D+01,
     &   0.10D+01,
     &   0.20D+01,
     &   0.20D+01,
     &   0.10D+01,
     &   0.10D+01 /
      data fx_vec /
     &  0.10000000000000000000D+02,
     &  0.45969769413186028260D+00,
     &  0.27267564329357957615D+00,
     &  0.17894056254885809051D+00,
     &  0.12402556531520681830D+00,
     &  0.88974396451575946519D-01,
     &  0.90393123848149944133D+00,
     &  0.81495684202992349481D+00,
     &  0.21887522421729849008D-01,
     &  0.17023439374069324596D-01 /
      data n_vec /
     &   0,
     &   1,
     &   2,
     &   3,
     &   4,
     &   5,
     &   5,
     &   5,
     &  10,
     &  11 /

      if ( n_data .lt. 0 ) then
        n_data = 0
      end if

      n_data = n_data + 1

      if ( n_max .lt. n_data ) then
        n_data = 0
        a = 0.0D+00
        b = 0.0D+00
        n = 0
        fx = 0.0D+00
      else
        a = a_vec(n_data)
        b = b_vec(n_data)
        n = n_vec(n_data)
        fx = fx_vec(n_data)
      end if

      return
      end
      subroutine slice ( dim_num, slice_num, piece_num )

c*********************************************************************72
c
cc SLICE: maximum number of pieces created by a given number of slices.
c
c  Discussion:
c
c    If we imagine slicing a pizza, each slice produce more pieces.  
c    The position of the slice affects the number of pieces created, but there
c    is a maximum.  
c
c    This function determines the maximum number of pieces created by a given
c    number of slices, applied to a space of a given dimension.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    12 August 2011
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Robert Banks,
c    Slicing Pizzas, Racing Turtles, and Further Adventures in 
c    Applied Mathematics,
c    Princeton, 1999,
c    ISBN13: 9780691059471,
c    LC: QA93.B358.
c
c  Parameters:
c
c    Input, integer DIM_NUM, the spatial dimension.
c
c    Input, integer SLICE_NUM, the number of slices.
c
c    Input, real PIECE_NUM, the maximum number of pieces that can
c    be created by the given number of slices applied in the given dimension.
c
      implicit none

      integer dim_num
      integer i4_choose
      integer j
      integer piece_num
      integer slice_num

      piece_num = 0
      do j = 0, min ( dim_num, slice_num )
        piece_num = piece_num + i4_choose ( slice_num, j )
      end do

      return
      end
      subroutine spherical_harmonic ( l, m, theta, phi, c, s )

c*********************************************************************72
c
cc SPHERICAL_HARMONIC evaluates spherical harmonic functions.
c
c  Discussion:
c
c    The spherical harmonic function Y(L,M,THETA,PHI,X) is the
c    angular part of the solution to Laplace's equation in spherical
c    coordinates.
c
c    Y(L,M,THETA,PHI,X) is related to the associated Legendre
c    function as follows:
c
c      Y(L,M,THETA,PHI,X) = FACTOR * P(L,M,cos(THETA)) * exp ( i * M * PHI )
c
c    Here, FACTOR is a normalization factor:
c
c      FACTOR = sqrt ( ( 2 * L + 1 ) * ( L - M )! / ( 4 * PI * ( L + M )! ) )
c
c    In Mathematica, a spherical harmonic function can be evaluated by
c
c      SphericalHarmonicY [ l, m, theta, phi ]
c
c    Note that notational tradition in physics requires that THETA
c    and PHI represent the reverse of what they would normally mean
c    in mathematical notation; that is, THETA goes up and down, and
c    PHI goes around.
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
c    Milton Abramowitz, Irene Stegun,
c    Handbook of Mathematical Functions,
c    National Bureau of Standards, 1964,
c    ISBN: 0-486-61272-4,
c    LC: QA47.A34.
c
c    Eric Weisstein,
c    CRC Concise Encyclopedia of Mathematics,
c    CRC Press, 2002,
c    Second edition,
c    ISBN: 1584883472,
c    LC: QA5.W45
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
c    Input, integer L, the first index of the spherical harmonic
c    function.  Normally, 0 <= L.
c
c    Input, integer M, the second index of the spherical harmonic 
c    function.  Normally, -L <= M <= L.
c
c    Input, double precision THETA, the polar angle, for which
c    0 <= THETA <= PI.
c
c    Input, double precision PHI, the longitudinal angle, for which
c    0 <= PHI <= 2*PI.
c
c    Output, double precision C(0:L), S(0:L), the real and imaginary
c    parts of the functions Y(L,0:L,THETA,PHI).
c
      implicit none

      integer l

      double precision c(0:l)
      integer i
      integer m
      integer m_abs
      double precision phi
      double precision plm(0:l)
      double precision s(0:l)
      double precision theta

      m_abs = abs ( m )

      call legendre_associated_normalized ( l, m_abs, cos ( theta ), 
     &  plm )

      do i = 0, l
        c(i) = plm(i) * cos ( dble ( m ) * phi )
        s(i) = plm(i) * sin ( dble ( m ) * phi )
      end do

      if ( m .lt. 0 ) then
        do i = 0, l
          c(i) = - c(i)
          s(i) = - s(i)
        end do
      end if

      return
      end
      subroutine spherical_harmonic_values ( n_data, l, m, theta, phi,
     &  yr, yi )

c*********************************************************************72
c
cc SPHERICAL_HARMONIC_VALUES returns values of spherical harmonic functions.
c
c  Discussion:
c
c    In Mathematica, the function can be evaluated by
c
c      SphericalHarmonicY [ l, m, theta, phi ]
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    25 March 2007
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
c    Eric Weisstein,
c    CRC Concise Encyclopedia of Mathematics,
c    CRC Press, 1998.
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
c    Input/output, integer N_DATA.
c    On input, if N_DATA is 0, the first test data is returned, and
c    N_DATA is set to the index of the test data.  On each subsequent
c    call, N_DATA is incremented and that test data is returned.  When
c    there is no more test data, N_DATA is set to 0.
c
c    Output, integer L, integer M, double precision THETA, PHI, the arguments
c    of the function.
c
c    Output, double precision YR, YI, the real and imaginary parts of
c    the function.
c
      implicit none

      integer n_max
      parameter ( n_max = 20 )

      integer l
      integer l_vec(n_max)
      integer m
      integer m_vec(n_max)
      integer n_data
      double precision phi
      double precision phi_vec(n_max)
      double precision theta
      double precision theta_vec(n_max)
      double precision yi
      double precision yi_vec(n_max)
      double precision yr
      double precision yr_vec(n_max)

      save l_vec
      save m_vec
      save phi_vec
      save theta_vec
      save yi_vec
      save yr_vec

      data l_vec /
     &   0,  1,  2,
     &   3,  4,  5,
     &   5,  5,  5,
     &   5,  4,  4,
     &   4,  4,  4,
     &   3,  3,  3,
     &   3,  3 /
      data m_vec /
     &   0,  0,  1,
     &   2,  3,  5,
     &   4,  3,  2,
     &   1,  2,  2,
     &   2,  2,  2,
     &  -1, -1, -1,
     &  -1, -1 /
      data phi_vec /
     &  0.1047197551196598D+01,
     &  0.1047197551196598D+01,
     &  0.1047197551196598D+01,
     &  0.1047197551196598D+01,
     &  0.1047197551196598D+01,
     &  0.6283185307179586D+00,
     &  0.6283185307179586D+00,
     &  0.6283185307179586D+00,
     &  0.6283185307179586D+00,
     &  0.6283185307179586D+00,
     &  0.7853981633974483D+00,
     &  0.7853981633974483D+00,
     &  0.7853981633974483D+00,
     &  0.7853981633974483D+00,
     &  0.7853981633974483D+00,
     &  0.4487989505128276D+00,
     &  0.8975979010256552D+00,
     &  0.1346396851538483D+01,
     &  0.1795195802051310D+01,
     &  0.2243994752564138D+01 /
      data theta_vec /
     &  0.5235987755982989D+00,
     &  0.5235987755982989D+00,
     &  0.5235987755982989D+00,
     &  0.5235987755982989D+00,
     &  0.5235987755982989D+00,
     &  0.2617993877991494D+00,
     &  0.2617993877991494D+00,
     &  0.2617993877991494D+00,
     &  0.2617993877991494D+00,
     &  0.2617993877991494D+00,
     &  0.6283185307179586D+00,
     &  0.1884955592153876D+01,
     &  0.3141592653589793D+01,
     &  0.4398229715025711D+01,
     &  0.5654866776461628D+01,
     &  0.3926990816987242D+00,
     &  0.3926990816987242D+00,
     &  0.3926990816987242D+00,
     &  0.3926990816987242D+00,
     &  0.3926990816987242D+00 /
      data yi_vec /
     &  0.0000000000000000D+00,
     &  0.0000000000000000D+00,
     & -0.2897056515173922D+00,
     &  0.1916222768312404D+00,
     &  0.0000000000000000D+00,
     &  0.0000000000000000D+00,
     &  0.3739289485283311D-02,
     & -0.4219517552320796D-01,
     &  0.1876264225575173D+00,
     & -0.3029973424491321D+00,
     &  0.4139385503112256D+00,
     & -0.1003229830187463D+00,
     &  0.0000000000000000D+00,
     & -0.1003229830187463D+00,
     &  0.4139385503112256D+00,
     & -0.1753512375142586D+00,
     & -0.3159720118970196D+00,
     & -0.3940106541811563D+00,
     & -0.3940106541811563D+00,
     & -0.3159720118970196D+00 /
      data yr_vec /
     &  0.2820947917738781D+00,
     &  0.4231421876608172D+00,
     & -0.1672616358893223D+00,
     & -0.1106331731112457D+00,
     &  0.1354974113737760D+00,
     &  0.5390423109043568D-03,
     & -0.5146690442951909D-02,
     &  0.1371004361349490D-01,
     &  0.6096352022265540D-01,
     & -0.4170400640977983D+00,
     &  0.0000000000000000D+00,
     &  0.0000000000000000D+00,
     &  0.0000000000000000D+00,
     &  0.0000000000000000D+00,
     &  0.0000000000000000D+00,
     &  0.3641205966137958D+00,
     &  0.2519792711195075D+00,
     &  0.8993036065704300D-01,
     & -0.8993036065704300D-01,
     & -0.2519792711195075D+00 /

      if ( n_data .lt. 0 ) then
        n_data = 0
      end if

      n_data = n_data + 1

      if ( n_max .lt. n_data ) then
        n_data = 0
        l = 0
        m = 0
        theta = 0.0D+00
        phi = 0.0D+00
        yr = 0.0D+00
        yi = 0.0D+00
      else
        l = l_vec(n_data)
        m = m_vec(n_data)
        theta = theta_vec(n_data)
        phi = phi_vec(n_data)
        yr = yr_vec(n_data)
        yi = yi_vec(n_data)
      end if

      return
      end
      subroutine stirling1 ( n, m, s1 )

c*********************************************************************72
c
cc STIRLING1 computes the Stirling numbers of the first kind.
c
c  Discussion:
c
c    The absolute value of the Stirling number S1(N,M) gives the number
c    of permutations on N objects having exactly M cycles, while the
c    sign of the Stirling number records the sign (odd or even) of
c    the permutations.  For example, there are six permutations on 3 objects:
c
c      A B C   3 cycles (A) (B) (C)
c      A C B   2 cycles (A) (BC)
c      B A C   2 cycles (AB) (C)
c      B C A   1 cycle  (ABC)
c      C A B   1 cycle  (ABC)
c      C B A   2 cycles (AC) (B)
c
c    There are 
c
c      2 permutations with 1 cycle, and S1(3,1) = 2
c      3 permutations with 2 cycles, and S1(3,2) = -3,
c      1 permutation with 3 cycles, and S1(3,3) = 1.
c
c    Since there are N! permutations of N objects, the sum of the absolute 
c    values of the Stirling numbers in a given row, 
c
c      sum ( 1 <= I <= N ) abs ( S1(N,I) ) = N!
c
c  First terms:
c
c    N/M:  1     2      3     4     5    6    7    8
c
c    1     1     0      0     0     0    0    0    0
c    2    -1     1      0     0     0    0    0    0
c    3     2    -3      1     0     0    0    0    0
c    4    -6    11     -6     1     0    0    0    0
c    5    24   -50     35   -10     1    0    0    0
c    6  -120   274   -225    85   -15    1    0    0
c    7   720 -1764   1624  -735   175  -21    1    0
c    8 -5040 13068 -13132  6769 -1960  322  -28    1
c
c  Recursion:
c
c    S1(N,1) = (-1)**(N-1) * (N-1)! for all N.
c    S1(I,I) = 1 for all I.
c    S1(I,J) = 0 if I < J.
c
c    S1(N,M) = S1(N-1,M-1) - (N-1) * S1(N-1,M)
c
c  Properties:
c
c    sum ( 1 <= K <= M ) S2(I,K) * S1(K,J) = Delta(I,J)
c
c    X_N = sum ( 0 <= K <= N ) S1(N,K) X^K
c    where X_N is the falling factorial function.
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
c    Input, integer N, the number of rows of the table.
c
c    Input, integer M, the number of columns of the table.
c
c    Output, integer S1(N,M), the Stirling numbers of the 
c    first kind.
c
      implicit none

      integer m
      integer n

      integer i
      integer j
      integer s1(n,m)

      if ( n .le. 0 ) then
        return
      end if

      if ( m .le. 0 ) then
        return
      end if

      s1(1,1) = 1
      do j = 2, m
        s1(1,j) = 0
      end do

      do i = 2, n

        s1(i,1) = - ( i - 1 ) * s1(i-1,1)

        do j = 2, m
          s1(i,j) = s1(i-1,j-1) - ( i - 1 ) * s1(i-1,j)
        end do

      end do
 
      return
      end
      subroutine stirling2 ( n, m, s2 )

c*********************************************************************72
c
cc STIRLING2 computes the Stirling numbers of the second kind.
c
c  Discussion:
c
c    S2(N,M) represents the number of distinct partitions of N elements
c    into M nonempty sets.  For a fixed N, the sum of the Stirling
c    numbers S2(N,M) is represented by B(N), called "Bell's number",
c    and represents the number of distinct partitions of N elements.
c
c    For example, with 4 objects, there are:
c
c    1 partition into 1 set:
c
c      (A,B,C,D)
c
c    7 partitions into 2 sets:
c
c      (A,B,C) (D)
c      (A,B,D) (C)
c      (A,C,D) (B)
c      (A) (B,C,D)
c      (A,B) (C,D)
c      (A,C) (B,D)
c      (A,D) (B,C)
c
c    6 partitions into 3 sets:
c
c      (A,B) (C) (D)
c      (A) (B,C) (D)
c      (A) (B) (C,D)
c      (A,C) (B) (D)
c      (A,D) (B) (C)
c      (A) (B,D) (C)
c
c    1 partition into 4 sets:
c
c      (A) (B) (C) (D)
c
c    So S2(4,1) = 1, S2(4,2) = 7, S2(4,3) = 6, S2(4,4) = 1, and B(4) = 15.
c
c    The Stirling numbers of the second kind S(N,1:N) are the coefficients of
c    the Bell polynomial B(N,X):
c
c      B(0,X) = 1
c      B(N,X) = sum ( 1 <= M <= N ) S(N,M) * X^M
c
c  First terms:
c
c    N/M: 1    2    3    4    5    6    7    8
c
c    1    1    0    0    0    0    0    0    0
c    2    1    1    0    0    0    0    0    0
c    3    1    3    1    0    0    0    0    0
c    4    1    7    6    1    0    0    0    0
c    5    1   15   25   10    1    0    0    0
c    6    1   31   90   65   15    1    0    0
c    7    1   63  301  350  140   21    1    0
c    8    1  127  966 1701 1050  266   28    1
c
c  Recursion:
c
c    S2(N,1) = 1 for all N.
c    S2(I,I) = 1 for all I.
c    S2(I,J) = 0 if I < J.
c
c    S2(N,M) = M * S2(N-1,M) + S2(N-1,M-1)
c
c  Properties:
c
c    sum ( 1 <= K <= M ) S2(I,K) * S1(K,J) = Delta(I,J)
c
c    X^N = sum ( 0 <= K <= N ) S2(N,K) X_K
c    where X_K is the falling factorial function.
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
c    Input, integer N, the number of rows of the table.
c
c    Input, integer M, the number of columns of the table.
c
c    Output, integer S2(N,M), the Stirling numbers of the 
c    second kind.
c
      implicit none

      integer m
      integer n

      integer i
      integer j
      integer s2(n,m)

      if ( n .le. 0 ) then
        return
      end if

      if ( m .le. 0 ) then
        return
      end if

      s2(1,1) = 1
      do j = 2, m
        s2(1,j) = 0
      end do

      do i = 2, n

        s2(i,1) = 1

        do j = 2, m
          s2(i,j) = j * s2(i-1,j) + s2(i-1,j-1)
        end do

      end do
 
      return
      end
      function tan_deg ( angle )

c*********************************************************************72
c
cc TAN_DEG returns the tangent of an angle given in degrees.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    06 July 2008
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, double precision ANGLE, the angle, in degrees.
c
c    Output, double precision TAN_DEG, the tangent of the angle.
c
      implicit none

      double precision angle
      double precision degrees_to_radians
      parameter 
     & ( degrees_to_radians = 3.141592653589793D+00 / 180.0D+00 )
      double precision tan_deg

      tan_deg = sin ( degrees_to_radians * angle ) 
     &        / cos ( degrees_to_radians * angle )

      return
      end
      subroutine tau ( n, taun )

c*********************************************************************72
c
cc TAU returns the value of TAU(N), the number of distinct divisors of N.
c
c  Discussion:
c
c    TAU(N) is the number of distinct divisors of N, including 1 and N.
c
c    If the prime factorization of N is
c
c      N = P1**E1 * P2**E2 * ... * PM**EM,
c
c    then
c
c      TAU(N) = ( E1 + 1 ) * ( E2 + 1 ) * ... * ( EM + 1 ).
c
c    One consequence of this fact is that TAU is odd if and only
c    if N is a perfect square.
c
c  First values:
c
c     N   TAU(N)
c
c     1    1
c     2    2
c     3    2
c     4    3
c     5    2
c     6    4
c     7    2
c     8    4
c     9    3
c    10    4
c    11    2
c    12    6
c    13    2
c    14    4
c    15    4
c    16    5
c    17    2
c    18    6
c    19    2
c    20    6
c    21    4
c    22    4
c    23    2
c    24    8
c    25    3
c    26    4
c    27    4
c    28    6
c    29    2
c    30    8
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    06 July 2008
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the value to be analyzed.  N must be 1 or
c    greater.
c
c    Output, integer TAUN, the value of TAU(N).  But if N is 0 or
c    less, TAUN is returned as 0, a nonsense value.  If there is
c    not enough room for factoring, TAUN is returned as -1.
c
      implicit none

      integer maxfactor
      parameter ( maxfactor = 20 )

      integer factor(maxfactor)
      integer i
      integer n
      integer nfactor
      integer nleft
      integer power(maxfactor)
      integer taun

      if ( n .le. 0 ) then
        taun = 0
        return
      end if

      if ( n .eq. 1 ) then
        taun = 1
        return
      end if
c
c  Factor N.
c
      call i4_factor ( n, maxfactor, nfactor, factor, power, nleft )

      if ( nleft .ne. 1 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'TAU - Fatal error!'
        write ( *, '(a)' ) '  Not enough factorization space.'
        taun = -1
        return
      end if

      taun = 1
      do i = 1, nfactor
        taun = taun * ( power(i) + 1 )
      end do

      return
      end
      subroutine tau_values ( n_data, n, c )

c*********************************************************************72
c
cc TAU_VALUES returns some values of the Tau function.
c
c  Discussion:
c
c    TAU(N) is the number of divisors of N, including 1 and N.
c
c    In Mathematica, the function can be evaluated by:
c
c      DivisorSigma[1,n]
c
c    If the prime factorization of N is
c
c      N = P1**E1 * P2**E2 * ... * PM**EM,
c
c    then
c
c      TAU(N) = ( E1 + 1 ) * ( E2 + 1 ) * ... * ( EM + 1 ).
c
c  First values:
c
c     N   TAU(N)
c
c     1    1
c     2    2
c     3    2
c     4    3
c     5    2
c     6    4
c     7    2
c     8    4
c     9    3
c    10    4
c    11    2
c    12    6
c    13    2
c    14    4
c    15    4
c    16    5
c    17    2
c    18    6
c    19    2
c    20    6
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    25 March 2007
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
c    Output, integer N, the argument of the Tau function.
c
c    Output, integer C, the value of the Tau function.
c
      implicit none

      integer n_max
      parameter ( n_max = 20 )

      integer c
      integer c_vec(n_max)
      integer n
      integer n_data
      integer n_vec(n_max)

      save c_vec
      save n_vec

      data c_vec /
     &  1,  2,  2,  3,  2,  4,  2,  4,  3,  4,
     &  2, 12, 12,  4, 18, 24,  2,  8, 14, 28 /
      data n_vec /
     &    1,   2,   3,   4,   5,   6,   7,   8,   9,  10,
     &   23,  72, 126, 226, 300, 480, 521, 610, 832, 960 /

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
      function tetrahedron_num ( n )

c*********************************************************************72
c
cc TETRAHEDRON_NUM returns the N-th tetrahedral number.
c
c  Discussion:
c
c    The N-th tetrahedral number T3(N) is formed by the sum of the first
c    N triangular numbers:
c
c      T3(N) = sum ( 1 <= I <= N ) T2(I)
c            = sum ( 1 <= I <= N ) sum ( 1 <= J < I ) J
c
c    By convention, T3(0) = 0.
c
c    The formula is:
c
c      T3(N) = ( N * ( N + 1 ) * ( N + 2 ) ) / 6
c
c  First Values:
c
c     0
c     1
c     4
c    10
c    20
c    35
c    56
c    84
c   120
c   165
c   220
c   275
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    06 July 2008
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the index of the desired number, which 
c    must be at least 0.
c
c    Output, integer TETRAHEDRON_NUM, the N-th tetrahedron number.
c
      implicit none

      integer n
      integer tetrahedron_num

      tetrahedron_num = ( n * ( n + 1 ) * ( n + 2 ) ) / 6

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
      function triangle_num ( n )

c*********************************************************************72
c
cc TRIANGLE_NUM returns the N-th triangular number.
c
c  Discussion:
c
c    The N-th triangular number T(N) is formed by the sum of the first
c    N integers:
c
c      T(N) = sum ( 1 <= I <= N ) I
c
c    By convention, T(0) = 0.
c
c    T(N) can be computed quickly by the formula:
c
c      T(N) = ( N * ( N + 1 ) ) / 2
c
c  First Values:
c
c     0
c     1
c     3
c     6
c    10
c    15
c    21
c    28
c    36
c    45
c    55
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    04 July 2008
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the index of the desired number, 
c    which must be at least 0.
c
c    Output, integer TRIANGLE_NUM, the N-th triangular number.
c
      implicit none

      integer n
      integer triangle_num

      triangle_num = ( n * ( n + 1 ) ) / 2

      return
      end
      subroutine triangle_to_i4 ( i, j, k )

c*********************************************************************72
c
cc TRIANGLE_TO_I4 converts a triangular coordinate to an integer.
c
c  Discussion:
c
c    Triangular coordinates are handy when storing a naturally triangular
c    array (such as the lower half of a matrix) in a linear array.
c
c    Thus, for example, we might consider storing 
c
c    (1,1)
c    (2,1) (2,2)
c    (3,1) (3,2) (3,3)
c    (4,1) (4,2) (4,3) (4,4)
c
c    as the linear array
c
c    (1,1) (2,1) (2,2) (3,1) (3,2) (3,3) (4,1) (4,2) (4,3) (4,4)    
c
c    Here, the quantities in parenthesis represent the natural row and
c    column indices of a single number when stored in a rectangular array.
c
c    Thus, our goal is, given the row I and column J of the data,
c    to produce the value K which indicates its position in the linear
c    array.
c
c    The triangular numbers are the indices associated with the
c    diagonal elements of the original array, T(1,1), T(2,2), T(3,3)
c    and so on.
c
c    The formula is:
c
c      K = J + ( (I-1) * I ) / 2
c
c  First Values:
c
c     I  J  K
c
c     0  0  0
c     1  1  1
c     2  1  2
c     2  2  3
c     3  1  4
c     3  2  5
c     3  3  6
c     4  1  7
c     4  2  8
c     4  3  9
c     4  4 10
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    04 July 2008
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer I, J, the row and column indices.  I and J must
c    be nonnegative, and J must not be greater than I.
c
c    Output, integer K, the linear index of the (I,J) element.
c
      implicit none

      integer i
      integer j
      integer k

      if ( i .lt. 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'TRIANGLE_TO_I4 - Fatal error!'
        write ( *, '(a)' ) '  I < 0.'
        write ( *, '(a,i8)' ) '  I = ', i
        stop
      else if ( j .lt. 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'TRIANGLE_TO_I4 - Fatal error!'
        write ( *, '(a)' ) '  J < 0.'
        write ( *, '(a,i8)' ) '  J = ', j
        stop
      else if ( i .lt. j ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'TRIANGLE_TO_I4 - Fatal error!'
        write ( *, '(a)' ) '  I < J.'
        write ( *, '(a,i8)' ) '  I = ', i
        write ( *, '(a,i8)' ) '  J = ', j
        stop
      end if

      k = j + ( ( i - 1 ) * i ) / 2

      return
      end
      subroutine vibonacci ( n, seed, v )

c*********************************************************************72
c
cc VIBONACCI computes the first N Vibonacci numbers.
c
c  Discussion:
c
c    The "Vibonacci numbers" are a generalization of the Fibonacci numbers:
c      V(N+1) = +/- V(N) +/- V(N-1)
c    where the signs are chosen randomly.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    18 July 2008
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Brian Hayes,
c    The Vibonacci Numbers,
c    American Scientist,
c    July-August 1999, Volume 87, Number 4.
c
c    Divakar Viswanath,
c    Random Fibonacci sequences and the number 1.13198824,
c    Mathematics of Computation,
c    1998.
c
c  Parameters:
c
c    Input, integer N, the highest number to compute.
c
c    Input/output, integer SEED, a seed for the random number 
c    generator.
c
c    Output, integer V(N), the first N Vibonacci numbers.  By 
c    convention, V(1) and V(2) are taken to be 1.
c
      implicit none

      integer n

      integer i
      integer i4_uniform
      integer j
      integer s1
      integer s2
      integer seed
      integer v(n)

      if ( n .le. 0 ) then
        return
      end if

      v(1) = 1

      if ( n .le. 1 ) then
        return
      end if

      v(2) = 1

      do i = 3, n
        
        j = i4_uniform ( 0, 1, seed )

        if ( j .eq. 0 ) then
          s1 = -1
        else
          s1 = +1
        end if

        j = i4_uniform ( 0, 1, seed )

        if ( j .eq. 0 ) then
          s2 = -1
        else
          s2 = +1
        end if

        v(i) = s1 * v(i-1) + s2 * v(i-2)

      end do
     
      return
      end
      subroutine zeckendorf ( n, m_max, m, i_list, f_list )

c*********************************************************************72
c
cc ZECKENDORF produces the Zeckendorf decomposition of a positive integer.
c
c  Discussion:
c
c    Zeckendorf proved that every positive integer can be represented
c    uniquely as the sum of non-consecutive Fibonacci numbers.
c
c    N = sum ( 1 <= I <= M ) F_LIST(I)
c
c  Example:
c
c     N    Decomposition
c
c    50    34 + 13 + 3
c    51    34 + 13 + 3 + 1
c    52    34 + 13 + 5
c    53    34 + 13 + 5 + 1
c    54    34 + 13 + 5 + 2
c    55    55
c    56    55 + 1
c    57    55 + 2
c    58    55 + 3
c    59    55 + 3 + 1
c    60    55 + 5
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    18 July 2008
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the positive integer to be decomposed.
c
c    Input, integer M_MAX, the maximum dimension of I_LIST 
c    and F_LIST.
c
c    Output, integer M, the number of parts in the decomposition.
c
c    Output, integer I_LIST(M_MAX), contains in entries 1 
c    through M the index of the Fibonacci numbers in the decomposition.
c
c    Output, integer F_LIST(M_MAX), contains in entries 1 
c    through M the value of the Fibonacci numbers in the decomposition.
c
      implicit none

      integer m_max

      integer f
      integer f_list(m_max)
      integer i
      integer i_list(m_max)
      integer j
      integer m
      integer n
      integer n_copy

      m = 0

      n_copy = n
c
c  Extract a sequence of Fibonacci numbers.
c
10    continue

      if ( 0 .lt. n_copy .and. m .lt. m_max ) then
        call fibonacci_floor ( n_copy, f, i )
        m = m + 1
        i_list(m) = i
        n_copy = n_copy - f
        go to 10
      end if
c
c  Replace any pair of consecutive indices ( I, I-1 ) by I+1.
c
      do i = m, 2, -1

        if ( i_list(i-1) .eq. i_list(i) + 1 ) then
          i_list(i-1) = i_list(i-1) + 1
          do j = i, m - 1
            i_list(j) = i_list(j+1)
          end do
          i_list(m) = 0
          m = m - 1
        end if

      end do
c
c  Fill in the actual values of the Fibonacci numbers.
c
      do i = 1, m
        call fibonacci_direct ( i_list(i), f_list(i) )
      end do

      return
      end
      subroutine zernike_poly ( m, n, rho, z )

!*********************************************************************72
!
!! ZERNIKE_POLY evaluates a Zernike polynomial at RHO.
!
!  Discussion:
!
!    This routine uses the facts that:
!
!    *) R^M_N = 0 if M < 0, or N < 0, or N < M.
!    *) R^M_M = RHO^M
!    *) R^M_N = 0 if mod ( N - M ) = 1.
!
!    and the recursion:
!
!    R^M_(N+2) = A * [ ( B * RHO * RHO - C ) * R^M_N - D * R^M_(N-2) ]
!
!    where 
!
!    A = ( N + 2 ) / ( ( N + 2 )^2 - M * M )
!    B = 4 * ( N + 1 )
!    C = ( N + M )^2 / N + ( N - M + 2 )^2 / ( N + 2 )
!    D = ( N^2 - M^2 ) / N
!
!    I wish I could clean up the recursion in the code, but for
!    now, I have to treat the case M = 0 specially.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    18 July 2008
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Eric Weisstein,
!    CRC Concise Encyclopedia of Mathematics,
!    CRC Press, 2002,
!    Second edition,
!    ISBN: 1584883472,
!    LC: QA5.W45
!
!  Parameters:
!
!    Input, integer M, the upper index.
!
!    Input, integer N, the lower index.
!
!    Input, double precision RHO, the radial coordinate.
!
!    Output, double precision Z, the value of the Zernike 
!    polynomial R^M_N at the point RHO.
!
      implicit none

      double precision a
      double precision b
      double precision c
      double precision d
      integer m
      integer n
      integer nn
      double precision rho
      double precision z
      double precision zm2
      double precision zp2
!
!  Do checks.
!
      if ( m .lt. 0 ) then
        z = 0.0D+00
        return
      end if

      if ( n .lt. 0 ) then
        z = 0.0D+00
        return
      end if

      if ( n .lt. m ) then
        z = 0.0D+00
        return
      end if

      if ( mod ( n - m, 2 ) .eq. 1 ) then
        z = 0.0D+00
        return
      end if

      zm2 = 0.0D+00
      z = rho**m

      if ( m .eq. 0 ) then

        if ( n .eq. 0 ) then
          return
        end if

        zm2 = z
        z = 2.0D+00 * rho * rho - 1.0D+00

        do nn = m + 2, n - 2, 2

          a = dble ( nn + 2 ) / dble ( ( nn + 2 )**2 - m**2 )
          b = dble ( 4 * ( nn + 1 ) )
          c = dble ( ( nn + m )**2 ) / dble ( nn ) 
     &      + dble ( ( nn - m + 2 )**2 ) / dble ( nn + 2 )
          d = dble ( nn**2 - m**2 ) / dble ( nn )

          zp2 = a * ( ( b * rho * rho - c ) * z - d * zm2 ) 
          zm2 = z
          z = zp2

        end do

      else

        do nn = m, n-2, 2

          a = dble ( nn + 2 ) / dble ( ( nn + 2 )**2 - m**2 )
          b = dble ( 4 * ( nn + 1 ) )
          c = dble ( ( nn + m )**2 ) / dble ( nn ) 
     &      + dble ( ( nn - m + 2 )**2 ) / dble ( nn + 2 )
          d = dble ( nn**2 - m**2 ) / dble ( nn )

          zp2 = a * ( ( b * rho * rho - c ) * z - d * zm2 ) 
          zm2 = z
          z = zp2

        end do

      end if

      return
      end
      subroutine zernike_poly_coef ( m, n, c )

c*********************************************************************72
c
cc ZERNIKE_POLY_COEF: coefficients of a Zernike polynomial.
c
c  Discussion:
c
c    With our coefficients stored in C(0:N), the
c    radial function R^M_N(RHO) is given by
c
c      R^M_N(RHO) = C(0) 
c                 + C(1) * RHO
c                 + C(2) * RHO^2
c                 + ...
c                 + C(N) * RHO^N
c
c    and the odd and even Zernike polynomials are
c
c      Z^M_N(RHO,PHI,odd)  = R^M_N(RHO) * sin(PHI)
c      Z^M_N(RHO,PHI,even) = R^M_N(RHO) * cos(PHI)
c
c    The first few "interesting" values of R are:
c
c    R^0_0 = 1
c
c    R^1_1 = RHO
c
c    R^0_2 = 2 * RHO^2 - 1
c    R^2_2 =     RHO^2
c
c    R^1_3 = 3 * RHO^3 - 2 * RHO
c    R^3_3 =     RHO^3
c
c    R^0_4 = 6 * RHO^4 - 6 * RHO^2 + 1
c    R^2_4 = 4 * RHO^4 - 3 * RHO^2
c    R^4_4 =     RHO^4
c
c    R^1_5 = 10 * RHO^5 - 12 * RHO^3 + 3 * RHO
c    R^3_5 =  5 * RHO^5 -  4 * RHO^3
c    R^5_5 =      RHO^5
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    18 July 2008
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Eric Weisstein,
c    CRC Concise Encyclopedia of Mathematics,
c    CRC Press, 2002,
c    Second edition,
c    ISBN: 1584883472,
c    LC: QA5.W45
c
c  Parameters:
c
c    Input, integer M, N, the parameters of the polynomial.
c    Normally, 0 <= M <= N and 0 <= N.
c  
c    Output, double precision C(0:N), the coefficients of the polynomial.
c
      implicit none

      integer n

      double precision c(0:n)
      integer i
      integer l
      integer m
      integer nm_minus
      integer nm_plus
      double precision r8_choose

      do i = 0, n
        c(i) = 0.0D+00
      end do

      if ( n .lt. 0 ) then
        return
      end if

      if ( m .lt. 0 ) then
        return
      end if
          
      if ( n .lt. m ) then
        return
      end if

      if ( mod ( n - m, 2 ) .eq. 1 ) then
        return
      end if

      nm_plus = ( m + n ) / 2
      nm_minus = ( n - m ) / 2

      c(n) = r8_choose ( n, nm_plus )

      do l = 0, nm_minus - 1

        c(n-2*l-2) = - dble ( ( nm_plus - l ) * ( nm_minus - l ) ) 
     &    * c(n-2*l) / dble ( ( n - l ) * ( l + 1 ) )

      end do

      return
      end
      function zeta ( p )

c*********************************************************************72
c
cc ZETA estimates the Riemann Zeta function.
c
c  Discussion:
c
c    For 1 < P, the Riemann Zeta function is defined as:
c
c      ZETA ( P ) = Sum ( 1 <= N < +oo ) 1 / N^P
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    03 July 2008
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Daniel Zwillinger, editor,
c    CRC Standard Mathematical Tables and Formulae,
c    30th Edition,
c    CRC Press, 1996.
c
c  Parameters:
c
c    Input, double precision P, the power to which the integers are raised.
c    P must be greater than 1.
c
c    Output, double precision ZETA, an approximation to the Riemann 
c    Zeta function.
c
      implicit none

      integer n
      double precision p
      double precision total
      double precision total_old
      double precision zeta

      if ( p .le. 1.0D+00 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'ZETA - Fatal error!'
        write ( *, '(a)' ) '  Exponent P <= 1.0.'
        zeta = -1.0D+00
        stop
      end if

      total = 0.0D+00
      n = 0

10    continue

        n = n + 1
        total_old = total
        total = total + 1.0D+00 / ( dble ( n ) )**p

        if ( total .le. total_old .or. 1000 .le. n ) then
          go to 20
        end if

      go to 10

20    continue

      zeta = total

      return
      end
      subroutine zeta_values ( n_data, n, zeta )

c*********************************************************************72
c
cc ZETA_VALUES returns some values of the Riemann Zeta function.
c
c  Discussion:
c
c    ZETA(N) = sum ( 1 <= I .lt. +oo ) 1 / I**N
c
c    In Mathematica, the function can be evaluated by:
c
c      Zeta[n]
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    24 March 2007
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
c    Output, integer N, the argument of the Zeta function.
c
c    Output, double precision ZETA, the value of the Zeta function.
c
      implicit none

      integer n_max
      parameter ( n_max = 15 )

      integer n
      integer n_data
      integer n_vec(n_max)
      double precision zeta
      double precision zeta_vec(n_max)

      save n_vec
      save zeta_vec

      data n_vec /
     &   2,
     &   3,
     &   4,
     &   5,
     &   6,
     &   7,
     &   8,
     &   9,
     &  10,
     &  11,
     &  12,
     &  16,
     &  20,
     &  30,
     &  40 /
      data zeta_vec /
     &  0.164493406684822643647D+01,
     &  0.120205690315959428540D+01,
     &  0.108232323371113819152D+01,
     &  0.103692775514336992633D+01,
     &  0.101734306198444913971D+01,
     &  0.100834927738192282684D+01,
     &  0.100407735619794433939D+01,
     &  0.100200839292608221442D+01,
     &  0.100099457512781808534D+01,
     &  0.100049418860411946456D+01,
     &  0.100024608655330804830D+01,
     &  0.100001528225940865187D+01,
     &  0.100000095396203387280D+01,
     &  0.100000000093132743242D+01,
     &  0.100000000000090949478D+01 /

      if ( n_data .lt. 0 ) then
        n_data = 0
      end if

      n_data = n_data + 1

      if ( n_max .lt. n_data ) then
        n_data = 0
        n = 0
        zeta = 0.0D+00
      else
        n = n_vec(n_data)
        zeta = zeta_vec(n_data)
      end if

      return
      end
