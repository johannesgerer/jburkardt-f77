      program main

c*********************************************************************72
c
cc MAIN is the main program for RANDOM_MPI.
c
c  Discussion:
c
c    This program demonstrates how P processors can generate the same
c    sequence of random numbers as 1 processor.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    22 May 2008
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    William Gropp, Ewing Lusk, Anthony Skjellum,
c    Using MPI: Portable Parallel Programming with the
c    Message-Passing Interface,
c    Second Edition,
c    MIT Press, 1999,
c    ISBN: 0262571323.
c
      include 'mpif.h'

      integer a
      integer an
      integer b
      integer bn
      integer c
      integer error
      integer id
      integer j
      integer k
      integer p
      integer u
      integer v
c
c  Initialize MPI.
c
      call MPI_Init ( error )
c
c  Get the number of processes.
c
      call MPI_Comm_size ( MPI_COMM_WORLD, p, error )
c
c  Get the rank of this processor.
c
      call MPI_Comm_rank ( MPI_COMM_WORLD, id, error )
c
c  Print a message.
c
      if ( id == 0 ) then
        call timestamp ( )
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'RANDOM_MPI - Master process:'
        write ( *, '(a)' ) '  FORTRAN77 version'
        write ( *, '(a,i8)' ) '  The number of processors is P = ', p
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 
     &  '  This program shows how a stream of random numbers'
        write ( *, '(a)' ) 
     &  '  can be computed "in parallel" in an MPI program.'
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 
     &  '  We assume we are using a linear congruential'
        write ( *, '(a)' ) 
     &  '  random number generator or "LCRG", which takes'
        write ( *, '(a)' ) 
     &  '  an integer input and returns a new integer output:'
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) '    U = ( A * V + B ) mod C'
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 
     &  '  We assume that we want the MPI program to produce'
        write ( *, '(a)' ) 
     &  '  the same sequence of random values as a sequential'
        write ( *, '(a)' ) 
     &  '  program would - but we want each processor to compute'
        write ( *, '(a)' ) '  one part of that sequence.'
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 
     &  '  We do this by computing a new LCRG which can compute'
        write ( *, '(a)' ) '  every P''th entry of the original one.'
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 
     &  '  Our LCRG works with integers, but it is easy to'
        write ( *, '(a)' ) 
     &    '  turn each integer into a real number between [0,1].'
      end if
c
c  A, B and C define the linear congruential random number generator.
c
      a = 16807
      b = 0
      c = 2147483647

      if ( id == 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) '  LCRG parameters:'
        write ( *, '(a)' ) ' '
        write ( *, '(a,i12)' ) '  A  = ', a
        write ( *, '(a,i12)' ) '  B  = ', b
        write ( *, '(a,i12)' ) '  C  = ', c
      end if

      k_hi = p * 10
c
c  Processor 0 generates 10 * P random values.
c
      if ( id == 0 ) then

        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 
     &    '  Let processor 0 generate the entire sequence.'
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) '     K    ID         Input        Output'
        write ( *, '(a)' ) ' '

        k = 0
        v = 12345
        write ( *, '(2x,i4,2x,i4,2x,12x,2x,i12)' ) k, id, v

        do k = 1, k_hi
          u = v
          call lcrg_evaluate ( a, b, c, u, v )
          write ( *, '(2x,i4,2x,i4,2x,i12,2x,i12)' ) k, id, u, v
        end do

      end if
c
c  Processor P now participates by computing the P-th part of the sequence.
c
      call lcrg_anbn ( a, b, c, p, an, bn )

      if ( id == 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) '  LCRG parameters for P processors:'
        write ( *, '(a)' ) ' '
        write ( *, '(a,i12)' ) '  AN = ', an
        write ( *, '(a,i12)' ) '  BN = ', bn
        write ( *, '(a,i12)' ) '  C  = ', c
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 
     &    '  Have ALL the processors participate in computing'
        write ( *, '(a)' ) '  the same random number sequence.'
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) '     K    ID         Input        Output'
        write ( *, '(a)' ) ' '
      end if
c
c  Use the basis LCRG to get the ID-th value in the sequence.
c
      v = 12345
      do j = 1, id
        u = v
        call lcrg_evaluate ( a, b, c, u, v )
      end do
      k = id

      write ( *, '(2x,i4,2x,i4,2x,12x,2x,i12)' ) k, id, v
c
c  Now use the "skipping" LCRG to compute the values with indices
c  ID, ID+P, ID+2P, ...,
c
      do k = id + p, k_hi, p
        u = v
        call lcrg_evaluate ( an, bn, c, u, v )
        write ( *, '(2x,i4,2x,i4,2x,i12,2x,i12)' ) k, id, u, v
      end do
c
c  Terminate MPI.
c
      call MPI_Finalize ( error )
c
c  Terminate.
c
      if ( id == 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'RANDOM_MPI - Master process:'
        write ( *, '(a)' ) '  Normal end of execution.'
        write ( *, '(a)' ) ' '
        call timestamp ( )
      end if

      stop
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
      else if ( a .ne. 0 .and. b .eq. 0 .and. c /= 0 ) then
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
      subroutine lcrg_anbn ( a, b, c, n, an, bn )

c*********************************************************************72
c
cc LCRG_ANBN computes the "N-th power" of a linear congruential generator.
c
c  Discussion:
c
c    We are considering a linear congruential random number generator.
c    The LCRG takes as input an integer value called SEED, and returns
c    an updated value of SEED, 
c
c      SEED(out) = ( a * SEED(in) + b ) mod c.
c
c    and an associated pseudorandom real value
c
c      U = SEED(out) / c.
c
c    In most cases, a user is content to call the LCRG repeatedly, with
c    the updating of SEED being taken care of automatically.
c
c    The purpose of this routine is to determine the values of AN and BN
c    that describe the LCRG that is equivalent to N applications of the
c    original LCRG.
c
c    One use for such a facility would be to do random number computations
c    in parallel.  If each of N processors is to compute many random values, 
c    you can guarantee that they work with distinct random values 
c    by starting with a single value of SEED, using the original LCRG to generate
c    the first N-1 "iterates" of SEED, so that you now have N "seed" values,
c    and from now on, applying the N-th power of the LCRG to the seeds.
c
c    If the K-th processor starts from the K-th seed, it will essentially
c    be computing every N-th entry of the original random number sequence,
c    offset by K.  Thus the individual processors will be using a random
c    number stream as good as the original one, and without repeating, and
c    without having to communicate.
c 
c    To evaluate the N-th value of SEED directly, we start by ignoring 
c    the modular arithmetic, and working out the sequence of calculations
c    as follows:
c
c      SEED(0)   =     SEED.
c      SEED(1)   = a * SEED      + b
c      SEED(2)   = a * SEED(1)   + b = a^2 * SEED           + a * b + b
c      SEED(3)   = a * SEED(2)   + b = a^3 * SEED + a^2 * b + a * b + b
c      ...
c      SEED(N-1) = a * SEED(N-2) + b 
c
c      SEED(N) = a * SEED(N-1) + b = a^N * SEED 
c                                    + ( a^(n-1) + a^(n-2) + ... + a + 1 ) * b
c
c    or, using the geometric series,
c
c      SEED(N) = a^N * SEED + ( a^N - 1) / ( a - 1 ) * b
c              = AN * SEED + BN
c
c    Thus, from any SEED, we can determine the result of N applications of the
c    original LCRG directly if we can solve
c
c      ( a - 1 ) * BN = ( a^N - 1 ) * b in modular arithmetic,
c
c    and evaluate:
c
c      AN = a^N
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    21 May 2008
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Barry Wilkinson, Michael Allen,
c    Parallel Programming:
c    Techniques and Applications Using Networked Workstations and Parallel Computers,
c    Prentice Hall, 
c    ISBN: 0-13-140563-2,
c    LC: QA76.642.W54.
c
c  Parameters:
c
c    Input, integer A, the multiplier for the LCRG.
c
c    Input, integer B, the added value for the LCRG.
c
c    Input, integer C, the base for the modular arithmetic.  
c    For 32 bit arithmetic, this is often 2^31 - 1, or 2147483647.  It is
c    required that 0 < C.
c
c    Input, integer N, the "index", or number of times that the
c    LCRG is to be applied.  It is required that 0 <= N.
c
c    Output, integer AN, BN, the multiplier and added value for
c    the LCRG that represent N applications of the original LCRG.
c
      implicit none

      integer a
      integer am1
      integer an
      integer anm1tb
      integer b
      integer bn
      integer c
      integer ierror
      integer n

      if ( n .lt. 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'LCRG_ANBN - Fatal error!'
        write ( *, '(a,i12)' ) '  Illegal input value of N = ', n
        stop
      end if

      if ( c .le. 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'LCRG_ANBN - Fatal error!'
        write ( *, '(a,i12)' ) '  Illegal input value of C = ', c
        stop
      end if

      if ( n .eq. 0 ) then
        an = 1
        bn = 0
      else if ( n .eq. 1 ) then
        an = a
        bn = b
      else
c
c  Compute A^N.
c
        call power_mod ( a, n, c, an )
c
c  Solve 
c    ( a - 1 ) * BN = ( a^N - 1 ) mod B
c  for BN.
c
        am1 = a - 1
        anm1tb = ( an - 1 ) * b

        call congruence ( am1, c, anm1tb, ierror, bn )

        if ( ierror .ne. 0 ) then
          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) 'LCRG_ANBN - Fatal error!'
          write ( *, '(a)' ) '  An error occurred in CONGRUENCE.'
          write ( *, '(a,i8)' ) '  The error code was IERROR = ', ierror
          stop
        end if

      end if

      return
      end
      subroutine lcrg_evaluate ( a, b, c, x, y )

c*********************************************************************72
c
cc LCRG_EVALUATE evaluates an LCRG, y = ( A * x + B ) mod C.
c
c  Discussion:
c
c    This routine cannot be recommended for production use.  Because we want
c    to do modular arithmetic, but the base is not a power of 2, we need to
c    use "double precision" integers to keep accuracy.
c
c    If we knew the base C, we could try to avoid overflow while not changing
c    precision.
c
c    If the base C was a power of 2, we could rely on the usual properties of 
c    integer arithmetic on computers, in which overflow bits, which are always 
c    ignored, don't actually matter.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    21 May 2008
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer A, the multiplier for the LCRG.
c
c    Input, integer B, the added value for the LCRG.
c
c    Input, integer C, the base for the modular arithmetic.  
c    For 32 bit arithmetic, this is often 2^31 - 1, or 2147483647.  It is
c    required that 0 < C.
c
c    Input, integer X, the value to be processed.
c
c    Output, integer Y, the processed value.
c
      implicit none

      integer a
      integer*8 a8
      integer b
      integer*8 b8
      integer c
      integer*8 c8
      integer x
      integer*8 x8
      integer y
      integer*8 y8
c
c  To avoid roundoff issues, we need to go to "double precision" integers.
c  (Not available on all planets.)
c
      a8 = a
      b8 = b
      c8 = c
      x8 = x

      y8 = mod ( a8 * x8 + b8, c8 )

      y = int ( y8 )

      if ( y < 0 ) then
        y = y + c
      end if

      return
      end
      subroutine power_mod ( a, n, m, x )

c*********************************************************************72
c
cc POWER_MOD computes mod ( A^N, M ).
c
c  Discussion:
c
c    Some programming tricks are used to speed up the computation, and to
c    allow computations in which the value A**N is much too large to
c    store in an integer word.
c
c    First, for efficiency, the power A**N is computed by determining
c    the binary expansion of N, then computing A, A^2, A^4, and so on
c    by repeated squaring, and multiplying only those factors that
c    contribute to A**N.
c
c    Secondly, the intermediate products are immediately "mod'ed", which
c    keeps them small.
c
c    For instance, to compute mod ( A^13, 11 ), we essentially compute
c
c       13 = 1 + 4 + 8
c
c       A^13 = A * A^4 * A^8
c
c       mod ( A^13, 11 ) = mod ( A, 11 ) * mod ( A^4, 11 ) * mod ( A^8, 11 ).
c
c    Fermat's little theorem says that if P is prime, and A is not divisible
c    by P, then ( A^(P-1) - 1 ) is divisible by P.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    22 May 2008
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer A, the base of the expression to be tested.
c    0 <= A is required.
c
c    Input, integer N, the power to which the base is raised.
c    0 <= N is required.
c
c    Input, integer M, the divisor against which the expression is tested.
c    0 < M is required.
c
c    Output, integer X, the remainder when A^N is divided by M.
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
c  A_SQUARE contains the successive squares of A.
c
      a_square2 = a
      x2 = 1
      m2 = m

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
c  Fix up X so that it is nonnegative.
c
20    continue

      if ( x2 .lt. 0 ) then
        x2 = x2 + m2
        go to 20
      end if

      x = int ( x2 )

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
