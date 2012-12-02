      function randlc ( x )

c*****************************************************************************80
c
cc RANDLC returns a uniform pseudorandom value.
c
c  Discussion:
c
c    The number returned is in the range (0, 1).  
c
c    The algorithm uses the linear congruential generator:
c
c      X(K+1) = A * X(K)  mod 2^46
c
c    This scheme generates 2^44 numbers before repeating.  
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    08 March 2010
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    David Bailey, Eric Barszcz, John Barton, D Browning, Robert Carter, 
c    Leonardo Dagum, Rod Fatoohi,
c    Samuel Fineberg, Paul Frederickson, Thomas Lasinski, Robert Schreiber, 
c    Horst Simon, V Venkatakrishnan, Sisira Weeratunga,
c    The NAS Parallel Benchmarks,
c    RNR Technical Report RNR-94-007,
c    March 1994.
c
c    Donald Knuth,
c    The Art of Computer Programming,
c    Volume 2, Seminumerical Algorithms,
c    Third Edition,
c    Addison Wesley, 1997,
c    ISBN: 0201896842,
c    LC: QA76.6.K64.
c
c  Parameters:
c
c    Input/output, double precision X, the seed.  X should be an 
c    odd integer such that 1 <= X <= 2^46.
c
c    Output, double precision RANDLC, the next pseudorandom value.
c
      implicit none

      double precision a
      double precision a1
      double precision a2
      integer i
      integer j
      integer ks
      double precision r23
      double precision r46
      double precision randlc
      double precision t1
      double precision t2
      double precision t23
      double precision t3
      double precision t4
      double precision t46
      double precision x
      double precision x1
      double precision x2
      double precision z

      save a
      save a1
      save a2
      save ks
      save r23
      save r46
      save t23
      save t46

      data a / 1220703125.0D+00 /
      data ks / 0 /
c
c  If this is the first call, compute 
c
c    R23 = 2 ^ -23, 
c    R46 = 2 ^ -46,
c    T23 = 2 ^ 23, 
c    T46 = 2 ^ 46.  
c
c  These are computed in loops, rather than by merely using the power operator, 
c  in order to insure that the results are exact on all systems.  
c
      if ( ks == 0 ) then

        r23 = 1.0D+00
        r46 = 1.0D+00
        t23 = 1.0D+00
        t46 = 1.0D+00

        do i = 1, 23
          r23 = 0.5D+00 * r23
          t23 = 2.0D+00 * t23
        end do

        do i = 1, 46
          r46 = 0.50D+00 * r46
          t46 = 2.0D+00 * t46
        end do
c
c  Break A into two parts such that A = 2^23 * A1 + A2.
c
        t1 = r23 * a
        a1 = dble ( int ( t1 ) )
        a2 = a - t23 * a1

        ks = 1

      end if
c
c  Deal with a 0 input value of X.
c
      if ( x .eq. 0.0D+00 ) then
        x = 314159265.0D+00
      end if
c
c  Deal somewhat arbitrarily with negative input X.
c
      if ( x < 0.0D+00 ) then
        x = - x
      end if
c
c  Break X into two parts X1 and X2 such that:
c
c    X = 2^23 * X1 + X2, 
c
c  then compute
c
c    Z = A1 * X2 + A2 * X1  (mod 2^23)
c    X = 2^23 * Z + A2 * X2  (mod 2^46).
c
      t1 = r23 * x
      x1 = dble ( int ( t1 ) )
      x2 = x - t23 * x1

      t1 = a1 * x2 + a2 * x1
      t2 = dble ( int ( r23 * t1 ) )
      z = t1 - t23 * t2

      t3 = t23 * z + a2 * x2
      t4 = dble ( int ( r46 * t3 ) )
      x = t3 - t46 * t4

      randlc = r46 * x

      return
      end
      function randlc_jump ( x, k )

c*********************************************************************72
c
cc RANDLC_JUMP returns the K-th element of a uniform pseudorandom sequence.
c
c  Discussion:
c
c    The sequence uses the linear congruential generator:
c
c      X(K+1) = A * X(K)  mod 2^46
c
c    The K-th element, which can be represented as
c
c      X(K) = A^K * X(0)  mod 2^46
c
c    is computed directly using the binary algorithm for exponentiation.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    11 March 2010
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    David Bailey, Eric Barszcz, John Barton, D Browning, Robert Carter, 
c    Leonardo Dagum, Rod Fatoohi,
c    Samuel Fineberg, Paul Frederickson, Thomas Lasinski, Robert Schreiber, 
c    Horst Simon, V Venkatakrishnan, Sisira Weeratunga,
c    The NAS Parallel Benchmarks,
c    RNR Technical Report RNR-94-007,
c    March 1994.
c
c    Donald Knuth,
c    The Art of Computer Programming,
c    Volume 2, Seminumerical Algorithms,
c    Third Edition,
c    Addison Wesley, 1997,
c    ISBN: 0201896842,
c    LC: QA76.6.K64.
c
c  Parameters:
c
c    Input, double precision X, the initial seed (with index 0).  
c
c    Input, integer K, the index of the desired value.
c
c    Output, double precision RANDLC_JUMP, the K-th value in the sequence.
c
      implicit none

      double precision a 
      parameter ( a = 1220703125.0D+00 )
      double precision a1
      save a1
      double precision a2
      save a2
      double precision b
      double precision b1
      double precision b2
      integer i
      integer j
      integer k
      integer k2
      integer ks
      save ks
      integer m
      double precision r23
      save r23
      double precision r46
      save r46
      double precision randlc_jump
      double precision t1
      double precision t2
      double precision t23
      save t23
      double precision t3
      double precision t4
      double precision t46
      save t46
      integer twom
      double precision x
      double precision x1
      double precision x2
      double precision xk
      double precision z

      data ks / 0 /
c
c  If this is the first call, compute 
c
c    R23 = 2 ^ -23, 
c    R46 = 2 ^ -46,
c    T23 = 2 ^ 23, 
c    T46 = 2 ^ 46.  
c
c  These are computed in loops, rather than by merely using the power operator, 
c  in order to insure that the results are exact on all systems.  
c
      if ( ks == 0 ) then

        r23 = 1.0D+00
        r46 = 1.0D+00
        t23 = 1.0D+00
        t46 = 1.0D+00

        do i = 1, 23
          r23 = 0.5D+00 * r23
          t23 = 2.0D+00 * t23
        end do

        do i = 1, 46
          r46 = 0.50D+00 * r46
          t46 = 2.0D+00 * t46
        end do
c
c  Break A into two parts such that A = 2^23 * A1 + A2.
c
        t1 = r23 * a
        a1 = real ( int ( t1 ), kind = 8 )
        a2 = a - t23 * a1

        ks = 1

      end if

      if ( k < 0 ) then

        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'RANDLC_JUMP - Fatal errorc'
        write ( *, '(a)' ) '  K < 0.'
        stop

      else if ( k == 0 ) then

        xk = x
c
c  Find M so that K < 2^M.
c
      else

        k2 = k
        xk = x

        m = 1
        twom = 2

10      continue

        if ( twom <= k ) then
          twom = twom * 2
          m = m + 1
          go to 10
        end if

        b = a
        b1 = a1
        b2 = a2

        do i = 1, m

          j = k2 / 2
c
c  Replace X by A * X, if appropriate.
c
          if ( 2 * j .ne. k2 ) then

            t1 = r23 * xk
            x1 = dble ( int ( t1 ) )
            x2 = xk - t23 * x1

            t1 = b1 * x2 + b2 * x1
            t2 = dble ( int ( r23 * t1 ) )
            z = t1 - t23 * t2

            t3 = t23 * z + b2 * x2
            t4 = dble ( int ( r46 * t3 ) )
            xk = t3 - t46 * t4

          end if
c
c  Replace A by A * A mod 2^46.
c
          t1 = r23 * b
          x1 = dble ( int ( t1 ) )
          x2 = b - t23 * x1

          t1 = b1 * x2 + b2 * x1
          t2 = dble ( int ( r23 * t1 ) )
          z = t1 - t23 * t2

          t3 = t23 * z + b2 * x2
          t4 = dble ( int ( r46 * t3 ) )
          b = t3 - t46 * t4
c
c  Update A1, A2.
c
          t1 = r23 * b
          b1 = dble ( int ( t1 ) )
          b2 = b - t23 * b1

          k2 = j

        end do

      end if

      randlc_jump = r46 * xk

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
