      subroutine ctof ( nc, uc, nf, uf )

c*********************************************************************72
c                                                    
cc CTOF transfers data from a coarse to a finer grid.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    25 November 2011
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    William Hager,
c    Applied Numerical Linear Algebra,
c    Prentice-Hall, 1988,
c    ISBN13: 978-0130412942,
c    LC: QA184.H33.
c
c  Parameters:
c
c    Input, integer NC, the number of coarse nodes.
c
c    Input, double precision UC(NC), the coarse correction data.
c
c    Input, integer NF, the number of fine nodes.
c
c    Input/output, double precision UF(NF), on input, the fine grid data.
c    On output, the data has been updated with prolonged coarse 
c    correction data.
c
      implicit none

      integer nc
      integer nf

      integer ic
      integer if
      double precision uc(nc)
      double precision uf(nf)

      do ic = 1, nc
        if = 2 * ic - 1
        uf(if) = uf(if) + uc(ic)
      end do

      do ic = 1, nc - 1
        if = 2 * ic
        uf(if) = uf(if) + 0.5D+00 * ( uc(ic) + uc(ic+1) )
      end do

      return
      end
      subroutine ftoc ( nf, uf, rf, nc, uc, rc )

c*********************************************************************72
c                                                    
cc FTOC transfers data from a fine grid to a coarser grid.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    07 December 2011
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    William Hager,
c    Applied Numerical Linear Algebra,
c    Prentice-Hall, 1988,
c    ISBN13: 978-0130412942,
c    LC: QA184.H33.
c
c  Parameters:
c
c    Input, integer NF, the number of fine nodes.
c
c    Input, double precision UF(NF), the fine data.
c
c    Input, double precision RF(NF), the right hand side for the fine grid.
c
c    Input, integer NC, the number of coarse nodes.
c
c    Output, double precision UC(NC), the coarse grid data, set to zero.
c
c    Output, double precision RC(NC), the right hand side for the coarse grid.
c
      implicit none

      integer nc
      integer nf

      integer ic
      integer if
      double precision rc(nc)
      double precision rf(nf)
      double precision uc(nc)
      double precision uf(nf)

      uc(1) = 0.0D+00
c     rc(1) = rf(1) - uf(1)
      rc(1) = 0.0D+00
      do ic = 2, nc - 1
        if = 2 * ic - 1
        uc(ic) = 0.0D+00
        rc(ic) = 4.0D+00 
     &    * ( rf(if) + uf(if-1) - 2.0D+00 * uf(if) + uf(if+1) )
      end do
      uc(nc) = 0.0D+00
c     rc(nc) = rf(nc) - uf(nc)
      rc(nc) = 0.0D+00

      return
      end
      subroutine gauss_seidel ( n, r, u, dif_l1 )

c*********************************************************************72
c                                                    
cc GAUSS_SEIDEL carries out one step of a Gauss-Seidel iteration.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    25 November 2011
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    William Hager,
c    Applied Numerical Linear Algebra,
c    Prentice-Hall, 1988,
c    ISBN13: 978-0130412942,
c    LC: QA184.H33.
c
c  Parameters:
c
c    Input, integer N, the number of unknowns.
c
c    Input, double precision R(N), the right hand side.
c
c    Input/output, double precision U(N), the estimated solution.
c
c    Output, double precision DIF_L1, the L1 norm of the difference between the
c    input and output solution estimates.
c
      implicit none

      integer n

      double precision dif_l1
      integer i
      double precision r(n)
      double precision u(n)
      double precision u_old

      dif_l1 = 0.0D+00
c
c  Setting U(1) = R(1), U(N) = R(N) seems right, but it makes the
c  code fail.
c
c     u(1) = r(1)
      do i = 2, n - 1
        u_old = u(i)
        u(i) = 0.5D+00 * ( u(i-1) + u(i+1) + r(i) )
        dif_l1 = dif_l1 + abs ( u(i) - u_old )
      end do
c     u(n) = r(n)

      return
      end
      function i4_log_2 ( i )

c*********************************************************************72
c
cc I4_LOG_2 returns the integer part of the logarithm base 2 of an I4.
c
c  Discussion:
c
c    For positive I4_LOG_2(I), it should be true that
c      2^I4_LOG_2(X) .le. |I| < 2^(I4_LOG_2(I)+1).
c    The special case of I4_LOG_2(0) returns -HUGE().
c
c    An I4 is an integer value.
c
c  Example:
c
c     I  I4_LOG_2
c
c     0  -1
c     1,  0
c     2,  1
c     3,  1
c     4,  2
c     5,  2
c     6,  2
c     7,  2
c     8,  3
c     9,  3
c    10,  3
c   127,  6
c   128,  7
c   129,  7
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    28 October 2007
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer I, the number whose logarithm base 2
c    is desired.
c
c    Output, integer I4_LOG_2, the integer part of the
c    logarithm base 2 of the absolute value of I.
c
      implicit none

      integer i
      integer i_abs
      integer i4_log_2
      integer i4_huge
      parameter ( i4_huge = 2147483647 )

      if ( i .eq. 0 ) then

        i4_log_2 = - i4_huge

      else

        i4_log_2 = 0

        i_abs = abs ( i )

10      continue

        if ( 2 .le. i_abs ) then
          i_abs = i_abs / 2
          i4_log_2 = i4_log_2 + 1
          go to 10
        end if

      end if

      return
      end
      subroutine monogrid_poisson_1d ( n, force, exact, it_num, u )

c*********************************************************************72
c                                                    
cc MONOGRID_POISSON_1D solves a 1D PDE, using the Gauss-Seidel method.
c
c  Discussion:
c
c    This routine solves a 1D boundary value problem of the form
c
c      - U''(X) = F(X) for A < X < B,
c
c    with boundary conditions U(A) = UA, U(B) = UB.
c
c    The Gauss-Seidel method is used. 
c
c    This routine is provided primarily for comparison with the
c    multigrid solver.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    25 November 2011
c
c  Author:
c
c    Original FORTRAN77 version by William Hager.
c    This FORTRAN77 version by John Burkardt.
c
c  Reference:
c
c    William Hager,
c    Applied Numerical Linear Algebra,
c    Prentice-Hall, 1988,
c    ISBN13: 978-0130412942,
c    LC: QA184.H33.
c
c  Parameters:
c
c    Input, integer N, the number of intervals.
c
c    Input, double precision external FORCE, the name of the function 
c    which evaluates the right hand side.
c
c    Input, double precision external EXACT, the name of the function which 
c    evaluates the exact solution.
c
c    Output, integer IT_NUM, the number of iterations.
c
c    Output, double precision U(N+1), the computed solution.
c
      implicit none

      integer n

      double precision d1
      double precision difmax
      external exact
      double precision exact
      external force
      double precision force
      double precision h
      integer i
      integer it_num
      integer k
      double precision r(n+1)
      double precision tol
      double precision u(n+1)
      double precision x
c
c  Initialization.
c
      tol = 0.0001D+00

      r(1) = 0.0D+00
      h = 1.0D+00 / dble ( n )
      do i = 2, n
        x = dble ( i - 1 ) / dble ( n )
        r(i) = h**2 * force ( x )
      end do
      r(n+1) = 0.0D+00

      do i = 1, n + 1
        u(i) = 0.0D+00
      end do

      it_num = 0
c
c  Gauss-Seidel iteration.
c
10    continue

      it_num = it_num + 1

      call gauss_seidel ( n + 1, r, u, d1 )

      if ( tol < d1 ) then
        go to 10
      end if

      return
      end
      subroutine multigrid_poisson_1d ( n, force, exact, it_num, u )

c*********************************************************************72
c                                                    
cc MULTIGRID_POISSON_1D solves a 1D PDE using the multigrid method.
c
c  Discussion:
c
c    This routine solves a 1D boundary value problem of the form
c
c      - U''(X) = F(X) for A < X < B,
c
c    with boundary conditions U(A) = UA, U(B) = UB.
c
c    The multigrid method is used. 
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    27 November 2011
c
c  Author:
c
c    Original FORTRAN77 version by William Hager.
c    This FORTRAN77 version by John Burkardt.
c
c  Reference:
c
c    William Hager,
c    Applied Numerical Linear Algebra,
c    Prentice-Hall, 1988,
c    ISBN13: 978-0130412942,
c    LC: QA184.H33.
c
c  Parameters:
c
c    Input, integer N, the number of intervals.
c    N must be a power of 2.
c
c    Input, double precision external FORCE, the name of the function 
c    which evaluates the right hand side.
c
c    Input, double precision external EXACT, the name of the function 
c    which evaluates the exact solution.
c
c    Output, integer IT_NUM, the number of iterations.
c
c    Output, integer U(N+1), the computed solution.
c
      implicit none

      integer n
      integer nl_max
      parameter ( nl_max = 518 )

      double precision d0
      double precision d1
      double precision difmax
      external exact
      double precision exact
      external force
      double precision force
      double precision h
      integer i
      integer i4_log_2
      integer it
      integer it_num
      integer j
      integer k
      integer l
      integer ll
      integer m
      integer nl
      double precision r(nl_max)
      double precision s
      double precision tol
      double precision u(n+1)
      double precision utol
      double precision uu(nl_max)
      double precision x
c
c  Determine if we have enough storage.
c
      k = i4_log_2 ( n )

      if ( n .ne. 2**k ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'MULTIGRID_POISSON_1D - Fatal error!'
        write ( *, '(a)' ) '  Interval number N is not a power of 2.'
        stop
      end if

      nl = n + n + k - 2

      if ( nl_max < nl ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'MULTIGRID_POISSON_1D - Fatal error!'
        write ( *, '(a)' ) '  Grid index K requires NL = ', nl
        write ( *, '(a)' ) '  Internal parameter NL_MAX = ', nl_max
        stop
      end if
c
c  Initialization.
c
      it = 4
      it_num = 0
      tol = 0.0001D+00
      utol = 0.7D+00
      m = n
c 
c  Set the right hand side.
c 
      r(1) = 0.0D+00
      h = 1.0D+00 / dble ( n )
      do i = 2, n
        x = dble ( i - 1 ) / dble ( n )
        r(i) = h**2 * force ( x )
      end do
      r(n+1) = 0.0D+00

      do i = 1, nl
        uu(i) = 0.0D+00
      end do
c
c  L points to first entry of solution
c  LL points to penultimate entry.
c
      l = 1
      ll = n
c 
c  Gauss-Seidel iteration
c
      d1 = 0.0D+00
      j = 0

10    continue

        d0 = d1
        j = j + 1
        call gauss_seidel ( n + 1, r(l), uu(l), d1 )
        it_num = it_num + 1
!
!  Must do at least 4 iterations at each level.
!
        if ( j .lt. it ) then

          go to 10
!
!  If enough iterations, and satisfactory decrease, and
!  on finest grid, exit.
!
        else if ( d1 .lt. tol .and. n .eq. m ) then

          go to 20
!
!  Enough iterations, satisfactory convergence, go finer.
!
        else if ( d1 .lt. tol ) then

          call ctof ( n + 1, uu(l), n + n + 1, uu(l-1-n-n) )

          n = n + n
          ll = l - 2
          l = l - 1 - n
          j = 0
!
!  Enough iterations, slow convergence, and 2 < n, go coarser.
!
        else if ( utol * d0 .le. d1 .and. 2 .lt. n ) then

          call ftoc ( n + 1, uu(l), r(l), (n/2)+1, uu(l+n+1), r(l+n+1) )

          n = n / 2
          l = ll + 2
          ll = ll + n + 1
          j = 0

        end if

      go to 10
c
c  Computation complete
c
20    continue

      do i = 1, n + 1
        u(i) = uu(i)
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
