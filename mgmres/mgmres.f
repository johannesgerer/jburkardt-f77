      subroutine atx_cr ( n, nz_num, ia, ja, a, x, w )

c*********************************************************************72
c
cc ATX_CR computes A'*x for a matrix stored in sparse compressed row form.
c
c  Discussion:
c
c    The Sparse Compressed Row storage format is used.
c
c    The matrix A is assumed to be sparse.  To save on storage, only
c    the nonzero entries of A are stored.  The vector JA stores the
c    column index of the nonzero value.  The nonzero values are sorted
c    by row, and the compressed row vector IA then has the property that
c    the entries in A and JA that correspond to row I occur in indices
c    IA[I] through IA[I+1]-1.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    30 August 2010
c
c  Author:
c
c    Original C version by Lili Ju.
c    FORTRAN77 version by John Burkardt.
c
c  Reference:
c
c    Richard Barrett, Michael Berry, Tony Chan, James Demmel,
c    June Donato, Jack Dongarra, Victor Eijkhout, Roidan Pozo,
c    Charles Romine, Henk van der Vorst,
c    Templates for the Solution of Linear Systems:
c    Building Blocks for Iterative Methods,
c    SIAM, 1994.
c    ISBN: 0898714710,
c    LC: QA297.8.T45.
c
c    Tim Kelley,
c    Iterative Methods for Linear and Nonlinear Equations,
c    SIAM, 2004,
c    ISBN: 0898713528,
c    LC: QA297.8.K45.
c
c    Yousef Saad,
c    Iterative Methods for Sparse Linear Systems,
c    Second Edition,
c    SIAM, 2003,
c    ISBN: 0898715342,
c    LC: QA188.S17.
c
c  Parameters:
c
c    Input, integer N, the order of the system.
c
c    Input, integer NZ_NUM, the number of nonzeros.
c
c    Input, integer IA(N+1), JA(NZ_NUM), the row and column
c    indices of the matrix values.  The row vector has been compressed.
c
c    Input, double precision A(NZ_NUM), the matrix values.
c
c    Input, double precision X(N), the vector to be multiplied by A'.
c
c    Output, double precision W(N), the value of A'*X.
c
      implicit none

      integer n
      integer nz_num

      double precision a(nz_num)
      integer i
      integer ia(n+1)
      integer ja(nz_num)
      integer k
      integer k1
      integer k2
      double precision w(n)
      double precision x(n)

      do i = 1, n
        w(i) = 0.0D+00
      end do

      do i = 1, n
        k1 = ia(i)
        k2 = ia(i+1) - 1
        do k = k1, k2
          w(ja(k)) = w(ja(k)) + a(k) * x(i)
        end do
      end do

      return
      end
      subroutine atx_st ( n, nz_num, ia, ja, a, x, w )

c*********************************************************************72
c
cc ATX_ST computes A'*x for a matrix stored in sparset triplet form.
c
c  Discussion:
c
c    The matrix A is assumed to be sparse.  To save on storage, only
c    the nonzero entries of A are stored.  For instance, the K-th nonzero
c    entry in the matrix is stored by:
c
c      A(K) = value of entry,
c      IA(K) = row of entry,
c      JA(K) = column of entry.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    30 August 2010
c
c  Author:
c
c    Original C version by Lili Ju.
c    FORTRAN77 version by John Burkardt.
c
c  Reference:
c
c    Richard Barrett, Michael Berry, Tony Chan, James Demmel,
c    June Donato, Jack Dongarra, Victor Eijkhout, Roidan Pozo,
c    Charles Romine, Henk van der Vorst,
c    Templates for the Solution of Linear Systems:
c    Building Blocks for Iterative Methods,
c    SIAM, 1994.
c    ISBN: 0898714710,
c    LC: QA297.8.T45.
c
c    Tim Kelley,
c    Iterative Methods for Linear and Nonlinear Equations,
c    SIAM, 2004,
c    ISBN: 0898713528,
c    LC: QA297.8.K45.
c
c    Yousef Saad,
c    Iterative Methods for Sparse Linear Systems,
c    Second Edition,
c    SIAM, 2003,
c    ISBN: 0898715342,
c    LC: QA188.S17.
c
c  Parameters:
c
c    Input, integer N, the order of the system.
c
c    Input, integer NZ_NUM, the number of nonzeros.
c
c    Input, integer IA(NZ_NUM), JA(NZ_NUM), the row and column
c    indices of the matrix values.
c
c    Input, double precision A(NZ_NUM), the matrix values.
c
c    Input, double precision X(N), the vector to be multiplied by A'.
c
c    Output, double precision W(N), the value of A'*X.
c
      implicit none

      integer n
      integer nz_num

      double precision a(nz_num)
      integer i
      integer ia(nz_num)
      integer j
      integer ja(nz_num)
      integer k
      double precision w(n)
      double precision x(n)

      do i = 1, n
        w(i) = 0.0D+00
      end do

      do k = 1, nz_num
        i = ia(k)
        j = ja(k)
        w(j) = w(j) + a(k) * x(i)
      end do

      return
      end
      subroutine ax_cr ( n, nz_num, ia, ja, a, x, w )

c*********************************************************************72
c
cc AX_CR computes A*x for a matrix stored in sparse compressed row form.
c
c  Discussion:
c
c    The Sparse Compressed Row storage format is used.
c
c    The matrix A is assumed to be sparse.  To save on storage, only
c    the nonzero entries of A are stored.  The vector JA stores the
c    column index of the nonzero value.  The nonzero values are sorted
c    by row, and the compressed row vector IA then has the property that
c    the entries in A and JA that correspond to row I occur in indices
c    IA[I] through IA[I+1]-1.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    30 August 2010
c
c  Author:
c
c    Original C version by Lili Ju.
c    FORTRAN77 version by John Burkardt.
c
c  Reference:
c
c    Richard Barrett, Michael Berry, Tony Chan, James Demmel,
c    June Donato, Jack Dongarra, Victor Eijkhout, Roidan Pozo,
c    Charles Romine, Henk van der Vorst,
c    Templates for the Solution of Linear Systems:
c    Building Blocks for Iterative Methods,
c    SIAM, 1994.
c    ISBN: 0898714710,
c    LC: QA297.8.T45.
c
c    Tim Kelley,
c    Iterative Methods for Linear and Nonlinear Equations,
c    SIAM, 2004,
c    ISBN: 0898713528,
c    LC: QA297.8.K45.
c
c    Yousef Saad,
c    Iterative Methods for Sparse Linear Systems,
c    Second Edition,
c    SIAM, 2003,
c    ISBN: 0898715342,
c    LC: QA188.S17.
c
c  Parameters:
c
c    Input, integer N, the order of the system.
c
c    Input, integer NZ_NUM, the number of nonzeros.
c
c    Input, integer IA(N+1), JA(NZ_NUM), the row and column
c    indices of the matrix values.  The row vector has been compressed.
c
c    Input, double precision A(NZ_NUM), the matrix values.
c
c    Input, double precision X(N), the vector to be multiplied by A.
c
c    Output, double precision W(N), the value of A*X.
c
      implicit none

      integer n
      integer nz_num

      double precision a(nz_num)
      integer i
      integer ia(n+1)
      integer ja(nz_num)
      integer k
      integer k1
      integer k2
      double precision w(n)
      double precision x(n)

      do i = 1, n
        w(i) = 0.0D+00
      end do

      do i = 1, n
        k1 = ia(i)
        k2 = ia(i+1) - 1
        do k = k1, k2
          w(i) = w(i) + a(k) * x(ja(k))
        end do
      end do

      return
      end
      subroutine ax_st ( n, nz_num, ia, ja, a, x, w )

c*********************************************************************72
c
cc AX_ST computes A*x for a matrix stored in sparset triplet form.
c
c  Discussion:
c
c    The matrix A is assumed to be sparse.  To save on storage, only
c    the nonzero entries of A are stored.  For instance, the K-th nonzero
c    entry in the matrix is stored by:
c
c      A(K) = value of entry,
c      IA(K) = row of entry,
c      JA(K) = column of entry.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    30 August 2010
c
c  Author:
c
c    Original C version by Lili Ju.
c    FORTRAN77 version by John Burkardt.
c
c  Reference:
c
c    Richard Barrett, Michael Berry, Tony Chan, James Demmel,
c    June Donato, Jack Dongarra, Victor Eijkhout, Roidan Pozo,
c    Charles Romine, Henk van der Vorst,
c    Templates for the Solution of Linear Systems:
c    Building Blocks for Iterative Methods,
c    SIAM, 1994.
c    ISBN: 0898714710,
c    LC: QA297.8.T45.
c
c    Tim Kelley,
c    Iterative Methods for Linear and Nonlinear Equations,
c    SIAM, 2004,
c    ISBN: 0898713528,
c    LC: QA297.8.K45.
c
c    Yousef Saad,
c    Iterative Methods for Sparse Linear Systems,
c    Second Edition,
c    SIAM, 2003,
c    ISBN: 0898715342,
c    LC: QA188.S17.
c
c  Parameters:
c
c    Input, integer N, the order of the system.
c
c    Input, integer NZ_NUM, the number of nonzeros.
c
c    Input, integer IA(NZ_NUM), JA(NZ_NUM), the row and column
c    indices of the matrix values.
c
c    Input, double precision A(NZ_NUM), the matrix values.
c
c    Input, double precision X(N), the vector to be multiplied by A.
c
c    Output, double precision W(N), the value of A*X.
c
      implicit none

      integer n
      integer nz_num

      double precision a(nz_num)
      integer i
      integer ia(nz_num)
      integer j
      integer ja(nz_num)
      integer k
      double precision w(n)
      double precision x(n)

      do i = 1, n
        w(i) = 0.0D+00
      end do

      do k = 1, nz_num
        i = ia(k)
        j = ja(k)
        w(i) = w(i) + a(k) * x(j)
      end do

      return
      end
      subroutine diagonal_pointer_cr ( n, nz_num, ia, ja, ua )

c*********************************************************************72
c
cc DIAGONAL_POINTER_CR finds diagonal entries.
c
c  Discussion:
c
c    The matrix A is assumed to be stored in compressed row format.  Only
c    the nonzero entries of A are stored.  The vector JA stores the
c    column index of the nonzero value.  The nonzero values are sorted
c    by row, and the compressed row vector IA then has the property that
c    the entries in A and JA that correspond to row I occur in indices
c    IA[I] through IA[I+1]-1.
c
c    The array UA can be used to locate the diagonal elements of the matrix.
c
c    It is assumed that every row of the matrix includes a diagonal element,
c    and that the elements of each row have been ascending sorted.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    30 August 2010
c
c  Author:
c
c    Original C version by Lili Ju.
c    FORTRAN77 version by John Burkardt.
c
c  Parameters:
c
c    Input, integer N, the order of the system.
c
c    Input, integer NZ_NUM, the number of nonzeros.
c
c    Input, integer IA(N+1), JA(NZ_NUM), the row and column 
c    indices of the matrix values.  The row vector has been compressed.  
c    On output, the order of the entries of JA may have changed because of
c    the sorting.
c
c    Output, integer UA(N), the index of the diagonal element 
c    of each row.
c
      implicit none

      integer n
      integer nz_num

      integer i
      integer ia(n+1)
      integer k
      integer ja(nz_num)
      integer ua(n)

      do i = 1, n
        ua(i) = - 1
      end do

      do i = 1, n
        do k = ia(i), ia(i+1) - 1
          if ( ja(k) .eq. i )  then
            ua(i) = k
          end if
        end do
      end do

      return
      end
      subroutine ilu_cr ( n, nz_num, ia, ja, a, ua, l )

c*********************************************************************72
c
cc ILU_CR computes the incomplete LU factorization of a matrix.
c
c  Discussion:
c
c    The matrix A is assumed to be stored in compressed row format.  Only
c    the nonzero entries of A are stored.  The vector JA stores the
c    column index of the nonzero value.  The nonzero values are sorted
c    by row, and the compressed row vector IA then has the property that
c    the entries in A and JA that correspond to row I occur in indices
c    IA(I) through IA(I+1)-1.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    30 August 2010
c
c  Author:
c
c    Original C version by Lili Ju.
c    FORTRAN77 version by John Burkardt.
c
c  Parameters:
c
c    Input, integer N, the order of the system.
c
c    Input, integer NZ_NUM, the number of nonzeros.
c
c    Input, integer IA(N+1), JA(NZ_NUM), the row and column 
c    indices of the matrix values.  The row vector has been compressed.
c
c    Input, double precision A(NZ_NUM), the matrix values.
c
c    Input, integer UA(N), the index of the diagonal element 
c    of each row.
c
c    Output, double precision L(NZ_NUM), the ILU factorization of A.
c
      implicit none

      integer n
      integer nz_num

      double precision a(nz_num)
      integer i
      integer ia(n+1)
      integer iw(n)
      integer j
      integer ja(nz_num)
      integer jj
      integer jrow
      integer jw
      integer k
      double precision l(nz_num)
      double precision tl
      integer ua(n)
c
c  Copy A.
c
      do i = 1, nz_num
        l(i) = a(i)
      end do

      do i = 1, n
c
c  IW points to the nonzero entries in row I.
c
        do j = 1, n
          iw(j) = -1
        end do

        do k = ia(i), ia(i+1) - 1
          iw(ja(k)) = k
        end do

        do j = ia(i), ia(i+1) - 1

          jrow = ja(j)

          if ( i .le. jrow ) then
            go to 10
          end if

          tl = l(j) * l(ua(jrow))
          l(j) = tl

          do jj = ua(jrow) + 1, ia(jrow+1) - 1
            jw = iw(ja(jj))
            if ( jw /= -1 ) then 
              l(jw) = l(jw) - tl * l(jj)
            end if
          end do

        end do

10      continue

        ua(i) = j

        if ( jrow .ne. i ) then
          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) 'ILU_CR - Fatal error!'
          write ( *, '(a)' ) '  JROW ~= I'
          write ( *, '(a,i8)' ) '  JROW = ', jrow
          write ( *, '(a,i8)' ) '  I    = ', i
          stop
        end if

        if ( l(j) .eq. 0.0D+00 ) then
          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) 'ILU_CR - Fatal error!'
          write ( *, '(a,i8)' ) '  Zero pivot on step I = ', i
          write ( *, '(a,i8,a)' ) '  L(', j, ') = 0.0'
          stop
        end if

        l(j) = 1.0D+00 / l(j)

      end do

      do i = 1, n 
        l(ua(i)) = 1.0D+00 / l(ua(i))
      end do

      return
      end
      subroutine lus_cr ( n, nz_num, ia, ja, l, ua, r, z )

c*********************************************************************72
c
cc LUS_CR applies the incomplete LU preconditioner.
c
c  Discussion:
c
c    The linear system M * Z = R is solved for Z.  M is the incomplete
c    LU preconditioner matrix, and R is a vector supplied by the user.
c    So essentially, we're solving L * U * Z = R.
c
c    The matrix A is assumed to be stored in compressed row format.  Only
c    the nonzero entries of A are stored.  The vector JA stores the
c    column index of the nonzero value.  The nonzero values are sorted
c    by row, and the compressed row vector IA then has the property that
c    the entries in A and JA that correspond to row I occur in indices
c    IA(I) through IA(I+1)-1.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    30 August 2010
c
c  Author:
c
c    Original C version by Lili Ju.
c    FORTRAN77 version by John Burkardt.
c
c  Parameters:
c
c    Input, integer N, the order of the system.
c
c    Input, integer NZ_NUM, the number of nonzeros.
c
c    Input, integer IA(N+1), JA(NZ_NUM), the row and column 
c    indices of the matrix values.  The row vector has been compressed.
c
c    Input, double precision L(NZ_NUM), the matrix values.
c
c    Input, integer UA(N), the index of the diagonal element 
c    of each row.
c
c    Input, double precision R(N), the right hand side.
c
c    Output, double precision Z(N), the solution of the system M * Z = R.
c
      implicit none

      integer n
      integer nz_num

      integer i
      integer ia(n+1)
      integer j
      integer ja(nz_num)
      double precision l(nz_num)
      double precision r(n)
      integer ua(n)
      double precision w(n)
      double precision z(n)
c
c  Copy R in.
c
      do i = 1, n
        w(i) = r(i)
      end do
c
c  Solve L * w = w where L is unit lower triangular.
c
      do i = 2, n
        do j = ia(i), ua(i) - 1
          w(i) = w(i) - l(j) * w(ja(j))
        end do
      end do
c
c  Solve U * w = w, where U is upper triangular.
c
      do i = n, 1, -1
        do j = ua(i) + 1, ia(i+1) - 1
          w(i) = w(i) - l(j) * w(ja(j))
        end do
        w(i) = w(i) / l(ua(i))
      end do
c
c  Copy Z out.
c
      do i = 1, n
        z(i) = w(i)
      end do

      return
      end
      subroutine mgmres_st ( n, nz_num, ia, ja, a, x, rhs, itr_max, mr, 
     &  tol_abs, tol_rel )

c*********************************************************************72
c
cc MGMRES_ST applies restarted GMRES to a sparse triplet matrix.
c
c  Discussion:
c
c    The linear system A*X=B is solved iteratively.
c
c    The matrix A is assumed to be stored in sparse triplet form.  Only
c    the nonzero entries of A are stored.  For instance, the K-th nonzero
c    entry in the matrix is stored by:
c
c      A(K) = value of entry,
c      IA(K) = row of entry,
c      JA(K) = column of entry.
c
c    Thanks to Jesus Pueblas Sanchez-Guerra for supplying two
c    corrections to the code on 31 May 2007.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    30 August 2010
c
c  Author:
c
c    Original C version by Lili Ju.
c    FORTRAN77 version by John Burkardt.
c
c  Reference:
c
c    Richard Barrett, Michael Berry, Tony Chan, James Demmel,
c    June Donato, Jack Dongarra, Victor Eijkhout, Roidan Pozo,
c    Charles Romine, Henk van der Vorst,
c    Templates for the Solution of Linear Systems:
c    Building Blocks for Iterative Methods,
c    SIAM, 1994.
c    ISBN: 0898714710,
c    LC: QA297.8.T45.
c
c    Tim Kelley,
c    Iterative Methods for Linear and Nonlinear Equations,
c    SIAM, 2004,
c    ISBN: 0898713528,
c    LC: QA297.8.K45.
c
c    Yousef Saad,
c    Iterative Methods for Sparse Linear Systems,
c    Second Edition,
c    SIAM, 2003,
c    ISBN: 0898715342,
c    LC: QA188.S17.
c
c  Parameters:
c
c    Input, integer N, the order of the linear system.
c
c    Input, integer NZ_NUM, the number of nonzero matrix values.
c
c    Input, integer IA(NZ_NUM), JA(NZ_NUM), the row and column
c    indices of the matrix values.
c
c    Input, double precision A(NZ_NUM), the matrix values.
c
c    Input/output, double precision X(N); on input, an approximation to
c    the solution.  On output, an improved approximation.
c
c    Input, double precision RHS(N), the right hand side of the linear system.
c
c    Input, integer ITR_MAX, the maximum number of (outer) 
c    iterations to take.
c
c    Input, integer MR, the maximum number of (inner) iterations 
c    to take.  0 < MR <= N.
c
c    Input, double precision TOL_ABS, an absolute tolerance applied to the
c    current residual.
c
c    Input, double precision TOL_REL, a relative tolerance comparing the
c    current residual to the initial residual.
c
      implicit none

      integer mr
      integer n
      integer nz_num

      double precision a(nz_num)
      double precision av
      double precision c(1:mr)
      double precision delta
      parameter ( delta = 1.0D-03 )
      double precision g(1:mr+1)
      double precision h(1:mr+1,1:mr)
      double precision htmp
      integer i
      integer ia(nz_num)
      integer itr
      integer itr_max
      integer itr_used
      integer j
      integer ja(nz_num)
      integer k
      integer k_copy
      double precision mu
      double precision r(1:n)
      double precision r8vec_dot_product
      double precision r8vec_norm
      double precision rho
      double precision rho_tol
      double precision rhs(1:n)
      double precision s(1:mr)
      double precision tol_abs
      double precision tol_rel
      double precision v(1:n,1:mr+1)
      logical verbose
      parameter ( verbose = .true. )
      double precision x(1:n)
      double precision y(1:mr+1)

      itr_used = 0

      if ( n .lt. mr ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'MGMRES_ST - Fatal error!'
        write ( *, '(a)' ) '  N < MR.'
        write ( *, '(a,i8)' ) '  N = ', n
        write ( *, '(a,i8)' ) '  MR = ', mr
        stop
      end if

      do itr = 1, itr_max

        call ax_st ( n, nz_num, ia, ja, a, x, r )

        do i = 1, n
          r(i) = rhs(i) - r(i)
        end do

        rho = r8vec_norm ( n, r )

        if ( verbose ) then
          write ( *, '(a,i8,a,g14.6)' ) '  ITR = ', itr, 
     &      '  Residual = ', rho 
        end if

        if ( itr .eq. 1 ) then
          rho_tol = rho * tol_rel
        end if

        do i = 1, n
          v(i,1) = r(i) / rho
        end do

        g(1) = rho
        do i = 2, mr + 1
          g(i) = 0.0D+00
        end do

        do j = 1, mr
          do i = 1, mr + 1
            h(i,j) = 0.0D+00
          end do
        end do

        do k = 1, mr

          k_copy = k

          call ax_st ( n, nz_num, ia, ja, a, v(1,k), v(1,k+1) )

          av = r8vec_norm ( n, v(1,k+1) )

          do j = 1, k
            h(j,k) = r8vec_dot_product ( n, v(1,k+1), v(1,j) )
            do i = 1, n
              v(i,k+1) = v(i,k+1) - h(j,k) * v(i,j)
            end do
          end do

          h(k+1,k) = r8vec_norm ( n, v(1,k+1) )

          if ( av + delta * h(k+1,k) .eq. av ) then

            do j = 1, k
              htmp = r8vec_dot_product ( n, v(1,k+1), v(1,j) )
              h(j,k) = h(j,k) + htmp
              do i = 1, n
                v(i,k+1) = v(i,k+1) - htmp * v(i,j)
              end do
            end do

            h(k+1,k) = r8vec_norm ( n, v(1,k+1) )

          end if

          if ( h(k+1,k) .ne. 0.0D+00 ) then
            do i = 1, n
              v(i,k+1) = v(i,k+1) / h(k+1,k)
            end do
          end if

          if ( 1 .lt. k ) then

            do i = 1, k + 1
              y(i) = h(i,k)
            end do

            do j = 1, k-1
              call mult_givens ( c(j), s(j), j, y )
            end do

            do i = 1, k + 1
              h(i,k) = y(i)
            end do

          end if

          mu = sqrt ( h(k,k)**2 + h(k+1,k)**2 )
          c(k) = h(k,k) / mu
          s(k) = - h(k+1,k) / mu
          h(k,k) = c(k) * h(k,k) - s(k) * h(k+1,k)
          h(k+1,k) = 0.0D+00
          call mult_givens ( c(k), s(k), k, g )
          rho = abs ( g(k+1) )

          itr_used = itr_used + 1

          if ( verbose ) then
            write ( *, '(a,i8,a,g14.6)' ) 
     &        '  K =   ', k, '  Residual = ', rho 
          end if

          if ( rho .le. rho_tol .and. rho .le. tol_abs ) then
            go to 10
          end if

        end do

10      continue

        k = k_copy - 1

        y(k+1) = g(k+1) / h(k+1,k+1)

        do i = k, 1, -1
          y(i) = g(i)
          do j = i + 1, k + 1
            y(i) = y(i) - h(i,j) * y(j)
          end do
          y(i) = y(i) / h(i,i)
        end do

        do i = 1, n
          do j = 1, k + 1
            x(i) = x(i) + v(i,j) * y(j)
          end do
        end do
             
        if ( rho .le. rho_tol .and. rho .le. tol_abs ) then
          go to 20
        end if

      end do

20    continue

      if ( verbose ) then
        write ( *, '(a)'       ) ' '
        write ( *, '(a)'       ) 'MGMRES_ST:'
        write ( *, '(a,i8)'    ) '  Iterations = ', itr_used
        write ( *, '(a,g14.6)' ) '  Final residual = ', rho
      end if

      return
      end
      subroutine mult_givens ( c, s, k, g )

c*********************************************************************72
c
cc MULT_GIVENS applies a Givens rotation to two successive entries of a vector.
c
c  Discussion:
c
c    In order to make it easier to compare this code with the Original C,
c    the vector indexing is 0-based.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    30 August 2010
c
c  Author:
c
c    Original C version by Lili Ju.
c    FORTRAN77 version by John Burkardt.
c
c  Reference:
c
c    Richard Barrett, Michael Berry, Tony Chan, James Demmel,
c    June Donato, Jack Dongarra, Victor Eijkhout, Roidan Pozo,
c    Charles Romine, Henk van der Vorst,
c    Templates for the Solution of Linear Systems:
c    Building Blocks for Iterative Methods,
c    SIAM, 1994.
c    ISBN: 0898714710,
c    LC: QA297.8.T45.
c
c    Tim Kelley,
c    Iterative Methods for Linear and Nonlinear Equations,
c    SIAM, 2004,
c    ISBN: 0898713528,
c    LC: QA297.8.K45.
c
c    Yousef Saad,
c    Iterative Methods for Sparse Linear Systems,
c    Second Edition,
c    SIAM, 2003,
c    ISBN: 0898715342,
c    LC: QA188.S17.
c
c  Parameters:
c
c    Input, double precision C, S, the cosine and sine of a Givens
c    rotation.
c
c    Input, integer K, indicates the location of the first 
c    vector entry.
c
c    Input/output, double precision G(1:K+1), the vector to be modified. 
c    On output, the Givens rotation has been applied to entries G(K) and G(K+1).
c
      implicit none

      integer k

      double precision c
      double precision g(k+1)
      double precision g1
      double precision g2
      double precision s

      g1 = c * g(k) - s * g(k+1)
      g2 = s * g(k) + c * g(k+1)

      g(k)   = g1
      g(k+1) = g2

      return
      end
      subroutine pmgmres_ilu_cr ( n, nz_num, ia, ja, a, x, rhs, 
     &  itr_max, mr, tol_abs, tol_rel )

c*********************************************************************72
c
cc PMGMRES_ILU_CR applies the preconditioned restarted GMRES algorithm.
c
c  Discussion:
c
c    The matrix A is assumed to be stored in compressed row format.  Only
c    the nonzero entries of A are stored.  The vector JA stores the
c    column index of the nonzero value.  The nonzero values are sorted
c    by row, and the compressed row vector IA then has the property that
c    the entries in A and JA that correspond to row I occur in indices
c    IA(I) through IA(I+1)-1.
c
c    This routine uses the incomplete LU decomposition for the
c    preconditioning.  This preconditioner requires that the sparse
c    matrix data structure supplies a storage position for each diagonal
c    element of the matrix A, and that each diagonal element of the
c    matrix A is not zero.
c
c    Thanks to Jesus Pueblas Sanchez-Guerra for supplying two
c    corrections to the code on 31 May 2007.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    29 August 2012
c
c  Author:
c
c    Original C version by Lili Ju.
c    FORTRAN77 version by John Burkardt.
c
c  Reference:
c
c    Richard Barrett, Michael Berry, Tony Chan, James Demmel,
c    June Donato, Jack Dongarra, Victor Eijkhout, Roidan Pozo,
c    Charles Romine, Henk van der Vorst,
c    Templates for the Solution of Linear Systems:
c    Building Blocks for Iterative Methods,
c    SIAM, 1994.
c    ISBN: 0898714710,
c    LC: QA297.8.T45.
c
c    Tim Kelley,
c    Iterative Methods for Linear and Nonlinear Equations,
c    SIAM, 2004,
c    ISBN: 0898713528,
c    LC: QA297.8.K45.
c
c    Yousef Saad,
c    Iterative Methods for Sparse Linear Systems,
c    Second Edition,
c    SIAM, 2003,
c    ISBN: 0898715342,
c    LC: QA188.S17.
c
c  Parameters:
c
c    Input, integer N, the order of the linear system.
c
c    Input, integer NZ_NUM, the number of nonzero matrix values.
c
c    Input, integer IA(N+1), JA(NZ_NUM), the row and column indices
c    of the matrix values.  The row vector has been compressed.
c
c    Input, double precision A(NZ_NUM), the matrix values.
c
c    Input/output, double precision X(N); on input, an approximation to
c    the solution.  On output, an improved approximation.
c
c    Input, double precision RHS(N), the right hand side of the linear system.
c
c    Input, integer ITR_MAX, the maximum number of (outer) 
c    iterations to take.
c
c    Input, integer MR, the maximum number of (inner) iterations 
c    to take.  MR must be less than N.
c
c    Input, double precision TOL_ABS, an absolute tolerance applied to the
c    current residual.
c
c    Input, double precision TOL_REL, a relative tolerance comparing the
c    current residual to the initial residual.
c
      implicit none

      integer mr
      integer n
      integer nz_num

      double precision a(nz_num)
      double precision av
      double precision c(mr+1)
      double precision delta
      parameter ( delta = 1.0D-03 )
      double precision g(mr+1)
      double precision h(mr+1,mr)
      double precision htmp
      integer i
      integer ia(n+1)
      integer itr
      integer itr_max
      integer itr_used
      integer j
      integer ja(nz_num)
      integer k
      integer k_copy
      double precision l(ia(n+1)+1)
      double precision mu
      double precision r(n)
      double precision r8vec_dot_product
      double precision rho
      double precision rho_tol
      double precision rhs(n)
      double precision s(mr+1)
      double precision tol_abs
      double precision tol_rel
      integer ua(n)
      double precision v(n,mr+1)
      logical verbose
      parameter ( verbose = .true. )
      double precision x(n)
      double precision y(mr+1)

      itr_used = 0

      call rearrange_cr ( n, nz_num, ia, ja, a )

      call diagonal_pointer_cr ( n, nz_num, ia, ja, ua )

      call ilu_cr ( n, nz_num, ia, ja, a, ua, l )

      if ( verbose ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'PMGMRES_ILU_CR'
        write ( *, '(a,i4)' ) '  Number of unknowns = ', n
      end if

      do itr = 1, itr_max

        call ax_cr ( n, nz_num, ia, ja, a, x, r )

        do j = 1, n
          r(j) = rhs(j) - r(j)
        end do

        call lus_cr ( n, nz_num, ia, ja, l, ua, r, r )

        rho = sqrt ( r8vec_dot_product ( n, r, r ) )

        if ( verbose ) then
          write ( *, '(a,i4,a,g14.6)' ) 
     &      '  ITR = ', itr, '  Residual = ', rho
        end if

        if ( itr .eq. 1 ) then
          rho_tol = rho * tol_rel
        end if

        do i = 1, n
          v(i,1) = r(i) / rho
        end do

        g(1) = rho
        do j = 2, mr + 1
          g(j) = 0.0D+00
        end do

        do j = 1, mr
          do i = 1, mr + 1
            h(i,j) = 0.0D+00
          end do
        end do

        do k = 1, mr

          k_copy = k

          call ax_cr ( n, nz_num, ia, ja, a, v(1,k), v(1,k+1) ) 

          call lus_cr ( n, nz_num, ia, ja, l, ua, v(1,k+1), v(1,k+1) )

          av = sqrt ( r8vec_dot_product ( n, v(1,k+1), v(1,k+1) ) )

          do j = 1, k
            h(j,k) = r8vec_dot_product ( n, v(1,k+1), v(1,j) )
            do i = 1, n
              v(i,k+1) = v(i,k+1) - v(i,j) * h(j,k)
            end do
          end do

          h(k+1,k) = sqrt ( 
     &      r8vec_dot_product ( n, v(1,k+1), v(1,k+1) ) )

          if ( ( av + delta * h(k+1,k)) .eq. av ) then
            do j = 1, k
              htmp = r8vec_dot_product ( n, v(1,k+1), v(1,j) )
              h(j,k) = h(j,k) + htmp
              do i = 1, n
                v(i,k+1) = v(i,k+1) - v(i,j) * htmp
              end do
            end do
            h(k+1,k) = sqrt ( 
     &        r8vec_dot_product ( n, v(1,k+1), v(1,k+1) ) )
          end if

          if ( h(k+1,k) .ne. 0.0D+00 ) then
            do i = 1, n
              v(i,k+1) = v(i,k+1) / h(k+1,k)
            end do
          end if

          if ( 1 .lt. k ) then
            do i = 1, k + 1
              y(i) = h(i,k)
            end do
            do j = 1, k - 1
              call mult_givens ( c(j), s(j), j, y )
            end do
            do i = 1, k + 1
              h(i,k) = y(i)
            end do
          end if

          mu = sqrt ( h(k,k)**2 + h(k+1,k)**2 )

          c(k) = h(k,k) / mu
          s(k) = -h(k+1,k) / mu
          h(k,k) = c(k) * h(k,k) - s(k) * h(k+1,k)
          h(k+1,k) = 0.0D+00
          call mult_givens ( c(k), s(k), k, g )

          rho = abs ( g(k+1) )

          itr_used = itr_used + 1

          if ( verbose ) then
            write ( *, '(a,i4,a,g14.6)' ) 
     &        '  K = ', k, '  Residual = ', rho
          end if

          if ( rho .le. rho_tol .and. rho .le. tol_abs ) then
            go to 10
          end if

        end do

10      continue

        k = k_copy - 1

        y(k+1) = g(k+1) / h(k+1,k+1)

        do i = k, 1, -1
          y(i) = g(i)
          do j = i + 1, k + 1
            y(i) = y(i) - h(i,j) * y(j)
          end do
          y(i) = y(i) / h(i,i)
        end do

        do i = 1, n
          do j = 1, k + 1
            x(i) = x(i) + v(i,j) * y(j)
          end do
        end do

        if ( rho .le. rho_tol .and. rho .le. tol_abs ) then
          go to 20
        end if

      end do

20    continue

      if ( verbose ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'PMGMRES_ILU_CR:'
        write ( *, '(a,i6)' ) '  Iterations = ', itr_used
        write ( *, '(a,g14.6)' ) '  Final residual = ', rho
      end if

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
      parameter ( i4_huge = 2147483647 )
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
          seed = seed + i4_huge
        end if

        r(i) = dble ( seed ) * 4.656612875D-10

      end do

      return
      end
      subroutine rearrange_cr ( n, nz_num, ia, ja, a )

c*********************************************************************72
c
cc REARRANGE_CR sorts a sparse compressed row matrix.
c
c  Discussion:
c
c    This routine guarantees that the entries in the CR matrix
c    are properly sorted.
c
c    After the sorting, the entries of the matrix are rearranged in such
c    a way that the entries of each column are listed in ascending order
c    of their column values.
c
c    The matrix A is assumed to be stored in compressed row format.  Only
c    the nonzero entries of A are stored.  The vector JA stores the
c    column index of the nonzero value.  The nonzero values are sorted
c    by row, and the compressed row vector IA then has the property that
c    the entries in A and JA that correspond to row I occur in indices
c    IA(I) through IA(I+1)-1.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    30 August 2010
c
c  Author:
c
c    Original C version by Lili Ju.
c    FORTRAN77 version by John Burkardt.
c
c  Reference:
c
c    Richard Barrett, Michael Berry, Tony Chan, James Demmel,
c    June Donato, Jack Dongarra, Victor Eijkhout, Roidan Pozo,
c    Charles Romine, Henk van der Vorst,
c    Templates for the Solution of Linear Systems:
c    Building Blocks for Iterative Methods,
c    SIAM, 1994.
c    ISBN: 0898714710,
c    LC: QA297.8.T45.
c
c    Tim Kelley,
c    Iterative Methods for Linear and Nonlinear Equations,
c    SIAM, 2004,
c    ISBN: 0898713528,
c    LC: QA297.8.K45.
c
c    Yousef Saad,
c    Iterative Methods for Sparse Linear Systems,
c    Second Edition,
c    SIAM, 2003,
c    ISBN: 0898715342,
c    LC: QA188.S17.
c
c  Parameters:
c
c    Input, integer N, the order of the system.
c
c    Input, integer NZ_NUM, the number of nonzeros.
c
c    Input, integer IA(N+1), the compressed row indices..
c
c    Input/output, integer JA(NZ_NUM), the column indices.  
c    On output, these may have been rearranged by the sorting.
c
c    Input/output, double precision A(NZ_NUM), the matrix values.  On output,
c    the matrix values may have been moved somewhat because of the sorting.
c
      implicit none

      integer n
      integer nz_num

      double precision a(nz_num)
      integer i
      integer ia(n+1)
      integer i4temp
      integer ja(nz_num)
      integer k
      integer l
      double precision r8temp

      do i = 1, n

        do k = ia(i), ia(i+1) - 2
          do l = k + 1, ia(i+1) - 1

            if ( ja(l) .lt. ja(k) ) then
              i4temp = ja(l)
              ja(l)  = ja(k)
              ja(k)  = i4temp

              r8temp = a(l)
              a(l)   = a(k)
              a(k)   = r8temp
            end if

          end do
        end do

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
