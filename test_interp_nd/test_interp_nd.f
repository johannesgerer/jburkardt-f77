      function csevl ( x, a, n )

c*********************************************************************72
c
cc CSEVL evaluates a Chebyshev series.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    09 March 2010
c
c  Author:
c
c    This FORTRAN77 version by John Burkardt.
c
c  Reference:
c
c    Roger Broucke,
c    Algorithm 446:
c    Ten Subroutines for the Manipulation of Chebyshev Series,
c    Communications of the ACM,
c    Volume 16, Number 4, April 1973, pages 254-256.
c
c  Parameters:
c
c    Input, double precision X, the evaluation point.
c
c    Input, double precision CS(N), the Chebyshev coefficients.
c
c    Input, integer N, the number of Chebyshev coefficients.
c
c    Output, double precision CSEVL, the Chebyshev series evaluated at X.
c
      implicit none

      integer n

      double precision a(n)
      double precision b0
      double precision b1
      double precision b2
      double precision csevl
      integer i
      double precision twox
      double precision x

      if ( n .lt. 1 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'CSEVL - Fatal error!'
        write ( *, '(a)' ) '  Number of terms <= 0.'
        stop
      end if

      if ( 1000 .lt. n ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'CSEVL - Fatal error!'
        write ( *, '(a)' ) '  Number of terms > 1000.'
        stop
      end if

      if ( x .lt. -1.1D+00 .or. 1.1D+00 .lt. x ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'CSEVL - Fatal error!'
        write ( *, '(a)' ) '  X outside (-1,+1)'
        write ( *, '(a,g14.6)' ) '  X = ', x
        stop
      end if

      twox = 2.0D+00 * x
      b1 = 0.0D+00
      b0 = 0.0D+00

      do i = n, 1, -1
        b2 = b1
        b1 = b0
        b0 = twox * b1 - b2 + a(i)
      end do

      csevl = 0.5D+00 * ( b0 - b2 )

      return
      end
      function i4vec_sum ( n, a )

c*********************************************************************72
c
cc I4VEC_SUM returns the sum of the entries of an I4VEC.
c
c  Discussion:
c
c    An I4VEC is a vector of I4's.
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
      function inits ( dos, nos, eta )

c*********************************************************************72
c
cc INITS initializes a Chebyshev series.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    09 March 2010
c
c  Author:
c
c    This FORTRAN77 version by John Burkardt.
c
c  Reference:
c
c    Roger Broucke,
c    Algorithm 446:
c    Ten Subroutines for the Manipulation of Chebyshev Series,
c    Communications of the ACM,
c    Volume 16, Number 4, April 1973, pages 254-256.
c
c  Parameters:
c
c    Input, double precision DOS(NOS), the Chebyshev coefficients.
c
c    Input, integer NOS, the number of coefficients.
c
c    Input, double precision ETA, the desired accuracy.
c
c    Output, integer INITS, the number of terms of the series needed
c    to ensure the requested accuracy.
c
      implicit none

      integer nos

      double precision dos(nos)
      double precision err 
      double precision eta
      integer i
      integer inits

      if ( nos .lt. 1 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'INITS - Fatal error!'
        write ( *, '(a)' ) '  Number of coefficients < 1.'
        stop
      end if

      err = 0.0D+00

      do i = nos, 1, -1
        err = err + dabs ( dos(i) )
        if ( eta .lt. err ) then
          inits = i
          return
        end if
      end do

      inits = nos
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'INITS - Warning!'
      write ( *, '(a)' ) '  ETA may be too small.'

      return
      end
      subroutine p00_c ( prob, m, seed, c )

c*********************************************************************72
c
cc P00_CW computes a random C parameter vector for any problem.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    28 August 2012
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer PROB, the problem number.
c
c    Input, integer M, the spatial dimension.
c
c    Input/output, integer SEED, a seed for the random 
c    number generator.
c
c    Output, double precision C(M), the parameter vector.
c
      implicit none

      integer m

      double precision b(6)
      double precision c(m)
      double precision c_sum
      integer i
      integer prob
      double precision r8vec_sum
      integer seed

      save b

      data b / 1.5D+00, 0.0D+00, 1.85D+00, 7.03D+00, 20.4D+00, 4.3D+00 /

      b(2) = dble ( m )

      call r8vec_uniform_01 ( m, seed, c )
      c_sum = r8vec_sum ( m, c )
      c(1:m) = b(prob) * c(1:m) / c_sum

      return
      end
      subroutine p00_d ( prob, m, id, c, w, n, x, d )

c*********************************************************************72
c
cc P00_D returns a derivative component of any function.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    28 August 2012
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer PROB, the index of the function.
c
c    Input, integer M, the spatial dimension.
c
c    Input, integer ID, the spatial coordinate to differentiate.
c
c    Input, double precision C(M), W(M), the problem parameters.
c
c    Input, integer N, the number of evaluation points.
c
c    Input, double precision X(M,N), the evalution points.
c
c    Output, double precision D(N), the ID-th derivative component.
c
      implicit none

      integer m
      integer n

      double precision c(m)
      double precision d(n)
      integer id
      integer prob
      double precision w(m)
      double precision x(m,n)  

      if ( id .lt. 0 .or. m .lt. id ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'P00_D - Fatal error!'
        write ( *, '(a,i4)' ) '  Illegal spatial coordinate ID = ', id
        stop
      end if

      if ( prob .eq. 1 ) then
        call p01_d ( m, id, c, w, n, x, d )
      else if ( prob .eq. 2 ) then
        call p02_d ( m, id, c, w, n, x, d )
      else if ( prob .eq. 3 ) then
        call p03_d ( m, id, c, w, n, x, d )
      else if ( prob .eq. 4 ) then
        call p04_d ( m, id, c, w, n, x, d )
      else if ( prob .eq. 5 ) then
        call p05_d ( m, id, c, w, n, x, d )
      else if ( prob .eq. 6 ) then
        call p06_d ( m, id, c, w, n, x, d )
      else
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'P00_D - Fatal error!'
        write ( *, '(a,i8)' ) '  Illegal function index PROB = ', prob
        stop
      end if

      return
      end
      subroutine p00_f ( prob, m, c, w, n, x, f )

c*********************************************************************72
c
cc P00_F returns the value of any function.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    28 August 2012
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer PROB, the index of the function.
c
c    Input, integer M, the spatial dimension.
c
c    Input, double precision C(M), W(M), the problem parameters.
c
c    Input, integer N, the number of evaluation points.
c
c    Input, double precision X(M,N), the evalution points.
c
c    Output, double precision F(N), the function values.
c
      implicit none

      integer m
      integer n

      double precision c(m)
      double precision f(n)
      integer prob
      double precision w(m)
      double precision x(m,n)  

      if ( prob .eq. 1 ) then
        call p01_f ( m, c, w, n, x, f )
      else if ( prob .eq. 2 ) then
        call p02_f ( m, c, w, n, x, f )
      else if ( prob .eq. 3 ) then
        call p03_f ( m, c, w, n, x, f )
      else if ( prob .eq. 4 ) then
        call p04_f ( m, c, w, n, x, f )
      else if ( prob .eq. 5 ) then
        call p05_f ( m, c, w, n, x, f )
      else if ( prob .eq. 6 ) then
        call p06_f ( m, c, w, n, x, f )
      else
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'P00_F - Fatal error!'
        write ( *, '(a,i8)' ) '  Illegal function index PROB = ', prob
        stop
      end if

      return
      end
      subroutine p00_prob_num ( prob_num )

c*********************************************************************72
c
cc P00_PROB_NUM returns the number of test functions available.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    28 August 2012
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c   Output, integer PROB_NUM, the number of test functions.
c
      implicit none

      integer prob_num

      prob_num = 6

      return
      end
      subroutine p00_q ( prob, m, c, w, q )

c*********************************************************************72
c
cc P00_Q returns the integral of any function.
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
c    John Burkardt
c
c  Parameters:
c
c    Input, integer PROB, the index of the function.
c
c    Input, integer M, the spatial dimension.
c
c    Input, double precision C(M), W(M), the problem parameters.
c
c    Output, double precision Q, the integral.
c
      implicit none

      integer m

      double precision c(m)
      integer prob
      double precision q
      double precision w(m)

      if ( prob .eq. 1 ) then
        call p01_q ( m, c, w, q )
      else if ( prob .eq. 2 ) then
        call p02_q ( m, c, w, q )
      else if ( prob .eq. 3 ) then
        call p03_q ( m, c, w, q )
      else if ( prob .eq. 4 ) then
        call p04_q ( m, c, w, q )
      else if ( prob .eq. 5 ) then
        call p05_q ( m, c, w, q )
      else if ( prob .eq. 6 ) then
        call p06_q ( m, c, w, q )
      else
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'P00_Q - Fatal error!'
        write ( *, '(a,i8)' ) '  Illegal function index PROB = ', prob
        stop
      end if

      return
      end
      subroutine p00_title ( prob, title )

c*********************************************************************72
c
cc P00_TITLE returns the title for any function.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    28 August 2012
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer PROB, the index of the function.
c
c    Output, character * ( * ) TITLE, the function title.
c
      implicit none

      integer prob
      character * ( * ) title

      if ( prob .eq. 1 ) then
        call p01_title ( title )
      else if ( prob .eq. 2 ) then
        call p02_title ( title )
      else if ( prob .eq. 3 ) then
        call p03_title ( title )
      else if ( prob .eq. 4 ) then
        call p04_title ( title )
      else if ( prob .eq. 5 ) then
        call p05_title ( title )
      else if ( prob .eq. 6 ) then
        call p06_title ( title )
      else
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'P00_TITLE - Fatal error!'
        write ( *, '(a,i8)' ) '  Illegal function index PROB = ', prob
        stop
      end if

      return
      end
      subroutine p00_w ( prob, m, seed, w )

c*********************************************************************72
c
cc P00_W computes a random W parameter vector for any problem.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    28 August 2012
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer PROB, the problem number.
c
c    Input, integer M, the spatial dimension.
c
c    Input/output, integer SEED, a seed for the random 
c    number generator.
c
c    Output, double precision W(M), the parameter vector.
c
      implicit none

      integer m

      integer prob
      integer seed
      double precision w(m)

      call r8vec_uniform_01 ( m, seed, w )

      return
      end
      subroutine p01_d ( m, id, c, w, n, x, d )

c*********************************************************************72
c
cc P01_D evaluates any derivative component for problem p01.
c
c  Discussion:
c
c    f(x) = cos ( 2 * pi * w(1) + sum ( c(1:m) * x(1:m) ) )
c
c    Default values are:
c
c    c(1:m) = 1/m
c    w(1) = 0.3
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    28 August 2012
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Alan Genz,
c    A Package for Testing Multiple Integration Subroutines,
c    in Numerical Integration: Recent Developments, Software
c    and Applications,
c    edited by Patrick Keast and Graeme Fairweather,
c    Reidel, 1987, pages 337-340,
c    ISBN: 9027725144,
c    LC: QA299.3.N38.
c
c  Parameters:
c
c    Input, integer M, the dimension of the argument.
c
c    Input, integer ID, the spatial coordinate to differentiate.
c
c    Input, double precision C(M), W(M), the problem parameters.
c
c    Input, integer N, the number of points.
c
c    Input, double precision X(M,N), the evaluation points.
c
c    Output, double precision D(N), the ID-th derivative component.
c
      implicit none

      integer m
      integer n

      double precision c(m)
      double precision d(n)
      integer i
      integer id
      integer j
      double precision pi
      parameter ( pi = 3.141592653589793D+00 )
      double precision w(m)
      double precision x(m,n)

      do j = 1, n
        d(j) = 2.0D+00 * pi * w(1)
      end do

      do i = 1, m
        do j = 1, n
          d(j) = d(j) + c(i) * x(i,j)
        end do
      end do

      do j = 1, n
        d(j) = - c(id) * sin ( d(j) )
      end do

      return
      end
      subroutine p01_f ( m, c, w, n, x, f )

c*********************************************************************72
c
cc P01_F evaluates the function for problem p01.
c
c  Discussion:
c
c    f(x) = cos ( 2 * pi * w(1) + sum ( c(1:m) * x(1:m) ) )
c
c    Default values are:
c
c    c(1:m) = 1/m
c    w(1) = 0.3
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    28 August 2012
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Alan Genz,
c    A Package for Testing Multiple Integration Subroutines,
c    in Numerical Integration: Recent Developments, Software
c    and Applications,
c    edited by Patrick Keast and Graeme Fairweather,
c    Reidel, 1987, pages 337-340,
c    ISBN: 9027725144,
c    LC: QA299.3.N38.
c
c  Parameters:
c
c    Input, integer M, the dimension of the argument.
c
c    Input, double precision C(M), W(M), the problem parameters.
c
c    Input, integer N, the number of points.
c
c    Input, double precision X(M,N), the evaluation points.
c
c    Output, double precision F(N), the function values.
c
      implicit none

      integer m
      integer n

      double precision c(m)
      double precision f(n)
      integer i
      integer j
      double precision pi
      parameter ( pi = 3.141592653589793D+00 )
      double precision w(m)
      double precision x(m,n)

      do j = 1, n
        f(j) = 2.0D+00 * pi * w(1)
      end do

      do i = 1, m
        do j = 1, n
          f(j) = f(j) + c(i) * x(i,j)
        end do
      end do

      do j = 1, n
        f(j) = cos ( f(j) )
      end do

      return
      end
      subroutine p01_q ( m, c, w, q )

c*********************************************************************72
c
cc P01_Q evaluates the integral for problem p01.
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
c    John Burkardt
c
c  Reference:
c
c    Alan Genz,
c    A Package for Testing Multiple Integration Subroutines,
c    in Numerical Integration: Recent Developments, Software
c    and Applications,
c    edited by Patrick Keast and Graeme Fairweather,
c    Reidel, 1987, pages 337-340,
c    ISBN: 9027725144,
c    LC: QA299.3.N38.
c
c  Parameters:
c
c    Input, integer M, the dimension of the argument.
c
c    Input, double precision C(M), W(M), the problem parameters.
c
c    Output, double precision Q, the integral.
c
      implicit none

      integer m

      double precision c(m)
      double precision c_prod
      double precision c_sum
      integer i
      double precision pi
      parameter ( pi = 3.141592653589793D+00 )
      double precision q
      double precision r8vec_sum
      double precision w(m)

      c_sum = r8vec_sum ( m, c )
      c_prod = 1.0D+00
      do i = 1, m
        c_prod = c_prod * sin ( 0.5D+00 * c(i) ) / c(i)
      end do

      q = 2.0D+00 ** m * 
     &  cos ( 2.0D+00 * pi * w(1) + 0.5D+00 * c_sum ) * c_prod

      return
      end
      subroutine p01_title ( title )

c*********************************************************************72
c
cc P01_TITLE returns the name of problem p01.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    28 August 2012
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, character * ( * ) TITLE, the title of the problem.
c
      implicit none

      character * ( * ) title

      title = 'Oscillatory'

      return
      end
      subroutine p02_d ( m, id, c, w, n, x, d )

c*********************************************************************72
c
cc P02_D evaluates an derivative component for problem p02.
c
c  Discussion:
c
c    f(x) = 1 / product ( c(1:m)^(-2) + ( x(1:m) - w(1:m) )^2 )
c
c    Default values are:
c
c    c(1:m) = 1
c    w(1:m) = 0.5
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    28 August 2012
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Alan Genz,
c    A Package for Testing Multiple Integration Subroutines,
c    in Numerical Integration: Recent Developments, Software
c    and Applications,
c    edited by Patrick Keast and Graeme Fairweather,
c    Reidel, 1987, pages 337-340,
c    ISBN: 9027725144,
c    LC: QA299.3.N38.
c
c  Parameters:
c
c    Input, integer M, the dimension of the argument.
c
c    Input, integer ID, the spatial coordinate to differentiate.
c
c    Input, double precision C(M), W(M), the problem parameters.
c
c    Input, integer N, the number of points.
c
c    Input, double precision X(M,N), the evaluation points.
c
c    Output, double precision D(N), the ID-th derivative component.
c
      implicit none

      integer m
      integer n

      double precision c(m)
      double precision d(n)
      integer i
      integer id
      integer j
      double precision w(m)
      double precision x(m,n)

      do j = 1, n
        d(j) = 1.0D+00
      end do

      do i = 1, m
        do j = 1, n
          d(j) = d(j) * ( c(i)**(-2) + ( x(i,j) - w(i) )**2 )
        end do
      end do

      do j = 1, n
        d(j) = - 2.0D+00 / d(j) * ( x(id,j) - w(id) ) / 
     &    ( c(id) ** ( -2 ) + ( x(id,j) - w(id) ) ** 2 )
      end do

      return
      end
      subroutine p02_f ( m, c, w, n, x, f )

c*********************************************************************72
c
cc P02_F evaluates the function for problem p02.
c
c  Discussion:
c
c    f(x) = 1 / product ( c(1:m)^(-2) + ( x(1:m) - w(1:m) )^2 )
c
c    Default values are:
c
c    c(1:m) = 1
c    w(1:m) = 0.5
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    28 August 2012
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Alan Genz,
c    A Package for Testing Multiple Integration Subroutines,
c    in Numerical Integration: Recent Developments, Software
c    and Applications,
c    edited by Patrick Keast and Graeme Fairweather,
c    Reidel, 1987, pages 337-340,
c    ISBN: 9027725144,
c    LC: QA299.3.N38.
c
c  Parameters:
c
c    Input, integer M, the dimension of the argument.
c
c    Input, double precision C(M), W(M), the problem parameters.
c
c    Input, integer N, the number of points.
c
c    Input, double precision X(M,N), the evaluation points.
c
c    Output, double precision F(N), the function values.
c
      implicit none

      integer m
      integer n

      double precision c(m)
      double precision f(n)
      integer i
      integer j
      double precision w(m)
      double precision x(m,n)

      do j = 1, n
        f(j) = 1.0D+00
      end do

      do i = 1, m
        do j = 1, n
          f(j) = f(j) * ( c(i)**(-2) + ( x(i,j) - w(i) )**2 )
        end do
      end do

      do j = 1, n
        f(j) = 1.0D+00 / f(j)
      end do

      return
      end
      subroutine p02_q ( m, c, w, q )

c*********************************************************************72
c
cc P02_Q evaluates the integral for problem p02.
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
c    John Burkardt
c
c  Reference:
c
c    Alan Genz,
c    A Package for Testing Multiple Integration Subroutines,
c    in Numerical Integration: Recent Developments, Software
c    and Applications,
c    edited by Patrick Keast and Graeme Fairweather,
c    Reidel, 1987, pages 337-340,
c    ISBN: 9027725144,
c    LC: QA299.3.N38.
c
c  Parameters:
c
c    Input, integer M, the dimension of the argument.
c
c    Input, double precision C(M), W(M), the problem parameters.
c
c    Output, double precision Q, the integral.
c
      implicit none

      integer m

      double precision c(m)
      integer i
      double precision q
      double precision w(m)

      q = 1.0D+00
      do i = 1, m
        q = q * 
     &    (   atan ( ( 1.0D+00 - w(i) ) * c(i) )
     &      + atan (             w(i)   * c(i) )
     &    ) * c(i)
      end do

      return
      end
      subroutine p02_title ( title )

c*********************************************************************72
c
cc P02_TITLE returns the title of problem p02.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    28 August 2012
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, character * ( * ) TITLE, the title of the problem.
c
      implicit none

      character * ( * ) title

      title = 'Product Peak'

      return
      end
      subroutine p03_d ( m, id, c, w, n, x, d )

c*********************************************************************72
c
cc P03_D evaluates any derivative component for problem p03.
c
c  Discussion:
c
c    f(x) = 1 / ( 1 + sum ( c(1:m) * x(1:m) ) ) ^ ( m + 1 )
c
c    Default values are:
c
c    c(1:m) = 1/m
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    28 August 2012
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Alan Genz,
c    A Package for Testing Multiple Integration Subroutines,
c    in Numerical Integration: Recent Developments, Software
c    and Applications,
c    edited by Patrick Keast and Graeme Fairweather,
c    Reidel, 1987, pages 337-340,
c    ISBN: 9027725144,
c    LC: QA299.3.N38.
c
c  Parameters:
c
c    Input, integer M, the dimension of the argument.
c
c    Input, integer ID, the spatial coordinate to differentiate.
c
c    Input, double precision C(M), W(M), the problem parameters.
c
c    Input, integer N, the number of points.
c
c    Input, double precision X(M,N), the evaluation points.
c
c    Output, double precision D(N), the ID-th derivative component.
c
      implicit none

      integer m
      integer n

      double precision c(m)
      double precision d(n)
      integer i
      integer id
      integer j
      double precision w(m)
      double precision x(m,n)

      do j = 1, n
        d(j) = 1.0D+00
      end do

      do i = 1, m
        do j = 1, n
          d(j) = d(j) + c(i) * x(i,j)
        end do
      end do

      do j = 1, n
        d(j) = - c(id) * dble ( m + 1 ) / d(j) ** ( m + 2 )
      end do

      return
      end
      subroutine p03_f ( m, c, w, n, x, f )

c*********************************************************************72
c
cc P03_F evaluates the function for problem p03.
c
c  Discussion:
c
c    f(x) = 1 / ( 1 + sum ( c(1:m) * x(1:m) ) ) ^ ( m + 1 )
c
c    Default values are:
c
c    c(1:m) = 1/m
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    28 August 2012
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Alan Genz,
c    A Package for Testing Multiple Integration Subroutines,
c    in Numerical Integration: Recent Developments, Software
c    and Applications,
c    edited by Patrick Keast and Graeme Fairweather,
c    Reidel, 1987, pages 337-340,
c    ISBN: 9027725144,
c    LC: QA299.3.N38.
c
c  Parameters:
c
c    Input, integer M, the dimension of the argument.
c
c    Input, double precision C(M), W(M), the problem parameters.
c
c    Input, integer N, the number of points.
c
c    Input, double precision X(M,N), the evaluation points.
c
c    Output, double precision F(N), the function values.
c
      implicit none

      integer m
      integer n

      double precision c(m)
      double precision f(n)
      integer i
      integer j
      double precision w(m)
      double precision x(m,n)

      do j = 1, n
        f(j) = 1.0D+00
      end do

      do i = 1, m
        do j = 1, n
          f(j) = f(j) + c(i) * x(i,j)
        end do
      end do

      do j = 1, n
        f(j) = 1.0D+00 / f(j) ** ( m + 1 )
      end do

      return
      end
      subroutine p03_q ( m, c, w, q )

c*********************************************************************72
c
cc P03_Q evaluates the integral for problem p03.
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
c    John Burkardt
c
c  Reference:
c
c    Alan Genz,
c    A Package for Testing Multiple Integration Subroutines,
c    in Numerical Integration: Recent Developments, Software
c    and Applications,
c    edited by Patrick Keast and Graeme Fairweather,
c    Reidel, 1987, pages 337-340,
c    ISBN: 9027725144,
c    LC: QA299.3.N38.
c
c  Parameters:
c
c    Input, integer M, the dimension of the argument.
c
c    Input, double precision C(M), W(M), the problem parameters.
c
c    Output, double precision Q, the integral.
c
      implicit none

      integer m

      double precision c(m)
      double precision c_prod
      integer i4vec_sum
      integer ivec(m)
      double precision q
      integer rank
      double precision r8_factorial
      double precision r8_mop
      double precision r8vec_i4vec_dot_product
      double precision r8vec_product
      integer s
      double precision w(m)
c
c  Here, we need to generate all possible DIM_NUM tuples with
c  values of 0 or 1.
c 
      q = 0.0D+00
      rank = 0

10    continue

        call tuple_next ( 0, 1, m, rank, ivec )

        if ( rank .eq. 0 ) then
          go to 20
        end if

        s = i4vec_sum ( m, ivec )

        q = q + r8_mop ( s ) 
     &    / ( 1.0D+00 + r8vec_i4vec_dot_product ( m, c, ivec ) )

      go to 10

20    continue

      c_prod = r8vec_product ( m, c )

      q = q / ( r8_factorial ( m ) * c_prod )

      return
      end
      subroutine p03_title ( title )

c*********************************************************************72
c
cc P03_TITLE returns the title of problem p03.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    28 August 2012
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, character * ( * ) TITLE, the title of the problem.
c
      implicit none

      character * ( * ) title

      title = 'Corner Peak'

      return
      end
      subroutine p04_d ( m, id, c, w, n, x, d )

c*********************************************************************72
c
cc P04_D evaluates any derivative component for problem p04.
c
c  Discussion:
c
c    f(x) = exp ( - sum ( c(1:m)^2 * ( x(1:m) - w(1:m) )^2 )
c
c    Default values are:
c
c    c(1:m) = 1 / m
c    w(1:m) = 0.5
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    28 August 2012
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Alan Genz,
c    A Package for Testing Multiple Integration Subroutines,
c    in Numerical Integration: Recent Developments, Software
c    and Applications,
c    edited by Patrick Keast and Graeme Fairweather,
c    Reidel, 1987, pages 337-340,
c    ISBN: 9027725144,
c    LC: QA299.3.N38.
c
c  Parameters:
c
c    Input, integer M, the dimension of the argument.
c
c    Input, integer ID, the spatial coordinate to differentiate.
c
c    Input, double precision C(M), W(M), the problem parameters.
c
c    Input, integer N, the number of points.
c
c    Input, double precision X(M,N), the evaluation points.
c
c    Output, double precision D(N), the ID-th derivative component.
c
      implicit none

      integer m
      integer n

      double precision c(m)
      double precision d(n)
      integer i
      integer id
      integer j
      double precision w(m)
      double precision x(m,n)

      do j = 1, n
        d(j) = 0.0D+00
      end do

      do i = 1, m
        do j = 1, n
          d(j) = d(j) + ( c(i) * ( x(i,j) - w(i) ) )**2
        end do
      end do

      do j = 1, n
        d(j) = exp ( - d(j) ) * c(id)**2 * ( - 2.0D+00 ) 
     &    * ( x(id,j) - w(id) )
      end do

      return
      end
      subroutine p04_f ( m, c, w, n, x, f )

c*********************************************************************72
c
cc P04_F evaluates the function for problem p04.
c
c  Discussion:
c
c    f(x) = exp ( - sum ( c(1:m)^2 * ( x(1:m) - w(1:m) )^2 )
c
c    Default values are:
c
c    c(1:m) = 1 / m
c    w(1:m) = 0.5
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    28 August 2012
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Alan Genz,
c    A Package for Testing Multiple Integration Subroutines,
c    in Numerical Integration: Recent Developments, Software
c    and Applications,
c    edited by Patrick Keast and Graeme Fairweather,
c    Reidel, 1987, pages 337-340,
c    ISBN: 9027725144,
c    LC: QA299.3.N38.
c
c  Parameters:
c
c    Input, integer M, the dimension of the argument.
c
c    Input, double precision C(M), W(M), the problem parameters.
c
c    Input, integer N, the number of points.
c
c    Input, double precision X(M,N), the evaluation points.
c
c    Output, double precision F(N), the function values.
c
      implicit none

      integer m
      integer n

      double precision c(m)
      double precision f(n)
      integer i
      integer j
      double precision w(m)
      double precision x(m,n)

      do j = 1, n
        f(j) = 0.0D+00
      end do

      do i = 1, m
        do j = 1, n
          f(j) = f(j) + ( c(i) * ( x(i,j) - w(i) ) )**2
        end do
      end do

      do j = 1, n
        f(j) = exp ( - f(j) )
      end do

      return
      end
      subroutine p04_q ( m, c, w, q )

c*********************************************************************72
c
cc P04_Q evaluates the integral for problem p04.
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
c    John Burkardt
c
c  Reference:
c
c    Alan Genz,
c    A Package for Testing Multiple Integration Subroutines,
c    in Numerical Integration: Recent Developments, Software
c    and Applications,
c    edited by Patrick Keast and Graeme Fairweather,
c    Reidel, 1987, pages 337-340,
c    ISBN: 9027725144,
c    LC: QA299.3.N38.
c
c  Parameters:
c
c    Input, integer M, the dimension of the argument.
c
c    Input, double precision C(M), W(M), the problem parameters.
c
c    Output, double precision Q, the integral.
c
      implicit none

      integer m

      double precision c(m)
      integer i
      double precision pi
      parameter ( pi = 3.141592653589793D+00 )
      double precision q
      double precision r8_error
      double precision w(m)

      q = 1.0D+00

      do i = 1, m

        q = q * sqrt ( pi ) 
     &    * ( r8_error ( c(i) * ( 1.0D+00 - w(i) ) ) 
     &      + r8_error ( c(i) *             w(i) ) ) 
     &    / ( 2.0D+00 * c(i) )

      end do

      return
      end
      subroutine p04_title ( title )

c*********************************************************************72
c
cc P04_TITLE returns the title of problem p04.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    28 August 2012
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, character * ( * ) TITLE, the title of the problem.
c
      implicit none

      character * ( * ) title

      title = 'Gaussian'

      return
      end
      subroutine p05_d ( m, id, c, w, n, x, d )

c*********************************************************************72
c
cc P05_D evaluates any derivative component for problem p05.
c
c  Discussion:
c
c    f(x) = exp ( - sum ( c(1:m) * abs ( x(1:m) - w(1:m) ) ) )
c
c    Default values are:
c
c    c(1:m) = 2.0
c    w(1:m) = 0.5
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    28 August 2012
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Alan Genz,
c    A Package for Testing Multiple Integration Subroutines,
c    in Numerical Integration: Recent Developments, Software
c    and Applications,
c    edited by Patrick Keast and Graeme Fairweather,
c    Reidel, 1987, pages 337-340,
c    ISBN: 9027725144,
c    LC: QA299.3.N38.
c
c  Parameters:
c
c    Input, integer M, the dimension of the argument.
c
c    Input, integer ID, the spatial coordinate to differentiate.
c
c    Input, double precision C(M), W(M), the problem parameters.
c
c    Input, integer N, the number of points.
c
c    Input, double precision X(M,N), the evaluation points.
c
c    Output, double precision D(N), the ID-th derivative component.
c
      implicit none

      integer m
      integer n

      double precision c(m)
      double precision d(n)
      integer i
      integer id
      integer j
      double precision w(m)
      double precision x(m,n)

      do j = 1, n
        d(j) = 0.0D+00
      end do

      do i = 1, m
        do j = 1, n
          d(j) = d(j) + c(i) * abs ( x(i,j) - w(i) )
        end do
      end do

      do j = 1, n
        d(1:n) = exp ( - d(j) )
      end do

      do j = 1, n
        if ( x(id,j) - w(id) .le. 0.0D+00 ) then
          d(j) = d(j) * c(id)
        else
          d(j) = - d(j) * c(id)
        end if
      end do

      return
      end
      subroutine p05_f ( m, c, w, n, x, f )

c*********************************************************************72
c
cc P05_F evaluates the function for problem p05.
c
c  Discussion:
c
c    f(x) = exp ( - sum ( c(1:m) * abs ( x(1:m) - w(1:m) ) ) )
c
c    Default values are:
c
c    c(1:m) = 2.0
c    w(1:m) = 0.5
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    28 August 2012
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Alan Genz,
c    A Package for Testing Multiple Integration Subroutines,
c    in Numerical Integration: Recent Developments, Software
c    and Applications,
c    edited by Patrick Keast and Graeme Fairweather,
c    Reidel, 1987, pages 337-340,
c    ISBN: 9027725144,
c    LC: QA299.3.N38.
c
c  Parameters:
c
c    Input, integer M, the dimension of the argument.
c
c    Input, double precision C(M), W(M), the problem parameters.
c
c    Input, integer N, the number of points.
c
c    Input, double precision X(M,N), the evaluation points.
c
c    Output, double precision F(N), the function values.
c
      implicit none

      integer m
      integer n

      double precision c(m)
      double precision f(n)
      integer i
      integer j
      double precision w(m)
      double precision x(m,n)

      do j = 1, n
        f(j) = 0.0D+00
      end do

      do i = 1, m
        do j = 1, n
          f(j) = f(j) + c(i) * abs ( x(i,j) - w(i) )
        end do
      end do

      do j = 1, n
        f(j) = exp ( - f(j) )
      end do

      return
      end
      subroutine p05_q ( m, c, w, q )

c*********************************************************************72
c
cc P05_Q evaluates the integral for problem p05.
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
c    John Burkardt
c
c  Reference:
c
c    Alan Genz,
c    A Package for Testing Multiple Integration Subroutines,
c    in Numerical Integration: Recent Developments, Software
c    and Applications,
c    edited by Patrick Keast and Graeme Fairweather,
c    Reidel, 1987, pages 337-340,
c    ISBN: 9027725144,
c    LC: QA299.3.N38.
c
c  Parameters:
c
c    Input, integer M, the dimension of the argument.
c
c    Input, double precision C(M), W(M), the problem parameters.
c
c    Output, double precision Q, the integral.
c
      implicit none

      integer m

      double precision c(m)
      integer i
      double precision q
      double precision w(m)

      q = 1.0D+00

      do i = 1, m
c
c  W < 0 < 1
c
c  | X - W | = X - W from 0 to 1.
c
        if ( w(i) .lt. 0.0D+00 ) then

          q = q * 
     &      ( exp ( - c(i) * (         - w(i) ) ) 
     &      - exp ( - c(i) * ( 1.0D+00 - w(i) ) ) ) / c(i)
c
c  0 < W < 1
c
c  | X - W | = W - X from 0 to Z, 
c            = X - W from      Z to 1.
c
        else if ( w(i) .lt. 1.0D+00 ) then

          q = q * ( 2.0D+00 
     &        - exp ( - c(i) * (           w(i) ) ) 
     &        - exp ( - c(i) * ( 1.0D+00 - w(i) ) ) ) / c(i)
c
c  0 < 1 < W
c
c  | X - W | = W - X from 0 to 1.
c
        else

          q = q * 
     &      ( exp ( - c(i) * ( w(i) - 1.0D+00 ) ) 
     &      - exp ( - c(i) * ( w(i)           ) ) ) / c(i)

        end if

      end do

      return
      end
      subroutine p05_title ( title )

c*********************************************************************72
c
cc P05_TITLE returns the title of problem p05.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    28 August 2012
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, character * ( * ) TITLE, the title of the problem.
c
      implicit none

      character * ( * ) title

      title = 'Continuous'

      return
      end
      subroutine p06_d ( m, id, c, w, n, x, d )

c*********************************************************************72
c
cc P06_D evaluates any derivative component for problem p06.
c
c  Discussion:
c
c    f(x) = exp ( c(1:m) * x(1:m) ) if x(1) <= w(1) and x(2) <= w(2).
c           0                          otherwise
c
c    Default values are:
c
c    c(1:m) = 0.5^(1/m)
c    w(1:2) = 0.5^(1/m)
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    28 August 2012
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Alan Genz,
c    A Package for Testing Multiple Integration Subroutines,
c    in Numerical Integration: Recent Developments, Software
c    and Applications,
c    edited by Patrick Keast and Graeme Fairweather,
c    Reidel, 1987, pages 337-340,
c    ISBN: 9027725144,
c    LC: QA299.3.N38.
c
c  Parameters:
c
c    Input, integer M, the dimension of the argument.
c
c    Input, integer ID, the spatial coordinate to differentiate.
c
c    Input, double precision C(M), W(M), the problem parameters.
c
c    Input, integer N, the number of points.
c
c    Input, double precision X(M,N), the evaluation points.
c
c    Output, double precision D(N), the ID-th derivative component.
c
      implicit none

      integer m
      integer n

      double precision c(m)
      double precision d(n)
      integer i
      integer id
      integer j
      double precision w(m)
      double precision x(m,n)

      if ( m .eq. 1 ) then

        do j = 1, n
          d(j) = c(1) * exp ( c(1) * x(1,j) )
        end do

        do j = 1, n
          if ( w(1) .lt. x(1,j) ) then
            d(j) = 0.0D+00
          end if
        end do

      else

        do j = 1, n
          d(j) = 0.0D+00
        end do

        do i = 1, m
          do j = 1, n
            d(j) = d(j) + c(i) * x(i,j)
          end do
        end do

        do j = 1, n
          d(j) = c(id) * exp ( d(j) )
        end do

        do j = 1, n
          if ( w(1) .lt. x(1,j) .or. w(2) .lt. x(2,j) ) then
            d(j) = 0.0D+00
          end if
        end do

      end if

      return
      end
      subroutine p06_f ( m, c, w, n, x, f )

c*********************************************************************72
c
cc P06_F evaluates the function for problem p06.
c
c  Discussion:
c
c    f(x) = exp ( c(1:m) * x(1:m) ) if x(1) <= w(1) and x(2) <= w(2).
c           0                          otherwise
c
c    Default values are:
c
c    c(1:m) = 0.5^(1/m)
c    w(1:2) = 0.5^(1/m)
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    28 August 2012
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Alan Genz,
c    A Package for Testing Multiple Integration Subroutines,
c    in Numerical Integration: Recent Developments, Software
c    and Applications,
c    edited by Patrick Keast and Graeme Fairweather,
c    Reidel, 1987, pages 337-340,
c    ISBN: 9027725144,
c    LC: QA299.3.N38.
c
c  Parameters:
c
c    Input, integer M, the dimension of the argument.
c
c    Input, double precision C(M), W(M), the problem parameters.
c
c    Input, integer N, the number of points.
c
c    Input, double precision X(M,N), the evaluation points.
c
c    Output, double precision F(N), the function values.
c
      implicit none

      integer m
      integer n

      double precision c(m)
      double precision f(n)
      integer i
      integer j
      double precision w(m)
      double precision x(m,n)

      if ( m .eq. 1 ) then

        do j = 1, n
          f(j) = exp ( c(1) * x(1,j) )
        end do

        do j = 1, n
          if ( w(1) .lt. x(1,j) ) then
            f(j) = 0.0D+00
          end if
        end do

      else

        do j = 1, n
          f(j) = 0.0D+00
        end do

        do i = 1, m
          do j = 1, n
            f(j) = f(j) + c(i) * x(i,j)
          end do
        end do

        do j = 1, n
          f(j) = exp ( f(j) )
        end do

        do j = 1, n
          if ( w(1) .lt. x(1,j) .or. w(2) .lt. x(2,j) ) then
            f(j) = 0.0D+00
          end if
        end do

      end if

      return
      end
      subroutine p06_q ( m, c, w, q )

c*********************************************************************72
c
cc P06_Q evaluates the integral for problem p06.
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
c    John Burkardt
c
c  Reference:
c
c    Alan Genz,
c    A Package for Testing Multiple Integration Subroutines,
c    in Numerical Integration: Recent Developments, Software
c    and Applications,
c    edited by Patrick Keast and Graeme Fairweather,
c    Reidel, 1987, pages 337-340,
c    ISBN: 9027725144,
c    LC: QA299.3.N38.
c
c  Parameters:
c
c    Input, integer M, the dimension of the argument.
c
c    Input, double precision C(M), W(M), the problem parameters.
c
c    Output, double precision Q, the integral.
c
      implicit none

      integer m

      double precision c(m)
      integer i
      double precision q
      double precision w(m)
c
c  To simplify the calculation, force W(3:M) to be at least 1.0.
c
      do i = 3, m
        w(i) = 1.0D+00
      end do

      q = 1.0D+00

      do i = 1, m

        if ( w(i) .le. 0.0D+00 ) then

          q = q * 0.0D+00

        else if ( w(i) .le. 1.0D+00 ) then

          if ( c(i) .eq. 0.0D+00 ) then
            q = q * w(i)
          else
            q = q * ( exp ( c(i) * w(i) ) - 1.0D+00 ) / c(i)
          end if

        else

          if ( c(i) .ne. 0.0D+00 ) then
            q = q * ( exp ( c(i) * w(i) ) - 1.0D+00 ) / c(i)
          end if

        end if

      end do

      return
      end
      subroutine p06_title ( title )

c*********************************************************************72
c
cc P06_TITLE returns the title of problem p06.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    28 August 2012
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, character * ( * ) TITLE, the title of the problem.
c
      implicit none

      character * ( * ) title

      title = 'Discontinuous'

      return
      end
      function r8_error ( x )

c*********************************************************************72
c
cc R8_ERROR evaluates the error function of an R8 argument.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    12 March 2010
c
c  Author:
c
c    Original FORTRAN77 version by Wayne Fullerton.
c    This FORTRAN77 version by John Burkardt.
c
c  Reference:
c
c    Wayne Fullerton,
c    Portable Special Function Routines,
c    in Portability of Numerical Software,
c    edited by Wayne Cowell,
c    Lecture Notes in Computer Science, Volume 57,
c    Springer 1977,
c    ISBN: 978-3-540-08446-4,
c    LC: QA297.W65.
c
c  Parameters:
c
c    Input, double precision X, the argument.
c
c    Output, double precision R8_ERROR, the error function of X.
c
      implicit none

      double precision erfcs(21)
      integer nterf
      double precision csevl
      double precision r8_error
      double precision r8_errorc
      integer inits
      double precision r8_mach
      double precision sqeps
      double precision sqrtpi
      double precision value
      double precision x
      double precision xbig
      double precision y

      save erfcs
      save nterf
      save sqeps
      save sqrtpi
      save xbig

      data erfcs(  1) / -0.49046121234691808039984544033376D-01 /
      data erfcs(  2) / -0.14226120510371364237824741899631D+00 /
      data erfcs(  3) / +0.10035582187599795575754676712933D-01 /
      data erfcs(  4) / -0.57687646997674847650827025509167D-03 /
      data erfcs(  5) / +0.27419931252196061034422160791471D-04 /
      data erfcs(  6) / -0.11043175507344507604135381295905D-05 /
      data erfcs(  7) / +0.38488755420345036949961311498174D-07 /
      data erfcs(  8) / -0.11808582533875466969631751801581D-08 /
      data erfcs(  9) / +0.32334215826050909646402930953354D-10 /
      data erfcs( 10) / -0.79910159470045487581607374708595D-12 /
      data erfcs( 11) / +0.17990725113961455611967245486634D-13 /
      data erfcs( 12) / -0.37186354878186926382316828209493D-15 /
      data erfcs( 13) / +0.71035990037142529711689908394666D-17 /
      data erfcs( 14) / -0.12612455119155225832495424853333D-18 /
      data erfcs( 15) / +0.20916406941769294369170500266666D-20 /
      data erfcs( 16) / -0.32539731029314072982364160000000D-22 /
      data erfcs( 17) / +0.47668672097976748332373333333333D-24 /
      data erfcs( 18) / -0.65980120782851343155199999999999D-26 /
      data erfcs( 19) / +0.86550114699637626197333333333333D-28 /
      data erfcs( 20) / -0.10788925177498064213333333333333D-29 /
      data erfcs( 21) / +0.12811883993017002666666666666666D-31 /

      data nterf / 0 /
      data sqeps / 0.0D+00 /
      data sqrtpi / 1.77245385090551602729816748334115D+00 /
      data xbig / 0.0D+00 /

      if ( nterf .eq. 0 ) then
        nterf = inits ( erfcs, 21, 0.1D+00 * r8_mach ( 3 ) )
        xbig = dsqrt ( - dlog ( sqrtpi * r8_mach ( 3 ) ) )
        sqeps = dsqrt ( 2.0D+00 * r8_mach ( 3 ) )
      end if

      y = dabs ( x )

      if ( y .le. sqeps ) then
        value = 2.0D+00 * x / sqrtpi
      else if ( y .le. 1.0D+00 ) then
        value = x * ( 1.0D+00 
     &    + csevl ( 2.0D+00 * x * x - 1.0D+00, erfcs, nterf ) )
      else if ( y .le. xbig ) then
        value = 1.0D+00 - r8_errorc ( y )
        if ( x .lt. 0.0D+00 ) then
          value = - value
        end if
      else
        value = 1.0D+00
        if ( x .lt. 0.0D+00 ) then
          value = - value
        end if
      end if

      r8_error = value

      return
      end
      function r8_errorc ( x )

c*********************************************************************72
c
cc R8_ERRORC evaluates the co-error function of an R8 argument.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    14 March 2010
c
c  Author:
c
c    Original FORTRAN77 version by Wayne Fullerton.
c    This FORTRAN77 version by John Burkardt.
c
c  Reference:
c
c    Wayne Fullerton,
c    Portable Special Function Routines,
c    in Portability of Numerical Software,
c    edited by Wayne Cowell,
c    Lecture Notes in Computer Science, Volume 57,
c    Springer 1977,
c    ISBN: 978-3-540-08446-4,
c    LC: QA297.W65.
c
c  Parameters:
c
c    Input, double precision X, the argument.
c
c    Output, double precision R8_ERRORC, the co-error function of X.
c
      implicit none

      double precision erc2cs(49)
      double precision erfccs(59)
      double precision erfcs(21)
      double precision eta
      integer nterc2
      integer nterf
      integer nterfc
      double precision csevl
      double precision r8_errorc
      integer inits
      double precision r8_mach
      double precision sqeps
      double precision sqrtpi
      double precision x
      double precision xmax
      double precision xsml
      double precision y

      save erfccs
      save erfcs
      save erc2cs
      save nterc2
      save nterf
      save nterfc
      save sqeps
      save sqrtpi
      save xmax
      save xsml

      data erfcs(  1) / -0.49046121234691808039984544033376D-01 /
      data erfcs(  2) / -0.14226120510371364237824741899631D+00 /
      data erfcs(  3) / +0.10035582187599795575754676712933D-01 /
      data erfcs(  4) / -0.57687646997674847650827025509167D-03 /
      data erfcs(  5) / +0.27419931252196061034422160791471D-04 /
      data erfcs(  6) / -0.11043175507344507604135381295905D-05 /
      data erfcs(  7) / +0.38488755420345036949961311498174D-07 /
      data erfcs(  8) / -0.11808582533875466969631751801581D-08 /
      data erfcs(  9) / +0.32334215826050909646402930953354D-10 /
      data erfcs( 10) / -0.79910159470045487581607374708595D-12 /
      data erfcs( 11) / +0.17990725113961455611967245486634D-13 /
      data erfcs( 12) / -0.37186354878186926382316828209493D-15 /
      data erfcs( 13) / +0.71035990037142529711689908394666D-17 /
      data erfcs( 14) / -0.12612455119155225832495424853333D-18 /
      data erfcs( 15) / +0.20916406941769294369170500266666D-20 /
      data erfcs( 16) / -0.32539731029314072982364160000000D-22 /
      data erfcs( 17) / +0.47668672097976748332373333333333D-24 /
      data erfcs( 18) / -0.65980120782851343155199999999999D-26 /
      data erfcs( 19) / +0.86550114699637626197333333333333D-28 /
      data erfcs( 20) / -0.10788925177498064213333333333333D-29 /
      data erfcs( 21) / +0.12811883993017002666666666666666D-31 /

      data erc2cs(  1) / -0.6960134660230950112739150826197D-01 /
      data erc2cs(  2) / -0.4110133936262089348982212084666D-01 /
      data erc2cs(  3) / +0.3914495866689626881561143705244D-02 /
      data erc2cs(  4) / -0.4906395650548979161280935450774D-03 /
      data erc2cs(  5) / +0.7157479001377036380760894141825D-04 /
      data erc2cs(  6) / -0.1153071634131232833808232847912D-04 /
      data erc2cs(  7) / +0.1994670590201997635052314867709D-05 /
      data erc2cs(  8) / -0.3642666471599222873936118430711D-06 /
      data erc2cs(  9) / +0.6944372610005012589931277214633D-07 /
      data erc2cs( 10) / -0.1371220902104366019534605141210D-07 /
      data erc2cs( 11) / +0.2788389661007137131963860348087D-08 /
      data erc2cs( 12) / -0.5814164724331161551864791050316D-09 /
      data erc2cs( 13) / +0.1238920491752753181180168817950D-09 /
      data erc2cs( 14) / -0.2690639145306743432390424937889D-10 /
      data erc2cs( 15) / +0.5942614350847910982444709683840D-11 /
      data erc2cs( 16) / -0.1332386735758119579287754420570D-11 /
      data erc2cs( 17) / +0.3028046806177132017173697243304D-12 /
      data erc2cs( 18) / -0.6966648814941032588795867588954D-13 /
      data erc2cs( 19) / +0.1620854541053922969812893227628D-13 /
      data erc2cs( 20) / -0.3809934465250491999876913057729D-14 /
      data erc2cs( 21) / +0.9040487815978831149368971012975D-15 /
      data erc2cs( 22) / -0.2164006195089607347809812047003D-15 /
      data erc2cs( 23) / +0.5222102233995854984607980244172D-16 /
      data erc2cs( 24) / -0.1269729602364555336372415527780D-16 /
      data erc2cs( 25) / +0.3109145504276197583836227412951D-17 /
      data erc2cs( 26) / -0.7663762920320385524009566714811D-18 /
      data erc2cs( 27) / +0.1900819251362745202536929733290D-18 /
      data erc2cs( 28) / -0.4742207279069039545225655999965D-19 /
      data erc2cs( 29) / +0.1189649200076528382880683078451D-19 /
      data erc2cs( 30) / -0.3000035590325780256845271313066D-20 /
      data erc2cs( 31) / +0.7602993453043246173019385277098D-21 /
      data erc2cs( 32) / -0.1935909447606872881569811049130D-21 /
      data erc2cs( 33) / +0.4951399124773337881000042386773D-22 /
      data erc2cs( 34) / -0.1271807481336371879608621989888D-22 /
      data erc2cs( 35) / +0.3280049600469513043315841652053D-23 /
      data erc2cs( 36) / -0.8492320176822896568924792422399D-24 /
      data erc2cs( 37) / +0.2206917892807560223519879987199D-24 /
      data erc2cs( 38) / -0.5755617245696528498312819507199D-25 /
      data erc2cs( 39) / +0.1506191533639234250354144051199D-25 /
      data erc2cs( 40) / -0.3954502959018796953104285695999D-26 /
      data erc2cs( 41) / +0.1041529704151500979984645051733D-26 /
      data erc2cs( 42) / -0.2751487795278765079450178901333D-27 /
      data erc2cs( 43) / +0.7290058205497557408997703680000D-28 /
      data erc2cs( 44) / -0.1936939645915947804077501098666D-28 /
      data erc2cs( 45) / +0.5160357112051487298370054826666D-29 /
      data erc2cs( 46) / -0.1378419322193094099389644800000D-29 /
      data erc2cs( 47) / +0.3691326793107069042251093333333D-30 /
      data erc2cs( 48) / -0.9909389590624365420653226666666D-31 /
      data erc2cs( 49) / +0.2666491705195388413323946666666D-31 /

      data erfccs(  1) / +0.715179310202924774503697709496D-01 /
      data erfccs(  2) / -0.265324343376067157558893386681D-01 /
      data erfccs(  3) / +0.171115397792085588332699194606D-02 /
      data erfccs(  4) / -0.163751663458517884163746404749D-03 /
      data erfccs(  5) / +0.198712935005520364995974806758D-04 /
      data erfccs(  6) / -0.284371241276655508750175183152D-05 /
      data erfccs(  7) / +0.460616130896313036969379968464D-06 /
      data erfccs(  8) / -0.822775302587920842057766536366D-07 /
      data erfccs(  9) / +0.159214187277090112989358340826D-07 /
      data erfccs( 10) / -0.329507136225284321486631665072D-08 /
      data erfccs( 11) / +0.722343976040055546581261153890D-09 /
      data erfccs( 12) / -0.166485581339872959344695966886D-09 /
      data erfccs( 13) / +0.401039258823766482077671768814D-10 /
      data erfccs( 14) / -0.100481621442573113272170176283D-10 /
      data erfccs( 15) / +0.260827591330033380859341009439D-11 /
      data erfccs( 16) / -0.699111056040402486557697812476D-12 /
      data erfccs( 17) / +0.192949233326170708624205749803D-12 /
      data erfccs( 18) / -0.547013118875433106490125085271D-13 /
      data erfccs( 19) / +0.158966330976269744839084032762D-13 /
      data erfccs( 20) / -0.472689398019755483920369584290D-14 /
      data erfccs( 21) / +0.143587337678498478672873997840D-14 /
      data erfccs( 22) / -0.444951056181735839417250062829D-15 /
      data erfccs( 23) / +0.140481088476823343737305537466D-15 /
      data erfccs( 24) / -0.451381838776421089625963281623D-16 /
      data erfccs( 25) / +0.147452154104513307787018713262D-16 /
      data erfccs( 26) / -0.489262140694577615436841552532D-17 /
      data erfccs( 27) / +0.164761214141064673895301522827D-17 /
      data erfccs( 28) / -0.562681717632940809299928521323D-18 /
      data erfccs( 29) / +0.194744338223207851429197867821D-18 /
      data erfccs( 30) / -0.682630564294842072956664144723D-19 /
      data erfccs( 31) / +0.242198888729864924018301125438D-19 /
      data erfccs( 32) / -0.869341413350307042563800861857D-20 /
      data erfccs( 33) / +0.315518034622808557122363401262D-20 /
      data erfccs( 34) / -0.115737232404960874261239486742D-20 /
      data erfccs( 35) / +0.428894716160565394623737097442D-21 /
      data erfccs( 36) / -0.160503074205761685005737770964D-21 /
      data erfccs( 37) / +0.606329875745380264495069923027D-22 /
      data erfccs( 38) / -0.231140425169795849098840801367D-22 /
      data erfccs( 39) / +0.888877854066188552554702955697D-23 /
      data erfccs( 40) / -0.344726057665137652230718495566D-23 /
      data erfccs( 41) / +0.134786546020696506827582774181D-23 /
      data erfccs( 42) / -0.531179407112502173645873201807D-24 /
      data erfccs( 43) / +0.210934105861978316828954734537D-24 /
      data erfccs( 44) / -0.843836558792378911598133256738D-25 /
      data erfccs( 45) / +0.339998252494520890627359576337D-25 /
      data erfccs( 46) / -0.137945238807324209002238377110D-25 /
      data erfccs( 47) / +0.563449031183325261513392634811D-26 /
      data erfccs( 48) / -0.231649043447706544823427752700D-26 /
      data erfccs( 49) / +0.958446284460181015263158381226D-27 /
      data erfccs( 50) / -0.399072288033010972624224850193D-27 /
      data erfccs( 51) / +0.167212922594447736017228709669D-27 /
      data erfccs( 52) / -0.704599152276601385638803782587D-28 /
      data erfccs( 53) / +0.297976840286420635412357989444D-28 /
      data erfccs( 54) / -0.126252246646061929722422632994D-28 /
      data erfccs( 55) / +0.539543870454248793985299653154D-29 /
      data erfccs( 56) / -0.238099288253145918675346190062D-29 /
      data erfccs( 57) / +0.109905283010276157359726683750D-29 /
      data erfccs( 58) / -0.486771374164496572732518677435D-30 /
      data erfccs( 59) / +0.152587726411035756763200828211D-30 /

      data nterc2 / 0 /
      data nterf / 0 /
      data nterfc / 0 /
      data sqeps / 0.0D+00 /
      data sqrtpi / 1.77245385090551602729816748334115D+00 /
      data xmax / 0.0D+00 /
      data xsml / 0.0D+00 /

      if ( nterf .eq. 0 ) then

        eta = 0.1D+00 * r8_mach ( 3 )
        nterf = inits ( erfcs, 21, eta )
        nterfc = inits ( erfccs, 59, eta )
        nterc2 = inits ( erc2cs, 49, eta )

        xsml = - dsqrt ( - dlog ( sqrtpi * r8_mach ( 3 ) ) )
        xmax = dsqrt (- dlog ( sqrtpi * r8_mach ( 1 ) ) )
        xmax = xmax - 0.5D+00 * dlog ( xmax ) / xmax - 0.01D+00
        sqeps = dsqrt ( 2.0D+00 * r8_mach ( 3 ) )

      end if

      if ( x .le. xsml ) then

        r8_errorc = 2.0D+00
        return

      end if

      if ( xmax .lt. x ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'R8_ERRORC - Warning!'
        write ( *, '(a)' ) '  X so big that ERFC underflows.'
        r8_errorc = 0.0D+00
        return
      end if

      y = dabs ( x )

      if ( y .lt. sqeps ) then
        r8_errorc = 1.0D+00 - 2.0D+00 * x / sqrtpi
        return
      else if ( y .le. 1.0D+00 ) then
        r8_errorc = 1.0D+00 - x * ( 1.0D+00 
     &    + csevl ( 2.0D+00 * x * x - 1.0D+00, erfcs, nterf ) )
        return
      end if

      y = y * y

      if ( y .le. 4.0D+00 ) then
        r8_errorc = dexp ( - y ) / dabs ( x ) * ( 0.5D+00 
     &    + csevl ( ( 8.0D+00 / y - 5.0D+00 ) / 3.0D+00, erc2cs, 
     &    nterc2 ) )
      else 
        r8_errorc = dexp ( - y ) / dabs ( x ) * ( 0.5D+00 
     &    + csevl ( 8.0D+00 / y - 1.0D+00, erfccs, nterfc ) )
      end if

      if ( x .lt. 0.0D+00 ) then
        r8_errorc = 2.0D+00 - r8_errorc
      end if

      return
      end
      function r8_mach ( i )

c*********************************************************************72
c
cc R8_MACH returns double precision real machine-dependent constants.
c
c  Discussion:
c
c    R8_MACH can be used to obtain machine-dependent parameters
c    for the local machine environment.  It is a function
c    with one input argument, and can be called as follows:
c
c      D = R8_MACH ( I )
c
c    where I=1,...,5.  The output value of D above is
c    determined by the input value of I:.
c
c    R8_MACH ( 1) = B^(EMIN-1), the smallest positive magnitude.
c    R8_MACH ( 2) = B^EMAX*(1 - B^(-T)), the largest magnitude.
c    R8_MACH ( 3) = B^(-T), the smallest relative spacing.
c    R8_MACH ( 4) = B^(1-T), the largest relative spacing.
c    R8_MACH ( 5) = LOG10(B)
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    25 April 2007
c
c  Author:
c
c    Original FORTRAN77 version by Phyllis Fox, Andrew Hall, Norman Schryer.
c    This FORTRAN77 version by John Burkardt.
c
c  Reference:
c
c    Phyllis Fox, Andrew Hall, Norman Schryer,
c    Algorithm 528:
c    Framework for a Portable Library,
c    ACM Transactions on Mathematical Software,
c    Volume 4, Number 2, June 1978, page 176-188.
c
c  Parameters:
c
c    Input, integer I, the index of the desired constant.
c
c    Output, double precision R8_MACH, the value of the constant.
c
      implicit none

      double precision r8_mach
      integer i

      if ( i .lt. 1 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'R8_MACH - Fatal error!'
        write ( *, '(a)' ) '  The input argument I is out of bounds.'
        write ( *, '(a)' ) '  Legal values satisfy 1 <= I <= 5.'
        write ( *, '(a,i12)' ) '  I = ', i
        r8_mach = 0.0D+00
        stop
      else if ( i .eq. 1 ) then
        r8_mach = 4.450147717014403D-308
      else if ( i .eq. 2 ) then
        r8_mach = 8.988465674311579D+307
      else if ( i .eq. 3 ) then
        r8_mach = 1.110223024625157D-016
      else if ( i .eq. 4 ) then
        r8_mach = 2.220446049250313D-016
      else if ( i .eq. 5 ) then
        r8_mach = 0.301029995663981D+000
      else if ( 5 .lt. i ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'R8_MACH - Fatal error!'
        write ( *, '(a)' ) '  The input argument I is out of bounds.'
        write ( *, '(a)' ) '  Legal values satisfy 1 <= I <= 5.'
        write ( *, '(a,i12)' ) '  I = ', i
        r8_mach = 0.0D+00
        stop
      end if

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
