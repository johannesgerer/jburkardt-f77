      subroutine basis_function_b_val ( tdata, tval, yval )

c*********************************************************************72
c
cc BASIS_FUNCTION_B_VAL evaluates the B spline basis function.
c
c  Discussion:
c
c    The B spline basis function is a piecewise cubic which
c    has the properties that:
c
c    * it equals 2/3 at TDATA(3), 1/6 at TDATA(2) and TDATA(4);
c    * it is 0 for TVAL .le. TDATA(1) or TDATA(5) .le. TVAL;
c    * it is strictly increasing from TDATA(1) to TDATA(3),
c      and strictly decreasing from TDATA(3) to TDATA(5);
c    * the function and its first two derivatives are continuous
c      at each node TDATA(I).
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    06 April 1999
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Alan Davies, Philip Samuels,
c    An Introduction to Computational Geometry for Curves and Surfaces,
c    Clarendon Press, 1996,
c    ISBN: 0-19-851478-6,
c    LC: QA448.D38.
c
c  Parameters:
c
c    Input, double precision TDATA(5), the nodes associated with the
c    basis function.  The entries of TDATA are assumed to be distinct
c    and increasing.
c
c    Input, double precision TVAL, a point at which the B spline basis
c    function is to be evaluated.
c
c    Output, double precision YVAL, the value of the function at TVAL.
c
      implicit none

      integer ndata
      parameter ( ndata = 5 )

      integer left
      integer right
      double precision tdata(ndata)
      double precision tval
      double precision u
      double precision yval

      if ( tval .le. tdata(1) .or. tdata(ndata) .le. tval ) then
        yval = 0.0D+00
        return
      end if
c
c  Find the interval [ TDATA(LEFT), TDATA(RIGHT) ] containing TVAL.
c
      call r8vec_bracket ( ndata, tdata, tval, left, right )
c
c  U is the normalized coordinate of TVAL in this interval.
c
      u = ( tval - tdata(left) ) / ( tdata(right) - tdata(left) )
c
c  Now evaluate the function.
c
      if ( tval .lt. tdata(2) ) then
        yval = u**3 / 6.0D+00
      else if ( tval .lt. tdata(3) ) then
        yval = ( ( (    - 3.0D+00   
     &              * u + 3.0D+00 ) 
     &              * u + 3.0D+00 ) 
     &              * u + 1.0D+00 ) / 6.0D+00
      else if ( tval .lt. tdata(4) ) then
        yval = ( ( (    + 3.0D+00   
     &              * u - 6.0D+00 ) 
     &              * u + 0.0D+00 ) 
     &              * u + 4.0D+00 ) / 6.0D+00
      else if ( tval .lt. tdata(5) ) then
        yval = ( 1.0D+00 - u )**3 / 6.0D+00
      end if

      return
      end
      subroutine basis_function_beta_val ( beta1, beta2, tdata, tval, 
     &  yval )

c*********************************************************************72
c
cc BASIS_FUNCTION_BETA_VAL evaluates the beta spline basis function.
c
c  Discussion:
c
c    With BETA1 = 1 and BETA2 = 0, the beta spline basis function
c    equals the B spline basis function.
c
c    With BETA1 large, and BETA2 = 0, the beta spline basis function
c    skews to the right, that is, its maximum increases, and occurs
c    to the right of the center.
c
c    With BETA1 = 1 and BETA2 large, the beta spline becomes more like
c    a linear basis function; that is, its value in the outer two intervals
c    goes to zero, and its behavior in the inner two intervals approaches
c    a piecewise linear function that is 0 at the second node, 1 at the
c    third, and 0 at the fourth.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    09 April 1999
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Alan Davies, Philip Samuels,
c    An Introduction to Computational Geometry for Curves and Surfaces,
c    Clarendon Press, 1996,
c    ISBN: 0-19-851478-6,
c    LC: QA448.D38.
c
c  Parameters:
c
c    Input, double precision BETA1, the skew or bias parameter.
c    BETA1 = 1 for no skew or bias.
c
c    Input, double precision BETA2, the tension parameter.
c    BETA2 = 0 for no tension.
c
c    Input, double precision TDATA(5), the nodes associated with the
c    basis function.  The entries of TDATA are assumed to be distinct
c    and increasing.
c
c    Input, double precision TVAL, a point at which the B spline
c    basis function is to be evaluated.
c
c    Output, double precision YVAL, the value of the function at TVAL.
c
      implicit none

      integer ndata
      parameter ( ndata = 5 )

      double precision a
      double precision b
      double precision beta1
      double precision beta2
      double precision c
      double precision d
      integer left
      integer right
      double precision tdata(ndata)
      double precision tval
      double precision u
      double precision yval

      if ( tval .le. tdata(1) .or. tdata(ndata) .le. tval ) then
        yval = 0.0D+00
        return
      end if
c
c  Find the interval [ TDATA(LEFT), TDATA(RIGHT) ] containing TVAL.
c
      call r8vec_bracket ( ndata, tdata, tval, left, right )
c
c  U is the normalized coordinate of TVAL in this interval.
c
      u = ( tval - tdata(left) ) / ( tdata(right) - tdata(left) )
c
c  Now evaluate the function.
c
      if ( tval .lt. tdata(2) ) then

        yval = 2.0D+00 * u**3

      else if ( tval .lt. tdata(3) ) then

        a = beta2 + 4.0D+00 * beta1 + 4.0D+00 * beta1 * beta1 
     &    + 6.0D+00 * ( 1.0D+00 - beta1 * beta1 ) 
     &    - 3.0D+00 * ( 2.0D+00 + beta2 + 2.0D+00 * beta1 ) 
     &    + 2.0D+00 * ( 1.0D+00 + beta2 + beta1 + beta1 * beta1 )

        b = - 6.0D+00 * ( 1.0D+00 - beta1 * beta1 ) 
     &      + 6.0D+00 * ( 2.0D+00 + beta2 + 2.0D+00 * beta1 ) 
     &      - 6.0D+00 * ( 1.0D+00 + beta2 + beta1 + beta1 * beta1 )

        c = - 3.0D+00 * ( 2.0D+00 + beta2 + 2.0D+00 * beta1 ) 
     &      + 6.0D+00 * ( 1.0D+00 + beta2 + beta1 + beta1 * beta1 )

        d = - 2.0D+00 * ( 1.0D+00 + beta2 + beta1 + beta1 * beta1 )

        yval = ( ( d * u + c ) * u + b ) * u + a

      else if ( tval .lt. tdata(4) ) then

        a = beta2 + 4.0D+00 * beta1 + 4.0D+00 * beta1 * beta1

        b = - 6.0D+00 * beta1 * ( 1.0D+00 - beta1 * beta1 )

        c = - 3.0D+00 * ( beta2 + 2.0D+00 * beta1**2 
     &    + 2.0D+00 * beta1**3 )

        d = 2.0D+00 * ( beta2 + beta1 + beta1**2 + beta1**3 )

        yval = ( ( d * u + c ) * u + b ) * u + a

      else if ( tval .lt. tdata(5) ) then

        yval = 2.0D+00 * beta1**3 * ( 1.0D+00 - u )**3

      end if

      yval = yval / ( 2.0D+00 + beta2 + 4.0D+00 * beta1 
     &  + 4.0D+00 * beta1**2 + 2.0D+00 * beta1**3 )

      return
      end
      subroutine basis_matrix_b_uni ( mbasis )

c*********************************************************************72
c
cc BASIS_MATRIX_B_UNI sets up the uniform B spline basis matrix.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    17 June 2006
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    James Foley, Andries vanDam, Steven Feiner, John Hughes,
c    Computer Graphics, Principles and Practice,
c    Second Edition,
c    Addison Wesley, 1995,
c    ISBN: 0201848406,
c    LC: T385.C5735.
c
c  Parameters:
c
c    Output, double precision MBASIS(4,4), the basis matrix.
c
      implicit none

      double precision mbasis(4,4)

      mbasis(1,1) = -1.0D+00 / 6.0D+00
      mbasis(2,1) =  3.0D+00 / 6.0D+00
      mbasis(3,1) = -3.0D+00 / 6.0D+00
      mbasis(4,1) =  1.0D+00 / 6.0D+00

      mbasis(1,2) =  3.0D+00 / 6.0D+00
      mbasis(2,2) = -6.0D+00 / 6.0D+00
      mbasis(3,2) =  0.0D+00 / 6.0D+00
      mbasis(4,2) =  4.0D+00 / 6.0D+00

      mbasis(1,3) = -3.0D+00 / 6.0D+00
      mbasis(2,3) =  3.0D+00 / 6.0D+00
      mbasis(3,3) =  3.0D+00 / 6.0D+00
      mbasis(4,3) =  1.0D+00 / 6.0D+00

      mbasis(1,4) = -1.0D+00 / 6.0D+00
      mbasis(2,3) =  0.0D+00 / 6.0D+00
      mbasis(3,4) =  0.0D+00 / 6.0D+00
      mbasis(4,4) =  0.0D+00 / 6.0D+00

      return
      end
      subroutine basis_matrix_beta_uni ( beta1, beta2, mbasis )

c*********************************************************************72
c
cc BASIS_MATRIX_BETA_UNI sets up the uniform beta spline basis matrix.
c
c  Discussion:
c
c    If BETA1 = 1 and BETA2 = 0, then the beta spline reduces to
c    the B spline.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    12 February 2004
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    James Foley, Andries vanDam, Steven Feiner, John Hughes,
c    Computer Graphics, Principles and Practice,
c    Second Edition,
c    Addison Wesley, 1995,
c    ISBN: 0201848406,
c    LC: T385.C5735.
c
c  Parameters:
c
c    Input, double precision BETA1, the skew or bias parameter.
c    BETA1 = 1 for no skew or bias.
c
c    Input, double precision BETA2, the tension parameter.
c    BETA2 = 0 for no tension.
c
c    Output, double precision MBASIS(4,4), the basis matrix.
c
      implicit none

      double precision beta1
      double precision beta2
      double precision delta
      integer i
      integer j
      double precision mbasis(4,4)

      mbasis(1,1) = - 2.0D+00 * beta1 * beta1 * beta1
      mbasis(1,2) =   2.0D+00 * beta2 
     &  + 2.0 * beta1 * ( beta1 * beta1 + beta1 + 1.0D+00 )
      mbasis(1,3) = - 2.0D+00 
     &  * ( beta2 + beta1 * beta1 + beta1 + 1.0D+00 )
      mbasis(1,4) =   2.0D+00

      mbasis(2,1) =   6.0D+00 * beta1 * beta1 * beta1
      mbasis(2,2) = - 3.0D+00 * beta2 
     &  - 6.0D+00 * beta1 * beta1 * ( beta1 + 1.0D+00 )
      mbasis(2,3) =   3.0D+00 * beta2 + 6.0D+00 * beta1 * beta1
      mbasis(2,4) =   0.0D+00

      mbasis(3,1) = - 6.0D+00 * beta1 * beta1 * beta1
      mbasis(3,2) =   6.0D+00 * beta1 * ( beta1 - 1.0D+00 ) 
     &  * ( beta1 + 1.0D+00 )
      mbasis(3,3) =   6.0D+00 * beta1
      mbasis(3,4) =   0.0D+00

      mbasis(4,1) =   2.0D+00 * beta1 * beta1 * beta1
      mbasis(4,2) =   4.0D+00 * beta1 * ( beta1 + 1.0D+00 ) + beta2
      mbasis(4,3) =   2.0D+00
      mbasis(4,4) =   0.0D+00

      delta = ( ( 2.0D+00   
     &  * beta1 + 4.0D+00 ) 
     &  * beta1 + 4.0D+00 ) 
     &  * beta1 + 2.0D+00 + beta2

      do j = 1, 4
        do i = 1, 4
          mbasis(i,j) = mbasis(i,j) / delta
        end do
      end do

      return
      end
      subroutine basis_matrix_bezier ( mbasis )

c*********************************************************************72
c
cc BASIS_MATRIX_BEZIER sets up the cubic Bezier spline basis matrix.
c
c  Discussion:
c
c    This basis matrix assumes that the data points are stored as
c    ( P1, P2, P3, P4 ).  P1 is the function value at T = 0, while
c    P2 is used to approximate the derivative at T = 0 by
c    dP/dt = 3 * ( P2 - P1 ).  Similarly, P4 is the function value
c    at T = 1, and P3 is used to approximate the derivative at T = 1
c    by dP/dT = 3 * ( P4 - P3 ).
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    05 April 1999
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    James Foley, Andries vanDam, Steven Feiner, John Hughes,
c    Computer Graphics, Principles and Practice,
c    Second Edition,
c    Addison Wesley, 1995,
c    ISBN: 0201848406,
c    LC: T385.C5735.
c
c  Parameters:
c
c    Output, double precision MBASIS(4,4), the basis matrix.
c
      implicit none

      double precision mbasis(4,4)

      mbasis(1,1) = -1.0D+00
      mbasis(1,2) =  3.0D+00
      mbasis(1,3) = -3.0D+00
      mbasis(1,4) =  1.0D+00

      mbasis(2,1) =  3.0D+00
      mbasis(2,2) = -6.0D+00
      mbasis(2,3) =  3.0D+00
      mbasis(2,4) =  0.0D+00

      mbasis(3,1) = -3.0D+00
      mbasis(3,2) =  3.0D+00
      mbasis(3,3) =  0.0D+00
      mbasis(3,4) =  0.0D+00

      mbasis(4,1) =  1.0D+00
      mbasis(4,2) =  0.0D+00
      mbasis(4,3) =  0.0D+00
      mbasis(4,4) =  0.0D+00

      return
      end
      subroutine basis_matrix_hermite ( mbasis )

c*********************************************************************72
c
cc BASIS_MATRIX_HERMITE sets up the Hermite spline basis matrix.
c
c  Discussion:
c
c    This basis matrix assumes that the data points are stored as
c    ( P1, P2, P1', P2' ), with P1 and P1' being the data value and
c    the derivative dP/dT at T = 0, while P2 and P2' apply at T = 1.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    06 April 1999
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    James Foley, Andries vanDam, Steven Feiner, John Hughes,
c    Computer Graphics, Principles and Practice,
c    Second Edition,
c    Addison Wesley, 1995,
c    ISBN: 0201848406,
c    LC: T385.C5735.
c
c  Parameters:
c
c    Output, double precision MBASIS(4,4), the basis matrix.
c
      implicit none

      double precision mbasis(4,4)

      mbasis(1,1) =  2.0D+00
      mbasis(1,2) = -2.0D+00
      mbasis(1,3) =  1.0D+00
      mbasis(1,4) =  1.0D+00

      mbasis(2,1) = -3.0D+00
      mbasis(2,2) =  3.0D+00
      mbasis(2,3) = -2.0D+00
      mbasis(2,4) = -1.0D+00

      mbasis(3,1) =  0.0D+00
      mbasis(3,2) =  0.0D+00
      mbasis(3,3) =  1.0D+00
      mbasis(3,4) =  0.0D+00

      mbasis(4,1) =  1.0D+00
      mbasis(4,2) =  0.0D+00
      mbasis(4,3) =  0.0D+00
      mbasis(4,4) =  0.0D+00

      return
      end
      subroutine basis_matrix_overhauser_nonuni ( alpha, beta, mbasis )

c*********************************************************************72
c
cc BASIS_MATRIX_OVERHAUSER_NONUNI: nonuniform Overhauser spline basis matrix.
c
c  Discussion:
c
c    This basis matrix assumes that the data points P1, P2, P3 and
c    P4 are not uniformly spaced in T, and that P2 corresponds to T = 0,
c    and P3 to T = 1.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    05 April 1999
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, double precision ALPHA, BETA.
c    ALPHA = || P2 - P1 || / ( || P3 - P2 || + || P2 - P1 || )
c    BETA  = || P3 - P2 || / ( || P4 - P3 || + || P3 - P2 || ).
c
c    Output, double precision MBASIS(4,4), the basis matrix.
c
      implicit none

      double precision alpha
      double precision beta
      double precision mbasis(4,4)

      mbasis(1,1) = - ( 1.0D+00 - alpha ) * ( 1.0D+00 - alpha ) / alpha
      mbasis(1,2) =   beta + ( 1.0D+00 - alpha ) / alpha
      mbasis(1,3) =   alpha - 1.0D+00 / ( 1.0D+00 - beta )
      mbasis(1,4) =   beta * beta / ( 1.0D+00 - beta )

      mbasis(2,1) =   2.0D+00 * ( 1.0D+00 - alpha ) 
     &  * ( 1.0D+00 - alpha ) / alpha
      mbasis(2,2) = ( - 2.0D+00 * ( 1.0D+00 - alpha ) - alpha * beta ) 
     &  / alpha
      mbasis(2,3) = ( 2.0D+00 * ( 1.0D+00 - alpha ) 
     &  - beta * ( 1.0D+00 - 2.0D+00 * alpha ) ) / ( 1.0D+00 - beta )
      mbasis(2,4) = - beta * beta / ( 1.0D+00 - beta )

      mbasis(3,1) = - ( 1.0D+00 - alpha ) * ( 1.0D+00 - alpha ) / alpha
      mbasis(3,2) =   ( 1.0D+00 - 2.0D+00 * alpha ) / alpha
      mbasis(3,3) =   alpha
      mbasis(3,4) =   0.0D+00

      mbasis(4,1) =   0.0D+00
      mbasis(4,2) =   1.0D+00
      mbasis(4,3) =   0.0D+00
      mbasis(4,4) =   0.0D+00

      return
      end
      subroutine basis_matrix_overhauser_nul ( alpha, mbasis )

c*********************************************************************72
c
cc BASIS_MATRIX_OVERHAUSER_NUL sets the nonuniform left Overhauser basis matrix.
c
c  Discussion:
c
c    This basis matrix assumes that the data points P1, P2, and
c    P3 are not uniformly spaced in T, and that P1 corresponds to T = 0,
c    and P2 to T = 1.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    27 August 1999
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, double precision ALPHA.
c    ALPHA = || P2 - P1 || / ( || P3 - P2 || + || P2 - P1 || )
c
c    Output, double precision MBASIS(3,3), the basis matrix.
c
      implicit none

      double precision alpha
      double precision mbasis(3,3)

      mbasis(1,1) =   1.0D+00 / alpha
      mbasis(1,2) = - 1.0D+00 / ( alpha * ( 1.0D+00 - alpha ) )
      mbasis(1,3) =   1.0D+00 / ( 1.0D+00 - alpha )

      mbasis(2,1) = - ( 1.0D+00 + alpha ) / alpha
      mbasis(2,2) =   1.0D+00 / ( alpha * ( 1.0D+00 - alpha ) )
      mbasis(2,3) = - alpha / ( 1.0D+00 - alpha )

      mbasis(3,1) =   1.0D+00
      mbasis(3,2) =   0.0D+00
      mbasis(3,3) =   0.0D+00

      return
      end
      subroutine basis_matrix_overhauser_nur ( beta, mbasis )

c*********************************************************************72
c
cc BASIS_MATRIX_OVERHAUSER_NUR: the nonuniform right Overhauser basis matrix.
c
c  Discussion:
c
c    This basis matrix assumes that the data points PN-2, PN-1, and
c    PN are not uniformly spaced in T, and that PN-1 corresponds to T = 0,
c    and PN to T = 1.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    27 August 1999
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, double precision BETA.
c    BETA = || P(N) - P(N-1) || 
c         / ( || P(N) - P(N-1) || + || P(N-1) - P(N-2) || )
c
c    Output, double precision MBASIS(3,3), the basis matrix.
c
      implicit none

      double precision beta
      double precision mbasis(3,3)

      mbasis(1,1) =   1.0D+00 / beta
      mbasis(1,2) = - 1.0D+00 / ( beta * ( 1.0D+00 - beta ) )
      mbasis(1,3) =   1.0D+00 / ( 1.0D+00 - beta )

      mbasis(2,1) = - ( 1.0D+00 + beta ) / beta
      mbasis(2,2) =   1.0D+00 / ( beta * ( 1.0D+00 - beta ) )
      mbasis(2,3) = - beta / ( 1.0D+00 - beta )

      mbasis(3,1) =   1.0D+00
      mbasis(3,2) =   0.0D+00
      mbasis(3,3) =   0.0D+00

      return
      end
      subroutine basis_matrix_overhauser_uni ( mbasis )

c*********************************************************************72
c
cc BASIS_MATRIX_OVERHAUSER_UNI sets the uniform Overhauser spline basis matrix.
c
c  Discussion:
c
c    This basis matrix assumes that the data points P1, P2, P3 and
c    P4 are uniformly spaced in T, and that P2 corresponds to T = 0,
c    and P3 to T = 1.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    05 April 1999
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    James Foley, Andries vanDam, Steven Feiner, John Hughes,
c    Computer Graphics, Principles and Practice,
c    Second Edition,
c    Addison Wesley, 1995,
c    ISBN: 0201848406,
c    LC: T385.C5735.
c
c  Parameters:
c
c    Output, double precision MBASIS(4,4), the basis matrix.
c
      implicit none

      double precision mbasis(4,4)

      mbasis(1,1) = - 1.0D+00 / 2.0D+00
      mbasis(1,2) =   3.0D+00 / 2.0D+00
      mbasis(1,3) = - 3.0D+00 / 2.0D+00
      mbasis(1,4) =   1.0D+00 / 2.0D+00

      mbasis(2,1) =   2.0D+00 / 2.0D+00
      mbasis(2,2) = - 5.0D+00 / 2.0D+00
      mbasis(2,3) =   4.0D+00 / 2.0D+00
      mbasis(2,4) = - 1.0D+00 / 2.0D+00

      mbasis(3,1) = - 1.0D+00 / 2.0D+00
      mbasis(3,2) =   0.0D+00
      mbasis(3,3) =   1.0D+00 / 2.0D+00
      mbasis(3,4) =   0.0D+00

      mbasis(4,1) =   0.0D+00
      mbasis(4,2) =   2.0D+00 / 2.0D+00
      mbasis(4,3) =   0.0D+00
      mbasis(4,4) =   0.0D+00

      return
      end
      subroutine basis_matrix_overhauser_uni_l ( mbasis )

c*********************************************************************72
c
cc BASIS_MATRIX_OVERHAUSER_UNI_L sets the left uniform Overhauser basis matrix.
c
c  Discussion:
c
c    This basis matrix assumes that the data points P1, P2, and P3
c    are not uniformly spaced in T, and that P1 corresponds to T = 0,
c    and P2 to T = 1.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    05 April 1999
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, double precision MBASIS(3,3), the basis matrix.
c
      implicit none

      double precision mbasis(3,3)

      mbasis(1,1) =   2.0D+00
      mbasis(1,2) = - 4.0D+00
      mbasis(1,3) =   2.0D+00

      mbasis(2,1) = - 3.0D+00
      mbasis(2,2) =   4.0D+00
      mbasis(2,3) = - 1.0D+00

      mbasis(3,1) =   1.0D+00
      mbasis(3,2) =   0.0D+00
      mbasis(3,3) =   0.0D+00

      return
      end
      subroutine basis_matrix_overhauser_uni_r ( mbasis )

c*********************************************************************72
c
cc BASIS_MATRIX_OVERHAUSER_UNI_R sets the right uniform Overhauser basis matrix.
c
c  Discussion:
c
c    This basis matrix assumes that the data points P(N-2), P(N-1),
c    and P(N) are uniformly spaced in T, and that P(N-1) corresponds to
c    T = 0, and P(N) to T = 1.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    05 April 1999
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, double precision MBASIS(3,3), the basis matrix.
c
      implicit none

      double precision mbasis(3,3)

      mbasis(1,1) =   2.0D+00
      mbasis(1,2) = - 4.0D+00
      mbasis(1,3) =   2.0D+00

      mbasis(2,1) = - 3.0D+00
      mbasis(2,2) =   4.0D+00
      mbasis(2,3) = - 1.0D+00

      mbasis(3,1) =   1.0D+00
      mbasis(3,2) =   0.0D+00
      mbasis(3,3) =   0.0D+00

      return
      end
      subroutine basis_matrix_tmp ( left, n, mbasis, ndata, tdata, 
     &  ydata, tval, yval )

c*********************************************************************72
c
cc BASIS_MATRIX_TMP computes Q = T * MBASIS * P
c
c  Discussion:
c
c    YDATA is a vector of data values, most frequently the values of some
c    function sampled at uniformly spaced points.  MBASIS is the basis
c    matrix for a particular kind of spline.  T is a vector of the
c    powers of the normalized difference between TVAL and the left
c    endpoint of the interval.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    10 February 2004
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer LEFT, indicats that TVAL is in the interval
c    [ TDATA(LEFT), TDATA(LEFT+1) ], or that this is the "nearest"
c    interval to TVAL.
c    For TVAL .lt. TDATA(1), use LEFT = 1.
c    For TDATA(NDATA) .lt. TVAL, use LEFT = NDATA - 1.
c
c    Input, integer N, the order of the basis matrix.
c
c    Input, double precision MBASIS(N,N), the basis matrix.
c
c    Input, integer NDATA, the dimension of the vectors TDATA
c    and YDATA.
c
c    Input, double precision TDATA(NDATA), the abscissa values.  This routine
c    assumes that the TDATA values are uniformly spaced, with an
c    increment of 1.0.
c
c    Input, double precision YDATA(NDATA), the data values to be
c    interpolated or approximated.
c
c    Input, double precision TVAL, the value of T at which the spline is to be
c    evaluated.
c
c    Output, double precision YVAL, the value of the spline at TVAL.
c
      implicit none

      integer maxn 
      parameter ( maxn = 4 )
      integer n
      integer ndata

      double precision arg
      double precision d
      integer first
      integer i
      integer j
      integer left
      double precision mbasis(n,n)
      double precision tdata(ndata)
      double precision tval
      double precision tvec(maxn)
      double precision ydata(ndata)
      double precision yval

      if ( left .eq. 1 ) then
        arg = 0.5D+00 * ( tval - tdata(left) )
        first = left
      else if ( left .lt. ndata - 1 ) then
        arg = tval - tdata(left)
        first = left - 1
      else if ( left .eq. ndata - 1 ) then
        arg = 0.5D+00 * ( 1.0D+00 + tval - tdata(left) )
        first = left - 1
      end if
c
c  TVEC(I) = ARG^(N-I).
c
      tvec(n) = 1.0D+00
      do i = n-1, 1, -1
        tvec(i) = arg * tvec(i+1)
      end do

      yval = 0.0D+00
      do j = 1, n
        d = 0.0D+00
        do i = 1, n
          d = d + tvec(i) * mbasis(i,j)
        end do
        yval = yval + d * ydata(first - 1 + j)
      end do

      return
      end
      subroutine bc_val ( n, t, xcon, ycon, xval, yval )

c*********************************************************************72
c
cc BC_VAL evaluates a parameterized N-th degree Bezier curve in 2D.
c
c  Discussion:
c
c    BC_VAL(T) is the value of a vector function of the form
c
c      BC_VAL(T) = ( X(T), Y(T) )
c
c    where
c
c      X(T) = sum ( 0 .le. I .le. N ) XCON(I) * BERN(I,N)(T)
c      Y(T) = sum ( 0 .le. I .le. N ) YCON(I) * BERN(I,N)(T)
c
c    BERN(I,N)(T) is the I-th Bernstein polynomial of order N
c    defined on the interval [0,1],
c
c    XCON(0:N) and YCON(0:N) are the coordinates of N+1 "control points".
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    12 February 2004
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    David Kahaner, Cleve Moler, Steven Nash,
c    Numerical Methods and Software,
c    Prentice Hall, 1989,
c    ISBN: 0-13-627258-4,
c    LC: TA345.K34.
c
c  Parameters:
c
c    Input, integer N, the degree of the Bezier curve.
c    N must be at least 0.
c
c    Input, double precision T, the point at which the Bezier curve should
c    be evaluated.  The best results are obtained within the interval
c    [0,1] but T may be anywhere.
c
c    Input, double precision XCON(0:N), YCON(0:N), the X and Y coordinates
c    of the control points.  The Bezier curve will pass through
c    the points ( XCON(0), YCON(0) ) and ( XCON(N), YCON(N) ), but
c    generally NOT through the other control points.
c
c    Output, double precision XVAL, YVAL, the X and Y coordinates of the point
c    on the Bezier curve corresponding to the given T value.
c
      implicit none

      integer n

      double precision bval(0:n)
      double precision r8vec_dot_product
      double precision t
      double precision xcon(0:n)
      double precision xval
      double precision ycon(0:n)
      double precision yval

      call bp01 ( n, t, bval )

      xval = r8vec_dot_product ( n + 1, xcon(0:n), bval(0:n) )
      yval = r8vec_dot_product ( n + 1, ycon(0:n), bval(0:n) )

      return
      end
      function bez_val ( n, x, a, b, y )

c*********************************************************************72
c
cc BEZ_VAL evaluates an N-th degree Bezier function at a point.
c
c  Discussion:
c
c    The Bezier function has the form:
c
c      BEZ(X) = sum ( 0 .le. I .le. N ) Y(I) * BERN(N,I)( (X-A)/(B-A) )
c
c    BERN(N,I)(X) is the I-th Bernstein polynomial of order N
c    defined on the interval [0,1],
c
c    Y(0:N) is a set of coefficients,
c
c    and if, for I = 0 to N, we define the N+1 points
c
c      X(I) = ( (N-I)*A + I*B) / N,
c
c    equally spaced in [A,B], the pairs ( X(I), Y(I) ) can be regarded as
c    "control points".
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    12 February 2004
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    David Kahaner, Cleve Moler, Steven Nash,
c    Numerical Methods and Software,
c    Prentice Hall, 1989,
c    ISBN: 0-13-627258-4,
c    LC: TA345.K34.
c
c  Parameters:
c
c    Input, integer N, the degree of the Bezier function.
c    N must be at least 0.
c
c    Input, double precision X, the point at which the Bezier function should
c    be evaluated.  The best results are obtained within the interval
c    [A,B] but X may be anywhere.
c
c    Input, double precision A, B, the interval over which the Bezier function
c    has been defined.  This is the interval in which the control
c    points have been set up.  Note BEZ(A) = Y(0) and BEZ(B) = Y(N),
c    although BEZ will not, in general pass through the other
c    control points.  A and B must not be equal.
c
c    Input, double precision Y(0:N), a set of data defining the Y coordinates
c    of the control points.
c
c    Output, double precision BEZ_VAL, the value of the Bezier function at X.
c
      implicit none

      integer n

      double precision a
      double precision b
      double precision bez_val
      double precision bval(0:n)
      double precision r8vec_dot_product
      double precision x
      double precision x01
      double precision y(0:n)

      if ( b - a .eq. 0.0D+00 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'BEZ_VAL - Fatal error!'
        write ( *, '(a,g14.6)' ) '  Null interval, A = B = ', a
        stop
      end if
c
c  X01 lies in [0,1], in the same relative position as X in [A,B].
c
      x01 = ( x - a ) / ( b - a )

      call bp01 ( n, x01, bval )

      bez_val = r8vec_dot_product ( n + 1, y(0:n), bval(0:n) )

      return
      end
      subroutine bp01 ( n, x, bern )

c*********************************************************************72
c
cc BP01 evaluates the Bernstein basis polynomials for [0,1] at a point.
c
c  Discussion:
c
c    For any N greater than or equal to 0, there is a set of N+1 Bernstein
c    basis polynomials, each of degree N, which form a basis for
c    all polynomials of degree N on [0,1].
c
c    BERN(N,I,X) = [N!/(I!*(N-I)!)] * (1-X)**(N-I) * X**I
c
c    N is the degree;
c
c    0 .le. I .le. N indicates which of the N+1 basis polynomials
c    of degree N to choose;
c
c    X is a point in [0,1] at which to evaluate the basis polynomial.
c
c  First values:
c
c    B(0,0,X) = 1
c
c    B(1,0,X) =      1-X
c    B(1,1,X) =                X
c
c    B(2,0,X) =     (1-X)**2
c    B(2,1,X) = 2 * (1-X)    * X
c    B(2,2,X) =                X**2
c
c    B(3,0,X) =     (1-X)**3
c    B(3,1,X) = 3 * (1-X)**2 * X
c    B(3,2,X) = 3 * (1-X)    * X**2
c    B(3,3,X) =                X**3
c
c    B(4,0,X) =     (1-X)**4
c    B(4,1,X) = 4 * (1-X)**3 * X
c    B(4,2,X) = 6 * (1-X)**2 * X**2
c    B(4,3,X) = 4 * (1-X)    * X**3
c    B(4,4,X) =                X**4
c
c  Special values:
c
c    B(N,I,1/2) = C(N,K) / 2**N
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    12 February 2004
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    David Kahaner, Cleve Moler, Steven Nash,
c    Numerical Methods and Software,
c    Prentice Hall, 1989,
c    ISBN: 0-13-627258-4,
c    LC: TA345.K34.
c
c  Parameters:
c
c    Input, integer N, the degree of the Bernstein basis 
c    polynomials.  N must be at least 0.
c
c    Input, double precision X, the evaluation point.
c
c    Output, double precision BERN(0:N), the values of the N+1 Bernstein basis
c    polynomials at X.
c
      implicit none

      integer n

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
          do j = i-1, 1, -1
            bern(j) = x * bern(j-1) + ( 1.0D+00 - x ) * bern(j)
          end do
          bern(0) = ( 1.0D+00 - x ) * bern(0)
        end do

      end if

      return
      end
      subroutine bpab ( n, a, b, x, bern )

c*********************************************************************72
c
cc BPAB evaluates the Bernstein basis polynomials for [A,B] at a point.
c
c  Discussion:
c
c    BERN(N,I,X) = [N!/(I!*(N-I)!)] * (B-X)**(N-I) * (X-A)**I / (B-A)**N
c
c    B(0,0,X) =   1
c
c    B(1,0,X) = (      B-X                ) / (B-A)
c    B(1,1,X) = (                 X-A     ) / (B-A)
c
c    B(2,0,X) = (     (B-X)**2            ) / (B-A)**2
c    B(2,1,X) = ( 2 * (B-X)    * (X-A)    ) / (B-A)**2
c    B(2,2,X) = (                (X-A)**2 ) / (B-A)**2
c
c    B(3,0,X) = (     (B-X)**3            ) / (B-A)**3
c    B(3,1,X) = ( 3 * (B-X)**2 * (X-A)    ) / (B-A)**3
c    B(3,2,X) = ( 3 * (B-X)    * (X-A)**2 ) / (B-A)**3
c    B(3,3,X) = (                (X-A)**3 ) / (B-A)**3
c
c    B(4,0,X) = (     (B-X)**4            ) / (B-A)**4
c    B(4,1,X) = ( 4 * (B-X)**3 * (X-A)    ) / (B-A)**4
c    B(4,2,X) = ( 6 * (B-X)**2 * (X-A)**2 ) / (B-A)**4
c    B(4,3,X) = ( 4 * (B-X)    * (X-A)**3 ) / (B-A)**4
c    B(4,4,X) = (                (X-A)**4 ) / (B-A)**4
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    12 February 2004
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    David Kahaner, Cleve Moler, Steven Nash,
c    Numerical Methods and Software,
c    Prentice Hall, 1989,
c    ISBN: 0-13-627258-4,
c    LC: TA345.K34.
c
c  Parameters:
c
c    Input, integer N, the degree of the Bernstein basis 
c    polynomials.  There is a set of N+1 Bernstein basis polynomials, each 
c    of degree N, which form a basis for polynomials of degree N on [A,B].  
c    N must be at least 0.
c
c    Input, double precision A, B, the endpoints of the interval on which the
c    polynomials are to be based.  A and B should not be equal.
c    Input, double precision X, the point at which the polynomials are to be
c    evaluated.  X need not lie in the interval [A,B].
c
c    Output, double precision BERN(0:N), the values of the N+1 Bernstein basis
c    polynomials at X.
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
          do j = i-1, 1, -1
            bern(j) = ( ( b - x ) * bern(j) 
     &        + ( x - a ) * bern(j-1) ) / ( b - a )
          end do
          bern(0) = ( b - x ) * bern(0) / ( b - a )
        end do

      end if

      return
      end
      subroutine bpab_approx ( n, a, b, ydata, xval, yval )

c*********************************************************************72
c
cc BPAB_APPROX evaluates the Bernstein polynomial approximant to F(X) on [A,B].
c
c  Formula:
c
c    BERN(F)(X) = sum ( 0 .le. I .le. N ) F(X(I)) * B_BASE(I,X)
c
c    where
c
c      X(I) = ( ( N - I ) * A + I * B ) / N
c      B_BASE(I,X) is the value of the I-th Bernstein basis polynomial at X.
c
c  Discussion:
c
c    The Bernstein polynomial BERN(F) for F(X) is an approximant, not an
c    interpolant; in other words, its value is not guaranteed to equal
c    that of F at any particular point.  However, for a fixed interval
c    [A,B], if we let N increase, the Bernstein polynomial converges
c    uniformly to F everywhere in [A,B], provided only that F is continuous.
c    Even if F is not continuous, but is bounded, the polynomial converges
c    pointwise to F(X) at all points of continuity.  On the other hand,
c    the convergence is quite slow compared to other interpolation
c    and approximation schemes.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    12 February 2004
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    David Kahaner, Cleve Moler, Steven Nash,
c    Numerical Methods and Software,
c    Prentice Hall, 1989,
c    ISBN: 0-13-627258-4,
c    LC: TA345.K34.
c
c  Parameters:
c
c    Input, integer N, the degree of the Bernstein polynomial
c    to be used.  N must be at least 0.
c
c    Input, double precision A, B, the endpoints of the interval on which the
c    approximant is based.  A and B should not be equal.
c
c    Input, double precision YDATA(0:N), the data values at N+1 equally
c    spaced points in [A,B].  If N = 0, then the evaluation point should
c    be 0.5 * ( A + B).  Otherwise, evaluation point I should be
c    ( (N-I)*A + I*B ) / N ).
c
c    Input, double precision XVAL, the point at which the Bernstein polynomial
c    approximant is to be evaluated.  XVAL does not have to lie in the
c    interval [A,B].
c
c    Output, double precision YVAL, the value of the Bernstein polynomial
c    approximant for F, based in [A,B], evaluated at XVAL.
c
      implicit none

      integer n

      double precision a
      double precision b
      double precision bvec(0:n)
      double precision r8vec_dot_product
      double precision xval
      double precision ydata(0:n)
      double precision yval
c
c  Evaluate the Bernstein basis polynomials at XVAL.
c
      call bpab ( n, a, b, xval, bvec )
c
c  Now compute the sum of YDATA(I) * BVEC(I).
c
      yval = r8vec_dot_product ( n + 1, ydata(0:n), bvec(0:n) )

      return
      end
      subroutine chfev ( x1, x2, f1, f2, d1, d2, ne, xe, fe, next, 
     &  ierr )

c*********************************************************************72
c
cc CHFEV evaluates a cubic polynomial given in Hermite form.
c
c  Discussion:
c
c    This routine evaluates a cubic polynomial given in Hermite form at an
c    array of points.  While designed for use by SPLINE_PCHIP_VAL, it may
c    be useful directly as an evaluator for a piecewise cubic
c    Hermite function in applications, such as graphing, where
c    the interval is known in advance.
c
c    The cubic polynomial is determined by function values
c    F1, F2 and derivatives D1, D2 on the interval [X1,X2].
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    25 June 2008
c
c  Author:
c
c    Original FORTRAN77 version by Fred Fritsch.
c    FORTRAN90 version by John Burkardt.
c
c  Reference:
c
c    Fred Fritsch, Ralph Carlson,
c    Monotone Piecewise Cubic Interpolation,
c    SIAM Journal on Numerical Analysis,
c    Volume 17, Number 2, April 1980, pages 238-246.
c
c    David Kahaner, Cleve Moler, Steven Nash,
c    Numerical Methods and Software,
c    Prentice Hall, 1989,
c    ISBN: 0-13-627258-4,
c    LC: TA345.K34.
c
c  Parameters:
c
c    Input, double precision X1, X2, the endpoints of the interval of
c    definition of the cubic.  X1 and X2 must be distinct.
c
c    Input, double precision F1, F2, the values of the function at X1 and
c    X2, respectively.
c
c    Input, double precision D1, D2, the derivative values at X1 and
c    X2, respectively.
c
c    Input, integer NE, the number of evaluation points.
c
c    Input, double precision XE(NE), the points at which the function is to
c    be evaluated.  If any of the XE are outside the interval
c    [X1,X2], a warning error is returned in NEXT.
c
c    Output, double precision FE(NE), the value of the cubic function
c    at the points XE.
c
c    Output, integer NEXT(2), indicates the number of
c    extrapolation points:
c    NEXT(1) = number of evaluation points to the left of interval.
c    NEXT(2) = number of evaluation points to the right of interval.
c
c    Output, integer IERR, error flag.
c    0, no errors.
c    -1, NE .lt. 1.
c    -2, X1 .eq. X2.
c
      implicit none

      integer ne

      double precision c2
      double precision c3
      double precision d1
      double precision d2
      double precision del1
      double precision del2
      double precision delta
      double precision f1
      double precision f2
      double precision fe(ne)
      double precision h
      integer i
      integer ierr
      integer next(2)
      double precision x
      double precision x1
      double precision x2
      double precision xe(ne)
      double precision xma
      double precision xmi

      if ( ne .lt. 1 ) then
        ierr = -1
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'CHFEV - Fatal error!'
        write ( *, '(a)' ) 
     &    '  Number of evaluation points is less than 1.'
        write ( *, '(a,i8)' ) '  NE = ', ne
        stop
      end if

      h = x2 - x1

      if ( h .eq. 0.0D+00 ) then
        ierr = -2
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'CHFEV - Fatal error!'
        write ( *, '(a)' ) '  The interval [X1,X2] is of zero length.'
        stop
      end if
c
c  Initialize.
c
      ierr = 0
      next(1) = 0
      next(2) = 0
      xmi = min ( 0.0D+00, h )
      xma = max ( 0.0D+00, h )
c
c  Compute cubic coefficients expanded about X1.
c
      delta = ( f2 - f1 ) / h
      del1 = ( d1 - delta ) / h
      del2 = ( d2 - delta ) / h
      c2 = -( del1 + del1 + del2 )
      c3 = ( del1 + del2 ) / h
c
c  Evaluation loop.
c
      do i = 1, ne

        x = xe(i) - x1
        fe(i) = f1 + x * ( d1 + x * ( c2 + x * c3 ) )
c
c  Count the extrapolation points.
c
        if ( x .lt. xmi ) then
          next(1) = next(1) + 1
        end if

        if ( xma .lt. x ) then
          next(2) = next(2) + 1
        end if

      end do

      return
      end
      subroutine data_to_dif ( ntab, xtab, ytab, diftab )

c*********************************************************************72
c
cc DATA_TO_DIF sets up a divided difference table from raw data.
c
c  Discussion:
c
c    Space can be saved by using a single array for both the DIFTAB and
c    YTAB dummy parameters.  In that case, the divided difference table will
c    overwrite the Y data without interfering with the computation.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    04 September 2004
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Carl deBoor,
c    A Practical Guide to Splines,
c    Springer, 2001,
c    ISBN: 0387953663.
c
c  Parameters:
c
c    Input, integer NTAB, the number of pairs of points
c    (XTAB(I),YTAB(I)) which are to be used as data.  The
c    number of entries to be used in DIFTAB, XTAB and YTAB.
c
c    Input, double precision XTAB(NTAB), the X values at which data was taken.
c    These values must be distinct.
c
c    Input, double precision YTAB(NTAB), the corresponding Y values.
c
c    Output, double precision DIFTAB(NTAB), the divided difference coefficients
c    corresponding to the input (XTAB,YTAB).
c
      implicit none

      integer ntab

      double precision diftab(ntab)
      integer i
      integer j
      logical r8vec_distinct
      double precision xtab(ntab)
      double precision ytab(ntab)

      if ( .not. r8vec_distinct ( ntab, xtab ) ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'DATA_TO_DIF - Fatal error!'
        write ( *, '(a)' ) '  Two entries of XTAB are equal!'
        return
      end if
c
c  Copy the data values into DIFTAB.
c
      do i = 1, ntab
        diftab(i) = ytab(i)
      end do
c
c  Compute the divided differences.
c
      do i = 2, ntab
        do j = ntab, i, -1

          diftab(j) = ( diftab(j) - diftab(j-1) ) 
     &      / ( xtab(j) - xtab(j+1-i) )

        end do
      end do

      return
      end
      subroutine dif_val ( ntab, xtab, diftab, xval, yval )

c*********************************************************************72
c
cc DIF_VAL evaluates a divided difference polynomial at a point.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    05 September 2004
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Carl deBoor,
c    A Practical Guide to Splines,
c    Springer, 2001,
c    ISBN: 0387953663.
c
c  Parameters:
c
c    Input, integer NTAB, the number of divided difference
c    coefficients in DIFTAB, and the number of points XTAB.
c
c    Input, double precision XTAB(NTAB), the X values upon which the
c    divided difference polynomial is based.
c
c    Input, double precision DIFTAB(NTAB), the divided difference
c    polynomial coefficients.
c
c    Input, double precision XVAL, the value where the polynomial
c    is to be evaluated.
c
c    Output, double precision YVAL, the value of the polynomial at XVAL.
c
      implicit none

      integer ntab

      double precision diftab(ntab)
      integer i
      double precision xtab(ntab)
      double precision xval
      double precision yval

      yval = diftab(ntab)
      do i = 1, ntab-1
        yval = diftab(ntab-i) + ( xval - xtab(ntab-i) ) * yval
      end do

      return
      end
      subroutine least_set_old ( ntab, xtab, ytab, ndeg, ptab, b, c, 
     &  d, eps, ierror )

c*********************************************************************72
c
cc LEAST_SET_OLD constructs the least squares polynomial approximation to data.
c
c  Discussion:
c
c    The least squares polynomial is not returned directly as a simple
c    polynomial.  Instead, it is represented in terms of a set of
c    orthogonal polynomials appopriate for the given data.  This makes
c    the computation more accurate, but means that the user can not
c    easily evaluate the computed polynomial.  Instead, the routine
c    LEAST_EVAL should be used to evaluate the least squares polynomial
c    at any point.  (However, the value of the least squares polynomial
c    at each of the data points is returned as part of this computation.)
c
c
c    A discrete unweighted inner product is used, so that
c
c      ( F(X), G(X) ) = sum ( 1 .le. I .le. NTAB ) F(XTAB(I)) * G(XTAB(I)).
c
c    The least squares polynomial is determined using a set of
c    orthogonal polynomials PHI.  These polynomials can be defined
c    recursively by:
c
c      PHI(0)(X) = 1
c      PHI(1)(X) = X - B(1)
c      PHI(I)(X) = ( X - B(I) ) * PHI(I-1)(X) - D(I) * PHI(I-2)(X)
c
c    The array B(1:NDEG) contains the values
c
c      B(I) = ( X*PHI(I-1), PHI(I-1) ) / ( PHI(I-1), PHI(I-1) )
c
c    The array D(2:NDEG) contains the values
c
c      D(I) = ( PHI(I-1), PHI(I-1) ) / ( PHI(I-2), PHI(I-2) )
c
c    Using this basis, the least squares polynomial can be represented as
c
c      P(X)(I) = sum ( 0 .le. I .le. NDEG ) C(I) * PHI(I)(X)
c
c    The array C(0:NDEG) contains the values
c
c      C(I) = ( YTAB(I), PHI(I) ) / ( PHI(I), PHI(I) )
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    15 May 2004
c
c  Author:
c
c    FORTRAN90 version by John Burkardt
c
c  Reference:
c
c    Gisela Engeln-Muellges, Frank Uhlig,
c    Numerical Algorithms with C,
c    Springer, 1996,
c    ISBN: 3-540-60530-4.
c
c  Parameters:
c
c    Input, integer NTAB, the number of data points.
c
c    Input, double precision XTAB(NTAB), the X data.  The values in XTAB
c    should be distinct, and in increasing order.
c
c    Input, double precision YTAB(NTAB), the Y data values corresponding
c    to the X data in XTAB.
c
c    Input, integer NDEG, the degree of the polynomial which the
c    program is to use.  NDEG must be at least 0, and less than or
c    equal to NTAB-1.
c
c    Output, double precision PTAB(NTAB), the value of the least
c    squares polynomial at the points XTAB(1:NTAB).
c
c    Output, double precision B(1:NDEG), C(0:NDEG), D(2:NDEG), arrays
c    needed to evaluate the polynomial.
c
c    Output, double precision EPS, the root-mean-square discrepancy of the
c    polynomial fit.
c
c    Output, integer IERROR, error flag.
c    zero, no error occurred;
c    nonzero, an error occurred, and the polynomial could not be computed.
c
      implicit none

      integer ndeg
      integer ntab

      double precision b(1:ndeg)
      double precision c(0:ndeg)
      double precision d(2:ndeg)
      double precision eps
      integer i
      integer i0l1
      integer i1l1
      integer ierror
      integer it
      integer k
      integer mdeg
      double precision ptab(ntab)
      double precision rn0
      double precision rn1
      double precision s
      double precision sum2
      double precision xtab(ntab)
      double precision y_sum
      double precision ytab(ntab)
      double precision ztab(2*ntab)

      ierror = 0
c
c  Check NDEG.
c
      if ( ndeg .lt. 0 ) then
        ierror = 1
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'LEAST_SET_OLD - Fatal error!'
        write ( *, '(a)' ) '  NDEG .lt. 0.'
        stop
      end if

      if ( ntab .le. ndeg ) then
        ierror = 1
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'LEAST_SET_OLD - Fatal error!'
        write ( *, '(a)' ) '  NTAB .le. NDEG.'
        stop
      end if
c
c  Check that the abscissas are strictly increasing.
c
      do i = 1, ntab-1
        if ( xtab(i+1) .le. xtab(i) ) then
          ierror = 1
          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) 'LEAST_SET_OLD - Fatal error!'
          write ( *, '(a)' ) '  XTAB must be strictly increasing, but'
          write ( *, '(a,i8,a,g14.6)' ) '  XTAB(', i, ') = ', xtab(i)
          write ( *, '(a,i8,a,g14.6)' ) 
     &      '  XTAB(', i+1, ') = ', xtab(i+1)
          stop
        end if
      end do

      i0l1 = 0
      i1l1 = ntab
c
c  The polynomial is of degree at least 0.
c
      y_sum = 0.0D+00
      do i = 1, ntab
        y_sum = y_sum + ytab(i)
      end do
      rn0 = ntab
      c(0) = y_sum / dble ( ntab )

      do i = 1, ntab
        ptab(i) = y_sum / dble ( ntab )
      end do

      if ( ndeg .eq. 0 ) then
        eps = 0.0D+00
        do i = 1, ntab
          eps = eps + ( ptab(i) - ytab(i) )**2
        end do
        eps = sqrt ( eps / dble ( ntab ) )
        return
      end if
c
c  The polynomial is of degree at least 1.
c
      b(1) = 0.0D+00
      do i = 1, ntab
        b(1) = b(1) + xtab(i)
      end do
      b(1) = b(1) / dble ( ntab )

      s = 0.0D+00
      sum2 = 0.0D+00
      do i = 1, ntab
        ztab(i1l1+i) = xtab(i) - b(1)
        s = s + ztab(i1l1+i)**2
        sum2 = sum2 + ztab(i1l1+i) * ( ytab(i) - ptab(i) )
      end do

      rn1 = s
      c(1) = sum2 / s

      do i = 1, ntab
        ptab(i) = ptab(i) + c(1) * ztab(i1l1+i)
      end do

      if ( ndeg .eq. 1 ) then
        eps = 0.0D+00
        do i = 1, ntab
          eps = eps + ( ptab(i) - ytab(i) )**2
        end do
        eps = sqrt ( eps / dble ( ntab ) )
        return
      end if

      do i = 1, ntab
        ztab(i) = 1.0D+00
      end do

      mdeg = 2
      k = 2

      do

        d(k) = rn1 / rn0

        sum2 = 0.0D+00
        do i = 1, ntab
          sum2 = sum2 + xtab(i) * ztab(i1l1+i) * ztab(i1l1+i)
        end do

        b(k) = sum2 / rn1

        s = 0.0D+00
        sum2 = 0.0D+00
        do i = 1, ntab
          ztab(i0l1+i) = ( xtab(i) - b(k) ) * ztab(i1l1+i) 
     &      - d(k) * ztab(i0l1+i)
          s = s + ztab(i0l1+i) * ztab(i0l1+i)
          sum2 = sum2 + ztab(i0l1+i) * ( ytab(i) - ptab(i) )
        end do

        rn0 = rn1
        rn1 = s

        c(k) = sum2 / rn1

        it = i0l1
        i0l1 = i1l1
        i1l1 = it

        do i = 1, ntab
          ptab(i) = ptab(i) + c(k) * ztab(i1l1+i)
        end do

        if ( ndeg .le. mdeg ) then
          exit
        end if

        mdeg = mdeg + 1
        k = k + 1

      end do
c
c  Compute the RMS error.
c
      eps = 0.0D+00
      do i = 1, ntab
        eps = eps + ( ptab(i) - ytab(i) )**2
      end do
      eps = sqrt ( eps / dble ( ntab ) )

      return
      end
      subroutine least_val_old ( x, ndeg, b, c, d, value )

c*********************************************************************72
c
cc LEAST_VAL_OLD evaluates a least squares polynomial defined by LEAST_SET_OLD.
c
c  Discussion:
c
c    This is an "old" version of the routine.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    15 May 2004
c
c  Author:
c
c    FORTRAN90 version by John Burkardt
c
c  Reference:
c
c    Gisela Engeln-Muellges, Frank Uhlig,
c    Numerical Algorithms with C,
c    Springer, 1996,
c    ISBN: 3-540-60530-4.
c
c  Parameters:
c
c    Input, double precision X, the point at which the polynomial is
c    to be evaluated.
c
c    Input, integer NDEG, the degree of the least squares
c    polynomial.
c
c    Input, double precision B(1:NDEG), C(0:NDEG), D(2:NDEG), arrays
c    defined by LEAST_SET_OLD, and needed to evaluate the polynomial.
c
c    Output, double precision VALUE, the value of the polynomial at X.
c
      implicit none

      integer ndeg

      double precision b(1:ndeg)
      double precision c(0:ndeg)
      double precision d(2:ndeg)
      integer k
      double precision sk
      double precision skp1
      double precision skp2
      double precision value
      double precision x

      if ( ndeg .le. 0 ) then

        value = c(0)

      else if ( ndeg .eq. 1 ) then

        value = c(0) + c(1) * ( x - b(1) )

      else

        skp2 = c(ndeg)
        skp1 = c(ndeg-1) + ( x - b(ndeg) ) * skp2

        do k = ndeg-2, 0, -1
          sk = c(k) + ( x - b(k+1) ) * skp1 - d(k+2) * skp2
          skp2 = skp1
          skp1 = sk
        end do

        value = sk

      end if

      return
      end
      subroutine least_set ( point_num, x, f, w, nterms, b, c, d )

c*********************************************************************72
c
cc LEAST_SET defines a least squares polynomial for given data.
c
c  Discussion:
c
c    This routine is based on ORTPOL by Conte and deBoor.
c
c    The polynomial may be evaluated at any point X by calling LEAST_VAL.
c
c    Thanks to Andrew Telford for pointing out a mistake in the form of
c    the check that there are enough unique data points, 25 June 2008.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    25 June 2008
c
c  Author:
c
c    Original FORTRAN77 version by Samuel Conte, Carl deBoor.
c    FORTRAN90 version by John Burkardt
c
c  Reference:
c
c    Samuel Conte, Carl deBoor,
c    Elementary Numerical Analysis,
c    Second Edition,
c    McGraw Hill, 1972,
c    ISBN: 07-012446-4,
c    LC: QA297.C65.
c
c  Parameters:
c
c    Input, integer POINT_NUM, the number of data values.
c
c    Input, double precision X(POINT_NUM), the abscissas of the data points.
c    At least NTERMS of the values in X must be distinct.
c
c    Input, double precision F(POINT_NUM), the data values at the points X(*).
c
c    Input, double precision W(POINT_NUM), the weights associated with
c    the data points.  Each entry of W should be positive.
c
c    Input, integer NTERMS, the number of terms to use in the
c    approximating polynomial.  NTERMS must be at least 1.
c    The degree of the polynomial is NTERMS-1.
c
c    Output, double precision B(NTERMS), C(NTERMS), D(NTERMS), are quantities
c    defining the least squares polynomial for the input data,
c    which will be needed to evaluate the polynomial.
c
      implicit none

      integer point_num
      integer nterms

      double precision b(nterms)
      double precision bsum
      double precision c(nterms)
      double precision d(nterms)
      double precision dsum
      double precision f(point_num)
      integer i
      integer j
      double precision p
      double precision pj(point_num)
      double precision pjm1(point_num)
      double precision s(nterms)
      double precision ssum
      double precision tol
      parameter ( tol = 0.0D+00 )
      integer unique_num
      double precision w(point_num)
      double precision x(point_num)
c
c  Make sure at least NTERMS X values are unique.
c
      call r8vec_unique_count ( point_num, x, tol, unique_num )

      if ( unique_num .lt. nterms ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'LEAST_SET - Fatal error!'
        write ( *, '(a)' ) '  The number of distinct X values must be'
        write ( *, '(a,i8)') '  at least NTERMS = ', nterms
        write ( *, '(a,i8)' ) 
     &    '  but the input data has only ', unique_num
        write ( *, '(a)' ) '  distinct entries.'
        return
      end if
c
c  Make sure all W entries are positive.
c
      do i = 1, point_num
        if ( w(i) .le. 0.0D+00 ) then
          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) 'LEAST_SET - Fatal error!'
          write ( *, '(a)' ) '  All weights W must be positive,'
          write ( *, '(a,i8)' ) '  but weight ', i
          write ( *, '(a,g14.6)' ) '  is ', w(i)
          return
        end if
      end do
c
c  Start inner product summations at zero.
c
      do i = 1, nterms
        b(i) = 0.0D+00
        c(i) = 0.0D+00
        d(i) = 0.0D+00
        s(i) = 0.0D+00
      end do
c
c  Set the values of P(-1,X) and P(0,X) at all data points.
c
      do i = 1, point_num
        pjm1(i) = 0.0D+00
        pj(i) = 1.0D+00
      end do
c
c  Now compute the value of P(J,X(I)) as
c
c    P(J,X(I)) = ( X(I) - B(J) ) * P(J-1,X(I)) - C(J) * P(J-2,X(I))
c
c  where
c
c    S(J) = .lt. P(J,X), P(J,X) >
c    B(J) = .lt. x*P(J,X), P(J,X) > / .lt. P(J,X), P(J,X) >
c    C(J) = S(J) / S(J-1)
c
c  and the least squares coefficients are
c
c    D(J) = .lt. F(X), P(J,X) > / .lt. P(J,X), P(J,X) >
c
      do j = 1, nterms

        dsum = 0.0D+00
        bsum = 0.0D+00
        ssum = 0.0D+00
        do i = 1, point_num
          dsum = dsum + w(i) * f(i) * pj(i)
          bsum = bsum + w(i) * x(i) * pj(i)**2
          ssum = ssum + w(i) * pj(i)**2
        end do

        d(j) = d(j) + dsum
        b(j) = b(j) + bsum
        s(j) = s(j) + ssum

        d(j) = d(j) / s(j)

        if ( j .eq. nterms ) then
          c(j) = 0.0D+00
          return
        end if

        b(j) = b(j) / s(j)

        if ( j .eq. 1 ) then
          c(j) = 0.0D+00
        else
          c(j) = s(j) / s(j-1)
        end if

        do i = 1, point_num
          p = pj(i)
          pj(i) = ( x(i) - b(j) ) * pj(i) - c(j) * pjm1(i)
          pjm1(i) = p
        end do

      end do

      return
      end
      subroutine least_val ( nterms, b, c, d, x, px )

c*********************************************************************72
c
cc LEAST_VAL evaluates a least squares polynomial defined by LEAST_SET.
c
c  Discussion:
c
c    The least squares polynomial is assumed to be defined as a sum
c
c      P(X) = sum ( 1 .le. I .le. NTERMS ) D(I) * P(I-1,X)
c
c    where the orthogonal basis polynomials P(I,X) satisfy the following
c    three term recurrence:
c
c      P(-1,X) = 0
c      P(0,X) = 1
c      P(I,X) = ( X - B(I-1) ) * P(I-1,X) - C(I-1) * P(I-2,X)
c
c    Therefore, the least squares polynomial can be evaluated as follows:
c
c    If NTERMS is 1, then the value of P(X) is D(1) * P(0,X) = D(1).
c
c    Otherwise, P(X) is defined as the sum of NTERMS > 1 terms.  We can
c    reduce the number of terms by 1, because the polynomial P(NTERMS,X)
c    can be rewritten as a sum of polynomials;  Therefore, P(NTERMS,X)
c    can be eliminated from the sum, and its coefficient merged in with
c    those of other polynomials.  Repeat this process for P(NTERMS-1,X)
c    and so on until a single term remains.
c    P(NTERMS,X) of P(NTERMS-1,X)
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    24 May 2005
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Samuel Conte, Carl deBoor,
c    Elementary Numerical Analysis,
c    Second Edition,
c    McGraw Hill, 1972,
c    ISBN: 07-012446-4,
c    LC: QA297.C65.
c
c  Parameters:
c
c    Input, integer NTERMS, the number of terms in the least
c    squares polynomial.  NTERMS must be at least 1.  The input value of NTERMS
c    may be reduced from the value given to LEAST_SET.  This will
c    evaluate the least squares polynomial of the lower degree specified.
c
c    Input, double precision B(NTERMS), C(NTERMS), D(NTERMS), the information
c    computed by LEAST_SET.
c
c    Input, double precision X, the point at which the least squares polynomial
c    is to be evaluated.
c
c    Output, double precision PX, the value of the least squares
c    polynomial at X.
c
      implicit none

      integer nterms

      double precision b(nterms)
      double precision c(nterms)
      double precision d(nterms)
      integer i
      double precision prev
      double precision prev2
      double precision px
      double precision x

      px = d(nterms)
      prev = 0.0D+00

      do i = nterms - 1, 1, -1

        prev2 = prev
        prev = px

        if ( i .eq. nterms-1 ) then
          px = d(i) + ( x - b(i) ) * prev
        else
          px = d(i) + ( x - b(i) ) * prev - c(i+1) * prev2
        end if

      end do

      return
      end
      subroutine least_val2 ( nterms, b, c, d, x, px, pxp )

c*********************************************************************72
c
cc LEAST_VAL2 evaluates a least squares polynomial defined by LEAST_SET.
c
c  Discussion:
c
c    This routine also computes the derivative of the polynomial.
c
c    The least squares polynomial is assumed to be defined as a sum
c
c      P(X) = sum ( 1 .le. I .le. NTERMS ) D(I) * P(I-1,X)
c
c    where the orthogonal basis polynomials P(I,X) satisfy the following
c    three term recurrence:
c
c      P(-1,X) = 0
c      P(0,X) = 1
c      P(I,X) = ( X - B(I-1) ) * P(I-1,X) - C(I-1) * P(I-2,X)
c
c    Therefore, the least squares polynomial can be evaluated as follows:
c
c    If NTERMS is 1, then the value of P(X) is D(1) * P(0,X) = D(1).
c
c    Otherwise, P(X) is defined as the sum of NTERMS > 1 terms.  We can
c    reduce the number of terms by 1, because the polynomial P(NTERMS,X)
c    can be rewritten as a sum of polynomials;  Therefore, P(NTERMS,X)
c    can be eliminated from the sum, and its coefficient merged in with
c    those of other polynomials.  Repeat this process for P(NTERMS-1,X)
c    and so on until a single term remains.
c    P(NTERMS,X) of P(NTERMS-1,X)
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    24 May 2005
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer NTERMS, the number of terms in the least 
c    squares polynomial.  NTERMS must be at least 1.  The value of NTERMS
c    may be reduced from the value given to LEAST_SET.
c    This will cause LEAST_VAL to evaluate the least squares polynomial
c    of the lower degree specified.
c
c    Input, double precision B(NTERMS), C(NTERMS), D(NTERMS), the information
c    computed by LEAST_SET.
c
c    Input, double precision X, the point at which the least squares polynomial
c    is to be evaluated.
c
c    Output, double precision PX, PXP, the value and derivative of the least
c    squares polynomial at X.
c
      implicit none

      integer nterms

      double precision b(nterms)
      double precision c(nterms)
      double precision d(nterms)
      integer i
      double precision px
      double precision pxm1
      double precision pxm2
      double precision pxp
      double precision pxpm1
      double precision pxpm2
      double precision x

      px = d(nterms)
      pxp = 0.0D+00
      pxm1 = 0.0D+00
      pxpm1 = 0.0D+00

      do i = nterms - 1, 1, -1

        pxm2 = pxm1
        pxpm2 = pxpm1
        pxm1 = px
        pxpm1 = pxp

        if ( i .eq. nterms - 1 ) then
          px = d(i) + ( x - b(i) ) * pxm1
          pxp = pxm1 + ( x - b(i) ) * pxpm1
        else
          px = d(i) + ( x - b(i) ) * pxm1 - c(i+1) * pxm2
          pxp = pxm1 + ( x - b(i) ) * pxpm1 - c(i+1) * pxpm2
        end if

      end do

      return
      end
      subroutine parabola_val2 ( dim_num, ndata, tdata, ydata, left, 
     &  tval, yval )

c*********************************************************************72
c
cc PARABOLA_VAL2 evaluates a parabolic interpolant through tabular data.
c
c  Discussion:
c
c    This routine is a utility routine used by OVERHAUSER_SPLINE_VAL.
c    It constructs the parabolic interpolant through the data in
c    3 consecutive entries of a table and evaluates this interpolant
c    at a given abscissa value.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    26 January 2004
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer DIM_NUM, the dimension of a single data point.
c    DIM_NUM must be at least 1.
c
c    Input, integer NDATA, the number of data points.
c    NDATA must be at least 3.
c
c    Input, double precision TDATA(NDATA), the abscissas of the data
c    points.  The values in TDATA must be in strictly ascending order.
c
c    Input, double precision YDATA(DIM_NUM,NDATA), the data points
c    corresponding to the abscissas.
c
c    Input, integer LEFT, the location of the first of the three
c    consecutive data points through which the parabolic interpolant
c    must pass.  1 .le. LEFT .le. NDATA - 2.
c
c    Input, double precision TVAL, the value of T at which the parabolic
c    interpolant is to be evaluated.  Normally, TDATA(1) .le. TVAL .le. T(NDATA),
c    and the data will be interpolated.  For TVAL outside this range,
c    extrapolation will be used.
c
c    Output, double precision YVAL(DIM_NUM), the value of the parabolic
c    interpolant at TVAL.
c
      implicit none

      integer ndata
      integer dim_num

      double precision dif1
      double precision dif2
      integer i
      integer left
      double precision t1
      double precision t2
      double precision t3
      double precision tval
      double precision tdata(ndata)
      double precision ydata(dim_num,ndata)
      double precision y1
      double precision y2
      double precision y3
      double precision yval(dim_num)
c
c  Check.
c
      if ( left .lt. 1 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'PARABOLA_VAL2 - Fatal error!'
        write ( *, '(a)' ) '  LEFT .lt. 1.'
        write ( *, '(a,i8)' ) '  LEFT = ', left
        stop
      end if

      if ( ndata-2 .lt. left ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'PARABOLA_VAL2 - Fatal error!'
        write ( *, '(a)' ) '  NDATA-2 .lt. LEFT.'
        write ( *, '(a,i8)' ) '  NDATA = ', ndata
        write ( *, '(a,i8)' ) '  LEFT =  ', left
        stop
      end if

      if ( dim_num .lt. 1 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'PARABOLA_VAL2 - Fatal error!'
        write ( *, '(a)' ) '  DIM_NUM .lt. 1.'
        stop
      end if
c
c  Copy out the three abscissas.
c
      t1 = tdata(left)
      t2 = tdata(left+1)
      t3 = tdata(left+2)

      if ( t2 .le. t1 .or. t3 .le. t2 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'PARABOLA_VAL2 - Fatal error!'
        write ( *, '(a)' ) '  T2 .le. T1 or T3 .le. T2.'
        stop
      end if
c
c  Construct and evaluate a parabolic interpolant for the data
c  in each dimension.
c
      do i = 1, dim_num

        y1 = ydata(i,left)
        y2 = ydata(i,left+1)
        y3 = ydata(i,left+2)

        dif1 = ( y2 - y1 ) / ( t2 - t1 )
        dif2 = ( ( y3 - y1 ) / ( t3 - t1 ) 
     &       - ( y2 - y1 ) / ( t2 - t1 ) ) / ( t3 - t2 )

        yval(i) = y1 + ( tval - t1 ) * ( dif1 + ( tval - t2 ) * dif2 )

      end do

      return
      end
      function pchst ( arg1, arg2 )

c*********************************************************************72
c
cc PCHST: PCHIP sign-testing routine.
c
c  Discussion:
c
c    This routine essentially computes the sign of ARG1 * ARG2.
c
c    The object is to do this without multiplying ARG1 * ARG2, to avoid
c    possible over/underflow problems.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    25 June 2008
c
c  Author:
c
c    Original FORTRAN77 version by Fred Fritsch.
c    FORTRAN90 version by John Burkardt.
c
c  Reference:
c
c    Fred Fritsch, Ralph Carlson,
c    Monotone Piecewise Cubic Interpolation,
c    SIAM Journal on Numerical Analysis,
c    Volume 17, Number 2, April 1980, pages 238-246.
c
c  Parameters:
c
c    Input, double precision ARG1, ARG2, two values to check.
c
c    Output, double precision PCHST,
c    -1.0, if ARG1 and ARG2 are of opposite sign.
c     0.0, if either argument is zero.
c    +1.0, if ARG1 and ARG2 are of the same sign.
c
      implicit none

      double precision arg1
      double precision arg2
      double precision pchst

      pchst = sign ( 1.0D+00, arg1 ) * sign ( 1.0D+00, arg2 )

      if ( arg1 .eq. 0.0D+00 .or. arg2 .eq. 0.0D+00 ) then
        pchst = 0.0D+00
      end if

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

      double precision r8_uniform_01
      integer k
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
        seed = seed + 2147483647
      end if

      r8_uniform_01 = dble ( seed ) * 4.656612875D-10

      return
      end
      subroutine r83_mxv ( n, a, x, b )

c*********************************************************************72
c
cc R83_MXV multiplies an R83 matrix times a vector.
c
c  Discussion:
c
c    The R83 storage format is used for a tridiagonal matrix.
c    The superdiagonal is stored in entries (1,2:N), the diagonal in
c    entries (2,1:N), and the subdiagonal in (3,1:N-1).  Thus, the
c    original matrix is "collapsed" vertically into the array.
c
c  Example:
c
c    Here is how an R83 matrix of order 5 would be stored:
c
c       *  A12 A23 A34 A45
c      A11 A22 A33 A44 A55
c      A21 A32 A43 A54  *
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    02 November 2003
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the order of the linear system.
c
c    Input, double precision A(3,N), the R83 matrix.
c
c    Input, double precision X(N), the vector to be multiplied by A.
c
c    Output, double precision B(N), the product A * x.
c
      implicit none

      integer n

      double precision a(3,n)
      double precision b(n)
      integer i
      double precision x(n)

      do i = 1, n
        b(i) = a(2,i) * x(i)
      end do

      do i = 1, n - 1
        b(i) = b(i) + a(1,i+1) * x(i+1)
      end do

      do i = 2, n
        b(i) = b(i) + a(3,i-1) * x(i-1)
      end do

      return
      end
      subroutine r83_np_fs ( n, a, b, x )

c*********************************************************************72
c
cc R83_NP_FS factors and solves an R83 system.
c
c  Discussion:
c
c    The R83 storage format is used for a tridiagonal matrix.
c    The superdiagonal is stored in entries (1,2:N), the diagonal in
c    entries (2,1:N), and the subdiagonal in (3,1:N-1).  Thus, the
c    original matrix is "collapsed" vertically into the array.
c
c    This algorithm requires that each diagonal entry be nonzero.
c    It does not use pivoting, and so can fail on systems that
c    are actually nonsingular.
c
c  Example:
c
c    Here is how an R83 matrix of order 5 would be stored:
c
c       *  A12 A23 A34 A45
c      A11 A22 A33 A44 A55
c      A21 A32 A43 A54  *
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    24 May 2009
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the order of the linear system.
c
c    Input/output, double precision A(3,N).
c    On input, the tridiagonal matrix.
c    On output, the data in these vectors has been overwritten
c    by factorization information.
c
c    Input, double precision B(N), the right hand side of the linear system.
c
c    Output, double precision X(N), the solution of the linear system.
c
      implicit none

      integer n

      double precision a(3,n)
      double precision b(n)
      integer i
      double precision x(n)
      double precision xmult
c
c  The diagonal entries can't be zero.
c
      do i = 1, n
        if ( a(2,i) .eq. 0.0D+00 ) then
          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) 'R83_NP_FS - Fatal error!'
          write ( *, '(a,i8,a)' ) '  A(2,', i, ') = 0.'
          return
        end if
      end do

      do i = 1, n
        x(i) = b(i)
      end do

      do i = 2, n
        xmult = a(3,i-1) / a(2,i-1)
        a(2,i) = a(2,i) - xmult * a(1,i)
        x(i)   = x(i)   - xmult * x(i-1)
      end do

      x(n) = x(n) / a(2,n)
      do i = n-1, 1, -1
        x(i) = ( x(i) - a(1,i+1) * x(i+1) ) / a(2,i)
      end do

      return
      end
      subroutine r83_uniform ( n, seed, a )

c*********************************************************************72
c
cc R83_UNIFORM randomizes an R83 matrix.
c
c  Discussion:
c
c    The R83 storage format is used for a tridiagonal matrix.
c    The superdiagonal is stored in entries (1,2:N), the diagonal in
c    entries (2,1:N), and the subdiagonal in (3,1:N-1).  Thus, the
c    original matrix is "collapsed" vertically into the array.
c
c    Here is how an R83 matrix of order 5 would be stored:
c
c       *  A12 A23 A34 A45
c      A11 A22 A33 A44 A55
c      A21 A32 A43 A54  *
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    19 September 2003
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the order of the linear system.
c
c    Input/output, integer SEED, a seed for the random number
c    generator.
c
c    Output, double precision A(3,N), the R83 matrix.
c
      implicit none

      integer n

      double precision a(3,n)
      integer seed

      a(1,1) = 0.0D+00
      call r8vec_uniform_01 ( n-1, seed, a(1,2:n) )

      call r8vec_uniform_01 ( n,   seed, a(2,1:n) )

      call r8vec_uniform_01 ( n-1, seed, a(3,1:n-1) )
      a(3,n) = 0.0D+00

      return
      end
      subroutine r8vec_bracket ( n, x, xval, left, right )

c*********************************************************************72
c
cc R8VEC_BRACKET searches a sorted array for successive brackets of a value.
c
c  Discussion:
c
c    An R8VEC is a vector of R8's.
c
c    If the values in the vector are thought of as defining intervals
c    on the real line, then this routine searches for the interval
c    nearest to or containing the given value.
c
c  Modified:
c
c    24 July 2008
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, length of input array.
c
c    Input, double precision X(N), an array that has been sorted into
c    ascending order.
c
c    Input, double precision XVAL, a value to be bracketed.
c
c    Output, integer LEFT, RIGHT, the results of the search.
c    Either:
c      XVAL < X(1), when LEFT = 1, RIGHT = 2;
c      X(N) < XVAL, when LEFT = N-1, RIGHT = N;
c    or
c      X(LEFT) <= XVAL <= X(RIGHT).
c
      implicit none

      integer n

      integer i
      integer left
      integer right
      double precision x(n)
      double precision xval

      do i = 2, n - 1

        if ( xval .lt. x(i) ) then
          left = i - 1
          right = i
          return
        end if

       end do

      left = n - 1
      right = n

      return
      end
      subroutine r8vec_bracket3 ( n, t, tval, left )

c*********************************************************************72
c
cc R8VEC_BRACKET3 finds the interval containing or nearest a given value.
c
c  Discussion:
c
c    An R8VEC is a vector of R8's.
c
c    The routine always returns the index LEFT of the sorted array
c    T with the property that either
c    *  T is contained in the interval [ T(LEFT), T(LEFT+1) ], or
c    *  T .lt. T(LEFT) = T(1), or
c    *  T > T(LEFT+1) = T(N).
c
c    The routine is useful for interpolation problems, where
c    the abscissa must be located within an interval of data
c    abscissas for interpolation, or the "nearest" interval
c    to the (extreme) abscissa must be found so that extrapolation
c    can be carried out.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    28 May 2009
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, length of the input array.
c
c    Input, double precision T(N), an array that has been sorted
c    into ascending order.
c
c    Input, double precision TVAL, a value to be bracketed by entries of T.
c
c    Input/output, integer LEFT.
c    On input, if 1 .le. LEFT .le. N-1, LEFT is taken as a suggestion for the
c    interval [ T(LEFT), T(LEFT+1) ] in which TVAL lies.  This interval
c    is searched first, followed by the appropriate interval to the left
c    or right.  After that, a binary search is used.
c    On output, LEFT is set so that the interval [ T(LEFT), T(LEFT+1) ]
c    is the closest to TVAL; it either contains TVAL, or else TVAL
c    lies outside the interval [ T(1), T(N) ].
c
      implicit none

      integer n

      integer high
      integer left
      integer low
      integer mid
      double precision t(n)
      double precision tval
c
c  Check the input data.
c
      if ( n .lt. 2 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'R8VEC_BRACKET3 - Fatal error!'
        write ( *, '(a)' ) '  N must be at least 2.'
        stop
      end if
c
c  If LEFT is not between 1 and N-1, set it to the middle value.
c
      if ( left .lt. 1 .or. n - 1 .lt. left ) then
        left = ( n + 1 ) / 2
      end if
c
c  CASE 1: TVAL .lt. T(LEFT):
c  Search for TVAL in [T(I), T(I+1)] for intervals I = 1 to LEFT-1.
c
      if ( tval .lt. t(left) ) then

        if ( left .eq. 1 ) then
          return
        else if ( left .eq. 2 ) then
          left = 1
          return
        else if ( t(left-1) .le. tval ) then
          left = left - 1
          return
        else if ( tval .le. t(2) ) then
          left = 1
          return
        end if
c
c  ...Binary search for TVAL in [T(I), T(I+1)] for intervals I = 2 to LEFT-2.
c
        low = 2
        high = left - 2

10      continue

          if ( low .eq. high ) then
            left = low
            return
          end if

          mid = ( low + high + 1 ) / 2

          if ( t(mid) .le. tval ) then
            low = mid
          else
            high = mid - 1
          end if

        go to 10
c
c  CASE2: T(LEFT+1) .lt. TVAL:
c  Search for TVAL in [T(I),T(I+1)] for intervals I = LEFT+1 to N-1.
c
      else if ( t(left+1) .lt. tval ) then

        if ( left .eq. n - 1 ) then
          return
        else if ( left .eq. n - 2 ) then
          left = left + 1
          return
        else if ( tval .le. t(left+2) ) then
          left = left + 1
          return
        else if ( t(n-1) .le. tval ) then
          left = n - 1
          return
        end if
c
c  ...Binary search for TVAL in [T(I), T(I+1)] for intervals I = LEFT+2 to N-2.
c
        low = left + 2
        high = n - 2

20      continue

          if ( low .eq. high ) then
            left = low
            return
          end if

          mid = ( low + high + 1 ) / 2

          if ( t(mid) .le. tval ) then
            low = mid
          else
            high = mid - 1
          end if

        go to 20
c
c  CASE3: T(LEFT) .le. TVAL .le. T(LEFT+1):
c  T is in [T(LEFT), T(LEFT+1)], as the user said it might be.
c
      else

      end if

      return
      end
      function r8vec_distinct ( n, a )

c*********************************************************************72
c
cc R8VEC_DISTINCT is true if the entries in an R8VEC are distinct.
c
c  Discussion:
c
c    An R8VEC is a vector of R8's.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    29 May 2009
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the number of entries in the vector.
c
c    Input, double precision A(N), the vector to be checked.
c
c    Output, logical R8VEC_DISTINCT is TRUE if the elements of A are distinct.
c
      implicit none

      integer n

      double precision a(n)
      integer i
      integer j
      logical r8vec_distinct

      r8vec_distinct = .false.

      do i = 2, n
        do j = 1, i - 1
          if ( a(i) .eq. a(j) ) then
            return
          end if
        end do
      end do

      r8vec_distinct = .true.

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
      subroutine r8vec_even ( n, alo, ahi, a )

c*********************************************************************72
c
cc R8VEC_EVEN returns an R8VEC of evenly spaced values.
c
c  Discussion:
c
c    An R8VEC is a vector of R8 values.
c
c    If N is 1, then the midpoint is returned.
c
c    Otherwise, the two endpoints are returned, and N-2 evenly
c    spaced points between them.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    09 December 2004
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the number of values.
c
c    Input, double precision ALO, AHI, the low and high values.
c
c    Output, double precision A(N), N evenly spaced values.
c    Normally, A(1) = ALO and A(N) = AHI.
c    However, if N = 1, then A(1) = 0.5*(ALO+AHI).
c
      implicit none

      integer n

      double precision a(n)
      double precision ahi
      double precision alo
      integer i

      if ( n .eq. 1 ) then

        a(1) = 0.5D+00 * ( alo + ahi )

      else

        do i = 1, n
          a(i) = ( dble ( n - i     ) * alo
     &           + dble (     i - 1 ) * ahi )
     &           / dble ( n     - 1 )
        end do

      end if

      return
      end
      subroutine r8vec_indicator ( n, a )

c*********************************************************************72
c
cc R8VEC_INDICATOR sets an R8VEC to the indicator vector.
c
c  Discussion:
c
c    An R8VEC is a vector of R8's.
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
      subroutine r8vec_order_type ( n, a, order )

c*********************************************************************72
c
cc R8VEC_ORDER_TYPE determines if R8VEC is (non)strictly ascending/descending.
c
c  Discussion:
c
c    An R8VEC is a vector of R8's.
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
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the number of entries of the array.
c
c    Input, double precision A(N), the array to be checked.
c
c    Output, integer ORDER, order indicator:
c    -1, no discernable order;
c    0, all entries are equal;
c    1, ascending order;
c    2, strictly ascending order;
c    3, descending order;
c    4, strictly descending order.
c
      implicit none

      integer n

      double precision a(n)
      integer i
      integer order
c
c  Search for the first value not equal to A(1).
c
      i = 1

10    continue

        i = i + 1

        if ( n .lt. i ) then
          order = 0
          return
        end if

        if ( a(1) .lt. a(i) ) then

          if ( i .eq. 2 ) then
            order = 2
          else
            order = 1
          end if

          go to 20

        else if ( a(i) .lt. a(1) ) then

          if ( i .eq. 2 ) then
            order = 4
          else
            order = 3
          end if

          go to 20

        end if

      go to 10

20    continue
c
c  Now we have a "direction".  Examine subsequent entries.
c
30    continue

      if ( i .lt. n ) then

        i = i + 1

        if ( order .eq. 1 ) then

          if ( a(i) .lt. a(i-1) ) then
            order = -1
            go to 40
          end if

        else if ( order .eq. 2 ) then

          if ( a(i) .lt. a(i-1) ) then
            order = -1
            go to 40
          else if ( a(i) .eq. a(i-1) ) then
            order = 1
          end if

        else if ( order .eq. 3 ) then

          if ( a(i-1) .lt. a(i) ) then
            order = -1
            go to 40
          end if

        else if ( order .eq. 4 ) then

          if ( a(i-1) .lt. a(i) ) then
            order = -1
            go to 40
          else if ( a(i) .eq. a(i-1) ) then
            order = 3
          end if

        end if

        go to 30

      end if

40    continue

      return
      end
      subroutine r8vec_print ( n, a, title )

c*********************************************************************72
c
cc R8VEC_PRINT prints an R8VEC.
c
c  Discussion:
c
c    An R8VEC is a vector of R8's.
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
c    Input, character * ( * ) TITLE, a title.
c
      implicit none

      integer n

      double precision a(n)
      integer i
      character ( len = * ) title

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) trim ( title )
      write ( *, '(a)' ) ' '
      do i = 1, n
        write ( *, '(2x,i8,a,1x,g16.8)' ) i, ':', a(i)
      end do

      return
      end
      subroutine r8vec_sort_bubble_a ( n, a )

c*********************************************************************72
c
cc R8VEC_SORT_BUBBLE_A ascending sorts an R8VEC using bubble sort.
c
c  Discussion:
c
c    An R8VEC is a vector of R8's.
c
c    Bubble sort is simple to program, but inefficient.  It should not
c    be used for large arrays.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    31 May 2009
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the number of entries in the array.
c
c    Input/output, double precision A(N).
c    On input, an unsorted array.
c    On output, the array has been sorted.
c
      implicit none

      integer n

      double precision a(n)
      integer i
      integer  j
      double precision t

      do i = 1, n - 1
        do j = i + 1, n
          if ( a(j) .lt. a(i) ) then
            t    = a(i)
            a(i) = a(j)
            a(j) = t
          end if
        end do
      end do

      return
      end
      subroutine r8vec_uniform_01 ( n, seed, r )

c*********************************************************************72
c
cc R8VEC_UNIFORM_01 returns a unit pseudorandom R8VEC.
c
c  Discussion:
c
c    An R8VEC is a vector of R8's.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    17 July 2006
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
      integer k
      integer seed
      double precision r(n)

      do i = 1, n

        k = seed / 127773

        seed = 16807 * ( seed - k * 127773 ) - k * 2836

        if ( seed .lt. 0 ) then
          seed = seed + 2147483647
        end if

        r(i) = dble ( seed ) * 4.656612875D-10

      end do

      return
      end
      subroutine r8vec_unique_count ( n, a, tol, unique_num )

c*********************************************************************72
c
cc R8VEC_UNIQUE_COUNT counts the unique elements in an unsorted R8VEC.
c
c  Discussion:
c
c    An R8VEC is a vector of R8 values.
c
c    Because the array is unsorted, this algorithm is O(N^2).
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    28 August 2008
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the number of elements of A.
c
c    Input, double precision A(N), the unsorted array to examine.
c
c    Input, double precision TOL, a nonnegative tolerance for equality.
c    Set it to 0.0 for the strictest test.
c
c    Output, integer UNIQUE_NUM, the number of unique elements
c    of A.
c
      implicit none

      integer n

      double precision a(n)
      integer i
      integer j
      integer unique_num
      double precision tol

      unique_num = 0

      do i = 1, n

        unique_num = unique_num + 1

        do j = 1, i - 1

          if ( abs ( a(i) - a(j) ) .le. tol ) then
            unique_num = unique_num - 1
            exit
          end if

        end do

      end do

      return
      end
      subroutine spline_b_val ( ndata, tdata, ydata, tval, yval )

c*********************************************************************72
c
cc SPLINE_B_VAL evaluates a cubic B spline approximant.
c
c  Discussion:
c
c    The cubic B spline will approximate the data, but is not
c    designed to interpolate it.
c
c    In effect, two "phantom" data values are appended to the data,
c    so that the spline will interpolate the first and last data values.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    11 February 2004
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Carl deBoor,
c    A Practical Guide to Splines,
c    Springer, 2001,
c    ISBN: 0387953663.
c
c  Parameters:
c
c    Input, integer NDATA, the number of data values.
c
c    Input, double precision TDATA(NDATA), the abscissas of the data.
c
c    Input, double precision YDATA(NDATA), the data values.
c
c    Input, double precision TVAL, a point at which the spline is
c    to be evaluated.
c
c    Output, double precision YVAL, the value of the function at TVAL.
c
      implicit none

      integer ndata

      double precision bval
      integer left
      integer right
      double precision tdata(ndata)
      double precision tval
      double precision u
      double precision ydata(ndata)
      double precision yval
c
c  Find the nearest interval [ TDATA(LEFT), TDATA(RIGHT) ] to TVAL.
c
      call r8vec_bracket ( ndata, tdata, tval, left, right )
c
c  Evaluate the 5 nonzero B spline basis functions in the interval,
c  weighted by their corresponding data values.
c
      u = ( tval - tdata(left) ) / ( tdata(right) - tdata(left) )
      yval = 0.0D+00
c
c  B function associated with node LEFT - 1, (or "phantom node"),
c  evaluated in its 4th interval.
c
      bval = ( ( (     - 1.0D+00   
     &             * u + 3.0D+00 ) 
     &             * u - 3.0D+00 ) 
     &             * u + 1.0D+00 ) / 6.0D+00

      if ( 0 .lt. left-1 ) then
        yval = yval + ydata(left-1) * bval
      else
        yval = yval + ( 2.0D+00 * ydata(1) - ydata(2) ) * bval
      end if
c
c  B function associated with node LEFT,
c  evaluated in its third interval.
c
      bval = ( ( (       3.0D+00   
     &             * u - 6.0D+00 ) 
     &             * u + 0.0D+00 ) 
     &             * u + 4.0D+00 ) / 6.0D+00

      yval = yval + ydata(left) * bval
c
c  B function associated with node RIGHT,
c  evaluated in its second interval.
c
      bval = ( ( (     - 3.0D+00   
     &             * u + 3.0D+00 ) 
     &             * u + 3.0D+00 ) 
     &             * u + 1.0D+00 ) / 6.0D+00

      yval = yval + ydata(right) * bval
c
c  B function associated with node RIGHT+1, (or "phantom node"),
c  evaluated in its first interval.
c
      bval = u**3 / 6.0D+00

      if ( right+1 .le. ndata ) then
        yval = yval + ydata(right+1) * bval
      else
        yval = yval + ( 2.0D+00 * ydata(ndata) - ydata(ndata-1) ) * bval
      end if

      return
      end
      subroutine spline_beta_val ( beta1, beta2, ndata, tdata, ydata, 
     &  tval, yval )

c*********************************************************************72
c
cc SPLINE_BETA_VAL evaluates a cubic beta spline approximant.
c
c  Discussion:
c
c    The cubic beta spline will approximate the data, but is not
c    designed to interpolate it.
c
c    If BETA1 = 1 and BETA2 = 0, the cubic beta spline will be the
c    same as the cubic B spline approximant.
c
c    With BETA1 = 1 and BETA2 large, the beta spline becomes more like
c    a linear spline.
c
c    In effect, two "phantom" data values are appended to the data,
c    so that the spline will interpolate the first and last data values.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    11 February 2004
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, double precision BETA1, the skew or bias parameter.
c    BETA1 = 1 for no skew or bias.
c
c    Input, double precision BETA2, the tension parameter.
c    BETA2 = 0 for no tension.
c
c    Input, integer NDATA, the number of data values.
c
c    Input, double precision TDATA(NDATA), the abscissas of the data.
c
c    Input, double precision YDATA(NDATA), the data values.
c
c    Input, double precision TVAL, a point at which the spline is
c    to be evaluated.
c
c    Output, double precision YVAL, the value of the function at TVAL.
c
      implicit none

      integer ndata

      double precision a
      double precision b
      double precision beta1
      double precision beta2
      double precision bval
      double precision c
      double precision d
      double precision delta
      integer left
      integer right
      double precision tdata(ndata)
      double precision tval
      double precision u
      double precision ydata(ndata)
      double precision yval
c
c  Find the nearest interval [ TDATA(LEFT), TDATA(RIGHT) ] to TVAL.
c
      call r8vec_bracket ( ndata, tdata, tval, left, right )
c
c  Evaluate the 5 nonzero beta spline basis functions in the interval,
c  weighted by their corresponding data values.
c
      u = ( tval - tdata(left) ) / ( tdata(right) - tdata(left) )

      delta = ( ( 2.0D+00   
     &  * beta1 + 4.0D+00 ) 
     &  * beta1 + 4.0D+00 ) 
     &  * beta1 + 2.0D+00 + beta2

      yval = 0.0D+00
c
c  Beta function associated with node LEFT - 1, (or "phantom node"),
c  evaluated in its 4th interval.
c
      bval = 2.0D+00 * ( beta1 * ( 1.0D+00 - u ) )**3  / delta

      if ( 0 .lt. left-1 ) then
        yval = yval + ydata(left-1) * bval
      else
        yval = yval + ( 2.0D+00 * ydata(1) - ydata(2) ) * bval
      end if
c
c  Beta function associated with node LEFT,
c  evaluated in its third interval.
c
      a = beta2 + ( 4.0D+00 + 4.0D+00 * beta1 ) * beta1

      b = - 6.0D+00 * beta1 * ( 1.0D+00 - beta1 ) * ( 1.0D+00 + beta1 )

      c = ( (     - 6.0D+00   
     &    * beta1 - 6.0D+00 ) 
     &    * beta1 + 0.0D+00 ) 
     &    * beta1 - 3.0D+00 * beta2

      d = ( (     + 2.0D+00   
     &    * beta1 + 2.0D+00 ) 
     &    * beta1 + 2.0D+00 ) 
     &    * beta1 + 2.0D+00 * beta2

      bval = ( a + u * ( b + u * ( c + u * d ) ) ) / delta

      yval = yval + ydata(left) * bval
c
c  Beta function associated with node RIGHT,
c  evaluated in its second interval.
c
      a = 2.0D+00

      b = 6.0D+00 * beta1

      c = 3.0D+00 * beta2 + 6.0D+00 * beta1 * beta1

      d = - 2.0D+00 * ( 1.0D+00 + beta2 + beta1 + beta1 * beta1 )

      bval = ( a + u * ( b + u * ( c + u * d ) ) ) / delta

      yval = yval + ydata(right) * bval
c
c  Beta function associated with node RIGHT+1, (or "phantom node"),
c  evaluated in its first interval.
c
      bval = 2.0D+00 * u**3 / delta

      if ( right + 1 .le. ndata ) then
        yval = yval + ydata(right+1) * bval
      else
        yval = yval + ( 2.0D+00 * ydata(ndata) - ydata(ndata-1) ) * bval
      end if

      return
      end
      subroutine spline_bezier_val ( dim_num, interval_num, data_val, 
     &  point_num, point_t, point_val )

c*********************************************************************72
c
cc SPLINE_BEZIER_VAL evaluates a cubic Bezier spline.
c
c  Discussion:
c
c    The cubic Bezier spline of N parts is defined by choosing
c    3*N+1 equally spaced T-abscissa values in the interval [0,N].
c    This defines N subintervals, each of length 1, and each containing
c    4 successives T abscissa values.
c
c    At each abscissa value, a DIM_NUM-dimensional Bezier control
c    value is assigned.  Over each interval, a Bezier cubic function
c    is used to define the value of the Bezier spline.  To the left of
c    the first interval, or to the right of the last interval,
c    extrapolation may be used to extend the spline definition to
c    the entire real line.
c
c    Note that the Bezier spline will pass through the 1st, 4th,
c    and in general 3*I+1 control values exactly.  The other control
c    values are not interpolating points.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    17 June 2006
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer DIM_NUM, the spatial dimension.
c
c    Input, integer INTERVAL_NUM, the number of intervals.
c
c    Input, double precision DATA_VAL(DIM_NUM,3*INTERVAL_NUM+1), the control
c    values.
c
c    Input, integer POINT_NUM, the number of sample points at
c    which the Bezier cubic spline is to be evaluated.
c
c    Input, double precision POINT_T(POINT_NUM), the "T" values associated
c    with the points.  A value of T between 0 and 1, for instance,
c    is associated with the first interval, and a value of T between
c    INTERVAL_NUM-1 and INTERVAL_NUM is in the last interval.
c
c    Output, double precision POINT_VAL(DIM_NUM,POINT_NUM), the value
c    of the Bezier cubic spline at the sample points.
c
      implicit none

      integer cubic 
      parameter ( cubic = 3 )
      integer interval_num
      integer dim_num
      integer point_num

      double precision bernstein_val(0:cubic)
      double precision d
      double precision data_val(dim_num,cubic*interval_num+1)
      integer dim
      integer interval
      integer j
      integer offset
      integer point
      double precision point_t(point_num)
      double precision point_val(dim_num,point_num)
      double precision t
      double precision t_01

      do point = 1, point_num

        t = point_t(point)

        interval = int ( t + 1 )

        interval = max ( interval, 1 )
        interval = min ( interval, interval_num )

        offset = 1 + ( interval - 1 ) * cubic

        t_01 = t - dble ( interval - 1 )

        call bp01 ( cubic, t_01, bernstein_val )

        do dim = 1, dim_num
          d = 0.0D+00
          do j = 0, cubic
            d = d + data_val(dim,offset+j) * bernstein_val(j)
          end do
          point_val(dim,point) = d
        end do

      end do

      return
      end
      subroutine spline_constant_val ( ndata, tdata, ydata, tval, yval )

c*********************************************************************72
c
cc SPLINE_CONSTANT_VAL evaluates a piecewise constant spline at a point.
c
c  Discussion:
c
c    NDATA-1 points TDATA define NDATA intervals, with the first
c    and last being semi-infinite.
c
c    The value of the spline anywhere in interval I is YDATA(I).
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    16 November 2001
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer NDATA, the number of data points defining
c    the spline.  NDATA must be at least 1.
c
c    Input, double precision TDATA(NDATA-1), the breakpoints.  The values
c    of TDATA should be distinct and increasing.
c
c    Input, double precision YDATA(NDATA), the values of the spline in
c    the intervals defined by the breakpoints.
c
c    Input, double precision TVAL, the point at which the spline is
c    to be evaluated.
c
c    Output, double precision YVAL, the value of the spline at TVAL.
c
      implicit none

      integer ndata

      integer i
      double precision tdata(ndata-1)
      double precision tval
      double precision ydata(ndata)
      double precision yval
c
c  Check NDATA.
c
      if ( ndata .lt. 1 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'SPLINE_CONSTANT_VAL - Fatal error!'
        write ( *, '(a)' ) '  NDATA .lt. 1.'
        stop
      end if

      do i = 1, ndata-1
        if ( tval .le. tdata(i) ) then
          yval = ydata(i)
          return
        end if
      end do

      yval = ydata(ndata)

      return
      end
      subroutine spline_cubic_set ( n, t, y, ibcbeg, ybcbeg, ibcend, 
     &  ybcend, ypp )

c*********************************************************************72
c
cc SPLINE_CUBIC_SET computes the second derivatives of a piecewise cubic spline.
c
c  Discussion:
c
c    For data interpolation, the user must call SPLINE_CUBIC_SET to
c    determine the second derivative data, passing in the data to be
c    interpolated, and the desired boundary conditions.
c
c    The data to be interpolated, plus the SPLINE_CUBIC_SET output,
c    defines the spline.  The user may then call SPLINE_CUBIC_VAL to
c    evaluate the spline at any point.
c
c    The cubic spline is a piecewise cubic polynomial.  The intervals
c    are determined by the "knots" or abscissas of the data to be
c    interpolated.  The cubic spline has continous first and second
c    derivatives over the entire interval of interpolation.
c
c    For any point T in the interval T(IVAL), T(IVAL+1), the form of
c    the spline is
c
c      SPL(T) = A(IVAL)
c             + B(IVAL) * ( T - T(IVAL) )
c             + C(IVAL) * ( T - T(IVAL) )**2
c             + D(IVAL) * ( T - T(IVAL) )**3
c
c    If we assume that we know the values Y(*) and YPP(*), which represent
c    the values and second derivatives of the spline at each knot, then
c    the coefficients can be computed as:
c
c      A(IVAL) = Y(IVAL)
c      B(IVAL) = ( Y(IVAL+1) - Y(IVAL) ) / ( T(IVAL+1) - T(IVAL) )
c        - ( YPP(IVAL+1) + 2 * YPP(IVAL) ) * ( T(IVAL+1) - T(IVAL) ) / 6
c      C(IVAL) = YPP(IVAL) / 2
c      D(IVAL) = ( YPP(IVAL+1) - YPP(IVAL) ) / ( 6 * ( T(IVAL+1) - T(IVAL) ) )
c
c    Since the first derivative of the spline is
c
c      SPL'(T) =     B(IVAL)
c              + 2 * C(IVAL) * ( T - T(IVAL) )
c              + 3 * D(IVAL) * ( T - T(IVAL) )**2,
c
c    the requirement that the first derivative be continuous at interior
c    knot I results in a total of N-2 equations, of the form:
c
c      B(IVAL-1) + 2 C(IVAL-1) * (T(IVAL)-T(IVAL-1))
c      + 3 * D(IVAL-1) * (T(IVAL) - T(IVAL-1))**2 = B(IVAL)
c
c    or, setting H(IVAL) = T(IVAL+1) - T(IVAL)
c
c      ( Y(IVAL) - Y(IVAL-1) ) / H(IVAL-1)
c      - ( YPP(IVAL) + 2 * YPP(IVAL-1) ) * H(IVAL-1) / 6
c      + YPP(IVAL-1) * H(IVAL-1)
c      + ( YPP(IVAL) - YPP(IVAL-1) ) * H(IVAL-1) / 2
c      =
c      ( Y(IVAL+1) - Y(IVAL) ) / H(IVAL)
c      - ( YPP(IVAL+1) + 2 * YPP(IVAL) ) * H(IVAL) / 6
c
c    or
c
c      YPP(IVAL-1) * H(IVAL-1) + 2 * YPP(IVAL) * ( H(IVAL-1) + H(IVAL) )
c      + YPP(IVAL) * H(IVAL)
c      =
c      6 * ( Y(IVAL+1) - Y(IVAL) ) / H(IVAL)
c      - 6 * ( Y(IVAL) - Y(IVAL-1) ) / H(IVAL-1)
c
c    Boundary conditions must be applied at the first and last knots.
c    The resulting tridiagonal system can be solved for the YPP values.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    27 January 2004
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Carl deBoor,
c    A Practical Guide to Splines,
c    Springer, 2001,
c    ISBN: 0387953663.
c
c  Parameters:
c
c    Input, integer N, the number of data points; N must be
c    at least 2.
c
c    Input, double precision T(N), the points where data is specified.
c    The values should be distinct, and increasing.
c
c    Input, double precision Y(N), the data values to be interpolated.
c
c    Input, integer IBCBEG, the left boundary condition flag:
c
c      0: the spline should be a quadratic over the first interval;
c      1: the first derivative at the left endpoint should be YBCBEG;
c      2: the second derivative at the left endpoint should be YBCBEG.
c
c    Input, double precision YBCBEG, the left boundary value, if needed.
c
c    Input, integer IBCEND, the right boundary condition flag:
c
c      0: the spline should be a quadratic over the last interval;
c      1: the first derivative at the right endpoint should be YBCEND;
c      2: the second derivative at the right endpoint should be YBCEND.
c
c    Input, double precision YBCEND, the right boundary value, if needed.
c
c    Output, double precision YPP(N), the second derivatives of
c    the cubic spline.
c
      implicit none

      integer n

      double precision a(3,n)
      integer i
      integer ibcbeg
      integer ibcend
      double precision t(n)
      double precision y(n)
      double precision ybcbeg
      double precision ybcend
      double precision ypp(n)
c
c  Check.
c
      if ( n .le. 1 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'SPLINE_CUBIC_SET - Fatal error!'
        write ( *, '(a)' ) '  The number of knots must be at least 2.'
        write ( *, '(a,i8)' ) '  The input value of N = ', n
        stop
      end if

      do i = 1, n-1
        if ( t(i+1) .le. t(i) ) then
          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) 'SPLINE_CUBIC_SET - Fatal error!'
          write ( *, '(a)' ) 
     &      '  The knots must be strictly increasing, but'
          write ( *, '(a,i8,a,g14.6)' ) '  T(',  i,') = ', t(i)
          write ( *, '(a,i8,a,g14.6)' ) '  T(',i+1,') = ', t(i+1)
          stop
        end if
      end do
c
c  Set the first equation.
c
      if ( ibcbeg .eq. 0 ) then
        ypp(1) = 0.0D+00
        a(2,1) = 1.0D+00
        a(1,2) = -1.0D+00
      else if ( ibcbeg .eq. 1 ) then
        ypp(1) = ( y(2) - y(1) ) / ( t(2) - t(1) ) - ybcbeg
        a(2,1) = ( t(2) - t(1) ) / 3.0D+00
        a(1,2) = ( t(2) - t(1) ) / 6.0D+00
      else if ( ibcbeg .eq. 2 ) then
        ypp(1) = ybcbeg
        a(2,1) = 1.0D+00
        a(1,2) = 0.0D+00
      else
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'SPLINE_CUBIC_SET - Fatal error!'
        write ( *, '(a)' ) 
     &    '  The boundary flag IBCBEG must be 0, 1 or 2.'
        write ( *, '(a,i8)' ) '  The input value is IBCBEG = ', ibcbeg
        stop
      end if
c
c  Set the intermediate equations.
c
      do i = 2, n-1
        ypp(i) = ( y(i+1) - y(i) ) / ( t(i+1) - t(i) ) 
     &         - ( y(i) - y(i-1) ) / ( t(i) - t(i-1) )
        a(3,i-1) = ( t(i) - t(i-1) ) / 6.0D+00
        a(2,i) = ( t(i+1) - t(i-1) ) / 3.0D+00
        a(1,i+1) = ( t(i+1) - t(i) ) / 6.0D+00
      end do
c
c  Set the last equation.
c
      if ( ibcend .eq. 0 ) then
        ypp(n) = 0.0D+00
        a(3,n-1) = -1.0D+00
        a(2,n) = 1.0D+00
      else if ( ibcend .eq. 1 ) then
        ypp(n) = ybcend - ( y(n) - y(n-1) ) / ( t(n) - t(n-1) )
        a(3,n-1) = ( t(n) - t(n-1) ) / 6.0D+00
        a(2,n) = ( t(n) - t(n-1) ) / 3.0D+00
      else if ( ibcend .eq. 2 ) then
        ypp(n) = ybcend
        a(3,n-1) = 0.0D+00
        a(2,n) = 1.0D+00
      else
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'SPLINE_CUBIC_SET - Fatal error!'
        write ( *, '(a)' ) 
     &    '  The boundary flag IBCEND must be 0, 1 or 2.'
        write ( *, '(a,i8)' ) '  The input value is IBCEND = ', ibcend
        stop
      end if
c
c  Special case:
c    N = 2, IBCBEG = IBCEND = 0.
c
      if ( n .eq. 2 .and. ibcbeg .eq. 0 .and. ibcend .eq. 0 ) then

        ypp(1) = 0.0D+00
        ypp(2) = 0.0D+00
c
c  Solve the linear system.
c
      else

        call r83_np_fs ( n, a, ypp, ypp )

      end if

      return
      end
      subroutine spline_cubic_val ( n, t, y, ypp, tval, yval, ypval, 
     &  yppval )

c*********************************************************************72
c
cc SPLINE_CUBIC_VAL evaluates a piecewise cubic spline at a point.
c
c  Discussion:
c
c    SPLINE_CUBIC_SET must have already been called to define the
c    values of YPP.
c
c    For any point T in the interval T(IVAL), T(IVAL+1), the form of
c    the spline is
c
c      SPL(T) = A
c             + B * ( T - T(IVAL) )
c             + C * ( T - T(IVAL) )**2
c             + D * ( T - T(IVAL) )**3
c
c    Here:
c      A = Y(IVAL)
c      B = ( Y(IVAL+1) - Y(IVAL) ) / ( T(IVAL+1) - T(IVAL) )
c        - ( YPP(IVAL+1) + 2 * YPP(IVAL) ) * ( T(IVAL+1) - T(IVAL) ) / 6
c      C = YPP(IVAL) / 2
c      D = ( YPP(IVAL+1) - YPP(IVAL) ) / ( 6 * ( T(IVAL+1) - T(IVAL) ) )
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    20 November 2000
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Carl deBoor,
c    A Practical Guide to Splines,
c    Springer, 2001,
c    ISBN: 0387953663.
c
c  Parameters:
c
c    Input, integer N, the number of data values.
c
c    Input, double precision T(N), the knot values.
c
c    Input, double precision Y(N), the data values at the knots.
c
c    Input, double precision YPP(N), the second derivatives of the
c    spline at the knots.
c
c    Input, double precision TVAL, a point, typically between T(1) and
c    T(N), at which the spline is to be evalulated.  If TVAL lies outside
c    this range, extrapolation is used.
c
c    Output, double precision YVAL, YPVAL, YPPVAL, the value of the spline, and
c    its first two derivatives at TVAL.
c
      implicit none

      integer n

      double precision dt
      double precision h
      integer left
      integer right
      double precision t(n)
      double precision tval
      double precision y(n)
      double precision ypp(n)
      double precision yppval
      double precision ypval
      double precision yval
c
c  Determine the interval [T(LEFT), T(RIGHT)] that contains TVAL.
c  Values below T(1) or above T(N) use extrapolation.
c
      call r8vec_bracket ( n, t, tval, left, right )
c
c  Evaluate the polynomial.
c
      dt = tval - t(left)
      h = t(right) - t(left)

      yval = y(left) 
     &     + dt * ( ( y(right) - y(left) ) / h 
     &            - ( ypp(right) / 6.0D+00 + ypp(left) / 3.0D+00 ) * h 
     &     + dt * ( 0.5D+00 * ypp(left) 
     &     + dt * ( ( ypp(right) - ypp(left) ) / ( 6.0D+00 * h ) ) ) )

      ypval = ( y(right) - y(left) ) / h 
     &     - ( ypp(right) / 6.0D+00 + ypp(left) / 3.0D+00 ) * h 
     &     + dt * ( ypp(left) 
     &     + dt * ( 0.5D+00 * ( ypp(right) - ypp(left) ) / h ) )

      yppval = ypp(left) + dt * ( ypp(right) - ypp(left) ) / h

      return
      end
      subroutine spline_cubic_val2 ( n, t, y, ypp, left, tval, yval, 
     &  ypval, yppval )

c*********************************************************************72
c
cc SPLINE_CUBIC_VAL2 evaluates a piecewise cubic spline at a point.
c
c  Discussion:
c
c    This routine is a modification of SPLINE_CUBIC_VAL; it allows the
c    user to speed up the code by suggesting the appropriate T interval
c    to search first.
c
c    SPLINE_CUBIC_SET must have already been called to define the
c    values of YPP.
c
c    In the LEFT interval, let RIGHT = LEFT+1.  The form of the spline is
c
c      SPL(T) =
c          A
c        + B * ( T - T(LEFT) )
c        + C * ( T - T(LEFT) )**2
c        + D * ( T - T(LEFT) )**3
c
c    Here:
c      A = Y(LEFT)
c      B = ( Y(RIGHT) - Y(LEFT) ) / ( T(RIGHT) - T(LEFT) )
c        - ( YPP(RIGHT) + 2 * YPP(LEFT) ) * ( T(RIGHT) - T(LEFT) ) / 6
c      C = YPP(LEFT) / 2
c      D = ( YPP(RIGHT) - YPP(LEFT) ) / ( 6 * ( T(RIGHT) - T(LEFT) ) )
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    24 February 2004
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Carl deBoor,
c    A Practical Guide to Splines,
c    Springer, 2001,
c    ISBN: 0387953663.
c
c  Parameters:
c
c    Input, integer N, the number of knots.
c
c    Input, double precision T(N), the knot values.
c
c    Input, double precision Y(N), the data values at the knots.
c
c    Input, double precision YPP(N), the second derivatives of the spline at
c    the knots.
c
c    Input/output, integer LEFT, the suggested T interval to 
c    search.  LEFT should be between 1 and N-1.  If LEFT is not in this range,
c    then its value will be ignored.  On output, LEFT is set to the
c    actual interval in which TVAL lies.
c
c    Input, double precision TVAL, a point, typically between T(1) and T(N), at
c    which the spline is to be evalulated.  If TVAL lies outside
c    this range, extrapolation is used.
c
c    Output, double precision YVAL, YPVAL, YPPVAL, the value of the spline, and
c    its first two derivatives at TVAL.
c
      implicit none

      integer n

      double precision dt
      double precision h
      integer left
      integer right
      double precision t(n)
      double precision tval
      double precision y(n)
      double precision ypp(n)
      double precision yppval
      double precision ypval
      double precision yval
c
c  Determine the interval [T(LEFT), T(RIGHT)] that contains TVAL.
c
c  What you want from R8VEC_BRACKET3 is that TVAL is to be computed
c  by the data in interval [T(LEFT), T(RIGHT)].
c
      call r8vec_bracket3 ( n, t, tval, left )
      right = left + 1
c
c  In the interval LEFT, the polynomial is in terms of a normalized
c  coordinate  ( DT / H ) between 0 and 1.
c
      dt = tval - t(left)
      h = t(right) - t(left)

      yval = y(left) + dt * ( ( y(right) - y(left) ) / h 
     &            - ( ypp(right) / 6.0D+00 + ypp(left) / 3.0D+00 ) * h 
     &     + dt * ( 0.5D+00 * ypp(left) 
     &     + dt * ( ( ypp(right) - ypp(left) ) / ( 6.0D+00 * h ) ) ) )

      ypval = ( y(right) - y(left) ) / h 
     &    - ( ypp(right) / 6.0D+00 + ypp(left) / 3.0D+00 ) * h 
     &    + dt * ( ypp(left) 
     &    + dt * ( 0.5D+00 * ( ypp(right) - ypp(left) ) / h ) )

      yppval = ypp(left) + dt * ( ypp(right) - ypp(left) ) / h

      return
      end
      subroutine spline_hermite_set ( ndata, tdata, ydata, ypdata, c )

c*********************************************************************72
c
cc SPLINE_HERMITE_SET sets up a piecewise cubic Hermite interpolant.
c
c  Discussion:
c
c    Once the array C is computed, then in the interval
c    (TDATA(I), TDATA(I+1)), the interpolating Hermite polynomial
c    is given by
c
c      SVAL(TVAL) =                 C(1,I)
c         + ( TVAL - TDATA(I) ) * ( C(2,I)
c         + ( TVAL - TDATA(I) ) * ( C(3,I)
c         + ( TVAL - TDATA(I) ) *   C(4,I) ) )
c
c    This is algorithm CALCCF of Conte and deBoor.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    11 February 2004
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Samuel Conte, Carl deBoor,
c    Elementary Numerical Analysis,
c    Second Edition,
c    McGraw Hill, 1972,
c    ISBN: 07-012446-4,
c    LC: QA297.C65.
c
c  Parameters:
c
c    Input, integer NDATA, the number of data points.
c    NDATA must be at least 2.
c
c    Input, double precision TDATA(NDATA), the abscissas of the data points.
c    The entries of TDATA are assumed to be strictly increasing.
c
c    Input, double precision Y(NDATA), YP(NDATA), the value of the
c    function and its derivative at TDATA(1:NDATA).
c
c    Output, double precision C(4,NDATA), the coefficients of the
c    Hermite polynomial.
c    C(1,1:NDATA) = Y(1:NDATA) and C(2,1:NDATA) = YP(1:NDATA).
c    C(3,1:NDATA-1) and C(4,1:NDATA-1) are the quadratic and cubic
c    coefficients.
c
      implicit none

      integer ndata

      double precision c(4,ndata)
      double precision divdif1
      double precision divdif3
      double precision dt
      integer i
      double precision tdata(ndata)
      double precision ydata(ndata)
      double precision ypdata(ndata)
c
c  Check NDATA.
c
      if ( ndata .lt. 2 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'SPLINE_HERMITE_SET - Fatal error!'
        write ( *, '(a)' ) '  NDATA .lt. 2.'
        stop
      end if

      do i = 1, ndata
        c(1,i) = ydata(i)
        c(2,i) = ypdata(i)
      end do

      do i = 1, ndata-1
        dt = tdata(i+1) - tdata(i)
        divdif1 = ( c(1,i+1) - c(1,i) ) / dt
        divdif3 = c(2,i) + c(2,i+1) - 2.0D+00 * divdif1
        c(3,i) = ( divdif1 - c(2,i) - divdif3 ) / dt
        c(4,i) = divdif3 / ( dt * dt )
      end do

      c(3,ndata) = 0.0D+00
      c(4,ndata) = 0.0D+00

      return
      end
      subroutine spline_hermite_val ( ndata, tdata, c, tval, sval, 
     &  spval )

c*********************************************************************72
c
cc SPLINE_HERMITE_VAL evaluates a piecewise cubic Hermite interpolant.
c
c  Discussion:
c
c    SPLINE_HERMITE_SET must be called first, to set up the
c    spline data from the raw function and derivative data.
c
c    In the interval (TDATA(I), TDATA(I+1)), the interpolating
c    Hermite polynomial is given by
c
c      SVAL(TVAL) =                 C(1,I)
c         + ( TVAL - TDATA(I) ) * ( C(2,I)
c         + ( TVAL - TDATA(I) ) * ( C(3,I)
c         + ( TVAL - TDATA(I) ) *   C(4,I) ) )
c
c    and
c
c      SVAL'(TVAL) =                    C(2,I)
c         + ( TVAL - TDATA(I) ) * ( 2 * C(3,I)
c         + ( TVAL - TDATA(I) ) *   3 * C(4,I) )
c
c    This is algorithm PCUBIC of Conte and deBoor.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    11 February 2004
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Samuel Conte, Carl deBoor,
c    Elementary Numerical Analysis,
c    Second Edition,
c    McGraw Hill, 1972,
c    ISBN: 07-012446-4,
c    LC: QA297.C65.
c
c  Parameters:
c
c    Input, integer NDATA, the number of data points.
c    NDATA must be at least 2.
c
c    Input, double precision TDATA(NDATA), the abscissas of the data points.
c    The entries of TDATA are assumed to be strictly increasing.
c
c    Input, double precision C(4,NDATA), the coefficient data computed by
c    SPLINE_HERMITE_SET.
c
c    Input, double precision TVAL, the point where the interpolant is to
c    be evaluated.
c
c    Output, double precision SVAL, SPVAL, the value of the interpolant
c    and its derivative at TVAL.
c
      implicit none

      integer ndata

      double precision c(4,ndata)
      double precision dt
      integer left
      integer right
      double precision spval
      double precision sval
      double precision tdata(ndata)
      double precision tval
c
c  Check NDATA.
c
      if ( ndata .lt. 2 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'SPLINE_HERMITE_VAL - Fatal error!'
        write ( *, '(a)' ) '  NDATA .lt. 2.'
        stop
      end if
c
c  Find the interval [ TDATA(LEFT), TDATA(RIGHT) ] that contains
c  or is nearest to TVAL.
c
      call r8vec_bracket ( ndata, tdata, tval, left, right )
c
c  Evaluate the cubic polynomial.
c
      dt = tval - tdata(left)

      sval = c(1,left) + dt * ( 
     &       c(2,left) + dt * ( 
     &       c(3,left) + dt * 
     &       c(4,left) ) )

      spval = c(2,left) + dt * ( 2.0D+00 * 
     &        c(3,left) + dt *   3.0D+00 * 
     &        c(4,left) )

      return
      end
      subroutine spline_linear_int ( ndata, tdata, ydata, a, b, 
     &  int_val )

c*********************************************************************72
c
cc SPLINE_LINEAR_INT evaluates the integral of a piecewise linear spline.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    01 November 2001
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer NDATA, the number of data points defining
c    the spline.  NDATA must be at least 2.
c
c    Input, double precision TDATA(NDATA), YDATA(NDATA), the values of
c    the independent and dependent variables at the data points.  The
c    values of TDATA should be distinct and increasing.
c
c    Input, double precision A, B, the interval over which the
c    integral is desired.
c
c    Output, double precision INT_VAL, the value of the integral.
c
      implicit none

      integer ndata

      double precision a
      double precision a_copy
      integer a_left
      integer a_right
      double precision b
      double precision b_copy
      integer b_left
      integer b_right
      integer i_left
      double precision int_val
      double precision tdata(ndata)
      double precision tval
      double precision ydata(ndata)
      double precision yp
      double precision yval
c
c  Check NDATA.
c
      if ( ndata .lt. 2 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'SPLINE_LINEAR_INT - Fatal error!'
        write ( *, '(a)' ) '  NDATA .lt. 2.'
        stop
      end if

      int_val = 0.0D+00

      if ( a .eq. b ) then
        return
      end if

      a_copy = min ( a, b )
      b_copy = max ( a, b )
c
c  Find the interval [ TDATA(A_LEFT), TDATA(A_RIGHT) ] that contains, or is
c  nearest to, A.
c
      call r8vec_bracket ( ndata, tdata, a_copy, a_left, a_right )
c
c  Find the interval [ TDATA(B_LEFT), TDATA(B_RIGHT) ] that contains, or is
c  nearest to, B.
c
      call r8vec_bracket ( ndata, tdata, b_copy, b_left, b_right )
c
c  If A and B are in the same interval...
c
      if ( a_left .eq. b_left ) then

        tval = ( a_copy + b_copy ) / 2.0D+00

        yp = ( ydata(a_right) - ydata(a_left) ) / 
     &       ( tdata(a_right) - tdata(a_left) )

        yval = ydata(a_left) + ( tval - tdata(a_left) ) * yp

        int_val = yval * ( b_copy - a_copy )

        return
      end if
c
c  Otherwise, integrate from:
c
c  A               to TDATA(A_RIGHT),
c  TDATA(A_RIGHT)  to TDATA(A_RIGHT+1),...
c  TDATA(B_LEFT-1) to TDATA(B_LEFT),
c  TDATA(B_LEFT)   to B.
c
c  Use the fact that the integral of a linear function is the
c  value of the function at the midpoint times the width of the interval.
c
      tval = ( a_copy + tdata(a_right) ) / 2.0D+00

      yp = ( ydata(a_right) - ydata(a_left) ) / 
     &     ( tdata(a_right) - tdata(a_left) )

      yval = ydata(a_left) + ( tval - tdata(a_left) ) * yp

      int_val = int_val + yval * ( tdata(a_right) - a_copy )

      do i_left = a_right, b_left - 1

        tval = ( tdata(i_left+1) + tdata(i_left) ) / 2.0D+00

        yp = ( ydata(i_left+1) - ydata(i_left) ) / 
     &       ( tdata(i_left+1) - tdata(i_left) )

        yval = ydata(i_left) + ( tval - tdata(i_left) ) * yp

        int_val = int_val + yval * ( tdata(i_left + 1) - tdata(i_left) )

      end do

      tval = ( tdata(b_left) + b_copy ) / 2.0D+00

      yp = ( ydata(b_right) - ydata(b_left) ) / 
     &     ( tdata(b_right) - tdata(b_left) )

      yval = ydata(b_left) + ( tval - tdata(b_left) ) * yp

      int_val = int_val + yval * ( b_copy - tdata(b_left) )

      if ( b .lt. a ) then
        int_val = - int_val
      end if

      return
      end
      subroutine spline_linear_intset ( n, int_x, int_v, data_x, 
     &  data_y )

c*********************************************************************72
c
cc SPLINE_LINEAR_INTSET: piecewise linear spline with given integral properties.
c
c  Discussion:
c
c    The user has in mind an interval, divided by INT_N+1 points into
c    INT_N intervals.  A linear spline is to be constructed,
c    with breakpoints at the centers of each interval, and extending
c    continuously to the left of the first and right of the last
c    breakpoints.  The constraint on the linear spline is that it is
c    required that it have a given integral value over each interval.
c
c    A tridiagonal linear system of equations is solved for the
c    values of the spline at the breakpoints.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    27 January 2004
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the number of intervals.
c
c    Input, double precision INT_X(N+1), the points that define the intervals.
c    Interval I lies between INT_X(I) and INT_X(I+1).
c
c    Input, double precision INT_V(N), the desired value of the integral of the
c    linear spline over each interval.
c
c    Output, double precision DATA_X(N), DATA_Y(N), the values of the
c    independent and dependent variables at the data points.  The values
c    of DATA_X are the interval midpoints.  The values of DATA_Y are
c    determined in such a way that the exact integral of the linear
c    spline over interval I is equal to INT_V(I).
c
      implicit none

      integer n

      double precision a(3,n)
      double precision data_x(n)
      double precision data_y(n)
      integer i
      double precision int_v(n)
      double precision int_x(n+1)
c
c  Set up the easy stuff.
c
      do i = 1, n
        data_x(i) = 0.5D+00 * ( int_x(i) + int_x(i+1) )
      end do
c
c  Set up C, D, E, the coefficients of the linear system.
c
      do i = 1, n - 2
        a(3,i) = 1.0D+00 
     &    - ( 0.5D+00 * ( data_x(i+1) + int_x(i+1) ) - data_x(i) ) 
     &    / ( data_x(i+1) - data_x(i) )
      end do
      a(3,n-1) = 0.0D+00
      a(3,n) = 0.0D+00

      a(2,1) = int_x(2) - int_x(1)

      do i = 2, n - 1
        a(2,i) = 1.0D+00 
     &    + ( 0.5D+00 * ( data_x(i) + int_x(i) ) 
     &    - data_x(i-1) ) 
     &    / ( data_x(i) - data_x(i-1) ) 
     &    - ( 0.5D+00 * ( data_x(i) + int_x(i+1) ) - data_x(i) ) 
     &    / ( data_x(i+1) - data_x(i) )
      end do

      a(2,n) = int_x(n+1) - int_x(n)

      a(1,1) = 0.0D+00
      a(1,2) = 0.0D+00

      do i = 3, n
        a(1,i) = ( 0.5D+00 * ( data_x(i-1) + int_x(i) ) 
     &    - data_x(i-1) ) / ( data_x(i) - data_x(i-1) )
      end do
c
c  Set up DATA_Y, which begins as the right hand side of the linear system.
c
      data_y(1) = int_v(1)
      do i = 2, n - 1
        data_y(i) = 2.0D+00 * int_v(i) / ( int_x(i+1) - int_x(i) )
      end do
      data_y(n) = int_v(n)
c
c  Solve the linear system.
c
      call r83_np_fs ( n, a, data_y, data_y )

      return
      end
      subroutine spline_linear_val ( ndata, tdata, ydata, tval, yval, 
     &  ypval )

c*********************************************************************72
c
cc SPLINE_LINEAR_VAL evaluates a piecewise linear spline at a point.
c
c  Discussion:
c
c    Because of the extremely simple form of the linear spline,
c    the raw data points ( TDATA(I), YDATA(I)) can be used directly to
c    evaluate the spline at any point.  No processing of the data
c    is required.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    06 April 1999
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer NDATA, the number of data points defining
c    the spline.  NDATA must be at least 2.
c
c    Input, double precision TDATA(NDATA), YDATA(NDATA), the values of
c    the independent and dependent variables at the data points.  The
c    values of TDATA should be distinct and increasing.
c
c    Input, double precision TVAL, the point at which the spline is
c    to be evaluated.
c
c    Output, double precision YVAL, YPVAL, the value of the spline and
c    its first derivative dYdT at TVAL.  YPVAL is not reliable if TVAL
c    is exactly equal to TDATA(I) for some I.
c
      implicit none

      integer ndata

      integer left
      integer right
      double precision tdata(ndata)
      double precision tval
      double precision ydata(ndata)
      double precision ypval
      double precision yval
c
c  Check NDATA.
c
      if ( ndata .lt. 2 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'SPLINE_LINEAR_VAL - Fatal error!'
        write ( *, '(a)' ) '  NDATA .lt. 2.'
        stop
      end if
c
c  Find the interval [ TDATA(LEFT), TDATA(RIGHT) ] that contains, or is
c  nearest to, TVAL.
c
      call r8vec_bracket ( ndata, tdata, tval, left, right )
c
c  Now evaluate the piecewise linear function.
c
      ypval = ( ydata(right) - ydata(left) ) 
     &  / ( tdata(right) - tdata(left) )

      yval = ydata(left) +  ( tval - tdata(left) ) * ypval

      return
      end
      subroutine spline_overhauser_nonuni_val ( ndata, tdata, ydata, 
     &  tval, yval )

c*********************************************************************72
c
cc SPLINE_OVERHAUSER_NONUNI_VAL evaluates the nonuniform Overhauser spline.
c
c  Discussion:
c
c    The nonuniformity refers to the fact that the abscissa values
c    need not be uniformly spaced.
c
c    Thanks to Doug Fortune for pointing out that the point distances
c    used to define ALPHA and BETA should be the Euclidean distances
c    between the points.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    08 May 2007
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer NDATA, the number of data points.
c    3 .le. NDATA is required.
c
c    Input, double precision TDATA(NDATA), the abscissas of the data points.
c    The values of TDATA are assumed to be distinct and increasing.
c
c    Input, double precision YDATA(NDATA), the data values.
c
c    Input, double precision TVAL, the value where the spline is to
c    be evaluated.
c
c    Output, double precision YVAL, the value of the spline at TVAL.
c
      implicit none

      integer ndata

      double precision alpha
      double precision beta
      double precision d21
      double precision d32
      double precision d43
      integer left
      double precision mbasis(4,4)
      double precision mbasis_l(3,3)
      double precision mbasis_r(3,3)
      integer right
      double precision tdata(ndata)
      double precision tval
      double precision ydata(ndata)
      double precision yval
c
c  Check NDATA.
c
      if ( ndata .lt. 3 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'SPLINE_OVERHAUSER_NONUNI_VAL - Fatal error!'
        write ( *, '(a)' ) '  NDATA .lt. 3.'
        stop
      end if
c
c  Find the nearest interval [ TDATA(LEFT), TDATA(RIGHT) ] to TVAL.
c
      call r8vec_bracket ( ndata, tdata, tval, left, right )
c
c  Evaluate the spline in the given interval.
c
      if ( left .eq. 1 ) then

        d21 = sqrt ( ( tdata(2) - tdata(1) )**2 
     &             + ( ydata(2) - ydata(1) )**2 )

        d32 = sqrt ( ( tdata(3) - tdata(2) )**2 
     &             + ( ydata(3) - ydata(2) )**2 )

        alpha = d21 / ( d32 + d21 )

        call basis_matrix_overhauser_nul ( alpha, mbasis_l )

        call basis_matrix_tmp ( left, 3, mbasis_l, ndata, tdata, ydata, 
     &    tval, yval )

      else if ( left .lt. ndata-1 ) then

        d21 = sqrt ( ( tdata(left) - tdata(left-1) )**2 
     &             + ( ydata(left) - ydata(left-1) )**2 )

        d32 = sqrt ( ( tdata(left+1) - tdata(left) )**2 
     &             + ( ydata(left+1) - ydata(left) )**2 )

        d43 = sqrt ( ( tdata(left+2) - tdata(left+1) )**2 
     &             + ( ydata(left+2) - ydata(left+1) )**2 )

        alpha = d21 / ( d32 + d21 )
        beta  = d32 / ( d43 + d32 )

        call basis_matrix_overhauser_nonuni ( alpha, beta, mbasis )

        call basis_matrix_tmp ( left, 4, mbasis, ndata, tdata, ydata, 
     &    tval, yval )

      else if ( left .eq. ndata-1 ) then

        d32 = sqrt ( ( tdata(ndata-1) - tdata(ndata-2) )**2 
     &             + ( ydata(ndata-1) - ydata(ndata-2) )**2 )

        d43 = sqrt ( ( tdata(ndata) - tdata(ndata-1) )**2 
     &             + ( ydata(ndata) - ydata(ndata-1) )**2 )

        beta  = d32 / ( d43 + d32 )

        call basis_matrix_overhauser_nur ( beta, mbasis_r )

        call basis_matrix_tmp ( left, 3, mbasis_r, ndata, tdata, ydata, 
     &    tval, yval )

      else

        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'SPLINE_OVERHAUSER_NONUNI_VAL - Fatal error!'
        write ( *, '(a,i8)' ) '  Nonsensical value of LEFT = ', left
        write ( *, '(a,i8)' ) '  but 0 .lt. LEFT .lt. NDATA = ', ndata
        write ( *, '(a)' ) '  is required.'
        stop

      end if

      return
      end
      subroutine spline_overhauser_uni_val ( ndata, tdata, ydata, 
     &  tval, yval )

c*********************************************************************72
c
cc SPLINE_OVERHAUSER_UNI_VAL evaluates the uniform Overhauser spline.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    06 April 1999
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer NDATA, the number of data points.
c    NDATA must be at least 3.
c
c    Input, double precision TDATA(NDATA), the abscissas of the data points.
c    The values of TDATA are assumed to be distinct and increasing.
c    This routine also assumes that the values of TDATA are uniformly
c    spaced; for instance, TDATA(1) = 10, TDATA(2) = 11, TDATA(3) = 12...
c
c    Input, double precision YDATA(NDATA), the data values.
c
c    Input, double precision TVAL, the value where the spline is to
c    be evaluated.
c
c    Output, double precision YVAL, the value of the spline at TVAL.
c
      implicit none

      integer ndata

      integer left
      double precision mbasis(4,4)
      double precision mbasis_l(3,3)
      double precision mbasis_r(3,3)
      integer right
      double precision tdata(ndata)
      double precision tval
      double precision ydata(ndata)
      double precision yval
c
c  Check NDATA.
c
      if ( ndata .lt. 3 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'SPLINE_OVERHAUSER_UNI_VAL - Fatal error!'
        write ( *, '(a)' ) '  NDATA .lt. 3.'
        stop
      end if
c
c  Find the nearest interval [ TDATA(LEFT), TDATA(RIGHT) ] to TVAL.
c
      call r8vec_bracket ( ndata, tdata, tval, left, right )
c
c  Evaluate the spline in the given interval.
c
      if ( left .eq. 1 ) then

        call basis_matrix_overhauser_uni_l ( mbasis_l )

        call basis_matrix_tmp ( left, 3, mbasis_l, ndata, tdata, ydata, 
     &    tval, yval )

      else if ( left .lt. ndata-1 ) then

        call basis_matrix_overhauser_uni ( mbasis )

        call basis_matrix_tmp ( left, 4, mbasis, ndata, tdata, ydata, 
     &    tval, yval )

      else if ( left .eq. ndata-1 ) then

        call basis_matrix_overhauser_uni_r ( mbasis_r )

        call basis_matrix_tmp ( left, 3, mbasis_r, ndata, tdata, 
     &    ydata, tval, yval )

      end if

      return
      end
      subroutine spline_overhauser_val ( dim_num, ndata, tdata, ydata, 
     &  tval, yval )

c*********************************************************************72
c
cc SPLINE_OVERHAUSER_VAL evaluates an Overhauser spline.
c
c  Discussion:
c
c    Over the first and last intervals, the Overhauser spline is a
c    quadratic.  In the intermediate intervals, it is a piecewise cubic.
c    The Overhauser spline is also known as the Catmull-Rom spline.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c   08 April 1999
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    JA Brewer, DC Anderson,
c    Visual Interaction with Overhauser Curves and Surfaces,
c    SIGGRAPH 77,
c    in Proceedings of the 4th Annual Conference on Computer Graphics
c    and Interactive Techniques,
c    ASME, July 1977, pages 132-137.
c
c    Edwin Catmull, Raphael Rom,
c    A Class of Local Interpolating Splines,
c    in Computer Aided Geometric Design,
c    edited by Robert Barnhill, Richard Reisenfeld,
c    Academic Press, 1974, pages 317-326,
c    ISBN: 0120790505.
c
c    David Rogers, Alan Adams,
c    Mathematical Elements of Computer Graphics,
c    Second Edition,
c    McGraw Hill, 1989,
c    ISBN: 0070535299.
c
c  Parameters:
c
c    Input, integer DIM_NUM, the dimension of a single data point.
c    DIM_NUM must be at least 1.
c
c    Input, integer NDATA, the number of data points.
c    NDATA must be at least 3.
c
c    Input, double precision TDATA(NDATA), the abscissas of the data
c    points.  The values in TDATA must be in strictly ascending order.
c
c    Input, double precision YDATA(DIM_NUM,NDATA), the data points
c    corresponding to the abscissas.
c
c    Input, double precision TVAL, the abscissa value at which the spline
c    is to be evaluated.  Normally, TDATA(1) .le. TVAL .le. T(NDATA), and
c    the data will be interpolated.  For TVAL outside this range,
c    extrapolation will be used.
c
c    Output, double precision YVAL(DIM_NUM), the value of the spline at TVAL.
c
      implicit none

      integer ndata
      integer dim_num

      integer i
      integer left
      integer order
      integer right
      double precision tdata(ndata)
      double precision tval
      double precision ydata(dim_num,ndata)
      double precision yl(dim_num)
      double precision yr(dim_num)
      double precision yval(dim_num)
c
c  Check.
c
      call r8vec_order_type ( ndata, tdata, order )

      if ( order .ne. 2 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'SPLINE_OVERHAUSER_VAL - Fatal error!'
        write ( *, '(a)' ) 
     &    '  The data abscissas are not strictly ascending.'
        stop
      end if

      if ( ndata .lt. 3 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'SPLINE_OVERHAUSER_VAL - Fatal error!'
        write ( *, '(a)' ) '  NDATA .lt. 3.'
        stop
      end if

      if ( dim_num .lt. 1 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'SPLINE_OVERHAUSER_VAL - Fatal error!'
        write ( *, '(a)' ) '  DIM_NUM .lt. 1.'
        stop
      end if
c
c  Locate the abscissa interval T(LEFT), T(LEFT+1) nearest to or
c  containing TVAL.
c
      call r8vec_bracket ( ndata, tdata, tval, left, right )
c
c  Evaluate the "left hand" quadratic defined at T(LEFT-1), T(LEFT), T(RIGHT).
c
      if ( 0 .lt. left-1 ) then
        call parabola_val2 ( dim_num, ndata, tdata, ydata, left-1, 
     &    tval, yl )
      end if
c
c  Evaluate the "right hand" quadratic defined at T(LEFT), T(RIGHT), T(RIGHT+1).
c
      if ( right+1 .le. ndata ) then
        call parabola_val2 ( dim_num, ndata, tdata, ydata, left, 
     &    tval, yr )
      end if
c
c  Average the quadratics.
c
      if ( left .eq. 1 ) then

        do i = 1, dim_num
          yval(i) = yr(i)
        end do

      else if ( right .lt. ndata ) then

        do i = 1, dim_num
          yval(i) =  
     &      ( ( tdata(right) - tval               ) * yl(i)   
     &      + (                tval - tdata(left) ) * yr(i) ) 
     &      / ( tdata(right)        - tdata(left) )
        end do

      else

        do i = 1, dim_num
          yval(i) = yl(i)
        end do

      end if

      return
      end
      subroutine spline_pchip_set ( n, x, f, d )

c*********************************************************************72
c
cc SPLINE_PCHIP_SET sets derivatives for a piecewise cubic Hermite interpolant.
c
c  Discussion:
c
c    This routine computes what would normally be called a Hermite
c    interpolant.  However, the user is only required to supply function
c    values, not derivative values as well.  This routine computes
c    "suitable" derivative values, so that the resulting Hermite interpolant
c    has desirable shape and monotonicity properties.
c
c    The interpolant will have an extremum at each point where
c    monotonicity switches direction.
c
c    The resulting piecewise cubic Hermite function may be evaluated
c    by SPLINE_PCHIP_VAL.
c
c    This routine was originally named "PCHIM".
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    15 December 2008
c
c  Author:
c
c    Original FORTRAN77 version by Fred Fritsch.
c    FORTRAN90 version by John Burkardt.
c
c  Reference:
c
c    Fred Fritsch, Ralph Carlson,
c    Monotone Piecewise Cubic Interpolation,
c    SIAM Journal on Numerical Analysis,
c    Volume 17, Number 2, April 1980, pages 238-246.
c
c    Fred Fritsch, Judy Butland,
c    A Method for Constructing Local Monotone Piecewise Cubic Interpolants,
c    SIAM Journal on Scientific and Statistical Computing,
c    Volume 5, Number 2, 1984, pages 300-304.
c
c  Parameters:
c
c    Input, integer N, the number of data points.  N must be
c    at least 2.
c
c    Input, double precision X(N), the strictly increasing independent
c    variable values.
c
c    Input, double precision F(N), dependent variable values to be
c    interpolated.  F(I) is the value corresponding to X(I).
c    This routine is designed for monotonic data, but it will work for any
c    F array.  It will force extrema at points where monotonicity switches
c    direction.
c
c    Output, double precision D(N), the derivative values at the
c    data points.  If the data are monotonic, these values will determine
c    a monotone cubic Hermite function.
c
      implicit none

      integer n

      double precision d(n)
      double precision del1
      double precision del2
      double precision dmax
      double precision dmin
      double precision drat1
      double precision drat2
      double precision dsave
      double precision f(n)
      double precision h1
      double precision h2
      double precision hsum
      double precision hsumt3
      integer i
      integer ierr
      integer nless1
      double precision pchst
      double precision temp
      double precision w1
      double precision w2
      double precision x(n)
c
c  Check the arguments.
c
      if ( n .lt. 2 ) then
        ierr = -1
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'SPLINE_PCHIP_SET - Fatal error!'
        write ( *, '(a)' ) '  Number of data points less than 2.'
        stop
      end if

      do i = 2, n
        if ( x(i) .le. x(i-1) ) then
          ierr = -3
          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) 'SPLINE_PCHIP_SET - Fatal error!'
          write ( *, '(a)' ) '  X array not strictly increasing.'
          stop
        end if
      end do

      ierr = 0
      nless1 = n - 1
      h1 = x(2) - x(1)
      del1 = ( f(2) - f(1) ) / h1
      dsave = del1
c
c  Special case N=2, use linear interpolation.
c
      if ( n .eq. 2 ) then
        d(1) = del1
        d(n) = del1
        return
      end if
c
c  Normal case, 3 .le. N.
c
      h2 = x(3) - x(2)
      del2 = ( f(3) - f(2) ) / h2
c
c  Set D(1) via non-centered three point formula, adjusted to be
c  shape preserving.
c
      hsum = h1 + h2
      w1 = ( h1 + hsum ) / hsum
      w2 = -h1 / hsum
      d(1) = w1 * del1 + w2 * del2

      if ( pchst ( d(1), del1 ) .le. 0.0D+00 ) then

        d(1) = 0.0D+00
c
c  Need do this check only if monotonicity switches.
c
      else if ( pchst ( del1, del2 ) .lt. 0.0D+00 ) then

         dmax = 3.0D+00 * del1

         if ( abs ( dmax ) .lt. abs ( d(1) ) ) then
           d(1) = dmax
         end if

      end if
c
c  Loop through interior points.
c
      do i = 2, nless1

        if ( 2 .lt. i ) then
          h1 = h2
          h2 = x(i+1) - x(i)
          hsum = h1 + h2
          del1 = del2
          del2 = ( f(i+1) - f(i) ) / h2
        end if
c
c  Set D(I)=0 unless data are strictly monotonic.
c
        d(i) = 0.0D+00

        temp = pchst ( del1, del2 )

        if ( temp .lt. 0.0D+00 ) then

          ierr = ierr + 1
          dsave = del2
c
c  Count number of changes in direction of monotonicity.
c
        else if ( temp .eq. 0.0D+00 ) then

          if ( del2 .ne. 0.0D+00 ) then
            if ( pchst ( dsave, del2 ) .lt. 0.0D+00 ) then
              ierr = ierr + 1
            end if
            dsave = del2
          end if
c
c  Use Brodlie modification of Butland formula.
c
        else

          hsumt3 = 3.0D+00 * hsum
          w1 = ( hsum + h1 ) / hsumt3
          w2 = ( hsum + h2 ) / hsumt3
          dmax = max ( abs ( del1 ), abs ( del2 ) )
          dmin = min ( abs ( del1 ), abs ( del2 ) )
          drat1 = del1 / dmax
          drat2 = del2 / dmax
          d(i) = dmin / ( w1 * drat1 + w2 * drat2 )

        end if

      end do
c
c  Set D(N) via non-centered three point formula, adjusted to be
c  shape preserving.
c
      w1 = -h2 / hsum
      w2 = ( h2 + hsum ) / hsum
      d(n) = w1 * del1 + w2 * del2

      if ( pchst ( d(n), del2 ) .le. 0.0D+00 ) then
        d(n) = 0.0D+00
      else if ( pchst ( del1, del2 ) .lt. 0.0D+00 ) then
c
c  Need do this check only if monotonicity switches.
c
        dmax = 3.0D+00 * del2

        if ( abs ( dmax ) .lt. abs ( d(n) ) ) then
          d(n) = dmax
        end if

      end if

      return
      end
      subroutine spline_pchip_val ( n, x, f, d, ne, xe, fe )

c*********************************************************************72
c
cc SPLINE_PCHIP_VAL evaluates a piecewise cubic Hermite function.
c
c  Description:
c
c    This routine may be used by itself for Hermite interpolation, or as an
c    evaluator for SPLINE_PCHIP_SET.
c
c    This routine evaluates the cubic Hermite function at the points XE.
c
c    Most of the coding between the call to CHFEV and the end of
c    the IR loop could be eliminated if it were permissible to
c    assume that XE is ordered relative to X.
c
c    CHFEV does not assume that X1 is less than X2.  Thus, it would
c    be possible to write a version of SPLINE_PCHIP_VAL that assumes a strictly
c    decreasing X array by simply running the IR loop backwards
c    and reversing the order of appropriate tests.
c
c    The present code has a minor bug, which I have decided is not
c    worth the effort that would be required to fix it.
c    If XE contains points in [X(N-1),X(N)], followed by points less than
c    X(N-1), followed by points greater than X(N), the extrapolation points
c    will be counted (at least) twice in the total returned in IERR.
c
c    The evaluation will be most efficient if the elements of XE are
c    increasing relative to X; that is, for all J .le. K,
c      X(I) .le. XE(J)
c    implies
c      X(I) .le. XE(K).
c
c    If any of the XE are outside the interval [X(1),X(N)],
c    values are extrapolated from the nearest extreme cubic,
c    and a warning error is returned.
c
c    This routine was originally called "PCHFE".
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    14 August 2005
c
c  Author:
c
c    Original FORTRAN77 version by Fred Fritsch.
c    FORTRAN90 version by John Burkardt.
c
c  Reference:
c
c    Fred Fritsch, Ralph Carlson,
c    Monotone Piecewise Cubic Interpolation,
c    SIAM Journal on Numerical Analysis,
c    Volume 17, Number 2, April 1980, pages 238-246.
c
c  Parameters:
c
c    Input, integer N, the number of data points.  N must be
c    at least 2.
c
c    Input, double precision X(N), the strictly increasing independent
c    variable values.
c
c    Input, double precision F(N), the function values.
c
c    Input, double precision D(N), the derivative values.
c
c    Input, integer NE, the number of evaluation points.
c
c    Input, double precision XE(NE), points at which the function is to
c    be evaluated.
c
c    Output, double precision FE(NE), the values of the cubic Hermite
c    function at XE.
c
      implicit none

      integer n
      integer ne

      double precision d(n)
      double precision f(n)
      double precision fe(ne)
      integer i
      integer ierc
      integer ierr
      integer ir
      integer j
      integer j_first
      integer j_new
      integer j_save
      integer next(2)
      integer nj
      double precision x(n)
      double precision xe(ne)
c
c  Check arguments.
c
      if ( n .lt. 2 ) then
        ierr = -1
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'SPLINE_PCHIP_VAL - Fatal error!'
        write ( *, '(a)' ) '  Number of data points less than 2.'
        stop
      end if

      do i = 2, n
        if ( x(i) .le. x(i-1) ) then
          ierr = -3
          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) 'SPLINE_PCHIP_VAL - Fatal error!'
          write ( *, '(a)' ) '  X array not strictly increasing.'
          stop
        end if
      end do

      if ( ne .lt. 1 ) then
        ierr = -4
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'SPLINE_PCHIP_VAL - Fatal error!'
        write ( *, '(a)' ) '  Number of evaluation points less than 1.'
        return
      end if

      ierr = 0
c
c  Loop over intervals.
c  The interval index is IL = IR-1.
c  The interval is X(IL) .le. X .lt. X(IR).
c
      j_first = 1
      ir = 2

      do
c
c  Skip out of the loop if have processed all evaluation points.
c
        if ( ne .lt. j_first ) then
          exit
        end if
c
c  Locate all points in the interval.
c
        j_save = ne + 1

        do j = j_first, ne
          if ( x(ir) .le. xe(j) ) then
            j_save = j
            if ( ir .eq. n ) then
              j_save = ne + 1
            end if
            exit
          end if
        end do
c
c  Have located first point beyond interval.
c
        j = j_save

        nj = j - j_first
c
c  Skip evaluation if no points in interval.
c
        if ( nj .ne. 0 ) then
c
c  Evaluate cubic at XE(J_FIRST:J-1).
c
          call chfev ( x(ir-1), x(ir), f(ir-1), f(ir), d(ir-1), d(ir), 
     &      nj, xe(j_first:j-1), fe(j_first:j-1), next, ierc )

          if ( ierc .lt. 0 ) then
            ierr = -5
            write ( *, '(a)' ) ' '
            write ( *, '(a)' ) 'SPLINE_PCHIP_VAL - Fatal error!'
            write ( *, '(a)' ) '  Error return from CHFEV.'
            stop
          end if
c
c  In the current set of XE points, there are NEXT(2) to the right of X(IR).
c
          if ( next(2) .ne. 0 ) then

            if ( ir .lt. n ) then
              ierr = -5
              write ( *, '(a)' ) ' '
              write ( *, '(a)' ) 'SPLINE_PCHIP_VAL - Fatal error!'
              write ( *, '(a)' ) '  IR .lt. N.'
              stop
            end if
c
c  These are actually extrapolation points.
c
            ierr = ierr + next(2)

          end if
c
c  In the current set of XE points, there are NEXT(1) to the left of X(IR-1).
c
          if ( next(1) .ne. 0 ) then
c
c  These are actually extrapolation points.
c
            if ( ir .le. 2 ) then
              ierr = ierr + next(1)
            else

              j_new = -1

              do i = j_first, j-1
                if ( xe(i) .lt. x(ir-1) ) then
                  j_new = i
                  exit
                end if
              end do

              if ( j_new .eq. -1 ) then
                ierr = -5
                write ( *, '(a)' ) ' '
                write ( *, '(a)' ) 'SPLINE_PCHIP_VAL - Fatal error!'
                write ( *, '(a)' ) '  Could not bracket the data point.'
                stop
              end if
c
c  Reset J.  This will be the new J_FIRST.
c
              j = j_new
c
c  Now find out how far to back up in the X array.
c
              do i = 1, ir-1
                if ( xe(j) .lt. x(i) ) then
                  exit
                end if
              end do
c
c  At this point, either XE(J) .lt. X(1) or X(i-1) .le. XE(J) .lt. X(I) .
c
c  Reset IR, recognizing that it will be incremented before cycling.
c
              ir = max ( 1, i-1 )

            end if

          end if

          j_first = j

        end if

        ir = ir + 1

        if ( n .lt. ir ) then
          exit
        end if

      end do

      return
      end
      subroutine spline_quadratic_val ( ndata, tdata, ydata, tval, 
     &  yval, ypval )

c*********************************************************************72
c
cc SPLINE_QUADRATIC_VAL evaluates a piecewise quadratic spline at a point.
c
c  Discussion:
c
c    Because of the simple form of a piecewise quadratic spline,
c    the raw data points ( TDATA(I), YDATA(I)) can be used directly to
c    evaluate the spline at any point.  No processing of the data
c    is required.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    24 October 1999
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer NDATA, the number of data points defining
c    the spline.  NDATA should be odd and at least 3.
c
c    Input, double precision TDATA(NDATA), YDATA(NDATA), the values of
c    the independent and dependent variables at the data points.  The
c    values of TDATA should be distinct and increasing.
c
c    Input, double precision TVAL, the point at which the spline is to
c    be evaluated.
c
c    Output, double precision YVAL, YPVAL, the value of the spline and
c    its first derivative dYdT at TVAL.  YPVAL is not reliable if TVAL
c    is exactly equal to TDATA(I) for some I.
c
      implicit none

      integer ndata

      double precision dif1
      double precision dif2
      integer left
      integer right
      double precision t1
      double precision t2
      double precision t3
      double precision tdata(ndata)
      double precision tval
      double precision y1
      double precision y2
      double precision y3
      double precision ydata(ndata)
      double precision ypval
      double precision yval

      if ( ndata .lt. 3 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'SPLINE_QUADRATIC_VAL - Fatal error!'
        write ( *, '(a)' ) '  NDATA .lt. 3.'
        stop
      end if

      if ( mod ( ndata, 2 ) .eq. 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'SPLINE_QUADRATIC_VAL - Fatal error!'
        write ( *, '(a)' ) '  NDATA must be odd.'
        stop
      end if
c
c  Find the interval [ TDATA(LEFT), TDATA(RIGHT) ] that contains, or is
c  nearest to, TVAL.
c
      call r8vec_bracket ( ndata, tdata, tval, left, right )
c
c  Force LEFT to be odd.
c
      if ( mod ( left, 2 ) .eq. 0 ) then
        left = left - 1
      end if
c
c  Copy out the three abscissas.
c
      t1 = tdata(left)
      t2 = tdata(left+1)
      t3 = tdata(left+2)

      if ( t2 .le. t1 .or. t3 .le. t2 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'SPLINE_QUADRATIC_VAL - Fatal error!'
        write ( *, '(a)' ) '  T2 .le. T1 or T3 .le. T2.'
        stop
      end if
c
c  Construct and evaluate a parabolic interpolant for the data
c  in each dimension.
c
      y1 = ydata(left)
      y2 = ydata(left+1)
      y3 = ydata(left+2)

      dif1 = ( y2 - y1 ) / ( t2 - t1 )

      dif2 = ( ( y3 - y1 ) / ( t3 - t1 ) 
     &     - ( y2 - y1 ) / ( t2 - t1 ) ) / ( t3 - t2 )

      yval = y1 + ( tval - t1 ) * ( dif1 + ( tval - t2 ) * dif2 )
      ypval = dif1 + dif2 * ( 2.0D+00 * tval - t1 - t2 )

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
