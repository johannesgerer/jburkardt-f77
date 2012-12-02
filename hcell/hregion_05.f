      program main

!***********************************************************************
! 
!! HCELL solves the Navier-Stokes equations in an H-cell region.
!
!  Discussion:
!
!    The fluid flow problem is formulated in terms of
!    primitive variables - u,v, and p.
!
!    u_t - laplacian u + (u.grad)u + grad p = f
!                                div u = 0
!
!    Boundary conditions:  (u,v)=(0,0) on top
!                          (u,v)=0 on left, right and bottom
!
!    This version uses finite element techniques
!    with piecewise linear functions on triangles to approximate
!    the pressure and quadratics on triangles for the velocity
!    (Taylor-Hood element), isoparametric element
!
!    Thanks to Nick Chepurniy for pointing out a defective call to UBDRY 
!    in NSTOKE, 8 August 2004.
!
!  Input files:
!
!    UP000.TXT contains the initial values of the solution coefficients.
!
!  Modified:
!
!    14 July 2004
!
!  Local parameters:
!
!    double precision A(MINUN,MAXUN), storage for the banded system matrix.
!
!    double precision AREA(NELEMN), the area of each element.
!
!    integer BC_TYPE selects the boundary conditions, by controlling the value of ALPHA.
!    Legal values are 1, for a step function, 2 for a "hat" function, 3 for a sinusoid.
!
!    character ( len = * ) FILE_NAME is the name of the input file containing the
!    initial values of the solution, and is also used as a template for the
!    output files that store the solution at each timestep.
!
!    integer INDX(MAXND,NUK), lists the indices of the U, V, and P
!    variables associated with the node.  It would be useful if the 
!    index associated with pressure also indicated whether there was
!    no pressure variable associated with the node, or that it was
!    prescribed.  This could be done by assigning INDX(NODE,3) = 0
!    for the midside nodes of the 6 node quadratic elements.
!
!    integer MAXEL is an overestimate of the number of elements.
!
!    integer MAXND is an overestimate of the number of nodes.
!
!    integer MAXUN is an overestimate of the number of unknowns.
!
!    integer MINUN is an estimate of the necessary bandwidth of the system matrix.
!
!    integer MX counts the nuber of columns of nodes.
!
!    integer MY counts the number of rows of nodes.
!
!    integer NBAND, the bandwidth for the finite element matrix.
!
!    integer NELEMN, the number of elements.
!
!    integer NEQN, the total number of unknowns.
!
!    integer NLBAND, NUBAND, the lowwer and upper half bandwidths
!    for the finite element matrix.
!
!    integer NNODES, the number of nodes per element.
!
!    integer NODE(MAXEL,NNODES), the nodes that make up each element.
!
!    integer NP, the number of nodes.
!
!    integer NQUAD, the number of quadrature points used in assembly.
!    (This is 3)
!
!    integer NUK, the maximum number of unknowns associated with one node.
!    (This is 3)
!
!    integer NX counts, not quite all the elements in the X direction, but the number
!    of elements plus 1.
!
!    integer NX1, NX2, NX3, count the elements in the X direction in the three subregions,
!
!    integer NY counts, not quite all the elements in the Y direction, but the number
!    of elements plus 1.
!
!    integer NY1, NY2, NY3, count the elements in the Y direction in the three subregions.
!
!    integer REGION_DENSITY_X(3), REGION_DENSITY_Y(3), specifies the
!    density of elements in the two coordinate directions.
!
!    double precision REGION_X(4), REGION_Y(4), the coordinates of 
!    breakpoints that define 9 logical subrectangles.
!
!    double precision XC(NP), YC(NP), the coordinates of the nodes.
!
!    double precision XM(MAXEL,3), YM(MAXEL,3), the coordinates
!    of quadrature points in each element.
!
      implicit none
!
!  Set the "master parameters".
!
      integer nrefine
      parameter ( nrefine = 2 )

      integer nx1
      integer nx2
      integer nx3
      integer ny1
      integer ny2
      integer ny3

      parameter ( nx1 = 45 * nrefine )
      parameter ( nx2 = 15 * nrefine )
      parameter ( nx3 = 45 * nrefine )
      parameter ( ny1 = 5 * nrefine )
      parameter ( ny2 = 1 * nrefine )
      parameter ( ny3 = 5 * nrefine )
!
!  Set parameters that depend on the master parameters.
!
      integer maxel
      integer maxnd
      integer maxun
      integer minun
      integer mx
      integer my
      integer n_time
      integer nnodes
      integer nquad
      integer nuk
      integer nx
      integer ny

      parameter ( nx = nx1 + nx2 + nx3 + 1 )
      parameter ( ny = ny1 + ny2 + ny3 + 1 )
      parameter ( mx = 2*nx-1 )
      parameter ( my = 2*ny-1 )
      parameter ( maxel =  2*(nx-1)*(ny-1) )
      parameter ( maxnd =  mx*my )
      parameter ( maxun =  2*mx*my+nx*ny )
      parameter ( minun =  27*ny )
      parameter ( nnodes = 6 )
      parameter ( nuk = 3 )
      parameter ( nquad = 3 )
      parameter ( n_time = 100 )
!
      double precision a(minun,maxun)
      double precision alpha
      double precision area(maxel)
      integer bc_type
      double precision beta(6)
      double precision deltat
      double precision f(maxun)
      double precision g(maxun)
      integer i
      integer indx(maxnd,nuk)
      integer ip
      integer ipivot(maxun)
      integer iter
      integer iu
      integer iv
      integer mrow1
      integer nband
      integer ncol1
      integer nelemn
      integer neqn
      integer nlband
      integer node(maxel,nnodes)
      character*80 node_file_name
      logical node_mask(maxnd)
      integer np
      integer nrow1
      integer nsim
      integer nsteps
      integer nuband
      character*20 p_file_name
      double precision pi
      double precision rdel
      integer region_density_x(3)
      integer region_density_y(3)
      double precision region_x(4)
      double precision region_y(4)
      character*80 tec_file_name
      character*80 title
      double precision tolns
      double precision tolopt
      double precision ubdry
      double precision uold(maxun)
      character*20 uv_file_name
      double precision uvel
      double precision visc_inv
      double precision vvel
      double precision xc(maxnd)
      double precision xlngth
      double precision xm(maxel,nquad)
      double precision yc(maxnd)
      double precision ylngth
      double precision ym(maxel,nquad)
!
!  Initialize some variables.
!
      bc_type = 1

      beta(1) = 1.0d0
      beta(2) = 1.0d0
      beta(3) = 1.0d0
      beta(4) = 1.0d0
      beta(5) = 1.0d0
      beta(6) = 1.0d0

      mrow1 = minun
      nsim = 3
      nsteps = 15
      p_file_name = 'p000.txt'
      pi = 4.d0 * datan ( 1.0d0 )

      region_density_x(1) = nx1 
      region_density_x(2) = nx2 
      region_density_x(3) = nx3 
      region_density_y(1) = ny1 
      region_density_y(2) = ny2 
      region_density_y(3) = ny3 

      region_x(1) =   0.0D+00
      region_x(2) =  45.0D+00
      region_x(3) =  60.0D+00
      region_x(4) = 105.0D+00

      region_y(1) =  0.0D+00
      region_y(2) =  5.0D+00
      region_y(3) =  6.0D+00
      region_y(4) = 11.0D+00

      tec_file_name = 'uv000.tec'
      tolns = 1.0d-6
      tolopt = 1.0d-6
      uv_file_name = 'uv000.txt'
      visc_inv = 100.0D+00

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'HREGION Version 5:'
      write ( *, '(a)' ) '  Solve the Navier Stokes fluid flow'
      write ( *, '(a)' ) '  equations in an H-shaped region,'
      write ( *, '(a)' ) '  using finite elements.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Modifications:'
      write ( *, '(a)' ) '    08 August 2004: '
      write ( *, '(a)' ) '      Corrected defective UBDRY call in NSTOKE.'
      write ( *, '(a)' ) '    25 May 2004: '
      write ( *, '(a)' ) '      Rename inverse viscosity VISC_INV.'
      write ( *, '(a)' ) '    22 May 2004: NREFINE is now a'
      write ( *, '(a)' ) '                 parameter, and seems OK.'
      write ( *, '(a)' ) '    21 May 2004: IMPLICIT NONE everywhere;'
      write ( *, '(a)' ) '                 Corrections to UBDRY.'
      write ( *, '(a)' ) '                 Correction to HCELL_DOF_SET.'
      write ( *, '(a)' ) '    20 May 2004: added TECPLOT output.'
      write ( *, '(a)' ) '    18 May 2004: correction to HCELL_DOF_SET.'
      write ( *, '(a)' ) ' '
      write ( *, '(a,i6)' ) '  Refinement factor = ', nrefine
      write ( *, '(a,g14.6)' ) '  Inverse viscosity = ', visc_inv
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  The X direction is divided into three'
      write ( *, '(a)' ) '  regions, with element densities:'
      write ( *, '(a)' ) ' '
      write ( *, '(4x,3i6)' ) ( region_density_x(i), i = 1, 3 )
      write ( *, '(a)' ) ' '
      write ( *, '(a,i6)' ) '  Corresponding NX = ', nx
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  The Y direction is divided into three'
      write ( *, '(a)' ) '  regions, with element densities:'
      write ( *, '(a)' ) ' '
      write ( *, '(4x,3i6)' ) ( region_density_y(i), i = 1, 3 )
      write ( *, '(a)' ) ' '
      write ( *, '(a,i6)' ) '  Corresponding NY = ', ny
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  The X subregions are demarked by 4 values:'
      write ( *, '(a)' ) ' '
      write ( *, '(4x,4f10.4)' ) ( region_x(i), i = 1, 4 )
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  The Y subregions are demarked by 4 values:'
      write ( *, '(a)' ) ' '
      write ( *, '(4x,4f10.4)' ) ( region_y(i), i = 1, 4 )

      write ( *, '(a)' ) ' '
      write ( *, '(a,i6)' ) '  Maximum number of nodes =    ', maxnd
      write ( *, '(a,i6)' ) '  Maximum number of elements = ', maxel
      write ( *, '(a,i6)' ) '  Maximum number of unknowns = ', maxun
      write ( *, '(a,i6)' ) '  Maximum matrix dimension 1 = ', minun
!
!  SETGRD constructs grid, numbers unknowns, calculates areas,
!  and points for midpoint quadrature rule, bandwidth and neqn
!
      call setgrd ( xc, yc, area, xm, ym, region_density_x,
     &   region_density_y, region_x, region_y,
     &   node, indx, nlband, nuband, nband, nelemn, np, 
     &   nnodes, nuk, nquad, neqn, maxnd, maxel )

      write ( *, '(a)' ) ' '
      write ( *, '(a,i6)' ) '  Number of nodes =    ', np
      write ( *, '(a,i6)' ) '  Number of elements = ', nelemn
      write ( *, '(a,i6)' ) '  Number of unknowns = ', neqn
!
!  Write the coordinates of the nodes of 6-node triangles to a file.
!
      call xy6_write ( 'xy6.txt', np, xc, yc )
!
!  Write the coordinates of the nodes of 3-node triangles to a file.
!
      call xy3_write ( 'xy3.txt', maxnd, np, indx, xc, yc )
!
!  Write the element node matrix to a file.
!
      if ( .false. ) then
        call element_node_write ( 'element_node.txt', maxel, nelemn,
     &    node )
      end if
!
!  Make a plot of the nodes.
!
      if ( .false. ) then

        do i = 1, np
          if ( i <= 100 ) then
            node_mask(i) = .true.
          else
            node_mask(i) = .false.
          end if
        end do

        title = 'H Cell Nodes'
        node_file_name = 'hcell_nodes.eps'

        call node_eps ( node_file_name, np, node_mask, xc, yc, title )

      end if

      nrow1 = nlband + nlband + nuband + 1
      ncol1 = neqn

      do i = 1, neqn
        f(i) = 0.d0
      end do
!
!  Read the initial value of the solution from a file.
!
!      if ( .false. ) then

         call uv_read ( uv_file_name, neqn, uold )
!         call p_read ( p_file_name, neqn, pold )

        deltat = 0.0002d0
        rdel = 1.d0 / deltat
        write ( *, * ) 'delta t = ', deltat
!      else
!
!  Or assume we are doing a steady state calculation.
!
!        do i = 1, neqn
!          uold(i) = 0.0D+00
!        end do
!
!        deltat = 0.0D+00 
!        rdel = 0.0D+00
!        write(*,*) " computing a steady state problem "

!      end if
!
!  Carry out the time iteration.
!
      do iter = 1, n_time

         write ( *, * )
         write ( *, * ) 'Time step ', iter
        
        if ( bc_type == 1 ) then

          if ( iter .le. 250 ) then
            alpha = 5.d0
          else 
            alpha = 1.d0
          end if

        else if ( bc_type == 2 ) then

          if ( iter .le. 250 ) then
            alpha = 80.d0*dble(iter)*deltat+1.d0
          else 
            alpha = -80.d0*dble(iter)*deltat+9.d0
          end if

        else if ( bc_type == 3 ) then

          alpha = 2.d0*sin(dble(iter)*0.01*pi)   

        end if

        do i = 1, neqn
          g(i) = f(i)
        end do

        do i = 1, neqn
          f(i) = 0.d0
        end do

        call nstoke ( xc, yc, area, xm, ym,
     &     a, f, g, uold, visc_inv, tolns, xlngth, ylngth,
     &     node, indx, ipivot, mrow1,
     &     nlband, nuband, nband, nrow1, ncol1,
     &     nelemn, np, nnodes, nuk, nquad, neqn,
     &     nsteps, nsim, maxnd, maxel, rdel, alpha, beta )
!
!  Save u=(gx,gy) to 'ue.dat' for 't.f'
!
        do i = 1, neqn
          uold(i) = f(i) 
        end do
!
!  Increment the file name, and 
!  save u=(gx,gy) to 'up???.dat'
!
         call file_name_inc ( uv_file_name )

         if ( iter .eq. n_time ) then

           call uv_write ( uv_file_name, maxnd, np, neqn, f, indx, 
     &        xc, yc, beta )

         end if
 
         call file_name_inc ( p_file_name )

         if ( iter .eq. n_time ) then

           call p_write ( p_file_name, maxnd, np, neqn, f, indx, 
     &        xc, yc )

        end if

         call file_name_inc ( tec_file_name )

         if ( iter .eq. n_time ) then

           call uv_write_tec ( tec_file_name, maxnd, np, 
     &        xc, yc, maxel, nelemn, node, 
     &        indx, neqn, f, beta )

        end if
     
      end do

      open ( unit = 2, file = 'data.txt' )

      write(2,*) " time dependent results for chemical attack"

      do ip = 1, np

        iu = indx(ip,1)

        if ( iu < 0 ) then
          uvel = ubdry ( 1, ip, iu, beta, xc, yc )
        else if ( iu == 0 ) then
          uvel = 0.0
        else if ( 0 < iu ) then
          uvel = f(iu)
        end if

        iv = indx(ip,2)

        if ( iv < 0 ) then
          vvel = ubdry ( 2, ip, iv, beta, xc, yc )
        else if ( iv == 0 ) then
          vvel = 0.0
        else if ( 0 < iv ) then
          vvel = f(iv)
        end if

        write(2,1212) ip, xc(ip), yc(ip), uvel, vvel

      end do 

1212  format("final values: node, x,y,u,v", i5, 2f12.5,10x,2e12.5)
      close ( unit = 2 )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'HCELL:'
      write ( *, '(a)' ) '  Normal end of execution.'

      stop
      end
      function ch_is_digit ( c )
!
!***********************************************************************
!
!! CH_IS_DIGIT returns .TRUE. if a character is a decimal digit.
!
!
!  Modified:
!
!    15 January 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character C, the character to be analyzed.
!
!    Output, logical CH_IS_DIGIT, .TRUE. if C is a digit, .FALSE. otherwise.
!
      implicit none
!
      character c
      logical ch_is_digit
!
      if ( lge ( c, '0' ) .and. lle ( c, '9' ) ) then
        ch_is_digit = .true.
      else
        ch_is_digit = .false.
      end if

      return
      end
      subroutine ch_to_digit ( c, digit )
!
!***********************************************************************
!
!! CH_TO_DIGIT returns the integer value of a base 10 digit.
!
!
!  Example:
!
!     C   DIGIT
!    ---  -----
!    '0'    0
!    '1'    1
!    ...  ...
!    '9'    9
!    ' '    0
!    'X'   -1
!
!  Modified:
!
!    04 August 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character C, the decimal digit, '0' through '9' or blank
!    are legal.
!
!    Output, integer DIGIT, the corresponding integer value.  If C was
!    'illegal', then DIGIT is -1.
!
      implicit none
!
      character c
      integer digit
!
      if ( lge ( c, '0' ) .and. lle ( c, '9' ) ) then

        digit = ichar ( c ) - 48

      else if ( c .eq. ' ' ) then

        digit = 0

      else

        digit = -1

      end if

      return
      end
      function d_huge ( )
!
!***********************************************************************
!
!! D_HUGE returns a "huge" double precision number.
!
!
!  Modified:
!
!    13 April 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, double precision D_HUGE, a huge number.
!
      implicit none
!
      double precision d_huge
!
      d_huge = 1.0D+30

      return
      end
      subroutine daxpy ( n, da, dx, incx, dy, incy )
!
!***********************************************************************
!
!! DAXPY adds a multiple of one vector to another.
!
!
!  Author:
!
!    Jack Dongarra
!
!  Parameters:
!
!    Input, integer N, the number of entries in the vectors.
!
!    Input, double precision DA, the multiplier.
!
!    Input, double precision DX(*), the first vector.
!
!    Input, integer INCX, the increment between entries in the first vector.
!
!    Input/output, double precision DY(*), the second vector, to which
!    a multiple of the first vector has been added on output.
!
!    Input, integer INCY, the increment between entries in the second vector.
!
      implicit none
!
      double precision da
      double precision dx(*)
      double precision dy(*)
      integer i
      integer incx
      integer incy
      integer ix
      integer iy
      integer m
      integer n
!
      if ( n .le. 0 ) then
        return
      end if

      if  ( da .eq. 0.0D0 ) then
        return
      end if

      if ( incx .eq. 1 .and. incy .eq. 1 ) then

        m = mod ( n, 4 )

        do i = 1, m
          dy(i) = dy(i) + da * dx(i)
        end do

        do i = m+1, n, 4
          dy(i)   = dy(i)   + da * dx(i)
          dy(i+1) = dy(i+1) + da * dx(i+1)
          dy(i+2) = dy(i+2) + da * dx(i+2)
          dy(i+3) = dy(i+3) + da * dx(i+3)
        end do

      else

        if ( incx .lt. 0 ) then
          ix = (-n+1) * incx + 1
        else
          ix = 1
        end if

        if ( incy .lt. 0 ) then
          iy = (-n+1) * incy + 1
        else
          iy = 1
        end if

        do i = 1, n
          dy(iy) = dy(iy) + da * dx(ix)
          ix = ix + incx
          iy = iy + incy
        end do

      end if

      return
      end
      function ddot ( n, dx, incx, dy, incy )
!
!***********************************************************************
!
!! DDOT forms the dot product of two vectors.
!
!
!  Author:
!
!    Jack Dongarra
!
!  Parameters:
!
!    Input, integer N, the number of entries in the vectors.
!
!    Input, double precision DX(*), the first vector.
!
!    Input, integer INCX, the increment between entries in the first vector.
!
!    Input, double precision DY(*), the second vector.
!
!    Input, integer INCY, the increment between entries in the second vector.
!
!    Output, double precision DDOT, the dot product of the two vectors.
!    
      implicit none
!
      double precision ddot
      double precision dtemp
      double precision dx(*)
      double precision dy(*)
      integer i
      integer incx
      integer incy
      integer ix
      integer iy
      integer m
      integer n
!
      ddot = 0.0D0
      dtemp = 0.0D0

      if ( n .le. 0 ) then
        return
      end if

      if ( incx .eq. 1 .and. incy .eq. 1 ) then

        m = mod ( n, 5 )

        do i = 1, m
          dtemp = dtemp + dx(i) * dy(i)
        end do

        do i = m+1, n, 5
          dtemp = dtemp 
     &       + dx(i)   * dy(i) 
     &       + dx(i+1) * dy(i+1) 
     &       + dx(i+2) * dy(i+2) 
     &       + dx(i+3) * dy(i+3) 
     &       + dx(i+4) * dy(i+4)
        end do

      else

        if ( incx .lt. 0 ) then
          ix = ( -n + 1 ) * incx + 1
        else
          ix = 1
        end if

        if ( incy .lt. 0 ) then
          iy = ( -n + 1 ) * incy + 1
        else
          iy = 1
        end if

        do i = 1, n
          dtemp = dtemp + dx(ix) * dy(iy)
          ix = ix + incx
          iy = iy + incy
        end do

      end if

      ddot = dtemp

      return
      end
      subroutine dgbfa ( abd, lda, n, ml, mu, ipvt, info )
!
!***********************************************************************
!
!! DGBFA factors a double precision band matrix by elimination.
!
!
!  Discussion:
!
!    If A is a band matrix, the following program segment
!    will set up the input.
!
!      ml = (band width below the diagonal)
!      mu = (band width above the diagonal)
!      m = ml + mu + 1
!      do j = 1, n
!        i1 = max(1, j-mu)
!        i2 = min(n, j+ml)
!        do i = i1, i2
!          k = i - j + m
!          abd(k,j) = a(i,j)
!        end do
!      end do
!
!    This uses rows  ml+1  through  2*ml+mu+1  of  abd .
!    in addition, the first  ml  rows in  abd  are used for
!    elements generated during the triangularization.
!    the total number of rows needed in  abd  is  2*ml+mu+1 .
!    the  ml+mu by ml+mu  upper left triangle and the
!    ml by ml  lower right triangle are not referenced.
!
!  Author:
!
!    Cleve Moler
!
!     on entry
!
!        abd     double precision(lda, n)
!                contains the matrix in band storage.  the columns
!                of the matrix are stored in the columns of  abd  and
!                the diagonals of the matrix are stored in rows
!                ml+1 through 2*ml+mu+1 of  abd .
!                see the comments below for details.
!
!        lda     integer
!                the leading dimension of the array  abd .
!                lda must be .ge. 2*ml + mu + 1 .
!
!        n       integer
!                the order of the original matrix.
!
!        ml      integer
!                number of diagonals below the main diagonal.
!                0 .le. ml .lt. n .
!
!        mu      integer
!                number of diagonals above the main diagonal.
!                0 .le. mu .lt. n .
!                more efficient if  ml .le. mu .
!     on return
!
!        abd     an upper triangular matrix in band storage and
!                the multipliers which were used to obtain it.
!                the factorization can be written  a = l*u  where
!                l  is a product of permutation and unit lower
!                triangular matrices and  u  is upper triangular.
!
!        ipvt    integer(n)
!                an integer vector of pivot indices.
!
!        info    integer
!                = 0  normal value.
!                = k  if  u(k,k) .eq. 0.0 .  this is not an error
!                     condition for this subroutine, but it does
!                     indicate that dgbsl will divide by zero if
!                     called.  use  rcond  in dgbco for a reliable
!                     indication of singularity.
!

      implicit none

      integer lda
      integer n

      double precision abd(lda,n)
      integer i
      integer idamax
      integer i0
      integer info
      integer ipvt(n)
      integer j
      integer ju
      integer jz
      integer j0
      integer j1
      integer k
      integer l
      integer lm
      integer m
      integer ml
      integer mm
      integer mu
      double precision t
!
      m = ml + mu + 1
      info = 0
!
!     zero initial fill-in columns
!
      j0 = mu + 2
      j1 = min ( n, m ) - 1

      do jz = j0, j1
        i0 = m + 1 - jz
        do i = i0, ml
          abd(i,jz) = 0.0d0
        end do
      end do

      jz = j1
      ju = 0
!
!     gaussian elimination with partial pivoting
!
      do k = 1, n-1
!
!        zero next fill-in column
!
         jz = jz + 1
         if ( jz <= n ) then
           do i = 1, ml
             abd(i,jz) = 0.0d0
           end do
         end if
!
!        find l = pivot index
!
         lm = min ( ml, n-k )
         l = idamax(lm+1,abd(m,k),1) + m - 1
         ipvt(k) = l + k - m
!
!        zero pivot implies this column already triangularized
!
         if (abd(l,k) .eq. 0.0d0) go to 100
!
!           interchange if necessary
!
            if (l .ne. m) then
              t = abd(l,k)
              abd(l,k) = abd(m,k)
              abd(m,k) = t
            end if
!
!  compute multipliers
!
            t = -1.0d0/abd(m,k)
            call dscal(lm,t,abd(m+1,k),1)
!
!  row elimination with column indexing
!
            ju = min ( max (ju,mu+ipvt(k)),n)
            mm = m

            do j = k+1, ju
               l = l - 1
               mm = mm - 1
               t = abd(l,j)
               if ( l .ne. mm ) then
                  abd(l,j) = abd(mm,j)
                  abd(mm,j) = t
               end if
               call daxpy(lm,t,abd(m+1,k),1,abd(mm+1,j),1)
            end do

         go to 110
  100    continue
            info = k
  110    continue

      end do

      ipvt(n) = n
      if (abd(m,n) .eq. 0.0d0) info = n

      return
      end
      subroutine dgbsl ( abd, lda, n, ml, mu, ipvt, b, job )
!
!***********************************************************************
!
!! DGBSL solves a double precision band linear system.
!
!
!  Discussion:
!
!    The linear system has the form
!
!      A * x = b  or  A' * x = b
!
!    The matrix has been factored by dgbco or dgbfa.
!
!    A division by zero will occur if the input factor contains a
!    zero on the diagonal.  technically this indicates singularity
!    but it is often caused by improper arguments or improper
!    setting of lda .  it will not occur if the subroutines are
!    called correctly and if dgbco has set rcond .gt. 0.0
!    or dgbfa has set info .eq. 0 .
!
!  Author:
!
!    Cleve Moler
!
!  Parameters:
!
!        abd     double precision(lda, n)
!                the output from dgbco or dgbfa.
!
!        lda     integer
!                the leading dimension of the array  abd .
!
!        n       integer
!                the order of the original matrix.
!
!        ml      integer
!                number of diagonals below the main diagonal.
!
!        mu      integer
!                number of diagonals above the main diagonal.
!
!        ipvt    integer(n)
!                the pivot vector from dgbco or dgbfa.
!
!        b       double precision(n)
!                the right hand side vector.
!
!        job     integer
!                = 0         to solve  a*x = b ,
!                = nonzero   to solve  trans(a)*x = b , where
!                            trans(a)  is the transpose.
!
!     on return
!
!        b       the solution vector  x .
!
      implicit none

      integer lda

      double precision abd(lda,*)
      double precision b(*)
      double precision ddot
      integer ipvt(*)
      integer job
      integer k
      integer kb
      integer l
      integer la
      integer lb
      integer lm
      integer m
      integer ml
      integer mu
      integer n
      double precision t
!
      m = mu + ml + 1

      if ( job .eq. 0 ) then

        if ( ml .ne. 0 ) then

          do k = 1, n - 1
            lm = min(ml,n-k)
            l = ipvt(k)
            t = b(l)
            if ( l .ne. k ) then
              b(l) = b(k)
              b(k) = t
            end if
            call daxpy(lm,t,abd(m+1,k),1,b(k+1),1)
          end do

        end if

        do kb = 1, n
          k = n + 1 - kb
          b(k) = b(k)/abd(m,k)
          lm = min(k,m) - 1
          la = m - lm
          lb = k - lm
          t = -b(k)
          call daxpy(lm,t,abd(la,k),1,b(lb),1)
        end do

      else

        do k = 1, n
          lm = min(k,m) - 1
          la = m - lm
          lb = k - lm
          t = ddot(lm,abd(la,k),1,b(lb),1)
          b(k) = (b(k) - t)/abd(m,k)
        end do

        if ( ml .ne. 0 ) then

          do kb = 1, n - 1
            k = n - kb
            lm = min (ml,n-k)
            b(k) = b(k) + ddot(lm,abd(m+1,k),1,b(k+1),1)
            l = ipvt(k)
            if ( l .ne. k ) then
              t = b(l)
              b(l) = b(k)
              b(k) = t
            end if
          end do

        end if

      end if 

      return
      end
      subroutine digit_inc ( c )
!
!***********************************************************************
!
!! DIGIT_INC increments a decimal digit.
!
!
!  Example:
!
!    Input  Output
!    -----  ------
!    '0'    '1'
!    '1'    '2'
!    ...
!    '8'    '9'
!    '9'    '0'
!    'A'    'A'
!
!  Modified:
!
!    04 August 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, character C, a digit to be incremented.
!
      implicit none
!
      character c
      integer digit

      call ch_to_digit ( c, digit )

      if ( digit .eq. -1 ) then
        return
      end if

      digit = digit + 1

      if ( digit .eq. 10 ) then
        digit = 0
      end if

      call digit_to_ch ( digit, c )

      return
      end
      subroutine digit_to_ch ( digit, c )
!
!***********************************************************************
!
!! DIGIT_TO_CH returns the character representation of a decimal digit.
!
!
!  Example:
!
!    DIGIT   C
!    -----  ---
!      0    '0'
!      1    '1'
!    ...    ...
!      9    '9'
!     17    '*'
!
!  Modified:
!
!    04 August 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer DIGIT, the digit value between 0 and 9.
!
!    Output, character C, the corresponding character, or '*' if DIGIT
!    was illegal.
!
      implicit none
!
      character c
      integer digit
!
      if ( 0 .le. digit .and. digit .le. 9 ) then

        c = char ( digit + 48 )

      else

        c = '*'

      end if

      return
      end
      subroutine dscal ( n, da, dx, incx )
!
!***********************************************************************
!
!! DSCAL scales a vector by a constant.
!
!
!  Author:
!
!    Jack Dongarra
!
!  Parameters:
!
!    Input, integer N, the number of entries in the vector.
!
!    Input, double precision DA, the scale factor.
!
!    Input/output, double precision DX(*), the vector to be scaled.
!
!    Input, integer INCX, the increment between successive entries.
!
      implicit none
!
      double precision da
      double precision dx(*)
      integer i
      integer incx
      integer m
      integer n
      integer nincx
!
      if ( n .le. 0 ) then
        return
      end if

      if ( incx .eq. 1 ) then

        m = mod ( n, 5 )

        do i = 1, m
          dx(i) = da * dx(i)
        end do

        do i = m+1, n, 5
          dx(i)   = da * dx(i)
          dx(i+1) = da * dx(i+1)
          dx(i+2) = da * dx(i+2)
          dx(i+3) = da * dx(i+3)
          dx(i+4) = da * dx(i+4)
        end do

      else

        nincx = n * incx
        do i = 1, nincx, incx
          dx(i) = da * dx(i)
        end do

      end if

      return
      end
      subroutine element_node_bandwidth ( maxnd, maxel, nnodes, 
     &   nelemn, node, neqn, np, indx, nlband, nuband, nband )
!
!***********************************************************************
!
!! ELEMENT_NODE_BANDWIDTH determines the bandwidth associated with the grid.
!
!
!  Modified:
!
!    31 March 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer MAXND, the maximum number of nodes
!
!    Input, integer MAXEL, the maximum number of elements.
!
!    Input, integer NNODES, the order of the elements.
!
!    Input, integer NELEMN, the number of elements.
!
!    Input, integer NODE(MAXEL,NNODES), the nodes in each element.
!
!    Input, integer NEQN, the number of degrees of freedom.
!
!    Input, integer NP, the number of nodes.
!
!    Input, integer INDX(MAXND,3), the nodal degrees of freedom.
!
!    Output, integer NLBAND, NUBAND, NBAND, the lower, upper and total bandwidths.
!
      implicit none
!
      integer maxnd
      integer maxel
      integer nelemn
      integer neqn
      integer nnodes
      integer np
!
      integer dof_max(neqn)
      integer dof_min(neqn)
      integer element
      integer i
      integer ieqn
      integer indx(maxnd,3)
      integer l1
      integer l2
      integer n1
      integer n2
      integer nband
      integer nlband
      integer node(maxel,nnodes)
      integer nuband
      integer p1
      integer p2
      integer u1
      integer u2
      integer v1
      integer v2
!
      do ieqn = 1, neqn
        dof_min(ieqn) = ieqn
        dof_max(ieqn) = ieqn
      end do

      do element = 1, nelemn

        do l1 = 1, nnodes

          n1 = node(element,l1)

          u1 = indx(n1,1)
          v1 = indx(n1,2)
          p1 = indx(n1,3)

          do l2 = 1, nnodes

            n2 = node(element,l2)

            u2 = indx(n2,1)
            v2 = indx(n2,2)
            p2 = indx(n2,3)

            if ( 1 <= u1 .and. 1 <= u2 ) then
              dof_min(u1) = min ( dof_min(u1), u2 )
              dof_max(u1) = max ( dof_max(u1), u2 )
            end if

            if ( 1 <= u1 .and. 1 <= v2 ) then
              dof_min(u1) = min ( dof_min(u1), v2 )
              dof_max(u1) = max ( dof_max(u1), v2 )
            end if

            if ( 1 <= u1 .and. 1 <= p2 ) then
              dof_min(u1) = min ( dof_min(u1), p2 )
              dof_max(u1) = max ( dof_max(u1), p2 )
            end if

            if ( 1 <= v1 .and. 1 <= u2 ) then
              dof_min(v1) = min ( dof_min(v1), u2 )
              dof_max(v1) = max ( dof_max(v1), u2 )
            end if

            if ( 1 <= v1 .and. 1 <= v2 ) then
              dof_min(v1) = min ( dof_min(v1), v2 )
              dof_max(v1) = max ( dof_max(v1), v2 )
            end if

            if ( 1 <= v1 .and. 1 <= p2 ) then
              dof_min(v1) = min ( dof_min(v1), p2 )
              dof_max(v1) = max ( dof_max(v1), p2 )
            end if

            if ( 1 <= p1 .and. 1 <= u2 ) then
              dof_min(p1) = min ( dof_min(p1), u2 )
              dof_max(p1) = max ( dof_max(p1), u2 )
            end if

            if ( 1 <= p1 .and. 1 <= v2 ) then
              dof_min(p1) = min ( dof_min(p1), v2 )
              dof_max(p1) = max ( dof_max(p1), v2 )
            end if

          end do
        end do
      end do

      nlband = 0
      nuband = 0
      do ieqn = 1, neqn
        nlband = max ( nlband, ieqn - dof_min(ieqn) )
        nuband = max ( nuband, dof_max(ieqn) - ieqn )
      end do

      nband = nlband + nuband + 1

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'ELEMENT_NODE_BANDWIDTH:'
      write ( *, '(a,i6)' ) '  Lower half bandwidth = ', nlband
      write ( *, '(a,i6)' ) '  Upper half bandwidth = ', nuband
      write ( *, '(a,i6)' ) '  Total bandwidth =      ', nband
    
      return
      end
      subroutine element_node_write ( file_name, element_max, 
     &  element_num, element_node )
!
!***********************************************************************
!
!! ELEMENT_NODE_WRITE writes the element node list to a file.
!
!
!  Modified:
!
!    21 April 2004
!
!  Parameters:
!
!    Input, character ( len = * ) FILE_NAME, the file to be written.
!
!    Input, integer ELEMENT_MAX, the maximum number of elements.
!
!    Input, integer ELEMENT_NUM, the number of elements.
!
!    Input, integer ELEMENT_NODE(ELEMENT_MAX,6), the nodes that make up each element.
!
      implicit none
!
      integer element_max
!
      integer element
      integer element_node(element_max,6)
      integer element_num
      character ( len = * ) file_name
      integer j
!
      open ( unit = 1, file = file_name, err = 10 )

      do element = 1, element_num
        write ( 1, * ) ( element_node(element,j), j = 1, 6 )
      end do

      close ( unit = 1 )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'ELEMENT_NODE_WRITE:'
      write ( *, '(a)' ) '  The 3-node triangle nodal coordinates were'
      write ( *, '(a)' ) '  written to the output file:'
      write ( *, '(a)' ) ' '
      write ( *, '(4x,a)' ) file_name
      return

10    continue

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'ELEMENT_NODE_WRITE - Fatal error!'
      write ( *, '(a)' ) '  Could not open the output file:'
      write ( *, '(a)' ) ' '
      write ( *, '(4x,a)' ) file_name

      stop
      end
      subroutine file_name_inc ( file_name )
!
!***********************************************************************
!
!! FILE_NAME_INC generates the next file name in a series.
!
!
!  Discussion:
!
!    It is assumed that the digits in the name, whether scattered or
!    connected, represent a number that is to be increased by 1 on
!    each call.  If this number is all 9's on input, the output number
!    is all 0's.  Non-numeric letters of the name are unaffected, and
!    if the name contains no digits, then nothing is done.
!
!  Examples:
!
!      Input          Output
!      -----          ------
!      a7to11.txt     a7to12.txt
!      a7to99.txt     a8to00.txt
!      a9to99.txt     a0to00.txt
!      cat.txt        cat.txt
!
!  Modified:
!
!    09 August 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, character ( len = * ) FILE_NAME.
!    On input, a character string to be incremented.
!    On output, the incremented string.
!
      implicit none
!
      character c
      logical ch_is_digit
      character*(*) file_name
      integer i
      integer lens
!
      lens = len ( file_name )

      do i = lens, 1, -1

        c = file_name(i:i)

        if ( ch_is_digit ( c ) ) then

          call digit_inc ( c )

          file_name(i:i) = c

          if ( c .ne. '0' ) then
            return
          end if

        end if

      end do

      return
      end
      subroutine hcell_dof_count ( region_density_x, region_density_y, 
     &  dof_num )
!
!*******************************************************************************
!
!! HCELL_DOF_COUNT determines the number of degrees of freedom in the region.
!
!
!  Discussion:
!
!    The count produced by this routine does not correspond to current practice,
!    because degrees of freedom associated with boundary conditions and othe
!    specifications are now not being counted.
!
!  Diagram:
!
!          +----------------------------+
!          |              :     :       |
!    row 3 |   (3,1)      :(3,2): (3,3) |
!          |              :     :       |
!          +--------------+.....+-------+
!                         |     |        
!    row 2     empty      |(2,2)|  empty
!                         |     |  
!          +--------------+.....+-------+
!          |              :     :       |
!    row 1 |   (1,1)      :(1,2): (1,3) |
!          |              :     :       |
!          +----------------------------+
!
!              col 1       col 2  col 3
!
!  Discussion:
!
!    The region is divided into a 3 by 3 grid.  Subregion ( I, J )
!    is divided into ELEMENT_DENSITY_X(J) * ELEMENT_DENSITY_Y(I) squares.
!    Then each square is split into two triangles, with the diagonal
!    going from the upper left to lower right.  
!
!    Each element, in turn, is made up of 6 nodes.  The corner
!    nodes have 3 degrees of freedom ( U, V, and P), while the
!    side nodes have 2 degrees of freedom (U and V only).
!
!  Modified:
!
!    17 April 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ELEMENT_DENSITY_X(3), the density of elements 
!    in the three columns.
!
!    Input, integer ELEMENT_DENSITY_Y(3), the density of elements 
!    in the three rows.
!
!    Output, integer DOF_NUM, the number of degrees of freedom.
!
      implicit none
!
      integer dof_num
      integer region_density_x(3)
      integer region_density_y(3)
      integer node_num_p
      integer node_num_u
!
      node_num_u = 
     &    ( 2 * region_density_x(1) + 1 ) 
     &  * ( 2 * region_density_y(1) + 1 ) 
     &  + ( 2 * region_density_x(1) + 1 ) 
     &  * ( 2 * region_density_y(3) + 1 ) 
     &  + ( 2 * region_density_x(2) - 1 ) 
     &  * ( 2 * region_density_y(1) + 1 ) 
     &  + ( 2 * region_density_x(2) + 1 ) 
     &  * ( 2 * region_density_y(2) - 1 ) 
     &  + ( 2 * region_density_x(2) - 1 ) 
     &  * ( 2 * region_density_y(3) + 1 ) 
     &  + ( 2 * region_density_x(3) + 1 ) 
     &  * ( 2 * region_density_y(1) + 1 ) 
     &  + ( 2 * region_density_x(3) + 1 ) 
     &  * ( 2 * region_density_y(3) + 1 )

      node_num_p = 
     &    ( region_density_x(1) + 1 ) 
     &  * ( region_density_y(1) + 1 ) 
     &  + ( region_density_x(1) + 1 ) 
     &  * ( region_density_y(3) + 1 ) 
     &  + ( region_density_x(2) - 1 ) 
     &  * ( region_density_y(1) + 1 ) 
     &  + ( region_density_x(2) + 1 ) 
     &  * ( region_density_y(2) - 1 ) 
     &  + ( region_density_x(2) - 1 ) 
     &  * ( region_density_y(3) + 1 ) 
     &  + ( region_density_x(3) + 1 ) 
     &  * ( region_density_y(1) + 1 ) 
     &  + ( region_density_x(3) + 1 ) 
     &  * ( region_density_y(3) + 1 ) 

      dof_num = 2 * node_num_u + node_num_p

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'HCELL_DOF_COUNT:'
      write ( *, '(a,i6)' ) '  Number of degrees of freedom = ', dof_num

      return
      end
      subroutine hcell_dof_set ( region_density_x, 
     &  region_density_y, maxnd, node_num, node_dof_index, dof )
!
!*******************************************************************************
!
!! HCELL_DOF_SET assigns degrees of freedom to each node.
!
!
!  Diagram:
!
!          +----------------------------+
!          |              :     :       |
!    row 3 |   (3,1)      :(3,2): (3,3) |
!          |              :     :       |
!          +--------------+.....+-------+
!                         |     |
!    row 2     empty      |(2,2)|  empty
!                         |     |
!          +--------------+.....+-------+
!          |              :     :       |
!    row 1 |   (1,1)      :(1,2): (1,3) |
!          |              :     :       |
!          +----------------------------+
!
!              col 1       col 2  col 3
!
!  Discussion:
!
!    The region is divided into a 3 by 3 grid.  Subregion ( I, J )
!    is divided into ELEMENT_DENSITY_X(J) * ELEMENT_DENSITY_Y(I) rectangles.
!    Then each rectangle is split into two triangles, with the diagonal
!    going from the upper left to lower right.
!
!  Modified:
!
!    19 April 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer REGION_DENSITY_X(3), the density of elements
!    in the three columns.
!
!    Input, integer REGION_DENSITY_Y(3), the density of elements
!    in the three rows.
!
!    Input, integer MAXND, the maximum number of nodes.
!
!    Input, integer NODE_NUM, the number of nodes.
!
!    Output, integer NODE_DOF_INDEX(NODE_NUM,3), the nodal degrees of freedom.
!
!  Local Parameters:
!
!    Integer DOF, contains the value of the most recently assigned degree of freedom.
!
      implicit none
!
      integer maxnd
!
      integer col
      integer dof
      logical found
      integer i
      integer j
      integer node
      integer node_dof_index(maxnd,3)
      integer node_num
      integer region_density_x(3)
      integer region_density_y(3)
      integer row
!
      node = 0
      dof = 0
      j = 0
!
!  Working in column 1, and NOT handling nodes in extreme right.
!
      do col = 1, 2 * region_density_x(1)
!
!  Working in row 1.
!
!  +--+--+--+
!  |        |
!  +--+  +--+
!     |  |
!  +--+  +--+
!  |11      |
!  +--+--+--+
!
        j = j + 1
        i = 0

        do row = 1, 2 * region_density_y(1) + 1

          i = i + 1
          node = node + 1

  
       if (col .eq. 1 ) then
         node_dof_index(node,1) = 0
         node_dof_index(node,2) = 0
         if ( row .ge. (2 * 2 * region_density_y(1) )/5 + 1
     &        .and.
     &        row .le. (3 * 2 * region_density_y(1) )/5 + 1 ) then 
            node_dof_index(node,1) = -3
         end if
       else
         if (row .eq. 1 ) then
            node_dof_index(node,1) = 0
            node_dof_index(node,2) = 0
         else
           if (row .eq. 2 * region_density_y(1) + 1 ) then
             node_dof_index(node,1) = 0
             node_dof_index(node,2) = 0
             if ( col .ge. (10 * 2 * region_density_x(1) )/45 + 1
     &            .and.
     &            col .le. (11 * 2 * region_density_x(1) )/45 + 1 ) then 
                node_dof_index(node,2) = -3
             end if
             if ( col .ge. (22 * 2 * region_density_x(1) )/45 + 1
     &            .and.
     &            col .le. (23 * 2 * region_density_x(1) )/45 + 1 ) then 
                node_dof_index(node,2) = -3
             end if
             if ( col .ge. (34 * 2 * region_density_x(1) )/45 + 1
     &            .and.
     &            col .le. (35 * 2 * region_density_x(1) )/45 + 1 ) then 
                node_dof_index(node,2) = -3
             end if
           else

             dof = dof + 1
             node_dof_index(node,1) = dof

             dof = dof + 1
             node_dof_index(node,2) = dof

           end if
         end if
       end if


          if ( mod ( j, 2 ) .eq. 1 .and. mod ( i, 2 ) .eq. 1 ) then
            dof = dof + 1
            node_dof_index(node,3) = dof
          else
            node_dof_index(node,3) = 0
          end if



        end do
!
!  Working in row 3.
!
!  +--+--+--+
!  |31     |
!  +--+  +--+
!     |  |
!  +--+  +--+
!  |        |
!  +--+--+--+
!
        i = i + 2 * region_density_y(2) - 1

        do row = 1, 2 * region_density_y(3) + 1

          i = i + 1
          node = node + 1


       if (col .eq. 1 ) then
         node_dof_index(node,1) = 0
         node_dof_index(node,2) = 0
         if ( row .ge. (2 * 2 * region_density_y(3) )/5 + 1
     &        .and.
     &        row .le. (3 * 2 * region_density_y(3) )/5 + 1 ) then 
            node_dof_index(node,1) = -1
         end if
       else
         if (row .eq. 2 * region_density_y(3) + 1 ) then
           node_dof_index(node,1) = 0
           node_dof_index(node,2) = 0
           if ( col .ge. (10 * 2 * region_density_x(1) )/45 + 1
     &          .and.
     &          col .le. (11 * 2 * region_density_x(1) )/45 + 1 ) then 
              node_dof_index(node,2) = -5
           end if
           if ( col .ge. (22 * 2 * region_density_x(1) )/45 + 1
     &          .and.
     &          col .le. (23 * 2 * region_density_x(1) )/45 + 1 ) then 
              node_dof_index(node,2) = -5
           end if
           if ( col .ge. (34 * 2 * region_density_x(1) )/45 + 1
     &          .and.
     &          col .le. (35 * 2 * region_density_x(1) )/45 + 1 ) then 
              node_dof_index(node,2) = -5
           end if
         else
           if (row .eq. 1 ) then
             node_dof_index(node,1) = 0
             node_dof_index(node,2) = 0
             if ( col .ge. (10 * 2 * region_density_x(1) )/45 + 1
     &            .and.
     &            col .le. (11 * 2 * region_density_x(1) )/45 + 1 ) then 
                node_dof_index(node,2) = -1
             end if
             if ( col .ge. (22 * 2 * region_density_x(1) )/45 + 1
     &            .and.
     &            col .le. (23 * 2 * region_density_x(1) )/45 + 1 ) then 
                node_dof_index(node,2) = -1
             end if
             if ( col .ge. (34 * 2 * region_density_x(1) )/45 + 1
     &            .and.
     &            col .le. (35 * 2 * region_density_x(1) )/45 + 1 ) then 
                node_dof_index(node,2) = -1
             end if
           else

             dof = dof + 1
             node_dof_index(node,1) = dof

             dof = dof + 1
             node_dof_index(node,2) = dof

           end if
         end if
       end if


          if ( mod ( j, 2 ) .eq. 1 .and. mod ( i, 2 ) .eq. 1 ) then
            dof = dof + 1
            node_dof_index(node,3) = dof
          else
            node_dof_index(node,3) = 0
          end if

        end do

      end do
!
!  Working in column 2, including nodes on extreme right and left.
!
      do col = 1, 2 * region_density_x(2) + 1
!
!  Working in row 1.
!
!  +--+--+--+
!  |        |
!  +--+  +--+
!     |  |
!  +--+  +--+
!  |   12   |
!  +--+--+--+
!
        j = j + 1
        i = 0

        do row = 1, 2 * region_density_y(1) + 1

          i = i + 1
          node = node + 1


       if (row .eq. 1 ) then
         node_dof_index(node,1) = 0
         node_dof_index(node,2) = 0
       else
         if (row .eq. 2 * region_density_y(1) + 1 
     &       .and.
     &       col .eq. 1 ) then
               node_dof_index(node,1) = 0
               node_dof_index(node,2) = 0
         else
           if (row .eq. 2 * region_density_y(1) + 1 
     &         .and.
     &         col .eq. 2 * region_density_x(2) + 1 ) then
                 node_dof_index(node,1) = 0
                 node_dof_index(node,2) = 0
           else

             dof = dof + 1
             node_dof_index(node,1) = dof

             dof = dof + 1
             node_dof_index(node,2) = dof

           end if
         end if
       end if


          if ( mod ( j, 2 ) .eq. 1 .and. mod ( i, 2 ) .eq. 1 ) then
            dof = dof + 1
            node_dof_index(node,3) = dof
          else
            node_dof_index(node,3) = 0
          end if

        end do
!
!  Working in row 2.
!
!  +--+--+--+
!  |        |
!  +--+  +--+
!     |22|
!  +--+  +--+
!  |        |
!  +--+--+--+
!
        do row = 2, 2 * region_density_y(2)

          i = i + 1
          node = node + 1

       if (col .eq. 1 ) then
          node_dof_index(node,1) = 0
          node_dof_index(node,2) = 0
       else
         if (col .eq. 2 * region_density_x(2) + 1 ) then
            node_dof_index(node,1) = 0
            node_dof_index(node,2) = 0
         else

             dof = dof + 1
             node_dof_index(node,1) = dof

             dof = dof + 1
             node_dof_index(node,2) = dof

         end if
       end if

          if ( mod ( j, 2 ) .eq. 1 .and. mod ( i, 2 ) .eq. 1 ) then
            dof = dof + 1
            node_dof_index(node,3) = dof
          else
            node_dof_index(node,3) = 0
          end if

        end do
!
!  Working in row 3.
!
!  +--+--+--+
!  |   32   |
!  +--+  +--+
!     |  |
!  +--+  +--+
!  |        |
!  +--+--+--+
!
        do row = 1, 2 * region_density_y(3) + 1

          i = i + 1
          node = node + 1

      
       if (row .eq. 1 
     &     .and.
     &     col .eq. 1 ) then
          node_dof_index(node,1) = 0
          node_dof_index(node,2) = 0
       else
         if (row .eq. 1 
     &       .and.
     &       col .eq. 2 * region_density_x(2) + 1 ) then
            node_dof_index(node,1) = 0
            node_dof_index(node,2) = 0
         else
           if (row .eq. 2 * region_density_y(3) + 1 ) then
             node_dof_index(node,1) = 0
             node_dof_index(node,2) = 0
             if (col .ge. (5 * 2 * region_density_x(2) )/15
     &           .and.
     &           col .le. (10 * 2 * region_density_x(2) )/15 ) then
               dof = dof + 1
!
!  18 May 2004, following line changed for MDG:
!
               node_dof_index(node,2) = dof
             end if
           else

             dof = dof + 1
             node_dof_index(node,1) = dof

             dof = dof + 1
             node_dof_index(node,2) = dof

           end if
         end if
       end if

          if ( mod ( j, 2 ) .eq. 1 .and. mod ( i, 2 ) .eq. 1 ) then
            dof = dof + 1
            node_dof_index(node,3) = dof
          else
            node_dof_index(node,3) = 0
          end if

        end do

      end do
!
!  Working in column 3, and NOT handling nodes in extreme left.
!
      do col = 2, 2 * region_density_x(3) + 1
!
!  Working in row 1.
!
!  +--+--+--+
!  |        |
!  +--+  +--+
!     |  |
!  +--+  +--+
!  |      13|
!  +--+--+--+
!
        j = j + 1
        i = 0

        do row = 1, 2 * region_density_y(1) + 1

          i = i + 1
          node = node + 1

       if (col .eq. 2 * region_density_x(3) + 1 ) then
         node_dof_index(node,1) = 0
         node_dof_index(node,2) = 0
         if ( row .ge. (2 * 2 * region_density_y(1) )/5 + 1
     &        .and.
     &        row .le. (3 * 2 * region_density_y(1) )/5 + 1 ) then 
            node_dof_index(node,1) = -4
         end if
       else
         if (row .eq. 1 ) then
            node_dof_index(node,1) = 0
            node_dof_index(node,2) = 0
         else
           if (row .eq. 2 * region_density_y(1) + 1 ) then
             node_dof_index(node,1) = 0
             node_dof_index(node,2) = 0
             if ( col .ge. (10 * 2 * region_density_x(3) )/45 + 1
     &            .and.
     &            col .le. (11 * 2 * region_density_x(3) )/45 + 1 ) then 
                node_dof_index(node,2) = -4
             end if
             if ( col .ge. (22 * 2 * region_density_x(3) )/45 + 1
     &            .and.
     &            col .le. (23 * 2 * region_density_x(3) )/45 + 1 ) then 
                node_dof_index(node,2) = -4
             end if
             if ( col .ge. (34 * 2 * region_density_x(3) )/45 + 1
     &            .and.
     &            col .le. (35 * 2 * region_density_x(3) )/45 + 1 ) then 
                node_dof_index(node,2) = -4
             end if
           else

             dof = dof + 1
             node_dof_index(node,1) = dof

             dof = dof + 1
             node_dof_index(node,2) = dof

           end if
         end if
       end if


          if ( mod ( j, 2 ) .eq. 1 .and. mod ( i, 2 ) .eq. 1 ) then
            dof = dof + 1
            node_dof_index(node,3) = dof
          else
            node_dof_index(node,3) = 0
          end if

        end do
!
!  Working in row 3.
!
!  +--+--+--+
!  |      33|
!  +--+  +--+
!     |  |
!  +--+  +--+
!  |        |
!  +--+--+--+
!
        i = i + 2 * region_density_y(2) - 1

        do row = 1, 2 * region_density_y(3) + 1

          i = i + 1
          node = node + 1

       if (col .eq. 2 * region_density_x(3) + 1 ) then
         node_dof_index(node,1) = 0
         node_dof_index(node,2) = 0
         if ( row .ge. (2 * 2 * region_density_y(3) )/5 + 1
     &        .and.
     &        row .le. (3 * 2 * region_density_y(3) )/5 + 1 ) then 
            node_dof_index(node,1) = -2
         end if
       else
         if (row .eq. 2 * region_density_y(3) + 1 ) then
           node_dof_index(node,1) = 0
           node_dof_index(node,2) = 0
           if ( col .ge. (10 * 2 * region_density_x(3) )/45 + 1
     &          .and.
     &          col .le. (11 * 2 * region_density_x(3) )/45 + 1 ) then 
              node_dof_index(node,2) = -6
           end if
           if ( col .ge. (22 * 2 * region_density_x(3) )/45 + 1
     &          .and.
     &          col .le. (23 * 2 * region_density_x(3) )/45 + 1 ) then 
              node_dof_index(node,2) = -6
           end if
           if ( col .ge. (34 * 2 * region_density_x(3) )/45 + 1
     &          .and.
     &          col .le. (35 * 2 * region_density_x(3) )/45 + 1 ) then 
              node_dof_index(node,2) = -6
           end if
         else
           if (row .eq. 1 ) then
             node_dof_index(node,1) = 0
             node_dof_index(node,2) = 0
             if ( col .ge. (10 * 2 * region_density_x(3) )/45 + 1
     &            .and.
     &            col .le. (11 * 2 * region_density_x(3) )/45 + 1 ) then 
                node_dof_index(node,2) = -2
             end if
             if ( col .ge. (22 * 2 * region_density_x(3) )/45 + 1
     &            .and.
     &            col .le. (23 * 2 * region_density_x(3) )/45 + 1 ) then 
                node_dof_index(node,2) = -2
             end if
             if ( col .ge. (34 * 2 * region_density_x(3) )/45 + 1
     &            .and.
     &            col .le. (35 * 2 * region_density_x(3) )/45 + 1 ) then 
                node_dof_index(node,2) = -2
             end if
           else

             dof = dof + 1
             node_dof_index(node,1) = dof

             dof = dof + 1
             node_dof_index(node,2) = dof

           end if
         end if
       end if


          if ( mod ( j, 2 ) .eq. 1 .and. mod ( i, 2 ) .eq. 1 ) then
            dof = dof + 1
            node_dof_index(node,3) = dof
          else
            node_dof_index(node,3) = 0
          end if

        end do

      end do

      return
      end
      subroutine hcell_element_count ( region_density_x, 
     &  region_density_y, element_num )
!
!*******************************************************************************
!
!! HCELL_ELEMENT_COUNT determines the number of elements in the region.
!
!
!  Diagram:
!
!          +----------------------------+
!          |              :     :       |
!    row 3 |   (3,1)      :(3,2): (3,3) |
!          |              :     :       |
!          +--------------+.....+-------+
!                         |     |        
!    row 2     empty      |(2,2)|  empty
!                         |     |  
!          +--------------+.....+-------+
!          |              :     :       |
!    row 1 |   (1,1)      :(1,2): (1,3) |
!          |              :     :       |
!          +----------------------------+
!
!              col 1       col 2  col 3
!
!  Discussion:
!
!    The region is divided into a 3 by 3 grid.  Subregion ( I, J )
!    is divided into ELEMENT_DENSITY_X(J) * ELEMENT_DENSITY_Y(I) squares.
!    Then each square is split into two triangles, with the diagonal
!    going from the upper left to lower right.  
!
!  Modified:
!
!    25 March 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ELEMENT_DENSITY_X(3), the density of elements 
!    in the three columns.
!
!    Input, integer ELEMENT_DENSITY_Y(3), the density of elements 
!    in the three rows.
!
!    Output, integer ELEMENT_NUM, the number of elements.
!
      implicit none
!
      integer element_num
      integer region_density_x(3)
      integer region_density_y(3)
!
      element_num = 
     &    2 * region_density_x(1) * region_density_y(1) 
     &  + 2 * region_density_x(1) * region_density_y(3) 
     &  + 2 * region_density_x(2) * region_density_y(1) 
     &  + 2 * region_density_x(2) * region_density_y(2) 
     &  + 2 * region_density_x(2) * region_density_y(3) 
     &  + 2 * region_density_x(3) * region_density_y(1) 
     &  + 2 * region_density_x(3) * region_density_y(3)

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'HCELL_ELEMENT_COUNT:'
      write ( *, '(a,i6)' ) '  Number of elements = ', element_num

      return
      end
      subroutine hcell_element_node ( region_density_x, 
     &  region_density_y, maxel, element_num, nnodes, element_node )
!
!*******************************************************************************
!
!! HCELL_ELEMENT_NODE determines the nodes that make up each element.
!
!
!  Diagram:
!
!          +----------------------------+
!          |              :     :       |
!    row 3 |   (3,1)      :(3,2): (3,3) |
!          |              :     :       |
!          +--------------+.....+-------+
!                         |     |        
!    row 2     empty      |(2,2)|  empty
!                         |     |  
!          +--------------+.....+-------+
!          |              :     :       |
!    row 1 |   (1,1)      :(1,2): (1,3) |
!          |              :     :       |
!          +----------------------------+
!
!              col 1       col 2  col 3
!
!  Discussion:
!
!    The region is divided into a 3 by 3 grid.  Subregion ( I, J )
!    is divided into ELEMENT_DENSITY_X(J) * ELEMENT_DENSITY_Y(I) squares.
!    Then each square is split into two triangles, with the diagonal
!    going from the upper left to lower right.  
!
!  Modified:
!
!    19 April 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ELEMENT_DENSITY_X(3), the density of elements 
!    in the three columns.
!
!    Input, integer ELEMENT_DENSITY_Y(3), the density of elements 
!    in the three rows.
!
!    Input, integer MAXEL, the maximum number of elements.
!
!    Input, integer ELEMENT_NUM, the number of elements.
!
!    Input, integer NNODES, the number of nodes per element.
!
!    Output, integer ELEMENT_NODE(MAXEL,NNODES), the nodes that make up
!    each element.
!
      implicit none
!
      integer maxel
      integer nnodes
!
      integer col
      integer col2
      integer element
      integer element_node(maxel,nnodes)
      integer element_num
      integer inc1
      integer inc2
      integer node_sw
      integer region_density_x(3)
      integer region_density_y(3)
      integer row
      integer row2
!
      element = 0

      do col = 1, 3
        do col2 = 1, region_density_x(col)

          do row = 1, 3

            if ( row /= 2 .or. col == 2 ) then

              if ( col == 1 ) then

                if ( col2 < region_density_x(1) ) then

                  if ( row == 1 ) then

                    if ( col2 == 1 ) then
                      node_sw = 1
                    else
                      node_sw = node_sw + inc1 + 1
                    end if

                  else

                    node_sw = node_sw + 1

                  end if

                  inc1 = ( 2 * region_density_y(1) + 1 ) 
     &                  + ( 2 * region_density_y(3) + 1 )

                  inc2 = ( 2 * region_density_y(1) + 1 ) 
     &                  + ( 2 * region_density_y(3) + 1 )

                else if ( row == 1 ) then

                  node_sw = node_sw + inc1 + 1

                   inc1 = ( 2 * region_density_y(1) + 1 ) 
     &                  + ( 2 * region_density_y(3) + 1 )

                   inc2 = ( 2 * region_density_y(1) + 1 ) 
     &                  + ( 2 * region_density_y(3) + 1 )

                else if ( row == 3 ) then

                  node_sw = node_sw + 1

                  inc1 = ( 2 * region_density_y(1) + 1 ) 
     &                  + ( 2 * region_density_y(3) + 1 )

                  inc2 = ( 2 * region_density_y(1) + 1 ) 
     &                  + ( 2 * region_density_y(2) - 1 ) 
     &                  + ( 2 * region_density_y(3) + 1 )

                end if

              else if ( col == 2 ) then

                if ( row == 1 ) then

                  node_sw = node_sw + inc1 + 1

                end if

                inc1 = ( 2 * region_density_y(1) + 1 ) 
     &                     + ( 2 * region_density_y(2) - 1 ) 
     &                + ( 2 * region_density_y(3) + 1 )

                inc2 = ( 2 * region_density_y(1) + 1 ) 
     &                 + ( 2 * region_density_y(2) - 1 ) 
     &                 + ( 2 * region_density_y(3) + 1 )

              else if ( col == 3 ) then

                if ( 1 < col2 ) then

                  if ( row == 1 ) then
                    node_sw = node_sw + inc1 + 1
                  else
                    node_sw = node_sw + 1
                  end if

                  inc1 = ( 2 * region_density_y(1) + 1 ) 
     &                   + ( 2 * region_density_y(3) + 1 )

                  inc2 = ( 2 * region_density_y(1) + 1 ) 
     &                   + ( 2 * region_density_y(3) + 1 )

                else if ( row == 1 ) then

                  node_sw = node_sw + inc1 + 1

                  inc1 = ( 2 * region_density_y(1) + 1 ) 
     &                   + ( 2 * region_density_y(2) - 1 ) 
     &                   + ( 2 * region_density_y(3) + 1 )

                  inc2 = ( 2 * region_density_y(1) + 1 ) 
     &                   + ( 2 * region_density_y(3) + 1 )

                else if ( row == 3 ) then

                  node_sw = node_sw 
     &               + ( 2 * region_density_y(2) - 1 ) + 1

                  inc1 = ( 2 * region_density_y(1) + 1 ) 
     &                  + ( 2 * region_density_y(3) + 1 )

                  inc2 = ( 2 * region_density_y(1) + 1 ) 
     &                  + ( 2 * region_density_y(3) + 1 )

                end if

              end if

              do row2 = 1, region_density_y(row)

                element = element + 1
                element_node(element,1) = node_sw
                element_node(element,2) = node_sw + inc1 + inc2
                element_node(element,3) = node_sw               + 2
                element_node(element,4) = node_sw + inc1
                element_node(element,5) = node_sw + inc1        + 1
                element_node(element,6) = node_sw               + 1

                element = element + 1
                element_node(element,1) = node_sw + inc1 + inc2 + 2
                element_node(element,2) = node_sw               + 2
                element_node(element,3) = node_sw + inc1 + inc2
                element_node(element,4) = node_sw + inc1        + 2
                element_node(element,5) = node_sw + inc1        + 1
                element_node(element,6) = node_sw + inc1 + inc2 + 1

                node_sw = node_sw + 2

              end do

            end if

          end do
        end do
      end do

      return
      end
      subroutine hcell_node_count ( region_density_x, region_density_y, 
     &  node_num )
!
!*******************************************************************************
!
!! HCELL_NODE_NUM determines the number of nodes in the region.
!
!
!  Diagram:
!
!          +----------------------------+
!          |              :     :       |
!    row 3 |   (3,1)      :(3,2): (3,3) |
!          |              :     :       |
!          +--------------+.....+-------+
!                         |     |        
!    row 2     empty      |(2,2)|  empty
!                         |     |  
!          +--------------+.....+-------+
!          |              :     :       |
!    row 1 |   (1,1)      :(1,2): (1,3) |
!          |              :     :       |
!          +----------------------------+
!
!              col 1       col 2  col 3
!
!  Discussion:
!
!    The region is divided into a 3 by 3 grid.  Subregion ( I, J )
!    is divided into ELEMENT_DENSITY_X(J) * ELEMENT_DENSITY_Y(I) squares.
!    Then each square is split into two triangles, with the diagonal
!    going from the upper left to lower right.  
!
!  Modified:
!
!    17 April 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ELEMENT_DENSITY_X(3), the density of elements 
!    in the three columns.
!
!    Input, integer ELEMENT_DENSITY_Y(3), the density of elements 
!    in the three rows.
!
!    Output, integer NODE_NUM, the number of nodes.
!
      implicit none
!
      integer node_num
      integer region_density_x(3)
      integer region_density_y(3)
! 
!  Count the nodes.
!
      node_num = 
     &    ( 2 * region_density_x(1) + 1 ) 
     &  * ( 2 * region_density_y(1) + 1 ) 
     &  + ( 2 * region_density_x(1) + 1 ) 
     &  * ( 2 * region_density_y(3) + 1 ) 
     &  + ( 2 * region_density_x(2) - 1 ) 
     &  * ( 2 * region_density_y(1) + 1 ) 
     &  + ( 2 * region_density_x(2) + 1 ) 
     &  * ( 2 * region_density_y(2) - 1 ) 
     &  + ( 2 * region_density_x(2) - 1 ) 
     &  * ( 2 * region_density_y(3) + 1 ) 
     &  + ( 2 * region_density_x(3) + 1 ) 
     &  * ( 2 * region_density_y(1) + 1 ) 
     &  + ( 2 * region_density_x(3) + 1 ) 
     &  * ( 2 * region_density_y(3) + 1 )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'HCELL_NODE_COUNT:'
      write ( *, '(a,i6)' ) '  Number of nodes = ', node_num

      return
      end
      subroutine hcell_node_xy ( region_density_x, region_density_y, 
     &  node_num, region_x, region_y, node_x, node_y )
!
!*******************************************************************************
!
!! HCELL_NODE_XY assigns coordinates to each node.
!
!
!  Diagram:
!
!          +----------------------------+
!          |              :     :       |
!    row 3 |   (3,1)      :(3,2): (3,3) |
!          |              :     :       |
!          +--------------+.....+-------+
!                         |     |        
!    row 2     empty      |(2,2)|  empty
!                         |     |  
!          +--------------+.....+-------+
!          |              :     :       |
!    row 1 |   (1,1)      :(1,2): (1,3) |
!          |              :     :       |
!          +----------------------------+
!
!              col 1       col 2  col 3
!
!  Discussion:
!
!    The region is divided into a 3 by 3 grid.  Subregion ( I, J )
!    is divided into ELEMENT_DENSITY_X(J) * ELEMENT_DENSITY_Y(I) squares.
!    Then each square is split into two triangles, with the diagonal
!    going from the upper left to lower right.  
!
!  Modified:
!
!    17 April 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer REGION_DENSITY_X(3), the density of elements 
!    in the three columns.
!
!    Input, integer REGION_DENSITY_Y(3), the density of elements 
!    in the three rows.
!
!    Input, integer NODE_NUM, the number of nodes.
!
!    Input, double precision REGION_X(4), REGION_Y(4), the coordinates of 
!    breakpoints that define 9 logical subrectangles.
!
!    Output, double precision NODE_X(NODE_NUM), NODE_Y(NODE_NUM),
!    the X and Y coordinates of the nodes.
!
      implicit none
!
      integer node_num
!
      integer col
      integer i
      integer j
      integer node
      double precision node_x(node_num)
      double precision node_y(node_num)
      integer region_density_x(3)
      integer region_density_y(3)
      double precision region_x(4)
      double precision region_y(4)
      integer row
!
      node = 0
      j = 0
!
!  Working in column 1, except for the extreme right.
!
      do col = 1, 2 * region_density_x(1)

        j = j + 1
        i = 0
!
!  Working in row 1.
!
!  +--+--+--+
!  |        |
!  +--+  +--+
!     |  |
!  +--+  +--+
!  |11      |
!  +--+--+--+
!
        do row = 1, 2 * region_density_y(1) + 1

          i = i + 1
          node = node + 1

          node_x(node) = 
     &      ( dble ( 2 * region_density_x(1) + 1 - col ) * region_x(1)
     &      + dble (                         - 1 + col ) * region_x(2) )
     &      / dble ( 2 * region_density_x(1)           )

          node_y(node) = 
     &      ( dble ( 2 * region_density_y(1) + 1 - row ) * region_y(1)
     &      + dble (                         - 1 + row ) * region_y(2) ) 
     &      / dble ( 2 * region_density_y(1)           )

        end do
!
!  Working in row 3.
!
!  +--+--+--+
!  |31      |
!  +--+  +--+
!     |  |
!  +--+  +--+
!  |        |
!  +--+--+--+
!
        i = i + 2 * region_density_y(2) - 1

        do row = 1, 2 * region_density_y(3) + 1

          i = i + 1
          node = node + 1

          node_x(node) = 
     &      ( dble ( 2 * region_density_x(1) + 1 - col ) * region_x(1)
     &      + dble (                         - 1 + col ) * region_x(2) )
     &      / dble ( 2 * region_density_x(1)           )

          node_y(node) = 
     &      ( dble ( 2 * region_density_y(3) + 1 - row ) * region_y(3) 
     &      + dble (                         - 1 + row ) * region_y(4) )
     &      / dble ( 2 * region_density_y(3)           )

        end do

      end do
!
!  Working in column 2, including extreme left and right.
!
      do col = 1, 2 * region_density_x(2) + 1
!
!  Working in row 1.
!
!  +--+--+--+
!  |        |
!  +--+  +--+
!     |  |
!  +--+  +--+
!  |   12   |
!  +--+--+--+
!
        j = j + 1
        i = 0

        do row = 1, 2 * region_density_y(1) + 1

          i = i + 1
          node = node + 1

          node_x(node) = 
     &      ( dble ( 2 * region_density_x(2) + 1 - col ) * region_x(2)   
     &      + dble (                         - 1 + col ) * region_x(3) ) 
     &      / dble ( 2 * region_density_x(2)           )

          node_y(node) = 
     &      ( dble ( 2 * region_density_y(1) + 1 - row ) * region_y(1)   
     &      + dble (                         - 1 + row ) * region_y(2) ) 
     &      / dble ( 2 * region_density_y(1)           )

        end do
!
!  Working in row 2.
!
!  +--+--+--+
!  |        |
!  +--+  +--+
!     |22|
!  +--+  +--+
!  |        |
!  +--+--+--+
!
        do row = 2, 2 * region_density_y(2)

          i = i + 1
          node = node + 1

          node_x(node) = 
     &      ( dble ( 2 * region_density_x(2) + 1 - col ) * region_x(2)   
     &      + dble (                         - 1 + col ) * region_x(3) ) 
     &      / dble ( 2 * region_density_x(2)           )

          node_y(node) = 
     &      ( dble ( 2 * region_density_y(2) + 1 - row ) * region_y(2)   
     &      + dble (                         - 1 + row ) * region_y(3) ) 
     &      / dble ( 2 * region_density_y(2)           )

        end do
!
!  Working in row 3.
!
!  +--+--+--+
!  |   32   |
!  +--+  +--+
!     |  |
!  +--+  +--+
!  |        |
!  +--+--+--+
!
        do row = 1, 2 * region_density_y(3) + 1

          i = i + 1
          node = node + 1

          node_x(node) = 
     &      ( dble ( 2 * region_density_x(2) + 1 - col ) * region_x(2)   
     &      + dble (                         - 1 + col ) * region_x(3) ) 
     &      / dble ( 2 * region_density_x(2)           )

          node_y(node) = 
     &      ( dble ( 2 * region_density_y(3) + 1 - row ) * region_y(3)   
     &      + dble (                         - 1 + row ) * region_y(4) ) 
     &      / dble ( 2 * region_density_y(3)           )

        end do

      end do
!
!  Working in column 3, except for extreme left.
!
      do col = 2, 2 * region_density_x(3) + 1
!
!  Working in row 1.
!
!  +--+--+--+
!  |        |
!  +--+  +--+
!     |  |
!  +--+  +--+
!  |      13|
!  +--+--+--+
!
        j = j + 1
        i = 0

        do row = 1, 2 * region_density_y(1) + 1

          i = i + 1
          node = node + 1

          node_x(node) = 
     &      ( dble ( 2 * region_density_x(3) + 1 - col ) * region_x(3)   
     &      + dble (                         - 1 + col ) * region_x(4) ) 
     &      / dble ( 2 * region_density_x(3)           )

          node_y(node) = 
     &      ( dble ( 2 * region_density_y(1) + 1 - row ) * region_y(1)   
     &      + dble (                         - 1 + row ) * region_y(2) ) 
     &      / dble ( 2 * region_density_y(1)           )

        end do

        i = i + 2 * region_density_y(2) - 1
!
!  Working in row 3.
!
!  +--+--+--+
!  |      33|
!  +--+  +--+
!     |  |
!  +--+  +--+
!  |        |
!  +--+--+--+
!
        do row = 1, 2 * region_density_y(3) + 1
          i = i + 1
          node = node + 1

          node_x(node) = 
     &      ( dble ( 2 * region_density_x(3) + 1 - col ) * region_x(3)   
     &      + dble (                         - 1 + col ) * region_x(4) ) 
     &      / dble ( 2 * region_density_x(3)           )

          node_y(node) = 
     &      ( dble ( 2 * region_density_y(3) + 1 - row ) * region_y(3)   
     &      + dble (                         - 1 + row ) * region_y(4) ) 
     &      / dble ( 2 * region_density_y(3)           )
        end do

      end do

      return
      end
      function idamax ( n, dx, incx )
!
!***********************************************************************
!
!! IDAMAX finds the vector element of largest magnitude.
!
!
!  Author:
!
!    Jack Dongarra
!
!  Parameters:
!
!    Input, integer N, the number of elements in the vector.
!
!    Input, double precision DX(*), the vector to examine.
!
!    Input, integer INCX, the increment between successive entries.
!
!    Output, integer IDAMAX, the index of the vector element of
!    largest magnitude.
!
      implicit none
!
      double precision dmax
      double precision dx(*)
      integer i
      integer idamax
      integer incx
      integer ix
      integer n
!
      idamax = 0
      if ( n < 1 ) then
        return
      end if

      idamax = 1
      if ( n == 1 ) then
        return
      end if

      if ( incx == 1 ) then

        dmax = dabs ( dx(1) )
        do i = 2, n
          if ( dmax < dabs ( dx(i) ) ) then
            idamax = i
            dmax = dabs ( dx(i) )
          end if
        end do

      else

        ix = 1
        dmax = dabs ( dx(1) )
        ix = ix + incx
        do i = 2, n
          if ( dmax < dabs ( dx(ix) ) ) then
            idamax = i
            dmax = dabs ( dx(ix) )
            ix = ix + incx
          end if
        end do

      end if

      return
      end
      subroutine node_eps ( node_file_name, np, node_mask, xc, yc, 
     &  title )
!
!***********************************************************************
!
!! NODE_EPS creates an EPS file containing an image of the nodes.
!
!
!  Modified:
!
!    20 April 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) NODE_FILE_NAME, the name of the file to create.
!
!    Input, integer NP, the number of nodes.
!
!    Input, logical NODE_MASK(NP), is TRUE for those nodes to be plotted.
!
!    Input, double precision XC(NP), YC(NP), the coordinates of the nodes.
!
!    Input, character ( len = * ) TITLE, a title for the plot.
!
      implicit none
!
      integer np
!
      double precision ave_x
      double precision ave_y
      integer circle_size
      double precision d_huge
      double precision dif
      integer eps_unit
      integer eps_x
      integer eps_y
      character ( len = * ) node_file_name
      integer i
      integer ios
      integer ip1
      integer j
      integer k
      integer node
      logical node_mask(np)
      integer node_unmasked_num
      double precision node_x_max
      double precision node_x_min
      double precision node_y_max
      double precision node_y_min
      double precision scale
      character ( len = 40 ) string
      character ( len = * ) title
      double precision xc(np)
      double precision yc(np)
!
!  Determine the range and number of the unmasked nodes.
!
      node_x_min =  d_huge ( )
      node_x_max = -d_huge ( )
      node_y_min =  d_huge ( )
      node_y_max = -d_huge ( )

      node_unmasked_num = 0

      do node = 1, np
        if ( node_mask(node) ) then
          node_unmasked_num = node_unmasked_num + 1
          node_x_min = min ( node_x_min, xc(node) )
          node_x_max = max ( node_x_max, xc(node) )
          node_y_min = min ( node_y_min, yc(node) )
          node_y_max = max ( node_y_max, yc(node) )
        end if
      end do

      if ( node_unmasked_num <= 200 ) then
        circle_size = 3
      else if ( node_unmasked_num <= 800 ) then
        circle_size = 2
      else
        circle_size = 1
      end if

      if ( node_y_max - node_y_min < node_x_max - node_x_min ) then
        scale = node_x_max - node_x_min
        dif = ( node_x_max - node_x_min ) - ( node_y_max - node_y_min )
        node_y_max = node_y_max + 0.5 * dif
        node_y_min = node_y_min - 0.5 * dif
      else
        scale = node_y_max - node_y_min
        dif = ( node_y_max - node_y_min ) - ( node_x_max - node_x_min )
        node_x_max = node_x_max + 0.5 * dif
        node_x_min = node_x_min - 0.5 * dif
      end if

      eps_unit = 1

      open ( unit = eps_unit, file = node_file_name, 
     &  status = 'unknown' )

      write ( eps_unit, '(a)' ) '%!PS-Adobe-3.0 EPSF-3.0'
      write ( eps_unit, '(a)' ) '%%Creator: node_eps(hcell.f)'
      write ( eps_unit, '(a,a)' ) '%%Title: ', node_file_name
      write ( eps_unit, '(a)' ) '%%Pages: 1'
      write ( eps_unit, '(a)' ) '%%BoundingBox:    36    36   576   756'
      write ( eps_unit, '(a)' ) '%%Document-Fonts: Times-Roman'
      write ( eps_unit, '(a)' ) '%%LanguageLevel: 1'
      write ( eps_unit, '(a)' ) '%%EndComments'
      write ( eps_unit, '(a)' ) '%%BeginProlog'
      write ( eps_unit, '(a)' ) '/inch {72 mul} def'
      write ( eps_unit, '(a)' ) '%%EndProlog'
      write ( eps_unit, '(a)' ) '%%Page:      1     1'
      write ( eps_unit, '(a)' ) 'save'
      write ( eps_unit, '(a)' ) '%'
      write ( eps_unit, '(a)' ) '% Set RGB line color.'
      write ( eps_unit, '(a)' ) '%'
      write ( eps_unit, '(a)' ) ' 0.9000 0.9000 0.9000 setrgbcolor'
      write ( eps_unit, '(a)' ) '%'
      write ( eps_unit, '(a)' ) '% Draw a gray border around the page.'
      write ( eps_unit, '(a)' ) '%'
      write ( eps_unit, '(a)' ) 'newpath'
      write ( eps_unit, '(a)' ) '    36   126 moveto'
      write ( eps_unit, '(a)' ) '   576   126 lineto'
      write ( eps_unit, '(a)' ) '   576   666 lineto'
      write ( eps_unit, '(a)' ) '    36   666 lineto'
      write ( eps_unit, '(a)' ) '    36   126 lineto'
      write ( eps_unit, '(a)' ) 'stroke'
      write ( eps_unit, '(a)' ) '%'
      write ( eps_unit, '(a)' ) '% Set RGB line color.'
      write ( eps_unit, '(a)' ) '%'
      write ( eps_unit, '(a)' ) ' 0.0000 0.0000 0.0000 setrgbcolor'

      write ( eps_unit, '(a)' ) '%'
      write ( eps_unit, '(a)' ) '%  Label the plot:'
      write ( eps_unit, '(a)' ) '%'
      write ( eps_unit, '(a)' ) ' 0.0000 0.0000 0.0000 setrgbcolor'
      write ( eps_unit, '(a)' ) 
     &   '/Times-Roman findfont 0.50 inch scalefont setfont'
      write ( eps_unit, '(a)' ) '    36   666 moveto'
      write ( eps_unit, '(a,a,a)' ) '(', title, ') show'

      write ( eps_unit, '(a)' ) '%'
      write ( eps_unit, '(a)' ) '% Define a clipping polygon'
      write ( eps_unit, '(a)' ) '%'
      write ( eps_unit, '(a)' ) '    36   126 moveto'
      write ( eps_unit, '(a)' ) '   576   126 lineto'
      write ( eps_unit, '(a)' ) '   576   666 lineto'
      write ( eps_unit, '(a)' ) '    36   666 lineto'
      write ( eps_unit, '(a)' ) '    36   126 lineto'
      write ( eps_unit, '(a)' ) 'clip newpath'

      write ( eps_unit, '(a)' ) '%'
      write ( eps_unit, '(a)' ) '%  Draw filled dots at each node:'
      write ( eps_unit, '(a)' ) '%'
      write ( eps_unit, '(a)' ) ' 0.0000 0.0000 1.0000 setrgbcolor'

      do node = 1, np

        if ( node_mask(node) ) then

          eps_x = int 
     &       ( ( node_x_max - xc(node)              ) *  61.0   
     &       + (            + xc(node) - node_x_min ) * 551.0 ) 
     &       / scale

          eps_y = int 
     &       ( ( node_y_max - yc(node)              ) * 151.0   
     &       + (              yc(node) - node_y_min ) * 641.0 ) 
     &       / scale

          write ( eps_unit, '(a,i4,2x,i4,2x,i4,a)' ) 
     &       'newpath  ', eps_x, eps_y, circle_size, 
     &       ' 0 360 arc closepath fill'

        end if

      end do
!
!  Label each node, but only if there aren't too many of them!
!
      if ( node_unmasked_num <= 200 ) then

        write ( eps_unit, '(a)' ) '%'
        write ( eps_unit, '(a)' ) '%  Label the nodes:'
        write ( eps_unit, '(a)' ) '%'
        write ( eps_unit, '(a)' ) ' 0.0000 0.0000 0.0000 setrgbcolor'
        write ( eps_unit, '(a)' ) 
     &    '/Times-Roman findfont 0.20 inch scalefont setfont'

        do node = 1, np

          if ( node_mask(node) ) then

            eps_x = int 
     &         ( ( node_x_max - xc(node)              ) *  61.0   
     &         + (            + xc(node) - node_x_min ) * 551.0 ) 
     &         / scale

            eps_y = int 
     &         ( ( node_y_max - yc(node)              ) * 151.0   
     &         + (              yc(node) - node_y_min ) * 641.0 ) 
     &         / scale
    
            write ( eps_unit, '(i4,2x,i4,a,i4,a)' ) eps_x, eps_y+5, 
     &         ' moveto (', node, ') show'

          end if

        end do

      end if

      write ( eps_unit, '(a)' ) '%'
      write ( eps_unit, '(a)' ) 'restore showpage'
      write ( eps_unit, '(a)' ) '%'
      write ( eps_unit, '(a)' ) '% End of page'
      write ( eps_unit, '(a)' ) '%'
      write ( eps_unit, '(a)' ) '%%Trailer'
      write ( eps_unit, '(a)' ) '%%EOF'

      close ( unit = eps_unit )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'NODE_EPS:'
      write ( *, '(a)' ) '  An encapsulated PostScript file was created'
      write ( *, '(a)' ) '  containing an image of the nodes.'
      write ( *, '(a)' ) '  The file name is'
      write ( *, '(4x,a)' ) node_file_name
    
      return
      end
      subroutine nstoke ( xc, yc, area, xm, ym,
     &   a, f, g, uold, visc_inv, tolns, xlngth, ylngth,
     &   node, indx, ipivot, mrow1,
     &   nlband, nuband, nband, nrow1, ncol1,
     &   nelemn, np, nnodes, nuk, nquad, neqn,
     &   nsteps, nsim, maxnd, maxel, rdel, alpha, beta )
!
!**********************************************************************
!
!! NSTOKE solves the Navier-Stokes equations using Taylor-Hood elements.
!
!
!  Modified:
!
!    08 August 2004
!
!  Author:
!
!    Hyung-Chun Lee,
!    Department of Mathematics,
!    Ajou University, Korea
!
!  Parameters:
!
!    Input, double precision XC(NP), YC(NP), the coordinates of the nodes.
!
!    Input, double precision AREA(NELEMN), the area of each element.
!
!    Input, double precision XM(MAXEL,3), YM(MAXEL,3), the coordinates
!    of quadrature points in each element.
!
!    Workspace, double precision A(MINUN,MAXUN), storage for the banded system matrix.
!
!    Input, integer NP, the number of nodes.
!
!    Input, integer NSIM, the number of simple iterations performed.
!
      implicit none
!
      integer maxel
      integer maxnd
      integer mrow1
      integer np

      double precision a(mrow1,*)
      double precision aij
      double precision alpha
      double precision ar
      double precision area(*)
      double precision arr
      double precision bb
      double precision bbb
      double precision bbbl
      double precision bbl
      double precision bbx
      double precision bby
      double precision beta(6)
      double precision bx
      double precision by
      double precision csim
      double precision det
      double precision diff
      double precision etax
      double precision etay
      double precision f(*)
      double precision g(*)
      integer i
      integer indx(maxnd,3)
      integer info
      integer ip
      integer ipivot(*)
      integer ipp
      integer iq
      integer iqq
      integer iquad
      integer it
      integer iter
      integer iuk
      integer iukk
      integer iun
      integer iuse
      integer j
      integer job
      integer kk
      integer nband
      integer ncol1
      integer nelemn
      integer neqn
      integer niter
      integer nlband
      integer nnodes
      integer node(maxel,*)
      integer nquad
      integer nrow1
      integer nsim
      integer nsteps
      integer nuband
      integer nuk
      double precision rdel
      double precision refbsp
      double precision tbbx
      double precision tbby
      double precision tbx
      double precision tby
      double precision tolns
      double precision ubc
      double precision ubdry
      double precision un(2)
      integer unk_u
      integer unk_v
      double precision unx(2)
      double precision uny(2)
      double precision uold(*)
      double precision uold_qp
      double precision visc
      double precision visc_inv
      double precision vold_qp
      double precision x
      double precision xc(np)
      double precision xix
      double precision xiy
      double precision xlngth
      double precision xm(maxel,*)
      double precision y
      double precision yc(np)
      double precision ylngth
      double precision ym(maxel,*)
!
      visc = 1.d0 / visc_inv
!
!  matrix assembly triangle by triangle
!
      do iter = 1, nsteps

        niter = iter

        if ( iter .lt. nsim ) then
          csim = 0.0D+00
        else
          csim = 1.0D+00
        end if

        do i = 1, nrow1
          do j = 1, ncol1
            a(i,j) = 0.d0
          end do
        end do

        do it = 1, nelemn 

          arr = area(it) / 3.d0

          do iquad = 1, nquad 

            y = ym(it,iquad)
            x = xm(it,iquad)
            call trans ( it, x, y, det, xix, xiy, etax, etay, xc,
     &          yc, node, maxel )
            ar = arr * det

            do kk = 1, 2
              un(kk) = 0.d0
              uny(kk) = 0.d0
              unx(kk) = 0.d0
            end do
            uold_qp = 0.0d0
            vold_qp = 0.0d0

            do iq = 1, nnodes

              call refqbf ( x, y, iq, bb, tbx, tby )
              bx = tbx * xix + tby * etax
              by = tbx * xiy + tby * etay
              ip = node(it,iq)

              do iuk = 1, 2

                iun = indx(ip,iuk)

                if ( 0 .lt. iun ) then
                  un(iuk)  = un(iuk)  + bb * g(iun)
                  unx(iuk) = unx(iuk) + bx * g(iun)
                  uny(iuk) = uny(iuk) + by * g(iun)
                  if ( iuk .eq. 1 ) uold_qp = uold_qp + bb * uold(iun)
                  if ( iuk .eq. 2 ) vold_qp = vold_qp + bb * uold(iun)
                else if ( iun .lt. 0 ) then
                  ubc = alpha * ubdry(iuk,ip,iun,beta,xc,yc)
                  un(iuk)  = un(iuk)  + bb * ubc
                  unx(iuk) = unx(iuk) + bx * ubc
                  uny(iuk) = uny(iuk) + by * ubc
                  if ( iuk .eq. 1 ) uold_qp = uold_qp + bb * ubc
                  if ( iuk .eq. 2 ) vold_qp = vold_qp + bb * ubc
                end if

              end do

            end do
                       
            do iq = 1, nnodes 

              ip=node(it,iq)
              call refqbf(x,y,iq,bb,tbx,tby)
              bx=tbx*xix+tby*etax
              by=tbx*xiy+tby*etay 

              if ( iq .le. 3 ) then
                bbl=refbsp(x,y,iq)
              end if

              do 210 iuk = 1, nuk

                i = indx(ip,iuk)

                if ( i .le. 0 ) go to 210

                if ( iuk .eq. 1 ) then

                  f(i) = f(i) + csim *
     &               ( (un(1)*unx(1)+un(2)*uny(1))*bb )*ar
     &               +rdel * uold_qp * bb*ar

                else if ( iuk .eq. 2 ) then

                  f(i) = f(i) + csim * 
     &               ( (un(1)*unx(2)+un(2)*uny(2))*bb )*ar
     &               +rdel * vold_qp * bb*ar

                end if

                do iqq = 1, nnodes

                  ipp=node(it,iqq)
                  call refqbf(x,y,iqq,bbb,tbbx,tbby)
                  bbx=tbbx*xix+tbby*etax
                  bby=tbbx*xiy+tbby*etay

                  if ( iqq .le. 3 ) then
                    bbbl = refbsp(x,y,iqq)
                  end if

                  do 190 iukk = 1, nuk

                    j=indx(ipp,iukk)

                    if ( j .eq. 0 ) go to 190
                    aij = 0.d0 
                    if ( i .eq. neqn ) go to 190

                    if ( iuk .eq. 1 ) then

                      if ( iukk .eq. 1 ) then
                        aij=visc*(by*bby+bx*bbx)
     &                   +(bbb*unx(1)*bb)*csim
     &                   +bb*bbx*un(1)
     &                   +bb*bby*un(2) + rdel*(bb*bbb)
                      else if ( iukk .eq. 2 ) then
                        aij=csim*(bb*bbb*uny(1))
                      else if ( iukk .eq. 3 ) then
                        aij=-bx*bbbl
                      end if

                    else if ( iuk .eq. 2 ) then

                      if ( iukk.eq.1 ) then
                        aij= csim*(bb*bbb*unx(2))
                      else if (iukk.eq.2) then
                        aij=(visc*(by*bby+bx*bbx)
     &                     +(bb*bbb*uny(2))*csim
     &                     +bb*bby*un(2)
     &                     +bb*bbx*un(1)) + rdel*(bb*bbb)
                      else if ( iukk .eq. 3 ) then
                        aij=-by*bbbl
                      end if

                    else if ( iuk .eq. 3 ) then

                      if ( iukk .eq. 1 ) then
                        aij=bbx*bbl
                      else if ( iukk .eq. 2 ) then
                        aij=bby*bbl
                      end if

                    end if
!
!  The coefficient AIJ is added to the matrix entry if J represents an unknown,
!  or is subtracted from the right hand side if J corresponds to a known "variable".
!
                  if ( 0 < j ) then
                    iuse=i-j+nband
                    a(iuse,j)=a(iuse,j)+aij*ar
                  else
                    f(i)=f(i)-ar*alpha*ubdry(iukk,ipp,j,beta,xc,yc)*aij 
                  end if

 190              continue
                end do
 210          continue
            end do
          end do
        end do
!
!  Replace last equation to set a reference pressure to 0.
!
        f(neqn)=0.d0
        do j = neqn-nlband, neqn-1
          i=neqn-j+nband
          a(i,j)=0.d0
        end do 
        a(nband,neqn)=1.0d0
!
!  Solve the system.
!
        job=0
        call dgbfa(a,mrow1,neqn,nlband,nuband,ipivot,info)

        call dgbsl(a,mrow1,neqn,nlband,nuband,ipivot,f,job)
!
!  Check for convergence
!
        diff = 0.d0

        do i = 1, np

          unk_u = indx(i,1)

          if ( 0 < unk_u ) then
            diff = diff + ( g(unk_u) - f(unk_u) )**2
          end if

          unk_v = indx(i,2)

          if ( 0 < unk_v ) then
            diff = diff + ( g(unk_v) - f(unk_v) )**2
          end if

        end do

        diff = sqrt ( diff )

        write ( 3, * ) '  Iteration ', iter, '  Difference is ', diff

        if ( diff .le. tolns ) then
          return
        end if

        do i = 1, neqn
          g(i) = f(i)
          f(i) = 0.d0
        end do

      end do

      return
      end
      subroutine p_write ( p_file_name, node_max, node_num, neqn, 
     &  c, indx, node_x, node_y )
!
!***********************************************************************
!
!! P_WRITE writes the pressures values for a given timestep to a file.
!
!
!  Discussion:
!
!    We HOPE that the INDX vector is set up so that INDX(NODE,3) = 0
!    exactly when NODE is a midside node.  We ASSUME that all corner
!    nodes have an INDX value that is 
!
!    * positive, if the value must be solved for;
!    * negative, if the value is prescribed in some other way.
!
!    I need to write U, V and P to files so I can plot and analyze the data.
!    The numerical scheme only computes P at some nodes.  The treatment of
!    boundary conditions may also result in special handling for the values
!    of velocity or pressure at other nodes.  
!
!    I would like to use the same simple format to INPUT solution data to
!    be used as startup, and to OUTPUT for analysis and plotting.  The
!    most sensible thing seems to be to determine the variables at all
!    nodes, either by interpolation, boundary specification, or the
!    finite element coefficients.  This code does not achieve that goal!
!
!  Modified:
!
!    21 April 2004
!
!  Parameters:
!
!    Input, character P_FILE_NAME, the name of the output file.
!
!    Input, integer NODE_MAX, the maximum number of nodes.
!
!    Input, integer NODE_NUM, the number of nodes.
!
!    Input, integer NEQN, the number of unknowns.
!
!    Input, double precision C(NEQN), the finite element coefficients.
!
!    Input, integer INDX(MAXND,NUK), lists the indices of the U, V, and P
!    variables associated with the node.
!
!    Input, double precision NODE_X(NODE_NUM), NODE_Y(NODE_NUM), the coordinates of the nodes.
!
      implicit none
!
      integer node_max
      integer node_num
      integer neqn
!
      double precision c(neqn)
      integer i
      integer indx(node_max,3)
      integer node
      double precision node_x(node_num)
      double precision node_y(node_num)
      double precision p
      character ( len = * ) p_file_name
!
      open ( unit = 1, file = p_file_name, err = 10 )

      do node = 1, node_num

        i = indx(node,3)

        if ( 0 == i ) then

        else if ( 0 < i ) then
          p = c(i)
          write ( 1, '(g14.6)' ) p
        else if ( i < 0 ) then
          p = 0.0D+00
          write ( 1, '(g14.6)' ) p
        end if

      end do

      close ( unit = 1 )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'P_WRITE:'
      write ( *, '(a)' ) '  Wrote the pressure output file:'
      write ( *, '(a)' ) ' '
      write (  *, '(4x,a)' ) p_file_name

      return

10    continue

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'P_WRITE - Fatal error!'
      write ( *, '(a)' ) '  Could not open the output file:'
      write ( *, '(a)' ) ' '
      write ( *, '(4x,a)' ) p_file_name
      stop
      end
      subroutine quad_a_set ( maxel, nelemn, nquad, area, xm, ym )
!
!***********************************************************************
!
!! QUAD_A_SET sets quadrature information for the assembly routine.
!
!
!  Modified:
!
!    20 April 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer MAXEL, the maximum number of elements.
!
!    Input, integer NELEMN, the number of elements.
!
!    Input, integer NQUAD, the number of quadrature points per element.
!
!    Output, double precision AREA(NELEMN), the nominal area of each element.
!
!    Output, double precision XM(MAXEL,3), YM(MAXEL,3), the coordinates of
!    the quadrature points.
!
      implicit none
!
      integer maxel
      integer nquad
!
      double precision area(maxel)
      integer element
      integer nelemn
      double precision xm(maxel,nquad)
      double precision ym(maxel,nquad)
!
      do element = 1, nelemn
        xm(element,1) = 0.5d0
        xm(element,2) = 0.5d0
        xm(element,3) = 0.0d0
      end do

      do element = 1, nelemn
        ym(element,1) = 0.0d0
        ym(element,2) = 0.5d0
        ym(element,3) = 0.5d0
      end do

      do element = 1, nelemn
        area(element) = 0.5d0
      end do

      return
      end
      function refbsp ( x, y, iq )

!***********************************************************************
!
!! REFBSP evaluates a linear basis functions on the reference triangle.
!
!  Diagram:
!
!    3
!    |\
!    | \
!    |  \
!    1---2
!
!  Modified:
!
!    18 June 2002
!
!  Parameters:
!
!    Input, double precision X, Y, the coordinates of a point in
!    the reference triangle.
!
!    Input, integer IQ, the index of a basis function in the reference triangle.
!
!    Output, double precision REFBSP, the value of the IQ-th basis function
!    at the point (X,Y) in the reference triangle.
!
      implicit none

      integer iq
      double precision refbsp
      double precision x
      double precision y

      if ( iq .eq. 1 ) then
        refbsp = 1.0D+00 - x - y
      else if ( iq .eq. 2 ) then
        refbsp = x
      else if ( iq .eq. 3 ) then
        refbsp = y
      else
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'REFBSP - Fatal error!'
        write ( *, '(a)' ) '  Illegal input value of IQ = ', iq
        stop
      end if

      return
      end
      subroutine refqbf ( x, y, in, bb, bx, by )
!
!************************************************************************
!
!! REFQBF evaluates quadratic basis functions on the reference triangle.
!
!
!  Modified:
!
!    18 June 2002
!
      implicit none

      double precision bb
      double precision bx
      double precision by
      integer in
      double precision x
      double precision y

      if ( in .eq. 1 ) then
        bb = (1.d0-x-y) * (1.d0-2.d0*x-2.d0*y)
        bx = -3.d0 + 4.d0 * x + 4.d0 * y
        by = -3.d0 + 4.d0 * x + 4.d0 * y
      else if ( in .eq. 2 ) then
        bb = x*(2.d0*x-1.d0)
        bx=4.d0*x-1.d0
        by=0.d0
      else if ( in .eq. 3 ) then
        bb=y*(2.d0*y-1.d0)
        bx=0.d0
        by=4.d0*y-1.d0
      else if ( in .eq. 4 ) then
        bb=4.d0*x*(1.d0-x-y)
        bx=4.d0*(1.d0-2.d0*x-y)
        by=-4.d0*x
      else if ( in .eq. 5 ) then
        bb=4.d0*x*y
        bx=4.d0*y
        by=4.d0*x
      else if ( in .eq. 6 ) then
        bb=4.d0*y*(1.d0-x-y)
        bx=-4.d0*y
        by=4.d0*(1.d0-x-2.d0*y)
      else
        bb = 0.0D+00
        bx = 0.0D+00
        by = 0.0D+00
      end if

      return
      end
      subroutine setgrd ( xc, yc, area, xm, ym, region_density_x,
     &  region_density_y, region_x, region_y,
     &  node, indx, nlband, nuband, nband, nelemn, np, 
     &  nnodes, nuk, nquad, neqn, maxnd, maxel )

!***********************************************************************
!
!! SETGRD sets up the grid for the problem.
!
!  Modified:
!
!    20 April 2004
!
!  Parameters:
!
!    Output, double precision XC(NP), YC(NP), the coordinates of the nodes.
!
!    Output, double precision AREA(NELEMN), the area of each element.
!
!    Output, double precision XM(MAXEL,NQUAD), YM(MAXEL,NQUAD), the coordinates
!    of quadrature points in each element.
!
!    Input, integer REGION_DENSITY_X(3), REGION_DENSITY_Y(3), specifies the
!    density of elements in the two coordinate directions.
!
!    Input, double precision REGION_X(4), REGION_Y(4), the coordinates of 
!    breakpoints that define 9 logical subrectangles.
!
!    Output, integer NODE(MAXEL,NNODES), the nodes that make up each element.
!
!    Output, integer INDX(MAXND,NUK), lists the indices of the U, V, and P
!    variables associated with the node.
!
!    Output, integer NLBAND, NUBAND, the lowwer and upper half bandwidths
!    for the finite element matrix.
!
!    Output, integer NBAND, the bandwidth for the finite element matrix.
!
!    Output, integer NELEMN, the number of elements.
!
!    Output, integer NP, the number of nodes.
!
!    Input, integer NNODES, the number of nodes per element.
!
!    Input, integer NUK, the maximum number of unknowns associated with one node.
!
!    Input, integer NQUAD, the number of quadrature points.
!
!    Output, integer NEQN, the total number of unknowns.
!
!    Input, integer MAXND, the maximum number of nodes.
!
!    Input, integer MAXEL, the maximum number of elements.
!
      implicit none

      integer maxel
      integer maxnd
      integer nnodes
      integer nquad
      integer nuk

      double precision area(maxel)
      integer dof_count
      integer indx(maxnd,nuk)
      integer nband
      integer nelemn
      integer neqn
      integer nlband
      integer node(maxel,nnodes)
      integer np
      integer nuband
      integer region_density_x(3)
      integer region_density_y(3)
      double precision region_x(4)
      double precision region_y(4)
      double precision xc(maxnd)
      double precision xlngth
      double precision xm(maxel,nquad)
      double precision yc(maxnd)
      double precision ylngth
      double precision ym(maxel,nquad)
!
!  Determine NP, the number of nodes.
!
      call hcell_node_count ( region_density_x, region_density_y, 
     &  np )
!
!  Assign coordinates to the nodes, XC and YC.
!
      call hcell_node_xy ( region_density_x, region_density_y, 
     &  np, region_x, region_y, xc, yc )
!
!  Determine NELEMN, the number of elements.
!
      call hcell_element_count ( region_density_x, 
     &  region_density_y, nelemn )
!
!  Assign nodes to elements in NODE.
!
      call hcell_element_node ( region_density_x, region_density_y, 
     &  maxel, nelemn, nnodes, node )
!
!  Determine NEQN, the number of degrees of freedom.
!  Assign degrees of freedom in INDX.
!
      call hcell_dof_set ( region_density_x, 
     &  region_density_y, maxnd, np, indx, neqn )
!
!  For the assembly routine, determine the quadrature data
!  NQUAD, AREA, XM and YM.
!
      call quad_a_set ( maxel, nelemn, nquad, area, xm, ym )
!
!  Compute the bandwidths NLBAND, NUBAND and NBAND.
!
      call element_node_bandwidth ( maxnd, maxel, nnodes, 
     &   nelemn, node, neqn, np, indx, nlband, nuband, nband )

      return
      end
      subroutine trans ( it, xq, yq, det, pj11, pj21, pj12, pj22,
     &  xc, yc, node, maxel )
!
! ***********************************************************************
!
!! TRANS transforms data between the reference and physical elements.
!
!
!  Modified:
!
!    18 June 2002
!
      implicit none
!
      integer maxel
!
      double precision det
      double precision f1x
      double precision f1y
      double precision f2x
      double precision f2y
      integer i1
      integer i2
      integer i3
      integer i4
      integer i5
      integer i6
      integer it
      integer node(maxel,*)
      double precision pj11
      double precision pj12
      double precision pj21
      double precision pj22
      double precision x1
      double precision x2
      double precision x3
      double precision x4
      double precision x5
      double precision x6
      double precision xc(*)
      double precision xq
      double precision y1
      double precision y2
      double precision y3
      double precision y4
      double precision y5
      double precision y6
      double precision yc(*)
      double precision yq
!
      i1 = node(it,1)
      i2 = node(it,2)
      i3 = node(it,3)
      i4 = node(it,4)
      i5 = node(it,5)
      i6 = node(it,6)

      x1 = xc(i1)
      y1 = yc(i1)
      x2 = xc(i2)
      y2 = yc(i2)
      x3 = xc(i3)
      y3 = yc(i3)
      x4 = xc(i4)
      y4 = yc(i4)
      x5 = xc(i5)
      y5 = yc(i5)
      x6 = xc(i6)
      y6 = yc(i6)
!
!  Compute partial derivatives at point (xq,yq)
!
      f1x=x1*(-3.d0+4.d0*xq+4.d0*yq)
     &       +x2*(4.d0*xq-1.d0)
     &       +x4*4.d0*(1.d0-2.d0*xq-yq)
     &       +x5*4.d0*yq + x6*4.d0*(-yq)

      f1y=x1*(-3.d0+4.d0*xq+4.d0*yq)
     &       +x3*(4.d0*yq-1.d0)
     &       +x4*4.d0*(-xq) + x5*4.d0*xq
     &       +x6*4.d0*(1.d0-xq-2.d0*yq)

      f2x=y1*(-3.d0+4.d0*xq+4.d0*yq)
     &       +y2*(4.d0*xq-1.d0)
     &       +y4*4.d0*(1.d0-2.d0*xq-yq)
     &       +y5*4.d0*yq + y6*4.d0*(-yq)

      f2y=y1*(-3.d0+4.d0*xq+4.d0*yq)
     &       +y3*(4.d0*yq-1.d0)
     &       +y4*4.d0*(-xq) + y5*4.d0*xq
     &       +y6*4.d0*(1.d0-xq-2.d0*yq)
!
!  Compute the determinant of the transformation evaluated at point (xq,yq)
!
      det=f1x*f2y-f1y*f2x
!
!  Compute j11, j22, j21, j22
!
      pj11 =  f2y / det
      pj12 = -f2x / det
      pj21 = -f1y / det
      pj22 =  f1x / det

      det = dabs ( det )

      return
      end
      function ubdry ( iuk, ip, iun, beta, xc, yc )
!
!***********************************************************************
!
!! UBDRY evaluates boundary conditions at a node.
!
!
!  Modified:
!
!    21 May 2004
!
!  Author:
!
!    Max Gunzburger
!
!  Parameters:
!
!    Input, integer IUK, indicates the class of variable being set.
!    * 1 = horizontal velocity;
!    * 2 = vertical velocity;
!    * 3 = pressure.
!
!    Input, integer IP, the index of the node.
!
!    Input, integer IUN, the value of INDX(IP,IUK), that is, the global
!    unknown number that was assigned to this "degree of freedom".
!
!    Input, double precision BETA(6), controls the strength of inlet and
!    outlet flows.
!
!    Input, double precision XC(*), YC(*), the coordinates of nodes.
!
      implicit none
!
      double precision beta(6)
      integer ip
      integer iuk
      integer iun
      double precision ubdry
      double precision xc(*)
      double precision yc(*)
!
      ubdry = 0.0D+00  
!
!  A nonzero horizontal velocity is specified at vent groups 1, 2, 3 and 4.
!
      if ( iuk .eq. 1 ) then
!
!  1: Zone 1 inlets
!
        if ( iun. eq. -1 ) then

          ubdry = beta(1) * 0.12D+00 * 4.0D+00 *
     &       ( yc(ip) - 8.0D+00 ) * ( 9.0D+00 - yc(ip) ) 
!
!  2: Zone 2 inlets.
!
        else if ( iun. eq. -2 ) then

          ubdry = - beta(2) * 0.125D+00 * 4.0D+00 *
     &       ( yc(ip) - 8.0D+00 ) * ( 9.0D+00 - yc(ip) ) 
!
!  3: Zone 3 inlets.
!
        else if ( iun. eq. -3 ) then

          ubdry = beta(3) * 0.11D+00 * 4.0D+00 *
     &       ( yc(ip) - 2.0D+00 ) * ( 3.0D+00 - yc(ip) )
!
!  4: Zone 4 inlets.
!
        else if ( iun. eq. -4 ) then

          ubdry = - beta(4) * 0.088D+00 * 4.0D+00 *
     &       ( yc(ip) - 2.0D+00 ) * ( 3.0D+00 - yc(ip) )

        else

          ubdry = 0.0D+00

        end if
!
!  A nonzero vertical velocity is specified at vent groups 1, 2, 3, 4, 5, 6.
!
      else if ( iuk .eq. 2 ) then
!
!  1: Zone 1 inlets.
!
        if ( iun. eq. -1 ) then

          if ( 10.0 .le. xc(ip) .and. xc(ip) .le. 11.0 ) then
            ubdry = beta(1) * 0.12D+00 * 4.0D+00 *
     &         ( xc(ip) - 10.0D+00 ) * ( 11.0D+00 - xc(ip) )
          else if ( 22.0 .le. xc(ip) .and. xc(ip) .le. 23.0 ) then
            ubdry = beta(1) * 0.12D+00 * 4.0D+00 *
     &         ( xc(ip) - 22.0D+00 ) * ( 23.0D+00 - xc(ip) )
          else if ( 34.0 .le. xc(ip) .and. xc(ip) .le. 35.0 ) then
            ubdry = beta(1) * 0.12D+00 * 4.0D+00 *
     &         ( xc(ip) - 34.0D+00 ) * ( 35.0D+00 - xc(ip) )
          end if
!
!  2: Zone 2 inlets.
!
        else if ( iun. eq. -2 ) then

          if ( 70.0 .le. xc(ip) .and. xc(ip) .le. 71.0 ) then
            ubdry = beta(2) * 0.125D+00 * 4.0D+00 *
     &         ( xc(ip) - 70.0D+00 ) * ( 71.0D+00 - xc(ip) )
          else if ( 82.0 .le. xc(ip) .and. xc(ip) .le. 83.0 ) then
            ubdry = beta(2) * 0.125D+00 * 4.0D+00 *
     &         ( xc(ip) - 82.0D+00 ) * ( 83.0D+00 - xc(ip) )
          else if ( 94.0 .le. xc(ip) .and. xc(ip) .le. 95.0 ) then
            ubdry = beta(2) * 0.125D+00 * 4.0D+00 *
     &         ( xc(ip) - 94.0D+00 ) * ( 95.0D+00 - xc(ip) )
          end if
!
!  3: Zone 3 inlets.
!
        else if ( iun. eq. -3 ) then

          if ( 10.0 .le. xc(ip) .and. xc(ip) .le. 11.0 ) then
            ubdry = - beta(3) * 0.11D+00 * 4.0D+00 *
     &         ( xc(ip) - 10.0D+00 ) * ( 11.0D+00 - xc(ip) )
          else if ( 22.0 .le. xc(ip) .and. xc(ip) .le. 23.0 ) then
            ubdry = - beta(3) * 0.11D+00 * 4.0D+00 *
     &         ( xc(ip) - 22.0D+00 ) * ( 23.0D+00 - xc(ip) )
          else if ( 34.0 .le. xc(ip) .and. xc(ip) .le. 35.0 ) then
            ubdry = - beta(3) * 0.11D+00 * 4.0D+00 *
     &         ( xc(ip) - 34.0D+00 ) * ( 35.0D+00 - xc(ip) )
          end if
!
!  4: Zone 4 inlets.
!
        else if ( iun. eq. -4 ) then

          if ( 70.0D+00 .le. xc(ip) .and. xc(ip) .le. 71.0D+00 ) then
            ubdry = - beta(4) * 0.088D+00 * 4.0D+00 *
     &         ( xc(ip) - 70.0D+00 ) * ( 71.0D+00 - xc(ip) )
          else if ( 82.0 .le. xc(ip) .and. xc(ip) .le. 83.0 ) then
            ubdry = - beta(4) * 0.088D+00 * 4.0D+00 *
     &         ( xc(ip) - 82.0D+00 ) * ( 83.0D+00 - xc(ip) )
          else if ( 94.0 .le. xc(ip) .and. xc(ip) .le. 95.0 ) then
            ubdry = - beta(4) * 0.088D+00 * 4.0D+00 *
     &         ( xc(ip) - 94.0D+00 ) * ( 95.0D+00 - xc(ip) )
          end if
!
!  5: Zone 1 outlets.
!
        else if ( iun. eq. -5 ) then

          if ( 10.0 .le. xc(ip) .and. xc(ip) .le. 11.0 ) then
            ubdry = beta(5) * 0.128D+00 * 4.0D+00 *
     &         ( xc(ip) - 10.0D+00 ) * ( 11.0D+00 - xc(ip) )
          else if ( 22.0 .le. xc(ip) .and. xc(ip) .le. 23.0 ) then
            ubdry = beta(5) * 0.128D+00 * 4.0D+00 *
     &         ( xc(ip) - 22.0D+00 ) * ( 23.0D+00 - xc(ip) )
          else if ( 34.0 .le. xc(ip) .and. xc(ip) .le. 35.0 ) then
            ubdry = beta(5) * 0.128D+00 * 4.0D+00 *
     &         ( xc(ip) - 34.0D+00 ) * ( 35.0D+00 - xc(ip) )
          end if
!
!  6: Zone 1 outlets.
!
        else if ( iun. eq. -6 ) then

          if ( 70.0 .le. xc(ip) .and. xc(ip) .le. 71.0 ) then
            ubdry = beta(6) * 0.224D+00 * 4.0D+00 *
     &         ( xc(ip) - 70.0D+00 ) * ( 71.0D+00 - xc(ip) )
          else if ( 82.0 .le. xc(ip) .and. xc(ip) .le. 83.0 ) then
            ubdry = beta(6) * 0.224D+00 * 4.0D+00 *
     &         ( xc(ip) - 82.0D+00 ) * ( 83.0D+00 - xc(ip) )
          else if ( 94.0 .le. xc(ip) .and. xc(ip) .le. 95.0 ) then
            ubdry = beta(6) * 0.224D+00 * 4.0D+00 *
     &         ( xc(ip) - 94.0D+00 ) * ( 95.0D+00 - xc(ip) )
          end if

        else

          ubdry = 0.0D+00

        end if

      end if

      return
      end
      subroutine uv_read ( uv_file_name, neqn, uold )
!
!***********************************************************************
!
!! UV_READ returns an initial value for the solution coefficient vector.
!
!
!  Modified:
!
!    20 April 2004
!
!  Parameters:
!
!    Input, character ( len = * ) UV_FILE_NAME, the file that contains the 
!    initial values.
!
!    Input, integer NEQN, the number of coefficients.
!
!    Output, double precision UOLD(NEQN), the initial value for the
!    solution coefficients.
!
      implicit none
!
      integer neqn
!
      character ( len = * ) uv_file_name
      integer i
      double precision uold(neqn)
!
!     open ( unit = 1, file = uv_file_name, err = 10 )

      do i = 1, neqn
!        read ( 1, * ) uold(i)
         uold(i)=0.0d0
      end do

!     close ( unit = 1 )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'UV_READ:'
      write ( *, '(a)' ) '  Read the file containing'
      write ( *, '(a)' ) '  the initial solution values:'
      write ( *, '(a)' ) ' '
      write ( *, '(4x,a)' ) uv_file_name

      return

10    continue

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'UV_READ - Fatal error!'
      write ( *, '(a)' ) '  Could not open the file containing'
      write ( *, '(a)' ) '  the initial solution values:'
      write ( *, '(a)' ) ' '
      write ( *, '(4x,a)' ) uv_file_name

      stop
      end
      subroutine uv_write ( uv_file_name, node_max, node_num, neqn, 
     &  c, indx, node_x, node_y, beta )
!
!***********************************************************************
!
!! UV_WRITE writes the velocity values for a given timestep to a file.
!
!
!  Discussion:
!
!    This code has been modified so that the value of velocity
!    is correctly computed by using either the solved-for coefficients,
!    or else the specified boundary values.
!
!  Modified:
!
!    18 May 2004
!
!  Parameters:
!
!    Input, character UV_FILE_NAME, the name of the output file.
!
!    Input, integer NODE_MAX, the maximum number of nodes.
!
!    Input, integer NODE_NUM, the number of nodes.
!
!    Input, integer NEQN, the number of unknowns.
!
!    Input, double precision C(NEQN), the finite element coefficients.
!
!    Input, integer INDX(MAXND,NUK), lists the indices of the U, V, and P
!    variables associated with the node.
!
!    Input, double precision NODE_X(NODE_NUM), NODE_Y(NODE_NUM), the 
!    coordinates of the nodes.
!
!    Input, double precision BETA(6), controls the strength of the
!    flow at the various inlets and outlets.
!
      implicit none
!
      integer node_max
      integer node_num
      integer neqn
!
      double precision beta(6)
      double precision c(neqn)
      integer i
      integer indx(node_max,3)
      integer node
      double precision node_x(node_num)
      double precision node_y(node_num)
      double precision u
      double precision ubdry
      character ( len = * ) uv_file_name
      double precision v
!
      open ( unit = 1, file = uv_file_name, err = 10 )

      do node = 1, node_num

        i = indx(node,1)

        if ( 0 < i ) then
          u = c(i)
        else if ( i < 0 ) then
          u = ubdry ( 1, node, i, beta, node_x, node_y )
        else
          u = 0.0D+00
        end if

        i = indx(node,2)

        if ( 0 < i ) then
          v = c(i)
        else if ( i < 0 ) then
          v = ubdry ( 2, node, i, beta, node_x, node_y )
        else
          v = 0.0D+00
        end if

        write ( 1, '(2g14.6)' ) u, v

      end do

      close ( unit = 1 )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'UV_WRITE:'
      write ( *, '(a)' ) '  Wrote the output file:'
      write ( *, '(a)' ) ' '
      write (  *, '(4x,a)' ) uv_file_name

      return

10    continue

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'UV_WRITE - Fatal error!'
      write ( *, '(a)' ) '  Could not open the output file:'
      write ( *, '(a)' ) ' '
      write ( *, '(4x,a)' ) uv_file_name
      stop
      end
      subroutine uv_write_tec ( tec_file_name, node_max, node_num, 
     &  node_x, node_y, element_max, element_num, element_node, 
     &  indx, neqn, c, beta )
!
!***********************************************************************
!
!! UV_WRITE_TEC writes velocity information to a TECPLOT file.
!
!
!  Modified:
!
!    20 May 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) TEC_FILE_NAME, the file to be written.
!
!    Input, integer NODE_MAX, the maximum number of nodes.
!
!    Input, integer NODE_NUM, the number of nodes.
!
!    Input, double precision NODE_X(NODE_NUM), NODE_Y(NODE_NUM), the
!    coordinates of the nodes.
!
!    Input, integer ELEMENT_MAX, the maximum number of elements.
!
!    Input, integer ELEMENT_NUM, the number of elements.
!
!    Input, integer ELEMENT_NODE(ELEMENT_MAX,6), the nodes that make 
!    up each element.
!
!    Input, integer INDX(NODE_MAX,3), lists the indices of the U, V, and P
!    variables associated with the node.
!
!    Input, integer NEQN, the number of finite element coefficients.
!
!    Input, double precision C(NEQN), the finite element coefficients.
!
!    Input, double precision BETA(6), controls the strength of the
!    flow at the various inlets and outlets.
!
      implicit none
!
      integer element_max
      integer neqn
      integer node_max
!
      double precision beta(6)
      double precision c(neqn)
      integer element
      integer element_node(element_max,6)
      integer element_num
      character*(*) tec_file_name
      integer i
      integer indx(node_max,3)
      integer node
      integer node_num
      double precision node_x(node_max)
      double precision node_y(node_max)
      double precision u
      double precision ubdry
      double precision v
!
      open ( unit = 1, file = tec_file_name, form = 'formatted', 
     &  access = 'sequential', status = 'replace' )

      write ( 1, * ) 'Title="HCELL data"'
      write ( 1, * ) 'Variables="X","Y","U","V"'
      write ( 1, * ) 'Zone N=', node_num, ', E=', 4 * element_num, 
     &  ', F=FEPOINT, ET=TRIANGLE'
!
!  Write out the data at each node.
!
      do node = 1, node_num

        i = indx(node,1)

        if ( 0 < i ) then
          u = c(i)
        else if ( i < 0 ) then
          u = ubdry ( 1, node, i, beta, node_x, node_y )
        else
          u = 0.0D+00
        end if

        i = indx(node,2)

        if ( 0 < i ) then
          v = c(i)
        else if ( i < 0 ) then
          v = ubdry ( 2, node, i, beta, node_x, node_y )
        else
          v = 0.0D+00
        end if

        write ( 1, '(4g15.6)' ) node_x(node), node_y(node), u, v

      end do
!
!  Write out the data that defines the elements.
!  Each 6 node quadratic element must be described as 4 linear elements.
!
!  3
!  |\
!  | \
!  6--5
!  |\ |\
!  | \| \
!  1--4--2
!
      do element = 1, element_num
        write ( 1, '(3i6)' ) 
     &    element_node(element,1), element_node(element,4), 
     &    element_node(element,6)
        write ( 1, '(3i6)' ) 
     &    element_node(element,2), element_node(element,5), 
     &    element_node(element,4)
        write ( 1, '(3i6)' ) 
     &    element_node(element,3), element_node(element,6), 
     &    element_node(element,5)
        write ( 1, '(3i6)' ) 
     &    element_node(element,4), element_node(element,5), 
     &    element_node(element,6)
      end do

      close ( unit = 1 )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'UV_WRITE_TEC:'
      write ( *, '(a)' ) '  Wrote the output file:'
      write ( *, '(a)' ) ' '
      write (  *, '(4x,a)' ) tec_file_name

      return
      end
      subroutine xy3_write ( file_name, maxnd, node_num, indx, node_x, 
     &  node_y )
!
!***********************************************************************
!
!! XY3_WRITE writes the coordinates of nodes of 3-node triangles to a file.
!
!
!  Discussion:
!
!    We HOPE that the INDX vector is set up so that INDX(NODE,3) = 0
!    exactly when NODE is a midside node.  We ASSUME that all corner
!    nodes have an INDX value that is 
!
!    * positive, if the value must be solved for;
!    * negative, if the value is prescribed in some other way.
!
!  Modified:
!
!    21 April 2004
!
!  Parameters:
!
!    Input, character ( len = * ) FILE_NAME, the file to be written.
!
!    Input, integer MAXND, the maximum number of nodes.
!
!    Input, integer NODE_NUM, the number of nodes.
!
!    Input, integer INDX(MAXND,3), lists the indices of the U, V, and P
!    variables associated with the node.  It would be useful if the 
!    index associated with pressure also indicated whether there was
!    no pressure variable associated with the node, or that it was
!    prescribed.  This could be done by assigning INDX(NODE,3) = 0
!    for the midside nodes of the 6 node quadratic elements.
!
!    Input, double precision NODE_X(NODE_NUM), NODE_Y(NODE_NUM), the
!    coordinates of the nodes.
!
      implicit none
!
      integer maxnd
      integer node_num
!
      character ( len = * ) file_name
      integer indx(maxnd,3)
      integer node
      double precision node_x(node_num)
      double precision node_y(node_num)
!
      open ( unit = 1, file = file_name, err = 10 )

      do node = 1, node_num
        if ( 0 /= indx(node,3) ) then
          write ( 1, * ) node_x(node), node_y(node)
        end if
      end do

      close ( unit = 1 )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'XY3_WRITE:'
      write ( *, '(a)' ) '  The 3-node triangle nodal coordinates were'
      write ( *, '(a)' ) '  written to the output file:'
      write ( *, '(a)' ) ' '
      write ( *, '(4x,a)' ) file_name
      return

10    continue

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'XY3_WRITE - Fatal error!'
      write ( *, '(a)' ) '  Could not open the output file:'
      write ( *, '(a)' ) ' '
      write ( *, '(4x,a)' ) file_name

      stop
      end
      subroutine xy6_write ( file_name, node_num, node_x, node_y )
!
!***********************************************************************
!
!! XY6_WRITE writes the coordinates of nodes of 6-node triangles to a file.
!
!
!  Modified:
!
!    20 April 2004
!
!  Parameters:
!
!    Input, character ( len = * ) FILE_NAME, the file to be written.
!
!    Input, integer NODE_NUM, the number of nodes.
!
!    Input, double precision NODE_X(NODE_NUM), NODE_Y(NODE_NUM), the
!    coordinates of the nodes.
!
      implicit none
!
      integer node_num
!
      character ( len = * ) file_name
      integer node
      double precision node_x(node_num)
      double precision node_y(node_num)
!
      open ( unit = 1, file = file_name, err = 10 )

      do node = 1, node_num
        write ( 1, * ) node_x(node), node_y(node)
      end do

      close ( unit = 1 )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'XY6_WRITE:'
      write ( *, '(a)' ) '  The 6-node triangle nodal coordinates were'
      write ( *, '(a)' ) '  written to the output file:'
      write ( *, '(a)' ) ' '
      write ( *, '(4x,a)' ) file_name
      return

10    continue

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'XY6_WRITE - Fatal error!'
      write ( *, '(a)' ) '  Could not open the output file:'
      write ( *, '(a)' ) ' '
      write ( *, '(4x,a)' ) file_name

      stop
      end
