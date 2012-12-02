      subroutine bin_preprocess ( ndim, box_min, box_max, n, 
     &  cell_generator, nbin, bin_start, bin_last, bin_next )

c*********************************************************************72
c
cc BIN_PREPROCESS organizes the preprocessing step for bins.
c
c  Discussion:
c
c    This routine is required to set up the bin data for use in the
c    nearest neighbor algorithm.
c
c    There are separate sets of calls for the 2D and 3D cases.  Although the
c    algorithms are essentially identical, the declarations of data
c    structures, in particular, of NBIN, were a little tricky to 
c    program generically.
c
c    To help with understanding the bin structure, the routine will
c    print out a few statistics about the bins on the first call only.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    22 July 2008
c
c  Author:
c
c    John Burkardt, Max Gunzburger, Janet Peterson
c
c  Parameters:
c
c    Input, integer NDIM, the spatial dimension.
c
c    Input, double precision BOX_MIN(NDIM), BOX_MAX(NDIM), the coordinates
c    of the two extreme corners of the bounding box.
c
c    Input, integer N, the number of Voronoi cells.
c
c    Input, double precision CELL_GENERATOR(NDIM,N), the Voronoi
c    cell generators.
c
c    Input, integer NBIN(3) is the number of bins to use in each direction.
c    For 3D problems, set NBIN(3) = 1.
c    For efficiency, these values should be set in such a way that the bins
c    are nearly square or cubical.
c
c    Output, integer BIN_START(NBIN(1),NBIN(2),NBIN(3)), the index of the first
c    cell generator in the bin, or -1 if none.
c
c    Output, integer BIN_LAST(NBIN(1),NBIN(2),NBIN(3)), the index of the last
c    cell generator in the bin, or -1 if none.
c
c    Output, integer BIN_NEXT(N), the index of the next cell generator in
c    the bin containing this cell generator.
c
      implicit none

      integer n
      integer ndim

      integer nbin(ndim)

      integer bin_last( nbin(1),nbin(2),nbin(3) )
      integer bin_next(n)
      integer bin_start( nbin(1),nbin(2),nbin(3) )
      double precision box_min(ndim)
      double precision box_max(ndim)
      double precision cell_generator(ndim,n)
      integer empty
      logical first_call
      integer i
      integer j
      integer k
      integer nonempty
      integer total

      save first_call

      data first_call / .true. /

      if ( ndim .eq. 2 ) then

        call r82vec_bin_even3 ( n, cell_generator, nbin, box_min, 
     &    box_max, bin_start, bin_last, bin_next )

        call r82vec_binned_reorder2 ( n, cell_generator, nbin, 
     &    bin_start, bin_last, bin_next )

        call r82vec_binned_sort_a2 ( n, cell_generator, nbin, 
     &    bin_start, bin_last )

      else if ( ndim .eq. 3 ) then

        call r83vec_bin_even3 ( n, cell_generator, nbin, box_min, 
     &    box_max, bin_start, bin_last, bin_next )

        call r83vec_binned_reorder2 ( n, cell_generator, nbin, 
     &    bin_start, bin_last, bin_next )

        call r83vec_binned_sort_a2 ( n, cell_generator, nbin, 
     &    bin_start, bin_last )

      end if
c
c  On the first call only, prepare a brief statistical report.
c
      if ( first_call ) then

        first_call = .false.

        total = 1
        do i = 1, ndim
          total = total * nbin(i)
        end do

        nonempty = 0
        empty = 0
     
        do i = 1, nbin(1)
          do j = 1, nbin(2)
            do k = 1, nbin(3)
              if ( 0 .lt. bin_start(i,j,k) ) then
                nonempty = nonempty + 1
              else
                empty = empty + 1
              end if
            end do
          end do
        end do

        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'BIN_PREPROCESS:'
        write ( *, '(a)' ) ' '
        write ( *, '(a,i8)' ) '  Number of points =     ', n
        write ( *, '(a,i8)' ) '  Total number of bins = ', total
        write ( *, '(a,i8)' ) '  Number of empty bins = ', empty
        write ( *, '(a,i8)' ) '          nonempy bins = ', nonempty
        write ( *, '(a)' ) ' '
        write ( *, '(a,f5.1)' ) 
     &    '  Percentage nonempty bins =         ', 
     &    dble ( 100 * nonempty ) / dble ( total )
        write ( *, '(a,g14.6)' ) 
     &    '  Number of points per bin =         ', 
     &    dble ( n ) / dble ( total )
        write ( *, '(a,g14.6)' ) 
     &    '  Number of points per nonempy bin = ', 
     &    dble ( n ) / dble ( nonempty )

      end if

      return
      end
      subroutine bin_to_r8_even2 ( nbin, bin, a, b, cmin, cmax )

c*********************************************************************72
c
cc BIN_TO_R8_EVEN2 returns the limits for a given "bin" in [A,B].
c
c  Discussion:
c
c    The interval from A to B is divided into NBIN equal subintervals or bins.
c
c  Example:
c
c    NBIN = 5, A = 10, B = 20
c
c    BIN      CMIN  CMAX
c
c    1         10    12
c    2         12    14
c    3         14    16
c    4         16    18
c    5         18    20
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    22 July 2008
c
c  Author:
c
c    John Burkardt, Max Gunzburger, Janet Peterson
c
c  Parameters:
c
c    Input, integer NBIN, the number of bins.
c
c    Input, integer BIN, the index of the bin to be considered.
c    If BIN is less than 1, or greater than NBIN, the user will get what
c    the user deserves.
c
c    Input, double precision A, B, the lower and upper limits of the bin
c    interval.  While A is expected to be less than B, the code should
c    return useful results if A is actually greater than B.
c
c    Output, double precision CMIN, CMAX, the minimum and maximum limits
c    on the bin.
c
      implicit none

      double precision a
      double precision b
      integer bin
      double precision cmax
      double precision cmin
      integer nbin
      double precision r8_huge
c
c  Compute the bin limits.
c
      if ( bin .lt. 1 ) then
        cmin = - r8_huge ( )
        cmax = a
      else if ( bin .le. nbin ) then
        cmin = ( dble ( nbin - bin + 1 ) * a 
     &         + dble (        bin - 1 ) * b ) 
     &         / dble ( nbin           )
        cmax = ( dble ( nbin - bin     ) * a 
     &         + dble (        bin     ) * b ) 
     &         / dble ( nbin           )
      else if ( nbin .le. bin ) then
        cmin = b
        cmax = r8_huge ( )
      end if

      return
      end
      subroutine bin_to_r8vec_even3 ( ndim, nbin, bin, a, b, cmin, 
     &  cmax )

c*********************************************************************72
c
cc BIN_TO_R8VEC_EVEN3 returns the limits for a given R8VEC "bin" in [A,B].
c
c  Discussion:
c
c    The interval from A(I) to B(I) is divided into NBIN(I) equal
c    subintervals or bins.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    22 July 2008
c
c  Author:
c
c    John Burkardt, Max Gunzburger, Janet Peterson
c
c  Parameters:
c
c    Input, integer NDIM, the dimension of the space.
c
c    Input, integer NBIN(NDIM), the number of bins in each dimension.
c
c    Input, integer BIN(NDIM), the index of the bin to be considered.
c
c    Input, double precision A(NDIM), B(NDIM), the lower and upper limits
c    of the bin interval.  While A(I) is expected to be less than B(I), the
c    code should return useful results if A(I) is actually greater than B(I).
c
c    Output, double precision CMIN(NDIM), CMAX(NDIM), the minimum and maximum
c    limits on the bin.
c
      implicit none

      integer ndim

      double precision a(ndim)
      double precision b(ndim)
      integer bin(ndim)
      double precision cmax(ndim)
      double precision cmin(ndim)
      integer i
      integer nbin(ndim)

      do i = 1, ndim
        call bin_to_r8_even2 ( nbin(i), bin(i), a(i), b(i), cmin(i), 
     &    cmax(i) )
      end do

      return
      end
      subroutine cvt_iteration ( ndim, box_min, box_max, n, 
     &  cell_generator, ns_cvt, use_diatom, use_bins, dr, updates, 
     &  nbin, bin_start, bin_last, bin_next )

c*********************************************************************72
c
cc CVT_ITERATION takes one step of the CVT iteration.
c
c  Discussion:
c
c    The routine is given a set of points, called "generators", which
c    define a tessellation of the region into Voronoi cells.  Each point
c    defines a cell.  Each cell, in turn, has a centroid, but it is
c    unlikely that the centroid and the generator coincide.
c
c    Each time this CVT iteration is carried out, an attempt is made
c    to modify the generators in such a way that they are closer and
c    closer to being the centroids of the Voronoi cells they generate.
c
c    A large number of sample points are generated, and the nearest generator
c    is determined.  A count is kept of how many points were nearest to each
c    generator.  Once the sampling is completed, the location of all the
c    generators is adjusted.  This step should decrease the discrepancy
c    between the generators and the centroids.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    22 July 2008
c
c  Author:
c
c    John Burkardt, Max Gunzburger, Janet Peterson
c
c  Reference:
c
c    John Burkardt, Max Gunzburger, Janet Peterson and Rebecca Brannon,
c    User Manual and Supporting Information for Library of Codes
c    for Centroidal Voronoi Placement and Associated Zeroth,
c    First, and Second Moment Determination,
c    Sandia National Laboratories Technical Report SAND2002-0099,
c    February 2002.
c
c  Parameters:
c
c    Input, integer NDIM, the spatial dimension.
c
c    Input, double precision BOX_MIN(NDIM), BOX_MAX(NDIM), the coordinates
c    of the two extreme corners of the bounding box.
c
c    Input, integer N, the number of Voronoi cells.
c
c    Input/output, double precision CELL_GENERATOR(NDIM,N), the Voronoi
c    cell generators.  On output, these have been modified
c
c    Input, integer NS_CVT, the number of sample points per generator.
c
c    Input, logical USE_DIATOM, is TRUE if DIATOM is to be called
c    to determine whether a point lies in the physical region; if it is
c    FALSE than a much simplified routine is used.
c
c    Input, logical USE_BINS, is TRUE if the bounding box is to be divided
c    up into bins to speed up the nearest neighbor search;
c    FALSE if the nearest neighbor seach is to be done naively.
c
c    Input, double precision DR, a tolerance used by DIATOM when testing
c    whether a point is within, outside of, or on the boundary of the
c    physical region.
c
c    Input/output, integer UPDATES(N), counts the number of times a cell
c    center has been updated.  Before the first call, all the entries of
c    UPDATES should be set to 1.  After each iteration, UPDATES will be
c    incremented by 1 for each cell generator that was updated.  Normally,
c    all of them will be so updated.
c
c    Input, integer NBIN(3) is the number of bins to use in each direction.
c    For 3D problems, set NBIN(3) = 1.
c    For efficiency, these values should be set in such a way that the bins
c    are nearly square or cubical.
c
c    Input, integer BIN_START(NBIN(1),NBIN(2),NBIN(3)), the index of the first
c    cell generator in the bin, or -1 if none.
c
c    Input, integer BIN_LAST(NBIN(1),NBIN(2),NBIN(3)), the index of the last
c    cell generator in the bin, or -1 if none.
c
c    Input, integer BIN_NEXT(N), the index of the next cell generator in
c    the bin containing this cell generator.
c
      implicit none

      integer n
      integer nbin(3)
      integer ndim

      integer bin_last( nbin(1),nbin(2),nbin(3) )
      integer bin_next(n)
      integer bin_start( nbin(1),nbin(2),nbin(3) )
      double precision box_max(ndim)
      double precision box_min(ndim)
      double precision cell_generator(ndim,n)
      double precision cell_generator2(ndim,n)
      integer count(n)
      logical debug 
      parameter ( debug = .false. )
      double precision dr
      integer i
      integer j
      integer nearest
      integer ngen
      integer ns_cvt
      integer random_sampler
      parameter ( random_sampler = 0 )
      logical reset
      integer updates(n)
      logical use_bins
      logical use_diatom
      double precision x(ndim)

      do j = 1, n
        do i = 1, ndim
          cell_generator2(i,j) = 0.0D+00
        end do
      end do

      do i = 1, n
        count(i) = 0
      end do

      reset = .false.

      do j = 1, ns_cvt * n
c
c  Generate a sampling point X.
c
        call region_sampler ( ndim, box_min, box_max, dr, x, 
     &    random_sampler, reset, use_diatom, ngen )
c
c  Find the nearest cell generator.
c
        if ( use_bins ) then

          if ( ndim .eq. 2 ) then

            call points_nearest_point_bins3_2d ( n, cell_generator, 
     &        nbin, box_min, box_max, bin_start, bin_last, bin_next, 
     &        x, nearest )

          else if ( ndim .eq. 3 ) then

            call points_nearest_point_bins3_3d ( n, cell_generator, 
     &        nbin, box_min, box_max, bin_start, bin_last, bin_next, 
     &        x, nearest )

          end if

        else

          call find_closest ( ndim, x, n, cell_generator, nearest )

        end if
c
c  Add X to the averaging data for CELL_GENERATOR(*,NEAREST).
c
        do i = 1, ndim
          cell_generator2(i,nearest) = 
     &      cell_generator2(i,nearest) + x(i)
        end do

        count(nearest) = count(nearest) + 1

      end do
c
c  Compute the new generators.
c
      do j = 1, n

        if ( count(j) .ne. 0 ) then

          do i = 1, ndim
            cell_generator(i,j) = cell_generator2(i,j) 
     &        / dble ( count(j) )
          end do

          updates(j) = updates(j) + 1

        end if

      end do

      return
      end
      subroutine find_closest ( ndim, x, n, cell_generator, nearest )

c*********************************************************************72
c
cc FIND_CLOSEST finds the Voronoi cell generator closest to a point X.
c
c  Discussion:
c
c    This routine finds the closest Voronoi cell generator by checking every
c    one.  For problems with many cells, this process can take the bulk
c    of the CPU time.  Other approaches, which group the cell generators into
c    bins, can run faster by a large factor.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    22 July 2008
c
c  Author:
c
c    John Burkardt, Max Gunzburger, Janet Peterson
c
c  Reference:
c
c    John Burkardt, Max Gunzburger, Janet Peterson and Rebecca Brannon,
c    User Manual and Supporting Information for Library of Codes
c    for Centroidal Voronoi Placement and Associated Zeroth,
c    First, and Second Moment Determination,
c    Sandia National Laboratories Technical Report SAND2002-0099,
c    February 2002.
c
c  Parameters:
c
c    Input, integer NDIM, the spatial dimension.
c
c    Input, double precision X(NDIM), the point to be checked.
c
c    Input, integer N, the number of cell generators.
c
c    Input, double precision CELL_GENERATOR(NDIM,N), the cell generators.
c
c    Output, integer NEAREST, the index of the nearest cell generators.
c
      implicit none

      integer n
      integer ndim

      double precision cell_generator(ndim,n)
      double precision distance
      double precision dist_sq
      integer i
      integer j
      integer nearest
      double precision r8_huge
      double precision x(ndim)

      nearest = 0
      distance = r8_huge ( )

      do j = 1, n

        dist_sq = 0.0D+00
        do i = 1, ndim
          dist_sq = dist_sq + ( cell_generator(i,j) - x(i) )**2
        end do

        if ( dist_sq .lt. distance ) then
          distance = dist_sq
          nearest = j
        end if

      end do

      distance = sqrt ( distance )

      return
      end
      subroutine generator_init ( ndim, box_min, box_max, n, 
     &  cell_generator, use_diatom, dr, random_generator )

c*********************************************************************72
c
cc GENERATOR_INIT initializes the Voronoi cell generators.
c
c  Discussion:
c
c    The points initialized here will be used to generate a tessellation
c    of the region into Voronoi cells.  Each generator point defines a
c    cell.  The CVT algorithm will try to modify these initial generators
c    in such a way that they are also the centroids of the cells they generate.
c
c    It is probably better to use Halton points for the centroids than
c    uniform random values, in the sense that the algorithm will probably
c    converge more quickly.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    22 July 2008
c
c  Author:
c
c    John Burkardt, Max Gunzburger, Janet Peterson
c
c  Reference:
c
c    John Burkardt, Max Gunzburger, Janet Peterson and Rebecca Brannon,
c    User Manual and Supporting Information for Library of Codes
c    for Centroidal Voronoi Placement and Associated Zeroth,
c    First, and Second Moment Determination,
c    Sandia National Laboratories Technical Report SAND2002-0099,
c    February 2002.
c
c  Parameters:
c
c    Input, integer NDIM, the spatial dimension.
c
c    Input, double precision BOX_MIN(NDIM), BOX_MAX(NDIM), the coordinates
c    of the two extreme corners of the bounding box.
c
c    Input, integer N, the number of Voronoi cells.
c
c    Output, double precision CELL_GENERATOR(NDIM,N), the Voronoi cell
c    generators.
c
c    Input, logical USE_DIATOM, is TRUE if DIATOM is to be called
c    to determine whether a point lies in the physical region; if it is
c    FALSE than a much simplified routine is used.
c
c    Input, double precision DR, a tolerance used by DIATOM when testing
c    whether a point is within, outside of, or on the boundary of the
c    physical region.
c
c    Input, integer RANDOM_GENERATOR, specifies how the
c    Voronoi cell generators are to be initialized.
c    0, use the F90 RANDOM_NUMBER routine;
c    1, use the Halton sequence.
c
      implicit none

      integer n
      integer ndim

      double precision box_max(ndim)
      double precision box_min(ndim)
      double precision cell_generator(ndim,n)
      double precision dr
      integer i
      integer ngen
      integer random_generator
      logical reset
      logical use_diatom

      reset = .true.

      do i = 1, n
        call region_sampler ( ndim, box_min, box_max, dr, 
     &    cell_generator(1,i), random_generator, reset, use_diatom, 
     &    ngen )
      end do

      return
      end
      subroutine i4_to_halton_vector ( seed, base, ndim, r )

c*********************************************************************72
c
cc I4_TO_HALTON_VECTOR computes an element of a vector Halton sequence.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    22 July 2008
c
c  Author:
c
c    John Burkardt, Max Gunzburger, Janet Peterson
c
c  Reference:
c
c    John Burkardt, Max Gunzburger, Janet Peterson and Rebecca Brannon,
c    User Manual and Supporting Information for Library of Codes
c    for Centroidal Voronoi Placement and Associated Zeroth,
c    First, and Second Moment Determination,
c    Sandia National Laboratories Technical Report SAND2002-0099,
c    February 2002.
c
c    John Halton,
c    On the efficiency of certain quasi-random sequences of points
c    in evaluating multi-dimensional integrals,
c    Numerische Mathematik,
c    Volume 2, pages 84-90.
c
c  Parameters:
c
c    Input, integer SEED, the index of the desired element.
c    Only the absolute value of SEED is considered.  SEED = 0 is allowed,
c    and returns R = 0.
c
c    Input, integer BASE(NDIM), the Halton bases, which should be
c    distinct prime numbers.  This routine only checks that each base
c    is greater than 1.
c
c    Input, integer NDIM, the dimension of the sequence.
c
c    Output, double precision R(NDIM), the SEED-th element of the Halton
c    sequence for the given bases.
c
      implicit none

      integer ndim

      logical all_zero
      integer base(ndim)
      double precision base_inv(ndim)
      integer digit(ndim)
      integer i
      double precision r(ndim)
      integer seed
      integer seed2(ndim)

      do i = 1, ndim
        seed2(i) = abs ( seed )
      end do

      do i = 1, ndim
        r(i) = 0.0D+00
      end do

      do i = 1, ndim

        if ( base(i) .le. 1 ) then
          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) 'I4_TO_HALTON_VECTOR - Fatal error!'
          write ( *, '(a)' ) '  An input base BASE is .LE. 1!'
          write ( *, '(2i8)' ) i, base(i)
          stop
        end if

      end do

      do i = 1, ndim
        base_inv(i) = 1.0D+00 / dble ( base(i) )
      end do

10    continue

      all_zero = .true.

      do i = 1, ndim
        if ( seed2(i) .ne. 0 ) then
          all_zero = .false.
        end if
      end do
       
      if ( .not. all_zero ) then
        do i = 1, ndim
          digit(i) = mod ( seed2(i), base(i) )
          r(i) = r(i) + dble ( digit(i) ) * base_inv(i)
          base_inv(i) = base_inv(i) / dble ( base(i) )
          seed2(i) = seed2(i) / base(i)
        end do
        go to 10
      end if

      return
      end
      subroutine index_box2_next_2d ( n1, n2, ic, jc, i, j, more )

c*********************************************************************72
c
cc INDEX_BOX2_NEXT_2D produces indices on the surface of a box in 2D.
c
c  Discussion:
c
c    The box has center at (IC,JC), and has half-widths N1 and N2.
c    The indices are exactly those which are between (IC-N1,JC-N2) and
c    (IC+N1,JC+N2) with the property that at least one of I and J
c    is an "extreme" value.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    22 July 2008
c
c  Author:
c
c    John Burkardt, Max Gunzburger, Janet Peterson
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
      subroutine index_box2_next_3d ( n1, n2, n3, ic, jc, kc, i, j, k, 
     &  more )

c*********************************************************************72
c
cc INDEX_BOX2_NEXT_3D produces indices on the surface of a box in 3D.
c
c  Discussion:
c
c    The box has a central cell of (IC,JC,KC), with a half widths of
c    (N1,N2,N3).  The indices are exactly those between (IC-N1,JC-N2,KC-N3)
c    and (IC+N1,JC+N2,KC+N3) with the property that at least one of I, J,
c    and K is an "extreme" value.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    22 July 2008
c
c  Author:
c
c    John Burkardt, Max Gunzburger, Janet Peterson
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

      if ( i .eq. ic + n1 .and. j .eq. jc + n2 .and. 
     &  k .eq. kc + n3 ) then
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
      subroutine points_nearest_point_bins3_2d ( nset, pset, nbin, 
     &  bin_min, bin_max, bin_start, bin_last, bin_next, ptest, 
     &  i_min )

c*********************************************************************72
c
cc POINTS_NEAREST_POINT_BINS3_2D finds the nearest point to a given point in 2D.
c
c  Discussion:
c
c    This code differs from POINTS_NEAREST_POINTS_BINS_2D by allowing the
c    user to specify the number of bins in each dimension.
c
c    A set of NSET points with coordinates PSET is assumed to lie within a
c    rectangle.  The limits of this rectangle are given in BIN_MIN and BIN_MAX.
c    The rectangle is divided up into an NBIN(1) by NBIN(2) regular grid of
c    cells.
c
c    The cells are ordered lexicographically, as suggested by the following
c    diagram for NBIN = (/ 5, 4 /)
c
c     *----*----*----*----*----*
c     |    |    |    |    |    |
c     | 16 | 17 | 18 | 19 | 20 |
c     |    |    |    |    |    |
c     *----*----*----*----*----*
c     |    |    |    |    |    |
c     | 11 | 12 | 13 | 14 | 15 |
c     |    |    |    |    |    |
c     *----*----*----*----*----*
c     |    |    |    |    |    |
c     |  6 |  7 |  8 |  9 | 10 |
c     |    |    |    |    |    |
c     *----*----*----*----*----*
c     |    |    |    |    |    |
c     |  1 |  2 |  3 |  4 |  5 |
c     |    |    |    |    |    |
c     *----*----*----*----*----*
c
c    The points in PSET are ordered by cell, and within a cell they
c    are ordered lexicographically.  Thus, for points P1 and P2,
c
c      P1 .lt. P2 implies that
c      * P1 is in a lower ordered cell than P2, or
c      * P1 is in the same cell as P2, but P1.X .lt. P2.X, or
c      * P1 is in the same cell as P2, P1.X = P2.X, but P1.Y .lt. P2.Y.
c
c    The arrays BIN_START(I,J) and BIN_LAST(I,J) are given the coordinates
c    I, J of a cell, and return the lowest and highest index in PSET of a
c    point in the I, J cell.  All indices in between also belong to
c    this cell.  If the cell has no points, then both arrays are -1.
c
c
c    After all this preprocessing, the algorithm for finding the nearest
c    point goes as follows:
c
c    1) for a point PTEST, determine the cell it belongs to.
c       Note that this algorithm will NOT be valid if PTEST lies outside
c       the bin limits.
c
c    2) Search for a cell that has at least one member of PSET in it.
c       We start at the cell containing PTEST, but if there are no members
c       there, we spiral out from the cell, one layer at a time.
c
c    3) Within this cell, find the point nearest to PTEST.  We now know that
c       we do not need to search any cell whose points will all be further
c       from PTEST than this distance.
c
c    4) Now, search in all other cells that could have a point closer to
c       PTEST than what we have found so far.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    22 July 2008
c
c  Author:
c
c    John Burkardt, Max Gunzburger, Janet Peterson
c
c  Reference:
c
c    Jon Bentley, Bruce Weide, Andrew Yao,
c    Optimal Expected Time Algorithms for Closest Point Problems,
c    ACM Transactions on Mathematical Software,
c    Volume 6, Number 4, December 1980, pages 563-580.
c
c  Parameters:
c
c    Input, integer NSET, the number of points in the set.
c
c    Input, double precision PSET(2,NSET), the coordinates of the points
c    in the set.
c
c    Input, integer NBIN(2), the number of cells in the horizontal and
c    vertical directions.
c
c    Input, double precision BIN_MIN(2), BIN_MAX(2), the minimum and maximum
c    bin values.
c
c    Input, integer BIN_START(NBIN(1),NBIN(2)), BIN_LAST(NBIN(1),NBIN(2)),
c    indicates the index of the first and last element in the bin, or -1
c    if none.
c
c    Input, integer BIN_NEXT(NSET), the index of the next element of the bin
c    containing this element.
c
c    Input, double precision PTEST(2), the coordinates of the test point.
c
c    Output, integer I_MIN, the index of the nearest point in PSET to PTEST.
c
      implicit none

      integer ndim
      parameter ( ndim = 2 )

      integer nbin(ndim)
      integer nset

      integer bin(ndim)
      integer bin_last(nbin(1),nbin(2))
      double precision bin_max(ndim)
      double precision bin_min(ndim)
      integer bin_start(nbin(1),nbin(2))
      integer bin_next(nset)
      double precision c_max(ndim)
      double precision c_min(ndim)
      double precision d_min
      double precision d_min_sq
      double precision d_sq
      integer i
      integer i_min
      integer ic
      integer j
      integer jc
      integer layer
      double precision layer_width
      logical more_bins
      integer node
      double precision pset(ndim,nset)
      double precision ptest(ndim)
      double precision r8_huge
      double precision search_radius
c
c  Special cases.
c
      if ( nset .le. 0 ) then
        d_min = r8_huge ( )
        i_min = 0
        return
      end if

      if ( nset .eq. 1 ) then
        d_min = 0.0
        do i = 1, ndim
          d_min = d_min + ( ptest(i) - pset(i,1) )**2
        end do
        d_min = sqrt ( d_min )
        i_min = 1
        return
      end if
c
c  The efficiency of the code will suffer if the data in the vector
c
c    ( bin_max(1:ndim) - bin_min(1:ndim) ) / dble ( nbin(1:ndim) )
c
c  varies significantly.
c
      layer_width = r8_huge ( )
      do i = 1, ndim
        layer_width = min ( layer_width, 
     &    ( bin_max(i) - bin_min(i) ) / dble ( nbin(i) ) )
      end do

      d_min_sq = r8_huge ( )
      i_min = 0
      search_radius = 0.0D+00
c
c  Determine the bin coordinates of the point P.
c
      call r8vec_to_bin_even3 ( ndim, nbin, bin_min, bin_max, 
     &  ptest, bin )
c
c  Determine the radius of the ball of space that will be completely
c  searched after this first bin, "layer 0", is done.
c
      call bin_to_r8vec_even3 ( ndim, nbin, bin, bin_min, bin_max, 
     &  c_min, c_max )
c
c  Set the central bin of the layers.
c
      ic = bin(1)
      jc = bin(2)

      layer = 0
c
c  Search all legal bins in layer LAYER.
c
10    continue

        more_bins = .false.

        call index_box2_next_2d ( layer, layer, ic, jc, i, j, 
     &    more_bins )
c
c  In layer LAYER, search each BIN I, J.
c
20      continue

          if ( 1 .le. i .and. i .le. nbin(1) .and. 
     &         1 .le. j .and. j .le. nbin(2) ) then

            node = bin_start(i,j)

30          continue

            if ( 0 .lt. node ) then

              d_sq = 0.0
              do i = 1, ndim
                d_sq = d_sq + ( ptest(i) - pset(i,node) )**2
              end do

              if ( d_sq .lt. d_min_sq ) then
                d_min_sq = d_sq
                i_min = node
              end if

              node = bin_next(node)

              go to 30

            end if

          end if
c
c  We have now searched all points in bin I, J.
c
c  Determine the next bin in this layer.
c
c  But if it lies outside the region, discard it, and go on to the next one.
c  Once you get past the last bin, exit because you are done the layer.
c
40        continue

            call index_box2_next_2d ( layer, layer, ic, jc, i, j, 
     &        more_bins )

            if ( .not. more_bins ) then
              go to 50
            end if

            if ( 1 .le. i .and. i .le. nbin(1) .and. 
     &           1 .le. j .and. j .le. nbin(2) ) then
              go to 50
            end if

          go to 40

50        continue

          if ( .not. more_bins ) then
            go to 60
          end if

        go to 20
c
c  We have completed layer LAYER.
c  Update the radius of the searched area.
c
60      continue

        if ( layer .eq. 0 ) then
          search_radius = r8_huge ( )
          do i = 1, ndim
            search_radius = min ( search_radius, 
     &        abs ( ptest(i) - c_min(i) ) )
          end do
          do i = 1, ndim
            search_radius = min ( search_radius, 
     &        abs ( ptest(i) - c_max(i) ) )
          end do
        else
          search_radius = search_radius + layer_width
        end if
c
c  We are done if:
c
c    * We have found at least one neighbor;
c    AND
c    * We have searched all legal boxes that contain the circle
c      with PTEST at the center and the nearest neighbor on the circumference.
c
        if ( i_min .ne. 0 ) then
          d_min = sqrt ( d_min_sq )
          if ( d_min .le. search_radius ) then
            go to 70
          end if
        end if

        layer = layer + 1

      go to 10

70    continue

      return
      end
      subroutine points_nearest_point_bins3_3d ( nset, pset, nbin, 
     &  bin_min, bin_max, bin_start, bin_last, bin_next, ptest, 
     &  i_min )

c*********************************************************************72
c
cc POINTS_NEAREST_POINT_BINS3_3D finds the nearest point to a given point in 3D.
c
c  Discussion:
c
c    This code differs from POINTS_NEAREST_POINTS_BINS_3D by allowing the
c    user to specify the number of bins in each dimension.
c
c    A set of NSET points with coordinates PSET is assumed to lie within a
c    box.  The limits of this rectangle are given in BIN_MIN and BIN_MAX.
c    The rectangle is divided up into an NBIN(1) by NBIN(2) by NBIN(3)
c    regular grid of cells.
c
c    The cells are ordered lexicographically, as suggested by the following
c    diagram for NBIN = (/ 5, 4, 2 /)
c
c             Z LAYER 1                       Z LAYER 2
c     *----*----*----*----*----*     *----*----*----*----*----*
c     |    |    |    |    |    |     |    |    |    |    |    |
c     | 16 | 17 | 18 | 19 | 20 |     | 36 | 37 | 38 | 39 | 40 |
c     |    |    |    |    |    |     |    |    |    |    |    |
c     *----*----*----*----*----*     *----*----*----*----*----*
c     |    |    |    |    |    |     |    |    |    |    |    |
c     | 11 | 12 | 13 | 14 | 15 |     | 31 | 32 | 33 | 34 | 35 |
c     |    |    |    |    |    |     |    |    |    |    |    |
c     *----*----*----*----*----*     *----*----*----*----*----*
c     |    |    |    |    |    |     |    |    |    |    |    |
c     |  6 |  7 |  8 |  9 | 10 |     | 26 | 27 | 28 | 29 | 30 |
c     |    |    |    |    |    |     |    |    |    |    |    |
c     *----*----*----*----*----*     *----*----*----*----*----*
c     |    |    |    |    |    |     |    |    |    |    |    |
c     |  1 |  2 |  3 |  4 |  5 |     | 21 | 22 | 23 | 24 | 25 |
c     |    |    |    |    |    |     |    |    |    |    |    |
c     *----*----*----*----*----*     *----*----*----*----*----*
c
c    The points in PSET are ordered by cell, and within a cell they
c    are ordered lexicographically.  Thus, for points P1 and P2,
c
c      P1 .lt. P2 implies that
c      * P1 is in a lower ordered cell than P2, or
c      * P1 is in the same cell as P2, but P1.X .lt. P2.X, or
c      * P1 is in the same cell as P2, P1.X = P2.X, but P1.Y .lt. P2.Y.
c      * P1 is in the same cell as P2, P1.X = P2.X, P1.Y = P2.Y,
c        but P1.Z .lt. P2.Z
c
c    The arrays BIN_START(I,J,K) and BIN_LAST(I,J,K) are given the coordinates
c    I, J, K of a cell, and return the lowest and highest index in PSET of a
c    point in the I, J, K cell.  All indices in between also belong to
c    this cell.  If the cell has no points, then both arrays are -1.
c
c
c    After all this preprocessing, the algorithm for finding the nearest
c    point goes as follows:
c
c    1) for a point PTEST, determine the cell it belongs to.
c       Note that this algorithm will NOT be valid if PTEST lies outside
c       the bin limits.
c
c    2) Search for a cell that has at least one member of PSET in it.
c       We start at the cell containing PTEST, but if there are no members
c       there, we spiral out from the cell, one layer at a time.
c
c    3) Within this cell, find the point nearest to PTEST.  We now know that
c       we do not need to search any cell whose points will all be further
c       from PTEST than this distance.
c
c    4) Now, search in all other cells that could have a point closer to
c       PTEST than what we have found so far.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    22 July 2008
c
c  Author:
c
c    John Burkardt, Max Gunzburger, Janet Peterson
c
c  Reference:
c
c    Jon Bentley, Bruce Weide, Andrew Yao,
c    Optimal Expected Time Algorithms for Closest Point Problems,
c    ACM Transactions on Mathematical Software,
c    Volume 6, Number 4, December 1980, pages 563-580.
c
c  Parameters:
c
c    Input, integer NSET, the number of points in the set.
c
c    Input, double precision PSET(3,NSET), the coordinates of the points
c    in the set.
c
c    Input, integer NBIN(3), the number of cells in the X, Y and Z directions.
c
c    Input, double precision BIN_MIN(3), BIN_MAX(3), the minimum and
c    maximum bin values.
c
c    Input, integer BIN_START(NBIN(1),NBIN(2),NBIN(3)),
c    BIN_LAST(NBIN(1),NBIN(2),NBIN(3)), the index of the first and last
c    element in the bin, or -1 if none.
c
c    Input, integer BIN_NEXT(NSET), the index of the next element of the bin
c    containing this element.
c
c    Input, double precision PTEST(3), the coordinates of the test points.
c
c    Output, integer I_MIN, the index of the nearest point in PSET to PTEST.
c
      implicit none

      integer ndim 
      parameter ( ndim = 3 )

      integer nbin(ndim)
      integer nset

      integer bin(ndim)
      integer bin_last(nbin(1),nbin(2),nbin(3))
      double precision bin_max(ndim)
      double precision bin_min(ndim)
      integer bin_start(nbin(1),nbin(2),nbin(3))
      integer bin_next(nset)
      double precision c_max(ndim)
      double precision c_min(ndim)
      double precision d_min
      double precision d_min_sq
      double precision d_sq
      integer i
      integer i_min
      integer ic
      integer j
      integer jc
      integer k
      integer kc
      integer layer
      double precision layer_width
      logical more_bins
      integer node
      double precision pset(ndim,nset)
      double precision ptest(ndim)
      double precision r8_huge
      double precision search_radius
c
c  Special cases.
c
      if ( nset .le. 0 ) then
        d_min = r8_huge ( )
        i_min = 0
        return
      end if

      if ( nset .eq. 1 ) then
        d_min = 0.0
        do i = 1, ndim
          d_min = d_min + ( ptest(i) - pset(i,1) )**2
        end do
        d_min = sqrt ( d_min )
        i_min = 1
        return
      end if
c
c  The efficiency of the code will suffer if the data in the vector
c
c    bin_max(1:ndim) - bin_min(1:ndim) / real ( nbin(1:ndim) )
c
c  varies significantly.
c
      layer_width = r8_huge ( )
      do i = 1, ndim
        layer_width = min ( layer_width, 
     &    ( bin_max(i) - bin_min(i) ) / dble ( nbin(i) ) )
      end do

      d_min_sq = r8_huge ( )
      i_min = 0
      search_radius = 0.0D+00
c
c  Determine the bin coordinates of the point P.
c
      call r8vec_to_bin_even3 ( ndim, nbin, bin_min, bin_max, 
     &  ptest(1), bin )
c
c  Determine the radius of the ball of space that will be completely
c  searched after this first bin, "layer 0", is done.
c
      call bin_to_r8vec_even3 ( ndim, nbin, bin, bin_min, bin_max, 
     &  c_min, c_max )
c
c  Set the central bin of the layers.
c
      ic = bin(1)
      jc = bin(2)
      kc = bin(3)

      layer = 0
c
c  Search all legal bins in layer LAYER.
c
10    continue

        more_bins = .false.

        call index_box2_next_3d ( layer, layer, layer, ic, jc, kc, 
     &    i, j, k, more_bins )
c
c  In layer LAYER, search each BIN I, J, K.
c
20      continue

          if ( 1 .le. i .and. i .le. nbin(1) .and. 
     &         1 .le. j .and. j .le. nbin(2) .and. 
     &         1 .le. k .and. k .le. nbin(3) ) then

            node = bin_start(i,j,k)

30          continue

            if ( 0 .lt. node ) then

              d_sq = 0.0
              do i = 1, ndim
                d_sq = d_sq + ( ptest(i) - pset(i,node) )**2
              end do

              if ( d_sq .lt. d_min_sq ) then
                d_min_sq = d_sq
                i_min = node
              end if

              node = bin_next(node)

              go to 30

            end if

          end if
c
c  We have now searched all points in bin I, J, K.
c
c  Determine the next bin in this layer.
c
c  But if it lies outside the region, discard it, and go on to the next one.
c  Once you get past the last bin, exit because you are done the layer.
c
40        continue

            call index_box2_next_3d ( layer, layer, layer, ic, jc, 
     &        kc, i, j, k, more_bins )

            if ( .not. more_bins ) then
              go to 50
            end if

            if ( 1 .le. i .and. i .le. nbin(1) .and. 
     &           1 .le. j .and. j .le. nbin(2) .and. 
     &           1 .le. k .and. k .le. nbin(3) ) then
              go to 50
            end if

          go to 40

50        continue

          if ( .not. more_bins ) then
            go to 60
          end if

        go to 20
c
c  We have completed layer LAYER.
c  Update the radius of the searched area.
c
60      continue

        if ( layer .eq. 0 ) then
          search_radius = r8_huge ( )
          do i = 1, ndim
            search_radius = min ( search_radius, 
     &        abs ( ptest(i) - c_min(i) ) )
          end do
          do i = 1, ndim
            search_radius = min ( search_radius, 
     &        abs ( ptest(i) - c_max(i) ) )
          end do
        else
          search_radius = search_radius + layer_width
        end if
c
c  We are done with PTEST if:
c
c    * We have found at least one neighbor;
c    AND
c    * We have searched all legal boxes that contain the circle
c      with PTEST at the center and the nearest neighbor on the circumference.
c
        if ( i_min .ne. 0 ) then
          d_min = sqrt ( d_min_sq )
          if ( d_min .le. search_radius ) then
            go to 70
          end if
        end if

        layer = layer + 1

      go to 10

70    continue
c
c  We are now done with all the layers.
c
      return
      end
      subroutine quality ( ndim, n, cell_moment, cell_volume,
     &  region_volume )

c*********************************************************************72
c
cc QUALITY computes some quality measures for a set of points in a region.
c
c  Discussion:
c
c    The quality measures report on how evenly spread the points are.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    22 July 2008
c
c  Author:
c
c    John Burkardt, Max Gunzburger, Janet Peterson
c
c  Reference:
c
c    John Burkardt, Max Gunzburger, Janet Peterson and Rebecca Brannon,
c    User Manual and Supporting Information for Library of Codes
c    for Centroidal Voronoi Placement and Associated Zeroth,
c    First, and Second Moment Determination,
c    Sandia National Laboratories Technical Report SAND2002-0099,
c    February 2002.
c
c  Parameters:
c
c    Input, integer NDIM, the spatial dimension.
c
c    Input, integer N, the number of cell generators.
c
c    Input, double precision CELL_MOMENT(NDIM,NDIM,N), the second moment
c    matrix for each Voronoi cell.
c
c    Input, double precision CELL_VOLUME(N), the Voronoi cell volumes.
c
c    Input, double precision REGION_VOLUME, the volume of the region,
c    as input by the user or estimated by VCM.
c
      implicit none

      integer n
      integer ndim

      double precision cell_det(n)
      double precision cell_moment(ndim,ndim,n)
      double precision cell_trace(n)
      double precision cell_volume(n)
      double precision dd_l1
      double precision dd_l2
      double precision dd_linf
      double precision r8mat_det_2d
      double precision r8mat_det_3d
      double precision ev
      double precision ev_l1
      double precision ev_l2
      double precision ev_linf
      integer i
      integer j
      integer k
      double precision matrix2(2,2)
      double precision matrix3(3,3)
      double precision r8_huge
      double precision region_volume
      double precision tr
      double precision tr_l1
      double precision tr_l2
      double precision tr_linf
c
c  Measure 1: the deviation of the cell volumes from the expected cell volume.
c
      ev = region_volume / dble ( n )

      ev_linf = - r8_huge ( )
      do i = 1, n
        ev_linf = max ( ev_linf, abs ( ev - cell_volume(i) ) )
      end do

      ev_l1 = 0.0
      do i = 1, n
        ev_l1 = ev_l1 + abs ( ev - cell_volume(i) )
      end do

      ev_l2 = 0.0
      do i = 1, n
        ev_l2 = ev_l2 + ( ev - cell_volume(i) )**2
      end do
      ev_l2 = sqrt ( ev_l2 )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'QUALITY'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Measure #1:'
      write ( *, '(a)' ) '    ( Cell_Volume - Expected Cell Volume )'
      write ( *, '(a)' ) ' '
      write ( *, '(a,g14.6)' ) '  Expected Cell Volume = ', ev
      write ( *, '(a)' ) ' '
      write ( *, '(a,g14.6)' ) '  L1 norm =              ', ev_l1
      write ( *, '(a,g14.6)' ) '  L2 norm =              ', ev_l2
      write ( *, '(a,g14.6)' ) '  L1 norm / N =          ', ev_l1 
     &  / dble ( n )
      write ( *, '(a,g14.6)' ) '  L2 norm / sqrt ( N ) = ', ev_l2 
     &  / sqrt ( dble ( n ) )
      write ( *, '(a,g14.6)' ) '  L-Inf norm =           ', ev_linf
c
c  Measure 2: the deviation of the traces of the cell second moment matrices
c  from the average.
c
      do i = 1, n
        cell_trace(i) = 0.0D+00
        do j = 1, ndim
          cell_trace(i) = cell_trace(i) + cell_moment(j,j,i)
        end do
      end do

      tr = 0.0
      do i = 1, n
        tr = tr + cell_trace(i)
      end do
      tr = tr / dble ( n )

      tr_linf = - r8_huge ( )
      do i = 1, n
        tr_linf = max ( tr_linf, abs ( tr - cell_trace(i) ) )
      end do

      tr_l1 = 0.0
      do i = 1, n
        tr_l1 = tr_l1 + abs ( tr - cell_trace(i) )
      end do

      tr_l2 = 0.0
      do i = 1, n
        tr_l2 = tr_l2 + ( tr - cell_trace(i) )**2
      end do
      tr_l2 = sqrt ( tr_l2 )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Measure #2:'
      write ( *, '(a)' ) '    ( Cell_Trace - Average Cell Trace )'
      write ( *, '(a)' ) ' '
      write ( *, '(a,g14.6)' ) '  Average Cell Trace = ', tr
      write ( *, '(a)' ) ' '
      write ( *, '(a,g14.6)' ) '  L1 norm =              ', tr_l1
      write ( *, '(a,g14.6)' ) '  L2 norm =              ', tr_l2
      write ( *, '(a,g14.6)' ) '  L1 norm / N =          ', tr_l1 
     &  / dble ( n )
      write ( *, '(a,g14.6)' ) '  L2 norm / sqrt ( N ) = ', tr_l2 
     &  / sqrt ( dble ( n ) )
      write ( *, '(a,g14.6)' ) '  L-Inf norm =           ', tr_linf
c
c  Measure 3: the determinant of the deviatoric matrix
c
      if ( ndim .eq. 2 ) then

        do i = 1, n

          do k = 1, ndim
            do j = 1, ndim
              matrix2(j,k) = cell_moment(j,k,i)
            end do
          end do

          do j = 1, ndim
            matrix2(j,j) = matrix2(j,j) - cell_trace(i) / dble ( ndim )
          end do

          cell_det(i) = r8mat_det_2d ( matrix2 )

        end do

      else if ( ndim .eq. 3 ) then

        do i = 1, n

          do k = 1, ndim
            do j = 1, ndim
              matrix3(j,k) = cell_moment(j,k,i)
            end do
          end do

          do j = 1, ndim
            matrix3(j,j) = matrix3(j,j) - cell_trace(i) / dble ( ndim )
          end do

          cell_det(i) = r8mat_det_3d ( matrix3 )

        end do

      end if

      dd_linf = - r8_huge ( )
      do i = 1, n
        dd_linf = max ( dd_linf, abs ( cell_det(i) ) )
      end do

      dd_l1 = 0.0
      do i = 1, n
        dd_l1 = dd_l1 + abs ( cell_det(i) )
      end do

      dd_l2 = 0.0
      do i = 1, n
        dd_l2 = dd_l2 + cell_det(i)**2
      end do
      dd_l2 = sqrt ( dd_l2 )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Measure #3:'
      write ( *, '(a)' ) 
     &  '    ( The determinant of the deviatoric matrix )'
      write ( *, '(a)' ) ' '
      write ( *, '(a,g14.6)' ) '  L1 norm = 	     ', dd_l1
      write ( *, '(a,g14.6)' ) '  L2 norm =              ', dd_l2
      write ( *, '(a,g14.6)' ) '  L1 norm / N =          ', dd_l1 
     &  / real ( n, kind = 8 )
      write ( *, '(a,g14.6)' ) '  L2 norm / sqrt ( N ) = ', dd_l2 
     &  / sqrt ( real ( n, kind = 8 ) )
      write ( *, '(a,g14.6)' ) '  L-Inf norm =           ', dd_linf

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
c    John Burkardt, Max Gunzburger, Janet Peterson
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
      subroutine r8_to_bin_even2 ( nbin, a, b, c, bin )

c*********************************************************************72
c
cc R8_TO_BIN_EVEN2 determines the appropriate "bin" for C in [A,B].
c
c  Discussion:
c
c    The interval from A to B is divided into NBIN equal subintervals or bins.
c
c  Example:
c
c    NBIN = 5, A = 5, B = 15
c
c    <-1-+-2-+-3-+-4-+-5->
c    5   7   9  11  13  15
c
c
c    C   BIN
c
c    1    1
c    3    1
c    4.9  1
c    5    1
c    6    1
c    7.1  2
c    8    2
c    9.5  3
c   12    4
c   14    5
c   15    5
c   15.1  5
c   99    5
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    22 July 2008
c
c  Author:
c
c    John Burkardt, Max Gunzburger, Janet Peterson
c
c  Parameters:
c
c    Input, integer NBIN, the number of bins.
c
c    Input, double precision A, B, the lower and upper limits of the bin
c    interval.  While A is expected to be less than B, the code should
c    return useful results if A is actually greater than B.
c
c    Input, double precision C, a value to be placed in a bin.
c
c    Output, integer BIN, the index of the bin to which C is assigned.
c
      implicit none

      double precision a
      double precision a2
      double precision b
      double precision b2
      integer bin
      double precision c
      integer nbin
      logical switch
c
c  Take care of special cases.
c
      if ( nbin .lt. 1 ) then
        bin = 0
        return
      end if

      if ( nbin .eq. 1 ) then
        bin = 1
        return
      end if

      if ( b .eq. a ) then
        bin = 1
        return
      end if
c
c  If the limits are descending, then we switch them now, and
c  unswitch the results at the end.
c
      if ( a .lt. b ) then
        switch = .false.
        a2 = a
        b2 = b
      else
        switch = .true.
        a2 = b
        b2 = a
      end if
c
c  Compute the bin.
c
      if ( c .le. a2 ) then
        bin = 1
      else if ( b2 .le. c ) then
        bin = nbin
      else
        bin = 1 + int ( dble ( nbin ) * ( c - a2 ) / ( b2 - a2 ) )
        bin = max ( bin, 1 )
        bin = min ( bin, nbin )
      end if
c
c  Reverse the switching.
c
      if ( switch ) then
        bin = nbin + 1 - bin
      end if

      return
      end
      subroutine r82vec_bin_even3 ( n, a, nbin, bin_min, bin_max, 
     &  bin_start, bin_last, bin_next )

c*********************************************************************72
c
cc R82VEC_BIN_EVEN3 bins an R82VEC into evenly spaced bins.
c
c  Discussion:
c
c    A different number of bins may be used in each dimension.
c
c    This is only a partial, indexed, sorting of the data.  To sort
c    the data, it is necessary to build a new array by extracting the
c    data for each bin, sorting that, and appending it to the array.
c
c    There are NBIN(1) 1D bins in the X direction, NBIN(2) for Y, making a
c    total of NBIN(1) * NBIN(2) 2D bins.  Each set of 1D bins begins and
c    ends at user specified mininum and maximum values.
c
c    The 2D bins are indexed by the X and Y bins that construct them,
c    and ordered lexicographically by these indices:
c
c      1,4 | 2,4 | 3,4 | 4,4 | 5,4
c      ----+-----+-----+-----+-----
c      1,3 | 2,3 | 3,3 | 4,3 | 5,3
c      ----+-----+-----+-----+-----
c      1,2 | 2,2 | 3,2 | 4,2 | 5,2
c      ----+-----+-----+-----+-----
c      1,1 | 2,1 | 3,1 | 4,1 | 5,1
c
c    Thus, the 2D bin sequence is (1,1), (1,2), (1,3), (1,4), (2,1),
c    ..., (5,4).
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    22 July 2008
c
c  Author:
c
c    John Burkardt, Max Gunzburger, Janet Peterson
c
c  Reference:
c
c    Thomas Cormen, Charles Leiserson, Ronald Rivest,
c    Introduction to Algorithms,
c    MIT Press, page 180.
c
c  Parameters:
c
c    Input, integer N, the number of points in the data set.
c
c    Input, double precision A(2,N), the D2 data to be binned.
c
c    Input, integer NBIN(2), the number of bins in each dimension.
c    NBIN must be at least 1.
c
c    Input, double precision BIN_MIN(2), BIN_MAX(2), the bin limits.
c
c    Output, integer BIN_START(NBIN(1),NBIN(2)), BIN_LAST(NBIN(1),NBIN(2)),
c    the index of the first and last elements of A that went into each bin,
c    or -1 if there are no entries in the bin.
c
c    Output, integer BIN_NEXT(N), contains the index of the next element of A
c    that follows this element in the same bin.  A value of 0 means this
c    is the last entry in the particular bin.
c
      implicit none

      integer ndim
      parameter ( ndim = 2 )

      integer n
      integer nbin(ndim)

      double precision a(ndim,n)
      integer bin(ndim)
      integer bin_last(nbin(1),nbin(2))
      integer bin_next(n)
      integer bin_start(nbin(1),nbin(2))
      double precision bin_max(ndim)
      double precision bin_min(ndim)
      integer i
      integer i1
      integer i2
      integer j
      integer k
c
c  Initialize the bin pointers to -1.
c
      do j = 1, nbin(2)
        do i = 1, nbin(1)
          bin_last(i,j) = -1
          bin_start(i,j) = -1
        end do
      end do

      do i = 1, n
        bin_next(i) = -1
      end do
c
c  Build up linked lists of entries that go into each bin.
c
      do j = 1, n

        call r8vec_to_bin_even3 ( ndim, nbin, bin_min, bin_max, 
     &    a(1,j), bin )

        i1 = bin(1)
        i2 = bin(2)

        if ( bin_start(i1,i2) .eq. -1 ) then
          bin_start(i1,i2) = j
        else
          k = bin_last(i1,i2)
          bin_next(k) = j
        end if

        bin_next(j) = 0

        bin_last(i1,i2) = j

      end do

      return
      end
      subroutine r82vec_binned_reorder2 ( n, a, nbin, bin_start,
     &  bin_last, bin_next )

c*********************************************************************72
c
cc R82VEC_BINNED_REORDER2 reorders a binned R82VEC.
c
c  Discussion:
c
c    This routine allows there to be a different number of bins in
c    each dimension.
c
c    The bin vectors define an implicit ordering of the data
c    vector.  This routine physically rearranges the data vector
c    so that items in the first bin come first, and so on.  The
c    data within a bin is not reordered.
c
c    On output, the BIN_START and BIN_NEXT arrays have also been updated
c    so that they still correspond to the (rearranged) vector A.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    22 July 2008
c
c  Author:
c
c    John Burkardt, Max Gunzburger, Janet Peterson
c
c  Parameters:
c
c    Input, integer N, the number of data points.
c
c    Input/output, double precision A(2,N), the D2 data to be sorted.
c
c    Input, integer NBIN(2), the number of bins in each direction.
c
c    Input/output, integer BIN_START(NBIN(1),NBIN(2)),
c    BIN_LAST(NBIN(1),NBIN(2)), contains the index of the first and last 
c    element of A that went into each bin, or -1 if there are no entries 
c    in the bin.
c
c    Input/output, integer BIN_NEXT(N), contains the index of the next
c    element of A that follows this element in the same bin.  A value of
c    0 means this is the last entry in the particular bin.
c
      implicit none

      integer ndim
      parameter ( ndim = 2 )

      integer n
      integer nbin(ndim)

      double precision a(ndim,n)
      double precision a2(ndim,n)
      integer bin_last(nbin(1),nbin(2))
      integer bin_next(n)
      integer bin_start(nbin(1),nbin(2))
      integer i
      integer i1
      integer i2
      integer j
      integer k
c
c  Bin by bin, copy the contents of A to A2.
c  The BIN_START array is also updated as we go.
c
      k = 0

      do i1 = 1, nbin(1)

        do i2 = 1, nbin(2)

          j = bin_start(i1,i2)

          if ( 0 .lt. j ) then
            bin_start(i1,i2) = k + 1
          end if

10        continue

          if ( 0 .lt. j ) then
            k = k + 1
            bin_last(i1,i2) = k
            do i = 1, ndim
              a2(i,k) = a(i,j)
            end do
            j = bin_next(j)
            go to 10
          end if

        end do

      end do
c
c  Copy A2 back into A.
c
      do j = 1, n
        do i = 1, ndim
          a(i,j) = a2(i,j)
        end do
      end do
c
c  Now update the BIN_NEXT array.
c
      do i = 1, n
        bin_next(i) = i + 1
      end do

      do i1 = 1, nbin(1)
        do i2 = 1, nbin(2)

          k = bin_last(i1,i2)

          if ( 0 .lt. k ) then
            bin_next(k) = 0
          end if

        end do
      end do

      return
      end
      subroutine r82vec_binned_sort_a2 ( n, a, nbin, bin_start, 
     &  bin_last )

c*********************************************************************72
c
cc R82VEC_BINNED_SORT_A2 sorts each bin of a binned R82VEC.
c
c  Discussion:
c
c    This routine allows a different number of bins in each dimension.
c
c    Presumably, the data vector was first binned by R82VEC_BIN_EVEN3,
c    then reordered by R82VEC_BINNED_REORDER2.  Now, each of the
c    bins of data is sorted one at a time.
c
c    The result is NOT a lexicographically sorted D2 vector.
c
c    What is true is that if I .lt. J, then either the I-th element of A occurs
c    in a lexicographically smaller bin than J, or they share a bin,
c    and the I-th element is lexicographically less than or equal to
c    the J-th element.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    22 July 2008
c
c  Author:
c
c    John Burkardt, Max Gunzburger, Janet Peterson
c
c  Parameters:
c
c    Input, integer N, the number of data points.
c
c    Input/output, double precision A(2,N), the R2 data to be sorted.
c
c    Input, integer NBIN(2), the number of bins in each dimension.
c
c    Input, integer BIN_START(NBIN(1),NBIN(2)), BIN_LAST(NBIN(1),NBIN(2)),
c    the index of the first and last element of A that went into each bin, or -1
c    if there are no entries in this bin.
c
      implicit none

      integer ndim
      parameter ( ndim = 2 )

      integer n
      integer nbin(ndim)

      double precision a(ndim,n)
      integer bin_last(nbin(1),nbin(2))
      integer bin_start(nbin(1),nbin(2))
      integer i1
      integer i2
      integer j1
      integer j2
      integer n1

      do i1 = 1, nbin(1)

        do i2 = 1, nbin(2)

          j1 = bin_start(i1,i2)

          if ( 0 .lt. j1 ) then

            j2 = bin_last(i1,i2)

            n1 = j2 + 1 - j1

            if ( 1 .lt. n1 ) then
              call r8col_sort_quick_a ( ndim, n1, a(1,j1) )
            end if

          end if

        end do

      end do

      return
      end
      subroutine r83vec_bin_even3 ( n, a, nbin, bin_min, bin_max, 
     &  bin_start, bin_last, bin_next )

c*********************************************************************72
c
cc R83VEC_BIN_EVEN3 bins an R83VEC into evenly spaced bins.
c
c  Discussion:
c
c    A different number of bins may be used in each dimension.
c
c    This is only a partial, indexed, sorting of the data.  To sort
c    the data, it is necessary to build a new array by extracting the
c    data for each bin, sorting that, and appending it to the array.
c
c    There are NBIN(1) 1D bins in the X direction, NBIN(2) for Y, and
c    NBIN(3) for Z, making a total of NBIN(1) * NBIN(2) * NBIN(3) 3D bins.
c    Each set of 1D bins begins and ends at user specified mininum and
c    maximum values.
c
c    The 3D bins are indexed by the X, Y and Z bins that construct them,
c    and ordered lexicographically by these indices.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    22 July 2008
c
c  Author:
c
c    John Burkardt, Max Gunzburger, Janet Peterson
c
c  Reference:
c
c    Thomas Cormen, Charles Leiserson, Ronald Rivest,
c    Introduction to Algorithms,
c    MIT Press, page 180.
c
c  Parameters:
c
c    Input, integer N, the number of points in the data set.
c
c    Input, double precision A(3,N), the R3 data to be binned.
c
c    Input, integer NBIN(3), the number of bins in each dimension.
c    NBIN must be at least 1.
c
c    Input, double precision BIN_MIN(3), BIN_MAX(3), the bin limits.
c
c    Output, integer BIN_START(NBIN(1),NBIN(2),NBIN(3)),
c    BIN_LAST(NBIN(1),NBIN(2),NBIN(3)), contains the
c    index of the first and last elements of A that went into each bin,
c    or -1 if there are no entries in the bin.
c
c    Output, integer BIN_NEXT(N), contains the index of the next element of A
c    that follows this element in the same bin.  A value of 0 means this
c    is the last entry in the particular bin.
c
      implicit none

      integer ndim
      parameter ( ndim = 3 )

      integer n
      integer nbin(ndim)

      double precision a(ndim,n)
      integer bin(ndim)
      integer bin_last(nbin(1),nbin(2),nbin(3))
      integer bin_next(n)
      integer bin_start(nbin(1),nbin(2),nbin(3))
      double precision bin_max(ndim)
      double precision bin_min(ndim)
      integer i
      integer i1
      integer i2
      integer i3
      integer j
      integer k
c
c  Initialize the bin pointers to -1.
c
      do k = 1, nbin(3)
        do j = 1, nbin(2)
          do i = 1, nbin(1)
            bin_last(i,j,k) = -1
            bin_start(i,j,k) = -1
          end do
        end do
      end do

      do i = 1, n
        bin_next(i) = -1
      end do
c
c  Build up linked lists of entries that go into each bin.
c
      do j = 1, n

        call r8vec_to_bin_even3 ( ndim, nbin, bin_min, bin_max, 
     &    a(1,j), bin )

        i1 = bin(1)
        i2 = bin(2)
        i3 = bin(3)

        if ( bin_start(i1,i2,i3) .eq. -1 ) then
          bin_start(i1,i2,i3) = j
        else
          k = bin_last(i1,i2,i3)
          bin_next(k) = j
        end if

        bin_next(j) = 0

        bin_last(i1,i2,i3) = j

      end do

      return
      end
      subroutine r83vec_binned_reorder2 ( n, a, nbin, bin_start, 
     &  bin_last, bin_next )

c*********************************************************************72
c
cc R83VEC_BINNED_REORDER2 reorders a binned R83VEC.
c
c  Discussion:
c
c    This routine allows there to be a different number of bins in
c    each dimension.
c
c    The bin vectors define an implicit ordering of the data
c    vector.  This routine physically rearranges the data vector
c    so that items in the first bin come first, and so on.  The
c    data within a bin is not reordered.
c
c    On output, the BIN_START and BIN_NEXT arrays have also been updated
c    so that they still correspond to the (rearranged) vector A.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    22 July 2008
c
c  Author:
c
c    John Burkardt, Max Gunzburger, Janet Peterson
c
c  Parameters:
c
c    Input, integer N, the number of data points.
c
c    Input/output, double precision A(3,N), the data to be sorted.
c
c    Input, integer NBIN(3), the number of bins in each dimension.
c
c    Input/output, integer BIN_START(NBIN(1),NBIN(2),NBIN(3)),
c    BIN_LAST(NBIN(1),NBIN(2),NBIN(3)), contains the
c    index of the first and last element of A that went into each bin,
c    or -1 if there are no entries in the bin.
c
c    Input/output, integer BIN_NEXT(N), contains the index of the next
c    element of A that follows this element in the same bin.  A value of 0
c    means this is the last entry in the particular bin.
c
      implicit none

      integer ndim
      parameter ( ndim = 3 )

      integer n
      integer nbin(ndim)

      double precision a(ndim,n)
      double precision a2(ndim,n)
      integer bin_last(nbin(1),nbin(2),nbin(3))
      integer bin_next(n)
      integer bin_start(nbin(1),nbin(2),nbin(3))
      integer i
      integer i1
      integer i2
      integer i3
      integer j
      integer k
c
c  Bin by bin, copy the contents of A to A2.
c  The BIN_START array is also updated as we go.
c
      k = 0

      do i1 = 1, nbin(1)

        do i2 = 1, nbin(2)

          do i3 = 1, nbin(3)

            j = bin_start(i1,i2,i3)

            if ( 0 .lt. j ) then
              bin_start(i1,i2,i3) = k + 1
            end if

10          continue

            if ( 0 .lt. j ) then
              k = k + 1
              bin_last(i1,i2,i3) = k
              do i = 1, ndim
                a2(i,k) = a(i,j)
              end do
              j = bin_next(j)
              go to 10
            end if

          end do

        end do

      end do
c
c  Copy A2 back into A.
c
      do j = 1, n
        do i = 1, ndim
          a(i,j) = a2(i,j)
        end do
      end do
c
c  Now update the BIN_NEXT array.
c
      do i = 1, n
        bin_next(i) = i + 1
      end do

      do i1 = 1, nbin(1)
        do i2 = 1, nbin(2)
          do i3 = 1, nbin(3)

            k = bin_last(i1,i2,i3)

            if ( 0 .lt. k ) then
              bin_next(k) = 0
            end if

          end do
        end do
      end do

      return
      end
      subroutine r83vec_binned_sort_a2 ( n, a, nbin, bin_start, 
     &  bin_last )

c*********************************************************************72
c
cc R83VEC_BINNED_SORT_A2 sorts each bin of a binned R83VEC.
c
c  Discussion:
c
c    This routine allows a different number of bins in each dimension.
c
c    Presumably, the data vector was first binned by R83VEC_BIN_EVEN3,
c    then reordered by R83VEC_BINNED_REORDER2.  Now, each of the
c    bins of data is sorted one at a time.
c
c    The result is NOT a lexicographically sorted D3 vector.
c
c    What is true is that if I .lt. J, then either the I-th element of A occurs
c    in a lexicographically smaller bin than J, or they share a bin,
c    and the I-th element is lexicographically less than or equal to
c    the J-th element.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    22 July 2008
c
c  Author:
c
c    John Burkardt, Max Gunzburger, Janet Peterson
c
c  Parameters:
c
c    Input, integer N, the number of data points.
c
c    Input/output, double precision A(3,N), the data to be sorted.
c
c    Input, integer NBIN(3), the number of bins in each dimension.
c
c    Input, integer BIN_START(NBIN(1),NBIN(2),NBIN(3)),
c    BIN_LAST(NBIN(1),NBIN(2),NBIN(3)), contains
c    the index of the first and last element of A that went into each bin,
c    or -1 if there are no entries in this bin.
c
      implicit none

      integer ndim 
      parameter ( ndim = 3 )

      integer n
      integer nbin(ndim)

      double precision a(ndim,n)
      integer bin_last(nbin(1),nbin(2),nbin(3))
      integer bin_start(nbin(1),nbin(2),nbin(3))
      integer i1
      integer i2
      integer i3
      integer j1
      integer j2
      integer n1

      do i1 = 1, nbin(1)

        do i2 = 1, nbin(2)

          do i3 = 1, nbin(3)

            j1 = bin_start(i1,i2,i3)

            if ( 0 .lt. j1 ) then

              j2 = bin_last(i1,i2,i3)

              n1 = j2 + 1 - j1

              if ( 1 .lt. n1 ) then
                call r8col_sort_quick_a ( ndim, n1, a(1,j1) )
              end if

            end if

          end do

        end do

      end do

      return
      end
      subroutine r8col_part_quick_a ( m, n, a, l, r )

c*********************************************************************72
c
cc R8COL_PART_QUICK_A reorders the columns of an array as part of a quick sort.
c
c  Discussion:
c
c    The routine reorders the columns of A.  Using A(1:M,1) as a
c    key, all entries of A that are less than or equal to the key will
c    precede the key, which precedes all entries that are greater than the key.
c
c  Example:
c
c    Input:
c
c      M = 2, N = 8
c      A = ( 2  8  6  0 10 10  0  5
c            4  8  2  2  6  0  6  8 )
c
c    Output:
c
c      L = 2, R = 4
c
c      A = (  0  0  2  8  6 10 10  4
c             2  6  4  8  2  6  0  8 )
c             ----     -------------
c             LEFT KEY     RIGHT
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    22 July 2008
c
c  Author:
c
c    John Burkardt, Max Gunzburger, Janet Peterson
c
c  Parameters:
c
c    Input, integer M, the row dimension of A, and the length of a column.
c
c    Input, integer N, the column dimension of A.
c
c    Input/output, double precision A(M,N).  On input, the array to be checked.
c    On output, A has been reordered as described above.
c
c    Output, integer L, R, the indices of A that define the three segments.
c    Let KEY = the input value of A(1:M,1).  Then
c    I .le. L                         A(1:M,I) .lt. KEY;
c           L .lt. I .lt. R           A(1:M,I) = KEY;
c                         R .le. I    A(1:M,I) .gt. KEY.
c
      implicit none

      integer m
      integer n

      double precision a(m,n)
      integer i
      integer j
      integer k
      double precision key(m)
      integer l
      integer r
      logical r8vec_eq
      logical r8vec_gt
      logical r8vec_lt

      if ( n .lt. 1 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'R8COL_PART_QUICK_A - Fatal error!'
        write ( *, '(a)' ) '  N .lt. 1.'
        return
      end if

      if ( n .eq. 1 ) then
        l = 0
        r = 2
        return
      end if

      do i = 1, m
        key(i) = a(i,1)
      end do

      k = 1
c
c  The elements of unknown size have indices between L+1 and R-1.
c
      l = 1
      r = n + 1

      do i = 2, n

        if ( r8vec_gt ( m, a(1,l+1), key ) ) then
          r = r - 1
          call r8vec_swap ( m, a(1,r), a(1,l+1) )
        else if ( r8vec_eq ( m, a(1,l+1), key ) ) then
          k = k + 1
          call r8vec_swap ( m, a(1,k), a(1,l+1) )
          l = l + 1
        else if ( r8vec_lt ( m, a(1,l+1), key ) ) then
          l = l + 1
        end if

      end do
c
c  Shift small elements to the left.
c
      do j = 1, l - k
        do i = 1, m
          a(i,j) = a(i,j+k)
        end do
      end do
c
c  Shift KEY elements to center.
c
      do j = l-k+1, l
        do i = 1, m
          a(i,j) = key(i)
        end do
      end do
c
c  Update L.
c
      l = l - k

      return
      end
      subroutine r8col_sort_quick_a ( m, n, a )

c*********************************************************************72
c
cc R8COL_SORT_QUICK_A ascending sorts the columns of a table using quick sort.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    22 July 2008
c
c  Author:
c
c    John Burkardt, Max Gunzburger, Janet Peterson
c
c  Parameters:
c
c    Input, integer M, the row order of A, and the length of a column.
c
c    Input, integer N, the number of columns of A.
c
c    Input/output, double precision A(M,N).
c    On input, the array to be sorted.
c    On output, the array has been sorted.
c
      implicit none

      integer MAXLEVEL
      parameter ( MAXLEVEL = 25 )

      integer m
      integer n

      double precision a(m,n)
      integer base
      integer l_segment
      integer level
      integer n_segment
      integer rsave(MAXLEVEL)
      integer r_segment

      if ( n .lt. 1 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'R8COL_SORT_QUICK_A - Fatal error!'
        write ( *, '(a)' ) '  N .lt. 1.'
        stop
      end if

      if ( n .eq. 1 ) then
        return
      end if

      level = 1
      rsave(level) = n + 1
      base = 1
      n_segment = n

10    continue
c
c  Partition the segment.
c
        call r8col_part_quick_a ( m, n_segment, a(1,base), l_segment, 
     &    r_segment )
c
c  If the left segment has more than one element, we need to partition it.
c
        if ( 1 .lt. l_segment ) then

          if ( MAXLEVEL .lt. level ) then
            write ( *, '(a)' ) ' '
            write ( *, '(a)' ) 'R8COL_SORT_QUICK_A - Fatal error!'
            write ( *, '(a,i8)' ) 
     &        '  Exceeding recursion maximum of ', MAXLEVEL
            stop
          end if

          level = level + 1
          n_segment = l_segment
          rsave(level) = r_segment + base - 1
c
c  The left segment and the middle segment are sorted.
c  Must the right segment be partitioned?
c
        else if ( r_segment .lt. n_segment ) then

          n_segment = n_segment + 1 - r_segment
          base = base + r_segment - 1
c
c  Otherwise, we back up a level if there is an earlier one.
c
        else

20        continue

            if ( level .le. 1 ) then
              go to 40
            end if

            base = rsave(level)
            n_segment = rsave(level-1) - rsave(level)
            level = level - 1

            if ( 0 .lt. n_segment ) then
              go to 30
            end if

          go to 20

30        continue

        end if

      go to 10

40    continue

      return
      end
      function r8mat_det_2d ( a )

c*********************************************************************72
c
cc R8MAT_DET_2D computes the determinant of a 2 by 2 matrix.
c
c  Discussion:
c
c    The determinant of a 2 by 2 matrix is
c
c      a11 * a22 - a12 * a21.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    22 July 2008
c
c  Author:
c
c    John Burkardt, Max Gunzburger, Janet Peterson
c
c  Parameters:
c
c    Input, double precision A(2,2), the matrix whose determinant is desired.
c
c    Output, double precision RMAT_DET_2D, the determinant of the matrix.
c
      implicit none

      double precision a(2,2)
      double precision r8mat_det_2d

      r8mat_det_2d = a(1,1) * a(2,2) - a(1,2) * a(2,1)

      return
      end
      function r8mat_det_3d ( a )

c*********************************************************************72
c
cc R8MAT_DET_3D computes the determinant of a 3 by 3 matrix.
c
c  Discussion:
c
c    The determinant of a 3 by 3 matrix is
c
c        a11 * a22 * a33 - a11 * a23 * a32
c      + a12 * a23 * a31 - a12 * a21 * a33
c      + a13 * a21 * a32 - a13 * a22 * a31
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    22 July 2008
c
c  Author:
c
c    John Burkardt, Max Gunzburger, Janet Peterson
c
c  Parameters:
c
c    Input, double precision A(3,3), the matrix whose determinant is desired.
c
c    Output, double precision RMAT_DET_3D, the determinant of the matrix.
c
      implicit none

      double precision a(3,3)
      double precision r8mat_det_3d

      r8mat_det_3d = 
     &       a(1,1) * ( a(2,2) * a(3,3) - a(2,3) * a(3,2) ) 
     &     + a(1,2) * ( a(2,3) * a(3,1) - a(2,1) * a(3,3) ) 
     &     + a(1,3) * ( a(2,1) * a(3,2) - a(2,2) * a(3,1) )

      return
      end
      function r8vec_eq ( n, a1, a2 )

c*********************************************************************72
c
cc R8VEC_EQ is true if every pair of entries in two vectors is equal.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    22 July 2008
c
c  Author:
c
c    John Burkardt, Max Gunzburger, Janet Peterson
c
c  Parameters:
c
c    Input, integer N, the number of entries in the vectors.
c
c    Input, double precision A1(N), A2(N), two vectors to compare.
c
c    Output, logical R8VEC_EQ.
c    R8VEC_EQ is .TRUE. if every pair of elements A1(I) and A2(I) are equal,
c    and .FALSE. otherwise.
c
      implicit none

      integer n

      double precision a1(n)
      double precision a2(n)
      integer i
      logical r8vec_eq

      r8vec_eq = .false.

      do i = 1, n
        if ( a1(i) .ne. a2(i) ) then
          return
        end if
      end do

      r8vec_eq = .true.

      return
      end
      function r8vec_gt ( n, a1, a2 )

c*********************************************************************72
c
cc R8VEC_GT == ( A1 > A2 ) for double precision vectors.
c
c  Discussion:
c
c    The comparison is lexicographic.
c
c    A1 > A2  <=>                              A1(1) > A2(1) or
c                 ( A1(1)     == A2(1)     and A1(2) > A2(2) ) or
c                 ...
c                 ( A1(1:N-1) == A2(1:N-1) and A1(N) > A2(N)
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    22 July 2008
c
c  Author:
c
c    John Burkardt, Max Gunzburger, Janet Peterson
c
c  Parameters:
c
c    Input, integer N, the dimension of the vectors.
c
c    Input, double precision A1(N), A2(N), the vectors to be compared.
c
c    Output, logical R8VEC_GT, is TRUE if and only if A1 > A2.
c
      implicit none

      integer n

      double precision a1(n)
      double precision a2(n)
      integer i
      logical r8vec_gt

      r8vec_gt = .false.

      do i = 1, n

        if ( a2(i) .lt. a1(i) ) then
          r8vec_gt = .true.
          return
        else if ( a1(i) .lt. a2(i) ) then
          return
        end if

      end do

      return
      end
      function r8vec_lt ( n, a1, a2 )

c*********************************************************************72
c
cc R8VEC_LT == ( A1 < A2 ) for double precision vectors.
c
c  Discussion:
c
c    The comparison is lexicographic.
c
c    A1 < A2  <=>                              A1(1) < A2(1) or
c                 ( A1(1)     == A2(1)     and A1(2) < A2(2) ) or
c                 ...
c                 ( A1(1:N-1) == A2(1:N-1) and A1(N) < A2(N)
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    22 July 2008
c
c  Author:
c
c    John Burkardt, Max Gunzburger, Janet Peterson
c
c  Parameters:
c
c    Input, integer N, the dimension of the vectors.
c
c    Input, double precision A1(N), A2(N), the vectors to be compared.
c
c    Output, logical R8VEC_LT, is TRUE if and only if A1 < A2.
c
      implicit none

      integer n

      double precision a1(n)
      double precision a2(n)
      integer i
      logical r8vec_lt

      r8vec_lt = .false.

      do i = 1, n

        if ( a1(i) .lt. a2(i) ) then
          r8vec_lt = .true.
          return
        else if ( a2(i) .lt. a1(i) ) then
          return
        end if

      end do

      return
      end
      subroutine r8vec_swap ( n, a1, a2 )

c*********************************************************************72
c
cc R8VEC_SWAP swaps two R8VEC's.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    22 July 2008
c
c  Author:
c
c    John Burkardt, Max Gunzburger, Janet Peterson
c
c  Parameters:
c
c    Input, integer N, the number of entries in the arrays.
c
c    Input/output, double precision A1(N), A2(N), the vectors to swap.
c
      implicit none

      integer n

      double precision a1(n)
      double precision a2(n)
      integer i
      double precision t

      do i = 1, n
        t     = a1(i)
        a1(i) = a2(i)
        a2(i) = t
      end do

      return
      end
      subroutine r8vec_to_bin_even3 ( ndim, nbin, a, b, c, bin )

c*********************************************************************72
c
cc R8VEC_TO_BIN_EVEN3 determines the appropriate "bin" for a R8VEC value.
c
c  Discussion:
c
c    The interval [A(I),B(I)] is divided into NBIN(I) equal subintervals
c    or bins.
c
c  Example:
c
c    NDIM = 3
c    NBIN = (/ 4, 5, 2 /),
c
c      A(1) = 1,  A(2) = 0,  A(3) = 8
c      B(1) = 17, B(2) = 20, B(3) = 10
c
c
c            8 < Z < 9                    9 < Z < 10
c
c   20 +     +     +     +     +     20 +     +     +     +     +
c        151 | 251 | 351 | 451            152 | 252 | 352 | 452
c   16 +-----+-----+-----+-----+     16 +-----+-----+-----+-----+
c        141 | 241 | 341 | 441            142 | 242 | 342 | 442
c   12 +-----+-----+-----+-----+     12 +-----+-----+-----+-----+
c        131 | 231 | 331 | 431            132 | 232 | 332 | 432
c    8 +-----+-----+-----+-----+      8 +-----+-----+-----+-----+
c        121 | 221 | 321 | 421            122 | 222 | 322 | 422
c    4 +-----+-----+-----+-----+      4 +-----+-----+-----+-----+
c        111 | 211 | 311 | 411            112 | 212 | 312 | 412
c    0 +     +     +     +     +      0 +     +     +     +     +
c      1     5     9    13    17        1     5     9    13    17
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    22 July 2008
c
c  Author:
c
c    John Burkardt, Max Gunzburger, Janet Peterson
c
c  Parameters:
c
c    Input, integer NDIM, the dimension of the space.
c
c    Input, integer NBIN(NDIM), the number of bins in each dimension.
c
c    Input, double precision A(NDIM), B(NDIM), the lower and upper limits of
c    the bin interval.  While A(I) is expected to be less than B(I), the code
c    should return useful results if A(I) is actually greater than B(I).
c
c    Input, double precision C(NDIM), a value to be placed in a bin.
c
c    Output, integer BIN(NDIM), the index of the bin to which C is assigned.
c
      implicit none

      integer ndim

      double precision a(ndim)
      double precision b(ndim)
      integer bin(ndim)
      double precision c(ndim)
      integer i
      integer nbin(ndim)

      do i = 1, ndim
        call r8_to_bin_even2 ( nbin(i), a(i), b(i), c(i), bin(i) )
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
c    John Burkardt, Max Gunzburger, Janet Peterson
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
      subroutine region_sampler ( ndim, box_min, box_max, dr, x, 
     &  random_function, reset, use_diatom, ngen )

c*********************************************************************72
c
cc REGION_SAMPLER returns a sample point in the physical region.
c
c  Discussion:
c
c    The calculations are done in NDIM dimensional space.
c
c    The physical region is enclosed in a bounding box.
c
c    A point is chosen in the bounding box, either by a uniform random
c    number generator, or from a vector Halton sequence.
c
c    If a user-supplied routine determines that this point is
c    within the physical region, this routine returns.  Otherwise,
c    a new random point is chosen.
c
c    The entries of the local vector HALTON_BASE should be distinct primes.
c    Right now, we are assuming NDIM is no greater than 3.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    22 July 2008
c
c  Author:
c
c    John Burkardt, Max Gunzburger, Janet Peterson
c
c  Reference:
c
c    John Burkardt, Max Gunzburger, Janet Peterson and Rebecca Brannon,
c    User Manual and Supporting Information for Library of Codes
c    for Centroidal Voronoi Placement and Associated Zeroth,
c    First, and Second Moment Determination,
c    Sandia National Laboratories Technical Report SAND2002-0099,
c    February 2002.
c
c  Parameters:
c
c    Input, integer NDIM, the spatial dimension.
c
c    Input, double precision BOX_MIN(NDIM), BOX_MAX(NDIM), the coordinates
c    of the two extreme corners of the bounding box.
c
c    Input, double precision DR, a tolerance used by DIATOM when testing
c    whether a point is within, outside of, or on the boundary of the
c    physical region.
c
c    Output, double precision X(NDIM), the random point.
c
c    Input, integer RANDOM_FUNCTION, specifies the random function.
c    0, uniform random numbers from F90 RANDOM_NUMBER.
c    1, Halton sequence.
c
c    Input/output, logical RESET.
c    If TRUE, the Halton sequence should be reset.
c
c    Input, logical USE_DIATOM, is TRUE if DIATOM is to be called 
c    to determine whether a point lies in the physical region; if it is
c    FALSE than a much simplified routine is used.
c
c    Output, integer NGEN, the number of points that were generated.
c    This is at least 1, but may be larger if some points were rejected.
c
      implicit none

      integer ndim

      double precision box_max(ndim)
      double precision box_min(ndim)
      double precision dr
      integer halton_base(3)
      integer halton_seed
      integer i
      integer ival
      double precision mdens
      integer ngen
      double precision r(ndim)
      integer random_function
      logical reset
      integer uniform_seed
      logical use_diatom
      double precision x(ndim)
      double precision zero
      parameter ( zero = 0.0D+00 )

      save halton_base
      save halton_seed
      save uniform_seed

      data halton_base / 2, 3, 5 /
      data halton_seed / 1 /
      data uniform_seed / 123456789 /

      ngen = 0

      if ( reset ) then
        halton_seed = 1
        reset = .false.
      end if

10    continue

        ngen = ngen + 1

        if ( 10000 .lt. ngen ) then
          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) 'REGION_SAMPLER - Fatal error!'
          write ( *, '(a,i12,a)') '  Generated ', ngen, 
     &      ' rejected points in a row.'
          write ( *, '(a)' ) 
     &      '  There may be a problem with the geometry definition.'
          write ( *, '(a)' ) ' '
          if ( random_function .eq. 0 ) then
            write ( *, '(a)' ) '  Using F90 RANDOM_NUMBER.'
          else if ( random_function .eq. 1 ) then
            write ( *, '(a)' ) '  Using Halton sequence.'
            write ( *, '(a,i12)' ) '  Current seed is ', halton_seed
          end if
          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) '  Current random value is:'
          write ( *, '(3g14.6)' ) ( r(i), i = 1, ndim )
          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) '  Current sample point is:'
          write ( *, '(3g14.6)' ) ( x(i), i = 1, ndim)
          stop
        end if
c
c  Generate a point at random using:
c  0: a uniformly distributed random value;
c  1: a Halton random value.
c
        if ( random_function .eq. 0 ) then

          call r8vec_uniform_01 ( ndim, uniform_seed, r )

        else if ( random_function .eq. 1 ) then

          call i4_to_halton_vector ( halton_seed, halton_base, 
     &      ndim, r )

          halton_seed = halton_seed + 1

        else

          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) 'REGION_SAMPLER - Fatal error!'
          write ( *, '(a,i8)' ) 
     &      '  Illegal value of RANDOM_FUNCTION = ', random_function
          stop

        end if
c
c  Determine a point in the bounding box.
c
        do i = 1, ndim
          x(i) = ( ( 1.0D+00 - r(i) ) * box_min(i) 
     &                       + r(i)   * box_max(i) )
        end do
c
c  Now determine if the point is in the region.
c
        if ( use_diatom ) then

          if ( ndim .eq. 2 ) then
            call diatom_point_test2 ( x(1), x(2), zero, dr, mdens, 
     &        ival )
          else if ( ndim .eq. 2 ) then
            call diatom_point_test2 ( x(1), x(2), x(3), dr, mdens, 
     &        ival )
          end if
c
c  Call the routine that has an analytic definition of the region:
c
        else

          call test_region ( x, ndim, ival )

        end if

        if ( ival .eq. 1 ) then
          go to 20
        end if

      go to 10

20    continue

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
c    John Burkardt, Max Gunzburger, Janet Peterson
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
      subroutine vcm ( ndim, box_min, box_max, n, cell_generator, 
     &  ns_mom, use_diatom, use_bins, region_volume_given, 
     &  region_volume, dr, nbin, bin_start, bin_last, bin_next, 
     &  cell_volume, cell_centroid, cell_moment )

c*********************************************************************72
c
cc VCM calculates Voronoi cell volumes, centroids and second moments.
c
c  Discussion:
c
c    A Monte Carlo sampling is used to estimate the quantities.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    22 July 2008
c
c  Author:
c
c    John Burkardt, Max Gunzburger, Janet Peterson
c
c  Reference:
c
c    John Burkardt, Max Gunzburger, Janet Peterson and Rebecca Brannon,
c    User Manual and Supporting Information for Library of Codes
c    for Centroidal Voronoi Placement and Associated Zeroth,
c    First, and Second Moment Determination,
c    Sandia National Laboratories Technical Report SAND2002-0099,
c    February 2002.
c
c  Parameters:
c
c    Input, integer NDIM, the spatial dimension.
c
c    Input, double precision BOX_MIN(NDIM), BOX_MAX(NDIM), the coordinates
c    of the two extreme corners of the bounding box.
c
c    Input, integer N, the number of cell generators.
c
c    Input, double precision CELL_GENERATOR(NDIM,N), the cell generators.
c
c    Input, integer NS_MOM, the number of sample points per cell generator.
c
c    Input, logical USE_DIATOM, is TRUE if DIATOM is to be called 
c    to determine whether a point lies in the physical region; if it is
c    FALSE than a much simplified routine is used.
c
c    Input, logical USE_BINS, is TRUE if the bounding box is to be divided
c    up into bins to speed up the nearest neighbor search;
c    FALSE if the nearest neighbor seach is to be done naively.
c
c    Input, logical REGION_VOLUME_GIVEN,
c    TRUE: the region volume is input in REGION_VOLUME.
c    FALSE: the region volume must be estimated by this routine.
c
c    Input/output, double precision REGION_VOLUME, the volume of the region.
c    If REGION_VOLUME_GIVEN is TRUE, then REGION_VOLUME is input by the user.
c    Otherwise, the volume is estimated and output by this routine.
c
c    Input, double precision DR, a tolerance used by DIATOM when testing
c    whether a point is within, outside of, or on the boundary of the
c    physical region.
c
c    Input, integer NBIN(3) is the number of bins to use in each direction.
c    For 3D problems, set NBIN(3) = 1.
c    For efficiency, these values should be set in such a way that the bins
c    are nearly square or cubical.
c
c    Input, integer BIN_START(NBIN(1),NBIN(2),NBIN(3)), the index of the first
c    cell generator in the bin, or -1 if none.
c
c    Input, integer BIN_LAST(NBIN(1),NBIN(2),NBIN(3)), the index of the last
c    cell generator in the bin, or -1 if none.
c
c    Input, integer BIN_NEXT(N), the index of the next cell generator in
c    the bin containing this cell generator.
c
c    Output, double precision CELL_VOLUME(N), the Voronoi cell volumes.
c
c    Output, double precision CELL_CENTROID(NDIM,N), the Voronoi cell
c    centroids.
c
c    Output, double precision CELL_MOMENT(NDIM,NDIM,N), the second moment
c    matrix for each Voronoi cell.
c
      implicit none

      integer nbin(3)
      integer ndim
      integer n

      integer bin_last(nbin(1),nbin(2),nbin(3))
      integer bin_next(n)
      integer bin_start(nbin(1),nbin(2),nbin(3))
      double precision box_max(ndim)
      double precision box_min(ndim)
      double precision box_volume
      double precision cell_centroid(ndim,n)
      double precision cell_generator(ndim,n)
      integer cell_hit(n)
      double precision cell_moment(ndim,ndim,n)
      double precision cell_volume(n)
      double precision dr
      integer i
      integer j
      integer k
      integer nearest
      integer ngen
      integer ns_mom
      integer ntries
      integer random_mom
      parameter ( random_mom = 0 )
      double precision region_volume
      double precision region_volume_estimate
      logical region_volume_given
      logical reset
      logical use_bins
      logical use_diatom
      double precision x(ndim)
c
c  Zero out the arrays.
c
      do j = 1, n
        do i = 1, ndim
          cell_centroid(i,j) = 0.0D+00
        end do
      end do

      do i = 1, n
        cell_hit(i) = 0
      end do

      do k = 1, n
        do j = 1, ndim
          do i = 1, ndim
            cell_moment(i,j,k) = 0.0D+00
          end do
        end do
      end do
c
c  Sample the region N * NS_MOM times, and keep track of which cell generator
c  is closest to each sampling point.
c
      ntries = 0
      reset = .true.

      do k = 1, n * ns_mom

        call region_sampler ( ndim, box_min, box_max, dr, x, 
     &    random_mom, reset, use_diatom, ngen )

        ntries = ntries + ngen

        if ( use_bins ) then

          if ( ndim .eq. 2 ) then

            call points_nearest_point_bins3_2d ( n, cell_generator, 
     &        nbin, box_min, box_max, bin_start, bin_last, bin_next, 
     &        x, nearest )

          else if ( ndim .eq. 3 ) then

            call points_nearest_point_bins3_3d ( n, cell_generator, 
     &        nbin, box_min, box_max, bin_start, bin_last, bin_next, 
     &        x, nearest )

          end if

        else

          call find_closest ( ndim, x, n, cell_generator, nearest )

        end if

        cell_hit(nearest) = cell_hit(nearest) + 1

        do i = 1, ndim
          cell_centroid(i,nearest) = cell_centroid(i,nearest) + x(i)
        end do

        do i = 1, ndim
          do j = 1, ndim
            cell_moment(i,j,nearest) = cell_moment(i,j,nearest) 
     &        + x(i) * x(j)
          end do
        end do

      end do
c
c  Estimate the area or volume if it was not given.
c
      box_volume = 1.0
      do i = 1, ndim
        box_volume = box_volume * ( box_max(i) - box_min(i) )
      end do

      region_volume_estimate = dble ( n * ns_mom ) * box_volume 
     &  / dble ( ntries  )

      write ( *, '(a)' ) ' '
      write ( *, '(a,g14.6)' ) 
     &  '  Volume of bounding box is     ', box_volume
      if ( region_volume_given ) then
        write ( *, '(a,g14.6)' ) 
     &    '  Given volume of region is     ', region_volume
      else
        region_volume = region_volume_estimate
      end if
      write ( *, '(a,g14.6)' ) '  Estimated volume of region is ', 
     &  region_volume_estimate
c
c  Estimate the geometric integrals for each Voronoi cell.
c
      do k = 1, n

        if ( 0 .lt. cell_hit(k) ) then

          cell_volume(k) = dble ( cell_hit(k) ) * region_volume 
     &      / dble ( n * ns_mom )

          do i = 1, ndim
            cell_centroid(i,k) = cell_centroid(i,k) 
     &        / dble ( cell_hit(k) )
          end do

          do j = 1, ndim
            do i = 1, ndim
              cell_moment(i,j,k) = cell_moment(i,j,k) 
     &          / dble ( cell_hit(k) )
            end do
          end do

          do i = 1, ndim
            do j = 1, ndim
              cell_moment(i,j,k) = cell_moment(i,j,k) 
     &          - cell_centroid(i,k) * cell_centroid(j,k)
            end do
          end do

        end if

      end do

      return
      end
