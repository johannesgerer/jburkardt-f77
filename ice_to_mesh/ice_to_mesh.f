      program main

c*********************************************************************72
c
cc ICE_TO_MESH reads ICE data from a NETCDF file and writes to a MESH file.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    23 November 2010
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
      implicit none

      integer arg_num
      integer dim
      integer edges
      character * ( 255 ) filename_mesh
      character * ( 255 ) filename_nc
      integer hexahedrons
      integer iarg
      character * ( 255 ) prefix
      integer quadrilaterals
      integer tetrahedrons
      integer triangles
      integer vertices

      call timestamp ( )
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'ICE_TO_MESH:'
      write ( *, '(a)' ) '  FORTRAN77 version'
      write ( *, '(a)' )
     &  '  Read ICE data from NETCDF file, write to MESH file.'
c
c  Get the number of command line arguments.
c
      arg_num = iargc ( )
c
c  Check the input argument.
c
      if ( 1 <= arg_num ) then

        iarg = 1
        call getarg ( iarg, prefix )

      else

        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) '  Enter the filename prefix:'
        read ( *, '(a)' ) prefix

      end if
c
c  Create the file names.
c
      filename_nc  = trim ( prefix ) // '.nc'
      filename_mesh = trim ( prefix ) // '.mesh'
c
c  Read sizes;
c
      call size_read ( filename_nc, dim, vertices, edges, triangles,
     &  quadrilaterals, tetrahedrons, hexahedrons )
c
c  Print sizes.
c
      call size_print ( dim, vertices, edges, triangles, quadrilaterals,
     &  tetrahedrons, hexahedrons )

      call ice_to_mesh_part2 ( filename_mesh, filename_nc, dim,
     &  vertices, edges, triangles, quadrilaterals, tetrahedrons,
     &  hexahedrons )
c
c  Terminate.
c
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'ICE_TO_MESH:'
      write ( *, '(a)' ) '  Normal end of execution.'
      write ( *, '(a)' ) ' '
      call timestamp ( )

      stop
      end
      subroutine ice_to_mesh_part2 ( filename_mesh, filename_nc, dim,
     &  vertices, edges, triangles, quadrilaterals, tetrahedrons,
     &  hexahedrons )

c*********************************************************************72
c
cc ICE_TO_MESH_PART2 continues the process, after creating automatic arrays.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    23 November 2010
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
      implicit none

      integer dim
      integer edges
      integer hexahedrons
      integer quadrilaterals
      integer tetrahedrons
      integer triangles
      integer vertices

      integer edge_label(edges)
      integer edge_vertex(2,edges)
      character * ( * ) filename_mesh
      character * ( * ) filename_nc
      integer hexahedron_label(hexahedrons)
      integer hexahedron_vertex(8,hexahedrons)
      integer quadrilateral_label(quadrilaterals)
      integer quadrilateral_vertex(4,quadrilaterals)
      integer tetrahedron_label(tetrahedrons)
      integer tetrahedron_vertex(4,tetrahedrons)
      integer triangle_label(triangles)
      integer triangle_vertex(3,triangles)
      double precision vertex_coordinate(dim,vertices)
      integer vertex_label(vertices)
c
c  Read the data.
c
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Reading "' // trim ( filename_nc ) // '".'

      call data_read ( filename_nc, dim, vertices, edges, triangles,
     &  quadrilaterals, tetrahedrons, hexahedrons, vertex_coordinate,
     &  vertex_label, edge_vertex, edge_label, triangle_vertex,
     &  triangle_label, quadrilateral_vertex, quadrilateral_label,
     &  tetrahedron_vertex, tetrahedron_label, hexahedron_vertex,
     &  hexahedron_label )
c
c  Print the data.
c
      if ( vertices .lt. 250 ) then

        call data_print ( dim, vertices, edges, triangles,
     &    quadrilaterals, tetrahedrons, hexahedrons, vertex_coordinate,
     &    vertex_label, edge_vertex, edge_label, triangle_vertex,
     &    triangle_label, quadrilateral_vertex, quadrilateral_label,
     &    tetrahedron_vertex, tetrahedron_label, hexahedron_vertex,
     &    hexahedron_label )

      end if
c
c  Write the data.
c
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Writing "' // trim ( filename_mesh ) // '".'

      call mesh_write ( filename_mesh, dim, vertices, edges, triangles,
     &  quadrilaterals, tetrahedrons, hexahedrons, vertex_coordinate,
     &  vertex_label, edge_vertex, edge_label,  triangle_vertex,
     &  triangle_label, quadrilateral_vertex, quadrilateral_label,
     &  tetrahedron_vertex, tetrahedron_label, hexahedron_vertex,
     &  hexahedron_label )

      write ( *, '(a)' ) '  Conversion completed.'

      return
      end
      subroutine data_print ( dim, vertices, edges, triangles,
     &  quadrilaterals, tetrahedrons, hexahedrons, vertex_coordinate,
     &  vertex_label, edge_vertex, edge_label, triangle_vertex,
     &  triangle_label, quadrilateral_vertex, quadrilateral_label,
     &  tetrahedron_vertex, tetrahedron_label, hexahedron_vertex,
     &  hexahedron_label )

c*********************************************************************72
c
cc DATA_PRINT prints the data of an ICE grid dataset.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    19 November 2010
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Pascal Frey,
c    MEDIT: An interactive mesh visualization software,
c    Technical Report RT-0253,
c    Institut National de Recherche en Informatique et en Automatique,
c    03 December 2001.
c
c  Parameters:
c
c    Input, integer DIM, the spatial dimension, which should be 2 or 3.
c
c    Input, integer VERTICES, the number of vertices.
c
c    Input, integer EDGES, the number of edges (may be 0).
c
c    Input, integer TRIANGLES, the number of triangles (may be 0).
c
c    Input, integer QUADRILATERALS, the number of quadrilaterals
c    (may be 0).
c
c    Input, integer TETRAHEDRONS, the number of tetrahedrons
c    (may be 0).
c
c    Input, integer HEXAHEDRONS, the number of hexahedrons
c    (may be 0).
c
c    Input, double VERTEX_COORDINATE(DIM,VERTICES), the coordinates
c    of each vertex.
c
c    Input, integer VERTEX_LABEL(VERTICES), a label for
c    each vertex.
c
c    Input, integer EDGE_VERTEX(2,EDGES), the vertices that form
c    each edge.
c
c    Input, integer EDGE_LABEL(EDGES), a label for each edge.
c
c    Input, integer TRIANGLE_VERTEX(3,TRIANGLES), the vertices
c    that form each triangle.
c
c    Input, integer TRIANGLE_LABEL(TRIANGLES), a label for
c    each triangle.
c
c    Input, integer QUADRILATERAL_VERTEX(4,QUADRILATERALS),
c    the vertices that form each quadrilateral.
c
c    Input, integer QUADRILATERAL_LABEL(QUADRILATERALS), a label
c    for each quadrilateral.
c
c    Input, integer TETRAHEDRON_VERTEX(4,TETRAHEDRONS), the
c    vertices that form each tetrahedron.
c
c    Input, integer TETRAHEDRON_LABEL(TETRAHEDRONS), a label for
c    each tetrahedron.
c
c    Input, integer HEXAHEDRON_VERTEX(8,HEXAHEDRONS), the vertices
c    that form each hexahedron.
c
c    Input, integer HEXAHEDRON_LABEL(HEXAHEDRONS), a label for
c    each hexahedron.
c
      implicit none

      integer dim
      integer edges
      integer hexahedrons
      integer quadrilaterals
      integer tetrahedrons
      integer triangles
      integer vertices

      integer edge_label(edges)
      integer edge_vertex(2,edges)
      integer hexahedron_label(hexahedrons)
      integer hexahedron_vertex(8,hexahedrons)
      integer i
      integer j
      integer quadrilateral_label(quadrilaterals)
      integer quadrilateral_vertex(4,quadrilaterals)
      integer tetrahedron_label(tetrahedrons)
      integer tetrahedron_vertex(4,tetrahedrons)
      integer triangle_label(triangles)
      integer triangle_vertex(3,triangles)
      real ( kind = 8 ) vertex_coordinate(dim,vertices)
      integer vertex_label(vertices)

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Vertices:'
      write ( *, '(a)' ) ' '
      if ( dim == 2 ) then
        do j = 1, vertices
          write ( *, '(2(2x,f10.4),2x,''('',i4,'')'')' )
     &      vertex_coordinate(1:dim,j), vertex_label(j)
        end do
      else if ( dim == 3 ) then
        do j = 1, vertices
          write ( *, '(3(2x,f10.4),2x,''('',i4,'')'')' )
     &      vertex_coordinate(1:dim,j), vertex_label(j)
        end do
      end if

      if ( 0 .lt. edges ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) '  Edges:'
        write ( *, '(a)' ) ' '
        do j = 1, edges
          write ( *, '(2(2x,i8),2x,''('',i4,'')'')' )
     &      edge_vertex(1:2,j), edge_label(j)
        end do
      end if

      if ( 0 .lt. triangles ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) '  Triangles:'
        write ( *, '(a)' ) ' '
        do j = 1, triangles
          write ( *, '(3(2x,i8),2x,''('',i4,'')'')' )
     &      triangle_vertex(1:3,j), triangle_label(j)
        end do
      end if

      if ( 0 .lt. quadrilaterals ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) '  Quadrilaterals:'
        write ( *, '(a)' ) ' '
        do j = 1, quadrilaterals
          write ( *, '(4(2x,i8),2x,''('',i4,'')'')' )
     &      quadrilateral_vertex(1:4,j), quadrilateral_label(j)
        end do
      end if

      if ( 0 .lt. tetrahedrons ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) '  Tetrahedrons:'
        write ( *, '(a)' ) ' '
        do j = 1, tetrahedrons
          write ( *, '(4(2x,i8),2x,''('',i4,'')'')' )
     &      tetrahedron_vertex(1:4,j), tetrahedron_label(j)
        end do
      end if

      if ( 0 .lt. hexahedrons ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) '  Hexahedrons:'
        write ( *, '(a)' ) ' '
        do j = 1, hexahedrons
          write ( *, '(8(2x,i8),2x,''('',i4,'')'')' )
     &      hexahedron_vertex(1:8,j), hexahedron_label(j)
        end do
      end if

      return
      end
      subroutine data_read ( filename, dim, vertices, edges, triangles,
     &  quadrilaterals, tetrahedrons, hexahedrons, vertex_coordinate,
     &  vertex_label, edge_vertex, edge_label, triangle_vertex,
     &  triangle_label, quadrilateral_vertex, quadrilateral_label,
     &  tetrahedron_vertex, tetrahedron_label, hexahedron_vertex,
     &  hexahedron_label )

c*********************************************************************72
c
cc DATA_READ reads ICE data from a NETCDF file.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    28 October 2010
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Pascal Frey,
c    MEDIT: An interactive mesh visualization software,
c    Technical Report RT-0253,
c    Institut National de Recherche en Informatique et en Automatique,
c    03 December 2001.
c
c  Parameters:
c
c    Input, character*(*) FILENAME, the name of the file to be read.
c    Ordinarily, the name should include the extension '.nc'.
c
c    Input, integer DIM, the spatial dimension, which should be 2 or 3.
c
c    Input, integer VERTICES, the number of vertices.
c
c    Input, integer EDGES, the number of edges (may be 0).
c
c    Input, integer TRIANGLES, the number of triangles (may be 0).
c
c    Input, integer QUADRILATERALS, the number of quadrilaterals (may be 0).
c
c    Input, integer TETRAHEDRONS, the number of tetrahedrons (may be 0).
c
c    Input, integer HEXAHEDRONS, the number of hexahedrons (may be 0).
c
c    Output, double precision VERTEX_COORDINATE(DIM,VERTICES), the coordinates
c   of each vertex.
c
c    Output, int VERTEX_LABEL(VERTICES), a label for each vertex.
c
c    Output, int EDGE_VERTEX(2,EDGES), the vertices that form each edge.
c
c    Output, int EDGE_LABEL(EDGES), a label for each edge.
c
c    Output, int TRIANGLE_VERTEX(3,TRIANGLES), the vertices that form
c    each triangle.
c
c    Output, int TRIANGLE_LABEL(TRIANGLES), a label for each triangle.c
c
c    Output, int QUADRILATERAL_VERTEX(4,QUADRILATERALS), the vertices that
c    form each quadrilateral.
c
c    Output, int QUADRILATERAL_LABEL(QUADRILATERALS), a label for
c    each quadrilateral.
c
c    Output, int TETRAHEDRON_VERTEX(4,TETRAHEDRONS), the vertices that
c    form each tetrahedron.
c
c    Output, int TETRAHEDRON_LABEL(TETRAHEDRONS), a label for
c    each tetrahedron.
c
c    Output, int HEXAHEDRON_VERTEX(8,HEXAHEDRONS), the vertices that form
c    each hexahedron.
c
c    Output, int HEXAHEDRON_LABEL(HEXAHEDRONS), a label for each hexahedron.
c
      implicit none

      include 'netcdf.inc'

      integer dim
      integer edges
      integer hexahedrons
      integer quadrilaterals
      integer tetrahedrons
      integer triangles
      integer vertices

      integer dim_dimension
      integer dim_edges
      integer dim_eight
      integer dim_four
      integer dim_hexahedrons
      integer dim_quadrilaterals
      integer dim_tetrahedrons
      integer dim_three
      integer dim_triangles
      integer dim_two
      integer dim_vertices
      integer dimid
      integer dimids(2)
      integer edge_label(edges)
      integer edge_vertex(2,edges)
      character*(*) filename
      integer hexahedron_label(hexahedrons)
      integer hexahedron_vertex(8,hexahedrons)
      integer i
      integer mode
      integer ncid
      integer ndims
      integer quadrilateral_label(quadrilaterals)
      integer quadrilateral_vertex(4,quadrilaterals)
      integer status
      integer tetrahedron_label(tetrahedrons)
      integer tetrahedron_vertex(4,tetrahedrons)
      integer triangle_label(triangles)
      integer triangle_vertex(3,triangles)
      integer var_edge_label
      integer var_edge_vertex
      integer var_hexahedron_label
      integer var_hexahedron_vertex
      integer var_quadrilateral_label
      integer var_quadrilateral_vertex
      integer var_tetrahedron_label
      integer var_tetrahedron_vertex
      integer var_triangle_label
      integer var_triangle_vertex
      integer var_vertex_coordinate
      integer var_vertex_label
      integer varid
      double precision vertex_coordinate(dim,vertices)
      integer vertex_label(vertices)
      integer xtype
c
c  Open the file.
c
      mode = NF_NOCLOBBER
      status = nf_open ( filename, mode, ncid )
c
c  Vertices.
c
      status = nf_inq_varid ( ncid, 'Vertex_Coordinate', varid )
      status = nf_get_var_double ( ncid, varid, vertex_coordinate )

      status = nf_inq_varid ( ncid, 'Vertex_Label', varid )
      status = nf_get_var_int ( ncid, varid, vertex_label )

c  Edges.
c
      if ( 0 .lt. edges ) then
        status = nf_inq_varid ( ncid, 'Edge_Vertex', varid )
        status = nf_get_var_int ( ncid, varid, edge_vertex )

        status = nf_inq_varid ( ncid, 'Edge_Label', varid )
        status = nf_get_var_int ( ncid, varid, edge_label )
      end if
c
c  Triangles.
c
      if ( 0 .lt. triangles ) then
        status = nf_inq_varid ( ncid, 'Triangle_Vertex', varid )
        status = nf_get_var_int ( ncid, varid, triangle_vertex )

        status = nf_inq_varid ( ncid, 'Triangle_Label', varid )
        status = nf_get_var_int ( ncid, varid, triangle_label )
      end if
c
c  Quadrilaterals.
c
      if ( 0 .lt. quadrilaterals ) then
        status = nf_inq_varid ( ncid, 'Quadrilateral_Vertex', varid )
        status = nf_get_var_int ( ncid, varid, quadrilateral_vertex )

        status = nf_inq_varid ( ncid, 'Quadrilateral_Label', varid )
        status = nf_get_var_int ( ncid, varid, quadrilateral_label )
      end if
c
c  Tetrahedrons.
c
      if ( 0 .lt. tetrahedrons ) then
        status = nf_inq_varid ( ncid, 'Tetrahedron_Vertex', varid )
        status = nf_get_var_int ( ncid, varid, tetrahedron_vertex )

        status = nf_inq_varid ( ncid, 'Tetrahedron_Label', varid )
        status = nf_get_var_int ( ncid, varid, tetrahedron_label )
      end if
c
c  Hexahedrons.
c
      if ( 0 .lt. hexahedrons ) then
        status = nf_inq_varid ( ncid, 'Hexahedron_Vertex', varid )
        status = nf_get_var_int ( ncid, varid, hexahedron_vertex )

        status = nf_inq_varid ( ncid, 'Hexahedron_Label', varid )
        status = nf_get_var_int ( ncid, varid, hexahedron_label )
      end if
c
c  Close the file.
c
      status = nf_close ( ncid )

      return
      end
      subroutine get_unit ( unit )

c*********************************************************************72
c
cc GET_UNIT returns a free FORTRAN unit number.
c
c  Discussion:
c
c    A "free" FORTRAN unit number is an integer between 1 and 99 which
c    is not currently associated with an I/O device.  A free FORTRAN unit
c    number is needed in order to open a file with the OPEN command.
c
c    If UNIT = 0, then no free FORTRAN unit could be found, although
c    all 99 units were checked (except for units 5, 6 and 9, which
c    are commonly reserved for console I/O).
c
c    Otherwise, UNIT is an integer between 1 and 99, representing a
c    free FORTRAN unit.  Note that GET_UNIT assumes that units 5 and 6
c    are special, and will never return those values.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    26 October 2008
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, integer UNIT, the free unit number.
c
      implicit none

      integer i
      integer unit

      unit = 0

      do i = 1, 99

        if ( i .ne. 5 .and. i .ne. 6 .and. i .ne. 9 ) then

          open ( unit = i, err = 10, status = 'scratch' )
          close ( unit = i )

          unit = i

          return
        end if

10      continue

      end do

      return
      end
      subroutine mesh_write ( filename, dim, vertices, edges,
     &  triangles, quadrilaterals, tetrahedrons, hexahedrons,
     &  vertex_coordinate, vertex_label, edge_vertex, edge_label,
     &  triangle_vertex, triangle_label, quadrilateral_vertex,
     &  quadrilateral_label, tetrahedron_vertex, tetrahedron_label,
     &  hexahedron_vertex, hexahedron_label )

c*********************************************************************72
c
cc MESH_WRITE writes sizes and data to a MESH file.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    03 December 2010
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Pascal Frey,
c    MEDIT: An interactive mesh visualization software,
c    Technical Report RT-0253,
c    Institut National de Recherche en Informatique et en Automatique,
c    03 December 2001.
c
c  Parameters:
c
c    Input, character * ( * ) FILENAME, the name of the file to be created.
c    Ordinarily, the name should include the extension ".mesh".
c
c    Input, integer DIM, the spatial dimension, which should be 2 or 3.
c
c    Input, integer VERTICES, the number of vertices.
c
c    Input, double precision VERTEX_COORDINATE(DIM,VERTICES), the coordinates
c    of each vertex.
c
c    Input, integer VERTEX_LABEL(VERTICES), a label for each vertex.
c
c    Input, integer EDGES, the number of edges (may be 0).
c
c    Input, integer EDGE_VERTEX(2,EDGES), the vertices that form each edge.
c
c    Input, integer EDGE_LABEL(EDGES), a label for each edge.
c
c    Input, integer TRIANGLES, the number of triangles (may be 0).
c
c    Input, integer TRIANGLE_VERTEX(3,TRIANGLES), the vertices that form
c    each triangle.
c
c    Input, integer TRIANGLE_LABEL(TRIANGLES), a label for each triangle.
c
c    Input, integer QUADRILATERALS, the number of quadrilaterals (may be 0).
c
c    Input, integer QUADRILATERAL_VERTEX(4,QUADRILATERALS), the vertices that
c    form each quadrilateral.
c
c    Input, integer QUADRILATERAL_LABEL(QUADRILATERALS), a label for
c    each quadrilateral.
c
c    Input, integer TETRAHEDRONS, the number of tetrahedrons (may be 0).
c
c    Input, integer TETRAHEDRON_VERTEX(4,TETRAHEDRONS), the vertices that
c    form each tetrahedron.
c
c    Input, integer TETRAHEDRON_LABEL(TETRAHEDRONS), a label for
c    each tetrahedron.
c
c    Input, integer HEXAHEDRONS, the number of hexahedrons (may be 0).
c
c    Input, integer HEXAHEDRON_VERTEX(8,HEXAHEDRONS), the vertices that form
c    each hexahedron.
c
c    Input, integer HEXAHEDRON_LABEL(HEXAHEDRONS), a label for each hexahedron.
c
      implicit none

      integer dim
      integer edges
      integer hexahedrons
      integer quadrilaterals
      integer tetrahedrons
      integer triangles
      integer vertices

      integer edge_label(edges)
      integer edge_vertex(2,edges)
      character * ( * ) filename
      integer fileunit
      integer hexahedron_label(hexahedrons)
      integer hexahedron_vertex(8,hexahedrons)
      integer i
      integer ios
      integer j
      integer quadrilateral_label(quadrilaterals)
      integer quadrilateral_vertex(4,quadrilaterals)
      integer tetrahedron_label(tetrahedrons)
      integer tetrahedron_vertex(4,tetrahedrons)
      integer triangle_label(triangles)
      integer triangle_vertex(3,triangles)
      double precision vertex_coordinate(dim,vertices)
      integer vertex_label(vertices)
c
c  Open the file.
c
      call get_unit ( fileunit )

      open ( unit = fileunit, file = filename, status = 'replace',
     &  iostat = ios )

      if ( ios .ne. 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'MESH_WRITE - Fatal errorc'
        write ( *, '(a)' ) '  Could not open file.'
        stop
      end if

      write ( fileunit, '(a)' ) 'MeshVersionFormatted 1'
      write ( fileunit, '(a)' ) '#  Created by mesh_write.f'
c
c  Dimension information.
c
      write ( fileunit, '(a)' ) ' '
      write ( fileunit, '(a)' ) 'Dimension'
      write ( fileunit, '(i8)' ) dim
c
c  Vertices.
c
      write ( fileunit, '(a)' ) ' '
      write ( fileunit, '(a)' ) 'Vertices'
      write ( fileunit, '(i8)' ) vertices
      if ( dim == 2 ) then
        do j = 1, vertices
          write ( fileunit, '(2(2x,f10.6),2x,i8)' )
     &      ( vertex_coordinate(i,j), i = 1, dim ), vertex_label(j)
        end do
      else if ( dim == 3 ) then
        do j = 1, vertices
          write ( fileunit, '(3(2x,f10.6),2x,i8)' )
     &      ( vertex_coordinate(i,j), i = 1, dim ), vertex_label(j)
        end do
      end if
c
c  Edges.
c
      if ( 0 .lt. edges ) then
        write ( fileunit, '(a)' ) ' '
        write ( fileunit, '(a)' ) 'Edges'
        write ( fileunit, '(i8)' ) edges
        do j = 1, edges
          write ( fileunit, '(2(2x,i8),2x,i8)' )
     &      ( edge_vertex(i,j), i = 1, 2 ), edge_label(j)
        end do
      end if
c
c  Triangles.
c
      if ( 0 .lt. triangles ) then
        write ( fileunit, '(a)' ) ' '
        write ( fileunit, '(a)' ) 'Triangles'
        write ( fileunit, '(i8)' ) triangles
        do j = 1, triangles
          write ( fileunit, '(3(2x,i8),2x,i8)' )
     &      ( triangle_vertex(i,j), i = 1, 3 ), triangle_label(j)
        end do
      end if
c
c  Quadrilaterals.
c
      if ( 0 .lt. quadrilaterals ) then
        write ( fileunit, '(a)' ) ' '
        write ( fileunit, '(a)' ) 'Quadrilaterals'
        write ( fileunit, '(i8)' ) quadrilaterals
        do j = 1, quadrilaterals
          write ( fileunit, '(4(2x,i8),2x,i8)' )
     &      ( quadrilateral_vertex(i,j), i = 1, 4 ),
     &      quadrilateral_label(j)
        end do
      end if
c
c  Tetrahedron.
c
      if ( 0 .lt. tetrahedrons ) then
        write ( fileunit, '(a)' ) ' '
        write ( fileunit, '(a)' ) 'Tetrahedra'
        write ( fileunit, '(i8)' ) tetrahedrons
        do j = 1, tetrahedrons
          write ( fileunit, '(4(2x,i8),2x,i8)' )
     &      ( tetrahedron_vertex(i,j), i = 1, 4 ), tetrahedron_label(j)
        end do
      end if
c
c  Hexahedron.
c
      if ( 0 .lt. hexahedrons ) then
        write ( fileunit, '(a)' ) ' '
        write ( fileunit, '(a)' ) 'Hexahedra'
        write ( fileunit, '(i8)' ) hexahedrons
        do j = 1, hexahedrons
          write ( fileunit, '(8(2x,i8),2x,i8)' )
     &      ( hexahedron_vertex(i,j), i = 1, 8 ), hexahedron_label(j)
        end do
      end if
c
c  End.
c
      write ( fileunit, '(a)' ) ' '
      write ( fileunit, '(a)' ) 'End'

      close ( unit = fileunit )

      return
      end
      subroutine size_print ( dim, vertices, edges, triangles,
     &  quadrilaterals, tetrahedrons, hexahedrons )

c*********************************************************************72
c
cc SIZE_PRINT prints the sizes of an ICE grid dataset.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    28 October 2010
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Pascal Frey,
c    MEDIT: An interactive mesh visualization software,
c    Technical Report RT-0253,
c    Institut National de Recherche en Informatique et en Automatique,
c    03 December 2001.
c
c  Parameters:
c
c    Input, integer DIM, the spatial dimension, which should be 2 or 3.
c
c    Input, integer VERTICES, the number of vertices.
c
c    Input, integer EDGES, the number of edges (may be 0).
c
c    Input, integer TRIANGLES, the number of triangles (may be 0).
c
c    Input, integer QUADRILATERALS, the number of quadrilaterals
c    (may be 0).
c
c    Input, integer TETRAHEDRONS, the number of tetrahedrons
c    (may be 0).
c
c    Input, integer HEXAHEDRONS, the number of hexahedrons
c    (may be 0).
c
      implicit none

      integer dim
      integer edges
      integer hexahedrons
      integer quadrilaterals
      integer tetrahedrons
      integer triangles
      integer vertices

      write ( *, '(a)' ) ' '
      write ( *, '(a,i8)' ) '  Number of dimensions = ', dim
      write ( *, '(a,i8)' ) '  Number of vertices = ', vertices
      write ( *, '(a,i8)' ) '  Number of edges = ', edges
      write ( *, '(a,i8)' ) '  Number of triangles = ', triangles
      write ( *, '(a,i8)' )
     &  '  Number of quadrilaterals = ', quadrilaterals
      write ( *, '(a,i8)' ) '  Number of tetrahedrons = ', tetrahedrons
      write ( *, '(a,i8)' ) '  Number of hexahedrons = ', hexahedrons

      return
      end
      subroutine size_read ( filename, dim, vertices, edges,
     &  triangles, quadrilaterals, tetrahedrons, hexahedrons )

c*********************************************************************72
c
cc SIZE_READ reads ICE sizes from a NETCDF file.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    28 October 2010
c
c  Author:
c
c   John Burkardt
c
c  Reference:
c
c    Pascal Frey,
c    MEDIT: An interactive mesh visualization software,
c    Technical Report RT-0253,
c    Institut National de Recherche en Informatique et en Automatique,
c    03 December 2001.
c
c  Parameters:
c
c    Input, character*(*) FILENAME, the name of the file to be read.
c    Ordinarily, the name should include the extension '.nc'.
c
c    Output, integer DIM, the spatial dimension, which should be 2 or 3.
c
c    Output, integer VERTICES, the number of vertices.
c
c    Output, integer EDGES, the number of edges (may be 0).
c
c    Output, integer TRIANGLES, the number of triangles (may be 0).
c
c    Output, integer QUADRILATERALS, the number of quadrilaterals
c    (may be 0).
c
c    Output, integer TETRAHEDRONS, the number of tetrahedrons
c    (may be 0).
c
c    Output, integer HEXAHEDRONS, the number of hexahedrons
c    (may be 0).
c
      implicit none

      include 'netcdf.inc'

      integer dim
      integer dimid
      integer edges
      character*(*) filename
      integer hexahedrons
      integer mode
      integer ncid
      integer quadrilaterals
      integer status
      integer tetrahedrons
      integer triangles
      integer vertices
c
c  Initialize everything to nothing.
c
      dim = 0
      vertices = 0
      edges = 0
      triangles = 0
      quadrilaterals = 0
      tetrahedrons = 0
      hexahedrons = 0
c
c  Open the file.
c
      mode = NF_NOWRITE
      status = nf_open ( filename, mode, ncid )
c
c  Get the dimension information.
c
      status = nf_inq_dimid ( ncid, 'Dimension', dimid )
      if ( status == NF_NOERR ) then
        status = nf_inq_dimlen ( ncid, dimid, dim )
      end if

      status = nf_inq_dimid ( ncid, 'Vertices', dimid )
      if ( status == NF_NOERR ) then
        status = nf_inq_dimlen ( ncid, dimid, vertices )
      end if

      status = nf_inq_dimid ( ncid, 'Edges', dimid )
      if ( status == NF_NOERR ) then
        status = nf_inq_dimlen ( ncid, dimid, edges )
      end if

      status = nf_inq_dimid ( ncid, 'Triangles', dimid )
      if ( status == NF_NOERR ) then
        status = nf_inq_dimlen ( ncid, dimid, triangles )
      end if

      status = nf_inq_dimid ( ncid, 'Quadrilaterals', dimid )
      if ( status == NF_NOERR ) then
        status = nf_inq_dimlen ( ncid, dimid, quadrilaterals )
      end if

      status = nf_inq_dimid ( ncid, 'Tetrahedra', dimid )
      if ( status == NF_NOERR ) then
        status = nf_inq_dimlen ( ncid, dimid, tetrahedrons )
      end if

      status = nf_inq_dimid ( ncid, 'Tetrahedrons', dimid )
      if ( status == NF_NOERR ) then
        status = nf_inq_dimlen ( ncid, dimid, tetrahedrons )
      end if

      status = nf_inq_dimid ( ncid, 'Hexahedra', dimid )
      if ( status == NF_NOERR ) then
        status = nf_inq_dimlen ( ncid, dimid, hexahedrons )
      end if

      status = nf_inq_dimid ( ncid, 'Hexahedrons', dimid )
      if ( status == NF_NOERR ) then
        status = nf_inq_dimlen ( ncid, dimid, hexahedrons )
      end if
c
c  Close the file.
c
      status = nf_close ( ncid )

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
