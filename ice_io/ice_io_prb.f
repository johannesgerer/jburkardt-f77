      program main

c*********************************************************************72
c
cc ICE_IO_PRB tests the ICE_IO library.
c
c  Discussion:
c
c    We begin by creating a file.
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
c    Russ Rew, Glenn Davis, Steve Emmerson, Harvey Davies, Ed Hartne,
c    The NETCDF User"s Guide,
c    Unidata Program Center, March 2009.
c
      implicit none

      character*255 filename

      call timestamp ( )
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'ICE_IO_PRB'
      write ( *, '(a)' ) '  FORTRAN77 version'
      write ( *, '(a)' ) '  Test the ICE_IO library.'
c
c  Create "hexahexa_2x2x2.nc"
c
      filename = 'hexahexa_2x2x2.nc'
      call test01 ( filename )
c
c  Read "hexahexa_2x2x2.nc"
c
      call test02 ( filename )
c
c  Create "cyl248.nc"
c
      filename = 'cyl248.nc'
      call test03 ( filename )
c
c  Read "cyl248.nc"
c
      call test04 ( filename )
c
c  Terminate.
c
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'ICE_IO_PRB'
      write ( *, '(a)' ) '  Normal end of execution.'
      write ( *, '(a)' ) ' '
      call timestamp ( )

      return
      end
      subroutine test01 ( filename )

c*********************************************************************72
c
cc TEST01 gets the HEXAHEXA_2X2X2 sizes and calls part 2.
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
c    Russ Rew, Glenn Davis, Steve Emmerson, Harvey Davies, Ed Hartne,
c    The NETCDF User"s Guide,
c    Unidata Program Center, March 2009.
c
c  Parameters:
c
c    Input, character*(*) FILENAME, the name of the file to be created.
c
      implicit none

      integer ( kind = 4 ) dim
      integer ( kind = 4 ) edges
      character*(*) filename
      integer ( kind = 4 ) hexahedrons
      integer ( kind = 4 ) quadrilaterals
      integer ( kind = 4 ) tetrahedrons
      integer ( kind = 4 ) triangles
      integer ( kind = 4 ) vertices

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST01:'
      write ( *, '(a)' ) '  Create an ICE grid dataset, print it,'
      write ( *, '(a)' ) '  and write it to an NETCDF file.'
c
c  Get sizes.
c
      call hexahexa_2x2x2_size ( dim, vertices, edges, triangles, 
     &  quadrilaterals, tetrahedrons, hexahedrons )
c
c  Print sizes;
c
      call size_print ( dim, vertices, edges, triangles, 
     &  quadrilaterals, tetrahedrons, hexahedrons )
c
c  I'm not going to try to use ALLOCATABLE arrays in a FORTRAN77 program.
c  But I am willing to try to get a similar result by calling a subroutine
c  and passing in sizes as parameters, relying on the subprogram to be
c  able to allocate memory at call time (NOT part of the FORTRAN77 standard).
c
      call test01_part2 ( filename, dim, vertices, edges, 
     &  triangles, quadrilaterals, tetrahedrons, hexahedrons )

      return
      end
      subroutine test01_part2 ( filename, dim, vertices, edges, 
     &  triangles, quadrilaterals, tetrahedrons, hexahedrons )

c*********************************************************************72
c
cc TEST01_PART2 creates the HEXAHEXA_2X2X2 dataset and writes it to NETCDF.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    23 October 2010
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Russ Rew, Glenn Davis, Steve Emmerson, Harvey Davies, Ed Hartne,
c    The NETCDF User"s Guide,
c    Unidata Program Center, March 2009.
c
c  Parameters:
c
c    Input, character*(*) FILENAME, the name of the file to be created.
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
      character*(*) filename
      integer hexahedron_label(hexahedrons)
      integer hexahedron_vertex(8,hexahedrons)
      integer quadrilateral_label(quadrilaterals)
      integer quadrilateral_vertex(4,quadrilaterals)
      integer tetrahedron_label(tetrahedrons)
      integer tetrahedron_vertex(4,tetrahedrons)
      integer triangle_label(triangles)
      integer triangle_vertex(3,triangles)
      double precision vertex_coordinate(3,vertices)
      integer vertex_label(vertices)
c
c  Get data.
c
      call hexahexa_2x2x2_data ( dim, vertices, edges, triangles, 
     &  quadrilaterals, tetrahedrons, hexahedrons, vertex_coordinate, 
     &  vertex_label, edge_vertex, edge_label, triangle_vertex, 
     &  triangle_label, quadrilateral_vertex, quadrilateral_label, 
     &  tetrahedron_vertex, tetrahedron_label, hexahedron_vertex, 
     &  hexahedron_label )
c
c  Print the data.
c
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Data to be written to "' 
     &  // trim ( filename ) // '".'

      call data_print ( dim, vertices, edges, triangles, 
     &  quadrilaterals, tetrahedrons, hexahedrons, vertex_coordinate, 
     &  vertex_label, edge_vertex, edge_label, triangle_vertex, 
     &  triangle_label, quadrilateral_vertex, quadrilateral_label, 
     &  tetrahedron_vertex, tetrahedron_label, hexahedron_vertex, 
     &  hexahedron_label )
c
c  Create the file.
c
      call ice_write ( filename, dim, vertices, edges, triangles,
     &  quadrilaterals, tetrahedrons, hexahedrons, vertex_coordinate, 
     &  vertex_label, edge_vertex, edge_label, triangle_vertex, 
     &  triangle_label, quadrilateral_vertex, quadrilateral_label,
     &  tetrahedron_vertex, tetrahedron_label, hexahedron_vertex, 
     &  hexahedron_label )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Created the file "' 
     &  // trim ( filename ) // '".'

      return
      end
      subroutine test02 ( filename )

c*********************************************************************72
c
cc TEST02 reads the HEXAHEXA_2X2X2 sizes and calls part 2.
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
c    Russ Rew, Glenn Davis, Steve Emmerson, Harvey Davies, Ed Hartne,
c    The NETCDF User"s Guide,
c    Unidata Program Center, March 2009.
c
c  Parameters:
c
c    Input, character*(*) FILENAME, the name of the file to be read.
c
      implicit none

      integer ( kind = 4 ) dim
      integer ( kind = 4 ) edges
      character*(*) filename
      integer ( kind = 4 ) hexahedrons
      integer ( kind = 4 ) quadrilaterals
      integer ( kind = 4 ) tetrahedrons
      integer ( kind = 4 ) triangles
      integer ( kind = 4 ) vertices

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST02:'
      write ( *, '(a)' ) 
     &  '  Read an ICE grid dataset from a NETCDF file,'
      write ( *, '(a)' ) '  and print the data.'
c
c  Get sizes.
c
      call size_read ( filename, dim, vertices, edges, triangles, 
     &  quadrilaterals, tetrahedrons, hexahedrons )
c
c  Print sizes;
c
      call size_print ( dim, vertices, edges, triangles, 
     &  quadrilaterals, tetrahedrons, hexahedrons )
c
c  Read the data within a subroutine, so we can try to get away
c  with automatic array allocation
c
      call test02_part2 ( filename, dim, vertices, edges, 
     &  triangles, quadrilaterals, tetrahedrons, hexahedrons )

      return
      end
      subroutine test02_part2 ( filename, dim, vertices, edges, 
     &  triangles, quadrilaterals, tetrahedrons, hexahedrons )

c*********************************************************************72
c
cc TEST02_PART2 reads the HEXAHEXA_2X2X2 dataset and prints it.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    29 October 2010
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Russ Rew, Glenn Davis, Steve Emmerson, Harvey Davies, Ed Hartne,
c    The NETCDF User"s Guide,
c    Unidata Program Center, March 2009.
c
c  Parameters:
c
c    Input, character*(*) FILENAME, the name of the file to be read.
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
      character*(*) filename
      integer hexahedron_label(hexahedrons)
      integer hexahedron_vertex(8,hexahedrons)
      integer quadrilateral_label(quadrilaterals)
      integer quadrilateral_vertex(4,quadrilaterals)
      integer tetrahedron_label(tetrahedrons)
      integer tetrahedron_vertex(4,tetrahedrons)
      integer triangle_label(triangles)
      integer triangle_vertex(3,triangles)
      double precision vertex_coordinate(3,vertices)
      integer vertex_label(vertices)
c
c  Read data.
c
      call data_read ( filename, dim, vertices, edges, triangles,
     &  quadrilaterals, tetrahedrons, hexahedrons, vertex_coordinate, 
     &  vertex_label, edge_vertex, edge_label, triangle_vertex, 
     &  triangle_label, quadrilateral_vertex, quadrilateral_label, 
     &  tetrahedron_vertex, tetrahedron_label, hexahedron_vertex, 
     &  hexahedron_label )
c
c  Print the data.
c
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Data read from "' 
     &  // trim ( filename ) // '".'

      call data_print ( dim, vertices, edges, triangles, 
     &  quadrilaterals, tetrahedrons, hexahedrons, vertex_coordinate, 
     &  vertex_label, edge_vertex, edge_label, triangle_vertex, 
     &  triangle_label, quadrilateral_vertex, quadrilateral_label, 
     &  tetrahedron_vertex, tetrahedron_label, hexahedron_vertex, 
     &  hexahedron_label )

      return
      end
      subroutine test03 ( filename )

c*********************************************************************72
c
cc TEST03 gets the CYL248 sizes and calls part 2.
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
c    Russ Rew, Glenn Davis, Steve Emmerson, Harvey Davies, Ed Hartne,
c    The NETCDF User"s Guide,
c    Unidata Program Center, March 2009.
c
c  Parameters:
c
c    Input, character*(*) FILENAME, the name of the file to be created.
c
      implicit none

      integer ( kind = 4 ) dim
      integer ( kind = 4 ) edges
      character*(*) filename
      integer ( kind = 4 ) hexahedrons
      integer ( kind = 4 ) quadrilaterals
      integer ( kind = 4 ) tetrahedrons
      integer ( kind = 4 ) triangles
      integer ( kind = 4 ) vertices

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST03:'
      write ( *, '(a)' ) '  Create an ICE grid dataset, print it,'
      write ( *, '(a)' ) '  and write it to an NETCDF file.'
c
c  Get sizes.
c
      call cyl248_size ( dim, vertices, edges, triangles, 
     &  quadrilaterals, tetrahedrons, hexahedrons )
c
c  Print sizes;
c
      call size_print ( dim, vertices, edges, triangles, 
     &  quadrilaterals, tetrahedrons, hexahedrons )
c
c  I'm not going to try to use ALLOCATABLE arrays in a FORTRAN77 program.
c  But I am willing to try to get a similar result by calling a subroutine
c  and passing in sizes as parameters, relying on the subprogram to be
c  able to allocate memory at call time (NOT part of the FORTRAN77 standard).
c
      call test03_part2 ( filename, dim, vertices, edges, 
     &  triangles, quadrilaterals, tetrahedrons, hexahedrons )

      return
      end
      subroutine test03_part2 ( filename, dim, vertices, edges, 
     &  triangles, quadrilaterals, tetrahedrons, hexahedrons )

c*********************************************************************72
c
cc TEST03_PART2 creates theCYL248 dataset and writes it to NETCDF.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    23 October 2010
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Russ Rew, Glenn Davis, Steve Emmerson, Harvey Davies, Ed Hartne,
c    The NETCDF User"s Guide,
c    Unidata Program Center, March 2009.
c
c  Parameters:
c
c    Input, character*(*) FILENAME, the name of the file to be created.
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
      character*(*) filename
      integer hexahedron_label(hexahedrons)
      integer hexahedron_vertex(8,hexahedrons)
      integer quadrilateral_label(quadrilaterals)
      integer quadrilateral_vertex(4,quadrilaterals)
      integer tetrahedron_label(tetrahedrons)
      integer tetrahedron_vertex(4,tetrahedrons)
      integer triangle_label(triangles)
      integer triangle_vertex(3,triangles)
      double precision vertex_coordinate(3,vertices)
      integer vertex_label(vertices)
c
c  Get data.
c
      call cyl248_data ( dim, vertices, edges, triangles, 
     &  quadrilaterals, tetrahedrons, hexahedrons, vertex_coordinate, 
     &  vertex_label, edge_vertex, edge_label, triangle_vertex, 
     &  triangle_label, quadrilateral_vertex, quadrilateral_label, 
     &  tetrahedron_vertex, tetrahedron_label, hexahedron_vertex, 
     &  hexahedron_label )
c
c  Print the data.
c
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Data to be written to "' 
     &  // trim ( filename ) // '".'

      call data_print ( dim, vertices, edges, triangles, 
     &  quadrilaterals, tetrahedrons, hexahedrons, vertex_coordinate, 
     &  vertex_label, edge_vertex, edge_label, triangle_vertex, 
     &  triangle_label, quadrilateral_vertex, quadrilateral_label, 
     &  tetrahedron_vertex, tetrahedron_label, hexahedron_vertex, 
     &  hexahedron_label )
c
c  Create the file.
c
      call ice_write ( filename, dim, vertices, edges, triangles,
     &  quadrilaterals, tetrahedrons, hexahedrons, vertex_coordinate, 
     &  vertex_label, edge_vertex, edge_label, triangle_vertex, 
     &  triangle_label, quadrilateral_vertex, quadrilateral_label,
     &  tetrahedron_vertex, tetrahedron_label, hexahedron_vertex, 
     &  hexahedron_label )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Created the file "' 
     &  // trim ( filename ) // '".'

      return
      end
      subroutine test04 ( filename )

c*********************************************************************72
c
cc TEST04 reads the CYL248 sizes and calls part 2.
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
c    Russ Rew, Glenn Davis, Steve Emmerson, Harvey Davies, Ed Hartne,
c    The NETCDF User"s Guide,
c    Unidata Program Center, March 2009.
c
c  Parameters:
c
c    Input, character*(*) FILENAME, the name of the file to be read.
c
      implicit none

      integer ( kind = 4 ) dim
      integer ( kind = 4 ) edges
      character*(*) filename
      integer ( kind = 4 ) hexahedrons
      integer ( kind = 4 ) quadrilaterals
      integer ( kind = 4 ) tetrahedrons
      integer ( kind = 4 ) triangles
      integer ( kind = 4 ) vertices

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST04:'
      write ( *, '(a)' ) 
     &  '  Read an ICE grid dataset from a NETCDF file,'
      write ( *, '(a)' ) '  and print the data.'
c
c  Get sizes.
c
      call size_read ( filename, dim, vertices, edges, triangles, 
     &  quadrilaterals, tetrahedrons, hexahedrons )
c
c  Print sizes;
c
      call size_print ( dim, vertices, edges, triangles, 
     &  quadrilaterals, tetrahedrons, hexahedrons )
c
c  Read the data within a subroutine, so we can try to get away
c  with automatic array allocation
c
      call test04_part2 ( filename, dim, vertices, edges, 
     &  triangles, quadrilaterals, tetrahedrons, hexahedrons )

      return
      end
      subroutine test04_part2 ( filename, dim, vertices, edges, 
     &  triangles, quadrilaterals, tetrahedrons, hexahedrons )

c*********************************************************************72
c
cc TEST04_PART2 reads the CYL248 dataset and prints it.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    29 October 2010
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Russ Rew, Glenn Davis, Steve Emmerson, Harvey Davies, Ed Hartne,
c    The NETCDF User"s Guide,
c    Unidata Program Center, March 2009.
c
c  Parameters:
c
c    Input, character*(*) FILENAME, the name of the file to be read.
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
      character*(*) filename
      integer hexahedron_label(hexahedrons)
      integer hexahedron_vertex(8,hexahedrons)
      integer quadrilateral_label(quadrilaterals)
      integer quadrilateral_vertex(4,quadrilaterals)
      integer tetrahedron_label(tetrahedrons)
      integer tetrahedron_vertex(4,tetrahedrons)
      integer triangle_label(triangles)
      integer triangle_vertex(3,triangles)
      double precision vertex_coordinate(3,vertices)
      integer vertex_label(vertices)
c
c  Read data.
c
      call data_read ( filename, dim, vertices, edges, triangles,
     &  quadrilaterals, tetrahedrons, hexahedrons, vertex_coordinate, 
     &  vertex_label, edge_vertex, edge_label, triangle_vertex, 
     &  triangle_label, quadrilateral_vertex, quadrilateral_label, 
     &  tetrahedron_vertex, tetrahedron_label, hexahedron_vertex, 
     &  hexahedron_label )
c
c  Print the data.
c
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Data read from "' 
     &  // trim ( filename ) // '".'

      call data_print ( dim, vertices, edges, triangles, 
     &  quadrilaterals, tetrahedrons, hexahedrons, vertex_coordinate, 
     &  vertex_label, edge_vertex, edge_label, triangle_vertex, 
     &  triangle_label, quadrilateral_vertex, quadrilateral_label, 
     &  tetrahedron_vertex, tetrahedron_label, hexahedron_vertex, 
     &  hexahedron_label )

      return
      end
