      program main

c*********************************************************************72
c
cc MESH_IO_PRB tests the MESH_IO library.
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
c  Reference:
c
c    Pascal Frey,
c    MEDIT: An interactive mesh visualization software,
c    Technical Report RT-0253,
c    Institut National de Recherche en Informatique et en Automatique,
c    03 December 2001.
c
      implicit none

      character * ( 255 ) filename

      call timestamp ( )
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'MESH_IO_PRB'
      write ( *, '(a)' ) '  FORTRAN77 version'
      write ( *, '(a)' ) '  Test the MESH_IO library.'
c
c  Create the file hexahexa_2x2x2.mesh
c
      call test01 ( )
c
c  Read and print the file hexahexa_2x2x2.mesh.
c
      filename = 'hexahexa_2x2x2.mesh'
      call test03 ( filename )
c
c  Create the file cyl248.mesh
c
      call test02 ( )
c
c  Read and print the sizes of file cyl248.mesh.
c
      filename = 'cyl248.mesh'
      call test03 ( filename )
c
c  Read and print the data in file cyl248.mesh.
c
      filename = 'cyl248.mesh'
      call test04 ( filename )
c
c  Terminate.
c
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'MESH_IO_PRB'
      write ( *, '(a)' ) '  Normal end of execution.'
      write ( *, '(a)' ) ' '
      call timestamp ( )

      return
      end
      subroutine test01 ( )

c*********************************************************************72
c
cc TEST01 creates a MESH dataset and writes it to a file.
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
      implicit none

      integer dim
      integer edges
      integer hexahedrons
      integer quadrilaterals
      integer tetrahedrons
      integer triangles
      integer vertices

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST01:'
      write ( *, '(a)' )
     &  '  Create a hexahedral mesh and write it to a file.'
c
c  Get sizes.
c
      call hexahexa_2x2x2_size ( dim, vertices, edges, triangles,
     & quadrilaterals, tetrahedrons, hexahedrons )
c
c  Call a subroutine and hope we can use automatic arrays:
c
      call test01_part2 ( dim, vertices, edges, triangles,
     & quadrilaterals, tetrahedrons, hexahedrons )

      return
      end
      subroutine test01_part2 ( dim, vertices, edges, triangles,
     & quadrilaterals, tetrahedrons, hexahedrons )

c*********************************************************************72
c
cc TEST01_PART2 creates automatic arrays and completes the test.
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
      character * ( 255 ) filename
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
c  Get the data.
c
      call hexahexa_2x2x2_data ( dim, vertices, edges, triangles,
     &  quadrilaterals, tetrahedrons, hexahedrons, vertex_coordinate,
     &  vertex_label, edge_vertex, edge_label, triangle_vertex,
     &  triangle_label, quadrilateral_vertex, quadrilateral_label,
     &  tetrahedron_vertex, tetrahedron_label, hexahedron_vertex,
     &  hexahedron_label )
c
c  Write the data.
c
      filename = 'hexahexa_2x2x2.mesh';

      call mesh_write ( filename, dim, vertices, edges, triangles,
     &  quadrilaterals, tetrahedrons, hexahedrons, vertex_coordinate,
     &  vertex_label, edge_vertex, edge_label, triangle_vertex,
     &  triangle_label, quadrilateral_vertex, quadrilateral_label,
     &  tetrahedron_vertex, tetrahedron_label, hexahedron_vertex,
     &  hexahedron_label )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' )
     &  '  Created the file "' // trim ( filename ) // '".'

      return
      end
      subroutine test02 ( )

c*********************************************************************72
c
cc TEST02 creates a MESH dataset and writes it to a file.
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
      implicit none

      integer dim
      integer edges
      integer hexahedrons
      integer quadrilaterals
      integer tetrahedrons
      integer triangles
      integer vertices

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST02:'
      write ( *, '(a)' )
     &  '  Create a hexahedral mesh and write it to a file.'
c
c  Get sizes.
c
      call cyl248_size ( dim, vertices, edges, triangles,
     & quadrilaterals, tetrahedrons, hexahedrons )
c
c  Call a subroutine and hope we can use automatic arrays:
c
      call test02_part2 ( dim, vertices, edges, triangles,
     & quadrilaterals, tetrahedrons, hexahedrons )

      return
      end
      subroutine test02_part2 ( dim, vertices, edges, triangles,
     & quadrilaterals, tetrahedrons, hexahedrons )

c*********************************************************************72
c
cc TEST02_PART2 creates automatic arrays and completes the test.
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
      character * ( 255 ) filename
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
c  Get the data.
c
      call cyl248_data ( dim, vertices, edges, triangles,
     &  quadrilaterals, tetrahedrons, hexahedrons, vertex_coordinate,
     &  vertex_label, edge_vertex, edge_label, triangle_vertex,
     &  triangle_label, quadrilateral_vertex, quadrilateral_label,
     &  tetrahedron_vertex, tetrahedron_label, hexahedron_vertex,
     &  hexahedron_label )
c
c  Write the data.
c
      filename = 'cyl248.mesh';

      call mesh_write ( filename, dim, vertices, edges, triangles,
     &  quadrilaterals, tetrahedrons, hexahedrons, vertex_coordinate,
     &  vertex_label, edge_vertex, edge_label, triangle_vertex,
     &  triangle_label, quadrilateral_vertex, quadrilateral_label,
     &  tetrahedron_vertex, tetrahedron_label, hexahedron_vertex,
     &  hexahedron_label )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' )
     &  '  Created the file "' // trim ( filename ) // '".'

      return
      end
      subroutine test03 ( filename )

c*********************************************************************72
c
cc MESH_IO_TEST03 reads and prints the sizes in a MESH dataset.
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
      implicit none

      integer dim
      integer edges
      character * ( * ) filename
      integer hexahedrons
      integer quadrilaterals
      integer tetrahedrons
      integer triangles
      integer vertices

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'MESH_IO_TEST03'
      write ( *, '(a)' ) '  Read a mesh file and print its sizes.'
c
c  Read sizes.
c
      call mesh_size_read ( filename, dim, vertices, edges, triangles,
     &  quadrilaterals, tetrahedrons, hexahedrons )
c
c  Print sizes.
c
      write ( *, '(a)' ) ' '
      write ( *, '(a)' )
     &  '  Header information for "' // trim ( filename ) // '".'

      call mesh_size_print ( dim, vertices, edges, triangles,
     &  quadrilaterals, tetrahedrons, hexahedrons )

      return
      end
      subroutine test04 ( filename )

c*********************************************************************72
c
cc MESH_IO_TEST04 reads a MESH dataset and prints its data.
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
      implicit none

      integer dim
      integer edges
      character * ( * ) filename
      integer hexahedrons
      integer quadrilaterals
      integer tetrahedrons
      integer triangles
      integer vertices

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'MESH_IO_TEST04'
      write ( *, '(a)' ) '  Read a mesh file and print its data.'
c
c  Read sizes.
c
      call mesh_size_read ( filename, dim, vertices, edges, triangles,
     &  quadrilaterals, tetrahedrons, hexahedrons )
c
c  Call a subroutine and hope we can use automatic arrays:
c
      call test04_part2 ( filename, dim, vertices, edges, triangles,
     &  quadrilaterals, tetrahedrons, hexahedrons )

      return
      end
      subroutine test04_part2 ( filename, dim, vertices, edges,
     &  triangles, quadrilaterals, tetrahedrons, hexahedrons )

c*********************************************************************72
c
cc MESH_IO_TEST04 reads a MESH dataset and prints its data.
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
      character * ( 255 ) filename
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
c  Read the data.
c
      call mesh_data_read ( filename, dim, vertices, edges, triangles,
     &  quadrilaterals, tetrahedrons, hexahedrons, vertex_coordinate,
     &  vertex_label, edge_vertex, edge_label, triangle_vertex,
     &  triangle_label, quadrilateral_vertex, quadrilateral_label,
     &  tetrahedron_vertex, tetrahedron_label, hexahedron_vertex,
     &  hexahedron_label )
c
c  Print the data.
c
      write ( *, '(a)' ) ' '
      write ( *, '(a)' )
     &  '  Data for file "' // trim ( filename ) // '".'

      call mesh_data_print ( dim, vertices, edges, triangles,
     &  quadrilaterals, tetrahedrons, hexahedrons, vertex_coordinate,
     &  vertex_label,  edge_vertex, edge_label, triangle_vertex,
     &  triangle_label, quadrilateral_vertex, quadrilateral_label,
     &  tetrahedron_vertex, tetrahedron_label, hexahedron_vertex,
     &  hexahedron_label )

      return
      end


