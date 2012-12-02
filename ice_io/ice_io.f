      subroutine cyl248_data ( dim, vertices, edges, triangles, 
     &  quadrilaterals, tetrahedrons, hexahedrons, vertex_coordinate, 
     &  vertex_label, edge_vertex, edge_label, triangle_vertex, 
     &  triangle_label, quadrilateral_vertex, quadrilateral_label, 
     &  tetrahedron_vertex, tetrahedron_label, hexahedron_vertex, 
     &  hexahedron_label )

c*********************************************************************72
c
cc CYL248_DATA defines the data for a 3D tetrahedral mesh.
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
c    Input, integer DIM, the spatial dimension, which should be 3.
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
c    Output, double precision VERTEX_COORDINATE(3,VERTICES), the XYZ coordinates
c    of each vertex.
c
c    Output, integer VERTEX_LABEL(VERTICES), a label for 
c    each vertex.
c
c    Output, integer EDGE_VERTEX(2,EDGES), the vertices that form 
c    each edge.
c
c    Output, integer EDGE_LABEL(EDGES), a label for each edge.
c
c    Output, integer TRIANGLE_VERTEX(3,TRIANGLES), the vertices 
c    that form each triangle.
c
c    Output, integer TRIANGLE_LABEL(TRIANGLES), a label for 
c    each triangle.
c
c    Output, integer QUADRILATERAL_VERTEX(4,QUADRILATERALS), the 
c    vertices that form each quadrilateral.
c
c    Output, integer QUADRILATERAL_LABEL(QUADRILATERALS), a label 
c    for each quadrilateral.
c
c    Output, integer TETRAHEDRON_VERTEX(4,TETRAHEDRONS), the 
c    vertices that form each tetrahedron.
c
c    Output, integer TETRAHEDRON_LABEL(TETRAHEDRONS), a label for 
c    each tetrahedron.
c
c    Output, integer HEXAHEDRON_VERTEX(8,HEXAHEDRONS), the 
c    vertices that form each hexahedron.
c
c    Output, integer HEXAHEDRON_LABEL(HEXAHEDRONS), a label for 
c    each hexahedron.
c
      implicit none

      integer dim
      integer dim_save 
      parameter ( dim_save = 3 )
      integer edges
      integer hexahedrons
      integer i
      integer quadrilaterals
      integer tetrahedrons
      integer tetrahedrons_save
      parameter ( tetrahedrons_save = 248 )
      integer triangles
      integer triangles_save
      parameter ( triangles_save = 154 )
      integer vertices
      integer vertices_save
      parameter ( vertices_save = 92 )

      integer edge_label(edges)
      integer edge_vertex(2,edges)
      integer hexahedron_label(hexahedrons)
      integer hexahedron_vertex(8,hexahedrons)
      integer quadrilateral_label(quadrilaterals)
      integer quadrilateral_vertex(4,quadrilaterals)
      integer tetrahedron_label(tetrahedrons)
      integer tetrahedron_vertex(4,tetrahedrons)
      integer tetrahedron_vertex_save(4,tetrahedrons_save) 
      integer triangle_label(triangles)
      integer triangle_label_save(triangles_save)
      integer triangle_vertex(3,triangles)
      integer triangle_vertex_save(3,triangles_save) 
      double precision vertex_coordinate(dim,vertices)
      double precision vertex_coordinate_save(3,vertices_save) 
      integer vertex_label(vertices)
      integer vertex_label_save(vertices_save)

      save tetrahedron_vertex_save
      save triangle_label_save
      save triangle_vertex_save
      save vertex_coordinate_save
      save vertex_label_save

      data tetrahedron_vertex_save /
     &  23, 1, 9, 8, 
     &  27, 9, 23, 1, 
     &  26, 8, 23, 9, 
     &  26, 9, 7, 8, 
     &  2, 9, 27, 1, 
     &  26, 9, 10, 7, 
     &  26, 28, 7, 10, 
     &  11, 29, 3, 2, 
     &  7, 6, 10, 28, 
     &  10, 6, 31, 28, 
     &  11, 29, 30, 3, 
     &  11, 30, 4, 3, 
     &  11, 30, 32, 4, 
     &  10, 6, 5, 31, 
     &  11, 5, 4, 32, 
     &  19, 33, 34, 20,  
     &  39, 22, 40, 16, 
     &  39, 17, 36, 22, 
     &  39, 22, 16, 17, 
     &  40, 22, 15, 16, 
     &  12, 19, 20, 33, 
     &  19, 20, 34, 18, 
     &  12, 33, 20, 35, 
     &  38, 37, 14, 21, 
     &  36, 22, 17, 18, 
     &  38, 14, 15, 21, 
     &  13, 14, 37, 21, 
     &  12, 20, 13, 35, 
     &  80, 32, 11, 30, 
     &  80, 28, 10, 31, 
     &  80, 31, 59, 28, 
     &  80, 58, 57, 26, 
     &  80, 28, 58, 26, 
     &  80, 59, 58, 28, 
     &  80, 28, 26, 10, 
     &  80, 10, 26, 9, 
     &  80, 9, 11, 10, 
     &  80, 9, 26, 23, 
     &  80, 23, 26, 57, 
     &  80, 23, 27, 9, 
     &  80, 23, 56, 27, 
     &  80, 30, 11, 29, 
     &  80, 5, 10, 11, 
     &  80, 5, 11, 32, 
     &  80, 5, 32, 31, 
     &  80, 31, 10, 5, 
     &  80, 2, 11, 9, 
     &  80, 29, 11, 2, 
     &  80, 2, 9, 27, 
     &  80, 27, 29, 2, 
     &  81, 40, 39, 22, 
     &  81, 22, 39, 36, 
     &  81, 18, 36, 34, 
     &  81, 34, 20, 18, 
     &  81, 22, 36, 18, 
     &  81, 20, 22, 18, 
     &  81, 37, 38, 21, 
     &  81, 20, 33, 35, 
     &  81, 13, 21, 20, 
     &  81, 13, 20, 35, 
     &  81, 13, 37, 21, 
     &  81, 35, 37, 13, 
     &  81, 20, 21, 22, 
     &  81, 34, 33, 20, 
     &  81, 21, 38, 15, 
     &  81, 38, 40, 15, 
     &  81, 22, 21, 15, 
     &  81, 15, 40, 22, 
     &  82, 60, 74, 59, 
     &  82, 74, 25, 59, 
     &  82, 73, 72, 58, 
     &  82, 25, 73, 58, 
     &  82, 59, 25, 58, 
     &  82, 58, 72, 57, 
     &  82, 57, 80, 58, 
     &  82, 58, 80, 59, 
     &  83, 71, 79, 70, 
     &  83, 70, 76, 78, 
     &  83, 79, 76, 70, 
     &  83, 79, 60, 76, 
     &  83, 82, 60, 74, 
     &  84, 54, 64, 55, 
     &  84, 64, 65, 55, 
     &  84, 65, 63, 55, 
     &  84, 65, 71, 63, 
     &  85, 29, 62, 30, 
     &  85, 80, 29, 30, 
     &  85, 29, 61, 62, 
     &  85, 78, 83, 76, 
     &  85, 78, 76, 30, 
     &  85, 62, 78, 30, 
     &  85, 76, 83, 60, 
     &  85, 76, 32, 30, 
     &  85, 32, 80, 30, 
     &  85, 32, 76, 60, 
     &  85, 27, 61, 29, 
     &  85, 80, 27, 29, 
     &  85, 83, 82, 60, 
     &  85, 77, 78, 62, 
     &  85, 60, 82, 59, 
     &  85, 59, 82, 80, 
     &  85, 32, 60, 31, 
     &  85, 80, 32, 31, 
     &  85, 60, 59, 31, 
     &  85, 59, 80, 31, 
     &  86, 51, 68, 52, 
     &  86, 69, 68, 51, 
     &  86, 68, 67, 52, 
     &  86, 52, 67, 53, 
     &  86, 67, 66, 53, 
     &  86, 53, 66, 54, 
     &  87, 50, 70, 49, 
     &  87, 71, 70, 50, 
     &  87, 63, 71, 50, 
     &  87, 63, 84, 71, 
     &  87, 70, 69, 49, 
     &  87, 71, 83, 70, 
     &  87, 49, 69, 51, 
     &  87, 69, 86, 51, 
     &  88, 64, 66, 73, 
     &  88, 72, 73, 66, 
     &  88, 72, 82, 73, 
     &  88, 24, 72, 66, 
     &  88, 64, 73, 25, 
     &  88, 73, 82, 25, 
     &  88, 66, 64, 54, 
     &  88, 84, 54, 64, 
     &  88, 87, 86, 84, 
     &  88, 67, 24, 66, 
     &  88, 66, 86, 67, 
     &  88, 64, 25, 65, 
     &  88, 65, 84, 64, 
     &  88, 25, 74, 65, 
     &  88, 25, 82, 74, 
     &  88, 83, 87, 71, 
     &  88, 71, 87, 84, 
     &  88, 82, 83, 74, 
     &  88, 74, 83, 71, 
     &  88, 65, 74, 71, 
     &  88, 71, 84, 65, 
     &  89, 86, 87, 84, 
     &  89, 39, 48, 44, 
     &  89, 44, 49, 43, 
     &  89, 44, 43, 36, 
     &  89, 44, 48, 50, 
     &  89, 48, 63, 50, 
     &  89, 86, 84, 54, 
     &  89, 51, 87, 86, 
     &  89, 44, 50, 49, 
     &  89, 50, 87, 49, 
     &  89, 43, 49, 51, 
     &  89, 49, 87, 51, 
     &  89, 39, 44, 36, 
     &  89, 36, 81, 39, 
     &  89, 63, 48, 47, 
     &  89, 47, 48, 40, 
     &  89, 46, 55, 47, 
     &  89, 38, 46, 47, 
     &  89, 55, 63, 47, 
     &  89, 55, 84, 63, 
     &  89, 43, 42, 34, 
     &  89, 43, 51, 42, 
     &  89, 45, 53, 54, 
     &  89, 53, 86, 54, 
     &  89, 45, 54, 46, 
     &  89, 42, 52, 41, 
     &  89, 41, 52, 53, 
     &  89, 52, 86, 53, 
     &  89, 42, 51, 52, 
     &  89, 51, 86, 52, 
     &  89, 46, 54, 55, 
     &  89, 54, 84, 55, 
     &  90, 56, 75, 61, 
     &  90, 24, 75, 56, 
     &  90, 27, 56, 61, 
     &  90, 61, 85, 27, 
     &  90, 75, 77, 61, 
     &  90, 80, 82, 57, 
     &  90, 85, 82, 80, 
     &  90, 57, 24, 56, 
     &  90, 72, 24, 57, 
     &  90, 57, 82, 72, 
     &  90, 80, 56, 27, 
     &  90, 85, 80, 27, 
     &  91, 85, 90, 77, 
     &  91, 86, 87, 69, 
     &  91, 78, 77, 69, 
     &  91, 83, 88, 82, 
     &  91, 90, 82, 88, 
     &  91, 67, 88, 86, 
     &  91, 88, 87, 86, 
     &  91, 87, 88, 83, 
     &  91, 83, 85, 78, 
     &  91, 78, 85, 77, 
     &  91, 77, 75, 68, 
     &  91, 77, 90, 75, 
     &  91, 69, 77, 68, 
     &  91, 68, 86, 69, 
     &  91, 68, 75, 67, 
     &  91, 67, 86, 68, 
     &  91, 24, 88, 67, 
     &  91, 90, 88, 24, 
     &  91, 69, 87, 70, 
     &  91, 87, 83, 70, 
     &  91, 75, 24, 67, 
     &  91, 75, 90, 24, 
     &  92, 89, 46, 45, 
     &  92, 41, 53, 45, 
     &  92, 89, 45, 53, 
     &  92, 89, 53, 41, 
     &  92, 89, 41, 42, 
     &  92, 35, 41, 45, 
     &  92, 33, 41, 35, 
     &  92, 35, 81, 33, 
     &  92, 35, 45, 37, 
     &  92, 81, 35, 37, 
     &  92, 34, 89, 42, 
     &  92, 81, 89, 34, 
     &  92, 33, 42, 41, 
     &  92, 37, 45, 46, 
     &  92, 37, 46, 38, 
     &  92, 81, 37, 38, 
     &  92, 33, 34, 42, 
     &  92, 33, 81, 34, 
     &  83, 74, 60, 71, 
     &  83, 60, 79, 71, 
     &  89, 39, 40, 48, 
     &  89, 39, 81, 40, 
     &  89, 36, 43, 34, 
     &  89, 34, 81, 36, 
     &  89, 63, 87, 50, 
     &  89, 84, 87, 63, 
     &  54, 88, 66, 86, 
     &  54, 88, 86, 84, 
     &  90, 72, 88, 24, 
     &  90, 82, 88, 72, 
     &  38, 47, 89, 40, 
     &  38, 89, 81, 40, 
     &  92, 46, 89, 38, 
     &  92, 89, 81, 38, 
     &  80, 23, 57, 56, 
     &  80, 57, 90, 56, 
     &  61, 85, 62, 77, 
     &  61, 90, 85, 77, 
     &  82, 85, 91, 83, 
     &  82, 90, 91, 85, 
     &  70, 91, 78, 83, 
     &  70, 78, 91, 69 /
      data triangle_label_save /
     &  4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 
     &  4, 4, 2, 2, 2, 2, 2, 2, 2, 2, 
     &  2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 
     &  3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 
     &  3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 
     &  3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 
     &  3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 
     &  3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 
     &  3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 
     &  3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 
     &  3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 
     &  3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 
     &  3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 
     &  3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 
     &  3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 
     &  3, 3, 3, 3 /
      data triangle_vertex_save /
     &  12, 20, 19, 
     &  12, 13, 20, 
     &  19, 20, 18, 
     &  20, 22, 18, 
     &  22, 17, 18, 
     &  22, 16, 17, 
     &  13, 21, 20, 
     &  13, 14, 21, 
     &  14, 15, 21, 
     &  22, 15, 16, 
     &  22, 21, 15, 
     &  20, 21, 22, 
     &  1, 9, 8, 
     &  2, 9, 1, 
     &  9, 7, 8, 
     &  2, 11, 9, 
     &  11, 2, 3, 
     &  11, 3, 4, 
     &  9, 10, 7, 
     &  7, 10, 6, 
     &  10, 5, 6, 
     &  11, 4, 5, 
     &  5, 10, 11, 
     &  9, 11, 10, 
     &  23, 1, 8, 
     &  26, 23, 8, 
     &  26, 8, 7, 
     &  27, 1, 23, 
     &  2, 1, 27, 
     &  26, 7, 28, 
     &  7, 6, 28, 
     &  27, 29, 2, 
     &  29, 3, 2, 
     &  29, 30, 3, 
     &  30, 4, 3, 
     &  6, 31, 28, 
     &  6, 5, 31, 
     &  5, 32, 31, 
     &  5, 4, 32, 
     &  12, 19, 33,  
     &  19, 34, 33, 
     &  19, 18, 34, 
     &  12, 33, 35, 
     &  12, 35, 13, 
     &  18, 36, 34, 
     &  36, 18, 17, 
     &  35, 37, 13, 
     &  13, 37, 14, 
     &  38, 14, 37, 
     &  38, 15, 14, 
     &  39, 36, 17, 
     &  39, 17, 16, 
     &  38, 40, 15, 
     &  40, 16, 15, 
     &  39, 16, 40, 
     &  33, 41, 35, 
     &  33, 42, 41, 
     &  33, 34, 42, 
     &  36, 43, 34, 
     &  43, 42, 34, 
     &  39, 44, 36, 
     &  44, 43, 36, 
     &  35, 45, 37, 
     &  35, 41, 45, 
     &  37, 46, 38, 
     &  37, 45, 46, 
     &  38, 47, 40, 
     &  38, 46, 47, 
     &  39, 48, 44, 
     &  39, 40, 48, 
     &  47, 48, 40, 
     &  44, 49, 43, 
     &  44, 50, 49, 
     &  44, 48, 50, 
     &  43, 51, 42, 
     &  43, 49, 51, 
     &  42, 52, 41, 
     &  42, 51, 52, 
     &  41, 53, 45, 
     &  41, 52, 53, 
     &  45, 54, 46, 
     &  45, 53, 54, 
     &  46, 55, 47, 
     &  46, 54, 55, 
     &  30, 32, 4, 
     &  23, 56, 27, 
     &  23, 57, 56, 
     &  23, 26, 57, 
     &  28, 58, 26, 
     &  58, 57, 26, 
     &  31, 59, 28, 
     &  59, 58, 28, 
     &  32, 60, 31, 
     &  60, 59, 31, 
     &  27, 61, 29, 
     &  27, 56, 61, 
     &  29, 62, 30, 
     &  29, 61, 62, 
     &  55, 63, 47, 
     &  63, 48, 47, 
     &  48, 63, 50, 
     &  54, 64, 55, 
     &  64, 65, 55, 
     &  65, 63, 55, 
     &  53, 66, 54, 
     &  66, 64, 54, 
     &  52, 67, 53, 
     &  67, 66, 53, 
     &  51, 68, 52, 
     &  68, 67, 52, 
     &  49, 69, 51, 
     &  69, 68, 51, 
     &  50, 70, 49, 
     &  70, 69, 49, 
     &  63, 71, 50, 
     &  71, 70, 50, 
     &  65, 71, 63, 
     &  64, 25, 65, 
     &  64, 73, 25, 
     &  64, 66, 73, 
     &  67, 24, 66, 
     &  24, 72, 66, 
     &  72, 73, 66, 
     &  68, 75, 67, 
     &  75, 24, 67, 
     &  69, 77, 68, 
     &  77, 75, 68, 
     &  70, 78, 69, 
     &  78, 77, 69, 
     &  62, 78, 30, 
     &  78, 76, 30, 
     &  76, 32, 30, 
     &  32, 76, 60, 
     &  61, 77, 62, 
     &  77, 78, 62, 
     &  56, 75, 61, 
     &  75, 77, 61, 
     &  57, 24, 56, 
     &  24, 75, 56, 
     &  58, 72, 57, 
     &  72, 24, 57, 
     &  59, 25, 58, 
     &  25, 73, 58, 
     &  73, 72, 58, 
     &  60, 74, 59, 
     &  74, 25, 59, 
     &  25, 74, 65, 
     &  65, 74, 71, 
     &  70, 76, 78, 
     &  71, 79, 70, 
     &  79, 76, 70, 
     &  79, 60, 76, 
     &  74, 60, 71, 
     &  60, 79, 71 /
      data vertex_coordinate_save /
     & 1.0,       0.2,         0.0, 
     & 1.0,       0.141421,    0.141421, 
     & 1.0,       0.0,         0.2, 
     & 1.0,      -0.141421,    0.141421, 
     & 1.0,      -0.2,         0.0, 
     & 1.0,      -0.141421,   -0.141421, 
     & 1.0,       0.0,        -0.2, 
     & 1.0,       0.141421,   -0.141421, 
     & 1.0,       0.066163,   -0.0302872, 
     & 1.0,      -0.0615154,  -0.0610739, 
     & 1.0,      -0.0306985,   0.0668017, 
     & 0.0,       0.2,         0.0, 
     & 0.0,       0.141421,   -0.141421, 
     & 0.0,       0.0,        -0.2, 
     & 0.0,      -0.141421,   -0.141421, 
     & 0.0,      -0.2,         0.0, 
     & 0.0,      -0.141421,    0.141421,  
     & 0.0,       0.0,         0.2,  
     & 0.0,       0.141421,    0.141421, 
     & 0.0,       0.0686748,   0.0255359, 
     & 0.0,       0.0,        -0.0865993, 
     & 0.0,      -0.0686749,   0.0255359, 
     & 0.8816,    0.185522,   -0.0747102, 
     & 0.642415,  0.187806,   -0.0687668, 
     & 0.627606, -0.0696445,  -0.187482, 
     & 0.876431,  0.0811908,  -0.182779, 
     & 0.881613,  0.186118,    0.0732131, 
     & 0.872048, -0.0699008,  -0.187387, 
     & 0.878318,  0.0844232,   0.181308, 
     & 0.845861, -0.0716063,   0.186742, 
     & 0.866503, -0.182493,   -0.0818307, 
     & 0.859402, -0.186751,    0.0715813, 
     & 0.131355,  0.18477,     0.0765501, 
     & 0.13317,   0.077694,    0.184292, 
     & 0.130862,  0.185301,   -0.0752567, 
     & 0.135181, -0.0749468,   0.185426, 
     & 0.130839,  0.0781729,  -0.18409, 
     & 0.131856, -0.0754694,  -0.185214, 
     & 0.135683, -0.184121,    0.0780993, 
     & 0.134207, -0.184959,   -0.0760928, 
     & 0.261923,  0.199982,    0.00264585, 
     & 0.263928,  0.144161,    0.138627, 
     & 0.268645,  0.00535339,  0.199928, 
     & 0.272346, -0.137646,    0.145098, 
     & 0.26108,   0.144683,   -0.138082, 
     & 0.260772,  0.00498797, -0.199938, 
     & 0.264253, -0.139152,   -0.143655, 
     & 0.270288, -0.199962,    0.00389323, 
     & 0.408181, -0.0730357,   0.186187, 
     & 0.411818, -0.184374,    0.0774991, 
     & 0.397539,  0.080738,    0.182979,  
     & 0.39192,   0.185619,    0.0744699, 
     & 0.392192,  0.184438,   -0.0773479, 
     & 0.389194,  0.0770141,  -0.184577, 
     & 0.38786,  -0.0747817,  -0.185493, 
     & 0.762413,  0.199986,   -0.0023425, 
     & 0.762987,  0.151152,   -0.13097, 
     & 0.741526,  0.0187858,  -0.199116, 
     & 0.746899, -0.128364,   -0.153371, 
     & 0.720076, -0.19917,    -0.0182053, 
     & 0.7628,    0.152219,    0.129728, 
     & 0.763882,  0.0434475,   0.195224, 
     & 0.399903, -0.1841,     -0.0781489, 
     & 0.506331, -0.00579066, -0.199916, 
     & 0.514514, -0.133894,   -0.148568, 
     & 0.526121,  0.135152,   -0.147424, 
     & 0.517967,  0.199953,   -0.0043215, 
     & 0.520585,  0.147847,    0.13469, 
     & 0.533956,  0.0124181,   0.199614, 
     & 0.558316, -0.136902,    0.145801, 
     & 0.549126, -0.199624,   -0.0122659, 
     & 0.657307,  0.117735,   -0.161674, 
     & 0.611189,  0.041829,   -0.195577, 
     & 0.631917, -0.164669,   -0.113508, 
     & 0.641444,  0.187001,    0.0709267, 
     & 0.720251, -0.155557,    0.125706, 
     & 0.647345,  0.0932963,   0.176906, 
     & 0.677484, -0.0430068,   0.195321, 
     & 0.635293, -0.188734,    0.0661777, 
     & 0.888023, -0.00868364, -0.00818647, 
     & 0.112146,  0.0,        -0.0118425, 
     & 0.676228,  0.0124197,  -0.0856487, 
     & 0.638436, -0.0639898,   0.0525795, 
     & 0.452586, -0.0410297,  -0.0704842, 
     & 0.762004, -0.0188614,   0.0693717, 
     & 0.463368,  0.0649048,   0.0262133, 
     & 0.473921, -0.0356443,   0.0388516, 
     & 0.557002,  0.0123705,  -0.0932599, 
     & 0.290986, -0.0200898,   0.00857934, 
     & 0.7038,    0.0856777,   0.0182744, 
     & 0.576134,  0.0436218,   0.0828782, 
     & 0.215187,  0.080855,   -0.0314946 /
      data vertex_label_save /
     &  3, 3, 3, 3, 3, 3, 3, 3, 2, 2, 
     &  2, 3, 3, 3, 3, 3, 3, 3, 3, 4, 
     &  4, 4, 3, 3, 3, 3, 3, 3, 3, 3, 
     &  3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 
     &  3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 
     &  3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 
     &  3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 
     &  3, 3, 3, 3, 3, 3, 3, 3, 3, 0, 
     &  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
     &  0, 0 /

      call r8mat_copy ( 3, vertices, vertex_coordinate_save, 
     &  vertex_coordinate )

      call i4vec_copy ( vertices, vertex_label_save, vertex_label )

      call i4mat_copy ( 3, triangles, triangle_vertex_save, 
     &  triangle_vertex )

      call i4vec_copy ( triangles, triangle_label_save, triangle_label )

      call i4mat_copy ( 4, tetrahedrons, tetrahedron_vertex_save, 
     &  tetrahedron_vertex )

      do i = 1, tetrahedrons
        tetrahedron_label(i) = 1
      end do

      return
      end
      subroutine cyl248_size ( dim, vertices, edges, triangles, 
     &  quadrilaterals, tetrahedrons, hexahedrons )

c*********************************************************************72
c
cc CYL248_SIZE defines the sizes for a 3D tetrahedral mesh.
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
c    Output, integer DIM, the spatial dimension, which should be 3.
c
c    Output, integer VERTICES, the number of vertices.
c
c    Output, integer EDGES, the number of edges (may be 0).
c
c    Output, integer TRIANGLES, the number of triangles (may be 0).
c
c    Output, integer QUADRILATERALS, the number of quadrilaterals (may be 0).
c
c    Output, integer TETRAHEDRONS, the number of tetrahedrons (may be 0).
c
c    Output, integer HEXAHEDRONS, the number of hexahedrons (may be 0).
c
      implicit none

      integer dim
      integer edges
      integer hexahedrons
      integer quadrilaterals
      integer tetrahedrons
      integer triangles
      integer vertices

      dim = 3
      vertices = 92
      edges = 0
      triangles = 154
      quadrilaterals = 0
      tetrahedrons = 248
      hexahedrons = 0

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
c    Input, integer DIM, the spatial dimension, which should be 3.
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
c    Input, double VERTEX_COORDINATE(3,VERTICES), the XYZ coordinates
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
      real ( kind = 8 ) vertex_coordinate(3,vertices)
      integer vertex_label(vertices)

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Vertices:'
      write ( *, '(a)' ) ' '
      do j = 1, vertices
        write ( *, '(3(2x,f10.4),2x,''('',i4,'')'')' ) 
     &    vertex_coordinate(1:3,j), vertex_label(j)
      end do

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
c    Input, integer DIM, the spatial dimension, which should be 3.
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
c    Output, double precision VERTEX_COORDINATE(3,VERTICES), the XYZ coordinates
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
      double precision vertex_coordinate(3,vertices)
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
      subroutine hexahexa_2x2x2_data ( dim, vertices, edges, triangles,
     &  quadrilaterals, tetrahedrons, hexahedrons, vertex_coordinate, 
     &  vertex_label, edge_vertex, edge_label, triangle_vertex, 
     &  triangle_label, quadrilateral_vertex, quadrilateral_label, 
     &  tetrahedron_vertex, tetrahedron_label, hexahedron_vertex, 
     &  hexahedron_label )

!*********************************************************************72
!
!! HEXAHEXA_2X2X2_DATA defines the data for a 3D hexahedral mesh.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    28 October 2010
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Pascal Frey,
!    MEDIT: An interactive mesh visualization software,
!    Technical Report RT-0253,
!    Institut National de Recherche en Informatique et en Automatique,
!    03 December 2001.
!
!  Parameters:
!
!    Input, integer DIM, the spatial dimension, which should be 3.
!
!    Input, integer VERTICES, the number of vertices.
!
!    Input, integer EDGES, the number of edges (may be 0).
!
!    Input, integer TRIANGLES, the number of triangles (may be 0).
!
!    Input, integer QUADRILATERALS, the number of quadrilaterals
!    (may be 0).
!
!    Input, integer TETRAHEDRONS, the number of tetrahedrons
!    (may be 0).
!
!    Input, integer HEXAHEDRONS, the number of hexahedrons
!    (may be 0).
!
!    Output, double precision VERTEX_COORDINATE(3,VERTICES), the XYZ 
!    coordinates of each vertex.
!
!    Output, integer VERTEX_LABEL(VERTICES), a label for
!    each vertex.
!
!    Output, integer EDGE_VERTEX(2,EDGES), the vertices that form
!    each edge.
!
!    Output, integer EDGE_LABEL(EDGES), a label for each edge.
!
!    Output, integer TRIANGLE_VERTEX(3,TRIANGLES), the vertices
!    that form each triangle.
!
!    Output, integer TRIANGLE_LABEL(TRIANGLES), a label for
!    each triangle.
!
!    Output, integer QUADRILATERAL_VERTEX(4,QUADRILATERALS), the
!    vertices that form each quadrilateral.
!
!    Output, integer QUADRILATERAL_LABEL(QUADRILATERALS), a label
!    for each quadrilateral.
!
!    Output, integer TETRAHEDRON_VERTEX(4,TETRAHEDRONS), the
!    vertices that form each tetrahedron.
!
!    Output, integer TETRAHEDRON_LABEL(TETRAHEDRONS), a label for
!    each tetrahedron.
!
!    Output, integer HEXAHEDRON_VERTEX(8,HEXAHEDRONS), the
!    vertices that form each hexahedron.
!
!    Output, integer HEXAHEDRON_LABEL(HEXAHEDRONS), a label for
!    each hexahedron.
!
      implicit none

      integer dim
      integer dim_save
      parameter ( dim_save = 3 )
      integer edges
      integer hexahedrons
      integer hexahedrons_save
      parameter ( hexahedrons_save = 8 )
      integer quadrilaterals
      integer quadrilaterals_save 
      parameter ( quadrilaterals_save = 24 )
      integer tetrahedrons
      integer triangles
      integer vertices
      integer vertices_save
      parameter ( vertices_save = 27 )

      integer edge_label(edges)
      integer edge_vertex(2,edges)
      integer hexahedron_label(hexahedrons)
      integer hexahedron_label_save(hexahedrons_save)
      integer hexahedron_vertex(8,hexahedrons)
      integer hexahedron_vertex_save(8,hexahedrons_save)
      integer quadrilateral_label(quadrilaterals)
      integer quadrilateral_label_save(quadrilaterals_save)
      integer quadrilateral_vertex(4,quadrilaterals)
      integer quadrilateral_vertex_save(4,quadrilaterals_save)
      integer tetrahedron_label(tetrahedrons)
      integer tetrahedron_vertex(4,tetrahedrons)
      integer triangle_label(triangles)
      integer triangle_vertex(3,triangles)
      double precision vertex_coordinate(dim,vertices)
      double precision vertex_coordinate_save(dim_save,vertices_save)
      integer vertex_label(vertices)
      integer vertex_label_save(vertices_save)

      save hexahedron_label_save
      save hexahedron_vertex_save
      save quadrilateral_label_save
      save quadrilateral_vertex_save
      save vertex_coordinate_save
      save vertex_label_save

      data hexahedron_label_save /
     &  1, 1, 1, 1, 1, 1, 1, 1 /

      data hexahedron_vertex_save /
     &  1,  2,  5,  4, 10, 11, 14, 13, 
     &  2,  3,  6,  5, 11, 12, 15, 14, 
     &  4,  5,  8,  7, 13, 14, 17, 16, 
     &  5,  6,  9,  8, 14, 15, 18, 17, 
     & 10, 11, 14, 13, 19, 20, 23, 22, 
     & 11, 12, 15, 14, 20, 21, 24, 23, 
     & 13, 14, 17, 16, 22, 23, 26, 25, 
     & 14, 15, 18, 17, 23, 24, 27, 26 /

      data quadrilateral_label_save /
     & 1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 
     & 3, 3, 4, 4, 4, 4, 5, 5, 5, 5, 
     & 6, 6, 6, 6 /

      data quadrilateral_vertex_save /
     &  1,  4,  5,  2, 
     &  2,  5,  6,  3, 
     &  4,  7,  8,  5, 
     &  5,  8,  9,  6, 
     &  1,  2, 11, 10, 
     &  2,  3, 12, 11, 
     & 10, 11, 20, 19, 
     & 11, 12, 21, 20, 
     &  3,  6, 15, 12, 
     &  6,  9, 18, 15, 
     & 12, 15, 24, 21, 
     & 15, 18, 27, 24, 
     &  7, 16, 17,  8, 
     &  8, 17, 18,  9, 
     & 16, 25, 26, 17, 
     & 17, 26, 27, 18, 
     &  1, 10, 13,  4, 
     &  4, 13, 16,  7, 
     & 10, 19, 22, 13, 
     & 13, 22, 25, 16, 
     & 19, 20, 23, 22, 
     & 20, 21, 24, 23, 
     & 22, 23, 26, 25, 
     & 23, 24, 27, 26 /

      data vertex_coordinate_save /
     &  0.0, 0.0, 0.0,
     &  0.5, 0.0, 0.0,
     &  1.0, 0.0, 0.0,
     &  0.0, 0.5, 0.0,
     &  0.5, 0.5, 0.0,
     &  1.0, 0.5, 0.0,
     &  0.0, 1.0, 0.0,
     &  0.5, 1.0, 0.0,
     &  1.0, 1.0, 0.0,
     &  0.0, 0.0, 0.5,
     &  0.5, 0.0, 0.5,
     &  1.0, 0.0, 0.5,
     &  0.0, 0.5, 0.5,
     &  0.5, 0.5, 0.5,
     &  1.0, 0.5, 0.5,
     &  0.0, 1.0, 0.5,
     &  0.5, 1.0, 0.5,
     &  1.0, 1.0, 0.5,
     &  0.0, 0.0, 1.0,
     &  0.5, 0.0, 1.0,
     &  1.0, 0.0, 1.0,
     &  0.0, 0.5, 1.0,
     &  0.5, 0.5, 1.0,
     &  1.0, 0.5, 1.0,
     &  0.0, 1.0, 1.0,
     &  0.5, 1.0, 1.0,
     &  1.0, 1.0, 1.0 /

      data vertex_label_save /
     &  5, 2, 3, 5, 1, 3, 5, 4, 4, 5,
     &  2, 3, 5, 0, 3, 5, 4, 4, 6, 6,
     &  6, 6, 6, 6, 6, 6, 6 /

      call r8vec_copy ( 3 * vertices, vertex_coordinate_save, 
     &  vertex_coordinate )
      call i4vec_copy ( vertices, vertex_label_save, vertex_label )
      call i4vec_copy ( 4 * quadrilaterals, quadrilateral_vertex_save,
     &  quadrilateral_vertex )
      call i4vec_copy ( quadrilaterals, quadrilateral_label_save,
     &  quadrilateral_label )
      call i4vec_copy ( 8 * hexahedrons, hexahedron_vertex_save, 
     &  hexahedron_vertex )
      call i4vec_copy ( hexahedrons, hexahedron_label_save, 
     &  hexahedron_label )

      return
      end
      subroutine hexahexa_2x2x2_size ( dim, vertices, edges, triangles, 
     &  quadrilaterals, tetrahedrons, hexahedrons )

c*********************************************************************72
c
cc HEXAHEXA_2x2x2_SIZE defines the sizes for a 3D hexahedral mesh.
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
c    Output, integer DIM, the spatial dimension, which should be 3.
c
c    Output, integer VERTICES, the number of vertices.
c
c    Output, integer EDGES, the number of edges (may be 0).
c
c    Output, integer TRIANGLES, the number of triangles (may be 0).
c
c    Output, integer QUADRILATERALS, the number of quadrilaterals (may be 0).
c
c    Output, integer TETRAHEDRONS, the number of tetrahedrons (may be 0).
c
c    Output, integer HEXAHEDRONS, the number of hexahedrons (may be 0).
c
      implicit none

      integer ( kind = 4 ) dim
      integer ( kind = 4 ) edges
      integer ( kind = 4 ) hexahedrons
      integer ( kind = 4 ) quadrilaterals
      integer ( kind = 4 ) tetrahedrons
      integer ( kind = 4 ) triangles
      integer ( kind = 4 ) vertices

      dim = 3
      vertices = 27
      edges = 0
      triangles = 0
      quadrilaterals = 24
      tetrahedrons = 0
      hexahedrons = 8

      return
      end
      subroutine i4mat_copy ( m, n, a1, a2 )

c*********************************************************************72
c
cc I4MAT_COPY copies an I4MAT.
c
c  Discussion:
c
c    An I4MAT is an array of I4's.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    04 June 2009
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer M, the number of rows.
c
c    Input, integer N, the number of columns.
c
c    Input, integer A1(M,N), the matrix to copy.
c
c    Output, integer A2(M,N), the copy.
c
      implicit none

      integer m
      integer n

      integer a1(m,n)
      integer a2(m,n)
      integer i
      integer j
 
      do j = 1, n
        do i = 1, m
          a2(i,j) = a1(i,j)
        end do
      end do

      return
      end
      subroutine i4vec_copy ( n, a1, a2 )

c*********************************************************************72
c
cc I4VEC_COPY copies an I4VEC.
c
c  Discussion:
c
c    An I4VEC is a vector of I4's.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    02 January 2007
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the length of the vectors.
c
c    Input, integer A1(N), the vector to be copied.
c
c    Output, integer A2(N), a copy of A1.
c
      implicit none

      integer n

      integer a1(n)
      integer a2(n)
      integer i

      do i = 1, n
        a2(i) = a1(i)
      end do

      return
      end
      subroutine ice_write ( filename, dim, vertices, edges, triangles,
     &  quadrilaterals, tetrahedrons, hexahedrons, vertex_coordinate, 
     &  vertex_label, edge_vertex, edge_label, triangle_vertex, 
     &  triangle_label, quadrilateral_vertex, quadrilateral_label, 
     &  tetrahedron_vertex, tetrahedron_label, hexahedron_vertex, 
     &  hexahedron_label )

c*********************************************************************72
c
cc ICE_WRITE writes 3D ICE sizes and data to a NETCDF file.
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
c    Input, character*(*) FILENAME, the name of the file to be created.
c    Ordinarily, the name should include the extension '.nc'.
c
c    Input, integer DIM, the spatial dimension, which should be 3.
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
c    Input, real VERTEX_COORDINATE(3,VERTICES), the XYZ coordinates
c    of each vertex.
c
c    Input, integer VERTEX_LABEL(VERTICES), a label for each vertex.
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
c    Input, integer QUADRILATERAL_VERTEX(4,QUADRILATERALS), the 
c    vertices that form each quadrilateral.
c
c    Input, integer QUADRILATERAL_LABEL(QUADRILATERALS), a label for 
c    each quadrilateral.
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

c
c  Disgruntled comment: 
c
c  In GFORTRAN for FORTRAN77, it seems that the
c  IMPLICIT must come BEFORE the "INCLUDE NETCDF" statement.
c
c  In GFORTRAN for FORTRAN90, it seems that the IMPLICIT must come 
c  AFTER the "USE NETCDF" statement.
c
c  Is it really necessary to be pointlessly picky about the
c  location of the IMPLICIT NONE statement?  Why can't I have seven
c  of them, scattered all through a routine, if I like?  When you
c  make a rule for no obvious benefit, all that's left over are
c  disadvantages and awkwardness.  And this becomes a serious matter
c  when you allow and promote the use of "INCLUDE" and "USE" statements
c  which often mean that you are including code you didn't write,
c  and can't easily examine.
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
      integer dimids(2)
      integer edge_label(edges)
      integer edge_vertex(2,edges) 
      character*(*) filename
      integer hexahedron_label(hexahedrons)
      integer hexahedron_vertex(8,hexahedrons)
      integer i
      integer j
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
      double precision vertex_coordinate(3,vertices)
      integer vertex_label(vertices)
      integer xtype
c
c  Create the file.  This automatically 'opens' it as well.
c
      mode = NF_CLOBBER
      status = nf_create ( filename, mode, ncid )
c
c  Put NETCDF into 'define' mode.
c
      status = nf_redef ( ncid )
c
c  Dimension information.
c
c  If a dimension has length 0, it seems to be taken to be the unlimited
c  dimension (not what you want) and then if you have two such dimensions,
c  you get a complaint that you have tried to define the unlimited dimension
c  twice.  The fix requires the programmer not to write anything whose
c  dimension is zero.
c
      status = nf_def_dim ( ncid, 'Dimension', dim, dim_dimension )

      status = nf_def_dim ( ncid, 'Vertices', vertices, dim_vertices )

      if ( 0 .lt. edges ) then
        status = nf_def_dim ( ncid, 'Edges', edges, dim_edges )
      end if

      if ( 0 .lt. triangles ) then
        status = nf_def_dim ( ncid, 'Triangles', triangles, 
     &    dim_triangles )
      end if

      if ( 0 .lt. quadrilaterals ) then
        status = nf_def_dim ( ncid, 'Quadrilaterals', quadrilaterals, 
     &    dim_quadrilaterals )
      end if

      if ( 0 .lt. tetrahedrons ) then
        status = nf_def_dim ( ncid, 'Tetrahedrons', tetrahedrons, 
     &    dim_tetrahedrons )
      end if

      if ( 0 .lt. hexahedrons ) then
        status = nf_def_dim ( ncid, 'Hexahedrons', hexahedrons, 
     &    dim_hexahedrons )
      end if

      status = nf_def_dim ( ncid, 'Two', 2, dim_two )
      status = nf_def_dim ( ncid, 'Three', 3, dim_three )
      status = nf_def_dim ( ncid, 'Four', 4, dim_four )
      status = nf_def_dim ( ncid, 'Eight', 8, dim_eight )
c
c  Define variables.
c
      ndims = 2
      dimids(1) = dim_three
      dimids(2) = dim_vertices
      status = nf_def_var ( ncid, 'Vertex_Coordinate', NF_DOUBLE, 
     &  ndims, dimids, var_vertex_coordinate )

      ndims = 1
      dimids(1) = dim_vertices
      status = nf_def_var ( ncid, 'Vertex_Label', NF_INT, ndims, 
     &  dimids, var_vertex_label )

      if ( 0 .lt. edges ) then
        ndims = 2
        dimids(1) = dim_two
        dimids(2) = dim_edges
        status = nf_def_var ( ncid, 'Edge_Vertex', NF_INT, ndims, 
     &    dimids, var_edge_vertex )

        ndims = 1
        dimids(1) = dim_edges
        status = nf_def_var ( ncid, 'Edge_Label', NF_INT, ndims, 
     &    dimids, var_edge_label )
      end if

      if ( 0 .lt. triangles ) then
        ndims = 2
        dimids(1) = dim_three
        dimids(2) = dim_triangles
        status = nf_def_var ( ncid, 'Triangle_Vertex', NF_INT, ndims, 
     &    dimids, var_triangle_vertex )

        ndims = 1
        dimids(1) = dim_triangles
        status = nf_def_var ( ncid, 'Triangle_Label', NF_INT, ndims, 
     &    dimids, var_triangle_label )
      end if

      if ( 0 .lt. quadrilaterals ) then
        ndims = 2
        dimids(1) = dim_four
        dimids(2) = dim_quadrilaterals
        status = nf_def_var ( ncid, 'Quadrilateral_Vertex', NF_INT, 
     &    ndims, dimids, var_quadrilateral_vertex )

        ndims = 1
        dimids(1) = dim_quadrilaterals
        status = nf_def_var ( ncid, 'Quadrilateral_Label', NF_INT, 
     &    ndims, dimids, var_quadrilateral_label )
      end if

      if ( 0 .lt. tetrahedrons ) then
        ndims = 2
        dimids(1) = dim_four
        dimids(2) = dim_tetrahedrons
        status = nf_def_var ( ncid, 'Tetrahedron_Vertex', NF_INT, 
     &    ndims, dimids, var_tetrahedron_vertex )

        ndims = 1
        dimids(1) = dim_tetrahedrons
        status = nf_def_var ( ncid, 'Tetrahedron_Label', NF_INT, 
     &    ndims, dimids, var_tetrahedron_label )
      end if

      if ( 0 .lt. hexahedrons ) then
        ndims = 2
        dimids(1) = dim_eight
        dimids(2) = dim_hexahedrons
        status = nf_def_var ( ncid, 'Hexahedron_Vertex', NF_INT, 
     &    ndims, dimids, var_hexahedron_vertex )

        ndims = 1
        dimids(1) = dim_hexahedrons
        status = nf_def_var ( ncid, 'Hexahedron_Label', NF_INT, 
     &    ndims, dimids, var_hexahedron_label )
      end if
c
c  Terminate the definition phase.
c
      status = nf_enddef ( ncid )
c
c  Write the data.
c
      status = nf_put_var_double ( ncid, var_vertex_coordinate, 
     &  vertex_coordinate )
      status = nf_put_var_int ( ncid, var_vertex_label, vertex_label )

      if ( 0 .lt. edges ) then
        status = nf_put_var_int ( ncid, var_edge_vertex, edge_vertex )
        status = nf_put_var_int ( ncid, var_edge_label, edge_label )
      end if

      if ( 0 .lt. triangles ) then
        status = nf_put_var_int ( ncid, var_triangle_vertex, 
     &    triangle_vertex )
        status = nf_put_var_int ( ncid, var_triangle_label, 
     &    triangle_label )
      end if

      if ( 0 .lt. quadrilaterals ) then
        status = nf_put_var_int ( ncid, var_quadrilateral_vertex, 
     &    quadrilateral_vertex )
        status = nf_put_var_int ( ncid, var_quadrilateral_label, 
     &    quadrilateral_label )
      end if

      if ( 0 .lt. tetrahedrons ) then
        status = nf_put_var_int ( ncid, var_tetrahedron_vertex, 
     &    tetrahedron_vertex )
        status = nf_put_var_int ( ncid, var_tetrahedron_label, 
     &    tetrahedron_label )
      end if

      if ( 0 .lt. hexahedrons ) then
        status = nf_put_var_int ( ncid, var_hexahedron_vertex, 
     &    hexahedron_vertex )
        status = nf_put_var_int ( ncid, var_hexahedron_label, 
     &    hexahedron_label )
      end if
c
c  Close the file.
c
      status = nf_close ( ncid )

      return
      end
      subroutine r8mat_copy ( m, n, a1, a2 )

c*********************************************************************72
c
cc R8MAT_COPY copies an R8MAT.
c
c  Discussion:
c
c    An R8MAT is an array of R8's.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    26 July 2008
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer M, N, the order of the matrix.
c
c    Input, double precision A1(M,N), the matrix to be copied.
c
c    Output, double precision A2(M,N), a copy of the matrix.
c
      implicit none

      integer m
      integer n

      double precision a1(m,n)
      double precision a2(m,n)
      integer i
      integer j

      do j = 1, n
        do i = 1, m
          a2(i,j) = a1(i,j)
        end do
      end do

      return
      end
      subroutine r8vec_copy ( n, a1, a2 )

c*********************************************************************72
c
cc R8VEC_COPY copies an R8VEC.
c
c  Discussion:
c
c    An R8VEC is a vector of R8 values.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    13 August 2008
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the length of the vectors.
c
c    Input, double precision A1(N), the vector to be copied.
c
c    Output, double precision A2(N), a copy of A1.
c
      implicit none

      integer n

      double precision a1(n)
      double precision a2(n)
      integer i

      do i = 1, n
        a2(i) = a1(i)
      end do

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
c    Input, integer DIM, the spatial dimension, which should be 3.
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
c    Output, integer DIM, the spatial dimension, which should be 3.
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
