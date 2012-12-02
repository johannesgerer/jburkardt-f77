      subroutine ch_cap ( ch )

c*********************************************************************72
c
cc CH_CAP capitalizes a single character.
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
c    Input/output, character CH, the character to capitalize.
c
      implicit none

      character ch
      integer itemp

      itemp = ichar ( ch )

      if ( 97 .le. itemp .and. itemp .le. 122 ) then
        ch = char ( itemp - 32 )
      end if

      return
      end
      function ch_eqi ( c1, c2 )

c*********************************************************************72
c
cc CH_EQI is a case insensitive comparison of two characters for equality.
c
c  Example:
c
c    CH_EQI ( 'A', 'a' ) is TRUE.
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
c    Input, character C1, C2, the characters to compare.
c
c    Output, logical CH_EQI, the result of the comparison.
c
      implicit none

      character c1
      character c1_cap
      character c2
      character c2_cap
      logical ch_eqi

      c1_cap = c1
      c2_cap = c2

      call ch_cap ( c1_cap )
      call ch_cap ( c2_cap )

      if ( c1_cap == c2_cap ) then
        ch_eqi = .true.
      else
        ch_eqi = .false.
      end if

      return
      end
      subroutine ch_to_digit ( c, digit )

c*********************************************************************72
c
cc CH_TO_DIGIT returns the integer value of a base 10 digit.
c
c  Example:
c
c     C   DIGIT
c    ---  -----
c    '0'    0
c    '1'    1
c    ...  ...
c    '9'    9
c    ' '    0
c    'X'   -1
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    04 August 1999
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, character C, the decimal digit, '0' through '9' or blank
c    are legal.
c
c    Output, integer DIGIT, the corresponding integer value.  If C was
c    'illegal', then DIGIT is -1.
c
      implicit none

      character c
      integer digit

      if ( lge ( c, '0' ) .and. lle ( c, '9' ) ) then

        digit = ichar ( c ) - 48

      else if ( c .eq. ' ' ) then

        digit = 0

      else

        digit = -1

      end if

      return
      end
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
c    Output, double precision VERTEX_COORDINATE(DIM,VERTICES), the coordinates
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

      call r8mat_copy ( dim, vertices, vertex_coordinate_save,
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
c    Output, integer DIM, the spatial dimension, which should be 2 or 3.
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
!    Input, integer DIM, the spatial dimension, which should be 2 or 3.
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
!    Output, double precision VERTEX_COORDINATE(DIM,VERTICES), the
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

      call r8vec_copy ( dim * vertices, vertex_coordinate_save,
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
c    Output, integer DIM, the spatial dimension, which should be 2 or 3.
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
      subroutine i4mat_zero ( m, n, a )

c*********************************************************************72
c
cc I4MAT_ZERO zeroes out an I4MAT.
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
c    22 November 2010
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer M, N, the row and column dimensions of the matrix.
c
c    Output, integer A(M,N), a matrix of zeroes.
c
      implicit none

      integer m
      integer n

      integer a(m,n)
      integer i
      integer j

      do j = 1, n
        do i = 1, m
          a(i,j) = 0
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
      subroutine i4vec_zero ( n, a )

c*********************************************************************72
c
cc I4VEC_ZERO sets the entries of an I4VEC to 0.
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
c    18 January 2007
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the number of entries in the vector.
c
c    Output, integer A(N), the vector, which has been set to zero.
c
      implicit none

      integer n

      integer a(n)
      integer i

      do i = 1, n
        a(i) = 0
      end do

      return
      end
      subroutine mesh_data_print ( dim, vertices, edges, triangles,
     &  quadrilaterals, tetrahedrons, hexahedrons, vertex_coordinate,
     &  vertex_label, edge_vertex, edge_label, triangle_vertex,
     &  triangle_label, quadrilateral_vertex, quadrilateral_label,
     &  tetrahedron_vertex, tetrahedron_label, hexahedron_vertex,
     &  hexahedron_label )

c*********************************************************************72
c
cc MESH_DATA_PRINT prints mesh data.
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
c    Input, double precision VERTEX_COORDINATE(DIM,VERTICES), the coordinates
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
      double precision vertex_coordinate(dim,vertices)
      integer vertex_label(vertices)

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Vertices:'
      write ( *, '(a)' ) ' '
      if ( dim == 2 ) then
        do j = 1, vertices
          write ( *, '(2(2x,f10.4),2x,''('',i4,'')'')' )
     &      vertex_coordinate(1:dim,j), vertex_label(j)
        end do
      else
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
      subroutine mesh_data_read ( filename, dim, vertices, edges,
     &  triangles, quadrilaterals, tetrahedrons, hexahedrons,
     &  vertex_coordinate, vertex_label, edge_vertex, edge_label,
     &  triangle_vertex, triangle_label, quadrilateral_vertex,
     &  quadrilateral_label, tetrahedron_vertex, tetrahedron_label,
     &  hexahedron_vertex, hexahedron_label )

c*********************************************************************72
c
cc MESH_DATA_READ reads data from a MESH file.
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
c  Parameters:
c
c    Input, character * ( * ) FILENAME, the name of the MESH file.
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
c    Output, double precision VERTEX_COORDINATE(DIM,VERTICES), the coordinates
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
c    Output, integer TRIANGLE_LABEL(TRIANGLES), a label for each
c    triangle.
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
c    Output, integer HEXAHEDRON_VERTEX(8,HEXAHEDRONS), the vertices
c    that form each hexahedron.
c
c    Output, integer HEXAHEDRON_LABEL(HEXAHEDRONS), a label for
c    each hexahedron.
c
      implicit none

      integer edges
      integer hexahedrons
      integer quadrilaterals
      integer tetrahedrons
      integer triangles
      integer vertices

      integer dim
      integer edge
      integer edge_label(edges)
      integer edge_vertex(2,edges)
      character * ( * ) filename
      integer fileunit
      integer hexahedron
      integer hexahedron_label(hexahedrons)
      integer hexahedron_vertex(8,hexahedrons)
      integer i
      integer i4vec(9)
      integer ierror
      integer ios
      integer j
      character * ( 80 ) keyword
      integer length
      integer line_num
      integer quadrilateral
      integer quadrilateral_label(quadrilaterals)
      integer quadrilateral_vertex(4,quadrilaterals)
      double precision r8vec(9)
      logical s_begin
      logical s_eqi
      integer s_len_trim
      integer tetrahedron
      integer tetrahedron_label(tetrahedrons)
      integer tetrahedron_vertex(4,tetrahedrons)
      character * ( 255 ) text
      integer triangle
      integer triangle_label(triangles)
      integer triangle_vertex(3,triangles)
      integer vertex
      double precision vertex_coordinate(dim,vertices)
      integer vertex_label(vertices)
c
c  Initialize everything to nothing.
c
      call r8mat_zero ( dim, vertices, vertex_coordinate )
      call i4vec_zero ( vertices, vertex_label )
      call i4mat_zero ( 2, edges, edge_vertex )
      call i4vec_zero ( edges, edge_label )
      call i4mat_zero ( 3, triangles, triangle_vertex )
      call i4vec_zero ( triangles, triangle_label )
      call i4mat_zero ( 4, quadrilaterals, quadrilateral_vertex )
      call i4vec_zero ( quadrilaterals, quadrilateral_label )
      call i4mat_zero ( 4, tetrahedrons, tetrahedron_vertex )
      call i4vec_zero ( tetrahedrons, tetrahedron_label )
      call i4mat_zero ( 8, hexahedrons, hexahedron_vertex )
      call i4vec_zero ( hexahedrons, hexahedron_label )
c
c  Open the file.
c
      call get_unit ( fileunit )

      open ( unit = fileunit, file = filename, status = 'old',
     &  iostat = ios )

      if ( ios .ne. 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'MESH_DATA_READ - Fatal errorc'
        write ( *, '(a)' ) '  Could not open file.'
        stop
      end if
c
c  Read lines til you get alphanumerics and determine a "mode"
c
      line_num = 0
      keyword = 'NONE'

10    continue

        read ( fileunit, '(a)', iostat = ios ) text

        if ( ios .ne. 0 ) then
          go to 20
        end if

        line_num = line_num + 1

        if ( s_len_trim ( text ) .eq. 0 ) then
          keyword = 'NONE'
          go to 10
        end if

        if ( text(1:1) .eq. '#' ) then
          go to 10
        end if
c
c  Remove initial blanks.
c
        call s_adjustl ( text )
c
c  Expecting a keyword.
c
            if ( s_eqi ( text, 'CORNERS' ) ) then

          keyword = 'SKIP'

        else if ( s_eqi ( text, 'DIMENSION' ) ) then

          keyword = 'DIMENSION'

        else if ( s_eqi ( text, 'EDGES' ) ) then

          keyword = 'EDGES'

        else if ( s_eqi ( text, 'END' ) ) then

          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) '  END statement encountered.'
          go to 20

        else if ( s_eqi ( text, 'HEXAHEDRA' ) .or.
     &            s_eqi ( text, 'HEXAHEDRONS' ) ) then

          keyword = 'HEXAHEDRONS'

        else if ( s_begin ( text, 'MESHVERSIONFORMATTED' ) ) then

        else if ( s_eqi ( text, 'NORMALATQUADRILATERALVERTICES' ) ) then

          keyword = 'SKIP'

        else if ( s_eqi ( text, 'NORMALATTRIANGLEVERTICES' ) ) then

          keyword = 'SKIP'

        else if ( s_eqi ( text, 'NORMALATVERTICES' ) ) then

          keyword = 'SKIP'

        else if ( s_eqi ( text, 'NORMALS' ) ) then

          keyword = 'SKIP'

        else if ( s_eqi ( text, 'QUADRILATERALS' ) ) then

          keyword = 'QUADRILATERALS'

        else if ( s_eqi ( text, 'REQUIREDEDGES' ) ) then

          keyword = 'SKIP'

        else if ( s_eqi ( text, 'REQUIREDVERTICES' ) ) then

          keyword = 'SKIP'

        else if ( s_eqi ( text, 'RIDGES' ) ) then

          keyword = 'SKIP'

        else if ( s_eqi ( text, 'TANGENTATEDGES' ) ) then

          keyword = 'SKIP'

        else if ( s_eqi ( text, 'TANGENTS' ) ) then

          keyword = 'SKIP'

        else if ( s_eqi ( text, 'TETRAHEDRA' ) .or.
     &            s_eqi ( text, 'TETRAHEDRONS' ) ) then

          keyword = 'TETRAHEDRONS'

        else if ( s_eqi ( text, 'TRIANGLES' ) ) then

          keyword = 'TRIANGLES'

        else if ( s_eqi ( text, 'VERTICES' ) ) then

          keyword = 'VERTICES'
c
c  Presumably, numeric data to be processed by keyword.
c
        else if ( s_eqi ( keyword, 'DIMENSION' ) ) then

          call s_to_i4 ( text, dim, ierror, length )

          keyword = 'NONE'

        else if ( s_eqi ( keyword, 'EDGES' ) ) then

          call s_to_i4 ( text, edges, ierror, length )

          keyword = 'EDGE_VERTEX'
          edge = 0

        else if ( s_eqi ( keyword, 'EDGE_VERTEX' ) ) then

          call s_to_i4vec ( text, 3, i4vec, ierror )
          edge = edge + 1
          do i = 1, 2
            edge_vertex(i,edge) = i4vec(i)
          end do
          edge_label(edge) = i4vec(3)

        else if ( s_eqi ( keyword, 'HEXAHEDRONS' ) ) then

          call s_to_i4 ( text, hexahedrons, ierror, length )

          keyword = 'HEXAHEDRON_VERTEX'
          hexahedron = 0

        else if ( s_eqi ( keyword, 'HEXAHEDRON_VERTEX' ) ) then

          call s_to_i4vec ( text, 9, i4vec, ierror )
          hexahedron = hexahedron + 1
          do i = 1, 8
            hexahedron_vertex(i,hexahedron) = i4vec(i)
          end do
          hexahedron_label(hexahedron) = i4vec(9)

        else if ( s_eqi ( keyword, 'QUADRILATERALS' ) ) then

          call s_to_i4 ( text, quadrilaterals, ierror, length )

          keyword = 'QUADRILATERAL_VERTEX'
          quadrilateral = 0

        else if ( s_eqi ( keyword, 'QUADRILATERAL_VERTEX' ) ) then

          call s_to_i4vec ( text, 5, i4vec, ierror )
          quadrilateral = quadrilateral + 1
          do i = 1, 4
            quadrilateral_vertex(i,quadrilateral) = i4vec(i)
          end do
          quadrilateral_label(quadrilateral) = i4vec(5)

        else if ( s_eqi ( keyword, 'TETRAHEDRONS' ) ) then

          call s_to_i4 ( text, tetrahedrons, ierror, length )

          keyword = 'TETRAHEDRON_VERTEX'
          tetrahedron = 0

        else if ( s_eqi ( keyword, 'TETRAHEDRON_VERTEX' ) ) then

          call s_to_i4vec ( text, 5, i4vec, ierror )
          tetrahedron = tetrahedron + 1
          do i = 1, 4
            tetrahedron_vertex(i,tetrahedron) = i4vec(i)
          end do
          tetrahedron_label(tetrahedron) = i4vec(5)

        else if ( s_eqi ( keyword, 'TRIANGLES' ) ) then

          call s_to_i4 ( text, triangles, ierror, length )

          keyword = 'TRIANGLE_VERTEX'
          triangle = 0

        else if ( s_eqi ( keyword, 'TRIANGLE_VERTEX' ) ) then

          call s_to_i4vec ( text, 4, i4vec, ierror )
          triangle = triangle + 1
          do i = 1, 3
            triangle_vertex(i,triangle) = i4vec(i)
          end do
          triangle_label(triangle) = i4vec(4)

        else if ( s_eqi ( keyword, 'VERTICES' ) ) then

          call s_to_i4 ( text, vertices, ierror, length )

          keyword = 'VERTEX_COORDINATE'
          vertex = 0

        else if ( s_eqi ( keyword, 'VERTEX_COORDINATE' ) ) then

          call s_to_r8vec ( text, dim + 1, r8vec, ierror )
          vertex = vertex + 1
          do i = 1, dim
            vertex_coordinate(i,vertex) = r8vec(i)
          end do
          vertex_label(vertex) = int ( r8vec(dim + 1) )

        else if ( s_eqi ( keyword, 'SKIP' ) ) then

        else

          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) 'MESH_DATA_READ - Fatal error!'
          write ( *, '(a,i8)' )
     &      '  Could not find keyword while reading line ', line_num
          write ( *, '(a)' ) '"' // trim ( text ) // '".'
          stop

        end if

      go to 10

20    continue
c
c  Close the file.
c
      close ( unit = fileunit )

      write ( *, '(a)' ) ' '
      write ( *, '(a,i8,a)' ) '  Read ', line_num,
     &  ' lines from "' // trim ( filename ) // '".'

      return
      end
      subroutine mesh_size_print ( dim, vertices, edges, triangles,
     &  quadrilaterals, tetrahedrons, hexahedrons )

c*********************************************************************72
c
cc MESH_SIZE_PRINT prints mesh sizes.
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
      subroutine mesh_size_read ( filename, dim, vertices, edges,
     &  triangles, quadrilaterals, tetrahedrons, hexahedrons )

c*********************************************************************72
c
cc MESH_SIZE_READ reads sizes from a MESH file.
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
c  Parameters:
c
c    Input, character * ( * ) FILENAME, the name of the MESH file.
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

      integer dim
      integer edges
      character * ( * ) filename
      integer fileunit
      integer hexahedrons
      integer ierror
      integer ios
      character * ( 80 ) keyword
      integer length
      integer line_num
      integer quadrilaterals
      logical s_begin
      logical s_eqi
      integer s_len_trim
      integer tetrahedrons
      character * ( 255 ) text
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
      call get_unit ( fileunit )

      open ( unit = fileunit, file = filename, status = 'old',
     &  iostat = ios )

      if ( ios .ne. 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'MESH_SIZE_READ - Fatal error!'
        write ( *, '(a)' ) '  Could not open file.'
        stop
      end if
c
c  Read lines til you get alphanumerics and determine a "mode"
c
      line_num = 0
      keyword = 'NONE'

10    continue

        read ( fileunit, '(a)', iostat = ios ) text

        if ( ios .ne. 0 ) then
          go to 20
        end if

        line_num = line_num + 1

        if ( s_len_trim ( text ) .eq. 0 ) then
          keyword = 'NONE'
          go to 10
        end if

        if ( text(1:1) .eq. '#' ) then
          go to 10
        end if
c
c  Remove initial blanks.
c
        call s_adjustl ( text )
c
c  Expecting a keyword.
c
            if ( s_eqi ( text, 'CORNERS' ) ) then

          keyword = 'SKIP'

        else if ( s_eqi ( text, 'DIMENSION' ) ) then

          keyword = 'DIMENSION'

        else if ( s_eqi ( text, 'EDGES' ) ) then

          keyword = 'EDGES'

        else if ( s_eqi ( text, 'END' ) ) then

          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) '  END statement encountered.'
          go to 20

        else if ( s_eqi ( text, 'HEXAHEDRA' ) .or.
     &            s_eqi ( text, 'HEXAHEDRONS' ) ) then

          keyword = 'HEXAHEDRONS'

        else if ( s_begin ( text, 'MESHVERSIONFORMATTED' ) ) then

        else if ( s_eqi ( text, 'NORMALATQUADRILATERALVERTICES' ) ) then

          keyword = 'SKIP'

        else if ( s_eqi ( text, 'NORMALATTRIANGLEVERTICES' ) ) then

          keyword = 'SKIP'

        else if ( s_eqi ( text, 'NORMALATVERTICES' ) ) then

          keyword = 'SKIP'

        else if ( s_eqi ( text, 'NORMALS' ) ) then

          keyword = 'SKIP'

        else if ( s_eqi ( text, 'QUADRILATERALS' ) ) then

          keyword = 'QUADRILATERALS'

        else if ( s_eqi ( text, 'REQUIREDEDGES' ) ) then

          keyword = 'SKIP'

        else if ( s_eqi ( text, 'REQUIREDVERTICES' ) ) then

          keyword = 'SKIP'

        else if ( s_eqi ( text, 'RIDGES' ) ) then

          keyword = 'SKIP'

        else if ( s_eqi ( text, 'TANGENTATEDGES' ) ) then

          keyword = 'SKIP'

        else if ( s_eqi ( text, 'TANGENTS' ) ) then

          keyword = 'SKIP'

        else if ( s_eqi ( text, 'TETRAHEDRA' ) .or.
     &            s_eqi ( text, 'TETRAHEDRONS' ) ) then

          keyword = 'TETRAHEDRONS'

        else if ( s_eqi ( text, 'TRIANGLES' ) ) then

          keyword = 'TRIANGLES'

        else if ( s_eqi ( text, 'VERTICES' ) ) then

          keyword = 'VERTICES'
c
c  Presumably, numeric data to be processed by keyword.
c
        else if ( s_eqi ( keyword, 'DIMENSION' ) ) then

          call s_to_i4 ( text, dim, ierror, length )

          keyword = 'NONE'

        else if ( s_eqi ( keyword, 'EDGES' ) ) then

          call s_to_i4 ( text, edges, ierror, length )

          keyword = 'EDGE_VERTEX'

        else if ( s_eqi ( keyword, 'EDGE_VERTEX' ) ) then

        else if ( s_eqi ( keyword, 'HEXAHEDRONS' ) ) then

          call s_to_i4 ( text, hexahedrons, ierror, length )

          keyword = 'HEXAHEDRON_VERTEX'

        else if ( s_eqi ( keyword, 'HEXAHEDRON_VERTEX' ) ) then

        else if ( s_eqi ( keyword, 'QUADRILATERALS' ) ) then

          call s_to_i4 ( text, quadrilaterals, ierror, length )

          keyword = 'QUADRILATERAL_VERTEX'

        else if ( s_eqi ( keyword, 'QUADRILATERAL_VERTEX' ) ) then

        else if ( s_eqi ( keyword, 'TETRAHEDRONS' ) ) then

          call s_to_i4 ( text, tetrahedrons, ierror, length )

          keyword = 'TETRAHEDRON_VERTEX'

        else if ( s_eqi ( keyword, 'TETRAHEDRON_VERTEX' ) ) then

        else if ( s_eqi ( keyword, 'TRIANGLES' ) ) then

          call s_to_i4 ( text, triangles, ierror, length )

          keyword = 'TRIANGLE_VERTEX'

        else if ( s_eqi ( keyword, 'TRIANGLE_VERTEX' ) ) then

        else if ( s_eqi ( keyword, 'VERTICES' ) ) then

          call s_to_i4 ( text, vertices, ierror, length )

          keyword = 'VERTEX_COORDINATE'

        else if ( s_eqi ( keyword, 'VERTEX_COORDINATE' ) ) then

        else if ( s_eqi ( keyword, 'SKIP' ) ) then

        else

          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) 'MESH_SIZE_READ - Fatal error!'
          write ( *, '(a,i8)' )
     &      '  Could not find keyword while reading line ', line_num
          write ( *, '(a)' ) '"' // trim ( text ) // '".'
          stop

        end if

      go to 10

20    continue
c
c  Close the file.
c
      close ( unit = fileunit )

      write ( *, '(a)' ) ' '
      write ( *, '(a,i8,a)' ) '  Read ', line_num,
     &  ' lines from "' // trim ( filename ) // '".'

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
c    Input, double precision VERTEX_COORRDINATE(DIM,VERTICES), the coordinates
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
      subroutine r8mat_zero ( m, n, a )

c*********************************************************************72
c
cc R8MAT_ZERO zeroes an R8MAT.
c
c  Discussion:
c
c    An R8MAT is an array of R8 values.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    18 July 2010
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer M, N, the number of rows and columns.
c
c    Output, double precision A(M,N), the matrix of zeroes.
c
      implicit none

      integer m
      integer n

      double precision a(m,n)
      integer i
      integer j

      do j = 1, n
        do i = 1, m
          a(i,j) = 0.0D+00
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
      subroutine s_adjustl ( s )

c*********************************************************************72
c
cc S_ADJUSTL flushes a string left.
c
c  Discussion:
c
c    Both blanks and tabs are treated as "white space".
c
c    This routine is similar to the FORTRAN90 ADJUSTL routine.
c
c  Example:
c
c    Input             Output
c
c    '     Hello'      'Hello     '
c    ' Hi therec  '    'Hi therec   '
c    'Fred  '          'Fred  '
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    03 Jun3 2008
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input/output, character * ( * ) S.
c    On input, S is a string of characters.
c    On output, any initial blank or tab characters have been cut.
c
      implicit none

      integer i
      integer nonb
      character * ( * ) s
      integer s_length
      character tab

      tab = char ( 9 )
c
c  Check the length of the string to the last nonblank.
c  If nonpositive, return.
c
      s_length = len_trim ( s )

      if ( s_length .le. 0 ) then
        return
      end if
c
c  Find NONB, the location of the first nonblank, nontab.
c
      nonb = 0

      do i = 1, s_length

        if ( s(i:i) .ne. ' ' .and. s(i:i) .ne. tab ) then
          nonb = i
          go to 10
        end if

      end do

10    continue

      if ( nonb .eq. 0 ) then
        s = ' '
        return
      end if
c
c  Shift the string left.
c
      if ( 1 .lt. nonb ) then
        do i = 1, s_length + 1 - nonb
          s(i:i) = s(i+nonb-1:i+nonb-1)
        end do
      end if
c
c  Blank out the end of the string.
c
      s(s_length+2-nonb:s_length) = ' '

      return
      end
      function s_begin ( s1, s2 )

c*********************************************************************72
c
cc S_BEGIN is TRUE if one string matches the beginning of the other.
c
c  Discussion:
c
c    The strings are compared, ignoring blanks, spaces and capitalization.
c
c  Example:
c
c     S1              S2      S_BEGIN
c
c    'Bob'          'BOB'     TRUE
c    '  B  o b '    ' bo b'   TRUE
c    'Bob'          'Bobby'   TRUE
c    'Bobo'         'Bobb'    FALSE
c    ' '            'Bob'     FALSE    (Do not allow a blank to match
c                                       anything but another blank string.)
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
c    Input, character * ( ) S1, S2, the strings to be compared.
c
c    Output, logical S_BEGIN, is TRUE if the strings match up to
c    the end of the shorter string, ignoring case.
c
      implicit none

      logical ch_eqi
      integer i1
      integer i2
      logical s_begin
      character * ( * )  s1
      integer s1_length
      character * ( * )  s2
      integer s2_length

      s1_length = len_trim ( s1 )
      s2_length = len_trim ( s2 )
c
c  If either string is blank, then both must be blank to match.
c  Otherwise, a blank string matches anything, which is not
c  what most people want.
c
      if ( s1_length .eq. 0 .or. s2_length .eq. 0 ) then

        if ( s1_length .eq. 0 .and. s2_length .eq. 0 ) then
          s_begin = .true.
        else
          s_begin = .false.
        end if

        return

      end if

      i1 = 0
      i2 = 0
c
c  Find the next nonblank in S1.
c
10    continue

20      continue

          i1 = i1 + 1

          if ( s1_length .lt. i1 ) then
            s_begin = .true.
            return
          end if

          if ( s1(i1:i1) .ne. ' ' ) then
            go to 30
          end if

        go to 20

30      continue
c
c  Find the next nonblank in S2.
c
40      continue

          i2 = i2 + 1

          if ( s2_length .lt. i2 ) then
            s_begin = .true.
            return
          end if

          if ( s2(i2:i2) .ne. ' ' ) then
            go to 50
          end if

        go to 40

50      continue
c
c  If the characters match, get the next pair.
c
        if ( .not. ch_eqi ( s1(i1:i1), s2(i2:i2) ) ) then
          go to 60
        end if

      go to 10

60    continue

      s_begin = .false.

      return
      end
      function s_eqi ( s1, s2 )

c*********************************************************************72
c
cc S_EQI is a case insensitive comparison of two strings for equality.
c
c  Example:
c
c    S_EQI ( 'Anjana', 'ANJANA' ) is TRUE.
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
c    Input, character*(*) S1, S2, the strings to compare.
c
c    Output, logical S_EQI, the result of the comparison.
c
      implicit none

      character c1
      character c2
      integer i
      integer lenc
      logical s_eqi
      character*(*) s1
      integer s1_length
      character*(*) s2
      integer s2_length

      s1_length = len ( s1 )
      s2_length = len ( s2 )
      lenc = min ( s1_length, s2_length )

      s_eqi = .false.

      do i = 1, lenc

        c1 = s1(i:i)
        c2 = s2(i:i)
        call ch_cap ( c1 )
        call ch_cap ( c2 )

        if ( c1 .ne. c2 ) then
          return
        end if

      end do

      do i = lenc + 1, s1_length
        if ( s1(i:i) .ne. ' ' ) then
          return
        end if
      end do

      do i = lenc + 1, s2_length
        if ( s2(i:i) .ne. ' ' ) then
          return
        end if
      end do

      s_eqi = .true.

      return
      end
      function s_len_trim ( s )

c*********************************************************************72
c
cc S_LEN_TRIM returns the length of a string to the last nonblank.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    05 March 2004
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, character*(*) S, a string.
c
c    Output, integer S_LEN_TRIM, the length of the string to the last nonblank.
c
      implicit none

      integer i
      character*(*) s
      integer s_len_trim

      do i = len ( s ), 1, -1

        if ( s(i:i) .ne. ' ' ) then
          s_len_trim = i
          return
        end if

      end do

      s_len_trim = 0

      return
      end
      subroutine s_to_i4 ( s, ival, ierror, length )

c*********************************************************************72
c
cc S_TO_I4 reads an I4 from a string.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    28 April 2008
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, character * ( * ) S, a string to be examined.
c
c    Output, integer IVAL, the integer value read from the string.
c    If the string is blank, then IVAL will be returned 0.
c
c    Output, integer IERROR, an error flag.
c    0, no error.
c    1, an error occurred.
c
c    Output, integer LENGTH, the number of characters of S
c    used to make IVAL.
c
      implicit none

      character c
      integer i
      integer ierror
      integer isgn
      integer istate
      integer ival
      integer length
      character * ( * ) s
      integer s_len_trim

      ierror = 0
      istate = 0
      isgn = 1
      ival = 0

      do i = 1, s_len_trim ( s )

        c = s(i:i)
c
c  Haven't read anything.
c
        if ( istate .eq. 0 ) then

          if ( c .eq. ' ' ) then

          else if ( c .eq. '-' ) then
            istate = 1
            isgn = -1
          else if ( c .eq. '+' ) then
            istate = 1
            isgn = + 1
          else if ( lle ( '0', c ) .and. lle ( c, '9' ) ) then
            istate = 2
            ival = ichar ( c ) - ichar ( '0' )
          else
            ierror = 1
            return
          end if
c
c  Have read the sign, expecting digits.
c
        else if ( istate .eq. 1 ) then

          if ( c .eq. ' ' ) then

          else if ( lle ( '0', c ) .and. lle ( c, '9' ) ) then
            istate = 2
            ival = ichar ( c ) - ichar ( '0' )
          else
            ierror = 1
            return
          end if
c
c  Have read at least one digit, expecting more.
c
        else if ( istate .eq. 2 ) then

          if ( lle ( '0', c ) .and. lle ( c, '9' ) ) then
            ival = 10 * ival + ichar ( c ) - ichar ( '0' )
          else
            ival = isgn * ival
            length = i - 1
            return
          end if

        end if

      end do
c
c  If we read all the characters in the string, see if we're OK.
c
      if ( istate .eq. 2 ) then
        ival = isgn * ival
        length = s_len_trim ( s )
      else
        ierror = 1
        length = 0
      end if

      return
      end
      subroutine s_to_i4vec ( s, n, ivec, ierror )

c*********************************************************************72
c
cc S_TO_I4VEC reads an I4VEC from a string.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    28 April 2008
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, character * ( * ) S, the string to be read.
c
c    Input, integer N, the number of values expected.
c
c    Output, integer IVEC(N), the values read from the string.
c
c    Output, integer IERROR, error flag.
c    0, no errors occurred.
c    -K, could not read data for entries -K through N.
c
      implicit none

      integer n

      integer i
      integer ierror
      integer ilo
      integer ivec(n)
      integer length
      character * ( * ) s

      i = 0
      ierror = 0
      ilo = 1

10    continue

      if ( i .lt. n ) then

        i = i + 1

        call s_to_i4 ( s(ilo:), ivec(i), ierror, length )

        if ( ierror .ne. 0 ) then
          ierror = -i
          go to 20
        end if

        ilo = ilo + length

      go to 10

      end if

20    continue

      return
      end
      subroutine s_to_r8 ( s, dval, ierror, length )

c*********************************************************************72
c
cc S_TO_R8 reads an R8 from a string.
c
c  Discussion:
c
c    The routine will read as many characters as possible until it reaches
c    the end of the string, or encounters a character which cannot be
c    part of the number.
c
c    Legal input is:
c
c       1 blanks,
c       2 '+' or '-' sign,
c       2.5 blanks
c       3 integer part,
c       4 decimal point,
c       5 fraction part,
c       6 'E' or 'e' or 'D' or 'd', exponent marker,
c       7 exponent sign,
c       8 exponent integer part,
c       9 exponent decimal point,
c      10 exponent fraction part,
c      11 blanks,
c      12 final comma or semicolon,
c
c    with most quantities optional.
c
c  Example:
c
c    S                 DVAL
c
c    '1'               1.0
c    '     1   '       1.0
c    '1A'              1.0
c    '12,34,56'        12.0
c    '  34 7'          34.0
c    '-1E2ABCD'        -100.0
c    '-1X2ABCD'        -1.0
c    ' 2E-1'           0.2
c    '23.45'           23.45
c    '-4.2E+2'         -420.0
c    '17d2'            1700.0
c    '-14e-2'         -0.14
c    'e2'              100.0
c    '-12.73e-9.23'   -12.73 * 10.0**(-9.23)
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    28 April 2008
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, character * ( * ) S, the string containing the
c    data to be read.  Reading will begin at position 1 and
c    terminate at the end of the string, or when no more
c    characters can be read to form a legal real.  Blanks,
c    commas, or other nonnumeric data will, in particular,
c    cause the conversion to halt.
c
c    Output, double precision DVAL, the value read from the string.
c
c    Output, integer IERROR, error flag.
c    0, no errors occurred.
c    1, 2, 6 or 7, the input number was garbled.  The
c    value of IERROR is the last type of input successfully
c    read.  For instance, 1 means initial blanks, 2 means
c    a plus or minus sign, and so on.
c
c    Output, integer LENGTH, the number of characters read
c    to form the number, including any terminating
c    characters such as a trailing comma or blanks.
c
      implicit none

      logical ch_eqi
      character c
      double precision dval
      integer ierror
      integer ihave
      integer isgn
      integer iterm
      integer jbot
      integer jsgn
      integer jtop
      integer length
      integer nchar
      integer ndig
      double precision rbot
      double precision rexp
      double precision rtop
      character * ( * ) s
      integer s_len_trim

      nchar = s_len_trim ( s )

      ierror = 0
      dval = 0.0D+00
      length = -1
      isgn = 1
      rtop = 0
      rbot = 1
      jsgn = 1
      jtop = 0
      jbot = 1
      ihave = 1
      iterm = 0

10    continue

        length = length + 1

        if ( nchar .lt. length+1 ) then
          go to 20
        end if

        c = s(length+1:length+1)
c
c  Blank character.
c
        if ( c .eq. ' ' ) then

          if ( ihave .eq. 2 ) then

          else if ( ihave .eq. 6 .or. ihave .eq. 7 ) then
            iterm = 1
          else if ( 1 .lt. ihave ) then
            ihave = 11
          end if
c
c  Comma.
c
        else if ( c .eq. ',' .or. c .eq. ';' ) then

          if ( ihave .ne. 1 ) then
            iterm = 1
            ihave = 12
            length = length + 1
          end if
c
c  Minus sign.
c
        else if ( c .eq. '-' ) then

          if ( ihave .eq. 1 ) then
            ihave = 2
            isgn = -1
          else if ( ihave .eq. 6 ) then
            ihave = 7
            jsgn = -1
          else
            iterm = 1
          end if
c
c  Plus sign.
c
        else if ( c .eq. '+' ) then

          if ( ihave .eq. 1 ) then
            ihave = 2
          else if ( ihave .eq. 6 ) then
            ihave = 7
          else
            iterm = 1
          end if
c
c  Decimal point.
c
        else if ( c .eq. '.' ) then

          if ( ihave .lt. 4 ) then
            ihave = 4
          else if ( 6 .le. ihave .and. ihave .le. 8 ) then
            ihave = 9
          else
            iterm = 1
          end if
c
c  Scientific notation exponent marker.
c
        else if ( ch_eqi ( c, 'E' ) .or. ch_eqi ( c, 'D' ) ) then

          if ( ihave .lt. 6 ) then
            ihave = 6
          else
            iterm = 1
          end if
c
c  Digit.
c
        else if ( ihave .lt. 11 .and. lle ( '0', c )
     &    .and. lle ( c, '9' ) ) then

          if ( ihave .le. 2 ) then
            ihave = 3
          else if ( ihave .eq. 4 ) then
            ihave = 5
          else if ( ihave .eq. 6 .or. ihave .eq. 7 ) then
            ihave = 8
          else if ( ihave .eq. 9 ) then
            ihave = 10
          end if

          call ch_to_digit ( c, ndig )

          if ( ihave .eq. 3 ) then
            rtop = 10.0D+00 * rtop + dble ( ndig )
          else if ( ihave .eq. 5 ) then
            rtop = 10.0D+00 * rtop + dble ( ndig )
            rbot = 10.0D+00 * rbot
          else if ( ihave .eq. 8 ) then
            jtop = 10 * jtop + ndig
          else if ( ihave .eq. 10 ) then
            jtop = 10 * jtop + ndig
            jbot = 10 * jbot
          end if
c
c  Anything else is regarded as a terminator.
c
        else
          iterm = 1
        end if
c
c  If we haven't seen a terminator, and we haven't examined the
c  entire string, go get the next character.
c
        if ( iterm .eq. 1 ) then
          go to 20
        end if

        go to 10

20    continue
c
c  If we haven't seen a terminator, and we have examined the
c  entire string, then we're done, and LENGTH is equal to NCHAR.
c
      if ( iterm .ne. 1 .and. length+1 .eq. nchar ) then
        length = nchar
      end if
c
c  Number seems to have terminated.  Have we got a legal number?
c  Not if we terminated in states 1, 2, 6 or 7.
c
      if ( ihave .eq. 1 .or. ihave .eq. 2 .or.
     &     ihave .eq. 6 .or. ihave .eq. 7 ) then
        ierror = ihave
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'S_TO_R8 - Serious error!'
        write ( *, '(a)' ) '  Illegal or nonnumeric input:'
        write ( *, '(a,a)' ) '    ', s
        return
      end if
c
c  Number seems OK.  Form it.
c
      if ( jtop .eq. 0 ) then
        rexp = 1.0D+00
      else
        if ( jbot .eq. 1 ) then
          rexp = 10.0D+00 ** ( jsgn * jtop )
        else
          rexp = 10.0D+00 ** ( dble ( jsgn * jtop ) / dble ( jbot ) )
        end if
      end if

      dval = dble ( isgn ) * rexp * rtop / rbot

      return
      end
      subroutine s_to_r8vec ( s, n, rvec, ierror )

c*********************************************************************72
c
cc S_TO_R8VEC reads an R8VEC from a string.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    28 April 2008
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, character * ( * ) S, the string to be read.
c
c    Input, integer N, the number of values expected.
c
c    Output, double precision RVEC(N), the values read from the string.
c
c    Output, integer IERROR, error flag.
c    0, no errors occurred.
c    -K, could not read data for entries -K through N.
c
      implicit none

      integer  n

      integer i
      integer ierror
      integer ilo
      integer lchar
      double precision rvec(n)
      character * ( * ) s

      i = 0
      ierror = 0
      ilo = 1

10    continue

      if ( i .lt. n ) then

        i = i + 1

        call s_to_r8 ( s(ilo:), rvec(i), ierror, lchar )

        if ( ierror .ne. 0 ) then
          ierror = -i
          go to 20
        end if

        ilo = ilo + lchar

        go to 10

      end if

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

