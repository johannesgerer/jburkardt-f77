      Program Funsub

c*********************************************************************72

      Real Area
      Real Area1
      Real Area2
      Real Length
      Real Width

      External Area

      Length = 2.0
      Width = 5.3
      Write (*,*) 'Length=', Length, ' Width=', Width

      Call Box(Length, Width, Area1)

      Write (*,*) 'Box computes the area as ', Area1

      Area2 = Area(Length, Width)

      Write (*,*) 'Area computes the area as ', Area2

      Stop
      End

      Subroutine Box(Length, Width, Area)

c*********************************************************************72

      Real Area
      Real Length
      Real Width

      Area = Length * Width

      Return
      End

      Function Area(Length, Width)

c*********************************************************************72

      Real Area
      Real Length
      Real Width

      Area = Length * Width

      Return
      End




