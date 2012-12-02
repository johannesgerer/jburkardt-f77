      Program Life

c*********************************************************************72
c
c  On a 10 by 10 board, place 0's and 1's randomly, leaving the border 0. 
c  1 is a live cell, 0 a dead cell.  On each step, a living cell ``dies'' if it 
c  has 0, or 1, or more than 3 living cell neighbors.
c  On each step, an dead cell comes alive if it has exactly 3 living
c  cell neighbors.

      Integer Nrow,Ncol
      Parameter (Nrow=10)
      Parameter (Ncol=10)
c
      Integer Cell(Nrow,Ncol)
      Integer I
      Integer Istep
      Integer Old(Nrow,Ncol)

      Call Init(Cell, Old, Nrow, Ncol)
      Istep=0
      Call Display(Cell, Istep, Nrow, Ncol)

      Do 10 I = 1, 5
        Call Compute(Cell, Old, Nrow, Ncol)
        Istep=I
        Call Display(Cell, Istep, Nrow, Ncol)
10    Continue

      Stop
      End

      Subroutine Init(Cell, Old, Nrow, Ncol)

c*********************************************************************72
c
c  Initialize the values in Cell to random 0 and 1 values.
c  Make sure the boundary is 0.

      Integer Nrow, Ncol
c
      Integer Cell(Nrow,Ncol)
      Integer I
      Integer Iseed, J, Old(Nrow,Ncol)
      Real Random
c
      External Random
c
      Iseed = 1993
c
      Do 20 I = 1, Nrow
        Do 10 J = 1, Ncol

          If (I .Eq. 1 .Or. I .Eq. Nrow .Or.
     &        J .Eq. 1 .Or. J .Eq .Ncol) Then
            Cell(I,J) = 0
          Else
            Cell(I,J) = Nint(Random(Iseed))
          Endif

          Old(I,J) = Cell(I,J)

10      Continue
20    Continue

      Return
      End

      Subroutine Compute(Cell, Old, Nrow, Ncol)

c*********************************************************************72
c
c  Update the cell information for the next step.

      Integer Nrow, Ncol
c
      Integer Newval, Cell(Nrow,Ncol)
      Integer I, J, Old(Nrow,Ncol)

      Do 20 I = 1, Nrow
        Do 10 J = 1, Ncol
          Old(I,J) = Cell(I,J)
10      Continue
20    Continue

      Do 40 I = 1, Nrow
        Do 30 J = 1, Ncol

          If (I .Eq. 1 .Or. I .Eq. Nrow .Or. 
     &        J .Eq. 1 .Or. J .Eq. Ncol) Then 

            Newval=0 

          Else

            Nabors=Old(I-1,J-1)+Old(I-1,J  )+Old(I-1,J+1)+
     &             Old(I,  J-1)             +Old(I,  J+1)+
     &             Old(I+1,J-1)+Old(I+1,J  )+Old(I+1,J+1)

            If (Nabors .Le. 1) Then
              Newval = 0
            Elseif (Nabors .Eq. 2) Then
              Newval = Old(I,J)
            Elseif (Nabors .Eq. 3) Then
              Newval = 1
            Else
              Newval = 0
            Endif

          Endif

          Cell(I,J) = Newval

30      Continue
40    Continue

      Return
      End

      Subroutine Display(Cell, Istep, Nrow, Ncol)

c*********************************************************************72
c
c  Display prints out the cell information.

      Integer Nrow, Ncol
c
      Integer Cell(Nrow,Ncol), I, Istep

      Write (*,*) ' '
      Write (*,*) 'Step ',Istep
      Write (*,*) ' '
      Do 10 I = 1, Nrow
        Write (*,'(1x,10i1)') (Cell(I,J), J = 1, Ncol)
10    Continue

      Return
      End

      Function Random(Iseed)

c*********************************************************************72
c
c  Random returns a "random" value.  Before the first call to
c  Random, set Iseed to some large, nonzero value.
c
      Integer Iseed
      Real Random
c
      Intrinsic Real
c
      Iseed = Iseed * 125

      Iseed = Iseed - (Iseed/2796203) * 2796203

      Random = Real(Iseed) / 2796203.0
c
      Return
      End
