      Program Table

c*********************************************************************72

      Integer I
      Integer J
      Integer X(3,4)
c
c  Set the values of the array.
c
      Do 20 I = 1, 3
        Do 10 J = 1, 4

          X(I,J) = 10*I + J

10      Continue
20    Continue

c
c  Print the values of the array.
c
      Do 30 I = 1, 3
        Write (*,*) (X(I,J), J = 1, 4)
30    Continue

      Stop
      End
