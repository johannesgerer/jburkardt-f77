      program cash

c*********************************************************************72

      Integer Fives
      Integer Ones
      Integer Sum
      Integer Tens
      Integer Value

      Ones = 0
      Fives = 0
      Tens = 0
      Sum = 0

10    Continue

      Write (*,*) 'Enter denomination (1, 5, or 10), or 0 to stop'
      Read (*,*) Value

      If (Value .Eq. 1) Then

        Ones = Ones + 1
        Sum = Sum + 1
        Go to 10

      Elseif (Value .Eq. 5) Then

        Fives = Fives + 1
        Sum = Sum + 5
        Go to 10

      Elseif (Value .Eq. 10) Then

        Tens = Tens + 1
        Sum = Sum + 10
        Go to 10

      Else

        Write (*,*) 'The current total is ', sum
        Stop

      Endif

      End
