      Program Count

c*********************************************************************72

      Integer Age
      Integer I
      Integer J
      Integer Number
      Integer Sum
c
      Age = 12

      Do 10 I = 1, Age

        Write (*,*) 'Happy birthday to you!'

10    Continue

      Sum = 0
      Number = 0
      Do 20 J = 1, Age

        Number = Number + 1
        Sum = Sum + Number

20    Continue

      Write (*,*) 'You''ve blown out ', Sum ,' candles in your life!'

      Stop
      End
