      Program Log2
c*********************************************************************72
c
c  LOG2 is a program which computes the integer part of the 
c  logarithm base two of the number X.  That is, LOG2 returns the 
c  number of times we can divide X by 2 before getting a value 
c  between 1 and 2.
c
c  If X is less than 1, we return a negative value, meaning we had 
c  to multiply rather than divide.
c
c  For example, LOG2(10)=3, because we can divide 10 by 2 just 3 
c  times before getting a result between 1 and 2.

      Integer Logtwo
      Real X

      Write (*,*) 'This program computes the integer part of'
      Write (*,*) 'the logarithm base 2 of any positive number.'

10    Continue

      Write (*,*) 'Enter a number X, or 0 to stop:'
      Read (*,*) X

      Logtwo = 0

      If (X .Gt. 0.0) Go To 20

c
c  If the user gave us a negative number, ask for a better one!
c
      If (X .Lt. 0.0) Then

        Write (*,*) 'LOG2 - Illegal input!'
        Write (*,*) 'The input value of X must be positive'
        Write (*,*) 'but your value was ', X
        Write (*,*) ' '
        Write (*,*) 'Try again!'
        Go to 10
c
c  But if the user gave us 0, quit.
c
      Else
        
        Write (*,*) 'The program is stopping because you typed 0.'
        Stop

      Endif

c
c  Here is where we work on getting the logarithm of a positive
c  number.
c

20    Continue

      If (X .Ge. 1.0 .And. X .Lt. 2.0) Then

        Write (*,*) 'LOG2 of X is ', Logtwo

        Go to 10

      Else
c
c  If a number is bigger than 2, divide it in half once.
c
        If (X .Ge. 2.0) Then

          X = X/2.0
          Logtwo = Logtwo + 1

          Go To 20

c
c  If a number is smaller than 1, double it.
c
        Elseif (X .Lt .1.0) Then

          X=2.0*X
          Logtwo=Logtwo-1
          Go To 20

        Endif

      Endif

      End
