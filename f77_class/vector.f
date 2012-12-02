      Program Vector

c*********************************************************************72
c
c  This program demonstrates some simple operations with vectors.  
c  A vector is simply a set or list of data, stored under a single 
c  variable name.  Each item in the list can be accessed by 
c  specifying its "index".

      Integer I, Jumble(10), Product, Sofar(10)

c  Here's how a DO loop may be used to move through a vector.  In 
c  this case, we would like to assign a value to each entry of the 
c  vector, and we do this by specifying a "general" entry 
c  Jumble(I), where I is to go from 1 to 10.  Our formula, 2*I-1, 
c  will actually store successive odd numbers in the Jumble vector.

      Do 10 I = 1, 10
        Jumble(I) = 2*I - 1
10    Continue

c  You can print out a short array this way:

      Write (*,*) 'Here''s what happens when we print Jumble using'
      Write (*,*) 'the statement  Write (*,*) Jumble:'
      Write (*,*) ' '
      Write (*,*) Jumble

c  You can print out an array one entry at a time using a DO loop:

      Write(*,*)' '
      Write(*,*)'Here''s what happens when we print Jumble using'
      Write(*,*)'the statement  Write (*,*) Jumble(I):'
      Write(*,*)' '
      Do 20 I = 1, 10
        Write (*,*) Jumble(I)
20    Continue

c  You can use each entry in a vector, one entry at a time, using 
c  a DO loop.

      Product = 1
      Do 30 I = 1,10
        Product = Product * Jumble(I)
30    Continue
      Write(*,*) 'The product of all the entries is ', Product

c  What happens if we add up entries of Jumble?  Let's store the 
c  first entry, then the sum of the first and second, then the sum 
c  of the first, second and third, and so on.  We'll save the 
c  values in another vector.

      Sum=0
      Do 40 I=1,10
        Sum = Sum + Jumble(I)
        Sofar(I) = Sum
40    Continue

c  Let's print the sums out now

      Write(*,*)' '
      Do 50 I = 1, 10
        Write(*,*)'The sum of the first ',I,' entries is ', Sofar(I)
50    Continue

      Stop
      End
