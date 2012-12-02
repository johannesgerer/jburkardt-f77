      Program Bill

c*********************************************************************72

      Real Tea, Coffee, Water, Chicken
      Real Beef, Cake, Pie
      Real Price1, Price2
c
c  Set the prices for food.
c
      Tea      = .50
      Coffee   = .45
      Water    = 0.0

      Chicken  = 1.75
      Beef     = 2.50 

      Cake     = 1.50
      Pie      = 1.25
      Icecream = 1.75     

      Call Waiter(Tea, Chicken, Pie, Price1)
      Write (*,*) 'First order will cost ', Price1

      Call Waiter(Water, Beef, Cake, Price2)
      Write (*,*) 'Second order will cost ', Price2

      Stop
      End

      Subroutine Waiter(Item1, Item2, Item3, Cost)

c*********************************************************************72

      Real Item1, Item2, Item3
      Real Total, Tip, Cost

      Total = Item1 + Item2 + Item3
      Tip = 0.15 * Total
      Cost = Total + Tip

      Return
      End
