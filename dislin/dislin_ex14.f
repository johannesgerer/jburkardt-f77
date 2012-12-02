      program main

c*********************************************************************72
c
cc DISLIN_EX14 demonstrates the use of TeX instructions for math formulas.
c
c  Discussion:
c
c    On Unix systems, the backslash character is interpreted as a special
c    symbol BEFORE THE COMPILER SEES IT, unless it is inside a double-quote
c    string.
c
c    If you know the program will be used on a UNIX system, you can 'simply'
c    type TWO backslash characters whenever you mean one.
c
c  Modified:
c
c    09 April 2011
c
c  Reference:
c
c    Helmut Michels,
c    The Data Plotting Software DISLIN - version 10.4,
c    Shaker Media GmbH, January 2010,
c    ISBN13: 978-3-86858-517-9.
c
      implicit none

      character*80 cstr
      integer nl
      integer nlmess

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'DISLIN_EX14:'
      write ( *, '(a)' ) '  FORTRAN77 version'
      write ( *, '(a)' ) '  Demonstrate the use of TeX instructions to'
      write ( *, '(a)' ) '  create mathematical formulas.'
c
c  Specify the format of the output file.
c
      call metafl ( 'png' )
c
c  Indicate that new data overwrites old data.
c
      call filmod ( 'delete' )
c
c  Specify the name of the output graphics file.
c
      call setfil ( 'dislin_ex14.png' )
c
c  Choose the page size and orientation.
c
      call setpag ( 'usap' )
c
c  For PNG output, reverse the default black background to white.
c
      call scrmod ( 'reverse' )
c
c  Open DISLIN.
c
      call disini ( )
c
c  Plot a border around the page.
c
      call pagera ( )
c
c  Use the COMPLEX font.
c
      call complx ( )
      call height ( 40 )

      cstr = 'TeX Instructions for Mathematical Formulas'
      nl = nlmess ( cstr )
      call messag ( cstr, (2100 - nl)/2, 100 ) 
  
      call texmod ( 'ON' )
      call messag ( '$\frac{1}{x+y}$', 150, 400 )
      call messag ( '$\frac{a^2 - b^2}{a+b} = a - b$', 1200, 400 )
  
      call messag ( '$r = \red{\sqrt{x^2 + y^2}}', 150, 700 )
      call messag ( '$\cos \phi = \frac{x}{\sqrt{x^2 + y^2}}$', 
     &  1200, 700 )

      call messag ( '$\Gamma(x) = \int_0^\infty e^{-t}t^{x-1}dt$', 
     &  150, 1000 )
      call messag ( '$\lim_{x \to \infty} (1 + \frac{1}{x})^x = e$', 
     &  1200, 1000 )

      call messag ( '$\mu = \sum_{i = 1}^n x_i p_i$', 150, 1300 )
      call messag ( '$\mu = \int_{-\infty}^ \infty x f(x) dx$', 
     &  1200, 1300 )

      call messag ( 
     &  '$\overline{x} = \frac{1}{n} \sum_{i = 1}^n x_i$', 
     &  150, 1600 )

      call messag ( '$s^2 = \frac{1}{n-1} \sum_{i = 1}^n' //
     &  '(x_i - \overline{x})^2$', 1200, 1600 )

      call messag ( '$\sqrt[n]{\frac{x^n - y^n}{1 + u^{2n}}}$', 
     &  150, 1900 )  
      call messag ( '$\sqrt[3]{-q + \sqrt{q^2 + p^3}}$', 1200, 1900 )

      call messag ( '$\int \frac{dx}{1+x^2} = \arctan x + C$', 
     &  150, 2200 )

      call messag ( 
     &  '$\int \frac{dx}{\sqrt{1+x^2}} = {\rm arcsinh} x + C$',
     &  1200, 2200 )

      call messag ( '$\overline{P_1P_2} = \sqrt{(x_2-x_1)^2 + ' //
     &  '(y_2-y_1)^2}$', 150, 2500 )
      call messag ( '$x = \frac{x_1 + \lambda x_2}{1 + \lambda}$', 
     &  1200, 2500 )
c
c  Close DISLIN.
c
      call disfin ( )
c
c  Terminate.
c
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'DISLIN_EX14:'
      write ( *, '(a)' ) '  Normal end of execution.'

      stop
      end
