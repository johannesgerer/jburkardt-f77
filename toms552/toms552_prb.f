      program main
c
c sample driver program for subroutine cl1.
c
c this program solves a k by n overdetermined system
c
c    ax=b
c
c in the l1 sense subject to l equality constraints
c
c    cx=d
c
c and m inequality constraints
c
c    ex.le.f.
c
c complete details of the parameters may be
c found in the documentation of the subroutine.
c
c the arrays are currently dimensioned to allow problems
c for which k+l+m .le. 100, n .le. 10.
c
c the program may be tested on the following data.
c
c     k = 8
c     l = 3
c     m = 2
c     n = 5
c
c      q = 2  0  1  3  1  7
c          7  4  4 15  7  4
c          9  4  7 20  6  7
c          2  2  1  5  3  4
c          9  3  2 14 10  0
c          4  5  0  9  9  4
c          4  4  9 17 -1  9
c          1  6  2  9  5  6
c          0  4  5  9 -1  5
c          3  2  7 12 -2  1
c          3  6 12 21 -3  6
c          0  3  6  9 -3  5
c          6  2  4 12  4  6
c
c     kode = 0
c     toler = 1.e-5
c     iter = 130
c
c
      dimension q(102,12), x(12), res(100), cu(2,110)
      integer iu(2,110), s(100)
      data klmd, klm2d, nklmd, n2d /100,102,110,12/
c input data.
      read (5,99999) k, l, m, n, kode, toler, iter
      klm = k + l + m
      n1 = n + 1
      do 10 i=1,klm
         read (5,99998) (q(i,j),j=1,n1)
         write (6,99994) (q(i,j),j=1,n1)
   10 continue
      call cl1(k, l, m, n, klmd, klm2d, nklmd, n2d, q,
     * kode, toler, iter, x, res, error, cu, iu, s)
c output kode, iteration count and error norm.
      write (6,99997) kode, iter, error
c output solution vector.
      write (6,99996) (i,x(i),i=1,n)
c output residual error at each point.
      write (6,99995) (i,res(i),i=1,klm)
c
c  Terminate.
c
      stop
99999 format (5i3, e10.0, i3)
99998 format (8f3.0)
99997 format (16h kode,iter,error, 2i10, e18.7)
99996 format (4h sol, i5, e18.7)
99995 format (6h error, i5, e18.7)
99994 format (2h  , 8f5.0)
      end
