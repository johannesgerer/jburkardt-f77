      subroutine altkos ( n, umspei, z1, unendl, izaehl, jzaehl, 
     &  alterk )

c*********************************************************************72
c
cc ALTKOS is used by QAP to determine alternative costs.
c
c  Discussion:
c
c    The routine determines alternative costs and the resulting
c    single assignment.
c
c  Modified:
c
c    10 December 2007
c
c  Author:
c
c    T Boenniger,
c    Rainer Burkard,
c    Karl-Heinz Stratmann
c
c  Reference:
c
c    Rainer Burkard, Ulrich Derigs,
c    Assignment and Matching Problems: Solution Methods with Fortran-Programs,
c    Springer, 1980,
c    ISBN: 3-540-10267-1,
c    LC: QA402.5 B86.
c
c  Parameters:
c
c    Input, integer N, ?
c
c    Input, integer UMSPEI(*), ?
c
c    Input, integer Z1(N), ?
c
c    Input, integer UNENDL, a large value.
c
c    Output, integer IZAEHL, JZAEHL, ?
c
c    Output, integer ALTERK, the alternative cost.
c
      implicit none

      integer n

      integer a
      integer alterk
      integer i
      integer i1
      integer izaehl
      integer j
      integer j1
      integer jj
      integer jzaehl
      integer min
      integer min1
      integer umspei(*)
      integer unendl
      integer z1(n)

      alterk = -1
      jj = 0

      do i = 1, n

        j = z1(i)
        j1 = j - n
        min = unendl

        do i1 = 1, n

          j1 = j1 + n

          if ( i1 .ne. i ) then
            a = umspei(j1)
            if ( a .lt. min ) then
              min = a
            end if
          end if

        end do

        min1 = min
        min = unendl

        do i1 = 1, n

          jj = jj + 1

          if ( i1 .ne. j ) then
            a = umspei(jj)
            if ( a .lt. min ) then
              min = a
            end if
          end if

        end do

        min = min + min1

        if ( alterk .lt. min ) then
          izaehl = i
          jzaehl = j
          alterk = min
        end if

      end do

      return
      end
      subroutine bmp ( n, cc, sup, nmatch, z, basis, mem, stack, sm,
     &  tma, tmb, p )

c*********************************************************************72
c
cc BMP solves the bottleneck matching problem.
c
c  Modified:
c
c    23 November 2007
c
c  Author:
c
c    Ulrich Derigs,
c    G Kazakidis
c
c  Reference:
c
c    Rainer Burkard, Ulrich Derigs,
c    Assignment and Matching Problems: Solution Methods with Fortran-Programs,
c    Springer, 1980,
c    ISBN: 3-540-10267-1,
c    LC: QA402.5 B86.
c
c  Parameters:
c
c    Input, integer N, the number of nodes.  N should be even.
c
c    Iput, integer CC(N*(N-1)/2), the strict upper triangle of the cost matrix,
c    stored column by column.
c
c    Input, integer SUP, a large integer.
c
c    Output, integer NMATCH(N), the optimal matching.
c
c    Output, integer Z, the cost of the optimal matching.
c
c    Workspace, integer BASIS(N), MEM(N), STACK(N), SM(N), TMA(N),
c    TMB(N), P(N).
c
      implicit none

      integer n

      integer basis(n)
      integer cc((n*(n-1))/2)
      integer cnsnt
      integer i
      integer i1
      integer ii
      integer iii
      integer ind
      integer is
      integer it
      integer iz
      integer j
      integer k
      integer max
      integer mem(n)
      integer min
      integer n1
      integer n2
      integer nb
      integer nb1
      integer nb2
      integer nb3
      integer nbs
      integer nbt
      integer ncard
      integer nhalb
      integer nk
      integer nk1
      integer nk2
      integer nka
      integer nkb
      integer nm
      integer nmatch(n)
      integer ns
      integer nt
      integer p(n)
      integer sm(n)
      integer stack(n)
      integer sup
      integer tma(n)
      integer tmb(n)
      integer top
      integer z
      integer zmin

      p(1) = 0
      p(2) = 0
      do i = 3, n
        p(i) = p(i-1) + i - 2
      end do

      do i = 1, n
        nmatch(i) = 0
      end do

      ncard = 0
      nhalb = n / 2
      top = n + 2
c
c  Start procedure.
c
      z = cc(1)

      do j = 3, n
        i1 = p(j) + 1
        if ( cc(i1) .lt. z ) then
          z = cc(i1)
        end if
      end do

      do i = 2, n

        i1 = p(i) + 1
        zmin = cc(i1)

        if ( zmin .le. z ) then
          go to 30
        end if

        ii = i - 1

        do j = 2, ii
          i1 = p(i) + j
          if ( cc(i1) .le. z ) then
            go to 30
          end if
          if ( cc(i1) .lt. zmin ) then
            zmin = cc(i1)
          end if
        end do

        ii = i + 1

        do j = ii, n
          i1 = p(j) + i
          if ( cc(i1) .le. z ) then
            go to 30
          end if
          if ( cc(i1) .lt. zmin ) then
            zmin = cc(i1)
          end if
        end do

        z = zmin

30      continue

      end do

      k = 0

      do j = 2, n
        i1 = p(j) + 1
        if ( cc(i1) .le. z ) then
          k = k + 1
        end if
      end do

      basis(1) = k
      mem(1) = 1

      do i = 2, n
        k = 0
        ii = i - 1
        do j = 1, ii
          i1 = p(i) + j
          if ( cc(i1) .le. z ) then
            k = k + 1
          end if
        end do

        ii = i + 1

        do j = ii, n
          i1 = p(j) + i
          if ( cc(i1) .le. z ) then
            k = k + 1
          end if
        end do

        basis(i) = k
        mem(i) = i

      end do

      call ssort ( basis, mem, n )

      do ii = 1, n

        i = mem(ii)

        if ( nmatch(i) .gt. 0 ) then
          go to 90
        end if

        iii = i - 1

        do j = 1, iii
          if ( nmatch(j) .eq. 0 ) then
            i1 = p(i) + j
            if ( cc(i1) .le. z ) then
              go to 85
            end if
          end if
        end do

        iii = i + 1

        do j = iii, n
          if ( nmatch(j) .eq. 0 ) then
            i1 = p(j) + i
            if ( cc(i1) .le. z ) then
              go to 85
            end if
          end if
        end do

        go to 90

85      continue

        nmatch(i) = j
        nmatch(j) = i
        ncard = ncard + 1

90      continue

      end do

      if ( ncard .ge. nhalb ) then
        return
      end if

110   continue

      zmin = sup
      it = 1

      do i = 1, n
        basis(i) = i
        mem(i) = i
        sm(i) = top
        tma(i) = top
        tmb(i) = top
        if ( nmatch(i) .eq. 0 ) then
          sm(i) = 0
        end if
      end do
c
c  A loop like this assumes that the value of I can
c  be "carried" out of the loop.
c
      do i = 1, n
        if ( nmatch(i) .eq. 0 ) then
          go to 130
        end if
      end do

130   continue

      stack(it) = i
      it = it + 1
      is = 1
c
c  Compute the shortest augmenting path starting from vertex I
c  by growing an alternating tree rooted at vertex I.
c
200   continue

      ns = stack(is)
      nbs = basis(ns)
      tma(nbs) = 0

      do nt = 1, n

        if ( ns .eq. nt ) then
          go to 240
        end if

        nbt = basis(nt)

        if ( tma(nbt) .ne. top ) then
          go to 240
        end if

        min = ns
        max = nt

        if ( ns .gt. nt ) then
          max = ns
          min = nt
        end if

        ind = p(max) + min
        cnsnt = cc(ind)

        if ( cnsnt .le. z ) then
          go to 219
        end if

        if ( cnsnt .ge. zmin ) then
          go to 240
        end if

        zmin = cnsnt
        go to 240

219     continue
c
c  Growing the alternating tree by adding two vertices and edges.
c
        if ( sm(nbt) .eq. top ) then

          tma(nbt) = ns
          tmb(nbt) = nbt
          nb = nmatch(nbt)
          sm(nb) = nbt
          stack(it) = nb
          it = it + 1
          go to 240

        end if

        nka = nt
        nkb = ns
        n1 = nbs
        nb1 = nbs
        n2 = nbt
        nb2 = n2

220     continue

        tma(nb1) = nb2
        nk = sm(nb1)

        if ( nk .eq. 0 ) then
          go to 225
        end if

        nb2 = basis(nk)
        nb1 = tma(nb2)
        nb1 = basis(nb1)
        go to 220

225     continue

        nb = nb1
        nb1 = n2
        nb2 = n1

230     continue

        if ( tma(nb1) .lt. top ) then
          go to 235
        end if

        tma(nb1) = nb2
        nk = sm(nb1)

        if ( nk .eq. 0 ) then
          go to 600
        end if

        nb2 = basis(nk)
        nb1 = tma(nb2)
        nb1 = basis(nb1)
        go to 230

235     continue
c
c  Shrinking a blossom.
c
        if ( nb1 .eq. nb ) then

400       continue

          nk1 = nb
          nk = mem(nb)

          if ( nb .ne. n2 ) then
            go to 436
          end if

435       continue

          n2 = n1
          nb2 = tma(nb)

436       continue

          mem(nk1) = nb2
          nm = nmatch(nb2)
          sm(nb2) = nm
          stack(it) = nb2
          it = it + 1
          nk1 = nb2

440       continue

          nk2 = nk1
          basis(nk2) = nb
          nk1 = mem(nk2)

          if ( nk1 .ne. nb2 ) then
            go to 440
          end if

          nb1 = basis(nm)
          mem(nk2) = nb1
          nk2 = nb1
      
445       continue

          nk1 = nk2
          basis(nk1) = nb
          nk2 = mem(nk1)

          if ( nk2 .ne. nb1 ) then
            go to 445
          end if

          if ( n2 .eq. nb1 ) then
            go to 450
          end if

          nb2 = tma(nb1)
          tma(nb1) = tmb(nb2)
          tmb(nb1) = tma(nb2)
          go to 436

450       continue

          if ( n2 .eq. n1 ) then
            go to 455
          end if

          tma(n2) = nkb
          tmb(n2) = nka

          if ( nb .ne. n1 ) then
            go to 435
          end if

          go to 460

455       continue

          tma(n1) = nka
          tmb(n1) = nkb

460       continue

          mem(nk1) = nk
          nbs = nb
          tma(nbs) = 0
          go to 240

        end if

        nk = tma(nb)
        tma(nb) = top
        nm = nmatch(nk)
        nb = basis(nm)
        go to 235

240     continue

      end do

      tma(nbs) = top
      is = is + 1

      if ( is .lt. it ) then
        go to 200
      end if

      z = zmin
      go to 110
c
c  Augmentation of the matching.
c  Exchange of the matching and nonmatching edges along
c  the augmenting path.
c
600   continue

      nb = nkb
      nb2 = nka

605   continue

      nb1 = nb

610   continue

      nmatch(nb1) = nb2
      nb3 = sm(nb1)

      if ( nb3 .eq. 0 ) then
        go to 620
      end if

      nb2 = tmb(nb3)

      if ( nb2 .ne. nb3 ) then
        go to 625
      end if

615   continue

      nb1 = tma(nb3)
      nmatch(nb2) = nb1
      go to 610

620   continue

      if ( nb .ne. nkb ) then
        go to 650
      end if

      nb = nka
      nb2 = nkb
      go to 605

625   continue

      iz = 1
      nk2 = nb3
      nk1 = nb2

630   continue

      stack(iz) = nk2
      nk = nk2

635   continue

      nk2 = sm(nk1)
      nk1 = tmb(nk2)

      if ( nk1 .eq. nk2 ) then
        go to 640
      end if

      iz = iz + 1
      go to 630

640   continue

      nk1 = tma(nk2)
      nmatch(nk2) = nk1
      nmatch(nk1) = nk2

      if ( nk1 .ne. nk ) then
        go to 635
      end if

645   continue

      if ( iz .eq. 1 ) then
        go to 615
      end if

      nk1 = tma(nk)
      nk2 = tmb(nk)
      nmatch(nk1) = nk2
      nmatch(nk2) = nk1
      iz = iz - 1
      nk = stack(iz)

      if ( nk1 .eq. nk ) then
        go to 645
      end if

      go to 635

650   continue

      ncard = ncard + 1

      if ( ncard .lt. nhalb ) then
        go to 110
      end if

      return
      end
      subroutine cmp ( n, m, nbl, index, nmatch, ncard, basis, mem, 
     &  stack, sm, tma, tmb )

c*********************************************************************72
c
cc CMP solves the cardinality matching problem.
c
c  Modified:
c
c    22 November 2007
c
c  Author:
c
c    Ulrich Derigs,
c    G Kazakidis
c
c  Reference:
c
c    Rainer Burkard, Ulrich Derigs,
c    Assignment and Matching Problems: Solution Methods with Fortran-Programs,
c    Springer, 1980,
c    ISBN: 3-540-10267-1,
c    LC: QA402.5 B86.
c
c  Parameters:
c
c    Input, integer N, the number of nodes.
c
c    Input, integer M, the number of edges.
c
c    Input, integer NBL(2*M), cumulated list of neighbors.
c
c    Input, integer INDEX(N+1), index of the first neighbor of node I
c    in vector NBL.  Note that the extra final entry of INDEX should point
c    to the first unused entry of NBL, namely 2*M+1!
c
c    Output, integer NMATCH(N), optimal matching.
c    0, if node I is unsaturated.
c    J, if edge (I,J) is the matching edge.
c
c    Output, integer NCARD, the cardinality of the optimal matching.
c
c    Workspace, integer BASIS(N), MEM(N), STACK(N), SM(N), TMA(N), TMB(N).
c
      implicit none

      integer m
      integer n

      integer anf
      integer basis(n)
      integer ende
      integer i
      integer ii
      integer index(n+1)
      integer ipe
      integer is
      integer it
      integer iz
      integer j
      integer jj
      integer mem(n)
      integer n1
      integer n2
      integer nb
      integer nb1
      integer nb2
      integer nb3
      integer nbl(2*m)
      integer nbs
      integer nbt
      integer ncard
      integer nexp
      integer ni
      integer nk
      integer nk1
      integer nk2
      integer nka
      integer nkb
      integer nm
      integer nmatch(n)
      integer ns
      integer nt
      integer sm(n)
      integer stack(n)
      integer tma(n)
      integer tmb(n)
      integer top

      do i = 1, n
        nmatch(i) = 0
      end do

      ncard = 0
      top = n + 2

      do i = 1, n
        basis(i) = index(i+1) - index(i)
        mem(i) = i
      end do

      call ssort ( basis, mem, n )

      do ii = 1, n

        i = mem(ii)

        if ( nmatch(i) .le. 0 ) then

          anf = index(i)
          ende = index(i+1) - 1

          do jj = anf, ende

            j = nbl(jj)

            if ( nmatch(j) .le. 0 ) then
              nmatch(j) = i
              nmatch(i) = j
              ncard = ncard + 1
              go to 55
            end if

          end do

        end if

55      continue

      end do

      nexp = n - 2 * ncard

      if ( nexp .lt. 2 ) then
        return
      end if

110   continue

      it = 1

      do i = 1, n

        basis(i) = i
        mem(i) = i
        sm(i) = top
        tma(i) = top
        tmb(i) = top

        if ( nmatch(i) .eq. 0 ) then
          sm(i) = 0
          stack(it) = i
          it = it + 1
        end if

      end do

      is = 1
c
c  Growing alternating trees from each unsaturated vertex.
c
200   continue

      ns = stack(is)
      nbs = basis(ns)
      tma(nbs) = 0
      i = index(ns)
      ii = ns + 1
      ii = index(ii) - 1

      do ni = i, ii

        nt = nbl(ni)
        nbt = basis(nt)

        if ( tma(nbt) .eq. top ) then
c
c  Growing one alternating tree by adding two vertices and edges.
c
          if ( sm(nbt) .eq. top ) then
            tma(nbt) = ns
            tmb(nbt) = nbt
            nb = nmatch(nbt)
            sm(nb) = nbt
            stack(it) = nb
            it = it + 1
            go to 240
          end if

          nka = nt
          nkb = ns
          n1 = nbs
          nb1 = n1
          n2 = nbt
          nb2 = n2

220       continue

          tma(nb1) = nb2
          nk = sm(nb1)

          if ( nk .eq. 0 ) then
            go to 225
          end if

          nb2 = basis(nk)
          nb1 = tma(nb2)
          nb1 = basis(nb1)

          go to 220

225       continue

          nb = nb1
          nb1 = n2
          nb2 = n1

230       continue

          if ( tma(nb1) .lt. top ) then
            go to 235
          end if

          tma(nb1) = nb2
          nk = sm(nb1)

          if ( nk .eq. 0 ) then
            go to 600
          end if

          nb2 = basis(nk)
          nb1 = tma(nb2)
          nb1 = basis(nb1)

          go to 230

235       continue
c
c  Shrinking a blossom.
c
          if ( nb1 .eq. nb ) then

            nk1 = nb
            nk = mem(nb)

            if ( nb .ne. n2 ) then
              go to 436
            end if

435         continue

            n2 = n1
            nb2 = tma(nb)

436         continue

            mem(nk1) = nb2
            nm = nmatch(nb2)
            sm(nb2) = nm
            stack(it) = nb2
            it = it + 1
            nk1 = nb2

440         continue

            nk2 = nk1
            basis(nk2) = nb
            nk1 = mem(nk2)

            if ( nk1 .ne. nb2 ) then
              go to 440
            end if

            nb1 = basis(nm)
            mem(nk2) = nb1
            nk2 = nb1

445         continue

            nk1 = nk2
            basis(nk1) = nb
            nk2 = mem(nk1)

            if ( nk2 .ne. nb1 ) then
              go to 445
            end if

            if ( n2 .eq. nb1 ) then
              go to 450
            end if

            nb2 = tma(nb1)
            tma(nb1) = tmb(nb2)
            tmb(nb1) = tma(nb2)
            go to 436

450         continue

            if ( n2 .eq. n1 ) then
              go to 455
            end if

            tma(n2) = nkb
            tmb(n2) = nka

            if ( nb .ne. n1 ) then
              go to 435
            end if

            go to 460

455         continue

            tma(n1) = nka
            tmb(n1) = nkb

460         continue

            mem(nk1) = nk
            nbs = nb
            tma(nbs) = 0
            go to 240

          end if

          nk = tma(nb)
          tma(nb) = top
          nm = nmatch(nk)
          nb = basis(nm)

          go to 235

        end if

240     continue

      end do

      tma(nbs) = top
      is = is + 1

      if ( is .lt. it ) then
        go to 200
      end if

      return
c
c  Augmentation of the matching.
c  Exchange of the matching and non-matching edges along the
c  augmenting path.
c
600   continue

      nb = nkb
      nb2 = nka

605   continue

      nb1 = nb

610   continue

      nmatch(nb1) = nb2
      nb3 = sm(nb1)

      if ( nb3 .eq. 0 ) then
        go to 620
      end if

      nb2 = tmb(nb3)

      if ( nb2 .ne. nb3 ) then
        go to 625
      end if

615   continue

      nb1 = tma(nb3)
      nmatch(nb2) = nb1
      go to 610

620   continue

      if ( nb .ne. nkb ) then
        go to 650
      end if

      nb = nka
      nb2 = nkb
      go to 605

625   continue

      iz = 1
      nk2 = nb3
      nk1 = nb2

630   continue

      stack(iz) = nk2
      nk = nk2

635   continue

      nk2 = sm(nk1)
      nk1 = tmb(nk2)

      if ( nk1 .eq. nk2 ) then
        go to 640
      end if

      iz = iz + 1
      go to 630

640   continue

      nk1 = tma(nk2)
      nmatch(nk2) = nk1
      nmatch(nk1) = nk2

      if ( nk1 .ne. nk ) then
        go to 635
      end if

645   continue

      if ( iz .eq. 1 ) then
        go to 615
      end if

      nk1 = tma(nk)
      nk2 = tmb(nk)
      nmatch(nk1) = nk2
      nmatch(nk2) = nk1
      iz = iz - 1
      nk = stack(iz)

      if ( nk1 .eq. nk ) then
        go to 645
      end if

      go to 635

650   continue

      nexp = nexp - 2
      ncard = ncard + 1

      if ( nexp .gt. 1 ) then
        go to 110
      end if

      return
      end
      subroutine connect ( n, nb, index, label, stack, icon )

c*********************************************************************72
c
cc CONNECT checks the connectivity of a graph.
c
c  Discussion:
c
c    For the CPP routine, it is essential that the graph is connected.
c    This routine is used to check that condition.
c
c  Modified:
c
c    05 December 2007
c
c  Author:
c
c    W Puetz
c
c  Reference:
c
c    Rainer Burkard, Ulrich Derigs,
c    Assignment and Matching Problems: Solution Methods with Fortran-Programs,
c    Springer, 1980,
c    ISBN: 3-540-10267-1,
c    LC: QA402.5 B86.
c
c  Parameters:
c
c    Input, integer N, the number of nodes.
c
c    Input, integer NB(2*M), the cumulated liswt of neighbors.  Here M
c    is the number of edges.  NB contains the neighbors of node 1, 
c    followed by the neighbors of node 2, and so on.
c
c    Input, integer INDEX(N+1), the index of the first neighbor of node
c    I in the neighbor list NB.
c
c    Workspace, integer LABEL(N).
c
c    Workspace, integer STACK(N).
c
c    Output, integer ICON, is 1 if the graph is connected, and
c    0 otherwise.
c
      implicit none

      integer n

      integer i
      integer i1
      integer i2
      integer icon
      integer index(n+1)
      integer j
      integer jj
      integer kstack
      integer label(n)
      integer nb(*)
      integer stack(n)

      label(1) = 0
      do i = 2, n
        label(i) = i
      end do

      icon = 1
      stack(1) = 1
      kstack = 1

15    continue

      i = stack(kstack)
      kstack = kstack - 1
      i1 = index(i)
      i2 = i + 1
      i2 = index(i2) - 1

      do j = i1, i2

        jj = nb(j)

        if ( label(jj) .ne. 0 ) then

          icon = icon + 1

          if ( icon .eq. n ) then
            icon = 1
            return
          end if

          label(jj) = 0
          kstack = kstack + 1
          stack(kstack) = jj

        end if

      end do

      if ( kstack .ne. 0 ) then
        go to 15
      end if

      icon = icon / n

      return
      end
      subroutine cpp ( n, m, kst, nb, kost, index, top, basis, mem,
     &  ka, kb, sm, tma, tmb, y1, y2, dplus, dminus, kurs, eps )

c*********************************************************************72
c
cc CPP solves the Chinese Postman Problem.
c
c  Modified:
c
c    07 December 2007
c
c  Author:
c
c    W Puetz
c
c  Reference:
c
c    Rainer Burkard, Ulrich Derigs,
c    Assignment and Matching Problems: Solution Methods with Fortran-Programs,
c    Springer, 1980,
c    ISBN: 3-540-10267-1,
c    LC: QA402.5 B86.
c
c  Parameters:
c
c    Input, integer N, the number of nodes.
c
c    Input, integer M, the number of edges.
c
c    Output, integer KST, the total cose of the duplicated edges.
c
c    Input, integer NB(2*M), the cumulated list of node neighbors.
c
c    Input/output, integer KOST(2*M).  On input KOST(I) is the cost of edge 
c    corresponding to edge to node NB(I).  On output, the cumulated next
c    node list of the optimal postman tour starting at node KURS.
c
c    Input/output, integer INDEX(N+1), the index of the first neighbor of
c    node I in the NB array.  On output, the index of the first neighbor
c    of node I in the next node list.
c
c    Input, integer TOP, a large integer.
c
c    Workspace, integer BASIS(N), MEM(N), KA(N), KB(N), SM(N), TMA(N),
c    TMB(N).
c
c    Workspace, real Y1(N), Y2(N), DPLUS(N), DMINUS(N).
c
c    Input, integer KURS, the starting node of the tour.
c
c    Input, real EPS, the machine accuracy.
c
      implicit none

      integer m
      integer n

      integer b
      integer b1
      integer b2
      integer b3
      integer base
      integer basis(n)
      integer bb
      integer bbb1
      integer bbest
      integer best
      real d
      real d1
      real dminus(n)
      real dplus(n)
      real eps
      real fltop
      integer i
      integer ii
      integer index(n+1)
      integer ka(n)
      integer kant
      integer kb(n)
      integer kn
      integer kn1
      integer kn2
      integer kn3
      integer kn4
      integer kost(2*m)
      integer kst
      integer kurs
      integer mem(n)
      integer nb(2*m)
      integer nm
      integer np
      integer sm(n)
      integer tma(n)
      integer tmb(n)
      integer top
      real y
      real y1(n)
      real y2(n)

      np = n + 1
      nm = - np

      do kn = 1, n
        dminus(kn) = top
      end do
c
c  Initial labeling.
c
      do kn = 1, n

        basis(kn) = kn
        mem(kn) = kn
        sm(kn) = nm
        tma(kn) = 0
        kb(kn) = kn
        y1(kn) = 0.0E+00
        y2(kn) = 0.0E+00
        i = index(kn)
        ii = index(kn+1)
        b = ii - i
        b1 = ( b / 2 ) * 2

        if ( b .ne. b1 ) then

          ii = ii - 1
          sm(kn) = 0
          dplus(kn) = 0.0E+00

          do kn1 = i, ii

            kn2 = nb(kn1)
            d1 = float ( kost(kn1) )

            if ( d1 .lt. dminus(kn2) ) then
              ka(kn2) = kn
              dminus(kn2) = d1
            end if

          end do

        end if

      end do
c
c  Examination of the labeling and decision of next step.
c
2001  continue

      d = top

      do b1 = 1, n

        if ( basis(b1) .eq. b1 ) then

          d1 = dminus(b1)

          if ( 0 .le. sm(b1) ) then

            d1 = 0.5E+00 * ( d1 + dplus(b1) )

            if ( d1 .le. d ) then
              d = d1
              best = b1
            end if

          else

            if ( 0 .lt. tma(b1) ) then
              d1 = d1 + y1(b1)
            end if

            if ( d1 .lt. d ) then
              d = d1
              best = b1
            end if

          end if

        end if

      end do
c
c  Check to see that at least once we picked up a better value of D.
c
      fltop = float ( top ) / 2.0E+00

      if ( fltop .le. d ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'CPP - Fatal error!'
        write ( *, '(a)' ) '  The input data is not correct, or'
        write ( *, '(a)' ) '  perhaps TOP is too small.'
        write ( *, '(a)' ) ' '
        write ( *, '(a,i12)' ) '  TOP = ', top
        write ( *, '(a,g14.6)' ) '  D = ', d
        stop
      end if

      if ( sm(best) .lt. 0 ) then
        go to 2010
      end if

      kn1 = ka(best)
      kn2 = kb(best)
      b = basis(kn1)
      b1 = best
      b2 = b

2020  continue

      tma(b1) = b2
      kn = sm(b1)

      if ( kn .ne. 0 ) then
        b2 = basis(kn)
        kn = tma(b2)
        b1 = basis(kn)
        go to 2020
      end if

      base = b1
      b1 = b
      b2 = best

2040  continue

      if ( tma(b1) .le. 0 ) then
        tma(b1) = b2
        kn = sm(b1)
        if ( kn .eq. 0 ) then
          go to 7000
        end if
        b2 = basis(kn)
        kn = tma(b2)
        b1 = basis(kn)
        go to 2040
      end if

2050  continue

      if ( b1 .eq. base ) then
        go to 4000
      end if

      b3 = tma(base)
      tma(base) = 0
      kn3 = - sm(b3)
      base = basis(kn3)
      go to 2050

2010  continue

      if ( tma(best) .le. 0 ) then
        go to 3000
      end if

      b = mem(best)

      if ( b .eq. best ) then
        go to 6000
      end if

      kn = tma(b)

      if ( 0 .lt. kn ) then
        go to 5000
      end if

      go to 6000
c
c  Growing an alternating tree.
c
3000  continue

      kn = - sm(best)

      if ( kn .le. n ) then

        tma(best) = ka(best)
        tmb(best) = kb(best)
        b = basis(kn)
        sm(b) = - sm(b)
        dminus(b) = top
        dplus(b) = d

        call pruef1 ( n, nb, kost, index, top, basis, mem, ka, kb,
     &    sm, y1, y2, dplus, dminus, b )

      else

        kn1 = ka(best)
        b = basis(kn1)
        sm(best) = sm(b)
        dminus(best) = top
        dplus(best) = d
        sm(b) = kn1
        y1(b) = y1(b) + d - dplus(b)
        kn4 = b

        call schr ( n, basis, mem, y1, y2, best, kn4 )

        kb(b) = kn4
        mem(best) = b
        mem(kn4) = best

        call pruef1 ( n, nb, kost, index, top, basis, mem, ka, kb,
     &    sm, y1, y2, dplus, dminus, best )

      end if

      go to 2001
c
c  Shrinking of a blossom
c
4000  continue

      y = y1(base) + d - dplus(base)
      y1(base) = 0.0E+00
      kn4 = base

4010  continue

      y2(kn4) = y2(kn4) + y
      kn4 = mem(kn4)

      if ( kn4 .ne. base ) then
        go to 4010
      end if

      kn3 = mem(base)

      if ( base .ne. b ) then
        go to 4020
      end if

4030  continue

      b = best
      b2 = tma(base)

4020  continue

      mem(kn4) = b2
      kn = - sm(b2)
      sm(b2) = kn
      y1(b2) = y1(b2) + dminus(b2) - d
      kn4 = b2

      call schr ( n, basis, mem, y1, y2, base, kn4 )

      kb(b2) = kn4
      b1 = basis(kn)
      mem(kn4) = b1
      y1(b1) = y1(b1) + d - dplus(b1)
      kn4 = b1

      call schr ( n, basis, mem, y1, y2, base, kn4 )

      kb(b1) = kn4

      if ( b .ne. b1 ) then
        b2 = tma(b1)
        tma(b1) = tmb(b2)
        tmb(b1) = tma(b2)
        go to 4020
      end if

      if ( b .ne. best ) then
        tma(b) = kn2
        tmb(b) = kn1
        if ( base .ne. best ) then
          go to 4030
        end if
      else
        tma(best) = kn1
        tmb(best) = kn2
      end if

      mem(kn4) = kn3
      b = mem(base)
      ka(b) = kn3
      dplus(b) = y
      tma(base) = 0
      dminus(base) = top
      dplus(base) = d

      call pruef1 ( n, nb, kost, index, top, basis, mem, ka, kb,
     &  sm, y1, y2, dplus, dminus, base )

      go to 2001
c
c  Expanding a blossom
c
5000  continue

      kn3 = ka(b)
      b1 = b
      call erw ( n, basis, mem, kb, y1, y2, b1, kn3 )
      y = dplus(b)
      y1(best) = y
      mem(best) = kn3

5010  continue

      y2(kn3) = y2(kn3) - y

      if ( kn3 .ne. best ) then
        kn3 = mem(kn3)
        go to 5010
      end if

      kn1 = - sm(best)
      b1 = basis(kn1)
      kn2 = sm(b1)
      base = basis(kn2)

      if ( base .eq. best ) then
        go to 5100
      end if

      b2 = base

5030  continue

      kn = tma(b2)
      b1 = basis(kn)

      if ( b1 .ne. best ) then
        kn = sm(b1)
        b2 = basis(kn)
        go to 5030
      end if

      tma(base) = tma(best)
      tma(best) = tmb(b2)
      tmb(base) = tmb(best)
      tmb(best) = kn
      kn3 = sm(base)
      b = basis(kn3)
      kn4 = sm(b)
      sm(base) = - kn1
      b1 = b

5050  continue

      kn1 = tma(b1)
      kn2 = tmb(b1)
      tma(b1) = kn4
      tmb(b1) = kn3
      sm(b1) = kn1
      b2 = basis(kn1)
      kn3 = sm(b2)
      sm(b2) = kn2

      if ( b2 .ne. best ) then
        b1 = basis(kn3)
        kn4 = sm(b1)
        tma(b2) = kn3
        tmb(b2) = kn4
      end if

5100  continue

      kn2 = tmb(base)
      b1 = basis(kn2)
      dminus(b1) = d

      if ( b1 .eq. base ) then
        go to 5200
      end if

      kn1 = tma(b1)
      b = basis(kn1)
      tma(b1) = tma(base)
      tmb(b1) = kn2

5110  continue

      kn = sm(b1)
      sm(b1) = - kn
      b2 = basis(kn)
      kn = tma(b2)
      tma(b2) = - kn
      dminus(b2) = top
      dplus(b2) = d
      b1 = basis(kn)
      dminus(b1) = d

      call pruef1 ( n, nb, kost, index, top, basis, mem, ka, kb,
     &  sm, y1, y2, dplus, dminus, b2 )

      if ( b1 .ne. base ) then
        go to 5110
      end if

      tma(base) = tmb(b2)
      tmb(base) = kn

      if ( b .eq. base ) then
        go to 2001
      end if

5200  continue

      b2 = b

5210  continue

      kn = sm(b2)
      sm(b2) = - kn
      b1 = basis(kn)
      tma(b2) = - b1
      kn = tma(b1)
      sm(b1) = - sm(b1)
      b2 = basis(kn)
      tma(b1) = - b2

      if ( b2 .ne. base ) then
        go to 5210
      end if

5220  continue

      b1 = - tma(b)

      call pruef2 ( n, nb, kost, index, top, basis, mem, ka, kb,
     &  sm, y1, y2, dplus, dminus, b, b )

      b = - tma(b1)

      call pruef2 ( n, nb, kost, index, top, basis, mem, ka, kb,
     &  sm, y1, y2, dplus, dminus, b1, b1 )

      if ( b .ne. base ) then
        go to 5220
      end if

      go to 2001
c
c  Modification of a blossom.
c
6000  continue

      dminus(best) = top
      dplus(best) = d
      i = 1
      y1(best) = 0.0E+00
      kn = - sm(best)
      b = basis(kn)
      kn1 = sm(b)

      if ( kn1 .ne. best ) then
        go to 6040
      end if

      i = 2
      sm(b) = kn
      kn3 = mem(best)
      mem(best) = b
      y1(b) = y1(b) + d - dplus(b)
      kn4 = b

      call schr ( n, basis, mem, y1, y2, best, kn4 )

      kb(b) = kn4
      mem(kn4) = kn3
      kn1 = tmb(best)

      if ( kn1 .ne. best ) then
        go to 6040
      end if

6110  continue

      kn = tma(best)
      b = basis(kn)
      sm(best) = sm(b)
      sm(b) = kn
      tma(best) = 0
      kn3 = mem(best)
      mem(best) = b
      y1(b) = y1(b) + d - dplus(b)
      kn4 = b

      call schr ( n, basis, mem, y1, y2, best, kn4 )

      kb(b) = kn4
      mem(kn4) = kn3

      call pruef1 ( n, nb, kost, index, top, basis, mem, ka, kb,
     &  sm, y1, y2, dplus, dminus, best )

      go to 2001

6040  continue

      kn2 = best
      b1 = mem(best)

6080  continue

      kn3 = b1
      kn4 = kb(b1)

6070  continue

      if ( kn3 .eq. kn1 ) then
        go to 6050
      end if

      if ( kn3 .ne. kn4 ) then
        kn3 = mem(kn3)
        go to 6070
      end if

6060  continue

      b1 = mem(kn4)
      kn2 = kn4
      go to 6080

6050  continue

      kn3 = mem(kn4)
      mem(kn2) = kn3

      call erw ( n, basis, mem, kb, y1, y2, b1, kn3 )

      dminus(b1) = d

      if ( i .eq. 2 ) then
        go to 6170
      end if

      i = 2
      tma(b1) = best
      tmb(b1) = sm(b1)
      sm(b1) = - kn
      kn1 = tmb(best)

      if ( kn1 .eq. best ) then
        go to 6110
      end if

      if ( basis(kn1) .eq. best ) then
        go to 6040
      end if

      tma(b1) = tma(best)
      tmb(b1) = kn1
      tma(best) = 0
      b1 = mem(best)

      if ( b1 .eq. best ) then

        sm(best) = nm
        bbest = best

        call pruef2 ( n, nb, kost, index, top, basis, mem, ka, kb,
     &    sm, y1, y2, dplus, dminus, bbest, best )

        go to 2001

      end if

      kn4 = kb(b1)
      kn3 = mem(kn4)
      mem(best) = kn3

      call erw ( n, basis, mem, kb, y1, y2, b1, kn3 )

      sm(best) = - sm(b1)
      sm(b1) = - best
      bbb1 = b1

      call pruef2 ( n, nb, kost, index, top, basis, mem, ka, kb,
     &  sm, y1, y2, dplus, dminus, bbb1, b1 )

      bbest = best

      call pruef2 ( n, nb, kost, index, top, basis, mem, ka, kb,
     &  sm, y1, y2, dplus, dminus, bbest, best )

      go to 2001

6170  continue

      tma(b1) = tma(best)
      tmb(b1) = kn1
      tma(best) = 0
      sm(best) = sm(b1)
      sm(b1) = - best

      call pruef1 ( n, nb, kost, index, top, basis, mem, ka, kb,
     &  sm, y1, y2, dplus, dminus, best )

      go to 2001
c
c  Augmentation.
c
7000  continue

      ii = 0

      do bb = 1, n

        if ( basis(bb) .eq. bb ) then

          kn3 = sm(bb)

          if ( 0 .le. kn3 ) then

            if ( kn3 .eq. 0 ) then
              ii = ii + 1
            end if

            d1 = d - dplus(bb)
            dplus(bb) = 0.0E+00
            y1(bb) = y1(bb) + d1
            sm(bb) = - kn3

          else

            kn3 = tma(bb)

            if ( 0 .lt. kn3 ) then
              d1 = dminus(bb) - d
              y1(bb) = y1(bb) + d1
              tma(bb) = - kn3
            end if

          end if

        end if

      end do

7030  continue

      if ( b1 .ne. b ) then
        b2 = tma(b1)
        tma(b1) = 0
        kn3 = - tma(b2)
        kn4 = tmb(b2)
        sm(b1) = - kn4
        kn = - sm(b2)
        sm(b2) = - kn3
        b1 = basis(kn)
        go to 7030
      end if

      if ( b .ne. best ) then
        tma(b) = 0
        sm(b) = - kn2
        sm(best) = - kn1
        b = best
        b1 = base
        go to 7030
      end if

      tma(best) = 0
      kn = 1

      if ( 2 .lt. ii ) then

        call pruef2 ( n, nb, kost, index, top, basis, mem, ka, kb,
     &    sm, y1, y2, dplus, dminus, kn, n )

        go to 2001

      end if
c
c  Generation of the original graph by expanding all pseudonodes.
c
8000  continue

      kst = 0

      do bb = 1, n

        if ( basis(bb) .eq. bb ) then

          kn1 = - sm(bb)

          if ( kn1 .ne. np ) then

            if ( 0 .lt. kn1 ) then

              b = basis(kn1)
              kn2 = - sm(b)

              call kasu ( nb, kost, index, top, kn1, kn2, kant )

              d = - float ( kost(kant) )
              d = d + y1(bb) + y1(b)
              d = d + y2(kn1) + y2(kn2)

              if ( eps .lt. abs ( d ) ) then
                go to 9995
              end if

              kst = kst + kost(kant)
              sm(b) = kn1
              sm(bb) = kn2

            end if
          end if
        end if

      end do

      do bb = 1, n

8040    continue

        if ( mem(bb) .ne. bb ) then

          base = basis(bb)
          b = mem(base)
          kn1 = tma(b)

          if ( 0 .lt. kn1 ) then
            go to 8050
          end if

          kn2 = sm(base)
          mem(base) = base
          y = y2(base)
          y1(base) = 0.0E+00
          y2(base) = 0.0E+00

8060      continue

          kn4 = kb(b)
          kn3 = mem(kn4)

          call erw ( n, basis, mem, kb, y1, y2, b, kn3 )

          b2 = basis(kn2)
          if ( b2 .eq. b ) then
            go to 8070
          end if

          kn1 = sm(b)

          call kasu ( nb, kost, index, top, base, kn1, kant )

          d = - float ( kost(kant) )
          kst = kst + kost(kant)
          d = d + y2(kn1) + y1(b) + y

          if ( eps .lt. abs ( d ) ) then
            go to 9995
          end if

          go to 8080

8070      continue

          sm(b) = kn2

8080      continue

          b = kn3

          if ( b .ne. base ) then
            go to 8060
          end if

          go to 8040

8050      continue

          kn3 = ka(b)
          b1 = b

          call erw ( n, basis, mem, kb, y1, y2, b1, kn3 )

          mem(base) = kn3
          y = dplus(b)
          y1(base) = y

8105      continue

          y2(kn3) = y2(kn3) - y

          if ( kn3 .ne. base ) then
            kn3 = mem(kn3)
            go to 8105
          end if

          kn1 = sm(base)
          b1 = basis(kn1)

          if ( b1 .eq. base ) then
            go to 8110
          end if

          b = tma(b1)
          b = basis(b)
          kn3 = sm(b1)
          sm(b1) = kn1

8120      continue

          b2 = basis(kn3)
          kn1 = tma(b2)
          kn2 = tmb(b2)
          b1 = basis(kn1)

          call kasu ( nb, kost, index, top, kn1, kn2, kant )

          d = - float ( kost(kant) )
          kst = kst + kost(kant)
          d = d + y1(b1) + y1(b2)
          d = d + y2(kn1) + y2(kn2)

          if ( eps .lt. abs ( d ) ) then
            go to 9995
          end if

          sm(b2) = kn2
          kn3 = sm(b1)
          sm(b1) = kn1

          if ( b1 .ne. base ) then
            go to 8120
          end if

8130      continue

          if ( b .eq. base ) then
            go to 8040
          end if

8110      continue

          kn3 = sm(b)
          b1 = basis(kn3)
          kn4 = sm(b1)

          call kasu ( nb, kost, index, top, kn3, kn4, kant )
 
          kst = kst + kost(kant)
          d = - float ( kost(kant) )
          d = d + y1(b) + y1(b1)
          d = d + y2(kn3) + y2(kn4)

          if ( eps .lt. abs ( d ) ) then
            go to 9995
          end if
        
          sm(b) = kn4
          sm(b1) = kn3
          kn2 = tma(b1)
          b = basis(kn2)
          go to 8130

        end if

      end do
c
c  Print the list of duplicated edges.
c
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  List of duplicated edges.'
      write ( *, '(a)' ) ' '
      i = index(2)

      do kn = 2, n

        ii = index(kn+1) - 1

        do kn2 = i, ii

          kn3 = nb(kn2)

          if ( kn3 .le. 0 ) then

            kn3 = - kn3

            if ( kn3 .le. kn ) then
              write ( *, '(a,i4,a,i4,a)' ) '  (', kn, ',', kn3, ')'
            end if

          end if

        end do

        i = ii + 1

      end do
c
c  Constructing the next node list of the optimal postman tour
c  starting at node KURS.
c
      call tour ( n, nb, kost, index, top, kb, basis, kurs )
      return

9995  continue

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'CPP - Fatal error!'
      write ( *, '(a)' ) '  The optimality conditions are violated.'
      stop
      end
      subroutine dreier ( n, a, b, olwert, perm, ai, bool )

c*********************************************************************72
c
cc DREIER carries out the triple exchange routine for QAPH2.
c
c  Modified:
c
c    28 November 2007
c
c  Author:
c
c    Ulrich Derigs
c
c  Reference:
c
c    Rainer Burkard, Ulrich Derigs,
c    Assignment and Matching Problems: Solution Methods with Fortran-Programs,
c    Springer, 1980,
c    ISBN: 3-540-10267-1,
c    LC: QA402.5 B86.
c
c  Parameters:
c
c    Input, integer N, dimension of the problem.
c
c    Input, integer A(N,N), the distance matrix.
c
c    Input, integer B(N,N), the connection matrix.
c
c    Input/output, integer OLWERT, the objective function value
c    corresponding to the permutation.
c
c    Input/output, integer PERM(N), the permutation.
c
c    Input, integer AI(N), the start sequence for the pairwise
c    exchange algorithm.
c
c    Workspace, logical BOOL(N).
c
      implicit none

      integer n

      integer a(n,n)
      integer a1
      integer a2
      integer a3
      integer adif1
      integer adif2
      integer adif3
      integer ai(n)
      integer as1
      integer as2
      integer as3
      integer asp1
      integer asp2
      integer asp3
      integer b(n,n)
      integer b1
      integer b2
      integer b3
      logical bool(n)
      integer bsp1
      integer bsp2
      integer bsp3
      integer delta
      integer delta1
      integer i
      integer i0
      integer i1
      integer i2
      integer i3
      integer it
      integer izaehl
      integer j
      integer j1
      integer j2
      integer j3
      integer k
      integer l
      integer l1
      integer olwert
      integer perm(n)

      do i = 1, n
        bool(i) = .false.
      end do

      it = 0

1310  continue

      it = it + 1
      izaehl = 0

      do i = 1, n - 2

        i1 = ai(i)
        j1 = perm(i1)
        bool(i1) = .true.

        do j = i + 1, n - 1

          i2 = ai(j)
          j2 = perm(i2)
          bool(i2) = .true.

          do k = j + 1, n

            i3 = ai(k)
            j3 = perm(i3)
            bool(i3) = .true.
!
!  Candidates I1, I2, I3 selected.
!
            delta = 0
            delta1 = 0

            do l = 1, n

              if ( .not. bool(l) ) then

                l1 = perm(l)

                a1 = a(l,i1)
                a2 = a(l,i2)
                a3 = a(l,i3)

                as1 = a(i1,l)
                as2 = a(i2,l)
                as3 = a(i3,l)

                adif1 = a3 - a1
                adif2 = a1 - a2
                adif3 = a2 - a3

                asp1 = as3 - as1
                asp2 = as1 - as2
                asp3 = as2 - as3

                b1 = b(l1,j1)
                b2 = b(l1,j2)
                b3 = b(l1,j3)

                bsp1 = b(j1,l1)
                bsp2 = b(j2,l1)
                bsp3 = b(j3,l1)

                delta = delta - adif1 * b1 - asp1 * bsp1
     &                        - adif2 * b2 - asp2 * bsp2
     &                        - adif3 * b3 - asp3 * bsp3

                delta1 = delta1 + adif1 * b3 + asp1 * bsp3
     &                          + adif2 * b1 + asp2 * bsp1
     &                          + adif3 * b2 + asp3 * bsp2

              end if

            end do

            a1 = a(i2,i3)
            a2 = a(i3,i1)
            a3 = a(i1,i2)

            as1 = a(i3,i2)
            as2 = a(i1,i3)
            as3 = a(i2,i1)

            adif1 = a1 - a2
            adif2 = a2 - a3
            adif3 = a3 - a1

            asp1 = as1 - as2
            asp2 = as2 - as3
            asp3 = as3 - as1

            b1 = b(j3,j1)
            b2 = b(j1,j2)
            b3 = b(j2,j3)

            bsp1 = b(j1,j3)
            bsp2 = b(j2,j1)
            bsp3 = b(j3,j2)

            delta = delta - adif1 * b1 - asp1 * bsp1
     &                    - adif2 * b2 - asp2 * bsp2
     &                    - adif3 * b3 - asp3 * bsp3

            delta1 = delta1 + adif1 * b3 + asp1 * bsp3
     &                      + adif2 * b1 + asp2 * bsp1
     &                      + adif3 * b2 + asp3 * bsp2

            a1 = a(i1,i1)
            a2 = a(i2,i2)
            a3 = a(i3,i3)

            adif1 = a3 - a1
            adif2 = a1 - a2
            adif3 = a2 - a3

            b1 = b(j1,j1)
            b2 = b(j2,j2)
            b3 = b(j3,j3)

            delta = delta - adif1 * b1
     &                    - adif2 * b2
     &                    - adif3 * b3

            delta1 = delta1 + adif1 * b3
     &                      + adif2 * b1
     &                      + adif3 * b2
!
!  If either DELTA is positive, we want to take the exchange.
!
            if ( 0 .lt. delta .or. 0 .lt. delta1 ) then

              if ( delta1 .le. delta ) then

                i0 = j1
                j1 = j2
                j2 = j3
                j3 = i0

                olwert = olwert - delta

              else

                i0 = j1
                j1 = j3
                j3 = j2
                j2 = i0

                olwert = olwert - delta1

              end if

              perm(i1) = j1
              perm(i2) = j2
              perm(i3) = j3

            end if
!
!  Release candidates I3, I2, I1.
!
            bool(i3) = .false.

          end do

          bool(i2) = .false.

        end do

        bool(i1) = .false.

      end do

      if ( 0 .lt. izaehl ) then
        go to 1310
      end if

      return
      end
      subroutine erw ( n, basis, mem, kb, y1, y2, bb, kkn3 )

c*********************************************************************72
c
cc ERW expands a blossom for CPP.
c
c  Modified:
c
c    06 December 2007
c
c  Author:
c
c    W Puetz
c
c  Reference:
c
c    Rainer Burkard, Ulrich Derigs,
c    Assignment and Matching Problems: Solution Methods with Fortran-Programs,
c    Springer, 1980,
c    ISBN: 3-540-10267-1,
c    LC: QA402.5 B86.
c
c  Parameters:
c
c    Unused, integer N, ?
c
c    Output, integer BASIS(*), ?
c
c    Input, integer MEM(*), ?
c
c    Input, integer KB(*), ?
c
c    Input, real Y1(*), ?
c
c    Input/output, real Y2(*), ?
c
c    Input, integer BB, ?
c
c    Input, integer KKN3, ?
c
      implicit none

      integer basis(*)
      integer bb
      integer kb(*)
      integer kkn2
      integer kkn3
      integer kkn4
      integer mem(*)
      integer n
      real y1(*)
      real y2(*)
      real yy

      kkn2 = bb

9800  continue

      bb = kkn2
      kkn4 = kb(bb)
      yy = y1(bb)

9810  continue

      basis(kkn2) = bb
      y2(kkn2) = y2(kkn2) - yy

      if ( kkn2 .ne. kkn4 ) then
        kkn2 = mem(kkn2)
        go to 9810
      end if

      kkn2 = mem(kkn4)
      mem(kkn4) = bb

      if ( kkn2 .ne. kkn3 ) then
        go to 9800
      end if

      return
      end
      subroutine heider ( n, a, b, olwert, perm, ai )

c*********************************************************************72
c
cc HEIDER examines the pairwise interchanges for QAPH2.
c
c  Modified:
c
c    28 November 2007
c
c  Author:
c
c    Ulrich Derigs
c
c  Reference:
c
c    Rainer Burkard, Ulrich Derigs,
c    Assignment and Matching Problems: Solution Methods with Fortran-Programs,
c    Springer, 1980,
c    ISBN: 3-540-10267-1,
c    LC: QA402.5 B86.
c
c  Parameters:
c
c    Input, integer N, dimension of the problem.
c
c    Input, integer A(N,N), the distance matrix.
c
c    Input, integer B(N,N), the connection matrix.
c
c    Input/output, integer OLWERT, the objective function value
c    corresponding to the permutation.
c
c    Input/output, integer PERM(N), the permutation.
c
c    Input, integer AI(N), the start sequence for the pairwise
c    exchange algorithm.
c
      implicit none

      integer n

      integer a(n,n)
      integer ai(n)
      integer b(n,n)
      integer delta
      integer i
      integer i0
      integer i1
      integer i2
      integer ibild
      integer it
      integer izaehl
      integer j
      integer j1
      integer j2
      integer jbild
      integer olwert
      integer perm(n)

      it = 0

1210  continue

      it = it + 1
      izaehl = 0

      do i = 1, n - 1

        i1 = ai(i)
        ibild = perm(i1)

        do j = i + 1, n

          i2 = ai(j)
          jbild = perm(i2)
          delta = 0

          do j1 = 1, n

            if ( j1 .ne. i1 .and. j1 .ne. i2 ) then

              j2 = perm(j1)

              delta = delta + 
     &          ( a(i1,j1) - a(i2,j1) ) * ( b(ibild,j2) - b(jbild,j2) )
     &        + ( a(j1,i1) - a(j1,i2) ) * ( b(j2,ibild) - b(j2,jbild) )

            end if

          end do

          delta = delta 
     &      + ( a(i1,i1) - a(i2,i2) ) 
     &      * ( b(ibild,ibild) - b(jbild,jbild) )
     &      + ( a(i1,i2) - a(i2,i1) )
     &      * ( b(ibild,jbild) - b(jbild,ibild) )
!
!  Make a swap if it will decrease the objective function.
!
          if ( 0 .lt. delta ) then

            olwert = olwert - delta
            i0 = ibild
            ibild = jbild
            jbild = i0

            perm(i1) = ibild
            perm(i2) = jbild

            izaehl = izaehl + 1

          end if

        end do

      end do
!
!  If at least one switch was made, it may be worth it to try again.
!
      if ( izaehl .gt. 0 ) then
        go to 1210
      end if

      return
      end
      subroutine kasu ( nb, kost, index, top, kn1, kn2, kant )

c*********************************************************************72
c
cc KASU duplicates matching edges for CPP.
c
c  Modified:
c
c    06 December 2007
c
c  Author:
c
c    W Puetz
c
c  Reference:
c
c    Rainer Burkard, Ulrich Derigs,
c    Assignment and Matching Problems: Solution Methods with Fortran-Programs,
c    Springer, 1980,
c    ISBN: 3-540-10267-1,
c    LC: QA402.5 B86.
c
c  Parameters:
c
c    Input, integer NB(*), ?
c
c    Unused, integer KOST(*), ?
c
c    Input, integer INDEX(*), ?
c
c    Unused, integer TOP, ?
c
c    Input, integer KN1, KN2, ?
c
c    Output, integer KANT, ?
c
      implicit none

      integer index(*)
      integer k1
      integer k2
      integer kant
      integer kn1
      integer kn2
      integer kost(*)
      integer nb(*)
      integer top

      k1 = kn1
      k2 = kn2

9100  continue

      kant = index(k1)

9120  continue

      if ( nb(kant) .ne. k2 ) then
        kant = kant + 1
        go to 9120
      end if

      nb(kant) = - k2

      if ( k1 .ne. kn2 ) then
        k1 = kn2
        k2 = kn1
        go to 9100
      end if

      return
      end
      subroutine lbap ( n, sup, c, z, zeile, spalte, dminus, dplus, 
     &  vor, vos, label )

c*********************************************************************72
c
cc LBAP solves the linear bottleneck assignment problem.
c
c  Modified:
c
c    19 November 2007
c
c  Author:
c
c    Ulrich Derigs
c
c  Reference:
c
c    Rainer Burkard, Ulrich Derigs,
c    Assignment and Matching Problems: Solution Methods with Fortran-Programs,
c    Springer, 1980,
c    ISBN: 3-540-10267-1,
c    LC: QA402.5 B86.
c
c  Parameters:
c
c    Input, integer N, the dimension of the cost matrix.
c
c    Input, integer SUP, a large machine number.
c
c    Input, integer C(N,N), the cost matrix.
c
c    Output, integer Z, the optimal value.
c
c    Workspace, integer ZEILE(N).
c
c    Output, integer SPALTE(N), the optimal assignment.
c
c    Workspace, integer DMINUS(N), DPLUS(N), VOR(N), VOS(N).
c
c    Workspace, logical LABEL(N).
c
      implicit none

      integer n

c     integer c(n,n)
      integer c(n*n)
      integer cij
      integer d
      integer dminus(n)
      integer dplus(n)
      integer i
      integer ii
      integer ind
      integer index
      integer is
      integer j
      integer jj
      integer k
      integer kj
      logical label(n)
      integer nfrei
      integer nfrei1
      integer s
      integer s1
      integer sj
      integer spalte(n)
      integer sup
      integer u
      integer u1
      integer us
      integer usi
      integer vgl
      integer vor(n)
      integer vos(n)
      integer w
      integer ws
      integer wsi
      integer z
      integer zeile(n)
      integer zmin
c
c  Threshhold totals method.
c
      do i = 1, n
        zeile(i) = 0
        dplus(i) = 0
        vos(i) = 0
        spalte(i) = 0
        label(i) = .false.
      end do
c
c  Determination of a lower bound.
c
      z = c(1)
      do j = 2, n
        if ( c(j) .lt. z ) then
          z = c(j)
        end if
      end do

      s = n
      do i = 2, n
        s1 = s + 1
        zmin = c(s+1)
        kj = 1
        if ( z .lt. zmin ) then

          do j = 2, n

            cij = c(s+j)
            if ( cij .le. z ) then
              label(j) = .true.
              go to 50
            end if

            if ( cij .lt. zmin ) then
              zmin = cij
              kj = j
            end if

          end do

          z = zmin

        end if

        label(kj) = .true.

50      continue

        s = s + n

      end do

      do j = 1, n

        if ( .not. label(j) ) then

          zmin = c(j)

          if ( z .lt. zmin ) then

            s = n

            do i = 2, n
              cij = c(s+j)
              if ( cij .le. z ) then
                go to 60
              end if
              if ( cij .lt. zmin ) then
                zmin = cij
              end if
              s = s + n
            end do

            z = zmin

60          continue

          end if

        end if

      end do
c
c  Sorting the rows and columns.
c
      s = 0

      do i = 1, n

        k = 0

        do j = 1, n

          if ( c(s+j) .le. z ) then
            k = k + 1
            dplus(j) = dplus(j) + 1
          end if

        end do

        dminus(i) = k
        vor(i) = i
        s = s + n

      end do

      call ssort ( dminus, vor, n )
      call ssort ( dplus, vos, n )
c
c  Construction of an initial partial assignment.
c
      do ii = 1, n

        i = vor(ii)
        is = ( i - 1 ) * n

        if ( spalte(i) .eq. 0 ) then

          do jj = 1, n
            j = vos(jj)
            if ( zeile(j) .ne. 0 ) go to 95
            if ( c(is+j) .gt. z ) go to 95
            zeile(j) = i
            spalte(i) = j
            go to 96
95          continue
          end do

        end if

96      continue

      end do

      nfrei = 0
      nfrei1 = n + 1

      do j = 1, n
        if ( zeile(j) .le. 0 ) then
          nfrei = nfrei + 1
          vos(nfrei) = j
        else
          nfrei1 = nfrei1 - 1
          vos(nfrei1) = j
        end if
      end do
c
c  Construction of the optimal assignment.
c
      do 1000 u = 1, n

        if ( spalte(u) .gt. 0 ) go to 1000
c
c  Shortest path computation.
c
        u1 = u - 1
        us = ( u - 1 ) * n

        do i = 1, n
          vor(i) = u
          label(i) = .false.
          dplus(i) = sup
          usi = us + i
          dminus(i) = c(usi)
        end do

        dplus(u) = 0

105     continue

        d = sup

        do ii = 1, nfrei

          i = vos(ii)

          if ( dminus(i) .lt. d ) then
            d = dminus(i)
            index = i
            ind = ii
            if ( d .le. z ) then
              go to 399
            end if
          end if

        end do

        do ii = nfrei1, n

          i = vos(ii)

          if ( .not. label(i) ) then

            if ( dminus(i) .lt. d ) then

              d = dminus(i)
              index = i
              ind = ii

              if ( d .le. z ) then
                go to 112
              end if

            end if

         end if

        end do

111     continue

        if ( zeile(index) .le. 0 ) then
          go to 399
        end if

112     continue

        label(index) = .true.
        w = zeile(index)
        ws = ( w - 1 ) * n
        dplus(w) = d

        do i = 1, n

          if ( .not. label(i) ) then

            vgl = c(ws+i)

            if ( vgl .lt. d ) then
              vgl = d
            end if

            if ( vgl .lt. dminus(i) ) then
              dminus(i) = vgl
              vor(i) = w
            end if

          end if

        end do

        go to 105
c
c  Augmentation
c
399     continue

        vos(ind) = vos(nfrei)
        vos(nfrei) = index
        nfrei1 = nfrei
        nfrei = nfrei - 1

400     continue

        w = vor(index)
        zeile(index) = w
        ind = spalte(w)
        spalte(w) = index

        if ( w .eq. u ) then
          go to 500
        end if

        index = ind

        go to 400
c
c  Update the lower bound.
c
500     continue

        if ( z .lt. d ) then
          z = d
        end if

1000  continue

      return
      end
      subroutine lsapi ( n, sup, c, z, zeile, spalte, dminus, dplus, 
     &  ys, yt, vor, label )

c*********************************************************************72
c
cc LSAPI solves the linear sum assignment problem with integer data.
c
c  Modified:
c
c    20 November 2007
c
c  Author:
c
c    Ulrich Derigs
c
c  Reference:
c
c    Rainer Burkard, Ulrich Derigs,
c    Assignment and Matching Problems: Solution Methods with Fortran-Programs,
c    Springer, 1980,
c    ISBN: 3-540-10267-1,
c    LC: QA402.5 B86.
c
c  Parameters:
c
c    Input, integer N, the dimension of the cost matrix.
c
c    Input, integer SUP, a large machine number.
c
c    Input, integer C(N,N), the cost matrix.
c
c    Output, integer Z, the optimal value.
c
c    Workspace, integer ZEILE(N).
c
c    Output, integer SPALTE(N), the optimal assignment.
c
c    Workspace, integer DMINUS(N), DPLUS(N), YS(N), YT(N), VOR(N).
c
c    Workspace, logical LABEL(N).
c
      implicit none

      integer n

c     integer c(n,n)
      integer c(n*n)
      integer cc
      integer d
      integer dminus(n)
      integer dplus(n)
      integer i
      integer ik
      integer ind
      integer index
      integer is
      integer isj
      integer j
      integer j0
      logical label(n)
      integer spalte(n)
      integer sup
      integer u
      integer ui
      integer us
      integer usi
      integer vgl
      integer vj
      integer vor(n)
      integer w
      integer ws
      integer wsi
      integer ys(n)
      integer yt(n)
      integer z
      integer zeile(n)
c
c  Construct an initial partial assignment.
c

      do i = 1, n
        zeile(i) = 0
        spalte(i) = 0
        vor(i) = 0
        ys(i) = 0
        yt(i) = 0
      end do

      ik = 0

      do i = 1, n

        do j = 1, n

          ik = ik + 1
          cc = c(ik)

          if ( j .eq. 1 ) then
            go to 4
          end if

          if ( ( cc - ui ) .ge. 0 ) then
            go to 3
          end if

4         continue

          ui = cc
          j0 = j

3         continue

        end do

        ys(i) = ui

        if ( zeile(j0) .eq. 0 ) then
          zeile(j0) = i
          spalte(i) = j0
        end if

      end do

      do j = 1, n
        yt(j) = 0
        if ( zeile(j) .eq. 0 ) then
          yt(j) = sup
        end if
      end do

      ik = 0

      do i = 1, n

        ui = ys(i)

        do j = 1, n

          ik = ik + 1
          vj = yt(j)

          if ( 0 .lt. vj ) then

            cc = c(ik) - ui

            if ( cc .lt. vj ) then
              yt(j) = cc
              vor(j) = i
            end if

          end if

        end do

      end do

      do j = 1, n

        i = vor(j)

        if ( i .ne. 0 ) then

          if ( spalte(i) .eq. 0 ) then
            spalte(i) = j
            zeile(j) = i
          end if

        end if

      end do

      do i = 1, n

        if ( spalte(i) .eq. 0 ) then

          ui = ys(i)
          ik = ( i - 1 ) * n

          do j = 1, n

            ik = ik + 1

            if ( zeile(j) .eq. 0 ) then

              cc = c(ik)

              if ( ( cc - ui - yt(j) ) .le. 0 ) then
                spalte(i) = j
                zeile(j) = i
              end if

            end if

          end do

        end if

      end do
c
c  Construct the optimal solution.
c
      do 1000 u = 1, n

        if ( spalte(u) .gt. 0 ) then
          go to 1000
        end if
c
c  Shortest path computation.
c
        us = ( u - 1 ) * n

        do i = 1, n
          vor(i) = u
          label(i) = .false.
          dplus(i) = sup
          usi = us + i
          dminus(i) = c(usi) - ys(u) - yt(i)
        end do

        dplus(u) = 0

105     continue

        d = sup
        index = 0

        do i = 1, n

          if ( .not. label(i) ) then
            if ( dminus(i) .lt. d ) then
              d = dminus(i)
              index = i
            end if
          end if

        end do

        if ( index .eq. 0 ) then
          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) 'LSAPI - Fatal error!'
          write ( *, '(a)' ) '  No unlabeled node with DMINUS < D.'
          stop
        end if

        write ( *, * ) 'ZEILE(', index, ' ) = ', zeile(index)
        write ( *, * ) 'Label(', index, ' ) = true'

        if ( zeile(index) .le. 0 ) then
          go to 400
        end if

        label(index) = .true.
        w = zeile(index)
        ws = ( w - 1 ) * n
        dplus(w) = d

        do i = 1, n

          if ( .not. label(i) ) then

            wsi = ws + i
            vgl = d + c(wsi) - ys(w) - yt(i)

            if ( vgl .lt. dminus(i) ) then
              dminus(i) = vgl
              vor(i) = w
            end if

          end if

        end do

        go to 105
c
c  Augmentation.
c
400     continue

        w = vor(index)
        zeile(index) = w
        ind = spalte(w)
        spalte(w) = index

        if ( w .ne. u ) then
          index = ind
          go to 400
        end if
c
c  Transformation.
c
        do i = 1, n

          if ( dplus(i) .ne. sup ) then
            ys(i) = ys(i) + d - dplus(i)
          end if

          if ( dminus(i) .lt. d ) then
            yt(i) = yt(i) + dminus(i) - d
          end if

        end do

1000  continue
c
c  Computation of the optimal value.
c
      z = 0
      do i = 1, n
        is = ( i - 1 ) * n
        j = spalte(i)
        z = z + c(is+j)
      end do

      return
      end
      subroutine lsapr ( n, isup, c, z, zeile, spalte, dminus, dplus, 
     &  ys, yt, vor, label, eps )

c*********************************************************************72
c
cc LSAPR solves the linear sum assignment problem with real data.
c
c  Modified:
c
c    22 November 2007
c
c  Author:
c
c    Ulrich Derigs
c
c  Reference:
c
c    Rainer Burkard, Ulrich Derigs,
c    Assignment and Matching Problems: Solution Methods with Fortran-Programs,
c    Springer, 1980,
c    ISBN: 3-540-10267-1,
c    LC: QA402.5 B86.
c
c  Parameters:
c
c    Input, integer N, the dimension of the cost matrix.
c
c    Input, integer ISUP, a large machine number.
c
c    Input, real C(N,N), the cost matrix.
c
c    Output, real Z, the optimal value.
c
c    Workspace, integer ZEILE(N).
c
c    Output, integer SPALTE(N), the optimal assignment.
c
c    Workspace, real DMINUS(N), DPLUS(N), YS(N), YT(N)
c
c    Workspace, integer VOR(N).
c
c    Workspace, logical LABEL(N).
c
c    Input, real EPS, machine accuracy.
c
      implicit none

      integer n

c     real c(n,n)
      real c(n*n)
      real cc
      integer d
      real dminus(n)
      real dplus(n)
      real eps
      integer i
      integer ik
      integer ind
      integer index
      integer is
      integer isj
      integer isup
      integer j
      integer j0
      logical label(n)
      integer spalte(n)
      real sup
      integer u
      real ui
      integer us
      integer usi
      real vgl
      real vj
      integer vor(n)
      integer w
      integer ws
      integer wsi
      real ys(n)
      real yt(n)
      real z
      integer zeile(n)
c
c  Construct an initial partial assignment.
c
      sup = real ( isup )

      do i = 1, n
        zeile(i) = 0
        spalte(i) = 0
        vor(i) = 0
        ys(i) = 0.0E+00
        yt(i) = 0.0E+00
      end do

      ik = 0

      do i = 1, n

        do j = 1, n

          ik = ik + 1
          cc = c(ik)

          if ( j .ne. 1 ) then

            if ( ( cc - ui ) .ge. eps ) then
              go to 3
            end if

          end if

          ui = cc
          j0 = j

3         continue

        end do

        ys(i) = ui

        if ( zeile(j0) .eq. 0 ) then
          zeile(j0) = i
          spalte(i) = j0
        end if

      end do

      do j = 1, n
        yt(j) = 0
        if ( zeile(j) .eq. 0 ) then
          yt(j) = sup
        end if
      end do

      ik = 0

      do i = 1, n

        ui = ys(i)

        do j = 1, n

          ik = ik + 1
          vj = yt(j)

          if ( eps .lt. vj ) then

            cc = c(ik) - ui

            if ( cc + eps .lt. vj ) then
              yt(j) = cc
              vor(j) = i
            end if

          end if

        end do

      end do

      do j = 1, n

        i = vor(j)

        if ( i .ne. 0 ) then

          if ( spalte(i) .eq. 0 ) then
            spalte(i) = j
            zeile(j) = i
          end if

        end if

      end do

      do i = 1, n

        if ( spalte(i) .eq. 0 ) then

          ui = ys(i)
          ik = ( i - 1 ) * n

          do j = 1, n

            ik = ik + 1

            if ( zeile(j) .eq. 0 ) then

              cc = c(ik)

              if ( ( cc - ui - yt(j) + eps ) .le. 0.0E+00 ) then
                spalte(i) = j
                zeile(j) = i
              end if

            end if

          end do

        end if

      end do
c
c  Construct the optimal solution.
c
      do 1000 u = 1, n

        if ( spalte(u) .gt. 0 ) then
          go to 1000
        end if
c
c  Shortest path computation.
c
        us = ( u - 1 ) * n

        do i = 1, n
          vor(i) = u
          label(i) = .false.
          dplus(i) = sup
          usi = us + i
          dminus(i) = c(usi) - ys(u) - yt(i)
        end do

        dplus(u) = 0.0E+00

105     continue

        d = sup
        index = 0

        do i = 1, n

          if ( .not. label(i) ) then
            if ( dminus(i) + eps .lt. d ) then
              d = dminus(i)
              index = i
            end if
          end if

        end do

        if ( index .eq. 0 ) then
          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) 'LSAPR - Fatal error!'
          write ( *, '(a)' ) '  No unlabeled node with DMINUS < D.'
          stop
        end if

        if ( zeile(index) .le. 0 ) then
          go to 400
        end if

        label(index) = .true.
        w = zeile(index)
        ws = ( w - 1 ) * n
        dplus(w) = d

        do i = 1, n

          if ( .not. label(i) ) then

            wsi = ws + i
            vgl = d + c(wsi) - ys(w) - yt(i)

            if ( vgl + eps .lt. dminus(i) ) then
              dminus(i) = vgl
              vor(i) = w
            end if

          end if

        end do

        go to 105
c
c  Augmentation.
c
400     continue

        w = vor(index)
        zeile(index) = w
        ind = spalte(w)
        spalte(w) = index

        if ( w .ne. u ) then
          index = ind
          go to 400
        end if
c
c  Transformation.
c
!500     continue

        do i = 1, n

          if ( dplus(i) .ne. sup ) then
            ys(i) = ys(i) + d - dplus(i)
          end if

          if ( dminus(i) + eps .lt. d ) then
            yt(i) = yt(i) + dminus(i) - d
          end if

        end do

1000  continue
c
c  Computation of the optimal value.
c
      z = 0.0E+00
      do i = 1, n
        is = ( i - 1 ) * n
        j = spalte(i)
        z = z + c(is+j)
      end do

      return
      end
      subroutine progno ( n, k, nmk, umspei, veksum, vekqu, weiter,
     &  aspei, bspei, cspei )

c*********************************************************************72
c
cc PROGNO is used by QAP to compute the linear subproblem cost matrix.
c
c  Discussion:
c
c    Computation of the final cost matrix UMSPEI for the linear subproblem
c    LSAPI on the K-th stage of the decision tree.
c
c  Modified:
c
c    15 December 2007
c
c  Author:
c
c    T Boenniger,
c    Rainer Burkard,
c    Karl-Heinz Stratmann
c
c  Reference:
c
c    Rainer Burkard, Ulrich Derigs,
c    Assignment and Matching Problems: Solution Methods with Fortran-Programs,
c    Springer, 1980,
c    ISBN: 3-540-10267-1,
c    LC: QA402.5 B86.
c
c  Parameters:
c
c    Input, integer N, ?
c
c    Input, integer K, ?
c
c    Input, integer NMK, ?
c
c    Output, integer UMSPEI(N*N), ?
c
c    Input, integer VEKSUM(N-2), ?
c
c    Input, integer VEKQU(N-2), ?
c
c    Input, logical WEITER, ?
c
c    Input, integer ASPEI(*), ?
c
c    Input, integer BSPEI(*), ?
c
c    Input/output, integer CSPEI(*), ?
c
      implicit none

      integer n

      integer aspei(*)
      integer bspei(*)
      integer cspei(*)
      integer i
      integer iz
      integer j
      integer j1
      integer j2
      integer j2iz
      integer j3
      integer j3iz
      integer ju
      integer k
      integer nmk
      integer nsum2
      integer t
      integer umspei(n*n)
      integer vekqu(n-2)
      integer veksum(n-2)
      logical weiter

      ju = 0

      if ( nmk .ne. n ) then
        go to 3600
      end if

      nsum2 = 0
      j1 = 0

      if ( weiter ) then
        go to 3680
      end if

      go to 3690

3600  continue

      if ( weiter ) then
        go to 3630
      end if

      j1 = vekqu(k)

3690  continue

      do i = 1, nmk
        do j = 1, nmk
          j1 = j1 + 1
          ju = ju + 1
          umspei(ju) = cspei(j1)
        end do
      end do

      return

3630  continue

      nsum2 = veksum(k+1)
      j1 = vekqu(k+1)

3680  continue

      j2 = nsum2

      do i = 1, nmk

        j3 = nsum2

        do j = 1, nmk

          j1 = j1 + 1
          ju = ju + 1

          t = 0
          do iz = 1, nmk - 1
            j2iz = j2 + iz
            j3iz = j3 + iz
            t = t + aspei(j2iz) * bspei(j3iz)
          end do

          j3 = j3 + nmk - 1
          cspei(j1) = t
          umspei(ju) = t

        end do

        j2 = j2 + nmk - 1

      end do
 
      return
      end
      subroutine pruef1 ( n, nb, kost, index, top, basis, mem, ka, kb,
     &  sm, y1, y2, dplus, dminus, bb )

c*********************************************************************72
c
cc PRUEF1 is used by CPP to scan node BB.
c
c  Modified:
c
c    06 December 2007
c
c  Author:
c
c    W Puetz
c
c  Reference:
c
c    Rainer Burkard, Ulrich Derigs,
c    Assignment and Matching Problems: Solution Methods with Fortran-Programs,
c    Springer, 1980,
c    ISBN: 3-540-10267-1,
c    LC: QA402.5 B86.
c
c  Parameters:
c
c    Unused, integer N.
c
c    Input, integer NB(*), ?
c
c    Input, integer KOST(*), ?
c
c    Input, integer INDEX(*), ?
c
c    Unused, integer TOP, ?
c
c    Input, integer BASIS(*), ?
c
c    Input, integer MEM(*), ?
c
c    Output, integer KA(*), ?
c
c    Output, integer KB(*), ?
c
c    Input, integer SM(*), ?
c
c    Input, real Y1(*), ?
c
c    Input, real Y2(*), ?
c
c    Input, real DPLUS(*), ?
c
c    Input/output, real DMINUS(*), ?
c
c    Input, integer BB, the node to scan.
c
      implicit none

      integer basis(*)
      integer bb
      integer bb1
      integer bb2
      real d1
      real d2
      real d3
      real dminus(*)
      real dplus(*)
      integer i
      integer index(*)
      integer ka(*)
      integer kant
      integer kb(*)
      integer kkn
      integer kkn1
      integer kkn2
      integer kost(*)
      integer mem(*)
      integer n
      integer nb(*)
      integer nn
      integer sm(*)
      integer top
      real y
      real y1(*)
      real y2(*)
      real yy

      d1 = dplus(bb) - y1(bb)
      kkn = bb
      kkn1 = sm(bb)
      bb1 = -1

      if ( 0 .lt. kkn1 ) then
        bb1 = basis(kkn1)
      end if

9330  continue

      i = index(kkn)
      nn = index(kkn+1) - 1
      y = y2(kkn)

      do kant = i, nn

        kkn2 = nb(kant)
        bb2 = basis(kkn2)

        if ( bb .ne. bb2 ) then
          if ( bb1 .ne. bb2 ) then

            d2 = dminus(bb2)
            yy = y1(bb2) + y2(kkn2)
            d3 = float ( kost(kant) )
            d3 = d3 - y - yy + d1

            if ( d3 .lt. d2 ) then
              dminus(bb2) = d3
              ka(bb2) = kkn
              kb(bb2) = kkn2
            end if

          end if
        end if

      end do

      kkn = mem(kkn)

      if ( kkn .ne. bb ) then
        go to 9330
      end if

      return
      end
      subroutine pruef2 ( n, nb, kost, index, top, basis, mem, ka, kb,
     &  sm, y1, y2, dplus, dminus, bb2, bbb )

c*********************************************************************72
c
cc PRUEF2 is used by CPP to scan nodes BB2 and BBB.
c
c  Modified:
c
c    06 December 2007
c
c  Author:
c
c    W Puetz
c
c  Reference:
c
c    Rainer Burkard, Ulrich Derigs,
c    Assignment and Matching Problems: Solution Methods with Fortran-Programs,
c    Springer, 1980,
c    ISBN: 3-540-10267-1,
c    LC: QA402.5 B86.
c
c  Parameters:
c
c    Unused, integer N.
c
c    Input, integer NB(*), ?
c
c    Input, integer KOST(*), ?
c
c    Input, integer INDEX(*), ?
c
c    Unused, integer TOP, ?
c
c    Input, integer BASIS(*), ?
c
c    Input, integer MEM(*), ?
c
c    Output, integer KA(*), ?
c
c    Output, integer KB(*), ?
c
c    Input, integer SM(*), ?
c
c    Input, real Y1(*), ?
c
c    Input, real Y2(*), ?
c
c    Input, real DPLUS(*), ?
c
c    Input/output, real DMINUS(*), ?
c
c    Input/output, integer BB2, a node to scan.
c
c    Input, integer BBB, a node to scan.
c
      implicit none

      integer basis(*)
      integer bb
      integer bb2
      integer bbb
      real d1
      real d2
      real d3
      real dminus(*)
      real dplus(*)
      integer i
      integer index(*)
      integer ka(*)
      integer kant
      integer kb(*)
      integer kkn
      integer kkn1
      integer kkn2
      integer kost(*)
      integer mem(*)
      integer n
      integer nb(*)
      integer nn
      integer sm(*)
      integer top
      real y
      real y1(*)
      real y2(*)
      real yy

9000  continue

      kkn2 = basis(bb2)

      if ( kkn2 .eq. bb2 ) then

        d2 = top
        yy = y1(bb2)

9040    continue

        i = index(kkn2)
        nn = index(kkn2+1) - 1
        y = y2(kkn2)

        do kant = i, nn

          kkn = nb(kant)
          bb = basis(kkn)

          if ( bb .ne. bb2 ) then
            if ( 0 .le. sm(bb) ) then

              d1 = dplus(bb) - y1(bb) - y2(kkn)
              d3 = float ( kost(kant) )
              d3 = d3 + d1 - yy - y

              if ( d3 .lt. d2 ) then
                d2 = d3
                ka(bb2) = kkn
                kb(bb2) = kkn2
              end if

            end if
          end if

        end do

        kkn2 = mem(kkn2)

        if ( kkn2 .ne. bb2 ) then
          go to 9040
        end if

        dminus(bb2) = d2

      end if

      bb2 = bb2 + 1

      if ( bb2 .le. bbb ) then
        go to 9000
      end if

      return
      end
      subroutine qap ( n, a, b, unendl, loesg, olwert, c, umspei, zul,
     &  u, v, dd, partpe, y, lab, z1, h1, meng, phiofm, menge, veksum,
     &  vekqu, alter, aspei, bspei, cspei, bool, bool1, hl1 )

c*********************************************************************72
c
cc QAP solves the quadratic assignment problem.
c
c  Discussion:
c
c    This routine, or those it calls, seems to have a bug.  It is not
c    able to compute the correct answer for the sample problem in the
c    reference.
c
c  Modified:
c
c    14 December 2007
c
c  Author:
c
c    T Boenniger,
c    Rainer Burkard,
c    Karl-Heinz Stratmann
c
c  Reference:
c
c    Rainer Burkard, Ulrich Derigs,
c    Assignment and Matching Problems: Solution Methods with Fortran-Programs,
c    Springer, 1980,
c    ISBN: 3-540-10267-1,
c    LC: QA402.5 B86.
c
c  Parameters:
c
c    Input, integer N, the dimension of the problem.
c
c    Input/output, integer A(N,N), the distance matrix.  On input, A
c    should contain nonnegative values.  Entries in A are changed 
c    during execution.
c
c    Input/output, integer B(N,N), the connection matrix.  On input, B
c    should contain nonnegative values.  Entries in B are changed
c    during execution.
c
c    Input, integer UNENDL, a large value.
c
c    Output, integer LOESG(N), the optimal assignment.
c
c    Output, integer OLWERT, the value of the objective function for
c    the optimal assignment.
c
c    Workspace, integer C(N,N).
c
c    Workspace, integer UMSPEI(N,N).
c
c    Workspace, integer ZUL(N,N).
c
c    Workspace, integer U(N).
c
c    Workspace, integer V(N).
c
c    Workspace, integer DD(N).
c
c    Workspace, integer PARTPE(N).
c
c    Workspace, integer Y(N).
c
c    Workspace, integer LAB(N).
c
c    Workspace, integer Z1(N).
c
c    Workspace, integer H1(N).
c
c    Workspace, integer MENG(N-2).
c
c    Workspace, integer PHIOFM(N-2).
c
c    Workspace, integer MENGE(N-2).
c
c    Workspace, integer VEKSUM(N-2).
c
c    Workspace, integer VEKQU(N-2).
c
c    Workspace, integer ALTER(N-2).
c
c    Workspace, integer ASPEI((N*(N+1)*(N-1))/3).
c
c    Workspace, integer BSPEI((N*(N+1)*(N-1))/3).
c
c    Workspace, integer CSPEI((N*(N+1)*(2*N+1))/6-1).
c
c    Workspace, logical BOOL(N).
c
c    Workspace, logical BOOL1(N).
c
c    Workspace, logical HL1(N).
c
      implicit none

      integer n

      integer a(n,n)
      integer alter(n-2)
      integer alterk
      integer am
      integer aspei((n*(n+1)*(n-1))/3)
      integer b(n,n)
      integer bm
      integer bmj
      logical bool(n)
      logical bool1(n)
      integer bspei((n*(n+1)*(n-1))/3)
      integer c(n,n)
      integer ca
      integer ca1
      integer ccc
      integer cspei((n*(n+1)*(2*n+1))/6-1)
      integer dd(n)
      integer h1(n)
      logical hl1(n)
      integer i
      integer i1
      integer i2
      integer ihalt
      integer ikap
      integer iz
      integer ize
      integer izaehl
      integer j
      integer j1
      integer j2
      integer jhalt
      integer jj
      integer jz
      integer jzaehl
      integer k
      integer kk
      integer lab(n)
      integer loesg(n)
      integer meng(n-2)
      integer menge(n-2)
      integer min
      integer min1
      integer min2
      integer nm2
      integer nmk
      integer olwert
      integer partpe(n)
      integer phiofm(n-2)
      integer ra
      integer raa
      integer rb
      integer rbb
      integer t
      integer t1
      integer t2
      integer u(n)
      integer umspei(n*n)
      integer unendl
      integer v(n)
      integer vekqu(n-2)
      integer veksum(n-2)
      logical weiter
      integer y(n)
      integer z1(n)
      integer zpart
      integer zstern
      integer zul(n,n)
c
c  Initialization.
c
      k = 0
      olwert = unendl
      nm2 = n - 2
      i2 = n + 1
      j = 0
      j1 = 0

      do i = 1, n - 2
        nmk = i2 - i
        j2 = nmk * nmk
        j = j + j2 - nmk
        veksum(i) = j
        j1 = j1 + j2
        vekqu(i) = j1
      end do
c
c  1: Reduction, C(I,J) = A(I,I) * B(J,J).
c
      do i = 1, n
        bool(i) = .false.
        bool1(i) = .false.
        zstern = a(i,i)
        do j = 1, n
          zul(i,j) = -1
          c(i,j) = zstern * b(j,j)
        end do
      end do

      do j = 1, n
        a(j,j) = unendl
        b(j,j) = unendl
      end do
c
c  2: Reduction of matrices A and B.
c
      do i = 1, n
        ra = a(i,1)
        rb = b(i,1)
        do j = 2, n
          raa = a(i,j)
          rbb = b(i,j)
          if ( raa .lt. ra ) then
            ra = raa
          end if
          if ( rbb .lt. rb ) then
            rb = rbb
          end if
        end do
        do j = 1, n
          a(i,j) = a(i,j) - ra
          b(i,j) = b(i,j) - rb
        end do
        u(i) = ra
        v(i) = rb
      end do
c
c  Determination of the matrix C.
c  Rowwise reduction of matrices A and B.
c
      do i = 1, n

        ca1 = u(i)

        do j = 1, n

          bmj = v(j)
          bm = ( n - 1 ) * bmj

          do kk = 1, n
            if ( kk .ne. j ) then
              bm = bm + b(j,kk)
            end if
          end do

          ca = ca1 * bm
          am = 0

          do kk = 1, n
            if ( kk .ne. i ) then
              am = am + a(i,kk)
            end if
          end do

          c(i,j) = ca + bmj * am + c(i,j)

        end do

      end do
c
c  Columnwise reduction of matrices A and B.
c
      do i = 1, n
        ra = a(1,i)
        rb = b(1,i)
        do j = 2, n
          raa = a(j,i)
          rbb = b(j,i)
          if ( raa .lt. ra ) then
            ra = raa
          end if
          if ( rbb .lt. rb ) then
            rb = rbb
          end if
        end do
        do j = 1, n
          a(j,i) = a(j,i) - ra
          b(j,i) = b(j,i) - rb
        end do
        u(i) = ra
        v(i) = rb
      end do

      do i = 1, n
        a(i,i) = 0
        b(i,i) = 0
        ca1 = u(i)
        do j = 1, n
          bmj = v(j)
          bm = ( n - 1 ) * bmj
          do kk = 1, n
            if ( kk .ne. j ) then
              bm = bm + b(kk,j)
            end if
          end do
          ca = ca1 * bm
          am = 0
          do kk = 1, n
            if ( kk .ne. i ) then
              am = am + a(kk,i)
            end if
          end do
          ccc = c(i,j) + ca + bmj * am
          c(i,j) = ccc
        end do
      end do

      zpart = 0
      nmk = n
      weiter = .true.
c
c  Improvement of the bounds C(I,J) by adding minimal scalar products
c  (Gilmore bounds).  The routines WEGSPE and PROGNO compute these
c  minimal scalar products.
c
      call wegspe ( n, k, nmk, izaehl, jzaehl, a, b, veksum, u, bool,
     &  bool1, aspei, bspei, h1 )

      call progno ( n, k, nmk, umspei, veksum, vekqu, weiter, aspei,
     &  bspei, cspei )

      weiter = .false.

5210  continue

      iz = 0
      i2 = 0
      j2 = 0
      ikap = 0

      do i = 1, n

        if ( .not. bool(i) ) then

          jj = iz * nmk
          iz = iz + 1
          jz = 0

          do j = 1, n

            if ( .not. bool1(j) ) then

              jz = jz + 1
              jj = jj + 1

              if ( 0 .le. zul(i,j) ) then

                umspei(jj) = unendl
                ikap = ikap + 1

                if ( ikap .lt. 2 ) then
                  ihalt = iz
                  jhalt = jz
                else if ( ikap .eq. 2 ) then
                  if ( iz .eq. ihalt ) then
                    i2 = iz
                  end if
                  if ( jz .eq. jhalt ) then
                    j2 = jz
                  end if
                end if

              else
                umspei(jj) = c(i,j) + umspei(jj)
              end if

            end if
          end do
        end if
      end do
c
c  Computation of a bound by solving LSAPI with cost matrix UMSPEI.
c
5180  continue

      call lsapi ( nmk, unendl, umspei, zstern, y, z1, lab, dd, u, v,
     &  h1, hl1 )

      jj = 0

      do i = 1, nmk
        do j = 1, nmk
          jj = jj + 1
          umspei(jj) = umspei(jj) - u(i) - y(j)
        end do
      end do
c
c  ZPART is the fixed fraction of the objective function value implied
c  by the present partial permutation.
c  The present bound is ZPART + ZSTERN.
c
      if ( olwert .le. zpart + zstern ) then
        go to 5250
      end if

      if ( weiter ) then
        go to 5220
      end if

5135  continue

      if ( ikap .eq. 1 ) then
        go to 5410
      end if

      if ( 1 .lt. ikap ) then
        go to 5420
      end if
c
c  Computation of the alternative costs.
c
      call altkos ( nmk, umspei, z1, unendl, izaehl, jzaehl, alterk )
      go to 5490
c
c  Computation of the next single assignment in a fixed row or column
c
5410  continue

      min = unendl
      j1 = z1(ihalt)
      jj = ( ihalt - 1 ) * nmk

      do j = 1, nmk
        jj = jj + 1
        t = umspei(jj)
        if ( t .lt. min .and. j .ne. j1 ) then
          min = t
        end if
      end do

      min1 = min
      min = unendl
      jj = j1

      do i = 1, nmk
        t = umspei(jj)
        jj = jj + nmk
        if ( t .lt. min .and. i .ne. ihalt ) then
          min = t
        end if
      end do

      min1 = min1 + min
      min = unendl
      jj = jhalt

      do i = 1, nmk
        t = umspei(jj)
        jj = jj + nmk
        if ( t .lt. min .and. jhalt .ne. z1(i) ) then
          min = t
        end if
      end do

      min2 = min

      do i = 1, nmk
        if ( z1(i) .eq. jhalt ) then
          i1 = i
          go to 5540
        end if
      end do

5540  continue

      jj = ( i1 - 1 ) * nmk
      min = unendl

      do j = 1, nmk
        jj = jj + 1
        t = umspei(jj)
        if ( t .lt. min .and. j .ne. jhalt ) then
          min = t
        end if
      end do

      if ( min1 .le. min + min2 ) then
        izaehl = i1
        jzaehl = jhalt
        alterk = min + min2
      else
        izaehl = ihalt
        jzaehl = j1
        alterk = min1
      end if

      go to 5490
c
c  Computation of the next single assignment in the previously fixed row
c  or previously fixed column.
c
5420  continue

      if ( i2 .ne. 0 ) then

        izaehl = i2
        jzaehl = z1(izaehl)

      else

        jzaehl = j2

        do i = 1, nmk
          if ( z1(i) .eq. jzaehl ) then
            izaehl = i
            go to 5485
          end if
        end do

5485    continue

      end if

      min = unendl
      jj = ( izaehl - 1 ) * nmk

      do i = 1, nmk
        jj = jj + 1
        t = umspei(jj)
        if ( t .lt. min .and. i .ne. jzaehl ) then
          min = t
        end if
      end do

      alterk = min
      min = unendl
      jj = jzaehl

      do j = 1, nmk
        t = umspei(jj)
        jj = jj + nmk
        if ( t .lt. min .and. j .ne. izaehl ) then
          min = t
        end if
      end do

      alterk = alterk + min

5490  continue

      iz = 0
      alter(k+1) = alterk + zstern

      do i = 1, n
        if ( .not. bool(i) ) then
          iz = iz + 1
          if ( izaehl .eq. iz ) then
            izaehl = i
            go to 5165
          end if
        end if
      end do

5165  continue

      iz = 0

      do j = 1, n
        if ( .not. bool1(j) ) then
          iz = iz + 1
          if ( jzaehl .eq. iz ) then
            jzaehl = j
            go to 5185
          end if
        end if
      end do

5185  continue

      zul(izaehl,jzaehl) = k
      weiter = .true.
      bool(izaehl) = .true.
      bool1(izaehl) = .true.
      nmk = n - k - 1
c
c  Computation of the cost matrix corresponding to the new partial
c  permutation.
c
      call wegspe ( n, k, nmk, izaehl, jzaehl, a, b, veksum, u, bool,
     &  bool1, aspei, bspei, h1 )

      call progno ( n, k, nmk, umspei, veksum, vekqu, weiter, aspei,
     &  bspei, cspei )

      iz = 0

      do i = 1, n
        if ( .not. bool(i) ) then
          t1 = a(i,izaehl)
          t2 = a(izaehl,i)
          jz = iz
          do j = 1, n
            if ( .not. bool1(j) ) then
              jz = jz + 1
              umspei(jz) = c(i,j) + t1 * b(j,jzaehl) + t2 * b(jzaehl,j)
     &          + umspei(jz)
            end if
          end do
          iz = iz + nmk
        end if
      end do

      zpart = zpart + c(izaehl,jzaehl)
      go to 5180
c
c  The bound for the new partial permutation is not less than a
c  previously bound objective function value.
c  Backtrack.
c
5250  continue

      if ( weiter ) then
        weiter = .false.
        k = k + 1
        go to 5230
      end if
c
c  Exit.  The solution tree is completely fathomed.
c
5255  continue

      if ( k .eq. 0 ) then
        return
      end if

      izaehl = menge(k)
      jzaehl = phiofm(k)

5220  continue

      do i = 1, n
        if ( .not. bool(i) ) then
          t1 = a(izaehl,i)
          t2 = a(i,izaehl)
          do j = 1, n
            if ( .not. bool1(j) ) then
              t = t1 * b(jzaehl,j) + t2 * b(j,jzaehl)
              if ( .not. weiter ) then
                t = - t
              end if
              c(i,j) = c(i,j) + t
            end if
          end do
        end if
      end do

      if ( weiter ) then
        partpe(izaehl) = jzaehl
        k = k + 1
        menge(k) = izaehl
        phiofm(k) = jzaehl
        if ( k .eq. n - 2 ) then
          go to 5270
        else
          ikap = 0
          go to 5135
        end if
      end if
c
c  Cancellation of the last single assignment.
c
5230  continue

      do i = 1, n
        if ( .not. bool(i) ) then
          do j = 1, n
            if ( .not. bool1(j) .and. zul(i,j) .eq. k ) then
              zul(i,j) = -1
            end if
          end do
        end if
      end do

      zpart = zpart - c(izaehl,jzaehl)
      bool(izaehl) = .false.
      bool1(jzaehl) = .false.
      k = k - 1
      nmk = n - k

      if ( olwert .le. alter(k+1) + zpart ) then
        go to 5255
      end if

      call progno ( n, k, nmk, umspei, veksum, vekqu, weiter, aspei,
     &  bspei, cspei )
      go to 5210
c
c  Computation of the objective function values for the remaining
c  two complete permutations.
c
5270  continue

      do i = 1, n
        if ( .not. bool(i) ) then
          iz = i
          go to 5285
        end if
      end do

5285  continue

      do i = 1, n
        if ( .not. bool1(i) ) then
          j = i
          go to 5295
        end if
      end do

5295  continue

      bool(iz) = .true.
      bool1(j) = .true.

      do i = 1, n
        if ( .not. bool(i) ) then
          i1 = i
          go to 5305
        end if
      end do

5305  continue

      do i = 1, n
        if ( .not. bool1(i) ) then
          j1 = i
          go to 5315
        end if
      end do

5315  continue

      weiter = .false.
      ize = 0

5330  continue

      zstern = c(iz,j) + c(i1,j1) + a(iz,i1) * b(j,j1) 
     &  + a(i1,iz) * b(j1,j)

      bool(iz) = .false.
      bool1(j) = .false.

      if ( zstern + zpart .lt. olwert ) then
        olwert = zstern + zpart
        do i = 1, n
          if ( bool(i) ) then
            loesg(i) = partpe(i)
          end if
        end do
        do i = 1, nm2
          meng(i) = menge(i)
        end do
        loesg(iz) = j
        loesg(i1) = j1
      end if

      if ( ize .ne. 0 ) then
        go to 5255
      end if

      ize = iz
      iz = i1
      i1 = ize
      go to 5330

      end
      subroutine qaph1 ( n, a, b, unendl, rep, perm, olwert, menge,
     &  phiofm, iperm, bool, bool1, seed )

c*********************************************************************72
c
cc QAPH1 is a heuristic solver for quadratic assignment problems.
c
c  Modified:
c
c    27 November 2007
c
c  Author:
c
c    T Boenniger,
c    Karl-Heinz Stratmann.
c
c  Reference:
c
c    Rainer Burkard, Ulrich Derigs,
c    Assignment and Matching Problems: Solution Methods with Fortran-Programs,
c    Springer, 1980,
c    ISBN: 3-540-10267-1,
c    LC: QA402.5 B86.
c
c  Parameters:
c
c    Input, integer N, the dimension of the problem.
c
c    Input, integer A(N,N), the distance matrix.  All entries should 
c    be nonnegative.
c
c    Input, integer B(N,N), the connection matrix.  All entries should
c    be nonnegative.
c
c    Input, integer UNENDL, a very large integer.
c
c    Input, integer REP, the number of restarts allowed.
c
c    Output, integer PERM(N), the best permutation found.
c
c    Output, integer OLWERT, the objective value.
c
c    Workspace, integer MENGE(N), PHIOFM(N), IPERM(N).
c
c    Workspace, logical BOOL(N), BOOL1(N).
c
c    Input/output, integer SEED, a seed for the random number generator.
c
      implicit none

      integer n

      integer a(n,n)
      integer b(n,n)
      logical bool(n)
      logical bool1(n)
      integer i
      integer i0
      integer ifest
      integer im
      integer iperm(n)
      integer ir
      integer j
      integer j0
      integer j1
      integer jfest
      integer jm
      integer k
      integer k1
      integer kfest
      integer km1
      integer lfest
      integer menge(n)
      integer mfest
      integer min
      integer nmk
      integer olwert
      integer perm(n)
      integer phiofm(n)
      integer rep
      integer seed
      integer sum
      integer sum1
      integer unendl
c
c  Determination of a random permutation.
c
      olwert = unendl

      do ir = 1, rep

        call zufall ( n, iperm, bool, seed )

        do i = 1, n
          bool(i) = .false.
          bool1(i) = .false.
        end do

        i0 = iperm(1)
        j0 = 1
        menge(1) = i0
        phiofm(1) = j0
        bool(i0) = .true.
        bool1(j0) = .true.
        k = 1

7250    continue

        min = unendl
        i0 = iperm(k+1)
        nmk = i0
c
c  Evaluation of the sum.
c
        do k1 = 1, k + 1

          if ( 1 .lt. k1 ) then

            ifest = i0
            i0 = menge(k1-1)
            menge(k1-1) = ifest

          end if

          sum = 0

          do i = 1, k - 1

            im = menge(i)
            jm = phiofm(i)
            sum = sum + a(im,im) * b(jm,jm)

            do j = i + 1, k
              lfest = menge(j)
              mfest = phiofm(j)
              sum = sum + a(im,lfest) * b(jm,mfest)
     &                  + a(lfest,im) * b(mfest,jm)
            end do

          end do

          if ( 1 .lt.  k ) then
            im = menge(k)
            jm = phiofm(k)
            sum = sum + a(im,im) * b(jm,jm)
          end if

          do i = 1, n

            if ( .not. bool1(i) ) then

              sum1 = a(i0,i0) * b(i,i)

              do j1 = 1, k
                im = menge(j1)
                jm = phiofm(j1)
                sum1 = sum1 + a(im,i0) * b(jm,i)
     &                      + a(i0,im) * b(i,jm)
              end do

              sum1 = sum1 + sum

              if ( sum1 .le. min ) then
                min = sum1
                kfest = k1
                jfest = i
              end if

            end if
            
          end do

          if ( 1 .lt. k1 ) then
            ifest = i0
            i0 = menge(k1-1)
            menge(k1-1) = ifest
          end if

        end do
c
c  Determination of the partial assignment with least costs.
c
        if ( kfest .ne. 1 ) then

          km1 = kfest - 1
          ifest = i0
          i0 = menge(km1)
          menge(km1) = ifest

        end if

        k = k + 1
        menge(k) = i0
        phiofm(k) = jfest
        bool(nmk) = .true.
        bool1(jfest) = .true.

        if ( k .lt. n ) then
          go to 7250
        end if
c
c  Determination of the final solution.
c
        if ( min .le. olwert ) then
          olwert = min
          do i = 1, n
            perm(menge(i)) = phiofm(i)
          end do
        end if

      end do

      return
      end
      subroutine qaph2 ( n, a, b, miter, rep, unendl, eps, ope, startp,
     &  ol, c, rest0, rest, h1, h2, h3, h4, ap, op, areal, lambda, u,
     &  v, dd, hr, bool, seed )

c*********************************************************************72
c
cc QAPH2 is a heuristic solver for quadratic assignment problems.
c
c  Discussion:
c
c    The method used involves the heuristic cutting plane and
c    exchange method.
c
c  Modified:
c
c    02 December 2007
c
c  Author:
c
c    T Boenniger
c
c  Reference:
c
c    Rainer Burkard, Ulrich Derigs,
c    Assignment and Matching Problems: Solution Methods with Fortran-Programs,
c    Springer, 1980,
c    ISBN: 3-540-10267-1,
c    LC: QA402.5 B86.
c
c  Parameters:
c
c    Input, integer N, the dimension of the problem.
c
c    Input, integer A(N,N), the distance matrix.
c
c    Input, integer B(N,N), the connection matrix.
c
c    Input, integer MITER, the number of iterations to take.
c
c    Input, integer REP, the number of restarts to use.
c
c    Input, integer UNENDL, a large number.
c
c    Input, real EPS, the machine accuracy.
c
c    Input/output, integer OPE(N).  On input with STARTP set TRUE, the
c    value stored in OPE is used as the starting permutation.  If STARTP
c    is FALSE, the input value is ignore.  On output, OPE contains the
c    best permutation that was found.
c
c    Input, logical STARTP, is TRUE if a starting permutation is
c    given in OPE.
c
c    Output, integer OL, the value of the objective function corresponding
c    to the output value of OPE.
c
c    Workspace, integer C(N,N), REST0(MITER), REST(MITER,N*N), H1(N*N),
c    H2(N), H3(N), H4(N),AP(N), OP(N).
c
c    Workspace, real AREAL(N*N), LAMBDA(N*N), U(N), V(N), DD(N), HR(N).
c
c    Workspace, logical BOOL(N).
c
c    Input/output, integer SEED, a seed for the random number generator.
c
      implicit none

      integer miter
      integer n

      integer a(n,n)
      integer alpha
      integer ap(n)
      real areal(n*n)
      integer b(n,n)
      real beta
      logical bool(n)
      integer c(n*n)
      real dd(n)
      real eps
      integer h1(n*n)
      integer h2(n)
      integer h3(n)
      integer h4(n)
      real hr(n)
      integer i
      integer ii
      integer ij
      integer irep
      integer izaehl
      integer j
      integer k
      real lambda(n*n)
      integer lbb
      integer ol
      integer op(n)
      integer ope(n)
      integer r
      integer rep
      integer rest(miter,n*n)
      integer rest0(miter)
      integer seed
      integer sp
      logical startp
      real u(n)
      integer ub
      integer ub1
      integer unendl
      real v(n)
      integer zstern

      r = 1
      irep = 1
      ub = unendl
c
c  Determination of the random starting permutation.
c
      if ( .not. startp ) then
        call zufall ( n, ope, bool, seed )
      end if

      do i = 1, n
        op(i) = ope(i)
      end do
c
c  Computation of the constants.
c  C(I,J) = N * max ( A(I,K) ) * max ( B(J,L) )
c
      do i = 1, n
        h1(i) = 0
        h2(i) = 0
        do k = 1, n
          h1(i) = max ( h1(i), a(i,k) )
          h2(i) = max ( h2(i), b(i,k) )
        end do
      end do

      izaehl = 0
      do i = 1, n
        do j = 1, n
          izaehl = izaehl + 1
          lambda(izaehl) = 0.0E+00
          c(izaehl) = n * h1(i) * h2(j)
        end do
      end do

      go to 1531

120   continue
c
c  Solution of the master problem for R = 1.
c
      ub1 = 0
      do i = 1, n
        do k = 1, n
          ub1 = ub1 + a(i,k) * b(op(i),op(k))
        end do
      end do

      if ( ub1 .lt. ub ) then
        ub = ub1
        do i = 1, n
          ope(i) = op(i)
        end do
      end if

      go to 150
c
c  Isolated chunk of code:
c
130   continue

      do i = 1, n
        op(i) = ap(i)
      end do

150   continue

      r = r + 1
      if ( mod ( r, miter ) .ne. 0 ) then
        go to 156
      end if

153   continue

      if ( irep .eq. rep ) then
        ol = ub
        return
      end if

      irep = irep + 1
      r = 1
      call zufall ( n, op, bool, seed )
c
c  Computation of the corresponding objective function.
c
1531  continue

      ub1 = 0
      izaehl = 0
      do i = 1, n
        do k = 1, n
          izaehl = izaehl + 1
          lambda(izaehl) = 0.0E+00
          ub1 = ub1 + a(i,k) * b(op(i),op(k))
        end do
      end do
c
c  Pair and triple exchange.
c
      call reihen ( n, a, b, op, h1, h2, h3, h4, bool )
      call heider ( n, a, b, ub1, op, h1 )
      call dreier ( n, a, b, ub1, op, h1, bool )

      if ( ub1 .lt. ub ) then
        ub = ub1
        do i = 1, n
          ope(i) = op(i)
        end do
      end if

156   continue
c
c  Determination of the new cut.
c
      alpha = 0
      izaehl = 0

      do i = 1, n

        do j = 1, n

          izaehl = izaehl + 1
          sp = 0

          do k = 1, n
            sp = sp + a(k,i) * b(op(k),j)
          end do

          if ( j .eq. op(i) ) then
            sp = sp + c(izaehl)
            alpha = alpha + c(izaehl)
          end if

          areal(izaehl) = float ( sp )
          rest(r,izaehl) = sp

        end do

      end do

      rest0(r) = alpha
c
c  Heurtistic solution of the restricted master problem.
c  Determination of BETA and of the vector LAMBDA(I).
c
      call lsapr ( n, unendl, areal, beta, h1, op, u, v, dd, hr,
     & h2, bool, eps )

      lbb = ifix ( beta ) - alpha
      beta = float ( lbb )
      beta = abs ( beta )
      if ( beta .eq. 0.0E+00 ) then
        beta = 1.0E+00
      end if

      do i = 1, n * n
        lambda(i) = lambda(i) + areal(i) / beta
      end do
c
c  For R = 1, LBB is the exact solution of the restricted master problem.
c
      if ( r .eq. 1 ) then
        go to 120
      end if

      call lsapr ( n, unendl, lambda, beta, h1, ap, u, v, dd, hr,
     &  h2, bool, eps )
c
c  Computation of the objective function value corresponding to AP(I).
c
      ub1 = 0
      do i = 1, n
        do k = 1, n
          ub1 = ub1 + a(i,k) * b(ap(i),ap(k))
        end do
      end do

      if ( ub1 .lt. ub ) then
         ub = ub1
         do i = 1, n
           ope(i) = op(i)
         end do
      end if
c
c  Improvement by pairwise exchange.
c
      do i = 1, n
        op(i) = ap(i)
      end do

      call reihen ( n, a, b, op, h1, h2, h3, h4, bool )
      call heider ( n, a, b, ub1, op, h1 )

      if ( ub1 .lt. ub ) then
        ub = ub1
        do i = 1, n
          ope(i) = op(i)
        end do
      end if
c
c  Determination of Z = max ( A(P)*X-ALPHA(P)) for X = OP.
c
      zstern = 0
      do k = 1, r
        sp = - rest0(k)
        ii = 0
        do i = 1, n
          ij = ii + op(i)
          sp = sp + rest(k,ij)
          ii = ii + n
        end do
        zstern = max ( zstern, sp )
      end do

      if ( zstern .lt. ub ) then
        go to 150
      end if
c
c  Determination of Z = max ( A(P)*X-ALPHA(P)) for X = AP.
c
      zstern = 0
      do k = 1, r
        sp = - rest0(k)
        ii = 0
        do i = 1, n
          ij = ii + ap(i)
          sp = sp + rest(k,ij)
          ii = ii + n
        end do
        zstern = max ( zstern, sp )
      end do

      if ( zstern .lt. ub ) then
        go to 130
      end if

      if ( irep .lt. rep ) then
        go to 153
      end if

      ol = ub

      return
      end
      function r4_uniform_01 ( seed )

c*********************************************************************72
c
cc R4_UNIFORM_01 returns a unit pseudorandom R4.
c
c  Discussion:
c
c    This routine implements the recursion
c
c      seed = 16807 * seed mod ( 2**31 - 1 )
c      r4_uniform_01 = seed / ( 2**31 - 1 )
c
c    The integer arithmetic never requires more than 32 bits,
c    including a sign bit.
c
c    If the initial seed is 12345, then the first three computations are
c
c      Input     Output      R4_UNIFORM_01
c      SEED      SEED
c
c         12345   207482415  0.096616
c     207482415  1790989824  0.833995
c    1790989824  2035175616  0.947702
c
c  Modified:
c
c    11 August 2004
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Paul Bratley, Bennett Fox, Linus Schrage,
c    A Guide to Simulation,
c    Second Edition,
c    Springer, 1987,
c    ISBN: 0387964673,
c    LC: QA76.9.C65.B73.
c
c    Bennett Fox,
c    Algorithm 647:
c    Implementation and Relative Efficiency of Quasirandom
c    Sequence Generators,
c    ACM Transactions on Mathematical Software,
c    Volume 12, Number 4, December 1986, pages 362-376.
c
c    Pierre L'Ecuyer,
c    Random Number Generation,
c    in Handbook of Simulation,
c    edited by Jerry Banks,
c    Wiley, 1998,
c    ISBN: 0471134031,
c    LC: T57.62.H37.
c
c    Peter Lewis, Allen Goodman, James Miller,
c    A Pseudo-Random Number Generator for the System/360,
c    IBM Systems Journal,
c    Volume 8, Number 2, 1969, pages 136-143.
c
c  Parameters:
c
c    Input/output, integer SEED, the "seed" value, which should NOT be 0.
c    On output, SEED has been updated.
c
c    Output, real R4_UNIFORM_01, a new pseudorandom variate,
c    strictly between 0 and 1.
c
      implicit none

      integer k
      integer seed
      real r4_uniform_01

      if ( seed .eq. 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'R4_UNIFORM_01 - Fatal error!'
        write ( *, '(a)' ) '  Input value of SEED = 0.'
        stop
      end if

      k = seed / 127773

      seed = 16807 * ( seed - k * 127773 ) - k * 2836

      if ( seed .lt. 0 ) then
        seed = seed + 2147483647
      end if
c
c  Although SEED can be represented exactly as a 32 bit integer,
c  it generally cannot be represented exactly as a 32 bit real number!
c
      r4_uniform_01 = real ( dble ( seed ) * 4.656612875D-10 )

      return
      end
      subroutine reihen ( n, a, b, perm, ai, bj, ai1, bj1, bool )

c*********************************************************************72
c
cc REIHEN choose a start sequence for QAPH2'S pairwise exchange algorithm.
c
c  Modified:
c
c    29 November 2007
c
c  Author:
c
c    Ulrich Derigs
c
c  Reference:
c
c    Rainer Burkard, Ulrich Derigs,
c    Assignment and Matching Problems: Solution Methods with Fortran-Programs,
c    Springer, 1980,
c    ISBN: 3-540-10267-1,
c    LC: QA402.5 B86.
c
c  Parameters:
c
c    Input, integer N, dimension of the problem.
c
c    Input, integer A(N,N), the distance matrix.
c
c    Input, integer B(N,N), the connection matrix.
c
c    Input, integer PERM(N), the permutation.
c
c    Output, integer AI(N), the start sequence for the pairwise
c    exchange algorithm.
c
c    Workspace, integer BJ(N).
c
c    Workspace, integer AI1(N).
c
c    Workspace, integer BJ1(N).
c
c    Unused, logical BOOL(N).
c
      implicit none

      integer n

      integer a(n,n)
      integer ai(n)
      integer ai1(n)
      integer b(n,n)
      integer bj(n)
      integer bj1(n)
      logical bool(n)
      integer i
      integer i1
      integer j
      integer j0
      integer j1
      integer k
      integer perm(n)
      integer sum

      do i = 1, n

        j1 = 0
        do j = 1, n
          if ( i .ne. j ) then
            j1 = j1 + 1
            bj(j1) = a(i,j) + a(j,i)
          end if
        end do

        do j = 1, i - 1
          ai1(j) = j
        end do

        do j = i, n - 1
          ai1(j) = j + 1
        end do

        call ssort ( bj, ai1, n - 1 )

        i1 = perm(i)
        j1 = 0

        do j = 1, n
          if ( i .ne. j ) then
            j1 = j1 + 1
            bj(j1) = b(i1,j) + b(j,i1)
          end if
        end do

        do j = 1, i1 - 1
          ai(j) = j
        end do

        do j = i1, n - 1
          ai(j) = j + 1
        end do

        call ssort ( bj, ai, n - 1 )

        sum = 0

        do i1 = 1, n - 1

          k = ai1(i1)
          j0 = perm(k)
          j1 = n

          do j = 1, n - 1
            j1 = j1 - 1
            if ( ai(j1) .eq. j0 ) then
              go to 6740
            end if
          end do

6740      continue

          sum = sum + ( n - 1 - i1 ) * ( n - 1 - j )

        end do

        bj1(i) = sum

      end do

      do j = 1, n
        ai(j) = j
      end do

      call ssort ( bj1, ai, n )

      return
      end
      subroutine scan1 ( nb1, n, sup, cc, p, basis, mem, ka, kb, sm, 
     & tma, tmb, y1, y2, dplus, dminus, m1 )

c*********************************************************************72
c
cc SCAN1 scans a node for the SMP algorithm.
c
c  Modified:
c
c    02 December 2007
c
c  Author:
c
c    Ulrich Derigs
c
c  Reference:
c
c    Rainer Burkard, Ulrich Derigs,
c    Assignment and Matching Problems: Solution Methods with Fortran-Programs,
c    Springer, 1980,
c    ISBN: 3-540-10267-1,
c    LC: QA402.5 B86.
c
c  Parameters:
c
c    Input, integer NB1, ?
c
c    Input, integer N, ?
c
c    Input, real SUP, a large real value.
c
c    Input, integer CC((N*(N-1))/2), the strict upper triangular part
c    of the cost matrix, stored by columns.
c
c    Input, integer P(*), ?
c
c    input, integer BASIS(N), ?
c
c    Input, integer MEM(N), ?
c
c    Output, integer KA(N), ?
c
c    Output, integer KB(N), ?
c
c    Not used, integer SM(*), ?
c
c    Input/output, integer TMA(*), ?
c
c    Input, integer TMB(*), ?
c
c    Input, real Y1(N), ?
c
c    Input, real Y2(N), ?
c
c    Input, real DPLUS(N), ?
c
c    Input/output, real DMINUS(N), ?
c
c    Input/output, integer M1(N), ?
c
      implicit none

      integer n

      integer basis(n)
      real c0
      integer cc(*)
      real d1
      real d2
      real dminus(n)
      real dplus(n)
      integer i1
      integer i2
      integer ind
      integer ka(n)
      integer kb(n)
      integer m1(n)
      integer max
      integer mem(n)
      integer min
      integer n1
      integer n2
      integer nb1
      integer nb2
      integer nc
      integer p(*)
      integer sm(*)
      real sup
      integer tma(*)
      integer tmb(*)
      integer top
      real y1(n)
      real y2(n)

      top = n + 2
      d1 = dplus(nb1) - y1(nb1)
      dminus(nb1) = sup
      d2 = d1 - y2(nb1)
      tma(nb1) = 0
      i1 = 0

      do n2 = 1, n

        nb2 = basis(n2)

        if ( top .le. tma(nb2) ) then

          i1 = i1 + 1
          m1(i1) = n2
          min = nb1
          max = nb2

          if ( nb1 .ne. n2 ) then

            if ( n2 .lt. nb1 ) then
              max = nb1
              min = n2
            end if

            ind = p(max) + min
            nc = cc(ind)
            c0 = float ( nc ) + d2
            c0 = c0 - y1(nb2) - y2(n2)

            if ( c0 .lt. dminus(nb2) ) then
              ka(nb2) = nb1
              kb(nb2) = n2
              dminus(nb2) = c0
            end if

          end if
        end if

      end do

      tma(nb1) = top
      n1 = nb1
      go to 315

305   continue

      d2 = d1 - y2(n1)

      do i2 = 1, i1

        n2 = m1(i2)
        nb2 = basis(n2)
        min = n1
        max = n2

        if ( n1 .ne. n2 ) then

          if ( n2 .lt. n1 ) then
            max = n1
            min = n2
          end if

          ind = p(max) + min
          nc = cc(ind)
          c0 = float ( nc ) + d2
          c0 = c0 - y1(nb2) - y2(n2)

          if ( c0 .lt. dminus(nb2) ) then
            ka(nb2) = n1
            kb(nb2) = n2
            dminus(nb2) = c0
          end if

        end if

      end do

315   continue

      if ( n1 .ne. nb1 ) then
        go to 305
      end if

      return
      end
      subroutine scan2 ( nb, n, sup, cc, p, basis, mem, ka, kb, sm, 
     & tma, tmb, y1, y2, dplus, dminus )

c*********************************************************************72
c
cc SCAN2 scans a node for the SMP algorithm.
c
c  Modified:
c
c    02 December 2007
c
c  Author:
c
c    Ulrich Derigs
c
c  Reference:
c
c    Rainer Burkard, Ulrich Derigs,
c    Assignment and Matching Problems: Solution Methods with Fortran-Programs,
c    Springer, 1980,
c    ISBN: 3-540-10267-1,
c    LC: QA402.5 B86.
c
c  Parameters:
c
c    Input, integer NB, ?
c
c    Input, integer N, ?
c
c    Input, real SUP, a large real value.
c
c    Input, integer CC((N*(N-1))/2), the strict upper triangular part
c    of the cost matrix, stored by columns.
c
c    Input, integer P(*), ?
c
c    input, integer BASIS(N), ?
c
c    Input, integer MEM(N), ?
c
c    Output, integer KA(N), ?
c
c    Output, integer KB(N), ?
c
c    Input, integer SM(*), ?
c
c    Input/output, integer TMA(*), ?
c
c    Input, integer TMB(*), ?
c
c    Input, real Y1(N), ?
c
c    Input, real Y2(N), ?
c
c    Input, real DPLUS(N), ?
c
c    Input/output, real DMINUS(N), ?
c
      implicit none

      integer n

      integer basis(n)
      real c0
      integer cc(*)
      real d
      real dminus(n)
      real dplus(n)
      integer ind
      integer ka(n)
      integer kb(n)
      integer max
      integer mem(n)
      integer min
      integer n1
      integer n2
      integer nb
      integer nb1
      integer nb2
      integer nc
      integer nka
      integer nkb
      integer p(*)
      integer sm(*)
      real sup
      integer tma(*)
      integer tmb(*)
      integer top
      real y1(n)
      real y1b
      real y2(n)
      real y2b

      top = n + 2

300   continue

      nb1 = nb
      nb = tmb(nb1)
      tmb(nb1) = top
      d = sup
      nka = 0
      nkb = 0
      n1 = nb1
      y1b = y1(nb1)

315   continue

      y2b = y2(n1)

      do n2 = 1, n

        nb2 = basis(n2)

        if ( sm(nb2) .lt. top ) then

          min = n1
          max = n2

          if ( n1 .ne. n2 ) then

            if ( n2 .lt. n1 ) then
              max = n1
              min = n2
            end if

            ind = p(max) + min
            nc = cc(ind)
            c0 = float ( nc ) - y1b - y2b
            c0 = c0 - y1(nb2) - y2(n2)
            c0 = c0 + dplus(nb2)

            if ( c0 .lt. d ) then
              nka = n2
              nkb = n1
              d = c0
            end if

          end if

        end if

      end do

      n1 = mem(n1)

      if ( n1 .ne. nb1 ) then
        go to 315
      end if

      ka(nb1) = nka
      kb(nb1) = nkb
      dminus(nb1) = d

      if ( nb .ne. 0 ) then
        go to 300
      end if

      return
      end
      subroutine schr ( n, basis, mem, y1, y2, bb, kkn )

c*********************************************************************72
c
cc SCHR shrinks a blossom for CPP.
c
c  Modified:
c
c    06 December 2007
c
c  Author:
c
c    W Puetz
c
c  Reference:
c
c    Rainer Burkard, Ulrich Derigs,
c    Assignment and Matching Problems: Solution Methods with Fortran-Programs,
c    Springer, 1980,
c    ISBN: 3-540-10267-1,
c    LC: QA402.5 B86.
c
c  Parameters:
c
c    Unused, integer N, ?
c
c    Output, integer BASIS(*), ?
c
c    Input, integer MEM(*), ?
c
c    Input, real Y1(*), ?
c
c    Input/output, real Y2(*), ?
c
c    Input, integer BB, ?
c
c    Input, integer KKN, ?
c
      implicit none

      integer basis(*)
      integer bb
      integer kkn
      integer kkn1
      integer kkn2
      integer mem(*)
      integer n
      real y1(*)
      real y2(*)
      real yy

      kkn1 = kkn
      yy = y1(kkn)

9700  continue

      basis(kkn) = bb
      y2(kkn) = y2(kkn) + yy
      kkn2 = mem(kkn)

      if ( kkn2 .ne. kkn1 ) then
        kkn = kkn2
        go to 9700
      end if

      return
      end
      subroutine smp ( n, zfw, nmatch, cc, p, basis, mem, ka, kb, sm, 
     &  tma, tmb, m1, y1, y2, dplus, dminus, sup, eps )

c*********************************************************************72
c
cc SMP solves the sum matching problem.
c
c  Modified:
c
c    03 December 2007
c
c  Author:
c
c    G Kazakidis
c
c  Reference:
c
c    Rainer Burkard, Ulrich Derigs,
c    Assignment and Matching Problems: Solution Methods with Fortran-Programs,
c    Springer, 1980,
c    ISBN: 3-540-10267-1,
c    LC: QA402.5 B86.
c
c  Parameters:
c
c    Input, integer N, the dimension of the vector.
c
c    Output, integer ZFW, the cost of the optimal matching.
c
c    Output, integer NMATCH(N), the optimal matching.
c
c    Input, integer CC((N*(N-1))/2), the strict upper triangular part
c    of the cost matrix, stored by columns.
c
c    Workspace, integer P(N), BASIS(N), MEM(N), KA(N), KB(N), SM(N),
c    TMA(N), TMB(N), M1(N).
c
c    Workspace, real Y1(N), Y2(N), DPLUS(N), DMINUS(N).
c
c    Input, real SUP, a large value.
c
c    Input, real EPS, the machine accuracy.
c
      implicit none

      integer n

      integer basis(n)
      real c0
      integer cc((n*(n-1))/2)
      real d
      real db
      real dbest
      real dminus(n)
      real dplus(n)
      real eps
      integer i
      integer ind
      integer ka(n)
      integer kb(n)
      integer m1(n)
      integer max
      integer mem(n)
      integer min
      integer n1
      integer n2
      integer n3
      integer nb
      integer nb1
      integer nb2
      integer nb3
      integer nbest
      integer nc
      integer ni
      integer nk
      integer nk1
      integer nk2
      integer nk3
      integer nk4
      integer nka
      integer nkb
      integer nm
      integer nmatch(n)
      integer nmb
      integer nn
      integer p(n)
      integer sm(n)
      real sup
      integer tma(n)
      integer tmb(n)
      integer top
      real y1(n)
      real y1b
      real y2(n)
      real y2b
      real yb
      integer zfw
c
c  Start.
c
      p(1) = 0
      p(2) = 0
      do i = 3, n
        p(i) = p(i-1) + i - 2
      end do

      top = n + 2

      do i = 1, n
        basis(i) = i
        mem(i) = i
        y1(i) = 0.0E+00
        y2(i) = 0.0E+00
        sm(i) = top
        tma(i) = top
        tmb(i) = top
        nmatch(i) = top
        dplus(i) = sup
        dminus(i) = sup
        ka(i) = 0
        kb(i) = i
      end do
c
c  Start procedure.
c
      do n1 = 1, n

        if ( nmatch(n1) .eq. top ) then

          nn = 0
          d = sup

          do n2 = 1, n

            min = n1
            max = n2

            if ( n1 .ne. n2 ) then

              if ( n2 .lt. n1 ) then
                max = n1
                min = n2
              end if

              ind = p(max) + min
              nc = cc(ind)
              c0 = float ( nc ) - y1(n2)

              if ( c0 .lt. d ) then
                d = c0
                nn = 0
                if ( nmatch(n2) .eq. top ) then
                  nn = n2
                end if
              else if ( c0 .eq. d ) then
                if ( nn .eq. 0 ) then
                  if ( nmatch(n2) .eq. top ) then
                    nn = n2
                  end if
                end if
              end if

            end if

          end do

          if ( nn .ne. 0 ) then
            y1(n1) = d
            nmatch(n1) = nn
            nmatch(nn) = n1
          end if

        end if

      end do
c
c  Initial labeling.
c
      nn = 0

      do ni = 1, n

        if ( nmatch(ni) .eq. top ) then

          nn = nn + 1
          sm(ni) = 0
          dplus(ni) = 0.0E+00
          y1b = y1(ni)

          do nk = 1, n

            min = ni
            max = nk

            if ( ni .ne. nk ) then
              if ( nk .lt. ni ) then
                max = ni
                min = nk
              end if

              ind = p(max) + min
              nc = cc(ind)
              c0 = float ( nc ) - y1b - y1(nk)

              if ( c0 .lt. dminus(nk) ) then
                dminus(nk) = c0
                ka(nk) = ni
              end if

            end if

          end do

        end if

      end do

      if ( nn .le. 1 ) then
        go to 700
      end if
c
c  Examination of the labeling, and decision for the next step.
c
200   continue

      dbest = sup

      do nb = 1, n

        if ( basis(nb) .eq. nb ) then

          d = dminus(nb)

          if ( sm(nb) .lt. top ) then

            d = 0.5E+00 * ( d + dplus(nb) )

            if ( d .le. dbest ) then
              nbest = nb
              dbest = d
            end if

          else if ( tma(nb) .lt. top ) then

            if ( mem(nb) .ne. nb ) then
              d = d + y1(nb)
              if ( d .lt. dbest ) then
                nbest = nb
                dbest = d
              end if
            end if

          else

            if ( db .lt. dbest ) then
              nbest = nb
              dbest = d
            end if

          end if

        end if

      end do

      if ( tma(nbest) .lt. top ) then
        go to 500
      end if

      if ( top .le. sm(nbest) ) then
        go to 300
      end if

      nka = ka(nbest)
      nkb = kb(nbest)
      n1 = nbest
      nb1 = n1
      n2 = basis(nka)
      nb2 = n2

220   continue

      tma(nb1) = nb2
      nk = sm(nb1)

      if ( nk .ne. 0 ) then
        nb2 = basis(nk)
        nb1 = tma(nb2)
        nb1 = basis(nb1)
        go to 220
      end if

      nb = nb1
      nb1 = n2
      nb2 = n1

230   continue

      if ( top .le. tma(nb1) ) then

        tma(nb1) = nb2
        nk = sm(nb1)

        if ( nk .eq. 0 ) then
          go to 600
        end if

        nb2 = basis(nk)
        nb1 = tma(nb2)
        nb1 = basis(nb1)
        go to 230

      end if

235   continue

      if ( nb1 .eq. nb ) then
        go to 400
      end if

      nk = tma(nb)
      tma(nb) = top
      nm = nmatch(nk)
      nb = basis(nm)
      go to 235
c
c  Growing an alternating tree by adding two edges.
c
300   continue

      tma(nbest) = ka(nbest)
      tmb(nbest) = kb(nbest)
      nm = nmatch(nbest)
      nmb = basis(nm)
      dplus(nmb) = dbest
      sm(nmb) = nmatch(nmb)

      call scan1 ( nmb, n, sup, cc, p, basis, mem, ka, kb, sm, tma,
     &  tmb, y1, y2, dplus, dminus, m1 )

      go to 200
c
c  Shrinking a blossom.
c
400   continue

      yb = y1(nb) + dbest - dplus(nb)
      y1(nb) = 0.0E+00
      nk1 = nb

430   continue

      y2(nk1) = y2(nk1) + yb
      nk1 = mem(nk1)

      if ( nk1 .ne. nb ) then
        go to 430
      end if

      nk = mem(nb)

      if ( nb .ne. n2 ) then
        go to 436
      end if

435   continue

      n2 = n1
      nb2 = tma(nb)

436   continue

      mem(nk1) = nb2
      nm = nmatch(nb2)
      sm(nb2) = nm
      y1b = y1(nb2) + dminus(nb2) - dbest
      nk1 = nb2

440   continue

      nk2 = nk1
      y2(nk2) = y2(nk2) + y1b
      basis(nk2) = nb
      nk1 = mem(nk2)

      if ( nk1 .ne. nb2 ) then
        go to 440
      end if

      kb(nb2) = nk2
      y1(nb2) = y1b
      nb1 = basis(nm)
      mem(nk2) = nb1
      y1b = y1(nb1) + dbest - dplus(nb1)
      nk2 = nb1

445   continue

      nk1 = nk2
      y2(nk1) = y2(nk1) + y1b
      basis(nk1) = nb
      nk2 = mem(nk1)

      if ( nk2 .ne. nb1 ) then
        go to 445
      end if

      kb(nb1) = nk1
      y1(nb1) = y1b

      if ( n2 .eq. nb1 ) then
        go to 450
      end if

      nb2 = tma(nb1)
      tma(nb1) = tmb(nb2)
      tmb(nb1) = tma(nb2)
      go to 436

450   continue

      if ( n2 .eq. nbest ) then
        go to 455
      end if

      tma(n2) = nkb
      tmb(n2) = nka

      if ( nb .ne. nbest ) then
        go to 435
      end if

      go to 460

455   continue

      tma(nbest) = nka
      tmb(nbest) = nkb

460   continue

      mem(nk1) = nk
      n1 = mem(nb)
      ka(n1) = nk
      dplus(n1) = yb
      tma(nb) = top
      dplus(nb) = dbest

      call scan1 ( nb, n, sup, cc, p, basis, mem, ka, kb, sm, tma,
     &  tmb, y1, y2, dplus, dminus, m1 )

      go to 200
c
c  Expansion of a T-labeled blossom.
c
500   continue

      n1 = mem(nbest)
      nb3 = n1
      nka = ka(n1)
      nk2 = n1

505   continue

      nk1 = nk2
      nkb = kb(nk1)
      y1b = y1(nk1)

510   continue

      basis(nk2) = nk1
      y2(nk2) = y2(nk2) - y1b

      if ( nk2 .eq. nkb ) then
        go to 515
      end if

      nk2 = mem(nk2)

      go to 510

515   continue

      nk2 = mem(nkb)
      mem(nkb) = nk1

      if ( nk2 .ne. nka ) then
        go to 505
      end if

      y1b = dplus(n1)
      y1(nbest) = y1b
      mem(nbest) = nka
      nk2 = nka

520   continue

      y2(nk2) = y2(nk2) - y1b

      if ( nk2 .ne. nbest ) then
        nk2 = mem(nk2)
        go to 520
      end if

!525  continue

      nk1 = nmatch(nbest)
      nb1 = basis(nk1)
      nk2 = sm(nb1)
      nb = basis(nk2)

      if ( nb .eq. nbest ) then
        go to 545
      end if

      nb2 = nb
530   continue

      nk = tma(nb2)
      nb1 = basis(nk)

      if ( nb1 .ne. nbest ) then
        nb2 = sm(nb1)
        nb2 = basis(nb2)
        go to 530
      end if

      tma(nb) = tma(nbest)
      tma(nbest) = tmb(nb2)
      tmb(nb) = tmb(nbest)
      tmb(nbest) = nk

      nk3 = sm(nb)
      nb3 = basis(nk3)
      nk4 = sm(nb3)
      sm(nb) = top
      nmatch(nb) = nk1
      nb1 = nb3

540   continue

      nk1 = tma(nb1)
      nk2 = tmb(nb1)
      tma(nb1) = nk4
      tmb(nb1) = nk3
      sm(nb1) = nk1
      nmatch(nb1) = nk1
      nb2 = basis(nk1)
      nmatch(nb2) = nk2
      nk3 = sm(nb2)
      sm(nb2) = nk2

      if ( nb2 .ne. nbest ) then
        nb1 = basis(nk3)
        nk4 = sm(nb1)
        tma(nb2) = nk3
        tmb(nb2) = nk4
        go to 540
      end if

545   continue

      nk2 = tmb(nb)
      nb1 = basis(nk2)
      dminus(nb1) = dbest
      n1 = 0

      if ( nb1 .eq. nb ) then
        go to 555
      end if

      nk1 = tma(nb1)
      nb3 = basis(nk1)
      tma(nb1) = tma(nb)
      tma(nb1) = nk2

550   continue

      nk = sm(nb1)
      sm(nb1) = top
      nb2 = basis(nk)
      nk = tma(nb2)
      tma(nb2) = top
      n2 = tmb(nb2)
      tmb(nb2) = n1
      n1 = nb2
      dplus(nb2) = dbest
      nb1 = basis(nk)
      dminus(nb1) = dbest

      if ( nb1 .ne. nb ) then
        go to 550
      end if

      tma(nb) = n2
      tmb(nb) = nk
      sm(nb) = top

      if ( nb3 .eq. nb ) then
        go to 570
      end if

555   continue

      nb1 = 0
      nb2 = nb3

560   continue

      nk = sm(nb2)
      sm(nb2) = top
      tma(nb2) = top
      tmb(nb2) = nb1
      nb1 = basis(nk)
      nk = tma(nb1)
      sm(nb1) = top
      tma(nb1) = top
      tmb(nb1) = nb2
      nb2 = basis(nk)

      if ( nb2 .ne. nb ) then
        go to 560
      end if

      call scan2 ( nb1, n, sup, cc, p, basis, mem, ka, kb, sm, tma, 
     &  tmb, y1, y2, dplus, dminus )

570   continue

575   continue

      if ( n1 .eq. 0 ) then
        go to 200
      end if

      nb = n1

      call scan1 ( nb, n, sup, cc, p, basis, mem, ka, kb, sm, tma,
     &  tmb, y1, y2, dplus, dminus, m1 )

      n1 = tmb(nb)
      tmb(nb) = top
      go to 575
c
c  Augmentation of the matching.
c  Exchange of the matching and non-matching edges along the
c  augmenting path.
c
600   continue

      nb = n1
      nk = nka

605   continue

      nb1 = nb

606   continue

      nmatch(nb1) = nk
      nk = sm(nb1)
      tma(nb1) = top

      if ( nk .ne. 0 ) then
        nb2 = basis(nk)
        nk1 = tma(nb2)
        nk = tmb(nb2)
        nb1 = basis(nk1)
        nmatch(nb2) = nk1
        go to 606
      end if

      if ( nb .eq. n1 ) then
        nb = n2
        nk = nkb
        go to 605
      end if
c
c  Removing all labels on  non-exposed base nodes.
c
      do nb = 1, n

        if ( basis(nb) .eq. nb ) then

          if ( sm(nb) .lt. top ) then

            d = dbest - dplus(nb)
            y1(nb) = y1(nb) + d
            sm(nb) = top

            if ( nmatch(nb) .ne. top ) then
              dplus(nb) = sup
              dminus(nb) = sup
            else
              sm(nb) = 0
              dplus(nb) = 0.0E+00
              dminus(nb) = sup
            end if

          else 

            if ( tma(nb) .lt. top ) then

              d = dminus(nb) - dbest
              y1(nb ) = y1(nb) + d
              tma(nb) = top
              tmb(nb) = top

            end if    

            dplus(nb) = sup
            dminus(nb) = sup

          end if

        end if

      end do

      nn = nn - 2

      if ( nn .le. 1 ) then
        go to 700
      end if
c
c  Determination of the new DMINUS values.
c
      do n1 = 1, n

        nb1 = basis(n1)

        if ( sm(nb1) .eq. 0 ) then

          y1b = y1(nb1)
          y2b = y2(n1)

          do n2 = 1, n

            nb2 = basis(n2)

            if ( nb1 .ne. nb2 ) then

              min = n1
              max = n2

              if ( n1 .ne. n2 ) then

                if ( n2 .lt. n1 ) then
                  max = n1
                  min = n2
                end if

                ind = p(max) + min
                nc = cc(ind)
                c0 = float ( nc ) - y1b - y2b
                c0 = c0 - y1(nb2) - y2(n2)

                if ( c0 .lt. dminus(nb2) ) then
                  ka(nb2) = n1
                  kb(nb2) = n2
                  dminus(nb2) = c0
                end if

              end if

            end if

          end do

        end if

      end do

      go to 200
c
c  Generation of the original graph by expansion of all shrunken blossoms.
c
700   continue

      zfw = 0

      do nb1 = 1, n

        if ( basis(nb1) .eq. nb1 ) then

          if ( 0 .le. sm(nb1) ) then

            n2 = nmatch(nb1)
            nb2 = basis(n2)
            n1 = nmatch(nb2)
            sm(nb1) = -1
            sm(nb2) = -1
            min = n1
            max = n2

            if ( n1 .ne. n2 ) then

              if ( n2 .lt. n1 ) then
                max = n1
                min = n2
              end if

              ind = p(max) + min
              nc = cc(ind)
              d = float ( nc ) - y1(nb1) - y1(nb2)
              d = d - y2(n1) - y2(n2)

              if ( eps .lt. abs ( d ) ) then
                write ( *, '(a)' ) ' '
                write ( *, '(a)' ) 'SMP - Fatal error!'
                write ( *, '(a)' ) '  Optimality conditions violated.'
                write ( *, '(a,i4,a,i4)' ) 
     &            '  Edge N1 = ', n1, ' N2 = ', n2
                stop
              end if

              zfw = zfw + nc

            end if

          end if
        end if
      end do

      do n1 = 1, n

705     continue

        nb = basis(n1)

        if ( nb .ne. n1 ) then

          nk2 = mem(nb)
          nka = ka(nk2)
          nb3 = nk2
          yb = dplus(nk2)

710       continue

          nk1 = nk2
          nkb = kb(nk1)
          y1b = y1(nk1)

715       continue

          basis(nk2) = nk1
          y2(nk2) = y2(nk2) - y1b

          if ( nk2 .ne. nkb ) then
            nk2 = mem(nk2)
            go to 715
          end if

          nk2 = mem(nkb)
          mem(nkb) = nk1

          if ( nk2 .ne. nka ) then
            go to 710
          end if

          y1(nb) = yb
          mem(nb) = nka
          nk2 = nka

725       continue

          y2(nk2) = y2(nk2) - yb

          if ( nk2 .ne. nb ) then
            nk2 = mem(nk2)
            go to 725
          end if

          nk = nmatch(nb)
          nk1 = basis(nk)
          nk1 = nmatch(nk1)
          nb1 = basis(nk1)

          if ( nb .eq. nb1 ) then
            go to 745
          end if

          nmatch(nb1) = nk
          nb3 = tma(nb1)
          nb3 = basis(nb3)

735       continue

          nk3 = sm(nb1)
          nb2 = basis(nk3)
          nk1 = tma(nb2)
          nk2 = tmb(nb2)
          nb1 = basis(nk1)
          nmatch(nb1) = nk2
          nmatch(nb2) = nk1
          min = nk1
          max = nk2

          if ( nk1 .eq. nk2 ) then
            go to 705
          end if

          if ( nk2 .lt. nk1 ) then
            max = nk1
            min = nk2
          end if

          ind = p(max) + min
          nc = cc(ind)
          d = float ( nc ) - y1(nb1) - y1(nb2)
          d = d - y2(nk1) - y2(nk2)

          if ( eps .lt. abs ( d ) ) then
            write ( *, '(a)' ) ' '
            write ( *, '(a)' ) 'SMP - Fatal error!'
            write ( *, '(a)' ) '  Optimality conditions violated.'
            write ( *, '(a,i4,a,i4)' ) 
     &        '  Edge NK1 = ', nk1, ' NK2 = ', nk2
            stop
          end if

          if ( nb1 .ne. nb ) then
            go to 735
          end if

740       continue

          if ( nb3 .eq. nb ) then
            go to 705
          end if

745       continue

          n2 = sm(nb3)
          nb2 = basis(n2)
          n3 = sm(nb2)
          min = n2
          max = n3

          if ( n2 .eq. n3 ) then
            go to 705
          end if

          if ( n3 .lt. n2 ) then
            max = n2
            min = n3
          end if

          ind = p(max) + min
          nc = cc(ind)
          d = float ( nc ) - y1(nb2) - y1(nb3)
          d = d - y2(n2) - y2(n3)

          if ( eps .lt. abs ( d ) ) then
            write ( *, '(a)' ) ' '
            write ( *, '(a)' ) 'SMP - Fatal error!'
            write ( *, '(a)' ) '  Optimality conditions violated.'
            write ( *, '(a,i4,a,i4)' ) 
     &        '  Edge N2 = ', n2, ' N3 = ', n3
            stop
          end if

          go to 740

        end if

      end do

      return
      end
      subroutine ssort ( a, b, n )

c*********************************************************************72
c
cc SSORT sorts an integer vector, and rearranges a second one.
c
c  Modified:
c
c    19 November 2007
c
c  Author:
c
c    Ulrich Derigs
c
c  Reference:
c
c    Rainer Burkard, Ulrich Derigs,
c    Assignment and Matching Problems: Solution Methods with Fortran-Programs,
c    Springer, 1980,
c    ISBN: 3-540-10267-1,
c    LC: QA402.5 B86.
c
c  Parameters:
c
c    Input, integer N, the dimension of the vector.
c
c    Input/output, integer A(N), on input, the vector to sort.
c    On output, the sorted vector.
c
c    Input/output, integer B(N), on input, the values 1 through N in order.
c    On output, the permutation vector.
c
      implicit none

      integer n

      integer a(n)
      integer ah
      integer b(n)
      integer bh
      integer f
      integer i
      integer is
      integer j
      integer js
      integer ls
      integer n2
      integer s
      integer t

      f = 1

      if ( n .le. f ) then
        return
      end if

      n2 = ( n - f + 1 ) / 2
      s = 1023

      do t = 1, 10

        if ( s .gt. n2 ) then
          go to 90
        end if

        ls = n - s

        do i = f, ls

          is = i + s
          ah = a(is)
          bh = b(is)
          j = i
          js = is

5         continue

          if ( ah .ge. a(j) ) then
            go to 10
          end if

          a(js) = a(j)
          b(js) = b(j)
          js = j
          j = j - s

          if ( j .ge. f ) then
            go to 5
          end if

10        continue

          a(js) = ah
          b(js) = bh

        end do

90      continue

        s = s / 2

      end do

      return
      end
      subroutine timestamp ( )

c*********************************************************************72
c
cc TIMESTAMP prints out the current YMDHMS date as a timestamp.
c
c  Discussion:
c
c    This FORTRAN77 version is made available for cases where the
c    FORTRAN90 version cannot be used.
c
c  Modified:
c
c    12 January 2007
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    None
c
      implicit none

      character * ( 8 ) ampm
      integer d
      character * ( 8 ) date
      integer h
      integer m
      integer mm
      character * ( 9 ) month(12)
      integer n
      integer s
      character * ( 10 ) time
      integer y

      save month

      data month /
     &  'January  ', 'February ', 'March    ', 'April    ', 
     &  'May      ', 'June     ', 'July     ', 'August   ', 
     &  'September', 'October  ', 'November ', 'December ' /

      call date_and_time ( date, time )

      read ( date, '(i4,i2,i2)' ) y, m, d
      read ( time, '(i2,i2,i2,1x,i3)' ) h, n, s, mm

      if ( h .lt. 12 ) then
        ampm = 'AM'
      else if ( h .eq. 12 ) then
        if ( n .eq. 0 .and. s .eq. 0 ) then
          ampm = 'Noon'
        else
          ampm = 'PM'
        end if
      else
        h = h - 12
        if ( h .lt. 12 ) then
          ampm = 'PM'
        else if ( h .eq. 12 ) then
          if ( n .eq. 0 .and. s .eq. 0 ) then
            ampm = 'Midnight'
          else
            ampm = 'AM'
          end if
        end if
      end if

      write ( *, 
     &  '(i2,1x,a,1x,i4,2x,i2,a1,i2.2,a1,i2.2,a1,i3.3,1x,a)' ) 
     &  d, month(m), y, h, ':', n, ':', s, '.', mm, ampm

      return
      end
      subroutine tour ( n, nb, kost, index, top, kb, ina, kurs )

c*********************************************************************72
c
cc TOUR determines an Eulerian tour for CPP.
c
c  Modified:
c
c    07 December 2007
c
c  Author:
c
c    W Puetz
c
c  Reference:
c
c    Rainer Burkard, Ulrich Derigs,
c    Assignment and Matching Problems: Solution Methods with Fortran-Programs,
c    Springer, 1980,
c    ISBN: 3-540-10267-1,
c    LC: QA402.5 B86.
c
c  Parameters:
c
c    Input, integer N, the number of nodes.
c
c    Input, integer NB(2*number of edges), cumulated list of node neighbors.
c
c    Output, integer KOST(2*numberof edges), cumulated next node list 
c    of the tour starting at KURS.
c
c    Input, integer INDEX(N+1), index of the first neighbor of node I in
c    the next node neighbor list.
c
c    Unused, integer TOP, a large integer.
c
c    Output, integer KB(N), index of the last neighbor of node I in the
c    next node neighbor list.
c
c    Workspace, integer INA(N).
c
c    Input, integer KURS, the starting node.
c
      implicit none

      integer n

      integer i
      integer ii
      integer ina(n)
      integer index(n+1)
      integer kb(n)
      integer kn
      integer kn1
      integer kn2
      integer kost(*)
      integer kurs
      integer m
      integer nb(*)
      integer nn
      integer top
      
      m = index(n+1)

      if ( kurs .le. 0 .or. n .lt. kurs ) then
        kurs = 1
      end if

      do kn = 1, n
        i = index(kn) - 1
        kb(kn) = i
        ina(kn) = i
      end do

      kn = kurs
 
9503  continue

      i = ina(kn)

9504  continue

      nn = index(kn+1) - 1

9505  continue

      i = i + 1

      if ( nn. lt. i ) then
        go to 9540
      end if

      kn1 = nb(i)

      if ( n .lt. kn1 ) then
        go to 9505
      end if

      if ( 0 .lt. kn1 ) then

        ii = ina(kn1)

9510    continue

        ii = ii + 1

        if ( nb(ii) .ne. kn ) then
          go to 9510
        end if

        nb(ii) = n + 1
        ii = kb(kn1) + 1
        kb(kn1) = ii
        kost(ii) = kn
        ina(kn) = i
        kn = kn1
        go to 9503

      end if

      kn2 = - kn
      kn1 = - kn1
      ii = ina(kn1)

9515  continue

      ii = ii + 1

      if ( nb(ii) .ne. kn2 ) then
        go to 9515
      end if

      nb(ii) = n + 1
      ii = kb(kn1) + 1
      kb(kn1) = ii
      kost(ii) = kn
      ii = kb(kn) + 1
      kb(kn) = ii
      kost(ii) = kn1

      go to 9505

9540  continue

      ina(kn) = m

      do kn = 1, n

        i = ina(kn)
        ii = kb(kn)

        if ( index(kn) .le. ii .and. i .lt. m ) then
          go to 9504
        end if

      end do

      return
      end
      subroutine wegspe ( n, k, nmk, izaehl, jzaehl, a, b, veksum,
     &  vekt, bool, bool1, aspei, bspei, h1 )

c*********************************************************************72
c
cc WEGSPE is used by QAP to order elements from the matrices A and B.
c
c  Modified:
c
c    10 December 2007
c
c  Author:
c
c    T Boenniger,
c    Rainer Burkard,
c    Karl-Heinz Stratmann
c
c  Reference:
c
c    Rainer Burkard, Ulrich Derigs,
c    Assignment and Matching Problems: Solution Methods with Fortran-Programs,
c    Springer, 1980,
c    ISBN: 3-540-10267-1,
c    LC: QA402.5 B86.
c
c  Parameters:
c
c    Input, integer N, ?
c
c    Input, integer K, ?
c
c    Input, integer NMK, ?
c
c    Input, integer IZAEHL, ?
c
c    Input, integer JZAEHL, ?
c
c    Input, integer A(N,N), the distance matrix.
c
c    Input, integer B(N,N), the connection matrix.
c
c    Input, integer VEKSUM(*), ?
c
c    Workspace, integer VEKT(*).
c
c    Input, logical BOOL(N), ?
c
c    Input, logical BOOL1(N), ?
c
c    Output, integer ASPEI(*), the ordered rows of A.
c
c    Output, integer BSPEI(*), the ordered rows of B.
c
c    Workspace, integer H1(N).
c
      implicit none

      integer n

      integer a(n,n)
      integer aspei(*)
      integer b(n,n)
      logical bool(n)
      logical bool1(n)
      integer bspei(*)
      integer h1(n)
      integer i
      integer iz
      integer izaehl
      integer j
      integer j1
      integer j2
      integer jzaehl
      integer k
      logical logi
      integer nmk
      integer nsum1
      integer t
      integer veksum(*)
      integer vekt(*)

      if ( nmk .ne. n ) then
        go to 2900
      end if
c
c  Sort the rows of A (minus diagonal entries), and copy into ASPEI.
c
      j1 = 1
      do i = 1, n
        iz = 0
        do j = 1, n
          if ( i .ne. j ) then
            iz = iz + 1
            vekt(iz) = a(i,j)
          end if
        end do
        call ssort ( vekt, h1, nmk - 1 )
        do j = 1, nmk - 1
          aspei(j1) = vekt(nmk-j)
          j1 = j1 + 1
        end do
      end do
c
c  Sort the rows of B (minus diagonal entries) and copy into BSPEI.
c
      j1 = 1
      do i = 1, n
        iz = 0
        do j = 1, n
          if ( i .ne. j ) then
            iz = iz + 1
            vekt(iz) = b(i,j)
          end if
        end do
        call ssort ( vekt, h1, nmk - 1 )
        do j = 1, nmk - 1
          bspei(j1) = vekt(j)
          j1 = j1 + 1
        end do
      end do

      return
c
c  If NMK is not equal to N...
c
2900  continue

      nsum1 = veksum(k+1)
      j1 = nsum1
      j2 = nsum1 - ( nmk + 1 ) * nmk

      do i = 1, n

        if ( .not. bool(i) ) then

          t = a(i,izaehl)
          logi = .true.

          do j = 1, nmk

            j2 = j2 + 1

            if ( aspei(j2) .ne. t .or. ( .not. logi ) ) then
              j1 = j1 + 1
              aspei(j1) = aspei(j2)
            else
              logi = .false.
            end if

          end do

        else

          if ( i .eq. izaehl ) then
            j2 = j2 + nmk
          end if

        end if

      end do

      j1 = nsum1
      j2 = nsum1 - ( nmk + 1 ) * nmk

      do i = 1, n

        if ( .not. bool1(i) ) then

          t = b(i,jzaehl)
          logi = .true.

          do j = 1, nmk

            j2 = j2 + 1

            if ( bspei(j2) .ne. t .or. ( .not. logi ) ) then
              j1 = j1 + 1
              bspei(j1) = bspei(j2)
            else
              logi = .false.
            end if

          end do

        else

          if ( i .eq. jzaehl ) then
            j2 = j2 + nmk
          end if

        end if

      end do

      return
      end
      subroutine zufall ( n, partpe, menge, seed )

c*********************************************************************72
c
cc ZUFALL determines a random permutation.
c
c  Modified:
c
c    27 November 2007
c
c  Author:
c
c    Ulrich Derigs
c
c  Reference:
c
c    Rainer Burkard, Ulrich Derigs,
c    Assignment and Matching Problems: Solution Methods with Fortran-Programs,
c    Springer, 1980,
c    ISBN: 3-540-10267-1,
c    LC: QA402.5 B86.
c
c  Parameters:
c
c    Input, integer N, the size of the set to be permuted.
c
c    Output, integer PARTPE(N), the permutation.
c
c    Workspace, logical MENGE(N).
c
c    Input/output, integer SEED, a seed for the random number generator.
c
      implicit none

      integer n

      real az
      integer i
      integer iz
      integer izaehl
      integer ize
      integer j
      logical menge(n)
      integer nmk
      integer partpe(n)
      real r4_uniform_01
      real ranf
      integer seed

      do i = 1, n
        menge(i) = .false.
      end do

      nmk = n

      do ize = 1, n

        az = float ( nmk ) * r4_uniform_01 ( seed )
!       az = float ( nmk ) * ranf ( az )
        j = int ( az ) + 1

        if ( 1 .lt. ize ) then
          izaehl = 0
          do iz = 1, n
            if ( .not. menge(iz) ) then
              izaehl = izaehl + 1
            end if
            if ( izaehl .eq. j ) then
              go to 3540
            end if
          end do
        end if

        iz = j

3540    continue

        partpe(ize) = iz
        menge(iz) = .true.
        nmk = n - ize

      end do

      return
      end
