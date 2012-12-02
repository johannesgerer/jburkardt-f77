      program main

c*********************************************************************72
c
cc MAIN is the main program for the BETIS boundary element code.
c
c  Discussion:
c
c    BETIS is a program for the solution of Laplace's equation,
c    using the boundary element method.
c
c    BETIS was developed by Federico Paris and Jose Canas, of
c    the department of elasticity and strength of materials,
c    in the Industrial Engineering school of the University of
c    Seville, Spain.
c
c    The program uses linear continuous elements and Dirichlet, Neumann
c    as well as mixed boundary conditions can be taken into
c    consideration.
c
c    The original version of BETIS was written in 1978.
c
c    Copies of BETIS may be found at
c
c      http://www.esi2.us.es/mmc/erm
c
c  Modified:
c
c    15 December 2007
c
c  Author:
c
c    Federico Paris, Jose Canas,
c
c  Reference:
c
c    Federico Paris, Jose Canas,
c    Boundary Element Method: Fundamentals and Applications,
c    Oxford, 1997,
c    ISBN: 0-19-856543-7,
c    LC: TA347.B69.P34.
c
c  Parameters:
c
c    Local, real A(2), stores the results of the
c    integrations of the kernels along an element.
c
c    Local, real B(2), stores the results of the
c    integrations of the kernels along an element.
c
c    Local, real CARGA(NPMAX), the load vector of the system of equations
c
c    Local, real COEF(NPMAX,NPMAX), the coefficients of the final system of
c    equations to be solved
c
c    Local, real DFI(NPMAX,2), the fluxes at the nodes of the boundary. It has
c    double the dimension of FI because two fluxes at each node
c    are stored.
c
c    Local, real FI(NPMAX), the potential values at the nodes of the boundary.
c
c    Local, real FINT(NPMAX-1), values of the potential at the internal points.
c
c    User input, character*60 IENC, the title of the problem.
c
c    User input, integer IGAUSS, the number of points to be used
c    in the Gauss quadrature rule.  A value of 4 is typical.
c    1 <= IGAUSS <= 50 is required.
c
c    Local, integer IMP, the unit number for the output data file.
c
c    Local, integer IPIVO(NPMAX), workspace used for pivoting during
c    the solution of the linear system.
c
c    Local, integer LEC, the unit number for the input data file.
c
c    User input, integer NANU, the index of the node where the potential is
c    going to be considered zero.  This is only for the Neumann problem, the
c    fluxes along the boundary being specified.
c
c    Local, integer NCOD(NPMAX), the codes of the nodes of the boundary used to
c    identify the boundary conditions.
c
c    User input, integer NP, the number of points (nodes) of the boundary where
c    the boundary integral equation is going to be applied. This number
c    represents also the number of elements used in the discretization of the
c    boundary.
c
c    User input, integer NPI, the number of internal points where the
c    potential is going to be calculated.
c
c    Local, integer parameter NPMAX, the maximum number of elements (plus one)
c    that can be used to model the boundary of a certain problem.
c
c    Local, real X(*), the X coordinates of the boundary points
c    used in the discretization.
c
c    Local, real XINT(NPMAX-1), the X coordinates of the internal points where
c    the potential is going to be calculated
c
c    Local, real Y(NPMAX), the Y coordinates of the boundary points
c    used in the discretization.
c
c    Local, real YINT(NPMAX-1), the Y coordinates of the internal points where
c    the potential is going to be calculated
c
      implicit none

      integer npmax
      parameter ( npmax = 101 )

      real a(2)
      real ax
      real ay
      real b(2)
      real carga(npmax)
      real cdiag
      real coef(npmax,npmax)
      real cose
      real den
      real dfi(npmax,2)
      real dfx
      real dfy
      real fi(npmax)
      real fint(npmax-1)
      integer i
      integer ielem
      character*60 ienc
      integer igauss
      integer imp
      integer ipivo(npmax)
      integer j
      integer jj
      integer k
      integer kk
      integer kn
      integer knma
      integer knme
      integer l
      integer lec
      integer na
      character*15 namedat
      character*15 namesal
      integer nanu
      integer nc
      integer ncod(npmax)
      integer nint
      integer nodo
      integer np
      integer npi
      real omeg(50)
      real x(npmax)
      real xi(50)
      real xint(npmax-1)
      real xnue
      real y(npmax)
      real yint(npmax-1)
      real ynue

      lec=10
      imp=11

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'BETIS'
      write ( *, '(a)' ) '  FORTRAN77 version'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  A 2D boundary element program'
      write ( *, '(a)' ) '  by Federico Paris and Jose Canas.'

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'Data File name?'
      read(*,'(a)') namedat
      open(lec,file=namedat,status='old',err=9988)

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'Result File name?'
      read(*,'(a)') namesal
      open(imp,file=namesal,status='new')

      write ( *, '(a)' ) 'Number of Gauss Points?'
      read(*,*) igauss

      call gauss_qn ( igauss, xi, omeg )
c
c  Reading of the basic variables of the program
c
c  The title of the example is written (maximum 50 characters) in the
c  first line of the data file and read as the variable ienc
c
      read(lec,'(a)') ienc
c
c  The second line of the data file includes the number of nodes of
c  the boundary (np), the number of internal points where the
c  potential is going to be calculated (npi) and the node where the
c  potential is going to be forced to be zero (nanu).
c
c  np cannot be zero. npi is zero when no information about
c  the domain is required. nanu is zero when the position of the
c  potential distribution is not undetermined, which is the case
c  of mixed or Dirichlet boundary conditions. In the case of Neumann
c  boundary conditions the variable must include the number of the
c  node where the potential is going to be taken as null.
c
c  The three variables must be written as integers separated by a
c  blank space
c
      read(lec,*) np,npi,nanu
      write(imp,1010) ienc
      write(imp,1001)np,npi,igauss,nanu
c
c  The main arrays are initialized to zero
c
      do i = 1, np
        fi(i)    = 0.0E+00
        dfi(i,1) = 0.0E+00
        dfi(i,2) = 0.0E+00
      end do
c
c  The coordinates of the internal points are read, if it has
c  previously been specified that this is required (npi not equal
c  zero).
c
c  As many lines as indicated by the variable npi must necessarily
c  be included. Each of these lines must have the coordinates x and y
c  of each internal point in real format separated by a blank space.
c
      if ( npi .ne. 0 ) then
        read(lec,*) (xint(i),yint(i),i=1,npi)
      end if
c
c  Reading and generation of the coordinates of the nodes of the
c  boundary.
c
c  Each line of this block must have the number of the node (integer)
c  and its two coordinates x and y (real), separated by a blank space.
c
c  The nodes of the boundary must be numbered anticlockwise,
c  starting at any point. The numbers of the nodes in the subsequent
c  lines must always increase although they are not required to be
c  consecutive.
c  If two lines have non consecutive node numbers the program
c  will automatically generate the coordinates of the intermediate
c  points dividing the straight line between the two specified nodes
c  into segments of equal length. The first node can also be used for
c  automatic generation of the last nodes of the boundary, giving it
c  the number np+1 and being obviously the last line of this block.
c
      l=0

    2 continue

      read(lec,*)i,x(i),y(i)

      if ( 0 < l ) then
        nint = i - l
        ax = ( x(i) - x(l) ) / nint
        ay = ( y(i) - y(l) ) / nint
      end if

    4 continue

      l = l + 1
      if(i-l) 5,6,7
    7 x(l) = x(l-1)+ax
      y(l) = y(l-1)+ay
      go to 4
    6 if(np-i) 8,9,2
    9 x(np+1) = x(1)
      y(np+1) = y(1)
    8 continue
c
c  Reading and generation of the boundary conditions of the problem
c
c  Each line of this block must include the number of the node
c  (integer), the code of the node (integer), the potential (real),
c  the flux before(real) and the flux after the node (real). Never, of
c  course, will the three real variables be known. When a variable is
c  unknown a zero value is placed, the code identifying which is/are
c  the known variable/s.
c  The following list will help with the preparation of the data.
c
c   Code     Variable/s to be specified
c
c     1   flux before the node [dfi(i,1)], flux after [dfi(i,2)]
c
c     2   potential [fi(i)], flux after the node[dfi(i,2)]
c
c     3   potential [fi(i)], flux before the node[dfi(i,1)]
c
c     4   potential [fi(i)], the external normal is discontinuous at i
c
c     5   potential [fi(i)], the external normal is continuous at i
c
c  A linear automatic generation similar to the one used for the nodes
c  coordinates is also employed here for the boundary conditions. It
c  can be used independently of the geometry, because the values of
c  the variables at the intermediate nodes are generated using the
c  number of the nodes and not with its geometrical position.
c
c  The same general features explained for the generation of the
c  coordinates are applicable in this case.
c
      l = 0

   10 continue

      read(lec,*) i,ncod(i),fi(i),dfi(i,1),dfi(i,2)
      if(i-l-1)11,12,13
   13 nint = i - l
      ax = ( fi(i) - fi(l) ) / nint
      ay = ( dfi(i,1) - dfi(l,2)) / nint
   12 l = l + 1
      if(i-l)11,14,15
   15 if ( ncod(i) .le. 3 .and. ncod(i) .ne. 2 ) go to 16
      ncod(l) = 5
      fi(l) = fi(l-1) + ax
      go to 12
   16 ncod(l) = 1
      dfi(l,1) = dfi(l-1,2) + ay
      dfi(l,2) = dfi(l,1)
      go to 12
   14 continue
      if(np-i)17,18,10
c
c  The first node is also used as the last plus one to simplify later
c  instructions.
c
   18 fi(np+1) = fi(1)
      dfi(np+1,1) = dfi(1,1)
      dfi(np+1,2) = dfi(1,2)
      ncod(np+1) = ncod(1)
   17 continue
c
c  Writing of input data in order to check them.
c
   19 write(imp,1002)
      write(imp,1003) (i,x(i),y(i),i=1,np)
      write(imp,1102)
      write(imp,1103) (i,ncod(i),fi(i),(dfi(i,k),k=1,2),i=1,np)
c
c   Starting of the application of the Boundary Element Method.
c
c  The variable NODO represents the node of the boundary where the
c  integral equation is being generated.
c
      do nodo = 1, np
c
c  The row of the coefficient matrix corresponding to the equation
c  to be generated and the corresponding position in the load vector
c  are initialized.
c
c  The variable CDIAG, also initialized, is used to store the sum of
c  the A coefficients corresponding to the point where the integral
c  equation is being applied.  In this way, the calculation of the free
c  term is avoided.
c
        do i = 1, np
          coef(nodo,i) = 0.0E+00
        end do

        cdiag = 0.0E+00
        carga(nodo) = 0.0E+00
c
c  The variable IELEM represents the element along which the
c  integration is going to be performed, from the node NODO.
c
        do ielem = 1, np

          if(nodo.eq.1.and.ielem.eq.np) go to 106
          if(ielem.eq.nodo)go to 105
          if((ielem+1).eq.nodo)go to 106
c
c  The node where the integral equation is applied does not belong
c  to the element along which the integration is performed, the
c  integrations being performed numerically in the subroutine HGNUM.
c
          call hgnum ( x(nodo), y(nodo), x(ielem), y(ielem),
     &      x(ielem+1), y(ielem+1), a(1), a(2), b(1), b(2), igauss,
     &      xi, omeg )
          go to 160
c
c  The node where the integral equation is applied belongs to the
c  element along which the integration is performed, an analytical
c  expression being used in the subroutine HGANA.
c
  105     call hgana ( x(nodo), y(nodo), x(ielem+1), y(ielem+1),
     &      a(1), a(2), b(1), b(2) )
          go to 160
  106     call hgana ( x(ielem), y(ielem), x(nodo), y(nodo), a(1),
     &      a(2), b(1), b(2) )
          ax = b(1)
          b(1) = b(2)
          b(2) = ax
c
c  The values of the coefficients A are stored in CDIAG in order to
c  calculate the free term, eqn 3.5.13, Section 3.5.2.3.
c
  160     cdiag = cdiag+a(1)+a(2)
c
c  The values of the integrations performed are going to be introduced
c  into the system of equations according to the codes of the nodes
c  of the elements along which the integrations have been performed.
c
c  The variable K is used here to refer to the initial (K=1) and final
c  node (K=2) of the element along which the integration has just been
c  performed.
c
          do k = 1, 2

            kn = ielem + k - 1
            if ( kn .eq. np + 1 ) then
              kn = 1
            end if
c
c  The system of equations is built up according to the boundary
c  conditions which are reflected in the codes.
c
c  Code 1: the values of the derivative of the function before and
c  after the node are known, case 1, Section 3.5.3.
c
            if ( ncod(kn) .eq. 1 ) then

              carga(nodo) = carga(nodo)+b(k)*dfi(kn,3-k)
              coef(nodo,kn) = coef(nodo,kn)+a(k)
c
c  Code 2: the value of the potential at the node and its derivative
c  after the node are known, case 2a, Section 3.5.3.
c
            else if ( ncod(kn) .eq. 2 ) then

              if ( k .eq. 2 ) then
                carga(nodo) = carga(nodo) - a(k) * fi(kn)
                coef(nodo,kn) = coef(nodo,kn) - b(k)
              else
                carga(nodo) = carga(nodo)+b(k)*dfi(kn,3-k)-a(k)*fi(kn)
              end if
c
c  Code 3: the value of the potential at the node and its derivative
c  before the node are known, case 2b, Section 3.5.3.
c
            else if ( ncod(kn) .eq. 3 ) then

              if ( k .eq. 2 ) then
                carga(nodo) = carga(nodo)+b(k)*dfi(kn,3-k)-a(k)*fi(kn)
              else
                carga(nodo) = carga(nodo)-a(k)*fi(kn)
                coef(nodo,kn) = coef(nodo,kn)-b(k)
              end if
c
c  Code 4: the value of the function at the node is known the normal
c  at the node being discontinuous. The approach of the gradient is
c  employed, case 3b, Section 3.5.3.
c
            else if ( ncod(kn) .eq. 4 ) then

              den = sqrt ( ( x(ielem) - x(ielem+1) )**2
     &                   + ( y(ielem) - y(ielem+1) )**2 )
              xnue = (y(ielem+1)-y(ielem))/den
              ynue = (x(ielem)-x(ielem+1))/den
              knma = kn+1
              knme = kn-1
              if(kn.eq.1)knme = np
              if(kn.eq.np)knma = 1
              den = (x(knma)-x(kn))*(y(knme)-y(kn))-(x(knme)-x(kn))*
     &          (y(knma)-y(kn))
c
c  DFX and DFY represent the components of the gradient according
c  to a linear assumption of the evolution of the function.
c
              dfx = ((y(kn)-y(knma))*fi(knme)+(y(knma)-y(knme))*fi(kn)+
     &          (y(knme)-y(kn))*fi(knma))/den
              dfy = ((x(knma)-x(kn))*fi(knme)+(x(knme)-x(knma))*fi(kn)+
     &          (x(kn)-x(knme))*fi(knma))/den
c
c  If the potential is constant in the neighborhood of the node under
c  study the direction of the gradient is non-determined. One is
c  arbitrarily chosen.
c
              if ( dfx .ne. 0.0E+00 .or. dfy .ne. 0.0E+00 ) then
                cose = (dfx*xnue+dfy*ynue)/sqrt(dfx**2+dfy**2)
              else
                cose = 0.7071E+00
              end if
c
c  The variable cose represents the constants D that appear in eqn
c  3.5.19. It is associated to the cosine of the angle that forms the
c  gradient whose components are DFX and DFY, with the flux at the
c  node of the element under consideration.
c  The variable cose is stored in DFI in order to avoid recalculating
c  after the solution of the system of equations.
c
              kk = 3-k
              dfi(kn,kk) = cose
              coef(nodo,kn) = coef(nodo,kn)-b(k)*cose
              carga(nodo) = carga(nodo)-a(k)*fi(kn)
c
c  Code 5: the value of the function at the node is known, the normal
c  being continuous at the node, case 3a, Section 3.5.3.
c  Code 5 is programmed using the same sentences as code 4, the
c  artificial variable COSE receiving the value 1.
c
            else if ( ncod(kn) .eq. 5 ) then

              cose = 1.0E+00
              kk = 3-k
              dfi(kn,kk) = cose
              coef(nodo,kn) = coef(nodo,kn)-b(k)*cose
              carga(nodo) = carga(nodo)-a(k)*fi(kn)

            end if

          end do
c
c  End of the integrations along all the elements from the node nodo.
c
        end do
c
c  The value of the free term is placed in its correct position
c  according to the code of the node.
c
        if ( ncod(nodo) .ne. 1 ) then
          carga(nodo) = carga(nodo) + cdiag * fi(nodo)
        else
          coef(nodo,nodo) = coef(nodo,nodo)-cdiag
        end if
c
c  End of the generation of the integral equations.
c
      end do
c
c  If some potential is specified along the boundary the system is
c  ready to be solved.  If all the derivatives are known along the
c  boundary the value of the potential must be forced to be zero at
c  one point (node nanu), equation 3.4.16, Section 3.4.3.
c
      na = 0

      if ( nanu .ne. 0 ) then

        na = 1
        do i = 1, np
          coef(np+1,i) = 0.0E+00
        end do

        coef(np+1,nanu) = 1.0E+00
        carga(np+1) = 0.0E+00

      end if
c
c  The system of equations is now ready to be solved.
c
  270 continue

      call pivo ( coef, np, carga, ipivo, npmax, imp, na )
c
c  The vector carga stores the solution of the system of equations,
c  which must now be placed in its correct position.
c
      do nodo = 1, np

        nc = ncod(nodo)
c
c  Code 1: CARGA stores the potential
c
        if ( nc .eq. 1 ) then

          fi(nodo) = carga(nodo)
c
c  Code 2 and Code 3: CARGA stores the value of the derivative
c  of the function before and after the node respectively.
c
        else if ( nc .eq. 2 .or. nc .eq. 3 ) then

          nc = nc - 1
          dfi(nodo,nc) = carga(nodo)
c
c  Code 4: CARGA stores the value of the gradient.  The derivatives are
c  obtained multiplying the gradient by the cosine of the angle that
c  it forms with the corresponding normals, eqns 3.5.18, Section
c  3.5.3. The variable COSE was stored in the variable DFI.
c
c  If the gradient forms 90 degrees with the normal the flux is zero
c  and coincides with the value of the cosine previously stored.
c
        else if ( nc .eq. 4 ) then

          do i=1,2

            if ( 1.0E-10 .lt. abs ( dfi(nodo,i) ) ) then
              dfi(nodo,i) = carga(nodo) * dfi(nodo,i)
            end if

          end do
c
c  Code 5: CARGA stores the value of the derivative of the potential
c  before and after the node.
c
        else if ( nc .eq. 5 ) then

          do i = 1, 2
            dfi(nodo,i) = carga(nodo)
          end do

        end if

      end do
c
c  The problem has been solved in the boundary, the vectors FI and DFI
c  containing the values of the potential and its external derivatives
c  along the boundary.
c
c  Now it is possible to evaluate the function inside the domain.  The
c  integral eqn 2.5.10, Section 2.5, discretized in the form 3.5.20,
c  Section 3.5.4, is applied to any point inside the domain previously
c  specified.
c
      do i = 1, npi

        fint(i) = 0.0E+00
c
c  Integrations along the boundary are performed from the point of
c  interest.  Numerical integration is now always used.
c
         do j = 1, np

           call hgnum ( xint(i),yint(i),x(j),y(j),x(j+1),y(j+1),
     &                a(1),a(2),b(1),b(2),igauss,xi,omeg)

           if ( j .eq. np ) then
             jj = 1
           else
             jj = j + 1
           end if

           fint(i) = fint(i) - a(1)*fi(j)-a(2)*fi(jj)+b(1)*dfi(j,2)+
     &             b(2)*dfi(jj,1)
         end do

         fint(i) = fint(i) / ( 2.0E+00 * 3.1415926E+00 )

       end do
c
c  The problem is finished and the results are going to be printed.
c
      write(imp,1005) ienc
      write(imp,1006) (i,fi(i),(dfi(i,j),j=1,2),i=1,np)

      if ( 0 < npi ) then
        write(imp,1007)
        write(imp,1008) (xint(i),yint(i),fint(i),i=1,npi)
      end if
c
c  Terminate.
c
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'BETIS:'
      write ( *, '(a)' ) '  Normal end of execution.'
      stop
c
c  Wrong data
c
    5 stop 5
   11 stop 11
9988  write(*,*) 'Input File does not exist. '
      stop 12
c
c   Printing formats
c
 1010 format(5x,a60,/,5x,30('**'))
 1001 format(//,5x,'General constants :',/
     @       /,5x,'Number of elements...............',i5,
     @       /,5x,'Number of internal points........',i5,
     @       /,5x,'Number of Gauss points...........',i5,
     @       /,5x,'Node where the potential is null.',i5,//)
 1002 format(5x,'Problem data:',
     @       //,3x,'node',12x,'x coor',9x,'y coor',/)
 1003 format(2x,i5,3x,2f15.4)
 1102 format(/,3x,'node',5x,'code',6x,'potential',4x
     @      ,'flux before',5x,'flux after',/)
 1103 format(2x,i5,3x,i5,1x,3f15.5)
 1005 format(///,5x,a60,/,5x,60('*'),//,5x,
     @       'SOLUTION ON THE BOUNDARY:',
     @       //,5x,'boundary points'//6x,'node'
     @      ,6x,'potential',4x,'flux before ',3x,'flux after',/)
 1006 format(5x,i5,3f15.5)
 1007 format(//,5X,'SOLUTION IN THE DOMAIN:',//,5x,
     @       'internal points',//,4x,'x coor',4x,'y coor',
     @       6x,'potential',/)
 1008 format(2f10.3,f15.3)
      end
      subroutine hgnum ( xp, yp, x1, y1, x2, y2, a1, a2, b1, b2,
     &  igauss, xi, omeg )

c*********************************************************************72
c
cc HGNUM calculates integrals when the point is not in the element.
c
c  Discussion:
c
c    The routine calculates numerically the integrations along
c    the elements of the boundary when the point where the integral
c    equation is applied does not belong to the element, the general
c    scheme represented by eqns 3.5.7, Section 3.5.2.1 being applied.
c
c  Modified:
c
c    15 December 2007
c
c  Author:
c
c    Federico Paris, Jose Canas,
c
c  Reference:
c
c    Federico Paris, Jose Canas,
c    Boundary Element Method: Fundamentals and Applications,
c    Oxford, 1997,
c    ISBN: 0-19-856543-7
c    LC: TA347.B69.P34.
c
c  Parameters:
c
c    Input, real XI(*), the natural coordinates of the Gauss
c    points used in the numerical integration.
c
c    Input, real OMEG(*), the weights associated with such points.
c
      implicit none

      real a1
      real a2
      real ax
      real ay
      real b1
      real b2
      real bx
      real by
      real dis
      real g
      real h
      integer i
      integer igauss
      real omeg(50)
      real r
      real sig
      real x1
      real x2
      real xc
      real xi(50)
      real xp
      real y1
      real y2
      real yc
      real yp

      ax = ( x2 - x1 ) / 2.0E+00
      ay = ( y2 - y1 ) / 2.0E+00
      bx = ( x2 + x1 ) / 2.0E+00
      by = ( y2 + y1 ) / 2.0E+00
c
c  DIS represents the distance from the point where the integral
c  equation is applied to the element along which the integration is
c  being performed.
c
      if ( ax .ne. 0.0E+00 ) then
        dis = abs ( ( ay * xp / ax - yp + y1 - ay * x1 / ax )
     &    / sqrt ( ( ay / ax )**2 + 1.0E+00 ) )
      else
        dis = abs ( xp - x1 )
      end if

      sig = ( x1 - xp ) * ( y2 - yp ) - ( x2 - xp ) * ( y1 - yp )

      if ( sig .lt. 0.0E+00 ) then
        dis = - dis
      end if

      a1 = 0.0E+00
      a2 = 0.0E+00
      b1 = 0.0E+00
      b2 = 0.0E+00
c
c  The kernel of the integral is now evaluated at the Gauss points.
c
      do i = 1, igauss

        xc = ax * xi(i) + bx
        yc = ay * xi(i) + by
c
c  R represents the distance from the point where the integral
c  equation is applied to the Gauss point under consideration.
c
        r = sqrt ( ( xp - xc )**2 + ( yp - yc )**2 )
        h = dis * omeg(i) * sqrt ( ax**2 + ay**2 ) / r**2
        g = alog ( 1.0E+00 / r ) * omeg(i) * sqrt ( ax**2 + ay**2 )
c
c  Equations 3.5.7.
c
        a1 = a1 + ( xi(i) - 1.0E+00 ) * h / 2.0E+00
        a2 = a2 - ( xi(i) + 1.0E+00 ) * h / 2.0E+00
        b1 = b1 - ( xi(i) - 1.0E+00 ) * g / 2.0E+00
        b2 = b2 + ( xi(i) + 1.0E+00 ) * g / 2.0E+00

      end do

      return
      end
      subroutine hgana ( x1, y1, x2, y2, a1, a2, b1, b2 )

c*********************************************************************72
c
cc HGANA calculates integrals along elements adjacent to the node.
c
c  Discussion:
c
c    Subroutine hgana is used to calculate analytically the integrations
c    along the elements adjacent to the node where the integral equation
c    is being applied.
c
c  Modified:
c
c    15 December 2007
c
c  Author:
c
c    Federico Paris, Jose Canas,
c
c  Reference:
c
c    Federico Paris, Jose Canas,
c    Boundary Element Method: Fundamentals and Applications,
c    Oxford, 1997,
c    ISBN: 0-19-856543-7
c    LC: TA347.B69.P34.
c
      implicit none

      real a1
      real a2
      real b1
      real b2
      real dis
      real x1
      real x2
      real y1
      real y2

      a1 = 0.0E+00
      a2 = 0.0E+00
c
c  The expressions 3.5.9 to 3.5.12, Section 3.5.2.2, are applied.
c
      dis = sqrt ( ( x2 - x1 )**2 + ( y2 - y1 )**2 )

      b1 = dis * ( 1.5E+00 - alog ( dis ) ) / 2.0E+00
      b2 = dis * ( 0.5E+00 - alog ( dis ) ) / 2.0E+00

      return
      end
      subroutine pivo ( a, n, c, ipivo, npmax, imp, na )

c*********************************************************************72
c
cc PIVO applies Gauss elimination to solve the linear system.
c
c  Modified:
c
c    15 December 2007
c
c  Author:
c
c    Federico Paris, Jose Canas,
c
c  Reference:
c
c    Federico Paris, Jose Canas,
c    Boundary Element Method: Fundamentals and Applications,
c    Oxford, 1997,
c    ISBN: 0-19-856543-7
c    LC: TA347.B69.P34.
c
      implicit none

      integer npmax

      real a(npmax,npmax)
      real aux
      real c(npmax)
      integer i
      integer imp
      integer ipivo(npmax)
      integer j
      integer k
      integer l
      real mx
      integer n
      integer na
      integer pmx

      do j = 1, n

        mx = a(j,j)
        pmx = j

        do i = j + 1, n + na
          if ( abs ( mx ) .lt. abs ( a(i,j) ) ) then
            mx = a(i,j)
            pmx = i
          end if
        end do

        if ( abs ( mx ) .lt. 1.0E-06 ) then

          write ( imp, '(a)' ) ' '
          write ( imp, '(a)' ) 'PIVO - Fatal error!'
          write ( imp, '(a)' ) '  The matrix is singular.'
          stop 1111

        else

          if ( pmx .ne. j ) then

            k = ipivo(pmx)
            ipivo(pmx) = ipivo(j)
            ipivo(j) = k

            do l = j, n
              aux = a(pmx,l)
              a(pmx,l) = a(j,l)
              a(j,l) = aux
            end do

            aux = c(pmx)
            c(pmx) = c(j)
            c(j) = aux

          end if

          do l = j + 1, n
            a(j,l) = a(j,l) / a(j,j)
          end do

          c(j) = c(j) / a(j,j)
          a(j,j) = 1.0E+00

          do i = 1, n + na
            if ( i .ne. j ) then
              do l = j + 1, n
                a(i,l) = a(i,l) - a(i,j) * a(j,l)
              end do
              c(i) = c(i) - a(i,j) * c(j)
              a(i,j) = 0.0E+00
            end if
          end do

        end if

      end do

      return
      end
      subroutine gauss_qn ( n, ep, om )

c*********************************************************************72
c
cc GAUSS_QN determines a Gauss quadrature rule.
c
c  Modified:
c
c    15 December 2007
c
c  Author:
c
c    Federico Paris, Jose Canas,
c
c  Reference:
c
c    Federico Paris, Jose Canas,
c    Boundary Element Method: Fundamentals and Applications,
c    Oxford, 1997,
c    ISBN: 0-19-856543-7
c    LC: TA347.B69.P34.
c
      implicit none

      real*8 coef(50)
      real ep(50)
      real*8 epsilon(50)
      integer expon
      integer i
      integer ii
      integer ir0
      integer ir1
      integer k
      integer n
      real om(50)
      real*8 omeg(50)
      real*8 root(0:50,2)
      real*8 tol

      tol = 1.0E-14

      if ( n .lt. 2 .or. 50 .lt. n ) then
        stop 1111
      end if
c
c  Initialization.
c
      ir1 = 1
      ir0 = 2
      root(0,1) = -1.0D+00
      root(1,1) =  0.0D+00
      root(2,1) =  1.0D+00
c
c  Computation of polynomial coefficients and roots
c
      do k = 2, n
        do i = 0, n / 2
          call coefic ( k, i, coef(i+1) )
        end do
        do i = 1, k
          call roots ( k, coef, root(i-1,ir1), root(i,ir1),
     &      tol, root(i,ir0) )
        end do
        root(k+1,ir0) =  1.0D+00
        root(0,ir0)   = -1.0D+00
        ii = ir0
        ir0 = ir1
        ir1 = ii
      end do
c
c  Computation of weights
c
      do i = 1, n
        epsilon(i) = root(i,ir1)
        omeg(i) = 0.0D+00
        do k = 0, n / 2
          expon = n - 2 * k
          if ( 1 .lt. expon ) then
            omeg(i) = omeg(i) + coef(k+1)*expon*epsilon(i)**(expon-1)
          else
            omeg(i) = omeg(i) + coef(k+1) * expon
          end if
        end do
        omeg(i) = omeg(i) / 2.0D+00**n
        omeg(i) = 2.0D+00 / ( omeg(i)**2 * ( 1.0D+00 - epsilon(i)**2))

        ep(i) = epsilon(i)
        om(i) = omeg(i)

      end do

      return
      end
      subroutine coefic ( n, k, value )

c*********************************************************************72
c
cc COEFIC determines the coefficients of the polynomial defining the Gauss rule.
c
c  Modified:
c
c    15 December 2007
c
c  Author:
c
c    Federico Paris, Jose Canas,
c
c  Reference:
c
c    Federico Paris, Jose Canas,
c    Boundary Element Method: Fundamentals and Applications,
c    Oxford, 1997,
c    ISBN: 0-19-856543-7
c    LC: TA347.B69.P34.
c
      implicit none

      integer i
      integer inter
      integer k
      integer n
      integer n1
      integer n2
      integer n3
      real*8 sig
      real*8 value

      value = 1.0D+00
      inter = k / 2

      if ( inter * 2 - k .ne. 0 ) then
        sig = -1.0D+00
      else
        sig =  1.0D+00
      end if

      n2 = n - k
      n1 = 2 * n2
      n3 = n - 2 * k

      if ( n1 .eq. 0 ) then
        value = 1.0D+00
      else
        do i = n2 + 1, n1
          value = value * dble ( i )
        end do
      end if

      if ( n3 .ne. 0 ) then
        do i = 1, n3
          value = value / dble ( i )
        end do
      end if

      if ( k .ne. 0 ) then
        do i = 1, k
          value = value / dble ( i )
        end do
      end if

      value = value * sig

      return
      end
      subroutine evalua ( ncoef, coef, x, value )

c*********************************************************************72
c
cc EVALUA evaluates the polynomial defining the Gauss rule.
c
c  Modified:
c
c    15 December 2007
c
c  Author:
c
c    Federico Paris, Jose Canas,
c
c  Reference:
c
c    Federico Paris, Jose Canas,
c    Boundary Element Method: Fundamentals and Applications,
c    Oxford, 1997,
c    ISBN: 0-19-856543-7
c    LC: TA347.B69.P34.
c
      implicit none

      integer ncoef

      real*8 coef(ncoef)
      integer expon
      integer k
      real*8 value
      real*8 x

      value = 0.0D+00

      do k = 0, ncoef / 2
        expon = ncoef - 2 * k
        if ( expon .eq. 0 ) then
          value = value + coef(k+1)
        else
          value = value + coef(k+1) * x**expon
        end if
      end do

      return
      end
      subroutine roots ( ncoef, coef, x00, xff, error, xi )

c*********************************************************************72
c
cc ROOTS seeks roots of the polynomial defining the Gauss rule.
c
c  Modified:
c
c    15 December 2007
c
c  Author:
c
c    Federico Paris, Jose Canas,
c
c  Reference:
c
c    Federico Paris, Jose Canas,
c    Boundary Element Method: Fundamentals and Applications,
c    Oxford, 1997,
c    ISBN: 0-19-856543-7
c    LC: TA347.B69.P34.
c
      implicit none

      integer ncoef

      real*8 coef(ncoef)
      real*8 errac
      real*8 error
      real*8 value0
      real*8 valuef
      real*8 valuei
      real*8 x0
      real*8 x00
      real*8 xf
      real*8 xff
      real*8 xi

      x0 = x00
      xf = xff
      call evalua ( ncoef, coef, x0, value0 )
      call evalua ( ncoef, coef, xf, valuef )
      errac = dabs ( xf - x0 )

      do while ( error .lt. errac )

        xi = ( xf + x0 ) / 2.0D+00

        call evalua ( ncoef, coef, xi, valuei )

        if ( dabs ( valuei ) .lt. error ) then
          return
        end if

        if ( valuei * value0 .lt. 0.0D+00 ) then
          xf = xi
          valuef = valuei
        else
          x0 = xi
          value0 = valuei
        end if

        errac = abs ( xf - x0 )

      end do

      return
      end

