      program main

c*********************************************************************72
c
cc MAIN is the main program for the SERBA boundary element code.
c
c  Discussion:
c
c    PROGRAM SERBA (original version: 1979)
c   
c    PROGRAM FOR THE SOLUTION OF THE ELASTIC PROBLEM  
c    BY THE BOUNDARY ELEMENT METHOD
c   
c    GROUP OF ELASTICITY AND STRENGTH OF MATERIALS
c    INDUSTRIAL ENGINEERING SCHOOL. UNIVERSITY OF SEVILLE
c   
c    The program uses linear continuous elements, and any kind of 
c    combination of boundary conditions in stresses and displacements
c    can be considered.
c
c    The main program is structured in three subroutines: serb1, serb2 
c    and serb3, which correspond to three segments of the original 
c    version of the program, developed when computer memory  
c    was much smaller.
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
      implicit real*8(a-h,o-z)

      character*50 ienc
      common imp,lec,ienc,np,npi,ianul,npmax,e,g,pois,dm,x(101), 
     @y(101),ncod(101),u(101),v(101),sga(101),sgd(101),taua(101),
     @xint(50), yint(50),nanu(5),ug(101),vg(101),coef(205,205),
     @carga(205),taud(101),igauss,isolve,omeg(50),xi(50)

      call serb1
      call gauss_qn (igauss,xi,omeg)
      call serb2 
      call serb3 

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'SERBA'
      write ( *, '(a)' ) '  Normal end of execution.'

      stop 
      end
      subroutine serb1 

c*********************************************************************72
c
cc SERB1 reads input, generates geometry and boundary conditions.
c
c  Discussion:
c
c    This routine reads the input data and generates the geometrical 
c    values and the boundary conditions of the problem.
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
c   list of the main variables used in the program:
c
c   lec:   variable used for the data input file 
c   imp:   variable used for the data output file
c   np:    number of points (nodes) of the boundary where the boundary 
c          integral equation is going to be applied. This number     
c          also represents the number of elements used in the 
c          discretization of the boundary
c   npi:   number of internal points where displacements and stresses 
c          are going to be calculated
c   ianul: number of displacements that are going to be fixed  
c         (independently of the boundary conditions), zero value,  
c          to apply support conditions
c   npmax: maximum dimension (2*number of elements that can be used to 
c          model the boundary maximum + 5 support conditions)
c   e:     elasticity modulus of the material
c   g:     shear modulus of the material
c   pois:  poisson ratio of the material
c   dm:    maximum dimension of the domain used to scale the system 
c             of equations
c
c   list of the main arrays used in this subroutine
c
c   ienc:  array used to store the title of the program
c   x,y:   include the coordinates x and y of the points of the 
c          boundary used in the discretization
c   ncod:  includes the codes of the nodes of the boundary used to
c          identify the boundary conditions
c   u,v:   store displacements in local coordinates of the nodes 
c          of the boundary. They are used as intermediate arrays
c   ug,    store displacements in global coordinates in x and y 
c   vg:    directions respectively
c   sga:   stores the normal stress before the nodes of the boundary
c   sgd:   stores the normal stress after the nodes of the boundary
c   taua:  stores the tangential stress before the nodes of the 
c          boundary
c   taud:  stores the tangential stress after the nodes of the  
c          boundary
c   nanu:  stores the positions of the displacements to be fixed as 
c          support conditions
c   xint,  store the coordinates x,y respectively of the internal
c   yint:  points where displacements and stresses are going to be  
c          calculated
c   coef:  stores the coefficients of the final system of 
c          equations to be solved
c   carga: stores the load vector of the system of equations
c
      implicit real*8(a-h,o-z)

      character*50 ienc
      character*15 namedat
      character*15 namesal

      common imp,lec,ienc,np,npi,ianul,npmax,e,g,pois,dm,x(101), 
     @y(101),ncod(101),u(101),v(101),sga(101),sgd(101),taua(101),
     @xint(50), yint(50),nanu(5),ug(101),vg(101),coef(205,205),
     @carga(205),taud(101),igauss,isolve,omeg(50),xi(50)
c
c   The function dis calculates the distance between two points.
c
      dis(x1,y1,x2,y2) = dsqrt((x1-x2)**2+(y1-y2)**2)
      lec = 10 
      imp = 11
      npmax = 205
      write (*,*) "PROGRAM SERBA ***************"
      write(*,*) "Data File name?"
      read(*,2000) namedat
      write(*,*) "Result File name?"
      read(*,2000) namesal
      open(lec,file = namedat,status = 'old',err = 9988)
      open(imp,file = namesal,status = 'new')
      write(*,*)"Number of Gauss points?"
      read(*,*) igauss
c
c     npara and ncontrol, later defined, are initialized.
c
      npara = 0
      ncontrol = 0
c
c   Reading of the basic variable of the problem.
c
c   The title of the example is written (maximum 50 characters) in the 
c   first line of the data file and read as the variable ienc.
c
      read(lec,2010) ienc
c
c     The second line of the data file includes the values of np, npi 
c     and ianul. np can not be zero. npi is zero when no 
c     information about the domain is required. ianul is zero when the 
c     boundary conditions determine the final position of the 
c     deformed body. When this is not the case the value of ianul 
c     represents the number of components of displacements to be 
c     fixed. The three variables must be written as integers separated 
c     by a blank space.
c
      read(lec,*)np,npi,ianul
c
c     ncontrol is a variable to detect wrong data.
c
      if(np.lt.0.or.np.gt.100) then
        ncontrol = 1
        write(imp,3000)
      end if
      if(npi.lt.0.or.npi.gt.50) then
        ncontrol = 1
        write(imp,3001)
      end if
      if(ianul.lt.0.or.ianul.gt.5) then
        ncontrol = 1
        write(imp,3002)
      end if
c
c     The main arrays are initialized.
c
      do i = 1,np+1 
        u(i) = 0.0 
        v(i) = 0.0 
        ug(i) = 0.0 
        vg(i) = 0.0
        sga(i) = 0.0
        taua(i) = 0.0  
        sgd(i) = 0.0
        taud(i) = 0.0
      end do
c
c     Reading of the properties of the material: e and pois.
c
c     e must have the same units as the stress boundary conditions 
c     and both must be consistent with the units used for the 
c     dimensions.
c
c     The variable ntip represents the kind of problem to be solved:
c            ntip = 0  plane strain 
c            ntip = 1  plane stress. 
c
c     The maximum dimension of the problem dm does not need to be 
c     an accurate value, and is only used for scaling purposes.
c
      read(lec,*)ntip,e,pois,dm
      if(ntip.lt.0.or.ntip.gt.1) then
        ncontrol = 1
        write(imp,3003)
      end if
      if(e.le.0) then
        ncontrol = 1
        write(imp,3004)
      end if
      if(pois.lt.0.or.pois.gt.0.5) then
        ncontrol = 1
        write(imp,3005)
      end if
      if(dm.le.0) then
         ncontrol = 1
         write(imp,3011)
      end if

      g = e/(2*(1+pois))  
      if(ntip.ne.0) pois = pois/(1+pois)
c
c     The coordinates of the internal points are read.
c
c     As many lines as are indicated by the variable npi must 
c     necessarily be included. Each of theses lines must have
c     the coordinates x and y of each internal point in real 
c     format separated by a blank space.
c
       do i = 1,npi
         read(lec,*)xint(i),yint(i)
       end do
c
c   Reading and generation of the coordinates of the nodes of the
c   boundary. Each line of this block must have the number of the node 
c   (integer) and its coordinates x and y (real), separated by a 
c   blank space. The nodes of the boundary must be numbered 
c   anticlockwise (for the interior problem), starting anywhere.
c   The numbers of the nodes in the subsequent lines must always
c   increase although they are not required to be consecutive. If two 
c   lines have non-consecutive node numbers the program will 
c   automatically generate the coordinates of the intermediate points
c   dividing the straight line between the two specified nodes into
c   segments of equal length. The first node can also be used for the
c   automatic generation of the last nodes of the boundary, being
c   given the number np+1 and obviously being the last line of this
c   block.
c
      l = 0 
   20 read(lec,*)i,x(i),y(i)
      if(l.eq.0.and.i.ne.1) then
        write(imp,3006)
        stop 1111
      end if
      if(l)25,30,25 
   25 nint = i-l  
      v1 = (x(i)-x(l))/nint 
      v2 = (y(i)-y(l))/nint 
   30 l = l+1
      if(i-l)35,40,45
   35  write(imp,3006)
      stop 1111
   45 x(l) = x(l-1)+v1  
      y(l) = y(l-1)+v2
      go to 30  
   40 if(np.gt.i) go to 20 
      x(np+1) = x(1)  
      y(np+1) = y(1)
c
c   Reading and generation of the boundary conditions of the problem.
c
c   Each line of this block must include the number of the node 
c   (integer), four integer numbers which identify the four real 
c   variables whose values represent the boundary conditions of the 
c   node. Then (see sentence labeled 50):
c
c     - i represents the number of a node whose boundary conditions
c       are going to be specified.
c       
c     - i1,i2,i3 and i4 (which may vary from 1 to 8) represent an 
c       enter code which identifies the values of the boundary 
c       conditions given by v1,v2,v3,v4.
c       
c     Both the displacements and the components of the stress 
c     vector are going to be specified in local coordinates: the 
c     normal to the boundary and the tangent direction anticlockwise.
c     Accordingly, eight possible variables are associated to a 
c     node, only four of them being known. The codes (values of i1,i2, 
c     i3 and i4) to identify them are:
c
c       1:  normal displacement at the node corresponding to the
c           element before the node
c       2:  tangential displacement at the node corresponding to the
c           element before the node
c       3:  normal displacement at the node corresponding to the
c           element after the node
c       4:  tangential displacement at the node corresponding to the
c           element after the node
c       5:  normal component of the stress vector at the node
c           corresponding to the element before the node
c       6:  tangential component of the stress vector at the node
c           corresponding to the element before the node
c       7:  normal component of the stress vector at the node
c           corresponding to the element after the node
c       8:  tangential component of the stress vector at the node
c           corresponding to the element after the node
c             
c       
c   A single internal code represented by an integer is used 
c   in the program and also in the output data. 
c   The following table summarizes all the possible cases with
c   their corresponding codes. Sometimes, as is shown, just
c   two of the four possible values of the boundary conditions
c   are sufficient to specify them. i1, i2, i3 and i4 must always  
c   be specified in increasing order.
c
c       case          line                                   code
c      -----------------------------------------------------------
c      
c              node  i1   v1   i2   v2   i3   v3   i4   v4
c
c       s-s           5   sb    6    tb   7   sa    8   ta    1
c
c       d-s           1   ub    2    vb   7   sa    8   ta    2
c       
c                     3   ua    4    va   5   sb    6   tb    3
c
c       ds-s          1   ub    6    tb   7   sa    8   ta    4
c       
c                     2   vb    5    sb   7   sa    8   ta    5
c       
c                     3   ua    5    sb   6   tb    8   ta    6
c       
c                     4   va    5    sb   6   tb    7   sa    7
c                                     
c       d-d(smooth)   1   ub    2    vb                       8
c       
c       d-d(corner)   1   ub    2    vb   3   ua    4   va    9
c       
c       d-ds          1   ub    2    vb   3   ua    8   ta    10
c       
c                     1   ub    2    vb   4   va    7   sa    11
c                             
c       ds-d          1   ub    3    ua   4   va    6   tb    12
c       
c                     2   vb    3    ua   4   va    5   sb    13
c
c       ds-ds(smt)    2 vb( = va) 5 sb( = sa)                     14
c             
c                     1 ub( = ua) 6 tb( = ta)                     15
c       
c       ds-ds(corn)   1   ub    3    ua   6   tb    8   ta    16
c       
c                     2   vb    4    va   5   sb    7   sa    17
c                             
c                     1   ub    4    va   6   tb    7   sa    18
c                             
c                     2   vb    3    ua   5   sb    8   ta    19
c                             
c
c     u and v represent the components of the displacements in the 
c     local coordinates defined. s and t represent the component 
c     normal (sigma) and tangential (tau) of the stress tensor.
c     b and a represent before and after the node respectively.
c
c     ds-d represents for instance that one displacement and one 
c     component of stress (ds) are known along one element beside the 
c     node and the displacements (d) are known along the other element 
c     beside the node. 
c
c     The remaining situations have a similar interpretation.
c       
c     Automatic generation of the intermediate boundary 
c     conditions is used when two nodes (whose numbers must
c     always be given in increasing order) are not consecutive.
c     The intermediate values are generated linearly with the
c     number of the nodes and not with their coordinates. The 
c     values after the first node and the values before the 
c     second node must be compatible for the intermediate
c     generation to be possible.
c
      l = 0
   50 read(lec,*)i,i1,v1,i2,v2,i3,v3,i4,v4
      if(l.eq.0.and.i.ne.1) then
        ncontrol = 1
        write(imp,3009)
      end if
      if(i1.eq.5.and.i2.eq.6.and.i3.eq.7.and.i4.eq.8) go to 100
      if(i1.eq.1.and.i2.eq.2.and.i3.eq.7.and.i4.eq.8) go to 200
      if(i1.eq.3.and.i2.eq.4.and.i3.eq.5.and.i4.eq.6) go to 210
      if(i1.eq.1.and.i2.eq.6.and.i3.eq.7.and.i4.eq.8) go to 300
      if(i1.eq.2.and.i2.eq.5.and.i3.eq.7.and.i4.eq.8) go to 310
      if(i1.eq.3.and.i2.eq.5.and.i3.eq.6.and.i4.eq.8) go to 320
      if(i1.eq.4.and.i2.eq.5.and.i3.eq.6.and.i4.eq.7) go to 330
      if(i1.eq.1.and.i2.eq.2.and.i3.eq.0.and.i4.eq.0) go to 400
      if(i1.eq.1.and.i2.eq.2.and.i3.eq.3.and.i4.eq.4) go to 500
      if(i1.eq.1.and.i2.eq.2.and.i3.eq.3.and.i4.eq.8) go to 600
      if(i1.eq.1.and.i2.eq.2.and.i3.eq.4.and.i4.eq.7) go to 610
      if(i1.eq.1.and.i2.eq.3.and.i3.eq.4.and.i4.eq.6) go to 620
      if(i1.eq.2.and.i2.eq.3.and.i3.eq.4.and.i4.eq.5) go to 630
      if(i1.eq.2.and.i2.eq.5.and.i3.eq.0.and.i4.eq.0) go to 700
      if(i1.eq.1.and.i2.eq.6.and.i3.eq.0.and.i4.eq.0) go to 800
      if(i1.eq.1.and.i2.eq.3.and.i3.eq.6.and.i4.eq.8) go to 900
      if(i1.eq.2.and.i2.eq.4.and.i3.eq.5.and.i4.eq.7) go to 910
      if(i1.eq.1.and.i2.eq.4.and.i3.eq.6.and.i4.eq.7) go to 920
      if(i1.eq.2.and.i2.eq.3.and.i3.eq.5.and.i4.eq.8) go to 930
      ncontrol = 1
           write(imp,3007)
c
c   The variable nc is used in all cases to store the code that 
c   would have to be assigned to the intermediate nodes that might
c   require to be generated.
c 
c   Case s-s, code 1 
c 
  100 ncod(i) = 1  
      sga(i) = v1 
      taua(i) = v2  
      sgd(i) = v3 
      taud(i) = v4  
      nc = 1 
      go to 1000
c 
c   Case d-s (the stresses after the node are known), code 2
c 
  200 ncod(i) = 2  
      u(i) = v1 
      v(i) = v2 
      sgd(i) = v3 
      taud(i) = v4  
      nc = 8 
c
c The global components of the displacements at the node are 
c generated.
c
      knme = i-1  
      if(i.eq.1)knme = np 
      den = dis(x(knme),y(knme),x(i),y(i))  
      xnue = (y(i)-y(knme))/den 
      ynue = (x(knme)-x(i))/den 
      ug(i) = xnue*u(i)-ynue*v(i) 
      vg(i) = ynue*u(i)+xnue*v(i)
      go to 1000  
c 
c  Case s-d (the stresses before the node are known), code 3  
c 
  210 ncod(i) = 3  
      u(i) = v1 
      v(i) = v2 
      sga(i) = v3 
      taua(i) = v4 
      nc = 1
c
c The global components of the displacements at the node are 
c generated.
c 
      knma = i+1  
      if(i.eq.np) knma = 1
      if(i.eq.np+1)knma = 2
      den = dis(x(knma),y(knma),x(i),y(i))  
      xnue = (y(knma)-y(i))/den 
      ynue = (x(i)-x(knma))/den 
      ug(i) = xnue*u(i)-ynue*v(i) 
      vg(i) = ynue*u(i)+xnue*v(i)
      go to 1000  
c 
c     cases ds-s, codes 4,5,6,7 
c 
c     code 4
c
  300 ncod(i) = 4  
      u(i) = v1 
      taua(i) = v2  
      sgd(i) = v3 
      taud(i) = v4  
      nc = 15 
      go to 1000  
c
c     Code 5
c
  310 ncod(i) = 5  
      u(i) = v1 
      sga(i) = v2 
      sgd(i) = v3 
      taud(i) = v4  
      nc = 14 
      go to 1000  
c
c     Code 6
c
  320 ncod(i) = 6  
      u(i) = v1 
      sga(i) = v2 
      taua(i) = v3  
      taud(i) = v4  
      nc = 1 
      go to 1000  
c
c     Code 7
c
  330 ncod(i) = 7  
      u(i) = v1 
      sga(i) = v2 
      taua(i) = v3  
      sgd(i) = v4 
      nc = 1 
      go to 1000    
c 
c     Case d-d smooth, code 8 
c 
  400 ncod(i) = 8
      u(i) = v1 
      v(i) = v2
      nc = 8 
c
c The global components of the displacements at the node are 
c generated.
c
      knme = i-1  
      if(i.eq.1)knme = np 
      den = dis(x(knme),y(knme),x(i),y(i))  
      xnue = (y(i)-y(knme))/den 
      ynue = (x(knme)-x(i))/den 
      ug(i) = xnue*u(i)-ynue*v(i) 
      vg(i) = ynue*u(i)+xnue*v(i)
      go to 1000  
c 
c     Case d-d corner, code 9 
c 
  500 ncod(i) = 9  
      u(i) = v1 
      v(i) = v2 
      nc = 8 
c
c The global components of the displacements at the node are 
c generated.
c
      knme = i-1  
      if(i.eq.1)knme = np 
      den = dis(x(knme),y(knme),x(i),y(i))  
      xnue = (y(i)-y(knme))/den 
      ynue = (x(knme)-x(i))/den 
      ug(i) = xnue*u(i)-ynue*v(i) 
      vg(i) = ynue*u(i)+xnue*v(i)
      go to 1000    
c
c     Cases d-ds, codes 10 and 11   
c
c     Code 10
c
  600 taud(i) = v4  
      ncod(i) = 10  
      go to 660
c
c     Code 11
c
  610 sgd(i) = v4 
      ncod(i) = 11  
  660 u(i) = v1 
      v(i) = v2 
      nc = 8 
c
c The global components of the displacements at the node are 
c generated.
c 
      knme = i-1  
      if(i.eq.1)knme = np 
      den = dis(x(knme),y(knme),x(i),y(i))  
      xnue = (y(i)-y(knme))/den 
      ynue = (x(knme)-x(i))/den 
      ug(i) = xnue*u(i)-ynue*v(i) 
      vg(i) = ynue*u(i)+xnue*v(i)
c
c   The local component of the displacement is stored in order to 
c   interpolate this displacement in intermediate nodes.
c
      u(i) = v3
      go to 1000   
c 
c     cases  ds-d, codes 12,13 
c
c     Code 12
c
  620 taua(i) = v4  
      ncod(i) = 12  
      nc = 15 
      go to 670
c
c     Code 13
c
  630 sga(i) = v4 
      ncod(i) = 13  
      nc = 14
c
c   The global components of the displacements at the node are 
c   generated.
c
  670 u(i) = v2 
      v(i) = v3 
      knma = i+1  
      if(i.eq.np) knma = 1
      if(i.eq.np+1)knma = 2
      den = dis(x(knma),y(knma),x(i),y(i))  
      xnue = (y(knma)-y(i))/den 
      ynue = (x(i)-x(knma))/den 
      ug(i) = xnue*u(i)-ynue*v(i) 
      vg(i) = ynue*u(i)+xnue*v(i)
c
c   The local component of the displacement is stored in order to 
c   interpolate this displacement in intermediate nodes.
c
      u(i) = v1
      go to 1000    
c 
c     Cases ds-ds smooth, codes 14,15 
c
c     Code 14
c
  700 ncod(i) = 14  
      sga(i) = v2 
      sgd(i) = v2 
      u(i) = v1 
      nc = 14
      go to 1000  
c
c     Code 15
c
  800 ncod(i) = 15  
      taua(i) = v2  
      taud(i) = v2  
      u(i) = v1 
      nc = 15
      go to 1000
c 
c     Cases ds-ds corner, codes 16,17,18,19 
c 
c
c     Code 16
c
  900 ncod(i) = 16  
      v(i) = v2
      u(i) = v1
      taua(i) = v3  
      taud(i) = v4  
      nc = 15 
      go to 1000
c
c     Code 17
c
  910 ncod(i) = 17 
        v(i) = v2
        u(i) = v1
      sga(i) = v3 
      sgd(i) = v4 
      nc = 14 
      go to 1000
c
c     Code 18
c
  920 ncod(i) = 18  
      v(i) = v2
      u(i) = v1
      taua(i) = v3  
      sgd(i) = v4 
      nc = 15 
      go to 1000
c
c     Code 19
c
  930 ncod(i) = 19  
         v(i) = v2
         u(i) = v1
      sga(i) = v3 
      taud(i) = v4  
      nc = 14 
      go to 1000

 1000 continue
c
c   The values of the functions (stresses and/or displacements)
c   at the intermediate points are generated.
c
c   Detection of internal errors
c
      if(i-l-1)1010,1020,1020
 1010 write(imp,3009)
      stop 1111
 1020 l = l+1 
      if(i-l)1030,1400,1040
 1030  write(imp,3008)
      stop 1111
 1040 continue
c
c   Message to prevent errors when using automatic generation
c   in curved boundaries.
c
      if(npara.eq.0) then
        npara = 1
        write(imp,2020)
      endif
c
c   The boundary conditions of the intermediate nodes are generated
c   according to their codes.
c
      ncod(l) = nc
      if(nc.eq.8) go to 1200
      if(nc.eq.14.or.nc.eq.15) go to 1300
      if(nc.eq.1)go to 1100
      write(imp,3008)
      stop 1111
c 
c      Case s-s 
c 
 1100 sga(l) = sgd(l-1)+(sga(i)-sgd(l-1))/(i-l+1) 
      taua(l) = taud(l-1)+(taua(i)-taud(l-1))/(i-l+1) 
      sgd(l) = sga(l) 
      taud(l) = taua(l) 
      go to 1020  
c 
c     Case d-d
c 
 1200 ug(l) = ug(l-1)+(ug(i)-ug(l-1))/(i-l+1) 
      vg(l) = vg(l-1)+(vg(i)-vg(l-1))/(i-l+1) 
      go to 1020  
c 
c     Case ds-ds  
c 
 1300 u(l) = u(l-1)+(u(i)-u(l-1))/(i-l+1)
c
c     Consideration of the case that a ds-ds node in a corner 
c     is used as a first node for automatic generation.
c
      if(ncod(l-1).le.15) go to 1310
      u(l) = v(l-1)+(u(i)-v(l-1))/(i-l+1)
 1310 continue
c
c     Generation of the stresses
c
      if(nc-14)1315,1320,1330 
 1315 write(imp,3008)
      stop 1111
 1320 sga(l) = sgd(l-1)+(sga(i)-sgd(l-1))/(i-l+1) 
      sgd(l) = sga(l) 
      go to 1020  
 1330 taua(l) = taud(l-1)+(taua(i)-taud(l-1))/(i-l+1) 
      taud(l) = taua(l) 
      go to 1020
c
c   End of sentences of intermediate generation.
c
 1400 if(np-i)1410,1410,50
c
c   The reading of the boundary condition has finished.
c
c   Generation of the boundary conditions of the first node used
c   as fictitious last one.
c
 1410 ncod(np+1) = ncod(1)  
      u(np+1) = u(1)  
      v(np+1) = v(1)  
      ug(np+1) = ug(1)  
      vg(np+1) = vg(1)  
      sga(np+1) = sga(1)  
      taua(np+1) = taua(1)  
      sgd(np+1) = sgd(1)  
      taud(np+1) = taud(1)
c
c  Printing of the generated input data.
c
      write(imp,2030)ienc,np,npi,e,pois,igauss
      if(ntip.eq.1) then
       write(imp,2040)
      else
       write(imp,2050)
      end if
      write(imp,2060)
      write(imp,2070)(i,x(i),y(i),i = 1,np+1)
      write(imp,2080)
      write(imp,2090)(i,ncod(i),u(i),v(i),ug(i),vg(i),i = 1,np+1)
      write(imp,2100)
      write(imp,2110)(i,ncod(i),sga(i),taua(i),sgd(i),
     @     taud(i),i = 1,np+1)
      if(ncontrol.eq.1) then
        write(imp,3010)
        stop 1111
      endif
      return
9988  write(*,*) 'Input File does not exist. Press any key to quit'
            pause
      stop 10000
c
c
c  printing formats
c
 2000 format(a15)
 2010 format(a50)
 2020 format(5x,'WARNING:',/,
     @ 2x,'The displacements and/or stresses at intermediate 
     @points are ',/,
     @ 2x,'generated linearly using the coordinates of the nodes. 
     @Consequently',/,
     @ 2x,'the automatic generation of the boundary conditions can 
     @not be ',/,
     @ 2x,'performed if the boundary is curved. However it can be used 
     @for the',/,
     @ 2x,'case of constant evolution.')
 2030 format(/5x,a50/,5x,25('**'),
     @///,5x,'General constants :',/,
     @ /,5x,'Number of elements........',6x,i3,
     @ /,5x,'Number of internal points.',6x,i3,
     @ /,5x,'Elasticity modulus........',f10.0,
     @ /,5x,'Poisson coefficient.......',9x,f4.3,
     @ /,5x,'Number of Gauss points....',6x,i3)
 2040 format(5x,'Type of problem...........  Plane Stress')
 2050 format(5x,'Type of problem...........  Plane Strain')
 2060 format(/,5x,'node',4x,'x coor',4x,'y coor',/)
 2070 format(4x,i5,2f10.4)
 2080 format(/,1x,'node',1x,' code',7x,'ul',10x,'vl',10x,'ug',10x,'vg'
     @,/)
 2090 format(i5,i7,4e12.4)
 2100 format(/,1x,'node',1x,'  code',6x,'nsb',8x,' tsb',9x,'nsa',8x,
     @' tsa',/)
 2110 format(i5,i7,4e12.4)
 3000 format(5x,'ERROR: ',1x,
     @      'The number of elements used is not allowed')
 3001 format(5x,'ERROR: ',1x,
     @      'The number of internal points used is not allowed')
 3002 format(5x,'ERROR: ',1x,
     @      'The number of displacements to be fixed is not allowed')
 3003 format(5x,'ERROR: ',1x,
     @      'The type of problem used is not valid')
 3004 format(5x,'ERROR: ',1x,
     @      'The Young modulus must be positive')
 3005 format(5x,'ERROR: ',1x,
     @      'The Poisson ratio must be >0 and ² 0.5')
 3006 format(5x,'ERROR: ',1x,'The numbering of the '
     @      'nodes in the generation of the coordinates is wrong')
 3007 format(5x,'ERROR: ',1x,
     @      'The codes used in the boundary conditions are wrong')
 3008 format(5x,'ERROR: ',1x,'Error in Program SERBA')
 3009 format(5x,'ERROR: ',1x,'The numbering of the nodes '
     @      'in the generation of the boundary conditions is wrong')
 3011 format(5x,'ERROR: ',1x,'Maximum dimension must be greater',
     @       ' than zero')
 3010 format(5x,'Wrong data (stop) ***** ')
      end
      subroutine serb2

c*********************************************************************72
c
cc SERB2 builds and solves the system of equations.
c
c  Discussion:
c
c    This routine builds the system of equations and solves it.
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
c  List of the additional main arrays used in this subroutine:
c
c   a,b:   arrays used to store the integration constants
c   ipivo: intermediate array used in the solver
c
      implicit real*8(a-h,o-z)

      real*8 a(2,4)
      real*8 b(2,4)
      character*50 ienc
      integer ipivo(205)

      common imp,lec,ienc,np,npi,ianul,npmax,e,g,pois,dm,x(101), 
     @y(101),ncod(101),u(101),v(101),sga(101),sgd(101),taua(101),
     @xint(50), yint(50),nanu(5),ug(101),vg(101),coef(205,205),
     @carga(205),taud(101),igauss,isolve,omeg(50),xi(50)
c
c     The function dis calculates the distance between two points.
c
      dis(x1,y1,x2,y2) = dsqrt((x1-x2)**2+(y1-y2)**2)          
c
c     Starting the application of the boundary element method.
c
c     The variable "nodo" represents the node of the boundary where 
c     the integral equation is being generated. Two linear equations 
c     corresponding to the node "nodo" used as collocation point are 
c     generated. 
c
      do 1000 nodo = 1,np
c
c     The rows of the matrix coefficient "coef" corresponding to the 
c     equation to be generated and the corresponding positions in the 
c     load vector "carga" are initialized.
c
      do k = 1,2*np 
        coef(nodo,k) = 0.0
        coef(nodo+np,k) = 0.0
      end do
      carga(nodo) = 0.0
      carga(nodo+np) = 0.0
c
c     Although the free term is calculated explicitly, the sum of the 
c     a coefficients is stored for checking purposes during the design 
c     phase of the program. cdigh stores the value corresponding 
c     to the x axis (horizontal) and cdigv the values corresponding to 
c     the y axis (vertical).
c 
      cdigh = 0.0
      cdigv = 0.
      knma = nodo+1 
      knme = nodo-1 
      if(nodo.eq.1) knme = np 
      if(nodo.eq.np) knma = 1 
c 
c     The components of the free term c(i,j) are calculated.  
c 
      den1 = dis(x(knme),y(knme),x(nodo),y(nodo)) 
      xnora = (y(nodo)-y(knme))/den1  
      ynora = (x(knme)-x(nodo))/den1  
      den2 = dis(x(nodo),y(nodo),x(knma),y(knma)) 
      xnord = (y(knma)-y(nodo))/den2  
      ynord = (x(nodo)-x(knma))/den2  
c     
c     al represents the angle subtended by the adjacent segments 
c     to the node "nodo" and cosal and senal its cosinus and sinus. 
c
      cosal = -(xnora*xnord+ynora*ynord)  
      ddd = 1.-cosal**2 
      if(abs(ddd).lt.0.00001)ddd = 0.0
      senal = dsqrt(ddd) 
      if(-ynora*xnord+xnora*ynord.lt.0.)senal = -senal  
      if(abs(cosal).ge.0.000001) then
         al = atan2(senal,cosal)
      else
         nsig = 1
         if(senal.lt.0) nsig = -1
         al = 3.14159265*nsig/2
      end if
      if(al.lt.0.)al = 2*3.14159265+al
c 
c     gam represents the angle between the positive direction of the 
c     axis x and the bisectrix of the solid angle at the node. cogam 
c     and segam are its cosinus and sinus. 
c   
      den = dsqrt((xnora+xnord)**2+(ynora+ynord)**2) 
      segam = (ynora+ynord)/den 
      cogam = (xnora+xnord)/den 
      xnora = (cogam**2-segam**2)*senal/(4*3.14159265*(1-pois))  
c
c     the next three variables correspond to the values of the free 
c     term, cft11 = c(1,1), cft22 = c(2,2), cft12 = c(1,2).
c 
      cft11 = 1-(al/(3.14159265*2))-xnora  
      cft22 = cft11+2*xnora 
      cft12 = -2*segam*cogam*senal/(4*3.14159265*(1-pois))
c
c     the variable "iel" represents the element along which the 
c     integration is going to be performed, from the node "nodo".
c
      do 1600 iel = 1,np
c
c     rlong represents the length of the element.
c
      rlong = dsqrt((y(iel+1)-y(iel))**2 +(x(iel+1)-x(iel))**2)
      xnue = (y(iel+1)-y(iel))/rlong
      ynue = -(x(iel+1)-x(iel))/rlong
c
c     The subroutine numer calculates the integration constants 
c     matrices a(2,4) and b(2,4) numerically if the collocation point
c     does not belong to the element. In the opposite case (the
c     collocation point belongs to the element) the subroutine ana
c     calculates the integration constants analytically.
c
      if(iel.eq.nodo) then
        call ana(rlong,-ynue,xnue,a,b,e,g,pois,dm,1)
        a(1,1) = a(1,1)+1-cft11
        a(1,3) = a(1,3)-cft12
        a(2,1) = a(2,1)-cft12
        a(2,3) = a(2,3)+1-cft22
      else
        if((iel+1.eq.nodo).or.(nodo.eq.1.and.iel.eq.np)) then
          call ana(rlong,-ynue,xnue,a,b,e,g,pois,dm,2)
        else
          call numer(x(nodo),y(nodo),x(iel),y(iel),x(iel+1),y(iel+1),  
     @               a,b,e,g,pois,dm,igauss,xi,omeg)  
        end if
      end if
c
c     The sum of the coefficients of each row is calculated in order
c     to apply for checking purposes a rigid body motion.
c 
      do kk = 1,4 
        cdigh = cdigh+a(1,kk) 
        cdigv = cdigv+a(2,kk)
      end do
c
c     The integration constants corresponding to each node of the 
c     element are going to be stored in the coefficient matrix.
c
c     k = 1 represents the case of the first node of the element,
c     k = 2 represents the case of the second node of the element.
c
      do 1650 k = 1,2
c
c     The variable kn represents the number of the node of the element
c     whose coefficients are going to be stored.
c
      kn = iel+k-1  
      if(kn.eq.np+1) kn = 1
c
c     knma and knme represent respectively the node anterior and 
c     posterior to kn.
c
      knma = kn+1 
      knme = kn-1 
      if(kn.eq.1) knme = np 
      if(kn.eq.np) knma = 1 
c
c     The code of the node of the element whose coefficients are going 
c     to be stored needs to be identified before proceeding to build  
c     up the system of equations.
c
      if(ncod(kn).eq.1)go to 100
      if(ncod(kn).eq.2.or.ncod(kn).eq.3)go to 200
      if(ncod(kn).gt.3.and.ncod(kn).lt.8)go to 300
      if(ncod(kn).eq.8)go to 400
      if(ncod(kn).eq.9)go to 500
      if(ncod(kn).gt.9.and.ncod(kn).lt.14)go to 600
      if(ncod(kn).eq.14)go to 700
      if(ncod(kn).eq.15)go to 800
      if(ncod(kn).gt.15)go to 900
c
  100 continue
c
c     Case s-s, code 1 
c 
      coef(nodo,kn) = coef(nodo,kn)+(k-1)*a(1,2)+(2-k)*a(1,1)
      coef(nodo,kn+np) = coef(nodo,kn+np)+(k-1)*a(1,4)+(2-k)*a(1,3)
      coef(nodo+np,kn) = coef(nodo+np,kn)+(k-1)*a(2,2)+(2-k)*a(2,1)
      coef(nodo+np,kn+np) = coef(nodo+np,kn+np)+(k-1)*a(2,4)
     @                    +(2-k)*a(2,3)
      carga(nodo) = carga(nodo)+(k-1)*(b(1,2)*sga(kn)+b(1,4)*taua(kn))/e
     @            +(2-k)*(b(1,1)*sgd(kn)+b(1,3)*taud(kn))/e  
      carga(nodo+np) = carga(nodo+np)+(k-1)*(b(2,2)*sga(kn)+b(2,4) * 
     @               taua(kn))/e+(2-k)*(b(2,1)*sgd(kn)+b(2,3)*
     @               taud(kn))/e
      go to 1650
c
  200 continue  
c 
c     Cases d-s, codes 2 and 3  
c
      if(ncod(kn).ne.2) go to 210
c 
c     Code 2
c
      carga(nodo) = carga(nodo)-(a(1,2)*(k-1)+a(1,1)*(2-k))*ug(kn)/dm
     @            -(a(1,4)*(k-1)+a(1,3)*(2-k))*vg(kn)/dm+  
     @            (2-k)*(b(1,1)*sgd(kn)+b(1,3)*taud(kn))/e 
      carga(nodo+np) = carga(nodo+np)-(a(2,2)*(k-1)+a(2,1)*(2-k))
     @               *ug(kn)/dm-(a(2,4)*(k-1)+a(2,3)*(2-k))*vg(kn)/dm 
     @               +(2-k)*(b(2,1)*sgd(kn)+b(2,3)*taud(kn))/e  
      coef(nodo,kn) = coef(nodo,kn)-(k-1)*b(1,2) 
      coef(nodo,kn+np) = coef(nodo,kn+np)-(k-1)*b(1,4) 
      coef(nodo+np,kn) = coef(nodo+np,kn)-(k-1)*b(2,2)
      coef(nodo+np,kn+np) = coef(nodo+np,kn+np)-(k-1)*b(2,4)
      go to 1650    
c
c     Code 3
c
  210 carga(nodo) = carga(nodo)-(a(1,2)*(k-1)+a(1,1)*(2-k))*ug(kn)/dm
     @            -(a(1,4)*(k-1)+a(1,3)*(2-k))*vg(kn)/dm+  
     @             (k-1)*(b(1,2)*sga(kn)+b(1,4)*taua(kn))/e 
      carga(nodo+np) = carga(nodo+np)-(a(2,2)*(k-1)+a(2,1)*(2-k))
     @               *ug(kn)/dm-(a(2,4)*(k-1)+a(2,3)*(2-k))*vg(kn)/dm 
     @               +(k-1)*(b(2,2)*sga(kn)+b(2,4)*taua(kn))/e  
      coef(nodo,kn) = coef(nodo,kn)-(2-k)*b(1,1) 
      coef(nodo,kn+np) = coef(nodo,kn+np)-(2-k)*b(1,3) 
      coef(nodo+np,kn) = coef(nodo+np,kn)-(2-k)*b(2,1)
      coef(nodo+np,kn+np) = coef(nodo+np,kn+np)-(2-k)*b(2,3)
      go to 1650
c
  300 continue
c 
c     Cases ds-s, codes 4,5,6 and 7  
c
c     The constants a are modified with the outward normal "nue"
c     to the element,before or after, depending on the code.
c
      if(ncod(kn).ge.6) go to 305
      den = dis(x(kn),y(kn),x(knme),y(knme))  
      xnue = (y(kn)-y(knme))/den
      ynue = -(x(kn)-x(knme))/den       
      go to 310
  305 den = dis(x(knma),y(knma),x(kn),y(kn))  
      xnue = (y(knma)-y(kn))/den 
      ynue = -(x(knma)-x(kn))/den
  310 continue
c
c The coefficients of a are expressed in the local system of 
c reference.
c
      if(k.eq.1)go to 315
      a12m = a(1,2)*xnue+a(1,4)*ynue
      a14m = -a(1,2)*ynue+a(1,4)*xnue
      a22m = a(2,2)*xnue+a(2,4)*ynue
      a24m = -a(2,2)*ynue+a(2,4)*xnue
      go to 320
c
  315 continue
      a11m = a(1,1)*xnue+a(1,3)*ynue
      a13m = -a(1,1)*ynue+a(1,3)*xnue
      a21m = a(2,1)*xnue+a(2,3)*ynue
      a23m = -a(2,1)*ynue+a(2,3)*xnue
c
  320 if(ncod(kn).ne.4) go to 325
c
c     Code 4
c
      coef(nodo,kn) = coef(nodo,kn)+(k-1)*a14m+(2-k)*a13m
      coef(nodo,kn+np) = coef(nodo,kn+np)-(k-1)*b(1,2)
      coef(nodo+np,kn) = coef(nodo+np,kn)+(k-1)*a24m+(2-k)*a23m
      coef(nodo+np,kn+np) = coef(nodo+np,kn+np)-(k-1)*b(2,2)    
      carga(nodo) = carga(nodo)+ (k-1)*(-a12m*u(kn)/dm+b(1,4)*taua(kn)
     @            /e)+(2-k)*(-a11m*u(kn)/dm+b(1,1)*sgd(kn)/e+b(1,3)*
     @            taud(kn)/e)  
      carga(nodo+np) = carga(nodo+np)+(k-1)*(-a22m*u(kn)/dm+b(2,4)*
     @               taua(kn)/e)+(2-k)*(-a21m*u(kn)/dm+b(2,1)*sgd(kn)
     @               /e+b(2,3)*taud(kn)/e)
      go to 1650
c
  325 if(ncod(kn).ne.5) go to 330
c
c     Code 5
c
      coef(nodo,kn) = coef(nodo,kn)+(k-1)*a12m+(2-k)*a11m
      coef(nodo,kn+np) = coef(nodo,kn+np)-(k-1)*b(1,4)
      coef(nodo+np,kn) = coef(nodo+np,kn)+(k-1)*a22m+(2-k)*a21m
      coef(nodo+np,kn+np) = coef(nodo+np,kn+np)-(k-1)*b(2,4)    
      carga(nodo) = carga(nodo)+(k-1)*(-a14m*u(kn)/dm+b(1,2)*sga(kn)/e)
     @            +(2-k)*(-a13m*u(kn)/dm+b(1,1)*sgd(kn)/e+b(1,3)
     @            *taud(kn)/e)       
      carga(nodo+np) = carga(nodo+np)+(k-1)*(-a24m*u(kn)/dm+b(2,2)*
     @               sga(kn)/e)+(2-k)*(-a23m*u(kn)/dm+b(2,1)*sgd(kn)
     @               /e+b(2,3)*taud(kn)/e)
      go to 1650
c
  330 if(ncod(kn).ne.6) go to 335
c
c     Code 6
c
      coef(nodo,kn) = coef(nodo,kn)+(k-1)*a14m+(2-k)*a13m
      coef(nodo,kn+np) = coef(nodo,kn+np)-(2-k)*b(1,1)
      coef(nodo+np,kn) = coef(nodo+np,kn)+(k-1)*a24m+(2-k)*a23m
      coef(nodo+np,kn+np) = coef(nodo+np,kn+np)-(2-k)*b(2,1)    
      carga(nodo) = carga(nodo)+(k-1)*(-a12m*u(kn)/dm+b(1,2)*sga(kn)/
     @            e+b(1,4)*taua(kn)/e)+(2-k)*(-a11m*u(kn)/dm+b(1,3)
     @            *taud(kn)/e)       
      carga(nodo+np) = carga(nodo+np)+(k-1)*(-a22m*u(kn)/dm+b(2,2)*
     @               sga(kn)/e+b(2,4)*taua(kn)/e)+(2-k)*(-a21m*u(kn)
     @               /dm+b(2,3)*taud(kn)/e)
      go to 1650
c
  335 continue
c
c     Code 7
c
      coef(nodo,kn) = coef(nodo,kn)+(k-1)*a12m+(2-k)*a11m
      coef(nodo,kn+np) = coef(nodo,kn+np)-(2-k)*b(1,3)
      coef(nodo+np,kn) = coef(nodo+np,kn)+(k-1)*a22m+(2-k)*a21m
      coef(nodo+np,kn+np) = coef(nodo+np,kn+np)-(2-k)*b(2,3)    
      carga(nodo) = carga(nodo)+(k-1)*(-a14m*u(kn)/dm+b(1,2)*sga(kn)
     @           /e+b(1,4)*taua(kn)/e)+(2-k)*(-a13m*u(kn)/dm+b(1,1)
     @           *sgd(kn)/e)   
      carga(nodo+np) = carga(nodo+np)+(k-1)*(-a24m*u(kn)/dm+b(2,2)
     @               *sga(kn)/e+b(2,4)*taua(kn)/e)+(2-k)*(-a23m*
     @               u(kn)/dm+b(2,1)*sgd(kn)/e)
      go to 1650
c
  400 continue
c 
c     Case d-d,smooth boundary, code 8
c 
      carga(nodo) = carga(nodo)-a(1,k)*ug(kn)/dm-a(1,k+2)*vg(kn)/dm 
      carga(np+nodo) = carga(np+nodo)-a(2,k)*ug(kn)/dm-a(2,k+2)*vg(kn)/ 
     @               dm 
      coef(nodo,kn) = coef(nodo,kn)-b(1,k) 
      coef(nodo,np+kn) = coef(nodo,np+kn)-b(1,k+2)
      coef(np+nodo,kn) = coef(np+nodo,kn)-b(2,k) 
      coef(np+nodo,kn+np) = coef(np+nodo,kn+np)-b(2,k+2)  
      go to 1650
c
  500 continue
c      
c     Case d-d, corner, code 9
c 
c     The subroutine cala calculates the values of the coefficients 
c     that relate the values of the stresses at both sides of the node 
c     with the principal stresses.
c
c     cosig1 and cosig2 relate sigma with the principal stresses I and 
c     II respectively. Similarly for cotan1 and cotan2 with respect 
c     to tau.
c
      call cala(x(knme),y(knme),x(kn),y(kn),x(knma),y(knma),ug(knme), 
     @          vg(knme),ug(kn),vg(kn),ug(knma),vg(knma),cosig1,
     @          cosig2,cotau1,cotau2,k)

      carga(nodo) = carga(nodo)-a(1,k)*ug(kn)/dm-a(1,k+2)*vg(kn)/dm 
      carga(np+nodo) = carga(np+nodo)-a(2,k)*ug(kn)/dm-a(2,k+2)*vg(kn)/ 
     @               dm 
      coef(nodo,kn) = coef(nodo,kn)-b(1,k)*cosig1-b(1,k+2)*cotau1 
      coef(nodo,np+kn) = coef(nodo,np+kn)-b(1,k)*cosig2-b(1,k+2)*cotau2  
      coef(np+nodo,kn) = coef(np+nodo,kn)-b(2,k)*cosig1-b(2,k+2)*cotau1 
      coef(np+nodo,kn+np) = coef(np+nodo,kn+np)-b(2,k)*cosig2-b(2,k+2)* 
     @                    cotau2  
      go to 1650

  600 continue
c 
c        cases d-ds, codes 10, 11, 12 and 13 
c
c   Calculation of the outward normals before (nuea) and 
c   after (nued) the node kn.
c
      den = dis(x(kn),y(kn),x(knme),y(knme))  
      xnuea = (y(kn)-y(knme))/den 
      ynuea = -(x(kn)-x(knme))/den
      den = dis(x(knma),y(knma),x(kn),y(kn))  
      xnued = (y(knma)-y(kn))/den 
      ynued = -(x(knma)-x(kn))/den
c
c   These codes require the use of the Cauchy relation, which
c   will be applied in the form:
c   xa*sga(kn)+ya*taua(kn) = xd*sgd(kn)+yd*taud(kn)
c   the coefficients xa, ya, xd and yd are now calculated.
c
      xa = xnued*xnuea+ynued*ynuea
      ya = -xnued*ynuea+ynued*xnuea
      xd = xa
      yd = -ya
c
c  nsmo is a variable used to designate whether the boundary is
c  smooth or not at the node under consideration. 
c  nsmo eq. 0 (ya = 0) means smooth boundary. 
c  nsmo eq. 1 means non smooth boundary.
c  nsmo eq. 2 is assigned to the particular case of a ninety
c  degree corner, which requires a particular treatment in the case
c  of codes 10 and 12.
c
        nsmo = 0
        if(abs(ya).ge.0.01)nsmo = 1
        if(abs(xa).le.0.01)nsmo = 2

      if(ncod(kn).ne.10)go to 625
c
c  code 10 
c
c  The case of a ninety degree corner (nsmo = 2) can not be considered
c  as a particular case of the general one. It will be considered,
c  together with the same situation corresponding to code 12, later 
c  on.
c
      if(nsmo.ne.2) go to 605
c
c  The continuity of the stress tensor (or the Cauchy relation)
c  leads directly in this case to the knowledge of the tangential 
c  components at both sides of the node under consideration.
c
      taua(kn) = -taud(kn)
      go to 690

  605 continue      
c
c  The Cauchy relation is applied for the remaining cases 
c  corresponding to code 10 in the form:
c  sgd(kn) = xsd1*sga(kn)+xsd2*taua(kn)-xsd3*taud(kn)
c  the coefficients xsd1, xsd2 and xsd3 are now calculated.
c
      if(nsmo.eq.0)go to 610
c
c   The boundary is not smooth and the coefficients of the Cauchy
c   relation are then calculated using the regular form.
c        
      xsd1 = xa/xd
      xsd2 = ya/xd
      xsd3 = yd/xd
      go to 615

  610 continue
c
c   The boundary has, at the node under consideration, a continuous 
c   normal.
c   The normal components of the stress vectors, before and after 
c   the node, are identified. On solving the problem, the remaining 
c   tangential component of the stress vectors is found. This will be 
c   compared with the known component to check whether the stress 
c   tensor is continuous or not at the node under consideration.
c
      xsd1 = 1.
      xsd2 = 0.
      xsd3 = 0.

  615 continue

      coef(nodo,kn) = coef(nodo,kn)-(k-1)*b(1,2)-(2-k)*b(1,1)*xsd1
      coef(nodo,kn+np) = coef(nodo,kn+np)-(k-1)*b(1,4)-(2-k)*b(1,1)*xsd2
      coef(nodo+np,kn) = coef(nodo+np,kn)-(k-1)*b(2,2)-(2-k)*b(2,1)*xsd1
      coef(nodo+np,kn+np) = coef(nodo+np,kn+np)-(k-1)*b(2,4)
     @                    -(2-k)*b(2,1)*xsd2
      carga(nodo) = carga(nodo)+
     @            (k-1)*(-a(1,2)*ug(kn)/dm-a(1,4)*vg(kn)/dm)
     @            +(2-k)*(-a(1,1)*ug(kn)/dm-a(1,3)*vg(kn)/dm+
     @            b(1,3)*taud(kn)/e-b(1,1)*xsd3*taud(kn)/e)
      carga(nodo+np) = carga(nodo+np)+
     @               (k-1)*(-a(2,2)*ug(kn)/dm-a(2,4)*vg(kn)/dm)
     @               +(2-k)*(-a(2,1)*ug(kn)/dm-a(2,3)*vg(kn)/dm+
     @               b(2,3)*taud(kn)/e-b(2,1)*xsd3*taud(kn)/e)

      go to 1650

  625 if(ncod(kn).ne.11) go to 650
c
c   code 11
c
c  The Cauchy relation is applied in the form:
c  taua(kn) = xta1*sgd(kn)+xta2*taud(kn)-xta3*sga(kn)
c  the coefficients xta1, xta2 and xta3 are now calculated.
c
      if(nsmo.eq.0)go to 630
c
c   The boundary is not smooth and the coefficients of the Cauchy
c   relation are then calculated using the regular form.
c        
      xta1 = xd/ya
      xta2 = yd/ya
      xta3 = xa/ya
      go to 635
c
  630 continue
c
c   The boundary has, at the node under consideration, a continuous 
c   normal.
c   The tangential components of the stress vectors, before and after 
c   the node, are identified. On solving the problem, the normal
c   components of the stress vectors are found. These will be compared 
c   with each other to check whether the stress tensor is continuous 
c   or not at the node under consideration.
c
      xta1 = 0.
      xta2 = 1.
      xta3 = 0.

  635 continue

      coef(nodo,kn) = coef(nodo,kn)+(k-1)*(-b(1,2)+b(1,4)*xta3)
      coef(nodo,kn+np) = coef(nodo,kn+np)-(k-1)*b(1,4)*xta2
     @                 -(2-k)*b(1,3)
      coef(nodo+np,kn) = coef(nodo+np,kn)+(k-1)*(-b(2,2)+b(2,4)*xta3)
      coef(nodo+np,kn+np) = coef(nodo+np,kn+np)-(k-1)*b(2,4)*xta2
     @                    -(2-k)*b(2,3)
      carga(nodo) = carga(nodo)+
     @            (k-1)*(-a(1,2)*ug(kn)/dm-a(1,4)*vg(kn)/dm+
     @            b(1,4)*xta1*sgd(kn)/e)
     @            +(2-k)*(-a(1,1)*ug(kn)/dm-a(1,3)*vg(kn)/
     @            dm+b(1,1)*sgd(kn)/e)     
      carga(nodo+np) = carga(nodo+np)+
     @              (k-1)*(-a(2,2)*ug(kn)/dm-a(2,4)*vg(kn)/dm+
     @              b(2,4)*xta1*sgd(kn)/e)
     @              +(2-k)*(-a(2,1)*ug(kn)/dm-a(2,3)*vg(kn)/
     @              dm+b(2,1)*sgd(kn)/e)

        go to 1650

  650 continue

      if(ncod(kn).ne.12)go to 675
c
c   code 12 
c
c  The case of a ninety degree corner can not be considered
c  as a particular case of the general one. It will be considered,
c  later on, together with the same situation corresponding to code 
c  10.
c
      if(nsmo.ne.2) go to 655
c
c  The continuity of the stress tensor (or the Cauchy relation)
c  leads directly in this case to the knowledge of the tangential 
c  components at both sides of the node under consideration.
c
      taud(kn) = -taua(kn)
      go to 690
c
  655 continue      
c
c  The Cauchy relation is applied for the remaining cases
c  corresponding to code 12 in the form:
c  sga(kn) = xsa1*sgd(kn)+xsa2*taud(kn)-xsa3*taua(kn)
c  the coefficients xsa1, xsa2 and xsa3 are now calculated.
c
      if(nsmo.eq.0)go to 660
c
c   The boundary is not smooth and the coefficients of the Cauchy
c   relation are then calculated using the regular form.
c        
      xsa1 = xd/xa
      xsa2 = yd/xa
      xsa3 = ya/xa
      go to 665
c
  660 continue
c
c   The boundary has, at the node under consideration, a continuous 
c   normal.
c   The normal components of the stress vectors, before and after 
c   the node, are identified. On solving the problem, the remaining 
c   tangential component of the stress vectors is found. This will be 
c   compared with the known component to check whether the stress 
c   tensor is continuous or not at the node under consideration.
c
      xsa1 = 1.
      xsa2 = 0.
      xsa3 = 0.

  665 continue

      coef(nodo,kn) = coef(nodo,kn)-(2-k)*b(1,1)-(k-1)*b(1,2)*xsa1
      coef(nodo,kn+np) = coef(nodo,kn+np)-(2-k)*b(1,3)-(k-1)*b(1,2)*xsa2
      coef(nodo+np,kn) = coef(nodo+np,kn)-(2-k)*b(2,1)-(k-1)*b(2,2)*xsa1
      coef(nodo+np,kn+np) = coef(nodo+np,kn+np)-(2-k)*b(2,3)
     @                    -(k-1)*b(2,2)*xsa2
      carga(nodo) = carga(nodo)+
     @            (k-1)*(-a(1,2)*ug(kn)/dm-a(1,4)*vg(kn)/dm+
     @            b(1,4)*taua(kn)/e-b(1,2)*xsa3*taua(kn)/e)
     @            +(2-k)*(-a(1,1)*ug(kn)/dm-a(1,3)*vg(kn)/dm)      
       carga(nodo+np) = carga(nodo+np)+
     @               (k-1)*(-a(2,2)*ug(kn)/dm-a(2,4)*vg(kn)/dm+
     @               b(2,4)*taua(kn)/e-b(2,2)*xsa3*taua(kn)/e)
     @               +(2-k)*(-a(2,1)*ug(kn)/dm-a(2,3)*vg(kn)/dm)

      go to 1650

  675 continue
c
c   code 13
c
c  The Cauchy relation is applied in the form:
c  taud(kn) = xtd1*sga(kn)+xtd2*taua(kn)-xtd3*sgd(kn)
c  the coefficients xta1, xta2 and xta3 are now calculated.
c
      if(nsmo.eq.0)go to 680
c
c   The boundary is not smooth and the coefficients of the Cauchy
c   relation are then calculated.
c
      xtd1 = xa/yd
      xtd2 = ya/yd
      xtd3 = xd/yd
      go to 685

  680 continue
c
c   The boundary has, at the node under consideration, a continuous 
c   normal.
c   The tangential components of the stress vectors before and after 
c   the node are identified. On solving the problem, the normal
c   components of the stress vectors are found, which will be compared 
c   with each other to check whether the stress tensor is continuous 
c   or not at the node under consideration.
c
      xtd1 = 0.
      xtd2 = 1.
      xtd3 = 0.
c
  685 continue
c
      coef(nodo,kn) = coef(nodo,kn)-(k-1)*b(1,4)-(2-k)*b(1,3)*xtd2
      coef(nodo,kn+np) = coef(nodo,kn+np)+(2-k)*(b(1,3)*xtd3-b(1,1))
      coef(nodo+np,kn) = coef(nodo+np,kn)-(k-1)*b(2,4)
     @                 -(2-k)*b(2,3)*xtd2
      coef(nodo+np,kn+np) = coef(nodo+np,kn+np)+(2-k)*(b(2,3)*xtd3
     @                    -b(2,1))
      carga(nodo) = carga(nodo)+
     @            (k-1)*(-a(1,2)*ug(kn)/dm-
     @            a(1,4)*vg(kn)/dm+b(1,2)*sga(kn)/e)
     @            +(2-k)*(-a(1,1)*ug(kn)/dm-a(1,3)*vg(kn)/dm+
     @            b(1,3)*xtd1*sga(kn)/e)   
      carga(nodo+np) = carga(nodo+np)+
     @               (k-1)*(-a(2,2)*ug(kn)/dm-
     @               a(2,4)*vg(kn)/dm+b(2,2)*sga(kn)/e)
     @               +(2-k)*(-a(2,1)*ug(kn)/dm-a(2,3)*vg(kn)/dm
     @               +b(2,3)*xtd1*sga(kn)/e)
c
      go to 1650
c
  690 continue
c
c  Cases of ninety degree corners corresponding to the codes 
c  10 and 12. The tangential stresses have been identified.
c
      coef(nodo,kn) = coef(nodo,kn)-(k-1)*b(1,2)
      coef(nodo,kn+np) = coef(nodo,kn+np)-(2-k)*b(1,1)
      coef(nodo+np,kn) = coef(nodo+np,kn)-(k-1)*b(2,2)
      coef(nodo+np,kn+np) = coef(nodo+np,kn+np)-(2-k)*b(2,1)
      carga(nodo) = carga(nodo)+
     @            (k-1)*(-a(1,2)*ug(kn)/dm-
     @            a(1,4)*vg(kn)/dm+b(1,4)*taua(kn)/e)
     @            +(2-k)*(-a(1,1)*ug(kn)/dm-
     @            a(1,3)*vg(kn)/dm+b(1,3)*taud(kn)/e)  
      carga(nodo+np) = carga(nodo+np)+
     @              (k-1)*(-a(2,2)*ug(kn)/dm
     @              -a(2,4)*vg(kn)/dm+b(2,4)*taua(kn)/e)
     @              +(2-k)*(-a(2,1)*ug(kn)/dm-
     @              a(2,3)*vg(kn)/dm+b(2,3)*taud(kn)/e)

        go to 1650

  700 continue
c
c     Case ds-ds, code 14  
c
      den = dis(x(kn),y(kn),x(knme),y(knme))  
      xnuea = (y(kn)-y(knme))/den 
      ynuea = -(x(kn)-x(knme))/den
      if(k.eq.1) go to 710

      a12m = a(1,2)*xnuea+a(1,4)*ynuea
      a14m = -a(1,2)*ynuea+a(1,4)*xnuea
      a22m = a(2,2)*xnuea+a(2,4)*ynuea
      a24m = -a(2,2)*ynuea+a(2,4)*xnuea
      go to 720

  710 a11m = a(1,1)*xnuea+a(1,3)*ynuea
      a13m = -a(1,1)*ynuea+a(1,3)*xnuea
      a21m = a(2,1)*xnuea+a(2,3)*ynuea
      a23m = -a(2,1)*ynuea+a(2,3)*xnuea

  720 continue

      coef(nodo,kn) = coef(nodo,kn)+(k-1)*a12m+(2-k)*a11m
      coef(nodo,kn+np) = coef(nodo,kn+np)-(k-1)*b(1,4)-(2-k)*b(1,3)
      coef(nodo+np,kn) = coef(nodo+np,kn)+(k-1)*a22m+(2-k)*a21m
      coef(nodo+np,kn+np) = coef(nodo+np,kn+np)-(k-1)*b(2,4)
     @                    -(2-k)*b(2,3)
      carga(nodo) = carga(nodo)+(k-1)*(-a14m*u(kn)/dm+b(1,2)*sga(kn)/e)
     @            +(2-k)*(-a13m*u(kn)/dm+b(1,1)*sga(kn)/e)   
      carga(nodo+np) = carga(nodo+np)+(k-1)*(-a24m*u(kn)/dm+b(2,2) 
     @             *sga(kn)/e)+(2-k)*(-a23m*u(kn)/dm+b(2,1)*sga(kn)/e)
      go to 1650

  800 continue
c 
c     Case ds-ds, code = 15 
c
      den = dis(x(kn),y(kn),x(knme),y(knme))  
      xnuea = (y(kn)-y(knme))/den 
      ynuea = -(x(kn)-x(knme))/den
      if(k.eq.1) go to 810

      a12m = a(1,2)*xnuea+a(1,4)*ynuea
      a14m = -a(1,2)*ynuea+a(1,4)*xnuea
      a22m = a(2,2)*xnuea+a(2,4)*ynuea
      a24m = -a(2,2)*ynuea+a(2,4)*xnuea
      go to 820

  810 a11m = a(1,1)*xnuea+a(1,3)*ynuea
      a13m = -a(1,1)*ynuea+a(1,3)*xnuea
      a21m = a(2,1)*xnuea+a(2,3)*ynuea
      a23m = -a(2,1)*ynuea+a(2,3)*xnuea

  820 continue

      coef(nodo,kn) = coef(nodo,kn)+(k-1)*a14m+(2-k)*a13m
      coef(nodo,kn+np) = coef(nodo,kn+np)-(k-1)*b(1,2)-(2-k)*b(1,1)
      coef(nodo+np,kn) = coef(nodo+np,kn)+(k-1)*a24m+(2-k)*a23m
      coef(nodo+np,kn+np) = coef(nodo+np,kn+np)-(k-1)*b(2,2)
     @                    -(2-k)*b(2,1)
      carga(nodo) = carga(nodo)+(k-1)*(-a12m*u(kn)/dm+b(1,4)*taua(kn)/e)
     @            +(2-k)*(-a11m*u(kn)/dm+b(1,3)*taua(kn)/e)  
      carga(nodo+np) = carga(nodo+np)+(k-1)*(-a22m*u(kn)/dm+b(2,4)* 
     @            taua(kn)/e)+(2-k)*(-a21m*u(kn)/dm+b(2,3)*taua(kn)/e)
      go to 1650
c
  900 continue
c 
c      Case ds-ds, corner, codes 16,17,18 and 19 
c 
      if(ncod(kn).eq.17.or.ncod(kn).eq.19) i1 = 2
      if(ncod(kn).eq.16.or.ncod(kn).eq.18) i1 = 1
      if(ncod(kn).eq.16.or.ncod(kn).eq.19) i2 = 3
      if(ncod(kn).eq.17.or.ncod(kn).eq.18) i2 = 4
c
c      The variables xnuea and ynuea represent the components of the 
c      outward normal before the node.  
c
      den = dis(x(kn),y(kn),x(knme),y(knme))  
      xnuea = (y(kn)-y(knme))/den 
      ynuea = -(x(kn)-x(knme))/den
c
c     The variables cdca1 and cdca2 represent the components of the
c     unit vector of the known displacement before the node.  
c
      cdca1 = xnuea
      cdca2 = ynuea
      if(i1.eq.1) go to 905
c
c  The known displacement before the node coincides with the direction
c  of the element.
c
      cdaux = cdca1  
      cdca1 = -cdca2  
      cdca2 = cdaux
  905 den = dis(x(knma),y(knma),x(kn),y(kn))
c
c   The variables xnued and ynued represent the components of the 
c   outward normal after the node. 
c
      xnued = (y(knma)-y(kn))/den 
      ynued = -(x(knma)-x(kn))/den
c
c   The variables cdcd1 and cdcd2 represent the components of the
c   unit vector of the known displacement after the node.  
c
      cdcd1 = xnued
      cdcd2 = ynued
      if(i2.eq.3) go to 910
c
c   The known displacement after the node coincides with the direction
c   of the element.
c
      cdaux = cdcd1  
      cdcd1 = -cdcd2  
      cdcd2 = cdaux
  910 continue
      prod = dabs(xnuea*xnued+ynuea*ynued)
      if(prod.ge.0.01)go to 915
c
c     Case of a corner of 90 degrees.
c
      if(ncod(kn).ne.16) go to 920
c
c     Case ds-ds 90 degrees code 16
c
      xnue = xnuea
      ynue = ynuea
      ug(kn) = xnue*u(kn)-ynue*v(kn)
      vg(kn) = ynue*u(kn)+xnue*v(kn)
c
      coef(nodo,kn) = coef(nodo,kn)-(k-1)*b(1,2)
      coef(nodo,kn+np) = coef(nodo,kn+np)-(2-k)*b(1,1)
      coef(nodo+np,kn) = coef(nodo+np,kn)-(k-1)*b(2,2)
      coef(nodo+np,kn+np) = coef(nodo+np,kn+np)-(2-k)*b(2,1)
      carga(nodo) = carga(nodo)+(k-1)*(-a(1,2)*ug(kn)/dm-a(1,4)* 
     @            vg(kn)/dm+b(1,4)*taua(kn)/e)+(2-k)*(-a(1,1)* 
     @            ug(kn)/dm-a(1,3)*vg(kn)/dm+b(1,3)*taud(kn)/e)    
       carga(nodo+np) = carga(nodo+np)+(k-1)*(-a(2,2)*ug(kn)/dm-
     @                a(2,4)*vg(kn)/dm+b(2,4)*taua(kn)/e)+(2-k)*
     @               (-a(2,1)*ug(kn)/dm-a(2,3)*vg(kn)/dm
     @               +b(2,3)*taud(kn)/e)    
        go to 1650
c
  920 if(ncod(kn).ne.17) go to 925
c
c     Case ds-ds 90 degrees code 17
c
      xnue = -ynuea
      ynue = xnuea
      ug(kn) = xnue*u(kn)-ynue*v(kn)
      vg(kn) = ynue*u(kn)+xnue*v(kn)
      coef(nodo,kn) = coef(nodo,kn)-(k-1)*b(1,4)
      coef(nodo,kn+np) = coef(nodo,kn+np)-(2-k)*b(1,3)
      coef(nodo+np,kn) = coef(nodo+np,kn)-(k-1)*b(2,4)
      coef(nodo+np,kn+np) = coef(nodo+np,kn+np)-(2-k)*b(2,3)
      carga(nodo) = carga(nodo)+(k-1)*(-a(1,2)*ug(kn)/dm-a(1,4)* 
     @            vg(kn)/dm+b(1,2)*sga(kn)/e)+(2-k)*(-a(1,1)* 
     @            ug(kn)/dm-a(1,3)*vg(kn)/dm+b(1,1)*sgd(kn)/e)     
       carga(nodo+np) = carga(nodo+np)+(k-1)*(-a(2,2)*ug(kn)/dm-
     @                a(2,4)*vg(kn)/dm+b(2,2)*sga(kn)/e)+(2-k)*
     @                (-a(2,1)*ug(kn)/dm-a(2,3)*vg(kn)/
     @                dm+b(2,1)*sgd(kn)/e)        
      go to 1650
c
  925 continue
c
c     The global displacements can not be calculated, a local axis 
c     being required. The matrix of integration constant a has to be 
c     expressed in this system.
c
      if(k.eq.1) go to 930
        a12m = a(1,2)*xnuea+a(1,4)*ynuea
        a14m = -a(1,2)*ynuea+a(1,4)*xnuea
        a22m = a(2,2)*xnuea+a(2,4)*ynuea
        a24m = -a(2,2)*ynuea+a(2,4)*xnuea
        go to 935
c
  930 a11m = a(1,1)*xnued+a(1,3)*ynued
        a13m = -a(1,1)*ynued+a(1,3)*xnued
        a21m = a(2,1)*xnued+a(2,3)*ynued
        a23m = -a(2,1)*ynued+a(2,3)*xnued
c
  935 continue
c
      if(ncod(kn).eq.19) go to 940
c
c     Case ds-ds 90 degrees code 18
c
      coef(nodo,kn) = coef(nodo,kn)+(k-1)*a14m+(2-k)*a11m
      coef(nodo,kn+np) = coef(nodo,kn+np)-(k-1)*b(1,2)
      coef(nodo+np,kn) = coef(nodo+np,kn)+(k-1)*a24m+(2-k)*a21m
      coef(nodo+np,kn+np) = coef(nodo+np,kn+np)-(k-1)*b(2,2)
      carga(nodo) = carga(nodo)+(k-1)*(-a12m*u(kn)/dm+b(1,4)*taua(kn)/e)
     @            +(2-k)*(a13m*u(kn)/dm+b(1,1)*sgd(kn)/e-b(1,3)* 
     @            taua(kn)/e)  
       carga(nodo+np) = carga(nodo+np)+(k-1)*(-a22m*u(kn)/dm+b(2,4)* 
     @                taua(kn)/e)+(2-k)*(a23m*u(kn)/dm+b(2,1) 
     @                *sgd(kn)/e-b(2,3)*taua(kn)/e)
      go to 1650
c
c     Case ds-ds 90 degrees code 19
c
  940 coef(nodo,kn) = coef(nodo,kn)+(k-1)*a12m-(2-k)*a13m
      coef(nodo,kn+np) = coef(nodo,kn+np)-(2-k)*b(1,1)
      coef(nodo+np,kn) = coef(nodo+np,kn)+(k-1)*a22m-(2-k)*a23m
      coef(nodo+np,kn+np) = coef(nodo+np,kn+np)-(2-k)*b(2,1)
      carga(nodo) = carga(nodo)+(k-1)*(-a14m*u(kn)/dm+b(1,2)*sga(kn)/e-
     @            b(1,4)*taud(kn)/e)+(2-k)*(-a11m*u(kn)/dm+b(1,3)* 
     @            taud(kn)/e)  
       carga(nodo+np) = carga(nodo+np)+(k-1)*(-a24m*u(kn)/dm+b(2,2)* 
     @                sga(kn)/e-b(2,4)*taud(kn)/e)+(2-k)*(-a21m*u(kn)/ 
     @                dm+b(2,3)*taud(kn)/e)
      go to 1650
c
c     Case ds-ds of a corner different from 90 degrees
c
c     The components of the global displacements can always be 
c     calculated.
c
  915 ug(kn) = (cdcd2*u(kn)-cdca2*v(kn))/(cdcd2*cdca1-cdcd1*cdca2)  
      vg(kn) = (-cdcd1*u(kn)+cdca1*v(kn))/(cdcd2*cdca1-cdcd1*cdca2)
      if(ncod(kn).ne.16) go to 945
c
c     Case ds-ds corner different from 90 degrees code 16
c
      coef(nodo,kn) = coef(nodo,kn)-(k-1)*b(1,2)
      coef(nodo,kn+np) = coef(nodo,kn+np)-(2-k)*b(1,1)
      coef(nodo+np,kn) = coef(nodo+np,kn)-(k-1)*b(2,2)
      coef(nodo+np,kn+np) = coef(nodo+np,kn+np)-(2-k)*b(2,1)
      carga(nodo) = carga(nodo)+(k-1)*(-a(1,2)*ug(kn)/dm-a(1,4)*
     @            vg(kn)/dm+b(1,4)*taua(kn)/e)+(2-k)*(-a(1,1) 
     @            *ug(kn)/dm-a(1,3)*vg(kn)/dm+b(1,3)*taud(kn)/e)   
      carga(nodo+np) = carga(nodo+np)+(k-1)*(-a(2,2)*ug(kn)/dm-
     @               a(2,4)*vg(kn)/dm+b(2,4)*taua(kn)/e)+(2-k)*
     @               (-a(2,1)*ug(kn)/dm-a(2,3)*vg(kn)/dm+b(2,3)* 
     @               taud(kn)/e)   
      go to 1650

  945 if(ncod(kn).ne.17) go to 950
c
c     Case ds-ds corner different from 90 degrees code 17
c
      coef(nodo,kn) = coef(nodo,kn)-(k-1)*b(1,4)
      coef(nodo,kn+np) = coef(nodo,kn+np)-(2-k)*b(1,3)
      coef(nodo+np,kn) = coef(nodo+np,kn)-(k-1)*b(2,4)
      coef(nodo+np,kn+np) = coef(nodo+np,kn+np)-(2-k)*b(2,3)
      carga(nodo) = carga(nodo)+(k-1)*(-a(1,2)*ug(kn)/dm-a(1,4)* 
     @            vg(kn)/dm+b(1,2)*sga(kn)/e)+(2-k)*(-a(1,1)*
     @            ug(kn)/dm-a(1,3)*vg(kn)/dm+b(1,1)*sgd(kn)/e)     
      carga(nodo+np) = carga(nodo+np)+(k-1)*(-a(2,2)*ug(kn)/dm-
     @               a(2,4)*vg(kn)/dm+b(2,2)*sga(kn)/e)+(2-k)*
     @               (-a(2,1)*ug(kn)/dm-a(2,3)*vg(kn)/dm+b(2,1)* 
     @               sgd(kn)/e)        
      go to 1650
c
  950 if(ncod(kn).ne.18) go to 955
c
c     Case ds-ds corner different from 90 degrees code 18
c
      coef(nodo,kn) = coef(nodo,kn)-(k-1)*b(1,2)
      coef(nodo,kn+np) = coef(nodo,kn+np)-(2-k)*b(1,3)
      coef(nodo+np,kn) = coef(nodo+np,kn)-(k-1)*b(2,2)
      coef(nodo+np,kn+np) = coef(nodo+np,kn+np)-(2-k)*b(2,3)
      carga(nodo) = carga(nodo)+(k-1)*(-a(1,2)*ug(kn)/dm-a(1,4)*
     @            vg(kn)/dm+b(1,4)*taua(kn)/e)+(2-k)*(-a(1,1)*
     @            ug(kn)/dm-a(1,3)*vg(kn)/dm+b(1,1)*sgd(kn)/e)     
       carga(nodo+np) = carga(nodo+np)+(k-1)*(-a(2,2)*ug(kn)/dm-
     @                a(2,4)*vg(kn)/dm+b(2,4)*taua(kn)/e)+(2-k)*
     @                (-a(2,1)*ug(kn)/dm-a(2,3)*vg(kn)/dm+b(2,1)
     @                *sgd(kn)/e) 

        go to 1650
c
c     Case ds-ds corner different from 90 degrees code 19
c
  955 coef(nodo,kn) = coef(nodo,kn)-(k-1)*b(1,4)
      coef(nodo,kn+np) = coef(nodo,kn+np)-(2-k)*b(1,1)
      coef(nodo+np,kn) = coef(nodo+np,kn)-(k-1)*b(2,4)
      coef(nodo+np,kn+np) = coef(nodo+np,kn+np)-(2-k)*b(2,1)
      carga(nodo) = carga(nodo)+(k-1)*(-a(1,2)*ug(kn)/dm-a(1,4)* 
     @            vg(kn)/dm+b(1,2)*sga(kn)/e)+(2-k)*(-a(1,1)*
     @            ug(kn)/dm-a(1,3)*vg(kn)/dm+b(1,3)*taud(kn)/e)    
      carga(nodo+np) = carga(nodo+np)+(k-1)*(-a(2,2)*ug(kn)/dm-a(2,4)
     @               *vg(kn)/dm+b(2,2)*sga(kn)/e)+(2-k)*(-a(2,1)*
     @               ug(kn)/dm-a(2,3)*vg(kn)/dm+b(2,3)*taud(kn)/e)
   
 1650 continue
c
c   End of the storage of the integration constants of the element 
c   "iel" integrating from the node "nodo".
c
 1600 continue  
c
c   End of the integrations along all the elements from the node 
c   "nodo".
c
c     The values of the variables cdigh and cdigv are printed when 
c     they are considered not to be close to zero, this fact 
c     indicating that the equation corresponding to a rigid body 
c     motion is not satisfactorily fulfilled.
c
      if(dabs(cdigh).gt.0.00001) then
          write(imp,3015)  nodo,cdigh
      end if
      if(dabs(cdigv).gt.0.00001) then
      write(imp,3020) nodo,cdigv
      end if

 1000 continue
c
c     End of the generation of the system of equations.
c
c 
c     The support conditions are now to be applied.
c
      if(ianul.eq.0)go to 2000 
      read(lec,*) (nanu(i),i = 1,ianul)
      write(imp,3025) ianul
      write(imp,3030) (nanu(i),i = 1,ianul)
        do 2016 i = 1,ianul
          if(nanu(i).gt.2*np.or.nanu(i).le.0) then
            write(imp,3040)
            stop 1111
          end if
          do j = 1,2*np
            coef(i+2*np,j) = 0.0
          end do
        coef(i+2*np,nanu(i)) = 1.
        carga(i+2*np) = 0.0
2016  continue
2000  continue  
c 
c     The subroutine pivo solves the system of equations.
c
      write(*,*)' Starting the solution of the system of equations'
      call pivo (coef,2*np,carga,ipivo,npmax,imp,ianul)
      write(*,*)' End of solution of the system of equations'
      return
c
c   Printing formats.
c
 3015 format(5x,'WARNING:',/,1x,
     @'The value of the sum of the integration constants at the node',
     @ 2x,i3,/' in the horizontal equation is',1x,e12.5)
 3020 format(5x,'WARNING:',/,1x,
     @'The value of the sum of the integration constants at the node',
     @ 2x,i3,/' in the vertical equation is',1x,e12.5)
 3025 format(/,5x,'The number of displacements to be prevented is'
     @,i5)
 3030 format(/,5x,'The positions of these displacements are'/,
     @(5x,i5))
3040  format(5x,'ERROR:',1x,'The codes given in the support'
     @ ' conditions are out of range')
      end 
      subroutine serb3 

c*********************************************************************72
c
cc SERB3 calculates and prints the displacements and stresses.
c
c  Discussion:
c
c    This routine reorders the results according to the code of each 
c    node, calculates displacements and stresses at internal points and 
c    finally prints the output results.
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
c    list of the additional main arrays used in this subroutine:
c
c    uint,  store the displacements along directions x and y
c    vint : respectively at the internal points
c    sigin: stores the stress tensor at the internal points with the
c           following equivalence:
c           sigin (ip,1):normal stress x at the internal point ip  
c           sigin (ip,2):tangential stress xy at the internal point ip  
c           sigin (ip,3):normal stress y at the internal point ip
c    hh,gg: store the integration constants when calculating the 
c           stresses at internal points
c
      implicit real*8(a-h,o-z)

      dimension a(2,4),b(2,4),uint(50),vint(50),sigin(50,3),
     @          hh(3,4),gg(3,4)
      character*50 ienc
      common imp,lec,ienc,np,npi,ianul,npmax,e,g,pois,dm,x(101), 
     @y(101),ncod(101),u(101),v(101),sga(101),sgd(101),taua(101),
     @xint(50),yint(50),nanu(5),ug(101),vg(101),coef(205,205),
     @carga(205),taud(101),igauss,isolve,omeg(50),xi(50) 
c
      dis(x1,y1,x2,y2) = dsqrt((x1-x2)**2+(y1-y2)**2)
      do 1000 nodo = 1,np  
c
c     The results are reordered according to the code of each node.
c
      if(ncod(nodo).eq.1)go to 100
      if(ncod(nodo).eq.2.or.ncod(nodo).eq.3)go to 200
      if(ncod(nodo).gt.3.and.ncod(nodo).lt.8)go to 300
      if(ncod(nodo).eq.8)go to 400
      if(ncod(nodo).eq.9)go to 500
      if(ncod(nodo).gt.9.and.ncod(nodo).lt.14)go to 600
      if(ncod(nodo).eq.14.or.ncod(nodo).eq.15)go to 700
      if(ncod(nodo).gt.15)go to 900
c 
c     Case s-s, code 1  
c 
  100 ug(nodo) = carga(nodo)*dm 
      vg(nodo) = carga(nodo+np)*dm  
      go to 1000
c
  200 continue
c 
c     Cases d-s, codes 2 and 3  
c 
      if(ncod(nodo).ne.2) go to 210
c
c     Code 2
c
      sga(nodo) = carga(nodo)*e 
      taua(nodo) = carga(nodo+np)*e 
      go to 1000
c
c     Code 3
c
  210 sgd(nodo) = carga(nodo)*e 
      taud(nodo) = carga(nodo+np)*e 
      go to 1000

  300 continue
c 
c     Cases ds-s, codes 4, 5, 6 and 7
c
      knma = nodo+1 
      knme = nodo-1 
      if(nodo.eq.1) knme = np 
      if(nodo.eq.np) knma = 1
c
c     The outward normal to the element along which the displacement
c     is known is calculated in order to be able to calculate the
c     displacements in global coordinates.
c
      if(ncod(nodo).ge.6) go to 310
      den = dis(x(nodo),y(nodo),x(knme),y(knme))  
      xnue = (y(nodo)-y(knme))/den
      ynue = -(x(nodo)-x(knme))/den     
      go to 320

  310 den = dis(x(knma),y(knma),x(nodo),y(nodo))  
      xnue = (y(knma)-y(nodo))/den 
      ynue = -(x(knma)-x(nodo))/den

  320 continue
      if(ncod(nodo).ne.4.and.ncod(nodo).ne.6)go to 330
c
c     The global displacements are calculated, codes 4 and 6.
c
      ug(nodo) = xnue*u(nodo)-ynue*carga(nodo)*dm
      vg(nodo) = ynue*u(nodo)+xnue*carga(nodo)*dm
c
c     Code 4
c
      if(ncod(nodo).eq.4)sga(nodo) = carga(nodo+np)*e
c
c     Code 6
c
      if(ncod(nodo).eq.6)sgd(nodo) = carga(nodo+np)*e
      go to 1000  
c
c     The global displacements are calculated, codes 5 and 7.
c
  330 ug(nodo) = xnue*carga(nodo)*dm-ynue*u(nodo)
      vg(nodo) = ynue*carga(nodo)*dm+xnue*u(nodo)
c
c     Code 5
c
      if(ncod(nodo).eq.5)taua(nodo) = carga(nodo+np)*e
c
c     Code 7
c
      if(ncod(nodo).eq.7)taud(nodo) = carga(nodo+np)*e
      go to 1000

  400 continue
c 
c     Case d-d smooth, code 8             
c 
      sga(nodo) = carga(nodo)*e 
      sgd(nodo) = sga(nodo) 
      taua(nodo) = carga(nodo+np)*e 
      taud(nodo) = taua(nodo) 
      go to 1000

  500 continue
c 
c     Case d-d corner, code 9          
c 
      knma = nodo+1 
      knme = nodo-1 
      if(nodo.eq.1) knme = np 
      if(nodo.eq.np) knma = 1 
c
c     Reordering before the node.
c
      k = 2
      call cala(x(knme),y(knme),x(nodo),y(nodo),x(knma),y(knma),  
     @          ug(knme),vg(knme),ug(nodo),vg(nodo),ug(knma),vg(knma), 
     @          cosig1,cosig2,cotau1,cotau2,k)   
      sga(nodo) = (cosig1*carga(nodo)+cosig2*carga(nodo+np))*e  
      taua(nodo) = (cotau1*carga(nodo)+cotau2*carga(nodo+np))*e
c
c     Reordering after the node.
c
      k = 1
      call cala(x(knme),y(knme),x(nodo),y(nodo),x(knma),y(knma),  
     @          ug(knme),vg(knme),ug(nodo),vg(nodo),ug(knma),vg(knma), 
     @          cosig1,cosig2,cotau1,cotau2,k)   
      sgd(nodo) = (cosig1*carga(nodo)+cosig2*carga(nodo+np))*e  
      taud(nodo) = (cotau1*carga(nodo)+cotau2*carga(nodo+np))*e
      go to 1000

  600 continue
c 
c   cases d-ds codes 10, 11, 12 and 13         
c
      knma = nodo+1 
      knme = nodo-1 
      if(nodo.eq.1) knme = np 
      if(nodo.eq.np) knma = 1
c
c   Calculation of the outward normals before (nuea)
c   and after (nued) the node "nodo".
c
      den = dis(x(nodo),y(nodo),x(knme),y(knme))  
      xnuea = (y(nodo)-y(knme))/den 
      ynuea = -(x(nodo)-x(knme))/den
      den = dis(x(knma),y(knma),x(nodo),y(nodo))  
      xnued = (y(knma)-y(nodo))/den 
      ynued = -(x(knma)-x(nodo))/den
c
c   Calculation of the coefficients of the Cauchy relation.
c
      xa = xnued*xnuea+ynued*ynuea
      ya = -xnued*ynuea+ynued*xnuea
      xd = xa
      yd = -ya
c
c  nsmo is a variable used to designate whether the boundary is
c  smooth or not at the node under consideration. 
c  nsmo eq. zero (ya = 0) means smooth boundary. 
c  nsmo eq. 1 means non smooth boundary.
c  nsmo eq. 2 is assigned to the particular case of a ninety
c  degree corner, which requires a particular treatment in the case
c  of codes 10 and 12.
c
        nsmo = 0
        if(abs(ya).ge.0.01)nsmo = 1
        if(abs(xa).le.0.01)nsmo = 2
c
      if(ncod(nodo).ne.10)go to 625
c
c  code 10 
c
c  The case of a ninety degree corner can not be considered
c  as a particular case of the general one. 
c
      if(nsmo.eq.2) go to 620
c
c   The Cauchy relation is applied in the form:
c   sgd(nodo) = xsd1*sga(nodo)+xsd2*taua(nodo)-xsd3*taud(nodo)
c   the coefficients xsd1, xsd2 and xsd3 are now calculated.
c
      if(nsmo.eq.0) go to 610
c
c   The boundary is not smooth and the coefficients of the Cauchy
c   relation are then calculated in the regular form.
c        
      xsd1 = xa/xd
      xsd2 = ya/xd
      xsd3 = yd/xd
      go to 615

  610 continue
c
c  The boundary is smooth.
c
      xsd1 = 1.
      xsd2 = 0.
      xsd3 = 0.

  615 continue

      sga(nodo) = carga(nodo)*e 
      taua(nodo) = carga(nodo+np)*e 
      sgd(nodo) = xsd1*sga(nodo)+xsd2*taua(nodo)-xsd3*taud(nodo)
c
c   The case of code = 10 with smooth boundary requires a particular
c   checking of the continuity of the stress tensor. The variable
c   nche is used for this purpose.
c
      nche = 0
      if(nsmo.eq.0) nche = 1
      go to 690

  620 continue

      sga(nodo) = carga(nodo)*e 
      sgd(nodo) = carga(nodo+np)*e
      nche = 0
      go to 690

  625 if(ncod(nodo).ne.11) go to 650
c
c   code 11
c
c  The Cauchy relation is applied in the form:
c  taua(nodo) = xta1*sgd(nodo)+xta2*taud(nodo)-xta3*sga(nodo)
c  the coefficients xta1, xta2 and xta3 are now calculated.
c
      if(nsmo.eq.0)go to 630
c
c   The boundary is not smooth and the coefficients of the Cauchy
c   relation are then calculated in the regular form.
c        
        xta1 = xd/ya
        xta2 = yd/ya
        xta3 = xa/ya
        go to 635

  630 continue
c
c  The boundary is smooth.
c
        xta1 = 0.
        xta2 = 1.
        xta3 = 0.

  635 continue

      sga(nodo) = carga(nodo)*e
      taud(nodo) = carga(nodo+np)*e 
      taua(nodo) = xta1*sgd(nodo)+xta2*taud(nodo)-xta3*sga(nodo)
      nche = 0
      go to 690

  650 continue

      if(ncod(nodo).ne.12)go to 675
c
c   code 12
c
c  The case of a ninety degree corner can not be considered
c  as a particular case of the general one. 
c
      if(nsmo.eq.2) go to 670
c
c  The Cauchy relation is applied in the form:
c  sga(nodo) = xsa1*sgd(nodo)+xsa2*taud(nodo)-xsa3*taua(nodo)
c  the coefficients xsa1, xsa2 and xsa3 are now calculated.
c
      if(nsmo.eq.0)go to 660
c
c   The boundary is not smooth and the coefficients of the Cauchy
c   relation are then calculated using the regular form.
c        
      xsa1 = xd/xa
      xsa2 = yd/xa
      xsa3 = ya/xa
      go to 665

  660 continue
c
c  The boundary is smooth.
c
      xsa1 = 1.
      xsa2 = 0.
      xsa3 = 0.

  665 continue

      sgd(nodo) = carga(nodo)*e
      taud(nodo) = carga(nodo+np)*e
      sga(nodo) = xsa1*sgd(nodo)+xsa2*taud(nodo)-xsa3*taua(nodo)
c
c   The case of code = 12 with smooth boundary requires a particular
c   checking of the continuity of the stress tensor. The variable
c   nche is used for this purpose.
c
      nche = 0
      if(nsmo.eq.0) nche = 1
      go to 690

  670 continue

      sga(nodo) = carga(nodo)*e 
      sgd(nodo) = carga(nodo+np)*e
      nche = 0
      go to 690

  675 continue
c
c   code 13
c
c  The Cauchy relation is applied in the form:
c  taud(nodo) = xtd1*sga(nodo)+xtd2*taua(nodo)-xtd3*sgd(nodo)
c  the coefficients xta1, xta2 and xta3 are now calculated.
c
      if(nsmo.eq.0)go to 680
c
c   The boundary is not smooth and the coefficients of the Cauchy
c   relation are then calculated.
c
      xtd1 = xa/yd
      xtd2 = ya/yd
      xtd3 = xd/yd

      go to 685

  680 continue
c
c   The boundary is smooth.
c
      xtd1 = 0.0
      xtd2 = 1.0
      xtd3 = 0.0

  685 continue

      taua(nodo) = carga(nodo)*e
      sgd(nodo) = carga(nodo+np)*e
      taud(nodo) = xtd1*sga(nodo)+xtd2*taua(nodo)-xtd3*sgd(nodo)
      nche = 0

  690 continue
c
c  The cases d-ds (codes 10, 11, 12 and 13) can be associated to
c  cases of a non continuous stress tensor. The Cauchy relation, used
c  to solve these cases, is checked, after having solved the 
c  system of equations.
c 
      sl = xa*sga(nodo)+ya*taua(nodo)
        if(nche.eq.1) sl = taua(nodo)
      sr = xd*sgd(nodo)+yd*taud(nodo)
        if(nche.eq.1) sr = taud(nodo)
      if(abs(sr).lt.1e-5.and.abs(sl).lt.1e-3)go to 1000
      if(abs(sr).lt.1e-3.and.abs(sl).lt.1e-5)go to 1000
      if(abs(sr).lt.1e-5) then
          xqu = sr/sl
      else
          xqu = sl/sr
      end if
      botlim = 0.99
      toplim = 1.01
      if(xqu.gt.botlim.and.xqu.lt.toplim) go to 1000

        write(imp,3010)nodo,nodo

        go to 1000

  700 continue
c
c     Cases ds-ds smooth, codes 14 and 15.
c
      knme = nodo-1
      if((nodo).eq.1)knme = np
      den = dis(x(nodo),y(nodo),x(knme),y(knme))
      xnuea = (y(nodo)-y(knme))/den
      ynuea = (x(knme)-x(nodo))/den
      if(ncod(nodo).eq.15)go to 710
c
c     Code 14
c
      v(nodo) = u(nodo)
      u(nodo) = carga(nodo)*dm
      taua(nodo) = carga(nodo+np)*e
      taud(nodo) = taua(nodo)
      go to 720
c
c     Code 15
c
  710 v(nodo) = carga(nodo)*dm
      sga(nodo) = carga(nodo+np)*e
      sgd(nodo) = sga(nodo)
c
c     Global displacements for codes 14 and 15.
c
  720 ug(nodo) = xnuea*u(nodo)-ynuea*v(nodo)
      vg(nodo) = ynuea*u(nodo)+xnuea*v(nodo)
      go to 1000
c
  900 continue
c 
c     Cases ds-ds corner, cases 16, 17, 18 and 19.  
c
      if (ncod(nodo).ne.16) go to 905
c
c     Code 16
c
      sga(nodo) = carga(nodo)*e 
      sgd(nodo) = carga(nodo+np)*e  
      go to 1000
c
  905 if(ncod(nodo).ne.17) go to 910
c
c     Code 17
c
      taua(nodo) = carga(nodo)*e  
      taud(nodo+np) = carga(nodo+np)*e  
      go to 1000
c
  910 if(ncod(nodo).ne.18) go to 920
c
c     Reordering of code 18 requires it to be known whether the 
c     internal angle is 90 degrees.
c
      knma = nodo+1 
      knme = nodo-1 
      if(nodo.eq.1) knme = np 
      if(nodo.eq.np) knma = 1
c
c     Calculation of the outward normals before and after the node.
c
      den = dis(x(nodo),y(nodo),x(knme),y(knme))
      xnuea = (y(nodo)-y(knme))/den
      ynuea = (x(knme)-x(nodo))/den
      den = dis(x(knma),y(knma),x(nodo),y(nodo))
      xnued = (y(knma)-y(nodo))/den
      ynued = (x(nodo)-x(knma))/den
      prod = dabs(xnuea*xnued+ynuea*ynued)
      if(prod.ge.0.01)go to 915
c
c     Code 18, 90 degrees.
c
      ug(nodo) = -ynuea*carga(nodo)*dm+xnuea*u(nodo)
      vg(nodo) = xnuea*carga(nodo)*dm+ynuea*u(nodo)
      sga(nodo) = carga(nodo+np)*e
      taud(nodo) = -taua(nodo)
      go to 1000
c
c     Code 18, different from 90 degrees.
c
  915 sga(nodo) = carga(nodo)*e 
      taud(nodo) = carga(nodo+np)*e 
      go to 1000
c
c     Reordering code 19 requires it to be known whether the 
c     internal angle is 90 degrees.
c
  920 knma = nodo+1 
      knme = nodo-1 
      if(nodo.eq.1) knme = np 
      if(nodo.eq.np) knma = 1 
c
c     Calculation of the outward normals before and after the node.
c
      den = dis(x(nodo),y(nodo),x(knme),y(knme))
      xnuea = (y(nodo)-y(knme))/den
      ynuea = (x(knme)-x(nodo))/den
      den = dis(x(knma),y(knma),x(nodo),y(nodo))
      xnued = (y(knma)-y(nodo))/den
      ynued = (x(nodo)-x(knma))/den
      prod = dabs(xnuea*xnued+ynuea*ynued)
      if(prod.ge.0.01)go to 925
c
c     Code 19, 90 degrees.
c
      ug(nodo) = xnuea*carga(nodo)*dm-ynuea*u(nodo)
      vg(nodo) = ynuea*carga(nodo)*dm+xnuea*u(nodo)
      sgd(nodo) = carga(nodo+np)*e
      taua(nodo) = -taud(nodo)
      go to 1000
c
c     Code 19, different from 90 degrees.
c
  925 taua(nodo) = carga(nodo)*e  
      sgd(nodo) = carga(nodo+np)*e  
      go to 1000

 1000 continue
c
c     The reordering of the results has been terminated.
c
c     Starting the calculations at internal points.
c
      ug(np+1) = ug(1)  
      vg(np+1) = vg(1)  
      taua(np+1) = taua(1)  
      taud(np+1) = taud(1)  
      sga(np+1) = sga(1)  
      sgd(np+1) = sgd(1)
c
c     Computation of the displacements u and v in global coordinates 
c     for the internal point i.
c
      do i = 1, npi

        uint(i) = 0.0  
        vint(i) = 0.0
c
c     Performing the integrations along the np elements of the 
c     boundary.
c
        do j = 1, np 
c
c     The subroutine lole calculates the values of the integration
c     constants of the integral expression of the displacements
c     integrating from the internal point i along the element j.
c 
          call numer(xint(i),yint(i),x(j),y(j),x(j+1),y(j+1),a,b,e,g,
     &    pois,dm,igauss,xi,omeg)  
c 
c     The values of the displacements at the internal point i are 
c     calculated.
c 
          uint(i) = uint(i)
     &    -a(1,1)*ug(j)-a(1,2)*ug(j+1)-a(1,3)*vg(j)-a(1,4)*  
     &    vg(j+1)+(b(1,1)*sgd(j)+b(1,2)*sga(j+1)+b(1,3)*taud(j)+  
     &    b(1,4)*taua(j+1))*dm/e
 
          vint(i) = vint(i)
     &    -a(2,1)*ug(j)-a(2,2)*ug(j+1)-a(2,3)*vg(j)-a(2,4)*  
     &    vg(j+1)+(b(2,1)*sgd(j)+b(2,2)*sga(j+1)+b(2,3)*taud(j)
     &    +b(2,4)*taua(j+1))*dm/e

        end do
c
c     The integration along the boundary from the internal point i 
c     has been performed.
c
      end do
c
c     The calculation of the displacements at the npi internal points 
c     has been performed.
c
c   Computation of the stress tensor in global coordinates 
c   for the internal point i.
c
      do i = 1, npi

        do j = 1, 3  
          sigin(i,j) = 0.0
        end do
c
c     Performing the integrations along the np elements of the 
c     boundary.
c
        do k = 1, np
c
c     The subroutine piyay calculates the values of the integration
c     constants of the integral expression of the stresses integrating 
c     from the internal point i along the element k.
c   
          call piyay(xint(i),yint(i),x(k),y(k),x(k+1),y(k+1),hh,gg,g,pois,
     &    igauss,xi,omeg)

          do j = 1, 3  
c 
c     The stress tensor at the internal point i is calculated.
c 
            sigin(i,j) = sigin(i,j) 
     &      + hh(j,1)*sgd(k)+hh(j,2)*sga(k+1)+hh(j,3)*  
     @      taud(k)+hh(j,4)*taua(k+1)-gg(j,1)*ug(k)-gg(j,2)*
     @      ug(k+1)-gg(j,3)*vg(k)-gg(j,4)*vg(k+1)
          end do
        end do
c
c     The integration from the internal point i along the boundary 
c     has been finished.
c
      end do
c
c     The calculation of the displacements at the npi internal points 
c     has been finished.
c 
c     Printing the results at the nodes of the boundary.
c 
      write(imp,2010) ienc  
      write(imp,2020)
      write(imp,2030)(i,ug(i),vg(i),i = 1,np)
      write(imp,2040)
      write(imp,2050)(i,sga(i),sgd(i),taua(i),taud(i),i = 1,np)
c
c     Printing the results at the internal points.
c
      write(imp,2060) 
      write(imp,2070)
      write(imp,2080)(xint(i),yint(i),uint(i),vint(i),i = 1,npi) 
      write(imp,2090) 
      write(imp,2100)(xint(i),yint(i),sigin(i,1),sigin(i,3),  
     @               sigin(i,2),i = 1,npi) 

      return
c
c     Printing formats.
c
 2010 format(/,5x,a50,/,5x,25('**'),//,5x,'Solution on the boundary:')
 2020 format(/,5x,'boundary points'//6x,'node',2x,'xdisplacement'
     @,2x,'ydisplacement',/)
 2030 format(5x,i5,2e15.7)
 2040 format(/,6x,'node',2x,'nor stress bf',2x,'nor stress af',2x,
     @'tan stress bf',2x,'tan stress af',/)
 2050 format(5x,i5,4e15.7)
 2060 format(/5x,'Solution in the domain:',//,5x,'internal points'//)
 2070 format(4x,'x coor',4x,'y coor',3x,'xdisplacement',3x,
     @'ydisplacement'/)  
 2080 format(2f10.3,1x,e15.8,1x,e15.8)
 2090 format(/4x,'x coor',4x,'y coor',7x,'stress xx'7x,'stress yy',7x,
     @'stress xy'/)  
 2100 format(2f10.3,1x,e15.8,1x,e15.8,1x,e15.8)
 3010 format(5x,'WARNING:',/,1x,
     @'The Cauchy relation has been used at node',i4,'.',/,1x,
     @'The results obtained suggest revising the hypothesis of',/,1x,
     @'definition of the stress tensor at this node',i4,' used for',/,       
     @1x,'the application of Cauchy relation.',/)

      end 
      subroutine numer ( xp, yp, x1, y1, x2, y2, a, b, e, g, pois,
     &  dm, igauss, xi, omeg ) 

c*********************************************************************72
c
cc NUMER numerical calculates integration constants of the fundamental equation.
c
c  Discussion:
c
c    The subroutine numer calculates the values of the integration
c    constants of the fundamental equation, integrating numerically
c    from the point of coordinates xp,yp along the element whose 
c    first node has coordinates x1,y1 and whose second node has x2,y2
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
c     list of the additional main arrays used in this subroutine:
c
c   xi:   stores the natural coordinates of the Gauss points used
c         for the numerical integration.
c
c   omeg: stores the weight functions for the Gauss points.
c
c   u:    represents the displacements of the Kelvin fundamental 
c         solution.
c
c   p:    represents the stresses of the Kelvin fundamental solution.
c
c   f:    represents the two linear interpolation functions.
c
c   a,b:  have the same meaning as in the main program.
c
      implicit real*8(a-h,o-z)

      dimension xi(50),omeg(50),a(2,4),b(2,4),f(2),u(2,2),p(2,2)

      imp = 11

      ax = (x2-x1)/2  
      ay = (y2-y1)/2  
      bx = (x2+x1)/2  
      by = (y2+y1)/2
c
c     dis represents the distance from the node under consideration
c     to the element along which the integration is going to be 
c     performed.
c 
      if(ax.eq.0.0)then  
c
c     Particular case of a vertical element.
c
        dis = dabs(xp-x1)
      else
        dis = dabs((ay*xp/ax-yp+y1-ay*x1/ax)/dsqrt((ay/ax)**2+1))
      end if
      sig = (x1-xp)*(y2-yp)-(x2-xp)*(y1-yp) 
      if(sig.lt.0) dis = -dis
c
c     The outward normal to the element is calculated.
c
      xnue =  ay/dsqrt(ax**2+ay**2) 
      ynue = -ax/dsqrt(ax**2+ay**2)
c
c     The arrays to be generated are initialized.
c
      do i = 1,2 
        do j = 1,4 
          a(i,j) = 0. 
          b(i,j) = 0.
        end do
      end do

      pi = 3.14159265*4*(1-pois)
      do 50 i = 1,igauss 
c
c     xc and yc are the coordinates of the integration points.
c 
      xc = ax*xi(i)+bx  
      yc = ay*xi(i)+by  
      r = dsqrt((xp-xc)**2+(yp-yc)**2) 
c
c     The values of the interpolation function at the integration 
c     point are calculated.
c
      f(1) = -(xi(i)-1)/2 
      f(2) = (xi(i)+1)/2  
c 
c     The values of the fundamental solution at the integration points
c     are calculated.
c 
      u(1,1) = ((3-4*pois)*dlog(1/r)+(xc-xp)**2/r**2)/(2*pi*g)*e/dm 
      u(1,2) = (xc-xp)*(yc-yp)/(2*pi*g*r**2)*e/dm 
      u(2,1) = u(1,2) 
      u(2,2) = ((3-4*pois)*dlog(1/r)+(yc-yp)**2/r**2)/(2*pi*g)*e/dm
      p(1,1) = -dis*((1-2*pois)+2*(xc-xp)**2/r**2)/(pi*r**2)  
      p(1,2) = -((2*dis*(xc-xp)*(yc-yp))/r**2+(1-2*pois)  
     @        *(xnue*(yc-yp)-ynue*(xc-xp)))/(pi*r**2) 
      p(2,1) = -((2*dis*(xc-xp)*(yc-yp))/r**2+(1-2*pois)*(ynue*(xc-xp)  
     @       -xnue*(yc-yp)))/(pi*r**2) 
      p(2,2) = -dis*((1-2*pois)+2*(yc-yp)**2/r**2)/(pi*r**2)

      do l = 1,2 
        do k = 1,2 
c
c     The values of the fundamental solution are multiplied by the 
c     interpolation function and by the weight function.
c 
          b(l,k) = b(l,k)+f(k)*u(l,1)*omeg(i)*dsqrt(ax**2+ay**2) 
          b(l,k+2) = b(l,k+2)+f(k)*u(l,2)*omeg(i)*dsqrt(ax**2+ay**2) 
          a(l,k) = a(l,k)+f(k)*p(l,1)*omeg(i)*dsqrt(ax**2+ay**2) 
          a(l,k+2) = a(l,k+2)+f(k)*p(l,2)*omeg(i)*dsqrt(ax**2+ay**2)

        end do
      end do
c
   50 continue
c
c     The array u is used as intermediate array.
c
      do 60 i = 1,2
      do 60 j = 1,2
        u(i,j) = b(i,j)
 60   continue
c
c     The matrix b is expressed in local coordinates.
c
      do 70 i = 1,2
      do 70 j = 1,2
        b(i,j) = b(i,j)*xnue+b(i,j+2)*ynue
        b(i,j+2) = -u(i,j)*ynue+b(i,j+2)*xnue
70    continue

      return  
      end
      subroutine ana ( rlong, coste, sinte, a, b, e, g, pois, dm, 
     &  icontrol ) 

c*********************************************************************72
c
cc ANA analytically calculates integration constants of the fundamental equation.
c
c  Discussion:
c
c    The subroutine ana calculates the values of the integration
c    constants of the fundamental equation, integrating analytically
c    from the node l or l+1 along the element l. 
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
c     List of the additional main arrays used in this subroutine:
c
c     a,b:  have the same meaning as in the main program
c
      implicit real*8(a-h,o-z)

      dimension a(2,4),b(2,4),u(2,2)
c
c     The arrays to be generated are initialized.
c
      do i = 1,2 
        do j = 1,4 
          a(i,j) = 0.0
          b(i,j) = 0.0
        end do
      end do

      pi = 3.14159265*4.*(1.-pois)
      c = (1-2*pois)/pi
      d = e/(2*pi*g)/dm
      if(icontrol.eq.1) then
c
c     Node 1
c
      a(1,3) = -c*(-1.0+dlog(rlong))
      a(1,4) = -c
      a(2,1) = -a(1,3)
      a(2,2) = -a(1,4)
      b(1,1) =  (d*rlong/2.)*((3.-4.*pois)*(1.5-dlog(rlong))+coste**2)
      b(1,2) =  (d*rlong/2.)*((3.-4.*pois)*(0.5-dlog(rlong))+coste**2)
      b(1,3) =  (d*rlong/2.)*coste*sinte
      b(1,4) = b(1,3)
      b(2,1) = b(1,4)
      b(2,2) = b(1,3)
      b(2,3) =  (d*rlong/2.)*((3.-4.*pois)*(1.5-dlog(rlong))+sinte**2)
      b(2,4) =  (d*rlong/2.)*((3.-4.*pois)*(0.5-dlog(rlong))+sinte**2)
      else
c
c     Node 2
c
      a(1,3) = c
      a(1,4) = c*(-1.0+dlog(rlong))
      a(2,1) = -a(1,3)
      a(2,2) = -a(1,4)
      b(1,1) =  (d*rlong/2.)*((3.-4.*pois)*(0.5-dlog(rlong))+coste**2)
      b(1,2) =  (d*rlong/2.)*((3.-4.*pois)*(1.5-dlog(rlong))+coste**2)
      b(1,3) =  (d*rlong/2.)*coste*sinte
      b(1,4) = b(1,3)
      b(2,1) = b(1,3)
      b(2,2) = b(1,4)
      b(2,3) =  (d*rlong/2.)*((3.-4.*pois)*(0.5-dlog(rlong))+sinte**2)
      b(2,4) =  (d*rlong/2.)*((3.-4.*pois)*(1.5-dlog(rlong))+sinte**2)
      end if
c
c     The array u is used as intermediate array.
c
      do i = 1,2
        do j = 1,2
          u(i,j) = b(i,j)
        end do
      end do
c
c     The matrix b is expressed in local coordinates.
c
      do i = 1,2
        do j = 1,2
          b(i,j) = b(i,j)*sinte-b(i,j+2)*coste
          b(i,j+2) = u(i,j)*coste+b(i,j+2)*sinte
        end do
      end do

      return  
      end
      subroutine cala ( xe, ye, x, y, xa, ya, ue, ve, u, v, ua, va, 
     &  a, b, c, d, k )  

c*********************************************************************72
c
cc CALA relates principal stresses to stress vectors in local coordinates.
c
c  Discussion:
c
c    The routine calculates the coefficients that relate
c    the stress vector in local coordinates with the principal 
c    stresses.
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
      implicit real*8(a-h,o-z)

      den = (xa-x)*(ye-y)-(xe-x)*(ya-y) 
c
c     epx, epy, epxy/2 are the components of the strain tensor.  
c 
      epx = ((y-ya)*ue+(ya-ye)*u+(ye-y)*ua)/den 
      epy = -((x-xa)*ve+(xa-xe)*v+(xe-x)*va)/den  
      epxy = ((y-ya)*ve+(ya-ye)*v+(ye-y)*va)/den-((x-xa)*ue+(xa-xe)*u
     @     +(xe-x)*ua)/den
c
      epxc = dabs(epx/1000.)
      epyc = dabs(epy/1000.)
      epxyc = dabs(epxy)
c
c     To detect whether the system xy is already the principal one,
c     the value of epxy is compared with the maximum value of the 
c     normal components of the strain tensor, an absolute comparison 
c     not being able to be performed due to the variability (though 
c     always small) of the values of the strain.
c
      if(epxyc.gt.epxc.or.epxyc.gt.epyc) go to 50
c
c     The strain tensor in the system xy is assumed to be the 
c     principal.
c
c     xnui,ynui and xnuii,ynuii represent the components x,y 
c     of the unit vectors of the principal directions I and II.
c 
      xnui = 1.0
      ynui = 0.0
      xnuii = 0.0
      ynuii = 1.0
      go to 75
c
c     ep1, ep2 are the principal strains
c 
   50 continue
      ep1 = (epx+epy+dsqrt((epx-epy)**2+epxy**2))/2  
      ep2 = (epx+epy-dsqrt((epx-epy)**2+epxy**2))/2
c
c     The unit vectors of the principal strains are calculated.
c
c     Axis I
c
      xnui = (epy-ep1)/dsqrt((epy-ep1)**2+(epxy/2)**2) 
      ynui = -epxy/2/dsqrt((epy-ep1)**2+(epxy/2)**2)
c
c     Axis II 
c
      xnuii = ynui/dsqrt(xnui**2+ynui**2)  
      ynuii = -xnui/dsqrt(xnui**2+ynui**2)
c
c     xnue,ynue are the components of the outward normal to the 
c     element.
c
  75  continue
c
c     epx and epy are now used as intermediate variables.
c
      if ( k .eq. 2 ) then
        epx = xe  
        epy = ye  
      else
        epx = xa  
        epy = ya  
      end if

  100 xnue = (-1)**k*(y-epy)/dsqrt((x-epx)**2+(y-epy)**2) 
      ynue = (-1)**k*(epx-x)/dsqrt((x-epx)**2+(y-epy)**2) 
c
c     a,b represent the coefficients that relate the normal stress
c     with the principal stresses I and II. c,d represent the 
c     coefficients that relate the tangential stress with
c     the principal stresses I and II.
c
      a = xnue*xnui+ynue*ynui 
      b = xnue*xnuii+ynue*ynuii 
      c = -a*b  
      d = -c  
      a = a**2  
      b = b**2

      return  
      end
      subroutine piyay ( xp, yp, x1, y1, x2, y2, hh, gg, g, pois, 
     &  igauss, xi, omeg )  

c*********************************************************************72
c
cc PIYAY calculates certain integration constants along boundary elements.
c
c  Discussion:
c
c    This subroutine calculates the values of the integration 
c    constants along the boundary elements corresponding to the 
c    integral expression of the stresses at internal points.
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
c     list of the additional main arrays used in this subroutine:
c
c   xi:   stores the natural coordinates of the Gauss points used
c         for the numerical integration
c   omeg: stores the weight functions for the Gauss points
c   d:    represents the kernel associated to the Kelvin fundamental 
c         solution that multiplies to the displacements in the 
c         integral expression of the stresses
c   s:    represents the kernel associated to the Kelvin fundamental 
c         solution that multiplies to the stress vector in the 
c         integral expression of the stresses
c   f:    represents the two linear interpolation functions
c   hh:   represents the integration constants matrix which multiplies 
c         to the displacements
c   gg:   represents the integration constants matrix which multiplies 
c         to the strees vector
c
      implicit real*8(a-h,o-z)

      real*8 d(3,2)
      real*8 f(2)
      real*8 gg(3,4)
      real*8 hh(3,4)
      integer i
      integer j
      real*8 omeg(50)
      real*8 s(3,2)
      real*8 xi(50)

      ax = (x2-x1)/2  
      ay = (y2-y1)/2  
      bx = (x2+x1)/2  
      by = (y2+y1)/2  
c
c     dis represents the distance from the point under consideration
c     to the element along which the integration is going to be 
c     performed.
c 
      if(ax.eq.0.0)then
c
c     Particular case of a vertical element.
c
        dis = dabs(xp-x1)
      else
        dis = dabs((ay*xp/ax-yp+y1-ay*x1/ax)/dsqrt((ay/ax)**2+1)) 
      end if

      sig = (x1-xp)*(y2-yp)-(x2-xp)*(y1-yp) 
      if(sig.lt.0) dis = -dis
c
c     The outward normal to the element is calculated.
c
      xnue = ay/dsqrt(ax**2+ay**2) 
      ynue = -ax/dsqrt(ax**2+ay**2)
c
c     The arrays to be generated are initialized.
c
      do 40 i = 1,3 
      do 40 j = 1,4 
        hh(i,j) = 0.  
        gg(i,j) = 0.
   40 continue
c
      pi1 = 3.14159265*4*(1-pois) 
      pi2 = pi1/g/2
c
      do 50 i = 1,igauss
c
c     xc and yc are the coordinates of the integration points.
c 
      xc = ax*xi(i)+bx  
      yc = ay*xi(i)+by  
      r = dsqrt((xp-xc)**2+(yp-yc)**2) 
      c1 = 1-2*pois 
      d1 = (xc-xp)/r  
      d2 = (yc-yp)/r  
c
c     The values of the interpolation function at the integration 
c     point are calculated.
c
      f(1) = -(xi(i)-1)/2 
      f(2) = (xi(i)+1)/2  
c 
c     The values of the fundamental solution at the integration points
c     are calculated.
c 
      d(1,1) = (c1*d1+2*d1**3)/(pi1*r)  
      d(3,2) = (c1*d2+2*d2**3)/(pi1*r)  
      d(1,2) = (-c1*d2+2*d1**2*d2)/(pi1*r)  
      d(3,1) = (-c1*d1+2*d1*d2**2)/(pi1*r)  
      d(2,1) = (c1*d2+2*d1**2*d2)/(pi1*r) 
      d(2,2) = (c1*d1+2*d1*d2**2)/(pi1*r) 
      s(1,1) = (2*dis/r*(d1-4*d1**3)+2*xnue*d1**2+xnue)/(pi2*r**2)  
      s(3,2) = (2*dis/r*(d2-4*d2**3)+2*ynue*d2**2+ynue)/(pi2*r**2)  
      s(1,2) = (2*dis/r*(c1*d2-4*d1**2*d2)+4*pois*xnue*d1*d2+c1*2 
     @       *ynue*d1**2-(1-4*pois)*ynue)/(pi2*r**2)  
      s(3,1) = (2*dis/r*(c1*d1-4*d1*d2**2)+4*pois*ynue*d1*d2+c1*2* 
     @       xnue*d2**2-(1-4*pois)*xnue)/(pi2*r**2)  
      s(2,1) = (2*dis/r*(pois*d2-4*d1**2*d2)+2*pois*ynue*d1**2+2*
     @        (1-pois)*xnue*d1*d2+c1*ynue)/(pi2*r**2)  
      s(2,2) = (2*dis/r*(pois*d1-4*d1*d2**2)+2*pois*xnue*d2**2+2*
     @       (1-pois)*ynue*d1*d2+c1*xnue)/(pi2*r**2)  
c
c     The values of the fundamental solution are multiplied by the 
c     interpolation function and by the weight function.
c 
      do 55 l = 1,3 
      do 55 k = 1,2 
        hh(l,k) = hh(l,k)+d(l,1)*f(k)*omeg(i)*dsqrt(ax**2+ay**2) 
        hh(l,k+2) = hh(l,k+2)+d(l,2)*f(k)*omeg(i)*dsqrt(ax**2+ay**2) 
        gg(l,k) = gg(l,k)+s(l,1)*f(k)*omeg(i)*dsqrt(ax**2+ay**2) 
        gg(l,k+2) = gg(l,k+2)+s(l,2)*f(k)*omeg(i)*dsqrt(ax**2+ay**2)
   55 continue

   50 continue
c
c     The array d is now used as intermediate array.
c
      do i = 1,3 
        do j = 1,2 
          d(i,j) = hh(i,j)
        end do
      end do
c
c     The matrix hh is expressed in local coordinates.
c 
      do 70 i = 1,3 
      do 70 j = 1,2 
        hh(i,j) = hh(i,j)*xnue+hh(i,j+2)*ynue 
        hh(i,j+2) = -d(i,j)*ynue+hh(i,j+2)*xnue
   70 continue  

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
c    ISBN: 0-19-856543-7,
c    LC: TA347.B69.P34.
c
      implicit none

      integer npmax,ipivo(npmax),n,imp,na,i,j,k,l,pmx
      real*8 a(npmax,npmax),c(npmax),mx,aux

      do 1 j = 1,n

        mx = a(j,j)
        pmx = j

        do i = j+1,n+na
          if(DABS(a(i,j)).gt.DABS(mx)) then
            mx = a(i,j)
            pmx = i
          end if
        end do

            if(DABS(mx).lt.1.e-6) then
                  write(imp,100)
 100        format(5x,' The matrix is Singular')
                  stop 1111
            else
                  if(pmx.ne.j) then
                        k = ipivo(pmx)
                        ipivo(pmx) = ipivo(j)
                        ipivo(j) = k
                        do 3 l = j,n
                              aux = a(pmx,l)
                              a(pmx,l) = a(j,l)
                              a(j,l) = aux
 3                      continue
                        aux = c(pmx)
                        c(pmx) = c(j)
                        c(j) = aux
                  end if
                  do 4 l = j+1,n
                        a(j,l) = a(j,l)/a(j,j)
 4                continue
                  c(j) = c(j)/a(j,j)
                  a(j,j) = 1.0
                  do 5 i = 1,n+na
                        if(i.ne.j) then
                              do 6 l = j+1,n
                         a(i,l) = a(i,l) - a(i,j) * a(j,l)
 6                            continue
                              c(i) = c(i)-a(i,j)*c(j)
                              a(i,j) = 0.0
                        end if
 5                continue
            end if
 1    continue

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
c    ISBN: 0-19-856543-7,
c    LC: TA347.B69.P34.
c
      implicit none

      integer n,k,i,ii,ir1,ir0,expon
      real*8 coef(50),root(0:50,2),tol,epsilon(50),omeg(50)
      real*8 ep(50),om(50)
      tol = 1./10.**14
      if(n.gt.50.or.n.lt.2) then
      write(11,3000)
      stop 1111
      end if
c
c     Initiations.
c
      ir1 = 1
      ir0 = 2
      root(0,1) = -1
      root(1,1) = 0.0
      root(2,1) = 1
c
c     Computation of polynomial coefficients and roots.
c
      do 1 k = 2,n

        do i = 0,n/2
          call coefic(k,i,coef(i+1))
        end do

        do 3 i = 1,k
          call roots(k,coef,root(i-1,ir1),root(i,ir1),
     @               tol,root(i,ir0))
3       continue
        root(k+1,ir0) = 1
        root(0,ir0) = -1
        ii = ir0
        ir0 = ir1
        ir1 = ii
1     continue
c
c     Computation of weights.
c
      do 4 i = 1,n
        epsilon(i) = root(i,ir1)
        omeg(i) = 0
        do 5 k = 0,n/2
          expon = n-2*k
          if(expon.gt.1)then 
            omeg(i) = omeg(i)+coef(k+1)*expon*epsilon(i)**(expon-1)
          else
            omeg(i) = omeg(i)+coef(k+1)*expon
          end if
5       continue
        omeg(i) = omeg(i)/2**n 
        omeg(i) = 2./(omeg(i)**2*(1-epsilon(i)**2))
        ep(i) = epsilon(i)
        om(i) = omeg(i) 
4     continue

      return
c
c     Printing formats.
c
 3000 format(5x,'ERROR: ',1x,
     @      'The number of Gauss points used is not allowed')
      end
      subroutine coefic ( n, k, value )

c*********************************************************************72
c
cc COEFIC sets the coefficients of a polynomial needed for the Gauss rule.
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
      implicit none

      integer inter,k,n,n1,n2,n3,i
      real*8 sig,value

      value = 1
      inter = k/2
      sig = 1.0
      if(inter*2-k.ne.0) sig = -1
      n2 = n-k
      n1 = 2*n2
      n3 = n-2*k
      if(n1.eq.0)then
        value = 1
      else
        do i = n2+1,n1
          value = value*i
        end do
      end if
      if(n3.ne.0) then
        do 2 i = 1,n3
          value = value/i
2       continue
      end if
      if(k.ne.0) then
        do 3 i = 1,k
          value = value/i
3       continue
      end if
      value = value*sig

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
c    ISBN: 0-19-856543-7,
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
c    ISBN: 0-19-856543-7,
c    LC: TA347.B69.P34.
c
      implicit none

      integer ncoef
      real*8 x00,xff,x0,xf,xi,error,coef(ncoef),value0,valuef,valuei,
     @       errac

      x0 = x00
      xf = xff
      call evalua(ncoef,coef,x0,value0)       
      call evalua(ncoef,coef,xf,valuef) 
      errac = dabs(xf-x0)

      do while (errac.gt.error)

        xi = (xf+x0)/2.
        call evalua(ncoef,coef,xi,valuei)    

        if(dabs(valuei).lt.error) then
          return
        end if

        if(valuei*value0.lt.0.0) then
          xf = xi
          valuef = valuei
        else
          x0 = xi
          value0 = valuei
        end if

        errac = abs(xf-x0)

      end do

      return
      end













