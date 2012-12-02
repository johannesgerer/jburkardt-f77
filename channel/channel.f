      program main

c*********************************************************************72
c
cc MAIN is the main program for CHANNEL.
c
c  Discussion:
c
c    CHANNEL is the main program for the channel flow problem.
c
c    Control for NAVIER-STOKES on channel using primitive variables.
c
c    This program solves a fluid flow problem on the unit square.
c
c    The fluid flow problem is formulated in terms of primitive variables 
c    - u, v, and p.
c
c    This code tries to match a downstream profile by altering
c    one parameter, the value of the inflow parameter at the mid-height
c    of the channel.
c
c    Piecewise linear functions on triangles approximate
c    the pressure and quadratics on triangles approximate the velocity.
c    This is the "Taylor-Hood" finite element basis.
c
c    The primitive variable formulation of the Navier Stokes equations
c    involves horizontal velocity U, vertical velocity V, and pressure P.  
c    The equations are:
c
c      U dUdx + V dUdy + dPdx - mu*(ddU/dxdx + ddU/dydy) = F1
c      U dVdx + V dVdy + dPdy - mu*(ddV/dxdx + ddV/dydy) = F2
c      dUdx + dVdy = 0
c
c    When reformulated into finite element form, with PHI(i) being the I-th
c    common basis function for U and V, and PSI(i) the I-th basis function 
c    for P, these equations become:
c
c      Integral (U dUdx PHI(i) + V dUdy PHI(I) - P dPHI(i)/dx 
c        + mu (dUdx dPHI(i)/dx + dUdy dPHI(i)/dy) ) = 
c        Integral (F1 * PHI(i))
c
c      Integral (U dVdx PHI(i) + V dVdy PHI(I) - P dPHI(i)/dy 
c        + mu (dVdx dPHI(i)/dx + dVdy dPHI(i)/dy) ) = 
c        Integral (F2 * PHI(i))
c
c      Integral (dUdx PSI(i) + dVdx PSI(i)) = 0
c
c  Change Record:
c
c    25 January 2007: Since I'm keeping this code around for reference,
c    I cleaned it up slightly, converting all routines to
c    using IMPLICIT NONE.
c
c    28 February 2006: There was a discrepancy between this code
c    and the FORTRAN90 version, seeming to come from the LINSYS routine. 
c    I worked on that for a while, and the FORTRAN77 code seems to be OK now.
c
c    4 September 1992: Added LONG to GDUMP output.
c
c    2 September 1992:  Set g to zero before secant iteration.
c    Made sure to update g BEFORE returning from NSTOKE.
c    New code converges in 3 secant steps, 13 Newton steps,
c    and 48 seconds.
c
c    27 August 1992: Corrected a mistake in XYDUMP that only occurred
c    for NX <= NY.
c
c    26 August 1992: broke out SETLIN, SETQUD and SETBAS for 
c    similarity with BUMP.
c
c    22 August 1992: Moved bandwidth calculations to SETBAN, mainly for
c    the benefit of the BUMP computation.
c
c    20 August 1992: Moved linear system setup and solve out of NSTOKE,
c    merged it with solvlin, and called it "LINSYS".
c
c    19 August 1992.  After I changed the double precision of the matrix
c    A in NSTOKE from (MAXROW,NEQN) to (NROW,NEQN) I found performance
c    was degraded.  The time went up from about 47 seconds to 57 seconds!
c    (I made other, also theoretically harmless changes at the same time).
c    Changing back to MAXROW brought the time down to 50 seconds.  What's
c    going on?  This is definitely a performance problem in LAPACK!
c
c    18 August 1992, the residual calculation has been healed, somehow!
c    That means that I can, if I wish, think about a jacobian.
c
c    13 August 1992, the way that the last pressure was set to 1 was not
c    working properly.  I don't know why, but I changed the code so that
c    the last pressure was set to zero, and the contour plots smoothed out,
c    and the secant iteration converged in three steps, rather than four.
c    So I'm also dropping the entire "PRESET" routine, which forced the
c    average pressure to zero.  What's the point?  Perhaps that's important
c    when you're modifying the grid, though.  So I won't actually delete 
c    the code.
c
c    11 August 1992, the original scheme for the last pressure equation
c    does NOT force the average pressure to be zero.  If you don't
c    believe me, repeat the calculation of PMEAN after the adjustment.
c    I've fixed that.
c
c  List of routines:
c
c  BSP     is given a point (X,Y), and evaluates, for a particular
c          triangle, the linear basis functions associated with
c          pressure.
c
c  CHANNEL is the main program.  It initializes data, and controls
c          an iteration which is intended to produce a value of the
c          parameter A which produces a desired internal flow profile.
c
c  DELETE  delete old copies of files before a new copy is opened.
c
c  GDUMP   writes information about a particular time step to a file for
c          use by a graphics display program.
c
c  GETG    returns a vector containing the values of a quantity at each
c          node where the internal profile is taken, by extracting those
c          values from a global list.
c
c  GRAM    computes the Gram matrix, GR(I,J)=INTEGRAL PHI(I)*PHI(J),
c          and the line integral of the internal velocity profile,
c          INTEGRAL UI*PHI(J).
c
c  IGETL   returns the local unknown number, given the global unknown 
c          number for any node on the line X=XZERO, where the internal
c          velocity profile is taken.
c
c  LINSYS  sets up and solves the Navier Stokes linear system.
c
c  NSTOKE  solves the Navier Stokes equations for a given set of
c          boundary conditions, including the inflow determined by
c          the parameter A.  The solution is done via Newton's method.
c
c  PVAL    computes a table of pressures for UVDUMP.
c
c  QBF     is given a point (X,Y), and evaluates
c          a particular basis function for a particular quadratic
c          triangular element at that point.
c
c  SETBAN  computes the bandwidth of the matrix.
c
c  SETBAS  evaluates and saves the value of the basis functions
c          at the quadrature points.
c
c  SETGRD  computes the location and numbering of nodes.
c          Constructs and numbers the elements, and assigns them nodes.
c          Calculates information needed for numerical quadrature.
c          Calculates the number of unknowns.
c
c  SETLIN  records the indices of the horizontal velocities that
c          lie along the profile sampling line.
c
c  SETQUD  computes the location of the quadrature points.
c
c  UBDRY   produces the value of the inflow at a given point as
c          determined by the value of the parameter A.
c
c  UVAL    evaluate the horizontal and vertical velocities, and their
c          derivatives with respect to X and Y, at a given point in an
c          element.
c
c  UVDUMP  dump velocites and pressures for PLOT3D.
c
c  XYDUMP  dump X and Y locations for PLOT3D.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    02 January 2009
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    DOUBLE PRECISION A(MAXROW,MAXEQN), the banded matrix 
c    used in NSTOKE and SOLVLIN.  It is stored according to
c    LINPACK/LAPACK general band storage mode.
c
c    DOUBLE PRECISION AREA(NELEMN), contains the area of each
c    element.
c
c    DOUBLE PRECISION F(MAXEQN).  After the call to NSTOKE,
c    F contains the solution to the current Navier Stokes problem.
c
c    CHARACTER*30 FILEG, the filename of the graphics data file
c    to be dumped for the DISPLAY program.
c
c    CHARACTER*30 FILEU, the filename of the graphics UV data file
c    to be dumped for the PLOT3D program.
c
c    CHARACTER*30 FILEX, the filename of the graphics XY data file
c    to be dumped for the PLOT3D program.
c
c    DOUBLE PRECISION G(MAXEQN).  After the call to SOLVLIN,
c    G contains the sensitivities.
c
c    INTEGER INDX(NP,2), records, for each node, whether
c    any velocity unknowns are associated with the node.
c    INDX(I,1) records information for horizontal velocities.
c    If it is 0, then no unknown is associated with the
c    node, and a zero horizontal velocity is assumed.
c    If it is -1, then no unknown is associated with the
c    node, but a velocity is specified via the UBDRY routine
c    If it is positive, then the value is the index in the
c    F and G arrays of the unknown coefficient.
c    INDX(I,2) records information for vertical velocities.
c    If it is 0, then no unknown is associated with the
c    node, and a zero vertical velocity is assumed.
c    If it is positive, then the value is the index in the
c    F and G arrays of the unknown coefficient.
c
c    INTEGER INSC(NP), records, for each node, whether an
c    unknown pressure is associated with the node.
c    If INSC(I) is zero, then no unknown is associated with the
c    node, and a pressure of 0 is assumed.  
c    If INSC(I) is positive, then the value is the index in the
c    F and G arrays of the unknown coefficient.
c
c    INTEGER IOUNIT, the FORTRAN unit number for the file
c    into which graphics data will be written.
c
c    INTEGER IPIVOT(MAXEQN), pivot information used by the linear
c    solver.
c
c    INTEGER ISOTRI(NELEMN), records whether or not a given 
c    element is isometric or not.
c    0, element is not isometric.
c    1, element is isometric.
c
c    Input, INTEGER IVUNIT, the output unit to which the data is to
c    be written.  The data file is unformated.
c
c    INTEGER IWRITE, controls the amount of output produced by the
c    program.
c    0, minimal output.
c    1, normal output, plus graphics data files created by
c       GDUMP, UVDUMP and XYDUMP.
c    2, copious output.
c
c    Input, INTEGER IXUNIT, the output unit to which the data is to
c    be written.  The data file is unformated.
c
c    LOGICAL LONG, 
c    .TRUE. if region is "long and thin", and
c    .FALSE. if region is "tall and skinny".
c    This determines how we number nodes, elements, and variables.
c
c    INTEGER MAXEQN, the maximum number of equations allowed.
c
c    INTEGER MAXNEW, the maximum number of Newton steps per iteration.
c
c    INTEGER MAXROW, the maximum row double precision of the coefficient
c    matrix that is allowed.
c
c    INTEGER MAXSEC, the maximum number of secant steps allowed.
c
c    INTEGER MX, MX = 2*NX-1, the total number of grid points on
c    the horizontal side of the region.
c
c    INTEGER MY, MY = 2*NY-1, the total number of grid points on
c    the vertical side of the region.
c
c    INTEGER NELEMN, the number of elements used.   
c
c    INTEGER NEQN, the number of equations or functions for
c    the full system.
c
c    INTEGER NLBAND, the number of diagonals below the main diagonal 
c    of the matrix A which are nonzero. 
c
c    INTEGER NNODES, the number of nodes per element, 6.
c
c    INTEGER NODE(NELEMN,6), records the global node numbers of the
c    6 nodes that make up each element.
c
c    INTEGER NODEX0, the lowest numbered node in the column of
c    nodes where the profile is measured.
c
c    INTEGER NP, the number of nodes.
c
c    INTEGER NQUAD, the number of quadrature points, currently
c    set to 3.
c
c    INTEGER NROW, the used row double precision of the coefficient
c    matrix.
c
c    INTEGER NUMNEW, total number of Newton iterations taken in
C    the NSTOKE routine during the entire run.
c
c    INTEGER NUMSEC, the number of secant steps taken.
c
c    INTEGER NX, the number of "main" grid points on the horizontal 
c    side of the region.
c
c    INTEGER NY, the number of "main" grid points on the vertical 
c    side of the region.
c
c    REAL PHI(NELEMN,NQUAD,NNODES,3).  Each entry of PHI contains
c    the value of a quadratic basis function or its derivative, 
c    evaluated at a quadrature point.
c    In particular, PHI(I,J,K,1) is the value of the quadratic basis
c    function associated with local node K in element I, evaluatated
c    at quadrature point J.
c    PHI(I,J,K,2) is the X derivative of that same basis function,
c    PHI(I,J,K,3) is the Y derivative of that same basis function.
c
c    REAL PSI(NELEMN,NQUAD,NNODES).  Each entry of PSI contains
c    the value of a linear basis function evaluated at a 
c    quadrature point.
c    PSI(I,J,K) is the value of the linear basis function associated
c    with local node K in element I, evaluated at quadrature point J.
c
c    REAL RES(MAXEQN), contains the residuals.
c
c    REAL REYNLD, the value of the Reynolds number.  In the
c    program's system of units, viscosity = 1 / REYNLD.
c
c    REAL RJPNEW, the derivative with respect to the
c    parameter A of the functional J.
c
c    REAL TARRAY(2), an array needed in order to store 
c    results of a call to the UNIX CPU timing routine ETIME.
c
c    REAL TOLNEW, the convergence tolerance for the Newton
c    iteration in NSTOKE.
c
c    REAL TOLSEC, the convergence tolerance for the secant
c    iteration in the main program.
c
c    REAL XC(NP), XC(I) is the X coordinate of node I.
c
c    REAL XLNGTH, the length of the region.
c
c    REAL XM(NELEMN,NQUAD), XM(IT,I) is the X coordinate of the 
c    I-th quadrature point in element IT.
c
c    REAL YC(NP), YC(I) is the Y coordinate of node I.
c
c    REAL YLNGTH, the height of the region.
c
c    REAL YM(NELEMN,NQUAD).  YM(IT,I) is the Y coordinate of the
c    I-th quadrature point in element IT.
c
      implicit none

      integer nx
      parameter (nx=21)

      integer ny
      parameter (ny=7)
c
c  This assignment should really read (maxrow=27*min(nx,ny))
c
      integer maxrow
      parameter (maxrow=27*ny)

      integer nelemn
      parameter (nelemn=2*(nx-1)*(ny-1))

      integer mx
      parameter (mx=2*nx-1)

      integer my
      parameter (my=2*ny-1)

      integer np
      parameter (np=mx*my)

      integer maxeqn
      parameter (maxeqn=2*mx*my+nx*ny)

      integer nnodes
      parameter (nnodes = 6)

      integer nquad
      parameter (nquad=3)

      double precision a(maxrow,maxeqn)
      double precision a2
      double precision abound
      double precision anew
      double precision aold
      double precision area(nelemn)
      double precision dcda(my)
      double precision f(maxeqn)
      character fileg*30
      character fileu*30
      character filex*30
      double precision g(maxeqn)
      double precision gr(my,my)
      integer i
      integer iline(my)
      integer indx(np,2)
      integer insc(np)
      integer iounit
      integer ipivot(maxeqn)
      integer isotri(nelemn)
      integer iter
      integer ivunit
      integer iwrite
      integer ixunit
      integer j
      logical long
      integer maxnew
      integer maxsec
      integer nband
      integer neqn
      integer nlband
      integer node(nelemn,nnodes)
      integer nodex0
      integer npara
      integer nrow
      integer numnew
      integer numsec
      double precision para
      double precision para2
      double precision phi(nelemn,nquad,nnodes,3)
      double precision psi(nelemn,nquad,nnodes)
      double precision r(my)
      double precision res(maxeqn)
      double precision reynld
      double precision rjpnew
      double precision rjpold
      double precision rtemp(my)
      real tarray(2)
      double precision temp
      double precision test
      double precision tolnew
      double precision tolsec
      double precision ui(my)
      double precision unew(my)
      real value
      double precision xc(np)
      double precision xlngth
      double precision xm(nelemn,nquad)
      double precision yc(np)
      double precision ylngth
      double precision ym(nelemn,nquad)

      call timestamp ( )
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'CHANNEL'
      write ( *, '(a)' ) '  FORTRAN77 version'
      write ( *, '(a)' ) '  Channel flow control problem'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Last modified:'
      write ( *, '(a)' ) '    4 September 1992.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Flow control problem:'
      write ( *, '(a)' ) '    Inflow controlled by one parameter.'
      write ( *, '(a)' ) '    Velocities measured along vertical line.'
      write ( *, '(a)' ) '    Try to match specified velocity profile.'
c
c  Set input data
c
      fileg = 'display.txt'
      fileu = 'uv.txt'
      filex = 'xy.txt'
      iounit = 2
      ivunit = 4
      iwrite = 2
      ixunit = 3
      maxnew = 10
      maxsec = 8
      npara = 1
      numnew = 0
      numsec = 0
      reynld = 1.0D+00
      rjpnew = 0.0D+00
      tolnew = 1.0D-04
      tolsec = 1.0D-06
      xlngth = 10.0D+00
      ylngth = 3.0D+00

      write ( *, * ) ' '
      write ( *, * ) '  NX = ', nx
      write ( *, * ) '  NY = ', ny
      write ( *, * ) '  Number of elements = ', nelemn
      write ( *, * ) '  Reynolds number = ', reynld
      write ( *, * ) '  Secant tolerance = ', tolsec
      write ( *, * ) '  Newton tolerance = ', tolnew
      write ( *, * ) ' '
c
c  SETGRD constructs grid, numbers unknowns, calculates areas,
c  and points for midpoint quadrature rule.
c
      call setgrd (indx,insc,isotri,iwrite,long,maxeqn,mx,my, 
     &  nelemn,neqn,nnodes,node,np,nx,ny)
c
c  Compute the bandwidth
c
      call setban(indx,insc,maxrow,nband,nelemn,nlband,nnodes, 
     &  node,np,nrow)
c
c  Record variable numbers along profile sampling line.
c
      call setlin(iline,indx,iwrite,long,mx,my,nodex0,np, 
     &  nx,ny,xlngth)
c
c  Set the coordinates of grid points.
c
      call setxy(iwrite,long,mx,my,np,nx,ny,xc,xlngth,yc,ylngth)
c
c  Set quadrature points
c
      call setqud(area,nelemn,nnodes,node,np,nquad,xc,xm,yc,ym)
c
c  Evaluate basis functions at quadrature points
c
      call setbas(nelemn,nnodes,node,np,nquad,phi,psi,xc,xm,yc,ym)
c
c  NSTOKE now solves the Navier Stokes problem for an inflow 
c  parameter of 1.0.

      para = 1.0D+00
      write (*,*) ' '
      write (*,*) 'Solve Navier Stokes problem with parameter = ',para
      write (*,*) 'for profile at x=',xc(nodex0)

      do i = 1, neqn
        g(i) = 1.0D+00
      end do

      call nstoke (a,area,f,g,indx,insc,ipivot,iwrite, 
     &  maxnew,maxrow,nelemn,neqn,nlband,nnodes,node, 
     &  np,nquad,nrow,numnew,para,phi,psi,reynld,tolnew,yc)
c
c  RESID computes the residual at the given solution
c
      if(iwrite.ge.1)then
        call resid (area,f,indx,insc,iwrite,nelemn,neqn, 
     &    nnodes,node,np,nquad,para,phi,psi,res,reynld,yc)
      end if
c
c  GETG computes the internal velocity profile at X=XC(NODEX0), which will
c  be used to measure the goodness-of-fit of the later solutions.
c
      call getg (f,iline,my,neqn,ui)

      if(iwrite.ge.1)then
        write (*,*) ' '
        write (*,*) 'U profile:'
        write (*,*) ' '
        write (*,'(5g14.6)') (ui(i),i=1,my)
      end if
c
c  GRAM generates the Gram matrix GR and the vector 
c  R=line integral of ui*phi
c
      call gram (gr,iline,indx,iwrite,my,nelemn,nnodes,node, 
     &   nodex0,np,para,r,ui,xc,yc)
c
c  GDUMP dumps information for graphics display by DISPLAY.
c
      if ( .false. )then
        write (*,*) 'Writing graphics data to file '//fileg
        call delete(fileg)
        open (unit=iounit,file=fileg,form='formatted',status='new', 
     &      err=50)
        rjpnew=0.0D+00
        call gdump (f,indx,insc,iounit,isotri,long,nelemn,neqn, 
     &      nnodes,node,np,npara,nx,ny,para,reynld,rjpnew,xc,yc)
      end if
c
c  XY_PLOT3D creates a PLOT3D "XY" file
c
      if ( .false. ) then
        call delete(filex)
        open(unit=ixunit,file=filex,form='formatted',status='new')
        call xy_plot3d (ixunit,long,np,nx,ny,xc,yc)
        close(unit=ixunit)
c
c  XY_TABLE creates a TABLE "XY" file
c
      else
        call delete ( filex )
        open ( unit = ixunit, file = filex, form = 'formatted',
     &    status = 'new' )
        call xy_table ( ixunit, np, xc, yc )
        close ( unit = ixunit )
      end if
c
c  UV_PLOT3D creates a PLOT3D "Q" file
c
      if ( .false. ) then
        call delete(fileu)
        open(unit=ivunit,file=fileu,form='formatted',status='new')
        call uv_plot3d (f,indx,insc,ivunit,long,mx,my, 
     &     nelemn,neqn,nnodes,node,np,para,a,reynld,yc)
        close(unit=ivunit)
c
c  UV_PLOT3D creates a PLOT3D "Q" file
c
      else
        call delete ( fileu )
        open ( unit = ivunit, file = fileu, form = 'formatted',
     &    status = 'new' )
        call uv_table ( f, indx, ivunit, neqn, np, para, yc )
        close ( unit = ivunit )
      end if
c
c  Destroy information about true solution.
c
      do i = 1, neqn
        f(i) = 0.0D+00
      end do

      do i = 1, neqn
        g(i) = 0.0D+00
      end do
c
c  Secant iteration loop
c
      aold = 0.0D+00
      rjpold = 0.0D+00
      anew = 0.1D+00

      do 30 iter = 1, maxsec

        numsec=numsec+1
        write (*,*) ' '
        write (*,*) 'Secant iteration ',iter
c
c  Solve for unew at new value of parameter anew
c
        write (*,*) ' '
        write (*,*) 'Solving Navier Stokes problem for parameter=',anew
c
c  Use solution F at previous value of parameter for starting point.
c
        do i = 1, neqn
          g(i) = f(i)
        end do

        para=anew

        call nstoke (a,area,f,g,indx,insc,ipivot,iwrite, 
     &      maxnew,maxrow,nelemn,neqn,nlband,nnodes,node, 
     &      np,nquad,nrow,numnew,para,phi,psi,reynld,tolnew,yc)
c
c  Get velocity profile
c
        call getg (f,iline,my,neqn,unew)

        if(iwrite.ge.1)then
          write (*,*) ' '
          write (*,*) 'Velocity profile:'
          write (*,*) ' '
          write (*,'(5g14.6)') (unew(i),i=1,my)
        end if
c
c  Solve linear system for du/da
c
        para=anew
        abound = 1.0D+00

        call linsys (a,area,g,f,indx,insc,ipivot, 
     &     maxrow,nelemn,neqn,nlband,nnodes,node, 
     &    np,nquad,nrow,para,abound,phi,psi,reynld,yc)
c
c  Output in DCDA
c
        call getg (g,iline,my,neqn,dcda)

        if(iwrite.ge.2)then
          write (*,*) ' '
          write (*,*) 'Sensitivities:'
          write (*,*) ' '
          write (*,'(5g14.6)') (dcda(i),i=1,my)
        end if
c
c  Evaluate J prime at current value of parameter where J is 
c  functional to be minimized.
c
c  JPRIME = 2.0 * DCDA(I) * (GR(I,J)*UNEW(J)-R(I))
c
        rjpnew = 0.0D+00
        do i = 1, my
          temp = -r(i)
          do j = 1, my
            temp = temp + gr(i,j) * unew(j)
          end do
          rjpnew = rjpnew + 2.0D+00 * dcda(i) * temp
        end do

        write (*,*) ' '
        write (*,*) 'Parameter =',anew,' J prime = ',rjpnew
c
c  Dump information for graphics
c
        if( .false. )then
          para = anew
          call gdump (f,indx,insc,iounit,isotri,long,nelemn,neqn, 
     &       nnodes,node,np,npara,nx,ny,para,reynld,rjpnew,xc,yc)
        end if
c
c  Update the estimate of the parameter using the secant step
c
        if (iter.eq.1) then
          a2 = 0.5D+00
        else
          a2 = aold-rjpold*(anew-aold)/(rjpnew-rjpold)
        end if

        aold = anew
        anew = a2
        rjpold = rjpnew
        test = abs ( anew - aold ) / abs ( anew )

        write (*,*) 'New value of parameter = ',anew
        write(*,*)'Convergence test=',test

        if (abs(anew-aold).lt.abs(anew)*tolsec) then
          write (*,*) 'Secant iteration converged.'
          go to 40
        end if

   30 continue

      write (*,*) 'Secant iteration failed to converge.'

   40 continue
c
c  Produce total CPU time used.
c
c     tarray(2)=0.0
c     tarray(1)=second()
c
      value = etime (tarray)
      tarray(1) = tarray(1)+tarray(2)
      tarray(2) = tarray(1) / 60.0D+00

      write (*,*) ' '
      write (*,*) 'Total execution time=',tarray(1),' seconds.'
      write (*,*) '                    =',tarray(2),' minutes.'
      write (*,*)'Number of secant steps = ',NUMSEC
      write (*,*)'Number of Newton steps = ',NUMNEW
c
c  Close graphics file
c
      if(.false.)then
        close (unit=iounit)
      end if

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'CHANNEL'
      write ( *, '(a)' ) '  Normal end of execution'

      write ( *, '(a)' ) ' '
      call timestamp ( )

      stop
c
c  Error opening graphics file
c
   50 continue
      write (*,*) 'CHANNEL could not open the graphics file!'
      stop
      end
      function bsp ( xq, yq, it, iq, id, nelemn, nnodes, node, np,
     & xc, yc )

c*********************************************************************72
c
cc BSP evaluates the linear basis functions associated with pressure.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    23 January 2007
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer nelemn
      integer nnodes
      integer np

      double precision bsp
      double precision d
      integer id
      integer i1
      integer i2
      integer i3
      integer iq
      integer iq1
      integer iq2
      integer iq3
      integer it
      integer node(nelemn,nnodes)
      double precision xc(np)
      double precision xq
      double precision yc(np)
      double precision yq

      iq1 = iq
      iq2 = mod(iq,3)+1
      iq3 = mod(iq+1,3)+1

      i1 = node(it,iq1)
      i2 = node(it,iq2)
      i3 = node(it,iq3)
      d = (xc(i2)-xc(i1))*(yc(i3)-yc(i1))
     &  -(xc(i3)-xc(i1))*(yc(i2)-yc(i1))

      if (id.eq.1) then
 
        bsp = 1.0D+00 + ((yc(i2)-yc(i3))
     &  *(xq-xc(i1))+(xc(i3)-xc(i2))*(yq-yc(i1))) / d
 
      elseif (id.eq.2) then
 
        bsp = (yc(i2)-yc(i3)) / d
 
      elseif (id.eq.3) then
 
        bsp = (xc(i3)-xc(i2)) / d
 
      else
 
        write (*,*) 'BSP - fatal error!'
        write (*,*) 'unknown value of id=',id
        stop

      end if

      return
      end
      subroutine daxpy ( n, da, dx, incx, dy, incy )

c*********************************************************************72
c
cc DAXPY: constant times a vector plus a vector.
c
c     uses unrolled loops for increments equal to one.
c     jack dongarra, linpack, 3/11/78.
c     modified 12/3/93, array(1) declarations changed to array(*)
c
      implicit none

      double precision dx(*),dy(*),da
      integer i,incx,incy,ix,iy,m,mp1,n

      if(n.le.0)return
      if (da .eq. 0.0d0) return
      if(incx.eq.1.and.incy.eq.1)go to 20
c
c        code for unequal increments or equal increments
c          not equal to 1
c
      ix = 1
      iy = 1
      if(incx.lt.0)ix = (-n+1)*incx + 1
      if(incy.lt.0)iy = (-n+1)*incy + 1
      do 10 i = 1,n
        dy(iy) = dy(iy) + da*dx(ix)
        ix = ix + incx
        iy = iy + incy
   10 continue
      return
c
c        code for both increments equal to 1
c
c
c        clean-up loop
c
   20 m = mod(n,4)
      if( m .eq. 0 ) go to 40
      do 30 i = 1,m
        dy(i) = dy(i) + da*dx(i)
   30 continue
      if( n .lt. 4 ) return
   40 mp1 = m + 1
      do i = mp1,n,4
        dy(i) = dy(i) + da*dx(i)
        dy(i + 1) = dy(i + 1) + da*dx(i + 1)
        dy(i + 2) = dy(i + 2) + da*dx(i + 2)
        dy(i + 3) = dy(i + 3) + da*dx(i + 3)
      end do

      return
      end
      double precision function ddot ( n, dx, incx, dy, incy )

c*********************************************************************72
c
cc DDOT forms the dot product of two vectors.
c
c     uses unrolled loops for increments equal to one.
c     jack dongarra, linpack, 3/11/78.
c     modified 12/3/93, array(1) declarations changed to array(*)
c
      implicit none

      double precision dx(*),dy(*),dtemp
      integer i,incx,incy,ix,iy,m,mp1,n

      ddot = 0.0d0
      dtemp = 0.0d0
      if(n.le.0)return
      if(incx.eq.1.and.incy.eq.1)go to 20
c
c        code for unequal increments or equal increments
c          not equal to 1
c
      ix = 1
      iy = 1
      if(incx.lt.0)ix = (-n+1)*incx + 1
      if(incy.lt.0)iy = (-n+1)*incy + 1
      do i = 1,n
        dtemp = dtemp + dx(ix)*dy(iy)
        ix = ix + incx
        iy = iy + incy
      end do
      ddot = dtemp
      return
c
c        code for both increments equal to 1
c
c
c        clean-up loop
c
   20 m = mod(n,5)
      if( m .eq. 0 ) go to 40
      do 30 i = 1,m
        dtemp = dtemp + dx(i)*dy(i)
   30 continue
      if( n .lt. 5 ) go to 60
   40 mp1 = m + 1
      do 50 i = mp1,n,5
        dtemp = dtemp + dx(i)*dy(i) + dx(i + 1)*dy(i + 1) +
     *   dx(i + 2)*dy(i + 2) + dx(i + 3)*dy(i + 3) + dx(i + 4)*dy(i + 4)
   50 continue
   60 ddot = dtemp
      return
      end
      subroutine delete (filnam)

c*********************************************************************72
c
cc DELETE deletes a file.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    23 January 2007
c
c  Author:
c
c    John Burkardt
c
      implicit none

      character filnam*(*)

      open (unit=99,file=filnam,status='old',err=10)
      close (unit=99,status='delete',err=10)
10    continue
      return
      end
      subroutine dgbfa(abd,lda,n,ml,mu,ipvt,info)

c*********************************************************************72
c
cc DGBFA factors a double precision band matrix.
c
c     dgbfa factors a double precision band matrix by elimination.
c
c     dgbfa is usually called by dgbco, but it can be called
c     directly with a saving in time if  rcond  is not needed.
c
c     on entry
c
c        abd     double precision(lda, n)
c                contains the matrix in band storage.  the columns
c                of the matrix are stored in the columns of  abd  and
c                the diagonals of the matrix are stored in rows
c                ml+1 through 2*ml+mu+1 of  abd .
c                see the comments below for details.
c
c        lda     integer
c                the leading dimension of the array  abd .
c                lda must be .ge. 2*ml + mu + 1 .
c
c        n       integer
c                the order of the original matrix.
c
c        ml      integer
c                number of diagonals below the main diagonal.
c                0 .le. ml .lt. n .
c
c        mu      integer
c                number of diagonals above the main diagonal.
c                0 .le. mu .lt. n .
c                more efficient if  ml .le. mu .
c     on return
c
c        abd     an upper triangular matrix in band storage and
c                the multipliers which were used to obtain it.
c                the factorization can be written  a = l*u  where
c                l  is a product of permutation and unit lower
c                triangular matrices and  u  is upper triangular.
c
c        ipvt    integer(n)
c                an integer vector of pivot indices.
c
c        info    integer
c                = 0  normal value.
c                = k  if  u(k,k) .eq. 0.0 .  this is not an error
c                     condition for this subroutine, but it does
c                     indicate that dgbsl will divide by zero if
c                     called.  use  rcond  in dgbco for a reliable
c                     indication of singularity.
c
c     band storage
c
c           if  a  is a band matrix, the following program segment
c           will set up the input.
c
c                   ml = (band width below the diagonal)
c                   mu = (band width above the diagonal)
c                   m = ml + mu + 1
c                   do 20 j = 1, n
c                      i1 = max0(1, j-mu)
c                      i2 = min0(n, j+ml)
c                      do 10 i = i1, i2
c                         k = i - j + m
c                         abd(k,j) = a(i,j)
c                10    continue
c                20 continue
c
c           this uses rows  ml+1  through  2*ml+mu+1  of  abd .
c           in addition, the first  ml  rows in  abd  are used for
c           elements generated during the triangularization.
c           the total number of rows needed in  abd  is  2*ml+mu+1 .
c           the  ml+mu by ml+mu  upper left triangle and the
c           ml by ml  lower right triangle are not referenced.
c
c     linpack. this version dated 08/14/78 .
c     cleve moler, university of new mexico, argonne national lab.
c
      implicit none

      integer lda,n,ml,mu,ipvt(1),info
      double precision abd(lda,1)
      double precision t
      integer i,idamax,i0,j,ju,jz,j0,j1,k,kp1,l,lm,m,mm,nm1

      m = ml + mu + 1
      info = 0
c
c  zero initial fill-in columns
c
      j0 = mu + 2
      j1 = min0(n,m) - 1
      if (j1 .lt. j0) go to 30
      do 20 jz = j0, j1
         i0 = m + 1 - jz
         do 10 i = i0, ml
            abd(i,jz) = 0.0d0
   10    continue
   20 continue
   30 continue
      jz = j1
      ju = 0
c
c  gaussian elimination with partial pivoting
c
      nm1 = n - 1
      if (nm1 .lt. 1) go to 130
      do 120 k = 1, nm1
         kp1 = k + 1
c
c  zero next fill-in column
c
         jz = jz + 1
         if (jz .gt. n) go to 50
         if (ml .lt. 1) go to 50
            do 40 i = 1, ml
               abd(i,jz) = 0.0d0
   40       continue
   50    continue
c
c        find l = pivot index
c
         lm = min0(ml,n-k)
         l = idamax(lm+1,abd(m,k),1) + m - 1
         ipvt(k) = l + k - m
c
c        zero pivot implies this column already triangularized
c
         if (abd(l,k) .eq. 0.0d0) go to 100
c
c           interchange if necessary
c
            if (l .eq. m) go to 60
               t = abd(l,k)
               abd(l,k) = abd(m,k)
               abd(m,k) = t
   60       continue
c
c           compute multipliers
c
            t = -1.0d0/abd(m,k)
            call dscal(lm,t,abd(m+1,k),1)
c
c           row elimination with column indexing
c
            ju = min0(max0(ju,mu+ipvt(k)),n)
            mm = m
            if (ju .lt. kp1) go to 90
            do 80 j = kp1, ju
               l = l - 1
               mm = mm - 1
               t = abd(l,j)
               if (l .eq. mm) go to 70
                  abd(l,j) = abd(mm,j)
                  abd(mm,j) = t
   70          continue
               call daxpy(lm,t,abd(m+1,k),1,abd(mm+1,j),1)
   80       continue
   90       continue
         go to 110
  100    continue
            info = k
  110    continue
  120 continue
  130 continue
      ipvt(n) = n
      if (abd(m,n) .eq. 0.0d0) info = n
      return
      end
      subroutine dgbsl(abd,lda,n,ml,mu,ipvt,b,job)

c*********************************************************************72
c
cc DGBSL solves a linear system given a matrix factored by DGBFA.
c
c
c     dgbsl solves the double precision band system
c     a * x = b  or  trans(a) * x = b
c     using the factors computed by dgbco or dgbfa.
c
c     on entry
c
c        abd     double precision(lda, n)
c                the output from dgbco or dgbfa.
c
c        lda     integer
c                the leading dimension of the array  abd .
c
c        n       integer
c                the order of the original matrix.
c
c        ml      integer
c                number of diagonals below the main diagonal.
c
c        mu      integer
c                number of diagonals above the main diagonal.
c
c        ipvt    integer(n)
c                the pivot vector from dgbco or dgbfa.
c
c        b       double precision(n)
c                the right hand side vector.
c
c        job     integer
c                = 0         to solve  a*x = b ,
c                = nonzero   to solve  trans(a)*x = b , where
c                            trans(a)  is the transpose.
c
c     on return
c
c        b       the solution vector  x .
c
c     error condition
c
c        a division by zero will occur if the input factor contains a
c        zero on the diagonal.  technically this indicates singularity
c        but it is often caused by improper arguments or improper
c        setting of lda .  it will not occur if the subroutines are
c        called correctly and if dgbco has set rcond .gt. 0.0
c        or dgbfa has set info .eq. 0 .
c
c     to compute  inverse(a) * c  where  c  is a matrix
c     with  p  columns
c           call dgbco(abd,lda,n,ml,mu,ipvt,rcond,z)
c           if (rcond is too small) go to ...
c           do 10 j = 1, p
c              call dgbsl(abd,lda,n,ml,mu,ipvt,c(1,j),0)
c        10 continue
c
c     linpack. this version dated 08/14/78 .
c     cleve moler, university of new mexico, argonne national lab.
c
      implicit none

      integer lda,n,ml,mu,ipvt(1),job
      double precision abd(lda,1),b(1)
      double precision ddot,t
      integer k,kb,l,la,lb,lm,m,nm1
c
      m = mu + ml + 1
      nm1 = n - 1
      if (job .ne. 0) go to 50
c
c        job = 0 , solve  a * x = b
c        first solve l*y = b
c
         if (ml .eq. 0) go to 30
         if (nm1 .lt. 1) go to 30
            do 20 k = 1, nm1
               lm = min0(ml,n-k)
               l = ipvt(k)
               t = b(l)
               if (l .eq. k) go to 10
                  b(l) = b(k)
                  b(k) = t
   10          continue
               call daxpy(lm,t,abd(m+1,k),1,b(k+1),1)
   20       continue
   30    continue
c
c        now solve  u*x = y
c
         do 40 kb = 1, n
            k = n + 1 - kb
            b(k) = b(k)/abd(m,k)
            lm = min0(k,m) - 1
            la = m - lm
            lb = k - lm
            t = -b(k)
            call daxpy(lm,t,abd(la,k),1,b(lb),1)
   40    continue
      go to 100
   50 continue
c
c        job = nonzero, solve  trans(a) * x = b
c        first solve  trans(u)*y = b
c
         do 60 k = 1, n
            lm = min0(k,m) - 1
            la = m - lm
            lb = k - lm
            t = ddot(lm,abd(la,k),1,b(lb),1)
            b(k) = (b(k) - t)/abd(m,k)
   60    continue
c
c        now solve trans(l)*x = y
c
         if (ml .eq. 0) go to 90
         if (nm1 .lt. 1) go to 90
            do 80 kb = 1, nm1
               k = n - kb
               lm = min0(ml,n-k)
               b(k) = b(k) + ddot(lm,abd(m+1,k),1,b(k+1),1)
               l = ipvt(k)
               if (l .eq. k) go to 70
                  t = b(l)
                  b(l) = b(k)
                  b(k) = t
   70          continue
   80       continue
   90    continue
  100 continue
      return
      end
      subroutine dscal ( n, da, dx, incx )

c*********************************************************************72
c
cc DSCAL scales a vector by a constant.
c
c     uses unrolled loops for increment equal to one.
c     jack dongarra, linpack, 3/11/78.
c     modified 3/93 to return if incx .le. 0.
c     modified 12/3/93, array(1) declarations changed to array(*)
c
      implicit none

      double precision da
      double precision dx(*)
      integer i
      integer incx
      integer m
      integer n
      integer nincx

      if ( n .le. 0 ) then
        return
      end if

      if ( incx .le. 0 ) then
        return
      end if

      if ( incx .eq. 1 )then

        m = mod ( n, 5 )

        do i = 1, m
          dx(i) = da * dx(i)
        end do

        do i = m+1, n, 5
          dx(i)   = da * dx(i)
          dx(i+1) = da * dx(i+1)
          dx(i+2) = da * dx(i+2)
          dx(i+3) = da * dx(i+3)
          dx(i+4) = da * dx(i+4)
        end do

      else

        nincx = n * incx
        do i = 1, nincx, incx
          dx(i) = da * dx(i)
        end do

      end if

      return
      end
      subroutine gdump ( f, indx, insc, iounit, isotri, long, nelemn,
     &  neqn, nnodes, node, np, npara, nx, ny, para, reynld, rjpnew,
     &  xc, yc )

c*********************************************************************72
c
cc GDUMP writes information to a file.
c
c  Discussion:
c
c    The information can be used to create
c    graphics images.  In order to keep things simple, exactly one
c    value, real or integer, is written per record.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    23 January 2007
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer NPARA, the number of parameters.  Fixed at 1
c    for now.
c
c    Input, DOUBLE PRECISION PARA(MAXPAR), the parameters.
c
      implicit none

      integer nelemn
      integer neqn
      integer nnodes
      integer np

      double precision f(neqn)
      double precision fval
      integer i
      integer indx(np,2)
      integer insc(np)
      integer iounit
      integer iset
      integer isotri(nelemn)
      integer j
      logical long
      integer node(nelemn,nnodes)
      integer npara
      integer nx
      integer ny
      double precision para
      double precision reynld
      double precision rjpnew
      double precision ubdry
      double precision xc(np)
      double precision yc(np)

      save iset

      data iset / 0 /

      iset = iset + 1

      write ( iounit, * ) long
      write ( iounit, * ) nelemn
      write ( iounit, * ) np
      write ( iounit, * ) npara
      write ( iounit, * ) nx
      write ( iounit, * ) ny
c
c  Pressures
c
      do i = 1, np
        j = insc(i)
        if (j.le.0) then
          fval = 0.0D+00
        else
          fval = f(j)
        end if
        write (iounit,*) fval
      end do
c
c  Horizontal velocities, U
c
      do i = 1, np
        j = indx(i,1)
        if (j.eq.0) then
          fval = 0.0D+00
        elseif (j.lt.0) then
          fval=ubdry(yc(i),para)
        else
          fval = f(j)
        end if
        write (iounit,*) fval
      end do
c
c  Vertical velocities, V
c
      do i = 1, np
        j = indx(i,2)
        if (j.le.0) then
          fval = 0.0D+00
        else
          fval = f(j)
        end if
        write (iounit,*) fval
      end do

      do i = 1, np
        write (iounit,*) indx(i,1)
        write (iounit,*) indx(i,2)
      end do

      do i = 1, np
        write (iounit,*) insc(i)
      end do

      do i = 1, nelemn
        write (iounit,*) isotri(i)
      end do

      do i = 1, nelemn
        do j = 1, 6
          write (iounit,*) node(i,j)
        end do
      end do

      write (iounit,*) para

      write (iounit,*) reynld

      write (iounit,*) rjpnew

      do i = 1, np
        write (iounit,*) xc(i)
      end do

      do i = 1, np
        write (iounit,*) yc(i)
      end do

      write (*,*) 'GDUMP wrote data set ',iset,' to file.'

      return
      end
      subroutine getg ( f, iline, my, neqn, u )

c*********************************************************************72
c
cc GETG outputs field values along the line X=XZERO.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    23 January 2007
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer neqn
      integer my

      double precision f(neqn)
      integer iline(my)
      integer j
      integer k
      double precision u(my)

      do j = 1, my

        k = iline(j)

        if ( 0 .lt. k ) then
          u(j) = f(k)
        else
          u(j) = 0.0D+00
        end if

      end do

      return
      end
      subroutine gram ( gr, iline, indx, iwrite, my, nelemn, nnodes,
     &  node, nodex0, np, para, r, ui, xc, yc )

c*********************************************************************72
c
cc GRAM computes the Gram matrix.
c
c  Discussion:
c
c    The routine computes:
c
c      the Gram matrix, GR(I,J) = INTEGRAL PHI(I)*PHI(J)
c
c      the vector R(I) = INTEGRAL UI*PHI(I).
c
c    The integrals are computed along the line where the profile is
c    specified.  The three point Gauss quadrature rule is used for the
c    line integral.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    23 January 2007
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer my
      integer nelemn
      integer nnodes
      integer np

      double precision ar
      double precision bb
      double precision bbb
      double precision bbx
      double precision bby
      double precision bma2
      double precision bx
      double precision by
      double precision gr(my,my)
      integer i
      integer igetl
      integer ii
      integer iline(my)
      integer indx(np,2)
      integer ip
      integer ipp
      integer iq
      integer iqq
      integer iquad
      integer it
      integer iun
      integer iwrite
      integer j
      integer jj
      integer k
      integer kk
      integer node(nelemn,nnodes)
      integer nodex0
      double precision para
      double precision r(my)
      double precision ubc
      double precision ubdry
      double precision ui(my)
      double precision uiqdpt
      double precision wt(3)
      double precision x
      double precision xc(np)
      double precision xzero
      double precision y
      double precision yq(3)
      double precision yc(np)
c
c  Input values for 3 point Gauss quadrature
c
      wt(1) = 5.0D+00 / 9.0D+00
      wt(2) = 8.0D+00 / 9.0D+00
      wt(3) = wt(1)
      yq(1) = -0.7745966692D+00
      yq(2) = 0.0D+00
      yq(3) = -yq(1)
c
c  zero arrays
c
      do i = 1, my
        r(i) = 0.0D+00
        do j = 1, my
          gr(i,j) = 0.0D+00
        end do
      end do
c
c  Compute line integral by looping over intervals along line
c  using three point Gauss quadrature
c
      xzero = xc(nodex0)
      do 70 it = 1, nelemn
c
c  Check to see if we are in a triangle with a side along line
c  x=xzero.  If not, skip out
c
        k = node(it,1)
        kk = node(it,2)
        if (abs(xc(k)-xzero).gt.1.0D-04) go to 70
        if (abs(xc(kk)-xzero).gt.1.0D-04) go to 70

        do 60 iquad = 1, 3

          bma2 = (yc(kk)-yc(k))/2.0D+00
          ar = bma2*wt(iquad)
          x = xzero
          y = yc(k)+bma2*(yq(iquad)+1.0D+00)
c
c  Compute u internal at quadrature points
c
          uiqdpt = 0.0D+00

          do 30 iq = 1, nnodes

            if (iq.gt.4) go to 30
            if (iq.eq.3) go to 30
            call qbf (x,y,it,iq,bb,bx,by,nelemn,nnodes,node,np,xc,yc)
            ip = node(it,iq)
            iun = indx(ip,1)

            if (iun.gt.0) then
              ii = igetl(iun,iline,my)
              uiqdpt = uiqdpt+bb*ui(ii)
            elseif (iun.lt.0) then
              ubc = ubdry(yc(ip),para)
              uiqdpt = uiqdpt+bb*ubc
            end if

   30     continue
c
c  Only loop over nodes lying on line x=xzero
c
          do 50 iq = 1, nnodes

            if(iq.eq.1.or.iq.eq.2.or.iq.eq.4)then
              ip = node(it,iq)
              call qbf (x,y,it,iq,bb,bx,by,nelemn,nnodes,node,np,xc,yc)
              i = indx(ip,1)
              if (i.le.0) go to 50
              ii = igetl(i,iline,my)
              r(ii) = r(ii)+bb*uiqdpt*ar

              do iqq = 1, nnodes
                if(iqq.eq.1.or.iqq.eq.2.or.iqq.eq.4)then
                  ipp = node(it,iqq)
                  call qbf (x,y,it,iqq,bbb,bbx,bby,nelemn,nnodes, 
     &                node,np,xc,yc)
                  j = indx(ipp,1)
                  if (j.ne.0) then
                    jj = igetl(j,iline,my)
                    gr(ii,jj) = gr(ii,jj)+bb*bbb*ar
                  end if
                end if
              end do

            end if
   50     continue
   60   continue
   70 continue

      if(iwrite.ge.2)then

        write(*,*)' '
        write(*,*)'Gram matrix:'
        write(*,*)' '

        do i=1,my
          do j=1,my
            write(*,*)i,j,gr(i,j)
          end do
        end do

        write(*,*)' '
        write(*,*)'R vector:'
        write(*,*)' '

        do i=1,my
          write(*,*)i,r(i)
        end do

      end if

      return
      end
      function igetl ( k, iline, my )

c*********************************************************************72
c
cc IGETL gets the local index for nodes along the line X=XZERO.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    23 January 2007
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer my

      integer igetl
      integer iline(my)
      integer j
      integer k

      do j = 1, my

        if ( iline(j) .eq. k ) then
          igetl = j
          return
        end if

      end do

      write (*,*) ' '
      write (*,*) 'IGETL - fatal error!'
      write (*,*) '  Unable to get local unknown number for '
      write (*,*) '  global variable number ',k

      igetl = 0

      stop
      end
      function idamax ( n, dx, incx )

c*********************************************************************72
c
cc IDAMAX finds the index of element having max. absolute value.
c
c     jack dongarra, linpack, 3/11/78.
c     modified 3/93 to return if incx .le. 0.
c     modified 12/3/93, array(1) declarations changed to array(*)
c
      implicit none

      double precision dmax
      double precision dx(*)
      integer i
      integer idamax
      integer incx
      integer ix
      integer n

      idamax = 0

      if ( n.lt.1 .or. incx .le. 0 ) then
        return
      end if

      idamax = 1

      if ( n .eq. 1 ) then
        return
      end if

      if ( incx .eq. 1 ) then

        dmax = dabs ( dx(1) )

        do i = 2, n

          if ( dmax .lt. dabs ( dx(i) ) ) then
            idamax = i
            dmax = dabs ( dx(i) )
          end if

        end do

      else

        ix = 1
        dmax = dabs ( dx(1) )
        ix = ix + incx

        do i = 2, n

          if ( dmax < dabs ( dx(ix) ) ) then
            idamax = i
            dmax = dabs ( dx(ix) )
          end if

          ix = ix + incx

        end do

      end if

      end
      subroutine linsys ( a, area, f, g, indx, insc, ipivot, 
     &  maxrow, nelemn, neqn, nlband, nnodes, node, 
     &  np, nquad, nrow, para1, para2, phi, psi, reynld, yc )

c*********************************************************************72
c
cc LINSYS sets up and solves the linear system.
c
c  Discussion:
c
c    The G array contains the previous solution.
c
c    The F array contains the right hand side initially and then the
c    current solution.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    28 February 2006
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer maxrow
      integer nelemn
      integer neqn
      integer nnodes
      integer np
      integer nquad
      integer nrow

      double precision a(maxrow,neqn)
      double precision ar
      double precision area(nelemn)
      double precision bb
      double precision bbb
      double precision bbbl
      double precision bbl
      double precision bbx
      double precision bby
      double precision bx
      double precision by
      double precision f(neqn)
      double precision g(neqn)
      integer i
      integer ihor
      integer indx(np,2)
      integer info
      integer insc(np)
      integer ioff
      integer ip
      integer ipivot(neqn)
      integer ipp
      integer iprs
      integer iq
      integer iqq
      integer iquad
      integer it
      integer iuse
      integer iver
      integer j
      integer job
      integer jp
      integer ju
      integer jv
      integer nlband
      integer node(nelemn,nnodes)
      double precision para1
      double precision para2
      double precision phi(nelemn,nquad,nnodes,3)
      double precision psi(nelemn,nquad,nnodes)
      double precision reynld
      double precision ubdry
      double precision un(2)
      double precision unx(2)
      double precision uny(2)
      double precision uu
      double precision visc
      double precision yc(np)

      ioff = nlband + nlband + 1
      visc = 1.0D+00 / reynld

      do i = 1, neqn
        f(i) = 0.0D+00
      end do

      do i = 1, nrow
        do j = 1, neqn
          a(i,j) = 0.0D+00
        end do
      end do
c
c  For each element,
c
      do it = 1, nelemn

        ar = area(it) / 3.0D+00
c
c  and for each quadrature point in the element,
c
        do iquad = 1, nquad
c
c  Evaluate the velocities at the quadrature point
c
          call uval ( g, indx, iquad, it, nelemn, neqn, nnodes, node, 
     &       np, nquad, para1, phi, un, unx, uny, yc )
c
c  For each basis function,
c
          do iq = 1, nnodes

            ip = node(it,iq)
            bb = phi(it,iquad,iq,1)
            bx = phi(it,iquad,iq,2)
            by = phi(it,iquad,iq,3)
            bbl = psi(it,iquad,iq)
            ihor = indx(ip,1)
            iver = indx(ip,2)
            iprs = insc(ip)

            if(ihor.gt.0)then
              f(ihor) = f(ihor)+ar*bb*(un(1)*unx(1)+un(2)*uny(1))
            end if

            if(iver.gt.0)then
              f(iver) = f(iver)+ar*bb*(un(1)*unx(2)+un(2)*uny(2))
            end if
c
c  For another basis function,
c
            do iqq = 1, nnodes

              ipp = node(it,iqq)
              bbb = phi(it,iquad,iqq,1)
              bbx = phi(it,iquad,iqq,2)
              bby = phi(it,iquad,iqq,3)
              bbbl = psi(it,iquad,iqq)
              ju = indx(ipp,1)
              jv = indx(ipp,2)
              jp = insc(ipp)
c
c  Horizontal velocity variable
c
              if(ju.gt.0)then

                if(ihor.gt.0)then
                  iuse = ihor-ju+ioff
                  a(iuse,ju) = a(iuse,ju) + ar * (
     &               visc*(by*bby+bx*bbx) 
     &               + bb*(bbb*unx(1)+bbx*un(1)+bby*un(2)) )
                end if

                if(iver.gt.0)then
                  iuse = iver-ju+ioff
                  a(iuse,ju) = a(iuse,ju)+ar*bb*bbb*unx(2)
                end if

                if(iprs.gt.0)then
                  iuse = iprs-ju+ioff
                  a(iuse,ju) = a(iuse,ju)+ar*bbx*bbl
                end if

              else if ( ju .lt. 0 ) then

                uu = ubdry(yc(ipp),para2)
                if ( 0 .lt. ihor ) then 
                  f(ihor) = f(ihor)-ar*uu*(visc*(by*bby+bx*bbx) 
     &               + bb*(bbb*unx(1)+bbx*un(1)+bby*un(2)))
                end if

                if ( 0 .lt. iver ) then
                  f(iver) = f(iver)-ar*uu*bb*bbb*unx(2)
                end if

                if ( 0 .lt. iprs ) then
                  f(iprs) = f(iprs)-ar*uu*bbx*bbl
                end if

              end if
c
c  Vertical velocity variable
c
              if(jv.gt.0)then

                if(ihor.gt.0)then
                  iuse = ihor-jv+ioff
                  a(iuse,jv) = a(iuse,jv)+ar*bb*bbb*uny(1)
                end if

                if(iver.gt.0)then
                  iuse = iver-jv+ioff
                  a(iuse,jv) = a(iuse,jv)+ar*(visc*(by*bby+bx*bbx) 
     &               +bb*(bbb*uny(2)+bby*un(2)+bbx*un(1)))
                end if

                if(iprs.gt.0)then
                  iuse = iprs-jv+ioff
                  a(iuse,jv) = a(iuse,jv)+ar*bby*bbl
                end if

              end if
c
c  Pressure variable
c
              if(jp.gt.0)then

                if(ihor.gt.0)then
                  iuse = ihor-jp+ioff
                  a(iuse,jp) = a(iuse,jp)-ar*bx*bbbl
                end if

                if(iver.gt.0)then
                  iuse = iver-jp+ioff
                  a(iuse,jp) = a(iuse,jp)-ar*by*bbbl
                end if

              end if

            end do
          end do
        end do
      end do
c
c  To avoid singularity of the pressure system, the last pressure
c  is simply assigned a value of 0.
c
      f(neqn) = 0.0D+00
      do j = neqn-nlband, neqn-1
        i=neqn-j+ioff
        a(i,j)=0.0D+00
      end do
      a(ioff,neqn) = 1.0D+00
c
c  Factor the matrix
c
      call dgbfa ( a, maxrow, neqn, nlband, nlband, ipivot, info )

      if ( info .ne. 0 ) then
        write (*,*) ' '
        write (*,*) 'LINSYS - fatal error!'
        write (*,*) 'DGBFA returns INFO=',info
        stop
      end if
c
c  Solve the linear system
c
      job = 0
      call dgbsl ( a, maxrow, neqn, nlband, nlband, ipivot, f, job )

      if (info.ne.0) then
        write (*,*) ' '
        write (*,*) 'LINSYS - fatal error!'
        write (*,*) 'DGBSL returns INFO=',info
        stop
      end if

      return
      end
      subroutine nstoke ( a, area, f, g, indx, insc, ipivot, iwrite, 
     &   maxnew, maxrow, nelemn, neqn, nlband, nnodes, node, 
     &   np, nquad, nrow, numnew, para, phi, psi, reynld, tolnew, yc )

c*********************************************************************72
c
cc NSTOKE solves the Navier Stokes equation using Taylor-Hood elements.
c
c  Discussion:
c
c    The G array contains the previous iterate.
c
c    The F array contains the right hand side initially and then the
c    current iterate.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    23 January 2007
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer maxrow
      integer nelemn
      integer neqn
      integer nnodes
      integer np
      integer nquad
      integer nrow

      double precision a(maxrow,neqn)
      double precision area(nelemn)
      double precision diff
      double precision f(neqn)
      double precision g(neqn)
      integer i
      integer idamax
      integer indx(np,2)
      integer insc(np)
      integer ipivot(neqn)
      integer iter
      integer iwrite
      integer maxnew
      integer nlband
      integer node(nelemn,nnodes)
      integer numnew
      double precision para
      double precision phi(nelemn,nquad,nnodes,3)
      double precision psi(nelemn,nquad,nnodes)
      double precision reynld
      double precision tolnew
      double precision yc(np)

      do iter = 1, maxnew

        numnew = numnew + 1

        call linsys (a,area,f,g,indx,insc,ipivot, 
     &     maxrow,nelemn,neqn,nlband,nnodes,node, 
     &     np,nquad,nrow,para,para,phi,psi,reynld,yc)
c
c  Check for convergence
c
        do i = 1, neqn
          g(i) = g(i) - f(i)
        end do

        diff = abs ( g(idamax(neqn,g,1)) )

        if ( 1 .le. iwrite ) then
          write(*,*)'NSTOKE iteration ',iter,' Mnorm=',diff
        end if

        do i = 1, neqn
          g(i) = f(i)
        end do

        if ( diff .le. tolnew ) then
          write (*,*) 'Navier Stokes iteration converged in ', 
     &       iter,' iterations.'
          return
        end if

      end do

      write (*,*) 'Navier Stokes solution did not converge!'

      return
      end
      subroutine pval (g,insc,long,mx,my,nelemn,neqn,nnodes,node,
     &  np,press)

c*********************************************************************72
c
cc PVAL computes a table of pressures.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    23 January 2007
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer nelemn
      integer neqn
      integer nnodes
      integer np

      double precision g(neqn)
      integer i
      integer insc(np)
      integer ip
      integer iq
      integer it
      integer ivar
      integer j
      logical long
      integer mx
      integer my
      integer node(nelemn,nnodes)
      double precision press(mx,my)

      do i=1,mx
        do j=1,my
          press(i,j) = 0.0D+00
        end do
      end do
c
c  Read the pressures where they are computed.  
c  These are "(odd, odd)" points.
c
      do it = 1, nelemn
        do iq = 1, 3

          ip=node(it,iq)
          ivar=insc(ip)

          if(long)then
            i=((ip-1)/my)+1
            j=mod(ip-1,my)+1
          else
            i=mod(ip-1,mx)+1
            j=((ip-1)/mx)+1
          end if

          if(ivar.gt.0)then
            press(i,j)=g(ivar)
          else
            press(i,j)=0.0
          end if

        end do
      end do
c
c  Interpolate the pressures at points (even, odd) and (odd, even).
c
      do i=2,mx-1,2
        do j=1,my,2
          press(i,j)=0.5*(press(i-1,j)+press(i+1,j))
        end do
      end do

      do j=2,my-1,2
        do i=1,mx,2
          press(i,j)=0.5*(press(i,j-1)+press(i,j+1))
        end do
      end do
c
c  Interpolate the pressures at points (even,even).
c
      do j=2,my-1,2
        do i=2,mx-1,2
          press(i,j)=0.5*(press(i-1,j-1)+press(i+1,j+1))
        end do
      end do

      return
      end
      subroutine qbf (x,y,it,in,bb,bx,by,nelemn,nnodes,node,np,xc,yc)

c*********************************************************************72
c
cc QBF evaluates a quadratic basis function in a triangle.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    23 January 2007
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer nelemn
      integer nnodes
      integer np

      double precision bb
      double precision bx
      double precision by
      double precision c
      double precision d
      integer i1
      integer i2
      integer i3
      integer in
      integer in1
      integer in2
      integer in3
      integer inn
      integer it
      integer j1
      integer j2
      integer j3
      integer node(nelemn,nnodes)
      double precision s
      double precision t
      double precision x
      double precision xc(np)
      double precision y
      double precision yc(np)

      if (in.le.3) then

        in1 = in
        in2 = mod(in,3)+1
        in3 = mod(in+1,3)+1

        i1 = node(it,in1)
        i2 = node(it,in2)
        i3 = node(it,in3)

        d = (xc(i2)-xc(i1))*(yc(i3)-yc(i1))
     &    -(xc(i3)-xc(i1))*(yc(i2)-yc(i1))
        t = 1.0+((yc(i2)-yc(i3))*(x-xc(i1))
     &    +(xc(i3)-xc(i2))*(y-yc(i1)))/d

        bb = t*(2.0*t-1.0)
        bx = (yc(i2)-yc(i3))*(4.0*t-1.0)/d
        by = (xc(i3)-xc(i2))*(4.0*t-1.0)/d

      else

        inn = in-3
        in1 = inn
        in2 = mod(inn,3)+1
        in3 = mod(inn+1,3)+1

        i1 = node(it,in1)
        i2 = node(it,in2)
        i3 = node(it,in3)

        j1 = i2
        j2 = i3
        j3 = i1

        d = (xc(i2)-xc(i1))*(yc(i3)-yc(i1))
     &    -(xc(i3)-xc(i1))*(yc(i2)-yc(i1))
        c = (xc(j2)-xc(j1))*(yc(j3)-yc(j1))
     &    -(xc(j3)-xc(j1))*(yc(j2)-yc(j1))
        t = 1.0 + ((yc(i2)-yc(i3))*(x-xc(i1))
     &    +(xc(i3)-xc(i2))*(y-yc(i1)))/d
        s = 1.0 + ((yc(j2)-yc(j3))*(x-xc(j1))
     &    +(xc(j3)-xc(j2))*(y-yc(j1)))/c

        bb = 4.0*s*t
        bx = 4.0*(t*(yc(j2)-yc(j3))/c+s*(yc(i2)-yc(i3))/d)
        by = 4.0*(t*(xc(j3)-xc(j2))/c+s*(xc(i3)-xc(i2))/d)

      end if

      return
      end
      subroutine resid (area,g,indx,insc,iwrite,nelemn,neqn,nnodes, 
     &   node,np,nquad,para,phi,psi,res,reynld,yc)

c*********************************************************************72
c
cc RESID computes the residual.
c
c  Discussion:
c
c    The G array contains the current iterate.
c    The RES array will contain the value of the residual.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    23 January 2007
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer nelemn
      integer neqn
      integer nnodes
      integer np
      integer nquad

      double precision aijpu
      double precision aijpv
      double precision aijup
      double precision aijuu
      double precision aijuv
      double precision aijvp
      double precision aijvu
      double precision aijvv
      double precision ar
      double precision area(nelemn)
      double precision bb
      double precision bbb
      double precision bbbl
      double precision bbl
      double precision bbx
      double precision bby
      double precision bx
      double precision by
      double precision g(neqn)
      integer i
      integer ibad
      integer ihor
      integer imax
      integer indx(np,2)
      integer insc(np)
      integer ip
      integer ipp
      integer iprs
      integer iq
      integer iqq
      integer iquad
      integer it
      integer iver
      integer iwrite
      integer j
      integer jp
      integer ju
      integer jv
      integer node(nelemn,nnodes)
      double precision para
      double precision phi(nelemn,nquad,nnodes,3)
      double precision psi(nelemn,nquad,nnodes)
      double precision res(neqn)
      double precision reynld
      double precision rmax
      double precision test
      double precision ubdry
      double precision un(2)
      double precision unx(2)
      double precision uny(2)
      double precision uu
      double precision visc
      double precision yc(np)
c
c  zero arrays
c
      visc = 1.0D+00 / reynld
      do i = 1, neqn
        res(i) = 0.0D+00
      end do
c
c  For each element,
c
      do it = 1, nelemn

        ar = area(it) / 3.0D+00
c
c  and for each quadrature point in the element,
c
        do iquad = 1, nquad
c
c  Evaluate velocities at quadrature point
c
          call uval (g,indx,iquad,it,nelemn,neqn,nnodes,node, 
     &       np,nquad,para,phi,un,unx,uny,yc)
c
c  For each basis function,
c
          do iq = 1, nnodes

            ip = node(it,iq)
            bb=phi(it,iquad,iq,1)
            bx=phi(it,iquad,iq,2)
            by=phi(it,iquad,iq,3)
            bbl = psi(it,iquad,iq)
            iprs = insc(ip)
            ihor=indx(ip,1)
            iver=indx(ip,2)

            if(ihor.gt.0)then
              res(ihor) = res(ihor)+(un(1)*unx(1)+un(2)*uny(1))*bb*ar
            end if

            if(iver.gt.0)then
              res(iver) = res(iver)+(un(1)*unx(2)+un(2)*uny(2))*bb*ar
            end if
c
c  For another basis function,
c
            do iqq = 1, nnodes

              ipp = node(it,iqq)
              bbb = phi(it,iquad,iqq,1)
              bbx = phi(it,iquad,iqq,2)
              bby = phi(it,iquad,iqq,3)
              bbbl = psi(it,iquad,iqq)
              ju = indx(ipp,1)
              jv = indx(ipp,2)
              jp = insc(ipp)
c
              if(ju.gt.0)then
                if(ihor.gt.0)then
                  aijuu = visc*(by*bby+bx*bbx) 
     &               + bb*(bbb*unx(1)+bbx*un(1)+bby*un(2))
                  res(ihor)=res(ihor)+aijuu*ar*g(ju)
                end if
                if(iver.gt.0)then
                  aijvu=bb*bbb*unx(2)
                  res(iver) = res(iver)+aijvu*ar*g(ju)
                end if
                if(iprs.gt.0)then
                  aijpu=bbx*bbl
                  res(iprs) = res(iprs)+aijpu*ar*g(ju)
                end if
              elseif(ju.lt.0)then
                uu=ubdry(yc(ipp),para)
                if(ihor.gt.0)then
                  aijuu = visc*(by*bby+bx*bbx) 
     &               + bb*(bbb*unx(1)+bbx*un(1)+bby*un(2))
                  res(ihor) = res(ihor)+ar*aijuu*uu
                end if
                if(iver.gt.0)then
                  aijvu=bb*bbb*unx(2)
                  res(iver) = res(iver)+ar*aijvu*uu
                end if
                if(iprs.gt.0)then
                  aijpu=bbx*bbl
                  res(iprs) = res(iprs)+ar*aijpu*uu
                end if
              end if

              if(jv.gt.0)then
                if(ihor.gt.0)then
                  aijuv = bb*bbb*uny(1)
                  res(ihor) = res(ihor)+aijuv*ar*g(jv)
                end if
                if(iver.gt.0)then
                  aijvv = visc*(by*bby+bx*bbx) 
     &               +bb*(bbb*uny(2)+bby*un(2)+bbx*un(1))
                  res(iver)=res(iver)+aijvv*ar*g(jv)
                  end if
                if(iprs.gt.0)then
                  aijpv=bby*bbl
                  res(iprs)=res(iprs)+aijpv*ar*g(jv)
                end if
              end if

              if(jp.gt.0)then
                if(ihor.gt.0)then
                  aijup = -bx*bbbl
                  res(ihor)=res(ihor)+aijup*ar*g(jp)
                end if
                if(iver.gt.0)then
                  aijvp=-by*bbbl
                  res(iver)=res(iver)+aijvp*ar*g(jp)
                end if
              end if

            end do
          end do
        end do
      end do

      res(neqn)=g(neqn)

      rmax=0.0D+00
      imax=0
      ibad=0

      do i=1,neqn

        test=abs(res(i))

        if(test.gt.rmax)then
          rmax=test
          imax=i
        end if

        if(test.gt.1.0D-03)then
          ibad=ibad+1
        end if

      end do

      if(iwrite.ge.1)then
        write(*,*)' '
        write(*,*)'RESIDUAL INFORMATION:'
        write(*,*)' '
        write(*,*)'Worst residual is number ',IMAX
        write(*,*)'of magnitude ',RMAX
        write(*,*)' '
        write(*,*)'Number of "bad" residuals is ',IBAD,' out of ',NEQN 
        write(*,*)' '
      end if

      if(iwrite.ge.2)then
        write(*,*)'Raw residuals:'
        write(*,*)' '
        i=0

        do j=1,np

          if(indx(j,1).gt.0)then
            i=i+1
            if(abs(res(i)).le.1.0D-03)then
              write(*,'(1x,a1,2i5,g14.6)')'U',i,j,res(i)
            else
              write(*,'(a1,a1,2i5,g14.6)')'*','U',i,j,res(i)
            end if
          end if

          if(indx(j,2).gt.0)then
            i=i+1
            if(abs(res(i)).le.1.0D-03)then
              write(*,'(1x,a1,2i5,g14.6)')'V',i,j,res(i)
            else
              write(*,'(a1,a1,2i5,g14.6)')'*','V',i,j,res(i)
            end if
          end if

          if(insc(j).gt.0)then
            i=i+1
            if(abs(res(i)).le.1.0D-03)then
              write(*,'(1x,a1,2i5,g14.6)')'P',i,j,res(i)
            else
              write(*,'(a1,a1,2i5,g14.6)')'*','P',i,j,res(i)
            end if
          end if

        end do

      end if

      return
      end
      subroutine setban(indx,insc,maxrow,nband,nelemn,nlband,nnodes, 
     &   node,np,nrow)

c*********************************************************************72
c
cc SETBAN computes the half band width
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    23 January 2007
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer nelemn
      integer nnodes
      integer np

      integer i
      integer indx(np,2)
      integer insc(np)
      integer ip
      integer ipp
      integer iq
      integer iqq
      integer it
      integer iuk
      integer iukk
      integer j
      integer maxrow
      integer nband
      integer nlband
      integer node(nelemn,nnodes)
      integer nrow

      nlband = 0
      do it = 1, nelemn
        do iq = 1, nnodes
          ip = node(it,iq)
          do iuk = 1, 3

            if (iuk.eq.3) then
              i = insc(ip)
            else
              i = indx(ip,iuk)
            end if

            if (i.gt.0)then
              do iqq = 1, nnodes
                ipp = node(it,iqq)
                do iukk = 1, 3
                  if ( iukk .eq. 3 ) then
                    j = insc(ipp)
                  else
                    j = indx(ipp,iukk)
                  end if
                  nlband = max(nlband,j-i)
                end do
              end do
            end if

          end do
        end do
      end do

      nband = nlband+nlband+1
      nrow = nlband+nlband+nlband+1

      write (*,*) 'Lower bandwidth=',nlband
      write (*,*) 'Total bandwidth=',nband
      write (*,*) 'NROW =',nrow
      if(nrow.gt.maxrow)then
        write(*,*)'SETBAN - NROW is too large!'
        write(*,*)'The maximum allowed is ',maxrow
        stop
      end if

      return
      end
      subroutine setbas(nelemn,nnodes,node,np,nquad,phi,psi,xc,xm,yc,ym)

c*********************************************************************72
c
cc SETBAS computes the basis functions at each integration point.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    23 January 2007
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer nelemn
      integer nnodes
      integer np
      integer nquad

      double precision bb
      double precision bsp
      double precision bx
      double precision by
      integer iq
      integer it
      integer j
      integer node(nelemn,nnodes)
      double precision phi(nelemn,nquad,nnodes,3)
      double precision psi(nelemn,nquad,nnodes)
      double precision x
      double precision xc(np)
      double precision xm(nelemn,nquad)
      double precision y
      double precision yc(np)
      double precision ym(nelemn,nquad)

      do it=1,nelemn
        do j=1,nquad
          x=xm(it,j)
          y=ym(it,j)
          do iq=1,6
            psi(it,j,iq)=bsp(x,y,it,iq,1,nelemn,nnodes,node,np,xc,yc)
            call qbf(x,y,it,iq,bb,bx,by,nelemn,nnodes,node,np,xc,yc)
            phi(it,j,iq,1)=bb
            phi(it,j,iq,2)=bx
            phi(it,j,iq,3)=by
          end do
        end do
      end do

      return
      end
      subroutine setgrd (indx,insc,isotri,iwrite,long,maxeqn,mx,my, 
     &    nelemn,neqn,nnodes,node,np,nx,ny)

c*********************************************************************72
c
cc SETGRD sets up the grid.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    23 January 2007
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer my
      integer nelemn
      integer nnodes
      integer np

      integer i
      integer ic
      integer icnt
      integer ielemn
      integer indx(np,2)
      integer insc(np)
      integer ip
      integer ip1
      integer ip2
      integer isotri(nelemn)
      integer it
      integer iwrite
      integer jc
      integer jcnt
      logical long
      integer maxeqn
      integer mx
      integer neqn
      integer node(nelemn,nnodes)
      integer nx
      integer ny
c
c  Determine whether region is long or skinny.  This will determine
c  how we number the nodes and elements.
c
      if(nx.gt.ny)then
        long=.true.
        write(*,*)'Using vertical ordering.'
      else
        long=.false.
        write(*,*)'Using horizontal ordering.'
      end if
c
c  Set parameters for Taylor Hood element
c
      write (*,*) ' '
      write (*,*) 'SETGRD: Taylor Hood element'
c
c  Construct grid coordinates, elements, and ordering of unknowns
c
      neqn = 0
      ielemn = 0

      do ip=1,np

        if ( long ) then
          ic=((ip-1)/my)+1
          jc=mod((ip-1),my)+1
        else
          ic=mod((ip-1),mx)+1
          jc=((ip-1)/mx)+1
        end if

        icnt = mod(ic,2)
        jcnt = mod(jc,2)
c
c  If both the row count and the column count are odd,
c  and we're not in the last row or top column, 
c  then we can define two new triangular elements based at the node.
c
c  For horizontal ordering, 
c  given the following arrangement of nodes, for instance:
c 
c    21 22 23 24 25
c    16 17 18 19 20
c    11 12 13 14 15
c    06 07 08 09 10
c    01 02 03 04 05
c
c  when we arrive at node 13, we will define 
c
c  element 7: (13, 23, 25, 18, 24, 19)
c  element 8: (13, 25, 15, 19, 20, 14)
c
c
c  For vertical ordering,
c  given the following arrangement of nodes, for instance:
c 
c    05 10 15 20 25
c    04 09 14 19 24
c    03 08 13 18 23
c    02 07 12 17 22
c    01 06 11 16 21
c
c  when we arrive at node 13, we will define 
c
c  element 7: (13, 25, 23, 19, 24, 18)
c  element 8: (13, 15, 25, 14, 20, 19)
c
        if((icnt.eq.1.and.jcnt.eq.1).and. 
     &      (ic.ne.mx).and.(jc.ne.my))then

          if(long)then
            ip1=ip+my
            ip2=ip+my+my
            ielemn=ielemn+1
            node(ielemn,1) = ip
            node(ielemn,2) = ip2+2
            node(ielemn,3) = ip2
            node(ielemn,4) = ip1+1
            node(ielemn,5) = ip2+1
            node(ielemn,6) = ip1
            isotri(ielemn)=0
            ielemn=ielemn+1
            node(ielemn,1) = ip
            node(ielemn,2) = ip+2
            node(ielemn,3) = ip2+2
            node(ielemn,4) = ip+1
            node(ielemn,5) = ip1+2
            node(ielemn,6) = ip1+1
            isotri(ielemn)=0

          else
            ip1=ip+mx
            ip2=ip+mx+mx
            ielemn=ielemn+1
            node(ielemn,1) = ip
            node(ielemn,2) = ip2
            node(ielemn,3) = ip2+2
            node(ielemn,4) = ip1
            node(ielemn,5) = ip2+1
            node(ielemn,6) = ip1+1
            isotri(ielemn)=0
            ielemn=ielemn+1
            node(ielemn,1) = ip
            node(ielemn,2) = ip2+2
            node(ielemn,3) = ip+2
            node(ielemn,4) = ip1+1
            node(ielemn,5) = ip1+2
            node(ielemn,6) = ip+1
            isotri(ielemn)=0
          end if
        end if
c
c  Consider whether velocity unknowns should be associated with this node.
c
c  If we are in column 1, horizontal velocities are specified, and 
c  vertical velocities are zero.
c
        if(ic.eq.1.and.1.lt.jc.and.jc.lt.my)then
          indx(ip,1) = -1
          indx(ip,2) = 0
c
c  If we are in column MX, horizontal velocities are unknown, and
c  vertical velocities are zero.
c
        else if (ic.eq.mx.and.1.lt.jc.and.jc.lt.my)then
          neqn = neqn+1
          indx(ip,1) = neqn
          indx(ip,2) = 0
c
c  Otherwise, if we are in row 1 or row MY, both horizontal and
c  vertical velocities are zero.
c
        else if (jc.eq.1.or.jc.eq.my)then
          indx(ip,1) = 0
          indx(ip,2) = 0    
c
c  Otherwise, we are at an interior node
c
        else
          neqn = neqn+2
          indx(ip,1) = neqn-1
          indx(ip,2) = neqn
        end if
c
c  Consider whether a pressure unknown should be associated with this node.
c  The answer is yes if both nodes are odd.
c
        if (jcnt.eq.1.and.icnt.eq.1) then
          neqn = neqn+1
          insc(ip) = neqn
        else
          insc(ip) = 0
        end if

      end do
c
c  If debugging is requested, print out data.
c
      if (iwrite.ge.2) then

        write(*,*)' '
        write(*,*)'    I      INDX 1 & 2, INSC'
        write(*,*)' '

        do i=1,np
          write (*,'(i5,3i5)') i,indx(i,1),indx(i,2),insc(i)
        end do

        write(*,*)' '
        write(*,*)'    IT    NODE(IT,I),I=1,6)'
        write(*,*)' '

        do it=1,nelemn
          write (*,'(7i6)') it,(node(it,i),i=1,6)
        end do

      end if

      write (*,*) 'Number of unknowns=',neqn
      if (neqn.gt.maxeqn)then
        write(*,*)'SETGRD - Too many unknowns!'
        write(*,*)'The maximum allowed is MAXEQN=',maxeqn
        write(*,*)'This problem requires NEQN=',neqn
        stop
      end if

      return
      end
      subroutine setlin(iline,indx,iwrite,long,mx,my,nodex0,np, 
     &   nx,ny,xlngth)

c*********************************************************************72
c
cc SETLIN determines the unknown numbers along the profile line.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    23 January 2007
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer my
      integer np

      integer i
      integer iline(my)
      integer indx(np,2)
      integer ip
      integer itemp
      integer iwrite
      logical long
      integer mx
      integer nodex0
      integer nx
      integer ny
      double precision xlngth
c
c  Determine the number of a node on the profile line
c
      itemp = nint((2.0D+00*(nx-1)*9.0D+00)/xlngth)

      if ( long ) then
        nodex0 = itemp*(2*ny-1)+1
      else
        nodex0 = itemp+1
      end if

      do i = 1, my
        if ( long ) then
          ip = nodex0+(i-1)
        else
          ip = nodex0+mx*(i-1)
        end if
        iline(i) = indx(ip,1)
      end do

      if ( 1 .le. iwrite ) then
        write(*,*)' '
        write (*,*) 'SETLIN: unknown numbers along line:'
        write(*,*)' '
        write (*,'(1X,15I5)') (iline(i),i=1,my)
        write(*,*)' '
      end if

      return
      end
      subroutine setqud(area,nelemn,nnodes,node,np,nquad,xc,xm,yc,ym)

c*********************************************************************72
c
cc SETQUD sets midpoint quadrature rule information.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    23 January 2007
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer nelemn
      integer nnodes
      integer np
      integer nquad

      double precision area(nelemn)
      integer ip1
      integer ip2
      integer ip3
      integer it
      integer node(nelemn,nnodes)
      double precision x1
      double precision x2
      double precision x3
      double precision xc(np)
      double precision xm(nelemn,nquad)
      double precision y1
      double precision y2
      double precision y3
      double precision yc(np)
      double precision ym(nelemn,nquad)

      do it = 1, nelemn

        ip1 = node(it,1)
        ip2 = node(it,2)
        ip3 = node(it,3)

        x1 = xc(ip1)
        x2 = xc(ip2)
        x3 = xc(ip3)

        y1 = yc(ip1)
        y2 = yc(ip2)
        y3 = yc(ip3)

        xm(it,1) = 0.5*(x1+x2)
        xm(it,2) = 0.5*(x2+x3)
        xm(it,3) = 0.5*(x3+x1)
        ym(it,1) = 0.5*(y1+y2)
        ym(it,2) = 0.5*(y2+y3)
        ym(it,3) = 0.5*(y3+y1)

        area(it) = 0.5*abs((y1+y2)*(x2-x1)+(y2+y3)*(x3-x2) 
     &     +(y3+y1)*(x1-x3))

      end do

      return
      end
      subroutine setxy(iwrite,long,mx,my,np,nx,ny,xc,xlngth,yc,ylngth)

c*********************************************************************72
c
cc SETXY sets the X, Y coordinates of grid points.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    23 January 2007
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer np

      integer i
      integer ic
      integer ip
      integer iwrite
      integer jc
      logical long
      integer mx
      integer my
      integer nx
      integer ny
      double precision xc(np)
      double precision xlngth
      double precision yc(np)
      double precision ylngth
c
c  Construct grid coordinates
c
      do ip = 1, np

        if(long)then
          ic=((ip-1)/my)+1
          jc=mod((ip-1),my)+1
        else
          ic=mod((ip-1),mx)+1
          jc=((ip-1)/mx)+1
        end if

        xc(ip) = (ic-1) * xlngth / (2*nx-2)
        yc(ip) = (jc-1) * ylngth / (2*ny-2)

      end do

      if ( 2 .le. iwrite ) then

        write(*,*)' '
        write(*,*)'    I      XC           YC'
        write(*,*)' '

        do i=1,np
          write (*,'(i5,2f12.5)')i,xc(i),yc(i)
        end do

      end if

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
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
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
      function ubdry ( y, para )

c*********************************************************************72
c
cc UBDRY sets the parabolic inflow in terms of the value of the parameter.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    23 January 2007
c
c  Author:
c
c    John Burkardt
c
      implicit none

      double precision para
      double precision ubdry
      double precision y

      ubdry = 4.0D+00 * para * y * ( 3.0D+00 - y ) / 9.0D+00

      return
      end
      subroutine uval (g,indx,iquad,it,nelemn,neqn,nnodes,node,
     &  np,nquad,para,phi,un,unx,uny,yc)

c*********************************************************************72
c
cc UVAL evaluates the velocities in a particular triangle.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    23 January 2007
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer nelemn
      integer neqn
      integer nnodes
      integer np
      integer nquad

      double precision bb
      double precision bx
      double precision by
      double precision g(neqn)
      integer indx(np,2)
      integer ip
      integer iq
      integer iquad
      integer it
      integer iuk
      integer iun
      integer kk
      integer node(nelemn,nnodes)
      double precision para
      double precision phi(nelemn,nquad,nnodes,3)
      double precision ubc
      double precision ubdry
      double precision un(2)
      double precision unx(2)
      double precision uny(2)
      double precision yc(np)

      do kk = 1, 2
        un(kk) = 0.0D+00
        uny(kk) = 0.0D+00
        unx(kk) = 0.0D+00
      end do

      do iq = 1, nnodes

        ip=node(it,iq)

        bb=phi(it,iquad,iq,1)
        bx=phi(it,iquad,iq,2)
        by=phi(it,iquad,iq,3)

        do iuk = 1, 2

          iun = indx(ip,iuk)

          if (iun.gt.0) then

            un(iuk) = un(iuk)+bb*g(iun)
            unx(iuk) = unx(iuk)+bx*g(iun)
            uny(iuk) = uny(iuk)+by*g(iun)

          else if (iun.lt.0) then

            ip = node(it,iq)
            ubc = ubdry(yc(ip),para)
            un(iuk) = un(iuk)+bb*ubc
            unx(iuk) = unx(iuk)+bx*ubc
            uny(iuk) = uny(iuk)+by*ubc

          end if

        end do

      end do

      return
      end
      subroutine uv_plot3d (f,indx,insc,ivunit,long,mx,my, 
     &    nelemn,neqn,nnodes,node,np,para,press,reynld,yc)

c*********************************************************************72
c
cc UV_PLOT3D creates a velocity file for use by PLOT3D.
c
c  Discussion:
c
c    Given the following set of nodes:
c
c    A  B  C
c    D  E  F
c    G  H  I
c
c    the file will have the form:
c
c    D, U(G), V(G), P
c    D, U(H), V(H), P
c    D, U(I), V(I), P
c    D, U(D), V(D), P
c    D, U(E), V(E), P
c    D, U(F), V(F), P
c    D, U(A), V(A), P
c    D, U(B), V(B), P
c    D, U(C), V(C), P
c
c    Here both D and P are set to 1 for now, representing dummy values
c    of density and pressure.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    23 January 2007
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, DOUBLE PRECISION PARA, the value of the parameter.
c
c    Workspace, DOUBLE PRECISION PRESS(MX,MY), used to hold the 
c    computed pressures.
c
      implicit none

      integer nelemn
      integer neqn
      integer nnodes
      integer np

      double precision alpha
      double precision dval
      double precision f(neqn)
      double precision fsmach
      integer i
      integer ii
      integer indx(np,2)
      integer insc(np)
      integer ip
      integer iset
      integer ivunit
      integer j
      integer k
      logical long
      integer mx
      integer my
      integer node(nelemn,nnodes)
      double precision para
      double precision press(mx,my)
      double precision reynld
      double precision time
      double precision ubdry
      double precision uval
      double precision vval
      double precision yc(np)

      save iset

      data iset /0/

      iset = iset + 1

      dval = 1.0D+00
      fsmach = 1.0D+00
      alpha = 1.0D+00
      time = 1.0D+00

      call pval (f,insc,long,mx,my,nelemn,neqn,nnodes,node,np,press)
c
c  If NX > NY, then nodes with a constant Y value are numbered consecutively.
c
      if(long)then

        write(ivunit,'(2I5)')mx,my
        write(ivunit,'(4G15.5)')fsmach,alpha,reynld,time

        do ii=1,4
          do j=1,my
            do i=1,mx

              ip=(i-1)*my+j

              if(ii.eq.1)then

                write(ivunit,'(G15.5)')dval

              elseif(ii.eq.2)then

                k = indx(ip,1)
                if (k.eq.0) then
                  uval = 0.0
                elseif (k.lt.0) then
                  uval=ubdry(yc(ip),para)
                else
                  uval = f(k)
                end if
                write(ivunit,'(G15.5)')uval

              elseif(ii.eq.3)then

                k = indx(ip,2)
                if (k.eq.0) then
                  vval = 0.0
                else
                  vval = f(k)
                end if
                write(ivunit,'(G15.5)')vval

              else

                write(ivunit,'(G15.5)')press(i,j)

              end if

            end do
          end do
        end do
c
c  If NX < NY, then nodes with a constant X value are numbered consecutively.
c
      else

        write(ivunit,'(2I5)')mx,my
        write(ivunit,'(4G15.5)')fsmach,alpha,reynld,time

        do ii=1,4
          do i=1,mx
            do j=1,my

              if(ii.eq.1)then

                 write(ivunit,'(G15.5)')dval

              elseif(ii.eq.2)then

                ip=(i-1)*my+j
                k = indx(ip,1)

                if (k.eq.0) then
                  uval = 0.0D+00
                elseif (k.lt.0) then
                  uval=ubdry(yc(i),para)
                else
                  uval = f(k)
                end if

                write(ivunit,'(G15.5)')uval

              elseif(ii.eq.3)then

                k = indx(ip,2)

                if (k.eq.0) then
                  vval = 0.0D+00
                else
                  vval = f(k)
                end if

                write(ivunit,'(G15.5)')vval

              else

                write(ivunit,'(G15.5)')press(i,j)

              end if

            end do
          end do
        end do

      end if

      write (*,*) 'UVDUMP wrote data set ',iset,' to file.'

      return
      end
      subroutine uv_table ( f, indx, ivunit, neqn, np, para, yc )

c*********************************************************************72
c
cc UV_TABLE creates a velocity table file.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    28 February 2006
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, double precision PARA, the value of the parameter.
c
      implicit none

      integer neqn
      integer np

      double precision f(neqn)
      integer i
      integer indx(np,2)
      integer ip
      integer ivunit
      integer k
      double precision para
      double precision ubdry
      double precision uval
      double precision vval
      double precision yc(np)

      do ip = 1, np

        k = indx(ip,1)

        if ( k .eq. 0 ) then
          uval = 0.0D+00
        elseif (k .lt. 0 ) then
          uval = ubdry ( yc(ip), para )
        else
          uval = f(k)
        end if

        k = indx(ip,2)

        if ( k .eq. 0 ) then
          vval = 0.0D+00
        else
          vval = f(k)
        end if

        write ( ivunit, '(2x,g14.6,2x,g14.6)' ) uval, vval

      end do

      return
      end
      subroutine xy_plot3d ( ixunit, long, np, nx, ny, xc, yc )

c*********************************************************************72
c
cc XY_PLOT3D creates a grid file for use by PLOT3D.
c
c  Discussion:
c
c    Given the following set of nodes:
c
c      A  B  C
c      D  E  F
c      G  H  I
c
c    the file will have the form:
c
c      X(G), X(H), X(I), X(D), X(E), X(F), X(A), X(B), X(C),
c      Y(G), Y(H), Y(I), Y(D), Y(E), Y(F), Y(A), Y(B), Y(C).
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    23 January 2007
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer np

      integer i
      integer ip
      integer iset
      integer ixunit
      integer j
      logical long
      integer mx
      integer my
      integer nx
      integer ny
      double precision xc(np)
      double precision yc(np)

      save iset

      data iset / 0 /

      iset = iset+1

      mx = 2*nx-1
      my = 2*ny-1
c
c  If NX > NY, then nodes with a constant Y value are numbered consecutively.
c
      if (long) then

        write(ixunit,'(2I15)')mx,my

        do i=1,my
          do j=1,mx
            ip=(j-1)*my+i
            write(ixunit,'(G15.5)')xc(ip)
          end do
        end do

        do i=1,my
          do  j=1,mx
            ip=(j-1)*my+i
            write(ixunit,'(G15.5)')yc(ip)
          end do
        end do
c
c  If NX < NY, then nodes with a constant X value are numbered consecutively.
c
      else

        write(ixunit,'(2I15)')my,mx

        do j=1,mx
          do i=1,my
            ip=(j-1)*my+i
            write(ixunit,'(G15.5)')xc(ip)
          end do
        end do

        do j = 1, mx
          do i = 1, my
            ip=(j-1)*my+i
            write(ixunit,'(G15.5)')yc(ip)
          end do
        end do

      end if

      write (*,*) 'XY_PLOT3D wrote data set ',iset,' to file.'

      return
      end
      subroutine xy_table ( ixunit, np, xc, yc )

c*********************************************************************72
c
cc XY_TABLE creates an (X,Y) table file.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    28 February 2006
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer np

      integer ip
      integer ixunit
      double precision xc(np)
      double precision yc(np)

      do ip = 1, np
        write ( ixunit, '(2x,g14.6,2x,g14.6)' ) xc(ip), yc(ip)
      end do

      return
      end
