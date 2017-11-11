!
!    ABAQUS format UEL subroutine
!
!    This file is compatible with both EN234_FEA and ABAQUS/Standard
!
!    The example implements a standard fully integrated 3D linear elastic continuum element
!
!    The file also needs the following subrouines:
!          abq_UEL_2D_integrationpoints           - defines integration points for 2D continuum elements
!          abq_UEL_2D_shapefunctions              - defines shape functions for 2D continuum elements
!=========================== ABAQUS format user element subroutine ===================

      SUBROUTINE UEL(RHS,AMATRX,SVARS,ENERGY,NDOFEL,NRHS,NSVARS,
     1     PROPS,NPROPS,COORDS,MCRD,NNODE,U,DU,V,A,JTYPE,TIME,DTIME,
     2     KSTEP,KINC,JELEM,PARAMS,NDLOAD,JDLTYP,ADLMAG,PREDEF,NPREDF,
     3     LFLAGS,MLVARX,DDLMAG,MDLOAD,PNEWDT,JPROPS,NJPROP,PERIOD)
    !
      INCLUDE 'ABA_PARAM.INC'
      !implicit none
    !
    !
      DIMENSION RHS(MLVARX,*),AMATRX(NDOFEL,NDOFEL),PROPS(*),
     1   SVARS(*),ENERGY(8),COORDS(MCRD,NNODE),U(NDOFEL),
     2   DU(MLVARX,*),V(NDOFEL),A(NDOFEL),TIME(2),PARAMS(*),
     3   JDLTYP(MDLOAD,*),ADLMAG(MDLOAD,*),DDLMAG(MDLOAD,*),
     4   PREDEF(2,NPREDF,NNODE),LFLAGS(*),JPROPS(*)

    !
    !       Variables that must be computed in this routine
    !       RHS(i)                     Right hand side vector.  In EN234_FEA the dimensions are always RHS(MLVARX,1)
    !       AMATRX(i,j)                Stiffness matrix d RHS(i)/ d DU(j)
    !       SVARS(1:NSVARS)            Element state variables.  Must be updated in this routine
    !       ENERGY(1:8)
    !                                  Energy(1) Kinetic Energy
    !                                  Energy(2) Elastic Strain Energy
    !                                  Energy(3) Creep Dissipation
    !                                  Energy(4) Plastic Dissipation
    !                                  Energy(5) Viscous Dissipation
    !                                  Energy(6) Artificial strain energy
    !                                  Energy(7) Electrostatic energy
    !                                  Energy(8) Incremental work done by loads applied to the element
    !       PNEWDT                     Allows user to control ABAQUS time increments.
    !                                  If PNEWDT<1 then time step is abandoned and computation is restarted with
    !                                  a time increment equal to PNEWDT*DTIME
    !                                  If PNEWDT>1 ABAQUS may increase the time increment by a factor PNEWDT
    !
    !       Variables provided for information
    !       NDOFEL                     Total # DOF for the element
    !       NRHS                       Dimension variable
    !       NSVARS                     Total # element state variables
    !       PROPS(1:NPROPS)            User-specified properties of the element
    !       NPROPS                     No. properties
    !       JPROPS(1:NJPROPS)          Integer valued user specified properties for the element
    !       NJPROPS                    No. integer valued properties
    !       COORDS(i,N)                ith coordinate of Nth node on element
    !       MCRD                       Maximum of (# coords,minimum of (3,#DOF)) on any node
    !       U                          Vector of DOF at the end of the increment
    !       DU                         Vector of DOF increments
    !       V                          Vector of velocities (defined only for implicit dynamics)
    !       A                          Vector of accelerations (defined only for implicit dynamics)
    !       JTYPE                      Integer identifying element type (the number n in the Un specification in the input file)
    !       TIME(1:2)                  TIME(1)   Current value of step time
    !                                  TIME(2)   Total time
    !       DTIME                      Time increment
    !       KSTEP                      Current step number (always 1 in EN234_FEA)
    !       KINC                       Increment number
    !       JELEM                      User assigned element number in ABAQUS (internally assigned in EN234_FEA)
    !       PARAMS(1:3)                Time increment parameters alpha, beta, gamma for implicit dynamics
    !       NDLOAD                     Number of user-defined distributed loads defined for this element
    !       JDLTYP(1:NDLOAD)           Integers n defining distributed load types defined as Un or (if negative) UnNU in input file
    !       ADLMAG(1:NDLOAD)           Distributed load magnitudes
    !       DDLMAG(1:NDLOAD)           Increment in distributed load magnitudes
    !       PREDEF(1:2,1:NPREDF,1:NNODE)   Predefined fields.
    !       PREDEF(1,...)              Value of predefined field
    !       PREDEF(2,...)              Increment in predefined field
    !       PREDEF(1:2,1,k)            Value of temperature/temperature increment at kth node
    !       PREDEF(1:2,2:NPREDF,k)     Value of user defined field/field increment at kth node (not used in EN234FEA)
    !       NPREDF                     Number of predefined fields (1 for en234FEA)
    !       LFLAGS                     Control variable
    !       LFLAGS(1)                  Defines procedure type
    !       LFLAGS(2)                  0 => small displacement analysis  1 => Large displacement (NLGEOM option)
    !       LFLAGS(3)                   1 => Subroutine must return both RHS and AMATRX (always true in EN234FEA)
    !                                   2 => Subroutine must return stiffness AMATRX = -dF/du
    !                                   3 => Subroutine must return daming matrix AMATRX = -dF/dudot
    !                                   4 => Subroutine must return mass matrix AMATRX = -dF/duddot
    !                                   5 => Define the RHS only
    !                                   6 => Define the mass matrix for the initial acceleration calculation
    !                                   100 => Define perturbation quantities for output
    !       LFLAGS(4)                   0 => General step   1 => linear perturbation step
    !       LFLAGS(5)                   0 => current approximation to solution based on Newton correction; 1 => based on extrapolation
    !       MLVARX                      Dimension variable (equal to NDOFEL in EN234FEA)
    !       PERIOD                      Time period of the current step
    !
    !
    ! Local Variables
      integer      :: i,j,n_points,kint
    !
      double precision  ::  xi(2,9)                          ! Area integration points
      double precision  ::  w(9)                             ! Area integration weights
      double precision  ::  N(9),Nstar(9)                    ! 2D shape functions
      double precision  ::  dNdxi(9,2),dNStardxi(9,2)        ! 2D shape function derivatives
      double precision  ::  dNdx(9,2) ,dNStardx(9,2)         ! Spatial derivatives
      double precision  ::  Nbar(9)
      double precision  ::  dNbardxi(9,2)
      double precision  ::  dNbardx(9,2)
      double precision  ::  dxdxi(2,2),lindxdxi(2,2)          ! Derivative of spatial coords wrt normalized coords

      double precision  ::  sol(9), dsol(9)                   ! Sol vector contains [mu, c, dmudx1, dmudx2, dcdx1, dcdx2]
      double precision  ::  q(9)                              ! q vector defined in class
      double precision  ::  D(9,9),Delastic(4,4),Dtemp(3,3)   ! D matrix defined in class
      double precision  ::  stress(4),strain(4),newStrain(4)  !Stress Strain measures
      double precision  ::  d11,d22,d44                       ! Dummy Dmatrix values
      double precision  ::  E,xnu                              !Elastic properties
      double precision  ::  B(9,24)                           ! p = B*U
      double precision  ::  dxidx(2,2),lindxidx(2,2),determinant,linDet    ! Jacobian inverse and determinant
      double precision  ::  diffusion_coeft,kappa,theta       ! Material properties
      double precision  ::  Omega,WEnergCost                  ! properties
      double precision  ::  c                                 ! concentration
      integer  ::k1,l1                                        ! Dummy integer
      double precision :: tempD0                              !Dummy float
      

    !
    !     Example ABAQUS UEL implementing 2D phase field model
      
      if (NNODE/=8) then
        write(6,*) ' The UEL should be used with only8 noded elements'
        write(6,*) ' The analysis will now stop'  
        stop
      else
          n_points = 4
      endif


      call abq_UEL_2D_integrationpoints(n_points, NNODE, xi, w)

      RHS(1:MLVARX,1) = 0.d0
      AMATRX(1:NDOFEL,1:NDOFEL) = 0.d0

      E = PROPS(1)
      xnu = PROPS(2)
      Omega = PROPS(3)
      WEnergCost = PROPS(4)
      kappa = PROPS(5)
      diffusion_coeft = PROPS(6)
      theta = PROPS(7)
      
      !Properties
      !  100.d0, 0.3d0      % E, nu
      !   0.0d0, 1.d0, 0.001d0, 1.d0   % Omega, W, kappa, diffusion coefft
      !   0.5d0                     % Theta

      
      !---------------------------------
      !Compute the D Matrix
      !--------------------------------
      
      
      Dtemp = 0.d0 ! Dtemp is a reduced 3x3 matrix corresponding to s11,s22,s12 = Dtemp * (e11,e22,e12)
      
      d44 = 0.5D0*E/(1+xnu)
      d11 = (1.D0-xnu)*E/( (1+xnu)*(1-2.D0*xnu) )      
      d12 = xnu*E/( (1+xnu)*(1-2.D0*xnu) )
      Dtemp(1:2,1:2)=d11
      Dtemp(1,2)=d12
      Dtemp(2,1)=d12      
      Dtemp(3,3) = d44 
      
      !--------------
      !Define D for a plane strain matrix
      !--------------
      Delastic = 0.d0 !The regular D Matrix for plane strain problem
      Delastic(1:3,1:3) = d12
      Delastic(1,1) = d11
      Delastic(2,2) = d11
      Delastic(3,3) = d11
      Delastic(4,4) = d44
      
      D=0.d0
      D(1:3,1:3) = Dtemp(1:3,1:3)
      D(4,1)=(-Omega/3.d0)*(Delastic(1,1)+Delastic(1,2)+Delastic(1,3))
      D(4,2)=(-Omega/3.d0)*(Delastic(2,1)+Delastic(2,2)+Delastic(2,3))
      D(4,4) = 1.d0
      D(1,5)=(-Omega/3.d0)*(Delastic(1,1)+Delastic(1,2)+Delastic(1,3))
      D(2,5)=(-Omega/3.d0)*(Delastic(2,1)+Delastic(2,2)+Delastic(2,3))
      D(8,6) = theta*diffusion_coeft
      D(9,7) = theta*diffusion_coeft
      D(7,9)=-kappa
      D(6,8)=-kappa
      D(5,5) = 1.d0/DTIME
      
      
      tempD0 = 0.d0
      
      do k1=1,3
          do l1=1,3
              tempD0=tempD0+Delastic(k1,l1)
          end do
      end do
      
      tempD0 = Omega*Omega*tempD0/9.d0
      
      
    !     --  Loop over integration points
      do kint = 1, n_points
          
        !---------------
        !Calculate regular shape function
        !--------------
        call abq_UEL_2D_shapefunctions(xi(1:2,kint),NNODE,N,dNdxi)
        dxdxi = matmul(coords(1:2,1:NNODE),dNdxi(1:NNODE,1:2))        
        determinant = dxdxi(1,1)*dxdxi(2,2)-dxdxi(1,2)*dxdxi(2,1)        
        dxidx(1,1:2) =  [dxdxi(2,2),-dxdxi(1,2)]/determinant
        dxidx(2,1:2) =  [-dxdxi(2,1),dxdxi(1,1) ]/determinant        
        dNdx(1:NNODE,1:2) = matmul(dNdxi(1:NNODE,1:2),dxidx)
        
        
        !---------------
        !Calculate shape function for 4 noded linear elements
        !--------------
        call abq_UEL_2D_shapefunctions(xi(1:2,kint),4,Nstar,
     +         dNstardxi)
         lindxdxi= matmul(coords(1:2,1:4),dNstardxi(1:4,1:2))        
         linDet=lindxdxi(1,1)*lindxdxi(2,2)-lindxdxi(1,2)*lindxdxi(2,1)        
         lindxidx(1,1:2) =  [ lindxdxi(2,2),-lindxdxi(1,2)]/linDet
         lindxidx(2,1:2) =  [-lindxdxi(2,1),lindxdxi(1,1) ]/linDet  
         dNStardx(1:4,1:2) = matmul(dNstardxi(1:4,1:2),lindxidx)
         
         
         
         
         
        !------------------
        !Compute the B Matrix
        !-----------------
        B = 0.d0
        
        B(1,1:2*NNODE:4) = dNdx(1:NNODE/2,1)
        B(1,2*NNODE+1:3*NNODE:2) = dNdx(1+NNODE/2:NNODE,1)
        
        B(2,2:2*NNODE:4) = dNdx(1:NNODE/2,2)
        B(2,2*NNODE+2:3*NNODE:2) = dNdx(1+NNODE/2:NNODE,2)
        
        B(3,1:2*NNODE:4) = dNdx(1:NNODE/2,2)
        B(3,2*NNODE+1:3*NNODE:2) = dNdx(1+NNODE/2:NNODE,2)
        
        B(3,2:2*NNODE:4) = dNdx(1:NNODE/2,1)
        B(3,2*NNODE+2:3*NNODE:2) = dNdx(1+NNODE/2:NNODE,1)
        
        B(4,3:2*NNODE:4) = Nstar(1:4)
        B(5,4:2*NNODE:4) = Nstar(1:4)
        
        B(6,3:2*NNODE:4) = dNStardx(1:4,1)
        B(7,3:2*NNODE:4) = dNStardx(1:4,2)
        
        B(8,4:2*NNODE:4) = dNStardx(1:4,1)
        B(9,4:2*NNODE:4) = dNStardx(1:4,2)
        
        
               
        !------------------
        !Compute the rest of D Matrix and assemble q vector
        !-----------------
        
        sol = matmul(B(1:9,1:3*NNODE),U(1:3*NNODE))      ! The p vector at the end of the step
        dsol = matmul(B(1:9,1:3*NNODE),DU(1:3*NNODE,1))    ! Increment in the p vector
        c = sol(5)
        
        strain= (/sol(1),sol(2),0.d0,sol(3)/)
        newStrain = strain-c*Omega*(/1.d0,1.d0,1.d0,0.d0/)
        stress=matmul(Delastic,newStrain)
        
        q=0.d0
        
        q(1) = stress(1)
        q(2) = stress(2)
        q(3) = stress(4)
       
        q(4) = sol(4)- (2.d0*WEnergCost*c*(c-1.d0)*(2.d0*c-1.d0))
     +          -((Omega/3.d0)*(stress(1)+stress(2)+stress(3)))
        
        q(5) = dsol(5)/DTIME
        q(6) = -kappa*sol(8)
        q(7) = -kappa*sol(9)
        q(8) = diffusion_coeft*(sol(6)+(theta-1.d0)*dsol(6))
        q(9) = diffusion_coeft*(sol(7)+(theta-1.d0)*dsol(7))
                
        D(4,5)=(-WEnergCost*(12.d0*c*c-12.d0*c+2.d0))+tempD0
        
        
        RHS(1:3*NNODE,1) = RHS(1:3*NNODE,1) 
     1   - matmul(transpose(B(1:9,1:3*NNODE)),q)*w(kint)*determinant

        AMATRX(1:3*NNODE,1:3*NNODE) = AMATRX(1:3*NNODE,1:3*NNODE) 
     1   + matmul(transpose(B(1:9,1:3*NNODE)),
     2             matmul(D,B(1:9,1:3*NNODE)))*w(kint)*determinant
     
        if (NSVARS>=n_points*4) then   ! Store stress at each integration point (if space was allocated to do so)            
            SVARS(4*kint-3:4*kint) = stress(1:4)
        end if 
            

      end do
  

      return

      END SUBROUTINE UEL

