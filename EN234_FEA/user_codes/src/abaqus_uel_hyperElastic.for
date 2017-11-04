!
!    ABAQUS format UEL subroutine
!
!    This file is compatible with both EN234_FEA and ABAQUS/Standard
!
!    The example implements a standard fully integrated 3D linear elastic continuum element
!
!    The file also contains the following subrouines:
!          abq_UEL_3D_integrationpoints           - defines integration ponits for 3D continuum elements
!          abq_UEL_3D_shapefunctions              - defines shape functions for 3D continuum elements
!          abq_UEL_invert3D                       - computes the inverse and determinant of a 3x3 matrix
!          abq_facenodes_3D                       - returns list of nodes on the face of a 3D element
!
!=========================== ABAQUS format user element subroutine ===================

      SUBROUTINE UEL(RHS,AMATRX,SVARS,ENERGY,NDOFEL,NRHS,NSVARS,
     1     PROPS,NPROPS,COORDS,MCRD,NNODE,U,DU,V,A,JTYPE,TIME,DTIME,
     2     KSTEP,KINC,JELEM,PARAMS,NDLOAD,JDLTYP,ADLMAG,PREDEF,NPREDF,
     3     LFLAGS,MLVARX,DDLMAG,MDLOAD,PNEWDT,JPROPS,NJPROP,PERIOD)
    !
      INCLUDE 'ABA_PARAM.INC'
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
      integer      :: i,j,n_points,kint, nfacenodes, ipoin
      integer      :: face_node_list(8)                       ! List of nodes on an element face
    !
      double precision  ::  xi(3,64)                          ! Volumetric Integration points
      double precision  ::  w(64)                             ! Integration weights
      double precision  ::  N(20)                             ! 3D Shape functions
      double precision  ::  dNdxi(20,3)                       ! 3D Shape function derivatives
      double precision  ::  dxdxi(3,3)                        ! Derivative of position wrt normalized coords
      double precision  ::  dNdx(20,3)                        ! Derivative of shape functions wrt spatial coords
    !
    !   Variables below are for computing integrals over element faces
      double precision  ::  face_coords(3,8)                  ! Coords of nodes on an element face
      double precision  ::  xi2(2,9)                          ! Area integration points
      double precision  ::  N2(9)                             ! 2D shape functions
      double precision  ::  dNdxi2(9,2)                       ! 2D shape function derivatives
      double precision  ::  norm(3)                           ! Normal to an element face
      double precision  ::  dxdxi2(3,2)                       ! Derivative of spatial coord wrt normalized areal coord
    !
      double precision  ::  strain(6)                         ! Strain vector contains [e11, e22, e33, 2e12, 2e13, 2e23]
      double precision  ::  stress(6)                         ! Stress vector contains [s11, s22, s33, s12, s13, s23]
      double precision  ::  D(6,6)                            ! stress = D*(strain)  (NOTE FACTOR OF 2 in shear strain)
      double precision,allocatable  ::  Bstar(:,:)                !strain = B*(dof_total)
      double precision  ::  dxidx(3,3), determinant           ! Jacobian inverse and determinant
      double precision  ::  E, xnu, D44, D11, D12             ! Material properties
      double precision  ::  F(3,3),Finv(3,3),Jac,F2(3,3)      !Deformation gradient tensor
      double precision :: B(100,100)
      double precision :: utemp(NNODE),dNdxTemp(NNODE)
      double precision :: PK2(3,3),qVec(9),PK2FT(3,3),cauchyStress(3,3)
      double precision :: H(6,9),YS(NNODE,3),YP(NNODE,3),YK(NNODE,NNODE)
      double precision :: Y(3*NNODE,3*NNODE)
      integer :: ie

    !
    !     Example ABAQUS UEL implementing 3D linear elastic elements
    !     El props are:

    !     PROPS(1)         Young's modulus
    !     PROPS(2)         Poisson's ratio

      if (NNODE == 4) n_points = 1               ! Linear tet
      if (NNODE == 10) n_points = 4              ! Quadratic tet
      if (NNODE == 8) n_points = 8               ! Linear Hex
      if (NNODE == 20) n_points = 27             ! Quadratic hex

      call abq_UEL_3D_integrationpoints(n_points, NNODE, xi, w)
      
      
      if (MLVARX<3*NNODE) then
        write(6,*) ' Error in abaqus UEL '
        write(6,*) ' Variable MLVARX must exceed 3*NNODE'
        write(6,*) ' MLVARX = ',MLVARX,' NNODE = ',NNODE
        stop
      endif
      
      !------------
      !Allocate necessary matrices
      !-----------
      allocate(Bstar(9,3*NNODE))
      Bstar= 0.d0
      

      RHS(1:MLVARX,1) = 0.d0
      AMATRX(1:NDOFEL,1:NDOFEL) = 0.d0
      
      

      ENERGY(1:8) = 0.d0

    !     --  Loop over integration points
      do kint = 1, n_points
      
        F=0.d0
        call abq_UEL_3D_shapefunctions(xi(1:3,kint),NNODE,N,dNdxi)
        dxdxi = matmul(coords(1:3,1:NNODE),dNdxi(1:NNODE,1:3))
        call abq_UEL_invert3d(dxdxi,dxidx,determinant)
        
        !------------------
        !Assemble the B Matrix
        !------------------
        
        dNdx(1:NNODE,1:3) = matmul(dNdxi(1:NNODE,1:3),dxidx)
        
        Bstar = 0.d0
        Bstar(1,1:3*NNODE-2:3) = dNdx(1:NNODE,1)
        Bstar(2,2:3*NNODE-1:3) = dNdx(1:NNODE,2)
        Bstar(3,3:3*NNODE:3)   = dNdx(1:NNODE,3)
        Bstar(4,1:3*NNODE-2:3) = dNdx(1:NNODE,2)        
        Bstar(5,2:3*NNODE-1:3) = dNdx(1:NNODE,1)        
        Bstar(6,1:3*NNODE-2:3) = dNdx(1:NNODE,3)
        Bstar(7,3:3*NNODE:3)   = dNdx(1:NNODE,1)
        Bstar(8,2:3*NNODE-1:3) = dNdx(1:NNODE,3)
        Bstar(9,3:3*NNODE:3)   = dNdx(1:NNODE,2)
        
        
        !---------
        !Calculate the F matrix
        !---------
        
        do i = 1,3
            ie = 3*(NNODE-1)+i
            F(i,1:3) = matmul(U(i:ie:3),dNdx(1:NNODE,1:3))
            F(i,i) = F(i,i) + 1.d0
        end do
        
        !------------------------
        !Calculate the stress from the material properties
        !------------------------
        call abq_UEL_invert3d(F,Finv,Jac)
        call hyper_Fung(PROPS(1:NPROPS),NPROPS,F,Jac,PK2,D)
        
        PK2FT = matmul(PK2,transpose(F))
        qVec =(/PK2FT(1,1),PK2FT(2,2),PK2FT(3,3),PK2FT(2,1),PK2FT(1,2),
     +          PK2FT(3,1),PK2FT(1,3),PK2FT(3,2),PK2FT(2,3)/)
        

        H = 0.d0
        H(1,1:7) = [F(1,1),0.d0,0.d0,0.d0,F(2,1),0.d0,F(3,1)]
        H(2,2) = F(2,2)
        H(2,4) = F(1,2)
        H(2,9) = F(3,2)
        H(3,3) = F(3,3)
        H(3,6) = F(1,3)
        H(3,8) = F(2,3)
        H(4,1:9) = [F(1,2),F(2,1),0.d0,F(1,1),F(2,2),0.d0,F(3,2),
     +             0.d0,F(3,1)]
        H(5,1:9) = [F(1,3),0.d0,F(3,1),0.d0,F(2,3),F(1,1),F(3,3),
     +             F(2,1),0.d0]
        H(6,1:9) = [0.d0,F(2,3),F(3,2),F(1,3),0.d0,F(1,2),0.d0,
     +             F(2,2),F(3,3)] 
      
        

      !-----------------
      !Calculating the Y matrix
      !-----------------
     
      
     
        !YS=0.d0
        !YK = 0.d0
        !YS = matmul(dNdx(1:NNODE,1:3),PK2)
        YK = matmul(dNdx(1:NNODE,1:3),
     +       matmul(PK2,transpose(dNdx(1:NNODE,1:3))))      
        Y=0.d0
        do i=1,NNODE
            do j=1,NNODE
                Y(i*3-2,j*3-2) = YK(i,j)
                Y(i*3-1,j*3-1) = YK(i,j)
                Y(i*3,j*3) = YK(i,j)
           end do
        end do
         
       !write(unit=222,FMT=*)"Y"
       !write(unit=222,FMT='(F8.5,",",f8.5)')((Y(i1,k1),k1=1,24),i1=1,24)
        
        cauchyStress = matmul(F,PK2FT)
        cauchyStress = cauchyStress/Jac
        stress(1) = cauchyStress(1,1)
        stress(2) = cauchyStress(2,2)
        stress(3) = cauchyStress(3,3)
        stress(4) = cauchyStress(1,2)
        stress(5) = cauchyStress(1,3)
        stress(6) = cauchyStress(2,3)
                
        
        RHS(1:3*NNODE,1) = RHS(1:3*NNODE,1)
     1   - matmul(transpose(Bstar),qVec)*w(kint)*determinant
                                              

        AMATRX(1:3*NNODE,1:3*NNODE) = AMATRX(1:3*NNODE,1:3*NNODE)
     1  + matmul(transpose(Bstar),(matmul(transpose(H),
     2    (matmul(D,matmul(H,Bstar))))))*w(kint)*determinant
     
        AMATRX(1:3*NNODE,1:3*NNODE) = AMATRX(1:3*NNODE,1:3*NNODE)
     +     +  Y*w(kint)*determinant


   !   ENERGY(2) = ENERGY(2)
   !1   + 0.5D0*dot_product(stress,strain)*w(kint)*determinant           ! Store the elastic strain energy

        if (NSVARS>=n_points*6) then   ! Store stress at each integration point (if space was allocated to do so)
            SVARS(6*kint-5:6*kint) = stress(1:6)
        endif
      end do
      
      
      PNEWDT = 1.d0 ! This leaves the timestep unchanged (ABAQUS will use its own algorithm to determine DTIME)
    
    
    !----------------------------------------------------------------------------------------------------------
    !Apply distributed loads    
    !Distributed loads are specified in the input file using the Un option in the input file.
    !n specifies the face number, following the ABAQUS convention
    !----------------------------------------------------------------------------------------------------------
    
      do j = 1,NDLOAD
      
        call abq_facenodes_3D(NNODE,iabs(JDLTYP(j,1)),
     1                                     face_node_list,nfacenodes)

        do i = 1,nfacenodes
            face_coords(1:3,i) = coords(1:3,face_node_list(i))
        end do

        if (nfacenodes == 3) n_points = 3
        if (nfacenodes == 6) n_points = 4
        if (nfacenodes == 4) n_points = 4
        if (nfacenodes == 8) n_points = 9

        call abq_UEL_2D_integrationpoints(n_points, nfacenodes, xi2, w)

        do kint = 1,n_points
            call abq_UEL_2D_shapefunctions(xi2(1:2,kint),
     1                        nfacenodes,N2,dNdxi2)
            dxdxi2 = matmul(face_coords(1:3,1:nfacenodes),
     1                           dNdxi2(1:nfacenodes,1:2))
            norm(1)=(dxdxi2(2,1)*dxdxi2(3,2))-(dxdxi2(2,2)*dxdxi2(3,1))
            norm(2)=(dxdxi2(1,1)*dxdxi2(3,2))-(dxdxi2(1,2)*dxdxi2(3,1))
            norm(3)=(dxdxi2(1,1)*dxdxi2(2,2))-(dxdxi2(1,2)*dxdxi2(2,1))

            do i = 1,nfacenodes
                ipoin = 3*face_node_list(i)-2
                RHS(ipoin:ipoin+2,1) = RHS(ipoin:ipoin+2,1)
     1                 - N2(1:nfacenodes)*adlmag(j,1)*norm(1:3)*w(kint)      ! Note determinant is already in normal
            end do
        end do
      end do

      deallocate(Bstar)
      return

      END SUBROUTINE UEL
      

      
      
      
      !-------------------------------------------
      !Subroutine to calculate D,PK2 from F and material constants
      !-------------------------------------------      
      
      subroutine hyper_Fung(elem_props,n_props,F,Jac,PK2,D)

           implicit none
     
           integer, intent(in)           :: n_props
           double precision, intent(in)  :: elem_props(n_props)
           double precision, intent(in)  :: F(3,3)
           double precision, intent(in)  :: Jac
           double precision, intent(out) :: PK2(3,3)
           double precision, intent(out) :: D(6,6)
           
     
           !double precision :: B(3,3)           
           double precision :: C(3,3),Cbar(3,3),Cbar_I(3,3)
           double precision :: Cinv(3,3),Cbar_I_vec(6)
           double precision :: Cvec(6),CstarVec(6),CbarStarVec(6)
           double precision :: Binvvec(6),CinvVec(6),PVec(6),newPvec(6)
           double precision :: eyevec(6),PK2Vec(6)
           
           double precision :: mu
           double precision :: K
           double precision :: ss
           double precision :: trB
           double precision :: Gvec(6),G(6,6),Omega(6,6),GCStarVec(6)
           double precision :: Q,P
           integer :: indexVec(6,2)
     
           integer :: i,j,temp1,temp2
           
           indexVec = RESHAPE((/1,2,3,1,1,2,1,2,3,2,3,3/),(/6,2/))
           mu = elem_props(1)
           K  = elem_props(2)
           G = 0.d0
           Gvec(1) = elem_props(3)
           Gvec(2) = elem_props(4)
           Gvec(3) = elem_props(5)
           Gvec(4) = elem_props(6)
           Gvec(5) = elem_props(6)
           Gvec(6) = elem_props(6)
           
           G(1,1) = elem_props(3)
           G(2,2) = elem_props(4)
           G(3,3) = elem_props(5)
           G(4,4) = elem_props(6)
           G(5,5) = elem_props(6)
           G(6,6) = elem_props(6)
           
           eyevec(1:3) = 1.d0
           eyevec(4:6) = 0.d0
           
           
           C = matmul(transpose(F),F)           
           call abq_UEL_invert3d(C,Cinv,ss)      ! ss is just a dummy variable here
           Cbar = C
           Cbar = Cbar*(Jac**(-2.d0/3.d0))
           Cbar_I = Cbar
           
           do i=1,3
              Cbar_I(i,i)= Cbar_I(i,i)-1.0
           end do
           Cbar_I_vec = (/Cbar_I(1,1),Cbar_I(2,2),Cbar_I(3,3),
     +                   Cbar_I(1,2),Cbar_I(1,3),Cbar_I(2,3)/)
           
           Q = .25d0*dot_product(Cbar_I_vec,(matmul(G,Cbar_I_vec)))
           
           
           
           !---------------
           !Calculate the stress
           !---------------
!           
!           do i=1,6
!           
!              temp1 = indexVec(i,1)
!              temp2 = indexVec(i,2)
!              
!              P=0.d0
!              do j=1,6
!                  P=P+C(indexVec(j,1),indexVec(j,2))
!     +              *G(j,j)*Cbar_I(indexVec(j,1),indexVec(j,2))
!              end do
!              
!              P=P*Cinv(temp1,temp2)
!              P=G(i,i)*Cbar_I(temp1,temp2)-(1.d0/3.d0)*P
!              P=P/(Jac**(-2.d0/3.d0))
!              PVec(i) = P
!              
!              PK2(temp1,temp2)=
!     +                 mu*exp(Q)*P+K*(Jac-1.d0)*Cinv(temp1,temp2)
!                  
!              PK2(temp1,temp2)=PK2(temp2,temp1)
!            end do
            
            
            !---------------
            !Calculating the PK2 stress
            !---------------
            
            Cvec=(/C(1,1),C(2,2),C(3,3),C(1,2),C(1,3),C(2,3)/)
            CstarVec=(/C(1,1),C(2,2),C(3,3),2*C(1,2),2*C(1,3),2*C(2,3)/)
            CbarStarVec= CstarVec*(Jac**(-2.d0/3.d0))
            CinvVec = (/Cinv(1,1),Cinv(2,2),Cinv(3,3),
     +                 Cinv(1,2),Cinv(1,3),Cinv(2,3)/)
     
            newPvec= 0.5d0*(Jac**(-2.d0/3.d0))*
     +              ((Gvec*(CbarStarVec-eyeVec))-
     +           (1.d0/3d0)*(dot_product(CstarVec,
     +           Gvec*(CbarStarVec-eyeVec)))*CinvVec)
     
               
     
            PK2Vec = mu*exp(Q)*newPvec+K*Jac*(Jac-1.d0)*CinvVec
            PK2(1,1) = PK2Vec(1)
            PK2(2,2) = PK2Vec(2)
            PK2(3,3) = PK2Vec(3)
            PK2(1,2) = PK2Vec(4)
            PK2(2,1) = PK2Vec(4)
            PK2(1,3) = PK2Vec(5)
            PK2(3,1) = PK2Vec(5)
            PK2(2,3) = PK2Vec(6)
            PK2(3,2) = PK2Vec(6)
     
            
     
           
            !------
            !Compute Omega
            !-------
     
            Omega = 0.d0
            do i=1,3
              do j=1,3
                  if (j<i) then
                      cycle
                  end if
                  Omega(i,j) = Cinv(i,j)*Cinv(i,j)
                  Omega(j,i) = Omega(i,j)
              end do
            end do
            
            Omega(1,4:6) = [Cinv(1,1)*Cinv(1,2),Cinv(1,1)*Cinv(1,3),
     +                      Cinv(1,2)*Cinv(1,3)]
            Omega(2,4:6) = [Cinv(2,1)*Cinv(2,2),Cinv(2,1)*Cinv(2,3),
     +                      Cinv(2,2)*Cinv(2,3)]
            Omega(3,4:6) = [Cinv(3,1)*Cinv(3,2),Cinv(3,1)*Cinv(3,3),
     +                      Cinv(3,2)*Cinv(3,3)]
     
            Omega(4,4:6)=[0.5*(Cinv(1,1)*Cinv(2,2)+Cinv(1,2)*Cinv(1,2)),
     +                    0.5*(Cinv(1,1)*Cinv(2,3)+Cinv(1,3)*Cinv(1,2)),
     +                    0.5*(Cinv(1,2)*Cinv(2,3)+Cinv(1,3)*Cinv(2,2))]
     
            Omega(5,5:6)=[0.5*(Cinv(1,1)*Cinv(3,3)+Cinv(1,3)*Cinv(1,3)),
     +                    0.5*(Cinv(1,2)*Cinv(3,3)+Cinv(1,3)*Cinv(2,3))]
     
            Omega(6,6)=  0.5*(Cinv(2,2)*Cinv(3,3)+Cinv(2,3)*Cinv(2,3))
            
            do i=1,6
              do j=1,6
                  if (j>=i) then
                      cycle
                  end if
                  Omega(i,j) = Omega(j,i)
               end do
            end do
            
          !write(unit=222,FMT='(F15.5,",",f15.5)')((G(i,k),k=1,6),i=1,6)
          !write(6,*) shape(G),shape(CbarStarVec),CbarStarVec-eyeVec
          !write(6,*) "Shape of G",Gvec,Gvec*(CbarStarVec-eyeVec)
            

            
        !--------------------
        !D Matrix calculation
        !--------------------
         D=0.d0
         GCstarVec= Gvec*CStarVec
         
         D = mu*exp(Q)*(Jac**(-4.d0/3.d0))*(G-(1.d0/3.d0)*(
     +       spread(GCstarVec,dim=2,ncopies=6)*
     +       spread(CinvVec,dim=1,ncopies=6)+
     +       spread(CinvVec,dim=2,ncopies=6)*
     +       spread(GCStarVec,dim=1,ncopies=6)))
    
        D = D- mu*exp(Q)*(Jac**(-2.d0/3.d0))*(1.d0/3.d0)*(
     + dot_product(CStarVec,Gvec*(CbarStarVec-eyeVec))*Omega)
     
        D = D+  mu*exp(Q)*(Jac**(-4.d0/3.d0))*(1.d0/9.d0)*(
     +      dot_product(CStarVec,GCStarVec)*  
     +spread(CinvVec,dim=2,ncopies=6)*spread(CinvVec,dim=1,ncopies=6))
     
        D = D+ mu*exp(Q)*spread(2*newPVec,dim=2,ncopies=6)*
     +spread(newPVec-(1.d0/3.d0)*CinvVec,dim=1,ncopies=6)
     
        D = D-mu*exp(Q)*(1.d0/3.d0)*(Jac**(-2.d0/3.d0))*(
     + spread(CinvVec,dim=2,ncopies=6)*spread(Gvec*(CbarStarVec-eyeVec),
     + dim=1,ncopies=6))
     
        D= D+K*Jac*(2.d0*Jac-1.d0)*(
     +spread(CinvVec,dim=2,ncopies=6)*spread(CinvVec,dim=1,ncopies=6))
     
        D=D+ 2.d0*K*Jac*(Jac-1.d0)*omega
     
        
        !write(unit=222,FMT='(F15.5,",",f15.5)')((D(i,k),k=1,6),i=1,6)
!     
       return
     
      end subroutine hyper_Fung