module CFDElementM

  use UtilitiesM
  use SparseKit

  use Triangle2D3NodeM
  use Triangle2D6NodeM
  use Quadrilateral2D4NodeM
  use Quadrilateral2D8NodeM

  use IntegratorPtrM

  use LeftHandSideM
  use ProcessInfoM

  use PointM
  use NodeM
  use NodePtrM

  use SourceM
  use SourcePtrM

  use ElementM

  use CFDMaterialM

  implicit none

  private
  public :: CFDElementDT, cfdElement, initGeometries

  type, extends(ElementDT) :: CFDElementDT
     class(CFDMaterialDT), pointer              :: material
     type(NodePtrDT), dimension(:), allocatable :: nodalPoints
     real(rkind)    , dimension(:), allocatable :: dNdx, dNdy
     real(rkind)    , dimension(3)              :: Tau
     real(rkind)                                :: nu
   contains
     procedure, public  :: init
     procedure, public  :: calculateLHS
     procedure, public  :: calculateRHS
     procedure, public  :: calculateLocalSystem
     procedure, public  :: calculateResults
     procedure, public  :: calculateDt
     procedure, public  :: calculateTauNu
     procedure, private :: setupIntegration
     procedure, private :: getValuedSource
     procedure, private :: vectorCalculus
  end type CFDElementDT

  interface cfdElement
     procedure :: constructor
  end interface cfdElement

  type(Triangle2D3NodeDT)     , target, save :: myTriangle2D3Node
  type(Triangle2D6NodeDT)     , target, save :: myTriangle2D6Node
  type(Quadrilateral2D4NodeDT), target, save :: myQuadrilateral2D4Node
  type(Quadrilateral2D8NodeDT), target, save :: myQuadrilateral2D8Node

contains

  type(CFDElementDT) function constructor(id, node, material)
    implicit none
    integer(ikind)                          , intent(in) :: id
    type(NodePtrDT)           , dimension(:), intent(in) :: node
    class(CFDMaterialDT), target            , intent(in) :: material
    call constructor%init(id, node, material)
  end function constructor

  subroutine init(this, id, node, material)
    implicit none
    class(CFDElementDT)                     , intent(inout) :: this
    integer(ikind)                          , intent(in)    :: id
    type(NodePtrDT)           , dimension(:), intent(in)    :: node
    class(CFDMaterialDT), target            , intent(in)    :: material
    this%id = id
    this%node = node
    this%material => material
    if(size(node) == 3) then
       this%geometry => myTriangle2D3Node
    else if(size(node) == 4) then
       this%geometry => myQuadrilateral2D4Node
    else if(size(node) == 6) then
       this%geometry => myTriangle2D6Node
    else if(size(node) == 8) then
       this%geometry => myQuadrilateral2D8Node
    end if
    allocate(this%source(1))
  end subroutine init

  subroutine initGeometries(nGauss)
    implicit none
    integer(ikind), intent(in) :: nGauss
    myTriangle2D3Node = triangle2D3Node(nGauss)
    myTriangle2D6Node = triangle2D6Node(nGauss)
    myQuadrilateral2D4Node = quadrilateral2D4Node(nGauss)
    myQuadrilateral2D8Node = quadrilateral2D8Node(nGauss)
  end subroutine initGeometries

  subroutine calculateLocalSystem(this, processInfo, lhs, rhs)
    implicit none
    class(CFDElementDT)                       , intent(inout) :: this
    type(ProcessInfoDT)                       , intent(inout) :: processInfo
    type(LeftHandSideDT)                      , intent(inout) :: lhs
    real(rkind)    , dimension(:), allocatable, intent(inout) :: rhs
    type(IntegratorPtrDT)                                     :: integrator
    integer(ikind)                                            :: i, j, k, elemID
    integer(ikind)                                            :: nNode, nDof, iNodeID
    real(rkind)    , dimension(:,:,:), allocatable            :: jacobian
    real(rkind)    , dimension(:,:)  , allocatable            :: U, theta
    real(rkind)    , dimension(:)    , allocatable            :: jacobianDet
    real(rkind)                                               :: gamma, V_sq 
    real(rkind)                                               :: bi, ci, v1, v2, e, rho
    real(rkind)                                               :: thetaK(4), AiUi(4), AiUi_theta(4)
    real(rkind)                                               :: U_k(4), theta_k(4), Ux(4), Uy(4)
    real(rkind)                                               :: A1AiUi_theta(4), A2AiUi_theta(4)
    real(rkind)                                               :: K1jUj(4), K2jUj(4), lambda, mu, temp
    real(rkind)                                               :: Cv, T_inf
    nNode = this%getnNode()
    nDof  = 4
    allocate(U(nDof,nNode))
    lhs = leftHandSide(0, 0, nNode*nDof)
    integrator  = this%getIntegrator()
    jacobian    = this%geometry%jacobianAtGPoints(this%nodalPoints)
    jacobianDet = this%geometry%jacobianDetAtGPoints(jacobian)
    
    lambda  = this%material%k
    mu      = this%material%mu
    temp    = 0.d0
    Cv      = this%material%Cv
    T_inf   = this%material%T_inf
    
    gamma   = this%material%gamma
    elemID  = this%getID()
    Ux      = 0.d0
    Uy      = 0.d0
    U_k     = 0.d0
    theta_k = 0.d0
    lhs%stiffness = 0.d0
    do i = 1, nNode
       do j = 1, nDof
          U(j,i) = this%node(i)%ptr%dof(j)%val
       end do
       Ux = Ux + U(:,i)*this%dNdx(i)
       Uy = Uy + U(:,i)*this%dNdy(i)
       temp = temp + (((U(4,i)/U(1,i))-0.5d0*((U(2,i)/U(1,i))**2+(U(3,i)/U(1,i))**2))/Cv)/nNode
    end do

    lambda = lambda*(temp/T_inf)**1.5d0*(T_inf + 194)/(temp + 194)
    
    do k = 1, integrator%getIntegTerms()
       do i = 1, nNode
          iNodeID = this%getNodeID(i)
          U_k     = U_k     + integrator%ptr%shapeFunc(k,i)*U(:,i)  
          theta_k = theta_k + integrator%ptr%shapeFunc(k,i)*processInfo%mat(:,iNodeID)
       end do
       rho  = U_k(1)
       v1   = U_k(2)/rho
       v2   = U_k(3)/rho
       e    = U_k(4)/rho
       V_sq = v1**2 + v2**2
       AiUi(1) = Ux(2) + Uy(3)
       AiUi(2) = (1.0d0/2.0d0)*Ux(1)*(V_sq*(gamma - 1) - 2*v1**2) - Ux(2)*v1*( &
            gamma - 3) - Ux(3)*v2*(gamma - 1) + Ux(4)*(gamma - 1) - Uy(1)*v1*v2 + &
            Uy(2)*v2 + Uy(3)*v1
       AiUi(3) = -Ux(1)*v1*v2 + Ux(2)*v2 + Ux(3)*v1 + (1.0d0/2.0d0)*Uy(1)*(V_sq*( &
            gamma - 1) - 2*v2**2) - Uy(2)*v1*(gamma - 1) - Uy(3)*v2*(gamma - 3) + &
            Uy(4)*(gamma - 1)
       AiUi(4) = Ux(1)*v1*(V_sq*(gamma - 1) - e*gamma) - 1.0d0/2.0d0*Ux(2)*(V_sq &
            *(gamma - 1) - 2*e*gamma + 2*v1**2*(gamma - 1)) - Ux(3)*v1*v2*( &
            gamma - 1) + Ux(4)*gamma*v1 + Uy(1)*v2*(V_sq*(gamma - 1) - e*gamma) - &
            Uy(2)*v1*v2*(gamma - 1) - 1.0d0/2.0d0*Uy(3)*(V_sq*(gamma - 1) - 2*e* &
            gamma + 2*v2**2*(gamma - 1)) + Uy(4)*gamma*v2
       AiUi_theta = theta_k + AiUi
       A1AiUi_theta(1) = AiUi_theta(2)
       A1AiUi_theta(2) = v1*(-gamma + 3)*AiUi_theta(2) - v2*(gamma - 1)* &
            AiUi_theta(3) + (gamma - 1)*AiUi_theta(4) + ((1.0d0/2.0d0)* &
            V_sq*(gamma - 1) - v1**2)*AiUi_theta(1)
       A1AiUi_theta(3) = -v1*v2*AiUi_theta(1) + v1*AiUi_theta(3) + v2* &
            AiUi_theta(2)
       A1AiUi_theta(4) = gamma*v1*AiUi_theta(4) - v1*v2*(gamma - 1)* &
            AiUi_theta(3) + v1*(V_sq*(gamma - 1) - e*gamma)*AiUi_theta(1) &
            + (-1.0d0/2.0d0*V_sq*(gamma - 1) + e*gamma - v1**2*(gamma - 1 &
            ))*AiUi_theta(2)
       A2AiUi_theta(1) = AiUi_theta(3)
       A2AiUi_theta(2) = -v1*v2*AiUi_theta(1) + v1*AiUi_theta(3) + v2* &
            AiUi_theta(2)
       A2AiUi_theta(3) = -v1*(gamma - 1)*AiUi_theta(2) + v2*(-gamma + 3)* &
            AiUi_theta(3) + (gamma - 1)*AiUi_theta(4) + ((1.0d0/2.0d0)* &
            V_sq*(gamma - 1) - v2**2)*AiUi_theta(1)
       A2AiUi_theta(4) = gamma*v2*AiUi_theta(4) - v1*v2*(gamma - 1)* &
            AiUi_theta(2) + v2*(V_sq*(gamma - 1) - e*gamma)*AiUi_theta(1) &
            + (-1.0d0/2.0d0*V_sq*(gamma - 1) + e*gamma - v2**2*(gamma - 1 &
            ))*AiUi_theta(3)
       A1AiUi_theta(1) = A1AiUi_theta(1)*this%Tau(1)
       A1AiUi_theta(2:3) = A1AiUi_theta(2:3)*this%Tau(2)
       A1AiUi_theta(4) = A1AiUi_theta(4)*this%Tau(3)
       A2AiUi_theta(1) = A2AiUi_theta(1)*this%Tau(1)    
       A2AiUi_theta(2:3) = A2AiUi_theta(2:3)*this%Tau(2)
       A2AiUi_theta(4) = A2AiUi_theta(4)*this%Tau(3)

       
       K1jUj(1) = 0.d0
       K1jUj(2) = (2.0d0/3.0d0)*mu*(-2*Ux(1)*v1 + 2*Ux(2) + Uy(1)*v2 - Uy(3))/rho
       K1jUj(3) = mu*(-Ux(1)*v2 + Ux(3) - Uy(1)*v1 + Uy(2))/rho
       K1jUj(4) = (1.0d0/3.0d0)*(Cv*mu*(-Uy(1)*v1*v2 + 3*Uy(2)*v2 - 2*Uy(3)*v1) &
            - Ux(1)*(Cv*mu*(3*V_sq + v1**2) - 3*lambda*(V_sq - e)) + Ux(2)*v1*(4*Cv*mu &
            - 3*lambda) + 3*Ux(3)*v2*(Cv*mu - lambda) + 3*Ux(4)*lambda)/(Cv*rho)
       K2jUj(1) = 0.d0
       K2jUj(2) = mu*(-Ux(1)*v2 + Ux(3) - Uy(1)*v1 + Uy(2))/rho
       K2jUj(3) = (2.0d0/3.0d0)*mu*(Ux(1)*v1 - Ux(2) - 2*Uy(1)*v2 + 2*Uy(3))/rho
       K2jUj(4) = (1.0d0/3.0d0)*(Cv*mu*(-Ux(1)*v1*v2 - 2*Ux(2)*v2 + 3*Ux(3)*v1) &
            - Uy(1)*(Cv*mu*(3*V_sq + v2**2) - 3*lambda*(V_sq - e)) + 3*Uy(2)*v1*(Cv*mu &
            - lambda) + Uy(3)*v2*(4*Cv*mu - 3*lambda) + 3*Uy(4)*lambda)/(Cv*rho)

       
       do i = 1, nNode
          iNodeID = this%getNodeID(i)
          lhs%stiffness(1:4,i) = lhs%stiffness(1:4,i) - (integrator%getShapeFunc(k,i)*AiUi &
               + this%dNdx(i)*A1AiUi_theta + this%dNdy(i)*A2AiUi_theta                     &
               + this%nu*(this%dNdx(i)*Ux + this%dNdy(i)*Uy)                               &
               + (this%dNdx(i)*K1jUj + this%dNdy(i)*K2jUj))                                &
               *jacobianDet(k)*integrator%getWeight(k)/processInfo%vect(iNodeID)
       end do
    end do
    deallocate(jacobian, jacobianDet, U)
  end subroutine calculateLocalSystem

  subroutine calculateDt(this, processInfo)
    implicit none
    class(CFDElementDT) , intent(inout)        :: this
    class(ProcessInfoDT), intent(inout)        :: processInfo 
    integer(ikind)                             :: i, nNode
    real(rkind), dimension(:,:,:), allocatable :: jacobian
    real(rkind), dimension(:)    , allocatable :: jacobianDet   
    real(rkind)                                :: alpha, deltaTU 
    real(rkind)                                :: V, deltaTC, dt 
    real(rkind)                                :: Vx, Vy, VxMAX, VyMAX, fSafe
    nNode       = this%getnNode()
    jacobian    = this%geometry%jacobianAtGPoints(this%nodalPoints)
    jacobianDet = this%geometry%jacobianDetAtGPoints(jacobian)
    fSafe       = processInfo%getConstants(1)
    VxMAX       = 0.d0
    VyMAX       = 0.d0
    do i = 1, nNode
       Vx = this%node(i)%ptr%dof(2)%val/this%node(i)%ptr%dof(1)%val
       Vy = this%node(i)%ptr%dof(3)%val/this%node(i)%ptr%dof(1)%val
       if (Vx.gt.VxMAX) VxMAX = Vx
       if (Vy.gt.VyMAX) VyMAX = Vy
    end do
    V       = (VxMAX**2+VyMAX**2)**0.5d0  
    alpha   = min((V*sqrt(jacobianDet(1)))/(2.d0*0.001d0)/3.d0,1.d0)    
    deltaTU = 1.d0/(4.d0*0.001d0/sqrt(jacobianDet(1))**2.d0+ALPHA*V/sqrt(jacobianDet(1)))    
    deltaTC = 1.d0/(4.d0*0.001d0/sqrt(jacobianDet(1))**2.d0)    
    dt      = fSafe/(1.d0/deltaTC+1.d0/deltaTU)   
    call processInfo%setMinimumDt(dt)
    deallocate(jacobian, jacobianDet)
  end subroutine calculateDt

  subroutine calculateTauNu(this, processInfo)
    implicit none
    class(cfdelementDT)        , intent(inout) :: this
    class(processinfoDT)       , intent(inout) :: processInfo
    type(integratorptrDT)                      :: integrator
    integer(ikind)                             :: i, j, nNode, nDof
    real(rkind), dimension(:,:,:), allocatable :: jacobian
    real(rkind), dimension(:,:)  , allocatable :: U
    real(rkind), dimension(:)    , allocatable :: jacobianDet, T
    real(rkind)                                :: rho, h_jgn, h_rgn, h_rgne  
    real(rkind)                                :: V2, deltaTC, dtmin, Tinf
    real(rkind)                                :: Vx, Vy, gamma, Cv, Vc
    real(rkind)                                :: fmu, smu, bi, ci, rrr, zzz, shoc 
    real(rkind)                                :: cte, rhoinf, Tau, t_sugn1, t_sugn2, t_sugn3
    real(rkind)                                :: term1, term2, tr1, tau_sung3, tau_sung3_e
    real(rkind)                                :: h_rgn1, h_rgn2, resumen, term_1, term_2
    real(rkind)                                :: R, temp, drx, dry, dr2, rjy, dtx, dty
    real(rkind)                                :: dt2, rtx, rty, dux, duy, du2, rux, ruy, rjx
    nNode = this%getnNode()
    nDof  = 4 
    allocate(U(4,nNode), T(nNode))
    integrator  = this%getIntegrator()
    jacobian    = this%geometry%jacobianAtGPoints(this%nodalPoints)
    jacobianDet = this%geometry%jacobianDetAtGpoints(jacobian)
    dtmin       = processInfo%getDt()
    cte         = processInfo%getConstants(2)
    R           = this%material%R
    Cv          = this%material%Cv
    gamma       = this%material%gamma
    Tinf        = this%material%T_inf
    rhoinf      = this%material%rho
    Tau         = 0.d0
    rho         = 0.d0
    Vx          = 0.d0
    Vy          = 0.d0
    temp        = 0.d0
    drx         = 0.d0
    dry         = 0.d0
    dTx         = 0.d0
    dTy         = 0.d0
    dUx         = 0.d0
    dUy         = 0.d0
    h_rgne      = 1.d-10
    h_rgn       = 1.d-10
    h_jgn       = 1.d-10
    do i = 1, nNode
       do j = 1, nDof
          U(j,i) = this%node(i)%ptr%dof(j)%val
       end do
    end do
    do i = 1, nNode
       rho = rho +  U(1,i)         
       Vx  = Vx  + (U(2,i)/U(1,i)) 
       Vy  = Vy  + (U(3,i)/U(1,i)) 
    end do
    rho = rho/nNode
    Vx  = Vx/nNode
    Vy  = Vy/nNode
    V2  = sqrt(Vx*Vx+Vy*Vy)
    do i = 1, nNode
       T(i) = ((U(4,i)/U(1,i))-0.5d0*((U(2,i)/U(1,i))**2+(U(3,i)/U(1,i))**2))/Cv
    end do
    do i = 1, nNode
       temp = temp + T(i)                                             
    end do
    temp = temp/nNode
    Vc   = sqrt(gamma*R*temp)
    do i = 1, nNode
       drx  = drx + U(1,i)*this%dNdx(i)
       dry  = dry + U(1,i)*this%dNdy(i)
    end do
    dr2  = sqrt(drx*drx+dry*dry)+1.d-20
    rjx  = drx/dr2
    rjy  = dry/dr2
    do i = 1, nNode
       dTx  = dTx + T(i)*this%dNdx(i)
       dTy  = dTy + T(i)*this%dNdy(i)
    end do
    dT2  = sqrt(dTx*dTx+dTy*dTy)+1.d-20
    rTx  = dTx/dT2
    rTy  = dTy/dT2    
    do i = 1, nNode
       dUx  = dUx + V2*this%dNdx(i)
       dUy  = dUy + V2*this%dNdy(i)
    end do
    dU2  = sqrt(dUx*dUx+dUy*dUy)+1.d-20
    rUx  = dUx/dU2
    rUy  = dUy/dU2 
    smu  = 110.d0
    fmu  = 0.41685d0*(temp/Tinf)**1.5d0*(Tinf+smu)/(temp+smu)   
    do i = 1, nNode
       term_1 = abs(Vx *this%dNdx(i)+Vy *this%dNdy(i))
       term_2 = abs(rjx*this%dNdx(i)+rjy*this%dNdy(i))
       h_rgn1 = abs(rtx*this%dNdx(i)+rty*this%dNdy(i)) 
       h_rgn2 = abs(rux*this%dNdx(i)+ruy*this%dNdy(i)) 
       Tau    = Tau+term_1+term_2*Vc       
       h_rgne = h_rgne+h_rgn1             
       h_rgn  = h_rgn+h_rgn2              
       h_jgn  = h_jgn+term_2              
    end do
    Tau     = 1.d0/Tau
    h_rgne  = 2.d0/h_rgne
    h_rgn   = 2.d0/h_rgn
    if (h_rgn .gt. 1.d3) h_rgn = 0.d0
    h_jgn   = 2.d0/h_jgn       
    if (h_jgn .gt. 1.d10) h_jgn = 0.d0 
    tr1     = dr2*h_jgn/rho
    zzz     = h_jgn/2.d0
    shoc    = (sqrt(tr1)+tr1**2)*Vc*2*zzz
    resumen = (1.d0/Tau)**2.d0 +(2.d0/dtmin)**2.d0
    rrr     = resumen**(-0.5d0)
    t_sugn1 = rrr
    t_sugn2 = rrr
    t_sugn3 = rrr
    if(fmu.ne.0.d0)then                 
       tau_sung3   = h_rgn**2.d0/(4.d0*fmu/rhoinf)
       tau_sung3_e = h_rgne**2.d0/(4.d0*fmu/rhoinf)
       t_sugn2     = (resumen+1.d0/tau_sung3**2.d0)**(-0.5d0)
       t_sugn3     = (resumen+1.d0/tau_sung3_e**2.d0)**(-0.5d0)
    end if
    this%Tau(1) = t_sugn1
    this%Tau(2) = t_sugn2
    this%Tau(3) = t_sugn3
    this%nu     = shoc*cte
    deallocate(jacobian, jacobianDet, U, T)
  end subroutine calculateTauNu

  subroutine calculateResults(this, processInfo, resultMat)
    implicit none
    class(CFDElementDT)                           , intent(inout) :: this
    type(ProcessInfoDT)                           , intent(inout) :: processInfo
    real(rkind)    , dimension(:,:,:), allocatable, intent(inout) :: resultMat
    type(IntegratorPtrDT)                                         :: integrator
    integer(ikind)                                                :: i, j, k, iNodeID
    integer(ikind)                                                :: nNode, nDof
    real(rkind)    , dimension(:,:,:), allocatable                :: jacobian
    real(rkind)    , dimension(:,:)  , allocatable                :: U
    real(rkind)    , dimension(:)    , allocatable                :: jacobianDet
    real(rkind)                                                   :: adder, bi, ci, gamma
    real(rkind)                                                   :: AiUi(4), U_loc(4)
    real(rkind)                                                   :: Ux(4), Uy(4)
    real(rkind)                                                   :: v1, v2, V_sq, e, rho
    if (processInfo%getProcess(1) == 1) then
       call this%vectorCalculus()
       nNode = this%getnNode()
       integrator = this%getIntegrator()
       jacobian = this%geometry%jacobianAtGPoints(this%nodalPoints)
       jacobianDet = this%geometry%jacobianDetAtGPoints(jacobian)
       do i = 1, nNode
          adder = 0.d0
          do j = 1, nNode
             do k = 1, integrator%getIntegTerms()
                adder = adder + (integrator%getWeight(k)*jacobianDet(k)&
                     *integrator%getShapeFunc(k,i)*integrator%getShapeFunc(k,j))
             end do
          end do
          iNodeID = this%getNodeID(i)
          processInfo%vect(iNodeID) = processInfo%vect(iNodeID) + adder
       end do       
       deallocate(jacobian, jacobianDet)
    end if
    if (processInfo%getProcess(2) == 1) then
       call this%calculateDt(processInfo)
    end if
    if (processInfo%getProcess(3) == 1) then
       nNode = this%getnNode()
       nDof  = 4
       allocate(U(nDof,nNode))
       integrator  = this%getIntegrator()
       jacobian    = this%geometry%jacobianAtGPoints(this%nodalPoints)
       jacobianDet = this%geometry%jacobianDetAtGPoints(jacobian)
       gamma       = this%material%gamma
       Ux          = 0.d0
       Uy          = 0.d0
       U_loc       = 0.d0
       do i = 1, nNode
          do j = 1, nDof
             U(j,i) = this%node(i)%ptr%dof(j)%val
          end do
          Ux = Ux + U(:,i)*this%dNdx(i)
          Uy = Uy + U(:,i)*this%dNdy(i)
       end do
       do k = 1, integrator%getIntegTerms()
          do i = 1, nNode
             U_loc = U_loc + integrator%ptr%shapeFunc(k,i)*U(:,i)  
          end do
          rho  = U_loc(1)
          v1   = U_loc(2)/rho
          v2   = U_loc(3)/rho
          e    = U_loc(4)/rho
          V_sq = v1**2 + v2**2
          AiUi(1) = Ux(2) + Uy(3) 
          AiUi(2) = (1.0d0/2.0d0)*Ux(1)*(V_sq*(gamma - 1) - 2*v1**2) - Ux(2)*v1*( &
               gamma - 3) - Ux(3)*v2*(gamma - 1) + Ux(4)*(gamma - 1) - Uy(1)*v1*v2 + &
               Uy(2)*v2 + Uy(3)*v1
          AiUi(3) = -Ux(1)*v1*v2 + Ux(2)*v2 + Ux(3)*v1 + (1.0d0/2.0d0)*Uy(1)*(V_sq*( &
               gamma - 1) - 2*v2**2) - Uy(2)*v1*(gamma - 1) - Uy(3)*v2*(gamma - 3) + &
               Uy(4)*(gamma - 1)
          AiUi(4) = Ux(1)*v1*(V_sq*(gamma - 1) - e*gamma) - 1.0d0/2.0d0*Ux(2)*(V_sq &
               *(gamma - 1) - 2*e*gamma + 2*v1**2*(gamma - 1)) - Ux(3)*v1*v2*( &
               gamma - 1) + Ux(4)*gamma*v1 + Uy(1)*v2*(V_sq*(gamma - 1) - e*gamma) - &
               Uy(2)*v1*v2*(gamma - 1) - 1.0d0/2.0d0*Uy(3)*(V_sq*(gamma - 1) - 2*e* &
               gamma + 2*v2**2*(gamma - 1)) + Uy(4)*gamma*v2
          do i = 1, nNode
             iNodeID = this%getNodeID(i)
             processInfo%mat(:,iNodeID) = processInfo%mat(:,iNodeID)               &
                  - integrator%getShapeFunc(k,i)*AiUi                              &
                  *jacobianDet(k)*integrator%getWeight(k)/processInfo%vect(iNodeID)
          end do
       end do
       call this%calculateTauNu(processInfo)
       deallocate(jacobian, jacobianDet, U)
    end if
  end subroutine calculateResults
  
  subroutine calculateLHS(this, processInfo, lhs)
    implicit none
    class(CFDElementDT)                                   , intent(inout) :: this
    type(ProcessInfoDT)                                   , intent(inout) :: processInfo
    type(LeftHandSideDT)                                  , intent(inout) :: lhs
  end subroutine calculateLHS

  subroutine calculateRHS(this, processInfo, rhs)
    implicit none
    class(CFDElementDT)                          , intent(inout) :: this
    type(ProcessInfoDT)                               , intent(inout) :: processInfo
    real(rkind)            , dimension(:)  , allocatable, intent(inout) :: rhs
  end subroutine calculateRHS

  subroutine setupIntegration(this, integrator, valuedSource, jacobianDet)
    implicit none
    class(CFDElementDT)                                 , intent(inout) :: this
    type(IntegratorPtrDT)                               , intent(in)    :: integrator
    real(rkind), dimension(4,integrator%getIntegTerms()), intent(out)   :: valuedSource
    real(rkind), dimension(integrator%getIntegTerms())  , intent(out)   :: jacobianDet
    integer(ikind)                                                      :: i, nNode
    type(NodePtrDT), dimension(:), allocatable                          :: nodalPoints
    nNode = this%getnNode()
    valuedSource = this%getValuedSource(integrator)
    jacobianDet = this%geometry%jacobianDetAtGPoints(this%nodalPoints)
  end subroutine setupIntegration

  function getValuedSource(this, integrator)
    implicit none
    class(CFDElementDT), intent(inout) :: this
    type(IntegratorPtrDT) , intent(in) :: integrator
    real(rkind), dimension(4,integrator%getIntegTerms()) :: getValuedSource
    integer(ikind) :: i, j, nNode
    real(rkind) :: x, y
    type(NodePtrDT), dimension(:), allocatable :: node
    nNode = this%getnNode()
    do i = 1, integrator%getIntegTerms()
       node = this%node
       x = 0
       y = 0
       do j = 1, nNode
          x = x + integrator%getShapeFunc(i,j)*node(j)%ptr%getx()
          y = y + integrator%getShapeFunc(i,j)*node(j)%ptr%gety()
       end do
       getValuedSource(1,i) = this%source(1)%evaluate(1, (/x,y/))
       getValuedSource(2,i) = this%source(1)%evaluate(2, (/x,y/))
       getValuedSource(3,i) = this%source(1)%evaluate(3, (/x,y/))
       getValuedSource(4,i) = this%source(1)%evaluate(4, (/x,y/))
    end do
  end function getValuedSource

  subroutine vectorCalculus(this)
    implicit none
    class(CFDElementDT), intent(inout) :: this
    type(IntegratorPtrDT)                          :: integrator
    type(NodePtrDT), dimension(:)    , allocatable :: nodalPoints
    real(rkind)    , dimension(:,:,:), allocatable :: jacobian
    real(rkind)    , dimension(:)    , allocatable :: jacobianDet
    real(rkind)                                    :: bi, ci
    integer(ikind)                                 :: i, j, nNode, nDof
    nNode = this%getnNode()
    nDof  = this%node(1)%getnDof()
    allocate(this%nodalPoints(nNode), this%dNdx(nNode), this%dNdy(nNode))
    do i = 1, nNode
       this%nodalPoints(i) = this%node(i)
    end do
    integrator    = this%getIntegrator()
    jacobian      = this%geometry%jacobianAtGPoints(this%nodalPoints)
    jacobianDet   = this%geometry%jacobianDetAtGPoints(jacobian)
    do i = 1, nNode
       bi = jacobian(1,2,2)*integrator%getDShapeFunc(1,1,i) &
            - jacobian(1,1,2)*integrator%getDShapeFunc(1,2,i)
       ci = jacobian(1,1,1)*integrator%getDShapeFunc(1,2,i) &
            - jacobian(1,2,1)*integrator%getDShapeFunc(1,1,i)
       this%dNdx(i) = bi/jacobianDet(1)
       this%dNdy(i) = ci/jacobianDet(1)
    end do
  end subroutine vectorCalculus
  
end module CFDElementM



