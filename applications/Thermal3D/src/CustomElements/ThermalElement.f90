module ThermalElementM

  use UtilitiesM

  use Tetrahedron3D4NodeM
  use Tetrahedron3D10NodeM
  use Hexahedron3D8NodeM
  use Hexahedron3D20NodeM

  use IntegratorPtrM

  use PointM
  use NodeM
  use NodePtrM

  use LeftHandSideM
  use ProcessInfoM

  use SourceM
  use SourcePtrM
  
  use ElementM

  use ThermalMaterialM

  implicit none

  private
  public :: ThermalElementDT, thermalElement, initGeometries

  type, extends(ElementDT) :: ThermalElementDT
     class(ThermalMaterialDT), pointer :: material
   contains
     procedure, public  :: init
     procedure, public  :: calculateLHS
     procedure, public  :: calculateRHS
     procedure, public  :: calculateLocalSystem
     procedure, public  :: calculateResults
     procedure, private :: setupIntegration
     procedure, private :: getValuedSource
  end type ThermalElementDT

  interface thermalElement
     procedure :: constructor
  end interface thermalElement

  type(Tetrahedron3D4NodeDT) , target, save :: myTetrahedron3D4Node
  type(Tetrahedron3D10NodeDT), target, save :: myTetrahedron3D10Node
  type(Hexahedron3D8NodeDT)  , target, save :: myHexahedron3D8Node
  type(Hexahedron3D20NodeDT) , target, save :: myHexahedron3D20Node

contains

  type(ThermalElementDT) function constructor(id, node, material)
    implicit none
    integer(ikind)                        , intent(in) :: id
    type(NodePtrDT)         , dimension(:), intent(in) :: node
    class(ThermalMaterialDT), target      , intent(in) :: material
    call constructor%init(id, node, material)
  end function constructor

  subroutine init(this, id, node, material)
    implicit none
    class(ThermalElementDT)               , intent(inout) :: this
    integer(ikind)                        , intent(in)    :: id
    type(NodePtrDT)         , dimension(:), intent(in)    :: node
    class(ThermalMaterialDT), target      , intent(in)    :: material
    this%id = id
    this%node = node
    this%material => material
    if(size(node) == 4) then
       this%geometry => myTetrahedron3D4Node
    else if(size(node) == 8) then
       this%geometry => myHexahedron3D8Node
    else if(size(node) == 10) then
       this%geometry => myTetrahedron3D10Node
    else if(size(node) == 20) then
       this%geometry => myHexahedron3D20Node
    end if
    allocate(this%source(1))
  end subroutine init

  subroutine initGeometries(nGauss)
    implicit none
    integer(ikind), intent(in) :: nGauss
    myTetrahedron3D4Node = tetrahedron3D4Node(nGauss)
    myTetrahedron3D10Node = tetrahedron3D10Node(nGauss)
    myHexahedron3D8Node = hexahedron3D8Node(nGauss)
    myHexahedron3D20Node = hexahedron3D20Node(nGauss)
  end subroutine initGeometries

  subroutine calculateLocalSystem(this, processInfo, lhs, rhs)
    implicit none
    class(ThermalElementDT)                               , intent(inout) :: this
    type(ProcessInfoDT)                                   , intent(inout) :: processInfo
    type(LeftHandSideDT)                                  , intent(inout) :: lhs
    real(rkind)            , dimension(:)    , allocatable, intent(inout) :: rhs
    integer(ikind)                                                        :: i, j, k, nNode
    real(rkind)                                                           :: dNidx, dNidy, dNidz
    real(rkind)                                                           :: dNjdx, dNjdy, dNjdz
    real(rkind)            , dimension(:,:,:), allocatable                :: jacobian
    real(rkind)            , dimension(:,:,:), allocatable                :: jacobianInv
    real(rkind)            , dimension(:)    , allocatable                :: jacobianDet
    real(rkind)                                                           :: val
    real(rkind)            , dimension(:)    , allocatable                :: valuedSource
    type(IntegratorPtrDT)                                                 :: integrator
    type(NodePtrDT)        , dimension(:)    , allocatable                :: nodalPoints
    nNode = this%getnNode()
    integrator = this%getIntegrator()
    lhs = leftHandSide(0, 0, nNode)
    allocate(rhs(nNode))
    allocate(nodalPoints(nNode))
    rhs = 0._rkind
    do i = 1, nNode
       nodalPoints(i) = this%node(i)
    end do
    jacobian = this%geometry%jacobianAtGPoints(nodalPoints)
    jacobianDet = this%geometry%jacobianDetAtGPoints(jacobian)
    allocate(jacobianInv(integrator%getIntegTerms(),3,3))
    do i = 1, integrator%getIntegTerms()
       jacobianInv(i,1:3,1:3) = matinv3(jacobian(i,1:3,1:3))
    end do
    do i = 1, nNode
       do j = 1, nNode
          lhs%stiffness(i,j) = 0._rkind
          do k = 1, integrator%getIntegTerms()
             dNidx = jacobianInv(k,1,1)*integrator%getDShapeFunc(k,1,i) &
                  + jacobianInv(k,1,2)*integrator%getDShapeFunc(k,2,i)  &
                  + jacobianInv(k,1,3)*integrator%getDShapeFunc(k,3,i)
             dNidy = jacobianInv(k,2,1)*integrator%getDShapeFunc(k,1,i) &
                  + jacobianInv(k,2,2)*integrator%getDShapeFunc(k,2,i)  &
                  + jacobianInv(k,2,3)*integrator%getDShapeFunc(k,3,i)
             dNidz = jacobianInv(k,3,1)*integrator%getDShapeFunc(k,1,i) &
                  + jacobianInv(k,3,2)*integrator%getDShapeFunc(k,2,i)  &
                  + jacobianInv(k,3,3)*integrator%getDShapeFunc(k,3,i)
             dNjdx = jacobianInv(k,1,1)*integrator%getDShapeFunc(k,1,j) &
                  + jacobianInv(k,1,2)*integrator%getDShapeFunc(k,2,j)  &
                  + jacobianInv(k,1,3)*integrator%getDShapeFunc(k,3,j)
             dNjdy = jacobianInv(k,2,1)*integrator%getDShapeFunc(k,1,j) &
                  + jacobianInv(k,2,2)*integrator%getDShapeFunc(k,2,j)  &
                  + jacobianInv(k,2,3)*integrator%getDShapeFunc(k,3,j)
             dNjdz = jacobianInv(k,3,1)*integrator%getDShapeFunc(k,1,j) &
                  + jacobianInv(k,3,2)*integrator%getDShapeFunc(k,2,j)  &
                  + jacobianInv(k,3,3)*integrator%getDShapeFunc(k,3,j)
             
             lhs%stiffness(i,j) = lhs%stiffness(i,j)                &
                  + integrator%getWeight(k)                         &
                  *(this%material%conductivity(1)*dNidx*dNjdx       &
                  + this%material%conductivity(2)*dNidy*dNjdy       &
                  + this%material%conductivity(3)*dNidz*dNjdz)      &
                  * jacobianDet(k)
          end do
       end do
       if(this%node(i)%hasSource()) then
          val = this%node(i)%ptr%source(1) &
               %evaluate((/this%node(i)%getx(), this%node(i)%gety(), this%node(i)%getz()/))
          rhs(i) = rhs(i) + val
       end if
    end do
    if(this%hasSource()) then
       allocate(valuedSource(integrator%getIntegTerms()))
       valuedSource = this%getValuedSource(integrator)
       do i = 1, nNode
          val = 0._rkind
          do j = 1, integrator%getIntegTerms()
             val = val + integrator%getWeight(j)*integrator%ptr%shapeFunc(j,i) &
                  *valuedSource(j)*jacobianDet(j)
          end do
          rhs(i) = rhs(i) + val
       end do
       deallocate(valuedSource)
       deallocate(jacobianDet)
    end if
    deallocate(jacobianInv)
  end subroutine calculateLocalSystem

  subroutine calculateLHS(this, processInfo, lhs)
    implicit none
    class(ThermalElementDT)                               , intent(inout) :: this
    type(ProcessInfoDT)                                   , intent(inout) :: processInfo
    type(LeftHandSideDT)                                  , intent(inout) :: lhs
    integer(ikind)                                                        :: i, j, k, nNode
    real(rkind)                                                           :: dNidx, dNidy, dNidz
    real(rkind)                                                           :: dNjdx, dNjdy, dNjdz
    real(rkind)            , dimension(:,:,:), allocatable                :: jacobian
    real(rkind)            , dimension(:,:,:), allocatable                :: jacobianInv
    real(rkind)            , dimension(:)    , allocatable                :: jacobianDet
    type(IntegratorPtrDT)                                                 :: integrator
    type(NodePtrDT)        , dimension(:)    , allocatable                :: nodalPoints
    nNode = this%getnNode()
    integrator = this%getIntegrator()
    lhs = leftHandSide(0, 0, nNode)
    allocate(nodalPoints(nNode))
    do i = 1, nNode
       nodalPoints(i) = this%node(i)
    end do
    jacobian = this%geometry%jacobianAtGPoints(nodalPoints)
    jacobianDet = this%geometry%jacobianDetAtGPoints(jacobian)
    allocate(jacobianInv(integrator%getIntegTerms(),3,3))
    do i = 1, integrator%getIntegTerms()
       jacobianInv(i,1:3,1:3) = matinv3(jacobian(i,1:3,1:3))
    end do
    do i = 1, nNode
       do j = 1, nNode
          lhs%stiffness(i,j) = 0._rkind
          do k = 1, integrator%getIntegTerms()
             dNidx = jacobianInv(k,1,1)*integrator%getDShapeFunc(k,1,i) &
                  + jacobianInv(k,1,2)*integrator%getDShapeFunc(k,2,i)  &
                  + jacobianInv(k,1,3)*integrator%getDShapeFunc(k,3,i)
             dNidy = jacobianInv(k,2,1)*integrator%getDShapeFunc(k,1,i) &
                  + jacobianInv(k,2,2)*integrator%getDShapeFunc(k,2,i)  &
                  + jacobianInv(k,2,3)*integrator%getDShapeFunc(k,3,i)
             dNidz = jacobianInv(k,3,1)*integrator%getDShapeFunc(k,1,i) &
                  + jacobianInv(k,3,2)*integrator%getDShapeFunc(k,2,i)  &
                  + jacobianInv(k,3,3)*integrator%getDShapeFunc(k,3,i)
             dNjdx = jacobianInv(k,1,1)*integrator%getDShapeFunc(k,1,j) &
                  + jacobianInv(k,1,2)*integrator%getDShapeFunc(k,2,j)  &
                  + jacobianInv(k,1,3)*integrator%getDShapeFunc(k,3,j)
             dNjdy = jacobianInv(k,2,1)*integrator%getDShapeFunc(k,1,j) &
                  + jacobianInv(k,2,2)*integrator%getDShapeFunc(k,2,j)  &
                  + jacobianInv(k,2,3)*integrator%getDShapeFunc(k,3,j)
             dNjdz = jacobianInv(k,3,1)*integrator%getDShapeFunc(k,1,j) &
                  + jacobianInv(k,3,2)*integrator%getDShapeFunc(k,2,j)  &
                  + jacobianInv(k,3,3)*integrator%getDShapeFunc(k,3,j)
             
             lhs%stiffness(i,j) = lhs%stiffness(i,j)                &
                  + integrator%getWeight(k)                         &
                  *(this%material%conductivity(1)*dNidx*dNjdx       &
                  + this%material%conductivity(2)*dNidy*dNjdy       &
                  + this%material%conductivity(3)*dNidz*dNjdz)      &
                  * jacobianDet(k)
          end do
       end do
    end do
  end subroutine calculateLHS

  subroutine calculateRHS(this, processInfo, rhs)
    implicit none
    class(ThermalElementDT)                           , intent(inout) :: this
    type(ProcessInfoDT)                               , intent(inout) :: processInfo
    real(rkind)            , dimension(:), allocatable, intent(inout) :: rhs
    integer(ikind)                                                    :: i, j, nNode
    real(rkind)                                                       :: val
    real(rkind)            , dimension(:), allocatable                :: valuedSource
    real(rkind)            , dimension(:), allocatable                :: jacobianDet
    type(IntegratorPtrDT)                                             :: integrator
    nNode = this%getnNode()
    allocate(rhs(nNode))
    rhs = 0._rkind
    do i = 1, nNode
       if(this%node(i)%hasSource()) then
          val = this%node(i)%ptr%source(1) &
               %evaluate((/this%node(i)%getx(), this%node(i)%gety(), this%node(i)%getz()/))
          rhs(i) = rhs(i) + val
       end if
    end do
    if(this%hasSource()) then
       integrator = this%getIntegrator()
       allocate(valuedSource(integrator%getIntegTerms()))
       allocate(jacobianDet(integrator%getIntegTerms()))
       call this%setupIntegration(integrator, valuedSource, jacobianDet)
       nNode = this%getnNode()
       do i = 1, nNode
          val = 0._rkind
          do j = 1, integrator%getIntegTerms()
             val = val + integrator%getWeight(j)*integrator%getShapeFunc(j,i) &
                  *valuedSource(j)*jacobianDet(j)
          end do
          rhs(i) = rhs(i) + val
       end do
       deallocate(valuedSource)
       deallocate(jacobianDet)
    end if
  end subroutine calculateRHS

  subroutine setupIntegration(this, integrator, valuedSource, jacobianDet)
    implicit none
    class(ThermalElementDT)                           , intent(inout) :: this
    type(IntegratorPtrDT)                             , intent(in)    :: integrator
    real(rkind), dimension(integrator%getIntegTerms()), intent(out)   :: valuedSource
    real(rkind), dimension(integrator%getIntegTerms()), intent(out)   :: jacobianDet
    integer(ikind)                                                    :: i, nNode
    type(NodePtrDT), dimension(:), allocatable                        :: nodalPoints
    nNode = this%getnNode()
    allocate(nodalPoints(nNode))
    valuedSource = this%getValuedSource(integrator)
    do i = 1, nNode
       nodalPoints(i) = this%node(i)
    end do
    jacobianDet = this%geometry%jacobianDetAtGPoints(nodalPoints)
  end subroutine setupIntegration

  function getValuedSource(this, integrator)
    implicit none
    class(ThermalElementDT), intent(inout) :: this
    type(IntegratorPtrDT) , intent(in) :: integrator
    real(rkind), dimension(integrator%getIntegTerms()) :: getValuedSource
    integer(ikind) :: i, j, nNode
    real(rkind) :: x, y, z
    type(NodePtrDT), dimension(:), allocatable :: node
    nNode = this%getnNode()
    do i = 1, integrator%getIntegTerms()
       node = this%node
       x = 0
       y = 0
       z = 0
       do j = 1, nNode
          x = x + integrator%getShapeFunc(i,j)*node(j)%ptr%getx()
          y = y + integrator%getShapeFunc(i,j)*node(j)%ptr%gety()
          z = z + integrator%getShapeFunc(i,j)*node(j)%ptr%getz()
       end do
       getValuedSource(i) = this%source(1)%evaluate((/x,y,z/))
    end do
  end function getValuedSource

  subroutine calculateResults(this, processInfo, resultMat)
    implicit none
    class(ThermalElementDT)                               , intent(inout) :: this
    type(ProcessInfoDT)                                   , intent(inout) :: processInfo
    real(rkind)            , dimension(:,:,:), allocatable, intent(inout) :: resultMat
    integer(ikind)                                                        :: i, iGauss, nNode
    real(rkind)                                                           :: qx, qy, qz, xi, eta, zeta
    real(rkind)                                                           :: kx, ky, kz
    real(rkind)                                                           :: dNidx, dNidy, dNidz
    real(rkind)            , dimension(:,:,:), allocatable                :: jacobian
    real(rkind)            , dimension(:,:,:), allocatable                :: jacobianInv
    real(rkind)            , dimension(:)    , allocatable                :: jacobianDet
    type(IntegratorPtrDT)                                                 :: integrator
    type(NodePtrDT)        , dimension(:)    , allocatable                :: nodalPoints
    integrator = this%getIntegrator()
    nNode = this%getnNode()
    allocate(nodalPoints(nNode))
    allocate(resultMat(1,integrator%getIntegTerms(),3))
    do i = 1, nNode
       nodalPoints(i) = this%node(i)
    end do
    jacobian = this%geometry%jacobianAtGPoints(nodalPoints)
    jacobianDet = this%geometry%jacobianDetAtGPoints(jacobian)
    allocate(jacobianInv(integrator%getIntegTerms(),3,3))
    do i = 1, integrator%getIntegTerms()
       jacobianInv(i,1:3,1:3) = matinv3(jacobian(i,1:3,1:3))
    end do
    do iGauss = 1, integrator%getIntegTerms()
       xi = integrator%getGPoint(iGauss,1)
       eta = integrator%getGPoint(iGauss,2)
       zeta = integrator%getGPoint(iGauss,3)
       qx = 0._rkind
       qy = 0._rkind
       qz = 0._rkind
       do i = 1, nNode
          dNidx = jacobianInv(iGauss,1,1)*integrator%getDShapeFunc(iGauss,1,i) &
               + jacobianInv(iGauss,1,2)*integrator%getDShapeFunc(iGauss,2,i)  &
               + jacobianInv(iGauss,1,3)*integrator%getDShapeFunc(iGauss,3,i)
          dNidy = jacobianInv(iGauss,2,1)*integrator%getDShapeFunc(iGauss,1,i) &
               + jacobianInv(iGauss,2,2)*integrator%getDShapeFunc(iGauss,2,i)  &
               + jacobianInv(iGauss,2,3)*integrator%getDShapeFunc(iGauss,3,i)
          dNidz = jacobianInv(iGauss,3,1)*integrator%getDShapeFunc(iGauss,1,i) &
               + jacobianInv(iGauss,3,2)*integrator%getDShapeFunc(iGauss,2,i)  &
               + jacobianInv(iGauss,3,3)*integrator%getDShapeFunc(iGauss,3,i)
          qx = qx + dNidx*this%node(i)%ptr%dof(1)%val
          qy = qy + dNidy*this%node(i)%ptr%dof(1)%val
          qz = qz + dNidz*this%node(i)%ptr%dof(1)%val
       end do
       kx = this%material%conductivity(1)
       ky = this%material%conductivity(2)
       kz = this%material%conductivity(3)
       resultMat(1,iGauss,1) = -1._rkind*kx*qx
       resultMat(1,iGauss,2) = -1._rkind*ky*qy
       resultMat(1,iGauss,3) = -1._rkind*kz*qz
    end do
  end subroutine calculateResults
  
end module ThermalElementM
