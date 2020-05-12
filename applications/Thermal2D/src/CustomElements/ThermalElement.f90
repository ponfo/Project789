module ThermalElementM

  use UtilitiesM

  use Triangle2D3NodeM
  use Triangle2D6NodeM
  use Quadrilateral2D4NodeM
  use Quadrilateral2D8NodeM

  use IntegratorPtrM

  use PointM
  use NodeM
  use NodePtrM

  use LeftHandSideM

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
     procedure, public  :: calculateDT
     procedure, private :: setupIntegration
     procedure, private :: getValuedSource
  end type ThermalElementDT

  interface thermalElement
     procedure :: constructor
  end interface thermalElement

  type(Triangle2D3NodeDT)     , target, save :: myTriangle2D3Node
  type(Triangle2D6NodeDT)     , target, save :: myTriangle2D6Node
  type(Quadrilateral2D4NodeDT), target, save :: myQuadrilateral2D4Node
  type(Quadrilateral2D8NodeDT), target, save :: myQuadrilateral2D8Node

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

  subroutine calculateLocalSystem(this, lhs, rhs)
    implicit none
    class(ThermalElementDT)                               , intent(inout) :: this
    type(LeftHandSideDT)                                  , intent(inout) :: lhs
    real(rkind)            , dimension(:)    , allocatable, intent(inout) :: rhs
    integer(ikind)                                                        :: i, j, k, nNode
    real(rkind)                                                           :: bi, bj, ci, cj
    real(rkind)            , dimension(:,:,:), allocatable                :: jacobian
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
    do i = 1, nNode
       do j = 1, nNode
          lhs%stiffness(i,j) = 0._rkind
          do k = 1, integrator%getIntegTerms()
             bi = jacobian(k,2,2)*integrator%getDShapeFunc(k,1,i) &
                  - jacobian(k,1,2)*integrator%getDShapeFunc(k,2,i)
             bj = jacobian(k,2,2)*integrator%getDShapeFunc(k,1,j) &
                  - jacobian(k,1,2)*integrator%getDShapeFunc(k,2,j)
             ci = jacobian(k,1,1)*integrator%getDShapeFunc(k,2,i) &
                  - jacobian(k,2,1)*integrator%getDShapeFunc(k,1,i)
             cj = jacobian(k,1,1)*integrator%getDShapeFunc(k,2,j) &
                  - jacobian(k,2,1)*integrator%getDShapeFunc(k,1,j)
             
             lhs%stiffness(i,j) = lhs%stiffness(i,j)     &
                  + integrator%getWeight(k)              &
                  *(this%material%conductivity(1)*bi*bj  &
                  + this%material%conductivity(2)*ci*cj) &
                  / jacobianDet(k)
          end do
       end do
       if(this%node(i)%hasSource()) then
          val = this%node(i)%ptr%source(1)%evaluate((/this%node(i)%getx(), this%node(i)%gety()/))
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
  end subroutine calculateLocalSystem

  subroutine calculateLHS(this, lhs)
    implicit none
    class(ThermalElementDT)                               , intent(inout) :: this
    type(LeftHandSideDT)                                  , intent(inout) :: lhs
    integer(ikind)                                                        :: i, j, k, nNode
    real(rkind)                                                           :: bi, bj, ci, cj
    real(rkind)            , dimension(:,:,:), allocatable                :: jacobian
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
    do i = 1, nNode
       do j = 1, nNode
          lhs%stiffness(i,j) = 0._rkind
          do k = 1, integrator%getIntegTerms()
             bi = jacobian(k,2,2)*integrator%getDShapeFunc(k,1,i) &
                  - jacobian(k,1,2)*integrator%getDShapeFunc(k,2,i)
             bj = jacobian(k,2,2)*integrator%getDShapeFunc(k,1,j) &
                  - jacobian(k,1,2)*integrator%getDShapeFunc(k,2,j)
             ci = jacobian(k,1,1)*integrator%getDShapeFunc(k,2,i) &
                  - jacobian(k,2,1)*integrator%getDShapeFunc(k,1,i)
             cj = jacobian(k,1,1)*integrator%getDShapeFunc(k,2,j) &
                  - jacobian(k,2,1)*integrator%getDShapeFunc(k,1,j)
             
             lhs%stiffness(i,j) = lhs%stiffness(i,j)     &
                  + integrator%getWeight(k)              &
                  *(this%material%conductivity(1)*bi*bj  &
                  + this%material%conductivity(2)*ci*cj) &
                  / jacobianDet(k)
          end do
       end do
    end do
  end subroutine calculateLHS

  subroutine calculateRHS(this, rhs)
    implicit none
    class(ThermalElementDT)                           , intent(inout) :: this
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
          val = this%node(i)%ptr%source(1)%evaluate((/this%node(i)%getx(), this%node(i)%gety()/))
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
       getValuedSource(i) = this%source(1)%evaluate((/x,y/))
    end do
  end function getValuedSource

  subroutine calculateResults(this, resultMat)
    implicit none
    class(ThermalElementDT)                               , intent(inout) :: this
    real(rkind)            , dimension(:,:,:), allocatable, intent(inout) :: resultMat
    integer(ikind)                                                        :: i, iGauss, nNode
    real(rkind)                                                           :: bi, ci, qx, qy, xi, eta
    real(rkind)                                                           :: dNidx, dNidy, kx, ky
    real(rkind)            , dimension(:,:,:), allocatable                :: jacobian
    real(rkind)            , dimension(:)    , allocatable                :: jacobianDet
    type(IntegratorPtrDT)                                                 :: integrator
    type(NodePtrDT)        , dimension(:)    , allocatable                :: nodalPoints
    integrator = this%getIntegrator()
    nNode = this%getnNode()
    allocate(nodalPoints(nNode))
    allocate(resultMat(1,integrator%getIntegTerms(),2))
    do i = 1, nNode
       nodalPoints(i) = this%node(i)
    end do
    jacobian = this%geometry%jacobianAtGPoints(nodalPoints)
    jacobianDet = this%geometry%jacobianDetAtGPoints(jacobian)
    do iGauss = 1, integrator%getIntegTerms()
       xi = integrator%getGPoint(iGauss,1)
       eta = integrator%getGPoint(iGauss,2)
       qx = 0._rkind
       qy = 0._rkind
       do i = 1, nNode
          bi = jacobian(iGauss,2,2)*integrator%getDShapeFunc(iGauss,1,i) &
               - jacobian(iGauss,1,2)*integrator%getDShapeFunc(iGauss,2,i)
          ci = jacobian(iGauss,1,1)*integrator%getDShapeFunc(iGauss,2,i) &
               - jacobian(iGauss,2,1)*integrator%getDShapeFunc(iGauss,1,i)
          dNidx = bi/jacobianDet(iGauss)
          dNidy = ci/jacobianDet(iGauss)
          qx = qx + dNidx*this%node(i)%ptr%dof(1)%val
          qy = qy + dNidy*this%node(i)%ptr%dof(1)%val
       end do
       kx = this%material%conductivity(1)
       ky = this%material%conductivity(2)
       resultMat(1,iGauss,1) = -1._rkind*kx*qx
       resultMat(1,iGauss,2) = -1._rkind*ky*qy
    end do
  end subroutine calculateResults

  subroutine calculateDT(this, dt)
    implicit none
    class(ThermalElementDT), intent(inout) :: this
    real(rkind)            , intent(inout) :: dt
  end subroutine calculateDT

end module ThermalElementM
