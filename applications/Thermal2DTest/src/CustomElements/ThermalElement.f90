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

  use SourceM
  
  use ElementM

  use ThermalMaterialM

  implicit none

  private
  public :: ThermalElementDT, thermalElement, initGeometries

  type, extends(ElementDT) :: ThermalElementDT
     type(ThermalMaterialDT), pointer :: material
   contains
     procedure, public  :: init
     procedure, public  :: calculateLHS
     procedure, public  :: calculateRHS
     procedure, public  :: calculateLocalSystem
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
    integer(ikind)                       , intent(in) :: id
    type(NodePtrDT)        , dimension(:), intent(in) :: node
    type(ThermalMaterialDT), target      , intent(in) :: material
    call constructor%init(id, node, material)
  end function constructor

  subroutine init(this, id, node, material)
    implicit none
    class(ThermalElementDT)              , intent(inout) :: this
    integer(ikind)                       , intent(in)    :: id
    type(NodePtrDT)        , dimension(:), intent(in)    :: node
    type(ThermalMaterialDT), target      , intent(in)    :: material
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
    class(ThermalElementDT)                             , intent(inout) :: this
    real(rkind)            , dimension(:,:), allocatable, intent(inout) :: lhs
    real(rkind)            , dimension(:)  , allocatable, intent(inout) :: rhs
    call this%calculateLHS(lhs)
    call this%calculateRHS(rhs)
  end subroutine calculateLocalSystem

  subroutine calculateLHS(this, lhs)
    implicit none
    class(ThermalElementDT)                               , intent(inout) :: this
    real(rkind)            , dimension(:,:)  , allocatable, intent(inout) :: lhs
    integer(ikind)                                                        :: i, j, k, nNode
    real(rkind)                                                           :: bi, bj, ci, cj
    real(rkind)            , dimension(:,:,:), allocatable                :: jacobian
    real(rkind)            , dimension(:)    , allocatable                :: jacobianDet
    type(IntegratorPtrDT)                                                 :: integrator
    type(PointDT)                                                         :: gaussPoint
    type(NodePtrDT)        , dimension(:)    , allocatable                :: nodalPoints
    nNode = this%getnNode()
    integrator%ptr => this%geometry%integrator
    allocate(lhs(nNode, nNode))
    allocate(jacobian(integrator%ptr%integTerms,2,2))
    allocate(jacobianDet(integrator%ptr%integTerms))
    allocate(nodalPoints(nNode))
    integrator%ptr => this%geometry%integrator
    gaussPoint = point(0._rkind, 0._rkind)
    do i = 1, nNode
       nodalPoints(i) = this%node(i)
    end do
    do i = 1, integrator%ptr%integTerms
       call gaussPoint%updatePoint(integrator%ptr%gPoint(i,1),integrator%ptr%gPoint(i,2))
       jacobian(i,1:2,1:2) = this%geometry%jacobian(gaussPoint, nodalPoints)
       jacobianDet(i) = this%geometry%jacobianDet(jacobian(i,1:2,1:2))
    end do
    do i = 1, nNode
       do j = 1, nNode
          lhs(i,j) = 0._rkind
          do k = 1, integrator%ptr%integTerms
             bi = jacobian(k,2,2)*integrator%ptr%dShapeFunc(k,1,i) &
                  - jacobian(k,1,2)*integrator%ptr%dShapeFunc(k,2,i)
             bj = jacobian(k,2,2)*integrator%ptr%dShapeFunc(k,1,j) &
                  - jacobian(k,1,2)*integrator%ptr%dShapeFunc(k,2,j)
             ci = jacobian(k,1,1)*integrator%ptr%dShapeFunc(k,2,i) &
                  - jacobian(k,2,1)*integrator%ptr%dShapeFunc(k,1,i)
             cj = jacobian(k,1,1)*integrator%ptr%dShapeFunc(k,2,j) &
                  - jacobian(k,2,1)*integrator%ptr%dShapeFunc(k,1,j)
             
             lhs(i,j) = lhs(i,j)                            &
                  + integrator%ptr%weight(k)                 &
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
       if(associated(this%node(i)%ptr%source)) then
          val = this%node(i)%ptr%source%func(1)%evaluate((/this%node(i)%getx(), this%node(i)%gety()/))
          rhs(i) = rhs(i) + val
       end if
    end do
    if(associated(this%source)) then
       integrator%ptr => this%geometry%integrator
       allocate(valuedSource(integrator%ptr%integTerms))
       allocate(jacobianDet(integrator%ptr%integTerms))
       call this%setupIntegration(integrator, valuedSource, jacobianDet)
       nNode = this%getnNode()
       do i = 1, nNode
          val = 0._rkind
          do j = 1, integrator%ptr%integTerms
             val = val + integrator%ptr%weight(j)*integrator%ptr%shapeFunc(j,i) &
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
    class(ThermalElementDT)                          , intent(inout) :: this
    type(IntegratorPtrDT)                            , intent(in)    :: integrator
    real(rkind), dimension(integrator%ptr%integTerms), intent(out)   :: valuedSource
    real(rkind), dimension(integrator%ptr%integTerms), intent(out)   :: jacobianDet
    integer(ikind)                                                   :: i, nNode
    real(rkind), dimension(2,2)                                      :: jacobian
    type(PointDT)                                                    :: gaussPoint
    type(NodePtrDT), dimension(:), allocatable                       :: nodalPoints
    nNode = this%getnNode()
    allocate(nodalPoints(nNode))
    valuedSource = this%getValuedSource(integrator)
    gaussPoint = point(0._rkind, 0._rkind)
    do i = 1, nNode
       nodalPoints(i) = this%node(i)
    end do
    do i = 1, integrator%ptr%integTerms
       call gaussPoint%updatePoint(integrator%ptr%gPoint(i,1),integrator%ptr%gPoint(i,2))
       jacobian = this%geometry%jacobian(gaussPoint, nodalPoints)
       jacobianDet(i) = this%geometry%jacobianDet(jacobian)
    end do
  end subroutine setupIntegration

  function getValuedSource(this, integrator)
    implicit none
    class(ThermalElementDT), intent(inout) :: this
    type(IntegratorPtrDT) , intent(in) :: integrator
    real(rkind), dimension(integrator%ptr%integTerms) :: getValuedSource
    integer(ikind) :: i, j, nNode
    real(rkind) :: x, y
    type(NodePtrDT), dimension(:), allocatable :: node
    nNode = this%getnNode()
    do i = 1, integrator%ptr%integTerms
       node = this%node
       x = 0
       y = 0
       do j = 1, nNode
          x = x + integrator%ptr%shapeFunc(i,j)*node(j)%ptr%getx()
          y = y + integrator%ptr%shapeFunc(i,j)*node(j)%ptr%gety()
       end do
       getValuedSource(i) = this%source%func(1)%evaluate((/x,y/))
    end do
  end function getValuedSource
          

end module ThermalElementM
