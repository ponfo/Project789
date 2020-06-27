module NormalVelocityM
  
  use UtilitiesM

  use PointM
  use NodeM
  use NodePtrM
  use GeometryM

  use LeftHandSideM

  use IntegratorPtrM

  use ProcessInfoM
  
  use ConditionM

  implicit none

  private
  public :: NormalVelocityDT, normalVelocity

  type, extends(ConditionDT) :: NormalVelocityDT
     integer(ikind), dimension(:), allocatable :: nodeIDList
   contains
     procedure, public :: init
     procedure, public :: calculateLocalSystem
     procedure, public :: calculateLHS
     procedure, public :: calculateRHS
  end type NormalVelocityDT

  interface normalVelocity
     procedure :: constructor
  end interface normalVelocity

contains

  type(NormalVelocityDT) function constructor(id, nodeIDList, node, geometry)
    implicit none
    integer(ikind)                           , intent(in) :: id
    integer(ikind)             , dimension(:), intent(in) :: nodeIDList
    type(NodePtrDT)            , dimension(:), intent(in) :: node
    class(GeometryDT)          , pointer     , intent(in) :: geometry
    call constructor%init(id, nodeIDList, node, geometry)
  end function constructor

  subroutine init(this, id, nodeIDList, node, geometry)
    implicit none
    class(NormalVelocityDT)                  , intent(inout) :: this
    integer(ikind)                           , intent(in)    :: id
    integer(ikind)             , dimension(:), intent(in)    :: nodeIDList
    type(NodePtrDT)            , dimension(:), intent(in)    :: node
    class(GeometryDT)          , pointer     , intent(in)    :: geometry
    this%id = id
    this%affectsLHS = .false.
    this%affectsRHS = .true.
    this%nodeIDList = nodeIDList
    this%node = node
    this%geometry => geometry
  end subroutine init

  subroutine calculateLocalSystem(this, processInfo, lhs, rhs)
    implicit none
    class(NormalVelocityDT)                          , intent(inout) :: this
    type(ProcessInfoDT)                              , intent(inout) :: processInfo
    type(LeftHandSideDT)                             , intent(inout) :: lhs
    real(rkind)         , dimension(:)  , allocatable, intent(inout) :: rhs
    call this%calculateRHS(processInfo, rhs)
  end subroutine calculateLocalSystem

  subroutine calculateLHS(this, processInfo, lhs)
    implicit none
    class(NormalVelocityDT)   , intent(inout) :: this
    type(ProcessInfoDT)       , intent(inout) :: processInfo
    type(LeftHandSideDT), intent(inout) :: lhs
    print*, 'No LHS component in normal velocity condition'
  end subroutine calculateLHS
  
 subroutine calculateRHS(this, processInfo, rhs)
    implicit none
    class(NormalVelocityDT)                             , intent(inout) :: this
    type(ProcessInfoDT)                                 , intent(inout) :: processInfo
    real(rkind)          , dimension(:)    , allocatable, intent(inout) :: rhs
    integer(ikind)                                                      :: i, j
    integer(ikind)                                                      :: nNode
    real(rkind)                                                         :: nx, ny
    real(rkind)                                                         :: int1, int2
    real(rkind)          , dimension(:,:,:), allocatable                :: jacobian
    real(rkind)          , dimension(:)    , allocatable                :: jacobianDet
    type(IntegratorPtrDT)                                               :: integrator
    type(NodePtrDT)      , dimension(:)    , allocatable                :: nodalPoints
    nNode = this%getnNode()
    integrator = this%getIntegrator()
    integrator%ptr => this%geometry%boundaryGeometry%integrator
    allocate(rhs(2*nNode))
    allocate(nodalPoints(nNode))
    rhs = 0._rkind
    do i = 1, nNode
       nodalPoints(i) = this%node(i)
    end do
    jacobian = this%geometry%boundaryGeometry%jacobianAtGPoints(nodalPoints)
    jacobianDet = this%geometry%boundaryGeometry%jacobianDetAtGPoints(jacobian)
    do i = 1, nNode
       int1 = 0._rkind
       int2 = 0._rkind
       do j = 1, integrator%getIntegTerms()
          nx = jacobian(j,1,1)/jacobianDet(j)
          ny = jacobian(j,1,2)/jacobianDet(j)
          int1 = int1 + integrator%getWeight(j)*integrator%getShapeFunc(j,i)*nx
          int2 = int2 + integrator%getWeight(j)*integrator%getShapeFunc(j,i)*ny
       end do
       rhs(2*i-1) = rhs(2*i-1) + int1
       rhs(2*i)   = rhs(2*i)   + int2
    end do
  end subroutine calculateRHS

end module NormalVelocityM
