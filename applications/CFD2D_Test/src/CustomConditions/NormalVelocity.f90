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
    real(rkind)          , dimension(:,:)  , allocatable                :: jacobian
    real(rkind)                                                         :: jacobianDet
    type(IntegratorPtrDT)                                               :: integrator
    type(PointDT)        , dimension(:)    , allocatable                :: pointToValue
    type(NodePtrDT)      , dimension(:)    , allocatable                :: nodalPoints
    nNode = this%getnNode()
    integrator = this%getIntegrator()
    integrator%ptr => this%geometry%boundaryGeometry%integrator
    allocate(rhs(2*nNode))
    allocate(nodalPoints(nNode))
    allocate(pointToValue(nNode))
    rhs = 0._rkind
    do i = 1, nNode
       nodalPoints(i) = this%node(i)
       pointToValue(i) = point(this%node(i)%getx(), this%node(i)%gety())
    end do
    
    do i = 1, nNode
       jacobian = this%geometry%boundaryGeometry%jacobian(pointToValue(i), nodalPoints)
       jacobianDet = this%geometry%boundaryGeometry%jacobianDet(jacobian)
       nx = jacobian(1,2)/jacobianDet
       ny = -jacobian(1,1)/jacobianDet
       !Que apunten para adentro:
       nx = -1*nx
       ny = -1*ny
       rhs(2*i-1) = nx
       rhs(2*i)   = ny
    end do
  end subroutine calculateRHS

end module NormalVelocityM
