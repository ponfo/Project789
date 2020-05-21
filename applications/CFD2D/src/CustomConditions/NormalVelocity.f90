module NormalVelocityM
  
  use UtilitiesM

  use PointM
  use NodeM
  use NodePtrM
  use GeometryM

  use LeftHandSideM

  use IntegratorPtrM
  
  use ConditionM

  implicit none

  private
  public :: NormalVelocityDT, normalVelocity

  type, extends(ConditionDT) :: NormalVelocityDT
     integer(ikind), dimension(:), allocatable :: nodeIDList
     real(rkind)                               :: velocity
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

  type(NormalVelocityDT) function constructor(id, nodeIDList, velocity, node, geometry)
    implicit none
    integer(ikind)                           , intent(in) :: id
    integer(ikind)             , dimension(:), intent(in) :: nodeIDList
    real(rkind)                              , intent(in) :: pressure
    type(NodePtrDT)            , dimension(:), intent(in) :: node
    class(GeometryDT)          , pointer     , intent(in) :: geometry
    call constructor%init(id, nodeIDList, pressure, node, geometry)
  end function constructor

  subroutine init(this, id, nodeIDList, velocity, node, geometry)
    implicit none
    class(NormalVelocityDT)                        , intent(inout) :: this
    integer(ikind)                           , intent(in)    :: id
    integer(ikind)             , dimension(:), intent(in)    :: nodeIDList
    real(rkind)                              , intent(in)    :: velocity
    type(NodePtrDT)            , dimension(:), intent(in)    :: node
    class(GeometryDT)          , pointer     , intent(in)    :: geometry
    this%id = id
    this%affectsLHS = .false.
    this%affectsRHS = .true.
    this%nodeIDList = nodeIDList
    this%velocity = velocity
    this%node = node
    this%geometry => geometry
  end subroutine init

  subroutine calculateLocalSystem(this, lhs, rhs)
    implicit none
    class(NormalVelocityDT)                                , intent(inout) :: this
    type(LeftHandSideDT)                             , intent(inout) :: lhs
    real(rkind)         , dimension(:)  , allocatable, intent(inout) :: rhs
    call this%calculateRHS(rhs)
  end subroutine calculateLocalSystem

  subroutine calculateLHS(this, lhs)
    implicit none
    class(NormalVelocityDT)   , intent(inout) :: this
    type(LeftHandSideDT), intent(inout) :: lhs
    print*, 'No LHS component in velocity condition'
  end subroutine calculateLHS
  
 subroutine calculateRHS(this, rhs)
    implicit none
    class(NormalVelocityDT)                             , intent(inout) :: this
    real(rkind)          , dimension(:)    , allocatable, intent(inout) :: rhs
    integer(ikind)                                                      :: i, j
    integer(ikind)                                                      :: nNode
    real(rkind)                                                         :: velocityx, velocityy
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
          velocityx = this%velocity*jacobian(j,1,2)/jacobianDet(j)
          velocityy = this%velocity*(-jacobian(j,1,1))/jacobianDet(j)
          int1 = int1 + integrator%getWeight(j)*integrator%getShapeFunc(j,i)  &
               * velocityx*jacobianDet(j)
          int2 = int2 + integrator%getWeight(j)*integrator%getShapeFunc(j,i)  &
               * velocityy*jacobianDet(j)
       end do
       rhs(2*i-1) = rhs(2*i-1) - int1
       rhs(2*i)   = rhs(2*i)   - int2
    end do
  end subroutine calculateRHS

end module NormalVelocityM
