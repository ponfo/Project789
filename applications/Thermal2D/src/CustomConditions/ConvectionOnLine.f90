module ConvectionOnLineM
  
  use UtilitiesM

  use PointM
  use NodeM
  use NodePtrM
  use GeometryM

  use LeftHandSideM
  use ProcessInfoM

  use IntegratorPtrM
  
  use ConditionM

  implicit none

  private
  public :: ConvectionOnLineDT, convectionOnLine

  type, extends(ConditionDT) :: ConvectionOnLineDT
     integer(ikind), dimension(:), allocatable :: nodeIDList
     real(rkind)                               :: coef
     real(rkind)                               :: temp
   contains
     procedure, public :: init

     procedure, public :: calculateLocalSystem
     procedure, public :: calculateLHS
     procedure, public :: calculateRHS
  end type ConvectionOnLineDT

  interface convectionOnLine
     procedure :: constructor
  end interface convectionOnLine

contains

  type(ConvectionOnLineDT) function constructor(id, nodeIDList, coef, temp, node, geometry)
    implicit none
    integer(ikind)                 , intent(in) :: id
    integer(ikind)   , dimension(:), intent(in) :: nodeIDList
    real(rkind)                    , intent(in) :: coef
    real(rkind)                    , intent(in) :: temp
    type(NodePtrDT)  , dimension(:), intent(in) :: node
    class(GeometryDT), pointer     , intent(in) :: geometry
    call constructor%init(id, nodeIDList, coef, temp, node, geometry)
  end function constructor

  subroutine init(this, id, nodeIDList, coef, temp, node, geometry)
    implicit none
    class(ConvectionOnLineDT)      , intent(inout) :: this
    integer(ikind)                 , intent(in)    :: id
    integer(ikind)   , dimension(:), intent(in)    :: nodeIDList
    real(rkind)                    , intent(in)    :: coef
    real(rkind)                    , intent(in)    :: temp
    type(NodePtrDT)  , dimension(:), intent(in)    :: node
    class(GeometryDT), pointer     , intent(in)    :: geometry
    this%id = id
    this%affectsLHS = .true.
    this%affectsRHS = .true.
    this%nodeIDList = nodeIDList
    this%coef = coef
    this%temp = temp
    this%node = node
    this%geometry => geometry
  end subroutine init

  subroutine calculateLocalSystem(this, processInfo, lhs, rhs)
    implicit none
    class(ConvectionOnLineDT)                             , intent(inout) :: this
    type(ProcessInfoDT)                                   , intent(inout) :: processInfo
    type(LeftHandSideDT)                                  , intent(inout) :: lhs
    real(rkind)              , dimension(:)  , allocatable, intent(inout) :: rhs
    integer(ikind)                                                        :: i, j, k
    integer(ikind)                                                        :: nNode
    real(rkind)                                                           :: int
    real(rkind)              , dimension(:)  , allocatable                :: jacobianDet
    type(IntegratorPtrDT)                                                 :: integrator
    type(NodePtrDT)          , dimension(:)  , allocatable                :: nodalPoints
    nNode = this%getnNode()
    integrator%ptr => this%geometry%boundaryGeometry%integrator
    lhs = leftHandSide(0, 0, nNode)
    allocate(rhs(nNode))
    allocate(nodalPoints(nNode))
    lhs%stiffness = 0._rkind
    rhs = 0._rkind
    do i = 1, nNode
       nodalPoints(i) = this%node(i)
    end do
    jacobianDet = this%geometry%boundaryGeometry%jacobianDetAtGPoints(nodalPoints)
    do i = 1, nNode
       do j = 1, nNode
          int = 0._rkind
          do k = 1, integrator%getIntegTerms()
             int = int + integrator%getWeight(k)*integrator%getShapeFunc(k,i) &
                  * integrator%getShapeFunc(k,j)*this%coef*jacobianDet(k)
          end do
          lhs%stiffness(i,j) = lhs%stiffness(i,j) + int
       end do
    end do
    do i = 1, nNode
       int = 0._rkind
       do j = 1, integrator%getIntegTerms()
          int = int + integrator%getWeight(j)*integrator%getShapeFunc(j,i) &
               * this%coef*this%temp*jacobianDet(j)
       end do
       rhs(i) = rhs(i) + int
    end do
  end subroutine calculateLocalSystem

  subroutine calculateLHS(this, processInfo, lhs)
    implicit none
    class(ConvectionOnLineDT)                             , intent(inout) :: this
    type(ProcessInfoDT)                                   , intent(inout) :: processInfo
    type(LeftHandSideDT)                                  , intent(inout) :: lhs
    integer(ikind)                                                        :: i, j, k
    integer(ikind)                                                        :: nNode
    real(rkind)                                                           :: int
    real(rkind)              , dimension(:)  , allocatable                :: jacobianDet
    type(IntegratorPtrDT)                                                 :: integrator
    type(NodePtrDT)        , dimension(:)    , allocatable                :: nodalPoints
    nNode = this%getnNode()
    integrator%ptr => this%geometry%boundaryGeometry%integrator
    lhs = leftHandSide(0, 0, nNode)
    allocate(nodalPoints(nNode))
    lhs%stiffness = 0._rkind
    do i = 1, nNode
       nodalPoints(i) = this%node(i)
    end do
    jacobianDet = this%geometry%boundaryGeometry%jacobianDetAtGPoints(nodalPoints)
    do i = 1, nNode
       do j = 1, nNode
          int = 0._rkind
          do k = 1, integrator%getIntegTerms()
             int = int + integrator%getWeight(k)*integrator%getShapeFunc(k,i) &
                  * integrator%getShapeFunc(k,j)*this%coef*jacobianDet(k)
          end do
          lhs%stiffness(i,j) = lhs%stiffness(i,j) + int
       end do
    end do
  end subroutine calculateLHS

  subroutine calculateRHS(this, processInfo, rhs)
    implicit none
    class(ConvectionOnLineDT)                             , intent(inout) :: this
    type(ProcessInfoDT)                                   , intent(inout) :: processInfo
    real(rkind)              , dimension(:)  , allocatable, intent(inout) :: rhs
    integer(ikind)                                                        :: i, j
    integer(ikind)                                                        :: nNode
    real(rkind)                                                           :: int
    real(rkind)              , dimension(:)  , allocatable                :: jacobianDet
    type(IntegratorPtrDT)                                                 :: integrator
    type(NodePTrDT)          , dimension(:)  , allocatable                :: nodalPoints
    nNode = this%getnNode()
    integrator%ptr => this%geometry%boundaryGeometry%integrator
    allocate(rhs(nNode))
    allocate(nodalPoints(nNode))
    rhs = 0._rkind
    do i = 1, nNode
       nodalPoints(i) = this%node(i)
    end do
    jacobianDet = this%geometry%boundaryGeometry%jacobianDetAtGPoints(nodalPoints)
    do i = 1, nNode
       int = 0._rkind
       do j = 1, integrator%getIntegTerms()
          int = int + integrator%getWeight(j)*integrator%getShapeFunc(j,i) &
               * this%coef*this%temp*jacobianDet(j)
       end do
       rhs(i) = rhs(i) + int
    end do
  end subroutine calculateRHS

end module ConvectionOnLineM
