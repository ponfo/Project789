module FluxOnLineM
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
  public :: FluxOnLineDT, fluxOnLine

  type, extends(ConditionDT) :: FluxOnLineDT
     integer(ikind), dimension(:), allocatable :: nodeIDList
     real(rkind)                               :: flux
   contains
     procedure, public :: init

     procedure, public :: calculateLocalSystem
     procedure, public :: calculateLHS
     procedure, public :: calculateRHS
  end type FluxOnLineDT

  interface fluxOnLine
     procedure :: constructor
  end interface fluxOnLine

contains

  type(FluxOnLineDT) function constructor(id, nodeIDList, flux, node, geometry)
    implicit none
    integer(ikind)                 , intent(in) :: id
    integer(ikind)   , dimension(:), intent(in) :: nodeIDList
    real(rkind)                    , intent(in) :: flux
    type(NodePtrDT)  , dimension(:), intent(in) :: node
    class(GeometryDT), pointer     , intent(in) :: geometry
    call constructor%init(id, nodeIDList, flux, node, geometry)
  end function constructor

  subroutine init(this, id, nodeIDList, flux, node, geometry)
    implicit none
    class(FluxOnLineDT)            , intent(inout) :: this
    integer(ikind)                 , intent(in)    :: id
    integer(ikind)   , dimension(:), intent(in)    :: nodeIDList
    real(rkind)                    , intent(in)    :: flux
    type(NodePtrDT)  , dimension(:), intent(in)    :: node
    class(GeometryDT), pointer     , intent(in)    :: geometry
    this%id = id
    this%affectsLHS = .false.
    this%affectsRHS = .true.
    this%nodeIDList = nodeIDList
    this%flux = flux
    this%node = node
    this%geometry => geometry
  end subroutine init

  subroutine calculateLocalSystem(this, lhs, rhs)
    implicit none
    class(FluxOnLineDT)                                   , intent(inout) :: this
    type(LeftHandSideDT)                                  , intent(inout) :: lhs
    real(rkind)              , dimension(:)  , allocatable, intent(inout) :: rhs
    call this%calculateRHS(rhs)
  end subroutine calculateLocalSystem

  subroutine calculateLHS(this, lhs)
    implicit none
    class(FluxOnLineDT)                                   , intent(inout) :: this
    type(LeftHandSideDT)                                  , intent(inout) :: lhs
    print*, 'No LHS component in FluxOnLine condition'
  end subroutine calculateLHS

  subroutine calculateRHS(this, rhs)
    implicit none
    class(FluxOnLineDT)                                   , intent(inout) :: this
    real(rkind)              , dimension(:)  , allocatable, intent(inout) :: rhs
    integer(ikind)                                                        :: i, j
    integer(ikind)                                                        :: nNode
    real(rkind)                                                           :: int
    real(rkind)              , dimension(:)  , allocatable                :: jacobianDet
    type(IntegratorPtrDT)                                                 :: integrator
    type(NodePtrDT)          , dimension(:)  , allocatable                :: nodalPoints
    nNode = this%getnNode()
    integrator = this%getIntegrator()
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
          int = int + integrator%ptr%weight(j)*integrator%getShapeFunc(j,i) &
               * this%flux*jacobianDet(j)
       end do
       int = int*(-1._rkind)
       rhs(i) = rhs(i) + int
    end do
  end subroutine calculateRHS

end module FluxOnLineM
    
