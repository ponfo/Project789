module FluxOnLineM
  use UtilitiesM

  use PointM
  use NodeM
  use NodePtrM
  use GeometryM

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
    real(rkind)              , dimension(:,:), allocatable, intent(inout) :: lhs
    real(rkind)              , dimension(:)  , allocatable, intent(inout) :: rhs
    call this%calculateRHS(rhs)
  end subroutine calculateLocalSystem

  subroutine calculateLHS(this, lhs)
    implicit none
    class(FluxOnLineDT)                                   , intent(inout) :: this
    real(rkind)              , dimension(:,:), allocatable, intent(inout) :: lhs
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
    type(PointDT)                                                         :: gaussPoint
    type(NodePtrDT)        , dimension(:)    , allocatable                :: nodalPoints
    nNode = this%getnNode()
    integrator%ptr => this%geometry%integrator
    allocate(rhs(nNode))
    allocate(jacobianDet(integrator%ptr%integTerms))
    allocate(nodalPoints(nNode))
    gaussPoint = point(0._rkind, 0._rkind)
    do i = 1, nNode
       nodalPoints(i) = this%node(i)
    end do
    do i = 1, integrator%ptr%integTerms
       call gaussPoint%updatePoint(integrator%ptr%gPoint(i,1),integrator%ptr%gPoint(i,2))
       jacobianDet(i) = this%geometry%jacobianDet(this%nodeIDList, gaussPoint, nodalPoints)
    end do
    do i = 1, nNode
       int = 0._rkind
       do j = 1, integrator%ptr%integTerms
          int = int + integrator%ptr%weight(j)*integrator%ptr%shapeFunc(j,this%nodeIDList(i)) &
               * this%flux*jacobianDet(j)
       end do
       int = int*(-1._rkind)
    end do
  end subroutine calculateRHS

end module FluxOnLineM
    
