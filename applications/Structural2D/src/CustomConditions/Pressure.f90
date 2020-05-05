module PressureM
  use UtilitiesM

  use PointM
  use NodeM
  use NodePtrM
  use GeometryM

  use IntegratorPtrM
  
  use ConditionM

  use StructuralMaterialM

  implicit none

  private
  public :: PressureDT, pressure

  type, extends(ConditionDT) :: PressureDT
     integer(ikind), dimension(:), allocatable :: nodeIDList
     real(rkind)                               :: pressure
     type(StructuralMaterialDT)  , pointer     :: material
   contains
     procedure, public :: init

     procedure, public :: calculateLocalSystem
     procedure, public :: calculateLHS
     procedure, public :: calculateRHS
  end type PressureDT

  interface pressure
     procedure :: constructor
  end interface pressure

contains

  type(PressureDT) function constructor(id, nodeIDList, pressure, node, geometry, material)
    implicit none
    integer(ikind)                          , intent(in) :: id
    integer(ikind)            , dimension(:), intent(in) :: nodeIDList
    real(rkind)                             , intent(in) :: pressure
    type(NodePtrDT)           , dimension(:), intent(in) :: node
    class(GeometryDT)         , pointer     , intent(in) :: geometry
    type(StructuralMaterialDT), target      , intent(in) :: material
    call constructor%init(id, nodeIDList, pressure, node, geometry, material)
  end function constructor

  subroutine init(this, id, nodeIDList, pressure, node, geometry, material)
    implicit none
    class(PressureDT)                       , intent(inout) :: this
    integer(ikind)                          , intent(in)    :: id
    integer(ikind)            , dimension(:), intent(in)    :: nodeIDList
    real(rkind)                             , intent(in)    :: pressure
    type(NodePtrDT)           , dimension(:), intent(in)    :: node
    class(GeometryDT)         , pointer     , intent(in)    :: geometry
    type(StructuralMaterialDT), pointer     , intent(in)    :: material
    this%id = id
    this%affectsLHS = .false.
    this%affectsRHS = .true.
    this%nodeIDList = nodeIDList
    this%pressure = pressure
    this%node = node
    this%geometry => geometry
    this%material => material
  end subroutine init

  subroutine calculateLocalSystem(this, lhs, rhs)
    implicit none
    class(PressureDT)                                     , intent(inout) :: this
    real(rkind)              , dimension(:,:), allocatable, intent(inout) :: lhs
    real(rkind)              , dimension(:)  , allocatable, intent(inout) :: rhs
    call this%calculateRHS(rhs)
  end subroutine calculateLocalSystem

  subroutine calculateLHS(this, lhs)
    implicit none
    class(PressureDT)                                     , intent(inout) :: this
    real(rkind)              , dimension(:,:), allocatable, intent(inout) :: lhs
    print*, 'No LHS component in pressure condition'
  end subroutine calculateLHS
  
 subroutine calculateRHS(this, rhs)
    implicit none
    class(PressureDT)                                   , intent(inout) :: this
    real(rkind)          , dimension(:)    , allocatable, intent(inout) :: rhs
    integer(ikind)                                                      :: i, j
    integer(ikind)                                                      :: nNode
    real(rkind)                                                         :: pressurex, pressurey
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
          pressurex = this%pressure*jacobian(j,1,2)*this%material%thickness/jacobianDet(j)
          pressurey = this%pressure*(-jacobian(j,1,1))*this%material%thickness/jacobianDet(j)
          int1 = int1 + integrator%getWeight(j)*integrator%getShapeFunc(j,i)  &
               * pressurex*jacobianDet(j)
          int2 = int2 + integrator%getWeight(j)*integrator%getShapeFunc(j,i)  &
               * pressurey*jacobianDet(j)
       end do
       rhs(2*i-1) = rhs(2*i-1) - int1
       rhs(2*i)   = rhs(2*i)   - int2
    end do
  end subroutine calculateRHS

end module PressureM
