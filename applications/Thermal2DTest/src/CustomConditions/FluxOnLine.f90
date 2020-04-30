module FluxOnLineM
  use ConditionM

  implicit none

  private
  public :: FluxOnLineDT, fluxOnLine

  type, extends(ConditionDT) :: FluxOnLineDT
     real(rkind) :: flux
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

  type(FluxOnLineDT) function constructor(id, flux, node, geometry)
    implicit none
    integer(ikind)                 , intent(in) :: id
    real(rkind)                    , intent(in) :: flux
    type(NodePtrDT)  , dimension(:), intent(in) :: node
    class(GeometryDT), target      , intent(in) :: geometry
    call constructor%init(id, flux, node, geometry)
  end function constructor

  subroutine init(this, id, flux, node, geometry)
    implicit none
    class(FluxOnLineDT)            , intent(inout) :: this
    integer(ikind)                 , intent(in)    :: id
    real(rkind)                    , intent(in)    :: flux
    type(NodePtrDT)  , dimension(:), intent(in)    :: node
    class(GeometryDT), target      , intent(in)    :: geometry
    this%id = id
    this%flux = flux
    this%node = node
    this%geometry => geometry
  end subroutine init

  subroutine calculateLocalSystem(this, lhs, rhs)
    implicit none
    class(FluxOnLineDT)                                   , intent(inout) :: this
    real(rkind)              , dimension(:,:), allocatable, intent(inout) :: lhs
    real(rkind)              , dimension(:)  , allocatable, intent(inout) :: rhs
    nNode = this%geometry%nNode
    allocate(lhs(nNode,nNode))
    allocate(rhs(nNode))
    print*, 'Reimplementación pendiente'
  end subroutine calculateLocalSystem

  subroutine calculateLHS(this, lhs)
    implicit none
    class(FluxOnLineDT)                                   , intent(inout) :: this
    real(rkind)              , dimension(:,:), allocatable, intent(inout) :: lhs
    print*, 'Reimplementación pendiente'
  end subroutine calculateLHS

  subroutine calculateRHS(this, rhs)
    implicit none
    class(FluxOnLineDT)                                   , intent(inout) :: this
    real(rkind)              , dimension(:)  , allocatable, intent(inout) :: rhs
    print*, 'Reimplementación pendiente'
  end subroutine calculateRHS

end module FluxOnLineM
    
