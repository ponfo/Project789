module ConvectionOnLineM
  use UtilitiesM

  use NodePtrM
  use GeometryM
  
  use ConditionM

  implicit none

  private
  public :: ConvectionOnLineDT, convectionOnLine

  type, extends(ConditionDT) :: ConvectionOnLineDT
     real(rkind) :: coef
     real(rkind) :: temp
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

  type(ConvectionOnLineDT) function constructor(id, coef, temp, node, geometry)
    implicit none
    integer(ikind)                 , intent(in) :: id
    real(rkind)                    , intent(in) :: coef
    real(rkind)                    , intent(in) :: temp
    type(NodePtrDT)  , dimension(:), intent(in) :: node
    class(GeometryDT), pointer     , intent(in) :: geometry
    call constructor%init(id, coef, temp, node, geometry)
  end function constructor

  subroutine init(this, id, coef, temp, node, geometry)
    implicit none
    class(ConvectionOnLineDT)      , intent(inout) :: this
    integer(ikind)                 , intent(in)    :: id
    real(rkind)                    , intent(in)    :: coef
    real(rkind)                    , intent(in)    :: temp
    type(NodePtrDT)  , dimension(:), intent(in)    :: node
    class(GeometryDT), pointer     , intent(in)    :: geometry
    this%id = id
    this%coef = coef
    this%temp = temp
    this%node = node
    this%geometry => geometry
  end subroutine init

  subroutine calculateLocalSystem(this, lhs, rhs)
    implicit none
    class(ConvectionOnLineDT)                             , intent(inout) :: this
    real(rkind)              , dimension(:,:), allocatable, intent(inout) :: lhs
    real(rkind)              , dimension(:)  , allocatable, intent(inout) :: rhs
    integer(ikind)                                                        :: nNode
    nNode = this%geometry%nNode
    allocate(lhs(nNode,nNode))
    allocate(rhs(nNode))
    print*, 'Reimplementación pendiente'
  end subroutine calculateLocalSystem

  subroutine calculateLHS(this, lhs)
    implicit none
    class(ConvectionOnLineDT)                             , intent(inout) :: this
    real(rkind)              , dimension(:,:), allocatable, intent(inout) :: lhs
    integer(ikind)                                                        :: nNode
    print*, 'Reimplementación pendiente'
  end subroutine calculateLHS

  subroutine calculateRHS(this, rhs)
    implicit none
    class(ConvectionOnLineDT)                             , intent(inout) :: this
    real(rkind)              , dimension(:)  , allocatable, intent(inout) :: rhs
    integer(ikind)                                                        :: nNode
    print*, 'Reimplementación pendiente'
  end subroutine calculateRHS

end module ConvectionOnLineM
