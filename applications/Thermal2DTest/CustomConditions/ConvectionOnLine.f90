module ConvectionOnLineM
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
    class(GeometryDT), target      , intent(in) :: geometry
    call constructor%init(id, coef, temp, node, geometry)
  end function constructor

  subroutine init(this, id, coef, temp, node, geometry)
    implicit none
    class(ConvectionOnLineDT)      , intent(inout) :: this
    integer(ikind)                 , intent(in)    :: id
    real(rkind)                    , intent(in)    :: coef
    real(rkind)                    , intent(in)    :: temp
    type(NodePtrDT)  , dimension(:), intent(in)    :: node
    class(GeometryDT), target      , intent(in)    :: geometry
    this%id = id
    this%coef = coef
    this%temp = temp
    this%node = node
    this%geometry => geometry
  end subroutine init

  subroutine calculateLocalSystem(this, lhs, rhs)
    implicit none
    class(ConvectionOnLineDT), intent(inout) :: this
    real(rkind), dimension(:,:), allocatable :: lhs
    real(rkind), dimension(:)  , allocatable :: rhs
    nNode = this%geometry%nNode
    allocate(lhs(nNode,nNode))
    allocate(rhs(nNode))
    print*, 'Reimplementación pendiente'
  end subroutine calculateLocalSystem

end module ConvectionOnLineM
