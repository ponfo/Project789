module ConditionPtrM
  use UtilitiesM
  use DebuggerM

  use NodeM
  use NodePtrM

  use GeometryM
  
  use ConditionM

  implicit none

  private
  public :: ConditionPtrDT

  type :: ConditionPtrDT
     class(ConditionDT), pointer :: ptr
   contains
     procedure, public :: assignGeometry
     procedure, public :: assignNode

     procedure, public :: getID
     procedure, public :: getnNode
     procedure, public :: getNode

     procedure, public :: calculateLocalSystem
     procedure, public :: calculateLHS
     procedure, public :: calculateRHS
  end type ConditionPtrDT

contains

  subroutine assignGeometry(this, geometry)
    implicit none
    class(ConditionPtrDT)        , intent(inout) :: this
    class(GeometryDT)    , target, intent(in)    :: geometry
    call this%ptr%assignGeometry(geometry)
  end subroutine assignGeometry

  subroutine assignNode(this, index, node)
    implicit none
    class(ConditionPtrDT)        , intent(inout) :: this
    integer(ikind)               , intent(in)    :: index
    type(NodeDT)         , target, intent(in)    :: node
    call this%ptr%assignNode(index, node)
  end subroutine assignNode

  integer(ikind) function getID(this)
    implicit none
    class(ConditionPtrDT), intent(inout) :: this
    getID = this%ptr%getID()
  end function getID

  integer(ikind) function getnNode(this)
    implicit none
    class(ConditionPtrDT), intent(inout) :: this
    getnNode = this%ptr%getnNode()
  end function getnNode

  type(NodePtrDT) function getNode(this, iNode)
    implicit none
    class(ConditionPtrDT), intent(inout) :: this
    integer(ikind)       , intent(in)    :: iNode
    getNode = this%ptr%getNode(iNode)
  end function getNode

  subroutine calculateLocalSystem(this, lhs, rhs)
    implicit none
    class(ConditionPtrDT)                             , intent(inout) :: this
    real(rkind)          , dimension(:,:), allocatable, intent(inout) :: lhs
    real(rkind)          , dimension(:)  , allocatable, intent(inout) :: rhs
    call this%ptr%calculateLocalSystem(lhs, rhs)
  end subroutine calculateLocalSystem

  subroutine calculateLHS(this, lhs)
    implicit none
    class(ConditionPtrDT)                             , intent(inout) :: this
    real(rkind)          , dimension(:,:), allocatable, intent(inout) :: lhs
    call this%ptr%calculateLHS(lhs)
  end subroutine calculateLHS

  subroutine calculateRHS(this, rhs)
    implicit none
    class(ConditionPtrDT)                           , intent(inout) :: this
    real(rkind)          , dimension(:), allocatable, intent(inout) :: rhs
    call this%ptr%calculateRHS(rhs)
  end subroutine calculateRHS

end module ConditionPtrM
