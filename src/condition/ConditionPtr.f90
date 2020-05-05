module ConditionPtrM
  use UtilitiesM
  use DebuggerM

  use NodeM
  use NodePtrM

  use GeometryM

  use IntegratorPtrM
  
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
     procedure, public :: getIntegrator
     procedure, public :: getNode
     procedure, public :: getNodeID
     procedure, public :: getAffectsLHS
     procedure, public :: getAffectsRHS

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

  integer(ikind) pure function getID(this)
    implicit none
    class(ConditionPtrDT), intent(in) :: this
    getID = this%ptr%getID()
  end function getID

  integer(ikind) pure function getnNode(this)
    implicit none
    class(ConditionPtrDT), intent(in) :: this
    getnNode = this%ptr%getnNode()
  end function getnNode

  type(IntegratorPtrDT) function getIntegrator(this)
    implicit none
    class(ConditionPtrDT), intent(in) :: this
    getIntegrator = this%ptr%getIntegrator()
  end function getIntegrator

  type(NodePtrDT) function getNode(this, iNode)
    implicit none
    class(ConditionPtrDT), intent(inout) :: this
    integer(ikind)       , intent(in)    :: iNode
    getNode = this%ptr%getNode(iNode)
  end function getNode

  integer(ikind) pure function getNodeID(this, iNode)
    implicit none
    class(ConditionPtrDT), intent(in) :: this
    integer(ikind)       , intent(in) :: iNode
    getNodeID = this%ptr%getNodeID(iNode)
  end function getNodeID

  logical pure function getAffectsLHS(this)
    implicit none
    class(ConditionPtrDT), intent(in) :: this
    getAffectsLHS = this%ptr%getAffectsLHS()
  end function getAffectsLHS

  logical pure function getAffectsRHS(this)
    implicit none
    class(ConditionPtrDT), intent(in) :: this
    getAffectsRHS = this%ptr%getAffectsRHS()
  end function getAffectsRHS

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
