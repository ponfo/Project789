module ElementPtrM
  use UtilitiesM
  use DebuggerM

  use NodeM
  use NodePtrM

  use GeometryM
  use GeometryObjectM

  use SourceM
  
  use ElementM

  implicit none

  private
  public :: ElementPtrDT

  type :: ElementPtrDT
     class(ElementDT), pointer :: ptr
   contains
     procedure, public :: assignGeometry
     procedure, public :: assignNode
     procedure, public :: assignSource

     procedure, public :: getID
     procedure, public :: getnNode
     procedure, public :: getNode
     procedure, public :: getNodeID

     procedure, public :: calculateLocalSystem
     procedure, public :: calculateLHS
     procedure, public :: calculateRHS
  end type ElementPtrDT

contains

  subroutine assignGeometry(this, geometry)
    implicit none
    class(ElementPtrDT)      , intent(inout) :: this
    class(GeometryDT), target, intent(in)    :: geometry
    call this%ptr%assignGeometry(geometry)
  end subroutine assignGeometry

  subroutine assignNode(this, index, node)
    implicit none
    class(ElementPtrDT)     , intent(inout) :: this
    integer(ikind)          , intent(in)    :: index
    type(NodeDT)    , target, intent(in)    :: node
    call this%ptr%assignNode(index, node)
  end subroutine assignNode

  subroutine assignSource(this, source)
    implicit none
    class(ElementPtrDT)     , intent(inout) :: this
    class(SourceDT) , target, intent(in)    :: source
    call this%ptr%assignSource(source)
  end subroutine assignSource

  integer(ikind) function getID(this)
    implicit none
    class(ElementPtrDT), intent(inout) :: this
    getID = this%ptr%getID()
  end function getID

  integer(ikind) function getnNode(this)
    implicit none
    class(ElementPtrDT), intent(inout) :: this
    getnNode = this%ptr%getnNode()
  end function getnNode

  type(NodePtrDT) function getNode(this, iNode)
    implicit none
    class(ElementPtrDT), intent(inout) :: this
    integer(ikind)     , intent(in)    :: iNode
    getNode = this%ptr%getNode(iNode)
  end function getNode

  integer(ikind) function getNodeID(this, iNode)
    implicit none
    class(ElementPtrDT), intent(inout) :: this
    integer(ikind)     , intent(in)    :: iNode
    getNodeID = this%ptr%getNodeID(iNode)
  end function getNodeID

  subroutine calculateLocalSystem(this, lhs, rhs)
    implicit none
    class(ElementPtrDT)                             , intent(inout) :: this
    real(rkind)        , dimension(:,:), allocatable, intent(inout) :: lhs
    real(rkind)        , dimension(:)  , allocatable, intent(inout) :: rhs
    call this%ptr%calculateLocalSystem(lhs, rhs)
  end subroutine calculateLocalSystem

  subroutine calculateLHS(this, lhs)
    implicit none
    class(ElementPtrDT)                             , intent(inout) :: this
    real(rkind)        , dimension(:,:), allocatable, intent(inout) :: lhs
    call this%ptr%calculateLHS(lhs)
  end subroutine calculateLHS

  subroutine calculateRHS(this, rhs)
    implicit none
    class(ElementPtrDT)                           , intent(inout) :: this
    real(rkind)        , dimension(:), allocatable, intent(inout) :: rhs
    call this%ptr%calculateRHS(rhs)
  end subroutine calculateRHS
    
end module ElementPtrM
