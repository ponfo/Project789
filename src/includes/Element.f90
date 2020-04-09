module ElementM
  use utilitiesM
  use DebuggerM

  use NodeM
  use NodePtrM

  use GeometryM

  implicit none

  private
  public :: ElementDT

  type, abstract :: ElementDT
     integer(ikind)                               :: id
     type(NodePtrDT)  , dimension(:), allocatable :: node
     class(GeometryDT)              , pointer     :: geometry
   contains
     procedure, public :: assignGeometry
     procedure, public :: assignNode
     procedure, public :: calculateLHS
     procedure, public :: calculateRHS
     procedure, public :: calculateResults
  end type ElementDT

contains

  subroutine assignGeometry(this, geometry)
    implicit none
    class(ElementDT)         , intent(inout) :: this
    class(GeometryDT), target, intent(in)    :: geometry
    this%geometry => geometry
    allocate(this%node(geometry%nNode))
  end subroutine assignGeometry

  subroutine assignNode(this, index, node)
    implicit none
    class(ElementDT)        , intent(inout) :: this
    integer(ikind)          , intent(in)    :: index
    type(NodeDT)    , target, intent(in)    :: node
    this%node(index) => node
  end subroutine assignNode

end module ElementM

  
     
