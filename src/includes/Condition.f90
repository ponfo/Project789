module ConditionM
  use UtilitiesM
  use DebuggerM

  use NodeM
  use NodePtrM

  use GeometryM

  implicit none

  private
  public :: ConditionDT

  type, abstract :: ConditionDT
     integer(ikind) :: id
     type(NodePtrDT)  , dimension(:), allocatable :: node
     class(GeometryDT)              , pointer     :: geometry
   contains
     procedure, public :: assignGeometry
     procedure, public :: assignNode
     procedure, public :: calculateLHS
     procedure, public :: calculateRHS
  end type ConditionDT

contains

  subroutine assignGeometry(this, geometry)
    implicit none
    class(ConditionDT)        , intent(inout) :: this
    class(GeometryDT) , target, intent(in)    :: geometry
    this%geometry => geometry
    allocate(this%node(geometry%nNode))
  end subroutine assignGeometry

  
