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
     integer(ikind)                               :: id
     type(NodePtrDT)  , dimension(:), allocatable :: node
     class(GeometryDT)              , pointer     :: geometry
   contains
     procedure, public :: assignGeometry
     procedure, public :: assignNode
     
     procedure, public :: calculateLHS
     procedure, public :: calculateRHS
     procedure, public :: calculateResults
  end type ConditionDT

contains

  subroutine assignGeometry(this, geometry)
    implicit none
    class(ConditionDT)        , intent(inout) :: this
    class(GeometryDT) , target, intent(in)    :: geometry
    this%geometry => geometry
    allocate(this%node(geometry%nNode))
  end subroutine assignGeometry

  subroutine assignNode(this, index, node)
    implicit none
    class(ConditionDT)        , intent(inout) :: this
    integer(ikind)            , intent(in)    :: index
    type(NodeDT)      , target, intent(in)    :: node
    this%node(index)%ptr => node
  end subroutine assignNode

  subroutine calculateLocalSystem(this, lhs, rhs)
    implicit none
    class(ConditionDT)                             , intent(inout) :: this
    real(rkind)       , dimension(:,:), allocatable, intent(out)   :: lhs
    real(rkind)       , dimension(:)  , allocatable, intent(out)   :: rhs
    print*, "** Condition's calculateLocalSystem not implemented **"
  end subroutine calculateLocalSystem

  function calculateLHS(this)
    implicit none
    class(ConditionDT), intent(inout)               :: this
    real(rkind)       , dimension(:,:), allocatable :: calculateLHS
    print*, "** Condition's calculateLHS not implemented **"
  end function calculateLHS

  function calculateRHS(this)
    implicit none
    class(ConditionDT), intent(inout)             :: this
    real(rkind)       , dimension(:), allocatable :: calculateRHS
    print*, "** Condition's calculateRHS not implemented **"
  end function calculateRHS

  subroutine calculateResults(this)
    implicit none
    class(ConditionDT), intent(inout) :: this
    print*, "** Condition's calculateResults not implemented **"
  end subroutine calculateResults

end module ConditionM

  
