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

     procedure, public :: getID
     procedure, public :: getnNode
     procedure, public :: getNode
     
     procedure(calculateLocalSystemInterf), deferred :: calculateLocalSystem
     procedure(calculateLHSInterf)        , deferred :: calculateLHS
     procedure(calculateRHSInterf)        , deferred :: calculateRHS
     procedure                                       :: calculateResults
  end type ConditionDT

  abstract interface
     subroutine calculateLocalSystem(this, lhs, rhs)
       use UtilitiesM
       import ConditionDT
       implicit none
       class(ConditionDT)                             , intent(inout) :: this
       real(rkind)       , dimension(:,:), allocatable, intent(inout) :: lhs
       real(rkind)       , dimension(:)  , allocatable, intent(inout) :: rhs
     end subroutine calculateLocalSystem
  end interface

  abstract interface
     subroutine calculateRHSInterf(this, rhs)
       use UtilitiesM
       import ConditionDT
       implicit none
       class(ConditionDT)                           , intent(inout) :: this
       real(rkind)       , dimension(:), allocatable, intent(inout) :: rhs
     end subroutine calculateRHSInterf
  end interface

  abstract interface
     subroutine calculateLHSInterf(this)
       use UtilitiesM
       import ConditionDT
       implicit none
       class(ConditionDT)                             , intent(inout) :: this
       real(rkind)       , dimension(:,:), allocatable, intent(inout) :: lhs
     end subroutine calculateLHSInterf
  end interface

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

  integer(ikind) function getID(this)
    implicit none
    class(ConditionDT), intent(inout) :: this
    getID = this%id
  end function getID

  integer(ikind) function getnNode(this)
    implicit none
    class(ConditionDT), intent(inout) :: this
    getnNode = size(this%node)
  end function getnNode

  type(NodePtrDT) function getNode(this, iNode)
    implicit none
    class(ConditionDT), intent(inout) :: this
    integer(ikind)    , intent(in)    :: iNode
    getNode = this%node(iNode)
  end function getNode

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

  
