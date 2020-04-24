module ElementM
  use UtilitiesM
  use DebuggerM

  use NodeM
  use NodePtrM

  use GeometryM
  use GeometryObjectM

  use SourceM

  implicit none

  private
  public :: ElementDT

  type, abstract :: ElementDT
     integer(ikind)                               :: id
     type(NodePtrDT)  , dimension(:), allocatable :: node
     class(GeometryDT)              , pointer     :: geometry
     type(SourceDT)                 , pointer     :: source
   contains
     procedure, public :: assignGeometry
     procedure, public :: assignNode
     procedure, public :: assignSource

     procedure                                   :: calculateLocalSystem
     procedure(calculateLHSInterf)    , deferred :: calculateLHS
     procedure(calculateRHSInterf)    , deferred :: calculateRHS
     procedure                                   :: calculateResults
  end type ElementDT

  abstract interface
     function calculateRHSInterf(this)
       use UtilitiesM
       import ElementDT
       implicit none
       class(ElementDT), intent(inout) :: this
       real(rkind)     , dimension(:), allocatable :: calculateRHSInterf
     end function calculateRHSInterf
  end interface

  abstract interface
     function calculateLHSInterf(this)
       use UtilitiesM
       import ElementDT
       implicit none
       class(ElementDT), intent(inout) :: this
       real(rkind)     , dimension(:,:), allocatable :: calculateLHSInterf
     end function calculateLHSInterf
  end interface

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
    this%node(index)%ptr => node
  end subroutine assignNode

  subroutine assignSource(this, source)
    implicit none
    class(ElementDT)        , intent(inout) :: this
    class(SourceDT) , target, intent(in)    :: source
    this%source => source
  end subroutine assignSource

  subroutine calculateLocalSystem(this, lhs, rhs)
    implicit none
    class(ElementDT)                               , intent(inout) :: this
    real(rkind)       , dimension(:,:), allocatable, intent(out)   :: lhs
    real(rkind)       , dimension(:)  , allocatable, intent(out)   :: rhs
    print*, "** Element's calculateLocalSystem not implemented **"
  end subroutine calculateLocalSystem

  subroutine calculateResults(this)
    implicit none
    class(ElementDT), intent(inout) :: this
    print*, "** Element's calculateResults not implemented **"
  end subroutine calculateResults

end module ElementM

  
     
