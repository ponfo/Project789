module ConditionM
  use UtilitiesM
  use DebuggerM

  use IntegratorPtrM

  use NodeM
  use NodePtrM

  use GeometryM

  use LeftHandSideM

  use ProcessInfoM

  implicit none

  private
  public :: ConditionDT

  type, abstract :: ConditionDT
     integer(ikind)                               :: id
     logical                                      :: affectsLHS
     logical                                      :: affectsRHS
     type(NodePtrDT)  , dimension(:), allocatable :: node
     class(GeometryDT)              , pointer     :: geometry
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
     
     procedure(calculateLocalSystemInterf), deferred :: calculateLocalSystem
     procedure(calculateLHSInterf)        , deferred :: calculateLHS
     procedure(calculateRHSInterf)        , deferred :: calculateRHS
  end type ConditionDT

  abstract interface
     subroutine calculateLocalSystemInterf(this, processInfo, lhs, rhs)
       use UtilitiesM
       import ConditionDT
       import ProcessInfoDT
       import LeftHandSideDT
       implicit none
       class(ConditionDT)                               , intent(inout) :: this
       type(ProcessInfoDT)                              , intent(inout) :: processInfo
       type(LeftHandSideDT)                             , intent(inout) :: lhs
       real(rkind)         , dimension(:)  , allocatable, intent(inout) :: rhs
     end subroutine calculateLocalSystemInterf
  end interface

  abstract interface
     subroutine calculateRHSInterf(this, processInfo, rhs)
       use UtilitiesM
       import ConditionDT
       import ProcessInfoDT
       implicit none
       class(ConditionDT)                           , intent(inout) :: this
       type(ProcessInfoDT)                          , intent(inout) :: processInfo
       real(rkind)       , dimension(:), allocatable, intent(inout) :: rhs
     end subroutine calculateRHSInterf
  end interface

  abstract interface
     subroutine calculateLHSInterf(this, processInfo, lhs)
       use UtilitiesM
       import ConditionDT
       import ProcessInfoDT
       import LeftHandSideDT
       implicit none
       class(ConditionDT)  , intent(inout) :: this
       type(ProcessInfoDT) , intent(inout) :: processInfo
       type(LeftHandSideDT), intent(inout) :: lhs
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

  integer(ikind) pure function getID(this)
    implicit none
    class(ConditionDT), intent(in) :: this
    getID = this%id
  end function getID

  integer(ikind) pure function getnNode(this)
    implicit none
    class(ConditionDT), intent(in) :: this
    getnNode = size(this%node)
  end function getnNode

  type(IntegratorPtrDT) function getIntegrator(this)
    implicit none
    class(ConditionDT), intent(in) :: this
    getIntegrator%ptr => this%geometry%integrator
  end function getIntegrator

  type(NodePtrDT) function getNode(this, iNode)
    implicit none
    class(ConditionDT), intent(inout) :: this
    integer(ikind)    , intent(in)    :: iNode
    getNode = this%node(iNode)
  end function getNode

  integer(ikind) pure function getNodeID(this, iNode)
    implicit none
    class(ConditionDT), intent(in) :: this
    integer(ikind)    , intent(in) :: iNode
    getNodeID = this%node(iNode)%ptr%getID()
  end function getNodeID

  logical pure function getAffectsLHS(this)
    implicit none
    class(ConditionDT), intent(in) :: this
    getAffectsLHS = this%affectsLHS
  end function getAffectsLHS

  logical pure function getAffectsRHS(this)
    implicit none
    class(ConditionDT), intent(in) :: this
    getAffectsRHS = this%affectsRHS
  end function getAffectsRHS

end module ConditionM

  
