module ElementM
  use UtilitiesM
  use DebuggerM

  use NodeM
  use NodePtrM

  use GeometryM

  use IntegratorPtrM

  use LeftHandSideM
  
  use SourceM
  use SourcePtrM

  implicit none

  private
  public :: ElementDT

  type, abstract :: ElementDT
     integer(ikind)                               :: id
     type(NodePtrDT)  , dimension(:), allocatable :: node
     class(GeometryDT)              , pointer     :: geometry
     type(SourcePTrDT), dimension(:), allocatable :: source
   contains
     procedure, public :: assignGeometry
     procedure, public :: assignNode
     procedure, public :: assignSourceOne
     procedure, public :: assignSourceMulti
     generic           :: assignSource => assignSourceOne, assignSourceMulti

     procedure, public :: getID
     procedure, public :: getnNode
     procedure, public :: getIntegrator
     procedure, public :: getNode
     procedure, public :: getNodeID
     procedure, public :: hasSourceOneSource
     procedure, public :: hasSourceMultiSource
     generic           :: hasSource => hasSourceOneSource, hasSourceMultiSource

     procedure(calculateLocalSystemInterf), deferred :: calculateLocalSystem
     procedure(calculateLHSInterf)        , deferred :: calculateLHS
     procedure(calculateRHSInterf)        , deferred :: calculateRHS
     procedure(calculateResultsInterf)    , deferred :: calculateResults
  end type ElementDT

  abstract interface
     subroutine calculateLocalSystemInterf(this, lhs, rhs)
       use UtilitiesM
       import ElementDT
       import LeftHandSideDT
       implicit none
       class(ElementDT)                                 , intent(inout) :: this
       type(LeftHandSideDT)                             , intent(inout) :: lhs
       real(rkind)         , dimension(:)  , allocatable, intent(inout) :: rhs
     end subroutine calculateLocalSystemInterf
  end interface

  abstract interface
     subroutine calculateRHSInterf(this, rhs)
       use UtilitiesM
       import ElementDT
       implicit none
       class(ElementDT)                           , intent(inout) :: this
       real(rkind)     , dimension(:), allocatable, intent(inout) :: rhs
     end subroutine calculateRHSInterf
  end interface

  abstract interface
     subroutine calculateLHSInterf(this, lhs)
       use UtilitiesM
       import ElementDT
       import LeftHandSideDT
       implicit none
       class(ElementDT)    , intent(inout) :: this
       type(LeftHandSideDT), intent(inout) :: lhs
     end subroutine calculateLHSInterf
  end interface

  !!In ResultMat(i,j,k) ->
  !!                       i: Number of result
  !!                       j: Number of point where the result lays
  !!                       k: Number of dimensions of the result
  abstract interface
     subroutine calculateResultsInterf(this, resultMat)
       use UtilitiesM
       import ElementDT
       implicit none
       class(ElementDT)                               , intent(inout) :: this
       real(rkind)     , dimension(:,:,:), allocatable, intent(inout) :: resultMat
     end subroutine calculateResultsInterf
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

  subroutine assignSourceOne(this, source)
    implicit none
    class(ElementDT)        , intent(inout) :: this
    class(SourceDT) , target, intent(in)    :: source
    call this%source(1)%associate(source)
  end subroutine assignSourceOne

  subroutine assignSourceMulti(this, iSource, source)
    implicit none
    class(ElementDT)        , intent(inout) :: this
    integer(ikind)          , intent(in)    :: iSource
    class(SourceDT) , target, intent(in)    :: source
    call this%source(iSource)%associate(source)
  end subroutine assignSourceMulti

  integer(ikind) pure function getID(this)
    implicit none
    class(ElementDT), intent(in) :: this
    getID = this%id
  end function getID

  integer(ikind) pure function getnNode(this)
    implicit none
    class(ElementDT), intent(in) :: this
    getnNode = size(this%node)
  end function getnNode

  type(IntegratorPtrDT) function getIntegrator(this)
    implicit none
    class(ElementDT), intent(inout) :: this
    getIntegrator%ptr => this%geometry%integrator
  end function getIntegrator

  type(NodePtrDT) function getNode(this, iNode)
    implicit none
    class(ElementDT), intent(inout) :: this
    integer(ikind)  , intent(in)    :: iNode
    getNode = this%node(iNode)
  end function getNode

  integer(ikind) pure function getNodeID(this, iNode)
    implicit none
    class(ElementDT), intent(in) :: this
    integer(ikind)  , intent(in) :: iNode
    getNodeID = this%node(iNode)%ptr%getID()
  end function getNodeID

  logical function hasSourceOneSource(this)
    implicit none
    class(ElementDT), intent(inout) :: this
    hasSourceOneSource = associated(this%source(1)%ptr)
  end function hasSourceOneSource

  logical function hasSourceMultiSource(this, iSource)
    implicit none
    class(ElementDT), intent(inout) :: this
    integer(ikind)  , intent(in)    :: iSource
    hasSourceMultiSource = associated(this%source(iSource)%ptr)
  end function hasSourceMultiSource

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

  
     
