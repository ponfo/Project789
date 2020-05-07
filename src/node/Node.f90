module NodeM
  use UtilitiesM

  use PointM

  use DofM

  use SourceM
  use SourcePtrM

  implicit none

  private
  public :: NodeDT, node

  type, extends(PointDT) :: NodeDT
     integer(ikind)                               :: id
     type(DofDT)      , dimension(:), allocatable :: dof
     type(SourcePtrDT), dimension(:), allocatable :: source
   contains
     procedure, public :: initNode1DOneSource
     procedure, public :: initNode1DMultiSource
     procedure, public :: initNode2DOneSource
     procedure, public :: initNode2DMultiSource
     procedure, public :: initNode3DOneSource
     procedure, public :: initNode3DMultiSource
     procedure, public :: assignSourceOne
     procedure, public :: assignSourceMulti
     procedure, public :: assignDof
     procedure, public :: fixDof
     procedure, public :: freeDof
     procedure, public :: getnDof
     procedure, public :: setID
     procedure, public :: getID
     procedure, public :: hasSourceOneSource
     procedure, public :: hasSourceMultiSource
     generic           :: assignSource => assignSourceOne, assignSourceMulti
     generic           :: hasSource => hasSourceOneSource, hasSourceMultiSource
  end type NodeDT

  interface node
     procedure :: constructor1DOneSource
     procedure :: constructor1DMultiSource
     procedure :: constructor2DOneSource
     procedure :: constructor2DMultiSource
     procedure :: constructor3DOneSource
     procedure :: constructor3DMultiSource
  end interface node
  
contains

  type(NodeDT) function constructor1DOneSource(id, nDof, x)
    implicit none
    integer(ikind), intent(in) :: id
    integer(ikind), intent(in) :: nDof 
    real(rkind)   , intent(in) :: x
    call constructor1DOneSource%initNode1DOneSource(id, nDof, x)
  end function constructor1DOneSource
  subroutine initNode1DOneSource(this, id, nDof, x)
    implicit none
    class(NodeDT) , intent(inout) :: this
    integer(ikind), intent(in)    :: id
    integer(ikind), intent(in)    :: nDof
    real(rkind)   , intent(in)    :: x
    call this%initPoint1D(x)
    call this%setID(id)
    allocate(this%dof(nDof))
    allocate(this%source(1))
  end subroutine initNode1DOneSource

  type(NodeDT) function constructor1DMultiSource(id, nDof, nSource, x)
    implicit none
    integer(ikind), intent(in) :: id
    integer(ikind), intent(in) :: nDof
    integer(ikind), intent(in) :: nSource
    real(rkind)   , intent(in) :: x
    call constructor1DMultiSource%initNode1DMultiSource(id, nDof, nSource, x)
  end function constructor1DMultiSource
  subroutine initNode1DMultiSource(this, id, nDof, nSource, x)
    implicit none
    class(NodeDT) , intent(inout) :: this
    integer(ikind), intent(in)    :: id
    integer(ikind), intent(in)    :: nDof
    integer(ikind), intent(in)    :: nSource
    real(rkind)   , intent(in)    :: x
    call this%initPoint1D(x)
    call this%setID(id)
    allocate(this%dof(nDof))
    allocate(this%source(nSource))
  end subroutine initNode1DMultiSource

  type(NodeDT) function constructor2DOneSource(id, nDof, x, y)
    implicit none
    integer(ikind), intent(in) :: id
    integer(ikind), intent(in) :: nDof
    real(rkind)   , intent(in) :: x
    real(rkind)   , intent(in) :: y
    call constructor2DOneSource%initNode2DOneSource(id, nDof, x, y)
  end function constructor2DOneSource
  subroutine initNode2DOneSource(this, id, nDof, x, y)
    implicit none
    class(NodeDT) , intent(inout) :: this
    integer(ikind), intent(in)    :: id
    integer(ikind), intent(in)    :: nDof
    real(rkind)   , intent(in)    :: x
    real(rkind)   , intent(in)    :: y
    call this%initPoint2D(x, y)
    call this%setID(id)
    allocate(this%dof(nDof))
    allocate(this%source(1))
  end subroutine initNode2DOneSource

  type(NodeDT) function constructor2DMultiSource(id, nDof, nSource, x, y)
    implicit none
    integer(ikind), intent(in) :: id
    integer(ikind), intent(in) :: nDof
    integer(ikind), intent(in) :: nSource
    real(rkind)   , intent(in) :: x
    real(rkind)   , intent(in) :: y
    call constructor2DMultiSource%initNode2DMultiSource(id, nDof, nSource, x, y)
  end function constructor2DMultiSource
  subroutine initNode2DMultiSource(this, id, nDof, nSource, x, y)
    implicit none
    class(NodeDT) , intent(inout) :: this
    integer(ikind), intent(in)    :: id
    integer(ikind), intent(in)    :: nDof
    integer(ikind), intent(in)    :: nSource
    real(rkind)   , intent(in)    :: x
    real(rkind)   , intent(in)    :: y
    call this%initPoint2D(x, y)
    call this%setID(id)
    allocate(this%dof(nDof))
    allocate(this%source(nSource))
  end subroutine initNode2DMultiSource

  type(NodeDT) function constructor3DOneSource(id, nDof, x, y, z)
    implicit none
    integer(ikind), intent(in) :: id
    integer(ikind), intent(in) :: nDof
    real(rkind)   , intent(in) :: x
    real(rkind)   , intent(in) :: y
    real(rkind)   , intent(in) :: z
    call constructor3DOneSource%initNode3DOneSource(id, nDof, x, y, z)
  end function constructor3DOneSource
  subroutine initNode3DOneSource(this, id, nDof, x, y, z)
    implicit none
    class(NodeDT) , intent(inout) :: this
    integer(ikind), intent(in)    :: id
    integer(ikind), intent(in)    :: nDof
    real(rkind)   , intent(in)    :: x
    real(rkind)   , intent(in)    :: y
    real(rkind)   , intent(in)    :: z
    call this%initPoint3D(x, y, z)
    call this%setID(id)
    allocate(this%dof(nDof))
    allocate(this%source(1))
  end subroutine initNode3DOneSource

  type(NodeDT) function constructor3DMultiSource(id, nDof, nSource, x, y, z)
    implicit none
    integer(ikind), intent(in) :: id
    integer(ikind), intent(in) :: nDof
    integer(ikind), intent(in) :: nSource
    real(rkind)   , intent(in) :: x
    real(rkind)   , intent(in) :: y
    real(rkind)   , intent(in) :: z
    call constructor3DMultiSource%initNode3DMultiSource(id, nDof, nSource, x, y, z)
  end function constructor3DMultiSource
  subroutine initNode3DMultiSource(this, id, nDof, nSource, x, y, z)
    implicit none
    class(NodeDT) , intent(inout) :: this
    integer(ikind), intent(in)    :: id
    integer(ikind), intent(in)    :: nDof
    integer(ikind), intent(in)    :: nSource
    real(rkind)   , intent(in)    :: x
    real(rkind)   , intent(in)    :: y
    real(rkind)   , intent(in)    :: z
    call this%initPoint3D(x, y, z)
    call this%setID(id)
    allocate(this%dof(nDof))
    allocate(this%source(nSource))
  end subroutine initNode3DMultiSource

  subroutine assignSourceOne(this, source)
    implicit none
    class(NodeDT)          , intent(inout) :: this
    class(SourceDT), target, intent(in)    :: source
    call this%source(1)%associate(source)
  end subroutine assignSourceOne

  subroutine assignSourceMulti(this, iSource, source)
    implicit none
    class(NodeDT)          , intent(inout) :: this
    integer(ikind)         , intent(in)    :: iSource
    class(SourceDT), target, intent(in)    :: source
    call this%source(iSource)%associate(source)
  end subroutine assignSourceMulti

  ! iDof  -> índice del dof en el nodo
  ! index -> índice del dof en el vector de dofs
  subroutine assignDof(this, iDof, dof)
    implicit none
    class(NodeDT) , intent(inout) :: this
    integer(ikind), intent(in)    :: iDof
    real(rkind)   , intent(in)    :: dof
    this%dof(iDof) = newDof(dof, .false.)
  end subroutine assignDof

  subroutine fixDof(this, iDof, fixedVal)
    implicit none
    class(NodeDT) , intent(inout) :: this
    integer(ikind), intent(in)    :: iDof
    real(rkind)   , intent(in)    :: fixedVal
    call this%dof(iDof)%fixDof(fixedVal)
  end subroutine fixDof

  subroutine freeDof(this, iDof)
    implicit none
    class(NodeDT) , intent(inout) :: this
    integer(ikind), intent(in)    :: iDof
    call this%dof(iDof)%freeDof()
  end subroutine freeDof

  integer(ikind) pure function getnDof(this)
    implicit none
    class(NodeDT), intent(in) :: this
    getnDof = size(this%dof)
  end function getnDof

  subroutine setID(this, id)
    implicit none
    class(NodeDT) , intent(inout) :: this
    integer(ikind), intent(in)    :: id
    this%id = id
  end subroutine setID

  integer(ikind) pure function getID(this)
    implicit none
    class(NodeDT), intent(in) :: this
    getID = this%id
  end function getID

  logical function hasSourceOneSource(this)
    implicit none
    class(NodeDT), intent(inout) :: this
    hasSourceOneSource = associated(this%source(1)%ptr)
  end function hasSourceOneSource

  logical function hasSourceMultiSource(this, iSource)
    implicit none
    class(NodeDT) , intent(inout) :: this
    integer(ikind), intent(in)    :: iSource
    hasSourceMultiSource = associated(this%source(iSource)%ptr)
  end function hasSourceMultiSource
  
end module NodeM
