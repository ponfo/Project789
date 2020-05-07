module NodePtrM
  use UtilitiesM
  
  use NodeM

  use PointM

  use DofM

  use SourceM

  implicit none

  private
  public :: NodePtrDT

  type :: NodePtrDT
     class(NodeDT), pointer :: ptr
   contains
     procedure, public :: associate
     procedure, public :: assignSource
     procedure, public :: assignDof
     procedure, public :: fixDof
     procedure, public :: freeDof
     procedure, public :: getnDof
     procedure, public :: setID
     procedure, public :: getID
     procedure, public :: setX
     procedure, public :: setY
     procedure, public :: setZ
     procedure, public :: getX
     procedure, public :: getY
     procedure, public :: getZ
     procedure, public :: getDimension
     procedure, public :: hasSourceOneSource
     procedure, public :: hasSourceMultiSource
     generic           :: hasSource => hasSourceOneSource, hasSourceMultiSource
  end type NodePtrDT

contains

  subroutine associate(this, node)
    implicit none
    class(NodePtrDT)        , intent(inout) :: this
    type(NodeDT)    , target, intent(in)    :: node
    this%ptr => node
  end subroutine associate

  subroutine assignSource(this, iSource, source)
    implicit none
    class(NodePtrDT)       , intent(inout) :: this
    integer(ikind)         , intent(in)    :: iSource
    class(SourceDT), target, intent(in)    :: source
    call this%ptr%assignSource(iSource, source)
  end subroutine assignSource

  subroutine assignDof(this, iDof, dof)
    implicit none
    class(NodePtrDT), intent(inout) :: this
    integer(ikind)  , intent(in)    :: iDof
    real(rkind)     , intent(in)    :: dof
    call this%ptr%assignDof(iDof, dof)
  end subroutine assignDof

  subroutine fixDof(this, iDof, fixedVal)
    implicit none
    class(NodePtrDT), intent(inout) :: this
    integer(ikind)  , intent(in)    :: iDof
    real(rkind)     , intent(in)    :: fixedVal
    call this%ptr%fixDof(iDof, fixedVal)
  end subroutine fixDof

  subroutine freeDof(this, iDof)
    implicit none
    class(NodePtrDT), intent(inout) :: this
    integer(ikind)  , intent(in)    :: iDof
    call this%ptr%freeDof(iDof)
  end subroutine freeDof

  integer(ikind) pure function getnDof(this)
    implicit none
    class(NodePtrDT), intent(in) :: this
    getnDof = this%ptr%getnDof()
  end function getnDof

  subroutine setX(this, x)
    implicit none
    class(NodePtrDT), intent(inout) :: this
    real(rkind)     , intent(in)    :: x
    call this%ptr%setX(x)
  end subroutine setX

  subroutine setY(this, y)
    implicit none
    class(NodePtrDT), intent(inout) :: this
    real(rkind)     , intent(in)    :: y
    call this%ptr%setY(y)
  end subroutine setY
  
  subroutine setZ(this, z)
    implicit none
    class(NodePtrDT), intent(inout) :: this
    real(rkind)     , intent(in)    :: z
    call this%ptr%setZ(z)
  end subroutine setZ
  
  real(rkind) pure function getX(this)
    implicit none
    class(NodePtrDT), intent(in) :: this
    getX = this%ptr%getX()
  end function getX

  real(rkind) pure function getY(this)
    implicit none
    class(NodePtrDT), intent(in) :: this
    getY = this%ptr%getY()
  end function getY

  real(rkind) pure function getZ(this)
    implicit none
    class(NodePtrDT), intent(in) :: this
    getZ = this%ptr%getZ()
  end function getZ

  integer(ikind) pure function getDimension(this)
    implicit none
    class(NodePtrDT), intent(in) :: this
    getDimension = this%ptr%getDimension()
  end function getDimension

  subroutine setID(this, id)
    implicit none
    class(NodePtrDT), intent(inout) :: this
    integer(ikind)  , intent(in)    :: id
    call this%ptr%setID(id)
  end subroutine setID

  integer(ikind) pure function getID(this)
    implicit none
    class(NodePtrDT), intent(in) :: this
    getID = this%ptr%getID()
  end function getID

  logical function hasSourceOneSource(this)
    implicit none
    class(NodePtrDT), intent(inout) :: this
    hasSourceOneSource = this%ptr%hasSourceOneSource()
  end function hasSourceOneSource

  logical function hasSourceMultiSource(this, iSource)
    implicit none
    class(NodePtrDT), intent(inout) :: this
    integer(ikind)  , intent(in)    :: iSource
    hasSourceMultiSource = this%ptr%hasSourceMultiSource(iSource)
  end function hasSourceMultiSource

end module NodePtrM
