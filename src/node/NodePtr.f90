module NodePtrM
  use NodeM

  implicit none

  private
  public :: NodePtrDT

  type :: NodePtrDT
     class(NodeDT), pointer :: ptr
   contains
     procedure, public :: assignSource
     procedure, public :: assignDof
     procedure, public :: fixDof
     procedure, public :: freeDof
     procedure, public :: getnDof
  end type NodePtrDT

contains

  subroutine assignSource(this, source)
    implicit none
    class(NodePtrDT)        , intent(inout) :: this
    class(SourceDT) , target, intent(in)    :: source
    call this%ptr%assignSource(source)
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

  integer(ikind) function getnDof(this)
    implicit none
    class(NodePtrDT), intent(inout) :: this
    getnDof = this%ptr%getnDof()
  end function getnDof

end module NodePtrM
