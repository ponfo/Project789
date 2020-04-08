module NodeM
  use UtilitiesM

  use PointM

  use DofM

  implicit none

  private
  public :: NodeDT, node

  type, extends(PointDT) :: NodeDT
     type(DofDT), dimension(:), allocatable :: dof
   contains
     procedure, public :: initNode1D
     procedure, public :: initNode2D
     procedure, public :: initNode3D
     procedure, public :: assignDof
     procedure, public :: fixDof
     procedure, public :: unfixDof
     procedure, public :: getnDof
  end type NodeDT

  interface node
     procedure :: constructor1D
     procedure :: constructor2D
     procedure :: constructor3D
  end interface node
  
contains

  type(NodeDT) function constructor1D(id, nDof, x)
    implicit none
    integer(ikind), intent(in) :: id
    integer(ikind), intent(in) :: nDof
    real(rkind)   , intent(in) :: x
    call constructor%initNode1D(id, nDof, x)
  end function constructor1D
  subroutine initNode1D(this, id, nDof, x)
    implicit none
    class(NodeDT) , intent(inout) :: this
    integer(ikind), intent(in)    :: id
    integer(ikind), intent(in)    :: nDof
    real(rkind)   , intent(in)    :: x
    call this%initPoint1D(id, x)
    allocate(this%dof(nDof))
  end subroutine initNode1D

  type(NodeDT) function constructor2D(id, nDof, x, y)
    implicit none
    integer(ikind), intent(in) :: id
    integer(ikind), intent(in) :: nDof
    real(rkind)   , intent(in) :: x
    real(rkind)   , intent(in) :: y
    call constructor%initNode2D(id, nDof, x, y)
  end function constructor2D
  subroutine initNode2D(this, id, nDof, x, y)
    implicit none
    class(NodeDT) , intent(inout) :: this
    integer(ikind), intent(in)    :: id
    integer(ikind), intent(in)    :: nDof
    real(rkind)   , intent(in)    :: x
    real(rkind)   , intent(in)    :: y
    call this%initPoint2D(id, x, y)
    allocate(this%dof(nDof))
  end subroutine initNode2D

  type(NodeDT) function constructor3D(id, nDof, x, y, z)
    implicit none
    integer(ikind), intent(in) :: id
    integer(ikind), intent(in) :: nDof
    real(rkind)   , intent(in) :: x
    real(rkind)   , intent(in) :: y
    real(rkind)   , intent(in) :: z
    call constructor%initNode3D(id, nDof, x, y, z)
  end function constructor3D
  subroutine initNode3D(this, id, nDof, x, y, z)
    implicit none
    class(NodeDT) , intent(inout) :: this
    integer(ikind), intent(in)    :: id
    integer(ikind), intent(in)    :: nDof
    real(rkind)   , intent(in)    :: x
    real(rkind)   , intent(in)    :: y
    real(rkind)   , intent(in)    :: z
    call this%initPoint3D(id, x, y, z)
    allocate(this%dof(nDof))
  end subroutine initNode3D

  ! iDof  -> índice del dof en el nodo
  ! index -> índice del dof en el vector de dofs
  subroutine assignDof(this, iDof, index, dofArray)
    implicit none
    class(NodeDT)              , intent(inout) :: this
    integer(ikind)             , intent(in)    :: iDof
    integer(ikind)             , intent(in)    :: index
    real(rkind)  , dimension(:), intent(in)    :: dofArray
    this%dof(iDof) = dof(index, dofArray, .false.)
  end subroutine assignDof

  subroutine fixDof(this, iDof)
    implicit none
    class(NodeDT) , intent(inout) :: this
    integer(ikind), intent(in)    :: iDof
    call this%dof(iDof)%fixDof()
  end subroutine fixDof

  subroutine unfixDof(this, iDof)
    implicit none
    class(NodeDT) , intent(inout) :: this
    integer(ikind), intent(in)    :: iDof
    call this%dof(iDof)%unfixDof()
  end subroutine unfixDof

  integer(ikind) function getnDof(this)
    implicit none
    class(NodeDT), intent(inout) :: this
    getnDof = size(this%dof)
  end function getnDof
  
end module NodeM
