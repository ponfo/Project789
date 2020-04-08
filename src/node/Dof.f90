module DofM
  use UtilitiesM

  implicit none

  private
  public :: DofDT, dof

  type DofDT
     real(rkind), pointer :: val
     logical              :: isFixed
   contains
     procedure, public :: init
     procedure, public :: fixDof
     procedure, public :: freeDof
  end type DofDT

  interface dof
     procedure :: constructor
  end interface dof

contains

  type(DofDT) function contructor(index, dofArray, isFixed)
    implicit none
    integer(ikind)                   , intent(in) :: index
    real(rkind), dimension(:), target, intent(in) :: dofArray
    logical                          , intent(in) :: isFixed
    call constructor%init(index, dofArray, isFixed)
  end function contructor

  subroutine init(this, index, dofArray, isFixed)
    implicit none
    class(DofDT)                     , intent(inout) :: this
    integer(ikind)                   , intent(in)    :: index
    real(rkind), dimension(:), target, intent(in)    :: dofArray
    logical                          , intent(in)    :: isFixed
    this%val     => dofArray(index)
    this%isFixed = isFixed
  end subroutine init

  subroutine fixDof(this)
    implicit none
    class(DofDT), intent(inout) :: this
    this%isFixed = .true.
  end subroutine fixDof

  subroutine freeDof(this)
    implicit none
    class(DofDT), intent(inout) :: this
    this%isFixed = .false.
  end subroutine freeDof
  
end module DofM
