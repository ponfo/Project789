module DofM
  use UtilitiesM

  implicit none

  private
  public :: DofDT, dof

  type DofDT
     real(rkind), pointer     :: val
     real(rkind), allocatable :: fixedVal
     logical                  :: isFixed
   contains
     procedure, public :: init
     procedure, public :: fixDof
     procedure, public :: freeDof
  end type DofDT

  interface dof
     procedure :: constructor
  end interface dof

contains

  type(DofDT) function contructor(dof, isFixed)
    implicit none
    real(rkind), target, intent(in) :: dof
    logical            , intent(in) :: isFixed
    call constructor%init(dof, isFixed)
  end function contructor

  subroutine init(this, dof, isFixed)
    implicit none
    class(DofDT)          , intent(inout) :: this
    real(rkind)   , target, intent(in)    :: dof
    logical               , intent(in)    :: isFixed
    this%val     => dof
    this%isFixed = isFixed
  end subroutine init

  subroutine fixDof(this, fixedVal)
    implicit none
    class(DofDT), intent(inout) :: this
    real(rkind) , intent(in)    :: fixedVal
    this%isFixed = .true.
    allocate(this%fixedVal, source = fixedVal)
  end subroutine fixDof

  subroutine freeDof(this)
    implicit none
    class(DofDT), intent(inout) :: this
    this%isFixed = .false.
    deallocate(this%fixedVal)
  end subroutine freeDof
  
end module DofM
