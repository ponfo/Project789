module DofM
  use UtilitiesM

  implicit none

  private
  public :: DofDT, newDof

  type DofDT
     real(rkind), pointer     :: val
     real(rkind), allocatable :: fixedVal
     logical                  :: isFixed
   contains
     procedure, public :: init
     procedure, public :: fixDof
     procedure, public :: freeDof
     procedure, public :: assignVal
     procedure, public :: getVal
  end type DofDT

  interface newDof
     procedure :: constructor
  end interface newDof

contains

  type(DofDT) function constructor(dof, isFixed)
    implicit none
    real(rkind), target, intent(in) :: dof
    logical            , intent(in) :: isFixed
    call constructor%init(dof, isFixed)
  end function constructor

  subroutine init(this, dof, isFixed)
    implicit none
    class(DofDT)          , intent(inout) :: this
    real(rkind)   , target, intent(in)    :: dof
    logical               , intent(in)    :: isFixed
    this%val     => dof
    this%isFixed = isFixed
  end subroutine init

  subroutine assignVal(this, dof)
    implicit none
    class(DofDT)        , intent(inout) ::this
    real(rkind) , target, intent(in)    :: dof
    this%val => dof
  end subroutine assignVal

  subroutine fixDof(this, fixedVal)
    implicit none
    class(DofDT), intent(inout) :: this
    real(rkind) , intent(in)    :: fixedVal
    this%isFixed = .true.
    if(allocated(this%fixedVal)) deallocate(this%fixedVal)
    allocate(this%fixedVal, source = fixedVal)
  end subroutine fixDof

  subroutine freeDof(this)
    implicit none
    class(DofDT), intent(inout) :: this
    this%isFixed = .false.
    if(allocated(this%fixedVal)) deallocate(this%fixedVal)
  end subroutine freeDof

  pure real(rkind) function getVal(this)
    implicit none
    class(DofDT), intent(in) :: this
    if(this%isFixed) then
       getVal = this%fixedVal
    else
       getVal = this%val
    end if
  end function getVal
  
end module DofM
