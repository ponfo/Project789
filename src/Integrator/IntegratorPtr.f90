module IntegratorPtrM
  use UtilitiesM
  use DebuggerM
  use IntegratorM
  implicit none
  private
  public :: IntegratorPtrDT
  type :: IntegratorPtrDT
     class(IntegratorDT), pointer :: ptr
   contains
     procedure, public :: allocate
  end type IntegratorPtrDT

contains

  subroutine allocate(this, integrator)
    implicit none
    class(IntegratorPtrDT), intent(inout) :: this
    type(IntegratorDT), target, intent(in) :: integrator
    this%ptr => integrator
  end subroutine allocate
  
end module IntegratorPtrM
