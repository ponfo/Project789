module ProcessM

  implicit none

  private
  public :: ProcessDT, NewProcessDT, SetNewProcess

  type, abstract :: NewProcessDT
  end type NewProcessDT

  interface SetNewProcess
     procedure :: constructor
  end interface SetNewProcess
  
  type ProcessDT
     class(NewProcessDT), allocatable :: process
   contains
     procedure :: setProcess
     procedure :: change
  end type ProcessDT

contains

  type(ProcessDT) function constructor(process)
    class(NewProcessDT), intent(in) :: process
    call constructor%setProcess(process)
  end function Constructor
  
  subroutine setProcess(this, process)
    implicit none
    class(ProcessDT)   , intent(inout) :: this
    class(NewProcessDT), intent(in)    :: process
    allocate(this%process, source = process)
  end subroutine setProcess

  subroutine change(this, newProcess)
    implicit none
    class(ProcessDT)   , intent(inout) :: this
    class(NewProcessDT), intent(in)    :: newProcess
    deallocate(this%process)
    allocate(this%process, source = newProcess)
  end subroutine change
  
end module ProcessM
  
