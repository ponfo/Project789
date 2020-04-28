module ProcessM

  implicit none

  private
  public :: ProcessDT, NewProcessDT

  type, abstract :: NewProcessDT
   contains
     procedure(NewProcess_Procedure), deferred :: useProcess
  end type NewProcessDT

  abstract interface
     subroutine NewProcess_Procedure(this)
       import NewProcessDT
       class(NewProcessDT), intent(inout) :: this
     end subroutine NewProcess_Procedure
  end interface
  
  type ProcessDT
     class(NewProcessDT), allocatable :: process
   contains
     procedure :: init
     procedure :: change
     procedure :: run
  end type ProcessDT

contains

  subroutine init(this, process)
    implicit none
    class(ProcessDT)   , intent(inout) :: this
    class(NewProcessDT), intent(in)    :: process
    allocate(this%process, source = process)
  end subroutine init

  subroutine change(this, newProcess)
    implicit none
    class(ProcessDT)   , intent(inout) :: this
    class(NewProcessDT), intent(in)    :: newProcess
    deallocate(this%process)
    allocate(this%process, source = newProcess)
  end subroutine change

  subroutine run(this)
    implicit none
    class(ProcessDT), intent(inout) :: this
    call this%process%useProcess()
  end subroutine run
  
end module ProcessM
  
