module ProcessInfoM
  use UtilitiesM
  use DebuggerM

  implicit none

  private
  public :: ProcessInfoDT

  type :: ProcessInfoDT
     real(rkind), allocatable :: t
     real(rkind), allocatable :: dt
   contains
     procedure, public :: setTime
     procedure, public :: setDT
     procedure, public :: getTime
     procedure, public :: getDT
  end type ProcessInfoDT

contains

  subroutine setTime(this, t)
    implicit none
    class(ProcessInfoDT), intent(inout) :: this
    real(rkind)         , intent(in)    :: t
    if(.not.allocated(this%t)) allocate(this%t)
    this%t = t
  end subroutine setTime

  subroutine setDT(this, dt)
    implicit none
    class(ProcessInfoDT), intent(inout) :: this
    real(rkind)         , intent(in)    :: dt
    if(.not.allocated(this%dt)) allocate(this%dt)
    this%dt = dt
  end subroutine setDT

  real(rkind) pure function getTime(this)
    implicit none
    class(ProcessInfoDT), intent(in) :: this
    getTime = this%t
  end function getTime

  real(rkind) pure function getDT(this)
    implicit none
    class(ProcessInfoDT), intent(in) :: this
    getDT = this%dt
  end function getDT

end module ProcessInfoM
     
