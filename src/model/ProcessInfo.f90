module ProcessInfoM
  use UtilitiesM
  use DebuggerM

  implicit none

  private
  public :: ProcessInfoDT

  type :: ProcessInfoDT
     real(rkind), allocatable :: t
     real(rkind), allocatable :: dt
     real(rkind), allocatable :: t0
     real(rkind), allocatable :: errorTol
     real(rkind), allocatable :: printStep
   contains
     procedure, public :: setTime
     procedure, public :: setDT
     procedure, public :: setMinimumDT
     procedure, public :: setErrorTol
     procedure, public :: setPrintStep
     procedure, public :: setT0
     procedure, public :: getTime
     procedure, public :: getDT
     procedure, public :: getErrorTol
     procedure, public :: getPrintStep
     procedure, public :: getT0
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

  subroutine setMinimumDT(this, dt)
    implicit none
    class(ProcessInfoDT), intent(inout) :: this
    real(rkind)         , intent(in)    :: dt
    if(.not.allocated(this%dt)) then
       allocate(this%dt)
       this%dt = dt
    else if(this%dt > dt) then
       this%dt = dt
    end if
  end subroutine setMinimumDT

  subroutine setPrintStep(this, printStep)
    implicit none
    class(ProcessInfoDT), intent(inout) :: this
    integer(ikind)      , intent(in)    :: printStep
    this%printStep = printStep
  end subroutine setPrintStep

    subroutine setT0(this, t0)
    implicit none
    class(ProcessInfoDT), intent(inout) :: this
    real(rkind)         , intent(in)    :: t0
    this%t0 = t0
  end subroutine setT0
  
  subroutine setErrorTol(this, errorTol)
    implicit none
    class(ProcessInfoDT), intent(inout) :: this
    real(rkind)         , intent(in)    :: errorTol
    this%errorTol = errorTol
  end subroutine setErrorTol
  
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

  integer(ikind) pure function getPrintStep(this)
    implicit none
    class(ProcessInfoDT), intent(in) :: this
    getPrintStep = this%printStep
  end function getPrintStep

  real(rkind) pure function getT0(this)
    implicit none
    class(ProcessInfoDT), intent(in) :: this
    getT0 = this%t0
  end function getT0

  real(rkind) pure function getErrorTol(this)
    implicit none
    class(ProcessInfoDT), intent(in) :: this
    getErrorTol = this%errorTol
  end function getErrorTol
  
end module ProcessInfoM
     
