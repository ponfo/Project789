module ProcessInfoM
  use UtilitiesM
  use DebuggerM

  implicit none

  private
  public :: ProcessInfoDT

  type :: ProcessInfoDT
     real(rkind)   , allocatable :: t
     real(rkind)   , allocatable :: dt
     real(rkind)   , allocatable :: t0
     real(rkind)   , allocatable :: errorTol
     real(rkind)   , allocatable :: printStep
     real(rkind)   , allocatable :: constants(:)
     integer(ikind), allocatable :: step
     integer(ikind), allocatable :: maxIter
   contains
     procedure, public :: setTime
     procedure, public :: setDT
     procedure, public :: setMinimumDT
     procedure, public :: setErrorTol
     procedure, public :: setPrintStep
     procedure, public :: setT0
     procedure, public :: setConstants
     procedure, public :: setStep
     procedure, public :: setMaxIter
     procedure, public :: getTime
     procedure, public :: getDT
     procedure, public :: getErrorTol
     procedure, public :: getPrintStep
     procedure, public :: getT0
     generic  , public :: getConstants => getOneConstant, getAllConstants
     procedure         :: getOneConstant
     procedure         :: getAllConstants
     procedure, public :: getStep
     procedure, public :: getMaxIter
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
    if(.not.allocated(this%printStep)) allocate(this%printStep)
    this%printStep = printStep
  end subroutine setPrintStep

    subroutine setT0(this, t0)
    implicit none
    class(ProcessInfoDT), intent(inout) :: this
    real(rkind)         , intent(in)    :: t0
    if(.not.allocated(this%t0)) allocate(this%t0)
    this%t0 = t0
  end subroutine setT0
  
  subroutine setErrorTol(this, errorTol)
    implicit none
    class(ProcessInfoDT), intent(inout) :: this
    real(rkind)         , intent(in)    :: errorTol
    if(.not.allocated(this%errorTol)) allocate(this%errorTol)
    this%errorTol = errorTol
  end subroutine setErrorTol

  subroutine setConstants(this, n, values)
    implicit none
    class(ProcessInfoDT)                  , intent(inout) :: this
    real(rkind), dimension(:), allocatable, intent(in)    :: values
    integer(ikind)                        , intent(in)    :: n
    allocate(this%constants(n))
    this%constants = values
  end subroutine setConstants

  subroutine setStep(this, step)
    implicit none
    class(ProcessInfoDT), intent(inout) :: this
    integer(ikind)      , intent(in)    :: step
    if(.not.allocated(this%step)) allocate(this%step)
    this%step = step
  end subroutine setStep

  subroutine setMaxIter(this, maxIter)
    implicit none
    class(ProcessInfoDT), intent(inout) :: this
    integer(ikind)      , intent(in)    :: maxIter
    if(.not.allocated(this%maxIter)) allocate(this%maxIter)
    this%maxIter = maxIter
  end subroutine setMaxIter
  
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

  real(rkind) pure function getOneConstant(this,i)
    implicit none
    class(ProcessInfoDT), intent(in) :: this
    integer(ikind)      , intent(in) :: i
    getOneConstant = this%constants(i)
  end function getOneConstant

  pure function getAllConstants(this)
    implicit none
    class(ProcessInfoDT), intent(in)            :: this
    real(rkind), dimension(size(this%constants)) :: getAllConstants
    getAllConstants = this%constants
  end function getAllConstants

  integer(ikind) pure function getStep(this)
    implicit none
    class(ProcessInfoDT), intent(in) :: this
    getStep = this%step
  end function getStep

  integer(ikind) pure function getMaxIter(this)
    implicit none
    class(ProcessInfoDT), intent(in) :: this
    getMaxIter = this%maxIter
  end function getMaxIter
  
end module ProcessInfoM
     
