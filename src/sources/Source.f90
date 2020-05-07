module SourceM
  use UtilitiesM
  use FortranParser
  
  implicit none
  
  private
  public :: SourceDT, source

  type :: SourceDT
     integer(ikind)                                  :: nDim
     type(EquationParser), dimension(:), allocatable :: func
   contains
     procedure, public :: init
     procedure, public :: evaluateFirstFunc
     procedure, public :: evaluateAnyFunc
     generic           :: evaluate => evaluateFirstFunc, evaluateAnyFunc
  end type SourceDT

  interface source
     procedure :: constructor
  end interface source

contains

  type(SourceDT) function constructor(nVar, nDim, var, func)
    implicit none
    integer(ikind)                   , intent(in) :: nVar
    integer(ikind)                   , intent(in) :: nDim
    character(len=*), dimension(nVar), intent(in) :: var
    character(len=*), dimension(nDim), intent(in) :: func
    call constructor%init(nVar, nDim, var, func)
  end function constructor

  subroutine init(this, nVar, nDim, var, func)
    implicit none
    class(SourceDT)                  , intent(inout) :: this
    integer(ikind)                   , intent(in)    :: nVar
    integer(ikind)                   , intent(in)    :: nDim
    character(len=*), dimension(nVar), intent(in)    :: var
    character(len=*), dimension(nDim), intent(in)    :: func
    integer(ikind)                                   :: i
    allocate(this%func(nDim))
    do i = 1, nDim
       this%func = EquationParser(func(i), var)
    end do
  end subroutine init

  real(rkind) function evaluateFirstFunc(this, val)
    implicit none
    class(SourceDT)              , intent(inout) :: this
    real(rkind)    , dimension(:), intent(in)    :: val
    evaluateFirstFunc = this%func(1)%evaluate(val)
  end function evaluateFirstFunc

  real(rkind) function evaluateAnyFunc(this, iFunc, val)
    implicit none
    class(SourceDT)              , intent(inout) :: this
    integer(ikind)               , intent(in)    :: iFunc
    real(rkind)    , dimension(:), intent(in)    :: val
    evaluateAnyFunc = this%func(iFunc)%evaluate(val)
  end function evaluateAnyFunc

end module SourceM

  
    
  
