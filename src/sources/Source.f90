module SourceM
  use UtilitiesM
  use FortranParser
  
  implicit none
  
  private
  public :: SourceDT, source

  type :: SourceDT
     integer(ikind)                       :: nDim
     type(FortranParser), dimension(nDim) :: func
   contains
     procedure, public :: init
  end type SourceDT

  interface source
     procedure :: constructor
  end interface source

contains

  type(SourceFuncDT) function constructor(nVar, nDim, var, func)
    implicit none
    integer(ikind)                   , intent(in) :: nVar
    integer(ikind)                   , intent(in) :: nDim
    character(len=*), dimension(nVar), intent(in) :: var
    character(len=*), dimension(nDim), intent(in) :: func
    call constructor%initSourceFunc(nVar, nDim, var, func)
  end function constructor

  subroutine init(this, nVar, nDim, var, func)
    implicit none
    class(SourceFuncDT)              , intent(inout) :: this
    integer(ikind)                   , intent(in)    :: nVar
    integer(ikind)                   , intent(in)    :: nDim
    character(len=*), dimension(nVar), intent(in)    :: var
    character(len=*), dimension(nDim), intent(in)    :: func
    integer(ikind)                                   :: i
    do i = 1, nDim
       this%func = EquationParser(func(i), var)
    end do
  end subroutine init

  !falta apply, si es que iría acá

end module SourceM

  
    
  
