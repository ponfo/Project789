module SourceM
  use UtilitiesM
  use FortranParser
  
  implicit none
  
  private
  public :: SourceFunc, sourceFunc, SourceDT, source

  type :: SourceFunc
     integer(ikind)                       :: nDim
     type(FortranParser), dimension(nDim) :: func
   contains
     procedure, public :: initSourceFunc
  end type SourceFunc
  
  type :: SourceDT
     integer(ikind)     , dimension(:), allocatable :: localNodes
     class(ElementDT)                 , pointer     :: element
     type(FortranParser)              , pointer     :: func
   contains
     procedure, public  :: initSource
     procedure, public  :: apply
     procedure, private :: setupIntegration
     procedure, private :: getValuedSource
  end type SourceDT

  interface sourceFunc
     procedure :: sourceFuncConstructor
  end interface sourceFunc

  interface source
     procedure :: sourceConstructor
  end interface source

contains

  type(SourceFuncDT) function sourceFuncConstructor(nVar, nDim, var, func)
    implicit none
    integer(ikind)                   , intent(in) :: nVar
    integer(ikind)                   , intent(in) :: nDim
    character(len=*), dimension(nVar), intent(in) :: var
    character(len=*), dimension(nDim), intent(in) :: func
    call sourceFuncConstructor%initSourceFunc(nVar, nDim, var, func)
  end function sourceFuncConstructor

  subroutine initSourceFunc(this, nVar, nDim, var, func)
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
  end subroutine initSourceFunc

  type(SourceDT) function sourceConstructor(localNodes, element, func)
    implicit none
    integer(ikind)      , dimension(:), intent(in) :: localNodes
    class(ElementDT)    , target      , intent(in) :: element
    class(FortranParser), target      , intent(in) :: func
    call sourceConstructor%initSource(localNodes, element, func)
  end function sourceConstructor

  subroutine initSource(this, localNodes, element, func)
    implicit none
    class(SourceDT)                   , intent(inout) :: this
    integer(ikind)      , dimension(:), intent(in)    :: localNodes
    class(ElementDT)    , target      , intent(in)    :: element
    class(FortranParser), target      , intent(in)    :: func
    this%localNodes = localNodes
    this%element => element
    this%func => func
  end subroutine initSource

  !falta apply, si es que iría acá

end module SourceM

  
    
  
