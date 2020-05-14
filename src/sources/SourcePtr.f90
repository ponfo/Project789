module SourcePtrM
  use UtilitiesM
  use SourceM

  implicit none

  private
  public :: SourcePtrDT

  type :: SourcePtrDT
     class(SourceDT), pointer :: ptr => null()
   contains
     procedure, public :: associate
     procedure, public :: evaluateFirstFunc
     procedure, public :: evaluateAnyFunc
     generic           :: evaluate => evaluateFirstFunc, evaluateAnyFunc
  end type SourcePtrDT

contains

  subroutine associate(this, source)
    implicit none
    class(SourcePtrDT)        , intent(inout) :: this
    type(SourceDT)    , target, intent(in)    :: source
    this%ptr => source
  end subroutine associate

  real(rkind) function evaluateFirstFunc(this, val)
    implicit none
    class(SourcePtrDT)              , intent(inout) :: this
    real(rkind)       , dimension(:), intent(in)    :: val
    evaluateFirstFunc = this%ptr%evaluateFirstFunc(val)
  end function evaluateFirstFunc

  real(rkind) function evaluateAnyFunc(this, iFunc, val)
    implicit none
    class(SourcePtrDT)              , intent(inout) :: this
    integer(ikind)                  , intent(in)    :: iFunc
    real(rkind)       , dimension(:), intent(in)    :: val
    evaluateAnyFunc = this%ptr%evaluateAnyFunc(iFunc, val)
  end function evaluateAnyFunc

end module SourcePtrM
