module PreconditionerMethodMod

  use tools
  use sparseKit

  use PreconditionerMod

  implicit none

  private
  public :: PreconditionerMethodTYPE

  type, extends(PreconditionerTYPE) :: PreconditionerMethodTYPE
   contains
     procedure :: usePreconditioner => method
  end type PreconditionerMethodTYPE

contains

  subroutine method(this, vector, matrix, solution, arg)
    implicit none
    class(PreconditionerMethodTYPE), intent(inout) :: this
    class(Sparse)                   , intent(inout) :: matrix
    real(rkind)   , dimension(:)    , intent(inout) :: vector
    real(rkind)   , dimension(:)    , intent(inout) :: solution
    integer(ikind), dimension(:)    , intent(inout) :: arg
    write(*,*) 'Preconditioner Method Implementation'
  end subroutine method

end module PreconditionerMethodMod
