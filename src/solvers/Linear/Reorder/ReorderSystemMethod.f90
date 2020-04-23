module ReorderSystemMethodMod

  use tools
  use sparseKit

  use ReorderSystemMod

  implicit none

  private
  public :: ReorderSystemMethodTYPE

  type, extends(ReorderSystemTYPE) :: ReorderSystemMethodTYPE
   contains
     procedure :: useReorder => method
  end type ReorderSystemMethodTYPE

contains

  subroutine method(this, vector, matrix, solution, arg)
    implicit none
    class(ReorderSystemMethodTYPE), intent(inout) :: this
    class(Sparse)                   , intent(inout) :: matrix
    real(rkind)   , dimension(:)    , intent(inout) :: vector
    real(rkind)   , dimension(:)    , intent(inout) :: solution
    integer(ikind), dimension(:)    , intent(inout) :: arg
    write(*,*) 'Reorder Method implementation'
  end subroutine method

end module ReorderSystemMethodMod
