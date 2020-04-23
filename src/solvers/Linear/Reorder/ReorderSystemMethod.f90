module ReorderSystemMethodM

  use UtilitiesM
  use SparseKit

  use ReorderSystemM

  implicit none

  private
  public :: ReorderSystemMethodDT

  type, extends(ReorderSystemDT) :: ReorderSystemMethodDT
   contains
     procedure :: useReorder => method
  end type ReorderSystemMethodDT

contains

  subroutine method(this, vector, matrix, solution, arg)
    implicit none
    class(ReorderSystemMethodDT), intent(inout) :: this
    class(Sparse)                   , intent(inout) :: matrix
    real(rkind)   , dimension(:)    , intent(inout) :: vector
    real(rkind)   , dimension(:)    , intent(inout) :: solution
    integer(ikind), dimension(:)    , intent(inout) :: arg
    write(*,*) 'Reorder Method implementation'
  end subroutine method

end module ReorderSystemMethodM
