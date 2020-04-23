module MethodM

  use UtilitiesM
  use SparseKit

  use NonLinearSolversM

  implicit none

  private
  public :: MethodDT

  type, extends(NonLinearSolversDT) :: MethodDT
   contains
     procedure :: useSolver => method
  end type MethodDT

contains

  subroutine method(this, matrix, vector, solution, arg)
    implicit none
    class(MethodDT)           , intent(inout) :: this
    class(Sparse)               , intent(inout) :: matrix
    real(rkind)   , dimension(:), intent(inout) :: vector
    real(rkind)   , dimension(:), intent(inout) :: solution
    integer(ikind), dimension(:), intent(inout) :: arg
    write(*,*) 'Non Linear Solver Method Implementation'
  end subroutine method

end module MethodM
