module MethodMod

  use tools
  use sparseKit

  use NonLinearSolversMod

  implicit none

  private
  public :: MethodTYPE

  type, extends(NonLinearSolversTYPE) :: MethodTYPE
   contains
     procedure :: useSolver => method
  end type MethodTYPE

contains

  subroutine method(this, matrix, vector, solution, arg)
    implicit none
    class(MethodTYPE)           , intent(inout) :: this
    class(Sparse)               , intent(inout) :: matrix
    real(rkind)   , dimension(:), intent(inout) :: vector
    real(rkind)   , dimension(:), intent(inout) :: solution
    integer(ikind), dimension(:), intent(inout) :: arg
    write(*,*) 'Non Linear Solver Method Implementation'
  end subroutine method

end module MethodMod
