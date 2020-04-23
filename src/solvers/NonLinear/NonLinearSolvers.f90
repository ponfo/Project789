module NonLinearSolversM

  use UtilitiesM
  use SparseKit

  implicit none

  private
  public :: NonLinearSolversDT

  type, abstract :: NonLinearSolversDT
   contains
     procedure(NonLinearSolvers_procedure), deferred  :: useSolver
  end type NonLinearSolversDT

  abstract interface
     subroutine NonLinearSolvers_procedure(this, matrix, vector, solution, arg)
       import NonLinearSolversDT, sparse, rkind, ikind
       class(NonLinearSolversDT)            , intent(inout) :: this
       class(Sparse)                        , intent(inout) :: matrix
       real(rkind)            , dimension(:), intent(inout) :: vector
       real(rkind)            , dimension(:), intent(inout) :: solution
       integer(ikind)         , dimension(:), intent(inout) :: arg
     end subroutine NonLinearSolvers_procedure
  end interface

end module NonLinearSolversM
