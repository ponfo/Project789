module NonLinearSolversMod

  use tools
  use sparseKit

  implicit none

  private
  public :: NonLinearSolversTYPE

  type, abstract :: NonLinearSolversTYPE
   contains
     procedure(NonLinearSolvers_procedure), deferred  :: useSolver
  end type NonLinearSolversTYPE

  abstract interface
     subroutine NonLinearSolvers_procedure(this, matrix, vector, solution, arg)
       import NonLinearSolversTYPE, sparse, rkind, ikind
       class(NonLinearSolversTYPE)            , intent(inout) :: this
       class(Sparse)                         , intent(inout) :: matrix
       real(rkind)            , dimension(:) , intent(inout) :: vector
       real(rkind)            , dimension(:) , intent(inout) :: solution
       integer(ikind)         , dimension(:) , intent(inout) :: arg
     end subroutine NonLinearSolvers_procedure
  end interface

end module NonLinearSolversMod
