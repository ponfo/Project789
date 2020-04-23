module IterativeLinearSolverM

  use UtilitiesM
  use Sparsekit

  use UsePreconditionerM
  
  implicit none

  private
  public :: IterativeLinearSolverDT

  type, abstract :: IterativeLinearSolverDT
     type(UsePreconditionerDT)                            :: preconditioner
   contains
     procedure(IterativeLinearSolver_procedure), deferred :: SolveSystem
  end type IterativeLinearSolverDT

  abstract interface
     subroutine IterativeLinearSolver_procedure(this, vector, matrix, solution, arg)
       import IterativeLinearSolverDT, sparse, rkind, ikind
       class(IterativeLinearSolverDT)        , intent(inout) :: this
       class(Sparse)                         , intent(inout) :: matrix
       real(rkind)             , dimension(:), intent(inout) :: vector
       real(rkind)             , dimension(:), intent(inout) :: solution
       integer(ikind)          , dimension(:), intent(inout) :: arg
     end subroutine IterativeLinearSolver_procedure
  end interface

end module IterativeLinearSolverM
