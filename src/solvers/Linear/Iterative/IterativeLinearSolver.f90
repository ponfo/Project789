module IterativeLinearSolverMod

  use tools
  use sparsekit

  use UsePreconditionerMod
  
  implicit none

  private
  public :: IterativeLinearSolverTYPE

  type, abstract :: IterativeLinearSolverTYPE
     type(UsePreconditionerTYPE)                            :: preconditioner
   contains
     procedure(IterativeLinearSolver_procedure)  , deferred :: SolveSystem
  end type IterativeLinearSolverTYPE

  abstract interface
     subroutine IterativeLinearSolver_procedure(this, vector, matrix, solution, arg)
       import IterativeLinearSolverTYPE, sparse, rkind, ikind
       class(IterativeLinearSolverTYPE)       , intent(inout) :: this
       class(Sparse)                          , intent(inout) :: matrix
       real(rkind)             , dimension(:) , intent(inout) :: vector
       real(rkind)             , dimension(:) , intent(inout) :: solution
       integer(ikind)          , dimension(:) , intent(inout) :: arg
     end subroutine IterativeLinearSolver_procedure
  end interface

end module IterativeLinearSolverMod
