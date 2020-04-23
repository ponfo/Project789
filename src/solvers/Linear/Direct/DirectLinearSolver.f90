module DirectLinearSolverMod

  use tools
  use sparsekit

  implicit none

  private
  public :: DirectLinearSolverTYPE

  type, abstract :: DirectLinearSolverTYPE
   contains
     procedure(DirectLinearSolver_procedure), deferred :: SolveSystem
  end type DirectLinearSolverTYPE

  abstract interface
     subroutine DirectLinearSolver_procedure(this, vector, matrix, solution, arg)
       import DirectLinearSolverTYPE, sparse, rkind, ikind
       class(DirectLinearSolverTYPE)           , intent(inout) :: this
       class(Sparse)                          , intent(inout)  :: matrix
       real(rkind)             , dimension(:) , intent(inout)  :: vector
       real(rkind)             , dimension(:) , intent(inout)  :: solution
       integer(ikind)          , dimension(:) , intent(inout)  :: arg
     end subroutine DirectLinearSolver_procedure
  end interface

end module DirectLinearSolverMod
