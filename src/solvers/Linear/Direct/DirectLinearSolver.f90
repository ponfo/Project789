module DirectLinearSolverM

  use UtilitiesM
  use Sparsekit

  implicit none

  private
  public :: DirectLinearSolverDT

  type, abstract :: DirectLinearSolverDT
   contains
     procedure(DirectLinearSolver_procedure), deferred :: SolveSystem
  end type DirectLinearSolverDT

  abstract interface
     subroutine DirectLinearSolver_procedure(this, vector, matrix, solution, arg)
       import DirectLinearSolverDT, sparse, rkind, ikind
       class(DirectLinearSolverDT)           , intent(inout) :: this
       class(Sparse)                         , intent(inout)  :: matrix
       real(rkind)             , dimension(:), intent(inout)  :: vector
       real(rkind)             , dimension(:), intent(inout)  :: solution
       integer(ikind)          , dimension(:), intent(inout)  :: arg
     end subroutine DirectLinearSolver_procedure
  end interface

end module DirectLinearSolverM
