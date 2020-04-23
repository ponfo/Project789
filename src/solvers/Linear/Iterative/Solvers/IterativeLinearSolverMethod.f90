module IterativeLinearSolverMethodMod

  use tools
  use sparseKit
    
  use IterativeLinearSolverMod
  
  implicit none
  
  private
  public :: IterativeLinearSolverMethodTYPE
  
  type, extends(IterativeLinearSolverTYPE) :: IterativeLinearSolverMethodTYPE
   contains
     procedure :: solveSystem => method
  end type IterativeLinearSolverMethodTYPE
  
contains
  
  subroutine method(this, vector, matrix, solution, arg)
    implicit none
    class(IterativeLinearSolverMethodTYPE) , intent(inout)  :: this
    class(Sparse)                          , intent(inout)  :: matrix
    real(rkind)             , dimension(:) , intent(inout)  :: vector
    real(rkind)             , dimension(:) , intent(inout)  :: solution
    integer(ikind)          , dimension(:) , intent(inout)  :: arg
    write(*,*), 'Iterative Method Implementation'
    return
  end subroutine method
  
end module IterativeLinearSolverMethodMod
  
