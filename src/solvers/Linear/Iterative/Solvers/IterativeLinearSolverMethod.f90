module IterativeLinearSolverMethodM

  use UtilitiesM
  use SparseKit
    
  use IterativeLinearSolverM
  
  implicit none
  
  private
  public :: IterativeLinearSolverMethodDT
  
  type, extends(IterativeLinearSolverDT) :: IterativeLinearSolverMethodDT
   contains
     procedure :: solveSystem => method
  end type IterativeLinearSolverMethodDT
  
contains
  
  subroutine method(this, vector, matrix, solution, arg)
    implicit none
    class(IterativeLinearSolverMethodDT)  , intent(inout)  :: this
    class(Sparse)                         , intent(inout)  :: matrix
    real(rkind)             , dimension(:), intent(inout)  :: vector
    real(rkind)             , dimension(:), intent(inout)  :: solution
    integer(ikind)          , dimension(:), intent(inout)  :: arg
    write(*,*), 'Iterative Method Implementation'
    return
  end subroutine method
  
end module IterativeLinearSolverMethodM
  
