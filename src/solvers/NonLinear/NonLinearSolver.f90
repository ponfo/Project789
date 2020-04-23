module NonLinearSolverMod

  use tools
  use sparseKit

  use NonLinearSolversMod

  implicit none

  private
  public :: NonLinearSolverTYPE, SetNonLinearSolver

  type NonLinearSolverTYPE
     class(NonLinearSolversTYPE), allocatable :: solver
   contains
     procedure :: init
     procedure :: changeSolver
     procedure :: solve
  end type NonLinearSolverTYPE

  interface SetNonLinearSolver
     procedure :: constructor
  end interface SetNonLinearSolver
  
contains

  type(NonLinearSolverTYPE) function constructor(solver)
    implicit none
    class(NonLinearSolversTYPE), intent(inout) :: solver
    call constructor%init(solver)
  end function constructor

  subroutine init(this, solver)
    implicit none
    class(NonLinearSolverTYPE) , intent(inout) :: this
    class(NonLinearSolversTYPE), intent(inout) :: solver
    allocate(this%solver, source = solver)
  end subroutine init

  subroutine changeSolver(this, newSolver)
    implicit none
    class(NonLinearSolverTYPE) , intent(inout) :: this
    class(NonLinearSolversTYPE), intent(inout) :: newSolver
    deallocate(this%solver)
    allocate(this%solver, source = newSolver)
  end subroutine changeSolver

  subroutine solve(this, matrix, vector, solution, arg)
    implicit none
    class(NonLinearSolverTYPE)  , intent(inout) :: this
    class(Sparse)               , intent(inout) :: matrix
    real(rkind)   , dimension(:), intent(inout) :: vector
    real(rkind)   , dimension(:), intent(inout) :: solution
    integer(ikind), dimension(:), intent(inout) :: arg
    call this%solver%useSolver(matrix, vector, solution, arg)
  end subroutine solve

end module NonLinearSolverMod
