module NonLinearSolverM

  use UtilitiesM
  use SparseKit

  use NonLinearSolversM

  implicit none

  private
  public :: NonLinearSolverDT, SetNonLinearSolver

  type NonLinearSolverDT
     class(NonLinearSolversDT), allocatable :: solver
   contains
     procedure :: init
     procedure :: changeSolver
     procedure :: solve
  end type NonLinearSolverDT

  interface SetNonLinearSolver
     procedure :: constructor
  end interface SetNonLinearSolver
  
contains

  type(NonLinearSolverDT) function constructor(solver)
    implicit none
    class(NonLinearSolversDT), intent(in) :: solver
    call constructor%init(solver)
  end function constructor

  subroutine init(this, solver)
    implicit none
    class(NonLinearSolverDT) , intent(inout) :: this
    class(NonLinearSolversDT), intent(in   ) :: solver
    allocate(this%solver, source = solver)
  end subroutine init

  subroutine changeSolver(this, newSolver)
    implicit none
    class(NonLinearSolverDT) , intent(inout) :: this
    class(NonLinearSolversDT), intent(inout) :: newSolver
    deallocate(this%solver)
    allocate(this%solver, source = newSolver)
  end subroutine changeSolver

  subroutine solve(this, matrix, vector, solution, arg)
    implicit none
    class(NonLinearSolverDT)    , intent(inout) :: this
    class(Sparse)               , intent(inout) :: matrix
    real(rkind)   , dimension(:), intent(inout) :: vector
    real(rkind)   , dimension(:), intent(inout) :: solution
    integer(ikind), dimension(:), intent(inout) :: arg
    call this%solver%useSolver(matrix, vector, solution, arg)
  end subroutine solve

end module NonLinearSolverM
