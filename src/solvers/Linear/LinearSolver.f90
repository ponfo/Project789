module LinearSolverM

  use UtilitiesM
  use SparseKit

  use DirectLinearSolverM
  use IterativeLinearSolverM

  use UseReorderSystemM

  implicit none

  private
  public :: LinearSolverDT, SetLinearSolver
  
  type, abstract :: LinearSolverDT
     class(DirectLinearSolverDT)   , allocatable :: directSolver
     class(IterativeLinearSolverDT), allocatable :: iterativeSolver
     type(UseReorderSystemDT)                    :: reorder
   contains
     generic   :: init            => iterativeInit, directInit
     procedure :: iterativeInit
     procedure :: directInit
     generic   :: changeSolver    => changeIterative, changeDirect
     procedure :: changeIterative
     procedure :: changeDirect
     procedure :: solve           => useSolver
  end type LinearSolverDT

  interface SetLinearSolver
     procedure :: iterativeConstructor
     procedure :: directConstructor
  end interface SetLinearSolver
  
contains

  type(LinearSolverDT) function iterativeConstructor(solver)
    implicit none
    class(IterativeLinearSolverDT), intent(in) :: solver
    call iterativeConstructor%init(solver)
  end function iterativeConstructor
  
  type(LinearSolverDT) function directConstructor(solver)
    implicit none
    class(DirectLinearSolverDT), intent(in) :: solver
    call directConstructor%init(solver)
  end function directConstructor

  subroutine iterativeInit(this, solver)
    implicit none
    class(LinearSolverDT)         , intent(inout) :: this
    class(IterativeLinearSolverDT), intent(in)    :: solver
    allocate(this%iterativeSolver, source = solver)
  end subroutine iterativeInit

  subroutine directInit(this, solver)
    implicit none
    class(LinearSolverDT)      , intent(inout) :: this
    class(DirectLinearSolverDT), intent(in)    :: solver
    allocate(this%directSolver, source = solver)
  end subroutine directInit
  
  subroutine changeIterative(this, newSolver)
    implicit none
    class(LinearSolverDT)         , intent(inout) :: this
    class(IterativeLinearSolverDT), intent(inout) :: newSolver
    if(allocated(this%directSolver)) then
       deallocate(this%directSolver)
    else if(allocated(this%iterativeSolver)) then
       deallocate(this%iterativeSolver)
    else
       write(*,*) '*** Solver Not Allocated ***'
       return
    end if
    allocate(this%iterativeSolver, source = newSolver)
  end subroutine changeIterative

  subroutine changeDirect(this, newSolver)
    implicit none
    class(LinearSolverDT)      , intent(inout) :: this
    class(DirectLinearSolverDT), intent(inout) :: newSolver
    if(allocated(this%iterativeSolver)) then
       deallocate(this%iterativeSolver)
    else if(allocated(this%directSolver)) then
       deallocate(this%directSolver)
       return
    end if
    allocate(this%directSolver, source = newSolver)
  end subroutine changeDirect
  
  subroutine useSolver(this, vector, matrix, solution, arg)
    implicit none
    class(LinearSolverDT)     , intent(inout) :: this
    class(Sparse)               , intent(inout) :: matrix
    real(rkind)   , dimension(:), intent(inout) :: vector
    real(rkind)   , dimension(:), intent(inout) :: solution
    integer(ikind), dimension(:), intent(inout) :: arg
    if(allocated(this%iterativeSolver)) then
       call this%iterativeSolver%solveSystem(vector, matrix, solution, arg)
    else if(allocated(this%directSolver)) then
       call this%directSolver%solveSystem(vector, matrix, solution, arg)
    else
       write(*,*) '*** Solver Not Allocated ***'
    end if
  end subroutine useSolver
  
end module LinearSolverM
