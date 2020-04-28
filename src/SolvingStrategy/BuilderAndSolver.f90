module BuilderAndSolverM

  implicit none

  private
  public :: NewBuilderAndSolverDT, BuilderAndSolverDT, SetBuilderAndSolver

  type, abstract :: NewBuilderAndSolverDT
  end type NewBuilderAndSolverDT

  interface SetBuilderAndSolver
     procedure :: constructor
  end interface SetBuilderAndSolver
  
  type BuilderAndSolverDT
     class(NewBuilderAndSolverDT), allocatable :: builderAndSolver
   contains
     procedure :: init
     procedure :: change
  end type BuilderAndSolverDT

contains

    type(BuilderAndSolverDT) function constructor(builderAndSolver)
    implicit none
    class(NewBuilderAndSolverDT), intent(in) :: builderAndSolver
    call constructor%init(builderAndSolver)
  end function constructor
  
  subroutine init(this, builderAndSolver)
    implicit none
    class(BuilderAndSolverDT)   , intent(inout) :: this
    class(NewBuilderAndSolverDT), intent(in)    :: builderAndSolver
    allocate(this%builderAndSolver, source = builderAndSolver)
  end subroutine init

  subroutine change(this, newBuilderAndSolver)
    implicit none
    class(BuilderAndSolverDT)   , intent(inout) :: this
    class(NewBuilderAndSolverDT), intent(in)    :: newBuilderAndSolver
    deallocate(this%builderAndSolver)
    allocate(this%builderAndSolver, source = newBuilderAndSolver)
  end subroutine change
  
end module BuilderAndSolverM
