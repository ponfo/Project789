module BuilderAndSolverM

  implicit none

  private
  public :: NewBuilderAndSolverDT, BuilderAndSolverDT, SetBuilderAndSolver

  type, abstract :: NewBuilderAndSolverDT
   contains
     procedure(BuilderAndSolver_procedure), deferred :: useBAndS
  end type NewBuilderAndSolverDT

  abstract interface
     subroutine BuilderAndSolver_Procedure(this)
       import NewBuilderAndSolverDT
       class(NewBuilderAndSolverDT), intent(inout) :: this
     end subroutine BuilderAndSolver_Procedure
  end interface

  interface SetBuilderAndSolver
     procedure :: constructor
  end interface SetBuilderAndSolver
  
  type BuilderAndSolverDT
     class(NewBuilderAndSolverDT), allocatable :: builderAndSolver
   contains
     procedure :: init
     procedure :: change
     procedure :: use
!!$     procedure :: buildAndSolve
!!$     procedure :: buildLHS
!!$     procedure :: buildRHS
!!$     procedure :: applyDirichletConditions
!!$     procedure :: systemSolve
     !procedure :: CalculateReactions
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

  subroutine use(this)
    implicit none
    class(BuilderAndSolverDT), intent(inout) :: this
    call this%builderAndSolver%useBandS()
  end subroutine use
  
end module BuilderAndSolverM
