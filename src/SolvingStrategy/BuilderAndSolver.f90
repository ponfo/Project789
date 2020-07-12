module BuilderAndSolverM

  use LinearSolverM
  use DirectLinearSolverM
  use IterativeLinearSolverM
  use NonLinearSolverM
  use NonLinearSolversM
  
  implicit none

  private
  public :: NewBuilderAndSolverDT, BuilderAndSolverDT, SetBuilderAndSolver
  
  type, abstract :: NewBuilderAndSolverDT
     type(LinearSolverDT)    :: userLinearSolver
     type(NonLinearSolverDT) :: userNonLinearSolver
   contains
     generic   :: setSolver => initNonLinearSolver&
          , initIterativeLinearSolver, initDirectLinearSolver
     procedure :: initIterativeLinearSolver
     procedure :: initDirectLinearSolver
     procedure :: initNonLinearSolver
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

  subroutine initIterativeLinearSolver(this, iterativeSolver)
    implicit none
    class(NewBuilderAndSolverDT  ), intent(inout) :: this
    class(IterativeLinearSolverDT), intent(in   ) :: iterativeSolver
    this%userLinearSolver = SetLinearSolver(iterativeSolver)
  end subroutine initIterativeLinearSolver
  
  subroutine initDirectLinearSolver(this, directSolver)
    implicit none
    class(NewBuilderAndSolverDT), intent(inout) :: this
    class(DirectLinearSolverDT ), intent(in   ) :: directSolver
    this%userLinearSolver = SetLinearSolver(directSolver)
  end subroutine initDirectLinearSolver
  
  subroutine initNonLinearSolver(this, NonLinearSolver)
    implicit none
    class(NewBuilderAndSolverDT), intent(inout) :: this
    class(NonLinearSolversDT   ), intent(in   ) :: nonLinearSolver
    this%userNonLinearSolver = SetNonLinearSolver(nonLinearSolver)
  end subroutine initNonLinearSolver
  
end module BuilderAndSolverM
