module SolvingStrategyM

  use ProcessM

  use BuilderAndSolverM
  use SchemeM

  implicit none

  private
  public :: NewSolvingStrategyDT, SolvingStrategyDT, SetSolvingStrategy

  type, abstract, extends(ProcessDT) :: NewSolvingStrategyDT
     class(BuilderAndSolverDT), allocatable :: builderAndSolver
     class(SchemeDT)          , allocatable :: scheme
  end type NewSolvingStrategyDT

  interface SetSolvingStrategy
     procedure :: constructor
  end interface SetSolvingStrategy
  
  type SolvingStrategyDT
     class(NewSolvingStrategyDT) , allocatable :: strategy
   contains
     procedure :: init
     procedure :: change
  end type SolvingStrategyDT

contains

  type(SolvingStrategyDT) function constructor(strategy)
    implicit none
    class(NewSolvingStrategyDT) , intent(in) :: strategy
    call constructor%init(strategy)
  end function constructor

  subroutine init(this, strategy)
    implicit none
    class(SolvingStrategyDT)    , intent(inout) :: this
    class(NewSolvingStrategyDT) , intent(in)    :: strategy
    allocate(this%strategy        , source = strategy)
  end subroutine init
  
  subroutine change(this, newStrategy)
    implicit none
    class(SolvingStrategyDT)   , intent(inout) :: this
    class(NewSolvingStrategyDT), intent(in)    :: newStrategy
    deallocate(this%strategy)
    allocate(this%strategy, source = newStrategy)
  end subroutine change
  
end module SolvingStrategyM
