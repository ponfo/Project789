module SolvingStrategyM

  use NewSolvingStrategyM
  use NewStrategyM
  use SchemeM
  use BuilderAndSolverM
  use StructuralModelM

  implicit none

  private
  public :: SolvingStrategyDT, InitSolvingStrategy

  type, extends(NewSolvingStrategyDT) :: SolvingStrategyDT
     type(BuilderAndSolverDT)          :: builderAndSolver
     type(SchemeDT)                    :: scheme
     class(StructuralModelDT), pointer :: structuralModel
   contains
  end type SolvingStrategyDT

  interface InitSolvingStrategy
     procedure :: constructor
  end interface InitSolvingStrategy

contains

  type(SolvingStrategyDT) function constructor(strategy, model)
    class(NewStrategyDT)            , intent(in) :: strategy
    class(StructuralModelDT), target, intent(in) :: model
    constructor%structuralModel => model
    call constructor%setStrategy(strategy)
  end function constructor
  
end module SolvingStrategyM
