module SolvingStrategyM

  use NewSolvingStrategyM
  use NewStrategyM
  use SchemeM
  use BuilderAndSolverM
  use ThermalModelM

  implicit none

  private
  public :: SolvingStrategyDT, InitSolvingStrategy

  type, extends(NewSolvingStrategyDT) :: SolvingStrategyDT
     type(BuilderAndSolverDT)       :: builderAndSolver
     type(SchemeDT)                 :: scheme
     class(ThermalModelDT), pointer :: thermalModel
   contains
  end type SolvingStrategyDT

  interface InitSolvingStrategy
     procedure :: constructor
  end interface InitSolvingStrategy

contains

  type(SolvingStrategyDT) function constructor(strategy, model)
    class(NewStrategyDT)         , intent(in) :: strategy
    class(ThermalModelDT), target, intent(in) :: model
    constructor%thermalModel => model
    call constructor%setStrategy(strategy)
  end function constructor
  
end module SolvingStrategyM
