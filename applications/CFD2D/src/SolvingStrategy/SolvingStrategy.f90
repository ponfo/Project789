module SolvingStrategyM

  use NewSolvingStrategyM
  use NewStrategyM
  use SchemeM
  use BuilderAndSolverM
  use CFDApplicationM

  implicit none

  private
  public :: SolvingStrategyDT, InitSolvingStrategy

  type, extends(NewSolvingStrategyDT) :: SolvingStrategyDT
     type(BuilderAndSolverDT)         :: builderAndSolver
     type(SchemeDT)                   :: scheme
     class(CFDApplicationDT), pointer :: application
   contains
  end type SolvingStrategyDT

  interface InitSolvingStrategy
     procedure :: constructor
  end interface InitSolvingStrategy 

contains

  type(SolvingStrategyDT) function constructor(strategy, application)
    class(NewStrategyDT)           , intent(in) :: strategy
    class(CFDApplicationDT), target, intent(in) :: application
    constructor%application => application
    call constructor%setStrategy(strategy)
  end function constructor
  
end module SolvingStrategyM
