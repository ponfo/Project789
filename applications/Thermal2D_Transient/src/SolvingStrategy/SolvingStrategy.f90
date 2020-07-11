module SolvingStrategyM

  use NewSolvingStrategyM
  use NewStrategyM
  use Thermal2DApplicationM

  implicit none

  private
  public :: SolvingStrategyDT, InitSolvingStrategy

  type, extends(NewSolvingStrategyDT) :: SolvingStrategyDT
     class(thermal2DApplicationDT), pointer :: application
   contains
  end type SolvingStrategyDT

  interface InitSolvingStrategy
     procedure :: constructor
  end interface InitSolvingStrategy

contains

  type(SolvingStrategyDT) function constructor(strategy, application)
    class(NewStrategyDT)                 , intent(in) :: strategy
    class(Thermal2DApplicationDT), target, intent(in) :: application
    constructor%application => application
    call constructor%setStrategy(strategy)
  end function constructor
  
end module SolvingStrategyM
