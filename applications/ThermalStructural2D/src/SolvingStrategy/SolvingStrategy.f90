module SolvingStrategyM

  use NewSolvingStrategyM
  use NewStrategyM
  use ModelM
  use ThermalModelM
  use StructuralModelM

  implicit none

  private
  public :: SolvingStrategyDT, InitSolvingStrategy

  type, extends(NewSolvingStrategyDT) :: SolvingStrategyDT
     class(ThermalModelDT)   , pointer :: thermalModel
     class(structuralModelDT), pointer :: structuralModel
   contains
  end type SolvingStrategyDT

  interface InitSolvingStrategy
     procedure :: structuralConstructor
     procedure :: thermalConstructor
  end interface InitSolvingStrategy
 
contains

  type(SolvingStrategyDT) function thermalConstructor(strategy, model)
    class(NewStrategyDT)         , intent(in) :: strategy
    class(ThermalModelDT), target, intent(in) :: model
    thermalConstructor%thermalModel => model
    call thermalConstructor%setStrategy(strategy)
  end function thermalConstructor

  type(SolvingStrategyDT) function structuralConstructor(strategy, model)
    class(NewStrategyDT)         , intent(in) :: strategy
    class(StructuralModelDT), target, intent(in) :: model
    structuralConstructor%structuralModel => model
    call structuralConstructor%setStrategy(strategy)
  end function structuralConstructor
  
end module SolvingStrategyM
