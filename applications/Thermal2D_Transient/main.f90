program main
  
  use DataInputM
  use Thermal2DApplicationM
  use ThermalStrategyM
  use SolvingStrategyM

  implicit none

  type(Thermal2DApplicationDT)   :: application
  type(ThermalStrategyDT)        :: strategy
  type(SolvingStrategyDT)        :: solvingStrategy

  call initFEM2D(application)
  solvingStrategy = InitSolvingStrategy(strategy, application)
  call solvingStrategy%useStrategy()
  
end program main
