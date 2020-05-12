program main
  
  use DataInputM
  use Thermal2DApplicationM
  use ThermalStrategyM

  implicit none

  type(Thermal2DApplicationDT)   :: application
  type(ThermalStrategyDT)        :: thermalStrategy

  call initFEM2D(application)
  call thermalStrategy%buildStrategyAndSolve(application%model)
  
end program main
