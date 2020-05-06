module ThermalStructuralNodeM
  use UtilitiesM

  implicit none

  private
  public :: ThermalStructuralNodeDT, thermalStructuralNode

  type, extends(NodeDT) :: ThermalStructuralNodeDT
     type(TemperatureLoadDT), pointer :: temperatureLoad
   contains
     procedure, public :: initThermalStructuralNode2D
  end type ThermalStructuralNodeDT

  interface thermalStructuralNode
     procedure :: contructor2D
  end interface thermalStructuralNode

contains

  type(ThermalStructuralNodeDT) function constructor2D(id, nDof, x, y, stableTemp, temperature)
