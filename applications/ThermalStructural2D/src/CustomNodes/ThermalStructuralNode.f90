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

  type(ThermalStructuralNodeDT) function constructor2D(id, nDof, nSource, x, y, stableTemp)
    implicit none
    integer(ikind), intent(in) :: id
    integer(ikind), intent(in) :: nDof
    integer(ikind), intent(in) :: nSource
    real(rkind)   , intent(in) :: x
    real(rkind)   , intent(in) :: y
    real(rkind)   , intent(in) :: stableTemp
    call constructor2D%initThermalStructuralNode2D(id, nDof, nSource, x, y, stableTemp)
  end function constructor2D

  subroutine initThermalStructuralNode2D(this, id, nDof, nSource, x, y, stableTemp)
    implicit none
    class(ThermalStructuralNodeDT), intent(inout) :: this
    integer(ikind)                , intent(in)    :: id
    integer(ikind)                , intent(in)    :: nDof
    integer(ikind)                , intent(in)    :: nSource
    real(rkind)                   , intent(in)    :: x
    real(rkind)                   , intent(in)    :: y
    real(rkind)                   , intent(in)    :: stableTemp
    call this%initNode2DMultiSource(id, nDof, nSource, x, y)

    !PARA QUE ESTO?
    
