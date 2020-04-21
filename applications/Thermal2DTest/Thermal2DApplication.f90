module Thermal2DApplicationM
  use ThermalElementM
  use ThermalConditionM

  implicit none

  private
  public :: Thermal2DApplicationDT, thermal2DApplication

  type :: Thermal2DApplicationDT
     type(NodeDT)            , dimension(:), allocatable :: node
     type(ThermalElementDT)  , dimension(:), allocatable :: element
     type(ThermalConditionDT), dimension(:), allocatable :: condition
     type(SourceFuncDT)      , dimension(:), allocatable :: sourceFunc
     type(ThermalSourceDT)   , dimension(:), allocatable :: source
     type(ThermalMaterial)   , dimension(:), allocatable :: material
     type(MeshM)                                         :: mesh
     type(ModelM)                                        :: model
   contains
     procedure, public :: assemble
     procedure, public :: solve
  end type Thermal2DApplicationDT

  interface thermal2DApplication
     procedure :: constructor
  end interface thermal2DApplication

contains

  type(Thermal2DApplicationDT) function  &
       constructor(nNode, nElement, nCondition, nSource, nMaterial, nGauss)
    implicit none
    integer(ikind), intent(in) :: nNode
    integer(ikind), intent(in) :: nElement
    integer(ikind), intent(in) :: nCondition
    integer(ikind), intent(in) :: nSource
    integer(ikind), intent(in) :: nMaterial
    integer(ikind), intent(in) :: nGauss
    call this%init(nNode, nElement, nCondition, nSource, nMaterial, nGauss)
  end function constructor

  subroutine init(this, nNode, nElement, nCondition, nSource, nMaterial, nGauss)
    implicit none
    class(Thermal2DApplicationDT), intent(inout) :: this
    integer(ikind)               , intent(in)    :: nNode
    integer(ikind)               , intent(in)    :: nElement
    integer(ikind)               , intent(in)    :: nCondition
    integer(ikind)               , intent(in)    :: nSource
    integer(ikind)               , intent(in)    :: nMaterial
    integer(ikind)               , intent(in)    :: nGauss
    allocate(this%node(nNode))
    allocate(this%element(nElement))
    allocate(this%condition(nCondition))
    allocate(this%source(nSource))
    allocate(this%material(nMaterial))
    call initGeometries(nGauss)
  end subroutine init

end module Thermal2DApplicationM
