module ThermalStruct2DApplicationM
  use UtilitiesM
  use DebuggerM

  use SourceM
  use NodeM

  use ThermalStructuralNodeM
  use ThermalElementM
  use ThermalStructuralElementM
  use ConvectionOnLineM
  use FluxOnLineM
  use PressureM
  use ThermalStructuralMaterialM
  use ThermalModelM
  use StructuralModelM

  implicit none

  private
  public :: ThermalStruct2DApplicationDT, thermalStruct2DApplication

  type :: ThermalStruct2DApplicationDT
     type(ThermalStructuralNodeDT)    , dimension(:), allocatable :: node
     type(ThermalElementDT)           , dimension(:), allocatable :: thermalElement
     type(ThermalStructuralElementDT) , dimension(:), allocatable :: structuralElement
     type(ConvectionOnLineDT)         , dimension(:), allocatable :: convectionOL
     type(FluxOnLineDT)               , dimension(:), allocatable :: normalFluxOL
     type(PressureDT)                 , dimension(:), allocatable :: pressure
     type(SourceDT)                   , dimension(:), allocatable :: source
     type(ThermalStructuralMaterialDT), dimension(:), allocatable :: material
     type(ThermalModelDT)                                         :: model
   contains
     procedure, public :: init
  end type ThermalStruct2DApplicationDT

  interface thermalStruct2DApplication
     procedure :: constructor
  end interface thermalStruct2DApplication

contains

  type(ThermalStruct2DApplicationDT) function  &
       constructor(nNode, nElement, nConvection, nNormalFlux, nPressure, nSource, nMaterial, nGauss)
    implicit none
    integer(ikind), intent(in) :: nNode
    integer(ikind), intent(in) :: nElement
    integer(ikind), intent(in) :: nConvection
    integer(ikind), intent(in) :: nNormalFlux
    integer(ikind), intent(in) :: nPressure
    integer(ikind), intent(in) :: nSource
    integer(ikind), intent(in) :: nMaterial
    integer(ikind), intent(in) :: nGauss
    call constructor%init &
         (nNode, nElement, nConvection, nNormalFlux, nPressure, nSource, nMaterial, nGauss)
  end function constructor

  subroutine init &
       (this, nNode, nElement, nConvection, nNormalFlux, nPressure, nSource, nMaterial, nGauss)
    implicit none
    class(ThermalStruct2DApplicationDT), intent(inout) :: this
    integer(ikind)                     , intent(in)    :: nNode
    integer(ikind)                     , intent(in)    :: nElement
    integer(ikind)                     , intent(in)    :: nConvection
    integer(ikind)                     , intent(in)    :: nNormalFlux
    integer(ikind)                     , intent(in)    :: nPressure
    integer(ikind)                     , intent(in)    :: nSource
    integer(ikind)                     , intent(in)    :: nMaterial
    integer(ikind)                     , intent(in)    :: nGauss
    allocate(this%node(nNode))
    allocate(this%thermalElement(nElement))
    allocate(this%structuralElement(nElement))
    allocate(this%convectionOL(nConvection))
    allocate(this%normalFluxOL(nNormalFlux))
    allocate(this%pressure(nPressure))
    allocate(this%source(nSource))
    allocate(this%material(nMaterial))
    call initGeometries(nGauss)
    this%model = thermalModel(                    &
           nDof = nNode                           &
         , nnz = nElement*64                      &
         , id = 1                                 &
         , nNode = nNode                          &
         , nElement = nElement                    &
         , nCondition = nConvection + nNormalFlux )
  end subroutine init

end module ThermalStruct2DApplicationM
