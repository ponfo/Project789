module Thermal3DApplicationM
  use UtilitiesM
  use DebuggerM

  use SourceM
  use NodeM
  
  use ThermalElementM
  use ConvectionOnSurfaceM
  use FluxOnSurfaceM
  use ThermalMaterialM
  use ThermalModelM

  implicit none

  private
  public :: Thermal3DApplicationDT, thermal3DApplication

  type :: Thermal3DApplicationDT
     type(NodeDT)               , dimension(:), allocatable :: node
     type(ThermalElementDT)     , dimension(:), allocatable :: element
     type(ConvectionOnSurfaceDT), dimension(:), allocatable :: convectionOS
     type(FluxOnSurfaceDT)      , dimension(:), allocatable :: normalFluxOS
     type(SourceDT)             , dimension(:), allocatable :: source
     type(ThermalMaterialDT)    , dimension(:), allocatable :: material
     type(ThermalModelDT)                                   :: model
   contains
     procedure, public :: init
  end type Thermal3DApplicationDT

  interface thermal3DApplication
     procedure :: constructor
  end interface thermal3DApplication

contains

  type(Thermal3DApplicationDT) function  &
       constructor(nNode, nElement, nConvection, nNormalFlux, nSource, nMaterial, nGauss)
    implicit none
    integer(ikind), intent(in) :: nNode
    integer(ikind), intent(in) :: nElement
    integer(ikind), intent(in) :: nConvection
    integer(ikind), intent(in) :: nNormalFlux
    integer(ikind), intent(in) :: nSource
    integer(ikind), intent(in) :: nMaterial
    integer(ikind), intent(in) :: nGauss
    call constructor%init(nNode, nElement, nConvection, nNormalFlux, nSource, nMaterial, nGauss)
  end function constructor

  subroutine init(this, nNode, nElement, nConvection, nNormalFlux, nSource, nMaterial, nGauss)
    implicit none
    class(Thermal3DApplicationDT), intent(inout) :: this
    integer(ikind)               , intent(in)    :: nNode
    integer(ikind)               , intent(in)    :: nElement
    integer(ikind)               , intent(in)    :: nConvection
    integer(ikind)               , intent(in)    :: nNormalFlux
    integer(ikind)               , intent(in)    :: nSource
    integer(ikind)               , intent(in)    :: nMaterial
    integer(ikind)               , intent(in)    :: nGauss
    allocate(this%node(nNode))
    allocate(this%element(nElement))
    allocate(this%convectionOS(nConvection))
    allocate(this%normalFluxOS(nNormalFlux))
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

end module Thermal3DApplicationM
