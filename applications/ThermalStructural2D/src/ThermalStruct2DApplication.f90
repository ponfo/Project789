module ThermalStruct2DApplicationM
  use UtilitiesM
  use DebuggerM

  use SourceM
  use NodeM

  use ThermalElementM
  use ThermalStructuralElementM
  use ConvectionOnLineM
  use FluxOnLineM
  use PressureM
  use StructuralMaterialM
  use ThermalModelM
  use StructuralModelM

  implicit none

  private
  public :: ThermalStruct2DApplicationDT, thermalStruct2DApplication

  type :: ThermalStruct2DApplicationDT
     integer(ikind)                                               :: nNode
     integer(ikind)                                               :: nElement
     type(NodeDT)                     , dimension(:), allocatable :: node
     type(ThermalElementDT)           , dimension(:), allocatable :: thermalElement
     type(ThermalStructuralElementDT) , dimension(:), allocatable :: structuralElement
     type(ConvectionOnLineDT)         , dimension(:), allocatable :: convectionOL
     type(FluxOnLineDT)               , dimension(:), allocatable :: normalFluxOL
     type(PressureDT)                 , dimension(:), allocatable :: pressure
     type(SourceDT)                   , dimension(:), allocatable :: source
     type(StructuralMaterialDT)       , dimension(:), allocatable :: material
     type(ThermalModelDT)                                         :: thermalModel
     type(StructuralModelDT)                                      :: structuralModel
   contains
     procedure, public :: init
     procedure, public :: transitionToStructural
  end type ThermalStruct2DApplicationDT

  interface thermalStruct2DApplication
     procedure :: constructor
  end interface thermalStruct2DApplication

contains

  type(ThermalStruct2DApplicationDT) function constructor &
       (nNode, nElement, nConvection, nNormalFlux, nPressure, nSource, nMaterial, nGauss)
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
         (nNode, nElement, nConvection, nNormalFlux, nPressure, nSource,  nMaterial, nGauss)
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
    this%nNode = nNode
    this%nElement = nElement
    allocate(this%node(nNode))
    allocate(this%thermalElement(nElement))
    allocate(this%structuralElement(this%nElement))
    allocate(this%convectionOL(nConvection))
    allocate(this%normalFluxOL(nNormalFlux))
    allocate(this%pressure(nPressure))
    allocate(this%source(nSource))
    allocate(this%material(nMaterial))
    call initGeometries(nGauss)
    call initGeometriesTS(nGauss)
    this%thermalModel = thermalModel(             &
           nDof = nNode                           &
         , nnz = nElement*64                      &
         , id = 1                                 &
         , nNode = nNode                          &
         , nElement = nElement                    &
         , nCondition = nConvection + nNormalFlux )
    call this%structuralModel%initWithoutSystem(   &
           id = 1                                 &
         , nNode = nNode                          &
         , nElement = nElement                    &
         , nCondition = nPressure                 )
  end subroutine init

  subroutine transitionToStructural(this)
    implicit none
    class(ThermalStruct2DApplicationDT), intent(inout) :: this
    integer(ikind)                                     :: i
    deallocate(this%thermalElement)
    deallocate(this%convectionOL)
    deallocate(this%normalFluxOL)
    call this%thermalModel%freeSystem()
    call this%structuralModel%initSystem(2*this%nNode, this%nElement*256)
    do i = 1, this%nNode
       call this%node(i)%assignDof(2, this%structuralModel%dof(i*2-1))
       call this%node(i)%assignDof(3, this%structuralModel%dof(i*2))
       call this%node(i)%assignName(2, this%structuralModel%displxDofName)
       call this%node(i)%assignName(3, this%structuralModel%displyDofName)
    end do
  end subroutine transitionToStructural

end module ThermalStruct2DApplicationM
