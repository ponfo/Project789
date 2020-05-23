module Structural3DApplicationM
  use UtilitiesM
  use DebuggerM

  use SourceM
  use NodeM
  
  use StructuralElementM
  use PressureM
  use StructuralMaterialM
  use StructuralModelM

  implicit none

  private
  public :: Structural3DApplicationDT, structural3DApplication

  type :: Structural3DApplicationDT
     type(NodeDT)              , dimension(:), allocatable :: node
     type(StructuralElementDT) , dimension(:), allocatable :: element
     type(PressureDT)          , dimension(:), allocatable :: pressure
     type(SourceDT)            , dimension(:), allocatable :: source
     type(StructuralMaterialDT), dimension(:), allocatable :: material
     type(StructuralModelDT)                               :: model
   contains
     procedure, public :: init
  end type Structural3DApplicationDT

  interface structural3DApplication
     procedure :: constructor
  end interface structural3DApplication

contains

  type(Structural3DApplicationDT) function  &
       constructor(nNode, nElement, nPressure, nSource, nMaterial, nGauss, nnz)
    implicit none
    integer(ikind), intent(in) :: nNode
    integer(ikind), intent(in) :: nElement
    integer(ikind), intent(in) :: nPressure
    integer(ikind), intent(in) :: nSource
    integer(ikind), intent(in) :: nMaterial
    integer(ikind), intent(in) :: nGauss
    integer(ikind), intent(in) :: nnz
    call constructor%init(nNode, nElement, nPressure, nSource, nMaterial, nGauss, nnz)
  end function constructor

  subroutine init(this, nNode, nElement, nPressure, nSource, nMaterial, nGauss, nnz)
    implicit none
    class(Structural3DApplicationDT), intent(inout) :: this
    integer(ikind)                  , intent(in)    :: nNode
    integer(ikind)                  , intent(in)    :: nElement
    integer(ikind)                  , intent(in)    :: nPressure
    integer(ikind)                  , intent(in)    :: nSource
    integer(ikind)                  , intent(in)    :: nMaterial
    integer(ikind)                  , intent(in)    :: nGauss
    integer(ikind)                  , intent(in)    :: nnz
    allocate(this%node(nNode))
    allocate(this%element(nElement))
    allocate(this%pressure(nPressure))
    allocate(this%source(nSource))
    allocate(this%material(nMaterial))
    call initGeometries(nGauss)
    this%model = structuralModel(                 &
           nDof = 3*nNode                         &
         , nnz = nnz                              &
         , id = 1                                 &
         , nNode = nNode                          &
         , nElement = nElement                    &
         , nCondition = nPressure                 )
  end subroutine init

end module Structural3DApplicationM
