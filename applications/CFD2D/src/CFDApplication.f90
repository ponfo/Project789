module CFDApplicationM
  use UtilitiesM
  use DebuggerM

  use SourceM
  use NodeM
  
  use CFDElementM
  use 
  use CFDModelM

  implicit none

  private
  public :: CFDApplicationDT, cfdApplication

  type :: CFDApplicationDT
     type(NodeDT)       , dimension(:), allocatable :: node
     type(CFDElementDT) , dimension(:), allocatable :: element
     type
     type(SourceDT)     , dimension(:), allocatable :: source
     type(CFDModelDT)                               :: model
   contains
     procedure, public :: init
  end type CFDApplicationDT

  interface cfdApplication
     procedure :: constructor
  end interface cfdApplication

contains

  type(CFDApplicationDT) function  &
       constructor(nNode, nElement, nPressure, nSource, nGauss)
    implicit none
    integer(ikind), intent(in) :: nNode
    integer(ikind), intent(in) :: nElement
    integer(ikind), intent(in) :: nPressure
    integer(ikind), intent(in) :: nSource
    integer(ikind), intent(in) :: nGauss
    call constructor%init(nNode, nElement, nPressure, nSource, nGauss)
  end function constructor

  subroutine init(this, nNode, nElement, nPressure, nSource, nGauss)
    implicit none
    class(CFDApplicationDT), intent(inout) :: this
    integer(ikind)                  , intent(in)    :: nNode
    integer(ikind)                  , intent(in)    :: nElement
    integer(ikind)                  , intent(in)    :: nPressure
    integer(ikind)                  , intent(in)    :: nSource
    integer(ikind)                  , intent(in)    :: nGauss
    allocate(this%node(nNode))
    allocate(this%element(nElement))
    allocate(this%pressure(nPressure))
    allocate(this%source(nSource))
    call initGeometries(nGauss)
    this%model = cfdModel(                &
           nDof = 2*nNode                &
         , nnz = nElement*256            &
         , id = 1                        &
         , nNode = nNode                 &
         , nElement = nElement           &
         , nCondition =                  )
  end subroutine init

end module CFDApplicationM
