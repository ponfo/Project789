module MeshM
  use UtilitiesM
  use DebuggerM

  use NodeM
  use NodePtrM
  
  use ElementM
  use ElementPtrM
  
  use ConditionM
  use ConditionPtrM

  implicit none

  private
  public :: MeshDT, mesh

  type :: MeshDT
     integer(ikind)                                  :: id
     type(NodePtrDT)     , dimension(:), allocatable :: node
     type(ElementPtrDT)  , dimension(:), allocatable :: element
     type(ConditionPtrDT), dimension(:), allocatable :: condition
   contains
     procedure, public :: init

     procedure, public :: addNode
     procedure, public :: addElement
     procedure, public :: addCondition

     procedure, public :: removeNode
     procedure, public :: removeElement
     procedure, public :: removeCondition

     procedure, public :: free
  end type MeshDT

  interface mesh
     procedure :: constructor
  end interface mesh

contains

  type(MeshDT) function constructor(id, nNode, nElement, nCondition)
    implicit none
    integer(ikind), intent(in) :: id
    integer(ikind), intent(in) :: nNode
    integer(ikind), intent(in) :: nElement
    integer(ikind), intent(in) :: nCondition
    call contructor%init(id, nNode, nElement, nCondition)
  end function constructor

  subroutine init(this, id, nNode, nElement, nCondition)
    implicit none
    class(MeshDT) , intent(inout) :: this
    integer(ikind), intent(in)    :: id
    integer(ikind), intent(in)    :: nNode
    integer(ikind), intent(in)    :: nElement
    integer(ikind), intent(in)    :: nCondition
    this%id = id
    allocate(this%node(nNode))
    allocate(this%element(nElemen))
    allocate(this%condition(nCondition))
  end subroutine init

  subroutine addNode(this, id, node)
    implicit none
    class(MeshDT)         , intent(inout) :: this
    integer(ikind)        , intent(in)    :: id
    class(NodeDT) , target, intent(in)    :: node
    call this%node(id)%allocate(node)
  end subroutine addNode

  subroutine addElement(this, id, element)
    implicit none
    class(MeshDT)           , intent(inout) :: this
    integer(ikind)          , intent(in)    :: id
    class(ElementDT), target, intent(in)    :: element
    call this%element(id)%allocate(element)
  end subroutine addElement

  subroutine addCondition(this, id, condition)
    implicit none
    class(MeshDT)             , intent(inout) :: this
    integer(ikind)            , intent(in)    :: id
    class(ConditionDT), target, intent(in)    :: condition
    call this%condition(id)%allocate(condition)
  end subroutine addCondition

  subroutine removeNode(this, id)
    implicit none
    class(MeshDT) , intent(inout) :: this
    integer(ikind), intent(in)    :: id
    call this%node(id)%deallocate()
  end subroutine removeNode

  subroutine removeElement(this, id)
    implicit none
    class(MeshDT) , intent(inout) :: this
    integer(ikind), intent(in)    :: id
    call this%element(id)%deallocate()
  end subroutine removeElement

  subroutine removeCondition(this, id)
    implicit none
    class(MeshDT) , intent(inout) :: this
    integer(ikind), intent(in)    :: id
    call this%condition(id)%deallocate()
  end subroutine removeCondition

  final subroutine free(this)
    implicit none
    class(MeshDT), intent(inout) :: this
    deallocate(this%node)
    deallocate(this%element)
    deallocate(this%condition)
  end subroutine free

end module MeshM

  
  
