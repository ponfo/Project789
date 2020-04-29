module ModelM

  use SparseKit
  
  use MeshM

  use SourceM

  use PropertyM

  implicit none

  private
  public :: ModelDT, model

  type ModelDT
     type(MeshDT), dimension(:), allocatable :: mesh
   contains
     procedure, public  :: initModel
     
     procedure, public  :: addMesh
     
     generic  , public  :: addNode         => addNodeOneMesh, addNodeMultiMesh
     generic  , public  :: addElement      => addElementOneMesh, addElementMultiMesh
     generic  , public  :: addCondition    => addConditionOneMesh, addConditionMultiMesh

     generic  , public  :: getnNode        => getnNodeOneMesh, getnNodeMultiMesh
     generic  , public  :: getnElement     => getnElementOneMesh, getnElementMultiMesh
     generic  , public  :: getnCondition   => getnConditionOneMesh, getnConditionMultiMesh
     generic  , public  :: getID           => getIDOneMesh, getIDMultiMesh
     generic  , public  :: getNode         => getNodeOneMesh, getNodeMultiMesh
     generic  , public  :: getElement      => getElementOneMesh, getElementMultiMesh
     generic  , public  :: getCondition    => getConditionOneMesh, getConditionMultiMesh

     generic  , public  :: removeNode      => removeNodeOneMesh, removeNodeMultiMesh
     generic  , public  :: removeElement   => removeElementOneMesh, removeElementMultiMesh
     generic  , public  :: removeCondition => removeConditionOneMesh, removeConditionMultiMesh

     procedure, private :: addNodeOneMesh
     procedure, private :: addNodeMultiMesh
     procedure, private :: addElementOneMesh
     procedure, private :: addElementMultiMesh
     procedure, private :: addConditionOneMesh
     procedure, private :: addConditionMultiMesh
     procedure, private :: getnNodeOneMesh
     procedure, private :: getnNodeMultiMesh
     procedure, private :: getnElementOneMesh
     procedure, private :: getnElementMultiMesh
     procedure, private :: getnConditionOneMesh
     procedure, private :: getnConditionMultiMesh
     procedure, private :: getIDOneMesh
     procedure, private :: getIDMultiMesh
     procedure, private :: getNodeOneMesh
     procedure, private :: getNodeMultiMesh
     procedure, private :: getElementOneMesh
     procedure, private :: getElementMultiMesh
     procedure, private :: getConditionOneMesh
     procedure, private :: getConditionMultiMesh
     procedure, private :: removeNodeOneMesh
     procedure, private :: removeNodeMultiMesh
     procedure, private :: removeElementOneMesh
     procedure, private :: removeElementMultiMesh
     procedure, private :: removeConditionOneMesh
     procedure, private :: removeConditionMultiMesh
   contains
     
  end type ModelDT

  interface model
     procedure :: constructor
  end interface model

contains

  type(ModelDT) function contructor(nMesh)
    implicit none
    integer(ikind), intent(in) :: nMesh
    call constructor%initModel(nMesh)
  end function contructor

  subroutine initModel(this, nMesh)
    implicit none
    class(ModelDT), intent(inout) :: this
    integer(ikind), intent(in)    :: nMesh
    allocate(this%mesh(nMesh))
  end subroutine initModel

  subroutine addMesh(this, id, nNode, nElement, nCondition)
    implicit none
    class(ModelDT), intent(inout) :: this
    integer(ikind), intent(in)    :: id
    integer(ikind), intent(in)    :: nNode
    integer(ikind), intent(in)    :: nElement
    integer(ikind), intent(in)    :: nCondition
    this%mesh(id) = mesh(id, nNode, nElement, nCondition)
  end subroutine addMesh

  subroutine addNodeOneMesh(this,  nodeID, node)
    implicit none
    class(ModelDT)        , intent(inout) :: this
    integer(ikind)        , intent(in)    :: nodeID
    class(NodeDT) , target, intent(in)    :: node
    call this%mesh(1)%addNode(nodeID, node)
  end subroutine addNodeOneMesh

  subroutine addNodeMultiMesh(this, meshID, nodeID, node)
    implicit none
    class(ModelDT)        , intent(inout) :: this
    integer(ikind)        , intent(in)    :: meshID
    integer(ikind)        , intent(in)    :: nodeID
    class(NodeDT) , target, intent(in)    :: node
    call this%mesh(meshID)%addNode(nodeID, node)
  end subroutine addNodeMultiMesh

  subroutine addElementOneMesh(this, elementID, element)
    implicit none
    class(ModelDT)          , intent(inout) :: this
    integer(ikind)          , intent(in)    :: elementID
    class(ElementDT), target, intent(in)    :: element
    call this%mesh(1)%addElement(elementID, element)
  end subroutine addElementOneMesh

  subroutine addElementMultiMesh(this, meshID, elementID, element)
    implicit none
    class(ModelDT)          , intent(inout) :: this
    integer(ikind)          , intent(in)    :: meshID
    integer(ikind)          , intent(in)    :: elementID
    class(ElementDT), target, intent(in)    :: element
    call this%mesh(meshID)%addElement(elementID, element)
  end subroutine addElementMultiMesh

  subroutine addConditionOneMesh(this, conditionID, condition)
    implicit none
    class(ModelDT)            , intent(inout) :: this
    integer(ikind)            , intent(in)    :: conditionID
    class(ConditionDT), target, intent(in)    :: condition
    call this%mesh(1)%addCondition(conditionID, condition)
  end subroutine addConditionOneMesh

  subroutine addConditionMultiMesh(this, meshID, conditionID, condition)
    implicit none
    class(ModelDT)            , intent(inout) :: this
    integer(ikind)            , intent(in)    :: meshID
    integer(ikind)            , intent(in)    :: conditionID
    class(ConditionDT), target, intent(in)    :: condition
    call this%mesh(meshID)%addCondition(conditionID, condition)
  end subroutine addConditionMultiMesh

  integer(ikind) function getnNodeOneMesh(this)
    implicit none
    class(ModelDT), intent(inout) :: this
    getnNodeOneMesh = this%mesh(1)%getnNode()
  end function getnNodeOneMesh

  integer(ikind) function getnNodeMultiMesh(this, meshID)
    implicit none
    class(ModelDT), intent(inout) :: this
    integer(ikind), intent(in)    :: meshID
    getnNodeMultiMesh = this%mesh(meshID)%getnNode()
  end function getnNodeMultiMesh
  
  integer(ikind) function getnElementOneMesh(this)
    implicit none
    class(ModelDT), intent(inout) :: this
    getnElementOneMesh = this%mesh(1)%getnElement()
  end function getnElementOneMesh

  integer(ikind) function getnElementMultiMesh(this, meshID)
    implicit none
    class(ModelDT), intent(inout) :: this
    integer(ikind), intent(in)    :: meshID
    getnElementMultiMesh = this%mesh(meshID)%getnElement()
  end function getnElementMultiMesh

  integer(ikind) function getnConditionOneMesh(this)
    implicit none
    class(ModelDT), intent(inout) :: this
    getnConditionOneMesh = this%mesh(1)%getnCondition()
  end function getnConditionOneMesh

  integer(ikind) function getnConditionMultiMesh(this, meshID)
    implicit none
    class(ModelDT), intent(inout) :: this
    integer(ikind), intent(in)    :: meshID
    getnConditionMultiMesh = this%mesh(meshID)%getnCondition()
  end function getnConditionMultiMesh

  integer(ikind) function getIDOneMesh(this)
    implicit none
    class(ModelDT), intent(inout) :: this
    getIDOneMesh = this%mesh(1)%getID()
  end function getIDOneMesh

  integer(ikind) function getIDMultiMesh(this, meshID)
    implicit none
    class(ModelDT), intent(inout) :: this
    integer(ikind), intent(in)    :: meshID
    getIDMultiMesh = this%mesh(meshID)%getID()
  end function getIDMultiMesh

  type(NodePtrDT) function getNodeOneMesh(this)
    implicit none
    class(ModelDT), intent(inout) :: this
    getNodeOneMesh = this%mesh(1)%getNode()
  end function getNodeOneMesh

  type(NodePtrDT) function getNodeMultiMesh(this, meshID)
    implicit none
    class(ModelDT), intent(inout) :: this
    integer(ikind), intent(in)    :: meshID
    getNodeMultiMesh = this%mesh(meshID)%getNode()
  end function getNodeMultiMesh

  type(ElementPtrDT) function getElementOneMesh(this)
    implicit none
    class(ModelDT), intent(inout) :: this
    getElementOneMesh = this%mesh(1)%getElement()
  end function getElementOneMesh

  type(ElementPtrDT) function getElementMultiMesh(this, meshID)
    implicit none
    class(ModelDT), intent(inout) :: this
    integer(ikind), intent(in)    :: meshID
    getElementMultiMesh = this%mesh(meshID)%getElement()
  end function getElementMultiMesh

  type(ConditionPtrDT) function getConditionOneMesh(this)
    implicit none
    class(ModelDT), intent(inout) :: this
    getConditionOneMesh = this%mesh(1)%getCondition()
  end function getConditionOneMesh

  type(ConditionPtrDT) function getConditionMultiMesh(this, meshID)
    implicit none
    class(ModelDT), intent(inout) :: this
    integer(ikind), intent(in)    :: meshID
    getConditionMultiMesh = this%mesh(meshID)%getCondition()
  end function getConditionMultiMesh

  subroutine removeNodeOneMesh(this, nodeID)
    implicit none
    class(ModelDT), intent(inout) :: this
    integer(ikind), intent(in)    :: nodeID
    call this%mesh(1)%removeNode(nodeID)
  end subroutine removeNodeOneMesh

  subroutine removeNodeMultiMesh(this, meshID, nodeID)
    implicit none
    class(ModelDT), intent(inout) :: this
    integer(ikind), intent(in)    :: meshID
    integer(ikind), intent(in)    :: nodeID
    call this%mesh(meshID)%removeNode(nodeID)
  end subroutine removeNodeMultiMesh

  subroutine removeElementOneMesh(this, elementID)
    implicit none
    class(ModelDT), intent(inout) :: this
    integer(ikind), intent(in)    :: elementID
    call this%mesh(1)%removeElement(elementID)
  end subroutine removeElementOneMesh

  subroutine removeElementMultiMesh(this, meshID, elementID)
    implicit none
    class(ModelDT), intent(inout) :: this
    integer(ikind), intent(in)    :: meshID
    integer(ikind), intent(in)    :: elementID
    call this%mesh(meshID)%removeElement(elementID)
  end subroutine removeElementMultiMesh

  subroutine removeConditionOneMesh(this, conditionID)
    implicit none
    class(ModelDT), intent(inout) :: this
    integer(ikind), intent(in)    :: conditionID
    call this%mesh(1)%removeCondition(conditionID)
  end subroutine removeConditionOneMesh

  subroutine removeConditionMultiMesh(this, meshID, conditionID)
    implicit none
    class(ModelDT), intent(inout) :: this
    integer(ikind), intent(in)    :: meshID
    integer(ikind), intent(in)    :: conditionID
    call this%mesh(meshID)%removeCondition(conditionID)
  end subroutine removeConditionMultiMesh
  
end module ModelM
