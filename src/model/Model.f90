module ModelM
  use UtilitiesM

  use SparseKit

  use NodeM
  use NodePtrM

  use ElementM
  use ElementPtrM

  use ConditionM
  use ConditionPtrM
  
  use MeshM
  use ProcessInfoM

  use SourceM

  use PropertyM

  implicit none

  private
  public :: ModelDT, model

  type ModelDT
     type(MeshDT)       , dimension(:), allocatable :: mesh
     type(ProcessInfoDT)                            :: processInfo
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
     generic  , public  :: getNodes        => getNodesOneMesh, getNodesMultiMesh
     generic  , public  :: getElements     => getElementsOneMesh, getElementsMultiMesh
     generic  , public  :: getConditions   => getConditionsOneMesh, getConditionsMultiMesh

     generic  , public  :: removeNode      => removeNodeOneMesh, removeNodeMultiMesh
     generic  , public  :: removeElement   => removeElementOneMesh, removeElementMultiMesh
     generic  , public  :: removeCondition => removeConditionOneMesh, removeConditionMultiMesh

     procedure, public  :: getMesh
     procedure, public  :: getProcessInfo

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
     procedure, private :: getNodesOneMesh
     procedure, private :: getNodesMultiMesh
     procedure, private :: getElementsOneMesh
     procedure, private :: getElementsMultiMesh
     procedure, private :: getConditionsOneMesh
     procedure, private :: getConditionsMultiMesh
     procedure, private :: removeNodeOneMesh
     procedure, private :: removeNodeMultiMesh
     procedure, private :: removeElementOneMesh
     procedure, private :: removeElementMultiMesh
     procedure, private :: removeConditionOneMesh
     procedure, private :: removeConditionMultiMesh
  end type ModelDT

  interface model
     procedure :: constructor
  end interface model

contains

  type(ModelDT) function constructor(nMesh)
    implicit none
    integer(ikind), intent(in) :: nMesh
    call constructor%initModel(nMesh)
  end function constructor

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

  subroutine addNodeOneMesh(this, nodeID, node)
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

  integer(ikind) pure function getnNodeOneMesh(this)
    implicit none
    class(ModelDT), intent(in) :: this
    getnNodeOneMesh = this%mesh(1)%getnNode()
  end function getnNodeOneMesh

  integer(ikind) pure function getnNodeMultiMesh(this, meshID)
    implicit none
    class(ModelDT), intent(in) :: this
    integer(ikind), intent(in) :: meshID
    getnNodeMultiMesh = this%mesh(meshID)%getnNode()
  end function getnNodeMultiMesh
  
  integer(ikind) pure function getnElementOneMesh(this)
    implicit none
    class(ModelDT), intent(in) :: this
    getnElementOneMesh = this%mesh(1)%getnElement()
  end function getnElementOneMesh

  integer(ikind) pure function getnElementMultiMesh(this, meshID)
    implicit none
    class(ModelDT), intent(in) :: this
    integer(ikind), intent(in) :: meshID
    getnElementMultiMesh = this%mesh(meshID)%getnElement()
  end function getnElementMultiMesh

  integer(ikind) pure function getnConditionOneMesh(this)
    implicit none
    class(ModelDT), intent(in) :: this
    getnConditionOneMesh = this%mesh(1)%getnCondition()
  end function getnConditionOneMesh

  integer(ikind) pure function getnConditionMultiMesh(this, meshID)
    implicit none
    class(ModelDT), intent(in) :: this
    integer(ikind), intent(in) :: meshID
    getnConditionMultiMesh = this%mesh(meshID)%getnCondition()
  end function getnConditionMultiMesh

  integer(ikind) pure function getIDOneMesh(this)
    implicit none
    class(ModelDT), intent(in) :: this
    getIDOneMesh = this%mesh(1)%getID()
  end function getIDOneMesh

  integer(ikind) pure function getIDMultiMesh(this, meshID)
    implicit none
    class(ModelDT), intent(in) :: this
    integer(ikind), intent(in) :: meshID
    getIDMultiMesh = this%mesh(meshID)%getID()
  end function getIDMultiMesh

  type(NodePtrDT) function getNodeOneMesh(this, nodeID)
    implicit none
    class(ModelDT), intent(inout) :: this
    integer(ikind), intent(in)    :: nodeID
    getNodeOneMesh = this%mesh(1)%getNode(nodeID)
  end function getNodeOneMesh

  type(NodePtrDT) function getNodeMultiMesh(this, meshID, nodeID)
    implicit none
    class(ModelDT), intent(inout) :: this
    integer(ikind), intent(in)    :: meshID
    integer(ikind), intent(in)    :: nodeID
    getNodeMultiMesh = this%mesh(meshID)%getNode(nodeID)
  end function getNodeMultiMesh

  type(ElementPtrDT) function getElementOneMesh(this, elementID)
    implicit none
    class(ModelDT), intent(inout) :: this
    integer(ikind), intent(in)    :: elementID
    getElementOneMesh = this%mesh(1)%getElement(elementID)
  end function getElementOneMesh

  type(ElementPtrDT) function getElementMultiMesh(this, meshID, elementID)
    implicit none
    class(ModelDT), intent(inout) :: this
    integer(ikind), intent(in)    :: meshID
    integer(ikind), intent(in)    :: elementID
    getElementMultiMesh = this%mesh(meshID)%getElement(elementID)
  end function getElementMultiMesh

  type(ConditionPtrDT) function getConditionOneMesh(this, conditionID)
    implicit none
    class(ModelDT), intent(inout) :: this
    integer(ikind), intent(in)    :: conditionID
    getConditionOneMesh = this%mesh(1)%getCondition(conditionID)
  end function getConditionOneMesh

  type(ConditionPtrDT) function getConditionMultiMesh(this, meshID, conditionID)
    implicit none
    class(ModelDT), intent(inout) :: this
    integer(ikind), intent(in)    :: meshID
    integer(ikind), intent(in)    :: conditionID
    getConditionMultiMesh = this%mesh(meshID)%getCondition(conditionID)
  end function getConditionMultiMesh

  function getNodesOneMesh(this)
    implicit none
    class(ModelDT) , intent(inout)                     :: this
    type(NodePtrDT), dimension(this%getnNodeOneMesh()) :: getNodesOneMesh
    getNodesOneMesh = this%mesh(1)%node
  end function getNodesOneMesh

  function getNodesMultiMesh(this, meshID)
    implicit none
    class(ModelDT) , intent(inout)                             :: this
    integer(ikind) , intent(in)                                :: meshID
    type(NodePtrDT), dimension(this%getnNodeMultiMesh(meshID)) :: getNodesMultiMesh
    getNodesMultiMesh = this%mesh(meshID)%node
  end function getNodesMultiMesh

  function getElementsOneMesh(this)
    implicit none
    class(ModelDT)    , intent(inout)                        :: this
    type(ElementPtrDT), dimension(this%getnElementOneMesh()) :: getElementsOneMesh
    getElementsOneMesh = this%mesh(1)%element
  end function getElementsOneMesh

  function getElementsMultiMesh(this, meshID)
    implicit none
    class(ModelDT)    , intent(inout)                                :: this
    integer(ikind)    , intent(in)                                   :: meshID
    type(ElementPtrDT), dimension(this%getnElementMultiMesh(meshID)) :: getElementsMultiMesh
    getElementsMultiMesh = this%mesh(meshID)%element
  end function getElementsMultiMesh

  function getConditionsOneMesh(this)
    implicit none
    class(ModelDT)      , intent(inout)                          :: this
    type(ConditionPtrDT), dimension(this%getnConditionOneMesh()) :: getConditionsOneMesh
    getConditionsOneMesh = this%mesh(1)%condition
  end function getConditionsOneMesh

  function getConditionsMultiMesh(this, meshID)
    implicit none
    class(ModelDT)      , intent(inout)                                  :: this
    integer(ikind)      , intent(in)                                     :: meshID
    type(ConditionPtrDT), dimension(this%getnConditionMultiMesh(meshID)) :: getConditionsMultiMesh
    getConditionsMultiMesh = this%mesh(meshID)%condition
  end function getConditionsMultiMesh

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
  
  type(MeshDT) function getMesh(this, meshID)
    implicit none
    class(ModelDT), intent(in) :: this
    integer(ikind), intent(in) :: meshID
    getMesh = this%mesh(meshID)
  end function getMesh

  type(ProcessInfoDT) pure function getProcessInfo(this)
    implicit none
    class(ModelDT), intent(in) :: this
    getProcessInfo = this%processInfo
  end function getProcessInfo
  
end module ModelM
