module ThermalElement2DM

  use ElementM

  implicit none

  private
  public :: ThermalElementDT, thermalElement, initGeometries

  type, extends(ElementDT) :: ThermalElementDT
     type(ThermalMaterialDT), pointer :: material
   contains
     procedure, public :: init
  end type ThermalElementDT

  interface thermalElement2D
     procedure :: constructor
  end interface thermalElement2D

  type(Triangle2D3NodeDT)     , save :: triangle2D3Node
  type(Triangle2D6NodeDT)     , save :: triangle2D6Node
  type(Quadrilateral2D4NodeDT), save :: Quadrilateral2D4Node
  type(Quadrilateral2D8NodeDT), save :: Quadrilateral2D8Node

contains

  type(ThermalElementDT) function contructor(id, node, material)
    implicit none
    integer(ikind)                       , intent(in) :: id
    type(NodePtrDT)        , dimension(:), intent(in) :: node
    type(ThermalMaterialDT), target      , intent(in) :: material
    call this%init(id, nGauss, node, material)
  end function contructor

  subroutine init(this, id, node, material)
    implicit none
    class(ThermalElementDT)              , intent(inout) :: this
    integer(ikind)                       , intent(in)    :: id
    type(NodePtrDT)        , dimension(:), intent(in)    :: node
    type(ThermalMaterialDT), target      , intent(in)    :: material
    this%id = id
    this%node = node
    this%material = material
    if(size(node) == 3) then
       this%geometry => triangle2D3Node
    else if(size(node) == 4) then
       this%geometry => quadrilateral2D4Node
    else if(size(node) == 6) then
       this%geometry => triangle2D6Node
    else if(size(node) == 8) then
       this%geometry => quadrilateral2D8Node
    end if
  end subroutine init

  subroutine initGeometries(nGauss)
    implicit none
    integer(ikind), intent(in) :: nGauss
    triangle2D3Node = triangle2D3Node(nGauss)
    triangle2D6Node = triangle2D6Node(nGauss)
    quadrilateral2D4Node = quadrilateral2D4Node(nGauss)
    quadrilateral2D8Node = quadrilateral2D8Node(nGauss)
  end subroutine initGeometries

end module ThermalElement2DM
