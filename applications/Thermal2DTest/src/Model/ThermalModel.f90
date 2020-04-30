module ThermalModelM
  use UtilitiesM

  use SparseKit

  use MeshM
  use ModelM
  
  implicit none

  private
  public :: ThermalModelDT, thermalModel

  type, extends(modelDT) :: ThermalModelDT
     type(Sparse)                           :: lhs
     real(rkind), dimension(:), allocatable :: rhs
     real(rkind), dimension(:), allocatable :: dof
   contains
     procedure, public :: init
  end type ThermalModelDT

  interface thermalModel
     procedure :: constructor
  end interface thermalModel
  
contains

  type(ThermalModelDT) function constructor(nDof, nnz, id, nNode, nElement, nCondition)
    implicit none
    integer(ikind), intent(in) :: nDof
    integer(ikind), intent(in) :: nnz
    integer(ikind), intent(in) :: id
    integer(ikind), intent(in) :: nNode
    integer(ikind), intent(in) :: nElement
    integer(ikind), intent(in) :: nCondition
    call constructor%init(nDof, nnz, id, nNode, nElement, nCondition)
  end function constructor

  subroutine init(this, nDof, nnz, id, nNode, nElement, nCondition)
    implicit none
    class(ThermalModelDT), intent(inout) :: this
    integer(ikind)       , intent(in)    :: nDof
    integer(ikind)       , intent(in)    :: nnz
    integer(ikind)       , intent(in)    :: id
    integer(ikind)       , intent(in)    :: nNode
    integer(ikind)       , intent(in)    :: nElement
    integer(ikind)       , intent(in)    :: nCondition
    this%lhs = sparse(nnz, nDof)
    allocate(this%rhs(nDof))
    allocate(this%dof(nDof))
    call this%initModel(1) !una sola malla en el modelo
    this%mesh(1) = mesh(id, nNode, nElement, nCondition)
  end subroutine init
  
end module ThermalModelM
