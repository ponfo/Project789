module CFDModelM
  use UtilitiesM

  use SparseKit

  use MeshM
  use ModelM

  use ResultsM
  
  implicit none

  private
  public :: CFDModelDT, cfdModel

  type, extends(modelDT) :: CFDModelDT
     type(Sparse)                           :: lhs
     type(Sparse)                           :: mass
     real(rkind), dimension(:), allocatable :: rhs
     real(rkind), dimension(:), allocatable :: dof
     type(ResultsDT)                        :: results
   contains
     procedure, public :: init
     procedure, public :: initWithoutSystem
     procedure, public :: initSystem
     procedure, public :: freeSystem
  end type CFDModelDT

  interface cfdModel
     procedure :: constructor
  end interface cfdModel
  
contains

  type(CFDModelDT) function constructor(nDof, nnz, id, nNode, nElement, nCondition)
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
    class(CFDModelDT), intent(inout) :: this
    integer(ikind)   , intent(in)    :: nDof
    integer(ikind)   , intent(in)    :: nnz
    integer(ikind)   , intent(in)    :: id
    integer(ikind)   , intent(in)    :: nNode
    integer(ikind)   , intent(in)    :: nElement
    integer(ikind)   , intent(in)    :: nCondition
    this%lhs  = sparse(nnz, nDof)
    this%mass = sparse(nnz, nDof)
    allocate(this%rhs(nDof))
    allocate(this%dof(nDof))
    call this%initModel(1) !una sola malla en el modelo
    this%mesh(1) = mesh(id, nNode, nElement, nCondition)
  end subroutine init

  subroutine initWithoutSystem(this, id, nNode, nElement, nCondition)
    implicit none
    class(CFDModelDT), intent(inout) :: this
    integer(ikind)   , intent(in)    :: id
    integer(ikind)   , intent(in)    :: nNode
    integer(ikind)   , intent(in)    :: nElement
    integer(ikind)   , intent(in)    :: nCondition
    call this%initModel(1)
    this%mesh(1) = mesh(id, nNode, nElement, nCondition)
  end subroutine initWithoutSystem

  subroutine initSystem(this, nDof, nnz)
    implicit none
    class(CFDModelDT), intent(inout) :: this
    integer(ikind)   , intent(in)    :: nDof
    integer(ikind)   , intent(in)    :: nnz
    this%lhs  = sparse(nnz, nDof)
    this%mass = sparse(nnz, nDof)
    allocate(this%rhs(nDof))
    allocate(this%dof(nDof))
  end subroutine initSystem
  
  subroutine freeSystem(this)
    implicit none
    class(CFDModelDT), intent(inout) :: this
    call this%lhs%free()
    call this%mass%free()
    if(allocated(this%rhs)) deallocate(this%rhs)
  end subroutine freeSystem
  
end module CFDModelM
