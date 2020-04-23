module UseReorderSystemMod

  use tools
  use sparseKit

  use ReorderSystemMod

  implicit none

  private
  public :: UseReorderSystemTYPE, SetReorderSystem

  type UseReorderSystemTYPE
     class(ReorderSystemTYPE), allocatable :: reorderMethod
   contains
     procedure :: init
     procedure :: changeMethod
     procedure :: use
  end type UseReorderSystemTYPE

  interface SetReorderSystem
     procedure :: constructor
  end interface SetReorderSystem
  
contains

  type(UseReorderSystemTYPE) function constructor(method)
    implicit none
    class(ReorderSystemTYPE), intent(inout) :: method
    call constructor%init(method)
  end function constructor

  subroutine init(this, method)
    implicit none
    class(UseReorderSystemTYPE), intent(inout) :: this
    class(ReorderSystemTYPE)   , intent(inout) :: method
    allocate(this%reorderMethod, source = method)
  end subroutine init

  subroutine changeMethod(this, newMethod)
    implicit none
    class(UseReorderSystemTYPE), intent(inout) :: this
    class(ReorderSystemTYPE)   , intent(inout) :: newMethod
    deallocate(this%reorderMethod)
    allocate(this%reorderMethod, source = newMethod)
  end subroutine changeMethod

  subroutine use(this, vector, matrix, solution, arg)
    implicit none
    class(UseReorderSystemTYPE)   , intent(inout) :: this
    class(Sparse)               , intent(inout) :: matrix
    real(rkind)   , dimension(:), intent(inout) :: vector
    real(rkind)   , dimension(:), intent(inout) :: solution
    integer(ikind), dimension(:), intent(inout) :: arg
    call this%reorderMethod%useReorder(vector, matrix, solution, arg)
  end subroutine use

end module UseReorderSystemMod
