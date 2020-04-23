module UsePreconditionerMod

  use tools
  use sparseKit

  use PreconditionerMod

  implicit none

  private
  public :: UsePreconditionerTYPE, SetPreconditioner

  type UsePreconditionerTYPE
     class(PreconditionerTYPE), allocatable :: preconditionerMethod
   contains
     procedure :: init
     procedure :: changePreconditioner
     procedure :: use
  end type UsePreconditionerTYPE

  interface SetPreconditioner
     procedure :: constructor
  end interface SetPreconditioner
  
contains

  type(UsePreconditionerTYPE) function constructor(preconditioner)
    implicit none
    class(PreconditionerTYPE), intent(inout) :: preconditioner
    call constructor%init(preconditioner)
  end function constructor

  subroutine init(this, method)
    implicit none
    class(UsePreconditionerTYPE), intent(inout) :: this
    class(PreconditionerTYPE)   , intent(inout) :: method
    allocate(this%preconditionerMethod, source = method)
  end subroutine init

  subroutine changePreconditioner(this, newMethod)
    implicit none
    class(UsePreconditionerTYPE), intent(inout) :: this
    class(PreconditionerTYPE)   , intent(inout) :: newMethod
    deallocate(this%preconditionerMethod)
    allocate(this%preconditionerMethod, source = newMethod)
  end subroutine changePreconditioner

  subroutine use(this, vector, matrix, solution, arg)
    implicit none
    class(UsePreconditionerTYPE)   , intent(inout) :: this
    class(Sparse)               , intent(inout) :: matrix
    real(rkind)   , dimension(:), intent(inout) :: vector
    real(rkind)   , dimension(:), intent(inout) :: solution
    integer(ikind), dimension(:), intent(inout) :: arg
    if(allocated(this%preconditionerMethod)) then
       call this%preconditionerMethod%usePreconditioner(vector, matrix, solution, arg)
    else
       write(*,*) '*** Preconditioner Not Allocated ***'
    end if
  end subroutine use

end module UsePreconditionerMod
