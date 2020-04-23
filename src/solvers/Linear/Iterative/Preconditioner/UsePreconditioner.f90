module UsePreconditionerM

  use UtilitiesM
  use SparseKit

  use PreconditionerM

  implicit none

  private
  public :: UsePreconditionerDT, SetPreconditioner

  type UsePreconditionerDT
     class(PreconditionerDT), allocatable :: preconditionerMethod
   contains
     procedure :: init
     procedure :: changePreconditioner
     procedure :: use
  end type UsePreconditionerDT

  interface SetPreconditioner
     procedure :: constructor
  end interface SetPreconditioner
  
contains

  type(UsePreconditionerDT) function constructor(preconditioner)
    implicit none
    class(PreconditionerDT), intent(inout) :: preconditioner
    call constructor%init(preconditioner)
  end function constructor

  subroutine init(this, method)
    implicit none
    class(UsePreconditionerDT), intent(inout) :: this
    class(PreconditionerDT)   , intent(inout) :: method
    allocate(this%preconditionerMethod, source = method)
  end subroutine init

  subroutine changePreconditioner(this, newMethod)
    implicit none
    class(UsePreconditionerDT), intent(inout) :: this
    class(PreconditionerDT)   , intent(inout) :: newMethod
    deallocate(this%preconditionerMethod)
    allocate(this%preconditionerMethod, source = newMethod)
  end subroutine changePreconditioner

  subroutine use(this, vector, matrix, solution, arg)
    implicit none
    class(UsePreconditionerDT)   , intent(inout) :: this
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

end module UsePreconditionerM
