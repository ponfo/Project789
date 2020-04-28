module PreconditionerMethodM

  use UtilitiesM
  use SparseKit

  use PreconditionerM

  implicit none

  private
  public :: PreconditionerMethodDT

  type, extends(PreconditionerDT) :: PreconditionerMethodDT
   contains
     procedure :: usePreconditioner => method
  end type PreconditionerMethodDT

contains

  subroutine method(this, vector, matrix, solution, arg)
    implicit none
    class(PreconditionerMethodDT), intent(inout) :: this
    class(Sparse)                   , intent(inout) :: matrix
    real(rkind)   , dimension(:)    , intent(inout) :: vector
    real(rkind)   , dimension(:)    , intent(inout) :: solution
    integer(ikind), dimension(:)    , intent(inout) :: arg
    write(*,*) 'Preconditioner Method Implementation'
  end subroutine method

end module PreconditionerMethodM
