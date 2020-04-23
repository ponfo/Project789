module PreconditionerM

  use UtilitiesM
  use SparseKit

  implicit none

  private
  public :: PreconditionerDT

  type, abstract :: PreconditionerDT
   contains
     procedure(Preconditioner_procedure), deferred  :: usePreconditioner
  end type PreconditionerDT

  abstract interface
     subroutine Preconditioner_procedure(this, vector, matrix, solution, arg)
       import PreconditionerDT, sparse, rkind, ikind
       class(PreconditionerDT)            , intent(inout) :: this
       class(Sparse)                         , intent(inout) :: matrix
       real(rkind)            , dimension(:) , intent(inout) :: vector
       real(rkind)            , dimension(:) , intent(inout) :: solution
       integer(ikind)         , dimension(:) , intent(inout) :: arg
     end subroutine Preconditioner_procedure
  end interface

end module PreconditionerM
