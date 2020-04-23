module PreconditionerMod

  use tools
  use sparseKit

  implicit none

  private
  public :: PreconditionerTYPE

  type, abstract :: PreconditionerTYPE
   contains
     procedure(Preconditioner_procedure), deferred  :: usePreconditioner
  end type PreconditionerTYPE

  abstract interface
     subroutine Preconditioner_procedure(this, vector, matrix, solution, arg)
       import PreconditionerTYPE, sparse, rkind, ikind
       class(PreconditionerTYPE)            , intent(inout) :: this
       class(Sparse)                         , intent(inout) :: matrix
       real(rkind)            , dimension(:) , intent(inout) :: vector
       real(rkind)            , dimension(:) , intent(inout) :: solution
       integer(ikind)         , dimension(:) , intent(inout) :: arg
     end subroutine Preconditioner_procedure
  end interface

end module PreconditionerMod
