module ReorderSystemMod

  use tools
  use sparseKit

  implicit none

  private
  public :: ReorderSystemTYPE

  type, abstract :: ReorderSystemTYPE
   contains
     procedure(ReorderSystem_procedure), deferred  :: useReorder
  end type ReorderSystemTYPE

  abstract interface
     subroutine ReorderSystem_procedure(this, vector, matrix, solution, arg)
       import ReorderSystemTYPE, sparse, rkind, ikind
       class(ReorderSystemTYPE)            , intent(inout) :: this
       class(Sparse)                         , intent(inout) :: matrix
       real(rkind)            , dimension(:) , intent(inout) :: vector
       real(rkind)            , dimension(:) , intent(inout) :: solution
       integer(ikind)         , dimension(:) , intent(inout) :: arg
     end subroutine ReorderSystem_procedure
  end interface

end module ReorderSystemMod
