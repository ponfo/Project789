module ReorderSystemM

  use UtilitiesM
  use SparseKit

  implicit none

  private
  public :: ReorderSystemDT

  type, abstract :: ReorderSystemDT
   contains
     procedure(ReorderSystem_procedure), deferred  :: useReorder
  end type ReorderSystemDT

  abstract interface
     subroutine ReorderSystem_procedure(this, vector, matrix, solution, arg)
       import ReorderSystemDT, sparse, rkind, ikind
       class(ReorderSystemDT)            , intent(inout) :: this
       class(Sparse)                         , intent(inout) :: matrix
       real(rkind)            , dimension(:) , intent(inout) :: vector
       real(rkind)            , dimension(:) , intent(inout) :: solution
       integer(ikind)         , dimension(:) , intent(inout) :: arg
     end subroutine ReorderSystem_procedure
  end interface

end module ReorderSystemM
