include 'mkl_pardiso.f90'
  
module mklPardisoM

  !***********************************************
  !*                 EXTERNS                     *
  !***********************************************
  use UtilitiesM
  use quickSortM
  use sparseKit
  use mkl_pardiso
  
  implicit none
  
  private
  public :: pardisoMKL

contains
  
  subroutine pardisoMKL(p, maxfct, mnum, mtype, phase, stiffness,&
       idum, nrhs, iparm, msglvl, rhs, dof, error)
    implicit none
    class(Sparse), intent(inout) :: stiffness
    type(mkl_pardiso_handle) :: pt(64)
    real(rkind) :: ddum(1), p(64), rhs(stiffness%n)
    real(rkind) :: dof(stiffness%n)
    integer(ikind) :: maxfct, mnum, mtype, phase, n, nnz, i, j
    integer(ikind) :: nrhs, iparm(64), msglvl, error, error1, idum(1)
    pt(1:64)%DUMMY = p(1:64)
    dof(1:stiffness%n) = rhs(1:stiffness%n)
    call pardiso (pt, maxfct, mnum, mtype, phase, stiffness%n, stiffness%a, &
         stiffness%ai, stiffness%aj, idum, nrhs, iparm, msglvl, rhs, dof, error)
    if (error /= 0) then
       write(*,'(a,i5)') 'the following error was detected: ', error
       stop
    end if
    phase = -1 ! release internal memory
    call pardiso (pt, maxfct, mnum, mtype, phase, n, ddum, idum, idum,&
         idum, nrhs, iparm, msglvl, ddum, ddum, error1)
    if (error1 /= 0) then
       write(*,'(a,i5)') 'the following error on release stage was detected: ', error1
       stop
    end if
    print'(a)', 'the system has been solved'
    return
  end subroutine pardisoMKL
  
end module MklPardisoM

