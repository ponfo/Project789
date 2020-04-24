include 'mkl_pardiso.f90'
module mklPardisoM

  use UtilitiesM
  use SparseKit
  
  use mkl_pardiso
  
  use DirectLinearSolverM
  
  implicit none
  
  private
  public :: MKLpardisoDT
  
  type, extends(DirectLinearSolverDT) :: MKLpardisoDT
   contains
     procedure :: solveSystem => pardisoMKL
  end type MKLpardisoDT
  
contains
  
  subroutine pardisoMKL(this, vector, matrix, solution, arg)
    implicit none
    class(MKLpardisoDT)                  , intent(inout)  :: this
    class(Sparse)                          , intent(inout)  :: matrix
    real(rkind)             , dimension(:) , intent(inout)  :: vector
    real(rkind)             , dimension(:) , intent(inout)  :: solution
    integer(ikind)          , dimension(:) , intent(inout)  :: arg
    type(mkl_pardiso_handle), dimension(64)                 :: pt
    real(rkind)             , dimension(1)                  :: ddum
    integer(ikind)                                          :: maxfct
    integer(ikind)                                          :: mnum
    integer(ikind)                                          :: mtype
    integer(ikind)                                          :: phase
    integer(ikind)                                          :: nrhs
    integer(ikind)          , dimension(64)                 :: iparm(64)
    integer(ikind)                                          :: msglvl
    integer(ikind)                                          :: error
    integer(ikind)          , dimension(1)                  :: idum
    integer(ikind)                                          :: i
    pt(1:64)%DUMMY = arg(1:64)
    maxfct        = arg(65)
    mnum          = arg(66)
    mtype         = arg(67)
    phase         = arg(68)
    idum          = arg(69)
    nrhs          = arg(70)
    do i = 1, 64
       iparm(i)   = arg(70+i)
    end do
    msglvl        = arg(136)
    error         = arg(137)
    solution      = vector
    call pardiso (pt, maxfct, mnum, mtype, phase, matrix%getn()  &
         , matrix%geta(), matrix%getai(), matrix%getaj(), idum, nrhs, iparm &
         , msglvl, vector, solution, error                  )
    if (error /= 0) then
       write(*,'(a,i5)') 'the following error was detected: ', error
       stop
    end if
    phase = -1 ! release internal memory
    call pardiso (pt, maxfct, mnum, mtype, phase, matrix%getn()  &
         , ddum, idum, idum, idum, nrhs, iparm, msglvl,     &
         ddum, ddum, error                                  )
    if (error /= 0) then
       write(*,'(a,i5)') 'the following error on release stage was detected: ', error
       stop
    end if
    print'(a)', 'the system has been solved'
    return
  end subroutine pardisoMKL
  
end module mklPardisoM
  
