module BiConjugateGradientM

  use UtilitiesM
  use SparseKit

  use IterativeLinearSolverM

  implicit none

  private
  public :: BiConjugateGradientDT

  type, extends(IterativeLinearSolverDT) :: BiConjugateGradientDT
   contains
     procedure :: solveSystem => method
  end type BiConjugateGradientDT

contains
       
  subroutine method(this, vector, matrix, solution, arg)
    implicit none
    class(BiConjugateGradientDT), intent(inout)  :: this
    class(Sparse)               , intent(inout)  :: matrix
    real(rkind)   , dimension(:), intent(inout)  :: vector
    real(rkind)   , dimension(:), intent(inout)  :: solution
    integer(ikind), dimension(:), intent(inout)  :: arg
    real(rkind)   , dimension(:), allocatable    :: y, p, r, z
    integer(ikind), dimension(:), allocatable    :: spRowptr, spIdx, x_fixIdx
    real(rkind)   , dimension(:), allocatable    :: spMtx, x
    real(rkind)   , dimension(:), allocatable    :: b, diagMtx, x_fix
    real(rkind)    :: tol, err_old, err_new, alfa, beta, py
    integer(ikind) :: npoin, nfix, k
    integer(ikind) :: i, j
    npoin = matrix%getn()
    spMtx = matrix%getA()
    spIdx = matrix%getAJ()
    spRowptr = matrix%getAI()
    allocate(diagMtx(npoin))
    x = solution
    nfix = 0
    do i = 1, npoin
       if (sqrt(x(i)**2) .ne. 0.d0) then
          nfix = nfix + 1
       end if
    end do
    allocate(x_fixIdx(nfix), x_fix(nfix))
    do i = 1, npoin
       if (sqrt(x(i)**2) .ne. 0.d0) then
          x_fix(i) = x(i)
          x_fixIdx(i) = i
       end if
    end do
    b = vector
    do i = 1, npoin
       do j = spRowptr(i), spRowptr(i+1)-1
          if (i == spIdx(j)) then
             diagMtx = spMtx(j) 
          end if
       end do
    end do
    allocate(y(npoin))
    allocate(p(npoin))
    allocate(r(npoin))
    allocate(z(npoin))
    k = 0
    tol = 1.d-10
    call copy1(nfix, 1.d0, x_fix, x_fixIdx, npoin, x)
    call SpMV(spMtx, spIdx, spRowptr, x, y, npoin, spRowptr(npoin+1))
    call copy2(npoin, 1.d30, x, nfix, x_fixIdx, y)
    call vecsum(npoin, -1.d0, y, b, r)
    call assign2(npoin, r, nfix, x_fixIdx, 0.d0)
    if(vecdot(npoin, r, r) < tol) return
    call vecdiv(npoin, r, diagMtx, p)
    err_new = vecdot(npoin, r, p)
    call SpMV(spMtx, spIdx, spRowptr, p, y, npoin, spRowptr(npoin+1))
    call copy2(npoin, 1.d30, p, nfix, x_fixIdx, y)
    py = vecdot(npoin, p, y)
    !DEBUGGGGG
    alfa = err_new/py
    call vecsum(npoin, alfa, p, x, x)
    err_old = err_new
    do while(dabs(err_old) > tol .and. k < 1000)
       k = k + 1
       call vecsum(npoin, -alfa, y, r, r)
       call vecdiv(npoin, r, diagMtx, z)
       err_new = vecdot(npoin, r, z)
       beta = err_new/err_old
       call vecsum(npoin, beta, p, z, p)
       call SpMV(spMtx, spIdx, spRowptr, p, y, npoin, spRowptr(npoin+1))
       call copy2(npoin, 1.d30, p, nfix, x_fixIdx, y)
       py = vecdot(npoin, p, y)
       !DEBUGGGGG
       alfa = err_new/py
       call vecsum(npoin, alfa, p, x, x)
       err_old = err_new
    end do
    solution = x
    return
  end subroutine method

  subroutine vecdiv(n, x, y, z)
    !---------------------------------
    !		Calcula z = x/y
    !---------------------------------
    implicit none
    integer(ikind) :: n, i
    real(rkind)    :: x(n), y(n), z(n)
    !$OMP PARALLEL DO
    do i = 1, n
       z(i) = x(i)/y(i)	
    end do
    !$OMP END PARALLEL DO
  end subroutine vecdiv

  subroutine vecsum(n, alfa, x, y, z)
    !-----------------------------------
    !		Calcula z = alfa*x + y
    !-----------------------------------
    implicit none
    integer(ikind) :: n, i
    real(rkind)    :: x(n), y(n), z(n), alfa
    !$OMP PARALLEL DO
    do i = 1, n
       z(i) = alfa*x(i) + y(i)
    end do
    !$OMP END PARALLEL DO
  end subroutine vecsum

  subroutine assign(n, y, scal)
    !-----------------------------
    !		y(:) = scal
    !-----------------------------
    implicit none
    integer(ikind) :: n, i
    real(rkind)    :: scal, y(n)
    !$OMP PARALLEL DO PRIVATE(i)
    do i = 1, n
       y(i) = scal
    end do
    !$OMP END PARALLEL DO
  end subroutine assign

  subroutine assign2(n, y, m, idx, scal)
    implicit none
    integer(ikind) :: n, m, i
    real(rkind)    :: y(n), scal
    integer(ikind) :: idx(m)
    !$OMP PARALLEL DO PRIVATE(i)
    do i = 1, m
       y(idx(i)) = scal
    end do
    !$OMP END PARALLEL DO
  end subroutine assign2

  subroutine copy1(m, alfa, x, idx, n, y)
    !-----------------------------------------
    !			y(idx(:)) = alfa*x
    !-----------------------------------------
    implicit none
    integer(ikind) :: n, m, i
    real(rkind)    :: x(m), y(n), alfa
    integer(ikind) :: idx(m)
    !$OMP PARALLEL DO PRIVATE(i)
    do i=1, m
       y(idx(i)) = alfa*x(i)
    end do
    !$OMP END PARALLEL DO
  end subroutine copy1

  subroutine copy2(n, alfa, x, m, idx, y)
    !-----------------------------------------
    !			y(idx(:)) = alfa*x(idx(:))
    !-----------------------------------------
    implicit none
    integer(ikind) :: n, m, i
    real(rkind)    :: x(n), y(n), alfa
    integer(ikind) :: idx(m)
    !$OMP PARALLEL DO PRIVATE(i)
    do i=1, m
       y(idx(i)) = alfa*x(idx(i))
    end do
    !$OMP END PARALLEL DO
  end subroutine copy2

  real(8) function vecdot(n, x, y)
    !-----------------------------------
    !	Devuelve producto punto <x,y>
    !-----------------------------------
    implicit none
    integer(ikind) :: n, i
    real(rkind)    :: x(n), y(n), res
    res = 0.d0
    !$OMP PARALLEL DO REDUCTION(+:res)
    do i = 1, n
       res = res + x(i)*y(i)
    end do
    !$OMP END PARALLEL DO
    vecdot = res
  end function vecdot

  subroutine SpMV(spMtx, spIdx, spRowptr, v, y, npoin, npos)
    !---------------------------------------------------------
    !			Multiplicacion sparse y = M.v
    !---------------------------------------------------------
    implicit none
    integer(ikind) :: npoin, npos
    integer(ikind) :: i, j
    real(rkind)    :: spMtx(npos), v(npoin), y(npoin), dot
    integer(ikind) :: spIdx(npos), spRowptr(npoin+1)
    !$OMP PARALLEL DO PRIVATE(i, j, dot)
    do i = 1, npoin
       dot = 0.d0
       do j = spRowptr(i)+1, spRowptr(i+1)
          dot = dot + spMtx(j)*v(spIdx(j))
       end do
       y(i) = dot
    end do
    !$OMP END PARALLEL DO 
  end subroutine SpMV

end module BiConjugateGradientM

