module Calculate_dtM

  use UtilitiesM

  use ProcessM

  use ElementPtrM

  use ThermalMaterialM

  use ThermalModelM

  implicit none

  private
  public :: Calculate_dtDT

  type, extends(NewProcessDT) :: Calculate_dtDT
   contains
     procedure :: useProcess => calculatedt
     procedure :: calculate
  end type Calculate_dtDT

contains

  real(rkind) function calculate(this, model)
    implicit none
    class(Calculate_dtDT), intent(inout)   :: this
    class(ThermalModelDT), intent(inout)   :: model
    type(ElementPtrDT)                     :: element
    real(rkind), dimension(:), allocatable :: dt
    integer(ikind)                         :: nElem
    integer(ikind)                         :: i
    integer(ikind)                         :: n
    nElem = model%getnElement()
    allocate(dt(nElem))
    do i = 1, nElem
       element = model%getElement(i)
       call element%calculateDt(dt(i))
    end do
    calculate = minval(dt)
    print*, '*** Calculate dt => dt = ', calculate,' ***'
  end function calculate

  subroutine calculatedt(this)
    implicit none
    class(Calculate_dtDT), intent(inout) :: this
  end subroutine calculatedt

end module Calculate_dtM
