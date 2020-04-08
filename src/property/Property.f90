module PropertyM
  use UtilitiesM
  use FortranParser

  implicit none

  private
  public :: PropertyDT, property

  type :: PropertyDT
     real(rkind)         , allocatable :: val
     type(EquationParser), allocatable :: f
   contains
     procedure, public :: init
  end type PropertyDT

  interface property
     procedure :: contructor
  end interface property

contains

  type(PropertyDT) function constructor(type)
    implicit NONE
    character(*), intent(in) :: type
    call constructor%init(type)
  end function constructor

  subroutine init(this, type)
    implicit none
    class(PropertyDT), intent(inout) :: this
    character(*)     , intent(in)    :: type
    if(trim(type) == 'constant') then
       allocate(this%val)
    else if(trim(type) == 'variable' .or. &
         trim(type) == 'function')   then
       allocate(this%f)
    end if
  end subroutine init

end module PropertyM
