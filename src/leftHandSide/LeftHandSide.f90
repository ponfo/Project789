module LeftHandSideM
  use UtilitiesM

  implicit none

  private
  public :: LeftHandSideDT, leftHandSide

  type :: LeftHandSideDT
     real(rkind), dimension(:,:), allocatable :: stiffness
     real(rkind), dimension(:,:), allocatable :: dampness
     real(rkind), dimension(:,:), allocatable :: mass
   contains
     procedure, public :: init
     procedure, public :: free
  end type LeftHandSideDT

  interface leftHandSide
     procedure :: constructor
  end interface leftHandSide

contains

  type(LeftHandSideDT) function constructor(massSize, dampnessSize, stiffnessSize)
    implicit none
    integer(ikind), intent(in) :: massSize
    integer(ikind), intent(in) :: dampnessSize
    integer(ikind), intent(in) :: stiffnessSize
    call constructor%init(massSize, dampnessSize, stiffnessSize)
  end function constructor

  subroutine init(this, massSize, dampnessSize, stiffnessSize)
    implicit none
    class(LeftHandSideDT), intent(inout) :: this
    integer(ikind)       , intent(in)    :: massSize
    integer(ikind)       , intent(in)    :: dampnessSize
    integer(ikind)       , intent(in)    :: stiffnessSize
    allocate(this%mass(massSize,massSize))
    allocate(this%dampness(dampnessSize,dampnessSize))
    allocate(this%stiffness(stiffnessSize,stiffnessSize))
  end subroutine init

  subroutine free(this)
    implicit none
    class(LeftHandSideDT), intent(inout) :: this
    if(allocated(this%stiffness)) deallocate(this%stiffness)
    if(allocated(this%dampness)) deallocate(this%dampness)
    if(allocated(this%mass)) deallocate(this%mass)
  end subroutine free

end module LeftHandSideM
    
