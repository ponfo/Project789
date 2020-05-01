module PointM
  use UtilitiesM

  implicit none

  private
  public :: PointDT, point
  
  type PointDT
     real(rkind), dimension(:), allocatable :: coord
   contains
     procedure, public :: initPoint1D
     procedure, public :: initPoint2D
     procedure, public :: initPoint3D
     generic  , public :: updatePoint    => updatePoint1D, updatePoint2D, updatePoint3D
     procedure, public :: updatePoint1D
     procedure, public :: updatePoint2D
     procedure, public :: updatePoint3D
     procedure, public :: setX
     procedure, public :: setY
     procedure, public :: setZ
     procedure, public :: getX
     procedure, public :: getY
     procedure, public :: getZ
     procedure, public :: getDimension
     procedure, public :: free
  end type PointDT

  interface point
     procedure :: constructor1D
     procedure :: constructor2D
     procedure :: constructor3D
  end interface point


contains

  type(PointDT) function constructor1D(x)
    implicit none
    real(rkind)   , intent(in) :: x
    call constructor1D%initPoint1D(x)
  end function constructor1D
  subroutine initPoint1D(this, x)
    implicit none
    class(PointDT), intent(inout) :: this
    real(rkind)   , intent(in)    :: x
    allocate(this%coord(1))
    this%coord(1) = x
  end subroutine initPoint1D

  type(PointDT) function constructor2D(x, y)
    implicit none
    real(rkind)   , intent(in) :: x
    real(rkind)   , intent(in) :: y
    call constructor2D%initPoint2D(x, y)
  end function constructor2D
  subroutine initPoint2D(this, x, y)
    implicit none
    class(PointDT), intent(inout) :: this
    real(rkind)   , intent(in)    :: x
    real(rkind)   , intent(in)    :: y
    allocate(this%coord(2))
    this%coord(1) = x
    this%coord(2) = y
  end subroutine initPoint2D

  type(PointDT) function constructor3D(x, y, z)
    implicit none
    real(rkind)   , intent(in) :: x
    real(rkind)   , intent(in) :: y
    real(rkind)   , intent(in) :: z
    call constructor3D%initPoint3D(x, y, z)
  end function constructor3D
  subroutine initPoint3D(this, x, y, z)
    implicit none
    class(PointDT), intent(inout) :: this
    real(rkind)   , intent(in)    :: x
    real(rkind)   , intent(in)    :: y
    real(rkind)   , intent(in)    :: z
    allocate(this%coord(3))
    this%coord(1) = x
    this%coord(2) = y
    this%coord(3) = z
  end subroutine initPoint3D

  subroutine updatePoint1D(this, x)
    implicit none
    class(PointDT), intent(inout) :: this
    real(rkind)   , intent(in)    :: x
    this%coord(1) = x
  end subroutine updatePoint1D

  subroutine updatePoint2D(this, x, y)
    implicit none
    class(PointDT), intent(inout) :: this
    real(rkind)   , intent(in)    :: x
    real(rkind)   , intent(in)    :: y
    this%coord(1) = x
    this%coord(2) = y
  end subroutine updatePoint2D

  subroutine updatePoint3D(this, x, y, z)
    implicit none
    class(PointDT), intent(inout) :: this
    real(rkind)   , intent(in)    :: x
    real(rkind)   , intent(in)    :: y
    real(rkind)   , intent(in)    :: z
    this%coord(1) = x
    this%coord(2) = y
    this%coord(3) = z
  end subroutine updatePoint3D

  subroutine setX(this, x)
    implicit none
    class(PointDT), intent(inout) :: this
    real(rkind)   , intent(in)    :: x
    this%coord(1) = x
  end subroutine setX

  subroutine setY(this, y)
    implicit none
    class(PointDT), intent(inout) :: this
    real(rkind)   , intent(in)    :: y
    this%coord(2) = y
  end subroutine setY
  
  subroutine setZ(this, z)
    implicit none
    class(PointDT), intent(inout) :: this
    real(rkind)   , intent(in)    :: z
    this%coord(3) = z
  end subroutine setZ
  
  real(rkind) pure function getX(this)
    implicit none
    class(PointDT), intent(in) :: this
    getX = this%coord(1)
  end function getX

  real(rkind) pure function getY(this)
    implicit none
    class(PointDT), intent(in) :: this
    getY = this%coord(2)
  end function getY

  real(rkind) pure function getZ(this)
    implicit none
    class(PointDT), intent(in) :: this
    getZ = this%coord(3)
  end function getZ

  integer(ikind) function getDimension(this)
    implicit none
    class(PointDT), intent(inout) :: this
    getDimension = size(this%coord)
  end function getDimension

  subroutine free(this)
    implicit none
    class(PointDT), intent(inout) :: this
    deallocate(this%coord)
  end subroutine free

end module PointM
