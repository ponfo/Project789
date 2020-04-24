module PointM
  use UtilitiesM

  implicit none

  private
  public :: PointDT, point
  
  type PointDT
     integer(ikind)                         :: id
     real(rkind), dimension(:), allocatable :: coord
   contains
     procedure, public :: initPoint1D
     procedure, public :: initPoint2D
     procedure, public :: initPoint3D
     procedure, public :: setID
     procedure, public :: setX
     procedure, public :: setY
     procedure, public :: setZ
     procedure, public :: getID
     procedure, public :: getX
     procedure, public :: getY
     procedure, public :: getZ
     procedure, public :: getDimension
  end type PointDT

  interface point
     procedure :: constructor1D
     procedure :: constructor2D
     procedure :: constructor3D
  end interface point

contains

  type(PointDT) function constructor1D(id, x)
    implicit none
    integer(ikind), intent(in) :: id
    real(rkind)   , intent(in) :: x
    call constructor1D%initPoint1D(id, x)
  end function constructor1D
  subroutine initPoint1D(this, id, x)
    implicit none
    class(PointDT), intent(inout) :: this
    integer(ikind), intent(in)    :: id
    real(rkind)   , intent(in)    :: x
    this%id = id
    allocate(this%coord(1))
    this%coord(1) = x
  end subroutine initPoint1D

  type(PointDT) function constructor2D(id, x, y)
    implicit none
    integer(ikind), intent(in) :: id
    real(rkind)   , intent(in) :: x
    real(rkind)   , intent(in) :: y
    call constructor2D%initPoint2D(id, x, y)
  end function constructor2D
  subroutine initPoint2D(this, id, x, y)
    implicit none
    class(PointDT), intent(inout) :: this
    integer(ikind), intent(in)    :: id
    real(rkind)   , intent(in)    :: x
    real(rkind)   , intent(in)    :: y
    this%id = id
    allocate(this%coord(2))
    this%coord(1) = x
    this%coord(2) = y
  end subroutine initPoint2D

  type(PointDT) function constructor3D(id, x, y, z)
    implicit none
    integer(ikind), intent(in) :: id
    real(rkind)   , intent(in) :: x
    real(rkind)   , intent(in) :: y
    real(rkind)   , intent(in) :: z
    call constructor3D%initPoint3D(id, x, y, z)
  end function constructor3D
  subroutine initPoint3D(this, id, x, y, z)
    implicit none
    class(PointDT), intent(inout) :: this
    integer(ikind), intent(in)    :: id
    real(rkind)   , intent(in)    :: x
    real(rkind)   , intent(in)    :: y
    real(rkind)   , intent(in)    :: z
    this%id = id
    allocate(this%coord(3))
    this%coord(1) = x
    this%coord(2) = y
    this%coord(3) = z
  end subroutine initPoint3D

  subroutine setID(this, id)
    implicit none
    class(PointDT), intent(inout) :: this
    integer(ikind), intent(in)    :: id
    this%id = id
  end subroutine setID

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

  integer(ikind) function getID(this)
    implicit none
    class(PointDT), intent(inout) :: this
    getID = this%id
  end function getID
  
  real(rkind) function getX(this)
    implicit none
    class(PointDT), intent(inout) :: this
    getX = this%coord(1)
  end function getX

  real(rkind) function getY(this)
    implicit none
    class(PointDT), intent(inout) :: this
    getY = this%coord(2)
  end function getY

  real(rkind) function getZ(this)
    implicit none
    class(PointDT), intent(inout) :: this
    getZ = this%coord(3)
  end function getZ

  integer(ikind) function getDimension(this)
    implicit none
    class(PointDT), intent(inout) :: this
    getDimension = size(this%coord)
  end function getDimension

end module PointM
