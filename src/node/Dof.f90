module DofM
  use UtilitiesM

  implicit none

  private
  public :: DofDT, newDof

  type DofDT
     real(rkind) , pointer     :: val
     real(rkind) , allocatable :: fixedVal
     character(:), pointer     :: name
     logical                   :: isFixed
   contains
     procedure, public :: init
     procedure, public :: initWithName
     procedure, public :: fixDof
     procedure, public :: freeDof
     procedure, public :: assignVal
     procedure, public :: assignName
     procedure, public :: getVal
     procedure, public :: getName
  end type DofDT

  interface newDof
     procedure :: constructor
     procedure :: constructorWithName
  end interface newDof

contains

  type(DofDT) function constructor(dof, isFixed)
    implicit none
    real(rkind), target, intent(in) :: dof
    logical            , intent(in) :: isFixed
    call constructor%init(dof, isFixed)
  end function constructor

  subroutine init(this, dof, isFixed)
    implicit none
    class(DofDT)          , intent(inout) :: this
    real(rkind)   , target, intent(in)    :: dof
    logical               , intent(in)    :: isFixed
    this%val     => dof
    this%name    => null()
    this%isFixed = isFixed
  end subroutine init

  type(DofDT) function constructorWithName(dof, name, isFixed)
    implicit none
    real(rkind) , target, intent(in) :: dof
    character(*), target, intent(in) :: name
    logical             , intent(in) :: isFixed
    call constructorWithName%initWithName(dof, name, isFixed)
  end function constructorWithName

  subroutine initWithName(this, dof, name, isFixed)
    implicit none
    class(DofDT)        , intent(inout) :: this
    real(rkind) , target, intent(in)    :: dof
    character(*), target, intent(in)    :: name
    logical             , intent(in)    :: isFixed
    this%val     => dof
    this%name    => name
    this%isFixed = isFixed
  end subroutine initWithName

  subroutine assignVal(this, dof)
    implicit none
    class(DofDT)        , intent(inout) :: this
    real(rkind) , target, intent(in)    :: dof
    this%val => dof
  end subroutine assignVal

  subroutine assignName(this, name)
    implicit none
    class(DofDT)        , intent(inout) :: this
    character(*), target, intent(in)    :: name
    this%name => name
  end subroutine assignName

  subroutine fixDof(this, fixedVal)
    implicit none
    class(DofDT), intent(inout) :: this
    real(rkind) , intent(in)    :: fixedVal
    this%isFixed = .true.
    if(allocated(this%fixedVal)) deallocate(this%fixedVal)
    allocate(this%fixedVal, source = fixedVal)
  end subroutine fixDof

  subroutine freeDof(this)
    implicit none
    class(DofDT), intent(inout) :: this
    this%isFixed = .false.
    if(allocated(this%fixedVal)) deallocate(this%fixedVal)
  end subroutine freeDof

  pure real(rkind) function getVal(this)
    implicit none
    class(DofDT), intent(in) :: this
    if(this%isFixed) then
       getVal = this%fixedVal
    else
       getVal = this%val
    end if
  end function getVal

  pure function getName(this)
    implicit none
    class(DofDT), intent(in)  :: this
    character(:), allocatable :: getName
    getName = this%name
  end function getName
  
end module DofM
