module SchemeM

  use UtilitiesM
  
  use ProcessM
  
  implicit none

  private
  public :: NewSchemeDT, SchemeDT, SetScheme

  type, abstract :: NewSchemeDT
   contains
     procedure(integrator_interface), nopass, deferred :: integrate
  end type NewSchemeDT

  abstract interface
     subroutine integrator_interface(this, dt)
       import :: NewProcessDT, rkind
       class(NewProcessDT), intent(inout) :: this
       real(rkind)        , intent(in)    :: dt
     end subroutine integrator_interface
  end interface
  
  interface SetScheme
     procedure :: constructor
  end interface SetScheme
  
  type SchemeDT
     class(NewSchemeDT), allocatable :: scheme
   contains
     procedure :: init
     procedure :: change
  end type SchemeDT

contains

  type(SchemeDT) function constructor(scheme)
    implicit none
    class(NewSchemeDT), intent(in) :: scheme
    call constructor%init(scheme)
  end function constructor
  
  subroutine init(this, scheme)
    implicit none
    class(SchemeDT)   , intent(inout) :: this
    class(NewSchemeDT), intent(in)    :: scheme
    allocate(this%scheme, source = scheme)
  end subroutine init

  subroutine change(this, newScheme)
    implicit none
    class(SchemeDT)   , intent(inout) :: this
    class(NewSchemeDT), intent(in)    :: newScheme
    deallocate(this%scheme)
    allocate(this%scheme, source = newScheme)
  end subroutine change
  
end module SchemeM
