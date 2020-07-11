module SchemeM
  
  implicit none

  private
  public :: NewSchemeDT, SchemeDT, SetScheme

  type, abstract :: NewSchemeDT
   contains
  end type NewSchemeDT
  
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
