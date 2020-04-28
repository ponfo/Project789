module SchemeM

  implicit none

  private
  public :: NewSchemeDT, SchemeDT, SetScheme

  type, abstract :: NewSchemeDT
   contains
     procedure(Scheme_Procedure), deferred :: useScheme
     !procedure(integrator_interface), deferred :: integrate
  end type NewSchemeDT

  abstract interface
     subroutine Scheme_Procedure(this)
       import NewSchemeDT
       class(NewSchemeDT), intent(inout) :: this
     end subroutine Scheme_Procedure
  end interface

  interface SetScheme
     procedure :: constructor
  end interface SetScheme
!!$abstract interface
!!$     subroutine integrator_interface(this,dt)
!!$       import :: ProcessDT
!!$       class(ProcessDT),intent(inout) :: this     ! integrand
!!$       real(rkind)     ,intent(in)    :: dt       ! time step size
!!$     end subroutine integrator_interface
!!$  end interface
  
  type SchemeDT
     class(NewSchemeDT), allocatable :: scheme
   contains
     procedure :: init
     procedure :: change
!!$     procedure :: InitializeElements
!!$     procedure :: InitializeSolutionStepStrategy
!!$     procedure :: FinalizeSolutionStep
!!$     procedure :: InitializeNonLinearIteration
!!$     procedure :: FinalizeNonLinearIteration
!!$     procedure :: Predict
!!$     procedure :: Update
     procedure :: use !outputData
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

  subroutine use(this)
    implicit none
    class(SchemeDT)   , intent(inout) :: this
    call this%scheme%useScheme()
  end subroutine use
  
end module SchemeM
