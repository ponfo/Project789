module BaseIntegrandM

  use UtilitiesM

  use ProcessM
  
  implicit none

  private
  public :: BaseIntegrandDT

  type, abstract :: BaseIntegrandDT
   contains
     procedure(integrator_interface), nopass, deferred :: integrate
  end type BaseIntegrandDT

  abstract interface
     subroutine integrator_interface(this, dt, multi_step)
       import :: NewProcessDT, rkind
       class(NewProcessDT), intent(inout) :: this
       real(rkind)        , intent(in)    :: dt
       logical            , intent(in)    :: multi_step
     end subroutine integrator_interface
  end interface
  
end module BaseIntegrandM
