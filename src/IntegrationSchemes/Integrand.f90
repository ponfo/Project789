module IntegrandM

  use UtilitiesM
  
  use ProcessM

  use SchemeM 

  implicit none 

  private
  public :: IntegrandDT
  
  type, abstract, extends(NewProcessDT) :: IntegrandDT
     class(NewSchemeDT), allocatable :: quadrature
   contains
     procedure, non_overridable :: integrate  
     procedure, non_overridable :: set_quadrature
     procedure, non_overridable :: get_quadrature
     procedure(time_derivative), deferred :: t 
     procedure(symmetric_operator), deferred :: add
     procedure(asymmetric_operator), deferred :: multiply
     procedure(symmetric_assignment), deferred :: assign
     generic :: operator(+) => add
     generic :: operator(*) => multiply
     generic :: assignment(=) => assign
  end type IntegrandDT
  
  abstract interface


     
     function time_derivative(this) result(dState_dt)
       import :: IntegrandDT
       class(IntegrandDT) ,intent(in) :: this
       class(IntegrandDT) ,allocatable :: dState_dt
     end function time_derivative
     
     function symmetric_operator(lhs,rhs) result(operator_result)
       import :: IntegrandDT
       class(IntegrandDT) ,intent(in) :: lhs,rhs
       class(IntegrandDT) ,allocatable :: operator_result
     end function symmetric_operator
     
     function asymmetric_operator(lhs,rhs) result(operator_result)  
       import :: IntegrandDT, rkind
       class(IntegrandDT) ,intent(in) :: lhs
       class(IntegrandDT) ,allocatable :: operator_result
       real(rkind),intent(in) :: rhs
     end function asymmetric_operator
     
     subroutine symmetric_assignment(lhs,rhs)
       import :: IntegrandDT
       class(IntegrandDT) ,intent(in) :: rhs
       class(IntegrandDT) ,intent(inout) :: lhs
     end subroutine symmetric_assignment
     
  end interface
  
contains
  
  subroutine set_quadrature (this, s)
    class(IntegrandDT), intent(inout) :: this
    class(NewSchemeDT), intent(in) :: s
    if (allocated(this%quadrature)) deallocate (this%quadrature)
    allocate (this%quadrature, source=s)
  end subroutine set_quadrature
  
  function get_quadrature (this) result (this_strategy)
    class(IntegrandDT), intent(in) :: this
    class(NewSchemeDT), allocatable :: this_strategy
    allocate (this_strategy, source=this%quadrature)
  end function get_quadrature
  
  subroutine integrate(model,dt)
    class(IntegrandDT) :: model
    real(rkind) ,intent(in) :: dt
    if (allocated(model%quadrature)) then
       call model%quadrature%integrate(model,dt)
    else
       stop 'integrate: no integration procedure available.'
    end if
  end subroutine integrate
  
end module IntegrandM