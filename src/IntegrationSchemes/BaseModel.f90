module BaseModelM

  use UtilitiesM
  
  use SchemeM
  
  use IntegrandM
  
  implicit none
  
  private
  public :: SetBaseModel, BaseModelDT
  
  type, extends(IntegrandDT) :: BaseModelDT
     real(rkind) :: constant
   contains
     procedure, public :: t          => dBaseModel_dt 
     procedure, public :: add        => add_BaseModel
     procedure, public :: multiply   => multiply_BaseModel
     procedure, public :: assign     => assign_BaseModel
     procedure, public :: useProcess => process
     procedure, public :: output    
  end type BaseModelDT
  
  interface SetBaseModel
     procedure constructor
  end interface SetBaseModel
  
contains

  subroutine process(this)
    implicit none
    class(BaseModelDT), intent(inout) :: this
  end subroutine process
  
  type(BaseModelDT) function constructor(initial_state,c,integrand)
    real(rkind), dimension(:) ,intent(in) :: initial_state
    real(rkind), intent(in) :: c
    class(NewSchemeDT),intent(in) :: integrand
    constructor%state=initial_state
    constructor%constant=c
    call constructor%set_quadrature(integrand)
  end function constructor
  
  function dBaseModel_dt(this, dof) result(dState_dt)
    class(BaseModelDT),intent(in) :: this
    class(IntegrandDT) ,allocatable :: dState_dt
    type(BaseModelDT),allocatable :: local_dState_dt
    real(rkind), dimension(:), intent(in) :: dof
    allocate(local_dState_dt)
    call local_dState_dt%set_quadrature(this%get_quadrature())
    allocate(local_dState_dt%state(size(this%state)))
    local_dState_dt%state(1) = this%state(2)-this%state(1)
    local_dState_dt%state(2) = this%state(1)*this%state(3)
    local_dState_dt%state(3) = this%state(1)-this%state(2)
    local_dState_dt%constant = 1.
    call move_alloc(local_dState_dt,dState_dt)
  end function dBaseModel_dt
  
  function add_BaseModel(lhs,rhs) result(sum)
    class(BaseModelDT),intent(in) :: lhs
    class(IntegrandDT),intent(in) :: rhs
    class(IntegrandDT),allocatable :: sum
    type(BaseModelDT),allocatable :: local_sum
    select type(rhs)
    class is (BaseModelDT)
       allocate(local_sum)
       call local_sum%set_quadrature(lhs%get_quadrature())
       local_sum%state = lhs%state + rhs%state
       local_sum%constant = lhs%constant + rhs%constant
    class default
       stop 'assig_BaseModel: unsupported class'
    end select
    call move_alloc(local_sum,sum)
  end function add_BaseModel
  
  function multiply_BaseModel(lhs,rhs) result(product)
    class(BaseModelDT),intent(in)   :: lhs
    real(rkind)       ,intent(in)   :: rhs
    class(IntegrandDT), allocatable :: product
    type(BaseModelDT) , allocatable :: local_product
    allocate(local_product)
    call local_product%set_quadrature(lhs%get_quadrature())
    local_product%state = rhs * lhs%state
    local_product%constant = rhs * lhs%constant
    call move_alloc(local_product,product)
  end function multiply_BaseModel
  
  subroutine assign_BaseModel(lhs,rhs)
    class(BaseModelDT),intent(inout) :: lhs
    class(IntegrandDT),intent(in) :: rhs
    select type(rhs)
    class is (BaseModelDT)
       lhs%state = rhs%state
       lhs%constant = rhs%constant
       call lhs%set_quadrature(rhs%get_quadrature())
       class default
       stop 'assign_BaseModel: unsupported class'
    end select
  end subroutine assign_BaseModel
  
  function output(this) result(coordinates)
    class(BaseModelDT),intent(in) :: this
    real(rkind),dimension(:) ,allocatable :: coordinates
    coordinates = [ this%constant, this%state ]
  end function output
  
end module BaseModelM
