module Poisson2DM

  use UtilitiesM
  use SparseKit
  
  use BaseIntegrandM

  use IntegrandM

  implicit none

  private
  public :: Poisson2DDT, SetPoisson2D

  type, extends(IntegrandDT) :: Poisson2DDT
     type(Sparse)                           :: stiffness
     type(Sparse)                           :: lumpedMassInverse
     real(rkind), dimension(:), allocatable :: rhs
   contains
     procedure :: t          => dpoisson2D_dt 
     procedure :: add        => add_poisson2D
     procedure :: multiply   => multiply_poisson2D
     procedure :: assign     => assign_poisson2D
     procedure :: getState
  end type Poisson2DDT

  interface SetPoisson2D
     procedure constructor
  end interface SetPoisson2D

contains

  type(poisson2DDT) function constructor(initial_state, stiffness&
       , rhs, lumpedMassInverse, this_strategy, step)
    class(BaseIntegrandDT)   , intent(in) :: this_strategy
    type(Sparse)             , intent(in) :: stiffness
    type(Sparse)             , intent(in) :: lumpedMassInverse
    real(rkind), dimension(:), intent(in) :: initial_state
    real(rkind), dimension(:), intent(in) :: rhs
    integer(ikind), intent(in), optional :: step
    constructor%state             = initial_state
    constructor%stiffness         = stiffness
    constructor%rhs               = rhs
    constructor%lumpedMassInverse = lumpedMassInverse
    if (present(step)) constructor%step = step
    call constructor%set_quadrature(this_strategy)
  end function constructor
  
  function dPoisson2D_dt(this, dof) result(dState_dt)
    class(Poisson2DDT), intent(in)                     :: this
    class(IntegrandDT), allocatable                    :: dState_dt
    type(Poisson2DDT) , allocatable                    :: local_dState_dt
    real(rkind), dimension(:), intent(in) :: dof
    integer(ikind) :: n
    allocate(local_dState_dt)
    call local_dState_dt%set_quadrature(this%get_quadrature())
    allocate(local_dState_dt%state(size(this%state)))
    local_dState_dt%state = (this%lumpedMassInverse)*(this%rhs-this%stiffness*dof)
    call move_alloc(local_dState_dt,dState_dt)
  end function dPoisson2D_dt
  
  function add_poisson2D(lhs,rhs) result(sum)
    class(Poisson2DDT), intent(in)  :: lhs
    class(IntegrandDT), intent(in)  :: rhs
    class(IntegrandDT), allocatable :: sum
    type(Poisson2DDT) , allocatable :: local_sum
    select type(rhs)
    class is (poisson2DDT)
       allocate(local_sum)
       call local_sum%set_quadrature(lhs%get_quadrature())
       local_sum%state = lhs%state + rhs%state
    class default
       stop 'assig_poisson2D: unsupported class'
    end select
    call move_alloc(local_sum,sum)
  end function add_poisson2D

  function multiply_poisson2D(lhs,rhs) result(product)
    class(Poisson2DDT), intent(in)  :: lhs
    class(IntegrandDT), allocatable :: product
    type(Poisson2DDT) , allocatable :: local_product
    real(rkind)       , intent(in)  :: rhs
    allocate(local_product)
    call local_product%set_quadrature(lhs%get_quadrature())
    local_product%state = lhs%state * rhs            
    call move_alloc(local_product,product)
  end function multiply_poisson2D

  subroutine assign_poisson2D(lhs,rhs)
    class(Poisson2DDT), intent(inout) :: lhs
    class(IntegrandDT), intent(in)    :: rhs
    select type(rhs)
    class is (Poisson2DDT)
       lhs%state = rhs%state
       lhs%step = rhs%step
       lhs%rhs = rhs%rhs
       lhs%stiffness = rhs%stiffness
       lhs%lumpedMassInverse = rhs%lumpedMassInverse
       call lhs%set_quadrature(rhs%get_quadrature())
    class default
       stop 'assign_poisson2D: unsupported class'
    end select
  end subroutine assign_poisson2D
  
  function getState(this) result(coordinates)
    class(Poisson2DDT)       , intent(in)  :: this
    real(rkind), dimension(:), allocatable :: coordinates
    coordinates = this%state
  end function getState
  
end module Poisson2DM
