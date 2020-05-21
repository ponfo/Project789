module NavierStokes2DM

  use UtilitiesM
  use SparseKit
  
  use SchemeM

  use IntegrandM

  implicit none

  private
  public :: NavierStokes2DDT, SetNavierStokes2D

  type, extends(IntegrandDT) :: NavierStokes2DDT
     type(Sparse)                           :: lhs
     type(Sparse)                           :: lumpedMassInverse
     real(rkind), dimension(:), allocatable :: rhs
   contains
     procedure :: t          => dnavierStokes2D_dt 
     procedure :: add        => add_navierStokes2D
     procedure :: multiply   => multiply_navierStokes2D
     procedure :: assign     => assign_navierStokes2D
     procedure :: getState
     procedure :: useProcess => process
  end type NavierStokes2DDT

  interface SetNavierStokes2D
     procedure constructor
  end interface SetNavierStokes2D

contains

  type(navierStokes2DDT) function constructor(initial_state, lhs&
       , rhs, lumpedMassInverse, this_strategy, step)
    class(NewSchemeDT)       , intent(in) :: this_strategy
    type(Sparse)             , intent(in) :: lhs
    type(Sparse)             , intent(in) :: lumpedMassInverse
    real(rkind), dimension(:), intent(in) :: initial_state
    real(rkind), dimension(:), intent(in) :: rhs
    integer(ikind), intent(in), optional :: step
    constructor%state             = initial_state
    constructor%lhs               = lhs
    constructor%rhs               = rhs
    constructor%lumpedMassInverse = lumpedMassInverse
    if (present(step)) constructor%step = step
    call constructor%set_quadrature(this_strategy)
  end function constructor
  
  function dNavierStokes2D_dt(this, dof) result(dState_dt)
    class(NavierStokes2DDT), intent(in)   :: this
    class(IntegrandDT), allocatable       :: dState_dt
    type(NavierStokes2DDT) , allocatable  :: local_dState_dt
    real(rkind), dimension(:), intent(in) :: dof
    integer(ikind) :: n
    allocate(local_dState_dt)
    call local_dState_dt%set_quadrature(this%get_quadrature())
    allocate(local_dState_dt%state(size(this%state)))
    local_dState_dt%state = (this%lumpedMassInverse)*(this%rhs-this%lhs*dof)
    call move_alloc(local_dState_dt,dState_dt)
  end function dNavierStokes2D_dt
  
  function add_navierStokes2D(lhs,rhs) result(sum)
    class(NavierStokes2DDT), intent(in)  :: lhs
    class(IntegrandDT), intent(in)  :: rhs
    class(IntegrandDT), allocatable :: sum
    type(NavierStokes2DDT) , allocatable :: local_sum
    select type(rhs)
    class is (navierStokes2DDT)
       allocate(local_sum)
       call local_sum%set_quadrature(lhs%get_quadrature())
       local_sum%state = lhs%state + rhs%state
    class default
       stop 'assig_navierStokes2D: unsupported class'
    end select
    call move_alloc(local_sum,sum)
  end function add_navierStokes2D

  function multiply_navierStokes2D(lhs,rhs) result(product)
    class(NavierStokes2DDT), intent(in)  :: lhs
    class(IntegrandDT), allocatable :: product
    type(NavierStokes2DDT) , allocatable :: local_product
    real(rkind)       , intent(in)  :: rhs
    allocate(local_product)
    call local_product%set_quadrature(lhs%get_quadrature())
    local_product%state = lhs%state * rhs            
    call move_alloc(local_product,product)
  end function multiply_navierStokes2D

  subroutine assign_navierStokes2D(lhs,rhs)
    class(NavierStokes2DDT), intent(inout) :: lhs
    class(IntegrandDT), intent(in)    :: rhs
    select type(rhs)
    class is (NavierStokes2DDT)
       lhs%state = rhs%state
       lhs%step = rhs%step
       lhs%rhs = rhs%rhs
       lhs%lhs = rhs%lhs
       lhs%lumpedMassInverse = rhs%lumpedMassInverse
       call lhs%set_quadrature(rhs%get_quadrature())
    class default
       stop 'assign_navierStokes2D: unsupported class'
    end select
  end subroutine assign_navierStokes2D
  
  function getState(this) result(coordinates)
    class(NavierStokes2DDT)       , intent(in)  :: this
    real(rkind), dimension(:), allocatable :: coordinates
    coordinates = this%state
  end function getState

  subroutine process(this)
    implicit none
    class(NavierStokes2DDT), intent(inout) :: this
  end subroutine process
  
end module NavierStokes2DM
