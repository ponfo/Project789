module CFDMaterialM
  
  use UtilitiesM
  use PropertyM
  
  implicit none
  
  private
  public :: CFDMaterialDT, cfdMaterial
  
  type :: CFDMaterialDT
     real(rkind) :: R
     real(rkind) :: gamma
     real(rkind) :: mu
     real(rkind) :: k
     real(rkind) :: Vx_inf
     real(rkind) :: Vy_inf
     real(rkind) :: T_inf
     real(rkind) :: rho
     real(rkind) :: P_inf
     real(rkind) :: Cv
     real(rkind) :: mach
     real(rkind) :: Vc
   contains
     procedure :: init
  end type CFDMaterialDT
  
  interface cfdMaterial
     procedure :: constructor
  end interface cfdMaterial
  
contains
  
  type(CFDMaterialDT) function constructor(R, gamma, mu, k, Vx, Vy, T, rho, mach, Cv, P, Vc)
    implicit none
    real(rkind), intent(in) :: R
    real(rkind), intent(in) :: gamma
    real(rkind), intent(in) :: mu
    real(rkind), intent(in) :: k
    real(rkind), intent(in) :: Vx
    real(rkind), intent(in) :: Vy
    real(rkind), intent(in) :: T
    real(rkind), intent(in) :: rho
    real(rkind), intent(in) :: mach
    real(rkind), intent(in) :: Cv
    real(rkind), intent(in) :: P
    real(rkind), intent(in) :: Vc
    call constructor%init(R, gamma, mu, k, Vx, Vy, T, rho, mach, Cv, P, Vc)
  end function constructor

  subroutine init(this, R, gamma, mu, k, Vx, Vy, T, rho, mach, Cv, P, Vc)
    implicit none
    class(CFDMaterialDT), intent(inout) :: this
    real(rkind)         , intent(in)    :: R
    real(rkind)         , intent(in)    :: gamma
    real(rkind)         , intent(in)    :: mu
    real(rkind)         , intent(in)    :: k
    real(rkind)         , intent(in)    :: Vx
    real(rkind)         , intent(in)    :: Vy
    real(rkind)         , intent(in)    :: T
    real(rkind)         , intent(in)    :: rho
    real(rkind)         , intent(in)    :: mach
    real(rkind)         , intent(in)    :: Cv
    real(rkind)         , intent(in)    :: P
    real(rkind)         , intent(in)    :: Vc
    this%R      = R
    this%gamma  = gamma
    this%mu     = mu
    this%k      = k
    this%Vx_inf = Vx
    this%Vy_inf = Vy
    this%T_inf  = T
    this%rho    = rho
    this%mach   = mach
    this%P_inf  = P
    this%Cv     = Cv
    this%Vc     = Vc
  end subroutine init
  
end module CFDMaterialM
