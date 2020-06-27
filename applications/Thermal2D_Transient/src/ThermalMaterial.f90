module ThermalMaterialM
  
  use UtilitiesM
  use PropertyM
  
  implicit none
  
  private
  public :: ThermalMaterialDT, thermalMaterial
  
  type :: ThermalMaterialDT
     real(rkind), dimension(2) :: conductivity
     real(rkind)               :: cp
     real(rkind)               :: rho
   contains
     procedure :: init
  end type ThermalMaterialDT
  
  interface thermalMaterial
     procedure :: constructor
  end interface thermalMaterial
  
contains
  
  type(ThermalMaterialDT) function constructor(kx,ky, cp, rho)
    implicit none
    real(rkind), intent(in) :: kx, ky, cp, rho
    call constructor%init(kx, ky, cp, rho)
  end function constructor
  
  subroutine init(this, kx, ky, cp, rho)
    implicit none
    real(rkind), intent(in) :: kx, ky, cp, rho
    class(ThermalMaterialDT), intent(inout) :: this
    this%conductivity = (/kx, ky/)
    this%cp           = cp
    this%rho          = rho
  end subroutine init
  
end module ThermalMaterialM
