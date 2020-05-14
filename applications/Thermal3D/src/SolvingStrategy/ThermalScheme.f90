module ThermalSchemeM

  use UtilitiesM
  use DebuggerM

  use IntegratorPtrM

  use ElementPtrM
  
  use ThermalmodelM
  
  use SchemeM

  use ProcessM

  implicit none

  private
  public :: ThermalSchemeDT

  type, extends(NewSchemeDT) :: ThermalSchemeDT
   contains
     procedure, public              :: calculateFlux
     procedure, nopass              :: integrate => integrator
  end type ThermalSchemeDT
  
contains

  subroutine calculateFlux(this, model)
    implicit none
    class(ThermalSchemeDT), intent(inout)                :: this
    class(ThermalModelDT), intent(inout)                 :: model
    logical                                              :: firstTetra = .true.
    logical                                              :: firstHexa = .true.
    integer(ikind)                                       :: iElem, iGauss, nNode, nTetra, nHexa
    integer(ikind)                                       :: tetraCounter, hexaCounter
    integer(ikind)                                       :: tetraPointCounter, hexaPointCounter
    integer(ikind)                                       :: nElem, nGauss, nPointsTetra, nPointsHexa
    real(rkind)          , dimension(:,:,:), allocatable :: localResultMat
    type(ElementPtrDT)                                   :: element
    type(IntegratorPtrDT)                                :: integrator
    write(*,*) '***  Direct Scheme ***'
    write(*,*) '*** Calculate Flux ***'
    nElem = model%getnElement()
    nPointsTetra = 0
    nTetra = 0
    nPointsHexa = 0
    nHexa = 0
    do iElem = 1, nElem
       element = model%getElement(iElem)
       nNode = element%getnNode()
       integrator = element%getIntegrator()
       nGauss = integrator%getIntegTerms()
       if(nNode == 4) then
          nPointsTetra = nPointsTetra + nGauss
          nTetra = nTetra + 1
          if(firstTetra) then
             model%heatFlux%tetraGPoint = integrator%getGPointFull()
             firstTetra = .false.
          end if
       else if(nNode == 8) then
          nPointsHexa = nPointsHexa + nGauss
          nHexa = nHexa + 1
          if(firstHexa) then
             model%heatFlux%hexaGPoint = integrator%getGPointFull()
             firstHexa = .false.
          end if
       end if
    end do
    allocate(model%heatFlux%tetraElemID(nTetra))
    allocate(model%heatFlux%tetraFlux(nPointsTetra,3))
    allocate(model%heatFlux%hexaElemID(nHexa))
    allocate(model%heatFlux%hexaFlux(nPointsHexa,3))
    tetraCounter = 0
    tetraPointCounter = 0
    hexaCounter = 0
    hexaPointCounter = 0
    do iElem = 1, nElem
       element = model%getElement(iElem)
       nNode = element%getnNode()
       call element%calculateResults(localResultMat)
       nGauss = size(localResultMat,2)
       if(nNode == 4) then
          tetraCounter = tetraCounter + 1
          model%heatFlux%tetraElemID(tetraCounter) = iElem
          do iGauss = 1, nGauss
             tetraPointCounter = tetraPointCounter + 1
             model%heatFlux%tetraFlux(tetraPointCounter,1) = localResultMat(1,iGauss,1)
             model%heatFlux%tetraFlux(tetraPointCounter,2) = localResultMat(1,iGauss,2)
             model%heatFlux%tetraFlux(tetraPointCounter,3) = localResultMat(1,iGauss,3)
          end do
       else if(nNode == 8) then
          hexaCounter = hexaCounter + 1
          model%heatFlux%hexaElemID(hexaCounter) = iElem
          do iGauss = 1, nGauss
             hexaPointCounter = hexaPointCounter + 1
             model%heatFlux%hexaFlux(hexaPointCounter,1) = localResultMat(1,iGauss,1)
             model%heatFlux%hexaFlux(hexaPointCounter,2) = localResultMat(1,iGauss,2)
             model%heatFlux%hexaFlux(hexaPointCounter,3) = localResultMat(1,iGauss,3)
          end do
       end if
       deallocate(localResultMat)
    end do
  end subroutine calculateFlux

  subroutine integrator(this, dt)
    implicit none
    class(NewProcessDT), intent(inout) :: this
    real(rkind)        , intent(in)    :: dt
  end subroutine integrator
  
end module ThermalSchemeM
