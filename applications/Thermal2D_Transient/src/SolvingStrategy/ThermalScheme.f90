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
    logical                                              :: firstTriang = .true.
    logical                                              :: firstQuad = .true.
    integer(ikind)                                       :: iElem, iGauss, nNode, nTriang, nQuad
    integer(ikind)                                       :: triangCounter, quadCounter
    integer(ikind)                                       :: triangPointCounter, quadPointCounter
    integer(ikind)                                       :: nElem, nGauss, nPointsTriang, nPointsQuad
    real(rkind)          , dimension(:,:,:), allocatable :: localResultMat
    type(ElementPtrDT)                                   :: element
    type(IntegratorPtrDT)                                :: integrator
    nElem = model%getnElement()
    nPointsTriang = 0
    nTriang = 0
    nPointsQuad = 0
    nQuad = 0
    do iElem = 1, nElem
       element = model%getElement(iElem)
       nNode = element%getnNode()
       integrator = element%getIntegrator()
       nGauss = integrator%getIntegTerms()
       if(nNode == 3 .or. nNode == 6) then
          nPointsTriang = nPointsTriang + nGauss
          nTriang = nTriang + 1
          if(firstTriang) then
             model%heatFlux%triangGPoint = integrator%getGPointFull()
             firstTriang = .false.
          end if
       else if(nNode == 4 .or. nNode == 8) then
          nPointsQuad = nPointsQuad + nGauss
          nQuad = nQuad + 1
          if(firstQuad) then
             model%heatFlux%quadGPoint = integrator%getGPointFull()
             firstQuad = .false.
          end if
       end if
    end do
!!$    if (allocated(model%heatFlux%triangElemID)) deallocate(model%heatFlux%triangElemID)
!!$    allocate(model%heatFlux%triangElemID(nTriang))
!!$    if (allocated(model%heatFlux%triangFlux)) deallocate(model%heatFlux%triangFlux)
!!$    allocate(model%heatFlux%triangFlux(nPointsTriang,2))
!!$    if (allocated(model%heatFlux%quadElemID)) deallocate(model%heatFlux%quadElemID)
!!$    allocate(model%heatFlux%quadElemID(nQuad))
!!$    if (allocated(model%heatFlux%quadFlux)) deallocate(model%heatFlux%quadFlux)
!!$    allocate(model%heatFlux%quadFlux(nPointsQuad,2))
    if(.not.allocated(model%heatFlux%triangElemID)) allocate(model%heatFlux%triangElemID(nTriang))
    if(.not.allocated(model%heatFlux%triangFlux)) allocate(model%heatFlux%triangFlux(nPointsTriang,2))
    if(.not.allocated(model%heatFlux%quadElemID)) allocate(model%heatFlux%quadElemID(nQuad))
    if(.not.allocated(model%heatFlux%quadFlux)) allocate(model%heatFlux%quadFlux(nPointsQuad,2))
    triangCounter = 0
    triangPointCounter = 0
    quadCounter = 0
    quadPointCounter = 0
    do iElem = 1, nElem
       element = model%getElement(iElem)
       nNode = element%getnNode()
       call element%calculateResults(model%processInfo, localResultMat)
       nGauss = size(localResultMat,2)
       if(nNode == 3 .or. nNode == 6) then
          triangCounter = triangCounter + 1
          model%heatFlux%triangElemID(triangCounter) = iElem
          do iGauss = 1, nGauss
             triangPointCounter = triangPointCounter + 1
             model%heatFlux%triangFlux(triangPointCounter,1) = localResultMat(1,iGauss,1)
             model%heatFlux%triangFlux(triangPointCounter,2) = localResultMat(1,iGauss,2)
          end do
       else if(nNode == 4 .or. nNode == 8) then
          quadCounter = quadCounter + 1
          model%heatFlux%quadElemID(quadCounter) = iElem
          do iGauss = 1, nGauss
             quadPointCounter = quadPointCounter + 1
             model%heatFlux%quadFlux(quadPointCounter,1) = localResultMat(1,iGauss,1)
             model%heatFlux%quadFlux(quadPointCounter,2) = localResultMat(1,iGauss,2)
          end do
       end if
       deallocate(localResultMat)
    end do
  end subroutine calculateFlux

  subroutine integrator(this, dt, multi_step)
    implicit none
    class(NewProcessDT), intent(inout) :: this
    real(rkind)        , intent(in)    :: dt
    logical            , intent(in) :: multi_step
  end subroutine integrator
  
end module ThermalSchemeM
