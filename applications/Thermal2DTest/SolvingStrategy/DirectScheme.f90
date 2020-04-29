module DirectSchemeM

  use UtilitiesM
  use DebuggerM

  use PointPtrM

  use IntegratorPtrM

  use Element1DPtrM
  use Element2DPtrM
  
  use ThermalmodelM
  
  use SchemeM

  implicit none

  private
  public :: DirectSchemeDT

  type, extends(NewSchemeDT) :: DirectSchemeDT
   contains
     procedure, public              :: calculateFlux
     procedure                      :: valueGPoints1D
     procedure                      :: valueGPointsTriang
     procedure                      :: valueGPointsQuad
  end type DirectSchemeDT

  procedure(addTriangFlux), pointer :: addFlux => null()
  integer(ikind)                    :: lineCount
  integer(ikind)                    :: triangCount
  integer(ikind)                    :: quadCount
  
contains

  subroutine calculateFlux(this, model)
    implicit none
    class(DirectSchemeDT), intent(inout)          :: this
    class(ThermalModelDT), intent(inout)          :: model
    type(Element1DPtrTYPE)                        :: element1D
    type(Element2DPtrTYPE)                        :: element2D
    type(IntegratorPtrTYPE)                       :: integrator
    type(PointPtrTYPE), dimension(:), allocatable :: point
    integer(ikind) :: iElem, iGauss, i, triangElemCount, quadElemCount
    integer(ikind) :: nElem, nPoint, nLine, nTriang, nQuad
    integer(ikind) :: nGauss, nGaussPointLine, nGaussPointTriang, nGaussPointQuad
    real(rkind) :: xi, eta, x, y,  bi, ci, dNidx, dNidy, jacobianDet
    real(rkind) :: k, kx, ky, q, qx, qy, jacobian1D
    real(rkind), dimension(:,:), allocatable :: dsf
    real(rkind), dimension(2,2) :: jacobian
    write(*,*) '***  Direct Scheme ***'
    write(*,*) '*** Calculate Flux ***'
    call this%valueGPoints1D(model)
    nLine = model%nLine
    if(nLine > 0) then
       nGauss = model%elementList1D%getGaussOrder()
    else
       nGauss = model%elementList2D%getGaussOrder()
    end if
    nGaussPointLine = getnGaussPointLine(nGauss, nLine)
    allocate(model%heatFlux%lineElemID(nLine))
    allocate(model%heatFlux%lineQ(nGaussPointLine))
    !Flux for one dimensional elements:
    lineCount = 0
    nElem = model%getnLine()
    do iElem = 1, nElem
       element1D = model%elementList1D%getElement(iElem)
       model%heatFlux%lineElemID(iElem) = element1D%getID()
       nPoint = element1D%getnPoint()
       integrator = element1D%getIntegrator()
       do iGauss = 1, integrator%ptr%integTerms
          xi = integrator%ptr%gPoint(iGauss,1)
          jacobian1D = element1D%jacobian(xi)
          q = 0.d0
          do i = 1, nPoint
             q = q + integrator%ptr%dShapeFunc(iGauss,1,i) &
                  *problem%dof(element1D%getPointID(i))/jacobian1D
          end do
          k = element1D%ptr%material%ptr%conductivity(1)
          call addLineFlux(this, -1.d0*k*q, model)
       end do
    end do
    !Flux for bidimensional elements:
    nTriang = model%nTriang
    nQuad = model%nQuad
    if(nTriang == 0 .and. nQuad == 0) return
    call this%valueGPointsTriang(model)
    call this%valueGPointsQuad(model)
    nGaussPointTriang = getnGaussPointTriang(nGauss, nTriang)
    nGaussPointQuad = getnGaussPointQuad(nGauss, nQuad)
    allocate(model%heatFlux%triangElemID(nTriang))
    allocate(model%heatFlux%triangQx(nGaussPointTriang))
    allocate(model%heatFlux%triangQy(nGaussPointTriang))
    allocate(model%heatFlux%quadElemID(nQuad))
    allocate(model%heatFlux%quadQx(nGaussPointQuad))
    allocate(model%heatFlux%quadQy(nGaussPointQuad))
    triangCount = 0
    quadCount = 0
    triangElemCount = 0
    quadElemCount = 0
    nElem = model%getnTriang() + model%getnQuad()
    do iElem = 1, nElem
       element2D = model%elementList2D%getElement(iElem)
       nPoint = element2D%getnPoint()
       integrator = element2D%getIntegrator()
       if(nPoint == 3 .or. nPoint == 6) then
          triangElemCount = triangElemCount + 1
          addFlux => addTriangFlux
          model%heatFlux%triangElemID(triangElemCount) = element2D%getID()
       else if(nPoint == 4 .or. nPoint == 8) then
          quadElemCount = quadElemCount + 1
          addFlux => addQuadFlux
          model%heatFlux%quadElemID(quadElemCount) = element2D%getID()
       else
          print'(A)', '** PostProcess ERROR1 **'
       end if
       do iGauss = 1, integrator%ptr%integTerms
          xi = integrator%ptr%gPoint(iGauss,1)
          eta = integrator%ptr%gPoint(iGauss,2)
          jacobian = element2D%jacobian(xi, eta)
          jacobianDet = element2D%jacobianDet(jacobian)
          allocate(dsf(2,nPoint))
          dsf = integrator%ptr%dShapeFunc(iGauss,:,:)
          qx = 0.d0
          qy = 0.d0
          do i = 1, nPoint
             bi = jacobian(2,2)*dsf(1,i) - jacobian(1,2)*dsf(2,i)
             ci = jacobian(1,1)*dsf(2,i) - jacobian(2,1)*dsf(1,i)
             dNidx = bi/jacobianDet
             dNidy = ci/jacobianDet
             qx = qx + dNidx*problem%dof(element2D%getPointID(i))
             qy = qy + dNidy*problem%dof(element2D%getPointID(i))
          end do
          kx = element2D%ptr%material%ptr%conductivity(1)
          ky = element2D%ptr%material%ptr%conductivity(2)
          call addFlux(this, -1.d0*kx*qx, -1.d0*ky*qy, model)
          deallocate(dsf)
       end do
    end do
    
  end subroutine calculateFlux

subroutine valueGPoints1D(this, model)
    implicit none
    class(DirectSchemeDT), intent(inout) :: this
    class(ThermalmodelDT), intent(inout) :: model
    type(IntegratorPtrDT) :: integrator
    integrator = model%elementList1D%getIntegrator()
    if(allocated(integrator%ptr%gPoint)) then
       allocate(model%heatFlux%lineGPoint(size(integrator%ptr%gPoint,1)))
       model%heatFlux%lineGPoint(:) = integrator%ptr%gPoint(:,1)
    end if
  end subroutine valueGPoints1D

  subroutine valueGPointsTriang(this, model)
    implicit none
    class(DirectSchemeDT), intent(inout) :: this
    class(ThermalModelDT), intent(inout) :: model
    type(IntegratorPtrDT) :: integrator
    integrator = model%elementList2D%getTriangIntegrator()
    if(allocated(integrator%ptr%gPoint)) then
       allocate(model%heatFlux%triangGPoint(size(integrator%ptr%gPoint,1),size(integrator%ptr%gPoint,2)))
       model%heatFlux%triangGPoint = integrator%ptr%gPoint
    end if
  end subroutine valueGPointsTriang

  subroutine valueGPointsQuad(this, model)
    implicit none
    class(DirectSchemeDT), intent(inout) :: this
    class(ThermalModelDT), intent(inout) :: model
    type(IntegratorPtrDT) :: integrator
    integrator = model%elementList2D%getQuadIntegrator()
    if(allocated(integrator%ptr%gPoint)) then
       allocate(model%heatFlux%quadGPoint(size(integrator%ptr%gPoint,1),size(integrator%ptr%gPoint,2)))
       model%heatFlux%quadGPoint = integrator%ptr%gPoint
    end if
  end subroutine valueGPointsQuad

  subroutine addLineFlux(this, q, model)
    implicit none
    class(DirectSchemeDT), intent(inout) :: this
    class(ThermalModelDT), intent(inout) :: model
    real(rkind), intent(in) :: q
    lineCount = lineCount + 1
    model%heatFlux%lineQ(lineCount) = q
  end subroutine addLineFlux

  subroutine addTriangFlux(this, qx, qy, model)
    implicit none
    class(DirectSchemeDT)   , intent(inout) :: this
    class(ThermalStrategyDT), intent(inout) :: model
    real(rkind), intent(in) :: qx
    real(rkind), intent(in) :: qy
    triangCount = triangCount + 1
    model%heatFlux%triangQx(triangCount) = qx
    model%heatFlux%triangQy(triangCount) = qy
  end subroutine addTriangFlux

  subroutine addQuadFlux(this, qx, qy, model)
    implicit none
    class(DirectSchemeDT), intent(inout) :: this
    class(ThermalModelDT), intent(inout) :: model
    real(rkind)          , intent(in)    :: qx
    real(rkind)          , intent(in)    :: qy
    quadCount = quadCount + 1
    model%heatFlux%quadQx(quadCount) = qx
    model%heatFlux%quadQy(quadCount) = qy
  end subroutine addQuadFlux

  integer(ikind) function getnGaussPointLine(nGauss, nLine)
    implicit none
    integer(ikind), intent(in) :: nGauss
    integer(ikind), intent(in) :: nLine
    getnGaussPointLine = nLine*nGauss
  end function getnGaussPointLine

  integer(ikind) function getnGaussPointTriang(nGauss, nTriang)
    implicit none
    integer(ikind), intent(in) :: nGauss
    integer(ikind), intent(in) :: nTriang
    if(nGauss == 1) then
       getnGaussPointTriang = nTriang
    else if(nGauss == 2) then
       getnGaussPointTriang = 3*nTriang
    else if(nGauss == 3) then
       getnGaussPointTriang = 4*nTriang
    else
       print*, '** Input Gauss Order not supported! **'
    end if
  end function getnGaussPointTriang
  
  integer(ikind) function getnGaussPointQuad(nGauss, nQuad)
    implicit none
    integer(ikind), intent(in) :: nGauss
    integer(ikind), intent(in) :: nQuad
    if(nGauss == 1) then
       getnGaussPointQuad = nQuad
    else if(nGauss == 2) then
       getnGaussPointQuad = 4*nQuad
    else if(nGauss == 3) then
       getnGaussPointQuad = 9*nQuad
    else
       print*, '** Input Gauss Order not supported! **'
    end if
  end function getnGaussPointQuad
  
end module DirectSchemeM
