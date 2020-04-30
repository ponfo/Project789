module ThermalElement2DM

  use ElementM

  implicit none

  private
  public :: ThermalElementDT, thermalElement, initGeometries

  type, extends(ElementDT) :: ThermalElementDT
     type(ThermalMaterialDT), pointer :: material
   contains
     procedure, public :: init
  end type ThermalElementDT

  interface thermalElement2D
     procedure :: constructor
  end interface thermalElement2D

  type(Triangle2D3NodeDT)     , save :: triangle2D3Node
  type(Triangle2D6NodeDT)     , save :: triangle2D6Node
  type(Quadrilateral2D4NodeDT), save :: Quadrilateral2D4Node
  type(Quadrilateral2D8NodeDT), save :: Quadrilateral2D8Node

contains

  type(ThermalElementDT) function contructor(id, node, material)
    implicit none
    integer(ikind)                       , intent(in) :: id
    type(NodePtrDT)        , dimension(:), intent(in) :: node
    type(ThermalMaterialDT), target      , intent(in) :: material
    call this%init(id, nGauss, node, material)
  end function contructor

  subroutine init(this, id, node, material)
    implicit none
    class(ThermalElementDT)              , intent(inout) :: this
    integer(ikind)                       , intent(in)    :: id
    type(NodePtrDT)        , dimension(:), intent(in)    :: node
    type(ThermalMaterialDT), target      , intent(in)    :: material
    this%id = id
    this%node = node
    this%material = material
    if(size(node) == 3) then
       this%geometry => triangle2D3Node
    else if(size(node) == 4) then
       this%geometry => quadrilateral2D4Node
    else if(size(node) == 6) then
       this%geometry => triangle2D6Node
    else if(size(node) == 8) then
       this%geometry => quadrilateral2D8Node
    end if
  end subroutine init

  subroutine initGeometries(nGauss)
    implicit none
    integer(ikind), intent(in) :: nGauss
    triangle2D3Node = triangle2D3Node(nGauss)
    triangle2D6Node = triangle2D6Node(nGauss)
    quadrilateral2D4Node = quadrilateral2D4Node(nGauss)
    quadrilateral2D8Node = quadrilateral2D8Node(nGauss)
  end subroutine initGeometries

  subroutine calculateLocalSystem(this, lhs, rhs)
    implicit none
    class(ThermalElementDT)                             , intent(inout) :: this
    real(rkind)            , dimension(:,:), allocatable, intent(inout) :: lhs
    real(rkind)            , dimension(:)  , allocatable, intent(inout) :: rhs
    call this%calculateLHS(lhs)
    call this%calculateRHS(rhs)
  end subroutine calculateLocalSystem

  subroutine calculateLHS(this, lhs)
    implicit none
    class(ThermalElementDT)                               , intent(inout) :: this
    real(rkind)            , dimension(:,:)  , allocatable, intent(inout) :: lhs
    integer(ikind)                                                        :: i, j, k, nNode
    real(rkind)                                                           :: bi, bj, ci, cj
    real(rkind)            , dimension(:,:,:), allocatable                :: jacobian
    real(rkind)            , dimension(:)    , allocatable                :: jacobianDet
    type(IntegratorPtrDT)                                                 :: integrator
    allocate(lhs(this%nNode,this%nNode))
    integrator%ptr => this%geometry%integrator
    allocate(jacobian(integrator%ptr%integTerms,2,2))
    allocate(jacobianDet(integrator%ptr%integTerms))
    do i = 1, integrator%ptr%integTerms
       jacobian(i,1:2,1:2) = this%jacobian(integrator%ptr%gPoint(i,1),integrator%ptr%gPoint(i,2))
       jacobianDet(i) = this%jacobianDet(jacobian(i,1:2,1:2))
    end do
    nNode = this%geometry%nNode
    do i = 1, nNode
       do j = 1, nNode
          lhs(i,j) = 0.d0
          do k = 1, integrator%ptr%integTerms
             bi = jacobian(k,2,2)*integrator%ptr%dShapeFunc(k,1,i) &
                  - jacobian(k,1,2)*integrator%ptr%dShapeFunc(k,2,i)
             bj = jacobian(k,2,2)*integrator%ptr%dShapeFunc(k,1,j) &
                  - jacobian(k,1,2)*integrator%ptr%dShapeFunc(k,2,j)
             ci = jacobian(k,1,1)*integrator%ptr%dShapeFunc(k,2,i) &
                  - jacobian(k,2,1)*integrator%ptr%dShapeFunc(k,1,i)
             cj = jacobian(k,1,1)*integrator%ptr%dShapeFunc(k,2,j) &
                  - jacobian(k,2,1)*integrator%ptr%dShapeFunc(k,1,j)
             lhs(i,j) = lhs(i,j)                            &
                  + integrator%ptr%weight(k)                 &
                  *(this%material%ptr%conductivity(1)*bi*bj  &
                  + this%material%ptr%conductivity(2)*ci*cj) &
                  / jacobianDet(k)
          end do
       end do
    end do
  end subroutine calculateLHS

  subroutine calculateRHS(this, rhs)
    implicit none
    class(ThermalElementDT)                           , intent(inout) :: this
    real(rkind)            , dimension(:), allocatable, intent(inout) :: rhs
    integer(ikind)                                                    :: i, j, nNode
    real(rkind)                                                       :: val
    type(IntegratorPtrDT)                                             :: integrator
    allocate(rhs(this%nNode))
    rhs = 0.d0
    do i = 1, this%nNode
       if(allocated(this%node(i)%ptr%source)) then
          val = this%node(i)%ptr%source%func(1)%evaluate((/this%node(i)%getx(), this%node(i)%gety()/))
          rhs(i) = rhs(i) + val
       end if
    end do
    if(allocated(this%source)) then
       integrator%ptr => this%geometry%integrator
       allocate(valuedSource(integrator%ptr%integTerms))
       allocate(jacobianDet(integrator%ptr%integTerms))
       call this%setupIntegration(integrator, valuedSource, jacobianDet)
       nNode = this%nNode
       do i = 1, nNode
          val = 0
          do j = 1, integrator%ptr%integTerms
             val = val + integrator%ptr%weight(j)*integrator%ptr%shapeFunc(j,i) &
                  *valuedSource(j)*jacobianDet(j)
          end do
          rhs(i) = rhs(i) + val
       end do
       deallocate(valuedSource)
       deallocate(jacobianDet)
    end if
  end subroutine calculateRHS

  subroutine setupIntegration(this, integrator, valuedSource, jacobianDet)
    implicit none
    class(ThermalElementDT)                          , intent(inout) :: this
    type(IntegratorPtrDT)                          , intent(in)      :: integrator
    real(rkind), dimension(integrator%ptr%integTerms), intent(out)   :: valuedSource
    real(rkind), dimension(integrator%ptr%integTerms), intent(out)   :: jacobianDet
    integer(ikind)                                                   :: i
    real(rkind), dimension(2,2)                                      :: jacobian
    valuedSource = this%getValuedSource(integrator)
    do i = 1, integrator%ptr%integTerms
       jacobian = this%jacobian(integrator%ptr%gPoint(i,1),integrator%ptr%gPoint(i,2))
       jacobianDet(i) = this%jacobianDet(jacobian)
    end do
  end subroutine setupIntegration

  function getValuedSource(this, integrator)
    implicit none
    class(ThermalElementDT), intent(inout) :: this
    type(IntegratorPtrDT) , intent(in) :: integrator
    real(rkind), dimension(integrator%ptr%integTerms) :: getValuedSource
    integer(ikind) :: i, j, nNode
    real(rkind) :: x, y
    type(NodePtrDT), dimension(:), allocatable :: node
    nNode = this%node
    do i = 1, integrator%ptr%integTerms
       node = this%node
       x = 0
       y = 0
       do j = 1, nNode
          x = x + integrator%ptr%shapeFunc(i,j)*node(j)%ptr%getx()
          y = y + integrator%ptr%shapeFunc(i,j)*node(j)%ptr%gety()
       end do
       getValuedSource(i) = this%source%func(1)%evaluate((/x,y/))
    end do
  end function getValuedSource
          

end module ThermalElement2DM
