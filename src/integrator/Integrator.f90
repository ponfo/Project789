module IntegratorM
  use UtilitiesM
  use DebuggerM
  implicit none
  private
  public :: IntegratorDT, integrator
  type :: IntegratorDT
     integer(ikind)                               :: gaussOrder
     integer(ikind)                               :: integTerms
     real(rkind), dimension(:)      , allocatable :: weight
     real(rkind), dimension(:,:)    , allocatable :: gPoint
     real(rkind), dimension(:,:)    , allocatable :: shapeFunc
     real(rkind), dimension(:,:,:)  , allocatable :: dShapeFunc
     real(rkind), dimension(:,:,:,:), allocatable :: ddShapeFunc
   contains 
     procedure, public  :: init
     procedure, public  :: valueGPoints
     procedure, private :: getG1D
     procedure, private :: getGTriangle
     procedure, private :: getGSquare
     procedure, private :: getGTetrahedron
     procedure, private :: getGHexahedron

     procedure, public :: getGaussOrder
     procedure, public :: getIntegTerms
     procedure, public :: getWeight
     procedure, public :: getGPoint
     procedure, public :: getShapeFunc
     procedure, public :: getDShapeFunc
     procedure, public :: getDDShapeFunc
     procedure, public :: getWeightFull
     procedure, public :: getGPointFull
     procedure, public :: getShapeFuncFull
     procedure, public :: getDShapeFuncFull
     procedure, public :: getDDShapeFuncFull
  end type IntegratorDT

  interface integrator
     procedure constructor
  end interface integrator
  
contains

  type(IntegratorDT) function constructor(gaussOrder, type)
    implicit none
    integer(ikind), intent(in) :: gaussOrder
    character(*), intent(in) :: type
    call constructor%init(gaussOrder, type)
  end function constructor

  subroutine init(this, gaussOrder, type)
    implicit none
    class(IntegratorDT), intent(inout) :: this
    integer(ikind), intent(in) :: gaussOrder
    character(*), intent(in) :: type
    this%gaussOrder = gaussOrder
    call this%valueGPoints(type)
  end subroutine init
  
  subroutine valueGPoints(this, type)
    implicit none
    class(IntegratorDT), intent(inout) :: this
    character(*), intent(in) :: type
    if     (trim(type) == 'Line'          .or. &
            trim(type) == 'line'          .or. &
            trim(type) == 'LINE')           then
       call this%getG1D()
    else if(trim(type) == 'Triangle'      .or. &
            trim(type) == 'triangle'      .or. &
            trim(type) == 'triang'        .or. & 
            trim(type) == 'TRIANGLE')       then
       call this%getGTriangle()
    else if(trim(type) == 'Quadrilateral' .or. &
            trim(type) == 'quadrilateral' .or. &
            trim(type) == 'QUADRILATERAL' .or. &
            trim(type) == 'quad')           then
       call this%getGSquare()
    else if(trim(type) == 'Tetrahedron'   .or. &
            trim(type) == 'tetrahedron'   .or. &
            trim(type) == 'TETRAHEDRON'   .or. &
            trim(type) == 'tetra')          then
       call this%getGTetrahedron()
    else if(trim(type) == 'Hexahedron'    .or. &
            trim(type) == 'hexahedron'    .or. &
            trim(type) == 'HEXAHEDRON'    .or. &
            trim(type) == 'hexa')           then
       call this%getGHexahedron()
    end if
  end subroutine valueGPoints

  subroutine getG1D(this)
    implicit none
    class(IntegratorDT), intent(inout) :: this
    real(rkind), parameter :: pi = dacos(-1._rkind)
    real(rkind) :: f, df, dx, r
    integer(ikind) :: i, j, iter, k, counter
    real(rkind), dimension(:), allocatable :: p0, p1, tmp
    real(rkind), dimension(:,:), allocatable :: G
    this%integTerms = this%gaussOrder
    allocate(G(2,this%gaussOrder))
    p0 = [1._rkind]
    p1 = [1._rkind, 0._rkind]
    do k = 2, this%gaussOrder
       tmp = ((2*k-1)*[p1,0._rkind]-(k-1)*[0._rkind, 0._rkind,p0])/k
       p0 = p1; p1 = tmp
    end do
    do i = 1, this%gaussOrder
       r = cos(pi*(i-0.25)/(this%gaussOrder+0.5))
       do iter = 1, 10
          f = p1(1); df = 0.
          do k = 2, size(p1)
             df = f + r*df
             f  = p1(k) + r * f
          end do
          dx =  f / df
          r = r - dx
          if (abs(dx)<10*epsilon(dx)) exit
       end do
       G(1,i) = r
       G(2,i) = 2/((1-r**2)*df**2)
    end do
    allocate(this%weight(this%integTerms))
    allocate(this%gPoint(this%integTerms,1))
    this%weight(:) = G(2,:)
    this%gPoint(:,1) = G(1,:)
    call debugLog('        Allocated weights: ', size(this%weight))
    call debugLog('        Allocated gPoints: ', size(this%gPoint))
    deallocate(G)
  end subroutine getG1D

  subroutine getGTriangle(this)
    implicit none
    class(IntegratorDT), intent(inout) :: this
    if(this%gaussOrder == 1) then
       this%integTerms  = 1
       allocate(this%weight(this%integTerms))
       allocate(this%gPoint(this%integTerms,2))
       this%weight(1)   = 1._rkind/2._rkind
       this%gPoint(1,1) = 1._rkind/3._rkind
       this%gPoint(1,2) = 1._rkind/3._rkind
    else if(this%gaussOrder == 2) then
       this%integTerms  = 3
       allocate(this%weight(this%integTerms))
       allocate(this%gPoint(this%integTerms,2))
       this%weight(1)   = 1._rkind/6._rkind
       this%weight(2)   = 1._rkind/6._rkind
       this%weight(3)   = 1._rkind/6._rkind
       this%gPoint(1,1) = 1._rkind/2._rkind
       this%gPoint(1,2) = 0._rkind
       this%gPoint(2,1) = 1._rkind/2._rkind
       this%gPoint(2,2) = 1._rkind/2._rkind
       this%gPoint(3,1) = 0._rkind
       this%gPoint(3,2) = 1._rkind/2._rkind
    else if(this%gaussOrder == 3) then
       this%integTerms  = 4
       allocate(this%weight(this%integTerms))
       allocate(this%gPoint(this%integTerms,2))
       this%weight(1)   = -27._rkind/96._rkind
       this%weight(2)   = 25._rkind/96._rkind
       this%weight(3)   = 25._rkind/96._rkind
       this%weight(4)   = 25._rkind/96._rkind
       this%gPoint(1,1) = 1._rkind/3._rkind
       this%gPoint(1,2) = 1._rkind/3._rkind
       this%gPoint(2,1) = 1._rkind/5._rkind
       this%gPoint(2,2) = 1._rkind/5._rkind
       this%gPoint(3,1) = 3._rkind/5._rkind
       this%gPoint(3,2) = 1._rkind/5._rkind
       this%gPoint(4,1) = 1._rkind/5._rkind
       this%gPoint(4,2) = 3._rkind/5._rkind
    else
       print'(A)', '** Input Gauss Order not supported for triangular elements! **'
    end if
    call debugLog('        Allocated weights: ', size(this%weight))
    call debugLog('        Allocated gPoints: ', size(this%gPoint))
  end subroutine getGTriangle

  subroutine getGSquare(this)
    implicit none
    class(IntegratorDT), intent(inout) :: this
    real(rkind), parameter :: pi = dacos(-1._rkind)
    real(rkind) :: f, df, dx, r
    integer(ikind) :: i, j, iter, k, counter
    real(rkind), dimension(:), allocatable :: p0, p1, tmp
    real(rkind), dimension(:,:), allocatable :: G
    this%integTerms = this%gaussOrder**2
    allocate(G(2,this%gaussOrder))
    p0 = [1._rkind]
    p1 = [1._rkind, 0._rkind]
    do k = 2, this%gaussOrder
       tmp = ((2*k-1)*[p1,0._rkind]-(k-1)*[0._rkind, 0._rkind,p0])/k
       p0 = p1; p1 = tmp
    end do
    do i = 1, this%gaussOrder
       r = cos(pi*(i-0.25)/(this%gaussOrder+0.5))
       do iter = 1, 10
          f = p1(1); df = 0.
          do k = 2, size(p1)
             df = f + r*df
             f  = p1(k) + r * f
          end do
          dx =  f / df
          r = r - dx
          if (abs(dx)<10*epsilon(dx)) exit
       end do
       G(1,i) = r
       G(2,i) = 2/((1-r**2)*df**2)
    end do
    allocate(this%weight(this%integTerms))
    allocate(this%gPoint(this%integTerms,2))
    counter = 0
    do i = 1, this%gaussOrder
       do j = 1, this%gaussOrder
          counter = counter + 1
          this%weight(counter) = G(2,i)*G(2,j)
          this%gPoint(counter,1) = G(1,i)
          this%gPoint(counter,2) = G(1,j)
       end do
    end do
    call debugLog('        Allocated weights: ', size(this%weight))
    call debugLog('        Allocated gPoints: ', size(this%gPoint))
    deallocate(G)
  end subroutine getGSquare

  subroutine getGTetrahedron(this)
    implicit none
    class(IntegratorDT), intent(inout) :: this
    if(this%gaussOrder == 1) then
       this%integTerms  = 1
       allocate(this%weight(this%integTerms))
       allocate(this%gPoint(this%integTerms,3))
       this%weight(1)   = 1
       this%gPoint(1,1) = .25
       this%gPoint(1,2) = .25
       this%gPoint(1,3) = .25
    else if(this%gaussOrder == 2) then
       this%integTerms  = 5
       allocate(this%weight(this%integTerms))
       allocate(this%gPoint(this%integTerms,3))
       this%weight(1)   = -.8
       this%weight(2)   = .45
       this%weight(3)   = .45
       this%weight(4)   = .45
       this%weight(5)   = .45
       this%gPoint(1,1) = .25
       this%gPoint(1,2) = .25
       this%gPoint(1,3) = .25
       this%gPoint(2,1) = .166666666666667
       this%gPoint(2,2) = .166666666666667
       this%gPoint(2,3) = .166666666666667
       this%gPoint(3,1) = .5
       this%gPoint(3,2) = .166666666666667
       this%gPoint(3,3) = .166666666666667
       this%gPoint(4,1) = .166666666666667
       this%gPoint(4,2) = .5
       this%gPoint(4,3) = .166666666666667
       this%gPoint(5,1) = .166666666666667
       this%gPoint(5,2) = .166666666666667
       this%gPoint(5,3) = .5
    else if(this%gaussOrder == 3) then
       this%integTerms  = 11
       allocate(this%weight(this%integTerms))
       allocate(this%gPoint(this%integTerms,3))
       this%weight(1)   = 0.013155555555555
       this%weight(2)   = 0.007622222222222
       this%weight(3)   = 0.007622222222222
       this%weight(4)   = 0.007622222222222
       this%weight(5)   = 0.007622222222222
       this%weight(6)   = 0.024888888888888
       this%weight(7)   = 0.024888888888888
       this%weight(8)   = 0.024888888888888
       this%weight(9)   = 0.024888888888888
       this%weight(10)  = 0.024888888888888
       this%weight(11)  = 0.024888888888888
       this%gPoint(1,1) = .25
       this%gPoint(1,2) = .25
       this%gPoint(1,3) = .25
       this%gPoint(2,1) = .0714285714285714
       this%gPoint(2,2) = .0714285714285714
       this%gPoint(2,3) = .785714285714286
       this%gPoint(3,1) = .0714285714285714
       this%gPoint(3,2) = .0714285714285714
       this%gPoint(3,3) = .0714285714285714
       this%gPoint(4,1) = .785714285714286
       this%gPoint(4,2) = .0714285714285714
       this%gPoint(4,3) = .0714285714285714
       this%gPoint(5,1) = .0714285714285714
       this%gPoint(5,2) = .785714285714286
       this%gPoint(5,3) = .0714285714285714
       this%gPoint(6,1) = 0.399403576166799
       this%gPoint(6,2) = 0.100596423833201
       this%gPoint(6,3) = 0.100596423833201
       this%gPoint(7,1) = 0.399403576166799
       this%gPoint(7,2) = 0.399403576166799
       this%gPoint(7,3) = 0.100596423833201
       this%gPoint(8,1) = 0.100596423833201
       this%gPoint(8,2) = 0.399403576166799
       this%gPoint(8,3) = 0.399403576166799
       this%gPoint(9,1) = 0.100596423833201
       this%gPoint(9,2) = 0.100596423833201
       this%gPoint(9,3) = 0.399403576166799
       this%gPoint(10,1) = 0.100596423833201
       this%gPoint(10,2) = 0.399403576166799
       this%gPoint(10,3) = 0.100596423833201
       this%gPoint(11,1) = 0.399403576166799
       this%gPoint(11,2) = 0.100596423833201
       this%gPoint(11,3) = 0.399403576166799
    else
       print'(A)', '** Input Gauss Order not supported for tetrahedral elements! **'
    end if
    call debugLog('        Allocated weights: ', size(this%weight))
    call debugLog('        Allocated gPoints: ', size(this%gPoint))
  end subroutine getGTetrahedron

  subroutine getGHexahedron(this)
    implicit none
    class(IntegratorDT), intent(inout) :: this
    real(rkind), parameter :: pi = dacos(-1._rkind)
    real(rkind) :: f, df, dx, r
    integer(ikind) :: i, j, iter, k, counter
    real(rkind), dimension(:), allocatable :: p0, p1, tmp
    real(rkind), dimension(:,:), allocatable :: G
    this%integTerms = this%gaussOrder**3
    allocate(G(2,this%gaussOrder))
    p0 = [1._rkind]
    p1 = [1._rkind, 0._rkind]
    do k = 2, this%gaussOrder
       tmp = ((2*k-1)*[p1,0._rkind]-(k-1)*[0._rkind, 0._rkind,p0])/k
       p0 = p1; p1 = tmp
    end do
    do i = 1, this%gaussOrder
       r = cos(pi*(i-0.25)/(this%gaussOrder+0.5))
       do iter = 1, 10
          f = p1(1); df = 0.
          do k = 2, size(p1)
             df = f + r*df
             f  = p1(k) + r * f
          end do
          dx =  f / df
          r = r - dx
          if (abs(dx)<10*epsilon(dx)) exit
       end do
       G(1,i) = r
       G(2,i) = 2/((1-r**2)*df**2)
    end do
    allocate(this%weight(this%integTerms))
    allocate(this%gPoint(this%integTerms,3))
    counter = 0
    do i = 1, this%gaussOrder
       do j = 1, this%gaussOrder
          do k = 1, this%gaussOrder
             counter = counter + 1
             this%weight(counter) = G(2,i)*G(2,j)*G(2,k)
             this%gPoint(counter,1) = G(1,i)
             this%gPoint(counter,2) = G(1,j)
             this%gPoint(counter,3) = G(1,k)
          end do
       end do
    end do
    call debugLog('        Allocated weights: ', size(this%weight))
    call debugLog('        Allocated gPoints: ', size(this%gPoint))
    deallocate(G)
  end subroutine getGHexahedron

  integer(ikind) pure function getGaussOrder(this)
    implicit none
    class(IntegratorDT), intent(in) :: this
    getGaussOrder = this%gaussOrder
  end function getGaussOrder

  integer(ikind) pure function getIntegTerms(this)
    implicit none
    class(IntegratorDT), intent(in) :: this
    getIntegTerms = this%integTerms
  end function getIntegTerms

  real(rkind) pure function getWeight(this, i)
    implicit none
    class(IntegratorDT), intent(in) :: this
    integer(ikind)     , intent(in) :: i
    getWeight = this%weight(i)
  end function getWeight

  real(rkind) pure function getGPoint(this, i, j)
    implicit none
    class(IntegratorDT), intent(in) :: this
    integer(ikind)     , intent(in) :: i
    integer(ikind)     , intent(in) :: j
    getGPoint = this%gPoint(i,j)
  end function getGPoint

  real(rkind) pure function getShapeFunc(this, i, j)
    implicit none
    class(IntegratorDT), intent(in) :: this
    integer(ikind)     , intent(in) :: i
    integer(ikind)     , intent(in) :: j
    getShapeFunc = this%shapeFunc(i,j)
  end function getShapeFunc

  real(rkind) pure function getDShapeFunc(this, i, j, k)
    implicit none
    class(IntegratorDT), intent(in) :: this
    integer(ikind)     , intent(in) :: i
    integer(ikind)     , intent(in) :: j
    integer(ikind)     , intent(in) :: k
    getDShapeFunc = this%dShapeFunc(i,j,k)
  end function getDShapeFunc

  real(rkind) pure function getDDShapeFunc(this, i, j, k, l)
    implicit none
    class(IntegratorDT), intent(in) :: this
    integer(ikind)     , intent(in) :: i
    integer(ikind)     , intent(in) :: j
    integer(ikind)     , intent(in) :: k
    integer(ikind)     , intent(in) :: l
    getDDShapeFunc = this%ddShapeFunc(i,j,k,l)
  end function getDDShapeFunc

  pure function getWeightFull(this)
    implicit none
    class(IntegratorDT)              , intent(in)  :: this
    real(rkind)        , dimension(:), allocatable :: getWeightFull
    getWeightFull = this%weight
  end function getWeightFull

  pure function getGPointFull(this)
    implicit none
    class(IntegratorDT)                , intent(in)  :: this
    real(rkind)        , dimension(:,:), allocatable :: getGPointFull
    getGPointFull = this%gPoint
  end function getGPointFull

  pure function getShapeFuncFull(this)
    implicit none
    class(IntegratorDT)                , intent(in)  :: this
    real(rkind)        , dimension(:,:), allocatable :: getShapeFuncFull
    getShapeFuncFull = this%shapeFunc
  end function getShapeFuncFull

  pure function getDShapeFuncFull(this)
    implicit none
    class(IntegratorDT)                  , intent(in)  :: this
    real(rkind)        , dimension(:,:,:), allocatable :: getDShapeFuncFull
    getDShapeFuncFull = this%dShapeFunc
  end function getDShapeFuncFull

  pure function getDDShapeFuncFull(this)
    implicit none
    class(IntegratorDT)                    , intent(in)  :: this
    real(rkind)        , dimension(:,:,:,:), allocatable :: getDDShapeFuncFull
    getDDShapeFuncFull = this%ddShapeFunc
  end function getDDShapeFuncFull
  
end module IntegratorM
