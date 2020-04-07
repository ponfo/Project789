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
    end if
  end subroutine valueGPoints

  subroutine getG1D(this)
    implicit none
    class(IntegratorDT), intent(inout) :: this
    real(rkind), parameter :: pi = dacos(-1.d0)
    real(rkind) :: f, df, dx, r
    integer(ikind) :: i, j, iter, k, counter
    real(rkind), dimension(:), allocatable :: p0, p1, tmp
    real(rkind), dimension(:,:), allocatable :: G
    this%integTerms = this%gaussOrder
    allocate(G(2,this%gaussOrder))
    p0 = [1.d0]
    p1 = [1.d0, 0.d0]
    do k = 2, this%gaussOrder
       tmp = ((2*k-1)*[p1,0.d0]-(k-1)*[0.d0, 0.d0,p0])/k
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
       this%weight(1)   = 1.d0/2.d0
       this%gPoint(1,1) = 1.d0/3.d0
       this%gPoint(1,2) = 1.d0/3.d0
    else if(this%gaussOrder == 2) then
       this%integTerms  = 3
       allocate(this%weight(this%integTerms))
       allocate(this%gPoint(this%integTerms,2))
       this%weight(1)   = 1.d0/6.d0
       this%weight(2)   = 1.d0/6.d0
       this%weight(3)   = 1.d0/6.d0
       this%gPoint(1,1) = 1.d0/2.d0
       this%gPoint(1,2) = 0.d0
       this%gPoint(2,1) = 1.d0/2.d0
       this%gPoint(2,2) = 1.d0/2.d0
       this%gPoint(3,1) = 0.d0
       this%gPoint(3,2) = 1.d0/2.d0
    else if(this%gaussOrder == 3) then
       this%integTerms  = 4
       allocate(this%weight(this%integTerms))
       allocate(this%gPoint(this%integTerms,2))
       this%weight(1)   = -27.d0/96.d0
       this%weight(2)   = 25.d0/96.d0
       this%weight(3)   = 25.d0/96.d0
       this%weight(4)   = 25.d0/96.d0
       this%gPoint(1,1) = 1.d0/3.d0
       this%gPoint(1,2) = 1.d0/3.d0
       this%gPoint(2,1) = 1.d0/5.d0
       this%gPoint(2,2) = 1.d0/5.d0
       this%gPoint(3,1) = 3.d0/5.d0
       this%gPoint(3,2) = 1.d0/5.d0
       this%gPoint(4,1) = 1.d0/5.d0
       this%gPoint(4,2) = 3.d0/5.d0
    else
       print'(A)', '** Input Gauss Order not supported for triangular elements! **'
    end if
    call debugLog('        Allocated weights: ', size(this%weight))
    call debugLog('        Allocated gPoints: ', size(this%gPoint))
  end subroutine getGTriangle

  subroutine getGSquare(this)
    implicit none
    class(IntegratorDT), intent(inout) :: this
    real(rkind), parameter :: pi = dacos(-1.d0)
    real(rkind) :: f, df, dx, r
    integer(ikind) :: i, j, iter, k, counter
    real(rkind), dimension(:), allocatable :: p0, p1, tmp
    real(rkind), dimension(:,:), allocatable :: G
    this%integTerms = this%gaussOrder**2
    allocate(G(2,this%gaussOrder))
    p0 = [1.d0]
    p1 = [1.d0, 0.d0]
    do k = 2, this%gaussOrder
       tmp = ((2*k-1)*[p1,0.d0]-(k-1)*[0.d0, 0.d0,p0])/k
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
    if(this%gaussOrder == 2) then
       this%integTerms  = 8
       allocate(this%weight(this%integTerms))
       allocate(this%gPoint(this%integTerms,3))
       this%weight(1)   = 0.001179673492382
       this%weight(2)   = 0.001179673492382
       this%weight(3)   = 0.004402601409914
       this%weight(4)   = 0.004402601409914
       this%weight(5)   = 0.016430731923420
       this%weight(6)   = 0.016430731923420
       this%weight(7)   = 0.061320326343747
       this%weight(8)   = 0.061320326343747    
       this%gPoint(1,1) = 0.009437387888358
       this%gPoint(1,2) = 0.035220811090087
       this%gPoint(1,3) = 0.166666666666667
       this%gPoint(2,1) = 0.035220811090087
       this%gPoint(2,2) = 0.009437387888358
       this%gPoint(2,3) = 0.166666666666667
       this%gPoint(3,1) = 0.035220811090087
       this%gPoint(3,2) = 0.131445856471988
       this%gPoint(3,3) = 0.044658198978444
       this%gPoint(4,1) = 0.131445856471988
       this%gPoint(4,2) = 0.035220811090087
       this%gPoint(4,3) = 0.044658198978444
       this%gPoint(5,1) = 0.035220810850163
       this%gPoint(5,2) = 0.131445855576580
       this%gPoint(5,3) = 0.622008467032738
       this%gPoint(6,1) = 0.131445855576580
       this%gPoint(6,2) = 0.035220810850163
       this%gPoint(6,3) = 0.622008467032738
       this%gPoint(7,1) = 0.131445855576580
       this%gPoint(7,2) = 0.490562611456158
       this%gPoint(7,3) = 0.166666666666667
       this%gPoint(8,1) = 0.490562611456158
       this%gPoint(8,2) = 0.131445855576580
       this%gPoint(8,3) = 0.166666666666667
    else if(this%gaussOrder == 3) then
       this%integTerms  = 27
       allocate(this%weight(this%integTerms))
       allocate(this%gPoint(this%integTerms,3))
       this%weight(1)   = 3.068198819728420e-5
       this%weight(2)   = 3.068198819728420e-5
       this%weight(3)   = 4.909118111565470e-5
       this%weight(4)   = 2.177926162424280e-4
       this%weight(5)   = 2.177926162424280e-4 
       this%weight(6)   = 3.484681859878840e-4 
       this%weight(7)   = 2.415587821057510e-4
       this%weight(8)   = 2.415587821057510e-4
       this%weight(9)   = 3.864940513692010e-4
       this%weight(10)  = 9.662351284230000e-4
       this%weight(11)  = 9.662351284230000e-4
       this%weight(12)  = 0.001545976205477
       this%weight(13)  = 0.006858710562414
       this%weight(14)  = 0.006858710562414
       this%weight(15)  = 0.010973936899863
       this%weight(16)  = 0.007607153074595
       this%weight(17)  = 0.007607153074595
       this%weight(18)  = 0.012171444919352
       this%weight(19)  = 0.001901788268649
       this%weight(20)  = 0.001901788268649
       this%weight(21)  = 0.003042861229838
       this%weight(22)  = 0.013499628508586
       this%weight(23)  = 0.013499628508586
       this%weight(24)  = 0.021599405613738
       this%weight(25)  = 0.014972747367084
       this%weight(26)  = 0.014972747367084
       this%weight(27)  = 0.023956395787334
       this%gPoint(1,1)  = 0.001431498841332
       this%gPoint(1,2)  = 0.011270166537926
       this%gPoint(1,3)  = 0.100000000000000
       this%gPoint(2,1)  = 0.011270166537926
       this%gPoint(2,2)  = 0.001431498841332
       this%gPoint(2,3)  = 0.100000000000000
       this%gPoint(3,1)  = 0.006350832689629
       this%gPoint(3,2)  = 0.006350832689629
       this%gPoint(3,3)  = 0.100000000000000
       this%gPoint(4,1)  = 0.006350832689629
       this%gPoint(4,2)  = 0.050000000000000
       this%gPoint(4,3)  = 0.056350832689629
       this%gPoint(5,1)  = 0.050000000000000
       this%gPoint(5,2)  = 0.006350832689629
       this%gPoint(5,3)  = 0.056350832689629
       this%gPoint(6,1)  = 0.028175416344815
       this%gPoint(6,2)  = 0.028175416344815
       this%gPoint(6,3)  = 0.056350832689629
       this%gPoint(7,1)  = 0.011270166537926
       this%gPoint(7,2)  = 0.088729833462074
       this%gPoint(7,3)  = 0.012701665379258
       this%gPoint(8,1)  = 0.088729833462074
       this%gPoint(8,2)  = 0.011270166537926
       this%gPoint(8,3)  = 0.012701665379258
       this%gPoint(9,1)  = 0.050000000000000
       this%gPoint(9,2)  = 0.050000000000000
       this%gPoint(9,3)  = 0.012701665379258
       this%gPoint(10,1) = 0.006350832689629
       this%gPoint(10,2) = 0.050000000000000
       this%gPoint(10,3) = 0.443649167310371
       this%gPoint(11,1) = 0.050000000000000
       this%gPoint(11,2) = 0.006350832689629
       this%gPoint(11,3) = 0.443649167310371
       this%gPoint(12,1) = 0.028175416344815
       this%gPoint(12,2) = 0.028175416344815
       this%gPoint(12,3) = 0.443649167310371
       this%gPoint(13,1) = 0.028175416344815
       this%gPoint(13,2) = 0.221824583655185
       this%gPoint(13,3) = 0.250000000000000
       this%gPoint(14,1) = 0.221824583655185
       this%gPoint(14,2) = 0.028175416344815
       this%gPoint(14,3) = 0.250000000000000
       this%gPoint(15,1) = 0.125000000000000
       this%gPoint(15,2) = 0.125000000000000
       this%gPoint(15,3) = 0.250000000000000
       this%gPoint(16,1) = 0.050000000000000
       this%gPoint(16,2) = 0.393649167310371
       this%gPoint(16,3) = 0.056350832689629
       this%gPoint(17,1) = 0.393649167310371
       this%gPoint(17,2) = 0.050000000000000
       this%gPoint(17,3) = 0.056350832689629
       this%gPoint(18,1) = 0.221824583655185
       this%gPoint(18,2) = 0.221824583655185
       this%gPoint(18,3) = 0.056350832689629
       this%gPoint(19,1) = 0.011270166537926
       this%gPoint(19,2) = 0.088729833462074
       this%gPoint(19,3) = 0.787298334620741
       this%gPoint(20,1) = 0.088729833462074
       this%gPoint(20,2) = 0.011270166537926
       this%gPoint(20,3) = 0.787298334620741
       this%gPoint(21,1) = 0.050000000000000
       this%gPoint(21,2) = 0.050000000000000
       this%gPoint(21,3) = 0.787298334620741
       this%gPoint(22,1) = 0.050000000000000
       this%gPoint(22,2) = 0.393649167310371
       this%gPoint(22,3) = 0.443649167310371
       this%gPoint(23,1) = 0.393649167310371
       this%gPoint(23,2) = 0.050000000000000
       this%gPoint(23,3) = 0.443649167310371
       this%gPoint(24,1) = 0.221824583655185
       this%gPoint(24,2) = 0.221824583655185
       this%gPoint(24,3) = 0.443649167310371
       this%gPoint(25,1) = 0.088729833462074
       this%gPoint(25,2) = 0.698568501158667
       this%gPoint(25,3) = 0.100000000000000
       this%gPoint(26,1) = 0.698568501158667
       this%gPoint(26,2) = 0.088729833462074
       this%gPoint(26,3) = 0.100000000000000
       this%gPoint(27,1) = 0.393649167310371
       this%gPoint(27,2) = 0.393649167310371
       this%gPoint(27,3) = 0.100000000000000
    else
       print'(A)', '** Input Gauss Order not supported for tetrahedral elements! **'
    end if
    call debugLog('        Allocated weights: ', size(this%weight))
    call debugLog('        Allocated gPoints: ', size(this%gPoint))
  end subroutine getGTetrahedron

end module IntegratorM
