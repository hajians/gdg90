!> a program to check the eigenvalues of the interior penalty
!! matrix for different penalization.

program test10

  use triangulate2d
  use alldgspace2d
  use dgopt2d
  use interpolate2d
  use mylinalgebra
  use func2mat
  use functional2mat
  use normdg2d
  use nestddmesh2d
  use intfspace2d

  use OptSchwarz

  implicit none

  type(tri2d)      :: ThC, ThF
  type(nestmesh2d) :: nest2d

  integer :: k

  type(primaldg2d) :: U, V

  real, allocatable :: A(:,:), rhsF(:)
  type(spmat) :: spA, spEnergy

  procedure(bilinearformPrimal), pointer :: bilinpr
  procedure(L2prod), pointer :: rhs

  !! element order: k
  k = 1

  !! reading the coarse and fine meshes
  call readtri2d('mesh/Hh/BOX-3.1',ThC)
  call buildtri2d(ThC)

  call readtri2d('mesh/Hh/BOX-3.2',ThF)
  call buildtri2d(ThF)

  call buildnestmesh2d(ThC, ThF, nest2d)

  call plottri2d('plot/Mesh2D',ThC)
  call plottri2d('plot/Mesh2DF',ThF)

  call plotnestinterface("plot/interface",nest2d)

  !! build primal DG spaces
  call buildprimaldg2d(ThF, k, U)
  call buildprimaldg2d(ThF, k, V)

  ! tuning the penalty scaling
  ! V%paramFlux%cp(1)%edgedof = sqrt(V%paramFlux%cp(1)%edgedof)
  ! V%paramFlux%cp(2)%edgedof = sqrt(V%paramFlux%cp(2)%edgedof)
  
  write(6,*) "start assembly"

  ! mapping bilinpr to aIP
  bilinpr => aIPH
  ! mapping
  rhs => L2prodF

  !! full matrix
  ! call bilinpr2mat(bilinpr,U,V,A)
  
  !! sparse matrix
  call bilinpr2spmat(bilinpr,U,V,spA)
  call writespmat('out/stiffness',spA)

  ! nullify bilinear pointer
  ! bilinpr => null()
  ! bilinpr => energy

  ! call bilinpr2spmat(bilinpr,U,V,spEnergy)
  ! call writespmat('energy',spEnergy)

  write(6,*) "end assembly"
  

contains

  !> bilinear form for IP
  subroutine aIPH(U,V,out)
    type(primaldg2d), intent(inout) :: U, V

    real, intent(inout) :: out

    out = 0.0

    U%flux = Dh( U%sc )
    call trace( U%sc,   U%trSc)
    call trace( U%flux, U%trflux)

    V%flux = Dh( V%sc )
    call trace( V%sc,   V%trSc)
    call trace( V%flux, V%trflux)

    out = out + dot0Th( U%flux, V%flux ) + dot0Th(U%sc, V%sc)

    call avetr( U%trFlux, U%singleFlux )
    call jumptr( V%trSc, V%singleFlux )
    out = out - dotproduct( U%singleFlux, V%singleFlux )

    call jumptr( U%trSc, U%singleFlux )
    call multiplySltrFlux( V%singleFlux, V%paramFlux, V%singleFlux )
    out = out + dotproduct( U%singleFlux, V%singleFlux )

    call avetr( V%trFlux, V%singleFlux )
    out = out - dotproduct( U%singleFlux, V%singleFlux )

    ! call jumptr( U%trFlux, U%JumpFlux )
    ! call jumptr( V%trFlux, V%JumpFlux )
    ! call multiplySltr( U%JumpFlux, U%MuMinusOne, U%JumpFlux)
    ! out = out - 0.5*dotproduct( U%JumpFlux, V%JumpFlux )

  end subroutine aIPH

    !> bilinear form for IP
  subroutine energy(U,V,out)
    type(primaldg2d), intent(inout) :: U, V

    real, intent(inout) :: out

    out = 0.0

    U%flux = Dh( U%sc )
    call trace( U%sc,   U%trSc)
    call trace( U%flux, U%trflux)

    V%flux = Dh( V%sc )
    call trace( V%sc,   V%trSc)
    call trace( V%flux, V%trflux)

    out = out + dot0Th( U%flux, V%flux ) + dot0Th(U%sc, V%sc)

    call jumptr( V%trSc, V%singleFlux )
    call jumptr( U%trSc, U%singleFlux )

    !call multiplySltrFlux( V%singleFlux, V%paramFlux, V%singleFlux )
    out = out + 0.1*dotproduct( U%singleFlux, V%singleFlux )

  end subroutine energy

  !> right hand side functional
  function L2prodF(V) result(out)
    type(dg2d), intent(in) :: V
    real :: out

    type(dg2d) :: Ihf

    procedure(fct), pointer :: func

    if ( V%init .eqv. .false. ) then
       write(6,*) "error in f(V):", V%init 
       stop
    endif

    Ihf = builddg2d(V%Th, V%lo1d, V%lo2d, V%k)

    func => f
    call Ih(func, Ihf)

    out = dot0Th( Ihf, V )
  end function L2prodF

  !> right-hand side function
  function f(coor) result(out)
    real, dimension(:,:), intent(in) :: coor
    real, dimension( size(coor,1) ) :: out

    real, dimension( size(coor,1) ) :: x, y
    real :: pi = acos( -1.0 )
    integer :: KK

    KK = 2
    x = coor(:,1)
    y = coor(:,2)
    out = (x**2-x)*(y**2-y)     

    out = 8*pi**2*sin( KK*pi*x ) * sin(KK*pi*y)
    out = -4*pi * ( cos(2*pi*x) - 2*pi*x*sin(2*pi*x) ) * sin(2*pi*y)
  end function f

  !! -------------------------------
  !> Displaying a Matrix in the shell
  !! @param [in] M matrix of size i,j
  subroutine ShowMat(M)
    real, intent(in), dimension(:,:) :: M
    integer :: dimx, dimy, i
    character(len=80) :: lb

    dimx = size(M, dim=1)
    dimy = size(M, dim=2)

    do i = 1, dimx
       lb = "('[',I2,']', 1000F8.2)"
       write(*,lb) i, M(i,:)
    enddo
    write(*,"(/)")

  end subroutine ShowMat

  subroutine ShowVec(V)
    real, intent(in), dimension(:) :: V
    integer :: dimx, i
    character(len=80) :: lb
    dimx = size(V, dim=1)

    lb = "('[',I4,']', 1000F8.4)"
    do i=1, dimx
       write(*,lb) i, V(i)
    enddo
    write(*,"(/)")
  end subroutine ShowVec
  !! -------------------------------
end program
