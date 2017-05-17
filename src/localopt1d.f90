!> @file locopt1d.f90

module localopt1d
  use jacobip
  use mylinalgebra
  implicit none

  type locopt1d

     logical :: init = .false.
     
     integer :: k
     integer :: Np

     real, allocatable :: refx(:)
     
     real, allocatable :: va1d(:,:)
     real, allocatable :: vr1d(:,:)
     real, allocatable :: dr1d(:,:)
     real, allocatable ::  m1d(:,:) ! mass matrix
     real, allocatable ::  s1d(:,:) ! gradient matrix
     real, allocatable :: gg1d(:,:) ! stiffness matrix

     real, allocatable :: va1dinv(:,:)
     real, allocatable ::  m1dinv(:,:)

  end type locopt1d

contains

  !> takes degree of polynomials and gives the coordinates
  !! of the Lagrange basis functions in the reference element.
  !> @param [in] k degree of polynomials
  !> @param [out] xy the coordinates of the nodal basis functions
  function nodes1d(k) result(x)
    integer, intent(in) :: k
    real, dimension(k+1) :: x

    integer :: i

    x = 0.0

    do i=1, k+1
       x(i) = -1 + 2.0/k * (i-1)
    end do

  end function nodes1d

  !> builds the local operators in 1d acting on the reference element
  !> @param [out] lo1d is an object containing all of these operators
  subroutine buildlocopt1d(k,lo1d)
    integer,        intent(in)    :: k
    type(locopt1d), intent(inout) :: lo1d

    if ( lo1d%init.eqv..true. ) then
       write(6,*) "error in buildlocopt1d: lo1d is already initialized"
       stop
    endif
    
    lo1d%k  = k
    lo1d%Np = k+1

    allocate( lo1d%refx( lo1d%Np ) )

    lo1d%refx = nodes1d(k)

    allocate( lo1d%va1d( lo1d%Np, lo1d%Np ) )
    allocate( lo1d%vr1d( lo1d%Np, lo1d%Np ) )
    allocate( lo1d%dr1d( lo1d%Np, lo1d%Np ) )
    allocate(  lo1d%m1d( lo1d%Np, lo1d%Np ) )
    allocate(  lo1d%s1d( lo1d%Np, lo1d%Np ) )
    allocate( lo1d%gg1d( lo1d%Np, lo1d%Np ) )

    allocate( lo1d%va1dinv( lo1d%Np, lo1d%Np ) )
    allocate(  lo1d%m1dinv( lo1d%Np, lo1d%Np ) )

    lo1d%va1d = VanderMonde( lo1d%refx )
    lo1d%vr1d = GradVaMonde( lo1d%refx )

    ! not completed!!!
    lo1d%va1dinv = inv( lo1d%va1d )

    lo1d%dr1d   = matmul( lo1d%vr1d, lo1d%va1dinv )
    lo1d%m1dinv = matmul( lo1d%va1d, transpose(lo1d%va1d) )
    lo1d%m1d    = inv( lo1d%m1dinv )
    lo1d%s1d    = matmul( lo1d%m1d, lo1d%dr1d )
    lo1d%gg1d   = matmul( transpose(lo1d%dr1d), matmul(lo1d%m1d, lo1d%dr1d) )

    !     I1D  = VX(-1.,1.,Deg)
    ! Va1D = VanderMonde(I1D)
    ! Vr1D = GradVaMonde(I1D)
    ! Call MIGS(Va1D, N_P, Va1DInv, IDS) ! MATRIX INVERSION; COPY TO Va1DInv

    ! Dr1D = MATMUL(Vr1D,Va1DInv)
    ! M1DInv  = MATMUL(Va1D, transpose(Va1D))
    ! Call MIGS(M1DInv, N_P, M1D, IDS)! MATRIX INVERSION; COPY TO M1D

    ! S1D  = MATMUL(M1D,DR1D)
    ! GradGrad1D = MATMUL( transpose(Dr1D), MATMUL(M1D, Dr1D) )
    lo1d%init = .true.
  end subroutine buildlocopt1d
  

end module localopt1d
