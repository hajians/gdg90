!> @file locopt2d.f90

module localopt2d
  use simplex2dp
  use mylinalgebra
  
  implicit none

  type locopt2d

     logical :: init = .false.
     
     integer :: k
     integer :: Np
     integer :: Npfaces
     
     real,    allocatable :: refxy(:,:)
     integer, allocatable :: nmask2d(:,:)

     real, allocatable :: va2d(:,:)

     real, allocatable :: grva2dr(:,:)
     real, allocatable :: grva2ds(:,:)
     
     real, allocatable :: dr2d(:,:)
     real, allocatable :: ds2d(:,:)

     real, allocatable :: m2d(:,:)
     real, allocatable :: sr2d(:,:)
     real, allocatable :: ss2d(:,:)

     real, allocatable :: stfrr(:,:)
     real, allocatable :: stfss(:,:)
     real, allocatable :: stfrs(:,:)
     real, allocatable :: stfsr(:,:)
  end type locopt2d

  private dmat
  private massmat
  private smat
  private stiffmat
  
contains

  !> takes degree of polynomials and gives the coordinates
  !! of the Lagrange basis functions in the reference element.
  !> @param [in] k degree of polynomials
  !> @param [out] xy the coordinates of the nodal basis functions
  function nodes2d(k) result(xy)
    integer, intent(in) :: k
    real, dimension( (k+1)*(k+2)/2, 2) :: xy

    real, dimension(3,2) :: refvtx

    real    :: l1, l2, l3
    integer :: i, j, m


    xy = 0.0

    refvtx(1,:) = (/ -1., -1. /)  ! vertices of the reference element
    refvtx(2,:) = (/ +1., -1. /)
    refvtx(3,:) = (/ -1., +1. /)

    do i=0, k
       do j=0, (k-i)
          l1 = real(i)/real(k)  ! barrycentric coordinates
          l3 = real(j)/real(k)
          l2 = 1. - l1 - l3

          m = j + (k+1)*i + 1 - i*(i-1)/2 ! map from 2d to 1d
          xy(m,:) = l1*refvtx(3,:) + l2*refvtx(1,:) + l3*refvtx(2,:)
       end do
    end do

  end function nodes2d

  !> builds a map between the face and face's nodes
  !> @param [in] k degree of polynomial
  function nodesmask2d(k)
    integer, intent(in) :: k
    integer, allocatable, dimension(:,:) :: nodesmask2d
    
    integer :: faces, Npfaces
    integer :: i, j
    
    faces   = 3
    Npfaces = k+1

    allocate( nodesmask2d( Npfaces, faces ) )

    nodesmask2d = 0

    nodesmask2d(:,1) = (/ (i, i=1,Npfaces) /)

    do i=0, k
       j = k-i
       nodesmask2d(i+1,2) = j + (k+1)*i + 1 - i*(i-1)/2
       j = 0
       nodesmask2d(i+1,3) = j + (k+1)*i + 1 - i*(i-1)/2
    end do
    
  end function nodesmask2d

  !> Dmatrices2D gives the Dr and Ds matrices for ref. element.
  !> @param [in] k
  !> @param [in] coorxy
  !> @param [in] va2d
  subroutine dmat(k, Np, coorxy, va2d,dr2d,ds2d)
    real, dimension(Np,2), intent(in)   :: coorxy
    integer,               intent(in)   :: k
    integer,               intent(in)   :: Np
    real, dimension(Np,Np), intent(in)  :: va2d

    real, dimension(Np,Np), intent(inout) :: dr2d, ds2d

    real, dimension(Np,Np) :: grva2dr, grva2ds, va2dinv

    call gradvandermonde2d(k, coorxy, grva2dr, grva2ds)

    va2dinv = inv( va2d )

    dr2d = matmul(grva2dr, va2dinv)
    ds2d = matmul(grva2ds, va2dinv)
    
  end subroutine dmat

  !> mass matrix in 2d on the reference element
  !> @param [in] va2d
  function massmat(Np,va2d) result(m2d)
    integer, intent(in)                :: Np
    real, dimension(Np,Np), intent(in) :: va2d
    
    real, dimension(Np,Np) :: m2d

    m2d = inv( matmul( va2d, transpose(va2d) ) )
  end function massmat

  !> smat produces gradient matrices in 2D. \int l_i d_{r} l_j
  subroutine smat(Np,m2d,dr2d,ds2d,sr2d,ss2d)
    integer, intent(in)                :: Np
    real, dimension(Np,Np), intent(in) :: m2d, dr2d,ds2d

    real, dimension(Np,Np), intent(inout) :: sr2d, ss2d

    sr2d = matmul( m2d, dr2d)
    ss2d = matmul( m2d, ds2d)
  end subroutine smat

  !> stiffmat computes stifness \int d_{r} l_i d_{r} l_j 
  subroutine stiffmat(Np,dr2d,ds2d,sr2d,ss2d,stfrr,stfss,stfrs,stfsr)
    integer, intent(in)                :: Np
    real, dimension(Np,Np), intent(in) :: dr2d,ds2d, sr2d, ss2d

    real, dimension(Np,Np), intent(inout) :: stfrr,stfss,stfrs,stfsr

    stfrr = matmul( transpose(dr2d), sr2d )
    stfss = matmul( transpose(ds2d), ss2d )

    stfrs = matmul( transpose(dr2d), ss2d )
    stfsr = matmul( transpose(ds2d), sr2d ) ! StfSR' = StfRS

  end subroutine stiffmat
  
  !> builds the local operators in 2d acting on the reference element
  !> @param [inout] lo2d is an object containing all of these operators
  subroutine buildlocopt2d(k,lo2d)
    integer,        intent(in)    :: k
    type(locopt2d), intent(inout) :: lo2d

    integer :: dim, faces

    if ( lo2d%init.eqv..true. ) then
       write(6,*) "error in buildlocopt2d: lo2d is already initialized"
       stop
    endif
    
    dim     = 2
    faces   = 3
    lo2d%k  = k
    lo2d%Np = (k+1)*(k+2)/2
    lo2d%Npfaces = (k+1)
    
    allocate( lo2d%refxy( lo2d%Np, dim ) )
    allocate( lo2d%nmask2d( lo2d%Npfaces, faces ) )
    
    lo2d%refxy   = nodes2d(k)
    lo2d%nmask2d = nodesmask2d(k)

    allocate( lo2d%va2d( lo2d%Np, lo2d%Np ) )
    allocate( lo2d%grva2dr( lo2d%Np, lo2d%Np ) )
    allocate( lo2d%grva2ds( lo2d%Np, lo2d%Np ) )

    allocate( lo2d%dr2d( lo2d%Np, lo2d%Np ) )
    allocate( lo2d%ds2d( lo2d%Np, lo2d%Np ) )

    allocate( lo2d%m2d( lo2d%Np, lo2d%Np ) )
    allocate( lo2d%sr2d( lo2d%Np, lo2d%Np ) )
    allocate( lo2d%ss2d( lo2d%Np, lo2d%Np ) )

    allocate( lo2d%stfrr ( lo2d%Np, lo2d%Np ) )
    allocate( lo2d%stfss( lo2d%Np, lo2d%Np ) )
    allocate( lo2d%stfrs( lo2d%Np, lo2d%Np ) )
    allocate( lo2d%stfsr( lo2d%Np, lo2d%Np ) )
    
    lo2d%va2d = vandermonde2d(k, lo2d%refxy)
    call dmat(k, lo2d%Np, lo2d%refxy, lo2d%va2d, lo2d%dr2d, lo2d%ds2d)

    lo2d%m2d = massmat(lo2d%Np, lo2d%va2d)

    call smat(lo2d%Np, lo2d%m2d, lo2d%dr2d, lo2d%ds2d, lo2d%sr2d, lo2d%ss2d)
    call stiffmat(lo2d%Np, lo2d%dr2d, lo2d%ds2d, lo2d%sr2d, lo2d%ss2d, &
         lo2d%stfrr, lo2d%stfss, lo2d%stfrs, lo2d%stfsr)

    lo2d%init = .true.
  end subroutine buildlocopt2d
  
  
end module localopt2d
