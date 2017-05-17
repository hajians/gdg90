  !> @file simplex2dp.f90

module simplex2dp
  use jacobip
  implicit none

contains

  !> maps r,s coordinates into a,b variables
  subroutine rs2ab(r,s,a,b)
    real, dimension(:),       intent(in)    :: r, s
    real, dimension(size(r)), intent(inout) :: a, b

    integer :: l, Np

    Np = size(r)

    do l=1, Np
       if ( abs(s(l) - 1.0)< 10.0**(-8) ) then
          a(l) = -1.0
       else
          a(l) = 2.0 * (1.0 + r(l))/(1.0 - s(l)) - 1.0
       end if
    end do

    b = s

  end subroutine rs2ab

  !> Simplex2DP gives the evaluation of an orthonormalized 
  !! basis function defined on the reference element.
  !> @param [in] coorxy is a list of coordinates at which the basis function
  !> @param [in] i index
  !> @param [in] j index
  !> @param [out] psi 
  function simpl2dp(coorxy,i,j) result(psi)
    integer, intent(in) :: i,j
    real, dimension(:,:), intent(in) :: coorxy

    real, allocatable, dimension(:) :: psi

    real, allocatable, dimension(:) :: r,s, a,b
    integer :: Np

    Np = size(coorxy,1)

    allocate( psi(Np) )
    allocate( r(Np), s(Np), a(Np), b(Np) )

    psi = 0.0
    r = coorxy(:,1); s = coorxy(:,2)

    call rs2ab(r,s,a,b)

    psi = sqrt(2.0) * JacP(a,0.0,0.0,i) * &
         JacP(b,2.0*i+1.0,0.0,j) * (1.0-b)**i

  end function simpl2dp

  !> GradSimplex2DP computes the D_r and D_s matrices.
  !! they are derivatives of \Psi evaluated at Node2D points.
  !> @param [in] a are arays that contains the transformed Node2D points.
  !> @param [in] b are arays that contains the transformed Node2D points.
  !> @param [in] i are integers that identify the orthonormal basis function.
  !> @param [in] j are integers that identify the orthonormal basis function.
  !> @param [out] dmodedr
  !> @param [out] dmodeds
  subroutine gradsimpl2dp(a,b,i,j,dmodedr,dmodeds)
    real,    intent(in), dimension(:) :: a,b
    integer, intent(in)               :: i,j

    real, intent(inout), dimension( size(a) ) :: dmodedr, dmodeds

    real, dimension( size(a) ) :: tmp
    real, dimension( size(a) ) :: fa, dfa
    real, dimension( size(a) ) :: gb, dgb

    if (size(a).ne.size(b)) then
       write(6,*) "error"
       stop
    end if

    fa  = JacP(a,0.,0.,i)
    dfa = GradJP(a,0.,0.,i)

    gb  = JacP(b, 2.*i+1., 0., j);
    dgb = GradJP(b, 2.*i+1., 0., j);

    dmodedr = dfa*gb            ! r-derivative
    if (i>0) then
       dmodedr = dmodedr * ( (0.5*(1.0-b))**(i-1.0) )
    end if

    dmodeds = dfa*(gb*(0.5*(1.0+a))) ! s-derivative
    if (i>0) then
       dmodeds = dmodeds*( (0.5*(1.0-b))**(i-1.0) )
    end if

    tmp = dgb*( (0.5*(1.0-b))**i ) ! 2nd term in d/ds
    if (i>0) then
       tmp = tmp - 0.5*i*gb * ( (0.5*(1.0-b))**(i-1.0) )
    end if
    dmodeds = dmodeds + fa*tmp

    dmodeds = 2.0**(i+0.5) * dmodeds
    dmodedr = 2.0**(i+0.5) * dmodedr

  end subroutine gradsimpl2dp

  !> Vandermonde2D gives the Vandermonde matrix evaluated at nodes
  !> @param [in] coorXY for polynomials of degree k.
  function vandermonde2d(k, coorxy) result(va2d)
    real, dimension(:,:), intent(in) :: coorxy
    integer,              intent(in) :: k

    real, dimension( size(coorxy,1) , size(coorxy,1) ) :: va2d

    integer :: Np, col, i, j

    Np = size(coorxy, 1)

    col = 1

    do i=0, k
       do j=0, k-i
          va2d(:,col) = simpl2dp(coorxy,i,j)
          col = col+1
       end do
    end do

  end function vandermonde2d

  !> GradVandermonde2D gives the derivatives orthonormal basis
  !! evaluated at Nodes2D with respect to r, s.
  !! Vandermonde2D gives the Vandermonde matrix evaluated at nodes
  !> @param [in] coorXY for polynomials of degree k. 
  subroutine gradvandermonde2d(k, coorxy, grva2dr, grva2ds)
    real, dimension(:,:), intent(in) :: coorxy
    integer,              intent(in) :: k

    real, dimension( size(coorxy,1), size(coorxy,1) ) :: grva2dr, grva2ds

    integer :: Np, col, i, j
    real, dimension( size(coorxy,1) ) :: r, s, a, b

    Np = size(coorxy, 1)
    
    if ( size(coorxy(:,1)).ne.size(coorxy(:,2)) ) then
       write(6,*) "error"
       stop
    end if

    r = coorxy(:,1)
    s = coorxy(:,2)

    call rs2ab(r,s,a,b)

    col = 1

    do i=0, k
       do j=0, k-i
          call gradsimpl2dp(a,b,i,j, grva2dr(:,col), grva2ds(:,col) )
          col = col+1
       end do
    end do


  end subroutine gradvandermonde2d

end module simplex2dp
