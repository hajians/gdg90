  !> @file dgfluxspace.f90


module dgfluxspace2d
  use dgspace2d
  implicit none

  type dgflux2d
     logical :: init = .false.
     
     type(dg2d), allocatable, dimension(:) :: cp ! components of dgflux
  end type dgflux2d

contains

  !> given Th and k it builds a null element
  !! of DG flux finite element space.
  !> @param [in] Th is the mesh object
  !> @param [in] k is the degree of polynomials
  function builddgflux2d(Th,lo1d, lo2d,k) result(sig)
    type(tri2d),         intent(in) :: Th
    type(locopt1d),      intent(in) :: lo1d
    type(locopt2d),      intent(in) :: lo2d
    integer,             intent(in) :: k

    type(dgflux2d) :: sig

    integer :: i
    allocate( sig%cp( Th%n_dim ) )

    do i=1, Th%n_dim
       sig%cp(i) = builddg2d(Th,lo1d,lo2d,k)
    end do

    sig%init = .true.
  end function builddgflux2d

  !>
  !!
  subroutine copydg2fluxdg(U, Sig)
    type(dg2d),     intent(in)    :: U
    type(dgflux2d), intent(inout) :: Sig

    integer :: dim
    if ( U%init .eqv. .false. .or. &
         Sig%init .eqv. .true. ) then
       write(6,*) "error in copydg2fluxdg",  Sig%init
       stop
    endif

    allocate( sig%cp( U%Th%n_dim ) )

    do dim=1, U%Th%n_dim
       call copydg2dg(U, sig%cp(dim) )
    enddo

    Sig%init = .true.
  end subroutine copydg2fluxdg
  
end module dgfluxspace2d
