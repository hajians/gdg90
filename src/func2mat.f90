  !> fctl2mat.f90
  !> @brief functional to matrices

module func2mat
  use dgopt2d
  implicit none

  !> abstract declaration
  abstract interface
     function bilinearformDG(U,V,DhU, DhV, trU,trV, trDu, trDv) result(out)
       import                 :: dg2d, dgflux2d, tr2d, trflux2d
       type(dg2d), intent(in) :: U, V
       type(dgflux2d)         :: DhU, DhV
       type(tr2d)             :: trU, trV
       type(trflux2d)         :: trDu, trDv
       real                   :: out
     end function bilinearformDG
  end interface

  abstract interface
     function specialformDG(U,V,DhU, DhV, trU,trV, trDu, trDv, index) result(out)
       import                 :: dg2d, dgflux2d, tr2d, trflux2d
       type(dg2d), intent(in) :: U, V
       type(dgflux2d)         :: DhU, DhV
       type(tr2d)             :: trU, trV
       type(trflux2d)         :: trDu, trDv
       real                   :: out
       integer, intent(in)    :: index(:)
     end function specialformDG
  end interface

  abstract interface
     function rhsDG(U) result(out)
       import :: dg2d
       type(dg2d), intent(in) :: U
       real                   :: out
     end function rhsDG
  end interface
contains

  !> build the right hand side vector
  subroutine rhsdg2vec(rhs_DG, Vorg, vec)
    procedure(rhsDG), pointer        :: rhs_DG
    type(dg2d), intent(in)           :: Vorg
    real, allocatable, intent(inout) :: vec(:)

    type(dg2d) :: V

    integer :: i, el, dof

    if ( Vorg%init.eqv..false. .or. &
         allocated( vec ).eqv..true. ) then
       write(6,*) "error in rhsdg2vec:", Vorg%init, allocated( vec )
       stop
    endif

    call copydg2dg(Vorg,V)

    allocate( vec( V%n_dof ) )
    vec   = 0.0

    V%dof = 0.0
    V%binary = .true.

    do el=1, V%Th%n_el
       V%elnotzero = el
       do i=1, V%Np
          dof = V%el2dof(el,i)         
          V%dof(dof) = 1.0
          vec( dof ) = rhs_DG(V)
          V%dof(dof) = 0.0
       enddo
    enddo

    V%binary    = .false.
    V%elnotzero = 0

  end subroutine rhsdg2vec

  !> build the matrix of the bilinear form
  subroutine bilindg2dmat(bilindg2d, Uorg, Vorg, mat)
    procedure(bilinearformDG), pointer    :: bilindg2d
    type(dg2d), intent(in)   :: Uorg, Vorg
    real, allocatable, dimension(:,:), intent(inout) :: mat

    integer :: i, j, el, face, Nface, adjel

    type(dg2d)     :: U, V
    type(tr2d)     :: trU, trV
    type(dgflux2d) :: DhU, DhV
    type(trflux2d) :: trDu, trDv

    integer, allocatable :: dofi(:), dofj(:)

    if ( Uorg%init .eqv. .false. .or. &
         Vorg%init .eqv. .false. .or. &
         allocated( mat )          &
         ) then
       write(6,*) "error in bilindg2dmat"
       stop
    endif

    call copydg2dg(Uorg,U)
    call copydg2dg(Vorg,V)
    call copydg2fluxdg(U, DhU)
    call copydg2fluxdg(V, DhV)

    call buildtr2d(U%Th, U%lo1d, U%k, trU)
    call buildtr2d(V%Th, V%lo1d, V%k, trV)

    call buildtrflux2d(U%Th, U%lo1d, U%k, trDu)
    call buildtrflux2d(V%Th, V%lo1d, V%k, trDv)

    call link(U, trU);    call link(V, trV)
    call link(DhU, trDu); call link(DhV, trDv)

    allocate( mat(U%n_dof, V%n_dof) )
    allocate( dofi(U%Np), dofj(U%Np) )

    U%dof = 0.0; V%dof = 0.0; mat = 0.0

    U%binary = .true.
    V%binary = .true.
    trU%binary = .true.
    trV%binary = .true.

    Nface = size( U%Th%el2el, 2 )

    do el=1, U%Th%n_el
       !
       dofi = U%el2dof(el,:)
       U%elnotzero = el
       trU%facenotzero = U%Th%el2face( el, : )
       ! element contribution
       dofj = U%el2dof(el,:)
       V%elnotzero = el
       trV%facenotzero = V%Th%el2face( el, : )
       !
       do i=1, U%Np
          U%dof( dofi(i) ) = 1.0
          do j=1, U%Np
             V%dof( dofj(j) ) = 1.0
             mat( dofi(i), dofj(j) ) =  &
                  bilindg2d(U,V,DhU, DhV,trU,trV, trDu, trDv)
             V%dof( dofj(j) ) = 0.0
          enddo
          U%dof( dofi(i) ) = 0.0
       enddo
       !
       ! neighbors contribution
       !
       do face=1, Nface
          adjel = U%Th%el2el(el,face)
          if (adjel.ne.-1) then
             V%elnotzero = adjel
             trV%facenotzero = V%Th%el2face( adjel, : )
             !
             dofj = U%el2dof(adjel,:)
             do i=1, U%Np
                U%dof( dofi(i) ) = 1.0
                do j=1, U%Np
                   V%dof( dofj(j) ) = 1.0
                   mat( dofi(i), dofj(j) ) =  &
                        bilindg2d(U,V,DhU, DhV,trU,trV, trDu, trDv)
                   V%dof( dofj(j) ) = 0.0
                enddo
                U%dof( dofi(i) ) = 0.0
             enddo
          endif
       enddo
       !
       !
       if ( mod(el,100).eq.0 ) then
          write(6,*) el
       endif
    enddo

    U%binary = .false.
    U%elnotzero = 0
    V%binary = .false.
    V%elnotzero = 0

  end subroutine bilindg2dmat

  subroutine specialdg2mat(formdg2d, U, V, IFindex, mat)
    procedure(specialformDG), pointer    :: formdg2d
    type(dg2d)               :: U, V
    integer, intent(in) :: IFindex(:)
    real, allocatable, dimension(:,:), intent(inout) :: mat

    integer :: i, j, NIFedge, ifedge, ifel, el, adjel, edge

    type(tr2d)     :: trU, trV
    type(dgflux2d) :: DhU, DhV
    type(trflux2d) :: trDu, trDv

    integer, allocatable :: dofi(:), dofj(:)

    if ( U%init .eqv. .false. .or. &
         V%init .eqv. .false. .or. &
         allocated( mat )          &
         ) then
       write(6,*) "error in specialdg2mat"
    endif

    call copydg2fluxdg(U, DhU)
    call copydg2fluxdg(V, DhV)

    call buildtr2d(U%Th, U%lo1d, U%k, trU)
    call buildtr2d(V%Th, V%lo1d, V%k, trV)

    call buildtrflux2d(U%Th, U%lo1d, U%k, trDu)
    call buildtrflux2d(V%Th, V%lo1d, V%k, trDv)

    call link(U, trU);    call link(V, trV)
    call link(DhU, trDu); call link(DhV, trDv)

    allocate( mat(U%n_dof, V%n_dof) )
    allocate( dofi(U%Np), dofj(U%Np) )

    U%dof = 0.0; V%dof = 0.0; mat = 0.0

    NIFedge = size( IFindex )

    do ifedge=1, NIFedge
       edge = IFindex(ifedge)
       do ifel=1, 2               ! two element per interface edge
          el    =  U%Th%edge2el(edge,ifel)
          adjel =  U%Th%edge2el(edge,3-ifel)
          dofi = U%el2dof(el,:)
          ! element contribution
          dofj = U%el2dof(el,:)
          do i=1, U%Np
             U%dof( dofi(i) ) = 1.0
             do j=1, U%Np
                V%dof( dofj(j) ) = 1.0
                mat( dofi(i), dofj(j) ) =  &
                     formdg2d(U,V,DhU, DhV, trU,trV, trDu, trDv, IFindex)
                V%dof( dofj(j) ) = 0.0
             enddo
             U%dof( dofi(i) ) = 0.0
          enddo
          !
       enddo
    enddo

    ! do i=1, U%n_dof
    !    U%dof(i) = 1.0
    !    do j=1, V%n_dof
    !       V%dof(j) = 1.0
    !       !
    !       mat(i,j) = formdg2d(U,V,DhU, DhV, trU,trV, trDu, trDv, index)
    !       !
    !       V%dof(j) = 0.0
    !    enddo
    !    !write(6,*) i
    !    U%dof(i) = 0.0
    ! enddo
  end subroutine specialdg2mat


end module func2mat
