module functional2mat

  use sparsemat
  use alldgspace2d
  implicit none

  !> abstract declaration
  abstract interface
     subroutine bilinearformPrimal(U,V,out)
       import :: primaldg2d

       type(primaldg2d), intent(inout) :: U, V
       real,             intent(inout) :: out
     end subroutine bilinearformPrimal
  end interface

  !> abstract declaration
  abstract interface
     function L2prod(U) result(out)
       import :: dg2d
       type(dg2d), intent(in) :: U
       real                   :: out
     end function L2prod
  end interface

contains

  !> convert a functional to a vector
  subroutine L2prod2vec(rhs_DG, Vorg, vec)
    procedure(L2prod), pointer       :: rhs_DG
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

  end subroutine L2prod2vec

  !> convert a bilinear form with primal spaces U, V
  !! into a sparse matrix
  subroutine bilinpr2spmat(bilinpr,U,V,spA)
    procedure(bilinearformPrimal), pointer :: bilinpr

    type(primaldg2d), intent(inout) :: U, V

    type(spmat), intent(inout) :: spA

    integer :: el, adjel, locface
    integer :: i, j

    integer :: est_nonzero

    integer, allocatable, dimension(:) :: dofi, dofj

    real :: out

    if ( U%init.eqv..false. .or. &
         V%init.eqv..false. .or. &
         spA%init.eqv..true. ) then
       write(6,*) "error in bilinpr2mat"
       stop
    endif

    est_nonzero = 3*(U%sc%Np**2)
    
    call buildspmat(spA, U%sc%n_dof, V%sc%n_dof, est_nonzero)

    allocate( dofi( U%sc%Np ), dofj( V%sc%Np ) )

    U%sc%dof = 0.0
    V%sc%dof = 0.0

    out = 0.0

    do el=1, U%sc%Th%n_el
       !
       dofi = U%sc%el2dof(el,:)
       dofj = V%sc%el2dof(el,:)
       !
       do i=1, U%sc%Np
          U%sc%dof(dofi(i)) = 1.0
          do j=1, V%sc%Np
             V%sc%dof( dofj(j) ) = 1.0
             !
             call bilinpr(U,V,out)
             call insert2spmat(spA, dofi(i), dofj(j), out)
             !
             V%sc%dof( dofj(j) ) = 0.0
          enddo
          U%sc%dof( dofi(i) ) = 0.0
       enddo
       !
       ! adjacent elements contribution
       do locface=1, size(U%sc%Th%el2el,2)

          adjel = U%sc%Th%el2el(el,locface)
          if (adjel.ne.-1) then
             !
             dofj = V%sc%el2dof(adjel,:)
             do i=1, U%sc%Np
                U%sc%dof(dofi(i)) = 1.0
                do j=1, V%sc%Np
                   V%sc%dof( dofj(j) ) = 1.0
                   !
                   call bilinpr(U,V,out)
                   call insert2spmat(spA, dofi(i), dofj(j), out)
                   !
                   V%sc%dof( dofj(j) ) = 0.0
                enddo
                U%sc%dof( dofi(i) ) = 0.0
             enddo
             !
          endif
       enddo
    enddo

  end subroutine bilinpr2spmat

  !> convert a bilinear form with primal spaces U, V
  !! into a matrix
  subroutine bilinpr2mat(bilinpr,U,V,mat)
    procedure(bilinearformPrimal), pointer :: bilinpr

    type(primaldg2d), intent(inout) :: U, V

    real, allocatable, dimension(:,:), intent(inout) :: mat

    integer :: el, adjel, locface
    integer :: i, j

    integer, allocatable, dimension(:) :: dofi, dofj

    real :: out

    if ( U%init.eqv..false. .or. &
         V%init.eqv..false. .or. &
         allocated(mat) ) then
       write(6,*) "error in bilinpr2mat"
       stop
    endif

    allocate( mat(U%sc%n_dof, V%sc%n_dof) )

    allocate( dofi( U%sc%Np ), dofj( V%sc%Np ) )

    U%sc%dof = 0.0
    V%sc%dof = 0.0

    mat = 0.0
    out = 0.0

    do el=1, U%sc%Th%n_el
       !
       dofi = U%sc%el2dof(el,:)
       dofj = V%sc%el2dof(el,:)
       !
       do i=1, U%sc%Np
          U%sc%dof(dofi(i)) = 1.0
          do j=1, V%sc%Np
             V%sc%dof( dofj(j) ) = 1.0
             !
             call bilinpr(U,V,out)
             mat(dofi(i), dofj(j)) = out
             !
             V%sc%dof( dofj(j) ) = 0.0
          enddo
          U%sc%dof( dofi(i) ) = 0.0
       enddo
       !
       ! adjacent elements contribution
       do locface=1, size(U%sc%Th%el2el,2)

          adjel = U%sc%Th%el2el(el,locface)
          if (adjel.ne.-1) then
             !
             dofj = V%sc%el2dof(adjel,:)
             do i=1, U%sc%Np
                U%sc%dof(dofi(i)) = 1.0
                do j=1, V%sc%Np
                   V%sc%dof( dofj(j) ) = 1.0
                   !
                   call bilinpr(U,V,out)
                   mat(dofi(i), dofj(j)) = out
                   !
                   V%sc%dof( dofj(j) ) = 0.0
                enddo
                U%sc%dof( dofi(i) ) = 0.0
             enddo
             !
          endif
       enddo
    enddo

  end subroutine bilinpr2mat
end module functional2mat
