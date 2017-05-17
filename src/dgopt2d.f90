!> @file dgopt2d.f90 contains operators which act on the dgspace2d

module dgopt2d
  use dgspace2d
  use dgfluxspace2d
  use brokentracespace2d
  use singletraceopt2d
  use interpolate2d
  use mylinalgebra

  implicit none

  !> interface for scalar product and vector dot product
  interface dot0Th
     module procedure dot0sc, dot0vec
  end interface dot0Th

  !> interface for \nabla_h
  interface Dh
     module procedure Dhdg2d
  end interface Dh

  !> interface for linkning trace and fem spaces
  interface link
     module procedure linkdg2d, linkfluxdg2d
  end interface link

  !> interface for trace
  interface trace
     module procedure tracedg2d, tracefluxdg2d
  end interface trace

  !> interface for average
  interface average
     module procedure avedg2d, aveflux2d
  end interface average

  !> interface for jump
  interface jump
     module procedure jumpdg2d, jumpflux2d
  end interface jump

  !> interface for 1d edge integral
  interface dot01d
     module procedure dot1dsc, dot1dvec
  end interface dot01d

  !> interface operator for dot
  interface operator (.dot.)
     module procedure paramdot1dsc, paramdot1dvec, consdot1dsc, consdot1dvec
  end interface operator (.dot.)

  interface operator (.add.)
     module procedure addtrace, addfluxtrace
  end interface operator (.add.)

  ! private
  private dot0sc
  private dot0vec

  private Dhdg2d

  private linkdg2d, linkfluxdg2d
  private tracedg2d, tracefluxdg2d
  private avedg2d, aveflux2d
  private jumpdg2d, jumpflux2d
  private dot1dsc, dot1dvec

  private paramdot1dsc, paramdot1dvec
  private consdot1dsc, consdot1dvec

  private addtrace, addfluxtrace
contains

  subroutine linkdg2d(V,trV)
    type(dg2d), intent(in)    :: V
    type(tr2d), intent(inout) :: trV

    logical, allocatable, dimension(:) :: faceflag
    integer :: el, adjel, face, adjface, Nface, globf, adjglobf
    integer, allocatable, dimension(:) :: facedof, adjfacedof

    integer :: i
    real    :: eps = 10.0**(-9)

    if (trV%init.eqv..false. .or. V%init.eqv..false.) then
       write(6,*) "error in linkdg2d: trV or V is not initialized", trV%init, V%init
       stop
    elseif ( V%k .ne. trV%k ) then
       write(6,*) "error in linkdg2d: trV and V do not have same polynomial deg."
       stop
    endif


    allocate( faceflag( V%Th%n_faces ) )
    faceflag = .false.
    allocate( facedof( trV%NpF ) )
    allocate( adjfacedof( trV%NpF ) )
    Nface = size( V%Th%el2face, 2 ) ! number of faces per element

    do el=1, V%Th%n_el
       do face=1, Nface             
          globf    = V%Th%el2face(el,face)
          !
          if ( faceflag(globf).eqv..false. ) then ! if not visited

             facedof = V%lo2d%nmask2d(:,face)

             trV%face2idof(globf,:) = V%el2dof(el,facedof)
             !trV%face2dof(globf,:)  = V%dof( trV%face2idof(globf,:) )
             faceflag(globf) = .true.

             adjglobf = V%Th%face2face(globf)
             !write(6,*) "globf", globf, "adj", adjglobf
             if (adjglobf.ne.-1) then ! if globf is not on the boundary
                !
                adjel   = V%Th%el2el(el,face)
                adjface = V%Th%face2locf(adjglobf)
                adjfacedof = V%lo2d%nmask2d(:,adjface)
                trV%face2idof(adjglobf,:) = V%el2dof(adjel,adjfacedof)
                if ( &
                     abs(V%dof2corr( trV%face2idof(globf,1), 1 ) - &
                     V%dof2corr( trV%face2idof(adjglobf,1), 1 )).le.eps &
                     .and. &
                     abs(V%dof2corr( trV%face2idof(globf,1), 2 ) - &
                     V%dof2corr( trV%face2idof(adjglobf,1), 2 )).le.eps &
                     ) then
                   trV%face2idof(adjglobf,:) = V%el2dof(adjel,adjfacedof)
                else

                   adjfacedof = adjfacedof( (/ (i,i=trV%NpF,1,-1) /) )
                   trV%face2idof(adjglobf,:) = V%el2dof(adjel,adjfacedof)

                endif

                !trV%face2dof(adjglobf,:)  = V%dof( trV%face2idof(adjglobf,:) )
                !
                faceflag(adjglobf) = .true.
             endif

          endif
          !
       enddo
    enddo

    trV%nx = V%nx
    trV%ny = V%ny

    trV%link = .true.
  end subroutine linkdg2d

  subroutine linkfluxdg2d(Sig, trSig)
    type(dgflux2d), intent(in)    :: Sig
    type(trflux2d), intent(inout) :: trSig

    integer :: dim

    if (trSig%init.eqv..false. .or. Sig%init.eqv..false.) then
       write(6,*) "error in linkfluxdg2d: trSig or Sig is not initialized", &
            trSig%init, Sig%init
       stop
    elseif ( Sig%cp(1)%k .ne. trSig%cp(1)%k ) then
       write(6,*) "error in linkfluxdg2d: trSig and Sig do not have same polynomial deg."
       stop
    endif

    do dim=1, Sig%cp(1)%Th%n_dim
       call linkdg2d( Sig%cp(dim), trSig%cp(dim) )
    enddo

    trSig%link = .true.
  end subroutine linkfluxdg2d


  subroutine tracedg2d(V,trV)
    type(dg2d), intent(in)    :: V
    type(tr2d), intent(inout) :: trV

    !integer, optional, intent(in) :: edge_index(:)

    integer :: globf

    !integer :: len_index, id
    integer :: face, Nface

    !real, allocatable :: locV(:)
    !real :: eps = 10.0**(-9)

    if (trV%init.eqv..false. .or. V%init.eqv..false.) then
       write(6,*) "error in tracedg2d: trV or V is not initialized"
       stop
    elseif ( V%k .ne. trV%k ) then
       write(6,*) "error in tracedg2d: trV and V do not have same polynomial deg."
       stop
    elseif (trV%link.eqv..false.) then
       write(6,*) "error in tracedg2d: trV is not linked to anyspace"
       stop
    endif

    if ( V%binary.eqv..true. ) then
       Nface = 3                ! three face per element
       trV%face2dof = 0.0
       trV%binary = .true.
       trV%facenotzero = V%Th%el2face( V%elnotzero, : )
       do face=1, Nface
          globf = V%Th%el2face( V%elnotzero, face )
          trV%face2dof(globf,:) = V%dof( trV%face2idof(globf,:) )
       enddo
    else
       do globf=1, V%Th%n_faces
          trV%face2dof(globf,:) = V%dof( trV%face2idof(globf,:) )
       enddo
    endif

  end subroutine tracedg2d

  !> tracefluxdg2d 
  subroutine tracefluxdg2d(Sig, trSig)
    type(dgflux2d), intent(in)    :: Sig
    type(trflux2d), intent(inout) :: trSig
    !integer, optional, intent(in) :: edge_index(:)

    integer :: dim

    if (trSig%init.eqv..false. .or. Sig%init.eqv..false.) then
       write(6,*) "error in tracefluxdg2d: trSig or Sig is not initialized"
       stop
    elseif ( Sig%cp(1)%k .ne. trSig%cp(1)%k ) then
       write(6,*) "error in tracefluxdg2d: trSig and Sig do not have same polynomial deg."
       stop
    elseif (trSig%link.eqv..false.) then
       write(6,*) "error in tracedg2d: trSig is not linked to anyspace"
       stop       
    endif

    do dim=1, Sig%cp(1)%Th%n_dim
       ! if ( present(edge_index) ) then
       !    call tracedg2d( Sig%cp(dim), trSig%cp(dim), edge_index )
       ! else
       call tracedg2d( Sig%cp(dim), trSig%cp(dim) )
       ! endif
    enddo

  end subroutine tracefluxdg2d

  !> computes average of trV and returns an object
  !! which is single valued on edges
  function avedg2d(trV) result(ave)
    type(tr2d), intent(in) :: trV
    type(tr2d)             :: ave

    integer :: globf, adjglobf

    integer :: face, Nface

    logical, dimension( trV%Th%n_faces ) :: faceflag 

    ! call buildtr2d(trV%Th,trV%lo1d,trV%k,ave)
    ! ave%face2idof = trV%face2idof
    ! ave%link = .true.
    call copytr2tr(trV, ave)

    ! allocate( faceflag( trV%Th%n_faces ) )

    faceflag = .false.
    if ( trV%binary .eqv. .true. ) then
       Nface = size( trV%facenotzero )
       ave%face2dof = 0.0
       !
       do face=1, Nface
          globf = trV%facenotzero(face)
          adjglobf = trV%Th%face2face(globf)
          !
          if ( faceflag(globf).eqv. .false. ) then
             !
             if (adjglobf.ne.-1) then
                ave%face2dof(globf,:) = &
                     0.5 * ( trV%face2dof(globf,:)    + &
                     trV%face2dof(adjglobf,:)   &
                     )
                ave%face2dof(adjglobf,:) = ave%face2dof(globf,:)
                faceflag(globf) = .false.
             else
                ave%face2dof(globf,:) = trV%face2dof(globf,:)
             endif
             !
             faceflag(globf) = .true.
          endif
          !
       enddo
       !
    else                        ! if trV is not binary
       do globf=1, trV%Th%n_faces
          adjglobf = trV%Th%face2face(globf)
          if ( faceflag(globf).eqv. .false. ) then
             !
             if (adjglobf.ne.-1) then
                ave%face2dof(globf,:) = &
                     0.5 * ( trV%face2dof(globf,:)    + &
                     trV%face2dof(adjglobf,:)   &
                     )
                ave%face2dof(adjglobf,:) = ave%face2dof(globf,:)
                faceflag(globf) = .false.
             else
                ave%face2dof(globf,:) = trV%face2dof(globf,:)
             endif
             !
             faceflag(globf) = .true.
          endif
       enddo
    endif

  end function avedg2d

  !> computes average of trSig and returns an object
  !! which is single valued on edges
  function aveflux2d(trSig) result(ave)
    type(trflux2d), intent(in) :: trSig
    type(trflux2d)             :: ave

    integer :: dim

    !call buildtrflux2d(trSig%cp(1)%Th,trSig%cp(1)%lo1d,trSig%cp(1)%k,ave)
    call copytr2trflux( trSig%cp(1), ave )

    do dim=1, trSig%cp(1)%Th%n_dim
       ave%cp(dim) = avedg2d(trSig%cp(dim))
    enddo

  end function aveflux2d


  !> computes the jump of trV and returns an object
  !! which is single valued on edges
  function jumpdg2d(trV) result(jump)
    type(tr2d), intent(in) :: trV
    type(trflux2d)         :: jump

    integer :: dim, globf, adjglobf

    integer :: face, Nface

    logical, dimension( trV%Th%n_faces ) :: faceflag 

    call copytr2trflux( trV, jump )


    !call buildtrflux2d(trV%Th,trV%lo1d,trV%k,jump)    
    ! do dim=1, trV%Th%n_dim
    !    jump%cp(dim)%face2idof = trV%face2idof
    !    jump%cp(dim)%link = .true.
    ! enddo
    faceflag = .false.

    if ( trV%binary .eqv. .true. ) then
       Nface = size( trV%facenotzero )
       dim = 1
       jump%cp(dim)%face2dof = 0.0
       dim = 2
       jump%cp(dim)%face2dof = 0.0
       !
       do face=1, Nface
          globf = trV%facenotzero(face)
          !
          if ( faceflag(globf) .eqv. .false.) then
             !
             adjglobf = trV%Th%face2face(globf)
             if (adjglobf.ne.-1) then
                dim = 1
                jump%cp(dim)%face2dof(globf,:) = &
                     trV%nx(globf)*trV%face2dof(globf,:) + &
                     trV%nx(adjglobf)*trV%face2dof(adjglobf,:)
                !
                jump%cp(dim)%face2dof(adjglobf,:) = &
                     jump%cp(dim)%face2dof(globf,:)
                !
                dim = 2
                jump%cp(dim)%face2dof(globf,:) = &
                     trV%ny(globf)*trV%face2dof(globf,:) + &
                     trV%ny(adjglobf)*trV%face2dof(adjglobf,:)
                !
                jump%cp(dim)%face2dof(adjglobf,:) = &
                     jump%cp(dim)%face2dof(globf,:)
                faceflag( adjglobf ) = .true.
             else
                dim = 1
                jump%cp(dim)%face2dof(globf,:) = &
                     trV%nx(globf)*trV%face2dof(globf,:) 
                dim = 2
                jump%cp(dim)%face2dof(globf,:) = &
                     trV%ny(globf)*trV%face2dof(globf,:) 
             endif
             !
             faceflag( globf ) = .true.
          endif
          !
       enddo
    else
       do globf=1, trV%Th%n_faces

          if ( faceflag(globf) .eqv. .false.) then
             !
             adjglobf = trV%Th%face2face(globf)
             if (adjglobf.ne.-1) then
                dim = 1
                jump%cp(dim)%face2dof(globf,:) = &
                     trV%nx(globf)*trV%face2dof(globf,:) + &
                     trV%nx(adjglobf)*trV%face2dof(adjglobf,:)
                !
                jump%cp(dim)%face2dof(adjglobf,:) = &
                     jump%cp(dim)%face2dof(globf,:)
                !
                dim = 2
                jump%cp(dim)%face2dof(globf,:) = &
                     trV%ny(globf)*trV%face2dof(globf,:) + &
                     trV%ny(adjglobf)*trV%face2dof(adjglobf,:)
                !
                jump%cp(dim)%face2dof(adjglobf,:) = &
                     jump%cp(dim)%face2dof(globf,:)
                faceflag( adjglobf ) = .true.
             else
                dim = 1
                jump%cp(dim)%face2dof(globf,:) = &
                     trV%nx(globf)*trV%face2dof(globf,:) 
                dim = 2
                jump%cp(dim)%face2dof(globf,:) = &
                     trV%ny(globf)*trV%face2dof(globf,:) 
             endif
             !
             faceflag( globf ) = .true.
          endif
       enddo
       !
    end if

  end function jumpdg2d

  !> computes the jump of trV and returns an object
  !! which is single valued on edges  
  function jumpflux2d(trSig) result(jump)
    type(trflux2d), intent(in) :: trSig
    type(tr2d)                 :: jump

    type(trflux2d) :: dummy
    integer        :: dim

    !call buildtrflux2d(trSig%cp(1)%Th,trSig%cp(1)%lo1d,trSig%cp(1)%k, dummy)

    !call buildtr2d(trSig%cp(1)%Th,trSig%cp(1)%lo1d,trSig%cp(1)%k, jump)
    !jump%face2idof = trSig%cp(1)%face2idof

    call copytr2tr( trSig%cp(1), jump )

    jump%face2dof  = 0.0

    call copytr2trflux( jump, dummy )

    do dim=1, trSig%cp(1)%Th%n_dim
       dummy =  jumpdg2d( trSig%cp(dim) )
       jump%face2dof = jump%face2dof + dummy%cp(dim)%face2dof
    enddo

  end function jumpflux2d

  !> compute Tau \dot \NOR can write it into trV
  subroutine normalproduct(trSig, trV)
    type(trflux2d), intent(in) :: trSig
    type(tr2d)                 :: trV

    integer :: globf, adjglobf

    logical, dimension( trV%Th%n_faces ) :: faceflag 

    faceflag = .false.

    trV%face2dof = 0.0

    do globf=1, trV%Th%n_faces

       if ( faceflag(globf) .eqv. .false.) then
          !
          adjglobf = trV%Th%face2face(globf)
          if (adjglobf.ne.-1) then

             trV%face2dof(globf,:) = &
                  trSig%cp(1)%nx(globf) * trSig%cp(1)%face2dof(globf,:) + &
                  trSig%cp(2)%ny(globf) * trSig%cp(2)%face2dof(globf,:)

             trV%face2dof(adjglobf,:) = &
                  trSig%cp(1)%nx(adjglobf) * trSig%cp(1)%face2dof(adjglobf,:) + &
                  trSig%cp(2)%ny(adjglobf) * trSig%cp(2)%face2dof(adjglobf,:)
          else
             trV%face2dof(globf,:) = &
                  trSig%cp(1)%nx(globf) * trSig%cp(1)%face2dof(globf,:) + &
                  trSig%cp(2)%ny(globf) * trSig%cp(2)%face2dof(globf,:)
          endif
          !
          faceflag( globf ) = .true.
       endif
    enddo

  end subroutine normalproduct

  !> computes the derivative of an object in dgopt2d
  function Dhdg2d(U) result(DhU)
    type(dg2d),        intent(in) :: U
    !integer, optional, intent(in) :: el_index(:)
    type(dgflux2d)             :: DhU

    real :: rx, ry, sx, sy
    real, dimension( U%Np, U%Np ) :: Dx, Dy

    integer :: el, i_cp
    !integer :: id, len_index
    !integer :: dim
    integer, dimension( U%Np ) :: locdof
    real,    dimension( U%Np ) :: locU
    real :: eps = 10.0**(-8)

    !DhU = builddgflux2d( U%Th, U%lo1d, U%lo2d, U%k )
    call copydg2fluxdg(U, DhU)

    if ( U%binary.eqv..true. ) then
       el = U%elnotzero
       !
       locdof = U%el2dof(el,:)
       locU   = U%dof( locdof )
       !
       rx = U%rx(1,el)
       ry = U%ry(1,el)
       sx = U%sx(1,el)
       sy = U%sy(1,el)

       Dx = rx * U%lo2d%dr2d + sx * U%lo2d%ds2d
       Dy = ry * U%lo2d%dr2d + sy * U%lo2d%ds2d

       locdof = U%el2dof(el,:)
       locU   = U%dof( locdof )

       i_cp   = 1
       DhU%cp(i_cp)%dof( locdof )  = matmul(Dx, locU )

       i_cp   = 2
       DhU%cp(i_cp)%dof( locdof )  = matmul(Dy, locU )

    else 
       do el=1, U%Th%n_el

          locdof = U%el2dof(el,:)
          locU   = U%dof( locdof )

          if ( sum(abs(locU)).ge.eps ) then
             rx = U%rx(1,el)
             ry = U%ry(1,el)
             sx = U%sx(1,el)
             sy = U%sy(1,el)

             Dx = rx * U%lo2d%dr2d + sx * U%lo2d%ds2d
             Dy = ry * U%lo2d%dr2d + sy * U%lo2d%ds2d

             locdof = U%el2dof(el,:)
             locU   = U%dof( locdof )

             i_cp   = 1
             DhU%cp(i_cp)%dof( locdof )  = matmul(Dx, locU )

             i_cp   = 2
             DhU%cp(i_cp)%dof( locdof )  = matmul(Dy, locU )
          else
          endif
       enddo
    endif
  end function Dhdg2d

  !> computes the L2-scalar product for two element of dgspace2d
  !> @param[in] U an element in dgspace2d
  !> @param[in] V an element in dgspace2d
  !> @param[out] dotv computes \int u v over Th
  function dot0sc(U,V)
    type(dg2d), intent(in) :: U, V
    !integer, optional, intent(in) :: el_index(:)
    real                   :: dot0sc

    real :: eps = 10.0**(-9)
    real :: locJ
    real, dimension(:), allocatable :: locU, locV

    integer :: el
    !integer :: id, len_index

    dot0sc = 0.0                ! init.

    ! cheking whether or not U and V belong to the same space. !COMPLETE
    if ( U%Th%n_dim.ne.V%Th%n_dim .or. &
         U%Th%n_el.ne.V%Th%n_el   .or. &
         U%k.ne.V%k                    &       
         ) then
       write(6,*) "error in dot0sc"
       stop
    end if
    !

    allocate( locU( U%Np ) )
    allocate( locV( V%Np ) )


    if ( (U%binary.eqv..true.) .and. (V%binary.eqv..true.) .and. &
         (U%elnotzero.eq.V%elnotzero) ) then
       el = U%elnotzero
       locU = U%dof( U%el2dof(el,:) )
       locV = V%dof( V%el2dof(el,:) )
       locJ = U%el2J(1,el)      ! assumed to have triangles

       dot0sc = dot0sc + &
            locJ * dot_product(locU,  matmul( U%lo2d%m2d, locV ) )
       !
    elseif ( (U%binary.eqv..true.) .and. (V%binary.eqv..true.) .and. &
         (U%elnotzero.ne.V%elnotzero) ) then
       !write(6,*) U%binary, V%binary, U%elnotzero, V%elnotzero
       dot0sc = 0.0
    else
       !
       do el=1, U%Th%n_el
          locU = U%dof( U%el2dof(el,:) )
          locV = V%dof( V%el2dof(el,:) )
          locJ = U%el2J(1,el)      ! assumed to have triangles

          if ( sum(abs(locU)).le.eps .or. &
               sum(abs(locV)).le.eps )    &
               then
          else
             dot0sc = dot0sc + &
                  locJ * dot_product(locU,  matmul( U%lo2d%m2d, locV ) )
             ! or V%lo2d%m2d
          end if
       end do
       !
    endif
  end function dot0sc

  function dot0vec(sig,tau)
    type(dgflux2d),    intent(in) :: sig, tau
    !integer, optional, intent(in) :: el_index(:)
    real                       :: dot0vec

    integer :: i, dim

    if ( sig%cp(1)%Th%n_dim .ne. tau%cp(1)%Th%n_dim ) then
       write(6,*) "error in dot0vec"
       stop
    end if

    dim = sig%cp(1)%Th%n_dim
    dot0vec = 0.0

    ! if ( present(el_index) ) then
    !    do i=1, dim
    !       dot0vec = dot0vec + dot0sc( sig%cp(i), tau%cp(i), el_index )
    !    enddo
    ! else
    do i=1, dim
       dot0vec = dot0vec + dot0sc( sig%cp(i), tau%cp(i) )
    enddo
    ! endif
  end function dot0vec

  !> dot product in 1D on a set of edges
  !> @param trU is supposed to be sing valued over edges
  !> @param trV
  function dot1dsc(trU, trV, edgeindex)
    type(tr2d), intent(in) :: trU, trV
    integer, optional, dimension(:), intent(in) :: edgeindex

    real :: dot1dsc

    real, dimension(:), allocatable    :: locU, locV
    integer :: dummy, n, globf, edge, n_edgeidx
    real    :: J1d
    integer, dimension(:), allocatable :: edgeidx

    real :: eps = 10.0**(-9)

    if ( trU%Th%n_dim .ne. trV%Th%n_dim .or.                     &
         trU%init .eqv. .false. .or. trV%init .eqv. .false. .or. &
         trU%link .eqv. .false. .or. trV%link .eqv. .false.      &
         ) then
       write(6,*) "error in dot1dsc", trU%link, trV%link
       stop
    end if

    if ( present(edgeindex) ) then
       n_edgeidx = size( edgeindex )
       allocate( edgeidx(n_edgeidx) )
       edgeidx   = edgeindex
    else
       n_edgeidx = trU%Th%n_edge
       allocate( edgeidx(n_edgeidx) )
       do dummy=1, n_edgeidx
          edgeidx(dummy) = dummy
       enddo
       !       edgeidx = (/ (n, n=1, n_edgeidx) /)       
    end if
    !
    dot1dsc = 0.0
    allocate( locU( trU%NpF ) )
    allocate( locV( trV%NpF ) )
    !
    do n=1, n_edgeidx
       edge  = edgeidx(n)
       globf = trU%Th%edge2face( edge, 1 ) ! trU, trV are single valued
       J1d   = trU%edge2J1d( edge )

       locU  = trU%face2dof( globf, : )
       locV  = trV%face2dof( globf, : )
       if ( sum(abs(locU)).le.eps .or. &
            sum(abs(locV)).le.eps )    &
            then
       else
          dot1dsc = dot1dsc + &
               J1d * dot_product( locU, matmul( trU%lo1d%m1d, locV ) )
          !J1d * sc_product( locU, trU%lo1d%m1d, locV )
       endif
    enddo

  end function dot1dsc

  !> dot product in 1D on a set of edges
  !> @param trU is supposed to be sing valued over edges
  !> @param trV
  function dot1dvec(trSig, trTau, edgeindex)
    type(trflux2d),    intent(in) :: trSig, trTau
    integer, optional, intent(in), dimension(:) :: edgeindex

    real :: dot1dvec

    !integer :: n
    integer :: dummy
    integer :: dim, n_edgeidx
    integer, dimension(:), allocatable :: edgeidx


    if ( trsig%cp(1)%Th%n_dim .ne. trtau%cp(1)%Th%n_dim ) then
       write(6,*) "error in dot1dvec"
       stop
    end if

    if ( present(edgeindex) ) then
       n_edgeidx = size( edgeindex )
       allocate( edgeidx(n_edgeidx) )
       edgeidx   = edgeindex
    else
       n_edgeidx = trSig%cp(1)%Th%n_edge
       allocate( edgeidx(n_edgeidx) )
       do dummy=1, n_edgeidx
          edgeidx(dummy) = dummy
       enddo
       !edgeidx = (/ (n, n=1, n_edgeidx) /)
    end if

    dot1dvec = 0.0

    do dim=1, trSig%cp(1)%Th%n_dim
       dot1dvec = dot1dvec + dot1dsc( trSig%cp(dim), trTau%cp(dim), edgeidx )
    enddo

  end function dot1dvec

  !> product of single valued constant with tr2d object
  function paramdot1dsc(param, trU) result(out)
    type(paramtr2d), intent(in) :: param
    type(tr2d)     , intent(in) :: trU

    type(tr2d) :: out

    integer :: globf

    if ( trU%init .eqv..false.  .or.  &
         param%init .eqv..false.      &
         ) then
       write(6,*) "error in paramdot1dsc"
       stop
    endif

    !call buildtr2d(trU%Th, trU%lo1d, trU%k, out)
    call copytr2tr( trU, out )

    out%face2idof = trU%face2idof

    do globf=1, trU%Th%n_faces
       out%face2dof(globf,:) = param%face2param(globf) * trU%face2dof(globf,:)
    enddo

  end function paramdot1dsc

  function paramdot1dvec(param, trSig) result(out)
    type(paramtr2d), intent(in) :: param
    type(trflux2d) , intent(in) :: trSig

    type(trflux2d) :: out

    integer :: dim 
    if ( trSig%cp(1)%init .eqv..false.  .or.  &
         trSig%link .eqv. .false.       .or.  &
         param%init .eqv..false.      &
         ) then
       write(6,*) "error in paramdot1dvec", trSig%cp(1)%init, trSig%link, param%init
       stop
    endif

    !call buildtrflux2d(trSig%cp(1)%Th, trSig%cp(1)%lo1d, trSig%cp(1)%k, out)

    call copytr2trflux( trSig%cp(1), out )

    do dim=1, trSig%cp(1)%Th%n_dim
       out%cp(dim) = paramdot1dsc(param, trSig%cp(dim))
    enddo

  end function paramdot1dvec

  !> multiply a constant to a trace object
  function consdot1dsc(cons, trU) result(out)
    real,       intent(in) :: cons
    type(tr2d), intent(in) :: trU

    type(tr2d) :: out

    if ( trU%init .eqv..false.  &
         ) then
       write(6,*) "error in constdot1dsc"
       stop
    endif

    call copytr2tr( trU, out )

    out%face2dof = trU%face2dof * cons

  end function consdot1dsc

  !> multiply a constant to a trace object
  function consdot1dvec(cons, trSig) result(out)
    real,       intent(in) :: cons
    type(trflux2d) , intent(in) :: trSig

    type(trflux2d) :: out

    integer :: dim 
    if ( trSig%cp(1)%init .eqv..false.  .or.  &
         trSig%link .eqv. .false.       &
         ) then
       write(6,*) "error in consdot1dvec"
       stop
    endif

    call copytr2trflux( trSig%cp(1), out )

    do dim=1, trSig%cp(1)%Th%n_dim
       out%cp(dim)%face2dof = cons * trSig%cp(dim)%face2dof
    enddo
  end function consdot1dvec

  !> add two dg trace objects 
  function addtrace(trU,trV) result(out)
    type(tr2d), intent(in) :: trU, trV

    type(tr2d) :: out

    if ( trU%init .eqv..false.  .or. &
         trV%init .eqv..false.       &
         ) then
       write(6,*) "error in addtrace"
       stop
    endif

    call copytr2tr(trU, out)

    out%face2dof = trU%face2dof + trV%face2dof
  end function addtrace

  !> add two dg flux traces
  function addfluxtrace(trSig, trTau) result(out)
    type(trflux2d) , intent(in) :: trSig, trTau

    type(trflux2d) :: out

    integer :: dim
    if ( trSig%init .eqv..false.  .or. &
         trSig%init .eqv..false.       &
         ) then
       write(6,*) "error in addfluxtrace"
       stop
    endif

    call copytr2trflux(trSig%cp(1), out)

    do dim=1, trSig%cp(1)%Th%n_dim
       out%cp(dim)%face2dof = trSig%cp(dim)%face2dof + trTau%cp(dim)%face2dof
    end do

  end function addfluxtrace

  subroutine plotdg2d(filename, U)
    character(len=*), intent(in) :: filename
    type(dg2d),       intent(in) :: U

    integer :: dof, myunit

    if ( U%init.eqv..false. ) then
       write(6,*) "error in plotdg2d"
       stop
    endif

    myunit = 11
    open(myunit, file=filename//'.gnu')

    do dof=1, U%n_dof
       write(myunit,*) U%dof2corr(dof,:), U%dof(dof)
    enddo

    close(myunit)

  end subroutine plotdg2d


end module dgopt2d
