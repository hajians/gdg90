!> ddopt.f90 contains domain decomposition operators

module ddopt2d
  use triangulate2d
  use ddgeo2d
  use dgspace2d
  use pnpoly

  implicit none

  !> restriction object
  type restr
     real, allocatable :: mat(:,:)
  end type restr

  !> domain decomposition object
  type dd2d
     logical :: init = .false.

     type(tri2d),  pointer :: Th => null()
     type(geom2d), pointer :: geo2d => null()
     type(dg2d),   pointer :: U  => null()

     integer :: n_subd

     logical, allocatable :: el2sub(:,:)
     integer, allocatable :: subN_dof(:)
     integer, allocatable :: subN_el(:)
     
     logical, allocatable :: edge2IF(:)
     integer              :: n_IFedge
     integer, allocatable :: IFedge(:)
     
     type(restr), allocatable :: res(:)
  end type dd2d

  interface Trp
     module procedure transposeRest
  end interface Trp

  interface operator (.dot.)
     module procedure resdotvec, resdotmat, matdotres
  end interface operator (.dot.)

  ! private
  private resdotvec
  private resdotmat
  private matdotres
  private transposeRest

contains

  !> determines whether or not a point is inside a convex hull
  subroutine builddd2d(Th, geo2d, U, ddopt)
    type(tri2d),  target, intent(in) :: Th
    type(geom2d), target, intent(in) :: geo2d
    type(dg2d),   target, intent(in) :: U
    type(dd2d),  intent(inout) :: ddopt

    integer :: sub, el, adjel, locf, face, vtx, edge
    integer :: elsub
    integer :: dummy

    !integer :: adjelsub

    logical :: inside
    logical, allocatable :: elflag(:), edgeflag(:)
    !real, allocatable    :: centroid(:,:)

    if ( Th%init.eqv..false.    .or. &
         geo2d%init.eqv..false. .or. &
         U%init.eqv..false.     .or. &
         ddopt%init.eqv..true. ) then
       write(6,*) "error in builddd2d", Th%init, U%init, ddopt%init
       stop
    endif

    ddopt%Th => Th
    ddopt%geo2d => geo2d
    ddopt%U  => U

    ddopt%n_subd = geo2d%n_subd

    allocate( ddopt%el2sub( Th%n_el, ddopt%n_subd ) )
    ddopt%el2sub = .false.
    ! allocate( centroid( Th%n_el, Th%n_dim ) )
    ! centroid = 0.0
    allocate( elflag( Th%n_el ), edgeflag( Th%n_edge ) )
    elflag   = .false.
    edgeflag = .false.

    vtx  = size( Th%el2v, 2 )    ! three vertex per element in Th
    face = size( Th%el2face, 2 ) ! number of faces of an element

    allocate( ddopt%res( ddopt%n_subd ) ) ! allocate restriction operators
    allocate( ddopt%subN_el( ddopt%n_subd ) ) ! # dof of each subdomain
    allocate( ddopt%subN_dof( ddopt%n_subd ) ) ! # dof of each subdomain
    allocate( ddopt%edge2IF( Th%n_edge ) )     ! edge to interface
    ddopt%subN_el  = 0
    ddopt%subN_dof = 0
    ddopt%edge2IF = .false.
    ddopt%n_IFedge = 0
!!!!! filling el2sub and subN_dof
    do sub=1, ddopt%n_subd
       do el=1, Th%n_el
          if ( elflag(el).eqv..false. ) then

             call polygon_contains_point_2d &
                  (geo2d%subd(sub)%n_vert, transpose(geo2d%subd(sub)%v2corr),&
                  Th%centroid(el,:),inside)
             !
             if ( inside.eqv..true. ) then ! if it is inside
                ddopt%el2sub(el,sub) = .true.
                ddopt%subN_el(sub) = ddopt%subN_el(sub) + 1
                elflag(el) = .true.
             endif
             !
          endif
       enddo
       ddopt%subN_dof(sub) = ddopt%subN_el(sub)*U%Np
    enddo
!!!!!
    do el=1, Th%n_el
       !
       do sub=1, ddopt%n_subd
          if ( ddopt%el2sub(el,sub).eqv..true.) then
             elsub = sub
          endif
       enddo
       !
       do locf=1, face
          adjel = Th%el2el(el,locf)
          if (adjel.ne.-1) then                             ! not on the the boundary
             !
             if ( ddopt%el2sub(adjel,elsub).eqv..true. ) then ! if in same subdomain
             else              ! if not in the same subdomain
!!!!!!!!!! works for only two subdomain !!!!!!!!!!
                if ( edgeflag( Th%el2edge(el,locf) ) .eqv. .false. ) then ! if not visited
                   edgeflag( Th%el2edge(el,locf) ) = .true.
                   ddopt%n_IFedge = ddopt%n_IFedge + 1      ! we found an interface edge
                   ddopt%edge2IF( Th%el2edge(el,locf) ) = .true.
                endif
                !
             endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          endif
       enddo
       !
    enddo

    allocate( ddopt%IFedge( ddopt%n_IFedge ) )
    ddopt%IFedge = 0
    dummy = 0
    do edge=1, Th%n_edge
       if ( ddopt%edge2IF(edge).eqv..true. ) then
          dummy = dummy + 1
          ddopt%IFedge(dummy) = edge
       endif
    enddo
    ! init the restriction matrix of each subdomain
    call buildrestrictiondg2d(U, ddopt)

    ddopt%init = .true.
  end subroutine builddd2d

  !> build restriction operators
  subroutine buildrestrictiondg2d(U, ddopt)
    type(dg2d),   target, intent(in) :: U
    type(dd2d),  intent(inout) :: ddopt

    integer :: el, sub
    integer :: subN_dof
    integer :: i, locel
    integer, allocatable :: locdof(:), globdof(:)

    real, allocatable :: ident(:,:)

    allocate( locdof(U%Np), globdof(U%Np), ident(U%Np,U%Np) )
    ident = 0.0
    do i=1, U%Np
       ident(i,i) = 1.0
    enddo


    do sub=1, ddopt%n_subd
       subN_dof = U%Np * ddopt%subN_el(sub)
       !
       allocate( ddopt%res(sub)%mat(subN_dof, U%n_dof) )
       ddopt%res(sub)%mat = 0.0

       locel = 0
       !
       do el=1, U%Th%n_el
          if ( ddopt%el2sub(el,sub) .eqv. .true. ) then
             locel = locel + 1

             globdof = U%el2dof(el,:)
             locdof  = (/ (i, i=(locel-1)*U%Np+1, (locel*U%Np) ) /)
             !
             ddopt%res(sub)%mat( locdof , globdof ) = ident
             !
          endif
       enddo
       !
    enddo
  end subroutine buildrestrictiondg2d

  !!
  !> operations for restrictions
  function transposeRest(res) result(resT)
    type(restr), intent(in) :: res
    type(restr)             :: resT

    allocate( resT%mat(size(res%mat,2), size(res%mat,1)) )

    resT%mat = transpose( res%mat )
  end function transposeRest

  !>
  function resdotvec(res,vec) result(out)
    type(restr), intent(in) :: res
    real,        intent(in) :: vec(:)

    real, dimension( size(res%mat,1) ) :: out

    out = matmul( res%mat, vec )
  end function resdotvec

  !>
  function resdotmat(res,mat) result(out)
    type(restr), intent(in) :: res
    real,        intent(in) :: mat(:,:)

    real, dimension( size(res%mat,1), size(mat,2) ) :: out

    out = matmul( res%mat, mat )
  end function resdotmat

  !>
  function matdotres(mat,res) result(out)
    type(restr), intent(in) :: res
    real,        intent(in) :: mat(:,:)

    real, dimension( size(mat,1), size(res%mat,2) ) :: out

    out = matmul( mat, res%mat )
  end function matdotres


end module ddopt2d
