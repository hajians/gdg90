!> nestddmesh2d.f90

module nestddmesh2d

  use triangulate2d
  use pnpoly

  implicit none

  type subdomain2el
     integer :: Nele_subd
     integer, allocatable :: ele(:)
  end type subdomain2el
  
  type nestmesh2d
     logical :: init = .false.

     type(tri2d), pointer :: ThCoarse => null()
     type(tri2d), pointer :: Thfine   => null()

     integer :: Nsubd
     integer :: N_edgeIntF

     logical, allocatable, dimension(:,:) :: el2subd
     integer, allocatable, dimension(:)   :: el2subidx

     type(subdomain2el), allocatable, dimension(:) :: subd2el

     logical, allocatable, dimension(:) :: edgeonIntF
     integer, allocatable, dimension(:) :: edgeidx

     logical, allocatable, dimension(:) :: edgenotIntF
     integer, allocatable, dimension(:) :: edgenotidx     

     integer, allocatable, dimension(:)   :: sub2NintFedge

     !> edge2subidx maps the interface edge index to the subdomain index
     integer, allocatable, dimension(:,:) :: edge2subidx


  end type nestmesh2d
contains

  !> 
  subroutine buildnestmesh2d(ThC, ThF, nest2d)
    type(tri2d), target, intent(in) :: ThC
    type(tri2d), target, intent(in) :: ThF

    type(nestmesh2d), intent(inout) :: nest2d

    integer :: dim, n_vtx, el, subd, subdel
    integer :: edge, edgeF, edgeC

    real, allocatable, dimension(:,:) :: ptsF, ptsC
    real, allocatable, dimension(:,:) :: vtx

    logical :: inside

    logical :: isin1, isin2

    if ( nest2d%init .eqv. .true. ) then
       write(6,*) "error in buildnestmesh2d: nest2d is initialized", nest2d%init
    endif

    nest2d%ThCoarse => ThC
    nest2d%Thfine   => ThF

    nest2d%Nsubd = ThC%n_el
    allocate( nest2d%el2subd( ThF%n_el, nest2d%Nsubd ) )
    nest2d%el2subd = .false.

    allocate( nest2d%el2subidx( Thf%n_el ) )
    nest2d%el2subidx = 0

    dim   = 2
    n_vtx =  size(ThC%el2v, 2) 

    allocate( vtx(dim,n_vtx) )

!!!
    do el=1, ThF%n_el
       inside = .false.
       subd   = 0
       !
       do while ( inside .eqv. .false. )      
          subd = subd + 1   
          vtx = transpose( ThC%v2corr(ThC%el2v(subd,:),:) )
          call polygon_contains_point_2d( n_vtx, vtx, ThF%centroid(el,:), inside )     
       end do
       !
       if ( inside .eqv. .true. ) then
          nest2d%el2subd(el,subd) = .true.
          nest2d%el2subidx(el) = subd
       endif
    end do
!!!

    ! build subdomain to elements
    allocate( nest2d%subd2el(ThC%n_el) )

    do subd=1, nest2d%Nsubd
       subdel = 0
       do el=1, ThF%n_el
          if ( nest2d%el2subd(el,subd).eqv..true. ) then
             subdel = subdel + 1
          endif
       enddo
       nest2d%subd2el(subd)%Nele_subd = subdel
       allocate( nest2d%subd2el(subd)%ele( subdel ) )
    enddo
    
    !! it is insured that the elements are sorted from smallest to
    !! bigest.
    do subd=1, nest2d%Nsubd
       subdel = 0
       do el=1, ThF%n_el
          if ( nest2d%el2subd(el,subd).eqv..true. ) then
             subdel = subdel + 1
             nest2d%subd2el(subd)%ele(subdel) = el
          endif
       enddo

    enddo
    !

    ! allocating the number of edges of a given subdomain
    allocate( nest2d%sub2NintFedge(nest2d%Nsubd) )
    nest2d%sub2NintFedge = 0

    allocate( nest2d%edgeonIntF( ThF%n_edge ) )
    nest2d%edgeonIntF = .false.

    ! finding edges on the interface
    nest2d%N_edgeIntF = 0

    allocate( ptsC( size(ThC%edge2v,2), dim ) )
    allocate( ptsF( size(ThF%edge2v,2), dim ) )

    ! finding interface edges of the ThF
    do edgeF=1, ThF%n_edge
       inside = .false.
       edgeC  = 0
       ptsF = ThF%v2corr( ThF%edge2v(edgeF,:), :)
       !
       do while ( inside .eqv. .false. )
          edgeC = edgeC + 1
          ptsC  = ThC%v2corr( ThC%edge2v(edgeC,:), :)

          call point_contains_line_2d(ptsC,ptsF(1,:),isin1)
          call point_contains_line_2d(ptsC,ptsF(2,:),isin2)

          if ( isin1 .and. isin2 .and. (ThF%edge2bd(edgeF).eq.-1) ) then ! found
             ! one
             inside = .true.
             nest2d%N_edgeIntF = nest2d%N_edgeIntF + 1
             nest2d%edgeonIntF(edgeF) = .true.

          endif

          if ( edgeC .eq. ThC%n_edge ) then ! if we havent found then it is not interface
             inside = .true.
          endif

       end do
       !
    end do

    allocate( nest2d%edge2subidx( nest2d%N_edgeIntF, 2 ) ) ! 2 subdomain per edge
    nest2d%edge2subidx = 0

    ! filling edgeIdx
    allocate( nest2d%edgeIdx( nest2d%N_edgeIntF ) )
    edge = 0

    do edgeF=1, ThF%n_edge
       if ( nest2d%edgeonIntF(edgeF) .eqv. .true. ) then
          edge = edge + 1
          !
          nest2d%edgeIdx(edge) = edgeF

          ! finding the number of edges per subdomain
          el = ThF%edge2el(edgeF,1)
          subd = nest2d%el2subidx(el)
          nest2d%sub2NintFedge(subd) = nest2d%sub2NintFedge(subd)+1
          nest2d%edge2subidx(edge,1) = subd

          ! finding the neighbor subdomain
          el = ThF%edge2el(edgeF,2)
          if (el.ne.-1) then
             subd = nest2d%el2subidx(el)
             nest2d%sub2NintFedge(subd) = nest2d%sub2NintFedge(subd)+1
             nest2d%edge2subidx(edge,2) = subd
          else
             nest2d%edge2subidx(edge,2) = -1
          endif
          !
       endif
    enddo

    ! edges which are not on the interface
    allocate( nest2d%edgenotIntF(ThF%n_edge) )
    nest2d%edgenotIntF = .true.

    allocate( nest2d%edgenotidx( ThF%n_edge - nest2d%N_edgeIntF ) )
    nest2d%edgenotidx = 0

    edge = 0
    do edgeF=1, ThF%n_edge
       if ( nest2d%edgeonIntF(edgeF).eqv..true. ) then
          nest2d%edgenotIntF(edgeF) = .false.
       else
          edge = edge + 1
          nest2d%edgenotidx(edge) = edgeF
       endif
    enddo

    nest2d%init = .true.

  end subroutine buildnestmesh2d

  subroutine plotnestinterface(filename,nest2d)
    character(len=*), intent(in) :: filename
    type(nestmesh2d), intent(in) :: nest2d

    integer :: myunit = 10
    integer :: n_pts
    integer :: edge, pts

    if (nest2d%init .eqv. .false.) then
       write(6,*) "error in plotnestinterface: nest2d is not init."
       stop
    endif
    !
    open(myunit, file=filename//'.edge')

    n_pts = 2                     ! two points per edge

    do edge=1, nest2d%N_edgeIntF
       do pts=1, n_pts
          write(myunit,*) &
               nest2d%ThFine%v2corr( nest2d%Thfine%edge2v(nest2d%edgeIdx(edge),pts), : )
       enddo
       write(myunit,*) 
    enddo
    close(myunit)

  end subroutine plotnestinterface


end module nestddmesh2d
