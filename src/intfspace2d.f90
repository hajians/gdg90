!> intfspace2d.f90 contains space of functions along interface which are
!! double valued.

module intfspace2d

  use nestddmesh2d
  use alldgspace2d
  
  implicit none

  !> intfsp2d type is the space of double valued functions along
  !! the interface between two subdomain.
  type intfsp2d

     logical :: init = .false.

     type(nestmesh2d), pointer :: nest2d => null()

     integer :: NintfEdge
     
     type(singletr2d), allocatable :: sltr(:)

     type(singletr2d) :: Mu
     type(singletr2d) :: MuMinusOne
  end type intfsp2d

contains

  subroutine buildinterface2d(ThF,V,nest2d,interface2d)
    type(tri2d),              intent(in) :: ThF
    type(primaldg2d),         intent(in) :: V
    type(nestmesh2d), target, intent(in) :: nest2d

    type(intfsp2d), intent(inout) :: interface2d

    ! number of adj subdomains of an interface edge
    integer, parameter :: Nside = 2

    !> interface edge
    !integer :: intf_edge

    !> fine edge index
    !integer :: fine_edge_idx

    integer :: side
    integer :: edge, globedge
    integer, allocatable :: locEdgedof(:), locdof(:)
    
    if ( interface2d%init .eqv. .true. ) then
       write(6,*) "error in interface2d"
    endif

    interface2d%nest2d => nest2d

    allocate( interface2d%sltr(Nside) )

    ! building the interface spaces
    do side=1, Nside
       call buildsingletr2d(ThF,V%sc%lo1d, V%sc%k, nest2d%edgeIdx, interface2d%sltr(side) )
    enddo
    ! build the corresponding space for penalty
    call buildsingletr2d(ThF,V%sc%lo1d,V%sc%k,nest2d%edgeIdx,interface2d%Mu)
    call buildsingletr2d(ThF,V%sc%lo1d,V%sc%k,nest2d%edgeIdx,interface2d%MuMinusOne)
    
    allocate( locEdgedof( V%trSc%NpF ), locdof( V%trSc%NpF ) )
    
    do edge=1, nest2d%N_edgeIntF
       !
       locEdgedof = interface2d%mu%edge2idof(edge,:)
       globedge = nest2d%edgeIdx(edge)
       locdof = V%paramSc%edge2idof(globedge,:)
       !
       interface2d%Mu%edgedof(locEdgedof) = V%paramSc%edgedof(locdof)
    enddo

    interface2d%MuMinusOne%edgedof = 1.0 / interface2d%Mu%edgedof
    
    interface2d%NintfEdge = nest2d%N_edgeIntF

    interface2d%init = .true.
  end subroutine buildinterface2d

  !> link interface single valued traces to the adjacent faces
  subroutine linkbrokentrace2interface(trV, interface2d)
    type(tr2d),     intent(in)    :: trV
    type(intfsp2d), intent(inout) :: interface2d

    integer, parameter :: Nadj_subd = 2   
    integer :: subd

    do subd=1, Nadj_subd
       call linktr(trV, interface2d%sltr(subd) )
    enddo

    call linktr(trV, interface2d%Mu )
    call linktr(trV, interface2d%MuMinusOne )
  end subroutine linkbrokentrace2interface

end module intfspace2d
