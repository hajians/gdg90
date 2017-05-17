  !> singletraceopt2d.f90

module singletraceopt2d
  use singletracespace2d
  use brokentracespace2d

  implicit none

  !> interface for trace links
  interface linktr
     module procedure linktr2d, linkfluxtr2d
  end interface linktr

  private linktr2d, linkfluxtr2d

  !> interface for averages
  interface avetr
     module procedure avetr2d, avefluxtr2d
  end interface avetr

  private avetr2d, avefluxtr2d

  !> interface for jump
  interface jumptr
     module procedure jumptr2d, jumpfluxtr2d
  end interface jumptr

  private jumptr2d, jumpfluxtr2d

  !> interface for dot_product
  interface dotproduct
     module procedure dotproductsc, dotproductflux
  end interface dotproduct

  private dotproductsc, dotproductflux
  
contains

  !> linking trace space to single valued dge spaces
  subroutine linktr2d(trV, SltrV)
    type(tr2d),       intent(in)    :: trV
    type(singletr2d), intent(inout) :: SltrV

    integer :: edge, globf, adjglobf, globEdge
    integer, allocatable, dimension(:) :: locedgedof

    if ( (trV%init .eqv. .false.)  .or. (trV%link .eqv. .false.) ) then
       write(6,*) "error in linktr2d: trV"
       stop
    elseif ( (SltrV%init .eqv. .false.) .or. (SltrV%link .eqv. .true.) ) then
       write(6,*) "error in linktr2d: SltrV"
       stop
    endif

    allocate( locedgedof( SltrV%NpF ) )

    do edge=1, SltrV%NintfEdge
       globEdge = SltrV%IntfEdgeIdx(edge)

       globf    = trV%Th%edge2face(globEdge,1) ! there are two faces
       adjglobf = trV%Th%edge2face(globEdge,2)

       locedgedof = SltrV%edge2idof(edge,:)

       if ( adjglobf .eq. -1  ) then
          SltrV%edgedof2facedof( locedgedof, 1 ) = &
               trV%face2idof(globf,:)
          SltrV%edgedof2facedof( locedgedof, 2 ) = -1
       else
          SltrV%edgedof2facedof( locedgedof, 1 ) = &
               trV%face2idof(globf,:)
          SltrV%edgedof2facedof( locedgedof, 2 ) = &
               trV%face2idof(adjglobf,:)
       endif

    enddo

    SltrV%link = .true.
  end subroutine linktr2d

  !> link trSig to Single trSig
  subroutine linkfluxtr2d(trSig,SltrSig)
    type(trflux2d),       intent(in)    :: trSig
    type(singletrflux2d), intent(inout) :: SltrSig

    integer :: dim

    if ( (trSig%init .eqv. .false.)  .or. (trSig%link .eqv. .false.) ) then
       write(6,*) "error in linkfluxtr2d: trSig"
       stop
    elseif ( (SltrSig%init .eqv. .false.) .or. (SltrSig%link .eqv. .true.) ) then
       write(6,*) "error in linkfluxtr2d: SltrSig"
       stop
    endif

    do dim=1, trSig%cp(1)%Th%n_dim
       call linktr2d(trSig%cp(dim), SltrSig%cp(dim))
    enddo

    SltrSig%link = .true.    
  end subroutine linkfluxtr2d

  !> lifts trace of trV on an edge from one of the two side
  subroutine tracelift2d(trV,whichface,SltrV)
    type(tr2d),       intent(in)    :: trV
    integer,          intent(in)    :: whichface
    type(singletr2d), intent(inout) :: SltrV

    integer :: edge, globf, globEdge
    integer, allocatable, dimension(:) :: locedgedof

    if ( (trV%init .eqv. .false.)  .or. (trV%link .eqv. .false.) ) then
       write(6,*) "error in avetr2d: trV"
       stop
    elseif ( (SltrV%init .eqv. .false.) .or. (SltrV%link .eqv. .false.) ) then
       write(6,*) "error in avetr2d: SltrV"
       stop
    endif

    allocate( locedgedof( SltrV%NpF ) )

    SltrV%edgedof = 0.0

    do edge=1, SltrV%NintfEdge
       globEdge = SltrV%IntfEdgeIdx(edge)

       locedgedof = SltrV%edge2idof(edge,:)

       globf = SltrV%Th%edge2face(globEdge,whichface)

       if (globf.ne.-1) then
          SltrV%edgedof(locedgedof) = trV%face2dof(globf,:)
       endif
    enddo
    
  end subroutine tracelift2d
  
  !> take average
  subroutine avetr2d(trV, SltrV)
    type(tr2d),       intent(in)    :: trV
    type(singletr2d), intent(inout) :: SltrV

    integer :: edge, face, globf, adjglobf, globEdge

    integer, allocatable, dimension(:) :: locedgedof

    if ( (trV%init .eqv. .false.)  .or. (trV%link .eqv. .false.) ) then
       write(6,*) "error in avetr2d: trV"
       stop
    elseif ( (SltrV%init .eqv. .false.) .or. (SltrV%link .eqv. .false.) ) then
       write(6,*) "error in avetr2d: SltrV"
       stop
    endif

    allocate( locedgedof( SltrV%NpF ) )

    do edge=1, SltrV%NintfEdge
       globEdge = SltrV%IntfEdgeIdx(edge)

       locedgedof = SltrV%edge2idof(edge,:)

       face  = 1
       globf = SltrV%Th%edge2face(globEdge,face)
       face  = 2
       adjglobf = SltrV%Th%edge2face(globEdge,face)       
       !
       ! take average
       if ( adjglobf.eq.-1 ) then
          SltrV%edgedof(locedgedof) = trV%face2dof(globf,:) 
       else
          SltrV%edgedof(locedgedof) = &
               0.5 * ( trV%face2dof(globf,:)    + &
               trV%face2dof(adjglobf,:)   &
               )
       endif
    enddo

  end subroutine avetr2d

  !> average for flux trace space
  subroutine avefluxtr2d(trSig, SltrSig)
    type(trflux2d),       intent(in)    :: trSig
    type(singletrflux2d), intent(inout) :: SltrSig

    integer :: dim

    if ( (trSig%init .eqv. .false.)  .or. (trSig%link .eqv. .false.) ) then
       write(6,*) "error in avefluxtr2d: trSig"
       stop
    elseif ( (SltrSig%init .eqv. .false.) .or. (SltrSig%link .eqv. .false.) ) then
       write(6,*) "error in avefluxtr2d: SltrSig"
       stop
    endif

    do dim=1, trSig%cp(1)%Th%n_dim
       call avetr2d( trSig%cp(dim), SltrSig%cp(dim) )
    end do
    !
  end subroutine avefluxtr2d

  !> jump of a scalar trace space
  subroutine jumptr2d(trV, SltrSig)
    type(tr2d),           intent(in)    :: trV
    type(singletrflux2d), intent(inout) :: SltrSig

    integer :: edge, globf, adjglobf, globEdge
    integer :: dim
    integer, allocatable, dimension(:) :: locedgedof

    !

    if ( (trV%init .eqv. .false.)  .or. (trV%link .eqv. .false.) ) then
       write(6,*) "error in jumptr2d: trSig"
       stop
    elseif ( (SltrSig%init .eqv. .false.) .or. (SltrSig%link .eqv. .false.) ) then
       write(6,*) "error in jumptr2d: SltrSig"
       stop
    endif

    allocate( locedgedof( SltrSig%cp(1)%NpF ) )

    do edge=1, SltrSig%cp(1)%NintfEdge
       globEdge = SltrSig%cp(1)%IntfEdgeIdx(edge)

       globf    = trV%Th%edge2face(globEdge,1) ! there are two faces
       adjglobf = trV%Th%edge2face(globEdge,2)

       ! x-component
       dim = 1
       locedgedof = SltrSig%cp(dim)%edge2idof(edge,:)

       if ( adjglobf.eq.-1 ) then
          SltrSig%cp(dim)%edgedof(locedgedof) = &
               trV%nx(globf)*trV%face2dof(globf,:) 
       else
          SltrSig%cp(dim)%edgedof(locedgedof) = &
               trV%nx(globf)*trV%face2dof(globf,:) + &
               trV%nx(adjglobf)*trV%face2dof(adjglobf,:)
       endif
       ! y-component
       dim = 2
       locedgedof = SltrSig%cp(dim)%edge2idof(edge,:)

       if ( adjglobf.eq.-1 ) then
          SltrSig%cp(dim)%edgedof(locedgedof) = &
               trV%ny(globf)*trV%face2dof(globf,:) 
       else
          SltrSig%cp(dim)%edgedof(locedgedof) = &
               trV%ny(globf)*trV%face2dof(globf,:) + &
               trV%ny(adjglobf)*trV%face2dof(adjglobf,:) 
       endif
       !
    enddo

  end subroutine jumptr2d

  !> build jump of trace flux
  subroutine jumpfluxtr2d(trSig, SltrV)
    type(trflux2d),   intent(in)    :: trSig
    type(singletr2d), intent(inout) :: SltrV

    integer :: edge, globf, adjglobf, globEdge

    integer, allocatable, dimension(:) :: locedgedof    

    if ( (trSig%init .eqv. .false.)  .or. (trSig%link .eqv. .false.) ) then
       write(6,*) "error in jumpfluxtr2d: trSig"
       stop
    elseif ( (SltrV%init .eqv. .false.) .or. (SltrV%link .eqv. .false.) ) then
       write(6,*) "error in jumpfluxtr2d: SltrV"
       stop
    endif

    allocate( locedgedof( SltrV%NpF ) )
    !
    do edge=1, SltrV%NintfEdge
       globEdge = SltrV%IntfEdgeIdx(edge)

       globf    = trSig%cp(1)%Th%edge2face(globEdge,1) ! there are two faces
       adjglobf = trSig%cp(1)%Th%edge2face(globEdge,2)

       locedgedof = SltrV%edge2idof(edge,:)

       if ( adjglobf.eq.-1) then
          SltrV%edgedof(locedgedof) = &
               trSig%cp(1)%nx(globf)*trSig%cp(1)%face2dof(globf,:) + &
               trSig%cp(2)%ny(globf)*trSig%cp(2)%face2dof(globf,:)       
       else
          SltrV%edgedof(locedgedof) = &
               trSig%cp(1)%nx(globf)*trSig%cp(1)%face2dof(globf,:)       + &
               trSig%cp(2)%ny(globf)*trSig%cp(2)%face2dof(globf,:)       + &
               trSig%cp(1)%nx(adjglobf)*trSig%cp(1)%face2dof(adjglobf,:) + &
               trSig%cp(2)%ny(adjglobf)*trSig%cp(2)%face2dof(adjglobf,:)                
       endif

    enddo
  end subroutine jumpfluxtr2d

  !> compute the dot product
  function dotproductsc(SltrU,SltrV) result(out)
    type(singletr2d), intent(in) :: SltrU, SltrV

    real :: out

    integer :: edge
    integer, allocatable, dimension(:) :: locedgedof

    out = 0.0
    allocate( locedgedof( SltrV%NpF ) )

    do edge=1, SltrV%NintfEdge
       locedgedof = SltrV%edge2idof(edge,:)
       out = out + SltrV%edge2J1D(edge) * &
            dot_product( SltrU%edgedof(locedgedof), &
            matmul( SltrV%lo1d%m1d, SltrV%edgedof(locedgedof) ) &
            )
    enddo

  end function dotproductsc

  !> compute the dot product of flux
  function dotproductflux(SltrSig,SltrTau) result(out)
    type(singletrflux2d), intent(in) :: SltrSig, SltrTau

    real :: out

    integer :: dim

    out = 0.0

    do dim=1, size( SltrSig%cp )
       out = out + dotproductsc( SltrSig%cp(dim), SltrTau%cp(dim) )
    enddo
  end function dotproductflux

  !> multiply two single valued trace function
  subroutine multiplysltr(SltrU, SltrV, SltrUV)
    type(singletr2d), intent(in) :: SltrU, SltrV

    type(singletr2d), intent(inout) :: SltrUV

    if ( (SltrU%init.eqv..false. .or. SltrU%link.eqv..false.) .and. &
         (SltrV%init.eqv..false. .or. SltrV%link.eqv..false.) .and. &
         (SltrUV%init.eqv..false. .or. SltrUV%link.eqv..false.) ) then
       write(6,*) "error in multiply"
       stop
    endif

    SltrUV%edgedof = SltrU%edgedof * SltrV%edgedof

  end subroutine multiplysltr

  !> multiply two single valued trace function for flux
  subroutine multiplySltrFlux(SltrSig, SltrTau, SltrSigTau)
    type(singletrflux2d), intent(in) :: SltrSig, SltrTau

    type(singletrflux2d), intent(inout) :: SltrSigTau

    integer :: dim

    do dim=1, SltrSig%cp(1)%Th%n_dim
       call multiplysltr(SltrSig%cp(dim), SltrTau%cp(dim), SltrSigTau%cp(dim) )
    enddo

  end subroutine multiplySltrFlux


end module singletraceopt2d
