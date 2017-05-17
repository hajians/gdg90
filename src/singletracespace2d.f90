  !> singletracespace2d.f90

module singletracespace2d

  use triangulate2d
  use localopt1d

  implicit none

  type singletr2d
     logical :: init = .false.
     logical :: link = .false.

     type(tri2d),    pointer :: Th   => null()
     type(locopt1d), pointer :: lo1d => null()

     integer :: NpF
     integer :: k
     integer :: NintfEdge
     integer :: Nedgedof

     integer, allocatable, dimension(:) :: intfEdgeIdx
     real,    allocatable, dimension(:) :: edge2J1d

     integer, allocatable, dimension(:,:) :: edge2idof
     real,    allocatable, dimension(:)   :: edgedof
     integer, allocatable, dimension(:,:) :: edgedof2facedof
     
     logical :: binary = .false.
     integer, allocatable :: edgenotzero(:)

  end type singletr2d

  type singletrflux2d
     logical :: init = .false.
     logical :: link = .false.

     type(singletr2d), allocatable, dimension(:) :: cp
  end type singletrflux2d

contains

  subroutine buildsingletr2d(Th,lo1d,k,edgeidx,sltr2d)
    type(tri2d),    target, intent(in) :: Th
    type(locopt1d), target, intent(in) :: lo1d
    integer,                intent(in) :: k
    integer, dimension(:),  intent(in) :: edgeidx

    type(singletr2d), intent(inout) :: sltr2d

    integer :: i, j, edge
    real, dimension(2,2) :: edgevtx2coor

    if (sltr2d%init .eqv. .true. ) then
       write(6,*) "error in buildsingletr2d: sltr2d is initialized", sltr2d%init
       stop
    endif

    sltr2d%k    = k
    sltr2d%NpF  = k+1
    sltr2d%Th   => Th
    sltr2d%lo1d => lo1d

    sltr2d%NintfEdge = size( edgeidx )

    sltr2d%Nedgedof = sltr2d%NpF*sltr2d%NintfEdge
    
    allocate( sltr2d%intfEdgeIdx( sltr2d%NintfEdge ) )
    allocate( sltr2d%edge2J1d( sltr2d%NintfEdge ) )
    allocate( sltr2d%edge2idof( sltr2d%NintfEdge, sltr2d%NpF ) )
    allocate( sltr2d%edgedof( sltr2d%NintfEdge*sltr2d%NpF ) )
    allocate( sltr2d%edgenotzero( sltr2d%NintfEdge ) )

    sltr2d%intfEdgeIdx = edgeidx

    do i=1, sltr2d%NintfEdge
       edge = edgeidx(i)      
       edgevtx2coor          = Th%v2corr( Th%edge2v(edge,:), :)
       sltr2d%edge2J1D(i) = &
            0.5*sqrt( sum( (edgevtx2coor(1,:)-edgevtx2coor(2,:))**2 ) )
    enddo

    sltr2d%edge2idof = 0

    do i=1, sltr2d%NintfEdge
       sltr2d%edge2idof(i,:) = (/ (j, j=( (i-1)*sltr2d%NpF + 1), i*sltr2d%NpF ) /)
    enddo

    sltr2d%edgedof = 0.0
    
    sltr2d%edgenotzero = 0
    !

    allocate( sltr2d%edgedof2facedof( sltr2d%NpF * sltr2d%NintfEdge, 2 ) ) ! 2 faces
    sltr2d%edgedof2facedof = 0

    sltr2d%init = .true.
  end subroutine buildsingletr2d

  !> build single valued trace for flux
  subroutine buildsinglefluxtr2d(Th,lo1d,k,edgeidx,sltrflux2d)
    type(tri2d),    target, intent(in) :: Th
    type(locopt1d), target, intent(in) :: lo1d
    integer,                intent(in) :: k
    integer, dimension(:),  intent(in) :: edgeidx

    type(singletrflux2d), intent(inout) :: sltrflux2d

    integer :: dim

    if ( sltrflux2d%init .eqv. .true. ) then
       write(6,*) "error in buildsinglefluxtr2d"
       stop
    endif

    allocate( sltrflux2d%cp( Th%n_dim ) )
    
    do dim=1, Th%n_dim
       call buildsingletr2d(Th,lo1d,k,edgeidx,sltrflux2d%cp(dim) )
    enddo

    sltrflux2d%init = .true.
    
  end subroutine buildsinglefluxtr2d


end module singletracespace2d
