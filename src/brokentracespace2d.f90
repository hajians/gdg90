  !> @file brokentracespace2d.f90

module brokentracespace2d
  use triangulate2d
  use localopt1d
  implicit none

  !> broken trace space in 2d
  type tr2d 
     logical :: init = .false.
     logical :: link = .false.  ! linking to a space definied over Th

     type(tri2d),    pointer :: Th => null()
     type(locopt1d), pointer :: lo1d => null()

     integer :: NpF
     integer :: k

     real, allocatable, dimension(:)      :: edge2J1d

     integer, allocatable, dimension(:,:) :: face2idof
     real, allocatable, dimension(:,:)    :: face2dof

     real, allocatable, dimension(:) :: nx
     real, allocatable, dimension(:) :: ny

     logical :: binary = .false.
     integer, allocatable :: facenotzero(:)
     integer :: elnotzero
  end type tr2d

  !> broken trace space for flux in 2d
  type trflux2d
     logical :: init = .false.
     logical :: link = .false.
     type(tr2d), allocatable, dimension(:) :: cp

  end type trflux2d


  type paramtr2d
     logical :: init = .false.

     type(tri2d), pointer :: Th => null()

     real, allocatable, dimension(:) :: face2param
  end type paramtr2d

contains

  subroutine buildtr2d(Th,lo1d,k,trV)
    type(tri2d),    target, intent(in) :: Th
    type(locopt1d), target, intent(in) :: lo1d
    integer,                intent(in) :: k

    type(tr2d), intent(inout) :: trV
    
    integer :: edge
    real, dimension(2,2) :: edgevtx2coor

    if (trV%init.eqv..true.) then
       write(6,*) "error in buildtrflux2d: trSig is initialized already"
       stop
    endif
    trV%k     = k
    trV%NpF   = k + 1
    trV%Th    => Th
    trV%lo1d  => lo1d

    allocate( trV%edge2J1D( Th%n_edge) )

    allocate( trV%face2dof( Th%n_faces, trV%NpF ) )
    allocate( trV%face2idof( Th%n_faces, trV%NpF ) )

    trV%face2dof  = 0.0
    trV%face2idof = 0

    do edge=1, Th%n_edge
       edgevtx2coor       = Th%v2corr( Th%edge2v(edge,:), :)
       trV%edge2J1D(edge) = &
            0.5*sqrt( sum( (edgevtx2coor(1,:)-edgevtx2coor(2,:))**2 ) )

    enddo

    allocate( trV%nx(Th%n_faces) )
    allocate( trV%ny(Th%n_faces) )

    allocate( trV%facenotzero( size(Th%el2face, 2)) ) ! three edge per element
    trV%facenotzero = 0
    trV%elnotzero   = 0
    trV%init = .true.
    
  end subroutine buildtr2d

  subroutine buildtrflux2d(Th,lo1d,k,trSig)
    type(tri2d),    target, intent(in) :: Th
    type(locopt1d), target, intent(in) :: lo1d
    integer,                intent(in) :: k

    type(trflux2d), intent(inout) :: trSig

    integer dim

    if (trSig%init.eqv..true.) then
       write(6,*) "error in buildtrflux2d: trSig is initialized already"
       stop
    endif

    allocate( trSig%cp( Th%n_dim ) )
    do dim=1, Th%n_dim
       call buildtr2d(Th,lo1d,k, trSig%cp(dim) )
    enddo

    trSig%init = .true.
  end subroutine buildtrflux2d


  subroutine buildparamtr2d(Th, param)
    type(tri2d), target, intent(in)    :: Th
    type(paramtr2d),     intent(inout) :: param
    
    if (param%init.eqv..true.) then
       write(6,*) "error in buildparamtr2d: param is initialized already"
       stop
    endif

    param%Th => Th

    allocate( param%face2param( Th%n_faces ) )

    param%face2param = 0.0

    param%init = .true.
  end subroutine buildparamtr2d

  subroutine copytr2tr(trU, trV)
    type(tr2d), intent(in)    :: trU
    type(tr2d), intent(inout) :: trV


    if ( trU%init .eqv. .false. .or. &
         trU%link .eqv. .false. .or. &
         trV%init .eqv. .true.  .or. &
         trV%link .eqv. .true.  ) then
       write(6,*) "error in copytr2tr"
       stop
    endif

    trV%Th => trU%Th
    trV%lo1d => trU%lo1d

    trV%NpF = trU%NpF
    trV%k   = trU%k

    allocate( trV%edge2J1D( trV%Th%n_edge ) )
    
    allocate( trV%face2idof( trV%Th%n_faces, trV%NpF ) )
    allocate( trV%face2dof( trV%Th%n_faces, trV%NpF ) )

    allocate( trV%nx( trV%Th%n_faces ) )
    allocate( trV%ny( trV%Th%n_faces ) )

    trV%edge2J1D = trU%edge2J1d

    trV%face2idof = trU%face2idof
    trV%face2dof  = trU%face2dof

    trV%nx = trU%nx
    trV%ny = trU%ny

    trV%binary = trU%binary
    allocate( trV%facenotzero(size(trU%facenotzero )) )
    trV%facenotzero = trU%facenotzero
    
    trV%init = .true.
    trV%link = .true.
  end subroutine copytr2tr

  !>
  !!
  subroutine copytr2trflux(trU, trSig)
    type(tr2d),     intent(in)    :: trU
    type(trflux2d), intent(inout) :: trSig

    integer :: dim 
    if ( trU%init .eqv. .false. .or. &
         trU%link .eqv. .false. .or. &
         trsig%init .eqv. .true.  .or. &
         trsig%link .eqv. .true.  ) then
       write(6,*) "error in copytr2trflux"
       stop
    endif

    allocate( trSig%cp( trU%Th%n_dim ) )

    do dim=1, trU%Th%n_dim
       call copytr2tr( trU, trSig%cp(dim) )
    enddo

    trSig%init = .true.
    trSig%link = .true.
    
  end subroutine copytr2trflux
  
    
  
end module brokentracespace2d

