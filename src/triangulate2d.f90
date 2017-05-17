!> @file triangulate2d.f90

module triangulate2d
  implicit none

  type tri2d

     logical :: init = .false.

     !> n_dim dimenstion of the domain
     integer :: n_dim

     !> n_el number of elements in triangulation
     integer :: n_el

     !> n_vtx number of vertices per element. It is contants for all.
     integer :: n_vtx

     
     !> n_faces total number of faces.
     integer :: n_faces

     !> n_egde total number of edges
     integer :: n_edge

     !> n_Bedge total number of boundary edges
     integer :: n_Bedge

     !> n_Iedge total number of interior edges
     integer :: n_Iedge


     
     !> el2v element to vertices index map
     integer, allocatable, dimension(:,:) :: el2v

     !> v2corr vertices index to coordinate map
     real,    allocatable, dimension(:,:) :: v2corr

     !> centroid centroid of an element
     real,    allocatable, dimension(:,:) :: centroid



     !> el2el adjacent element to element map
     integer, allocatable, dimension(:,:) :: el2el
     
     !> el2face element to global face index map
     integer, allocatable, dimension(:,:) :: el2face

     !> face2v global face index to vertices map
     integer, allocatable, dimension(:,:) :: face2v

     !> face2face global face to face map
     integer, allocatable, dimension(:) :: face2face

     !> face2locf global face to local face of an element map
     integer, allocatable, dimension(:) :: face2locf


     
     !> edge2el global edge index to element map
     integer, allocatable, dimension(:,:) :: edge2el

     !> el2edge global element to global edge index
     integer, allocatable, dimension(:,:) :: el2edge

     !> edge2face global edge to global face index
     integer, allocatable, dimension(:,:) :: edge2face
     
     !> edge2v global edge to global vertices index
     integer, allocatable, dimension(:,:) :: edge2v
     
     !> edge2bd shows if a global edge on the boundary
     integer, allocatable, dimension(:)   :: edge2bd
     
     !> inedge_idx global indices of the interior edges
     integer, allocatable, dimension(:)   :: inedge_idx
     
     !> bdedge_idx global indices of the boundary edges
     integer, allocatable, dimension(:)   :: bdedge_idx
     
  end type tri2d

  ! operator overloading for assignment
  ! interface assignment(=)
  !   module procedure copymesh
  ! end interface assignment(=)
  
  
  
  
  ! private subroutines
  private buildfacetri2d
  private buildel2eltri2d
  private buildedgestri2d
  private copymesh

contains

  !> Read a mesh in 2-dimenstion.
  !> @param[in] filename name of mesh file
  !> @param[inout] Th triangulation
  subroutine readtri2d(filename,Th)
    character(len=*), intent(in)    :: filename
    type(tri2d),      intent(inout) :: Th

    integer :: myunit = 10
    integer :: n_attr, vtx, bdd
    integer :: i, dummy

    !write(6,*) "start reading mesh"
    ! read .ele file
    open(myunit, file=filename//'.ele')
    read(myunit,*) Th%n_el, vtx, n_attr
    allocate( Th%el2v( Th%n_el, vtx ) ) ! initiate the element to vertices

    do i=1, Th%n_el
       read(myunit,*) dummy, Th%el2v(i,:)
    enddo
        
    close(myunit)

    ! read .node file
    open(myunit, file=filename//'.node')
    read(myunit,*) Th%n_vtx, Th%n_dim, n_attr, bdd
    allocate( Th%v2corr( Th%n_vtx, Th%n_dim ) )

    do i=1, Th%n_vtx
       read(myunit,*) dummy, Th%v2corr(i,:), dummy
    enddo
    
    close(myunit)
    !write(6,*) "end reading mesh"
  end subroutine readtri2d

  !> Plot a mesh in 2-dimension.
  !> @param[in] filename name of mesh file
  !> @param[inout] Th triangulation
  subroutine plottri2d(filename,Th)
    character(len=*), intent(in) :: filename
    type(tri2d),      intent(in) :: Th

    integer :: myunit = 10
    integer :: i, j, vtx
    character(len=80) :: str
    
    vtx = size(Th%el2v, 2)

    open(myunit, file=filename//'.gnu')

    do i=1, Th%n_el
       do j=1, vtx
          write(myunit,*) Th%v2corr(Th%el2v(i,j), :)
       enddo
       write(myunit,*) Th%v2corr(Th%el2v(i,1), :)
       write(myunit,*) 
    enddo

    close(myunit)

    vtx = size(Th%edge2v, 2)
    open(myunit, file=filename//'.edge')
    open(myunit+1, file=filename//'.Bedge')
    open(myunit+2, file=filename//'.ele')
    
    do i=1, Th%n_edge
       do j=1, vtx
          write(myunit,*) Th%v2corr( Th%edge2v(i,j), :)
          if ( Th%edge2bd(i).ne.-1 ) then
             write(myunit+1,*) Th%v2corr( Th%edge2v(i,j), :)
          end if
       end do
       !
       write(myunit,*)          ! make blank lines
       write(myunit+1,*)        ! make blank lines
    end do

    str = "('set label ',A1,i2,A1,' at first', f5.2, ', first', f5.2)"
    do i=1, Th%n_el
       write(myunit+2,str) "''",i,"''", Th%centroid(i,1), Th%centroid(i,2)
    enddo
    
    close(myunit)
    close(myunit+1)
    close(myunit+2)

  end subroutine plottri2d

  !> build faces' map
  subroutine buildfacetri2d(Th)
    type(tri2d), intent(inout) :: Th

    integer :: faces ! number of faces of an element
    integer :: vtx
    integer :: el, locf, globf, adjglobf, i

    integer, dimension(:,:), allocatable :: locf2locv
    
    faces = 3
    vtx   = size( Th%el2v, 2 )

    allocate( locf2locv(faces, 2) ) ! two vertices per face

    locf2locv(1,:) = (/ 1 , 2 /)
    locf2locv(2,:) = (/ 2 , 3 /)
    locf2locv(3,:) = (/ 3 , 1 /)

    Th%n_faces = Th%n_el*faces
    
    allocate( Th%el2face( Th%n_el, faces ) )
    allocate( Th%face2face( Th%n_faces ) )
    allocate( Th%face2locf( Th%n_faces ) )
    allocate( Th%face2v( Th%n_faces, 2 ) )

    ! build el2f, f2v
    i = 1
    do el=1, Th%n_el
       do locf=1, faces
          Th%el2face(el,locf) = (el-1)*faces + locf
          Th%face2v(i,:) = Th%el2v( el, locf2locv(locf,:) )
          !write(6,*) locf2locv(locf,:)
          Th%face2locf(i) = locf
          i = i + 1
       enddo
    enddo

    ! build face2face
    Th%face2face = -1
    do globf=1, Th%n_faces
       if ( Th%face2face(globf) .eq. -1 ) then

          do adjglobf=globf+1, Th%n_faces
             if ( (Th%face2v(globf,1).eq.Th%face2v(adjglobf,1) .and. &
                   Th%face2v(globf,2).eq.Th%face2v(adjglobf,2) )     &
                   .or.                                              &
                  (Th%face2v(globf,1).eq.Th%face2v(adjglobf,2) .and. &
                   Th%face2v(globf,2).eq.Th%face2v(adjglobf,1) )     &
                  ) then
                !(Th%face2v(globf,:) .eq. Th%face2v(adjglobf,(/2,1/)) ) &
                Th%face2face(globf)    = adjglobf
                Th%face2face(adjglobf) = globf
             end if
          end do
          
       end if
    end do
    
  end subroutine buildfacetri2d

  !> builds the map of adjacent elements
  subroutine buildel2eltri2d(Th)
    type(tri2d), intent(inout) :: Th

    integer :: faces, vtx ! number of faces of an element

    integer :: el, locf, globf, adjglobf, adjel

    faces = 3; vtx = 3
    
    allocate( Th%el2el( Th%n_el, faces ) )
    allocate( Th%centroid( Th%n_el, Th%n_dim ) )

    Th%el2el = -1               ! init

    do el=1, Th%n_el
       do locf=1, faces
          globf    = Th%el2face(el,locf)
          adjglobf = Th%face2face(globf)

          if (adjglobf.ne.-1) then
             adjel = int( ceiling( real(adjglobf)/faces ) )
             ! from definition
             Th%el2el(el,locf) = adjel
          end if
       end do
    end do

    do el=1, Th%n_el
       Th%centroid(el,1) = sum( Th%v2corr( Th%el2v(el,:), 1 ) )/real(vtx)
       Th%centroid(el,2) = sum( Th%v2corr( Th%el2v(el,:), 2 ) )/real(vtx)
    enddo

  end subroutine buildel2eltri2d

  !> build edges 
  subroutine buildedgestri2d(Th)
    type(tri2d), intent(inout) :: Th

    integer :: faces ! number of faces of an element
    integer :: vtx

    integer :: locf, adjlocf, globf, adjglobf, globedge, bdedge, inedge
    integer :: el, adjel
    integer, dimension(Th%n_faces) :: faceflag

    globedge = 0
    bdedge   = 0
    inedge   = 0
    faceflag = -1

    faces = size( Th%el2face, 2 )
    vtx   = size( Th%face2v, 2 )

    Th%n_edge  = 0
    Th%n_Bedge = 0
    Th%n_Iedge = 0

    do globf=1, Th%n_faces
       if ( Th%face2face(globf).eq.-1 ) then
          Th%n_Bedge = Th%n_Bedge + 1
       else
          Th%n_Iedge = Th%n_Iedge + 1
       end if
    end do

    if ( mod(Th%n_Iedge,2).eq.0 ) then
       Th%n_Iedge = Th%n_Iedge/2
    else
       write(6,*) "error in triangulate2d: buildedgestri2d"
       stop
    end if

    Th%n_edge = Th%n_Iedge + Th%n_Bedge

    allocate( Th%edge2el( Th%n_edge, 2 ) ) ! max two adjacent elements
    allocate( Th%el2edge( Th%n_el, faces ) ) ! 
    allocate( Th%edge2face( Th%n_edge, 2 ) ) ! max two adjacent faces
    allocate( Th%edge2v( Th%n_edge, vtx ) )
    allocate( Th%edge2bd( Th%n_edge ) )
    allocate( Th%inedge_idx( Th%n_Iedge ) )
    allocate( Th%bdedge_idx( Th%n_Bedge ) )

    Th%edge2el = -1             ! init
    Th%el2edge = -1
    Th%edge2bd = -1             ! init

    Th%inedge_idx = -1
    Th%bdedge_idx = -1

    do globf=1, Th%n_faces
       if ( faceflag(globf).eq.-1 ) then

          globedge        = globedge + 1
          faceflag(globf) = 0

          Th%edge2v(globedge,:) = Th%face2v(globf,:)

          adjglobf = Th%face2face( globf )

          ! 
          el = int( ceiling( real(globf)/faces ) )
          locf = mod( globf, faces )
          if (locf.eq.0) then
             locf = faces
          endif
          Th%edge2el(globedge,1) = el
          Th%el2edge(el,locf) = globedge
          Th%edge2face(globedge,1) = globf
          ! 

          if ( adjglobf.ne.-1 ) then
             inedge = inedge + 1
             faceflag(adjglobf) = 0
             ! 
             adjel = int( ceiling( real(adjglobf)/faces ) )
             adjlocf = mod( adjglobf, faces )
             if (adjlocf.eq.0) then
                adjlocf = faces
             endif
             Th%edge2el(globedge,2) = adjel
             Th%el2edge(adjel,adjlocf) = globedge
             Th%edge2face(globedge,2) = adjglobf
             Th%inedge_idx( inedge ) = globedge
             ! 
          else                  ! if it is on the boundary
             bdedge = bdedge + 1
             Th%edge2face(globedge,2) = -1
             Th%edge2bd(globedge) = 0
             Th%bdedge_idx( bdedge ) = globedge
          end if

       end if
    end do

  end subroutine buildedgestri2d
  
  

  !> builds the maps for triangulation.
  !> It is necessary to have maps for elements to vertices
  !> and vectices to coordinates
  subroutine buildtri2d(Th)
    type(tri2d), intent(inout) :: Th

    call buildfacetri2d(Th)
    call buildel2eltri2d(Th)
    call buildedgestri2d(Th)

    Th%init = .true.
  end subroutine buildtri2d
  
  !> copy mesh T1h into T2h. T1h components has to be filled
  subroutine copymesh(T2h, T1h)
    type(tri2d), intent(in)    :: T1h
    type(tri2d), intent(out) :: T2h

    T2h%n_dim = T1h%n_dim
    T2h%n_el  = T1h%n_el
    T2h%n_vtx = T1h%n_vtx

    T2h%n_faces = T1h%n_faces

    T2h%n_edge   = T1h%n_edge
    T2h%n_Bedge  = T1h%n_Bedge
    T2h%n_Iedge  = T1h%n_Iedge    

    allocate( T2h%el2v( size(T1h%el2v,1) , size(T1h%el2v,2) ) )
    allocate( T2h%v2corr( size(T1h%v2corr,1) , size(T1h%v2corr,2) ) )

    allocate( T2h%el2el( size(T1h%el2el,1) , size(T1h%el2el,2) ) )

    allocate( T2h%el2face( size(T1h%el2face,1) , size(T1h%el2face,2) ) )
    allocate( T2h%face2v( size(T1h%face2v,1) , size(T1h%face2v,2) ) )
    allocate( T2h%face2face( size(T1h%face2face) ) )

    allocate( T2h%edge2el( size(T1h%edge2el,1) , size(T1h%edge2el,2) ) )
    allocate( T2h%edge2v( size(T1h%edge2v,1) , size(T1h%edge2v,2) ) )
    allocate( T2h%edge2bd( size(T1h%edge2bd) ) )

    T2h%el2v   = T1h%el2v
    T2h%v2corr = T1h%v2corr

    T2h%el2el   = T1h%el2el

    T2h%el2face   = T1h%el2face
    T2h%face2v    = T1h%face2v
    T2h%face2face = T1h%face2face

    T2h%edge2el = T1h%edge2el
    T2h%edge2v  = T1h%edge2v
    T2h%edge2bd = T1h%edge2bd
  end subroutine copymesh
  
  
end module triangulate2d
