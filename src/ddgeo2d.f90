!> @file ddgeo2d.f90 domain decomposition geometry in 2-d

module ddgeo2d

  implicit none
  
  !> convex polygon
  type cpoly

     logical :: init = .false.
     
     integer :: n_vert
     integer :: n_dim

     real, allocatable :: v2corr(:,:)
          
  end type cpoly

  !> geometry and domain decomposition into convex polygons in 2-d
  type geom2d

     logical :: init = .false.

     integer :: n_subd          ! number of subdomains

     type(cpoly), allocatable :: subd(:) ! subdomains
     type(cpoly)              :: domain  ! domain
     
  end type geom2d

  
contains

  !> build a convex polygon object
  subroutine buildcpoly(poly, dim, n_v, v2corr)
    type(cpoly), intent(inout) :: poly

    integer, intent(in) :: dim, n_v
    real,    intent(in) :: v2corr(n_v,dim)

    if (poly%init.eqv..true.) then
       write(6,*) "error in buildcpoly: poly is initialized already"
       stop
    endif

    poly%n_dim  = dim
    poly%n_vert = n_v

    allocate( poly%v2corr(n_v, dim) )
    poly%v2corr = v2corr

    poly%init = .true.
  end subroutine buildcpoly

  !> read the geometry
  subroutine readgeo2d(filename, geo2d)
    character(len=*), intent(in)    :: filename
    type(geom2d),     intent(inout) :: geo2d

    integer :: myunit = 11

    integer :: tot_vtx, dim, dummy_vtx
    real,    allocatable :: v2corr(:,:), dummycorr(:)
    integer, allocatable :: dummyv(:)

    integer :: i, j

    if (geo2d%init.eqv..true.) then
       write(6,*) "error in readgeo2d: geo2d is already initi."
    endif

    i = 0
    j = 0

    open(myunit, file=filename//'.geo')

    read(myunit,*) tot_vtx, dim
    allocate( v2corr(tot_vtx,dim), dummycorr(dim) )

    do i=1, tot_vtx             ! read all vertices
       read(myunit,*) j, dummycorr
       v2corr(j,:) = dummycorr
    enddo

!!!!!
    read(myunit,*) dummy_vtx   ! vertices of the domain
    allocate( dummyv( dummy_vtx ) )
    read(myunit,*) dummyv       ! read vertices constructing the domain
    call buildcpoly(geo2d%domain, dim, dummy_vtx, v2corr(dummyv,:) )
    deallocate( dummyv )
!!!!!

    read(myunit,*) geo2d%n_subd ! number of subdomain
    allocate( geo2d%subd( geo2d%n_subd ) )

    do i=1, geo2d%n_subd
       read(myunit,*) dummy_vtx   ! vertices of the subdomain
       allocate( dummyv( dummy_vtx ) )
       read(myunit,*) dummyv       ! read vertices constructing the subdomain
       call buildcpoly(geo2d%subd(i), dim, dummy_vtx, v2corr(dummyv,:) )
       deallocate( dummyv )
    enddo

    close(myunit)

    geo2d%init = .true.

  end subroutine readgeo2d
  
  subroutine plotgeom2d(filename,geo2d)
    character(len=*), intent(in) :: filename
    type(geom2d),     intent(in) :: geo2d

    integer :: i, sub
    integer :: myunit = 11

    if (geo2d%init.eqv..false.) then
       write(6,*) "error in plotgeom2d: geo2d is not initi."
       stop
    endif

    open(myunit, file=filename//'D.gnu') ! Domain geometry file

    do i=1, geo2d%domain%n_vert
       write(myunit,*) geo2d%domain%v2corr(i,:)
    enddo
    write(myunit,*) geo2d%domain%v2corr(1,:)

    close(myunit)

    open(myunit, file=filename//'S.gnu') ! subdomain geometry file
    do sub=1, geo2d%n_subd
       do i=1, geo2d%subd(sub)%n_vert
          write(myunit,*) geo2d%subd(sub)%v2corr(i,:)
       enddo
       write(myunit,*) geo2d%subd(sub)%v2corr(1,:)
       write(myunit,*) 
    enddo

    close(myunit)
    
  end subroutine plotgeom2d
  
  
end module ddgeo2d

     
  
