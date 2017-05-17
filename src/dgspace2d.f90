  !> @file dgspace2d.f90
  !> @brief contains objects which represent scalar discontinuous
  !> Galerkin finite element in 2-dimension.contains

module dgspace2d
  use triangulate2d
  use localopt1d
  use localopt2d
  
  implicit none

  type dg2d

     logical :: init = .false.
     
     type(tri2d),    pointer :: Th   => null()
     type(locopt1d), pointer :: lo1d => null()
     type(locopt2d), pointer :: lo2d => null()
     
     integer :: k               ! deg of polynomial

     integer       :: Np
     integer       :: n_dof
     real, allocatable :: dof(:)

     integer, allocatable :: el2dof(:,:)
     real,    allocatable :: dof2corr(:,:)

     real, allocatable :: rx(:,:)
     real, allocatable :: sx(:,:)
     real, allocatable :: ry(:,:)
     real, allocatable :: sy(:,:)

     real, allocatable :: el2J(:,:)

     real, allocatable :: nx(:)
     real, allocatable :: ny(:)

     logical :: binary    = .false. ! if the nodal values is zero except one element
     integer :: elnotzero = 0 

  end type dg2d

  private builddof
  private geofac2d
contains

  !> build the map from element index to dofs
  !! representing an object in dg2d space
  !> @param [in] lo2d
  !> @param [inout] V is an object in the dg2d space
  subroutine builddof(lo2d, V)
    type(dg2d),     intent(inout) :: V
    type(locopt2d), intent(in)    :: lo2d
    
    integer :: i, el
    integer :: faces, vtx

    integer, allocatable, dimension(:) :: locvtx

    real, allocatable, dimension(:) :: r, s

    faces  = size(V%Th%el2face, 2)
    vtx    = size(V%Th%el2v, 2)

    allocate( locvtx(faces) )

    V%Np    = (V%k+1)*(V%k+2)/2    ! dofs of a triangle
    V%n_dof = V%Np * V%Th%n_el

    allocate( r(V%Np), s(V%Np) )
    r = lo2d%refxy(:,1)
    s = lo2d%refxy(:,2)
    
    allocate( V%dof( V%n_dof ) )
    allocate( V%el2dof( V%Th%n_el, V%Np) )
    allocate( V%dof2corr( V%n_dof, V%Th%n_dim ) )

    V%dof = 0.0

    do el=1, V%Th%n_el
       V%el2dof(el,:) = (/ (i, i=(el-1)*V%Np+1,el*V%Np) /)

       locvtx = V%Th%el2v(el,:)

       V%dof2corr( V%el2dof(el,:), 1 ) = &
            0.5*( -(r+s) * V%Th%v2corr( locvtx(1), 1) + &
            (1+r) * V%Th%v2corr( locvtx(2), 1) + &
            (1+s) * V%Th%v2corr( locvtx(3), 1) )
       V%dof2corr( V%el2dof(el,:), 2 ) = &
            0.5*( -(r+s) * V%Th%v2corr( locvtx(1), 2) + &
            (1+r) * V%Th%v2corr( locvtx(2), 2) + &
            (1+s) * V%Th%v2corr( locvtx(3), 2) )
    end do

  end subroutine builddof

  subroutine geofac2d(lo2d,V)
    type(locopt2d), intent(in)    :: lo2d
    type(dg2d),     intent(inout) :: V

    real, dimension( V%Np ) :: xr, xs, yr, ys
    real, dimension( V%Np ) :: x, y

    integer :: el, face, globf
    real    :: norm

    allocate( V%rx( V%Np, V%Th%n_el ) )
    allocate( V%sx( V%Np, V%Th%n_el ) )
    allocate( V%ry( V%Np, V%Th%n_el ) )
    allocate( V%sy( V%Np, V%Th%n_el ) )

    allocate( V%el2J( V%Np, V%Th%n_el ) )

    allocate( V%nx( V%Th%n_faces ) )
    allocate( V%ny( V%Th%n_faces ) )

    do el=1, V%Th%n_el
       x = V%dof2corr( V%el2dof(el,:), 1 )
       y = V%dof2corr( V%el2dof(el,:), 2 )

       xr = matmul(lo2d%dr2d, x)
       xs = matmul(lo2d%ds2d, x)
       yr = matmul(lo2d%dr2d, y)
       ys = matmul(lo2d%ds2d, y)

       V%el2J(:,el) = -xs*yr + xr*ys

       if (V%el2J(1,el)<0.0 ) then
          write(6,*) "negative J"
          stop
       end if

       V%rx(:,el) =  ys / V%el2J(:,el)
       V%sx(:,el) = -yr / V%el2J(:,el)
       V%ry(:,el) = -xs / V%el2J(:,el)
       V%sy(:,el) =  xr / V%el2J(:,el)

       face  = 1
       globf = V%Th%el2face(el,face) 
       V%nx(globf) = -V%sx(1,el)
       V%ny(globf) = -V%sy(1,el)

       norm = sqrt(V%nx(globf)**2 + V%ny(globf)**2)
       V%nx(globf) = V%nx(globf)/norm
       V%ny(globf) = V%ny(globf)/norm
       

       face  = 2
       globf = V%Th%el2face(el,face) 
       V%nx(globf) = V%sx(1,el) + V%rx(1,el)
       V%ny(globf) = V%sy(1,el) + V%ry(1,el)

       norm = sqrt(V%nx(globf)**2 + V%ny(globf)**2)
       V%nx(globf) = V%nx(globf)/norm
       V%ny(globf) = V%ny(globf)/norm


       face  = 3
       globf = V%Th%el2face(el,face) 
       V%nx(globf) = -V%rx(1,el)
       V%ny(globf) = -V%ry(1,el)

       norm = sqrt(V%nx(globf)**2 + V%ny(globf)**2)
       V%nx(globf) = V%nx(globf)/norm
       V%ny(globf) = V%ny(globf)/norm
    end do

  end subroutine geofac2d
  
  !> given Th and k it builds a null element
  !! of DG finite element space.
  !> @param [in] Th is the mesh object
  !> @param [in] k is the degree of polynomials
  function builddg2d(Th, lo1d, lo2d, k) result(V)
    type(tri2d), target,    intent(in) :: Th
    type(locopt1d), target, intent(in) :: lo1d
    type(locopt2d), target, intent(in) :: lo2d
    integer,                intent(in) :: k

    type(dg2d) :: V

    integer :: faces, vtx
    !    integer :: Np, NpEdge, NpInEdge, NpIn

    V%Th   => Th                ! assignment should be rigorous
    V%lo1d => lo1d
    V%lo2d => lo2d
    V%k  = k

    faces  = size(V%Th%el2face, 2)
    vtx    = size(V%Th%el2v, 2)

    call builddof(lo2d, V)

    call geofac2d(lo2d, V)

    V%init = .true.
  end function builddg2d

  subroutine copydg2dg(U,V)
    type(dg2d), intent(in)    :: U
    type(dg2d), intent(inout) :: V

    if ((V%init .eqv. .true.) .or. &
         (U%init .eqv. .false.) ) then
       write(6,*) "error in copydg2dg:", V%init, U%init
       stop
    endif

    V%Th   => U%Th
    V%lo1d => U%lo1d
    V%lo2d => U%lo2d

    V%k = U%k

    V%Np = U%Np
    V%n_dof = U%n_dof

    allocate( V%dof( V%n_dof ) )
    allocate( V%el2dof( V%Th%n_el, V%Np) )
    allocate( V%dof2corr( V%n_dof, V%Th%n_dim ) )

    allocate( V%rx( V%Np, V%Th%n_el ) )
    allocate( V%sx( V%Np, V%Th%n_el ) )
    allocate( V%ry( V%Np, V%Th%n_el ) )
    allocate( V%sy( V%Np, V%Th%n_el ) )

    allocate( V%el2J( V%Np, V%Th%n_el ) )

    allocate( V%nx( V%Th%n_faces ) )
    allocate( V%ny( V%Th%n_faces ) )

    V%dof = U%dof

    V%el2dof = U%el2dof
    V%dof2corr = U%dof2corr

    V%rx = U%rx
    V%sx = U%sx
    V%ry = U%ry
    V%sy = U%sy

    V%el2J = U%el2J

    V%nx = U%nx
    V%ny = U%ny

    V%binary    = U%binary
    V%elnotzero = U%elnotzero
    
    V%init = .true.

  end subroutine copydg2dg
  

end module dgspace2d
