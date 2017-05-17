!> @file femspace2d.f90
!> @brief contains objects which represent scalar conforming finite element
!> in 2-dimension.

module femspace2d
  use triangulate2d
  use localopt2d
  
  implicit none

  
  type fem2d
     type(tri2d) :: Th

     integer :: k               ! deg of polynomial

     integer       :: Np
     integer       :: n_dof
     real, pointer :: dof(:)

     integer, pointer :: el2dof(:,:)
     real,    pointer :: dof2corr(:,:)
     
  end type fem2d
  
contains

  !> build the map from element index to dofs
  !! representing an object in fem2d space
  !> @param [inout] V is an object in the fem2d space
  subroutine builddof(lo2d, V)
    type(locopt2d), intent(in)    :: lo2d 
    type(fem2d),    intent(inout) :: V

    integer*2, allocatable, dimension(:) :: flagdof
    
    integer :: faces, vtx
    integer :: Np, NpEdge, NpInEdge, NpIn

    faces  = size(V%Th%el2face, 2)
    vtx    = size(V%Th%el2v, 2)
    
    Np       = (V%k+1)*(V%k+2)/2    ! dofs of a triangle
    NpEdge   = (V%k+1)            ! dofs on an edge
    NpInEdge = NpEdge - 2       ! dofs inside an edge
    NpIn     = Np - (vtx + faces*NpInEdge)

    V%Np    = Np
    V%n_dof = V%Th%n_vtx + V%Th%n_edge*NpInEdge + V%Th%n_el*NpIn

    allocate( V%dof( V%n_dof ) )
    allocate( V%el2dof( V%Th%n_el, Np) )
    allocate( V%dof2corr( V%n_dof, V%Th%n_dim ) )

    allocate( flagdof( V%n_dof ) )

    V%dof   = 0.0
    flagdof = -1

    
    
  end subroutine builddof
  

  !> given Th and k it builds a null element
  !! of conforming finite element space.
  !> @param [in] Th is the mesh object
  !> @param [in] k is the degree of polynomials
  function buildfem2d(Th, lo2d, k) result(V)
    type(tri2d),    intent(in) :: Th
    type(locopt2d), intent(in) :: lo2d
    integer,        intent(in) :: k

    type(fem2d) :: V

    integer :: faces, vtx
!    integer :: Np, NpEdge, NpInEdge, NpIn
    
    V%Th = Th                   ! assignment should be rigorous
    V%k  = k

    faces  = size(V%Th%el2face, 2)
    vtx    = size(V%Th%el2v, 2)
    
    call builddof(lo2d, V)
    
  end function buildfem2d
  
end module femspace2d

