!> OSM.f90 provides tools necessary to assemble OSM for IPH method
!!  

module OptSchwarz

  use alldgspace2d
  use intfspace2d
  use mylinalgebra
  
  implicit none

  !> restriction operator
  type Rest
     integer, allocatable :: dof(:)
  end type Rest
  
  !> a type which contains tools to assemble the global matrix for OSM
  type OSMtools

     logical :: init = .false.

     !! number of dofs of each component
     integer :: n_dof_U
     integer :: n_dof_R
     integer :: n_dof_tot

     !! an array with the indeces of u and r on the big matrix
     integer, allocatable :: index_u(:)
     integer, allocatable :: index_r(:)

     !! a map from local dofs of the auxiliary variable along the
     !! interface into the global matrix dofs
     integer, allocatable :: AUXloc2glob(:,:,:)

     !! Restriction operator for each subdomain dofs. It acts on dofs
     !! of the primal variable: u.
     type(Rest), allocatable :: Res(:)

  end type OSMtools

contains

  subroutine buildOSM(V, auxXi, nest2d, osm)
    type(primaldg2d), intent(in) :: V
    type(intfsp2d),   intent(in) :: auxXi
    type(nestmesh2d), intent(in) :: nest2d

    type(OSMtools),  intent(inout) :: osm

    integer :: n_dof_subd, subd, el, glob_el, i

    integer, allocatable :: locdof(:)

    integer :: Np, NpF
    integer :: Nside
    integer :: NintfEdge

    if ( osm%init.eqv..true. ) then
       write(6,*) "error in osm%init", osm%init
    elseif ( V%init.eqv..false. ) then
       write(6,*) "error in V%init", V%init
    elseif ( auxXi%init.eqv..false. ) then
       write(6,*) "error in auxXi%init", auxXi%init
    endif

    Np  = V%sc%Np
    NpF = V%trSc%NpF

    Nside     = 2               ! two side per subdomain
    NintfEdge = auxXi%NintfEdge ! number of interface edges

    allocate( osm%AUXloc2glob(NintfEdge, Nside, NpF) )

    osm%AUXloc2glob(:,1,:) = auxXi%sltr(1)%edge2idof ! first side indices
    osm%AUXloc2glob(:,2,:) = osm%AUXloc2glob(:,1,:) + NpF*NintfEdge ! second
    ! side
    ! with
    ! shift
    osm%n_dof_U = V%sc%n_dof
    osm%n_dof_R = Nside * NpF * NintfEdge

    osm%n_dof_tot = osm%n_dof_R + osm%n_dof_U

    allocate( osm%index_u(osm%n_dof_U), osm%index_r(osm%n_dof_R) )

    osm%index_u = (/ (i, i=1, osm%n_dof_U) /)
    osm%index_r = (/ (i, i=1, osm%n_dof_R) /) + osm%n_dof_U
    
    allocate( osm%Res(nest2d%Nsubd) )
    allocate( locdof( Np ) )

    !! building dofs of each subdomain
    do subd=1, nest2d%Nsubd
       n_dof_subd = Np*nest2d%subd2el(subd)%Nele_subd
       allocate( osm%Res(subd)%dof(n_dof_subd) )
       !
       do el=1, nest2d%subd2el(subd)%Nele_subd
          glob_el = nest2d%subd2el(subd)%ele(el)
          locdof = (/ (i, i=((el-1)*Np+1), el*Np ) /)
          osm%Res(subd)%dof( locdof ) = V%sc%el2dof(glob_el,:)
       enddo
       !
    enddo

    osm%init = .true.

  end subroutine buildOSM

  !! build the block diagonals of the original matrix A
  subroutine blockdiag(A, osm, out)
    real, dimension(:,:), intent(in)    :: A
    type(osmtools),       intent(in)    :: osm

    real, allocatable,    intent(inout) :: out(:,:)

    integer :: subd

    if ( allocated( out ) ) then
       write(6,*) "error in blockdiag"
       stop
    endif

    allocate( out(size(A,1), size(A,2)) )
    out = 0.0

    write(6,*) size( osm%Res )
    
    do subd=1, size( osm%Res )
       out( osm%Res(subd)%dof, osm%Res(subd)%dof ) = &
            out( osm%Res(subd)%dof, osm%Res(subd)%dof ) + &
            A( osm%Res(subd)%dof, osm%Res(subd)%dof )
    enddo
    

  end subroutine blockdiag

  
  !> performing blockJacobi
  subroutine blockJacobi(Adiag, Aoff, rhs, v, Nite)
    real, dimension(:,:), intent(in) :: Adiag, Aoff
    real, dimension(:), intent(in) :: rhs
    real, dimension(:), intent(inout) :: v
    integer, intent(inout) :: Nite
  
    real, parameter :: eps = 10.0**(-1)

    real :: initErr

    real, dimension( size(Adiag,1) ) :: Un, Um, Ainvrhs

    real, dimension( size(Adiag,1), size(Adiag,2) ) :: Ainv

    Ainv = inv(Adiag)
    
    Un = 0.0
    Um = 0.0

    Nite = 0

    initErr = sqrt( sum(V**2) )

    Ainvrhs = matmul(Ainv,rhs)
    
    do while ( sqrt( sum((Um-V)**2) )/initErr.ge.eps )
       Un = Ainvrhs - matmul( matmul(Ainv, Aoff), Um )
       Um = Un
       Nite = Nite + 1
       write(6,*) Nite, sqrt( sum((Um-V)**2) )/initErr
    enddo

    v = Um
  end subroutine blockJacobi

end module OptSchwarz
