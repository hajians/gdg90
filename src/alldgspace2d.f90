!> alldgspace2d.f90 provides essential components
!! to be used in a bilinear forms. The components will be
!! initialized before hand and therefore computing
!! a bilinear form will be fast

module alldgspace2d

  use dgspace2d
  use dgfluxspace2d
  use brokentracespace2d
  use singletracespace2d
  
  use dgopt2d

  implicit none

  !> primal dg space which contains all components
  type primaldg2d

     logical :: init

     integer :: k

     type(locopt1d) :: lo1d
     type(locopt2d) :: lo2d

     type(dg2d)     :: sc
     type(dgflux2d) :: flux

     type(tr2d)     :: trSc
     type(trflux2d) :: trFlux
   
     type(singletr2d)     :: singleSc
     type(singletrflux2d) :: singleFlux

     type(singletr2d)     :: paramSc
     type(singletrflux2d) :: paramFlux

     !> auxiliary trace - for dummy computations
     type(tr2d)     :: auxtrSc

     !> extra space for the jump of \nabla u. it is defined only on
     !! the interior edges.
     type(singletr2d) :: JumpFlux
     type(singletr2d) :: MuMinusOne
     
  end type primaldg2d


contains

  !> build primal dg space using triangulation Th and
  !! polynomial degree k
  subroutine buildprimaldg2d(Th, k, prdg)
    type(tri2d), intent(in) :: Th
    integer,     intent(in) :: k

    type(primaldg2d), intent(inout) :: prdg

    integer :: i, edge, globedge, dim
    integer, allocatable, dimension(:) :: locedgedof

    if ( (k.le.0) .or. (prdg%init.eqv..true.) ) then
       write(6,*) "error in buildprimaldg2d", k
       stop
    endif

    prdg%k = k

    call buildlocopt1d(k,prdg%lo1d)
    call buildlocopt2d(k,prdg%lo2d)

    ! build sc space
    prdg%sc = builddg2d(Th, prdg%lo1d, prdg%lo2d, k)
    ! build flux space
    prdg%flux = builddgflux2d(Th, prdg%lo1d, prdg%lo2d, k)

    ! build trace space for scalar space
    call buildtr2d(Th, prdg%lo1d, k, prdg%trSc)
    ! build trace space for flux space
    call buildtrflux2d(Th, prdg%lo1d, k, prdg%trFlux)

    ! link scalar with the trace of scalar
    call link(prdg%sc, prdg%trSc)
    ! link Flux with trace of Flux
    call link(prdg%flux, prdg%trFlux)

    ! build auxiliary trace space for scalar space
    call buildtr2d(Th, prdg%lo1d, k, prdg%auxtrSc)
    ! link scalar with the trace of scalar
    call link(prdg%sc, prdg%auxtrSc)

    !> build single-valued trace spaces
    call buildsingletr2d(Th,prdg%lo1d,k, (/ (i,i=1,Th%n_edge) /), prdg%singleSc)
    call buildsinglefluxtr2d(Th,prdg%lo1d,k, (/ (i,i=1,Th%n_edge) /), prdg%singleFlux)

    call buildsingletr2d(Th,prdg%lo1d,k, (/ (i,i=1,Th%n_edge) /), prdg%paramSc)
    call buildsinglefluxtr2d(Th,prdg%lo1d,k, (/ (i,i=1,Th%n_edge) /), prdg%paramFlux)

    !> link single-valued trace spaces
    call linktr(prdg%trSc, prdg%singleSc)
    call linktr(prdg%trFlux, prdg%singleFlux)

    call linktr(prdg%trSc, prdg%paramSc)
    call linktr(prdg%trFlux, prdg%paramFlux)

    !> for IPH only
    call buildsingletr2d(Th, prdg%lo1d,k, Th%inedge_idx, prdg%JumpFlux)
    call linktr(prdg%trSc, prdg%JumpFlux)
    call buildsingletr2d(Th, prdg%lo1d,k, Th%inedge_idx, prdg%MuMinusOne)
    call linktr(prdg%trSc, prdg%MuMinusOne)

    allocate( locedgedof( prdg%singleSc%NpF ) )

    do edge=1, Th%n_edge
       locedgedof = prdg%paramSc%edge2idof(edge,:)

       prdg%paramSc%edgedof(locedgedof) =50.0 / prdg%paramSc%edge2J1d(edge)
       do dim=1, Th%n_dim
          prdg%paramFlux%cp(dim)%edgedof(locedgedof) = &
               prdg%paramSc%edgedof(locedgedof)               
       enddo
    enddo


    do edge=1, Th%n_Iedge
       globedge = Th%inedge_idx(edge)
       locedgedof = prdg%MuMinusOne%edge2idof(edge,:)
       prdg%MuMinusOne%edgedof( locedgedof ) = 1.0 / &
            prdg%paramSc%edgedof( prdg%paramSc%edge2idof(globedge,:) )
    enddo

    prdg%init = .true.
  end subroutine buildprimaldg2d

end module alldgspace2d
