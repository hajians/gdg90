  !> normdg2d.f90

module normdg2d
  use dgspace2d
  use dgfluxspace2d
  use dgopt2d

  implicit none

contains

  !> L2norm of U
  function L2norm(U) result(out)
    type(dg2d), intent(in) :: U
    real :: out

    out = dot0Th(U,U)

    out = sqrt(out)
  end function L2norm

  !> DG norm
  function DGnorm(U) result(out)
    type(dg2d), intent(in) :: U
    real :: out
    
    type(dgflux2d)  :: Du
    type(tr2d)      :: trU
    type(paramtr2d) :: mu 

    call copydg2fluxdg(U,Du)
    call buildtr2d(U%Th, U%lo1d, U%k, trU)
    call link(U,trU)
    
    Du = Dh(U)
    call trace(U,trU)

    mu = penaltyparam(trU)
    
    out = dot0Th(Du,Du) + dot01d( mu.dot.jump(trU), jump(trU) )

    out = sqrt(out)
  end function DGnorm

  !> penalty param
  function penaltyparam(trV) result(mu)
    type(tr2d), intent(in) :: trV
    type(paramtr2d) :: mu

    integer :: edge, globf

    if ( trV%init .eqv. .false. ) then
       write(6,*) "error in penalty"
       stop
    endif

    call buildparamtr2d(trV%Th,mu)

    do edge=1, trV%Th%n_edge

       globf = trV%Th%edge2face(edge,1)
       if (globf > 0 ) then
          mu%face2param(globf) = 1.0/trV%edge2J1D(edge)
       endif

       globf = trV%Th%edge2face(edge,2)
       if (globf > 0 ) then
          mu%face2param(globf) = 1.0/trV%edge2J1D(edge)
       endif

    enddo

    mu%face2param = mu%face2param * 4*(trV%k+1)*(trV%k+2)
    mu%init = .true.
  end function penaltyparam
end module normdg2d
