  !> @file interpolate2d.f90

module interpolate2d
  use dgspace2d
  implicit none

  
  abstract interface
     function fct(x) result(out)
       real, dimension(:,:), intent(in) :: x
       real, dimension( size(x,1) ) :: out
     end function fct
  end interface

  !> interface for interpolation
  interface Ih
     module procedure interplolatedg2d
  end interface Ih
    
contains

  subroutine interplolatedg2d(func, V)
    procedure(fct), pointer   :: func
    type(dg2d), intent(inout) :: V

    integer :: el
    integer, allocatable, dimension(:) :: gdof
    
    if (V%init.eqv..false.) then
       write(6,*) "error in interpolatedg2d: V is not initialized already"
       stop
    endif
    allocate( gdof( V%Np ) )
    
    do el=1, V%Th%n_el
       gdof = V%el2dof(el,:)
       V%dof(gdof) = func( V%dof2corr( gdof, : ) )
    enddo
        
  end subroutine interplolatedg2d
  
  
end module interpolate2d
