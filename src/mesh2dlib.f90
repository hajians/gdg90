!> @file mesh2d.f90

module mesh2dlib
  implicit none
  
  type triangulate2d
     integer :: n_dim
     integer :: n_el
     integer :: n_vtx
     
     integer :: n_edge
     integer :: n_Bedge
     integer :: n_Iedge
     
     integer, pointer :: el2v(:,:)
     real,    pointer :: v2corr(:,:)
     
     integer, pointer :: el2el(:)
     integer, pointer :: el2edge(:,:)

     integer, pointer :: edge2bd(:)
     
  end type triangulate2d
  
end module mesh2dlib
