# Introduction

GDG90 contains libraries to implement different discontinuous Galerkin
(DG) finite element methods in an object-oriented fashion while
maintaining the speed during assembling of the underlying matrices.

The user only needs to define the bilinear forms and/or norms
associated to the DG method. For instance in order to compute the
L2-norm of a discontinuous function, one writes:

```fortran
function L2norm(U) result(out)
	type(dg2d), intent(in) :: U ! a function in DG space
	real :: out                 ! a real value for norm

	out = dot0Th(U,U)           ! computes L2-norm by 
                            ! performing inner product
	out = sqrt(out)
end function L2norm
```


# Example

Compile the main.f90 using:
```bash
make main
```

This file implements the stiffness matrix associated to interior
penalty method (which is one of the popular DG methods) and writes it
as a sparse matrix into a file. 
