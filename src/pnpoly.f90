module pnpoly
  implicit none
contains
  
subroutine polygon_contains_point_2d ( n, v, p, inside )

  !*****************************************************************************80
  !
  !! POLYGON_CONTAINS_POINT_2D finds if a point is inside a simple polygon in 2D.
  !
  !  Discussion:
  !
  !    A simple polygon is one whose boundary never crosses itself.
  !    The polygon does not need to be convex.
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license. 
  !
  !  Modified:
  !
  !    09 May 2005
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Reference:
  !
  !    M Shimrat,
  !    ACM Algorithm 112,
  !    Position of Point Relative to Polygon,
  !    Communications of the ACM,
  !    Volume 5, Number 8, page 434, August 1962.
  !
  !  Parameters:
  !
  !    Input, integer ( kind = 4 ) N, the number of nodes or vertices in 
  !    the polygon.  N must be at least 3.
  !
  !    Input, real ( kind = 8 ) V(2,N), the vertices of the polygon.
  !
  !    Input, real ( kind = 8 ) P(2), the coordinates of the point to be tested.
  !
  !    Output, logical INSIDE, is TRUE if the point is inside the polygon.
  !
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ), parameter :: dim_num = 2

  integer ( kind = 4 ) i
  logical inside
  integer ( kind = 4 ) ip1
  real ( kind = 8 ) p(dim_num)
  real ( kind = 8 ) v(dim_num,n)

  inside = .false.

  do i = 1, n

     if ( i < n ) then
        ip1 = i + 1
     else
        ip1 = 1
     end if

     if ( ( v(2,i)   <  p(2) .and. p(2) <= v(2,ip1)   ) .or. &
          ( p(2) <= v(2,i)   .and. v(2,ip1)   < p(2) ) ) then
        if ( ( p(1) - v(1,i) ) - ( p(2) - v(2,i) ) &
             * ( v(1,ip1) - v(1,i) ) / ( v(2,ip1) - v(2,i) ) < 0 ) then
           inside = .not. inside
        end if
     end if

  end do

  return
end subroutine polygon_contains_point_2d

! is a point contains between a line passing two points
subroutine point_contains_line_2d(v,p,inside)
  real, dimension(2,2), intent(in) :: v
  real, dimension(2),   intent(in) :: p

  logical, intent(inout) :: inside

  real, dimension(2) :: vec1, vec2
  real :: d1, d2
  real :: eps = 10.0**(-9)

  inside = .false.

  vec1 = v(1,:) - p; d1 = sqrt(dot_product( vec1, vec1 ));
  vec2 = v(2,:) - p; d2 = sqrt(dot_product( vec2, vec2 ));

  if ( (abs(d1).le.eps) .or. (abs(d2).le.eps) ) then
     inside = .true.
  else
     vec1 = vec1 / d1
     vec2 = vec2 / d2

     if ( abs( dot_product(vec1,vec2) - (-1.0) ).le.eps ) then
        inside = .true.
     endif
  endif
end subroutine point_contains_line_2d
end module pnpoly
