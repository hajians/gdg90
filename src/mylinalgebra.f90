! Returns the inverse of a matrix calculated by finding the LU
! decomposition.  Depends on LAPACK.

module mylinalgebra
  implicit none
contains
  
function inv(A) result(Ainv)
  real, dimension(:,:), intent(in) :: A
  real, dimension(size(A,1),size(A,2)) :: Ainv

  real, dimension(size(A,1)) :: work  ! work array for LAPACK
  integer, dimension(size(A,1)) :: ipiv   ! pivot indices
  integer :: n, info

  ! External procedures defined in LAPACK
  external DGETRF
  external DGETRI

  ! Store A in Ainv to prevent it from being overwritten by LAPACK
  Ainv = A
  n = size(A,1)

  ! DGETRF computes an LU factorization of a general M-by-N matrix A
  ! using partial pivoting with row interchanges.
  call DGETRF(n, n, Ainv, n, ipiv, info)

  if (info /= 0) then
     stop 'Matrix is numerically singular!'
  end if

  ! DGETRI computes the inverse of a matrix using the LU factorization
  ! computed by DGETRF.
  call DGETRI(n, Ainv, n, ipiv, work, n, info)

  if (info /= 0) then
     stop 'Matrix inversion failed!'
  end if
end function inv

function optmatmul(A, v) result(out)
  real, intent(in), dimension(:,:) :: A
  real, intent(in), dimension( size(A,2) ) :: v

  real, dimension( size(A,1) ) :: out

  integer :: j
  real    :: eps = 10.0**(-9)

  out = 0.0
  
  do j=1, size(A,2)
     if ( abs(v(j)).ge.eps ) then
       out = out + A(:,j) * v(j)
     endif
  enddo

end function optmatmul

!>
function sc_product(u,A,v) result(out)
  real, intent(in), dimension(:,:) :: A
  real, intent(in), dimension( size(A,2) ) :: v
  real, intent(in), dimension( size(A,1) ) :: u

  integer :: isU, isV, idU, idV
  real :: out

  out = 0.0

  call islocbinary(u,isU,idU)
  call islocbinary(v,isV,idV)

  if ( (isU.eq.1) .and. (isV.eq.1) ) then
     out = A(idU,idV)
     !write(6,*) "yes"
  else
     out = dot_product( u, matmul(A,v) )
  endif

end function sc_product

!>
subroutine islocbinary(v,flag,idx)
  real, intent(in) :: v(:)

  integer, intent(inout) :: flag
  integer, intent(inout) :: idx

  integer :: len, i
  real :: epsP = +10.0**(-8)
  real :: epsN = -10.0**(-8)

  flag = 0
  len = size(v)

  do i=1, len
     if ( (v(i).ge.epsP) .or. (v(i).ge.epsN) ) then
        flag = flag+1
        idx  = i
     endif
  enddo

end subroutine islocbinary


subroutine findev(A, ev, B)
  real, dimension(:,:),           intent(in)    :: A
  real, optional, dimension(:,:), intent(in)     :: B
  real, dimension( size(A,1) ),   intent(inout) :: ev

  integer :: len, lda, ldb, info 

  integer :: lwork, itype
  real, allocatable :: work(:)

  real, allocatable :: dummA(:,:)
  real, allocatable :: dummB(:,:)

  EXTERNAL         DSYEV
  EXTERNAL         DSYGV

  allocate( dummA(size(A,1), size(A,2)) )
  dummA = A

  len = size(A,1)
  lda = len
  ldb = len
  ev  = 0.0

  lwork = 3*len;
  allocate( work( lwork ) )

  itype = 1                     ! find eigenvalues of A*v = lam*B*v
  !write(6,*) size(work)
  if ( present(B) ) then
     allocate( dummB(size(B,1), size(B,2)) )
     dummB = B
     CALL DSYGV(itype,'N','U',len,dummA,LDA,dummB,LDB, ev,WORK,LWORK,INFO)
  else
     CALL DSYEV('N','U',len,dummA,LDA,ev,WORK,LWORK,INFO)
  endif

  if (info.ne.0) then
     write(6,*) "error in findev:", INFO
  endif

end subroutine findev

end module mylinalgebra

