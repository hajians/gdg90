!> @file jacobip.f90
!> version
!> there was a bug

module JacobiP
  implicit none
  ! a module that contains the ...


contains

  function JacP(x, alpha, beta, n) RESULT(jacobi)

    real, intent(in), dimension(:) :: x
    real, intent(in) :: alpha, beta
    integer :: n

    real, dimension(size(x)) ::  jacobi
    real, dimension(:), allocatable :: p0, p1, pnew
    real :: aold, anew, bold
    integer :: length, i

    length = size(x)
    allocate( p0(length), p1(length), pnew(length) )

    p0 = sqrt( 0.5**(alpha+beta+1)*gamma(alpha+beta+2)/&
         ( gamma(alpha+1)*gamma(beta+1) ) )
    p1 = 0.5*p0*sqrt( (alpha+beta+3)/( (alpha+1)*(beta+1) ) ) * &
         ( (alpha+beta+2)*x + (alpha-beta) )
    aold = 2.*sqrt( (1+alpha)*(1+beta)/(alpha+beta+3) )/(2+alpha+beta)

    if (n > 1) then
       do i=2, n
          anew = 2.*sqrt( i*(i+alpha+beta)*(i+alpha)*(i+beta)/&
               (2.*i+alpha+beta-1)/(2.*i+alpha+beta+1) )/(2*i+alpha+beta)
          bold = -1.*(alpha**2 - beta**2)/(2.*(i-1)+alpha+beta)/&
               (2*(i-1)+alpha+beta+2) ! B old!
          pnew = (1./anew)*( (x-bold)*P1 - aold*P0)

          aold = anew
          p0 = p1
          p1 = pnew
       end do
    endif

    if (n==0) then
       jacobi = p0
    else
       jacobi = p1
    endif
  end function JacP

  function VanderMonde(x) RESULT(Van)
    real, intent(in), dimension(:) :: x
    real, dimension(size(x),size(x)) :: Van
    real, dimension(size(x)) :: row
    integer :: length, i
    length = size(x)

    do i=1, length
       row = JacP(x,0.,0.,i-1)
       Van(:,i) = row
    end do
  end function VanderMonde

  function VanderMonde1(x) RESULT(Van)
    real, intent(in), dimension(:) :: x
    real, dimension(size(x),size(x)) :: Van
    integer :: length, i
    length = size(x)
    do i=1, length
       Van(:,i) = JacP(x,0.,0.,i-1)
    enddo

  end function VanderMonde1

  FUNCTION GradJP(x,alpha,beta,n) RESULT(DP)
    real, intent(in), dimension(:) :: x
    real, intent(in) :: alpha, beta
    integer, intent(in) :: n

    real, dimension(size(x)) ::  DP
    if (n==0) then
       DP = 0.0
    else
       DP = sqrt( n*(n+alpha+beta+1) )*JacP(x,alpha+1,beta+1,n-1)
    end if

  END FUNCTION GradJP

  !> Integral of Jacobi ???
  FUNCTION IntJP(x,n) RESULT(IP)
    real, intent(in), dimension(:) :: x
    integer, intent(in) :: n
    real, dimension( size(x) ) :: IP

    IF (n==0) THEN
       IP = (x + 1.)/sqrt(2.)
    ELSE
       IP = (x**2 - 1.)/(n*(n+1)) * GradJP(x, 0., 0., n) 
    ENDIF
  END FUNCTION IntJP

  FUNCTION GradVaMonde(x) RESULT(GradVM)
    real, intent(in), dimension(:) :: x
    real, dimension(size(x),size(x)) :: GradVM 
    integer :: length, i
    length = size(x)
    do i=1, length
       GradVM(:,i) = GradJP(x,0.,0.,i-1)
    enddo

  END FUNCTION GradVaMonde



end module JacobiP

  
    
    
