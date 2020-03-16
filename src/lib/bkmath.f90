!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! BKMATH Module
! Author: Brian Krauth, Simon Fraser University
!
! Mathematical procedures related to SMLE program.  All real-valued functions 
! work in either single or double precision; results will be of same type as
! inputs.
!
! Public procedures available (by category) are:
! 
!   Probability and statistics:
!
!	PDFN(x):	Standard normal PDF of x. For PDFN, CDFN,
!			and CDFINVN, x can be scalar, vector, or
!			2-dimensional matrix.  Result will be in 
!			same type as x.  
!	CDFN(x):	Standard normal CDF of x
! 	CDFINVN(x):	Inverse standard normal CDF of x
!	CHOL(x):	Cholesky decomposition of x, assuming
!			x is symmetric and positive definite.
!	OLS(x,y):	OLS coefficients for a regression of y on x.
!	RANDINT(n):	Vector of n random integers between 1 and n.
!	
!
!   Combinatorics:
!
!	FACTORIAL(x):	Factorial of x (i.e., x!)
!	NCHOOSEK(n,k):	(n!)/(k!(n-k)!)
!
!   Halton sequences (Halton sequences give an improvement over random numbers for simulation-based
!	estimation.  See Kenneth Train's book for more information):
!
!	FIND_PRIMES(n):		Gives the first n prime numbers.  DOES IT INCLUDE 2?
!	HMAT(r,c):		Gives an r-by-c matrix of Halton sequences
!	HALTON_SEQUENCE(n,s):	Gives a Halton sequence.
!	RHALT(r,c):		Gives an r-by-c matrix of randomized Halton sequences
!			
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module bkmath  
use bklib, only : SP, DP, runif
implicit none
private
private :: pdfns0,pdfnd0,pdfns1,pdfns2,pdfnd1,pdfnd2, &
     cdfn0,cdfn1,cdfn2,cdfn0d,cdfn1d,cdfn2d, &
     cdfinvn0,cdfinvn1,cdfinvn2,cdfinvn0d,cdfinvn1d,cdfinvn2d, &
     chols,chold,olss,olsd,factorial0,factorial1,nchoosek0,nchoosek1, &
     swap,outerprod,outerand,inversed
public :: pdfn,cdfn,cdfinvn,chol,ols,randint,factorial,nchoosek,find_primes, &
     hmat,halton_sequence,rhalt,inverse
interface pdfn
   module procedure pdfns0,pdfnd0,pdfns1,pdfns2,pdfnd1,pdfnd2
end interface
interface cdfn
   module procedure cdfn0,cdfn1,cdfn2,cdfn0d,cdfn1d,cdfn2d
end interface
interface cdfinvn
   module procedure cdfinvn0,cdfinvn1,cdfinvn2,cdfinvn0d,cdfinvn1d,cdfinvn2d
end interface
interface chol
   module procedure chols,chold
end interface
interface ols
   module procedure olss,olsd
end interface
interface factorial
   module procedure factorial0,factorial1
end interface
interface nchoosek
   module procedure nchoosek0,nchoosek1
end interface
interface inverse
   module procedure inversed
end interface

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
contains
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine inversed(a)
  real(kind=DP), dimension(:,:), intent(inout) :: a
  integer, dimension(size(a,1)) :: ipiv,indxr,indxc
  logical, dimension(size(a,1)) :: lpiv
  real(kind=DP) :: pivinv
  real(kind=DP), dimension(size(a,1)) :: dumc
  integer, target, dimension(2) :: irc
  integer :: i,l,n
  integer, pointer :: irow,icol
  n=size(a,1)
  irow => irc(1)
  icol => irc(2)
  ipiv=0
  do i=1,n
     lpiv = (ipiv==0)
     irc = maxloc(abs(a),outerand(lpiv,lpiv))
     ipiv(icol)=ipiv(icol)+1
     if (ipiv(icol) > 1) stop "Singular matrix"
     if (irow /= icol) then
        call swap(a(irow,:),a(icol,:))
     end if
     indxr(i) = irow
     indxc(i) = icol
     if (a(icol,icol)==0.0_dp) stop "Singular matrix"
     pivinv=1.0_dp/a(icol,icol)
     a(icol,icol)=1.0_dp
     a(icol,:)=a(icol,:)*pivinv
     dumc=a(:,icol)
     a(:,icol)=0.0_dp
     a(icol,icol)=pivinv
     a(1:(icol-1),:)=a(1:(icol-1),:)-outerprod(dumc(1:(icol-1)),a(icol,:))
     a(icol+1:,:)=a(icol+1:,:)-outerprod(dumc(icol+1:),a(icol,:))
  end do
  do l=n,1,-1
     call swap(a(:,indxr(l)),a(:,indxc(l)))
  end do
end subroutine inversed


subroutine swap(a,b)
  real(kind=DP), dimension(:), intent(inout) :: a,b
  real(kind=DP), dimension(size(a))  :: dum
  dum=a
  a=b
  b=dum
end subroutine swap
  
function outerprod(a,b) result (outer)
  real(kind=DP), dimension(:), intent(in) :: a,b
  real(kind=DP), dimension(size(a),size(b)) :: outer
  outer = spread(a,dim=2,ncopies=size(b))*spread(b,dim=1,ncopies=size(a))
end function outerprod

function outerand(a,b) result (outer)
  logical, dimension(:), intent(in) :: a,b
  logical, dimension(size(a),size(b)) :: outer
  outer = spread(a,dim=2,ncopies=size(b)) .and. spread(b,dim=1,ncopies=size(a))
end function outerand


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! PDFN function
!
! Format: pdfn(x)
!
! Gives standard normal PDF of x.
! 
! X can be either single or double precision, and can be a scalar,
! vector, or 2-dimensional array.  The output
! vector will be the same type as X.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
pure function pdfns0(a) result (pdfn)
  real(kind=SP), intent(in) :: a
  real(kind=SP) :: pdfn
  real(kind=SP), parameter :: pi=3.14159265 ! 3589793238462643
  pdfn = exp((-a**2)/2.0)/((2.0*pi)**0.5)
end function pdfns0
pure function pdfnd0(a) result (pdfn)
  real(kind=DP), intent(in) :: a
  real(kind=DP) :: pdfn
  real(kind=DP), parameter :: pi=3.14159265358979323_dp ! 8462643d0
  pdfn = exp((-a**2)/2.0_dp)/((2.0_dp*pi)**0.5_dp)
end function pdfnd0
pure function pdfns1(a) result (pdfn)
  real(kind=SP), dimension(:), intent(in) :: a
  real(kind=SP), dimension(size(a)) :: pdfn
  real(kind=SP), parameter :: pi=3.14159265 ! 3589793238462643
  pdfn = exp((-a**2)/2.0)/((2.0*pi)**0.5)
end function pdfns1
pure function pdfnd1(a) result (pdfn)
  real(kind=DP), dimension(:), intent(in) :: a
  real(kind=DP), dimension(size(a)) :: pdfn
  real(kind=DP), parameter :: pi=3.14159265358979323_dp ! 8462643d0
  pdfn = exp((-a**2)/2.0_dp)/((2.0_dp*pi)**0.5_dp)
end function pdfnd1
pure function pdfns2(a) result (pdfn)
  real(kind=SP), dimension(:,:), intent(in) :: a
  real(kind=SP), dimension(size(a,1),size(a,2)) :: pdfn
  real(kind=SP), parameter :: pi=3.14159265 ! 3589793238462643
  pdfn = exp((-a**2)/2.0)/((2.0*pi)**0.5)
end function pdfns2
pure function pdfnd2(a) result (pdfn)
  real(kind=DP), dimension(:,:), intent(in) :: a
  real(kind=DP), dimension(size(a,1),size(a,2)) :: pdfn
  real(kind=DP), parameter :: pi=3.14159265358979323_dp ! 8462643d0
  pdfn = exp((-a**2)/2.0_dp)/((2.0_dp*pi)**0.5_dp)
end function pdfnd2


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! CDFN function
!
! Format: cdfn(x)
!
! Gives standard normal CDF of x.
! 
! X can be either single or double precision, and can be a scalar,
! vector, or 2-dimensional array.  The output
! vector will be the same type as X.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
function cdfn0(x) result(cdfn)
  real, intent(in) :: x
  real :: cdfn
  cdfn = real(cdfn0d(real(x,kind=DP)),kind=SP)
end function cdfn0
function cdfn1(x) result(cdfn)
  real, dimension(:), intent(in) :: x
  real, dimension(size(x)) :: cdfn
  cdfn=real(cdfn1d(real(x,kind=DP)),kind=SP)
end function cdfn1
function cdfn2(x) result(cdfn)
  real, dimension(:,:), intent(in) :: x
  real, dimension(Size(x,1),Size(x,2)) :: cdfn
  real, dimension(Size(x,1)*Size(x,2)) :: tmpx
  tmpx = reshape(x,(/ size(x,1)*size(x,2) /))
  tmpx = cdfn1(tmpx)
  cdfn = reshape(tmpx,(/ size(x,1), size(x,2) /))
end function cdfn2
pure function cdfn0d(x) result(cdfn)
  real(kind=DP), intent(in) :: x
  real(kind=DP) :: cdfn
  real(kind=DP), parameter :: a1 = 0.398942280444_dp, a2 = 0.399903438504_dp, &
       a3 = 5.75885480458_dp, a4 = 29.8213557808_dp, a5 = 2.62433121679_dp, &
       a6 = 48.6959930692_dp, a7 = 5.92885724438_dp, b0 = 0.398942280385_dp, &
       b1 = 3.8052e-08_dp, b2 = 1.00000615302_dp, b3 = 3.98064794e-04_dp, &
       b4 = 1.98615381364_dp, b5 = 0.151679116635_dp, b6 = 5.29330324926_dp, &
       b7 = 4.8385912808_dp, b8 = 15.1508972451_dp, b9 = 0.742380924027_dp, &
       b10 = 30.789933034_dp, b11 = 3.99019417011_dp
  real(kind=DP) ::  y
  logical :: smallx
  y = 0.5_dp * x**2
  cdfn=0.0_dp
  smallx = (abs(x) < 1.28_dp)
  if (smallx) cdfn = 0.5_dp - abs(x)  * ( a1 - a2 * y / ( y + a3 - a4 / ( y + a5 &
      + a6 / ( y + a7 ) ) ) )
  if ((.not.smallx).and.(abs(x) < 12.7_dp)) cdfn = exp ( - y ) * b0 / ( abs(x)  - b1 &
      + b2 / ( abs(x)  + b3 &
      + b4 / ( abs(x)  - b5 &
      + b6 / ( abs(x)  + b7 &
      - b8 / ( abs(x)  + b9 &
      + b10 / ( abs(x)  + b11 ) ) ) ) ) )
  if ( x > 0.0_dp ) cdfn = 1.0_dp-cdfn
end function cdfn0d
!function cdfn0d(x)
!  real(kind=DP), intent(in) :: x
!  real(kind=DP) :: cdfn0d
!  real(kind=DP) :: xtmp(1)
!  xtmp(1)=x
!  xtmp=cdfn1d(xtmp)
!  cdfn0d=xtmp(1)
!end function cdfn0d
pure function cdfn1d(x) result(cdfn)
  real(kind=DP), dimension(:), intent(in) :: x
  real(kind=DP), dimension(size(x)) :: cdfn
  real(kind=DP), parameter :: a1 = 0.398942280444_dp, a2 = 0.399903438504_dp, &
       a3 = 5.75885480458_dp, a4 = 29.8213557808_dp, a5 = 2.62433121679_dp, &
       a6 = 48.6959930692_dp, a7 = 5.92885724438_dp, b0 = 0.398942280385_dp, &
       b1 = 3.8052e-08_dp, b2 = 1.00000615302_dp, b3 = 3.98064794e-04_dp, &
       b4 = 1.98615381364_dp, b5 = 0.151679116635_dp, b6 = 5.29330324926_dp, &
       b7 = 4.8385912808_dp, b8 = 15.1508972451_dp, b9 = 0.742380924027_dp, &
       b10 = 30.789933034_dp, b11 = 3.99019417011_dp
  real(kind=DP), dimension(size(x)) ::  y
  logical, dimension(size(x)) :: smallx
  y = 0.5_dp * x**2
  cdfn=0.0_dp
  smallx = (abs(x) < 1.28_dp)
  where (smallx) cdfn = 0.5_dp - abs(x)  * ( a1 - a2 * y / ( y + a3 - a4 / ( y + a5 &
      + a6 / ( y + a7 ) ) ) )
  where ((.not.smallx).and.(abs(x) < 12.7_dp)) cdfn = exp ( - y ) * b0 / ( abs(x)  - b1 &
      + b2 / ( abs(x)  + b3 &
      + b4 / ( abs(x)  - b5 &
      + b6 / ( abs(x)  + b7 &
      - b8 / ( abs(x)  + b9 &
      + b10 / ( abs(x)  + b11 ) ) ) ) ) )
  where ( x > 0.0_dp ) cdfn = 1.0_dp-cdfn
end function cdfn1d
function cdfn2d(x) result(cdfn)
  real(kind=DP), dimension(:,:), intent(in) :: x
  real(kind=DP), dimension(Size(x,1),Size(x,2)) :: cdfn
  real(kind=DP), dimension(Size(x,1)*Size(x,2)) :: tmpx
  tmpx = reshape(x,(/ size(x,1)*size(x,2) /))
  tmpx = cdfn1d(tmpx)
  cdfn = reshape(tmpx,(/ size(x,1), size(x,2) /))
end function cdfn2d

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! CDFINVN function
!
! Format: cdfinvn(x)
!
! Gives inverse of standard normal CDF of x.
! 
! X can be either single or double precision, and can be a scalar,
! vector, or 2-dimensional array.  The output
! vector will be the same type as X.  In order for the inverse 
! CDF to exist, X must be strictly between zero and one.  This
! function does not check!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
function cdfinvn0(p) result(cdfinvn)
  real, intent(in) :: p
  real :: cdfinvn
  real, dimension(1) :: tmp
  tmp(1) = p
  tmp = cdfinvn1(tmp)
  cdfinvn=tmp(1)
end function cdfinvn0
function cdfinvn1(p) result(cdfinvn)
  real, dimension(:), intent(in) :: p
  real, dimension(Size(p)) :: cdfinvn
  logical, dimension(Size(p)) :: maskgt
  real, dimension(Size(p)) :: y,xp
  real, parameter :: lim=1.0e20, p0=-0.32223243, p1=-1.0, &
                     p2=-0.34224208, p3=-0.02042312, &
                     p4=-0.45364221e-4, q0=0.09934846, &
                     q1=0.58858157, q2=0.53110346, &
                     q3=0.10353775, q4=0.38560700e-2
  maskgt = (p > 0.5)
  xp = p
  where (maskgt) xp = 1-xp
  y = sqrt(-2.0*log(xp))
  xp = y + ((((p4*y + p3) * y + p2)* y + p1)*y + p0)/ &
       ((((q4*y + q3) * y + q2) * y + q1) * y + q0)
  where (maskgt) xp = -xp
  where (p == 0.5) xp = 0.0
  cdfinvn = -xp
end function cdfinvn1
function cdfinvn2(p) result(cdfinvn)
  real, dimension(:,:), intent(in) :: p
  real, dimension(Size(p,1),Size(p,2)) :: cdfinvn
  real, dimension(Size(p,1)*Size(p,2)) :: tmp
  tmp = reshape(p,(/ size(p,1)*size(p,2) /))
  tmp = cdfinvn1(tmp)
  cdfinvn = reshape(tmp,(/ size(p,1), size(p,2) /))
end function cdfinvn2
pure function cdfinvn0d(p)  result(cdfinvn)
  real(kind=DP), intent(in) :: p
  real(kind=DP) :: cdfinvn
  logical :: maskgt
  real(kind=DP) :: y
  real(kind=DP), parameter :: p0=-0.322232431088_dp, p1=-1.0_dp, & 
                     p2=-0.342242088547_dp, p3=-0.0204231210245_dp, &
                     p4=-0.453642210148e-4_dp, q0=0.0993484626060_dp, &
                     q1=0.588581570495_dp, q2=0.531103462366_dp, &
                     q3=0.103537752850_dp, q4=0.38560700634e-2_dp, &
                     minp=1.0e-300_dp
  maskgt = (p > 0.5_dp)
  cdfinvn = p
  if (maskgt) cdfinvn = 1.0_dp-cdfinvn
  if (cdfinvn < minp) then
     cdfinvn = minp
  end if
  y = sqrt(-2.0_dp*log(cdfinvn))
  cdfinvn = y + ((((p4*y + p3) * y + p2)* y + p1)*y + p0)/ &
       ((((q4*y + q3) * y + q2) * y + q1) * y + q0)
  if (.not.maskgt) cdfinvn=-cdfinvn
  if (p == 0.5_dp) cdfinvn = 0.0_dp
end function cdfinvn0d
pure function cdfinvn1d(p)  result(cdfinvn)
  real(kind=DP), dimension(:), intent(in) :: p
  real(kind=DP), dimension(size(p)) :: cdfinvn
  logical, dimension(size(p)) :: maskgt
  real(kind=DP), dimension(size(p)) :: y
  real(kind=DP), parameter :: p0=-0.322232431088_dp, p1=-1.0_dp, & 
                     p2=-0.342242088547_dp, p3=-0.0204231210245_dp, &
                     p4=-0.453642210148e-4_dp, q0=0.0993484626060_dp, &
                     q1=0.588581570495_dp, q2=0.531103462366_dp, &
                     q3=0.103537752850_dp, q4=0.38560700634e-2_dp, &
                     minp=1.0e-300_dp
  maskgt = (p > 0.5_dp)
  cdfinvn = p
  where (maskgt) cdfinvn = 1.0_dp-cdfinvn
  where (cdfinvn < minp) cdfinvn=minp
  y = sqrt(-2.0_dp*log(cdfinvn))
  cdfinvn = y + ((((p4*y + p3) * y + p2)* y + p1)*y + p0)/ &
       ((((q4*y + q3) * y + q2) * y + q1) * y + q0)
  where (.not.maskgt) cdfinvn=-cdfinvn
  where (p == 0.5_dp) cdfinvn = 0.0_dp
end function cdfinvn1d
function cdfinvn2d(p) result(cdfinvn)
  real(kind=DP), dimension(:,:), intent(in) :: p
  real(kind=DP), dimension(Size(p,1),Size(p,2)) :: cdfinvn
  real(kind=DP), dimension(Size(p,1)*Size(p,2)) :: tmp
  tmp = reshape(p,(/ size(p,1)*size(p,2) /))
  tmp = cdfinvn1d(tmp)
  cdfinvn = reshape(tmp,(/ size(p,1), size(p,2) /))
end function cdfinvn2d

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! CHOL function
!
! Format: chol(a)
!
! Gives Cholesky decomposition of the symmetric positive definite
! matrix A, i.e., the matrix L such that LL'=a.  If A is not
! a symmetric positive definite matrix, then this function
! will give back garbage.  A can be either single or double
! precision, with the result being of the same type as A.
! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
function chols(a) result(chol)
  real(kind=SP), dimension(:,:), intent(in) :: a
  real(kind=SP), dimension(size(a,1),size(a,2)) :: chol
  real(kind=SP), dimension(size(a,1)) :: p
  integer :: i,n
  real(kind=SP) :: sumn
  n = size(a,1) ! should check to make sure a is square
  chol = a
  do i=1,n
     sumn = chol(i,i) - dot_product(chol(1:i-1,i),chol(1:i-1,i))
     p(i) = sqrt(sumn)
     chol(i,i+1:n) = (chol(i+1:n,i)-matmul(chol(1:i-1,i),chol(1:i-1,i+1:n)))/p(i)
  end do
  do i = 1,n
     chol(i,i)=p(i)
     chol(i+1:n,i) = 0.0
  end do
end function chols
function chold(a) result(chol)
  real(kind=DP), dimension(:,:), intent(in) :: a
  real(kind=DP), dimension(size(a,1),size(a,1)) :: chol
  real(kind=DP), dimension(size(a,1)) :: p
  integer :: i,n
  real(kind=DP) :: sumn
  n = size(a,1) ! should check to make sure a is square
  chol = a
  do i=1,n
     sumn = chol(i,i) - dot_product(chol(1:i-1,i),chol(1:i-1,i))
     p(i) = sqrt(sumn)
     chol(i,i+1:n) = (chol(i+1:n,i)-matmul(chol(1:i-1,i),chol(1:i-1,i+1:n)))/p(i)
  end do
  do i = 1,n
     chol(i,i)=p(i)
     chol(i+1:n,i) = 0.0_dp
  end do
end function chold




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! OLS function
!
! Format: ols(x,y)
!
! Gives OLS coefficient vector for regression of x on y.
! 
! X and Y can be either single or double precision.  The output
! vector will be the same type as X and Y.  X should have a 
! column of ones if you want an intercept.
! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
function olss(x,y) result(olsout)
  real(kind=SP), dimension(:,:), intent(in) :: x
  real(kind=SP), dimension(:), intent(in) :: y
  real(kind=SP), dimension(size(x,2)) :: olsout,xprimey,p
  real(kind=SP), dimension(size(x,2),size(x,2)) :: xprimex
  real(kind=SP) :: summ
  integer :: i,n
  n=size(x,2)
  xprimex = matmul(transpose(x),x)
  xprimey = matmul(transpose(x),y)
  do i=1,n
     summ = xprimex(i,i)-dot_product(xprimex(i,1:i-1),xprimex(i,1:i-1))
     p(i) = sqrt(summ)
     xprimex(i+1:n,i) = (xprimex(i,i+1:n)-matmul(xprimex(i+1:n,1:i-1),xprimex(i,1:i-1)))/p(i)
  end do
  do i=1,n
     olsout(i)=(xprimey(i)-dot_product(xprimex(i,1:i-1),olsout(1:i-1)))/p(i)
  end do
  do i=n,1, -1
     olsout(i)=(olsout(i)-dot_product(xprimex(i+1:n,i),olsout(i+1:n)))/p(i)
  end do
end function olss
function olsd(x,y) result(ols)
  real(kind=DP), dimension(:,:), intent(in) :: x
  real(kind=DP), dimension(:), intent(in) :: y
  real(kind=DP), dimension(size(x,2)) :: ols,xprimey,p
  real(kind=DP), dimension(size(x,2),size(x,2)) :: xprimex
  real(kind=DP) :: summ
  integer :: i,n
  n=size(x,2)
  xprimex = matmul(transpose(x),x)
  xprimey = matmul(transpose(x),y)
  do i=1,n
     summ = xprimex(i,i)-dot_product(xprimex(i,1:i-1),xprimex(i,1:i-1))
     p(i) = sqrt(summ)
     xprimex(i+1:n,i) = (xprimex(i,i+1:n)-matmul(xprimex(i+1:n,1:i-1),xprimex(i,1:i-1)))/p(i)
  end do
  do i=1,n
     ols(i)=(xprimey(i)-dot_product(xprimex(i,1:i-1),ols(1:i-1)))/p(i)
  end do
  do i=n,1, -1
     ols(i)=(ols(i)-dot_product(xprimex(i+1:n,i),ols(i+1:n)))/p(i)
  end do
end function olsd

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! RANDINT subroutine (formerly a function)
!
! Format: randint(n)
!
! Gives a list of n random integers between 1 and n.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine randint(randintout)
  integer, dimension(:), intent(out) :: randintout
  real, dimension(size(randintout)) :: rnd_num
  integer :: n
  n=size(randintout)
  call runif(rnd_num)
  randintout = 1+floor(real(n,SP)*rnd_num)
end subroutine randint




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! FACTORIAL function
!
! Format: factorial(x)
!
! Given the integer (or vector of integers) x, returns 
!	x! or x*(x-1)*(x-2)*...*3*2*1
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
pure function factorial0(n) result(factorial)
  integer, intent(in) :: n
  integer :: factorial
  integer :: i
  factorial=1
  if (n > 1) then
     do i=2,n
        factorial=factorial*i
     end do
  end if
end function factorial0
pure function factorial1(n) result(factorial)
  integer, dimension(:), intent(in) :: n
  integer, dimension(size(n)) :: factorial
  integer :: i,maxn
  factorial=1
  maxn=maxval(n)
  if (maxn > 1) then
     do i=2,maxn
        where (n >= i) factorial=factorial*i
     end do
  end if
end function factorial1

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! NCHOOSEK function
!
! Format: nchoosek(n,k)
!
! Gives the number of k-length combinations of n items, or
!       n!
!     ----------
!     k!(n-k)!
! 
! N and K can both be scalars, or both be vectors, but
! must be integers.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
pure function nchoosek0(n,k) result(nchoosek)
  integer, intent(in) :: n,k
  integer :: nchoosek
  nchoosek = factorial(n)/(factorial(k)*factorial(n-k))
end function nchoosek0
pure function nchoosek1(n,k) result(nchoosek)
  integer, dimension(:), intent(in) :: n
  integer, dimension(size(n)), intent(in) :: k
  integer, dimension(size(n)) :: nchoosek
  nchoosek = factorial(n)/(factorial(k)*factorial(n-k))
end function nchoosek1


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! FIND_PRIMES function
!
! Format: find_primes(n)
!
! Returns the first n primes in a vector of length n, through 
! a very brute-force algorithm.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
function find_primes(n) result(primes)
  integer, intent(in) :: n
  integer, dimension(n) :: primes
  integer :: i,j,nprimes
  logical :: isaprime
  primes(1)=2
  nprimes=1
  do i=3,huge(1)
     isaprime = .true.
     do j=1,nprimes
        if (modulo(i,primes(j)) == 0) then
           isaprime = .false.
           exit
        end if
     end do
     if (isaprime) then
        nprimes=nprimes+1
        primes(nprimes) = i
        if (nprimes == n) return
     end if
  end do
end function find_primes

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! RHALT function
!
! Format: rhalt(r,c,shift)
!
! Gives a randomized Halton matrix with r rows and c columns.
! "shift" is the random number
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
function rhalt(r,c,shift) result(rhaltout)
  integer, intent(in) :: r,c
  real(kind=DP), intent(in) :: shift
  real(kind=DP), dimension(r,c) :: rhaltout
  rhaltout=hmat(r,c)+shift
  where (rhaltout > 1.0_dp) rhaltout=rhaltout-1.0_dp
end function rhalt

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! HMAT function
!
! Format: hmat(r,c)
!
! Gives a Halton matrix with r rows and c columns.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
function hmat(r,c) result(hmatout)
  integer, intent(in) :: r,c
  real(kind=DP), dimension(r,c) :: hmatout
  integer, dimension(c) :: prim
  integer :: droppit,i
  real(kind=DP), dimension(:), allocatable :: tmp
  prim=find_primes(c)
  droppit = max(10,prim(c))
  allocate(tmp(r+droppit))
  do i=1,c
     tmp = halton_sequence(r+droppit,prim(i))
     hmatout(:,i) = tmp((droppit+1):(droppit+r))
  end do
  deallocate(tmp)
end function hmat
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! HALTON_SEQUENCE function
!
! Format: halton_sequence(n,s)
!
! Gives a Halton sequence
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
function halton_sequence(n,s) result(hseqout)
  integer, intent(in) :: n,s
  real(kind=DP), dimension(n) :: hseqout
  integer :: i,j,k,len
  real(kind=DP), dimension(:), allocatable :: x
  k = ceiling(log(real(n+1))/log(real(s)))
  allocate(x(1+s**k))
  x(1) = 0.0_dp
  len=1
  do i=1,k
     do j=1,(s-1)
        x((j*len+1):((j+1)*len))=x(1:len)+real(j,DP)/real(s**i,DP)
     end do
     len=len*s
  end do
  hseqout  = x(2:(n+1))
  deallocate(x)
end function halton_sequence


end module bkmath
