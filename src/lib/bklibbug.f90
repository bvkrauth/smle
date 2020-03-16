!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! BKLIB Module
! Author: Brian Krauth, Simon Fraser University
!
! Machine-specific code for SMLE program.  This particular
! file is for SFU's Bugaboo machine.
!
! Public procedures defined in this module are:
!
!   runif,runifd:	Draw random numbers from uniform distribution
!			runif is single precision, runifd is double.
!			runif() gives a single number, runif(n) gives
!			a vector and runif(r,c) give a r-by-c matrix.
!   setseed: 		Reseed random number generator based on system 
!			clock.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module bklib ! contains any machine-specific code
implicit none
private
private :: runif0s,runif1s,runif2s,runif0d,runif1d,runif2d, &
     assert_eq2,assert_eq3,assert_eq4
public :: runif,setseed,sp,dp,assert_eq
integer, parameter, public :: SP=KIND(1.0),DP=selected_real_kind(9,99)
interface runif
   module procedure runif0s,runif1s,runif2s,runif0d,runif1d,runif2d
end interface
interface assert_eq
   module procedure assert_eq2,assert_eq3,assert_eq4
end interface
external g05faf,g05ccf  ! These functions are in the NAG (?) library

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! ASSERT_EQ function
! 
! Takes up to 4 arguments, stops program if they are not equal.
! Returns the (common) value if they are.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
function assert_eq4(a,b,c,d) result (n)
  integer, intent(in) :: a,b,c,d
  integer, intent(out) :: n
  n=a
  if (a /= b) stop "Error in assert_eq"
  if (a /= c) stop "Error in assert_eq"
  if (a /= d) stop "Error in assert_eq"
end function assert_eq4
function assert_eq3(a,b,c) result (n)
  integer, intent(in) :: a,b,c
  integer, intent(out) :: n
  n=a
  if (a /= b) stop "Error in assert_eq"
  if (a /= c) stop "Error in assert_eq"
end function assert_eq3
function assert_eq2(a,b) result (n)
  integer, intent(in) :: a,b
  integer, intent(out) :: n
  n=a
  if (a /= b) stop "Error in assert_eq"
end function assert_eq2


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! SETSEED subroutine
! 
! Reseeds random number generator based on system clock.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine setseed()
  call g05ccf 
end subroutine setseed


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! RUNIF, RUNIFD functions
! 
! Draw pseudorandom number(s) from uniform(0,1) distribution.
! 
! RUNIF gives single-precision numbers, RUNIFD gives double.
! 
! Format: RUNIF() gives a single number, RUNIF(n) gives a
! vector of length n, and RUNIF(r,c) gives an r-by-c matrix.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine runif0s(runifout)
  real(kind=SP), intent(inout) :: runifout
  real(kind=DP), dimension(1) :: x
  call g05faf(0.0_dp,1.0_dp,1,x)
  runifout=real(x(1),kind=SP)
end subroutine runif0s
subroutine runif1s(runifout)
  real(kind=SP), dimension(:), intent(inout) :: runifout
  real(kind=DP), dimension(size(runifout)) :: x
  call g05faf(0.0_dp,1.0_dp,size(runifout),x)
  runifout=real(x,kind=SP)
end subroutine runif1s
subroutine runif2s(runifout)
  real(kind=SP), dimension(:,:), intent(inout) :: runifout
  real(kind=DP), dimension(size(runifout,1)*size(runifout,2)) :: x
  integer :: r,c
  call g05faf(0.0_dp,1.0_dp,size(x),x)
  runifout=reshape(real(x,kind=SP),(/size(runifout,1),size(runifout,2)/))
end subroutine runif2s

subroutine runif0d(runifout)
  real(kind=DP), intent(inout) :: runifout
  real(kind=DP), dimension(1) :: x
  call g05faf(0.0_dp,1.0_dp,1,x)
  runifout=x(1)
end subroutine runif0d
subroutine runif1d(runifout)
  real(kind=DP), dimension(:), intent(inout) :: runifout
  call g05faf(0.0_dp,1.0_dp,size(runifout),runifout)
end subroutine runif1d
subroutine runif2d(runifout)
  real(kind=DP), dimension(:,:), intent(inout) :: runifout
  real(kind=DP), dimension(size(runifout,1)*size(runifout,2)) :: x
  integer :: r,c
  call g05faf(0.0_dp,1.0_dp,size(x),x)
  runifout=reshape(x,(/size(runifout,1),size(runifout,2)/))
end subroutine runif2d

end module bklib
