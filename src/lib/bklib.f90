!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! SMLIB Module
! Author: Brian Krauth, Simon Fraser University
!
! Machine-specific code for SMLE program.  By default, the
! program uses the standard random number generator associated
! with the random_number subroutine.  The standard random number
! generator is not always all that good; so this file can be altered
! to substitute a preferred algorithm.
!
! Public procedures defined in this module are:
!
!   runif:	        Draw random numbers from uniform distribution
!   setseed: 		Reseed random number generator based on system 
!			clock.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module bklib ! contains any machine-specific code
implicit none
public :: setseed, runif
private :: runif0s,runif1s,runif2s,runif0d,runif1d,runif2d
integer, parameter, public :: SP=KIND(1.0),DP=selected_real_kind(9,99)
interface runif
   module procedure runif0s,runif1s,runif2s,runif0d,runif1d,runif2d
end interface
! external g05faf,g05ccf ! These functions are in the NAG library, which is unavailable on this machine

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! SETSEED subroutine
! 
! Reseeds random number generator based on system clock.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine setseed()
  integer :: i,n
  real :: j
  call system_clock(n) 
  n = modulo(n,40000)
  do i=1,n
     call random_number(j)
  end do
!   call g05ccf 
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
subroutine runif0s(x) 
  real(kind=SP), intent(inout) :: x
  call random_number(x)
end subroutine runif0s
subroutine runif1s(x)
  real(kind=SP), dimension(:), intent(inout) :: x
  call random_number(x)
end subroutine runif1s
subroutine runif2s(x)
  real(kind=SP), dimension(:,:), intent(inout) :: x
  call random_number(x)
end subroutine runif2s
subroutine runif0d(x) 
  real(kind=DP), intent(inout) :: x
  call random_number(x)
end subroutine runif0d
subroutine runif1d(x)
  real(kind=DP), dimension(:), intent(inout) :: x
  call random_number(x)
end subroutine runif1d
subroutine runif2d(x)
  real(kind=DP), dimension(:,:), intent(inout) :: x
  call random_number(x)
end subroutine runif2d


end module bklib
