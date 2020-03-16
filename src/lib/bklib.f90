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
public :: setseed, runif, read_buffer, get_parmfile
private :: runif0s,runif1s,runif2s,runif0d,runif1d,runif2d
integer, parameter, public :: SP=KIND(1.0),DP=selected_real_kind(9,99), strlen=32767
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


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! READ_BUFFER subroutine
!
! Format: read_buffer(buffer)
!
! Reads a line of data into a character buffer. This subroutine
! was written to work around a problem in the IBM XLFortran
! compiler: when reading in a line of text from a DOS/Windows
! text file (i.e., one with lines ending in a CR-LF rather than a LF)
! the program would include the CR as a character rather than
! dropping it.  This leads to all sorts of problems, so we
! just read in each line into a big string buffer, and change
! any CR at the end of the record into a space.
!
! Note that this workaround is only necessary for the IBM 
! compiler; other compilers seem to strip the CR automatically.
! It also comes at the cost of imposing a fixed limit on
! the size of a record (right now I have it at 1000 characters).
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine read_buffer(unit,buffer)
    integer, intent(in) :: unit
    character(len=*), intent(inout) :: buffer
    integer :: ios
    ! Read the data into the text buffer
    read (unit=unit,iostat=ios,fmt="(a)") buffer 
    ! If the last character is anything but a blank, then we might
    ! have records that are too big.  
    if (len(buffer) == len_trim(buffer)) then 
       stop "Error: Records in data file are larger than maximum buffer size"
    else
       ! Here's the important part.  If the lat non-blank character is a CR, replace it with a space.
       if (iachar(buffer(len_trim(buffer):len_trim(buffer))) == 13) then
          buffer(len_trim(buffer):len_trim(buffer)) = " "
       end if
    end if 
  end subroutine read_buffer

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! GET_PARMFILE subroutine
!
! Format: get_parmfile(parmfile)
! 
! Retrieves the name of the input file if provided in the command
! arguments.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine get_parmfile(parmfile)
    character(len=*), intent(inout) :: parmfile
    character(len=strlen) :: buffer
    integer :: status
    status = command_argument_count()
    if (status > 0) then
      call get_command_argument(1,buffer)
      parmfile = trim(buffer)
    end if
  end subroutine get_parmfile
  	
end module bklib
