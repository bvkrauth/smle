!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! SMGLOB MODULE
! Author: Brian Krauth, Simon Fraser University
!
! File of global variables for SMLE program.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module smglob 
use bklib, only : DP
implicit none
private
public :: save_dfp_checkpoint, save_sa_checkpoint, load_checkpoint, load_dfp_checkpoint, &
     load_sa_checkpoint


type, public :: smcontrols
   integer :: nobs,nvar,ngroups,maxgroupsize,nsim,restarts,numagg,istart,run_number
   logical :: fix_gamma
   real(kind=DP) :: fixed_rho,fixed_gamma,dfpstop,gamma_start,gamma_min,gamma_max
   character(len=12) :: logfile,resultfile
   character(len=1) :: covmat_type,equilibrium_type,search_method,rho_type
end type smcontrols

type, public :: sacontrols
   integer :: ns,nt
   real(kind=DP) :: t,rt,eps
end type sacontrols


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Control parameters, set in PARM.DAT
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
type(smcontrols), public :: smglobal
type(sacontrols), public :: saglobal
logical, public :: resume_from_checkpoint=.false.
character(len=*), parameter, public :: checkpointfile="check.dat",lockfile="check.lock",version_number="1.1"
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Data 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
real(kind=DP), public :: fstar=1.0e50_dp,dfp_fp
real(kind=DP), dimension(:), allocatable, public :: bstar
real(kind=DP), dimension(:,:), allocatable, public :: u,x,dfp_hessin,z,covmat
real(kind=DP), dimension(:), allocatable, public :: y,ygroup,dfp_p,dfp_g,bx
integer, dimension(:), allocatable, public :: original_groupid,groupid,groupsize
integer, public :: dfp_itstart
character(len=80), public :: STARTTIME,ENDTIME
real(kind=DP), parameter, public :: maxrho=0.99_dp,mingam=0.001_dp

contains

subroutine save_dfp_checkpoint(its,p,fp,g,hessin)
  integer, intent(in) :: its
  real(kind=DP), dimension(:), intent(in) :: p,g
  real(kind=DP), intent(in) :: fp
  real(kind=DP), dimension(:,:), intent(in) :: hessin
  integer :: ios
  logical :: file_found
  inquire(file=LOCKFILE,exist=file_found) 
  if (file_found) then
     open (unit=1,file=LOCKFILE,iostat=ios,action="write",form="unformatted",position="rewind",status="replace")
  else
     open (unit=1,file=LOCKFILE,iostat=ios,action="write",form="unformatted",status="new")
  end if
  if (ios == 0) then
     write (unit=1,iostat=ios) its
     close (unit=1,iostat=ios)
  end if
  inquire(file=CHECKPOINTFILE,exist=file_found) 
  if (file_found) then
     open (unit=1,file=CHECKPOINTFILE,iostat=ios,action="write",form="unformatted",position="rewind",status="replace")
  else
     open (unit=1,file=CHECKPOINTFILE,iostat=ios,action="write",form="unformatted",status="new")
  end if
  if (ios == 0) then 
     write (unit=1,iostat=ios) SMGLOBAL
     write (unit=1,iostat=ios) FSTAR,BSTAR,U,X,Z,COVMAT,Y,YGROUP,BX,ORIGINAL_GROUPID,GROUPID,GROUPSIZE
     write (unit=1,iostat=ios) its,p,fp,g,hessin
     close (unit=1,iostat=ios)
  end if
end subroutine save_dfp_checkpoint

subroutine save_sa_checkpoint(t,p)
  real(kind=DP), dimension(:), intent(in) :: p
  real(kind=DP), intent(in) :: t
  integer :: ios
  logical :: file_found
  inquire(file=LOCKFILE,exist=file_found) 
  if (file_found) then
     open (unit=1,file=LOCKFILE,iostat=ios,action="write",form="unformatted",position="rewind",status="replace")
  else
     open (unit=1,file=LOCKFILE,iostat=ios,action="write",form="unformatted",status="new")
  end if
  if (ios == 0) then  
     write (unit=1,iostat=ios) t
     close (unit=1,iostat=ios)
  end if
  inquire(file=CHECKPOINTFILE,exist=file_found) 
  if (file_found) then
     open (unit=1,file=CHECKPOINTFILE,iostat=ios,action="write",form="unformatted",position="rewind",status="replace")
  else
     open (unit=1,file=CHECKPOINTFILE,iostat=ios,action="write",form="unformatted",status="new")
  end if
  if (ios == 0) then
     write (unit=1,iostat=ios) SMGLOBAL
     write (unit=1,iostat=ios) FSTAR,BSTAR,U,X,Z,COVMAT,Y,YGROUP,BX,ORIGINAL_GROUPID,GROUPID,GROUPSIZE
     write (unit=1,iostat=ios) t,p
     close (unit=1,iostat=ios)
  end if
end subroutine save_sa_checkpoint

subroutine load_checkpoint()
  integer :: k,nobs,ios
  open (unit=1,file=CHECKPOINTFILE,iostat=ios,action="read",form="unformatted",position="rewind",status="old")
  if (ios /= 0) then
     stop "Error: checkpoint file exists, but cannot be opened."
  end if
  read (unit=1,iostat=ios) SMGLOBAL  
  if (ios /= 0) then
     stop "Error: checkpoint file exists, but is corrupted."
  end if
  nobs=SMGLOBAL%NOBS
  k=SMGLOBAL%NVAR+3
  allocate(BSTAR(k),U(SMGLOBAL%MAXGROUPSIZE,SMGLOBAL%NSIM),X(nobs,SMGLOBAL%NVAR-SMGLOBAL%NUMAGG),&
       COVMAT(k,k),Y(nobs),YGROUP(nobs),BX(nobs), &
       ORIGINAL_GROUPID(nobs),GROUPID(nobs),GROUPSIZE(SMGLOBAL%NGROUPS), &
       DFP_P(k),DFP_G(k),DFP_HESSIN(k,k))
  if (SMGLOBAL%NUMAGG > 0) then
     allocate(z(smglobal%ngroups,smglobal%numagg))
  else 
     allocate(z(1,1))
  end if
  read (unit=1,iostat=ios) fstar,bstar,u,x,z,covmat,y,ygroup,bx,original_groupid,groupid,groupsize
  if (ios /= 0) then
     stop "Error: checkpoint file exists, but is corrupted."
  end if
  if (SMGLOBAL%SEARCH_METHOD == "D") then
     read (unit=1,iostat=ios) dfp_itstart,dfp_p,dfp_fp,dfp_g,dfp_hessin
     if (ios /= 0) then
        stop "Error: checkpoint file exists, but is corrupted."
     end if
  elseif (SMGLOBAL%SEARCH_METHOD == "S") then
     read (unit=1,iostat=ios) dfp_fp,dfp_p
     if (ios /= 0) then
        stop "Error: checkpoint file exists, but is corrupted."
     end if
  else
     stop "Error: checkpoint file exists, but is corrupted"
  end if
  close(unit=1,iostat=ios)
  RESUME_FROM_CHECKPOINT=.true.
end subroutine load_checkpoint

subroutine load_dfp_checkpoint(its,p,fp,g,hessin)
  integer, intent(out) :: its
  real(kind=DP), dimension(:), intent(out) :: p,g
  real(kind=DP), intent(out) :: fp
  real(kind=DP), dimension(:,:), intent(out) :: hessin
  its = dfp_itstart
  p = dfp_p
  fp = dfp_fp
  g = dfp_g
  hessin = dfp_hessin
  deallocate(dfp_p,dfp_g,dfp_hessin)
end subroutine load_dfp_checkpoint

subroutine load_sa_checkpoint(t,p)
  real(kind=DP), dimension(:), intent(out) :: p
  real(kind=DP), intent(out) :: t
  t=dfp_fp
  p = dfp_p
  deallocate(dfp_p,dfp_g,dfp_hessin)
end subroutine load_sa_checkpoint

end module smglob



