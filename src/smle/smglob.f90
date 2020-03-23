!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! SMGLOB MODULE
! Author: Brian Krauth, Simon Fraser University
!
! File of global variables for SMLE program.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module smglob 
use bklib, only : DP,strlen
implicit none
private
public :: save_dfp_checkpoint,save_sa_checkpoint,load_checkpoint, &
     load_dfp_checkpoint,load_sa_checkpoint

type, public :: smcontrols
 integer :: nobs,nvar,maxgroupsize,nsim,restarts,numagg,fixedeffects,istart,run_number
 logical :: underreporting_correction,bootstrap,fix_gamma
 real(kind=DP) :: fixed_rho,fixed_gamma,dfpstop,gamma_start,gamma_min,gamma_max
 character(len=strlen) :: logfile,resultfile
 character(len=1) :: equilibrium_type,search_method,rho_type,simulator_type
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
character(len=*), parameter, public:: checkpointfile="check.dat",lockfile="check.lock", &
     version_number="1.2.1"
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Data 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
real(kind=DP), public :: fstar=1.0e50_dp,dfp_fp
real(kind=DP), dimension(:), allocatable, public :: bstar
real(kind=DP), dimension(:,:), allocatable, public :: u,x,dfp_hessin
real(kind=DP), dimension(:), allocatable, public :: y,ygroup,reporting_rate,dfp_p,dfp_g 
integer, dimension(:), allocatable, public :: nfriends,yint,ygint
integer, public :: dfp_itstart
character(len=80), public :: starttime,endtime

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Bigmatsize contains the maximum number of equilibria (sort of)
! for a given number of agents.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
integer, parameter, dimension(11), public :: bigmatsize=(/0,2,4,8,15,28,51,92,164,290,509/)
real(kind=DP), parameter, public :: maxrho=0.99_dp,mingam=0.001_dp


contains


subroutine save_dfp_checkpoint(its,p,fp,g,hessin)
  integer, intent(in) :: its
  real(kind=DP), dimension(:), intent(in) :: p,g
  real(kind=DP), intent(in) :: fp
  real(kind=DP), dimension(:,:), intent(in) :: hessin
  integer :: ios
  open (unit=1,file=LOCKFILE,iostat=ios,action="write",form="unformatted",status="replace")
  if (ios == 0) then
     write (unit=1,iostat=ios) its
     close (unit=1,iostat=ios)
  else
     continue
  end if
  open (unit=1,file=CHECKPOINTFILE,iostat=ios,action="write",form="unformatted",status="replace")
  if (ios == 0) then
     write (unit=1,iostat=ios) SMGLOBAL
     write (unit=1,iostat=ios) FSTAR,BSTAR,U,X,Y,YGROUP,REPORTING_RATE,NFRIENDS,YINT,YGINT
     write (unit=1,iostat=ios) its,p,fp,g,hessin 
     close (unit=1,iostat=ios)
  else
     continue
  end if
end subroutine save_dfp_checkpoint

subroutine save_sa_checkpoint(t,p)
  real(kind=DP), dimension(:), intent(in) :: p
  real(kind=DP), intent(in) :: t
  integer :: ios
  open (unit=1,file=LOCKFILE,iostat=ios,action="write",form="unformatted",status="replace")
  if (ios == 0) then
     write (unit=1,iostat=ios) t
     close (unit=1,iostat=ios)
  else
     continue
  end if
  open (unit=1,file=CHECKPOINTFILE,iostat=ios,action="write",form="unformatted",status="replace")
  if (ios == 0) then
     write (unit=1,iostat=ios) SMGLOBAL
     write (unit=1,iostat=ios) FSTAR,BSTAR,U,X,Y,YGROUP,REPORTING_RATE,NFRIENDS,YINT,YGINT
     write (unit=1,iostat=ios) t,p
     close (unit=1,iostat=ios)
  else
     continue
  end if
end subroutine save_sa_checkpoint

subroutine load_checkpoint()
  integer :: k,nobs,ios
  logical :: file_found
  open (unit=1,file=checkpointfile,iostat=ios,action="read",form="unformatted",position="rewind",status="old")  
  if (ios == 0) then
     read (unit=1,iostat=ios) smglobal
     if (ios /= 0) then
        stop "Error: Checkpoint file is corrupted.  Please delete and start over."
     end if
     k=smglobal%nvar+3
     nobs=SMGLOBAL%NOBS
     allocate(ygroup(nobs),y(nobs),x(nobs,smglobal%nvar),nfriends(nobs), &
          yint(nobs),ygint(nobs),reporting_rate(nobs),bstar(k), &
          u(smglobal%maxgroupsize,smglobal%nsim),dfp_p(k),dfp_g(k),dfp_hessin(k,k))
     read (unit=1,iostat=ios) fstar,bstar,u,x,y,ygroup,reporting_rate,nfriends,yint,ygint
     if (ios /= 0) then
        stop "Error: Checkpoint file is corrupted.  Please delete and start over."
     end if
     if (SMGLOBAL%SEARCH_METHOD == "D") then
        read (unit=1,iostat=ios) dfp_itstart,dfp_p,dfp_fp,dfp_g,dfp_hessin
        if (ios /= 0) then
           stop "Error: Checkpoint file is corrupted.  Please delete and start over."
        end if
     elseif (SMGLOBAL%SEARCH_METHOD == "S") then
        read (unit=1,iostat=ios) dfp_fp,dfp_p
        if (ios /= 0) then
           stop "Error: Checkpoint file is corrupted.  Please delete and start over."
        end if
     else
        stop "Error: Checkpoint file is corrupted.  Please delete and start over."
     end if
     close(unit=1,iostat=ios)
  else
     stop "Error: There is a checkpoint file, but it cannot be opened."
  end if
  inquire(file=SMGLOBAL%LOGFILE,exist=file_found)
  if (file_found) then
     open(unit=1,file=SMGLOBAL%LOGFILE,iostat=ios,action="write",position="append",status="old") 
  else
     open(unit=1,file=SMGLOBAL%LOGFILE,iostat=ios,action="write",status="new") 
  end if
  if (ios == 0) then
     write (unit=1,iostat=ios,fmt=*) "SMLE ",version_number
     write (unit=1,iostat=ios,fmt=*) "Author: Brian Krauth, Simon Fraser University"
     write (unit=1,iostat=ios,fmt=*) "Time restarted from checkpoint:   ",starttime(1:2),":",starttime(3:4),":",starttime(5:10)
     close (unit=1,iostat=ios)
  else
     stop "Error: Unable to open log file"
  end if
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
  t = dfp_fp
  p = dfp_p
  deallocate(dfp_p,dfp_g,dfp_hessin)
end subroutine load_sa_checkpoint


end module smglob



