!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! SMLE PROGRAM 
! Author: Brian Krauth, Simon Fraser University
!
! Estimates model developed in my paper "Simulation-based
! estimation of peer effects", and used in my papers
! "Peer effects and selection effects in youth smoking"
! and others.
!
! The program reads in the text file PARM.DAT, which contains
! all of the user-selected options.  It then reads in data
! from a file, estimates the model, and reports (appends) the results in 
! another file.  Intermediate results are stored in a log file.
! The exact names of the data file, result file, and log file
! are specified in PARM.DAT.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
program smle
  use bklib, only : dp,setseed  ! machine-specific code, e.g. code which uses LAPACK functions
  use smglob, only : starttime,smglobal,fstar,bstar ! global variables
  use smutil, only : load_data,estimate,cleanup ! other functions and subroutines
! SMC  use monte, only : montecarlo
  implicit none
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  VARIABLE DECLARATIONS
!
!  Most global variables are defined in the module smglob.
!  The program runs primarily in double precision.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  integer :: iter
  integer, parameter :: maxiter=1000
 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! START SMLE EXECUTABLE CODE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! SET UP
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  call date_and_time(time=STARTTIME) ! get current time
  call setseed() ! reset random number generator - machine specific 
  do iter=1,maxiter
! SMC call montecarlo()
     call load_data() ! loads all data from files 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! MAIN ESTIMATION CODE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     if ((SMGLOBAL%EQUILIBRIUM_TYPE == "B").or.(SMGLOBAL%EQUILIBRIUM_TYPE == "M")) then
        call likelihood_bounds(5)
     elseif (SMGLOBAL%EQUILIBRIUM_TYPE == "P") then
        call likelihood_plot(4.0_dp,10)
     else ! the usual case
        if (SMGLOBAL%RHO_TYPE == "I") then     ! User has requested interval estimates 
           call interval_estimate(0.9_dp,12)
        else
           call estimate()                        ! ESTIMATE is where the model is estimated
           call report_results()
        end if
     end if
  

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! CLEAN UP AND EXIT
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  call cleanup(STARTTIME)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! FUNCTIONS AND SUBROUTINES, IF ANY
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! SMC cycle  ! if monte carlo, repeat
  if (.not.SMGLOBAL%BOOTSTRAP) exit ! if not monte carlo, repeat only if bootstrapping
end do


contains

  subroutine report_results()
    logical :: file_found
    integer :: ios,pos
    real(kind=DP) :: loglikelihood,rhox,rhoe,gam
    real(kind=DP), dimension(SMGLOBAL%NVAR+1) :: beta
    loglikelihood = -FSTAR
    rhox = BSTAR(1)
    pos=2
    if ((SMGLOBAL%RHO_TYPE == "E").or.((SMGLOBAL%RHO_TYPE == "I").and.(SMGLOBAL%RUN_NUMBER == 2))) then
       rhoe = BSTAR(pos)
       pos = pos + 1
    elseif ((SMGLOBAL%RHO_TYPE == "X").or.((SMGLOBAL%RHO_TYPE=="I").and.(SMGLOBAL%RUN_NUMBER==1))) then
       rhoe=rhox
    else
       rhoe=SMGLOBAL%FIXED_RHO
    end if
    if ((SMGLOBAL%FIX_GAMMA).or.((SMGLOBAL%RHO_TYPE=="I").and.(SMGLOBAL%RUN_NUMBER==2))) then
       gam=SMGLOBAL%FIXED_GAMMA
    else
       gam = BSTAR(pos)
       pos = pos + 1
    end if
    beta = bstar(pos:size(bstar))
    inquire(file=SMGLOBAL%RESULTFILE,exist=file_found) ! see if resultfile already exists
    if (file_found) then 
       open(unit=1,file=SMGLOBAL%RESULTFILE,iostat=ios,action="write",position="append",status="old") ! append results if it exists
    else
       open(unit=1,file=SMGLOBAL%RESULTFILE,iostat=ios,action="write",status="new") ! if it doesn't exist, create a new one
    end if
    if (ios == 0) then
!       if (SMGLOBAL%RHO_TYPE == "I") then
!          if (SMGLOBAL%RUN_NUMBER == 1) then
!             write (unit=1,iostat=ios,fmt=*) -FSTAR, BSTAR(1), BSTAR 
!          elseif (SMGLOBAL%RUN_NUMBER == 2) then
!             write (unit=1,iostat=ios,fmt=*) -FSTAR, BSTAR(1:2), SMGLOBAL%FIXED_GAMMA, BSTAR(3:size(BSTAR))
!          else
!             write (unit=1,iostat=ios,fmt=*) -FSTAR, BSTAR(1), SMGLOBAL%FIXED_RHO, BSTAR(2:size(BSTAR))
!          end if
!       else
!          write (unit=1,iostat=ios,fmt=*) -FSTAR, BSTAR 
!       end if
       write (unit=1,iostat=ios,fmt=*) loglikelihood, rhox, rhoe, gam, beta
       close(unit=1,iostat=ios)
    else
       open(unit=1,file=SMGLOBAL%LOGFILE,iostat=ios,action="write",position="append",status="old")
       write (unit=1,iostat=ios,fmt=*) "WARNING: Unable to write results to file ", SMGLOBAL%RESULTFILE
       close(unit=1,iostat=ios)
    end if
    open(unit=1,file=SMGLOBAL%LOGFILE,iostat=ios,action="write",position="append",status="old")
    if (ios == 0) then
       write (unit=1,iostat=ios,fmt=*) "b is: ",BSTAR 
       write (unit=1,iostat=ios,fmt=*) "loglikelihood is: ",-FSTAR 
       close(unit=1,iostat=ios)
    else
       continue
    end if
  end subroutine report_results

  subroutine interval_estimate(rho_max,nsteps)
    real(kind=DP), intent(in) :: rho_max
    integer, intent(in) :: nsteps
    real(kind=DP) :: step
    integer :: i,istart,ios
    step=rho_max/real(nsteps-3,kind=DP)
    istart=SMGLOBAL%RUN_NUMBER
    SMGLOBAL%FIXED_GAMMA=0.001_dp
    SMGLOBAL%FIXED_RHO=0.0_dp
    do i=istart,nsteps
       SMGLOBAL%RUN_NUMBER=i
       if (i==1) then
          continue
       elseif (i==2) then
          SMGLOBAL%FIXED_GAMMA=0.001_dp
          open(unit=1,file=SMGLOBAL%LOGFILE,iostat=ios,action="write",position="append",status="old")
          if (ios == 0) then
             write (unit=1,iostat=ios,fmt=*) "Interval estimation: gamma=0.0"
             close(unit=1,iostat=ios)
          else
             continue
          end if
       else 
          open(unit=1,file=SMGLOBAL%LOGFILE,iostat=ios,action="write",position="append",status="old")
          if (ios == 0) then
             write (unit=1,iostat=ios,fmt=*) "Interval estimation: rho=",SMGLOBAL%FIXED_RHO
             close(unit=1,iostat=ios)
          else
             continue
          end if
       end if
       call estimate()
       call report_results()
       if (i > 2) then
          SMGLOBAL%FIXED_RHO=SMGLOBAL%FIXED_RHO+step
       end if
    end do
  end subroutine interval_estimate

  subroutine likelihood_bounds(nsteps)
    integer, intent(in) :: nsteps
    integer :: i,k,istart,istop,ios
    logical :: file_found
    istart=SMGLOBAL%RUN_NUMBER
    if (istart > 1) then
       k=size(BSTAR)
       deallocate(BSTAR)
       allocate(BSTAR(k-1))
    end if
    if (SMGLOBAL%EQUILIBRIUM_TYPE == "B") then
       istop=2+2*nsteps
    else
       istop=1+nsteps
    end if
    do i=istart,istop
       SMGLOBAL%RUN_NUMBER=i
       if (i==1) then 
          open(unit=1,file=SMGLOBAL%LOGFILE,iostat=ios,action="write",position="append",status="old")
          if (ios == 0) then
             write (unit=1,iostat=ios,fmt=*) "Finding lower bound for maximum of likelihood function"
             close(unit=1,iostat=ios)
          else
             continue
          end if
          call estimate()
          SMGLOBAL%DFPSTOP=FSTAR
          open(unit=1,file=SMGLOBAL%LOGFILE,iostat=ios,action="write",position="append",status="old")
          if (ios == 0) then
             write (unit=1,iostat=ios,fmt=*) "Lower bound for maximum of likelihood function = ", -SMGLOBAL%DFPSTOP
             close(unit=1,iostat=ios)
          else
             continue
          end if
          if (SMGLOBAL%RHO_TYPE == "E") then
             SMGLOBAL%GAMMA_START=BSTAR(3)
          else
             SMGLOBAL%GAMMA_START=BSTAR(2)
          end if
          SMGLOBAL%GAMMA_MAX=SMGLOBAL%GAMMA_START
          SMGLOBAL%GAMMA_MIN=0.0_dp
          SMGLOBAL%FIX_GAMMA=.true.
       elseif (i < nsteps+2) then
          SMGLOBAL%FIXED_GAMMA=(SMGLOBAL%GAMMA_MIN+SMGLOBAL%GAMMA_MAX)/2.0_dp       
          open(unit=1,file=SMGLOBAL%LOGFILE,iostat=ios,action="write",position="append",status="old")
          if (ios == 0) then
             write (unit=1,iostat=ios,fmt=*) "Searching for lower bound on gamma, currently considering gamma=",SMGLOBAL%FIXED_GAMMA
             close(unit=1,iostat=ios)
          else
             continue
          end if
          call estimate()
          if (FSTAR < SMGLOBAL%DFPSTOP) then
             SMGLOBAL%GAMMA_MAX=SMGLOBAL%FIXED_GAMMA
          else
             SMGLOBAL%GAMMA_MIN=SMGLOBAL%FIXED_GAMMA
          end if
          if (i== nsteps+1) then
             open(unit=1,file=SMGLOBAL%LOGFILE,iostat=ios,action="write",position="append",status="old")
             if (ios == 0) then
                write (unit=1,iostat=ios,fmt=*) "Lower bound found"
                close(unit=1,iostat=ios)
             else
                continue
             end if
             inquire(file=SMGLOBAL%RESULTFILE,exist=file_found)
             if (file_found) then              
                open(unit=1,file=SMGLOBAL%RESULTFILE,iostat=ios,action="write",position="append",status="old")
                if (ios /= 0) then
                   open(unit=1,file=SMGLOBAL%LOGFILE,iostat=ios,action="write",position="append",status="old")
                   if (ios == 0) then
                      write(unit=1,iostat=ios,fmt=*) "WARNING: Unable to write results to file ",SMGLOBAL%RESULTFILE
                   end if
                end if
             else
                open(unit=1,file=SMGLOBAL%RESULTFILE,iostat=ios,action="write",status="new") 
                if (ios /= 0) then
                   open(unit=1,file=SMGLOBAL%LOGFILE,iostat=ios,action="write",position="append",status="old")
                   if (ios == 0) then
                      write(unit=1,iostat=ios,fmt=*) "WARNING: Unable to write results to file ",SMGLOBAL%RESULTFILE
                   end if
                end if
             end if
             write (unit=1,iostat=ios,fmt=*) (SMGLOBAL%GAMMA_MAX+SMGLOBAL%GAMMA_MIN)/2.0_dp
             close (unit=1,iostat=ios)
          end if
       elseif (i == nsteps+2) then
          SMGLOBAL%GAMMA_MIN=SMGLOBAL%GAMMA_START
          SMGLOBAL%GAMMA_MAX=4.0_dp
          SMGLOBAL%FIXED_GAMMA=4.0_dp
          open(unit=1,file=SMGLOBAL%LOGFILE,iostat=ios,action="write",position="append",status="old")
          if (ios == 0) then
             write (unit=1,iostat=ios,fmt=*) "Searching for upper bound on gamma"
             close(unit=1,iostat=ios)
          else
             continue
          end if
          call estimate()
          if (FSTAR < SMGLOBAL%DFPSTOP) then
             open(unit=1,file=SMGLOBAL%LOGFILE,iostat=ios,action="write",position="append",status="old")
             if (ios == 0) then
                write (unit=1,iostat=ios,fmt=*) "Upper bound for gamma is > 4.0, stopping calculation there"
                close(unit=1,iostat=ios)
             else
                continue
             end if
             open(unit=1,file=SMGLOBAL%RESULTFILE,iostat=ios,action="write",position="append",status="old") 
             if (ios == 0) then
                write (unit=1,fmt=*) 4.0_dp
                close (unit=1)
             else
                open(unit=1,file=SMGLOBAL%LOGFILE,iostat=ios,action="write",position="append",status="old")
                if (ios == 0) then
                   write (unit=1,iostat=ios,fmt=*) "WARNING: Unable to write results to file ",SMGLOBAL%RESULTFILE
                   close (unit=1,iostat=ios)
                else
                   continue
                end if
             end if
             return
          end if
       else
          SMGLOBAL%FIXED_GAMMA=(SMGLOBAL%GAMMA_MIN+SMGLOBAL%GAMMA_MAX)/2.0_dp
          open(unit=1,file=SMGLOBAL%LOGFILE,iostat=ios,action="write",position="append",status="old")
          if (ios == 0) then
             write (unit=1,iostat=ios,fmt=*) "Searching for upper bound on gamma, currently considering gamma=",SMGLOBAL%FIXED_GAMMA
             close(unit=1,iostat=ios)
          else
             continue
          end if
          call estimate()
          if (FSTAR < SMGLOBAL%DFPSTOP) then
             SMGLOBAL%GAMMA_MIN=SMGLOBAL%FIXED_GAMMA
          else
             SMGLOBAL%GAMMA_MAX=SMGLOBAL%FIXED_GAMMA
          endif
          if (i == (2+2*nsteps)) then
             open(unit=1,file=SMGLOBAL%LOGFILE,iostat=ios,action="write",position="append",status="old")
             if (ios == 0) then
                write (unit=1,iostat=ios,fmt=*) "Upper bound found"
                close(unit=1,iostat=ios)
             else
                continue
             end if
             open(unit=1,file=SMGLOBAL%RESULTFILE,iostat=ios,action="write",position="append",status="old") 
             if (ios == 0) then
                write (unit=1,iostat=ios,fmt=*) (SMGLOBAL%GAMMA_MAX+SMGLOBAL%GAMMA_MIN)/2.0_dp
                close (unit=1,iostat=ios)
             else
                open(unit=1,file=SMGLOBAL%LOGFILE,iostat=ios,action="write",position="append",status="old")
                if (ios == 0) then
                   write (unit=1,iostat=ios,fmt=*) "WARNING: Unable to write results to file ",SMGLOBAL%RESULTFILE
                   close (unit=1,iostat=ios)
                else
                   continue
                end if
             end if
          end if
       end if
       open(unit=1,file=SMGLOBAL%LOGFILE,iostat=ios,action="write",position="append",status="old")
       if (ios == 0) then
          write (unit=1,iostat=ios,fmt=*) "b is: ",BSTAR 
          if (i == 1) then
             k=size(BSTAR)
             deallocate(BSTAR)
             allocate(BSTAR(k-1))
          else
             write (unit=1,iostat=ios,fmt=*) "Gamma is: ", SMGLOBAL%FIXED_GAMMA
          end if
          write (unit=1,iostat=ios,fmt=*) "loglikelihood is: ",-FSTAR 
          close (unit=1,iostat=ios)
       else
          stop "Error: Unable to write results to log file"
       end if
    end do
  end subroutine likelihood_bounds




  subroutine likelihood_plot(gamma_max,nsteps)
    real(kind=DP), intent(in) :: gamma_max
    integer, intent(in) :: nsteps
    real(kind=DP) :: step
    integer :: i,istart,istop,ios
    logical :: file_found
    step = gamma_max/real(nsteps-1,kind=DP)
    SMGLOBAL%FIX_GAMMA=.true.
    SMGLOBAL%FIXED_GAMMA=0.0_dp
    istart=SMGLOBAL%RUN_NUMBER
    istop=2*nsteps
    do i=istart,istop
       SMGLOBAL%RUN_NUMBER=i
       call estimate()
       inquire(file=SMGLOBAL%RESULTFILE,exist=file_found) ! see if resultfile already exists
       if (file_found) then 
          open(unit=1,file=SMGLOBAL%RESULTFILE,iostat=ios,action="write",position="append",status="old") ! append results if it exists
       else
          open(unit=1,file=SMGLOBAL%RESULTFILE,iostat=ios,action="write",status="new") ! if it doesn't exist, create a new one
       end if
       if (ios == 0) then
          if (modulo(i,2) > 0) then
             write (unit=1,iostat=ios,fmt=*) SMGLOBAL%FIXED_GAMMA, -FSTAR
          else
             write (unit=1,iostat=ios,fmt=*) -FSTAR          
             SMGLOBAL%FIXED_GAMMA = SMGLOBAL%FIXED_GAMMA+step
          end if
          close (unit=1,iostat=ios)
       else
          stop "Error: Unable to write results to result file"
       end if
    end do
  end subroutine likelihood_plot
  



end program smle










