!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! SMLE2 PROGRAM 
! Author: Brian Krauth, Simon Fraser University
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module smle2
  use bklib, only : DP,runif  ! machine-specific code
  use bkmath, only : pdfn,cdfinvn,ols,rhalt,inverse  ! useful math functions
  use smglob ! global variables
  use loglik, only : loglikelihood, dloglikelihood, loglikevec ! log-likelihood
  use dfpmin, only : dfp
  use simann, only : sa
  implicit none
  private
  private :: uppercase,get_startb
  public :: estimate,cleanup,load_data

contains




  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! ESTIMATE subroutine
!
! Format: estimate()
!
! Estimate the model, returning the resulting vector of parameter
! estimates in BSTAR and the maximized log-likelihood in 
! FSTAR.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine estimate()
    real(kind=DP), dimension(size(BSTAR)) :: startb,b,c,lb,ub,shift
    real(kind=DP), dimension(size(X,2)) :: mu,sigmasq
!    real(kind=DP), dimension(SMGLOBAL%NGROUPS) :: fhat
    real(kind=DP), dimension(SMGLOBAL%NGROUPS) :: fhat,fhat2
    real(kind=DP), dimension(SMGLOBAL%NGROUPS,size(BSTAR)) :: ghat
    real(kind=DP) :: fret,n,ftmp
    real(kind=DP), parameter :: eps=0.01_dp
    integer :: istart,i,iter,j,nacc,nfcnev,nobds,ier,intercept,ios
!---------------------------------------------------------------
! Standardize x for estimation - this makes the numerical optimization
! algorithm (DFP) less likely to get lost
!---------------------------------------------------------------
    n=real(size(X,1),kind=DP)
    do i=1,size(mu)
       mu(i)=sum(X(:,i))/n
       sigmasq(i)=sum((X(:,i)-mu(i))**2)/n
       X(:,i)=(X(:,i)-mu(i))/sqrt(sigmasq(i))
    end do
    if (SMGLOBAL%NUMAGG > 0) then
       do i=1,SMGLOBAL%NUMAGG
          z(:,i)=(z(:,i)-mu(i))/sqrt(sigmasq(i))
       end do
    end if 
!---------------------------------------------------------------
! Get a good initial guess on the parameter vector.
!---------------------------------------------------------------
    startb = get_startb(Y,YGROUP,X) ! MUST BE UPDATED
    open(unit=1,file=SMGLOBAL%LOGFILE,iostat=ios,action="write",position="append",status="old")
    if (ios == 0) then
       write (unit=1,iostat=ios,fmt=*) "startb: ", startb
       close(unit=1,iostat=ios)
    end if
    b = startb
    if (RESUME_FROM_CHECKPOINT) then
       ! if this is the case, we are actually restarting a partially-completed run of the program
       ! based on information in a "checkpoint" file (CHECK.DAT) that has been saved
       ! in the current directory.  The last saved values of
       ! bstar and fstar have already been retrieved from the checkpoint file
       istart=SMGLOBAL%ISTART
    else
       BSTAR = startb ! start out with the initial guess
!       FSTAR = loglikelihood(BSTAR) ! figure out what that gives us...
       call loglikelihood(BSTAR,FSTAR) ! new 
       istart = 1
    end if
!---------------------------------------------------------------
! Then adjust the parameter vector to maximize the log-likelihood
!---------------------------------------------------------------
    if (SMGLOBAL%SEARCH_METHOD == "S") then
       c=2.0_dp
       lb=-2.0_dp
       ub=2.0_dp
! Limits on rho_x
       lb(1)=-1.0_dp/real(SMGLOBAL%MAXGROUPSIZE,kind=DP)
       ub(1)=0.99_dp
! Limits on rho_e, if estimated
       i=2
       if ((SMGLOBAL%RHO_TYPE=="E").or.((SMGLOBAL%RHO_TYPE=="I").and.(SMGLOBAL%RUN_NUMBER==2))) then
          lb(2)=-1.0_dp/real(SMGLOBAL%MAXGROUPSIZE,kind=DP)
          ub(2)=0.99_dp
          i=i+1
       end if
! Limits on gamma, if estimated
       if ((SMGLOBAL%FIX_GAMMA).or.((SMGLOBAL%RHO_TYPE=="I").and.(SMGLOBAL%RUN_NUMBER==2))) then
          continue  
       else
          lb(i)=-0.1_dp
          ub(i)=4.0_dp
          i=i+1
       endif
! Limits on intercept
       lb(i)=-3.0_dp
       ub(i)=3.0_dp
       nacc=0
       nfcnev=0
       nobds=0
       ier=0
       if (RESUME_FROM_CHECKPOINT) then
          call load_sa_checkpoint(SAGLOBAL%T,startb)
       end if
       call sa(startb,.false.,SAGLOBAL%RT,SAGLOBAL%EPS,SAGLOBAL%NS,SAGLOBAL%NT,4,100000,lb,ub,c,2, &
            SAGLOBAL%T,bstar,fstar,nacc,nfcnev,nobds,ier,loglikelihood)
    elseif (SMGLOBAL%SEARCH_METHOD == "D") then
       do i = istart,SMGLOBAL%RESTARTS ! Start from a few different points to avoid finding only a local optimum
          SMGLOBAL%ISTART=i ! save the current iteration number 
          fret=FSTAR ! this is just to make sure that fret is initialized
          call dfp(b,0.00001_dp,iter,fret,loglikelihood,dloglikelihood)  ! Main optimization procedure
          if (fret < FSTAR) then ! If this iteration of DFP has improved on our highest likelihood so far...
             FSTAR = fret ! ...then save that result
             BSTAR = b
          end if
          call runif(shift)
          b = 2.0_dp*shift*b ! Generate a random (but not crazy) starting point for restart.
       end do
    else
       open(unit=1,file=SMGLOBAL%LOGFILE,iostat=ios,action="write",position="append",status="old")
       if (ios == 0) then
          write (unit=1,iostat=ios,fmt=*) "Error in ESTIMATE: search_method not defined"
          close(unit=1,iostat=ios) 
       end if
       stop "Error in ESTIMATE: search_method not defined"
    end if
!---------------------------------------------------------------
! At this point BSTAR contains the ML coefficients for the 
! standardized X variables.  We need to transform the coefficients
! so that they apply to the original X variables.
!---------------------------------------------------------------
    intercept = size(BSTAR)-size(mu)
    do i=1,size(mu)
       BSTAR(i+intercept)=BSTAR(i+intercept)/sqrt(sigmasq(i)) ! unstandardize the coefficients   
       BSTAR(intercept)=BSTAR(intercept)-BSTAR(i+intercept)*mu(i)     ! then the intercept        
       X(:,i)=(sqrt(sigmasq(i))*X(:,i))+mu(i) ! then unstandardize the X variable
    end do
    if (SMGLOBAL%NUMAGG > 0) then
       do i=1,SMGLOBAL%NUMAGG
          Z(:,i)=sqrt(sigmasq(i))*Z(:,i)+mu(i)
       end do
    end if
!---------------------------------------------------------------
! Then estimate the covariance matrix by the inverse hessian 
!---------------------------------------------------------------
    if (SMGLOBAL%COVMAT_TYPE == "H") then
       do i=1,size(b)
          b=BSTAR
          b(i)=b(i)+eps
!          startb(i)=loglikelihood(b)
          call loglikelihood(b,ftmp) ! new
          startb(i)=ftmp ! new
          do j=i,size(b)
             b(j)=b(j)+eps
!             COVMAT(i,j)=loglikelihood(b)
             call loglikelihood(b,ftmp)
             COVMAT(i,j)=ftmp
             b(j)=b(j)-eps
          end do
       end do
       do i=1,size(b)
          do j=i,size(b)
             COVMAT(i,j)=COVMAT(i,j)+fret-startb(i)-startb(j)
             COVMAT(j,i)=COVMAT(i,j)
          end do
       end do
       COVMAT=COVMAT/(eps*eps)
       call inverse(COVMAT)
!---------------------------------------------------------------
! Or estimate the covariance matrix by the BHHH/OPG method 
!---------------------------------------------------------------
    elseif (SMGLOBAL%COVMAT_TYPE == "O") then
!       fhat=loglikevec(BSTAR)
       call loglikevec(BSTAR,fhat)
       do i=1,size(b)
          b=BSTAR
          b(i)=b(i)+eps
!          ghat(:,i)=(loglikevec(b)-fhat)/eps
          call loglikevec(b,fhat2) ! new
          ghat(:,i)=(fhat2-fhat)/eps ! new
       end do
       COVMAT=matmul(transpose(ghat),ghat)
       call inverse(COVMAT)
    end if
    SMGLOBAL%ISTART=1 ! reset in case ESTIMATE is called more than once in program
  end subroutine estimate

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! GET_STARTB function
!
! Format: get_startb(y,ygroup,x)
!
! Get a good starting guess for the parameter vector.  Our
! guess is proportional to the results of 
! a simple linear probability model.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
  function get_startb(y,ygroup,x) result(startb)
    real(kind=DP), dimension(:), intent(in) :: y, ygroup
    real(kind=DP), dimension(:,:), intent(in) :: x
    real(kind=DP), dimension(size(BSTAR)) :: startb 
    real(kind=DP), dimension(size(y)) :: ytmp
    real(kind=DP), dimension(size(x,1),size(X,2)+2) :: z
    integer :: i
    z(:,1) = ygroup
    z(:,2) = 1.0_dp
    z(:,3:Size(z,2)) = x
    ytmp=y-0.5_dp
    i=1+size(startb)-size(z,2)
    startb(i:size(startb)) = ols(z,ytmp)/pdfn(cdfinvn(sum(y)/real(size(y),DP)))
    startb(1) = 0.2_dp  ! rho_x
    i=2
    if (SMGLOBAL%RHO_TYPE=="X") then
       continue
    elseif (SMGLOBAL%RHO_TYPE=="E") then
       startb(i) = 0.2_dp
       i=i+1
    elseif (SMGLOBAL%RHO_TYPE=="F") then
       continue
    elseif (SMGLOBAL%RHO_TYPE=="I") then
       if (SMGLOBAL%RUN_NUMBER == 2) then
          startb(i) = 0.2_dp
          i=i+1
       end if
    end if
    if ((SMGLOBAL%FIX_GAMMA).or.((SMGLOBAL%RHO_TYPE=="I").and.(SMGLOBAL%RUN_NUMBER == 2))) then
       continue
    else
       startb(i) = startb(i)/2.0_dp 
    end if
  end function get_startb


  subroutine cleanup()
    logical :: file_found
    character(len=80) :: endtime
    integer :: ios
    open(unit=1,file=SMGLOBAL%LOGFILE,iostat=ios,action="write",position="append",status="old") 
    if (ios == 0) then
       write(unit=1,iostat=ios,fmt=*) "IN CLEANUP"
       write(unit=1,iostat=ios,fmt=*) "Checkpointfile: ",CHECKPOINTFILE
       write(unit=1,iostat=ios,fmt=*) "Lockfile: ",LOCKFILE
       close(unit=1,iostat=ios)
    end if
    inquire(file=CHECKPOINTFILE,exist=file_found)
    if (file_found) then
       open(unit=1,file=CHECKPOINTFILE,iostat=ios,form="unformatted",action="write",position="append",status="old")
       close(unit=1,iostat=ios,status="delete")
    end if
    inquire(file=LOCKFILE,exist=file_found)
    if (file_found) then
       open(unit=1,file=LOCKFILE,iostat=ios,form="unformatted",action="write",position="append",status="old")
       close(unit=1,iostat=ios,status="delete")
    end if
    call date_and_time(time=endtime)  
    open(unit=1,file=smglobal%logfile,iostat=ios,action="write",position="append",status="old") 
    if (ios == 0) then
       write (unit=1,iostat=ios,fmt=*) "Time started:   ",STARTTIME(1:2),":",STARTTIME(3:4),":",STARTTIME(5:10)
       write (unit=1,iostat=ios,fmt=*) "Time completed: ",endtime(1:2),":",endtime(3:4),":",endtime(5:10)
       close (unit=1,iostat=ios)
    end if
  end subroutine cleanup


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! LOAD_DATA subroutine
!
! Format: load_data
!
! Loads data and user options from the appropriate files, and 
! stores it in various global variables.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine load_data()
    character(len=12) :: datafile,ufile,parmfile="parm.dat    ",resultfile,logfile, &
         covtypelong="opg         ",eqtypelong="low         ",searchmethlong="dfp         ", &
         rhotypelong=".false.     "
    character(len=80) :: toss
    character(len=1) :: covmat_type,equilibrium_type,search_method,rho_type
    logical :: load_u,file_found,locked
    integer :: i,j,nobs,nvar,ngroups,maxgroupsize,nsim,restarts,numagg,istart=1, &
         nofriends,ns,nt,bstar_size,ios
    logical :: fix_gamma
    real(kind=DP) :: fixed_rho,fixed_gamma,dfpstop=-1.0e50_dp,t,rt,eps,shift
    real(kind=DP), dimension(:,:), allocatable :: dat
    integer, dimension(:), allocatable :: tmpsize
    logical, dimension(:), allocatable :: msk
!---------------------------------------------------------------
! First check to see if a previous run of the program was 
! interrupted.  If so, restart that run and skip everything
! else.
!---------------------------------------------------------------
    inquire(file=checkpointfile,exist=resume_from_checkpoint)
    if (resume_from_checkpoint) then
       inquire(file=lockfile,exist=locked)
       if (locked) stop "ERROR: A copy of SMLE is already running in this directory"
       call load_checkpoint() 
       inquire(file=SMGLOBAL%LOGFILE,exist=file_found)
       if (file_found) then
          open(unit=1,file=smglobal%logfile,iostat=ios,action="write",position="append",status="old") 
       else
          open(unit=1,file=smglobal%logfile,iostat=ios,action="write",status="new") 
       end if
       if (ios == 0) then
          write (unit=1,iostat=ios,fmt=*) "Time restarted:   ",starttime(1:2),":",starttime(3:4),":",starttime(5:10)
          close (unit=1,iostat=ios)
       end if
    else
!---------------------------------------------------------------
! Read in parameter file (always PARM.DAT).
! The parameter file has a specific line-by-line format; 
! for example, the program always looks in line #5 for 
! NVAR.
!---------------------------------------------------------------
       inquire (file=parmfile,exist=file_found)
       if (.not.file_found) then 
          stop "Error: parm.dat file not found"
       end if
       open(unit=1,file=parmfile,iostat=ios,action="read",position="rewind",status="old")
       if (ios /= 0) then
          stop "Error: parm.dat file exists but cannot be opened"
       else
          read (unit=1,fmt=*,iostat=ios) toss ! TOSS is used for comment lines; the program doesn't use it
          read (unit=1,fmt=*,iostat=ios) toss
          read (unit=1,fmt=*,iostat=ios) datafile ! DATAFILE = name of file where data is stored
          read (unit=1,fmt=*,iostat=ios) toss
          read (unit=1,fmt=*,iostat=ios) resultfile ! RESULT file = name of file where program should append results
          read (unit=1,fmt=*,iostat=ios) toss
          read (unit=1,fmt=*,iostat=ios) logfile ! LOGFILE = name of file where program should write log information
          read (unit=1,fmt=*,iostat=ios) toss
          read (unit=1,fmt=*,iostat=ios) nobs ! NOBS = number of observations
          read (unit=1,fmt=*,iostat=ios) toss
          read (unit=1,fmt=*,iostat=ios) nvar ! NVAR = number of explanatory variables
          read (unit=1,fmt=*,iostat=ios) toss
          read (unit=1,fmt=*,iostat=ios) numagg ! NUMAGG = number of explanatory variables that are aggregates
          read (unit=1,fmt=*,iostat=ios) toss
          read (unit=1,fmt=*,iostat=ios) nsim ! NSIM = number of simulations to use for GHK
          read (unit=1,fmt=*,iostat=ios) toss
          read (unit=1,fmt=*,iostat=ios) searchmethlong ! "DFP", or "SA"
          read (unit=1,fmt=*,iostat=ios) toss
          read (unit=1,fmt=*,iostat=ios) restarts ! RESTARTS = number of times to restart search with a random vector
          read (unit=1,fmt=*,iostat=ios) toss
          read (unit=1,fmt=*,iostat=ios) eqtypelong ! EQUILIBRIUM_TYPE: high,low,random
          read (unit=1,fmt=*,iostat=ios) toss
          read (unit=1,fmt=*,iostat=ios) covtypelong ! "OPG", "HESSIAN", or "NONE" 
          read (unit=1,fmt=*,iostat=ios) toss
          read (unit=1,fmt=*,iostat=ios) rhotypelong ! RHO_TYPE = .true. if you want to set correlation in unobservables
          read (unit=1,fmt=*,iostat=ios) toss
          read (unit=1,fmt=*,iostat=ios) fixed_rho ! FIXED_RHO = correlation in unobservables (ignored if RHO_TYPE!="FIXED")
          read (unit=1,fmt=*,iostat=ios) toss
          read (unit=1,fmt=*,iostat=ios) fix_gamma ! FIX_GAMMA = .true. if you want to set peer effect
          read (unit=1,fmt=*,iostat=ios) toss
          read (unit=1,fmt=*,iostat=ios) fixed_gamma ! FIXED_GAMMA = peer effect (ignored if FIX_GAMMA=.false.)
          read (unit=1,fmt=*,iostat=ios) toss
          read (unit=1,fmt=*,iostat=ios) load_u ! LOAD_U = .true. if you want to use random numbers from a file
          read (unit=1,fmt=*,iostat=ios) toss
          read (unit=1,fmt=*,iostat=ios) ufile ! UFILE = if (LOAD_U=.true.), the name of the file where the random
                                               ! numbers can be found.  If (LOAD_U=.false.), the name of the file
                                               ! where program should write the random numbers it generates.
          close(unit=1,iostat=ios) ! We're done with PARM.DAT  
!---------------------------------------------------------------
! Perform basic validation and write out information to log file.
!---------------------------------------------------------------
       end if
       open(unit=1,file=logfile,iostat=ios,action="write",status="replace") ! opens the log file, writes start time
       if (ios == 0) then
          write (unit=1,iostat=ios,fmt=*) "SMLE2 Version ",version_number
          write (unit=1,iostat=ios,fmt=*) "Author: Brian Krauth, Simon Fraser University"
          write (unit=1,iostat=ios,fmt=*) "Time started:   ",starttime(1:2),":",starttime(3:4),":",starttime(5:10)
          write (unit=1,iostat=ios,fmt=*) "Data:"
          write (unit=1,iostat=ios,fmt=*) "  Data in file ",datafile
          write (unit=1,iostat=ios,fmt=*) "  There are",nobs," observations"
          write (unit=1,iostat=ios,fmt=*) "  There are",nvar," explanatory variables"
          if (numagg > 0) then
             write (unit=1,iostat=ios,fmt=*) "  Of them,",numagg," are aggregate variables"
          end if
          if (load_u) then
             write (unit=1,iostat=ios,fmt=*) "  Random numbers will be read from file ",ufile
          end if
          write (unit=1,iostat=ios,fmt=*) "Output:"
          write (unit=1,iostat=ios,fmt=*) "  Results will be written to file ",resultfile
          write (unit=1,iostat=ios,fmt=*) "  Log file is ",logfile
          if (.not.load_u) then
             write (unit=1,iostat=ios,fmt=*) "  Random numbers will be written to file ",ufile
          end if
          write (unit=1,iostat=ios,fmt=*) "Estimation parameters:"
          write (unit=1,iostat=ios,fmt=*) "  Likelihood will be estimated using  ",nsim," simulations"
          write (unit=1,iostat=ios,fmt=*) "  Search algorithm will run",restarts," times"
          close (unit=1,iostat=ios)
       end if
       if (numagg < 0) then
          open(unit=1,file=logfile,iostat=ios,action="write",position="append",status="old") 
          if (ios == 0) then
             write (unit=1,iostat=ios,fmt=*) "Fatal error in parm.dat: numagg < 0 not allowed"
             close (unit=1,iostat=ios)
          end if
          stop "Fatal error in parm.dat: numagg < 0 not allowed"
       end if
       if (nvar < (numagg+1)) then
          open(unit=1,file=logfile,iostat=ios,action="write",position="append",status="old") 
          if (ios == 0) then
             write (unit=1,iostat=ios,fmt=*) "Fatal error in parm.dat: numagg > nvar-1 not allowed"
             close (unit=1,iostat=ios)
          end if
          stop "Fatal error in parm.dat: numagg > nvar-1 not allowed"
       end if
       if (nobs < nvar) then
          open(unit=1,file=logfile,iostat=ios,action="write",position="append",status="old") 
          if (ios == 0) then
             write (unit=1,iostat=ios,fmt=*) "Fatal error in parm.dat: nobs < nvar not allowed"
             close (unit=1,iostat=ios)
          end if
          stop "Fatal error in parm.dat: nobs < nvar not allowed"
       end if
       if (nsim < 1) then 
          open(unit=1,file=logfile,iostat=ios,action="write",position="append",status="old") 
          if (ios == 0) then
             write (unit=1,iostat=ios,fmt=*) "Fatal error in parm.dat: nsim < 1 not allowed"
             close (unit=1,iostat=ios)
          end if
          stop "Fatal error in parm.dat: nsim < 1 not allowed"
       end if
       if (restarts < 1) then 
          open(unit=1,file=logfile,iostat=ios,action="write",position="append",status="old") 
          if (ios == 0) then
             write (unit=1,iostat=ios,fmt=*) "Fatal error in parm.dat: restarts < 1 not allowed"
             close (unit=1,iostat=ios)
          end if
          stop "Fatal error in parm.dat: restarts < 1 not allowed"
       end if
       covmat_type=uppercase(covtypelong(1:1))
       if (covmat_type == "N") then
          open(unit=1,file=logfile,iostat=ios,action="write",position="append",status="old") 
          if (ios == 0) then
             write (unit=1,iostat=ios,fmt=*) "  No covariance matrix to be estimated"
             close (unit=1,iostat=ios)
          end if
       elseif (covmat_type == "H") then
          open(unit=1,file=logfile,iostat=ios,action="write",position="append",status="old") 
          if (ios == 0) then
             write (unit=1,iostat=ios,fmt=*) "  Covariance matrix to be estimated by Hessian method"
             close (unit=1,iostat=ios)
          end if
       elseif (covmat_type == "O") then
          open(unit=1,file=logfile,iostat=ios,action="write",position="append",status="old") 
          if (ios == 0) then
             write (unit=1,iostat=ios,fmt=*) "  Covariance matrix to be estimated by OPG/BHHH method"
             close (unit=1,iostat=ios)
          end if
       else
          open(unit=1,file=logfile,iostat=ios,action="write",position="append",status="old") 
          if (ios == 0) then
             write (unit=1,iostat=ios,fmt=*) "Fatal error in parm.dat: covmat_type=",covtypelong," not allowed"
             close (unit=1,iostat=ios)
          end if
          stop "Fatal error in parm.dat: chosen covmat_type not allowed"
       end if
       rho_type=uppercase(rhotypelong(1:1))
       if (rho_type == ".") then                                 ! Added for backward compatibility
          if (uppercase(rhotypelong(2:2)) == "T") then
             rho_type="F"
          elseif (uppercase(rhotypelong(2:2)) == "F") then
             rho_type="X"
          end if
       end if
       if (rho_type == "X") then ! X
          open(unit=1,file=logfile,iostat=ios,action="write",position="append",status="old") 
          if (ios == 0) then
             write (unit=1,iostat=ios,fmt=*) "  Correlation in unobservables set to equal correlation in observables (AET)"
             close (unit=1,iostat=ios)
          end if
       elseif (rho_type == "F") then ! FIXED
          if (fixed_rho > MAXRHO) then
             fixed_rho = MAXRHO
          end if
          open(unit=1,file=logfile,iostat=ios,action="write",position="append",status="old") 
          if (ios == 0) then
             write (unit=1,iostat=ios,fmt=*) "  Correlation in unobservables set to:", fixed_rho
             close (unit=1,iostat=ios)
          end if
       elseif (rho_type == "E") then ! ESTIMATED
          open(unit=1,file=logfile,iostat=ios,action="write",position="append",status="old") 
          if (ios == 0) then
             write (unit=1,iostat=ios,fmt=*) "  Correlation in unobservables to be estimated (warning: weakly identified)"
             close (unit=1,iostat=ios)
          end if
       elseif (rho_type == "I") then ! INTERVAL
          open(unit=1,file=logfile,iostat=ios,action="write",position="append",status="old") 
          if (ios == 0) then
             write (unit=1,iostat=ios,fmt=*) "  Correlation in unobservables to be varied for interval estimation"
             close (unit=1,iostat=ios)
          end if
          if (fix_gamma) then
             open(unit=1,file=logfile,iostat=ios,action="write",position="append",status="old") 
             if (ios == 0) then
                write (unit=1,iostat=ios,fmt=*) "Warning: rho_type=I conflicts with fix_gamma=.true."
                write (unit=1,iostat=ios,fmt=*) "fix_gamma set to .false."
                close (unit=1,iostat=ios)
             end if
             fix_gamma=.false.
          end if
       else  ! SOMETHING ELSE (NOT ALLOWED)
          open(unit=1,file=logfile,iostat=ios,action="write",position="append",status="old") 
          if (ios == 0) then
             write (unit=1,iostat=ios,fmt=*) "Fatal error in parm.dat: RHO_TYPE=",rhotypelong," not allowed"
             close (unit=1,iostat=ios)
          end if
          stop "Fatal error in parm.dat: value of RHO_TYPE is not allowed"
       end if
       if (fix_gamma) then
          if (fixed_gamma < mingam) then
             fixed_gamma=mingam
          end if
          open(unit=1,file=logfile,iostat=ios,action="write",position="append",status="old") 
          if (ios == 0) then
             write (unit=1,iostat=ios,fmt=*) "  Peer effect set to:", fixed_gamma
             close (unit=1,iostat=ios)
          end if
       else
          open(unit=1,file=logfile,iostat=ios,action="write",position="append",status="old") 
          if (ios == 0) then
             write (unit=1,iostat=ios,fmt=*) "  Peer effect to be estimated"
             if (rho_type=="E") then
                write (unit=1,iostat=ios,fmt=*) "WARNING: parameters are only weakly identified!"
             end if
             close (unit=1,iostat=ios)
          end if
       end if
       equilibrium_type=uppercase(eqtypelong(1:1))
       if (equilibrium_type == "L") then
          open(unit=1,file=logfile,iostat=ios,action="write",position="append",status="old") 
          if (ios == 0) then
             write (unit=1,iostat=ios,fmt=*) "  Equilibrium selection rule: always play low-activity equilibrium"
             close (unit=1,iostat=ios)
          end if
       elseif (equilibrium_type == "H") then
          open(unit=1,file=logfile,iostat=ios,action="write",position="append",status="old") 
          if (ios == 0) then
             write (unit=1,iostat=ios,fmt=*) "  Equilibrium selection rule: always play high-activity equilibrium"
             close (unit=1,iostat=ios)
          end if
       elseif (equilibrium_type == "R") then
          open(unit=1,file=logfile,iostat=ios,action="write",position="append",status="old") 
          if (ios == 0) then
             write (unit=1,iostat=ios,fmt=*) "  Equilibrium selection rule: play randomly selected equilibrium"
             close (unit=1,iostat=ios)
          end if
       elseif ((equilibrium_type == "B").or.(equilibrium_type == "M")) then
          open(unit=1,file=logfile,iostat=ios,action="write",position="append",status="old") 
          if (ios == 0) then
             write (unit=1,iostat=ios,fmt=*) "  No equilibrium selection rule: estimate likelihood bounds"
             close (unit=1,iostat=ios)
          end if
          if (rho_type == "I") then 
             open(unit=1,file=logfile,iostat=ios,action="write",position="append",status="old") 
             if (ios == 0) then
                write (unit=1,iostat=ios,fmt=*) "Warning: conflict between equilibrium_type=", &
                     equilibrium_type," and rho_type=I"
                write (unit=1,iostat=ios,fmt=*) "rho_type set to X"             
                close (unit=1,iostat=ios)
             end if
             rho_type = "X"
          end if
          fix_gamma=.false.
          covmat_type="N"
       elseif (equilibrium_type == "P") then
          open(unit=1,file=logfile,iostat=ios,action="write",position="append",status="old") 
          if (ios == 0) then
             write (unit=1,iostat=ios,fmt=*) "  No equilibrium selection rule: plot likelihood bounds"
             close (unit=1,iostat=ios) 
          end if
          if (rho_type == "I") then 
             open(unit=1,file=logfile,iostat=ios,action="write",position="append",status="old") 
             if (ios == 0) then
                write (unit=1,iostat=ios,fmt=*) "Warning: conflict between equilibrium_type=P and rho_type=I"
                write (unit=1,iostat=ios,fmt=*) "rho_type set to X"             
                close (unit=1,iostat=ios)
             end if
             rho_type = "X"
          end if
          fix_gamma=.true.
          covmat_type="N"
       else
          open(unit=1,file=logfile,iostat=ios,action="write",position="append",status="old") 
          if (ios == 0) then
             write (unit=1,iostat=ios,fmt=*) "Fatal error in parm.dat: equilibrium_type=",eqtypelong," not allowed"
             close (unit=1,iostat=ios)
          end if
          stop "Fatal error in parm.dat: chosen equilibrium_type not allowed"
       end if       
       search_method=uppercase(searchmethlong(1:1))
       if (search_method == "D") then
          open(unit=1,file=logfile,iostat=ios,action="write",position="append",status="old") 
          if (ios == 0) then
             write (unit=1,iostat=ios,fmt=*) "  Optimization by Davidson-Fletcher-Powell (DFP) algorithm"
             close (unit=1,iostat=ios)
          end if
       elseif (search_method == "S") then
          open(unit=1,file=logfile,iostat=ios,action="write",position="append",status="old") 
          if (ios == 0) then
             write (unit=1,iostat=ios,fmt=*) "  Optimization by simulated annealing (SA) algorithm"
             close (unit=1,iostat=ios)
          end if
       else
          open(unit=1,file=logfile,iostat=ios,action="write",position="append",status="old") 
          if (ios == 0) then
             write (unit=1,iostat=ios,fmt=*) "Fatal error in parm.dat: search_method=",searchmethlong," not allowed"
             close (unit=1,iostat=ios)
          end if
          stop "Fatal error in parm.dat: chosen search_method not allowed"
       end if       
!---------------------------------------------------------------
! Read in data file, and store the results in various global
! variables
!---------------------------------------------------------------
       inquire (file=datafile,exist=file_found)
       if (file_found) then 
          open (unit=1,file=datafile,iostat=ios,form="formatted",action="read",position="rewind",status="old")
       else
          open(unit=1,file=logfile,iostat=ios,action="write",position="append",status="old") 
          if (ios == 0) then
             write (unit=1,iostat=ios,fmt=*) "Fatal error: datafile ",datafile," not found"
             close (unit=1,iostat=ios)
          end if
          stop "Error: datafile not found"
       end if
!---------------------------------------------------------------
! We didn't know how big our data was until reading NOBS
! and NVAR from the PARM.DAT file.  So now we need to 
! allocate these variables.
!---------------------------------------------------------------
       bstar_size=nvar+3
       if (fix_gamma) then
          bstar_size=bstar_size-1
       end if
       if (rho_type=="E") then
          bstar_size=bstar_size+1
       end if
       if (.not.allocated(ygroup)) then
          allocate(ygroup(nobs),y(nobs),x(nobs,nvar),original_groupid(nobs),groupid(nobs), &
               bx(nobs),bstar(bstar_size),covmat(bstar_size,bstar_size))
       end if
       allocate(tmpsize(nobs),msk(nobs),dat(nobs,nvar+2))
       bx=0.0_dp
       bstar=0.0_dp
       covmat=0.0_dp
!---------------------------------------------------------------
! Read the data file into the array DAT.
!---------------------------------------------------------------
       do i=1,nobs
          read (unit=1,fmt=*,iostat=ios) dat(i,:)
          if (ios /= 0) then 
             close(unit=1,iostat=ios)
             open(unit=1,file=logfile,iostat=ios,action="write",position="append",status="old") 
             write (unit=1,iostat=ios,fmt=*) "Fatal error: datafile has too few records"
             close (unit=1,iostat=ios)
             stop "Fatal error: datafile has too few records"
          end if
       end do
       close(unit=1)
!---------------------------------------------------------------
! Create the main data matrices (global variables): 
!   X (explanatory individual-level variables)
!   Y (endogenous binary choice)
!---------------------------------------------------------------
       y = dat(:,2)
       x=dat(:,3:size(dat,2))
!---------------------------------------------------------------
! Next, we need to create the list of peer groups.  In the 
! input data, groups are identified with some integer ID number.
! For convenience we will redefine the ID number on a scale
! of 1 to NGROUPS (where NGROUPS is the number of distinct 
! groups).  We will create the global variables:
! 
!  ORIGINAL_GROUPID    For each observation, original ID number 
!                      of peer group from data
!  GROUPID             For each observation, internally used ID
!                      number of peer group
!  GROUPSIZE           For each multiple-member 
!                      peer group, number of members
!  NGROUPS             Number of distinct multiple-member
!                      peer groups
!  MAXGROUPSIZE        Size of largest peer group
!
! One complication is that the data may have individuals who
! are the only members of their peer group.  I code these
! respondents as being members of "group zero".  Group zero
! is not included in any group-level calculation, and is thus
! basically dropped from the analysis.  I recommend not
! including them in the input data at all, and the program
! issues a warning whenever there are such respondents.
! The number of respondents who are in one-member groups
! is stored in the local variable NOFRIENDS.
!---------------------------------------------------------------
       original_groupid = nint(dat(:,1))
       groupid=0
       ygroup=0.0_dp
       j=0
       nofriends=0
       do i=minval(original_groupid),maxval(original_groupid)
          msk = (original_groupid == i)
          if (any(msk)) then
             if (count(msk) > 1) then 
                j=j+1
                tmpsize(j) = count(msk)
                where (msk) 
                   groupid=j
                   ygroup=(sum(y,mask=msk)-y)/real(count(msk)-1,kind=DP)
                end where
             else
                nofriends=nofriends+1
             end if
          end if
       end do
       if (nofriends > 0) then
          open(unit=1,file=logfile,iostat=ios,action="write",position="append",status="old") 
          if (ios == 0) then
             write (unit=1,iostat=ios,fmt=*) "Warning: ",nofriends, "observations are the only members of their "
             write (unit=1,iostat=ios,fmt=*) "peer group; not included in estimation. "
             close (unit=1,iostat=ios)
          end if
       end if
       ngroups=j
       if (.not.allocated(groupsize)) then
          allocate(groupsize(ngroups))
       end if
       groupsize=tmpsize(1:ngroups)
       maxgroupsize=maxval(groupsize)
!---------------------------------------------------------------
! Next we create (if needed) the matrix of aggregate
! variables Z.  One issue with this is that the aggregate
! variables are reported separately by each group member.
! In other words, they should be identical across members
! of the same group, but currently the program has no 
! mechanism for checking this.
!---------------------------------------------------------------
       if (numagg > 0) then
          if (.not.allocated(z)) then
             allocate(z(ngroups,numagg))
          end if
          do i=1,ngroups
             do j=1,numagg
                z(i,j) = sum(x(:,j),mask=(groupid==i))/real(groupsize(i),kind=DP)
             end do
          end do
       else
          if (.not.allocated(z)) then
             allocate(z(1,1))
          end if
          z=0.0_dp
       end if
!---------------------------------------------------------------
! Finally we either generate (if LOAD_U=.false.) or read in
! (if LOAD_U=.true.) the random numbers to drive the simulator.
!---------------------------------------------------------------
       if (.not.allocated(u)) then
          allocate(u(maxgroupsize,nsim))
       end if
       if (load_u) then ! load in random numbers from the file named in UFILE
          inquire(file=ufile,exist=file_found)
          if (file_found) then
             open (unit=1,file=ufile,iostat=ios,form="formatted",action="read",position="rewind",status="old")
             if (ios == 0) then
                do i=1,nsim
                   read (unit=1,iostat=ios,fmt=*) u(:,i)
                   if (ios /= 0) then
                      stop "Error: ufile does not have enough records"
                   end if
                end do
                close(unit=1,iostat=ios)
             else
                stop "Error: ufile exists but cannot be opened"
             end if
          else
             stop "Error: ufile not found"
          end if
       else ! generate random numbers with rhalt
          call runif(shift)
          u=transpose(rhalt(nsim,maxgroupsize,shift))
          open(unit=1,file=ufile,iostat=ios,form="formatted",action="write",status="replace")
          if (ios == 0) then
             do i=1,nsim
                write (unit=1,iostat=ios,fmt=*) u(:,i)
             end do
             close(unit=1,iostat=ios)
          end if
       end if
       smglobal=smcontrols(nobs,nvar,ngroups,maxgroupsize,nsim,restarts,numagg,istart,1, &
            fix_gamma,fixed_rho,fixed_gamma,dfpstop,0.0_dp,0.0_dp,0.0_dp,logfile,resultfile, &
            covmat_type,equilibrium_type,search_method,rho_type)
       deallocate(dat,msk,tmpsize)
    end if
    if (SMGLOBAL%SEARCH_METHOD == "S") then
       ns=10
       nt=10
       t=2.0_dp
       rt=0.5_dp
       eps=0.01_dp
       inquire (file="parmsa.dat",exist=file_found)
       if (file_found) then 
          open (unit=1,file="parmsa.dat",iostat=ios,form="formatted",action="read",position="rewind",status="old")
          if (ios == 0) then
             read (unit=1,iostat=ios,fmt=*) toss ! TOSS is used for comment lines; the program doesn't use it
             read (unit=1,iostat=ios,fmt=*) toss
             read (unit=1,iostat=ios,fmt=*) ns
             read (unit=1,iostat=ios,fmt=*) toss
             read (unit=1,iostat=ios,fmt=*) nt
             read (unit=1,iostat=ios,fmt=*) toss
             read (unit=1,iostat=ios,fmt=*) t
             read (unit=1,iostat=ios,fmt=*) toss
             read (unit=1,iostat=ios,fmt=*) rt
             read (unit=1,iostat=ios,fmt=*) toss
             read (unit=1,iostat=ios,fmt=*) eps
             close (unit=1,iostat=ios)
          end if
       else
          continue
       end if
       saglobal=sacontrols(ns,nt,t,rt,eps)
    end if
end subroutine load_data



function uppercase(str) result (strup)
  character(len=1), intent(in) :: str
  character(len=1) :: strup
  strup=str
  if ((ichar(str) >= ichar("a")) .and. (ichar(str) <= ichar("z"))) then 
     strup = char(ichar("A")+ichar(str)-ichar("a"))
  end if
end function uppercase

end module smle2










