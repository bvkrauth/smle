!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! SMUTIL Module
! Author: Brian Krauth, Simon Fraser University
!
! Miscellaneous procedures for SMLE program.  
!
! Public procedures available are:
!
!  	LOAD_DATA		Load user options from PARM.DAT and
!				data from data file.
!
!	RESAMPLE(dat)		Given a data set of size N (dat),
!				sample N times with replacement and
!				return the result. 		
!
!	UNDERREPORTING		Figure out if underreporting correction
! 				is needed and calculate the reporting
! 				rate.
!
!	ESTIMATE	        Estimate the model, returning the results
!				in bstar and fstar.
!
!       GET_STARTB(y,ygroup,x)	Get a good starting value for the
!				parameter vector.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module smutil
  use bklib, only : dp,runif 
  use bkmath, only : pdfn,cdfinvn,ols,randint,rhalt
  use smglob
  use loglik, only : loglikelihood,dloglikelihood
  use dfpmin, only : dfp
  use simann, only : sa
  implicit none
  private
  private :: uppercase
  public :: get_startb,load_data,resample,underreporting,estimate,cleanup 


contains


  subroutine cleanup(starttime)
    logical :: file_found
    character(len=80), intent(in) :: starttime
    character(len=80) :: endtime
    integer :: ios
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
       write (unit=1,iostat=ios,fmt=*) "Time started:   ",starttime(1:2),":",starttime(3:4),":",starttime(5:10)
       write (unit=1,iostat=ios,fmt=*) "Time completed: ",endtime(1:2),":",endtime(3:4),":",endtime(5:10)
       close (unit=1,iostat=ios)
    end if
    deallocate(bstar,u,x,y,ygroup,reporting_rate,nfriends,yint,ygint)
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
    character(len=12) :: datafile,ufile,bootfile,parmfile="parm.dat    ",resultfile,logfile,eqtypelong,rhotypelong,simtypelong
    character(len=80) :: toss
    character(len=1) :: equilibrium_type,rho_type,search_method,simulator_type
    logical :: load_u,file_found,locked
    real(kind=DP), dimension(:,:), allocatable :: dat
    integer :: i,nobs,nvar,maxgroupsize,nsim,restarts,numagg,fixedeffects=0,istart=1,bstarsize,ns,nt,ios
    logical :: underreporting_correction,bootstrap,fix_gamma
    real(kind=DP) :: fixed_rho,fixed_gamma,dfpstop=-1.0e50_dp,t,rt,eps,shift
!---------------------------------------------------------------
! First check to see if a previous run of the program was 
! interrupted.  If so, restart that run.
!---------------------------------------------------------------
    inquire(file=checkpointfile,exist=resume_from_checkpoint)
    if (resume_from_checkpoint) then
       inquire(file=lockfile,exist=locked)
       if (locked) stop "ERROR: A copy of SMLE is already running in this directory"
       call load_checkpoint()
    else
!---------------------------------------------------------------
! Start reading in parameter file (always PARM.DAT).
! The parameter file has a specific line-by-line format; 
! for example, the program always looks in line #5 for 
! NVAR.
!---------------------------------------------------------------
       inquire (file=parmfile,exist=file_found)
       if (.not.file_found) then 
          stop "Error: parm.dat file not found"
       end if
       open(unit=1,file=parmfile,iostat=ios,action="read",position="rewind",status="old")
       if (ios == 0) then
          read (unit=1,fmt=*,iostat=ios) toss ! TOSS is used for comment lines; the program doesn't use it
          read (unit=1,fmt=*,iostat=ios) toss
          read (unit=1,fmt=*,iostat=ios) toss
          read (unit=1,fmt=*,iostat=ios) toss 
          read (unit=1,fmt=*,iostat=ios) nvar ! NVAR = number of explanatory variables
          read (unit=1,fmt=*,iostat=ios) toss 
          read (unit=1,fmt=*,iostat=ios) nobs ! NOBS = number of observations
          read (unit=1,fmt=*,iostat=ios) toss 
          read (unit=1,fmt=*,iostat=ios ) numagg ! NUMAGG = number of explanatory variables that are aggregates
          read (unit=1,fmt=*,iostat=ios) toss 
          read (unit=1,fmt=*,iostat=ios) nsim ! NSIM = number of simulations to use for GHK
          read (unit=1,fmt=*,iostat=ios) toss 
          read (unit=1,fmt=*,iostat=ios) restarts ! RESTARTS = number of times to restart search with a random vector
          read (unit=1,fmt=*,iostat=ios) toss 
          read (unit=1,fmt=*,iostat=ios) simtypelong ! formerly GTOL!
          read (unit=1,fmt=*,iostat=ios) toss 
          read (unit=1,fmt=*,iostat=ios) eqtypelong ! EQUILIBRIUM_TYPE = not operational (yet!)
          read (unit=1,fmt=*,iostat=ios) toss 
          read (unit=1,fmt=*,iostat=ios) underreporting_correction ! UNDERREPORTING_CORRECTION = .true. if correction wanted
          read (unit=1,fmt=*,iostat=ios) toss 
          read (unit=1,fmt=*,iostat=ios) bootstrap ! BOOTSTRAP = .true. if you want to resample data
          read (unit=1,fmt=*,iostat=ios) toss 
          read (unit=1,fmt=*,iostat=ios) load_u ! LOAD_U = .true. if you want to use random numbers from a file
          read (unit=1,fmt=*,iostat=ios) toss 
          read (unit=1,fmt=*,iostat=ios) rhotypelong ! RHO_TYPE = X, F(ixed), E(stimated), or I(nterval)
          read (unit=1,fmt=*,iostat=ios) toss 
          read (unit=1,fmt=*,iostat=ios) fixed_rho ! FIXED_RHO = correlation in unobservables (ignored if FIX_RHO=.false.)
          read (unit=1,fmt=*,iostat=ios) toss 
          read (unit=1,fmt=*,iostat=ios) fix_gamma ! FIX_GAMMA = .true. if you want to set peer effect
          read (unit=1,fmt=*,iostat=ios) toss 
          read (unit=1,fmt=*,iostat=ios) fixed_gamma ! FIXED_GAMMA = peer effect (ignored if FIX_GAMMA=.false.)
          read (unit=1,fmt=*,iostat=ios) toss 
          read (unit=1,fmt=*,iostat=ios) datafile ! DATAFILE = name of file where data is stored
          read (unit=1,fmt=*,iostat=ios) toss 
          read (unit=1,fmt=*,iostat=ios) logfile ! LOGFILE = name of file where program should write log information
          read (unit=1,fmt=*,iostat=ios) toss 
          read (unit=1,fmt=*,iostat=ios) resultfile ! RESULT file = name of file where program should append results
          read (unit=1,fmt=*,iostat=ios) toss 
          read (unit=1,fmt=*,iostat=ios) ufile ! UFILE = if (LOAD_U=.true.), the name of the file where the random
          read (unit=1,fmt=*,iostat=ios) toss 
          ! numbers can be found.  If (LOAD_U=.false.), the name of the file
          ! where program should write the random numbers it generates.
          read (unit=1,fmt=*,iostat=ios) toss 
          read (unit=1,fmt=*,iostat=ios) bootfile ! BOOTFILE = file for bootstraps
          read (unit=1,fmt=*,iostat=ios) toss 
          read (unit=1,fmt=*,iostat=ios) fixedeffects ! FIXEDEFFECTS = number of aggregate variables used for underreporting
          close(unit=1,iostat=ios) ! We're done with PARM.DAT       
       else
          stop "Error: parm.dat file exists but cannot be opened"
       end if
!---------------------------------------------------------------
! A little bit of validation:
!---------------------------------------------------------------
       open(unit=1,file=logfile,iostat=ios,action="write",status="replace") ! opens the log file, writes start time
       if (ios /= 0) then
          stop "Error: Cannot open log file"
       else
          write (unit=1,iostat=ios,fmt=*) "SMLE ",version_number
          write (unit=1,iostat=ios,fmt=*) "Author: Brian Krauth, Simon Fraser University"
          write (unit=1,iostat=ios,fmt=*) "Time started:   ",starttime(1:2),":",starttime(3:4),":",starttime(5:10)
          if (nvar < 1) then
             write (unit=1,iostat=ios,fmt=*) "Error in parm.dat: nvar < 1 not allowed"
             stop "Error in parm.dat: nvar < 1 not allowed"
          end if
          if (nobs < 1) then
             write (unit=1,iostat=ios,fmt=*) "Error in parm.dat: nobs < 1 not allowed"
             stop "Error in parm.dat: nobs < 1 not allowed"
          end if
          if (numagg < 0) then
             write (unit=1,iostat=ios,fmt=*) "Error in parm.dat: numagg < 0 not allowed"
             stop "Error in parm.dat: numagg < 0 not allowed"
          elseif (numagg >= nvar) then
             write (unit=1,iostat=ios,fmt=*) "Error in parm.dat: numagg >= nvar not allowed"
             stop "Error in parm.dat: numagg >= nvar not allowed"
          end if
          if (nsim < 1) then
             write (unit=1,iostat=ios,fmt=*) "Error in parm.dat: nsim < 1 not allowed"
             stop "Error in parm.dat: nsim < 1 not allowed"
          end if
          simulator_type = uppercase(simtypelong(1:1))
          if (restarts < 1) then
             search_method="S"
          else 
             search_method="D"
          end if
          if (fixedeffects > numagg) then 
             write (unit=1,iostat=ios,fmt=*) "Error in parm.dat: fixedeffects>numagg not allowed"
             stop "Error in parm.dat: fixedeffects>numagg not allowed"
          end if
          write (unit=1,iostat=ios,fmt=*) "Data: "
          write (unit=1,iostat=ios,fmt=*) "     Data is in file ",datafile
          write (unit=1,iostat=ios,fmt=*) "     There are ",nobs, " observations"
          write (unit=1,iostat=ios,fmt=*) "     There are ",nvar," explanatory variables "
          if (numagg > 0) then
             write (unit=1,iostat=ios,fmt=*) "     Of these,",numagg," are aggregate variables"
          end if
          if (load_u) then
             write (unit=1,iostat=ios,fmt=*) "    Random numbers will be loaded from file ",ufile
             if (bootstrap) then
                write (unit=1,iostat=ios,fmt=*) "    The bootstrap sample will be loaded from file ",bootfile
             end if
          end if
          write (unit=1,iostat=ios,fmt=*) "Output: "
          write (unit=1,iostat=ios,fmt=*) "    Results will be reported in file ",resultfile
          write (unit=1,iostat=ios,fmt=*) "    The log file is ",logfile
          if (.not.load_u) then
             write (unit=1,iostat=ios,fmt=*) "    Random numbers will be saved in file ",ufile," for later replication"
             if (bootstrap) then
                write (unit=1,iostat=ios,fmt=*) "    The bootstrap sample will be stored in file ",bootfile
             end if
          end if
          write (unit=1,iostat=ios,fmt=*) "Calculation/simulation parameters:"
          if (simulator_type == "G") then
             write (unit=1,iostat=ios,fmt=*) "    The GHK method will be used for estimating selection probabilities"
          elseif (simulator_type == "H") then
             write (unit=1,iostat=ios,fmt=*) "    The GHK/CFS hybrid method will be used for estimating selection probabilities"
          else
             write (unit=1,iostat=ios,fmt=*) "WARNING: Program does not recognize option chosen for simulator_type:"
             write (unit=1,iostat=ios,fmt=*) "         ",simtypelong
             write (unit=1,iostat=ios,fmt=*) "         This is probably because of a PARM.DAT file that doesn't match"
             write (unit=1,iostat=ios,fmt=*) "         the current version of the program.  Rather than stopping here,"
             write (unit=1,iostat=ios,fmt=*) "         the program will choose the default simulator."
             simulator_type = "G"
          end if
          write (unit=1,iostat=ios,fmt=*) "    The likelihood function will be estimated using ",nsim," simulation runs"
          if (search_method == "D") then
             write (unit=1,iostat=ios,fmt=*) "    Optimization will be done by the Davidson-Fletcher-Powell (DFP) algorithm"
             write (unit=1,iostat=ios,fmt=*) "    with ",restarts," restarts of search for maximum"
          else !  (search_method == "S") 
             write (unit=1,iostat=ios,fmt=*) "    Optimization will be done by the simulated annealing (SA) algorithm"
          end if
          if (bootstrap) then
             write (unit=1,iostat=ios,fmt=*) "    Bootstrap resample of original data taken"
          end if
          write (unit=1,iostat=ios,fmt=*) "Model characteristics:"
          if (underreporting_correction) then
             write (unit=1,iostat=ios,fmt=*) "    Underreporting correction to be made"
          end if
          bstarsize=nvar+3
          rho_type=uppercase(rhotypelong(1:1))
          if (rho_type == ".") then
             if (uppercase(rhotypelong(2:2)) == "T") then
                rho_type="F"
             else
                rho_type="X"
             end if
          end if
          if (rho_type == "X") then
             write (unit=1,iostat=ios,fmt=*) "    rho_e (within-group correlation in unobservables)"
             write (unit=1,iostat=ios,fmt=*) "    set equal to rho_x (correlation in obserables)"
          elseif (rho_type == "E") then
             write (unit=1,iostat=ios,fmt=*) "    rho_e (within-group correlation in unobservables) estimated"
             bstarsize=bstarsize+1
          elseif (rho_type == "F") then
             if (fixed_rho > MAXRHO) then
                write(unit=1,iostat=ios,fmt=*) "WARNING: fixed_rho > MAXRHO, so fixed_rho reset to MAXRHO=",MAXRHO
                fixed_rho = MAXRHO
             end if
             write (unit=1,iostat=ios,fmt=*) "    rho_e (within-group correlation in unobservables) fixed at: ",fixed_rho
          elseif (rho_type == "I") then
             write (unit=1,iostat=ios,fmt=*) "    rho_e (within-group correlation in unobservables) varied for interval estimation"
             if (fix_gamma) then
                write (unit=1,iostat=ios,fmt=*) "WARNING: rho_type=Interval not consistent ", &
                     "with fix_gamma=.true.; fix_gamma reset to .false."
             end if
          else
             write (unit=1,iostat=ios,fmt=*) "    Fatal error in parm.dat: Illegal value for rho_type: ",rhotypelong
             stop "Fatal error in parm.dat: Illegal value for rho_type"
          end if
          if (fix_gamma) then
             bstarsize=bstarsize-1
             write (unit=1,iostat=ios,fmt=*) "    gamma (peer effect) fixed at: ",fixed_gamma
             if (fixed_gamma < MINGAM) then
                fixed_gamma=MINGAM
             end if
          else
             write (unit=1,iostat=ios,fmt=*) "    gamma (peer effect) to be estimated"
          end if
          equilibrium_type=uppercase(eqtypelong(1:1))
          if (equilibrium_type == "L") then
             write (unit=1,iostat=ios,fmt=*) "    Equilibrium selection rule: always play low-activity equilibrium"
          elseif (equilibrium_type == "H") then
             write (unit=1,iostat=ios,fmt=*) "    Equilibrium selection rule: always play high-activity equilibrium"
             if (simulator_type == "G") then
                simulator_type = "H"
                write (unit=1,iostat=ios,fmt=*) "WARNING: Chosen simulator type (GHK) not yet available for this selection rule"
                write (unit=1,iostat=ios,fmt=*) "Simulator type switched to GHK/CFS hybrid simulator"
             end if
          elseif (equilibrium_type == "R") then
             write (unit=1,iostat=ios,fmt=*) "    Equilibrium selection rule: play randomly selected equilibrium"
             if (simulator_type == "G") then
                simulator_type = "H"
                write (unit=1,iostat=ios,fmt=*) "WARNING: Chosen simulator type (GHK) not yet available for this selection rule"
                write (unit=1,iostat=ios,fmt=*) "Simulator type switched to GHK/CFS hybrid simulator"
             end if
          elseif (equilibrium_type == "P") then
             write (unit=1,iostat=ios,fmt=*) "    Equilibrium selection rule: plot likelihood bounds"
             if (simulator_type == "G") then
                simulator_type = "H"
                write (unit=1,iostat=ios,fmt=*) "WARNING: Chosen simulator type (GHK) not yet available for this selection rule"
                write (unit=1,iostat=ios,fmt=*) "Simulator type switched to GHK/CFS hybrid simulator"
             end if
          elseif (equilibrium_type == "M") then
             write (unit=1,iostat=ios,fmt=*) "    Equilibrium selection rule: use likelihood bounds method (lower bound only)"
             if (simulator_type == "G") then
                simulator_type = "H"
                write (unit=1,iostat=ios,fmt=*) "WARNING: Chosen simulator type (GHK) not yet available for this selection rule"
                write (unit=1,iostat=ios,fmt=*) "Simulator type switched to GHK/CFS hybrid simulator"
             end if
          elseif (equilibrium_type == "B") then
             write (unit=1,iostat=ios,fmt=*) "    Equilibrium selection rule: use likelihood bounds method"
             if (simulator_type == "G") then
                simulator_type = "H"
                write (unit=1,iostat=ios,fmt=*) "WARNING: Chosen simulator type (GHK) not yet available for this selection rule"
                write (unit=1,iostat=ios,fmt=*) "Simulator type switched to GHK/CFS hybrid simulator"
             end if
          elseif (equilibrium_type == ".") then
             write (unit=1,iostat=ios,fmt=*) "    Equilibrium selection rule: always play low-activity equilibrium"
             equilibrium_type = "L"
          else
             write (unit=1,iostat=ios,fmt=*) "Fatal error in parm.dat: equilibrium_type=",eqtypelong," not allowed"
             stop "Fatal error in parm.dat: chosen equilibrium_type not allowed"
          end if
          close (unit=1,iostat=ios)
       end if
!---------------------------------------------------------------
! Read in data file, and store the results in various global
! variables
!---------------------------------------------------------------
       inquire (file=datafile,exist=file_found)
       if (.not.file_found) then 
          stop "datafile not found"
       end if
       open (unit=1,file=datafile,iostat=ios,form="formatted",action="read",position="rewind",status="old")
!---------------------------------------------------------------
! We didn't know how big our data was until reading NOBS
! and NVAR from the PARM.DAT file.  So now we need to 
! allocate these variables.
!---------------------------------------------------------------
       if (ios /= 0) then
          stop "Error: data file found, but cannot be opened"
       else
          allocate(ygroup(nobs),y(nobs),x(nobs,nvar),nfriends(nobs), &
               yint(nobs),ygint(nobs),reporting_rate(nobs),bstar(bstarsize))
!          if (maxgroupsize > 0) then
!             allocate(dat(nobs,nvar+2))
!          else
          allocate(dat(nobs,nvar+3))
!          end if
!---------------------------------------------------------------
! Read in the data file, into the array DAT.
!---------------------------------------------------------------
          do i=1,nobs
             read (unit=1,iostat=ios,fmt=*) dat(i,:)
             if (ios /= 0) then
                stop "Error: Data file has too few records"
             end if
          end do
          close(unit=1,iostat=ios)
       end if
!---------------------------------------------------------------
! Resample DAT if we're bootstrapping.
!---------------------------------------------------------------
       if (bootstrap) then
!          dat = resample(dat,load_u,bootfile)
          call resample(dat,load_u,bootfile)
       end if
!---------------------------------------------------------------
! Y is the binary choice of the individual, YGROUP is the 
! average choice of the group.  Both are real numbers.
!---------------------------------------------------------------
       y = dat(:,1)
       ygroup = dat(:,2)
!---------------------------------------------------------------
! This is code to deal with the case that MAXGROUPSIZE is set in
! the PARM.DAT file
!---------------------------------------------------------------
!       if (maxgroupsize > 0) then
!          nfriends=maxgroupsize-1
!          x=dat(:,3:(nvar+2))
!       else
          nfriends=nint(dat(:,3))
          x=dat(:,4:(nvar+3))
          maxgroupsize=maxval(nfriends)+1
!       end if
       if (maxgroupsize > 11) then
          stop "maxgroupsize > 11, not allowed"
       end if
!---------------------------------------------------------------
! YINT is the integer version of Y, and YGINT is an integer equal
! to the number (not fraction) of friends who choose Y=1.
! It is handy to have integers for calculations.
!---------------------------------------------------------------
       yint = nint(y)
       ygint = nint(ygroup*real(nfriends,kind(ygroup)))
!---------------------------------------------------------------
! Finally we either generate (if LOAD_U=.false.) or read in
! (if LOAD_U=.true.) the random numbers to drive the simulator.
!---------------------------------------------------------------
       allocate(u(maxgroupsize,nsim))
       if (load_u) then ! load in random numbers from the file named in UFILE
          inquire(file=ufile,exist=file_found)
          if (file_found) then
             open (unit=1,file=ufile,iostat=ios,form="formatted",action="read",position="rewind",status="old")
             if (ios == 0) then
                do i=1,nsim
                   read (unit=1,iostat=ios,fmt=*) u(:,i)
                   if (ios /= 0) then
                      stop "Error: ufile has too few records"
                   end if
                end do
                close(unit=1,iostat=ios)
             else
                stop "Error: ufile found but could not be opened"
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
          else
             continue
          end if
       end if
       smglobal=smcontrols(nobs,nvar,maxgroupsize,nsim,restarts,numagg,fixedeffects,istart,1, &
            underreporting_correction,bootstrap,fix_gamma, &
            fixed_rho,fixed_gamma,dfpstop,0.0_dp,0.0_dp,0.0_dp, &
            logfile,resultfile,equilibrium_type,search_method,rho_type,simulator_type)
       deallocate(dat)
    end if
    if (smglobal%search_method == "S") then
       ns = 10
       nt = 10
       t = 2.0_dp
       rt = 0.5_dp
       eps = 0.01_dp
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
          else
             continue
          end if
       else
          continue
       end if
       saglobal=sacontrols(ns,nt,t,rt,eps)
    end if
  end subroutine load_data
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! RESAMPLE subroutine (formerly function)
!
! Format: resample(dat)
!
! Given the N-by-K array DAT, create a new N-by-K array by 
! randomly sampling (with replacement) rows of DAT.
! Used to do bootstrapping.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine resample(dat,loadit,bootfile) 
    real(kind=DP), dimension(:,:), intent(inout) :: dat
    logical, intent(in) :: loadit
    character(len=12), intent(in) :: bootfile
    real(kind=DP), dimension(size(dat,1),size(dat,2)) :: datout
    integer :: i,n,ios
    integer, dimension(size(dat,1)) :: swr
    logical :: file_found
    n=size(swr)
    if (loadit) then
       inquire(file=bootfile,exist=file_found)
       if (file_found) then
          open (unit=1,file=bootfile,iostat=ios,form="formatted",action="read",position="rewind",status="old")
          if (ios == 0) then
             do i=1,n
                read (unit=1,iostat=ios,fmt=*) swr(i)
                if (ios /=0) then
                   stop "Error: bootfile does not contain enough records"
                end if
             end do
             close(unit=1,iostat=ios)
          else
             stop "Error: bootfile found but could not be opened"
          end if
       else
          stop "Error: bootfile not found"
       end if
    else
!       swr=randint(n)
       call randint(swr)
       open(unit=1,file=bootfile,iostat=ios,form="formatted",action="write",status="replace")
       if (ios == 0) then
          do i=1,n
             write (unit=1,iostat=ios,fmt=*) swr(i)
          end do
          close(unit=1,iostat=ios)
       else
          continue
       end if
    end if
    do i=1,n
       datout(i,:)=dat(swr(i),:)
    end do
    dat = datout
  end subroutine resample

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! UNDERREPORTING subroutine
!
! Format: underreporting
!
! Sets up the underreporting correction.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine underreporting()
    real(kind=DP) :: reprate
    real(kind=DP), dimension(smglobal%nobs,smglobal%fixedeffects+1) :: z
    real(kind=DP), dimension(smglobal%fixedeffects+1) :: b,bgroup
    integer :: ios
    if (smglobal%underreporting_correction) then
       reprate = sum(y)/sum(ygroup)
       open(unit=1,file=smglobal%logfile,iostat=ios,action="write",position="append",status="old")
       if (ios == 0) then
          write (unit=1,iostat=ios,fmt=*) "Underreporting correction:"
          write (unit=1,iostat=ios,fmt=*) "    Overall reporting rate: ",reprate
!---------------------------------------------------------------
! If REPORTING_RATE > 1, then there is no underreporting to
! correct, so shut the whole thing off.  Because the underreporting
! correction takes up a lot of CPU time, this is much better
! than just setting REPORTING_RATE to 1.
!---------------------------------------------------------------
          if (reprate > 1.0_dp) then
             smglobal%underreporting_correction = .false.
             write (unit=1,fmt=*) "underreporting correction turned off"
          else 
             z(:,1)=1.0_dp
             if (smglobal%fixedeffects > 0) then
                z(:,2:smglobal%fixedeffects+1)=x(:,1:smglobal%fixedeffects)
             end if
             b=ols(z,y)
             write (unit=1,iostat=ios,fmt=*) "    Self-reported avg. choice",b
             bgroup=ols(z,ygroup)
             write (unit=1,iostat=ios,fmt=*) "    Peer-reported avg. choice",bgroup
             !          reporting_rate = reprate
             reporting_rate=matmul(z,b)/matmul(z,bgroup)
             if (smglobal%fixedeffects > 0) then 
                write (unit=1,iostat=ios,fmt=*) "    Reporting rate by group:", reporting_rate
             end if
          end if
          close (unit=1,iostat=ios)
       else
          stop "Error: Unable to write to log file"
       end if
    end if
  end subroutine underreporting


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
    real(kind=DP), dimension(size(bstar)) :: startb,b,c,lb,ub
    real(kind=DP), dimension(size(x,2)) :: mu,sigmasq
    real(kind=DP) :: fret,n
    integer :: istart,i,iter,nacc,nfcnev,nobds,ier,intercept,ios
    call underreporting()
!---------------------------------------------------------------
! Standardize x for estimation (prevents numerical problems?) added experimentally 1/23
!---------------------------------------------------------------
    n=real(size(x,1),kind=DP)
    do i=1,size(mu)
       mu(i)=sum(x(:,i))/n
       sigmasq(i)=sum((x(:,i)-mu(i))**2)/n
       x(:,i)=(x(:,i)-mu(i))/sqrt(sigmasq(i))
    end do
!---------------------------------------------------------------
! Get a good initial guess on the parameter vector.
!---------------------------------------------------------------
    startb = get_startb(y,ygroup,x) 
    b = startb
    if (.not.resume_from_checkpoint) then
       bstar = startb
       fstar = 1.0e50_dp
    end if
    istart=smglobal%istart
!---------------------------------------------------------------
! Start looking
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
       do i = istart,smglobal%restarts ! Start from a few different points to avoid finding only local optimum
          smglobal%istart=i
          call dfp(b,0.00001_dp,iter,fret,loglikelihood,dloglikelihood)
          if (fret < fstar) then ! If this iteration of DFP has improved on our highest likelihood so far...
             fstar = fret ! ...then save that result
             bstar = b
          end if
          call runif(c)
          b = 2.0_dp*c*b ! Generate a random (but not crazy) starting point for restart.
       end do
    else
       open(unit=1,file=SMGLOBAL%LOGFILE,iostat=ios,action="write",position="append",status="old")
       if (ios == 0) then
          write (unit=1,iostat=ios,fmt=*) "Error in ESTIMATE: search_method not defined"
          close(unit=1,iostat=ios)
       else
          continue
       end if
       stop "Error in ESTIMATE: search_method not defined"
    end if
    intercept=size(BSTAR)-size(mu)
    do i=1,size(mu)
       bstar(i+intercept)=bstar(i+intercept)/sqrt(sigmasq(i))
       bstar(intercept)=bstar(intercept)-bstar(i+intercept)*mu(i)
       x(:,i)=(sqrt(sigmasq(i))*x(:,i))+mu(i)
    end do
    smglobal%istart=1
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
    real(kind=DP), dimension(Size(BSTAR)) :: startb
    real(kind=DP), dimension(size(y)) :: ytmp
    real(kind=DP), dimension(Size(x,1),Size(x,2)+2) :: z
    integer :: i
    z(:,1) = ygroup
    z(:,2) = 1.0_dp
    z(:,3:Size(z,2)) = x
    ytmp=y-0.5_dp
    i=1+size(startb)-size(z,2)
    startb(i:size(startb)) = ols(z,ytmp)/pdfn(cdfinvn(sum(y)/real(size(y),DP)))
    startb(1) =  0.2_dp ! rho_x
    i=2
    if ((SMGLOBAL%RHO_TYPE == "E").or.((SMGLOBAL%RHO_TYPE=="I").and.(SMGLOBAL%RUN_NUMBER == 2))) then
       startb(i) = 0.2_dp
       i=i+1
    end if   
    if ((SMGLOBAL%FIX_GAMMA).or.((SMGLOBAL%RHO_TYPE=="I").and.(SMGLOBAL%RUN_NUMBER == 2))) then
       continue
    else
       startb(i) = startb(i)/2.0_dp 
    end if
  end function get_startb



function uppercase(str) result (strup)
  character(len=1), intent(in) :: str
  character(len=1) :: strup
  strup=str
  if ((ichar(str) >= ichar("a")) .and. (ichar(str) <= ichar("z"))) then 
     strup = char(ichar("A")+ichar(str)-ichar("a"))
  end if
end function uppercase




end module smutil

