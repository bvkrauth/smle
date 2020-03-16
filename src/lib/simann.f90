module simann
  use bklib, only : dp,setseed,runif
  use smglob, only : smglobal, save_sa_checkpoint
  implicit none
  private 
  private :: exprep,prt1,prt2,prt3,prt4,prt5,prt6,prt7,prt8,prt9,prt10,prtvec
  public :: sa
   
contains

!  Version: 3.2
!  Date: 1/22/94.
!  Differences compared to Version 2.0:
!     1. If a trial is out of bounds, a point is randomly selected
!        from LB(i) to UB(i). Unlike in version 2.0, this trial is
!        evaluated and is counted in acceptances and rejections.
!        All corresponding documentation was changed as well.
!  Differences compared to Version 3.0:
!     1. If VM(i) > (UB(i) - LB(i)), VM is set to UB(i) - LB(i).
!        The idea is that if T is high relative to LB & UB, most
!        points will be accepted, causing VM to rise. But, in this
!        situation, VM has little meaning; particularly if VM is
!        larger than the acceptable region. Setting VM to this size
!        still allows all parts of the allowable region to be selected.
!  Differences compared to Version 3.1:
!     1. Test made to see if the initial temperature is positive.
!     2. WRITE statements prettied up.
!     3. References to paper updated.
!
!  Synopsis:
!  This routine implements the continuous simulated annealing global
!  optimization algorithm described in Corana et al.'s article
!  "Minimizing Multimodal Functions of Continuous Variables with the
!  "Simulated Annealing" Algorithm" in the September 1987 (vol. 13,
!  no. 3, pp. 262-280) issue of the ACM Transactions on Mathematical
!  Software.
!
!  A very quick (perhaps too quick) overview of SA:
!     SA tries to find the global optimum of an N dimensional function.
!  It moves both up and downhill and as the optimization process
!  proceeds, it focuses on the most promising area.
!     To start, it randomly chooses a trial point within the step length
!  VM (a vector of length N) of the user selected starting point. The
!  function is evaluated at this trial point and its value is compared
!  to its value at the initial point.
!     In a maximization problem, all uphill moves are accepted and the
!  algorithm continues from that trial point. Downhill moves may be
!  accepted; the decision is made by the Metropolis criteria. It uses T
!  (temperature) and the size of the downhill move in a probabilistic
!  manner. The smaller T and the size of the downhill move are, the more
!  likely that move will be accepted. If the trial is accepted, the
!  algorithm moves on from that point. If it is rejected, another point
!  is chosen instead for a trial evaluation.
!     Each element of VM periodically adjusted so that half of all
!  function evaluations in that direction are accepted.
!     A fall in T is imposed upon the system with the RT variable by
!  T(i+1) = RT*T(i) where i is the ith iteration. Thus, as T declines,
!  downhill moves are less likely to be accepted and the percentage of
!  rejections rise. Given the scheme for the selection for VM, VM falls.
!  Thus, as T declines, VM falls and SA focuses upon the most promising
!  area for optimization.
!
!  The importance of the parameter T:
!     The parameter T is crucial in using SA successfully. It influences
!  VM, the step length over which the algorithm searches for optima. For
!  a small intial T, the step length may be too small; thus not enough
!  of the function might be evaluated to find the global optima. The user
!  should carefully examine VM in the intermediate output (set IPRINT =
!  1) to make sure that VM is appropriate. The relationship between the
!  initial temperature and the resulting step length is function
!  dependent.
!     To determine the starting temperature that is consistent with
!  optimizing a function, it is worthwhile to run a trial run first. Set
!  RT = 1.5 and T = 1.0. With RT > 1.0, the temperature increases and VM
!  rises as well. Then select the T that produces a large enough VM.
!
!  For modifications to the algorithm and many details on its use,
!  (particularly for econometric applications) see Goffe, Ferrier
!  and Rogers, "Global Optimization of Statistical Functions with
!  Simulated Annealing," Journal of Econometrics, vol. 60, no. 1/2, 
!  Jan./Feb. 1994, pp. 65-100.
!  For more information, contact 
!              Bill Goffe
!              Department of Economics and International Business
!              University of Southern Mississippi 
!              Hattiesburg, MS  39506-5072 
!              (601) 266-4484 (office)
!              (601) 266-4920 (fax)
!              bgoffe@whale.st.usm.edu (Internet)
!
!  As far as possible, the parameters here have the same name as in
!  the description of the algorithm on pp. 266-8 of Corana et al.
!
!  In this description, SP is single precision, DP is double precision,
!  INT is integer, L is logical and (N) denotes an array of length n.
!  Thus, DP(N) denotes a double precision array of length n.
!
!  Input Parameters:
!    Note: The suggested values generally come from Corana et al. To
!          drastically reduce runtime, see Goffe et al., pp. 90-1 for
!          suggestions on choosing the appropriate RT and NT.
!    X - The starting values for the variables of the function to be
!        optimized. (DP(N))
!    MAX - Denotes whether the function should be maximized or
!          minimized. A true value denotes maximization while a false
!          value denotes minimization. Intermediate output (see IPRINT)
!          takes this into account. (L)
!    RT - The temperature reduction factor. The value suggested by
!         Corana et al. is .85. See Goffe et al. for more advice. (DP)
!    EPS - Error tolerance for termination. If the final function
!          values from the last neps temperatures differ from the
!          corresponding value at the current temperature by less than
!          EPS and the final function value at the current temperature
!          differs from the current optimal function value by less than
!          EPS, execution terminates and IER = 0 is returned. (EP)
!    NS - Number of cycles. After NS*N function evaluations, each
!         element of VM is adjusted so that approximately half of
!         all function evaluations are accepted. The suggested value
!         is 20. (INT)
!    NT - Number of iterations before temperature reduction. After
!         NT*NS*N function evaluations, temperature (T) is changed
!         by the factor RT. Value suggested by Corana et al. is
!         MAX(100, 5*N). See Goffe et al. for further advice. (INT)
!    NEPS - Number of final function values used to decide upon termi-
!           nation. See EPS. Suggested value is 4. (INT)
!    MAXEVL - The maximum number of function evaluations. If it is
!             exceeded, IER = 1. (INT)
!    LB - The lower bound for the allowable solution variables. (DP(N))
!    UB - The upper bound for the allowable solution variables. (DP(N))
!         If the algorithm chooses X(I) .LT. LB(I) or X(I) .GT. UB(I),
!         I = 1, N, a point is from inside is randomly selected. This
!         This focuses the algorithm on the region inside UB and LB.
!         Unless the user wishes to concentrate the search to a par-
!         ticular region, UB and LB should be set to very large positive
!         and negative values, respectively. Note that the starting
!         vector X should be inside this region. Also note that LB and
!         UB are fixed in position, while VM is centered on the last
!         accepted trial set of variables that optimizes the function.
!    C - Vector that controls the step length adjustment. The suggested
!        value for all elements is 2.0. (DP(N))
!    IPRINT - controls printing inside SA. (INT)
!             Values: 0 - Nothing printed.
!                     1 - Function value for the starting value and
!                         summary results before each temperature
!                         reduction. This includes the optimal
!                         function value found so far, the total
!                         number of moves (broken up into uphill,
!                         downhill, accepted and rejected), the
!                         number of out of bounds trials, the
!                         number of new optima found at this
!                         temperature, the current optimal X and
!                         the step length VM. Note that there are
!                         N*NS*NT function evalutations before each
!                         temperature reduction. Finally, notice is
!                         is also given upon achieveing the termination
!                         criteria.
!                     2 - Each new step length (VM), the current optimal
!                         X (XOPT) and the current trial X (X). This
!                         gives the user some idea about how far X
!                         strays from XOPT as well as how VM is adapting
!                         to the function.
!                     3 - Each function evaluation, its acceptance or
!                         rejection and new optima. For many problems,
!                         this option will likely require a small tree
!                         if hard copy is used. This option is best
!                         used to learn about the algorithm. A small
!                         value for MAXEVL is thus recommended when
!                         using IPRINT = 3.
!             Suggested value: 1
!             Note: For a given value of IPRINT, the lower valued
!                   options (other than 0) are utilized.
!
!  Input/Output Parameters:
!    T - On input, the initial temperature. See Goffe et al. for advice.
!        On output, the final temperature. (DP)
!
!  Output Parameters:
!    XOPT - The variables that optimize the function. (DP(N))
!    FOPT - The optimal value of the function. (DP)
!    NACC - The number of accepted function evaluations. (INT)
!    NFCNEV - The total number of function evaluations. In a minor
!             point, note that the first evaluation is not used in the
!             core of the algorithm; it simply initializes the
!             algorithm. (INT).
!    NOBDS - The total number of trial function evaluations that
!            would have been out of bounds of LB and UB. Note that
!            a trial point is randomly selected between LB and UB.
!            (INT)
!    IER - The error return number. (INT)
!          Values: 0 - Normal return; termination criteria achieved.
!                  1 - Number of function evaluations (NFCNEV) is
!                      greater than the maximum number (MAXEVL).
!                  2 - The starting value (X) is not inside the
!                      bounds (LB and UB).
!                  3 - The initial temperature is not positive.
!                  99 - Should not be seen; only used internally.
!
!
!  Required Functions (included):
!    EXPREP - Replaces the function EXP to avoid under- and overflows.
!             It may have to be modified for non IBM-type main-
!             frames. (DP)
!
!  Required Subroutines (included):
!    PRTVEC - Prints vectors.
!    PRT1 ... PRT10 - Prints intermediate output.
!    FCN - Function to be optimized. The form is
!            SUBROUTINE FCN(N,X,F)
!            INTEGER N
!            DOUBLE PRECISION  X(N), F
!            ...
!            function code with F = F(X)
!            ...
!            RETURN
!            END
!          Note: This is the same form used in the multivariable
!          minimization algorithms in the IMSL edition 10 library.
!
!  Machine Specific Features:
!    1. EXPREP may have to be modified if used on non-IBM type main-
!       frames. Watch for under- and overflows in EXPREP.
!    2. Some FORMAT statements use G25.18; this may be excessive for
!       some machines.


  subroutine sa(x,max,rt,eps,ns,nt,neps,maxevl,lb,ub,c,iprint, &
       t,xopt,fopt,nacc,nfcnev,nobds,ier,fcn)


!  Type all external variables.
    real(kind=DP), dimension(:), intent(inout) :: x
    real(kind=DP), dimension(size(x)), intent(in) :: lb,ub,c
    real(kind=DP), dimension(size(x)), intent(out) :: xopt
    real(kind=DP), intent(inout) :: t
    real(kind=DP), intent(in) :: eps, rt
    real(kind=DP), intent(out) :: fopt
    integer, intent(in) ::  ns, nt, neps,  maxevl, iprint
    integer, intent(out) :: nacc,nfcnev,nobds,ier
    logical, intent(in) :: max

!  Type all internal variables.
    real(kind=DP), dimension(neps) ::  fstar 
    real(kind=DP), dimension(size(x)) ::  xp,vm 
    integer, dimension(size(x)) :: nacp 
    real(kind=DP) :: f, fp, p, pp, ratio, shift
    integer :: n, nup, ndown, nrej, nnew, lnobds, h, i, j, m, ios
    logical :: quit

!  Type all functions.
!      DOUBLE PRECISION  EXPREP
!      REAL  RANMAR
    interface
!       function fcn(p) result (fout)
       subroutine fcn(p,fout) ! new
         use bklib
         real(kind=dp), dimension(:), intent(in) :: p
!         real(kind=dp) :: fout
         real(kind=dp), intent(inout) :: fout ! new
!       end function fcn
       end subroutine fcn ! new
    end interface

!  Initialize the random number generator RANMAR. (deleted)
!  call rmarin(iseed1,iseed2)
    call setseed()
  
!  Set initial values.
    n = size(x)
    nacc = 0
    nobds = 0
    nfcnev = 0
    ier = 99
    xopt = x
    nacp = 0
    fstar = 1.0e20_dp
    vm = ub-lb

!  If the initial temperature is not positive, notify the user and 
!  return to the calling routine.  
    if (t <= 0.0_dp) then
       open (unit=1,file=SMGLOBAL%LOGFILE,iostat=ios,action="write",position="append",status="old")
       if (ios == 0) then
          write(unit=1,iostat=ios,fmt=*) "  THE INITIAL TEMPERATURE IS NOT POSITIVE. "
          write(unit=1,iostat=ios,fmt=*) "  RESET THE VARIABLE T. "
          close(unit=1,iostat=ios)
       else
          continue
       end if
       ier = 3
       return
    end if

!  If the initial value is out of bounds, notify the user and return
!  to the calling routine.
    if (any(x > ub) .or. any(x < lb)) then
       call prt1()
       ier = 2
       return
    end if

!  Evaluate the function with input X and return value as F.
!    f = fcn(x)
    call fcn(x,f) ! new

!  If the function is to be minimized, switch the sign of the function.
!  Note that all intermediate and final output switches the sign back
!  to eliminate any possible confusion for the user.
    if(.not. max) f = -f
    nfcnev = nfcnev + 1
    fopt = f
    fstar(1) = f
    if(iprint >= 1) call prt2(max,n,x,f)

!  Start the main loop. Note that it terminates if (i) the algorithm
!  succesfully optimizes the function or (ii) there are too many
!  function evaluations (more than MAXEVL).
    do
       nup = 0
       nrej = 0
       nnew = 0
       ndown = 0
       lnobds = 0
       do m = 1, nt
          do j = 1, ns
             do h = 1, n
              
!  Generate XP, the trial value of X. Note use of VM to choose XP.
                do i = 1, n
                   if (i == h) then
                      call runif(shift)
                      xp(i) = x(i) + (shift*2.0_dp- 1.0_dp) * vm(i)
                   else
                      xp(i) = x(i)
                   end if

!  If XP is out of bounds, select a point in bounds for the trial.
                   if((xp(i) < lb(i)) .or. (xp(i) > ub(i))) then
                      call runif(shift)
                      xp(i) = lb(i) + (ub(i) - lb(i))*shift
                      lnobds = lnobds + 1
                      nobds = nobds + 1
                      if(iprint >= 3) then
                         call prt3(max,n,xp,x,fp,f)
                      end if
                   end if
                end do

!  Evaluate the function with the trial point XP and return as FP.
!                fp= fcn(xp)
                call fcn(xp,fp)
                if(.not. max) fp = -fp
                nfcnev = nfcnev + 1
                if(iprint >= 3) call prt4(max,n,xp,x,fp,f)

!  If too many function evaluations occur, terminate the algorithm.
                if(nfcnev >= maxevl) then
                   call prt5()
                   if (.not. max) then
                      fopt = -fopt
                   end if
                   ier = 1
                   return
                end if

!  Accept the new point if the function value increases.
                if(fp >= f) then
                   if(iprint >= 3) then
                      open (unit=1,file=SMGLOBAL%LOGFILE,iostat=ios,action="write",position="append",status="old")
                      if (ios == 0) then
                         write(unit=1,iostat=ios,fmt=*) "  POINT ACCEPTED"
                         close (unit=1,iostat=ios)
                      else
                         continue   
                      end if
                   end if
                   x=xp
                   f = fp
                   nacc = nacc + 1
                   nacp(h) = nacp(h) + 1
                   nup = nup + 1

!  If greater than any other point, record as new optimum.
                   if (fp > fopt) then
                      if(iprint >= 3) then
                         open (unit=1,file=SMGLOBAL%LOGFILE,iostat=ios,action="write",position="append",status="old")
                         if (ios == 0) then
                            write(unit=1,iostat=ios,fmt=*) "  NEW OPTIMUM"
                            close (unit=1,iostat=ios)
                         else
                            continue
                         end if
                      end if
                      xopt=xp
                      fopt = fp
                      nnew = nnew + 1
                   end if

!  If the point is lower, use the Metropolis criteria to decide on
!  acceptance or rejection.
                else
                   p = exprep((fp - f)/t)
                   call runif(pp)
                   if (pp < p) then
                      if(iprint >= 3) call prt6(max)
                      x=xp
                      f = fp
                      nacc = nacc + 1
                      nacp(h) = nacp(h) + 1
                      ndown = ndown + 1
                   else
                      nrej = nrej + 1
                      if(iprint >= 3) then
                         call prt7(max)
                      end if
                   end if
                end if
             end do
          end do


!  Adjust VM so that approximately half of all evaluations are accepted.
          do i = 1, n
             ratio = real(nacp(i),kind=DP) /real(ns,kind=DP)
             if (ratio > 0.6_dp) then
                vm(i) = vm(i)*(1.0_dp + c(i)*(ratio - 0.6_dp)/0.4_dp)
             else if (ratio < 0.4_dp) then
                vm(i) = vm(i)/(1.0_dp + c(i)*((0.4_dp - ratio)/0.4_dp))
             end if
             if (vm(i) > (ub(i)-lb(i))) then
                vm(i) = ub(i) - lb(i)
             end if
          end do
        
          if(iprint >= 2) then
             call prt8(n,vm,xopt,fopt,x)
          end if
        
          nacp=0
        
       end do
     
       if(iprint >= 1) then
          call prt9(max,n,t,xopt,vm,fopt,nup,ndown,nrej,lnobds,nnew)
       end if
     
!  Check termination criteria.
       quit = .false.
       fstar(1) = f
       if ((fopt - fstar(1)) <= eps) quit = .true.
       do i = 1, neps
          if (abs(f - fstar(i)) > eps) quit = .false.
       end do
       if (fopt > (-SMGLOBAL%DFPSTOP)) then 
          open (unit=1,file=SMGLOBAL%LOGFILE,iostat=ios,action="write",position="append",status="old")
          if (ios == 0) then
             write (unit=1,iostat=ios,fmt=*) "fret < dfpstop"
             close (unit=1,iostat=ios)
          else
             continue
          end if
          quit = .true.
       end if

!  Terminate SA if appropriate.
       if (quit) then
          x=xopt
          ier = 0
          if (.not. max) fopt = -fopt
          if(iprint >= 1) call prt10()
          return
       end if

!  If termination criteria is not met, prepare for another loop.
       t = rt*t
       call save_sa_checkpoint(t,xopt)
       do i = neps, 2, -1
          fstar(i) = fstar(i-1)
       end do

       f = fopt
       x=xopt

!  Loop again.
    end do
  
  end subroutine sa

function  exprep(rdum) result (exprepout)
!  This function replaces exp to avoid under- and overflows and is
!  designed for IBM 370 type machines. It may be necessary to modify
!  it for other machines. Note that the maximum and minimum values of
!  EXPREP are such that they has no effect on the algorithm.

  real(kind=DP), intent(in) :: rdum
  real(kind=DP) :: exprepout

  if (rdum > 174.0_dp) then
     exprepout = 3.69e+75_dp
  else if (rdum < -180.0_dp) then
     exprepout = 0.0_dp
  else
     exprepout = exp(rdum)
  end if

end function exprep

!subroutine rmarin(ij,kl)
!  This subroutine and the next function generate random numbers. See
!  the comments for SA for more information. The only changes from the
!  orginal code is that (1) the test to make sure that RMARIN runs first
!  was taken out since SA assures that this is done (this test didn't
!  compile under IBM's VS Fortran) and (2) typing ivec as integer was
!  taken out since ivec isn't used. With these exceptions, all following
!  lines are original.

! This is the initialization routine for the random number generator
!     RANMAR()
! NOTE: The seed variables can have values between:    0 <= IJ <= 31328
!                                                      0 <= KL <= 30081
!  real :: U(97), C, CD, CM
!  integer :: I97, J97
!  common /raset1/ U, C, CD, CM, I97, J97
!  if( IJ < 0  .or.  IJ .gt. 31328  .or. KL .lt. 0  .or.  KL .gt. 30081 ) then
!     print "(A)", " The first random number seed must have a value between 0 and 31328"
!     print "(A)"," The second seed must have a value between 0 and 30081"
!     stop
!  endif
!  i = mod(IJ/177, 177) + 2
!  j = mod(IJ    , 177) + 2
!  k = mod(KL/169, 178) + 1
!  l = mod(KL,     169)
!  do ii = 1, 97
!     s = 0.0
!     t = 0.5
!     do jj = 1, 24
!        m = mod(mod(i*j, 179)*k, 179)
!        i = j
!        j = k
!        k = m
!        l = mod(53*l+1, 169)
!        if (mod(l*m, 64) >= 32) then
!           s = s + t
!        endif
!        t = 0.5 * t
!     end do
!     U(ii) = s
!  end do
!  C = 362436.0 / 16777216.0
!  CD = 7654321.0 / 16777216.0
!  CM = 16777213.0 /16777216.0
!  I97 = 97
!  J97 = 33
!end subroutine rmarin

!function ranmar()
!  real :: U(97), C, CD, CM
!  integer :: I97, J97
!  common /raset1/ U, C, CD, CM, I97, J97
!  uni = U(I97) - U(J97)
!  if( uni < 0.0 ) uni = uni + 1.0
!  U(I97) = uni
!  I97 = I97 - 1
!  if(I97 == 0) I97 = 97
!  J97 = J97 - 1
!  if(J97 == 0) J97 = 97
!  C = C - CD
!  if( C .lt. 0.0 ) C = C + CM
!  uni = uni - C
!  if( uni .lt. 0.0 ) uni = uni + 1.0
!  ranmar = uni
!end function ranmar

subroutine prt1()
!  This subroutine prints intermediate output, as does PRT2 through
!  PRT10. Note that if SA is minimizing the function, the sign of the
!  function value and the directions (up/down) are reversed in all
!  output to correspond with the actual function optimization. This
!  correction is because SA was written to maximize functions and
!  it minimizes by maximizing the negative a function.
  integer :: ios
  open (unit=1,file=SMGLOBAL%LOGFILE,iostat=ios,action="write",position="append",status="old")
  if (ios == 0) then
     write (unit=1,iostat=ios,fmt=*) "  THE STARTING VALUE (X) IS OUTSIDE THE BOUNDS "
     write (unit=1,iostat=ios,fmt=*) "  (LB AND UB). EXECUTION TERMINATED WITHOUT ANY"
     write (unit=1,iostat=ios,fmt=*) "  OPTIMIZATION. RESPECIFY X, UB OR LB SO THAT  "
     write (unit=1,iostat=ios,fmt=*) "  LB(I) < X(I) < UB(I), I = 1, N. "
     close (unit=1,iostat=ios)
  else
     continue
  end if

end subroutine prt1

subroutine prt2(max,n,x,f)
  real(kind=DP), dimension(:), intent(in)  ::  x
  real(kind=DP), intent(in) ::  f
  integer, intent(in)  :: n
  logical, intent(in) ::  max
  integer :: ios

  call prtvec(x,n,"INITIAL X")
  open (unit=1,file=SMGLOBAL%LOGFILE,iostat=ios,action="write",position="append",status="old")
  if (ios == 0) then
     if (max) then
        write(unit=1,iostat=ios,fmt=*) "  INITIAL F: ", f
     else
        write(unit=1,iostat=ios,fmt=*) "  INITIAL F: ", -f
     end if
     close(unit=1,iostat=ios)
  else
     continue
  end if

end subroutine prt2

subroutine prt3(max,n,xp,x,fp,f)

  real(kind=DP), dimension(:), intent(in) ::  xp, x
  real(kind=DP), intent(in) ::  fp, f
  integer, intent(in) ::  n
  logical, intent(in) ::  max
  integer :: ios

  call prtvec(x,n,"CURRENT X")
  open (unit=1,file=SMGLOBAL%LOGFILE,iostat=ios,action="write",position="append",status="old")
  if (ios == 0) then
     if (max) then
        write(unit=1,iostat=ios,fmt=*) "  CURRENT F: ", f
     else
        write(unit=1,iostat=ios,fmt=*) "  CURRENT F: ", -f
     end if
     close(unit=1,iostat=ios)
  else
     continue
  end if
  call prtvec(xp,n,"TRIAL X")
  open (unit=1,file=SMGLOBAL%LOGFILE,iostat=ios,action="write",position="append",status="old")
  if (ios == 0) then
     write(unit=1,iostat=ios,fmt=*) "  POINT REJECTED SINCE OUT OF BOUNDS"
     close(unit=1,iostat=ios)
  else
     continue
  end if

end subroutine prt3

subroutine prt4(max,n,xp,x,fp,f)

  real(kind=DP), dimension(:), intent(in)  ::  xp, x 
  real(kind=DP), intent(in)   ::  fp, f
  integer, intent(in)   :: n
  logical, intent(in)   :: max
  integer :: ios

  call prtvec(x,n,"CURRENT X")
  if (max) then
     open (unit=1,file=SMGLOBAL%LOGFILE,iostat=ios,action="write",position="append",status="old")
     if (ios == 0) then
        write(unit=1,iostat=ios,fmt=*) "  CURRENT F: ",f
        close (unit=1,iostat=ios)
     else
        continue
     end if
     call prtvec(xp,n,"TRIAL X")
     open (unit=1,file=SMGLOBAL%LOGFILE,iostat=ios,action="write",position="append",status="old")
     if (ios == 0) then
        write(unit=1,iostat=ios,fmt=*) "  RESULTING F: ",fp
        close (unit=1,iostat=ios)
     else
        continue
     end if
  else
     open (unit=1,file=SMGLOBAL%LOGFILE,iostat=ios,action="write",position="append",status="old")
     if (ios == 0) then
        write(unit=1,iostat=ios,fmt=*) "  CURRENT F: ",-f
        close (unit=1,iostat=ios)
     else
        continue
     end if
     call prtvec(xp,n,"TRIAL X")
     open (unit=1,file=SMGLOBAL%LOGFILE,iostat=ios,action="write",position="append",status="old")
     if (ios == 0) then
        write(unit=1,iostat=ios,fmt=*) "  RESULTING F: ",-fp
        close (unit=1,iostat=ios)
     else
        continue
     end if
  end if

end subroutine prt4

subroutine prt5()

  integer :: ios

  open (unit=1,file=SMGLOBAL%LOGFILE,iostat=ios,action="write",position="append",status="old")
  if (ios == 0) then
     write(unit=1,iostat=ios,fmt=*) "  TOO MANY FUNCTION EVALUATIONS; CONSIDER "
     write(unit=1,iostat=ios,fmt=*) "  INCREASING MAXEVL OR EPS, OR DECREASING "
     write(unit=1,iostat=ios,fmt=*) "  NT OR RT. THESE RESULTS ARE LIKELY TO BE "
     write(unit=1,iostat=ios,fmt=*) "  POOR."
     close(unit=1,iostat=ios)
  else
     continue
  end if
    
end subroutine prt5


subroutine prt6(max)
  logical, intent(in)   ::  max
  integer :: ios

  open (unit=1,file=SMGLOBAL%LOGFILE,iostat=ios,action="write",position="append",status="old")
  if (ios == 0) then
     if (max) then
        write(unit=1,iostat=ios,fmt=*) "  THOUGH LOWER, POINT ACCEPTED"
     else
        write(unit=1,iostat=ios,fmt=*) "  THOUGH HIGHER, POINT ACCEPTED"
     end if
     close (unit=1,iostat=ios)
  else
     continue
  end if

end subroutine prt6

subroutine prt7(max)

  logical, intent(in)   :: max
  integer :: ios

  open (unit=1,file=SMGLOBAL%LOGFILE,iostat=ios,action="write",position="append",status="old")
  if (ios == 0) then
     if (max) then
        write(unit=1,iostat=ios,fmt=*) "  LOWER POINT REJECTED"
     else
        write(unit=1,iostat=ios,fmt=*) "  HIGHER POINT REJECTED"
     end if
     close(unit=1,iostat=ios)
  else
     continue
  end if

end subroutine prt7

subroutine prt8(n,vm,xopt,fopt,x)

  real(kind=DP), dimension(:), intent(in) :: vm, xopt, x
  real(kind=DP), intent(in) :: fopt
  integer, intent(in) :: n
  integer :: ios

  open (unit=1,file=SMGLOBAL%LOGFILE,iostat=ios,action="write",position="append",status="old")
  if (ios == 0) then
     write(unit=1,iostat=ios,fmt=*) " INTERMEDIATE RESULTS AFTER STEP LENGTH ADJUSTMENT"
     close(unit=1,iostat=ios)
  else
     continue
  end if
  call prtvec(vm,n,"NEW STEP LENGTH (VM)")
  call prtvec(xopt,n,"CURRENT OPTIMAL X")
  open (unit=1,file=SMGLOBAL%LOGFILE,iostat=ios,action="write",position="append",status="old")
  if (ios == 0) then
     write(unit=1,iostat=ios,fmt=*) " CURRENT OPTIMAL F(X): ",fopt
     close(unit=1,iostat=ios)  
  else
     continue
  end if
  call prtvec(x,n,"CURRENT X")

end subroutine prt8

subroutine prt9(max,n,t,xopt,vm,fopt,nup,ndown,nrej,lnobds,nnew)

  real(kind=DP), dimension(:), intent(in) ::  xopt, vm
  real(kind=DP), intent(in) :: t, fopt
  integer, intent(in) :: n, nup, ndown, nrej, lnobds, nnew
  logical, intent(in) :: max
  integer :: totmov,ios

  totmov = nup + ndown + nrej

  open (unit=1,file=SMGLOBAL%LOGFILE,iostat=ios,action="write",position="append",status="old")
  if (ios == 0) then
     write(unit=1,iostat=ios,fmt=*) " INTERMEDIATE RESULTS BEFORE NEXT TEMPERATURE REDUCTION"
     write(unit=1,iostat=ios,fmt=*) "  CURRENT TEMPERATURE:  ",t
     if (max) then
        write(unit=1,iostat=ios,fmt=*) "  MAX FUNCTION VALUE SO FAR: ",fopt
        write(unit=1,iostat=ios,fmt=*) "  TOTAL MOVES:               ", totmov
        write(unit=1,iostat=ios,fmt=*) "    UPHILL:                  ", nup
        write(unit=1,iostat=ios,fmt=*) "    ACCEPTED DOWNHILL:       ", ndown
        write(unit=1,iostat=ios,fmt=*) "    REJECTED DOWNHILL:       ", nrej
        write(unit=1,iostat=ios,fmt=*) "    OUT OF BOUNDS TRIALS:    ", lnobds
        write(unit=1,iostat=ios,fmt=*) "  NEW MAXIMA THIS TEMPERATURE:", nnew
     else
        write(unit=1,iostat=ios,fmt=*) "  MIN FUNCTION VALUE SO FAR: ",-fopt
        write(unit=1,iostat=ios,fmt=*) "  TOTAL MOVES:               ", totmov
        write(unit=1,iostat=ios,fmt=*) "    DOWNHILL:                ", nup
        write(unit=1,iostat=ios,fmt=*) "    ACCEPTED UPHILL:       ", ndown
        write(unit=1,iostat=ios,fmt=*) "    REJECTED UPHILL:       ", nrej
        write(unit=1,iostat=ios,fmt=*) "    OUT OF BOUNDS TRIALS:    ", lnobds
        write(unit=1,iostat=ios,fmt=*) "  NEW MINIMA THIS TEMPERATURE:", nnew
     end if
     close(unit=1,iostat=ios)
  else
     continue
  end if
  call prtvec(xopt,n,"CURRENT OPTIMAL X")
  call prtvec(vm,n,"STEP LENGTH (VM)")
end subroutine prt9

subroutine prt10()

  integer :: ios

  open (unit=1,file=SMGLOBAL%LOGFILE,iostat=ios,action="write",position="append",status="old")
  if (ios == 0) then
     write(unit=1,iostat=ios,fmt=*) "  SA ACHIEVED TERMINATION CRITERIA. IER = 0. "
     close(unit=1,iostat=ios)
  else
     continue
  end if


end subroutine prt10

subroutine prtvec(vector,ncols,name)
!  This subroutine prints the double precision vector named VECTOR.
!  Elements 1 thru NCOLS will be printed. NAME is a character variable
!  that describes VECTOR. Note that if NAME is given in the call to
!  PRTVEC, it must be enclosed in quotes. If there are more than 10
!  elements in VECTOR, 10 elements will be printed on each line.

  integer, intent(in) :: ncols
  real(kind=DP), dimension(ncols), intent(in) :: vector
!  character *(*) name
  character(len=*), intent(in) :: name
  integer :: lines,i,ll,j,ios

  open (unit=1,file=SMGLOBAL%LOGFILE,iostat=ios,action="write",position="append",status="old")
  if (ios == 0) then
     write(unit=1,iostat=ios,fmt=*) name
     if (ncols > 10) then  
        lines = floor(real(ncols)/10.0) !  formerly   lines = int(ncols/10.)
        do  i = 1, lines
           ll = 10*(i - 1)
           write(unit=1,iostat=ios,fmt=*) (VECTOR(J),J = 1+LL, 10+LL)
        end do
        write(unit=1,iostat=ios,fmt=*) (VECTOR(J),J = 11+LL, NCOLS)
     else
        write(unit=1,iostat=ios,fmt=*) (VECTOR(J),J = 1, NCOLS)
     end if
     close(unit=1,iostat=ios)
  else
     continue
  end if


! 1000 format( 10(g12.5,1x))

end subroutine prtvec

end module simann
