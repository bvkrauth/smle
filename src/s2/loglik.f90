!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! LOGLIK Module
! Author: Brian Krauth, Simon Fraser University
!
! Calculates log-likelihood function for SMLE program.  
!
! Public procedures available are:
! 
!   	LOGLIKELIHOOD(p):	Log-likelihod of the data associated with the 
!				parameter vector p.
!	DLOGLIKELIHOOD(p):	Numerical approximation of the first
!				derivative of the log-likelihood function.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module loglik
use bklib, only : DP
use bkmath, only : cdfn,cdfinvn,chol
use smglob, only : SMGLOBAL, Y, Z, X, BX, U, GROUPSIZE, GROUPID, MAXRHO, MINGAM
implicit none
private
private :: probi,check_constraints
public :: loglikelihood,dloglikelihood,loglikevec


contains
  

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! LOGLIKELIHOOD function
!
! Format: loglikelihood(p)
!
! Gives loglikelihod of data given parameter vector p.  P must
! be double-precision.  The calcuation itself uses a number
! of additional variables which are defined in the module
! SMGLOB.
!
! Most of the work is done in the LOGLIKEVEC function
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!function loglikelihood(p) result (like)
subroutine loglikelihood(p,like) ! new
  real(kind=DP), dimension(:), intent(in) :: p
!  real(kind=DP) :: like
  real(kind=DP), intent(inout) :: like ! new
  real(kind=DP), dimension(SMGLOBAL%NGROUPS) :: likevec ! new
  integer :: ios
!  like=sum(loglikevec(p))
  call loglikevec(p,likevec) ! new
  like=sum(likevec) ! new
!---------------------------------------------------------------
! Then we write the result to the log file
!---------------------------------------------------------------
  open (unit=1,file=SMGLOBAL%LOGFILE,iostat=ios,action="write",position="append",status="old")
  if (ios == 0) then
     write (unit=1,iostat=ios,fmt="(7f12.3)") like,p 
     close (unit=1,iostat=ios)
  end if
!---------------------------------------------------------------
! Finally, the function actually returns the negative of the
! loglikelihood.  This is because our optimizer DFPMIN
! minimizes the function is given, and we want to maximize 
! the likelihood.
!---------------------------------------------------------------
  like=-like ! since dfpmin minimizes
!end function loglikelihood
end subroutine loglikelihood

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! LOGLIKEVEC function
!
! Format: loglikevec(p)
!
! Gives loglikelihod of data by observation given parameter vector p.  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!function loglikevec(p) result(likevec)
subroutine loglikevec(p,likevec) ! new
  real(kind=DP), dimension(:), intent(in) :: p
  real(kind=DP) :: like
  real(kind=DP), dimension(SMGLOBAL%NOBS) :: bxsq
!  real(kind=DP), dimension(SMGLOBAL%NGROUPS) :: bz,likevec,npeer
  real(kind=DP), dimension(SMGLOBAL%NGROUPS) :: bz,npeer ! new
  real(kind=DP), dimension(SMGLOBAL%NGROUPS), intent(inout) :: likevec ! new
  real(kind=DP), dimension(SMGLOBAL%NGROUPS,2) :: bxtmp  
  real(kind=DP), dimension(SMGLOBAL%MAXGROUPSIZE,SMGLOBAL%MAXGROUPSIZE) :: c
  real(kind=DP) :: rhox,rhoe,gam
  real(kind=DP), parameter :: minlike=1.0e-300_dp
  integer :: i,intercept
!---------------------------------------------------------------
! There are various constraints on parameter values; for example,
! gamma cannot be negative.  We implement this by starting with 
! a loglikelihood of zero, then applying penalties for any violations
! of constraints.   This allows us to use a standard optimization 
! method.
!---------------------------------------------------------------
  like=0.0_dp
!---------------------------------------------------------------
! RHOX is the correlation in BX among peers, which is a parameter
! to be estimated.  RHOE is the correlation in unobservables
! among peers.  GAM is the peer effect parameter.  The program
! allows the user to set either RHOE or GAM (not both) 
! directly, or they can be estimated under the RHOX=RHOE 
! identifying assumption.
!---------------------------------------------------------------
  rhox = p(1) 
  intercept=2
  if (SMGLOBAL%RHO_TYPE=="X") then
     rhoe=rhox
  elseif ((SMGLOBAL%RHO_TYPE=="E").or.((SMGLOBAL%RHO_TYPE=="I").and.(SMGLOBAL%RUN_NUMBER==2))) then
     rhoe=p(intercept)
     intercept=intercept+1
  else ! either RHO_TYPE=F(ixed) or RHO_TYPE=I(nterval) and RUN_NUMBER /= 2
     rhoe=SMGLOBAL%FIXED_RHO
  end if
  if ((SMGLOBAL%FIX_GAMMA).or.((SMGLOBAL%RHO_TYPE=="I").and.(SMGLOBAL%RUN_NUMBER==2))) then
     gam=SMGLOBAL%FIXED_GAMMA
  else
     gam=p(intercept)
     intercept=intercept+1
  end if
!---------------------------------------------------------------
! BX is a vector of length NOBS (# of observations);
! BZ (if it exists) is a vector of length NGROUPS (# of groups)
! The exogenous
! explanatory variables (in X) are classified as individual-level 
! or aggregate.  The first NUMAGG columns of X are assumed to
! be aggregate, and the rest are assumed to be individual level.
! Agent i's private utility from choosing one is BZ[GROUPID(i)]+BX[i],
! where BZ[GROUPID(i)] is the product of i's aggregate variables and their
! coefficients, and BX[i] is the product of i's individual
! variables and their coefficients.
! 
! The distinction between the two is that it is assumed that 
! agent i's friends have the same value of BZ, while
! his/her friends' values of BX are random, with correlation 
! RHO_X.
!---------------------------------------------------------------
  if (SMGLOBAL%NUMAGG > 0) then
     bz = matmul(Z,p(intercept+1:(intercept+SMGLOBAL%NUMAGG))) 
  else
     bz=0.0_dp
  end if
  BX=p(intercept) + matmul(X(:,(1+SMGLOBAL%NUMAGG):size(X,2)),p((intercept+1+SMGLOBAL%NUMAGG):size(p)))
!---------------------------------------------------------------
! The procedure CHECK_CONSTRAINTS checks to see if the parameter
! values satisfy the various constraints, and subtracts a penalty
! from LOGLIKELIHOOD if they are not.
!---------------------------------------------------------------
  call check_constraints(like,rhox,rhoe,gam,BX,bz)
!---------------------------------------------------------------
! Now, given parameter values the distribution of private 
! preferences is multivariate normal with mean 
!     BZ[GROUPID(i)] + BX[i]
! and a covariance matrix which has ones on the diagonal
! and RHOE off the diagonal.
! Because the covariance  matrix does not vary across observations, 
! we calculate it here directly, and store its Cholesky 
! decomposition in the matrix C.
!---------------------------------------------------------------
  c=rhoe
  forall (i=1:size(c,1)) c(i,i)=1.0_dp
  c=transpose(chol(c))
! !HPF$ INDEPENDENT
  do i=1,size(likevec)
     likevec(i) = probi(i,GROUPSIZE(i),bz(i),gam,c)
  end do
!---------------------------------------------------------------
! Having calculated the likelihood of each observation, we 
! take the natural log and add across observations to get
! the log-likelihood.
!---------------------------------------------------------------
  where (likevec < minlike) likevec=minlike ! avoid floating point error
  likevec=log(likevec)
  bxsq= BX-(sum(BX)/real(size(BX),kind=DP))
  bxsq = bxsq/sqrt(sum(bxsq**2)/real(size(bxsq),kind=DP))
  forall (i=1:SMGLOBAL%NGROUPS)
     bxtmp(i,1) = sum(bxsq,mask=(GROUPID==i))**2
  end forall
  bxsq = bxsq**2
  forall (i=1:SMGLOBAL%NGROUPS) 
     bxtmp(i,2) = sum(bxsq,mask=(GROUPID==i))
  end forall
  npeer=real(GROUPSIZE-1,kind=DP)
  likevec=likevec-0.5_dp*(npeer*log(1.0_dp-rhox) + log(1.0_dp+npeer*rhox) + &
       bxtmp(:,2)/(1.0_dp-rhox)  - (rhox/((1.0_dp-rhox)*(1.0_dp+npeer*rhox)))*bxtmp(:,1))
  likevec(1)=like+likevec(1)
!end function loglikevec
end subroutine loglikevec

function probi(gidloc,gsloc,bz,gam,c) result (probiout)
  integer, intent(in) :: gidloc,gsloc
  real(kind=DP), intent(in) :: bz,gam
  real(kind=DP), dimension(:,:), intent(in) :: c
  real(kind=DP), dimension(gsloc) :: yloc,bxloc
  real(kind=DP), dimension(gsloc) :: mine,maxe,ypeer,cj
  real(kind=DP), dimension(SMGLOBAL%NSIM) :: ta,tb,wgt,selected
  real(kind=DP), dimension(gsloc,SMGLOBAL%NSIM)  :: tt
  real(kind=DP), parameter :: bignum=1.0e3_dp
  real(kind=DP) :: probiout
  integer :: i,j,ytot
  integer, dimension(gsloc,SMGLOBAL%NSIM) :: cutoff
  logical, dimension(gsloc+1,SMGLOBAL%NSIM) :: equilibrium
!
! First we pull data from the global variables
! YLOC is the vector of choices by members of group GIDLOC
! BXLOC is the vector of their characteristics
! YPEER is the vector of leave-yourself-out averages for behavior
!
  yloc=pack(Y,mask=(GROUPID==gidloc))
  bxloc=pack(BX,mask=(GROUPID==gidloc))
  ypeer=(sum(yloc)-yloc)/real(gsloc-1,kind=DP)
!
! Now, in order for YLOC to be an equilibrium, the vector of preferences
! must lie between MINE and MAXE as calculated below.
!
  where (yloc < 0.1_dp) 
     mine = -bignum
     maxe = -(bz+bxloc+gam*ypeer)
  elsewhere 
     mine = -(bz+bxloc+gam*ypeer)
     maxe = bignum
  end where
!
! This is the start of the GHK simulation
! The purpose is to calculate the probability that the vector
! of preferences is such that YLOC is *an* equlibrium, i.e.,
! that it lies between MINE and MAXE.  In addition, a side
! effect of the GHK simulator is that it returns a randomly
! drawn vector from the interval [MINE,MAXE].  This will 
! be used shortly to calculate the probability that YLOC
! will be the *selected* equilibrium, conditional on being
! *an* equilibrium.
!
  ta = cdfn(mine(1)) ! formerly cdfn(((mine(1))/c(1,1)))
  tb = cdfn(maxe(1)) ! formerly cdfn(((maxe(1))/c(1,1)))
  wgt = tb-ta
  do j=2,gsloc
     cj=c(j,1:gsloc) 
     tt((j-1),:) = cdfinvn(tb+U((j-1),:)*(ta-tb))
     ta = cdfn(((mine(j))-matmul(cj(1:(j-1)),tt(1:(j-1),:)))/c(j,j))
     tb = cdfn(((maxe(j))-matmul(cj(1:(j-1)),tt(1:(j-1),:)))/c(j,j))
     wgt = wgt*(tb-ta)
  end do
  probiout = sum(wgt)/real(size(wgt),DP)
  tt((gsloc),:) = cdfinvn(tb+U((gsloc),:)*(ta-tb))
!
! Now we figure out if our equilibrium was selected
! Step 1 is to calculate, for each simulation run,
! a row of cutoffs for each individual.  For example
! CUTOFF[i,j] is the highest number of friends choosing
! y=1 such that individual j in simulation i will
! prefer to choose y=0.
!
  cutoff=floor((real(1-gsloc,kind=DP)/gam)*(matmul(c(1:gsloc,1:gsloc),tt)+spread(bxloc+bz,dim=2,ncopies=size(tt,2))))
  selected=1.0_dp
!
! The next step is to calculate the set of equilibrium,
! then apply the selection rule to determine the probability
! that YLOC is the selected equilibrium.
!
  if (SMGLOBAL%EQUILIBRIUM_TYPE== "L") then
     ytot=sum(nint(yloc))
     if (ytot > 1) then
        do j=1,(ytot-1)
           where ((count(cutoff < j-2,dim=1)==(j-1)).and.(count(cutoff == (j-2),dim=1)==0)) selected=0.0_dp
        end do
     end if
  elseif (SMGLOBAL%EQUILIBRIUM_TYPE=="H") then 
     ytot=sum(nint(yloc))
     if (ytot < gsloc) then
        do j=ytot+2,(gsloc+1)
           where ((count(cutoff < j-2,dim=1)==(j-1)).and.(count(cutoff == (j-2),dim=1)==0)) selected=0.0_dp
        end do
     end if
  elseif (SMGLOBAL%EQUILIBRIUM_TYPE=="R") then
     do i=1,size(equilibrium,1) 
        do j=1,size(equilibrium,2)
           equilibrium(i,j) = (count(cutoff(:,j) < i-2)==(i-1)).and.(count(cutoff(:,j)==(i-2))==0)
        end do
     end do
     do i=1,size(selected)
        selected(i)=1.0_dp/real(count(equilibrium(:,i)),kind=DP)
     end do
  else ! SMGLOBAL%EQUILILIBRIUM_TYPE == "M", "B", or "P"
     if ((SMGLOBAL%RUN_NUMBER==1).or.((SMGLOBAL%EQUILIBRIUM_TYPE == "P").and.(modulo(SMGLOBAL%RUN_NUMBER,2)>0))) then
        do i=1,size(equilibrium,1) 
           do j=1,size(equilibrium,2)
              equilibrium(i,j) = (count(cutoff(:,j) < i-2)==(i-1)).and.(count(cutoff(:,j)==(i-2))==0)
           end do
        end do
        do i=1,size(selected)
           if (count(equilibrium(:,i)) > 1) then ! there is a faster way of doing this with "dim"
              selected(i)=0.0_dp
           end if
        end do
     end if
  end if
  probiout=probiout*(sum(selected)/real(size(selected),kind=DP)) 
end function probi








!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! CHECK_CONSTRAINTS subroutine (PRIVATE, used by LOGLIKELIHOOD)
!
! Format: check_constraints(loglikelihood,rho,rho2,gam,bx,bz)
!
! Checks to see if the parameters rho,rho2,gam, and the 
! vectors bx and bz satisfy various constraints.  If they do not,
! we change them so they satisfy the constraints exactly, then
! add a penalty to LOGLIKELIHOOD.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine check_constraints(loglikelihood,rho,rho2,gam,bx,bz)
  real(kind=DP), intent(inout) :: loglikelihood 
  real(kind=DP), intent(inout) :: rho,rho2,gam
  real(kind=DP), dimension(:), intent(inout) :: bx,bz
  real(kind=DP), parameter :: penalty=50.0_dp,maxbx=5.0_dp,penalty2=50.0_dp
  real(kind=DP) :: minrho
!---------------------------------------------------------------
! Constraint #1: RHO must be greater than -1/MAXGROUPSIZE, or else
! the covariance matrix fails to be positive definite, and
! thus a valid covariance matrix.
!---------------------------------------------------------------
  minrho = -MAXRHO/real(SMGLOBAL%MAXGROUPSIZE,kind=DP)
  if (rho < minrho) then
     loglikelihood= loglikelihood-penalty*(rho-minrho)**2
     rho = minrho
!---------------------------------------------------------------
! Constraint #2: RHO must be less than (or equal to) one, since 
! it is a correlation. 
!---------------------------------------------------------------
  else if (rho > MAXRHO) then
     loglikelihood = loglikelihood-penalty*(rho-MAXRHO)**2
     rho = MAXRHO
  end if
!---------------------------------------------------------------
! Constraints #3 and #4: RHO2 faces the same constraints as RHO.
!---------------------------------------------------------------
  if (rho2 < minrho) then
     loglikelihood= loglikelihood-penalty*(rho2-minrho)**2
     rho2 = minrho
  else if (rho2 > MAXRHO) then
     loglikelihood = loglikelihood-penalty*(rho2-MAXRHO)**2
     rho2 = MAXRHO
  end if
!---------------------------------------------------------------
! Constraint #5: GAM must be nonnegative. For computational
! purposes, it is easiest if GAM is always greater than or equal to
! some very timy positive number.
!---------------------------------------------------------------
  if (gam < MINGAM) then
     loglikelihood=loglikelihood-penalty*(gam-MINGAM)**2
     gam=MINGAM
  end if
!---------------------------------------------------------------
! Constraint #6: to avoid wasting time with obviously crazy 
! parameter values, we also put a maximum on  |BX| and |BZ|. 
!---------------------------------------------------------------
  if ((maxval(abs(bx)) > maxbx) .or. (maxval(abs(bz)) > maxbx)) then
     loglikelihood = loglikelihood - penalty2*sum((bx-maxbx)**2,(abs(bx) > maxbx))
     where (bx > maxbx) bx = maxbx
     where (bx < -maxbx) bx= -maxbx
     loglikelihood = loglikelihood - penalty2*sum((bz-maxbx)**2,(abs(bz) >maxbx))
     where (bz > maxbx) bz = maxbx
     where (bz < -maxbx) bz=-maxbx
  end if
end subroutine check_constraints





!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! DLOGLIKELIHOOD function
!
! Format: dloglikelihood(p)
!
! Gives approximate first derivative (at p) of 
! LOGLIKELIHOOD(p).
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!function dloglikelihood(p,fp) result(fout)
subroutine dloglikelihood(p,fp,fout) ! new
  real(kind=DP), dimension(:), intent(in) :: p
  real(kind=DP), intent(in) :: fp
!  real(kind=DP), dimension(size(p)) :: fout,ptmp,gradp
  real(kind=DP), dimension(size(p)), intent(inout) :: fout ! new
  real(kind=DP), dimension(size(p)) :: ptmp,gradp ! new
  real(kind=DP) :: ftmp ! new
  integer :: i
  ptmp=abs(p)
  where (ptmp < 1.0e-2_dp) ptmp = 1.0e-2_dp
!
! Changed by BK 5/10/2004 - seems like a big change.  It *looks* like I had the 
! sign function backwards, but things mostly seemed to work.  If things start
! to act up it should probably be changed back
!  gradp = 1.0e-8_dp*ptmp*sign(p,1.0_dp)
  gradp = 1.0e-8_dp*ptmp*sign(1.0_dp,p) ! new 
  ptmp = p+gradp
  gradp = ptmp-p
  do i=1,Size(p)
     ptmp=p
     ptmp(i) = ptmp(i) + gradp(i)
!     fout(i) = (loglikelihood(ptmp)-fp)/gradp(i)
     call loglikelihood(ptmp,ftmp) ! new
     fout(i) = (ftmp-fp)/gradp(i) ! new
  end do
!end function dloglikelihood
end subroutine dloglikelihood



end module loglik







