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
use bklib, only : sp,dp
use bkmath, only : cdfn,cdfinvn,chol,nchoosek
use smglob, only : smglobal,x,y,yint,ygint,nfriends,reporting_rate,u,bigmatsize,mingam,maxrho
implicit none
private
private :: ghkd,ghkdu,rintersects,rintersectd,rintersecti, &
     nonemptys,nonemptyd,nonemptyi,smprob_hybrid,smprob_ghk, &
     check_constraints,eqregion,eqrange,makew
public :: loglikelihood,dloglikelihood,rintersect,ghk,nonempty
interface ghk
   module procedure ghkd,ghkdu
end interface
interface rintersect
   module procedure rintersects,rintersectd,rintersecti
end interface
interface nonempty
   module procedure nonemptys,nonemptyd,nonemptyi
end interface



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
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!function loglikelihood(p) result(like)
subroutine loglikelihood(p,like) ! new
  real(kind=DP), dimension(:), intent(in) :: p
!  real(kind=DP) :: like
  real(kind=DP), intent(inout) :: like ! new
  integer :: i,intercept,ios
  real(kind=DP), dimension(SMGLOBAL%NOBS) :: bz,bx,likevec,underreporting_rate
  real(kind=DP), dimension(SMGLOBAL%MAXGROUPSIZE,SMGLOBAL%MAXGROUPSIZE) :: c
  real(kind=DP) :: mu,sigmasq,rhox,rhoe,gam
!---------------------------------------------------------------
! First we include some HPF (high-performance fortran) directives
! so that the HPF compiler knows how best to parallelize.
!---------------------------------------------------------------
!HPF$ DISTRIBUTE LIKEVEC(BLOCK)
!HPF$ ALIGN YINT(I) WITH LIKEVEC(I)
!HPF$ ALIGN YGINT(I) WITH LIKEVEC(I)
!HPF$ ALIGN BX(I) WITH LIKEVEC(I)
!HPF$ ALIGN BZ(I) WITH LIKEVEC(I)
!HPF$ ALIGN NFRIENDS(I) WITH LIKEVEC(I)
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
! allows the user to set either RHO2 or GAM (not both) 
! directly, or they can be estimated under the RHOX=RHOE 
! identifying assumption.
!
!---------------------------------------------------------------
  rhox=p(1)
  intercept=2
  if ((SMGLOBAL%RHO_TYPE=="X").or.((SMGLOBAL%RHO_TYPE=="I").and.(SMGLOBAL%RUN_NUMBER==1))) then
     rhoe=rhox
  elseif ((SMGLOBAL%RHO_TYPE=="E").or.((SMGLOBAL%RHO_TYPE=="I").and.(SMGLOBAL%RUN_NUMBER==2))) then
     rhoe=p(intercept)
     intercept=intercept+1
  else ! either RHO_TYPE=F(ixed) or RHO_TYPE=I(nterval) and RUN_NUMBER >  2
     rhoe=SMGLOBAL%FIXED_RHO
  end if
  if ((SMGLOBAL%FIX_GAMMA).or.((SMGLOBAL%RHO_TYPE=="I").and.(SMGLOBAL%RUN_NUMBER==2))) then
     gam=SMGLOBAL%FIXED_GAMMA
  else
     gam=p(intercept)
     intercept=intercept+1
  end if
!---------------------------------------------------------------
! BZ and BX are both vectors of length NOBS.  The exogenous
! explanatory variables (in X) are classified as individual-level 
! or aggregate.  The first NUMAGG columns of X are assumed to
! be aggregate, and the rest are assumed to be individual level.
! Agent i's private utility from choosing one is BZ[i]+BX[i],
! where BZ[i] is the product of i's aggregate variables and their
! coefficients, and BX[i] is the product of i's individual
! variables and their coefficients.
! 
! The distinction between the two is that it is assumed that 
! agent i's friends have the same value of BZ, while
! his/her friends' values of BX are random, with correlation 
! RHO_X.
!---------------------------------------------------------------
  if (SMGLOBAL%NUMAGG > 0) then
     bz = matmul(X(:,1:SMGLOBAL%NUMAGG),p(intercept+1:(intercept+SMGLOBAL%NUMAGG)))
  else
     bz=0.0_dp
  end if
  bx=p(intercept) + matmul(X(:,(1+SMGLOBAL%NUMAGG):size(X,2)),p((intercept+SMGLOBAL%NUMAGG+1):size(p))) 
!---------------------------------------------------------------
! The procedure CHECK_CONSTRAINTS checks to see if the parameter
! values satisfy the various constraints, and subtracts a penalty
! from LOGLIKELIHOOD if they are not.
!---------------------------------------------------------------
  call check_constraints(like,rhox,rhoe,gam,bx,bz)
!---------------------------------------------------------------
! MU is the average of BX in the population, and SIGMASQ is
! its variance.   
!---------------------------------------------------------------
  mu=sum(bx)/real(size(bx,1),kind=DP)
  sigmasq=sum((bx-mu)**2)/real(size(bx,1),kind=DP)
!---------------------------------------------------------------
! Now, given parameter values the distribution of private 
! preferences is multivariate normal with mean 
!     BZ[i] + MU*(1-RHO)+RHO*BX[i]
! and a covariance matrix which is calculated in the function
! MAKEW.  In order to simulate a multivariate normal random
! vector with mean zero and covariance matrix W, you take a 
! vector of independent normal random variables and multiply
! it by the Cholesky decomposition of W.  Because the covariance
! matrix does not vary across observations, we calculate it 
! here directly, and store its Cholesky decomposition in
! the matrix C.
!---------------------------------------------------------------
  c = transpose(chol(makew(sigmasq,rhox,rhoe,gam)))
!---------------------------------------------------------------
! The two FORALL constructs in the next few lines are where the
! program spends 99.99% of its CPU time, and where any parallel
! processing is done.  For each observation i, we calculate 
! the probability of observing that particular value of YINT[i]
! (the binary choice of agent i) and YGINT[i] (the number of 
! i's friends that choose one), conditional on the assumed
! parameter vector and the value of the observed exogenous
! variables.  This probability is calculated in the function 
! SMPROB, and stored in the vector LIKEVEC.
!---------------------------------------------------------------
  if (SMGLOBAL%SIMULATOR_TYPE=="G") then
!HPF$ INDEPENDENT
     forall (i=1:size(likevec)) 
        likevec(i)=smprob_ghk(YINT(i),YGINT(i),bx(i),bz(i),mu,rhox,gam,1+NFRIENDS(i),c)
     end forall
  else
!HPF$ INDEPENDENT
     forall (i=1:size(likevec)) 
        likevec(i)=smprob_hybrid(YINT(i),YGINT(i),bx(i),bz(i),mu,rhox,gam,1+NFRIENDS(i),c)
     end forall
  end if
  if (SMGLOBAL%UNDERREPORTING_CORRECTION) then
     likevec= likevec*(1.0_dp + y*(REPORTING_RATE-1.0_dp))
     underreporting_rate=1.0_dp-REPORTING_RATE
     if (SMGLOBAL%SIMULATOR_TYPE=="G") then
!HPF$ INDEPENDENT
        forall (i=1:size(likevec))
           likevec(i)=likevec(i)+ (underreporting_rate(i))*smprob_ghk(YINT(i)+1, &
                YGINT(i),bx(i),bz(i),mu,rhox,gam,1+NFRIENDS(i),c)
        end forall
     else
!HPF$ INDEPENDENT
        forall (i=1:size(likevec))
           likevec(i)=likevec(i)+ (underreporting_rate(i))*smprob_hybrid(YINT(i)+1, &
                YGINT(i),bx(i),bz(i),mu,rhox,gam,1+NFRIENDS(i),c)
        end forall
     end if
  end if
!---------------------------------------------------------------
! Having calculated the likelihood of each observation, we 
! take the natural log and add across observations to get
! the log-likelihood.
!---------------------------------------------------------------
  like=like+sum(log(likevec))
!---------------------------------------------------------------
! Then we write the result to the log file
!---------------------------------------------------------------
  open (unit=1,file=SMGLOBAL%LOGFILE,iostat=ios,action="write",position="append",status="old")
  if (ios == 0) then
     write (unit=1,iostat=ios,fmt="(7f12.3)") like,p
     close (unit=1,iostat=ios)
  else
     continue
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





pure function smprob_hybrid(y,ygroup,bx,bz,mu,rho,gam,groupsize,c) result(smprobout)
  integer, intent(in) :: y,ygroup,groupsize
  real(kind=DP), intent(in) :: bx,bz,mu,rho,gam
  real(kind=DP), dimension(:,:), intent(in) :: c
  real(kind=DP) :: smprobout
  real(kind=DP), dimension(groupsize) :: vmu
  real(kind=DP), dimension(size(c,1),size(c,2)) :: ctmp
  real(kind=DP), dimension(groupsize) :: yloc,ypeer,mincutoff,maxcutoff,cj
  real(kind=DP), dimension(groupsize,SMGLOBAL%NSIM)  :: tt,cutofftmp1,cutofftmp2
  real(kind=DP), dimension(SMGLOBAL%NSIM) :: selected,ta,tb,wgt
  real(kind=DP), parameter :: bignum=1.0e3_dp
  integer :: i,j,ytot
  integer, dimension(groupsize,SMGLOBAL%NSIM) :: cutoff
  logical, dimension(groupsize+1,SMGLOBAL%NSIM) :: equilibrium
  if (y > 1) then
     smprobout=0.0_dp
     return
  end if
  ctmp=real(groupsize-1,kind=DP)*c
  vmu=(bz+mu*(1-rho)+rho*bx)
  vmu(1)=bz+bx
  vmu = vmu*(real(1-groupsize,kind=DP)/gam)
  yloc = 0.0_dp
  yloc(1) = real(y,kind=DP)
  if (ygroup > 0) then
     yloc(groupsize+1-ygroup:groupsize) = 1.0_dp
  end if
  ypeer = sum(yloc)-yloc
  where (yloc < 0.1_dp)
     mincutoff=ypeer
     maxcutoff=bignum
  elsewhere
     mincutoff=-bignum
     maxcutoff=ypeer
  end where
  ta = cdfn((mincutoff(1)-vmu(1))/ctmp(1,1)) ! formerly cdfn(((mine(1))/c(1,1)))
  tb = cdfn((maxcutoff(1)-vmu(1))/ctmp(1,1)) ! formerly cdfn(((maxe(1))/c(1,1)))
  wgt = tb-ta
  do j=2,groupsize
     cj=ctmp(j,1:groupsize) 
     tt((j-1),:) = cdfinvn(tb+U((j-1),:)*(ta-tb))
     ta = cdfn(((mincutoff(j)-vmu(j))-matmul(cj(1:(j-1)),tt(1:(j-1),:)))/ctmp(j,j))
     tb = cdfn(((maxcutoff(j)-vmu(j))-matmul(cj(1:(j-1)),tt(1:(j-1),:)))/ctmp(j,j))
     wgt = wgt*(tb-ta)
  end do
  smprobout = real(nchoosek(groupsize-1,ygroup),kind=DP)*sum(wgt)/real(size(wgt),DP)
  tt(groupsize,:) = cdfinvn(tb+U((groupsize),:)*(ta-tb))
  cutofftmp1=matmul(ctmp(1:groupsize,1:groupsize),tt)
  cutofftmp2=spread(vmu,dim=2,ncopies=size(tt,2))
  cutoff=floor(cutofftmp1+cutofftmp2)
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
     if (ytot < groupsize) then
        do j=ytot+2,(groupsize+1)
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
  smprobout=max(smprobout*(sum(selected)/real(size(selected),kind=DP)),1.0e-50_dp)
end function smprob_hybrid


pure function smprob_ghk(y,ygroup,bx,bz,mu,rho,gam,groupsize,c) result(smprobout)
  integer, intent(in) :: y,ygroup,groupsize
  real(kind=DP), intent(in) :: bx,bz,mu,rho,gam
  real(kind=DP), dimension(:,:), intent(in) :: c
  real(kind=DP) :: smprobout
  real(kind=DP), dimension(groupsize) :: vmu
  real(kind=DP), dimension(BIGMATSIZE(groupsize)) :: like
  real(kind=DP), dimension(BIGMATSIZE(groupsize),2*groupsize+1) :: bigmat
  real(kind=DP), dimension(size(c,1),size(c,2)) :: ctmp
  integer :: j
!---------------------------------------------------------------
! This is a trick related to the underreporting correction.
!---------------------------------------------------------------
  if (y > 1) then
     smprobout=0.0_dp
     return
  end if
!---------------------------------------------------------------
! Finish calculating the covariance matrix.  This is a trick.
!---------------------------------------------------------------
  ctmp=real(groupsize-1,kind=DP)*c
!---------------------------------------------------------------
! Each person's personal preference for y=1 can be described
! by a random variable V[i] where person i will choose y=1
! if and only if at least V[i] peers choose the same.
! The vector V is distributed multivariate normal with mean
! vector VMU and a covariance matrix whose Cholesky decomposition
! is given by CTMP. 
!---------------------------------------------------------------
  vmu=(bz+mu*(1-rho)+rho*bx)
  vmu(1)=bz+bx
  vmu = vmu*(real(1-groupsize,kind=DP)/gam)
!---------------------------------------------------------------
! The matrix BIGMAT describes the set of values for V that are
! associated with Y and YGROUP being the equilibrium behavior
! of the individual and his/her peer group. 
!
! EXPLANATION????????
!
!
!---------------------------------------------------------------
  bigmat = real(eqregion(y,ygroup,groupsize),kind=DP)
!---------------------------------------------------------------
! For each row in BIGMAT, we need to calculate the probability
! that V will be in the region described by that row.
!---------------------------------------------------------------
  do j=1,size(like)
     if (bigmat(j,1)==0.0_dp) then
        like(j) = 0.0_dp
     else 
        like(j) = bigmat(j,1)*ghk(vmu,ctmp,bigmat(j,2:(groupsize+1)), &
          bigmat(j,(groupsize+2):size(bigmat,2)))
     end if
  end do
!---------------------------------------------------------------
! Then finally, we add up these probabilities.  Because we will
! be taking logs in LOGLIKELIHOOD, we put a floor on the probability
! of a very small positive number (1.0d-50).  This keeps us from
! getting overflow errors.
!---------------------------------------------------------------
  smprobout=max(real(nchoosek(groupsize-1,ygroup),kind=DP)*sum(like),1.0e-50_dp)
end function smprob_ghk





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
  real(kind=DP), parameter :: penalty=50.0_dp,maxbx=5.0_dp, penalty2=50.0_dp
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
! EQREGION function (PRIVATE, used by SMPROB and thus LOGLIKELIHOOD)
!
! Format: eqregion(y,ygroup,groupsize)
!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
pure function eqregion(y,ygroup,groupsize) result(eqregionout)
  integer, intent(in) :: y,ygroup,groupsize
  integer, dimension(BIGMATSIZE(groupsize),2*groupsize+1) :: eqregionout
  integer :: nrow
  integer, dimension(2*groupsize-3,2*groupsize+1) :: requilibria
  integer, dimension(2*groupsize) :: actual,tmp
  integer :: i,j,k,jmax
! Calculate equilibria
  requilibria=0
  actual = eqrange(y,ygroup,groupsize)
  k=2*groupsize+1
  j=0
  if (ygroup > 1) then
     do i = 0,(ygroup-2)
        tmp = rintersect(eqrange(y,i,groupsize),actual) 
        if (nonempty(tmp)) then
           j=j+1
           requilibria(j,2:k) = tmp
           requilibria(j,1) = nchoosek(ygroup,i)
        end if
     end do
  end if
  if ((y > 0) .and. (ygroup > 0)) then
     do i=0,(ygroup-1)
        tmp = rintersect(eqrange(0,i,groupsize),actual)
        if (nonempty(tmp)) then
           j=j+1
           requilibria(j,2:k)=tmp
           requilibria(j,1) = nchoosek(ygroup,i)
        end if
     end do
  end if
! Then take their union
  eqregionout(1,1) = -1
  eqregionout(1,2:size(eqregionout,2)) = eqrange(y,ygroup,groupsize)
  eqregionout(3:size(eqregionout,1),1)=0 
  eqregionout(2,:)=requilibria(1,:)
  nrow = 2
  do i=2,size(requilibria,1)
     if (requilibria(i,1)==0) then
        exit 
     end if
     nrow=nrow+1
     eqregionout(nrow,:)=requilibria(i,:)
     jmax=nrow-1
     do j=2,jmax
        tmp = rintersect(requilibria(i,2:k),eqregionout(j,2:k))
        if (nonempty(tmp)) then
           nrow=nrow+1
           eqregionout(nrow,2:k)=tmp
           eqregionout(nrow,1)=-(eqregionout(j,1)*requilibria(i,1))
        end if
     end do
  end do          
  eqregionout(:,1)=-eqregionout(:,1)
end function eqregion


pure function eqrange(y,ygroup,groupsize) result(eqrangeout)
  integer, intent(in) :: y,ygroup,groupsize
  integer, dimension(2*groupsize) :: eqrangeout
  integer, parameter :: infty=huge(1)
  integer, dimension(groupsize) :: equil,nequil
  equil(1) = y
  equil(2:(groupsize-ygroup))=0
  if (ygroup > 0) then
     equil((groupsize-ygroup+1):groupsize)=1
  end if
  nequil = 1 - equil
  eqrangeout(1:groupsize) = -infty*equil+(y+ygroup)*nequil
  eqrangeout((groupsize+1):(2*groupsize)) = (y+ygroup-1)*equil+infty*nequil
end function eqrange


pure function rintersectd(r1,r2) result(rintersectout)
  real(kind=DP), dimension(:), intent(in) :: r1,r2
  real(kind=DP), dimension(size(r1)) :: rintersectout 
  integer :: k,i
  k = size(r1)/2
  do i=1,k
     rintersectout(i) = max(r1(i),r2(i))
     rintersectout(i+k) = min(r1(i+k),r2(i+k))
  end do
end function rintersectd
pure function rintersects(r1,r2) result(rintersectout)
  real(kind=SP), dimension(:), intent(in) :: r1,r2
  real(kind=SP), dimension(size(r1)) :: rintersectout
  integer :: k,i
  k = size(r1)/2
  do i=1,k
     rintersectout(i) = max(r1(i),r2(i))
     rintersectout(i+k) = min(r1(i+k),r2(i+k))
  end do
end function rintersects
pure function rintersecti(r1,r2) result(rintersectout)
  integer, dimension(:), intent(in) :: r1,r2
  integer, dimension(size(r1)) :: rintersectout
  integer :: k,i
  k = size(r1)/2
  do i=1,k
     rintersectout(i) = max(r1(i),r2(i))
     rintersectout(i+k) = min(r1(i+k),r2(i+k))
  end do
end function rintersecti

pure function nonemptys(rect) result(nonemptyout)
  real(kind=SP), dimension(:), intent(in) :: rect
  logical :: nonemptyout
  integer :: k,i
  k=size(rect)/2
  nonemptyout = .true.
  do i=1,k
     if (rect(i) >= rect(i+k)) nonemptyout = .false.
  end do
end function nonemptys
pure function nonemptyd(rect) result(nonemptyout)
  real(kind=DP), dimension(:), intent(in) :: rect
  logical :: nonemptyout
  integer :: k,i
  k=size(rect)/2
  nonemptyout = .true.
  do i=1,k
     if (rect(i) >= rect(i+k)) nonemptyout = .false.
  end do
end function nonemptyd
pure function nonemptyi(rect) result(nonemptyout)
  integer, dimension(:), intent(in) :: rect
  logical :: nonemptyout
  integer :: k,i
  k=size(rect)/2
  nonemptyout = .true.
  do i=1,k
     if (rect(i) >= rect(i+k)) nonemptyout = .false.
  end do
end function nonemptyi

function makew(sigmasq,rho,rho2,gam) result(w)
  real(kind=DP), intent(in) :: sigmasq,rho,rho2,gam
  real(kind=DP), dimension(SMGLOBAL%MAXGROUPSIZE,SMGLOBAL%MAXGROUPSIZE) :: w
  real(kind=DP), dimension(SMGLOBAL%MAXGROUPSIZE) :: s
  integer :: i
  s=rho2 
  s(1)=1.0_dp
  w =rho*sigmasq + rho2 - sigmasq*rho**2
  do i=1,SMGLOBAL%MAXGROUPSIZE
     w(i,i) = sigmasq+1-sigmasq*rho**2
  end do
  w(1,:)=s
  w(:,1)=s
  w=w/(gam**2)
end function makew


function ghkd(mu,c,a,b,u) result(ghkout)
  real(kind=DP), dimension(:), intent(in) :: mu
  real(kind=DP), dimension(size(mu)), intent(in) :: a,b
  real(kind=DP), dimension(size(mu),size(mu)), intent(in) :: c
  real(kind=DP), dimension(:,:), intent(in) :: u
  real(kind=DP) :: ghkout
  integer :: r,m,j,ii
  real(kind=DP), dimension(size(u,2)) :: ones,ta,tb,wgt
  real(kind=DP), dimension(size(mu)) :: cj
  real(kind=DP), dimension(size(mu),size(u,2)) :: tt
  real(kind=DP), parameter :: tinynum=1.0e-100_dp
  m = Size(mu)
  r = Size(u,2)
  j = 1
  ii = 1
  ones=1.0_dp
  ta = cdfn(((a(1)-mu(1))/(c(1,1)+tinynum))*ones) 
  tb = cdfn(((b(1)-mu(1))/(c(1,1)+tinynum))*ones)
  tt=0.0_dp
  tt(1,:) = cdfinvn(u(1,:)*ta+(ones-u(1,:))*tb)
  wgt = tb-ta
  do j=2,m
     cj=c(j,:) 
     ta = cdfn(((a(j)-mu(j))*ones-matmul(cj(1:j),tt(1:j,:)))/(c(j,j)+tinynum))
     tb = cdfn(((b(j)-mu(j))*ones-matmul(cj(1:j),tt(1:j,:)))/(c(j,j)+tinynum))
     tt(j,:) = cdfinvn(u(j,:)*ta+(ones-u(j,:))*tb)
     wgt = wgt*(tb-ta)
  end do
  ghkout = sum(wgt)/real(r,kind(1.0_dp))
end function ghkd
pure function ghkdu(mu,c,a,b) result(ghkout)
  real(kind=DP), dimension(:), intent(in) :: mu,a,b
  real(kind=DP), dimension(:,:), intent(in) :: c
  real(kind=DP) :: ghkout
  integer :: j
  real(kind=DP), dimension(smglobal%nsim) :: ta,tb,wgt
  real(kind=DP), dimension(Size(a-1),smglobal%nsim)  :: tt
  real(kind=DP), dimension(size(c,1)) :: cj
  real(kind=DP) :: r
  r = real(SMGLOBAL%NSIM,DP)
  ta = cdfn(((a(1)-mu(1))/c(1,1)))
  tb = cdfn(((b(1)-mu(1))/c(1,1)))
  wgt = tb-ta
  do j=2,size(a)
     cj=c(j,:) 
     tt((j-1),:) = cdfinvn(U((j-1),:)*ta+(1.0_dp-U((j-1),:))*tb)
     ta = cdfn(((a(j)-mu(j))-matmul(cj(1:(j-1)),tt(1:(j-1),:)))/c(j,j))
     tb = cdfn(((b(j)-mu(j))-matmul(cj(1:(j-1)),tt(1:(j-1),:)))/c(j,j))
     wgt = wgt*(tb-ta)
  end do
  ghkout = sum(wgt)/r
end function ghkdu



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! DLOGLIKELIHOOD function
!
! Format: dloglikelihood(p)
!
! Gives approximate first derivative (at p) of 
! LOGLIKELIHOOD(p).
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! function dloglikelihood(p,fp) result(fout)
subroutine dloglikelihood(p,fp,fout) ! new
  real(kind=DP), dimension(:), intent(in) :: p
  real(kind=DP), intent(in) :: fp
!  real(kind=DP), dimension(size(p)) :: fout
  real(kind=DP), dimension(size(p)), intent(inout) :: fout ! new
  real(kind=DP), dimension(size(p)) :: ptmp,gradp
  real(kind=DP) :: ftmp ! new
  integer :: i
  ptmp=abs(p)
  where (ptmp < 1.0e-2_dp) ptmp = 1.0e-2_dp
!
! Changed by BK 6/24/2004 - seems like a big change.  It *looks* like I had the 
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
end subroutine dloglikelihood ! new


end module loglik







