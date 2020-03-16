!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! PSIM PROGRAM: MONTE module 
! Author: Brian Krauth, Simon Fraser University
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module monte
  use bklib, only : dp,sp,runif ! ,runifd  ! machine-specific code
  use bkmath, only : cdfinvn,chol ! useful math functions
  implicit none
  private
  private :: uppercase
  public :: montecarlo, probit_simulate

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  VARIABLE DECLARATIONS
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! none

contains


subroutine montecarlo()
  character(len=12) :: parmfile="parmonte.dat",eqtypelong,xtypelong
  character(len=80) :: toss
  character(len=1) :: eqtype,xtype 
  logical :: file_found
  integer :: ngroup,maxgroupsize,nvar,numagg,ios
  real(kind=DP), dimension(4) :: reporting_rate
  real(kind=DP), dimension(:), allocatable :: b
  real(kind=DP), dimension(2) :: b2
  b2=0.0_dp
  inquire (file=parmfile,exist=file_found)
  if (.not.file_found) then 
     stop "Error: parmonte.dat file not found"
  end if
  open(unit=1,file=parmfile,iostat=ios,action="read",position="rewind",status="old")
  if (ios /= 0) then
     stop "Error: parmonte.dat file found, but could not be opened"
  end if 
  read (unit=1,fmt=*,iostat=ios) toss ! TOSS is used for comment lines; the program doesn't use it
  read (unit=1,fmt=*,iostat=ios) toss
  read (unit=1,fmt=*,iostat=ios) ngroup
  read (unit=1,fmt=*,iostat=ios) toss
  read (unit=1,fmt=*,iostat=ios) maxgroupsize
  read (unit=1,fmt=*,iostat=ios) toss
  read (unit=1,fmt=*,iostat=ios) nvar
  read (unit=1,fmt=*,iostat=ios) toss
  read (unit=1,fmt=*,iostat=ios) numagg
  read (unit=1,fmt=*,iostat=ios) toss
  read (unit=1,fmt=*,iostat=ios) eqtypelong  
  read (unit=1,fmt=*,iostat=ios) toss
  read (unit=1,fmt=*,iostat=ios) xtypelong  
  allocate(b(nvar+5)) 
  read (unit=1,fmt=*,iostat=ios) toss
  read (unit=1,fmt=*,iostat=ios) b
  read (unit=1,fmt=*,iostat=ios) toss
  read (unit=1,fmt=*,iostat=ios) reporting_rate
  read (unit=1,fmt=*,iostat=ios) toss
  read (unit=1,fmt=*,iostat=ios) b2
  close (unit=1,iostat=ios)
  eqtype=uppercase(eqtypelong(1:1))
  xtype=uppercase(xtypelong(1:1))
  if (ngroup < 1) then
    stop "Error in parmonte.dat: ngroup < 1 not allowed"
  elseif (maxgroupsize < 2) then 
    stop "Error in parmonte.dat: maxgroupsize < 2 not allowed"
  elseif (nvar < 1) then
    stop "Error in parmonte.dat: nvar < 1 not allowed"
  elseif (numagg < 0) then
    stop "Error in parmonte.dat: numagg < 0 not allowed"
  elseif ((numagg+1) > nvar) then
    stop "Error in parmonte.dat: numagg must be less than nvar"
  elseif ((eqtype /= "L").and.(eqtype /= "H").and.(eqtype /= "R")) then
    stop "Error in parmonte.dat: specified eqtype not allowed"
  elseif ((xtype /= "N").and.(xtype /= "B")) then
    stop "Error in parmonte.dat: specified xtype not allowed"
  elseif (b(1) < -1.0_dp/real(maxgroupsize-1,kind=DP)) then 
    stop "Error in parmonte.dat: rho_x too low"
  elseif (b(2) < -1.0_dp/real(maxgroupsize-1,kind=DP)) then 
    stop "Error in parmonte.dat: rho_e too low"
  elseif (b(3) < 0.0_dp) then
    stop "Error in parmonte.dat: gamma < 0 not allowed"
! temporarily deleted
!  elseif (reporting_rate < 0.0_dp) then
!    stop "Error in parmonte.dat: reporting rate < 0 not allowed"
!  elseif (reporting_rate > 1.0_dp) then 
!    stop "Error in parmonte.dat; reporting rate > 1 not allowed"
  end if
write (*,*) "Done loading settings"
  call probit_simulate(ngroup,maxgroupsize,numagg,eqtype,xtype,b,reporting_rate,b2) 
  deallocate(b)
end subroutine montecarlo



function uppercase(lower) result (upper) 
  character(len=1), intent(in)  :: lower
  character(len=1) :: upper
  if ((ichar(lower) >= ichar("a")) .and. (ichar(lower) <= ichar("z"))) then 
     upper = char(ichar("A")+ichar(lower)-ichar("a"))
  else
     upper = lower
  end if
end function uppercase

  




!---------------------------------------------------------------
! PROBIT_SIMULATE subroutine
! 
! Purpose: Creates a simulated data set.  Writes the results to
! two files:
!    ONEOBS.TXT:  Contains a single observation from each group,
!                 to replicate each random sample.
!    ALLOBS.TXT:  Contains all observations, to replicate a group
!                 based sample.
!
! Inputs:
!    NGROUP       Number of peer groups to simulate (integer)
!    MAXGROUPSIZE Size of each group (integer)
!    B            Vector of parameters (5-vector of DP reals)
!                 B(1)=correlation in x's
!                 B(2)=correlation in u's
!                 B(3)=peer effect parameter
!                 B(4)=contextual effect
!                 B(5)=intercept
!                 B(6:...)=coefficient on X
! Outputs:
!    None.  Creates the two files described above
!---------------------------------------------------------------
subroutine probit_simulate(ngroup,maxgroupsize,numagg,eqtype,xtype,b,reporting_rate,b2)
  integer, intent(in) :: ngroup,maxgroupsize,numagg
  real(kind=DP), dimension(:), intent(in) :: b,b2
  real(kind=DP), dimension(4), intent(in) :: reporting_rate
  character(len=1), intent(in) :: eqtype,xtype 
  ! Changed 10/2008 to solve a stack overflow problem.  These arrays can turn out
  ! to be quite big.  When they are declared locally, the memory used to store
  ! them will come from the stack.  But there is a lot less stack memory
  ! than heap memory.  So the program can handle much larger data if we make
  ! them allocatable instead (since the memory will come from the heap).
  ! There is probably some better way to do this, but this will do.
  real(kind=DP), dimension(:,:), allocatable :: u,bx,p1,p2,equil5,h
  real(kind=DP), dimension(:,:,:), allocatable :: x
  integer, dimension(:,:), allocatable :: y,cutoff,cutoff2,equil1,equil2,equil3,equil4
!  real(kind=DP), dimension(ngroup,maxgroupsize) :: u,bx,p1
!  real(kind=DP), dimension(ngroup,maxgroupsize,size(b)-5)  :: x 
!  real(kind=DP), dimension(ngroup,size(b)-5) :: p2
!  integer, dimension(ngroup,maxgroupsize) :: y,cutoff,cutoff2
!  integer, dimension(ngroup,maxgroupsize+1) :: equil1,equil2,equil3,equil4
!  real(kind=DP), dimension(ngroup,maxgroupsize+1) :: equil5
  real(kind=DP), dimension(1) :: equil6
!  real(kind=DP), dimension(maxgroupsize,maxgroupsize) :: h 
  real(kind=DP) :: mu,rhox,rhoe,gam,lambda,eps,k,tmp1,tmp2,sigmax,rxe,sxe
  real(kind=DP), parameter :: about_zero=1.0e-20_dp
  integer :: i,j,ios
  ! Now we need to allocate all the big variables.
  allocate( u(ngroup,maxgroupsize), bx(ngroup,maxgroupsize), p1(ngroup,maxgroupsize), &
       x(ngroup,maxgroupsize,size(b)-5), p2(ngroup,size(b)-5), y(ngroup,maxgroupsize), &
       cutoff(ngroup,maxgroupsize), cutoff2(ngroup,maxgroupsize), equil1(ngroup,maxgroupsize+1), &
       equil2(ngroup,maxgroupsize+1), equil3(ngroup,maxgroupsize+1), equil4(ngroup,maxgroupsize+1), &
       equil5(ngroup,maxgroupsize+1), h(maxgroupsize,maxgroupsize) ,stat=ios)
  if (ios /= 0) then
     write (*,*) "Error allocating variables in PROBIT_SIMULATE, ios'=",ios
  end if
!
! Step 1: First we create the x and u matrices
!
  print *, "loading data"
!
! We start with IID N(0,1) random variables in u and x
!
  call runif(p1)
  call runif(p2)
  call runif(u)
! Warning: The function call below will give a segmentation fault if u is too big.
! This is because my CDFINVN function has the same problem we fixed above 
! (trying to make a big array from stack memory), and 
! is the main constraint on the size of data set we can generate.
  u=cdfinvn(u) 
  do i=1,size(x,3)
!     x(:,:,i)=runifd(size(x,1),size(x,2))
     call runif(x(:,:,i)) ! new
  end do
!
! If there are aggregate variables, make them identical
! across all members of the same group
!
  if (numagg > 0) then                                
     do i = 1,numagg
        do j=1,ngroup
           x(j,:,i)=x(j,1,i)
        end do
     end do
  end if
!
! Break up coefficient vector B into its components
!
  rhox = b(1)
  rhoe = b(2)
  gam = b(3)
  lambda = b(4)
  sxe = b2(1)
  rxe = b2(2)
!
! Step 2: The parameter vector B indicates a particular
! covariance structure for both x and u, so we transform
! the matrices to have that covariance.
!
! This makes use of the handy trick: If X~N(0,I) then MU + Cholesky(H)*X ~ N(MU,H)
!
  h=rhox                                             
  do i=1,size(h,1)                                   
     h(i,i)=1.0_dp                                   
  end do
  h=chol(h)
  if (xtype == "N") then  
     do i=1,size(x,3) 
        x(:,:,i) = cdfinvn(x(:,:,i)) 
     end do 
     do i=(1+numagg),size(x,3)
        x(:,:,i)=matmul(x(:,:,i),h)
     end do
  elseif (xtype == "B") then 
     eps = sqrt(rhox)/2.0_dp 
     where (p2 < 0.5_dp)  
        p2 = 0.5_dp + eps 
     elsewhere
        p2 = 0.5_dp - eps 
     end where 
     do i=1,size(p2,1) 
        do j=1,size(p2,2)
           where (x(i,:,j) < p2(i,j)) 
              x(i,:,j) = 1.0_dp 
           elsewhere 
              x(i,:,j) = 0.0_dp 
           end where
        end do
     end do
  else 
    stop "Illegal xtype" 
  end if 
! Create vector BX, describes preferences
  bx = b(5)
  do i=1,size(x,3)
     bx=bx+b(i+5)*x(:,:,i) 
  end do
!
! Now we will do the same for u, reusing the H matrix.  
!
  k=real(maxgroupsize-1,kind=DP)
  sigmax = sum(bx**2)/real(size(bx,1)*size(bx,2),kind=DP)-(sum(bx)/real(size(bx,1)*size(bx,2),kind=DP))**2
  tmp1 = rhoe + ((((k-1.0_dp)-k*rhox)*rxe**2)+2.0_dp*rxe*sxe-rhox*sxe)/(sigmax*(k*rhox**2-(k-1.0_dp)*rhox-1.0_dp))
  h = tmp1
  tmp2 = 1.0_dp + (k*rxe-2.0_dp*k*rhox*rxe*sxe+sxe**2+(k-1.0_dp)*rhox*sxe**2)/(sigmax*(k*rhox**2-(k-1.0_dp)*rhox-1.0_dp))
  do i=1,size(h,1)
     h(i,i)=tmp2
  end do
  h=chol(h)
  u=matmul(u,h)
! fiddle with u
  if ((b2(1) > 0.0_dp).or.(b2(2) > 0.0_dp)) then
     tmp1 = (k*rhox*rxe-sxe-(k-1.0_dp)*rhox*sxe)/(sigmax*(k*rhox**2-(k-1.0_dp)*rhox-1.0_dp))
     tmp2 = (rhox*sxe-rxe)/(sigmax*(k*rhox**2-(k-1.0_dp)*rhox-1.0_dp))
     do i=1,ngroup 
        mu = sum(bx(i,:)) ! /real(maxgroupsize,kind=DP) 
        do j=1,maxgroupsize
           u(i,j)=u(i,j)+tmp1*bx(i,j)+tmp2*mu
        end do
     end do
  end if
! Add contextual effects !!!!!!!!!!!!!!!MAY HAVE HAD BUG!!!!!!!!!!!!!!!!!
  do i=1,ngroup 
     mu = lambda*sum(bx(i,:))/real(maxgroupsize,kind=DP) 
     bx(i,:)= (1.0_dp-lambda/real(maxgroupsize,kind=DP))*bx(i,:)+mu 
  end do 
  bx=bx+u
!
! Step 3: Now we solve for equilibrium choices
!
  y=0                                                ! Initialize all the variables
  equil1=0
  equil2=0
  equil3=0
  equil4=0
  if (gam > about_zero) then                             ! Don't want to divide by zero!
     bx = -real(maxgroupsize-1,kind=DP)*bx/gam       ! BX = # of friends (real) choosing y=1 at which you are indifferent
     cutoff = floor(bx)                              ! CUTOFF= max # of friends choosing y=1 at which you prefer y=0
!
! Avoid messing with anything for the next few lines,
! It's very hard to explain and very easy to screw up!
! The overall goal is to get (in EQUIL4) a matrix of zeros
! and ones such that EQUIL4(i,j) = 1 if and only if it is a Nash 
! equilibrium for exactly j-1 players in group i to choose y=1
!
     do i=1,size(equil1,2)                           
        equil1(:,i)=i-1                              
     end do                                          
     do i=2,size(equil2,2)                           
        where (cutoff < i-2)                         
           cutoff2=1               
        elsewhere
           cutoff2=0
        end where
        do j=1,size(equil2,1)
           equil2(j,i) = sum(cutoff2(j,:))
        end do
     end do
     equil2(:,1)=equil2(:,2)
     do i=2,size(equil3,2)
        where (cutoff < i-1) 
           cutoff2=1
        elsewhere
           cutoff2=0
        end where
        do j=1,size(equil3,1)
           equil3(j,i) = sum(cutoff2(j,:))
        end do
     end do
     where (equil1 == equil2) equil4=1 
     where (equil1 /= equil3) equil4=0               
!
! At this point EQUIL4 has coded the Nash equilibria.
! It has NGROUP rows and MAXGROUPSIZE+1 columns,
! EQUIL4[i,j] =    1 if it is an equilibrium in group i for j-1 people to choose y=1
!                  0 otherwise.
!
!
! Next we apply the equilibrium selection rule.
!
     if (eqtype == "L") then
        equil4=equil4*(equil1+1)                       
        where (equil4 == 0) equil4=maxgroupsize+1    
        do i=1,size(cutoff2,1)
           cutoff2(i,:)=minval(equil4(i,:))-1
        end do
     elseif (eqtype == "H") then
        equil4=equil4*(equil1+1)
        do i=1,size(cutoff2,1)
           cutoff2(i,:)=maxval(equil4(i,:))-1
        end do
     elseif (eqtype=="R") then
        call runif(equil5)
        equil5 = real(equil4,kind=DP)*equil5
        do i=1,size(cutoff2,1)
           equil6 = maxloc(equil5(i,:))
           cutoff2(i,:)= equil6(1)-1
        end do
     else
        stop "Error: illegal equilibrium type selected"
     end if
! Now we set Y at its selected equilibrium value
     where (cutoff < cutoff2) 
        y=1
     elsewhere
        y=0
     end where
  else
     where (bx > 0.0_dp)                          ! This is in here to avoid dividing by zero when gam=0 (i.e., b(3)=0)
        y=1
     elsewhere
        y=0
     end where
  end if
!
! Step 3: Now we write the simulated data to a file
!
  print *, "Saving"
  open(unit=1,file="allobs.txt",iostat=ios,action="write",position="rewind",status="replace")
  do i=1,size(x,1)
     do j=1,size(x,2)
        write(unit=1,fmt=*,iostat=ios) i, y(i,j), x(i,j,:)
     end do
  end do
  close(unit=1,iostat=ios)
  open(unit=1,file="oneobs.txt",iostat=ios,action="write",position="rewind",status="replace")
  where (y(:,1) < 0.5_dp)
     where (p1(:,1) < reporting_rate(1))
        y(:,1) = 1.0_dp
     end where
  elsewhere
     where (p1(:,1) > reporting_rate(2))
        y(:,1) = 0.0_dp
     end where
  end where
  where (y(:,2:size(y,2)) < 0.5_dp)
     where (p1(:,2:size(y,2)) < reporting_rate(3))
        y(:,2:size(y,2)) = 1.0_dp
     end where
  elsewhere
     where (p1(:,2:size(y,2)) > reporting_rate(4))
        y(:,2:size(y,2)) = 0.0_dp
     end where
  end where
  do i=1,size(x,1)
     write (unit=1,fmt=*,iostat=ios) y(i,1), real(sum(y(i,:))-y(i,1),kind=SP)/real(maxgroupsize-1,kind=SP), maxgroupsize-1, x(i,1,:)
  end do
  close(unit=1,iostat=ios)
 deallocate( u,bx,p1,x,p2,y,cutoff,cutoff2,equil1,equil2,equil3,equil4,equil5,h)
end subroutine probit_simulate


end module monte










