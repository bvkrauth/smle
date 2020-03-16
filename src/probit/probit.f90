program probit
use bklib, only : setseed,dp
use bkmath, only : pdfn,cdfn,inverse
! PROBITMC use monte, only : montecarlo
implicit none
real(kind=DP), dimension(:), allocatable :: y
real(kind=DP), dimension(:,:), allocatable :: x
! PROBITMC real(kind=DP), dimension(100) :: gamma
real(kind=DP) :: gammai,reprate
logical :: underreporting_correction
character(len=12) :: resultfile
! PROBITMC integer :: i

! PROBITMC do i=1,100
! PROBITMC   print *, i
   call setseed()
! PROBITMC   call montecarlo()
   call load_data()
!   gamma(i) = estimate(y,x)
   call estimate(y,x,reprate,gammai)
! PROBITMC gamma(i)=gammai
! PROBITMC deallocate(x,y)
! PROBITMC end do
! PROBITMC print *, "Avg gamma: ", sum(gamma)/100.0_dp
! PROBITMC print *, "Std. Dev.:", sqrt(sum((gamma-sum(gamma)/100.0_dp)**2)/100.0_dp)


contains


subroutine estimate(y,x,reprate,gamma)
  real(kind=DP), dimension(:), intent(in) :: y
  real(kind=DP), dimension(:,:), intent(in) :: x
  real(kind=DP), dimension(size(x,2)) :: b,db,g
  real(kind=DP), dimension(size(x,1)) :: pdf,cdf,d,xb,gfac
  real(kind=DP), dimension(size(x,2),size(x,2)) :: h
  real(kind=DP), dimension(size(x,1),size(x,2)) :: xtmp
  real(kind=DP), intent(in) :: reprate
  real(kind=DP), intent(inout) :: gamma
  integer :: i,j,n,ios
  logical :: file_found
  real(kind=DP) :: eps=100.0_dp,tol=0.000001_dp
  n=size(y)
  h=matmul(transpose(x),x)
  call inverse(h)
  b=matmul(h,matmul(transpose(x),y))
  do i=1,1000
     xb=matmul(x,b)
     pdf=pdfn(xb)*reprate
     cdf=cdfn(xb)*reprate
     gfac=y*(pdf/cdf)-(1.0_dp-y)*(pdf/(1.0_dp-cdf))
     do j=1,size(x,2)
        g(j)=sum(x(:,j)*gfac)
     end do     
     d=pdf*((y*(pdf+(xb*cdf))/cdf**2) + ((1.0_dp-y)*(pdf-xb*(1.0_dp-cdf))/((1.0_dp-cdf)**2)))
     do j=1,size(x,2)
        xtmp(:,j)=x(:,j)*d
     end do
     h=matmul(transpose(xtmp),x)
     call inverse(h)
     db=matmul(h,g)
     b=b+db
     eps=sum(db**2)/real(size(db),kind=DP)
     if (eps < tol) then
        exit
     end if
  end do
  inquire (file=resultfile,exist=file_found)
  if (file_found) then
     open(unit=1,file=resultfile,iostat=ios,action="write",position="append",status="old")
  else
     open(unit=1,file=resultfile,iostat=ios,action="write",position="rewind",status="new")
  end if
  if (ios == 0) then
     write (unit=1,iostat=ios,fmt=*) sum(y*log(cdf)+(1.0_dp-y)*log(1.0_dp-cdf)), b
     close(unit=1,iostat=ios)
  else
     stop "Error: Unable to open result file"
  end if
  gamma=b(2)
end subroutine estimate


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
    character(len=12) :: datafile,parmfile="parm.dat  ",logfile
    character(len=80) :: toss,samptypelong
    character(len=1) :: sample_type
    logical :: file_found
    real(kind=DP), dimension(:,:), allocatable :: dat
    integer, dimension(:), allocatable :: groupid
    logical, dimension(:), allocatable :: msk
    integer :: i,j,nobs,nvar,ios
    inquire (file=parmfile,exist=file_found)
    if (.not.file_found) then 
       stop "Error: parm.dat file not found"
    end if
    open(unit=1,file=parmfile,iostat=ios,action="read",position="rewind",status="old")
    if (ios /=0) then 
       stop "Error: parm.dat file found but could not be opened"
    end if 
    read (unit=1,fmt=*,iostat=ios) toss ! TOSS is used for comment lines; the program doesn't use it
    read (unit=1,fmt=*,iostat=ios) toss
    read (unit=1,fmt=*,iostat=ios) nvar 
    read (unit=1,fmt=*,iostat=ios) toss
    read (unit=1,fmt=*,iostat=ios) nobs 
    read (unit=1,fmt=*,iostat=ios) toss
    read (unit=1,fmt=*,iostat=ios) samptypelong
    read (unit=1,fmt=*,iostat=ios) toss
    read (unit=1,fmt=*,iostat=ios) underreporting_correction 
    read (unit=1,fmt=*,iostat=ios) toss
    read (unit=1,fmt=*,iostat=ios) datafile
    read (unit=1,fmt=*,iostat=ios) toss
    read (unit=1,fmt=*,iostat=ios) logfile 
    read (unit=1,fmt=*,iostat=ios) toss
    read (unit=1,fmt=*,iostat=ios) resultfile 
    close (unit=1,iostat=ios)
    sample_type=uppercase(samptypelong(1:1))
    allocate(x(nobs,nvar+2),y(nobs))
    open(unit=1,file=logfile,iostat=ios,action="write",position="rewind",status="replace")
    write(unit=1,iostat=ios,fmt=*) "Probit estimation"
    write(unit=1,iostat=ios,fmt=*) "Number of variables:",nvar 
    write(unit=1,iostat=ios,fmt=*) "Number of observations:",nobs
    write(unit=1,iostat=ios,fmt=*) "Sample type:",samptypelong
    write(unit=1,iostat=ios,fmt=*) "Underreporting corrrection?:",underreporting_correction
    write(unit=1,iostat=ios,fmt=*) "Data read from file:",datafile
    write(unit=1,iostat=ios,fmt=*) "Log file:",logfile
    write(unit=1,iostat=ios,fmt=*) "Results written to file:",resultfile
    close(unit=1,iostat=ios)
    if (sample_type == "I") then
       allocate(dat(nobs,nvar+3))
    elseif (sample_type=="G") then
       allocate(dat(nobs,nvar+2),groupid(nobs),msk(nobs))
    else
       stop "Illegal sample type"
    end if
    open (unit=1,file=datafile,iostat=ios,form="formatted",action="read",position="rewind",status="old")    
    if (ios == 0) then
       do i=1,nobs
          read (unit=1,fmt=*,iostat=ios) dat(i,:)
          if (ios /= 0) then
             stop "Error: not enough observations in data file"
          end if
       end do
    else
       stop "Cannot open data file"
    end if
    close(unit=1,iostat=ios)
    x(:,1)=1.0_dp
    x(:,3:size(x,2))=dat(:,(3+size(dat,2)-size(x,2)):size(dat,2))
    if (sample_type == "I") then
       y=dat(:,1)
       x(:,2)=dat(:,2)
    else
       y=dat(:,2)
       groupid=nint(dat(:,1))
       do i=1,nobs
          j=groupid(i)
          msk=(groupid==j)
          msk(i)=.false.
          x(i,2)=sum(y,msk)/real(count(msk),kind=DP)
       end do
    end if
    deallocate(dat)
    if (sample_type == "G") then
       deallocate(msk,groupid)
    end if
    if (underreporting_correction) then
       reprate=sum(y)/sum(x(:,2))
       if (reprate > 1.0_dp) then
          reprate = 1.0_dp
       end if
    else
       reprate=1.0_dp
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

end program probit
