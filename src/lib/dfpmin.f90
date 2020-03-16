!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! DFPMIN Module
! Author: Brian Krauth, Simon Fraser University
!
! Optimization code for SMLE program.  
!
! Public procedures defined in this module are:
!
!   dfp:	Optimize a function using the method of Davidson,
!		Fletcher and Powell.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module dfpmin
use bklib, only : SP,DP
use smglob, only : save_dfp_checkpoint, SMGLOBAL, RESUME_FROM_CHECKPOINT, load_dfp_checkpoint
implicit none
private
private :: dfpd,outerprod,outerprods,outerprodd,vabs,vabss,vabsd, &
     unit_matrix,unit_matrixs,unit_matrixd,linmib
public :: dfp 
interface dfp
   module procedure dfpd
end interface
interface outerprod
   module procedure outerprods,outerprodd
end interface
interface vabs
   module procedure vabss,vabsd
end interface
interface unit_matrix
   module procedure unit_matrixs,unit_matrixd
end interface


contains


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! DFP subroutine
!
! Format: dfp(p,gtol,iter,fret,func,dfunc)
!
! Finds the minimum of function FUNC, given a starting value
! P.  Returns the minimizer in P, the minimum in FRET, and the
! number of iterations required in ITER.  GTOL is the convergence
! tolerance, and DFUNC is a function that returns the first 
! derivative (gradient) of FUNC.
! 
! The algorithm used is the Davidson-Fletcher-Powell algorithm
! and the code is adapted from "Numerical recipes in Fortran 90"
! See that book for detailed documentation.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE dfpd(p,gtol,iter,fret,func,dfunc)
integer, intent(out) :: iter
real(kind=DP), intent(in) :: gtol
real(kind=DP), intent(out) :: fret
real(kind=DP), dimension(:), intent(inout) :: p
interface
!   function func(p) result (fout)
   subroutine func(p,fout) ! new 
     use bklib
     use smglob
     real(kind=dp), dimension(:), intent(in) :: p
!     real(kind=dp) :: fout
     real(kind=dp), intent(inout) :: fout ! new
!   end function func
   end subroutine func ! new
!   function dfunc(p,fp) result(fout)
   subroutine dfunc(p,fp,fout) ! new
     use bklib
     use smglob
     real(kind=dp), dimension(:), intent(in) :: p
     real(kind=dp), intent(in) :: fp
!     real(kind=dp), dimension(size(p)) :: fout
     real(kind=dp), dimension(size(p)), intent(inout) :: fout ! new
!   END FUNCTION dfunc
   end subroutine dfunc ! new
END INTERFACE
! constants differ from gauss
INTEGER, PARAMETER  :: ITMAX=200
! REAL(kind=DP), PARAMETER :: STPMX=100.0_dp,EPS=0.0000000001_dp,ftol=0.00001_dp
REAL(kind=DP), PARAMETER :: STPMX=5.0_dp,EPS=0.0000000001_dp,ftol=0.00001_dp
! description: does BFGS search
INTEGER :: its,itstart,ios
!LOGICAL :: check
REAL(kind=DP) :: den,fac,fad,fae,fp,stpmax,sumdg,sumxi
REAL(kind=DP), DIMENSION(size(p)) :: dg,g,hdg,pnew,xi
REAL(kind=DP), DIMENSION(size(p),size(p)) :: hessin
!
! Based on the function DFPMIN  in Numerical Recipes in Fortran 90, pages 1215-1216
!
! Given a starting point P that is a vector of length N, the Broyden-Fletcher-Goldfarb-Shanno
! variant of the Davidon-Fletcher-Powell algorithm is performed on a function FUNC, using
! its gradient as calculated by the routine DFUNC. The convergence requirement on zeroing
! the gradient is input as GTOL.  Returned quantities are P (the location of the minumum)
! ITER (the number of iterations that were performed) and FRET (the minimum value of the 
! function).  ITMAX is the maximum allowed number of iterations; STPMX is the scaled
! maximum step length allowed in line searches; EPS is the machine precision; TOLX is
! the convergence criterion on x values
!
if (RESUME_FROM_CHECKPOINT) then  
   call load_dfp_checkpoint(itstart,p,fp,g,hessin)
else 
   itstart=1 
!   fp = func(p) ! Calculate starting value and gradient
   call func(p,fp) ! Calculate starting value and gradient ! new
!   g = dfunc(p,fp)
   call dfunc(p,fp,g)
   call unit_matrix(hessin) ! Initialize inverse hessian to the unit matrix
end if
! deleted by BK 11/28/03  
! xi=-g ! initial line direction, moved by BK inside loop
stpmax=STPMX*max(vabs(p),real(size(p),dp)) ! Used by lnsrch, but not linmib
open (unit=1,file=smglobal%logfile,iostat=ios,action="write",position="append",status="old")
if (ios==0) then
   write (unit=1,iostat=ios,fmt=*) "STPMX: ", STPMX
   write (unit=1,iostat=ios,fmt=*) "vabs(p): ", vabs(p)
   write (unit=1,iostat=ios,fmt=*) "real(size(p),dp): ", real(size(p),dp)
   write (unit=1,iostat=ios,fmt=*) "stpmax: ", stpmax
   close (unit=1,iostat=ios)
else
   continue ! no need for error handling here
end if
do its=itstart,ITMAX ! main search loop
  open (unit=1,file=smglobal%logfile,iostat=ios,action="write",position="append",status="old")
  if (ios==0) then
     write (unit=1,iostat=ios,fmt=*) "Iteration: ", its
     write (unit=1,iostat=ios,fmt=*) "Current optimal b: ", p
     write (unit=1,iostat=ios,fmt=*) "Current value of function: ",-fp
     write (unit=1,iostat=ios,fmt=*) "g: ", g
     close (unit=1,iostat=ios)
  else
     continue ! no need for error handling here
  end if
  call save_dfp_checkpoint(its,p,fp,g,hessin) ! save current status to disk in case program is interrupted
  iter=its
  xi=-matmul(hessin,g) ! set line direction; moved up from bottom of loop by BK 11/28/03
!  call lnsrch(p,fp,g,xi,pnew,fret,stpmax,check,func)
  call linmib(p,fp,xi,pnew,fret,func) ! use LINMIB to do linear search for minimum in direction XI
  if (2.0_dp*abs(fp-fret) < ftol*(abs(fp)+abs(fret)+eps)) then ! function has converged
     open (unit=1,file=smglobal%logfile,iostat=ios,action="write",position="append",status="old")
     if (ios == 0) then
        write (unit=1,iostat=ios,fmt=*) "relative function convergence"
        close (unit=1,iostat=ios)
     else
        continue
     end if
     return
  end if
! added by BK - only used when likelihood bounds are being used, should probably be removed
  if (fret < SMGLOBAL%DFPSTOP) then 
     open (unit=1,file=smglobal%logfile,iostat=ios,action="write",position="append",status="old")
     if (ios == 0) then
        write (unit=1,iostat=ios,fmt=*) "fret < dfpstop"
        close (unit=1,iostat=ios)
     else
        continue
     end if
     return
  end if
! back to original code
  fp=fret ! update the current minimum
!   xi=pnew-p ! direction (THIS APPEARS TO BE A BUG IN THE NUMERICAL RECIPES CODE; deleted) 
  p=pnew ! and minimizer
! NOT IN GAUSS VERSION
!   if (maxval(abs(xi)*max(abs(p),1.0_dp)/den) < gtol) then ! test for convergence on dX
!      return
!   end if
  dg=g ! save the old gradient
!  g=dfunc(p,fp) ! and get the new gradient
  call dfunc(p,fp,g) ! new
  den=max(fret,1.0_dp)
!   if(maxval(abs(g)*max(abs(p),1.0_dp)/den) < gtol) then ! test for convergence on zero gradient
!      return
!   end if
  dg=g-dg ! Compute difference of gradients
  hdg=matmul(hessin,dg) ! and difference times current matrix
  fac=dot_product(dg,xi) ! calculate dot products for the denominators
  fae=dot_product(dg,hdg)
  sumdg=dot_product(dg,dg)
  sumxi=dot_product(xi,xi)
  open (unit=1,file=smglobal%logfile,iostat=ios,action="write",position="append",status="old")
  if (ios == 0) then
     write (unit=1,iostat=ios,fmt=*) "STPMX: ", sumdg,sumxi,eps,fac,fae
     close (unit=1,iostat=ios)
  else
     continue
  end if
  if (fac > sqrt(eps*sumdg*sumxi)) then ! skip update if FAC not sufficiently positive
     fac=1.0_dp/fac
     fad=1.0_dp/fae
     dg=fac*xi-fad*hdg ! vector that makes BFGS different from DFP
     hessin=hessin+fac*outerprod(xi,xi)-fad*outerprod(hdg,hdg)+fae*outerprod(dg,dg) ! the BFGS updating formula
  else
     open (unit=1,file=smglobal%logfile,iostat=ios,action="write",position="append",status="old")
     if (ios == 0) then
        write (unit=1,iostat=ios,fmt=*) "fac nearly zero, fac= ",fac," EPS=",EPS, & 
             "sumdg=",sumdg,"sumxi=",sumxi
        close (unit=1,iostat=ios)
     else
        continue
     end if
     return
  end if
! Moved to beginning by BK 11/28/03
!  xi=-matmul(hessin,g)
end do ! start the next iteration
stop "Error: too many iterations in DFP algorithm"
END SUBROUTINE dfpd

FUNCTION outerprods(a,b) result(outerprodout)
REAL(kind=SP), DIMENSION(:), INTENT(IN) :: a,b
REAL(kind=SP), DIMENSION(size(a),size(b)) :: outerprodout
outerprodout = spread(a,dim=2,ncopies=size(b))*spread(b,dim=1,ncopies=size(a))
END FUNCTION outerprods
FUNCTION outerprodd(a,b) result(outerprodout)
REAL(kind=DP), DIMENSION(:), INTENT(IN) :: a,b
REAL(kind=DP), DIMENSION(size(a),size(b)) :: outerprodout
outerprodout = spread(a,dim=2,ncopies=size(b))*spread(b,dim=1,ncopies=size(a))
END FUNCTION outerprodd

FUNCTION vabss(v) result(vabsout)
REAL(kind=SP), DIMENSION(:), INTENT(IN) :: v
REAL(kind=SP) :: vabsout
vabsout=sqrt(dot_product(v,v))
END FUNCTION vabss
FUNCTION vabsd(v) result(vabsout)
REAL(kind=DP), DIMENSION(:), INTENT(IN) :: v
REAL(kind=DP) :: vabsout
vabsout=sqrt(dot_product(v,v))
END FUNCTION vabsd

SUBROUTINE unit_matrixs(mat) 
REAL(kind=SP), DIMENSION(:,:), INTENT(OUT) :: mat
INTEGER :: i,n
n=min(size(mat,1),size(mat,2))
mat(:,:)=0.0_sp
do i=1,n
   mat(i,i)=1.0_sp
end do
END SUBROUTINE unit_matrixs
SUBROUTINE unit_matrixd(mat) 
REAL(kind=DP), DIMENSION(:,:), INTENT(OUT) :: mat
INTEGER :: i,n
n=min(size(mat,1),size(mat,2))
mat(:,:)=0.0_dp
do i=1,n
   mat(i,i)=1.0_dp
end do
END SUBROUTINE unit_matrixd



subroutine linmib(pst,fp,xi,pnew,fret,func)
! implicit none
real(kind=DP), dimension(:), intent(in) :: pst,xi
real(kind=DP), intent(in) :: fp
real(kind=DP), dimension(:), intent(inout) :: pnew
real(kind=DP), intent(inout) :: fret
real(kind=DP) :: a1,a2,f1,f2,am,fm,ad,a,b,u,v,w,x,e,fx,xm,tol1,tol2,p,q,r,etemp, &
     d,fu,fw,fv
integer :: i
real(kind=DP), parameter :: cgold=0.381966_dp,zeps=1.0e-10_dp,tol=0.00001_dp
! real(DP), parameter :: cgold=0.381966d0,zeps=1d-5,tol=0.001d0
integer, parameter :: maxit=100
interface
!   function func(p) result(fout)
   subroutine func(p,fout) ! new
     use bklib
     use smglob
     implicit none
!     real(kind=DP) :: fout
     real(kind=DP), intent(inout) :: fout ! new
     real(kind=DP), dimension(:), intent(in) :: p
!   end function func
   end subroutine func
end interface
a1=0.0_dp
a2=0.01_dp
! new stuff - see if it works better
if (xi(1) > 0.0_dp) then
   b=(0.9_dp-pst(1))/xi(1)
   a=(-0.5_dp-pst(1))/xi(1)
else
   a=(0.9_dp-pst(1))/xi(1)
   b=(-0.5_dp-pst(1))/xi(1)  
end if
!f1 = func(pst+a*xi)
call func(pst+a*xi,f1) ! new
!f2 = func(pst+b*xi)
call func(pst+b*xi,f2) ! new
if ((f1 > fp).and.(f2 > fp)) then
   v=0.0_dp
   x=v
   w=v
   e=v
   fx=fp
   fv=fx
   fw=fx
   goto 5
end if
! end of new stuff
f1=fp
do i=1,10
   if (i==10) then
      stop "Error in linmib: function is too flat "
   end if
!   f2=func(pst+a2*xi)
   call func(pst+a2*xi,f2) ! new
   if (f1 /= f2) exit 
   a2=a2+1.0_dp
end do
if (f1 <= f2) then
   am=a1
   fm=f1
   f1=f2
   a1=a2
else
   am=a2
   fm=f2
end if
ad=5.0_dp*(am-a1)
a2=am+ad
!f2 = func(pst+a2*xi)
 call func(pst+a2*xi,f2) ! new
do i=1,maxit
   if (f2 > fm) exit
   ad=2.0_dp*ad
   a2=am+ad
!   f2=func(pst+a2*xi)
   call func(pst+a2*xi,f2) ! new
end do
if (a1 < a2) then
   a=a1
   b=a2
else
   a=a2
   b=a1
end if
v=am
x=v
w=v
e=0.0_dp
fx=fm
fv=fx
fw=fx
5 continue
do i=1,maxit
   xm=(a+b)*0.5_dp
   tol1=tol*abs(x)+zeps 
   tol2=2.0_dp*tol1
   if (abs(x-xm) <= (tol2-0.5_dp*(b-a))) goto 30
   if (abs(e) > tol1) then
      r=(x-w)*(fx-fv)
      q=(x-v)*(fx-fw)
      p=(x-v)*q-(x-w)*r
      q=2.0_dp*(q-r)
      if (q > 0) p=-p
      q=abs(q)
      etemp=e
      e=d
      if ((abs(p)>=abs(0.5_dp*q*etemp)).or.(p<=q*(a-x)).or.(p>=(q*(b-x)))) goto 10
      d=p/q
      u=x+d
      if (((u-a)<tol2).or.((b-u)<tol2)) then
         d=abs(tol1)
         if (xm<x) d=-d
      end if
      goto 20
   end if
10 continue
   if(x>=xm) then
      e=a-x
   else
      e=b-x
   end if 
   d=cgold*e
20 continue
   if (abs(d) >= tol1) then
      u=x+d
   else
      if (d>0.0_dp) then
         u=x+abs(tol1)
      else
         u=x-abs(tol1)
      end if
   end if
!   fu=func(pst+u*xi)
   call func(pst+u*xi,fu) !new
   if (fu <= fx) then
      if (u >= x) then
         a=x
      else
         b=x
      end if
      v=w
      fv=fw
      w=x
      fw=fx
      x=u
      fx=fu
   else
      if (u<x) then
         a=u
      else
         b=u
      end if
      if ((fu<= fw).or.(w==x)) then
         v=w
         fv=fw
         w=u
         fw=fu
      else 
         if ((fu <= fv).or.(v == x).or.(v == w)) then
            v=u
            fv=fu
         end if
      end if
   end if
end do
stop "Error in linmib: max number of iterations exceeded"
30 continue 
  pnew=pst+x*xi
  fret=fx
end subroutine linmib







end module dfpmin
