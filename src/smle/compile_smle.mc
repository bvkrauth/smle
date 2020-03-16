perl makemc.pl
hpf -O2 -Mautopar -o smc -llapack -lblas -lnag ../lib/bklibbug.f90 smglob.f90 ../lib/bkmath.f90 loglik.f90 ../lib/dfpmin.f90 ../lib/simann.f90 ../psim/monte.f90 smutil.f90 smc.f90 
rm *.o *.mod smc.f90
# 

