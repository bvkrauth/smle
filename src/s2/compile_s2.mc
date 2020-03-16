perl makemc.pl
hpf -fast -o s2mc -llapack -lblas -lnag ../lib/bklibbug.f90 ../lib/bkmath.f90 s2glob.f90 loglik.f90 ../lib/dfpmin.f90 ../lib/simann.f90 smle2.f90 ../psim/monte.f90 s2mc.f90
rm *.o *.mod *~ s2mc.f90 
# 

