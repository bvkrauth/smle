g95  -o probit ../lib/bklib.f90 ../lib/bkmath.f90 probit.f90
del *.mod *.o
\progfile\perl\bin\perl.exe makemc.pl
g95 -o probitmc ../lib/bklib.f90 ../lib/bkmath.f90 ../psim/monte.f90 probitmc.f90
del *.mod *.o probitmc.f90

