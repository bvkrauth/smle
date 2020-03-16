@rem This file is used for compiling the SMLE program on an MS Windows system using the F compiler
@rem F is a subset of Fortran 90.  A free F compiler is available from the Fortran company
@rem at http://www.fortran.com
@echo Compiling...
g95 -Wall -O3 -o smle.exe ../lib/bklib.f90 ../lib/bkmath.f90 smglob.f90 ../lib/dfpmin.f90 ../lib/simann.f90 loglik.f90 smutil.f90 smle.f90
@rem del *.mod
@rem \Perl\bin\perl.exe makemcrep.pl
@rem F -O3 -o smlerep.exe -tempdir . ../lib/bklib.f90 ../lib/bkmath.f90 smglob.f90 ../lib/dfpmin.f90 ../lib/simann.f90 loglik.f90 smutil.f90 smlerep.f90
@rem del *.mod smlerep.f90
@rem F -O3 -o smc.exe -tempdir . ../lib/bklib.f90 ../lib/bkmath.f90 smglob.f90 ../lib/dfpmin.f90 ../lib/simann.f90 loglik.f90 ../psim/monte.f90 smutil.f90 smc.f90
@rem del *.mod smc.f90 
