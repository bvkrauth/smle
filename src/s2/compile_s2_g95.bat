@rem -----------------------------------------------------
@rem  COMPILE_S2.BAT
@rem 
@rem  This batch file is used to compile the program
@rem  S2.EXE.  It requires the F compiler, but can
@rem  be modified to use any Fortran 90 compiler.
@rem  An F compiler is available free from 
@rem        http://www.fortran.com
@rem -----------------------------------------------------
@rem
g95 -Wall -O3 -o s2.exe  ../lib/bklib.f90 ../lib/bkmath.f90 s2glob.f90 loglik.f90 ../lib/dfpmin.f90 ../lib/simann.f90 smle2.f90 s2.f90
del *.mod
