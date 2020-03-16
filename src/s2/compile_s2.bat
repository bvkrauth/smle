@rem -----------------------------------------------------
@rem  COMPILE_S2.BAT
@rem 
@rem  This batch file is used to compile the program
@rem  S2.EXE.  It uses the Intel Fortran compiler, but can
@rem  be modified to use any Fortran 90 compiler.
@rem -----------------------------------------------------
@rem
@rem These are the compiler options that seems to give the best performance in tests
@rem ifort /O2 /Qipo /Qprec-div- /QaxHost /arch:IA32 -o s2.exe ../lib/bklib.f90 ../lib/bkmath.f90 s2glob.f90 loglik.f90 ../lib/dfpmin.f90 ../lib/simann.f90 smle2.f90 s2.f90
@rem This is an experimental version
@rem ifort /MT /Qopenmp-link:static /O2 /Qipo /Qprec-div- /Qparallel /QxHost -o s2a.exe ../lib/bklib.f90 ../lib/bkmath.f90 s2glob.f90 loglik.f90 ../lib/dfpmin.f90 ../lib/simann.f90 smle2.f90 s2.f90 
ifort /O2 /Qipo /Qprec-div- /QaxHost /arch:IA32 -o s2b.exe ../lib/bklib.f90 ../lib/bkmath.f90 s2glob.f90 loglik.f90 ../lib/dfpmin.f90 ../lib/simann.f90 smle2.f90 s2.f90
 