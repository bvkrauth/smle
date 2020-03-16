@rem This file is used for compiling the SMLE program on an MS Windows system using the F compiler
@rem F is a subset of Fortran 90.  A free F compiler is available from the Fortran company
@rem at http://www.fortran.com
@echo Compiling...
@rem F  -O3 -o smle.exe -tempdir . ../lib/bklib.f90 ../lib/bkmath.f90 smglob.f90 ../lib/dfpmin.f90 ../lib/simann.f90 loglik.f90 smutil.f90 smle.f90
@rem del *.mod
@rem \Perl\bin\perl.exe makemcrep.pl
@rem F -O3 -o smlerep.exe -tempdir . ../lib/bklib.f90 ../lib/bkmath.f90 smglob.f90 ../lib/dfpmin.f90 ../lib/simann.f90 loglik.f90 smutil.f90 smlerep.f90
@rem del *.mod smlerep.f90
@rem F -O3 -o smc.exe -tempdir . ../lib/bklib.f90 ../lib/bkmath.f90 smglob.f90 ../lib/dfpmin.f90 ../lib/simann.f90 loglik.f90 ../psim/monte.f90 smutil.f90 smc.f90
@rem del *.mod smc.f90 


@rem Baseline version - 7:52
@rem ifort -o smle0.exe ../lib/bklib.f90 ../lib/bkmath.f90 smglob.f90 ../lib/dfpmin.f90 ../lib/simann.f90 loglik.f90 smutil.f90 smle.f90
@rem Optimized version - 7:28
@rem ifort /O2 /Qipo /Qprec-div- -o smle1.exe ../lib/bklib.f90 ../lib/bkmath.f90 smglob.f90 ../lib/dfpmin.f90 ../lib/simann.f90 loglik.f90 smutil.f90 smle.f90
@rem Uses SSE3 instruction set (Intel Core processors only) - 6:06
@rem ifort /O2 /QxSSE3 /Qipo /Qprec-div- -o smle2.exe ../lib/bklib.f90 ../lib/bkmath.f90 smglob.f90 ../lib/dfpmin.f90 ../lib/simann.f90 loglik.f90 smutil.f90 smle.f90
@rem Uses SSE3 instruction set (Intel Core processor) if available, will default to Pentium if not - 6:05
@rem ifort /O2 /QaxSSE3 /arch:IA32 /Qipo /Qprec-div- -o smle3.exe ../lib/bklib.f90 ../lib/bkmath.f90 smglob.f90 ../lib/dfpmin.f90 ../lib/simann.f90 loglik.f90 smutil.f90 smle.f90
@rem Optimized version with O3 - 7:54
@rem ifort /O3 /Qipo /Qprec-div- -o smle4.exe ../lib/bklib.f90 ../lib/bkmath.f90 smglob.f90 ../lib/dfpmin.f90 ../lib/simann.f90 loglik.f90 smutil.f90 smle.f90
@rem Parallelized- 6:42
@rem ifort /O2 /Qparallel /Qipo /Qprec-div- -o smle5.exe ../lib/bklib.f90 ../lib/bkmath.f90 smglob.f90 ../lib/dfpmin.f90 ../lib/simann.f90 loglik.f90 smutil.f90 smle.f90
@rem Parallelized with SSE3 - 7:05
@rem ifort /O2 /QxSSE3 /Qparallel /Qipo /Qprec-div- -o smle6.exe ../lib/bklib.f90 ../lib/bkmath.f90 smglob.f90 ../lib/dfpmin.f90 ../lib/simann.f90 loglik.f90 smutil.f90 smle.f90
@rem original 10:43

@rem this is what seems to work best.  It gives a 43% reduction in calculation time.
ifort /O2 /QaxSSE3 /arch:IA32 /Qipo /Qprec-div- -o smle.exe ../lib/bklib.f90 ../lib/bkmath.f90 smglob.f90 ../lib/dfpmin.f90 ../lib/simann.f90 loglik.f90 smutil.f90 smle.f90


del *.mod *.obj
