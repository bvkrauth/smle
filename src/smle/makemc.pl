#!/usr/bin/perl
######################################################
# This perl script simply goes into the file SMLE.F90,
# finds lines beginning with "! SMC"
# and removes the "! SMC" part from the line.
# The output is written to the new file SMC.F90
######################################################
open(SMCF90,">smc.f90") || die "cannot open smc.f90" ;
open(INDATA,"./smle.f90") || die "cannot open smle.f90";
while (<INDATA>){
    chop;
    $tmp=substr($_,0,5);
    if ($tmp eq "! SMC"){
	print SMCF90 (substr($_,5), "\n");
    } else {
	print SMCF90 ($_,"\n");
    }
}
close (INDATA);
close(SMCF90);



