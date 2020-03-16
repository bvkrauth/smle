#!/usr/bin/perl
######################################################
# This perl script simply goes into the file S2.F90,
# finds lines beginning with "! S2MC"
# and removes the "! S2MC" part from the line.
# The output is written to the new file S2MC.F90
######################################################
open(S2MCF90,">s2mc.f90") || die "cannot open s2mc.f90" ;
open(INDATA,"./s2.f90") || die "cannot open s2.f90";
while (<INDATA>){
    chop;
    $tmp=substr($_,0,6);
    if ($tmp eq "! S2MC"){
	print S2MCF90 (substr($_,6), "\n");
    } else {
	print S2MCF90 ($_,"\n");
    }
}
close (INDATA);
close(S2MCF90);


