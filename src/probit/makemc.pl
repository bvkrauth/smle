#!/usr/bin/perl
open(PROBITMCF90,">probitmc.f90") || die "cannot open s2mc.f90" ;
open(INDATA,"./probit.f90") || die "cannot open s2.f90";
while (<INDATA>){
    chop;
    $tmp=substr($_,0,10);
    if ($tmp eq "! PROBITMC"){
	print PROBITMCF90 (substr($_,10), "\n");
    } else {
	print PROBITMCF90 ($_,"\n");
    }
}
close (INDATA);
close(PROBITMCF90);


