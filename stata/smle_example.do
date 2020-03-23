/* Stata do-file for demonstrating and testing the smle package */
discard
clear
sysuse census, clear

/* Generate some fake data */
gen own = (popurban > 3328253)
gen npeers = region
set seed 339487731
gen peeravg = round(runiform()*region)/region
gen zero = 0
gen one = 1
gen oneandahalf = 1.5
quietly probit own peeravg pop
estimates store probit

/* If this do-file is run with no arguments, Stata will skip the estimation step for most examples to save time */
/* If the "execute" argument is included, Stata will actually estimate all of the models */
if (strpos("`0'","execute") == 0) {
	local execute "noexecute"
}

/* If the "clear" argument is included, Stata will delete some local intermediate files to ensure clean execution */
if (strpos("`0'","clear") > 0) {
	capture erase "parm.txt"
	capture erase "data.txt"
	capture erase "smle_tmp\parm.txt"
	capture erase "smle_tmp\data.txt"
}

/* Confirm that SMLE calculation is stable */
/* This requires that we have the data file .\testing\indu.txt */
capture confirm file ".\testing\indu.txt"
if {_rc == 0) {
	smle own pop , peeravg(peeravg) npeers(npeers) restarts(1) ufile(".\testing\indu.txt") 
	estimates store indstable
	qui {
	mat T_b = J(1,5,0)
	mat T_b[1,1] =   .108166016638279
	mat T_b[1,2] =  2.80907215711e-07
	mat T_b[1,3] = -1.633952379226685
	mat T_b[1,4] =  .0111995115876198
	mat T_b[1,5] =  .0111995115876198
	}
	matrix C_b = e(b)
	assert mreldif( C_b , T_b ) < 1E-8
}

/* Confirm that S2 calculation is stable */
/* This requires that we have the data file .\testing\grpu.txt */
capture confirm file ".\testing\grpu.txt"
if {_rc == 0) {
	smle own pop , groupid(region) restarts(1) ufile(".\testing\grpu.txt")
	estimates store grstable
	qui {
	mat T_b = J(1,5,0)
	mat T_b[1,1] =  .3889549374580383
	mat T_b[1,2] =  2.98495479001e-07
	mat T_b[1,3] =  -2.06155252456665
	mat T_b[1,4] = -.0618750564754009
	mat T_b[1,5] = -.0618750564754009
	}
	matrix C_b = e(b)
	assert mreldif( C_b , T_b ) < 1E-8
}

/* Basic usage, individual-based sample */
smle own pop , peeravg(peeravg) npeers(npeers) replace `execute'
estimates store basic

/* Basic usage, group-based sample */
smle own pop, groupid(region) replace `execute'
estimates store group

/* Basic usage, no execution */
smle own pop , peeravg(peeravg) npeers(npeers) replace noexecute
confirm file "parm.txt"
confirm file "data.txt"

/* Specify location for files */
smle own pop , peeravg(peeravg) npeers(npeers) replace noexecute save("smle_tmp")
confirm file "smle_tmp\parm.txt"
confirm file "smle_tmp\data.txt"

/* Basic usage, error messages */
/* Require "replace" if the parameter and/or data file already exists */
rcof `"noisily smle own pop, peeravg(peeravg) npeers(npeers) noexecute"' == 602
/* There should be at least 2 variables */
rcof "noisily smle own , peeravg(peeravg) npeers(npeers)" == 102
/* Must specify either groupid OR peeravg and npeers */
rcof "noisily smle own pop " == 198
rcof "noisily smle own pop, peeravg(peeravg)" == 198
rcof "noisily smle own pop, npeers(npeers)" == 198
rcof "noisily smle own pop, groupid(region) peeravg(peeravg) npeers(npeers)" == 198
/* First variable (own choice) should be binary and should vary */
rcof "noisily smle peeravg pop, peeravg(peeravg) npeers(npeers)" == 450
smle zero pop , peeravg(peeravg) npeers(npeers) replace noexecute /* THIS CRASHES THE FORTRAN PROGRAM IF WE EXECUTE */
smle one pop , peeravg(peeravg) npeers(npeers) replace noexecute /* THIS CRASHES THE FORTRAN PROGRAM IF WE EXECUTE */
/* Peer average (peeravg) should vary over the range 0 to 1 */
rcof "noisily smle own pop , peeravg(pop) npeers(npeers)" == 459
smle own pop , peeravg(zero) npeers(npeers) replace noexecute/* THIS CRASHES THE FORTRAN PROGRAM IF WE EXECUTE */
/* Number of peers (npeers) should be a strictly positive integer. It does not need to vary */
rcof "noisily smle own pop, peeravg(peeravg) npeers(zero)" == 459
smle own pop, peeravg(peeravg) npeers(oneandahalf) replace noexecute /* THIS SHOULD GET REJECTED BUT WE DON'T CHECK WHETHER IT'S AN INTEGER */
smle own pop, peeravg(peeravg) npeers(one) replace noexecute
/* Group ID (groupid) should be an integer or long. */
smle own pop, groupid(peeravg) replace noexecute /* THIS CRASHES THE FORTRAN PROGRAM IF WE EXECUTE */
/* All other variables should vary and be linearly independent */
rcof "noisily smle own one, peeravg(peeravg) npeers(npeers)" == 459
rcof "noisily smle own pop pop, peeravg(peeravg) npeers(npeers)" == 459  

/* Add aggregate variables */
/* Aggregate variables must vary and be linearly independent of individual variables */
rcof "noisily smle own pop, peeravg(peeravg) npeers(npeers) aggregate(one)" == 459  
rcof "noisily smle own pop, peeravg(peeravg) npeers(npeers) aggregate(pop)" == 459  
smle own pop , peeravg(peeravg) npeers(npeers) replace aggregate(popurban) `execute'
estimates store numagg_1
smle own pop , groupid(region) replace aggregate(popurban) `execute'
estimates store gnumagg_1

/* Set rhotype and rho */
/* Valid options: x, fixed, estimate, interval */
rcof `"noisily smle own pop , peeravg(peeravg) npeers(npeers) replace rhotype(invalid)"' == 198
/* rhotype(x) is the default */
smle own pop , peeravg(peeravg) npeers(npeers) replace rhotype(x) noexecute
smle own pop , groupid(region) replace rhotype(x) noexecute
/* rhotype(fixed) uses the value given for rho */
/* Rho is a correlation so it needs to be between -1 and 1 */
/* NOTE: when rho > 0.99, the Fortran program automatically sets it to 0.99 */
rcof "noisily smle own pop , peeravg(peeravg) npeers(npeers) replace rho(2)" == 198
rcof "noisily smle own pop , peeravg(peeravg) npeers(npeers) replace rho(-2)" == 198
smle own pop , peeravg(peeravg) npeers(npeers) replace rhotype(fixed) rho(0.2) noexecute
smle own pop , peeravg(peeravg) npeers(npeers) replace rhotype(fixed) `execute'
estimates store rho_0
smle own pop , groupid(region) replace rhotype(fixed) `execute'
estimates store grho_0
/* rhotype(estimate) estimates rho directly */
smle own pop , peeravg(peeravg) npeers(npeers) replace rhotype(estimate) `execute'
estimates store rho_estimate
smle own pop , groupid(region) replace rhotype(estimate) `execute'
estimates store grho_estimate
/* rhotype(interval) is not yet supported with execution */
rcof `"noisily smle own pop , peeravg(peeravg) npeers(npeers) replace rhotype(interval) execute"' == 198
rcof `"noisily smle own pop , groupid(region) replace rhotype(interval) execute"' == 198
smle own pop , peeravg(peeravg) npeers(npeers) replace rhotype(interval) noexecute /* RUNS BUT ONLY READS FIRST ROW */
smle own pop , groupid(region) replace rhotype(interval) noexecute /* RUNS BUT ONLY READS FIRST ROW */

/* Set fixgamma and gamma */
/* No restrictions on gamma, but when gamma = 0, the Fortran program makes it 0.001 to avoid divide-by-zero errors */
smle own pop , peeravg(peeravg) npeers(npeers) replace fixgamma gamma(0) `execute'
est store gamma0
smle own pop , groupid(region) replace fixgamma gamma(0) `execute'
est store ggamma0

/* Set equilibrium type */
/* Valid options: low (default, high, random, bounds, minimum, plot */
rcof `"noisily smle own pop , peeravg(peeravg) npeers(npeers) replace equilibrium(invalid)"' == 198
/* equilibrium(low) is the default */
smle own pop , peeravg(peeravg) npeers(npeers) replace equilibrium(low) noexecute 
smle own pop , groupid(region) replace equilibrium(low) noexecute 
/* equilibrium(high) assumes high equilibrium is selected */
smle own pop , peeravg(peeravg) npeers(npeers) replace equilibrium(high) `execute' 
estimates store high
smle own pop , groupid(region) replace equilibrium(high) `execute' 
estimates store ghigh
/* equilibrium(random) randomly selects an equilibrium */
smle own pop , peeravg(peeravg) npeers(npeers) replace equilibrium(random) `execute'
estimates store random
smle own pop , groupid(region) replace equilibrium(random) `execute'
estimates store grandom
/* Options bounds, minimum and plot are not yet supported with execution */
rcof `"noisily smle own pop , peeravg(peeravg) npeers(npeers) replace equilibrium(bounds) execute"' == 198
rcof `"noisily smle own pop , groupid(region) replace equilibrium(bounds) execute"' == 198
smle own pop , peeravg(peeravg) npeers(npeers) replace equilibrium(bounds) noexecute 
smle own pop , groupid(region) replace equilibrium(bounds) noexecute 
rcof `"noisily smle own pop , peeravg(peeravg) npeers(npeers) replace equilibrium(minimum) execute"' == 198
rcof `"noisily smle own pop , groupid(region) replace equilibrium(minimum) execute"' == 198
smle own pop , peeravg(peeravg) npeers(npeers) replace equilibrium(minimum) noexecute 
smle own pop , groupid(region) replace equilibrium(minimum) noexecute 
rcof `"noisily smle own pop , peeravg(peeravg) npeers(npeers) replace equilibrium(plot) execute"' == 198
rcof `"noisily smle own pop , groupid(region) replace equilibrium(plot) execute"' == 198
smle own pop , peeravg(peeravg) npeers(npeers) replace equilibrium(plot) noexecute 
smle own pop , groupid(region) replace equilibrium(plot) noexecute 

/* Use underreporting correction */
smle own pop , peeravg(peeravg) npeers(npeers) replace underreporting `execute'
estimates store under
/* Not available for group-based sample */
rcof "noisily smle own pop , groupid(region) replace underreporting noexecute" == 198

/* Set simulator */
/* Valid options: ghk, hybrid - I'm pretty sure GHK only works with Low + individual-based sample */
rcof `"noisily smle own pop , peeravg(peeravg) npeers(npeers) replace simulator(invalid)"' == 198
/* ghk is the default */
smle own pop , peeravg(peeravg) npeers(npeers) replace simulator(ghk) noexecute
smle own pop , groupid(region) replace simulator(ghk) noexecute
/* hybrid */
smle own pop , peeravg(peeravg) npeers(npeers) replace simulator(hybrid) `execute'
estimates store hybrid
smle own pop , groupid(region) replace simulator(hybrid) `execute'
estimates store ghybrid

/* Set nsim (integer, minimum 1)*/
rcof "noisily smle own pop , peeravg(peeravg) npeers(npeers) nsim(0)" == 198
rcof "noisily smle own pop , peeravg(peeravg) npeers(npeers) nsim(2.5)" == 198
smle own pop , peeravg(peeravg) npeers(npeers) replace nsim(1) `execute'
estimates store insim_1
smle own pop , groupid(region) replace nsim(1) `execute'
estimates store gnsim_1

/* Set optimizer */
/* Valid options: dfp, sa */
rcof `"noisily smle own pop , peeravg(peeravg) npeers(npeers) replace optimizer(invalid)"' == 198
/* dfp is the default */
smle own pop , peeravg(peeravg) npeers(npeers) replace optimizer(dfp) noexecute
smle own pop , groupid(region) replace optimizer(dfp) noexecute
/* sa is simulated annealing */
smle own pop , peeravg(peeravg) npeers(npeers) replace optimizer(sa) `execute'
estimates store sa
smle own pop , groupid(region) replace optimizer(sa) `execute'
estimates store gsa

/* Set restarts (integer, minimum 0) */
rcof "noisily smle own pop , peeravg(peeravg) npeers(npeers) restarts(-1)" == 198
rcof "noisily smle own pop , peeravg(peeravg) npeers(npeers) restarts(2.5)" == 198
/* We already used restarts(1) above, so no need to test it */

est table *, b(%9.3f) se title(`"execute = "`execute'" "') equations(1)

exit

/* Using the bootstrap */
bootstrap _b , reps(5): smle own pop , peeravg(peeravg) npeers(npeers) replace `execute'
est store boot
estat bootstrap, all

est table *, b(%9.3f) se

