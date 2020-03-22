discard
clear
sysuse census, clear

gen own = (popurban > 3328253)
gen npeers = region
gen peeravg = round(runiform()*region)/region
gen zero = 0
gen one = 1
gen oneandahalf = 1.5
quietly probit own peeravg pop
estimates store probit

if (strpos("`0'","execute") == 0) {
	local execute "noexecute"
}

if (strpos("`0'","clear") > 0) {
	capture erase "parm.txt"
	capture erase "data.txt"
	capture erase "smle_tmp\parm.txt"
	capture erase "smle_tmp\data.txt"
}

/* Basic usage */
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
/* First variable (own choice) should be binary and should vary */
rcof "noisily smle peeravg pop, peeravg(peeravg) npeers(npeers)" == 450
smle zero pop , peeravg(peeravg) npeers(npeers) replace noexecute /* THIS CRASHES THE FORTRAN PROGRAM IF WE EXECUTE */
smle one pop , peeravg(peeravg) npeers(npeers) replace noexecute /* THIS CRASHES THE FORTRAN PROGRAM IF WE EXECUTE */
/* Peer average should vary over the range 0 to 1 */
rcof "noisily smle own pop , peeravg(pop) npeers(npeers)" == 459
smle own pop , peeravg(zero) npeers(npeers) replace noexecute/* THIS CRASHES THE FORTRAN PROGRAM IF WE EXECUTE */
/* Number of peers should be a strictly positive integer. It does not need to vary */
rcof "noisily smle own pop, peeravg(peeravg) npeers(zero)" == 459
smle own pop, peeravg(peeravg) npeers(oneandahalf) replace noexecute /* THIS SHOULD GET REJECTED BUT WE DON'T CHECK WHETHER IT'S AN INTEGER */
smle own pop, peeravg(peeravg) npeers(one) replace noexecute
/* All other variables should vary and be linearly independent */
rcof "noisily smle own one, peeravg(peeravg) npeers(npeers)" == 459
rcof "noisily smle own pop pop, peeravg(peeravg) npeers(npeers)" == 459  

/* Set nsim (integer, minimum 1)*/
rcof "noisily smle own pop , peeravg(peeravg) npeers(npeers) nsim(0)" == 198
rcof "noisily smle own pop , peeravg(peeravg) npeers(npeers) nsim(2.5)" == 198
smle own pop , peeravg(peeravg) npeers(npeers) replace nsim(1) `execute'
estimates store nsim_1

/* Set restarts (integer, minimum 0) */
rcof "noisily smle own pop , peeravg(peeravg) npeers(npeers) restarts(-1)" == 198
rcof "noisily smle own pop , peeravg(peeravg) npeers(npeers) restarts(2.5)" == 198
smle own pop , peeravg(peeravg) npeers(npeers) replace restarts(0) `execute'
estimates store restarts_0

/* Add aggregate variables */
/* Aggregate variables must vary and be linearly independent of individual variables */
rcof "noisily smle own pop, peeravg(peeravg) npeers(npeers) aggregate(one)" == 459  
rcof "noisily smle own pop, peeravg(peeravg) npeers(npeers) aggregate(pop)" == 459  
smle own pop , peeravg(peeravg) npeers(npeers) replace aggregate(popurban) `execute'
estimates store numagg_1

/* Set simulator type */
/* Valid options: G(HK) H(ybrid - only first letter matters*/
rcof `"noisily smle own pop , peeravg(peeravg) npeers(npeers) replace simulator(invalid)"' == 198
smle own pop , peeravg(peeravg) npeers(npeers) replace simulator(g) noexecute
smle own pop , peeravg(peeravg) npeers(npeers) replace simulator(hybrid) `execute'
estimates store hybrid

/* Set equilibrium type */
rcof `"noisily smle own pop , peeravg(peeravg) npeers(npeers) replace equilibrium(invalid)"' == 198
smle own pop , peeravg(peeravg) npeers(npeers) replace equilibrium(high) `execute' 
estimates store high
smle own pop , peeravg(peeravg) npeers(npeers) replace equilibrium(random) `execute'
estimates store random
/* These options are not yet supported with execution */
rcof `"noisily smle own pop , peeravg(peeravg) npeers(npeers) replace equilibrium(bounds) execute"' == 198
smle own pop , peeravg(peeravg) npeers(npeers) replace equilibrium(bounds) noexecute 
rcof `"noisily smle own pop , peeravg(peeravg) npeers(npeers) replace equilibrium(plot) execute"' == 198
smle own pop , peeravg(peeravg) npeers(npeers) replace equilibrium(plot) noexecute 
rcof `"noisily smle own pop , peeravg(peeravg) npeers(npeers) replace equilibrium(minimum) execute"' == 198
smle own pop , peeravg(peeravg) npeers(npeers) replace equilibrium(minimum) noexecute 

/* Set rho value/type */
/* Valid options: X, F(ixed), E(stimate), I(nterval) - only first letter matters*/
rcof `"noisily smle own pop , peeravg(peeravg) npeers(npeers) replace rhotype(notvalid)"' == 198
/* Rho is a correlation so it needs to be between -1 and 1 */
/* NOTE: when rho > 0.99, the Fortran program automatically sets it to 0.99 */
rcof "noisily smle own pop , peeravg(peeravg) npeers(npeers) replace rho(2)" == 198
smle own pop , peeravg(peeravg) npeers(npeers) replace rhotype(fixed) rho(0.2) noexecute
smle own pop , peeravg(peeravg) npeers(npeers) replace rhotype(fixed) `execute'
estimates store rho_0
smle own pop , peeravg(peeravg) npeers(npeers) replace rhotype(estimate) `execute'
estimates store rho_estimate
/* These options are not yet supported with execution */
smle own pop , peeravg(peeravg) npeers(npeers) replace rhotype(interval) noexecute /* RUNS BUT ONLY READS FIRST ROW */
rcof `"noisily smle own pop , peeravg(peeravg) npeers(npeers) replace rhotype(interval) execute"' == 198

/* Set gamma value */
/* No restrictions on gamma, but when gamma = 0, the Fortran program makes it 0.001 to avoid divide-by-zero errors */
smle own pop , peeravg(peeravg) npeers(npeers) replace fixgamma gamma(0) `execute'
est store gamma0

/* Use underreporting correction */
smle own pop , peeravg(peeravg) npeers(npeers) replace underreporting `execute'
estimates store under

est table *, b(%9.3f) se title(`"execute = "`execute'" "') equations(1)

exit

/* Using the bootstrap */
bootstrap _b , reps(5): smle own pop , peeravg(peeravg) npeers(npeers) replace `execute'
est store boot
estat bootstrap, all

est table *, b(%9.3f) se

