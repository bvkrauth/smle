capture program drop smle
* SMLE command: Estimates SMLE model
program define smle, eclass
	/* If we are "replaying" (i.e. the user typed SMLE with no arguments, or ESTIMATES REPLAY) just display the most recent results */
	if (replay()){;
		di_smle `0'
		exit
	}
	syntax varlist(min = 2) [if] [in] , [ AGGregate(varlist) GRoupid(varname) NPeers(varname) PEERavg(varname) /*
		integers: */ RESTarts(integer 3) NSim(integer 100) /*
		reals: */ rho(real 0.0) GAMma(real 0.0) /*
		optionally_on: */ FIXGamma loadu UNDerreporting_correction EXEcute replace /*
		strings: */ EQuilibrium(string) RHOType(string) SIMulator_type(string) OPTimizer(string) COVmat(string)]
	** Add defaults 
	*** Optionally on 
	foreach opt in fixgamma loadu underreporting_correction {
		if "``opt''" == "" {
			local `opt' ".false." 
		}
		else {
			local `opt' ".true."
		}
	}
	*** Optional strings
	if "`simulator_type'" == "" {
		local simulator_type "GHK"
	}
	if "`optimizer'" == "" {
		local optimizer "DFP"
	}
	if "`covmat'" == "" {
		local covmat "None"
	}
	if "`equilibrium'" == "" {
		local equilibrium "Low"
	}
	if "`rhotype'" == "" {
		local rhotype "X"
	}
	** Check arguments for validity
	*** Numerical arguments 
	if (`restarts' < 0) { 
		di as error "Invalid value for restarts: `restarts' (must be >= 0)"
		error 198
	}
	if (`nsim' < 1) { 
		di as error "Invalid value for nsim: `nsim' (must be >= 1)"
		error 198
	}
	if !inrange(`rho',-1,1) { 
		di as error "Invalid value for rho: `rho' (must be between -1 and 1)" /* TODO: Find out actual limits */
		error 198
	}
	if !inrange(`gamma',.,.) { 
		di as error "Invalid value for gamma: `gamma' (must be finite)" /* TODO: Find out actual limits */
		error 198
	}
	*** String arguments 
	if !inlist(upper(substr("`equilibrium'",1,1)), "L","H","R","B","P","M") { 
		di as error "Invalid value for equilibrium: `equilibrium' (should be Low, High, Random, Bounds, Plot or Minimum)"
		error 198
	}
	if (("`execute'" == "execute") & inlist(upper(substr("`equilibrium'",1,1)), "B","P","M")) {
		di as error "Cannot execute for equilibrium = `equilibrium'"
		error 198
	}
	if !inlist(upper(substr("`rhotype'",1,1)), "X","F","E","I") { 
		di as error "Invalid value for rhotype: `rhotype' (should be X, Fixed, Estimate, or Interval)"
		error 198
	}
	if (("`execute'" == "execute") & inlist(upper(substr("`rhotype'",1,1)), "I")) {
		di as error "Cannot execute for rhotype = `rhotype'"
		error 198
	}
	if !inlist(upper(substr("`simulator_type'",1,1)), "G", "H") { 
		di as error "Invalid value for simulator_type: `simulator_type' (should be GHK or Hybrid)"
		error 198
	}
	if !inlist(upper(substr("`optimizer'",1,1)), "D", "S") { 
		di as error "Invalid value for optimizer: `optimizer' (should be DFP or SA)"
		error 198
	}
	*** Variable lists 
	if ("`groupid'" == "" & ("`npeers'" == "" | "`peeravg'" == "")) {
		di as error "Cannot determine peer choice:you must specify groupid OR npeers and peeravg
		error 198
	}
	** Find/create folder
	local folder "smle_tmp"
	capture mkdir "`folder'"
	local pwd = c(pwd)
	local folder "`pwd'\\`folder'"
	di "folder = `folder'"
	** Get file names
	local parmfile "parm.txt"
	local datafile "data.txt"
	local resultfile "results.txt"
	local logfile "log.txt"
	local ufile "u.txt"
	local bootfile "boot.txt"
	if (c(os) == "Windows") {
		quietly findfile "smle.exe"
	}
	if (c(os) == "Unix") {
		quietly findfile "smle"
	}
	local smle_exe = r(fn)
	di "smle executable in `smle_exe'"
	if (c(os) == "Windows") {
		quietly findfile "s2.exe"
	}
	if (c(os) == "Unix") {
		quietly findfile "s2"
	}
	local s2_exe = r(fn)
	di "s2 executable in `s2_exe'"
	** Confirm that files can be read/written
	if ("`replace'" != "replace") {
		confirm new file "`folder'\`parmfile'"
		confirm new file "`folder'\`datafile'"
		confirm new file "`folder'\`resultfile'"
		confirm new file "`folder'\`logfile'"
	}
	if ("`loadu'"==".false." ) {
		if ("`replace'" != "replace") {
			confirm new file "`folder'\`ufile'"
			confirm new file "`folder'\`bootfile'"
		}
	}
	else {
		confirm file "`folder'\`ufile'"
		confirm file "`folder'\`bootfile'"
	}
	** Prepare data for export
	local numagg : word count `aggregate'
	di "aggregate = `aggregate', numagg = `numagg'"
	quietly count
	local nraw = r(N)
	marksample touse
	markout `touse' `aggregate' `groupid' `npeers' `peeravg'
	quietly count if `touse'
	local nobs = r(N)
	di "`nraw' total observations in data, `nobs' will be used"
	local nvar : word count `varlist'
	local nvar = `nvar' +`numagg'- 1
	di "`nvar' explanatory variables, `numagg' are aggregate"
	gettoken own_choice xvars : varlist
	capture assert inlist(`own_choice',0,1)
	if _rc {
		di as error "Own choice variable `own_choice' must be 0/1"
		error 450
	}
	if ("`groupid'" == "") {
		capture assert inrange(`peeravg',0,1)
		if _rc {
			di as error "Peer choice variable `peeravg' must be between 0 and 1"
			error 459
		}
		capture assert `npeers' > 0
		if _rc {
			di as error "Number of peers `npeers' must be > 0"
			error 459 
		}
	}
	quietly correlate `aggregate' `xvars' if `touse', cov
	if (det(r(C)) == 0) {
		di as error "Explanatory variables are perfectly collinear, please drop some."
		error 459
	}
	** Create data file
	if "`groupid'" == "" {
		outsheet `own_choice' `peeravg' `npeers' `aggregate' `xvars' using "`folder'/`datafile'" if `touse', nonames nolabel `replace'
	}
	else {
		outsheet `groupid' `own_choice' `aggregate' `xvars' using "`folder'/`datafile'" if `touse', nonames nolabel `replace'
	}
	** Create parameter file
	if "`groupid'" == "" {
		file open fh using "`folder'/`parmfile'", write text `replace'
		file write fh "Parameter file `parmfile': Created by Stata command:" _newline
		file write fh `"smle `0' "' _newline 
		file write fh "." _newline 
		file write fh "NVAR: number of exogenous explanatory variables" _newline
		file write fh "`nvar'" _newline
		file write fh "NOBS: number of observations" _newline
		file write fh "`nobs'" _newline  
		file write fh "NUMAGG: number of exogenous explanatory variables that are aggregates" _newline
		file write fh "`numagg'" _newline
		file write fh "NSIM: number of simulations to use in calculating estimated loglikelihood" _newline
		file write fh "`nsim'" _newline
		file write fh "RESTARTS: number of times to run search algorithm" _newline
		file write fh "`restarts'" _newline
		file write fh "SIMULATOR_TYPE: simulator to use in calculating normal rectangle probabilities" _newline
		file write fh "`simulator_type'" _newline  
		file write fh "EQUILIBRIUM_TYPE: equilibrium selection rule, either low, random, or high" _newline
		file write fh "`equilibrium'" _newline
		file write fh "UNDERREPORTING_CORRECTION: correct for underreporting" _newline
		file write fh "`underreporting_correction'" _newline  
		file write fh "BOOTSTRAP: " _newline
		file write fh ".false." _newline  
		file write fh "LOAD_U: .true. if you want random numbers loaded from UFILE, .false. if you want new random numbers" _newline
		file write fh "`loadu'" _newline
		file write fh "RHO_TYPE:" _newline
		file write fh "`rhotype'" _newline
		file write fh "FIXED_RHO: value to fix rho_e at if RHO_TYPE=Fixed" _newline
		file write fh "`rho'" _newline
		file write fh "FIX_GAMMA: normally gamma is estimated, but it is fixed if this is .true." _newline
		file write fh "`fixgamma'" _newline
		file write fh "FIXED_GAMMA: value to fix gamma at if FIX_GAMMA=.true." _newline
		file write fh "`gamma'" _newline
		file write fh "DATAFILE: name of file where data is located" _newline
		file write fh "`folder'\\`datafile'" _newline
		file write fh "LOGFILE: name of file to send logging information" _newline
		file write fh "`logfile'" _newline
		file write fh "RESULTFILE: name of file to which results should be appended" _newline
		file write fh "`resultfile'" _newline
		file write fh "UFILE: name of file to which random numbers should be written (if LOAD_U=.false.) or read (if LOAD_U=.true.)" _newline
		file write fh "`ufile'" _newline
		file write fh "BOOTFILE: " _newline
		file write fh "`bootfile'" _newline 
		file write fh "FIXEDEFFECTS: " _newline
		file write fh "0" _newline 
		file close fh
		di "Data file `datafile' and parameter file `parmfile' created.  
		di "To estimate model, put these two files in the same directory as smle.exe, and execute the command" 
		di as result " smle `parmfile'"
	} 
	else {
*		file open fh using "`folder'/`parmfile'", write text `replace' /* WORKAROUND */
		file open fh using "`folder'/parm.dat", write text `replace' /* WORKAROUND */
		file write fh `"Parameter file `parmfile' created by Stata command: smle"' _newline
		file write fh "DATAFILE: name of file where data is located" _newline
*		file write fh "`folder'\\`datafile'" _newline /* WORKAROUND */
		file write fh "`datafile'" _newline /* WORKAROUND */
		file write fh "RESULTFILE: name of file to which results should be appended" _newline
		file write fh "`resultfile'" _newline
		file write fh "LOGFILE: name of file to send logging information" _newline
		file write fh "`logfile'" _newline
		file write fh "NOBS: number of observations" _newline
		file write fh "`nobs'" _newline  
		file write fh "NVAR: number of exogenous explanatory variables" _newline
		file write fh "`nvar'" _newline
		file write fh "NUMAGG: number of exogenous explanatory variables that are aggregates" _newline
		file write fh "`numagg'" _newline
		file write fh "NSIM: number of simulations to use in calculating estimated loglikelihood" _newline
		file write fh "`nsim'" _newline
		file write fh "SEARCH_METHOD: DFP or SA" _newline
		file write fh "DFP" _newline /* UPDATE THIS! */
		file write fh "RESTARTS: number of times to run search algorithm" _newline
		file write fh "`restarts'" _newline
		file write fh "EQUILIBRIUM_TYPE: equilibrium selection rule, either low, random, or high" _newline
		file write fh "`equilibrium'" _newline
		file write fh "COVMAT_TYPE: OPG, Hessian or None" _newline
		file write fh "None" _newline /* UPDATE THIS! */  
		file write fh "RHO_TYPE:" _newline
		file write fh "`rhotype'" _newline
		file write fh "FIXED_RHO: value to fix rho_e at if RHO_TYPE=Fixed" _newline
		file write fh "`rho'" _newline
		file write fh "FIX_GAMMA: normally gamma is estimated, but it is fixed if this is .true." _newline
		file write fh "`fixgamma'" _newline
		file write fh "FIXED_GAMMA: value to fix gamma at if FIX_GAMMA=.true." _newline
		file write fh "`gamma'" _newline
		file write fh "LOAD_U: .true. if you want random numbers loaded from UFILE, .false. if you want new random numbers" _newline
		file write fh "`loadu'" _newline
		file write fh "UFILE: name of file to which random numbers should be written (if LOAD_U=.false.) or read (if LOAD_U=.true.)" _newline
		file write fh "`ufile'" _newline
		file close fh
		di "Data file `datafile' and parameter file `parmfile' created.  
		di "To estimate model, put these two files in the same directory as smle.exe, and execute the command" 
		di as result " smle `parmfile'"
	}
	** Execute
	if "`execute'" == "execute" {
		capture erase "`resultfile'"
		di "Attempting to execute (this could take a long time)"
		if ("`groupid'" == "") {
			winexec `smle_exe' ".\smle_tmp\parm.txt"
		}
		else {
			cd smle_tmp
			capture erase "`resultfile'" /* WORKAROUND */
			winexec ../s2.exe 
		}
		/* The following lines of code pauses the Stata program until the RCR program has ended.  */
		sleep 1000
		/* Check to see if the output_file exists */
		capture confirm file "`resultfile'"
		/* Then repeat the following until it exists */
		while _rc != 0 {
			di "." _continue
			/* Wait five seconds */
			sleep 5000
			/* 	Check to see if the file is present */
			capture confirm file "`resultfile'"
		}
		di "DONE." _newline 
		preserve
		quietly infile loglik rhox rhoe gamma cons `aggregate' `xvars' using "`resultfile'", clear
		if ("`groupid'" == "") {
			rename gamma `peeravg'
		}
		else {
			local peeravg "gamma"
			quietly cd ..
		}
		mkmat loglik `peeravg' `aggregate' `xvars' cons rhox rhoe, matrix(b)
		restore
		scalar loglik = b[1,1]
		matrix b = b[1,2..colsof(b)]
		local cons = colsof(b) - 2
		matname b _cons , columns(`cons')
		local bnames : colnames b
		matrix V = J(colsof(b),colsof(b),0)
		matrix rownames V = `bnames'
		matrix colnames V = `bnames'
		ereturn post b V, depname(`own_choice') obs(`nobs') esample(`touse')
		ereturn local title "SMLE estimates"
		ereturn local simulator_type "`simulator_type'"
		ereturn local equilibrium_type "`equilibrium_type'"
		ereturn local rhotype "`rhotype'"
		ereturn scalar loglik = loglik
		ereturn scalar nsim = `nsim'
		ereturn scalar restarts = `restarts'
		ereturn scalar nvar = `nvar'
		ereturn scalar numagg = `numagg'
		if (upper(substr("`rhotype'",1,1)) == "F") {
			ereturn scalar fixed_rho = `rho'
		}
		if (upper("`fixgamma'") == ".TRUE.") {
			ereturn scalar fixed_gamma = `gamma'
		}
		if (upper("`underreporting_correction'") == ".TRUE.") {
			ereturn scalar report_rate = .
		}
		ereturn local cmdline "smle `0'"
		/* This should be the very last thing added to e() */
		ereturn local cmd "smle"
		di_smle
	}
end

program di_smle
		di _n as text `"`e(title)'"' _col(50) `"Number of obs  = "' as result %12.0g e(N)
		di _col(50) as text "Number of sims = " as result %12.0g e(nsim) 
		di as text "log-likelihood = " as result %9.0g e(loglik) 
		ereturn display
		di "To get standard errors, call" _newline ". bootstrap _b : " `"`e(cmdline)'"'  
end

