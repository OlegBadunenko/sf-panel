*! version 1.0.0  2Jul2025 adopted from fron_p.ado to be used with margins

//  see https://friosavila.github.io/playingwithstata/main_mlols.html
// So we just need to write a program that defines this log-likelihood function. This is stright forward:
//
// . program myols_mle
//   1.         args lnf xb lnsigma
//   2.         local sigma exp(`lnsigma')
//   3.         qui:replace `lnf' = log(normalden($ML_y1,`xb',`sigma'))
//   4. end
// Notice that this program has 3 arguments.
//
// lnf: Will store the LL for each observation
// xb: Will store the linear combination for all 
// lnsigma: Will store the 
// Notice here that we are not estimating 
//  directly because it is numerically more stable to estimate its log.
//
// Now the interesting part. Create our "predict" program:
//
// . program myols_mle_p
//   1.         syntax newvarname [if] [in] , [ mean sigma emean *]
//   2.         if "`mean'`sigma'`emean'" =="" {
//   3.             ml_p `0'
//   4.         }
//   5.         marksample touse, novarlist
//   6.         if "`mean'" !=""  {
//   7.             tempvar xb
//   8.             _predict double `xb' , eq(#1)
//   9.                 gen `typlist' `varlist' = `xb' if `touse'
//  10.                 label var `varlist' "E(y|X)"
//  11.         }
//  12.         else if "`sigma'" !=""  {
//  13.             tempvar xb
//  14.             _predict double `xb' , eq(#2)
//  15.                 gen `typlist' `varlist' = exp(`xb') if `touse'
//  16.                 label var `varlist' "E(sigma|X)"
//  17.         }
//  18.         else if "`emean'"!="" {
//  19.             tempvar xb lns
//  20.                 _predict double `xb' , eq(#1)
//  21.                 _predict double `lns' , eq(#2)
//  22.                 local sigma (exp(`lns'))
//  23.                 gen `typlist' `varlist' = exp(`xb')*exp(0.5*`sigma'^2) if `t
// > ouse'
//  24.                 label var `varlist' "E(exp(Y)|X)"
//  25.         }
//  26. end
// This program, here named -myols_mle_p-, can be used to estimate 3 things:
//
// mean: This is the standard outcome. Just looking into the linear combination of X and betas
// sigma: When this option is used, you will obtain 
//  .
// emean: And this will be something new. -emean- could be used if, say, your outcome of interest was "wages", but you were modeling "log(wages)"

capture program drop xtsf_p
program define xtsf_p
	version 8 

	local myopts "persistent transient xb"
	display "myopts"
	display "|`myopts'|"

		/* Step 2:
			call _propts, exit if done, 
			else collect what was returned.
		*/

	_pred_se "`myopts'" `0'
	if `s(done)' {
		exit
	}
	local vtyp  `s(typ)'
	local varn `s(varn)'
	local 0 `"`s(rest)'"'

		/* Step 3:
			Parse your syntax.
		*/

	syntax [if] [in] [, `myopts']

	marksample touse
	
// 	display "varn"
// 	display "|`varn'|"

	local sumt = ("`persistent'"!="") + ("`transient'"!="") + ("`xb'"!="")
	if `sumt' > 1 {
		di as err "only one statistic may be specified"
		exit 198
	}
	else { 
		if `sumt' == 0 {
			local stat "xb"
			di as txt "(option {bf:xb} assumed; fitted values)"
		}
		else 	local stat "`xb'`persistent'`transient'"
	}

	display "stat"
	display "|`stat'|"

	xtsf_p_myCalc `"`vtyp'"' `varn' `stat' `touse' 
end

capture program drop xtsf_p_myCalc
program define xtsf_p_myCalc
	args vtyp varn stat cond
			/* vtyp: type of new variable
			   varn: name of new variable
	 		   cond: if & in
		 	   stat: option specified
			*/

// 	local y `e(depvar)'
// 	display "xtsf_p_myCalc: stat"
// 	display "|`stat'|"	
	tempname mysample
	generate `mysample' = e(sample)
	quietly summarize `mysample' if `mysample' == 0
	if r(N) > 0 {
		display "{text} (" as result r(N) "{text} missing values generated)"
	}	
	
	quietly generate `varn' = .
	
	if "`stat'"=="xb" {
		display "stat in xtsf_p_myCalc"
		display "|`stat'|"
		tempvar xb
		quietly _predict double `xb' if `cond', xb
		drop `varn'
		generate `vtyp' `varn'=`xb' if `mysample'
		label variable `varn' `"`:var label `xb''"'
	}

	if "`stat'"=="persistent" {
		display "stat in xtsf_p_myCalc"
		display "|`stat'|"
		mata: st_store(., st_local("varn"), st_local("mysample"), st_matrix("e(eff_persistent)"))
		label variable `varn' "Persistent Efficiency"
	}
	if "`stat'"=="transient" {
		display "stat in xtsf_p_myCalc"
		display "|`stat'|"
		mata: st_store(., st_local("varn"), st_local("mysample"), st_matrix("e(eff_transient)"))
		label variable `varn' "Transient Efficiency"
	}

end
