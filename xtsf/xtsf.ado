*! version 1.0.0  18Jun2025
*! version 1.1.0  16Jul2025
*! version 1.2.0  23Feb2025
*! author Oleg Badunenko


// if(c(MP)){
//  	set processors 1
// }

capture program drop xtsf
program define xtsf, eclass
	version 14

  if !replay() {

		syntax varlist(numeric fv min=2) [if] [in]                  ///
			[, noCONStant ROBUST COST MODEL(string)  ///
			UITLNVariance(string) VITLNVariance(string)  ///
			UILNVariance(string) VILNVariance(string)  ///
			FROM(string) INITSCAF(real 1) INITSCAU(real 0.25) INITSCAV(real 0.25) ///
			simtype(string) R(integer 199)  ///
			HALTONBase(integer 2)  BASESApart(integer 1)  ///
			ITERate(integer 102) TRACElevel(string) ///
			LEVEL(string) NOLOG ] [NTHREADS(integer 1)] ///
			[TECHNIQUE(string)] ///
			[noCI] [noPValues] [noOMITted] [noEMPTYcells] ///
			[VSQUISH] [BASElevels] [ALLBASElevels] [noFVLABel] ///
			[fvwrap(passthru)] [fvwrapon(passthru)] ///
			[CFORMAT(passthru)] [PFORMAT(passthru)] ///
			[SFORMAT(passthru)] [nolstretch] 

		tempname b V n1 nt1 R2 R2adj aic bic llopt  ///
			rez fitted shat RSS mypanelvar mytimevar convstatus ///
			u_mean u_mode eff_mean eff_mode eff_bc eff_lb eff_ub ///
			mylevel iter dist uizeromean modelcoefs modelfn ///
			myprod function cnames vitlnvarianceN uilnvarianceN uimeanN ///
			eff_pers eff_tran theta_re


		// handle specifications of uimean uilnvariance vitlnvariance
		foreach opt in uitlnvariance vitlnvariance uilnvariance vilnvariance {
			if "``opt''" != "" {
				local `opt'opt ``opt''
				tokenize "``opt'opt'", parse(,)
				fvunab `opt' : `1'
				if "`s(fvops)'" == "" {
					confirm numeric var ``opt''
				}
				/* noCONStant */
				if "`3'" != "" {
					local l = length(`"`3'"')
					if `"`3'"' == bsubstr("noconstant", /*
					*/ 1,max(6,`l')) {
						local `opt'nocns "noconstant"
					}
					else {
						di as err "`3' invalid"
						exit 198
					}
				}
			}
		}

		// handle from
		if "`from'" == "" { 
			local init_supplied = 0
		}
		else {
			local init_supplied = 1
			matrix _from_192837465 = `from'
		}

		// handle simtype
		if "`simtype'" == "" | "`simtype'" == "halton"{ 
			local mySimtype = "1"
		}
		else {
			local mySimtype = "2"
		}

		// 	display "from"
		// 	display "|`from'|"

		// handle level
		if ("`level'" == "") {
			local mylevel 95
		}
		else 	if `level' < 10 |  `level' > 99.99 {
			display "{p 1 1 7}{error}level() must be between 10 and 99.99 inclusive{p_end}"
			exit 198
		}
		else {
			local mylevel `level'
		}

		// handle iterate
		if `iterate' < 0  {
			display as error " iterate() must be greater than 0"
			exit 198
		}

		// handle tracelevel
		if "`nolog'" == "" {
			if "`tracelevel'" == ""{
				local tracelevel "value"
			}
			if "`tracelevel'" != "none" & "`tracelevel'" != "value"  & ///
				"`tracelevel'" != "tolerance" & "`tracelevel'" != "step" & ///
				"`tracelevel'" != "paramdiffs" & "`tracelevel'" != "params" & ///
				"`tracelevel'" != "gradient" & "`tracelevel'" != "hessian" {
				di as error "tracelevel() can be specified as follows"
				di as input "none, value, tolerance, step, paramdiffs, params, gradient, or hessian"
				di as text "see help for " as result "mf_optimize##i_tracelevel" as text " for details"
				exit 198
			}
		}
		else {
			local tracelevel "none"
		}

		// handle technique
		if "`technique'" == ""{
			local technique "nr" 
		} 
		else if "`technique'" == "bhhh" | "`technique'" == "bfgs" | "`technique'" == "nr" | "`technique'" == "dfp"{
      // good
		} 
		else {
			display as text "see help for " as result "'help mf_optimize##i_technique'" as text " for details"
			// 		display as error "technique() is not appropriate"
			// 		exi 198
		}

		// handle model
		local model = lower("`model'")
		// 	display "`model'"
		if "`model'" == "" {
			local model "gtremsle"
		}
		if "`model'" != "gtremsle" & "`model'" != "gtrefull" & "`model'" != "tre" {
			di as error "model() is not appropriate"
			exit 198
		}
		else if "`model'" == "gtremsle" | "`model'" == "" {
			local model 41
			scalar _vu0zero_192837465 = 0
			scalar _vv0zero_192837465 = 0
			scalar _vuizero_192837465 = 0
			local modelfn "4comp MSLE"
		}
		else if "`model'" == "gtrefull" {
			local model 42
			scalar _vu0zero_192837465 = 0
			scalar _vv0zero_192837465 = 0
			scalar _vuizero_192837465 = 0
			local modelfn "4comp Full"
		}
		else if "`model'" == "tre" {
			local model 3
			scalar _vu0zero_192837465 = 1
			scalar _vv0zero_192837465 = 0
			scalar _vuizero_192837465 = 0
			//     local modelcoefs "beta[t]:gamma beta[t]:delta"
			//     local modelfn "beta[t] = 1 + gamma*(t-T_i) + delta*(t-T_i)^2"
		}

		// handle R
		if "`r'" == "" { 
			scalar _r_192837465 = 250
		}
		else {
			scalar _r_192837465 = `r'		
		}
		// handle base
		if "`haltonbase'" == "" { 
			scalar _haltonb_192837465 = 1
		}
		else {
			scalar _haltonb_192837465 = `haltonbase'
		}
		// base apart
		if "`basesapart'" == "" { 
			scalar _basesa_192837465 = 1
		}
		else {
			scalar _basesa_192837465 = `basesapart'
		}

		// 	  di 5
		// Multiprocessing
		//   scalar `_nthreads' = "`nthreads'"
		scalar _nthreads_192837465 = `nthreads'
		// 	scalar list _nthreads_192837465

		// handle production/cost function
		if "`cost'" == "" { 
			local myprod = 1 
			local function "production"
			scalar _prod_192837465 = 1
		}
		else {
			local myprod = -1
			local function "cost"
			scalar _prod_192837465 = -1
		}

		marksample touse
		// handle the lists
		gettoken depvar indepvars : varlist
		_fv_check_depvar `depvar'
		_rmcoll `indepvars' if `touse', expand `constant'
		local indepvars `r(varlist)'

		markout `touse' `vitlnvariance' `uilnvariance' `uimean'

		quietly count if `touse'==1
		scalar _n_192837465 = r(N)
		if r(N) == 0 {
			error 2000
		}

		quietly xtset
		local mypanelvar `r(panelvar)'
		local mytimevar `r(timevar)'

		quietly xtreg `depvar' `indepvars' if `touse'
		// 	matrix `theta_re'  = e(b)
		// 	matrix list `theta_re'
		// 	matrix _theta_re_192837465 = `theta_re'
		matrix _theta_re_192837465 = e(b)


		//   regressors
		if "`constant'" == "" {
			local cnames "`indepvars' _cons"
		}
		else {
			local cnames "`indepvars'"
		}

		// handle error components

		foreach errorcomp in ui uit vi vit {
			if "`errcomplnvariance'" != "" {
        local errcomplnvarianceN ""
        foreach lname of local errcomplnvariance {
					local errcomplnvarianceN "`errcomplnvarianceN' ln[var(`errorcomp')]:`lname'"
        }
        if "`errcomplnvariancenocns'" == "" {
					local cnames_`errorcomp' "`cnames_`errorcomp'' `errcomplnvarianceN' ln[var(`errorcomp')]:_cons"
        }
        else {
					local cnames_`errorcomp' "`cnames_`errorcomp'' `errcomplnvarianceN'"
        }
			}
			else {
        local cnames_`errorcomp' "`cnames_`errorcomp'' ln[var(`errorcomp')]:_cons"
			}
		}

		local cnames "`cnames' `cnames_vi' `cnames_ui' `cnames_vit' `cnames_uit'"
		// 	display "`cnames'"

		//   matrix _lnls_192837465 = J(_n_192837465,1,17)

		mata: xtsf_vykonaty("`depvar'", /// 1
			"`indepvars'", /// 2
			"`touse'", /// 3
			"`constant'", /// 4
			"`mypanelvar'", /// 5
			"`mytimevar'", /// 6
			"`model'", /// 7
			"`robust'", /// 8
			"`uitlnvariance'", /// 9
			"`uitlnvariancenocns'", /// 10
			"`vitlnvariance'", /// 11
			"`vitlnvariancenocns'", /// 12
			"`uilnvariance'", /// 13
			"`uilnvariancenocns'", /// 14
			"`vilnvariance'", /// 15
			"`vilnvariancenocns'", /// 16
			"`b'", /// 17
			"`V'", /// 18
			"`n1'", /// 19
			"`nt1'", /// 20
			"`eff_pers'", /// 21
			"`eff_tran'", /// 22
			"`technique'", /// 23
			"`iterate'", /// 24
			"`initscaf'", /// 25
			"`initscau'", /// 26
			"`initscav'", /// 27
			"`tracelevel'", /// 28
			"`convstatus'", /// 29
			"`mylevel'", /// 30
			"`rez'", /// 31
			"`llopt'" /// 32
			)

		// 	mata: mata drop tymV
		// 	mata: mata drop tymU

		// di 11

		matrix colnames `b'        = `cnames'
		// di 12
		// mat list `b'
		// di 13
		matrix rownames `V'        = `cnames'
		// 	di 14
		matrix colnames `V'        = `cnames'
		// 	di 15
		ereturn post `b' `V', esample(`touse') buildfvinfo depname("`depvar'")
		//   matrix colnames `eff_jlms'   = "Efficiency"
		// di 16
		ereturn matrix eff_persistent    = `eff_pers'
		// 	di 17
		ereturn matrix eff_transient    = `eff_tran'
		// 	di 18
		ereturn scalar N           = `n1'
		// 	di 19
		ereturn scalar sumTi       = `nt1'
		// 	di 20
		ereturn scalar converged   = `convstatus'
		// 	di 21
		ereturn scalar level       = `mylevel'
		// 	di 22
		ereturn scalar ll          = `llopt'

		//   ereturn matrix residuals   = `rez'
		//   ereturn matrix xb          = `fitted'

		ereturn local predict "xtsf_p"
		ereturn local cmd   "xtsf"
		ereturn local cmdline "`0'"

  }
	if replay() {
		//     display "replay here"
    syntax, [LEVel(real `c(level)')] [noCI] [noPValues] [noOMITted] [noEMPTYcells] [VSQUISH] [BASElevels] [ALLBASElevels] [noFVLABel] [fvwrap(passthru)] [fvwrapon(passthru)] ///
			[CFORMAT(passthru)] [PFORMAT(passthru)] [SFORMAT(passthru)] [nolstretch]
  }
  if "`nolog'" == "" {
    display
    display as result "Sample:" as input "{hline 22}
    display as input " Number of obs    " as text "= " as result `nt1'
    display as input " Number of groups " as text "= " as result `n1'
    display as result "Diagnostics:" as input "{hline 17}
    display as input " Simulated logL   " as text "= " as result  %-15.4f `llopt'
    display as input "{hline 29}"

    display as input _newline " The 4-components estimator of"
    display as text " the " as input "`function' " as text "stochastic frontier model for panel data"
		//     display as text " where effciency is" as input " time-varying"
    display

		ereturn display, level(`mylevel')	`ci' `pvalues' `omitted' `emptycells' `vsquish' `baselevels' `allbaselevels' `fvlabel' `fvwrap' `fvwrapon' `cformat' `pformat' `sformat' `lstretch'
		//     display as result "Note, the function of inefficiency change over time is"
		//     display as input "`modelfn'" as result ": "
		//     display as text " the larger is the beta[t] the larger is the inefficiency"
	}

end

**# xtsf_vykonaty


capture mata mata drop xtsf_vykonaty()

mata:

	void xtsf_vykonaty( string scalar depvar, /// 1
		string scalar indepvars, /// 2
		string scalar touse, /// 3
		string scalar constant, /// 4
		string scalar mypanelvar, /// 5
		string scalar mytimevar, /// 6
		string scalar model0, /// 7
		string scalar robust, /// 8
		string scalar uitlnvariance, /// 9
		string scalar uitlnvariancenocns, /// 10
		string scalar vitlnvariance, /// 11
		string scalar vitlnvariancenocns, /// 12
		string scalar uilnvariance, /// 13
		string scalar uilnvariancenocns, /// 14
		string scalar vilnvariance, /// 15
		string scalar vilnvariancenocns, /// 16
		string scalar bsname, /// 17
		string scalar vsname, /// 18
		string scalar n1name, /// 19
		string scalar nt1name, /// 20
		string scalar eff_pers_name, /// 21
		string scalar eff_tran_name, /// 22
		string scalar technique, /// 23
		string scalar iter0, /// 24
		string scalar initscaf0, /// 25
		string scalar initscau0, /// 26
		string scalar initscav0, /// 27
		string scalar tracelevel, /// 28
		string scalar convstatus, /// 29
		string scalar level0, /// 30
		string scalar rezname, /// 31
		string scalar LLname /// 32
		)
	{

		real scalar init_supplied
		real scalar mySimtype


		model     = strtoreal(model0)
		// 	"MODEL"
		// 	model
		level     = strtoreal(level0)
		alpha1    = (100-level)/100
		//   100


		ids0      = st_data(., mypanelvar, touse)
		ids       = panelsetup(ids0, 1)

		unique_ids = uniqrows(ids0) // all unique values in mypanelvar (not necessarily sorted or in panel order)
		ids       = ids, ids[,2] - ids[,1] :+ 1
		nobs      = rows(ids)

		maxTi = max(ids[,3]) // Max number of obs per ID
		// 	maxTi

		//   101

		Ti0   = st_data(., mytimevar, touse)
		Ti    = Ti0 :- min(Ti0) :+ 1
		//   maxTi = max(Ti)


		nt        = rows(Ti)
		tLESSmaxT = J(nt,1,.)
		for (i = 1; i <= nobs; i++) {
			tymch1  = panelsubmatrix(Ti,  i, ids)
			tLESSmaxT[ids[i,1]..ids[i,2]] = tymch1 :- max( tymch1 )
		}
		// 102

		// ## Data

		yit       = st_data(., depvar, touse)

		xit       = st_data(., indepvars, touse)
		if (constant == "") {
			xit = xit, J(nt, 1, 1)
		}
		k         = cols(xit)

		//   103
		// --- vit component ---
		if (vitlnvariance == "") {
			zvit    = J(nt, 1, 1)
			kv      = 1
		}
		else {
			zvit    = st_data(., vitlnvariance, touse)
			if (vitlnvariancenocns == "") {
				zvit  = zvit, J(nt, 1, 1)
			}
			kv      = cols(zvit)
		}

		// --- uit component ---
		if (uitlnvariance == "") {
			zuit = J(nt, 1, 1)
			ku   = 1
		} else {
			zuit = st_data(., uitlnvariance, touse)
			if (uitlnvariancenocns == "") {
        zuit = zuit, J(nt, 1, 1)
			}
			ku = cols(zuit)
		}



		//   104


		// --- ui component ---
		if (uilnvariance == "") {
			zui  = J(nt, 1, 1)
			ku0  = 1
		} else {
			zui = st_data(., uilnvariance, touse)
			if (uilnvariancenocns == "") {
        zui = zui, J(nt, 1, 1)
			}
			ku0 = cols(zui)
		}

		// --- vi component ---
		if (vilnvariance == "") {
			zvi  = J(nt, 1, 1)
			kv0  = 1
		} else {
			zvi = st_data(., vilnvariance, touse)
			if (vilnvariancenocns == "") {
        zvi = zvi, J(nt, 1, 1)
			}
			kv0 = cols(zvi)
		}

		//   zui

		//   105
		//   if (cost == ""){
		//     s = -1
		//   }
		//   else {
		//     s = 1
		//   }
		//	
		// 	1051
		// 	s
		//
		// 	s = (cost == "") ? -1 : 1
		//
		// 	1052
		// 	s

		//   107

		// ## starting values

		initscaf = strtoreal(initscaf0)
		initscau = strtoreal(initscau0)
		initscav = strtoreal(initscav0)



		init_supplied = strtoreal(st_local("init_supplied"))
		mySimtype = strtoreal(st_local("mySimtype"))
		// 	init_supplied

		if (init_supplied == 1) {
			// 		printf("init_supplied is 1\n")
			stata("mat list _from_192837465")
			theta0 = st_matrix("_from_192837465")
			theta0
		} else if (init_supplied == 0) {
			// 		printf("init_supplied is 0\n")
		  xxi       = cross(xit,xit)
			ols_b     = invsym( xxi ) * cross(xit,yit)

			// panel RE
			re_b = st_matrix("_theta_re_192837465")
			ols_b = re_b'



			// exit()


			ols_res   = yit - xit * ols_b
			ols_res_m = J(nobs,1,.)
			for (i = 1; i <= nobs; i++) {
				ols_res_m[i] = mean( panelsubmatrix(ols_res,  i, ids) )
			}
			ols_res_m_abs = abs(ols_res_m)

			olsZv_b     = invsym( cross(zvit, zvit) ) * cross(zvit,(log(ols_res:^2)))
			olsZu_b     = invsym( cross(zuit, zuit) ) * cross(zuit,(log(abs(ols_res):^2)))
			//   olsZv_b
			// zvi
			// zui
			// ols_res_m
			//   olsZv0_b     = invsym( cross(zvi, zvi) ) * cross(zvi,(log(ols_res_m:^2)))
			//   olsZu0_b     = invsym( cross(zui, zui) ) * cross(zui,(log(ols_res_m_abs:^2)))
			olsZv0_b     = invsym( cross(zvi, zvi) ) * cross(zvi,(log(ols_res:^2)))
			olsZu0_b     = invsym( cross(zui, zui) ) * cross(zui,(log(abs(ols_res):^2)))
			//   olsZu_b

			// Kb+Kv0+Ku0+Kvi+Kui

			// 	1 v_0i
			// 	2 u_0i
			// 	3 v_it
			// 	4 u_it

			// 	myScale = 0.25
			// 	myScale2 = 0.25

			// initscaf
			// initscav
			// initscau

			theta0  = initscaf * ols_b \ initscav * olsZv0_b \ initscau * olsZu0_b \ initscav * olsZv_b \ initscau * olsZu_b

			theta0 = theta0'

			// 	theta0
		} else {
			printf("init_supplied is something else: %f\n", init_supplied)
		}

		// exit()


		//   1072

		//   theta0s


		//   1071
		Ktheta    = length(theta0)

		iter = strtoreal(iter0)

		//   iter

		if (model == 41) {

			// 	Halton draws
			// 	108
			// 	st_numscalar("_r_192837465")
			// 	nobs
			// 	st_numscalar("_haltonb_192837465")

			if (mySimtype == 1) {
				// 			printf("mySimtype is 1\n")
				tymV = generate_halton_matrix(st_numscalar("_r_192837465"), nobs, st_numscalar("_haltonb_192837465"))
				tymU = generate_halton_matrix(st_numscalar("_r_192837465"), nobs, st_numscalar("_haltonb_192837465") + st_numscalar("_basesa_192837465"))
			} else if (mySimtype == 2) {
				// 			printf("mySimtype is 2\n")
				tym = runiform(st_numscalar("_r_192837465"), 2*nobs)
				tymV = tym[ , 1..nobs]         // First N columns
				tymU = tym[ , (nobs+1)..(2*nobs)] // Last N columns
			} else {
				printf("simtype is something else: %f\n", mySimtype)
			}

			// 	tymV
			tymV = invnormal(tymV')
			st_matrix("_V_192837465", tymV)

			// 	tymU
			tymU = abs(invnormal(tymU'))
			st_matrix("_U_192837465", tymU)
			// 	stata("matrix list _U_192837465")

			// Create pointers to matrices
			// pY = &yit
			// pX = &xit
			YX = (yit, xit) 
			pYX = &YX    // combined matrix
			pV = &tymV
			pU = &tymU
			pzvi = &zvi
			pzui = &zui
			pzvit = &zvit
			pzuit = &zuit



			//	put into stata matrices foreach

			st_matrix("_y_192837465", yit)
			//   11122
			st_matrix("_x_192837465", xit)
			//   11123
			st_matrix("_zvi_192837465", zvi)
			//   11124
			st_matrix("_zui_192837465", zui)
			//   11123
			st_matrix("_zvit_192837465", zvit)
			//   11124
			st_matrix("_zuit_192837465", zuit)

			st_matrix("_ids_192837465", unique_ids)
			st_matrix("_lnls_192837465", J(nobs, 1, 17))

			// 	111241
			st_matrix("_idvar_192837465", ids0)	
			// 	111242

			//   11126
			st_numscalar("_kb_192837465", k)
			//   11127
			st_numscalar("_kv0_192837465", kv0)
			//   11128
			st_numscalar("_ku0_192837465", ku0)
			//   11129
			st_numscalar("_kv_192837465", kv)
			//   11128
			st_numscalar("_ku_192837465", ku)

			st_numscalar("_nt_192837465", nt)
			st_numscalar("_n_192837465", nobs)
			st_numscalar("_idlenmax_192837465", maxTi)

			st_matrix("_t1_192837465", theta0)
			//   11121
			// length(theta0)
			st_matrix("_grad_192837465", J(nobs, length(theta0), 17))
			// 	111211

			// ## optimize

			// 1233

			S  = optimize_init()
			optimize_init_evaluator(S, &xtsf2optimizeC())
			optimize_init_valueid(S, "simulated log-likelihood")
			optimize_init_iterid(S, "iter")
			optimize_init_evaluatortype(S, "gf0")
			optimize_init_technique(S, technique)
			// 	1234
			optimize_init_nmsimplexdeltas(S,J(1,length(theta0),1))
			//   optimize_init_tracelevel(S, "none")
			//   optimize_init_conv_warning(S, "off")
			//   optimize_init_verbose(S, 0)
			optimize_init_tracelevel(S, tracelevel)
			//   1235
			optimize_init_conv_maxiter(S, iter)
			//   1236
			//   optimize_init_singularHmethod(S, "m-marquardt")
// 			optimize_init_singularHmethod(S, "hybrid")
			//   1237
				optimize_init_deriv_usemin(S, "on")
			// 	1238
			//   optimize_init_argument(S, 1, y)
			// 	optimize_init_argument(S, 2, x)
			// 	optimize_init_argument(S, 3, zsv)
			// 	optimize_init_argument(S, 4, zsk)
			// 	optimize_init_argument(S, 5, zsu)
			// 	optimize_init_argument(S, 6, myscalars)
			optimize_init_params(S, theta0)

			// 	1239

// 			"The Model is MSLE"

		}
		else if (model == 42) {
			"The Model is Full MLE"
			myscalars = k, kv0, ku0, kv, ku, Ktheta, nobs, nt, ///
			st_numscalar("_prod_192837465"), 
			st_numscalar("_r_192837465"),
			st_numscalar("_haltonb_192837465")
			// 		myscalars

			if (mySimtype == 1) {
				printf("mySimtype is 1\n")
				mydraws = generate_halton_matrix(st_numscalar("_r_192837465"), ///
					nt, ///
					st_numscalar("_haltonb_192837465")) // this generates RxNT, I need transpose
			} else if (mySimtype == 2) {
				printf("mySimtype is 2\n")
				mydraws = runiform(st_numscalar("_r_192837465"), nt)
			} else {
				printf("simtype is something else: %f\n", mySimtype)
			}

			mydraws = mydraws'


			S  = optimize_init()
			optimize_init_evaluator(S, &gtrefull2optimize())
			optimize_init_valueid(S, "log-likelihood")
			optimize_init_iterid(S, "iter")
			optimize_init_evaluatortype(S, "gf0")
			optimize_init_technique(S, "nr")
			optimize_init_deriv_usemin(S, "on")
			optimize_init_tracelevel(S, tracelevel)
			optimize_init_conv_maxiter(S, iter)
			//   optimize_init_singularHmethod(S, "hybrid")
			optimize_init_argument(S, 1, yit)
			optimize_init_argument(S, 2, xit)
			optimize_init_argument(S, 3, zvi)
			optimize_init_argument(S, 4, zui)
			optimize_init_argument(S, 5, zvit)
			optimize_init_argument(S, 6, zuit)
			optimize_init_argument(S, 7, ids)
			optimize_init_argument(S, 8, mydraws)
			optimize_init_argument(S, 9, myscalars)
			optimize_init_params(S, theta0)
		}

		if (tracelevel != "none"){
			printf("\n{input:Log-likelihood maximization using} {result:optimize()} \n\n")
		}

		// 		exit()

		//   211
		bh        = optimize(S)
		// 	bh
		vh        = optimize_result_V_oim(S)
		//   hess = optimize_result_Hessian(S)
		//   vh = invsym(-hess)
		if(robust == ""){  
			vh = optimize_result_V(S)
		}
		else {
			vh = optimize_result_V_robust(S)
		}


		//   211
		st_matrix(bsname,        bh)
		// // 	212
		//   //bcov
		st_matrix(vsname,        vh)

		// ## eff

		// 	1 v_0i
		// 	2 u_0i
		// 	3 v_it
		// 	4 u_it

		bet       = bh[1,1..k]
		// 	bet
		gvi0      = bh[(k+1)..(k+kv0)]
		// 	"gvi0"
		// 	gvi0
		gui0      = bh[(k+kv0+1)..(k+kv0+ku0)]
		// 	"gui0"
		// 	gui0
		gvit      = bh[(k+kv0+ku0+1)..(k+kv0+ku0+kv)]
		// 	"gvit"
		// 	gvit
		guit      = bh[(k+kv0+ku0+kv+1)..(k+kv0+ku0+kv+ku)]
		// 	"guit"
		// 	guit


		xb        = xit * bet'
		eit       = yit - xb
		// 	"sv02"
		if(st_numscalar("_vv0zero_192837465") == 1){
			sv02 = J(nt,1,0)
		}
		else {
			sv02 = exp( zvi * gvi0' )
		}

		if(st_numscalar("_vu0zero_192837465") == 1){
			su02 = J(nt,1,0)
		}
		else {
			su02 = exp( zui * gui0' )
		}	
		if(st_numscalar("_vuizero_192837465") == 1){
			sui2 = J(nt,1,0)
		}
		else {
			sui2 = exp( zuit * guit' )
		}
		// 	"svi2"
		svi2 = exp( zvit * gvit' )

		// ### the loop for eff

		// "te_pers_vec = J(nt, 1, .)"
		te_pers_vec = J(nt, 1, .)
		te_tran_vec = J(nt, 1, .)

		for (i = 1; i <= nobs; i++) {
			// 	"i"
			// 	i
			ei      = panelsubmatrix(eit,  i, ids)
			sv02i   = panelsubmatrix(sv02,  i, ids)
			su02i   = panelsubmatrix(su02,  i, ids)
			sui2i   = panelsubmatrix(sui2,  i, ids)
			svi2i   = panelsubmatrix(svi2,  i, ids)

			//     result = e_exp_tu(ei, i, sui2i, svi2i, sv02i[1], su02i[1],
			//                       1, "halton", 12345, 500, 2, epsilon(1), 0)





			// result = condmean_exp(sui2i, svi2i, su02i[1], sv02i[1], ei,
			// 	st_numscalar("_r_192837465"), st_numscalar("_prod_192837465"))			

			result = e_exp_tu(ei, i, sui2i, svi2i, sv02i[1], su02i[1],
			1, "halton", 12345, 
			st_numscalar("_r_192837465"), 
			st_numscalar("_haltonb_192837465"), 
			epsilon(1), 0)	

			Ti = rows(ei)
			te_pers_vec[ids[i,1]..ids[i,2]] = J(Ti, 1, result[1])
			// 		"te_pers_vec"
			// 		result[1]s
			te_tran_vec[ids[i,1]..ids[i,2]] = result[2..(Ti+1)]
			// 		"te_tran_vec"
			// 		(result[2..(Ti+1)])'

	}


		// exit()


		//   220
		//   te_jlms_mean
		//   st_matrix(u1_name,        u_mean)
		//   221
		//   st_matrix(u2_name,        u_mode)
		//   222
		// te_pers_vec
		st_matrix(eff_pers_name,        te_pers_vec)
		//   221
		st_matrix(eff_tran_name,        te_tran_vec)
		//   st_matrix(eff3_name,        te_bc)
		//   223
		//   st_matrix(eff_lb,        te_l)
		//   224
		//   st_matrix(eff_ub,        te_u)
		//   225

		// Diagnostics
		R2        = 1 - cross(eit,eit) / cross(yit,yit)
		R2a       = 1 - (1-R2) * (nt-1) / (nt-k-nobs)

		// Score
		shat2     = variance(eit)
		shat      = sqrt(shat2)
		//   998
		RSS       = cross(eit, eit)
		//   999
		aic       = log((nt-1)/nt*shat2)+1+2*(k+1)/nt
		//   9991
		bic       = log((nt-1)/nt*shat2)+1+(k+1)*log(nt)/nt

		st_numscalar(n1name,     nobs)
		// 	217
		st_numscalar(nt1name,    nt)
		// 	218
		// 	st_numscalar(R2name,     R2)
		// 	219
		// 	st_numscalar(R2adjname,  R2a)
		// 	222
		// 	st_numscalar(aicname,    aic)
		// 	223
		// 	st_numscalar(bicname,    bic)

		//   st_matrix(xbname,        xb)
		// 	236
		// 	st_numscalar(shatname,   shat)
		//   237
		//   st_numscalar(RSSname,    RSS)
		st_numscalar(LLname,     optimize_result_value(S))
		//   238
		st_numscalar(convstatus, optimize_result_converged(S))
		//   239
		st_matrix(rezname,       eit)
		//   240

	}
end

**# xtsf2optimizeC

capture mata mata drop xtsf2optimizeC()
mata:
	void xtsf2optimizeC(real scalar todo, real vector theta,                 ///
		ll, grad, H)
	{ 
		//   stata("mat list _t1_192837465")
		st_matrix("_t1_192837465", theta)
		// 	stata("mat list _t1_192837465")
		// 	stata("mat dir")
		// 	stata("sca dir")
// 		stata(`"plugin call xtsfC, _nthreads_192837465 _prod_192837465 _V_192837465 _U_192837465 _y_192837465 _x_192837465 _zvi_192837465 _vv0zero_192837465 _zui_192837465 _vu0zero_192837465 _zvit_192837465 _zuit_192837465 _vuizero_192837465 _ids_192837465 _idvar_192837465 _t1_192837465 _nt_192837465 _n_192837465 _r_192837465 _idlenmax_192837465 _kb_192837465 _kv0_192837465 _ku0_192837465 _kv_192837465 _ku_192837465 _grad_192837465 _lnls_192837465 "')
		stata(`"plugin call xtsfC, _nthreads_192837465 _prod_192837465 _V_192837465 _U_192837465 _y_192837465 _x_192837465 _zvi_192837465 _vv0zero_192837465 _zui_192837465 _vu0zero_192837465 _zvit_192837465 _zuit_192837465 _vuizero_192837465 _ids_192837465 _idvar_192837465 _t1_192837465 _nt_192837465 _n_192837465 _r_192837465 _idlenmax_192837465 _kb_192837465 _kv0_192837465 _ku0_192837465 _kv_192837465 _ku_192837465 _lnls_192837465 "')
		// 		stata("mat list _lnls_192837465")

		ll    = st_matrix("_lnls_192837465")
		if (todo>=1) {
			grad = st_matrix("_grad_192837465")
			//     "ll"
			//     ll[1..10,]
			//     "gradient"
			//     grad[1..10,]
			//     colsum(grad)
			//     st_matrix("_t1_192837465")
			//     st_matrix("_t1_5_192837465")
		}
	}

end

* Generate Halton draws in RxN matrix format using MATA
* Find N primes starting from prime A using prime checking algorithm


**# is_prime

capture mata mata drop is_prime()
mata:
	// Function to check if a number is prime (equivalent to C code)
	real scalar is_prime(real scalar num) {
    if (num < 2) return(0)
    for (i = 2; i * i <= num; i++) {
			if (mod(num, i) == 0) return(0)
    }
    return(1)
	}
end

**# generate_primes

capture mata mata drop generate_primes()
mata:
	// Function to generate N primes starting from prime A
	real rowvector generate_primes(real scalar A, real scalar N) {
    primes = J(1, N, .)
    count = 0
    candidate = A

    while (count < N) {
			if (is_prime(candidate)) {
				count++
				primes[count] = candidate
			}
			candidate++
    }

    return(primes)
	}

end

**# generate_halton_matrix

capture mata mata drop generate_halton_matrix()
mata:
	// Function to generate RxN Halton matrix
	real matrix generate_halton_matrix(real scalar R, real scalar N, real scalar A) {

    // Generate N primes starting from A
		//     printf("Finding %f primes starting from %f...\n", N, A)
    primes = generate_primes(A, N)

		//     printf("Generated primes: ")
    if (N <= 1) {
			primes'  // Show all if N <= 20
    } else {
			//         printf("First 10: ")
			//         primes[1..10]  // Show first 10 if N > 20
			//         printf("Last 10: ")
			//         primes[(N-9)..N]  // Show last 10 if N > 20
    }

    // Initialize the RxN matrix
    halton_matrix = J(R, N, .)

    // Generate each column using halton(R, 1, prime)
    for (n = 1; n <= N; n++) {
			current_prime = primes[n]

			// Generate R-dimensional Halton sequence with 1 draw using current prime
			halton_column = halton(R, 1, current_prime)

			// Fill the column in our matrix
			halton_matrix[., n] = halton_column

			//         if (n <= 5) printf("Column %f: Using prime %f\n", n, current_prime)
    }

    return(halton_matrix)
	}
end

**# generate_halton_matrix

mata:
	// Set parameters
	R = 10      // Number of rows (dimensions)  
	N = 38     // Number of columns (draws)
	A = 1      // Starting prime

	// Generate the Halton matrix
	// printf("Generating %fx%f Halton matrix starting from prime %f...\n", R, N, A)
	H = generate_halton_matrix(R, N, A)
end

capture program drop xtsfC
capture program drop xtsfll

local os = upper(substr(c(os), 1, 3))
display "`os'"

local mach = upper(substr(c(machine_type), 1, 3))
display "`mach'"

if "`mach'" == "MAC"{
	local mac_mach = upper(substr(c(machine_type), 12, 3))
	display "`mac_mach'" 
}

if "`os'" == "MAC"{
	if "`mac_mach'" == "INT" {
		// 		di "not ready yet"
		display "Mac Intel"
// 		cd "~/research/coding/c/C_in_Stata"
		program xtsfC, plugin using ("xtsf-mac-intel.plugin")

		}
		else if "`mac_mach'" == "SIL" {
			display "Mac Silicon"
// 			cd "~/research/coding/c/C_in_Stata"
			program xtsfC, plugin using ("xtsf-mac-arm.plugin")
				// 		program xtsfll, plugin using ("xtsf-ll-mac.plugin")
			}
			else {
				display "smth is wrong"
			}
		}
		else if "`os'" == "UNI" {
			if "`mach'" == "MAC"{
				if "`mac_mach'" == "INT" {
					display "Mac Intel"
// 					cd "~/research/coding/c/C_in_Stata"
					program xtsfC, plugin using ("xtsf-mac-intel.plugin")
					}
					else if "`mac_mach'" == "SIL" {
						display "Mac Silicon"
// 						cd "~/research/coding/c/C_in_Stata"
						program xtsfC, plugin using ("xtsf-mac-arm.plugin")
							// 		program xtsfll, plugin using ("xtsf-ll-mac.plugin")
						}
						else {
							display "smth is wrong"
						}
					}
					else {
						cd "~/research/coding/c/C_in_Stata"
						capture program xtsfC, plugin using ("xtsf-linux.plugin")
					}
				}
				else if "`os'" == "WIN"{
					display "Windows"
// 					cd "~/research/coding/c/C_in_Stata"
					program xtsfC, plugin using ("xtsf-win.plugin")
				} 

capture mata mata drop e_exp_tu()
mata:

	real vector e_exp_tu(real vector ei, real scalar id, real vector sui2, real vector svi2,
	real scalar sv02, real scalar su02, real scalar prod,
	string scalar simtype, real scalar seed, real scalar R,
	real scalar halton_base, real scalar inv_tol, real scalar print_level) {

		real scalar Ti, do_prod, su02_zero, ghkReLmd, te_max, range_star
		real matrix A, Sig, V_1, Sig_1, Lmd, R_, CholLmd, mydraws
		real vector Re, te_it, te_it_alt, te_resi, result
		real scalar t00, t11, t22, temp

		Ti = rows(ei)
		do_prod = (prod == 1 ? 1 : -1)
		A = -1 * do_prod * (I(Ti) , J(Ti, 1, 1))

		Sig = diag(svi2) + sv02 * J(Ti, Ti, 1)
		su02_zero = (su02 < sqrt(epsilon(1)))
		if (su02_zero) su02 = epsilon(1)
		sui2 = select(sui2, sui2 :> 0) + (sui2 :== 0) :* epsilon(1)

		V_1 = diag(1 :/ (su02 \ sui2))
		Sig_1 = invsym(Sig)
		Lmd = invsym(V_1 + A' * Sig_1 * A)

		// Lmd = Lmd + 1e-8 * I(rows(Lmd))

		// 		"Lmd"
		// 		Lmd
		R_ = Lmd * A' * Sig_1
		Re = R_ * ei
		// 		"Re"
		// 		Re
		//     CholLmd = cholesky(Lmd)'

		// CholLmd = cholesky(Lmd)
		// CholLmd = CholLmd'

		// 		"CholLmd"
		// 		CholLmd

		if (simtype == "halton") {
			//         primes = generate_primes(halton_base, Ti)
			mydraws = generate_halton_matrix(R, Ti, 2)
		} else {
			mydraws = runiform(R, Ti)
		}
		// 		"mydraws[1..10, 1..5]"
		// 		mydraws[1..10, 1..5]

		//     ghkReLmd = ghk(CholLmd, Re, mydraws, R)

		// Calculate GHK for Re and Lambda
		//     from_vec = J(Ti+1, 1, -1e16)  // Approximation of -Inf
		//     to_vec = Re
		//     ghkReLmd = ghk(CholLmd, from_vec, to_vec, mydraws, R)

		// 		"ghkReLmd - ghk"
		// 		ghkReLmd


		ghkReLmd = ghk_Richard_Gates(Lmd, Re, mydraws, R)

		// 		"ghkReLmd - ghk_Richard_Gates"
		// 		ghkReLmd

		te_it = J(Ti + 1, 1, .)

		// 		printf("Re: %g x %g\n", rows(Re), cols(Re))
		// printf("Lmd: %g x %g\n", rows(Lmd), cols(Lmd))
		// printf("CholLmd: %g x %g\n", rows(CholLmd), cols(CholLmd))
		//
		// printf("Re: %g x %g, Lmd[ , 1]: %g x %g\n", rows(Re), cols(Re), rows(Lmd[ , 1]), cols(Lmd[ , 1]))

		te_it = J(Ti + 1, 1, .)

		for (i = 1; i <= Ti + 1; i++) {

			// printf("i = %g, Lmd[ , i]: %g x %g\n", i, rows(Lmd[ , i]), cols(Lmd[ , i]))

			t00 = -Re[i] + 0.5 * Lmd[i, i]
			t11 = exp(t00)
			// 				t11
			// 				t22 = ghk(CholLmd, Re :- Lmd[ , i], mydraws, R)


			// t22 = ghk(CholLmd, from_vec, Re :- vec(Lmd[ , i]), mydraws, R)
			t22 = ghk_Richard_Gates(Lmd, Re :- vec(Lmd[ , i]), mydraws, R)
			// t22

			te_it[i] = t11 :* t22 :/ ghkReLmd
		}
		// 		te_it
		// 		"loop done"

		// 		finite_idx = selectindex(te_it :== te_it :& te_it :!= .)
		// 		"1"
		// 		nonfinite_idx = selectindex(!(te_it :== te_it :& te_it :!= .))
		//		
		// printf("finite_idx = ")
		// finite_idx
		// printf("nonfinite_idx = ")
		// nonfinite_idx
		// "te_it :== te_it"
		// te_it :== te_it
		// "te_it :!= ."
		// te_it :!= .
		// "te_it :== te_it :& te_it :!= ."
		// te_it :== te_it :& te_it :!= .
		// "!(te_it :== te_it :& te_it :!= .)"
		// !(te_it :== te_it :& te_it :!= .)

		// 		"2"
		// 		te_it[nonfinite_idx] = .
		// 		"3"
		te_max = max(te_it)
		// 		"4"
		te_it_alt = (exp(-prod * sqrt(2 * su02 / pi())) \ exp(-prod * sqrt(2 * sui2 :/ pi())))
		// 		"5"

		if (su02_zero) te_it[1] = 1
		// 		"6"
		te_it[2..(Ti+1)] = select(te_it[2..(Ti+1)], sui2 :> 0) + (sui2 :== 0) :* 1
		// "7"
		return(te_it)
	}

end

**# ghk ghk_Richard_Gates

capture mata mata drop ghk_Richard_Gates()
mata:
	// In the code snippet below, assume that V= Σ and that x is a vector of length m containing the upper bounds of integration. 

	real scalar ghk_Richard_Gates(real matrix V, real vector x, real matrix draws, real scalar n) {
		// here "x" is "to", and "n" is "R"
		real scalar m
		real matrix z, a, p

		m = cols(V)

		// 	"n"
		// 	n		
		// 	"x"
		// 	x	
		// 	"V"
		// 	V

		z = J(n,m-1,0)
		// 	11
		a = J(n,1,x[1])
		// 	12
		p = J(n,1,1)
		// 	13
		T = cholesky(V)'
		// 	14
		for (j=1; j<=m; j++) {
			if (j > 1) a = J(n,1,x[j]) - z[,1::(j-1)]*T[1::(j-1),j]
			a = normal(a:/T[j,j])
			p = p:*a
			// 		if (j < m) z[.,j] = invnormal(uniform(n,1):*a)
			if (j < m) z[.,j] = invnormal(draws[ ,j]:*a)
		}
		// 	21
		pr = sum(p)/n
		return(pr)	
	}
end

**# ghk

capture mata mata drop ghk()
mata:
	// GHK simulator function
	real scalar ghk(real matrix C, real vector to, real matrix draws, real scalar R) {
		real scalar K, k
		real matrix z, a1, b1, pi, mysum

		K = cols(C)
		z = J(R, K - 1, 0)
		pi = J(R, 1, 1)

		for (k = 1; k <= K; k++) {
			if (k == 1) {
				a1 = J(R, 1, to[1])
				b1 = J(R, 1, to[1])
			} else {
				mysum = z[ , 1..(k - 1)] * C[k, 1..(k - 1)]'
				a1 = J(R, 1, to[k]) :- mysum
				b1 = J(R, 1, to[k]) :- mysum
			}

			a1 = normal(a1 :/ C[k, k])
			b1 = normal(b1 :/ C[k, k])
			pi = pi :* (b1 :- a1)

			if (k < K) {
				z[ , k] = invnormal(a1 :+ draws[ , k] :* (b1 :- a1))
			}
		}

		return(mean(pi))
	}
end

//# multivariate normal density
capture mata mata drop mvn_density()
capture mata mata drop mvn_density_log()

mata

	real scalar mvn_density(real vector x, real vector mu, real matrix Sigma)
	{
		real scalar k, detSigma, norm_const, exponent
		real vector diff
		real matrix invSigma

		k = rows(x)
		diff = x - mu
		invSigma = invsym(Sigma)
		detSigma = det(Sigma)

		norm_const = 1 / ((2 * pi())^(k/2) * sqrt(detSigma))
		exponent = -0.5 * diff' * invSigma * diff

		return(norm_const * exp(exponent))
	}

	real scalar mvn_density_log(real vector x, real vector mu, real matrix invSigma, real matrix cholSigma, | real scalar log)
	{
		real scalar k, norm_const, exponent
		real vector diff

		if (args() < 4) log = 1  // default to log = true
		k = rows(x)
		diff = x - mu
		//     invSigma = invsym(Sigma)
		//     detSigma = det(Sigma)

		real scalar log_detSigma
		log_detSigma = 2 * sum(ln(diagonal(cholSigma)))


		//     norm_const = -0.5 * k * ln(2 * pi()) - 0.5 * ln(detSigma)
		norm_const = -0.5 * k * ln(2 * pi()) - 0.5 * log_detSigma
		exponent = -0.5 * diff' * invSigma * diff

		if (log) {
			return(norm_const + exponent)
		} else {
			return(exp(norm_const + exponent))
		}
	}


end

**# gtre_ll_het C
capture mata mata drop gtre_ll_het0()
mata:
	real scalar gtre_ll_het0(real scalar prod,
	real matrix V, real matrix U,
	real colvector Y, real matrix X,
	real matrix Zv0, real scalar vv0zero,
	real matrix Zu0, real scalar vu0zero,
	real matrix Zvi,
	real matrix Zui, real scalar vuizero,
	real colvector ids, real colvector idvar,
	real colvector theta,
	real scalar NT, real scalar N, real scalar R,
	real scalar Kb, real scalar Kv0, real scalar Ku0,
	real scalar Kvi, real scalar Kui)
	{
		real colvector YXb, Vv0, Vu0, Vvi, Vui, S, L
		real colvector Pi
		real scalar q, w, i, r, j
		real scalar es, els, Pitr, Pir, lnls

		// Initialize arrays
		YXb = J(NT, 1, 0)
		Vv0 = J(NT, 1, 0)
		Vu0 = J(NT, 1, 0)
		Vvi = J(NT, 1, 0)
		Vui = J(NT, 1, 0)
		S   = J(NT, 1, 0)
		L   = J(NT, 1, 0)
		Pi  = J(N, 1, 0)

		// Calculate YXb = Y - X*b and variances
		for(q = 1; q <= NT; q++) {
			// XB calculation
			YXb[q] = Y[q]
			for(w = 1; w <= Kb; w++) {
				YXb[q] = YXb[q] - X[q, w] * theta[w]
			}

			// Vv0 calculation
			if(vv0zero == 1) {
				Vv0[q] = 0
			}
			else {
				Vv0[q] = 0
				for(w = 1; w <= Kv0; w++) {
					Vv0[q] = Vv0[q] + Zv0[q, w] * theta[Kb + w]
				}
				Vv0[q] = exp(Vv0[q])
			}

			// Vu0 calculation
			if(vu0zero == 1) {
				Vu0[q] = 0
			}
			else {
				Vu0[q] = 0
				for(w = 1; w <= Ku0; w++) {
					Vu0[q] = Vu0[q] + Zu0[q, w] * theta[Kb + Kv0 + w]
				}
				Vu0[q] = exp(Vu0[q])
			}

			// Vvi calculation
			Vvi[q] = 0
			for(w = 1; w <= Kvi; w++) {
				Vvi[q] = Vvi[q] + Zvi[q, w] * theta[Kb + Kv0 + Ku0 + w]
			}
			Vvi[q] = exp(Vvi[q])

			// Vui calculation
			if(vuizero == 1) {
				Vui[q] = 0
			}
			else {
				Vui[q] = 0
				for(w = 1; w <= Kui; w++) {
					Vui[q] = Vui[q] + Zui[q, w] * theta[Kb + Kv0 + Ku0 + Kvi + w]
				}
				Vui[q] = exp(Vui[q])
			}

			// Sum of Vui and Vvi and Lambda_it
			S[q] = Vui[q] + Vvi[q]
			L[q] = sqrt(Vui[q] / Vvi[q])
		}

		// Initialize log-likelihood
		lnls = 0

		// Main likelihood calculation loop
		for(i = 1; i <= N; i++) {
			Pi[i] = 0

			for(r = 1; r <= R; r++) {
				Pir = 1

				for(j = 1; j <= NT; j++) {
					if(idvar[j] == ids[i]) {
                    // Main calculation
                    es = (YXb[j] - sqrt(Vv0[j]) * V[i, r] + sqrt(Vu0[j]) * U[i, r] * prod) * (S[j]^(-0.5))
                    els = prod * es * L[j]

                    // Calculate Pitr using standard normal density and cumulative distribution
                    Pitr = 2 * (S[j]^(-0.5)) * normalden(es) * normal(-els)

                    // Product over t
                    Pir = Pir * Pitr
					}
				}

				// Sum of Pir
				Pi[i] = Pi[i] + Pir
			}

			// Mean of Pir
			Pi[i] = Pi[i] / R
		}

		// Calculate log-likelihood by summing all the logged elements
		for(i = 1; i <= N; i++) {
			lnls = lnls + ln(Pi[i])
		}

		return(lnls)
	}
end

**# gtre_ll_het1
capture mata mata drop gtre_ll_het1()
mata:
	real colvector gtre_ll_het1(real scalar prod,
	real matrix V, real matrix U,
	real colvector Y, real matrix X,
	real matrix zvi, real scalar vv0zero,
	real matrix zui, real scalar vu0zero,
	real matrix zvit,
	real matrix zuit, real scalar vuizero,
	real matrix ids, real colvector idvar,
	real colvector theta,
	real scalar NT, real scalar N, real scalar R,
	real scalar k, real scalar kv0, real scalar ku0,
	real scalar kv, real scalar ku)
	{
		real rowvector bet, gvi0, gui0, gvit, guit
		real colvector eit, sv02, su02, svi2, sui2
		real colvector ei, svi2i, sui2i, S, L
		real scalar sv02i, su02i
		real colvector Pi, lnls
		real scalar i, r, j, t
		real scalar es, els, Pitr, Pir

		// Extract parameter vectors using Mata-specific syntax
		bet  = theta[1..k]'
		gvi0 = theta[(k+1)..(k+kv0)]'
		gui0 = theta[(k+kv0+1)..(k+kv0+ku0)]'
		gvit = theta[(k+kv0+ku0+1)..(k+kv0+ku0+kv)]'
		guit = theta[(k+kv0+ku0+kv+1)..(k+kv0+ku0+kv+ku)]'

		// Calculate residuals
		eit = Y - X * bet

		// Calculate variance components
		if(vv0zero == 1) {
			sv02 = J(NT, 1, 0)
		}
		else {
			sv02 = exp(zvi * gvi0)
		}

		if(vu0zero == 1) {
			su02 = J(NT, 1, 0)
		}
		else {
			su02 = exp(zui * gui0)
		}

		svi2 = exp(zvit * gvit)

		if(vuizero == 1) {
			sui2 = J(NT, 1, 0)
		}
		else {
			sui2 = exp(zuit * guit)
		}

		// Initialize Pi vector and lnls vector
		Pi = J(N, 1, 0)
		lnls = J(N, 1, 0)

		// Main likelihood calculation loop using panelsubmatrix
		for(i = 1; i <= N; i++) {
			// Extract panel-specific data
			ei    = panelsubmatrix(eit, i, ids)
			sv02i = panelsubmatrix(sv02, i, ids)[1,1]
			su02i = panelsubmatrix(su02, i, ids)[1,1]
			svi2i = panelsubmatrix(svi2, i, ids)
			sui2i = panelsubmatrix(sui2, i, ids)

			// Calculate S and L for this panel
			S = sui2i + svi2i
			L = sqrt(sui2i :/ svi2i)

			Pi[i] = 0

			for(r = 1; r <= R; r++) {
				Pir = 1

				// Loop over time periods for this individual
				for(t = 1; t <= rows(ei); t++) {
					// Main calculation
					es = (ei[t] - sqrt(sv02i) * V[i, r] + sqrt(su02i) * U[i, r] * prod) * (S[t]^(-0.5))
					els = prod * es * L[t]

					// Calculate Pitr using standard normal density and cumulative distribution
					Pitr = 2 * (S[t]^(-0.5)) * normalden(es) * normal(-els)

					// Product over t
					Pir = Pir * Pitr
				}

				// Sum of Pir over r
				Pi[i] = Pi[i] + Pir
			}

			// Mean of Pir
			Pi[i] = Pi[i] / R
		}

		// Calculate log-likelihood for each individual
		for(i = 1; i <= N; i++) {
			lnls[i] = ln(Pi[i])
		}

		return(lnls)
	}
end

**# gtre_ll_het pointers
capture mata mata drop gtre_ll_het_p()
mata:
	void gtre_ll_het_p(real scalar todo, real vector theta, ///
		pointer(real matrix) scalar pYX,    // 1
	pointer(real matrix) scalar pzvi,    // 2
	pointer(real matrix) scalar pzui,    // 3
	pointer(real matrix) scalar pzvit,    // 4
	pointer(real matrix) scalar pzuit,    // 5
	pointer(real matrix) scalar pV,    // 6
	pointer(real matrix) scalar pU,    // 7
	real matrix ids,    // 8
	real rowvector scalars, /// 9
		llf, grad, H)
	{
		real rowvector bet, gvi0, gui0, gvit, guit
		real colvector eit, sv02, su02, svi2, sui2
		real colvector ei, svi2i, sui2i, S, L
		real scalar sv02i, su02i
		real colvector Pi
		real scalar i, r, t
		real scalar es, els, Pitr, Pir

		k         = scalars[1]
		kv0       = scalars[2]
		ku0       = scalars[3]
		kv        = scalars[4]
		ku        = scalars[5]
		//   Ktheta    = scalars[6]
		nobs      = scalars[7]
		nt        = scalars[8]
		prod      = scalars[9]
		R         = scalars[10]
		vv0zero         = scalars[11]
		vu0zero         = scalars[12]
		vuizero         = scalars[13]

		// Extract parameter vectors using Mata-specific syntax
		// 	1
		// 	theta
		// 	11
		// 	scalars
		// 	12
		bet  = theta[1,1..k]'
		// 	13
		// 	bet
		gvi0 = theta[(k+1)..(k+kv0)]'
		// 	14
		// 	gvi0
		gui0 = theta[(k+kv0+1)..(k+kv0+ku0)]'
		// 	15
		// 	gui0
		gvit = theta[(k+kv0+ku0+1)..(k+kv0+ku0+kv)]'
		// 	16
		// 	gvit
		guit = theta[(k+kv0+ku0+kv+1)..(k+kv0+ku0+kv+ku)]'
		// 	17
		// 	guit
		// 	1
		// Calculate residuals
		bet1 = (1\ -bet)
		// 	bet1
		// 	11
		// 	cols(*pYX)
		// 	rows(*pYX)
		eit = *pYX * bet1
		//     eit = *pY - *pX * bet
		//     2
		// Calculate variance components
		if(vv0zero == 1) {
			sv02 = J(nt, 1, 0)
		}
		else {
			sv02 = exp(*pzvi * gvi0)
		}
		//     3
		if(vu0zero == 1) {
			su02 = J(nt, 1, 0)
		}
		else {
			su02 = exp(*pzui * gui0)
		}
		//     4
		svi2 = exp(*pzvit * gvit)
		//     5
		if(vuizero == 1) {
			sui2 = J(nt, 1, 0)
		}
		else {
			sui2 = exp(*pzuit * guit)
		}
		//     6
		// Initialize Pi vector and lnls vector
		Pi = J(nobs, 1, 0)
		llf = J(nobs, 1, 0)
		//     7
		// Main likelihood calculation loop using panelsubmatrix
		for(i = 1; i <= nobs; i++) {
			// Extract panel-specific data
			ei    = panelsubmatrix(eit, i, ids)
			sv02i = panelsubmatrix(sv02, i, ids)[1,1]
			su02i = panelsubmatrix(su02, i, ids)[1,1]
			svi2i = panelsubmatrix(svi2, i, ids)
			sui2i = panelsubmatrix(sui2, i, ids)

			// Calculate S and L for this panel
			S = sui2i + svi2i
			L = sqrt(sui2i :/ svi2i)

			Pi[i] = 0

			for(r = 1; r <= R; r++) {
				Pir = 1

				// Loop over time periods for this individual
				for(t = 1; t <= rows(ei); t++) {
					// Main calculation
					es = (ei[t] - sqrt(sv02i) * (*pV)[i, r] + sqrt(su02i) * (*pU)[i, r] * prod) * (S[t]^(-0.5))
					els = prod * es * L[t]

					// Calculate Pitr using standard normal density and cumulative distribution
					Pitr = 2 * (S[t]^(-0.5)) * normalden(es) * normal(-els)

					// Product over t
					Pir = Pir * Pitr
				}

				// Sum of Pir over r
				Pi[i] = Pi[i] + Pir
			}

			// Mean of Pir
			Pi[i] = Pi[i] / R
		}

		// Calculate log-likelihood for each individual
		for(i = 1; i <= nobs; i++) {
			llf[i] = ln(Pi[i])
		}
		// 		printf("llf = %f\n", sum(llf))

	}
end

**# gtre_ll_het no pointers
capture mata mata drop gtre_ll_het()
mata:
	void gtre_ll_het(real scalar todo, real vector theta, ///
		real matrix YX,    // 1
	real matrix zvi,    // 2
	real matrix zui,    // 3
	real matrix zvit,    // 4
	real matrix zuit,    // 5
	real matrix V,    // 6
	real matrix U,    // 7
	real matrix ids,    // 8
	real rowvector scalars, /// 9
		llf, grad, H)
	{
		real rowvector bet, gvi0, gui0, gvit, guit
		real colvector eit, sv02, su02, svi2, sui2
		real colvector ei, svi2i, sui2i, S, L
		real scalar sv02i, su02i
		real colvector Pi
		real scalar i, r, t
		real scalar es, els, Pitr, Pir

		k         = scalars[1]
		kv0       = scalars[2]
		ku0       = scalars[3]
		kv        = scalars[4]
		ku        = scalars[5]
		//   Ktheta    = scalars[6]
		nobs      = scalars[7]
		nt        = scalars[8]
		prod      = scalars[9]
		R         = scalars[10]
		vv0zero         = scalars[11]
		vu0zero         = scalars[12]
		vuizero         = scalars[13]

		// Extract parameter vectors using Mata-specific syntax
		// 	1
		// 	theta
		// 	11
		// 	scalars
		// 	12
		bet  = theta[1,1..k]'
		// 	13
		// 	bet
		gvi0 = theta[(k+1)..(k+kv0)]'
		// 	14
		// 	gvi0
		gui0 = theta[(k+kv0+1)..(k+kv0+ku0)]'
		// 	15
		// 	gui0
		gvit = theta[(k+kv0+ku0+1)..(k+kv0+ku0+kv)]'
		// 	16
		// 	gvit
		guit = theta[(k+kv0+ku0+kv+1)..(k+kv0+ku0+kv+ku)]'
		// 	17
		// 	guit
		// 	1
		// Calculate residuals
		bet1 = (1\ -bet)
		// 	bet1
		// 	11
		// 	cols(YX)
		// 	rows(YX)
		eit = YX * bet1
		//     eit = *pY - *pX * bet
		//     2
		// Calculate variance components
		if(vv0zero == 1) {
			sv02 = J(nt, 1, 0)
		}
		else {
			sv02 = exp(zvi * gvi0)
		}
		//     3
		if(vu0zero == 1) {
			su02 = J(nt, 1, 0)
		}
		else {
			su02 = exp(zui * gui0)
		}
		//     4
		svi2 = exp(zvit * gvit)
		//     5
		if(vuizero == 1) {
			sui2 = J(nt, 1, 0)
		}
		else {
			sui2 = exp(zuit * guit)
		}
		//     6
		// Initialize Pi vector and lnls vector
		Pi = J(nobs, 1, 0)
		llf       = J(nobs,1,.)

		//     7
		// Main likelihood calculation loop using panelsubmatrix
		for(i = 1; i <= nobs; i++) {
			// 		printf("doing i = %f\n", i)

			// Extract panel-specific data
			ei    = panelsubmatrix(eit, i, ids)
			sv02i = panelsubmatrix(sv02, i, ids)[1,1]
			su02i = panelsubmatrix(su02, i, ids)[1,1]
			svi2i = panelsubmatrix(svi2, i, ids)
			sui2i = panelsubmatrix(sui2, i, ids)

			// Calculate S and L for this panel
			S = sui2i + svi2i
			L = sqrt(sui2i :/ svi2i)

			Pi[i] = 0

			for(r = 1; r <= R; r++) {
				Pir = 1

				// Loop over time periods for this individual
				for(t = 1; t <= rows(ei); t++) {
					// Main calculation
					es = (ei[t] - sqrt(sv02i) * V[i, r] + sqrt(su02i) * U[i, r] * prod) * (S[t]^(-0.5))
					els = prod * es * L[t]

					// Calculate Pitr using standard normal density and cumulative distribution
					Pitr = 2 * (S[t]^(-0.5)) * normalden(es) * normal(-els)

					// Product over t
					Pir = Pir * Pitr
				}

				// Sum of Pir over r
				Pi[i] = Pi[i] + Pir
			}

			// Mean of Pir
			Pi[i] = Pi[i] / R
		}

		// 	1717

		// Calculate log-likelihood for each individual
		for(i = 1; i <= nobs; i++) {
			llf[i] = ln(Pi[i])
		}
		// 	printf("llf = %f\n", sum(llf))
		// 	printf("i = %d done \n", i)

	}
end
