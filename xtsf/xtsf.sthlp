{smcl}
{right:version 1.2  23Feb2026}
{cmd:help xtsf}
{hline}

{marker title}{...}
{title:Title}

{p2colset 5 20 22 2}{...} {phang} {bf:xtsf} {hline 2} A four‑component estimator for stochastic frontier panel‑data models. All four error‑term components can include determinants, and the method fully supports unbalanced panels.
{p_end} {p2colreset}{...}

{marker syntax}{...}
{title:Syntax}

{p 10 17 2} {cmd:xtsf} {it:{help varlist:depvar}} {it:{help varlist:indepvars}} {ifin}, [{cmd:}{it:{help xtsf##options:options}}] 

{synoptset 17 tabbed}{...}
{marker Specification}{...}
{synopthdr:Specification}
{synoptline}
{syntab :Model}
{synopt :{it:{help varname:depvars}}}left-hand-side variable{p_end}
{synopt :{it:{help varname:indepvars}}}right-hand-side variables. {it:indepvars} may contain factor variables; see {help fvvarlist}{p_end}

{synoptset 17 tabbed}{...}
{synopthdr :options}
{synoptline}

{syntab:Model}
{synopt:{cmd:noconstant}}suppress constant term{p_end}
{synopt:{cmd:robust}}robust covariance estimator{p_end}
{synopt:{opt cost}}fit cost frontier model; default is production frontier{p_end}
{synopt:{opt model(string)}}specify functional form or model label{p_end}


{syntab:Variance components}
{synopt:{opt uitlnvariance(string)}}determinants of ln var(u_it){p_end}
{synopt:{opt vitlnvariance(string)}}determinants of ln var(v_it){p_end}
{synopt:{opt uilnvariance(string)}}determinants of ln var(u_i){p_end}
{synopt:{opt vilnvariance(string)}}determinants of ln var(v_i){p_end}


{syntab:Initial values}
{synopt:{opt from(string)}}initial values vector{p_end}
{synopt:{opt initscaf(real)}}initial scale factor for a(f); default 1{p_end}
{synopt:{opt initscau(real)}}initial scale factor for u; default 0.25{p_end}
{synopt:{opt initscav(real)}}initial scale factor for v; default 0.25{p_end}

{syntab:Simulation}
{synopt:{opt simtype(string)}}simulation type for simulated likelihood{p_end}
{synopt:{opt r(integer)}}number of draws; default 199{p_end}
{synopt:{opt haltonbase(integer)}}Halton base; default 2{p_end}
{synopt:{opt basesapart(integer)}}skip factor for Halton sequence; default 1{p_end}

{syntab:Optimization}
{synopt:{opt iterate(integer)}}maximum number of iterations; default 102{p_end}
{synopt:{opt tracelevel(string)}}controls amount of tracing output{p_end}
{synopt:{opt nolog}}suppress log of iteration progress{p_end}

{syntab:Reporting}
{synopt:{opt level(#)}}set confidence level; default as set by {cmd:set level}{p_end}
{synopt:{help ereturn##display_options:{it:display_options}}}display options for coefficient table{p_end}

{synoptline}
{p2colreset}{...}

{marker description}{...}
{title:Description}

{pstd} 
{cmd:xtsf} fits the four‑component estimator of the stochastic frontier model for panel data. 

{pstd}
It allows using factor variables (see {help fvvarlist}). Unbalanced panels are supported.

{pstd}See {help xtsf_postestimation:xtsf postestimation} for features available after estimation.

{marker options}{...}
{title:Options}

{dlgtab:Model}

{phang}
{cmd:noconstant} suppresses the constant term in the frontier model.

{phang}
{cmd:robust} requests a robust covariance estimator.

{phang}
{opt cost} specifies that {cmd:xtsf} fit a cost frontier model.
The default is a production frontier.

{phang}
{opt model(string)} specifies the functional form or a named model choice.


{dlgtab:Variance components}

{phang}
{opt uitlnvariance(string)} specifies determinants of ln var(u_it), the time‑varying inefficiency component.

{phang}
{opt vitlnvariance(string)} specifies determinants of ln var(v_it), the time‑varying noise component.

{phang}
{opt uilnvariance(string)} specifies determinants of ln var(u_i), the individual‑specific inefficiency component.

{phang}
{opt vilnvariance(string)} specifies determinants of ln var(v_i), the individual‑specific noise component.


{dlgtab:Initial values}

{phang}
{opt from(string)} supplies an initial values vector for optimization.

{phang}
{opt initscaf(real)} sets the initial scale factor for a(f); the default is 1.

{phang}
{opt initscau(real)} sets the initial scale factor for u; the default is 0.25.

{phang}
{opt initscav(real)} sets the initial scale factor for v; the default is 0.25.


{dlgtab:Simulation}

{phang}
{opt simtype(string)} specifies the type of simulation used in the simulated‑likelihood evaluation.

{phang}
{opt r(integer)} sets the number of simulation draws; the default is 199.

{phang}
{opt haltonbase(integer)} sets the Halton base used in generating quasi‑random sequences; the default is 2.

{phang}
{opt basesapart(integer)} specifies the skip factor (“spacing apart”) for the Halton sequence; the default is 1.


{dlgtab:Multiprocessing}

{phang}
{opt nthreads(integer)} specifies the number of processor threads to use.
The default is {cmd:nthreads(1)} (single‑threaded execution).
Increasing the number of threads enables parallel computation where supported.


{dlgtab:Optimization}

{phang}
{opt iterate(integer)} sets the maximum number of iterations; the default is 102.

{phang}
{opt tracelevel(string)} controls the amount of tracing information printed during optimization.

{phang}
{opt nolog} suppresses the iteration log.


{dlgtab:Reporting}

{phang}
{opt level(#)} sets the confidence level; the default is the current setting of {cmd:set level}.

{phang}
{it:display_options} are standard reporting options for coefficient tables;
see {helpb ereturn##display_options:[R] estimation options}.

{title:Example}

{pstd}Load data{p_end}

{phang2}{cmd:. webuse xtfrontier1, clear}{p_end}

{pstd}A four‑component estimator with default options.{p_end}

{phang2}{cmd:. xtsf lnwidgets lnmachines lnworkers t}{p_end}

{pstd}Using factor variables for interaction{p_end}

{phang2}{cmd:. xtsf lnwidgets c.lnmachines##c.lnworkers c.t##c.t, level(99) cformat(%9.4f)}{p_end}

{pstd}Variance in transient inefficiency{p_end}

{phang2}{cmd:. timer on 1}{p_end}
{phang2}{cmd:. xtsf lnwidgets c.lnmachines##c.lnworkers c.t##c.t, level(99) cformat(%9.4f)}{p_end}
{phang2}{cmd:. timer off 1}{p_end}

{pstd}Multiprocessing{p_end}

{phang2}{cmd:. timer on 2}{p_end}
{phang2}{cmd:. xtsf lnwidgets c.lnmachines##c.lnworkers c.t##c.t, nthreads(10)}{p_end}
{phang2}{cmd:. timer off 2}{p_end}

{phang2}{cmd:. timer list}{p_end}

{pstd}List the stored results{p_end}

{phang2}{cmd:. ereturn list}{p_end}

{title:Saved results}

{pstd}
{cmd:xtsf} saves the following in {cmd:e()}:

{synoptset 20 tabbed}{...}
{p2col 5 20 24 2: Scalars}{p_end}
{synopt:{cmd:e(N)}}number of observations{p_end}
{synopt:{cmd:e(sumTi})}the sum of T_i, i = 1,...,N{p_end}
{synopt:{cmd:e(converged})}1 if converged{p_end}
{synopt:{cmd:e(level})}confidence level{p_end}
{synopt:{cmd:e(ll})}simulated log-likelihood value{p_end}
{synoptset 20 tabbed}{...} {p2col 5 20 24 2: Macros}{p_end}
{synopt:{cmd:e(cmd)}}{cmd:xtsf}{p_end}
{synopt:{cmd:e(cmdline)}}command as typed{p_end}
{synopt:{cmd:e(depvar)}}name of dependent variable{p_end}
{synopt:{cmd:e(predict)}}program used to implement {opt predict}{p_end}
{synopt:{cmd:e(properties)}}{opt b V}{p_end}
{synoptset 20 tabbed}{...} {p2col 5 20 24 2: Matrices}{p_end}
{synopt:{cmd:e(b)}}vector of estimated coefficients{p_end}
{synopt:{cmd:e(V)}}estimated variance-covariance matrix{p_end}
{synopt:{cmd:e(residuals)}}residuals{p_end}
{synopt:{cmd:e(eff_persistent)}}Measures of persistent technical efficiency, sumTi x 1, where sumTi is sum of T_i, i = 1,...,N{p_end}
{synopt:{cmd:e(eff_transient)}}Measures of transient technical efficiency, sumTi x 1, where sumTi is sum of T_i, i = 1,...,N{p_end}

{synoptset 20 tabbed}{...}{p2col 5 20 24 2: Functions}{p_end}
{synopt:{cmd:e(sample)}}marks estimation sample{p_end}
{p2colreset}{...}

{marker references}{...}
{title:References}

{phang}
Oleg Badunenko and Subal C. Kumbhakar (2017) "Economies of Scale, Technical Change and Persistent and Time-Varying Cost Efficiency in Indian Banking: Do Ownership, Regulation and Heterogeneity Matter?" {it:European Journal of Operational Research}, DOI https://doi.org/10.1016/j.ejor.2017.01.025

{marker codes}{...}
{title:Codes, Programming Kit}

{pstd} {cmd:xtsf} is an implementation of the C code {p_end}

{title:Author}

{psee} Oleg Badunenko{p_end}{psee} Brunel University London{p_end}{psee}E-mail: oleg.badunenko@brunel.ac.uk {p_end}

{title:Disclaimer}
 
{pstd} This software is provided "as is" without warranty of any kind, either expressed or implied. The entire risk as to the quality and 
performance of the program is with you. Should the program prove defective, you assume the cost of all necessary servicing, repair or 
correction. In no event will the copyright holders or their employers, or any other party who may modify and/or redistribute this software, 
be liable to you for damages, including any general, special, incidental or consequential damages arising out of the use or inability to 
use the program.{p_end}

{title:Also see}

{p 7 14 2}Help: {help xtsf1g}, {help xtsf2g}, {help xtsf2gbi}, {helpb xtsf3gkss}, {helpb xtsf3gpss1}, {help xtsf3gpss2}, {helpb xtsf3gpss3} (if installed){p_end}

