{smcl}
{right:version 1.3  18Sep2020}
{cmd:help xtsf3gkss}
{hline}

{marker title}{...}
{title:Title}

{p2colset 5 20 22 2}{...} {phang} {bf:xtsf3gkss} {hline 2} KSS estimator of the stochastic frontier model for panel data, where arbitrary temporal heterogeneity is allowed. Unbalanced panels are supported{p_end} {p2colreset}{...}

{marker syntax}{...}
{title:Syntax}

{p 10 17 2} {cmd:xtsf3gkss} {it:{help varlist:depvar}} {it:{help varlist:indepvars}} {ifin}, [{cmd:}{it:{help xtsf3gkss##options:options}}] 

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
{syntab :Cost frontier}
{synopt :{opt cost}}fit cost frontier model; default is production frontier
model{p_end}

{syntab:Smoothing parameter}
{synopt :{opt gr0(#)}}grid begin{p_end}
{synopt :{opt gr1(#)}}grid end{p_end}
{synopt :{opt gri(#)}}grid step/increment{p_end}
{synopt :{opt li(#)}}minimum dimension choice{p_end}
{synopt :{opt ls(#)}}maximum dimension choice{p_end}

{syntab :De-meaning}
{synopt :{opt gmean(#)}}elimination of global mean{p_end}
{synopt :{opt tmean(#)}}elimination of individual effects{p_end}
{synopt :{opt tmean(#)}}elimination of time fixed effects{p_end}

{syntab :Reporting}
{synopt :{opt lev:el(#)}}set confidence level; default as set by set level{p_end}
{synopt :{opt {ul on}nolog{ul off}}}suppress display of a log{p_end}
{synopt :{help ereturn##display_options :{it:display_options}}}further options for displaying output{p_end}
{synoptline}
{p2colreset}{...}

{marker description}{...}
{title:Description}

{pstd} 
{cmd:xtsf3gkss} fits the KSS estimator of the stochastic frontier model for panel data, where arbitrary temporal heterogeneity is allowed. 

{pstd}
It allows using factor variables (see {help fvvarlist}). Unbalanced panels are supported.

{pstd}See {help xtsf3gkss_postestimation:xtsf3gkss postestimation} for features available after estimation.

{marker options}{...}
{title:Options}

{dlgtab:Main}

{phang} {opt cost} specifies that frontier fit a cost frontier model

{dlgtab:Smoothing parameter}

{phang} {opt gr0(#)}, {opt gr1(#)} and {opt gri(#)} specify the lower limit, upper limit, and increment of the grid for search, respectively. The default values are 0.1, 1, and 0.1, which implies that {cmd:xtsf3gkss} will go over 0.1, 0.2, 0.3,..., 1.0. 

{phang} {opt li(#)} and {opt ls(#)} specify the minimum and maximum choice of dimensions.

{dlgtab:De-meaning}

{phang} {opt gmean(#)}}, {opt tmean(#)}}, and/or {opt tmean(#)}} allow elimination of global mean, individual effects, and/or time fixed effects

{dlgtab:Reporting}

{phang} {opt level(#)}; see {helpb estimation options##level():[R] estimation options}.

{phang} {opt nolog} prevents {cmd:xtsf3gkss} from displaying any output on the screen (akin to using {cmd:quietly}).


{title:Example}

{pstd}Load data{p_end}

{phang2}{cmd:. webuse xtfrontier1, clear}{p_end}

{pstd}KSS estimator when optimal the smoothing parameter is chosen on the [0.1 1.0] grid with interval of 0.1. In each grid point, dimentions 3 through 7 are checked.{p_end}

{phang2}{cmd:. xtsf3gkss lnwidgets lnmachines lnworkers t, gr0(0.1) gr1(1.0) gri(0.1) li(3) ls(7) imean tmean level(99) cformat(%9.4f)}{p_end}

{pstd}Using factor variables for interaction{p_end}

{phang2}{cmd:. xtsf3gkss lnwidgets c.lnmachines##c.lnworkers c.t##c.t, gr0(0.1) gr1(1.0) gri(0.1) li(3) ls(7) imean tmean level(99) cformat(%9.4f)}{p_end}

{pstd}List the stored results{p_end}

{phang2}{cmd:. ereturn list}{p_end}

{pstd}Using if{p_end}

{phang2}{cmd:. quietly summarize lnwidgets}{p_end}

{phang2}{cmd:. xtsf3gkss lnwidgets c.lnmachines##c.lnworkers if lnwidgets < r(mean), gr0(0.1) gr1(0.6) gri(0.05) reps(9)}{p_end}

{title:Saved results}

{pstd}
{cmd:xtsf3gkss} saves the following in {cmd:e()}:

{synoptset 20 tabbed}{...}
{p2col 5 20 24 2: Scalars}{p_end}
{synopt:{cmd:e(N)}}number of observations{p_end}
{synopt:{cmd:e(sumTi})}the sum of T_i, i = 1,...,N{p_end}
{synopt:{cmd:e(bandwidth)}}optimal L{p_end}
{synopt:{cmd:e(r2)}}R-squared{p_end}
{synopt:{cmd:e(r2_a)}}R-squared adjusted{p_end}
{synopt:{cmd:e(aic)}}AIC{p_end}
{synopt:{cmd:e(bic)}}BIC{p_end}
{synopt:{cmd:e(cp})}Mallows's Cp = RSS/shat^2-NT+2*(k+1), where k is the number of regressors{p_end}
{synopt:{cmd:e(shat})}standard error of the regression{p_end}
{synopt:{cmd:e(RSS})}RSS{p_end}
{synoptset 20 tabbed}{...} {p2col 5 20 24 2: Macros}{p_end}
{synopt:{cmd:e(cmd)}}{cmd:xtsf3gkss}{p_end}
{synopt:{cmd:e(cmdline)}}command as typed{p_end}
{synopt:{cmd:e(depvar)}}name of dependent variable{p_end}
{synopt:{cmd:e(predict)}}program used to implement {opt predict}{p_end}
{synopt:{cmd:e(properties)}}{opt b V}{p_end}
{synoptset 20 tabbed}{...} {p2col 5 20 24 2: Matrices}{p_end}
{synopt:{cmd:e(b)}}vector of estimated coefficients{p_end}
{synopt:{cmd:e(V)}}estimated variance-covariance matrix{p_end}
{synopt:{cmd:e(residuals)}}residuals{p_end}
{synopt:{cmd:e(alpha)}}Individual effects, sumTi x 1, where sumTi is sum of T_i, i = 1,...,N{p_end}
{synopt:{cmd:e(eff)}}Measures of technical efficiency, sumTi x 1, where sumTi is sum of T_i, i = 1,...,N{p_end}

{synoptset 20 tabbed}{...}{p2col 5 20 24 2: Functions}{p_end}
{synopt:{cmd:e(sample)}}marks estimation sample{p_end}
{p2colreset}{...}

{marker references}{...}
{title:References}

{phang}
Kneip, Alois and Sickles, Robin C. and Song, Wonho 2012, "A new Panel Data Treatment for Heterogeneity In Time Trends," {it:Econometric Theory}, 28, 590–628.

{marker codes}{...}
{title:Codes, Programming Kit}

{pstd} {cmd:xtsf3gkss} is an implementation of the MATLAB code {p_end}

{pmore2}Written by Wonho Song{p_end}
{pmore2}Last Update: September 1, 2017{p_end}
{pmore2}E-mail: whsong@cau.ac.kr, whsong73@hotmail.com{p_end}

{title:Author}

{psee} Oleg Badunenko{p_end}{psee} Brunel University London{p_end}{psee}E-mail: oleg.badunenko@brunel.ac.uk {p_end}

{title:Disclaimer}
 
{pstd} This software is provided "as is" without warranty of any kind, either expressed or implied. The entire risk as to the quality and 
performance of the program is with you. Should the program prove defective, you assume the cost of all necessary servicing, repair or 
correction. In no event will the copyright holders or their employers, or any other party who may modify and/or redistribute this software, 
be liable to you for damages, including any general, special, incidental or consequential damages arising out of the use or inability to 
use the program.{p_end}

{title:Also see}

{p 7 14 2}Help: {help xtsf1g}, {help xtsf2g}, {help xtsf2gbi}, {helpb xtsf3gpss1}, {help xtsf3gpss2}, {helpb xtsf3gpss3} (if installed){p_end}

