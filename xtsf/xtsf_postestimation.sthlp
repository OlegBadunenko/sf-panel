{smcl}
{right:version 1.2  18Sep2020}
{cmd:help xtsf3gkss} postestimation
{hline}

{marker title}{...}
{title:Title}

{p2colset 1 32 34 2}{...}
{phang}{bf:xtsf3gkss postestimation} {hline 2} Postestimation tools for {help xtsf3gkss:xtsf3gkss}{p_end} {p2colreset}{...} 


{marker description}{...}
{title:Postestimation commands}

{pstd}
The following postestimation commands are available after {opt xtsf3gkss}:

{synoptset 17}{...}
{p2coldent :Command}Description{p_end}
{synoptline}
{synopt :{helpb lincom}}Linear combinations of parameters{p_end}
{synopt :{helpb nlcom}}Nonlinear combinations of estimators{p_end}
{synopt :{helpb predictnl}}Obtain nonlinear predictions, standard errors, etc., after estimation{p_end}
{synopt :{helpb test}}Test linear hypotheses after estimation{p_end}
{synopt :{helpb testnl}}Test nonlinear hypotheses after estimation{p_end}
{synopt :{helpb xtsf3gkss postestimation##predict:predict}}predictions, residuals, efficiency measures, effects{p_end}
{synoptline}
{p2colreset}{...}

{marker syntax_predict}{...}
{marker predict}{...}
{title:Syntax for predict}

{p 8 16 2}
{cmd:predict}
{newvar}
{ifin}
[{cmd:,} {it:statistic}]

{synoptset 17 tabbed}{...}
{synopthdr :statistic}
{synoptline}
{syntab :Main}
{synopt :{opt xb}}linear prediction; the default{p_end}
{synopt :{opt resid:uals}}residuals{p_end}
{synopt :{opt te}}estimates of the time-constant efficiency{p_end}
{synopt :{opt alpha}}estimates of the time-constant effects, alpha{p_end}
{synoptline}
{p2colreset}{...}
{p 4 6 2}
These statistics are available only in sample; by default {cmd:e(sample)} is used, which defines the estimation sample.


{marker des_predict}{...}
{title:Description for predict}

{pstd}
{cmd:predict} creates a new variable containing predictions such as
linear predictions, residuals, and estimates of technical efficiency.


{marker options_predict}{...}
{title:Options for predict}

{dlgtab:Main}

{phang}
{opt xb}, calculates the linear prediction.

{phang}
{opt resid:uals} calculates the residuals.

{phang}
{opt te} produces estimates of time-constant efficiency measures for both within and GLS models.

{phang}
{opt alpha} produces estimates of time-constant individual effects for both within and GLS models.

{marker examples}{...}
{title:Examples}

{pstd}Load data{p_end}

{phang2}{cmd:. webuse xtfrontier1, clear}{p_end}

{pstd}KSS estimator when optimal the smoothing parameter is chosen on the [0.1 1.0] grid with interval of 0.1. In each grid point, dimentions 3 through 7 are checked.{p_end}

{phang2}{cmd:. xtsf3gkss lnwidgets lnmachines lnworkers t, gr0(0.1) gr1(1.0) gri(0.1) li(3) ls(7) imean tmean level(99) cformat(%9.4f)}{p_end}

{pstd}Using factor variables for interaction{p_end}

{phang2}{cmd:. xtsf3gkss lnwidgets c.lnmachines##c.lnworkers c.t##c.t, gr0(0.1) gr1(1.0) gri(0.1) li(3) ls(7) imean tmean level(99) cformat(%9.4f)}{p_end}

{pstd}Estimate technical efficiency{p_end}
{phang2}{cmd:. predict efficiency, te}{p_end}

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


