# xtsf Illustration

In this article, the functionality In this article, the functionality of the command **xtsf** is showcased. The subset of the PWT V10 data is used:

``` stata
use ../../pwt10_sbsmpl, clear
```

# Define specifications

Here are the specifications of the translog and Cobb-Douglas production functions.

``` stata
global spePWTCD "lnY lnK lnHL"
global spePWTTL "lnY (c.lnK c.lnHL)##(c.lnK c.lnHL)"
global CSyear 2000
global itermax = 1000
```

# Translog specification

## Attempt 1

`xtsf` fits the 4-components estimator of the stochastic frontier
model for panel data. It allows using factor variables (see `fvvarlist`). Unbalanced panels
are supported.

The options `initscaf` (scaling factor of the OLS estimates for the *frontier* starting values), `initscau` (scaling factor of the initial value of the variance of $u$), and `initscav` (scaling factor of the initial value of the variance of $v$) should be kept at their default values, but varying them may help convergence. The option `haltonbase` defines the base for the Halton draws for the noise term for the first observation. The option `basesapart` governs how many primes apart the draws for the inefficiency term lie. Even with $R=199$ the estimation is computer intensive

``` stata
timer clear 41
timer on 41
xtsf $spePWTTL if sample_panel, simtype(halton) haltonbase(3) basesapart(2) technique(nr) r(199) iter(102) model(gtreMSLE) initscaf(0.99) initscau(0.5) initscav(0.25)
timer off 41

...
iter 102: Simulated log-likelihood =  881.46667  (not concave)
convergence not achieved

Sample:----------------------
 Number of obs    = 750
 Number of groups = 38
Diagnostics:-----------------
 Simulated logL   = 881.4667       
-----------------------------

 The 4-components estimator of
 the production stochastic frontier model for panel data

-------------------------------------------------------------------------------
          lnY | Coefficient  Std. err.      z    P>|z|     [95% conf. interval]
--------------+----------------------------------------------------------------
          lnK |   2.805423   .0034264   818.76   0.000     2.798707    2.812139
         lnHL |  -1.151636          .        .       .            .           .
              |
  c.lnK#c.lnK |  -.1008305   .0000526 -1915.58   0.000    -.1009337   -.1007274
              |
 c.lnK#c.lnHL |     .14631   .0003042   480.91   0.000     .1457137    .1469063
              |
c.lnHL#c.lnHL |  -.0534775   .0007695   -69.49   0.000    -.0549857   -.0519693
              |
        _cons |  -8.854076          .        .       .            .           .
--------------+----------------------------------------------------------------
ln[var(vi)]   |
        _cons |  -3.330546    .076926   -43.30   0.000    -3.481318   -3.179774
--------------+----------------------------------------------------------------
ln[var(ui)]   |
        _cons |  -7.939136   1.293957    -6.14   0.000    -10.47525   -5.403027
--------------+----------------------------------------------------------------
ln[var(vit)]  |
        _cons |  -5.471688   .0952354   -57.45   0.000    -5.658346    -5.28503
--------------+----------------------------------------------------------------
ln[var(uit)]  |
        _cons |  -9.757716   16.45169    -0.59   0.553    -42.00244    22.48701
-------------------------------------------------------------------------------

. timer clear 1

. timer on 1

. timer list 41
  41:    264.93 /        1 =     264.9260
. display r(t41)/60 " minutes"
4.4154333 minutes
```

## Attempt 2

We can try different scaling value for the variance of $u$ ($0.5$ to $0.25$) and increase the number of permitted iterations to $1020$:


``` stata
timer clear 42
timer on 42
xtsf $spePWTTL if sample_panel, simtype(halton) haltonbase(3) basesapart(2) technique(nr) r(199) iter(1020) model(gtreMSLE) initscaf(0.99) initscau(0.25) initscav(0.25)
timer off 42

...

iter 1017: Simulated log-likelihood =  813.69897  (not concave)
iter 1018: Simulated log-likelihood =  813.69929  (not concave)
iter 1019: Simulated log-likelihood =  813.69961  (not concave)
iter 1020: Simulated log-likelihood =  813.69992  (not concave)
convergence not achieved

Sample:----------------------
 Number of obs    = 750
 Number of groups = 38
Diagnostics:-----------------
 Simulated logL   = 813.6999       
-----------------------------

 The 4-components estimator of
 the production stochastic frontier model for panel data

-------------------------------------------------------------------------------
          lnY | Coefficient  Std. err.      z    P>|z|     [95% conf. interval]
--------------+----------------------------------------------------------------
          lnK |   2.793017   .0018322  1524.39   0.000     2.789426    2.796608
         lnHL |  -1.125678   .0232385   -48.44   0.000    -1.171225   -1.080131
              |
  c.lnK#c.lnK |  -.1003921   .0000654 -1535.14   0.000    -.1005203   -.1002639
              |
 c.lnK#c.lnHL |   .1448535   .0004709   307.64   0.000     .1439307    .1457764
              |
c.lnHL#c.lnHL |  -.0425306   .0025249   -16.84   0.000    -.0474793    -.037582
              |
        _cons |  -9.393699          .        .       .            .           .
--------------+----------------------------------------------------------------
ln[var(vi)]   |
        _cons |  -1.019085   .0326685   -31.19   0.000    -1.083114   -.9550555
--------------+----------------------------------------------------------------
ln[var(ui)]   |
        _cons |  -1.436433   .0690178   -20.81   0.000    -1.571705   -1.301161
--------------+----------------------------------------------------------------
ln[var(vit)]  |
        _cons |  -6.182829   .1626949   -38.00   0.000    -6.501705   -5.863952
--------------+----------------------------------------------------------------
ln[var(uit)]  |
        _cons |  -5.000942   .1657858   -30.17   0.000    -5.325876   -4.676008
-------------------------------------------------------------------------------
\end{lstlisting}
taking prohibitively a lot of time
\begin{lstlisting}[language = stata]
. timer list 42
  42:   2785.68 /        1 =    2785.6850
. display r(t42)/60 " minutes"
46.428083 minutes
```

Unfortunately the convergence is still not achieved. 

## Attempt 3

In this choice of options
```stata
xtsf $spePWTTL if sample_panel, simtype(halton) haltonbase(2) basesapart(1) technique(bfgs) r(199) iter(102) model(gtreMSLE)
\end{lstlisting}
```
we encounter another issue:
> error: "could not calculate numerical derivatives -- flat or discontinuous region encountered"

which could indicate bad initial values that led optimization into flat region. 

## Attempt 4

Possible solution could be to find some starting values with a low $R$:
```stata
xtsf $spePWTTL if sample_panel, simtype(halton) haltonbase(2) basesapart(1) nthreads(10) technique(nr) r(199) iter(100) model(gtreMSLE) from(init_theta)
```
This approach still does not lead to satisfactory results:

> convergence not achieved

# Cobb-Douglas specification

The model seems too complex. (Hint at the possibility of no inefficiency in the PWT data once random effects are taken into account). We wish to draw attention to difficulty estimating the models involving random effects. Estimating a simpler model, the Cobb-Douglas specification is successful and quick

``` stata
timer clear 40
timer on 40
xtsf $spePWTCD if sample_panel, simtype(halton) technique(nr) r(249) iter(102) model(gtreMSLE)
timer off 40
timer list 40
display r(t40)/60 " minutes"

  The Model is MSLE

Log-likelihood maximization using optimize() 

iter 0:  Simulated log-likelihood = -493.98842  (not concave)
iter 1:  Simulated log-likelihood =  543.32339  (not concave)
...
iter 25: Simulated log-likelihood =  822.23872  

Sample:----------------------
 Number of obs    = 750
 Number of groups = 38
Diagnostics:-----------------
 Simulated logL   = 822.2387       
-----------------------------

 The 4-components estimator of
 the production stochastic frontier model for panel data

------------------------------------------------------------------------------
         lnY | Coefficient  Std. err.      z    P>|z|     [95% conf. interval]
-------------+----------------------------------------------------------------
         lnK |    .355902    .005228    68.08   0.000     .3456553    .3661487
        lnHL |   .6764163   .0090621    74.64   0.000      .658655    .6941776
       _cons |   5.915055   .0721282    82.01   0.000     5.773686    6.056423
-------------+----------------------------------------------------------------
ln[var(vi)]  |
       _cons |  -2.707901   .1012418   -26.75   0.000    -2.906332   -2.509471
-------------+----------------------------------------------------------------
ln[var(ui)]  |
       _cons |  -7.861818   1.028573    -7.64   0.000    -9.877785   -5.845851
-------------+----------------------------------------------------------------
ln[var(vit)] |
       _cons |  -5.716886   .2527519   -22.62   0.000    -6.212271   -5.221502
-------------+----------------------------------------------------------------
ln[var(uit)] |
       _cons |  -5.329486   .4937253   -10.79   0.000     -6.29717   -4.361802
------------------------------------------------------------------------------


. timer off 40

. timer list 40
  40:     52.33 /        1 =      52.3260

. display r(t40)/60 " minutes"
.8721 minutes
```

# Efficiencies

We can now obtain efficiencies from the GTRE model using the `predict` command. We start with persistent efficiency. `predict` with the `persistent` option generates efficiencies for all time periods, which is useful to calculate the overall efficiency, however if we wish to look at the summary, which need to count persistent efficiency once per country:

## Persistent efficiency

``` stata
capture drop te_pers
predict te_pers if sample_panel, persistent

* last year available
capture drop year_last_available
bysort id: egen year_last_available = max(year) if sample_panel

capture drop te_pers_initial
generate te_pers_initial = te_pers if year == year_last_available & sample_panel

summarize te_pers_initial

    Variable |        Obs        Mean    Std. dev.       Min        Max
-------------+---------------------------------------------------------
te_pers_in~l |         38    .9825579    .0043489   .9717373   .9910768
```

## Transient efficiency

There is not much persistent inefficiency. This is expected from a low estimate `ln[var(ui)]` of $-7.861818$ implying $\sigma^2_{u0i}=.00038$. We obtain transient efficiency estimates using the `predict` command with the `transient` option:

``` stata
capture drop te_tran
predict te_tran if sample_panel, transient
capture drop fitted
predict fitted if sample_panel

summarize te_pers te_tran

    Variable |        Obs        Mean    Std. dev.       Min        Max
-------------+---------------------------------------------------------
     te_pers |        750    .9825288    .0043116   .9717373   .9910768
     te_tran |        750    .9503249    .0219938   .8470802   .9945518
```

That is it for now.
