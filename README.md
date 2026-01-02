# [README: sf-panel](https://olegbadunenko.github.io/sf-panel/index.html)

<a rel="license" href="http://creativecommons.org/licenses/by/4.0/"><img alt="Creative Commons License" style="border-width:0" src="https://i.creativecommons.org/l/by/4.0/88x31.png" /></a><br />This work is licensed under a <a rel="license" href="http://creativecommons.org/licenses/by/4.0/">Creative Commons Attribution 4.0 International License</a>.

> Cite codes usage as shown [here](https://olegbadunenko.github.io/sf-panel/authors.html).

Stata commands to estimate the panel data stochastic frontier models of
the first and second generation written and maintained by Oleg Badunenko
(<oleg.badunenko@brunel.ac.uk>). The details are discussed in

Stochastic frontier analysis in Stata: using existing and coding new
commands in “Handbook of Research Methods and Applications in Empirical
Microeconomics” (edited by Nigar Hashimzade and Michael A. Thornton),
2021, Chapter 17, ***Edward Elgar Publishing***, [DOI
<img src="man/figures/doi.png"  width="12" height="12">](https://doi.org/10.4337/9781788976480.00027)

> Codes to replicate the results in the chapter can be found [here](https://olegbadunenko.github.io/sf-panel/book-examples.html).


# Installation

First, the tool to install from GitHub is required. It can be installed by typing in the Stata command line:

```stata
net install github, from("https://haghish.github.io/github/")
```

typical search is performed by typing

```stata
github search sf
```
which provides tools for installation. The package can be installed directly by typing

```stata
github install OlegBadunenko/sf-panel
```

The package can be installed without the *github* command:

```stata
net install xtsf3g, from("https://raw.githubusercontent.com/OlegBadunenko/sf-panel/main")
```

# Help

Help files after installation can be called as usual:

```stata
help xtsf
```

# Uses

This [article](https://olegbadunenko.github.io/sf-panel/book-examples.html) guides through the code to replicate the results presented in the chapter.
