# BEAVARs.jl Documentation

This is a personal package implementing various Bayesian VARs for economic analysis and forecasting. 

Current VARs:
 - Bayesian VAR implemented with dummy observations as in [Banbura, M., Giannone, D., and Reichlin, L. (2010), Large Bayesian vecotr auto regressions, _Journal of Applied Econometrics_, Vol 1(1), doi.org/10.1002/jae.1137](https://doi.org/10.1002/jae.1137).
 - Bayesian VAR with Minnesota prior and fixed variance-covariance matrix as in [Chan, J.C.C. (2020), Large Bayesian Vecotrautoregressions, P. Fuleky (Eds), _Macroeconomic Forecasting in the Era of Big Data_, 95-125, Springer, Cham](https://joshuachan.org/papers/large_BVAR.pdf), see also [joshuachan.org](https://joshuachan.org) 
 - Bayesian VAR with Minnesota prior and common stochastic volatility (CSV) as in [Chan, J.C.C. (2020), Large Bayesian Vecotrautoregressions, P. Fuleky (Eds), _Macroeconomic Forecasting in the Era of Big Data_, 95-125, Springer, Cham](https://joshuachan.org/papers/large_BVAR.pdf), see also [joshuachan.org](https://joshuachan.org)

The codes can be optimized.

```@contents
```


# References
Banbura, Marta, Giannone, Domenico and Reichlin, Lucrezia, (2010), Large Bayesian vector auto regressions, _Journal of Applied Econometrics_, 25, issue 1, p. 71-92. [https://doi.org/10.1002/jae.1137](https://doi.org/10.1002/jae.1137). 