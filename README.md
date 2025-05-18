# BEAVARs.jl: Bayesian Econometric Analysis with VAR models

![Credit Mikelde Ferro: https://www.youtube.com/watch?v=WIYQWK4pkqg&ab_channel=MikedelFerro-Music](BeaverMunchingCabbage.png)
Still taken from a video by [Mikel Ferro](https://www.youtube.com/watch?v=WIYQWK4pkqg&ab_channel=MikedelFerro-Music)

This is a personal package implementing various Bayesian VARs for economic analysis and forecasting. 



Currently implemented VARs:
 - BVAR with dummy observations as in Banbura, M., Giannone, D., and Reichlin, L. (2010), Large Bayesian vector auto regressions, _Journal of Applied Econometrics_, Vol 25(1), [doi.org/10.1002/jae.1137](https://doi.org/10.1002/jae.1137).

 - BVAR with classical  Minnesota prior (homoscedastic fixed variance-covariance matrix) as in Chan, J.C.C. (2020), Large Bayesian Vecotrautoregressions, P. Fuleky (Eds), _Macroeconomic Forecasting in the Era of Big Data_, 95-125, Springer, Cham, [https://doi.org/10.1007/978-3-030-31150-6](https://doi.org/10.1007/978-3-030-31150-6), see also [joshuachan.org](https://joshuachan.org) and his [pdf](https://joshuachan.org/papers/large_BVAR.pdf).

 - BVAR with Minnesota prior and common stochastic volatility (csv) as in Chan, J.C.C. (2020), Large Bayesian Vecotrautoregressions, P. Fuleky (Eds), _Macroeconomic Forecasting in the Era of Big Data_, 95-125, Springer, Cham, [https://doi.org/10.1007/978-3-030-31150-6](https://doi.org/10.1007/978-3-030-31150-6), see also [joshuachan.org](https://joshuachan.org) and his [pdf](https://joshuachan.org/papers/large_BVAR.pdf).

Each model is contained in a separate function. See the documetnation for details. Note that notation follows the original reference. Consequently variable and parameter names are different across functions (e.g. $\lambda_1$ in one paper can be $c_1$ in another). 

Some codes have been translated from Matlab, so there is a lot of room for optimization. 


# References
Banbura, Marta, Giannone, Domenico and Reichlin, Lucrezia, (2010), Large Bayesian vector auto regressions, _Journal of Applied Econometrics_, 25, issue 1, p. 71-92. [https://doi.org/10.1002/jae.1137](https://doi.org/10.1002/jae.1137). 

Chan, J.C.C. (2020), Large Bayesian Vecotrautoregressions, P. Fuleky (Eds), _Macroeconomic Forecasting in the Era of Big Data_, 95-125, Springer, Cham, [https://doi.org/10.1007/978-3-030-31150-6](https://doi.org/10.1007/978-3-030-31150-6)

# Requests
This package is begin developed for my personal use mainly to learn Julia in my sparse spare time and I cannot add a particular paper or method. If you choose to do so and open an issue it might take a while for me to review it but feel free to do so.

# Notes on the name
The name BEAVARs is an obvious play of words with a misspelled version of my favourite animal.

It is also a nod to the [BEAR Toolbox](https://www.ecb.europa.eu/press/research-publications/working-papers/html/bear-toolbox.en.html) - Bayesian  Estimation, Analysis and Regression, which is a powerful Matlab toolbox for estimating various VAR, BVAR, and Panel VAR models. This is not an attempt to reach the size and scope of BEAR in the Julia ecosystem. There are, however, some similarities in the idea of easy estimation of different models.

It does not conform to the widely accepted convention of naming Julia packages (capital letter followed by all lowercase) but it [doesn't break any rules either](https://pkgdocs.julialang.org/v1/creating-packages/#Package-naming-rules). It isn't the only package with more than one capital letter, e.g. FFT, CUDA, CSV etc. Also, because it's misspelled on purpose, the name Beaver.jl remains open, and if someone does that one could still distinguish the packages BEAVARs.jl and Beavers.jl easily.


#  Acknowledgmenets
I would like to thank [Guillaume Dalle](https://github.com/gdalle]). He is not associated with this package but went out of his way to help me get my first steps in Github and Julia optimization. Also, many users in the Julia discourse helped me when I was struggling - this community is great.

