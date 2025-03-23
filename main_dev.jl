include("devPkgs.jl")
using(BEAVARs)
using DelimitedFiles
using LinearAlgebra
using Statistics

YY = readdlm("data/FRED_Chan2020_LBA.txt",',');
Y_small = YY[:,[1,2,4,8]];
p = 1;
n = 4;
intercept = 1;
n_irf = 5;
nsave = 500;
store_beta, store_sigma = Banbura2010(Y_small,lags=p,nsave=nsave);

includet("src/irfs.jl")

IRF_median, IRF_68_low, IRF_68_high, IRF_4d = irf_chol_overDraws(store_beta,store_sigma,n,p,intercept,n_irf,nsave);

