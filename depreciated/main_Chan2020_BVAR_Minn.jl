# this must be depreciated with multiple dispatch

nclude("devPkgs.jl")
using BEAVARs
using DelimitedFiles

Y_20 = readdlm("data/FRED_Chan2020_LBA.txt", ',');
YY = Y_20[:,[1,7,8,12]];

setup_str, hyper_str = makeSetup(YY,"Chan2020_LBA_Minn")
store_beta, store_Σ =  Chan2020_LBA_Minn(YY, setup_str,hyper_str)