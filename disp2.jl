include("devPkgs.jl")
using Parameters
using BEAVARs
using TimeSeries
using DelimitedFiles
# using Plots
using Statistics
using LinearAlgebra
using Distributions


YY20 = readdlm("data/FRED_Chan2020_LBA.txt", ',');
YY = YY20; n = size(YY20,2);
# YY = YY20[:,[1,4,5,6]]
@btime out_struct, var_strcut, hyp_strct = beavar("Chan2020iniw",YY,n_save = 10,n_burn=10);
out_struct2, var_strcut, hyp_strct = beavar("Chan2020iniw2",YY,n_save = 1000,n_burn=1000);
@unpack p,nburn,nsave = var_strcut
@unpack nu0 = hyp_strct
p = 4

