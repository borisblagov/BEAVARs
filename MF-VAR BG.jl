include("devPkgs.jl")
using BEAVARs
using TimeSeries
using Parameters
using Plots
# using DelimitedFiles
# using Statistics
# using LinearAlgebra
# using Distributions
# using SparseArrays

include("src/plot_functions.jl")
trans = 0;
dataHF_tab, dataLF_tab, varList = BEAVARs.readSpec("bg_julL","data/Specifications_mfvar.xlsx");
out_strct, varSetup,hypSetup = beavar("Blagov2025",dataHF_tab,dataLF_tab,varList,0,n_burn=100,n_save=100);
fanChart(out_strct.store_YY[:,1,:])
Yfor3D = BEAVARs.forecast(out_strct,varSetup);
fanChart(Yfor3D[:,1,:])


YY_HF_med = percentile_mat(out_strct.store_YY,0.5,dims=3);
@unpack M_zsp, store_YY, z_vec, Sm_bit = out_strct
plot(M_zsp*YY_HF_med'[Sm_bit])
plot!(z_vec)

