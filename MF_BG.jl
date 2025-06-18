include("devPkgs.jl")
using BEAVARs
using TimeSeries
using Parameters
# using DelimitedFiles
using Plots
# using Statistics
# using LinearAlgebra


dataHF_tab, dataLF_tab, varList = BEAVARs.readSpec("bg_julL","data/Specifications_mfvar.xlsx");


out_strct, varSetup,hypSetup = beavar("Blagov2025",dataHF_tab,dataLF_tab,varList,0,n_burn=10,n_save=10);

@unpack M_zsp, store_YY, z_vec, Sm_bit = out_strct

yy1 = dropdims(median(store_YY,dims=3),dims=3);
yy_low =  percentile_mat(store_YY, 0.05; dims=3);
yy_high =  percentile_mat(store_YY, 0.95; dims=3);
plot(yy1[:,1]); plot!(yy_low[:,1]); plot!(yy_high[:,1])
plot(M_zsp*yy1[Sm_bit'])
plot!(z_vec)

Yfor3D = BEAVARs.forecast(out_strct,varSetup);
