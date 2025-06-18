include("devPkgs.jl")
using BEAVARs
using TimeSeries
using Parameters
# using DelimitedFiles
using Plots
# using Statistics
# using LinearAlgebra


dataNRW_full = readtimearray("data/nrw_2025_05_03.csv"; format="dd/mm/yyyy", delim=',');
varNames = colnames(dataNRW_full)
varList_HF = varNames[1:9];
dataM_raw = dataNRW_full[varList_HF]

dataM_tab = map((timestamp, values) -> (timestamp, log.(values)), dataM_raw[varNames[[1:4;6:9]]]);
dataM_tab = [dataM_tab dataM_raw[:ifo]./100];

dataQ_tab = log.(collapse(dataNRW_full[:gdpNRW], year, first))


dataLF_tab = dataQ_tab[1:end-1];
dataHF_tab = dataM_tab;

# dataLF_tab = diff(dataQ_tab[1:end-1]).*100
# dataHF_tab = diff(dataM_tab).*100

varList   = [colnames(dataLF_tab); colnames(dataHF_tab)];

trans = 1
out_strct, varSetup,hypSetup = beavar("CPZ2024",dataHF_tab,dataLF_tab,varList,trans,n_burn=5,n_save=5,p=5);
out_strct, varSetup,hypSetup = beavar("CPZ2024",dataHF_tab,dataLF_tab,varList,trans,n_burn=5000,n_save=5000,p=5);

@unpack M_zsp, store_YY, z_vec, Sm_bit = out_strct

yy1 = dropdims(median(store_YY,dims=3),dims=3);
yy_low =  percentile_mat(store_YY, 0.05; dims=3);
yy_high =  percentile_mat(store_YY, 0.95; dims=3);
plot(yy1[:,1]); plot!(yy_low[:,1]); plot!(yy_high[:,1])
plot(M_zsp*yy1[Sm_bit'])
plot!(z_vec)

Yfor3D = BEAVARs.forecast(out_strct,varSetup);

fc = dropdims(median(out_strct.store_YY,dims=3),dims=3)

include("fcast_plot.jl")