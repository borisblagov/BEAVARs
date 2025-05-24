include("devPkgs.jl")
using BEAVARs
using TimeSeries
using Parameters
# using DelimitedFiles
using Plots
# using Statistics
# using LinearAlgebra



dataM_raw = readtimearray("data/YM_de_opt.csv"; format="dd/mm/yyyy", delim=',');
varList_HF = colnames(dataM_raw)

dataM1_tab = map((timestamp, values) -> (timestamp, log.(values)), dataM_raw[varList_HF[[1:5;7:8;10]]]);
dataM_tab = [dataM1_tab dataM_raw[varList_HF[[6;9]]]];

dataQ_raw = log.(readtimearray("data/YQ_de_opt.csv"; format="dd/mm/yyyy", delim=','))
dataQ_tab = dataQ_raw[:ygdp]


dataLF_tab = dataQ_tab;
dataHF_tab = dataM_tab;

# dataLF_tab = diff(dataQ_tab[1:end-1]).*100
# dataHF_tab = diff(dataM_tab).*100

varList   = [colnames(dataLF_tab); colnames(dataHF_tab)];

trans = 1
out_strct, varSetup,hypSetup = beavar("CPZ2024",dataHF_tab,dataLF_tab,varList,trans,n_burn=50,n_save=50);

@unpack M_zsp, store_YY, z_vec, Sm_bit = out_strct

yy1 = dropdims(median(store_YY,dims=3),dims=3);
yy_low =  percentile_mat(store_YY, 0.05; dims=3);
yy_high =  percentile_mat(store_YY, 0.95; dims=3);
plot(yy1[:,1]); plot!(yy_low[:,1]); plot!(yy_high[:,1])
plot(M_zsp*yy1[Sm_bit'])
plot!(z_vec)

Yfor3D = BEAVARs.forecast(out_strct,varSetup);
include("fcast_plot.jl")




YY = dropdims(median(out_strct.store_YY,dims=3),dims=3);
out_strct, varSetup,hypSetup = beavar("Chan2020_LBA_Minn",YY)
Yfor3D = BEAVARs.forecast(out_strct,varSetup);
include("fcast_plot.jl")