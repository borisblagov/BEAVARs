include("devPkgs.jl")
using BEAVARs
using TimeSeries
using Parameters
# using DelimitedFiles
using Plots
# using Statistics
# using LinearAlgebra


YY = rand(35,4);
YYlist = [:Y1; :Y2; :Y3; :Y4]


# out_strct, varSetup,hypSetup = 
out_strct, varSetup,hypSetup = beavar("Chan2020minn",YY,n_save=10,n_burn=10);
out_strct, varSetup,hypSetup = beavar("Chan2020iniw",YY,n_save=10,n_burn=10);
out_strct, varSetup,hypSetup = beavar("Chan2020csv",YY,n_save=10,n_burn=10);
out_strct2, varSetup2,hypSetup2 = beavar("Chan2020csv",YY,n_save=10,n_burn=10);
out_strct2, varSetup2,hypSetup2 = beavar("BGR2010",YY,n_save=1000,n_burn=5000);

Yfor3D = BEAVARs.forecast(out_strct2,varSetup2)


Yfor_low1 = percentile_mat(Yfor3D,0.05,dims=3);
Yfor_low = percentile_mat(Yfor3D,0.16,dims=3);
Yfor_med = percentile_mat(Yfor3D,0.5,dims=3)
Yfor_hih = percentile_mat(Yfor3D,0.84,dims=3);
Yfor_hih1 = percentile_mat(Yfor3D,0.97,dims=3);

ik = 1
plot(Yfor_med[:,ik],w=0;ribbon=(Yfor_med[:,ik]-Yfor_low1[:,ik],Yfor_hih1[:,ik]-Yfor_med[:,ik]),fillalpha = 0.1,color=1,legend=false)
plot!(Yfor_med[:,ik],w=2;ribbon = (Yfor_med[:,ik]-Yfor_low[:,ik],Yfor_hih[:,ik]-Yfor_med[:,ik]),fillalpha=0.05,color=1)

dataM_bg_full = readtimearray("data/bg_julL.csv"; format="dd/mm/yyyy", delim=',');
dataQ_bg_full = readtimearray("data/dataQ_BG.csv"; format="dd/mm/yyyy", delim=',')
varNamesM_full = colnames(dataM_bg_full)

dataM_bg_raw = dataM_bg_full
dataQ_bg_raw = dataQ_bg_full


varNamesM_full = colnames(dataM_bg_full)
varNamesHF = [:survIndustryBG];
varNamesLF = [:gdpBG]
varList   = [varNamesLF; varNamesHF];

# select only the needed data and transform it if needed
dataM_bg_tab = dataM_bg_raw[varNamesHF];
dataQ_bg_tab = percentchange(dataQ_bg_raw[varNamesLF])

dataLF_tab = dataQ_bg_tab;
dataHF_tab = dataM_bg_tab;

out_strct, varSetup,hypSetup = beavar("CPZ2024",dataHF_tab,dataLF_tab,varList,0,n_burn=200,n_save=100);

@unpack M_zsp, store_YY, z_vec, Sm_bit = out_strct

yy1 = dropdims(median(store_YY,dims=3),dims=3);
yy_low =  percentile_mat(store_YY, 0.05; dims=3);
yy_high =  percentile_mat(store_YY, 0.95; dims=3);
plot(yy1[:,1]); plot!(yy_low[:,1]); plot!(yy_high[:,1])
plot(M_zsp*yy1[Sm_bit'])
plot!(z_vec)

Yfor3D = BEAVARs.forecast(out_strct,varSetup);
