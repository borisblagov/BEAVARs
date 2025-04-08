include("devPkgs.jl")
using BEAVARs
# using DelimitedFiles
# using TimeSeries
# using Plots
# using Parameters
# using Statistics
# using LinearAlgebra

# data_ta_full = readtimearray("data/data_tpu.csv"; format="dd/mm/yyyy", delim=',')

# data_de = data_ta_full[[:tpuCaldara, :pmiDE, :gdpDE, :invDE, :hicpDE, :euribor]];
# var_names = colnames(data_de)
# YY = values(data_de);

# # YY20 = readdlm("data/FRED_Chan2020_LBA.txt", ',');
# # YY = YY20[:,[1,4,5,6]]
# setup_str, hyper_str, model_type = makeSetup(YY,"Chan2020_LBA_csv",nburn=120,nsave=500,n_fcst=20)
# store_beta, store_h, store_Σ, s2_h_store, store_ρ, store_σ_h2, eh_store, XX = Chan2020_LBA_csv(YY,setup_str,hyper_str);

YY = rand(35,4);
hyp_mine = hypBanbura2010(lambda=0.5)
VAR_str, Hyp_str = beavar(YY,"Chan2020_LBA_csv",p=1,nburn=100,nsave=500);
VAR_str, Hyp_str = beavar(YY,"Chan2020_LBA_Minn",p=1,nburn=100,nsave=500);
VAR_str, Hyp_str = beavar(YY,"Banbura2010",p=1,nburn=100,nsave=500,hyp=hyp_mine);


@time Chan2020_LBA_csv(YY,VAR_str,Hyp_str);
@btime Chan2020_LBA_Minn(YY,VAR_str,Hyp_str);