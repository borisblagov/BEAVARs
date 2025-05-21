include("devPkgs.jl")
using BEAVARs
using TimeSeries
# using DelimitedFiles
# using Plots
 using Parameters
# using Statistics
# using LinearAlgebra


YY = rand(35,4);
YYlist = [:Y1; :Y2; :Y3; :Y4]
out_strct, varSetup,hypSetup = beavar("BGR2010",YY)

# YY_tup, model_type, hyp_strct = BEAVARs.beavar2("BGR2010",YY)
# YY_tup, model_type, hyp_strct, set_strct, store_beta, store_sigma = BEAVARs.beavar2("CPZ2024",dataHF_tab,dataLF_tab,varList,p=1,n_burn=100,n_save=100)

# VAR_str, Hyp_str = beavar(YY,"Chan2020_LBA_csv",p=1,nburn=100,nsave=500);
# VAR_str, Hyp_str = beavar(YY,"Chan2020_LBA_Minn",p=1,nburn=100,nsave=500);


# hyp_mine = hypBGR2010(lambda=0.5)
# VAR_str, Hyp_str = beavar(YY,"BGR2010",p=1,nburn=100,nsave=500,hyp=hyp_mine);



dataM_bg_full = readtimearray("data/bg_julL.csv"; format="dd/mm/yyyy", delim=',')
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

@time out_strct, varSetup,hypSetup = beavar("CPZ2024",dataHF_tab,dataLF_tab,varList,n_burn=2000,n_save=2000);

@unpack M_zsp, store_YY, z_vec, Sm_bit = out_strct

yy1 = dropdims(median(store_YY,dims=3),dims=3);
yy_low =  percentile_mat(store_YY, 0.05; dims=3);
yy_high =  percentile_mat(store_YY, 0.95; dims=3);
plot(yy1[:,1]); plot!(yy_low[:,1]); plot!(yy_high[:,1])
plot(M_zsp*yy1[Sm_bit'])
plot!(z_vec)