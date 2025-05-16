include("devPkgs.jl")
using BEAVARs
using DelimitedFiles
using LinearAlgebra
using TimeSeries
using SparseArrays
using Plots
using Distributions


p = 5; # number of lags
dataM_bg_full = readtimearray("data/bg_julL.csv"; format="dd/mm/yyyy", delim=',')
dataQ_bg_full = readtimearray("data/dataQ_BG.csv"; format="dd/mm/yyyy", delim=',')
varNamesM_full = colnames(dataM_bg_full)

dataM_bg_raw = dataM_bg_full
dataQ_bg_raw = dataQ_bg_full


varNamesM_full = colnames(dataM_bg_full)

varNamesHF = [:survIndustryBG];
varNamesLF = [:gdpBG]
varOrder   = [:gdpBG,:survIndustryBG,:ipBG];
# varNamesHF = varNamesM_full;
varOrder = [:gdpBG; varNamesHF]

# select only the needed data and transform it if needed
dataM_bg_tab = dataM_bg_raw[varNamesHF];
dataQ_bg_tab = percentchange(dataQ_bg_raw[varNamesLF])


dataLF_tab = dataQ_bg_tab;
dataHF_tab = dataM_bg_tab;


# beyond this point we shouldn't need any more input from the user, i.e. we switch from Q and M (if you have monthly and quarterly data) or A and Q (if you have annual and quarterly) to HF and LF


## prepare missing data draw
nsave = 1000;
nburn = 1000;
model_str = "CPZ2024";

setup_str, model_type = makeSetup(varOrder,model_str,p,16,8,nsave,nburn)
# hyperparameter setup
hyper_str=hypChan2020()

#---------------------------------------
store_YY,store_beta, store_Σt, M_zsp, z_vec, Sm_bit = BEAVARs.CPZ_loop!(dataLF_tab,dataHF_tab,varOrder,setup_str,hyper_str);


yy1 = median(store_YY,dims=3);
plot(yy1[:,1])

# plot(M_zsp*YYt[Sm_bit])
plot(M_zsp*yy1[Sm_bit'])
plot!(z_vec)
# #---------------------------------------

