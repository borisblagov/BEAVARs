include("devPkgs.jl")
using BEAVARs
# using DelimitedFiles
using TimeSeries
# using Plots
# using Parameters
# using Statistics
# using LinearAlgebra


YY = rand(35,4);

YYinp, model_type, nargs, hyp_strct, var_strct = BEAVARs.beavar2("CPZ2024",YY)

VAR_str, Hyp_str = beavar(YY,"Chan2020_LBA_csv",p=1,nburn=100,nsave=500);
VAR_str, Hyp_str = beavar(YY,"Chan2020_LBA_Minn",p=1,nburn=100,nsave=500);


hyp_mine = hypBGR2010(lambda=0.5)
VAR_str, Hyp_str = beavar(YY,"BGR2010",p=1,nburn=100,nsave=500,hyp=hyp_mine);



dataM_bg_full = readtimearray("data/bg_julL.csv"; format="dd/mm/yyyy", delim=',')
dataQ_bg_full = readtimearray("data/dataQ_BG.csv"; format="dd/mm/yyyy", delim=',')
varNamesM_full = colnames(dataM_bg_full)

dataM_bg_raw = dataM_bg_full
dataQ_bg_raw = dataQ_bg_full


varNamesM_full = colnames(dataM_bg_full)
varNamesHF = [:survIndustryBG];
varNamesLF = [:gdpBG]
varList   = [:gdpBG,:survIndustryBG,:ipBG];
# varOrder = [:gdpBG; varNamesHF]

# select only the needed data and transform it if needed
dataM_bg_tab = dataM_bg_raw[varNamesHF];
dataQ_bg_tab = percentchange(dataQ_bg_raw[varNamesLF])

dataLF_tab = dataQ_bg_tab;
dataHF_tab = dataM_bg_tab;

YYtup = (dataHF_tab,dataLF_tab,varList)

YYtest = BEAVARs.beavar2("CPZ2024",YY)