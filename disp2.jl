include("devPkgs.jl")
using Parameters
using BEAVARs
using TimeSeries
using DelimitedFiles
# using Plots
using Statistics
using LinearAlgebra
using Distributions
# using Cthulhu


dataHF_tab, dataLF_tab, varList = BEAVARs.readSpec("bg_L250703","data/Specifications_mfvar.xlsx");
var_strct = makeSetup("CPZ2024",dataHF_tab,dataLF_tab,0,n_burn=10;n_save=10)

BEAVARs.beavar2(var_strct);

@code_warntype BEAVARs.beavar2(var_strct);

beavar(var_strct);
@code_warntype beavar(var_strct);


@btime beavar(var_strct);
@btime BEAVARs.beavar2(var_strct);


out_strct, varSetup,hypSetup = beavar("CPZ2024",dataHF_tab,dataLF_tab,varList,0,n_burn=10,n_save=10);

@btime beavar("CPZ2024",dataHF_tab,dataLF_tab,varList,0,n_burn=10,n_save=10);
@code_warntype beavar("CPZ2024",dataHF_tab,dataLF_tab,varList,0,n_burn=10,n_save=10);



set_strct = BEAVARs.VARSetup(4,10,10,16,8,1);
BEAVARs.CPZ2024(dataHF_tab,dataLF_tab,varList,set_strct,hypSetup,0);
@btime BEAVARs.CPZ2024(dataHF_tab,dataLF_tab,varList,set_strct,hypSetup,0);
@code_warntype BEAVARs.CPZ2024(dataHF_tab,dataLF_tab,varList,set_strct,hypSetup,0);


@unpack model_type, p,nburn,nsave,n_irf,n_fcst,dataHF_tab,dataLF_tab, aggMix,hyp_strct = var_strct
    varList = [colnames(dataHF_tab); colnames(dataLF_tab)]
@code_warntype dispatchModel(model_type,dataHF_tab,dataLF_tab,varList,aggMix, hyp_strct, p,nburn,nsave,n_irf,n_fcst)