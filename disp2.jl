include("devPkgs.jl")
using BEAVARs
# using Parameters
using TimeSeries
# using DelimitedFiles
# using Plots
# using Statistics
# using LinearAlgebra
# using Distributions


dataHF_tab, dataLF_tab, varList = BEAVARs.readSpec("bg_L250703","data/Specifications_mfvar.xlsx");

model_type, hyp_strct, set_strct = makeSetup("CPZ2024",n_burn=10;n_save=10)

data_strct = BEAVARs.makeDataSetup(model_type,dataHF_tab, dataLF_tab,0)

# beavar(model_type,dataHF_tab,dataLF_tab,varList,aggWgh, hyp_strct, set_strct);
beavar(model_type,data_strct, hyp_strct, set_strct);
@btime beavar(model_type,data_strct, hyp_strct, set_strct);
@code_warntype beavar(model_type,data_strct, hyp_strct, set_strct);


@code_warntype beavar(model_type,dataHF_tab,dataLF_tab,varList,aggWgh, hyp_strct, set_strct)
@btime beavar(model_type,dataHF_tab,dataLF_tab,varList,aggWgh, hyp_strct, set_strct);






