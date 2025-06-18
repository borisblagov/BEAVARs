include("devPkgs.jl")

using BEAVARs
using TimeSeries
using Parameters
# using CSV
# using Dates
 using Plots
 using XLSX

#  data_ta_full = readtimearray("data/data_tpu.csv"; format="dd/mm/yyyy", delim=',')
# data_de = data_ta_full[[:tpuCaldara, :pmiDE, :gdpDE, :invDE, :hicpDE, :euribor]];
# var_names = colnames(data_de)
# YY = values(data_de);


    
modelstr = "bg_julL";                                    # chosen modelstr
data_path = "data/Specifications_mfvar.xlsx"

hf_ta, lf_ta, varListF = BEAVARs.readSpec("bg_julL","data/Specifications_mfvar.xlsx")
