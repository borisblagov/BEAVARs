include("devPkgs.jl")

using BEAVARs
# using DelimitedFiles
using TimeSeries
using Parameters
# using CSV
# using Dates
 using Plots

 data_ta_full = readtimearray("data/data_tpu.csv"; format="dd/mm/yyyy", delim=',')

data_de = data_ta_full[[:tpuCaldara, :pmiDE, :gdpDE, :invDE, :hicpDE, :euribor]];
var_names = colnames(data_de)
YY = values(data_de);


