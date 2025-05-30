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
xf = XLSX.readxlsx("data/Specifications_mfvar.xlsx");
sh_names = XLSX.sheetnames(xf)

sh_ref = xf["setup"];
sh_mat = sh_ref[:]

model_ind = findall(sh_mat[1,:].==modelstr)[1]
vb_ind = findall(sh_mat[:,1].=="lastRow")[1] + 1;       # where the variables start
varListA_str = sh_mat[vb_ind:end,1];                    # strings of variables
varListA_sym = Symbol.(varListA_str[:]);                    # symbols of variables


vm_bit = .!iszero.(sh_mat[vb_ind:end,model_ind]);       # boolean list of ALL variables
vm_trans = sh_mat[vb_ind:end,model_ind];                # list of variables with transformations for ALL vars

varListF = varListA_sym[vm_bit]
trans_dictA = Dict(varListA_sym .=> vm_trans);          # transformation dictionary vector for ALL vars. 

datasheet_str = "datasheet" * string(sh_mat[sh_mat[:,1].=="datasheet",model_ind][1]) * "_M";    # the high-frequency datasheet
data_mat = xf[datasheet_str][:];

nan_mat =  fill(NaN, size(data_mat,1),size(data_mat,2));
data_mat[ismissing.(data_mat)]=nan_mat[ismissing.(data_mat)];   # replace all missing with NaN
values_mat = convert(Array{Float64},data_mat[2:end,2:end])      # convert values to numbers
date_vec = DateTime.(data_mat[2:end,1])                         # convert first column to DateTime
# date_vecstr = Dates.DateFormat.(date_vec,"yyyy-mm-dd");

dataf_TA = TimeArray(date_vec,values_mat,Symbol.(data_mat[1,2:end]))

varList_HF = intersect(varListF,colnames(dataf_TA));      # looks for which variables are required and which are found
dataHF_TA = dataf_TA[varList_HF];                         # selects the variables found in this TA

transHF_vec = Vector{Int}();                                # vector of transformations for the found high-freq. vars
push!(transHF_vec,(trans_dictA[i] for i in varList_HF)...);  # fill it


# transformations

# if 1, leave it be, initialize TA
hf_ta = map((timestamp, values) -> (timestamp, (values)), dataHF_TA[varList_HF[transHF_vec.==1]])
#if 2, divide by 100
hf_ta_temp = map((timestamp, values) -> (timestamp, (values)./100.0), dataHF_TA[varList_HF[transHF_vec.==2]])
hf_ta = merge(hf_ta,hf_ta_temp)
# if 3, take logs
hf_ta_temp = map((timestamp, values) -> (timestamp, log.(values)), dataHF_TA[varList_HF[transHF_vec.==3]])
hf_ta = merge(hf_ta,hf_ta_temp)
# if 7 take percentage changes
hf_ta_temp = percentchange(dataHF_TA[varList_HF[transHF_vec.==7]]);
hf_ta = merge(hf_ta,hf_ta_temp)
# if 8, multiply by 100
hf_ta_temp = map((timestamp, values) -> (timestamp, (values).*100.0), dataHF_TA[varList_HF[transHF_vec.==6]])
hf_ta = merge(hf_ta,hf_ta_temp)
# if 9, take exp
hf_ta_temp = map((timestamp, values) -> (timestamp, exp.(values)), dataHF_TA[varList_HF[transHF_vec.==9]])
hf_ta = merge(hf_ta,hf_ta_temp)

hf_ta = hf_ta[varList_HF];  # reorder back the variables