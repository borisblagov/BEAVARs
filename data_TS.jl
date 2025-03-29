include("devPkgs.jl")

using BEAVARs
using DelimitedFiles
using TimeSeries
using CSV
using Dates
using DataFrames
using XLSX
using Plots
using BEAVARs

#t_csv = CSV.File("data/data_tpu.csv")
data_ta_full = readtimearray("data/data_tpu.csv"; format="dd/mm/yyyy", delim=',')
#TimeArray(t_csv; timestamp=:date)

#tt = XLSX.readtable("data/data_tpu.xlsx", "Sheet1")
#df=DataFrame(tt)

data_de = data_ta_full[[:tpuCaldara, :pmiDE, :gdpDE, :invDE, :hicpDE, :euribor]];
var_names = colnames(data_de)
YY = values(data_de);

n= size(data_de,2);
p = 4;
n_irf = 16;
intercept = -1;


hyper = hypChan2020_CSV();
store_A, store_h, store_Σ, store_s2_h = Chan2020_LBA_CSV(YY;hyp=hyper,nsave = 5000,nburn = 5000);

IRF_median, IRF_68_low, IRF_68_high = irf_chol_overDraws_csv(store_A,store_Σ,store_h,n,p,intercept,n_irf;shSize = "stdev");

i_var = 4;
i_shock = 1;
plot(IRF_median[:,i_var,i_shock],ribbon = (IRF_median[:,i_var,i_shock].-IRF_68_low[:,i_var,1],IRF_68_high[:,i_var,i_shock].-IRF_median[:,i_var,1]),fillalpha=0.2,label="",title=var_names[i_var])

plot(plot(0:10),plot(0:10))

plt = plot(layout=(3,2))
# for i_var = 1:n
    for i_var = 1:n;
        plot!(plt,IRF_median[:,i_var,i_shock],ribbon = (IRF_median[:,i_var,i_shock].-IRF_68_low[:,i_var,1],IRF_68_high[:,i_var,i_shock].-IRF_median[:,i_var,1]),fillalpha=0.2,label="",title=var_names[i_var],titlefontsize=8,subplot=i_var)
    end
# end
display(plt)