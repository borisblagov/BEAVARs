include("devPkgs.jl")

using BEAVARs
# using DelimitedFiles
using TimeSeries
using Parameters
# using CSV
# using Dates
 using Plots

#t_csv = CSV.File("data/data_tpu.csv")
data_ta_full = readtimearray("data/data_tpu.csv"; format="dd/mm/yyyy", delim=',')
#TimeArray(t_csv; timestamp=:date)

#tt = XLSX.readtable("data/data_tpu.xlsx", "Sheet1")
#df=DataFrame(tt)

data_de = data_ta_full[[:tpuCaldara, :pmiDE, :gdpDE, :invDE, :hicpDE, :euribor]];
var_names = colnames(data_de)
YY = values(data_de);



hyper = hypChan2020_csv();
store_A, store_h, store_Σ, s2_h_store, ρ_store, σ_h2_store, eh_store = Chan2020_LBA_csv(YY;hyp=hyper);
@btime Chan2020_LBA_csv(YY;hyp=hyper,nsave=100,nburn=100);

VARsetup = makeSetup(YY,nsave=1000,nburn=1000)
@unpack n,p, n_irf, const_loc = VARsetup
const_loc = 1
store_A, store_h, store_Σ, s2_h_store, ρ_store, σ_h2_store, eh_store = Chan2020_LBA_csv_strct(YY;hyp=hyper,VARSetup = VARsetup);
@btime Chan2020_LBA_csv_strct(YY;hyp=hyper,VARSetup = VARsetup);

IRF_median, IRF_68_low, IRF_68_high = irf_chol_overDraws_csv(store_A,store_Σ,store_h,n,p,const_loc,n_irf;shSize = "stdev");

i_var = 4;
i_shock = 1;
plot(IRF_median[:,i_var,i_shock],ribbon = (IRF_median[:,i_var,i_shock].-IRF_68_low[:,i_var,1],IRF_68_high[:,i_var,i_shock].-IRF_median[:,i_var,1]),fillalpha=0.2,label="",title=var_names[i_var])


plt = plot(layout=(3,2))
# for i_var = 1:n
    for i_var = 1:n;
        plot!(plt,IRF_median[:,i_var,i_shock],ribbon = (IRF_median[:,i_var,i_shock].-IRF_68_low[:,i_var,1],IRF_68_high[:,i_var,i_shock].-IRF_median[:,i_var,1]),fillalpha=0.2,label="",title=var_names[i_var],titlefontsize=8,subplot=i_var)
    end
# end
display(plt)