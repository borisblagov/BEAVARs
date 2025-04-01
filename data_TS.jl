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


VARsetup, hyp_str = makeSetup(YY,"Chan2020_LBA_csv",n_irf=32,nsave=1000,nburn=1000)
store_beta, store_h, store_Σ, s2_h_store, ρ_store, σ_h2_store, eh_store = Chan2020_LBA_csv(YY;VARSetup = VARsetup,hyp=hyp_str);
# @btime Chan2020_LBA_csv(YY;hyp=hyper,VARSetup = VARsetup);

IRF_median, IRF_68_low, IRF_68_high = irf_chol_overDraws_csv(beta_store,store_Σ,store_h,VARsetup;shSize = "stdev");

store_beta, store_Σ =  Chan2020_LBA_Minn(YY;hyp=hypChan2020,VARSetup = VARsetup)

IRF_median, IRF_68_low, IRF_68_high = irf_chol_overDraws(beta_store,store_Σ,VARsetup;shSize = "stdev");

i_var = 3;
i_shock = 1;
plot(IRF_median[:,i_var,i_shock],ribbon = (IRF_median[:,i_var,i_shock].-IRF_68_low[:,i_var,1],IRF_68_high[:,i_var,i_shock].-IRF_median[:,i_var,1]),fillalpha=0.2,label="",title=var_names[i_var])

n=6

plt = plot(layout=(3,2))
# for i_var = 1:n
    for i_var = 1:n;
        plot!(plt,IRF_median[:,i_var,i_shock],ribbon = (IRF_median[:,i_var,i_shock].-IRF_68_low[:,i_var,1],IRF_68_high[:,i_var,i_shock].-IRF_median[:,i_var,1]),fillalpha=0.2,label="",title=var_names[i_var],titlefontsize=8,subplot=i_var)
    end
# end
display(plt)