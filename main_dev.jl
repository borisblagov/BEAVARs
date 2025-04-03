include("devPkgs.jl")
using BEAVARs
using DelimitedFiles
using Plots
using Parameters
using Statistics
using TimeSeries



data_ta_full = readtimearray("data/data_tpu.csv"; format="dd/mm/yyyy", delim=',')

data_de = data_ta_full[[:tpuCaldara, :pmiDE, :gdpDE, :invDE, :hicpDE, :euribor]];
var_names = colnames(data_de)
YY = values(data_de);

# YY20 = readdlm("data/FRED_Chan2020_LBA.txt", ',');
# YY = YY20[:,[1,4,5,6]]
setup_str, hyper_str = makeSetup(YY,"Chan2020_LBA_csv",nburn=2000,nsave=1000)
store_beta, store_h, store_Σ, s2_h_store, ρ_store, σ_h2_store, eh_store, XX = Chan2020_LBA_csv(YY,setup_str,hyper_str);

# Using Impulse responses
IRF_median, IRF_68_low, IRF_68_high = irf_chol_overDraws_csv(store_beta,store_Σ,store_h,setup_str);

# store_beta, store_Σ =  Chan2020_LBA_Minn(YY;hyp=hypChan2020,VARSetup = setup_str);
# IRF_median, IRF_68_low, IRF_68_high = irf_chol_overDraws(store_beta,store_Σ,setup_str;shSize = "stdev");

i_var = 3;
i_shock = 1;            # which variable you shock
plot(IRF_median[:,i_var,i_shock],ribbon = (IRF_median[:,i_var,i_shock].-IRF_68_low[:,i_var,1],IRF_68_high[:,i_var,i_shock].-IRF_median[:,i_var,1]),fillalpha=0.2) # ,label="",title=var_names[i_var]

@unpack n,p=setup_str
A = reshape(median(store_beta[:,:],dims=2),n*p+1,n)

global Yfit = XX*A;
global Yact = @views YY[p+1:end,:]
plot([Yfit[:,1],Yact[:,1]])

# n = setup_str.n
plt = plot(layout=(3,2))
for i_var = 1:n;
    plot!(plt,[Yfit[:,i_var] Yact[:,i_var]],subplot=i_var)
end
display(plt)

# n = setup_str.n
# plt = plot(layout=(3,2))
#     for i_var = 1:n;
#         plot!(plt,IRF_median[:,i_var,i_shock],ribbon = (IRF_median[:,i_var,i_shock].-IRF_68_low[:,i_var,1],IRF_68_high[:,i_var,i_shock].-IRF_median[:,i_var,1]),fillalpha=0.2,label="",title=var_names[i_var],titlefontsize=8,subplot=i_var)
#     end
# # end
# display(plt)