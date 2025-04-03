include("devPkgs.jl")
using BEAVARs
using DelimitedFiles
using Plots
using Parameters


YY = readdlm("data/FRED_Chan2020_LBA.txt", ',');
setup_str, hyper_str = makeSetup(YY,"Chan2020_LBA_csv")
store_beta, store_h, store_Σ, s2_h_store, ρ_store, σ_h2_store, eh_store = Chan2020_LBA_csv(YY,setup_str,hyper_str);

# Using Impulse responses
IRF_median, IRF_68_low, IRF_68_high = irf_chol_overDraws_csv(store_beta,store_Σ,store_h,setup_str);

# store_beta, store_Σ =  Chan2020_LBA_Minn(YY;hyp=hypChan2020,VARSetup = setup_str);
# IRF_median, IRF_68_low, IRF_68_high = irf_chol_overDraws(store_beta,store_Σ,setup_str;shSize = "stdev");

i_var = 3;
plot(IRF_median[:,i_var,i_shock],ribbon = (IRF_median[:,i_var,i_shock].-IRF_68_low[:,i_var,1],IRF_68_high[:,i_var,i_shock].-IRF_median[:,i_var,1]),fillalpha=0.2) # ,label="",title=var_names[i_var]

@unpack n,p=setup_str
A = reshape(store_beta[:,1],n*p+1,n)

# n = setup_str.n
i_shock = 1;            # which variable you shock
# plt = plot(layout=(3,2))
#     for i_var = 1:n;
#         plot!(plt,IRF_median[:,i_var,i_shock],ribbon = (IRF_median[:,i_var,i_shock].-IRF_68_low[:,i_var,1],IRF_68_high[:,i_var,i_shock].-IRF_median[:,i_var,1]),fillalpha=0.2,label="",title=var_names[i_var],titlefontsize=8,subplot=i_var)
#     end
# # end
# display(plt)