include("devPkgs.jl")

using BEAVARs
using DelimitedFiles


YY_20 = readdlm("data/FRED_Chan2020_LBA.txt",',');
YY = YY_20[:,[1,2,4,8]];

setup_str, hyper_str = makeSetup(YY,"Banbura2010")
store_beta, store_sigma = Banbura2010(YY,setup_str,hyper_str);

IRF_median, IRF_68_low, IRF_68_high, IRF_4d = irf_chol_overDraws(store_beta,store_sigma,setup_str);



# using Plots
# i_var = 1;
# i_shock = 1;
# plot(IRF_median[:,i_var,i_shock],ribbon = (IRF_median[:,i_var,i_shock].-IRF_68_low[:,i_var,1],IRF_68_high[:,i_var,i_shock].-IRF_median[:,i_var,1]),fillalpha=0.2)

# n = setup_str.n
# plt = plot(layout=(3,2))
# # for i_var = 1:n
#     for i_var = 1:n;
#         plot!(plt,IRF_median[:,i_var,i_shock],ribbon = (IRF_median[:,i_var,i_shock].-IRF_68_low[:,i_var,1],IRF_68_high[:,i_var,i_shock].-IRF_median[:,i_var,1]),fillalpha=0.2,titlefontsize=8,subplot=i_var)
#     end
# # end
# display(plt)
