include("devPkgs.jl")
using BEAVARs
using DelimitedFiles
using Plots
using Parameters
using Statistics
using TimeSeries
using LinearAlgebra
using XLSX

data_ta_full = readtimearray("data/data_tpu.csv"; format="dd/mm/yyyy", delim=',')
data_ta_full = readtimearray("data/data_for_csv.csv"; format="dd/mm/yyyy", delim=',')

xf = XLSX.readxlsx("data/results_new.xlsx")

cumGDP_effects = zeros(16,6);
cumINV_effects = zeros(16,6);
#-------------------------------------------------#
#
#           Germany Caldara
#
#-------------------------------------------------#

# data_de = data_ta_full[1:100,[:tpuCaldara, :pmiDE, :gdpDE, :invDE, :hicpDE, :euribor]]
# var_names = colnames(data_de)
# YY = values(data_de);

# setup_str, hyper_str = makeSetup(YY,"Chan2020_LBA_Minn",p=2,nburn=5000,nsave=5000,n_fcst=200)
# store_beta, store_Σ = Chan2020_LBA_Minn(YY,setup_str,hyper_str);

# # Using Impulse responses
# IRF_median, IRF_68_low, IRF_68_high = irf_chol_overDraws(store_beta,store_Σ,setup_str,shSize=0.4);



# @unpack n_fcst,n,p,nsave=setup_str

# # n = setup_str.n
# i_shock = 1;            # which variable you shock
# plt = plot(layout=(3,2))
#     for i_var = 1:n;
#         plot!(plt,IRF_median[:,i_var,i_shock],ribbon = (IRF_median[:,i_var,i_shock].-IRF_68_low[:,i_var,1],IRF_68_high[:,i_var,i_shock].-IRF_median[:,i_var,1]),fillalpha=0.2,label="",title=var_names[i_var],titlefontsize=8,subplot=i_var)
#     end
# # end
# display(plt)


#-------------------------------------------------#
#
#           Germany Baker
#
#-------------------------------------------------#

data_deB = data_ta_full[1:80,[:tpuBaker, :pmiDE, :gdpDE, :invDE, :hicpDE, :euribor]]
var_names = colnames(data_deB)
YYdeB = values(data_deB);

setup_str, hyper_str = makeSetup(YYdeB,"Chan2020_LBA_Minn",p=2,nburn=5000,nsave=5000,n_fcst=200)
store_beta, store_Σ = Chan2020_LBA_Minn(YYdeB,setup_str,hyper_str);

# Using Impulse responses
IRF_median, IRF_68_low, IRF_68_high = irf_chol_overDraws(store_beta,store_Σ,setup_str,shSize=1.4);


@unpack n_fcst,n,p,nsave=setup_str

# -----------------------
# Full plot
# -----------------------
# i_shock = 1;            # which variable you shock
# plt = plot(layout=(3,2))
#     for i_var = 1:n;
#         plot!(plt,IRF_median[:,i_var,i_shock],ribbon = (IRF_median[:,i_var,i_shock].-IRF_68_low[:,i_var,1],IRF_68_high[:,i_var,i_shock].-IRF_median[:,i_var,1]),fillalpha=0.2,label="",title=var_names[i_var],titlefontsize=8,subplot=i_var)
#     end
# # end
# display(plt)

#-------------------------
# Paper plot
#-------------------------


plot_add = 2;
i_shock = 1;            # which variable you shock
plt = plot(layout=(6,2), size=(600,950))
    for i_var = 3:4;
        plot!(plt,IRF_median[:,i_var,i_shock],ribbon = (IRF_median[:,i_var,i_shock].-IRF_68_low[:,i_var,1],IRF_68_high[:,i_var,i_shock].-IRF_median[:,i_var,1]),fillalpha=0.2,label="",title=var_names[i_var],titlefontsize=8,subplot=i_var-2)
        if i_var == 3
        # plot!(plt,gdp_mat[:,2],ribbon = (gdp_mat[:,2].-gdp_mat[:,1],gdp_mat[:,3].-gdp_mat[:,2]),fillalpha=0.2,legend=false,subplot=i_var-plot_add)
        else
        # plot!(plt,inv_mat[:,2],ribbon = (inv_mat[:,2].-inv_mat[:,1],inv_mat[:,3].-inv_mat[:,2]),fillalpha=0.2,legend=false,subplot=i_var-plot_add)
        end
    end
# end
# display(plt)

cumGDP_effects[:,1] = cumsum(IRF_median[:,3,1]);
cumINV_effects[:,1] = cumsum(IRF_median[:,4,1]);

#-------------------------------------------------#
#
#           Germany CSV Baker
#
#-------------------------------------------------#

# data_deB = data_ta_full[1:100,[:tpuBaker, :pmiDE, :gdpDE, :invDE, :hicpDE, :euribor]]
# var_names = colnames(data_deB)
# YYdeB = values(data_deB);


# # ---- CSV model
# store_beta, store_h, store_Σ, s2_h_store, ρ_store, σ_h2_store, eh_store = Chan2020_LBA_csv(YY,setup_str,hyper_str);
# # Using Impulse responses
# IRF_median2, IRF_68_low2, IRF_68_high2 = irf_chol_overDraws_csv(store_beta,store_Σ,store_h,setup_str,shSize=1.4);

# i_shock = 1;            # which variable you shock
# #plt = plot(layout=(3,2))
#     for i_var = 1:n;
#         plot!(plt,IRF_median2[:,i_var,i_shock],ribbon = (IRF_median2[:,i_var,i_shock].-IRF_68_low2[:,i_var,1],IRF_68_high2[:,i_var,i_shock].-IRF_median2[:,i_var,1]),fillalpha=0.2,label="",title=var_names[i_var],titlefontsize=8,subplot=i_var)
#     end
# # end
# display(plt)


#-------------------------------------------------#
#
#           France Caldara
#
#-------------------------------------------------#


# data_fr = data_ta_full[1:80,[:tpuCaldara, :pmiFR, :gdpFR, :invFR, :hicpFR, :euribor]]
# var_names = colnames(data_fr)
# YYfr = values(data_fr);

# setup_str, hyper_str = makeSetup(YYfr,"Chan2020_LBA_Minn",p=2,nburn=5000,nsave=5000,n_fcst=200)
# store_beta, store_Σ = Chan2020_LBA_Minn(YYfr,setup_str,hyper_str);

# # Using Impulse responses
# IRF_median, IRF_68_low, IRF_68_high = irf_chol_overDraws(store_beta,store_Σ,setup_str,shSize=0.4);


# @unpack n_fcst,n,p,nsave=setup_str

# # n = setup_str.n
# i_shock = 1;            # which variable you shock
# plt = plot(layout=(3,2),plot_title="France",titlefontsize=10)
#     for i_var = 1:n;
#         plot!(plt,IRF_median[:,i_var,i_shock],ribbon = (IRF_median[:,i_var,i_shock].-IRF_68_low[:,i_var,1],IRF_68_high[:,i_var,i_shock].-IRF_median[:,i_var,1]),fillalpha=0.2,label="",title=var_names[i_var],titlefontsize=8,subplot=i_var)
#     end
# # end
# display(plt)


#-------------------------------------------------#
#
#           France Baker
#
#-------------------------------------------------#


data_frB = data_ta_full[1:80,[:tpuBaker, :pmiFR, :gdpFR, :invFR, :hicpFR, :euribor]]
var_names = colnames(data_frB)
YYfrB = values(data_frB);

setup_str, hyper_str = makeSetup(YYfrB,"Chan2020_LBA_Minn",p=2,nburn=5000,nsave=5000,n_fcst=200)
store_beta, store_Σ = Chan2020_LBA_Minn(YYfrB,setup_str,hyper_str);

# Using Impulse responses
IRF_median, IRF_68_low, IRF_68_high = irf_chol_overDraws(store_beta,store_Σ,setup_str,shSize=1.4);


@unpack n_fcst,n,p,nsave=setup_str

# n = setup_str.n
# i_shock = 1;            # which variable you shock
# # plt = plot(layout=(3,2))
#     for i_var = 1:n;
#         plot!(plt,IRF_median[:,i_var,i_shock],ribbon = (IRF_median[:,i_var,i_shock].-IRF_68_low[:,i_var,1],IRF_68_high[:,i_var,i_shock].-IRF_median[:,i_var,1]),fillalpha=0.2,label="",title=var_names[i_var],titlefontsize=8,subplot=i_var)
#     end
# # end
# display(plt)

#-------------------------
# Paper plot
#-------------------------


gdp_sh = xf["gdpFR"]
inv_sh = xf["invFR"]

gdp_mat = convert(Array{Float64},gdp_sh["B2:D17"])
inv_mat = convert(Array{Float64},inv_sh["B2:D17"])

plot_add = 0; # adds a number to specify whihc plot we need
i_shock = 1;            # which variable you shock
# plt = plot(layout=(6,2), size=(600,950))
    for i_var = 3:4;
        plot!(plt,IRF_median[:,i_var,i_shock],ribbon = (IRF_median[:,i_var,i_shock].-IRF_68_low[:,i_var,1],IRF_68_high[:,i_var,i_shock].-IRF_median[:,i_var,1]),fillalpha=0.2,label="",title=var_names[i_var],titlefontsize=8,subplot=i_var-plot_add)
        if i_var == 3
        # plot!(plt,gdp_mat[:,2],ribbon = (gdp_mat[:,2].-gdp_mat[:,1],gdp_mat[:,3].-gdp_mat[:,2]),fillalpha=0.2,legend=false,subplot=i_var-plot_add)
        else
        # plot!(plt,inv_mat[:,2],ribbon = (inv_mat[:,2].-inv_mat[:,1],inv_mat[:,3].-inv_mat[:,2]),fillalpha=0.2,legend=false,subplot=i_var-plot_add)
        end
    end
# end
# display(plt)

cumGDP_effects[:,2] = cumsum(IRF_median[:,3,1]);
cumINV_effects[:,2] = cumsum(IRF_median[:,4,1]);
# for i_var = 3:4;
#     # plot!(plt,IRF_median[:,i_var,i_shock],ribbon = (IRF_median[:,i_var,i_shock].-IRF_68_low[:,i_var,1],IRF_68_high[:,i_var,i_shock].-IRF_median[:,i_var,1]),fillalpha=0.2,label="",title=var_names[i_var],titlefontsize=8,subplot=i_var-2)
#     if i_var == 3
#     plot!(plt_sv,gdp_mat[:,2],ribbon = (gdp_mat[:,2].-gdp_mat[:,1],gdp_mat[:,3].-gdp_mat[:,2]),fillalpha=0.2,legend=false,subplot=i_var-plot_add,title=var_names[i_var],titlefontsize=8)
#     else
#     plot!(plt_sv,inv_mat[:,2],ribbon = (inv_mat[:,2].-inv_mat[:,1],inv_mat[:,3].-inv_mat[:,2]),fillalpha=0.2,legend=false,subplot=i_var-plot_add,title=var_names[i_var],titlefontsize=8)
#     end
# end
# display(plt_sv)


#-------------------------------------------------#
#
#           Italy Caldara
#
#-------------------------------------------------#


# data_it = data_ta_full[1:80,[:tpuCaldara, :pmiIT, :gdpIT, :invIT, :hicpIT, :euribor]]
# var_names = colnames(data_it)
# YYit = values(data_it);

# setup_str, hyper_str = makeSetup(YYit,"Chan2020_LBA_Minn",p=2,nburn=5000,nsave=5000,n_fcst=200)
# store_beta, store_Σ = Chan2020_LBA_Minn(YYit,setup_str,hyper_str);

# # Using Impulse responses
# IRF_median, IRF_68_low, IRF_68_high = irf_chol_overDraws(store_beta,store_Σ,setup_str,shSize=0.4);


# @unpack n_fcst,n,p,nsave=setup_str

# # n = setup_str.n
# i_shock = 1;            # which variable you shock
# plt = plot(layout=(3,2))
#     for i_var = 1:n;
#         plot!(plt,IRF_median[:,i_var,i_shock],ribbon = (IRF_median[:,i_var,i_shock].-IRF_68_low[:,i_var,1],IRF_68_high[:,i_var,i_shock].-IRF_median[:,i_var,1]),fillalpha=0.2,label="",title=var_names[i_var],titlefontsize=8,subplot=i_var)
#     end
# # end
# display(plt)



#-------------------------------------------------#
#
#           Italy Baker
#
#-------------------------------------------------#


data_itB = data_ta_full[1:80,[:tpuBaker, :pmiIT, :gdpIT, :invIT, :hicpIT, :euribor]]
var_names = colnames(data_itB)
YYitB = values(data_itB);

setup_str, hyper_str = makeSetup(YYitB,"Chan2020_LBA_Minn",p=2,nburn=5000,nsave=5000,n_fcst=200)
store_beta, store_Σ = Chan2020_LBA_Minn(YYitB,setup_str,hyper_str);

# Using Impulse responses
IRF_median, IRF_68_low, IRF_68_high = irf_chol_overDraws(store_beta,store_Σ,setup_str,shSize=1.4);


@unpack n_fcst,n,p,nsave=setup_str

# n = setup_str.n
# i_shock = 1;            # which variable you shock
# # plt = plot(layout=(3,2))
#     for i_var = 1:n;
#         plot!(plt,IRF_median[:,i_var,i_shock],ribbon = (IRF_median[:,i_var,i_shock].-IRF_68_low[:,i_var,1],IRF_68_high[:,i_var,i_shock].-IRF_median[:,i_var,1]),fillalpha=0.2,label="",title=var_names[i_var],titlefontsize=8,subplot=i_var)
#     end
# # end
# display(plt)


cumGDP_effects[:,3] = cumsum(IRF_median[:,3,1]);
cumINV_effects[:,3] = cumsum(IRF_median[:,4,1]);
#-------------------------
# Paper plot
#-------------------------


gdp_sh = xf["gdpIT"]
inv_sh = xf["invIT"]

gdp_mat = convert(Array{Float64},gdp_sh["B2:D17"])
inv_mat = convert(Array{Float64},inv_sh["B2:D17"])

plot_add = -2; # adds a number to specify whihc plot we need
i_shock = 1;            # which variable you shock
# plt = plot(layout=(6,2), size=(600,950))
    for i_var = 3:4;
        plot!(plt,IRF_median[:,i_var,i_shock],ribbon = (IRF_median[:,i_var,i_shock].-IRF_68_low[:,i_var,1],IRF_68_high[:,i_var,i_shock].-IRF_median[:,i_var,1]),fillalpha=0.2,label="",title=var_names[i_var],titlefontsize=8,subplot=i_var-plot_add)
        if i_var == 3
        # plot!(plt,gdp_mat[:,2],ribbon = (gdp_mat[:,2].-gdp_mat[:,1],gdp_mat[:,3].-gdp_mat[:,2]),fillalpha=0.2,legend=false,subplot=i_var-plot_add)
        else
        # plot!(plt,inv_mat[:,2],ribbon = (inv_mat[:,2].-inv_mat[:,1],inv_mat[:,3].-inv_mat[:,2]),fillalpha=0.2,legend=false,subplot=i_var-plot_add)
        end
    end
# end
# display(plt)




#-------------------------------------------------#
#
#           Spain Caldara
#
#-------------------------------------------------#


# data_es = data_ta_full[1:80,[:tpuCaldara, :pmiES, :gdpES, :invES, :hicpES, :euribor]]
# var_names = colnames(data_es)
# YYes = values(data_es);

# setup_str, hyper_str = makeSetup(YYes,"Chan2020_LBA_Minn",p=2,nburn=5000,nsave=5000,n_fcst=200)
# store_beta, store_Σ = Chan2020_LBA_Minn(YYes,setup_str,hyper_str);

# # Using Impulse responses
# IRF_median, IRF_68_low, IRF_68_high = irf_chol_overDraws(store_beta,store_Σ,setup_str,shSize=0.4);


# @unpack n_fcst,n,p,nsave=setup_str

# # n = setup_str.n
# i_shock = 1;            # which variable you shock
# plt = plot(layout=(3,2))
#     for i_var = 1:n;
#         plot!(plt,IRF_median[:,i_var,i_shock],ribbon = (IRF_median[:,i_var,i_shock].-IRF_68_low[:,i_var,1],IRF_68_high[:,i_var,i_shock].-IRF_median[:,i_var,1]),fillalpha=0.2,label="",title=var_names[i_var],titlefontsize=8,subplot=i_var)
#     end
# # end
# display(plt)



#-------------------------------------------------#
#
#           Spain Baker
#
#-------------------------------------------------#


data_esB = data_ta_full[1:80,[:tpuBaker, :pmiES, :gdpES, :invES, :hicpES, :euribor]]
var_names = colnames(data_esB)
YYesB = values(data_esB);

setup_str, hyper_str = makeSetup(YYesB,"Chan2020_LBA_Minn",p=2,nburn=25000,nsave=5000,n_fcst=200)
store_beta, store_Σ = Chan2020_LBA_Minn(YYesB,setup_str,hyper_str);

# Using Impulse responses
IRF_median, IRF_68_low, IRF_68_high = irf_chol_overDraws(store_beta,store_Σ,setup_str,shSize=1.4);


@unpack n_fcst,n,p,nsave=setup_str

# # n = setup_str.n
# i_shock = 1;            # which variable you shock
# #plt = plot(layout=(3,2))
#     for i_var = 1:n;
#         plot!(plt,IRF_median[:,i_var,i_shock],ribbon = (IRF_median[:,i_var,i_shock].-IRF_68_low[:,i_var,1],IRF_68_high[:,i_var,i_shock].-IRF_median[:,i_var,1]),fillalpha=0.2,label="",title=var_names[i_var],titlefontsize=8,subplot=i_var)
#     end
# # end
# display(plt)

cumGDP_effects[:,4] = cumsum(IRF_median[:,3,1]);
cumINV_effects[:,4] = cumsum(IRF_median[:,4,1]);

#-------------------------
# Paper plot
#-------------------------


gdp_sh = xf["gdpES"]
inv_sh = xf["invES"]

gdp_mat = convert(Array{Float64},gdp_sh["B2:D17"])
inv_mat = convert(Array{Float64},inv_sh["B2:D17"])

plot_add = -4; # adds a number to specify whihc plot we need
i_shock = 1;            # which variable you shock
# plt = plot(layout=(6,2), size=(600,950))
    for i_var = 3:4;
        plot!(plt,IRF_median[:,i_var,i_shock],ribbon = (IRF_median[:,i_var,i_shock].-IRF_68_low[:,i_var,1],IRF_68_high[:,i_var,i_shock].-IRF_median[:,i_var,1]),fillalpha=0.2,label="",title=var_names[i_var],titlefontsize=8,subplot=i_var-plot_add)
        if i_var == 3
        # plot!(plt,gdp_mat[:,2],ribbon = (gdp_mat[:,2].-gdp_mat[:,1],gdp_mat[:,3].-gdp_mat[:,2]),fillalpha=0.2,legend=false,subplot=i_var-plot_add)
        else
        # plot!(plt,inv_mat[:,2],ribbon = (inv_mat[:,2].-inv_mat[:,1],inv_mat[:,3].-inv_mat[:,2]),fillalpha=0.2,legend=false,subplot=i_var-plot_add)
        end
    end
# end
# display(plt)


#-------------------------------------------------#
#
#           Spain CSV Baker
#
#-------------------------------------------------#

# data_esB = data_ta_full[1:100,[:tpuBaker, :pmiES, :gdpES, :invES, :hicpES, :euribor]]
# var_names = colnames(data_esB)
# YY = values(data_esB);

# # ---- CSV model
# store_beta, store_h, store_Σ, s2_h_store, ρ_store, σ_h2_store, eh_store = Chan2020_LBA_csv(YY,setup_str,hyper_str);
# # Using Impulse responses
# IRF_median2, IRF_68_low2, IRF_68_high2 = irf_chol_overDraws_csv(store_beta,store_Σ,store_h,setup_str,shSize=1.4);

# i_shock = 1;            # which variable you shock
# #plt = plot(layout=(3,2))
#     for i_var = 1:n;
#         plot!(plt,IRF_median2[:,i_var,i_shock],ribbon = (IRF_median2[:,i_var,i_shock].-IRF_68_low2[:,i_var,1],IRF_68_high2[:,i_var,i_shock].-IRF_median2[:,i_var,1]),fillalpha=0.2,label="",title=var_names[i_var],titlefontsize=8,subplot=i_var)
#     end
# # end
# display(plt)


#-------------------------------------------------#
#
#           EA Caldara
#
#-------------------------------------------------#


# data_ea = data_ta_full[1:80,[:tpuCaldara, :pmiEA, :gdpEA, :invEA, :hicpEA, :euribor]]
# var_names = colnames(data_ea)
# YY = values(data_ea);

# setup_str, hyper_str = makeSetup(YY,"Chan2020_LBA_Minn",p=2,nburn=5000,nsave=5000,n_fcst=200)
# store_beta, store_Σ = Chan2020_LBA_Minn(YY,setup_str,hyper_str);

# # Using Impulse responses
# IRF_median, IRF_68_low, IRF_68_high = irf_chol_overDraws(store_beta,store_Σ,setup_str,shSize=0.4);


# @unpack n_fcst,n,p,nsave=setup_str

# # n = setup_str.n
# i_shock = 1;            # which variable you shock
# plt = plot(layout=(3,2))
#     for i_var = 1:n;
#         plot!(plt,IRF_median[:,i_var,i_shock],ribbon = (IRF_median[:,i_var,i_shock].-IRF_68_low[:,i_var,1],IRF_68_high[:,i_var,i_shock].-IRF_median[:,i_var,1]),fillalpha=0.2,label="",title=var_names[i_var],titlefontsize=8,subplot=i_var)
#     end
# # end
# display(plt)



#-------------------------------------------------#
#
#           EA Baker
#
#-------------------------------------------------#


# data_eaB = data_ta_full[1:80,[:tpuBaker, :pmiEA, :gdpEA, :invEA, :hicpEA, :euribor]]
# var_names = colnames(data_eaB)
# YY = values(data_eaB);

# setup_str, hyper_str = makeSetup(YY,"Chan2020_LBA_Minn",nburn=5000,nsave=5000,n_fcst=200)
# store_beta, store_Σ = Chan2020_LBA_Minn(YY,setup_str,hyper_str);

# # Using Impulse responses
# IRF_median, IRF_68_low, IRF_68_high = irf_chol_overDraws(store_beta,store_Σ,setup_str,shSize=1.4);


# @unpack n_fcst,n,p,nsave=setup_str

# # n = setup_str.n
# i_shock = 1;            # which variable you shock
# # plt = plot(layout=(3,2))
#     for i_var = 1:n;
#         plot!(plt,IRF_median[:,i_var,i_shock],ribbon = (IRF_median[:,i_var,i_shock].-IRF_68_low[:,i_var,1],IRF_68_high[:,i_var,i_shock].-IRF_median[:,i_var,1]),fillalpha=0.2,label="",title=var_names[i_var],titlefontsize=8,subplot=i_var)
#     end
# # end
# display(plt)




#-------------------------------------------------#
#
#           netherlands Caldara
#
#-------------------------------------------------#

# data_nl = data_ta_full[1:80,[:tpuCaldara, :pmiNL, :gdpNL, :invNL, :hicpNL, :euribor]]
# var_names = colnames(data_nl)
# YY = values(data_nl);

# setup_str, hyper_str = makeSetup(YY,"Chan2020_LBA_Minn",nburn=5000,nsave=5000,n_fcst=200)
# store_beta, store_Σ = Chan2020_LBA_Minn(YY,setup_str,hyper_str);

# # Using Impulse responses
# IRF_median, IRF_68_low, IRF_68_high = irf_chol_overDraws(store_beta,store_Σ,setup_str,shSize=0.4);


# @unpack n_fcst,n,p,nsave=setup_str

# # n = setup_str.n
# i_shock = 1;            # which variable you shock
# plt = plot(layout=(3,2))
#     for i_var = 1:n;
#         plot!(plt,IRF_median[:,i_var,i_shock],ribbon = (IRF_median[:,i_var,i_shock].-IRF_68_low[:,i_var,1],IRF_68_high[:,i_var,i_shock].-IRF_median[:,i_var,1]),fillalpha=0.2,label="",title=var_names[i_var],titlefontsize=8,subplot=i_var)
#     end
# # end
# display(plt)


#-------------------------------------------------#
#
#           netherlands Baker
#
#-------------------------------------------------#

data_nl = data_ta_full[1:80,[:tpuBaker, :pmiNL, :gdpNL, :invNL, :hicpNL, :euribor]]
var_names = colnames(data_nl)
YY = values(data_nl);

setup_str, hyper_str = makeSetup(YY,"Chan2020_LBA_Minn",nburn=5000,nsave=5000,n_fcst=200)
store_beta, store_Σ = Chan2020_LBA_Minn(YY,setup_str,hyper_str);

# Using Impulse responses
IRF_median, IRF_68_low, IRF_68_high = irf_chol_overDraws(store_beta,store_Σ,setup_str,shSize=1.4);


@unpack n_fcst,n,p,nsave=setup_str

# i_shock = 1;            # which variable you shock
# # plt = plot(layout=(3,2))
#     for i_var = 1:n;
#         plot!(plt,IRF_median[:,i_var,i_shock],ribbon = (IRF_median[:,i_var,i_shock].-IRF_68_low[:,i_var,1],IRF_68_high[:,i_var,i_shock].-IRF_median[:,i_var,1]),fillalpha=0.2,label="",title=var_names[i_var],titlefontsize=8,subplot=i_var)
#     end
# # end
# display(plt)

cumGDP_effects[:,5] = cumsum(IRF_median[:,3,1]);
cumINV_effects[:,5] = cumsum(IRF_median[:,4,1]);

#-------------------------
# Paper plot
#-------------------------


gdp_sh = xf["gdpNL"]
inv_sh = xf["invNL"]

gdp_mat = convert(Array{Float64},gdp_sh["B2:D17"])
inv_mat = convert(Array{Float64},inv_sh["B2:D17"])

plot_add = -6; # adds a number to specify whihc plot we need
i_shock = 1;            # which variable you shock
# plt = plot(layout=(6,2), size=(600,950))
    for i_var = 3:4;
        plot!(plt,IRF_median[:,i_var,i_shock],ribbon = (IRF_median[:,i_var,i_shock].-IRF_68_low[:,i_var,1],IRF_68_high[:,i_var,i_shock].-IRF_median[:,i_var,1]),fillalpha=0.2,label="",title=var_names[i_var],titlefontsize=8,subplot=i_var-plot_add)
        if i_var == 3
        # plot!(plt,gdp_mat[:,2],ribbon = (gdp_mat[:,2].-gdp_mat[:,1],gdp_mat[:,3].-gdp_mat[:,2]),fillalpha=0.2,legend=false,subplot=i_var-plot_add)
        else
        # plot!(plt,inv_mat[:,2],ribbon = (inv_mat[:,2].-inv_mat[:,1],inv_mat[:,3].-inv_mat[:,2]),fillalpha=0.2,legend=false,subplot=i_var-plot_add)
        end
    end
# end
# display(plt)




#-------------------------------------------------#
#
#           United Kingdom Caldara
#
#-------------------------------------------------#

# data_uk = data_ta_full[1:80,[:tpuCaldara, :pmiUK, :gdpUK, :invUK, :cpiUK, :polrate]]
# var_names = colnames(data_uk)
# YY = values(data_uk);

# setup_str, hyper_str = makeSetup(YY,"Chan2020_LBA_Minn",p=2,nburn=5000,nsave=5000,n_fcst=200)
# store_beta, store_Σ = Chan2020_LBA_Minn(YY,setup_str,hyper_str);

# # Using Impulse responses
# IRF_median, IRF_68_low, IRF_68_high = irf_chol_overDraws(store_beta,store_Σ,setup_str,shSize=0.4);


# @unpack n_fcst,n,p,nsave=setup_str

# # n = setup_str.n
# i_shock = 1;            # which variable you shock
# plt = plot(layout=(3,2))
#     for i_var = 1:n;
#         plot!(plt,IRF_median[:,i_var,i_shock],ribbon = (IRF_median[:,i_var,i_shock].-IRF_68_low[:,i_var,1],IRF_68_high[:,i_var,i_shock].-IRF_median[:,i_var,1]),fillalpha=0.2,label="",title=var_names[i_var],titlefontsize=8,subplot=i_var)
#     end
# # end
# display(plt)

#-------------------------------------------------#
#
#           United Kingdom Baker
#
#-------------------------------------------------#

data_uk = data_ta_full[1:80,[:tpuBaker, :pmiUK, :gdpUK, :invUK, :cpiUK, :polrate]]
var_names = colnames(data_uk)
YY = values(data_uk);

setup_str, hyper_str = makeSetup(YY,"Chan2020_LBA_Minn",p=2,nburn=5000,nsave=5000,n_fcst=200)
store_beta, store_Σ = Chan2020_LBA_Minn(YY,setup_str,hyper_str);

# Using Impulse responses
IRF_median, IRF_68_low, IRF_68_high = irf_chol_overDraws(store_beta,store_Σ,setup_str,shSize=1.4);


cumGDP_effects[:,6] = cumsum(IRF_median[:,3,1]);
cumINV_effects[:,6] = cumsum(IRF_median[:,4,1]);

@unpack n_fcst,n,p,nsave=setup_str

# i_shock = 1;            # which variable you shock
# # plt = plot(layout=(3,2))
#     for i_var = 1:n;
#         plot!(plt,IRF_median[:,i_var,i_shock],ribbon = (IRF_median[:,i_var,i_shock].-IRF_68_low[:,i_var,1],IRF_68_high[:,i_var,i_shock].-IRF_median[:,i_var,1]),fillalpha=0.2,label="",title=var_names[i_var],titlefontsize=8,subplot=i_var)
#     end
# # end
# display(plt)



gdp_sh = xf["gdpUK"]
inv_sh = xf["invUK"]

gdp_mat = convert(Array{Float64},gdp_sh["B2:D17"])
inv_mat = convert(Array{Float64},inv_sh["B2:D17"])

plot_add = -8; # adds a number to specify whihc plot we need
i_shock = 1;            # which variable you shock
# plt = plot(layout=(6,2), size=(600,950))
    for i_var = 3:4;
        plot!(plt,IRF_median[:,i_var,i_shock],ribbon = (IRF_median[:,i_var,i_shock].-IRF_68_low[:,i_var,1],IRF_68_high[:,i_var,i_shock].-IRF_median[:,i_var,1]),fillalpha=0.2,label="",title=var_names[i_var],titlefontsize=8,subplot=i_var-plot_add)
        if i_var == 3
        # plot!(plt,gdp_mat[:,2],ribbon = (gdp_mat[:,2].-gdp_mat[:,1],gdp_mat[:,3].-gdp_mat[:,2]),fillalpha=0.2,legend=false,subplot=i_var-plot_add)
        else
        # plot!(plt,inv_mat[:,2],ribbon = (inv_mat[:,2].-inv_mat[:,1],inv_mat[:,3].-inv_mat[:,2]),fillalpha=0.2,legend=false,subplot=i_var-plot_add)
        end
    end
# end
display(plt)


#-------------------------
# plot of the cumulative effects
labs = ["Germany" "France" "Italy" "Spain" "Netherlands" "UK"]
p1 = plot(cumGDP_effects,title = "GDP",titlefontsize=10, labels=labs,legend=false)
p2 = plot(cumINV_effects,title = "Investment",titlefontsize=10, labels=labs,legend=false)
legendh = plot([0 0 0 0 0 0], label=labs, legendcolumns=3,legendfontsize=8,showaxis = false, grid = false)
plt_cum = plot(p1,p2,legendh, size=(600,400),layout = @layout([A B; C]))



#-------------------------
# SV plot
#-------------------------


gdp_cunt = "gdpDE"; inv_cunt = "invDE";
gdp_sh = xf[gdp_cunt]
inv_sh = xf[inv_cunt]

gdp_mat = convert(Array{Float64},gdp_sh["B2:D17"])
inv_mat = convert(Array{Float64},inv_sh["B2:D17"])

plot_add = 2;
plt_sv = plot(layout=(6,2), size=(600,900))
plot!(plt_sv,gdp_mat[:,2],ribbon = (gdp_mat[:,2].-gdp_mat[:,1],gdp_mat[:,3].-gdp_mat[:,2]),fillalpha=0.2,legend=false,subplot=3-plot_add,title=gdp_cunt,titlefontsize=8)
plot!(plt_sv,inv_mat[:,2],ribbon = (inv_mat[:,2].-inv_mat[:,1],inv_mat[:,3].-inv_mat[:,2]),fillalpha=0.2,legend=false,subplot=4-plot_add,title=inv_cunt,titlefontsize=8)



gdp_cunt = "gdpFR"; inv_cunt = "invFR";
plot_add = 0;

gdp_sh = xf[gdp_cunt]
inv_sh = xf[inv_cunt]

gdp_mat = convert(Array{Float64},gdp_sh["B2:D17"])
inv_mat = convert(Array{Float64},inv_sh["B2:D17"])

plot!(plt_sv,gdp_mat[:,2],ribbon = (gdp_mat[:,2].-gdp_mat[:,1],gdp_mat[:,3].-gdp_mat[:,2]),fillalpha=0.2,legend=false,subplot=3-plot_add,title=gdp_cunt,titlefontsize=8)
plot!(plt_sv,inv_mat[:,2],ribbon = (inv_mat[:,2].-inv_mat[:,1],inv_mat[:,3].-inv_mat[:,2]),fillalpha=0.2,legend=false,subplot=4-plot_add,title=inv_cunt,titlefontsize=8)



gdp_cunt = "gdpIT"; inv_cunt = "invIT";
plot_add = -2;

gdp_sh = xf[gdp_cunt]
inv_sh = xf[inv_cunt]

gdp_mat = convert(Array{Float64},gdp_sh["B2:D17"])
inv_mat = convert(Array{Float64},inv_sh["B2:D17"])

plot!(plt_sv,gdp_mat[:,2],ribbon = (gdp_mat[:,2].-gdp_mat[:,1],gdp_mat[:,3].-gdp_mat[:,2]),fillalpha=0.2,legend=false,subplot=3-plot_add,title=gdp_cunt,titlefontsize=8)
plot!(plt_sv,inv_mat[:,2],ribbon = (inv_mat[:,2].-inv_mat[:,1],inv_mat[:,3].-inv_mat[:,2]),fillalpha=0.2,legend=false,subplot=4-plot_add,title=inv_cunt,titlefontsize=8)


gdp_cunt = "gdpES"; inv_cunt = "invES";
plot_add = -4;

gdp_sh = xf[gdp_cunt]
inv_sh = xf[inv_cunt]

gdp_mat = convert(Array{Float64},gdp_sh["B2:D17"])
inv_mat = convert(Array{Float64},inv_sh["B2:D17"])

plot!(plt_sv,gdp_mat[:,2],ribbon = (gdp_mat[:,2].-gdp_mat[:,1],gdp_mat[:,3].-gdp_mat[:,2]),fillalpha=0.2,legend=false,subplot=3-plot_add,title=gdp_cunt,titlefontsize=8)
plot!(plt_sv,inv_mat[:,2],ribbon = (inv_mat[:,2].-inv_mat[:,1],inv_mat[:,3].-inv_mat[:,2]),fillalpha=0.2,legend=false,subplot=4-plot_add,title=inv_cunt,titlefontsize=8)


gdp_cunt = "gdpNL"; inv_cunt = "invNL";
plot_add = -6;

gdp_sh = xf[gdp_cunt]
inv_sh = xf[inv_cunt]

gdp_mat = convert(Array{Float64},gdp_sh["B2:D17"])
inv_mat = convert(Array{Float64},inv_sh["B2:D17"])

plot!(plt_sv,gdp_mat[:,2],ribbon = (gdp_mat[:,2].-gdp_mat[:,1],gdp_mat[:,3].-gdp_mat[:,2]),fillalpha=0.2,legend=false,subplot=3-plot_add,title=gdp_cunt,titlefontsize=8)
plot!(plt_sv,inv_mat[:,2],ribbon = (inv_mat[:,2].-inv_mat[:,1],inv_mat[:,3].-inv_mat[:,2]),fillalpha=0.2,legend=false,subplot=4-plot_add,title=inv_cunt,titlefontsize=8)


gdp_cunt = "gdpUK"; inv_cunt = "invUK";
plot_add = -8;

gdp_sh = xf[gdp_cunt]
inv_sh = xf[inv_cunt]

gdp_mat = convert(Array{Float64},gdp_sh["B2:D17"])
inv_mat = convert(Array{Float64},inv_sh["B2:D17"])

plot!(plt_sv,gdp_mat[:,2],ribbon = (gdp_mat[:,2].-gdp_mat[:,1],gdp_mat[:,3].-gdp_mat[:,2]),fillalpha=0.2,legend=false,subplot=3-plot_add,title=gdp_cunt,titlefontsize=8)
plot!(plt_sv,inv_mat[:,2],ribbon = (inv_mat[:,2].-inv_mat[:,1],inv_mat[:,3].-inv_mat[:,2]),fillalpha=0.2,legend=false,subplot=4-plot_add,title=inv_cunt,titlefontsize=8)



display(plt_sv)
