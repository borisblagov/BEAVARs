include("devPkgs.jl")
using BEAVARs
using DelimitedFiles
using Plots
using Parameters
using Statistics
using TimeSeries
using LinearAlgebra



data_ta_full = readtimearray("data/data_tpu.csv"; format="dd/mm/yyyy", delim=',')

#-------------------------------------------------#
#
#           Germany Caldara
#
#-------------------------------------------------#


data_de = data_ta_full[1:80,[:tpuCaldara, :pmiDE, :gdpDE, :invDE, :hicpDE, :euribor]]
var_names = colnames(data_de)
YY = values(data_de);


setup_str, hyper_str = makeSetup(YY,"Chan2020_LBA_Minn",nburn=5000,nsave=5000,n_fcst=200)
store_beta, store_Σ = Chan2020_LBA_Minn(YY,setup_str,hyper_str);

# Using Impulse responses
IRF_median, IRF_68_low, IRF_68_high = irf_chol_overDraws(store_beta,store_Σ,setup_str,shSize=0.4);



@unpack n_fcst,n,p,nsave=setup_str

# n = setup_str.n
i_shock = 1;            # which variable you shock
plt = plot(layout=(3,2))
    for i_var = 1:n;
        plot!(plt,IRF_median[:,i_var,i_shock],ribbon = (IRF_median[:,i_var,i_shock].-IRF_68_low[:,i_var,1],IRF_68_high[:,i_var,i_shock].-IRF_median[:,i_var,1]),fillalpha=0.2,label="",title=var_names[i_var],titlefontsize=8,subplot=i_var)
    end
# end
display(plt)


# ---- CSV model
store_beta, store_h, store_Σ, s2_h_store, ρ_store, σ_h2_store, eh_store = Chan2020_LBA_csv(YY,setup_str,hyper_str);
# Using Impulse responses
IRF_median2, IRF_68_low2, IRF_68_high2 = irf_chol_overDraws_csv(store_beta,store_Σ,store_h,setup_str,shSize=0.4);

i_shock = 1;            # which variable you shock
#plt = plot(layout=(3,2))
    for i_var = 1:n;
        plot!(plt,IRF_median2[:,i_var,i_shock],ribbon = (IRF_median2[:,i_var,i_shock].-IRF_68_low2[:,i_var,1],IRF_68_high2[:,i_var,i_shock].-IRF_median2[:,i_var,1]),fillalpha=0.2,label="",title=var_names[i_var],titlefontsize=8,subplot=i_var)
    end
# end
display(plt)


#-------------------------------------------------#
#
#           Germany Baker
#
#-------------------------------------------------#

data_deB = data_ta_full[1:80,[:tpuBaker, :pmiDE, :gdpDE, :invDE, :hicpDE, :euribor]]
var_names = colnames(data_deB)
YYdeB = values(data_deB);

setup_str, hyper_str = makeSetup(YYdeB,"Chan2020_LBA_Minn",nburn=5000,nsave=5000,n_fcst=200)
store_beta, store_Σ = Chan2020_LBA_Minn(YYdeB,setup_str,hyper_str);

# Using Impulse responses
IRF_median, IRF_68_low, IRF_68_high = irf_chol_overDraws(store_beta,store_Σ,setup_str,shSize=1.1);


@unpack n_fcst,n,p,nsave=setup_str

# n = setup_str.n
i_shock = 1;            # which variable you shock
plt = plot(layout=(3,2))
    for i_var = 1:n;
        plot!(plt,IRF_median[:,i_var,i_shock],ribbon = (IRF_median[:,i_var,i_shock].-IRF_68_low[:,i_var,1],IRF_68_high[:,i_var,i_shock].-IRF_median[:,i_var,1]),fillalpha=0.2,label="",title=var_names[i_var],titlefontsize=8,subplot=i_var)
    end
# end
display(plt)

#-------------------------------------------------#
#
#           France Caldara
#
#-------------------------------------------------#


data_fr = data_ta_full[1:80,[:tpuCaldara, :pmiFR, :gdpFR, :invFR, :hicpFR, :euribor]]
var_names = colnames(data_fr)
YYfr = values(data_fr);

setup_str, hyper_str = makeSetup(YYfr,"Chan2020_LBA_Minn",nburn=5000,nsave=5000,n_fcst=200)
store_beta, store_Σ = Chan2020_LBA_Minn(YYfr,setup_str,hyper_str);

# Using Impulse responses
IRF_median, IRF_68_low, IRF_68_high = irf_chol_overDraws(store_beta,store_Σ,setup_str,shSize=0.4);


@unpack n_fcst,n,p,nsave=setup_str

# n = setup_str.n
i_shock = 1;            # which variable you shock
plt = plot(layout=(3,2),plot_title="France",titlefontsize=10)
    for i_var = 1:n;
        plot!(plt,IRF_median[:,i_var,i_shock],ribbon = (IRF_median[:,i_var,i_shock].-IRF_68_low[:,i_var,1],IRF_68_high[:,i_var,i_shock].-IRF_median[:,i_var,1]),fillalpha=0.2,label="",title=var_names[i_var],titlefontsize=8,subplot=i_var)
    end
# end
display(plt)


# ---- CSV model
store_beta, store_h, store_Σ, s2_h_store, ρ_store, σ_h2_store, eh_store = Chan2020_LBA_csv(YYfr,setup_str,hyper_str);
# Using Impulse responses
IRF_median2, IRF_68_low2, IRF_68_high2 = irf_chol_overDraws_csv(store_beta,store_Σ,store_h,setup_str,shSize=0.4);

i_shock = 1;            # which variable you shock
plt = plot(layout=(3,2))
    for i_var = 1:n;
        plot!(plt,IRF_median2[:,i_var,i_shock],ribbon = (IRF_median2[:,i_var,i_shock].-IRF_68_low2[:,i_var,1],IRF_68_high2[:,i_var,i_shock].-IRF_median2[:,i_var,1]),fillalpha=0.2,label="",title=var_names[i_var],titlefontsize=8,subplot=i_var)
    end
# end
display(plt)

#-------------------------------------------------#
#
#           France Baker
#
#-------------------------------------------------#


data_frB = data_ta_full[1:80,[:tpuBaker, :pmiFR, :gdpFR, :invFR, :hicpFR, :euribor]]
var_names = colnames(data_frB)
YYfrB = values(data_frB);

setup_str, hyper_str = makeSetup(YYfrB,"Chan2020_LBA_Minn",nburn=5000,nsave=5000,n_fcst=200)
store_beta, store_Σ = Chan2020_LBA_Minn(YYfrB,setup_str,hyper_str);

# Using Impulse responses
IRF_median, IRF_68_low, IRF_68_high = irf_chol_overDraws(store_beta,store_Σ,setup_str,shSize=1.1);


@unpack n_fcst,n,p,nsave=setup_str

# n = setup_str.n
i_shock = 1;            # which variable you shock
plt = plot(layout=(3,2))
    for i_var = 1:n;
        plot!(plt,IRF_median[:,i_var,i_shock],ribbon = (IRF_median[:,i_var,i_shock].-IRF_68_low[:,i_var,1],IRF_68_high[:,i_var,i_shock].-IRF_median[:,i_var,1]),fillalpha=0.2,label="",title=var_names[i_var],titlefontsize=8,subplot=i_var)
    end
# end
display(plt)


#-------------------------------------------------#
#
#           Italy Caldara
#
#-------------------------------------------------#


data_it = data_ta_full[1:80,[:tpuCaldara, :pmiIT, :gdpIT, :invIT, :hicpIT, :euribor]]
var_names = colnames(data_it)
YYit = values(data_it);

setup_str, hyper_str = makeSetup(YYit,"Chan2020_LBA_Minn",nburn=5000,nsave=5000,n_fcst=200)
store_beta, store_Σ = Chan2020_LBA_Minn(YYit,setup_str,hyper_str);

# Using Impulse responses
IRF_median, IRF_68_low, IRF_68_high = irf_chol_overDraws(store_beta,store_Σ,setup_str,shSize=0.4);


@unpack n_fcst,n,p,nsave=setup_str

# n = setup_str.n
i_shock = 1;            # which variable you shock
plt = plot(layout=(3,2))
    for i_var = 1:n;
        plot!(plt,IRF_median[:,i_var,i_shock],ribbon = (IRF_median[:,i_var,i_shock].-IRF_68_low[:,i_var,1],IRF_68_high[:,i_var,i_shock].-IRF_median[:,i_var,1]),fillalpha=0.2,label="",title=var_names[i_var],titlefontsize=8,subplot=i_var)
    end
# end
display(plt)



#-------------------------------------------------#
#
#           Italy Baker
#
#-------------------------------------------------#


data_itB = data_ta_full[1:80,[:tpuBaker, :pmiIT, :gdpIT, :invIT, :hicpIT, :euribor]]
var_names = colnames(data_itB)
YYitB = values(data_itB);

setup_str, hyper_str = makeSetup(YYitB,"Chan2020_LBA_Minn",nburn=5000,nsave=5000,n_fcst=200)
store_beta, store_Σ = Chan2020_LBA_Minn(YYitB,setup_str,hyper_str);

# Using Impulse responses
IRF_median, IRF_68_low, IRF_68_high = irf_chol_overDraws(store_beta,store_Σ,setup_str,shSize=1.1);


@unpack n_fcst,n,p,nsave=setup_str

# n = setup_str.n
i_shock = 1;            # which variable you shock
plt = plot(layout=(3,2))
    for i_var = 1:n;
        plot!(plt,IRF_median[:,i_var,i_shock],ribbon = (IRF_median[:,i_var,i_shock].-IRF_68_low[:,i_var,1],IRF_68_high[:,i_var,i_shock].-IRF_median[:,i_var,1]),fillalpha=0.2,label="",title=var_names[i_var],titlefontsize=8,subplot=i_var)
    end
# end
display(plt)




#-------------------------------------------------#
#
#           Spain Caldara
#
#-------------------------------------------------#


data_es = data_ta_full[1:80,[:tpuCaldara, :pmiES, :gdpES, :invES, :hicpES, :euribor]]
var_names = colnames(data_es)
YYes = values(data_es);

setup_str, hyper_str = makeSetup(YYes,"Chan2020_LBA_Minn",nburn=5000,nsave=5000,n_fcst=200)
store_beta, store_Σ = Chan2020_LBA_Minn(YYes,setup_str,hyper_str);

# Using Impulse responses
IRF_median, IRF_68_low, IRF_68_high = irf_chol_overDraws(store_beta,store_Σ,setup_str,shSize=0.4);


@unpack n_fcst,n,p,nsave=setup_str

# n = setup_str.n
i_shock = 1;            # which variable you shock
plt = plot(layout=(3,2))
    for i_var = 1:n;
        plot!(plt,IRF_median[:,i_var,i_shock],ribbon = (IRF_median[:,i_var,i_shock].-IRF_68_low[:,i_var,1],IRF_68_high[:,i_var,i_shock].-IRF_median[:,i_var,1]),fillalpha=0.2,label="",title=var_names[i_var],titlefontsize=8,subplot=i_var)
    end
# end
display(plt)



#-------------------------------------------------#
#
#           Spain Baker
#
#-------------------------------------------------#


data_esB = data_ta_full[1:80,[:tpuBaker, :pmiES, :gdpES, :invES, :hicpES, :euribor]]
var_names = colnames(data_esB)
YYesB = values(data_esB);

setup_str, hyper_str = makeSetup(YYesB,"Chan2020_LBA_Minn",nburn=5000,nsave=5000,n_fcst=200)
store_beta, store_Σ = Chan2020_LBA_Minn(YYesB,setup_str,hyper_str);

# Using Impulse responses
IRF_median, IRF_68_low, IRF_68_high = irf_chol_overDraws(store_beta,store_Σ,setup_str,shSize=1.1);


@unpack n_fcst,n,p,nsave=setup_str

# n = setup_str.n
i_shock = 1;            # which variable you shock
plt = plot(layout=(3,2))
    for i_var = 1:n;
        plot!(plt,IRF_median[:,i_var,i_shock],ribbon = (IRF_median[:,i_var,i_shock].-IRF_68_low[:,i_var,1],IRF_68_high[:,i_var,i_shock].-IRF_median[:,i_var,1]),fillalpha=0.2,label="",title=var_names[i_var],titlefontsize=8,subplot=i_var)
    end
# end
display(plt)