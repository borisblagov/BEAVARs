include("devPkgs.jl")
using BEAVARs
using DelimitedFiles
using LinearAlgebra
using TimeSeries
using SparseArrays
using Plots
using Distributions


p = 5; # number of lags
dataM_bg_full = readtimearray("data/bg_julL.csv"; format="dd/mm/yyyy", delim=',')
dataQ_bg_full = readtimearray("data/dataQ_BG.csv"; format="dd/mm/yyyy", delim=',')
varNamesM_full = colnames(dataM_bg_full)


#adding second GDP for testing
# bgTest_tab = dataQ_bg_full[:gdpBG];
# bgTest_tab = rename(bgTest_tab,:gdp2)
# dataQ_bg_full = merge(dataQ_bg_full,bgTest_tab)

dataM_bg_raw = dataM_bg_full
dataQ_bg_raw = dataQ_bg_full


varNamesM_full = colnames(dataM_bg_full)
# varNamesQ_full = colnames(dataQ_bg_full)

varNamesHF = [:survIndustryBG];
varNamesLF = [:gdpBG]
varOrder   = [:gdpBG,:survIndustryBG,:ipBG];
varNamesHF = varNamesM_full;
varOrder = [:gdpBG; varNamesHF]

# select only the needed data and transform it if needed
dataM_bg_tab = dataM_bg_raw[varNamesHF];
dataQ_bg_tab = percentchange(dataQ_bg_raw[varNamesLF])


# dataLF_tab = dataQ_bg_tab[1:4];
# dataHF_tab = dataM_bg_tab[1:16];
dataLF_tab = dataQ_bg_tab;
dataHF_tab = dataM_bg_tab;


# beyond this point we shouldn't need any more input from the user, i.e. we switch from Q and M (if you have monthly and quarterly data) or A and Q (if you have annual and quarterly) to HF and LF


## prepare missing data draw
nsave = 100;
nburn = 100;
model_stre = "CPZ2024";

setup_str, model_type = makeSetup(varOrder,model_stre,p,16,8,nsave,nburn)
# hyperparameter setup
hyper_str=hypChan2020()



fdataHF_tab, z_tab, freq_mix_tp, datesHF, varNamesLF, fvarNames = BEAVARs.CPZ_prep_TimeArrays(dataLF_tab,dataHF_tab,varOrder)

YYwNA = values(fdataHF_tab);
YY = deepcopy(YYwNA)
Tf,n = size(YY);

B_draw, structB_draw, Σt_inv, b0 = BEAVARs.initParamMatrices(n,p,1) 

YYt, Y0, longyo, nm, H_B, H_B_CI, strctBdraw_LI, Σ_invsp, Σt_LI, Σp_invsp, Σpt_ind, Xb, cB, cB_b0_LI, Smsp, Sosp, Sm_bit, Gm, Go, GΣ, Kym = BEAVARs.CPZ_initMatrices(YY,structB_draw,b0,Σt_inv,p);

M_zsp, z_vec, T_z, MOiM, MOiz = BEAVARs.CPZ_makeM_inter(z_tab,YYt,Sm_bit,datesHF,varNamesLF,fvarNames,freq_mix_tp,nm,Tf);
# YYt = BEAVARs.CPZ_draw_wz3!(YYt,longyo,Y0,cB,B_draw,structB_draw,strctBdraw_LI,Σt_inv,Σt_LI,Xb,cB_b0_LI,Σ_invsp,p,n,Sm_bit,Smsp,Sosp,nm,MOiM,MOiz,Gm,Go,H_B,GΣ,Kym,H_B_CI);
YYt = BEAVARs.CPZ_draw_wz3!(YYt,longyo,Y0,cB,B_draw,structB_draw,strctBdraw_LI,Σt_inv,Σt_LI,Xb,cB_b0_LI,Σ_invsp,p,n,Sm_bit,Smsp,Sosp,nm,MOiM,MOiz,Gm,Go,H_B,GΣ,Kym,H_B_CI);
     
upd_these_vec = sum(Sm_bit,dims=2).>T*0.25;
Y, X, T, deltaP, sigmaP, mu_prior, V_Minn_inv, V_Minn_inv_elview, XtΣ_inv_den, XtΣ_inv_X, Xsur_den, Xsur_CI, X_CI, k, K_β, beta = BEAVARs.CPZ_initMinn(YY,p);
beta,b0,B_draw,Σt_inv,structB_draw, XtΣ_inv_X = BEAVARs.CPZ_Minn5!(YY,p,hyper_str,n,k,b0,B_draw,Σt_inv,structB_draw,Σp_invsp,Σpt_ind,Y,X,T,mu_prior,deltaP,sigmaP,intercept,Xsur_den,Xsur_CI,X_CI,XtΣ_inv_den,XtΣ_inv_X,V_Minn_inv,V_Minn_inv_elview,upd_these_vec,K_β,beta);

@time BEAVARs.CPZ_Minn5!(YY,p,hyper_str,n,k,b0,B_draw,Σt_inv,structB_draw,Σp_invsp,Σpt_ind,Y,X,T,mu_prior,deltaP,sigmaP,intercept,Xsur_den,Xsur_CI,X_CI,XtΣ_inv_den,XtΣ_inv_X,V_Minn_inv,V_Minn_inv_elview,upd_these_vec,K_β,beta);
@time BEAVARs.CPZ_Minn4!(YY,p,hyper_str,n,k,b0,B_draw,Σt_inv,structB_draw,Σp_invsp,Σpt_ind,Y,X,T,mu_prior,deltaP,sigmaP,intercept,Xsur_den,Xsur_CI,X_CI,XtΣ_inv_den,XtΣ_inv_X,V_Minn_inv,V_Minn_inv_elview);


#---------------------------------------
store_YY, store_beta, store_Σt = BEAVARs.CPZ_loop!(dataLF_tab,dataHF_tab,varOrder,setup_str,hyper_str);
store_YY2, store_beta2 = BEAVARs.CPZ_loop2!(YYwNA,setup_str,hyper_str,MOiM, MOiz);
# @time store_YY, store_beta, store_Σt = BEAVARs.CPZ_loop!(dataLF_tab,dataHF_tab,varOrder,setup_str,hyper_str);
# @time store_YY2, store_beta2 = BEAVARs.CPZ_loop2!(YYwNA,setup_str,hyper_str,MOiM, MOiz);

yy1 = median(store_YY,dims=3);
plot(yy1[:,1])
yy2 = median(store_YY2,dims=3);
plot!(yy2[:,1])

# plot(M_zsp*YYt[Sm_bit])
plot(M_zsp*yy1[Sm_bit'])
plot!(z_vec)
# #---------------------------------------

