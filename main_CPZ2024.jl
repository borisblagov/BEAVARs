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
varOrder   = [:gdpBG,:survIndustryBG];
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

n = size(varOrder,1)

# hyperparameter setup
hyper_str=hypChan2020()
k = n*p+1; 


## prepare missing data draw
fdataHF_tab, z_tab, freq_mix_tp, datesHF, varNamesLF, fvarNames = BEAVARs.CPZ_prep_TA(dataLF_tab,dataHF_tab,varOrder)

YYwNA = values(fdataHF_tab);


nsave = 50;
nburn = 50;
YY = deepcopy(YYwNA)
Tf,n = size(YY);
model_stre = "CPZ2024";
setup_str, model_type = makeSetup(YY,model_stre,p,16,8,nsave,nburn)

B_draw = [zeros(n,)  1.0*I(n) zeros(n,n*(p-1))]; b0 = B_draw[:,1];
structB_draw = [-1.0*I(n) B_draw[:,2:end]]
Σt  = .001*Matrix(I,n,n).+0.00000000001; Σt_inv = inv(Σt)

YYt, Y0, H_Bsp, sBd_ind, Σ_invsp, Σt_ind, Σp_invsp, Σpt_ind, cB_b0_ind, Xb, cB, Smsp, Sosp, Sm_bit, longyo, nm, inputs_str, H_B,H_B_CI = BEAVARs.CPZ_initMatrices(YY,structB_draw,b0,Σt_inv,p);

M_zsp, z_vec, T_z, MOiM, MOiz = BEAVARs.CPZ_makeM_inter(z_tab,YYt,Sm_bit,datesHF,varNamesLF,fvarNames,freq_mix_tp,nm,Tf);


#---------------------------------------
intercept = 1;
store_YY1, store_beta = BEAVARs.CPZ_loop5!(YYwNA,setup_str,hyper_str,k,Tf,intercept,MOiM, MOiz);
# @btime store_YY1, store_beta = BEAVARs.CPZ_loop5!(YYwNA,setup_str,hyper_str,k,Tf,intercept,MOiM, MOiz);
store_YY2, store_beta = BEAVARs.CPZ_loop6!(YYwNA,setup_str,hyper_str,k,Tf,intercept,MOiM, MOiz);
# @btime store_YY2, store_beta = BEAVARs.CPZ_loop6!(YYwNA,setup_str,hyper_str,k,Tf,intercept,MOiM, MOiz);
store_YY3, store_beta2 = BEAVARs.CPZ_loop7!(YYwNA,setup_str,hyper_str,k,Tf,intercept,MOiM, MOiz);
@btime store_YY3, store_beta = BEAVARs.CPZ_loop7!(YYwNA,setup_str,hyper_str,k,Tf,intercept,MOiM, MOiz);

yy1 = median(store_YY1,dims=3)
yy2 = median(store_YY2,dims=3)
yy3 = median(store_YY3,dims=3)
plot(yy1[:,1])
plot!(yy2[:,1])
plot!(yy3[:,1])
# #---------------------------------------


# @time store_YY, store_beta = BEAVARs.CPZ_loop!(YY, setup_str, hypSetup,nu0,k,b0,B_draw,Σt_inv,structB_draw,YYt,longyo,Y0,cB,sBd_ind,Σt_ind,Xb,cB_b0_ind,H_Bsp,Σ_invsp,Sm_bit,Smsp,Sosp,nm,MOiM,MOiz,Tf,Xsur_ind,kronI_V_invsp,Xsur,Σp_invsp,Σpt_ind);
# @time store_YY, store_beta = BEAVARs.CPZ_loop2!(YY, setup_str, hypSetup,nu0,k,b0,B_draw,Σt_inv,structB_draw,YYt,longyo,Y0,cB,sBd_ind,Σt_ind,Xb,cB_b0_ind,H_Bsp,Σ_invsp,Sm_bit,Smsp,Sosp,nm,MOiM,MOiz,Tf,Xsur_ind,kronI_V_invsp,Xsur,Σp_invsp,Σpt_ind,Y,X,T,mu_prior,deltaP,sigmaP,intercept);
# @btime store_YY, store_beta = BEAVARs.CPZ_loop3!(YY, setup_str, hypSetup,nu0,k,b0,B_draw,Σt_inv,structB_draw,YYt,longyo,Y0,cB,sBd_ind,Σt_ind,Xb,cB_b0_ind,H_Bsp,Σ_invsp,Sm_bit,Smsp,Sosp,nm,MOiM,MOiz,Tf,kronI_V_invsp,Σp_invsp,Σpt_ind,Y,X,T,mu_prior,deltaP,sigmaP,intercept,Xsur_den,Xsur_CI,X_CI,XtΣ_inv_den,XtΣ_inv_X);
# @btime store_YY2, store_beta2 = BEAVARs.CPZ_loop4!(YY, setup_str, hypSetup,nu0,k,b0,B_draw,Σt_inv,structB_draw,YYt,longyo,Y0,cB,sBd_ind,Σt_ind,Xb,cB_b0_ind,H_Bsp,Σ_invsp,Sm_bit,Smsp,Sosp,nm,MOiM,MOiz,Tf,kronI_V_invsp,Σp_invsp,Σpt_ind,Y,X,T,mu_prior,deltaP,sigmaP,intercept,Xsur_den,Xsur_CI,X_CI,XtΣ_inv_den,XtΣ_inv_X,V_Minn_inv,V_Minn_inv_elview);
# # @btime store_YY, store_beta = BEAVARs.CPZ_loop_old!(YY,setup_str,hypSetup,nu0,k,b0,B_draw,Σt_inv,structB_draw,YYt,longyo,Y0,cB,sBd_ind,Σt_ind,Xb,cB_b0_ind,H_Bsp,Σ_invsp,Sm_bit,Smsp,Sosp,nm,MOiM,MOiz,nburn,nsave,Tf);

#----------------------------


#  YYt = BEAVARs.CPZ_draw_wz!(YYt,longyo,Y0,cB,B_draw,structB_draw,sBd_ind,Σt_inv,Σt_ind,Xb,cB_b0_ind,H_Bsp,Σ_invsp,p,n,Sm_bit,Smsp,Sosp,nm,MOiM,MOiz);

# H_B = Matrix(H_Bsp);
# Gm = H_B*Smsp; Go = H_B*Sosp; GΣ = Gm'*Σ_invsp; Kym = GΣ*Gm; # we can initialize all these and then mutate with mul!()
# # Initialize matrices for updating the parameter draws from CPZ_Minn  
# # ------------------------------------

# (deltaP, sigmaP, mu_prior) = trainPriors(YY,p);                                     # do OLS to initialize priors
# V_Minn_inv          = 1.0*Matrix(I,n*k,n*k);                                        # prior matrix
# V_Minn_inv_elview   = @view(V_Minn_inv[diagind(V_Minn_inv)]);                       # will be used to update the diagonal



# Y, X, T             = mlagL(YY,p);
# Yt                  = vec(Y)

# XtΣ_inv_den         = zeros(k*n,T*n);                   # this is X' ( I(T) ⊗ Σ-1 )   from page 6 in Chan 2020 LBA
# XtΣ_inv_X           = zeros(n*k,n*k);                   # this is X' ( I(T) ⊗ Σ-1 ) X from page 6 in Chan 2020 LBA    
# Xsur_den, Xsur_CI, X_CI = BEAVARs.SUR_form_dense(X,n);  # prepares the SUR form and the indices of the parameters for updating


# BEAVARs.CPZ_draw_wz3!(YYt,longyo,Y0,cB,B_draw,structB_draw,sBd_ind,Σt_inv,Σt_ind,Xb,cB_b0_ind,H_Bsp,Σ_invsp,p,n,Sm_bit,Smsp,Sosp,nm,MOiM,MOiz,Gm,Go,H_B,GΣ,Kym,H_B_CI);

# BEAVARs.CPZ_draw_wz2!(YYt,longyo,Y0,cB,B_draw,structB_draw,sBd_ind,Σt_inv,Σt_ind,Xb,cB_b0_ind,H_Bsp,Σ_invsp,p,n,Sm_bit,Smsp,Sosp,nm,MOiM,MOiz,Gm,Go,H_B)
    
# BEAVARs.CPZ_draw_wz!(YYt,longyo,Y0,cB,B_draw,structB_draw,sBd_ind,Σt_inv,Σt_ind,Xb,cB_b0_ind,H_Bsp,Σ_invsp,p,n,Sm_bit,Smsp,Sosp,nm,MOiM,MOiz)
    