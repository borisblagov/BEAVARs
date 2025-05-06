include("devPkgs.jl")
using BEAVARs
using DelimitedFiles
using LinearAlgebra
using TimeSeries
using SparseArrays
using Plots
using Distributions


p = 5; # number of lags
dataM_bg_full = readtimearray("data/dataM_BG.csv"; format="dd/mm/yyyy", delim=',')
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
dataM_bg_tab = dataM_bg_raw[varNamesHF]./100
dataQ_bg_tab = percentchange(dataQ_bg_raw[varNamesLF])


# dataLF_tab = dataQ_bg_tab[1:4];
# dataHF_tab = dataM_bg_tab[1:16];
dataLF_tab = dataQ_bg_tab;
dataHF_tab = dataM_bg_tab;


# beyond this point we shouldn't need any more input from the user, i.e. we switch from Q and M (if you have monthly and quarterly data) or A and Q (if you have annual and quarterly) to HF and LF


n = size(varOrder,1)

# hyperparameter setup
hypSetup=hypChan2020()
nu0 = 3; # hyperparameter on the inverse-wishart dist
b0 = zeros(n,)
B0 = -1.0*I(n)
B1 = 1.0*I(n);# B1 = [0.26 0; 0 0.96]
B2 = zeros(size(B1))
B_draw = [b0 B1 B2 B2 B2 B2]
structB_draw = [B0 B_draw[:,2:end]]
Σt  = .0010*Matrix(I,n,n).+0.00000000001; Σt_inv = inv(Σt)
# Σt_inv = inv([0.00014309 0; 0 0.00047340])
k = n*p+1; 


# Y0[1:5,1]=[0.01404399 ;0.01404399 ;0.01404399 ;0.01404399 ;0.01404399 ]


# cB[n*p-n+1+n : Tfn]=b0[cB_b0_ind]


## prepare missing data draw
fdataHF_tab, z_tab, freq_mix_tp, datesHF, varNamesLF, fvarNames = BEAVARs.CPZ_prep_TA(dataLF_tab,dataHF_tab,varOrder)

YYwNA = values(fdataHF_tab);
YYwNA[:,2:end]=randn(242,16)
YY = deepcopy(YYwNA)
model_stre = "CPZ2024";
setup_str, model_type = makeSetup(YY,model_stre,p,16,8,100,100)



Tf,n = size(YY);
YYt, Y0, H_Bsp, sBd_ind, Σ_invsp, Σt_ind, Σp_invsp, Σpt_ind, cB_b0_ind, Xb, cB, Smsp, Sosp, Sm_bit, longyo, nm, inputs_str = BEAVARs.CPZ_initMat(YY,structB_draw,b0,Σt_inv,p);

M_zsp, z_vec, T_z, MOiM, MOiz = BEAVARs.CPZ_makeM_inter(z_tab,YYt,Sm_bit,datesHF,varNamesLF,fvarNames,freq_mix_tp,nm,Tf);

# Initialize estimation stage
Y, X, T = mlag_r(YY,p);
Yt = vec(Y)
(deltaP, sigmaP, mu_prior) = trainPriors(YY,p);

(idx_kappa1,idx_kappa2, V_Minn, beta_Minn) = prior_Minn(n,p,sigmaP,hypSetup);

# Σhat_sp  = sparse(1:n,1:n,1.0); Σhat_sp.nzval[:] = 1.0./sigmaP; 
T_speye  = sparse(Matrix(1.0I, T, T)); # for updating XiSig. Ins't this just Σ_invsp without the t=0,...t=-p??? 
kronIT_Σhat_sp = kron(T_speye,Σt_inv);
Xsur = SUR_form(X,n)
Xsur_ind = repeat(1:T*(n*p+1),n); # indices from X to match to Xsur for updating Xsur
kronI_V_invsp = sparse(1:n*k,1:n*k,1.0);    # for updating prior in K_β

# draw of the missing values
BEAVARs.CPZ_draw_wz!(YYt,longyo,Y0,cB,B_draw,structB_draw,sBd_ind,Σt_inv,Σt_ind,Xb,cB_b0_ind,H_Bsp,Σ_invsp,p,n,Sm_bit,Smsp,Sosp,nm,MOiM,MOiz);
# @btime BEAVARs.CPZ_draw_wz!(YYt,longyo,Y0,cB,B_draw,structB_draw,sBd_ind,Σt_inv,Σt_ind,Xb,cB_b0_ind,H_Bsp,Σ_invsp,p,n,Sm_bit,Smsp,Sosp,nm,MOiM,MOiz);


# draw of the parameters
# beta,b0,B_draw,Σt_inv,structB_draw = BEAVARs.CPZ_Minn!(YY,p,hypSetup,nu0,n,k,b0,B_draw,Σt_inv,structB_draw);

Y, X, T = mlag_r(YY,p)
Yt = vec(Y)
(deltaP, sigmaP, mu_prior) = trainPriors(YY,p);

(idx_kappa1,idx_kappa2, V_Minn, beta_Minn) = prior_Minn(n,p,sigmaP,hypSetup);

# Σhat_sp  = sparse(1:n,1:n,1.0); Σhat_sp.nzval[:] = 1.0./sigmaP; 
# T_speye  = sparse(Matrix(1.0I, T, T)); # for updating XiSig. Ins't this just Σ_invsp without the t=0,...t=-p??? 
# Σhat_sp = Σt_inv;
# kronIT_Σhat_sp = kron(T_speye,Σt_inv);
Xsur = SUR_form(X,n)
Xsur.nzval[:] .= X[Xsur_ind];            # updating Xsur

Σp_invsp.nzval[:] = Σt_inv[Σpt_ind];
XtΣ_inv = Xsur'*Σp_invsp;

kronI_V_invsp.nzval[:] = 1.0./V_Minn;
K_β = kronI_V_invsp + XtΣ_inv*Xsur;
cholK_β = cholesky(Hermitian(K_β));
beta_hat = cholK_β.UP\(cholK_β.PtL\(beta_Minn./V_Minn + XtΣ_inv * Yt))

beta = beta_hat + cholK_β.UP\randn(k*n);

B_draw[:,:] = reshape(beta,k,n)'
b0[:] = B_draw[:,1]
structB_draw[:,n+1:end] = B_draw[:,2:end]


# plot(YY)

nsave = 100;
nburn = 100;
setup_str, model_type = makeSetup(YY,model_stre,p,16,8,nsave,nburn)
@btime store_YY, store_beta = BEAVARs.CPZ_loop!(YY, setup_str, hypSetup,nu0,k,b0,B_draw,Σt_inv,structB_draw,YYt,longyo,Y0,cB,sBd_ind,Σt_ind,Xb,cB_b0_ind,H_Bsp,Σ_invsp,Sm_bit,Smsp,Sosp,nm,MOiM,MOiz,Tf,Xsur_ind,kronI_V_invsp,Xsur,Σp_invsp,Σpt_ind);
@btime store_YY, store_beta = BEAVARs.CPZ_loop_old!(YY,setup_str,hypSetup,nu0,k,b0,B_draw,Σt_inv,structB_draw,YYt,longyo,Y0,cB,sBd_ind,Σt_ind,Xb,cB_b0_ind,H_Bsp,Σ_invsp,Sm_bit,Smsp,Sosp,nm,MOiM,MOiz,nburn,nsave,Tf);


# nsave = 5000;
# nburn = 5000;
# @btime store_YY, store_beta = BEAVARs.CPZ_loop!(YY,p,hypSetup,nu0,n,k,b0,B_draw,Σt_inv,structB_draw,YYt,longyo,Y0,cB,sBd_ind,Σt_ind,Xb,cB_b0_ind,H_Bsp,Σ_invsp,Sm_bit,Smsp,Sosp,nm,MOiM,MOiz,nburn,nsave,Tf);
