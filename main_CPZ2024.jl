include("devPkgs.jl")
using BEAVARs
using DelimitedFiles
using LinearAlgebra
using TimeSeries
using SparseArrays
using Plots
using Parameters
using Distributions


p = 5; # number of lags
dataM_bg_full = readtimearray("data/bg_julL.csv"; format="dd/mm/yyyy", delim=',')
dataQ_bg_full = readtimearray("data/dataQ_BG.csv"; format="dd/mm/yyyy", delim=',')
varNamesM_full = colnames(dataM_bg_full)

dataM_bg_raw = dataM_bg_full
dataQ_bg_raw = dataQ_bg_full


varNamesM_full = colnames(dataM_bg_full)

varNamesHF = [:survIndustryBG];
varNamesLF = [:gdpBG]
varOrder   = [:gdpBG,:survIndustryBG,:ipBG];
# varNamesHF = varNamesM_full;
varList = [:gdpBG; varNamesHF]

# select only the needed data and transform it if needed
dataM_bg_tab = dataM_bg_raw[varNamesHF];
dataQ_bg_tab = percentchange(dataQ_bg_raw[varNamesLF])


dataLF_tab = dataQ_bg_tab;
dataHF_tab = dataM_bg_tab;


# beyond this point we shouldn't need any more input from the user, i.e. we switch from Q and M (if you have monthly and quarterly data) or A and Q (if you have annual and quarterly) to HF and LF


## prepare missing data draw
n_save = 10;
n_burn = 10;
model_str = "CPZ2024";

out_strct, varSetup,hypSetup = beavar("CPZ2024",dataHF_tab,dataLF_tab,varList,0,n_burn=100,n_save=10);

# hyperparameter setup
hyper_str=hypChan2020()

set_strct = VARSetup(p,n_save,n_burn,16,8,1);

@unpack p, nburn,nsave, const_loc = set_strct
ndraws = nsave+nburn;
nmdraws = 10;               # given a draw from the parameters to draw multiple time from the distribution of the missing data for better confidence intervals
trans = 0;
fdataHF_tab, z_tab, freq_mix_tp, datesHF, varNamesLF, fvarNames = BEAVARs.CPZ_prep_TimeArrays(dataLF_tab,dataHF_tab,varList,trans)

YYwNA = values(fdataHF_tab);
YY = deepcopy(YYwNA);
Tf,n = size(YY);

B_draw, structB_draw, Σt_inv, b0 = BEAVARs.initParamMatrices(n,p,const_loc) 

YYt, Y0, longyo, nm, H_B, H_B_CI, strctBdraw_LI, Σ_invsp, Σt_LI, Σp_invsp, Σpt_ind, Xb, cB, cB_b0_LI, Smsp, Sosp, Sm_bit, Gm, Go, GΣ, Kym = BEAVARs.CPZ_initMatrices(YY,structB_draw,b0,Σt_inv,p);

M_zsp, z_vec, T_z, MOiM, MOiz = BEAVARs.CPZ_makeM_inter(z_tab,YYt,Sm_bit,datesHF,varNamesLF,fvarNames,freq_mix_tp,nm,Tf);

# YY has missing values so we need to draw them once to be able to initialize matrices and prior values
YYt = BEAVARs.CPZ_draw_wz!(YYt,longyo,Y0,cB,B_draw,structB_draw,strctBdraw_LI,Σt_inv,Σt_LI,Xb,cB_b0_LI,Σ_invsp,p,n,Sm_bit,Smsp,Sosp,nm,MOiM,MOiz,Gm,Go,H_B,GΣ,Kym,H_B_CI,nmdraws);

# we will be updating the priors for variables with many missing observations (>25%)
updP_vec = sum(Sm_bit,dims=2).>size(Sm_bit,2)*0.25;

# Initialize matrices for updating the parameter draws from CPZ_Minn  
# ------------------------------------
Y, X, T, deltaP, sigmaP, mu_prior, V_Minn_inv, V_Minn_inv_elview, XtΣ_inv_den, XtΣ_inv_X, Xsur_den, Xsur_CI, X_CI, k, K_β, beta, intercept = BEAVARs.CPZ_initMinn(YY,p)


# prepare matrices for storage
store_YY    = zeros(Tf,n,nsave);
store_β     = zeros(n^2*p+n,nsave);
store_Σt_inv= zeros(n,n,nsave);
store_Σt    = zeros(n,n,nsave);

# @showprogress for ii in 1:ndraws
    ii=1
    # draw of the missing values
    BEAVARs.CPZ_draw_wz!(YYt,longyo,Y0,cB,B_draw,structB_draw,strctBdraw_LI,Σt_inv,Σt_LI,Xb,cB_b0_LI,Σ_invsp,p,n,Sm_bit,Smsp,Sosp,nm,MOiM,MOiz,Gm,Go,H_B,GΣ,Kym,H_B_CI,nmdraws);
    
    # draw of the parameters
    beta,b0,B_draw,Σt_inv,structB_draw,Σt = BEAVARs.CPZ_Minn!(YY,p,hypSetup,n,k,b0,B_draw,Σt_inv,structB_draw,Σp_invsp,Σpt_ind,Y,X,T,mu_prior,deltaP,sigmaP,const_loc,Xsur_den,Xsur_CI,X_CI,XtΣ_inv_den,XtΣ_inv_X,V_Minn_inv,V_Minn_inv_elview,updP_vec,K_β,beta);
    # beta,b0,B_draw,Σt_inv,structB_draw = BEAVARs.CPZ_Minn4!(YY,p,hypSetup,n,k,b0,B_draw,Σt_inv,structB_draw,Σp_invsp,Σpt_ind,Y,X,T,mu_prior,deltaP,sigmaP,const_loc,Xsur_den,Xsur_CI,X_CI,XtΣ_inv_den,XtΣ_inv_X,V_Minn_inv,V_Minn_inv_elview);

    Y, X = mlagL!(YY,Y,X,p,n)
    (deltaP, sigmaP, mu_prior) = BEAVARs.updatePriors3!(Y,X,n,mu_prior,deltaP,sigmaP,intercept,updP_vec);

    (idx_kappa1,idx_kappa2, V_Minn_vec, beta_Minn) = prior_Minn(n,p,sigmaP,hypSetup);
    V_Minn_vec_inv = 1.0./V_Minn_vec;
    Σp_invsp.nzval[:] = Σt_inv[Σpt_ind];
   
    
    Xsur_den[Xsur_CI] = X[X_CI]; 
    mul!(XtΣ_inv_den,Xsur_den',Σp_invsp);                 # X' ( I(T) ⊗ Σ-1 )
    mul!(XtΣ_inv_X,XtΣ_inv_den,Xsur_den);                   # XtΣ_inv_X = X' ( I(T) ⊗ Σ-1 ) X
    V_Minn_inv_elview[:] = V_Minn_vec_inv;  # update the diagonal of V_Minn_inv, i.e. V_Minn^-1
    K_β[:,:] .= V_Minn_inv .+ XtΣ_inv_X;
    cholK_β = cholesky(Hermitian(K_β));    
    
    prior_mean = V_Minn_inv*beta_Minn;                   # V^-1_Minn * beta_Minn 
    mul!(prior_mean,XtΣ_inv_den,  vec(Y),1.0,1.0);        # (V^-1_Minn * beta_Minn) + X' ( I(T) ⊗ Σ-1 ) y


    beta_hat = ldiv!(cholK_β.U,ldiv!(cholK_β.L,prior_mean));

    helper_vec = randn(k*n,);
    beta = beta_hat + ldiv!(cholK_β.U,helper_vec);

    B_draw[:,:] = reshape(beta,k,n)'
    b0[:] = B_draw[:,1]
    structB_draw[:,n+1:end] = B_draw[:,2:end]


    # errors 
    fit = zeros(size(Y))
    ee = Y-mul!(fit,X,B_draw');
    # Σ_t = rand(InverseWishart(T+hypSetup.nu0,Diagonal(sigmaP)+ee'*ee));
    Σ_t = rand(InverseWishart(T+hypSetup.nu0+n, ee'*ee));
    Σt_inv[:,:] = Σ_t\I;

    # if ii>nburn
    #     store_β[:,ii-nburn]  = beta;
    #     store_YY[:,:,ii-nburn]  = YY;
    #     store_Σt_inv[:,:,ii-nburn]    = Σt_inv;
    #     store_Σt[:,:,ii-nburn] = Σt;
    # end
# end