include("devPkgs.jl")
using Parameters
using BEAVARs
using TimeSeries
using DelimitedFiles
using Plots
using Statistics
using LinearAlgebra
using Distributions


YY20 = readdlm("data/FRED_Chan2020_LBA.txt", ',');
YY = YY20; n = size(YY20,2);
# YY = YY20[:,[1,4,5,6]]
@btime out_struct, var_strcut, hyp_strct = beavar("Chan2020iniw",YY,n_save = 10,n_burn=10);
@btime out_struct2, var_strcut, hyp_strct = beavar("Chan2020iniw2",YY,n_save = 10,n_burn=10);
@unpack p,nburn,nsave = var_strcut
@unpack nu0 = hyp_strct
p = 4

Y, X, T, n, sigmaP, S_0, Σt_inv, Vβminn_inv, Vβminn_inv_elview, Σ_invsp, Σt_LI, XtΣ_inv_den, XtΣ_inv_X, Xsur_den, Xsur_CI, X_CI, k, K_β, beta, intercept = BEAVARs.initMinn(YY,p);

(idx_kappa1,idx_kappa2, Vβminn, βminn) = prior_Minn(n,p,sigmaP,hyp_strct)


Vβminn_inv_elview[:] = 1.0./Vβminn;

Xsur_den[Xsur_CI] = X[X_CI];
Σ_invsp.nzval[:] = Σt_inv[Σt_LI];

mul!(XtΣ_inv_den,Xsur_den',Σ_invsp);            #  X'*( I(T) ⊗ Σ^{-1} )
mul!(XtΣ_inv_X,XtΣ_inv_den,Xsur_den);           #  X'*( I(T) ⊗ Σ^{-1} )*X
K_β[:,:] .= Vβminn_inv .+ XtΣ_inv_X;            #  K_β = V^{-1} + X'*( I(T) ⊗ Σ^{-1} )*X
prior_mean = Vβminn_inv*βminn;                  #  V^-1 * βminn 
mul!(prior_mean,XtΣ_inv_den, vec(Y'),1.0,1.0);  # (V^-1_Minn * beta_Minn) + X' ( I(T) ⊗ Σ-1 ) y
cholK_β = cholesky(Hermitian(K_β));
beta_hat = ldiv!(cholK_β.U,ldiv!(cholK_β.L,prior_mean));
beta = beta_hat + ldiv!(cholK_β.U,randn(k*n,));

U = reshape(vec(Y') - Xsur_den*beta,n,T);
Σt = rand(InverseWishart(nu0+n+T,S_0+U*U'));
Σt_inv = Σt\I;
