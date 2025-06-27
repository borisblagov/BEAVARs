include("devPkgs.jl")
using BEAVARs
using TimeSeries
using Parameters
using DelimitedFiles
using Plots
using Statistics
using LinearAlgebra
using Distributions
using SparseArrays

include("src/plot_functions.jl")
trans = 0;
dataHF_tab, dataLF_tab, varList = BEAVARs.readSpec("bg_julL","data/Specifications_mfvar.xlsx");
out_strct, varSetup,hypSetup = beavar("Blagov2025",dataHF_tab,dataLF_tab,varList,0,n_burn=100,n_save=100);
fanChart(out_strct.store_YY[:,1,:])
Yfor3D = BEAVARs.forecast(out_strct,varSetup);
fanChart(Yfor3D[:,1,:])


YY_HF_med = percentile_mat(out_strct.store_YY,0.5,dims=3);
@unpack M_zsp, store_YY, z_vec, Sm_bit = out_strct
plot(M_zsp*YY_HF_med'[Sm_bit])
plot!(z_vec)


varSetup, hypSetup =  BEAVARs.beavar_debug("Blagov2025",n_save=10,n_burn=10)


out = BEAVARs.Blagov2025(dataHF_tab,dataLF_tab,varList,varSetup,hypSetup,trans);


@unpack p, nburn,nsave, const_loc = varSetup

@unpack ρ, σ_h2, v_h0, S_h0, ρ_0, V_ρ = hypSetup

# nmdraws = 10;               # given a draw from the parameters to draw multiple time from the distribution of the missing data for better confidence intervals

ndraws = nsave+nburn;

trans = 0;

fdataHF_tab, z_tab, freq_mix_tp, datesHF, varNamesLF, fvarNames = BEAVARs.CPZ_prep_TimeArrays(dataLF_tab,dataHF_tab,varList,trans)

YYwNA = values(fdataHF_tab);
YY = deepcopy(YYwNA);
Tf,n = size(YY);

B_draw, structB_draw, Σt_inv, b0 = BEAVARs.initParamMatrices(n,p,const_loc) 

YYt, Y0, longyo, nm, H_B, H_B_CI, strctBdraw_LI, Σ_invsp, Σt_LI, Σp_invsp, Σpt_ind, Xb, cB, cB_b0_LI, Smsp, Sosp, Sm_bit, Gm, Go, GΣ, Kym = BEAVARs.CPZ_initMatrices(YY,structB_draw,b0,Σt_inv,p);

M_zsp, z_vec, T_z, MOiM, MOiz = BEAVARs.CPZ_makeM_inter(z_tab,YYt,Sm_bit,datesHF,varNamesLF,fvarNames,freq_mix_tp,nm,Tf);


# YY has missing values so we need to draw them once to be able to initialize matrices and prior values
cB,H_B,Σ_invsp  = BEAVARs.Blagov2025_updCPZ!(cB,H_B,Σ_invsp,B_draw,structB_draw,Σt_inv,Y0,cB_b0_LI,p,n,H_B_CI,strctBdraw_LI,Σt_LI);
YYt             = BEAVARs.Blagov2025_draw_wz!(YYt,longyo,cB,Σ_invsp,Sm_bit,Smsp,Sosp,nm,MOiM,MOiz,Gm,Go,H_B,GΣ,Kym);

# BEAVARs.CPZ_draw_wz!(YYt,longyo,Y0,cB,B_draw,structB_draw,strctBdraw_LI,Σt_inv,Σt_LI,Xb,cB_b0_LI,Σ_invsp,p,n,Sm_bit,Smsp,Sosp,nm,MOiM,MOiz,Gm,Go,H_B,GΣ,Kym,H_B_CI,nmdraws);


# Initialize matrices for updating the parameter draws from CPZ_iniv  
# ------------------------------------
Y, X, T, k, sigmaP, S_0, Σ, A_0, V_Ainv, v_0, H_ρ,h,eh,Ωinv, dg_ind_Ωinv, VAinvDA0, AVAinvA   = BEAVARs.Blagov2025_initcsv(YY,p,hypSetup);


# prepare matrices for storage
store_YY    = zeros(Tf,n,nsave);


Y, X = mlagL!(YY,Y,X,p,n);
A, cholΣU, Σ, s2_h, U, Ωinv = BEAVARs.Chan2020_drawA(Y,X,n,k,T,v_0,h,Ωinv,dg_ind_Ωinv,V_Ainv,S_0,VAinvDA0,AVAinvA);
h = BEAVARs.Chan2020_draw_h!(h,s2_h,ρ,σ_h2,n,H_ρ,T);
ρ, σ_h2, eh = BEAVARs.Chan2020_draw_ρ!(ρ,h,eh,v_h0,S_h0,ρ_0,V_ρ,T);

Σt_inv = inv(Σ) 
B_draw[:,:] = A';
b0[:] = B_draw[:,1];
structB_draw[:,n+1:end] = B_draw[:,2:end];

cB,H_B,Σ_invsp  = BEAVARs.Blagov2025_updCPZ!(cB,H_B,Σ_invsp,B_draw,structB_draw,Σt_inv,Y0,cB_b0_LI,p,n,H_B_CI,strctBdraw_LI,Σt_LI);
Σ_invsp         = BEAVARs.Blagov2025_updCPZcsv!(Σ_invsp,p,n,Tf,Ωinv);
YYt             = BEAVARs.Blagov2025_draw_wz!(YYt,longyo,cB,Σ_invsp,Sm_bit,Smsp,Sosp,nm,MOiM,MOiz,Gm,Go,H_B,GΣ,Kym);