include("devPkgs.jl")
using BEAVARs
using DelimitedFiles
using LinearAlgebra
using TimeSeries
using SparseArrays
using Plots


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
varNamesHF = varNamesM_full;
varNamesLF = [:gdpBG]
varOrder   = [:gdpBG,:survIndustryBG];
varOrder = [:gdpBG; varNamesHF]

# select only the needed data and transform it if needed
dataM_bg_tab = dataM_bg_raw[varNamesHF]./100
dataQ_bg_tab = percentchange(dataQ_bg_raw[varNamesLF])


# dataLF_tab = dataQ_bg_tab[1:4];
# dataHF_tab = dataM_bg_tab[1:16];
dataLF_tab = dataQ_bg_tab;
dataHF_tab = dataM_bg_tab;


# beyond this point we shouldn't need any more input from the user, i.e. we switch from Q and M (if you have monthly and quarterly data) or A and Q (if you have annual and quarterly) to HF and LF


n = 17
b0 = zeros(n,)
B0 = -1.0*I(n)
B1 = 1.0*I(n);# B1 = [0.26 0; 0 0.96]
B2 = zeros(size(B1))
B_draw = [b0 B1 B2 B2 B2 B2]
structB_draw = [B0 B_draw[:,2:end]]
Σt  = 1000*Matrix(I,n,n); Σt_inv = inv(Σt)
# Σt_inv = inv([0.00014309 0; 0 0.00047340])


# Y0[1:5,1]=[0.01404399 ;0.01404399 ;0.01404399 ;0.01404399 ;0.01404399 ]


# cB[n*p-n+1+n : Tfn]=b0[cB_b0_ind]



fdataHF_tab, z_tab, freq_mix_tp, datesHF, varNamesLF, fvarNames = BEAVARs.CPZ_prep_TA(dataLF_tab,dataHF_tab,varOrder)

YY = values(fdataHF_tab);
Tf,n = size(YY);
YYt, Y0, H_Bsp, sBd_ind, Σ_invsp, Σt_ind, cB_b0_ind, Xb, cB, Smsp, Sosp, Sm_bit, longyo, nm = BEAVARs.CPZ_initMat(YY,structB_draw,b0,Σt_inv,p);

M_zsp, z_vec, T_z, MOiM, MOiz = BEAVARs.CPZ_makeM_inter(z_tab,YYt,Sm_bit,datesHF,varNamesLF,fvarNames,freq_mix_tp,nm,Tf);



# updating cB
BEAVARs.CPZ_update_cB!(cB,B_draw[:,2:end],b0,Y0,cB_b0_ind,p,n)

# updating H_B
H_Bsp.nzval[:] = -structB_draw[sBd_ind];

# updating Σ_invsp
Σ_invsp.nzval[:] = Σt_inv[Σt_ind];

Gm = H_Bsp*Smsp
Go = H_Bsp*Sosp
Kym     = Gm'*Σ_invsp*Gm
CL = cholesky(Hermitian(Kym))
μ_y = CL.UP\(CL.PtL\Gm'*Σ_invsp)*(Xb*cB-Go*longyo)

# YYt[Sm_bit] = μ_y + CL.UP\randn(nm,)

# BEAVARs.CPZ_draw!(YYt,longyo,Y0,cB,B_draw,structB_draw,sBd_ind,Σt_inv,Σt_ind,Xb,cB_b0_ind,H_Bsp,Σ_invsp,p,n,Sm_bit,Smsp,Sosp,nm)

KymBar = MOiM + Kym;

CLBar = cholesky(Hermitian(KymBar))
# C = CLBar.PtL; # Ct = CLBar.UP;
μ_yBar = (CLBar.UP)\(CLBar.PtL\(MOiz + Kym*μ_y))

YYt[Sm_bit] = μ_yBar +  CLBar.UP\randn(nm,)

BEAVARs.CPZ_draw_wz!(YYt,longyo,Y0,cB,B_draw,structB_draw,sBd_ind,Σt_inv,Σt_ind,Xb,cB_b0_ind,H_Bsp,Σ_invsp,p,n,Sm_bit,Smsp,Sosp,nm,MOiM,MOiz);
# @btime BEAVARs.CPZ_draw_wz!(YYt,longyo,Y0,cB,B_draw,structB_draw,sBd_ind,Σt_inv,Σt_ind,Xb,cB_b0_ind,H_Bsp,Σ_invsp,p,n,Sm_bit,Smsp,Sosp,nm,MOiM,MOiz);

plot(YY)


hypSetup=hypChan2020()
Y, X, T = mlag_r(YY,p)

Yt = vec(Y)

(deltaP, sigmaP, mu_prior) = trainPriors(YY,p);

(idx_kappa1,idx_kappa2, V_Minn, beta_Minn) = prior_Minn(n,p,sigmaP,hypSetup);
k = n*p+1;
nk_speye = sparse(1:n*k,1:n*k,1.0); nk_speye.nzval[:] = 1.0./V_Minn;    # for updating prior in K_β
T_speye  = sparse(Matrix(1.0I, T, T)); # for updating XiSig. Ins't this just Σ_invsp without the t=0,...t=-p??? 
Xsur = SUR_form(X,n)
XiSig = Xsur'*kron(T_speye,sparse(1:n,1:n,1.0./sigmaP))

nk_speye.nzval[:] = 1.0./V_Minn;
K_β = nk_speye + XiSig*Xsur;
# K_β = sparse(1:n*k,1:n*k,1.0./V_Minn) + XiSig*Xsur;
cholK_β = cholesky(Hermitian(K_β));
beta_hat = cholK_β.UP\(cholK_β.PtL\(beta_Minn./V_Minn + XiSig * Yt))


# ndraws = nsave+nburn;
# store_beta=zeros(n^2*p+n,nsave)
beta = beta_hat + cholK_β.UP\randn(k*n);

k_nc = n*p; 
B_draw = reshape(beta,k_nc+1,n)'

b0[:] = B_draw[:,1]
blk[:,n+1:end] = B_draw[:,2:end]


xx = sprand(10,hel)

t1 = sparse(Matrix(1.0I, 5, 5))
t2 = sparse(Matrix(1.0I, 5, 5))
XC1.PtL\t2
XC1 = cholesky(t1)
XC2 = cholesky(t2)

ldiv!(XC1.PtL,t1)
ldiv!(XC1.PtL,rand(5,1))
ldiv!(XC1.PtL,rand(5,))

Y = zeros(5,)
ldiv!(Y,XC1.PtL,rand(5,))