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


# varNamesM_full = colnames(dataM_bg_full)
# varNamesQ_full = colnames(dataQ_bg_full)

varNamesHF = [:survIndustryBG];
varNamesLF = [:gdpBG]
varOrder   = [:gdpBG,:survIndustryBG]

# select only the needed data and transform it if needed
dataM_bg_tab = dataM_bg_raw[varNamesHF]./100
dataQ_bg_tab = percentchange(dataQ_bg_raw[varNamesLF])


# dataLF_tab = dataQ_bg_tab[1:4];
# dataHF_tab = dataM_bg_tab[1:16];
dataLF_tab = dataQ_bg_tab;
dataHF_tab = dataM_bg_tab;


# beyond this point we shouldn't need any more input from the user, i.e. we switch from Q and M (if you have monthly and quarterly data) or A and Q (if you have annual and quarterly) to HF and LF


n = 2
b0 = zeros(n,)
B0 = -1.0*I(n)
B1 = 1.0*I(n); B1 = [0.26 0; 0 0.96]
B2 = zeros(size(B1))
blk = [B0 B1 B2 B2*0.1 B2*0.01 B2*0.01]
Σt_inv  = 100*Matrix(I,n,n); Σt_inv = inv([0.00014309 0; 0 0.00047340])


# Y0[1:5,1]=[0.01404399 ;0.01404399 ;0.01404399 ;0.01404399 ;0.01404399 ]


# cB[n*p-n+1+n : Tfn]=b0[cB_b0_ind]


fdataHF_tab, z_tab, freq_mix_tp, datesHF, varNamesLF, fvarNames = BEAVARs.CPZ_prep_TAfrequencies(dataLF_tab,dataHF_tab,varOrder)


YY = values(fdataHF_tab)
Y0 = @views YY[1:p,:]

(Tf,n) = size(YY); # full time span (with initial conditions)
k = n*(p+1); kn = k*n
Tfn = n*Tf;

YYt = (YY');
vYYt = vec(YYt);


Sm_bit = isnan.(YYt)
So_bit = .!isnan.(YYt)
longyo = vYYt[vec(So_bit)];

indC_nan_wide = findall(Sm_bit) #  Cartesian indices of missing values
# indC_non_wide = findall(!isnan,YYt)  # Cartesian indices of not missing values
indC_non_wide = findall(So_bit)  # Cartesian indices of not missing values

# convert between linear and cartesian indices
indL_all = LinearIndices(YYt);
indL_nan_wide = indL_all[indC_nan_wide] # are the linear indices of NaN values
indL_non_wide = indL_all[indC_non_wide] # are the linear indices of non NaN values



nm = sum(Sm_bit)
S_full = I(Tf*n);
Sm = S_full[:,indL_nan_wide]
So = S_full[:,indL_non_wide]
Smsp = sparse(Sm);
Sosp = sparse(So);

# Initialize matrices
H_Bsp, blk_ind = BEAVARs.makeBlkDiag(Tfn,n,p,-blk);
Σsp_inv, Σt_ind = BEAVARs.makeBlkDiag(Tfn,n,0,Σt_inv);
cB_b0_ind = repeat(1:n,div(Tfn-n*p-n+1+n,2));  # this repeats [1:n] so that we can update cB[indicesAfter Y_0,Y_{-1}, ..., Ymp] = b0[cB_b0_ind]
Xb = sparse(Matrix(1.0I, Tfn, Tfn))
cB = repeat(b0,Tf);


Gm = H_Bsp*Smsp
Go = H_Bsp*Sosp


Kym     = Gm'*Σsp_inv*Gm
CL = cholesky(Hermitian(Kym))
μ_y = CL.UP\(CL.PtL\Gm'*Σsp_inv)*(Xb*cB-Go*longyo)

# YYt[Sm_bit] = μ_y + CL.UP\randn(nm,)


M_zsp, z_vec, T_z = BEAVARs.CPZ_makeM_inter(z_tab,YYt,Sm_bit,datesHF,varNamesLF,fvarNames,(3,12,1),nm,Tf);
O_zsp = sparse(I,T_z,T_z)*0.0000002;
MOiM = M_zsp'*(O_zsp\M_zsp)
MOiz = M_zsp'*(O_zsp\z_vec)
KymBar = MOiM + Kym;


CLBar = cholesky(Hermitian(KymBar))
# C = CLBar.PtL; # Ct = CLBar.UP;
μ_yBar = (CLBar.UP)\(CLBar.PtL\(MOiz + Kym*μ_y))

YYt[Sm_bit] = μ_yBar +  CLBar.UP\randn(nm,)

plot(YY)


hypSetup=hypChan2020()
Y, X, T, n = mlag_r(YY,p)

Yt = vec(Y')

(deltaP, sigmaP, mu_prior) = trainPriors(YY,p);

(idx_kappa1,idx_kappa2, V_Minn, beta_Minn) = prior_Minn(n,p,sigmaP,hypSetup);
k = n*p+1;
Xsur = SUR_form(X,n)
XiSig = Xsur'*kron(sparse(Matrix(1.0I, T, T)),sparse(1:n,1:n,1.0./sigmaP));
K_β = sparse(1:n*k,1:n*k,1.0./V_Minn) + XiSig*Xsur;
cholK_β = cholesky(Hermitian(K_β));
beta_hat = cholK_β.UP\(cholK_β.PtL\(beta_Minn./V_Minn + XiSig * Yt))


# ndraws = nsave+nburn;
# store_beta=zeros(n^2*p+n,nsave)
beta = beta_hat + cholK_β.UP\randn(k*n);

k_nc = n*p; 
B_draw = reshape(beta,k_nc+1,n)'

b0[:] = B_draw[:,1]
blk[:,n+1:end] = B_draw[:,2:end]
# updating cB
BEAVARs.CPZ_update_cB!(cB,B_draw[:,2:end],b0,Y0,cB_b0_ind,p,n)

# updating H_B
H_Bsp.nzval[:] = -blk[blk_ind]