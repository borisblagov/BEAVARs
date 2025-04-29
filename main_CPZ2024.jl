include("devPkgs.jl")
using BEAVARs
using DelimitedFiles
using LinearAlgebra
# using FillArrays
# using BlockArrays
# using BlockBandedMatrices
using TimeSeries
using SparseArrays


dataM_bg_full = readtimearray("data/dataM_BG.csv"; format="dd/mm/yyyy", delim=',')
dataQ_bg_full = readtimearray("data/dataQ_BG.csv"; format="dd/mm/yyyy", delim=',')
varNamesM_full = colnames(dataM_bg_full)
#data_de = data_ta_full[[:tpuCaldara, :pmiDE, :gdpDE, :invDE, :hicpDE, :euribor]];
#var_names = colnames(data_de)
#YY = values(data_de);

#adding second GDP for testing
# bgTest_tab = dataQ_bg_full[:gdpBG];
# bgTest_tab = rename(bgTest_tab,:gdp2)
# dataQ_bg_full = merge(dataQ_bg_full,bgTest_tab)

dataM_bg_raw = dataM_bg_full
dataQ_bg_raw = dataQ_bg_full


varNamesM_full = colnames(dataM_bg_full)
varNamesQ_full = colnames(dataQ_bg_full)

varNamesHF = [:survIndustryBG];
varNamesLF = [:gdpBG]
varOrder   = [:gdpBG,:survIndustryBG]

# select only the needed data and transform it if needed
dataM_bg_tab = dataM_bg_raw[varNamesHF]
dataQ_bg_tab = percentchange(dataQ_bg_raw[varNamesLF])


# dataLF_tab = dataQ_bg_tab[1:4];
# dataHF_tab = dataM_bg_tab[1:16];
dataLF_tab = dataQ_bg_tab;
dataHF_tab = dataM_bg_tab;


# beyond this point we shouldn't need any more input from the user, i.e. we switch from Q and M (if you have monthly and quarterly data) or A and Q (if you have annual and quarterly) to HF and LF

# create the final z_tab
z_tab = dataLF_tab;
# add the z_tab as NaN values in the high-frequency tab
fdataHF_tab = merge(dataHF_tab,map((timestamp, values) -> (timestamp, values.*NaN), dataLF_tab[varNamesLF]),method=:outer)
fdataHF_tab = fdataHF_tab[varOrder]              # ordering the variables as the user wants them
fvarNames = colnames(fdataHF_tab)                # full list of the variable names


YY = values(fdataHF_tab)
datesHF = timestamp(fdataHF_tab)
datesLF = timestamp(dataLF_tab)
freqL_date = Month(datesLF[2])-Month(datesLF[1])
freqH_date = Month(datesHF[2])-Month(datesHF[1])

# tuple showing the specification: 1, 3, 12 are monthly quarterly, annually and 0,1 is growth rates or log-levels
freq_mix_tp = (convert(Int,freqH_date/Month(1)), convert(Int,freqL_date/Month(1)),0) # tuple with the high and low frequencies. 1 is monthly, 3 is quarterly, 12 is annually

p = 5; # number of lags
# n = 3; # number of vars
(Tf,n) = size(YY); # full time span (with initial conditions)
k = n*(p+1); kn = k*n
Tfn = n*Tf;
# YY = rand(Tf,n);
# YY[2,1] = NaN; YY[4,1] = NaN; YY[3,3] = NaN;
# YY[:,z_var] = fill(NaN, Tf,1);

YYt = (YY')
Sm_bit = isnan.(YYt)
So_bit = .!isnan.(YYt)
indC_nan_wide = findall(Sm_bit) #  Cartesian indices of missing values
indC_non_wide = findall(!isnan,YYt);  # Cartesian indices of not missing values
vYYt = vec(YYt);

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

# Go = HB So ; Go yo = HB So yo = HB YY[So_bit]
B0 = -1.0*I(n)
B1 = 1.0*I(n); B1 = [0.26 0; 0 0.96]
B2 = zeros(size(B1))


blk = [B0 B1 B2 B2*0.1 B2*0.01 B2*0.01]


Σt_inv  = 100*Matrix(I,n,n); Σt_inv = inv([0.000143 0; 0 0.000473])
H_Bsp, blk_ind = BEAVARs.makeBlkDiag(Tfn,n,p,-blk);
Σsp_inv, Σt_ind = BEAVARs.makeBlkDiag(Tfn,n,0,Σt_inv);

b0 = zeros(n,)
cB = repeat(b0,Tf);
X = sparse(Matrix(1.0I, Tfn, Tfn))
Gm = H_Bsp*Smsp
Go = H_Bsp*Sosp
longyo = vYYt[vec(So_bit)];


Kym     = Gm'*Σsp_inv*Gm
CL = cholesky(Hermitian(Kym))
μ_y = CL.U\(CL.L\Gm'*Σsp_inv)*(X*cB-Go*longyo)

# YYt[Sm_bit] = μ_y + CL.U\randn(nm,)


M_zsp, z_vec, T_z = BEAVARs.makeMinter(z_tab,YYt,Sm_bit,datesHF,varNamesLF,fvarNames,freq_mix_tp,nm,Tf);
O_zsp = sparse(I,T_z,T_z)*0.000002;
MOiM = M_zsp'*(O_zsp\M_zsp)
MOiz = M_zsp'*(O_zsp\z_vec)
KymBar = MOiM + Kym;

CLBar = cholesky(Hermitian(KymBar))
μ_yBar = CLBar.U\(CLBar.L\(MOiz + Kym*μ_y))
YYt[Sm_bit] = μ_yBar +  CLBar.U\randn(nm,)






# Old code# lin = LinearIndices(blk)

# bbm = BlockBandedMatrix(ones(Float64,Tfn,Tfn),n*ones(Int,Tf,),n*ones(Int,Tf,),(p,0))


# # create a matrix in the shape of H_B which hast the linear indices of the coefficient matrix. It will be used to select the coefficients
# bbm_lin = BlockBandedMatrix(ones(Int,Tfn,Tfn),n*ones(Int,Tf,),n*ones(Int,Tf,),(p,0))
#
# # this one iterates over the Blocks of the block-banded matrix bbm and fills them with the linear indices of the coefficients to generate the indices that select from the coefficient matrix 
# for ij = 0:p
#     for ii in 1:div(Tfn,n)-ij
#         bbm_lin[Block(ii+ij,ii)]=lin[ 1 + (ij-0)*n : n + (ij-0)*n ,:]
#     end
# end
# blklin_ind = deepcopy(bbm_lin.data)

# # Initialize sparse H_B
# H_Bsp=sparse(bbm)
# H_Bsp.nzval[:] = @view blk[blklin_ind];
# global Σ_lin = zeros(Float64,Tfn,Tfn)


# for ij = 0:p
#     for ii = 1:Tf-ij
#         Σ_lin[ ij*n + (ii-1)*n + 1 : n + (ii-1)*n +  ij*n, (ii-1)*n + 1 : n + (ii-1)*n] = ones(Float64,n,n)
#     end
# end
# Σ_lin_sp = BEAVARs.makeBlkDiag(Tfn,n,0,type=Float64)
# global Σ_linInt_sp = BEAVARs.makeBlkDiag(Tfn,n,0,type=Int)

# Σt_lin = LinearIndices(Σt)


# for ij = 0:0
#     for ii = 1:div(Tfn,n)-ij
#         Σ_linInt_sp[ ij*n + (ii-1)*n + 1 : n + (ii-1)*n +  ij*n, (ii-1)*n + 1 : n + (ii-1)*n] = Σt_lin
#     end
# end
# lin_ind = deepcopy(Σ_linInt_sp.nzval)


# T_z, n_z = size(z_tab);    # number of z vars
# z_mat = values(z_tab)
# global M_z = zeros(T_z*n_z,nm)
# global z_vec = zeros(T_z*n_z,)
# for ii_z = 1:n_z # iterator going through each variable in z_tab (along the columns)
#     datesLF_ii = timestamp(z_tab[varNamesLF[ii_z]])
#     iter = CartesianIndices(YYt)
#     ym_ci = iter[Sm_bit] # a vector of y_m with cartesian indices of the missing values in YYt
#     z_ci = CartesianIndices((z_var_pos[ii_z]:z_var_pos[ii_z],1:Tf))
#     z_Mind_vec_ii = vec(sum(ym_ci.==z_ci,dims=2)) # alternative z_Mind_vec=vec(indexin(ym_ci,z_ci)).!==nothing

#     M_inter_ii = zeros(T_z,nm)
#     M_z_ii = @views M_inter_ii[:,z_Mind_vec_ii.==1]

#     if size(datesHF,1)!==size(M_z_ii,2)
#         error("The size of M does not match the number of dates available in z_tab. Maybe the low-frequency data is longer? The problem is with variable number ", z_var_pos[ii_z])
#     end

#     # we need to watch out with the dates due to how the intertemporal constraint works Take for example growth rates Q and M
#     # y_t = 1/3 y_t - 2/3 y_{t-1} \dots - - 2/3 y_{t-3} - 1/3 y_{t-5}
#     # Intuitively, Q1 quarterly GDP (e.g. 01.01.2000) is the weighted sum of the monthly March, February, January, December, November, and October
#     # if y_t^Q is 01.01.2000, we need +2 and -2 months for the weights
#     if freq_mix_tp==(1,3,0)
#         hfWeights = [1/3; 2/3; 3/3; 2/3; 1/3]; n_hfw = size(hfWeights,1); #number of weights, depends on the variable transformation and frequency
#     elseif freq_mix_tp==(3,12,0)
#         # quarterly and yearly data with growth rates
#         hfWeights = [1/4; 2/4; 3/4; 1; 3/4; 2/4; 1/4]; n_hfw = size(hfWeights,1); #number of weights, depends on the variable transformation and frequency
#     else
#         error("This combination of frequencies and transformation has not been implemented")
#     end

#     for ii_zi in eachindex(datesLF_ii) # iterator going through each time point in datesHF
#         ii_M = findall(datesHF.==datesLF_ii[ii_zi])[1]       # find the low-frequency index that corresponds to the high-frequency missing value
#         # M_z_ii[ii_zi, findall(datesHF.==datesLF_ii[ii_zi])[1]-n_hfw+1:findall(datesHF.==datesLF_ii[ii_zi])[1]] = hfWeights # if shifted above
#         M_z_ii[ii_zi,ii_M-div((n_hfw-1),2): ii_M+div((n_hfw-1),2)]=hfWeights; # +2 and - 2 months for the weights or +3 and -3
#     end
#     M_z[(ii_z-1)*T_z + 1:T_z + (ii_z-1)*T_z,:] = M_inter_ii;
#     z_vec[(ii_z-1)*T_z + 1:T_z + (ii_z-1)*T_z,]  = values(z_tab[varNamesLF[ii_z]]);
# end
