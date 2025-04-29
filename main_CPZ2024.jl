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
#data_de = data_ta_full[[:tpuCaldara, :pmiDE, :gdpDE, :invDE, :hicpDE, :euribor]];
#var_names = colnames(data_de)
#YY = values(data_de);


z_var = 3;
z_vec = randn(4,);
p = 2; # number of lags
n = 3; # number of vars
Tf = 10; # full time span (with initial conditions)
k = n*(p+1); kn = k*n
Tfn = n*Tf;
YY = rand(Tf,n);
YY[2,1] = NaN; YY[4,2] = NaN; YY[3,3] = NaN;
YY[:,z_var] = fill(NaN, Tf,1);

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
B0 = 1.0*I(n)
B1 = rand(n,n)
B2 = rand(n,n)*0.1


blk = [B0 B1 B2 B2*0.1 B2*0.01 B2*0.01]
# lin = LinearIndices(blk)

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


xx = rand(n,n);
Σt_inv  = xx'*xx;
H_Bsp, blk_ind = BEAVARs.makeBlkDiag(Tfn,n,p,blk);
Σsp_inv, Σt_ind = BEAVARs.makeBlkDiag(Tfn,n,0,Σt_inv);

b0 = zeros(n,)
cB = repeat(b0,Tf);
X = sparse(Matrix(1.0I, Tfn, Tfn))
Gm = H_Bsp*Smsp
Go = H_Bsp*Sosp
longyo = vYYt[vec(So_bit)];


Kym     = Gm'*Σsp_inv*Gm
CL = cholesky(Hermitian(Kym))
μ_y = CL.U\(CL.L\Gm'*Σsp_inv)*(X*cB-Go*longyo);

YYt[Sm_bit] = μ_y + CL.U\randn(nm,)

iter = CartesianIndices(YYt)
ym_ci = iter[Sm_bit] # a vector of y_m with cartesian indices of the missing values in YYt
z_ci = CartesianIndices((z_var:z_var,1:Tf))
z_Mind_vec = vec(sum(ym_ci.==z_ci,dims=2))

M_inter = zeros(size(z_vec,1),nm)
M_z = @views M_inter[:,z_Mind_vec.==1]

hfWeights = [1/3; 2/3; 3/3; 2/3; 1/3];
M_z[1,1:5]=hfWeights;
M_z[2,2:6]=hfWeights;
M_z[3,3:7]=hfWeights;
M_z[4,4:8]=hfWeights


M_inter