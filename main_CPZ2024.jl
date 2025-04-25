include("devPkgs.jl")
using BEAVARs
using DelimitedFiles
using LinearAlgebra
#using FillArrays
using BlockArrays
using BlockBandedMatrices
using SparseArrays

p = 2;
n = 3;
T = 5;
k = n*(p+1); kn = k*3
Tf = n*T;
YY = rand(T,n);
YY[2,1] = NaN; YY[4,2] = NaN; YY[3,3] = NaN;

YYt = YY'
Sm_bit = isnan.(YYt)
So_bit = .!isnan.(YYt)
indC_nan_wide = findall(Sm_bit) #  Cartesian indices of missing values
indC_non_wide = findall(!isnan,YYt);  # Cartesian indices of not missing values
vYYt = vec(YYt);

# convert between linear and cartesian indices
indL_all = LinearIndices(YYt);
indL_nan_wide = indL_all[indC_nan_wide] # are the linear indices of NaN values
indL_non_wide = indL_all[indC_non_wide] # are the linear indices of non NaN values


S_full = I(T*n);
Sm = S_full[:,indL_nan_wide]
So = S_full[:,indL_non_wide]


# Go = HB So ; Go yo = HB So yo = HB YY[So_bit]
B0 = 1.0*I(3)
B1 = rand(3,3)
B2 = rand(3,3)*0.1

blk = [B0;B1;B2]
lin = LinearIndices(blk)


bbm = BlockBandedMatrix(ones(Float64,Tf,Tf),n*ones(Int,T,),n*ones(Int,T,),(p,0))


# create a matrix in the shape of H_B which hast the linear indices of the coefficient matrix. It will be used to select the coefficients
bbm_lin = BlockBandedMatrix(ones(Int,Tf,Tf),n*ones(Int,T,),n*ones(Int,T,),(p,0))
 
# this one iterates over the Blocks of the block-banded matrix bbm and fills them with the linear indices of the coefficients to generate the indices that select from the coefficient matrix 
for ij = 0:p
    for ii in 1:div(Tf,n)-ij
        bbm_lin[Block(ii+ij,ii)]=lin[ 1 + (ij-0)*n : n + (ij-0)*n ,:]
    end
end


blklin_ind = deepcopy(bbm_lin.data)


# Initialize sparse H_B
H_Bsp=sparse(bbm)
H_Bsp.nzval[:] = @view blk[blklin_ind];

