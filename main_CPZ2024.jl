include("devPkgs.jl")
using BEAVARs
using DelimitedFiles
using LinearAlgebra

n = 3;
T = 12;
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