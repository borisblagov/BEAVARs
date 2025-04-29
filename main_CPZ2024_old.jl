
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
