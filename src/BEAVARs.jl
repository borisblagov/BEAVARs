module BEAVARs
using LinearAlgebra, Distributions, SparseArrays

# from init_functions.jl
export mlag, mlag_r, ols, precentile_mat

# from Banbura2010
export makeDummiesMinn!, makeDummiesSumOfCoeff!, getBeta!, getSigma!, gibbs_beta_sigma,trainPriors, Banbura2010

export irf_chol

include("init_functions.jl")
include("Banbura2010.jl")

# functions for irfs
#includet("irfs.jl")
"""
    irf_chol(beta_vec,sigma_vec,n,p,intercept,n_irf)

    Calculates impulse responses for a set of parameters in a vector beta_vec and variance-covariance matrix in a vector sigma_vec
    Uses the Cholesky decomposition and the companion form.
        - beta_vec needs to have the parameters for each equation stacked on top of each other
        - intercept: 
            - the function needs to know whether the intercept is at the left or the right to drop it correctly
            - use 1 if it is on the right, 0 if there is none, and -1 if it is on the left
"""
function irf_chol(beta_vec,sigma_vec,n,p,intercept,n_irf)

    k_nc = n*p;
    IRF_mat = zeros(n_irf,n,n);
    for i_var = 1:n
        B_draw = reshape(beta_vec,k_nc+1,n)'
        if intercept == 1
            Bnc_draw = B_draw[:,1:k_nc]
        end
    
        Σ_draw = reshape(sigma_vec,n,n)
    
        F = [Bnc_draw; 1.0I(n*(p-1)) zeros(n*(p-1),n)]
    
        A_chol = cholesky(Σ_draw)
    
        responses = zeros(n,n_irf);
        shocks    = zeros(n,);     # shocks vector (either with unity shocks or std.dev shocks)
    
        shocks[i_var,] = 1.0
    
        responses[:,1] = A_chol.L*shocks;
        impulses = [responses[:,1]; zeros(n*(p-1))]
    
        for i_irf = 2:n_irf
            F_n = F^(i_irf-1)
            resp_all = F_n*impulses
            responses[:,i_irf] = resp_all[1:n]
        end
    
        IRF_mat[:,:,i_var] = transpose(responses)
    
    end
        
    return IRF_mat
    end




end

