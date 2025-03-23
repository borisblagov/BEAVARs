
"""
    mlag(Yfull::Matrix{Float64},p::Integer)
    Creates lags of a matrix for a VAR representation with a constant on the right
        Yfull: a matrix of dimensions T+p x N returns a matrix Y with dimensions TxN and X with dimenions Tx(N*p+1)
"""
function mlag(Yfull::Matrix{Float64},p::Integer)
    (Tf, n) = size(Yfull)
    T = Tf-p;
    X = zeros(T,n*p+1)
    for i = 1:p
        X[:,1+n*(i-1):n+n*(i-1)] = Yfull[p-i+1:end-i,:]
    end
    X[:,end] = ones(T,1)
    Y = Yfull[p+1:end,:]
    return Y, X, T, n
end


"""
    mlag_r(Yfull::Matrix{Float64},p::Integer)
    Creates lags of a matrix for a VAR representation with a constant on the left
        Yfull: a matrix of dimensions T+p x N returns a matrix Y with dimensions TxN and X with dimenions Tx(N*p+1)
"""
function mlag_r(Yfull::Matrix{Float64},p::Integer)
    (Tf, n) = size(Yfull)
    T = Tf-p;
    X = zeros(T,n*p+1)
    X[:,1] = ones(T,1)
    for i = 1:p
        X[:,2+n*(i-1):1+n+n*(i-1)] = Yfull[p-i+1:end-i,:]
    end
    Y = Yfull[p+1:end,:]
    return Y, X, T, n
end

"""
    ols(Y,X)

Performs standard linear regression on two matrices Y and X,
returning β as a vector, the vector of residuals ε and the variance σ_sq
"""
function ols(Y,X)
    n = size(Y,2)
    β = vec(X\Y);
    ε = Y-X*reshape(β,size(X,2),n);
    σ_sq = ε'*ε/size(Y,1);
    return β, ε, σ_sq
end




@doc raw"""
    trainPriors(Z0::Matrix{Float64},p::Int64)

Independent AR(p) regressions with constant to estimate prior values for further Bayesian estimation

For a training sample `Z0` with `n` variables and `p` lags the function will do column-wise
`n` linear regressions of order p and return a matrix 

deltaP has the constant on the bottom and the lags (1) to (p) in rows [1:end-1,:]

"""
function trainPriors(Z0::Matrix{Float64},p::Int64)
    Y, X, T, n = mlag(Z0,p)     # mlag takes constant from the right
    mu_prior = vec(mean(Y, dims=1));
    deltaP = zeros(p+1,n)
    sigmaP = vec(zeros(n,1))

    # Do univariate AR(p) linear regressions with constant,
    # Assumes that if there is a constant in X is on the right
    for ii = 1:n
        b,res,sig = ols(Y[:,ii],[X[:,ii:n:end-1] ones(T,1)])
        deltaP[:,ii] = b
        sigmaP[ii,:] = sig
    end
    return deltaP, sigmaP, mu_prior
end


"""
    precentile_mat(A, p; dims) 

Calculates percentiles in a matrix across a specific dimension. Taken from https://github.com/JuliaStats/Statistics.jl/issues/23
and written by [https://github.com/holomorphism](holomorphism). I couldn't find such a function in Julia

"""
function precentile_mat(A, p; dims) 
    prctile_mat = mapslices(x->quantile(x, p), A; dims);
    prctileSlim_mat = dropdims(prctile_mat, dims = dims)
    return prctileSlim_mat
end