
"""
    mlag(Yfull::Matrix{Float64},p::Integer)
    Creates lags of a matrix for a VAR representation with a constant on the right of X
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
    const_loc = 0;          # means the constant is on the right of X, i.e. the bottom of beta in X*β   
    Y = @views Yfull[p+1:end,:]
    return Y, X, T, n, const_loc
end

"""
    mlagL!(YY,Y,X,p,n;intercept=1)
    Creates lagged matrices Y and X by mutating. Assumes constant is on the left of X
"""
function mlagL!(YY,Y,X,p,n)
    for i = 1:p
        X[:,1 + 1 +n*(i-1):n + (1) + n*(i-1)] = YY[p-i+1:end-i,:]
    end
    Y[:,:] = YY[p+1:end,:];
    return Y, X
end


"""
    mlagL(Yfull::Matrix{Float64},p::Integer)
    Creates lags of a matrix for a VAR representation with a constant on the left
        Yfull: a matrix of dimensions T+p x N returns a matrix Y with dimensions TxN and X with dimenions Tx(N*p+1)
"""
function mlagL(Yfull::Matrix{Float64},p::Integer)
    (Tf, n) = size(Yfull)
    T = Tf-p;
    X = zeros(T,n*p+1)
    X[:,1] = ones(T,1)
    for i = 1:p
        X[:,2+n*(i-1):1+n+n*(i-1)] = Yfull[p-i+1:end-i,:]
    end
    Y = Yfull[p+1:end,:]
    const_loc = 1;
    return Y, X, T, n, const_loc
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
    # ε = Y-X*β;
    σ_sq = ε'*ε/(size(Y,1)-size(X,2));
    return β, ε, σ_sq
end

"""
    ols(Y,X)

Performs standard linear regression on two matrices Y and X,
returning β as a vector, the vector of residuals ε and the variance σ_sq
"""
function ols2(Y,X,β,ε,T,n,σ_sq)
    # β = vec(X\Y);
    β[:,:] = X\Y;
    # ε = Y-X*reshape(β,size(X,2),n);
    ε[:,:] = Y-X*β;
    σ_sq[:,:] = mul!(σ_sq,ε',ε)./(T-n);
    return β, ε, σ_sq
end



@doc raw"""
    trainPriors(Z0::Matrix{Float64},p::Int64)

Independent AR(p) regressions with constant to estimate prior values for further Bayesian estimation

For a training sample `Z0` with `n` variables and `p` lags the function will do column-wise
`n` linear regressions of order p and return a matrix 

deltaP has the constant on the bottom and the lags (1) to (p) in rows [1:end-1,:]

"""
function trainPriors(YY::Matrix{Float64},p::Int64)
    Y, X, T, n, const_loc = mlag(YY,p)     # mlag takes constant from the right
    mu_prior = vec(mean(Y, dims=1));
    deltaP = zeros(p+1,n)
    sigmaP = vec(zeros(n,1))

    # Do univariate AR(p) linear regressions with constant,
    # Assumes that if there is a constant in X is on the right
    for ii = 1:n
        b,res,sig = ols(Y[:,ii],[X[:,ii:n:end-1] X[:,end]])
        deltaP[:,ii] = b
        sigmaP[ii,:] = sig
    end
    return deltaP, sigmaP, mu_prior
end




@doc raw"""
    trainPriors(Z0::Matrix{Float64},p::Int64)

Independent AR(p) regressions with constant to estimate prior values for further Bayesian estimation

For a training sample `Z0` with `n` variables and `p` lags the function will do column-wise
`n` linear regressions of order p and return a matrix 

deltaP has the constant on the bottom and the lags (1) to (p) in rows [1:end-1,:]

"""
function updatePriors!(Y,X,n::Int,mu_prior,deltaP,sigmaP,intercept)
    # Do univariate AR(p) linear regressions with constant,
    # if intercept = 1, constant is on the left of X (bottom of beta in X*beta)
    for ii = 1:n
        b,res,sig = ols(Y[:,ii],[ones(size(X,1),) X[:,intercept+ii:n:end-1+intercept]])
        deltaP[:,ii] = b
        sigmaP[ii,:] = sig
    end
    return deltaP, sigmaP, mu_prior
end


"""
    percentile_mat(A, p; dims) 

Calculates percentiles in a matrix across a specific dimension. Taken from https://github.com/JuliaStats/Statistics.jl/issues/23
and written by [https://github.com/holomorphism](holomorphism). I couldn't find such a function in Julia

"""
function percentile_mat(A, p; dims) 
    prctile_mat = mapslices(x->quantile(x, p), A; dims);
    prctileSlim_mat = dropdims(prctile_mat, dims = dims)
    return prctileSlim_mat
end


