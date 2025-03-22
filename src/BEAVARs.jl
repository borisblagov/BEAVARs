module BEAVARs
export mlag

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

end
