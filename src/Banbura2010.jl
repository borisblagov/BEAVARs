@doc raw"""
# makeDummiesMinn!(sigma::Vector{Float64},delta,lambda,n::Integer,p::Integer,Y_d1))

Fills a matrix $Y_d$ and $X_d$ following eq. (5) in [Banbura, Giannone, Reichling (2010), JAE,
Large Bayesian Autoregressions](https://doi.org/10.1002/jae.1137). 

```math
Y_d =   \begin{bmatrix}
         diag(σ_1*δ_1 \dots σ_N*δ_n)/λ\\
         \mathbf{0}_{n*(p-1) \times n}\\
         diag(σ_1, \dots, σ_n)\\
         \mathbf{0}_{1,n \times n}
         \end{bmatrix}\\
```

```math
X_d = \begin{bmatrix}
        diag(1,\dots,p) \otimes diag(σ_1, \dots σ_n)./ λ \quad \mathbf{0}_{n*p,1}\\
        \mathbf{0}_{n,n*p+1}\\
        \mathbf{0}_{1,n*p}  ε\\
        \end{bmatrix}
```

Instad of creating the matrix every time, the function uses mutation on matrices called `Y_d1` and `X_d1`.

For `Y_d1` it is the the diagonal of the first `(1:n,1:n)` block and in the diagonal of the third block `n+n*(p-1)+1:n+n*(p-1)+n`

For `X_d1` it is the diagonal of the first (n * p x n * p), esesntially `X_d1` is populated along its diagonal and only the constant is added at the end

```lang=julia
function makeDummiesMinn!(sigma::Vector{Float64},delta::Vector{Float64},lambda,n::Integer,p::Integer,epsi,Y_d1,X_d1)
    CI1_Yd1 = CartesianIndex.(1:n,1:n)  # These are the diagonal indices for the top block 
    CI2_Yd1 = CartesianIndex.(n+n*(p-1)+1:n+n*(p-1)+n,1:n)  # These are the diagonal indices for the block in the middle
    CI_Xd1 = CartesianIndex.(1:n*p,1:n*p)
    Y_d1[CI1_Yd1] = sigma.*delta./lambda;
    Y_d1[CI2_Yd1] = sigma;
    X_d1[end,end] = epsi;
    X_d1[CI_Xd1]  = repeat(sigma./lambda,p).*repeat(1:p,inner=n)
    return Y_d1, X_d1;
end
```

"""
function makeDummiesMinn!(sigma::Vector{Float64},delta::Vector{Float64},lambda,n::Integer,p::Integer,epsi,Y_d1,X_d1)
    CI1_Yd1 = CartesianIndex.(1:n,1:n)  # These are the diagonal indices for the top block 
    CI2_Yd1 = CartesianIndex.(n+n*(p-1)+1:n+n*(p-1)+n,1:n)  # These are the diagonal indices for the block in the middle
    CI_Xd1 = CartesianIndex.(1:n*p,1:n*p)
    Y_d1[CI1_Yd1] = sigma.*delta./lambda;
    Y_d1[CI2_Yd1] = sigma;
    X_d1[end,end] = epsi;
    X_d1[CI_Xd1]  = repeat(sigma./lambda,p).*repeat(1:p,inner=n)
    return Y_d1, X_d1;
end

"""
    makeDummiesSumOfCoeff!()

TBW
"""
function makeDummiesSumOfCoeff!(delta,mu,tau,n,p,Yd_2,Xd_2)
    CI1_Yd2 = CartesianIndex.(1:n,1:n)
    Yd_2[CI1_Yd2] = delta.*mu./tau;
    #xd2=[kron(ones(1,p),yd2) zeros(n,1)];
    
    CI1_Xd2 = CartesianIndex.(repeat(1:n,p),1:n*p);
    Xd_2[CI1_Xd2] = repeat(delta.*mu./tau,p)
    return Yd_2, Xd_2
end


"""
    getBeta!(beta::Vector{Float64},Sigma_b::Hermitian)

    Draws from the multivariate normal with mean β and var-covar Σ_b and mutates β

"""
function getBeta!(beta_mean::Vector{Float64},beta::Vector{Float64},Sigma_b::Hermitian)
    d=MvNormal(beta_mean,Sigma_b);
    beta[:]=vec(rand(d,1));
end

function getSigma!(Tstar,P,Sigma_tilde)
    d=InverseWishart(Tstar,P);
    Sigma_tilde[:,:]=rand(d);
    return Sigma_tilde
end



"""
    gibbs_beta_sigma(Tstar,n,p,Ystar,Xstar,nburn,nsave)

    Implements the gibbs sampler for the Bayesian VAR using dummy variables 
     as in Banbura, Giannone, and Reichlin (2010), Large Bayesian Vectorautoregressions, _Journal of Applied Econometrics_
"""
function gibbs_beta_sigma(Tstar,n,p,Ystar,Xstar,nburn,nsave)
    k = n*(n*p+1)
    ndraws = nsave+nburn;
    store_beta=zeros(k,nsave);
    store_sigma=zeros(n*n,nsave);
    beta_mean = vec(Xstar\Ystar)
    beta = similar(beta_mean)
    iXX = inv(Xstar'*Xstar)

    Sigma_tilde = zeros(n,n); Sigma_tilde[CartesianIndex.(1:n,1:n)] = ones(n);

    for ii = 1:ndraws
        Sigma_b = Hermitian(kron(Sigma_tilde,iXX))
            
        getBeta!(beta_mean,beta,Sigma_b)

        res = Ystar-Xstar*reshape(beta,n*p+1,n);
        P   = res'*res

        getSigma!(Tstar,P,Sigma_tilde)

        if ii>nburn
            store_beta[:,ii-nburn] = beta;
            store_sigma[:,ii-nburn] = vec(Sigma_tilde);
        end
    end
    return store_beta, store_sigma
end



function Banbura2010(Z::Matrix{Float64};lags::Integer=1,lambda::Float64=0.1,epsi::Float64=0.001,nburn::Integer=1000,nsave::Integer=2000)
#function BVAR_dummies(Z::Matrix{Float64},lags::Integer,param_str::PriorSetup)
#    @unpack_PriorSetup param_str
    p = lags;
    deltaP_mat, sigmaP_vec, mu_prior = trainPriors(Z,1)
    delta = deltaP_mat[1,:];

    Y,X,T,n = mlag(Z,p);

    
    tau = 10*lambda;
    Tstar = T+n+n*(p-1)+n+1+n;

    Xstar = zeros(Tstar,n*p+1)
    Ystar = zeros(Tstar,n);                   # holds the data with the prior. This is Y_star in eq. (6) of Banbura, Giannone, and Reichlin (2010), JAE,


    Ystar[1:T,:] = Y;                                   # The top part of Ystar is the data
    Xstar[1:T,:] = X;
    Yd1          = @view Ystar[T+1:T+n+n*(p-1)+n+1,:];  # views into Ystar and defines part of it as Yd1 which will be mutated later in makeDummiesMinn!
    Yd2          = @view Ystar[T+n+n*(p-1)+n+2:Tstar,:]; 

    Xd1          = @view Xstar[T+1:T+n+n*(p-1)+n+1,:];
    Xd2          = @view Xstar[T+n+n*(p-1)+n+1+1:Tstar,:];

    Yd1, Xd1 = makeDummiesMinn!(sigmaP_vec,delta,lambda,n,p,epsi,Yd1,Xd1);
    Yd2, Xd2 = makeDummiesSumOfCoeff!(delta,mu_prior,tau,n,p,Yd2,Xd2);


    store_beta, store_sigma =  gibbs_beta_sigma(Tstar,n,p,Ystar,Xstar,nburn,nsave);

    # display("done")

    return store_beta, store_sigma;
end