module BEAVARs
using   Revise, 
        LinearAlgebra, 
        Distributions, 
        SparseArrays,
        TimeSeries, 
        Parameters

# from init_functions.jl
export mlag, mlag_r, ols, percentile_mat, mlag_linked!

# from Banbura2010
export makeDummiesMinn!, makeDummiesSumOfCoeff!, getBeta!, getSigma!, gibbs_beta_sigma,trainPriors, Banbura2010, hypBanbura2010

# from irfs.jl
export irf_chol, irf_chol_overDraws, irf_chol_overDraws_csv

# from Chan_2020
export prior_Minn, SUR_form, Chan2020_LBA_Minn, draw_h_csv!, Chan2020_LBA_csv, prior_NonConj, draw_h_csv_opt!
export hypChan2020, Chan2020_LBA_csv_keywords, fcastChan2020_LBA_csv

export makeSetup, fortschr!, beavar, dispatchModel, makeOutput

# Structures, to be uncommented later
export modelSetup, Chan2020_LBA_csv_type, Chan2020_LBA_Minn_type, modelHypSetup, hypDefault_strct, outChan2020_LBA_csv, VARModelType


function fortschr!(Yfor,p,A)
    for i_for = 1:8
        tclass = @views vec(reverse(Yfor[1+i_for-1:p+i_for-1,:],dims=1)')
        tclass = [1;tclass];
        Yfor[p+i_for,:]=tclass'*A;
    end
end


# Structures for multiple dispatch across models

# types for models
abstract type VARModelType end
struct Chan2020_LBA_Minn_type <: VARModelType end
struct Chan2020_LBA_csv_type <: VARModelType end
struct Banbura2010_type <: VARModelType end


# types for parameters (VAR setup parameters, hyperparameters, and output storage)
abstract type modelSetup end
abstract type modelHypSetup end

# empty structure for initialising the model
struct hypDefault_strct <: modelHypSetup end

# structure initializing the VAR
@with_kw struct VARSetup <: modelSetup
    n::Int          # number of variables (will be overwritten)
    p::Int          # number of lags
    nsave::Int      # gibbs to save
    nburn::Int      # gibbs to burn
    n_irf::Int      # number of impulse responses
    n_fcst::Int     # number of forecast periods
    const_loc::Int  # location of the constant
end



# Constructor
@doc raw"""
# makeSetup(YY::Array{Float64},model=model_str; p::Int,nburn::Int,nsave::Int,n_irf::Int,n_fcst::Int)

Populates the constructor VARSetup with default and/or custom values.

Mandatory arguments (write directly the argument):

    YY:     your data in a matrix form
    model:  the model string you want to use. Currently supported are "Banbura2010", "Chan2020_LBA_Minn", "Chan2020_LBA_csv".

Keyword arguments (e.g. p = 4). They don't have to be specified.

    p:      number of lags, default is 4
    nburn:  number of burn-in draws that will be discarded, default is 2000
    nsave:  number of retained draws (total is then nburn + nsave), default is 1000
    n_irf:  horizon of impulse responses, default is 16
    n_fcst: horizon of forecasting periods, default is 8

Outputs

    VARSetup: the setup structure for the BEAVARs
    hypSetup: the setup structure for hyperparameters with default values. If you want to modify these, disregard the output and do it outside
    
Examples

```lang=julia
julia> YY = randn(10,3); makeSetup(YY::Array{Float64},"Banbura2010")
(BEAVARs.VARSetup
  n: Int64 3
  p: Int64 4
  nsave: Int64 1000
  nburn: Int64 1000
  n_irf: Int64 16
  n_fcst: Int64 8
  const_loc: Int64 0
, hypBanbura2010
  lambda: Float64 0.1
  epsi: Float64 0.001
)
```

```lang=julia
julia> YY = randn(10,5); makeSetup(YY::Array{Float64},"Banbura2010",p=1,nburn=500,nsave=100,n_irf=4,n_fcst=12)
(BEAVARs.VARSetup
  n: Int64 5
  p: Int64 1
  nsave: Int64 100
  nburn: Int64 500
  n_irf: Int64 4
  n_fcst: Int64 12
  const_loc: Int64 0
, hypBanbura2010
  lambda: Float64 0.1
  epsi: Float64 0.001
)
```
"""
function makeSetup(YY::Array{Float64},model_str::String,p::Int,n_irf::Int,n_fcst::Int,nburn::Int,nsave::Int)
    T,n = size(YY);
    if model_str == "Chan2020_LBA_csv"
        # hypSetup = hypChan2020()
        const_loc = 1;
        model_type = Chan2020_LBA_csv_type()
    elseif model_str == "Chan2020_LBA_Minn"
        # hypSetup = hypChan2020()
        const_loc = 1;
        model_type = Chan2020_LBA_Minn_type()
    elseif model_str == "Banbura2010"
        # hypSetup = hypBanbura2010()
        const_loc = 0;
        model_type = Banbura2010_type()
    else
        error("Model not found, make sure the spelling is completely correct, upper and lowercase matters!\n Possible models are: \n    Banbura2010 \n    Chan2020_LBA_Minn\n    Chan2020_LBA_csv\n")
    end
    return VARSetup(n,p,nsave,nburn,n_irf,n_fcst,const_loc), model_type
end


function makeHypSetup(::Chan2020_LBA_csv_type)
    return hypChan2020()
end
function makeHypSetup(::Chan2020_LBA_Minn_type)
    return hypChan2020()
end
function makeHypSetup(::Banbura2010_type)
    return hypBanbura2010()
end


function beavar(YY::Array{Float64},model_str=model_name::String;p::Int=4,nburn::Int=1000,nsave::Int=1000,n_irf::Int=16,n_fcst::Int = 8,hyp::modelHypSetup=hypDefault_strct())
    setup_str, model_type = makeSetup(YY,model_str,p,n_irf,n_fcst,nburn,nsave);
    
    # checking if user supplied the hyperparameter structure
    if isa(hyp,hypDefault_strct)
        # if not supplied, make a default one
        hyp_strct = makeHypSetup(model_type)
        # println("using the default hyperparameters")
    else        
        # else use supplied
        hyp_strct = hyp
        # println("using the supplied parameters")
    end
    dispatchModel(YY,model_type,setup_str, hyp_strct);
    
    return setup_str, hyp_strct
end

function dispatchModel(YY,::Chan2020_LBA_Minn_type,setup_str, hyper_str)
    # println("Hello Minn")
    store_beta, store_sigma = Chan2020_LBA_Minn(YY,setup_str,hyper_str);
    return store_beta, store_sigma
end


function dispatchModel(YY,::Chan2020_LBA_csv_type,setup_str, hyper_str)
    println("Hello csv")
    store_beta, store_h, store_Σ, s2_h_store, store_ρ, store_σ_h2, eh_store = Chan2020_LBA_csv(YY,setup_str,hyper_str);
    return store_beta, store_h, store_Σ, s2_h_store, store_ρ, store_σ_h2, eh_store
end

function dispatchModel(YY,::Banbura2010_type,setup_str, hyper_str)
    println("Hello Banbura2010")
    store_beta, store_sigma = Banbura2010(YY,setup_str,hyper_str);
end

#----------------------------------------
# Chan 2020 LBA functions

@with_kw struct hypChan2020 <: modelHypSetup
    c1::Float64     = 0.04; # hyperparameter on own lags
    c2::Float64     = 0.01; # hyperparameter on other lags
    c3::Float64     = 100;  # hyperparameter on the constant
    ρ::Float64      = 0.8;
    σ_h2::Float64   = 0.1;
    v_h0::Float64    = 5.0; 
    S_h0::Float64    = 0.01*(v_h0-1.0); 
    ρ_0::Float64     = .9; 
    V_ρ::Float64     = 0.04;
    q::Float64       = 0.5;
end


@doc raw"""
# Chan2020_LBA_Minn(YY,VARSetup,hypSetup)

Implements the classic homoscedastic Minnesota prior with a SUR form following Chan (2020)

"""
function Chan2020_LBA_Minn(YY,VARSetup::modelSetup,hypSetup::modelHypSetup)
    @unpack p,nburn,nsave = VARSetup
    Y, X, T, n = mlag_r(YY,p)

    Yt = vec(Y')

    (deltaP, sigmaP, mu_prior) = trainPriors(YY,4);

    (idx_kappa1,idx_kappa2, V_Minn, beta_Minn) = prior_Minn(n,p,sigmaP,hypSetup);
    k = n*p+1;
    Xsur = SUR_form(X,n)
    XiSig = Xsur'*kron(sparse(Matrix(1.0I, T, T)),sparse(1:n,1:n,1.0./sigmaP));
    K_β = sparse(1:n*k,1:n*k,1.0./V_Minn) + XiSig*Xsur;
    cholK_β = cholesky(Hermitian(K_β));
    beta_hat = cholK_β.UP\(cholK_β.PtL\(beta_Minn./V_Minn + XiSig * Yt))
    # CSig = sparse(1:n,1:n,sqrt.(sigmaP));
    
    ndraws = nsave+nburn;
    store_beta=zeros(n^2*p+n,nsave)
    for ii = 1:ndraws 
        beta = beta_hat + cholK_β.UP\randn(k*n);
        if ii>nburn
            store_beta[:,ii-nburn] = beta;
        end
    end

    Σ = Matrix(1.0I, n, n); Σ[diagind(Σ)]=diag(Σ).*sigmaP;
    store_sigma = repeat(vec(Σ),1,nsave);
    return store_beta, store_sigma
end


@doc raw"""
    (idx_kappa1,idx_kappa2, C, beta_Minn) = prior_Minn(n,p,sigmaP,hypSetup)

Impements a Minnesota prior with scaling for the off-diagonal elements as in Chan, J.C.C. (2020). Large Bayesian Vector Autoregressions. In: P. Fuleky (Eds),
Macroeconomic Forecasting in the Era of Big Data, 95-125, Springer, Cham. If you remove the hyperparameters ``c_1``, ``c_2``, ``c_3`` (or set them equal to 1),
you would end up with the getC function from Chan, J.C.C. (2021). Minnesota-Type Adaptive Hierarchical Priors for 
Large Bayesian VARs, International Journal of Forecasting, 

Outputs: 
- C is the variance of the Minnesota prior defined as
```math
C_{n^2p+n \times 1} = 
\begin{cases}
        \dfrac{c_1}{l^2}, & \text{for the coefficient on the $l$-th lag of variable $i$}\\
        \dfrac{c_2 s^2_i}{l^2 s^2_j}, & \text{for the coefficient on the $l$-th lag of variable $j$, $j\neq i$}\\
        c_3, & \text{for the intercept}\\
\end{cases}
```
It is commonly known as `V_Minn` or `V_theta` in other codes.

The hyperparameters have default values of ``c_1 = 0.04``; ``c_2 = 0.01``; ``c_3 = 100``. In the Chan (2021) IJF paper they
are referred to as ``\kappa_1``, ``\kappa_2``, and ``\kappa_4`` (``\kappa_3`` is used there for prior on contemporaneous relationships)

"""
function prior_Minn(n::Integer,p::Integer,sigmaP_vec::Vector{Float64},hypSetup)    
    @unpack c1,c2,c3 = hypSetup
    C = zeros(n^2*p+n,);
    beta_Minn = zeros(n^2*p+n);
    np1 = n*p+1 # number of parameters per equation
    idx_kappa1 =  Vector{Int}()
    idx_kappa2 =  Vector{Int}()
    idx_count = 1    
    Ci = zeros(np1,1)     # vector for equation i

    for ii = 1:n
      for j = 1:n*p+1       # for j=1:n*p+1 
        l = ceil((j-1)/n)       # Here we need a float, as afterwards we will divide by l 
        idx = mod(j-1,n);       # this will count if its own lag, non-own lag, or constant
        if idx==0
            idx = n;
        end
        if j == 1
            Ci[j] = c3;
        elseif idx == ii
            Ci[j] = c1/l^2;
            push!(idx_kappa1,idx_count)
        else
            Ci[j] = c2*sigmaP_vec[ii]/(l^2*sigmaP_vec[idx]);
            push!(idx_kappa2,idx_count)
        end
        idx_count += 1
    end

    C[(ii-1)*np1+1:ii*np1] = Ci

    end
    return idx_kappa1,idx_kappa2, C, beta_Minn

end


@doc raw"""
    (idx_kappa1,idx_kappa2, C, beta_Minn) = prior_NonConj(n,p,sigmaP,hyp)

Impements a Minnesota prior for a non-conjugate case as in Chan, J.C.C. (2020). Large Bayesian Vector Autoregressions. In: P. Fuleky (Eds),
Macroeconomic Forecasting in the Era of Big Data, 95-125, Springer, Cham.

Outputs: 
- C is the variance of the Minnesota prior defined as
```math
C_{n*p+1 \times 1} = 
\begin{cases}
        \dfrac{c_1}{l^2 s^2_i}, & \text{for the coefficient on the $l$-th lag of variable $i$}\\
        c_3, & \text{for the intercept}\\
\end{cases}
```

Note that C is now (n*p+1 x 1) and not n^2*p+n as [`prior_Minn(x)`](@ref)

"""
function prior_NonConj(n::Integer,p::Integer,sigmaP_vec::Vector{Float64},hypSetup)    
    @unpack c1,c3 = hypSetup;
    C = zeros(n^2*p+n,);
    beta_Minn = zeros(n^2*p+n);
    np1 = n*p+1 # number of parameters per equation
    idx_kappa1 =  Vector{Int}()
    idx_kappa2 =  Vector{Int}()
    idx_count = 1    
    C = zeros(np1,)     # vector for equation i


    for j = 1:n*p+1       # for j=1:n*p+1 
    l = ceil((j-1)/n)       # Here we need a float, as afterwards we will divide by l 
    idx = mod(j-1,n);       # this will count if its own lag, non-own lag, or constant
    if idx==0
        idx = n;
    end
    if j == 1
        C[j] = c3;
    else
        C[j] = c1/(l^2*sigmaP_vec[idx]);
        push!(idx_kappa1,idx_count)
    end

    end
    return idx_kappa1,idx_kappa2, C, beta_Minn

end

@doc raw"""
    draw_h_csv!(h,s2_h,ρ,σ_h2,n,H_ρ)

Draws the log-volatilities h by mutating the array from the previous draw using an MH algorithm

"""
function draw_h_csv!(h,s2_h,ρ,σ_h2,n,H_ρ)

    accept = false;
    T_h = size(s2_h,1);
    # H_ρ = sparse(Matrix(1.0I, T_h, T_h)) - sparse(ρ*diagm(-1=>repeat(1:1,T_h-1)));
    H_inv = H_ρ'*sparse(1:T_h,1:T_h,[(1-ρ^2)/σ_h2; 1/σ_h2*ones(T_h-1,)])*H_ρ
    ϵ_h = 1.0;
    ht = similar(h); ht[:,] = h[:,];
    
    local(K_h)
    while ϵ_h > 10^(-3)
        eht = exp.(ht);
        sieht = s2_h[:,]./eht
        f_h = -n/2.0 .+ 0.5*sieht;
        G_h = 0.5*sieht
        K_h = H_inv + sparse(1:T_h,1:T_h,G_h);
        new_ht = K_h\(f_h + G_h.*ht)
        ϵ_h = maximum(abs.(new_ht-ht))
        ht[:,] = new_ht[:,]
    end
    C_K_h = cholesky(K_h)
    
    hstar = ht
    logc = -0.5*hstar'*H_inv*hstar - (n/2.0) *sum(hstar) .- 0.5*exp.(-hstar)'*s2_h[:,] .+ log(3)
    
    flag = false
    local hc = similar(h)
    local alpARc = 0.0;
    while flag == 0
        hc = ht + C_K_h.U\randn(T_h,);        
        alpARc = -0.5*hc'*H_inv*hc .- (n/2.0)*sum(hc) .- 0.5*exp.(-hc)'*s2_h[:,] .+ 0.5*(hc-ht)'*K_h*(hc-ht) .- logc;
        
        if alpARc > log(rand())
            flag = true;
        end
    end        
    alpAR = -0.5*h'*H_inv*h .- (n/2.0)*sum(h) .- 0.5*exp.(-h)'*s2_h[:,] .+ 0.5*(h-ht)'*K_h*(h-ht) .- logc;
    if alpAR < 0
        alpMH = 1.0;
    elseif alpARc < 0
        alpMH = - alpAR;
    else
        alpMH = alpARc - alpAR;
    end    
    if alpMH > log(rand())
        h[:,] = hc[:,]; 
        accept = true;
    end

    
    return h, accept
end # end draw_h_csv!





function draw_h_csv_opt!(h::Vector{Float64},s2_h,ρ,σ_h2,n,H_ρ)

    accept = false;
    T_h = size(s2_h,1);
    # H_ρ = sparse(Matrix(1.0I, T_h, T_h)) - sparse(ρ*diagm(-1=>repeat(1:1,T_h-1)));
    H_inv = H_ρ'*sparse(1:T_h,1:T_h,[(1.0-ρ^2)/σ_h2; 1.0/σ_h2*ones(T_h-1,)])*H_ρ
    ϵ_h = 1.0;
    # ht = similar(h); 
    # ht[:,] = h[:,];
    ht = h;
    nd2 = n/2.0;

    local K_h
    while ϵ_h > 10^(-3)
    #     eht = exp.(ht);
    #     sieht = s2_h[:,]./eht
        # f_h = -nd2 .+ 0.5*sieht;
        # G_h = 0.5*sieht
        G_h = 0.5.*s2_h[:,]./exp.(ht)
        K_h = H_inv + sparse(1:T_h,1:T_h,G_h);
        new_ht = K_h\(-nd2 .+ G_h + G_h.*ht)
        ϵ_h = maximum(abs.(new_ht-ht))
        ht[:,] = new_ht[:,]
    end
    C_K_h = cholesky(K_h)
    
    hstar = ht
    logc = -0.5*hstar'*H_inv*hstar - nd2 *sum(hstar) .- 0.5*exp.(-hstar)'*s2_h .+ log(3)
    
    flag = false
    local hc = similar(hstar)
    # f(x,y) = -0.5*x'*H_inv*x .- nd2*sum(x) .- 0.5*exp.(-x)'*s2_h .+ 0.5*(x-y)'*K_h*(x-y) .- logc;
    alpARc = 0.0;
    while flag == 0
        hc = ht + C_K_h.U\randn(T_h,);   
        
        alpARc = -0.5*hc'*H_inv*hc .- nd2*sum(hc) .- 0.5*exp.(-hc)'*s2_h .+ 0.5*(hc-ht)'*K_h*(hc-ht) .- logc;
        # alpARc = f(hc,ht);
        if alpARc > log(rand())
            flag = true;
        end
    end        
        alpAR = -0.5*h'*H_inv*h .- nd2*sum(h) .- 0.5*exp.(-h)'*s2_h .+ 0.5*(h-ht)'*K_h*(h-ht) .- logc;
        # alpAR = f(h,ht);
    if alpAR < 0.0
        alpMH = 1.0;
    elseif alpARc < 0.0
        alpMH = - alpAR;
    else
        alpMH = alpARc - alpAR;
    end    
    if alpMH > log(rand())
        h[:,] = hc[:,]; 
        accept = true;
    end

    
    return h, accept
end # end draw_h_csv!


function SUR_form(X,n)
    repX = kron(X,ones(n,1));
    r,c = size(X);
    idi = repeat(1:r*n,inner=c);
    idj=repeat(1:c*n,r);
    Xout = sparse(idi,idj,vec(repX'));

    return Xout
end


@doc raw"""
# Chan2020_LBA_csv(YY,VARSetup,hypSetup)

Implements the BVAR with Minnesota prior with a SUR form and common stochastic volatilty (csv) following Chan (2020)

"""
function Chan2020_LBA_csv(YY::Array{Float64},VARSetup::modelSetup,hypSetup::modelHypSetup)
    @unpack ρ, σ_h2, v_h0, S_h0, ρ_0, V_ρ = hypSetup
    @unpack p, nsave, nburn = VARSetup

    Y, Z, T, n = mlag_r(YY,p)
    (deltaP, sigmaP, mu_prior) = trainPriors(YY,4)  # the 4 here is hard coded so that you can compare models with the same initial conditions
    np1 = n*p+1; # number of parameters per equation
    # VARsetup.n = 20;
    # print((n))
    
    (idx_kappa1,idx_kappa2, V_Minn, beta_Minn) = prior_NonConj(n,p,sigmaP,hypSetup);
    
    A_0 = reshape(beta_Minn,np1,n);
    V_Ainv = sparse(1:np1,1:np1,1.0./V_Minn)
    
    S_0 = Diagonal(sigmaP);
    Σ = S_0;
    v_0 = n+3;
    h = zeros(T,)
    H_ρ = sparse(Matrix(1.0I, T, T)) - sparse(ρ*diagm(-1=>repeat(1:1,T-1)));
    # dv = ones(T,); ev = -ρ.*ones(T-1,);
    # H_ρ = Bidiagonal(dv,ev,:L)
    
    # This part follows page 19 of Chan, J. (2020)
    ndraws = nsave+nburn;
    store_beta = zeros(np1*n,nsave);
    store_h = zeros(T,nsave);
    store_Σ = zeros(n,n,nsave);
    store_s2_h = zeros(T,nsave);
    ρ_store = zeros(nsave,);
    σ_h2_store = zeros(nsave,); 
    eh_store = zeros(T,nsave);
    
    eh = similar(h);

    Ωinv = sparse(1:T,1:T,exp.(-h));
    dg_ind_Ωinv = diagind(Ωinv);
    for ii = 1:ndraws 
        Ωinv[dg_ind_Ωinv]=exp.(-h);
        ZtΩinv = Z'*Ωinv;
        
        K_A = V_Ainv + ZtΩinv*Z;
        A_hat = K_A\(V_Ainv\A_0 + ZtΩinv*Y);
        S_hat = S_0 + A_0'*V_Ainv*A_0 + Y'*Ωinv*Y - A_hat'*K_A*A_hat;
        S_hat = (S_hat+S_hat')/2;
    
        Σ = rand(InverseWishart(v_0+T,S_hat));
        cholΣ  = cholesky(Σ).U; # if we get the upper we don't need constant transpose
        A = A_hat + (cholesky(Hermitian(K_A)).U\randn(np1,n))*cholΣ;
    
        # Errors
        U = Y - Z*A
        s2_h = sum((U/cholΣ).^2,dims=2)
    
        draw_h_csv!(h,s2_h,ρ,σ_h2,n,H_ρ)
        eh[1,] = h[1]*sqrt(1-ρ^2);
        eh[2:end,] = h[2:end,].-ρ.*h[1:end-1,]
        σ_h2 = 1.0./rand(Gamma(v_h0+T/2.0, 1/( S_h0 + sum(eh.^2)./2.0 ) ) )
    
        K_rho = 1.0/V_ρ + sum(h[1:T-1,].^2)./σ_h2;
        ρ_hat = K_rho\(ρ_0/V_ρ + h[1:T-1,]'*h[2:T,]./σ_h2);
        ρ_c = ρ_hat + sqrt(K_rho)\randn();
    
        g_ρ(x) = -0.5*log(σ_h2./(1-x.^2))-0.5*(1.0-x.^2)/σ_h2*h[1]^2;
        if abs(ρ_c)<.999
            alpMH = exp(g_ρ(ρ_c)-g_ρ(ρ));
            if alpMH>rand()
                ρ = ρ_c;            
            end
        end    
        H_ρ[diagind(H_ρ,-1)]=fill(-ρ,T-1);
        # H_ρ.ev.=-ρ.*ones(T-1,)

        if ii>nburn
            store_beta[:,ii-nburn] = vec(A);
            store_h[:,ii-nburn] = h;
            store_Σ[:,:,ii-nburn] = Σ;
            store_s2_h[:,ii-nburn] = s2_h;
            ρ_store[ii-nburn,] = ρ;
            σ_h2_store[ii-nburn,] = σ_h2;
            eh_store[:,ii-nburn] = eh;
        end
    end

    return store_beta, store_h, store_Σ, store_s2_h, ρ_store, σ_h2_store, eh_store
    
end # end function Chan2020_LBA_csv

#--------------------------------------
#--- FORECASTING BLOCK
#----------------
function fcastChan2020_LBA_csv(YY,VARSetup, store_beta, store_h,store_Σ, store_ρ, store_σ_h2)
    @unpack n_fcst,n,p,nsave = VARSetup

    Yfor3D    = fill(NaN,(p+n_fcst,n,nsave))
    hfor3D    = fill(NaN,(p+n_fcst,nsave)); 
    
    
    Yfor3D[1:p,:,:] .= @views YY[end-p+1:end,:];
    hfor3D[1:p,:] = @views store_h[end-p+1:end,:];
    
    for i_draw = 1:nsave
    
        hfor = @views hfor3D[:,i_draw];
        Yfor = @views Yfor3D[:,:,i_draw];
        A_draw = @views reshape(store_beta[:,i_draw],n*p+1,n);
        ρ_draw = @view store_ρ[i_draw];
        σ_h2_draw = @views store_σ_h2[i_draw];
        Σ_draw = @views store_Σ[:,:,i_draw];
                
        for i_for = 1:n_fcst
            hfor[p+i_for,] = ρ_draw.*hfor[p+i_for-1,] + sqrt(σ_h2_draw).*randn()
            tclass = @views vec(reverse(Yfor[1+i_for-1:p+i_for-1,:],dims=1)')
            tclass = [1;tclass];
            Yfor[p+i_for,:]=tclass'*A_draw  .+ (exp.(hfor[p+i_for,]./2.0)*cholesky(Σ_draw).U*randn(n,1))';    
        end
    end
    return Yfor3D, hfor3D


end # end function fcastChan2020_LBA_csv()


#----------------------------------------
#
# CPZ 2024 functions
#
#----------------------------------------




@doc raw"""
# blkDiagMat_sp,blockMatInd_vec = makeBlkDiag(Tfn::Int,n::Int,p::Int,blockMat)

Initializes a block-diagonal matrix with p optional blocks below the main diagonal matrix and populates it.

After you have done so, you can later update the entries of the matrix by simply using 

```lang=julia
    blkDiagMat_sp.nzval[:] = @views blockMat[blkDiagMat_sp]
```

It can be used to generate the shape of H_B or Σ from Chan, Poon, and Zhu (2024) which hast the linear indices of the coefficient matrix.

"""
function makeBlkDiag(Tfn::Int,n::Int,p::Int,blockMat)

    type = eltype(blockMat)

    # initialize first the large matrix
    blkDiagMat = zeros(type,Tfn,Tfn)

    # The loop first populates the main block-diagonal (ij = 0)
    # then it populates the first off-diagonal (ij = 1), 
    # but stops before the last bottom right block, which is on the main diagonal and under which there is no entry
    for ij = 0:p
        for ii = 1:div(Tfn,n)-ij
            blkDiagMat[ ij*n + (ii-1)*n + 1 : n + (ii-1)*n +  ij*n, (ii-1)*n + 1 : n + (ii-1)*n] = ones(type,n,n)
        end
    end
    blkDiagMat_sp = sparse(blkDiagMat)

    # Now we will use linear indices to map the values in blockMat to their positions in blkDiagMat
    # - first initialize a matrix with Int as indices cannot be other than Int numbers
    blkDiagMatInt_sp = convert.(Int,blkDiagMat_sp)
    # - take the linear indices of blockMat
    blockMatInd = LinearIndices(blockMat)

    # - populate the diagonal of the Int block-diag matrix with the linear indices of blockMat
    for ij = 0:p
        for ii = 1:div(Tfn,n)-ij
            blkDiagMatInt_sp[ ij*n + (ii-1)*n + 1 : n + (ii-1)*n +  ij*n, (ii-1)*n + 1 : n + (ii-1)*n] = blockMatInd[ :, 1 + (ij-0)*n : n + (ij-0)*n]
        end
    end

    # - copy those linear indices to be used
    blockMatInd_vec = deepcopy(blkDiagMatInt_sp.nzval)

       # we can now update
    blkDiagMat_sp.nzval[:] = blockMat[blockMatInd_vec]

    return blkDiagMat_sp,blockMatInd_vec


    # if you would like to use the BlockBandedMatrices package you can do so. The code for the loops above would be
    # blk = blockMat  
    # lin = LinearIndices(blk)
    # bbm = BlockBandedMatrix(ones(Float64,Tfn,Tfn),n*ones(Int,Tf,),n*ones(Int,Tf,),(p,0))
    # bbm_lin = BlockBandedMatrix(ones(Int,Tfn,Tfn),n*ones(Int,Tf,),n*ones(Int,Tf,),(p,0))
    # for ij = 0:p
    #     for ii in 1:div(Tfn,n)-ij
    #         bbm_lin[Block(ii+ij,ii)]=lin[ 1 + (ij-0)*n : n + (ij-0)*n ,:]
    #     end
    # end
    # blklin_ind = deepcopy(bbm_lin.data)
    # blkDiagMat_sp=sparse(bbm)
    # blkDiagMat_sp.nzval[:] = @view blk[blklin_ind];

end




@doc raw"""

"""
function makeMinter(z_tab,YYt,Sm_bit,datesHF,varNamesLF,fvarNames,freq_mix_tp,nm,Tf)
    
    z_var_pos  = indexin(varNamesLF,fvarNames); # positions of the variables in z
    T_z, n_z = size(z_tab);    # number of z vars
    
    M_z = zeros(T_z*n_z,nm)
    z_vec = zeros(T_z*n_z,)
    for ii_z = 1:n_z # iterator going through each variable in z_tab (along the columns)
        datesLF_ii = timestamp(z_tab[varNamesLF[ii_z]])
        iter = CartesianIndices(YYt)
        ym_ci = iter[Sm_bit] # a vector of y_m with cartesian indices of the missing values in YYt
        z_ci = CartesianIndices((z_var_pos[ii_z]:z_var_pos[ii_z],1:Tf))
        z_Mind_vec_ii = vec(sum(ym_ci.==z_ci,dims=2)) # alternative z_Mind_vec=vec(indexin(ym_ci,z_ci)).!==nothing
    
        M_inter_ii = zeros(T_z,nm)
        M_z_ii = @views M_inter_ii[:,z_Mind_vec_ii.==1]
    
        if size(datesHF,1)!==size(M_z_ii,2)
            error("The size of M does not match the number of dates available in z_tab. Maybe the low-frequency data is longer? The problem is with variable number ", z_var_pos[ii_z])
        end
    
        # we need to watch out with the dates due to how the intertemporal constraint works Take for example growth rates Q and M
        # y_t = 1/3 y_t - 2/3 y_{t-1} \dots - - 2/3 y_{t-3} - 1/3 y_{t-5}
        # Intuitively, Q1 quarterly GDP (e.g. 01.01.2000) is the weighted sum of the monthly March, February, January, December, November, and October
        # if y_t^Q is 01.01.2000, we need +2 and -2 months for the weights
        if freq_mix_tp==(1,3,0)
            hfWeights = [1/3; 2/3; 3/3; 2/3; 1/3]; n_hfw = size(hfWeights,1); #number of weights, depends on the variable transformation and frequency
        elseif freq_mix_tp==(3,12,0)
            # quarterly and yearly data with growth rates
            hfWeights = [1/4; 2/4; 3/4; 1; 3/4; 2/4; 1/4]; n_hfw = size(hfWeights,1); #number of weights, depends on the variable transformation and frequency
        else
            error("This combination of frequencies and transformation has not been implemented")
        end
    
        for ii_zi in eachindex(datesLF_ii) # iterator going through each time point in datesHF
            ii_M = findall(datesHF.==datesLF_ii[ii_zi])[1]       # find the low-frequency index that corresponds to the high-frequency missing value
            # M_z_ii[ii_zi, findall(datesHF.==datesLF_ii[ii_zi])[1]-n_hfw+1:findall(datesHF.==datesLF_ii[ii_zi])[1]] = hfWeights # if shifted above
            M_z_ii[ii_zi,ii_M-div((n_hfw-1),2): ii_M+div((n_hfw-1),2)]=hfWeights; # +2 and - 2 months for the weights or +3 and -3
        end
        M_z[(ii_z-1)*T_z + 1:T_z + (ii_z-1)*T_z,:] = M_inter_ii;
        z_vec[(ii_z-1)*T_z + 1:T_z + (ii_z-1)*T_z,]  = values(z_tab[varNamesLF[ii_z]]);
    end
    M_zsp = sparse(M_z);
    return M_zsp, z_vec, T_z
end


include("init_functions.jl")
include("Banbura2010.jl")
include("irfs.jl")


#-------------------------------------
end # END OF MODULE
#-------------------------------------

