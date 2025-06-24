include("devPkgs.jl")
using BEAVARs
using TimeSeries
using Parameters
using DelimitedFiles
using Plots
using Statistics
using LinearAlgebra
using Distributions
using SparseArrays


VARSetup, hypSetup =  BEAVARs.beavar_debug("Chan2020csv",n_save=10,n_burn=10)
    

YY = rand(35,2);
YYlist = [:Y1; :Y2; :Y3; :Y4]

@unpack ρ, σ_h2, v_h0, S_h0, ρ_0, V_ρ = hypSetup
@unpack p, nsave, nburn = VARSetup

Y, X, T, n, const_loc = mlagL(YY,p)
(deltaP, sigmaP, mu_prior) = trainPriors(YY,p)
np1 = n*p+1; # number of parameters per equation
    
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
store_β = zeros(np1*n,nsave);
store_h = zeros(T,nsave);
store_Σ = zeros(n,n,nsave);
store_s2_h = zeros(T,nsave);
ρ_store = zeros(nsave,);
σ_h2_store = zeros(nsave,); 
eh_store = zeros(T,nsave);

eh = similar(h);

Ωinv = sparse(1:T,1:T,exp.(-h));
dg_ind_Ωinv = diagind(Ωinv);
# for ii = 1:ndraws 
    ii = 1
    Ωinv[dg_ind_Ωinv]=exp.(-h);
    ZtΩinv = X'*Ωinv;
    
    K_A = V_Ainv + ZtΩinv*X;
    A_hat = K_A\(V_Ainv\A_0 + ZtΩinv*Y);
    S_hat = S_0 + A_0'*V_Ainv*A_0 + Y'*Ωinv*Y - A_hat'*K_A*A_hat;
    S_hat = (S_hat+S_hat')/2;

    Σ = rand(InverseWishart(v_0+T,S_hat));
    cholΣ  = cholesky(Σ).U; # if we get the upper we don't need constant transpose
    A = A_hat + (cholesky(Hermitian(K_A)).U\randn(np1,n))*cholΣ;

    # Errors
    U = Y - X*A
    s2_h = sum((U/cholΣ).^2,dims=2)

    BEAVARs.draw_h_csv!(h,s2_h,ρ,σ_h2,n,H_ρ)
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
        store_β[:,ii-nburn] = vec(A);
        store_h[:,ii-nburn] = h;
        store_Σ[:,:,ii-nburn] = Σ;
        store_s2_h[:,ii-nburn] = s2_h;
        ρ_store[ii-nburn,] = ρ;
        σ_h2_store[ii-nburn,] = σ_h2;
        eh_store[:,ii-nburn] = eh;
    end
# end
