include("devPkgs.jl")
using BEAVARs
using DelimitedFiles
# using Statistics
# using Plots
using Distributions
using Parameters
using SparseArrays
using LinearAlgebra

YY = readdlm("data/FRED_Chan2020_LBA.txt", ',');

hyper = hypChan2020_CSV();
# A_store,h_store,s2_h_store,ρ_store, σ_h2_store, eh_store = Chan2020_LBA_CSV(YY,hypChan2020_CSV = hyper,nsave=500,nburn=100);
# @time( Chan2020_LBA_CSV(YY,hypChan2020_CSV=hyper,nsave=5000,nburn=100));

@unpack c1, c2, c3, ρ, σ_h2, v_h0, S_h0, ρ_0, V_ρ, q = hyper

Y, Z, T, n = mlag_r(YY,4);
(deltaP, sigmaP, mu_prior) = trainPriors(YY,4);
p=4;
np1 = n*p+1; # number of parameters per equation
nburn = 100
nsave = 500

(idx_kappa1,idx_kappa2, V_Minn, beta_Minn) = prior_NonConj(n,p,sigmaP,c1,c3);

A_0 = reshape(beta_Minn,np1,n);
V_Ainv = sparse(1:np1,1:np1,1.0./V_Minn)

S_0 = Diagonal(sigmaP);
Σ = S_0;
v_0 = n+3;
h = zeros(T,)
H_ρ = sparse(Matrix(1.0I, T, T)) - sparse(ρ*diagm(-1=>repeat(1.0:1.0,T-1)));
# dv = ones(T,); ev = -ρ.*ones(T-1,);
# H_ρ = Bidiagonal(dv,ev,:L)

# This part follows page 19 of Chan, J. (2020)
ndraws = nsave+nburn;
A_store = zeros(np1,n,nsave);
h_store = zeros(T,nsave);
s2_h_store = zeros(T,nsave);
ρ_store = zeros(nsave,);
σ_h2_store = zeros(nsave,); 
eh_store = zeros(T,nsave);

eh = similar(h);

Ωinv = sparse(1:T,1:T,exp.(-h));
dg_ind_Ωinv = diagind(Ωinv);
ii = 1;
# for ii = 1:ndraws 
Ωinv[dg_ind_Ωinv]=exp.(-h);
ZtΩinv = Z'*Ωinv;

K_A = V_Ainv + ZtΩinv*Z;
A_hat = K_A\(V_Ainv\A_0 + ZtΩinv*Y);
S_hat = S_0 + A_0'*V_Ainv*A_0 + Y'*Ωinv*Y - A_hat'*K_A*A_hat;
S_hat = (S_hat+S_hat')/2;

Σ = rand(InverseWishart(v_0+T,S_hat));
CSig_t = cholesky(Σ).U; # if we get the upper we don't need constant transpose
A = A_hat + (cholesky(Hermitian(K_A)).U\randn(np1,n))*CSig_t;

# Errors
U = Y - Z*A
s2_h = sum((U/CSig_t).^2,dims=2); # this is a matrix, must be explicitly converted

draw_h_csv!(h,s2_h,ρ,σ_h2,n,H_ρ);
draw_h_csv_opt!(h,s2_h[:,],ρ,σ_h2,n,H_ρ);


@time draw_h_csv!(h,s2_h,ρ,σ_h2,n,H_ρ);
@time draw_h_csv_opt!(h,s2_h[:,],ρ,σ_h2,n,H_ρ);

@btime draw_h_csv!(h,s2_h,ρ,σ_h2,n,H_ρ);
@btime draw_h_csv_opt!(h,s2_h[:,],ρ,σ_h2,n,H_ρ);


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


# if ii>
    A_store[:,:,ii-nburn] = A;
    h_store[:,ii-nburn] = h;
    s2_h_store[:,ii-nburn] = s2_h;
    ρ_store[ii-nburn,] = ρ;
    σ_h2_store[ii-nburn,] = σ_h2;
    eh_store[:,ii-nburn] = eh;
# end
# end