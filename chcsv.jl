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

include("src/plot_functions.jl")

VARSetup, hypSetup =  BEAVARs.beavar_debug("Chan2020csv",n_save=10,n_burn=10)
    

YY = rand(35,2); 
varList = [:Y1; :Y2]

@unpack ρ, σ_h2, v_h0, S_h0, ρ_0, V_ρ = hypSetup
@unpack p, nsave, nburn = VARSetup

Y, X, T, n, k, sigmaP, S_0, Σ, A_0, V_Ainv, v_0, H_ρ,h,eh,Ωinv, dg_ind_Ωinv,ZtΩinv,K_A,A_hat, S_hat, VAinvDA0, AVAinvA = BEAVARs.initcsv(YY,p,hypSetup);

# This part follows page 19 of Chan, J. (2020)
ndraws = nsave+nburn;
store_β = zeros(k*n,nsave);
store_h = zeros(T,nsave);
store_Σ = zeros(n,n,nsave);
store_s2_h = zeros(T,nsave);
ρ_store = zeros(nsave,);
σ_h2_store = zeros(nsave,); 
eh_store = zeros(T,nsave);

# for ii = 1:ndraws 
    ii = 1
    A, cholΣU, Σ = BEAVARs.Chan2020_drawA(Y,X,n,k,T,v_0,h,Ωinv,dg_ind_Ωinv,V_Ainv,S_0,VAinvDA0,AVAinvA);
    h, s2_h = BEAVARs.Chan2020_draw_h!(h,A,Y,X,cholΣU,ρ,σ_h2,n,H_ρ,T);
    ρ, σ_h2, eh = BEAVARs.Chan2020_draw_ρ!(ρ,h,eh,v_h0,S_h0,ρ_0,V_ρ,T);
   
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


YY20 = readdlm("data/FRED_Chan2020_LBA.txt", ',');
YY = YY20[:,[1,4,5,6]]
varList = [:Y1; :Y2;:Y3;:Y4]
YY = YY20;
varList = [:Y1; :Y2;:Y3;:Y4;:Y1; :Y2;:Y3;:Y4;:Y1; :Y2;:Y3;:Y4;:Y1; :Y2;:Y3;:Y4;:Y1; :Y2;:Y3;:Y4]
out_strct2, varSetup2,hypSetup2 = beavar("Chan2020csv2",YY,n_save=10,n_burn=10);
out_strct2, varSetup2,hypSetup2 = beavar("Chan2020csv2",YY,n_save=10000,n_burn=10000);

forc2 = BEAVARs.forecast(out_strct2,varSetup2);

fanChart(forc2)
out_strct, varSetup,hypSetup = beavar("Chan2020csv",YY,n_save=10,n_burn=10);
out_strct, varSetup,hypSetup = beavar("Chan2020csv",YY,n_save=10000,n_burn=10000);

forc = BEAVARs.forecast(out_strct,varSetup);
fanChart(forc)

fanChart(out_strct2.store_h)
fanChart(out_strct.store_h)
plot!(median(out_strct2.store_h,dims=2))

@time out_strct3, varSetup3,hypSetup = beavar("Chan2020iniw",YY,n_save=10,n_burn=10);
out_strct3, varSetup3,hypSetup = beavar("Chan2020iniw",YY,n_save=100,n_burn=100);

Yfor3D = BEAVARs.forecast(out_strct3,varSetup3);
include("fcast_plot.jl")