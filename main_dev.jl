include("devPkgs.jl")
using BEAVARs
using DelimitedFiles
using Plots
using Parameters
using Statistics
using TimeSeries
using LinearAlgebra



data_ta_full = readtimearray("data/data_tpu.csv"; format="dd/mm/yyyy", delim=',')

data_de = data_ta_full[[:tpuCaldara, :pmiDE, :gdpDE, :invDE, :hicpDE, :euribor]];
var_names = colnames(data_de)
YY = values(data_de);

# YY20 = readdlm("data/FRED_Chan2020_LBA.txt", ',');
# YY = YY20[:,[1,4,5,6]]
setup_str, hyper_str = makeSetup(YY,"Chan2020_LBA_csv",nburn=1200,nsave=5000,n_fcst=200)
store_beta, store_h, store_Σ, s2_h_store, store_ρ, store_σ_h2, eh_store = Chan2020_LBA_csv(YY,setup_str,hyper_str);

# Using Impulse responses
IRF_median, IRF_68_low, IRF_68_high = irf_chol_overDraws_csv(store_beta,store_Σ,store_h,setup_str);

# store_beta, store_Σ =  Chan2020_LBA_Minn(YY;hyp=hypChan2020,VARSetup = setup_str);
# IRF_median, IRF_68_low, IRF_68_high = irf_chol_overDraws(store_beta,store_Σ,setup_str;shSize = "stdev");

i_var = 3;
i_shock = 1;            # which variable you shock
plot(IRF_median[:,i_var,i_shock],ribbon = (IRF_median[:,i_var,i_shock].-IRF_68_low[:,i_var,1],IRF_68_high[:,i_var,i_shock].-IRF_median[:,i_var,1]),fillalpha=0.2) # ,label="",title=var_names[i_var]

@unpack n_fcst,n,p,nsave=setup_str
Amed = reshape(median(store_beta[:,:],dims=2),n*p+1,n);

global Yfit = XX*Amed;
global Yact = @views YY[p+1:end,:]
plot([Yfit[:,1],Yact[:,1]])

# n = setup_str.n
plt = plot(layout=(3,2))
for i_var = 1:n;
    plot!(plt,[Yfit[:,i_var] Yact[:,i_var]],subplot=i_var)
end
display(plt)



# Forecasting
i_for = 1
i_draw = 1;
Yfor3D    = fill(NaN,(p+n_fcst,n,nsave))
hfor3D    = fill(NaN,(p+n_fcst,nsave)); 


Yfor3D[1:p,:,:] .= @views YY[end-p+1:end,:];
hfor3D[1:p,:] = @views store_h[end-p+1:end,:];


hfor = @views hfor3D[:,i_draw];
Yfor = @views Yfor3D[:,:,i_draw];
A_draw = @views reshape(store_beta[:,i_draw],n*p+1,n);
ρ_draw = @view store_ρ[i_draw];
σ_h2_draw = @views store_σ_h2[i_draw];
Σ_draw = @views store_Σ[:,:,i_draw];
        
@time for i_for = 1:n_fcst
    hfor[p+i_for,] = ρ_draw.*hfor[p+i_for-1,] + sqrt(σ_h2_draw).*randn()
    tclass = @views vec(reverse(Yfor[1+i_for-1:p+i_for-1,:],dims=1)')
    tclass = [1;tclass];
    Yfor[p+i_for,:]=tclass'*A_draw  .+ (exp.(hfor[p+i_for,]./2.0)*cholesky(Σ_draw).U*randn(n,1))';    
end


Yfor3D, hfor3D = fcastChan2020_LBA_csv(YY,setup_str, store_beta, store_h,store_Σ, store_ρ, store_σ_h2);


Yfor_med = median(Yfor3D,dims=3)

plot(Yfor_med[:,4])