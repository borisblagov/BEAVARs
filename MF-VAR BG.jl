include("devPkgs.jl")
using BEAVARs
using TimeSeries
using Parameters
using Plots
using DelimitedFiles
# using Statistics
# using LinearAlgebra
# using Distributions
# using SparseArrays

include("src/plot_functions.jl")
trans = 0;
dataHF_tab, dataLF_tab, varList = BEAVARs.readSpec("bg_L250911","data/Specifications_mfvar.xlsx");
out_strct, varSetup,hypSetup = beavar("CPZ2024",dataHF_tab,dataLF_tab,varList,0,n_burn=10000,n_save=10000);
Yfor_low1, Yfor_low, Yfor_med, Yfor_hih, Yfor_hih1 = fanChart(out_strct.store_YY[:,1,:].*100, timestamp(dataHF_tab));
Plots.ylabel!("percent"); Plots.title!("Monthly GDP growth")

out_strct2, varSetup2,hypSetup2 = beavar("Blagov2025",dataHF_tab,dataLF_tab,varList,0,n_burn=10000,n_save=10000);
Yfor2_low1, Yfor_low2, Yfor_med2, Yfor_hih2, Yfor_hih12 = fanChart(out_strct2.store_YY[:,1,:].*100, timestamp(dataHF_tab));
Plots.ylabel!("percent"); Plots.title!("Monthly GDP growth")


data = (datetime = timestamp(dataHF_tab), col0 = Yfor_low1, col1 = Yfor_low, col2 = Yfor_med,col3 = Yfor_hih, col4 = Yfor_hih1);
ta = TimeArray(data; timestamp = :datetime, meta = "Example");
plot(ta.*100)

# out_strct2, varSetup,hypSetup = beavar("Blagov2025",dataHF_tab,dataLF_tab,varList,0,n_burn=100,n_save=100);


Yfor3D = BEAVARs.forecast(out_strct,varSetup);
# fanChart(Yfor3D[:,1,:])


YY_HF_med = percentile_mat(out_strct.store_YY,0.5,dims=3);
# YY_HF_med2 = percentile_mat(out_strct2.store_YY,0.5,dims=3);
@unpack M_zsp, store_YY, z_vec, Sm_bit = out_strct
plot(M_zsp*YY_HF_med'[Sm_bit])
plot!(z_vec)



@unpack p = varSetup;
Y,X,T,n = BEAVARs.mlagL(YY_HF_med,p);
Amed = reshape(percentile_mat(out_strct.store_β,0.5,dims=2),n*p+1,n)
global Yfit = X*Amed;
global Yact = @views Y
plot([Yfit[:,1],Yact[:,1]]*100)


@unpack p = varSetup;
Y,X,T,n = BEAVARs.mlagL(YY_HF_med2,p);
Amed = reshape(percentile_mat(out_strct2.store_β,0.5,dims=2),n*p+1,n)
global Yfit2 = X*Amed;
global Yact2 = @views Y
plot([Yfit2[:,1],Yact2[:,1]]*100)

# n = setup_str.n
plt = plot(layout=(5,3))
for i_var = 1:n;
    plot!(plt,[Yfit[:,i_var] Yact[:,i_var]],subplot=i_var)
end
display(plt)

percentliles_mat =[Yfor_low1 Yfor_low Yfor_med Yfor_hih Yfor_hih1]