include("devPkgs.jl")
using BEAVARs
using TimeSeries
using Parameters
# using DelimitedFiles
using Plots
# using Statistics
# using LinearAlgebra


dataNRW_full = readtimearray("data/nrw_2025_05_03.csv"; format="dd/mm/yyyy", delim=',');
varNames = colnames(dataNRW_full)
varList_HF = varNames[1:9];
dataM_raw = dataNRW_full[varList_HF]

dataM_tab = map((timestamp, values) -> (timestamp, log.(values)), dataM_raw[varNames[[1:4;6:9]]]);
dataM_tab = [dataM_tab dataM_raw[:ifo]./100];

dataQ_tab = log.(collapse(dataNRW_full[:gdpNRW], year, first))


dataLF_tab = dataQ_tab;
dataHF_tab = dataM_tab;
varList   = [colnames(dataLF_tab); colnames(dataHF_tab)];


fdataHF_tab, z_tab, freq_mix_tp, datesHF, varNamesLF, fvarNames = BEAVARs.CPZ_prep_TimeArrays(dataLF_tab,dataHF_tab,varList,trans)

YYwNA = values(fdataHF_tab);
YY = deepcopy(YYwNA);
Tf,n = size(YY);
p = 4; const_loc = 1;
B_draw, structB_draw, Σt_inv, b0 = BEAVARs.initParamMatrices(n,p,const_loc) 

YYt, Y0, longyo, nm, H_B, H_B_CI, strctBdraw_LI, Σ_invsp, Σt_LI, Σp_invsp, Σpt_ind, Xb, cB, cB_b0_LI, Smsp, Sosp, Sm_bit, Gm, Go, GΣ, Kym = BEAVARs.CPZ_initMatrices(YY,structB_draw,b0,Σt_inv,p);

M_zsp, z_vec, T_z, MOiM, MOiz = BEAVARs.CPZ_makeM_inter(z_tab,YYt,Sm_bit,datesHF,varNamesLF,fvarNames,freq_mix_tp,nm,Tf);


z_var_pos  = indexin(varNamesLF,fvarNames); # positions of the variables in z
T_z, n_z = size(z_tab);    # number of z vars

global M_z = zeros(T_z*n_z,nm)
global z_vec = zeros(T_z*n_z,)
for ii_z = 1:n_z # iterator going through each variable in z_tab (along the columns)
    datesLF_ii = timestamp(z_tab[varNamesLF[ii_z]])
    iter = CartesianIndices(YYt)
    ym_ci = iter[Sm_bit] # a vector of y_m with cartesian indices of the missing values in YYt
    z_ci = CartesianIndices((z_var_pos[ii_z]:z_var_pos[ii_z],1:Tf))
    z_Mind_vec_ii = vec(sum(ym_ci.==z_ci,dims=2)) # alternative z_Mind_vec=vec(indexin(ym_ci,z_ci)).!==nothing

    M_inter_ii = zeros(T_z,nm)
    M_z_ii = @views M_inter_ii[:,z_Mind_vec_ii.==1]

    if size(datesHF,1)!==size(M_z_ii,2)
        # error("The size of M does not match the number of dates available in z_tab. Maybe the low-frequency data is longer? The problem is with variable number ", z_var_pos[ii_z])
    end

    # we need to watch out with the dates due to how the intertemporal constraint works Take for example growth rates Q and M
    # y_t = 1/3 y_t - 2/3 y_{t-1} \dots - - 2/3 y_{t-3} - 1/3 y_{t-5}
    # Intuitively, Q1 quarterly GDP (e.g. 01.01.2000) is the weighted sum of the monthly March, February, January, December, November, and October
    # if y_t^Q is 01.01.2000, we need +2 and -2 months for the weights
    if freq_mix_tp==(1,3,0)
        hfWeights = [1/3; 2/3; 3/3; 2/3; 1/3]; n_hfw = size(hfWeights,1); #number of weights, depends on the variable transformation and frequency
        hf_num1 = 1; hf_num2 = 1;  # this solves the range below ii_M-div((n_hfw-hf_num1),2): ii_M+div((n_hfw-hf_num2),2). This should give the indices -2, -1, 0, +1, +2
    elseif freq_mix_tp==(3,12,0)
        # quarterly and yearly data with growth rates
        hfWeights = [1/4; 2/4; 3/4; 1; 3/4; 2/4; 1/4]; n_hfw = size(hfWeights,1); #number of weights, depends on the variable transformation and frequency
        hf_num1 = 1; hf_num2 = 1;  # this solves the range below ii_M-div((n_hfw-hf_num1),2): ii_M+div((n_hfw-hf_num2),2). This should give the indices -3, -2, -1, 0, +1, +2, +3
    elseif freq_mix_tp==(1,3,1)
        hfWeights = [1/3; 1/3; 1/3]; n_hfw = size(hfWeights,1); #number of weights, depends on the variable transformation and frequency
        hf_num1 = 3; hf_num2 = -1;  # this solves the range below ii_M-div((n_hfw-hf_num),2): ii_M+div((n_hfw-hf_num),2). This should give the indices -0, +1, +2
    elseif freq_mix_tp==(3,12,1)
        hfWeights = [1/4; 1/4; 1/4; 1/4]; n_hfw = size(hfWeights,1); #number of weights, depends on the variable transformation and frequency
        hf_num1 = 4; hf_num2 = -3;  # this solves the range below ii_M-div((n_hfw-hf_num),2): ii_M+div((n_hfw-hf_num),2). This should give the indices -0, +1, +2
    else
        error("This combination of frequencies and transformation has not been implemented")
    end

    for ii_zi in eachindex(datesLF_ii) # iterator going through each time point in datesHF
        ii_M = findall(datesHF.==datesLF_ii[ii_zi])[1]       # find the low-frequency index that corresponds to the high-frequency missing value
        # M_z_ii[ii_zi, findall(datesHF.==datesLF_ii[ii_zi])[1]-n_hfw+1:findall(datesHF.==datesLF_ii[ii_zi])[1]] = hfWeights # if shifted above
        M_z_ii[ii_zi,ii_M-div((n_hfw-hf_num1),2): ii_M+div((n_hfw-hf_num2),2)]=hfWeights; # +2 and - 2 months for the weights or +3 and -3
    end
    M_z[(ii_z-1)*T_z + 1:T_z + (ii_z-1)*T_z,:] = M_inter_ii;
    z_vec[(ii_z-1)*T_z + 1:T_z + (ii_z-1)*T_z,]  = values(z_tab[varNamesLF[ii_z]]);
end

out_strct, varSetup,hypSetup = beavar("CPZ2024",dataHF_tab,dataLF_tab,varList,1,n_burn=20,n_save=20);

@unpack M_zsp, store_YY, z_vec, Sm_bit = out_strct

yy1 = dropdims(median(store_YY,dims=3),dims=3);
yy_low =  percentile_mat(store_YY, 0.05; dims=3);
yy_high =  percentile_mat(store_YY, 0.95; dims=3);
plot(yy1[:,1]); plot!(yy_low[:,1]); plot!(yy_high[:,1])
plot(M_zsp*yy1[Sm_bit'])
plot!(z_vec)

Yfor3D = BEAVARs.forecast(out_strct,varSetup);
