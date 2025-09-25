include("devPkgs.jl")
using   BEAVARs, TimeSeries, Parameters, Plots
# using DelimitedFiles
# using Plots
# using Statistics
# using LinearAlgebra

data_ta_full = readtimearray("data/data_tpu.csv"; format="dd/mm/yyyy", delim=',')

YYlist = [:tpuCaldara, :pmiDE, :gdpDE, :invDE, :hicpDE, :euribor];
data_de = data_ta_full[YYlist];
var_names = colnames(data_de)
YY = values(data_de);

YY = rand(30,3)

model_type, hyp_strct, set_strct = makeSetup("BGR2010",n_burn=1000;n_save=1000)
data_strct = BEAVARs.makeDataSetup(model_type,data_de,var_list=var_names)
out_strct_minn, varSetup = beavar(model_type, set_strct, hyp_strct, data_strct);



model_type, hyp_strct, set_strct = makeSetup("Chan2020iniw",n_burn=10;n_save=10)


nsave = 100; nburn = 100;
out_strct_minn, varSetup,hypSetup = beavar("Chan2020minn",YY,n_save=nsave,n_burn=nburn);
out_strct_iniw, varSetup,hypSetup = beavar("Chan2020iniw",YY,n_save=nsave,n_burn=nburn);
out_strct_csv2, varSetup,hypSetup = beavar("Chan2020csv2",YY,n_save=nsave,n_burn=nburn);
out_strct_csv, varSetup,hypSetup = beavar("Chan2020csv",YY,n_save=nsave,n_burn=nburn);
out_strct_bgr, varSetup2,hypSetup2 = beavar("BGR2010",YY,n_save=nsave,n_burn=nburn);


###
dataHF_tab, dataLF_tab, varList = BEAVARs.readSpec("bg_L250703","data/Specifications_mfvar.xlsx");
model_type, hyp_strct, set_strct = makeSetup("CPZ2023",n_burn=10;n_save=10)
data_strct = BEAVARs.makeDataSetup(model_type,dataHF_tab, dataLF_tab,0)
out_strct_cpz, varSetup = beavar(model_type, set_strct, hyp_strct, data_strct);
beavar()



###

var_no = 2
Yfit, Yact = BEAVARs.modelFit(out_strct_cpz,varSetup);
plot([Yfit[:,var_no],Yact[:,var_no]]*100)
Yfit, Yact = BEAVARs.modelFit(out_strct_csv,varSetup);
plot([Yfit[:,var_no],Yact[:,var_no]]*100)

Yfor3D = BEAVARs.forecast(out_strct_cpz,varSetup)


Yfor_low1 = percentile_mat(Yfor3D,0.05,dims=3);
Yfor_low = percentile_mat(Yfor3D,0.16,dims=3);
Yfor_med = percentile_mat(Yfor3D,0.5,dims=3)
Yfor_hih = percentile_mat(Yfor3D,0.84,dims=3);
Yfor_hih1 = percentile_mat(Yfor3D,0.97,dims=3);

ik = 1
plot(Yfor_med[:,ik],w=0;ribbon=(Yfor_med[:,ik]-Yfor_low1[:,ik],Yfor_hih1[:,ik]-Yfor_med[:,ik]),fillalpha = 0.1,color=1,legend=false)
plot!(Yfor_med[:,ik],w=2;ribbon = (Yfor_med[:,ik]-Yfor_low[:,ik],Yfor_hih[:,ik]-Yfor_med[:,ik]),fillalpha=0.05,color=1)

dataM_bg_full = readtimearray("data/bg_julL.csv"; format="dd/mm/yyyy", delim=',');
dataQ_bg_full = readtimearray("data/dataQ_BG.csv"; format="dd/mm/yyyy", delim=',')
varNamesM_full = colnames(dataM_bg_full)

dataM_bg_raw = dataM_bg_full
dataQ_bg_raw = dataQ_bg_full


varNamesM_full = colnames(dataM_bg_full)
varNamesHF = [:survIndustryBG];
varNamesLF = [:gdpBG]
varList   = [varNamesLF; varNamesHF];

# select only the needed data and transform it if needed
dataM_bg_tab = dataM_bg_raw[varNamesHF];
dataQ_bg_tab = percentchange(dataQ_bg_raw[varNamesLF])

dataLF_tab = dataQ_bg_tab;
dataHF_tab = dataM_bg_tab;

out_strct, varSetup,hypSetup = beavar("CPZ2024",dataHF_tab,dataLF_tab,varList,0,n_burn=200,n_save=100);

@unpack M_zsp, store_YY, z_vec, Sm_bit = out_strct

yy1 = dropdims(median(store_YY,dims=3),dims=3);
yy_low =  percentile_mat(store_YY, 0.05; dims=3);
yy_high =  percentile_mat(store_YY, 0.95; dims=3);
plot(yy1[:,1]); plot!(yy_low[:,1]); plot!(yy_high[:,1])
plot(M_zsp*yy1[Sm_bit'])
plot!(z_vec)

Yfor3D = BEAVARs.forecast(out_strct,varSetup);
