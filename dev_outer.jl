using BEAVARs, TimeSeries
using Dates     # these are required only for the example, your data may already have a time-series format

# Chan2020minn
model_type, hyp_strct, set_strct = makeSetup("Chan2020minn",n_burn=20;n_save=50,p=2)
data = TimeArray(DateTime(2020,1,1):Quarter(1):DateTime(2027,4,1),rand(30,3));
data_strct = BEAVARs.makeDataSetup(model_type,data);

loop_strct = LoopSetup(model_type,set_strct,hyp_strct,data_strct)

vint_dict = Dict{String,BEAVARs.BVARmodelLoopSetup}()
vint_dict["v2"] = loop_strct
vint_dict["v1"] = loop_strct

vint_out_dict, fcast_out_dict = beavars(vint_dict)

out_strct = vint_out_dict["v1"]
fcast_mat = fcast_out_dict["v1"]


YYfcastErr_tph = dropdims(trueYY.-mean(YYfcast3D_mat,dims=3),dims=3)

YYfcastErr_mat = repeat(YYfcastErr_tph,inner=(1,1,2))
YYfcastErr_mat[9,:,2]= fill(NaN,3,)

