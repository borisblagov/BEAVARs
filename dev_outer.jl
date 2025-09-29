using BEAVARs, TimeSeries
using Dates     # these are required only for the example, your data may already have a time-series format

# Chan2020minn
model_type, hyp_strct, set_strct = makeSetup("Chan2020minn",n_burn=20;n_save=50,p=2)
data = TimeArray(DateTime(2020,1,1):Quarter(1):DateTime(2027,4,1),rand(30,3));
data_strct = BEAVARs.makeDataSetup(model_type,data);
out_strct = beavar(model_type, set_strct, hyp_strct, data_strct)

propertynames(out_strct)

trueYY = rand(10,3)
trueYY[10,:] = fill(NaN,3,)
YYfcast3D_mat = BEAVARs.forecast(out_strct,set_strct)

YYfcastErr_tph = dropdims(trueYY.-mean(YYfcast3D_mat,dims=3),dims=3)

YYfcastErr_mat = repeat(YYfcastErr_tph,inner=(1,1,2))
YYfcastErr_mat[9,:,2]= fill(NaN,3,)

_nanfunc(f, A, ::Colon) = f(filter(!isnan, A))
_nanfunc(f, A, dims) = mapslices(a->_nanfunc(f,a,:), A, dims=dims)
nanfunc(f, A; dims=:) = _nanfunc(f, A, dims)