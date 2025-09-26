using BEAVARs
using TimeSeries # required to create a TimeArray
using Dates     # these are required only for the example, your data may already have a time-series format

# Chan2020minn
model_type, hyp_strct, set_strct = makeSetup("Chan2020minn",n_burn=20;n_save=50,p=2)
data = TimeArray(DateTime(2020,1,1):Quarter(1):DateTime(2027,4,1),rand(30,3));
data_strct = BEAVARs.makeDataSetup(model_type,data);
out_strct, varSetup = beavar(model_type, set_strct, hyp_strct, data_strct);

# Chan2020iniw
model_type, hyp_strct, set_strct = makeSetup("Chan2020iniw",n_burn=20;n_save=50,p=2)
data = TimeArray(DateTime(2020,1,1):Quarter(1):DateTime(2027,4,1),rand(30,3));
data_strct = BEAVARs.makeDataSetup(model_type,data);
out_strct, varSetup = beavar(model_type, set_strct, hyp_strct, data_strct);

# Chan2020csv
model_type, hyp_strct, set_strct = makeSetup("Chan2020csv",n_burn=20;n_save=50,p=2)
data = TimeArray(DateTime(2020,1,1):Quarter(1):DateTime(2027,4,1),rand(30,3));
data_strct = BEAVARs.makeDataSetup(model_type,data);
out_strct, varSetup = beavar(model_type, set_strct, hyp_strct, data_strct);


