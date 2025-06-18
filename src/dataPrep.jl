
function TAtrans(dataHF_TA,varList_HF,trans_dictA)
transHF_vec = getindex.(Ref(trans_dictA),varList_HF);       # vector of transformations for the vars

    # transformations
# if 1, leave it be, initialize TA
hf_ta = map((timestamp, values) -> (timestamp, (values)), dataHF_TA[varList_HF[transHF_vec.==1]])
#if 2, divide by 100
hf_ta_temp = map((timestamp, values) -> (timestamp, (values)./100.0), dataHF_TA[varList_HF[transHF_vec.==2]])
hf_ta = merge(hf_ta,hf_ta_temp)
# if 3, take logs
hf_ta_temp = map((timestamp, values) -> (timestamp, log.(values)), dataHF_TA[varList_HF[transHF_vec.==3]])
hf_ta = merge(hf_ta,hf_ta_temp)
# if 7 take percentage changes
hf_ta_temp = percentchange(dataHF_TA[varList_HF[transHF_vec.==7]]);
hf_ta = merge(hf_ta,hf_ta_temp)
# if 8, multiply by 100
hf_ta_temp = map((timestamp, values) -> (timestamp, (values).*100.0), dataHF_TA[varList_HF[transHF_vec.==6]])
hf_ta = merge(hf_ta,hf_ta_temp)
# if 9, take exp
hf_ta_temp = map((timestamp, values) -> (timestamp, exp.(values)), dataHF_TA[varList_HF[transHF_vec.==9]])
hf_ta = merge(hf_ta,hf_ta_temp)

hf_ta = hf_ta[varList_HF];  # reorder back the variables

return hf_ta
end


function xlsx2ta(data_mat)
    data_mat = data_mat[.!ismissing.(data_mat[:,1]),:];             # sometimes XLSX reads empty rows below the last one (especially conditional formatting), here we remove them

    nan_mat =  fill(NaN, size(data_mat,1),size(data_mat,2));
    data_mat[ismissing.(data_mat)]=nan_mat[ismissing.(data_mat)];   # replace all missing with NaN
    values_mat = convert(Array{Float64},data_mat[2:end,2:end])      # convert values to numbers
    date_vec = DateTime.(data_mat[2:end,1])                         # convert first column to DateTime
    # date_vecstr = Dates.DateFormat.(date_vec,"yyyy-mm-dd");

    dataf_TA = TimeArray(date_vec,values_mat,Symbol.(data_mat[1,2:end]))
    return dataf_TA
end


function readSpec(modelstr,data_path)
    xf = XLSX.readxlsx(data_path);
    sh_names = XLSX.sheetnames(xf)

    sh_ref = xf["setup"];
    sh_mat = sh_ref[:]

    # Reads the setup sheet
    model_ind = findall(sh_mat[1,:].==modelstr)[1]
    vb_ind = findall(sh_mat[:,1].=="lastRow")[1] + 1;       # where the variables start
    varListA_str = sh_mat[vb_ind:end,1];                    # strings of variables
    varListA_sym = Symbol.(varListA_str[:]);                # symbols of variables

    vm_bit = .!iszero.(sh_mat[vb_ind:end,model_ind]);       # boolean list of ALL variables
    vm_trans = sh_mat[vb_ind:end,model_ind];                # list of variables with transformations for ALL vars

    varListF = varListA_sym[vm_bit]
    trans_dictA = Dict(varListA_sym .=> vm_trans);          # transformation dictionary vector for ALL vars. 



    # Reads the high-frequency
    datasheetHF_str = "datasheet" * string(sh_mat[sh_mat[:,1].=="datasheet",model_ind][1]) * "_HF";    # the high-frequency datasheet
    data_mat = xf[datasheetHF_str][:];
    dataf_HF_TA = BEAVARs.xlsx2ta(data_mat)
    varList_HF = intersect(varListF,colnames(dataf_HF_TA));      # looks for which variables are required and which are found
    dataHF_TA = dataf_HF_TA[varList_HF];                         # selects the variables found in this TA

    hf_ta = BEAVARs.TAtrans(dataHF_TA,varList_HF,trans_dictA)


    # Reads the low-frequency
    datasheetLF_str = "datasheet" * string(sh_mat[sh_mat[:,1].=="datasheet",model_ind][1]) * "_LF";    # the low-frequency datasheet
    data_mat = xf[datasheetLF_str][:];
    dataf_LF_TA = BEAVARs.xlsx2ta(data_mat)
    varList_LF = intersect(varListF,colnames(dataf_LF_TA));      # looks for which variables are required and which are found
    dataLF_TA = dataf_LF_TA[varList_LF];                         # selects the variables found in this TA

    # transHF_vec = Vector{Int}();                                # vector of transformations for the found high-freq. vars
    # push!(transHF_vec,(trans_dictA[i] for i in varList_HF)...);  # fill it
    # transHF_vec = getindex.(Ref(trans_dictA),varList_HF);       # vector of transformations for the vars

    lf_ta = BEAVARs.TAtrans(dataLF_TA,varList_LF,trans_dictA)
    
    return hf_ta, lf_ta, varListF

end