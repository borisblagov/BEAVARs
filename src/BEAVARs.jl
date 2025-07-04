module BEAVARs
using   Revise, 
        LinearAlgebra, 
        Distributions, 
        SparseArrays,
        TimeSeries, 
        Parameters,
        ProgressMeter,
        XLSX

# from init_functions.jl
export mlag, mlagL, mlagL!, ols, percentile_mat

# from Banbura, Giannone, Reichling 2010
export makeDummiesMinn!, makeDummiesSumOfCoeff!, getBeta!, getSigma!, gibbs_beta_sigma,trainPriors, BGR2010, hypBGR2010

# from irfs.jl
export irf_chol, irf_chol_overDraws, irf_chol_overDraws_csv

# from Chan_2020
export prior_Minn, Chan2020minn, Chan2020csv, prior_NonConj
export hypChan2020

export makeSetup, beavar, dispatchModel, makeOutput

# Structures, to be uncommented later
# export modelSetup, modelOutput, Chan2020csv_type, Chan2020minn_type, modelHypSetup, hypDefault_strct, outChan2020csv, VARModelType, VARSetup




# Types for multiple dispatch across models
abstract type VARModelType end      # types for models
abstract type modelSetup end        # type for VAR setup parameters
abstract type modelOutput end       # type for output storage
abstract type modelHypSetup end     # types for hyperparameters


# Structures for multiple dispatch across models
struct Chan2020minn_type <: VARModelType end
struct Chan2020iniw_type <: VARModelType end
struct Chan2020iniw_type2 <: VARModelType end
struct Chan2020csv_type <: VARModelType end
struct Chan2020csv_type2 <: VARModelType end
struct BGR2010_type <: VARModelType end
struct CPZ2024_type <: VARModelType end
struct Blagov2025_type <: VARModelType end
struct hypDefault_strct <: modelHypSetup end    # empty structure for initialising the hyperparameters


function selectModel(model_str::String)
    if model_str == "Chan2020csv"
        model_type = Chan2020csv_type()
    elseif model_str == "Chan2020csv2"
        model_type = Chan2020csv_type2()
    elseif model_str == "Chan2020minn"
        model_type = Chan2020minn_type()
    elseif model_str == "Chan2020iniw"
        model_type = Chan2020iniw_type()
    elseif model_str == "Chan2020iniw2"
        model_type = Chan2020iniw_type2()
    elseif model_str == "BGR2010"
        model_type = BGR2010_type()
    elseif model_str == "CPZ2024"
        model_type = CPZ2024_type()
    elseif model_str == "Blagov2025"
        model_type = Blagov2025_type()
    else
        error("Model not found, make sure the spelling is completely correct, upper and lowercase matters!\n Possible models are: \n    BGR2010 \n    Chan2020minn\n    Chan2020csv\n    Chan2020iniw2\n")
    end
    return model_type
end


# structure initializing the VAR
# Constructor
@doc raw"""
# set_strct = VARSetup(p,n_save,n_burn,n_irf,n_fcst,const_loc);
Populates the constructor VARSetup with default and/or custom values.

    p:      number of lags, default is 4
    nburn:  number of burn-in draws that will be discarded, default is 2000
    nsave:  number of retained draws (total is then nburn + nsave), default is 1000
    n_irf:  horizon of impulse responses, default is 16
    n_fcst: horizon of forecasting periods, default is 8

Outputs

    VARSetup: the setup structure for the BEAVARs
```
"""
@with_kw struct VARSetup <: modelSetup
    p::Int          # number of lags
    nsave::Int      # gibbs to save
    nburn::Int      # gibbs to burn
    n_irf::Int      # number of impulse responses
    n_fcst::Int     # number of forecast periods
    const_loc::Int  # location of the constant
end


# @with_kw struct setupCPZ <: modelSetup
#     p::Int          # number of lags
#     nsave::Int      # gibbs to save
#     nburn::Int      # gibbs to burn
#     n_irf::Int      # number of impulse responses
#     n_fcst::Int     # number of forecast periods
#     const_loc::Int  # location of the constant
#     dataHF_tab:TimeArray     #
#     dataLF_tab:TimeArray     #
#     aggWgh::Int:     #
# end

# ------------------------
# MAIN FUNCTION 
function beavar(model_str=model_name::String,YY_tup... ;p::Int=4,n_burn::Int=1000,n_save::Int=1000,n_irf::Int=16,n_fcst::Int = 8,hyp::modelHypSetup=hypDefault_strct())
    model_type = BEAVARs.selectModel(model_str)
    
    # checking if user supplied the hyperparameter structure
    if isa(hyp,hypDefault_strct)                    # if not supplied, make a default one
        hyp_strct = BEAVARs.makeHypSetup(model_type); # println("using the default hyperparameters")
    else                                            # else use supplied    
        hyp_strct = hyp; # println("using the supplied parameters")
    end
        
    out_strct, set_strct = dispatchModel(model_type,YY_tup, hyp_strct,p,n_burn,n_save,n_irf,n_fcst);
    return out_strct, set_strct, hyp_strct
end


function beavar_debug(model_str=model_name::String,YY_tup... ;p::Int=4,n_burn::Int=1000,n_save::Int=1000,n_irf::Int=16,n_fcst::Int = 8,hyp::modelHypSetup=hypDefault_strct())
    model_type = BEAVARs.selectModel(model_str)
    
    # checking if user supplied the hyperparameter structure
    if isa(hyp,hypDefault_strct)                    # if not supplied, make a default one
        hyp_strct = BEAVARs.makeHypSetup(model_type); # println("using the default hyperparameters")
    else                                            # else use supplied    
        hyp_strct = hyp; # println("using the supplied parameters")
    end
        
    set_strct = VARSetup(p,n_save,n_burn,n_irf,n_fcst,1)
    return set_strct, hyp_strct
end

function makeVARstrct(model_str=model_name::String,dataHF_tab::TimeArray,dataLF_tab::TimeArray,aggWgh::Int;p::Int=4,n_burn::Int=1000,n_save::Int=1000,n_irf::Int=16,n_fcst::Int = 8,hyp::modelHypSetup=hypDefault_strct())
    model_type = BEAVARs.selectModel(model_str)
    
    # checking if user supplied the hyperparameter structure
    if isa(hyp,hypDefault_strct)                    # if not supplied, make a default one
        hyp_strct = BEAVARs.makeHypSetup(model_type); # println("using the default hyperparameters")
    else                                            # else use supplied    
        hyp_strct = hyp; # println("using the supplied parameters")
    end
        
    set_strct = VARSetup(p,n_save,n_burn,n_irf,n_fcst,1)
end

include("dataPrep.jl")
include("init_functions.jl")
include("BGR2010.jl")
include("irfs.jl")
include("Chan2020.jl")
include("Chan2020minn.jl")
include("Chan2020iniw.jl")
include("Chan2020csv.jl")
include("CPZ2024.jl")
include("Blagov2025.jl")



#-------------------------------------
end # END OF MODULE
#-------------------------------------

