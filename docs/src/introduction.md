# General introduction

## Installation
Installing the package follows the typical Julia scheme
```julia
julia> ]
pkg> add BEAVARs
pkg>
```
and, after pressing backspace  to get back to the julia terminal, typing
```julia
julia> using BEAVARs
```


## Usage
The main function of the package is
```julia
beavar(model_type, set_strct, hyp_strct, data_strct)
```
which calls the relevant models and performs the estimation. Using the package boils down to the correct specification of these arguments, for which special functions exist to help you create them. Before jumping in the details let's give a brief overview.

- `model_type`: An object that allows Julia to know which model you want to use and call the relevant functions. 
- `set_strct`:  A structure containing the general VAR setup such as number of lags, number of draws, etc.
- `hyp_strct`:  A structure for setting the hyperparameters for the Bayesian estimation.
- `data_strct`: A structure containing your data.

The first three objects are generated using a helper function `make_setup`. The fourth object is generated using the helper function `makeDataSetup()`. Let us showcase each of these by estimating the following model: [Chan2020minn](@ref).