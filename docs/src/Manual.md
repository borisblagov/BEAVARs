# Manual

## Installation
Installing the package follows the typical Julia scheme
```
pkg> add BEAVARs
```

## Usage
The main function of the package is the function 
```
beavar(model_type, data_strct, set_strct, hyp_strct)
```
which calls the relevant models and performs the estimation. Using the package boils down to the correct specification of its arguments. Therefore, before jumping in details we will sketch what they are and how they are defined.

- `model_type`: A custom type of the package, i.e. a special object that allows Julia to know which function to call. 
- `hyp_strct`:  A structure with hyperparameter values for the Bayesian estimation.
- `set_strct`: A structure with general VAR setup such as number of lags, number of draws, etc.
- `data_strct`: A structure containing your data.