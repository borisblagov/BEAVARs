using Revise
using BayesianVAR
using BenchmarkTools
using DelimitedFiles
using CSV
using Parameters
using SparseArrays
using LinearAlgebra
using Plots


Y_20 = readdlm("data/FRED_Chan2019.txt", ',');
Y_4 = Y_20[:,[1,7,8,12]];

hyper=hypChan2020()
store_beta = Chan2020_LBA(Y_4)
Chan2020_LBA(Y_20,p=1,nburn = 100, nsave = 5000)
Chan2020_LBA(Y_4,p=1,nburn = 100, nsave = 5000)