# Additional functions

```@docs
makeSetup(YY::Array{Float64},model=model_str; p::Int,nburn::Int,nsave::Int,n_irf::Int,n_fcst::Int)
```

```@docs
ols(Y,X)
```

```@docs
mlag(Yfull::Matrix{Float64},p::Integer)
```

```@docs
mlag_r(Yfull::Matrix{Float64},p::Integer)
```


```@docs
trainPriors(Z0::Matrix{Float64},p::Int64)
```

```@percentile_mat
percentile_mat(A, p; dims)
```