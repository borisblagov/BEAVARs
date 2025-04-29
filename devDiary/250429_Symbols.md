To use data structures effectively I stumbled upon Symbols in Julia. These are my notes for future reference and how I may need them in coding.

My use case is that I want to have a dataset with all possible variables in there and then manipulate that list, for example by referencing only some of them, reordering, or transforming them. This is from the perspective of TimeSeries.jl package but I imagine it is usefuf for others such as DataFrames.

Searching the internet I found some useful references, mainly 
- [This stackoverflow thread](https://stackoverflow.com/questions/23480722/what-is-a-symbol-in-julia/23482257#23482257)
- [Metaprogramming documentation](https://docs.julialang.org/en/v1/manual/metaprogramming/)

My use case
```julia
    using TimeSeries.jl
    dataM_bg_full = readtimearray("data/dataM_BG.csv"; format="dd/mm/yyyy", delim=',')
    varnamesM_full = colnames(dataM_bg_full)
```

`varnamesM_full` has the Symbol type

```julia
julia> varnamesM_full = colnames(dataM_bg_full)
16-element Vector{Symbol}:
 :fdiBalBG
 ...
 :survConsBG
 ```

 The symbol is "bound" to a variable. You can use `dump()` to see what the symbol refers to. An example where we bind the string `"bar"` to the symbol `:foo` and store that in an expression `ex`

```
julia> ex = :(foo = "bar")
:(foo = "bar")

julia> dump(ex)
Expr
  head: Symbol =
  args: Array{Any}((2,))
    1: Symbol foo
    2: String "bar"

julia> eval(ex)
"bar"
```

For my use

```julia
 julia> dump(varnamesM_full)
Array{Symbol}((16,))
  1: Symbol fdiBalBG
  2: Symbol bnbBalBG
  ...
  15: Symbol arrivalsBG
  16: Symbol survConsBG
```
