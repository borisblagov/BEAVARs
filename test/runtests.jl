using BEAVARs
using Test
@testset "BEAVARs.jl" begin
    # Write your tests here.
    TimeArray(DateTime(2020,1,1):Quarter(1):DateTime(2027,4,1),rand(30,3))
end

