using BayesianAnimalMovementModels
using Distributions, Random
using LinearAlgebra, PDMats
using Impute
using Test

const tests = [
    "functions",
    "models",
]

printstyled("Running tests:\n", color=:blue)

for t in tests
    @testset "Test $t" begin
        Random.seed!(345679)
        include("$t.jl")
    end
end
