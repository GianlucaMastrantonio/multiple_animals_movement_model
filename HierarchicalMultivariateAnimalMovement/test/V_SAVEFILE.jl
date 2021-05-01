using Revise
using BayesianAnimalMovementModels
using Distributions, Random
using LinearAlgebra, PDMats
using Impute
using Test
# Pkg.test("BayesianAnimalMovementModels")
# @code_llvm
# @code_lowered
# @code_native
# @code_typed
# @code_warntype



IndexA   = Int32[2;4]
#IndexA   = Int32[1;2;3;4;5;6]
CovInv = rand(InverseWishart(20, Matrix{Float64}(I,6,6)  ) )
Mean    = [1.0;2.0;3.0;-0.4;0.1;-10.9]
Obs     = rand(6)

CondMean, CondVar, CondVarInv = BayesianAnimalMovementModels.compute_CondMeanAndVariance_MvNormal(IndexA,Mean, CovInv,Obs)

@testset "zeta" begin
    IndexB = [1;3;5;6]
    Cov = inv(CovInv)
    CondMean2 = Mean[IndexA]+Cov[IndexA,IndexB]*inv(Cov[IndexB,IndexB])*(Obs[IndexB]-Mean[IndexB]  )
    CondVar2 = Cov[IndexA,IndexA]-Cov[IndexA,IndexB]*inv(Cov[IndexB,IndexB])*Cov[IndexB,IndexA]
    @test CondMean ≈ CondMean2
    @test CondVar  ≈ CondVar2
end

#### #### #### #### #### #### #### #### #### #### #### ####
#### #### #### #### #### #### #### #### #### #### #### ####
#### #### #### FUNCTIONs
#### #### #### #### #### #### #### #### #### #### #### ####
#### #### #### #### #### #### #### #### #### #### #### ####

#### MISSING

GenMod = BayesianAnimalMovementModels
# create dataset
Dataset = Matrix{Union{Float64,Missing}}(undef,5,4)
Dataset[1,:] = [1 2 3 4]
Dataset[2,:] = [1 2 3 4].+2
Dataset[3,:] = [1 2 3 4].+5
Dataset[4,:] = [1 2 3 4].+2
Dataset[5,:] = [1 2 3 4].+1
Dataset[2,1:2] = [missing missing]
Dataset[4,4]   = missing

MCMCLikelihood = OptionsLikelihood(
    # numer of regimes
    kmax                = Int16(3),
    # data
    data                = Dataset,
    # model type
    likelihood_type     = "OU",

    # Missing
    update_missing      = true,
    # mu0
    update_mu0          = true,
    inits_mu0           = [zeros(size(Dataset,2),3)[:,i] for i in 1:size(zeros(size(Dataset,2),3),2)],
    prior_mu0           = Dict("name"=>"MvNormal","μ"=>zeros(4), "Σ"=>Matrix{Float64}(I,4,4)),
    # psi
    update_psi          = true,
    inits_psi           = [zeros(Integer(size(Dataset,2)/2),3)[:,i] for i in 1:size(zeros(Integer(size(Dataset,2)/2),3),2)],
    prior_psi           = Dict("name"=>"MvUniform","a"=>zeros(2), "b"=>ones(2)),
    # Sigma
    update_sigma        = true,
    inits_sigma         = [Matrix{Float64}(I,4,4), Matrix{Float64}(I,4,4), Matrix{Float64}(I,4,4)],
    prior_sigma         = Dict("name"=>"InverseWishart","df"=>Float64(10), "Ψ"=>Matrix{Float64}(I,4,4)),
    # zeta
    update_zeta         = true,
    inits_zeta          = Int16.(sample(1:3, 5, replace = true))
)

MCMCClusterization = OptionsClusterization(
    Likelihood =           MCMCLikelihood,
    clustering_type     = "HMM"::String, # one of "HMM"
    # pi
    update_pi           = true::Bool,
    inits_pi            = [ones(3)/3.0,ones(3)/3.0,ones(3)/3.0]::Vector{Vector{Float64}},
    prior_pi            = Dict("name"=>"Dirichlet","alpha"=>ones(3)/3.0)
    )

MCMCout = OptionsMCMC(MCMCLikelihood,MCMCClusterization;
        iterations = Int64(100),
        burnin = Int64(50),
        thin = Int64(10),

        )

# @testset "zeta" begin
#     MCMCClusterization.clusterization.zeta[1] = 1
#     @test MCMCClusterization.clusterization.zeta[1]    ==  MCMCLikelihood.clusterization.zeta[1]
#     @test MCMCClusterization.clusterization.zeta       === MCMCLikelihood.clusterization.zeta
# end

#### #### #### #### #### #### #### ####
#### TEST
#### #### #### #### #### #### #### ####

using BayesianAnimalMovementModels
ModelOUT = MCMCalgorithm(
    MCMCout,
    MCMCLikelihood,
    MCMCClusterization
)
#
# appBurninOrThin = MCMCout.burnin
# for saveMCMCindex in 1:MCMCout.nsamplesave
#
#     for burnithinMCMCindex in 1:appBurninOrThin
#
#         ### ### ### ### ### ### ###
#         ### ### LIKELIHOOD
#         ### ### ### ### ### ### ###
#
#         ### missing
#         TODELETE.sample_missing!(MCMCLikelihood)
#
#         ### mu0
#         TODELETE.sample_mu0!(MCMCLikelihood,MCMCLikelihood.mu0)
#
#         ### sigma
#         TODELETE.sample_sigma!(MCMCLikelihood, MCMCLikelihood.sigma)
#
#         ### psi
#         TODELETE.sample_psi!(MCMCLikelihood, MCMCLikelihood.psi)
#
#         ### ### ### ### ### ### ###
#         ### ### ZETA
#         ### ### ### ### ### ### ###
#
#         # zeta
#         TODELETE.sample_zeta!(MCMCLikelihood, MCMCClusterization)
#
#         ### ### ### ### ### ### ###
#         ### ### CLUSTERIZATION
#         ### ### ### ### ### ### ###
#
#         ### pi
#         TODELETE.sample_pi!(MCMCClusterization, MCMCClusterization.pi)
#
#     end # second for MCMC
#     appBurninOrThin = MCMCindices.thin
#
#     TODELETE.save_posteriorsamples!(MCMCout,MCMCLikelihood, MCMCClusterization)
#
# end # first for MCMC

#end
#using BayesianAnimalMovementModels
