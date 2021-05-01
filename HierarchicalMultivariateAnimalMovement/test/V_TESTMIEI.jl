##### Packages
using Revise
using BayesianAnimalMovementModels
using Distributions, Random
using LinearAlgebra, PDMats
using Test

##### JULIA PARAMETERS
rng     = MersenneTwister(123);
##### Simulation Parameters

nt = 1500    # number observations
K  = 3       # number total groups
nanim = 2    # numner of animales

# coordinates of mu0 (K,Animale, (x,y))
mu0 = Vector{Vector}(
    [ [0.0; 0.0; 10.0;  -10.0], [ 10.0 ; 10.0; 0.0; 0.0], [-10.0; 10.0;  10.0;  10.0] ]
        )

# covariances
rng         = MersenneTwister(123);
SIGMApar    = InverseWishart((nanim*2-1)+5,Matrix{Float64}(I, nanim*2,nanim*2 ))
#Sigma = Vector{PDMat{Float64,Array{Float64,2}}}(UndefInitializer(),K)
Sigma = Vector(UndefInitializer(),K)
for k in 1:K
    Sigma[k] = PDMat(rand(rng,SIGMApar))
end

# means
vec0 = Vector{Float64}(undef,nanim*2)
for i in 1:(nanim*2)
    vec0[i] = 0.0::Float64
end

# ProbMatrix
ProbMat = Vector{Vector{Float64}}(undef,K)
ProbMat[1] = [0.8,0.1,0.1]
ProbMat[2] = [0.2,0.6,0.2]
ProbMat[3] = [0.15,0.05,0.8]

# Psi
Psi = Vector{Float64}([0.1,0.6,0.2])
#### Data similation
DataCoords              = Matrix{Float64}(undef,nt,nanim*2)  # data
LatentClassification    = Vector{Int64}(undef,nt)



# latent classification
LatentClassification[1] = 1
for i in 2:nt
    z = LatentClassification[i-1]
    LatentClassification[i] = rand(rng,Categorical(ProbMat[z]))
end
# the data
for i in 1:(nanim*2)
    DataCoords[1,i] = 0.0
end
for i in 2:nt
    z = LatentClassification[i-1]
    DataCoords[i,:] = mu0[z]+Psi[z]*(DataCoords[i-1,:]-mu0[z]  ) +  rand(rng,MultivariateNormal(vec0,Sigma[z]))
end




#### The model

nt          = 1500
kmax        = 3
nanim       = 2
nc          = nanim*2

InitSigma   = Vector{Matrix{Float64}}()
for k in 1:kmax
    push!(InitSigma ,Matrix{Float64}(I,nc,nc))
end



DataCoords2 = Matrix{Union{Float64,Missing} }(DataCoords)

MCMCLikelihood = OptionsLikelihood(
    # numer of regimes
    kmax                = Int16(kmax),
    # data
    data                = DataCoords2,
    # model type
    likelihood_type     = "OU",

    # Missing
    update_missing      = true,
    # mu0
    update_mu0          = true,
    inits_mu0           = [zeros(nc,kmax)[:,i] for i in 1:size(zeros(nc,kmax),2)],
    prior_mu0           = Dict("name"=>"MvNormal","μ"=>zeros(nc), "Σ"=>Matrix{Float64}(I,nc,nc)),
    # psi
    update_psi          = true,
    inits_psi           = [zeros(Integer(nanim),kmax)[:,i] for i in 1:size(zeros(Integer(nanim),kmax),2)],
    prior_psi           = Dict("name"=>"MvUniform","a"=>zeros(nanim), "b"=>ones(nanim)),
    # Sigma
    update_sigma        = true,
    inits_sigma         = InitSigma,
    prior_sigma         = Dict("name"=>"InverseWishart","df"=>Float64(nc+1), "Ψ"=>Matrix{Float64}(I,nc,nc)),
    # zeta
    update_zeta         = true,
    inits_zeta          = Int16.(sample(1:kmax, nt, replace = true))
)

InitPi   = Vector{Vector{Float64}}()
for k in 1:kmax
    push!(InitPi ,ones(kmax)/3.0)
end
MCMCClusterization = OptionsClusterization(
    Likelihood =           MCMCLikelihood,
    clustering_type     = "HMM"::String, # one of "HMM"
    # pi
    update_pi           = true::Bool,
    inits_pi            = InitPi,
    prior_pi            = Dict("name"=>"Dirichlet","alpha"=>(ones(kmax)/kmax))
    )

MCMCout = OptionsMCMC(MCMCLikelihood,MCMCClusterization;
        iterations = Int64(10000),
        burnin = Int64(5000),
        thin = Int64(2),

        )


using BayesianAnimalMovementModels
TODELETE = BayesianAnimalMovementModels
# using Profile
# using ProfileView
# using Traceur
timings = zeros(6)
appBurninOrThin = 500

appBurninOrThin = MCMCout.burnin
for saveMCMCindex in 1:MCMCout.nsamplesave

    for burnithinMCMCindex in 1:appBurninOrThin

        ### ### ### ### ### ### ###
        ### ### LIKELIHOOD
        ### ### ### ### ### ### ###

        ### missing
        timings[1] += @elapsed TODELETE.sample_missing!(MCMCLikelihood)

        # @trace      TODELETE.sample_missing!(MCMCLikelihood)
        # @profview   TODELETE.sample_missing!(MCMCLikelihood)
        # @time       TODELETE.sample_missing!(MCMCLikelihood)

        ### mu0
        timings[2] += @elapsed TODELETE.sample_mu0!(MCMCLikelihood,MCMCLikelihood.mu0)

        # @trace      TODELETE.sample_mu0!(MCMCLikelihood,MCMCLikelihood.mu0)
        # @profile       TODELETE.sample_mu0!(MCMCLikelihood,MCMCLikelihood.mu0);Profile.print()
        # @time       TODELETE.sample_mu0!(MCMCLikelihood,MCMCLikelihood.mu0)

        ### sigma
        timings[3] += @elapsed TODELETE.sample_sigma!(MCMCLikelihood, MCMCLikelihood.sigma)

        # @trace      TODELETE.sample_sigma!(MCMCLikelihood, MCMCLikelihood.sigma)
        # @profview   TODELETE.sample_sigma!(MCMCLikelihood, MCMCLikelihood.sigma)
        # @time       TODELETE.sample_sigma!(MCMCLikelihood, MCMCLikelihood.sigma)

        ### psi
        timings[4] += @elapsed TODELETE.sample_psi!(MCMCLikelihood, MCMCLikelihood.psi)

        # @trace      TODELETE.sample_psi!(MCMCLikelihood, MCMCLikelihood.psi)
        # @profview   TODELETE.sample_psi!(MCMCLikelihood, MCMCLikelihood.psi)
        # @time       TODELETE.sample_psi!(MCMCLikelihood, MCMCLikelihood.psi)

        ### ### ### ### ### ### ###
        ### ### ZETA
        ### ### ### ### ### ### ###

        # zeta
        timings[5] += @elapsed TODELETE.sample_zeta!(MCMCLikelihood, MCMCClusterization)

        # @trace      TODELETE.sample_zeta!(MCMCLikelihood, MCMCClusterization)
        # @profview   TODELETE.sample_zeta!(MCMCLikelihood, MCMCClusterization)
        # @time       TODELETE.sample_zeta!(MCMCLikelihood, MCMCClusterization)

        ### ### ### ### ### ### ###
        ### ### CLUSTERIZATION
        ### ### ### ### ### ### ###

        ### pi
        timings[6] += @elapsed TODELETE.sample_pi!(MCMCClusterization, MCMCClusterization.pi)


    end # second for MCMC
    appBurninOrThin = MCMCindices.thin

    TODELETE.save_posteriorsamples!(MCMCout,MCMCLikelihood, MCMCClusterization)

end
