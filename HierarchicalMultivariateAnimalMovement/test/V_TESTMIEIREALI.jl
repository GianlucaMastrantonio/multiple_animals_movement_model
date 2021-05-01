##### Packages
using Revise
using BayesianAnimalMovementModels
using Distributions, Random
using LinearAlgebra, PDMats
using Test
using RCall
##### DIRECTORIES and NAMES
NAME 	= "RealDataFirstModel"

DIR  	= "/Users/gianlucamastrantonio/Dropbox (Politecnico di Torino Staff)/lavori/sheepdog/"
DIRDATA = string(DIR,"Data/")
DIRCODE = string(DIR,"analisi/")
DIROUT  = string(DIR,"analisi/plot/")

##### JULIA PARAMETERS
#rng     = MersenneTwister(123);

#### LOAD DATA FROM R

@rput DIRDATA
@rput DIROUT
R" load(paste(DIRDATA,'Data3.Rdata'  ,sep='')) "
@rget DataAn3

R"""

	DataCoords = DataAn3[,1:14]

	SDtot1 = sd(DataCoords[,seq(1,ncol(DataCoords),by=2 )],na.rm=T)
	SDtot2 = sd(DataCoords[,seq(2,ncol(DataCoords),by=2 )],na.rm=T)

	Meantot1 = mean(DataCoords[,seq(1,ncol(DataCoords),by=2 )],na.rm=T)
	Meantot2 = mean(DataCoords[,seq(2,ncol(DataCoords),by=2 )],na.rm=T)

	DataCoords[,seq(1,ncol(DataCoords),by=2 )] = (DataCoords[,seq(1,ncol(DataCoords),by=2 )]-Meantot1)/SDtot1
	DataCoords[,seq(2,ncol(DataCoords),by=2 )] = (DataCoords[,seq(2,ncol(DataCoords),by=2 )]-Meantot2)/SDtot2


    # pdf(paste(DIROUT,"Data.pdf",sep=""))
    # nCov = (ncol(DataCoords))
    # nanim = nCov/2
    # for(i in 1:nanim)
    # {
    #     plot(range(c(  DataCoords[,seq(1,nCov, by=2)] ),na.rm =T),range(c(  DataCoords[,1+seq(1,nCov, by=2)] ),na.rm =T), type="n")
    #     points(DataCoords[,1:2+(i-1)*2], type="l")
    #     points(DataCoords[,1:2+(i-1)*2], pch=20, cex=1)
	#
    # }
    # par(mfrow=c(3,3))
    # for(i in 1:(nCov-1))
    # {
    #     for(j in (i+1):nCov)
    #     {
    #         plot(DataCoords[,c(i,j)], pch=20, cex=1)
    #     }
    # }
    # dev.off()
"""


@rget DataCoords
NANindex 				= isnan.(DataCoords)
DataCoords  			= Matrix{Union{Float64,Missing}}(DataCoords)
DataCoords[NANindex]   .= missing

Dataset = DataCoords

#### The model
kmax        = 7
NAMETOT     = string(NAME,"K=",kmax)
nt          = size(Dataset,1)

nc          = size(Dataset,2)
nanim       = Int32(nc/2)


InitSigma   = Vector{Matrix{Float64}}()
for k in 1:kmax
    push!(InitSigma ,Matrix{Float64}(I,nc,nc))
end
rand(InverseWishart(nc*3.0,Matrix{Float64}(I,nc,nc)*25))

#DataCoords2 = Matrix{Union{Float64,Missing} }(DataCoords)

MCMCLikelihood = OptionsLikelihood(
    # numer of regimes
    kmax                = Int16(kmax),
    # data
    data                = Dataset,
    # model type
    likelihood_type     = "OU",

    # Missing
    update_missing      = true,
    # mu0
    update_mu0          = true,
    inits_mu0           = [rand(nc,kmax)[:,i] for i in 1:kmax],
    prior_mu0           = Dict("name"=>"MvNormal","μ"=>zeros(nc), "Σ"=>Matrix{Float64}(I,nc,nc)*5),
    # psi
    update_psi          = true,
    inits_psi           = [0.5*ones(Integer(nanim),kmax)[:,i] for i in 1:kmax],
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
    push!(InitPi ,ones(kmax)/kmax)
end

MCMCClusterization = OptionsClusterization(
    Likelihood =           MCMCLikelihood,
    clustering_type     = "HMM"::String, # one of "HMM"
    # pi
    update_pi           = true::Bool,
    inits_pi            = InitPi,
    prior_pi            = Dict("name"=>"Dirichlet","alpha"=>ones(kmax))
    )

MCMCout = OptionsMCMC(MCMCLikelihood,MCMCClusterization;
        iterations = Int64(1000),
        burnin = Int64(500),
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

		MCMCLikelihood.mu0.prior.Σ.mat
		MCMCLikelihood.mu0.prior.μ
		MCMCLikelihood.mu0.parameteracc
		MCMCLikelihood.data.data
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
