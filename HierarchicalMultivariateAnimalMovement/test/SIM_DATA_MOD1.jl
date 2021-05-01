##### Packages
using Revise
using HierarchicalMultivariateAnimalMovement
using Distributions, Random
using LinearAlgebra, PDMats
using Test
TODELETE = HierarchicalMultivariateAnimalMovement
##### JULIA PARAMETERS
rng     = MersenneTwister(123);
##### Simulation Parameters

nt = 1500    # number observations
K  = 3       # number total groups
nanim = 2    # numner of animales

# coordinates of mu0 (K,Animale, (x,y))
mu = Vector{Vector{Vector}}(undef, nanim)
mu[1] = [ [0.0; 0.0], [10.0 ; 10.0] , [-10.0; 10.0] ]
mu[2] = [ [0.0;  0.0], [ 0.0; 0.0] , [10.0;  10.0] ]

eta = Vector{Vector{Vector}}(undef, nanim)
eta[1] = [ [0.0; 1.0], [1.0 ; 1.0] , [0.0; 0.0] ]
eta[2] = [ [0.0; 1.0], [0.0; -1.0] , [-1.0;  -1.0] ]

rho = Vector{Vector{Vector}}(undef, nanim)
rho[1] = [ [0.0], [0.5] , [1.0] ]
rho[2] = [ [0.5], [0.7] , [0.7] ]

nu = Vector{Vector{Vector}}(undef, nanim)
nu[1] = [ [0.9], [0.7] , [0.5] ]
nu[2] = [ [0.1], [0.5] , [0.7] ]

Sigma = Vector{Vector{PDMat{Float64,Array{Float64,2}}}}(undef, nanim)
Sigma[1] = [PDMat([	1.0 0;
					0 2.0]);
			PDMat([	2.0  0.8;
					0.8 1.0]);
			PDMat([	2.0  0.8;
					0.8 1.0])
			]
Sigma[2] = [PDMat([	2.0 0.5;
					0.5 2.0]);
			PDMat([	1.0  0.5;
					0.5 2.0]);
			PDMat([	2.0  0.8;
					0.8 1.0])
			]

# means
vec0 = Vector{Float64}(undef,2)
for i in 1:(2)
    vec0[i] = 0.0::Float64
end

# ProbMatrix
ProbMat = Vector{Matrix{Float64}}(undef,nanim)
ProbMat[1]  = [  [0.8 0.1 0.1];
			    [0.2 0.6 0.2];
				[0.15 0.05 0.8]
			  ]

ProbMat[2]  = [ [0.4 0.3 0.3];
				[0.1 0.7 0.2];
				[0.1 0.3 0.6]
			  ]


#### Data similation
DataCoords				 = Vector{Matrix{Float64}}(undef,nanim)
DataCoords[1]            = Matrix{Float64}(undef,nt,2)
DataCoords[2]            = Matrix{Float64}(undef,nt,2)

LatentClassification     = Vector{Vector{Int64}}(undef, nanim)
LatentClassification[1]  = Vector{Int64}(undef,nt)
LatentClassification[2]  = Vector{Int64}(undef,nt)



# latent classification
for ian in 1:nanim
	LatentClassification[ian][1] = 1
	for i in 2:nt
	    z = LatentClassification[ian][i-1]
	    LatentClassification[ian][i] = rand(rng,Categorical(ProbMat[ian][z,:]))
	end
end

# the data
for ian in 1:nanim
	DataCoords[ian][1,1] = 0.0
	DataCoords[ian][1,2] = 0.0
end


Angle = [0.0]
for ian in 1:nanim
	global Angle
	Angle = [0.0]
	for i in 2:nt

	    z = LatentClassification[ian][i-1]

		R = [cos.(rho[ian][z][].*Angle[]) -sin.(rho[ian][z][].*Angle[]); sin.(rho[ian][z][].*Angle[]) cos.(rho[ian][z][].*Angle[])]

		DataCoords[ian][i,:] =  DataCoords[ian][i-1,:]+(1.0 .-rho[ian][z]) .* nu[ian][z] .*(mu[ian][z,]-DataCoords[ian][i-1,:])
		DataCoords[ian][i,:] += rho[ian][z][]*R*eta[ian][z,]
		DataCoords[ian][i,:]      += R*rand(rng,MultivariateNormal(vec0,Sigma[ian][z]))

		Angle = atan(DataCoords[ian][i,:][2]-DataCoords[ian][i-1,:][2],DataCoords[ian][i,:][1]-DataCoords[ian][i-1,:][1])
	end
end





#### The model

nt          = 1500
kmax        = 10
nanim       = 2
nc          = nanim*2

InitSigma   = Vector{Vector{Matrix{Float64}}}(undef, nanim)
for ian in 1:nanim
	InitSigma[ian] = Vector{Matrix{Float64}}()
	for k in 1:kmax
		push!(InitSigma[ian], Matrix{Float64}(I,2,2))
	end
end


DataCoords2 = Vector{Matrix{Union{Float64,Missing} }}(DataCoords)
#using BayesianAnimalMovementModels

DataCoords2[1][1,:] .= missing
DataCoords2[1][100,1:2] .= missing

DataCoords2[2][10,:] .= missing
DataCoords2[2][110,1:2] .= missing



##########
#MCMCLikelihood.clusterization.zeta
MCMCLikelihood = TODELETE.OptionsLikelihood(
    # numer of regimes
    kmax                = Int16(kmax)::Int16,
    # data
    data                = DataCoords2::Vector{Matrix{Union{Float64,Missing}}},

    # Missing
    update_missing      = true::Bool,
    # mu0
    update_mu          = true::Bool,
    inits_mu           = [	[zeros(2,kmax)[:,i] for i in 1:size(zeros(2,kmax),2)],
						   	[zeros(2,kmax)[:,i] for i in 1:size(zeros(2,kmax),2)]]::Vector{Vector{Vector{Float64}}},
    prior_mu           = Dict("name"=>"MvNormal","μ"=>zeros(2), "Σ"=>Matrix{Float64}(I,2,2))::Dict,
    # psi
    update_nu          = true::Bool,
    inits_nu           = [[zeros(Integer(1),kmax)[:,i] for i in 1:size(zeros(Integer(1),kmax),2)] ,
						  [zeros(Integer(1),kmax)[:,i] for i in 1:size(zeros(Integer(1),kmax),2) ] ]::Vector{Vector{Vector{Float64}}},
    prior_nu           = Dict("name"=>"MvUniform","a"=>zeros(1), "b"=>ones(1))::Dict,
	# eta
    update_eta          = true::Bool,
    inits_eta           = [	[zeros(2,kmax)[:,i] for i in 1:size(zeros(2,kmax),2)],
						   	[zeros(2,kmax)[:,i] for i in 1:size(zeros(2,kmax),2)]]::Vector{Vector{Vector{Float64}}},
    prior_eta           = Dict("name"=>"MvNormal","μ"=>zeros(2), "Σ"=>Matrix{Float64}(I,2,2))::Dict,
    # rho
    update_rho          = true::Bool,
    inits_rho           = [[ones(Integer(1),kmax)[:,i].*0.5 for i in 1:size(ones(Integer(1),kmax),2)] ,
						  [ones(Integer(1),kmax)[:,i].*0.5 for i in 1:size(ones(Integer(1),kmax),2) ] ]::Vector{Vector{Vector{Float64}}},
    prior_rho           = Dict("name"=>"MvUniform","a"=>zeros(1), "b"=>ones(1))::Dict,
    # Sigma
    update_sigma        = true::Bool,
    inits_sigma         = InitSigma::Vector{Vector{Matrix{Float64}}},
    prior_sigma         = Dict("name"=>"InverseWishart","df"=>Float64(2+1), "Ψ"=>Matrix{Float64}(I,2,2))::Dict,

	# angle
	inits_angle 		= [[zeros(Integer(1),kmax)[:,i] for i in 1:size(zeros(Integer(1),kmax),2)] ,
						  [zeros(Integer(1),kmax)[:,i] for i in 1:size(zeros(Integer(1),kmax),2) ] ]::Vector{Vector{Vector{Float64}}},
	prior_angle 		= Dict("name"=>"MvUniform","a"=> -pi*ones(Float64,Integer(1)), "b"=> pi*ones(Float64,Integer(1))),
	# zeta
    update_zeta         = true::Bool,
    inits_zeta          = [Int16.(sample(1:kmax, nt, replace = true)),
							Int16.(sample(1:kmax, nt, replace = true))]::Vector{Vector{Int16}}

)


InitPi   = Vector{Vector{Vector{Float64}}}(undef, nanim)
for ian in 1:nanim
	InitPi[ian]   = Vector{Vector{Float64}}()
	for k in 1:kmax
	    push!(InitPi[ian] ,ones(kmax)/3.0)
	end
end


MCMCClusterization = TODELETE.OptionsClusterization(
	nanim      = 2,
    Likelihood =           MCMCLikelihood,
    clustering_type     = "HDP-HMM-Divided"::String, # one of "HMM"
    # pi
    update_pi           = true::Bool,
    inits_pi            = InitPi,
    prior_pi            = Dict("name"=>"Dirichlet","alpha"=>(ones(kmax)/kmax))
    )

MCMCout = Vector{TODELETE.MCMCutils}(undef,nanim)
for ian in 1:nanim
	MCMCout[ian] = TODELETE.OptionsMCMC(MCMCLikelihood[ian],MCMCClusterization[ian];
	        iterations = Int64(10000),
	        burnin = Int64(5000),
	        thin = Int64(2),

	        )
end



using BayesianAnimalMovementModels


timings = zeros(6)
appBurninOrThin = 500

appBurninOrThin = MCMCout[1].burnin
for saveMCMCindex in 1:MCMCout.nsamplesave

    for burnithinMCMCindex in 1:appBurninOrThin

        ### ### ### ### ### ### ###
        ### ### LIKELIHOOD
        ### ### ### ### ### ### ###

        ### missing
        TODELETE.sample_missing!(MCMCLikelihood,MCMCLikelihood[1].miss)

        ### mu0 (mu)
		TODELETE.sample_mu0!(MCMCLikelihood,MCMCLikelihood[1].mu)

        ### sigma
        TODELETE.sample_sigma!(MCMCLikelihood, MCMCLikelihood[1].sigma)

        ### psi
        TODELETE.sample_psi!(MCMCLikelihood, MCMCLikelihood[1].nu)

		# eta
		TODELETE.sample_muC!(MCMCLikelihood,MCMCLikelihood[1].eta)

		# rho
		TODELETE.sample_rho!(MCMCLikelihood,MCMCLikelihood[1].rho)

        ### ### ### ### ### ### ###
        ### ### ZETA
        ### ### ### ### ### ### ###

        # zeta
        TODELETE.sample_zeta!(MCMCLikelihood, MCMCClusterization)


        ### ### ### ### ### ### ###
        ### ### CLUSTERIZATION
        ### ### ### ### ### ### ###

        ### pi
        TODELETE.sample_pi!(MCMCClusterization, MCMCClusterization[1].pi)

		TODELETE.compute_m!(MCMCClusterization, MCMCClusterization[1].pi)
		TODELETE.sample_beta!(MCMCClusterization, MCMCClusterization[1].pi)
		TODELETE.sample_ak!(MCMCClusterization, MCMCClusterization[1].pi)
		TODELETE.sample_gamma!(MCMCClusterization, MCMCClusterization[1].pi)
		TODELETE.sample_rhodp!(MCMCClusterization, MCMCClusterization[1].pi)


    end # second for MCMC
    appBurninOrThin = MCMCindices.thin

    TODELETE.save_posteriorsamples!(MCMCout,MCMCLikelihood, MCMCClusterization)

end
