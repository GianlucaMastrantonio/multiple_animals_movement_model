
### Parameter to set
DIR_CODE = "/Users/gianlucamastrantonio/Dropbox (Politecnico di Torino Staff)/lavori/HDP_Animal/ToSumbit/codes/multiple_animals_movement_model/"

#### #### #### #### #### #### ####
#### #### Code
#### #### #### #### #### #### ####

### Packages
using Revise
using HierarchicalMultivariateAnimalMovement
using Distributions, Random
using LinearAlgebra, PDMats
using RCall

imod = 1
IndAn = 1:12
rngseed = 310

rng    	= MersenneTwister(rngseed);
Random.seed!(rngseed)# 1234 1

##### DIRECTORIES and NAMES
NAME 	= ["Standard", "Prop", "singleHDP"][imod]
NAME    = string(NAME,"_",rngseed,"_1_1_")

DIR  	= "/Users/gianlucamastrantonio/Dropbox (Politecnico di Torino Staff)/lavori/HDP_Animal/ToSumbit/codes/multiple_animals_movement_model/"
DIRDATA = string(DIR,"data/")
DIROUT  = string(DIR,"output/")

#### LOAD DATA FROM R
@rput DIRDATA
@rput DIROUT
@rput IndAn
R" load(paste(DIRDATA,'Data2.Rdata'  ,sep='')) "
@rget DataAn2

R"""

	#DataCoords = DataAn3[,1:14]
	DataCoords = DataAn2[,IndAn]
	print(c(IndAn[1],IndAn[2]))
	SDtot1   = sd(DataCoords[,seq(1,ncol(DataCoords),by=2 )],na.rm=T)
	SDtot2   = sd(DataCoords[,seq(2,ncol(DataCoords),by=2 )],na.rm=T)

	Meantot1 = mean(DataCoords[,seq(1,ncol(DataCoords),by=2 )],na.rm=T)
	Meantot2 = mean(DataCoords[,seq(2,ncol(DataCoords),by=2 )],na.rm=T)

	SDtot    = mean(c(SDtot1,SDtot2))
	Meantot  = mean(c(Meantot1,Meantot2))

	DataCoords[,seq(1,ncol(DataCoords),by=2 )] = (DataCoords[,seq(1,ncol(DataCoords),by=2 )]-Meantot1)/SDtot
	DataCoords[,seq(2,ncol(DataCoords),by=2 )] = (DataCoords[,seq(2,ncol(DataCoords),by=2 )]-Meantot2)/SDtot



    nCov = (ncol(DataCoords))
    nanim = nCov/2

	set.seed(1)


"""
@rget DataCoords

NANindex 				= isnan.(DataCoords)

Dataset = DataCoords



###
nanim = Int16(size(Dataset,2)/2)
DataVec = Vector{Matrix{Union{Float64,Missing}}}(undef,nanim )
for ian in 1:nanim
	DataVec[ian] = Matrix{Union{Float64,Missing}}(Dataset[:,  [1+(ian-1)*2,2+(ian-1)*2] ])
	NANindex 	 = isnan.(DataVec[ian][:,1])
	# DataVec[ian][NANindex,1]   .= missing
	# DataVec[ian][NANindex,2]   .= missing
	# DataVec[ian][IndexNA,:]  .= missing
end


#### The model

kmax        = 100
NAMETOT     = string(NAME,"K=",kmax)
nt          = size(Dataset,1)

nc          = size(Dataset,2)
nanim       = nanim


InitSigma   = Vector{Vector{Matrix{Float64}}}(undef, nanim)
for ian in 1:nanim
	InitSigma[ian] = Vector{Matrix{Float64}}()
	for k in 1:kmax
		push!(InitSigma[ian], Matrix{Float64}(I,2,2))
	end
end

InitMu  = Vector{Vector{Vector{Float64}}}(undef, nanim)
InitNu  = Vector{Vector{Vector{Float64}}}(undef, nanim)
InitEta = Vector{Vector{Vector{Float64}}}(undef, nanim)
InitRho = Vector{Vector{Vector{Float64}}}(undef, nanim)
InitZeta = Vector{Vector{Int16}}(undef, nanim)
InitAngle = Vector{Vector{Vector{Float64}}}(undef, nanim)
for ian in 1:nanim
	InitMu[ian]  = [zeros(2,kmax)[:,i] for i in 1:size(zeros(2,kmax),2)]
	InitNu[ian]  = [0.5*ones(Integer(1),kmax)[:,i] for i in 1:size(zeros(Integer(1),kmax),2) ]
	InitEta[ian] = 	[zeros(2,kmax)[:,i] for i in 1:size(zeros(2,kmax),2)]
	InitRho[ian] = [ones(Integer(1),kmax)[:,i].*0.5 for i in 1:size(ones(Integer(1),kmax),2) ]
	InitZeta[ian] = Int16.(sample(1:min(kmax,15), nt, replace = true))
	InitAngle[ian] = [zeros(Integer(1),kmax)[:,i] for i in 1:size(zeros(Integer(1),kmax),2) ]
end


NAMETOTapp = deepcopy(NAMETOT)



global NAMETOTapp
NAMETOT     = string(NAMETOTapp)
MCMCLikelihood = OptionsLikelihood(
	# numer of regimes
	kmax                = Int16(kmax),
	# data
	data                = DataVec::Vector{Matrix{Union{Float64,Missing}}},

	# Missing
    update_missing      = true::Bool,

	# mu0

    update_mu          = true::Bool,
    inits_mu           = InitMu,
    prior_mu           = Dict("name"=>"MvNormal","μ"=>zeros(2), "Σ"=>Matrix{Float64}(I,2,2).*100)::Dict,

    # nu

    update_nu          = true::Bool,
    inits_nu           = InitNu,
    prior_nu           = Dict("name"=>"MvUniform","a"=>zeros(1), "b"=>ones(1))::Dict,

	# eta
    update_eta          = true::Bool,
    inits_eta           = InitEta,
    prior_eta           = Dict("name"=>"MvNormal","μ"=>zeros(2), "Σ"=>Matrix{Float64}(I,2,2).*100)::Dict,

    # rho

	update_rho          = true::Bool,
    inits_rho           = InitRho,
    prior_rho           = Dict("name"=>"MvUniform","a"=>zeros(1), "b"=>ones(1))::Dict,

    # Sigma

    update_sigma        = true::Bool,
    inits_sigma         = InitSigma::Vector{Vector{Matrix{Float64}}},
    prior_sigma         = Dict("name"=>"InverseWishart","df"=>Float64(2+1), "Ψ"=>Matrix{Float64}(I,2,2))::Dict,

	# angle
	inits_angle 		= InitAngle::Vector{Vector{Vector{Float64}}},
	prior_angle 		= Dict("name"=>"MvUniform","a"=> -pi*ones(Float64,Integer(1)), "b"=> pi*ones(Float64,Integer(1))),

	#zeta

    update_zeta         = true::Bool,
    inits_zeta          = InitZeta

)


InitPi   = Vector{Vector{Vector{Float64}}}(undef, nanim)
for ian in 1:nanim
	InitPi[ian]   = Vector{Vector{Float64}}()
	for k in 1:kmax
	    push!(InitPi[ian] ,ones(kmax)/3.0)
	end
end


MCMCClusterization = OptionsClusterization(
	nanim      = Int64(nanim),
    Likelihood =           MCMCLikelihood,
    clustering_type     = "HDP-HMM-Divided"::String, # one of "HMM"

    # pi
    update_pi           = true::Bool,
    inits_pi            = InitPi,
    prior_pi            = Dict("name"=>"Dirichlet","alpha"=>(ones(kmax)/kmax)),

	prior_ak 			= [1.1; 0.1],
    prior_gamma 		= [1.1; 0.1],
    prior_rho 			= [1.0; 1.0]
    )









HierarchicalParameters = OptionHierarchicalParameters(MCMCLikelihood)

IterMolt = Int16(15)
MCMCout = Vector{MCMCutils}(undef,nanim)
for ian in 1:nanim
	MCMCout[ian] = OptionsMCMC(
		MCMCLikelihood[ian],
		MCMCClusterization[ian],
		HierarchicalParameters;
		iterations = Int64(5000*IterMolt),
		burnin = Int64(2500*IterMolt),
		thin = Int64(1*IterMolt)

	        )
end

#imod = 1
if imod == 1
	ModelOut  = MCMCalgorithm(
		MCMCout,
		MCMCLikelihood,
	    MCMCClusterization,
		rng
	)
end

if imod == 2
	ModelOut  = MCMCalgorithm(
		MCMCout,
		MCMCLikelihood,
	    MCMCClusterization,
		HierarchicalParameters,
		rng;
		Iter_lev2 = Int32(0)
	)
end

if imod == 3
	ModelOut  = MCMCalgorithm_singleHDP(
		MCMCout,
		MCMCLikelihood,
	    MCMCClusterization,
		HierarchicalParameters,
		rng;

end






NAMEOUT = NAMETOT
R"""
	library(GenFunc)
	library(lattice)
	MCMCMODEL_OUT = list()
"""
for ian in 1:nanim

	ModelOUT = ModelOut["PosteriorSamples"][ian].mcmc_out
	@rput ModelOUT
	@rput NAMEOUT
	@rput ian
	R"""
		MCMCMODEL_OUT[[ian]] = ModelOUT

	"""
end

@rput DataVec
R"""
	save.image(paste(DIROUT,"MODELOUT","_",NAMEOUT,".Rdata",sep=""))
"""
