
#### #### #### #### #### #### #### #### #### #### #### ####
#### #### #### #### #### LIKELIHOOD
#### #### #### #### #### #### #### #### #### #### #### ####
function OptionsLikelihood(;
    # numer of regimes
    kmax::Int16,
    # data
    data::Vector{Matrix{Union{Float64,Missing}}},

    # Missing
    update_missing::Bool,
    # mu0
    update_mu::Bool,
    inits_mu::Vector{Vector{Vector{Float64}}},
    prior_mu::Dict,
    # psi
    update_nu::Bool,
    inits_nu::Vector{Vector{Vector{Float64}}},
    prior_nu::Dict,
	# eta
    update_eta::Bool,
    inits_eta::Vector{Vector{Vector{Float64}}},
    prior_eta::Dict,
    # rho
    update_rho::Bool,
    inits_rho::Vector{Vector{Vector{Float64}}},
    prior_rho::Dict,
    # Sigma
    update_sigma::Bool,
    inits_sigma::Vector{Vector{Matrix{Float64}}},
    prior_sigma::Dict,

	# angle
	inits_angle::Vector{Vector{Vector{Float64}}},
	prior_angle::Dict,
	# zeta
    update_zeta::Bool,
    inits_zeta::Vector{Vector{Int16}},
)

###### missing
datacopy     = deepcopy(data)
nanim 		 = size(datacopy)[1]
nt           = size(datacopy[1])[1]
ncol         = size(datacopy[1])[2]

mcmc_missing = Vector{MissingDoUpdate}(undef, nanim)
mcmc_data    = Vector{CoordinatesDataset}(undef, nanim)
mcmc_zeta    = Vector{AbstractZeta}(undef, nanim)
mcmc_mu      = Vector{AbstractVecPar}(undef, nanim)
mcmc_eta     = Vector{AbstractVecPar}(undef, nanim)
mcmc_nu      = Vector{AbstractVecPar}(undef, nanim)
mcmc_rho     = Vector{AbstractVecPar}(undef, nanim)
mcmc_sigma   = Vector{AbstractPosDefMatPar}(undef, nanim)
mcmc_angle   = Vector{AbstractVecPar}(undef, nanim)
Likelihood   = Vector{Likelihood_OU_CircLinmodel}(undef, nanim)
for ian in 1:nanim

	app_missing = size(findall(ismissing, datacopy[ian]))[1] == 0

	if !update_missing
		mcmc_missing[ian] = MissingDoNotUpdate(datacopy[ian])
		if app_missing

		else
			printstyled("ATTENTION: missing values have been interpolated and there will be not estimated \n", color=:blue)
		end
		#mcmc_missing = MissingDoNotUpdate(datacopy)
	else
		mcmc_missing[ian] = MissingDoUpdate(datacopy[ian])
	end

	# dataset
	mcmc_data[ian]   = CoordinatesDataset(datacopy[ian]::Matrix{Union{Missing,Float64}})

	# zeta
	if size(inits_zeta[ian],1)!=nt error("length inits_zeta must be ", nt) end
	if update_zeta==false
		mcmc_zeta[ian] = ZetaDoNotUpdate(inits_zeta[ian],kmax)
	else
		mcmc_zeta[ian] = ZetaDoUpdate(inits_zeta[ian],kmax)
	end

	# data model
	mcmc_mu[ian]    = VecParDoNotUpdate()
	mcmc_rho[ian]   = VecParDoNotUpdate()

	# mu
	if update_mu==false
		mcmc_mu[ian] = VecParDoNotUpdate()
	else
		checkPriorName  = size(findall(["MvNormal"] .== prior_mu["name"]),1)== 0
		if checkPriorName error("mu prior must be one of MvNormal") end

		namecheck       = string("check_Mixture",prior_mu["name"])
		nameprior       = string(prior_mu["name"],"_mu")
		checkfunc       = getfield(HierarchicalMultivariateAnimalMovement,Symbol(namecheck))
		p_aftercheck    = checkfunc(nameprior,inits_mu[ian], prior_mu, kmax, ncol)

		structobject    = getfield(HierarchicalMultivariateAnimalMovement,Symbol(string("VecPar",  prior_mu["name"])))
		mcmc_mu[ian]    = structobject(p_aftercheck)

	end

	# eta
	if update_eta==false
		mcmc_eta[ian] = VecParDoNotUpdate()
	else
		checkPriorName  = size(findall(["MvNormal"] .== prior_eta["name"]),1)== 0
		if checkPriorName error("eta prior must be one of MvNormal") end

		namecheck       = string("check_Mixture",prior_eta["name"])
		nameprior       = string(prior_eta["name"],"_eta")
		checkfunc       = getfield(HierarchicalMultivariateAnimalMovement,Symbol(namecheck))
		p_aftercheck    = checkfunc(nameprior,inits_eta[ian], prior_eta, kmax, ncol)

		structobject    = getfield(HierarchicalMultivariateAnimalMovement,Symbol(string("VecPar",  prior_eta["name"])))
		mcmc_eta[ian]        = structobject(p_aftercheck)

	end


	# nu
	if update_nu==false
		mcmc_nu[ian] = VecParDoNotUpdate()
	else
		checkPriorName  = size(findall(["MvBeta_NotSupported_Anymore","MvUniform"] .== prior_nu["name"]),1)== 0
		if checkPriorName error("nu prior must be one of MvBeta") end

		namecheck       = string("check_Mixture",prior_nu["name"])
		nameprior       = string(prior_nu["name"],"_nu")
		checkfunc       = getfield(HierarchicalMultivariateAnimalMovement,Symbol(namecheck))
		p_aftercheck    = checkfunc(nameprior,inits_nu[ian], prior_nu, kmax, Integer(ncol/2))

		structobject    = getfield(HierarchicalMultivariateAnimalMovement,Symbol(string("VecPar",  prior_nu["name"])))
		mcmc_nu[ian]    = structobject(p_aftercheck)

	end

	# rho
	if update_rho==false

		namecheck       = string("check_Mixture",prior_rho["name"])
		nameprior       = string("MvUniform","_rho")
		checkfunc       = getfield(HierarchicalMultivariateAnimalMovement,Symbol(namecheck))
		p_aftercheck    = checkfunc(nameprior,inits_rho[ian], prior_rho, kmax, Integer(ncol/2))
		#print(p_aftercheck[3])

		mcmc_rho[ian]        = VecParMvNoPrior(p_aftercheck[3])

	else
		checkPriorName  = size(findall(["MvBeta_NotSupported_Anymore","MvUniform"] .== prior_rho["name"]),1)== 0
		if checkPriorName error("rho prior must be one of MvBeta") end

		namecheck       = string("check_Mixture",prior_rho["name"])
		nameprior       = string(prior_rho["name"],"_rho")
		checkfunc       = getfield(HierarchicalMultivariateAnimalMovement,Symbol(namecheck))
		p_aftercheck    = checkfunc(nameprior,inits_rho[ian], prior_rho, kmax, Integer(ncol/2))

		structobject    = getfield(HierarchicalMultivariateAnimalMovement,Symbol(string("VecPar",  prior_rho["name"])))
		mcmc_rho[ian]   = structobject(p_aftercheck)

	end

	# Sigma
	if update_sigma==false
		mcmc_sigma[ian] = PosDefMatParDoNotUpdate()
	else

		checkPriorName  = size(findall(["InverseWishart"] .== prior_sigma["name"]),1)== 0
		if checkPriorName error("sigma prior must be one of InverseWishart") end

		namecheck       = string("check_Mixture",prior_sigma["name"])
		nameprior       = string(prior_sigma["name"],"_sigma")
		checkfunc       = getfield(HierarchicalMultivariateAnimalMovement,Symbol(namecheck))
		p_aftercheck    = checkfunc(nameprior,inits_sigma[ian], prior_sigma, kmax, Integer(ncol))

		structobject    = getfield(HierarchicalMultivariateAnimalMovement,Symbol(string("PosDefMat",  prior_sigma["name"])))
		mcmc_sigma[ian]        = structobject(p_aftercheck)

	end

	# angle
	update_angle = true

	# if update_angle==false
	# 	mcmc_angle[ian] = VecParDoNotUpdate()
	# else

		checkPriorName  = size(findall(["MvBeta_NotSupported_Anymore","MvUniform"] .== prior_angle["name"]),1)== 0
		if checkPriorName error("angle prior must be one of MvBeta") end

		namecheck       = string("check_Mixture",prior_angle["name"])
		nameprior       = string(prior_angle["name"],"_angle")
		checkfunc       = getfield(HierarchicalMultivariateAnimalMovement,Symbol(namecheck))
		p_aftercheck    = checkfunc(nameprior,inits_angle[ian], prior_angle, kmax, Integer(ncol/2))

		structobject    = getfield(HierarchicalMultivariateAnimalMovement,Symbol(string("VecPar",  prior_angle["name"])))
		mcmc_angle[ian]      = structobject(p_aftercheck)

	# end

	Likelihood[ian] = Likelihood_OU_CircLinmodel(mcmc_data[ian],mcmc_missing[ian],mcmc_mu[ian],mcmc_nu[ian],mcmc_sigma[ian],mcmc_zeta[ian],kmax,mcmc_eta[ian],mcmc_rho[ian],mcmc_angle[ian], true)
end
datacopy = nothing
return Likelihood

end

#### #### #### #### #### #### #### #### #### #### #### ####
#### #### #### #### #### LATENT CLASSIFICATION
#### #### #### #### #### #### #### #### #### #### #### ####
function OptionsClusterization(;
	nanim::Int,
	Likelihood::Vector{Tl},
    clustering_type::String, # one of HDP-HMM-Divided, HDP-HMM-Joint
    # pi
    update_pi::Bool,
    inits_pi::Vector{Vector{Vector{Float64}}},
    prior_pi::Dict,
    prior_ak = [1.0; 1.0]::Vector{Float64},
    prior_gamma = [1.0; 1.0]::Vector{Float64},
    prior_rho = [1.0 ;1.0]::Vector{Float64}
    ) where {Tl<:AbstractLikelihood}

    kmax = Likelihood[1].kmax

	mcmc_pi 		         = Vector{AbstractVecPar}(undef, nanim)
	mcmc_initpi              = Vector{typeof(1.0/(Float64(kmax))*ones(Float64,kmax))}(undef,nanim)
	mcmc_initz               = Vector{typeof(ones(Int16,1))}(undef,nanim)
	mcmc_beta                = Vector{typeof(1.0/kmax*ones(Float64,kmax))}(undef,nanim)
	mcmc_rho                 = Vector{typeof([0.5])}(undef,nanim)
	mcmc_gamma               = Vector{typeof([1.0])}(undef,nanim)
	mcmc_ak                  = Vector{typeof([1.0])}(undef,nanim)
	#Clusterization			 = Vector{AbstractClusterization}(undef,nanim)

	if clustering_type == "HDP-HMM-Divided"
		Clusterization = Vector{Clusterization_HDPHMM}(undef,nanim)
	elseif clustering_type == "HDP-HMM-Joint"
		Clusterization = Vector{Clusterization_HDPHMM_Joint}(undef,nanim)
	end

	for ian in 1:nanim
		if clustering_type == "HDP-HMM-Divided"
	        # pi
	        if update_pi==false
	            mcmc_pi[ian] = PiDoNotUpdate()
	        else

	            checkPriorName  = size(findall(["Dirichlet"] .== prior_pi["name"]),1)== 0
	            if checkPriorName error("pi prior must be one of Dirichlet") end

	            namecheck       = string("check_Clusterization",prior_pi["name"])
	            nameprior       = string(prior_pi["name"],"_pi")
	            checkfunc       = getfield(HierarchicalMultivariateAnimalMovement,Symbol(namecheck))
	            p_aftercheck    = checkfunc(nameprior,inits_pi[ian], prior_pi, kmax)

	            structobject    = getfield(HierarchicalMultivariateAnimalMovement,Symbol(string("VecPar",  prior_pi["name"])))
	            mcmc_pi[ian]        = structobject(p_aftercheck)
	        end

	        mcmc_initpi[ian]              = 1.0/(Float64(kmax))*ones(Float64,kmax)
	        mcmc_initz[ian]               = ones(Int16,1)
	        mcmc_beta[ian]                = 1.0/kmax*ones(Float64,kmax)
	        mcmc_rho[ian]                 = [0.5]
	        mcmc_gamma[ian]               = [1.0]
	        mcmc_ak[ian]                  = [1.0]



	        Clusterization[ian] = Clusterization_HDPHMM(Likelihood[ian].clusterization,mcmc_pi[ian],mcmc_initpi[ian],mcmc_initz[ian],mcmc_beta[ian],mcmc_rho[ian],mcmc_gamma[ian],  mcmc_ak[ian],prior_ak,prior_gamma,prior_rho)

	    elseif clustering_type == "HDP-HMM-Joint"
			# pi
	        if update_pi==false
	            mcmc_pi[ian] = PiDoNotUpdate()
	        else

	            checkPriorName  = size(findall(["Dirichlet"] .== prior_pi["name"]),1)== 0
	            if checkPriorName error("pi prior must be one of Dirichlet") end

	            namecheck       = string("check_Clusterization",prior_pi["name"])
	            nameprior       = string(prior_pi["name"],"_pi")
	            checkfunc       = getfield(HierarchicalMultivariateAnimalMovement,Symbol(namecheck))
	            p_aftercheck    = checkfunc(nameprior,inits_pi[ian], prior_pi, kmax)

	            structobject    = getfield(HierarchicalMultivariateAnimalMovement,Symbol(string("VecPar",  prior_pi["name"])))
	            mcmc_pi[ian]        = structobject(p_aftercheck)
	        end

	        mcmc_initpi[ian]              = 1.0/(Float64(kmax))*ones(Float64,kmax)
	        mcmc_initz[ian]               = ones(Int16,1)
	        mcmc_beta[ian]                = 1.0/kmax*ones(Float64,kmax)
	        mcmc_rho[ian]                 = [0.9]
	        mcmc_gamma[ian]               = [0.1]
	        mcmc_ak[ian]                  = [1.0]



	        Clusterization[ian] = Clusterization_HDPHMM_Joint(Likelihood[ian].clusterization,mcmc_pi[ian],mcmc_initpi[ian],mcmc_initz[ian],mcmc_beta[ian],mcmc_rho[ian],mcmc_gamma[ian],  mcmc_ak[ian],prior_ak,prior_gamma,prior_rho)
	    else
			error("clustering_type must be one of  ...")
	    end
	end



    return Clusterization

end

#### #### #### #### #### #### #### #### #### #### #### ####
#### #### #### #### #### MCMC interations
#### #### #### #### #### #### #### #### #### #### #### ####
