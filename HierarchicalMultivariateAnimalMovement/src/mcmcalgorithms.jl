
#### GENERAL
function MCMCalgorithm(
    MCMCout::MCMCutils,
    MCMCLikelihood::AbstractLikelihood,
    MCMCClusterization::AbstractClusterization
)
    error("MCMC non implemented for Likelihood ",typeof(Likelihood), " and Clusterization ",typeof(Clusterization))
end

#### MODEL1

function MCMCalgorithm(
    MCMCout::Vector{MCMCutils},
    MCMCLikelihood::Vector{Likelihood_OU_CircLinmodel},
    #MCMCClusterization::Vector{AbstractClusterization}
    MCMCClusterization::Vector{T1},
	rng::MersenneTwister
)::Dict where {T1<:AbstractClusterization}

nanim = size(MCMCout)[1]

for ian in 1:nanim
MCMCout[ian].indexsave[1]    = 1
end


appBurninOrThin         = MCMCout[1].burnin
iterations = Int64(0)
for saveMCMCindex in 1:MCMCout[1].nsamplesave

  	for burnithinMCMCindex in 1:appBurninOrThin
		#print(rand(),"\n")
		if mod(iterations,50) == 0 print(string("Iterations=",iterations,"\n"))  end
		if mod(iterations,50) == 0
		for ian in 1:nanim
		    print(string("Clusters=",MCMCClusterization[ian].clusterization.n_nonemptyC[1],"\n"))
		end
		end
		iterations += one(iterations)

		### ### ### ### ### ### ###
		### ### ZETA
		### ### ### ### ### ### ###

		# zeta
		sample_zeta!(rng,MCMCLikelihood, MCMCClusterization,MCMCLikelihood[1].clusterization)
		sample_zeta_mergesplit!(rng,MCMCLikelihood, MCMCClusterization,MCMCLikelihood[1].clusterization)
		#print(rand(rng), " zeta \n")
		sample_pi!(rng,MCMCClusterization, MCMCClusterization[1].pi)
		#print(rand(rng), " pi \n")


		### ### ### ### ### ### ###
		### ### CLUSTERIZATION
		### ### ### ### ### ### ###

		#print("iterations ",iterations," ",rand(rng), "\n")
		### pi

		compute_m!(rng,MCMCClusterization, MCMCClusterization[1].pi)
		#print(rand(rng), " m \n")


		sample_beta!(rng, MCMCClusterization, MCMCClusterization[1].pi)
		#print(rand(rng), " beta \n")



		sample_ak!(rng,MCMCClusterization, MCMCClusterization[1].pi)
		#print(rand(rng), " ak \n")
		sample_gamma!(rng,MCMCClusterization, MCMCClusterization[1].pi)
		#print(rand(rng), " gamm< \n")
		sample_rhodp!(rng,MCMCClusterization, MCMCClusterization[1].pi)
		#print(rand(rng), " rhodp \n")




		### ### ### ### ### ### ###
		### ### LIKELIHOOD
		### ### ### ### ### ### ###

		### missing
		sample_missing!(rng,MCMCLikelihood,MCMCLikelihood[1].miss)
		#print(rand(rng), " miss \n")

		### mu0 (mu)
		sample_mu0!(rng,MCMCLikelihood,MCMCLikelihood[1].mu)
		#print(rand(rng), " mu \n")

		### sigma
		sample_sigma!(rng,MCMCLikelihood, MCMCLikelihood[1].sigma)
		#print(rand(rng), " sigma \n")

		### psi
		sample_psi!(rng,MCMCLikelihood, MCMCLikelihood[1].nu)
		#print(rand(rng), " psi \n")

		# eta
		sample_muC!(rng,MCMCLikelihood,MCMCLikelihood[1].eta)
		#print(rand(rng), " muC \n")

		# rho
		sample_rho!(rng,MCMCLikelihood,MCMCLikelihood[1].rho)
		#print(rand(rng), " rho \n")



	end # second for MCMC
	appBurninOrThin = MCMCout[1].thin

	save_posteriorsamples!(MCMCout,MCMCLikelihood, MCMCClusterization)

end # first for MCMC

# for saveMCMCindex in 1:MCMCout[1].nsamplesave
#
#   	for burnithinMCMCindex in 1:appBurninOrThin
# 		#print(rand(),"\n")
# 		if mod(iterations,50) == 0 print(string("Iterations=",iterations,"\n"))  end
# 		if mod(iterations,50) == 0
# 		for ian in 1:nanim
# 		    print(string("Clusters=",MCMCClusterization[ian].clusterization.n_nonemptyC[1],"\n"))
# 		end
# 		end
# 		iterations += one(iterations)
#
# 		### ### ### ### ### ### ###
# 		### ### CLUSTERIZATION
# 		### ### ### ### ### ### ###
#
# 		print("iterations ",iterations," ",rand(rng), "\n")
# 		### pi
#
# 		compute_m!(rng,MCMCClusterization, MCMCClusterization[1].pi)
# 		print(rand(rng), " m \n")
#
#
# 		sample_beta!(rng, MCMCClusterization, MCMCClusterization[1].pi)
# 		print(rand(rng), " beta \n")
#
#
# 		sample_pi!(rng,MCMCClusterization, MCMCClusterization[1].pi)
# 		print(rand(rng), " pi \n")
# 		sample_ak!(rng,MCMCClusterization, MCMCClusterization[1].pi)
# 		print(rand(rng), " ak \n")
# 		sample_gamma!(rng,MCMCClusterization, MCMCClusterization[1].pi)
# 		print(rand(rng), " gamm< \n")
# 		sample_rhodp!(rng,MCMCClusterization, MCMCClusterization[1].pi)
# 		print(rand(rng), " rhodp \n")
#
# 		### ### ### ### ### ### ###
# 		### ### ZETA
# 		### ### ### ### ### ### ###
#
# 		# zeta
# 		sample_zeta!(rng,MCMCLikelihood, MCMCClusterization,MCMCLikelihood[1].clusterization)
# 		print(rand(rng), " zeta \n")
#
#
#
#
# 		### ### ### ### ### ### ###
# 		### ### LIKELIHOOD
# 		### ### ### ### ### ### ###
#
# 		### missing
# 		sample_missing!(rng,MCMCLikelihood,MCMCLikelihood[1].miss)
# 		print(rand(rng), " miss \n")
#
# 		### mu0 (mu)
# 		sample_mu0!(rng,MCMCLikelihood,MCMCLikelihood[1].mu)
# 		print(rand(rng), " mu \n")
#
# 		### sigma
# 		sample_sigma!(rng,MCMCLikelihood, MCMCLikelihood[1].sigma)
# 		print(rand(rng), " sigma \n")
#
# 		### psi
# 		sample_psi!(rng,MCMCLikelihood, MCMCLikelihood[1].nu)
# 		print(rand(rng), " psi \n")
#
# 		# eta
# 		sample_muC!(rng,MCMCLikelihood,MCMCLikelihood[1].eta)
# 		print(rand(rng), " muC \n")
#
# 		# rho
# 		sample_rho!(rng,MCMCLikelihood,MCMCLikelihood[1].rho)
# 		print(rand(rng), " rho \n")
#
# 		if iterations >40
# 			error(" ")
# 		end
#
# 	end # second for MCMC
# 	appBurninOrThin = MCMCout[1].thin
#
# 	save_posteriorsamples!(MCMCout,MCMCLikelihood, MCMCClusterization)
#
# end


return Dict("PosteriorSamples"=>MCMCout, "Likelihood"=>MCMCLikelihood, "Clusterization"=>MCMCClusterization)
end






#### MODEL1

function MCMCalgorithm(
    MCMCout::Vector{MCMCutils},
    MCMCLikelihood::Vector{Likelihood_OU_CircLinmodel},
    MCMCClusterization::Vector{T1},
    MCMCHierarchicalParameters::OptionHierarchicalParameters,
	rng::MersenneTwister;
	Iter_lev2 = 10::Integer,
	Iter_print = 50::Integer
)::Dict where {T1<:AbstractClusterization}

    nanim = size(MCMCout)[1]

    for ian in 1:nanim
        MCMCout[ian].indexsave[1]    = 1
    end


    appBurninOrThin         = MCMCout[1].burnin
    iterations = Int64(0)
    for saveMCMCindex in 1:MCMCout[1].nsamplesave

        for burnithinMCMCindex in 1:appBurninOrThin

            if mod(iterations,Iter_print) == 0 print(string("Iterations_DHV2=",iterations,"\n"))  end
            if mod(iterations,Iter_print) == 0
                for ian in 1:nanim
                    print(string("Clusters=",MCMCClusterization[ian].clusterization.n_nonemptyC[1],"\n"))
                end
             end
            iterations += one(iterations)





	        ### ### ### ### ### ### ###
	        ### ### LIKELIHOOD
	        ### ### ### ### ### ### ###

			### missing
	        sample_missing!(rng,MCMCLikelihood,MCMCLikelihood[1].miss)

			if iterations>Iter_lev2

				sample_muC_lv2!(rng,MCMCLikelihood,MCMCLikelihood[1].eta,MCMCHierarchicalParameters, MCMCClusterization)
				#sample_muC_split_merge!(rng,MCMCLikelihood,MCMCLikelihood[1].eta,MCMCHierarchicalParameters, MCMCClusterization)


				sample_mu0_lv2!(rng,MCMCLikelihood,MCMCLikelihood[1].mu,MCMCHierarchicalParameters, MCMCClusterization)
				#sample_mu0_split_merge!(rng,MCMCLikelihood,MCMCLikelihood[1].mu,MCMCHierarchicalParameters, MCMCClusterization)

				sample_rho_lv2!(rng,MCMCLikelihood,MCMCLikelihood[1].rho,MCMCHierarchicalParameters, MCMCClusterization)
				#sample_rho_lv2!(rng,MCMCLikelihood,MCMCLikelihood[1].rho,MCMCHierarchicalParameters, MCMCClusterization)

				sample_psi_lv2!(rng,MCMCLikelihood,MCMCLikelihood[1].nu,MCMCHierarchicalParameters, MCMCClusterization)
				#sample_psi_split_merge!(rng,MCMCLikelihood,MCMCLikelihood[1].nu,MCMCHierarchicalParameters, MCMCClusterization)

				sample_sigma_lv2!(rng,MCMCLikelihood,MCMCLikelihood[1].sigma,MCMCHierarchicalParameters, MCMCClusterization)
				#sample_sigma_split_merge!(rng,MCMCLikelihood,MCMCLikelihood[1].sigma,MCMCHierarchicalParameters, MCMCClusterization)

			end

			#eta
			sample_muC!(rng,MCMCLikelihood,MCMCLikelihood[1].eta,MCMCHierarchicalParameters)

	        #mu0 (mu)
			sample_mu0!(rng,MCMCLikelihood,MCMCLikelihood[1].mu,MCMCHierarchicalParameters)

	        # sigma
	        sample_sigma!(rng,MCMCLikelihood, MCMCLikelihood[1].sigma,MCMCHierarchicalParameters)

	        # ### psi
	        sample_psi!(rng,MCMCLikelihood, MCMCLikelihood[1].nu,MCMCHierarchicalParameters)
			# rho
			sample_rho!(rng,MCMCLikelihood,MCMCLikelihood[1].rho,MCMCHierarchicalParameters)


			### ### ### ### ### ### ###
			### ### ZETA
			### ### ### ### ### ### ###

			# zeta
			#sample_zeta!(rng,MCMCLikelihood, MCMCClusterization,MCMCLikelihood[1].clusterization)

			if iterations>50
				if rand(rng,Uniform(0.0,1.0))<0.9
					sample_zeta!(rng,MCMCLikelihood, MCMCClusterization,MCMCLikelihood[1].clusterization)
				else
					sample_zeta_mergesplit!(rng,MCMCLikelihood, MCMCClusterization,MCMCLikelihood[1].clusterization)
				end
			else
				sample_zeta_app!(rng,MCMCLikelihood, MCMCClusterization,MCMCLikelihood[1].clusterization)
			end



			sample_pi!(rng,MCMCClusterization, MCMCClusterization[1].pi)



            ### ### ### ### ### ### ###
			### ### CLUSTERIZATION - Second Level
			### ### ### ### ### ### ###

			compute_m!(rng,MCMCClusterization, MCMCClusterization[1].pi)
			sample_prob_lv2!(rng,MCMCClusterization, MCMCClusterization[1].pi,MCMCHierarchicalParameters)
			sample_beta!(rng,MCMCClusterization, MCMCClusterization[1].pi,MCMCHierarchicalParameters)
			### ### ### ### ### ### ###
			### ### CLUSTERIZATION
			### ### ### ### ### ### ###

			### pi




			sample_ak!(rng,MCMCClusterization, MCMCClusterization[1].pi,MCMCHierarchicalParameters)
			sample_rhodp!(rng,MCMCClusterization, MCMCClusterization[1].pi,MCMCHierarchicalParameters)
			#sample_gamma!(rng,MCMCClusterization, MCMCClusterization[1].pi,MCMCHierarchicalParameters)





        end # second for MCMC
        appBurninOrThin = MCMCout[1].thin

        save_posteriorsamples!(MCMCout,MCMCLikelihood, MCMCClusterization,MCMCHierarchicalParameters)


    end # first for MCMC
    return Dict("PosteriorSamples"=>MCMCout, "Likelihood"=>MCMCLikelihood, "Clusterization"=>MCMCClusterization)
end




#### MODEL1

function MCMCalgorithm_joint(
    MCMCout::Vector{MCMCutils},
    MCMCLikelihood::Vector{Likelihood_OU_CircLinmodel},
    MCMCClusterization::Vector{T1},
    MCMCHierarchicalParameters::OptionHierarchicalParameters,
	rng::MersenneTwister;
	Iter_lev2 = 0::Integer,
	Iter_print = 50::Integer
)::Dict where {T1<:AbstractClusterization}

    nanim = size(MCMCout)[1]

    for ian in 1:nanim
        MCMCout[ian].indexsave[1]    = 1
    end
	print("\n \n \n \n \n \n")
	print("no Split e Merge")
	print("\n \n \n \n \n \n")

    appBurninOrThin         = MCMCout[1].burnin
    iterations = Int64(0)
    for saveMCMCindex in 1:MCMCout[1].nsamplesave

        for burnithinMCMCindex in 1:appBurninOrThin

            if mod(iterations,Iter_print) == 0 print(string("Iterations_H2=",iterations,"\n"))  end
            if mod(iterations,Iter_print) == 0
                for ian in 1:nanim
                    print(string("Clusters=",MCMCClusterization[ian].clusterization.n_nonemptyC[1],"\n"))
                end
             end
            iterations += one(iterations)




	        ### ### ### ### ### ### ###
	        ### ### LIKELIHOOD
	        ### ### ### ### ### ### ###

			### missing
	        sample_missing!(rng,MCMCLikelihood,MCMCLikelihood[1].miss)

			if iterations>Iter_lev2

				sample_muC_lv2_type2!(rng,MCMCLikelihood,MCMCLikelihood[1].eta,MCMCHierarchicalParameters, MCMCClusterization)
				#sample_muC_split_merge!(rng,MCMCLikelihood,MCMCLikelihood[1].eta,MCMCHierarchicalParameters, MCMCClusterization)


				sample_mu0_lv2_type2!(rng,MCMCLikelihood,MCMCLikelihood[1].mu,MCMCHierarchicalParameters, MCMCClusterization)
				#sample_mu0_split_merge!(rng,MCMCLikelihood,MCMCLikelihood[1].mu,MCMCHierarchicalParameters, MCMCClusterization)

				sample_rho_lv2_type2!(rng,MCMCLikelihood,MCMCLikelihood[1].rho,MCMCHierarchicalParameters, MCMCClusterization)
				#sample_rho_lv2!(rng,MCMCLikelihood,MCMCLikelihood[1].rho,MCMCHierarchicalParameters, MCMCClusterization)

				sample_psi_lv2_type2!(rng,MCMCLikelihood,MCMCLikelihood[1].nu,MCMCHierarchicalParameters, MCMCClusterization)
				#sample_psi_split_merge!(rng,MCMCLikelihood,MCMCLikelihood[1].nu,MCMCHierarchicalParameters, MCMCClusterization)

				sample_sigma_lv2_type2!(rng,MCMCLikelihood,MCMCLikelihood[1].sigma,MCMCHierarchicalParameters, MCMCClusterization)
				#sample_sigma_split_merge!(rng,MCMCLikelihood,MCMCLikelihood[1].sigma,MCMCHierarchicalParameters, MCMCClusterization)

			end

			#eta
			sample_muC!(rng,MCMCLikelihood,MCMCLikelihood[1].eta,MCMCHierarchicalParameters)

	        #mu0 (mu)
			sample_mu0!(rng,MCMCLikelihood,MCMCLikelihood[1].mu,MCMCHierarchicalParameters)

	        # sigma
	        sample_sigma!(rng,MCMCLikelihood, MCMCLikelihood[1].sigma,MCMCHierarchicalParameters)

	        # ### psi
	        sample_psi!(rng,MCMCLikelihood, MCMCLikelihood[1].nu,MCMCHierarchicalParameters)
			# rho
			sample_rho!(rng,MCMCLikelihood,MCMCLikelihood[1].rho,MCMCHierarchicalParameters)


			### ### ### ### ### ### ###
			### ### ZETA
			### ### ### ### ### ### ###

			# zeta
			sample_zeta_type2!(rng,MCMCLikelihood, MCMCClusterization,MCMCLikelihood[1].clusterization)

			if rand(rng,Uniform(0.0,1.0))<1.1
				sample_zeta_mergesplit!(rng,MCMCLikelihood, MCMCClusterization,MCMCLikelihood[1].clusterization)
			end
			sample_pi_type2!(rng,MCMCClusterization, MCMCClusterization[1].pi)
            ### ### ### ### ### ### ###
			### ### CLUSTERIZATION - Second Level
			### ### ### ### ### ### ###

			sample_prob_lv2_type2!(rng,MCMCClusterization, MCMCClusterization[1].pi,MCMCHierarchicalParameters)
			sample_beta_type2!(rng,MCMCClusterization, MCMCClusterization[1].pi,MCMCHierarchicalParameters)
			### ### ### ### ### ### ###
			### ### CLUSTERIZATION
			### ### ### ### ### ### ###

			### pi
			compute_m_type2!(rng,MCMCClusterization, MCMCClusterization[1].pi)



			sample_ak_type2!(rng,MCMCClusterization, MCMCClusterization[1].pi,MCMCHierarchicalParameters)
			sample_rhodp_type2!(rng,MCMCClusterization, MCMCClusterization[1].pi,MCMCHierarchicalParameters)
			sample_gamma_type2!(rng,MCMCClusterization, MCMCClusterization[1].pi,MCMCHierarchicalParameters)






        end # second for MCMC
        appBurninOrThin = MCMCout[1].thin

        save_posteriorsamples!(MCMCout,MCMCLikelihood, MCMCClusterization,MCMCHierarchicalParameters)


    end # first for MCMC
    return Dict("PosteriorSamples"=>MCMCout, "Likelihood"=>MCMCLikelihood, "Clusterization"=>MCMCClusterization)
end











#### #### #### #### #### #### ####
#### MODEL SINGLE HDP
#### #### #### #### #### #### ####

function MCMCalgorithm_singleHDP(
    MCMCout::Vector{MCMCutils},
    MCMCLikelihood::Vector{Likelihood_OU_CircLinmodel},
    MCMCClusterization::Vector{T1},
    MCMCHierarchicalParameters::OptionHierarchicalParameters,
	rng::MersenneTwister;
	Iter_lev2 = 0::Integer,
	Iter_print = 50::Integer
)::Dict where {T1<:AbstractClusterization}

    nanim = size(MCMCout)[1]

    for ian in 1:nanim
        MCMCout[ian].indexsave[1]    = 1
    end




	nanim_v2    = size(MCMCLikelihood)[1]
    nanim       = nanim_v2
    kmax        = MCMCLikelihood[1].kmax

	for ian in 1:nanim_v2
		for k in 1:kmax
			MCMCHierarchicalParameters.h_mu.clust[k,ian] = k
			MCMCHierarchicalParameters.h_eta.clust[k,ian] = k
			MCMCHierarchicalParameters.h_nu.clust[k,ian] = k
			MCMCHierarchicalParameters.h_rho.clust[k,ian] = k
			MCMCHierarchicalParameters.h_sigma.clust[k,ian] = k
		end
	end

	#
	#
	# app = Float64(0.0)
	# j1 = parHier.h_mu.clust[k,ian]
	# app += log(parHier.h_mu.prob[j1])
	# j1 = parHier.h_eta.clust[k,ian]
	# app += log(parHier.h_eta.prob[j1])
	# j1 = parHier.h_nu.clust[k,ian]
	# app += log(parHier.h_nu.prob[j1])
	# j1 = parHier.h_rho.clust[k,ian]
	# app += log(parHier.h_rho.prob[j1])
	# j1 = parHier.h_sigma.clust[k,ian]
	# app += log(parHier.h_sigma.prob[j1])

    appBurninOrThin         = MCMCout[1].burnin
    iterations = Int64(0)
    for saveMCMCindex in 1:MCMCout[1].nsamplesave

        for burnithinMCMCindex in 1:appBurninOrThin

            if mod(iterations,Iter_print) == 0 print(string("Iterations_H=",iterations,"\n"))  end
            if mod(iterations,Iter_print) == 0
                for ian in 1:nanim
                    print(string("Clusters=",MCMCClusterization[ian].clusterization.n_nonemptyC[1],"\n"))
                end
             end
            iterations += one(iterations)





	        ### ### ### ### ### ### ###
	        ### ### LIKELIHOOD
	        ### ### ### ### ### ### ###

			### missing
	        sample_missing!(rng,MCMCLikelihood,MCMCLikelihood[1].miss)

			# if iterations>Iter_lev2
			#
			# 	sample_muC_lv2!(rng,MCMCLikelihood,MCMCLikelihood[1].eta,MCMCHierarchicalParameters, MCMCClusterization)
			# 	#sample_muC_split_merge!(rng,MCMCLikelihood,MCMCLikelihood[1].eta,MCMCHierarchicalParameters, MCMCClusterization)
			#
			#
			# 	sample_mu0_lv2!(rng,MCMCLikelihood,MCMCLikelihood[1].mu,MCMCHierarchicalParameters, MCMCClusterization)
			# 	#sample_mu0_split_merge!(rng,MCMCLikelihood,MCMCLikelihood[1].mu,MCMCHierarchicalParameters, MCMCClusterization)
			#
			# 	sample_rho_lv2!(rng,MCMCLikelihood,MCMCLikelihood[1].rho,MCMCHierarchicalParameters, MCMCClusterization)
			# 	#sample_rho_lv2!(rng,MCMCLikelihood,MCMCLikelihood[1].rho,MCMCHierarchicalParameters, MCMCClusterization)
			#
			# 	sample_psi_lv2!(rng,MCMCLikelihood,MCMCLikelihood[1].nu,MCMCHierarchicalParameters, MCMCClusterization)
			# 	#sample_psi_split_merge!(rng,MCMCLikelihood,MCMCLikelihood[1].nu,MCMCHierarchicalParameters, MCMCClusterization)
			#
			# 	sample_sigma_lv2!(rng,MCMCLikelihood,MCMCLikelihood[1].sigma,MCMCHierarchicalParameters, MCMCClusterization)
			# 	#sample_sigma_split_merge!(rng,MCMCLikelihood,MCMCLikelihood[1].sigma,MCMCHierarchicalParameters, MCMCClusterization)
			#
			# end

			#eta
			sample_muC!(rng,MCMCLikelihood,MCMCLikelihood[1].eta,MCMCHierarchicalParameters)

	        #mu0 (mu)
			sample_mu0!(rng,MCMCLikelihood,MCMCLikelihood[1].mu,MCMCHierarchicalParameters)

	        # sigma
	        sample_sigma!(rng,MCMCLikelihood, MCMCLikelihood[1].sigma,MCMCHierarchicalParameters)

	        # ### psi
	        sample_psi!(rng,MCMCLikelihood, MCMCLikelihood[1].nu,MCMCHierarchicalParameters)
			# rho
			sample_rho!(rng,MCMCLikelihood,MCMCLikelihood[1].rho,MCMCHierarchicalParameters)


			### ### ### ### ### ### ###
			### ### ZETA
			### ### ### ### ### ### ###

			# zeta
			sample_zeta!(rng,MCMCLikelihood, MCMCClusterization,MCMCLikelihood[1].clusterization)
			if rand(rng,Uniform(0.0,1.0))<1.1
				#sample_zeta_mergesplit!(rng,MCMCLikelihood, MCMCClusterization,MCMCLikelihood[1].clusterization)
				#sample_zeta_mergesplit_v2!(rng,MCMCLikelihood, MCMCClusterization,MCMCLikelihood[1].clusterization, MCMCHierarchicalParameters)
			end
			sample_pi!(rng,MCMCClusterization, MCMCClusterization[1].pi)

            ### ### ### ### ### ### ###
			### ### CLUSTERIZATION - Second Level
			### ### ### ### ### ### ###

			#compute_m!(rng,MCMCClusterization, MCMCClusterization[1].pi)
			#sample_prob_lv2!(rng,MCMCClusterization, MCMCClusterization[1].pi,MCMCHierarchicalParameters)
			#sample_beta!(rng,MCMCClusterization, MCMCClusterization[1].pi,MCMCHierarchicalParameters)

			compute_m_singleHDP!(rng,MCMCClusterization, MCMCClusterization[1].pi)
			sample_prob_lv2_singleHDP!(rng,MCMCClusterization, MCMCClusterization[1].pi,MCMCHierarchicalParameters)
			sample_beta_singleHDP!(rng,MCMCClusterization, MCMCClusterization[1].pi,MCMCHierarchicalParameters)
			### ### ### ### ### ### ###
			### ### CLUSTERIZATION
			### ### ### ### ### ### ###

			### pi




			sample_ak!(rng,MCMCClusterization, MCMCClusterization[1].pi,MCMCHierarchicalParameters)
			sample_rhodp!(rng,MCMCClusterization, MCMCClusterization[1].pi,MCMCHierarchicalParameters)
			#sample_gamma!(rng,MCMCClusterization, MCMCClusterization[1].pi,MCMCHierarchicalParameters)






        end # second for MCMC
        appBurninOrThin = MCMCout[1].thin

        save_posteriorsamples!(MCMCout,MCMCLikelihood, MCMCClusterization,MCMCHierarchicalParameters)


    end # first for MCMC
    return Dict("PosteriorSamples"=>MCMCout, "Likelihood"=>MCMCLikelihood, "Clusterization"=>MCMCClusterization)
end
