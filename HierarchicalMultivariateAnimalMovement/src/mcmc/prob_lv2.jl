function sample_prob_lv2!(rng::MersenneTwister,Clusterization::Vector{Clusterization_HDPHMM}, piPar::VecParDirichlet, parHier::OptionHierarchicalParameters)

    #sample_beta!(Clusterization[1], Clusterization[1].pi)

    kmax  = Clusterization[1].clusterization.kmax
    nanim = size(Clusterization)[1]
    ktot  = kmax

    # #Clusterization = MCMCClusterization
    # # piPar = MCMCClusterization.pi

    # ktot        = Clusterization.clusterization.kmax
    # nt          = size(Clusterization.clusterization.zeta,1)
    # # MatM        = zeros(Float64,ktot,ktot)
    # # MatMbar     = zeros(Float64,ktot,ktot)
    # nnonempty   = Clusterization.clusterization.n_nonemptyC[1]
    #
    # gammaPar    = Clusterization.mcmc_gamma[1]
    # rhoPar      = Clusterization.mcmc_rho[1]
    # kappaPar    = Clusterization.mcmc_rho[1]*Clusterization.mcmc_ak[1]
    # alphaPar    = Clusterization.mcmc_ak[1]-kappaPar
    # betaPar     = Clusterization.mcmc_beta
    #
    # logrj = 0.0
    # sj    = 0.0
    #
    #      DirPar[kk1] += sum(MatMbar[:,kk1])
    #     end
    #
    #     betaPar[1:ktot]  = rand(Dirichlet(DirPar))
    #
    gammaPar        = Clusterization[1].mcmc_gamma[1]
    DirPar_mu       = gammaPar*Float64.(ones(ktot)/ktot)
    DirPar_nu       = gammaPar*Float64.(ones(ktot)/ktot)
    DirPar_eta      = gammaPar*Float64.(ones(ktot)/ktot)
    DirPar_sigma    = gammaPar*Float64.(ones(ktot)/ktot)
    DirPar_rho      = gammaPar*Float64.(ones(ktot)/ktot)

    for ian in 1:nanim
        MatMsum         = sum(Clusterization[ian].MatM,dims=1)

        #for k in 1:kmax
        for k1 in 1:Clusterization[ian].clusterization.n_nonemptyC[1]
            k = Clusterization[ian].clusterization.nonemptyC[k1]

            j = parHier.h_mu.clust[k,ian]
            DirPar_mu[j] += MatMsum[k]

            j = parHier.h_nu.clust[k,ian]
            DirPar_nu[j] += MatMsum[k]

            j = parHier.h_eta.clust[k,ian]
            DirPar_eta[j] += MatMsum[k]

            j = parHier.h_sigma.clust[k,ian]
            DirPar_sigma[j] += MatMsum[k]

            j = parHier.h_rho.clust[k,ian]
            DirPar_rho[j] += MatMsum[k]

        end

    end
    parHier.h_mu.prob[1:kmax]  = rand(rng,Dirichlet(DirPar_mu))
    parHier.h_mu.prob[parHier.h_mu.prob.<1.0e-100] .= 1.0e-100

    parHier.h_nu.prob[1:kmax]  = rand(rng,Dirichlet(DirPar_nu))
    parHier.h_nu.prob[parHier.h_nu.prob.<1.0e-100] .= 1.0e-100

    parHier.h_eta.prob[1:kmax]  = rand(rng,Dirichlet(DirPar_eta))
    parHier.h_eta.prob[parHier.h_eta.prob.<1.0e-100] .= 1.0e-100

    parHier.h_sigma.prob[1:kmax]  = rand(rng,Dirichlet(DirPar_sigma))
    parHier.h_sigma.prob[parHier.h_sigma.prob.<1.0e-100] .= 1.0e-100

    parHier.h_rho.prob[1:kmax]  = rand(rng,Dirichlet(DirPar_rho))
    parHier.h_rho.prob[parHier.h_rho.prob.<1.0e-100] .= 1.0e-100

    return nothing
end




function sample_prob_lv2_type2!(rng::MersenneTwister,Clusterization::Vector{Clusterization_HDPHMM}, piPar::VecParDirichlet, parHier::OptionHierarchicalParameters)

    #sample_beta!(Clusterization[1], Clusterization[1].pi)

    kmax  = Clusterization[1].clusterization.kmax
    nanim = size(Clusterization)[1]
    ktot  = kmax

    # #Clusterization = MCMCClusterization
    # # piPar = MCMCClusterization.pi

    # ktot        = Clusterization.clusterization.kmax
    # nt          = size(Clusterization.clusterization.zeta,1)
    # # MatM        = zeros(Float64,ktot,ktot)
    # # MatMbar     = zeros(Float64,ktot,ktot)
    # nnonempty   = Clusterization.clusterization.n_nonemptyC[1]
    #
    # gammaPar    = Clusterization.mcmc_gamma[1]
    # rhoPar      = Clusterization.mcmc_rho[1]
    # kappaPar    = Clusterization.mcmc_rho[1]*Clusterization.mcmc_ak[1]
    # alphaPar    = Clusterization.mcmc_ak[1]-kappaPar
    # betaPar     = Clusterization.mcmc_beta
    #
    # logrj = 0.0
    # sj    = 0.0
    #
    #      DirPar[kk1] += sum(MatMbar[:,kk1])
    #     end
    #
    #     betaPar[1:ktot]  = rand(Dirichlet(DirPar))
    #
    gammaPar        = Clusterization[1].mcmc_gamma[1]
    DirPar_mu       = gammaPar*Float64.(ones(ktot)/ktot)
    DirPar_nu       = gammaPar*Float64.(ones(ktot)/ktot)
    DirPar_eta      = gammaPar*Float64.(ones(ktot)/ktot)
    DirPar_sigma    = gammaPar*Float64.(ones(ktot)/ktot)
    DirPar_rho      = gammaPar*Float64.(ones(ktot)/ktot)

    for ian in 1:1
        MatMsum         = sum(Clusterization[ian].MatM,dims=1)

        #for k in 1:kmax
        for k1 in 1:Clusterization[ian].clusterization.n_nonemptyC[1]
            k = Clusterization[ian].clusterization.nonemptyC[k1]

            j = parHier.h_mu.clust[k,ian]
            DirPar_mu[j] += MatMsum[k]

            j = parHier.h_nu.clust[k,ian]
            DirPar_nu[j] += MatMsum[k]

            j = parHier.h_eta.clust[k,ian]
            DirPar_eta[j] += MatMsum[k]

            j = parHier.h_sigma.clust[k,ian]
            DirPar_sigma[j] += MatMsum[k]

            j = parHier.h_rho.clust[k,ian]
            DirPar_rho[j] += MatMsum[k]

        end

    end
    parHier.h_mu.prob[1:kmax]  = rand(rng,Dirichlet(DirPar_mu))
    parHier.h_mu.prob[parHier.h_mu.prob.<1.0e-100] .= 1.0e-100

    parHier.h_nu.prob[1:kmax]  = rand(rng,Dirichlet(DirPar_nu))
    parHier.h_nu.prob[parHier.h_nu.prob.<1.0e-100] .= 1.0e-100

    parHier.h_eta.prob[1:kmax]  = rand(rng,Dirichlet(DirPar_eta))
    parHier.h_eta.prob[parHier.h_eta.prob.<1.0e-100] .= 1.0e-100

    parHier.h_sigma.prob[1:kmax]  = rand(rng,Dirichlet(DirPar_sigma))
    parHier.h_sigma.prob[parHier.h_sigma.prob.<1.0e-100] .= 1.0e-100

    parHier.h_rho.prob[1:kmax]  = rand(rng,Dirichlet(DirPar_rho))
    parHier.h_rho.prob[parHier.h_rho.prob.<1.0e-100] .= 1.0e-100

    return nothing
end
