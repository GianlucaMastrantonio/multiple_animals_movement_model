

function sample_beta!(rng::MersenneTwister,Clusterization::Vector{Clusterization_HDPHMM}, piPar::VecParDirichlet, parHier::OptionHierarchicalParameters)

    #sample_beta!(Clusterization[1], Clusterization[1].pi)

    kmax  = Clusterization[1].clusterization.kmax
    nanim = size(Clusterization)[1]



    for ian in 1:nanim

        for k in 1:kmax
            app = Float64(0.0)
            j1 = parHier.h_mu.clust[k,ian]
            app += log(parHier.h_mu.prob[j1])
            j1 = parHier.h_eta.clust[k,ian]
            app += log(parHier.h_eta.prob[j1])
            j1 = parHier.h_nu.clust[k,ian]
            app += log(parHier.h_nu.prob[j1])
            j1 = parHier.h_rho.clust[k,ian]
            app += log(parHier.h_rho.prob[j1])
            j1 = parHier.h_sigma.clust[k,ian]
            app += log(parHier.h_sigma.prob[j1])

            Clusterization[ian].mcmc_beta[k] = exp(app)
        end
        Clusterization[ian].mcmc_beta[Clusterization[ian].mcmc_beta.<1.0e-100] .= 1.0e-100
    end

end



function sample_beta_singleHDP!(rng::MersenneTwister,Clusterization::Vector{Clusterization_HDPHMM}, piPar::VecParDirichlet, parHier::OptionHierarchicalParameters)

    #sample_beta!(Clusterization[1], Clusterization[1].pi)

    kmax  = Clusterization[1].clusterization.kmax
    nanim = size(Clusterization)[1]



    for ian in 1:nanim

        for k in 1:kmax

            Clusterization[ian].mcmc_beta[k] = parHier.h_mu.prob[k]
        end
        Clusterization[ian].mcmc_beta[Clusterization[ian].mcmc_beta.<1.0e-100] .= 1.0e-100
    end

end



function sample_beta_type2!(rng::MersenneTwister,Clusterization::Vector{Clusterization_HDPHMM}, piPar::VecParDirichlet, parHier::OptionHierarchicalParameters)

    #sample_beta!(Clusterization[1], Clusterization[1].pi)

    kmax  = Clusterization[1].clusterization.kmax
    nanim = size(Clusterization)[1]

    ian   = Int16(1)
    for k in 1:kmax
        app = Float64(0.0)
        j1 = parHier.h_mu.clust[k,ian]
        app += log(parHier.h_mu.prob[j1])
        j1 = parHier.h_eta.clust[k,ian]
        app += log(parHier.h_eta.prob[j1])
        j1 = parHier.h_nu.clust[k,ian]
        app += log(parHier.h_nu.prob[j1])
        j1 = parHier.h_rho.clust[k,ian]
        app += log(parHier.h_rho.prob[j1])
        j1 = parHier.h_sigma.clust[k,ian]
        app += log(parHier.h_sigma.prob[j1])

        Clusterization[ian].mcmc_beta[k] = exp(app)
    end
    Clusterization[ian].mcmc_beta[Clusterization[ian].mcmc_beta.<1.0e-100] .= 1.0e-100



    for ian in 1:nanim
        for k in 1:kmax
            Clusterization[ian].mcmc_beta[k] = Clusterization[1].mcmc_beta[k]
        end
    end

end


function sample_beta!(rng::MersenneTwister,Clusterization::Vector{Clusterization_HDPHMM_Joint}, piPar::VecParDirichlet)

    sample_beta!(Clusterization[1], Clusterization[1].pi)

    nanim = size(Clusterization)[1]
    if nanim>1
        for ian in 2:nanim
            for i1 in 1:size(Clusterization[1].mcmc_beta)[1]
                Clusterization[ian].mcmc_beta[i1] = Clusterization[1].mcmc_beta[i1]
            end
        end
    end

end

function sample_beta!(rng::MersenneTwister,Clusterization::Clusterization_HDPHMM_Joint, piPar::VecParDirichlet)

    #
    #Clusterization = MCMCClusterization[2]
    #piPar = MCMCClusterization.pi
    MatM        = Clusterization.MatM
    MatMbar     = Clusterization.MatMbar
    ktot        = Clusterization.clusterization.kmax
    nt          = size(Clusterization.clusterization.zeta,1)
    # MatM        = zeros(Float64,ktot,ktot)
    # MatMbar     = zeros(Float64,ktot,ktot)
    nnonempty   = Clusterization.clusterization.n_nonemptyC[1]

    gammaPar    = Clusterization.mcmc_gamma[1]
    rhoPar      = Clusterization.mcmc_rho[1]
    kappaPar    = Clusterization.mcmc_rho[1]*Clusterization.mcmc_ak[1]
    alphaPar    = Clusterization.mcmc_ak[1]-kappaPar
    betaPar     = Clusterization.mcmc_beta

    logrj = 0.0
    sj    = 0.0




    DirPar  = gammaPar*Float64.(ones(ktot)/ktot)

    for k1 in 1:nnonempty
        kk1 = Clusterization.clusterization.nonemptyC[k1]
        DirPar[kk1] += sum(MatMbar[:,kk1])
    end

    betaPar[1:ktot]  = rand(rng,Dirichlet(DirPar))
    betaPar[betaPar.<1.0e-100] .= 1.0e-100


    return nothing
end


function sample_beta!(rng::MersenneTwister,Clusterization::Vector{Clusterization_HDPHMM}, piPar::VecParDirichlet)

    nanim = size(Clusterization)[1]
    for ian in 1:nanim
        sample_beta!(rng,Clusterization[ian], Clusterization[ian].pi)
    end

end

function sample_beta!(rng::MersenneTwister,Clusterization::Clusterization_HDPHMM, piPar::VecParDirichlet)

    #
    #Clusterization = MCMCClusterization
    # piPar = MCMCClusterization.pi
    MatM        = Clusterization.MatM
    MatMbar     = Clusterization.MatMbar
    ktot        = Clusterization.clusterization.kmax
    nt          = size(Clusterization.clusterization.zeta,1)
    # MatM        = zeros(Float64,ktot,ktot)
    # MatMbar     = zeros(Float64,ktot,ktot)
    nnonempty   = Clusterization.clusterization.n_nonemptyC[1]

    gammaPar    = Clusterization.mcmc_gamma[1]
    rhoPar      = Clusterization.mcmc_rho[1]
    kappaPar    = Clusterization.mcmc_rho[1]*Clusterization.mcmc_ak[1]
    alphaPar    = Clusterization.mcmc_ak[1]-kappaPar
    betaPar     = Clusterization.mcmc_beta

    logrj = 0.0
    sj    = 0.0




    DirPar  = gammaPar*Float64.(ones(ktot)/ktot)

    for k1 in 1:nnonempty
        kk1 = Clusterization.clusterization.nonemptyC[k1]
        DirPar[kk1] += sum(MatMbar[:,kk1])
    end

    betaPar[1:ktot]  = rand(rng,Dirichlet(DirPar))
    betaPar[betaPar.<1.0e-100] .= 1.0e-100


    return nothing
end
