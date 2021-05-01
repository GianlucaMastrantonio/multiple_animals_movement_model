



function sample_rhodp!(rng::MersenneTwister,Clusterization::Vector{Clusterization_HDPHMM}, piPar::VecParDirichlet,parHier::OptionHierarchicalParameters)

    nanim = size(Clusterization)[1]


    a = Clusterization[1].prior_rho[1]
    b = Clusterization[1].prior_rho[2]
    sumM    = Float64(0.0)
    rhoapp = 0.0
    for ian in 1:nanim



        MatM        = Clusterization[ian].MatM
        MatMbar     = Clusterization[ian].MatMbar
        ktot        = Clusterization[ian].clusterization.kmax
        nt          = size(Clusterization[ian].clusterization.zeta,1)
        MatM        = Clusterization[ian].MatM
        MatMbar     = Clusterization[ian].MatMbar
        nnonempty   = Clusterization[ian].clusterization.n_nonemptyC[1]

        gammaPar    = Clusterization[ian].mcmc_gamma[1]
        rhoPar      = Clusterization[ian].mcmc_rho[1]
        kappaPar    = Clusterization[ian].mcmc_rho[1]*Clusterization[ian].mcmc_ak[1]
        alphaPar    = Clusterization[ian].mcmc_ak[1]-kappaPar
        betaPar     = Clusterization[ian].mcmc_beta

        sumM        += 1.0*sum(MatM)


        if sum(MatM) != 0
            for j in 1:sum(MatM)
                rhoapp += Float64(rand(rng,  Bernoulli( rhoPar  ) ))
            end
        end

    end
    if sumM != 0
        Clusterization[1].mcmc_rho[1] = rand(rng,Beta( rhoapp+a  ,sumM-rhoapp+b   ))
    else
        Clusterization[1].mcmc_rho[1] = rand(rng,Beta( a  ,b   ))
    end

    if nanim>1
        for ian in 2:nanim
            Clusterization[ian].mcmc_rho[1] = Clusterization[1].mcmc_rho[1]
        end
    end

    return nothing
end


function sample_rhodp_type2!(rng::MersenneTwister,Clusterization::Vector{Clusterization_HDPHMM}, piPar::VecParDirichlet,parHier::OptionHierarchicalParameters)

    sample_rhodp!(rng,Clusterization[1], Clusterization[1].pi)

    nanim = size(Clusterization)[1]
    if nanim>1
        for ian in 2:nanim
            Clusterization[ian].mcmc_rho[1] = Clusterization[1].mcmc_rho[1]
        end
    end

end


function sample_rhodp!(rng::MersenneTwister,Clusterization::Vector{Clusterization_HDPHMM_Joint}, piPar::VecParDirichlet)

    sample_rhodp!(rng,Clusterization[1], Clusterization[1].pi)

    nanim = size(Clusterization)[1]
    if nanim>1
        for ian in 2:nanim
            Clusterization[ian].mcmc_rho[1] = Clusterization[1].mcmc_rho[1]
        end
    end

end

function sample_rhodp!(rng::MersenneTwister,Clusterization::Clusterization_HDPHMM_Joint, piPar::VecParDirichlet)

    #
    #Clusterization = MCMCClusterization
    # piPar = MCMCClusterization.pi
    MatM        = Clusterization.MatM
    MatMbar     = Clusterization.MatMbar
    ktot        = Clusterization.clusterization.kmax
    nt          = size(Clusterization.clusterization.zeta,1)
    MatM        = Clusterization.MatM
    MatMbar     = Clusterization.MatMbar
    nnonempty   = Clusterization.clusterization.n_nonemptyC[1]

    gammaPar    = Clusterization.mcmc_gamma[1]
    rhoPar      = Clusterization.mcmc_rho[1]
    kappaPar    = Clusterization.mcmc_rho[1]*Clusterization.mcmc_ak[1]
    alphaPar    = Clusterization.mcmc_ak[1]-kappaPar
    betaPar     = Clusterization.mcmc_beta

    a = Clusterization.prior_rho[1]
    b = Clusterization.prior_rho[2]
    rhoapp = 0.0
    if sum(MatM) != 0
        for j in 1:sum(MatM)
            rhoapp += Float64(rand(rng,  Bernoulli( rhoPar  ) ))
        end
        Clusterization.mcmc_rho[1] = rand(rng,Beta( rhoapp+a  ,sum(MatM)-rhoapp+b   ))
    else
        Clusterization.mcmc_rho[1] = rand(rng,Beta( a  ,b   ))
    end


    return nothing
end



function sample_rhodp!(rng::MersenneTwister,Clusterization::Vector{Clusterization_HDPHMM}, piPar::VecParDirichlet)

    nanim = size(Clusterization)[1]
    for ian in 1:nanim
        #print("rho ", rand(rng),"\n")
        sample_rhodp!(rng,Clusterization[ian], Clusterization[ian].pi)
    end

end

function sample_rhodp!(rng::MersenneTwister,Clusterization::Clusterization_HDPHMM, piPar::VecParDirichlet)

    #
    #Clusterization = MCMCClusterization
    # piPar = MCMCClusterization.pi
    MatM        = Clusterization.MatM
    MatMbar     = Clusterization.MatMbar
    ktot        = Clusterization.clusterization.kmax
    nt          = size(Clusterization.clusterization.zeta,1)
    MatM        = Clusterization.MatM
    MatMbar     = Clusterization.MatMbar
    nnonempty   = Clusterization.clusterization.n_nonemptyC[1]

    gammaPar    = Clusterization.mcmc_gamma[1]
    rhoPar      = Clusterization.mcmc_rho[1]
    kappaPar    = Clusterization.mcmc_rho[1]*Clusterization.mcmc_ak[1]
    alphaPar    = Clusterization.mcmc_ak[1]-kappaPar
    betaPar     = Clusterization.mcmc_beta

    a = Clusterization.prior_rho[1]
    b = Clusterization.prior_rho[2]
    rhoapp = 0.0
    if sum(MatM) != 0
        for j in 1:sum(MatM)
            rhoapp += Float64(rand(rng,  Bernoulli( rhoPar  ) ))
        end
        Clusterization.mcmc_rho[1] = rand(rng,Beta( rhoapp+a  ,sum(MatM)-rhoapp+b   ))
    else
        Clusterization.mcmc_rho[1] = rand(rng,Beta( a  ,b   ))
    end


    return nothing
end
