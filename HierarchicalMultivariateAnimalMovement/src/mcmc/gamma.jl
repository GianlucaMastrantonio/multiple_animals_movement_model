
function sample_gamma!(rng::MersenneTwister,Clusterization::Vector{Clusterization_HDPHMM}, piPar::VecParDirichlet,parHier::OptionHierarchicalParameters )


    nanim = size(Clusterization)[1]
    ktot  = Clusterization[1].clusterization.kmax
    # gamma
    a = Clusterization[1].prior_gamma[1]
    b = Clusterization[1].prior_gamma[2]

    app1 = Float64(0.0)
    app2 = Float64(0.0)
    nnon = Float64(0.0)
    for ian in 1:nanim
        #Clusterization = MCMCClusterization
        # piPar = MCMCClusterization.pi
        MatM        = Clusterization[ian].MatM
        MatMbar     = Clusterization[ian].MatMbar

        nt          = size(Clusterization[ian].clusterization.zeta,1)
        MatM        = Clusterization[ian].MatM
        MatMbar     = Clusterization[ian].MatMbar
        nnonempty   = Clusterization[ian].clusterization.n_nonemptyC[1]

        gammaPar    = Clusterization[ian].mcmc_gamma[1]
        rhoPar      = Clusterization[ian].mcmc_rho[1]
        kappaPar    = Clusterization[ian].mcmc_rho[1]*Clusterization[ian].mcmc_ak[1]
        alphaPar    = Clusterization[ian].mcmc_ak[1]-kappaPar
        betaPar     = Clusterization[ian].mcmc_beta


        if sum(MatMbar)!=0
            app1 += log(rand(Beta(gammaPar/ktot+1, sum(MatMbar) )))
            app2 += Float64(rand(Bernoulli(sum(MatMbar)/(sum(MatMbar)+gammaPar/ktot) )))

            nnon += 1.0*Clusterization[ian].clusterization.n_nonemptyC[1]


        end
    end
    Clusterization[1].mcmc_gamma[1] = rand(Gamma(a+nnon -app2, 1.0/( b-app1 )  ))
    if nanim>1
        for ian in 2:nanim
            Clusterization[ian].mcmc_gamma[1] = Clusterization[1].mcmc_gamma[1]
        end
    end

    #


    return nothing
end


function sample_gamma_type2!(rng::MersenneTwister,Clusterization::Vector{Clusterization_HDPHMM}, piPar::VecParDirichlet,parHier::OptionHierarchicalParameters)

    sample_gamma!(rng::MersenneTwister,Clusterization[1], Clusterization[1].pi)

    nanim = size(Clusterization)[1]
    if nanim>1
        for ian in 2:nanim
            Clusterization[ian].mcmc_gamma[1] = Clusterization[1].mcmc_gamma[1]
        end
    end

end




function sample_gamma!(rng::MersenneTwister,Clusterization::Vector{Clusterization_HDPHMM_Joint}, piPar::VecParDirichlet)

    sample_gamma!(rng::MersenneTwister,Clusterization[1], Clusterization[1].pi)

    nanim = size(Clusterization)[1]
    if nanim>1
        for ian in 2:nanim
            Clusterization[ian].mcmc_gamma[1] = Clusterization[1].mcmc_gamma[1]
        end
    end

end

function sample_gamma!(rng::MersenneTwister,Clusterization::Clusterization_HDPHMM_Joint, piPar::VecParDirichlet)

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

    # gamma
    a = Clusterization.prior_gamma[1]
    b = Clusterization.prior_gamma[2]
    if sum(MatMbar)!=0
        app1 = log(rand(Beta(gammaPar/ktot+1, sum(MatMbar) )))
        app2 = Float64(rand(Bernoulli(sum(MatMbar)/(sum(MatMbar)+gammaPar/ktot) )))

        Clusterization.mcmc_gamma[1] = rand(Gamma(a+Clusterization.clusterization.n_nonemptyC[1]-app2, ktot/( ktot*b-app1 )  ))
    else
        Clusterization.mcmc_gamma[1] = rand(Gamma(a+Clusterization.clusterization.n_nonemptyC[1], 1.0/b  ))
    end
    #


    return nothing
end


function sample_gamma!(rng::MersenneTwister,Clusterization::Vector{Clusterization_HDPHMM}, piPar::VecParDirichlet)

    nanim = size(Clusterization)[1]
    for ian in 1:nanim
        sample_gamma!(rng::MersenneTwister,Clusterization[ian], Clusterization[ian].pi)
    end

end

function sample_gamma!(rng::MersenneTwister,Clusterization::Clusterization_HDPHMM, piPar::VecParDirichlet)

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

    # gamma
    a = Clusterization.prior_gamma[1]
    b = Clusterization.prior_gamma[2]
    if sum(MatMbar)!=0
        app1 = log(rand(Beta(gammaPar/ktot+1, sum(MatMbar) )))
        app2 = Float64(rand(Bernoulli(sum(MatMbar)/(sum(MatMbar)+gammaPar/ktot) )))

        Clusterization.mcmc_gamma[1] = rand(Gamma(a+Clusterization.clusterization.n_nonemptyC[1]-app2, ktot/( ktot*b-app1 )  ))
    else
        Clusterization.mcmc_gamma[1] = rand(Gamma(a+Clusterization.clusterization.n_nonemptyC[1], 1.0/b  ))
    end
    #


    return nothing
end
