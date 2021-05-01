

function sample_ak!(rng::MersenneTwister,Clusterization::Vector{Clusterization_HDPHMM}, piPar::VecParDirichlet,parHier::OptionHierarchicalParameters)

    nanim = size(Clusterization)[1]

    logrj = 0.0
    sj    = 0.0
    a       = Clusterization[1].prior_ak[1]
    b       = Clusterization[1].prior_ak[2]
    sumM    = Float64(0.0)
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



        # ak

        for k1 in 1:nnonempty
            kk1 = Clusterization[ian].clusterization.nonemptyC[k1]
            if sum(Clusterization[ian].clusterization.n_itojC[kk1] ) != 0
                logrj += log(rand(rng,Beta(alphaPar +kappaPar+1.0, sum(Clusterization[ian].clusterization.n_itojC[kk1] )  )))

                sj += Float64(rand(rng,Bernoulli(sum(Clusterization[ian].clusterization.n_itojC[kk1] )/(sum(Clusterization[ian].clusterization.n_itojC[kk1] )+alphaPar +kappaPar)   )))
            end


        end
        sumM += 1.0*sum(MatM)
    end
    #Clusterization.mcmc_ak[1] = rand(Gamma(a+sum(MatM)  -sj, 1.0/( b-logrj )  ))
    Clusterization[1].mcmc_ak[1] = rand(rng,Gamma(a+sumM  -sj, 1.0/(b-logrj )  ))

    if nanim>1
        for ian in 2:nanim
            Clusterization[ian].mcmc_ak[1] = Clusterization[1].mcmc_ak[1]
        end
    end


    #


    return nothing



end




function sample_ak_type2!(rng::MersenneTwister,Clusterization::Vector{Clusterization_HDPHMM}, piPar::VecParDirichlet,parHier::OptionHierarchicalParameters)

    sample_ak!(rng,Clusterization[1], Clusterization[1].pi)

    nanim = size(Clusterization)[1]
    if nanim>1
        for ian in 2:nanim
            Clusterization[ian].mcmc_ak[1] = Clusterization[1].mcmc_ak[1]
        end
    end

end



function sample_ak!(rng::MersenneTwister,Clusterization::Vector{Clusterization_HDPHMM_Joint}, piPar::VecParDirichlet)

    sample_ak!(Clusterization[1], Clusterization[1].pi)

    nanim = size(Clusterization)[1]
    if nanim>1
        for ian in 2:nanim
            Clusterization[ian].mcmc_ak[1] = Clusterization[1].mcmc_ak[1]
        end
    end

end

function sample_ak!(rng::MersenneTwister,Clusterization::Clusterization_HDPHMM_Joint, piPar::VecParDirichlet)

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

    logrj = 0.0
    sj    = 0.0

    # ak
    a = Clusterization.prior_ak[1]
    b = Clusterization.prior_ak[2]
    for k1 in 1:nnonempty
        kk1 = Clusterization.clusterization.nonemptyC[k1]
        if sum(Clusterization.clusterization.n_itojC[kk1] ) != 0
            logrj += log(rand(rng,Beta(alphaPar +kappaPar+1.0, sum(Clusterization.clusterization.n_itojC[kk1] )  )))

            sj += Float64(rand(rng,Bernoulli(sum(Clusterization.clusterization.n_itojC[kk1] )/(sum(Clusterization.clusterization.n_itojC[kk1] )+alphaPar +kappaPar)   )))
        end


    end
    # print(a)
    # print("\n")
    # print(b)
    # print("\n")
    # print(sj)
    # print("\n")
    # print(logrj)
    # print("\n")
    # print(sum(MatM))
    # print("\n")

    Clusterization.mcmc_ak[1] = rand(rng,Gamma(a+sum(MatM)  -sj, 1.0/( b-logrj )  ))
    kappaPar = Clusterization.mcmc_rho[1]*Clusterization.mcmc_ak[1]
    alphaPar = Clusterization.mcmc_ak[1]-kappaPar
    #


    return nothing
end


function sample_ak!(rng::MersenneTwister,Clusterization::Vector{Clusterization_HDPHMM}, piPar::VecParDirichlet)

    nanim = size(Clusterization)[1]
    for ian in 1:nanim
        sample_ak!(rng,Clusterization[ian], Clusterization[ian].pi)
    end

end

function sample_ak!(rng::MersenneTwister,Clusterization::Clusterization_HDPHMM, piPar::VecParDirichlet)

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

    logrj = 0.0
    sj    = 0.0

    # ak
    a = Clusterization.prior_ak[1]
    b = Clusterization.prior_ak[2]
    for k1 in 1:nnonempty
        kk1 = Clusterization.clusterization.nonemptyC[k1]
        if sum(Clusterization.clusterization.n_itojC[kk1] ) != 0
            logrj += log(rand(rng,Beta(alphaPar +kappaPar+1.0, sum(Clusterization.clusterization.n_itojC[kk1] )  )))

            sj += Float64(rand(rng,Bernoulli(sum(Clusterization.clusterization.n_itojC[kk1] )/(sum(Clusterization.clusterization.n_itojC[kk1] )+alphaPar +kappaPar)   )))
        end


    end
    # print(a)
    # print("\n")
    # print(b)
    # print("\n")
    # print(sj)
    # print("\n")
    # print(logrj)
    # print("\n")
    # print(sum(MatM))
    # print("\n")

    Clusterization.mcmc_ak[1] = rand(rng,Gamma(a+sum(MatM)  -sj, 1.0/( b-logrj )  ))
    kappaPar = Clusterization.mcmc_rho[1]*Clusterization.mcmc_ak[1]
    alphaPar = Clusterization.mcmc_ak[1]-kappaPar
    #


    return nothing
end
