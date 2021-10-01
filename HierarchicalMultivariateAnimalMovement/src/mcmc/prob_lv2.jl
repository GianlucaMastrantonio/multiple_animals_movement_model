function sample_prob_lv2!(rng::MersenneTwister,Clusterization::Vector{Clusterization_HDPHMM}, piPar::VecParDirichlet, parHier::OptionHierarchicalParameters)
    ## My model

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
    DirPar_mu       = gammaPar*Float64.(zeros(ktot)/ktot)
    DirPar_nu       = gammaPar*Float64.(zeros(ktot)/ktot)
    DirPar_eta      = gammaPar*Float64.(zeros(ktot)/ktot)
    DirPar_sigma    = gammaPar*Float64.(zeros(ktot)/ktot)
    DirPar_rho      = gammaPar*Float64.(zeros(ktot)/ktot)

    for ian in 1:nanim
        MatMsum         = sum(Clusterization[ian].MatMbar,dims=1)

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

    ### sample gamma
    kmu = Float64(0)
    keta = Float64(0)
    ksigma = Float64(0)
    knu = Float64(0)
    krho = Float64(0)
    summu = Float64(0)
    sumeta = Float64(0)
    sumsigma = Float64(0)
    sumnu = Float64(0)
    sumrho = Float64(0)
    for k1 in 1:ktot
        if DirPar_mu[k1]>0
            kmu = kmu+Float64(1)
            summu = summu+DirPar_mu[k1]
        end
        if DirPar_eta[k1]>0
            keta = keta+Float64(1)
            sumeta = sumeta+DirPar_eta[k1]
        end
        if DirPar_sigma[k1]>0
            ksigma = ksigma+Float64(1)
            sumsigma = sumsigma+DirPar_sigma[k1]
        end
        if DirPar_nu[k1]>0
            knu = knu+Float64(1)
            sumnu = sumnu+DirPar_nu[k1]
        end
        if DirPar_rho[k1]>0
            krho = krho+Float64(1)
            sumrho = sumrho+DirPar_rho[k1]
        end
    end


    # eta = rand(Beta(gamma_par[1]+1, sum(MatMbar) ))
    # gamma_par
    # gamma_mu = rand(Gamma(a+kmu, 1/( b-app1 )  ))
    a = Clusterization[1].prior_gamma[1]
    b = Clusterization[1].prior_gamma[2]

    if summu!=0
        app1 = log(rand(Beta(Clusterization[1].gamma_par[1]+1, summu )))
        app2 = Float64(rand(Bernoulli(summu/(summu+Clusterization[1].gamma_par[1]) )))

        Clusterization[1].gamma_par[1] = rand(Gamma(a+kmu-app2, 1.0/( b-app1 )  ))
    else
        Clusterization[1].gamma_par[1] = rand(Gamma(a+kmu, 1.0/b  ))
    end
    if sumeta!=0
        app1 = log(rand(Beta(Clusterization[1].gamma_par[2]+1, sumeta )))
        app2 = Float64(rand(Bernoulli(sumeta/(sumeta+Clusterization[1].gamma_par[2]) )))

        Clusterization[1].gamma_par[2] = rand(Gamma(a+keta-app2, 1.0/( b-app1 )  ))
    else
        Clusterization[1].gamma_par[2] = rand(Gamma(a+keta, 1.0/b  ))
    end
    if sumsigma!=0
        app1 = log(rand(Beta(Clusterization[1].gamma_par[3]+1, sumsigma )))
        app2 = Float64(rand(Bernoulli(sumsigma/(sumsigma+Clusterization[1].gamma_par[3]) )))

        Clusterization[1].gamma_par[3] = rand(Gamma(a+ksigma-app2, 1.0/( b-app1 )  ))
    else
        Clusterization[1].gamma_par[3] = rand(Gamma(a+ksigma, 1.0/b  ))
    end

    if sumnu!=0
        app1 = log(rand(Beta(Clusterization[1].gamma_par[4]+1, sumnu )))
        app2 = Float64(rand(Bernoulli(sumnu/(sumnu+Clusterization[1].gamma_par[4]) )))

        Clusterization[1].gamma_par[4] = rand(Gamma(a+knu-app2, 1.0/( b-app1 )  ))
    else
        Clusterization[1].gamma_par[4] = rand(Gamma(a+knu, 1.0/b  ))
    end

    if sumrho!=0
        app1 = log(rand(Beta(Clusterization[1].gamma_par[5]+1, sumrho )))
        app2 = Float64(rand(Bernoulli(sumrho/(sumrho+Clusterization[1].gamma_par[5]) )))

        Clusterization[1].gamma_par[5] = rand(Gamma(a+krho-app2, 1.0/( b-app1 )  ))
    else
        Clusterization[1].gamma_par[5] = rand(Gamma(a+krho, 1.0/b  ))
    end

    for k1 in 1:ktot
        DirPar_mu[k1] +=    Clusterization[1].gamma_par[1]*Float64.(1.0/ktot)
        DirPar_eta[k1] +=   Clusterization[1].gamma_par[2]*Float64.(1.0/ktot)
        DirPar_sigma[k1] += Clusterization[1].gamma_par[3]*Float64.(1.0/ktot)
        DirPar_nu[k1] +=    Clusterization[1].gamma_par[4]*Float64.(1.0/ktot)
        DirPar_rho[k1] +=   Clusterization[1].gamma_par[5]*Float64.(1.0/ktot)
    end

    #### end sample gamma


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


# function sample_prob_lv2!(rng::MersenneTwister,Clusterization::Vector{Clusterization_HDPHMM}, piPar::VecParDirichlet, parHier::OptionHierarchicalParameters)
#  ## APPROCCIO MIO FUNZIOANNTE
#     #sample_beta!(Clusterization[1], Clusterization[1].pi)
#
#     kmax  = Clusterization[1].clusterization.kmax
#     nanim = size(Clusterization)[1]
#     ktot  = kmax
#
#     # #Clusterization = MCMCClusterization
#     # # piPar = MCMCClusterization.pi
#
#     # ktot        = Clusterization.clusterization.kmax
#     # nt          = size(Clusterization.clusterization.zeta,1)
#     # # MatM        = zeros(Float64,ktot,ktot)
#     # # MatMbar     = zeros(Float64,ktot,ktot)
#     # nnonempty   = Clusterization.clusterization.n_nonemptyC[1]
#     #
#     # gammaPar    = Clusterization.mcmc_gamma[1]
#     # rhoPar      = Clusterization.mcmc_rho[1]
#     # kappaPar    = Clusterization.mcmc_rho[1]*Clusterization.mcmc_ak[1]
#     # alphaPar    = Clusterization.mcmc_ak[1]-kappaPar
#     # betaPar     = Clusterization.mcmc_beta
#     #
#     # logrj = 0.0
#     # sj    = 0.0
#     #
#     #      DirPar[kk1] += sum(MatMbar[:,kk1])
#     #     end
#     #
#     #     betaPar[1:ktot]  = rand(Dirichlet(DirPar))
#     #
#     gammaPar        = Clusterization[1].mcmc_gamma[1]
#     DirPar_mu       = gammaPar*Float64.(ones(ktot)/ktot)
#     DirPar_nu       = gammaPar*Float64.(ones(ktot)/ktot)
#     DirPar_eta      = gammaPar*Float64.(ones(ktot)/ktot)
#     DirPar_sigma    = gammaPar*Float64.(ones(ktot)/ktot)
#     DirPar_rho      = gammaPar*Float64.(ones(ktot)/ktot)
#
#     for ian in 1:nanim
#         MatMsum         = sum(Clusterization[ian].MatMbar,dims=1)
#
#         #for k in 1:kmax
#         for k1 in 1:Clusterization[ian].clusterization.n_nonemptyC[1]
#             k = Clusterization[ian].clusterization.nonemptyC[k1]
#
#             j = parHier.h_mu.clust[k,ian]
#             DirPar_mu[j] += MatMsum[k]
#
#             j = parHier.h_nu.clust[k,ian]
#             DirPar_nu[j] += MatMsum[k]
#
#             j = parHier.h_eta.clust[k,ian]
#             DirPar_eta[j] += MatMsum[k]
#
#             j = parHier.h_sigma.clust[k,ian]
#             DirPar_sigma[j] += MatMsum[k]
#
#             j = parHier.h_rho.clust[k,ian]
#             DirPar_rho[j] += MatMsum[k]
#
#         end
#
#     end
#     parHier.h_mu.prob[1:kmax]  = rand(rng,Dirichlet(DirPar_mu))
#     parHier.h_mu.prob[parHier.h_mu.prob.<1.0e-100] .= 1.0e-100
#
#     parHier.h_nu.prob[1:kmax]  = rand(rng,Dirichlet(DirPar_nu))
#     parHier.h_nu.prob[parHier.h_nu.prob.<1.0e-100] .= 1.0e-100
#
#     parHier.h_eta.prob[1:kmax]  = rand(rng,Dirichlet(DirPar_eta))
#     parHier.h_eta.prob[parHier.h_eta.prob.<1.0e-100] .= 1.0e-100
#
#     parHier.h_sigma.prob[1:kmax]  = rand(rng,Dirichlet(DirPar_sigma))
#     parHier.h_sigma.prob[parHier.h_sigma.prob.<1.0e-100] .= 1.0e-100
#
#     parHier.h_rho.prob[1:kmax]  = rand(rng,Dirichlet(DirPar_rho))
#     parHier.h_rho.prob[parHier.h_rho.prob.<1.0e-100] .= 1.0e-100
#
#     return nothing
# end

function sample_prob_lv2_singleHDP!(rng::MersenneTwister,Clusterization::Vector{Clusterization_HDPHMM}, piPar::VecParDirichlet, parHier::OptionHierarchicalParameters)

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
    DirPar_mu       = gammaPar*Float64.(zeros(ktot)/ktot)
    # DirPar_nu       = gammaPar*Float64.(ones(ktot)/ktot)
    # DirPar_eta      = gammaPar*Float64.(ones(ktot)/ktot)
    # DirPar_sigma    = gammaPar*Float64.(ones(ktot)/ktot)
    # DirPar_rho      = gammaPar*Float64.(ones(ktot)/ktot)

    for ian in 1:nanim
        MatMsum         = sum(Clusterization[ian].MatMbar,dims=1)

        #for k in 1:kmax
        for k1 in 1:Clusterization[ian].clusterization.n_nonemptyC[1]
            k = Clusterization[ian].clusterization.nonemptyC[k1]

            #j = parHier.h_mu.clust[k,ian]
            DirPar_mu[k] += MatMsum[k]

            # #j = parHier.h_nu.clust[k,ian]
            # DirPar_nu[k] += MatMsum[k]
            #
            # #j = parHier.h_eta.clust[k,ian]
            # DirPar_eta[k] += MatMsum[k]
            #
            # #j = parHier.h_sigma.clust[k,ian]
            # DirPar_sigma[k] += MatMsum[k]
            #
            # #j = parHier.h_rho.clust[k,ian]
            # DirPar_rho[k] += MatMsum[k]

        end

    end

    ### sample gamma
    kmu = Float64(0)

    summu = Float64(0)
    for k1 in 1:ktot
        if DirPar_mu[k1]>0
            kmu = kmu+Float64(1)
            summu = summu+DirPar_mu[k1]
        end
    end

    a = Clusterization[1].prior_gamma[1]
    b = Clusterization[1].prior_gamma[2]

    if summu!=0
        app1 = log(rand(Beta(Clusterization[1].gamma_par[1]+1, summu )))
        app2 = Float64(rand(Bernoulli(summu/(summu+Clusterization[1].gamma_par[1]) )))

        Clusterization[1].gamma_par[1] = rand(Gamma(a+kmu-app2, 1.0/( b-app1 )  ))
    else
        Clusterization[1].gamma_par[1] = rand(Gamma(a+kmu, 1.0/b  ))
    end


    for k1 in 1:ktot
        DirPar_mu[k1] +=    Clusterization[1].gamma_par[1]*Float64.(1.0/ktot)
    end



    parHier.h_mu.prob[1:kmax]  = rand(rng,Dirichlet(DirPar_mu))
    parHier.h_mu.prob[parHier.h_mu.prob.<1.0e-100] .= 1.0e-100

    # parHier.h_nu.prob[1:kmax]  = parHier.h_mu.prob[1:kmax]
    #
    # parHier.h_eta.prob[1:kmax]  = parHier.h_mu.prob[1:kmax]
    #
    # parHier.h_sigma.prob[1:kmax]  = parHier.h_mu.prob[1:kmax]
    #
    # parHier.h_rho.prob[1:kmax]  = parHier.h_mu.prob[1:kmax]

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
        MatMsum         = sum(Clusterization[ian].MatMbar,dims=1)

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
