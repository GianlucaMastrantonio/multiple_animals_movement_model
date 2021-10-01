
function compute_m!(rng::MersenneTwister,Clusterization::Vector{Clusterization_HDPHMM_Joint}, piPar::VecParDirichlet)

    compute_m!(rng::MersenneTwister,Clusterization[1], Clusterization[1].pi)

    nanim = size(Clusterization)[1]
    if nanim>1
        for ian in 2:nanim

            for i1 in 1:size(Clusterization[1].MatM)[1]
                for i2 in 1:size(Clusterization[1].MatM)[2]
                    Clusterization[ian].MatM[i1,i2] = Clusterization[1].MatM[i1,i2]
                    Clusterization[ian].MatMbar[i1,i2] = Clusterization[1].MatMbar[i1,i2]
                end
            end
        end
    end

end

function compute_m!(rng::MersenneTwister,Clusterization::Clusterization_HDPHMM_Joint, piPar::VecParDirichlet)


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

    MatM    .= 0
    MatMbar .= 0
    for k1 in 1:nnonempty
        kk1 = Clusterization.clusterization.nonemptyC[k1]
        for k2 in 1:nnonempty

            kk2 = Clusterization.clusterization.nonemptyC[k2]
            napp = 0
            for j in 1:Clusterization.clusterization.n_itojC[kk1][kk2]

                u = rand(rng,Uniform(0.0,1.0))

                if u<  ( alphaPar*betaPar[kk2]   +kappaPar*(k1==k2) )/( alphaPar*betaPar[kk2]   +kappaPar*(k1==k2) + napp   )
                    MatM[kk1,kk2] += one(MatM[kk1,kk2])
                    MatMbar[kk1,kk2] += one(MatMbar[kk1,kk2])
                end

                napp += one(napp)
            end
        end

    end

    for k1 in 1:nnonempty
        kk1 = Clusterization.clusterization.nonemptyC[k1]
        for j in 1:Int(MatM[kk1,kk1])

            u = rand(rng,Uniform(0.0,1.0))

            if u < ( rhoPar/(rhoPar + betaPar[kk1]*(1.0-rhoPar)) )
                MatMbar[kk1,kk1] -= one(MatMbar[kk1,kk1])
            end
        end
    end


    return nothing
end



function compute_m_type2!(rng::MersenneTwister,Clusterization::Vector{Clusterization_HDPHMM}, piPar::VecParDirichlet)


    compute_m!(rng,Clusterization[1], Clusterization[1].pi)

    nanim = size(Clusterization)[1]
    if nanim>1
        for ian in 2:nanim

            for i1 in 1:size(Clusterization[1].MatM)[1]
                for i2 in 1:size(Clusterization[1].MatM)[2]
                    Clusterization[ian].MatM[i1,i2] = Clusterization[1].MatM[i1,i2]
                    Clusterization[ian].MatMbar[i1,i2] = Clusterization[1].MatMbar[i1,i2]
                end
            end
        end
    end

end


function compute_m!(rng::MersenneTwister,Clusterization::Vector{Clusterization_HDPHMM}, piPar::VecParDirichlet)

    nanim = size(Clusterization)[1]
    for ian in 1:nanim
        compute_m!(rng,Clusterization[ian], Clusterization[ian].pi)
    end

end

function compute_m!(rng::MersenneTwister,Clusterization::Clusterization_HDPHMM, piPar::VecParDirichlet)


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

    MatM    .= 0
    MatMbar .= 0
    for k1 in 1:nnonempty
        kk1 = Clusterization.clusterization.nonemptyC[k1]
        for k2 in 1:nnonempty

            kk2 = Clusterization.clusterization.nonemptyC[k2]
            napp = 0
            for j in 1:Clusterization.clusterization.n_itojC[kk1][kk2]

                u = rand(rng,Uniform(0.0,1.0))

                if u<  ( alphaPar*betaPar[kk2]   +kappaPar*(k1==k2) )/( alphaPar*betaPar[kk2]   +kappaPar*(k1==k2) + napp   )
                    MatM[kk1,kk2] += one(MatM[kk1,kk2])
                    MatMbar[kk1,kk2] += one(MatMbar[kk1,kk2])
                end

                napp += one(napp)
            end
        end

    end

    for k1 in 1:nnonempty
        kk1 = Clusterization.clusterization.nonemptyC[k1]
        for j in 1:Int(MatM[kk1,kk1])

            u = rand(rng,Uniform(0.0,1.0))

            if u < ( rhoPar/(rhoPar + betaPar[kk1]*(1.0-rhoPar)) )
                MatMbar[kk1,kk1] -= one(MatMbar[kk1,kk1])
            end
        end
    end


    return nothing
end




function compute_m_singleHDP!(rng::MersenneTwister,Clusterization::Vector{Clusterization_HDPHMM}, piPar::VecParDirichlet)

    nanim = size(Clusterization)[1]
    for ian in 1:nanim
        compute_m_singleHDP!(rng,Clusterization[ian], Clusterization[ian].pi)
    end

end

function compute_m_singleHDP!(rng::MersenneTwister,Clusterization::Clusterization_HDPHMM, piPar::VecParDirichlet)


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

    MatM    .= 0
    MatMbar .= 0
    for k1 in 1:nnonempty
        kk1 = Clusterization.clusterization.nonemptyC[k1]
        for k2 in 1:nnonempty

            kk2 = Clusterization.clusterization.nonemptyC[k2]
            napp = 0
            for j in 1:Clusterization.clusterization.n_itojC[kk1][kk2]

                u = rand(rng,Uniform(0.0,1.0))

                if u<  ( alphaPar*betaPar[kk2]   +kappaPar*(k1==k2) )/( alphaPar*betaPar[kk2]   +kappaPar*(k1==k2) + napp   )
                    MatM[kk1,kk2] += one(MatM[kk1,kk2])
                    MatMbar[kk1,kk2] += one(MatMbar[kk1,kk2])
                end

                napp += one(napp)
            end
        end

    end

    for k1 in 1:nnonempty
        kk1 = Clusterization.clusterization.nonemptyC[k1]
        for j in 1:Int(MatM[kk1,kk1])

            u = rand(rng,Uniform(0.0,1.0))

            if u < ( rhoPar/(rhoPar + betaPar[kk1]*(1.0-rhoPar)) )
                MatMbar[kk1,kk1] -= one(MatMbar[kk1,kk1])
            end
        end
    end


    return nothing
end
