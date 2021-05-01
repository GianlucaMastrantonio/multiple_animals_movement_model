
function sampleDP(rng::MersenneTwister,par1::Vector{Float64}, par2::Vector{Float64}, order::Vector{Int16}, K::Int16)

    ret = zeros(Float64,K)

    j = order[1]
    ret[j] = rand(Beta(par1[j],par2[j]))
    sum   = 1.0-ret[j]
    if K>2
        for k in 2:(K-1)

            j = order[k]
            #print(par1[j], " ",par2[j], "\n")
            ret[j] = rand(Beta(par1[j],par2[j]))*sum
            sum   -= ret[j]
        end
    end
    #print("\n")
    j = order[K]
    ret[j] = sum


    return ret
end

function sample_pi!(rng::MersenneTwister,Clusterization::AbstractClusterization, piPar::AbstractVecPar)
    error(string("sample_mu0 not defined for Likelihood ", typeof(Clusterization), " and mu0 ", typeof(mu0Par)) )
end






function sample_pi!(rng::MersenneTwister,Clusterization::Vector{Clusterization_HDPHMM_Joint}, piPar::VecParDirichlet)

    sample_pi!(rng::MersenneTwister,Clusterization[1], Clusterization[1].pi)

    nanim = size(Clusterization)[1]
    if nanim>1
        for ian in 2:nanim
            for i1 in 1:size(Clusterization[1].pi.parameteracc)[1]
                for i2 in 1:size(Clusterization[1].pi.parameteracc[1])[1]
                    Clusterization[ian].pi.parameteracc[i1][i2] = Clusterization[1].pi.parameteracc[i1][i2]
                end
            end

        end
    end


end

function sample_pi!(rng::MersenneTwister,Clusterization::Clusterization_HDPHMM_Joint, piPar::VecParDirichlet)

    #
    #Clusterization = MCMCClusterization
    # piPar = MCMCClusterization.pi
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

    # pi
    for k in 1:Clusterization.clusterization.kmax

        n_app                   = Clusterization.clusterization.n_itojC[k]
        alpha_P                 = betaPar*alphaPar.+ Float64.(n_app)
        alpha_P[k]              += kappaPar
        piPar.parameteracc[k]   =  rand(Dirichlet(alpha_P))
    end

    return nothing
end




function sample_pi_type2!(rng::MersenneTwister,Clusterization::Vector{Clusterization_HDPHMM}, piPar::VecParDirichlet)

    sample_pi!(rng::MersenneTwister,Clusterization[1], Clusterization[1].pi)

    nanim = size(Clusterization)[1]
    if nanim>1
        for ian in 2:nanim
            for i1 in 1:size(Clusterization[1].pi.parameteracc)[1]
                for i2 in 1:size(Clusterization[1].pi.parameteracc[1])[1]
                    Clusterization[ian].pi.parameteracc[i1][i2] = Clusterization[1].pi.parameteracc[i1][i2]
                end
            end

        end
    end


end

function sample_pi!(rng::MersenneTwister,Clusterization::Vector{Clusterization_HDPHMM}, piPar::VecParDirichlet)

    nanim = size(Clusterization)[1]
    for ian in 1:nanim
        sample_pi!(rng::MersenneTwister,Clusterization[ian], Clusterization[ian].pi)
    end

end
function sample_pi!(rng::MersenneTwister,Clusterization::Clusterization_HDPHMM, piPar::VecParDirichlet)

    #
    #Clusterization = MCMCClusterization[2]
    #piPar = MCMCClusterization[2].pi
    ktot        = Clusterization.clusterization.kmax
    kmax = ktot
    nt          = size(Clusterization.clusterization.zeta,1)
    MatM        = Clusterization.MatM
    MatMbar     = Clusterization.MatMbar
    nnonempty   = Clusterization.clusterization.n_nonemptyC[1]

    gammaPar    = Clusterization.mcmc_gamma[1]
    rhoPar      = Clusterization.mcmc_rho[1]
    kappaPar    = Clusterization.mcmc_rho[1]*Clusterization.mcmc_ak[1]
    alphaPar    = Clusterization.mcmc_ak[1]-kappaPar
    betaPar     = Clusterization.mcmc_beta

    # pi
    for k in 1:Clusterization.clusterization.kmax

        n_app                   = Clusterization.clusterization.n_itojC[k]
        alpha_P                 = betaPar*alphaPar.+ Float64.(n_app)
        alpha_P[k]              += kappaPar

        OrderVec                = zeros(Int16,kmax)

        par1                    = zeros(Float64, kmax)
        par1[1:kmax]            = alpha_P[1:kmax]
        par2                    = zeros(Float64, kmax)

        for iord in 1:Clusterization.clusterization.n_nonemptyC[1]
            OrderVec[iord] = Clusterization.clusterization.nonemptyC[iord]
        end

        for iord in 1:Clusterization.clusterization.n_emptyC[1]
            OrderVec[Clusterization.clusterization.n_nonemptyC[1]+iord] = Clusterization.clusterization.emptyC[iord]
        end

        SumTot  = Float64(max(1.0,sum(betaPar))*alphaPar+kappaPar+sum(n_app))
        #SumTot  = max( alphaPar+kappaPar+Float64(sum(n_app)), Float64(sum(par1[1:kmax])) )
        #print("PAR \n")
        for iord in 1:kmax
            j = OrderVec[iord]

            SumTot -= par1[j]
            par2[j] = SumTot

            #print(par1[j]," ", par2[j], " ",sum(betaPar)," ",OrderVec[iord]," \n")
        end
        par2[par2.<1.0e-100] .= 1.0e-100


        piPar.parameteracc[k][1:kmax] = sampleDP(rng,par1, par2,  OrderVec, Int16(kmax))[1:kmax]
        #piPar.parameteracc[k]   =  rand(Dirichlet(alpha_P))
    end

    return nothing
end

#
# function sample_pi!(Clusterization::Clusterization_HDPHMM, piPar::VecParDirichlet)
#
#     #
#     #Clusterization = MCMCClusterization
#     # piPar = MCMCClusterization.pi
#     ktot        = Clusterization.clusterization.kmax
#     nt          = size(Clusterization.clusterization.zeta,1)
#     MatM        = zeros(Float64,ktot,ktot)
#     MatMbar     = zeros(Float64,ktot,ktot)
#     nnonempty   = Clusterization.clusterization.n_nonemptyC[1]
#
#     gammaPar    = Clusterization.mcmc_gamma[1]
#     rhoPar      = Clusterization.mcmc_rho[1]
#     kappaPar    = Clusterization.mcmc_rho[1]*Clusterization.mcmc_ak[1]
#     alphaPar    = Clusterization.mcmc_ak[1]-kappaPar
#     betaPar     = Clusterization.mcmc_beta
#
#     logrj = 0.0
#     sj    = 0.0
#
#
#     for k1 in 1:nnonempty
#         kk1 = Clusterization.clusterization.nonemptyC[k1]
#         for k2 in 1:nnonempty
#
#             kk2 = Clusterization.clusterization.nonemptyC[k2]
#             napp = 0
#             for j in 1:Clusterization.clusterization.n_itojC[kk1][kk2]
#
#                 u = rand(Uniform(0.0,1.0))
#
#                 if u<  ( alphaPar*betaPar[kk2]   +kappaPar*(k1==k2) )/( alphaPar*betaPar[kk2]   +kappaPar*(k1==k2) + napp   )
#                     MatM[kk1,kk2] += one(MatM[kk1,kk2])
#                 end
#
#                 napp += one(napp)
#             end
#         end
#
#     end
#
#     DirPar  = gammaPar*Float64.(ones(ktot)/ktot)
#     MatMbar = deepcopy(MatM)
#
#     for k1 in 1:nnonempty
#         kk1 = Clusterization.clusterization.nonemptyC[k1]
#         for j in 1:Int(MatM[kk1,kk1])
#
#             u = rand(Uniform(0.0,1.0))
#
#             if u < ( rhoPar/(rhoPar + betaPar[kk1]*(1.0-rhoPar)) )
#                 MatMbar[kk1,kk1] -= one(MatMbar[kk1,kk1])
#             end
#         end
#         DirPar[kk1] += sum(MatMbar[:,kk1])
#     end
#
#     betaPar[1:ktot]  = rand(Dirichlet(DirPar))
#     betaPar[betaPar.<1.0e-100] .= 1.0e-100
#     # pi
#     for k in 1:Clusterization.clusterization.kmax
#
#         n_app                   = Clusterization.clusterization.n_itojC[k]
#         alpha_P                 = betaPar*alphaPar.+ Float64.(n_app)
#         alpha_P[k]              += kappaPar
#         piPar.parameteracc[k]   =  rand(Dirichlet(alpha_P))
#     end
#
#     # ak
#     a = Clusterization.prior_ak[1]
#     b = Clusterization.prior_ak[2]
#     for k1 in 1:nnonempty
#         kk1 = Clusterization.clusterization.nonemptyC[k1]
#         logrj += log(rand(Beta(alphaPar +kappaPar+1.0, sum(Clusterization.clusterization.n_itojC[kk1] )  )))
#
#         sj += Float64(rand(Bernoulli(sum(Clusterization.clusterization.n_itojC[kk1] )/(sum(Clusterization.clusterization.n_itojC[kk1] )+alphaPar +kappaPar)   )))
#
#     end
#
#     Clusterization.mcmc_ak[1] = rand(Gamma(a+sum(MatM)  -sj, 1.0/( b-logrj )  ))
#     kappaPar = Clusterization.mcmc_rho[1]*Clusterization.mcmc_ak[1]
#     alphaPar = Clusterization.mcmc_ak[1]-kappaPar
#     #
#
#
#
#     # gamma
#     a = Clusterization.prior_gamma[1]
#     b = Clusterization.prior_gamma[2]
#     if sum(MatMbar)!=0
#         app1 = log(rand(Beta(gammaPar/ktot+1, sum(MatMbar) )))
#         app2 = Float64(rand(Bernoulli(sum(MatMbar)/(sum(MatMbar)+gammaPar/ktot) )))
#
#         Clusterization.mcmc_gamma[1] = rand(Gamma(a+Clusterization.clusterization.n_nonemptyC[1]-app2, ktot/( ktot*b-app1 )  ))
#     else
#         Clusterization.mcmc_gamma[1] = rand(Gamma(a+Clusterization.clusterization.n_nonemptyC[1], 1.0/b  ))
#     end
#
#     #
#     #
#     #
#     a = Clusterization.prior_rho[1]
#     b = Clusterization.prior_rho[2]
#     rhoapp = 0.0
#     if sum(MatM) != 0
#         for j in 1:sum(MatM)
#             rhoapp += Float64(rand(  Bernoulli( rhoPar  ) ))
#         end
#         Clusterization.mcmc_rho[1] = rand(Beta( rhoapp+a  ,sum(MatM)-rhoapp+b   ))
#     else
#         Clusterization.mcmc_rho[1] = rand(Beta( a  ,b   ))
#     end
#
#
#
#     return nothing
# end
