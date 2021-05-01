


####################################
##### PROPOSED MODELS
####################################

function sample_zeta!(rng::MersenneTwister,Likelihood::Vector{Likelihood_OU_CircLinmodel}, Clusterization::Vector{Clusterization_HDPHMM},doupdate::ZetaDoNotUpdate)



end



function sample_zeta!(rng::MersenneTwister,Likelihood::Vector{Likelihood_OU_CircLinmodel}, Clusterization::Vector{Clusterization_HDPHMM},doupdate::ZetaDoUpdate)

    nanim = size(Likelihood)[1]
    for ian in 1:nanim
        sample_zeta!(rng,Likelihood[ian], Clusterization[ian])
    end

end
function sample_zeta!(rng::MersenneTwister,Likelihood::AbstractLikelihood, Clusterization::AbstractClusterization_Divided)

    u = rand(rng,Uniform(0.0,1.0))
    if u<0.95
        sample_zeta_beam!(rng,Likelihood,Clusterization)
    else
        sample_zeta_marg!(rng,Likelihood,Clusterization)
    end

end


###### BEAM SAMPLER

function sample_zeta_beam!(rng::MersenneTwister,Likelihood::AbstractLikelihood, Clusterization::Clusterization_HDPHMM)

    # Likelihood      = MCMCLikelihood
    # Clusterization  = MCMCClusterization
Moltbeam = Float64(1.0)
#    Clusterization.pi.parameteracc
    kmax        = Likelihood.kmax
    nc          = Likelihood.data.ncol
    nanim       = Likelihood.data.nanimals
    nt          = Likelihood.data.nt

    zetaprex = 1
    veckacc = zeros(Int16,kmax)
    probsVect       = zeros(Float64,kmax)


    Clust = Clusterization.clusterization
    kappaPar = Clusterization.mcmc_ak[1]*Clusterization.mcmc_rho[1]
    alphaPar = Clusterization.mcmc_ak[1]-kappaPar




    un    = zeros(Float64,nt)
    zetaprex = 1
    for i  in 1:(nt-1)
        @inbounds k1          = Likelihood.clusterization.zeta[i]
        @inbounds un[i]       = rand(rng,Uniform(0.0, Clusterization.pi.parameteracc[zetaprex][k1]/Moltbeam))
        @inbounds zetaprex    = k1
    end

    zetaprex = 1
    @inbounds for i  in 1:(nt-2)
        #print("i=",i,"\n")
        nk              = 0
         @inbounds zetasux        = Likelihood.clusterization.zeta[i+1]
         @inbounds kold           = Likelihood.clusterization.zeta[i]
        for k in 1:kmax

            if (Clusterization.pi.parameteracc[zetaprex][k]>un[i]) && (Clusterization.pi.parameteracc[k][zetasux]>un[i+1])
                nk            += one(nk)
                @inbounds veckacc[nk]   = k

                @inbounds probsVect[nk] = compute_logpdf(Likelihood, Int16(k), Int16(i+1)) #comp




            end
        end

        @inbounds sampvec = probsVect[1:nk]
        @inbounds Likelihood.clusterization.zeta[i] = veckacc[sample_discretevar(rng,sampvec)]
        @inbounds zetaprex = Likelihood.clusterization.zeta[i]

    end
    i       = nt-1
    nk      = 0
    kold           = Likelihood.clusterization.zeta[i]
    for k in 1:kmax

        if Clusterization.pi.parameteracc[zetaprex][k]>un[i]
            nk           += one(nk)
            veckacc[nk]   = k
            @inbounds probsVect[nk] = compute_logpdf(Likelihood, Int16(k), Int16(i+1)) # # compute_logpdf(Likelihood, Int16(k), Int16(i+1))



        end

    end
    sampvec = probsVect[1:nk]
    @inbounds Likelihood.clusterization.zeta[i] = veckacc[sample_discretevar(rng,sampvec)]

    Update_ObsInClust!(Likelihood.clusterization)
    return nothing

end

##### MARGINALIZATION

function sample_zeta_marg!(rng::MersenneTwister,Likelihood::AbstractLikelihood, Clusterization::Clusterization_HDPHMM)

    #Likelihood      = MCMCLikelihood
    # Clusterization  = MCMCClusterization

    Clust = Clusterization.clusterization
#    Clusterization.pi.parameteracc
    kmax        = Likelihood.kmax
    nc          = Likelihood.data.ncol
    nanim       = Likelihood.data.nanimals
    nt          = Likelihood.data.nt

    zetaprex = 1
    veckacc = zeros(Int16,kmax)
    probsVect       = zeros(Float64,kmax)



    # LogPdfMatrix = zeros(Float64,nt,kmax)
    # for i  in 1:(nt-1)
    #     for k in 1:kmax
    #         @inbounds LogPdfMatrix[i+1,k] = compute_logpdf(Likelihood, Int16(k), Int16(i+1))
    #     end
    # end

    #
    # un    = zeros(Float64,nt)
    # zetaprex = 1
    # for i  in 1:(nt-1)
    #     @inbounds k1          = Likelihood.clusterization.zeta[i]
    #     @inbounds un[i]       = rand(Uniform(0.0, Clusterization.pi.parameteracc[zetaprex][k1]))
    #     @inbounds zetaprex    = k1
    # end
    kappaPar = Clusterization.mcmc_ak[1]*Clusterization.mcmc_rho[1]
    alphaPar = Clusterization.mcmc_ak[1]-kappaPar
    zetaprex = 1
    @inbounds for i  in 1:(nt-2)

        @inbounds zetasux        = Likelihood.clusterization.zeta[i+1]
        @inbounds kold           = Likelihood.clusterization.zeta[i]
        @inbounds Clust.n_itojC[zetaprex][kold] -= 1
        @inbounds Clust.n_itojC[kold][zetasux]  -= 1
        for k in 1:kmax

            @inbounds probsVect[k] = compute_logpdf(Likelihood, Int16(k), Int16(i+1))

            if k == zetaprex
                @inbounds probsVect[k] += log(alphaPar*Clusterization.mcmc_beta[k]+ Clust.n_itojC[zetaprex][k] + kappaPar)
            else
                @inbounds probsVect[k] += log(alphaPar*Clusterization.mcmc_beta[k]+ Clust.n_itojC[zetaprex][k])
            end

            if k == zetasux
                @inbounds probsVect[k] += log(alphaPar*Clusterization.mcmc_beta[zetasux]+ Clust.n_itojC[k][zetasux] + kappaPar)
            else
                @inbounds probsVect[k] += log(alphaPar*Clusterization.mcmc_beta[zetasux]+ Clust.n_itojC[k][zetasux])
            end
        end
    #    print("probsVect \n", probsVect ,"\n")
        @inbounds Likelihood.clusterization.zeta[i] = sample_discretevar(rng,probsVect)

        @inbounds Clust.n_itojC[zetaprex][Likelihood.clusterization.zeta[i]] += 1
        @inbounds Clust.n_itojC[Likelihood.clusterization.zeta[i]][zetasux]  += 1

        @inbounds zetaprex = Likelihood.clusterization.zeta[i]

    end
    i       = nt-1
    nk      = 0
    kold           = Likelihood.clusterization.zeta[i]
    Clust.n_itojC[zetaprex][kold] -= 1
    for k in 1:kmax

        @inbounds probsVect[k] = compute_logpdf(Likelihood, Int16(k), Int16(i+1))

        if k == zetaprex
            @inbounds probsVect[k] += log(alphaPar*Clusterization.mcmc_beta[k]+ Clust.n_itojC[zetaprex][k] + kappaPar)
        else
            @inbounds probsVect[k] += log(alphaPar*Clusterization.mcmc_beta[k]+ Clust.n_itojC[zetaprex][k])
        end


    end

    Likelihood.clusterization.zeta[i] = sample_discretevar(rng,probsVect)

    Update_ObsInClust!(Likelihood.clusterization)
    return nothing

end






function sample_zeta_mergesplit!(rng::MersenneTwister,Likelihood::Vector{Likelihood_OU_CircLinmodel}, Clusterization::Vector{Clusterization_HDPHMM},doupdate::ZetaDoUpdate)

    #print("splitmerge","\n")
    nanim = size(Likelihood)[1]
    for ian in 1:nanim
        sample_zeta_mergesplit!(rng,Likelihood[ian], Clusterization[ian])
    end

end


function sample_zeta_mergesplit!(rng::MersenneTwister,Likelihood::AbstractLikelihood, Clusterization::Clusterization_HDPHMM)
    #print("splitmerge","\n")
    #Likelihood      = MCMCLikelihood
    # Clusterization  = MCMCClusterization

    Clust = Clusterization.clusterization
#    Clusterization.pi.parameteracc
    kmax        = Likelihood.kmax
    nc          = Likelihood.data.ncol
    nanim       = Likelihood.data.nanimals
    nt          = Likelihood.data.nt

    zetaprex = 1
    veckacc = zeros(Int16,kmax)
    probsVect       = zeros(Float64,kmax)

    kappaPar = Clusterization.mcmc_ak[1]*Clusterization.mcmc_rho[1]
    alphaPar = Clusterization.mcmc_ak[1]-kappaPar
    zetaprex = 1


    ksamp1_app     = Int16(rand(1:Likelihood.clusterization.n_nonemptyC[1],1)[1])
    ksamp1         = Int16(Likelihood.clusterization.nonemptyC[ksamp1_app])
    domerge        = Int16(0)

    ProbPropSplit = Float64(0.5)
    ProbPropMerge = Float64(1.0-ProbPropSplit)

    proposal_split = Float64(0.0);
    proposal_merge = Float64(0.0);

    uu = rand(rng, Uniform(0.0,1.0))
    if uu<ProbPropMerge
        u = Int16(1)
    else
        u = Int16(2)
    end


    if Likelihood.clusterization.n_nonemptyC[1] == 1
        u = 2
        proposal_merge += log(Float64(ProbPropMerge))
        proposal_split += log(Float64(1.0))

    elseif Likelihood.clusterization.n_nonemptyC[1] == kmax
        u = 1
        proposal_merge += log(Float64(1.0))
        proposal_split += log(Float64(ProbPropSplit))

    elseif u == 1
        proposal_merge += log(Float64(ProbPropMerge))
        proposal_split += log(Float64(ProbPropSplit))
    else
        proposal_merge += log(Float64(ProbPropMerge))
        proposal_split += log(Float64(ProbPropSplit))
    end

    if u==1
        ## MERGE
        domerge        = Int16(1)
        ksamp2_app     = Int16(rand(1:(Likelihood.clusterization.n_nonemptyC[1]-1),1)[1])
        if ksamp2_app >= ksamp1_app
            ksamp2_app = ksamp2_app+one(ksamp2_app)
        end
        ksamp2         = Int16(Likelihood.clusterization.nonemptyC[ksamp2_app])

        t_index_split1 = Int16.(findall(x->x == ksamp1, Likelihood.clusterization.zeta[1:(nt-1)]))
        t_index_split2 = Int16.(findall(x->x == ksamp2, Likelihood.clusterization.zeta[1:(nt-1)]))

        t_index_merge = sort([t_index_split1;t_index_split2])::Vector{Int16}


        # proposal
        proposal_merge += -log(Float64(Likelihood.clusterization.n_nonemptyC[1]))-log(Float64(Likelihood.clusterization.n_nonemptyC[1]-1))

        proposal_split += -log(Float64(Likelihood.clusterization.n_nonemptyC[1]-1))-log(Float64(Likelihood.clusterization.n_emptyC[1]+1))

        proposal_split += -Float64(Likelihood.clusterization.n_nonemptyC[1]-1)*log(2.0)

        ## n i to j
        nitoj_split1 = deepcopy(Clust.n_itojC[ksamp1])
        nitoj_split2 = deepcopy(Clust.n_itojC[ksamp2])

        nitoj_merge     = deepcopy(Clust.n_itojC[ksamp1])
        nitoj_merge_app = deepcopy(Clust.n_itojC[ksamp2])

        nitoj_merge[ksamp1]     += nitoj_merge_app[ksamp2]
        nitoj_merge[ksamp1]     += nitoj_merge_app[ksamp1]
        nitoj_merge[ksamp1]     += nitoj_merge[ksamp2]
        nitoj_merge_app[ksamp2]  = 0
        nitoj_merge_app[ksamp1]  = 0
        nitoj_merge[ksamp2]      = 0
        nitoj_merge             .+= nitoj_merge_app


    else
        ## split
        ksamp2_app     = Int16(rand(1:Likelihood.clusterization.n_emptyC[1],1)[1])
        ksamp2         = Int16(Likelihood.clusterization.emptyC[ksamp2_app])

        t_index_merge = Int16.(findall(x->x == ksamp1, Likelihood.clusterization.zeta[1:(nt-1)]))
        sets_vec      = find_const(Likelihood.clusterization.zeta[1:(nt-1)],ksamp1)

        t_index_split1 = zeros(Int16,length(t_index_merge))
        t_index_split2 = zeros(Int16,length(t_index_merge))

        nsplit1        = Int16(0)
        nsplit2        = Int16(0)

        ## n i to j
        nitoj_merge     = deepcopy(Clust.n_itojC[ksamp1])

        nitoj_split1    = zeros(Int16, kmax)
        nitoj_split2    = zeros(Int16, kmax)

        from_nitoj_split1    = zeros(Int16, kmax)
        from_nitoj_split2    = zeros(Int16, kmax)

        zprec = Int16(1)
        for ii in 1:size(sets_vec,1)

            if sets_vec[ii,1] != 1
                zprec = sets_vec[ii,1]-Int16(1)
            end

            u = Int16(rand(rng,1:2,1)[1])
            if u == 1

                nitoj_split1[ksamp1] += sets_vec[ii,2]-sets_vec[ii,1]
                if sets_vec[ii,2]< (nt-1)
                    nitoj_split1[Likelihood.clusterization.zeta[sets_vec[ii,2]+1]] += Int16(1)
                end

                for jj in sets_vec[ii,1]:sets_vec[ii,2]
                    nsplit1 += one(nsplit1)
                    t_index_split1[nsplit1] = jj

                end

                from_nitoj_split1[Likelihood.clusterization.zeta[zprec]] += Int16(1)

            else
                nitoj_split2[ksamp2] += sets_vec[ii,2]-sets_vec[ii,1]
                if sets_vec[ii,2]< (nt-1)
                    nitoj_split2[Likelihood.clusterization.zeta[sets_vec[ii,2]+1]] += Int16(1)
                end

                for jj in sets_vec[ii,1]:sets_vec[ii,2]
                    nsplit2 += one(nsplit2)
                    t_index_split2[nsplit2] = jj
                end

                from_nitoj_split2[Likelihood.clusterization.zeta[zprec]] += Int16(1)
            end

            #print("\n\n",sum(sets_vec),"\n",sets_vec,"\n\n")
        end
        #print("\n n", nsplit1," ",nsplit2," ",nsplit1+nsplit2," \n")
        t_index_split1 = t_index_split1[1:nsplit1]
        t_index_split2 = t_index_split2[1:nsplit2]
        # proposal
        proposal_merge += -log(Float64(Likelihood.clusterization.n_nonemptyC[1]+1))-log(Float64(Likelihood.clusterization.n_nonemptyC[1]))

        proposal_split += -log(Float64(Likelihood.clusterization.n_nonemptyC[1]))-log(Float64(Likelihood.clusterization.n_emptyC[1]))

        proposal_split += -Float64(Likelihood.clusterization.n_nonemptyC[1])*log(2.0)





    end
    #print("\n ", domerge," ", sum(nitoj_merge)==(sum(nitoj_split1)+sum(nitoj_split2))," ", sum(nitoj_merge)," ",sum(nitoj_split1)," ",sum(nitoj_split2), "\n")
    prob_split = 0.0;
    prob_merge = 0.0;
    # for t in t_index_split1
    #     @inbounds prob_split += compute_logpdf(Likelihood, ksamp1, Int16(t+1))
    # end
    # prob_merge = deepcopy(prob_split)

    for t in t_index_split2
        @inbounds prob_split += compute_logpdf(Likelihood, ksamp2, Int16(t+1))
        @inbounds prob_merge += compute_logpdf(Likelihood, ksamp1, Int16(t+1))
    end

    #print("\n pre ",domerge," ", prob_split," ",prob_merge,"\n")


    for kk in 1:kmax

        if kk!=ksamp1
            prob_split += loggamma(alphaPar*Clusterization.mcmc_beta[kk]+ nitoj_split1[kk])
            prob_merge += loggamma(alphaPar*Clusterization.mcmc_beta[kk]+ nitoj_merge[kk])
        else
            prob_split += loggamma(alphaPar*Clusterization.mcmc_beta[kk]+ nitoj_split1[kk]+ kappaPar)
            prob_merge += loggamma(alphaPar*Clusterization.mcmc_beta[kk]+ nitoj_merge[kk]+ kappaPar)

        end
        if kk!=ksamp2
            prob_split += loggamma(alphaPar*Clusterization.mcmc_beta[kk]+ nitoj_split2[kk])
            prob_merge += loggamma(alphaPar*Clusterization.mcmc_beta[kk])
        else
            prob_split += loggamma(alphaPar*Clusterization.mcmc_beta[kk]+ nitoj_split2[kk]+ kappaPar)
            prob_merge += loggamma(alphaPar*Clusterization.mcmc_beta[kk]+ kappaPar)
        end


    end

    prob_split -= loggamma(alphaPar+ sum(nitoj_split1)+ kappaPar)
    prob_split -= loggamma(alphaPar+ sum(nitoj_split2)+ kappaPar)
    prob_merge -= loggamma(alphaPar+ sum(nitoj_merge)+ kappaPar)
    prob_merge -= loggamma(alphaPar+  kappaPar)

    if domerge == 1
        for kk in 1:kmax

            if (kk!=ksamp1) && (kk!=ksamp2)
                prob_split += loggamma(alphaPar*Clusterization.mcmc_beta[ksamp1]+ Clust.n_itojC[kk][ksamp1])
                prob_split += loggamma(alphaPar*Clusterization.mcmc_beta[ksamp2]+ Clust.n_itojC[kk][ksamp2])

                prob_merge += loggamma(alphaPar*Clusterization.mcmc_beta[ksamp1]+ Clust.n_itojC[kk][ksamp2]+ Clust.n_itojC[kk][ksamp1])
                prob_merge += loggamma(alphaPar*Clusterization.mcmc_beta[ksamp2])
            end


        end
    else
        for kk in 1:kmax

            if (kk!=ksamp1) && (kk!=ksamp2)
                prob_split += loggamma(alphaPar*Clusterization.mcmc_beta[ksamp1]+ from_nitoj_split1[kk])
                prob_split += loggamma(alphaPar*Clusterization.mcmc_beta[ksamp2]+ from_nitoj_split2[kk])

                prob_merge += loggamma(alphaPar*Clusterization.mcmc_beta[ksamp1]+ Clust.n_itojC[kk][ksamp1])
                prob_merge += loggamma(alphaPar*Clusterization.mcmc_beta[ksamp2])

                #print("\n DD ", " ",kk," ",Clust.n_itojC[kk][ksamp1]== from_nitoj_split1[kk]+from_nitoj_split2[kk], Clust.n_itojC[kk][ksamp1], " ", from_nitoj_split1[kk]," ",from_nitoj_split2[kk])
            end

        end
    end


    if domerge == 1
        #merge
        MH = (prob_merge+proposal_split)-(prob_split+proposal_merge)

        #print(MH," P ",prob_split," ",proposal_merge," ",prob_merge," " ,proposal_split, "\n")
        u = rand(rng,Uniform( 0.0,1.0))

        if u<exp(MH)
            for t in t_index_merge
                @inbounds Likelihood.clusterization.zeta[t] = ksamp1
            end
            Update_ObsInClust!(Likelihood.clusterization)
        end

    else
        #split
        MH = (prob_split+proposal_merge)-(prob_merge+proposal_split)
        #print("\n","Split",ksamp1," ",ksamp2,"\n")
        #print(MH," P ",prob_split," ",proposal_merge, " ",prob_merge," " ,proposal_split, "\n")
        u = rand(rng,Uniform( 0.0,1.0))

        if u<exp(MH)
            for t in t_index_split1
                @inbounds Likelihood.clusterization.zeta[t] = ksamp1
            end
            for t in t_index_split2
                @inbounds Likelihood.clusterization.zeta[t] = ksamp2
            end
            Update_ObsInClust!(Likelihood.clusterization)
        end
    end
    #print(exp(MH),"\n")
end





function sample_zeta_mergesplit_v2!(rng::MersenneTwister,Likelihood::Vector{Likelihood_OU_CircLinmodel}, Clusterization::Vector{Clusterization_HDPHMM},doupdate::ZetaDoUpdate,parHier::OptionHierarchicalParameters)

    #print("splitmerge","\n")
    nanim = size(Likelihood)[1]
    for ian in 1:nanim
        sample_zeta_mergesplit_v2!(rng,Likelihood[ian], Clusterization[ian],parHier,ian)
    end

end


function sample_zeta_mergesplit_v2!(rng::MersenneTwister,Likelihood::AbstractLikelihood, Clusterization::Clusterization_HDPHMM,parHier::OptionHierarchicalParameters,ian)


    Clust = Clusterization.clusterization
#    Clusterization.pi.parameteracc
    kmax        = Likelihood.kmax
    nc          = Likelihood.data.ncol
    nanim       = Likelihood.data.nanimals
    nt          = Likelihood.data.nt

    zetaprex = 1
    veckacc = zeros(Int16,kmax)
    probsVect       = zeros(Float64,kmax)

    kappaPar = Clusterization.mcmc_ak[1]*Clusterization.mcmc_rho[1]
    alphaPar = Clusterization.mcmc_ak[1]-kappaPar
    zetaprex = 1


    ksamp1_app     = Int16(rand(1:Likelihood.clusterization.n_nonemptyC[1],1)[1])
    ksamp1         = Int16(Likelihood.clusterization.nonemptyC[ksamp1_app])
    domerge        = Int16(0)

    ProbPropSplit = Float64(0.5)
    ProbPropMerge = Float64(1.0-ProbPropSplit)

    proposal_split = Float64(0.0);
    proposal_merge = Float64(0.0);

    uu = rand(rng, Uniform(0.0,1.0))
    if uu<ProbPropMerge
        u = Int16(1)
    else
        u = Int16(2)
    end


    if Likelihood.clusterization.n_nonemptyC[1] == 1
        u = 2
        proposal_merge += log(Float64(ProbPropMerge))
        proposal_split += log(Float64(1.0))

    elseif Likelihood.clusterization.n_nonemptyC[1] == kmax
        u = 1
        proposal_merge += log(Float64(1.0))
        proposal_split += log(Float64(ProbPropSplit))

    elseif u == 1
        proposal_merge += log(Float64(ProbPropMerge))
        proposal_split += log(Float64(ProbPropSplit))
    else
        proposal_merge += log(Float64(ProbPropMerge))
        proposal_split += log(Float64(ProbPropSplit))
    end

    if u==1
        ## MERGE
        domerge        = Int16(1)
        ksamp2_app     = Int16(rand(1:(Likelihood.clusterization.n_nonemptyC[1]-1),1)[1])
        if ksamp2_app >= ksamp1_app
            ksamp2_app = ksamp2_app+one(ksamp2_app)
        end
        ksamp2         = Int16(Likelihood.clusterization.nonemptyC[ksamp2_app])

        t_index_split1 = Int16.(findall(x->x == ksamp1, Likelihood.clusterization.zeta[1:(nt-1)]))
        t_index_split2 = Int16.(findall(x->x == ksamp2, Likelihood.clusterization.zeta[1:(nt-1)]))

        t_index_merge = sort([t_index_split1;t_index_split2])::Vector{Int16}


        # proposal
        proposal_merge += -log(Float64(Likelihood.clusterization.n_nonemptyC[1]))-log(Float64(Likelihood.clusterization.n_nonemptyC[1]-1))

        proposal_split += -log(Float64(Likelihood.clusterization.n_nonemptyC[1]-1))

        proposal_split += -Float64(Likelihood.clusterization.n_nonemptyC[1]-1)*log(2.0)


        sumSim = 0
        sumSim += parHier.h_mu.clust[ksamp2,ian] == parHier.h_mu.clust[ksamp1,ian]
        sumSim += parHier.h_eta.clust[ksamp2,ian] == parHier.h_eta.clust[ksamp1,ian]
        sumSim += parHier.h_nu.clust[ksamp2,ian] == parHier.h_nu.clust[ksamp1,ian]
        sumSim += parHier.h_rho.clust[ksamp2,ian] == parHier.h_rho.clust[ksamp1,ian]
        sumSim += parHier.h_sigma.clust[ksamp2,ian] == parHier.h_sigma.clust[ksamp1,ian]


        proposal_split += -5.0*log(2.0)-Float64(sumSim)*log(Float64(kmax))

        ## n i to j
        nitoj_split1 = deepcopy(Clust.n_itojC[ksamp1])
        nitoj_split2 = deepcopy(Clust.n_itojC[ksamp2])

        nitoj_merge     = deepcopy(Clust.n_itojC[ksamp1])
        nitoj_merge_app = deepcopy(Clust.n_itojC[ksamp2])

        nitoj_merge[ksamp1]     += nitoj_merge_app[ksamp2]
        nitoj_merge[ksamp1]     += nitoj_merge_app[ksamp1]
        nitoj_merge[ksamp1]     += nitoj_merge[ksamp2]
        nitoj_merge_app[ksamp2]  = 0
        nitoj_merge_app[ksamp1]  = 0
        nitoj_merge[ksamp2]      = 0
        nitoj_merge             .+= nitoj_merge_app


    else
        ## split

        #ksamp2_app     = Int16(rand(1:Likelihood.clusterization.n_emptyC[1],1)[1])
        ksamp2_app     = Int16(1)
        ksamp2         = Int16(Likelihood.clusterization.emptyC[ksamp2_app])

        #### NEW
        Aindex      = Matrix{Int16}(undef,kmax,5)
        @inbounds Aindex[:,1] = deepcopy(parHier.h_mu.clust[:,ian])
        @inbounds Aindex[:,2] = deepcopy(parHier.h_eta.clust[:,ian])
        @inbounds Aindex[:,3] = deepcopy(parHier.h_nu.clust[:,ian])
        @inbounds Aindex[:,4] = deepcopy(parHier.h_rho.clust[:,ian])
        @inbounds Aindex[:,5] = deepcopy(parHier.h_sigma.clust[:,ian])

        @inbounds Likelihood.eta.parameteracc[ksamp2]       =  deepcopy(Likelihood.eta.parameteracc[ksamp1])
        @inbounds Likelihood.mu.parameteracc[ksamp2]        =  deepcopy(Likelihood.mu.parameteracc[ksamp1])
        @inbounds Likelihood.rho.parameteracc[ksamp2]       =  deepcopy(Likelihood.rho.parameteracc[ksamp1])
        @inbounds Likelihood.nu.parameteracc[ksamp2]        =  deepcopy(Likelihood.nu.parameteracc[ksamp1])
        @inbounds Likelihood.sigma.parameteracc[ksamp2]     =  deepcopy(Likelihood.sigma.parameteracc[ksamp1])
        @inbounds Likelihood.sigma.parameteraccinv[ksamp2]  =  deepcopy(Likelihood.sigma.parameteraccinv[ksamp1])


        @inbounds Likelihood.eta.parameterprop[ksamp2]       =  deepcopy(Likelihood.eta.parameterprop[ksamp1])
        @inbounds Likelihood.mu.parameterprop[ksamp2]        =  deepcopy(Likelihood.mu.parameterprop[ksamp1])
        @inbounds Likelihood.rho.parameterprop[ksamp2]       =  deepcopy(Likelihood.rho.parameterprop[ksamp1])
        @inbounds Likelihood.nu.parameterprop[ksamp2]        =  deepcopy(Likelihood.nu.parameterprop[ksamp1])
        @inbounds Likelihood.sigma.parameterprop[ksamp2]     =  deepcopy(Likelihood.sigma.parameterprop[ksamp1])
        @inbounds Likelihood.sigma.parameterpropinv[ksamp2]  =  deepcopy(Likelihood.sigma.parameterpropinv[ksamp1])



        @inbounds Aindex[ksamp2,1] = Aindex[ksamp1,1]
        @inbounds Aindex[ksamp2,2] = Aindex[ksamp1,2]
        @inbounds Aindex[ksamp2,3] = Aindex[ksamp1,3]
        @inbounds Aindex[ksamp2,4] = Aindex[ksamp1,4]
        @inbounds Aindex[ksamp2,5] = Aindex[ksamp1,5]

        @inbounds parHier.h_mu.clust[ksamp2,ian] = parHier.h_mu.clust[ksamp1,ian]
        @inbounds parHier.h_eta.clust[ksamp2,ian] = parHier.h_eta.clust[ksamp1,ian]
        @inbounds parHier.h_nu.clust[ksamp2,ian] = parHier.h_nu.clust[ksamp1,ian]
        @inbounds parHier.h_rho.clust[ksamp2,ian] = parHier.h_rho.clust[ksamp1,ian]
        @inbounds parHier.h_sigma.clust[ksamp2,ian] = parHier.h_sigma.clust[ksamp1,ian]

        #ru = Int16(rand(rng,1:5))
        ru_1 = Int16(rand(rng,1:2))
        ru_2 = Int16(rand(rng,1:2))
        ru_3 = Int16(rand(rng,1:2))
        ru_4 = Int16(rand(rng,1:2))
        ru_5 = Int16(rand(rng,1:2))


        kbeta_1     = Int16(rand(rng,1:kmax,1)[1])
        kbeta_2     = Int16(rand(rng,1:kmax,1)[1])
        kbeta_3     = Int16(rand(rng,1:kmax,1)[1])
        kbeta_4     = Int16(rand(rng,1:kmax,1)[1])
        kbeta_5     = Int16(rand(rng,1:kmax,1)[1])

        if ru_1 == 1
            Aindex[ksamp2,1] = kbeta_1
        end
        if ru_2 == 1
            Aindex[ksamp2,2] = kbeta_2
        end
        if ru_3 == 1
            Aindex[ksamp2,3] = kbeta_3
        end
        if ru_4 == 1
            Aindex[ksamp2,4] = kbeta_4
        end
        if ru_5 == 1
            Aindex[ksamp2,5] = kbeta_5
        end



        iWhile        = Int16(0)
        while (check_unique(Aindex,ksamp2,kmax)!= true) && (iWhile<100)
            iWhile += one(iWhile)

            ru_1 = Int16(rand(rng,1:2))
            ru_2 = Int16(rand(rng,1:2))
            ru_3 = Int16(rand(rng,1:2))
            ru_4 = Int16(rand(rng,1:2))
            ru_5 = Int16(rand(rng,1:2))

            kbeta_1     = Int16(rand(rng,1:kmax,1)[1])
            kbeta_2     = Int16(rand(rng,1:kmax,1)[1])
            kbeta_3     = Int16(rand(rng,1:kmax,1)[1])
            kbeta_4     = Int16(rand(rng,1:kmax,1)[1])
            kbeta_5     = Int16(rand(rng,1:kmax,1)[1])

            if ru_1 == 1
                Aindex[ksamp2,1] = kbeta_1
            end
            if ru_2 == 1
                Aindex[ksamp2,2] = kbeta_2
            end
            if ru_3 == 1
                Aindex[ksamp2,3] = kbeta_3
            end
            if ru_4 == 1
                Aindex[ksamp2,4] = kbeta_4
            end
            if ru_5 == 1
                Aindex[ksamp2,5] = kbeta_5
            end
        end

        if iWhile >=100
            error("iWhile splitmerge_v2")
        end

        if ru_1 == 1
            parHier.h_mu.clust[ksamp2,ian] = kbeta_1
            @inbounds Likelihood.mu.parameteracc[ksamp2]    = deepcopy(parHier.h_mu.par.parameteracc[kbeta_1])
            @inbounds Likelihood.mu.parameterprop[ksamp2]    = deepcopy(parHier.h_mu.par.parameteracc[kbeta_1])
        end
        if ru_2 == 1
            parHier.h_eta.clust[ksamp2,ian] = kbeta_2
            @inbounds Likelihood.eta.parameteracc[ksamp2]    = deepcopy(parHier.h_eta.par.parameteracc[kbeta_2])
            @inbounds Likelihood.eta.parameterprop[ksamp2]    = deepcopy(parHier.h_eta.par.parameteracc[kbeta_2])
        end
        if ru_3 == 1
            parHier.h_nu.clust[ksamp2,ian] = kbeta_3
            @inbounds Likelihood.nu.parameteracc[ksamp2]    = deepcopy(parHier.h_nu.par.parameteracc[kbeta_3])
            @inbounds Likelihood.nu.parameterprop[ksamp2]    = deepcopy(parHier.h_nu.par.parameteracc[kbeta_3])
        end
        if ru_4 == 1
            parHier.h_rho.clust[ksamp2,ian] = kbeta_4
            @inbounds Likelihood.rho.parameteracc[ksamp2]    = deepcopy(parHier.h_rho.par.parameteracc[kbeta_4])
            @inbounds Likelihood.rho.parameterprop[ksamp2]    = deepcopy(parHier.h_rho.par.parameteracc[kbeta_4])
        end
        if ru_5 == 1
            parHier.h_sigma.clust[ksamp2,ian] = kbeta_5
            @inbounds Likelihood.sigma.parameteracc[ksamp2]    = deepcopy(parHier.h_sigma.par.parameteracc[kbeta_5])
            @inbounds Likelihood.sigma.parameterprop[ksamp2]    = deepcopy(parHier.h_sigma.par.parameteracc[kbeta_5])

            @inbounds Likelihood.sigma.parameteraccinv[ksamp2]    = deepcopy(parHier.h_sigma.par.parameteraccinv[kbeta_5])
            @inbounds Likelihood.sigma.parameterpropinv[ksamp2]    = deepcopy(parHier.h_sigma.par.parameteraccinv[kbeta_5])

        end

        app = Float64(0.0)
        j1 = parHier.h_mu.clust[ksamp2,ian]
        app += log(parHier.h_mu.prob[j1])
        j1 = parHier.h_eta.clust[ksamp2,ian]
        app += log(parHier.h_eta.prob[j1])
        j1 = parHier.h_nu.clust[ksamp2,ian]
        app += log(parHier.h_nu.prob[j1])
        j1 = parHier.h_rho.clust[ksamp2,ian]
        app += log(parHier.h_rho.prob[j1])
        j1 = parHier.h_sigma.clust[ksamp2,ian]
        app += log(parHier.h_sigma.prob[j1])

        Clusterization.mcmc_beta[ksamp2] = max(1.0e-100,exp(app))

        ### END NEW



        t_index_merge = Int16.(findall(x->x == ksamp1, Likelihood.clusterization.zeta[1:(nt-1)]))
        sets_vec      = find_const(Likelihood.clusterization.zeta[1:(nt-1)],ksamp1)

        t_index_split1 = zeros(Int16,length(t_index_merge))
        t_index_split2 = zeros(Int16,length(t_index_merge))

        nsplit1        = Int16(0)
        nsplit2        = Int16(0)

        ## n i to j
        nitoj_merge     = deepcopy(Clust.n_itojC[ksamp1])

        nitoj_split1    = zeros(Int16, kmax)
        nitoj_split2    = zeros(Int16, kmax)

        from_nitoj_split1    = zeros(Int16, kmax)
        from_nitoj_split2    = zeros(Int16, kmax)

        zprec = Int16(1)
        for ii in 1:size(sets_vec,1)

            if sets_vec[ii,1] != 1
                zprec = sets_vec[ii,1]-Int16(1)
            end

            u = Int16(rand(rng,1:2,1)[1])
            if u == 1

                nitoj_split1[ksamp1] += sets_vec[ii,2]-sets_vec[ii,1]
                if sets_vec[ii,2]< (nt-1)
                    nitoj_split1[Likelihood.clusterization.zeta[sets_vec[ii,2]+1]] += Int16(1)
                end

                for jj in sets_vec[ii,1]:sets_vec[ii,2]
                    nsplit1 += one(nsplit1)
                    t_index_split1[nsplit1] = jj

                end

                from_nitoj_split1[Likelihood.clusterization.zeta[zprec]] += Int16(1)

            else
                nitoj_split2[ksamp2] += sets_vec[ii,2]-sets_vec[ii,1]
                if sets_vec[ii,2]< (nt-1)
                    nitoj_split2[Likelihood.clusterization.zeta[sets_vec[ii,2]+1]] += Int16(1)
                end

                for jj in sets_vec[ii,1]:sets_vec[ii,2]
                    nsplit2 += one(nsplit2)
                    t_index_split2[nsplit2] = jj
                end

                from_nitoj_split2[Likelihood.clusterization.zeta[zprec]] += Int16(1)
            end

            #print("\n\n",sum(sets_vec),"\n",sets_vec,"\n\n")
        end
        #print("\n n", nsplit1," ",nsplit2," ",nsplit1+nsplit2," \n")
        t_index_split1 = t_index_split1[1:nsplit1]
        t_index_split2 = t_index_split2[1:nsplit2]
        # proposal
        proposal_merge += -log(Float64(Likelihood.clusterization.n_nonemptyC[1]+1))-log(Float64(Likelihood.clusterization.n_nonemptyC[1]))

        proposal_split += -log(Float64(Likelihood.clusterization.n_nonemptyC[1]))
        proposal_split += -Float64(Likelihood.clusterization.n_nonemptyC[1])*log(2.0)

        sumSim = 0
        sumSim += parHier.h_mu.clust[ksamp2,ian] == parHier.h_mu.clust[ksamp1,ian]
        sumSim += parHier.h_eta.clust[ksamp2,ian] == parHier.h_eta.clust[ksamp1,ian]
        sumSim += parHier.h_nu.clust[ksamp2,ian] == parHier.h_nu.clust[ksamp1,ian]
        sumSim += parHier.h_rho.clust[ksamp2,ian] == parHier.h_rho.clust[ksamp1,ian]
        sumSim += parHier.h_sigma.clust[ksamp2,ian] == parHier.h_sigma.clust[ksamp1,ian]


        proposal_split += -5.0*log(2.0)-Float64(sumSim)*log(Float64(kmax))




    end
    #print("\n ", domerge," ", sum(nitoj_merge)==(sum(nitoj_split1)+sum(nitoj_split2))," ", sum(nitoj_merge)," ",sum(nitoj_split1)," ",sum(nitoj_split2), "\n")
    prob_split = 0.0;
    prob_merge = 0.0;
    # for t in t_index_split1
    #     @inbounds prob_split += compute_logpdf(Likelihood, ksamp1, Int16(t+1))
    # end
    # prob_merge = deepcopy(prob_split)

    for t in t_index_split2
        @inbounds prob_split += compute_logpdf(Likelihood, ksamp2, Int16(t+1))
        @inbounds prob_merge += compute_logpdf(Likelihood, ksamp1, Int16(t+1))
    end

    #print("\n pre ",domerge," ", prob_split," ",prob_merge,"\n")


    for kk in 1:kmax

        if kk!=ksamp1
            prob_split += loggamma(alphaPar*Clusterization.mcmc_beta[kk]+ nitoj_split1[kk])
            prob_merge += loggamma(alphaPar*Clusterization.mcmc_beta[kk]+ nitoj_merge[kk])
        else
            prob_split += loggamma(alphaPar*Clusterization.mcmc_beta[kk]+ nitoj_split1[kk]+ kappaPar)
            prob_merge += loggamma(alphaPar*Clusterization.mcmc_beta[kk]+ nitoj_merge[kk]+ kappaPar)

        end
        if kk!=ksamp2
            prob_split += loggamma(alphaPar*Clusterization.mcmc_beta[kk]+ nitoj_split2[kk])
            prob_merge += loggamma(alphaPar*Clusterization.mcmc_beta[kk])
        else
            prob_split += loggamma(alphaPar*Clusterization.mcmc_beta[kk]+ nitoj_split2[kk]+ kappaPar)
            prob_merge += loggamma(alphaPar*Clusterization.mcmc_beta[kk]+ kappaPar)
        end


    end

    prob_split -= loggamma(alphaPar+ sum(nitoj_split1)+ kappaPar)
    prob_split -= loggamma(alphaPar+ sum(nitoj_split2)+ kappaPar)
    prob_merge -= loggamma(alphaPar+ sum(nitoj_merge)+ kappaPar)
    prob_merge -= loggamma(alphaPar+  kappaPar)

    if domerge == 1
        for kk in 1:kmax

            if (kk!=ksamp1) && (kk!=ksamp2)
                prob_split += loggamma(alphaPar*Clusterization.mcmc_beta[ksamp1]+ Clust.n_itojC[kk][ksamp1])
                prob_split += loggamma(alphaPar*Clusterization.mcmc_beta[ksamp2]+ Clust.n_itojC[kk][ksamp2])

                prob_merge += loggamma(alphaPar*Clusterization.mcmc_beta[ksamp1]+ Clust.n_itojC[kk][ksamp2]+ Clust.n_itojC[kk][ksamp1])
                prob_merge += loggamma(alphaPar*Clusterization.mcmc_beta[ksamp2])
            end


        end
    else
        for kk in 1:kmax

            if (kk!=ksamp1) && (kk!=ksamp2)
                prob_split += loggamma(alphaPar*Clusterization.mcmc_beta[ksamp1]+ from_nitoj_split1[kk])
                prob_split += loggamma(alphaPar*Clusterization.mcmc_beta[ksamp2]+ from_nitoj_split2[kk])

                prob_merge += loggamma(alphaPar*Clusterization.mcmc_beta[ksamp1]+ Clust.n_itojC[kk][ksamp1])
                prob_merge += loggamma(alphaPar*Clusterization.mcmc_beta[ksamp2])

                #print("\n DD ", " ",kk," ",Clust.n_itojC[kk][ksamp1]== from_nitoj_split1[kk]+from_nitoj_split2[kk], Clust.n_itojC[kk][ksamp1], " ", from_nitoj_split1[kk]," ",from_nitoj_split2[kk])
            end

        end
    end


    if domerge == 1
        #merge
        MH = (prob_merge+proposal_split)-(prob_split+proposal_merge)

        #print(MH," P ",prob_split," ",proposal_merge," ",prob_merge," " ,proposal_split, "\n")
        u = rand(rng,Uniform( 0.0,1.0))

        if u<exp(MH)
            for t in t_index_merge
                @inbounds Likelihood.clusterization.zeta[t] = ksamp1
            end
            Update_ObsInClust!(Likelihood.clusterization)
        end

    else
        #split
        MH = (prob_split+proposal_merge)-(prob_merge+proposal_split)
        #print("\n","Split",ksamp1," ",ksamp2,"\n")
        #print(MH," P ",prob_split," ",proposal_merge, " ",prob_merge," " ,proposal_split, "\n")
        u = rand(rng,Uniform( 0.0,1.0))

        if u<exp(MH)
            for t in t_index_split1
                @inbounds Likelihood.clusterization.zeta[t] = ksamp1
            end
            for t in t_index_split2
                @inbounds Likelihood.clusterization.zeta[t] = ksamp2
            end
            Update_ObsInClust!(Likelihood.clusterization)
        end
    end
    #print(exp(MH),"\n")
end





####SECT HMM Joint
function sample_zeta!(rng::MersenneTwister,Likelihood::Vector{Likelihood_OU_CircLinmodel}, Clusterization::Vector{Clusterization_HDPHMM_Joint},doupdate::ZetaDoNotUpdate)



end

function sample_zeta!(rng::MersenneTwister,Likelihood::Vector{Likelihood_OU_CircLinmodel}, Clusterization::Vector{Clusterization_HDPHMM_Joint},doupdate::ZetaDoUpdate)

    u = rand(rng,Uniform(0.0,1.0))
    if u<0.95
        sample_zeta_beam!(rng,Likelihood[1],Clusterization[1],Likelihood)
    else
        sample_zeta_marg!(rng,Likelihood[1],Clusterization[1],Likelihood)
    end

end
function sample_zeta_beam!(rng::MersenneTwister,Likelihood::Likelihood_OU_CircLinmodel, Clusterization::Clusterization_HDPHMM_Joint,Likelihood_JOINT::Vector{Likelihood_OU_CircLinmodel})

    Moltbeam = Float64(1.0)
    # Likelihood      = MCMCLikelihood
    # Clusterization  = MCMCClusterization

    nanim_v2 = size(Likelihood_JOINT)[1]
#    Clusterization.pi.parameteracc
    kmax        = Likelihood.kmax
    nc          = Likelihood.data.ncol
    nanim       = Likelihood.data.nanimals
    nt          = Likelihood.data.nt

    zetaprex = 1
    veckacc = zeros(Int16,kmax)
    probsVect       = zeros(Float64,kmax)



    un    = zeros(Float64,nt)
    zetaprex = 1
    for i  in 1:(nt-1)
        @inbounds k1          = Likelihood.clusterization.zeta[i]
        @inbounds un[i]       = rand(rng,Uniform(0.0, Clusterization.pi.parameteracc[zetaprex][k1]/Moltbeam))
        @inbounds zetaprex    = k1
    end

    zetaprex = 1
    @inbounds for i  in 1:(nt-2)
        nk              = 0
         @inbounds zetasux        = Likelihood.clusterization.zeta[i+1]
         @inbounds kold           = Likelihood.clusterization.zeta[i]
        for k in 1:kmax

            if (Clusterization.pi.parameteracc[zetaprex][k]>un[i]) && (Clusterization.pi.parameteracc[k][zetasux]>un[i+1])
                nk            += one(nk)
                @inbounds veckacc[nk]   = k
                #probsVect[nk] = LogPdfMatrix[i+1,k] #comp
                @inbounds probsVect[nk] = compute_logpdf(Likelihood_JOINT, Int16(k), Int16(i+1)) #comp
                #probsVect[nk] += log( max(Clusterization.pi.parameteracc[zetaprex][k]*1000,1.0e-300))
                #probsVect[nk] += log( max(Clusterization.pi.parameteracc[k][zetasux]*1000,1.0e-300))
            end
        end

        @inbounds sampvec = probsVect[1:nk]
        @inbounds Likelihood.clusterization.zeta[i] = veckacc[sample_discretevar(rng,sampvec)]
        @inbounds zetaprex = Likelihood.clusterization.zeta[i]

        if nanim_v2>1
            for ian in 2:nanim_v2
                @inbounds Likelihood_JOINT[ian].clusterization.zeta[i] = Likelihood.clusterization.zeta[i]
            end
        end
    end
    i       = nt-1
    nk      = 0
    kold           = Likelihood.clusterization.zeta[i]
    for k in 1:kmax

        if Clusterization.pi.parameteracc[zetaprex][k]>un[i]
            nk           += one(nk)
            veckacc[nk]   = k
            @inbounds probsVect[nk] = compute_logpdf(Likelihood_JOINT, Int16(k), Int16(i+1)) # # compute_logpdf(Likelihood, Int16(k), Int16(i+1))
            @inbounds probsVect[nk] += log( max( Clusterization.pi.parameteracc[zetaprex][k]*1000,1.0e-300))
        end

    end
    sampvec = probsVect[1:nk]
    @inbounds Likelihood.clusterization.zeta[i] = veckacc[sample_discretevar(rng,sampvec)]

    Update_ObsInClust!(Likelihood.clusterization)

    nanim_v2 = size(Likelihood_JOINT)[1]
    #print(nanim_v2)
    if nanim_v2>1
        for ian in 2:nanim_v2
            for i in 1:(nt)
                Likelihood_JOINT[ian].clusterization.zeta[i] = Likelihood.clusterization.zeta[i]
            end

            Update_ObsInClust!(Likelihood_JOINT[ian].clusterization)
        end
    end

    return nothing

end

##### MARGINALIZATION

function sample_zeta_marg!(rng::MersenneTwister,Likelihood::Likelihood_OU_CircLinmodel, Clusterization::Clusterization_HDPHMM_Joint,Likelihood_JOINT::Vector{Likelihood_OU_CircLinmodel})

    #Likelihood      = MCMCLikelihood
    # Clusterization  = MCMCClusterization
    nanim_v2 = size(Likelihood_JOINT)[1]
    Clust = Clusterization.clusterization
#    Clusterization.pi.parameteracc
    kmax        = Likelihood.kmax
    nc          = Likelihood.data.ncol
    nanim       = Likelihood.data.nanimals
    nt          = Likelihood.data.nt

    zetaprex = 1
    veckacc = zeros(Int16,kmax)
    probsVect       = zeros(Float64,kmax)



    # LogPdfMatrix = zeros(Float64,nt,kmax)
    # for i  in 1:(nt-1)
    #     for k in 1:kmax
    #         @inbounds LogPdfMatrix[i+1,k] = compute_logpdf(Likelihood, Int16(k), Int16(i+1))
    #     end
    # end

    #
    # un    = zeros(Float64,nt)
    # zetaprex = 1
    # for i  in 1:(nt-1)
    #     @inbounds k1          = Likelihood.clusterization.zeta[i]
    #     @inbounds un[i]       = rand(Uniform(0.0, Clusterization.pi.parameteracc[zetaprex][k1]))
    #     @inbounds zetaprex    = k1
    # end
    kappaPar = Clusterization.mcmc_ak[1]*Clusterization.mcmc_rho[1]
    alphaPar = Clusterization.mcmc_ak[1]-kappaPar
    zetaprex = 1
    @inbounds for i  in 1:(nt-2)

        zetasux        = Likelihood.clusterization.zeta[i+1]
        kold           = Likelihood.clusterization.zeta[i]
        Clust.n_itojC[zetaprex][kold] -= 1
        Clust.n_itojC[kold][zetasux]  -= 1
        for k in 1:kmax

            probsVect[k] = compute_logpdf(Likelihood_JOINT, Int16(k), Int16(i+1))

            if k == zetaprex
                probsVect[k] += log(alphaPar*Clusterization.mcmc_beta[k]+ Clust.n_itojC[zetaprex][k] + kappaPar)
            else
                probsVect[k] += log(alphaPar*Clusterization.mcmc_beta[k]+ Clust.n_itojC[zetaprex][k])
            end

            if k == zetasux
                probsVect[k] += log(alphaPar*Clusterization.mcmc_beta[zetasux]+ Clust.n_itojC[k][zetasux] + kappaPar)
            else
                probsVect[k] += log(alphaPar*Clusterization.mcmc_beta[zetasux]+ Clust.n_itojC[k][zetasux])
            end
        end

        Likelihood.clusterization.zeta[i] = sample_discretevar(rng,probsVect)
        if nanim_v2>1
            for ian in 2:nanim_v2
                @inbounds Likelihood_JOINT[ian].clusterization.zeta[i] = Likelihood.clusterization.zeta[i]
            end
        end
        Clust.n_itojC[zetaprex][Likelihood.clusterization.zeta[i]] += 1
        Clust.n_itojC[Likelihood.clusterization.zeta[i]][zetasux]  += 1

        zetaprex = Likelihood.clusterization.zeta[i]

    end
    i       = nt-1
    nk      = 0
    kold           = Likelihood.clusterization.zeta[i]
    Clust.n_itojC[zetaprex][kold] -= 1
    for k in 1:kmax

        probsVect[k] = compute_logpdf(Likelihood_JOINT, Int16(k), Int16(i+1))

        if k == zetaprex
            probsVect[k] += log(alphaPar*Clusterization.mcmc_beta[k]+ Clust.n_itojC[zetaprex][k] + kappaPar)
        else
            probsVect[k] += log(alphaPar*Clusterization.mcmc_beta[k]+ Clust.n_itojC[zetaprex][k])
        end


    end

    Likelihood.clusterization.zeta[i] = sample_discretevar(rng,probsVect)

    Update_ObsInClust!(Likelihood.clusterization)


    if nanim_v2>1
        for ian in 2:nanim_v2
            for i in 1:(nt)
                Likelihood_JOINT[ian].clusterization.zeta[i] = Likelihood.clusterization.zeta[i]
            end

            Update_ObsInClust!(Likelihood_JOINT[ian].clusterization)
        end
    end


    return nothing

end








#### ALTRO METODO
####SECT HMM Joint
# function sample_zeta!(rng::MersenneTwister,Likelihood::Vector{Likelihood_OU_CircLinmodel}, Clusterization::Vector{Clusterization_HDPHMM_Joint},doupdate::ZetaDoNotUpdate)
#
#
#
# end

# function sample_zeta!(rng::MersenneTwister,Likelihood::Vector{Likelihood_OU_CircLinmodel}, Clusterization::Vector{Clusterization_HDPHMM_Joint},doupdate::ZetaDoUpdate)
#
#     u = rand(rng,Uniform(0.0,1.0))
#     if u<0.95
#         sample_zeta_beam!(rng,Likelihood[1],Clusterization[1],Likelihood)
#     else
#         sample_zeta_marg!(rng,Likelihood[1],Clusterization[1],Likelihood)
#     end
#
# end
# function sample_zeta_beam!(rng::MersenneTwister,Likelihood::Likelihood_OU_CircLinmodel, Clusterization::Clusterization_HDPHMM_Joint,Likelihood_JOINT::Vector{Likelihood_OU_CircLinmodel})
#
#     Moltbeam = Float64(1.0)
#     # Likelihood      = MCMCLikelihood
#     # Clusterization  = MCMCClusterization
#
#     nanim_v2 = size(Likelihood_JOINT)[1]
# #    Clusterization.pi.parameteracc
#     kmax        = Likelihood.kmax
#     nc          = Likelihood.data.ncol
#     nanim       = Likelihood.data.nanimals
#     nt          = Likelihood.data.nt
#
#     zetaprex = 1
#     veckacc = zeros(Int16,kmax)
#     probsVect       = zeros(Float64,kmax)
#
#
#
#     un    = zeros(Float64,nt)
#     zetaprex = 1
#     for i  in 1:(nt-1)
#         @inbounds k1          = Likelihood.clusterization.zeta[i]
#         @inbounds un[i]       = rand(rng,Uniform(0.0, Clusterization.pi.parameteracc[zetaprex][k1]/Moltbeam))
#         @inbounds zetaprex    = k1
#     end
#
#     zetaprex = 1
#     @inbounds for i  in 1:(nt-2)
#         nk              = 0
#          @inbounds zetasux        = Likelihood.clusterization.zeta[i+1]
#          @inbounds kold           = Likelihood.clusterization.zeta[i]
#         for k in 1:kmax
#
#             if (Clusterization.pi.parameteracc[zetaprex][k]>un[i]) && (Clusterization.pi.parameteracc[k][zetasux]>un[i+1])
#                 nk            += one(nk)
#                 @inbounds veckacc[nk]   = k
#                 #probsVect[nk] = LogPdfMatrix[i+1,k] #comp
#                 @inbounds probsVect[nk] = compute_logpdf(Likelihood_JOINT, Int16(k), Int16(i+1)) #comp
#                 #probsVect[nk] += log( max(Clusterization.pi.parameteracc[zetaprex][k]*1000,1.0e-300))
#                 #probsVect[nk] += log( max(Clusterization.pi.parameteracc[k][zetasux]*1000,1.0e-300))
#             end
#         end
#
#         @inbounds sampvec = probsVect[1:nk]
#         @inbounds Likelihood.clusterization.zeta[i] = veckacc[sample_discretevar(rng,sampvec)]
#         @inbounds zetaprex = Likelihood.clusterization.zeta[i]
#
#         if nanim_v2>1
#             for ian in 2:nanim_v2
#                 @inbounds Likelihood_JOINT[ian].clusterization.zeta[i] = Likelihood.clusterization.zeta[i]
#             end
#         end
#     end
#     i       = nt-1
#     nk      = 0
#     kold           = Likelihood.clusterization.zeta[i]
#     for k in 1:kmax
#
#         if Clusterization.pi.parameteracc[zetaprex][k]>un[i]
#             nk           += one(nk)
#             veckacc[nk]   = k
#             @inbounds probsVect[nk] = compute_logpdf(Likelihood_JOINT, Int16(k), Int16(i+1)) # # compute_logpdf(Likelihood, Int16(k), Int16(i+1))
#             @inbounds probsVect[nk] += log( max( Clusterization.pi.parameteracc[zetaprex][k]*1000,1.0e-300))
#         end
#
#     end
#     sampvec = probsVect[1:nk]
#     @inbounds Likelihood.clusterization.zeta[i] = veckacc[sample_discretevar(rng,sampvec)]
#
#     Update_ObsInClust!(Likelihood.clusterization)
#
#     nanim_v2 = size(Likelihood_JOINT)[1]
#     #print(nanim_v2)
#     if nanim_v2>1
#         for ian in 2:nanim_v2
#             for i in 1:(nt)
#                 Likelihood_JOINT[ian].clusterization.zeta[i] = Likelihood.clusterization.zeta[i]
#             end
#
#             Update_ObsInClust!(Likelihood_JOINT[ian].clusterization)
#         end
#     end
#
#     return nothing
#
# end
#
# ##### MARGINALIZATION
#
# function sample_zeta_marg!(rng::MersenneTwister,Likelihood::Likelihood_OU_CircLinmodel, Clusterization::Clusterization_HDPHMM_Joint,Likelihood_JOINT::Vector{Likelihood_OU_CircLinmodel})
#
#     #Likelihood      = MCMCLikelihood
#     # Clusterization  = MCMCClusterization
#     nanim_v2 = size(Likelihood_JOINT)[1]
#     Clust = Clusterization.clusterization
# #    Clusterization.pi.parameteracc
#     kmax        = Likelihood.kmax
#     nc          = Likelihood.data.ncol
#     nanim       = Likelihood.data.nanimals
#     nt          = Likelihood.data.nt
#
#     zetaprex = 1
#     veckacc = zeros(Int16,kmax)
#     probsVect       = zeros(Float64,kmax)
#
#
#
#     # LogPdfMatrix = zeros(Float64,nt,kmax)
#     # for i  in 1:(nt-1)
#     #     for k in 1:kmax
#     #         @inbounds LogPdfMatrix[i+1,k] = compute_logpdf(Likelihood, Int16(k), Int16(i+1))
#     #     end
#     # end
#
#     #
#     # un    = zeros(Float64,nt)
#     # zetaprex = 1
#     # for i  in 1:(nt-1)
#     #     @inbounds k1          = Likelihood.clusterization.zeta[i]
#     #     @inbounds un[i]       = rand(Uniform(0.0, Clusterization.pi.parameteracc[zetaprex][k1]))
#     #     @inbounds zetaprex    = k1
#     # end
#     kappaPar = Clusterization.mcmc_ak[1]*Clusterization.mcmc_rho[1]
#     alphaPar = Clusterization.mcmc_ak[1]-kappaPar
#     zetaprex = 1
#     @inbounds for i  in 1:(nt-2)
#
#         zetasux        = Likelihood.clusterization.zeta[i+1]
#         kold           = Likelihood.clusterization.zeta[i]
#         Clust.n_itojC[zetaprex][kold] -= 1
#         Clust.n_itojC[kold][zetasux]  -= 1
#         for k in 1:kmax
#
#             probsVect[k] = compute_logpdf(Likelihood_JOINT, Int16(k), Int16(i+1))
#
#             if k == zetaprex
#                 probsVect[k] += log(alphaPar*Clusterization.mcmc_beta[k]+ Clust.n_itojC[zetaprex][k] + kappaPar)
#             else
#                 probsVect[k] += log(alphaPar*Clusterization.mcmc_beta[k]+ Clust.n_itojC[zetaprex][k])
#             end
#
#             if k == zetasux
#                 probsVect[k] += log(alphaPar*Clusterization.mcmc_beta[zetasux]+ Clust.n_itojC[k][zetasux] + kappaPar)
#             else
#                 probsVect[k] += log(alphaPar*Clusterization.mcmc_beta[zetasux]+ Clust.n_itojC[k][zetasux])
#             end
#         end
#
#         Likelihood.clusterization.zeta[i] = sample_discretevar(rng,probsVect)
#         if nanim_v2>1
#             for ian in 2:nanim_v2
#                 @inbounds Likelihood_JOINT[ian].clusterization.zeta[i] = Likelihood.clusterization.zeta[i]
#             end
#         end
#         Clust.n_itojC[zetaprex][Likelihood.clusterization.zeta[i]] += 1
#         Clust.n_itojC[Likelihood.clusterization.zeta[i]][zetasux]  += 1
#
#         zetaprex = Likelihood.clusterization.zeta[i]
#
#     end
#     i       = nt-1
#     nk      = 0
#     kold           = Likelihood.clusterization.zeta[i]
#     Clust.n_itojC[zetaprex][kold] -= 1
#     for k in 1:kmax
#
#         probsVect[k] = compute_logpdf(Likelihood_JOINT, Int16(k), Int16(i+1))
#
#         if k == zetaprex
#             probsVect[k] += log(alphaPar*Clusterization.mcmc_beta[k]+ Clust.n_itojC[zetaprex][k] + kappaPar)
#         else
#             probsVect[k] += log(alphaPar*Clusterization.mcmc_beta[k]+ Clust.n_itojC[zetaprex][k])
#         end
#
#
#     end
#
#     Likelihood.clusterization.zeta[i] = sample_discretevar(rng,probsVect)
#
#     Update_ObsInClust!(Likelihood.clusterization)
#
#
#     if nanim_v2>1
#         for ian in 2:nanim_v2
#             for i in 1:(nt)
#                 Likelihood_JOINT[ian].clusterization.zeta[i] = Likelihood.clusterization.zeta[i]
#             end
#
#             Update_ObsInClust!(Likelihood_JOINT[ian].clusterization)
#         end
#     end
#
#
#     return nothing
#
# end
#




####################################
##### PROPOSED MODELS - SAME Z
####################################


function sample_zeta_type2!(rng::MersenneTwister,Likelihood::Vector{Likelihood_OU_CircLinmodel}, Clusterization::Vector{Clusterization_HDPHMM},doupdate::ZetaDoNotUpdate)



end

function sample_zeta_type2!(rng::MersenneTwister,Likelihood::Vector{Likelihood_OU_CircLinmodel}, Clusterization::Vector{Clusterization_HDPHMM},doupdate::ZetaDoUpdate)

    u = rand(rng,Uniform(0.0,1.0))
    if u<0.95
        sample_zeta_beam_type2!(rng,Likelihood[1],Clusterization[1],Likelihood)
    else
        sample_zeta_marg_type2!(rng,Likelihood[1],Clusterization[1],Likelihood)
    end

end
function sample_zeta_beam_type2!(rng::MersenneTwister,Likelihood::Likelihood_OU_CircLinmodel, Clusterization::Clusterization_HDPHMM,Likelihood_JOINT::Vector{Likelihood_OU_CircLinmodel})

    Moltbeam = Float64(1.0)
    # Likelihood      = MCMCLikelihood
    # Clusterization  = MCMCClusterization

    nanim_v2 = size(Likelihood_JOINT)[1]
#    Clusterization.pi.parameteracc
    kmax        = Likelihood.kmax
    nc          = Likelihood.data.ncol
    nanim       = Likelihood.data.nanimals
    nt          = Likelihood.data.nt

    zetaprex = 1
    veckacc = zeros(Int16,kmax)
    probsVect       = zeros(Float64,kmax)



    un    = zeros(Float64,nt)
    zetaprex = 1
    for i  in 1:(nt-1)
        @inbounds k1          = Likelihood.clusterization.zeta[i]
        @inbounds un[i]       = rand(rng,Uniform(0.0, Clusterization.pi.parameteracc[zetaprex][k1]/Moltbeam))
        @inbounds zetaprex    = k1
    end

    zetaprex = 1
    @inbounds for i  in 1:(nt-2)
        nk              = 0
         @inbounds zetasux        = Likelihood.clusterization.zeta[i+1]
         @inbounds kold           = Likelihood.clusterization.zeta[i]
        for k in 1:kmax

            if (Clusterization.pi.parameteracc[zetaprex][k]>un[i]) && (Clusterization.pi.parameteracc[k][zetasux]>un[i+1])
                nk            += one(nk)
                @inbounds veckacc[nk]   = k
                #probsVect[nk] = LogPdfMatrix[i+1,k] #comp
                @inbounds probsVect[nk] = compute_logpdf(Likelihood_JOINT, Int16(k), Int16(i+1)) #comp
                #probsVect[nk] += log( max(Clusterization.pi.parameteracc[zetaprex][k]*1000,1.0e-300))
                #probsVect[nk] += log( max(Clusterization.pi.parameteracc[k][zetasux]*1000,1.0e-300))
            end
        end

        @inbounds sampvec = probsVect[1:nk]
        @inbounds Likelihood.clusterization.zeta[i] = veckacc[sample_discretevar(rng,sampvec)]
        @inbounds zetaprex = Likelihood.clusterization.zeta[i]

        if nanim_v2>1
            for ian in 2:nanim_v2
                @inbounds Likelihood_JOINT[ian].clusterization.zeta[i] = Likelihood.clusterization.zeta[i]
            end
        end
    end
    i       = nt-1
    nk      = 0
    kold           = Likelihood.clusterization.zeta[i]
    for k in 1:kmax

        if Clusterization.pi.parameteracc[zetaprex][k]>un[i]
            nk           += one(nk)
            veckacc[nk]   = k
            @inbounds probsVect[nk] = compute_logpdf(Likelihood_JOINT, Int16(k), Int16(i+1)) # # compute_logpdf(Likelihood, Int16(k), Int16(i+1))
            @inbounds probsVect[nk] += log( max( Clusterization.pi.parameteracc[zetaprex][k]*1000,1.0e-300))
        end

    end
    sampvec = probsVect[1:nk]
    @inbounds Likelihood.clusterization.zeta[i] = veckacc[sample_discretevar(rng,sampvec)]

    Update_ObsInClust!(Likelihood.clusterization)

    nanim_v2 = size(Likelihood_JOINT)[1]
    #print(nanim_v2)
    if nanim_v2>1
        for ian in 2:nanim_v2
            for i  in 1:(nt)
                Likelihood_JOINT[ian].clusterization.zeta[i] = Likelihood.clusterization.zeta[i]
            end
            Update_ObsInClust!(Likelihood_JOINT[ian].clusterization)
        end
    end

    return nothing

end

##### MARGINALIZATION

function sample_zeta_marg_type2!(rng::MersenneTwister,Likelihood::Likelihood_OU_CircLinmodel, Clusterization::Clusterization_HDPHMM,Likelihood_JOINT::Vector{Likelihood_OU_CircLinmodel})

    #Likelihood      = MCMCLikelihood
    # Clusterization  = MCMCClusterization
    nanim_v2 = size(Likelihood_JOINT)[1]
    Clust = Clusterization.clusterization
#    Clusterization.pi.parameteracc
    kmax        = Likelihood.kmax
    nc          = Likelihood.data.ncol
    nanim       = Likelihood.data.nanimals
    nt          = Likelihood.data.nt

    zetaprex = 1
    veckacc = zeros(Int16,kmax)
    probsVect       = zeros(Float64,kmax)



    # LogPdfMatrix = zeros(Float64,nt,kmax)
    # for i  in 1:(nt-1)
    #     for k in 1:kmax
    #         @inbounds LogPdfMatrix[i+1,k] = compute_logpdf(Likelihood, Int16(k), Int16(i+1))
    #     end
    # end

    #
    # un    = zeros(Float64,nt)
    # zetaprex = 1
    # for i  in 1:(nt-1)
    #     @inbounds k1          = Likelihood.clusterization.zeta[i]
    #     @inbounds un[i]       = rand(Uniform(0.0, Clusterization.pi.parameteracc[zetaprex][k1]))
    #     @inbounds zetaprex    = k1
    # end
    kappaPar = Clusterization.mcmc_ak[1]*Clusterization.mcmc_rho[1]
    alphaPar = Clusterization.mcmc_ak[1]-kappaPar
    zetaprex = 1
    @inbounds for i  in 1:(nt-2)

        zetasux        = Likelihood.clusterization.zeta[i+1]
        kold           = Likelihood.clusterization.zeta[i]
        Clust.n_itojC[zetaprex][kold] -= 1
        Clust.n_itojC[kold][zetasux]  -= 1
        for k in 1:kmax

            probsVect[k] = compute_logpdf(Likelihood_JOINT, Int16(k), Int16(i+1))

            if k == zetaprex
                probsVect[k] += log(alphaPar*Clusterization.mcmc_beta[k]+ Clust.n_itojC[zetaprex][k] + kappaPar)
            else
                probsVect[k] += log(alphaPar*Clusterization.mcmc_beta[k]+ Clust.n_itojC[zetaprex][k])
            end

            if k == zetasux
                probsVect[k] += log(alphaPar*Clusterization.mcmc_beta[zetasux]+ Clust.n_itojC[k][zetasux] + kappaPar)
            else
                probsVect[k] += log(alphaPar*Clusterization.mcmc_beta[zetasux]+ Clust.n_itojC[k][zetasux])
            end
        end

        Likelihood.clusterization.zeta[i] = sample_discretevar(rng,probsVect)
        if nanim_v2>1
            for ian in 2:nanim_v2
                @inbounds Likelihood_JOINT[ian].clusterization.zeta[i] = Likelihood.clusterization.zeta[i]
            end
        end
        Clust.n_itojC[zetaprex][Likelihood.clusterization.zeta[i]] += 1
        Clust.n_itojC[Likelihood.clusterization.zeta[i]][zetasux]  += 1

        zetaprex = Likelihood.clusterization.zeta[i]

    end
    i       = nt-1
    nk      = 0
    kold           = Likelihood.clusterization.zeta[i]
    Clust.n_itojC[zetaprex][kold] -= 1
    for k in 1:kmax

        probsVect[k] = compute_logpdf(Likelihood_JOINT, Int16(k), Int16(i+1))

        if k == zetaprex
            probsVect[k] += log(alphaPar*Clusterization.mcmc_beta[k]+ Clust.n_itojC[zetaprex][k] + kappaPar)
        else
            probsVect[k] += log(alphaPar*Clusterization.mcmc_beta[k]+ Clust.n_itojC[zetaprex][k])
        end


    end

    Likelihood.clusterization.zeta[i] = sample_discretevar(rng,probsVect)

    Update_ObsInClust!(Likelihood.clusterization)


    if nanim_v2>1
        for ian in 2:nanim_v2
            for i  in 1:(nt)
                Likelihood_JOINT[ian].clusterization.zeta[i] = Likelihood.clusterization.zeta[i]
            end
            Update_ObsInClust!(Likelihood_JOINT[ian].clusterization)
        end
    end


    return nothing

end





function sample_zeta_mergesplit!(rng::MersenneTwister,Likelihood::Likelihood_OU_CircLinmodel, Clusterization::Clusterization_HDPHMM,Likelihood_JOINT::Vector{Likelihood_OU_CircLinmodel})

    print("splitmerge non omplementato bene, copiare con opportune modifiche quello sopra","\n")
    error("")
    #Likelihood      = MCMCLikelihood
    # Clusterization  = MCMCClusterization

    Clust = Clusterization.clusterization
#    Clusterization.pi.parameteracc
    kmax        = Likelihood.kmax
    nc          = Likelihood.data.ncol
    nanim       = Likelihood.data.nanimals
    nt          = Likelihood.data.nt

    zetaprex = 1
    veckacc = zeros(Int16,kmax)
    probsVect       = zeros(Float64,kmax)

    kappaPar = Clusterization.mcmc_ak[1]*Clusterization.mcmc_rho[1]
    alphaPar = Clusterization.mcmc_ak[1]-kappaPar
    zetaprex = 1


    ksamp1_app     = Int16(rand(1:Likelihood.clusterization.n_nonemptyC[1],1)[1])
    ksamp1         = Int16(Likelihood.clusterization.nonemptyC[ksamp1_app])
    ksamp2_app     = Int16(rand(1:Likelihood.clusterization.n_emptyC[1],1)[1])
    ksamp2         = Int16(Likelihood.clusterization.emptyC[ksamp2_app])
    domerge        = Int16(0)


    proposal_split = Float64(0.0);
    proposal_merge = Float64(0.0);
    u = Int16(rand(rng,1:2,1)[1])

    if Likelihood.clusterization.n_nonemptyC[1] == 1
        u = 2
        proposal_merge += log(Float64(0.5))
    end

    if Likelihood.clusterization.n_nonemptyC[1] == kmax
        u = 1
        proposal_split += log(Float64(0.5))
    end

    if u==1
        ## MERGE
        domerge        = Int16(1)
        ksamp2_app     = Int16(rand(1:(Likelihood.clusterization.n_nonemptyC[1]-1),1)[1])
        if ksamp2_app >= ksamp1_app
            ksamp2_app = ksamp2_app+one(ksamp2_app)
        end
        ksamp2         = Int16(Likelihood.clusterization.nonemptyC[ksamp2_app])

        t_index_split1 = Int16.(findall(x->x == ksamp1, Likelihood.clusterization.zeta[1:(nt-1)]))
        t_index_split2 = Int16.(findall(x->x == ksamp2, Likelihood.clusterization.zeta[1:(nt-1)]))

        t_index_merge = sort([t_index_split1;t_index_split2])::Vector{Int16}


        # proposal
        proposal_merge += -log(Float64(Likelihood.clusterization.n_nonemptyC[1]))-log(Float64(Likelihood.clusterization.n_nonemptyC[1]-1))
        proposal_split += -log(Float64(Likelihood.clusterization.n_nonemptyC[1]-1))-log(Float64(Likelihood.clusterization.n_emptyC[1]+1))

        ## n i to j
        nitoj_split1 = deepcopy(Clust.n_itojC[ksamp1])
        nitoj_split2 = deepcopy(Clust.n_itojC[ksamp2])

        nitoj_merge     = deepcopy(Clust.n_itojC[ksamp1])
        nitoj_merge_app = deepcopy(Clust.n_itojC[ksamp2])

        nitoj_merge[ksamp1]     += nitoj_merge_app[ksamp2]
        nitoj_merge[ksamp1]     += nitoj_merge_app[ksamp1]
        nitoj_merge[ksamp1]     += nitoj_merge[ksamp2]
        nitoj_merge_app[ksamp2]  = 0
        nitoj_merge_app[ksamp1]  = 0
        nitoj_merge[ksamp2]      = 0
        nitoj_merge             .+= nitoj_merge_app


    else
        ## split
        t_index_merge = Int16.(findall(x->x == ksamp1, Likelihood.clusterization.zeta[1:(nt-1)]))
        sets_vec      = find_const(Likelihood.clusterization.zeta[1:(nt-1)],ksamp1)

        t_index_split1 = zeros(Int16,length(t_index_merge))
        t_index_split2 = zeros(Int16,length(t_index_merge))

        nsplit1        = Int16(0)
        nsplit2        = Int16(0)

        ## n i to j
        nitoj_merge     = deepcopy(Clust.n_itojC[ksamp1])

        nitoj_split1    = zeros(Int16, kmax)
        nitoj_split2    = zeros(Int16, kmax)

        from_nitoj_split1    = zeros(Int16, kmax)
        from_nitoj_split2    = zeros(Int16, kmax)

        zprec = Int16(1)
        for ii in 1:size(sets_vec,1)

            if sets_vec[ii,1] != 1
                zprec = sets_vec[ii,1]-Int16(1)
            end

            u = Int16(rand(rng,1:2,1)[1])
            if u == 1

                nitoj_split1[ksamp1] += sets_vec[ii,2]-sets_vec[ii,1]
                if sets_vec[ii,2]< (nt-1)
                    nitoj_split1[Likelihood.clusterization.zeta[sets_vec[ii,2]+1]] += Int16(1)
                end

                for jj in sets_vec[ii,1]:sets_vec[ii,2]
                    nsplit1 += one(nsplit1)
                    t_index_split1[nsplit1] = jj

                end

                from_nitoj_split1[Likelihood.clusterization.zeta[zprec]] += Int16(1)

            else
                nitoj_split2[ksamp2] += sets_vec[ii,2]-sets_vec[ii,1]
                if sets_vec[ii,2]< (nt-1)
                    nitoj_split2[Likelihood.clusterization.zeta[sets_vec[ii,2]+1]] += Int16(1)
                end

                for jj in sets_vec[ii,1]:sets_vec[ii,2]
                    nsplit2 += one(nsplit2)
                    t_index_split2[nsplit2] = jj
                end

                from_nitoj_split2[Likelihood.clusterization.zeta[zprec]] += Int16(1)
            end

            #print("\n\n",sum(sets_vec),"\n",sets_vec,"\n\n")
        end

        #print("\n n", nsplit1," ",nsplit2," ",nsplit1+nsplit2," \n")
        t_index_split1 = t_index_split1[1:nsplit1]
        t_index_split2 = t_index_split2[1:nsplit2]
        # proposal
        proposal_merge += -log(Float64(Likelihood.clusterization.n_nonemptyC[1]+1))-log(Float64(Likelihood.clusterization.n_nonemptyC[1]))
        proposal_split += -log(Float64(Likelihood.clusterization.n_nonemptyC[1]))-log(Float64(Likelihood.clusterization.n_emptyC[1]))





    end
    print("\n ", domerge," ", sum(nitoj_merge)==(sum(nitoj_split1)+sum(nitoj_split2))," ", sum(nitoj_merge)," ",sum(nitoj_split1)," ",sum(nitoj_split2), "\n")
    prob_split = 0.0;
    prob_merge = 0.0;
    # for t in t_index_split1
    #     @inbounds prob_split += compute_logpdf(Likelihood, ksamp1, Int16(t+1))
    # end
    # prob_merge = deepcopy(prob_split)

    for t in t_index_split2
        @inbounds prob_split += compute_logpdf(Likelihood_JOINT, ksamp2, Int16(t+1))
        @inbounds prob_merge += compute_logpdf(Likelihood_JOINT, ksamp1, Int16(t+1))
    end




    for kk in 1:kmax

        if kk!=ksamp1
            prob_split += loggamma(alphaPar*Clusterization.mcmc_beta[kk]+ nitoj_split1[kk])
            prob_merge += loggamma(alphaPar*Clusterization.mcmc_beta[kk]+ nitoj_merge[kk])
        else
            prob_split += loggamma(alphaPar*Clusterization.mcmc_beta[kk]+ nitoj_split1[kk]+ kappaPar)
            prob_merge += loggamma(alphaPar*Clusterization.mcmc_beta[kk]+ nitoj_merge[kk]+ kappaPar)
        end
        if kk!=ksamp2
            prob_split += loggamma(alphaPar*Clusterization.mcmc_beta[kk]+ nitoj_split2[kk])
        else
            prob_split += loggamma(alphaPar*Clusterization.mcmc_beta[kk]+ nitoj_split2[kk]+ kappaPar)
        end
    end

    prob_split -= loggamma(alphaPar+ sum(nitoj_split1)+ kappaPar)
    prob_split -= loggamma(alphaPar+ sum(nitoj_split2)+ kappaPar)
    prob_merge -= loggamma(alphaPar+ sum(nitoj_merge)+ kappaPar)


    if domerge == 1
        for kk in 1:kmax

            prob_split += loggamma(alphaPar*Clusterization.mcmc_beta[ksamp1]+ Clust.n_itojC[kk][ksamp1])
            prob_split += loggamma(alphaPar*Clusterization.mcmc_beta[ksamp2]+ Clust.n_itojC[kk][ksamp2])

            prob_merge += loggamma(alphaPar*Clusterization.mcmc_beta[ksamp1]+ Clust.n_itojC[kk][ksamp2]+ Clust.n_itojC[kk][ksamp1])

        end
    else
        for kk in 1:kmax

            prob_split += loggamma(alphaPar*Clusterization.mcmc_beta[ksamp1]+ from_nitoj_split1[kk])
            prob_split += loggamma(alphaPar*Clusterization.mcmc_beta[ksamp2]+ from_nitoj_split2[kk])

            prob_merge += loggamma(alphaPar*Clusterization.mcmc_beta[ksamp1]+ Clust.n_itojC[kk][ksamp1])

        end
    end

    if u==1
        #merge
        MH = (prob_merge+proposal_split)-(prob_split+proposal_merge)

        u = rand(rng,Uniform( 0.0,1.0))

        if u<exp(MH)
            for t in t_index_merge
                @inbounds Likelihood.clusterization.zeta[t] = ksamp1
            end
            Update_ObsInClust!(Likelihood.clusterization)
        end

    else
        #split
        MH = (prob_split+proposal_merge)-(prob_merge+proposal_split)

        u = rand(rng,Uniform( 0.0,1.0))

        if u<exp(MH)
            for t in t_index_split1
                @inbounds Likelihood.clusterization.zeta[t] = ksamp1
            end
            for t in t_index_split2
                @inbounds Likelihood.clusterization.zeta[t] = ksamp2
            end
            Update_ObsInClust!(Likelihood.clusterization)
        end
    end

    if nanim_v2>1
        for ian in 2:nanim_v2
            for t in t_index_merge
                Likelihood_JOINT[ian].clusterization.zeta[t] = Likelihood.clusterization.zeta[t]
            end
            Update_ObsInClust!(Likelihood_JOINT[ian].clusterization)
        end
    end

end
